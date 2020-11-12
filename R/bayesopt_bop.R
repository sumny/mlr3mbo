#' @export
bayesop_bop = function(instance, acq_function, acq_optimizer, n_design = 4 * instance$search_space$length) {
  # FIXME: maybe do not have this here, but have a general assert helper
  assert_r6(instance, "OptimInstanceQDO")
  assert_r6(acq_function, "AcqFunctionQDO")
  assert_r6(acq_optimizer, "AcqOptimizer")
  archive = instance$archive

  # FIXME: maybe do not have this here, but have a general init helper
  if (archive$n_evals == 0) {
    design = if(instance$search_space$has_deps) {
      generate_design_random(instance$search_space, n_design)$data
    } else {
      generate_design_lhs(instance$search_space, n_design)$data
    }
    instance$eval_batch(design)
  }

  # FIXME: better way to pass the feature surrogate
  acq_function$setup(archive, feature_function_eval_dt = instance$feature$feature_function$eval_dt, feature_surrogate_predict = instance$feature$surrogate$predict, niches = instance$feature$niches)  # setup necessary to determine the domain, codomain (for opt direction) of acq function, niche boundaries and surrogate predict

  ps_ch = names(which(instance$objective$domain$storage_type == "character"))
  ps_lgl = names(which(instance$objective$domain$storage_type == "logical"))


  repeat {
    # FIXME:
    xydt = xydt[, c(archive$cols_x, archive$cols_g, archive$cols_y), with = FALSE]
    xydt = setDT(imap(xydt, function(x, name) {
      if (name %in% ps_ch) {
        factor(x, levels = instance$objective$domain$params[[name]]$levels)
      } else if (name %in% ps_lgl) {
        x + 0L
      } else {
        x
      }
    }))

    if (instance$feature$model_feature_function) {
      instance$feature$surrogate$update(xydt = xydt[, c(archive$cols_x, archive$cols_g), with = FALSE], y_cols = archive$cols_g)  # update feature function surrogate model with new data
      cat(c("PE", "R2", "\n"),
        c(instance$feature$surrogate$model$model$regr.ranger$model$prediction.error, instance$feature$surrogate$model$model$regr.ranger$model$r.squared),
        "\n"
      )
    }
    acq_function$surrogate$update(xydt = xydt[, c(archive$cols_x, archive$cols_y), with = FALSE], y_cols = archive$cols_y)  # update surrogate model with new data
    cat(c("PE", "R2", "\n"),
      c(acq_function$surrogate$model$model$regr.ranger$model$prediction.error, acq_function$surrogate$model$model$regr.ranger$model$r.squared),
      "\n"
    )

    acq_function$update(archive)  # NOTE: necessary, see bayesop_soo
    xdt = acq_optimizer$optimize(acq_function)
    instance$eval_batch(xdt)
    if (instance$is_terminated || instance$terminator$is_terminated(archive)) break
  }

  return(instance$archive)
}

if (FALSE) {
  # single feature function example
  set.seed(1)
  devtools::load_all("../bbotk")
  devtools::load_all()
  library(paradox)
  library(mlr3learners)

  obfun = ObjectiveRFun$new(
    fun = function(xs) -sin(xs[[1L]] / 2L) + cos(3L * xs[[1L]]),
    domain = ParamSet$new(list(ParamDbl$new("x", 0, 10))),
    id = "y"
  )

  ffun = ObjectiveRFun$new(
    fun = function(xs) xs[[1L]]^2 / 100,
    domain = ParamSet$new(list(ParamDbl$new("x", 0, 10))),
    codomain = ParamSet$new(list(ParamDbl$new("g", tags = "minimize"))),  # FIXME: minimize / maximize makes no sense for feature functions but is required
    id = "g"
  )

  nb1 = NicheBoundaries$new("niche1", niche_boundaries = list(g = c(0, 0.5)))
  nb2 = NicheBoundaries$new("niche2", niche_boundaries = list(g = c(0.5, 0.9)))

  nb = NichesBoundaries$new("test", niches_boundaries = list(niche1 = nb1, niche2 = nb2))
  #nb$get_niche_dt(ffun$eval_dt(data.table(x = c(0, 1.4, 3, 8))))

  surrogate = lrn("regr.km")
  surrogate$param_set$values = list(covtype = "matern3_2", optim.method = "gen", jitter = 0)

  ftfun = Feature$new("test", ffun, nb, SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))

  terminator = trm("evals", n_evals = 20)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    feature = ftfun,
    terminator = terminator
  )

  acq_function = AcqFunctionEJIE$new(SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))
  acq_optimizer = AcqOptimizerRandomSearch$new()
  # n_design = 4 * instance$search_space$length

  bayesop_bop(instance, acq_function, acq_optimizer)

  x = seq(from = 0.01, to = 9.99, length.out = 1001)
  y = obfun$eval_dt(data.table(x = x))[[1L]]
  g = ffun$eval_dt(data.table(x = x))[[1L]]
  niche = ftfun$niches$get_niche_dt(data.table(g))[[1L]]

  cols = c("red", "green")
  niches = c("niche1", "niche2")
  y_ = acq_function$surrogate$predict(data.table(x = x))[[1]]
  g_ = instance$feature$surrogate$predict(data.table(x = x))[[1]]

  plot(x, y, type = "l")
  points(x, y_, type = "l", col = "red")
  
  points(x, g, type = "l")
  points(x, g_, type = "l", col = "green")

  plot(x, y, type = "n")
  points(x[is.na(niche)], y[is.na(niche)], col = "black", cex = 0.5, type = "l")
  for (i in 1:2) {
    points(x[niche == niches[i]], y[niche == niches[i]], col = cols[i], cex = 0.5, type = "l")
    best = instance$archive$best(j = niches[i])
    points(best[[1L]], best[[2L]], col = "brown", pch = 19)
  }
}

if (FALSE) {
  # multiple feature functions example
  set.seed(1)
  devtools::load_all("../bbotk")
  devtools::load_all()
  library(paradox)
  library(mlr3learners)

  obfun = ObjectiveRFun$new(
    fun = function(xs) -sin(xs[[1L]] / 2L) + cos(3L * xs[[1L]]),
    domain = ParamSet$new(list(ParamDbl$new("x", 0, 10))),
    id = "y"
  )

  ffun = ObjectiveRFun$new(
    fun = function(xs) list(g1 = xs[[1L]]^2 / 100, g2 = sqrt(xs[[1L]])),
    domain = ParamSet$new(list(ParamDbl$new("x", 0, 10))),
    codomain = ParamSet$new(list(ParamDbl$new("g1", tags = "minimize"), ParamDbl$new("g2", tags = "minimize"))),
    id = "g"
  )

  nb1 = NicheBoundaries$new("niche1", niche_boundaries = list(g1 = c(0, 0.5), g2 = c(0, 2)))
  nb2 = NicheBoundaries$new("niche2", niche_boundaries = list(g1 = c(0.5, 0.9), g2 = c(2, 3)))

  nb = NichesBoundaries$new("test", niches_boundaries = list(niche1 = nb1, niche2 = nb2))
  #nb$get_niche_dt(ffun$eval_dt(data.table(x = c(0, 1.4, 3, 8))))

  surrogate = lrn("regr.km")
  surrogate$param_set$values = list(covtype = "matern3_2", optim.method = "gen", jitter = 0)

  ftfun = Feature$new("test", ffun, nb, SurrogateMultiCritLearners$new(list(surrogate$clone(deep = TRUE), surrogate$clone(deep = TRUE))))

  terminator = trm("evals", n_evals = 20)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    feature = ftfun,
    terminator = terminator
  )

  acq_function = AcqFunctionEJIE$new(SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))
  acq_optimizer = AcqOptimizerRandomSearch$new()
  # n_design = 4 * instance$search_space$length

  bayesop_bop(instance, acq_function, acq_optimizer)
  
  x = seq(from = 0.01, to = 9.99, length.out = 1001)
  y = obfun$eval_dt(data.table(x = x))[[1L]]
  g = ffun$eval_dt(data.table(x = x))
  niche = ftfun$niches$get_niche_dt(data.table(g))[[1L]]

  cols = c("red", "green")
  niches = c("niche1", "niche2")
  y_ = acq_function$surrogate$predict(data.table(x = x))[[1]]
  g_ = instance$feature$surrogate$predict(data.table(x = x))

  df = data.frame(x = x, y = y, col = ifelse(is.na(niche), yes = "black", no = ifelse(niche == "niche1", yes = "red", no = "green")))

  best_x = instance$archive$best()[["x"]]
  best_y = obfun$eval_dt(data.table(x = best_x))[[1L]]
  
  df_point = data.frame(x = best_x, y = best_y)

  ggplot(df, aes(x = x, y = y)) + 
    geom_line(aes(colour = col, group = 1)) + 
    scale_colour_identity() +
    geom_point(aes(x = x, y = y), df_point)
}

if (FALSE) {
  # ROC example
  set.seed(1)
  devtools::load_all("../bbotk")
  devtools::load_all()
  library(paradox)
  library(mlr3pipelines)
  library(mlr3learners)
  library(mlr3tuning)
  library(mlr3viz)
  library(precrec)

  lgr::get_logger("mlr3")$set_threshold("warn")

  task = tsk("sonar")
  learner = lrn("classif.ranger", predict_type = "prob")
  resampling = rsmp("cv", folds = 10L)
  measure = msr("classif.auc")

  # FIXME: needs all Tuning classes later
  #obfun = ObjectiveTuning$new(task = task, learner = learner,
  #  resampling = resampling, measures = list(measure),
  #  store_benchmark_result = TRUE,
  #  store_models = TRUE, check_values = TRUE)

  domain = ParamSet$new(list(
    ParamInt$new("mtry", lower = 1L, upper = 60L),
    ParamInt$new("num.trees", lower = 1L, upper = 1000L)))

  #refinstance = TuningInstanceSingleCrit$new(task, learner, resampling, measure, domain, trm("evals", n_evals = 20))
  #tnr("random_search")$optimize(refinstance) # 5, 401, 0.9380148

  obfun = ObjectiveRFun$new(
    fun = function(xs) {
      learner$param_set$values = xs
      rr = resample(task, learner, resampling)
      rr$aggregate(measure)
    },
    domain = domain,
    codomain = ParamSet$new(list(ParamDbl$new("y", tags = "maximize"))),
    id = "y"
  )

  # FIXME: Need an ObjectiveFeature to be more liberal
  ffun = ObjectiveRFun$new(
    fun = function(xs) {
      learner$param_set$values = xs
      rr = resample(task, learner, resampling)
      mod = evalmod(as_precrec(rr), raw_curves = FALSE)
      dmod = as.data.table(mod)[type == "ROC", c("x", "y")]
      #plot(mod, "ROC")
      data.table(roc = list(gx = dmod[["x"]], gy = dmod[["y"]]))
    },
    domain = domain,
    codomain = ParamSet$new(list(ParamDbl$new("roc", tags = "minimize"))),
    id = "groc",
    check_values = FALSE
  )

  e1 = Ellipsoid2D$new("e1", center = c(0., 0.7), radii = c(0.025, 0.025))
  #e2 = Ellipsoid2D$new("e2", center = c(0.4, 0.8), radii = c(0.05, 0.05))
  e3 = Ellipsoid2D$new("e3", center = c(0.3, 1), radii = c(0.025, 0.025))
  #e4 = Ellipsoid2D$new("e4", center = c(0.4, 0.9), radii = c(0.05, 0.1))

  niches = NichesROC$new("test", ellipsoids = list(niche1 = list("e1" = e1), niche2 = list("e3" = e3)))

  surrogate = lrn("regr.ranger")

  ftfun = Feature$new("test", ffun, niches, NULL)

  terminator = trm("evals", n_evals = 1)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    feature = ftfun,
    terminator = terminator
  )

  acq_function = AcqFunctionEJIE$new(SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))
  acq_optimizer = AcqOptimizerRandomSearch$new()
  acq_optimizer$param_set$values$iters = 100
  #n_design = 4 * instance$search_space$length

  bayesop_bop(instance, acq_function, acq_optimizer)

  par(mfrow = c(1, 2))
  learner$param_set$values = list(mtry = 18, num.trees = 896)
  rr = resample(task, learner, resampling)
  plot(evalmod(as_precrec(rr)), "ROC")
  learner$param_set$values = list(mtry = 8, num.trees = 617)
  rr = resample(task, learner, resampling)
  plot(evalmod(as_precrec(rr)), "ROC")
}

if (FALSE) {
  # runtime example
  set.seed(1)
  devtools::load_all("../bbotk")
  devtools::load_all()
  library(paradox)
  library(mlr3pipelines)
  library(mlr3learners)
  library(mlr3tuning)
  library(mlr3viz)

  #lgr::get_logger("mlr3")$set_threshold("warn")

  task = tsk("sonar")
  learner = lrn("classif.ranger", predict_type = "prob")
  resampling = rsmp("cv", folds = 3L)
  measure = msr("classif.acc")

  # FIXME: needs all Tuning classes later
  #obfun = ObjectiveTuning$new(task = task, learner = learner,
  #  resampling = resampling, measures = list(measure),
  #  store_benchmark_result = TRUE,
  #  store_models = TRUE, check_values = TRUE)

  domain = ParamSet$new(list(
    ParamInt$new("mtry", lower = 1L, upper = 60L),
    ParamInt$new("num.trees", lower = 1L, upper = 10000L),
    ParamDbl$new("sample.fraction", lower = 0.1, upper = 1)))

  obfun = ObjectiveRFun$new(
    fun = function(xs) {
      learner$param_set$values = xs
      rr = resample(task, learner, resampling)
      rr$aggregate(measure)
    },
    domain = domain,
    codomain = ParamSet$new(list(ParamDbl$new("y", tags = "maximize"))),
    id = "y"
  )

  # FIXME: Need an ObjectiveFeature to be more liberal
  ffun = ObjectiveRFun$new(
    fun = function(xs) {
      learner$param_set$values = xs
      learner$train(task)      
      learner$timings["train"]
    },
    domain = domain,
    codomain = ParamSet$new(list(ParamDbl$new("time", tags = "minimize"))),
    id = "gtime",
    check_values = FALSE
  )

  # FIXME: allow for quantile definiton of niche storing times in a container
  nb1 = NicheBoundaries$new("niche1", niche_boundaries = list(time = c(0, 1)))
  nb2 = NicheBoundaries$new("niche2", niche_boundaries = list(time = c(1, 10)))

  nb = NichesBoundaries$new("test", niches_boundaries = list(niche1 = nb1, niche2 = nb2))

  surrogate = lrn("regr.km")
  surrogate$param_set$values = list(covtype = "matern3_2", optim.method = "gen", jitter = 0)

  ftfun = Feature$new("test", ffun, nb, SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))

  terminator = trm("evals", n_evals = 50)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    feature = ftfun,
    terminator = terminator
  )

  acq_function = AcqFunctionEJIE$new(SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))
  acq_optimizer = AcqOptimizerRandomSearch$new()
  #n_design = 4 * instance$search_space$length

  bayesop_bop(instance, acq_function, acq_optimizer)

  # niche2 objective
  p1 = acq_function$surrogate$predict(data.table(mtry = 1:60, num.trees = 6017, sample.fraction = 0.7519022))
  p1$lwr = p1$mean - 1.96 * p1$se
  p1$upr = p1$mean + 1.96 * p1$se
  p1$mtry = 1:60
  ggplot(p1, aes(mtry, mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha =0.3) +
    labs(x = "mtry") +
    labs(y = "Mean classif.acc, num.trees = 3707, sample.fraction = 0.883")

  p2 = acq_function$surrogate$predict(data.table(mtry = 24, num.trees = 1:10000, sample.fraction = 0.7519022))
  p2$lwr = p2$mean - 1.96 * p2$se
  p2$upr = p2$mean + 1.96 * p2$se
  p2$num.trees = 1:10000
  ggplot(p2, aes(num.trees, mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha =0.3) +
    labs(x = "num.trees") +
    labs(y = "Mean Fitting Time, mtry = 54, sample.fraction = 0.883")

  p3 = acq_function$surrogate$predict(data.table(mtry = 24, num.trees = 9068, sample.fraction = seq(0, 1, length.out = 1000)))
  p3$lwr = p3$mean - 1.96 * p3$se
  p3$upr = p3$mean + 1.96 * p3$se
  p3$sample.fraction = seq(0, 1, length.out = 1000)
  ggplot(p3, aes(sample.fraction, mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha =0.3) +
    labs(x = "mtry") +
    labs(y = "Mean Fitting Time, mtry = 54, num.trees = 3707")

  # niche2 feature
  p1 = instance$feature$surrogate$predict(data.table(mtry = 1:60, num.trees = 9068, sample.fraction = 0.7519022))
  p1$lwr = p1$mean - 1.96 * p1$se
  p1$upr = p1$mean + 1.96 * p1$se
  p1$mtry = 1:60
  ggplot(p1, aes(mtry, mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha =0.3) +
    labs(x = "mtry") +
    labs(y = "Mean Fitting Time, num.trees = 3707, sample.fraction = 0.883")

  p2 = instance$feature$surrogate$predict(data.table(mtry = 54, num.trees = 1:10000, sample.fraction = 0.8825619))
  p2$lwr = p2$mean - 1.96 * p2$se
  p2$upr = p2$mean + 1.96 * p2$se
  p2$num.trees = 1:10000
  ggplot(p2, aes(num.trees, mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha =0.3) +
    labs(x = "num.trees") +
    labs(y = "Mean Fitting Time, mtry = 54, sample.fraction = 0.883")

  p3 = instance$feature$surrogate$predict(data.table(mtry = 54, num.trees = 3707, sample.fraction = seq(0, 1, length.out = 1000)))
  p3$lwr = p3$mean - 1.96 * p3$se
  p3$upr = p3$mean + 1.96 * p3$se
  p3$sample.fraction = seq(0, 1, length.out = 1000)
  ggplot(p3, aes(sample.fraction, mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha =0.3) +
    labs(x = "mtry") +
    labs(y = "Mean Fitting Time, mtry = 54, num.trees = 3707")
}

if (FALSE) {
  # nasbench example
  set.seed(1)
  devtools::load_all("../bbotk")
  devtools::load_all()
  library(paradox)
  library(mlr3pipelines)
  library(mlr3learners)
  library(mlr3tuning)
  library(mlr3learners.lightgbm)

  source("/home/user/nb_reticulate.R")

  future::plan(future::multicore)

  #lgr::get_logger("mlr3")$set_threshold("warn")

  #FIXME: surrogate

  domain = ps_tune

  # FIXME: parallel
  obfun = ObjectiveRFun$new(
    fun = function(xs) {
      psvals = insert_named(xs, map(ps$params[ps$tags == "constant"], "default"))
      #psvals = Filter(Negate(is.null), psvals)
      names(psvals) = fix_name(names(psvals), mode = "r_to_py")
      configspace_config = py$ConfigSpace$Configuration(py$configspace, psvals)
      py$performance_model$predict(config = configspace_config, representation = "configspace", with_noise = TRUE)
    },
    domain = ps_tune,
    codomain = ParamSet$new(list(ParamDbl$new("performance", tags = "maximize"))),
    id = "yperformance",
    check_values = FALSE
  )

  # FIXME: Need an ObjectiveFeature to be more liberal
  ffun = ObjectiveRFun$new(
    fun = function(xs) {
      psvals = insert_named(xs, map(ps$params[ps$tags == "constant"], "default"))
      #psvals = Filter(Negate(is.null), psvals)
      names(psvals) = fix_name(names(psvals), mode = "r_to_py")
      configspace_config = py$ConfigSpace$Configuration(py$configspace, psvals)
      log(py$runtime_model$predict(config = configspace_config, representation = "configspace"))
    },
    domain = ps_tune,
    codomain = ParamSet$new(list(ParamDbl$new("runtime", tags = "minimize"))),
    id = "gruntime",
    check_values = FALSE
  )

  # FIXME: allow for quantile definiton of niche storing times in a container
  nb1 = NicheBoundaries$new("niche1", niche_boundaries = list(time = c(log(1), log(3000))))
  nb2 = NicheBoundaries$new("niche2", niche_boundaries = list(time = c(log(3000), log(5000))))
  nb3 = NicheBoundaries$new("niche3", niche_boundaries = list(time = c(log(5000), log(10000))))

  nb = NichesBoundaries$new("test", niches_boundaries = list(niche1 = nb1, niche2 = nb2, niche3 = nb3))

  surrogate = GraphLearner$new(po("imputeoor")  %>>% lrn("regr.ranger"))
  surrogate$param_set$values$regr.ranger.se.method = "jack"
  surrogate$param_set$values$regr.ranger.keep.inbag = TRUE

  # FIXME: RegrAVG SE aggregation
  #surrogate_rt = GraphLearner$new(po("imputeoor") %>>% po("encode") %>>% lrn("regr.lightgbm"))

  ftfun = Feature$new("test", ffun, nb, SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))

  terminator = trm("evals", n_evals = 3000)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    feature = ftfun,
    terminator = terminator
  )

  acq_function = AcqFunctionEJIE$new(SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))
  acq_optimizer = AcqOptimizerRandomSearch$new()
  acq_optimizer$param_set$values$iters = 10000
  #n_design = 4 * instance$search_space$length

  bayesop_bop(instance, acq_function, acq_optimizer)
 
  # FIXME: rs, other qd-algos
  #rs = OptimizerRandomSearch$new()
  #rs$optimize(instance)
}

