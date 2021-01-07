#' @export
bayesopt_bop = function(instance, acq_function, acq_optimizer, n_design = 4 * instance$search_space$length) {
  # FIXME: maybe do not have this here, but have a general assert helper
  assert_r6(instance, "OptimInstanceQDO")
  assert_r6(acq_function, "AcqFunction")
  #assert_r6(acq_optimizer, "AcqOptimizer")
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

  acq_function$setup(archive)  # setup necessary to determine the domain, codomain (for opt direction) of acq function

  repeat {
    acq_function$surrogate$update(xydt = char_to_fct(archive_xyg(archive)), y_cols = c(archive$cols_y, archive$cols_g))  # update surrogate model with new data
    acq_function$update(archive)
    xdt = acq_optimizer$optimize(acq_function)
    #xdt = acq_optimizer$optimize(acq_function, archive = archive)
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

  nb1 = NicheBoundaries$new("niche1", niche_boundaries = list(g = c(0, 0.5)))
  nb2 = NicheBoundaries$new("niche2", niche_boundaries = list(g = c(0.5, 0.9)))
  nb = NichesBoundaries$new("test", niches_boundaries = list(niche1 = nb1, niche2 = nb2))

  obfun = ObjectiveRFun$new(
    fun = function(xs) {
      g = xs[[1L]]^2 / 100
      list(
        y = -sin(xs[[1L]] / 2L) + cos(3L * xs[[1L]]),
        g = g,
        niche = nb$get_niche(g)
      )
    },
    domain = ParamSet$new(list(ParamDbl$new("x", 0, 10))),
    codomain = ParamSet$new(list(
      ParamDbl$new("y", tags = "minimize"),
      ParamDbl$new("g", tags = "feature"),
      ParamFct$new("niche", levels = c("niche1", "niche2"), special_vals = list(NA_character_), tags = "niche"))
    ),
    properties = "single-crit",
    id = "test"
  )
  #obfun$eval_dt(data.table(x = c(0, 1.4, 3, 8)))

  terminator = trm("evals", n_evals = 20)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    terminator = terminator
  )

  surrogate = lrn("regr.km")
  surrogate$param_set$values = list(covtype = "matern3_2", optim.method = "gen", nugget.stability = 10^-8)

  acq_function = AcqFunctionEJIE$new(SurrogateMultiCritLearners$new(list(surrogate$clone(deep = TRUE), surrogate$clone(deep = TRUE))), niches = nb)
  acq_optimizer = AcqOptimizerMIES_old$new()
  #n_design = 4 * instance$search_space$length

  bayesopt_bop(instance, acq_function, acq_optimizer)

  x = seq(from = 0.01, to = 9.99, length.out = 1001)
  dt = obfun$eval_dt(data.table(x = x))
  y = dt[["y"]]
  g = dt[["g"]]
  niche = dt[["niche"]]

  cols = c("red", "green")
  niches = c("niche1", "niche2")
  p = acq_function$surrogate$predict(data.table(x = x))
  y_ = p[["y"]][["mean"]]
  g_ = p[["g"]][["mean"]]

  plot(x, y, type = "l")
  points(x, y_, type = "l", col = "red")

  points(x, g, type = "l")
  points(x, g_, type = "l", col = "green")

  plot(x, y, type = "n")
  points(x[is.na(niche)], y[is.na(niche)], col = "black", cex = 0.5, type = "l")
  for (i in 1:2) {
    points(x[niche == niches[i]], y[niche == niches[i]], col = cols[i], cex = 0.5, type = "l")
    best = instance$archive$best(niche = niches[i])
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

  nb1 = NicheBoundaries$new("niche1", niche_boundaries = list(g1 = c(0, 0.5), g2 = c(0, 2)))
  nb2 = NicheBoundaries$new("niche2", niche_boundaries = list(g1 = c(0.5, 0.9), g2 = c(2, 3)))
  nb = NichesBoundaries$new("test", niches_boundaries = list(niche1 = nb1, niche2 = nb2))

  obfun = ObjectiveRFun$new(
    fun = function(xs) {
      g1 = xs[[1L]]^2 / 100
      g2 = sqrt(xs[[1L]])
      list(
        y = -sin(xs[[1L]] / 2L) + cos(3L * xs[[1L]]),
        g1 = g1,
        g2 = g2,
        niche = nb$get_niche(list(g1 = g1, g2 = g2))
      )
    },
    domain = ParamSet$new(list(ParamDbl$new("x", 0, 10))),
    codomain = ParamSet$new(list(
      ParamDbl$new("y", tags = "minimize"),
      ParamDbl$new("g1", tags = "feature"),
      ParamDbl$new("g2", tags = "feature"),
      ParamFct$new("niche", levels = c("niche1", "niche2"), special_vals = list(NA_character_), tags = "niche"))
    ),
    properties = "single-crit",
    id = "test"
  )
  #obfun$eval_dt(data.table(x = c(0, 1.4, 3, 8)))

  terminator = trm("evals", n_evals = 20)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    terminator = terminator
  )

  surrogate = lrn("regr.km")
  surrogate$param_set$values = list(covtype = "matern3_2", optim.method = "gen", nugget.stability = 10^-8)

  acq_function = AcqFunctionEJIE$new(SurrogateMultiCritLearners$new(list(surrogate$clone(deep = TRUE), surrogate$clone(deep = TRUE), surrogate$clone(deep = TRUE))), niches = nb)
  acq_optimizer = AcqOptimizer$new(opt("random_search", batch_size = 1000), trm("evals", n_evals = 1000))
  #n_design = 4 * instance$search_space$length

  bayesopt_bop(instance, acq_function, acq_optimizer)

  x = seq(from = 0.01, to = 9.99, length.out = 1001)
  dt = obfun$eval_dt(data.table(x = x))
  y = dt[["y"]]
  g1 = dt[["g1"]]
  g2 = dt[["g2"]]
  niche = dt[["niche"]]

  cols = c("red", "green")
  niches = c("niche1", "niche2")
  p = acq_function$surrogate$predict(data.table(x = x))
  y_ = p[["y"]][["mean"]]
  g1_ = p[["g1"]][["mean"]]
  g2_ = p[["g2"]][["mean"]]

  df = data.frame(x = x, y = y, col = ifelse(is.na(niche), yes = "black", no = ifelse(niche == "niche1", yes = "red", no = "green")))

  best_x = instance$archive$best()[["x"]]
  best_y = instance$archive$best()[["y"]]

  df_point = data.frame(x = best_x, y = best_y)

  ggplot(df, aes(x = x, y = y)) +
    geom_line(aes(colour = col, group = 1)) +
    scale_colour_identity() +
    geom_point(aes(x = x, y = y), df_point)
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
  #library(mlr3learners.lightgbm)

  source("/home/lps/nb_reticulate.R")

  future::plan(future::multicore)

  #lgr::get_logger("mlr3")$set_threshold("warn")

  # FIXME: allow for quantile definiton of niche storing times in a container
  nb1 = NicheBoundaries$new("niche1", niche_boundaries = list(time = c(log(1), log(3000))))
  nb2 = NicheBoundaries$new("niche2", niche_boundaries = list(time = c(log(3000), log(5000))))
  nb3 = NicheBoundaries$new("niche3", niche_boundaries = list(time = c(log(5000), log(10000))))
  nb = NichesBoundaries$new("test", niches_boundaries = list(niche1 = nb1, niche2 = nb2, niche3 = nb3))

  # FIXME: parallel
  obfun = ObjectiveRFun$new(
    fun = function(xs) {
      psvals = insert_named(xs, map(ps$params[ps$tags == "constant"], "default"))[ps$ids()]
      psvals = Filter(Negate(is.null), psvals)
      names(psvals) = fix_name(names(psvals), mode = "r_to_py")
      configspace_config = py$ConfigSpace$Configuration(py$configspace, psvals)
      time = log(py$runtime_model$predict(config = configspace_config, representation = "configspace"))
      list(
        performance = py$performance_model$predict(config = configspace_config, representation = "configspace", with_noise = TRUE),
        time = time,
        niche = nb$get_niche(time)
      )
    },
    domain = ps_tune,
    codomain = ParamSet$new(list(
      ParamDbl$new("performance", tags = "maximize"),
      ParamDbl$new("time", tags = "feature"),
      ParamFct$new("niche", levels = c("niche1", "niche2", "niche3"), special_vals = list(NA_character_), tags = "niche"))
    ),
    properties = "single-crit",
    check_values = FALSE,
    id = "test"
  )

  terminator = trm("evals", n_evals = 300)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    terminator = terminator
  )

  surrogate = GraphLearner$new(po("imputeoor")  %>>% lrn("regr.ranger", num.trees = 500L, keep.inbag = TRUE, se.method = "jack"))

  # FIXME: RegrAVG SE aggregation
  #surrogate_rt = GraphLearner$new(po("imputeoor") %>>% po("encode") %>>% lrn("regr.lightgbm"))

  acq_function = AcqFunctionEJIE$new(SurrogateMultiCritLearners$new(list(surrogate$clone(deep = TRUE), surrogate$clone(deep = TRUE))), niches = nb)
  acq_optimizer = AcqOptimizerRandomSearch_old$new()
  #acq_optimizer = AcqOptimizer$new(opt("random_search", batch_size = 1000), trm("evals", n_evals = 1000))  # FIXME: debug
  #acq_optimizer = AcqOptimizerMutateCrossover$new()
  #n_design = 4 * instance$search_space$length
  n_design = 50L

  bayesopt_bop(instance, acq_function, acq_optimizer, n_design = 50L)

  # FIXME: rs, other qd-algos
  rs = OptimizerRandomSearch$new()
  rs$param_set$values$batch_size = 1L
  rs$optimize(instance)
}

