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
  acq_function$setup(archive, instance$feature$niche_boundaries$niche_boundaries, instance$feature$surrogate$predict)  # setup necessary to determine the domain, codomain (for opt direction) of acq function

  repeat {
    xydt = archive$data()
    # FIXME: surrogate can also use the archive data now
    instance$feature$surrogate$update(xydt = cbind(xydt[, archive$cols_x, with = FALSE], instance$feature$feature_function$eval_dt(xydt[, archive$cols_x, with = FALSE])), y_cols = instance$feature$feature_function$codomain$ids()) # update feature function surrogate model with new data
    acq_function$surrogate$update(xydt = xydt[, c(archive$cols_x, archive$cols_y), with = FALSE], y_cols = archive$cols_y) # update surrogate model with new data

    acq_function$update(archive) # NOTE: necessary, see bayesop_soo
    xdt = acq_optimizer$optimize(acq_function)
    instance$eval_batch(xdt)
    if (instance$is_terminated || instance$terminator$is_terminated(archive)) break
  }

  return(instance$archive)
}

if (FALSE) {
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

  nb1 = NicheBoundary$new("1", niche_boundary = list(g = c(0, 0.5)))
  nb2 = NicheBoundary$new("2", niche_boundary = list(g = c(0.5, 1)))

  nb = NicheBoundaries$new("test", niche_boundaries = list(nb1, nb2))

  surrogate = lrn("regr.km")
  surrogate$param_set$values <- list(covtype = "matern3_2", optim.method = "gen", jitter = 0)

  ftfun = Feature$new("test", ffun, nb, SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))

  terminator = trm("evals", n_evals = 20)

  instance = OptimInstanceQDOSingleCrit$new(
    objective = obfun,
    feature = ftfun,
    terminator = terminator
  )

  acq_function = AcqFunctionEJIE$new(SurrogateSingleCritLearner$new(surrogate$clone(deep = TRUE)))
  acq_optimizer = AcqOptimizerRandomSearch$new()

  bayesop_bop(instance, acqfun, acqopt)

  x = seq(from = 0.01, to = 9.99, length.out = 1001)
  y = obfun$eval_dt(data.table(x = x))[[1L]]
  g = ffun$eval_dt(data.table(x = x))[[1L]]
  niche = ftfun$niche_boundaries$get_niche_dt(data.table(g))[[1L]]

  cols = c("red", "green")
  niches = as.character(1:2)
  y_ = acq_function$surrogate$predict(data.table(x = x))[[1]]
  g_ = instance$feature$surrogate$predict(data.table(x = x))[[1]]

  #plot(x, y, type = "l")
  #points(x, y_, type = "l", col = "red")
  
  #points(x, g, type = "l")
  #points(x, g_, type = "l", col = "red")

  plot(x, y, type = "n")
  for (i in 1:2) {
    points(x[niche == niches[i]], y[niche == niches[i]], col = cols[i], cex = 0.5, type = "l")
    best = instance$archive$best(c = niches[i])
    points(best[[1L]], best[[2L]], pch = 17)
  }
}

