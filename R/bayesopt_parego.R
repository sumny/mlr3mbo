bayesop_parego = function(instance, acq_function, acq_optimizer, q = 1, s = 100, rho = 0.05) {
  #FIXME maybe do not have this here, but have a general assert helper
  assert_r6(instance, "OptimInstanceMultiCrit")
  assert_r6(acq_function, "AcqFunction")
  assert_r6(acq_optimizer, "AcqOptimizer")
  assert_int(q, lower = 1)

  archive = instance$archive
  #FIXME maybe do not have this here, but have a general init helper
  if (archive$n_evals == 0) {
    design = generate_design_lhs(instance$search_space, 4 * instance$search_space$length)$data
    instance$eval_batch(design)
  }

  #acq_function$setup(archive) #FIXME: Will fail internal update call for surrogate, or y_min if acq_fun is EI because archive as multiple y columns
  # Alternatives
  # - only setup with domain, codomain, ycolumns (not cool because often we will pass things we do not need and blow up signature)
  # - Implement ArchiveScalarized: cool?
  # - ?
  #manual fix: #FIXMEL Write ParEGO Infill Crit?
  acq_function$update()
  acq_function$search_space = instance$search_space
  acq_function$codomain = generate_acq_codomain(codomain, direction = "minimize", id = acq_function$id)
  acq_function$mult_max_to_min = 1 # we always minimize
  
  d = archive$codomain$length

  repeat {

    ydt = Map("*", ydt, mult_max_to_min(archive$codomain))

    xdt = map_dtr(seq_len(q), function(i) {
      # FIXME: Wrong way to calculate lambda
      lambda = runif(d)
      lambda = lambda / sum(lambda)
      lambda = c(1,0)
      
      mult = Map('*', ydt, lambda)
      y_scal = do.call('+', mult)
      archive$set_column("parego_y", y_scal)
      acq_function$surrogate$update(archive) #update surrogate model with new data
      #acq_function$update(archive) #FIXME e.g. EI will fail w/o update
      ydt = archive$ydt(archive)
      acq_optimizer$optimize(acq_function)  
    })

    instance$eval_batch(xdt)
    if (instance$is_terminated || instance$terminator$is_terminated(instance$archive)) break
  }

  return(instance$archive)
}

if (FALSE) {
  set.seed(1)
  devtools::load_all()
  library(bbotk)
  library(paradox)
  library(mlr3learners)

  FUN_2D_2D = function(xs) {
    list(y1 = xs[[1]]^2, y2 = -xs[[2]]^2)
  }
  PS_2D = ParamSet$new(list(
    ParamDbl$new("x1", lower = -1, upper = 1),
    ParamDbl$new("x2", lower = -1, upper = 1)
  ))
  FUN_2D_2D_CODOMAIN = ParamSet$new(list(
    ParamDbl$new("y1", tags = "minimize"),
    ParamDbl$new("y2", tags = "maximize")
  ))
  obfun = ObjectiveRFun$new(fun = FUN_2D_2D, domain = PS_2D,
    codomain = FUN_2D_2D_CODOMAIN, properties = "multi-crit")

  terminator = trm("evals", n_evals = 20)

  instance = MboInstanceMultiCrit$new(
    objective = obfun, 
    terminator = terminator
  )

  surrogate = SurrogateLearner$new(learner = lrn("regr.km"), columns = "parego_y")
  acqfun = AcqFunctionEI$new(surrogate = surrogate)
  acqopt = AcqOptimizerRandomSearch$new()

  bayesop_parego(instance, acqfun, acqopt, q = 2)

  archdata = instance$archive$data()
  archdata = cbind(archdata[1:18,], instance$archive$temp_cols)
  yhat = surrogate$predict(instance$archive$xydt(y = character()))
  archdata = cbind(archdata, yhat[[1]][1:18, ])

  xgrid = generate_design_grid(instance$search_space, 100)$data
  preds = cbind(xgrid, acqfun$surrogate$predict(xgrid)$parego_y)
  preds = cbind(preds, acqfun$eval_dt(xgrid))

  library(ggplot2)
  g = ggplot(archdata, aes_string(x = "x1", y = "y1", color = "batch_nr"))
  g + geom_point()
}




