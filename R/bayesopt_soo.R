#FIXME: should we pass the arvive here, too?

# most important controls:
# - the acqf, with its settings
# - the terminator, thats in objective
# - the design?
# - the optimizer?

# and why do we pass some tings in "control" and some separately?

bayesop_soo = function(instance, acq_function, acq_optimizer) {
  #FIXME maybe do not have this here, but have a general assert helper
  assert_r6(instance, "OptimInstanceSingleCrit")
  assert_r6(acq_function, "AcqFunction")
  assert_r6(acq_optimizer, "AcqOptimizer")
  archive = instance$archive

  #FIXME maybe do not have this here, but have a general init helper
  if (archive$n_evals == 0) {
    design = generate_design_lhs(instance$search_space, 4 * instance$search_space$length)$data
    instance$eval_batch(design)
  }

  acq_function$setup(archive) #setup necessary to determine the domain, codomain (for opt direction) of acq function
  acq_function$surrogate$setup(archive) #maybe not necessary

  repeat {
    xydt = archive$data()
    acq_function$surrogate$update(archive) #update surrogate model with new data
    
    acq_function$update(archive) # NOTE: necessary becaue we have to dertermine e.g. y_best for ei, there are possible other costy calculations that we just want to do once for each state. We might not want to do these calculation in acq_function$eval_dt() because this can get called several times during the optimization.
    # one more costy example would be AEI, where we ask the surrogate for the mean prediction of the points in the design
    # alternatively the update could be called by the AcqOptimizer (but he should not need to know about the archive, so then the archive also has to live in the AcqFun), 
    xdt = acq_optimizer$optimize(acq_function)
    instance$eval_batch(xdt)
    if (instance$is_terminated || instance$terminator$is_terminated(archive)) break
  }

  return(instance$archive)
}

if (FALSE) {
  set.seed(1)
  library(bbotk)
  devtools::load_all()
  library(paradox)
  library(mlr3learners)

  obfun = ObjectiveRFun$new(
    fun = function(xs) sum(unlist(xs)^2),
    domain = ParamSet$new(list(ParamDbl$new("x", -5, 5))),
    id = "test"
  )

  terminator = trm("evals", n_evals = 10)

  instance = MboInstanceSingleCrit$new(
    objective = obfun, 
    terminator = terminator
  )

  surrogate = SurrogateLearner$new(learner = lrn("regr.km", nugget = 0.001), columns = "y")
  acqfun = AcqFunctionEI$new(surrogate = surrogate)
  acqopt = AcqOptimizerRandomSearch$new()

  bayesop_soo(instance, acqfun, acqopt)

  data = instance$archive$data()
  plot(y~batch_nr, data[batch_nr>1,], type = "b")

  xgrid = generate_design_grid(instance$search_space, 100)$data
  preds = cbind(xgrid, acqfun$surrogate$predict(xgrid)$y)
  preds = cbind(preds, acqfun$eval_dt(xgrid))

  acqopt$optimize(acqfun)

  library(ggplot2)
  g = ggplot(data, aes(x = x, y = y, col = batch_nr))
  g = g + geom_line() + geom_point()
  g = g + geom_line(data = preds, aes(x = x, y = mean), col = "blue")
  g = g + geom_line(data = preds, aes(x = x, y = acq_ei), col = "red")
  g

}



