if (FALSE) {
  set.seed(1)
  library(bbotk)
  devtools::load_all()
  library(paradox)
  library(mlr3learners)

  obfun = ObjectiveRFun$new(
    fun = function(xs) list(y = sum(unlist(xs)^2), time = 1),
    domain = ParamSet$new(list(ParamDbl$new("x", -5, 5))),
    codomain = ParamSet$new(list(ParamDbl$new("y"))),
    id = "test"
  )

  terminator = trm("evals", n_evals = 20)

  instance = MboInstanceSingleCrit$new(
    objective = obfun, 
    terminator = terminator
  )

  surrogate = SurrogateCollection$new(list(
    SurrogateLearner$new(learner = lrn("regr.km"), columns = "y"),
    SurrogateLearner$new(learner = lrn("regr.km"), columns = "time")
  ))
  acqfun = AcqFunctionEIPS$new(surrogate = surrogate)
  acqopt = AcqOptimizerRandomSearch$new()

  bayesop_soo(instance, acqfun, acqopt)
}



