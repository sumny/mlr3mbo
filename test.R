library(mlr3)
library(paradox)
# library(bbotk)
# roxygenize()
if (fs::dir_exists("~/cos/bbotk")) {
  load_all("~/cos/bbotk")
} else if (fs::dir_exists("~/gits/bbotk")) {
  load_all("~/gits/bbotk")
}
load_all()

PS_2D = ParamSet$new(list(
  ParamDbl$new("x1", lower = -1, upper = 1),
  ParamDbl$new("x2", lower = -1, upper = 1)
))
FUN_2D = function(xs) {
  y = sum(as.numeric(xs)^2)
  list(y = y)
}
OBJ_2D = ObjectiveRFun$new(fun = FUN_2D, domain = PS_2D)

term = TerminatorEvals$new()
term$param_set$values$n_evals = 9

inst = OptimInstance$new(objective = OBJ_2D, param_set = PS_2D, terminator = term)


ll = lrn("regr.rpart")
acqf = AcqFunctionMean$new()
acqf_optim = AcqOptimizerRandom$new()

a = bayesop_soo(inst, ll, acqf)
print(a)

