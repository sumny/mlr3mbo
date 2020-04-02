AcqOptimizer = R6Class("AcqOptimizer", 
  public = list(
    
    param_set = NULL,

    initialize = function(param_set = ParamSet$new()) {
      assert_param_set(param_set)
      self$param_set = param_set
    },

    optimize = function(acqf) {
      stop("not implemented")
    }
  )
)

AcqOptimizerRandom = R6Class("AcqOptimizerRandom",
  
  inherit = AcqOptimizer,

  public = list(
    initialize = function() {
      param_set = ParamSet$new(list(
        ParamInt$new("n_evals", lower = 0)
      ))
      super$initialize(param_set = param_set)
    },

    optimize = function(acqf) {
      assert_r6(acqf, "AcqFunction")
      n_evals = assert_int(n_evals)
      d = generate_design_random(acqf$domain, self$param_set$values$n_evals)
      ydt = acqf$eval(d$data)
      # FIXME:
      which_best = which_min
      j = which_best(ydt$y)
      d$data[j,]
    }
  )

)
