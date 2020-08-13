#' @export
AcqOptimizerRandomSearch = R6Class("AcqOptimizerRandomSearch",
  inherit = AcqOptimizer,

  public = list(

    param_set = NULL,

    initialize = function() {
      self$param_set = ParamSet$new(list(
        ParamInt$new("iters", lower = 1L, default = 1000L)
      ))
      self$param_set$values$iters = 1000L
    },

    optimize = function(acqf) {
      xdt = generate_design_random(acqfun$search_space, self$param_set$values$iters)$data
      ydt = acqfun$eval_dt(xdt)
      xydt = cbind(xdt, ydt)
      setorderv(xydt, acqf$codomain$ids(), order = -1 * acqf$mult_max_to_min, na.last = TRUE)
      xydt[1, acqf$search_space$ids(), with = FALSE]
    }
))
