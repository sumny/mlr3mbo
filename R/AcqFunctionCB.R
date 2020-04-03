#' @title Acquisition function: Confidence Bound
#'
#' @usage NULL
#' @format [R6::R6Class] object.
#'
#' @section Construction:
#'
#' @section Fields: See [AcqFunction]
#' @section Methods: See [AcqFunction]
#' @export
AcqFunctionCB = R6Class("AcqFunctionCB", 
  inherit = AcqFunction,

  public = list(
    initialize = function(surrogate, domain, codomain) {
      param_set = ParamSet$new(list(
        ParamDbl$new("lambda", lower = 0, default = 2)
      ))
      param_set$values$lambda = 2
      
      super$initialize(
        id = "AcqCB", 
        param_set = param_set,
        surrogate = surrogate,
      )
    },

    set_up = function(optim_instance) {
      super$set_up(optim_instance)
      self$minimize = map_lgl(self$codomain$tags, identical, y = "minimize")
      self$mult
    }

    eval = function(xdt) {
      p = self$surrogate$predict_newdata(task = NULL, newdata = xdt)
      res = p$mean + self$mult_max_to_min * self$settings$lambda * p$se
      data.table(acq_y = res)
    }
  )
)
