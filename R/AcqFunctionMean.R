#' @title Acquisition function: Posterior Mean
#'
#' @usage NULL
#' @format [R6::R6Class] object.
#'
#' @section Construction:
#'
#' @section Fields: See [AcqFunction]
#' @section Methods: See [AcqFunction]
#' @export
AcqFunctionMean = R6Class("AcqFunctionMean", 
  inherit = AcqFunction,
  
  public = list(

    initialize = function(surrogate, domain, codomain) {
      super$initialize(
        id = "AcqMean", 
        surrogate = surrogate,
        domain = domain,
        codomain = codomain
      )
    },

    eval = function(xdt) {
      p = self$surrogate$predict_newdata(xdt)
      data.table(y = p$response)
    }
  )
)

