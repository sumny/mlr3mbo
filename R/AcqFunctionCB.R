#' @title Acquisition Function Confidence Bound
#'
#' @description
#'  Confidence Bound.
#'
#' @section Parameters:
#' * `lambda` (`numeric(1)`)\cr
#'   TODO DESCRIPTION and Reference
#'
#'
#' @family Acquisition Function
#'
#' @export
AcqFunctionCB = R6Class("AcqFunctionCB",
  inherit = AcqFunction,

  public = list(

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param surrogate [SurrogateSingleCrit].
    initialize = function(surrogate) {
      assert_r6(surrogate, "SurrogateSingleCrit")

      constants = ParamSet$new(list(
        ParamDbl$new("lambda", lower = 0, default = 2)
      ))
      constants$values$lambda = 2

      fun = function(xdt) {
        p = self$surrogate$predict(xdt)
        res = p$mean - self$surrogate_max_to_min * self$constants$values$lambda * p$se
        # FIXME: what do we return here? do we want to see se, mean, too?
        data.table(acq_cb = res)
      }

      super$initialize("acq_cb", constants, surrogate, direction = "same", fun = fun)
    }
  )
)
