#' @title Acquisition function: Expected Improvement
#'
#' @usage NULL
#' @format [R6::R6Class] object.
#'
#' @section Construction:
#'
#' @section Fields: See [AcqFunction]
#' @section Methods: See [AcqFunction]
#' @export
AcqFunctionEIPS = R6Class("AcqFunctionTEI",
  inherit = AcqFunction,
  public = list(

    y_best = NULL, 

    initialize = function(surrogate) {
      param_set = ParamSet$new()
      assert_r6(surrogate, "SurrogateCollection") #FIXME: Write assert_surrogate function
      assert_set_equal(surrogate$columns, c("y", "time")) #FIXME: y should be codomain$ids 
      super$initialize("acq_eips", param_set, surrogate, direction = "maximize")
    },

    eval_dt = function(xdt) {
      p = self$surrogate$predict(xdt)
      mu = p[["y"]]$mean #FIXME: y should be codomain$ids
      se = p[["y"]]$se
      mu_t = p[["time"]]
      d = self$y_best - self$mult_max_to_min * mu
      d_norm = d / se
      ei = d * pnorm(d_norm) + se + dnorm(d_norm)
      eips = ei / mu_t
      eips = ifelse(se < 1e-20 | mu_t < 1e-20, 0, eips)
      data.table(acq_ei = ei)
    },

    update = function(archive) {
      super$update(archive)
      self$y_best = archive$best()[[archive$cols_y]]
    }
  )
)
