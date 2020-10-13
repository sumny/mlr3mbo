#' @title Acquisition Function Expected Joint Improvement of ELites
#'
#' @description
#' Expected Improvement.
#'
#' @export
AcqFunctionEJIE = R6Class("AcqFunctionEJIE",
  inherit = AcqFunctionQDO,
  public = list(

    #' @field y_bests (`numeric()`).
    y_bests = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param surrogate [SurrogateSingleCrit].
    initialize = function(surrogate) {
      param_set = ParamSet$new()
      assert_r6(surrogate, "SurrogateSingleCrit")
      super$initialize("acq_ei", param_set, surrogate, direction = "maximize")
    },

    #' @description
    #' Evaluates all input values in `xdt`.
    #'
    #' @param xdt [data.table::data.table]
    #'
    #' @return `data.table`
    eval_dt = function(xdt) {
      # FIXME: currently only 1d features
      p = self$surrogate$predict(xdt)
      mu = p$mean
      se = p$se

      ei_c = map(self$niches, function(niche) {
        best = self$y_bests[[niche]]
        if (is.na(best)) best = self$surrogate_max_to_min * min(unlist(self$y_bests), na.rm = TRUE)
        d = best - self$surrogate_max_to_min * mu
        d_norm = d / se
        ei_c = d * pnorm(d_norm) + se * dnorm(d_norm)
      })

      p_g = self$feature_surrogate_predict(xdt)
      mu_g = p_g$mean
      se_g = p_g$se

      prob_c = map(self$niches, function(niche) {
        boundary = self$niche_boundaries[[niche]]$niche_boundary[[1L]]
        pnorm((mu_g - boundary[1L]) / se_g) - pnorm((mu_g - boundary[2L]) / se_g)
      })

      ei = Reduce("+", pmap(list(ei_c, prob_c), function(ec, pc) ec * pc))

      ei = ifelse(se < 1e-20, 0, ei)
      data.table(acq_ei = ei)
    },

    #' @description
    #' Updates acquisition function and sets `y_best`.
    #'
    #' @param archive [bbotk::Archive]
    update = function(archive) {
      super$update(archive)
      self$y_bests = map(self$niches, function(niche) archive$best(c = niche)[[archive$cols_y]])
      names(self$y_bests) = self$niches
    }
  )
)
