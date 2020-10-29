#' @title Acquisition Function Expected Joint Improvement of Elites
#'
#' @description
#' Expected Joint Improvement of Elites.
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
    #' @param surrogate ([SurrogateSingleCrit]).
    initialize = function(surrogate) {
      param_set = ParamSet$new()
      assert_r6(surrogate, "SurrogateSingleCrit")
      super$initialize("acq_ejie", param_set, surrogate, direction = "maximize")
    },

    #' @description
    #' Evaluates all input values in `xdt`.
    #'
    #' @param xdt [data.table::data.table]
    #'
    #' @return `data.table`
    eval_dt = function(xdt) {

      p = self$surrogate$predict(xdt)
      mu = p$mean
      se = p$se

      ei_j = map(names(self$niches$niches), function(niche) {
        best = self$y_bests[[niche]]
        if (!length(best)) best = self$surrogate_max_to_min * min(unlist(self$y_bests), na.rm = TRUE)
        if (!is.finite(best)) best = 0  # FIXME:
        d = best - self$surrogate_max_to_min * mu
        d_norm = d / se
        ei_j = d * pnorm(d_norm) + se * dnorm(d_norm)
      })

      niche = self$niches$get_niche_dt(self$feature_function_eval_dt(xdt))

      prob_j = if (is.null(self$feature_surrogate_predict)) {
        map(transpose_list(niche), function(x) {
          p_j = rep(0, length(self$niches$niches))
          niche_match = match(x, names(self$niches$niches), nomatch = 0)
          p_j[niche_match] = 1
          p_j
        })
      } else {
        # FIXME: sum to 1
        p_g = self$feature_surrogate_predict(xdt)
        # FIXME: feature_function_ids
        feature_function_ids = names(self$niches$niches[[1L]]$niche_boundaries)
        mu_g = if (is.list(p_g) & !is.data.table(p_g)) map(p_g, "mean") else setNames(list(p_g$mean), nm = feature_function_ids)
        se_g = if (is.list(p_g) & !is.data.table(p_g)) map(p_g, "se") else setNames(list(p_g$se), nm = feature_function_ids)

        if (test_r6(self$niches, classes = "NichesBoundaries")) {
          map(names(self$niches$niches), function(niche) {
            boundaries = self$niches$niches[[niche]]$niche_boundaries
            pj = 1
            for (id in feature_function_ids) {
              pj = pj * (pnorm((mu_g[[id]] - boundaries[[id]][1L]) / se_g[[id]]) - pnorm((mu_g[[id]] - boundaries[[id]][2L]) / se_g[[id]]))
            }
            pj
          })
        }
      }

      ejie = Reduce("+", pmap(list(ei_j, prob_j), function(ec, pc) ec * pc))

      # NOTE: we have to check if a point is in a niche, if not, we do not want to propose it at all, i.e, return ejie = 0      
      ejie[se < 1e-20 | is.na(niche)] = 0
      data.table(acq_ejie = ejie)
    },

    #' @description
    #' Updates acquisition function and sets `y_best`.
    #'
    #' @param archive [bbotk::ArchiveQDO]
    update = function(archive) {
      super$update(archive)
      self$y_bests = setNames(map(names(self$niches$niches), function(niche) archive$best(j = niche)[[archive$cols_y]]), nm = names(self$niches$niches))
    }
  )
)
