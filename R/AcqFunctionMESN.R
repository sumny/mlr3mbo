#' @title Acquisition Function Max-value Entropy Search for Niches
#'
#' @description
#' Max-value Entropy Search for Niches
#'
#' TODO DESCRIPTION and Reference
#'
#' @family Acquisition Function
#'
#' @export
AcqFunctionMESN = R6Class("AcqFunctionMESN",
  inherit = AcqFunction,
  public = list(

    #' @field niches ([Niches]).
    niches = NULL,

    #' @field bests ([data.table::data.table]).
    bests = NULL,

    #' @field archive_data ([data.table::data.table]).
    archive_data = NULL,

    #' @field grid [data.table::data.table]
    grid = NULL,

    #' @field maxesn (`list()`).
    maxesn = NULL,

    #' @field cols_x (`character()`).
    #' #FIXME: Name
    cols_x = NULL,

    #' @field cols_y (`character()`).
    #' #FIXME: Name
    cols_y = NULL,

    #' @field cols_g (`character()`).
    #' #FIXME: Name
    #' #FIXME: must match ids in NichesBoundaries
    cols_g = NULL,

    #' @field cols_niche (`character()`).
    #' #FIXME: Name
    #' #FIXME: must match ids in NichesBoundaries
    cols_niche = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param surrogate [SurrogateSingleCrit].
    #' @param niches ([Niches]).
    initialize = function(surrogate, niches) {
      assert_r6(surrogate, "SurrogateMultiCrit")
      assert_r6(niches, "Niches")
      self$niches = niches


      fun = function(xdt) {
        if (is.null(self$maxesn)) {
          stop("maxesn is not set. Missed to call $update(archive)?")
        }
        p = self$surrogate$predict(xdt)
        mesn = map_dbl(seq_len(NROW(p[[self$cols_y]])), function(i) {
          mu = p[[self$cols_y]]$mean[i]
          se = p[[self$cols_y]]$se[i]
          sum(map_dbl(self$maxesn, function(maxes) {
            if (is.null(maxes)) return(0)  # early exit
            gamma = (maxes - (- self$surrogate_max_to_min[[self$cols_y]] * mu)) / se
            p_gamma = pnorm(gamma)
            mean(((gamma * dnorm(gamma)) / (2 * p_gamma)) - log(p_gamma), na.rm = TRUE)
          }))
        })
        mesn[is.na(mesn)] = 0  # FIXME: check NAs
        data.table(acq_mesn = mesn)
      }

      super$initialize("acq_mesn", surrogate = surrogate, direction = "maximize", fun = fun)
    },

    #' @description
    #' Sets up the acquisition function.
    #'
    #' @param archive [bbotk::Archive].
    setup = function(archive) {
      # FIXME: Should we allow alternative search_space as additional argument?

      # here we can change the optim direction of the codomain for the acq function
      self$codomain = generate_acq_codomain(archive$codomain, id = self$id, direction = self$direction)

      self$surrogate_max_to_min = mult_max_to_min(archive$codomain)

      self$domain = archive$search_space$clone(deep = TRUE)
      #self$domain$trafo = NULL # FIXME is it okay to do this?

      self$cols_x = archive$cols_x
      self$cols_y = archive$cols_y
      self$cols_g = archive$cols_g
      self$cols_niche = archive$cols_niche

      if (is.null(self$grid)) {
        self$grid = generate_design_lhs(self$domain, n = 1000L)$data  # FIXME: this does scale extremely poorly
      }
    },

    #' @description
    #' Updates acquisition function and sets `maxes`.
    #'
    #' @param archive [bbotk::Archive]
    update = function(archive) {
      super$update(archive)
      self$bests = archive$best()
      self$archive_data = archive$data
      self$maxesn = get_maxesn(x = archive$data[, archive$cols_x, with = FALSE], grid = self$grid, surrogate = self$surrogate, surrogate_max_to_min = self$surrogate_max_to_min, cols_y = self$cols_y, cols_g = self$cols_g, niches = self$niches)
    }
  )
)



# FIXME: AcqFunction ParamSet with at least gridsize and nk
get_maxesn = function(nK = 10000L, x, grid, surrogate, surrogate_max_to_min, cols_y, cols_g, niches) {
  xgrid = rbind(grid, x)
  p = surrogate$predict(xgrid)
  mu = p[[cols_y]]$mean
  se = p[[cols_y]]$se
  se[se < .Machine$double.eps] = .Machine$double.eps  # FIXME:

  # FIXME: also used in ejie --> own function
  prob_j = {
  # FIXME: not to model feature function via tags?
  #if (is.null(self$feature_surrogate_predict)) {
  #  niche = self$niches$get_niche_dt(self$feature_function_eval_dt(xdt))
  #  map(transpose_list(niche), function(x) {
  #    p_j = rep(0, length(self$niches$niches))
  #    niche_match = match(x, names(self$niches$niches), nomatch = 0)
  #    p_j[niche_match] = 1
  #    p_j
  #  })
  #} else {
    p_g = p[cols_g]
    mu_g = if (is.list(p_g) & !is.data.table(p_g)) map(p_g, "mean") else setNames(list(p_g$mean), nm = cols_g)
    se_g = if (is.list(p_g) & !is.data.table(p_g)) map(p_g, "se") else setNames(list(p_g$se), nm = cols_g)

    if (test_r6(niches, classes = "NichesBoundaries")) {
      map(names(niches$niches), function(niche) {
        boundaries = niches$niches[[niche]]$niche_boundaries
        pj = 1
        for (id in cols_g) {
          pj = pj * (pnorm((mu_g[[id]] - boundaries[[id]][1L]) / se_g[[id]]) - pnorm((mu_g[[id]] - boundaries[[id]][2L]) / se_g[[id]]))
        }
        pj
      })
    }
  }
  names(prob_j) = names(niches$niches)
  setDT(prob_j)
  argmax_prob_j = c(names(prob_j), NA_character_)[apply(prob_j, MARGIN = 1L, FUN = function(x) {
    tmp = which(x > 0.99)
    if (length(tmp) == 0L) length(x) + 1L else tmp
  })]

  setNames(map(names(niches$niches), function(niche) {

    ids = which(argmax_prob_j == niche)
    if (length(ids) == 0L) {
      return(NULL)  # early exit
    }

    mu_max = max(- surrogate_max_to_min[[cols_y]] * mu[ids])

    left = mu_max
    leftprob = probf(left, mu = mu[ids], se = se[ids], surrogate_max_to_min = surrogate_max_to_min[[cols_y]], prob_j = prob_j[[niche]][ids])
    # FIXME:
    ntry = 0L
    while (leftprob > 0.50 && ntry < 3L) {
      left = if (left > 0.01) left / 2 else 2 * left - 0.05
      leftprob = probf(left, mu = mu[ids], se = se[ids], surrogate_max_to_min = surrogate_max_to_min[[cols_y]], prob_j = prob_j[[niche]][ids])
      ntry = ntry + 1L
    }

    ntry = 0L
    right = max(- surrogate_max_to_min[[cols_y]] * (mu[ids] - (surrogate_max_to_min[[cols_y]] * 5 * se[ids])))
    rightprob = probf(right, mu = mu[ids], se = se[ids], surrogate_max_to_min = surrogate_max_to_min[[cols_y]], prob_j = prob_j[[niche]][ids])
    while (rightprob < 0.50 && ntry < 3L) {
      right = right + right - left
      rightprob = probf(right, mu = mu[ids], se = se[ids], surrogate_max_to_min = surrogate_max_to_min[[cols_y]], prob_j = prob_j[[niche]][ids])
      ntry = ntry + 1L
    }

    if (leftprob > 0.50 || rightprob < 0.50) {
      return(mu_max + runif(nK, min = 0, max = 1))
    }

    mgrid = seq(from = left, to = right, length.out = 100L)

    prob = map_dbl(mgrid, function(mg) {
      prod(pnorm((mg - (- surrogate_max_to_min[[cols_y]] * mu[ids])) / se[ids]))
    })

    if (sum(prob > 0.05 & prob < 0.95) == 0L) {
      return(mu_max + runif(nK, min = 0, max = 1))
    }

    # Gumbel sampling
    q1 = optimize(function(x) abs(probf(x, mu = mu[ids], se = se[ids], surrogate_max_to_min = surrogate_max_to_min[[cols_y]], prob_j = prob_j[[niche]][ids]) - 0.25), interval = range(mgrid))$minimum
    q2 = optimize(function(x) abs(probf(x, mu = mu[ids], se = se[ids], surrogate_max_to_min = surrogate_max_to_min[[cols_y]], prob_j = prob_j[[niche]][ids]) - 0.5), interval = range(mgrid))$minimum
    q3 = optimize(function(x) abs(probf(x, mu = mu[ids], se = se[ids], surrogate_max_to_min = surrogate_max_to_min[[cols_y]], prob_j = prob_j[[niche]][ids]) - 0.75), interval = range(mgrid))$minimum
    beta = (q1 - q3) / (log(log(4 / 3)) - log(log(4)))  # FIXME: assert beta > 0
    alpha = q2 + beta * log(log(2))

    -log(-log(runif(nK, min = 0, max = 1))) * beta + alpha
    # FIXME: maxes that are <= mu_max + eps should be replaced by mu_max + eps
  }), nm = names(niches$niches))
}



probf = function(mu_, mu, se, surrogate_max_to_min, prob_j) {
  prod(pnorm((mu_ - (- surrogate_max_to_min * mu)) / se) * prob_j)
}

