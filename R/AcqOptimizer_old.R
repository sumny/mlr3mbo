#' @title Acquisition Optimizer
#'
#' @description
#' Optimizer for [AcqFunction] objects.
#'
#' @export
AcqOptimizer_old = R6Class("AcqOptimizer_old",
  public = list(

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
    },

    #' @description
    #' Optimize the acquisition function.
    #'
    #' @param acq_function [AcqFunction].
    #' @param archive [bbotk::Archive].
    optimize = function(acq_function, archive) {
      stop("abstract")
    }
  )
)



#' @title Acquisition Optimizer Random Search
#'
#' @description
#' `AcqOptimizerRandomSearch_old` class that implements a random search for the
#' optimization of acquisition functions.
#'
#' @export
AcqOptimizerRandomSearch_old = R6Class("AcqOptimizerRandomSearch_old",
  inherit = AcqOptimizer_old,

  public = list(

    #' @field param_set ([paradox::ParamSet]).
    param_set = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      self$param_set = ParamSet$new(list(
        ParamInt$new("iters", lower = 1L),
        ParamLgl$new("trafo"))
      )
      self$param_set$values = list(iters = 1000L, trafo = FALSE)
    },

    #' @description
    #' Optimize the acquisition function.
    #'
    #' @param acq_function [AcqFunction].
    #' @param archive [bbotk::Archive].
    optimize = function(acq_function, archive) {
      xdt = generate_design_random(acq_function$domain, self$param_set$values$iters)
      txdt = if (self$param_set$values$trafo) {
        as.data.table(do.call(rbind, xdt$transpose()))
      } else {
        xdt$data
      }
      ydt = acq_function$eval_dt(txdt) * mult_max_to_min(acq_function$codomain)
      best = which(ydt[[1L]] == min(ydt[[1L]]))
      if (length(best) > 1L) {
        best = sample(best, 1L)
      }
      xdt$data[best, ]
    }
))




#' @title Acquisition Optimizer MIES
#'
#' @description
#' `AcqOptimizerMIES_old` class that uses [CEGO::optimMIES] for the
#' optimization of acquisition functions.
#'
#'
#' @export
AcqOptimizerMIES_old = R6Class("AcqOptimizerMIES_old",
  inherit = AcqOptimizer_old,

  public = list(

    #' @field param_set ([paradox::ParamSet]).
    param_set = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      self$param_set = ParamSet$new(list(
        ParamLgl$new("classical"),
        ParamInt$new("budget", lower = 1L),
        ParamInt$new("popsize", lower = 1L),
        ParamLgl$new("niches"))
      )
      self$param_set$values = list(classical = FALSE, budget = 1000L, popsize = 100L, niches = TRUE)
    },

    #' @description
    #' Optimize the acquisition function.
    #'
    #' @param acq_function [AcqFunction].
    #' @param archive [bbotk::Archive].
    optimize = function(acq_function, archive) {
      fun = function(x) {
        mult_max_to_min(acq_function$codomain) * unlist(acq_function$eval_dt(setNames(transpose(as.data.table(x)), nm = acq_function$domain$ids())))
      }
      # FIXME: CEGO with custom start design is buggy af
      get_sigma0 = function(ps) {
        t(as.matrix(map_dbl(ps$params, .f = function(param) {
          switch(class(param)[1L],
            "ParamDbl" = (param$upper - param$lower) * 0.1,
            "ParamInt" = (param$upper - param$lower) * (1 / 3),
            "ParamFct" = 0.1
          )
        })))
      }

      # FIXME: currently always uses best_niches + 1 random
      if (self$param_set$values$classical) {
        mies_res = CEGO::optimMIES(x = NULL, fun = fun,
          control = list(budget = self$param_set$values$budget,
            popsize = self$param_set$values$popsize, vectorized = TRUE,
            generations = Inf,
            types = acq_function$domain$storage_type,
            lower = unlist(map(acq_function$domain$params, "lower")),
            upper = unlist(map(acq_function$domain$params, "upper"))))
        setNames(transpose(as.data.table(mies_res$xbest)), nm = acq_function$domain$ids())
      } else {
        # FIXME: currently always uses best_niches + 1 random
        best_niches = if (self$param_set$values$niches) {
          acq_function$bests[, acq_function$cols_x, with = FALSE]
        } else {
          acq_function$archive_data[, acq_function$cols_x, with = FALSE]
        }
        xdt = cbind(best_niches, sigma0 = get_sigma0(acq_function$domain))
        x = as.list(transpose(xdt))

        mies_res = CEGO::optimMIES(x = x, fun = fun,
          control = list(budget = self$param_set$values$budget,
            popsize = length(x), vectorized = TRUE,
            generations = Inf,
            types = acq_function$domain$storage_type,
            lower = unlist(map(acq_function$domain$params, "lower")),
            upper = unlist(map(acq_function$domain$params, "upper"))))
        setNames(transpose(as.data.table(mies_res$xbest)), nm = acq_function$domain$ids())
      }
    }
))



#' @title Acquisition Optimizer Mutation Bananas
#'
#' @description
#' `AcqOptimizerMutateBananas` class that implements a mutation
#' algorithm for the optimization of acquisition functions.
#'
#' @export
AcqOptimizerMutateBananas = R6Class("AcqOptimizerMutateBananas",
  inherit = AcqOptimizer_old,

  public = list(

    #' @field param_set ([paradox::ParamSet]).
    param_set = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      self$param_set = ParamSet$new(list(
        ParamInt$new("n", lower = 1L),
        ParamLgl$new("trafo"))
      )
      self$param_set$values = list(n = 100, trafo = FALSE)
    },

    #' @description
    #' Optimize the acquisition function.
    #'
    #' @param acq_function [AcqFunction].
    #' @param archive [bbotk::Archive].
    optimize = function(acq_function, archive) {

      transpose = function(data, ps, filter_na = TRUE, trafo = TRUE) {
        assert_flag(filter_na)
        assert_flag(trafo)
        xs = transpose_list(data)
        if (filter_na) {
          xs = map(xs, function(x) Filter(Negate(is_scalar_na), x))
        }
        if (ps$has_trafo && trafo) {
          xs = map(xs, function(x) ps$trafo(x, ps))
        }
        return(xs)
      }

      data = archive$best()[, archive$cols_x, with = FALSE]
      xdt = map_dtr(seq_len(self$param_set$values$n), .f = function(x) mutate(data, acq_function))
      xdt = unique(xdt)
      check = apply(xdt, MARGIN = 1L, FUN = function(x) x == data)
      check[is.na(check)] = TRUE
      drop = which(colSums(check) == NCOL(data))
      xdt = xdt[-drop, ]

      if (NROW(xdt) == 0L) {
        return(generate_design_random(acq_function$domain, n = 1L)$data)
      }

      txdt = if (self$param_set$values$trafo) {
        as.data.table(do.call(rbind, transpose(xdt, acq_function$domain)))
      } else {
        xdt
      }

      ydt = acq_function$eval_dt(txdt) * mult_max_to_min(acq_function$codomain)
      best = which(ydt[[1L]] == min(ydt[[1L]]))
      if (length(best) > 1L) {
        best = sample(best, 1L)
      }
      xdt[best, ]
    }
))

mutate = function(data, acq_function) {
  # uniform mutation
  ind = sample(setdiff(seq_len(NCOL(data)), which(is.na(data))), size = 1L)
  qunif = runif(1, min = 0, max = 1)
  data[[ind]] = acq_function$domain$params[[colnames(data)[ind]]]$qunif(qunif)
  data

  # params that are NA need their default here
  data = setDT(imap(data, .f = function(value, name) {
    if (is.na(value)) {
      acq_function$domain$params[[name]]$default
    } else {
      value
    }
  }))

  for (i in seq_len(NROW(acq_function$domain$deps))) {
    dep = acq_function$domain$deps[i, ]

    if (any(map_lgl(dep[["cond"]], .f = function(cond) cond$test(data[[dep[["on"]]]])) == FALSE)) {
      data[[dep[["id"]]]] = switch(acq_function$domain$storage_type[[dep[["id"]]]], "integer" = NA_integer_, "double" = NA_real_, "character" = NA_character_)
    }
  }

  data
}
