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
    optimize = function(acq_function) {
      stop("abstract")
    }
  )
)



#' @title Acquisition Optimizer Random Search Bananas Trafo
#'
#' @description
#' `AcqOptimizerRandomSearch_old` class that implements a random search for the
#' optimization of acquisition functions.
#'
#' @export
AcqOptimizerRandomSearch_old_trafo = R6Class("AcqOptimizerRandomSearch_old_trafo",
  inherit = AcqOptimizer_old,

  public = list(

    #' @field param_set ([paradox::ParamSet]).
    param_set = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      self$param_set = ParamSet$new(list(
        ParamInt$new("iters", lower = 1L)
      ))
      self$param_set$values$iters = 1000L
    },

    #' @description
    #' Optimize the acquisition function.
    #'
    #' @param acq_function [AcqFunction].
    optimize = function(acq_function) {
      xdt = generate_design_random(acq_function$domain, self$param_set$values$iters)
      txdt = as.data.table(do.call(rbind, xdt$transpose()))
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
    optimize = function(acq_function) {
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



#' @title Acquisition Optimizer Mutation Crossover
#'
#' @description
#' `AcqOptimizerMutationCrossover` class that implements a mutation crossover
#' (QDO) algorithm for the optimization of acquisition functions.
#'
#' @export
AcqOptimizerMutateCrossover_old = R6Class("AcqOptimizerMutateCrossover_old",
  inherit = AcqOptimizer_old,

  public = list(

    #' @field param_set ([paradox::ParamSet]).
    param_set = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      self$param_set = ParamSet$new(list(
        ParamInt$new("iters", lower = 1L),
        ParamLgl$new("niches")
      ))
      self$param_set$values = list(iters = 1000L, niches = TRUE)
    },

    #' @description
    #' Optimize the acquisition function.
    #'
    #' @param acq_function [AcqFunction].
    optimize = function(acq_function) {
      best_niches = if (self$param_set$values$niches) {
        acq_function$bests[, acq_function$cols_x, with = FALSE]
      } else {
        acq_function$archive_data[, acq_function$cols_x, with = FALSE]
      }
      # resolve dependencies by setting up a Design
      #xdt = Design$new(acq_function$search_space,
        #map_dtr(seq_len(self$param_set$values$iters), .f = function(x) mutate_niches(best_niches, acq_function)), remove_dupl = FALSE)$data
      xdt = map_dtr(seq_len(self$param_set$values$iters), .f = function(x) mutate_niches(best_niches, acq_function))
      ydt = acq_function$eval_dt(xdt) * mult_max_to_min(acq_function$codomain)
      best = which(ydt[[1L]] == min(ydt[[1L]]))
      if (length(best) > 1L) {
        best = sample(best, 1L)
      }
      xdt[best, ]
    }
))

mutate_niches = function(best_niches, acq_function) {
  #checkmate::assert_data_table(best_niches, min.rows = 1L, min.cols = 1, null.ok = FALSE)
  number_of_niches = NROW(best_niches)

  # uniform mutation
  best_niches = setDT(imap(best_niches, .f = function(value, name) {
    mutation_prob = runif(number_of_niches, min = 0, max = 1)
    qunif = runif(number_of_niches, min = 0, max = 1)
    mutate = mutation_prob > 0.5
    value[mutate] = acq_function$domain$params[[name]]$qunif(qunif[mutate])
    # FIXME: paradox issue 318
    if (acq_function$domain$params[[name]]$storage_type == "integer") {
      as.integer(value)
    } else {
      value
    }
  }))

  if (number_of_niches > 1L) {
    # uniform crossover
    best_niches = setDT(map(best_niches, .f = function(value) {
      value[which.max(runif(number_of_niches, min = 0, max = 1))]
    }))
  }

  # FIXME: params that are NA need their default here
  best_niches = setDT(imap(best_niches, .f = function(value, name) {
    if (is.na(value)) {
      acq_function$domain$params[[name]]$default
    } else {
      value
    }
  }))

  for(i in seq_len(NROW(acq_function$domain$deps))) {
    dep = acq_function$domain$deps[i, ]

    if (any(map_lgl(dep[["cond"]], .f = function(cond) cond$test(best_niches[[dep[["on"]]]])) == FALSE)) {
      best_niches[[dep[["id"]]]] = switch(acq_function$domain$storage_type[[dep[["id"]]]], "integer" = NA_integer_, "double" = NA_real_, "character" = NA_character_)
    }
  }

  best_niches
}
