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



#' @title Acquisition Optimizer Random Search
#'
#' @description
#' `AcqOptimizerRandomSearch` class that implements a random search for the
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
        ParamInt$new("iters", lower = 1L, default = 1000L)
      ))
      self$param_set$values$iters = 1000L
    },

    #' @description
    #' Optimize the acquisition function.
    #'
    #' @param acq_function [AcqFunction].
    optimize = function(acq_function) {
      xdt = generate_design_random(acq_function$domain, self$param_set$values$iters)$data
      ydt = acq_function$eval_dt(char_to_fct(xdt)) * mult_max_to_min(acq_function$codomain)
      best = which(ydt[[1]] == min(ydt[[1]]))
      if (length(best) > 1) {
        best = sample(best, 1)
      }
      xdt[best, ]
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

  # params that are NA need their default here
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
