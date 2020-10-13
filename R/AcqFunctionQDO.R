#' @title Acquisition Function QDO
#'
#' @description
#' Based on a surrogate model, the acquisition function encodes the preference to evaluate
#' a new point for evaluation. QDO.
#'
#' @export
AcqFunctionQDO = R6Class("AcqFunctionQDO",
  public = list(

    #' @field id (`character(1)`).
    id = NULL,

    #' @field surrogate [Surrogate].
    surrogate = NULL,

    #' @field param_set ([paradox::ParamSet]).
    param_set = NULL,

    #' @field search_space ([paradox::ParamSet]).
    search_space = NULL,

    #' @field codomain ([paradox::ParamSet]).
    codomain = NULL,

    #' @field direction (`character(1)`).
    direction = NULL, # optim direction of the acq function

    #' @field surrogate_max_to_min (`numeric(1)`).
    surrogate_max_to_min = NULL, # optim direction of the obj function 1 for min, -1 for max, maybe it makes sense to make this private so it is clear that this is not meant to turn the acq into a minimization problem

    #' @field niches.
    #' FIXME:
    niches = NULL, 

    #' @field feature_surrogate_pointer 
    #' FIXME: actually having a pointer here to the feature surrogate would be perfect
    feature_surrogate_predict = NULL,

    #' @field niche_boundaries
    niche_boundaries = NULL,
 
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param id (`character(1)`).
    #' @param param_set ([paradox::ParamSet]).
    #' @param surrogate [Surrogate].
    #' @param direction (`character(1)`).
    initialize = function(id, param_set, surrogate, direction) {
      self$id = assert_string(id)
      self$param_set = assert_param_set(param_set)
      self$surrogate = assert_r6(surrogate, "Surrogate")
      self$direction = assert_choice(direction, c("same", "minimize", "maximize"))
    },

    #' @description
    #' Evaluates all input values in `xdt`.
    #'
    #' @param xdt [data.table::data.table]
    #'
    #' @return `data.table` \cr
    #' The column has to have the same name as the id of the acq_fun, because we
    #' renamed the id of the codomain
    eval_dt = function(xdt) {
      stop("abstract")
    },

    #' @description
    #' Sets up the acquisition function
    #'
    #' @param archive [bbotk::Archive].
    setup = function(archive, niche_boundaries, feature_surrogate_predict) {
      # FIXME: Should we allow alternative search_space as additional argument?

      # here we can change the optim direction of the codomain for the acq function
      self$codomain = generate_acq_codomain(archive$codomain, id = self$id, direction = self$direction)

      self$surrogate_max_to_min = mult_max_to_min(archive$codomain)

      self$search_space = archive$search_space

      self$niches = archive$niches

      self$niche_boundaries = niche_boundaries

      self$feature_surrogate_predict = feature_surrogate_predict
    },

    #' @description
    #' Update the acquisition function
    #'
    #' @param archive [bbotk::Archive].
    update = function(archive) {
      # it's okay to do nothing here
    },

    #' @description
    #' Generates an objective function for [bbotk::Optimizer].
    #'
    #' @return [bbotk::ObjectiveRFunDt]
    generate_objective = function() {
      bbotk::ObjectiveRFunDt$new(
        fun = self$eval_dt,
        domain = self$search_space,
        codomain = self$codomain,
        id = self$id
      )
    }
  )
)
