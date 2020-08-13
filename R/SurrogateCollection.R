#' @title Surrogate Model from mlr3 Learner
#'
#' @description
#'
#' @export
SurrogateCollection = R6Class("SurrogateCollection",
  inherit = Surrogate,
  public = list(

    #' @field Stores surrogate model
    surrogates = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param model Model
    initialize = function(surrogates) {
      self$surrogates = assert_list(surrogates, "Surrogate")
      super$initialize(map_chr(surrogates, "columns"))
    },

    #' @description
    #' Train model with new points.
    #'
    #' @return `NULL`
    update = function(archive) {
      for (surrogate in self$surrogates) {
        surrogate$update(archive)
      }
    },

    #' @description
    #' Possible setup routine of the surrogate
    #'
    #' @return `NULL`
    setup = function(archive) {
      for (surrogate in self$surrogates) {
        surrogate$setup(archive)
      }
    },


    #' @description
    #' Returns mean response and standard error
    #'
    #' @return list of [data.table::data.table] objects
    predict = function(xdt) {
      preds = lapply(self$surrogates, function(surrogate) surrogate$predict(xdt))
      unlist(preds, recursive = FALSE)
    }
  )
)

