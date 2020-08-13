#' @title Surrogate Model from mlr3 Learner
#'
#' @description
#'
#' @export
SurrogateLearner = R6Class("SurrogateLearner",
  inherit = Surrogate,
  public = list(

    learner = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param learner [mlr3::LearnerRegr]
    initialize = function(learner, columns) {
      super$initialize(columns)
      self$learner = assert_learner(learner)
      if("se" %in% self$learner$predict_types) {
        self$learner$predict_type = "se"
      }
    },

    #' @description
    #' Train model with new points.
    #'
    update = function(archive) {
      xydt = archive$xydt(y = self$columns)
      task = TaskRegr$new(id = "surrogate_task", backend = xydt, target = self$columns)
      self$learner$train(task)
    },

    #' @description
    #' Returns mean response and standard error
    #'
    #' @return named list with one [data.table::data.table]
    predict = function(xdt) {
      pred = self$learner$predict_newdata(newdata = xdt)
      if(self$learner$predict_type == "se") {
        ydt = data.table(mean = pred$response, se = pred$se)
      } else {
        ydt = data.table(mean = pred$response)
      }
      res = list(ydt)
      names(res) = self$columns
      return(res)
    }
  )
)
