#' @title Surrogate Model for MultiCriteria response surfaces
#'
#' @description
#' Multi Criteria response surfaces modeled by multiple regression [mlr3::Learner] objects.
#'
#' @export
SurrogateMultiCritLearners = R6Class("SurrogateMultiCritLearners",
  inherit = SurrogateMultiCrit,
  public = list(

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param learners (list of [mlr3::LearnerRegr]).
    initialize = function(learners) {
      self$model = assert_learners(learners)
      for (model in self$model) {
        if (model$predict_type != "se" && "se" %in% model$predict_types) {
          model$predict_type = "se"
        }
      }

      ps = ParamSet$new(list(
        ParamLgl$new("calc_insample_perf"),
        ParamUty$new("perf_measures", custom_check = function(x) check_list(x, types = "MeasureRegr", any.missing = FALSE, len = length(self$model))),  # FIXME: actually want check_measures
        ParamUty$new("perf_thresholds", custom_check = function(x) check_double(x, lower = -Inf, upper = Inf, any.missing = FALSE, len = length(self$model))))
      )
      ps$values = list(calc_insample_perf = FALSE, perf_measures = replicate(length(self$model), msr("regr.rsq")), perf_thresholds = rep(0, length(self$model)))
      ps$add_dep("perf_measures", on = "calc_insample_perf", cond = CondEqual$new(TRUE))
      ps$add_dep("perf_thresholds", on = "calc_insample_perf", cond = CondEqual$new(TRUE))
      private$.param_set = ps
    },

    #' @description
    #' Returns mean response and standard error
    #'
    #' @param xdt [data.table::data.table()]\cr
    #' New data.
    #'
    #' @return [data.table::data.table()] with the columns `mean` and `se`.
    predict = function(xdt) {
      assert_xdt(xdt)

      preds = lapply(self$model, function(model) {
        pred = model$predict_newdata(newdata = xdt)
        if (model$predict_type == "se") {
          data.table(mean = pred$response, se = pred$se)
        } else {
          data.table(mean = pred$response)
        }
      })
      names(preds) = names(self$model)
      return(preds)
    }
  ),

  active = list(

    #' @field k (`integer(1)`)\cr
    #' Returns the number of models.
    k = function() {
      length(self$model)
    },

    #' @field assert_insample_perf (`numeric()`) \cr
    #' Asserts whether the current insample performance meets the performance threshold.
    assert_insample_perf = function(rhs) {
      if (!missing(rhs)) {
        stopf("Field/Binding is read-only.")
      }

      if (!self$param_set$values$calc_insample_perf) {
        return(invisible(self$insample_perf))
      }

      check = all(pmap_lgl(
        list(
          insample_perf = self$insample_perf,
          perf_threshold = self$param_set$values$perf_thresholds,
          perf_measure = self$param_set$values$perf_measures
        ),
        .f = function(insample_perf, perf_threshold, perf_measure) {
          if (perf_measure$minimize) {
            insample_perf < perf_threshold
          } else {
            insample_perf > perf_threshold
          }
        })
      )

      if (!check) {
        stopf("Current insample performance of the Surrogate Model does not meet the performance threshold")
      }
      invisible(self$insample_perf)
    }
  ),

  private = list(

    # Train model with new points.
    # Also calculates the insample performance based on the `perf_measures` hyperparameter if `calc_insample_perf = TRUE`.
    .update = function(xydt, y_cols) {
      assert_xydt(xydt, y_cols)

      backend = as_data_backend(xydt)
      features = setdiff(names(xydt), y_cols)

      tasks = lapply(y_cols, function(y_col) {
        # If this turns out to be a bottleneck, we can also operate on a
        # single task here
        task = TaskRegr$new(
          id = paste0("surrogate_task_", y_col),
          backend = backend,
          target = y_col)
        task$col_roles$feature = features
        task
      })
      pmap(list(model = self$model, task = tasks), .f = function(model, task) {
        model$train(task)
        NULL
      })
      names(self$model) = y_cols

      if (self$param_set$values$calc_insample_perf) {
        private$.insample_perf = setNames(pmap_dbl(list(model = self$model, task = tasks, perf_measure = self$param_set$values$perf_measures),
          .f = function(model, task, perf_measure) {
            assert_measure(perf_measure, task = task, learner = model)
            model$predict(task)$score(perf_measure, task = task, learner = model)
          }
        ), nm = map_chr(self$param_set$values$perf_measures, "id"))
        self$assert_insample_perf
      }
    }
  )
)
