#' @title Single Point Proposal Generator
#' @export
#'
ProposalGeneratorSingle = R6Class("ProposalGeneratorSingle",
  inherit = ProposalGenerator,

  public = list(

    acq_function = NULL,
    acq_optimizer = NULL,

    initialize = function(acq_function, acq_optimizer) {
      self$acq_function = assert_r6(acq_function, "AcqFunction")
      self$acq_optimizer = assert_r6(acq_optimizer, "AcqOptimizer")
    },

    #' @return data.table \cr
    #'   data.table with columns of domain$ids() and possible extras
    propose = function() {
      self$acq_optimizer$optimize(self$acq_function)
    },

    setup = function(archive) {
      super$setup(archive)
      self$acq_function$setup(archive)
      #self$acq_optimizer$setup()
    },

    update = function() {
      self$acq_function$update(self$archive)
      #self$acq_optimizer$update()
    }


  )
)
