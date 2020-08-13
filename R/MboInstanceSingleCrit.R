#' @title Single Criteria MBO Instance
#'

#' @export
MboInstanceSingleCrit = R6Class("MboInstanceSingleCrit",
  inherit = OptimInstanceSingleCrit,
  public = list(
    initialize = function(objective, search_space = NULL, terminator) {
      super$initialize(objective, search_space = NULL, terminator)
      self$archive = MboArchive$new(search_space = self$search_space,
        codomain = objective$codomain)
    }
  )
)
