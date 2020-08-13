MboArchive = R6Class("MboArchive",
  inherit = Archive,
  public = list(
    
    temp_cols = NULL,

    # allows to set new columns, if new column has same name as existing columns in the archive, the existing columns will remain but will be ignored for the sake of the new column
    set_column = function(name, vals) {
      assert_numeric(vals, len = self$n_evals)
      assert_string(name) #FIXME, check for valid column name
      if (is.null(self$temp_cols) || nrow(self$temp_cols) != length(vals)) { #new data.table
        self$temp_cols = data.table(vals)
        names(self$temp_cols) = name
      } else { #add column
        self$temp_cols[, (name) := vals]
      }
    },

    set_y_column = function()

    xydt = function(x = self$cols_x, y = self$cols_y) {
      if (!is.null(self$temp_cols) && any(y %in% colnames(self$temp_cols))) {
        # first we respect columns that are in temp_cols
        y_temp_cols = intersect(y, colnames(self$temp_cols))
        # the rest will be taken from self$data()
        # effect: if columns exist in temp_cols AND data() we will take the ones in temp_cols
        if (length(y_temp_cols) > 0 && nrow(self$temp_cols) != self$n_evals) {
          stop("Data in MboArchive$temp_cols is outdated!")
        }
        y_data = setdiff(y, y_temp_cols)
        cbind(self$data()[, c(x, y_data), with = FALSE], self$temp_cols[, y_temp_cols, with = FALSE])  
      } else {
        self$data()[, c(x, y), with = FALSE]
      } 
    },

    ydt = function(y = self$cols_y) {
      self$xydt(x = character(), y = y)
    }
  )
)