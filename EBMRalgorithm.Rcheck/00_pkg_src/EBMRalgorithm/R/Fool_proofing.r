#' Check Data for Missing Response Values
#'
#' This method checks if the outcome variable has missing values and creates a response indicator column (`r`) if necessary.
#' If the outcome variable does not have missing values, it checks if the response indicator column (`r`) already exists.
#'
#' @param y_names A vector of column names for the outcome variable(s).
#' @param data A data frame containing the data.
#' @return A data frame with response indicator column `r`.
#' @examples
#' \dontrun{
#' data = data.frame(outcome = c(1, NA, 3, 4, NA), covariate = c(2, 3, 4, 5, 6))
#' y_names = "outcome"
#' cleaned_data = ebmr_algorithm$private$check_data(y_names, data)
#' }
#'
#' @keywords internal

check_data = function(y_names, data) {
  if(sum(is.na(data[y_names])) > 0){
    data$r = ifelse(is.na(data[y_names]), 0, 1)
  }else{
    if(!("r" %in% colnames(data))){
      stop("there must be a column of response indicator with column name \"r\".")
    }
  }
  return(data)
}

