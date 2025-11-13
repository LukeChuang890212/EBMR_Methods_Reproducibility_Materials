#' Parse a Formula into Response Indicator and Predictor Components
#'
#' This method extracts key components from a model formula, including:
#' the response indicator variable (`r`), the variable inside `o()`,
#' and the remaining predictor variables.
#'
#' @import stringr
#'
#' @param formula A formula specifying the model structure (e.g., `r ~ x1 + o(y) + x2`),
#' where the missing outcome variables should be placed inside `o()`.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{r_names}}{The name of the response indicator variable.}
#'   \item{\code{y_names}}{The variable(s) enclosed in `o()`, representing missing outcome variables, if present.}
#'   \item{\code{x_names}}{The names of the remaining predictor variables.}
#' }
#'
#' @details
#' This function processes the formula by extracting:
#' \itemize{
#'   \item The response indicator variable on the left-hand side of `~`.
#'   \item The missing outcome variable(s) inside `o()`, if any.
#'   \item The remaining predictors, excluding variables in `o()`.
#' }
#'
#' @examples
#' \dontrun{
#' ebmr <- EBMRAlgorithm$new(data = data)
#' parsed_formula <- ebmr$parse_formula(r ~ x1 + o(y) + x2)
#' print(parsed_formula)
#' }
#'
#' @keywords internal

parse_formula = function(formula) {
  # Convert the formula to a string
  formula <- as.character(formula)

  # Extract the left-hand side of the formula (everything before the ~)
  lhs <- formula[2]

  # Isolate r
  r_names <- lhs

  # Extract the right-hand side of the formula (everything after the ~)
  rhs <- formula[3]

  # Extract the variables inside o() using regular expressions
  y_names <- str_extract(rhs, "o\\(([^)]+)\\)")   # Variables in o()

  if(is.na(y_names)){
    y_names = character(0)
  }else{
    # Clean up: Remove 'o()' and split the variables (there's only one in this case)
    y_names <- gsub("o\\(|\\)", "", y_names)
    y_names <- unlist(strsplit(y_names, split = "\\+"))
  }


  # Extract the rest of the variables (excluding o())
  x_names <- gsub("o\\([^)]+\\)\\s*\\+?", "", rhs)

  # Clean up: Split the rest of the variables by '+' and trim whitespace
  x_names <- unlist(strsplit(x_names, split = "\\+"))
  x_names <- trimws(x_names)

  # Return a list with r_names, y_names, x_names
  return(list(r_names = r_names, y_names = y_names, x_names = x_names))
}

#' Separate Continuous and Discrete Variables from a Data Frame
#'
#' This method separates continuous and discrete variables from a given data frame
#' based on their characteristics. Continuous variables are numeric columns with more than 10 unique values,
#' while discrete variables are either non-numeric or numeric columns with fewer than or equal to 10 unique values.
#'
#' @param x A data frame containing the variables to be separated.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{x1}}{A data frame with the discrete variables.}
#'   \item{\code{x2}}{A data frame with the continuous variables.}
#'   \item{\code{x1_names}}{A character vector with the names of the discrete variables.}
#'   \item{\code{x2_names}}{A character vector with the names of the continuous variables.}
#' }
#'
#' @details
#' The function separates variables into two categories:
#' \itemize{
#'   \item Continuous variables: Numeric variables with more than 10 unique values.
#'   \item Discrete variables: Non-numeric variables or numeric variables with 10 or fewer unique values.
#' }
#' The function returns the separated variables in two data frames: one for continuous variables (`x2`)
#' and one for discrete variables (`x1`), along with the corresponding names of these variables.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   age = c(25, 30, 35, 40),
#'   gender = c('Male', 'Female', 'Female', 'Male'),
#'   salary = c(50000, 60000, 55000, 70000)
#' )
#' ebmr <- EBMRAlgorithm$new(data = data)
#' result <- ebmr$separate_variable_types(data)
#' print(result)
#' }
#'
#' @keywords internal

separate_variable_types = function(x) {
  # Initialize empty lists to store continuous and discrete columns
  x1 <- list()
  x2 <- list()

  col_names = colnames(x)
  x1_names = c()
  x2_names = c()

  # Iterate through each column in the matrix 'x'
  for (i in 1:ncol(x)) {
    column <- x[, i]

    # Check if the column is numeric
    if (is.numeric(column)) {
      # If the column has more than 10 unique values, it's treated as continuous
      if (length(unique(column)) > 10) {
        x2[[colnames(x)[i]]] <- column
        x2_names = c(x2_names, col_names[i])
      } else {
        x1[[colnames(x)[i]]] <- column
        x1_names = c(x1_names, col_names[i])
      }
    } else {
      # For non-numeric columns, treat them as discrete
      x1[[colnames(x)[i]]] <- column
      x1_names = c(x1_names, col_names[i])
    }
  }

  # Convert the lists to matrices (if you want matrices, otherwise keep them as lists)
  x1 <- as.data.frame(x1)
  x2 <- as.data.frame(x2)

  # Return the results
  return(list(x1 = x1, x2 = x2, x1_names = x1_names, x2_names = x2_names))
}
