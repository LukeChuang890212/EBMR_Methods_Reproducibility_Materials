#' Parse a Formula and Identify Outcome-Related Terms
#'
#' This function parses a model formula and automatically identifies which terms
#' involve the outcome variable (which may be subject to missingness).
#'
#' @param formula A formula specifying the propensity score model.
#' Examples:
#' \itemize{
#'   \item Simple MNAR: `r ~ y + health + father`
#'   \item MNAR with interactions: `r ~ y + y:health + y:father + health + father`
#'   \item MAR model: `r ~ health + father`
#' }
#' @param outcome Character string specifying the outcome variable name (e.g., "y").
#' Terms containing this variable are treated as MNAR terms.
#' If NULL, the model is treated as MAR.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{r_name}}{The name of the response indicator variable.}
#'   \item{\code{outcome}}{The outcome variable name.}
#'   \item{\code{y_terms}}{Character vector of terms involving the outcome.}
#'   \item{\code{x_terms}}{Character vector of terms not involving the outcome.}
#'   \item{\code{all_terms}}{All terms in the formula (excluding intercept).}
#'   \item{\code{is_mnar}}{Logical, TRUE if any outcome-related terms exist.}
#' }
#'
#' @details
#' The function uses R's native formula parsing via `terms()` to properly handle
#' interactions (`:`, `*`) and other formula operators. It then identifies which
#' terms contain the outcome variable.
#'
#' @examples
#' \dontrun{
#' # MNAR model with interactions
#' result <- parse_formula(
#'   r ~ y + y:health + y:father + health + father,
#'   outcome = "y"
#' )
#' # result$y_terms: c("y", "y:health", "y:father")
#' # result$x_terms: c("health", "father")
#'
#' # MAR model
#' result <- parse_formula(r ~ health + father, outcome = NULL)
#' # result$y_terms: character(0)
#' # result$x_terms: c("health", "father")
#' }
#'
#' @keywords internal

parse_formula <- function(formula, outcome = NULL) {

  # Get the response variable name (left-hand side)
  r_name <- as.character(formula[[2]])

  # Use terms() to properly parse the formula
  term_obj <- terms(formula)
  term_labels <- attr(term_obj, "term.labels")

  # If no outcome specified, treat all terms as x_terms (MAR model)
  if (is.null(outcome)) {
    return(list(
      r_name = r_name,
      outcome = NULL,
      y_terms = character(0),
      x_terms = term_labels,
      all_terms = term_labels,
      is_mnar = FALSE
    ))
  }

  # Identify which terms contain the outcome variable
  # A term contains outcome if the outcome variable name appears as a component
  # For "y:health", the components are "y" and "health"
  # For "y", the component is just "y"

  # Function to check if a term involves the outcome
  involves_outcome <- function(term, outcome_var) {
    # Split term by : to get components
    components <- strsplit(term, ":")[[1]]
    outcome_var %in% components
  }

  y_mask <- sapply(term_labels, involves_outcome, outcome_var = outcome)
  y_terms <- term_labels[y_mask]
  x_terms <- term_labels[!y_mask]

  return(list(
    r_name = r_name,
    outcome = outcome,
    y_terms = y_terms,
    x_terms = x_terms,
    all_terms = term_labels,
    is_mnar = length(y_terms) > 0
  ))
}


#' Build Design Matrix from Parsed Formula
#'
#' Constructs the design matrix for the propensity score model using the
#' parsed formula components.
#'
#' @param data A data frame containing all variables.
#' @param parsed_formula Output from `parse_formula()`.
#' @param include_intercept Logical, whether to include an intercept column.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{design}}{The full design matrix.}
#'   \item{\code{y_design}}{Design matrix columns for outcome-related terms.}
#'   \item{\code{x_design}}{Design matrix columns for non-outcome terms.}
#'   \item{\code{n_y}}{Number of outcome-related columns.}
#'   \item{\code{n_x}}{Number of non-outcome columns (excluding intercept).}
#' }
#'
#' @keywords internal

build_design_matrix <- function(data, parsed_formula, include_intercept = TRUE) {
  n <- nrow(data)

  # Build y_design (outcome-related terms)
  if (length(parsed_formula$y_terms) > 0) {
    y_formula <- as.formula(
      paste("~", paste(parsed_formula$y_terms, collapse = " + "), "- 1")
    )
    y_design <- model.matrix(y_formula, data = data)
  } else {
    y_design <- matrix(nrow = n, ncol = 0)
  }

  # Build x_design (non-outcome terms)
  if (length(parsed_formula$x_terms) > 0) {
    x_formula <- as.formula(
      paste("~", paste(parsed_formula$x_terms, collapse = " + "), "- 1")
    )
    x_design <- model.matrix(x_formula, data = data)
  } else {
    x_design <- matrix(nrow = n, ncol = 0)
  }

  # Combine: intercept (optional), y_design, x_design
  if (include_intercept) {
    design <- cbind(1, y_design, x_design)
    colnames(design)[1] <- "(Intercept)"
  } else {
    design <- cbind(y_design, x_design)
  }

  return(list(
    design = design,
    y_design = y_design,
    x_design = x_design,
    n_y = ncol(y_design),
    n_x = ncol(x_design)
  ))
}


#' Separate Continuous and Discrete Variables from a Data Frame
#'
#' This method separates continuous and discrete variables from a given data
#' frame based on their characteristics.
#'
#' @param x A data frame containing the variables to be separated.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{x1}}{A data frame with the discrete variables.}
#'   \item{\code{x2}}{A data frame with the continuous variables.}
#'   \item{\code{x1_names}}{Names of the discrete variables.}
#'   \item{\code{x2_names}}{Names of the continuous variables.}
#' }
#'
#' @details
#' Continuous variables: numeric with more than 10 unique values.
#' Discrete variables: non-numeric or numeric with 10 or fewer unique values.
#'
#' @keywords internal

separate_variable_types <- function(x) {
  x1 <- list()
  x2 <- list()
  col_names <- colnames(x)
  x1_names <- c()
  x2_names <- c()

  for (i in seq_len(ncol(x))) {
    column <- x[, i]
    if (is.numeric(column)) {
      if (length(unique(column)) > 10) {
        x2[[col_names[i]]] <- column
        x2_names <- c(x2_names, col_names[i])
      } else {
        x1[[col_names[i]]] <- column
        x1_names <- c(x1_names, col_names[i])
      }
    } else {
      x1[[col_names[i]]] <- column
      x1_names <- c(x1_names, col_names[i])
    }
  }

  x1 <- as.data.frame(x1)
  x2 <- as.data.frame(x2)

  return(list(x1 = x1, x2 = x2, x1_names = x1_names, x2_names = x2_names))
}
