% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Preprocessor.r
\name{parse_formula}
\alias{parse_formula}
\title{Parse a Formula into Response Indicator and Predictor Components}
\usage{
parse_formula(formula)
}
\arguments{
\item{formula}{A formula specifying the model structure (e.g., \code{r ~ x1 + o(y) + x2}),
where the missing outcome variables should be placed inside \code{o()}.}
}
\value{
A list containing:
\describe{
\item{\code{r_names}}{The name of the response indicator variable.}
\item{\code{y_names}}{The variable(s) enclosed in \code{o()}, representing missing outcome variables, if present.}
\item{\code{x_names}}{The names of the remaining predictor variables.}
}
}
\description{
This method extracts key components from a model formula, including:
the response indicator variable (\code{r}), the variable inside \code{o()},
and the remaining predictor variables.
}
\details{
This function processes the formula by extracting:
\itemize{
\item The response indicator variable on the left-hand side of \code{~}.
\item The missing outcome variable(s) inside \code{o()}, if any.
\item The remaining predictors, excluding variables in \code{o()}.
}
}
\examples{
\dontrun{
ebmr <- EBMRAlgorithm$new(data = data)
parsed_formula <- ebmr$parse_formula(r ~ x1 + o(y) + x2)
print(parsed_formula)
}

}
\keyword{internal}
