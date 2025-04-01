pkgname <- "EBMRalgorithm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "EBMRalgorithm-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('EBMRalgorithm')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("EBMRAlgorithm")
### * EBMRAlgorithm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: EBMRAlgorithm
### Title: EBMRAlgorithm Class
### Aliases: EBMRAlgorithm

### ** Examples

## Not run: 
##D # Define propensity score specifications (for illustration)
##D ps_specifications = list(
##D   formula.list = list(formula1, formula2),
##D   h_x_names.list = list(h_x1, h_x2),
##D   inv_link = inv_link_function
##D )
##D 
##D # Initialize the EBMRAlgorithm class
##D ebmr = EBMRAlgorithm$new(y_names = "outcome", ps_specifications = ps_specifications, data = dataset)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("EBMRAlgorithm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("EBMR_IPW")
### * EBMR_IPW

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: EBMR_IPW
### Title: Compute the proposed Inverse Probability Weighting (IPW)
###   Estimator
### Aliases: EBMR_IPW

### ** Examples

## Not run: 
##D ebmr <- EBMRAlgorithm$new(data = data)
##D ipw_estimates <- ebmr$EBMR_IPW(h_x_names = covariates, true_ps = true_ps_data)
##D print(ipw_estimates)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("EBMR_IPW", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("EBMR_IPW_with_locally_misspecified_model")
### * EBMR_IPW_with_locally_misspecified_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: EBMR_IPW_with_locally_misspecified_model
### Title: Compute the IPW Estimator with Locally Misspecified Model
### Aliases: EBMR_IPW_with_locally_misspecified_model

### ** Examples

## Not run: 
##D ebmr <- EBMRAlgorithm$new(data = data)
##D ipw_sensitivity <- ebmr$EBMR_IPW_with_locally_misspecified_model(
##D   ps.matrix = estimated_ps_matrix,
##D   perturb_ps = 2,
##D   exp_tilt = function(y, data) exp(y * data$covariate),
##D   exp_tilt_x_names = c("covariate")
##D )
##D print(ipw_sensitivity)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("EBMR_IPW_with_locally_misspecified_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("WangShaoKim2014")
### * WangShaoKim2014

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: WangShaoKim2014
### Title: Wang, Shao, and Kim (2014) Method for Propensity Score
###   Estimation
### Aliases: WangShaoKim2014

### ** Examples

## Not run: 
##D data <- data.frame(x = rnorm(100), y = rnorm(100), r = rbinom(100, 1, 0.5))
##D ebmr <- EBMRAlgorithm$new(data = data)
##D result <- ebmr$WangShaoKim2014(
##D   formula = r ~ x,
##D   h_x_names = c("x"),
##D   inv_link = function(eta) 1 / (1 + exp(-eta))
##D )
##D print(result$coefficients)
## End(Not run)

## Not run: 
##D data <- data.frame(x = rnorm(100), y = rnorm(100), r = rbinom(100, 1, 0.5))
##D ebmr <- EBMRAlgorithm$new(data = data)
##D result <- ebmr$WangShaoKim2014(
##D   formula = r ~ x,
##D   h_x_names = c("x"),
##D   inv_link = function(eta) 1 / (1 + exp(-eta))
##D )
##D print(result$coefficients)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("WangShaoKim2014", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("check_data")
### * check_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: check_data
### Title: Check Data for Missing Response Values
### Aliases: check_data
### Keywords: internal

### ** Examples

## Not run: 
##D data = data.frame(outcome = c(1, NA, 3, 4, NA), covariate = c(2, 3, 4, 5, 6))
##D y_names = "outcome"
##D cleaned_data = ebmr_algorithm$private$check_data(y_names, data)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("check_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ensemble")
### * ensemble

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ensemble
### Title: Estimate the Coefficients for the Propensity Score Model
### Aliases: ensemble
### Keywords: internal

### ** Examples

## Not run: 
##D ebmr <- EBMRAlgorithm$new(data = data)
##D nu_estimates <- ebmr$estimate_nu(ps.matrix = ps_data, h_x = covariates)
##D print(nu_estimates)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ensemble", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("parse_formula")
### * parse_formula

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: parse_formula
### Title: Parse a Formula into Response Indicator and Predictor Components
### Aliases: parse_formula
### Keywords: internal

### ** Examples

## Not run: 
##D ebmr <- EBMRAlgorithm$new(data = data)
##D parsed_formula <- ebmr$parse_formula(r ~ x1 + o(y) + x2)
##D print(parsed_formula)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("parse_formula", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("separate_variable_types")
### * separate_variable_types

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: separate_variable_types
### Title: Separate Continuous and Discrete Variables from a Data Frame
### Aliases: separate_variable_types
### Keywords: internal

### ** Examples

## Not run: 
##D data <- data.frame(
##D   age = c(25, 30, 35, 40),
##D   gender = c('Male', 'Female', 'Female', 'Male'),
##D   salary = c(50000, 60000, 55000, 70000)
##D )
##D ebmr <- EBMRAlgorithm$new(data = data)
##D result <- ebmr$separate_variable_types(data)
##D print(result)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("separate_variable_types", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
