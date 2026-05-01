# 測試更robust的formula解析方式
library(stringr)

# 想要支援的formula格式
test_formulas <- list(
  as.formula("r ~ o(y) + health"),
  as.formula("r ~ o(y) + o(y:health) + health + father"),
  as.formula("r ~ o(y) + o(y*health) + health + father"),
  as.formula("r ~ o(y) + o(y:health) + o(y:father) + health + father + health:father"),
  as.formula("r ~ health + father")
)

parse_formula_v2 <- function(formula) {
  formula_str <- as.character(formula)
  r_names <- formula_str[2]
  rhs <- formula_str[3]

  # 提取所有o()內容
  o_matches <- str_extract_all(rhs, "o[(]([^)]+)[)]")[[1]]

  if (length(o_matches) == 0) {
    y_terms <- character(0)
  } else {
    y_terms <- gsub("o[(]|[)]", "", o_matches)
    y_terms <- trimws(y_terms)
  }

  # 移除o(...)得到x部分 - use fixed patterns to avoid escaping issues
  x_rhs <- rhs
  for (m in o_matches) {
    x_rhs <- sub(m, "", x_rhs, fixed = TRUE)
  }
  # Clean up: remove leading/trailing + and multiple +
  x_rhs <- gsub("^[[:space:]]*[+][[:space:]]*", "", x_rhs)
  x_rhs <- gsub("[[:space:]]*[+][[:space:]]*$", "", x_rhs)
  x_rhs <- gsub("[[:space:]]*[+][[:space:]]*[+][[:space:]]*", " + ", x_rhs)
  x_rhs <- trimws(x_rhs)

  if (x_rhs == "" || is.na(x_rhs)) {
    x_formula <- NULL
  } else {
    x_formula <- as.formula(paste("~", x_rhs))
  }

  if (length(y_terms) > 0) {
    y_formula <- as.formula(paste("~", paste(y_terms, collapse = " + "), "- 1"))
  } else {
    y_formula <- NULL
  }

  return(list(
    r_names = r_names,
    y_terms = y_terms,
    y_formula = y_formula,
    x_formula = x_formula
  ))
}

cat("Testing parse_formula_v2:\n\n")

for (f in test_formulas) {
  result <- parse_formula_v2(f)

  cat("Formula:", deparse(f), "\n")
  cat("  r_names:", result$r_names, "\n")
  cat("  y_terms:", if(length(result$y_terms)==0) "(none)" else paste(result$y_terms, collapse="; "), "\n")
  cat("  y_formula:", if(is.null(result$y_formula)) "(none)" else deparse(result$y_formula), "\n")
  cat("  x_formula:", if(is.null(result$x_formula)) "(none)" else deparse(result$x_formula), "\n")
  cat("\n")
}

# 測試用model.matrix建構design matrix
cat("\n========================================\n")
cat("Testing model.matrix construction:\n")
cat("========================================\n\n")

set.seed(123)
dat <- data.frame(
  y = rbinom(100, 1, 0.3),
  health = rnorm(100),
  father = rbinom(100, 1, 0.5),
  r = rbinom(100, 1, 0.7)
)

f <- as.formula("r ~ o(y) + o(y:health) + o(y:father) + health + father + health:father")
result <- parse_formula_v2(f)

cat("Formula:", deparse(f), "\n\n")

# 建構 y design matrix (o() terms)
if (!is.null(result$y_formula)) {
  y_design <- model.matrix(result$y_formula, data = dat)
  cat("y_design (from o() terms):\n")
  cat("  Columns:", paste(colnames(y_design), collapse = ", "), "\n")
  cat("  Dim:", nrow(y_design), "x", ncol(y_design), "\n\n")
}

# 建構 x design matrix
if (!is.null(result$x_formula)) {
  x_design <- model.matrix(result$x_formula, data = dat)
  cat("x_design (non-o() terms):\n")
  cat("  Columns:", paste(colnames(x_design), collapse = ", "), "\n")
  cat("  Dim:", nrow(x_design), "x", ncol(x_design), "\n\n")
}

# 完整的 design matrix: cbind(1, y_design, x_design without intercept)
full_design <- cbind(1, y_design, x_design[, -1, drop=FALSE])
cat("Full design matrix:\n")
cat("  Columns:", paste(colnames(full_design), collapse = ", "), "\n")
cat("  Dim:", nrow(full_design), "x", ncol(full_design), "\n")
