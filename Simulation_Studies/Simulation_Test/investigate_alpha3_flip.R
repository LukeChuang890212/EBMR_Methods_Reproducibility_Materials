setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EBMRalgorithmFast4")
})
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
dat <- all_data[1:n_val, ]  # Rep 1

# Model 3: r ~ y + u2 + z2, h_alpha = (u1, u2, z1, z2)
formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]
cat(sprintf("Formula: %s\n", deparse(formula_j)))
cat(sprintf("h_alpha: %s\n", paste(h_alpha_vars, collapse=", ")))

r_vec <- dat[["r"]]
design_mat <- model.matrix(formula_j, data=dat)
h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
h_full <- cbind(1, h_alpha_mat)
n <- n_val

cat(sprintf("design_mat cols: %s\n", paste(colnames(design_mat), collapse=", ")))
cat(sprintf("h_full cols: intercept, %s\n", paste(h_alpha_vars, collapse=", ")))
cat(sprintf("param_dim=%d, esteq_dim=%d\n\n", ncol(design_mat), ncol(h_full)))

model_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  1 / (1 + exp(-eta))
}

Phi_alpha <- function(alpha) {
  pi_hat <- model_fn(alpha)
  (r_vec / pi_hat - 1) * h_full
}

gmm_obj <- function(alpha, W_mat=NULL) {
  g_mat <- Phi_alpha(alpha)
  G <- colMeans(g_mat)
  if (is.null(W_mat)) W_mat <- solve(crossprod(g_mat)/n)
  as.numeric(t(G) %*% W_mat %*% G)
}

# GLM init
glm_fit <- glm(formula_j, data=dat, family=binomial(link="logit"))
alpha_glm <- coef(glm_fit)
cat(sprintf("GLM init: (%s)\n", paste(round(alpha_glm, 4), collapse=", ")))

# GMM converged values from previous test
alpha_pos <- c(2.0425, -3.7987, 1.0111, 0.3857)  # from max_outer=200, rep 1
alpha_neg <- -alpha_pos

cat(sprintf("alpha_pos: (%s)\n", paste(round(alpha_pos, 4), collapse=", ")))
cat(sprintf("alpha_neg: (%s)\n", paste(round(alpha_neg, 4), collapse=", ")))

# Check: are pi(alpha) and pi(-alpha) related?
pi_pos <- model_fn(alpha_pos)
pi_neg <- model_fn(alpha_neg)
cat(sprintf("\npi(alpha_pos): mean=%.4f, range=[%.4f, %.4f]\n",
    mean(pi_pos), min(pi_pos), max(pi_pos)))
cat(sprintf("pi(alpha_neg): mean=%.4f, range=[%.4f, %.4f]\n",
    mean(pi_neg), min(pi_neg), max(pi_neg)))
cat(sprintf("pi_pos + pi_neg: mean=%.4f (should be 1 if symmetric)\n", mean(pi_pos + pi_neg)))
cat(sprintf("max|pi_pos + pi_neg - 1| = %.6f\n", max(abs(pi_pos + pi_neg - 1))))

# GMM objectives at both solutions
cat("\n=== GMM objectives ===\n")
# With identity W
cat(sprintf("Obj(alpha_pos, W=I): %.8f\n", gmm_obj(alpha_pos, diag(ncol(h_full)))))
cat(sprintf("Obj(alpha_neg, W=I): %.8f\n", gmm_obj(alpha_neg, diag(ncol(h_full)))))
cat(sprintf("Obj(alpha_glm, W=I): %.8f\n", gmm_obj(alpha_glm, diag(ncol(h_full)))))

# With optimal W at each
g_pos <- Phi_alpha(alpha_pos); W_pos <- solve(crossprod(g_pos)/n)
g_neg <- Phi_alpha(alpha_neg); W_neg <- solve(crossprod(g_neg)/n)
g_glm <- Phi_alpha(alpha_glm); W_glm <- solve(crossprod(g_glm)/n)

cat(sprintf("\nObj(alpha_pos, W(alpha_pos)): %.8f\n", gmm_obj(alpha_pos, W_pos)))
cat(sprintf("Obj(alpha_neg, W(alpha_neg)): %.8f\n", gmm_obj(alpha_neg, W_neg)))
cat(sprintf("Obj(alpha_glm, W(alpha_glm)): %.8f\n", gmm_obj(alpha_glm, W_glm)))

# Cross-evaluate: does W from one solution favor the other?
cat(sprintf("\nObj(alpha_pos, W(alpha_neg)): %.8f\n", gmm_obj(alpha_pos, W_neg)))
cat(sprintf("Obj(alpha_neg, W(alpha_pos)): %.8f\n", gmm_obj(alpha_neg, W_pos)))

# Check moment conditions at each
G_pos <- colMeans(g_pos)
G_neg <- colMeans(g_neg)
G_glm <- colMeans(g_glm)
cat(sprintf("\nG(alpha_pos): (%s)\n", paste(round(G_pos, 6), collapse=", ")))
cat(sprintf("G(alpha_neg): (%s)\n", paste(round(G_neg, 6), collapse=", ")))
cat(sprintf("G(alpha_glm): (%s)\n", paste(round(G_glm, 6), collapse=", ")))

# Now trace the iterative GMM to see the cycling
cat("\n\n=== Tracing iterative GMM from GLM init ===\n")
alpha_t <- alpha_glm
for (t in 1:30) {
  g_t <- Phi_alpha(alpha_t)
  G_t <- colMeans(g_t)
  W_t <- tryCatch(solve(crossprod(g_t)/n), error=function(e) diag(ncol(h_full)))
  obj_t <- as.numeric(t(G_t) %*% W_t %*% G_t)

  # Inner optim with this fixed W
  obj_fn <- function(a) {
    g_m <- Phi_alpha(a)
    G_v <- colMeans(g_m)
    as.numeric(t(G_v) %*% W_t %*% G_v)
  }
  opt <- optim(alpha_t, obj_fn, method="L-BFGS-B", control=list(maxit=1000))
  alpha_new <- opt$par

  cat(sprintf("  t=%2d: alpha=(%s), obj=%.6f, ||step||=%.4f\n",
      t, paste(round(alpha_t, 3), collapse=", "), obj_t,
      sqrt(sum((alpha_new - alpha_t)^2))))

  alpha_t <- alpha_new
}

# Check: what does the design matrix look like?
cat("\n\n=== Design matrix properties ===\n")
cat(sprintf("Correlation between design_mat columns:\n"))
print(round(cor(design_mat), 3))

cat(sprintf("\nCorrelation between h_full columns:\n"))
print(round(cor(h_full), 3))
