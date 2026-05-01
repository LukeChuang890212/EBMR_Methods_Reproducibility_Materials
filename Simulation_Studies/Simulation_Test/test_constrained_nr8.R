setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200
COND_THRESH <- 1e8
TRUST_RADIUS <- 2.0

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

# Track details for non-converged reps
nonconv_details <- list()

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
    design_mat <- model.matrix(formula_j, data=dat)
    h_full <- cbind(1, as.matrix(dat[, h_alpha_vars, drop=FALSE]))
    n <- n_val; p <- ncol(design_mat); h_dim <- ncol(h_full)

    model_fn <- function(alpha) { eta <- as.vector(design_mat %*% alpha); 1/(1+exp(-eta)) }
    Phi_alpha <- function(alpha) { pi_hat <- model_fn(alpha); (r_vec/pi_hat - 1) * h_full }
    Gamma_fn <- function(alpha) {
      pi_hat <- model_fn(alpha)
      crossprod(h_full * (-r_vec*(1-pi_hat)/pi_hat), design_mat) / n
    }
    W_func <- function(g_mat) solve(crossprod(g_mat) / nrow(g_mat))

    nr_inner <- function(start, W_hat) {
      alpha <- start
      for (k in 1:500) {
        g_mat <- Phi_alpha(alpha); G_vec <- colMeans(g_mat)
        Gamma_hat <- Gamma_fn(alpha)
        g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
        if (max(abs(g_vec)) < 1e-8) break
        H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
        direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
        sn <- sqrt(sum(direction^2))
        if (sn > TRUST_RADIUS) direction <- direction * (TRUST_RADIUS / sn)
        obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
        step <- 1.0; accepted <- FALSE
        for (ls in 1:30) {
          cand <- alpha + step * direction
          G_new <- colMeans(Phi_alpha(cand))
          obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
          if (is.finite(obj_new) && obj_new < obj_cur - 1e-4*step*sum(g_vec*direction)) {
            Gamma_c <- Gamma_fn(cand)
            H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
            eig_c <- eigen(H_c, symmetric=TRUE, only.values=TRUE)$values
            cc <- max(eig_c)/max(min(eig_c), 1e-15)
            if (cc < COND_THRESH) { accepted <- TRUE; break }
          }
          step <- step * 0.5
        }
        if (accepted) { alpha <- cand } else { break }
      }
      alpha
    }

    estimates <- nr_inner(rep(0, p), diag(h_dim))
    best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE
    grad_trace <- numeric(0)
    final_t <- 0

    for (t in 1:500) {  # extended to 500
      g_mat <- Phi_alpha(estimates); G_vec <- colMeans(g_mat)
      W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
      Gamma_hat <- Gamma_fn(estimates)
      conv_grad <- 2*as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
      grad_norm <- max(abs(conv_grad))
      grad_trace <- c(grad_trace, grad_norm)

      if (grad_norm < 1e-6) { converged <- TRUE; final_t <- t; break }
      if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
      } else { no_improve <- no_improve+1L
        if (no_improve >= 20L) { estimates <- best_est; final_t <- t; break }
      }
      estimates <- nr_inner(estimates, W_hat)
      final_t <- t
    }

    if (!converged) {
      nonconv_details[[length(nonconv_details)+1]] <- list(
        rep=rep_i, final_grad=best_grad, final_t=final_t,
        no_improve=no_improve, grad_trace=grad_trace,
        reason=if(no_improve >= 20L) "stall" else if(final_t >= 500) "max_iter" else "other"
      )
    }
  }, error = function(e) {
    nonconv_details[[length(nonconv_details)+1]] <<- list(rep=rep_i, reason="error", msg=conditionMessage(e))
  })
}

cat(sprintf("Total non-converged: %d / %d\n\n", length(nonconv_details), n_reps))

# Summarize reasons
reasons <- sapply(nonconv_details, `[[`, "reason")
cat("Reasons:\n")
print(table(reasons))

cat("\nDetails of non-converged reps:\n")
for (d in nonconv_details) {
  if (d$reason == "error") {
    cat(sprintf("  Rep %3d: ERROR - %s\n", d$rep, d$msg))
  } else {
    cat(sprintf("  Rep %3d: reason=%s, final_grad=%.2e, iters=%d, best_grad=%.2e\n",
        d$rep, d$reason, tail(d$grad_trace, 1), d$final_t, d$final_grad))
    # Show last 10 grad values
    gt <- d$grad_trace
    n_show <- min(10, length(gt))
    cat(sprintf("           last %d grads: %s\n", n_show,
        paste(sprintf("%.2e", tail(gt, n_show)), collapse=", ")))
  }
}

# Check: how many would converge with looser criterion?
cat("\nConvergence at different thresholds:\n")
for (tol in c(1e-6, 1e-5, 1e-4, 1e-3)) {
  n_conv <- sum(sapply(nonconv_details, function(d) {
    if (d$reason == "error") return(FALSE)
    d$final_grad < tol
  }))
  cat(sprintf("  tol=%.0e: %d additional would converge (total = %d)\n",
      tol, n_conv, 200 - length(nonconv_details) + n_conv))
}
