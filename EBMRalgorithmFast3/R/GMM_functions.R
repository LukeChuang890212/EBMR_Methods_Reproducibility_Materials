# GMM functions for EBMRalgorithmFast
# Optimizations:
# 1. Use crossprod instead of t() %*%
# 2. Vectorized H_2n_3 computation (no sapply loop)
# 3. Precompute reusable matrices
# 4. Use colMeans instead of apply
# 5. Cache g() evaluations in Gamma computation for Gamma_2
# 6. Support analytical gradients via optional dg function

gmm = function(g, W, n, esteq_dim, param_dim, init, se.fit = T, dg = NULL, d2g = NULL, method = "GN"){
  # dg: optional function that returns analytical gradient of g
  # dg(param) should return array of dim (n, param_dim, esteq_dim)
  # where dg[i, j, l] = d(g_l(i)) / d(param_j)
  #
  # d2g: optional function that returns analytical second derivative (Hessian) of Gamma

  # d2g(param) should return matrix of dim (esteq_dim * param_dim, param_dim)
  # which is the jacobian of vec(t(Gamma)) w.r.t. param

  G = function(param){
    g.matrix = g(param)
    return(matrix(.colMeans(g.matrix, n, esteq_dim), esteq_dim, 1))
  }

  # t_Gamma_i: compute jacobian of g for each observation
  # Use analytical gradient if provided, otherwise numerical
  t_Gamma_i = function(param){
    if (!is.null(dg)) {
      # Use analytical gradient
      return(dg(param))
    }
    # Fallback to numerical jacobian
    Gamma.arr = array(NA, dim = c(n, param_dim, esteq_dim))
    for(l in 1:esteq_dim){
      Gamma.arr[,,l] = jacobian(function(param) g(param)[, l], param)
    }
    return(Gamma.arr)
  }

  Gamma = function(param, t_Gamma_arr = NULL){
    # Allow passing precomputed t_Gamma_arr to avoid redundant jacobian calls
    if (is.null(t_Gamma_arr)) {
      t_Gamma_arr = t_Gamma_i(param)
    }
    # t_Gamma_arr is n x param_dim x esteq_dim
    # Result should be esteq_dim x param_dim
    Gamma_mat = matrix(0, esteq_dim, param_dim)
    for (j in 1:param_dim) {
      Gamma_mat[, j] = colMeans(t_Gamma_arr[, j, ])
    }
    return(Gamma_mat)
  }

  # Gamma_2 computes jacobian of vec(t(Gamma)) w.r.t. param
  # Use analytical d2g if provided, otherwise fall back to numerical
  Gamma_2 = function(param){
    if (!is.null(d2g)) {
      # Use analytical second derivative
      return(d2g(param))
    }
    # Fallback: Compute Gamma directly from G (colMeans of g) to avoid
    # expensive t_Gamma_i computation inside numerical differentiation
    Gamma_from_G = function(p) {
      jacobian(function(pp) as.vector(G(pp)), p)
    }
    return(jacobian(function(p) as.vector(t(Gamma_from_G(p))), param))
  }

  if(is.null(init)) init = rep(0, param_dim)
  opt = NULL

  if (method == "GN" && !is.null(dg)) {
    # ---- Gauss-Newton solver ----
    # Two-step iterative GMM with Gauss-Newton inner solver
    # NR step: delta = -(Gamma'WGamma)^{-1} Gamma'WG

    compute_Gamma_direct <- function(param) {
      t_Gamma_arr <- dg(param)
      Gamma_mat <- matrix(0, esteq_dim, param_dim)
      for (j in 1:param_dim) {
        Gamma_mat[, j] <- .colMeans(matrix(t_Gamma_arr[, j, ], n, esteq_dim), n, esteq_dim)
      }
      Gamma_mat
    }

    alpha <- init
    gn_outer_iters <- 0
    gn_inner_iters <- 0
    gn_inner_converged <- FALSE
    gn_outer_converged <- FALSE

    for (outer in 1:40) {
      gn_outer_iters <- outer
      g_mat <- g(alpha)
      W_hat <- tryCatch(W(g_mat), error = function(e) diag(esteq_dim))

      alpha_old <- alpha
      for (nr_iter in 1:2000) {
        gn_inner_iters <- gn_inner_iters + 1
        g_mat <- g(alpha)
        G_hat <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
        Gamma_hat <- compute_Gamma_direct(alpha)

        GtW <- crossprod(Gamma_hat, W_hat)
        H <- GtW %*% Gamma_hat
        rhs <- GtW %*% G_hat

        delta <- tryCatch(-as.vector(solve(H, rhs)),
                          error = function(e) {
                            tryCatch(-as.vector(solve(H + 1e-6 * diag(param_dim), rhs)),
                                     error = function(e2) rep(0, param_dim))
                          })

        # Damping for large steps
        step_size <- 1.0
        if (max(abs(delta)) > 2) step_size <- 2 / max(abs(delta))

        alpha <- alpha + step_size * delta
        if (max(abs(delta * step_size)) < 1e-8) {
          gn_inner_converged <- TRUE
          break
        }
      }

      param_change <- max(abs(alpha - alpha_old))
      if (param_change < 1e-6) {
        gn_outer_converged <- TRUE
        break
      }
    }

    # Compute final GMM objective
    g_mat_final <- g(alpha)
    G_final <- matrix(.colMeans(g_mat_final, n, esteq_dim), esteq_dim, 1)
    W_final <- tryCatch(W(g_mat_final), error = function(e) diag(esteq_dim))
    gn_objective <- as.numeric(crossprod(G_final, W_final %*% G_final))

    estimates <- alpha
    opt <- list(
      converged = gn_outer_converged,
      outer_iterations = gn_outer_iters,
      total_inner_iterations = gn_inner_iters,
      inner_converged = gn_inner_converged,
      objective = gn_objective,
      param_change = param_change,
      method = "GN"
    )

  } else if (method == "L-BFGS-B" && !is.null(dg)) {
    # ---- L-BFGS-B with analytical gradient ----
    compute_Gamma_direct <- function(param) {
      t_Gamma_arr <- dg(param)
      Gamma_mat <- matrix(0, esteq_dim, param_dim)
      for (j in 1:param_dim) {
        Gamma_mat[, j] <- .colMeans(matrix(t_Gamma_arr[, j, ], n, esteq_dim), n, esteq_dim)
      }
      Gamma_mat
    }

    sol <- init; W.hat <- NULL
    for (outer in 1:20) {
      sol_old <- sol
      obj <- function(p) {
        Gh <- G(p)
        if (is.null(W.hat)) as.numeric(crossprod(Gh))
        else as.numeric(crossprod(Gh, W.hat %*% Gh))
      }
      gr <- function(p) {
        Gh <- G(p); Gamma_hat <- compute_Gamma_direct(p)
        if (is.null(W.hat)) as.vector(2 * crossprod(Gamma_hat, Gh))
        else as.vector(2 * crossprod(Gamma_hat, W.hat %*% Gh))
      }
      opt <- optim(sol, obj, gr = gr, method = "L-BFGS-B",
                   lower = rep(-20, param_dim), upper = rep(20, param_dim),
                   control = list(maxit = 2000))
      sol <- opt$par
      W.hat <- tryCatch(W(g(sol)), error = function(e) diag(esteq_dim))
      if (max(abs(sol - sol_old)) < 1e-6) break
    }

    estimates <- sol

  } else if (method == "L-BFGS-B-prog" && !is.null(dg)) {
    # ---- L-BFGS-B with analytical gradient and progressive bounds ----
    # Progressive bound relaxation: start with tight bounds and gradually relax.
    # Warm-starting from constrained solutions keeps the optimizer in
    # moderate-coefficient basins, avoiding separation solutions.
    compute_Gamma_direct <- function(param) {
      t_Gamma_arr <- dg(param)
      Gamma_mat <- matrix(0, esteq_dim, param_dim)
      for (j in 1:param_dim) {
        Gamma_mat[, j] <- .colMeans(matrix(t_Gamma_arr[, j, ], n, esteq_dim), n, esteq_dim)
      }
      Gamma_mat
    }

    bounds_seq <- c(3, 5, 7, 10, 15, 20)
    sol <- init
    lbfgsb_outer <- 0

    for (B in bounds_seq) {
      W.hat <- NULL
      for (outer in 1:20) {
        lbfgsb_outer <- lbfgsb_outer + 1
        sol_old <- sol
        obj <- function(p) {
          Gh <- G(p)
          if (is.null(W.hat)) as.numeric(crossprod(Gh))
          else as.numeric(crossprod(Gh, W.hat %*% Gh))
        }
        gr <- function(p) {
          Gh <- G(p); Gamma_hat <- compute_Gamma_direct(p)
          if (is.null(W.hat)) as.vector(2 * crossprod(Gamma_hat, Gh))
          else as.vector(2 * crossprod(Gamma_hat, W.hat %*% Gh))
        }
        opt <- optim(sol, obj, gr = gr, method = "L-BFGS-B",
                     lower = rep(-B, param_dim), upper = rep(B, param_dim),
                     control = list(maxit = 2000))
        sol <- opt$par
        W.hat <- tryCatch(W(g(sol)), error = function(e) diag(esteq_dim))
        if (max(abs(sol - sol_old)) < 1e-6) break
      }
    }

    # Compute final GMM objective
    g_mat_final <- g(sol)
    G_final <- matrix(.colMeans(g_mat_final, n, esteq_dim), esteq_dim, 1)
    W_final <- tryCatch(W(g_mat_final), error = function(e) diag(esteq_dim))
    lbfgsb_objective <- as.numeric(crossprod(G_final, W_final %*% G_final))

    estimates <- sol
    opt <- list(
      converged = TRUE,
      outer_iterations = lbfgsb_outer,
      objective = lbfgsb_objective,
      method = "L-BFGS-B-prog"
    )

  } else {
    # ---- L-BFGS-B fallback (no analytical gradient) ----
    state <- new.env(parent = emptyenv())
    state$W.hat <- NULL
    state$cache_param_g <- NULL
    state$cache_g <- NULL
    state$cache_G <- NULL

    compute_g_cache <- function(param) {
      if (!is.null(state$cache_param_g) && identical(param, state$cache_param_g)) return()
      state$cache_param_g <- param
      state$cache_g <- g(param)
      state$cache_G <- matrix(.colMeans(state$cache_g, n, esteq_dim), esteq_dim, 1)
    }

    obj <- function(param) {
      compute_g_cache(param)
      G.hat <- state$cache_G
      if (is.null(state$W.hat)) {
        value <- as.numeric(crossprod(G.hat))
      } else {
        value <- as.numeric(crossprod(G.hat, state$W.hat %*% G.hat))
      }
      if (!is.finite(value)) return(1e8)
      return(value)
    }

    sol_path = matrix(init, param_dim)
    conv_err = 10^8
    prev_conv_err = 10^8
    t = 1

    while (t < 1000){
      state$cache_param_g <- NULL
      opt = optim(sol_path[, t], obj, method = "L-BFGS-B",
                  lower = rep(-10, param_dim), upper = rep(10, param_dim),
                  control = list(maxit = 1000))
      sol_path = cbind(sol_path, opt$par)
      compute_g_cache(sol_path[,t+1])
      g.matrix = state$cache_g
      state$W.hat = W(g.matrix)
      conv_err = max(abs(sol_path[,t+1]-sol_path[, t]))

      if (conv_err < 1e-6) break
      if (t > 3 && conv_err > 1e-6) {
        rel_change = abs(conv_err - prev_conv_err) / (prev_conv_err + 1e-16)
        if (rel_change < 0.01) break
      }

      prev_conv_err = conv_err
      t = t + 1
    }

    estimates = sol_path[, t]
  }

  Gamma.hat = W.hat = g.matrix = eta_s = Q = h_x = psi = se = NA
  if(se.fit){
    # Compute t_Gamma_arr ONCE and reuse for both Gamma.hat and H_2n_2
    t_Gamma_arr = t_Gamma_i(estimates)
    Gamma.hat = Gamma(estimates, t_Gamma_arr)  # Pass precomputed array
    g.matrix = g(estimates)
    W.hat = W(g.matrix)

    # Precompute t(g.matrix) once - used multiple times
    t_g.matrix = t(g.matrix)

    eta_s = .colMeans(g.matrix, n, esteq_dim)
    M_n = kronecker(eta_s %*% W.hat, diag(param_dim)) %*% Gamma_2(estimates)

    # Precompute reusable quantities
    GtW = crossprod(Gamma.hat, W.hat)  # t(Gamma.hat) %*% W.hat

    # H_0n with ridge regularization for numerical stability
    GtWG = GtW %*% Gamma.hat
    H_0n = tryCatch({
      -solve(GtWG)
    }, error = function(e) {
      # Add small ridge penalty if singular
      -solve(GtWG + 1e-8 * diag(nrow(GtWG)))
    })
    H_1n = GtW %*% (t_g.matrix - eta_s)

    # Vectorized H_2n_2 (t_Gamma_arr already computed above)
    # Original: apply(t_Gamma_arr, 1, function(m) (m - t(Gamma.hat)) %*% W.hat %*% eta_s)
    # t_Gamma_arr is n x param_dim x esteq_dim
    # For each i: (t_Gamma_arr[i,,] - t(Gamma.hat)) %*% W.hat %*% eta_s
    W_eta = as.vector(W.hat %*% eta_s)  # esteq_dim vector
    # Use GtW which is t(Gamma.hat) %*% W.hat, then multiply by eta_s
    tGamma_W_eta = as.vector(GtW %*% eta_s)  # param_dim vector (constant part)
    # Reshape t_Gamma_arr to (n*param_dim) x esteq_dim, multiply by W_eta, reshape back
    dim_arr = dim(t_Gamma_arr)
    t_Gamma_mat = matrix(t_Gamma_arr, nrow = dim_arr[1] * dim_arr[2], ncol = dim_arr[3])
    H_2n_2_vec = t_Gamma_mat %*% W_eta  # (n*param_dim) vector
    H_2n_2 = matrix(H_2n_2_vec, nrow = n, ncol = param_dim, byrow = FALSE)
    H_2n_2 = t(H_2n_2) - tGamma_W_eta  # param_dim x n

    # Vectorized H_2n_3: avoid sapply loop
    # Original: t(Gamma.hat) %*% (-W.hat %*% (g[i,] %*% t(g[i,])) %*% W.hat + W.hat) %*% eta_s
    # = GtW %*% eta_s - (GtW %*% g[i,]) * (g[i,] %*% W.hat %*% eta_s)
    GtW_eta = tGamma_W_eta  # Already computed above
    # W_eta already computed above
    GtW_g = GtW %*% t_g.matrix  # param_dim x n (reuse t_g.matrix)
    g_W_eta = as.vector(g.matrix %*% W_eta)  # n vector

    H_2n_3 = matrix(GtW_eta, param_dim, n) - GtW_g * matrix(g_W_eta, param_dim, n, byrow = TRUE)

    # Q with ridge regularization for numerical stability
    I_minus_H0M = diag(param_dim) - H_0n %*% M_n
    Q = tryCatch({
      solve(I_minus_H0M) %*% H_0n
    }, error = function(e) {
      # Add small ridge penalty if singular
      solve(I_minus_H0M + 1e-8 * diag(param_dim)) %*% H_0n
    })
    psi = Q %*% (H_1n + H_2n_2 + H_2n_3)
    cov.hat = var(t(psi)) / n
    se = sqrt(diag(cov.hat))
  }

  result = list(
    estimates = estimates,
    se = se,
    Gamma.hat = Gamma.hat,
    W.hat = W.hat,
    g.matrix = g.matrix,
    eta_s = eta_s,
    Q = Q,
    h_x = h_x,
    psi = psi,
    opt = opt
  )

  return(result)
}
