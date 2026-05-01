# GMM functions for EBMRalgorithmFast
# Optimizations:
# 1. Use crossprod instead of t() %*%
# 2. Vectorized H_2n_3 computation (no sapply loop)
# 3. Precompute reusable matrices
# 4. Use colMeans instead of apply
# 5. Cache g() evaluations in Gamma computation for Gamma_2
# 6. Support analytical gradients via optional dg function

gmm = function(g, W, n, esteq_dim, param_dim, init, se.fit = T, dg = NULL, d2g = NULL){
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

  # Use an environment to hold state (W.hat and cached values)
  state <- new.env(parent = emptyenv())
  state$W.hat <- NULL
  state$cache_param_g <- NULL  # Separate cache keys for g and Gamma
  state$cache_param_Gamma <- NULL
  state$cache_g <- NULL
  state$cache_G <- NULL
  state$cache_Gamma <- NULL

  # Helper to compute and cache g and G for a given param (for objective function)
  compute_g_cache <- function(param) {
    if (!is.null(state$cache_param_g) && identical(param, state$cache_param_g)) {
      return()  # Cache is valid
    }
    state$cache_param_g <- param
    state$cache_g <- g(param)
    state$cache_G <- matrix(.colMeans(state$cache_g, n, esteq_dim), esteq_dim, 1)
  }

  # Helper to compute and cache Gamma for a given param (for gradient - only when needed)
  compute_Gamma_cache <- function(param) {
    if (!is.null(state$cache_param_Gamma) && identical(param, state$cache_param_Gamma)) {
      return()  # Cache is valid
    }
    if (!is.null(dg)) {
      state$cache_param_Gamma <- param
      t_Gamma_arr <- dg(param)  # n x param_dim x esteq_dim
      Gamma_mat <- matrix(0, esteq_dim, param_dim)
      for (j in 1:param_dim) {
        Gamma_mat[, j] <- .colMeans(matrix(t_Gamma_arr[, j, ], n, esteq_dim), n, esteq_dim)
      }
      state$cache_Gamma <- Gamma_mat
    }
  }

  # Single objective function that uses state$W.hat
  # Optimized: use if/else instead of ifelse (faster)
  obj <- function(param) {
    compute_g_cache(param)  # Only compute g and G, not Gamma
    G.hat <- state$cache_G
    if (is.null(state$W.hat)) {
      value <- as.numeric(crossprod(G.hat))
    } else {
      value <- as.numeric(crossprod(G.hat, state$W.hat %*% G.hat))
    }
    if (!is.finite(value)) return(1e8)
    return(value)
  }

  # Analytical gradient of objective function: d(obj)/d(param) = 2 * Gamma' W G
  # where Gamma = d(G)/d(param) is esteq_dim x param_dim
  # Only use analytical gradient if dg is provided
  grad_obj <- if (!is.null(dg)) {
    function(param) {
      compute_g_cache(param)     # Ensure G is computed
      compute_Gamma_cache(param) # Compute Gamma only when gradient is needed
      G.hat <- state$cache_G
      Gamma_mat <- state$cache_Gamma

      if (is.null(state$W.hat)) {
        grad <- 2 * as.vector(crossprod(Gamma_mat, G.hat))
      } else {
        grad <- 2 * as.vector(crossprod(Gamma_mat, state$W.hat %*% G.hat))
      }
      return(grad)
    }
  } else {
    NULL  # No analytical gradient available
  }

  if(is.null(init)) init = rep(0, param_dim)
  sol_path = matrix(init, param_dim)
  conv_err = 10^8
  prev_conv_err = 10^8
  t = 1

  # Improved convergence: stop when either:
  # 1. conv_err < 1e-8 (absolute tolerance)
  # 2. conv_err is not improving significantly (relative change < 1%)
  # 3. t >= 1000 (max iterations)
  #
  # Note: Using L-BFGS-B for optimization. It provides significant speedup
  # especially for binary outcomes (4-24x faster than BFGS for setting14).
  # Bounds of [-10, 10] prevent extreme coefficients while being wide enough
  # for typical propensity score parameters.
  while (t < 1000){
    # Clear cache before optimization (new W.hat)
    state$cache_param_g <- NULL
    state$cache_param_Gamma <- NULL
    # L-BFGS-B with bounds - faster for binary y, similar for continuous y
    opt = optim(sol_path[, t], obj, gr = grad_obj, method = "L-BFGS-B",
                lower = rep(-100, param_dim), upper = rep(100, param_dim),
                control = list(maxit = 1000))

    # Original BFGS (uncomment if L-BFGS-B causes issues)
    # opt = optim(sol_path[, t], obj, gr = grad_obj, method = "BFGS",
    #             control = list(maxit = 1000))
    sol_path = cbind(sol_path, opt$par)
    # Use cached g.matrix if available from last obj/grad call, otherwise compute
    compute_g_cache(sol_path[,t+1])
    g.matrix = state$cache_g
    state$W.hat = W(g.matrix)  # Update W.hat in state
    conv_err = max(abs(sol_path[,t+1]-sol_path[, t]))

    # Check convergence
    if (conv_err < 1e-6) break  # Absolute convergence

    # Check if conv_err stopped improving (stagnation detection)
    if (t > 3 && conv_err > 1e-6) {
      rel_change = abs(conv_err - prev_conv_err) / (prev_conv_err + 1e-16)
      if (rel_change < 0.01) break  # Less than 1% improvement
    }

    prev_conv_err = conv_err
    t = t + 1
  }

  estimates = sol_path[, t]

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

    H_0n = -solve(GtW %*% Gamma.hat)
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

    Q = solve(diag(param_dim) - H_0n %*% M_n) %*% H_0n
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
    opt = opt,
    obj = obj
  )

  return(result)
}
