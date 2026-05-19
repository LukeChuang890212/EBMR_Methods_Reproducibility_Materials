# GMM functions for EPS
#
# Optimizations:
# 1. Use crossprod instead of t() %*%
# 2. Vectorized H_2n_3 computation (no sapply loop)
# 3. Precompute reusable matrices
# 4. Use colMeans instead of apply
# 5. Cache g() evaluations in Gamma computation for Gamma_2
# 6. Support analytical gradients via optional dg function

gmm = function(g, W, n, esteq_dim, param_dim, init, se.fit = T, dg = NULL, d2g = NULL, lower = -Inf, upper = Inf, stall_limit = 20L, max_outer_override = NULL, Gamma_direct = NULL, optimizer = c("L-BFGS-B", "constrained_nr"), cond_threshold = 1e8, trust_radius = 2.0){
  optimizer <- match.arg(optimizer)
  # Gamma_direct: optional fast path for Gamma = E[dg/dparam]
  # Gamma_direct(param) should return matrix of dim (esteq_dim, param_dim)
  # This avoids the full 3D array allocation in dg() when only the mean is needed.
  # dg is still required for SE computation (t_Gamma_i needs the full 3D array).

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

  # ===========================================================================
  # Main optimization: Iteratively update W and minimize G'W_tG
  # ===========================================================================

  # Gamma_fast: compute Gamma matrix efficiently (used by both optimizers)
  Gamma_fast <- if (!is.null(Gamma_direct)) {
    Gamma_direct
  } else if (!is.null(dg)) {
    function(param) {
      dg_arr <- dg(param)
      Gm <- matrix(0, esteq_dim, param_dim)
      for (j in 1:param_dim) Gm[, j] <- .colMeans(dg_arr[, j, , drop = FALSE], n, esteq_dim)
      Gm
    }
  } else { NULL }

  max_outer = if (!is.null(max_outer_override)) max_outer_override else 200L
  tol = 1e-6

  if (optimizer == "L-BFGS-B") {
    # =========================================================================
    # L-BFGS-B optimizer (baseline)
    # =========================================================================

    # Fixed-W objective and gradient (reused across outer iterations)
    state <- new.env(parent = emptyenv())
    state$W.hat <- NULL
    state$cache_param <- NULL
    state$cache_g <- NULL
    state$cache_G <- NULL

    cache <- function(param) {
      if (!is.null(state$cache_param) && identical(param, state$cache_param)) return()
      state$cache_param <- param
      state$cache_g <- g(param)
      state$cache_G <- matrix(.colMeans(state$cache_g, n, esteq_dim), esteq_dim, 1)
    }

    obj <- function(param) {
      cache(param)
      v <- as.numeric(crossprod(state$cache_G, state$W.hat %*% state$cache_G))
      if (!is.finite(v)) 1e8 else v
    }

    # Inner-loop gradient: 2 * Gamma(alpha)' W G(alpha), W fixed from outer iteration
    grad <- if (!is.null(Gamma_fast)) {
      function(param) {
        cache(param)
        Gamma_mat <- Gamma_fast(param)
        2 * as.vector(crossprod(Gamma_mat, state$W.hat %*% state$cache_G))
      }
    } else { NULL }

    # Convergence criterion: 2 * Gamma(alpha)' W(alpha) G(alpha), W at current alpha.
    conv_grad <- if (!is.null(Gamma_fast)) {
      function(param) {
        g_mat <- g(param)
        G_vec <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
        W_hat <- tryCatch(W(g_mat), error = function(e) diag(esteq_dim))
        2 * as.vector(crossprod(Gamma_fast(param), W_hat %*% G_vec))
      }
    } else { NULL }

    # Fast version for the outer loop: reuses state$cache_G and state$W.hat
    conv_grad_fast <- if (!is.null(Gamma_fast)) {
      function(param) {
        2 * as.vector(crossprod(Gamma_fast(param), state$W.hat %*% state$cache_G))
      }
    } else { NULL }

    # Projected gradient: zero out components where alpha is at bound and gradient pushes outward
    proj_grad_norm <- function(grad_vec, param) {
      pg <- grad_vec
      for (j in seq_along(param)) {
        if (is.finite(lower[j]) && param[j] <= lower[j] + 1e-10 && pg[j] > 0) pg[j] <- 0
        if (is.finite(upper[j]) && param[j] >= upper[j] - 1e-10 && pg[j] < 0) pg[j] <- 0
      }
      max(abs(pg))
    }

    # Step 1: W = I, minimize G(alpha)'G(alpha)
    state$W.hat <- diag(esteq_dim)
    opt_init = optim(init, obj, gr = grad, method = "L-BFGS-B",
                     lower = lower, upper = upper,
                     control = list(maxit = 1000))
    estimates = opt_init$par

    # Step 2: Iterative GMM with W updates
    alpha_change <- Inf
    for (t in 1:max_outer) {
      cache(estimates)
      state$W.hat = tryCatch(W(state$cache_g), error = function(e) diag(esteq_dim))

      # Check convergence: |alpha_t - alpha_{t-1}| < tol
      if (alpha_change < tol) break

      state$cache_param <- NULL
      estimates_old <- estimates
      opt_t = optim(estimates, obj, gr = grad, method = "L-BFGS-B",
                    lower = lower, upper = upper,
                    control = list(maxit = 1000))
      estimates = opt_t$par
      alpha_change <- max(abs(estimates - estimates_old))
    }
    iter_total = t
    final_grad_norm = alpha_change

    # ----- Fallback: if L-BFGS-B did not converge, try constrained_nr -----
    if (!is.na(final_grad_norm) && final_grad_norm >= tol && !is.null(Gamma_fast)) {
      lbfgsb_estimates <- estimates
      lbfgsb_grad_norm <- final_grad_norm

      # constrained_nr inner optimizer (fast: degeneracy guard via cond check on H_gn)
      nr_inner_fb <- function(start, W_hat) {
        alpha <- start
        for (k in 1:100) {
          Gamma_hat <- Gamma_fast(alpha)
          g_mat <- g(alpha)
          G_vec <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
          g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
          if (max(abs(g_vec)) < 1e-8) break
          H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
          direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
          sn <- sqrt(sum(direction^2))
          if (sn > trust_radius) direction <- direction * (trust_radius / sn)
          obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
          step <- 1.0
          accepted <- FALSE
          for (ls in 1:15) {
            cand <- alpha + step * direction
            g_cand <- g(cand)
            G_cand <- matrix(.colMeans(g_cand, n, esteq_dim), esteq_dim, 1)
            obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
            if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
              # Fast degeneracy check: use rcond instead of eigen
              Gamma_c <- Gamma_fast(cand)
              H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
              if (rcond(H_c) > 1/cond_threshold) { accepted <- TRUE; break }
            }
            step <- step * 0.5
          }
          if (accepted) { alpha <- cand } else { break }
        }
        alpha
      }

      # Run constrained_nr from fresh init
      est_nr <- nr_inner_fb(init, diag(esteq_dim))  # Step 1: W=I
      alpha_change_nr <- Inf

      for (t_nr in 1:max_outer) {
        if (alpha_change_nr < tol) break
        g_mat_nr <- g(est_nr)
        W_nr <- tryCatch(W(g_mat_nr), error = function(e) diag(esteq_dim))
        est_nr_old <- est_nr
        est_nr <- nr_inner_fb(est_nr, W_nr)
        alpha_change_nr <- max(abs(est_nr - est_nr_old))
      }

      # Check final grad norm
      g_nr_f <- g(est_nr)
      G_nr_f <- matrix(.colMeans(g_nr_f, n, esteq_dim), esteq_dim, 1)
      W_nr_f <- tryCatch(W(g_nr_f), error = function(e) diag(esteq_dim))
      nr_grad_norm <- max(abs(2 * as.vector(crossprod(Gamma_fast(est_nr), W_nr_f %*% G_nr_f))))

      # Accept if constrained_nr did better
      if (nr_grad_norm < lbfgsb_grad_norm) {
        estimates <- est_nr
        final_grad_norm <- nr_grad_norm
      }
    }

  } else {
    # =========================================================================
    # Constrained NR optimizer
    # =========================================================================
    # Custom Gauss-Newton inner optimizer with:
    #   - Trust region (clamp step norm to trust_radius)
    #   - Degeneracy guard: reject steps where cond(Gamma'W Gamma) > cond_threshold
    # This prevents the optimizer from entering the flat plateau where pi -> 1,
    # which causes zero gradient AND zero Hessian (spurious convergence).

    if (is.null(Gamma_fast)) {
      stop("constrained_nr optimizer requires Gamma_direct or dg")
    }

    # Inner optimizer: fully minimize Q(alpha) = G(alpha)' W G(alpha) for fixed W,
    # subject to cond(Gamma'W Gamma) < cond_threshold
    nr_inner <- function(start, W_hat) {
      alpha <- start
      for (k in 1:500) {
        g_mat <- g(alpha)
        G_vec <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
        Gamma_hat <- Gamma_fast(alpha)
        g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
        if (max(abs(g_vec)) < 1e-8) break

        # Gauss-Newton direction
        H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
        direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)

        # Trust region clamp
        sn <- sqrt(sum(direction^2))
        if (sn > trust_radius) direction <- direction * (trust_radius / sn)

        # Backtracking line search with degeneracy guard (rcond)
        obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
        step <- 1.0
        accepted <- FALSE
        for (ls in 1:30) {
          cand <- alpha + step * direction
          g_cand <- g(cand)
          G_cand <- matrix(.colMeans(g_cand, n, esteq_dim), esteq_dim, 1)
          obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
          # Armijo condition + degeneracy guard
          if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
            Gamma_c <- Gamma_fast(cand)
            H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
            if (rcond(H_c) > 1/cond_threshold) { accepted <- TRUE; break }
          }
          step <- step * 0.5
        }
        if (accepted) { alpha <- cand } else { break }
      }
      alpha
    }

    # Step 1: W = I
    estimates <- nr_inner(init, diag(esteq_dim))

    # Step 2: Iterative GMM with W updates
    alpha_change <- Inf
    for (t in 1:max_outer) {
      if (alpha_change < tol) break

      g_mat <- g(estimates)
      W_hat <- tryCatch(W(g_mat), error = function(e) diag(esteq_dim))

      estimates_old <- estimates
      estimates <- nr_inner(estimates, W_hat)
      alpha_change <- max(abs(estimates - estimates_old))
    }
    iter_total = t
    final_grad_norm <- alpha_change

    # ----- Fallback: if constrained_nr did not converge, try L-BFGS-B -----
    if (final_grad_norm >= tol) {
      cnr_estimates <- estimates
      cnr_grad_norm <- final_grad_norm

      # L-BFGS-B iterative GMM subroutine
      run_lbfgsb_gmm <- function(start_point) {
        st <- new.env(parent = emptyenv())
        st$W.hat <- NULL; st$cache_param <- NULL
        st$cache_g <- NULL; st$cache_G <- NULL

        cache_fn <- function(param) {
          if (!is.null(st$cache_param) && identical(param, st$cache_param)) return()
          st$cache_param <- param
          st$cache_g <- g(param)
          st$cache_G <- matrix(.colMeans(st$cache_g, n, esteq_dim), esteq_dim, 1)
        }
        obj_fn <- function(param) {
          cache_fn(param)
          v <- as.numeric(crossprod(st$cache_G, st$W.hat %*% st$cache_G))
          if (!is.finite(v)) 1e8 else v
        }
        grad_fn <- function(param) {
          cache_fn(param)
          2 * as.vector(crossprod(Gamma_fast(param), st$W.hat %*% st$cache_G))
        }

        est <- start_point
        # Step 1: W = I
        st$W.hat <- diag(esteq_dim)
        opt0 <- optim(est, obj_fn, gr = grad_fn, method = "L-BFGS-B",
                       control = list(maxit = 1000))
        est <- opt0$par

        # Step 2: Iterative W updates
        for (t_fb in 1:max_outer) {
          cache_fn(est)
          st$W.hat <- tryCatch(W(st$cache_g), error = function(e) diag(esteq_dim))
          fg <- max(abs(grad_fn(est)))
          if (fg < tol) break
          st$cache_param <- NULL
          opt_t <- optim(est, obj_fn, gr = grad_fn, method = "L-BFGS-B",
                         control = list(maxit = 1000))
          est <- opt_t$par
        }
        g_f <- g(est)
        G_f <- matrix(.colMeans(g_f, n, esteq_dim), esteq_dim, 1)
        W_f <- tryCatch(W(g_f), error = function(e) diag(esteq_dim))
        list(estimates = est,
             grad_norm = max(abs(2 * as.vector(crossprod(Gamma_fast(est), W_f %*% G_f)))))
      }

      # Try L-BFGS-B from two starting points: constrained_nr's best and init (zero)
      fb1 <- run_lbfgsb_gmm(cnr_estimates)
      fb2 <- run_lbfgsb_gmm(init)

      # Pick the better fallback
      fb <- if (fb1$grad_norm <= fb2$grad_norm) fb1 else fb2

      # Accept if fallback converged better than constrained_nr
      if (fb$grad_norm < cnr_grad_norm) {
        estimates <- fb$estimates
        final_grad_norm <- fb$grad_norm
      }
    }
  }

  # ===========================================================================
  # Compute objective and convergence info
  # ===========================================================================
  g_mat_final = g(estimates)
  G_final = matrix(.colMeans(g_mat_final, n, esteq_dim), esteq_dim, 1)
  W_final = tryCatch(W(g_mat_final), error = function(e) diag(esteq_dim))
  objective = as.numeric(crossprod(G_final, W_final %*% G_final))
  # Compute cond(GtWG) at solution for diagnostics
  Gamma_final <- if (!is.null(Gamma_fast)) Gamma_fast(estimates) else NULL
  solution_cond <- NA_real_
  if (!is.null(Gamma_final)) {
    H_final <- crossprod(Gamma_final, W_final %*% Gamma_final)
    eig_final <- eigen(H_final, symmetric = TRUE, only.values = TRUE)$values
    solution_cond <- max(eig_final) / max(min(eig_final), 1e-15)
  }

  opt = list(
    converged = !is.na(final_grad_norm) && final_grad_norm < tol,
    objective = objective,
    iterations = iter_total,
    final_grad_norm = final_grad_norm,
    solution_cond = solution_cond
  )

  # ===========================================================================
  # Standard error computation
  # ===========================================================================
  Gamma.hat = W.hat = g.matrix = eta_s = Q = h_x = psi = se = NA
  if(se.fit){
    # Compute t_Gamma_arr ONCE and reuse for both Gamma.hat and SE
    t_Gamma_arr = t_Gamma_i(estimates)
    Gamma.hat = Gamma(estimates, t_Gamma_arr)  # Pass precomputed array
    g.matrix = g(estimates)
    W.hat = W(g.matrix)

    # Precompute reusable quantities
    GtW = crossprod(Gamma.hat, W.hat)  # t(Gamma.hat) %*% W.hat
    GtWG = GtW %*% Gamma.hat
    eta_s = .colMeans(g.matrix, n, esteq_dim)

    # # ---- SE1: one-step GMM standard error (commented out) ----
    # t_g.matrix = t(g.matrix)
    # M_n = kronecker(eta_s %*% W.hat, diag(param_dim)) %*% Gamma_2(estimates)
    # H_0n = tryCatch({
    #   -solve(GtWG)
    # }, error = function(e) {
    #   -solve(GtWG + 1e-8 * diag(nrow(GtWG)))
    # })
    # H_1n = GtW %*% (t_g.matrix - eta_s)
    # W_eta = as.vector(W.hat %*% eta_s)
    # tGamma_W_eta = as.vector(GtW %*% eta_s)
    # dim_arr = dim(t_Gamma_arr)
    # t_Gamma_mat = matrix(t_Gamma_arr, nrow = dim_arr[1] * dim_arr[2], ncol = dim_arr[3])
    # H_2n_2_vec = t_Gamma_mat %*% W_eta
    # H_2n_2 = matrix(H_2n_2_vec, nrow = n, ncol = param_dim, byrow = FALSE)
    # H_2n_2 = t(H_2n_2) - tGamma_W_eta
    # GtW_eta = tGamma_W_eta
    # GtW_g = GtW %*% t_g.matrix
    # g_W_eta = as.vector(g.matrix %*% W_eta)
    # H_2n_3 = matrix(GtW_eta, param_dim, n) - GtW_g * matrix(g_W_eta, param_dim, n, byrow = TRUE)
    # I_minus_H0M = diag(param_dim) - H_0n %*% M_n
    # Q = tryCatch({
    #   solve(I_minus_H0M) %*% H_0n
    # }, error = function(e) {
    #   solve(I_minus_H0M + 1e-8 * diag(param_dim)) %*% H_0n
    # })
    # psi1 = Q %*% (H_1n + H_2n_2 + H_2n_3)
    # cov.hat = var(t(psi1)) / n
    # se1 = sqrt(diag(cov.hat))

    # ---- SE (CUE standard error) ----
    # H = Gamma'WGamma + (G'W ⊗ I_k) R(theta) - (G'W ⊗ Gamma'W) S(theta)
    # R(theta) = d vec(Gamma') / d theta'  [k*h x k]
    # S(theta) = d vec(W^{-1}) / d theta'  [h^2 x k]  -- computed analytically
    # Q_i = Gamma'W g_i + Gamma_i'WG - (Gamma'W g_i)(g_i'WG)  [k x 1]
    # psi_i = -H^{-1} Q_i
    dim_arr = dim(t_Gamma_arr)
    t_Gamma_mat = matrix(t_Gamma_arr, nrow = dim_arr[1] * dim_arr[2], ncol = dim_arr[3])
    tryCatch({
      # ---- R(theta): d vec(Gamma') / d theta' [(k*h) x k] ----
      # Use d2g directly if available (= Gamma_2, already the Jacobian of vec(Gamma')).
      # Fall back to finite differences on Gamma function (2k evaluations).
      R_mat <- if (!is.null(d2g)) {
        d2g(estimates)   # (k*h) x k  -- same object as Gamma_2(estimates)
      } else {
        eps2 <- 1e-5
        Gamma_func2 <- if (!is.null(Gamma_direct)) {
          Gamma_direct
        } else if (!is.null(dg)) {
          function(p) {
            dg_arr2 <- dg(p)
            Gm2 <- matrix(0, esteq_dim, param_dim)
            for (jj in 1:param_dim)
              Gm2[, jj] <- .colMeans(dg_arr2[, jj, , drop = FALSE], n, esteq_dim)
            Gm2
          }
        } else {
          function(p) jacobian(function(pp) as.vector(G(pp)), p)
        }
        Rm <- matrix(0, param_dim * esteq_dim, param_dim)
        for (j in 1:param_dim) {
          ej <- rep(0, param_dim); ej[j] <- eps2
          Rm[, j] <- (as.vector(t(Gamma_func2(estimates + ej))) -
                      as.vector(t(Gamma_func2(estimates - ej)))) / (2 * eps2)
        }
        Rm
      }

      # ---- S(theta): d vec(W^{-1}) / d theta' [h^2 x k] ----
      # Detect whether W depends on theta:
      #   - If W = I (identity or other constant), dW/dtheta = 0, S = 0
      #   - If W = solve(g'g/n), W^{-1} = (1/n) sum g_i g_i'
      #   - If W = solve(var(g)), W^{-1} = (1/(n-1)) sum (g_i-gbar)(g_i-gbar)'
      # Both parameter-dependent cases require analytical S_mat.
      W_inv <- tryCatch(solve(W.hat), error = function(e) NULL)
      W_type <- "constant"  # default: W doesn't depend on theta
      if (!is.null(W_inv)) {
        gg_n <- crossprod(g.matrix) / n        # (1/n) g'g — uncentered second moment
        var_g <- cov(g.matrix)                 # var(g) — centered covariance
        rel_tol <- 1e-6
        denom <- max(max(abs(W_inv)), 1e-10)
        if (max(abs(W_inv - gg_n)) / denom < rel_tol) {
          W_type <- "second_moment"            # W^{-1} = g'g/n
        } else if (max(abs(W_inv - var_g)) / denom < rel_tol) {
          W_type <- "covariance"               # W^{-1} = var(g)
        }
      }

      S_mat <- matrix(0, esteq_dim^2, param_dim)
      if (W_type == "second_moment") {
        # W^{-1} = (1/n) sum g_i g_i'
        # d(W^{-1})/dtheta_j = (1/n) sum [dg_i/dtheta_j g_i' + g_i dg_i/dtheta_j']
        for (j in 1:param_dim) {
          dg_j <- t_Gamma_arr[, j, ]           # n x esteq_dim
          cross <- (crossprod(dg_j, g.matrix) + crossprod(g.matrix, dg_j)) / n
          S_mat[, j] <- as.vector(cross)
        }
      } else if (W_type == "covariance") {
        # W^{-1} = var(g) = (1/(n-1)) sum (g_i - gbar)(g_i - gbar)'
        # d(W^{-1})/dtheta_j = (1/(n-1)) sum [d(g_i-gbar)/dtheta_j (g_i-gbar)'
        #                                    + (g_i-gbar) d(g_i-gbar)/dtheta_j']
        # where d(g_i-gbar)/dtheta_j = dg_i/dtheta_j - (1/n) sum_k dg_k/dtheta_j
        g_centered <- scale(g.matrix, center = TRUE, scale = FALSE)  # g_i - gbar
        for (j in 1:param_dim) {
          dg_j <- t_Gamma_arr[, j, ]           # n x esteq_dim: dg_i/dtheta_j
          dg_j_centered <- scale(dg_j, center = TRUE, scale = FALSE)
          cross <- (crossprod(dg_j_centered, g_centered) +
                    crossprod(g_centered, dg_j_centered)) / (n - 1)
          S_mat[, j] <- as.vector(cross)
        }
      }
      # else W_type == "constant": S_mat stays zero

      # ---- H = Gamma'W Gamma + (G'W ⊗ I_k) R - (G'W ⊗ Gamma'W) S  [k x k] ----
      GW_row      <- as.vector(t(eta_s) %*% W.hat)                              # length h (G'W)
      GW_kron_Ik  <- kronecker(matrix(GW_row, 1, esteq_dim), diag(param_dim))   # k x (h*k)
      GW_kron_GtW <- kronecker(matrix(GW_row, 1, esteq_dim), GtW)               # k x (h*h)
      H_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat

      # ---- Q_i = Gamma'W g_i + Gamma_i'WG - (Gamma'W g_i)(g_i'WG)  [k x n] ----
      WG         <- as.vector(W.hat %*% eta_s)                              # length h
      GtW_g      <- GtW %*% t(g.matrix)                                    # k x n
      g_WG       <- as.vector(g.matrix %*% WG)                             # n-vector

      # Gamma_i'WG for each i: t_Gamma_mat (n*k x h) %*% WG → (n*k) x 1
      # Reshape to n x k, transpose to k x n
      GammaT_WG <- t(matrix(t_Gamma_mat %*% WG, n, param_dim))             # k x n

      Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)             # k x n

      # ---- psi_i = -H^{-1} Q_i; cov from sample variance ----
      # Guard against ill-conditioned H_mat (large condition number → unstable SE)
      H_eig    <- eigen(H_mat, symmetric = FALSE, only.values = TRUE)$values
      cond_num <- max(Mod(H_eig)) / max(min(Mod(H_eig)), 1e-15)
      if (!is.finite(cond_num) || cond_num > 1e15) {
        warning(sprintf("SE: H_mat ill-conditioned (cond=%.2e); se set to NA", cond_num))
      } else {
        H_inv <- tryCatch(solve(H_mat),
                          error = function(e) solve(H_mat + 1e-8 * diag(param_dim)))
        psi <- -(H_inv %*% Q_mat)                                          # k x n
        cov.hat <- var(t(psi)) / n
        se  <- sqrt(diag(cov.hat))
      }
    }, error = function(e) {
      warning("SE computation failed: ", conditionMessage(e))
    })
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
    psi2 = psi,  # backward compatibility alias
    opt = opt
  )

  return(result)
}
