# ======================================================================
# FILE 1 / 3: ARLSA setup
# - Global settings / utility functions / DGP
# - rho estimation / benchmark estimators
# - ARLSA Algorithm 1
# - run_one_seed()
#
# 이 파일은 설정과 함수만 정의합니다.
# Monte Carlo 반복 자체는 실행하지 않습니다.
# ======================================================================

# ----------------------------------------------------------------------
# 0. Packages / global settings
# ----------------------------------------------------------------------

required_pkgs <- c("parallel", "doParallel", "foreach")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg, "' is required. Install it first with: ",
        "install.packages('", pkg, "')"
      )
    )
  }
}

library(parallel)
library(doParallel)
library(foreach)

# ---------------- USER SETTINGS ----------------
total_seeds <- 300

n_source <- 500
n_target <- 500

lambda_rho <- 0.001      # Tikhonov penalty inside rho estimation
lambda_arlsa <- 0.001    # Tikhonov penalty in ARLSA Algorithm 1, lambda_rho랑 같은 parameter는 아님 

rho_floor <- 0.01 #너무 작은 rho방지
rho_tilde_cap <- 20 #rho tilda의 극단값 제한
rho_soft_cap <- 10 # rho hat의 극단값을 부드럽게 -> numerical setting

d_floor <- 1e-4          # numerical positivity floor for D(X) projection
projection_ridge <- 1e-6 #hat matrix 보정
final_ridge <- 1e-6 #theta hat에서 inverse를 numerical하게 불안정 보정

use_parallel <- TRUE #병렬계산(true)
max_cores <- 8 #8코어
# ------------------------------------------------


# ----------------------------------------------------------------------
# 1. Basic utilities
# ----------------------------------------------------------------------

K_gauss <- function(d, bw) {
  dnorm(d, mean = 0, sd = bw) # d:두 값사이의 거리, bw:bandwidth
}

true_rho_func <- function(y) { #simulation 평가용 truth
  # q_Y(y) / p_Y(y)
  dnorm(y, mean = 1, sd = 1) /
    dnorm(y, mean = 0, sd = sqrt(2))
}

row_normalize <- function(M, eps = 1e-12) { #kernel 값들을 확률가중치처럼 만듬
  rs <- rowSums(M) #행별 합을 1로 만듬
  rs[!is.finite(rs) | rs < eps] <- eps
  M / rs
}

safe_linear_solve <- function(A, B, ridge = 0) {
  A <- as.matrix(A)
  B <- as.matrix(B)

  if (ridge > 0) {
    A <- A + ridge * diag(nrow(A)) #diag(nrow(A)):identity matrix
  }

  out <- tryCatch( #역행렬계산이 에러가 생겨도 계속 실행
    solve(A, B), #AX=B에서 matrix X 계산
    error = function(e) NULL
  )

  if (is.null(out) || any(!is.finite(out))) {
    return(NULL)
  }

  out
}


# ----------------------------------------------------------------------
# 2. Data generation; label shift상황 생성
# ----------------------------------------------------------------------

generate_data <- function(
    seed,
    n1 = n_source,
    n0 = n_target,
    alpha = c(-0.5, 0.5, 1)
) {
  set.seed(seed)

  N <- n1 + n0

  # Source and target Y differ -> label shift
  y_s <- rnorm(n1, mean = 0, sd = sqrt(2))
  y_t_hidden <- rnorm(n0, mean = 1, sd = 1) #target에서 Y 숨기기

  Y_all_hidden <- c(y_s, y_t_hidden)

  # Same X | Y in source and target
  X_all <- matrix(0, nrow = N, ncol = 3)

  X_all[, 1] <- alpha[1] * Y_all_hidden + rnorm(N, 0, 1)
  X_all[, 2] <- alpha[2] * Y_all_hidden + rnorm(N, 0, 1)
  X_all[, 3] <- alpha[3] * Y_all_hidden + rnorm(N, 0, 1)

  x_s <- X_all[seq_len(n1), , drop = FALSE]
  x_t <- X_all[n1 + seq_len(n0), , drop = FALSE]

  list(
    y_s = y_s,
    y_t_hidden = y_t_hidden,
    x_s = x_s,
    x_t = x_t,
    x_all = X_all,
    n1 = n1,
    n0 = n0,
    N = N,
    pi = n1 / N
  )
}


# Theoretical target OLS parameter:
# Y_T ~ N(1,1), X = alpha Y + eps
true_target_beta <- function(alpha = c(-0.5, 0.5, 1)) { #theta0 계산
  var_y <- 1
  mu_y <- 1

  Sigma_x <- diag(length(alpha)) + var_y * tcrossprod(alpha)
  cov_xy <- var_y * alpha

  slopes <- as.numeric(solve(Sigma_x, cov_xy))
  intercept <- mu_y - sum(slopes * alpha * mu_y)

  c(
    b0 = intercept,
    b1 = slopes[1],
    b2 = slopes[2],
    b3 = slopes[3]
  )
}


# ----------------------------------------------------------------------
# 3. Efficient Inference: progressive rho estimation
# ----------------------------------------------------------------------

build_ml_weight_matrix <- function(X_all, X_source, h_ml) { #Nadaraya-Watson으로 weight matrix
  N <- nrow(X_all)
  n1 <- nrow(X_source)

  W <- matrix(1, nrow = N, ncol = n1)

  for (k in seq_len(ncol(X_all))) {
    Dk <- outer(
      X_all[, k],
      X_source[, k],
      "-"
    )

    W <- W * K_gauss(Dk, h_ml)
  }

  row_normalize(W)
}


estimate_rho_efficient <- function(
    y_s,
    X_all,
    n1,
    n0,
    pi_val,
    lambda_reg = lambda_rho,
    y_grid = seq(-2.5, 2.5, length.out = 50),
    rho_floor_local = rho_floor,
    rho_tilde_cap_local = rho_tilde_cap,
    rho_soft_cap_local = rho_soft_cap
) {

  N <- n1 + n0
  X_source <- X_all[seq_len(n1), , drop = FALSE]

  # Bandwidths from Efficient Inference simulation section
  h_ml <- 3 * n1^(-1/7) #E_p( |x)
  l_bw <- 1.5 * n1^(-1/3) #E( |Y)
  h_bw <- 0.5 * n1^(-1/16) #K_h(Y-y0)

  # E_s(. | X) NW weights
  W_ML <- build_ml_weight_matrix(
    X_all = X_all,
    X_source = X_source,
    h_ml = h_ml
  )

  # 바깥 E(. | Y=y) NW weights: Nadaraya-Watson 근사
  W_Y_unnorm <- outer(
    y_s,
    y_s,
    function(y1, y2) K_gauss(y1 - y2, l_bw)
  )
  W_Y <- row_normalize(W_Y_unnorm)

  # K_h(Y_i - y0) over a grid
  V_mat_grid <- outer(
    y_s,
    y_grid,
    function(y_i, y_0) K_gauss(y_i - y_0, h_bw)
  )

  # Misspecified working model rho_star used in the Efficient Inference paper:휴리스틱하게
  rho_star_raw_s <- true_rho_func(y_s) *
    exp(0.2 * y_s + 0.1 * y_s^2)

  c_star <- 1 / mean(rho_star_raw_s)

  rho_star_grid <- c_star *
    true_rho_func(y_grid) *
    exp(0.2 * y_grid + 0.1 * y_grid^2)

  solve_rho_grid <- function(rho_guess_grid) { #rho star->rho tilde->rho hat

    rho_guess_s <- approx(
      x = y_grid,
      y = rho_guess_grid,
      xout = y_s,
      rule = 2
    )$y

    # rho는 positivity 보정
    rho_guess_s[!is.finite(rho_guess_s)] <- 1
    rho_guess_s <- pmax(rho_floor_local, rho_guess_s)

    denom_vec <- as.vector(
      W_ML %*% (
        rho_guess_s^2 +
          (pi_val / (1 - pi_val)) * rho_guess_s #rho_guess : rho(Y_i)
      )
    )

    denom_vec[!is.finite(denom_vec) | denom_vec < 1e-10] <- 1e-10
    w <- 1 / denom_vec

    # Discretized Fredholm operator
    M <- W_Y %*%
      diag(w[seq_len(n1)]) %*%
      W_ML[seq_len(n1), , drop = FALSE] %*%
      diag(rho_guess_s)

    lhs <- crossprod(M) + lambda_reg * diag(n1) #500X500
    rhs <- crossprod(M, V_mat_grid) #500X50

    A_mat <- safe_linear_solve(lhs, rhs) #500X50

    if (is.null(A_mat)) {
      stop("rho estimation failed at the Tikhonov solve.")
    }

    # E_s[a(Y) rho(Y) | X]
    E_term <- W_ML %*% #1000X50
      diag(rho_guess_s) %*%
      A_mat

    # Target direct term
    term1 <- colMeans(
      w[n1 + seq_len(n0)] *
        E_term[n1 + seq_len(n0), , drop = FALSE] #target Y는 사용X
    )

    # Source rectifier
    term2 <- colMeans(
      rho_guess_s *
        (
          V_mat_grid -
            w[seq_len(n1)] *
            E_term[seq_len(n1), , drop = FALSE] #target X의 정보+source data
        )
    )

    p_hat <- colMeans(V_mat_grid)

    (term1 + term2) / (p_hat + 1e-8)
  }

  # Stage 1: consistent rho_tilde
  rho_tilde_grid <- solve_rho_grid(rho_star_grid) #rho star -> rho tilda
  rho_tilde_grid[!is.finite(rho_tilde_grid)] <- 1 #0.01<=rhotilda<=20
  rho_tilde_grid <- pmax(
    rho_floor_local,
    pmin(rho_tilde_grid, rho_tilde_cap_local)
  )

  # Stage 2: refined / efficient rho_hat
  rho_hat_grid <- solve_rho_grid(rho_tilde_grid) #rho tilda -> rho hat
  rho_hat_grid[!is.finite(rho_hat_grid)] <- 1

  # Evaluate rho_hat at source labels
  rho_raw <- approx( #rho_raw:500
    x = y_grid,
    y = rho_hat_grid,
    xout = y_s,
    rule = 2
  )$y

  rho_raw[!is.finite(rho_raw)] <- 1
  rho_positive <- pmax(rho_floor_local, rho_raw) #최소 0.01

  # Numerical control:
  # smooth cap + renormalization, finite-sample safeguard : numerical stabilization을 위해서 논문X
  rho_soft <- rho_soft_cap_local *
    tanh(rho_positive / rho_soft_cap_local)

  rho_soft <- pmax(rho_floor_local, rho_soft)

  rho_final <- rho_soft / mean(rho_soft) #rho_final=(rho1 hat,...,rho500 hat)'

  #최종 output
  list(
    rho = rho_final,
    rho_raw = rho_raw,
    rho_star_grid = rho_star_grid,
    rho_tilde_grid = rho_tilde_grid,
    rho_hat_grid = rho_hat_grid,
    y_grid = y_grid,
    h_ml = h_ml,
    l_bw = l_bw,
    h_bw = h_bw
  )
}


# ----------------------------------------------------------------------
# 4. Benchmark estimators
# ----------------------------------------------------------------------

#label shift 보장X
fit_naive <- function(y_s, x_s) { 
  as.numeric(coef(lm(y_s ~ x_s)))
}

#IPW
fit_ipw <- function(y_s, x_s, rho) {
  as.numeric(coef(lm(y_s ~ x_s, weights = rho)))
}

#유니콘
fit_oracle <- function(y_t_hidden, x_t) {
  as.numeric(coef(lm(y_t_hidden ~ x_t)))
}


# ----------------------------------------------------------------------
# 5. ARLSA Algorithm 1
# ----------------------------------------------------------------------

make_U1_U2 <- function(y_s, Phi) { #phi:500X4, phi<-cbind(1, x_s)
  n <- nrow(Phi)
  p <- ncol(Phi)

  # OLS score:
  # U(y,x,theta) = x*y - x*x^T theta
  U1 <- Phi * y_s #500X4

  U2 <- matrix(0, nrow = n, ncol = p^2) #500X16, p=4

  for (i in seq_len(n)) {
    xxT <- tcrossprod(Phi[i, ])
    U2[i, ] <- as.vector(xxT)
  }

  list(U1 = U1, U2 = U2, p = p)
}


fit_arlsa <- function(
    y_s, #500X1
    x_s, #500X3
    x_all, #1000X3
    rho, #500
    pi_val, #pi=0.5
    lambda_reg = lambda_arlsa, #lambda_arlsa=0.001
    d_floor_local = d_floor,
    projection_ridge_local = projection_ridge,
    final_ridge_local = final_ridge
) {

  n1 <- length(y_s) #n1=500
  N <- nrow(x_all) #N=1000

  # Uploaded reference code uses this simple linear basis.
  # This is also an allowed simple choice in the ARLSA draft.
  Phi <- cbind(1, x_s) #500X4
  Phi_all <- cbind(1, x_all) #1000X4

  p_dim <- ncol(Phi)

  # Score matrices U1, U2
  score_obj <- make_U1_U2(y_s, Phi)
  U1 <- score_obj$U1
  U2 <- score_obj$U2

  # ----------------------------------------------------------
  # Step 1: source pre-computation
  # ----------------------------------------------------------
  XtX <- crossprod(Phi) #4X4

  G <- safe_linear_solve( #4X4
    XtX,
    diag(ncol(Phi)),
    ridge = projection_ridge_local
  )

  if (is.null(G)) {
    return(rep(NA_real_, p_dim))
  }

  H <- Phi %*% G %*% t(Phi) #500X500
  H_all <- Phi_all %*% G %*% t(Phi) #1000X500

  # NW kernel matrix K_h over source Y
  h_y <- 1.06 * sd(y_s) * n1^(-1/5)

  dist_y <- as.matrix(dist(y_s)) #500X500
  K_h <- exp(-(dist_y^2) / (2 * h_y^2)) #Gaussian kernel
  K_h <- row_normalize(K_h) #500X500

  P <- diag(rho) #500X500
  P2 <- P %*% P #500X500

  # d = H(P rho + pi/(1-pi) rho)
  d <- as.vector(
    H %*% (
      P %*% rho +
        (pi_val / (1 - pi_val)) * rho
    )
  )

  d_all <- as.vector(
    H_all %*% (
      P %*% rho +
        (pi_val / (1 - pi_val)) * rho
    )
  )

  # Same finite-sample positivity protection as the uploaded reference code
  d_safe <- pmax(d, d_floor_local) #numerical safeguard
  d_all_safe <- pmax(d_all, d_floor_local)

  W <- diag(1 / d_safe) #500X500
  W_all <- diag(1 / d_all_safe) #1000X1000

  # ----------------------------------------------------------
  # Step 2: nuisance functions A1, A2
  # ----------------------------------------------------------
  L_h <- K_h %*% W %*% H %*% P #500X500
                               #rho -> E_s( |x) -> D(x)^-1 -> Y방향 smoothing

  # Common Tikhonov operator
  lhs_L <- crossprod(L_h) + lambda_reg * diag(n1) #500X500

  M_h_lambda <- safe_linear_solve( #500X500
    lhs_L,
    t(L_h)
  )

  if (is.null(M_h_lambda)) {
    return(rep(NA_real_, p_dim))
  }

  R1 <- K_h %*% #500X4
    (diag(n1) - W %*% H %*% P2) %*%
    U1

  R2 <- K_h %*% #500X16
    (diag(n1) - W %*% H %*% P2) %*%
    U2

  A1 <- M_h_lambda %*% R1 #500X4
  A2 <- M_h_lambda %*% R2 #500X16

  # ----------------------------------------------------------
  # Step 3: global projection B1, B2; B_k=source+target
  # ----------------------------------------------------------
  B1 <- W_all %*%
    H_all %*%
    (P2 %*% U1 + P %*% A1) #1000X4

  B2 <- W_all %*%
    H_all %*%
    (P2 %*% U2 + P %*% A2) #1000X16

  # ----------------------------------------------------------
  # Step 4: Psi1, Psi2 and closed-form theta
  # ----------------------------------------------------------
  n0 <- N - n1

  Psi1 <- matrix(0, nrow = N, ncol = p_dim) #1000X4
  Psi2_vec <- matrix(0, nrow = N, ncol = p_dim^2) #1000X16

  # Source r=1
  Psi1[seq_len(n1), ] <-
    (1 / pi_val) *
    rho *
    (U1 - B1[seq_len(n1), , drop = FALSE])

  Psi2_vec[seq_len(n1), ] <-
    (1 / pi_val) *
    rho *
    (U2 - B2[seq_len(n1), , drop = FALSE])

  # Target r=0; target Y never appears
  target_idx <- n1 + seq_len(n0)

  Psi1[target_idx, ] <-
    (1 / (1 - pi_val)) *
    B1[target_idx, , drop = FALSE]

  Psi2_vec[target_idx, ] <-
    (1 / (1 - pi_val)) *
    B2[target_idx, , drop = FALSE]

  S1 <- colSums(Psi1) #4X1

  # vec ordering is consistent with as.vector(xxT)
  S2 <- matrix( #4X4
    colSums(Psi2_vec),
    nrow = p_dim,
    ncol = p_dim
  )

  theta_hat <- safe_linear_solve(
    S2,
    S1,
    ridge = final_ridge_local
  )

  if (is.null(theta_hat)) {
    return(rep(NA_real_, p_dim))
  }

  as.numeric(theta_hat)
}


# ----------------------------------------------------------------------
# 6. One Monte Carlo replication: theta naive,ipw,arlsa,oracle 
# ----------------------------------------------------------------------

run_one_seed <- function(seed, return_full = FALSE) {

  dat <- generate_data(
    seed = seed,
    n1 = n_source,
    n0 = n_target
  )

  rho_obj <- estimate_rho_efficient(
    y_s = dat$y_s,
    X_all = dat$x_all,
    n1 = dat$n1,
    n0 = dat$n0,
    pi_val = dat$pi
  )

  rho <- rho_obj$rho

  theta_naive <- fit_naive(
    dat$y_s,
    dat$x_s
  )

  theta_ipw <- fit_ipw(
    dat$y_s,
    dat$x_s,
    rho
  )

  theta_arlsa <- fit_arlsa(
    y_s = dat$y_s,
    x_s = dat$x_s,
    x_all = dat$x_all,
    rho = rho,
    pi_val = dat$pi
  )

  theta_oracle <- fit_oracle(
    dat$y_t_hidden,
    dat$x_t
  )

  rho_true_s <- true_rho_func(dat$y_s)

  result <- data.frame(
    Seed = seed,

    b0_Naive = theta_naive[1],
    b0_IPW = theta_ipw[1],
    b0_ARLSA = theta_arlsa[1],
    b0_Oracle = theta_oracle[1],

    b1_Naive = theta_naive[2],
    b1_IPW = theta_ipw[2],
    b1_ARLSA = theta_arlsa[2],
    b1_Oracle = theta_oracle[2],

    b2_Naive = theta_naive[3],
    b2_IPW = theta_ipw[3],
    b2_ARLSA = theta_arlsa[3],
    b2_Oracle = theta_oracle[3],

    b3_Naive = theta_naive[4],
    b3_IPW = theta_ipw[4],
    b3_ARLSA = theta_arlsa[4],
    b3_Oracle = theta_oracle[4],

    Rho_RMSE = sqrt(mean((rho - rho_true_s)^2)),
    Rho_Min_Raw = min(rho_obj$rho_raw, na.rm = TRUE),
    Rho_Max_Raw = max(rho_obj$rho_raw, na.rm = TRUE),
    Rho_Min_Final = min(rho, na.rm = TRUE),
    Rho_Max_Final = max(rho, na.rm = TRUE)
  )

  if (!return_full) {
    return(result)
  }

  list(
    result = result,
    data = dat,
    rho = rho_obj
  )
}
