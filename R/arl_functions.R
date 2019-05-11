#' Estimate average run length (ARL) for monitoring by tdpca
#'
#' @param threshold A numeric specifying the threshold value for when a change
#'   is declared.
#' @param x_train The training data (not lag-extended).
#' @param lag The number of lags used in TDPCA.
#' @param axes A vector specifying which axes/projections to use.
#' @param n The number of observations to monitor for an estimate of the ARL.
#'   See details.
#' @param w The window size. Number of recent time-points to consider for a
#'   change.
#' @param n_sim The number of simulations to base the estimate on.
#'
#' @return An estimate of the average run length.
#'
#' @export

tdpca_arl <- function(threshold, x_train, lag, axes, n, w, n_sim) {
  m <- ncol(x_train)

  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  run_lengths <- foreach::foreach(b = 1:n_sim, .combine = 'c') %dopar% {
    boot_train_obj <- boot_z_train_np(m, x_train, lag, axes)
    z_boot <- boot_z_mon_np(n, x_train, boot_train_obj, lag, axes)
    mixture_rl(threshold, boot_train_obj$z, z_boot, 1, w)
  }
  stop_parallel(comp_cluster)
  est_arl(run_lengths, n - lag, n_sim)
}

boot_z_train_np <- function(n_boot, x_train, lag, axes) {
  x_boot <- block_boot(n_boot, x_train, lag)
  x_boot_e <- lag_extend(x_boot, lag)
  mu_hat <- rowMeans(x_boot_e)
  sigma_hat <- rowSds(x_boot_e)
  x_boot_e <- (x_boot_e - mu_hat) / sigma_hat

  m <- ncol(x_boot_e)
  cor_mat_hat <- 1 / (m - 1) * x_boot_e %*% t(x_boot_e)
  pca_obj <- tpca::pca(cor_mat_hat, axes = axes)
  V <- pca_obj$vectors
  lambda <- pca_obj$values

  z_boot <- V %*% x_boot_e / sqrt(lambda)
  list('z'= z_boot, 'sigma' = sigma_hat, 'mu' = mu_hat, 'V' = V, 'lambda' = lambda)
}

boot_z_mon_np <- function(n_boot, x_train, boot_train_obj, lag, axes) {
  x_boot <- block_boot(n_boot, x_train, lag)
  x_boot_e <- lag_extend(x_boot, lag)

  mu_hat_train <- boot_train_obj$mu
  sigma_hat_train <- boot_train_obj$sigma
  V_train <- boot_train_obj$V
  lambda_train <- boot_train_obj$lambda

  x_boot_e_norm <- (x_boot_e - mu_hat_train) / sigma_hat_train
  z_boot <- V_train %*% x_boot_e_norm / sqrt(lambda_train)
  z_boot
}


block_boot <- function(n_boot, x, lag) {
  l <- lag + 1
  b <- ceiling(n_boot/l)
  ind_vec <- 1:(n_boot - l)
  block_ind <- rep(0, b)
  block_ind[1] <- sample(ind_vec, 1)
  for (j in 2:b) {
    prev <- block_ind[j - 1]
    avail_ind <- ind_vec[-min(prev + l - 1, n_boot - l)]
    block_ind[j] <- sample(avail_ind, 1)
  }
  x_boot <- lapply(block_ind, function(i) x[, i:(i + l - 1)])
  x_boot <- do.call('cbind', x_boot)
  x_boot[, 1:n_boot]
}

lag_extend <- function(x, lag) {
  # Matrix: (x_{t - lag}, ...,  x_{t})^T, I.e., time t is in the lower rows.
  n <- ncol(x)
  p <- nrow(x)
  x_extended <- matrix(0, nrow = p * (lag + 1), ncol = n - lag)
  vapply((lag + 1):n, function(i) as.vector(x[, (i - lag):i]),
         numeric(p * (lag + 1)))
}


est_arl <- function(run_lengths, n, n_sim) {
  arl_est <- n / mean(as.numeric(run_lengths < n))
  if (is.infinite(arl_est)) arl_est <- 5 * n / (1 / n_sim)
  round(arl_est)
}
