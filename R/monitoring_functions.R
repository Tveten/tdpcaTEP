#' @useDynLib TEPExample
#' @importFrom Rcpp sourceCpp
NULL

mixture_rl <- function(threshold, x_train, x, p0, w) {
  m <- ncol(x_train)
  n_max <- ncol(x)

  t <- 1
  sums <- init_sums(x_train, x[, t], n_max)
  detection_stat <- 0
  while ((detection_stat < threshold) & (t < n_max)) {
    t <- t + 1
    sums <- update_sumsC(sums, x[, t], m, t)
    log_liks <- mixture_log_liksC(sums, m, t, w, p0)
    detection_stat <- max(log_liks)
  }
  t
}

init_sums <- function(x_train, x_1, n_max) {
  d <- nrow(x_train)
  m <- ncol(x_train)
  sums <- list('u' = matrix(0, nrow = d, ncol = n_max + 1),
               'v' = matrix(0, nrow = d, ncol = n_max + 1))
  sums$u[, 1] <- rowSums(x_train)
  sums$v[, 1] <- (m - 1) * rowVars(x_train)
  update_sumsC(sums, x_1, m, 1)
}

bartlett_corr <- function(m, t, k) {
  1 / 2 * (-(m + t) * log(m + t) +
            (m + k) * log(m + k) +
            (t - k) * log(t - k) +
            (m + t) * digamma((m + t - 1) / 2) -
            (m + k) * digamma((m + k - 1) / 2) -
            (t - k) * digamma((t - k - 1) / 2))
}

mixture_log_liksC <- function(sums, m, t, w, p0) {
  log_liks <- log_liks_matC(sums, m, t, w, p0)
  ks <- max(0, t - (w + 1)):(t - 2)
  c <- bartlett_corr(m, t, ks)
  corr_log_liks <- log_liks %*% diag(1 / c, nrow = length(c))
  global_log_liks <- colSums(log(1 - p0 + p0 * exp(1 / 2 * corr_log_liks)))
  global_log_liks
}
