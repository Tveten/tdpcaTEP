find_dpca_thresholds <- function(dpca_model, n, alpha, type = 'val') {
  load_TEP_data_globally(no_test = TRUE)
  if (type == 'val') {
    x <- get_TEP_train(101:150, return_all = FALSE)
    x <- lapply(x, lag_extend, lag = dpca_model$lag)
    x <- do.call('cbind', x)
  } else if (type == 'train') {
    x <- dpca_model$x_extended
  } else stop('Invalid "type" argument. Must be "val" or "train".')

  T2 <- calc_T2(x, dpca_model)
  Q <- calc_Q(x, dpca_model)
  ind <- ceiling(alpha / (2 * n))
  thresholds <- list('T2' = sort(T2, decreasing = TRUE)[ind],
                     'Q'   = sort(Q, decreasing = TRUE)[ind])

  thresholds
}

get_analytical_dpca_thresholds <- function(dpca_model, n, alpha) {
  k <- length(dpca_model$lambda)
  alpha_n <- alpha / (2 * n)
  T2_thresh <- qchisq(1 - alpha_n, k)

  theta <- vapply(1:3, function(i) sum(dpca_model$lambda_rest^i), numeric(1))
  h0 <- 1 - 2 * theta[1] * theta[3] / (3 * theta[2]^2)
  c <- qnorm(1 - alpha_n)
  Q_thresh <- theta[1] / h0 * (c * sqrt(2 * theta[2] * h0^2) / theta[1] +
                               1 + theta[2] * h0 * (h0 - 1) / theta[1]^2)

  list('T2' = T2_thresh, 'Q' = Q_thresh)
}

train_dpca_model <- function(x_train, lag, cpv = 0.95) {
  x_train_e <- lag_extend(x_train, lag)
  d <- nrow(x_train_e)
  m <- ncol(x_train_e)

  # Training:
  mu_x <- rowMeans(x_train_e)
  sigma_x <- rowSds(x_train_e)
  x_train_e <- (x_train_e - mu_x) / sigma_x
  cor_mat_x <- 1 / (m - 1) * x_train_e %*% t(x_train_e)
  pca_obj <- tpca::pca(cor_mat_x)
  V <- pca_obj$vectors
  lambda <- pca_obj$values

  cumprop_lambda <- cumsum(lambda) / sum(lambda)
  k <- min(which(cumprop_lambda >= cpv))
  V_k <- V[1:k, ]
  lambda_k <- lambda[1:k]

  list('x_extended'  = sigma_x * x_train_e + mu_x,
       'cor_mat'     = cor_mat_x,
       'lag'         = lag,
       'V'           = V_k,
       'lambda'      = lambda_k,
       'mu'          = mu_x,
       'sigma'       = sigma_x,
       'lambda_rest' = lambda[(k + 1):d])
}

calc_T2 <- function(x, dpca_model) {
  V <- dpca_model$V
  lambda <- dpca_model$lambda
  mu_x <- dpca_model$mu
  sigma_x <- dpca_model$sigma
  x <- (x - mu_x) / sigma_x

  y_mat <- V %*% x
  diag_lambda_inv <- diag(1 / lambda)
  T2 <- vapply(1:ncol(y_mat), function(i) {
    t(y_mat[, i]) %*% diag_lambda_inv %*% y_mat[, i]
  }, numeric(1))
  T2
}

calc_Q <- function(x, dpca_model) {
  V <- dpca_model$V
  mu_x <- dpca_model$mu
  sigma_x <- dpca_model$sigma
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  d <- nrow(x)
  x <- (x - mu_x) / sigma_x

  pred_mat <- t(V) %*% V
  Id <- diag(rep(1, d))
  Q <- vapply(1:ncol(x), function(i) {
    t(x[, i]) %*% (Id - pred_mat) %*% x[, i]
  }, numeric(1))
  Q
}

#' Run the DPCA method on the TEP data
#'
#' This function reproduces the run length results of DPCA on the TEP data.
#' Results are stored in results/TEP_rl_dpca_<threshold>.txt.
#'
#' @param fault_nr A numeric vector of desired fault numbers between 1 and 20.
#' @param threshold 'analytical' or 'val' to determine whether the analytical
#' thresholds for DPCA should be used, ot thresholds determined by a validation
#' set.
#'
#' @export
run_TEP_dpca <- function(fault_nr = 1:20, threshold = 'analytical') {
  lag <- 5
  load_TEP_data_globally()
  training_sets <- get_TEP_train(2, return_all = FALSE)[[1]]
  dpca_model <- train_dpca_model(training_sets, lag)

  kappa <- 160 - lag
  alpha <- 0.01
  if (threshold == 'analytical')
    thresholds <- get_analytical_dpca_thresholds(dpca_model, kappa, alpha)
  if (threshold == 'val')
    thresholds <- find_dpca_thresholds(dpca_model, kappa, alpha)
  print(thresholds)

  n_test_sets <- 500
  # run_lengths <- matrix(0, nrow = n_test_sets, ncol = length(fault_nr))
  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  run_lengths <- foreach::foreach(i = seq_along(fault_nr), .combine = 'cbind') %dopar% {
  # for (i in seq_along(fault_nr)) {
    test_sets <- get_TEP_test(1:n_test_sets, fault_nr[i], return_all = FALSE)
    test_sets <- lapply(test_sets, lag_extend, lag = dpca_model$lag)
    # run_lengths[, i] <- sim_rl_dpca(test_sets, dpca_model, thresholds, kappa, fault_nr[i])
    sim_rl_dpca(test_sets, dpca_model, thresholds, kappa, fault_nr[i])
  }

  colnames(run_lengths) <- paste0('Fault', fault_nr)
  file_name <- paste0('./results/TEP_rl_dpca_', threshold, '.txt')
  write.table(run_lengths, file_name, row.names = FALSE)
}

sim_rl_dpca <- function(test_sets, dpca_model, thresholds, kappa, fault_nr) {
  mon_length <- ncol(test_sets[[1]])

  run_lengths <- rep(0, length(test_sets))
  for (i in seq_along(test_sets)) {
    T2 <- calc_T2(test_sets[[i]], dpca_model)
    Q <- calc_Q(test_sets[[i]], dpca_model)
    T2_stop <- suppressWarnings(min(mon_length + 1, min(which(T2 > thresholds$T2))))
    Q_stop <- suppressWarnings(min(mon_length + 1, min(which(Q > thresholds$Q))))
    run_lengths[i] <- min(T2_stop, Q_stop)
  }

  run_lengths
}

est_edd_dpca <- function(lag = 5, threshold = 'analytical', show_plot = FALSE) {
  kappa <- 160 - lag
  file_name <- paste0('./results/TEP_rl_dpca_', threshold, '.txt')
  run_lengths <- read.table(file_name, sep = ' ', header = TRUE)

  if (show_plot) {
    par(mfrow = c(3, 3))
    lapply(1:9, function(i) {
      plot(density(run_lengths[[i]], adjust = 0.3))
      abline(v = kappa, lwd = 3)
    })
  }

  pfa_edd_list <- lapply(run_lengths, function(rl) {
    pfa <- sum(rl <= kappa) / length(rl)
    edd <- mean(rl[rl > kappa] - kappa)
    c(pfa, edd)
  })
  pfa_edd_df <- as.data.frame(do.call('rbind', pfa_edd_list))
  names(pfa_edd_df) <- c('pfa', 'edd')

  pfa_edd_df
}

plot_T_Q <- function(T2, Q, thresholds) {
  par(mfrow = c(2, 1))
  plot(T2, type = 'l')
  abline(h = thresholds$T2, col = 'red')
  plot(Q, type = 'l')
  abline(h = thresholds$Q, col = 'red')
}

test_dpca <- function(fault_nr) {
  n_test_sets <- 100
  lag <- 5
  load_TEP_data_globally()
  training_sets <- get_TEP_train(2, return_all = FALSE)[[1]]
  dpca_model <- train_dpca_model(training_sets, lag)

  kappa <- 160 - lag
  alpha <- 0.01
  thresholds <- get_analytical_dpca_thresholds(dpca_model, kappa, alpha)
  print(thresholds)

  test_sets <- get_TEP_test(1:n_test_sets, fault_nr, return_all = FALSE)
  test_sets <- lapply(test_sets, lag_extend, lag = dpca_model$lag)
  T2 <- calc_T2(test_sets[[1]], dpca_model)
  Q <- calc_Q(test_sets[[1]], dpca_model)
  plot_T_Q(T2, Q, thresholds)
  run_length <- sim_rl_dpca(test_sets, dpca_model, thresholds, kappa, fault_nr)
  print(run_length)
}
