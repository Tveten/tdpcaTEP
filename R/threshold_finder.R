#' Estimating the threshold for tpca changepoint detection
#'
#' Description
#'
#' \code{n} and \code{alpha} governs the false alarm rate by the relation
#'   P(T < n | H_0) <= \code{alpha}.
#' The corresponding average run length (\code{arl}) is approximately given by n / alpha.
#'
#' \code{rel_tol} and \code{thresh_alpha} governs the number of simulations used in
#' each step of the algorithm towards a more and more certain estimate.
#' At each step, the number of simulations is chosen so that
#' \code{[(1 - rel_tol)arl, (1 + rel_tol)arl]} approximately covers the true
#' average run length at confidence level \code{thresh_alpha}.
#' For example, when the last \code{rel_tol} is 0.025, it means that the final
#' estimated threshold corresponds to an average run length of approximately
#' \code{arl} +- 0.025 * \code{arl} at confidence level \code{thresh_alpha}.
#' The algorithm should start with a large relative error tolerance,
#' and then narrow it down for the quickest convergence.
#'
#' @param x_train The d x m training data matrix, where d is the dimension of the data
#'   and m the number of training samples.
#' @param lag The number of lags in TDPCA and the block bootstrap.
#' @param n The length of the segment to monitor for false alarms. See details.
#' @param alpha Probability of type I error (false alarm) within the time window . See details.
#' @param axes Indices of the principal axes to be used in simulations.
#'   MUST be supplied if mon_type == 'tpca'.
#' @param p0 The mixture probability of a dimension being affected by a change.
#'   MUST be supplied if mon_type == 'mixture'.
#' @param w The window size (an integer).
#' @param rel_tol A vector with the sequence of relative error tolerances
#'   allowed at each step towards convergence in the algorithm. See details.
#' @param thresh_alpha 1 - \code{thresh_alpha} is the confidence level. See details.
#' @param init_thresh Use if a custom initial value for the threshold is wanted.
#' @param learning_coef The learning rate is defined as \code{rel_tol} /
#' \code{learning_coef}. The default \code{learning_coef} is 3, which has
#' been found a good choice after a lot of experimenting.
#' #' @param file_id A string that identifies the files associated with the particular
#' threshold beyond the data dimension d, axes/p0, training size m or
#' probability of false alarm alpha.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{threshold}}{The final estimate }
#' }
#' \code{tpca_threshold} also creates a .txt file with each line showing
#' summaries of each step in the estimation procedure. The last line corresponds
#' to the final estimate.
#'
#' @export

threshold_finder_tdpca <- function(x_train, lag, n, alpha,
                             axes          = NULL,
                             w             = 200,
                             rel_tol       = c(0.2, 0.1, 0.05, 0.025),
                             thresh_alpha  = 0.05,
                             init_thresh   = NULL,
                             learning_coef = NULL,
                             file_id       = NULL) {

  set_file_name <- function(file_type) {
    get_n_equal_files <- function(root_name, files_in_wdir) {
      split_files <- strsplit(files_in_wdir, '[(.]')
      stripped_files <- vapply(split_files, `[`, character(1), 1)
      n_equals <- sum(root_name == stripped_files)
      n_equals
    }

    files_in_wdir <- list.files()
    dir <- './thresholds/threshold_results/'
    root_name <- paste0('tpca_threshold_', file_id, file_type, '_',
                        'n', n, 'alpha', substr_right(alpha, 2),
                        'm', m, 'd', d,
                        'ax', paste(axes, collapse = '-'))

    n_equal_files <- get_n_equal_files(root_name, files_in_wdir)
    if (n_equal_files > 0)
      root_name <- paste0(root_name, '(', as.character(n_equal_files + 1), ')')
    paste0(dir, root_name, '.txt')
  }

  set_log_name <- function() {
    set_file_name('log')
  }

  set_results_name <- function() {
    set_file_name('results')
  }

  log_next_stage <- function(rel_tol) {
    write(sprintf('=========================================================='),
          file = log_file, append = TRUE)
    write(sprintf('Relative error tolerance: %.3f',
                  rel_tol),
          file = log_file, append = TRUE)
    write(sprintf(' '),
          file = log_file, append = TRUE)
  }

  log_current <- function(n_sim, threshold) {
    write(sprintf('Currently: threshold = %.4f with %.0f simulations',
                  threshold, n_sim),
          file = log_file, append = TRUE)
  }

  log_results <- function(arl_est, threshold) {
    curr_time <- proc.time()[3] / 60
    write(sprintf('Result:    ARL       = %.0f, %.1f min used.',
                  arl_est, curr_time - start_time),
          file = log_file, append = TRUE)
  }

  init_log_file <- function() {
    write(sprintf('Data dimension             = %.0f', d), file = log_file, append = TRUE)
    write(sprintf('Reduced dimension          = %.0f', r), file = log_file, append = TRUE)
    write(paste0('Axes                       = ', paste(axes, collapse = ' ')), file = log_file, append = TRUE)
    write(sprintf('Number of training samples = %.0f', m), file = log_file, append = TRUE)
    write(sprintf('Window length              = %.0f', w), file = log_file, append = TRUE)
  }

  init_results_file <- function() {
    header_items <- 'data_dim'
    header_items <- c(header_items, 'axes')
    header_items <- c(header_items,  'n_training',
                      'window_size', 'n_sim', 'threshold',
                      'arl_est', 'arl_lower', 'arl_upper')
    write(paste(header_items, collapse = ' '),
          file = results_file, append = TRUE, ncolumns = length(header_items))
  }

  store_results <- function(n_sim, threshold, arl_est, arl_conf_int) {
    stored_values <- d
    axes_string <- paste(axes, collapse = '-')
    stored_values <- c(stored_values, axes_string)
    stored_values <- c(stored_values, m, w, n_sim, round(threshold, 5),
                       arl_est, arl_conf_int)
    write(stored_values, file = results_file,
          append = TRUE, ncolumns = length(stored_values))
  }

  get_n_sim <- function(thresh_alpha, rel_tol) {
    z <- qnorm(1 - thresh_alpha / 2)
    n_sim <- (z / rel_tol)^2
    return(round(n_sim))
  }

  init_threshold <- function(r) {
    if (is.null(init_thresh)) return(9 + 1.2 * r)
    else return(init_thresh)
  }

  update_threshold <- function(threshold, rel_tol, arl, arl_est) {
    if (is.null(learning_coef)) learning_rate <- rel_tol / 3
    else learning_rate <- rel_tol / learning_coef
    threshold_change <- (arl - arl_est) / arl * learning_rate
    threshold * (1 + threshold_change)
  }

  get_arl_func <- function() {
    function(threshold, n_sim)
      tdpca_arl(threshold, x_train, lag, axes, n, w, n_sim)
  }

  # TODO: Error handling: Either axes or p0 must be supplied.
  # TODO: Update log and result functions.

  m <- ncol(x_train) - lag
  d <- nrow(x_train) * (lag + 1)

  r <- length(axes)
  n_sim <- get_n_sim(thresh_alpha, rel_tol)
  threshold <- init_threshold(r)
  arl <- n / alpha
  arl_func <- get_arl_func()

  log_file <- set_log_name()
  results_file <- set_results_name()
  init_log_file()
  init_results_file()
  start_time <- proc.time()[3] / 60

  # When rel_tol is low, the algorithm might jump back and forth over the true
  # value. Thus, it should be stopped at some point. Too many runs on the final
  # stage is not worth it.
  max_retries <- max(8, round(20 - 1 / 4 * r))
  n_retries <- 1
  for (k in 1:length(n_sim)) {
    arl_conf_int <- c(0, 0)
    log_next_stage(rel_tol[k])
    while (!is_in_interval(arl, arl_conf_int) && (n_retries <= max_retries)) {
      log_current(n_sim[k], threshold)

      arl_est <- arl_func(threshold, n_sim[k])
      arl_conf_int <- geom_conf_int(arl_est, n_sim[k], thresh_alpha)

      log_results(arl_est, threshold)
      store_results(n_sim[k], threshold, arl_est, arl_conf_int)

      if (!is_in_interval(arl, arl_conf_int))
        threshold <- update_threshold(threshold, rel_tol[k], arl, arl_est)
      if (k == length(n_sim)) n_retries <- n_retries + 1
    }
  }
  threshold
}

conf_int <- function(x, alpha = 0.05) {
  ci <- quantile(x, probs = c(alpha / 2, (1 - alpha) / 2))
}

geom_conf_int <- function(mean_est, n, thresh_alpha) {
  p_est <- 1 / mean_est
  sd_est <- sqrt((1 - p_est) / p_est^2)
  percentile <- qnorm((1 - thresh_alpha / 2))
  lower <- mean_est - percentile * sd_est / sqrt(n)
  upper <- mean_est + percentile * sd_est / sqrt(n)
  conf_int <- c(floor(lower), ceiling(upper))
  return(conf_int)
}
