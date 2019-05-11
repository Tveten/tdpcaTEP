plot_TEP_proj <- function(D = nrow(x_train), sim_nr = 2, fault_nr = 11) {
  x_train <- get_TEP_train(1, return_all = FALSE)
  x_train <- lag_extend(x_train, 3)
  x_train <- (x_train - rowMeans(x_train)) / rowSds(x_train)
  m <- ncol(x_train)
  Sigma0_hat <- 1 / (m - 1) * x_train %*% t(x_train)
  pca_obj <- tpca::pca(Sigma0_hat)
  V <- pca_obj$vectors
  lambda <- pca_obj$values
  z_train <- V %*% x_train / sqrt(lambda)

  x_test <- get_TEP_test(sim_nr, fault_nr, return_all = FALSE)[[1]]
  x_test <- lag_extend(x_test, 3)
  x_test <- (x_test - rowMeans(x_test)) / rowSds(x_test)
  z_test <- V %*% x_test / sqrt(lambda)

  kappa <- 160
  n_before <- 2 * kappa
  n_after <- 160
  for (d in D) {
    print(d)
    z_d <- c(z_train[d, ], z_test[d, ])
    par(mfrow = c(2, 1))
    plot(z_d[(m + kappa - n_before):(m + kappa + n_after)], type = 'l',
         main = paste0('Projection ', d))
    abline(v = n_before - kappa, lwd = 2)
    abline(v = n_before, lwd = 2, col = 'red')
    acf(z_train[d, ])
    Sys.sleep(2)
  }
}

get_rl_file <- function(lag, fault_nr, change_distr = 'semisparse_uniform') {
  dir <- './results/'
  n <- 160 - lag
  unlagged_dim <- 52
  d <- unlagged_dim * (lag + 1)
  alpha_str <- '01'
  m <- 500 - lag
  paste0(dir, 'TEP_rl_fault', fault_nr, '_', change_distr,
         '_n', n, 'alpha', alpha_str, 'm', m, 'd', d, 'l', lag, '.txt')
}

get_rl_table <- function(lag, fault_nr, keep_names = NULL,
                         change_distr = c('semisparse_uniform',
                                          'semisparse_sd_only',
                                          'semisparse_mean_only')) {
  n_distr <- length(change_distr)
  rl_file <- get_rl_file(lag, fault_nr, change_distr[1])
  rl_df <- read.table(rl_file, header = TRUE, sep = ' ')
  tpca_ind <- grepl('tpca', names(rl_df))
  names(rl_df)[tpca_ind] <-  paste0(names(rl_df)[tpca_ind], ' ', change_distr[1])
  if (n_distr > 1) {
    for (i in 2:n_distr) {
      rl_file <- get_rl_file(lag, fault_nr, change_distr[i])
      rl_df_temp <- read.table(rl_file, header = TRUE, sep = ' ')
      rl_df_temp <- rl_df_temp[, grepl('tpca', names(rl_df_temp))]
      names(rl_df_temp) <- paste0(names(rl_df_temp), ' ', change_distr[i])
      rl_df <- cbind(rl_df, rl_df_temp)
    }
  }

  dpca_rl <- read.table('./results/TEP_rl_dpca_analytical.txt', sep = ' ', header = TRUE)
  rl_df$dpca <- dpca_rl[, fault_nr]
  if (!is.null(keep_names)) rl_df <- rl_df[keep_names]

  rl_df
}

rl_hist <- function(lag, fault_nr,
                    change_distr = c('semisparse_uniform',
                                     'semisparse_sd_only',
                                     'semisparse_mean_only')) {
  n_distr <- length(change_distr)
  rl_df <- get_rl_table(lag, fault_nr, change_distr = change_distr)

  par(mfrow = c(2 + n_distr, 3))
  for (i in 1:length(rl_df)) {
    hist(rl_df[[i]], breaks = 20, probability = TRUE,
         xlab = 'run length', main = names(rl_df)[i])
  }
}

setup_legend <- function(method_names) {
  legend_info <- list('max_pca'              = c('Max DPCA', 'darkgoldenrod'),
                      'min_pca'              = c('Min DPCA', 'cornflowerblue'),
                      'semisparse_uniform'   = c('TDPCA(unif)', 'firebrick3'),
                      'semisparse_mean_only' = c('TDPCA(mean)', 'mediumpurple3'),
                      'semisparse_sd_only'   = c('TDPCA(var)', 'darkgreen'),
                      'dpca'                 = c('DPCA', 'black'))

  col <- rep('', along = method_names)
  label <- rep('', along = method_names)
  for (i in seq_along(legend_info)) {
    ind <- which(grepl(names(legend_info)[i], method_names))
    label[ind] <- legend_info[[i]][1]
    col[ind] <- legend_info[[i]][2]
  }

  return(list('col' = col, 'label' = label))
}

#' Plot density estimates of the run lengths
#'
#' WARNING: Before running this function
#' your must either have run \code{run_TEP_sim} yourself, or have downloaded
#' the results/ and thresholds/ directories from https://github.com/Tveten/tdpcaTEP
#' and put them in your working directory.
#'
#' @param lag The lag used in TDPCA.
#' @param fault_nr A number between 1 and 20 indicating the fault number.
#' @param keep_names A vector of names indicating which methods should be kept for
#' display (for internal use).
#' @param change_distr Which change distributions TDPCA to show results for
#' (see \code{?tpca::tpca}).
#' @param xlim A vector with lower and upper bounds for the x-axis.
#'
#' @export
rl_density <- function(lag, fault_nr, keep_names = NULL,
                       change_distr = c('semisparse_uniform',
                                        'semisparse_sd_only',
                                        'semisparse_mean_only'),
                       xlim = c(50, 300)) {
  n_distr <- length(change_distr)
  rl_df <- get_rl_table(lag, fault_nr, keep_names, change_distr)
  kappa <- 160 - lag

  lwd_lines <- 1.3
  legend_setup <- setup_legend(names(rl_df))
  col <- legend_setup$col
  label <- legend_setup$label
  adj <- 1
  dens_ests <- lapply(rl_df, function(rl_vec) density(rl_vec, adjust = adj))
  max_y <- max(unlist(lapply(dens_ests, function(dens_est) max(dens_est$y))))

  plot(dens_ests[[1]], xlab = 'Run length', col = col[1],
       xlim = xlim, ylim = c(0, max_y), lwd = lwd_lines,
       main = paste0('Fault ', fault_nr))
  for (i in 2:length(rl_df)) {
    lines(dens_ests[[i]], col = col[i], lwd = lwd_lines)
  }
  abline(v = kappa, lwd = 3, lty = 2)
  legend('topleft', unique(label), col = unique(col), lty = 1, lwd = lwd_lines)
}

p_false <- function(lag, fault_nr = 11, keep_names = NULL,
                    change_distr = c('semisparse_uniform',
                                     'semisparse_sd_only',
                                     'semisparse_mean_only')) {
  kappa <- 160 - lag
  rl_df <- get_rl_table(lag, fault_nr, keep_names, change_distr)
  apply(rl_df, 2, function(x) round(sum(x <= kappa) / length(x), 3))
}

est_edd <- function(lag, fault_nr = 11, keep_names = NULL,
                    change_distr = c('semisparse_uniform',
                                     'semisparse_sd_only',
                                     'semisparse_mean_only')) {
  kappa <- 160 - lag
  rl_df <- get_rl_table(lag, fault_nr, keep_names, change_distr)
  apply(rl_df, 2, function(x) round(mean(x[x > kappa]) - kappa, 1))
}

#' Make a table of EDD results and PFA for all tested methods.
#'
#' WARNING: Before running this function
#' your must either have run \code{run_TEP_sim} yourself, or have downloaded
#' the results/ and thresholds/ directories from https://github.com/Tveten/tdpcaTEP
#' and put them in your working directory.
#'
#' @param lag The lag used in TDPCA.
#' @param fault_nr A number between 1 and 20 indicating the fault number.
#' @param keep_names A vector of names indicating which methods should be kept for
#' display (for internal use).
#' @param change_distr Which change distributions TDPCA to show results for
#' (see \code{?tpca::tpca}).
#' @param sort_by If not NULL, either 'edd' or 'pfa', to sort the output table.
#' @param show_plot Logical. Should a density plot also be shown?
#'
#' @export
TEP_summary <- function(lag, fault_nr, keep_names = NULL,
                        change_distr = c('semisparse_uniform',
                                         'semisparse_sd_only',
                                         'semisparse_mean_only'),
                        sort_by = NULL, show_plot = TRUE, ...) {
  if (show_plot) rl_density(lag, fault_nr, keep_names, change_distr, ...)
  pfa <- p_false(lag, fault_nr, keep_names, change_distr)
  edd <- est_edd(lag, fault_nr, keep_names, change_distr)
  summary_table <- data.frame('method'  = names(edd),
                              'edd'     = edd,
                              'pfa' = pfa)
  row.names(summary_table) <- NULL
  if (!is.null(sort_by))
    summary_table <- summary_table[order(summary_table[[sort_by]]), ]
  summary_table
}

#' Make a table summarizing the EDD results for all faults.
#'
#' This function reproduces the TEP table in the paper.
#' WARNING: Before running this function
#' your must either have run \code{run_TEP_sim} yourself, or have downloaded
#' the results/ and thresholds/ directories from https://github.com/Tveten/tdpcaTEP
#' and put them in your working directory.
#'
#' Prints the PFAs and sorts the columns of the output table in the order of PFA.
#'
#' @export
pfa_edd_table <- function() {
  lag <- 5
  fault_nr <- 1:20
  keep_names <- c('max_pca20', 'min_pca20', 'tpca0.9 semisparse_uniform',
                  'tpca0.9 semisparse_mean_only', 'tpca0.99 semisparse_sd_only',
                  'dpca')

  labels <- setup_legend(keep_names)$label
  edd_table <- as.data.frame(matrix(0, ncol = length(labels),
                                       nrow = length(fault_nr)))
  names(edd_table) <- labels
  for (i in fault_nr) {
    edd_table[i, ] <- est_edd(lag, fault_nr[i], keep_names)
  }

  pfa <- p_false(lag, 1, keep_names)
  print('Probabilities of false alarms:')
  print(sort(pfa))
  print('EDD table:')
  edd_table[order(pfa)]
}

joint_TEP_summary <- function(lag,
                              change_distr = c('semisparse_uniform',
                                               'semisparse_sd_only',
                                               'semisparse_mean_only'),
                              sort_by = NULL) {
  table4 <- TEP_summary(lag, 4, show_plot = FALSE)
  table11 <- TEP_summary(lag, 11, show_plot = FALSE)
  joint_table <- cbind(table4, table11[, 2:3])
  names(joint_table) <- c('method', 'edd4', 'p_false4', 'edd11', 'p_false11')
  joint_table$method <- c('Max PCA (J = 1)', 'Max PCA (J = 5)', 'Max PCA (J = 20)',
                          'Min PCA (J = 1)', 'Min PCA (J = 5)', 'Min PCA (J = 20)',
                          'TPCA (uniform, c = 0.8)', 'TPCA (uniform, c = 0.99)', 'TPCA (uniform, c = 0.999)',
                          'TPCA (var, c = 0.8)', 'TPCA (var, c = 0.99)', 'TPCA (var, c = 0.999)',
                          'TPCA (mean, c = 0.8)', 'TPCA (mean, c = 0.99)', 'TPCA (mean, c = 0.999)')
  rownames(joint_table) <- NULL
  if (!is.null(sort_by))
    joint_table <- joint_table[order(joint_table[[sort_by]]), ]
  xtable::xtable(joint_table, digits = c(0, 0, 1, 3, 1, 3))
}

save_density_plots <- function(lag) {
  base_width <- 5
  base_height <- 5
  n_col <- 2
  n_row <- 1
  width <- n_col * base_width
  height <- n_row * base_height

  fault_nr <- c(10, 11)
  dir <- './results/figures/'
  file_name <- paste0('rl_densities_TEP_fault', paste(fault_nr, collapse = '-'), '.png')
  png(file = paste0(dir, file_name), width = width, height = height, units = 'in', res = 480)
  par(mfrow = c(1, 2))
  xlim <- list(c(60, 250), c(60, 250))
  keep_names <- c('max_pca20', 'min_pca20', 'tpca0.9 semisparse_uniform',
                  'tpca0.9 semisparse_mean_only', 'tpca0.99 semisparse_sd_only',
                  'dpca')
  rl_density(lag, fault_nr = fault_nr[1], keep_names, xlim = xlim[[1]])
  rl_density(lag, fault_nr = fault_nr[2], keep_names, xlim = xlim[[2]])
  dev.off()
}

get_TEP_thresholds <- function(change_distr) {
  get_tpca_train_info <- function(cov_mat_nr, change_distr) {
    dir <- './thresholds/axes/'
    axes_location <- paste0(dir, 'tpca_axes_', cov_mat_nr, '_m', m, 'd', d, '.txt')
    all_lines_vec <- readLines(axes_location)
    all_lines_split <- strsplit(all_lines_vec, ' ')
    axes_list <- list()
    cutoffs <- list()
    for (i in seq_along(all_lines_split)) {
      line <- all_lines_split[[i]]
      if (line[1] == change_distr) {
        axes_list[[length(axes_list) + 1]] <- as.numeric(line[3:length(line)])
        cutoffs[[length(cutoffs) + 1]] <- as.numeric(line[2])
      }
    }
    cutoffs <- unlist(cutoffs)
    names(axes_list) <- as.character(cutoffs)
    list('axes' = axes_list, 'cutoffs' = cutoffs)
  }

  get_tpca_threshold <- function(m, d, cov_mat_nr, axes) {
    dir <- './thresholds/'
    alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
    threshold_file <- paste0(dir, 'tpca_thresholds_', cov_mat_nr,
                             '_n', n, 'alpha', alpha_str,
                             'm', m, 'd', d, '_FINAL.txt')
    threshold_df <- read.table(threshold_file, header = TRUE, sep = ' ')
    threshold_df$axes <- as.character(threshold_df$axes)
    axes_in_df <- lapply(strsplit(threshold_df$axes, '-'), as.numeric)
    axes_in_df <- lapply(axes_in_df, sort)
    axes <- sort(axes)
    ind <- which(vapply(axes_in_df, is_equal_vectors, logical(1), y = axes))
    threshold <- threshold_df$threshold[ind]
    if (length(threshold) == 0)
      stop(paste0('Threshold for m = ', m, ', d = ', d, ', axes = (',
                  paste(axes, collapse = ', '), ') does not exist.'))
    if (length(threshold) > 1) {
      threshold <- threshold[1]
      warning('More than one threshold returned from file. Picking the first one.')
    }
    threshold
  }

  m <- 497
  d <- 132
  alpha <- 0.01
  n <- 160
  cov_mat_nr <- 2
  thresholds_list <- list()

  tpca_info <- get_tpca_train_info(cov_mat_nr, change_distr)
  tpca_threshold_df <- data.frame('cutoffs' = tpca_info$cutoffs)
  tpca_threshold_df$axes <- lapply(tpca_info$axes, function(axis) {
    paste(axis, collapse = '-')
  })
  tpca_threshold_df$threshold <- lapply(tpca_info$axes, function(axis) {
    get_tpca_threshold(m, d, cov_mat_nr, axis)
  })
  thresholds_list$tpca <- tpca_threshold_df

  r_pca <- c(1, 2, 3, 5, 10, 20)
  axes_max_pca <- lapply(r_pca, function(j) 1:j)
  axes_min_pca <- lapply(r_pca, function(j) (d - j + 1):d)
  min_pca_threshold_df <- data.frame('J' = r_pca)
  max_pca_threshold_df <- data.frame('J' = r_pca)
  min_pca_threshold_df$threshold <- lapply(axes_min_pca, function(axis) {
    get_tpca_threshold(m, d, cov_mat_nr, axis)
  })
  max_pca_threshold_df$threshold <- lapply(axes_max_pca, function(axis) {
    get_tpca_threshold(m, d, cov_mat_nr, axis)
  })

  thresholds_list$min_pca <- min_pca_threshold_df
  thresholds_list$max_pca <- max_pca_threshold_df
  thresholds_list
}
