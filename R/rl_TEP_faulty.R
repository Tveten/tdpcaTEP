rl_TEP_faulty <- function(test_sets, train_obj, lag, change_distr, edd_file,
                          w, r_pca, threshold_settings) {

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
    alpha_str <- strsplit(as.character(threshold_settings$alpha), '[.]')[[1]][2]
    threshold_file <- paste0(dir, 'tpca_thresholds_', cov_mat_nr,
                             '_n', threshold_settings$n, 'alpha', alpha_str,
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

  get_rl_df <- function(all_run_lengths, method, method_params) {
    rl_mat <- lapply(all_run_lengths, function(x) x[[method]])
    rl_mat <- do.call('rbind', rl_mat)
    rl_df <- as.data.frame(rl_mat)
    colnames(rl_df) <- vapply(method_params, function(x) paste0(method, x), character(1))
    rl_df
  }

  store_results <- function(all_run_lengths) {
    methods <- names(all_run_lengths[[1]])
    method_params <- list(r_pca, r_pca, c_tpca)
    rl_dfs <- lapply(seq_along(methods), function(i) {
      get_rl_df(all_run_lengths, methods[i], method_params[[i]])
    })
    entire_rl_df <- do.call('cbind', rl_dfs)
    write.table(entire_rl_df, edd_file)
  }

  x_train <- train_obj$x
  cov_mat_nr <- train_obj$nr
  x_train_e <- lag_extend(x_train, lag)
  d <- nrow(x_train_e)
  m <- ncol(x_train_e)
  mu_x <- rowMeans(x_train_e)
  sigma_x <- rowSds(x_train_e)
  x_train_e <- (x_train_e - mu_x) / sigma_x

  Sigma0_hat <- 1 / (m - 1) * x_train_e %*% t(x_train_e)
  pca_obj <- tpca::pca(Sigma0_hat)
  V <- pca_obj$vectors
  lambda <- pca_obj$values
  z_train <- V %*% x_train_e / sqrt(lambda)

  n_sim <- length(test_sets)

  # Method parameters
  axes_max_pca <- lapply(r_pca, function(j) 1:j)
  axes_min_pca <- lapply(r_pca, function(j) (d - j + 1):d)
  tpca_train_info <- get_tpca_train_info(cov_mat_nr, change_distr)
  c_tpca <- tpca_train_info$cutoffs
  axes_tpca <- tpca_train_info$axes

  sim_tpca_rl <- function(axes, z) {
    threshold <- get_tpca_threshold(m, d, cov_mat_nr, axes)
    z_train_sub <- z_train[axes, , drop = FALSE]
    z_sub <- z[axes, , drop = FALSE]
    mixture_rl(threshold, z_train_sub, z_sub, 1, w)
  }

  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  all_run_lengths <- foreach::foreach(l = 1:n_sim) %dopar% {
    x <- test_sets[[l]]
    x_e <- lag_extend(x, lag)
    x_e_norm <- (x_e - mu_x) / sigma_x
    z <- V %*% x_e_norm / sqrt(lambda)
    rl_max_pca <- vapply(axes_max_pca, sim_tpca_rl, numeric(1), z = z)
    rl_min_pca <- vapply(axes_min_pca, sim_tpca_rl, numeric(1), z = z)
    rl_tpca <- vapply(axes_tpca, sim_tpca_rl, numeric(1), z = z)
    rl_list <- list('max_pca' = rl_max_pca,
                    'min_pca' = rl_min_pca,
                    'tpca'    = rl_tpca)
    rl_list
  }
  stop_parallel(comp_cluster)

  store_results(all_run_lengths)
}
