find_tdpca_thresholds <- function(train_obj, lag, n, alpha, rel_tol, r,
                                 cutoff, change_distrs) {
  unique_vectors <- function(list_of_vectors) {
    list_of_vectors[!duplicated(lapply(list_of_vectors, sort))]
  }

  log_tpca <- function(change_distr, cutoff, axes) {
    out <- c(change_distr, cutoff, axes)
    dir <- './thresholds/axes/'
    tpca_log_file <- paste0(dir, 'tpca_axes_', cov_mat_nr,
                            '_m', m, 'd', d, '.txt')
    write(out, file = tpca_log_file, ncolumns = length(out), append = TRUE)
  }

  which_axes <- function(prop_max, keep_prop, max_axes) {
    order_axes <- order(prop_max, decreasing = TRUE)
    cum_prop <- cumsum(prop_max[order_axes])
    n_keep <- min(sum(cum_prop < keep_prop) + 1, max_axes)
    order_axes[1:n_keep]
  }

  get_axes_list <- function(x_train, max_axes) {
    x_train_e <- lag_extend(x_train, lag)
    x_train_e <- (x_train_e - rowMeans(x_train_e)) / rowSds(x_train_e)
    m2 <- ncol(x_train_e)
    cov_mat <- 1 / (m2 - 1) * x_train_e %*% t(x_train_e)

    axes_pca_max <- lapply(r, function(j) 1:j)
    axes_pca_min <- lapply(r, function(j) (d - j + 1):d)
    axes_list <- list(axes_pca_max, axes_pca_min)
    for (i in seq_along(change_distrs)) {
      prop_max <- tpca::tdpca(cov_mat, lag = lag, n_sim = 10^4,
                              change_distr = change_distrs[i],
                              cutoff = cutoff[1],
                              max_axes = max_axes)$prop_axes_max
      which_axes_list <- list()
      for (j in seq_along(cutoff)) {
        which_axes_list[[j]] <- which_axes(prop_max, cutoff[j], max_axes)
        log_tpca(change_distrs[i], cutoff[j], which_axes_list[[j]])
      }
      axes_list[[2 + i]] <- which_axes_list
    }
    axes_list <- unlist(axes_list, recursive = FALSE)
    unique_vectors(axes_list)
  }

  get_axes_list_from_file <- function(cov_mat_nr) {
    axes_pca_max <- lapply(r, function(j) 1:j)
    axes_pca_min <- lapply(r, function(j) (d - j + 1):d)
    axes_list <- c(axes_pca_max, axes_pca_min)
    dir <- './thresholds/axes/'
    axes_location <- paste0(dir, 'tpca_axes_', cov_mat_nr, '_m', m, 'd', d, '.txt')
    all_lines_vec <- readLines(axes_location)
    all_lines_split <- strsplit(all_lines_vec, ' ')
    for (i in seq_along(all_lines_split)) {
      line <- all_lines_split[[i]]
      axes_list[[length(axes_list) + 1]] <- as.numeric(line[3:length(line)])
    }
    unique_vectors(axes_list)
  }

  write_global_log(paste0('Finding tpca thresholds for TEP training set nr. ', train_obj$nr, '.'))

  x_train <- train_obj$x
  cov_mat_nr <- train_obj$nr
  m <- ncol(x_train) - lag
  d <- nrow(x_train) * (lag + 1)
  max_axes <- min(round(d / 5), 40)
  axes_list <- get_axes_list(x_train, max_axes)
  # axes_list <- get_axes_list_from_file(2)

  for (i in seq_along(axes_list)) {
    write_global_log(paste0('Finding tpca thresholds for axis ', i,
                            ' of ', length(axes_list), '.'))
    threshold_finder(x_train, lag, n, alpha, axes = axes_list[[i]],
                     rel_tol = rel_tol,
                     file_id = paste0(cov_mat_nr, '_'))
  }
}

find_all_thresholds <- function(training_sets, lag,  n, alpha, rel_tol, r,
                                cutoff, change_distrs) {
  lapply(training_sets, find_tdpca_thresholds, lag = lag, n = n, alpha = alpha,
         rel_tol = rel_tol, r = r, cutoff = cutoff, change_distrs = change_distrs)
}

extract_final_thresholds <- function(d, m, n, alpha) {
  get_final_threshold <- function(curr_file, id) {
    threshold_df <- read.table(paste0(path2, curr_file), header = TRUE, sep = ' ')
    split_id <- strsplit(id, 'n|alpha|m|d')
    n <- as.numeric(split_id[[1]][2])
    alpha <- as.numeric(paste0('0.', split_id[[1]][3]))
    target_arl <- n / alpha
    threshold_sub_df <- subset(threshold_df, n_sim == max(threshold_df$n_sim))
    final_threshold_ind <- which.min(abs(threshold_sub_df$arl_est - target_arl))
    threshold_sub_df$threshold[final_threshold_ind]
  }

  alpha <- strsplit(as.character(alpha), '[.]')[[1]][2]
  path1 <- './thresholds/'
  path2 <- paste0(path1, 'threshold_results/')

  files <- list_files_matching(path2, 'tpca', 'results', paste0('d', d),
                               paste0('m', m), paste0('n', n), paste0('alpha', alpha))
  split_files <- strsplit(files, '_|p0|ax|[.]')
  id_ind <- c(3, 5)
  simulation_id <- vapply(split_files, function(x) {
      paste(x[id_ind], collapse = '_')
    }, character(1))
  unique_id <- unique(simulation_id)

  final_files <- paste0(path1, 'tpca_thresholds_', unique_id, '_FINAL.txt')
  for (i in seq_along(unique_id)) {
    which_files <- simulation_id == unique_id[i]
    curr_files <- files[which_files]
    curr_split_files <- split_files[which_files]

    axes_id <- vapply(curr_split_files, function(split_file) split_file[6], character(1))
    write(c('axes', 'threshold'), final_files[i], ncolumns = 2)
    for (j in seq_along(curr_files)) {
      final_threshold <- get_final_threshold(curr_files[j], unique_id[i])
      write(c(axes_id[j], final_threshold), final_files[i], ncolumns = 2, append = TRUE)
    }
  }
}
