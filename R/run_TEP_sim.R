#' Load the TEP data to the global environment.
#'
#' This function loads the TEP data and puts it in the global environment.
#'
#' For this function to work, make sure you have downloaded the TEP data from
#' https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6C3JR1
#' and placed the fault-free training and faulty testing data in your working directory.
#'
#' The objects loaded are named "fault_free_training" and "faulty_testing".
#'
#' @param only_training Logical. If TRUE, only the training data is loaded.
#' The test data is a large file, so it can be beneficial to not load it if not needed.
#'
#' @export
load_TEP_data_globally <- function(only_training = FALSE) {
  # https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6C3JR1
  training_set <- "fault_free_training.RData"
  load(training_set, envir = .GlobalEnv)

  if (!only_training) {
    test_set <- "faulty_testing.RData"
    load(test_set, envir = .GlobalEnv)
  }
}

get_TEP_train <- function(sim_nr, return_all = TRUE) {
  # Depends on load_TEP_data_globally() having been run.
  training_sets <- lapply(sim_nr, function(i) {
    x_train <- subset(fault_free_training, simulationRun == i)
    x_train <- adjust_TEP_set(x_train)
    list('x'  = x_train, 'nr' = sim_nr)
  })
  if (return_all)
    return(training_sets)
  else
    return(lapply(training_sets, `[[`, 1))
}

get_TEP_test <- function(sim_nr, fault_nr, return_all = TRUE) {
  # Depends on load_TEP_data_globally() having been run.
  test_sets <- lapply(sim_nr, function(i) {
    x_test <- subset(faulty_testing, simulationRun == i & faultNumber == fault_nr)
    x_test <- adjust_TEP_set(x_test)
    list('x'  = x_test, 'nr' = sim_nr)
  })
  if (return_all)
    return(test_sets)
  else
    return(lapply(test_sets, `[[`, 1))
}

adjust_TEP_set <- function(x) {
  # Removes the dimensions that are not observations and transposes the data.
  rm_dim <- 1:3
  t(as.matrix(x[, -rm_dim]))
}

write_global_log <- function(text) {
  # text indicates where the simulations are at.
  curr_time <- Sys.time()
  write(paste0(text, ' Current time: ', curr_time, '. Time used so far: ',
               round(difftime(curr_time, start_time_G, units = 'hour'), 2), ' hours.'),
        log_file_G, append = TRUE)
}

make_dirs <- function() {
  # Creates the following directories in the current directory:
  # ./thresholds
  #     /axes
  #     /threshold_results
  # ./results
  #     /figures

  dir.create('thresholds', showWarnings = FALSE)
  dir.create('thresholds/axes', showWarnings = FALSE)
  dir.create('thresholds/threshold_results', showWarnings = FALSE)
  dir.create('results', showWarnings = FALSE)
  dir.create('results/figures', showWarnings = FALSE)
}

#' Reproduce the TEP simulation study.
#'
#' This function reproduces the results in
#' "Online Detection of Sparse Changes in High-Dimensional Data Streams Using Tailored Projections"
#' for the TEP data (Section 5). On an average computer it takes about two days to run on three cores.
#' Thus, the files with the results of running \code{run_TEP_sim} used in the paper
#' are included in the package or found at github.com/Tveten/tdpcaTEP.
#'
#' The TEP data are automatically loaded from the github repository to the
#' global environement.
#'
#' The following directories are created in the current directory:
#' ./thresholds,
#' ./thresholds/axes,
#' .thresholds/threshold_results,
#' ./results,
#' ./results/figures.
#' Files created to save information about thresholds are saved in the threshold
#' directory, while EDD results are stored in "results".
#' The results we obtained by running \code{run_TEP_sim} are stored in the
#' same structure of directories on Github.
#'
#' @export
run_TEP_sim <- function() {
  init_global_log <- function() {
    alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
    log_file_G <<- paste0('overall_log_', 'n', n, 'alpha', alpha_str,
                               'm', m, 'd', d, 'txt')
    start_time_G <<- Sys.time()
    write(paste0('Simulation start at ', start_time_G), log_file_G, append = TRUE)

    # Log every
    #   - find_mixture_thresholds
    #   - find_tpca_thresholds
    #   - In edd_sim:
    #      * run_mean_change_simulations
    #      * run_sd_change_simulations
    #      * run_cor_change_simulations
  }

  set_faulty_rl_file <- function(fault_nr, change_distr) {
    alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
    dir <- './results/'
    paste0(dir, 'TEP_rl_fault', fault_nr, '_', change_distr,
           '_n', n, 'alpha', alpha_str, 'm', m, 'd', d, 'l', lag, '.txt')
  }

  load_TEP_data_globally()
  lag <- 5
  TEP_train <- get_TEP_train(2)
  m <- ncol(TEP_train[[1]]$x) - lag
  d <- nrow(TEP_train[[1]]$x) * (lag + 1)

  # Method parameters
  r <- c(1, 2, 3, 5, 10, 20, 40)
  cutoff <- c(0.8, 0.9, 0.95, 0.99, 0.995, 0.999)
  change_distrs <- c('semisparse_mean_only',
                     'semisparse_sd_only',
                     'semisparse_uniform')

  # Threshold handling
  kappa <- 160 - lag
  n <- kappa
  alpha <- 0.01
  rel_tol <- c(0.2, 0.1, 0.05)
  init_global_log()
  make_dirs()
  find_all_thresholds(TEP_train, lag, n, alpha, rel_tol, r, cutoff, change_distrs)
  extract_final_thresholds(d, m, n, alpha)

  # sim_nr ranges from 1 to 500.
  # fault_nr from 1 to 20.
  # Faults inserted at sample nr. 160 in test sets.

  n_test_sets <- 500
  fault_nr <- 1:20
  w <- 200

  for (i in seq_along(fault_nr)) {
    TEP_test_sets <- get_TEP_test(1:n_test_sets, fault_nr[i], return_all = FALSE)
    for (j in seq_along(change_distrs)) {
      write_global_log(paste0('EDD sims for fault nr. ', fault_nr[i],
                              ' and change distribution ', change_distrs[j], '.'))
      rl_TEP_faulty(TEP_test_sets,
                    train_obj    = TEP_train[[1]],
                    lag          = lag,
                    change_distr = change_distrs[j],
                    edd_file     = set_faulty_rl_file(fault_nr[i], change_distrs[j]),
                    w            = w,
                    r_pca        = r,
                    threshold_settings = list('n' = n, 'alpha' = alpha))
    }
  }

  run_TEP_dpca(threshold = 'analytical')
  run_TEP_dpca(threshold = 'val')
}

test_tpca_TEP <- function(x_train) {
  m <- ncol(x_train)
  Sigma0_hat <- 1 / (m - 1) * x_train %*% t(x_train)
  D <- nrow(x_train)
  J <- 4
  A <- tpca::pca(Sigma0_hat, axes = (D - J + 1):D)$vectors
  y <- A %*% x_train

  par(mfrow = c(J, 1), mar = c(2.5, 2.5, 1.1, 2.5))
  invisible(lapply(1:J, function(i) plot(y[i, ], type = 'l')))
  Sys.sleep(10)
  par(mfrow = c(J, 1), mar = c(2.5, 2.5, 1.1, 2.5))
  invisible(lapply(1:J, function(i) acf(y[i, ])))
  Sys.sleep(10)
  par(mfrow = c(J, 1), mar = c(2.5, 2.5, 1.1, 2.5))
  invisible(lapply(1:J, function(i) {
    qqnorm(y[i, ])
    qqline(y[i, ])
  }))
}
