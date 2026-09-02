# ======================================================================
# FILE 2 / 3: simulation
# - 01_setup.R 불러오기
# - Monte Carlo 반복 실행
# - replication-level raw 결과 저장
# ======================================================================

source("01_setup.R")

# ----------------------------------------------------------------------
# 7. Monte Carlo
# ----------------------------------------------------------------------

cat("=============================================================\n")
cat("Algorithm 1 simulation starts\n")
cat(sprintf("Source n = %d, Target n = %d\n", n_source, n_target))
cat(sprintf("Monte Carlo replications = %d\n", total_seeds))
cat("rho: Efficient Inference progressive estimation\n")
cat("=============================================================\n")


if (use_parallel) {

  num_cores <- min(
    max_cores,
    max(1, parallel::detectCores() - 1)
  )

  cat(sprintf("Parallel cores: %d\n", num_cores))

  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  results_combined <- foreach(
    s = seq_len(total_seeds),
    .combine = rbind,
    .packages = character(0)
  ) %dopar% {

    tryCatch(
      run_one_seed(s, return_full = FALSE),
      error = function(e) NULL
    )
  }

  parallel::stopCluster(cl)

} else {

  results_list <- lapply(
    seq_len(total_seeds),
    function(s) {
      tryCatch(
        run_one_seed(s, return_full = FALSE),
        error = function(e) NULL
      )
    }
  )

  results_list <- Filter(Negate(is.null), results_list)
  results_combined <- do.call(rbind, results_list)
}


if (is.null(results_combined) || nrow(results_combined) == 0) {
  stop("All simulation replications failed.")
}

cat(sprintf(
  "\nCompleted replications: %d / %d\n",
  nrow(results_combined),
  total_seeds
))

# ----------------------------------------------------------------------
# Save raw Monte Carlo results
# ----------------------------------------------------------------------

write.csv(
  results_combined,
  "EfficientRho_MC_raw.csv",
  row.names = FALSE
)

cat("\nSaved raw Monte Carlo results:\n")
cat("  -EfficientRho_MC_raw.csv\n")