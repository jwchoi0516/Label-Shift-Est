# ======================================================================
# FILE 3 / 3: ARLSA results
# - 01_ARLSA_setup.R 불러오기
# - 02_ARLSA_simulation.R에서 저장한 raw 결과 읽기
# - Bias / SD / RMSE 계산
# - boxplot 및 rho diagnostic plot
# - summary 결과 저장
# ======================================================================

source("01_ARLSA_setup.R")

raw_result_file <- "ARLSA_EfficientRho_MC_raw.csv"

if (!file.exists(raw_result_file)) {
  stop(
    paste0(
      "'", raw_result_file, "' does not exist. ",
      "Run 02_ARLSA_simulation.R first."
    )
  )
}

results_combined <- read.csv(
  raw_result_file,
  stringsAsFactors = FALSE
)

cat(sprintf(
  "Loaded Monte Carlo replications: %d\n",
  nrow(results_combined)
))

# ----------------------------------------------------------------------
# 8. Numerical summary
# ----------------------------------------------------------------------

truth <- true_target_beta()

summarize_parameter <- function(param, truth_value) {

  cols <- paste0(
    param,
    c("_Naive", "_IPW", "_ARLSA", "_Oracle")
  )

  method_names <- c(
    "Naive",
    "IPW",
    "ARLSA",
    "Oracle"
  )

  out <- lapply(seq_along(cols), function(j) {

    x <- results_combined[[cols[j]]]
    x <- x[is.finite(x)]

    data.frame(
      Parameter = param,
      Method = method_names[j],
      N_success = length(x),
      Mean = mean(x),
      Bias = mean(x - truth_value),
      SD = sd(x),
      RMSE = sqrt(mean((x - truth_value)^2))
    )
  })

  do.call(rbind, out)
}


summary_table <- rbind(
  summarize_parameter("b0", truth["b0"]),
  summarize_parameter("b1", truth["b1"]),
  summarize_parameter("b2", truth["b2"]),
  summarize_parameter("b3", truth["b3"])
)

cat("\n=============================================================\n")
cat("Monte Carlo summary: Bias / SD / RMSE\n")
cat("=============================================================\n")

summary_print <- summary_table
numeric_cols <- vapply(summary_print, is.numeric, logical(1))
summary_print[, numeric_cols] <- lapply(
  summary_print[, numeric_cols, drop = FALSE],
  round,
  digits = 4
)
print(summary_print, row.names = FALSE)


cat("\n=============================================================\n")
cat("Density ratio diagnostics\n")
cat("=============================================================\n")

cat(sprintf(
  "Average rho RMSE          : %.4f\n",
  mean(results_combined$Rho_RMSE, na.rm = TRUE)
))

cat(sprintf(
  "Average raw rho range     : [%.4f, %.4f]\n",
  mean(results_combined$Rho_Min_Raw, na.rm = TRUE),
  mean(results_combined$Rho_Max_Raw, na.rm = TRUE)
))

cat(sprintf(
  "Average final rho range   : [%.4f, %.4f]\n",
  mean(results_combined$Rho_Min_Final, na.rm = TRUE),
  mean(results_combined$Rho_Max_Final, na.rm = TRUE)
))


# ----------------------------------------------------------------------
# 9. Boxplots: Naive / IPW / ARLSA / Oracle
# ----------------------------------------------------------------------

cat("\nDrawing boxplots...\n")

par(
  mfrow = c(2, 2),       # b0 b1 / b2 b3
  mar = c(5, 6, 4, 2)    # 왼쪽 method 이름 공간 조금 넓게
)

truth_vec <- c(
  b0 = truth["b0"],
  b1 = truth["b1"],
  b2 = truth["b2"],
  b3 = truth["b3"]
)

for (param in c("b0", "b1", "b2", "b3")) {

  cols <- paste0(
    param,
    c("_Naive", "_IPW", "_ARLSA", "_Oracle")
  )

  data_subset <- results_combined[, cols, drop = FALSE]

  # horizontal boxplot에서는 첫 번째 항목이 아래에 위치하므로
  # 위에서부터 Naive, IPW, ARLSA, Oracle이 되도록 역순으로 배치
  data_plot <- data_subset[, 4:1, drop = FALSE]

  boxplot(
    data_plot,
    names = c("Oracle", "ARLSA", "IPW", "Naive"),
    col = c("lightblue", "lightgreen", "orange", "salmon"),
    main = param,
    xlab = "Estimated Value",
    outline = TRUE,
    horizontal = TRUE,
    ylim=c(-1,1)
  )

  grid(nx = NULL, ny = NA)

  # horizontal이므로 true value는 세로선
  abline(
    v = truth_vec[param],
    lty = 2,
    lwd = 2
  )

  # 각 method의 Monte Carlo 평균
  means_subset <- colMeans(
    data_plot,
    na.rm = TRUE
  )

  # horizontal이므로 x = 추정값, y = method 위치
  points(
    x = means_subset,
    y = 1:4,
    pch = 23,
    bg = "white",
    cex = 1.3
  )
}

par(mfrow = c(1, 1))


# ----------------------------------------------------------------------
# 10. rho diagnostic plot for seed 1
# ----------------------------------------------------------------------

diag_seed1 <- run_one_seed(
  1,
  return_full = TRUE
)

rho_diag <- diag_seed1$rho

plot(
  rho_diag$y_grid,
  true_rho_func(rho_diag$y_grid),
  type = "l",
  lwd = 2,
  xlab = "y",
  ylab = expression(rho(y)),
  main = "Density ratio diagnostic (Seed 1)"
)

lines(
  rho_diag$y_grid,
  rho_diag$rho_star_grid,
  lty = 2
)

lines(
  rho_diag$y_grid,
  rho_diag$rho_tilde_grid,
  lty = 3,
  lwd = 2
)

lines(
  rho_diag$y_grid,
  rho_diag$rho_hat_grid,
  lty = 4,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "True rho",
    "Working rho*",
    "rho_tilde",
    "rho_hat"
  ),
  lty = c(1, 2, 3, 4),
  lwd = c(2, 1, 2, 2),
  bty = "n"
)

# ----------------------------------------------------------------------
# 11. Save summary results
# ----------------------------------------------------------------------

write.csv(
  summary_table,
  "ARLSA_EfficientRho_MC_summary.csv",
  row.names = FALSE
)

cat("\n=============================================================\n")
cat("Result analysis complete.\n")
cat("Saved file:\n")
cat("  - ARLSA_EfficientRho_MC_summary.csv\n")
cat("=============================================================\n")


# ----------------------------------------------------------------------
# X distribution diagnostic: Source vs Target
# ----------------------------------------------------------------------

x_diag <- generate_data(
  seed = 1,
  n1 = n_source,
  n0 = n_target
)

par(
  mfrow = c(1, 3),
  mar = c(5, 4, 4, 1)
)

for (j in 1:3) {
  
  x_s_j <- x_diag$x_s[, j]
  x_t_j <- x_diag$x_t[, j]
  
  # Source와 Target에서 동일한 구간 사용
  common_breaks <- pretty(
    range(c(x_s_j, x_t_j)),
    n = 20
  )
  
  hist(
    x_s_j,
    breaks = common_breaks,
    probability = TRUE,
    col = rgb(1, 0, 0, 0.35),
    border = "white",
    main = paste0("Distribution of X", j),
    xlab = paste0("X", j),
    ylim = c(
      0,
      max(
        hist(x_s_j, breaks = common_breaks, plot = FALSE)$density,
        hist(x_t_j, breaks = common_breaks, plot = FALSE)$density
      ) * 1.15
    )
  )
  
  hist(
    x_t_j,
    breaks = common_breaks,
    probability = TRUE,
    col = rgb(0, 0, 1, 0.35),
    border = "white",
    add = TRUE
  )
  
  lines(
    density(x_s_j),
    lwd = 2,
    lty = 1
  )
  
  lines(
    density(x_t_j),
    lwd = 2,
    lty = 2
  )
  
  legend(
    "topright",
    legend = c("Source", "Target"),
    fill = c(
      rgb(1, 0, 0, 0.35),
      rgb(0, 0, 1, 0.35)
    ),
    bty = "n"
  )
}

par(mfrow = c(1, 1))