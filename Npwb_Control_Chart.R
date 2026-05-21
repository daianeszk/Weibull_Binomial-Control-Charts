## ============================================
## Step 00 - Parameters
## ============================================

options(digits = 16)

# Sample sizes to test
n_vec <- c(5, 10, 15, 30, 50)

# p0 = false alarm probability
p0     <- 1 / 370
alpha0 <- p0

# Weibull in-control
k0      <- 5
lambda0 <- 4


## ============================================
## Step 01 - Helper functions
## ============================================

.validate_inputs <- function(Y, n_sample, alpha) {
  if (length(n_sample) != 1 || !is.finite(n_sample) || n_sample <= 0 || n_sample != as.integer(n_sample)) {
    stop("n_sample must be a positive integer.")
  }
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }
  if (any(!is.finite(Y)) || any(Y != as.integer(Y)) || any(Y < 1) || any(Y > n_sample)) {
    stop("Y must contain integers in the interval [1, n_sample].")
  }
  invisible(TRUE)
}

## Find p1 such that P(N >= Y | Bin(n, p1)) = alpha
## Exact inversion via Beta: p1 = qbeta(alpha; Y, n-Y+1)
find_p1 <- function(Y, n_sample, alpha) {
  .validate_inputs(Y, n_sample, alpha)
  qbeta(alpha, shape1 = Y, shape2 = n_sample - Y + 1)
}

## Detailed binomial table for a fixed p, rows indexed by N
make_binom_table <- function(p, n_sample) {
  if (length(p) != 1 || !is.finite(p) || p < 0 || p > 1) {
    stop("p must be in [0,1].")
  }
  if (n_sample != as.integer(n_sample) || n_sample <= 0) {
    stop("invalid n_sample.")
  }
  
  N <- 0:n_sample
  
  data.frame(
    N       = N,
    P_eq_N  = dbinom(N, size = n_sample, prob = p),
    P_leq_N = pbinom(N, size = n_sample, prob = p, lower.tail = TRUE),
    P_gt_N  = pbinom(N, size = n_sample, prob = p, lower.tail = FALSE),
    P_geq_N = pbinom(N - 1, size = n_sample, prob = p, lower.tail = FALSE)
  )
}

## Direction helper
get_mean_shift_direction <- function(mean_1, mean_0, tol = 1e-12) {
  direction <- rep("same", length(mean_1))
  direction[mean_1 < mean_0 - tol] <- "decrease"
  direction[mean_1 > mean_0 + tol] <- "increase"
  direction
}


## ============================================
## Step 02 - In-control Weibull moments
## ============================================

weibull_mean0 <- lambda0 * gamma(1 + 1 / k0)
weibull_var0  <- lambda0^2 * (gamma(1 + 2 / k0) - gamma(1 + 1 / k0)^2)

tol_mean <- 1e-12 * max(1, abs(weibull_mean0))


## ============================================
## Step 03 - Out-of-control scenarios
## ============================================

k1_vec <- c(3, 4, 5, 6, 7)

lambda1_vec <- c(
  2.47, 3.29, 3.91, 4.32, 4.94, 5.76,   # k1 = 3
  2.43, 3.24, 3.85, 4.25, 4.86, 5.67,   # k1 = 4
  2.40, 3.20, 3.80, 4.20, 4.80, 5.60,   # k1 = 5
  2.38, 3.17, 3.76, 4.16, 4.75, 5.54,   # k1 = 6
  2.36, 3.14, 3.73, 4.12, 4.71, 5.497   # k1 = 7
)

scenario_grid <- expand.grid(
  k1      = k1_vec,
  lambda1 = lambda1_vec
)


## ============================================
## Step 04 - Function to run one sample size
## ============================================

run_np_weibull_for_n <- function(n_sample) {
  
  cat("\n============================================\n")
  cat("Running n_sample =", n_sample, "\n")
  cat("============================================\n")
  
  Y_vec <- seq_len(n_sample)
  
  ## --------------------------------------------
  ## Design table
  ## --------------------------------------------
  
  p1_vec <- find_p1(Y_vec, n_sample, alpha0)
  
  alert_point_weibull_low_vec <- qweibull(
    p1_vec,
    shape = k0,
    scale = lambda0
  )
  
  alert_point_weibull_high_vec <- qweibull(
    1 - p1_vec,
    shape = k0,
    scale = lambda0
  )
  
  design_table <- data.frame(
    n_sample                 = n_sample,
    Y                        = Y_vec,
    p1                       = p1_vec,
    p0                       = alpha0,
    alert_point_weibull_low  = alert_point_weibull_low_vec,
    alert_point_weibull_high = alert_point_weibull_high_vec,
    mean_weibull_0           = weibull_mean0,
    var_weibull_0            = weibull_var0
  )
  
  design_table$check_p1_low_from_weibull0 <- pweibull(
    design_table$alert_point_weibull_low,
    shape = k0,
    scale = lambda0
  )
  
  design_table$check_p1_high_from_weibull0 <- pweibull(
    design_table$alert_point_weibull_high,
    shape = k0,
    scale = lambda0,
    lower.tail = FALSE
  )
  
  design_table$check_alpha_low_from_binom <- pbinom(
    design_table$Y - 1,
    size = n_sample,
    prob = design_table$check_p1_low_from_weibull0,
    lower.tail = FALSE
  )
  
  design_table$check_alpha_high_from_binom <- pbinom(
    design_table$Y - 1,
    size = n_sample,
    prob = design_table$check_p1_high_from_weibull0,
    lower.tail = FALSE
  )
  
  design_table$err_p1_low     <- design_table$check_p1_low_from_weibull0 - design_table$p1
  design_table$err_p1_high    <- design_table$check_p1_high_from_weibull0 - design_table$p1
  design_table$err_alpha_low  <- design_table$check_alpha_low_from_binom - design_table$p0
  design_table$err_alpha_high <- design_table$check_alpha_high_from_binom - design_table$p0
  
  
  ## --------------------------------------------
  ## Full binomial table
  ## --------------------------------------------
  
  binom_list <- lapply(seq_along(Y_vec), function(i) {
    
    y  <- Y_vec[i]
    p1 <- p1_vec[i]
    
    binom_table <- make_binom_table(p1, n_sample)
    
    P_signal <- pbinom(
      y - 1,
      size = n_sample,
      prob = p1,
      lower.tail = FALSE
    )
    
    data.frame(
      n_sample          = n_sample,
      Y                 = y,
      p1                = p1,
      P_signal_geq_Y    = P_signal,
      signal_if_N_geq_Y = binom_table$N >= y,
      binom_table
    )
  })
  
  full_results <- do.call(rbind, binom_list)
  
  
  ## --------------------------------------------
  ## Performance table: ARL1
  ## --------------------------------------------
  
  performance_list <- vector("list", nrow(scenario_grid))
  
  for (i in seq_len(nrow(scenario_grid))) {
    
    k1_i      <- scenario_grid$k1[i]
    lambda1_i <- scenario_grid$lambda1[i]
    
    mean_1_i <- lambda1_i * gamma(1 + 1 / k1_i)
    var_1_i  <- lambda1_i^2 * (gamma(1 + 2 / k1_i) - gamma(1 + 1 / k1_i)^2)
    
    mean_shift_direction <- get_mean_shift_direction(
      mean_1 = mean_1_i,
      mean_0 = weibull_mean0,
      tol    = tol_mean
    )
    
    if (mean_shift_direction == "increase") {
      
      alert_point_used <- design_table$alert_point_weibull_high
      
      p_ooc <- pweibull(
        alert_point_used,
        shape = k1_i,
        scale = lambda1_i,
        lower.tail = FALSE
      )
      
      chart_side_used <- "high"
      signal_rule <- "N_high >= Y, where N_high = # T >= t_high"
      
    } else {
      
      alert_point_used <- design_table$alert_point_weibull_low
      
      p_ooc <- pweibull(
        alert_point_used,
        shape = k1_i,
        scale = lambda1_i
      )
      
      chart_side_used <- "low"
      
      if (mean_shift_direction == "decrease") {
        signal_rule <- "N_low >= Y, where N_low = # T <= t_low"
      } else {
        signal_rule <- "same mean; low-side convention used"
      }
    }
    
    alpha1 <- pbinom(
      design_table$Y - 1,
      size = n_sample,
      prob = p_ooc,
      lower.tail = FALSE
    )
    
    ARL1 <- ifelse(alpha1 > 0, 1 / alpha1, Inf)
    
    tmp <- design_table
    
    tmp$k1                   <- k1_i
    tmp$lambda1              <- lambda1_i
    tmp$mean_weibull_1       <- mean_1_i
    tmp$var_weibull_1        <- var_1_i
    tmp$mean_shift_direction <- mean_shift_direction
    tmp$chart_side_used      <- chart_side_used
    tmp$signal_rule          <- signal_rule
    tmp$alert_point_used     <- alert_point_used
    tmp$p_ooc_weibull        <- p_ooc
    tmp$alpha1               <- alpha1
    tmp$ARL1                 <- ARL1
    
    performance_list[[i]] <- tmp
  }
  
  performance_table <- do.call(rbind, performance_list)
  
  
  ## --------------------------------------------
  ## Summary table
  ## --------------------------------------------
  
  summary_table <- performance_table[, c(
    "n_sample",
    "Y",
    "k1",
    "lambda1",
    "mean_weibull_1",
    "var_weibull_1",
    "mean_shift_direction",
    "chart_side_used",
    "ARL1",
    "p1",
    "p0",
    "alert_point_used",
    "p_ooc_weibull",
    "alpha1"
  )]
  
  
  ## --------------------------------------------
  ## Matrix table for Excel
  ## --------------------------------------------
  
  library(dplyr)
  library(tidyr)
  
  tcrit_low_vec  <- alert_point_weibull_low_vec
  tcrit_high_vec <- alert_point_weibull_high_vec
  
  arl_long <- scenario_grid %>%
    mutate(
      n_sample = n_sample,
      
      mean_weibull_1 = lambda1 * gamma(1 + 1 / k1),
      var_weibull_1  = lambda1^2 * (gamma(1 + 2 / k1) - gamma(1 + 1 / k1)^2),
      
      mean_shift_direction = get_mean_shift_direction(
        mean_1 = mean_weibull_1,
        mean_0 = weibull_mean0,
        tol    = tol_mean
      ),
      
      chart_side_used = ifelse(
        mean_shift_direction == "increase",
        "high",
        "low"
      ),
      
      signal_rule = ifelse(
        chart_side_used == "high",
        "N_high >= Y, where N_high = # T >= t_high",
        "N_low >= Y, where N_low = # T <= t_low"
      ),
      
      signal_rule = ifelse(
        mean_shift_direction == "same",
        "same mean; low-side convention used",
        signal_rule
      )
    ) %>%
    tidyr::crossing(Y = Y_vec) %>%
    mutate(
      tcrit_low  = tcrit_low_vec[Y],
      tcrit_high = tcrit_high_vec[Y],
      
      p_h1_low = pweibull(
        tcrit_low,
        shape = k1,
        scale = lambda1
      ),
      
      p_h1_high = pweibull(
        tcrit_high,
        shape = k1,
        scale = lambda1,
        lower.tail = FALSE
      ),
      
      tcrit = ifelse(
        chart_side_used == "high",
        tcrit_high,
        tcrit_low
      ),
      
      p_h1 = ifelse(
        chart_side_used == "high",
        p_h1_high,
        p_h1_low
      ),
      
      P_signal = pbinom(
        Y - 1,
        size = n_sample,
        prob = p_h1,
        lower.tail = FALSE
      ),
      
      ARL1 = ifelse(P_signal > 0, 1 / P_signal, NA_real_)
    ) %>%
    select(
      n_sample,
      mean_weibull_1,
      var_weibull_1,
      k1,
      lambda1,
      mean_shift_direction,
      chart_side_used,
      signal_rule,
      Y,
      ARL1
    )
  
  matrix_table <- arl_long %>%
    pivot_wider(
      names_from   = Y,
      values_from  = ARL1,
      names_prefix = "Y_",
      names_sort   = TRUE
    ) %>%
    arrange(k1, lambda1)
  
  
  return(list(
    n_sample          = n_sample,
    design_table      = design_table,
    full_results      = full_results,
    performance_table = performance_table,
    summary_table     = summary_table,
    arl_long          = arl_long,
    matrix_table      = matrix_table
  ))
}


## ============================================
## Step 05 - Run all sample sizes
## ============================================

results_list <- lapply(n_vec, run_np_weibull_for_n)
names(results_list) <- paste0("n_", n_vec)


## ============================================
## Step 06 - Combine and export CSV files
## ============================================

design_all <- do.call(rbind, lapply(results_list, function(x) x$design_table))
full_results_all <- do.call(rbind, lapply(results_list, function(x) x$full_results))
performance_all <- do.call(rbind, lapply(results_list, function(x) x$performance_table))
summary_all <- do.call(rbind, lapply(results_list, function(x) x$summary_table))
arl_long_all <- do.call(rbind, lapply(results_list, function(x) x$arl_long))

write.csv(
  design_all,
  "np_vs_weibull_all_n.csv",
  row.names = FALSE
)

write.csv(
  full_results_all,
  "np_probabilities_all_n.csv",
  row.names = FALSE
)

write.csv(
  performance_all,
  "np_weibull_performance_ARL1_all_n.csv",
  row.names = FALSE
)

write.csv(
  summary_all,
  "np_weibull_summary_ARL1_all_n.csv",
  row.names = FALSE
)

write.csv(
  arl_long_all,
  "np_weibull_ARL1_long_all_n.csv",
  row.names = FALSE
)


## ============================================
## Step 07 - Excel matrix report
## ============================================

library(openxlsx)

wb <- createWorkbook()

# Styles
style_4d <- createStyle(numFmt = "0.0000")
style_k1 <- createStyle(numFmt = "0")
style_l1 <- createStyle(numFmt = "0.00")

for (nm in names(results_list)) {
  
  matrix_table <- results_list[[nm]]$matrix_table
  sheet_name <- paste0("matrix_", nm)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, matrix_table)
  
  rows_data <- 2:(nrow(matrix_table) + 1)
  
  col_k1 <- which(names(matrix_table) == "k1")
  col_l1 <- which(names(matrix_table) == "lambda1")
  
  num_cols <- which(sapply(matrix_table, is.numeric))
  num_cols_4d <- setdiff(num_cols, c(col_k1, col_l1))
  
  addStyle(
    wb,
    sheet_name,
    style_4d,
    rows = rows_data,
    cols = num_cols_4d,
    gridExpand = TRUE
  )
  
  addStyle(
    wb,
    sheet_name,
    style_k1,
    rows = rows_data,
    cols = col_k1,
    gridExpand = TRUE
  )
  
  addStyle(
    wb,
    sheet_name,
    style_l1,
    rows = rows_data,
    cols = col_l1,
    gridExpand = TRUE
  )
  
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(matrix_table), widths = "auto")
}

# Add combined long table
addWorksheet(wb, "ARL1_long_all_n")
writeData(wb, "ARL1_long_all_n", arl_long_all)
freezePane(wb, "ARL1_long_all_n", firstRow = TRUE)
setColWidths(wb, "ARL1_long_all_n", cols = 1:ncol(arl_long_all), widths = "auto")

# Add combined summary table
addWorksheet(wb, "summary_all_n")
writeData(wb, "summary_all_n", summary_all)
freezePane(wb, "summary_all_n", firstRow = TRUE)
setColWidths(wb, "summary_all_n", cols = 1:ncol(summary_all), widths = "auto")

# Add design table
addWorksheet(wb, "design_all_n")
writeData(wb, "design_all_n", design_all)
freezePane(wb, "design_all_n", firstRow = TRUE)
setColWidths(wb, "design_all_n", cols = 1:ncol(design_all), widths = "auto")

saveWorkbook(
  wb,
  "np_weibull_matrix_all_n.xlsx",
  overwrite = TRUE
)


## ============================================
## Step 08 - Quick checks
## ============================================

cat("\n--- Design table, first rows ---\n")
print(head(design_all, 20))

cat("\n--- Summary table, first rows ---\n")
print(head(summary_all, 20))

cat("\nFiles exported:\n")
cat("- np_vs_weibull_all_n.csv\n")
cat("- np_probabilities_all_n.csv\n")
cat("- np_weibull_performance_ARL1_all_n.csv\n")
cat("- np_weibull_summary_ARL1_all_n.csv\n")
cat("- np_weibull_ARL1_long_all_n.csv\n")
cat("- np_weibull_matrix_all_n.xlsx\n")

# ============================================================
# Step 09 - Wide ARL comparison fixing each baseline n
# Baselines: n = 5, 15, 30, 50
# ============================================================

library(dplyr)
library(tidyr)
library(readr)

# ------------------------------------------------------------
# Check required object
# ------------------------------------------------------------

if (!exists("arl_long_all")) {
  stop("Object 'arl_long_all' was not found. Run Steps 00 to 08 first.")
}

# ------------------------------------------------------------
# Prepare base
# ------------------------------------------------------------

arl_compare_base <- arl_long_all |>
  mutate(
    direction = mean_shift_direction,
    n = n_sample,
    UCL = Y,
    mean_shift = (mean_weibull_1 - weibull_mean0) / weibull_mean0,
    var_shift  = (var_weibull_1 - weibull_var0) / weibull_var0,
    mean_shift_label = paste0(round(100 * mean_shift), "%"),
    var_shift_label  = paste0(round(100 * var_shift), "%"),
    inspected_units = n * ARL1
  ) |>
  filter(direction %in% c("decrease", "increase"))

# ------------------------------------------------------------
# Best ARL1 for each scenario and each n
# ------------------------------------------------------------

best_by_n <- arl_compare_base |>
  group_by(
    direction,
    n,
    k1,
    lambda1,
    mean_weibull_1,
    var_weibull_1,
    mean_shift_label,
    var_shift_label
  ) |>
  slice_min(ARL1, n = 1, with_ties = TRUE) |>
  summarise(
    best_ARL1 = min(ARL1, na.rm = TRUE),
    best_UCL = paste(UCL, collapse = ", "),
    best_inspected_units = min(inspected_units, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Wide base: one row per scenario, columns for each n
# ------------------------------------------------------------

best_by_n_wide <- best_by_n |>
  select(
    direction,
    k1,
    lambda1,
    mean_weibull_1,
    var_weibull_1,
    mean_shift_label,
    var_shift_label,
    n,
    best_ARL1,
    best_UCL,
    best_inspected_units
  ) |>
  pivot_wider(
    names_from = n,
    values_from = c(best_ARL1, best_UCL, best_inspected_units),
    names_prefix = "n_"
  ) |>
  arrange(direction, k1, lambda1)

# ------------------------------------------------------------
# Function to create gains fixing one baseline n
# ------------------------------------------------------------

make_gain_wide <- function(baseline_n, data_wide) {
  
  n_values <- c(5, 10, 15, 30, 50)
  comparison_n <- setdiff(n_values, baseline_n)
  
  baseline_arl_col <- paste0("best_ARL1_n_", baseline_n)
  baseline_units_col <- paste0("best_inspected_units_n_", baseline_n)
  
  out <- data_wide
  
  for (n_comp in comparison_n) {
    
    comp_arl_col <- paste0("best_ARL1_n_", n_comp)
    comp_units_col <- paste0("best_inspected_units_n_", n_comp)
    
    gain_arl_col <- paste0("gain_ARL_n", n_comp, "_vs_n", baseline_n)
    gain_units_col <- paste0("gain_units_n", n_comp, "_vs_n", baseline_n)
    
    out[[gain_arl_col]] <-
      100 * (out[[baseline_arl_col]] - out[[comp_arl_col]]) /
      out[[baseline_arl_col]]
    
    out[[gain_units_col]] <-
      100 * (out[[baseline_units_col]] - out[[comp_units_col]]) /
      out[[baseline_units_col]]
  }
  
  out |>
    mutate(baseline_n = baseline_n) |>
    relocate(baseline_n, .after = direction)
}

# ------------------------------------------------------------
# Create one wide comparison table for each fixed baseline
# ------------------------------------------------------------

gain_vs_n5  <- make_gain_wide(5,  best_by_n_wide)
gain_vs_n15 <- make_gain_wide(15, best_by_n_wide)
gain_vs_n30 <- make_gain_wide(30, best_by_n_wide)
gain_vs_n50 <- make_gain_wide(50, best_by_n_wide)

# Optional: combine all baselines in one file
gain_all_baselines_wide <- bind_rows(
  gain_vs_n5,
  gain_vs_n15,
  gain_vs_n30,
  gain_vs_n50
)

# ------------------------------------------------------------
# View results
# ------------------------------------------------------------

print(gain_vs_n5)
print(gain_vs_n15)
print(gain_vs_n30)
print(gain_vs_n50)

# ------------------------------------------------------------
# Export
# ------------------------------------------------------------

dir.create("outputs", showWarnings = FALSE)

write_csv(best_by_n, "outputs/best_by_n_long.csv")
write_csv(best_by_n_wide, "outputs/best_by_n_wide.csv")

write_csv(gain_vs_n5,  "outputs/gain_vs_n5_wide.csv")
write_csv(gain_vs_n15, "outputs/gain_vs_n15_wide.csv")
write_csv(gain_vs_n30, "outputs/gain_vs_n30_wide.csv")
write_csv(gain_vs_n50, "outputs/gain_vs_n50_wide.csv")

write_csv(gain_all_baselines_wide, "outputs/gain_all_baselines_wide.csv")

cat("\nFiles exported:\n")
cat("- outputs/best_by_n_long.csv\n")
cat("- outputs/best_by_n_wide.csv\n")
cat("- outputs/gain_vs_n5_wide.csv\n")
cat("- outputs/gain_vs_n15_wide.csv\n")
cat("- outputs/gain_vs_n30_wide.csv\n")
cat("- outputs/gain_vs_n50_wide.csv\n")
cat("- outputs/gain_all_baselines_wide.csv\n")
