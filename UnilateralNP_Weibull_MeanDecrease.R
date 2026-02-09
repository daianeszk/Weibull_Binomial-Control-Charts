## ============================================
## Step 00 - Parameters
## ============================================

options(digits = 16)

n_sample <- 10

# p0 = false alarm probability
p0     <- 1 / 370
alpha0 <- p0

# Weibull in-control
k0      <- 5
lambda0 <- 4

# Set of Y values (np-chart limits) to explore
Y_vec <- seq_len(n_sample)


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

## Detailed binomial table for a fixed p (rows indexed by N)
make_binom_table <- function(p, n_sample) {
  if (length(p) != 1 || !is.finite(p) || p < 0 || p > 1) stop("p must be in [0,1].")
  if (n_sample != as.integer(n_sample) || n_sample <= 0) stop("invalid n_sample.")
  
  N <- 0:n_sample
  
  data.frame(
    N       = N,
    P_eq_N  = dbinom(N, size = n_sample, prob = p),                        # P(X = N)
    P_leq_N = pbinom(N, size = n_sample, prob = p, lower.tail = TRUE),     # P(X <= N)
    P_gt_N  = pbinom(N, size = n_sample, prob = p, lower.tail = FALSE),    # P(X > N)
    P_geq_N = pbinom(N - 1, size = n_sample, prob = p, lower.tail = FALSE) # P(X >= N)
  )
}


## ============================================
## Step 02 - Convert p1 to Weibull "alert point"
## ============================================

weibull_mean0 <- lambda0 * gamma(1 + 1 / k0)
weibull_var0  <- lambda0^2 * (gamma(1 + 2 / k0) - gamma(1 + 1 / k0)^2)

# Critical p1 for each Y
p1_vec <- find_p1(Y_vec, n_sample, alpha0)

# Alert time (truncated life): t_y such that P(T <= t_y) = p1
alert_point_weibull_vec <- qweibull(p1_vec, shape = k0, scale = lambda0)

design_table <- data.frame(
  Y                   = Y_vec,
  p1                  = p1_vec,
  p0                  = alpha0,
  alert_point_weibull = alert_point_weibull_vec,
  mean_weibull_0      = weibull_mean0,
  var_weibull_0       = weibull_var0,
  n_sample            = n_sample
)

# Automatic checks (sanity checks)
design_table$check_p1_from_weibull0 <- pweibull(
  design_table$alert_point_weibull,
  shape = k0, scale = lambda0
)

design_table$check_alpha_from_binom <- pbinom(
  design_table$Y - 1,
  size = n_sample,
  prob = design_table$p1,
  lower.tail = FALSE
)

design_table$err_p1    <- design_table$check_p1_from_weibull0 - design_table$p1
design_table$err_alpha <- design_table$check_alpha_from_binom - design_table$p0


# Full binomial table (N=0..n_sample for each Y)
binom_list <- lapply(seq_along(Y_vec), function(i) {
  y  <- Y_vec[i]
  p1 <- p1_vec[i]
  
  binom_table <- make_binom_table(p1, n_sample)
  
  P_signal <- pbinom(y - 1, size = n_sample, prob = p1, lower.tail = FALSE)
  
  data.frame(
    Y                  = y,
    p1                 = p1,
    P_signal_geq_Y     = P_signal,
    signal_if_N_geq_Y  = (binom_table$N >= y),
    binom_table
  )
})

full_results <- do.call(rbind, binom_list)


## ============================================
## Step 03 - Export tables
## ============================================

cat("\n--- np_probabilities (full binomial table) ---\n")
print(head(full_results, 20))

write.csv(
  full_results,
  "np_probabilities.csv",
  row.names = FALSE
)

cat("\n--- np_vs_weibull (chart design) ---\n")
print(design_table)

write.csv(
  design_table,
  "np_vs_weibull.csv",
  row.names = FALSE
)


## ============================================
## Step 04 - Vary lambda and k (out-of-control) and compute ARL1
## ============================================

k1_vec      <- c(3, 4, 5, 6, 7)
lambda1_vec <- c(
  1.25, 1.5, 1.75, 2, 2.25, 2.36, 2.38, 2.4, 2.43, 2.47, 2.5, 2.75,
  3, 3.14, 3.17, 3.2, 3.24, 3.25, 3.29, 3.5, 3.73, 3.75, 3.76, 3.8,
  3.85, 3.91, 4, 4.12, 4.16, 4.2, 4.25, 4.32
)

scenario_grid <- expand.grid(
  k1      = k1_vec,
  lambda1 = lambda1_vec
)

performance_list <- vector("list", nrow(scenario_grid))

for (i in seq_len(nrow(scenario_grid))) {
  
  k1_i      <- scenario_grid$k1[i]
  lambda1_i <- scenario_grid$lambda1[i]
  
  mean_1_i <- lambda1_i * gamma(1 + 1 / k1_i)
  var_1_i  <- lambda1_i^2 * (gamma(1 + 2 / k1_i) - gamma(1 + 1 / k1_i)^2)
  
  # Out-of-control Weibull CDF evaluated at in-control alert points
  p_ooc <- pweibull(
    design_table$alert_point_weibull,
    shape = k1_i,
    scale = lambda1_i
  )
  
  # Signal probability under Bin(n, p_ooc)
  alpha1 <- pbinom(
    design_table$Y - 1,
    size = n_sample,
    prob = p_ooc,
    lower.tail = FALSE
  )
  
  ARL1 <- ifelse(alpha1 > 0, 1 / alpha1, Inf)
  
  tmp <- design_table
  tmp$k1             <- k1_i
  tmp$lambda1        <- lambda1_i
  tmp$mean_weibull_1 <- mean_1_i
  tmp$var_weibull_1  <- var_1_i
  tmp$p_ooc_weibull  <- p_ooc
  tmp$alpha1         <- alpha1
  tmp$ARL1           <- ARL1
  
  performance_list[[i]] <- tmp
}

performance_table <- do.call(rbind, performance_list)

cat("\n--- Chart performance for different (k1, lambda1) ---\n")
print(head(performance_table, 20))

write.csv(
  performance_table,
  "np_weibull_performance_ARL1.csv",
  row.names = FALSE
)


## ============================================
## Step 05 - Summary table
## ============================================

summary_table <- performance_table[, c(
  "n_sample",
  "Y",
  "k1",
  "lambda1",
  "ARL1",
  "p1",
  "p0",
  "alert_point_weibull"
)]

cat("\n--- FINAL TABLE (summary) ---\n")
print(head(summary_table, 20))

write.csv(
  summary_table,
  "np_weibull_summary_ARL1.csv",
  row.names = FALSE
)


## ============================================
## Step 06 - Excel matrix report 
## ============================================

library(dplyr)
library(tidyr)
library(openxlsx)

param_grid <- tidyr::crossing(k1 = k1_vec, lambda1 = lambda1_vec)

# Critical times under H0 (one per Y)
tcrit_vec <- alert_point_weibull_vec  # length = n_sample

# Long: (k1, lambda1, Y) -> ARL1
arl_long <- param_grid %>%
  mutate(
    mean_weibull_1 = lambda1 * gamma(1 + 1 / k1),
    var_weibull_1  = lambda1^2 * (gamma(1 + 2 / k1) - gamma(1 + 1 / k1)^2)
  ) %>%
  tidyr::crossing(Y = Y_vec) %>%
  mutate(
    tcrit    = tcrit_vec[Y],
    p_h1     = pweibull(tcrit, shape = k1, scale = lambda1),
    P_signal = pbinom(Y - 1, size = n_sample, prob = p_h1, lower.tail = FALSE),
    ARL1     = ifelse(P_signal > 0, 1 / P_signal, NA_real_)
  ) %>%
  select(mean_weibull_1, var_weibull_1, k1, lambda1, Y, ARL1)

# Wide (matrix) like the screenshot
matrix_table <- arl_long %>%
  pivot_wider(
    names_from  = Y,
    values_from = ARL1,
    names_sort  = TRUE
  ) %>%
  arrange(k1, lambda1)

# ----------------------------
# Export to Excel with formats
# ----------------------------
wb <- createWorkbook()
addWorksheet(wb, "matrix")

writeData(wb, "matrix", matrix_table)

# Styles (Excel number formats)
style_4d <- createStyle(numFmt = "0.0000")
style_k1 <- createStyle(numFmt = "0")     
style_l1 <- createStyle(numFmt = "0.00")

# Rows that contain data (exclude header row)
rows_data <- 2:(nrow(matrix_table) + 1)

# Column indices
col_k1 <- which(names(matrix_table) == "k1")
col_l1 <- which(names(matrix_table) == "lambda1")

# Apply 4 decimals to all numeric cols EXCEPT k1 and lambda1
num_cols <- which(sapply(matrix_table, is.numeric))
num_cols_4d <- setdiff(num_cols, c(col_k1, col_l1))

addStyle(wb, "matrix", style_4d, rows = rows_data, cols = num_cols_4d, gridExpand = TRUE)
addStyle(wb, "matrix", style_k1, rows = rows_data, cols = col_k1,     gridExpand = TRUE)
addStyle(wb, "matrix", style_l1, rows = rows_data, cols = col_l1,     gridExpand = TRUE)

saveWorkbook(wb, "np_weibull_matrix.xlsx", overwrite = TRUE)
