## ============================================================
## Sweep study:
## Double-sampling np-Weibull control chart
## Vary n_a, n_b, alpha_a, Y_a, Y_b
## Target: global ARL0 close to 370
##
## Revised version:
## - Double sampling only
## - n_b <= 15
## - ASS restriction: 5 <= ASS <= 10
## - Exports simple LaTeX tables for Stage A and Stage B
## ============================================================

options(digits = 16)


## ============================================================
## Step 00 - Global parameters
## ============================================================

ARL0_target <- 370
alpha_global_h0 <- 1 / ARL0_target

ASS_min <- 5
ASS_max <- 10

n_a_vec <- 1:ASS_max

n_b_max <- 15
n_b_vec <- 1:n_b_max

alpha_a_vec <- c(0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50)

k0      <- 5
lambda0 <- 4

weibull_mean0 <- lambda0 * gamma(1 + 1 / k0)

weibull_var0 <- lambda0^2 * (
  gamma(1 + 2 / k0) - gamma(1 + 1 / k0)^2
)

tol_mean <- 1e-12 * max(1, abs(weibull_mean0))

## Selected design for individual LaTeX tables
selected_n_a <- 6
selected_alpha_a <- 0.05

selected_n_b <- 6
selected_alpha_b <- alpha_global_h0 / selected_alpha_a

cat("\n============================================================\n")
cat("GLOBAL SETTINGS\n")
cat("============================================================\n")
cat("ARL0_target     =", ARL0_target, "\n")
cat("alpha_global_h0 =", alpha_global_h0, "\n")
cat("ASS_min         =", ASS_min, "\n")
cat("ASS_max         =", ASS_max, "\n")
cat("n_b_max         =", n_b_max, "\n")
cat("k0              =", k0, "\n")
cat("lambda0         =", lambda0, "\n")
cat("Weibull mean H0 =", weibull_mean0, "\n")
cat("Weibull var H0  =", weibull_var0, "\n")


## ============================================================
## Step 01 - Helper functions
## ============================================================

.validate_inputs <- function(Y, n_sample, alpha) {
  
  if (
    length(n_sample) != 1 ||
    !is.finite(n_sample) ||
    n_sample <= 0 ||
    n_sample != as.integer(n_sample)
  ) {
    stop("n_sample must be a positive integer.")
  }
  
  if (
    length(alpha) != 1 ||
    !is.finite(alpha) ||
    alpha <= 0 ||
    alpha >= 1
  ) {
    stop("alpha must be in (0,1).")
  }
  
  if (
    any(!is.finite(Y)) ||
    any(Y != as.integer(Y)) ||
    any(Y < 1) ||
    any(Y > n_sample)
  ) {
    stop("Y must contain integers in the interval [1, n_sample].")
  }
  
  invisible(TRUE)
}


find_p1 <- function(Y, n_sample, alpha) {
  
  .validate_inputs(
    Y        = Y,
    n_sample = n_sample,
    alpha    = alpha
  )
  
  qbeta(
    alpha,
    shape1 = Y,
    shape2 = n_sample - Y + 1
  )
}


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


get_mean_shift_direction <- function(mean_1, mean_0, tol = 1e-12) {
  
  direction <- rep("same", length(mean_1))
  direction[mean_1 < mean_0 - tol] <- "decrease"
  direction[mean_1 > mean_0 + tol] <- "increase"
  
  direction
}


make_alpha_label <- function(x) {
  
  out <- formatC(x, format = "fg", digits = 8)
  out <- gsub("\\.", "p", out)
  out <- gsub("-", "m", out)
  
  out
}


safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}


safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}


safe_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}


safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}


make_stage_design <- function(stage, n_sample, alpha, k0, lambda0) {
  
  Y_vec <- seq_len(n_sample)
  
  p1_vec <- find_p1(
    Y        = Y_vec,
    n_sample = n_sample,
    alpha    = alpha
  )
  
  tcrit_low_vec <- qweibull(
    p1_vec,
    shape = k0,
    scale = lambda0
  )
  
  tcrit_high_vec <- qweibull(
    1 - p1_vec,
    shape = k0,
    scale = lambda0
  )
  
  out <- data.frame(
    stage       = stage,
    n_sample    = n_sample,
    alpha_stage = alpha,
    Y           = Y_vec,
    p1          = p1_vec,
    tcrit_low   = tcrit_low_vec,
    tcrit_high  = tcrit_high_vec,
    stringsAsFactors = FALSE
  )
  
  out$check_p1_low_from_weibull0 <- pweibull(
    out$tcrit_low,
    shape = k0,
    scale = lambda0
  )
  
  out$check_p1_high_from_weibull0 <- pweibull(
    out$tcrit_high,
    shape = k0,
    scale = lambda0,
    lower.tail = FALSE
  )
  
  out$check_alpha_low_from_binom <- pbinom(
    out$Y - 1,
    size = out$n_sample,
    prob = out$check_p1_low_from_weibull0,
    lower.tail = FALSE
  )
  
  out$check_alpha_high_from_binom <- pbinom(
    out$Y - 1,
    size = out$n_sample,
    prob = out$check_p1_high_from_weibull0,
    lower.tail = FALSE
  )
  
  out$err_alpha_low <- out$check_alpha_low_from_binom - out$alpha_stage
  out$err_alpha_high <- out$check_alpha_high_from_binom - out$alpha_stage
  
  out
}


make_design_pairs <- function(design_a, design_b) {
  
  design_a_perf <- data.frame(
    Y_a          = design_a$Y,
    p1_a         = design_a$p1,
    tcrit_low_a  = design_a$tcrit_low,
    tcrit_high_a = design_a$tcrit_high,
    stringsAsFactors = FALSE
  )
  
  design_b_perf <- data.frame(
    Y_b          = design_b$Y,
    p1_b         = design_b$p1,
    tcrit_low_b  = design_b$tcrit_low,
    tcrit_high_b = design_b$tcrit_high,
    stringsAsFactors = FALSE
  )
  
  idx <- expand.grid(
    row_a = seq_len(nrow(design_a_perf)),
    row_b = seq_len(nrow(design_b_perf))
  )
  
  out <- cbind(
    design_a_perf[idx$row_a, , drop = FALSE],
    design_b_perf[idx$row_b, , drop = FALSE]
  )
  
  rownames(out) <- NULL
  
  out
}


summarise_by_rule <- function(df, rule_cols) {
  
  if (nrow(df) == 0) {
    warning("Input data frame has zero rows.")
    return(data.frame())
  }
  
  split_list <- split(df, df$rule_id)
  
  out <- do.call(
    rbind,
    lapply(names(split_list), function(id) {
      
      d <- split_list[[id]]
      base <- d[1, rule_cols, drop = FALSE]
      
      data.frame(
        base,
        scenario_count = nrow(d),
        
        ARL1_mean      = safe_mean(d$ARL1),
        ARL1_median    = safe_median(d$ARL1),
        ARL1_max       = safe_max(d$ARL1),
        ARL1_min       = safe_min(d$ARL1),
        
        Avg_Sample1_mean   = safe_mean(d$Avg_Sample1_per_cycle),
        Avg_Sample1_median = safe_median(d$Avg_Sample1_per_cycle),
        Avg_Sample1_max    = safe_max(d$Avg_Sample1_per_cycle),
        Avg_Sample1_min    = safe_min(d$Avg_Sample1_per_cycle),
        
        P_stop_after_A_h0_mean = safe_mean(d$P_stop_after_A_h0),
        P_stop_after_A_h1_mean = safe_mean(d$P_stop_after_A_h1),
        
        tcritsum_mean      = safe_mean(d$tcritsum),
        tcritsum_median    = safe_median(d$tcritsum),
        tcritsum_max       = safe_max(d$tcritsum),
        tcritsum_min       = safe_min(d$tcritsum),
        
        Expected_Alert_Time_H0_mean =
          safe_mean(d$Expected_Alert_Time_H0_per_cycle),
        
        Expected_Alert_Time_H1_mean =
          safe_mean(d$Expected_Alert_Time_H1_per_cycle),
        
        Expected_Alert_Time_H1_median =
          safe_median(d$Expected_Alert_Time_H1_per_cycle),
        
        Expected_Alert_Time_H1_max =
          safe_max(d$Expected_Alert_Time_H1_per_cycle),
        
        Expected_Alert_Time_H1_min =
          safe_min(d$Expected_Alert_Time_H1_per_cycle),
        
        Expected_Alert_Time_To_Signal_H1_mean =
          safe_mean(d$Expected_Alert_Time_To_Signal_H1),
        
        Expected_Alert_Time_To_Signal_H1_median =
          safe_median(d$Expected_Alert_Time_To_Signal_H1),
        
        Expected_Alert_Time_To_Signal_H1_max =
          safe_max(d$Expected_Alert_Time_To_Signal_H1),
        
        Expected_Alert_Time_To_Signal_H1_min =
          safe_min(d$Expected_Alert_Time_To_Signal_H1),
        
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  
  out
}


get_best_per_scenario <- function(df) {
  
  if (nrow(df) == 0) {
    warning("Input data frame has zero rows.")
    return(data.frame())
  }
  
  split_list <- split(df, df$scenario_id)
  
  out <- do.call(
    rbind,
    lapply(names(split_list), function(id) {
      
      d <- split_list[[id]]
      
      d <- d[order(
        d$ARL1,
        d$Expected_Alert_Time_H1_per_cycle,
        d$Expected_Alert_Time_To_Signal_H1,
        d$tcritsum,
        d$Avg_Sample1_per_cycle
      ), ]
      
      d[1, , drop = FALSE]
    })
  )
  
  rownames(out) <- NULL
  
  out
}


make_simple_stage_latex_table <- function(
    stage_df,
    stage_letter,
    side = c("decrease", "increase"),
    caption = NULL,
    label = NULL
) {
  
  side <- match.arg(side)
  
  stage_df <- stage_df[order(stage_df$UCL), ]
  
  n_stage <- unique(stage_df$n_stage)
  alpha_stage <- unique(stage_df$alpha_stage)
  
  if (length(n_stage) != 1 || length(alpha_stage) != 1) {
    stop("stage_df must contain exactly one n_stage and one alpha_stage.")
  }
  
  if (side == "decrease") {
    ucl_symbol <- paste0("$UCL_{y", stage_letter, "}$")
    p_symbol   <- paste0("$p_{y", stage_letter, "}$")
    t_symbol   <- paste0("$t_{y", stage_letter, "}$")
    t_values   <- stage_df$t_y
    side_text  <- "mean decrease"
  } else {
    ucl_symbol <- paste0("$UCL_{w", stage_letter, "}$")
    p_symbol   <- paste0("$p_{w", stage_letter, "}$")
    t_symbol   <- paste0("$t_{w", stage_letter, "}$")
    t_values   <- stage_df$t_w
    side_text  <- "mean increase"
  }
  
  n_symbol <- ifelse(stage_letter == "a", "n_a", "n_b")
  alpha_symbol <- ifelse(stage_letter == "a", "\\alpha_a", "\\alpha_b")
  stage_name <- ifelse(stage_letter == "a", "first stage", "second stage")
  
  if (is.null(caption)) {
    caption <- paste0(
      "In-control design parameters of the proposed $np^{2s}_{wb}$ chart for the ",
      stage_name,
      ". Out-of-control case: ",
      side_text,
      "."
    )
  }
  
  if (is.null(label)) {
    label <- paste0(
      "tab:np2s_",
      stage_letter,
      "_",
      side,
      "_n",
      n_stage,
      "_alpha",
      make_alpha_label(alpha_stage)
    )
  }
  
  header_values <- paste0("\\textbf{", stage_df$UCL, "}", collapse = " & ")
  p_values <- paste0(sprintf("%.5f", stage_df$p_y), collapse = " & ")
  t_values_fmt <- paste0(sprintf("%.4f", t_values), collapse = " & ")
  
  n_cols <- nrow(stage_df) + 1
  col_spec <- paste0("c|", paste(rep("c", nrow(stage_df)), collapse = ""))
  
  latex <- paste0(
    "\\begin{table}[H]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n\n",
    "\\scriptsize\n",
    "\\setlength{\\tabcolsep}{3pt}\n",
    "\\renewcommand{\\arraystretch}{1.15}\n\n",
    "\\setlength{\\aboverulesep}{0pt}\n",
    "\\setlength{\\belowrulesep}{0pt}\n\n",
    "\\begin{adjustbox}{max width=\\textwidth}\n",
    "\\begin{tabular}{", col_spec, "}\n",
    "\\toprule\n",
    "\\multicolumn{", n_cols, "}{c}{$",
    n_symbol, "=", n_stage, "$, $",
    alpha_symbol, "=", sprintf("%.6f", alpha_stage), "$} \\\\\n",
    "\\midrule\n",
    ucl_symbol, " & ", header_values, " \\\\\n",
    "\\midrule\n",
    p_symbol, " & ", p_values, " \\\\\n",
    t_symbol, " & ", t_values_fmt, " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{adjustbox}\n\n",
    "\\end{table}\n\n"
  )
  
  latex
}


write_all_simple_stage_tables <- function(stage_summary, stage_letter, side, file_name) {
  
  groups <- unique(stage_summary[, c("n_stage", "alpha_stage")])
  groups <- groups[order(groups$n_stage, groups$alpha_stage), ]
  
  tables <- vector("list", nrow(groups))
  
  for (i in seq_len(nrow(groups))) {
    
    n_i <- groups$n_stage[i]
    alpha_i <- groups$alpha_stage[i]
    
    d_i <- stage_summary[
      stage_summary$n_stage == n_i &
        abs(stage_summary$alpha_stage - alpha_i) < 1e-14,
    ]
    
    tables[[i]] <- make_simple_stage_latex_table(
      stage_df = d_i,
      stage_letter = stage_letter,
      side = side
    )
  }
  
  writeLines(unlist(tables), con = file_name)
}


write_selected_simple_stage_table <- function(
    stage_summary,
    stage_letter,
    side,
    n_selected,
    alpha_selected,
    file_name
) {
  
  d <- stage_summary[
    stage_summary$n_stage == n_selected &
      abs(stage_summary$alpha_stage - alpha_selected) < 1e-14,
  ]
  
  if (nrow(d) == 0) {
    warning("Selected stage table not found: ", file_name)
    return(invisible(NULL))
  }
  
  latex <- make_simple_stage_latex_table(
    stage_df = d,
    stage_letter = stage_letter,
    side = side
  )
  
  writeLines(latex, con = file_name)
}


## ============================================================
## Step 02 - Out-of-control scenarios
## ============================================================

k1_vec <- c(3, 4, 5, 6, 7)

lambda1_values <- c(
  2.47, 3.29, 3.91, 4.32, 4.94, 5.76,
  2.43, 3.24, 3.85, 4.25, 4.86, 5.67,
  2.40, 3.20, 3.80, 4.20, 4.80, 5.60,
  2.38, 3.17, 3.76, 4.16, 4.75, 5.54,
  2.36, 3.14, 3.73, 4.12, 4.71, 5.497
)

n_lambda_per_k <- length(lambda1_values) / length(k1_vec)

if (n_lambda_per_k != as.integer(n_lambda_per_k)) {
  stop("lambda1_values must have the same number of values for each k1.")
}

scenario_grid <- data.frame(
  k1      = rep(k1_vec, each = n_lambda_per_k),
  lambda1 = lambda1_values
)

scenario_grid$scenario_id <- paste0(
  "k", scenario_grid$k1,
  "_lambda", gsub("\\.", "p", as.character(scenario_grid$lambda1))
)


## ============================================================
## Step 03 - Double-sampling parameter grid
## ============================================================

double_param_grid <- expand.grid(
  n_a     = n_a_vec,
  n_b     = n_b_vec,
  alpha_a = alpha_a_vec
)

double_param_grid$alpha_b <- alpha_global_h0 / double_param_grid$alpha_a

double_param_grid$valid_alpha_b <- (
  is.finite(double_param_grid$alpha_b) &
    double_param_grid$alpha_b > 0 &
    double_param_grid$alpha_b < 1
)

double_param_grid <- double_param_grid[double_param_grid$valid_alpha_b, ]

double_param_grid$Avg_Sample0_per_cycle <- (
  double_param_grid$n_a +
    double_param_grid$n_b * double_param_grid$alpha_a
)

double_param_grid$ASS0_in_range <- (
  double_param_grid$Avg_Sample0_per_cycle >= ASS_min &
    double_param_grid$Avg_Sample0_per_cycle <= ASS_max
)

double_param_grid_all <- double_param_grid

double_param_grid <- double_param_grid[
  double_param_grid$ASS0_in_range,
]

double_param_grid$ARL0_global <- 1 / (
  double_param_grid$alpha_a * double_param_grid$alpha_b
)

double_param_grid$n_max_per_cycle <- (
  double_param_grid$n_a + double_param_grid$n_b
)

double_param_grid$design_family_id <- paste0(
  "DS_na", double_param_grid$n_a,
  "_nb", double_param_grid$n_b,
  "_aa", make_alpha_label(double_param_grid$alpha_a),
  "_ab", make_alpha_label(double_param_grid$alpha_b)
)

write.csv(
  double_param_grid_all,
  "np_weibull_DS_parameter_grid_all_expected_alert.csv",
  row.names = FALSE
)

write.csv(
  double_param_grid,
  "np_weibull_DS_parameter_grid_ASS0_5_10_expected_alert.csv",
  row.names = FALSE
)


## ============================================================
## Step 04 - Double-sampling sweep
## ============================================================

double_performance_list <- list()
double_stage_design_catalog <- list()

counter <- 1
design_counter <- 1

for (g in seq_len(nrow(double_param_grid))) {
  
  n_a_g     <- double_param_grid$n_a[g]
  n_b_g     <- double_param_grid$n_b[g]
  alpha_a_g <- double_param_grid$alpha_a[g]
  alpha_b_g <- double_param_grid$alpha_b[g]
  
  Avg_Sample0_g <- double_param_grid$Avg_Sample0_per_cycle[g]
  n_max_per_cycle_g <- double_param_grid$n_max_per_cycle[g]
  design_family_id_g <- double_param_grid$design_family_id[g]
  
  design_a_g <- make_stage_design(
    stage    = "A_first_sample",
    n_sample = n_a_g,
    alpha    = alpha_a_g,
    k0       = k0,
    lambda0  = lambda0
  )
  
  design_b_g <- make_stage_design(
    stage    = "B_confirmation_sample",
    n_sample = n_b_g,
    alpha    = alpha_b_g,
    k0       = k0,
    lambda0  = lambda0
  )
  
  design_a_export <- data.frame(
    design_family_id = design_family_id_g,
    n_a              = n_a_g,
    n_b              = n_b_g,
    alpha_a          = alpha_a_g,
    alpha_b          = alpha_b_g,
    n_max_per_cycle  = n_max_per_cycle_g,
    design_a_g,
    stringsAsFactors = FALSE
  )
  
  design_b_export <- data.frame(
    design_family_id = design_family_id_g,
    n_a              = n_a_g,
    n_b              = n_b_g,
    alpha_a          = alpha_a_g,
    alpha_b          = alpha_b_g,
    n_max_per_cycle  = n_max_per_cycle_g,
    design_b_g,
    stringsAsFactors = FALSE
  )
  
  double_stage_design_catalog[[design_counter]] <- rbind(
    design_a_export,
    design_b_export
  )
  
  design_counter <- design_counter + 1
  
  design_pairs_g <- make_design_pairs(
    design_a = design_a_g,
    design_b = design_b_g
  )
  
  design_pairs_g$rule_id <- paste0(
    design_family_id_g,
    "_Ya", design_pairs_g$Y_a,
    "_Yb", design_pairs_g$Y_b
  )
  
  for (s in seq_len(nrow(scenario_grid))) {
    
    k1_s      <- scenario_grid$k1[s]
    lambda1_s <- scenario_grid$lambda1[s]
    
    mean_1_s <- lambda1_s * gamma(1 + 1 / k1_s)
    
    var_1_s <- lambda1_s^2 * (
      gamma(1 + 2 / k1_s) - gamma(1 + 1 / k1_s)^2
    )
    
    mean_shift_direction_s <- get_mean_shift_direction(
      mean_1 = mean_1_s,
      mean_0 = weibull_mean0,
      tol    = tol_mean
    )
    
    tmp <- design_pairs_g
    
    if (mean_shift_direction_s == "increase") {
      
      chart_side_used_s <- "high"
      
      signal_rule_s <- paste(
        "Stage A: N_high_a >= Y_a;",
        "Stage B: N_high_b >= Y_b;",
        "final signal only if both stages signal"
      )
      
      p_h1_a <- pweibull(
        tmp$tcrit_high_a,
        shape = k1_s,
        scale = lambda1_s,
        lower.tail = FALSE
      )
      
      p_h1_b <- pweibull(
        tmp$tcrit_high_b,
        shape = k1_s,
        scale = lambda1_s,
        lower.tail = FALSE
      )
      
      tcrit_a <- tmp$tcrit_high_a
      tcrit_b <- tmp$tcrit_high_b
      
    } else {
      
      chart_side_used_s <- "low"
      
      if (mean_shift_direction_s == "decrease") {
        signal_rule_s <- paste(
          "Stage A: N_low_a >= Y_a;",
          "Stage B: N_low_b >= Y_b;",
          "final signal only if both stages signal"
        )
      } else {
        signal_rule_s <- paste(
          "same mean; low-side convention used;",
          "final signal only if both stages signal"
        )
      }
      
      p_h1_a <- pweibull(
        tmp$tcrit_low_a,
        shape = k1_s,
        scale = lambda1_s
      )
      
      p_h1_b <- pweibull(
        tmp$tcrit_low_b,
        shape = k1_s,
        scale = lambda1_s
      )
      
      tcrit_a <- tmp$tcrit_low_a
      tcrit_b <- tmp$tcrit_low_b
    }
    
    alpha1_a <- pbinom(
      tmp$Y_a - 1,
      size = n_a_g,
      prob = p_h1_a,
      lower.tail = FALSE
    )
    
    alpha1_b <- pbinom(
      tmp$Y_b - 1,
      size = n_b_g,
      prob = p_h1_b,
      lower.tail = FALSE
    )
    
    alpha1_global_h1 <- alpha1_a * alpha1_b
    
    ARL1 <- ifelse(
      alpha1_global_h1 > 0,
      1 / alpha1_global_h1,
      Inf
    )
    
    Avg_Sample1_per_cycle <- n_a_g + n_b_g * alpha1_a
    
    P_stop_after_A_h0 <- 1 - alpha_a_g
    P_stop_after_A_h1 <- 1 - alpha1_a
    
    tcritsum <- tcrit_a + tcrit_b
    
    Expected_Alert_Time_H0_per_cycle <- tcrit_a + alpha_a_g * tcrit_b
    Expected_Alert_Time_H1_per_cycle <- tcrit_a + alpha1_a * tcrit_b
    
    Expected_Alert_Time_To_Signal_H1 <- (
      ARL1 * Expected_Alert_Time_H1_per_cycle
    )
    
    tmp$model                 <- "double_sampling"
    tmp$design_family_id      <- design_family_id_g
    tmp$scenario_id           <- scenario_grid$scenario_id[s]
    tmp$k1                    <- k1_s
    tmp$lambda1               <- lambda1_s
    tmp$mean_weibull_0        <- weibull_mean0
    tmp$var_weibull_0         <- weibull_var0
    tmp$mean_weibull_1        <- mean_1_s
    tmp$var_weibull_1         <- var_1_s
    tmp$mean_shift_direction  <- mean_shift_direction_s
    tmp$chart_side_used       <- chart_side_used_s
    tmp$signal_rule           <- signal_rule_s
    
    tmp$n_a                   <- n_a_g
    tmp$n_b                   <- n_b_g
    tmp$n_max_per_cycle       <- n_max_per_cycle_g
    tmp$alpha_a_h0            <- alpha_a_g
    tmp$alpha_b_h0            <- alpha_b_g
    tmp$alpha_global_h0       <- alpha_global_h0
    tmp$ARL0_global           <- ARL0_target
    tmp$Avg_Sample0_per_cycle <- Avg_Sample0_g
    
    tmp$tcrit_a               <- tcrit_a
    tmp$tcrit_b               <- tcrit_b
    tmp$tcritsum              <- tcritsum
    
    tmp$P_stop_after_A_h0     <- P_stop_after_A_h0
    tmp$P_stop_after_A_h1     <- P_stop_after_A_h1
    
    tmp$Expected_Alert_Time_H0_per_cycle <- Expected_Alert_Time_H0_per_cycle
    tmp$Expected_Alert_Time_H1_per_cycle <- Expected_Alert_Time_H1_per_cycle
    tmp$Expected_Alert_Time_To_Signal_H1 <- Expected_Alert_Time_To_Signal_H1
    
    tmp$p_h1_a                <- p_h1_a
    tmp$p_h1_b                <- p_h1_b
    tmp$alpha1_a              <- alpha1_a
    tmp$alpha1_b              <- alpha1_b
    tmp$alpha1_global_h1      <- alpha1_global_h1
    tmp$ARL1                  <- ARL1
    tmp$Avg_Sample1_per_cycle <- Avg_Sample1_per_cycle
    
    tmp$ASS0_in_range <- (
      tmp$Avg_Sample0_per_cycle >= ASS_min &
        tmp$Avg_Sample0_per_cycle <= ASS_max
    )
    
    tmp$ASS1_in_range <- (
      tmp$Avg_Sample1_per_cycle >= ASS_min &
        tmp$Avg_Sample1_per_cycle <= ASS_max
    )
    
    double_performance_list[[counter]] <- tmp
    counter <- counter + 1
  }
}

double_performance_table <- do.call(rbind, double_performance_list)
double_stage_design_catalog <- do.call(rbind, double_stage_design_catalog)


## ============================================================
## Step 05 - np2s vs Weibull stage design
## ============================================================

np2s_vs_weibull_stage_design <- double_stage_design_catalog

np2s_vs_weibull_stage_design$stage_type <- ifelse(
  np2s_vs_weibull_stage_design$stage == "A_first_sample",
  "Stage A",
  "Stage B"
)

np2s_vs_weibull_stage_design$P_signal_geq_Y_h0 <- pbinom(
  np2s_vs_weibull_stage_design$Y - 1,
  size = np2s_vs_weibull_stage_design$n_sample,
  prob = np2s_vs_weibull_stage_design$p1,
  lower.tail = FALSE
)

np2s_vs_weibull_stage_design$ARL_stage_if_alone <- ifelse(
  np2s_vs_weibull_stage_design$P_signal_geq_Y_h0 > 0,
  1 / np2s_vs_weibull_stage_design$P_signal_geq_Y_h0,
  Inf
)

np2s_vs_weibull_stage_design <- np2s_vs_weibull_stage_design[, c(
  "design_family_id",
  "stage_type",
  "stage",
  "n_a",
  "n_b",
  "n_max_per_cycle",
  "alpha_a",
  "alpha_b",
  "n_sample",
  "alpha_stage",
  "Y",
  "p1",
  "tcrit_low",
  "tcrit_high",
  "check_p1_low_from_weibull0",
  "check_p1_high_from_weibull0",
  "check_alpha_low_from_binom",
  "check_alpha_high_from_binom",
  "P_signal_geq_Y_h0",
  "ARL_stage_if_alone",
  "err_alpha_low",
  "err_alpha_high"
)]

write.csv(
  np2s_vs_weibull_stage_design,
  "np2s_vs_weibull_stage_design.csv",
  row.names = FALSE
)


## ============================================================
## Step 06 - Simple summaries for Stage A and Stage B
## ============================================================

stage_A_simple <- np2s_vs_weibull_stage_design[
  np2s_vs_weibull_stage_design$stage == "A_first_sample",
]

stage_A_simple <- unique(data.frame(
  n_stage     = stage_A_simple$n_a,
  alpha_stage = stage_A_simple$alpha_a,
  UCL         = stage_A_simple$Y,
  p_y         = stage_A_simple$p1,
  t_y         = stage_A_simple$tcrit_low,
  t_w         = stage_A_simple$tcrit_high,
  stringsAsFactors = FALSE
))

stage_A_simple <- stage_A_simple[order(
  stage_A_simple$n_stage,
  stage_A_simple$alpha_stage,
  stage_A_simple$UCL
), ]

stage_B_simple <- np2s_vs_weibull_stage_design[
  np2s_vs_weibull_stage_design$stage == "B_confirmation_sample",
]

stage_B_simple <- unique(data.frame(
  n_stage     = stage_B_simple$n_b,
  alpha_stage = stage_B_simple$alpha_b,
  UCL         = stage_B_simple$Y,
  p_y         = stage_B_simple$p1,
  t_y         = stage_B_simple$tcrit_low,
  t_w         = stage_B_simple$tcrit_high,
  stringsAsFactors = FALSE
))

stage_B_simple <- stage_B_simple[order(
  stage_B_simple$n_stage,
  stage_B_simple$alpha_stage,
  stage_B_simple$UCL
), ]

write.csv(
  stage_A_simple,
  "np2s_stage_A_simple_summary.csv",
  row.names = FALSE
)

write.csv(
  stage_B_simple,
  "np2s_stage_B_simple_summary.csv",
  row.names = FALSE
)



## ============================================================
## Step 08 - np2s rule pairs
## ============================================================

np2s_vs_weibull_rule_pairs <- unique(double_performance_table[, c(
  "design_family_id",
  "rule_id",
  "n_a",
  "n_b",
  "n_max_per_cycle",
  "Y_a",
  "Y_b",
  "alpha_a_h0",
  "alpha_b_h0",
  "alpha_global_h0",
  "ARL0_global",
  "p1_a",
  "p1_b",
  "tcrit_low_a",
  "tcrit_low_b",
  "tcrit_high_a",
  "tcrit_high_b",
  "Avg_Sample0_per_cycle"
)])

np2s_vs_weibull_rule_pairs$P_signal_A_h0 <- (
  np2s_vs_weibull_rule_pairs$alpha_a_h0
)

np2s_vs_weibull_rule_pairs$P_signal_B_h0 <- (
  np2s_vs_weibull_rule_pairs$alpha_b_h0
)

np2s_vs_weibull_rule_pairs$P_global_signal_h0 <- (
  np2s_vs_weibull_rule_pairs$P_signal_A_h0 *
    np2s_vs_weibull_rule_pairs$P_signal_B_h0
)

np2s_vs_weibull_rule_pairs$ARL0_check <- (
  1 / np2s_vs_weibull_rule_pairs$P_global_signal_h0
)

np2s_vs_weibull_rule_pairs$err_alpha_global_h0 <- (
  np2s_vs_weibull_rule_pairs$P_global_signal_h0 -
    np2s_vs_weibull_rule_pairs$alpha_global_h0
)

np2s_vs_weibull_rule_pairs <- np2s_vs_weibull_rule_pairs[order(
  np2s_vs_weibull_rule_pairs$n_a,
  np2s_vs_weibull_rule_pairs$n_b,
  np2s_vs_weibull_rule_pairs$alpha_a_h0,
  np2s_vs_weibull_rule_pairs$Y_a,
  np2s_vs_weibull_rule_pairs$Y_b
), ]

write.csv(
  np2s_vs_weibull_rule_pairs,
  "np2s_vs_weibull_rule_pairs.csv",
  row.names = FALSE
)


## ============================================================
## Step 09 - Full binomial probability tables
## ============================================================

np2s_probabilities_list <- vector(
  "list",
  nrow(np2s_vs_weibull_stage_design)
)

for (i in seq_len(nrow(np2s_vs_weibull_stage_design))) {
  
  row_i <- np2s_vs_weibull_stage_design[i, ]
  
  binom_table <- make_binom_table(
    p        = row_i$p1,
    n_sample = row_i$n_sample
  )
  
  P_signal <- pbinom(
    row_i$Y - 1,
    size = row_i$n_sample,
    prob = row_i$p1,
    lower.tail = FALSE
  )
  
  np2s_probabilities_list[[i]] <- data.frame(
    design_family_id  = row_i$design_family_id,
    stage_type        = row_i$stage_type,
    stage             = row_i$stage,
    n_a               = row_i$n_a,
    n_b               = row_i$n_b,
    n_max_per_cycle   = row_i$n_max_per_cycle,
    alpha_a           = row_i$alpha_a,
    alpha_b           = row_i$alpha_b,
    n_sample          = row_i$n_sample,
    alpha_stage       = row_i$alpha_stage,
    Y                 = row_i$Y,
    p1                = row_i$p1,
    tcrit_low         = row_i$tcrit_low,
    tcrit_high        = row_i$tcrit_high,
    P_signal_geq_Y    = P_signal,
    signal_if_N_geq_Y = binom_table$N >= row_i$Y,
    binom_table,
    stringsAsFactors = FALSE
  )
}

np2s_probabilities_stage_tables <- do.call(
  rbind,
  np2s_probabilities_list
)

write.csv(
  np2s_probabilities_stage_tables,
  "np2s_probabilities_stage_tables.csv",
  row.names = FALSE
)


## ============================================================
## Step 10 - ASS filtered performance
## ============================================================

double_performance_table_ASS_5_10 <- double_performance_table[
  double_performance_table$ASS1_in_range,
]

write.csv(
  double_stage_design_catalog,
  "np_weibull_DS_stage_design_catalog_expected_alert.csv",
  row.names = FALSE
)

write.csv(
  double_performance_table,
  "np_weibull_DS_sweep_performance_all_expected_alert.csv",
  row.names = FALSE
)

write.csv(
  double_performance_table_ASS_5_10,
  "np_weibull_DS_sweep_performance_ASS_5_10_expected_alert.csv",
  row.names = FALSE
)


## ============================================================
## Step 11 - Mapped scenarios with ASS in [5, 10]
## ============================================================

mapped_scenarios_ASS_5_10 <- double_performance_table_ASS_5_10[, c(
  "scenario_id",
  "k1",
  "lambda1",
  "mean_weibull_1",
  "var_weibull_1",
  "mean_shift_direction",
  "chart_side_used",
  "n_a",
  "n_b",
  "n_max_per_cycle",
  "Y_a",
  "Y_b",
  "alpha_a_h0",
  "alpha_b_h0",
  "Avg_Sample0_per_cycle",
  "Avg_Sample1_per_cycle",
  "P_stop_after_A_h0",
  "P_stop_after_A_h1",
  "ARL1",
  "alpha1_global_h1",
  "Expected_Alert_Time_H1_per_cycle",
  "Expected_Alert_Time_To_Signal_H1",
  "tcrit_a",
  "tcrit_b",
  "tcritsum"
)]

mapped_scenarios_ASS_5_10 <- mapped_scenarios_ASS_5_10[order(
  mapped_scenarios_ASS_5_10$scenario_id,
  mapped_scenarios_ASS_5_10$ARL1,
  mapped_scenarios_ASS_5_10$Expected_Alert_Time_H1_per_cycle,
  mapped_scenarios_ASS_5_10$Avg_Sample1_per_cycle
), ]

write.csv(
  mapped_scenarios_ASS_5_10,
  "np_weibull_DS_mapped_scenarios_ASS_5_10_expected_alert.csv",
  row.names = FALSE
)


## ============================================================
## Step 12 - Double-sampling ranking
## ============================================================

double_rule_cols <- c(
  "model",
  "rule_id",
  "design_family_id",
  "n_a",
  "n_b",
  "n_max_per_cycle",
  "Y_a",
  "Y_b",
  "alpha_a_h0",
  "alpha_b_h0",
  "alpha_global_h0",
  "ARL0_global",
  "Avg_Sample0_per_cycle"
)

double_ranking <- summarise_by_rule(
  df        = double_performance_table_ASS_5_10,
  rule_cols = double_rule_cols
)

if (nrow(double_ranking) > 0) {
  
  double_ranking <- double_ranking[order(
    double_ranking$ARL1_mean,
    double_ranking$Expected_Alert_Time_H1_mean,
    double_ranking$Expected_Alert_Time_To_Signal_H1_mean,
    double_ranking$Avg_Sample1_mean,
    double_ranking$tcritsum_mean
  ), ]
  
  rownames(double_ranking) <- NULL
}

write.csv(
  double_ranking,
  "np_weibull_DS_sweep_ranking_ASS_5_10_expected_alert.csv",
  row.names = FALSE
)


## ============================================================
## Step 13 - Best double-sampling rule per scenario
## ============================================================

best_double_per_scenario_ASS_5_10 <- get_best_per_scenario(
  double_performance_table_ASS_5_10
)

write.csv(
  best_double_per_scenario_ASS_5_10,
  "np_weibull_DS_best_rule_per_scenario_ASS_5_10_expected_alert.csv",
  row.names = FALSE
)


## ============================================================
## Step 14 - Optional Excel export
## ============================================================

if (requireNamespace("openxlsx", quietly = TRUE)) {
  
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "global_settings")
  openxlsx::addWorksheet(wb, "DS_param_all")
  openxlsx::addWorksheet(wb, "DS_param_ASS0_5_10")
  openxlsx::addWorksheet(wb, "stage_A_simple")
  openxlsx::addWorksheet(wb, "stage_B_simple")
  openxlsx::addWorksheet(wb, "np2s_vs_weibull_stage")
  openxlsx::addWorksheet(wb, "np2s_vs_weibull_pairs")
  openxlsx::addWorksheet(wb, "np2s_probabilities")
  openxlsx::addWorksheet(wb, "DS_perf_all")
  openxlsx::addWorksheet(wb, "DS_perf_ASS1_5_10")
  openxlsx::addWorksheet(wb, "DS_mapped_ASS_5_10")
  openxlsx::addWorksheet(wb, "DS_ranking_ASS_5_10")
  openxlsx::addWorksheet(wb, "best_DS_scenario")
  
  global_settings <- data.frame(
    ARL0_target      = ARL0_target,
    alpha_global_h0  = alpha_global_h0,
    ASS_min          = ASS_min,
    ASS_max          = ASS_max,
    n_a_min          = min(n_a_vec),
    n_a_max          = max(n_a_vec),
    n_b_min          = min(n_b_vec),
    n_b_max          = max(n_b_vec),
    selected_n_a     = selected_n_a,
    selected_alpha_a = selected_alpha_a,
    selected_n_b     = selected_n_b,
    selected_alpha_b = selected_alpha_b,
    k0               = k0,
    lambda0          = lambda0,
    weibull_mean0    = weibull_mean0,
    weibull_var0     = weibull_var0,
    n_double_configs_all = nrow(double_param_grid_all),
    n_double_configs_ASS0_5_10 = nrow(double_param_grid),
    n_double_rules_all = length(unique(double_performance_table$rule_id)),
    n_double_rules_ASS1_5_10 =
      length(unique(double_performance_table_ASS_5_10$rule_id)),
    n_scenarios      = nrow(scenario_grid),
    rows_DS_perf_all = nrow(double_performance_table),
    rows_DS_perf_ASS1_5_10 = nrow(double_performance_table_ASS_5_10),
    rows_stage_A_simple = nrow(stage_A_simple),
    rows_stage_B_simple = nrow(stage_B_simple)
  )
  
  openxlsx::writeData(wb, "global_settings", global_settings)
  openxlsx::writeData(wb, "DS_param_all", double_param_grid_all)
  openxlsx::writeData(wb, "DS_param_ASS0_5_10", double_param_grid)
  openxlsx::writeData(wb, "stage_A_simple", stage_A_simple)
  openxlsx::writeData(wb, "stage_B_simple", stage_B_simple)
  openxlsx::writeData(wb, "np2s_vs_weibull_stage", np2s_vs_weibull_stage_design)
  openxlsx::writeData(wb, "np2s_vs_weibull_pairs", np2s_vs_weibull_rule_pairs)
  openxlsx::writeData(wb, "np2s_probabilities", np2s_probabilities_stage_tables)
  openxlsx::writeData(wb, "DS_perf_all", double_performance_table)
  openxlsx::writeData(wb, "DS_perf_ASS1_5_10", double_performance_table_ASS_5_10)
  openxlsx::writeData(wb, "DS_mapped_ASS_5_10", mapped_scenarios_ASS_5_10)
  openxlsx::writeData(wb, "DS_ranking_ASS_5_10", double_ranking)
  openxlsx::writeData(wb, "best_DS_scenario", best_double_per_scenario_ASS_5_10)
  
  sheet_names <- names(wb)
  
  for (sh in sheet_names) {
    openxlsx::freezePane(wb, sh, firstRow = TRUE)
    openxlsx::setColWidths(wb, sh, cols = 1:100, widths = "auto")
  }
  
  openxlsx::saveWorkbook(
    wb,
    "np_weibull_DS_sweep_ASS_5_10_expected_alert_report.xlsx",
    overwrite = TRUE
  )
  
  cat("\nExcel file written: np_weibull_DS_sweep_ASS_5_10_expected_alert_report.xlsx\n")
  
} else {
  
  cat("\nPackage 'openxlsx' is not installed.\n")
  cat("Excel export skipped, but all CSV and TEX files were written.\n")
}


## ============================================================
## Step 15 - Final checks
## ============================================================

cat("\n============================================================\n")
cat("FINAL CHECKS\n")
cat("============================================================\n")

cat("ASS restriction:", ASS_min, "<= ASS <=", ASS_max, "\n")
cat("n_b constraint: n_b <=", n_b_max, "\n")
cat("Double configurations - all:", nrow(double_param_grid_all), "\n")
cat("Double configurations - ASS0 in [5,10]:", nrow(double_param_grid), "\n")
cat("Double rules - all:", length(unique(double_performance_table$rule_id)), "\n")
cat("Double rules - ASS1 in [5,10]:",
    length(unique(double_performance_table_ASS_5_10$rule_id)), "\n")
cat("Scenarios:", nrow(scenario_grid), "\n")
cat("Double performance rows - all:", nrow(double_performance_table), "\n")
cat("Double performance rows - ASS1 in [5,10]:",
    nrow(double_performance_table_ASS_5_10), "\n")

cat("\nSimple summaries:\n")
cat("Stage A simple rows:", nrow(stage_A_simple), "\n")
cat("Stage B simple rows:", nrow(stage_B_simple), "\n")

cat("\nSelected tables:\n")
cat("Stage A: n_a =", selected_n_a, ", alpha_a =", selected_alpha_a, "\n")
cat("Stage B: n_b =", selected_n_b, ", alpha_b =", selected_alpha_b, "\n")

cat("\nFiles written:\n")
cat("- np_weibull_DS_parameter_grid_all_expected_alert.csv\n")
cat("- np_weibull_DS_parameter_grid_ASS0_5_10_expected_alert.csv\n")
cat("- np2s_vs_weibull_stage_design.csv\n")
cat("- np2s_vs_weibull_rule_pairs.csv\n")
cat("- np2s_probabilities_stage_tables.csv\n")
cat("- np2s_stage_A_simple_summary.csv\n")
cat("- np2s_stage_B_simple_summary.csv\n")
cat("- np2s_stage_A_tables_mean_decrease.tex\n")
cat("- np2s_stage_B_tables_mean_decrease.tex\n")
cat("- np2s_stage_A_tables_mean_increase.tex\n")
cat("- np2s_stage_B_tables_mean_increase.tex\n")
cat("- np2s_selected_stage_A_decrease.tex\n")
cat("- np2s_selected_stage_B_decrease.tex\n")
cat("- np2s_selected_stage_A_increase.tex\n")
cat("- np2s_selected_stage_B_increase.tex\n")
cat("- np_weibull_DS_stage_design_catalog_expected_alert.csv\n")
cat("- np_weibull_DS_sweep_performance_all_expected_alert.csv\n")
cat("- np_weibull_DS_sweep_performance_ASS_5_10_expected_alert.csv\n")
cat("- np_weibull_DS_mapped_scenarios_ASS_5_10_expected_alert.csv\n")
cat("- np_weibull_DS_sweep_ranking_ASS_5_10_expected_alert.csv\n")
cat("- np_weibull_DS_best_rule_per_scenario_ASS_5_10_expected_alert.csv\n")
cat("- np_weibull_DS_sweep_ASS_5_10_expected_alert_report.xlsx, only if openxlsx is installed\n")

