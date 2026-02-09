options(digits = 16)

## ============================================
## Step 00 - Parâmetros
## ============================================

n_sample    <- 50
ARL0        <- 370
alpha_total <- 1 / ARL0
alpha_upper <- alpha_total / 2

# Weibull sob controle (H0)
k0      <- 5
lambda0 <- 4

# Varreduras
LCL_vec <- 0:(n_sample - 1)   # LCL = n é impossível (P(N<=n)=1)
UCL_vec <- 1:n_sample         # UCL = 0 é trivial (P(N>=0)=1)

# H1
k1_vec      <- c(3, 4, 5, 6, 7)
lambda1_vec <- c(1.25, 1.5, 1.75, 2, 2.25, 2.36, 2.38, 2.4, 2.43, 2.47, 2.5, 2.75,
                 3, 3.14, 3.17, 3.2, 3.24, 3.25, 3.29, 3.5, 3.73, 3.75, 3.76, 3.8,
                 3.85, 3.91, 4, 4.12, 4.16, 4.2, 4.25, 4.32, 4.5, 4.71,
                 4.75, 4.8, 4.86, 4.94, 5, 5.25, 5.497, 5.5, 5.54, 5.6,
                 5.67, 5.75, 5.76, 6)

## ============================================
## Step 01 - Funções p1 (LCL e UCL)
## ============================================

# LCL-driven (LOWER): objetivo P(N <= LCL) = alpha_upper
# identidade: P(X <= k) = I_{1-p}(n-k, k+1)
# => 1-p = qbeta(alpha, n-k, k+1)  => p = 1 - qbeta(...)
find_p1_lcl_from_LCL <- function(LCL, n, alpha_upper) {
  1 - qbeta(alpha_upper, shape1 = n - LCL, shape2 = LCL + 1)
}

# UCL-driven (UPPER): objetivo P(N >= UCL) = alpha_upper
# => P(N >= UCL) = alpha  <=> p = qbeta(alpha; UCL, n-UCL+1)
find_p1_ucl_from_UCL <- function(UCL, n, alpha_upper) {
  qbeta(alpha_upper, shape1 = UCL, shape2 = n - UCL + 1)
}

## ============================================
## Step 02 - Gerar relatório np_probabilities.csv
## ============================================

out_list <- vector("list", length(LCL_vec) * length(UCL_vec))
idx <- 1L

for (l in seq_along(LCL_vec)) {
  
  LCL    <- as.integer(LCL_vec[l])
  p1_lcl <- find_p1_lcl_from_LCL(LCL, n_sample, alpha_upper)
  
  # N serve pros dois blocos
  N <- 0:n_sample
  
  ## ----- bloco LCL (LOWER): sinal se N <= LCL
  P_signal_LCL <- pbinom(LCL, size = n_sample, prob = p1_lcl, lower.tail = TRUE)  # ~ alpha_upper
  signal_LCL   <- (N <= LCL)
  
  P_eq_N_LCL  <- dbinom(N, size = n_sample, prob = p1_lcl)
  P_leq_N_LCL <- pbinom(N, size = n_sample, prob = p1_lcl, lower.tail = TRUE)
  P_gt_N_LCL  <- pbinom(N, size = n_sample, prob = p1_lcl, lower.tail = FALSE)
  P_geq_N_LCL <- pbinom(N - 1, size = n_sample, prob = p1_lcl, lower.tail = FALSE)
  
  Always_Signal_Lower <- (N < LCL)
  
  # Checagem (tolerância numérica)
  if (abs(P_signal_LCL - alpha_upper) > 1e-10) {
    warning(sprintf("LCL=%d: P_signal_LCL=%.12g difere de alpha_upper=%.12g (dif=%.3e)",
                    LCL, P_signal_LCL, alpha_upper, P_signal_LCL - alpha_upper))
  }
  
  ## ----- varre UCL dentro de cada LCL
  for (u in seq_along(UCL_vec)) {
    
    UCL    <- as.integer(UCL_vec[u])
    p1_UCL <- find_p1_ucl_from_UCL(UCL, n_sample, alpha_upper)
    
    ## ----- bloco UCL (UPPER): sinal se N >= UCL
    P_signal_UCL <- pbinom(UCL - 1, size = n_sample, prob = p1_UCL, lower.tail = FALSE)  # ~ alpha_upper
    signal_UCL   <- (N >= UCL)
    
    P_eq_N_UCL  <- dbinom(N, size = n_sample, prob = p1_UCL)
    P_leq_N_UCL <- pbinom(N, size = n_sample, prob = p1_UCL, lower.tail = TRUE)
    P_gt_N_UCL  <- pbinom(N, size = n_sample, prob = p1_UCL, lower.tail = FALSE)
    P_geq_N_UCL <- pbinom(N - 1, size = n_sample, prob = p1_UCL, lower.tail = FALSE)
    
    # checagem UCL
    if (abs(P_signal_UCL - alpha_upper) > 1e-10) {
      warning(sprintf("UCL=%d: P_signal_UCL=%.12g difere de alpha_upper=%.12g (dif=%.3e)",
                      UCL, P_signal_UCL, alpha_upper, P_signal_UCL - alpha_upper))
    }
    
    out_list[[idx]] <- data.frame(
      n_sample = n_sample,
      p1_lcl   = p1_lcl,
      
      # LCL block
      P_signal_LCL        = P_signal_LCL,
      signal_LCL          = signal_LCL,
      LCL                 = LCL,
      N                   = N,
      P_eq_N_LCL          = P_eq_N_LCL,
      P_leq_N_LCL         = P_leq_N_LCL,
      P_gt_N_LCL          = P_gt_N_LCL,
      P_geq_N_LCL         = P_geq_N_LCL,
      Always_Signal_Lower = Always_Signal_Lower,
      
      # UCL block
      UCL          = UCL,
      p1_UCL       = p1_UCL,
      P_signal_UCL = P_signal_UCL,
      signal_UCL   = signal_UCL,
      P_eq_N_UCL   = P_eq_N_UCL,
      P_leq_N_UCL  = P_leq_N_UCL,
      P_gt_N_UCL   = P_gt_N_UCL,
      P_geq_N_UCL  = P_geq_N_UCL,
      
      check.names = FALSE
    )
    
    idx <- idx + 1L
  }
}

np_probabilities <- do.call(rbind, out_list)

# Ordenar por LCL, UCL, N
np_probabilities <- np_probabilities[order(np_probabilities$LCL,
                                           np_probabilities$UCL,
                                           np_probabilities$N), ]

write.csv(np_probabilities, "np_probabilities.csv", row.names = FALSE)

cat("\n--- OK: gerado np_probabilities.csv ---\n")
cat(sprintf("n=%d | alpha_upper=%.12g\n", n_sample, alpha_upper))


## ============================================
## Step 03 - Converter p1_lcl e p1_ucl em tempos Weibull (H0)
## Saída: np_vs_weibull_LCL_UCL.csv
## ============================================

mean_weibull_0 <- lambda0 * gamma(1 + 1 / k0)
var_weibull_0  <- lambda0^2 * (gamma(1 + 2 / k0) - gamma(1 + 1 / k0)^2)

design_step2 <- unique(np_probabilities[, c("n_sample", "LCL", "UCL", "p1_lcl", "p1_UCL")])
names(design_step2)[names(design_step2) == "p1_UCL"] <- "p1_ucl"

# Pares bilaterais válidos
design_step2 <- design_step2[design_step2$LCL < design_step2$UCL, ]

# Y
design_step2$Y <- design_step2$LCL

# Proteção numérica para qweibull
eps <- 1e-15
p_lcl_safe <- pmin(pmax(design_step2$p1_lcl, eps), 1 - eps)
p_ucl_safe <- pmin(pmax(design_step2$p1_ucl, eps), 1 - eps)

design_step2$lower_alert_point_weib <- qweibull(p_lcl_safe, shape = k0, scale = lambda0)
design_step2$upper_alert_point_weib <- qweibull(p_ucl_safe, shape = k0, scale = lambda0)

design_step2$mean_weibull_0 <- mean_weibull_0
design_step2$var_weibull_0  <- var_weibull_0

design_step2 <- design_step2[order(design_step2$LCL, design_step2$UCL), ]

design_step2 <- design_step2[, c(
  "n_sample","Y","LCL","UCL",
  "p1_lcl","lower_alert_point_weib",
  "p1_ucl","upper_alert_point_weib",
  "mean_weibull_0","var_weibull_0"
)]

write.csv(design_step2, "np_vs_weibull_LCL_UCL.csv", row.names = FALSE)

cat("\n--- OK: gerado np_vs_weibull_LCL_UCL.csv ---\n")
print(head(design_step2, 10))


## ============================================
## Step 04 - Desempenho sob H1 (k1, lambda1)
## ============================================

scenarios <- expand.grid(k1 = k1_vec, lambda1 = lambda1_vec)
perf_list <- vector("list", nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  
  k1_i      <- scenarios$k1[i]
  lambda1_i <- scenarios$lambda1[i]
  
  mean1 <- lambda1_i * gamma(1 + 1 / k1_i)
  var1  <- lambda1_i^2 * (gamma(1 + 2 / k1_i) - gamma(1 + 1 / k1_i)^2)
  
  # Weibull CDF sob H1 nos tempos do lower e do upper 
  p_oc_lower <- pweibull(design_step2$lower_alert_point_weib, shape = k1_i, scale = lambda1_i)
  p_oc_upper <- pweibull(design_step2$upper_alert_point_weib, shape = k1_i, scale = lambda1_i)
  
  # Binomial tails sob H1
  alpha_lcl <- pbinom(design_step2$LCL,     size = design_step2$n_sample, prob = p_oc_lower, lower.tail = TRUE)
  alpha_ucl <- pbinom(design_step2$UCL - 1, size = design_step2$n_sample, prob = p_oc_upper, lower.tail = FALSE)
  
  alpha_total <- alpha_lcl + alpha_ucl
  ARL1_total  <- ifelse(alpha_total > 0, 1 / alpha_total, Inf)
  
  tmp <- data.frame(
    # colunas do Step 02
    n_sample               = design_step2$n_sample,
    Y                      = design_step2$Y,
    LCL                    = design_step2$LCL,
    UCL                    = design_step2$UCL,
    p1_lcl                 = design_step2$p1_lcl,
    lower_alert_point_weib = design_step2$lower_alert_point_weib,
    p1_ucl                 = design_step2$p1_ucl,
    upper_alert_point_weib = design_step2$upper_alert_point_weib,
    mean_weibull_0         = design_step2$mean_weibull_0,
    var_weibull_0          = design_step2$var_weibull_0,
    
    # cenário H1
    k1      = k1_i,
    lambda1 = lambda1_i,
    mean1   = mean1,
    var1    = var1,
    
    # desempenho
    p_oc_upper  = p_oc_upper,
    p_oc_lower  = p_oc_lower,
    alpha_lcl   = alpha_lcl,
    alpha_ucl   = alpha_ucl,
    alpha_total = alpha_total,
    ARL1_total  = ARL1_total,
    
    check.names = FALSE
  )
  
  perf_list[[i]] <- tmp
}

desempenho <- do.call(rbind, perf_list)
desempenho <- desempenho[order(desempenho$k1, desempenho$lambda1, desempenho$LCL, desempenho$UCL), ]

write.csv(desempenho, "desempenho_np_vs_weibull_LCL_UCL.csv", row.names = FALSE)

cat("\n--- OK: gerado desempenho_np_vs_weibull_LCL_UCL.csv ---\n")
print(head(desempenho, 10))

## ============================================
## Step 05 - Relatório estilo "planilha" (5 casas decimais)
## Saída: relatorio_np_weibull.xlsx
## ============================================

# Pacotes (instala se precisar)
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("dplyr", quietly = TRUE))    install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE))    install.packages("tidyr")

library(openxlsx)
library(dplyr)
library(tidyr)

generate_excel_report <- function(design_step2,
                                  desempenho,
                                  ARL0,
                                  file = "relatorio_np_weibull.xlsx",
                                  digits = 5,
                                  LCL_only = NULL) {
  # quais LCLs vão entrar
  LCLs <- sort(unique(design_step2$LCL))
  if (!is.null(LCL_only)) {
    LCLs <- LCLs[LCLs %in% LCL_only]
  }
  
  wb <- createWorkbook()
  
  # estilos
  fmt_dec <- paste0("0.", paste(rep("0", digits), collapse = ""))
  style_dec <- createStyle(numFmt = fmt_dec)
  style_int <- createStyle(numFmt = "0")
  style_hdr <- createStyle(textDecoration = "bold")
  style_lbl <- createStyle(textDecoration = "bold", halign = "right")
  
  for (LCL_target in LCLs) {
    
    # ---- topo (design) para este LCL
    dL <- design_step2 %>%
      filter(LCL == LCL_target) %>%
      arrange(UCL)
    
    if (nrow(dL) == 0) next
    
    ucls <- dL$UCL
    m <- length(ucls)
    
    # valores (constantes no LCL)
    p1_lcl_val   <- unique(dL$p1_lcl)
    low_t_val    <- unique(dL$lower_alert_point_weib)
    
    if (length(p1_lcl_val) != 1 || length(low_t_val) != 1) {
      warning(sprintf("LCL=%d: p1_lcl ou lower_alert_point_weib não são únicos.", LCL_target))
      p1_lcl_val <- p1_lcl_val[1]
      low_t_val  <- low_t_val[1]
    }
    
    # matriz do topo (com 4 colunas vazias + bloco em colunas)
    top_mat <- rbind(
      c("", "", "", "LCL",                   rep(LCL_target, m)),
      c("", "", "", "p1_lcl",                rep(p1_lcl_val, m)),
      c("", "", "", "lower_alert_point_weib",rep(low_t_val, m)),
      c("", "", "", "UCL",                   ucls),
      c("", "", "", "p1_ucl",                dL$p1_ucl),
      c("", "", "", "upper_alert_point_weib",dL$upper_alert_point_weib),
      c("", "", "", "ARL0_total",            rep(ARL0, m))
    )
    
    # ---- matriz ARL1_total (cenários x UCL) para este LCL
    wide <- desempenho %>%
      filter(LCL == LCL_target) %>%
      select(mean1, var1, k1, lambda1, UCL, ARL1_total) %>%
      mutate(UCL = as.character(UCL)) %>%
      pivot_wider(names_from = UCL, values_from = ARL1_total) %>%
      arrange(k1, lambda1)
    
    # arredondar valores (opcional, além do formato)
    num_cols <- names(wide)[sapply(wide, is.numeric)]
    wide[num_cols] <- lapply(wide[num_cols], function(x) round(x, digits))
    
    # ---- escrever na aba
    sheet_name <- paste0("LCL_", LCL_target)
    addWorksheet(wb, sheetName = sheet_name)
    
    # Topo
    writeData(wb, sheet = sheet_name, x = top_mat, startRow = 1, startCol = 1, colNames = FALSE)
    
    # Estilo: labels em negrito (coluna 4)
    addStyle(wb, sheet = sheet_name, style = style_lbl,
             rows = 1:nrow(top_mat), cols = 4, gridExpand = TRUE, stack = TRUE)
    
    # Estilo: UCL (linha 4 do topo) inteiro
    addStyle(wb, sheet = sheet_name, style = style_int,
             rows = 4, cols = 5:(4 + m), gridExpand = TRUE, stack = TRUE)
    
    # Estilo: linhas decimais do topo (p1/tempos/ARL0)
    dec_rows <- c(2, 3, 5, 6, 7)
    addStyle(wb, sheet = sheet_name, style = style_dec,
             rows = dec_rows, cols = 5:(4 + m), gridExpand = TRUE, stack = TRUE)
    
    # Bloco de baixo (matriz)
    start_row <- nrow(top_mat) + 2
    writeData(wb, sheet = sheet_name, x = wide, startRow = start_row, startCol = 1, colNames = TRUE)
    
    # Cabeçalho do bloco de baixo em negrito
    addStyle(wb, sheet = sheet_name, style = style_hdr,
             rows = start_row, cols = 1:ncol(wide), gridExpand = TRUE, stack = TRUE)
    
    # Formato 5 casas no bloco de baixo (todas as células numéricas, exceto cabeçalho)
    addStyle(wb, sheet = sheet_name, style = style_dec,
             rows = (start_row + 1):(start_row + nrow(wide)),
             cols = 1:ncol(wide), gridExpand = TRUE, stack = TRUE)
    
    # Ajuste visual
    #freezePane(wb, sheet = sheet_name, firstRow = TRUE, firstCol = TRUE)
    setColWidths(wb, sheet = sheet_name, cols = 1:ncol(wide), widths = "auto")
  }
  
  saveWorkbook(wb, file = file, overwrite = TRUE)
  cat(sprintf("\n--- OK: relatório gerado em '%s' ---\n", file))
}

# Gerar para TODOS os LCLs existentes em design_step2
generate_excel_report(design_step2, desempenho, ARL0 = ARL0, file = "relatorio_np_weibull.xlsx", digits = 4)

