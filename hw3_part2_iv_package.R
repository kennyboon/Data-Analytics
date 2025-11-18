#!/usr/bin/env Rscript

## Homework 3 – Part 2: IV practice with ivreg package

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(broom)
  library(ivreg)
  library(purrr)
})

data_path <- file.path(getwd(), "HW3_IV.xlsx")
if (!file.exists(data_path)) {
  stop("Cannot find HW3_IV.xlsx at: ", data_path)
}

df <- read_xlsx(data_path)

cat("Q2.1 – Data summary:\n")
summary_tbl <- df %>%
  summarise(across(
    everything(),
    list(mean = mean, sd = sd, min = min, max = max),
    .names = "{.col}_{.fn}"
  ))
print(summary_tbl)

model_form <- as.formula(y ~ x1 + x2)
instr_z1 <- as.formula(~ x1 + z1)
instr_z12 <- as.formula(~ x1 + z1 + z2)

ols_mod <- lm(model_form, data = df)
iv_z1 <- ivreg(model_form | instr_z1, data = df)
iv_z12 <- ivreg(model_form | instr_z12, data = df)

extract_tidy <- function(fit, label) {
  tidy(fit) %>%
    mutate(model = label) %>%
    select(model, term, estimate, std.error, statistic, p.value)
}

results_tbl <- bind_rows(
  extract_tidy(ols_mod, "OLS"),
  extract_tidy(iv_z1, "IV: z1"),
  extract_tidy(iv_z12, "IV: z1+z2")
)

cat("\nQ2.2–2.5 – Coefficient comparison:\n")
print(results_tbl)

hausman_test <- function(ols_fit, iv_fit, param = "x2") {
  beta_ols <- coef(ols_fit)[param]
  beta_iv <- coef(iv_fit)[param]
  vcov_ols <- vcov(ols_fit)[param, param, drop = FALSE]
  vcov_iv <- vcov(iv_fit)[param, param, drop = FALSE]
  diff <- beta_iv - beta_ols
  diff_var <- vcov_iv - vcov_ols
  stat <- as.numeric(diff^2 / diff_var)
  p <- 1 - pchisq(stat, df = 1)
  tibble(
    parameter = param,
    statistic = stat,
    df = 1,
    p.value = p,
    critical_5pct = qchisq(0.95, 1)
  )
}

haus_tbl <- hausman_test(ols_mod, iv_z1)
cat("\nQ2.4 – Durbin-Wu-Hausman test (x2 exogenous?):\n")
print(haus_tbl)

cat("\nWeak-IV and over-ID diagnostics from ivreg:\n")
print(summary(iv_z1, diagnostics = TRUE))
print(summary(iv_z12, diagnostics = TRUE))

first_stage <- lm(x2 ~ x1 + z1, data = df)
cat("\nFirst-stage regression (x2 on x1 + z1):\n")
print(summary(first_stage))

sargan_stat <- function(iv_fit) {
  k_instr <- length(iv_fit$instruments)
  k_reg <- length(coef(iv_fit))
  df_J <- k_instr - k_reg
  if (df_J <= 0) {
    return(tibble(statistic = NA_real_, df = df_J, p.value = NA_real_))
  }
  res <- residuals(iv_fit)
  inst_mat <- model.matrix(iv_fit$instruments, data = df)
  Pz <- inst_mat %*% solve(crossprod(inst_mat)) %*% t(inst_mat)
  sigma2 <- sum(res^2) / iv_fit$df.residual
  stat <- as.numeric(t(res) %*% Pz %*% res / sigma2)
  tibble(
    statistic = stat,
    df = df_J,
    p.value = 1 - pchisq(stat, df_J),
    critical_5pct = qchisq(0.95, df_J)
  )
}

sargan_tbl <- sargan_stat(iv_z12)
cat("\nQ2.5 – Over-identification test for instruments (z1, z2):\n")
print(sargan_tbl)

cat("\nScript complete. Use printed tables for the written answers.\n")
