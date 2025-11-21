#!/usr/bin/env Rscript
# Part 2.1-2.5 analysis for HW3 using HW3_IV.xlsx

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ivreg)
  library(broom)
})

data_path <- "HW3_IV.xlsx"
stopifnot(file.exists(data_path))
d <- read_excel(data_path)

cat("=========== 2.1 Data Summary ===========\n")
print(summary(d))
cat("\nPairwise correlations:\n")
print(cor(d))

cat("\n=========== 2.2 OLS: y ~ x1 + x2 ===========\n")
ols_fit <- lm(y ~ x1 + x2, data = d)
print(tidy(ols_fit, conf.int = TRUE))

cat("\n=========== 2.3 IV (z1 instruments x2) ===========\n")
iv_z1 <- ivreg(y ~ x1 + x2 | x1 + z1, data = d)
print(tidy(iv_z1, conf.int = TRUE))
cat("\nModel comparison (OLS vs IV with z1):\n")
print(glance(ols_fit))
print(glance(iv_z1))

cat("\n=========== 2.4 Durbin-Wu-Hausman Test ===========\n")
beta_diff <- coef(iv_z1)["x2"] - coef(ols_fit)["x2"]
var_diff <- vcov(iv_z1)["x2", "x2"] - vcov(ols_fit)["x2", "x2"]
hausman_stat <- as.numeric(beta_diff^2 / var_diff)
hausman_p <- 1 - pchisq(hausman_stat, df = 1)
cat(sprintf("Hausman statistic (df=1): %.3f, p-value: %.4g\n", hausman_stat, hausman_p))

cat("\n=========== 2.5 IV (z1 & z2 instruments) ===========\n")
iv_z12 <- ivreg(y ~ x1 + x2 | x1 + z1 + z2, data = d)
print(tidy(iv_z12, conf.int = TRUE))

cat("\nFirst-stage diagnostics for instruments:\n")
print(summary(iv_z12, diagnostics = TRUE))
