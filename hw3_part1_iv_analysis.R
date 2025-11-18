#!/usr/bin/env Rscript

## Homework 3 – Part 1: Instrumental Variables on the soft drinks data

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

safe_solve <- function(mat, label = "matrix") {
  tryCatch(
    solve(mat),
    error = function(e) {
      stop(sprintf("Matrix inversion failed for %s: %s", label, e$message), call. = FALSE)
    }
  )
}

run_ols <- function(X, Y) {
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  beta <- safe_solve(XtX, "X'X") %*% XtY
  resid <- Y - X %*% beta
  df <- nrow(X) - ncol(X)
  sigma2 <- as.numeric(crossprod(resid) / df)
  cov <- sigma2 * safe_solve(XtX, "Var(OLS)")
  se <- sqrt(diag(cov))
  tidy <- tibble(
    term = colnames(X),
    estimate = as.numeric(beta),
    std.error = se,
    statistic = estimate / std.error,
    p.value = 2 * pt(-abs(statistic), df = df)
  )
  list(beta = beta, residuals = resid, sigma2 = sigma2, cov = cov, df = df, tidy = tidy)
}

projection_matrix <- function(W) {
  solve_W <- safe_solve(crossprod(W), "W'W")
  W %*% solve_W %*% t(W)
}

run_iv <- function(X, Y, W, P = NULL) {
  if (is.null(P)) {
    P <- projection_matrix(W)
  }
  XtPX <- t(X) %*% P %*% X
  XtPY <- t(X) %*% P %*% Y
  beta <- safe_solve(XtPX, "X'PX") %*% XtPY
  resid <- Y - X %*% beta
  df <- nrow(X) - ncol(X)
  sigma2 <- as.numeric(crossprod(resid) / df)
  cov <- sigma2 * safe_solve(XtPX, "Var(IV)")
  se <- sqrt(diag(cov))
  tidy <- tibble(
    term = colnames(X),
    estimate = as.numeric(beta),
    std.error = se,
    statistic = estimate / std.error,
    p.value = 2 * pt(-abs(statistic), df = df)
  )
  list(beta = beta, residuals = resid, sigma2 = sigma2, cov = cov, df = df, tidy = tidy, P = P)
}

hausman_test <- function(beta_ols, beta_iv, cov_ols, cov_iv, df = 1) {
  diff <- beta_iv - beta_ols
  middle <- cov_iv - cov_ols
  stat <- tryCatch(
    as.numeric(t(diff) %*% safe_solve(middle, "Hausman middle") %*% diff),
    error = function(e) {
      warning(sprintf("Hausman statistic not available: %s", e$message))
      NA_real_
    }
  )
  list(
    statistic = stat,
    df = df,
    p.value = ifelse(is.na(stat), NA_real_, 1 - pchisq(stat, df)),
    critical_5pct = qchisq(0.95, df)
  )
}

overid_test <- function(resid, P, sigma2, df) {
  stat <- as.numeric(t(resid) %*% P %*% resid / sigma2)
  list(
    statistic = stat,
    df = df,
    p.value = 1 - pchisq(stat, df),
    critical_5pct = qchisq(0.95, df)
  )
}

data_path <- file.path(getwd(), "soft drinks.xls")
if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path)
}

soft <- read_xls(
  path = data_path,
  sheet = "soft drinks data",
  skip = 1,
  col_names = TRUE,
  .name_repair = "minimal"
) %>%
  mutate(across(everything(), as.numeric))

if (anyNA(soft)) {
  dropped <- sum(!complete.cases(soft))
  warning(sprintf("Dropping %d rows that contain missing values.", dropped))
  soft <- soft %>% drop_na()
}

names(soft) <- make.names(names(soft), unique = TRUE)

N <- nrow(soft)
cat(sprintf("Loaded %d observations and %d variables from %s\n", N, ncol(soft), data_path))

cat("\nQ1 – Why lagged prices might work (or fail) as IVs:\n")
cat(paste(
  "- Prices in packaged goods tend to adjust slowly because trade promotions are scheduled weeks ahead, so lagged prices predict current prices well (relevance).\n",
  "- Weekly demand shocks (weather, local events) are short-lived relative to four-week lags, so past prices predate the current error term if shocks are transient (exclusion).\n",
  "- Violations arise when promotions span multiple weeks or when retailers react to anticipated demand shocks, which induces correlation between lagged prices and the current error.\n",
  "- Serially correlated demand shocks or competitive pricing responses can therefore break the exclusion restriction, motivating the Hausman and over-identification tests below.\n",
  sep = ""
))

y_col <- make.names("log(Sales)")
own_price_col <- make.names("own log(price/ medianprice)")
lag_cols <- make.names(c(
  "own price of past 4 weeks",
  "own price of past 3 weeks",
  "own price of past 2 weeks",
  "own price of past week"
))

exo_idx <- c(2:28, 30:43)
X_exog <- as.matrix(soft[, exo_idx])
colnames(X_exog) <- make.names(colnames(soft)[exo_idx], unique = TRUE)

own_price <- as.matrix(soft[, own_price_col, drop = FALSE])
colnames(own_price) <- "own_log_price"

lag_prices <- as.matrix(soft[, lag_cols])
colnames(lag_prices) <- c("lag_price_4w", "lag_price_3w", "lag_price_2w", "lag_price_1w")

Y <- as.matrix(soft[, y_col, drop = FALSE])
colnames(Y) <- "log_sales"

X <- cbind(X_exog, own_price, intercept = 1)
colnames(X)[ncol(X)] <- "intercept"

W <- cbind(X_exog, lag_prices, intercept = 1)
colnames(W)[ncol(W)] <- "intercept"

ols <- run_ols(X, Y)
iv <- run_iv(X, Y, W)

comparison <- ols$tidy %>%
  select(term, estimate_ols = estimate, std.error_ols = std.error, p.value_ols = p.value) %>%
  left_join(
    iv$tidy %>%
      select(term, estimate_iv = estimate, std.error_iv = std.error, p.value_iv = p.value),
    by = "term"
  )

cat("\nQ2 – Coefficient comparison (OLS vs. IV):\n")
print(comparison)

hausman <- hausman_test(ols$beta, iv$beta, ols$cov, iv$cov, df = 1)
cat("\nQ3 – Hausman (Durbin-Wu) test for exogeneity of regressors:\n")
print(hausman)

overid <- overid_test(iv$residuals, iv$P, iv$sigma2, df = ncol(W) - ncol(X))
cat("\nQ3 – Over-identification test for the instrument set (incl. lagged prices):\n")
print(overid)

stage1_df <- bind_cols(
  tibble(own_log_price = as.numeric(own_price)),
  as_tibble(X_exog),
  as_tibble(lag_prices)
)

lag_names <- colnames(lag_prices)
exo_terms <- setdiff(names(stage1_df)[-1], lag_names)
full_formula <- reformulate(termlabels = names(stage1_df)[-1], response = "own_log_price")
exo_formula <- reformulate(termlabels = exo_terms, response = "own_log_price")
stage1_full <- lm(full_formula, data = stage1_df)
stage1_restricted <- lm(exo_formula, data = stage1_df)
anova_stage1 <- anova(stage1_restricted, stage1_full)

partial_F <- anova_stage1[["F"]][2]
partial_p <- anova_stage1[["Pr(>F)"]][2]
first_stage_R2 <- summary(stage1_full)$r.squared

cat("\nFirst-stage relevance diagnostics for lagged prices:\n")
cat(sprintf("Partial F (lagged prices | other covariates): %.2f (p = %.4f)\n", partial_F, partial_p))
cat(sprintf("Full first-stage R-squared: %.3f\n", first_stage_R2))

lag_only_model <- lm(reformulate(lag_names, response = "own_log_price"), data = stage1_df)
lag_only_F <- summary(lag_only_model)$fstatistic
cat(sprintf(
  "Lag-only first-stage F: %.2f (df1 = %.0f, df2 = %.0f)\n",
  lag_only_F["value"], lag_only_F["numdf"], lag_only_F["dendf"]
))

cat("\n2SLS using only the four lagged prices as instruments (diagnostic):\n")
lag_only_W <- cbind(lag_prices, intercept = 1)
colnames(lag_only_W)[ncol(lag_only_W)] <- "intercept"

lag_only_attempt <- tryCatch(
  {
    P_lag <- projection_matrix(lag_only_W)
    XtPX_lag <- t(X) %*% P_lag %*% X
    beta_lag <- safe_solve(XtPX_lag, "lag-only X'PX") %*% t(X) %*% P_lag %*% Y
    list(success = TRUE, beta = beta_lag)
  },
  error = function(e) list(success = FALSE, message = e$message)
)

if (isTRUE(lag_only_attempt$success)) {
  warning("Lag-only IV estimation unexpectedly succeeded; inspect condition numbers manually.")
} else {
  cat("As expected, lag-only IV fails because X'P_lag X is near-singular:\n")
  cat(lag_only_attempt$message, "\n")
}

cat("\nScript complete. Use the printed diagnostics to interpret elasticity estimates and instrument validity.\n")
