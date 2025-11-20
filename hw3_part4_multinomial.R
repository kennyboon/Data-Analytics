#!/usr/bin/env Rscript

## Homework 3 – Part 4: Multinomial Logit and Probit

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(mlogit)
  library(mprobit)
  library(broom)
})

data_path <- file.path(getwd(), "hw_multi.xlsx")
if (!file.exists(data_path)) {
  stop("Cannot find hw_multi.xlsx at: ", data_path)
}

raw <- read_xlsx(data_path, sheet = 1, col_names = TRUE) %>%
  transmute(
    id = bifid,
    age = as.numeric(age),
    stype = factor(stype, levels = c("CT", "OJT", "JSA", "Other")),
    ed1 = as.integer(ed1),
    ed2 = as.integer(ed2),
    ed3 = as.integer(ed3),
    black = as.integer(black),
    hisp = as.integer(hisp),
    nvrwrk = as.integer(nvrwrk)
  )

cat("Q4.1 – Sample counts by service type:\n")
print(count(raw, stype))

covars <- c("age", "ed2", "ed3", "black", "hisp", "nvrwrk")

ml_data <- mlogit.data(
  data = raw,
  choice = "stype",
  shape = "wide",
  id.var = "id",
  alt.levels = levels(raw$stype)
)

ml_formula <- stype ~ 0 | age + ed2 + ed3 + black + hisp + nvrwrk
mlogit_fit <- mlogit(ml_formula, data = ml_data, reflevel = "CT")

cat("\nQ4.1 – Multinomial logit coefficients (baseline = CT):\n")
print(summary(mlogit_fit))

mlogit_tidy <- tidy(mlogit_fit) %>%
  mutate(model = "mlogit")

mp_formula <- stype ~ age + ed2 + ed3 + black + hisp + nvrwrk
mprobit_fit <- mprobit(
  formula = mp_formula,
  data = raw,
  reflevel = "CT",
  link = "probit",
  draws = 2000
)

cat("\nQ4.2 – Multinomial probit coefficients:\n")
print(summary(mprobit_fit))

mprobit_tidy <- tidy(mprobit_fit) %>%
  mutate(model = "mprobit")

coef_tbl <- bind_rows(mlogit_tidy, mprobit_tidy)
write.csv(coef_tbl, "part4_r_coef_summary.csv", row.names = FALSE)

pred_profiles <- tibble(
  age = c(25, 35, 45),
  ed2 = c(1, 0, 0),
  ed3 = c(0, 1, 0),
  black = 0,
  hisp = c(0, 0, 1),
  nvrwrk = c(0, 0, 1),
  id = 999001:999003
) %>%
  crossing(stype = levels(raw$stype))

pred_data <- mlogit.data(pred_profiles, choice = "stype", shape = "long", chid.var = "id", alt.levels = levels(raw$stype))
pred_probs <- predict(mlogit_fit, newdata = pred_data)
cat("\nIllustrative predicted probabilities from mlogit for three profiles:\n")
print(pred_probs)

cat("\nScript complete – see console output and CSV summaries for reporting.\n")
