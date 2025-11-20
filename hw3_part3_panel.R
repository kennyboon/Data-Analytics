#!/usr/bin/env Rscript

## Homework 3 – Part 3: Manager Pay Panel (FE vs. RE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(plm)
  library(broom)
  library(purrr)
})

data_path <- file.path(getwd(), "ManagerPay.xls")
if (!file.exists(data_path)) {
  stop("Cannot find ManagerPay.xls at: ", data_path)
}

raw <- read_xls(data_path, sheet = 1, col_names = TRUE)

rename_years <- setNames(
  c("Year2004", "Year2005", "Year2006", "Year2007"),
  names(raw)[1:4]
)

rename_sic <- setNames(
  paste0("SIC", sprintf("%02d", 1:45)),
  names(raw)[5:49]
)

rename_tail <- setNames(
  c("FirmNumber", "CEO_Bonus", "CEO_Options", "CEO_Salary", "CEO_Total", "MarketValue"),
  names(raw)[50:55]
)

d <- raw %>%
  rename(any_of(rename_years)) %>%
  rename(any_of(rename_sic)) %>%
  rename(any_of(rename_tail)) %>%
  mutate(
    Year = Year2004 * 2004 + Year2005 * 2005 + Year2006 * 2006 + Year2007 * 2007,
    mv_scaled = `MarketValue` / 1e7
  )

cat("Q3.1 – Sample overview:\n")
cat(sprintf("Observations: %d; Firms: %d; Years: %s\n",
            nrow(d),
            dplyr::n_distinct(d$FirmNumber),
            paste(range(d$Year), collapse = "-")))
print(
  d %>%
    summarise(
      CEO_Total_mean = mean(`CEO_Total`),
      CEO_Total_sd = sd(`CEO_Total`),
      mv_scaled_mean = mean(mv_scaled),
      mv_scaled_sd = sd(mv_scaled)
    )
)

formula_panel <- as.formula(`CEO_Total` ~ Year2005 + Year2006 + Year2007 + mv_scaled)
pdata <- pdata.frame(d, index = c("FirmNumber", "Year"))

fixed_mod <- plm(formula_panel, data = pdata, model = "within", effect = "individual")
random_mod <- plm(formula_panel, data = pdata, model = "random", effect = "individual")

tidy_with_ci <- function(model, label) {
  tidy(model, conf.int = TRUE) %>%
    mutate(model = label) %>%
    select(model, everything())
}

results_tbl <- bind_rows(
  tidy_with_ci(fixed_mod, "Fixed Effects"),
  tidy_with_ci(random_mod, "Random Effects")
)

cat("\nQ3.1 & Q3.2 – Coefficient estimates:\n")
print(results_tbl)

haus <- phtest(fixed_mod, random_mod)
cat("\nQ3.2 – Hausman test (RE vs. FE):\n")
print(haus)

cat("\nModel diagnostics:\n")
print(summary(fixed_mod))
print(summary(random_mod))

cat("\nInterpretation notes:\n")
cat("- Fixed effects absorb firm-level heterogeneity and attribute inter-temporal variation in CEO pay to year dummies and market value changes.\n")
cat("- Random effects assumes firm-specific effects are uncorrelated with covariates; Hausman test evaluates this.\n")
cat("- Compare the pay-performance slope (mv_scaled) under both estimators to gauge alignment with agency theory’s unit elasticity benchmark.\n")

cat("\nScript complete—use printed tables for the write-up.\n")
