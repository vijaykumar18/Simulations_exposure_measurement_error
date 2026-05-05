################################################################################
# CreateDataFinal.R
# Purpose : Generate the two "true" datasets used in all simulation experiments:
#             1. truedta_null_Final.csv  -- null association (true RR = 1.00)
#             2. truedta_ass_Final.csv   -- harmful association (true RR = 1.115
#                                          per 10 ug/m3 PM2.5)
# Inputs  : pm25_js_zip_2000_2016_ny_daily.csv  (daily 1-km PM2.5, NYC)
#           zipall_1719.sas7bdat                 (ZIP-level PTB counts 2017-2019)
# Outputs : truedta_null_Final.csv
#           truedta_ass_Final.csv
# Authors : Vijay Kumar
# Updated : March 2026
################################################################################


# ── Libraries ─────────────────────────────────────────────────────────────────
library(haven)
library(dplyr)
library(readr)
library(MASS)
library(lubridate)
library(lme4)
rm(list = ls())

# ── Working directory ─────────────────────────────────────────────────────────
# Set this to the root of your project folder before running.
# e.g. setwd("~/Downloads/CodeReview/CodeReview_Alan")
# Alternatively, open the project in RStudio and it will be set automatically.
setwd("~/Downloads/CodeReview_Alan")

# ── 1. Load PM2.5 data ────────────────────────────────────────────────────────
# Daily ZIP-code-level PM2.5 predictions for NYC (2000-2016), 1-km ensemble model
# Source: Di et al. (2019), Harvard Dataverse

daily_pm25 <- read_csv("pm25_js_zip_2000_2016_ny_daily.csv")
daily_pm25$date2 <- as.Date(as.character(daily_pm25$date), format = "%Y%m%d")
daily_pm25$date <- daily_pm25$date2
daily_pm25 <- daily_pm25%>%dplyr::select(date,ZIP,pm25)
daily_pm25 <- as.data.frame(daily_pm25)
head(daily_pm25)


# ── 2. Load preterm birth (PTB) data ─────────────────────────────────────────
# ZIP-code-level aggregated PTB counts from NY State Vital Statistics 2017-2019.
# We assume a constant PTB rate over the 11-year simulation window (2006-2016).
Mom_zip_all <- read_sas("zipall_1719.sas7bdat")
Mom_zip_all$Maternal_all <- rowSums(Mom_zip_all[, 2:10], na.rm = TRUE)
Mom_zip_all$Preterm_all <- rowSums(Mom_zip_all[, 80:83], na.rm = TRUE)
Mom_zip <- Mom_zip_all[, c("ZIP", "Maternal_all", "Preterm_all")]
summary(Mom_zip)
hist(Mom_zip$Preterm_all)

# Filter for zip codes with births
zip_data <- Mom_zip %>% filter(Maternal_all > 0)
head(zip_data)


################################################################################
# STEP 1 — NULL ASSOCIATION DATA (true RR = 1.00)
# PTB counts are generated independently of PM2.5, so no association exists
# by design.  The overall NYC PTB prevalence (~8.4%) is preserved.
################################################################################

# ── 3a. Simulate daily birth and PTB counts for each ZIP code ─────────────────
# Daily birth counts: drawn from N(mean_daily_births, 10% CV)
# Daily PTB counts: drawn from NegBin(mu = overall_prob * daily_births, k)
# Counts constrained to [0, daily_births].

# ── 3. Simulation parameters ──────────────────────────────────────────────────
set.seed(212)  # <-- explicit seed here ensures reproducibility independent

# Filter for zip codes with births
zip_data <- Mom_zip %>% filter(Maternal_all > 0)

# ── NOTE ON OVERDISPERSION CALIBRATION ────────────────────────────────────────
# phi = 1.21 was determined through the following calibration procedure,
# conducted prior to finalizing this script:
#
# STAGE 1 — Initial estimation:
#   A provisional null dataset was simulated using phi = 1.57 (initial guess).
#   A negative binomial GLM was then fitted to the simulated daily null PTB
#   counts to estimate phi from the data:
#
#     fit_phi <- glm.nb(preterm_births ~ offset(log(all_births)),
#                       data = truedta_null_provisional)
#     theta_hat <- fit_phi$theta          # = 3.35
#     phi_hat   <- 1 + 1/theta_hat       # = 1.30
#
# STAGE 2 — Grid search:
#   phi was varied over seq(1.20, 1.60, by = 0.01) and both null and
#   association datasets were simulated at each value. Performance was
#   evaluated on four criteria:
#     - Null RR close to 1.000
#     - Association RR close to 1.115
#     - Dispersion ~ 1.0 in both datasets
#     - Consistent phi across null and association steps
#
# STAGE 3 — Multi-seed validation:
#   The best candidate (phi = 1.21) was validated across 5 random seeds
#   (212, 123, 456, 789, 1000):
#     - Mean null RR = 1.0005 (SD = 0.003)
#     - Mean ass  RR = 1.1153 (SD = 0.005)
#     - Mean null dispersion = 0.968
#     - Mean ass  dispersion = 0.964
#   Results were stable across all seeds, confirming phi = 1.21 as the
#   final value used in both null and association data generation.
# ─────────────────────────────────────────────────────────────────────────────
overdispersion <- 1.21  # single source of truth — used throughout


# TRUE NULL function with overdispersion
simulate_true_null_births <- function(maternal_all, n_years_sim = 11, 
                                      n_years_obs = 3,
                                      overdispersion = 1.21) {  # <-- argument with default
  total_days <- round(n_years_sim * 365.25)
  
  # Use OVERALL average preterm probability for NULL data
  overall_preterm_prob <- sum(zip_data$Preterm_all) / sum(zip_data$Maternal_all)
  
  # maternal_all is total births over 3 observed years (2017-2019)
  # Daily mean = (3-year total) / (days in 3 years)
  daily_mean <- maternal_all / (365.25 * n_years_obs)
  
  # Generate daily births using the correct daily rate
  daily_births <- round(rnorm(total_days, daily_mean, 0.1 * daily_mean))
  daily_births <- pmax(daily_births, 0)
  
  # Simulate daily PTB counts with overdispersion via negative binomial
  mu         <- daily_births * overall_preterm_prob
  size_param <- 1 / (overdispersion - 1)  # uses argument, not hardcoded value
  
  simulated_daily_preterm <- rnbinom(total_days, mu = mu, size = size_param)
  simulated_daily_preterm <- pmin(pmax(round(simulated_daily_preterm), 0), daily_births)
  
  return(cbind(daily_births, simulated_daily_preterm))
}

# Create TRUE NULL data
simulated_data <- data.frame()
common_zips    <- unique(zip_data$ZIP)

for (zip_code in common_zips) {
  zip_subset <- zip_data %>% filter(ZIP == zip_code)
  
  total_days       <- round(11 * 365.25)
  simulated_births <- simulate_true_null_births(
    maternal_all   = zip_subset$Maternal_all,
    n_years_sim    = 11,
    n_years_obs    = 3,
    overdispersion = overdispersion  # <-- passed from single source at top
  )
  
  zip_simulated_data <- data.frame(
    ZIP            = rep(zip_code, total_days),
    date           = seq(as.Date("2006-01-01"), by = "day", length.out = total_days),
    preterm_births = simulated_births[, 2],
    all_births     = simulated_births[, 1]
  )
  
  simulated_data <- bind_rows(simulated_data, zip_simulated_data)
}

# ── 3b. Merge with PM2.5 and finalise null dataset ────────────────────────────
# PM2.5 is scaled to per-10-ug/m3 units so the regression coefficient directly
# represents the RR per 10 ug/m3 (consistent with the target RR = 1.115).
# Merge with daily PM2.5 data
truedta_null <- merge(simulated_data, daily_pm25, by = c("ZIP", "date"))
truedta_null <- truedta_null %>% filter(all_births > 0)
truedta_null$pm25<-truedta_null$pm25/10 # rescale: 1 unit = 10 ug/m3

# Save null data
write_csv(truedta_null, "truedta_null_Final.csv")
# Verify TRUE NULL
# we want to make sure the simulated null data has True RR = 1, also prevalence is same than null prevalence
mod_null <- glm(preterm_births ~ pm25 + offset(log(all_births)),
                data = truedta_null, family = quasipoisson)

cat("===NULL DATA VERIFICATION ===\n")
cat("Sample size:", nrow(truedta_null), "observations\n")
cat("Overall preterm probability used:", round(sum(zip_data$Preterm_all) / sum(zip_data$Maternal_all), 4), "\n")
cat("Simulated prevalence:", round(sum(truedta_null$preterm_births) / sum(truedta_null$all_births), 4), "\n")
cat("Null RR:", round(exp(coef(mod_null)["pm25"]), 4), "\n")
cat("Null RR 95% CI:", round(exp(confint(mod_null)["pm25",]), 4), "\n")
cat("Dispersion parameter:", summary(mod_null)$dispersion, "\n")

############################CREATE ASSOCIATON DATA###############

################################################################################
# STEP 2 — HARMFUL ASSOCIATION DATA (true RR = 1.115 per 10 ug/m3)
# Uses the CORRECTED null data (with proper daily birth rates)
################################################################################

# ── 4a. Parameters ────────────────────────────────────────────────────────────
base_prev <- 0.09
target_rr <- 1.115
beta1 <- log(target_rr)

# Use corrected null data as base
truedta_ass <- truedta_null %>%  # This now has CORRECT birth counts
  dplyr::select(-preterm_births)

# Calculate intercept
mean_pm25 <- mean(truedta_ass$pm25, na.rm = TRUE)
beta0 <- log(base_prev) - beta1 * mean_pm25

# ── 4b. Generate PTB counts with PM2.5 association ────────────────────────────
mu <- truedta_ass$all_births * exp(beta0 + beta1 * truedta_ass$pm25)

# use overdispersion using negative binomial
size_param <- 1/(overdispersion - 1)

set.seed(212)  # <-- explicit seed here ensures reproducibility independent
outc <- rnbinom(nrow(truedta_ass), mu = mu, size = size_param)
outc <- pmin(pmax(round(outc), 0), truedta_ass$all_births)

# Create final association dataset
truedta_ass <- truedta_ass %>%
  mutate(preterm_births = outc)

# Save association data  
write_csv(truedta_ass, "truedta_ass_Final.csv")

# Verify association data
mod_ass <- glm(preterm_births ~ pm25 + offset(log(all_births)),
               data = truedta_ass, family = quasipoisson)

cat("=== Association DATA VERIFICATION ===\n")
cat("Sample size:", nrow(truedta_ass), "observations\n")
cat("Overall preterm probability used (null):", 
    round(sum(zip_data$Preterm_all) / sum(zip_data$Maternal_all), 4), "\n")
cat("Simulated null prevalence (from corrected data):", 
    round(sum(truedta_null$preterm_births) / sum(truedta_null$all_births), 4), "\n")
cat("Simulated association prevalence:", 
    round(sum(truedta_ass$preterm_births) / sum(truedta_ass$all_births), 4), "\n")
cat("Association RR:", round(exp(coef(mod_ass)["pm25"]), 4), "\n")
cat("Association RR 95% CI:", round(exp(confint(mod_ass)["pm25",]), 4), "\n")
cat("Dispersion parameter:", summary(mod_ass)$dispersion, "\n")


# ── Session info for reproducibility ─────────────────────────────────────────
cat("\n=== SESSION INFO ===\n")
sessionInfo()