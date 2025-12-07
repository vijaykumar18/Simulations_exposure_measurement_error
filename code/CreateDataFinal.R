####Association data########
# Set working directory and load libraries
setwd("~/Downloads/CodeReview/CodeReview_Heather")
library(survival)
library(ggplot2)
library(broom)
library(haven)
library(dplyr)
library(tidyr)
library(readr)
library(MASS)
library(lubridate)
library(lme4)
rm(list = ls())



# Load PM2.5 data
daily_pm25 <- read_csv("pm25_js_zip_2000_2016_ny_daily.csv")
daily_pm25$date2 <- as.Date(as.character(daily_pm25$date), format = "%Y%m%d")
daily_pm25$date <- daily_pm25$date2
daily_pm25 <- daily_pm25%>%dplyr::select(date,ZIP,pm25)
daily_pm25 <- as.data.frame(daily_pm25)
head(daily_pm25)

# Load maternal and preterm birth data
Mom_zip_all <- read_sas("zipall_1719.sas7bdat")
Mom_zip_all$Maternal_all <- rowSums(Mom_zip_all[, 2:10], na.rm = TRUE)
Mom_zip_all$Preterm_all <- rowSums(Mom_zip_all[, 80:83], na.rm = TRUE)
Mom_zip <- Mom_zip_all[, c("ZIP", "Maternal_all", "Preterm_all")]
summary(Mom_zip)
hist(Mom_zip$Preterm_all)

# Filter for zip codes with births
zip_data <- Mom_zip %>% filter(Maternal_all > 0)
head(zip_data)

set.seed(212)

# Filter for zip codes with births
zip_data <- Mom_zip %>% filter(Maternal_all > 0)

# Filter for zip codes with births
zip_data <- Mom_zip %>% filter(Maternal_all > 0)

set.seed(212)

# TRUE NULL function with overdispersion
simulate_true_null_births <- function(maternal_all, n_years = 11) {
  total_days <- round(n_years * 365.25)
  
  # Use OVERALL average preterm probability for NULL data
  overall_preterm_prob <- sum(zip_data$Preterm_all) / sum(zip_data$Maternal_all)
  
  daily_births <- round(rnorm(total_days, maternal_all/365.25, 0.1*(maternal_all/365.25)))
  daily_births <- pmax(daily_births, 0)
  
  # Add overdispersion using negative binomial
  mu <- daily_births * overall_preterm_prob
  overdispersion <- 1.57  # Same as your association data
  size_param <- 1/(overdispersion - 1)
  
  simulated_daily_preterm <- rnbinom(total_days, mu = mu, size = size_param)
  simulated_daily_preterm <- pmin(pmax(round(simulated_daily_preterm), 0), daily_births)
  
  return(cbind(daily_births, simulated_daily_preterm))
}

# Create TRUE NULL data
simulated_data <- data.frame()
common_zips <- unique(zip_data$ZIP)

for (zip_code in common_zips) {
  zip_subset <- zip_data %>% filter(ZIP == zip_code)
  
  total_days <- round(11 * 365.25)
  simulated_births <- simulate_true_null_births(zip_subset$Maternal_all, 11)
  
  zip_simulated_data <- data.frame(
    ZIP = rep(zip_code, total_days),
    date = seq(as.Date("2006-01-01"), by = "day", length.out = total_days),
    preterm_births = simulated_births[,2],
    all_births = simulated_births[,1]
  )
  
  simulated_data <- bind_rows(simulated_data, zip_simulated_data)
}

# Merge with daily PM2.5 data
truedta_null <- merge(simulated_data, daily_pm25, by = c("ZIP", "date"))
truedta_null <- truedta_null %>% filter(all_births > 0)
truedta_null$pm25<-truedta_null$pm25/10

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

# Parameters for association data
base_prev <- 0.09
target_rr <- 1.115
beta1 <- log(target_rr)

# Use the same base structure but create association
truedta_ass <- truedta_null %>%  # Start with the same structure
  dplyr::select(-preterm_births)

# Calculate intercept
mean_pm25 <- mean(truedta_ass$pm25, na.rm = TRUE)
beta0 <- log(base_prev) - beta1 * mean_pm25

# Generate outcomes with association and overdispersion
mu <- truedta_ass$all_births * exp(beta0 + beta1 * truedta_ass$pm25)

# Add overdispersion using negative binomial
overdispersion <- 1.57
size_param <- 1/(overdispersion - 1)

set.seed(212)
outc <- rnbinom(nrow(truedta_ass), mu = mu, size = size_param)
outc <- pmin(pmax(round(outc), 0), truedta_ass$all_births)

# Create final association dataset
truedta_ass <- truedta_ass %>%
  mutate(preterm_births = outc)

# Save association data  
write_csv(truedta_ass, "truedta_ass_Final.csv")

#Verify association data
# we want to make sure the simulated association data has True RR = 1.115, also association prevalence is higher 
#than null prevalence

mod_ass <- glm(preterm_births ~ pm25 + offset(log(all_births)),
               data = truedta_ass, family = quasipoisson)

cat("=== Association DATA VERIFICATION ===\n")
cat("Sample size:", nrow(truedta_ass), "observations\n")
cat("Overall preterm probability used:", round(sum(zip_data$Preterm_all) / sum(zip_data$Maternal_all), 4), "\n")
cat("Simulated null prevalence:", round(sum(truedta_null$preterm_births) / sum(truedta_null$all_births), 4), "\n")
cat("Simulated association prevalence:", round(sum(truedta_ass$preterm_births) / sum(truedta_ass$all_births), 4), "\n")
cat("Association RR:", round(exp(coef(mod_ass)["pm25"]), 4), "\n")
cat("Association RR 95% CI:", round(exp(confint(mod_ass)["pm25",]), 4), "\n")
cat("Dispersion parameter:", summary(mod_ass)$dispersion, "\n")
