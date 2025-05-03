####Association data########
# Set working directory and load libraries
setwd("/scratch/user/u.vk230134/Downloads/DiffMeasError/CodeReview")
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
daily_pm25 <- daily_pm25[, c(2, 3, 5)]
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

# Function to simulate daily preterm births
simulate_daily_preterm_births <- function(maternal_all, preterm_all, n_years = 11) {
  # Total number of days in the simulation period
  total_days <- round(n_years * 365.25)
  
  # Simulate daily preterm birth probabilities based on the ratio of preterm to total births
  preterm_prob        <- preterm_all / maternal_all
  daily_preterm_probs <- rep(preterm_prob, each = total_days)
  daily_births        <- round(rnorm(total_days, maternal_all/365.25, 0.1*(maternal_all/365.25)))
  
  # Simulate daily preterm births
  simulated_daily_preterm <- rbinom(total_days, size = daily_births, prob = daily_preterm_probs)
  
  return(cbind(daily_births, simulated_daily_preterm))
}



# Create a data frame for simulation results
simulated_data <- data.frame()
common_zips    <- unique(zip_data$ZIP)

# Loop through each ZIP code
for (zip_code in common_zips) {
  zip_subset <- zip_data %>% filter(ZIP == zip_code)
  
  # Define simulation years
  n_years_zip <- 11
  total_days <- round(n_years_zip * 365.25)
  
  # Simulate daily preterm births
  simulated_preterm <- simulate_daily_preterm_births(zip_subset$Maternal_all, zip_subset$Preterm_all, n_years_zip)
  
  # Corrected data frame: Ensure 1 row per date
  zip_simulated_data <- data.frame(
    ZIP = rep(zip_code, total_days),  # Only 1 row per date
    date = seq(as.Date("2006-01-01"), by = "day", length.out = total_days),
    preterm_births = simulated_preterm[,2],
    all_births = simulated_preterm[,1]
  )
  
  # Append the results
  simulated_data <- bind_rows(simulated_data, zip_simulated_data)
}


head(simulated_data)
# Check the first rows of the simulated data
head(simulated_data)

# Merge with daily PM2.5 data, keeping only rows with non-zero births
simulated_preterm_pm <- merge(simulated_data, daily_pm25, by = c("ZIP", "date"))
simulated_preterm_pm <- simulated_preterm_pm %>% filter(all_births > 0)

# Check merged data
head(simulated_preterm_pm)
truedta_null<-simulated_preterm_pm

# Model without assumed association
mod1 <- glm(preterm_births ~ pm25 + offset(log(all_births)), data = truedta_null, family = quasipoisson)
summary(mod1)
exp(coef(mod1)["pm25"])  # Final estimated RR
# Write the final simulated data to CSV
write_csv(truedta_null, "truedta_null_rpois.csv")

# Clear the environment
rm(list = ls())

#################### Associations
# Read the null dataset

truedta_null <- read_csv("/scratch/user/u.vk230134/Downloads/DiffMeasError/CodeReview/truedta_null_rpois.csv")

# Filter out zero all_births rows
truedta_null <- truedta_null %>% filter(all_births > 0)

# Scale PM2.5 by 10
truedta_null <- truedta_null %>% mutate(pm25 = pm25 / 10)
set.seed(212)
# Calculate lambda based on observed rates
lambda <- mean(truedta_null$preterm_births / truedta_null$all_births, na.rm = TRUE)

# Define beta1 (log RR for a 10-unit increase in PM2.5)
beta1 <- log(1.12)  # RR = 1.12

simulated_counts <- replicate(100, {
  rpois(nrow(truedta_null), lambda * (truedta_null$all_births * exp(beta1 * truedta_null$pm25)))
})

# Compute mean simulated counts and round to integer
final_outc <- round(rowMeans(simulated_counts))

# Ensure preterm births do not exceed total births
final_outc <- pmin(final_outc, truedta_null$preterm_births)

# Create final dataset
truedta_ass <- data.frame(truedta_null[, c("ZIP", "date", "all_births", "pm25")], preterm_births = final_outc)

# Fit model
mod2 <- glm(preterm_births ~ pm25 + offset(log(all_births)), data = truedta_ass, family = quasipoisson)
summary(mod2)
exp(coef(mod2)["pm25"])  # Final estimated RR
write_csv(truedta_ass, "truedta_ass_rpois.csv")
