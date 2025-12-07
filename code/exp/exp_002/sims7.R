################new code for running PTB simulations############

#####################Vijay's simulations##############
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
library(parallel)
library(data.table)

rm(list=ls())

## 1. bring in "true data"

truedta_null <- read_csv("truedta_null_Final.csv")

truedta_ass <- read_csv("truedta_ass_Final.csv")


truedta_null$term_births <- truedta_null$all_births - truedta_null$preterm_births
truedta_ass$term_births  <- truedta_ass$all_births - truedta_ass$preterm_births

###################PTB scenario############################

set.seed(212)
n.sims <- 200
err.lev <- c(0, 0.5, 1, 2)

##################################################################
##################################################################
##################################################################

####################### DIFFERENTIAL ERROR Sims7.R #######################

##################################################################
##################################################################
##################################################################
### harmful association

# less error among those high-risk

diff.harm.lessErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.diff.harm.lessErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims) {
  set.seed(212 + s)
  
  for (e in 1:length(err.lev)) {
    
    # Generate differential errors (less error for high-risk)
    pm.err.norisk <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), 0, err.lev[e] * sd(truedta_ass$pm25))
    
    pm.err.hirisk <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), 0.25, (err.lev[e] / 10) * sd(truedta_ass$pm25))
    
    # Weighted exposure by outcome prevalence
    pm.err <- ifelse(
      truedta_ass$all_births > 0,
      (truedta_ass$term_births    * pm.err.norisk +
         truedta_ass$preterm_births * pm.err.hirisk) /
        truedta_ass$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )
    
    # Add variable to dataset
    truedta_ass$pm_err <- pm.err
    
    # -------- RANDOM INTERCEPT quasi-Poisson model --------
    mod <- glmmPQL(
      fixed  = preterm_births ~ pm_err + offset(log(all_births)),
      random = ~ 1 | ZIP,
      family = quasipoisson(link = "log"),
      data   = truedta_ass,
      verbose = FALSE
    )
    # -------------------------------------------------------
    
    # Store coefficient and standard error
    diff.harm.lessErr[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
    diff.harm.lessErr[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
    
    # Correlation
    cor.diff.harm.lessErr[s, e] <- cor(pm.err, truedta_ass$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
}

# Proper column names
diff.harm.lessErr <- as.data.frame(diff.harm.lessErr)
names(diff.harm.lessErr) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                                   rep(err.lev, each = 2))

cor.diff.harm.lessErr <- as.data.frame(cor.diff.harm.lessErr)
names(cor.diff.harm.lessErr) <- paste0("corr", err.lev)

# Save files
write_csv(diff.harm.lessErr, "./exp/exp_002/Results/diff_harm_lessErr.csv")
write_csv(cor.diff.harm.lessErr, "./exp/exp_002/Results/cor_diff_harm_lessErr.csv")

rm(diff.harm.lessErr, cor.diff.harm.lessErr)