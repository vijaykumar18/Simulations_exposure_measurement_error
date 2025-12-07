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

##################### NON-DIFFERENTIAL ERROR Sims2.R #####################

##################################################################
##################################################################
##################################################################

### harmful association
nondiff.harm <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.nondiff.harm <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){
  set.seed(212 + s)
  
  for (e in 1:length(err.lev)){
    
    # Mismeasured exposure
    pm.err <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), -0.25, err.lev[e] * sd(truedta_ass$pm25))
    
    # Add to dataset
    truedta_ass$pm_err <- pm.err
    
    # -------- RANDOM INTERCEPT QUASI-POISSON MODEL --------
    mod <- glmmPQL(
      fixed  = preterm_births ~ pm_err + offset(log(all_births)),
      random = ~ 1 | ZIP,
      family = quasipoisson(link = "log"),
      data   = truedta_ass,
      verbose = FALSE
    )
    # ------------------------------------------------------
    
    # Extract coefficient + SE (same indexing as glm version)
    nondiff.harm[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]  # beta
    nondiff.harm[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]  # se
    
    # Correlation
    cor.nondiff.harm[s, e] <- cor(pm.err, truedta_ass$pm25)
    
    rm(mod, pm.err)
  }
}

# Column names (same structure as before)
nondiff.harm <- as.data.frame(nondiff.harm)
names(nondiff.harm) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                              rep(err.lev, each = 2))

cor.nondiff.harm <- as.data.frame(cor.nondiff.harm)
names(cor.nondiff.harm) <- paste0("corr", err.lev)

# Save
write_csv(nondiff.harm, "./exp/exp_003/Results/nondiff_harm.csv")
write_csv(cor.nondiff.harm, "./exp/exp_003/Results/cor_nondiff_harm.csv")


rm(nondiff.harm, cor.nondiff.harm)