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

####################### DIFFERENTIAL ERROR Sims5.R #######################

##################################################################
##################################################################
##################################################################

# ## less error among those high-risk & 50% of the rest

diff.null.lessErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.diff.null.lessErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){
  set.seed(212 + s)
  
  for (e in 1:length(err.lev)){
    
    # Generate exposure error for different risk groups
    pm.err.norisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, err.lev[e] * sd(truedta_null$pm25))
    
    pm.err.hirisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0.25, (err.lev[e] / 10) * sd(truedta_null$pm25))
    
    # Weighted combination (50/50 mix for term births)
    pm.err <- ifelse(
      truedta_null$all_births > 0,
      (0.5 * truedta_null$term_births    * pm.err.norisk +
         0.5 * truedta_null$term_births    * pm.err.hirisk +
         truedta_null$preterm_births * pm.err.hirisk) /
        truedta_null$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )
    
    # Add to dataset
    truedta_null$pm_err <- pm.err
    
    # -------- RANDOM INTERCEPT quasi-Poisson model --------
    mod <- glmmPQL(
      fixed  = preterm_births ~ pm_err + offset(log(all_births)),
      random = ~ 1 | ZIP,
      family = quasipoisson(link = "log"),
      data   = truedta_null,
      verbose = FALSE
    )
    # -------------------------------------------------------
    
    # Extract beta and SE
    diff.null.lessErr50[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
    diff.null.lessErr50[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
    
    # Correlation
    cor.diff.null.lessErr50[s, e] <- cor(pm.err, truedta_null$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
}

# Column names
diff.null.lessErr50 <- as.data.frame(diff.null.lessErr50)
names(diff.null.lessErr50) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                                     rep(err.lev, each = 2))

cor.diff.null.lessErr50 <- as.data.frame(cor.diff.null.lessErr50)
names(cor.diff.null.lessErr50) <- paste0("corr", err.lev)

# Save output
write_csv(diff.null.lessErr50, "./exp/exp_002/Results/diff_null_lessErr50.csv")
write_csv(cor.diff.null.lessErr50, "./exp/exp_002/Results/cor_diff_null_lessErr50.csv")

rm(diff.null.lessErr50, cor.diff.null.lessErr50)