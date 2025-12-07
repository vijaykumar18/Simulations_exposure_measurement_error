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

##################### NON-DIFFERENTIAL ERROR Sims1.R #####################

##################################################################
##################################################################
##################################################################
# non diff null
nondiff.null <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.nondiff.null <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){
  set.seed(212 + s)
  
  for (e in 1:length(err.lev)){
    
    # error-prone exposure
    pm.err <- truedta_null$pm25 + 
      rnorm(nrow(truedta_null), 0, err.lev[e] * sd(truedta_null$pm25))
    
    truedta_null$pm_err <- pm.err   # add to dataset
    
    # -------- RANDOM INTERCEPT MODEL (quasi-Poisson) --------
    mod <- glmmPQL(
      fixed  = preterm_births ~ pm_err + offset(log(all_births)),
      random = ~ 1 | ZIP,
      family = quasipoisson(link = "log"),
      data   = truedta_null,
      verbose = FALSE
    )
    # --------------------------------------------------------
    
    # coefficient and standard error (same position as glm)
    nondiff.null[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]  # beta
    nondiff.null[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]  # SE
    
    # correlation
    cor.nondiff.null[s, e] <- cor(pm.err, truedta_null$pm25)
    
    rm(mod, pm.err)
  }
}

# Same naming scheme as your glm version
nondiff.null <- as.data.frame(nondiff.null)
names(nondiff.null) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                              rep(err.lev, each = 2))

cor.nondiff.null <- as.data.frame(cor.nondiff.null)
names(cor.nondiff.null) <- paste0("corr", err.lev)

# Save results
write_csv(nondiff.null, "./exp/exp_001/Results/nondiff_null.csv")
write_csv(cor.nondiff.null, "./exp/exp_001/Results/cor_nondiff_null.csv")

rm(nondiff.null, cor.nondiff.null)