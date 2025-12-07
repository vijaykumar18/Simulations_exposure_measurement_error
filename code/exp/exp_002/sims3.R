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

####################### DIFFERENTIAL ERROR Sims3.R #######################

##################################################################
##################################################################
##################################################################

## less error among those high-risk


diff.null.lessErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.diff.null.lessErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){
  set.seed(212 + s)
  
  for (e in 1:length(err.lev)){
    
    # Generate exposure error separately by risk group
    pm.err.norisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, err.lev[e] * sd(truedta_null$pm25))
    
    pm.err.hirisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0.25, (err.lev[e] / 10) * sd(truedta_null$pm25))
    
    # Weighted combination
    pm.err <- ifelse(
      truedta_null$all_births > 0,
      (truedta_null$term_births    * pm.err.norisk +
         truedta_null$preterm_births * pm.err.hirisk) /
        truedta_null$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )
    
    # Add to data
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
    
    # beta and se (same indexing as glm version)
    diff.null.lessErr[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
    diff.null.lessErr[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
    
    # correlation
    cor.diff.null.lessErr[s, e] <- cor(pm.err, truedta_null$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
}

# Same naming format as before
diff.null.lessErr <- as.data.frame(diff.null.lessErr)
names(diff.null.lessErr) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                                   rep(err.lev, each = 2))

cor.diff.null.lessErr <- as.data.frame(cor.diff.null.lessErr)
names(cor.diff.null.lessErr) <- paste0("corr", err.lev)

# Save
write_csv(diff.null.lessErr, "./exp/exp_002/Results/diff_null_lessErr.csv")
write_csv(cor.diff.null.lessErr, "./exp/exp_002/Results/cor_diff_null_lessErr.csv")

rm(diff.null.lessErr, cor.diff.null.lessErr)