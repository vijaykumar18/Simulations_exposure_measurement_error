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
      rnorm(nrow(truedta_null), 0.25, err.lev[e] * sd(truedta_null$pm25))

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
write_csv(nondiff.null, "./exp/exp_002/Results/nondiff_null.csv")
write_csv(cor.nondiff.null, "./exp/exp_002/Results/cor_nondiff_null.csv")

rm(nondiff.null, cor.nondiff.null)

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
      rnorm(nrow(truedta_ass), 0.25, err.lev[e] * sd(truedta_ass$pm25))

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
write_csv(nondiff.harm, "./exp/exp_002/Results/nondiff_harm.csv")
write_csv(cor.nondiff.harm, "./exp/exp_002/Results/cor_nondiff_harm.csv")


rm(nondiff.harm, cor.nondiff.harm)



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




##################################################################
##################################################################
##################################################################

####################### DIFFERENTIAL ERROR Sims4.R #######################

##################################################################
##################################################################
##################################################################

# ## more error among those high-risk
diff.null.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.diff.null.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){
  set.seed(212 + s)

  for (e in 1:length(err.lev)){

    # Generate exposure error separately by risk group
    pm.err.norisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, err.lev[e] * sd(truedta_null$pm25))

    pm.err.hirisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0.25, (err.lev[e] * 10) * sd(truedta_null$pm25))

    # Weighted combination
    pm.err <- ifelse(
      truedta_null$all_births > 0,
      (truedta_null$term_births    * pm.err.norisk +
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

    # beta and se
    diff.null.moreErr[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
    diff.null.moreErr[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]

    # correlation
    cor.diff.null.moreErr[s, e] <- cor(pm.err, truedta_null$pm25)

    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
}

# Naming (same as glm version)
diff.null.moreErr <- as.data.frame(diff.null.moreErr)
names(diff.null.moreErr) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                                   rep(err.lev, each = 2))

cor.diff.null.moreErr <- as.data.frame(cor.diff.null.moreErr)
names(cor.diff.null.moreErr) <- paste0("corr", err.lev)

# Save files
write_csv(diff.null.moreErr, "./exp/exp_002/Results/diff_null_moreErr.csv")
write_csv(cor.diff.null.moreErr, "./exp/exp_002/Results/cor_diff_null_moreErr.csv")

rm(diff.null.moreErr, cor.diff.null.moreErr)

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


##################################################################
##################################################################
##################################################################

####################### DIFFERENTIAL ERROR Sims6.R #######################

##################################################################
##################################################################
##################################################################

## more error among those high-risk & 50% of the rest

diff.null.moreErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.diff.null.moreErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){
  set.seed(212 + s)

  for (e in 1:length(err.lev)){

    # Generate exposure error
    pm.err.norisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, err.lev[e] * sd(truedta_null$pm25))

    pm.err.hirisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0.25, (err.lev[e] * 10) * sd(truedta_null$pm25))

    # Weighted combination (50/50 mixture for term births, hi-risk gets 10× error)
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

    # Extract coefficient and SE
    diff.null.moreErr50[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
    diff.null.moreErr50[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]

    # Correlation
    cor.diff.null.moreErr50[s, e] <- cor(pm.err, truedta_null$pm25)

    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
}

# Name columns exactly as original
diff.null.moreErr50 <- as.data.frame(diff.null.moreErr50)
names(diff.null.moreErr50) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                                     rep(err.lev, each = 2))

cor.diff.null.moreErr50 <- as.data.frame(cor.diff.null.moreErr50)
names(cor.diff.null.moreErr50) <- paste0("corr", err.lev)

# Save
write_csv(diff.null.moreErr50, "./exp/exp_002/Results/diff_null_moreErr50.csv")
write_csv(cor.diff.null.moreErr50, "./exp/exp_002/Results/cor_diff_null_moreErr50.csv")

rm(diff.null.moreErr50, cor.diff.null.moreErr50)


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

##################################################################
##################################################################
##################################################################

####################### DIFFERENTIAL ERROR Sims8.R #######################

##################################################################
##################################################################
##################################################################

## more error among those high-risk

diff.harm.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.harm.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){ 
  set.seed(212 + s)

  for (e in 1:length(err.lev)){

    # Generate differential *more error* for high-risk
    pm.err.norisk <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), 0, err.lev[e] * sd(truedta_ass$pm25))

    pm.err.hirisk <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), 0.25, (err.lev[e] * 10) * sd(truedta_ass$pm25))

    # Weighted PM by birth counts
    pm.err <- ifelse(
      truedta_ass$all_births > 0,
      (truedta_ass$term_births    * pm.err.norisk +
       truedta_ass$preterm_births * pm.err.hirisk) /
         truedta_ass$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )

    # Add to dataset
    truedta_ass$pm_err <- pm.err

    # -------- RANDOM INTERCEPT QUASI-POISSON ----------
    mod <- glmmPQL(
      fixed  = preterm_births ~ pm_err + offset(log(all_births)),
      random = ~ 1 | ZIP,
      family = quasipoisson(link = "log"),
      data   = truedta_ass,
      verbose = FALSE
    )
    # ---------------------------------------------------

    # Store coefficient and SE
    diff.harm.moreErr[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
    diff.harm.moreErr[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]

    # Correlation
    cor.diff.harm.moreErr[s, e] <- cor(pm.err, truedta_ass$pm25)

    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
}

# Column names
diff.harm.moreErr <- as.data.frame(diff.harm.moreErr)
names(diff.harm.moreErr) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                                   rep(err.lev, each = 2))

cor.diff.harm.moreErr <- as.data.frame(cor.diff.harm.moreErr)
names(cor.diff.harm.moreErr) <- paste0("corr", err.lev)

# Save results
write_csv(diff.harm.moreErr, "./exp/exp_002/Results/diff_harm_moreErr.csv")
write_csv(cor.diff.harm.moreErr, "./exp/exp_002/Results/cor_diff_harm_moreErr.csv")
rm(diff.harm.moreErr, cor.diff.harm.moreErr)


##################################################################
##################################################################
##################################################################

####################### DIFFERENTIAL ERROR Sims9.R #######################

##################################################################
##################################################################
##################################################################

## less error among those high-risk & 50% of the rest

diff.harm.lessErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.harm.lessErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){ 
  set.seed(212 + s)
  
  for (e in 1:length(err.lev)){

    # Generate differential errors
    pm.err.norisk <- truedta_ass$pm25 + 
      rnorm(nrow(truedta_ass), 0, err.lev[e] * sd(truedta_ass$pm25))
    
    pm.err.hirisk <- truedta_ass$pm25 + 
      rnorm(nrow(truedta_ass), 0.25, (err.lev[e] / 10) * sd(truedta_ass$pm25))
    
    # Weighted combination by birth counts
    pm.err <- ifelse(
      truedta_ass$all_births > 0,
      (0.5 * truedta_ass$term_births * pm.err.norisk +
       0.5 * truedta_ass$term_births * pm.err.hirisk +
       truedta_ass$preterm_births * pm.err.hirisk) / truedta_ass$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )

    # Add to dataset
    truedta_ass$pm_err <- pm.err

    # -------- RANDOM INTERCEPT QUASI-POISSON ----------
    mod <- glmmPQL(
      fixed  = preterm_births ~ pm_err + offset(log(all_births)),
      random = ~ 1 | ZIP,
      family = quasipoisson(link = "log"),
      data   = truedta_ass,
      verbose = FALSE
    )
    # ---------------------------------------------------

    # Store coefficient and SE
    diff.harm.lessErr50[s, (e*2-1)] <- summary(mod)$tTable["pm_err", 1]
    diff.harm.lessErr50[s, (e*2)]   <- summary(mod)$tTable["pm_err", 2]

    # Store correlation
    cor.diff.harm.lessErr50[s, e] <- cor(pm.err, truedta_ass$pm25)

    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
}

# Column names
diff.harm.lessErr50 <- as.data.frame(diff.harm.lessErr50)
names(diff.harm.lessErr50) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                                    rep(err.lev, each = 2))

cor.diff.harm.lessErr50 <- as.data.frame(cor.diff.harm.lessErr50)
names(cor.diff.harm.lessErr50) <- paste0("corr", err.lev)

# Save results
write_csv(diff.harm.lessErr50, "./exp/exp_002/Results/diff_harm_lessErr50.csv")
write_csv(cor.diff.harm.lessErr50, "./exp/exp_002/Results/cor_diff_harm_lessErr50.csv")
rm(diff.harm.lessErr50, cor.diff.harm.lessErr50)


#################################################################
##################################################################
##################################################################

####################### DIFFERENTIAL ERROR Sims10.R #######################

##################################################################
##################################################################
##################################################################

## more error among those high-risk & 50% of the rest

diff.harm.moreErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.harm.moreErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){ 
  set.seed(212 + s)
  
  for (e in 1:length(err.lev)){

    # Generate differential errors
    pm.err.norisk <- truedta_ass$pm25 + 
      rnorm(nrow(truedta_ass), 0, err.lev[e] * sd(truedta_ass$pm25))
    
    pm.err.hirisk <- truedta_ass$pm25 + 
      rnorm(nrow(truedta_ass), 0.25, (err.lev[e] * 10) * sd(truedta_ass$pm25))
    
    # Weighted combination by birth counts
    pm.err <- ifelse(
      truedta_ass$all_births > 0,
      (0.5 * truedta_ass$term_births * pm.err.norisk +
       0.5 * truedta_ass$term_births * pm.err.hirisk +
       truedta_ass$preterm_births * pm.err.hirisk) / truedta_ass$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )

    truedta_ass$pm_err <- pm.err

    # -------- RANDOM INTERCEPT QUASI-POISSON ----------
    mod <- glmmPQL(
      fixed  = preterm_births ~ pm_err + offset(log(all_births)),
      random = ~ 1 | ZIP,
      family = quasipoisson(link = "log"),
      data   = truedta_ass,
      verbose = FALSE
    )
    # ---------------------------------------------------

    # Store coefficient and SE
    diff.harm.moreErr50[s, (e*2-1)] <- summary(mod)$tTable["pm_err", 1]
    diff.harm.moreErr50[s, (e*2)]   <- summary(mod)$tTable["pm_err", 2]

    # Store correlation
    cor.diff.harm.moreErr50[s, e] <- cor(pm.err, truedta_ass$pm25)

    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
}

# Column names
diff.harm.moreErr50 <- as.data.frame(diff.harm.moreErr50)
names(diff.harm.moreErr50) <- paste0(rep(c("b_", "se_"), times = length(err.lev)),
                                     rep(err.lev, each = 2))

cor.diff.harm.moreErr50 <- as.data.frame(cor.diff.harm.moreErr50)
names(cor.diff.harm.moreErr50) <- paste0("corr", err.lev)

# Save results
write_csv(diff.harm.moreErr50, "./exp/exp_002/Results/diff_harm_moreErr50.csv")
write_csv(cor.diff.harm.moreErr50, "./exp/exp_002/Results/cor_diff_harm_moreErr50.csv")

rm(diff.harm.moreErr50, cor.diff.harm.moreErr50)
