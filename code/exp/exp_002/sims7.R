#####################Vijay's simulations##############
setwd("~/Downloads/DiffMeasError")

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

truedta_null <- read_csv("truedta_null_rpois.csv")

truedta_ass <- read_csv("truedta_ass_rpois.csv")

# mod2 <- glm(preterm_births ~ pm25 + offset(loglog(all_births)), 
#             data = truedta_null, family=quasipoisson)
# summary(mod2)


truedta_null$term_births <- truedta_null$all_births - truedta_null$preterm_births
truedta_ass$term_births  <- truedta_ass$all_births - truedta_ass$preterm_births

set.seed(212)
n.sims <- 200
err.lev <- c(0, 0.5, 1, 2)



##################################################################
##################################################################
##################################################################

####################### DIFFERENTIAL ERROR #######################

##################################################################
##################################################################
##################################################################
### harmful association

## less error among those high-risk

diff.harm.lessErr     <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.diff.harm.lessErr <- matrix(NA, nrow = n.sims, ncol=length(err.lev))

for (s in 1:n.sims) {  # Iterate over the number of simulations
  set.seed(212 + s)    # Set a unique seed for each simulation
  
  index <- 1
  for (e in 1:length(err.lev)) { # Iterate over the error levels
    # Generate errors for no-risk and high-risk groups
    pm.err.norisk <- rnorm(nrow(truedta_ass), truedta_ass$pm25, err.lev[e] * sd(truedta_ass$pm25))
    pm.err.hirisk <- rnorm(nrow(truedta_ass), truedta_ass$pm25+0.25, (err.lev[e] / 10) * sd(truedta_ass$pm25))
    
    # Combine errors using weighting scheme
    pm.err <- ifelse(truedta_ass$all_births > 0,
                     (truedta_ass$term_births * pm.err.norisk + truedta_ass$preterm_births * pm.err.hirisk) / 
                       (truedta_ass$all_births),
                     mean(c(pm.err.norisk, pm.err.hirisk)))
    
    # Fit the model
    mod <- glm(preterm_births ~ pm.err + offset(log(all_births)),
               data = truedta_ass, family = quasipoisson)
    
    # Extract results
    diff.harm.lessErr[s, c(index, (index + 1))] <- summary(mod)$coefficients[2, 1:2]
    cor.diff.harm.lessErr[s, e] <- cor(pm.err, truedta_ass$pm25)
    
    # Clean up and increment index
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
    index <- index + 2
  }
}


diff.harm.lessErr        <- as.data.frame(diff.harm.lessErr)
names(diff.harm.lessErr) <- paste0(rep(c("b_", "se_"), each = length(err.lev)), err.lev)
write_csv(diff.harm.lessErr, "./exp/exp_002/Results/diff_harm_lessErr.csv")

cor.diff.harm.lessErr        <- as.data.frame(cor.diff.harm.lessErr)
names(cor.diff.harm.lessErr) <- paste0("corr", err.lev)
write_csv(cor.diff.harm.lessErr, "./exp/exp_002/Results/cor_diff_harm_lessErr.csv")

rm(diff.harm.lessErr, cor.diff.harm.lessErr)
