#####################Vijay's simulations##############
setwd("~/Downloads/DiffMeasError/CodeReview")

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

## more error among those high-risk

diff.harm.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.diff.harm.moreErr <- matrix(NA, nrow = n.sims, ncol=length(err.lev))

for (s in 1:n.sims){ 
  index <- 1
  for (e in 1:length(err.lev)){ # e=1
    pm.err.norisk <- rnorm(nrow(truedta_ass), truedta_ass$pm25, err.lev[e]*sd(truedta_ass$pm25))
    pm.err.hirisk <- rnorm(nrow(truedta_ass), truedta_ass$pm25-0.25, (err.lev[e]*10)*sd(truedta_ass$pm25))
    
    pm.err <- ifelse(truedta_ass$all_births > 0,
                     (truedta_ass$term_births*pm.err.norisk + truedta_ass$preterm_births*pm.err.hirisk)/(truedta_ass$all_births),
                     mean(c(pm.err.norisk, pm.err.hirisk)))
    
    mod <- glm(preterm_births ~ pm.err + offset(log(all_births)), 
               data = truedta_ass, family=quasipoisson)
    
    diff.harm.moreErr[s, c(index, (index+1))] <- summary(mod)$coefficients[2,1:2]
    cor.diff.harm.moreErr[s,e] <- cor(pm.err, truedta_ass$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
    index <- index + 2
  }}

diff.harm.moreErr        <- as.data.frame(diff.harm.moreErr)
names(diff.harm.moreErr) <- paste0(rep(c("b_", "se_"), each = length(err.lev)), err.lev)
write_csv(diff.harm.moreErr, "./exp/exp_003/Results/diff_harm_moreErr.csv")

cor.diff.harm.moreErr        <- as.data.frame(cor.diff.harm.moreErr)
names(cor.diff.harm.moreErr) <- paste0("corr", err.lev)
write_csv(cor.diff.harm.moreErr, "./exp/exp_003/Results/cor_diff_harm_moreErr.csv")

rm(diff.harm.moreErr, cor.diff.harm.moreErr)
