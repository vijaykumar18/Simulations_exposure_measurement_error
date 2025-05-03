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

##################### NON-DIFFERENTIAL ERROR #####################

##################################################################
##################################################################
##################################################################

#
# ### null association
#
nondiff.null<- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)  # Store b and se
cor.nondiff.null <- matrix(NA, nrow = n.sims, ncol = length(err.lev))
#
for (s in 1:n.sims){
  index <- 1
  for (e in 1:length(err.lev)){ # e=1
    pm.err <- rnorm(nrow(truedta_null), truedta_null$pm25-0.25, err.lev[e]*sd(truedta_null$pm25))
    mod    <- glm(preterm_births ~ pm.err + offset(log(all_births)),
                  data = truedta_null, family=quasipoisson)
    nondiff.null[s, c(index, (index+1))] <- summary(mod)$coefficients[2,1:2]
    cor.nondiff.null[s,e] <- cor(pm.err, truedta_null$pm25)
    rm(mod, pm.err)
    index <- index + 2
  }}

nondiff.null        <- as.data.frame(nondiff.null)
names(nondiff.null) <- paste0(rep(c("b_", "se_"), each = length(err.lev)), err.lev)
write_csv(nondiff.null, "./exp/exp_003/Results/nondiff_null.csv")

cor.nondiff.null        <- as.data.frame(cor.nondiff.null)
names(cor.nondiff.null) <- paste0("corr", err.lev)
write_csv(cor.nondiff.null, "./exp/exp_003/Results/cor_nondiff_null.csv")

rm(nondiff.null, cor.nondiff.null)
