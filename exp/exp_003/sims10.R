################################################################################
# sims10.R  —  Scenario 10: Differential error (MORE error in at-risk + 50% non-risk), HARMFUL
#
# Purpose : Same structure as sims6.R (more error + 50% split) but applied to
#           the harmful association dataset (true RR = 1.115).
#           This is the most extreme differential error scenario and produces
#           the largest attenuation bias (up to -99%).
#
# Error structure (Section 2.2, Table 1, E1):
#   At-risk (preterm births):      v ~ N(0, (err.lev * 10 * SD(X))^2)
#   50% of non-at-risk (term):     v ~ N(0, (err.lev * 10 * SD(X))^2)
#   50% of non-at-risk (term):     v ~ N(0, (err.lev       * SD(X))^2)
#
# Usage   : Rscript sims10.R  /  source("sims10.R")
# Output  : exp/exp_003/Results/diff_harm_moreErr50.csv
#           exp/exp_003/Results/cor_diff_harm_moreErr50.csv
################################################################################

##################### Vijay's simulations ##############################

rm(list=ls())

setwd("~/Downloads/CodeReview_Alan")

library(readr)
library(MASS)
library(dplyr)

##################### Logging ##########################################

log_file <- "./exp/exp_003/sims10_log.txt"

start_time <- Sys.time()

cat("Script: sims10.R (exp_003)\n", file = log_file)
cat("Start time:", start_time, "\n\n", file = log_file, append = TRUE)

########################################################################
## 1. bring in "true data"

truedta_null <- read_csv("truedta_null_Final.csv")
truedta_ass  <- read_csv("truedta_ass_Final.csv")

truedta_null$term_births <- truedta_null$all_births - truedta_null$preterm_births
truedta_ass$term_births  <- truedta_ass$all_births - truedta_ass$preterm_births

################### PTB scenario ######################################

set.seed(212)
n.sims <- 5
err.lev <- c(0, 0.5, 1, 2)

########################################################################
####################### DIFFERENTIAL ERROR Sims10.R ####################
########################################################################

## more error among those high-risk & 50% of the rest

diff.harm.moreErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.harm.moreErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

for (s in 1:n.sims){
  
  progress <- round(100 * s / n.sims, 1)
  
  msg <- paste0(
    "Simulation ", s, "/", n.sims,
    " (", progress, "%) started at ", Sys.time(), "\n"
  )
  
  cat(msg)
  cat(msg, file = log_file, append = TRUE)
  
  set.seed(212 + s)
  
  for (e in 1:length(err.lev)){
    
    # Generate differential errors
    pm.err.norisk <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), 0,
            err.lev[e] * sd(truedta_ass$pm25))
    
    pm.err.hirisk <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), -0.25,
            (err.lev[e] * 10) * sd(truedta_ass$pm25))
    
    # Weighted exposure
    pm.err <- ifelse(
      truedta_ass$all_births > 0,
      (0.5 * truedta_ass$term_births * pm.err.norisk +
         0.5 * truedta_ass$term_births * pm.err.hirisk +
         truedta_ass$preterm_births * pm.err.hirisk) /
        truedta_ass$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )
    
    truedta_ass$pm_err <- pm.err
    
    mod <- try(
      glmmPQL(
        fixed  = preterm_births ~ pm_err + offset(log(all_births)),
        random = ~1 | ZIP,
        family = quasipoisson(link = "log"),
        data   = truedta_ass,
        verbose = FALSE
      ),
      silent = TRUE
    )
    
    if(!inherits(mod, "try-error")){
      
      diff.harm.moreErr50[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
      diff.harm.moreErr50[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
      
    } else {
      
      diff.harm.moreErr50[s, (e*2 - 1)] <- NA
      diff.harm.moreErr50[s, (e*2)]     <- NA
      
      fail_msg <- paste(
        "Model failed at simulation", s,
        "error level", err.lev[e], "\n"
      )
      
      cat(fail_msg)
      cat(fail_msg, file = log_file, append = TRUE)
      
    }
    
    cor.diff.harm.moreErr50[s, e] <- cor(pm.err, truedta_ass$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
    
  }
  
  done_msg <- paste("Simulation", s, "completed at", Sys.time(), "\n")
  
  cat(done_msg)
  cat(done_msg, file = log_file, append = TRUE)
  
}

##################### Save results ###################################

diff.harm.moreErr50 <- as.data.frame(diff.harm.moreErr50)

names(diff.harm.moreErr50) <- paste0(
  rep(c("b_", "se_"), times = length(err.lev)),
  rep(err.lev, each = 2)
)

cor.diff.harm.moreErr50 <- as.data.frame(cor.diff.harm.moreErr50)

names(cor.diff.harm.moreErr50) <- paste0("corr", err.lev)

write_csv(diff.harm.moreErr50,
          "./exp/exp_003/Results/diff_harm_moreErr50.csv")

write_csv(cor.diff.harm.moreErr50,
          "./exp/exp_003/Results/cor_diff_harm_moreErr50.csv")

##################### End logging ####################################

end_time <- Sys.time()
runtime <- end_time - start_time

cat("\nEnd time:", end_time, "\n", file = log_file, append = TRUE)
cat("Total runtime:", runtime, "\n", file = log_file, append = TRUE)

cat("\nScript finished.\nTotal runtime:", runtime, "\n")

rm(diff.harm.moreErr50, cor.diff.harm.moreErr50)