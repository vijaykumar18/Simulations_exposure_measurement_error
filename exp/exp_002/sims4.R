################################################################################
# sims4.R  —  Scenario 4: Differential error (MORE error in at-risk), NULL
#
# Purpose : Run 200 simulations where the at-risk (preterm) population has
#           MORE measurement error than the non-at-risk population.
#           Applied to the null dataset (true RR = 1.00).
#
# Error structure (Section 2.2, Table 1, E1):
#   Non-at-risk (term births):  v ~ N(0, (err.lev      * SD(X))^2)
#   At-risk (preterm births):   v ~ N(0, (err.lev * 10 * SD(X))^2)
#
# Aggregation: Birth-count-weighted ZIP-day exposure (see sims3.R for details).
#
# Usage   : Rscript sims4.R
#           source("sims4.R")  (from Sims_all.R)
#
# Output  : exp/exp_002/Results/diff_null_moreErr.csv
#           exp/exp_002/Results/cor_diff_null_moreErr.csv
################################################################################

##################### Vijay's simulations ##############################

rm(list=ls())

setwd("~/Downloads/CodeReview_Alan")

library(readr)
library(MASS)
library(dplyr)
library(lme4)

##################### Logging ##########################################

log_file <- "./exp/exp_002/sims4_log.txt"

start_time <- Sys.time()

cat("Script: sims4.R (exp_002)\n", file = log_file)
cat("Start time:", start_time, "\n\n", file = log_file, append = TRUE)

########################################################################
## 1. bring in "true data"

truedta_null <- read_csv("truedta_null_Final.csv")
truedta_ass  <- read_csv("truedta_ass_Final.csv")

truedta_null$term_births <- truedta_null$all_births - truedta_null$preterm_births
truedta_ass$term_births  <- truedta_ass$all_births - truedta_ass$preterm_births

################### PTB scenario ######################################

set.seed(212)
n.sims  <- 5
err.lev <- c(0, 0.5, 1, 2)

########################################################################
################ DIFFERENTIAL ERROR Sims4.R ###########################
########################################################################

diff.null.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.null.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

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
    
    #  Non-at-risk: baseline error
    
    pm.err.norisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0,
            err.lev[e] * sd(truedta_null$pm25))
    
    # At-risk: 10x more error (err.lev * 10)
    # Represents high-risk individuals having less accurate exposure estimates
    # (e.g., residing in areas with sparser monitoring)
    pm.err.hirisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0.25,
            (err.lev[e] * 10) * sd(truedta_null$pm25))
    
    # Weighted exposure
    
    pm.err <- ifelse(
      truedta_null$all_births > 0,
      (truedta_null$term_births * pm.err.norisk +
         truedta_null$preterm_births * pm.err.hirisk) /
        truedta_null$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )
    
    truedta_null$pm_err <- pm.err
    
    mod <- try(
      glmmPQL(
        fixed  = preterm_births ~ pm_err + offset(log(all_births)),
        random = ~1 | ZIP,
        family = quasipoisson(link = "log"),
        data   = truedta_null,
        verbose = FALSE
      ),
      silent = TRUE
    )
    
    if(!inherits(mod, "try-error")){
      
      diff.null.moreErr[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
      diff.null.moreErr[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
      
    } else {
      
      diff.null.moreErr[s, (e*2 - 1)] <- NA
      diff.null.moreErr[s, (e*2)]     <- NA
      
      fail_msg <- paste(
        "Model failed at simulation", s,
        "error level", err.lev[e], "\n"
      )
      
      cat(fail_msg)
      cat(fail_msg, file = log_file, append = TRUE)
      
    }
    
    cor.diff.null.moreErr[s, e] <- cor(pm.err, truedta_null$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
    
  }
  
  done_msg <- paste("Simulation", s, "completed at", Sys.time(), "\n")
  
  cat(done_msg)
  cat(done_msg, file = log_file, append = TRUE)
  
}

##################### Save results ###################################

diff.null.moreErr <- as.data.frame(diff.null.moreErr)

names(diff.null.moreErr) <- paste0(
  rep(c("b_", "se_"), times = length(err.lev)),
  rep(err.lev, each = 2)
)

cor.diff.null.moreErr <- as.data.frame(cor.diff.null.moreErr)

names(cor.diff.null.moreErr) <- paste0("corr", err.lev)

write_csv(diff.null.moreErr,
          "./exp/exp_002/Results/diff_null_moreErr.csv")

write_csv(cor.diff.null.moreErr,
          "./exp/exp_002/Results/cor_diff_null_moreErr.csv")

##################### End logging ####################################

end_time <- Sys.time()
runtime  <- end_time - start_time

cat("\nEnd time:", end_time, "\n", file = log_file, append = TRUE)
cat("Total runtime:", runtime, "\n", file = log_file, append = TRUE)

cat("\nScript finished.\nTotal runtime:", runtime, "\n")

rm(diff.null.moreErr, cor.diff.null.moreErr)