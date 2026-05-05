################################################################################
# sims1.R  —  Scenario 1: Non-differential error, NULL association
#
# Purpose : Run 200 simulations adding classical (non-differential) measurement
#           error to the null dataset.  Error is independent of PTB outcome,
#           so no bias is expected (true RR = 1.00).
#
# Error structure (Section 2.2, E1):
#   Z_it = X_it + v_it,  v_it ~ N(0, (err.lev * SD(X))^2)
#   Error is the same for all births regardless of PTB status => non-differential
#
# Usage   : Rscript sims1.R          (standalone screen session)
#           source("sims1.R")        (called from Sims_all.R — data pre-loaded)
#
# Output  : exp/exp_003/Results/nondiff_null.csv
#           exp/exp_003/Results/cor_nondiff_null.csv
################################################################################

##################### Vijay's simulations ##############################

rm(list=ls())

setwd("~/Downloads/CodeReview_Alan")

library(readr)
library(MASS)
library(dplyr)

##################### Logging ##########################################

log_file <- "./exp/exp_003/sims1_log.txt"

start_time <- Sys.time()

cat("Script: sims1.R (exp_003)\n", file = log_file)
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
##################### NON-DIFFERENTIAL ERROR Sims1.R ###################
########################################################################
### non diff null

nondiff.null <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.nondiff.null <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

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
    
    # error-prone exposure
    pm.err <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), -0.25,
            err.lev[e] * sd(truedta_null$pm25))
    
    truedta_null$pm_err <- pm.err
    
    # Quasi-Poisson GLMM with random intercept by ZIP (Equation 1, Section 2.3)
    # Random intercept accounts for within-ZIP clustering over time.
    # Quasi-Poisson handles residual overdispersion in daily PTB counts.
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
      
      nondiff.null[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
      nondiff.null[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
      
    } else {
      
      nondiff.null[s, (e*2 - 1)] <- NA
      nondiff.null[s, (e*2)]     <- NA
      
      fail_msg <- paste(
        "Model failed at simulation", s,
        "error level", err.lev[e], "\n"
      )
      
      cat(fail_msg)
      cat(fail_msg, file = log_file, append = TRUE)
      
    }
    
    cor.nondiff.null[s, e] <- cor(pm.err, truedta_null$pm25)
    
    rm(mod, pm.err)
    
  }
  
  done_msg <- paste("Simulation", s, "completed at", Sys.time(), "\n")
  
  cat(done_msg)
  cat(done_msg, file = log_file, append = TRUE)
  
}

##################### Save results ###################################

nondiff.null <- as.data.frame(nondiff.null)

names(nondiff.null) <- paste0(
  rep(c("b_", "se_"), times = length(err.lev)),
  rep(err.lev, each = 2)
)

cor.nondiff.null <- as.data.frame(cor.nondiff.null)

names(cor.nondiff.null) <- paste0("corr", err.lev)

write_csv(nondiff.null,
          "./exp/exp_003/Results/nondiff_null.csv")

write_csv(cor.nondiff.null,
          "./exp/exp_003/Results/cor_nondiff_null.csv")

##################### End logging ####################################

end_time <- Sys.time()
runtime  <- end_time - start_time

cat("\nEnd time:", end_time, "\n", file = log_file, append = TRUE)
cat("Total runtime:", runtime, "\n", file = log_file, append = TRUE)

cat("\nScript finished.\nTotal runtime:", runtime, "\n")

rm(nondiff.null, cor.nondiff.null)