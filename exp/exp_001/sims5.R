################################################################################
# sims5.R  —  Scenario 5: Differential error (LESS error in at-risk + 50% non-risk), NULL
#
# Purpose : Run 200 simulations where the at-risk population AND 50% of the
#           non-at-risk population have LESS error; the remaining 50% of the
#           non-at-risk population have the baseline error.
#           Applied to the null dataset (true RR = 1.00).
#
# Error structure (Section 2.2, Table 1, E1):
#   At-risk (preterm births):      v ~ N(0, (err.lev / 10 * SD(X))^2)
#   50% of non-at-risk (term):     v ~ N(0, (err.lev / 10 * SD(X))^2)  [less]
#   50% of non-at-risk (term):     v ~ N(0, (err.lev       * SD(X))^2)  [baseline]
#
# Aggregation: The 50/50 split for term births is implemented by mixing the two
#   error processes with equal weight in the weighted average:
#     pm.err = (0.5*term*Z_norisk + 0.5*term*Z_hirisk + preterm*Z_hirisk) / all_births
#   This is equivalent to: 50% of term births get baseline error, 50% get less
#   error, and all preterm births get less error.
#
# Usage   : Rscript sims5.R
#           source("sims5.R")  (from Sims_all.R)
#
# Output  : exp/exp_001/Results/diff_null_lessErr50.csv
#           exp/exp_001/Results/cor_diff_null_lessErr50.csv
################################################################################

##################### Vijay's simulations ##############################
project.folder = paste0(print(here::here()),'/') 
project.folder
setwd(project.folder)
getwd()

library(readr)
library(MASS)
library(dplyr)

##################### Logging setup ###################################

log_file <- "./exp/exp_001/sims5_log.txt"

start_time <- Sys.time()

cat("Script: sims5.R\n", file = log_file)
cat("Start time:", start_time, "\n\n", file = log_file, append = TRUE)

#######################################################################
## 1. bring in "true data"

truedta_null <- read_csv("truedta_null_Final.csv")
truedta_ass  <- read_csv("truedta_ass_Final.csv")

truedta_null$term_births <- truedta_null$all_births - truedta_null$preterm_births
truedta_ass$term_births  <- truedta_ass$all_births - truedta_ass$preterm_births

################### PTB scenario ######################################

set.seed(212)
n.sims  <- 200
err.lev <- c(0, 0.5, 1, 2)

##################################################################
#################### DIFFERENTIAL ERROR Sims5.R ##################
##################################################################

# less error among those high-risk & 50% of the rest

diff.null.lessErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.null.lessErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

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
    
    # Generate exposure error for different risk groups
    pm.err.norisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, err.lev[e] * sd(truedta_null$pm25))
    
    pm.err.hirisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, (err.lev[e] / 10) * sd(truedta_null$pm25))
    
    # Weighted combination (50/50 mix for term births)
    pm.err <- ifelse(
      truedta_null$all_births > 0,
      (0.5 * truedta_null$term_births * pm.err.norisk +
         0.5 * truedta_null$term_births * pm.err.hirisk +
         truedta_null$preterm_births * pm.err.hirisk) /
        truedta_null$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )
    
    # Add to dataset
    truedta_null$pm_err <- pm.err
    
    # -------- RANDOM INTERCEPT quasi-Poisson model --------
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
    # ------------------------------------------------------
    
    if(!inherits(mod, "try-error")){
      
      diff.null.lessErr50[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
      diff.null.lessErr50[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
      
    } else {
      
      diff.null.lessErr50[s, (e*2 - 1)] <- NA
      diff.null.lessErr50[s, (e*2)]     <- NA
      
      fail_msg <- paste(
        "Model failed at simulation", s,
        "error level", err.lev[e], "\n"
      )
      
      cat(fail_msg)
      cat(fail_msg, file = log_file, append = TRUE)
      
    }
    
    # Correlation
    cor.diff.null.lessErr50[s, e] <- cor(pm.err, truedta_null$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
    
  }
  
  done_msg <- paste("Simulation", s, "completed at", Sys.time(), "\n")
  
  cat(done_msg)
  cat(done_msg, file = log_file, append = TRUE)
  
}

##################### Save results ###################################

diff.null.lessErr50 <- as.data.frame(diff.null.lessErr50)

names(diff.null.lessErr50) <- paste0(
  rep(c("b_", "se_"), times = length(err.lev)),
  rep(err.lev, each = 2)
)

cor.diff.null.lessErr50 <- as.data.frame(cor.diff.null.lessErr50)

names(cor.diff.null.lessErr50) <- paste0("corr", err.lev)

write_csv(diff.null.lessErr50, "./exp/exp_001/Results/diff_null_lessErr50.csv")
write_csv(cor.diff.null.lessErr50, "./exp/exp_001/Results/cor_diff_null_lessErr50.csv")

##################### End logging ####################################

end_time <- Sys.time()
runtime  <- end_time - start_time

cat("\nEnd time:", end_time, "\n", file = log_file, append = TRUE)
cat("Total runtime:", runtime, "\n", file = log_file, append = TRUE)

cat("\nScript finished.\nTotal runtime:", runtime, "\n")

####################### clear memory #################################

rm(diff.null.lessErr50, cor.diff.null.lessErr50)