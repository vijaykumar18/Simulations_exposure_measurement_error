################################################################################
# sims3.R  —  Scenario 3: Differential error (LESS error in at-risk), NULL
#
# Purpose : Run 200 simulations where the at-risk (preterm) population has
#           LESS measurement error than the non-at-risk population.
#           Applied to the null dataset (true RR = 1.00).
#
# Error structure (Section 2.2, Table 1, E1):
#   Non-at-risk (term births):  v ~ N(0, (err.lev       * SD(X))^2)
#   At-risk (preterm births):   v ~ N(0, (err.lev / 10  * SD(X))^2)
#
# Aggregation: The ZIP-day-level observed exposure is a birth-count-weighted
#   average of the two group-specific error-prone exposures:
#     Z_it = (term_births * Z_norisk + preterm_births * Z_hirisk) / all_births
#   This reflects the fact that the observed ZIP-level exposure is an aggregate
#   of individual-level exposures across births with different error variances.
#
# Usage   : Rscript sims3.R
#           source("sims3.R")  (from Sims_all.R)
#
# Output  : exp/exp_001/Results/diff_null_lessErr.csv
#           exp/exp_001/Results/cor_diff_null_lessErr.csv
################################################################################



##################### Vijay's simulations ##############################

setwd("~/Downloads/CodeReview_Alan")

library(readr)
library(MASS)
library(dplyr)

##################### Logging setup ###################################

log_file <- "./exp/exp_001/sims3_log.txt"

start_time <- Sys.time()

cat("Script: sims3.R\n", file = log_file)
cat("Start time:", start_time, "\n\n", file = log_file, append = TRUE)

#######################################################################
## 1. bring in "true data"

truedta_null <- read_csv("truedta_null_Final.csv")
truedta_ass  <- read_csv("truedta_ass_Final.csv")

truedta_null$term_births <- truedta_null$all_births - truedta_null$preterm_births
truedta_ass$term_births  <- truedta_ass$all_births - truedta_ass$preterm_births

################### PTB scenario ######################################

set.seed(212)
n.sims  <- 5
err.lev <- c(0, 0.5, 1, 2)

##################################################################
#################### DIFFERENTIAL ERROR Sims3.R ##################
##################################################################

## less error among those high-risk

diff.null.lessErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.null.lessErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

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
    
    # Generate exposure error separately by risk group
    pm.err.norisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, err.lev[e] * sd(truedta_null$pm25))
    
    pm.err.hirisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, (err.lev[e] / 10) * sd(truedta_null$pm25))
    
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
    
    # Quasi-Poisson GLMM with random intercept by ZIP (Equation 1, Section 2.3)
    # Random intercept accounts for within-ZIP clustering over time.
    # Quasi-Poisson handles residual overdispersion in daily PTB counts.
    
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
      
      diff.null.lessErr[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
      diff.null.lessErr[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
      
    } else {
      
      diff.null.lessErr[s, (e*2 - 1)] <- NA
      diff.null.lessErr[s, (e*2)]     <- NA
      
      fail_msg <- paste(
        "Model failed at simulation", s,
        "error level", err.lev[e], "\n"
      )
      
      cat(fail_msg)
      cat(fail_msg, file = log_file, append = TRUE)
    }
    
    # correlation
    cor.diff.null.lessErr[s, e] <- cor(pm.err, truedta_null$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
  }
  
  done_msg <- paste("Simulation", s, "completed at", Sys.time(), "\n")
  
  cat(done_msg)
  cat(done_msg, file = log_file, append = TRUE)
}

##################### Save results ###################################

diff.null.lessErr <- as.data.frame(diff.null.lessErr)

names(diff.null.lessErr) <- paste0(
  rep(c("b_", "se_"), times = length(err.lev)),
  rep(err.lev, each = 2)
)

cor.diff.null.lessErr <- as.data.frame(cor.diff.null.lessErr)

names(cor.diff.null.lessErr) <- paste0("corr", err.lev)

write_csv(diff.null.lessErr, "./exp/exp_001/Results/diff_null_lessErr.csv")
write_csv(cor.diff.null.lessErr, "./exp/exp_001/Results/cor_diff_null_lessErr.csv")

##################### End logging ####################################

end_time <- Sys.time()
runtime  <- end_time - start_time

cat("\nEnd time:", end_time, "\n", file = log_file, append = TRUE)
cat("Total runtime:", runtime, "\n", file = log_file, append = TRUE)

cat("\nScript finished.\nTotal runtime:", runtime, "\n")

####################### clear memory #################################

rm(diff.null.lessErr, cor.diff.null.lessErr)
