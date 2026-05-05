################################################################################
# sims8.R  —  Scenario 8: Differential error (MORE error in at-risk), HARMFUL
#
# Error structure: At-risk gets 10x MORE error than non-at-risk.
# Applied to harmful association dataset (true RR = 1.115).
#
# Usage   : Rscript sims8.R  /  source("sims8.R")
# Output  : exp/exp_001/Results/diff_harm_moreErr.csv
#           exp/exp_001/Results/cor_diff_harm_moreErr.csv
################################################################################



##################### Vijay's simulations ##############################

setwd("~/Downloads/CodeReview_Alan")

library(readr)
library(MASS)
library(dplyr)

##################### Logging setup ###################################

log_file <- "./exp/exp_001/sims8_log.txt"

start_time <- Sys.time()

cat("Script: sims8.R\n", file = log_file)
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
#################### DIFFERENTIAL ERROR Sims8.R ##################
##################################################################

## more error among those high-risk

diff.harm.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.harm.moreErr <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

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
    
    # Generate differential more error for high-risk
    pm.err.norisk <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), 0, err.lev[e] * sd(truedta_ass$pm25))
    
    pm.err.hirisk <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), 0, (err.lev[e] * 10) * sd(truedta_ass$pm25))
    
    # Weighted PM by birth counts
    pm.err <- ifelse(
      truedta_ass$all_births > 0,
      (truedta_ass$term_births * pm.err.norisk +
         truedta_ass$preterm_births * pm.err.hirisk) /
        truedta_ass$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )
    
    truedta_ass$pm_err <- pm.err
    
    # -------- RANDOM INTERCEPT quasi-Poisson model --------
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
    # ------------------------------------------------------
    
    if(!inherits(mod, "try-error")){
      
      diff.harm.moreErr[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
      diff.harm.moreErr[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
      
    } else {
      
      diff.harm.moreErr[s, (e*2 - 1)] <- NA
      diff.harm.moreErr[s, (e*2)]     <- NA
      
      fail_msg <- paste(
        "Model failed at simulation", s,
        "error level", err.lev[e], "\n"
      )
      
      cat(fail_msg)
      cat(fail_msg, file = log_file, append = TRUE)
      
    }
    
    cor.diff.harm.moreErr[s, e] <- cor(pm.err, truedta_ass$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
    
  }
  
  done_msg <- paste("Simulation", s, "completed at", Sys.time(), "\n")
  
  cat(done_msg)
  cat(done_msg, file = log_file, append = TRUE)
  
}

##################### Save results ###################################

diff.harm.moreErr <- as.data.frame(diff.harm.moreErr)

names(diff.harm.moreErr) <- paste0(
  rep(c("b_", "se_"), times = length(err.lev)),
  rep(err.lev, each = 2)
)

cor.diff.harm.moreErr <- as.data.frame(cor.diff.harm.moreErr)

names(cor.diff.harm.moreErr) <- paste0("corr", err.lev)

write_csv(diff.harm.moreErr, "./exp/exp_001/Results/diff_harm_moreErr.csv")
write_csv(cor.diff.harm.moreErr, "./exp/exp_001/Results/cor_diff_harm_moreErr.csv")

##################### End logging ####################################

end_time <- Sys.time()
runtime  <- end_time - start_time

cat("\nEnd time:", end_time, "\n", file = log_file, append = TRUE)
cat("Total runtime:", runtime, "\n", file = log_file, append = TRUE)

cat("\nScript finished.\nTotal runtime:", runtime, "\n")

####################### clear memory #################################

rm(diff.harm.moreErr, cor.diff.harm.moreErr)