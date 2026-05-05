################################################################################
# sims6.R  —  Scenario 6: Differential error (MORE error in at-risk + 50% non-risk), NULL
#
# Error structure (Section 2.2, Table 1, E1):
#   At-risk (preterm births):      v ~ N(0, (err.lev * 10 * SD(X))^2)
#   50% of non-at-risk (term):     v ~ N(0, (err.lev * 10 * SD(X))^2)  [more]
#   50% of non-at-risk (term):     v ~ N(0, (err.lev       * SD(X))^2)  [baseline]
#
# See sims5.R for details on the 50/50 weighting scheme.
#
# Usage   : Rscript sims6.R  /  source("sims6.R")
# Output  : exp/exp_001/Results/diff_null_moreErr50.csv
#           exp/exp_001/Results/cor_diff_null_moreErr50.csv
################################################################################

###################### Vijay's simulations ##############################

project.folder = paste0(print(here::here()),'/') 
project.folder
setwd(project.folder)
getwd()

library(readr)
library(MASS)
library(dplyr)

##################### Logging setup ###################################

log_file <- "./exp/exp_001/sims6_log.txt"

start_time <- Sys.time()

cat("Script: sims6.R\n", file = log_file)
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
#################### DIFFERENTIAL ERROR Sims6.R ##################
##################################################################

## more error among those high-risk & 50% of the rest

diff.null.moreErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.diff.null.moreErr50 <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

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
    
    # Generate exposure error
    pm.err.norisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, err.lev[e] * sd(truedta_null$pm25))
    
    pm.err.hirisk <- truedta_null$pm25 +
      rnorm(nrow(truedta_null), 0, (err.lev[e] * 10) * sd(truedta_null$pm25))
    
    # Weighted combination (50/50 mixture for term births)
    pm.err <- ifelse(
      truedta_null$all_births > 0,
      (0.5 * truedta_null$term_births * pm.err.norisk +
         0.5 * truedta_null$term_births * pm.err.hirisk +
         truedta_null$preterm_births * pm.err.hirisk) /
        truedta_null$all_births,
      mean(c(pm.err.norisk, pm.err.hirisk))
    )
    
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
      
      diff.null.moreErr50[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
      diff.null.moreErr50[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
      
    } else {
      
      diff.null.moreErr50[s, (e*2 - 1)] <- NA
      diff.null.moreErr50[s, (e*2)]     <- NA
      
      fail_msg <- paste(
        "Model failed at simulation", s,
        "error level", err.lev[e], "\n"
      )
      
      cat(fail_msg)
      cat(fail_msg, file = log_file, append = TRUE)
      
    }
    
    cor.diff.null.moreErr50[s, e] <- cor(pm.err, truedta_null$pm25)
    
    rm(mod, pm.err, pm.err.norisk, pm.err.hirisk)
    
  }
  
  done_msg <- paste("Simulation", s, "completed at", Sys.time(), "\n")
  
  cat(done_msg)
  cat(done_msg, file = log_file, append = TRUE)
  
}

##################### Save results ###################################

diff.null.moreErr50 <- as.data.frame(diff.null.moreErr50)

names(diff.null.moreErr50) <- paste0(
  rep(c("b_", "se_"), times = length(err.lev)),
  rep(err.lev, each = 2)
)

cor.diff.null.moreErr50 <- as.data.frame(cor.diff.null.moreErr50)

names(cor.diff.null.moreErr50) <- paste0("corr", err.lev)

write_csv(diff.null.moreErr50, "./exp/exp_001/Results/diff_null_moreErr50.csv")
write_csv(cor.diff.null.moreErr50, "./exp/exp_001/Results/cor_diff_null_moreErr50.csv")

##################### End logging ####################################

end_time <- Sys.time()
runtime  <- end_time - start_time

cat("\nEnd time:", end_time, "\n", file = log_file, append = TRUE)
cat("Total runtime:", runtime, "\n", file = log_file, append = TRUE)

cat("\nScript finished.\nTotal runtime:", runtime, "\n")

####################### clear memory #################################

rm(diff.null.moreErr50, cor.diff.null_moreErr50)