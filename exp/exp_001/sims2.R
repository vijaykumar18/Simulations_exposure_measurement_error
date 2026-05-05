################################################################################
# sims2.R  —  Scenario 2: Non-differential error, HARMFUL association
#
# Purpose : Run 200 simulations adding classical (non-differential) measurement
#           error to the harmful association dataset (true RR = 1.115).
#           Non-differential error is expected to attenuate estimates toward
#           the null (attenuation bias).
#
# Error structure (Section 2.2, E1):
#   Z_it = X_it + v_it,  v_it ~ N(0, (err.lev * SD(X))^2)
#   Same error applied to all observations regardless of PTB status.
#
# Usage   : Rscript sims2.R          (standalone screen session)
#           source("sims2.R")        (called from Sims_all.R — data pre-loaded)
#
# Output  : exp/exp_001/Results/nondiff_harm.csv
#           exp/exp_001/Results/cor_nondiff_harm.csv
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

log_file <- "./exp/exp_001/sims2_log.txt"

start_time <- Sys.time()

cat("Script: sims2.R\n", file = log_file)
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
##################### NON-DIFFERENTIAL ERROR Sims2.R #############
##################################################################

### harmful association
nondiff.harm <- matrix(NA, nrow = n.sims, ncol = length(err.lev) * 2)
cor.nondiff.harm <- matrix(NA, nrow = n.sims, ncol = length(err.lev))

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
    
    # Add classical (non-differential) measurement error to true PM2.5
    # v_it ~ N(0, (err.lev[e] * SD(X))^2)
    # Mismeasured exposure
    pm.err <- truedta_ass$pm25 +
      rnorm(nrow(truedta_ass), 0, err.lev[e] * sd(truedta_ass$pm25))
    
    # Add to dataset
    truedta_ass$pm_err <- pm.err
    
    # Quasi-Poisson GLMM with random intercept by ZIP (Equation 1, Section 2.3)
    # Random intercept accounts for within-ZIP clustering over time.
    # Quasi-Poisson handles residual overdispersion in daily PTB counts.
    
    # -------- RANDOM INTERCEPT QUASI-POISSON MODEL --------
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
    # --------------------------------------------------------
    # coefficient and standard error (same position as glm)
    # check if model failed
    if(!inherits(mod, "try-error")){
      
      nondiff.harm[s, (e*2 - 1)] <- summary(mod)$tTable["pm_err", 1]
      nondiff.harm[s, (e*2)]     <- summary(mod)$tTable["pm_err", 2]
      
    } else {
      
      nondiff.harm[s, (e*2 - 1)] <- NA
      nondiff.harm[s, (e*2)]     <- NA
      
      fail_msg <- paste(
        "Model failed at simulation", s,
        "error level", err.lev[e], "\n"
      )
      
      cat(fail_msg)
      cat(fail_msg, file = log_file, append = TRUE)
      
    }
    # Pearson correlation between error-prone and true exposure (Table A1)
    # Correlation
    cor.nondiff.harm[s, e] <- cor(pm.err, truedta_ass$pm25)
    
    rm(mod, pm.err)
  }
  
  done_msg <- paste("Simulation", s, "completed at", Sys.time(), "\n")
  
  cat(done_msg)
  cat(done_msg, file = log_file, append = TRUE)
}

##################### Save results ###################################

nondiff.harm <- as.data.frame(nondiff.harm)

names(nondiff.harm) <- paste0(
  rep(c("b_", "se_"), times = length(err.lev)),
  rep(err.lev, each = 2)
)

cor.nondiff.harm <- as.data.frame(cor.nondiff.harm)

names(cor.nondiff.harm) <- paste0("corr", err.lev)

write_csv(nondiff.harm, "./exp/exp_001/Results/nondiff_harm.csv")
write_csv(cor.nondiff.harm, "./exp/exp_001/Results/cor_nondiff_harm.csv")

##################### End logging ####################################

end_time <- Sys.time()
runtime  <- end_time - start_time

cat("\nEnd time:", end_time, "\n", file = log_file, append = TRUE)
cat("Total runtime:", runtime, "\n", file = log_file, append = TRUE)

cat("\nScript finished.\nTotal runtime:", runtime, "\n")

####################### clear memory #################################

rm(nondiff.harm, cor.nondiff.harm)

