################################################################################
# Sims_all.R  —  Run all 10 simulation scenarios for Experiment 1 (E1)
#
# Purpose : Convenience script that runs all 10 scenarios sequentially in a
#           single R session.  Data are loaded ONCE here and shared across all
#           scenario scripts, avoiding redundant file I/O.
#
# Individual scripts (sims1.R – sims10.R) can also be run independently in
# separate screen sessions for parallel execution (~2.5 hrs total vs ~10-12 hrs
# sequentially).  Each script guards against re-loading data if objects are
# already present in the environment (see the if(!exists(...)) blocks).
#
# Usage (single screen session, ~10-12 hrs on M2 Mac / 16 GB RAM):
#   screen -S sims_all
#   cd exp/exp_001
#   Rscript Sims_all.R
#   [Ctrl+A, D to detach]
#
# Usage (parallel screen sessions, ~2.5 hrs total):
#   Run sims1.R through sims10.R each in its own screen session (see readme).
#
# Outputs : All CSV results saved to ./exp/exp_001/Results/
################################################################################

# ── Working directory ─────────────────────────────────────────────────────────
# Set this to the project root before running.
# e.g. setwd("~/Downloads/CodeReview/CodeReview_Alan")

# ── Libraries ─────────────────────────────────────────────────────────────────
library(dplyr)
library(readr)
library(MASS)    # glmmPQL
library(lme4)

rm(list = ls())

# ── Load data ONCE ────────────────────────────────────────────────────────────
# Data are loaded here so each sourced scenario script skips its own load block.
# This avoids reading ~500 MB from disk 10 times when running all scenarios.
cat("Loading data...\n")

truedta_null <- read_csv("truedta_null_Final.csv")
truedta_ass  <- read_csv("truedta_ass_Final.csv")

truedta_null$term_births <- truedta_null$all_births - truedta_null$preterm_births
truedta_ass$term_births  <- truedta_ass$all_births  - truedta_ass$preterm_births

cat("Rows (null):", nrow(truedta_null), "| Rows (ass):", nrow(truedta_ass), "\n")

# ── Shared simulation parameters ─────────────────────────────────────────────
set.seed(212)
n.sims  <- 200
err.lev <- c(0, 0.5, 1, 2)

getwd()

# ── Run all 10 scenarios ─────────────────────────────────────────────────────
# Each script detects that truedta_null / truedta_ass already exist and skips
# its standalone data-loading block, using the shared objects above instead.

cat("\n[1/10] Scenario 1: Non-differential error, Null\n")
source("exp/exp_001/sims1.R")

cat("\n[2/10] Scenario 2: Non-differential error, Harmful\n")
source("exp/exp_001/sims2.R")

cat("\n[3/10] Scenario 3: Differential less error (at-risk only), Null\n")
source("exp/exp_001/sims3.R")

cat("\n[4/10] Scenario 4: Differential more error (at-risk only), Null\n")
source("exp/exp_001/sims4.R")

cat("\n[5/10] Scenario 5: Differential less error (at-risk + 50% non-risk), Null\n")
source("exp/exp_001/sims5.R")

cat("\n[6/10] Scenario 6: Differential more error (at-risk + 50% non-risk), Null\n")
source("exp/exp_001/sims6.R")

cat("\n[7/10] Scenario 7: Differential less error (at-risk only), Harmful\n")
source("exp/exp_001/sims7.R")

cat("\n[8/10] Scenario 8: Differential more error (at-risk only), Harmful\n")
source("exp/exp_001/sims8.R")

cat("\n[9/10] Scenario 9: Differential less error (at-risk + 50% non-risk), Harmful\n")
source("exp/exp_001/sims9.R")

cat("\n[10/10] Scenario 10: Differential more error (at-risk + 50% non-risk), Harmful\n")
source("exp/exp_001/sims10.R")

cat("\nAll scenarios complete. Results saved to ./exp/exp_001/Results/\n")
sessionInfo()
