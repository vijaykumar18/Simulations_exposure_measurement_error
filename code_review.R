## Code Review for following paper
##Title: The impact of differential exposure measurement error on the PM2.5 -- preterm birth association: A simulation study
## Code Review: Alan Domínguez
## Updated: March 20, 2026

# Comment: Please add your comments

#A. Comment on each RScript

##################################################################################################################

#1. CreateDataFinal.R 

I suggest to install/load packages and create folders separetely (this is just in terms of reproducibility)

project.folder = paste0(print(here::here()),'/') 
source(paste0(project.folder,'setup_directory_structure.R'))
source(paste0(functions.folder,'script_initiate.R'))

# File: packages_to_use.R ----------------------------------------------------------------------

# List with all the packages used 
list.packages <- c('usethis','tidyverse', 'tidylog', 'sf', 'usmap','excessmort')

new.packages <- list.packages[!(list.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load all packages
lapply(list.packages, require, character.only = TRUE)

# File: script_initiate.R ---------------------------------------------------------------------

# sourcing the packages we use
source(paste0(packages.folder,'packages_to_use.R')) 

# load functions used in many scripts
source(paste0(functions.folder,'functions.R'))

I check the number of spatial units in the daily_pm25 dataset and are 574 ZIP.

Then for the simulation you just pick the ones that contain birth data (181 ZIP). 

For the STEP 2 — HARMFUL ASSOCIATION DATA (true RR = 1.115 per 10 ug/m3) I check "Association of Air Pollution and Heat Exposure With Preterm Birth,
Low Birth Weight, and Stillbirth in the US: A Systematic Review", the RR was correctly transfer.

Maybe you can add this into the code. 
  
increase_percent <- 11.5 # this is the number reported in the systematic review. 
RR <- 1 + increase_percent/100
RR # 1.115 

##################################################################################################################

I checked both ways to run the simulation study. 

#2. sims_all.R (I test wih a few number of simulations and it worked well)
sims_all.R is running withouth any issue. I think we can remove this part of the code: 
  
# ── Shared simulation parameters ─────────────────────────────────────────────
set.seed(212)
n.sims  <- 200
err.lev <- c(0, 0.5, 1, 2)

The simulation parameters were already setted in each script.  Just a minor addition to make fully reproducible the simuliton.

Add this to set the project folder (this depends if you upload a Rproject in github or if you just add the code). 

project.folder = paste0(print(here::here()),'/') 
project.folder
setwd(project.folder)
getwd() # this should match with project.folder name 


##########################################
### --- EXP 001, EXP 002, EXP 003 --- ###
########################################

# 2. sims1.R to sims10.R in screen session (I run the entire simulation using this apporach)

I really like the simulation progress addition in each of the simsx.R scripts (definitely necessary in simulation studies!!!). This was one of my previous comment glad 
that you already add this in the second version.

In each of the scripts, I made this tiny change since I was working in a Rproject "simulation-ap_PTB_EM" i used this code for rooting the folder (again this is purely style)

project.folder = paste0(print(here::here()),'/') 
project.folder
setwd(project.folder)


There were no major issues in EXP001 AND EXP002. I was able to reproduce all results

# EXP 003 
Simulation 94 failed in my computer "sims10.R" leading to different results in panel F 

Simulation 93/200 (46.5%) started at 2026-03-23 22:47:19.697884
Simulation 93 completed at 2026-03-23 22:47:45.783524 
Simulation 94/200 (47%) started at 2026-03-23 22:47:45.783896
Model failed at simulation 94 error level 2 
Simulation 94 completed at 2026-03-23 22:48:12.181075 
Simulation 95/200 (47.5%) started at 2026-03-23 22:48:12.181329
Simulation 95 completed at 2026-03-23 22:48:38.426851 

This can be checked in exp/exp_003/sims10_log.txt

##################################################################################################################

#3.VisualizeResults.R

Does it make sense to add a Figure of NYC? a shapefile with ZIPs and air pollution levels and PTB distributions? just as an overall picture of the data that we use 
in the simulation study? 

A general comment in all the scripts for visuallizing the results is that for loading the data is a little bit difficult to follow, 
I would comment a little bit more wich dataset are we using for generate the figures, since most of the datasets are called the same,
I suggest to use absolute paths to the experiments folders "exp/exp_001_results/data.csv" or change the name of the csv files (probably changing the paths is more easy)
  
I would change how the constant is displayed in panel E and F for both figures I believe we should displayed as -C if not is a litle bit misleading. Additionally, 
I suggest to add a subtitle for A, C, E (null association) and for B, D, F (Harmfull association) also in the right side of each figure you could specify, 
the experiment E1: random error; E2: Random error + c; E3: Random error - c. This will help the reader to understand the scenarios without the
necesite of checking the footnotes. 

I was not able to fully reproduce FIGURE 2 panel F and Figure 3 panel F. 

After meeting with Vijay, we fixed this by  adding na.rm = TRUE       

# old version 
for (m in seq_along(c.dta)) {
  dta <- as.data.frame(get(c.dta[m]))
  
  for (e in seq_along(err.lev)) {
    b_col <- which(names(dta) == paste0("b_", err.lev[e]))
    se_col <- which(names(dta) == paste0("se_", err.lev[e]))
    
    if (length(b_col) > 0 && length(se_col) > 0) {
      beta <- mean(dta[, b_col])                       
      vbar <- mean(dta[, se_col]^2)
      capb <- var(dta[, b_col])
      

# new version (fixed)
for (m in seq_along(c.dta)) {
  dta <- as.data.frame(get(c.dta[m]))
        
    for (e in seq_along(err.lev)) {
      b_col <- which(names(dta) == paste0("b_", err.lev[e]))
      se_col <- which(names(dta) == paste0("se_", err.lev[e]))
          
       if (length(b_col) > 0 && length(se_col) > 0) {
        beta <- mean(dta[, b_col], na.rm = TRUE)                       
        vbar <- mean(dta[, se_col]^2, na.rm = TRUE)
        capb <- var(dta[, b_col], na.rm = TRUE)  

##################################################################################################################

# CI_check_A1-A5.Rmd

updated necessary packages to run all chunks

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
library(glmmTMB)
library(ggplot2)
library(knitr) # this package was added
library(kableExtra) # this package was added 

        
When genererating the tables I found that tables dont follow the same order as manuscript.
Numbers were correctly transfered to the new version of the manuscript. 



##################################################################################################################

#C. Comments on overall code

Really nice and tidy code, most of my initial suggestion were changed in the second version. 

I suggest some small changes (mostly on how root the project)-this are pure style comments and it will depend on how you want to set the github project.

Also for the github I would name it something like simulation_PM25_PTB_exposure_measurement_error is longer but more specific to the manuscript. 

Finally, numbers in the manuscript should be double check (abstract) to ensure that all changes were transfered correctly. 






