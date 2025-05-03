#Author: Vijay Kumar
#Project: Exposure measurement error
#Code Reviever: Atharv Jayprakash
##Feb 2025
## Code Review
## The impact of differential exposure measurement error on the PM2.5-- preterm birth association: a simulation study
## @Atharv Jayprakash
## Updated: Feb 1, 2025

#First you need to install these packages..
install.packages(c("survival", "ggplot2", "broom", "haven", "dplyr", "tidyr", 
                   "readr", "MASS", "lubridate", "lme4", "cowplot"))



#@Atharv: fill this with any error you get

#To start of code review you need 1) daily PM2.5 zip code level (pm25_js_zip_2000_2016_ny_daily.csv)  data 2) Pretermbirth aggregated data (zipall_1719.sas7bdat)
# You need to set working directory properly...like mine is setwd("~/Downloads/DiffMeasError/CodeReview") inside code review we have data and code for data generation and then exp subfolder with simulations


#Details of each experiment. author, exp_xxx, detil about simulations

#Step-0 Create data; Run CreateData.R it will create two csv files namely truedta_null_rpois.csv and truedta_ass_rpois.csv

#Step-1 Preprocess.R: It creates folder structure copying from last exp folder (not required in code review)

#Step-2 Add details to each simX.R like (not required in code review)

#Step-3 Run the simsx.R in each exp_xxxx folder in screen session to make process faster by screen-S sims1 to open screen and cd directory where sims1.R is located and Rscript sims1.R in screen sims1. Repeat it for ten scripts in 10 screens. 

#Step-4 Run VisualizeResults.R in results folder of each experiment to get plots for the experiments 

#Note Directory: DiffMeasError is the main folder but everything else in inside the exp subfolder 

################################################Please fill the details##########################


#exp_001, vijay, first experiment with multiple scripts runnin in screen sessions, scenario      (Paper Scenaro-1)
#exp_002, vijay, added +0.25 in mean of PM25 (Paper Scenaro-2)
#exp_003, vijay, added -0.25 in mean of PM25 (Paper Scenaro-3)
