# Details for code reviewer
Use the code_review.R file to add your comments in case anything is added to the code

This file contains details of each experiment. author, exp_xxx, detil about simulations

To start of code review, you need;

1) Daily PM2.5 zip code level data  (pm25_js_zip_2000_2016_ny_daily.csv) from US-wide Daily 1 KM PM2.5 model by Di, Q.(2020) https://doi.org/10.7927/9yp5-hz11
2) Preterm birth aggregated data (zipall_1719.sas7bdat) from New York State Health https://www.health.ny.gov/statistics/vital_statistics

A) You need to set the working directory properly...like mine is setwd("~/Downloads/DiffMeasError/Code/codeReview"). Inside code review, we have data and code for data generation, and then an exp subfolder with simulations

There are three major steps of the code; 

1) Data creation: In this step, you need to create null and association data

2) Run simulations in each experiment in the exp folder, where three experiments are listed: exp_001, exp_002, and exp_003. There are 10 different simulations in each experiment with sims1.R, sims2.R, and up to sims10.R. Ideally, you can do it in a screen session because it will take around 3 hours on my Mac Pro M2, 16GB RAM, and 19-core GPU. You can run each R script manually if you want, but it will take too much time. You can open screen by screen -S sims1 for let say one screen, and inside the screen run Rscript sims1.R and repeat for all 10 screens. 

3) The results of all simulations will be saved in Results, and there is a visualization R script. Run that to get the required plots and compare them to the figures in the paper.

Data is already created using step #1. I would suggest going through step #3, and once everything looks good, generating the data by yourself using step #1 and then doing steps #2 to 3 again to make sure everything looks good. 

Figure 1--2, Table A2--A4 are the results from all three experiments.

# details of experiments
# there were around fifty experiments in total, but we have included the top three, which are used in the paper

details of each experiment: author, exp_xxx, details about simulations.
# exp_001, vijay, first experiment with multiple scripts running in screen sessions, scenario (Paper Scenario E-1)
# exp_002, vijay, added c with c= +2.5 in mean of PM25 (Paper Scenario E-2)
# exp_003, vijay, added c with c= -2.5 in mean of PM25 (Paper Scenario E-3)
