# Details for code reviewer
Use code_review.R file to aadd your comments in case anything is added to code.

# Code
This file contains details of each experiment. author, exp_xxx, detil about simulations

To start of code review you need;

1) daily PM2.5 zip code level data  (pm25_js_zip_2000_2016_ny_daily.csv)  
2) Pretermbirth aggregated data (zipall_1719.sas7bdat)

# You need to set working directory properly...like mine is setwd("~/Downloads/DiffMeasError/Code") inside code review we have data and code for data generation and then exp subfolder with simulations

There are three major steps of the code; 

1) Data creation: in this step, you need to create null and association data

2) Run simulations in each experiment in the exp folder where three experiments are listed: exp_001, exp_002, and exp_003. There are 10 different simulations in each experiment with sims1.R, sims2.R, and up to sims10.R. Ideally, you can do it in a screen session because it will take around 3 hours on my Mac Pro M2, 16GB RAM, and 19-core GPU. You can run each R script manually if you want, but it will take too much time. You can open screen by screen -S sims1 for let say one screen, and inside the screen run Rscript sims1.R and repeat for all 10 screens. 

3) The results of all simulations will be saved in Results, and there is a visualization R script. Run that to get the required plots and compare them to the figures in the paper.

I am also attaching data I have already created using step #1. I would suggest going through step #3, and once everything looks good, generating the data by yourself using step #1 and then doing step #2 to 3 again to make sure everything looks good. 

The Figure 2 to 5 in paper is coming from experiment rxp_001, Figure 6 to 9 in paper is coming from experiment exp_002, and Figure 10 to 13 in paper is coming from experiment exp_003,

#Note Directory: DiffMeasError is the main folder but everything else in inside the exp subfolder 

# details of experiments
# Note there were around fifty experiments in total but we have included top three which are used in paper.

Note details of each experiment. author, exp_xxx, detil about simulations.
#exp_001, vijay, first experiment with multiple scripts runnin in screen sessions, scenario      (Paper Scenaro-1)
#exp_002, vijay, added +0.25 in mean of PM25 (Paper Scenaro-2)
#exp_003, vijay, added -0.25 in mean of PM25 (Paper Scenaro-3)
