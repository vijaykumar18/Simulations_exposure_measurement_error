setwd("~/Downloads/CodeReview/CodeReview_Heather/exp/exp_002/Results")
library(survival)
library(ggplot2)
library(broom)
library(haven)
library(dplyr)
#library(tidyr)
library(readr)
library(MASS)
library(lubridate)
library(lme4)
library(cowplot) # for combining plots

rm(list=ls())

# This scenario combines 0 error from baseline E1 experiment and +c from E2 experiment. 
# Process exp_043 (baseline - ALL error levels 0, 0.5, 1, 2 but WITHOUT +c suffix)
setwd("/scratch/user/u.vk230134/Downloads/DiffMeasError/CodeReview/exp/exp_043/Results")

# Load baseline datasets
nondiff_null      <- read_csv("nondiff_null.csv")
diff_null_moreErr <- read_csv("diff_null_moreErr.csv")
diff_null_lessErr <- read_csv("diff_null_lessErr.csv")
diff_null_lessErr50 <- read_csv("diff_null_lessErr50.csv")
diff_null_moreErr50 <- read_csv("diff_null_moreErr50.csv")

nondiff_harm      <- read_csv("nondiff_harm.csv")
diff_harm_moreErr <- read_csv("diff_harm_moreErr.csv")
diff_harm_lessErr <- read_csv("diff_harm_lessErr.csv")
diff_harm_lessErr50 <- read_csv("diff_harm_lessErr50.csv")
diff_harm_moreErr50 <- read_csv("diff_harm_moreErr50.csv")

# List of baseline dataset names
baseline.dta <- ls()

# Process exp_044 (+c datasets - all error levels)
setwd("/scratch/user/u.vk230134/Downloads/DiffMeasError/CodeReview/exp/exp_044/Results")

# Load +c datasets
nondiff_null_c      <- read_csv("nondiff_null.csv")
diff_null_moreErr_c <- read_csv("diff_null_moreErr.csv")
diff_null_lessErr_c <- read_csv("diff_null_lessErr.csv")
diff_null_lessErr50_c <- read_csv("diff_null_lessErr50.csv")
diff_null_moreErr50_c <- read_csv("diff_null_moreErr50.csv")

nondiff_harm_c      <- read_csv("nondiff_harm.csv")
diff_harm_moreErr_c <- read_csv("diff_harm_moreErr.csv")
diff_harm_lessErr_c <- read_csv("diff_harm_lessErr.csv")
diff_harm_lessErr50_c <- read_csv("diff_harm_lessErr50.csv")
diff_harm_moreErr50_c <- read_csv("diff_harm_moreErr50.csv")

# List of +c dataset names
c.dta <- ls()[!ls() %in% baseline.dta]

# Error levels for both baseline and +c datasets
err.lev <- c(0, 0.5, 1, 2)

# Initialize results data frames
results <- data.frame()
bias_results <- data.frame()


# PROCESS BASELINE (exp_043) - ONLY ERROR LEVEL 0
for (m in seq_along(baseline.dta)) {
  dta <- as.data.frame(get(baseline.dta[m]))
  
  # ONLY process error level 0 for baseline
  e <- 1  # This corresponds to error level 0
  b_col <- which(names(dta) == paste0("b_", err.lev[e]))
  se_col <- which(names(dta) == paste0("se_", err.lev[e]))
  
  if (length(b_col) > 0 && length(se_col) > 0) {
    beta <- mean(dta[, b_col])
    vbar <- mean(dta[, se_col]^2)
    capb <- var(dta[, b_col])
    sterr <- sqrt(vbar + (nrow(dta) + 1) / nrow(dta) * capb)
    rr <- exp(beta)
    lci <- exp(beta - 1.96 * sterr)
    uci <- exp(beta + 1.96 * sterr)
    
    temp.res <- data.frame(beta, sterr, rr, lci, uci)
    temp.res$scenario <- paste(baseline.dta[m], err.lev[e], sep = "_")
    
    results <- bind_rows(results, temp.res)
    
    # Calculate bias
    beta_values <- dta[, b_col]
    if (grepl("harm", baseline.dta[m])) {
      result_type <- "Association"
      bias_values <- (beta_values - log(1.115)) / log(1.115) * 100
      # COVERAGE: For association, check if true beta (log(1.115)) is within CI
      true_beta <- log(1.115)
    } else {
      result_type <- "Null Association" 
      bias_values <- beta_values
      # COVERAGE: For null, check if 0 is within CI
      true_beta <- 0
    }
    
    bias_mean <- mean(bias_values)
    bias_2.5 <- quantile(bias_values, 0.025)
    bias_97.5 <- quantile(bias_values, 0.975)
    
    # CALCULATE COVERAGE
    ci_lower <- beta_values - 1.96 * dta[, se_col]
    ci_upper <- beta_values + 1.96 * dta[, se_col]
    coverage <- mean(ci_lower <= true_beta & true_beta <= ci_upper) * 100
    
    temp.bias <- data.frame(
      scenario = paste(baseline.dta[m], err.lev[e], sep = "_"),
      bias_mean = bias_mean,
      bias_2.5 = bias_2.5,
      bias_97.5 = bias_97.5,
      coverage = coverage,  # ADD COVERAGE
      result_type = result_type
    )
    
    bias_results <- bind_rows(bias_results, temp.bias)
  }
}

# PROCESS +C DATASETS (exp_045) - ALL ERROR LEVELS (0, 0.5, 1, 2) WITH +c
for (m in seq_along(c.dta)) {
  dta <- as.data.frame(get(c.dta[m]))
  
  for (e in seq_along(err.lev)) {
    b_col <- which(names(dta) == paste0("b_", err.lev[e]))
    se_col <- which(names(dta) == paste0("se_", err.lev[e]))
    
    if (length(b_col) > 0 && length(se_col) > 0) {
      beta <- mean(dta[, b_col])
      vbar <- mean(dta[, se_col]^2)
      capb <- var(dta[, b_col])
      sterr <- sqrt(vbar + (nrow(dta) + 1) / nrow(dta) * capb)
      rr <- exp(beta)
      lci <- exp(beta - 1.96 * sterr)
      uci <- exp(beta + 1.96 * sterr)
      
      temp.res <- data.frame(beta, sterr, rr, lci, uci)
      temp.res$scenario <- paste0(gsub("_c$", "", c.dta[m]), "_", err.lev[e], "+c")
      
      results <- bind_rows(results, temp.res)
      
      # Calculate bias
      beta_values <- dta[, b_col]
      if (grepl("harm", c.dta[m])) {
        result_type <- "Association"
        bias_values <- (beta_values - log(1.115)) / log(1.115) * 100
        true_beta <- log(1.115)
      } else {
        result_type <- "Null Association"
        bias_values <- beta_values
        true_beta <- 0
      }
      
      bias_mean <- mean(bias_values)
      bias_2.5 <- quantile(bias_values, 0.025)
      bias_97.5 <- quantile(bias_values, 0.975)
      
      # CALCULATE COVERAGE
      ci_lower <- beta_values - 1.96 * dta[, se_col]
      ci_upper <- beta_values + 1.96 * dta[, se_col]
      coverage <- mean(ci_lower <= true_beta & true_beta <= ci_upper) * 100
      
      temp.bias <- data.frame(
        scenario = paste0(gsub("_c$", "", c.dta[m]), "_", err.lev[e], "+c"),
        bias_mean = bias_mean,
        bias_2.5 = bias_2.5,
        bias_97.5 = bias_97.5,
        coverage = coverage,  # ADD COVERAGE
        result_type = result_type
      )
      
      bias_results <- bind_rows(bias_results, temp.bias)
    }
  }
}

coverage_results <- bias_results %>% 
  dplyr::select(scenario, coverage, result_type)

summary(results)
summary(bias_results)
summary(coverage_results)

# Save final results
write_csv(results, "results.csv")
write_csv(bias_results, "bias_results.csv")
write_csv(coverage_results, "coverage_results.csv")
###########################################
##E2 null Results


# Get E1 null scenarios (error level 0 only)
e1.null.res <- results[grep("null.*_0$", results$scenario),]  # Ends with "_0" (no +c)

# Get E2 null scenarios (all +c error levels)  
e2.null.res <- results[grep("null.*\\+c", results$scenario),]

# Combine E1 and E2
e2.combined.res <- bind_rows(e1.null.res, e2.null.res)

# Define the order by SCENARIO TYPE first, then error level
desired_order <- c(
  # Diff Less (5 points)
  "diff_null_lessErr_0", "diff_null_lessErr_0+c", "diff_null_lessErr_0.5+c", "diff_null_lessErr_1+c", "diff_null_lessErr_2+c",
  # Diff Less +50% (5 points)
  "diff_null_lessErr50_0", "diff_null_lessErr50_0+c", "diff_null_lessErr50_0.5+c", "diff_null_lessErr50_1+c", "diff_null_lessErr50_2+c",
  # Diff More (5 points)
  "diff_null_moreErr_0", "diff_null_moreErr_0+c", "diff_null_moreErr_0.5+c", "diff_null_moreErr_1+c", "diff_null_moreErr_2+c",
  # Diff More +50% (5 points)
  "diff_null_moreErr50_0", "diff_null_moreErr50_0+c", "diff_null_moreErr50_0.5+c", "diff_null_moreErr50_1+c", "diff_null_moreErr50_2+c",
  # Nondiff (5 points)
  "nondiff_null_0", "nondiff_null_0+c", "nondiff_null_0.5+c", "nondiff_null_1+c", "nondiff_null_2+c"
)


# Reorder the data
e2.combined.res <- e2.combined.res[match(desired_order, e2.combined.res$scenario), ]

# Apply styling
e2.combined.res$sym <- ifelse(grepl("nondiff", e2.combined.res$scenario), 15, 19)
e2.combined.res$col <- ifelse(grepl("less", e2.combined.res$scenario) & !grepl("more", e2.combined.res$scenario), "blue", 
                              ifelse(grepl("more", e2.combined.res$scenario), "cadetblue", "black"))
e2.combined.res$lty <- ifelse(grepl("Err50", e2.combined.res$scenario), 2, 1)
e2.combined.res$lwd <- ifelse(grepl("Err50", e2.combined.res$scenario), 1.5, 1)

# Extract error levels for x-axis labels (simplified)
e2.combined.res$err <- ifelse(
  grepl("_0$", e2.combined.res$scenario), "0",
  ifelse(grepl("_0\\+c", e2.combined.res$scenario), "0+c",
         ifelse(grepl("_0\\.5\\+c", e2.combined.res$scenario), "0.5+c",
                ifelse(grepl("_1\\+c", e2.combined.res$scenario), "1+c", "2+c")))
)

summary(e2.null.res)

pdf("null_results.pdf", width = 6, height = 6)
par(mar = c(5, 5, 4, 2) + 0.01)

# For log scale, we need to adjust the y-axis range and use log="y"
y_min <- 0.99
y_max <- 1.16

# Use log scale on y-axis - add log="y" parameter
plot(c(1, nrow(e2.combined.res)), c(y_min, y_max), pch = '', las = 1,
     xlab = "", ylab = "RR & 95% CI", cex.lab = 1.5, cex.axis = 1.25, xaxt = "n",
     log = "y")  # This enables log scale on y-axis

abline(h = 1, col = "darkgray")

# X-axis labels (error levels)
axis(side = 1, at = 1:nrow(e2.combined.res), labels = e2.combined.res$err, cex.axis = 1.1, las = 2)

# Add vertical lines to separate scenario groups
group_breaks <- c(5, 10, 15, 20)  # After each scenario group (5 points each)
#abline(v = group_breaks + 0.5, col = "gray", lty = 3)

# Plot points and CIs (no changes needed here - they'll automatically use log scale)
points(1:nrow(e2.combined.res), e2.combined.res$rr, 
       pch = ifelse(grepl("Err50", e2.combined.res$scenario), 23, e2.combined.res$sym),
       col = e2.combined.res$col, 
       cex = 1.5)

segments(
  1:nrow(e2.combined.res), e2.combined.res$lci,
  1:nrow(e2.combined.res), e2.combined.res$uci,
  col = e2.combined.res$col, 
  lty = e2.combined.res$lty, 
  lwd = e2.combined.res$lwd
)

# Add scenario group labels
scenario_labels <- c("Diff Less", "Diff Less\n+50%", "Diff More", "Diff More\n+50%", "Nondiff")
label_positions <- c(3, 8, 13, 18, 23)  # Middle of each group


# Legend
legend("topright", 
       pch = c(19, 23, 19, 23, 15), 
       col = c("blue", "blue", "cadetblue", "cadetblue", "black"), 
       legend = c("Diff Less", "Diff Less + 50%", "Diff More", "Diff More + 50%", "Nondiff"), 
       bty = "n")

legend("bottomleft", lty = c(1, 2), 
       legend = c("High-Risk", "High-Risk + 50% Rest"), 
       bty = 'n')

mtext("Error Level", side = 1, line = 3.5, cex = 1.5)
dev.off()

## E2 Harm Results - GROUPED BY SCENARIO TYPE
# Get E1 harm scenarios (error level 0 only)
e1.harm.res <- results[grep("harm.*_0$", results$scenario),]  # Ends with "_0" (no +c)

# Get E2 harm scenarios (all +c error levels)  
e2.harm.res <- results[grep("harm.*\\+c", results$scenario),]

# Combine E1 and E2
e2.combined.res <- bind_rows(e1.harm.res, e2.harm.res)

# Define the order by SCENARIO TYPE first, then error level
desired_order <- c(
  # Diff Less (5 points)
  "diff_harm_lessErr_0", "diff_harm_lessErr_0+c", "diff_harm_lessErr_0.5+c", "diff_harm_lessErr_1+c", "diff_harm_lessErr_2+c",
  # Diff Less +50% (5 points)
  "diff_harm_lessErr50_0", "diff_harm_lessErr50_0+c", "diff_harm_lessErr50_0.5+c", "diff_harm_lessErr50_1+c", "diff_harm_lessErr50_2+c",
  # Diff More (5 points)
  "diff_harm_moreErr_0", "diff_harm_moreErr_0+c", "diff_harm_moreErr_0.5+c", "diff_harm_moreErr_1+c", "diff_harm_moreErr_2+c",
  # Diff More +50% (5 points)
  "diff_harm_moreErr50_0", "diff_harm_moreErr50_0+c", "diff_harm_moreErr50_0.5+c", "diff_harm_moreErr50_1+c", "diff_harm_moreErr50_2+c",
  # Nondiff (5 points)
  "nondiff_harm_0", "nondiff_harm_0+c", "nondiff_harm_0.5+c", "nondiff_harm_1+c", "nondiff_harm_2+c"
)

# Reorder the data
e2.combined.res <- e2.combined.res[match(desired_order, e2.combined.res$scenario), ]

# Apply styling
e2.combined.res$sym <- ifelse(grepl("nondiff", e2.combined.res$scenario), 15, 19)
e2.combined.res$col <- ifelse(grepl("less", e2.combined.res$scenario) & !grepl("more", e2.combined.res$scenario), "blue", 
                              ifelse(grepl("more", e2.combined.res$scenario), "cadetblue", "black"))
e2.combined.res$lty <- ifelse(grepl("Err50", e2.combined.res$scenario), 2, 1)
e2.combined.res$lwd <- ifelse(grepl("Err50", e2.combined.res$scenario), 1.5, 1)

# Extract error levels for x-axis labels (simplified)
e2.combined.res$err <- ifelse(
  grepl("_0$", e2.combined.res$scenario), "0",
  ifelse(grepl("_0\\+c", e2.combined.res$scenario), "0+c",
         ifelse(grepl("_0\\.5\\+c", e2.combined.res$scenario), "0.5+c",
                ifelse(grepl("_1\\+c", e2.combined.res$scenario), "1+c", "2+c")))
)

pdf("harm_results.pdf", width = 6, height = 6)
par(mar = c(5, 5, 4, 2) + 0.01)

# For log scale, we need to adjust the y-axis range and use log="y"
y_min <- 0.99
y_max <- 1.27

# Use log scale on y-axis - add log="y" parameter
plot(c(1, nrow(e2.combined.res)), c(y_min, y_max), pch = '', las = 1,
     xlab = "", ylab = "RR & 95% CI", cex.lab = 1.5, cex.axis = 1.25, xaxt = "n",
     log = "y")  # This enables log scale on y-axis

abline(h = 1, col = "darkgray")
abline(h = 1.115, col = "hotpink2", lty = 4)

# X-axis labels (error levels)
axis(side = 1, at = 1:nrow(e2.combined.res), labels = e2.combined.res$err, cex.axis = 1.1, las = 2)


# Add vertical lines to separate scenario groups
group_breaks <- c(5, 10, 15, 20)  # After each scenario group (5 points each)
#abline(v = group_breaks + 0.5, col = "gray", lty = 3)

# Plot points and CIs (no changes needed here - they'll automatically use log scale)
points(1:nrow(e2.combined.res), e2.combined.res$rr, 
       pch = ifelse(grepl("Err50", e2.combined.res$scenario), 23, e2.combined.res$sym),
       col = e2.combined.res$col, 
       cex = 1.5)

segments(
  1:nrow(e2.combined.res), e2.combined.res$lci,
  1:nrow(e2.combined.res), e2.combined.res$uci,
  col = e2.combined.res$col, 
  lty = e2.combined.res$lty, 
  lwd = e2.combined.res$lwd
)

# Add scenario group labels
scenario_labels <- c("Diff Less", "Diff Less\n+50%", "Diff More", "Diff More\n+50%", "Nondiff")
label_positions <- c(3, 8, 13, 18, 23)  # Middle of each group


# Legend
legend("topright", 
       pch = c(19, 23, 19, 23, 15), 
       col = c("blue", "blue", "cadetblue", "cadetblue", "black"), 
       legend = c("Diff Less", "Diff Less + 50%", "Diff More", "Diff More + 50%", "Nondiff"), 
       bty = "n")

legend("bottomleft", lty = c(1, 2), 
       legend = c("High-Risk", "High-Risk + 50% Rest"), 
       bty = 'n')

mtext("Error Level", side = 1, line = 3.5, cex = 1.5)

dev.off()



# ############################Bias Heat maps ####################

bias_results <- read_csv("bias_results.csv")

# Add a column to indicate the result type

bias_assoc<-bias_results%>%dplyr::select(scenario,bias_mean,result_type)%>%filter(result_type=="Association")
bias_null<-bias_results%>%dplyr::select(scenario,bias_mean,result_type)%>%filter(result_type=="Null Association")


# Modify bias_type levels and prepare the data
bias_assoc <- bias_assoc %>%
  mutate(bias_type = factor(gsub("_\\d.*", "", scenario),
                            levels = c("nondiff_harm", 
                                       "diff_harm_lessErr", 
                                       "diff_harm_lessErr50",
                                       "diff_harm_moreErr",  
                                       "diff_harm_moreErr50")),
         error_level = factor(gsub(".*_", "", scenario), 
                              levels = c("0","0+c","0.5+c", "1+c", "2+c")))


bias_null <- bias_null %>%
  mutate(bias_type = factor(gsub("_\\d.*", "", scenario),
                            levels = c("nondiff_null", 
                                       "diff_null_lessErr", 
                                       "diff_null_lessErr50",
                                       "diff_null_moreErr",  
                                       "diff_null_moreErr50")),
         error_level = factor(gsub(".*_", "", scenario), 
                              levels = c("0","0+c","0.5+c", "1+c", "2+c")))




summary(bias_null)
summary(bias_assoc)
# Define color limits based on the type of bias
assoc_limits <- c(-134, 115)    # Adjust these limits for association
null_limits <- c(-1, 1) # Adjust these limits for null association

custom_scale_fill <- function(result_type) {
  if (result_type == "Association") {
    scale_fill_gradient2(
      name = "Bias (%)",
      low = "purple",   # lower end
      mid = "white",   # midpoint (0 bias)
      high = "red",     # higher end
      midpoint = 0,     # force mid color at 0
      limits = assoc_limits,
      breaks = seq(assoc_limits[1], assoc_limits[2], length.out = 5),
      guide = guide_colorbar(reverse = FALSE)
    )
  } else {
    scale_fill_gradient2(
      name = "Bias (absolute)",
      low = "lightblue",   # lower end
      mid = "white",   # midpoint (0 bias)
      high = "#5E82C9",     # higher end - lightish blue
      midpoint = 0,     # force mid color at 0
      limits = null_limits,  # Added limits for consistency
      breaks = seq(null_limits[1], null_limits[2], length.out = 5),
      guide = guide_colorbar(reverse = FALSE)
    )
  }
}


# Apply the custom scales manually and combine the plots
library(cowplot) # for combining plots

# Create the heatmap for bias_assoc with consistent color scale
p1 <- ggplot(bias_assoc, aes(x = bias_type, y = error_level, fill = bias_mean)) +
  geom_tile(color = "black", linewidth = 0.05)  +
  custom_scale_fill("Association") +
  geom_text(aes(
    label = sprintf("%.1f%%", pmin(pmax(bias_mean, assoc_limits[1]), assoc_limits[2])),
    color = ifelse(bias_mean > -50, "black", "white")  # Adjust this condition based on your scale
  ), size = 5) +  # Apply limits to the labels
  labs(title = "Percentage Bias by Scenario - Association",
       x = "Scenario", y = "Error Level") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 16, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 16, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 14, face = "plain", colour = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 14, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black",  angle = 30, hjust = 1),
        legend.title  = element_text(size = 14, face = "plain", colour = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  scale_color_identity()  # Ensure that the text color is respected

p1

# Create the heatmap for bias_null with consistent color scale
p2 <- ggplot(bias_null, aes(x = bias_type, y = error_level, fill = bias_mean)) +
  geom_tile(color = "black", linewidth = 0.05) +
  geom_text(aes(label = sprintf("%.2f", bias_mean)), size = 5, color = "black") +  # Removed % sign
  custom_scale_fill("Null Association") +
  labs(title = "Absolute Bias by Scenario - Null Association",  # Changed title
       x = "Scenario", y = "Error Level") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 16, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 16, face = "plain", colour = "black"),    
        axis.text.x  = element_text(size = 14, face = "plain", colour = "black", angle = 30, hjust = 1),
        axis.text.y  = element_text(size = 14, face = "plain", colour = "black"),
        legend.text  = element_text(size = 10, face = "plain", colour = "black",   angle = 30, hjust = 1),
        legend.title = element_text(size = 14, face = "plain", colour = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "top")

p2
# Combine the plots
pdf("ass_bias_perc_heatmap.pdf", width = 8, height = 6)
p1
dev.off()

pdf("null_bias_perc_heatmap.pdf", width = 8, height = 6)
p2
dev.off()
rm(list=ls())
