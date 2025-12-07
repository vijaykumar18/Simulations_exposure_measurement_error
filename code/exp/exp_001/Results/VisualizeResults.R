setwd("~/Downloads/CodeReview/CodeReview_Heather/exp/exp_001/Results")
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

# List of dataset names
my.dta <- ls()

# Add 0 to error levels
err.lev <- c(0, 0.5, 1, 2)
results <- data.frame()
bias_results <- data.frame()

# Main loop - WITH COVERAGE
for (m in seq_along(my.dta)) {
  dta <- as.data.frame(get(my.dta[m]))
  
  for (e in seq_along(err.lev)) {
    # Calculate column positions directly for each error level
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
      temp.res$scenario <- paste(my.dta[m], err.lev[e], sep = "_")
      
      results <- bind_rows(results, temp.res)
      
      # Calculate bias AND COVERAGE
      beta_values <- dta[, b_col]
      
      # Determine result type and calculate bias
      if (grepl("harm", my.dta[m])) {
        result_type <- "Association"
        bias_values <- (beta_values - log(1.115)) / log(1.115) * 100
        true_beta <- log(1.115)
      } else {
        result_type <- "Null Association"
        bias_values <- beta_values 
        true_beta <- 0
      }
      
      # CALCULATE COVERAGE
      ci_lower <- beta_values - 1.96 * dta[, se_col]
      ci_upper <- beta_values + 1.96 * dta[, se_col]
      coverage <- mean(ci_lower <= true_beta & true_beta <= ci_upper) * 100
      
      # Calculate percentiles
      bias_mean <- mean(bias_values)
      bias_2.5 <- quantile(bias_values, 0.025)
      bias_97.5 <- quantile(bias_values, 0.975)
      
      temp.bias <- data.frame(
        scenario = paste(my.dta[m], err.lev[e], sep = "_"),
        bias_mean = bias_mean,
        bias_2.5 = bias_2.5,
        bias_97.5 = bias_97.5,
        coverage = coverage,  # ADD COVERAGE
        result_type = result_type
      )
      
      bias_results <- bind_rows(bias_results, temp.bias)
      
      rm(temp.res, beta, vbar, capb, sterr, rr, lci, uci)
    } else {
      warning("Columns not found for dataset: ", my.dta[m], " error level: ", err.lev[e])
    }
  }
}

# Create coverage_results
coverage_results <- bias_results %>% 
  dplyr::select(scenario, coverage, result_type)

# Output the results
summary(results)
summary(bias_results)
summary(coverage_results)

write_csv(results, "results.csv")
write_csv(bias_results, "bias_results.csv")
write_csv(coverage_results, "coverage_results.csv")
#################################
## Null Results

null.res <- results[grep("null", results$scenario),]
null.res$sym <- ifelse(null.res$scenario %in% null.res$scenario[grep("nondiff", null.res$scenario)], 15, 19)
null.res$col <- ifelse(null.res$scenario %in% null.res$scenario[grep("less", null.res$scenario)], "blue", 
                       ifelse(null.res$scenario %in% null.res$scenario[grep("more", null.res$scenario)], "cadetblue",
                              "black"))
null.res$lty <- ifelse(null.res$scenario %in% null.res$scenario[grep("Err50", null.res$scenario)], 2, 1)
null.res$lwd <- ifelse(null.res$scenario %in% null.res$scenario[grep("Err50", null.res$scenario)], 1.5, 1)
null.res$err <- rep(c("0", "0.5", "1", "2"), length.out = nrow(null.res))

# Ensure the 'err' vector is correctly assigned
summary(null.res)

pdf("null_results.pdf", width = 6, height = 6)
par(mar = c(5, 5, 4, 2) + 0.01)

plot(
  c(1, nrow(null.res)), c(0.99,1.01),
  pch = '', las = 1,
  xlab = "",
  ylab = "RR & 95% CI",
  cex.lab = 1.5, cex.axis = 1.25, xaxt = 'n',
log="y")

abline(h = 1, col = "darkgray")

axis(
  side = 1,
  at = 1:dim(null.res)[1],
  labels = null.res$err,
  cex.axis = 1.2,
  las = 2
)

mtext("Error Level", side = 1, line = 3.5, cex = 1.5)

# Keep original colors, but use diamond shape for Err50 scenarios
points(1:dim(null.res)[1], null.res$rr, 
       pch = ifelse(grepl("Err50", null.res$scenario), 23, null.res$sym),  # Diamond for Err50, original shape for others
       col = null.res$col, 
       cex = 1.5)

segments(
  1:dim(null.res)[1], null.res$lci,
  1:dim(null.res)[1], null.res$uci,
  col = null.res$col, 
  lty = null.res$lty, 
  lwd = null.res$lwd
)

# Update legend to show diamond for Err50 scenarios
legend("topright", 
       pch = c(19, 23, 19, 23, 15),  # Circle, Diamond, Circle, Diamond, Square
       col = c("blue", "blue", "cadetblue", "cadetblue", "black"), 
       legend = c("Diff Less", "Diff Less + 50%", "Diff More", "Diff More + 50%", "Nondiff"), 
       bty = "n")

legend("bottomleft", lty = c(1, 2), 
       legend = c("High-Risk", "High-Risk + 50% Rest"), 
       bty = 'n')

dev.off()

## Harmful Results

harm.res     <- results[grep("harm", results$scenario),]
harm.res$sym <- ifelse(harm.res$scenario %in% harm.res$scenario[grep("nondiff", harm.res$scenario)], 15, 19)
harm.res$col <- ifelse(harm.res$scenario %in% harm.res$scenario[grep("less", harm.res$scenario)], "blue", 
                       ifelse(harm.res$scenario %in% harm.res$scenario[grep("more", harm.res$scenario)], "cadetblue",
                              "black"))
harm.res$lty <- ifelse(harm.res$scenario %in% harm.res$scenario[grep("Err50", harm.res$scenario)], 2, 1)
harm.res$lwd <- ifelse(harm.res$scenario %in% harm.res$scenario[grep("Err50", harm.res$scenario)], 1.5, 1)
harm.res$err <- rep(c("0", "0.5", "1", "2"), length.out = nrow(harm.res))

harm.res$err <- as.character(harm.res$err)
summary(harm.res)

pdf("harm_results.pdf", width = 6, height = 6)
par(mar = c(5, 5, 4, 2) + 0.01)


plot(c(1, dim(harm.res)[1]), c(0.99,1.124), pch = '', las = 1,
     xlab = "", ylab = "RR & 95% CI", cex.lab = 1.5, cex.axis = 1.25, xaxt = "n",
     log = "y")  # This enables log scale on y-axis

abline(h = 1, col = "darkgray")
abline(h = 1.115, col = "hotpink2", lty = 4)

axis(side = 1, at = 1:length(harm.res$err), labels = as.character(harm.res$err), cex.axis = 1.2, las = 2)

mtext("Error Level", side = 1, line = 3.5, cex = 1.5)

points(1:dim(harm.res)[1], harm.res$rr, 
       pch = ifelse(grepl("Err50", harm.res$scenario), 23, harm.res$sym),  # Diamond for Err50
       col = harm.res$col, 
       cex = 1.5)

segments(
  1:dim(harm.res)[1], harm.res$lci,
  1:dim(harm.res)[1], harm.res$uci,
  col = harm.res$col, 
  lty = harm.res$lty, 
  lwd = harm.res$lwd
)

# Update legend
legend("topright", 
       pch = c(19, 23, 19, 23, 15), 
       col = c("blue", "blue", "cadetblue", "cadetblue", "black"), 
       legend = c("Diff Less", "Diff Less + 50%", "Diff More", "Diff More + 50%", "Nondiff"), 
       bty = "n")
#legend("topright", pch = c(19, 19, 15), col = c("blue", "cadetblue", "Black"), legend = c("Diff Less", "Diff More", "Nondiff"), bty = "n")
legend("bottomleft", lty = c(1, 2), legend = c("High-Risk", "High-Risk + 50% Rest"), bty = 'n')

dev.off()




rm(list=ls())

#############
# Add a column to indicate the result type

bias_results <- read_csv("bias_results.csv")

bias_assoc<-bias_results%>%dplyr::select(scenario,bias_mean,result_type)%>%filter(result_type=="Association")
bias_null<-bias_results%>%dplyr::select(scenario,bias_mean,result_type)%>%filter(result_type=="Null Association")

# Combine the two datasets
#combined_bias <- bind_rows(bias_assoc, bias_null)

# Modify bias_type levels and prepare the data
bias_assoc <- bias_assoc %>%
  mutate(bias_type = factor(gsub("_\\d.*", "", scenario),
                            levels = c("nondiff_harm", 
                                       "diff_harm_lessErr", 
                                       "diff_harm_lessErr50",
                                       "diff_harm_moreErr",  
                                       "diff_harm_moreErr50")),
         error_level = factor(gsub(".*_", "", scenario), 
                              levels = c("0","0.5", "1", "2")))


bias_null <- bias_null %>%
  mutate(bias_type = factor(gsub("_\\d.*", "", scenario),
                            levels = c("nondiff_null", 
                                       "diff_null_lessErr", 
                                       "diff_null_lessErr50",
                                       "diff_null_moreErr",  
                                       "diff_null_moreErr50")),
         error_level = factor(gsub(".*_", "", scenario), 
                              levels = c("0","0.5", "1", "2")))




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
      name = "Bias (%)",
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

