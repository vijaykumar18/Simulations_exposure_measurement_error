
setwd("/scratch/user/u.vk230134/Downloads/DiffMeasError/CodeReview/exp/exp_002/Results")

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
library(glmmTMB)
library(stringr)



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
# List of dataset names
my.dta <- ls()

# Add 0 to error levels
err.lev <- c(0, 0.5, 1, 2)
results <- data.frame()


# Main loop
for (m in seq_along(my.dta)) {
  dta <- as.data.frame(get(my.dta[m]))
  index <- 1
  for (e in seq_along(err.lev)) {
    # Ensure index stays within bounds of the dataset's columns
    if (index + 1 <= ncol(dta)) {
      beta <- mean(dta[, index])
      vbar <- mean(dta[, index + 1])
      capb <- var(dta[, index])
      sterr <- sqrt(vbar + (nrow(dta) + 1) / nrow(dta) * capb)
      rr <- exp(beta)
      lci <- exp(beta - 1.96 * sterr)
      uci <- exp(beta + 1.96 * sterr)
      
      temp.res <- data.frame(beta, sterr, rr, lci, uci)
      # Add '+c' to the scenario name
      temp.res$scenario <- paste0(my.dta[m], "_", err.lev[e], "+c")
      
      results <- bind_rows(results, temp.res)
      
      rm(temp.res, beta, vbar, capb, sterr, rr, lci, uci)
    } else {
      # Handle case where index+1 is out of bounds (optional)
      warning("Skipping iteration for dataset: ", my.dta[m], " due to index bounds issue.")
    }
    
    index <- index + 2
  }
}

# Output the results
results
results_001 <- read_csv("/scratch/user/u.vk230134/Downloads/DiffMeasError/CodeReview/exp/exp_001/Results/results_001.csv")


# Filter rows where 'scenario' contains "diff_harm_lessErr"
df <- results %>% filter(str_detect(scenario, "diff_harm_lessErr_"))
df1<-results_001%>%filter(scenario=="diff_harm_lessErr_0")
df1<-rbind(df1,df)
df <- results %>% filter(str_detect(scenario, "diff_harm_lessErr50_"))
df2<-results_001%>%filter(scenario=="diff_harm_lessErr50_0")
df2<-rbind(df2,df)
df <- results %>% filter(str_detect(scenario, "diff_harm_moreErr_"))
df3<-results_001%>%filter(scenario=="diff_harm_moreErr_0")
df3<-rbind(df3,df)
df <- results %>% filter(str_detect(scenario, "diff_harm_moreErr50_"))
df4<-results_001%>%filter(scenario=="diff_harm_moreErr50_0")
df4<-rbind(df4,df)
df <- results %>% filter(str_detect(scenario, "diff_null_lessErr_"))
df5<-results_001%>%filter(scenario=="diff_null_lessErr_0")
df5<-rbind(df5,df)
df <- results %>% filter(str_detect(scenario, "diff_null_lessErr50_"))
df6<-results_001%>%filter(scenario=="diff_null_lessErr50_0")
df6<-rbind(df6,df)

df <- results %>% filter(str_detect(scenario, "diff_null_moreErr_"))
df7<-results_001%>%filter(scenario=="diff_null_moreErr_0")
df7<-rbind(df7,df)
df <- results %>% filter(str_detect(scenario, "diff_null_moreErr50_"))
df8<-results_001%>%filter(scenario=="diff_null_moreErr50_0")
df8<-rbind(df8,df)

df <- results %>% filter(str_detect(scenario, "nondiff_harm_"))
df9<-results_001%>%filter(scenario=="nondiff_harm_0")
df9<-rbind(df9,df)
df <- results %>% filter(str_detect(scenario, "nondiff_null_"))
df10<-results_001%>%filter(scenario=="nondiff_null_0")
df10<-rbind(df10,df)


results<- rbind(df1,df2,df3,df4,df5, df6, df7,df8,df9, df10)

results

# Calculating percent bias for association
results$perc.bias.assoc <- ifelse(results$scenario %in% results$scenario[grep("harm", results$scenario)],
                                  (results$beta - log(1.12))/log(1.12)*100, NA)

# Calculating percent bias for null associations
results$perc.bias.null <- ifelse(results$scenario %in% results$scenario[grep("null", results$scenario)],
                                 results$beta * 100, NA)

# Separate data frames for association and null association
results_assoc <- results[!is.na(results$perc.bias.assoc), c("scenario", "perc.bias.assoc")]
results_null  <- results[!is.na(results$perc.bias.null),  c("scenario", "perc.bias.null")]


# Save percent bias with association to a text file
sink("perc_bias_assoc.txt")
print(results_assoc)
sink()
# Save the data frame as a CSV file
write.csv(results_assoc, "perc_bias_assoc.csv", row.names = FALSE)


# Save percent bias for null associations to a different text file
sink("perc_bias_null.txt")
print(results_null)
sink()
write.csv(results_null, "perc_bias_null.csv", row.names = FALSE)

# Filter the data for null scenarios
null.res <- results[grep("null", results$scenario),]

# Define symbols for plotting based on scenario
null.res$sym <- ifelse(null.res$scenario %in% null.res$scenario[grep("nondiff", null.res$scenario)], 15, 19)

# Define colors for plotting based on scenario
null.res$col <- ifelse(null.res$scenario %in% null.res$scenario[grep("less", null.res$scenario)], "blue", 
                       ifelse(null.res$scenario %in% null.res$scenario[grep("more", null.res$scenario)], "cadetblue", 
                              "black"))

# Define line types based on scenario
null.res$lty <- ifelse(null.res$scenario %in% null.res$scenario[grep("Err50", null.res$scenario)], 2, 1)

# Define line widths based on scenario
null.res$lwd <- ifelse(null.res$scenario %in% null.res$scenario[grep("Err50", null.res$scenario)], 1.5, 1)

# Clean the 'err' column to remove unwanted parts before extracting error size
null.res$err <- gsub("diff_null_lessErr|diff_null_moreErr|nondiff_null", "", null.res$scenario)  # Remove prefixes like "diff_null_lessErr"
null.res$err <- gsub("50_", "", null.res$err)  # Remove '_50' from error size if it exists
null.res$err <- gsub("_", "", null.res$err)  # Remove underscore if it exists




pdf("null_results.pdf", width = 6, height = 6)  # Fixed dimensions

# Set plot margins (bottom, left, top, right)
par(mar = c(5, 5, 4, 2) + 0.1)
# Adjust margins to give more space for x-axis labels
#par(mar = c(10, 4, 2, 2))   # Increase bottom margin

# Set up the plot with custom axis
plot(
  c(1, dim(null.res)[1]),
  c(0.96, 1.04),
  pch = '', las = 1,
  xlab = "",
  ylab = "RR & 95% CI",
  cex.lab = 1.5, cex.axis = 1.25, xaxt = 'n'
)

# Add horizontal reference line at RR=1
abline(h = 1, col = "darkgray")

# Customize x-axis with rotated labels
axis(
  side = 1,
  at = 1:dim(null.res)[1],
  labels = null.res$err,
  cex.axis = 1.2,
  las = 2
)

# Reposition the x-axis label below the rotated tick labels
mtext("Error Level", side = 1, line = 3.8, cex = 1.5)

# Plot points and error bars
points(1:dim(null.res)[1], null.res$rr, pch = null.res$sym, col = null.res$col, cex = 1.5)
segments(
  1:dim(null.res)[1], null.res$lci,
  1:dim(null.res)[1], null.res$uci,
  col = null.res$col, cex = 1.5, lty = null.res$lty, lwd = null.res$lwd
)

# Add plot title
#title("Null Scenarios")


# Add legends
legend("topright", pch = c(19, 19, 15), col = c("blue", "cadetblue", "black"), 
       legend = c("Diff Less", "Diff More", "Nondiff"), bty = "n")
legend("bottomright", lty = c(1, 2), 
       legend = c("High-Risk", "High-Risk + 50% Rest"), bty = 'n')

# Close the PDF device
dev.off()


## Harmful Results

harm.res <- results[grep("harm", results$scenario),]
harm.res$sym <- ifelse(harm.res$scenario %in% harm.res$scenario[grep("nondiff", harm.res$scenario)], 15, 19)
harm.res$col <- ifelse(harm.res$scenario %in% harm.res$scenario[grep("less", harm.res$scenario)], "blue", 
                       ifelse(harm.res$scenario %in% harm.res$scenario[grep("more", harm.res$scenario)], "cadetblue",
                              "black"))
harm.res$lty <- ifelse(harm.res$scenario %in% harm.res$scenario[grep("Err50", harm.res$scenario)], 2, 1)
harm.res$lwd <- ifelse(harm.res$scenario %in% harm.res$scenario[grep("Err50", harm.res$scenario)], 1.5, 1)
harm.res$err <- rep(c("0", "0.5", "1", "2"), length.out = nrow(harm.res))

# Clean the 'err' column to remove unwanted parts before extracting error size
harm.res$err <- gsub("diff_harm_lessErr|diff_harm_moreErr|nondiff_harm", "", harm.res$scenario)  # Remove prefixes like "diff_null_lessErr"
harm.res$err <- gsub("50_", "", harm.res$err)  # Remove '_50' from error size if it exists
harm.res$err <- gsub("_", "", harm.res$err)  # Remove underscore if it exists


# Ensure 'err' is a factor and contains the correct levels
harm.res$err <- as.character(harm.res$err)

# Create the plot with custom x-axis
pdf("harm_results.pdf", width = 6, height = 6)  # Fixed dimensions

# Set plot margins (bottom, left, top, right)
par(mar = c(5, 5, 4, 2) + 0.1)


# Plot with a large range for Y to ensure no clipping of segments
plot(c(1, dim(harm.res)[1]), c(0.94,1.30), pch = '', las = 1,
     xlab = "", ylab = "RR & 95% CI", cex.lab = 1.5, cex.axis = 1.25, xaxt = "n")

# Add horizontal reference lines at RR=1 and RR=1.115
abline(h = 1, col = "darkgray")
abline(h = 1.115, col = "hotpink2", lty = 4)

# Customize x-axis with rotated labels
axis(side = 1, at = 1:length(harm.res$err), labels = as.character(harm.res$err), cex.axis = 1.2, las = 2)

# Reposition the x-axis label below the rotated tick labels
mtext("Error Level", side = 1, line = 3.8, cex = 1.5)

# Plot points and segments
points(1:length(harm.res$err), harm.res$rr, pch = harm.res$sym, col = harm.res$col, cex = 1.5)
segments(1:length(harm.res$err), harm.res$lci, 1:length(harm.res$err), harm.res$uci, col = harm.res$col, cex = 1.5, 
         lty = harm.res$lty, lwd = harm.res$lwd)

# Add title and legends
#title("Association Scenarios")
legend("topright", pch = c(19, 19, 15), col = c("blue", "cadetblue", "Black"), legend = c("Diff Less", "Diff More", "Nondiff"), bty = "n")
legend("bottomright", lty = c(1, 2), legend = c("High-Risk", "High-Risk + 50% Rest"), bty = 'n')

# Close the PDF device
dev.off()




#############
# Add a column to indicate the result type
# Read tabular data into R
bias <- read.csv("perc_bias_assoc.txt", header = FALSE)

# Read the CSV file into a data frame
bias <- read.csv("perc_bias_assoc.csv")

bias_null<-read.csv("perc_bias_null.csv")
bias_assoc <- bias %>%
  mutate(result_type = "Association")%>%rename(perc.bias=perc.bias.assoc)

bias_null <- bias_null %>%
  mutate(result_type = "Null Association")%>%rename(perc.bias=perc.bias.null)

# Combine the two datasets
combined_bias <- bind_rows(bias_assoc, bias_null)

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
bias_assoc
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
assoc_limits <- c(-100, 40) # Adjust these limits for association
null_limits <- c(-0.24, 0.24)      # Adjust these limits for null association

# Define custom scale based on the type
custom_scale_fill <- function(result_type) {
  if (result_type == "Association") {
    scale_fill_distiller(
      palette = "YlOrRd",
      direction = -1,
      name = "Bias (%)",
      limits = assoc_limits,
      breaks = seq(assoc_limits[1], assoc_limits[2], length.out = 5),
      guide = guide_colorbar(reverse = TRUE)
    )
  } else {
    scale_fill_distiller(
      palette = "Blues",
      direction = -1,
      name = "Bias (%)",
      limits = null_limits,
      breaks = seq(null_limits[1], null_limits[2], length.out = 5),
      guide = guide_colorbar(reverse = TRUE)
    )
  }
}


# Apply the custom scales manually and combine the plots
library(cowplot) # for combining plots

# Create the heatmap for bias_assoc with consistent color scale
p1 <- ggplot(bias_assoc, aes(x = bias_type, y = error_level, fill = perc.bias)) +
  geom_tile(color = "white") +
  custom_scale_fill("Association") +
  geom_text(aes(
    label = sprintf("%.1f%%", pmin(pmax(perc.bias, assoc_limits[1]), assoc_limits[2])),
    color = ifelse(perc.bias > -50, "black", "white")  # Adjust this condition based on your scale
  ), size = 5) +  # Apply limits to the labels
  labs(title = "Percentage Bias by Scenario - Association",
       x = "Scenario", y = "Error Level") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 16, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 16, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 14, face = "plain", colour = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 14, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        legend.title  = element_text(size = 14, face = "plain", colour = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  scale_color_identity()  # Ensure that the text color is respected

p1

# Create the heatmap for bias_null with consistent color scale
p2 <- ggplot(bias_null, aes(x = bias_type, y = error_level, fill = perc.bias)) +
  geom_tile(color = "white") +
  custom_scale_fill("Null Association") +
  geom_text(aes(label = sprintf("%.1f%%", pmin(pmax(perc.bias, null_limits[1]), null_limits[2]))), 
            size = 5) +  # Apply limits to the labels
  labs(title = "Percentage Bias by Scenario - Null Association",
       x = "Scenario", y = "Error Level") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 16, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 16, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 14, face = "plain", colour = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 14, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        legend.title  = element_text(size = 14, face = "plain", colour = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "top")
p2

p2
# Combine the plots
pdf("ass_bias_perc_heatmap_044.pdf", width = 8, height = 6)
p1
dev.off()

pdf("null_bias_perc_heatmap_044.pdf", width = 8, height = 6)
p2
dev.off()
rm(list=ls())

