# Load necessary libraries
library(glmmTMB)
library(ggplot2)
library(pwr)
library(effectsize)
library(tidyverse)

# Load data
data <- read.table("CoancestryInput.txt", sep="\t", header=F)

# Run this to load inbreeding coefs without loci removal
inbreed_coef <- read.table("CoancestryOutput.txt", sep=" ", header=T)

# Add survivorship column
data["survivorship"] <- c(0,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1,1)

# Create trimmed data with inbreeding coefficient and survivorship
trimmed_data <- data[c("V1", "survivorship")]
trimmed_data["inbreed_coef"] <- inbreed_coef["TrioML"]
trimmed_data["group"] <- factor(rep(c("AF", "AG", "AL", "BI", "BU", "BW"), 
                             times=c(4, 5, 2, 5, 4, 6)))

#trimmed_data <- trimmed_data[-c(7,8,9,18,20),] #remove heterozygous dead birds
#trimmed_data <- trimmed_data[-c(5,6,7,8,9),] #remove ag group birds
#trimmed_data <- trimmed_data[-c(17,18,19,20),] #remove bu group birds
trimmed_data <- trimmed_data[-c(5,6,7,8,9,17,18,19,20),] #remove ag and bu group birds


# T-test
t.test(inbreed_coef ~ survivorship, data = trimmed_data,alternative = "greater")
d <- (cohens_d(inbreed_coef ~ survivorship, data = trimmed_data))
d_value <- as.numeric(d$Cohens_d)

g <- (hedges_g(inbreed_coef ~ survivorship, data = trimmed_data))
g_value <- as.numeric(g$Hedges_g)

# Find current power given our sample size
pwr.t.test(d = g_value,
           n = 12,
           sig.level = 0.05,
           type = "two.sample",
           alternative = "greater")

# Find sample size needed for 0.8 power
pwr.t.test(d = g_value,
           power = 0.8,
           sig.level = 0.05,
           type = "two.sample",
           alternative = "greater")


# Fit a logistic regression model using glmmTMB
model_glmmTMB <- glmmTMB(survivorship ~ inbreed_coef+(1|group), data = trimmed_data, family = binomial)

# Summary of the model
summary(model_glmmTMB)

# Get model summary
model_summary <- summary(model_glmmTMB)

# Extract p-values for the fixed effects
p_values <- model_summary$coefficients$cond[, "Pr(>|z|)"]  # p-values for fixed effects

# Get the p-value for inbreed_coef (for example)
p_value_inbreed_coef <- p_values["inbreed_coef"]

# Generate a sequence of inbreeding coefficients for predictions
inbreeding_seq <- seq(min(trimmed_data$inbreed_coef), max(trimmed_data$inbreed_coef), length.out = 100)
newdata <- data.frame(inbreed_coef = inbreeding_seq)
newdata$group <- factor(rep(levels(trimmed_data$group)[1], nrow(newdata)))

# Get predicted values and confidence intervals
predictions <- predict(model_glmmTMB, newdata = newdata, type = "link", se.fit = TRUE)

# Convert to probability scale (logit to probability)
newdata$fit <- plogis(predictions$fit)  # Predicted probabilities
newdata$lwr <- plogis(predictions$fit - 1.96 * predictions$se.fit)  # Lower CI
newdata$upr <- plogis(predictions$fit + 1.96 * predictions$se.fit)  # Upper CI

# Plot logistic regression curve with confidence interval
ggplot(trimmed_data, aes(x = inbreed_coef, y = survivorship, color = group)) +
  geom_point(size = 2) +  # Black points without jitter
  scale_y_continuous(breaks = seq(0, 1, by = 1)) +
  geom_line(data = newdata, aes(x = inbreed_coef, y = fit), color = "black", size = 1) +  # Black logistic curve
  geom_ribbon(data = newdata, aes(x = inbreed_coef, ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkgray", inherit.aes = FALSE) +  # Gray confidence interval shading
  labs(
    title = "",
    x = "Individual Inbreeding Coefficient",
    y = "Hatching Success (0 = unhatched, 1 = hatched)"
  ) +
  theme(
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16)   # Increase y-axis label size
  )
#  + annotate("text", x = 0.5, y = 0.9, label = paste("p=", round(p_value_inbreed_coef, 4)), size = 3, hjust = 0)  # Add p-value text


ggplot(data=trimmed_data, aes(group, inbreed_coef, colour=group)) + 
  geom_boxplot() + 
  labs(
    title = "",
    x = "Breeding Groups",
    y = "Individual Inbreeding Coefficient"
  )+
  theme(
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16)   # Increase y-axis label size
  )


with(trimmed_data, pairwise.t.test(inbreed_coef, group, p.adj="bonferroni"))

# Calculate group mean and standard deviation
summary_df <- trimmed_data %>%
  group_by(group) %>%
  summarize(
    mean_value = mean(inbreed_coef, na.rm = TRUE),
    sd_value = sd(inbreed_coef, na.rm = TRUE)
  )

# View the result
print(summary_df)



