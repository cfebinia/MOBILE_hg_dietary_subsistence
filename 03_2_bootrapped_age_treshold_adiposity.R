rm(list=ls())


library(lmerTest)
library(boot)
library(tidyverse)

###############
### Obese x Adiposity Age treshold prediction
################
# load data and iput
df <- read.delim("out/anthrop_data_median_imputed_processed.tsv")
final_model <- readRDS("out/adiposity_lmer_final.rds")

# Extract coefficients
coeffs <- fixef(final_model)
b0 <- coeffs["(Intercept)"]
b_bmi <- coeffs["bmi"]
b_age <- coeffs["age"]
b_sexM <- coeffs["sexM"]

# Calculate point estimates for Age
# Formula: Age = (0.51 - b0 - (b_bmi * 23) - (b_sex * Sex)) / b_age
age_female_est <- (0.51 - b0 - (b_bmi * 23)) / b_age
age_male_est <- (0.51 - b0 - (b_bmi * 23) - b_sexM) / b_age

# Parametric bootstrapping for 95% Confidence Intervals
boot_calc <- function(model) {
  b <- fixef(model)
  # Calculate age for Female (sexM = 0) and Male (sexM = 1)
  a_f <- (0.51 - b[1] - (b[2] * 23)) / b[3]
  a_m <- (0.51 - b[1] - (b[2] * 23) - b[4]) / b[3]
  return(c(Female = a_f, Male = a_m))
}

set.seed(123)
boot_results <- bootMer(final_model, boot_calc, nsim = 1000) # simulate parametric bootstrap
intervals <- confint(boot_results, type = "perc")

# Combine results into a formatted table
results_table <- data.frame(
  Sex = c("Female", "Male"),
  Predicted_Age = c(age_female_est, age_male_est),
  Lower_CI = intervals[, 1],
  Upper_CI = intervals[, 2]
)

print(results_table)

#############################
## Usig BOOT (original data)
#############################
library(boot)

# 1. Define the function for the boot operator
# 'data' is the original dataframe, 'indices' is what boot uses to resample
boot_age_logic <- function(data, indices) {
  d <- data[indices, ]
  
  # Refit the model on resampled data
  # Use try() to handle potential convergence issues during resampling
  fit <- try(lmer(WHtR ~ bmi + age + sex + (1 | population), data = d), silent = TRUE)
  
  if (inherits(fit, "try-error")) return(c(NA, NA))
  
  # Extract coefficients
  b <- fixef(fit)
  
  # Calculate Age for Females (sexM = 0) and Males (sexM = 1)
  # Target WHtR = 0.51, Target BMI = 23
  age_f <- (0.51 - b["(Intercept)"] - (b["bmi"] * 23)) / b["age"]
  age_m <- (0.51 - b["(Intercept)"] - (b["bmi"] * 23) - b["sexM"]) / b["age"]
  
  return(c(Female = age_f, Male = age_m))
}

# 2. Run the bootstrap
set.seed(123)
boot_results <- boot(data = df, statistic = boot_age_logic, R = 1000)

# 3. View results and Confidence Intervals
print(boot_results)

# Get CI
combined_ci <- do.call(rbind, lapply(1:2, function(i) {
  res <- boot.ci(boot_results, type = "perc", index = i)
  
  data.frame(
    Sex = ifelse(i == 1, "Female", "Male"),
    Lower_CI = res$percent[4],
    Upper_CI = res$percent[5]
  )
}))

final_summary <- data.frame(
  Sex = c("Female", "Male"),
  Predicted_Age = boot_results$t0
) %>%
  cbind(combined_ci[, c("Lower_CI", "Upper_CI")])

print(final_summary)

# writeout
write.table(final_summary, "out/adiposity_age_tresholds_boot_nonParam.tsv", quote = F, sep = "\t", row.names = F)
write.table(results_table, "out/adiposity_age_tresholds_boot_Param.tsv", quote = F, sep = "\t", row.names = F)
