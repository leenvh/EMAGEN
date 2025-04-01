
#### haplotype correlation analysis Alex 

# Load necessary libraries
library(dplyr)
library(psych)
library(officer)
library(flextable)
library(tidyr)

# reading the data
hap_dat <- read.csv("data/data_haplotype.csv")

# 
dat1 <- hap_dat %>%
  select(pvivaxtopfalciparumproportions, apiannualparasiteindex, AICNI, AICNI_ISAKAA, AICNI.ISGEAA, AIRNI, AIRNI_ISGEAA, CVIET, CVMNK_Wt, ISAKAA_Wt, ISGEAA, NFSND)

# To calculate the  Spearman correlation and p-values
cor_results <- corr.test(dat1, method = "spearman")

# Extract the correlation coefficients and p-values
cor_matrix <- cor_results$r
p_matrix <- cor_results$p

#   data frame for the results
variables <- rownames(cor_matrix)[-1:-2]  # to remove the first two rows (self-correlations)
results <- data.frame(
  Variable = rep(variables, 2),
  Metric = c(rep("Correlation_with_pvivaxtopfalciparumproportions", length(variables)),
             rep("Correlation_with_apiannualparasiteindex", length(variables))),
  Value = c(cor_matrix[-1:-2, "pvivaxtopfalciparumproportions"],
            cor_matrix[-1:-2, "apiannualparasiteindex"]),
  P_Value = c(p_matrix[-1:-2, "pvivaxtopfalciparumproportions"],
              p_matrix[-1:-2, "apiannualparasiteindex"])
)


# Export the results 
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(flextable(results))

print(doc, target = "Spearman_Correlation_Results_alex.docx")

### bootstrapping 
# Function to perform bootstrap for each variable with p vivax to p falciparum proportions and api annual parasite index
bootstrap_results <- function(data, var1, var2) {
  bootstrap_result <- boot(data = hap_dat, statistic = cor_func, R = 1000, var1 = var1, var2 = var2)
  
  # Calculate confidence intervals using the percentile method
  bootstrap_ci <- boot.ci(bootstrap_result, type = "perc")
  
  # Return the results
  return(list(
    estimate = bootstrap_result$t0, 
    ci_lower = bootstrap_ci$perc[4], 
    ci_upper = bootstrap_ci$perc[5]
  ))
}

# library
library(dplyr)
library(psych)
library(boot)

# Select the  the variables 
dat2 <- hap_dat %>%
  select(pvivaxtopfalciparumproportions, apiannualparasiteindex, AICNI, AICNI_ISAKAA, AICNI.ISGEAA, AIRNI, AIRNI_ISGEAA, CVIET, CVMNK_Wt, ISAKAA_Wt, ISGEAA, NFSND)

# Define a function to compute Spearman correlation
cor_func <- function(data, indices, var1, var2) {
  # Resample the data (with replacement)
  resampled_data <- data[indices, ]
  
  # Compute Spearman correlation
  cor_value <- cor(resampled_data[[var1]], resampled_data[[var2]], method = "spearman")
  return(cor_value)
}

bootstrap_results <- function(data, var1, var2) {
  # Apply bootstrap for Spearman correlation
  bootstrap_result <- boot(data = data, statistic = cor_func, R = 1000, var1 = var1, var2 = var2)
  
  bootstrap_ci <- boot.ci(bootstrap_result, type = "perc")
  
  # Return the results
  return(list(
    estimate = bootstrap_result$t0, 
    ci_lower = bootstrap_ci$perc[4], 
    ci_upper = bootstrap_ci$perc[5]
  ))
}

# List of other variables
other_variables <- colnames(dat2)[-c(1, 2)]  

# Initialize a results data frame to store the outcomes
results <- data.frame(
  Variable = character(),
  SpearmanCorrelation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each variable and apply the bootstrap function
for (var in other_variables) {
  # To calculate p correlation for p vivax to p falciparum proportions  
  result_pvivax <- bootstrap_results(dat2, "pvivaxtopfalciparumproportions", var)
  
  # To calculate  correlation for api annual parasite index  
  result_api <- bootstrap_results(dat2, "apiannualparasiteindex", var)
  
  # Append the results for p vivax to p falciparum proportions
  results <- rbind(results, data.frame(
    Variable = paste("pvivaxtopfalciparumproportions", "vs", var),
    SpearmanCorrelation = result_pvivax$estimate,
    CI_Lower = result_pvivax$ci_lower,
    CI_Upper = result_pvivax$ci_upper
  ))
  
  # Append the results for api annual parasite index
  results <- rbind(results, data.frame(
    Variable = paste("apiannualparasiteindex", "vs", var),
    SpearmanCorrelation = result_api$estimate,
    CI_Lower = result_api$ci_lower,
    CI_Upper = result_api$ci_upper
  ))
}


# Create a Word doc
doc <- read_docx()

doc <- doc %>%
  body_add_par("Bootstrap Spearman Correlation Results", style = "heading 1") %>%
  body_add_flextable(flextable(results))

# Save 
print(doc, target = "Bootstrap_Spearman_Correlation_Results_alex.docx")

