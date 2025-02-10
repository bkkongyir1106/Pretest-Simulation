library(nortest)  
library(tseries)  
library(moments)  

# Outliers
calculate_outliers <- function(samples) {
  qnt <- quantile(samples, probs = c(0.25, 0.75))
  H <- 1.5 * IQR(samples)
  outliers <- sum(samples < (qnt[1] - H) | samples > (qnt[2] + H))
  return(outliers)
}

# Function to perform the normality tests and calculations
determine_normality <- function(samples) {
  shapiro_result <- shapiro.test(samples)
  ad_result <- ad.test(samples)
  ks_result <- ks.test(samples, "pnorm", mean(samples), sd(samples))
  jb_result <- jarque.bera.test(samples)
  cvm_result <- nortest::cvm.test(samples)
  Outliers <- calculate_outliers(samples)
  

  sample_skewness <- skewness(samples)
  sample_kurtosis <- kurtosis(samples) - 3 
  skewness_check <- abs(sample_skewness) <= 0.3
  kurtosis_check <- abs(sample_kurtosis) <= 3
  
  # Check the p-values of the normality tests (threshold = 0.05)
  shapiro_check <- shapiro_result$p.value > 0.05
  ad_check <- ad_result$p.value > 0.05
  ks_check <- ks_result$p.value > 0.05
  jb_check <- jb_result$p.value > 0.05
  cvm_check <- cvm_result$p.value > 0.05
  number_outlier <- Outliers < 0.003 * length(samples)
  
  # Count the number of satisfied conditions
  total_conditions <- 8  
  satisfied_conditions <- sum(skewness_check, kurtosis_check, shapiro_check, ad_check, ks_check, jb_check, cvm_check, number_outlier)
  
  threshold <- 0.50 * total_conditions
  classification <- ifelse(satisfied_conditions < threshold, "Non-Normally Distributed", "Normally Distributed")
  
  # Return results as a list
  list(
    Shapiro_P_Value = shapiro_result$p.value,
    AD_P_Value = ad_result$p.value,
    KS_P_Value = ks_result$p.value,
    JB_P_Value = jb_result$p.value,
    CVM_P_Value = cvm_result$p.value,
    Skewness = sample_skewness,
    Kurtosis = sample_kurtosis,
    Number_outlier = Outliers,
    Classification = classification,
    Satisfied_Conditions = satisfied_conditions,
    Total_Conditions = total_conditions,
    Threshold = threshold
  )
}

# Example usage
set.seed(123)  # For reproducibility
example_samples <- rgamma(100, shape = 1, rate = 1) # Generate a random sample from a chi-squared distribution
result <- determine_normality(example_samples)
print(result)

