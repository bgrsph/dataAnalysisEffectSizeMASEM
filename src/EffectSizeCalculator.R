# This script involves custom functions of Cohen's d, Hedge's g, and transformations between these effect sizes and correlation coefficients. 


# Function to compute Cohen's d
#Note: treatment - control because it was mentioned in the paper or somwehere. 
getCohensD <- function(mean_t, SD_t, n_t, mean_c, SD_c, n_c) {
  # Check if any required values are missing
  if (any(is.na(c(mean_t, SD_t, n_t, mean_c, SD_c, n_c)))) {
    return(NA)  # Return NA if any input is missing
  }
  
  # Compute pooled standard deviation
  Spooled <- sqrt(((SD_t^2) * (n_t - 1) + (SD_c^2) * (n_c - 1)) / (n_t + n_c - 2))
  
  # Compute Cohen's d
  d <- (mean_t - mean_c) / Spooled
  
  return(d)
}


# Function to compute Hedges' g
getHedgesG <- function(mean_t, SD_t, n_t, mean_c, SD_c, n_c) {
  # Check if any required values are missing
  if (any(is.na(c(mean_t, SD_t, n_t, mean_c, SD_c, n_c)))) {
    return(NA)  # Return NA if any input is missing
  }
  
  # Compute Cohen's d using the previously defined function
  d <- getCohensD(mean_t, SD_t, n_t, mean_c, SD_c, n_c)
  
  # Compute correction factor for small sample sizes
  correction_factor <- 1 - (3 / (4 * (n_t + n_c) - 9))
  
  # Compute Hedges' g
  g <- d * correction_factor
  
  return(g)
}