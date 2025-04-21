################################################################################
# Script Name: effectSizeFunctions.R
# Author: Bugra Sipahioglu
# Date: 12/03/2025
# Version: 0.1
# 
# Description: This script provides custom functions to compute effect sizes (Cohen's d, Hedges' g)
#' and transform them into point-biserial correlations using Aaron et al.'s (1998) formula.
#
# Functions
#' - `getCohensD()`: Computes Cohen's d for two independent groups.
#' - `getHedgesG()`: Computes Hedges' g by applying a correction factor to Cohen's d.
#' - `getPointBiserialCorrelation()`: Converts effect sizes (Cohen's d or Hedges' g) into point-biserial correlations.
#
# Dependencies:
#   - R version: 4.4.2 (2024-10-31) -- "Pile of Leaves"
#   - No external packages are required as the functions work with base R. 
#
# Contact:
#   - Email: b.sipahioglu@umail.leidenuniv.nl
#   - GitHub: https://github.com/bgrsph
################################################################################

#' Compute Cohen's d for Two Independent Groups
#'
#' This function calculates Cohen's d, a standardized effect size measure, 
#' for comparing two independent groups (treatment and control).
#'
#' @param mean_t Numeric. Mean of the treatment group.
#' @param SD_t Numeric. Standard deviation of the treatment group.
#' @param n_t Integer. Sample size of the treatment group.
#' @param mean_c Numeric. Mean of the control group.
#' @param SD_c Numeric. Standard deviation of the control group.
#' @param n_c Integer. Sample size of the control group.
#'
#' @return Numeric. Cohen's d value. Returns `NA` if any input is missing.
#' 
#' @examples
#' getCohensD(mean_t = 30, SD_t = 5, n_t = 50, mean_c = 28, SD_c = 6, n_c = 50)
#'
#' @export
getCohensD <- function(mean_t, SD_t, n_t, mean_c, SD_c, n_c) {

  if (any(is.na(c(mean_t, SD_t, n_t, mean_c, SD_c, n_c)))) {   # Check if any required values are missing
    print("\nInput is missing to calculate Cohen's d. Check the dataset.\n")
    return(NA) 
  }
  
  # Compute pooled standard deviation
  Spooled <- sqrt(((SD_t^2) * (n_t - 1) + (SD_c^2) * (n_c - 1)) / (n_t + n_c - 2))
  
  # Compute Cohen's d
  d <- (mean_t - mean_c) / Spooled
  
  return(d)
}



#' Compute Hedges' g for Two Independent Groups
#'
#' This function calculates Hedges' g, a standardized effect size measure, 
#' for comparing two independent groups (treatment and control).
#'
#' @param mean_t Numeric. Mean of the treatment group.
#' @param SD_t Numeric. Standard deviation of the treatment group.
#' @param n_t Integer. Sample size of the treatment group.
#' @param mean_c Numeric. Mean of the control group.
#' @param SD_c Numeric. Standard deviation of the control group.
#' @param n_c Integer. Sample size of the control group.
#'
#' @return Numeric. Hedges' g value. Returns `NA` if any input is missing.
#' 
#' @examples
#' getHedgesG(mean_t = 30, SD_t = 5, n_t = 50, mean_c = 28, SD_c = 6, n_c = 50)
#'
#' @export
getHedgesG <- function(mean_t, SD_t, n_t, mean_c, SD_c, n_c) {
 
  if (any(is.na(c(mean_t, SD_t, n_t, mean_c, SD_c, n_c)))) {  # Check if any required values are missing
    print("\nInput is missing to calculate Hedges' g. Check the dataset.\n")
    return(NA)
  }
  
  # Compute Cohen's d using the previously defined function
  d <- getCohensD(mean_t, SD_t, n_t, mean_c, SD_c, n_c)
  
  # Compute correction factor for small sample sizes
  correction_factor <- 1 - (3 / (4 * (n_t + n_c) - 9))
  
  # Compute Hedges' g
  g <- d * correction_factor
  
  return(g)
}


#' Compute Point-Biserial Correlation from Effect Size
#'
#' This function calculates the point-biserial correlation (r_pb) 
#' from a given effect size (Cohen's d or Hedges' g) using Aaron et al.'s (1998) formula.
#'
#' @param d Numeric. Effect size (e.g., Cohen's d or Hedges' g).
#' @param N Integer. Total sample size.
#'
#' @return Numeric. The computed point-biserial correlation. Returns `NA` if inputs are invalid.
#'
#' @examples
#' compute_rpb(d = 0.5, N = 100)
#'
#' @export
getPointBiserialCorrelation <- function(d, N) {
  # Check for valid inputs
  if (is.na(d) || is.na(N) || N <= 0) {
    print("\nInput is missing to calculate point-biserial correlation. Check the dataset.\n")
    return(NA)
  }
  
  # Compute point-biserial correlation using Aaron et al.'s (1998) formula
  r_pb <- d / sqrt(d^2 + 4 * (8/N)) 
  
  return(r_pb)
}