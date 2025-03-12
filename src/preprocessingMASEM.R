################################################################################
# Script Name: preprocessingMASEM.R
# Author: Bugra Sipahioglu
# Date: 12/03/2025
# Version: 0.1
# 
# Description: This script imports effect size data, identifies relevant CBT 
#              studies, calculates Cohen's d and Hedges' g, converts them to point-biserial correlation coefficients
#              and preprocesses data for Meta-Analytic Structural Equation Modeling (MASEM).
#
# Datasets Used:
#   - ESdata.csv
#   - unimeta.CSV
#   - MeasureCharacteristics.xlsx
#   - StudyCharacteristics.xlsx
#
# Dependencies:
#   - R version: 4.4.2 (2024-10-31) -- "Pile of Leaves"
#   - Required Libraries: dplyr, readxl, writexl, tidyr, gt (and any others used)
#   - Custom Functions: effectSizeFunctions.R
#
# Outputs:
#   - 2 WebMASEM-compatible Excel files that include all the correlation coefficients of the variables of the mediation model
#
# Usage:
#   - Make sure that all required datasets and the function file are in the correct directory.
#   - Execute this script in RStudio or run source("preprocessingMASEM.R").
#
# Notes:
#   - Later on, this script will focus on preprocessing the 2 Excel files following WebMASEM Shiny App Tutorial
#
# Contact:
#   - Email: b.sipahioglu@umail.leidenuniv.nl
#   - GitHub: https://github.com/bgrsph
################################################################################

# Load libraries
library(dplyr)
library(readxl)
library(writexl)
library(tidyr)
library(gt)
source("~/Desktop/repo/dataAnalysisEffectSizeMASEM/src/effectSizeFunctions.R")

# Import datasets
ESdata<-read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/ESdata.csv")
unimeta <- read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Univariate/unimeta.CSV")
measureCharacteristicsData<-read_xlsx("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/MeasureCharacteristics.xlsx")
studyCharacteristicsData <- read_xlsx("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/StudyCharacteristics.xlsx")

# Extract the CBT study IDs from the unimeta.csv
cbtStudyIDs <- unimeta %>%
  filter(IPT == 0) %>%  # Select only CBT trials (Study 1)
  pull(studyid) %>%  # Extract study IDs
  unique()  # Keep unique study IDs

# Extract the relevant effect size data of CBT studies from ESdata.csv
# - Mediation paths:
#   - No time-point CBT (0) → Post-treatment Negative Cognition (3)
#   - Post-treatment Negative Cognition (3) → Post-treatment Depression (3)
#   - No time-point CBT (0) → Post-treatments Depression (3)
ESDataCBTStudies <- ESdata %>%
  filter(
    studyid %in% cbtStudyIDs &
      (
        (var1type == "1" & var2type == "MA" & Var1time == 0 & Var2time == 3) |    # Path a: CBT Treatment (1) → Negative Cognition (MA)
          (var1type == "MA" & var2type == "O" & Var1time == 3 & Var2time == 3) |  # Path b: Negative Cognition (MA) → Depression Severity (O)
          (var1type == "1" & var2type == "O" & Var1time == 0 & Var2time == 3)     # Path c': CBT Treatment (1) → Depression Severity
      )
  ) %>%
  select(studyid, var1type, Var1time, var2type, Var2time, # Include only parameters required to compute SMDs
         mean_t, SD_t, n_t, mean_c, SD_c, n_c, r)


# Compute Cohen's d and point-biserial correlation (r_pb) for each row
ESDataCBTStudiesCohensD <- ESDataCBTStudies %>%
  mutate(
    Cohen_d = mapply(getCohensD, mean_t, SD_t, n_t, mean_c, SD_c, n_c),  
    r_pb = mapply(getPointBiserialCorrelation, Cohen_d, n_t + n_c)  #
  )

# Compute Hedges' g and point-biserial correlation (r_pb) for each row
ESDataCBTStudiesHedgesG <- ESDataCBTStudies %>%
  mutate(
    Hedges_g = mapply(getHedgesG, mean_t, SD_t, n_t, mean_c, SD_c, n_c),
    r_pb = mapply(getPointBiserialCorrelation, Hedges_g, n_t + n_c)  
  )

# Next: preprocess datasets for WebMASEM



