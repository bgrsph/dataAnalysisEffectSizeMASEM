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


# Compute Cohen's d and Hedges' g and convert them into point-biserial correlation (r_pb)
CBTStudiesCohensD <- ESDataCBTStudies %>%
  mutate(
    Cohen_d = mapply(getCohensD, mean_t, SD_t, n_t, mean_c, SD_c, n_c),  
    r_pb = mapply(getPointBiserialCorrelation, Cohen_d, n_t + n_c)  
  )

CBTStudiesHedgesG <- ESDataCBTStudies %>%
  mutate(
    Hedges_g = mapply(getHedgesG, mean_t, SD_t, n_t, mean_c, SD_c, n_c),
    r_pb = mapply(getPointBiserialCorrelation, Hedges_g, n_t + n_c)  
  )

# If a study includes more than one data per relationship, aggregate the effect sizes or correlation coefficients by averaging
CBTStudiesCohensD <- CBTStudiesCohensD %>%
  group_by(studyid, var1type, var2type) %>%
  summarise(
    r = ifelse(all(is.na(r)), NA, mean(r, na.rm = TRUE)),  # Keep NA if all are NA
    Cohen_d = ifelse(all(is.na(Cohen_d)), NA, mean(Cohen_d, na.rm = TRUE)),
    r_pb = ifelse(all(is.na(r_pb)), NA, mean(r_pb, na.rm = TRUE))
  ) %>%
  ungroup()

CBTStudiesHedgesG <- CBTStudiesHedgesG %>%
  group_by(studyid, var1type, var2type) %>%
  summarise(
    r = ifelse(all(is.na(r)), NA, mean(r, na.rm = TRUE)),  # Keep NA if all are NA
    Hedges_g = ifelse(all(is.na(Hedges_g)), NA, mean(Hedges_g, na.rm = TRUE)),
    r_pb = ifelse(all(is.na(r_pb)), NA, mean(r_pb, na.rm = TRUE))
  ) %>%
  ungroup()


# Preprocess the data sets to make them webMASEM-compatible

# Add Y_X, Y_M, and M_X columns
CBTStudiesCohensD <- CBTStudiesCohensD %>%
  mutate(
    Y_X = case_when(var1type == "1" & var2type == "O" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_),
    Y_M = case_when(var1type == "MA" & var2type == "O" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_),
    M_X = case_when(var1type == "1" & var2type == "MA" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_)
  )

CBTStudiesHedgesG <- CBTStudiesHedgesG %>%
  mutate(
    Y_X = case_when(var1type == "1" & var2type == "O" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_),
    Y_M = case_when(var1type == "MA" & var2type == "O" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_),
    M_X = case_when(var1type == "1" & var2type == "MA" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_)
  )


# Collapse data so that each studyid has one row with all available variable pairs
webMASEMCohensD <- CBTStudiesCohensD %>%
  group_by(studyid) %>%
  summarise(
    Y_X = if (all(is.na(Y_X))) NA else mean(Y_X, na.rm = TRUE),  # Keep NA if all values are NA
    Y_M = if (all(is.na(Y_M))) NA else mean(Y_M, na.rm = TRUE),
    M_X = if (all(is.na(M_X))) NA else mean(M_X, na.rm = TRUE)
  ) %>%
  ungroup()

webMASEMHedgesG <- CBTStudiesHedgesG %>%
  group_by(studyid) %>%
  summarise(
    Y_X = if (all(is.na(Y_X))) NA else mean(Y_X, na.rm = TRUE),  # Keep NA if all values are NA
    Y_M = if (all(is.na(Y_M))) NA else mean(Y_M, na.rm = TRUE),
    M_X = if (all(is.na(M_X))) NA else mean(M_X, na.rm = TRUE)
  ) %>%
  ungroup()

# Write the data sets into Excel files
write_xlsx(webMASEMCohensD, "~/Desktop/repo/dataAnalysisEffectSizeMASEM/output/webMASEM_CohensD.xlsx")
write_xlsx(webMASEMHedgesG, "~/Desktop/repo/dataAnalysisEffectSizeMASEM/output/webMASEM_HedgesG.xlsx")
