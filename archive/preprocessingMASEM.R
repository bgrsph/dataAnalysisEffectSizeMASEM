################################################################################
# Script Name: preprocessingMASEM.R
# Author: Bugra Sipahioglu
# Date: 21/04/2025
# Version: 1.0
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
#   - 1 WebMASEM-compatible Excel files that include all the correlation coefficients of the variables of the mediation model
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
ESdata <- read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/ESdata.csv")
unimeta <- read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Univariate/unimeta.CSV")

# Extract the CBT study IDs from the unimeta.csv
cbtStudyIDs <- unimeta %>%
  filter(IPT == 0) %>%  # Select only CBT trials (Study 1)
  pull(studyid) %>%  # Extract study IDs
  unique()  # Keep unique study IDs

# Extract the relevant effect size data of CBT studies from ESdata.csv
# - Mediation paths (as shown in Table 1): 
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
  select(studyid, var1type, Var1time, var2type, Var2time, mean_t, SD_t, n_t, mean_c, SD_c, n_c, r, NumOppDir) # Only include parameters required to compute Cohen's d & Hedges' g

# Compute Cohen's d & Hedges' g
# Reverse them and r if one of the variables is reverse-coded compared to the theoretical expectation
# TODO: should I reverse r? try it out. 
CBTStudiesCohensD <- ESDataCBTStudies %>%
  mutate(
    Cohen_d = mapply(getCohensD, mean_t, SD_t, n_t, mean_c, SD_c, n_c),
    Cohen_d = ifelse(NumOppDir == 1, -Cohen_d, Cohen_d),
    r = ifelse(NumOppDir == 1 , -r, r),
    N_total = n_t + n_c
  )  %>%
  select(studyid, var1type, var2type, Cohen_d, r, N_total) 

CBTStudiesHedgesG <- ESDataCBTStudies %>%
  mutate(
    Hedges_g = mapply(getHedgesG, mean_t, SD_t, n_t, mean_c, SD_c, n_c),
    Hedges_g = ifelse(NumOppDir == 1, -Hedges_g, Hedges_g),
    r = ifelse(NumOppDir == 1 , -r, r),
    N_total = n_t + n_c
  ) %>%
  select(studyid, var1type, var2type, Hedges_g, r, N_total)
# 
# 
# # Aggregate effect sizes, correlation coefficients and sample sizes
# # Sample sizes are aggregated here for the first time:
# #   - This gives us one average sample size per mediation path within each study,
# #     in case multiple effect size estimates (e.g., different measures) were reported for the same path.
# 
# CBTStudiesCohensDAggregated <- CBTStudiesCohensD %>%
#   group_by(studyid, var1type, var2type) %>%
#   summarise(
#     Cohen_d = if (all(is.na(Cohen_d))) NA else mean(Cohen_d, na.rm = TRUE),
#     r = if (all(is.na(r))) NA else mean(r, na.rm = TRUE),
#     N_total = if (all(is.na(N_total))) NA else mean(N_total, na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# CBTStudiesHedgesGAggregated <- CBTStudiesHedgesG %>%
#   group_by(studyid, var1type, var2type) %>%
#   summarise(
#     Hedges_g = if (all(is.na(Hedges_g))) NA else mean(Hedges_g, na.rm = TRUE),
#     r = if (all(is.na(r))) NA else mean(r, na.rm = TRUE),
#     N_total = if (all(is.na(N_total))) NA else mean(N_total, na.rm = TRUE)
#   ) %>%
#   ungroup()


# Group-based reversal logic (mimicking meta-analysis rules)
#   Apply reverse rules for r
#     Reverse r
#       If var1type=1 AND var2type=Group A: O, MA, ML2, ML29  
#       If var1type=O, MA, ML2, or ML29 AND var2type=MB, MC, MD, ME
#       If var1type=MB, MC, MD, ME AND var2type=O, MA, ML2, or ML29
#     Same r
#       If var1type=1 AND var2type=Group B: MB, MC, MD, ME, 1 
#       If var1type=var2type 
#       If both var group A var1type=O, MA, ML2, or ML29 AND var2type=O, MA, ML2, or ML29
#       If both var group B var1type=MB, MC, MD, ME, 1 AND var2type=MB, MC, MD, ME, 1

CBTStudiesCohensDAggregated <- CBTStudiesCohensD %>%
  mutate(
    # === Reversal Cases ===
    reversalRule1 = ifelse(var1type == "1" & var2type %in% c("O", "MA", "ML2", "ML29"), TRUE, FALSE),
    reversalRule2 = ifelse(var1type %in% c("O", "MA", "ML2", "ML29") & var2type %in% c("MB", "MC", "MD", "ME"), TRUE, FALSE),
    reversalRule3 = ifelse(var1type %in% c("MB", "MC", "MD", "ME") & var2type %in% c("O", "MA", "ML2", "ML29"), TRUE, FALSE),
    
    # # === Same Direction Cases ===
    # sameDirectionRule1 = ifelse(var1type == "1" & var2type %in% c("MB", "MC", "MD", "ME", "1"), TRUE, FALSE),
    # sameDirectionRule2 = ifelse(var1type == var2type, TRUE, FALSE),
    # sameDirectionRule3 = ifelse(var1type %in% c("O", "MA", "ML2", "ML29") & var2type %in% c("O", "MA", "ML2", "ML29"), TRUE, FALSE),
    # sameDirectionRule4 = ifelse(var1type %in% c("MB", "MC", "MD", "ME", "1") & var2type %in% c("MB", "MC", "MD", "ME", "1"), TRUE, FALSE),
    
    # Combine reversal logic (TRUE if any reversal rule is triggered)
    isReverseDirection = reversalRule1 | reversalRule2 | reversalRule3,
    
    # Apply reversal
    Cohen_d = ifelse(isReverseDirection, -Cohen_d, Cohen_d),
    r = ifelse(isReverseDirection, -r, r)
  )

CBTStudiesHedgesGAggregated <- CBTStudiesHedgesG %>%
  mutate(
    # === Reversal Cases ===
    reversalRule1 = ifelse(var1type == "1" & var2type %in% c("O", "MA", "ML2", "ML29"), TRUE, FALSE),
    reversalRule2 = ifelse(var1type %in% c("O", "MA", "ML2", "ML29") & var2type %in% c("MB", "MC", "MD", "ME"), TRUE, FALSE),
    reversalRule3 = ifelse(var1type %in% c("MB", "MC", "MD", "ME") & var2type %in% c("O", "MA", "ML2", "ML29"), TRUE, FALSE),
    
    # # === Same Direction Cases ===
    # sameDirectionRule1 = ifelse(var1type == "1" & var2type %in% c("MB", "MC", "MD", "ME", "1"), TRUE, FALSE),
    # sameDirectionRule2 = ifelse(var1type == var2type, TRUE, FALSE),
    # sameDirectionRule3 = ifelse(var1type %in% c("O", "MA", "ML2", "ML29") & var2type %in% c("O", "MA", "ML2", "ML29"), TRUE, FALSE),
    # sameDirectionRule4 = ifelse(var1type %in% c("MB", "MC", "MD", "ME", "1") & var2type %in% c("MB", "MC", "MD", "ME", "1"), TRUE, FALSE),
    
    # Combine reversal logic (TRUE if any reversal rule is triggered)
    isReverseDirection = reversalRule1 | reversalRule2 | reversalRule3,
    
    # Apply reversal
    Hedges_g = ifelse(isReverseDirection, -Hedges_g, Hedges_g),
    r = ifelse(isReverseDirection, -r, r)
  )



# Transform Cohen's d & Hedges' g into point-biserial correlation via Aaron et al. (1998) formula
CBTStudiesCohensDAggregated <- CBTStudiesCohensDAggregated %>%
  mutate(
    r_pb = mapply(getPointBiserialCorrelation, Cohen_d, N_total)
  )

CBTStudiesHedgesGAggregated <- CBTStudiesHedgesGAggregated %>%
  mutate(
    r_pb = mapply(getPointBiserialCorrelation, Hedges_g, N_total)
  )



# CBTStudiesCohensDAggregated <- CBTStudiesCohensDAggregated %>%
#   group_by(studyid, var1type, var2type) %>%
#   summarise(
#     Cohen_d = if (all(is.na(Cohen_d))) NA else mean(Cohen_d, na.rm = TRUE),
#     r = if (all(is.na(r))) NA else mean(r, na.rm = TRUE),
#     r_pb = if (all(is.na(r_pb))) NA else mean(r_pb, na.rm = TRUE),
#     N_total = if (all(is.na(N_total))) NA else mean(N_total, na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# CBTStudiesHedgesGAggregated <- CBTStudiesHedgesGAggregated %>%
#   group_by(studyid, var1type, var2type) %>%
#   summarise(
#     Hedges_g = if (all(is.na(Hedges_g))) NA else mean(Hedges_g, na.rm = TRUE),
#     r = if (all(is.na(r))) NA else mean(r, na.rm = TRUE),
#     r_pb = if (all(is.na(r_pb))) NA else mean(r_pb, na.rm = TRUE),
#     N_total = if (all(is.na(N_total))) NA else mean(N_total, na.rm = TRUE)
#   ) %>%
#   ungroup()



# Preprocess the data sets to make them webMASEM-compatible

# Add Y_X, Y_M, and M_X columns. If the relationship is between continious variables, the r_pb returns null and "r" is used instead. 
CBTStudiesCohensDAggregated <- CBTStudiesCohensDAggregated %>%
  mutate(
    Y_X = case_when(var1type == "1" & var2type == "O" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_),
    Y_M = case_when(var1type == "MA" & var2type == "O" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_),
    M_X = case_when(var1type == "1" & var2type == "MA" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_)
  )

CBTStudiesHedgesGAggregated <- CBTStudiesHedgesGAggregated %>%
  mutate(
    Y_X = case_when(var1type == "1" & var2type == "O" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_),
    Y_M = case_when(var1type == "MA" & var2type == "O" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_),
    M_X = case_when(var1type == "1" & var2type == "MA" ~ ifelse(!is.na(r_pb), r_pb, r), TRUE ~ NA_real_)
  )


# Collapse data so that each studyid has one row with all available variable pairs
webMASEMCohensD <- CBTStudiesCohensDAggregated %>%
  group_by(studyid) %>%
  summarise(
    M_X = if (all(is.na(M_X))) NA else mean(M_X, na.rm = TRUE),
    Y_X = if (all(is.na(Y_X))) NA else mean(Y_X, na.rm = TRUE),  # Keep NA if all values are NA
    Y_M = if (all(is.na(Y_M))) NA else mean(Y_M, na.rm = TRUE)
  ) %>%
  ungroup()

webMASEMHedgesG <- CBTStudiesHedgesGAggregated %>%
  group_by(studyid) %>%
  summarise(
    M_X = if (all(is.na(M_X))) NA else mean(M_X, na.rm = TRUE),
    Y_X = if (all(is.na(Y_X))) NA else mean(Y_X, na.rm = TRUE),  # Keep NA if all values are NA
    Y_M = if (all(is.na(Y_M))) NA else mean(Y_M, na.rm = TRUE)
  ) %>%
  ungroup()


# Sample sizes are aggregated a second time here:
#   - This step computes one overall average sample size per study,
#     across all mediation paths, to represent the total N in the final WebMASEM input.
CBTStudiesSampleSize <- ESDataCBTStudies %>%
  filter(!is.na(n_t), !is.na(n_c)) %>%
  mutate(total_N = n_t + n_c) %>%
  group_by(studyid) %>%
  summarise(N = mean(total_N)) %>% 
  ungroup()

# Merge the sample size data with the WebMASEM-compatible data frame
webMASEMCohensD <- webMASEMCohensD %>%
  left_join(CBTStudiesSampleSize, by = "studyid")

webMASEMHedgesG2 <- webMASEMHedgesG %>%
  left_join(CBTStudiesSampleSize, by = "studyid")


# Write the data sets into Excel files
write.csv(webMASEMCohensD, "~/Desktop/repo/dataAnalysisEffectSizeMASEM/output/webMASEM_CohensD.csv")
write.csv(webMASEMHedgesG2, "~/Desktop/repo/dataAnalysisEffectSizeMASEM/output/webMASEM_HedgesG.csv")
