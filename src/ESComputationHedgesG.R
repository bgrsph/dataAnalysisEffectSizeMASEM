library(dplyr)

# 1- Import datasets
ESdata <- read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/ESdata.csv")
unimeta <- read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Univariate/unimeta.CSV")

# 2- Extract the CBT study IDs from the unimeta.csv (i.e., IPT column = 0)
cbtStudyIDs <- unimeta %>%
  filter(IPT == 0) %>%  # Select only CBT trials (Study 1)
  pull(studyid) %>%  # Extract study IDs
  unique()  # Keep unique study IDs

# 3- Extract the relevant effect size data of CBT studies from ESdata.csv based on mediation path and time of measurement
#     -Mediation paths (as shown in Table 1): 
#       - No time-point CBT (0) → Post-treatment Negative Cognition (3)
#       - Post-treatment Negative Cognition (3) → Post-treatment Depression (3)
#       - No time-point CBT (0) → Post-treatments Depression (3)
ESDataCBTStudies <- ESdata %>%
  filter(
    studyid %in% cbtStudyIDs &
      (
        (var1type == "1" & var2type == "MA" & Var1time == 0 & Var2time == 3) |    # Path a: CBT Treatment (1) → Negative Cognition (MA)
          (var1type == "MA" & var2type == "O" & Var1time == 3 & Var2time == 3) |  # Path b: Negative Cognition (MA) → Depression Severity (O)
          (var1type == "1" & var2type == "O" & Var1time == 0 & Var2time == 3)     # Path c': CBT Treatment (1) → Depression Severity
      )
  ) %>%
  select(studyid, var1type, Var1time, var2type, Var2time, mean_t, SD_t, n_t, mean_c, SD_c, n_c, r, NumOppDir,n_r) # Only include parameters required to compute Cohen's d & Hedges' g

# 4- Compute Hedges g & reverse it's sign if the relationship has 1 variable that is in the opposite direction, else same direction. 
ESDataHedgesG <- mutate (ESDataCBTStudies,
                             
                             #Compute Hedges'g
                             hedges_g = (1 - (3 / (4 * (n_t + n_c) - 9))) * ((mean_t - mean_c) / (sqrt(((SD_t^2) * (n_t - 1) + (SD_c^2) * (n_c - 1)) / (n_t + n_c - 2)))),
                          
                             #Reverse score hedges_g that are comprised of 1 variable that is in the opposite direction, otherwise same direction
                             hedges_g_opp = ifelse((NumOppDir==1), -hedges_g, hedges_g),
                         
                             #Reverse score r that is comprised of 1 variable that is in the opposite direction, otherwise same direction
                             r = ifelse((NumOppDir==1), -r, r)
)


#5- Based on the code of meta-analysis (see ESComputation.R), reverse r if one of the conditions below applies
#   -If var1type=1 AND var2type=Group A: O, MA, ML2, ML29
#   -If var1type=O, MA, ML2, or ML29 AND var2type=MB, MC, MD, ME
#   -If var1type=MB, MC, MD, ME AND var2type=O, MA, ML2, or ML29
ESDataHedgesGSignsReversed <- mutate (ESDataHedgesG,
                               var1grp = ifelse((var1type=="O"|var1type=="MA"|var1type=="ML2"|var1type=="ML29"), "A", "B"),
                               var2grp = ifelse((var2type=="O"|var2type=="MA"|var2type=="ML2"|var2type=="ML29"), "A", "B"),
                               hedges_g = ifelse((var1grp==var2grp), hedges_g_opp, -hedges_g_opp), 
                               r = ifelse((var1grp==var2grp), r, -r),
)


# 6- Aggregate sample size and Hedges'g in case same relationship within a study exists with different effect sizes and/or sample sizes
ESDataAggregated <- ESDataHedgesGSignsReversed %>%
  group_by(studyid, var1type, var2type) %>%
  summarise(
    hedges_g_avg = mean(hedges_g, na.rm = TRUE), # Average Hedges' g (for relationships that may appear more than once per study)
    r_avg = mean(r, na.rm = TRUE),               # Average raw r (fallback when Hedges' g is missing)
    n_total_avg = mean(n_t + n_c, na.rm = TRUE)  # Average total sample size (used to fill in missing sample sizes later)
  ) %>%
  ungroup() %>%
  # Step: Derive rpb from hedges_g_avg when available; otherwise fallback to r_avg
  mutate(
    # Compute rpb from Hedges' g when it's available (Aaron et al.)
    rpb = ifelse(!is.na(hedges_g_avg), (hedges_g_avg / sqrt(hedges_g_avg^2 + 4 - (8 / n_total_avg))), r_avg) # if Hedges'g not reported, fallback to reported r (which is the case for continious relationship of M <--> Y)
)  

# 6- Compute webMASEM-compatible dataset with one row per study
webMASEMHedgesG <- ESDataAggregated %>%
  mutate(
    # Recode effect paths
    M_X = ifelse(var1type == "1"  & var2type == "MA", rpb, NA),  # CBT → Neg. Cognition
    Y_X = ifelse(var1type == "1"  & var2type == "O",  rpb, NA),   # CBT → Depression
    Y_M = ifelse(var1type == "MA" & var2type == "O",  rpb, NA)  # Neg. Cognition → Depression
  ) %>%
  # Now group by studyid and compute average per path. Here, taking the mean does not takes the mean, it just collapses all the aggregated data into one row. 
  group_by(studyid) %>%
  summarise(
    M_X = if (all(is.na(M_X))) NA else mean(M_X, na.rm = TRUE),
    Y_X = if (all(is.na(Y_X))) NA else mean(Y_X, na.rm = TRUE),
    Y_M = if (all(is.na(Y_M))) NA else mean(Y_M, na.rm = TRUE),
    N   = mean(n_total_avg, na.rm = TRUE)
  ) %>%
  ungroup()

write.csv(webMASEMHedgesG, "~/Desktop/repo/dataAnalysisEffectSizeMASEM/output/webMASEMHedgesG.csv")


