
#Raw and computed ES 

library(dplyr)
source("~/Desktop/repo/dataAnalysisEffectSizeMASEM/src/effectSizeFunctions.R")

#Load dataset--select ESdata.cvs when window opens
#ESdata<-read.csv(file.choose())
ESdata <- read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/ESdata.csv")

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
  select(studyid, var1type, Var1time, var2type, Var2time, mean_t, SD_t, n_t, mean_c, SD_c, n_c, r, NumOppDir,n_r, TCid, var1mul, var2mul) # Only include parameters required to compute Cohen's d & Hedges' g


ESdata1 <- mutate (ESDataCBTStudies,
                   
                   #1. Compute pooled SD and then Hedges' g from ns, means, and SDs
                   Spooled1 = sqrt((SD_t*SD_t*(n_t-1) + SD_c*SD_c*(n_c-1))/(n_t+n_c-2)), 
                   g1 = (mean_t - mean_c)/Spooled1, 
                   
                   #2. Apply small sample correction
                   ga = g1 - (3*g1 / (4*(n_t + n_c) - 9) ),
                   
                   #3. Compute ns for rs that are converted from means and SDs and ns
                   n_r1 = ifelse((is.na(Spooled1)), n_r, (n_t+n_c)),
                   
                   #4. Convert corrected Hedges' g into rs, if ga is missing, then input entered r, otherwise convert ga to r                                       
                   rga = ifelse((is.na(ga)), r, ga/(sqrt((ga*ga) + ((n_r*n_r)/(n_t*n_c)) ) )   ), 
                   
                   #5. Reverse score rs that are comprised of 1 variable that is in the opposite direction, otherwise same direction
                   rgaopp = ifelse((NumOppDir==1), rga*-1, rga),
                   
                   #6. Reverse score gs that are comprised of 1 variable that is in the opposite direction, otherwise same direction
                   gaopp = ifelse((NumOppDir==1), ga*-1, ga),
                   
                   #7. Compute variance of corrected Hedges' g
                   Vga  =       ((n_t + n_c) / (n_t*n_c) ) + ((ga*ga) / (2*(n_t + n_c)) ),
                   
                   #8. Compute SE of un/corrected Hedges' g
                   SEga = sqrt( ((n_t + n_c) / (n_t*n_c) ) + ((ga*ga) / (2*(n_t + n_c)) ) ),
                   
                   #9. BUGRA
                   hedges_g = (1 - (3 / (4 * (n_t + n_c) - 9))) * ((mean_t - mean_c) / (sqrt(((SD_t^2) * (n_t - 1) + (SD_c^2) * (n_c - 1)) / (n_t + n_c - 2)))),
                   
                   #10.BUGRA: convert Hedges' g into r using Aaron et al. formula. If needed values are NA, use the reported r (which is the case for continous bivariate relationship M-Y)
                   rpb_converted = ifelse((is.na(hedges_g)), r, hedges_g / (sqrt(hedges_g^2 + 4 - (8/(n_t + n_c))))),
              
                   #11. BUGRA
                   rpb_converted_opp = ifelse((NumOppDir==1), -rpb_converted, rpb_converted),
                   
                   #12. BUGRA
                   hedges_g_opp = ifelse((NumOppDir==1), -hedges_g, hedges_g)
                   
)


# Collapse multiple measures of the same category at the same timepoint within study/tx-control comparison

# ESdata2 <- ESdata1 %>% 
#   group_by(studyid, TCid, Var1time, var1type, var1mul, Var2time, var2type, var2mul) %>% 
#   summarise(rpb_converted_opp=mean(rpb_converted_opp), hedges_g_opp = mean(hedges_g_opp), rgaopp=mean(rgaopp), gaopp=mean(gaopp), SEga=mean(SEga), Vga=mean(Vga), n_r1=mean(n_r1), n_t=mean(n_t), n_c=mean(n_c)) 
# 
# ESdata3 <- ESdata2 %>% 
#   group_by(studyid, TCid, Var1time, var1type, var1mul, Var2time, var2type, var2mul) %>% 
#   summarise(rpb_converted_opp=mean(rpb_converted_opp), hedges_g_opp = mean(hedges_g_opp), rgaopp=mean(rgaopp), gaopp=mean(gaopp), SEga=mean(SEga), Vga=mean(Vga), n_r1=mean(n_r1), n_t=mean(n_t), n_c=mean(n_c)) 


#Apply reverse rules for r
#Reverse r
#If var1type=1 AND var2type=Group A: O, MA, ML2, ML29  
#If var1type=O, MA, ML2, or ML29 AND var2type=MB, MC, MD, ME
#If var1type=MB, MC, MD, ME AND var2type=O, MA, ML2, or ML29

#Same r
#If var1type=1 AND var2type=Group B: MB, MC, MD, ME, 1 
#If var1type=var2type 
#If both var group A var1type=O, MA, ML2, or ML29 AND var2type=O, MA, ML2, or ML29
#If both var group B var1type=MB, MC, MD, ME, 1 AND var2type=MB, MC, MD, ME, 1

ESdata4 <- mutate (ESdata1,
                   var1grp = ifelse((var1type=="O"|var1type=="MA"|var1type=="ML2"|var1type=="ML29"), "A", "B"),
                   var2grp = ifelse((var2type=="O"|var2type=="MA"|var2type=="ML2"|var2type=="ML29"), "A", "B"),
                   rga = ifelse((var1grp==var2grp), rgaopp, rgaopp*-1),
                   ga = ifelse((var1grp==var2grp), gaopp, gaopp*-1),
                   rpb_converted = ifelse((var1grp==var2grp), rpb_converted_opp, rpb_converted_opp*-1),
                   hedges_g = ifelse((var1grp==var2grp), hedges_g_opp, hedges_g_opp*-1)
) 




#Use to create data files for univariate meta-analysis and MASEM.
webMASEM <- ESdata4 %>%
  mutate(
    M_X = ifelse(var1type == "1" & var2type == "MA", rpb_converted, NA),
    Y_M = ifelse(var1type == "MA" & var2type == "O", rpb_converted, NA),
    Y_X = ifelse(var1type == "1" & var2type == "O", rpb_converted, NA)
  ) %>%
  group_by(studyid) %>%
  summarise(
    M_X = if (all(is.na(M_X))) NA else mean(M_X, na.rm = TRUE),
    Y_X = if (all(is.na(Y_X))) NA else mean(Y_X, na.rm = TRUE),
    Y_M = if (all(is.na(Y_M))) NA else mean(Y_M, na.rm = TRUE),
    N = mean(n_t + n_c, na.rm = TRUE)
  ) %>%
  ungroup()


write.csv(webMASEM, "~/Desktop/repo/dataAnalysisEffectSizeMASEM/output/webMASEM.csv")

cat("Final row counts — ESdata1:", nrow(ESdata1), 
    "| ESdata4 (after reversal):", nrow(ESdata4), 
    "| webMASEM:", nrow(webMASEM), "\n")