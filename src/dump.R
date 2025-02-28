#Import libraries
library(dplyr)
library(readxl)
library(writexl)


#Import datasheets
ESdata<-read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/ESdata.csv")
unimeta <- read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Univariate/unimeta.CSV")
measureCharacteristicsData<-read_xlsx("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/MeasureCharacteristics.xlsx")

#Extract the data for the selected mediation model where X (CBT or not), M (negative cognition), Y(depression outcome)

# Extract Path a (CBT → Negative Cognition)
pathA <- unimeta %>%
  filter(IPT == 0, VarPair == "0_1_3_MA") %>%
  mutate(Path = "a") 

# Extract Path b (Negative Cognition → Depression Severity)
pathB <- unimeta %>%
  filter(IPT == 0, VarPair == "3_MA_3_O") %>%
  mutate(Path = "b") 

# Extract Path c (CBT → Depression Severity) - Handling Multiple Possible VarPairs
pathC <- unimeta %>%
  filter(IPT == 0, VarPair == "0_1_3_O") %>%
  mutate(Path = "c") 

extractedDataset <-  bind_rows(pathA, pathB, pathC)
write.csv(extractedDataset, "ExtractedData.csv", row.names = FALSE)












# Load extracted dataset
data <- read.csv("ExtractedData.csv")

Path_a <- data %>%
  filter(VarPair == "0_1_3_MA") %>%
  select(studyid, X_M = rga, n_r1)

# Extract Path b (Negative Cognition → Outcome)
Path_b <- data %>%
  filter(VarPair == "3_MA_3_O") %>%
  select(studyid, M_Y = rga)

# Extract Path c (CBT → Outcome)
Path_c <- data %>%
  filter(VarPair == "0_1_3_O") %>%
  select(studyid, X_Y = rga)

# Merge datasets (keep only studies that appear in all paths)
MASEM_ready <- Path_a %>%
  inner_join(Path_b, by = "studyid") %>%
  inner_join(Path_c, by = "studyid")

write.csv(MASEM_ready, "MASEMData.csv", row.names = FALSE)




cbt_mediation_es <- cbt_mediation_es %>%
  mutate(cohens_d = mapply(getCohensD, mean_t, SD_t, n_t, mean_c, SD_c, n_c))

cbt_mediation_es <- cbt_mediation_es %>%
  mutate(hedges_g = mapply(getHedgesG, mean_t, SD_t, n_t, mean_c, SD_c, n_c))

#Add an dummy column for transformedr for now
# Append an empty column to cbt_mediation_es
cbt_mediation_es <- cbt_mediation_es %>%
  mutate(r_transformed = NA)

# Calculating sample sizes (control and treatment) for information

# Check summary statistics for sample sizes
sample_size_summary <- cbt_mediation_es %>%
  summarise(
    Mean_n_t = mean(n_t, na.rm = TRUE),
    Median_n_t = median(n_t, na.rm = TRUE),
    Min_n_t = min(n_t, na.rm = TRUE),
    Max_n_t = max(n_t, na.rm = TRUE),
    Q1_n_t = quantile(n_t, 0.25, na.rm = TRUE),
    Q3_n_t = quantile(n_t, 0.75, na.rm = TRUE),
    
    Mean_n_c = mean(n_c, na.rm = TRUE),
    Median_n_c = median(n_c, na.rm = TRUE),
    Min_n_c = min(n_c, na.rm = TRUE),
    Max_n_c = max(n_c, na.rm = TRUE),
    Q1_n_c = quantile(n_c, 0.25, na.rm = TRUE),
    Q3_n_c = quantile(n_c, 0.75, na.rm = TRUE)
  )

# Print the summary statistics
print(sample_size_summary)



# Compute imbalance metrics
imbalance_metrics <- cbt_mediation_es %>%
  mutate(
    R = n_t / n_c,  # Ratio of group sizes
    VIF = (n_t + n_c)^2 / (4 * n_t * n_c),  # Variance Inflation Factor
    I = abs(n_t - n_c) / (n_t + n_c)  # Relative Imbalance Index
  ) %>%
  summarise(
    Mean_R = mean(R, na.rm = TRUE),
    Median_R = median(R, na.rm = TRUE),
    Min_R = min(R, na.rm = TRUE),
    Max_R = max(R, na.rm = TRUE),
    
    Mean_VIF = mean(VIF, na.rm = TRUE),
    Median_VIF = median(VIF, na.rm = TRUE),
    Min_VIF = min(VIF, na.rm = TRUE),
    Max_VIF = max(VIF, na.rm = TRUE),
    
    Mean_I = mean(I, na.rm = TRUE),
    Median_I = median(I, na.rm = TRUE),
    Min_I = min(I, na.rm = TRUE),
    Max_I = max(I, na.rm = TRUE)
  )

# Print results
print(imbalance_metrics)

