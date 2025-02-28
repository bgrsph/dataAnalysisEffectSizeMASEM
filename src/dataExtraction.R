#Import (custom) libraries
library(dplyr)
library(readxl)
library(writexl)
library(tidyr)
library(gt)
source("~/Desktop/repo/dataAnalysisEffectSizeMASEM/src/EffectSizeCalculator.R")


#Import datasets
ESdata<-read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/ESdata.csv")
unimeta <- read.csv("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Univariate/unimeta.CSV")
measureCharacteristicsData<-read_xlsx("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/MeasureCharacteristics.xlsx")
studyCharacteristicsData <- read_xlsx("~/Desktop/repo/dataAnalysisEffectSizeMASEM/data/Files/Study Measure ES Data/StudyCharacteristics.xlsx")

# Identify Study 1 (CBT-only) study IDs
cbt_study_ids <- unimeta %>%
  filter(IPT == 0) %>%  # Select only CBT trials (Study 1)
  pull(studyid) %>%  # Extract study IDs
  unique()  # Keep unique study IDs

# Filter ESdata for CBT trials from Study 1
cbt_mediation_es <- ESdata %>%
  filter(
    studyid %in% cbt_study_ids &  # Use automatically identified CBT study IDs
      (
        (var1type == "1" & var2type == "MA") |  # Path a: CBT → Negative Cognition
          (var1type == "MA" & var2type == "O") |  # Path b: Negative Cognition → Depression
          (var1type == "1" & var2type == "O")     # Path c: CBT → Depression
      )
  ) %>%
  select(studyid, var1type, Var1time, var2type, Var2time, 
         mean_t, SD_t, n_t, mean_c, SD_c, n_c, r, n_r)


# Filter the studies where mediation path is chronologically preserved. 
#     var1Time: 0 (pre-treatment), 1 (immediate post-treatment), 2(short-term follow-up), 3(long-term follow up)

#     pre-treatment X (CBT) --> post-treatment M (negative cognition)
#     post-treatment M (negative cognition) --> post-treatment/follow-up Y (depression)
#     pre-treatment X (CBT) --> post-treatment/follow-up Y (depression)


cbt_mediation_es_filtered <- cbt_mediation_es %>%
  filter(
    (var1type == "1" & var2type == "MA" & Var1time == 0 & Var2time == 3) |  # Path A: CBT → Negative Cognition
      (var1type == "MA" & var2type == "O" & Var1time == 3 & Var2time %in% c(3,4)) |  # Path B: Negative Cognition → Depression
      (var1type == "1" & var2type == "O" & Var1time == 0 & Var2time %in% c(3,4))    # Path C: CBT → Depression
  )


# Count all the possible timings with all possible mediation path to choose the ones with the most data
mediation_paths <- expand.grid(
  var1type = c("1", "MA"), 
  var2type = c("MA", "O")
) %>%
  mutate(
    Path = case_when(
      var1type == "1"  & var2type == "MA" ~ "X → M",
      var1type == "MA" & var2type == "O"  ~ "M → Y",
      var1type == "1"  & var2type == "O"  ~ "X → Y",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Path))

# Define all possible Var1time → Var2time combinations
time_combinations <- expand.grid(
  Var1time = unique(cbt_mediation_es$Var1time),
  Var2time = unique(cbt_mediation_es$Var2time)
) %>%
  mutate(Timing = paste(Var1time, "→", Var2time))

# Generate full table with zero-filled missing combinations
mediation_summary <- cbt_mediation_es %>%
  mutate(
    Path = case_when(
      var1type == "1"  & var2type == "MA" ~ "X → M",
      var1type == "MA" & var2type == "O"  ~ "M → Y",
      var1type == "1"  & var2type == "O"  ~ "X → Y",
      TRUE ~ NA_character_
    ),
    Timing = paste(Var1time, "→", Var2time)
  ) %>%
  filter(!is.na(Path)) %>%
  count(Path, Timing) %>%
  complete(Path = mediation_paths$Path, Timing = time_combinations$Timing, fill = list(n = 0)) %>%
  pivot_wider(names_from = Timing, values_from = n, values_fill = 0)

# Create a nice table plot
mediation_summary %>%
  gt() %>%
  tab_header(
    title = "Mediation Paths by Timing Combinations",
    subtitle = "Counts of Observations for Each Mediation Path and Measurement Timing"
  ) %>%
  cols_label(Path = "Mediation Path") %>%
  data_color(
    columns = -Path,
    colors = scales::col_numeric(
      palette = "Greens",
      domain = NULL
    )
  ) %>%
  fmt_number(
    columns = -Path,
    decimals = 0
  ) %>%
  tab_options(
    table.font.size = px(14),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(14)
  )


# Calculate average sample size
cbt_mediation_es_filtered %>%
  mutate(avg_sample_size = (n_t + n_c)) %>%
  summarise(mean_sample_size = mean(avg_sample_size, na.rm = TRUE))



# in the measure characteristics: varID: (studyid).(mediator)(category letter “a-e” or “L” for low priority).(2-digit number, 1st one 01,
#   2nd one 02, 3rd one 03, no particular order required)

# here you can extract either only parent, child, or aggregate the multiple effect sizes in ONE study. 
