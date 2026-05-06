rm(list=ls())
library(dplyr)

path <- getwd()
setwd(path) 
list_analyses <- c(
  "00_2_food_data_preprocess.R",
  "01_descriptive_food_source_revised.R",
  "02_1_diversity_food_frequency.R",
  "02_2_meal_composition.R",
  "03_1_bmi_and_adiposity.R",
  "03_2_bootrapped_age_treshold_adiposity.R",
  "03_3_GDD_data.R",
  "04_1_sugar.R",
  "04_2_smoking.R",
  "04_4_medicine.R",
  "05_1_build_fig3_and_fig4.R",
  "05_2_non_food_compile.R",
  "06_reported_stats.R"
)

for(file in list_analyses){
  print(file)
  source(file)
  setwd(path) 
}
