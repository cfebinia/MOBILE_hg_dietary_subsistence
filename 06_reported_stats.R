rm(list=ls()) 

library(tidyverse)

# helper functions
####
calculate_missingness <- function(data, vars) {
  subset_data <- data[vars]
  
  missing_counts <- colSums(is.na(subset_data))
  total_missing <- sum(missing_counts)
  
  num_observations <- nrow(data)
  total_possible <- num_observations * length(vars)
  missing_proportion <- total_missing / total_possible
  
  results <- list(
    n_obs = num_observations,
    total_na = total_missing,
    proportion_na = missing_proportion
  )
  
  return(results)
}

### Sample Match

full_list <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>%
  select(sampleid, age, sex) %>%
  mutate(population=substr(sampleid,1,3)) %>%
  mutate(group=ifelse(population %in% c("BTU","ORT","ORS"),"earlyTransition",
                      ifelse(population=="LDY","Agriculture",
                      "lateTransition"))) %>%
  mutate(population = factor(population, c("BTU","ORT","ORS","APT","TBU","BSP","LDY")),
         group = factor(group, c("earlyTransition","lateTransition","Agriculture")),
         sex = factor(sex, c("F","M")))

file_list = list.files("input_files/", 
                       pattern = ".*sample_list\\.rds$", 
                       full.names = TRUE)
dataset_names <- gsub("^input_files/|\\_sample_list.rds$", "", file_list)

subset_samples <- lapply(file_list, readRDS)
names(subset_samples) <- dataset_names

for(dat in dataset_names){
  full_list[full_list$sampleid %in% subset_samples[[dat]],dat] <- subset_samples[[dat]]
}

# count samples
countA <- full_list %>%
  select(total = sampleid, sex, all_of(dataset_names)) %>%
  summarise(
    across(-sex, list(
      val1 = ~ as.character(n_distinct(.x, na.rm = TRUE)),
      val2 = ~ paste0("(", n_distinct(.x[sex == "M"], na.rm = TRUE), 
                      "|", n_distinct(.x[sex == "F"], na.rm = TRUE), ")")
    ))
  ) %>%
  mutate(group="all") %>%
  pivot_longer(
    cols = -group, 
    names_to = c("dataset", ".value"), 
    names_pattern = "(.*)_val(.+)"
  ) %>%
  mutate(n = paste(`1`, `2`)) %>%
  mutate(n = gsub("^0 \\(0\\|0\\)$", "—", n)) %>%
  select(group, dataset, n) %>%
  pivot_wider(
    names_from = dataset, 
    values_from = n
  ) %>%
  rename(cohort=group)

countB <- full_list %>%
  select(total = sampleid, group, sex, all_of(dataset_names)) %>%
  group_by(group) %>%
  summarise(
    across(-sex, list(
      val1 = ~ as.character(n_distinct(.x, na.rm = TRUE)),
      val2 = ~ paste0("(", n_distinct(.x[sex == "M"], na.rm = TRUE), 
                      "|", n_distinct(.x[sex == "F"], na.rm = TRUE), ")")
    ))) %>%
  pivot_longer(
    cols = -group, 
    names_to = c("dataset", ".value"), 
    names_pattern = "(.*)_val(.+)"
  ) %>%
  mutate(n = paste(`1`, `2`)) %>%
  mutate(n = gsub("^0 \\(0\\|0\\)$", "—", n)) %>%
  select(group, dataset, n) %>%
  pivot_wider(
    names_from = dataset, 
    values_from = n
  ) %>%
  rename(cohort=group)
  

countC <- full_list %>%
  select(total = sampleid, population, sex, all_of(dataset_names)) %>%
  group_by(population) %>%
  summarise(
    across(-sex, list(
      val1 = ~ as.character(n_distinct(.x, na.rm = TRUE)),
      val2 = ~ paste0("(", n_distinct(.x[sex == "M"], na.rm = TRUE), 
                      "|", n_distinct(.x[sex == "F"], na.rm = TRUE), ")")
    ))) %>%
  pivot_longer(
    cols = -population, 
    names_to = c("dataset", ".value"), 
    names_pattern = "(.*)_val(.+)"
  ) %>%
  mutate(n = paste(`1`, `2`)) %>%
  mutate(n = gsub("^0 \\(0\\|0\\)$", "—", n)) %>%
  select(population, dataset, n) %>%
  pivot_wider(
    names_from = dataset, 
    values_from = n
  ) %>% 
  rename(cohort=population)


consolidated_sample_count <- rbind(countA, countB, countC) %>%
  select(cohort, source, ffq, meal_compos, obesity, medicine, smoking, total) %>%
  rename(food_source=source)

write.table(consolidated_sample_count, 
            "out/consolidated_sample_count.tsv", quote = F, sep = "\t", row.names = F)

##################
# Reported Stats
#################
# load datasets
df.src <- read.delim("input_files/Food_Source_Compiled_final - IDPROP.tsv") %>% select(-population)
df.src.imp <- read.delim("out/food_source_data_processed.tsv")  %>% select(-population,-group)
df.idds <- read.delim("out/idds_data_processed.tsv")
df.ffq <- read.delim("input_files/ffq_corrected_imputed_median.tsv")
df.mealcompos <- read.delim("out/meal_compos_data_processed.tsv")
df.bio <- read.delim("out/anthrop_data_processed.tsv") %>% rename(sampleid=sid)
df.bio.imp <- read.delim("out/anthrop_data_median_imputed_processed.tsv") %>% rename(sampleid=sid) # if using imputed
df.gddbio <- read.delim("out/gdd_biomarkers_processed.tsv") %>% rename(sampleid=sid)
df.sugar <- read.delim("out/sugar_data_processed.tsv")
df.smoking <- read.delim("out/smoking_data_processed.tsv") %>% rename(sampleid=sid)
df.medicine <- read.delim("out/medicine_data_processed.tsv") %>% rename(sampleid=sid)


pop.ord <- c("BTU","ORT","ORS","APT","TBU","BSP","LDY")

###########
# Food Source
###########
df.src <- df.bio %>%
  select(sampleid, sex, population, group) %>%
  full_join(df.src.imp, by="sampleid", keep = F) %>%
  filter(sampleid %in% full_list$source)

length(unique(df.src$sampleid))

df.src$population <- factor(df.src$population, pop.ord)

df.src %>%
  group_by(population) %>%
  filter(source=="wild") %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(proportion), 
            sd = sd(proportion), 
            med = median(proportion),
            q25 = quantile(proportion, 0.25),
            q75 = quantile(proportion, 0.75),
            q95 = quantile(proportion, 0.95),
            max = max(proportion, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

df.src %>%
  filter(!is.na(proportion)) %>%
  group_by(population) %>%
  filter(source=="market") %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(proportion), 
            sd = sd(proportion), 
            med = median(proportion),
            q25 = quantile(proportion, 0.25),
            q75 = quantile(proportion, 0.75),
            q95 = quantile(proportion, 0.95),
            max = max(proportion, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

df.src %>%
  filter(!is.na(proportion)) %>%
  group_by(population) %>%
  filter(source=="garden") %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(proportion), 
            sd = sd(proportion), 
            med = median(proportion),
            q25 = quantile(proportion, 0.25),
            q75 = quantile(proportion, 0.75),
            q95 = quantile(proportion, 0.95),
            max = max(proportion, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

df.src %>%
  filter(!is.na(proportion)) %>%
  group_by(group) %>%
  filter(source=="garden") %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(proportion), 
            sd = sd(proportion), 
            med = median(proportion),
            q25 = quantile(proportion, 0.25),
            q75 = quantile(proportion, 0.75),
            q95 = quantile(proportion, 0.95),
            max = max(proportion, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

###########
# Food Diversity
###########
length(unique(df.idds$sampleid))

df.idds %>%
filter(!is.na(idds)) %>%
  group_by(group) %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(idds), 
            sd = sd(idds), 
            med = median(idds),
            q25 = quantile(idds, 0.25),
            q75 = quantile(idds, 0.75),
            q95 = quantile(idds, 0.95),
            max = max(idds, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

df.idds %>%
  filter(population != "LDY") %>%
  summarise(
    p_val = wilcox.test(idds ~ group)$p.value
  )

###########
# Food Frequency
###########
head(df.ffq)
length(unique(df.ffq$sampleid))

###########
# Anthropometry
###########
df.bio %>%
  filter(!is.na(bmi)) %>%
  group_by(population) %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(bmi), 
            sd = sd(bmi), 
            med = median(bmi),
            q25 = quantile(bmi, 0.25),
            q75 = quantile(bmi, 0.75),
            q95 = quantile(bmi, 0.95),
            max = max(bmi, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 1)))

select_sample <- df.bio.imp$sampleid

df.bio <- df.bio %>%
  filter(sampleid %in% select_sample)

length(unique(df.bio$sampleid))
apply(df.bio, 2, function(x) sum(is.na(x))) 
select_var <- c("bmi", "WHtR")
calculate_missingness(df.bio, select_var)

no_missing <- unique(c(df.bio$sampleid[!is.na(df.bio$bmi)], df.bio$sampleid[!is.na(df.bio$WHtR)]))
select_sample[!select_sample %in% no_missing]

table(df.bio$group)

df.bio%>%
  filter(!is.na(bmi)) %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(bmi), 
            sd = sd(bmi), 
            med = median(bmi),
            q25 = quantile(bmi, 0.25),
            q75 = quantile(bmi, 0.75),
            q95 = quantile(bmi, 0.95),
            max = max(bmi, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 1)))

df.bio %>%
  filter(!is.na(bmi)) %>%
  group_by(group) %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(bmi), 
            sd = sd(bmi), 
            med = median(bmi),
            q25 = quantile(bmi, 0.25),
            q75 = quantile(bmi, 0.75),
            q95 = quantile(bmi, 0.95),
            max = max(bmi, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 1)))

df.bio %>%
  filter(!is.na(bmi)) %>%
  group_by(group, sex) %>%
  summarise(n=n_distinct(sampleid),
            mean = mean(bmi), 
            sd = sd(bmi), 
            med = median(bmi),
            q25 = quantile(bmi, 0.25),
            q75 = quantile(bmi, 0.75),
            q95 = quantile(bmi, 0.95),
            max = max(bmi, 0.25))%>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 1)))

df.bio %>%
  filter(!is.na(is_overweight)) %>%
  group_by(group, sex, is_overweight) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  filter(is_overweight==1)

df.bio %>%
  filter(!is.na(is_overweight)) %>%
  filter(!is.na(sex)) %>%
  group_by(sex) %>%
  summarise(
    p_val = fisher.test(table(group, is_overweight))$p.value,
    .groups = "drop"
  ) %>%
  mutate(fdr_adj_p = p.adjust(p_val, method = "fdr"))

df.bio %>%
  filter(!is.na(is_underweight)) %>%
  group_by(group, sex, is_underweight) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  filter(is_underweight==1)

df.bio %>%
  filter(!is.na(is_underweight)) %>%
  filter(!is.na(sex)) %>%
  group_by(sex) %>%
  summarise(
    p_val = fisher.test(table(group, is_underweight))$p.value,
    .groups = "drop"
  ) %>%
  mutate(fdr_adj_p = p.adjust(p_val, method = "fdr"))

lmm_list <- intersect(df.bio$sampleid, df.ffq$sampleid)

df.bio %>% 
  filter(sampleid %in% lmm_list) %>%
  group_by(group) %>%
  summarise(n=n_distinct(sampleid))

###########
# GDD Biomarkers
###########
length(unique(df.gddbio$sampleid))
table(df.gddbio$sex)
table(df.gddbio$group)
table(df.gddbio$population)

table(df.gddbio$group, df.gddbio$sex)
table(df.gddbio$population, df.gddbio$sex)

lipids <- c("cholesterol","ldl","hdl","triglyceride","blood_sugar")

df.gddbio %>%
  group_by(sex) %>%
  summarise(
    across(all_of(lipids), 
           ~ list(cor.test(bmi, .x, method = "spearman", exact = FALSE))),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(lipids),
    names_to = "lipid",
    values_to = "test_obj"
  ) %>%
  mutate(
    rho = sapply(test_obj, function(x) round(x$estimate,2)),
    p_val = sapply(test_obj, function(x) x$p.value),
    fdr_p = p.adjust(p_val, method = "fdr")
  ) %>%
  select(-test_obj) %>%
  mutate(across(c(rho, p_val, fdr_p), ~ round(.x, 4)))

df.gddbio %>%
  filter(!is.na(is_hypertension)) %>%
  group_by(sex,group, is_hypertension) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  filter(is_hypertension==1)

contingency_table <- df.gddbio %>%
  filter(!is.na(is_hypertension)) %>%
  filter(sex == "M") %>%
  select(group, is_hypertension) %>%
  table()
contingency_table
fisher.test(contingency_table)

longitudinal <- read.delim("input_files/Longitudinal Samples - Matching.tsv")
longitudinal$sampleid <- gsub("-","", longitudinal$X2022_SampleID)
longitudinal <- longitudinal %>% filter(sampleid %in% full_list$sampleid)
length(unique(longitudinal$sampleid))

###########
# sugar
###########
length(unique(df.sugar$sampleid)) / n_distinct(full_list$ffq)
table(df.sugar$group)
df.sugar %>%
  filter(population != "LDY") %>%
  summarise(n=n_distinct(sampleid), mean=mean(daily_g), sd=sd(daily_g))

df.sugar %>%
  filter(population != "LDY") %>%
  mutate(is_high = ifelse(daily_g > 50,1,0)) %>%
  group_by(is_high) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  filter(is_high==1)

df.sugar %>%
  mutate(is_high = ifelse(daily_g > 50,1,0)) %>%
  group_by(group, is_high) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  filter(is_high==1)

df.sugar %>%
  group_by(group) %>%
  summarise(mean = mean(daily_g), 
            sd = sd(daily_g), 
            q95 = quantile(daily_g, 0.95), 
            max = max(daily_g))

###########
# smoking 
###########
length(unique(df.smoking$sampleid))
table(df.smoking$group)

table(df.smoking$Smoker)
table(df.smoking$Smoker)/ n_distinct(full_list$sampleid)

df.smoking %>%
  group_by(sex, Smoker) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  ungroup()

df.smoking %>%
  group_by(group, Smoker) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  ungroup()


df.smoking %>%
  group_by(sex, group, Smoker) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  ungroup()


df.smoking %>%
  filter(!is.na(intensity)) %>%
  group_by(Smoker, sex) %>%
  summarise(n = n(), 
            mean = mean(intensity),
            SD = sd(intensity),
            med = median(intensity),
            q25 = quantile(intensity, 0.25),
            q75 = quantile(intensity, 0.75),
            q95 = quantile(intensity, 0.95),
            max = max(intensity, 0.25))

df.smoking %>%
  filter(!is.na(intensity)) %>%
  filter(Smoker=="Smoker") %>%
  group_by(Smoker, sex, group) %>%
  summarise(n = n(), 
            mean = mean(intensity),
            SD = sd(intensity),
            med = median(intensity),
            q25 = quantile(intensity, 0.25),
            q75 = quantile(intensity, 0.75),
            q95 = quantile(intensity, 0.95),
            max = max(intensity, 0.25))

###########
# medicine 
###########
length(unique(df.medicine$sampleid))

table(df.medicine$all_medicine_list)

med_cat <- c("antibiotic", "analgesic", 
"respiratory", "antihyperlipidemics", 
"antigout", "antihistamine", 
"gastrointestinal", "herbal", 
"antihypertension", "other")

df.medicine %>%
  rowwise() %>%
  mutate(usage = ifelse(sum(c_across(all_of(med_cat)), na.rm = TRUE) > 0, 1, 0)) %>%
  ungroup() -> df.medicine

table(df.medicine$usage)
table(df.medicine$usage) / n_distinct(full_list$sampleid)

df.medicine %>%
  filter(usage == 1) %>%
  group_by(group, population, usage) %>%
  summarise(n_distinct(sampleid))
  
df.medicine %>%
  group_by(group, usage) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  ungroup()

df.medicine %>%
  group_by(group, population, analgesic) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  filter(analgesic == 1) %>%
  arrange(group, proportion) 

df.medicine %>%
  group_by(group, population, respiratory) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  filter(respiratory == 1) %>%
  arrange(group, proportion) 

df.medicine %>%
  group_by(group, sex, population, antibiotic) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group,2)) %>%
  filter(antibiotic == 1) %>%
  arrange(sex, group, desc(proportion)) 

unique(df.medicine$all_medicine_list[df.medicine$antibiotic==1])
       
df.medicine %>% 
  filter(antihypertension==1 | antihyperlipidemics == 1) %>%
  group_by(group, population) %>%
  summarise(n=n_distinct(sampleid),
            F=sum(sex=="F"),
            M=sum(sex=="M")) %>%
  filter(n>0)
    
       