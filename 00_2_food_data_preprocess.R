rm(list=ls())
library(dplyr)

# utils
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

### RUN
df.src <- read.delim("input_files/Food_Source_Compiled_final - IDPROP.tsv")
dff <- read.delim("input_files/ffq_corrected.tsv") # using ffq data
bio <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>% 
  select(sampleid, population, sex) %>%
  mutate(population=substr(sampleid,1,3))

# composite Sugar
dff$sugar_g <- round((dff$coffee_tea * dff$added_sugar * 12.5) / 30,1)

###
# Food Source
###
select_var <- c("garden","market","wild")

colSums(is.na(df.src[select_var]))

apply(df.src,2, function(x) sum(is.na(x)))

nNA <- calculate_missingness(df.src, select_var) # %NA
saveRDS(nNA, "out/propNA_df.src.rds")
print(nNA)

df.src <- bio %>%
  select(-population) %>%
  left_join(df.src, by="sampleid",unmatched = "drop") %>%
  filter(if_any(all_of(select_var), ~ .x != 0))# remove samples if diet is 0 across

df.src <- df.src %>%
  select(sampleid, population, select_var) %>%
  group_by(population) %>%
  # Apply population-specific median imputation across all food variables
  mutate(across(all_of(select_var), ~ replace_na(.x, round(median(.x, na.rm = TRUE), 0)))) %>%
  ungroup()

write.table(df.src, "input_files/fsource_corrected_imputed_median.tsv", quote = F, row.names = F, sep="\t")

###
# Food Frequency
###
groupA <- c("fish","white_meat","red_meat")
groupB <- c("rice","legumes")
groupC <- c("tuber","vegetables","fruits")
groupD <- c("eggs","noodles","cooking_oil")
no_var <- c("coffee_tea", "added_sugar", "lard", "sago")
select_var <- c(groupA,"honey",groupB,groupC,groupD,"sugar_g",no_var)

dff <- bio %>%
  left_join(dff, by="sampleid",unmatched = "drop") %>%
  filter(if_any(all_of(select_var), ~ .x != 0))# remove samples if diet is 0 across

nNA <- calculate_missingness(dff, select_var) # %NA
saveRDS(nNA, "out/propNA_dff.rds")
print(nNA)

dff_imp <- dff %>%
  select(sampleid, population, sex, select_var) %>%
  group_by(population) %>%
  # Apply population-specific median imputation across all food variables
  mutate(across(all_of(select_var), ~ replace_na(.x, round(median(.x, na.rm = TRUE), 0)))) %>%
  mutate(sugar_g=round((coffee_tea * added_sugar * 12.5 / 30),1)) %>%
  ungroup()

imputed_samples <- apply(dff[,select_var], 1, function(x) any(is.na(x)))
imputed_samples <- dff$sampleid[imputed_samples]

apply(dff_imp[,-c(1,4)],2, function(x) sum(is.na(x)))

write.table(dff_imp, "input_files/ffq_corrected_imputed_median.tsv", quote = F, row.names = F, sep="\t")

# convert to z
df_freq_Z <- dff_imp %>%
  filter(if_any(all_of(select_var), ~ .x != 0)) %>% # remove samples if diet is 0 across
  mutate(across(all_of(select_var), ~ as.vector(scale(.))))

# PCA
library(vegan)
bio <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>% 
  select(sampleid, population, sex) %>%
  mutate(population=substr(sampleid,1,3))

PCA_var <- c(groupA,groupB,groupC,groupD)
pca_data <- dff %>%
  select(-c(population,sex)) %>%
  right_join(bio, by = "sampleid") %>%
  filter(population != "LDY") %>%
  group_by(population) %>%
  mutate(across(all_of(PCA_var), ~ replace_na(.x, round(median(.x, na.rm = TRUE), 0)))) %>%
  ungroup()%>%
  mutate(across(all_of(PCA_var), ~ as.vector(scale(.)))) %>%
  arrange(sampleid)
dim(pca_data)

select_ids <- pca_data$sampleid
pca_data <- pca_data %>%
  select(all_of(PCA_var)) %>%
  as.matrix()
rownames(pca_data) <- select_ids

ord <- rda(pca_data)
ev <- eigenvals(ord)
total_inertia <- sum(ev)
pc1_explained <- round(100*(ev[1] / total_inertia),1)
pc2_explained <- round(100*(ev[2] / total_inertia),1)
print(c(pc1_explained,pc2_explained))

pca_data <- data.frame(sampleid=rownames(pca_data), pca_data)

# writeout
write.table(data.frame(dff %>% select(-population, -sex)), "input_files/ffq_corrected.tsv", quote = F, row.names = F, sep="\t")
write.table(data.frame(dff_imp %>% select(-population, -sex)), "input_files/ffq_corrected_imputed_median.tsv", quote = F, row.names = F, sep="\t")
write.table(data.frame(df_freq_Z %>% select(-population, -sex)), "input_files/ffq_corrected_imputed_median_Zscore.tsv", quote = F, row.names = F, sep="\t")
write.table(data.frame(pca_data), "input_files/ffq_pca_data.tsv", quote = F, row.names = F, sep="\t")
saveRDS(imputed_samples, "input_files/ffq_imputed_samples.ord")
saveRDS(ord, "input_files/ffq_pca.ord")

# sample counts
dff %>%
  summarise(n=n_distinct(sampleid))

dff %>%
  group_by(population) %>%
  summarise(n=n_distinct(sampleid))

saveRDS(unique(dff$sampleid),"input_files/ffq_sample_list.rds")
saveRDS(unique(df_freq_Z$sampleid),"input_files/ffqZ_sample_list.rds")
saveRDS(unique(df_freq_Z$sampleid),"input_files/ffqZ_sample_list.rds")

# reload
#ffq <- read.delim("input_files/ffq_corrected_imputed.tsv") # older imputation
#ffq <- read.delim("input_files/ffq_corrected_imputed_median.tsv")
