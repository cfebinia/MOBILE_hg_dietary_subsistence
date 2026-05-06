rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyverse)
library(lmerTest)
library(performance)
library(ggeffects)
library(patchwork)

###############
sex_base <- c("F","M")

# Helper functions
identify_outliers <- function(value_vec, id_vec, coef = 1.5) {
  qs <- quantile(value_vec, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qs[2] - qs[1]
  
  upper_fence <- qs[2] + (coef * iqr)
  lower_fence <- qs[1] - (coef * iqr)
  
  is_outlier <- value_vec > upper_fence | value_vec < lower_fence
  
  return(ifelse(!is.na(is_outlier) & is_outlier, as.character(id_vec), NA))
}

export_mixed_models <- function(model_list, folder_path = ".", base_name = "model_analysis") {
  
  if (!requireNamespace("performance", quietly = TRUE)) stop("Package 'performance' is required.")
  if (!requireNamespace("broom.mixed", quietly = TRUE)) stop("Package 'broom.mixed' is required.")
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required.")
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  if (is.null(names(model_list))) {
    names(model_list) <- paste0("Model_", seq_along(model_list))
  }
  
  estimates_path <- file.path(folder_path, paste0(base_name, "_estimates.tsv"))
  fit_scores_path <- file.path(folder_path, paste0(base_name, "_fit_indices.tsv"))
  
  fixed_effects <- do.call(rbind, lapply(names(model_list), function(m_name) {
    res <- broom.mixed::tidy(model_list[[m_name]], effects = "fixed", conf.int = TRUE)
    stat_label <- if (inherits(model_list[[m_name]], "glmerMod")) "z_score" else "t_score"
    colnames(res)[colnames(res) == "statistic"] <- stat_label
    
    res$p_adj <- p.adjust(as.numeric(res$p.value), method = "fdr")
    
    res$model_id <- m_name
    return(res)
  }))
  
  select_cols <- intersect(c("estimate", "std.error", "t_score", "z_score", "df", "conf.low", "conf.high"), colnames(fixed_effects))
  fixed_effects[, select_cols] <- apply(fixed_effects[, select_cols], 2, function(x) round(as.numeric(x), 2))
  
  fit_scores <- do.call(rbind, lapply(names(model_list), function(m_name) {
    m <- model_list[[m_name]]
    
    r2_res <- tryCatch(performance::r2(m), error = function(e) NULL)
    icc_res <- tryCatch(performance::icc(m), error = function(e) NULL)
    is_singular <- lme4::isSingular(m)
    
    # Check if results are lists and contain the expected elements
    m_r2 <- if (is.list(r2_res) && "R2_marginal" %in% names(r2_res)) round(r2_res$R2_marginal, 2) else 0
    c_r2 <- if (is.list(r2_res) && "R2_conditional" %in% names(r2_res)) round(r2_res$R2_conditional, 2) else 0
    adj_icc <- if (is.list(icc_res) && "ICC_adjusted" %in% names(icc_res)) round(icc_res$ICC_adjusted, 3) else 0
    
    data.frame(
      model_id = m_name,
      is_singular = is_singular,
      marginal_R2 = m_r2,
      conditional_R2 = c_r2,
      variance_random = adj_icc,
      AIC = round(AIC(m), 0),
      n_obs = nobs(m)
    )
  }))
  
  write.table(fixed_effects, estimates_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(fit_scores, fit_scores_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("Success: Files saved to ", folder_path)
}

# Colours
population_colors <- readRDS("input_files/mycols.rds")
popcols <- population_colors$col
names(popcols) <- population_colors$population

groupCols <- c("lateTransition" = "maroon", "earlyTransition" = "forestgreen")
groupShape <- c("lateTransition" = 15, "earlyTransition" = 16) # this one need flat colours
groupshape2 <- c("lateTransition" = 22, "earlyTransition" = 21) # this one can have borders


# populations
pop.ord <- c("BTU","ORT","ORS","APT","TBU","BSP") # LDY omitted
earlyTransition <- c("BTU","ORS","ORT")
lateTransition <- c("APT","TBU","BSP")
pop.ord <- c("BTU","ORT","ORS","APT","TBU","BSP")
pop.select <- c("BTU","ORT","ORS","APT","TBU","BSP") # sometimes this can be different
HGgroup <- data.frame(population=c(earlyTransition, lateTransition), 
                      group=c(rep("earlyTransition", length(earlyTransition)),
                              rep("lateTransition", length(lateTransition))))
HGgroup$group <- factor(HGgroup$group, c("earlyTransition","lateTransition"))
HGgroup$ethnic <- ifelse(HGgroup$population %in% c("ORT","ORS"), "OrangRimba", "Others")
HGgroup$ethnic <- factor(HGgroup$ethnic, c("Others","OrangRimba"))

###############
# LOAD DATA
# this "RAW" data has been corrected (e.g. converted when inch was recorded as cm)
# See this table for the list of corrections: "2023_MOBILE_AnthropometryData - Corrections"
#df3 <- read.delim("input_files/2023_MOBILE_AnthropometryData(RAW).tsv") %>%
df0 <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>%
  mutate(population = substr(sampleid, 1,3)) %>%
  left_join(HGgroup, by="population")

# REPORTED STATS
df3 <- df0 %>% 
  filter(population %in% pop.select) %>%
  mutate(sex = factor(sex, sex_base),
         population = factor(population, levels = pop.select))

# IMPUTATION and FILTERING
# normality test
apply(df3[,c("age", "w_kg", "h_cm", "wc_cm", "hc_cm")],2, function(x) shapiro.test(na.omit(x[df3$sex=="M"])))
apply(df3[,c("age", "w_kg", "h_cm", "wc_cm", "hc_cm")],2, function(x) shapiro.test(na.omit(x[df3$sex=="F"])))
message("Only height is normally distributed. We'll universally impute with median (by population&sex).")

# record imputed sampleids
imputed_sid <- unique(unlist(apply(df3[,c("age", "w_kg", "h_cm", "wc_cm", "hc_cm")],2, function(x) df3$sampleid[is.na(x)])))
length(imputed_sid)
table(substr(imputed_sid,1,3))

# impute
df3_imputed <- df3 %>%
  mutate(population=substr(sampleid,1,3)) %>%
  filter(population %in% pop.ord) %>%
  group_by(population, sex) %>%
  mutate(across(
    .cols = c(age, w_kg, h_cm, wc_cm, hc_cm),
    .fns = ~ ifelse(is.na(.), median(., na.rm = TRUE), .)
  )) %>%
  ungroup() %>%
  select(sampleid, population, sex, age, w_kg, h_cm, wc_cm, hc_cm)

# TRESHOLDING & RENAME COLUMNS
df3_imputed <- df3_imputed %>%
  mutate(BMI=round(w_kg/((h_cm/100)^2),1),
       WHR=round(wc_cm/hc_cm,2),
       WHtR=round(wc_cm/h_cm,2),
       BMI_category=case_when(
         BMI < 18.5 ~ "underweight",
         BMI >= 18.5 & BMI < 23.0 ~ "normal",
         BMI >= 23.0 & BMI < 27.5 ~ "overweight",
         BMI >= 27.5 ~ "obese", # normally 27.5 or 30
         TRUE ~ NA_character_)) %>%
  mutate(
    is_obese = ifelse(BMI >= 27.5,1,0),
    is_overweight = ifelse(BMI >= 23,1,0),
    is_underweight = ifelse(BMI < 18.5,1,0)
    ) %>% 
  rename(sid=sampleid, 
         height_cm=h_cm,
         weight_kg=w_kg,
         bmi=BMI,
         cOb=BMI_category)

apply(df3, 2, function(x) sum(is.na(x))) 
apply(df3_imputed, 2, function(x) sum(is.na(x))) 

# SUBSETTING and MERGE WITH GROUP DATA
dfob <- df3_imputed %>%
  filter(population %in% pop.select) %>%
  left_join(HGgroup, by="population") # %>%
  # filter(age<61 & age>=18)

dfob %>% 
  group_by(population) %>%
  summarise(n=n_distinct(sid),
            M=sum(sex=="M"),
            F=sum(sex=="F"))

###########
## MALNUTRITION
##########
fitdata <- dfob[,c("sid","population","age","sex", "cOb", "is_overweight","is_underweight" )] %>%
  left_join(HGgroup, by="population") %>%
  mutate(
    population = factor(population, pop.ord),
    cOb = factor(cOb, c("underweight","normal","overweight","obese")))

legend_labs <- c("Underweight (BMI<18.5)", "Normal", "Overweight (BMI:23 - 27.4)", "Obese (BMI≥27.5)")

fitdata_labeled <- fitdata %>%
  group_by(sex, population) %>%
  mutate(n = n(),
         pop_sex_label = paste0(population, "\n(", n, ")")) %>%
  ungroup() %>%
  arrange(sex) %>% arrange(population) %>% 
  mutate(pop_sex_label=factor(pop_sex_label, unique(pop_sex_label)))

fitdata_labeled <- fitdata_labeled %>% mutate(sex=factor(sex, c("M","F")))

plot_outcomes <- ggplot(fitdata_labeled, aes(x = pop_sex_label, fill = cOb)) +
  geom_bar(colour = "black", position = "fill", width = 1) +
  labs(
    title = "Weight Status Proportions by Population and Sex",
    x = "Population",
    y = "Proportion",
    fill = "Weight Status"
  ) +
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  theme_bw(base_size = 10) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "underweight" = "#4575b4",
    "normal"      = "#e0e0e0",
    "overweight"  = "#f98e00",
    "obese"       = "#d73027"
  ), labels = legend_labs) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(face="bold"),
    legend.position = "bottom", 
    legend.text = element_text(size = 9), 
    legend.key.size = unit(4,"mm"),
    legend.title = element_blank(),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_rect(fill = "#ffb380"),
    strip.text = element_text(face = "bold", size=9)
  )

plot_outcomes

figout=file.path(getwd(),"fig_out","04.biometrics/bmi","weightStats_preval")
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = plot_outcomes,filename = outpath, device = dev, width = 8.5,height = 3.5,dpi = 300, scale=1)
}

fitdata_labeled <- fitdata %>%
  group_by(sex, group) %>%
  mutate(n = n(),
         group_sex_label = paste0(group, "\n(", n, ")")) %>%
  ungroup() %>%
  arrange(sex) %>% arrange(group) %>% 
  mutate(group_sex_label = gsub("earlyT","Early-t",group_sex_label)) %>%
  mutate(group_sex_label = gsub("lateT","Late-t",group_sex_label)) %>%
  mutate(group_sex_label=factor(group_sex_label, unique(group_sex_label)))

fitdata_labeled <- fitdata_labeled %>% mutate(sex=factor(sex, c("M","F")))
plot_outcomes2 <- ggplot(fitdata_labeled, aes(x = group_sex_label, fill = cOb)) +
  geom_bar(colour = "black", position = "fill", width = 1) +
  labs(
    title = "Weight Status Proportions by Population and Sex",
    x = "Population",
    y = "Proportion",
    fill = "Weight Status"
  ) +
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  theme_bw(base_size = 10) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "underweight" = "#4575b4",
    "normal"      = "#e0e0e0",
    "overweight"  = "#f98e00",
    "obese"       = "#d73027"
  ), labels = legend_labs) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(face="bold"),
    legend.position = "bottom", 
    legend.text = element_text(size = 9), 
    legend.key.size = unit(4,"mm"),
    legend.title = element_blank(),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_rect(fill = "#ffb380"),
    strip.text = element_text(face = "bold", size=9)
  )

plot_outcomes2

figout=file.path(getwd(),"fig_out","04.biometrics/bmi","weightStats_preval_group")
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = plot_outcomes2,filename = outpath, device = dev, width = 3.5,height = 3.5,dpi = 300, scale=1)
}

###############
### BMI x Adiposity
################
fitdata <- dfob[,c("sid","population","age","sex", "cOb", "is_overweight","is_underweight", "bmi", "WHtR")] %>%
  left_join(HGgroup, by="population") %>%
  filter(population!="LDY") %>%
  mutate(
    sex=factor(sex, sex_base),
    population = factor(population, pop.ord),
    cOb = factor(cOb, c("underweight","normal","overweight","obese")))

fit0 <- lmer(WHtR ~ bmi + age + sex + (1|population), data=fitdata)
fit1 <- lmer(WHtR ~ bmi + sex + age + ethnic + (1|population), data=fitdata)
fit2 <- lmer(WHtR ~ bmi + age + sex + group + (1|population), data=fitdata)
fit3 <- lmer(WHtR ~ bmi + age + sex + ethnic + group + (1|population), data=fitdata)
fit4 <- lmer(WHtR ~ bmi + age + sex + group * sex + (1|population), data=fitdata)

fit0BMI <- lmer(bmi ~ age * sex + ethnic + group * sex + (1|population), data=fitdata)

model_list <- list(Model0=fit0, Model1=fit3)
lapply(model_list, function(x) summary(x))

export_mixed_models(
  model_list,
  folder_path = "out/", 
  base_name = "LMM_Results_WHTR_vs_BMI")

final_model <- fit0
summary(fit0)

r2_values <- model_performance(final_model)
marginal_r2 <- round(r2_values$R2_marginal, 2)

saveRDS(final_model, "out/adiposity_lmer_final.rds")

#plot
fitdata$pred <- predict(
  final_model, 
  newdata = transform(fitdata, age = mean(fitdata$age, na.rm = TRUE)), 
  re.form = NA
)

# 2. Create the plot
fitdata <- fitdata %>% mutate(sex=factor(sex, c("M","F")))
plot_WHTR <- ggplot(fitdata, aes(y = bmi, x = WHtR, color = group, shape=group)) + 
  geom_vline(xintercept = 0.5, linetype = "dashed", col = "grey") + 
  geom_hline(yintercept = 23, linetype = "dashed", col = "grey") + 
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  #geom_text(size=2, col="grey20", aes(label=sid), nudge_x = 1)+
  #geom_smooth(aes(group=group), method = "lm", se = FALSE, linewidth = 1, col="white") +
  #geom_smooth(method = "lm", linewidth = 1, se = FALSE) +
  geom_line(aes(x = pred, group=interaction(group,age)), linewidth = 1, col="black") +
  geom_point(alpha = 0.4, size=2) +
  scale_color_manual(values = groupCols)+
  scale_shape_manual(values = groupShape)+
  labs(
    y = "BMI (kg/m²)", 
    x = "Waist (cm) / Height (cm)", 
    title = "Central Adiposity by BMI and Sex",
    subtitle = paste0("Marginal R-squared = ", marginal_r2),
    caption = "—— LMM: WHtR ~ BMI + age + sex + (1|community)") +
  theme_bw(base_size = 10) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    plot.caption = element_text(face = "plain", size=9),
    plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(face="bold"),
    axis.text = element_text(size = 9),
    legend.position = "bottom", 
    legend.text = element_text(size = 9), 
    legend.key.size = unit(4,"mm"),
    legend.title = element_blank(),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_rect(fill = "#ffb380"),
    strip.text = element_text(face = "bold", size=9)
  )

plot_WHTR

figout="fig_out/04.biometrics/FigureS4_raw"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = plot_WHTR,filename = outpath, device = dev, width = 11,height = 7,dpi = 300, scale=0.6)
}


###############
### BMI x Diet
################
library(vegan)

# load
ord <- readRDS("input_files/ffq_pca.ord")
dff <- read.delim("input_files/ffq_pca_data.tsv")

envfit <- envfit(ord, dff, perm=999)

select_var <- colnames(dff)[-1]

dff <- df3 %>% 
  filter(population != "LDY") %>%
  mutate(group = case_when(
    population %in% c("BTU", "ORT", "ORS") ~ "earlyTransition",
    population %in% c("APT", "TBU", "BSP") ~ "lateTransition",
    TRUE ~ NA_character_
  )) %>%
  select(sampleid, population, group, age, sex) %>%
  left_join(dff, unmatched = "drop") %>%
  filter(if_any(all_of(colnames(dff)[-1]), ~ .x != 0))

no_diet_ids <- unique(df3$sampleid[!df3$sampleid %in% dff$sampleid])

# group definitions
groupcols <- c(earlyTransition="forestgreen",
               lateTransition="maroon",
               Agri="grey50")

groupshapes <- c(earlyTransition=21,
                 lateTransition=22,
                 Agri=24)

foodcols <- c("Group A" = "forestgreen",
              "Group B"="maroon",
              "Group C"="gold2", 
              "Other"="grey")

groupA <- c("fish","white_meat","red_meat","honey")
groupB <- c("rice","legumes","tuber","vegetables","fruits")
groupC <- c("eggs","noodles","cooking_oil")
select_var <- c(groupA,groupB,groupC)

#########
# get scores to plot
species_scores <- as.data.frame(scores(ord, display = "species")) *0.6
write.table(data.frame(species_scores %>%
                         mutate(food=rownames(species_scores))),
            "fig_out/04.biometrics/species_scores.tsv",quote = F,sep = "\t",row.names = F)

site_scores <- as.data.frame(scores(ord, display = "sites"))
site_scores$sampleid <- rownames(site_scores)
site_scores <- dff[,c("sampleid","population","group")] %>%
  left_join(site_scores, by="sampleid") %>%
  filter(! sampleid %in% no_diet_ids)


df_env <- as.data.frame(scores(envfit, display = "vectors"))
df_env$r2 <- envfit$vectors$rsq
df_env$p <- envfit$vectors$pvals
df_env$var <- rownames(df_env)


limit_val <- max(
  abs(site_scores$PC1), abs(site_scores$PC2),
  abs(species_scores$PC1), abs(species_scores$PC2)
) * 1.15

# Variance Explained
ev <- eigenvals(ord)
total_inertia <- sum(ev)
pc1_explained <- round(100*(ev[1] / total_inertia),1)
pc2_explained <- round(100*(ev[2] / total_inertia),1)
xlab <- paste0("PC1 (",pc1_explained,"%)")
ylab <- paste0("PC2 (",pc2_explained,"%)")

# PC loadings
loadings <- scores(ord, display = "species", choices = c(1, 2), scaling = 0)
loadings <- data.frame(
  Variable = rownames(loadings),
  PC1 = loadings[, 1],
  PC2 = loadings[, 2],
  Abs_PC1 = abs(loadings[, 1]),
  Abs_PC2 = abs(loadings[, 2])
)
pc1_ranked <- loadings[order(loadings$Abs_PC1, decreasing = TRUE), c("Variable", "PC1")]
pc2_ranked <- loadings[order(loadings$Abs_PC2, decreasing = TRUE), c("Variable", "PC2")]

outfile=paste0(file.path(getwd(),"fig_out/04.biometrics/bmi/","PCA_FFQ_Zscore_loadings"),"_loadings.txt")
sink(outfile)
print(pc1_ranked)
print(pc2_ranked)
sink()

loadings_df <- pc1_ranked %>% left_join(pc2_ranked, by="Variable") %>%
  mutate(food_group=case_when(
    Variable %in% groupA ~ "Group A",
    Variable %in% groupB ~ "Group B",
    Variable %in% groupC ~ "Group C",
    TRUE ~ "Other"
  )) %>%
  mutate(food_group = factor(food_group, c("Group A","Group B", "Group C", "Group D", "Other")))

pc1 <- loadings_df %>%
  mutate(Variable = factor(Variable, levels = Variable[order(PC1, decreasing = TRUE)]),
         sign = PC1 > 0) %>%
  ggplot() +  
  geom_col(aes(x = PC1, y = Variable, fill = food_group)) +
  scale_fill_manual(
    values = foodcols,
    labels = c(
      "Group A" = "Animal-based",
      "Group B" = "Plant-based",
      "Group C" = "Market"
    )
  ) +
  labs(x = NULL, y = NULL, 
       title = "PC1 Loadings",
       subtitle = paste0("PC1 (",pc1_explained,"% Variance Explained)")) +
  theme_bw(base_size = 10) +
  theme(
    plot.margin = margin(5, 5, 5, 5, "pt"),
    title = element_text(face="bold"),
    plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom", 
    legend.text = element_text(size = 7), 
    legend.key.size = unit(4,"mm"),
    legend.title = element_blank())
pc1

figout=file.path(getwd(),"fig_out","04.biometrics/bmi","PC1_loadings")
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,"_diet.",dev)
  ggsave(plot = pc1,filename = outpath, device = dev, width = 4.2,height = 6)
}

pc2 <- loadings_df %>%
  mutate(Variable = factor(Variable, levels = Variable[order(PC2, decreasing = TRUE)]),
         sign = PC2 > 0) %>%
  ggplot() +  
  geom_col(aes(x = PC2, y = Variable, fill = food_group)) +
  scale_fill_manual(
    values = foodcols,
    labels = c(
      "Group A" = "Animal-based",
      "Group B" = "Plant-based",
      "Group C" = "Market"
    )
  ) +
  labs(x = NULL, y = NULL, 
       title = "PC2 Loadings",
       subtitle = paste0("PC2 (",pc2_explained,"% Variance Explained)")) +
  theme_bw(base_size = 10) +
  theme(
    plot.margin = margin(5, 5, 5, 5, "pt"),
    plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40"),
    title = element_text(face="bold"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom", 
    legend.text = element_text(size = 7), 
    legend.key.size = unit(4,"mm"),
    legend.title = element_blank())
pc2

figout=file.path(getwd(),"fig_out","04.biometrics/bmi","PC2_loadings")
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,"_diet.",dev)
  ggsave(plot = pc2,filename = outpath, device = dev, width = 4.2,height = 6)
}

# plot PCA
p <- ggplot() +  theme_bw(base_size = 14) +
  theme(plot.margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"),
        panel.grid = element_blank())+
  geom_vline(xintercept = 0, linetype="dashed", col="grey")+ 
  geom_hline(yintercept = 0, linetype="dashed", col="grey")+
  #scale_x_continuous(limits = c(-limit_val, limit_val)) +
  #scale_y_continuous(limits = c(-limit_val, limit_val)) +
  coord_fixed(ratio = 1)+
  labs(x = xlab, y = ylab, title = paste("PCA -", "FFQ_Zscore"))

for(sc in c("population","lifestyle")){
  if(sc == "population"){
    q <- p + 
      geom_point(data = site_scores, alpha=0.6, 
                 aes(x = PC1, y = PC2, col=population, fill=population, shape=group), size=2)+
      #geom_segment(data = species_scores, 
      #             aes(x = 0, y = 0, xend = PC1, yend = PC2),
      #             arrow = arrow(length = unit(0.2, "cm")), colour = "red") +
      geom_segment(data = df_env, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2),
                   arrow = arrow(length = unit(0.2, "cm")), 
                   color = "red", linewidth = 0.7) +
      geom_text(data = species_scores, size=3,
                aes(x = PC1, y = PC2, label = rownames(species_scores)), 
                colour = "red", vjust = -0.5)+
      scale_shape_manual(values = groupshapes) +
      scale_colour_manual(values = popcols) +
      scale_fill_manual(values = popcols)
    
    figout=file.path(getwd(),"fig_out","04.biometrics/bmi","PCA_FFQ_Zscore")
    for(dev in c("svg","png","pdf")){
      outpath=paste0(figout,"_population.",dev)
      ggsave(plot = q,filename = outpath, device = dev, width = 6,height = 6)
    }
    pca_population <- q
    
  } else {
    q <- p + 
      geom_point(data = site_scores, alpha=0.6, 
                 aes(x = PC1, y = PC2, col=group, fill=group, shape=group), size=2)+
      geom_segment(data = species_scores, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2),
                   arrow = arrow(length = unit(0.2, "cm")), colour = "red") +
      geom_text(data = species_scores, size=3,
                aes(x = PC1, y = PC2, label = rownames(species_scores)), 
                colour = "red", vjust = -0.5) +
      scale_colour_manual(values = groupcols) +
      scale_fill_manual(values = groupcols) +
      scale_shape_manual(values = groupshapes)
    
    print(q)
    
    figout=file.path(getwd(),"fig_out","04.biometrics/bmi","PCA_FFQ_Zscore")
    for(dev in c("svg","png","pdf")){
      outpath=paste0(figout,"_group.",dev)
      ggsave(plot = q,filename = outpath, device = dev, width = 6,height = 6)
    }
    pca_group <- q
    
    }}

outfile=paste0(file.path(getwd(),"fig_out","04.biometrics/bmi","PCA_FFQ_Zscore"),"_scores.tsv")
write.table(site_scores, file=outfile,quote = F,row.names = F,sep = "\t")

### lmer
pc_scores <- site_scores %>%
  select(sampleid, PC1, PC2)

fitdata0 <- dfob[,c("sid","population","age","sex", "cOb", "is_overweight","is_underweight", "bmi", "WHtR")] %>%
  rename(sampleid=sid) %>%
  left_join(HGgroup, by="population") %>%
  filter(population!="LDY") %>%
  mutate(
    sex=factor(sex, sex_base),
    population = factor(population, pop.ord))

fitdata <- fitdata0 %>%
  left_join(pc_scores, by="sampleid") %>%
  filter(!is.na(PC1)) 

dim(fitdata0)
dim(fitdata)
table(fitdata$population)

fitdata$sex <- factor(fitdata$sex, c("F","M"))
fit0 <- lmer(bmi ~ group + sex + age + group:sex + sex: age + + (1|population), data=fitdata0) 
fit2 <- lmer(bmi ~ PC1 + PC2 + sex + age + (1|population), data=fitdata)
fit3 <- lmer(bmi ~ PC1 + PC2 + sex + age + PC1:sex + (1|population), data=fitdata)

fit_wtht <- lmer(WHtR ~ bmi + sex + age + (1|population), data=fitdata0)

model_list <- list(Model0 = fit0, Model2=fit_wtht, Model3=fit2, Model4=fit3) 
lapply(model_list, function(x) summary(x))
dir.create("fig_out/04.biometrics/bmi/", showWarnings = F, recursive = T)
export_mixed_models(
  model_list,
  folder_path = "out/", 
  base_name = "BMI_LMM_Results_FFQ")

cor.test(fitdata$bmi, fitdata$PC1)
cor.test(fitdata$bmi, fitdata$PC2)

# write table
write.table(fitdata, file = "fig_out/04.biometrics/bmi/bmi_lmerPCA_data.tsv", sep = "\t",quote = F,row.names = F)
write.table(fitdata0, file="fig_out/04.biometrics/lmer_bmi_data.tsv",quote = F,row.names = F,sep = "\t")

#####################
## Plot BMI vs PCs
ggdata <- fitdata %>% 
  select(bmi, sex, group, PC1) 

bmi_model <- fit3
library(performance)
r2_values <- model_performance(bmi_model)
marginal_r2 <- round(r2_values$R2_marginal, 2)
icall <- paste(
  "-- LMM: BMI ~",
  paste0(colnames(coefficients(bmi_model)[[1]])[-1], collapse = " + "),
  "+ (1|community)")

pred_full <- ggpredict(bmi_model, 
                       terms = c("PC1", "sex"),
                       condition = c(age = mean(fitdata$age, na.rm = TRUE))) %>%
  as.data.frame() %>%
  rename(PC1 = x, 
         sex = group) %>%
  mutate(sex = factor(sex, sex_base))

ggdata <- ggdata %>% mutate(sex=factor(sex, c("M","F")))

pc1_bmi <- ggplot(ggdata, aes(y = bmi, x = PC1, color = group)) + 
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  geom_hline(yintercept = 23, linetype = "dashed", col = "grey") + 
  geom_point(size=2, alpha = 0.4, aes(shape=group, group=interaction(group, sex))) +
  #geom_smooth(aes(group=group), method = "lm", se = FALSE, linewidth = 1.2, col="white") +
  #geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  # geom_smooth(col="black", method = "lm", se = FALSE, linewidth = 1, linetype="dotted") +
  geom_line(data = pred_full, 
            aes(y = predicted, x = PC1), 
            color = "black", 
            linewidth = 1)+
  scale_color_manual(values = groupCols)+
  scale_shape_manual(values = groupShape)+
  labs(
    x = "PC1 Score", 
    y = "BMI (kg/m2)", 
    title = "BMI Associations with Food Intake (PC1)",
    subtitle = paste0("Marginal R-squared = ", marginal_r2),
    caption = icall
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    plot.caption = element_text(face = "plain", size=9),
    plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40"),
    title=element_text(face="bold"),
    axis.text = element_text(size=10),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(face="bold"),
    legend.position = "bottom", 
    legend.text = element_text(size = 9), 
    legend.key.size = unit(4,"mm"),
    legend.title = element_blank(),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_rect(fill = "#ffb380"),
    strip.text = element_text(face = "bold", size=9)
  )

pc1_bmi

figout="fig_out/04.biometrics/bmi/bmi_vs_PC1"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = pc1_bmi,filename = outpath, device = dev, width = 11,height = 5,dpi = 300, scale=0.8)
}


## PC2
ggdata <- fitdata %>% 
  select(bmi, sex, ethnic, group, PC2)

bmi_model <- fit3
r2_values <- model_performance(bmi_model)
marginal_r2 <- round(r2_values$R2_marginal, 2)

pred_full <- ggpredict(bmi_model, 
                       terms = c("PC2", "sex"),
                       condition = c(age = mean(fitdata$age, na.rm = TRUE))) %>%
  as.data.frame() %>%
  rename(PC2 = x, 
         sex = group) %>%
  mutate(sex = factor(sex, sex_base))

ggdata <- ggdata %>% mutate(sex=factor(sex, c("M","F")))
pc2_bmi <- ggplot(ggdata, aes(y = bmi, x = PC2, color = group)) + 
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  geom_hline(yintercept = 23, linetype = "dashed", col = "grey") + 
  geom_point(size=2, alpha = 0.4, aes(shape=group, group=interaction(group, sex))) +
  #geom_smooth(aes(group=group), method = "lm", se = FALSE, linewidth = 1.2, col="white") +
  #geom_smooth(col="black", method = "lm", se = FALSE, linewidth = 1, linetype="dotted") +
  geom_line(data = pred_full, 
            aes(y = predicted, x = PC2), group=1, 
            color = "black", 
            linewidth = 1)+
  scale_color_manual(values = groupCols)+
  scale_shape_manual(values = groupShape)+
  labs(
    x = "PC2 Score", 
    y = "BMI (kg/m2)", 
    title = "BMI Associations with Food Intake (PC2)",
    subtitle = paste0("Marginal R-squared = ", marginal_r2),
    caption = icall
  )  +
  theme_bw(base_size = 10) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    plot.caption = element_text(face = "plain", size=9),
    plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40"),
    title=element_text(face="bold"),
    axis.text = element_text(size=10),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(face="bold"),
    legend.position = "bottom", 
    legend.text = element_text(size = 9), 
    legend.title = element_blank(),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_rect(fill = "#ffb380"),
    strip.text = element_text(face = "bold", size=9)
  )

pc2_bmi

figout="fig_out/04.biometrics/bmi/bmi_vs_PC2"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = pc2_bmi,filename = outpath, device = dev, width = 11,height = 5,dpi = 300, scale=0.8)
}

# Save Plots for Figure 3
figout=file.path(getwd(),"fig_out","04.biometrics/bmi","PCA_FFQ_Zscore")
saveRDS(pca_population, paste0(figout,"_group.rds"))
saveRDS(pca_group, paste0(figout,"_group.rds"))

figout=file.path(getwd(),"fig_out","03.dietary_assesment","Fig3B.rds")
plot_3B <- pca_group + 
  theme(legend.position = "bottom", title=element_text(size=8, face="bold"))
saveRDS(plot_3B, figout)

# Save Plots for Figure 4
plot_A1 <- plot_outcomes2 + labs(title = NULL, x=NULL)
plot_A2 <- plot_outcomes + labs(title = NULL, x=NULL, y=NULL) + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())
plot_B <- pc1_bmi + theme(title=element_text(size=8, face="bold"))
plot_C <- pc1 + theme(title=element_text(size=8), face="bold")
plot_D <- pc2_bmi + theme(title=element_text(size=8, face="bold"))
plot_E <- pc2 + theme(title=element_text(size=8, face="bold"))


to_save=list(
  Fig4A1=plot_A1,
  Fig4A2=plot_A2,
  Fig4B=plot_B,
  Fig4C=plot_C,
  Fig4D=plot_D,
  Fig4E=plot_E
)

for(fig in names(to_save)){
  outfile=paste0(fig,".rds")
  saveRDS(to_save[[fig]], file = file.path("fig_out/04.biometrics/",outfile))
}

# sample list
length(dfob$sid)

dfob %>% 
  group_by(group) %>%
  summarise(n=n_distinct(sid))

dfob %>% 
  group_by(sex) %>%
  summarise(n=n_distinct(sid))

dfob %>% 
  mutate(population=factor(population, pop.ord)) %>%
  group_by(population) %>%
  summarise(n=n_distinct(sid))


saveRDS(unique(dfob$sid),file = "input_files/obesity_sample_list.rds")

# save dataset
df0 <- df0 %>% 
  mutate(BMI=round(w_kg/((h_cm/100)^2),1),
               WHR=round(wc_cm/hc_cm,2),
               WHtR=round(wc_cm/h_cm,2),
               BMI_category=case_when(
                 BMI < 18.5 ~ "underweight",
                 BMI >= 18.5 & BMI < 23.0 ~ "normal",
                 BMI >= 23.0 & BMI < 27.5 ~ "overweight",
                 BMI >= 27.5 ~ "obese", # normally 27.5 or 30
                 TRUE ~ NA_character_)) %>%
  mutate(
    is_obese = ifelse(BMI >= 27.5,1,0),
    is_overweight = ifelse(BMI >= 23,1,0),
    is_underweight = ifelse(BMI < 18.5,1,0)
  ) %>% 
  rename(sid=sampleid, 
         height_cm=h_cm,
         weight_kg=w_kg,
         bmi=BMI,
         cOb=BMI_category)

df0 <- df0 %>%
  mutate(group = case_when(
    population %in% c("BTU", "ORT", "ORS") ~ "earlyTransition",
    population %in% c("APT", "TBU", "BSP") ~ "lateTransition",
    population == "LDY" ~ "Agriculture",
    TRUE ~ NA_character_
  ))

dfob <- dfob %>%
  mutate(group = case_when(
    population %in% c("BTU", "ORT", "ORS") ~ "earlyTransition",
    population %in% c("APT", "TBU", "BSP") ~ "lateTransition",
    population == "LDY" ~ "Agriculture",
    TRUE ~ NA_character_
  ))

write.table(df0, "out/anthrop_data_processed.tsv", quote = F, sep = "\t", row.names = F)  
write.table(dfob, "out/anthrop_data_median_imputed_processed.tsv", quote = F, sep = "\t", row.names = F)