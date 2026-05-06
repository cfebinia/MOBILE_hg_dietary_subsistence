rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyverse)
library(lmerTest)
library(ggeffects)

# colours
population_colors <- readRDS("input_files/mycols.rds")
popcols <- population_colors$col
names(popcols) <- population_colors$population
names(popcols)[8] <- "HUL"

pop.ord <- c("BTU","HUL","APT","TBU","BSP")

df1 <- read.delim("input_files/Selected_GDD_anthro+clinical - compiled.tsv") %>%   rename(w_kg=weight_kg, h_cm=height_cm)

df1 <- df1 %>%
  filter(population %in% pop.ord) %>%
  mutate(BMI=round(w_kg/((h_cm/100)^2),1),
         WHR=round(wc_cm/hc_cm,3),
         WHtR=round(wc_cm/h_cm,3),
         BMI_category=case_when(
           BMI < 18.5 ~ "underweight",
           BMI >= 18.5 & BMI < 23.0 ~ "normal",
           BMI >= 23.0 & BMI < 27.5 ~ "overweight",
           BMI >= 27.5 ~ "obese", # normally 27.5 or 30
           TRUE ~ NA_character_)) %>%
  mutate(
    is_obese = ifelse(BMI >= 27.5,1,0),
    is_overweight = ifelse(BMI >= 23,1,0),
    is_underweight = ifelse(BMI < 18.5,1,0),
    is_hypertension = ifelse(systole>=140 | diastole>=90, 1, 0)
  ) %>% 
  mutate(
    sex=factor(sex, c("M","F")),
    population=factor(population, pop.ord)
  ) %>%
  rename(height_cm=h_cm,
         weight_kg=w_kg,
         bmi=BMI,
         cOb=BMI_category)

# sample number
table(df1$population)
table(df1$population, df1$sex)

# plot
lipids <- c("cholesterol","ldl","hdl","triglyceride","blood_sugar")

ggdata <- df1 %>%
  select(all_of(c("sid", "population", "sex", "bmi", "is_obese", lipids))) %>%
  pivot_longer(
    cols = all_of(c("bmi",lipids)), 
    names_to = "metrics", 
    values_to = "value", 
    values_drop_na = TRUE
  ) %>%
  filter(!is.na(is_obese)) %>%
  mutate(is_obese=ifelse(is_obese==1, "obese","lean"),
         metrics=factor(metrics, c("bmi",lipids)))

ggdata %>% 
  group_by(sex) %>%
  summarise(n=n_distinct(sid))


is_outlier <- ggdata %>%
  group_by(metrics) %>%
  filter(
    value <= quantile(value, 0.001, na.rm = TRUE) | 
      value >= quantile(value, 0.999, na.rm = TRUE)
  ) %>%
  ungroup()

ggdata_M <- ggdata %>% 
  group_by(metrics) %>%
  filter(
    value >= quantile(value, 0.001, na.rm = TRUE) &
      value <= quantile(value, 0.999, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(sex == "M") %>%
  group_by(population) %>%
  filter(n_distinct(sid) >= 5) %>%
  ungroup() %>%
  mutate(metrics=factor(metrics, c("bmi",lipids)))


ggdata_F <- ggdata %>% 
  group_by(metrics) %>%
  filter(
    value >= quantile(value, 0.001, na.rm = TRUE) &
      value <= quantile(value, 0.999, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(sex == "F") %>%
  group_by(population) %>%
  filter(n_distinct(sid) >= 5) %>%
  ungroup() %>%
  mutate(metrics=factor(metrics, c("bmi",lipids)))


# count samples
ggdata_M %>%
  filter(!is.na(value)) %>%
  group_by(population) %>%
  summarise(n=n_distinct(sid)) %>%
  ungroup() %>%
  mutate(xlab=paste0(population,"\n(",n,")")) -> count_M

group_labelM <- count_M$xlab
names(group_labelM) <- count_M$population

ggdata_F %>%
  filter(!is.na(value)) %>%
  group_by(population) %>%
  summarise(n=n_distinct(sid)) %>%
  ungroup() %>%
  mutate(xlab=paste0(population,"\n(",n,")")) -> count_F

group_labelF <- count_F$xlab
names(group_labelF) <- count_F$population


plot_base <-   theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(size=9),
        title=element_text(face="bold", size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.8, "lines"),
        legend.title = element_text(size=6, face="bold"),
        panel.grid.minor.x = element_blank(),
        panel.spacing.y = unit(0.1, "lines"), 
        panel.spacing.x = unit(0.1, "lines"))

plot_M <- ggplot(ggdata_M,
            aes(x = population, y = value, fill = population, colour = population)) +
  facet_wrap(~ metrics, scales = "free_y") +
  geom_point(
    aes(shape = is_obese), alpha = 0.6, size = 2,
    position = position_jitterdodge(seed = 1, jitter.width = 0.25, dodge.width = 0.75, jitter.height = 0)) +
  geom_boxplot(fill=NA, outlier.shape = NA, width=0.75) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_x_discrete(limits = pop.ord, label=group_labelM) +
  scale_colour_manual(values = popcols) +
  labs(title="Men") +
  plot_base 
  
plot_F <- ggplot(ggdata_F,
       aes(x = population, y = value, fill = population, colour = population)) +
  facet_wrap(~ metrics, scales = "free_y") +
  geom_point(
    aes(shape = is_obese), alpha = 0.6, size = 2,
    position = position_jitterdodge(seed = 1, jitter.width = 0.25, dodge.width = 0.75, jitter.height = 0)) +
  geom_boxplot(fill=NA, outlier.shape = NA, width=0.75) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  scale_x_discrete(limits = pop.ord, label=group_labelF) +
  scale_colour_manual(values = popcols) +
  labs(title="Women") +
  plot_base 

wilcox_data <- ggdata %>%
  group_by(population, sex, metrics) %>%
  filter(n_distinct(sid) >= 5) %>%
  ungroup() %>%
  group_by(metrics,sex) %>%
  filter(n_distinct(population) >= 2) %>%
  ungroup()

metrics_list <- as.character(unique(wilcox_data$metrics))

wilcox_res <- list("M" = list(), "F" = list())

for(x in metrics_list) {
  for(s in c("M", "F")) {
    dat <- wilcox_data %>% filter(metrics == x & sex == s)
    
    resp <- dat[["value"]]
    grp <- dat[["population"]]
    
    out_matrix <- pairwise.wilcox.test(
      x = resp, g = grp, 
      p.adjust.method = "fdr", 
      paired = FALSE, exact = FALSE
    )$p.value
    
    out_df <- out_matrix %>% 
      as.table() %>% 
      as.data.frame() %>% 
      filter(!is.na(Freq)) %>% 
      rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
      filter(p_adj < 0.05) %>%
      mutate(p_adj=round(p_adj, 4))
    
    if(nrow(out_df)>0){
      out_final <- data.frame(sex = s, metrics = x, out_df)
      wilcox_res[[s]][[x]] <- out_final}
  }
}

# Combine all list elements into one master data frame
final_results <- bind_rows(lapply(wilcox_res, bind_rows))

wilcox_annot <- final_results %>%
  mutate(metrics=factor(metrics, c("bmi",lipids))) %>%
  mutate(group1=factor(group1, pop.ord),
         group2=factor(group2, pop.ord),
         x1=as.numeric(group1),
         x2=as.numeric(group2),
         xlab=(x2+x1)/2)

wilcox_annot_M <- wilcox_annot %>%
    filter(sex == "M") %>%
    group_by(metrics) %>%
    mutate(
      # Logic: Place labels above the highest point in each facet
      y_base = max(ggdata_M$value[ggdata_M$metrics == unique(metrics)], na.rm = TRUE),
      y_pos = y_base + (row_number() * (y_base * 0.1)),
      label = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    ) %>%
    ungroup()
  
wilcox_annot_F <- wilcox_annot %>%
    filter(sex == "F") %>%
    group_by(metrics) %>%
    mutate(
      # Logic: Place labels above the highest point in each facet
      y_base = max(ggdata_F$value[ggdata_F$metrics == unique(metrics)], na.rm = TRUE),
      y_pos = y_base + (row_number() * (y_base * 0.1)),
      label = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    ) %>%
    ungroup()


plot_M <- plot_M +
  geom_segment(
    data = wilcox_annot_M,
    aes(x = group1, xend = group2, y = y_pos, yend = y_pos),
    inherit.aes = FALSE, colour = "black"
  ) +
  geom_text(
    data = wilcox_annot_M, size=3,
    aes(x = xlab, y = y_pos, label = label),
    inherit.aes = FALSE, colour = "black", vjust = 0
  ) 

plot_M  


plot_F <- plot_F +
    geom_segment(
      data = wilcox_annot_F,
      aes(x = group1, xend = group2, y = y_pos, yend = y_pos),
      inherit.aes = FALSE, colour = "black"
    ) +
    geom_text(
      data = wilcox_annot_F, size=3,
      aes(x = xlab, y = y_pos, label = label),
      inherit.aes = FALSE, colour = "black", vjust = 0
    )


plot_F


##
library(lmerTest)
cordata <- df1 %>%
  select(sid, population, sex, bmi, all_of(lipids)) %>%
  mutate(group = ifelse(population %in% c("BTU", "HUL"), "earlyTransition", "lateTransition")) %>%
  mutate(group = factor(group, levels = c("earlyTransition", "lateTransition")))

subset_M <- cordata %>% 
  filter(sex == "M") %>% 
  filter(if_all(c(bmi, all_of(lipids)), ~ !is.na(.x)))

subset_F <- cordata %>% 
  filter(sex == "F") %>% 
  filter(if_all(c(bmi, all_of(lipids)), ~ !is.na(.x)))

cor_stats <- cordata %>%
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
    rho = sapply(test_obj, function(x) x$estimate),
    p_val = sapply(test_obj, function(x) x$p.value),
    fdr_p = p.adjust(p_val, method = "fdr")
  ) %>%
  select(-test_obj) %>%
  mutate(across(c(rho, p_val, fdr_p), ~ round(.x, 4)))

print(cor_stats)

## 
cordata$obese <- ifelse(cordata$bmi >= 27.5,1,0)

library(lmerTest)
fitdata <- subset(cordata, sex=="M") %>%
  left_join(df1[,c("sid","age")], by="sid") %>%
  mutate(state=ifelse(population %in% c("APT","TBU","BSP"), "lateTransition","earlyTransition")) %>%
  mutate(state=factor(group, c("earlyTransition","lateTransition")))

fit0 <- lmer(cholesterol ~ state + age + (1|population), data=fitdata )
fit <- lmer(cholesterol ~ state + obese + age + (1|population), data=fitdata )
lapply(list(fit0,fit), function(x) summary(x))


fit0 <- lmer(ldl ~ state + age + (1|population), data=fitdata )
fit <- lmer(ldl ~ state + obese + age + (1|population), data=fitdata )
lapply(list(fit0,fit), function(x) summary(x))


fit0 <- lmer(hdl ~ state + age + (1|population), data=fitdata )
fit <- lmer(hdl ~ state + obese + age + (1|population), data=fitdata )
lapply(list(fit0,fit), function(x) summary(x))


fit0 <- lmer(triglyceride ~ state + age + (1|population), data=fitdata )
fit <- lmer(triglyceride ~ state + obese + age + (1|population), data=fitdata )
lapply(list(fit0,fit), function(x) summary(x))


## Hypertension
df_hyp <- df1 %>%
  select(sid, population, sex, is_hypertension) %>%
  mutate(group=ifelse(population %in% c("BTU","HUL"), "earlyTransition","lateTransition"))

df_hyp %>%
  group_by(sex) %>%
  summarise(n=n_distinct(sid))

df_hyp %>%
  filter(sex=="M") %>%
  group_by(population) %>%
  summarise(n=n_distinct(sid)) %>%
  ungroup() %>%
  mutate(xlab=paste0(population," (",n,")")) -> count_Hyp

Hyp_label <- count_Hyp$xlab
names(Hyp_label) <- count_Hyp$population

# 1. Prepare data for plotting (Prevalence Calculation)
plot_data <- df_hyp %>%
  filter(sex == "M") %>%
  group_by(population) %>%
  summarise(
    n_total = length(na.omit(is_hypertension)),
    n_hyper = sum(is_hypertension == 1, na.rm = T),
    prevalence = (n_hyper / n_total) * 100
  ) %>%
  ungroup()

# 2. Perform Fisher's Exact Test
plot_data_group <- df_hyp %>%
  filter(sex == "M") %>%
  group_by(group) %>%
  summarise(
    n_total = length(na.omit(is_hypertension)),
    n_hyper = sum(is_hypertension == 1, na.rm = TRUE),
    prevalence = (n_hyper / n_total) * 100
  ) %>%
  ungroup()

contingency_table <- df_hyp %>%
  filter(!is.na(is_hypertension)) %>%
  filter(sex == "M") %>%
  select(group, is_hypertension) %>%
  table()

fisher_res <- fisher.test(contingency_table)
p_val <- fisher_res$p.value

stat_label <- paste0(
  "Fisher's Exact Test\np = ", 
  ifelse(p_val < 0.001, "< 0.001", round(p_val, 3))
)

pops <- unique(as.character(df_hyp$population[df_hyp$sex == "M"]))
pairs <- combn(pops, 2, simplify = FALSE)

# 2. Iterate through pairs with reinforced data types
pairwise_results <- lapply(pairs, function(p) {
  
  # Logic: Ensure 'population' is character and 'p' is unlisted
  pair_data <- df_hyp %>% 
    filter(sex == "M") %>%
    filter(as.character(population) %in% as.character(p))
  
  # Check if both groups in the pair actually have data
  if(length(unique(pair_data$population)) < 2) return(NULL)
  
  # Create contingency table
  tab <- table(pair_data$population, pair_data$is_hypertension)
  
  # Execute Fisher's Exact Test
  p_val <- fisher.test(tab)$p.value
  
  data.frame(
    group1 = p[1], 
    group2 = p[2], 
    p_raw = p_val, 
    stringsAsFactors = FALSE
  )
})

# 3. Combine and adjust p-values
pairwise_df <- bind_rows(pairwise_results) %>%
  mutate(p_adj = p.adjust(p_raw, method = "fdr"))

plot_data <- df_hyp %>%
  filter(sex == "M") %>%
  group_by(population) %>%
  summarise(prevalence = mean(is_hypertension == 1, na.rm = TRUE) * 100)

annot_df <- pairwise_df %>%
  mutate(
    x1 = as.numeric(factor(group1, levels = pop.ord)),
    x2 = as.numeric(factor(group2, levels = pop.ord)),
    xlab = (x1 + x2) / 2,
    # Logic: Stack brackets above the highest bar
    y_pos = max(plot_data$prevalence) + (row_number() * 5)
  )

popcols2 <- c(rep("forestgreen",2), rep("maroon",3))
names(popcols2) <- c("BTU","HUL","APT","TBU","BSP")

plot_hyp <- ggplot(plot_data, aes(x = population, y = prevalence, fill = population)) +
  geom_bar(stat = "identity", colour = "black", width = 0.7) +
  geom_text(
    aes(label = paste0(round(prevalence, 1), "%")), 
    vjust = -0.5, size = 4
  ) +
  annotate(
    "text", 
    # Logic: -Inf targets the far left of the plotting area
    x = -Inf, y = Inf, 
    label = stat_label, 
    # Logic: hjust = 0 ensures the text starts AT the x-coordinate (left-aligned)
    hjust = -0.1, vjust = 1.5, 
    fontface = "italic", size = 4
  ) +
  labs(
    title = "Hypertension Prevalence in Men",
    y = "Prevalence (%)",
    x = "Population"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_manual(values = popcols2) +
  plot_base+
  scale_x_discrete(limit=pop.ord, label=Hyp_label)+
  scale_y_continuous(
    labels = scales::label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.2))
  )

plot_hyp



# stack
library(patchwork)

narrow_margin <- theme(
  plot.margin = margin(t = 3, r = 3, b = 3, l = 3, unit = "pt"),
  panel.spacing = unit(0.5, "lines")
)

combined_plot <- (plot_M + narrow_margin) / (plot_F + narrow_margin) / (plot_hyp + narrow_margin)


combined_plot <- combined_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(2, 2, 1))& 
  theme(plot.margin = margin(1, 1, 1, 1, "pt"),
        plot.tag = element_text(face = "bold", size = 18))

combined_plot

# save plots

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","04.biometrics","blood_markers_byPop"),".",dev)
  ggsave(plot = combined_plot,filename = outpath, device = dev, width = 20,height = 27, scale=1.2, units = "cm")
}


# count samples
df1 %>% 
  mutate(group=ifelse(population%in%c("BTU","HUL"), "earlyTransition", "lateTransition")) -> df1

length(unique(df1$sid))

df1 %>% 
  group_by(sex) %>%
  summarise(n=n_distinct(sid))

table(df1$group, df1$sex)

df1 %>%
  group_by(group) %>%
  summarise(n=n_distinct(sid))

df1 %>%
  group_by(population) %>%
  summarise(n=n_distinct(sid))

outdata <- merge.data.frame(cordata, df_hyp[,c("sid","is_hypertension")], by="sid")
write.table(outdata, "out/gdd_biomarkers_processed.tsv", quote = F, sep = "\t", row.names = F)
