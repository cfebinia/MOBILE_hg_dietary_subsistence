rm(list=ls())
library(tidyverse)
library(ggplot2)
library(lmerTest)
library(ggeffects)
library(broom.mixed)

#################
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

custom_theme <- theme_bw(base_size = 10) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = unit(0.2, "mm"), colour = "grey", linetype = "dashed"),
    panel.spacing.y = unit(0.2, "lines"), 
    panel.spacing.x = unit(0.2, "lines"),
    strip.background = element_rect(fill = "#ffb380"),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right", 
    legend.text = element_text(size = 9), 
    legend.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10)
  )

# facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
# guides(fill = guide_legend(ncol = 1), col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))+

#################
# Colours
population_colors <- readRDS("input_files/mycols.rds")
popcols <- population_colors$col
names(popcols) <- population_colors$population

groupCols <- c("lateTransition" = "maroon", "earlyTransition" = "forestgreen", "Agricultural"="grey40")
groupShape <- c("lateTransition" = 15, "earlyTransition" = 16, "Agricultural"=17) # this one need flat colours
groupshpae2 <- c("lateTransition" = 22, "earlyTransition" = 21, "Agricultural"=24) # this one can have borders

# populations
pop.ord <- c("BTU","ORT","ORS","APT","TBU","BSP", "LDY")
earlyTransition <- c("BTU","ORS","ORT")
lateTransition <- c("APT","TBU","BSP")
Agricultural <- "LDY"
pop.select <- c("BTU","ORT","ORS","APT","TBU","BSP", "LDY")
HGgroup <- data.frame(population=c(earlyTransition, lateTransition, Agricultural), 
                      group=c(rep("earlyTransition", length(earlyTransition)),
                              rep("lateTransition", length(lateTransition)),
                              "Agricultural"))
HGgroup$group <- factor(HGgroup$group, c("earlyTransition","lateTransition","Agricultural"))
HGgroup$ethnic <- ifelse(HGgroup$population %in% c("ORT","ORS"), "OrangRimba", "Others")
HGgroup$ethnic <- factor(HGgroup$ethnic, c("Others","OrangRimba"))

#################
## load data
# df_obesity <- read.delim("input_files/2023_MOBILE_AnthropometryData(RAW).tsv") %>%
df_obesity <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>%
  mutate(population = substr(sampleid, 1,3)) %>%
  mutate(bmi=w_kg/((h_cm/100)^2)) %>%
  filter(population %in% pop.select) %>%
  select(sampleid, population, age, sex, bmi) %>%
  mutate(sex = factor(sex, levels = c("M","F")),
         population = factor(population, levels = pop.select))

# df_obesity <- read.delim("fig_out/04.biometrics/lmer_bmi_data.tsv") %>% rename(sampleid=sid)

sugar_intake <- c("coffee_tea","added_sugar")

df_sugar <-  #read.delim("input_files/ffq_corrected_imputed.tsv") %>%
  read.delim("input_files/ffq_corrected_imputed_median.tsv") %>%
  select(sampleid, all_of(sugar_intake)) %>%
  mutate(population=substr(sampleid,1,3)) %>%
  filter(!is.na(coffee_tea) & !is.na(added_sugar)) %>%
  left_join(df_obesity[,c("sampleid","age","sex")],by = "sampleid") %>%
  group_by(population, sex) %>%
  mutate(coffee_tea = ifelse(is.na(coffee_tea), 
                             median(coffee_tea, na.rm = TRUE), 
                             coffee_tea),
         added_sugar = ifelse(is.na(added_sugar), 
                              median(added_sugar, na.rm = TRUE), 
                              added_sugar)) %>%
  ungroup() %>%
  mutate(monthly_kg=round(12.5*coffee_tea*added_sugar/1000,3),
         daily_g=round(12.5*coffee_tea*added_sugar/30,0)) %>%
  select(sampleid,daily_g)


### REPORTED STATS
df <- df_sugar %>%
  filter(!is.na(daily_g)) %>%
  mutate(population=substr(sampleid, 1,3)) %>%
  left_join(HGgroup, by="population") %>%
  left_join(df_obesity[,c("sampleid","age","sex")], by="sampleid") %>%
  filter(sampleid %in% df_sugar$sampleid)

n_distinct(df$sampleid)

shapiro.test(df$daily_g)
mean(df$daily_g)
sd(df$daily_g)
median(df$daily_g)
quantile(df$daily_g)
min(df$daily_g)
max(df$daily_g)

df %>% 
  group_by(group) %>%
  summarise(
    n=n_distinct(sampleid),
    daily_intake=mean(daily_g),
    SD=sd(daily_g))

x <- sort(df$daily_g[df$group=="earlyTransition"])
quantile(x,.95)


### Define data

fitdata <- df_obesity %>%
  select(sampleid, population, sex, age, bmi) %>%
  right_join(HGgroup, by="population") %>%
  left_join(df_sugar, by="sampleid") %>%
  filter(!is.na(daily_g)) %>%
  mutate(sugar_Z = as.numeric(scale(as.numeric(daily_g)))) %>%
  group_by(population, sex) %>%
  mutate(across(
    .cols = c(age),
    .fns = ~ ifelse(is.na(.), median(., na.rm = TRUE), .)
  )) %>%
  ungroup() %>%
  mutate(population=factor(population, pop.ord),
         sex=factor(sex, c("F","M")))

n_distinct(fitdata$sampleid)

round(prop.table(table(fitdata$daily_g>50, fitdata$group),margin = 2),2)

fitM <- subset(fitdata, sex =="M")
round(prop.table(table(fitM$daily_g>50, fitM$group),margin = 2),2)

fitF <- subset(fitdata, sex =="F")
round(prop.table(table(fitF$daily_g>50, fitF$group),margin = 2),2)


# plot
p_sugar_pop <- 
  fitdata %>%
  mutate(sex=factor(sex, c("M","F"))) %>%
  ggplot(aes(x=population, y=daily_g, shape=group, colour = group)) +
  geom_point(alpha=0.4, size=2,
             position=position_jitter(width = 0.2, seed=1))+
  geom_boxplot(fill=NA, outlier.shape = NA, linewidth=0.5)+
  scale_shape_manual(values = groupShape) +
  scale_color_manual(values = groupCols) +
  custom_theme

# boxplot by group
fitdata_labelled <- fitdata %>%
  mutate(sex=factor(sex, c("M","F"))) %>%
  group_by(sex, group) %>%
  mutate(n = n(),
         group_sex_label = paste0(group, "\n(", n, ")")) %>%
  ungroup() %>%
  arrange(sex, group) %>% 
  mutate(group_sex_label = factor(group_sex_label, levels = unique(group_sex_label)))

p_sugar <- 
  ggplot(fitdata_labelled, aes(x = group_sex_label, y = daily_g, shape = group, colour = group)) +
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  geom_point(alpha = 0.4, size = 2,
             position = position_jitter(width = 0.2, seed = 1)) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.5) +
  geom_hline(yintercept = 50, linewidth = 1, linetype = "dashed", col = "red") +
  scale_shape_manual(values = groupShape) +
  scale_color_manual(values = groupCols) +
  custom_theme +
  guides(fill = guide_legend(ncol = 2), col = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 9, angle=45, hjust=1)) +
  labs(title = "Daily Sugar Intake", x = NULL, y = "Sugar (Grams/Day)")

p_sugar

wilcox_data <- fitdata %>%
  mutate(sex=factor(sex, c("M","F"))) %>%
  group_by(group, sex) %>%
  filter(n_distinct(sampleid) >= 5) %>%
  ungroup() %>%
  group_by(sex) %>%
  filter(n_distinct(group) >= 2) %>%
  ungroup()

metrics_list <- "daily_g"

wilcox_res <- list("M" = list(), "F" = list())

for(s in c("M", "F")) {
    dat <- wilcox_data %>% filter(sex == s)
    
    resp <- dat[["daily_g"]]
    grp <- dat[["group"]]
    
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
      out_final <- data.frame(sex = s, out_df)
      wilcox_res[[s]] <- out_final}
}

# Combine all list elements into one master data frame
final_results <- bind_rows(lapply(wilcox_res, bind_rows))

wilcox_annot <- final_results %>%
  mutate(sex=factor(sex, c("M","F")),
         group1=factor(group1, levels(fitdata$group)),
         group2=factor(group2, levels(fitdata$group)),
         x1=as.numeric(group1),
         x2=as.numeric(group2),
         xlab=(x2+x1)/2) %>%
  group_by(sex) %>%
  mutate(
    y_base = max(fitdata$daily_g, na.rm = TRUE) - 100 ,
    y_pos = y_base + (row_number() * (y_base * 0.25)),
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )) %>%
  ungroup()

p_sugar <- p_sugar +
  geom_segment(
    data = wilcox_annot,
    aes(x = x1, xend = x2, y = y_pos, yend = y_pos),
    inherit.aes = F, colour = "black"
  ) +
  geom_text(
    data = wilcox_annot, size=5,
    aes(x = xlab, y = y_pos-5, label = label),
    inherit.aes = F, colour = "black", vjust = 0
  )


# save figs
figout="fig_out/05.non_food/sugar_popSex"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = p_sugar_pop,filename = outpath, device = dev, width = 6,height = 4,dpi = 300, scale=1.2)
}

figout="fig_out/05.non_food/sugar_groupSex"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = p_sugar,filename = outpath, device = dev, width = 3,height = 4,dpi = 300, scale=1.2)
}

# stats
max(fitdata$daily_g)
mean(fitdata$daily_g)
fitdata %>% group_by(group) %>% summarise(n=length(!is.na(daily_g)),
                                          max=max(daily_g, na.rm = T),
                                          mean=mean(daily_g, na.rm=T),
                                          se=sd(daily_g,na.rm=T)/sqrt(length(!is.na(daily_g))))
fitdata %>%
  filter(group == "earlyTransition") %>%
  group_by(sex) %>%
  summarise(q95 = quantile(daily_g, probs = 0.95),
            q100 = max(daily_g, na.rm = T))

fitdata %>%
  filter(group == "lateTransition") %>%
  group_by(sex) %>%
  summarise(q95 = quantile(daily_g, probs = 0.95),
            q100 = max(daily_g, na.rm = T))

# LMER
fitdata %>%
  filter(!is.na(bmi)) %>%
  filter(population != "LDY") %>%
  #filter(
  #  daily_g >= quantile(daily_g, 0.001, na.rm = TRUE) &
  #    daily_g <= quantile(daily_g, 0.999, na.rm = TRUE)
  #) %>%
  #filter(daily_g >0) %>%
  mutate(
    high_sugar=ifelse(daily_g>50,1,0),
    sex=factor(sex, c("F","M")),
    group=factor(group, c("earlyTransition","lateTransition"))) -> fit_subset

# categorical
fit0 <- lmer(bmi ~ high_sugar + group + sex + age + (1|population), data=fit_subset)
fit1 <- lmer(bmi ~ high_sugar + group * sex + age + (1|population), data=fit_subset)
fit2 <- lmer(bmi ~ high_sugar * sex + sex + age + (1 | population), data=fit_subset)
model_list <- list(Model0=fit0, Model1=fit1, Model2=fit2)
lapply(model_list, function(x) summary(x))

# continuous
fit0 <- lmer(bmi ~ sugar_Z + group +  sex + age + (1|population), data=fit_subset)
fit1 <- lmer(bmi ~ sugar_Z + group * sex + age + (1|population), data=fit_subset)
fit2 <- lmer(bmi ~ sugar_Z * sex + sugar_Z * group + age + (1 | population), data=fit_subset)
fit3 <- lmer(bmi ~ sugar_Z * group + group * sex + age + (1 | population), data=fit_subset)
fit4 <- lmer(bmi ~ sugar_Z * sex + group * sex + age + (1 | population), data=fit_subset)
model_list <- list(Model0=fit0, Model1=fit1, Model2=fit2, Model3=fit3, Model4=fit4)
lapply(model_list, function(x) summary(x))

export_mixed_models(
  model_list,
  folder_path = "out/", 
  base_name = "LMM_Results_Sugar_vs_BMI_raw")


bmi_model <- fit1
summary(bmi_model)
library(performance)
r2_values <- model_performance(bmi_model)
marginal_r2 <- round(r2_values$R2_marginal, 2)

pred_full <- ggpredict(bmi_model, 
                       terms = c("sugar_Z [all]", "group", "sex"),
                       condition = c(age = mean(fit_subset$age, na.rm = TRUE)))

pred_full <- as.data.frame(pred_full)

pred_full <- rename(pred_full, 
                    bmi = predicted,
                    sugar_Z = x,
                    sex = facet)

pred_full$sex <- factor(pred_full$sex, levels = levels(fitdata$sex))
pred_full$group <- factor(pred_full$group, levels = levels(fitdata$group))

pred_full <- pred_full %>%
  group_by(sex) %>%
  filter(sugar_Z < quantile(sugar_Z, 0.999, na.rm = TRUE)) %>%
  ungroup()

icaption <- summary(bmi_model)
icaption <- gsub("population", "community", as.character(icaption$call$formula))[3]
icaption <- paste0("-- LMM: BMI ~ ", icaption)

p_sugar_bmi <- fit_subset %>%
  filter(daily_g < 350) %>%
  mutate(sex=factor(sex, c("M","F"))) %>%
  ggplot(aes(x=sugar_Z, y=bmi, colour = group)) +
  geom_point(aes(shape=group), alpha=0.4, size=2)+
  scale_shape_manual(values = groupShape) +
  scale_color_manual(values = groupCols) +
  #geom_smooth(method="lm") +
  geom_line(data=pred_full, inherit.aes = F, linewidth=1,
            aes(x=sugar_Z, y=bmi, colour = group))+
  facet_grid(~sex, labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  custom_theme +
  theme(axis.text = element_text(size=10))+
  #scale_x_continuous(limits = c(0,300), breaks = seq(0,400, by=100))+
  labs(title="BMI Association with Sugar Intake", 
       subtitle = paste0("Marginal R-squared = ", marginal_r2),
       x="Sugar (Z-score)", 
       y="BMI (kg/m2)",       
       caption = icaption) +
  theme(plot.caption = element_text(face = "plain"),
        plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40"))
  
p_sugar_bmi


figout="fig_out/05.non_food/sugar_lmBMI"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = p_sugar_bmi,filename = outpath, device = dev, width = 6,height = 4,dpi = 300, scale=1.2)
}

# to combine
panel_A <- p_sugar
panel_B <- p_sugar_bmi

saveRDS(panel_A, "fig_out/05.non_food/Panel_A.rds")
saveRDS(panel_B, "fig_out/05.non_food/Panel_B.rds")




# save sample list
length(df$sampleid)

df %>% 
  group_by(group) %>%
  summarise(n=n_distinct(sampleid))

df %>% 
  group_by(sex) %>%
  summarise(n=n_distinct(sampleid))

df %>% 
  mutate(population=factor(population, pop.ord)) %>%
  group_by(population) %>%
  summarise(n=n_distinct(sampleid))


saveRDS(unique(df$sampleid),file = "input_files/sugar_sample_list.rds")
write.table(df, "out/sugar_data_processed.tsv", quote = F, sep = "\t", row.names = F)
