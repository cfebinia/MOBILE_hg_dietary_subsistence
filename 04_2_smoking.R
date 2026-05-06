rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggeffects)


##########
# Helper Functions
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
# scale_shape_manual(values = c("Smoker"=17, "Non-smoker"=16)) +
# scale_color_manual(values = groupCols) +

##########
# load data
smoking <- read.delim("input_files/Smoking_compiled_final - COMPILED_FINAL.tsv") %>% 
  select(sid, Yes_No, pcs_per_day)
age_sex <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>%
  rename(sid=sampleid) %>%
  mutate(population=substr(sid,1,3))%>%
  select(sid, population, sex, age)
smoking <- age_sex %>% right_join(smoking, by="sid", unmatched = "drop")

trg_HG <- c("BTU","ORT","ORS")
trd_HG <- c("APT","TBU","BSP")
pop.ord <- c(trg_HG, trd_HG)

smoking$population <- factor(smoking$population, pop.ord)
table(smoking$population, smoking$sex)


# Reported stats
df <- smoking %>%
  mutate(group=ifelse(population %in% trg_HG, "earlyTransition","lateTransition")) %>%
  mutate(group=factor(group, c("earlyTransition","lateTransition")),
         sex=factor(sex, c("F","M")),
         Smoker=factor(ifelse(Yes_No==0,"Non-smoker","Smoker"))) %>%
  rename(intensity=pcs_per_day) %>%
  filter(!is.na(Yes_No)) %>%
  select(-Yes_No)

nrow(df)

wilcox.test(intensity ~ group, 
            data = df[!is.na(df$intensity),], 
            subset = sex == "M")

wilcox.test(intensity ~ group, 
            data = df[!is.na(df$intensity),], 
            subset = sex == "F")

df %>%
  filter(!is.na(intensity)) %>%
  filter(Smoker=="Smoker") %>%
  group_by(Smoker, sex, group) %>%
  summarise(n = n(), 
            median = round(median(intensity),0),
            SD = round(sd(intensity),0),
            q25 = quantile(intensity, 0.25),
            q75 = quantile(intensity, 0.75))

df %>%
  group_by(group, Smoker) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_by_group = sum(n),
         proportion = round(n / total_by_group)) %>%
  ungroup()

groupCols <- c("lateTransition" = "maroon", "earlyTransition" = "forestgreen")
groupShape <- c("lateTransition" = 15, "earlyTransition" = 16) # this one need flat colours
groupshape2 <- c("lateTransition" = 22, "earlyTransition" = 21) # this one can have borders


# PREVALENCE
# ggdata <- data.frame(table(smoking$population, smoking$sex, smoking$Yes_No))
# colnames(ggdata) <- c("population","sex","smoking","count")

ggdata <- smoking
gr <- c(rep("earlyTransition",length(trg_HG)),
        rep("lateTransition", length(trd_HG)))
names(gr) <- c(trg_HG, trd_HG)

ggdata$group <- gr[ggdata$population]
ggdata$group <- factor(ggdata$group, c("earlyTransition","lateTransition"))

ggdata %>% filter(!is.na(Yes_No)) %>%
  group_by(sex, population, group) %>% 
  summarise(smoking=sum(Yes_No==1),
            total=length(Yes_No),
            preval=round(100*smoking/total,2)) %>%
  mutate(sex=factor(sex, c("M","F"))) -> ggdata

mytitle <- paste0("Smoking, n=", nrow(smoking))

p <- ggplot(ggdata) + custom_theme +
  geom_bar(aes(x=population, y=preval, fill=group), stat = "identity") +
  scale_fill_manual(values = c(earlyTransition = "forestgreen", lateTransition = "maroon"))+
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  guides(fill = guide_legend(ncol = 1), col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))+
  scale_y_continuous(limits=c(0,120), breaks = seq(0,100, by=50), labels = paste0(seq(0,100, by=50), "%"))+
  ggtitle(mytitle) + ylab("prevalence")

group_stats <- ggdata %>%
  group_by(sex, group) %>%
  summarise(
    total_smoking = sum(smoking),
    total_n = sum(total),
    .groups = "drop"
  ) %>%
  group_by(sex) %>%
  # Perform the proportion test for each sex
  summarise(
    p_val = prop.test(x = total_smoking, n = total_n)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    # Create the significance stars
    label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    # Position the star in the middle of the plot area
    y_pos = 110 
  )

p <- p + 
  geom_text(data = group_stats, size = 5,
            aes(x = 3.5, y = 110, label = label), 
            inherit.aes = FALSE, 
            fontface = "bold", color = "black") +
  # Expand limits so the stars don't get cut off
  coord_cartesian(ylim = c(0, 120))

p

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","05.non_food","smoking_prevalence"),".",dev)
  ggsave(plot = p,filename = outpath, device = dev, width = 4,height = 4, scale=1.3)
}


# INTENSITY
# Age bins
age_bins <- smoking[,c("sid","age")] %>%
  filter(!is.na(age)) %>%
  mutate(age_bin = case_when(
    age < 20 ~ 20,
    age >= 60 ~ 60,
    TRUE ~(age %/% 20) * 20)
  ) %>%
  mutate(age_lab=case_when(
    age_bin==20 ~ "18-39",
    #age_bin==30 ~ "30-39",
    age_bin==40 ~ "40-59",
    #age_bin==50 ~ "50-59",
    age_bin==60 ~ "≥60"
  )) %>%
  mutate(age_lab=factor(age_lab, levels = c("18-39","40-59","≥60")))

table(age_bins$age_lab)
table(substr(age_bins$sid,1,3),age_bins$age_lab)  

age_bin_labs <- distinct(age_bins[,c("age_bin","age_lab")])
tmp <- as.character(age_bin_labs$age_lab)
names(tmp) <- as.numeric(as.character(age_bin_labs$age_bin))
age_bin_labs <- tmp[order(names(tmp),decreasing = F)]
age_bin_labs

smoking$group <- gr[smoking$population]

smoking <- smoking %>%
  select(-age) %>%
  left_join(age_bins, by="sid")

table(smoking$group, smoking$age_bin, smoking$sex)

smoking %>% 
  filter(Yes_No==1) %>%
  filter(!is.na(pcs_per_day)) %>%
  filter(!is.na(age)) -> wilcox_data

library(ggplot2)
library(ggpubr)

wilcox_data$group <- factor(wilcox_data$group, c("earlyTransition","lateTransition"))
wilcox_data$sex <- factor(wilcox_data$sex, c("M","F"))

p <-
  ggplot(wilcox_data) +
  custom_theme +
  facet_grid(~ sex, scales = "free_x", space = "free_x") +
  geom_point(alpha=0.6,size=2,
             aes(x=age_bin, y=pcs_per_day, shape=group, fill=group, colour = group, group = interaction(group, age_bin)),
             position = position_jitterdodge(seed=1, jitter.width = 5, dodge.width = 12,jitter.height = 0))+
  geom_boxplot(aes(x=age_bin, y = pcs_per_day, colour = group, group = interaction(group, age_bin)), 
               fill = NA, width = 12, outlier.shape = NA)+
  scale_x_continuous(labels = age_bin_labs, breaks = seq(20,60, by=20), expand = c(0.1,0)) +
  scale_shape_manual(values = groupShape) +
  scale_color_manual(values = groupCols) +
  labs(
    title = "Adult Smoking Habit",
    subtitle = paste0("(n=", n_distinct(wilcox_data$sid),")"),
    caption = "p-value from Wilcoxon Test",
    x = "Age (Years)",
    y = "Intensity (Pcs/Day)",
    colour = "HG Group"
  ) 

signif_wilcox <- wilcox_data %>%
  group_by(age_bin, sex, group) %>%
  filter(n_distinct(sid) >= 3) %>%
  group_by(age_bin, sex) %>%
  filter(n_distinct(group) == 2) %>%
  summarise(
    p_val = wilcox.test(pcs_per_day ~ group, exact=F)$p.value,
    max_y = max(pcs_per_day, na.rm = TRUE) +3,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_val, method = "fdr"), 
    labs = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    x_center = age_bin, 
    x1 = x_center - 3, 
    x2 = x_center + 3
  )

p <- p + 
  guides(fill = guide_legend(nrow = 2, byrow = T)) +
  guides(col = guide_legend(nrow = 2, byrow = T)) +
  geom_segment(data = signif_wilcox, 
               aes(x = x1, xend = x2, y = max_y + 2, yend = max_y + 2), 
               inherit.aes = FALSE) +
  geom_segment(data = signif_wilcox, 
               aes(x = x1, xend = x1, y = max_y + 2, yend = max_y + 1), 
               inherit.aes = FALSE) +
  geom_segment(data = signif_wilcox, 
               aes(x = x2, xend = x2, y = max_y + 2, yend = max_y + 1), 
               inherit.aes = FALSE) +
  geom_text(data = signif_wilcox, size = 5,
            aes(x = x_center, y = max_y+0.7, label = labs), 
            inherit.aes = FALSE)

p

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","05.non_food","smoking_intensity"),".",dev)
  ggsave(plot = p,filename = outpath, device = dev, width = 6,height = 6, scale=0.8)
}


# box
smoking %>% 
  filter(!is.na(pcs_per_day)) %>%
  #filter(!is.na(age)) %>%
  mutate(group=ifelse(population %in% c("BTU","ORS","ORT"), "earlyTransition","lateTransition"))%>%
  mutate(sex=factor(sex, c("M","F")),
         group=factor(group, levels = c( "earlyTransition","lateTransition")))-> ggdata

ggdata0<- ggdata
ggdata <- ggdata0 %>% mutate(Smoking=ifelse(Yes_No==1, "Smoker","Non-smoker"))
# %>% filter(Yes_No==1)

ggdata_labelled <- ggdata %>%
  filter(Smoking == "Smoker") %>%
  group_by(sex, group) %>%
  mutate(n = n(),
         group_sex_label = paste0(group, "\n(", n, ")")) %>%
  ungroup() %>%
  arrange(sex, group) %>% 
  mutate(group_sex_label = factor(group_sex_label, unique(group_sex_label)))

label_lookup <- ggdata_labelled %>%
  distinct(sex, group, group_sex_label)

p_box <- ggplot(ggdata_labelled,
                aes(x = group_sex_label, y = pcs_per_day, colour = group)) +
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  geom_point(aes(shape = Smoking),
             alpha = 0.4, size = 2, 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.7, seed = 1)) +
  geom_boxplot(aes(group = interaction(group_sex_label, Smoking)),
               fill = NA, outlier.shape = NA, linewidth = 0.2, width=0.8,
               position = position_dodge(width = 0.7)) +
  scale_shape_manual(values = c("Smoker" = 17, "Non-smoker" = 16)) +
  scale_color_manual(values = groupCols) +
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Smoking Intensity", x = "", y = "Pcs/Day")

# Statistical Analysis
wilcox_data <- ggdata_labelled %>%
  group_by(group, sex) %>%
  filter(n_distinct(sid) >= 5) %>%
  ungroup() %>%
  group_by(sex) %>%
  filter(n_distinct(group) >= 2) %>%
  ungroup()

wilcox_res <- list()

for(s in c("M", "F")) {
  dat <- wilcox_data %>% filter(sex == s)
  if(nrow(dat) == 0) next
  
  out_matrix <- pairwise.wilcox.test(
    x = dat$pcs_per_day, g = dat$group, 
    p.adjust.method = "fdr", 
    paired = FALSE, exact = FALSE
  )$p.value
  
  if(!is.null(out_matrix)){
    out_df <- out_matrix %>% 
      as.table() %>% 
      as.data.frame() %>% 
      filter(!is.na(Freq)) %>% 
      rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
      filter(p_adj < 0.05) %>%
      mutate(p_adj = round(p_adj, 4), sex = s)
    
    wilcox_res[[s]] <- out_df
  }
}

final_results <- bind_rows(wilcox_res)

# Annotation Alignment
wilcox_annot <- final_results %>%
  inner_join(label_lookup, by = c("sex", "group1" = "group")) %>%
  rename(label1 = group_sex_label) %>%
  inner_join(label_lookup, by = c("sex", "group2" = "group")) %>%
  rename(label2 = group_sex_label) %>%
  mutate(sex = factor(sex, levels = c("M", "F")),
         group1 = factor(group1, unique(gr)),
         group2 = factor(group2, unique(gr))) %>%
  group_by(sex) %>%
  mutate(
    x1 = as.numeric(group1), 
    x2 = as.numeric(group2),
    y_base = max(ggdata$pcs_per_day, na.rm = TRUE),
    y_pos = y_base + (row_number() * (y_base * 0.1)),
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )) %>%
  ungroup()

# Add Annotations to Plot
p_box <- p_box +
  geom_segment(
    data = wilcox_annot,
    aes(x = x1, xend = x2, 
        y = y_pos, 
        yend = y_pos),
    inherit.aes = FALSE, colour = "black"
  ) +
  geom_text(
    data = wilcox_annot, size = 5,
    aes(x = (x1 + x2)/2, 
        y = y_pos - 5, 
        label = label),
    inherit.aes = FALSE, colour = "black", vjust = -0.5,
    stat = "unique") %>% 
  { # Re-calculate x-midpoint manually to avoid factor issues
    .$data$xlab <- (as.numeric(.$data$label1) + as.numeric(.$data$label2)) / 2
    .
  }
p_box

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","05.non_food","smoking_intensity_box"),".",dev)
  ggsave(plot = p_box,filename = outpath, device = dev, width = 3,height = 4, scale=1.2)
}

ggdata %>% group_by(sex, group) %>% summarise(n=length(!is.na(pcs_per_day)))


n_distinct(ggdata0$sid)
table(ggdata0$population)
table(ggdata0$Yes_No)
table(ggdata0$sex)
ggdata0 %>% group_by(sex) %>% filter(Yes_No==1) %>% summarise(n=n_distinct(sid))
table(ggdata0$group, ggdata0$sex)
ggdata0 %>% group_by(sex, group) %>% filter(Yes_No==1) %>% summarise(n=n_distinct(sid))

ggdata0 %>% group_by(sex, population) %>% filter(Yes_No==1) %>% summarise(n=n_distinct(sid)) %>% filter(population=="BTU")
ggdata0 %>% group_by(sex, population) %>% summarise(n=n_distinct(sid)) %>% filter(population=="BTU")

ggdata %>% group_by(sex, group) %>% filter(Yes_No==1) %>% summarise(median=median(pcs_per_day, na.rm=T),
                                                                    q25=quantile(pcs_per_day, 0.25),
                                                                    q75=quantile(pcs_per_day, 0.75))


### vs BMI
#df_obesity <- read.delim("input_files/2023_MOBILE_AnthropometryData(RAW).tsv") %>%
df_obesity <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>%
  mutate(population = substr(sampleid, 1,3)) %>%
  mutate(bmi=round(w_kg/((h_cm/100)^2),1)) %>%
  filter(population %in% pop.ord) %>%
  select(sampleid, bmi) %>%
  rename(sid=sampleid)

boxdata <- smoking  %>%
  select(sid, population, group, sex, age, Yes_No)  %>%
  left_join(df_obesity, by="sid") %>%
  filter(!is.na(bmi)) %>%
  mutate(Smoking=ifelse(Yes_No==1, "Smoker", "Non-smoker")) %>%
  mutate(sex=factor(sex, c("M","F")))

boxdata_labelled <- boxdata %>%
  group_by(sex, group) %>%
  mutate(n = n(),
         group_sex_label = paste0(group, "\n(", n, ")")) %>%
  ungroup() %>%
  arrange(sex, group) %>% 
  mutate(group_sex_label = factor(group_sex_label, unique(group_sex_label)))

label_lookup <- boxdata_labelled %>%
  distinct(sex, group, group_sex_label)

p_smoking_bmi <- boxdata_labelled %>%
  ggplot(aes(x = group_sex_label, y = bmi, colour = group)) +
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  geom_point(aes(shape = Smoking, group = interaction(Smoking, group)),
             alpha = 0.6, size = 2,
             position = position_jitterdodge(seed = 1, jitter.width = 0.2, dodge.width = 0.75)) +
  geom_boxplot(aes(group = interaction(Smoking, group)), 
               fill = NA, width = 0.7, outlier.shape = NA,
               position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c("Smoker"=17, "Non-smoker"=16)) +
  scale_color_manual(values = groupCols) +
  scale_y_continuous(breaks = seq(5, 50, by = 5)) +
  custom_theme +
  labs(
    title = "BMI of Smokers and Non-smokers",
    # subtitle = paste0("(n=", n_distinct(boxdata_labelled$sid), ")"),
    caption = "p-value from Wilcoxon Test",
    x = "Transition State",
    y = "BMI (kg/m2)",
    colour = "HG Group"
  )

p_smoking_bmi

stats_data <- boxdata_labelled %>%
  group_by(sex, group, Smoking) %>%
  filter(n_distinct(sid) > 2) %>%
  ungroup() %>%
  group_by(sex, group) %>%
  summarise(
    n = n_distinct(sid),
    p_val = wilcox.test(bmi ~ Smoking, exact=F)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_val, method = "fdr"),
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    y_start = max(boxdata$bmi) - 2,
    y_text = y_start + 0.2,
    x0 = as.numeric(factor(group))
  )

d_width <- 0.8

plot_out0 <- p_smoking_bmi + 
  scale_y_continuous(limits = c(min(boxdata_labelled$bmi),2+max(stats_data$y_text)))+
  guides(shape = guide_legend(nrow = 2, byrow = T)) +
  guides(col = guide_legend(nrow = 2, byrow = T)) +
  geom_segment(data = stats_data, inherit.aes = F,
               aes(x = x0 - d_width/4, 
                   xend = x0 + d_width/4, 
                   y = y_start, yend = y_start),
               colour = "black") +
  geom_text(data = stats_data, inherit.aes = F, size = 5, 
            aes(x = x0, y = y_text , label = label),
            vjust = 0)

plot_out0
plot_out <- p_smoking_bmi

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","05.non_food","smoking_vs_bmi"),".",dev)
  ggsave(plot = plot_out0,filename = outpath, device = dev, width = 5, height = 4, scale=1.2)
}

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","05.non_food","smoking_vs_bmi_noAnnot"),".",dev)
  ggsave(plot = p_smoking_bmi,filename = outpath, device = dev, width = 5, height = 4, scale=1.2)
}
stats_data

##############
## LMM
library(lmerTest)
table(boxdata$sex, boxdata$group)

fitdata <- merge.data.frame(smoking, df_obesity[,c("sid","bmi")], by="sid") %>%
  mutate(
    Smoker=Yes_No,
    sex=factor(sex, c("F","M")),
    group=factor(group, unique(gr)))
dim(fitdata)
table(fitdata$sex, fitdata$group)

# by smoker vs non-smoker
fit0 <- lmer(bmi ~ Smoker + group + sex + age + (1|population), data=fitdata)
summary(fit0)

fit1 <- lmer(bmi ~ Smoker + group + sex +  age + Smoker:group + (1|population), data=fitdata)
summary(fit1)

fit2 <- lmer(bmi ~ Smoker  + group + sex + age + Smoker:sex + group:sex + (1|population), data=fitdata)
summary(fit2)

model_list <- list(Model0 = fit0, Model1 = fit1, Model2=fit2)
lapply(model_list, function (x) summary(x))

lmm_path="out/"
export_mixed_models(
  model_list,
  folder_path = lmm_path,
  base_name = "LMM_bmi_vs_smoking")

# by intensity
fit0 <- lmer(bmi ~ pcs_per_day + sex + age + (1|population), data=fitdata)
summary(fit0)

fit1 <- lmer(bmi ~ pcs_per_day + sex + age + pcs_per_day:sex + (1|population), data=fitdata)
summary(fit1)

fit2 <- lmer(bmi ~ pcs_per_day + group + sex + group + age + pcs_per_day:sex + group:sex + (1|population), data=fitdata)
summary(fit2)

model_list <- list(Model0 = fit0, Model1 = fit1, Model2=fit2)

lmm_path="out"
export_mixed_models(
  model_list,
  folder_path = lmm_path,
  base_name = "LMM_bmi_vs_smoking_intensity")



# plot
# Ensure Smoker is a factor
fitdata$Smoker <- factor(fitdata$Smoker, levels = c(0, 1))

# Refit model and generate predictions as before
bmi_model <- lmer(bmi ~ Smoker * sex + group * sex + age + (1|population), data = fitdata)
pred_full <- ggpredict(bmi_model, terms = c("Smoker", "group", "sex"),
                       condition = c(age = mean(fitdata$age, na.rm = TRUE)))
pred_full <- as.data.frame(pred_full)
pred_full <- rename(pred_full, bmi = predicted, Smoker = x, sex = facet)

library(performance)
r2_values <- model_performance(bmi_model)
marginal_r2 <- round(r2_values$R2_marginal, 2)

# Plot with updated x-axis labels
p_smoking_bmi2 <- fitdata %>%
  filter(!is.na(bmi), !is.na(Smoker), !is.na(sex)) %>%
  mutate(sex = factor(sex, levels = c("M", "F"))) %>%
  ggplot(aes(x = Smoker, y = bmi, colour = group)) +
  # Points: set dodge.width to 0.8 for clarity
  geom_point(aes(shape = Smoker), alpha = 0.4, size = 2,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8, seed = 1)) +
  # Boxplots: dodge.width matches the points
  geom_boxplot(aes(group = interaction(Smoker, group)), 
               linewidth = 0.5, width = 0.5, 
               outlier.shape = NA, fill = NA,
               position = position_dodge(width = 0.8)) +
  # Predicted lines: dodge.width matches the points and boxplots
  geom_line(data = pred_full, linewidth = 1,
            aes(x = Smoker, y = bmi, colour = group, group = group),
            position = position_dodge(width = 0.8)) +
  scale_x_discrete(labels = c("0" = "Non-smoker", "1" = "Smoker")) +
  scale_shape_manual(values = c("0" = 16, "1" = 17), 
                     labels = c("0" = "Non-smoker", "1" = "Smoker")) +
  scale_color_manual(values = groupCols) +
  facet_grid(~sex, labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  custom_theme + 
  labs(title = "BMI Association with Smoking Status", 
       subtitle = paste0("Marginal R-squared = ", marginal_r2),
       x = "Smoking Status", 
       y = "BMI (kg/m2)",
       caption = "—— LMM: BMI ~ Smoker + sex + state + age + Smoker:sexM + state:sexM + (1|community)") +
  theme(plot.caption = element_text(face = "plain"),
        plot.subtitle = element_text(face = "italic", size = 9, colour = "grey40"))

p_smoking_bmi2

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","05.non_food","LMM_smoking_vs_bmi"),".",dev)
  ggsave(plot = p_smoking_bmi2,filename = outpath, device = dev, width = 5, height = 4, scale=1.2)
}


# to combine
panel_C <- p_box
panel_D <- plot_out0

saveRDS(panel_C, "fig_out/05.non_food/Panel_C.rds")
saveRDS(panel_D, "fig_out/05.non_food/Panel_D.rds")


# save sample list
length(df$sid)

df %>% 
  group_by(group) %>%
  summarise(n=n_distinct(sid))

df %>% 
  group_by(sex) %>%
  summarise(n=n_distinct(sid))

df %>% 
  mutate(population=factor(population, pop.ord)) %>%
  group_by(population) %>%
  summarise(n=n_distinct(sid))


saveRDS(unique(df$sid),file = "input_files/smoking_sample_list.rds")
write.table(df, "out/smoking_data_processed.tsv", quote = F, sep = "\t", row.names = F)

