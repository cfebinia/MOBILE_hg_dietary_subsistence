rm(list=ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

population_colors <- readRDS("input_files/mycols.rds")
popcols <- population_colors$col
names(popcols) <- population_colors$population


#######################
# Food Diversity Score
#######################
# this is a modified scoring. As long one reported consuming one item in the food groups, it's considered as True.
# further guidelines from FAO: https://openknowledge.fao.org/handle/20.500.14283/cb3434en (https://doi.org/10.4060/cb3434en)
df.con <- read.delim("input_files/consumption_corrected.tsv") %>% select(-population)

bio <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>% 
  select(sampleid, population, sex) %>%
  mutate(population=substr(sampleid,1,3))

df.con <- bio %>%
  left_join(df.con, by="sampleid",unmatched = "drop") %>%
  filter(if_any(all_of(colnames(df.con)[-1]), ~ .x != 0))# remove samples if variables is 0 across
  
food.groups <- read.csv("input_files/food_groups.csv") %>%
  filter(!is.na(food_groups1))

pop.ord <- c("BTU", "ORT", "ORS", "APT", "TBU", "BSP", "LDY")
pop.select <- pop.ord
mdd.cutoff <- round(length(unique(food.groups$food_groups1)) * 0.75, 0)

idds <- df.con %>%
  as.data.frame() %>%
  # Filter populations and replace NA
  mutate(population = substr(sampleid, 1, 3)) %>%
  filter(population %in% pop.ord) %>%
  mutate(across(all_of(food.groups$food), ~ replace_na(.x, 0))) %>%
  # Reshape to long for group-wise consumption check
  pivot_longer(
    cols = all_of(food.groups$food),
    names_to = "food_item",
    values_to = "value"
  ) %>%
  left_join(food.groups, by = c("food_item" = "food")) %>%
  group_by(sampleid, population, food_groups1) %>%
  summarise(any_consumed = any(value > 0), .groups = "drop") %>%
  # Reshape to wide matrix format
  pivot_wider(names_from = food_groups1, values_from = any_consumed, values_fill = FALSE) %>%
  # Calculate Scores and Categories
  mutate(
    idds = rowSums(across(where(is.logical))),
    group = case_when(
      population %in% c("BTU", "ORT", "ORS") ~ "earlyTransition",
      population %in% c("APT", "TBU", "BSP") ~ "lateTransition",
      population == "LDY" ~ "Agri",
      TRUE ~ NA_character_
    ),
    mdd = if_else(idds >= mdd.cutoff, "High DDS", "Low DDS")
  ) %>%
  mutate(population = factor(population, levels = pop.select),
         group=factor(group, levels = c("earlyTransition","lateTransition","Agri"))) %>%
  filter(idds > 0)

food_category_list <- idds %>% 
  select(where(is.logical)) %>% 
  colnames()

idds_stats <- idds %>%
  group_by(population) %>%
  summarise(
    n = n(),
    med = median(idds),
    ave = mean(idds),
    se = sd(idds) / sqrt(n()),
    .groups = "drop"
  )

outliers <- idds %>%
  group_by(population) %>%
  group_modify(~ {
    out <- boxplot.stats(.x$idds)$out
    .x %>% mutate(is_outlier = idds %in% out)
  }) %>%
  ungroup() %>%
  mutate(outlier_label = if_else(is_outlier, sampleid, ""))

# xlab
idds %>%
  group_by(population) %>%
  summarise(n=n_distinct(sampleid)) %>%
  mutate(xlab=paste0(population," (", n, ")")) -> n_pop

xlab <- n_pop$xlab
names(xlab) <- n_pop$population

myshape <- c(earlyTransition = 16,
             lateTransition = 15,
             Agri = 17)

groupcol <- c(earlyTransition = "forestgreen",
              lateTransition = "maroon",
              Agri = "grey50")

p <- ggplot(idds) + theme_bw(base_size = 12) + 
  theme(axis.title.x = element_blank(), 
        panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45,hjust = 1)) +
  geom_vline(xintercept = c(3.5,6.5), linetype=2, linewidth=1.2, col="grey60")+
  geom_point(aes(x=population, y=idds, col=group, fill=group, shape=group),
             size=2, alpha=0.5, position = position_jitter(height = 0.1,width = 0.2,seed = 100)) +
  geom_boxplot(aes(x=population, y=idds, col=group),
               fill=NA, outlier.shape = NA, linewidth=1, alpha=0.1, width=0.8) + 
  scale_color_manual(values = groupcol)+
  scale_shape_manual(values = myshape)+
  scale_x_discrete(labels=xlab)+
  labs(title="Individual Food Diversity Score (IDDS)")
p

pairwise.wilcox.test(idds$idds,g = idds$population,p.adjust.method = "fdr",paired = F, exact=F) -> idds_pwilcox
pairwise.wilcox.test(idds$idds,g = idds$group,p.adjust.method = "fdr",paired = F, exact=F) -> idds_pwilcox_gr
print(idds_pwilcox)
print(idds_pwilcox_gr)

to_lab <- data.frame(x1=c(1,4,2,2,5), x2=c(3,6,5,2,5), y1=c(11.5,11.5,12,11.5,11.5), y2=c(11.5,11.5,12,12,12))
plab <- paste0("q = ",round(idds_pwilcox_gr$p.value["lateTransition","earlyTransition"],3),"")

p <- p +   
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), breaks = seq(0,12,by=2))+
  geom_segment(data = to_lab, aes(x = x1, xend = x2, y = y1, yend = y2)) +
  annotate("text", x = 3.5, y = 12.7, label = plab, vjust = 1.5)

print(p)

plot_idds <- p

figout=file.path(getwd(),"fig_out","03.dietary_assesment","IDDS_boxplot")
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = plot_idds,filename = outpath, device = dev, width = 2.5,height = 3.5,dpi = 600, scale=1.2)
}

idds %>%
  summarise(n=n_distinct(sampleid))

sink("fig_out/03.dietary_assesment/food_diversity_index_pw_wilcox.txt")
print(idds_pwilcox)
print(idds_pwilcox_gr)
sink()


################
# Food Frequency
###############
# load dataset
dff <- read.delim("input_files/ffq_corrected.tsv")
dff_imp <- read.delim("input_files/ffq_corrected_imputed_median.tsv")
dff_Z <- read.delim("input_files/ffq_corrected_imputed_median_Zscore.tsv")

bio <- read.delim("input_files/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>% 
  select(sampleid, population, sex) %>%
  mutate(population=substr(sampleid,1,3))

ffq_list <- list(dff=dff, dff_imp=dff_imp, dff_Z=dff_Z)

for(dat in names(ffq_list)){
  tmp <- ffq_list[[dat]]
  tmp <- bio %>%
    left_join(tmp, by="sampleid",unmatched = "drop") %>%
    filter(if_any(all_of(colnames(tmp)[-1]), ~ .x != 0))
  ffq_list[[dat]] <- tmp
}

dff <- ffq_list[["dff"]]
dff_imp <- ffq_list[["dff_imp"]]
dff_Z <- ffq_list[["dff_Z"]]

# deal with honey
#dff_imp %>%
#  group_by(population) %>%
#  summarise(mean_honey=mean(honey, na.rm=T)) %>%
#  arrange(desc(mean_honey))
# dff_imp$honey <- ifelse(dff_imp$honey>0,1,0)

# groups:
groupA <- c("fish","white_meat","red_meat")
groupB <- c("rice","legumes")
groupC <- c("tuber","vegetables","fruits")
groupD <- c("eggs","noodles","cooking_oil")
no_var <- c("coffee_tea", "added_sugar", "lard", "sago")
select_var <- c(groupA,"honey",groupB,groupC,groupD,"sugar_g")

food_lvl <- c("Animal-based food","Plant-based food","Market food")

zscore_group <- c(rep("Animal-based food",4), rep("Plant-based food",5),
                  rep("Market food",4))
names(zscore_group) <- select_var  

# Z-Score plot
ggdata <- dff_Z %>%
  select(sampleid, population, all_of(select_var)) %>%
  pivot_longer(
    cols = all_of(select_var), 
    names_to = "variable", 
    values_to = "Zscore"
  ) %>%
  # group_by(variable) %>% mutate(Zscore = as.numeric(scale(freq))) %>%
  ungroup() %>%
  mutate(
    type = factor(zscore_group[variable], levels = unique(zscore_group)),
    group = case_when(
      population %in% c("BTU", "ORT", "ORS") ~ "earlyTransition",
      population %in% c("APT", "TBU", "BSP") ~ "lateTransition",
      population == "LDY" ~ "Agri",
      TRUE ~ NA_character_
    )) %>%
  mutate(
    group = factor(group, levels = c("earlyTransition","lateTransition","Agri")),
    population = factor(population, levels = pop.ord),
    variable = factor(variable, levels = select_var),
    type = factor(type, levels = food_lvl))

# 4. Statistical Summaries
# Consolidated calculation of CI and SE
resort_food <- names(sort(factor(zscore_group,levels = food_lvl)))

ggdata %>%
  group_by(population) %>%
  summarise(n=n_distinct(sampleid)) %>%
  mutate(xlab=paste0(population,"\n(", n, ")")) -> n_pop

xlab <- n_pop$xlab
names(xlab) <- n_pop$population

ggdata_mean <- ggdata %>%
  group_by(variable, type, population) %>%
  summarise(
    n = n(),
    mean = mean(Zscore, na.rm = TRUE),
    se = sd(Zscore, na.rm = TRUE) / sqrt(n),
    tcrit = qt(0.975, n - 1), # Corrected for two-tailed 95% CI
    #lowerCI = mean - (tcrit * se),
    #upperCI = mean + (tcrit * se),
    lowerCI = mean - se,
    upperCI = mean + se,
    .groups = "drop"
  ) %>%
  mutate(
    population = factor(population, levels = pop.ord),
    variable = factor(variable, levels = resort_food),
    #type = factor(type, levels = c("Staples","Animals","Plants","Other")
    type = factor(type, levels = food_lvl
    ))

ylim <- max(abs(c(ggdata_mean$upperCI, ggdata_mean$lowerCI)))

# 5. Visualisation
library(RColorBrewer)

set.seed(1)
groupA_cols <- brewer.pal(7, "Set1") %>% sample(length(groupA)+1)
groupB_cols <- brewer.pal(5, "YlGn") %>% sample(length(groupB))
groupC_cols <- brewer.pal(6, "Dark2") %>% sample(length(groupC)+1)
groupD_cols <- brewer.pal(6, "Set2") %>% sample(length(groupD))

foodcols <- c(groupA_cols,groupB_cols,groupC_cols,groupD_cols)
names(foodcols) <- names(zscore_group)
foodcols["sugar_g"] <- "#528B8B"
foodcols["noodles"] <- "#00CDCD"
foodcols["honey"] <- "#333"
foodcols["cooking_oil"] <- "#FFB90F"

diet.freqlines <- 
  ggplot(ggdata_mean, 
         aes(x = population, y = mean, col = variable, group = variable)) +
  geom_vline(xintercept = c(3.5,6.5), linetype=2, linewidth=0.8, col="grey60")+
  theme_bw(base_size = 12) +
  facet_wrap(~type, nrow = 1) +
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2, linewidth=0.8) +
  geom_line(data=subset(ggdata_mean, population != "LDY"),
            linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = foodcols) +
  theme(
    panel.border = element_rect(linewidth = 0.5, colour = "black"),
    panel.spacing.x = unit(0, "lines"),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size=10),
    strip.background = element_blank(),
    strip.text.x = element_text(face="bold"),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.5)
  ) +
  #scale_y_continuous(limits = c(-1.5,1.5))+
  scale_x_discrete(labels=xlab)+
  labs(y = "Z-score", caption = "NAs imputed with population mean.")
diet.freqlines


walk(c("svg", "png", "pdf"), ~ {
  ggsave(
    filename = paste0("fig_out/03.dietary_assesment/food_freq_lines.", .x),
    plot = diet.freqlines, width = 8, height = 4, scale=1
  )
})


# barplot
bardata <- dff %>% # not imputed
  filter(population != "LDY") %>%
  select(sampleid, population, all_of(select_var)) %>%
  pivot_longer(
    cols = all_of(select_var), 
    names_to = "variable", 
    values_to = "freq"
  ) %>%
  group_by(variable) %>%
  mutate(
    Zscore = as.numeric(scale(freq))
  ) %>%
  ungroup() %>%
  mutate(
    type = factor(zscore_group[variable], levels = unique(zscore_group)),
    group = case_when(
      population %in% c("BTU", "ORT", "ORS") ~ "earlyTransition",
      population %in% c("APT", "TBU", "BSP") ~ "lateTransition",
      population == "LDY" ~ "Agri",
      TRUE ~ NA_character_
    )) %>%
  mutate(
    group = factor(group, levels = c("earlyTransition","lateTransition","Agri")),
    population = factor(population, levels = pop.ord),
    variable = factor(variable, levels = select_var),
    type = factor(type, levels = food_lvl)) %>%
  mutate(variable = factor(variable, levels = select_var),
        type = factor(type, levels = food_lvl),
         group = factor(group, c("earlyTransition","lateTransition")))

bardata %>%
  group_by(group) %>%
  summarise(n=n_distinct(sampleid))

ggdata_summary <- bardata %>%
  group_by(variable, group, type) %>%
  summarise(
    mean_freq = mean(freq, na.rm = TRUE),
    se_freq = sd(freq, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

dodge_val <- position_dodge(width = 0.8)

nGroup <- bardata %>%
  group_by(group) %>%
  summarise(n=n_distinct(sampleid)) %>%
  ungroup() %>%
  mutate(lab=paste0(group," (", n,")"))

legend_labs <- nGroup$lab
names(legend_labs) <- nGroup$group

pflip <- bardata %>%
  mutate(variable=factor(variable, rev(names(foodcols)))) %>%
  filter(!is.na(freq)) %>%
  ggplot(aes(x = variable, y = freq, fill = group)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.background = element_rect(fill = "white", linewidth = 0),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "grey70", linewidth = 0.5, linetype = "dashed"),
        strip.text = element_text(face = "bold"),
        panel.spacing.y = unit(0, "lines"), 
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  stat_summary(fun = mean, geom = "bar", 
               position = dodge_val, width = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = dodge_val, width = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = c(earlyTransition = "forestgreen", 
                               lateTransition = "maroon"), 
                    labels = legend_labs) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  coord_flip() +
  labs(
    title = "Mean Food Frequency by Group",
    y = "Consumption*",
    x = "Food Variable",
    fill = "Transition Group"
  ) -> pflip

# Filter and prepare label data
stat_test <- bardata %>%
  group_by(variable, type) %>%
  summarise(
    freq_early = median(freq[group == "earlyTransition"], na.rm = TRUE),
    freq_late = median(freq[group == "lateTransition"], na.rm = TRUE), 
    Zscore_early = median(Zscore[group == "earlyTransition"], na.rm = TRUE),
    Zscore_late = median(Zscore[group == "lateTransition"], na.rm = TRUE), 
    deltaZ = Zscore_late -  Zscore_early,
    # The tryCatch ensures that if a test fails (e.g., zero variance), 
    # the script returns NA instead of stopping.
    p_val = tryCatch(
      wilcox.test(freq ~ group)$p.value,
      error = function(e) return(NA)
    ),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_val, method = "fdr", n = length(select_var))) %>%
  mutate(p_label = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE ~ ""
  ),
  bracket = case_when(
    p_adj < 0.001 ~ "]",
    p_adj < 0.01  ~ "]",
    p_adj < 0.05  ~ "]",
    TRUE ~ ""
  ))

label_data <- bardata %>%
  group_by(variable, group, type) %>%
  summarise(
    avg = mean(freq, na.rm = TRUE),
    se = sd(freq, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  group_by(variable) %>%
  # Identify the furthest right point for significance placement
  slice_max(avg, n = 1, with_ties = FALSE) %>%
  mutate(max_pos = avg + se) %>%
  ungroup() %>%
  left_join(stat_test, by = c("variable", "type")) %>%
  filter(p_label != "")

# Add significance layers to the plot
p <- pflip + 
  # Significance Stars (***)
  geom_text(data = label_data, 
            aes(y = max_pos + 2, x = variable, label = p_label),
            inherit.aes = FALSE, 
            angle = 90,
            #hjust = 0,    
            #vjust = 0.5,
            size = 3.2,
            col = "black") + 
  # Brackets (])
  geom_text(data = label_data, 
            aes(y = max_pos + 1, x = variable, label = bracket),
            angle = 0,      # Kept at 0 as the ']' is already vertical
            inherit.aes = FALSE, 
            #hjust = 0, 
            #vjust = 0.5, 
            size = 4,
            col = "black")


plot_ffqbar <- p

walk(c("svg", "png", "pdf"), ~ {
  ggsave(
    filename = paste0("fig_out/03.dietary_assesment/ffq_twosided_barplot.", .x),
    plot = p, width = 4, height = 3.5, scale = 1.2
  )
})

p + coord_flip()

# Save Plots for Figure 3
library(patchwork)
plot_A <- plot_idds + theme(title=element_text(size=8, face="bold"), legend.position = "bottom")
plot_C <- plot_ffqbar + theme(title=element_text(size=8, face="bold"), legend.position = "bottom")
plot_D <- diet.freqlines + theme(title=element_text(size=8, face="bold"), legend.position = "bottom")

to_save=list(
  Fig3A=plot_A,
  Fig3C=plot_C,
  Fig3D=plot_D
)

for(fig in names(to_save)){
  outfile=paste0(fig,".rds")
  saveRDS(to_save[[fig]], file = file.path("fig_out/03.dietary_assesment/",outfile))
}

# save_sample_list
saveRDS(unique(idds$sampleid), "input_files/idds_sample_list.rds")


write.table(idds,file="out/idds_data_processed.tsv",quote = F,row.names = F,sep = "\t")
write.table(stat_test,file="out/ffZ_wilcox_out.tsv",quote = F,row.names = F,sep = "\t")
