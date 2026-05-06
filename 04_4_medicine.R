rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)


###########
# custom theme
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

###########
# load data


medicine_raw <- read.delim("input_files/Medicine_use_compiled - medicine_type_matrix_corrected.tsv")

# fix labels
fixed_labels <- colnames(medicine_raw)
fixed_labels <- gsub("anticholesterol","antihyperlipidemics", fixed_labels)
fixed_labels <- gsub("hypertension","antihypertension", fixed_labels)
colnames(medicine_raw) <- fixed_labels

age_sex <- read.csv("input_files/Anthropometry_Compiled_final - IMPUTED.csv") %>% 
  select(sampleid,population,age,sex) %>% rename(sid=sampleid)

# define population order
trg_HG <- c("BTU","ORT","ORS")
trd_HG <- c("APT","TBU","BSP")
pop.ord <- c(trg_HG, trd_HG)

# define age bins
age_sex <- subset(age_sex, age_sex$population %in% pop.ord)
bin_size <- 20
max_age  <- max(age_sex$age[age_sex$sid %in% medicine_raw$sid], na.rm = TRUE)
bins     <- c(18, seq(30, max_age + bin_size, by = bin_size))
labels    <- paste(head(bins, -1), head(bins, -1) + bin_size - 1, sep = "-")
labels[1] <- "18-29"

# merge
medicine <- age_sex %>%
  filter(population != "LDY") %>%
  left_join(medicine_raw, by=c("sid","population")) %>%
  # filter(!is.na(all_medicine_list)) %>%
  mutate(across(c(antibiotic:other), ~ replace_na(.x, 0)),
         all_medicine_list = replace_na(all_medicine_list, "no_report")) %>%
  mutate(group=ifelse(population %in% c("BTU","ORT","ORS"),"earlyTransition","lateTransition")) %>%
  mutate(age_bin = cut(
    age,
    breaks = bins,
    labels = labels,
    include.lowest = TRUE,
    right = FALSE
  )) %>%
  mutate(group=factor(group, c("earlyTransition","lateTransition")),
         sex=factor(sex, c("M","F")),
         population=factor(population, pop.ord))

# rank by total
medicine_rank <- colSums(medicine[,-c(1:5,16,17)])
medicine_rank <- medicine_rank[order(medicine_rank,decreasing = T)]
medicine_rank <- names(medicine_rank)
medicine_rank <- c(medicine_rank[medicine_rank!="other"],"other")

# generate colours
my_colors <- c("forestgreen","#70A470","#376B37", "maroon","#B7778F","#7E3E56")
names(my_colors) <- pop.ord

ggdata <- medicine %>%
  select(!all_medicine_list) %>%
  pivot_longer(cols=c("antibiotic","analgesic","respiratory","antihyperlipidemics","antigout","antihistamine","gastrointestinal","herbal","antihypertension","other"), 
               names_to = "medicine_type", values_to = "consumption") %>%
  filter(consumption == 1) %>%
  mutate(medicine_type=factor(medicine_type, medicine_rank))%>%
  filter(medicine_type != "other") %>%
  group_by(medicine_type,group,population,sex) %>%
  summarise(n=n_distinct(sid), .groups = "drop") 

pop_total <- age_sex %>%
  filter(sid %in% medicine$sid) %>%
  group_by(population,sex) %>%
  summarise(pop_total=n_distinct(sid), .groups = "drop")

ggdata <- ggdata %>%
  left_join(pop_total, by=c("population","sex")) %>%
  mutate(prop=n/pop_total) %>%
  mutate(population=factor(population, pop.ord),
         sex=factor(sex, c("M","F")))



p <- ggplot(ggdata, aes(x = medicine_type, y = prop, fill = population, group = population)) +
  guides(fill = guide_legend(nrow = 2, byrow = T)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9, preserve = "total")) +
  facet_grid(~sex, scales = "free_x", labeller = as_labeller(c("M" = "Men", "F" = "Women"))) +
  labs(
    x = NULL, 
    y = "Proportion", 
    title = "Medicine Usage") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = my_colors)
p

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","05.non_food","medicine_usage"),".",dev)
  ggsave(plot = p,filename = outpath, device = dev, width = 10,height = 5)
}

group_totals <- medicine %>%
  group_by(sex, group) %>%
  summarise(total_participants = n_distinct(sid), .groups = "drop")

stats_data <- medicine %>%
  select(!all_medicine_list) %>%
  pivot_longer(
    cols = c("antibiotic", "analgesic", "respiratory", "antihyperlipidemics", 
             "antigout", "antihistamine", "gastrointestinal", "herbal", 
             "antihypertension", "other"), 
    names_to = "medicine_type", 
    values_to = "consumption"
  ) %>%
  filter(consumption == 1, medicine_type != "other") %>%
  mutate(medicine_type = factor(medicine_type, levels = medicine_rank)) %>%
  group_by(medicine_type, group, sex) %>%
  summarise(n = n_distinct(sid), .groups = "drop")


fisher_results <- stats_data %>%
  left_join(group_totals, by = c("sex", "group")) %>%
  mutate(not_n = total_participants - n) %>%
  group_by(medicine_type, sex) %>%
  filter(n_distinct(group) == 2) %>%
  summarise(
    earlyTransition_yes  = n[1],
    earlyTransition_no   = not_n[1],
    lateTransition_yes  = n[2],
    lateTransition_no   = not_n[2],
    p_val = fisher.test(matrix(c(earlyTransition_yes, earlyTransition_no, lateTransition_yes, lateTransition_no), nrow = 2))$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_val, method = "fdr"))

fisher_results

write.table(fisher_results,file = "out/medicine_usage_fisher.tsv",quote = F,sep = "\t",row.names = FALSE)

message("note that fisher test of transitioning vs transitioned all not significant")


# to combine
panel_E <- p
saveRDS(panel_E, "fig_out/05.non_food/Panel_E.rds")



# save sample list
df <- medicine
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


saveRDS(unique(df$sid),file = "input_files/medicine_sample_list.rds")
write.table(df, "out/medicine_data_processed.tsv", quote = F, sep = "\t", row.names = F)