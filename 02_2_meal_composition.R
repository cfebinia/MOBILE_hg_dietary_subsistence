rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

meal_raw <- read.delim("input_files/Meal_Time_Composition - meal_matrix.tsv") %>% 
  rename(population = populaiton)
age_sex <- read.csv("input_files/Anthropometry_Compiled_final - IMPUTED.csv") %>% 
  select(sampleid,population,age,sex) %>% rename(sid=sampleid)

# define population order
trg_HG <- c("BTU","ORT","ORS")
trd_HG <- c("APT","TBU","BSP")
pop.ord <- c(trg_HG, trd_HG)

# define age bins
age_sex <- subset(age_sex, age_sex$population %in% pop.ord)
bin_size <- 20
max_age  <- max(age_sex$age[age_sex$sid %in% meal_raw$sid], na.rm = TRUE)
bins     <- c(18, seq(30, max_age + bin_size, by = bin_size))
labels    <- paste(head(bins, -1), head(bins, -1) + bin_size - 1, sep = "-")
labels[1] <- "18-29"

# merge
essentials <- c("rice","noodle","white_meat","red_meat","fish","eggs","vegetables","fruit","legumes")
meal <- age_sex %>%
  filter(population != "LDY") %>%
  left_join(meal_raw[,c("sid","population",essentials)], by=c("sid","population")) %>%
  filter(sid %in% meal_raw$sid) %>%
  mutate(total_items=rowSums(select(., all_of(essentials)), na.rm = TRUE)) %>%
  filter(total_items > 0) %>%
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

ggdata <- meal %>%
  select(sid, group, population, sex, total_items)

gcol <- c(earlyTransition="forestgreen", lateTransition="maroon")

ggplot(ggdata, aes(x = population, y = total_items, col = group)) +
  theme_bw() +
  # geom_boxplot(aes(alpha = sex), position = position_dodge(width = 0.8), outlier.shape = NA, fill=NA) +
  geom_violin(
    aes(fill=group, group=interaction(population, sex)), alpha=0.2,
    position = position_dodge(width = 0.8), 
    trim = TRUE,    # Keeps the violin within the data range
    scale = "width" # Makes all violins the same maximum width
  ) +
  geom_point(aes(shape = sex), size=2.5, 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5, jitter.height = 0.1), 
             alpha = 0.5) +
  scale_alpha_manual(values = c("M" = 1, "F" = 0.6)) + # Slight transparency for females to distinguish
  labs(y = "Total Essential Items Consumed", x = "Population Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = gcol) +
  scale_color_manual(values = gcol)

set.seed(2)
p <- ggplot(ggdata, aes(x = group, y = total_items, col = group)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")+
  geom_violin(
    aes(fill=group, group=interaction(group, sex)), alpha=0.2,
    position = position_dodge(width = 0.8), 
    trim = TRUE,    # Keeps the violin within the data range
    scale = "width" # Makes all violins the same maximum width
  ) +
 geom_point(aes(shape = sex), size=2.5, 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5, jitter.height = 0.1), 
             alpha = 0.5) +
  scale_alpha_manual(values = c("M" = 1, "F" = 0.6)) + # Slight transparency for females to distinguish
  labs(title = "Items Consumed per Meal", x = "", y = "") +
  scale_fill_manual(values = gcol) +
  scale_color_manual(values = gcol) +
  scale_shape_manual(values = c(M=15, F=16))

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","03.dietary_assesment","meal_composition","items_per_meal"),".",dev)
  ggsave(plot = p,filename = outpath, device = dev, width = 6.4,height = 4.5,dpi = 300, scale=0.8)
}


# most_paired
ggdata2 <- meal %>%
  pivot_longer(cols=all_of(essentials), 
               names_to = "food_type", values_to = "consumption") %>%
  filter(consumption == 1) %>%
  group_by(sid,group,population,age,sex) %>%
  summarise(meal_composition = paste(food_type, collapse = ","),
            n_items = n(),
            .groups = "drop")

meal_composition_freq <- ggdata2 %>%
  group_by(meal_composition) %>%
  summarise(freq=n(), .groups = "drop") %>%
  arrange(desc(freq)) 
meal_composition_rank <- meal_composition_freq$meal_composition
top_meal_comp <- meal_composition_rank[1:15]

total_population <- age_sex %>%
  filter(sid %in% meal_raw$sid) %>%
  group_by(population,sex) %>%
  summarise(pop_total=n_distinct(sid))

sample_list <- unique(ggdata2$sid)
ggdata2 %>%
  group_by(population) %>%
  summarise(n=n_distinct(sid)) %>%
  mutate(xlab=paste0(population," (", n, ")")) -> n_pop
xlab <- n_pop$xlab
names(xlab) <- n_pop$population

meal_compos <- ggdata2 %>%
  group_by(group,population,sex,meal_composition) %>%
  summarise(freq=n_distinct(sid), .groups = "drop") %>%
  left_join(total_population, by=c("population","sex")) %>%
  mutate(prop=freq/pop_total) %>%
  mutate(meal_composition=factor(meal_composition,meal_composition_rank),
         population=factor(population, pop.ord),
         sex=factor(sex, c("M","F"))) 

ggdata3 <- meal_compos %>%
  filter(meal_composition %in% top_meal_comp)

my_colors <- c("forestgreen","#70A470","#376B37", "maroon","#B7778F","#7E3E56")
names(my_colors) <- pop.ord

p <- ggplot(ggdata3, aes(x = meal_composition, y = prop, fill = population, group = population)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")+
  geom_bar(stat = "identity", position = position_dodge(width = 0.8,preserve = "total")) +
  facet_grid(group~ sex) +
  guides(fill = guide_legend(nrow = 1, byrow = T)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(labels = scales::percent)+
  labs(title = "Common Meal Compositions (top 15)", y = "proportion", x = "")
p  

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","03.dietary_assesment","meal_composition","meal_composition_top"),".",dev)
  ggsave(plot = p,filename = outpath, device = dev, width = 10,height = 7,dpi = 300, scale=0.8)
}

meal_composition_freq <- ggdata2 %>%
  group_by(meal_composition, population, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = c(population, sex), 
    values_from = n, 
    values_fill = 0,
    names_sep = "_"
  ) %>%
  mutate(
    Total_earlyTransition=rowSums(select(., starts_with(trg_HG))),
    Total_transitioninedHG=rowSums(select(., starts_with(trd_HG))),
    BTU_Both = rowSums(select(., starts_with("BTY_"))),
    ORT_Both = rowSums(select(., starts_with("ORT_"))),
    ORS_Both = rowSums(select(., starts_with("ORS_"))),
    APT_Both = rowSums(select(., starts_with("APT_"))),
    TBU_Both = rowSums(select(., starts_with("TBU_"))),
    BSP_Both = rowSums(select(., starts_with("BSP_"))),
    Total_Male   = rowSums(select(., ends_with("_M"))),
    Total_Female = rowSums(select(., ends_with("_F"))),
    Grand_Total  = Total_Male + Total_Female
  ) %>%
 arrange(desc(Grand_Total))

write.table(meal_composition_freq,file = "out/meal_composition_byPopSex.tsv",quote = F,sep = "\t",row.names = FALSE)

# By rank visualisation
meal_compos_ranked <- meal_compos %>%
  group_by(population, sex) %>%
  mutate(rank = dense_rank(desc(freq))) %>%
  arrange(population, sex, rank) %>%
  ungroup()

meal_compos_ranked_tops <- 
  meal_compos_ranked %>%
  group_by(population, sex) %>%
  # Get top ranked (5)
  slice_max(order_by = freq, n = 5, with_ties = FALSE) %>%  
  arrange(population, sex, desc(freq))
tops_list <- unique(as.character(meal_compos_ranked_tops$meal_composition))
meal_compos_ranked_tops <-
  meal_compos_ranked %>%
  filter(meal_composition %in% tops_list)
not_listed <- as.character(meal_compos$meal_composition)[!as.character(meal_compos$meal_composition) %in% tops_list]
head(meal_compos_ranked_tops)
print(not_listed)

# All

plotdata <- meal_compos_ranked %>% # or can change to the meal_compos_ranked_tops for the top5
  filter(meal_composition %in% meal_composition_rank[1:30])
  
plot_meal_rank <- plotdata %>%
  ggplot(aes(x = sex, y = meal_composition)) + 
  # geom_point(aes(size = prop, color = sex), alpha=0.6) +
  # scale_color_manual(values = c("M" = "#0072B2", "F" = "#CC79A7"))
  geom_point(aes(size = prop, color = group), alpha=0.6) +
  scale_colour_manual(values=gcol) +
  facet_grid(~sex+population, scale="free_x") +
  scale_y_discrete(limits = rev, expand = expansion(mult = c(0.05, 0.05))) +
  scale_size_continuous(range = c(1, 10)) + 
  theme_bw(base_size = 14) +
  labs(title = "Meal Compositions Across Communities",
       x = "Sex",
       y = "Meal Composition",
       size = "Proportion") +
  theme(panel.grid.major = element_line(linewidth = 0.5, color = "grey90"),
        panel.spacing.x = unit(0, "lines"),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

plot_meal_rank

length(unique(as.character(meal_compos_ranked$meal_composition)))

for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","03.dietary_assesment","meal_composition","meal_composition_prop_ranked"),".",dev)
  ggsave(plot = plot_meal_rank,filename = outpath, device = dev, width = 10,height = 9,dpi = 300, scale=0.8)
}

#######
# MDS
library(tidyverse)
library(vegan)
library(ggdendro)

sig.level = 0.1 # max pvalue for significant assoc

for(db in c("rank","proportion")){
# distance_by <- "proportion" # options: "rank" or "proportion"
  distance_by <- db 
  
if(distance_by == "rank") {
  penalty_val <- max(meal_compos_ranked$rank, na.rm = TRUE) + 1
  value_col <- "rank"
} else {
  raw_min <- min(meal_compos_ranked$prop, na.rm = TRUE)
  penalty_val <- 10^floor(log10(raw_min))
  value_col <- "prop"
}

process_gender_data <- function(df, g, v_col, p_val) {
  mat <- df %>%
    filter(sex == g) %>% 
    select(population, meal_composition, !!sym(v_col)) %>%
    pivot_wider(names_from = population, values_from = !!sym(v_col), 
                values_fill = list(rank = p_val, prop = p_val)) %>%
    filter(meal_composition%in%meal_composition_rank[1:30])
  
  dists <- mat %>%
    select(-meal_composition) %>%
    as.matrix() %>%
    t() %>%
    dist(method = "euclidean")
  
  return(list(matrix = mat, dists = dists))
}

resM <- process_gender_data(meal_compos_ranked, "M", value_col, penalty_val)
resF <- process_gender_data(meal_compos_ranked, "F", value_col, penalty_val)

distM <- resM$dists
distF <- resF$dists

plot_hc <- function(dists, lab_title, direction="right") {
  hc <- hclust(dists, method = "ward.D2")
  dend_data <- dendro_data(hc)
  
  tmp <- dend_data$labels
  tmp$group <- ifelse(tmp$label %in% c("BTU","ORT","ORS"), "earlyTransition","lateTransition")
  dend_data$labels <- tmp
  
    if(direction =="right"){
    segmt <- segment(dend_data)
  } else if(direction == "left"){
    segmt <- segment(dend_data)
    segmt[,c("y","yend")] <- segmt[,c("y","yend")] * -1
  }

  figout <- ggplot(segmt) + 
    theme_void() +
    theme(plot.background = element_rect(fill="white", linewidth = 0),
          legend.position = "none")+
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    scale_colour_manual(values = gcol)+
    labs(title = lab_title)
  
  if(direction == "right"){
  figout <- figout +
    geom_text(data = label(dend_data), aes(x = x, y = y-0.02, label = label, colour = group), 
              hjust = 1, angle = 0, size = 4, fontface="bold") +
    coord_flip()
  } else if(direction == "left"){
    figout <- figout +
      geom_text(data = label(dend_data), aes(x = x, y = y+0.02, label = label, colour = group), 
                hjust = 0, angle = 0, size = 4, fontface="bold") +
      coord_flip()
  
    return(figout)
  }
  }

p_dendro_M <- plot_hc(distM, direction = "left", paste0("Meal Composition Similarity in Men  (", distance_by, ")"))
p_dendro_M

figname=paste("meal_composition_dendro",distance_by,"M",sep="_")
for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","03.dietary_assesment","meal_composition",figname),".",dev)
  ggsave(plot = p_dendro_M,filename = outpath, device = dev, width = 4,height = 2)
}

p_dendro_F <- plot_hc(distF, direction = "right", paste0("Meal Composition Similarity in Women (", distance_by, ")"))
p_dendro_F

figname=paste("meal_composition_dendro",distance_by,"F",sep="_")
for(dev in c("svg","png","pdf")){
  outpath=paste0(file.path(getwd(),"fig_out","03.dietary_assesment","meal_composition",figname),".",dev)
  ggsave(plot = p_dendro_F,filename = outpath, device = dev, width = 4,height = 2)
}

if(db=="rank") {hc_plots <- list(M=p_dendro_M, F=p_dendro_F)}

## MDS 
ordM <- monoMDS(distM, stress = 1)
ordF <- monoMDS(distF, stress = 1)
all_points <- rbind(scores(ordM, display="sites"), scores(ordF, display="sites"))
max_val <- max(abs(all_points)) * 1.05 # Add 20% padding for labels
plot_limit <- c(-max_val, max_val)

for(i in c("M","F")){
  if(i == "M"){
    imat <- resM$matrix
    idist <- distM
    ititle <- paste0("Meal Composition Similarity in Men (", distance_by, ")")
  } else {
    imat <- resF$matrix
    idist <- distF
    ititle <- paste0("Meal Composition Similarity in Women (", distance_by, ")")
  }
  
  # set ordination
  set.seed(1)
  ord <- monoMDS(idist, stress = 1)
  site_coords <- as.data.frame(scores(ord, display = "sites"))
  site_coords$population <- rownames(site_coords)
  site_coords$group <- ifelse(site_coords$population %in% c("BTU","ORT","ORS"),
                              "earlyTransition", "lateTransition")
  
  # get variance explained
  fit <- goodness(ord)
  dist_orig <- as.vector(idist)
  dist_nmds1 <- dist(scores(ord, display="sites")[,1])
  r2_axis1 <- cor(dist_orig, as.vector(dist_nmds1))^2
  dist_nmds2 <- dist(scores(ord, display="sites")[,2])
  r2_axis2 <- cor(dist_orig, as.vector(dist_nmds2))^2
  xlab <- paste0("NMDS1 (Proxy: ", round(r2_axis1 * 100, 1), "% variance)")
  ylab <- paste0("NMDS2 (Proxy: ", round(r2_axis2 * 100, 1), "% variance)")
  
  # get environmental fitting
  env_data <- imat %>%
    select(-meal_composition) %>%
    as.data.frame()
  rownames(env_data) <- imat$meal_composition
  en <- envfit(ord, as.data.frame(t(env_data)), permutations = 9999, na.rm = TRUE)
  vector_data <- as.data.frame(scores(en, display = "vectors"))
  vector_data$p_val <- en$vectors$pvals
  vector_data$r2 <- en$vectors$r
  vector_data$taster <- rownames(vector_data)
  
  sig_vectors <- subset(vector_data, p_val < sig.level)
  
  # Safety check for scaling
  if(nrow(sig_vectors) > 0) {
    scaling_factor <- 0.7
    sig_vectors$NMDS1 <- sig_vectors$NMDS1 * sqrt(sig_vectors$r2) * scaling_factor
    sig_vectors$NMDS2 <- sig_vectors$NMDS2 * sqrt(sig_vectors$r2) * scaling_factor
  }
  
  p <- ggplot() +
    theme_bw(base_size = 14) +
    coord_fixed(xlim = plot_limit, ylim = plot_limit)+
    theme(panel.grid = element_blank(), 
          legend.position = "bottom",
          plot.margin = margin(10, 10, 10, 10)) +
    geom_vline(xintercept = 0, colour="grey70", linewidth=0.5, linetype="dashed") +
    geom_hline(yintercept = 0, colour="grey70", linewidth=0.5, linetype="dashed") +
    geom_point(data = site_coords, aes(x = NMDS1, y = NMDS2, colour = group), 
               size = 4, alpha = 0.7) +
    geom_text(data = site_coords, aes(x = NMDS1, y = NMDS2, label = population), 
              vjust = 1.5, size = 3, color = "black") +
    scale_color_manual(values=c(earlyTransition="forestgreen", lateTransition="maroon")) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(title = ititle,
         subtitle = "Displaying only associated variables (p<0.1, 9999 permutations)",
         caption = paste("Stress:", round(ord$stress, 3)),
         x = xlab, y = ylab)
  
  # Check if we have significant vectors to plot
  if(nrow(sig_vectors) > 0) {
    p <- p + 
      geom_segment(data = sig_vectors, 
                   aes(x = 0, y = 0, xend = -NMDS1, yend = -NMDS2),
                   arrow = arrow(length = unit(0.25, "cm")), 
                   color = "firebrick", linewidth = 0.8) +
      geom_text(data = sig_vectors, size = 3,
                aes(x = -NMDS1, y = -NMDS2, label = taster), 
                color = "firebrick", fontface = "bold", vjust = -0.5, hjust = 1.1) +
      labs(subtitle = paste0("Displaying significant meal drivers: p<",sig.level))
    
  } else {
    # If no variables are significant, add a message to the plot
    p <- p + 
      labs(subtitle = "No variables met the significance threshold")
    
    # Also print a message to the console for the operations lead
    message(paste("Notice: No significant environmental variables found for gender:", i))
  }
  
  print(p)
  
  if(distance_by == "rank" & i == "M") {
    mds_plots <- list()
    mds_plots[[i]] <- p + coord_fixed(xlim = plot_limit, ylim = plot_limit, ratio = 0.7)
  } else if(distance_by == "rank" & i == "F") {
    mds_plots[[i]] <- p + coord_fixed(xlim = plot_limit, ylim = plot_limit, ratio = 0.7)
    }
  
  figname=paste("meal_composition_NMDS",distance_by,i,sep="_")
  for(dev in c("svg","png","pdf")){
    outpath=paste0(file.path(getwd(),"fig_out","03.dietary_assesment","meal_composition",figname),".",dev)
    ggsave(plot = p ,filename = outpath, device = dev, width = 8,height = 8)
  }
}
}

## Compare rank
meal_compos_ranked$geo <- ifelse(meal_compos_ranked$population %in% c("ORT","ORS"),
                                 "Sumatra","Borneo")

# by sex
wilcox_data <- meal_compos_ranked %>%
  group_by(meal_composition, population) %>%
  filter(n_distinct(sex) == 2) %>%
  ungroup() %>%
  group_by(meal_composition) %>%
  filter(n_distinct(group) == 2) %>%
  ungroup()

var <- unique(as.character(wilcox_data$meal_composition))
res_list <- list()
for(i in var){
  df <- subset(wilcox_data, wilcox_data$meal_composition==i)
  df <- df[,c("prop","sex")]
  colnames(df) <- c("x","g")
  res <- wilcox.test(df$x ~ df$g, exact=F)
  res_list[[i]] <- res
}



distance_by <- db 
penalty_val <- max(meal_compos_ranked$rank, na.rm = TRUE) + 1
value_col <- "rank"

resM <- process_gender_data(meal_compos_ranked, "M", value_col, penalty_val)
resF <- process_gender_data(meal_compos_ranked, "F", value_col, penalty_val)

compare_rankM <- 
  resM$matrix %>%
  pivot_longer(cols=all_of(as.character(pop.ord)),names_to = "population",values_to = "rank") %>%
  mutate(group = ifelse(population %in% c("BTU","ORS","ORT"),
                        "earlyTransition","lateTransition"),
         geo = ifelse(population %in% c("ORT","ORS"),
                      "Sumatra","Borneo"))

compare_rankF <- 
  resF$matrix %>%
  pivot_longer(cols=all_of(as.character(pop.ord)),names_to = "population",values_to = "rank") %>%
  mutate(group = ifelse(population %in% c("BTU","ORS","ORT"),
                        "earlyTransition","lateTransition"),
         geo = ifelse(population %in% c("ORT","ORS"),
                      "Sumatra","Borneo"))

compare_rankSum <- function(data, group, signif.level=0.05){
  var <- unique(data$meal_composition)
  res_list <- list()
  for(i in var){
    df <- subset(data, data$meal_composition==i)
    df <- df[,c("rank",group)]
    colnames(df) <- c("x","g")
    res <- wilcox.test(df$x ~ df$g, exact=F)
    is_signif <- (res$p.value < signif.level)
    
    if(isTRUE(is_signif)){res_list[i] <- res$p.value}
  }
  
  if(length(res_list)>0){
    return(res_list)}else{message("no significant difference")}
}

compare_rankSum(compare_rankM, group = "group", signif.level = 0.1)
compare_rankSum(compare_rankF, group = "group", signif.level = 0.1)

compare_rankSum(compare_rankM, group = "geo", signif.level = 0.1)
compare_rankSum(compare_rankF, group = "geo", signif.level = 0.1)


compare_sex <-
  rbind(data.frame(compare_rankM, sex="M"),
      data.frame(compare_rankF, sex="F")) %>%
  group_by(meal_composition, population) %>%
  # Filter for combinations present in both Male and Female groups
  filter(n_distinct(sex) == 2) %>%
  ungroup()

compare_rankSum(subset(compare_sex, geo=="Borneo"), group = "sex", signif.level = 0.1)
compare_rankSum(subset(compare_sex, geo=="Sumatra"), group = "sex", signif.level = 0.1)
compare_rankSum(subset(compare_sex, group=="earlyTransition"), group = "sex", signif.level = 0.1)
compare_rankSum(subset(compare_sex, group=="lateTransition"), group = "sex", signif.level = 0.1)


# Build Figure S2
library(patchwork)

plot_A <- plot_meal_rank 
plot_B1 <- hc_plots[["M"]] + theme(legend.position = "none") + labs(title="") + scale_y_continuous(expand = expansion(mult = c(0, 0.7)))
plot_B2 <- mds_plots[["M"]] + 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) + 
  labs(title="Men", subtitle = NULL)
plot_B3 <- mds_plots[["F"]] + 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) + 
  labs(title="Women", subtitle = NULL)
plot_B4 <- hc_plots[["F"]] + theme(legend.position = "none")  + labs(title="") + scale_y_continuous(expand = expansion(mult = c(0.7, 0)))


row1 <- plot_A
row2 <- (plot_B1 | plot_B2 | plot_B3 | plot_B4) + 
  plot_layout(widths = c(0.1, 0.4,0.4, 0.1), guides = "collect") &
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))
  
figureS2 <- (wrap_elements(row1) / row2) + 
  plot_layout(heights = c(7, 3)) + 
  plot_annotation(tag_levels = list(c("A", "B"))) & 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.size = unit(4, "mm"),
    plot.tag = element_text(size = 20, face = "bold")
  )


figureS2

figout="fig_out/03.dietary_assesment/FigureS2_mealcomposition"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = figureS2, filename = outpath, device = dev, 
         width = 210, height = 260, units = "mm", dpi = 300,
         scale = 1.2)
}

# save_sample_list
saveRDS(unique(ggdata2$sid), "input_files/meal_compos_sample_list.rds")
saveRDS(unique(ggdata2$sid[ggdata2$sex=="M"]), "input_files/meal_compos_M_sample_list.rds")
saveRDS(unique(ggdata2$sid[ggdata2$sex=="F"]), "input_files/meal_compos_F_sample_list.rds")

write.table(ggdata2, file="out/meal_compos_data_processed.tsv",quote = F,sep = "\t",row.names = FALSE)

