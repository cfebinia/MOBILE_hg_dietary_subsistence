rm(list=ls())
library(dplyr)
library(tidyverse)
library(ggplot2)


# load data
df.src <- read.delim("input_files/food_source_corrected.tsv")
#src.pc.byid <- read.delim("input_files/Food_Source_Compiled_final - IDPROP.tsv")
src.pc.byid <- read.delim("input_files/fsource_corrected_imputed_median.tsv") 
dff <- read.delim("input_files/ffq_corrected_imputed.tsv") # using ffq imputed data

pop.ord <- c("BTU","ORT","ORS","APT","TBU","BSP")
var <- colnames(dff)[-c(1,2)]
food.var <- var[!var%in%c("water","coffee_tea","added_sugar")]

# load colour themes
mycols <- readRDS("input_files/mycols.rds")
popcols <- mycols$col
names(popcols) <- mycols$population

foodcols <- rainbow(length(colnames(dff)[-c(1,2)]),s = 0.6,v = 0.8)
set.seed(1) ; foodcols <- sample(foodcols)
names(foodcols) <- colnames(dff)[-c(1,2)]

# sample count
samples <- df.src$sampleid
length(samples)
table(substr(samples,1,3))

##############
## individual boxplot
boxdata <-src.pc.byid

n_pop <- boxdata %>%
  group_by(population) %>%
  summarise(n=n_distinct(sampleid)) %>%
  mutate(xlab=paste0(population,"\n(",n,")"))
xlab <- n_pop$xlab
names(xlab) <- n_pop$population

boxdata <- boxdata %>%
  pivot_longer(cols=c(wild,market,garden), names_to="source", values_to="proportion") %>%
  mutate(group=case_when(
    population %in% c("BTU","ORT","ORS")  ~ "Early-transition",
    population %in% c("APT","TBU","BSP") ~ "Late-transition",
    population == "LDY" ~ "Agri.",
    TRUE ~ NA_character_)) %>%
  mutate(source=factor(source, c("wild","market","garden")),
         population=factor(population, c(pop.ord,"LDY")),
         group=factor(group, c("Early-transition","Late-transition", "Agri.")))

p <- ggplot(boxdata, aes(y = proportion, x = population, colour = population)) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  geom_point(aes(shape=group), size = 2, position = position_jitter(seed = 1, width = 0.2)) +
  geom_vline(xintercept = c(3.5, 6.5), linewidth = 0.8, col = "grey40", linetype = 2) +
  facet_wrap(~source, ncol = 3) +
  scale_color_manual(values = popcols, labels=xlab) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = scales::percent) +
  scale_shape_manual(values = c(16,15,17))+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face="bold"),
        strip.background = element_blank(),
        legend.position = "top") +
  guides(colour = guide_legend(nrow = 1, reverse = TRUE)) +
  labs(title = "Frequency of wild, garden, and market goods in individual diet",
       caption = "freq of wild + garden + market = 1, across 14 items",
       x = "Population",
       y = "Proportion")

p

figout <- "fig_out/02.food_source/source_proportion_boxplot"
for(dev in c("png","svg","pdf")){
  fileout <- paste0(figout,".",dev)
  ggsave(plot=p, filename = fileout, device = dev,width = 8,height = 3.8,units = "in",dpi = 600,scale = 1.2)
}


boxdata %>%
  group_by(group, population) %>%
  filter(source=="wild") %>%
  summarise(median=median(proportion[proportion>0], na.rm=T),
            #iqr = IQR(proportion, na.rm = T),
            q1 = quantile(proportion, 0.25, na.rm = TRUE),
            q3 = quantile(proportion, 0.75, na.rm = TRUE),
            .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2)))

boxdata %>%
  group_by(group) %>%
  filter(source=="garden") %>%
  summarise(median=median(proportion[proportion>0], na.rm=T),
            #iqr = IQR(proportion, na.rm = T),
            q1 = quantile(proportion, 0.25, na.rm = TRUE),
            q3 = quantile(proportion, 0.75, na.rm = TRUE),
            .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2)))

boxdata %>%
  filter(group != "Agri.") %>%
  filter(source == "garden") %>%
  wilcox.test(proportion ~ group, data = .)

boxdata %>%
  summarise(n=n_distinct(sampleid))

boxdata %>% 
  group_by(group) %>%
  summarise(n=n_distinct(sampleid))

boxdata %>% 
  group_by(population) %>%
  summarise(n=n_distinct(sampleid))

# Ternary Plot
library("Ternary")
library("animation")
pop.ord0 <- c(pop.ord,"LDY")

df3d <- src.pc.byid

transpcol <- adjustcolor(popcols, alpha.f = 1)
names(transpcol) <- names(popcols)

dftr <- data.frame(df3d, col=transpcol[as.character(df3d$population)])
dftr$population<-factor(dftr$population,pop.ord0)
dftr <- dftr[order(dftr$population, decreasing = F),]
dftr$group <-  case_when(dftr$population %in% c("BTU","ORT","ORS")  ~ "ET",
                         dftr$population %in% c("APT","TBU","BSP") ~ "LT",
                         dftr$population == "LDY" ~ "AG")
dftr$group <- factor(dftr$group, c("ET","LT","AG"))  
myshape <- c(ET=16,LT=15,AG=17)
popshape <- c(16,16,16,15,15,15,17)
names(popshape) <- levels(dftr$population)
tr_data <- list()
saveGIF({
  i<-""
  outfile <- paste0("fig_out/food_source_ternary",i,collapse = "")
  #png(filename = paste0(outfile,".png"),res=600, width=7, height=7,units = "in")
  #svg(filename = paste0(outfile,".svg"), width=7, height=7)
  par(mar = c(1, 1, 3, 1)) 
  TernaryPlot( cex=1.5,
               atip = "Market", btip = "Garden", ctip = "Wild",
               alab = "%", blab = "%", clab = "%",
               grid.lines = 5, grid.lty = "dotted",
               grid.minor.lines = 1, grid.minor.lty = "dotted"
  )
  for(g in levels(dftr$group)){
    iShape <- myshape[g]
    ipoints <- dftr[dftr$group == g, c("market", "garden", "wild")]
    icols <- dftr$col[dftr$group == g]
    TernaryPoints(ipoints,
                  pch = iShape, 
                  cex = 2,
                  col = icols)
  }
  legend("topright",
         legend = names(popcols)[names(popcols)%in%dftr$population],
         col = popcols[names(popcols)%in%dftr$population],
         pch = popshape,
         pt.cex = 1.8,cex = 0.8,
         bty = "n") 
  title("Food Source", cex.main=1.5)
  #dev.off()
  
  # per pop
  for(i in levels(dftr$population)){
    outfile <- paste0("fig_out/food_source_ternary",i,collapse = "")
    idat <- dftr[dftr$population==i,]
    len <- length(idat$sampleid)
    icol <- transpcol[i]
    
    # save data as list
    # ni <- which(levels(dftr$population) == i)
    tr_data[[i]] <- idat
    
    # plot
    par(mar = c(1, 1, 3, 1)) 
    TernaryPlot( cex=1.5,
                 atip = "Market", btip = "Garden", ctip = "Wild",
                 alab = "%", blab = "%", clab = "%",
                 grid.lines = 5, grid.lty = "dotted",
                 grid.minor.lines = 1, grid.minor.lty = "dotted"
    )
    
    iShape=popshape[i]
    TernaryPoints(idat[,c("market","garden", "wild")],
                  pch = iShape, cex=2,
                  col = icol)
    title(cex.main=1.5, main=paste0(i, "\n(n=",len,")"))
    #dev.off()
  }
}, movie.name=file.path(getwd(),"fig_out/food_source_trenary.gif"), 
interval=1, ani.width=600, ani.height=600)

#### Set output for static plot
outfile=paste0("fig_out/02.food_source/food_source_ternary_","centroid",collapse = "")
pdf(file = paste0(outfile,".pdf"),width=7.5, height=7.5)
# png and svg does not import well
# png(filename = paste0(outfile,".png"),res=600, width=7.5, height=7.5,units = "in")
# svg(filename = paste0(outfile,".svg"), width=7.5, height=7.5)
par(mar = c(1, 1, 3, 1)) 

# Base plot
select_pop <- pop.ord0
TernaryPlot( cex=1.5,
             atip = "Market", btip = "Garden", ctip = "Wild",
             alab = "%", blab = "%", clab = "%",
             grid.lines = 5, grid.lty = "dotted",
             grid.minor.lines = 1, grid.minor.lty = "dotted"
)
legend(
  "topright",
  legend = select_pop,
  col = popcols[select_pop],
  y.intersp = 1.2,
  x.intersp = 1.2,
  pt.cex = 2.5, bty = "n", 
  pch = popshape, border = NA,fill = NA)

popshape2 <- c(21,21,21,22,22,22,24)
names(popshape2) <- levels(dftr$population)

for(i in select_pop){
  idat <- tr_data[[i]]
  icol <- popcols[i]
  iShape <- popshape2[i]
  fill_col <- adjustcolor(icol, alpha.f = 0.3)
  TernaryPoints(idat[,c("market","garden", "wild")],
                pch = iShape, cex = 1.2,
                col = icol, bg = fill_col)
}

# Centroids
n <- length(select_pop)
connect_centres <- matrix(0,nrow = n, ncol = 3,
                          dimnames = list(select_pop,
                                          c("market","garden", "wild")))

for(i in select_pop){
  idat <- tr_data[[i]]
  icol <- popcols[i]
  centres <- colMeans(idat[,c("market","garden", "wild")])
  connect_centres[i,] <- centres
}

# Connecting Lines
set1 <-connect_centres[-which(rownames(connect_centres) %in% c("BTU","ORS","ORT","LDY")),]
TernaryArrows(
  from = set1[-nrow(set1), c("market","garden", "wild")],
  to   = set1[-1, c("market","garden", "wild")],
  length = 0,
  lwd = 5,
  col="maroon"
)

set2 <-connect_centres[which(rownames(connect_centres) %in% c("BTU","ORS","ORT")),]
TernaryArrows(
  from = set2[-nrow(set2), c("market","garden", "wild")],
  to   = set2[-1, c("market","garden", "wild")],
  length = 0,
  lwd = 5,
  col="forestgreen"
)

# Points
for(i in select_pop){
  idat <- tr_data[[i]]
  icol <- popcols[i]
  iShape <- popshape2[i]
  centres <- colMeans(idat[,c("market","garden", "wild")])
  connect_centres[i,] <- centres
  TernaryPoints(centres,
                pch = iShape, cex = 3,
                bg = "white")
  TernaryPoints(centres,
                pch = iShape, cex = 2.2,
                col = "black", bg = icol)
}

dev.off()


#####################
### By Item procurement
# food order
food_reorder <- c("fish","white_meat","red_meat","honey",
                  "rice","legumes",
                  "eggs","noodles","cooking_oil",
                  "tuber","vegetables","fruits",
                  "lard","sago")

# plot raw data
src_bar <- df.src %>%
  pivot_longer(cols=all_of(food.var), names_to = "food", values_to = "source")%>%
  mutate(population=substr(sampleid,1,3),
         group=case_when(
           population %in% c("BTU","ORT","ORS")  ~ "Early-transition",
           population %in% c("APT","TBU","BSP") ~ "Late-transition",
           population == "LDY" ~ "Agri.",
           TRUE ~ NA_character_),
         island=case_when(
           population %in% c("ORT","ORS") ~ "Sumatra",
           population %in% c("BTU","APT","TBU","BSP","LDY") ~ "Borneo",
           TRUE ~ NA_character_))
  
src_bar %>%
  summarise(n=n_distinct(sampleid))

sample_list <- unique(src_bar$sampleid)

src_bar %>%
  group_by(group) %>%
  summarise(n=n_distinct(sampleid))

src_bar %>%
  group_by(population) %>%
  summarise(n=n_distinct(sampleid))


src_bar <- src_bar%>%   
  group_by(group, island, population, food) %>%
  summarise(n=n_distinct(sampleid),
            market=mean(!is.na(source) & source=="market"),
            garden=mean(!is.na(source) & source=="garden"),
            wild=mean(!is.na(source) & source=="wild"),
            not_consumed=mean(!is.na(source) & source=="not_consumed"),
            isNA=sum(is.na(source)/n)) %>%
  ungroup() %>%
  pivot_longer(cols=c(wild,market,garden, not_consumed, isNA), names_to="source", values_to="proportion") %>%
  mutate(xlab=paste0(population, "(",n,")")) %>%
  mutate(food=factor(food, food_reorder),
         source=factor(source, c("wild","market","garden","not_consumed","isNA")),
         population=factor(population, c(pop.ord,"LDY")),
         group=factor(group, c("Early-transition","Late-transition", "Agri.")))

distinct(src_bar[,c("population","xlab")]) -> xlab_df
population_labels <- setNames(as.character(xlab_df$xlab), xlab_df$population)        

p <- ggplot(src_bar) +
  theme_bw(base_size=10) +
  geom_bar(aes(x=population, y=proportion, fill=source), stat="identity", width = 1)+
  facet_grid(food~group, scales="free_x", space = "free_x")+
  scale_fill_manual(values=c(market="gold2",
                             garden="royalblue1",
                             wild="springgreen3",
                             not_consumed="grey50",
                             isNA="grey80"))+
  scale_x_discrete(expand = c(0,0), labels = population_labels) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.5,1), labels = scales::percent)+
  theme(panel.grid = element_blank())+
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(0, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.text.x = element_text(face="bold"),
        strip.text.y = element_text(angle=-90, hjust=0.5, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1,size=7))+
  theme(legend.position = "top")+
  geom_hline(yintercept = 0.5, linetype="dashed", col="black")+
  labs(title="Proportion of food sources", x="Population", y="Proportion")

print(p)  
  
figout=file.path(getwd(),"fig_out","02.food_source","multinom","sourceBar_raw_byfood+pop")
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = p,filename = outpath, device = dev, width = 2.2,height = 7,dpi = 300, scale=1.5)
}


# exclude NAs and not-consumed
select_food <- c("fish","white_meat","red_meat", "honey",
                 "rice","legumes",
                 "tuber","vegetables","fruits",
                 "eggs","noodles","cooking_oil")

src_bar_noNA <- df.src %>%
  pivot_longer(cols=all_of(food.var), names_to = "food", values_to = "source")%>%
  filter(source %in% c("market","garden","wild","isNA")) %>%
  filter(food %in% select_food) %>%
  mutate(population=substr(sampleid,1,3),
         group=case_when(
           population %in% c("BTU","ORT","ORS")  ~ "Early-transition",
           population %in% c("APT","TBU","BSP") ~ "Late-transition",
           population == "LDY" ~ "Agri.",
           TRUE ~ NA_character_),
         island=case_when(
           population %in% c("ORT","ORS") ~ "Sumatra",
           population %in% c("BTU","APT","TBU","BSP","LDY") ~ "Borneo",
           TRUE ~ NA_character_)) 

src_bar_noNA %>% 
  summarise(n=n_distinct(sampleid))

src_bar_noNA %>% 
  group_by(group) %>%
  summarise(n=n_distinct(sampleid))

src_bar_noNA %>% 
  group_by(population) %>%
  summarise(n=n_distinct(sampleid)) %>%
  mutate(xlab=paste0(population, "(",n,")")) -> ncount
ncount

population_labels <- setNames(as.character(ncount$xlab), ncount$population)        

src_bar_noNA <- src_bar_noNA %>%
  group_by(group, island, population, food) %>%
  summarise(n=n_distinct(sampleid),
            market=mean(!is.na(source) & source=="market"),
            garden=mean(!is.na(source) & source=="garden"),
            wild=mean(!is.na(source) & source=="wild"),
            .groups = "drop") %>%
  ungroup() %>%
  pivot_longer(cols=c(wild,market,garden), names_to="source", values_to="proportion") %>%
  mutate(food=factor(food, food_reorder),
         source=factor(source, c("wild","market","garden")),
         population=factor(population, c(pop.ord,"LDY")),
         group=factor(group, c("Early-transition","Late-transition", "Agri.")))

p <- ggplot(src_bar_noNA) +
  theme_bw(base_size=10) +
  geom_bar(aes(x=population, y=proportion, fill=source), stat="identity", width = 1)+
  facet_grid(food~group, scales="free_x", space = "free_x")+
  scale_fill_manual(values=c(market="gold2",
                             garden="royalblue1",
                             wild="springgreen3",
                             not_consumed="grey50",
                             isNA="grey80"))+
  scale_x_discrete(expand = c(0,0), labels = population_labels) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.5,1), labels = scales::percent)+
  theme(panel.grid = element_blank())+
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(0, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.text.x = element_text(face="bold"),
        strip.text.y = element_text(angle=-90, hjust=0.5, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1,size=7))+
  theme(legend.position = "top")+
  geom_hline(yintercept = 0.5, linetype="dashed", col="black")+
  labs(title="Proportion of food sources", x="Population", y="Proportion")

print(p)  

figout=file.path(getwd(),"fig_out","02.food_source","multinom","sourceBar_noNA_byfood+pop")
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = p,filename = outpath, device = dev, width = 2.5,height = 6,dpi = 300, scale=1.5)
}



# save_sample_list
saveRDS(unique(boxdata$sampleid), "input_files/source_sample_list.rds")
saveRDS(sample_list, "input_files/sourcebyfood_sample_list.rds")

write.table(boxdata, file="out/food_source_data_processed.tsv",quote = F,row.names = F,sep = "\t")
write.table(src_bar, file="out/food_source_byFood_data_processed.tsv",quote = F,row.names = F,sep = "\t")




