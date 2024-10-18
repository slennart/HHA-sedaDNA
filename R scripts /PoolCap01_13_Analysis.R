

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readxl)
library(extrafont)

################################################################################
################################################################################

###import metadata
PoolCap01_13_metadata <- fread("metadata_files/PoolCap01_13_metadata.csv")

###import read count per sample data
PoolCap01_13_data <- fread("PoolCap01_13_ss98_03Sep24.tsv")

###add metadata to flat file
PoolCap01_13 <- merge(PoolCap01_13_data, PoolCap01_13_metadata, by="lib_id")

################################################################################
################################################################################

#read in damage data

damage_files <- grep(".tsv", list.files("damage_files/"), value = T)
filtered_damage_all <- data.frame()
for(file in damage_files ){
  metaDMG <- fread(paste("damage_files/", file, sep = ""))
  metaDMG[End=="5p", deamination_5p := `C>T`/C]
  #metaDMG[End=="3p", deamination_3p := `C>T`/C]
  metaDMG[,count := mean(Total), by=Chr]
  ## double stranded SE ##
  metaDMG[End=="3p", deamination_3p := `G>A`/G]
  metaDMG_2A<- rbind(metaDMG[Pos<=3 & End=="5p", .(mean=mean(deamination_5p)), by=c("Chr", 'End')],
                     metaDMG[Pos<=3 & End=="3p", .(mean=mean(deamination_3p)), by=c("Chr", 'End')])
  metaDMG_2A[,End:=paste(End, "_mean", sep = "")]
  metaDMG_2B<- rbind(metaDMG[Pos==1 & End=="5p", .(mean=mean(deamination_5p)), by=c("Chr", 'End')],
                     metaDMG[Pos==1 & End=="3p", .(mean=mean(deamination_3p)), by=c("Chr", 'End')])
  metaDMG_2B[,End:=paste(End, "_end", sep = "")]
  
  metaDMG_2 <- merge(data.table(dcast(metaDMG_2A, Chr~End)), 
                     data.table(dcast(metaDMG_2B, Chr~End)), by = "Chr")
  
  filtered_damage <- metaDMG_2
  colnames(filtered_damage) <- c("taxa_ID", "mean_damage_3p", "mean_damage_5p",  "end_damage_3p", "end_damage_5p")
  filtered_damage$lib_id <- gsub(".tsv", "", file)
  filtered_damage_all <- rbind(filtered_damage_all, filtered_damage)}

#rename taxa_ID to taxid
colnames(filtered_damage_all)[1] <- "taxid"
###add damage data to flat file
PoolCap01_13 <- merge(PoolCap01_13, filtered_damage_all, by=c("lib_id", "taxid"))

################################################################################
################################################################################

#minimum read filter
PoolCap01_13_filt <- PoolCap01_13[read_count_sum>=3]

### make datatables

#### summarize taxa, read_count_sum per library
PoolCap01_13_filt_w <- PoolCap01_13_filt %>% 
  filter(core %in% c("LK21_ICst26", "Ryder19_12_GC","blank", "DA17_st7_73G", "Ryder19_24_PC")) %>%
  filter(rank == "species") %>%
  select (lib_id, core, depth_in_core, age_median, taxon_name, read_count_sum) %>%
  group_by(taxon_name, lib_id) %>%
  reframe(read_count_sum = sum(read_count_sum), across()) %>%
  slice(which.max(read_count_sum), .by = c(taxon_name,lib_id)) %>%
  pivot_wider(names_from = "taxon_name", values_from = "read_count_sum", values_fill = 0) %>% 
  arrange(core, depth_in_core) %>%
  data.frame()

#calculate relative abundance of reads (based on unique read counts & number of clean reads per sample)
PoolCap01_13_filt <- mutate(PoolCap01_13_filt, rel_read_count = read_count_sum/cleaned_reads)

################################################################################
################################################################################
##overall stats

PoolCap01_13_filt_target <- filter(PoolCap01_13_filt, family_name %in% c("Phocidae", "Balaenidae", "Monodontidae",
                                                                         "Balaenopteridae", "Delphinidae"))

#overall on_target read counts (family level)
sum(PoolCap01_13_filt_target$read_count)

#per sample read counts (family level)
agg_tbl <- PoolCap01_13_filt_target %>% group_by(lib_id) %>% 
  summarise(sum_counts=sum(read_count),
            .groups = 'drop')

#mean read count per sample + sd (family level)
mean(agg_tbl$sum_counts)
sd(agg_tbl$sum_counts)

################################################################################
################################################################################

#extract age model information 
sample_ages_all <- PoolCap01_13_metadata %>% distinct(age_median, .keep_all = TRUE)
sample_ages_all$age_median <- as.numeric(sample_ages_all$age_median)

################################################################################

#set up alternating background bars;
rects <- data.frame(ymin = seq(0,11,1), ymax = seq(1,12,1), col = rep(c(1, -1), length.out = 12))
rects$col <- ifelse(rects$col==-1, NA, 'gray70')

geom_rects <-   geom_rect(data=rects, aes(xmin=-Inf, xmax=Inf, ymin=ymin,
                                            ymax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2)
#import font:
windowsFonts(Calibri=windowsFont("Calibri"))

################################################################################
################################################################################

cores <- unique(PoolCap01_13_filt$core)

#set x-axis limit (adjust for each core?):
x_axis_max_age <- 13

#select the core for which plots should be generated, e.g. only LK21_ICst26 
i <- 1

#set up core-specific colour palette: 
colour_palette <- c("#b31629ff", "#d85f4cff", "#f7a482ff", "#fedbc7ff") #LK21st26

################################################################################

#only detections for one core, e.g. ICst26 
data_filtered <- filter(PoolCap01_13_filt, core == cores[i])

#filter metadata table for DNA concentrations and prepare for plotting
DNA_conc_data <- filter(PoolCap01_13_metadata, core == cores[i])
DNA_conc_data <- filter(DNA_conc_data, !is.na(DNA_conc_raw_extracts))
DNA_conc_data <- distinct(DNA_conc_data, age_median, .keep_all = T)


#change class of age_median column from character to numeric
data_filtered$age_median <- as.numeric(data_filtered$age_median)

#exclude family level IDs
data_filtered <- filter(data_filtered, rank =="species")

#filter for on-target reads
data_filtered <- filter(data_filtered, order_name %in% c("Carnivora", "Artiodactyla"))

# create a new binned read count variable from read_count
data_filtered <- data_filtered %>% 
  mutate(read_count_binned=cut(read_count_sum, breaks=c(2, 10, 50, 100, max(read_count_sum, na.rm=T)),
                                                labels=c("3-10", "10-50", "50-100", ">100"))) %>%
  # change level order
  mutate(read_count_binned=factor(as.character(read_count_binned), levels=rev(levels(read_count_binned)))) 

################################################################################
#ages
sample_ages <- filter(sample_ages_all, core == cores[i])

#convert to ky
sample_ages <- mutate(sample_ages, age_min = (age_min/1000), 
                             age_max= (age_max/1000), age_median=(age_median/1000))
#sort by age_median
sample_ages <- sample_ages[order(sample_ages$age_median),]

#sort by age_median
sample_ages <- filter(sample_ages, !is.na(cleaned_reads))

#prepare x-axis labels (position: median, top:min, bottom:max)
x_axis <- sample_ages$age_median

################################################################################
################################################################################

##explore damage
age_damage_3p <- select(data_filtered, taxon_name, read_count_sum, age_median, mean_damage_3p, read_count_binned)
colnames(age_damage_3p)[4] <- "Mean damage"
age_damage_5p <- select(data_filtered, taxon_name, read_count_sum, age_median, mean_damage_5p, read_count_binned)
colnames(age_damage_5p)[4] <- "Mean damage"
age_damage <- new <- rbind(age_damage_3p, age_damage_5p)

age_damage_all <- rbind(age_damage, age_damage_sg)
  
p4 <- ggplot() +
  geom_point(data=age_damage_all, aes(x=age_median, y=`Mean damage`, colour = read_count_binned), size = 2) + 
  theme_classic() + ylim(-0.01, 0.15) + ylab("Mean deamination across the last 3bp") + 
  xlab("estimated median age (cal ka BP)") + ggtitle("North-East Greenland 73G") +
  scale_color_manual(values=c(colour_palette[1], colour_palette[2],  colour_palette[3], colour_palette[4]))

pp1 <- list(p1,p2,p3,p4)
composite_plot <- plot_grid(plotlist=pp1, ncol=2, align='v', axis = "lr", rel_heights = c(3,3,3,3))
composite_plot

################################################################################
#set order of x axis ticks (remember to use inverse order)
data_filtered$taxon_name <- factor(data_filtered$taxon_name, 
                                       levels = c("Pusa hispida", "Cystophora cristata","Halichoerus grypus",
                                                  "Phoca groenlandica","Erignathus barbatus","Monodon monoceros",
                                                  "Delphinapterus leucas", "Orcinus orca",  "Balaena mysticetus",
                                                  "Balaenoptera acutorostrata","Balaenoptera physalus"))
#replace scientific names by common names
names <- scale_x_discrete(labels=c("Pusa hispida" = "Ringed seal","Phoca vitulina" = "Harbor seal",
                          "Cystophora cristata" ="Hooded seal", "Halichoerus grypus" = "Grey seal",
                          "Phoca groenlandica" = "Harp seal", "Erignathus barbatus" = "Bearded seal",
                          "Monodon monoceros" = "Narwhal","Delphinapterus leucas" = "Beluga",
                          "Orcinus orca" = "Orca", "Balaenoptera acutorostrata" = "Minke whale",
                          "Balaena mysticetus" = "Bowhead whale", "Balaenoptera physalus" = "Fin whale")) 


################################################################################
### "length.out" HAS TO BE ADJUSTED FOR EACH CORE
### "rep()" HAS TO BE ADJUSTED FOR EACH CORE

timeline <-ggplot(sample_ages, aes(x=age_min, xend=age_max, 
                                  y=rep(c(2, 1, 0, -1 ,-2, -3), length.out = 42), 
                                  yend=rep(c(2, 1, 0, -1 ,-2, -3), length.out = 42))) +
  geom_rect(data=rects, aes(ymin=-Inf, ymax=Inf, xmin=ymin,
                          xmax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2) + 
  geom_segment(linewidth=0.6, colour= "grey40") +  scale_y_reverse(limits = c(2.5,-3.5)) + 
  scale_x_reverse(limits = c(x_axis_max_age,-0.1), n.breaks= 12) + 
  theme(text = element_text(family = "Calibri", colour= "#1A242F"),
        panel.grid = element_line(), panel.background = element_rect(fill = 'white'), 
        axis.text.y=element_blank(), axis.ticks.y = element_blank(), 
        axis.ticks.x = element_line(linewidth=0.5, colour = "#474F58"), 
        axis.line.x = element_line(linewidth=0.5, colour = "#474F58"),
        axis.title.y = element_blank(), axis.text.x=element_text(size=14, vjust=-0.05), 
        axis.title.x=element_text(size=16, vjust=-0.1)) + 
  xlab("Calibrated thousands of years BP") +   
  geom_point(aes(x = x_axis_2,y= rep(c(1.5, 0.5,-0.5, -1.5, -2.5, -3.5), length.out = 42)), 
             size =1,shape=25, colour = "grey40", fill= "grey40") 

################################################################################
##### DNA concentrations
DNA_conc <-ggplot(DNA_conc_data, aes(x=(age_median/1000), y=log10(DNA_conc_raw_extracts+1.1))) +
  geom_rect(data=rects, aes(ymin=-0.2, ymax=1.4, xmin=ymin,xmax=ymax), 
            fill=rects$col, inherit.aes=FALSE, alpha=0.2) +
  geom_bar(stat="identity", fill=colour_palette[2], alpha=1, width = 0.08) + 
  scale_x_reverse(limits = c(x_axis_max_age,-0.1), labels=NULL) +
  theme(text = element_text(family = "Calibri", colour= "#474F58"), 
        panel.background = element_blank(), panel.spacing = element_blank(),
        axis.text.y=element_text(size=9), axis.ticks.y = element_line(linewidth=0.5), 
        axis.title.y = element_text(size= 14, vjust= 4), 
        axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
        axis.title.x=element_blank(), axis.line.y = element_line(linewidth=0.5, colour = "#474F58")) +
  geom_point(aes(x = (DNA_conc_data$age_median/1000),y= -0.1), size =1,shape=25, colour = "grey40", fill= "grey40") +
  ylab("DNA conc. \n log10(ng uL⁻¹)")

################################################################################
#define theme for heatmap
theme <-   theme(text = element_text(family = "Calibri", colour= "#1A242F"),
                 panel.grid = element_line(), panel.background = element_rect(fill = 'white'),
                 axis.text.y=element_text(size=14), #colour= y_label_colours),
                 axis.text.x=element_text(size=14),axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(), axis.title.x = element_blank(),
                 legend.key=element_blank(), legend.text = element_text(size=14),
                 legend.title = element_text(size=14), legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-1,-1,-1,-1), plot.margin = unit(c(0,0,0,0), "cm"))


################################################################################
##### detection bubble heat map

bubble_heatmap <- ggplot() + geom_rects +
  coord_flip(clip = "off") + scale_y_reverse(breaks =x_axis, labels = NULL, limits = c(x_axis_max_age,-0.1)) +
  geom_point(data= sample_ages, aes(y=age_median, x = 1), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 2), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 3), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 4), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 5), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 6), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 7), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 8), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 9), size = 0.5, colour = "grey70") +
  geom_point(data= sample_ages, aes(y=age_median, x = 10), size = 0.5, colour = "grey70") +
  geom_point(data = data_filtered, aes(x=taxon_name, y=(age_median/1000), size = rel_read_count, colour=read_count_binned), 
             alpha=1, stroke = 0, shape=16) + 
  scale_size_continuous(range = c(4, 12), breaks = c(min(data_filtered$rel_read_count), mean(data_filtered$rel_read_count), 
                                                     max(data_filtered$rel_read_count)), name ='Relative sequence \n abundance') + 
  labs(colour ='Unique sequences') + names + theme +  guides(colour = guide_legend(override.aes = list(size=10))) +
  scale_color_manual(values=c(colour_palette[1],colour_palette[2], colour_palette[3], colour_palette[4])) 


################################################################################

pp1 <- list(sea_ice, productivity, ocean_temp,  
  DNA_conc, bubble_heatmap, timeline)
set_null_device(cairo_pdf)
composite_heatmap <- plot_grid(plotlist=pp1, ncol=1, align='v', axis = "lr", 
                               rel_heights = c(1.8, 1.8, 1.8, 
                                               1.7, 4, 1.8))

composite_heatmap


