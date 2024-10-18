

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readxl)
library(extrafont)

################################################################################
################################################################################

#set up alternating background bars;
rects <- data.frame(ymin = seq(0,11,1), ymax = seq(1,12,1), col = rep(c(1, -1), length.out = 12))
rects$col <- ifelse(rects$col==-1, NA, 'gray70')

geom_rects <-   geom_rect(data=rects, aes(xmin=-Inf, xmax=Inf, ymin=ymin,
                                          ymax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2)

windowsFonts(Calibri = windowsFont("Calibri"))

sg_theme <-   theme(text = element_text(family = "Calibri", colour= "#1A242F"), panel.grid = element_line(), panel.background = element_rect(fill = 'white'),
                    axis.text.y=element_text(size=14), #colour= y_label_colours),
                    axis.line.x = element_blank(),
                    axis.text.x=element_text(size=14),axis.ticks.x = element_blank(),
                    axis.title.y = element_blank(), axis.title.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.key=element_blank(), legend.text = element_text(size=14),
                    legend.title = element_text(size=14), legend.margin=margin(0,0,0,0),
                    plot.margin = unit(c(0,0,0,0), "cm"))

################################################################################
################################################################################

###import metadata
PoolScreen01_metadata <- fread("Analysis/R/ngsLCA/metadata_files/PoolScreen01_metadata.csv")

################################################################################
################################################################################

###import read count per sample data
PoolScreen01_data <- fread("PoolScreen01_03Sep24.tsv")
###correct lib_id
PoolScreen01_data$lib_id <- substr(PoolScreen01_data$lib_id, 0, 14)

###add metadata to flat file
PoolScreen01 <- merge(PoolScreen01_data, PoolScreen01_metadata, by="lib_id")

################################################################################
################################################################################

#extract age model information 
sg_sample_ages_all <- PoolScreen01 %>% distinct(age_median, .keep_all = TRUE)
sg_sample_ages_all$age_median <- as.numeric(sg_sample_ages_all$age_median)

################################################################################
################################################################################
#read in damage data

damage_files <- grep(".tsv", list.files("damage_files/"), value = T)
filtered_damage_all <- data.frame()
for(file in damage_files ){
  metaDMG <- fread(paste("damage_files/", file, sep = ""))
  metaDMG[End=="5p", deamination_5p := `C>T`/C]
  #metaDMG[End=="3p", deamination_3p := `C>T`/C]
  ####CURRENTLY ONLY CONSIDERING 3P DAMAGE OF DOUBLE-STRANDED LIBRARIES
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

#remove _mapDamage suffix
filtered_damage_all$lib_id <- gsub("_EKDL220009290-1A_HWYJCDSX3_L3", "", filtered_damage_all$lib_id)
#rename taxa_ID to taxid
colnames(filtered_damage_all)[1] <- "taxid"

###add damage data to flat file
PoolScreen01 <- merge(PoolScreen01, filtered_damage_all, by=c("lib_id", "taxid"))

################################################################################
################################################################################
##data filtering

#minimum read filter
PoolScreen01_filt <- PoolScreen01[read_count>=10]
#sum(PoolScreen01_filt$read_count) #303543

PoolScreen01_filt_w_bac_sum <- PoolScreen01_filt %>%  filter(core %in% c("LK21_ICst26", "Ryder19_12_GC","blank", "DA17_st7_73G", "Ryder19_24_PC")) %>% dplyr::group_by(lib_id) %>% dplyr::summarize(read_count_lib = sum(read_count))

#remove Bacteria
PoolScreen01_filt <- filter(PoolScreen01_filt, !phylum_name %in% c("", "Pseudomonadota", "Actinomycetota", "Thermodesulfobacteriota", 
                                                                   "Campylobacterota", "Bacteroidota", "Bacillota", "Spirochaetota",
                                                                   "Nitrososphaerota", "Campylobacterota", "Chloroflexota", "Nitrospinota",
                                                                   "Nitrospirota", "Ignavibacteriota", "Planctomycetota", "Verrucomicrobiota",
                                                                   "Myxococcota", "Gemmatimonadota", "Acidobacteriota", "Atribacterota",
                                                                   "Euryarchaeota",  "Uroviricota", "Rhodothermota", "Cyanobacteriota", "Deinococcota"))

PoolScreen01_filt <- filter(PoolScreen01_filt, !grepl('Candidatus', phylum_name))
#sum(PoolScreen01_filt$read_count) # 117588
PoolScreen01_filt_no_bac_sum <- PoolScreen01_filt %>%  filter(core %in% c("LK21_ICst26", "Ryder19_12_GC","blank", "DA17_st7_73G", "Ryder19_24_PC")) %>% dplyr::group_by(lib_id) %>% dplyr::summarize(read_count_lib = sum(read_count))

#stats
PoolScreen01_filt_perc_bac <- merge(PoolScreen01_filt_w_bac_sum, PoolScreen01_filt_no_bac_sum, by="lib_id", all=T)
PoolScreen01_filt_perc_bac <- mutate(PoolScreen01_filt_perc_bac, perc = read_count_lib.y/read_count_lib.x)
PoolScreen01_filt_perc_bac[is.na(PoolScreen01_filt_perc_bac)] <- 0
PoolScreen01_filt_perc_bac <- mutate(PoolScreen01_filt_perc_bac, perc = 1-perc)
mean(PoolScreen01_filt_perc_bac$read_count_lib.x)
sd(PoolScreen01_filt_perc_bac$read_count_lib.x)

min(PoolScreen01_filt_perc_bac$perc)
max(PoolScreen01_filt_perc_bac$perc)
mean(PoolScreen01_filt_perc_bac$perc)

################################################################################
### make datatables

#### summarize taxa, read_count_sum per library
PoolScreen01_filt_w <- PoolScreen01_filt %>% 
  filter(core %in% c("LK21_ICst26", "Ryder19_12_GC","blank", "DA17_st7_73G", "Ryder19_24_PC")) %>%
  select (lib_id, core, depth_in_core, age_median, family_name, read_count_sum) %>%
  filter(!family_name == "") %>%
  group_by(family_name, lib_id) %>%
  reframe(read_count_sum = sum(read_count_sum), across()) %>%
  slice(which.max(read_count_sum), .by = c(family_name,lib_id)) %>%
  pivot_wider(names_from = "family_name", values_from = "read_count_sum", values_fill = 0) %>% 
  data.frame()

################################################################################
#calculate relative abundance of reads (based on unique read counts & number of clean reads per sample)
PoolScreen01_filt <- mutate(PoolScreen01_filt, rel_read_count = read_count/cleaned_reads)

################################################################################
################################################################################

cores <- unique(PoolScreen01_filt$core)

#select the core for which plots should be genereated, e.g. only LK21_ICst26 
i <- c(3)
x_axis_max_age <- 12

#filter by core 
sg_data_filtered <- filter(PoolScreen01_filt, core %in% cores[i])

#filter out samples without age estimate 
sg_data_filtered <- sg_data_filtered %>% filter(!is.na(age_median))

#change class of age_median column from character to numeric
sg_data_filtered$age_median <- as.numeric(sg_data_filtered$age_median)

### only  keep family level and lower tax classfications 
sg_data_filtered <- filter(sg_data_filtered, !family_name =="")
sg_data_filtered <- sg_data_filtered %>% filter(!is.na(family_name))

### merge by family & lib_id
sg_data_filtered <- sg_data_filtered %>%
  slice(which.max(read_count_sum), .by = c(lib_id, family_name)) 

# create a new binned read count variable from read_count
sg_data_filtered <- sg_data_filtered %>% 
  mutate(read_count_binned=cut(read_count_sum, breaks=c(9, 50, 100, max(read_count_sum, na.rm=T)),
                               labels=c("10-50", "50-100", ">100"))) %>%
  # change level order
  mutate(read_count_binned=factor(as.character(read_count_binned), levels=rev(levels(read_count_binned)))) 

################################################################################

#### filter by at least 3 occurrences:
sg_data_filtered <- sg_data_filtered %>%
  group_by(family_name) %>%
  filter(n()>2)

################################################################################
#ages
sg_sample_ages <- filter(sg_sample_ages_all, core %in% cores[i])

#convert to ky
sg_sample_ages <- mutate(sg_sample_ages, age_min = (age_min/1000), 
                             age_max= (age_max/1000), age_median=(age_median/1000))
#sort by age_median
sg_sample_ages_2 <- sg_sample_ages[order(sg_sample_ages$age_median),]

#prepare x-axis labels (position: median, top:min, bottom:max)
sg_x_axis_2 <- sg_sample_ages_2$age_median


################################################################################
#set order of x axis ticks (remember to use inverse order)
sg_data_filtered$family_name <- factor(sg_data_filtered$family_name, 
                                     levels = c("Balaenopteridae", "Delphinidae", "Phocidae", "Canidae", "Hominidae",
                                                "Gadidae","Liparidae","Priapulidae","Chaetocerotaceae", "Cytherideidae", "Bacillariaceae", 
                                                "Pinaceae","Poaceae","Mustelidae","Zoarcidae",
                                                "Chrysopetalidae", "Hippopodiidae","Forskaliidae","Thalassiosiraceae"))

#dirty fix to make the layering work (otherwise it will complain about discrete/continuous scales)
sg_names <- scale_x_discrete(labels=c("Balaenopteridae"="Balaenopteridae", "Delphinidae" = "Delphinidae", "Phocidae"="Phocidae",
                                   "Gadidae"="Gadidae","Chaetocerotaceae"="Chaetocerotaceae", "Canidae" = "Canidae",
                                   "Hominidae" ="Hominidae","Cytherideidae"="Cytherideidae", "Bacillariaceae"="Bacillariaceae",
                                   "Pinaceae"="Pinaceae","Priapulidae"="Priapulidae","Poaceae"="Poaceae","Mustelidae"="Mustelidae",
                                   "Zoarcidae" ="Zoarcidae", "Chrysopetalidae"= "Chrysopetalidae", "Hippopodiidae"="Hippopodiidae",
                                   "Forskaliidae"="Forskaliidae","Thalassiosiraceae"="Thalassiosiraceae", "Liparidae"="Liparidae"))


################################################################################
################################################################################
### "length.out" HAS TO BE ADJUSTED FOR EACH CORE
### "rep()" HAS TO BE ADJUSTED FOR EACH CORE

timeline_sg_2 <-ggplot(sg_sample_ages_2, aes(x=age_min, xend=age_max, 
                                       y=rep(c(2, 1, 0, -1 ,-2, -3), length.out = 23), 
                                       yend=rep(c(2, 1, 0, -1 ,-2, -3), length.out = 23))) +
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
  geom_point(aes(x = sg_x_axis_2,y= rep(c(1.5, 0.5,-0.5, -1.5, -2.5, -3.5), length.out = 23)), 
             size =1,shape=25, colour = "grey40", fill= "grey40") 


################################################################################

bubble_heatmap_shotgun_2 <- ggplot() + geom_rects + 
  coord_flip(clip = "off") + scale_y_reverse(breaks= NULL, limits = c(x_axis_max_age,-0.1)) + 
  geom_point(data= sg_sample_ages, aes(y=age_median, x = 1), size = 0.5, colour = "grey70") +
  geom_point(data= sg_sample_ages, aes(y=age_median, x = 2), size = 0.5, colour = "grey70") +
  geom_point(data= sg_sample_ages, aes(y=age_median, x = 3), size = 0.5, colour = "grey70") +
  geom_point(data= sg_sample_ages, aes(y=age_median, x = 4), size = 0.5, colour = "grey70") +
  geom_point(data= sg_sample_ages, aes(y=age_median, x = 5), size = 0.5, colour = "grey70") +
  geom_point(data = sg_data_filtered, aes(x=family_name, y=(age_median/1000), 
              size = rel_read_count, colour=read_count_binned), alpha=1, stroke = 2, shape=1) +
  scale_size_continuous(range = c(1, 9), breaks = c(min(sg_data_filtered$rel_read_count), 
                                                     mean(sg_data_filtered$rel_read_count), max(sg_data_filtered$rel_read_count))) + 
  labs(colour ='Unique sequences') +  sg_names + sg_theme + 
  labs(size ='Relative sequence \n abundance') + guides(colour = guide_legend(override.aes = list(size=10))) +
  scale_color_manual(values=c(colour_palette[1], colour_palette[2], colour_palette[3])) 

################################################################################

pp1 <- list(sea_ice, productivity, ocean_temp,  
            bubble_heatmap_shotgun_2, timeline_sg_2)
set_null_device(cairo_pdf)
composite_heatmap <- plot_grid(plotlist=pp1, ncol=1, align='v', axis = "lr", 
                               rel_heights = c(1.8, 1.8, 1.8, 
                                               3, 1.8))

composite_heatmap
