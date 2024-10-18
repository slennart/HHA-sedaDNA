
library(tidyr)
library(vegan)
library(psych)
library(ggpubr)
library(ggord)
library(reshape2)

#################################################################################################################################################
### Detection matrix all cores

PoolCap01_13_filt_all_w <- PoolCap01_13_filt_target %>% 
  filter(core %in% c("LK21_ICst26", "Ryder19_12_GC","DA17_st7_73G", "Ryder19_24_PC")) %>% 
  filter(rank == "species") %>% 
  group_by(age_median, taxon_name) %>%
  reframe(read_count_sum = sum(read_count_sum), across()) %>%
  slice(which.max(read_count_sum), .by = c(age_median, taxon_name))  %>%
  select (taxon_name, read_count_sum, age_median) %>% 
  pivot_wider(names_from = "taxon_name", values_from = "read_count_sum", values_fill = 0) %>% 
  data.frame()


### only capture data - all cores
# add zero count samples
hyb_data_all_w <- PoolCap01_13_filt_all_w  %>% add_row(age_median = setdiff(unique(PoolCap01_13_metadata$age_median), PoolCap01_13_filt_all_w $age_median))
hyb_data_all_w[is.na(hyb_data_all_w)] <- 0

hyb_data_all_w_ordered <- hyb_data_all_w[order(hyb_data_all_w$age_median),]

row.names(hyb_data_all_w_ordered) <- hyb_data_all_w_ordered$age_median
hyb_data_all_w_ordered[1] = NULL
hyb_data_all_w_ordered_0rm <- subset(hyb_data_all_w_ordered, rowSums(hyb_data_all_w_ordered)!=0) 

### both capture & shotgun data - all cores
PoolScreen01_filt_all_w <- PoolScreen01_filt %>% 
  filter(core %in% c("LK21_ICst26", "Ryder19_12_GC","DA17_st7_73G", "Ryder19_24_PC")) %>%
  filter(!family_name  == "") %>% #only keep family level and lower taxonomic assignments
  group_by(family_name, age_median) %>%
  reframe(read_count_sum = sum(read_count), across()) %>% 
  slice(which.max(read_count_sum), .by = c(family_name,age_median)) %>%
  ungroup()%>%
  group_by(family_name, core) %>%
  filter(n()>2) %>%
  ungroup() %>%
  select (family_name, read_count_sum, age_median) %>% 
  pivot_wider(names_from = "family_name", values_from = "read_count_sum", values_fill = 0) %>% 
  data.frame()

sg_hyb_data_all_w <- merge(PoolCap01_13_filt_all_w, PoolScreen01_filt_all_w, by="age_median")

# add zero count samples
sg_hyb_data_all_w <- sg_hyb_data_all_w %>% add_row(age_median = setdiff(unique(PoolCap01_13_metadata$age_median), sg_hyb_data_all_w$age_median))
sg_hyb_data_all_w[is.na(sg_hyb_data_all_w)] <- 0

sg_hyb_data_all_w_ordered <- sg_hyb_data_all_w[order(sg_hyb_data_all_w$age_median),]

row.names(sg_hyb_data_all_w_ordered) <- sg_hyb_data_all_w_ordered$age_median
sg_hyb_data_all_w_ordered[1] = NULL
sg_hyb_data_all_w_ordered_0rm <- subset(sg_hyb_data_all_w_ordered, rowSums(sg_hyb_data_all_w_ordered)!=0) 

################################################################################################################################################
################################################################################################################################################

#Correlation analysis

#pairwise correlations including interpolated proxy data with Benjamini-Hochberg p-value correction

# first filter for species occurring in more than one sample
res <- sg_hyb_data_all_w_ordered_0rm[, apply(sg_hyb_data_all_w_ordered_0rm, 2, function(x) sum(x > 0)) > 1]
res <- subset(res, rowSums(res)!=0) 


interpolations <- all_metadata %>%
                    filter(core %in% c("LK21_ICst26", "DA17_st7_73G","LK21_ICst26", "Ryder19_24_PC")) %>%
                    select(age_median, temp_interpol, ip25_interpol, brass_interpol, dino_interpol, forams_interpol) %>%
                    data.frame()

interpolations<- interpolations[order(interpolations$age_median),]
interpolations<- interpolations %>% filter(!age_median == -8) %>% distinct(age_median, .keep_all = T)
row.names(interpolations) <- interpolations$age_median
interpolations[1] = NULL
interpolations.z <- decostand(interpolations, method = "standardize")
#### use standardized!


res.ct.int = subset(res, (row.names(res) %in% row.names(interpolations.z)))
interpolations.z = subset(interpolations.z, (row.names(interpolations.z) %in% row.names(res)))
res.ct.int.ma=cbind(interpolations.z, res.ct.int)
res.ct.int.ma=as.matrix(res.ct.int.ma)

corr.ma=corr.test(res.ct.int.ma, y = NULL, use = "pairwise",method="spearman",
                       adjust="BH",
                       alpha=.1,ci=TRUE,minlength=5)

#visualize

cormat <- corr.ma$r
cormat.p <- corr.ma$p

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)
melted_cormat <- mutate(melted_cormat, p = melted_cormat.p$value)
melted_cormat <- filter(melted_cormat, p < 0.1)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)


corr.ma.coeff.ma=corr.ma$r
corr.ma.p=corr.ma$p

corr.ma.coeff.ma[ lower.tri(corr.ma.coeff.ma, diag=TRUE) ]<- 0
corr.ma.p[ lower.tri(corr.ma.p, diag=TRUE) ]<- 1

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#filter for correlations above the threshold of the coefficient 
#rho>0.4 and p<0.1
flat.matrix.ma=flattenCorrMatrix(corr.ma.coeff.ma, corr.ma.p)
flat.matrix.ma.s<-subset(flat.matrix.ma, p < 0.1)
flat.matrix.ma.s.cor<-subset(flat.matrix.ma.s, abs(cor) > 0.6)

################################################################################################################################################
################################################################################################################################################

### ORDINATION ANALYSIS

###         RDA 

ages <-unique((filter(PoolCap01_13_metadata, core %in% c("LK21_ICst26", "DA17_st7_73G")))$age_median)
###  ADJUST CORE NAMES!
input_data <- hyb_data_all_w_ordered_0rm[rownames(hyb_data_all_w_ordered_0rm) %in% ages,]
#input_data <- input_data[-1,]

vare.dca <- decorana(input_data)
vare.dca
# -> length of DCA1 <4 means homogeneous dataset where rda is more appropriate than cca 

data_hel <- decostand(input_data, method="hellinger")

interpolations <- all_metadata %>%
  filter(core %in% c("LK21_ICst26", "DA17_st7_73G")) %>% # "Ryder19_12_GC" , "Ryder19_24_PC", "DA17_st7_73G"
  select(age_median, temp_interpol, ip25_interpol, brass_interpol, dino_interpol, forams_interpol) %>%
  data.frame()

interpolations <- filter(interpolations, !is.na(temp_interpol))
interpolations<- interpolations[order(interpolations$age_median),]
interpolations<- interpolations %>% distinct(age_median, .keep_all = T)
row.names(interpolations) <- interpolations$age_median
interpolations = subset(interpolations, (row.names(interpolations) %in% row.names(data_hel)))
data_hel = subset(data_hel, (row.names(data_hel) %in% row.names(interpolations)))
interpolations[1] = NULL
interpolations.z <- decostand(interpolations, method = "standardize")

vare.rda <- rda(data_hel~temp_interpol+ ip25_interpol + brass_interpol+ dino_interpol+ forams_interpol, data = interpolations.z, na.action = na.exclude) 
head(summary(vare.rda))
plot(vare.rda)
anova(vare.rda, step = 1000, by = "terms")

## extract % explained by the first 2 axes
perc <- round(100*(summary(vare.rda)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(vare.rda, display="sites", choices=c(1,2), scaling=2)
sc_sp <- scores(vare.rda, display="species", choices=c(1,2), scaling=2)
sc_bp <- scores(vare.rda, display="bp", choices=c(1, 2), scaling=2)


### CODE FROM https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html 

###  ADJUST CORE NAMES!
pdf("raw_figures/ordination_analysis/73G_26G_hyb_all_rda_v2.pdf", width=8, height=6)
# Set up a blank plot with scaling, axes, and labels
plot(vare.rda,
     scaling = 2, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-1,1), 
     ylim = c(-1,1),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 2",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1.2) # size
# add points for species scores
points(sc_sp, 
       pch = 22, # set shape (here, square with a fill colour)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_sp), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)
# add arrows for effects of the explanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "grey2", 
       lwd = 3)
# add text labels for arrows
text(x = sc_bp[,1] + 0.2, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.1, 
     labels = rownames(sc_bp), 
     col = "grey2", 
     cex = 1, 
     font = 2)

dev.off()
