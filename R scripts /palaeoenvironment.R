################################################################################

library(reshape)
library(ggplot2)
library(readr)
library(dplyr)
library(cowplot)
library(readxl)
library(data.table)
library(extrafont)

#setwd("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment")

### define theme & set x axis maximum
##################################################################################
windowsFonts(Calibri = windowsFont("Calibri"))

proxy_theme <-  theme(text = element_text(family = "Calibri", colour= "#474F58"), panel.grid = element_line(), 
                      panel.background = element_rect(fill = 'white'), panel.spacing = element_blank(),
                      axis.text.x=element_blank(), axis.title.x=element_blank(),
                      axis.ticks.y = element_line(linewidth=0.5), axis.line.y = element_line(linewidth=0.5, colour = "#474F58"),
                      axis.ticks.x = element_blank(), axis.title.y=element_text(size = 14, vjust = 0.5),
                      axis.title.y.right = element_text(size = 14, vjust = 0.5),
                      legend.key=element_blank(), axis.text.y=element_text(size = 9)) 

##################################################################################

### set up alternating background bars: define rects
##################################################################################

rects <- data.frame(ymin = seq(0,10,1), ymax = seq(1,11,1), col = rep(c(1, -1), length.out = 11))
rects$col <- ifelse(rects$col==-1, 'white', 'gray70')

##################################################################################

### timeslices (needed for linear interpolations)

###import metadata
PoolCap01_13_metadata <- fread("C:/Users/vcl943/Documents/Analysis/R/ngsLCA/metadata_files/PoolCap01_13_metadata.csv")
timeslices <- unique(PoolCap01_13_metadata$age_median)

PoolScreen01_metadata <- fread("C:/Users/vcl943/Documents/Analysis/R/ngsLCA/metadata_files/PoolScreen01_metadata.csv")
PoolScreen01_metadata_red <- PoolScreen01_metadata %>% select(!sample_ID)


all_metadata <- rbind(PoolCap01_13_metadata, PoolScreen01_metadata_red) 
timeslices <- unique(all_metadata$age_median)
##################################################################################


### air temp reconstruction based on Lecavalier et al. 2017
##################################################################################

agassiz_temp_data <- read_excel("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/Temperature_Reconstruction.xlsx", sheet = "data_Lecavalier")
#adjust this number if necessary!

#filter temp_data to younger than 11.66 ky BP
agassiz_temp_data <- filter(agassiz_temp_data,`Age [ka BP]` < 11.66)

#plot agassiz temp reconstruction:
agassiz_temp_curve <- ggplot(data=agassiz_temp_data,aes(`Age [ka BP]`, MAAT_preindustrial)) +
  geom_rect(data=rects, aes(ymin=-3, ymax=9, xmin=ymin,xmax=ymax), 
            fill=rects$col, inherit.aes=FALSE, alpha=0.2) +
  geom_ribbon(data=agassiz_temp_data, aes(x=`Age [ka BP]`, ymin=MAAT_preindustrial_2_sigma_lower, 
                                          ymax=MAAT_preindustrial_2_sigma_upper),fill="grey70", alpha=0.4) +
  geom_line(color="#474F58",linewidth=0.8) +
  scale_y_continuous(name ="MAAT anomaly (°C)", breaks=c(-2,0,2,4,6,8)) +
  geom_segment(y= 0,yend=0, color="grey70", x=-x_axis_max_age, xend=0.1) +
  theme(legend.position="none")   + scale_x_reverse(limits = c(x_axis_max_age,-0.1)) +  proxy_theme 

#prepare for linear interpolation
temp.df<- agassiz_temp_data %>% select(c(`Age [ka BP]`,`MAAT [°C]`))
temp.df$`Age [ka BP]`<-temp.df$`Age [ka BP]`*1000
colnames(temp.df)[1] <- "age_bp"
colnames(temp.df)[2] <- "maat_c"
temp.df=as.data.frame(temp.df)

# linear temperature interpolations for each core 

{y <- temp.df$maat_c
 x <- temp.df$age_bp
 y_interpol <- approx(x, y, xout=timeslices, method="linear")$y
 temp_interpol <- data.frame(Age_a=timeslices, temp_interpol=y_interpol)}

# merge temp_interpol with PoolCap01_13_metadata

PoolCap01_13_metadata$temp_interpol=temp_interpol$temp_interpol[match(PoolCap01_13_metadata$age_median, temp_interpol$Age_a)]

# merge temp_interpol with all_metadata
all_metadata$temp_interpol=temp_interpol$temp_interpol[match(all_metadata$age_median, temp_interpol$Age_a)]

##################################################################################


### PALAEOPROXIES NORTH-EAST GREENLAND 

### NEG: foram data (Pados-Dibattista)
DA17_st7_calc_form <- read_delim("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/proxies/Pados-Dibattista-etal_2021/datasets/DA17-NG-ST7-73G_calc_foram_modified.tab", 
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

DA17_st7_calc_form <- DA17_st7_calc_form %>%
                      mutate(atl_forams = `C. neoteretis [#/g]`+ `P. bulloides [#/g]`)

### NEG: dinosterol brassicasterol (Brassicasterol & Dinosterol) PS093-025 (Syring et al 2020)
sterol_data_PS93_025 <- read_delim("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/proxies/SyringN_2019/datasets/PS93-025_Steroles_edit.tab", 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)

### NEG: sterol sea ice reconstructions (IP25 & HBI III) PS93-025 (Syring et al 2020)
ice_data_PS93_025 <- read_delim("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/proxies/SyringN_2019/datasets/PS93-025_HBI_edit.tab", 
                                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)

################################################################################

# ip25 prepare for linear interpolation
NEG_ip25.df<- ice_data_PS93_025 %>% select(c(`Age [ka BP]`,`IP25/TOC [µg/g]`))
NEG_ip25.df$`Age [ka BP]`<-NEG_ip25.df$`Age [ka BP]`*1000
colnames(NEG_ip25.df)[1] <- "age_bp"
colnames(NEG_ip25.df)[2] <- "ip25"
NEG_ip25.df=as.data.frame(NEG_ip25.df)
NEG_ip25.df <- na.omit(NEG_ip25.df)

# brassicasterol prepare for linear interpolation
NEG_brass.df<- sterol_data_PS93_025 %>% select(c(`Age [ka BP]`,`Brassicasterol/TOC [µg/g]`))
NEG_brass.df$`Age [ka BP]`<-NEG_brass.df$`Age [ka BP]`*1000
colnames(NEG_brass.df)[1] <- "age_bp"
colnames(NEG_brass.df)[2] <- "brass"
NEG_brass.df=as.data.frame(NEG_brass.df)
NEG_brass.df <- na.omit(NEG_brass.df)

# dinosterol prepare for linear interpolation
NEG_dino.df<- sterol_data_PS93_025 %>% select(c(`Age [ka BP]`,`Dinosterol/TOC [µg/g]`))
NEG_dino.df$`Age [ka BP]`<-NEG_dino.df$`Age [ka BP]`*1000
colnames(NEG_dino.df)[1] <- "age_bp"
colnames(NEG_dino.df)[2] <- "dino"
NEG_dino.df=as.data.frame(NEG_dino.df)
NEG_dino.df <- na.omit(NEG_dino.df)

# forams prepare for linear interpolation
NEG_forams.df<- DA17_st7_calc_form %>% select(c(`Age [ka BP]`,`atl_forams`))
NEG_forams.df$`Age [ka BP]`<-NEG_forams.df$`Age [ka BP]`*1000
colnames(NEG_forams.df)[1] <- "age_bp"
NEG_forams.df=as.data.frame(NEG_forams.df)
NEG_forams.df <- na.omit(NEG_forams.df)

### perform interpolation for NEG

{ y <- NEG_ip25.df$ip25
  x <- NEG_ip25.df$age_bp
  y_ip25_interpol <- approx(x, y, xout=timeslices, method="linear")$y
  y <- NEG_brass.df$brass
  x <- NEG_brass.df$age_bp
  y_brass_interpol <- approx(x, y, xout=timeslices, method="linear")$y
  y <- NEG_dino.df$dino
  x <- NEG_dino.df$age_bp
  y_dino_interpol <- approx(x, y, xout=timeslices, method="linear")$y
  y <- NEG_forams.df$atl_forams
  x <- NEG_forams.df$age_bp
  y_forams_interpol <- approx(x, y, xout=timeslices, method="linear")$y
  NEG_proxy_interpol <- data.frame(Age_a=timeslices, 
                                   ip25_interpol=y_ip25_interpol, 
                                   brass_interpol= y_brass_interpol,
                                   dino_interpol = y_dino_interpol,
                                   forams_interpol = y_forams_interpol)}


# merge NEG_proxy_interpol with PoolCap01_13_metadata

PoolCap01_13_metadata$ip25_interpol <- ifelse(PoolCap01_13_metadata$core == "DA17_st7_73G", NEG_proxy_interpol$ip25_interpol[match(PoolCap01_13_metadata$age_median, NEG_proxy_interpol$Age_a)], NA)
PoolCap01_13_metadata$brass_interpol <- ifelse(PoolCap01_13_metadata$core == "DA17_st7_73G", NEG_proxy_interpol$brass_interpol[match(PoolCap01_13_metadata$age_median, NEG_proxy_interpol$Age_a)], NA)
PoolCap01_13_metadata$dino_interpol <- ifelse(PoolCap01_13_metadata$core == "DA17_st7_73G", NEG_proxy_interpol$dino_interpol[match(PoolCap01_13_metadata$age_median, NEG_proxy_interpol$Age_a)], NA)
PoolCap01_13_metadata$forams_interpol <- ifelse(PoolCap01_13_metadata$core == "DA17_st7_73G", NEG_proxy_interpol$forams_interpol[match(PoolCap01_13_metadata$age_median, NEG_proxy_interpol$Age_a)], NA)

# merge NEG_proxy_interpol with all_metadata
all_metadata$ip25_interpol <- ifelse(all_metadata$core == "DA17_st7_73G", NEG_proxy_interpol$ip25_interpol[match(all_metadata$age_median, NEG_proxy_interpol$Age_a)], NA)
all_metadata$brass_interpol <- ifelse(all_metadata$core == "DA17_st7_73G", NEG_proxy_interpol$brass_interpol[match(all_metadata$age_median, NEG_proxy_interpol$Age_a)], NA)
all_metadata$dino_interpol <- ifelse(all_metadata$core == "DA17_st7_73G", NEG_proxy_interpol$dino_interpol[match(all_metadata$age_median, NEG_proxy_interpol$Age_a)], NA)
all_metadata$forams_interpol <- ifelse(all_metadata$core == "DA17_st7_73G", NEG_proxy_interpol$forams_interpol[match(all_metadata$age_median, NEG_proxy_interpol$Age_a)], NA)

################################################################################


### PALAEOPROXIES Melville Bay 

### Saini et al 2020: GeoB19927
GeoB19927 <- read_delim("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/proxies/Saini-etal_2020/datasets/Saini-etal_2020_edit.tab", 
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#filter dataframe for NA rows
GeoB19927_noNA <- filter(GeoB19927, !is.na(`IP25/TOC [µg/g]`))

### LK21st26 foram counts:
st26_proxies <- read_excel("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/proxies/LK21-IC-st26-GC1_algae_counts.xlsx", sheet = "Foraminifera #")
st26_proxies <- subset(st26_proxies, select = c("median_age", "chilled_atl_forams/gram", "polar_arct_forams/gram"))
st26_proxies <- filter(st26_proxies, median_age >0)
st26_proxies <- mutate(st26_proxies, median_age = round(median_age/1000,1))

# ip25 prepare for linear interpolation
MB_ip25.df<- GeoB19927_noNA %>% select(c(`Age [ka BP]`,`IP25/TOC [µg/g]`))
MB_ip25.df$`Age [ka BP]`<-MB_ip25.df$`Age [ka BP]`*1000
colnames(MB_ip25.df)[1] <- "age_bp"
colnames(MB_ip25.df)[2] <- "ip25"
MB_ip25.df=as.data.frame(MB_ip25.df)
MB_ip25.df <- na.omit(MB_ip25.df)

# brassicasterol prepare for linear interpolation
MB_brass.df<- GeoB19927_noNA %>% select(c(`Age [ka BP]`,`Brassicasterol/TOC [µg/g]`))
MB_brass.df$`Age [ka BP]`<-MB_brass.df$`Age [ka BP]`*1000
colnames(MB_brass.df)[1] <- "age_bp"
colnames(MB_brass.df)[2] <- "brass"
MB_brass.df=as.data.frame(MB_brass.df)
MB_brass.df <- na.omit(MB_brass.df)

# dinosterol prepare for linear interpolation
MB_dino.df<- GeoB19927_noNA %>% select(c(`Age [ka BP]`,`Dinosterol/TOC [µg/g]`))
MB_dino.df$`Age [ka BP]`<-MB_dino.df$`Age [ka BP]`*1000
colnames(MB_dino.df)[1] <- "age_bp"
colnames(MB_dino.df)[2] <- "dino"
MB_dino.df=as.data.frame(MB_dino.df)
MB_dino.df <- na.omit(MB_dino.df)

# forams prepare for linear interpolation
MB_forams.df<- st26_proxies %>% select(c(`median_age`,`chilled_atl_forams/gram`))
MB_forams.df$median_age<-MB_forams.df$median_age*1000
colnames(MB_forams.df)[1] <- "age_bp"
colnames(MB_forams.df)[2] <- "chilled_atl_forams"
MB_forams.df=as.data.frame(MB_forams.df)
MB_forams.df <- na.omit(MB_forams.df)

### perform interpolation for MB

{ y <- MB_ip25.df$ip25
  x <- MB_ip25.df$age_bp
  y_ip25_interpol <- approx(x, y, xout=timeslices, method="linear")$y
  y <- MB_brass.df$brass
  x <- MB_brass.df$age_bp
  y_brass_interpol <- approx(x, y, xout=timeslices, method="linear")$y
  y <- MB_dino.df$dino
  x <- MB_dino.df$age_bp
  y_dino_interpol <- approx(x, y, xout=timeslices, method="linear")$y
  y <- MB_forams.df$chilled_atl_forams
  x <- MB_forams.df$age_bp
  y_forams_interpol <- approx(x, y, xout=timeslices, method="linear")$y
  MB_proxy_interpol <- data.frame(Age_a=timeslices, 
                                   ip25_interpol=y_ip25_interpol, 
                                   brass_interpol= y_brass_interpol,
                                   dino_interpol = y_dino_interpol,
                                   forams_interpol = y_forams_interpol)}


# merge MB_proxy_interpol with PoolCap01_13_metadata

PoolCap01_13_metadata$ip25_interpol <- ifelse(PoolCap01_13_metadata$core == "LK21_ICst26", MB_proxy_interpol$ip25_interpol[match(PoolCap01_13_metadata$age_median, MB_proxy_interpol$Age_a)], PoolCap01_13_metadata$ip25_interpol)
PoolCap01_13_metadata$brass_interpol <- ifelse(PoolCap01_13_metadata$core == "LK21_ICst26", MB_proxy_interpol$brass_interpol[match(PoolCap01_13_metadata$age_median, MB_proxy_interpol$Age_a)], PoolCap01_13_metadata$brass_interpol)
PoolCap01_13_metadata$dino_interpol <- ifelse(PoolCap01_13_metadata$core == "LK21_ICst26", MB_proxy_interpol$dino_interpol[match(PoolCap01_13_metadata$age_median, MB_proxy_interpol$Age_a)], PoolCap01_13_metadata$dino_interpol)
PoolCap01_13_metadata$forams_interpol <- ifelse(PoolCap01_13_metadata$core == "LK21_ICst26", MB_proxy_interpol$forams_interpol[match(PoolCap01_13_metadata$age_median, MB_proxy_interpol$Age_a)], PoolCap01_13_metadata$forams_interpol)

# merge MB_proxy_interpol with PoolCap01_13_metadata

all_metadata$ip25_interpol <- ifelse(all_metadata$core == "LK21_ICst26", MB_proxy_interpol$ip25_interpol[match(all_metadata$age_median, MB_proxy_interpol$Age_a)], all_metadata$ip25_interpol)
all_metadata$brass_interpol <- ifelse(all_metadata$core == "LK21_ICst26", MB_proxy_interpol$brass_interpol[match(all_metadata$age_median, MB_proxy_interpol$Age_a)], all_metadata$brass_interpol)
all_metadata$dino_interpol <- ifelse(all_metadata$core == "LK21_ICst26", MB_proxy_interpol$dino_interpol[match(all_metadata$age_median, MB_proxy_interpol$Age_a)], all_metadata$dino_interpol)
all_metadata$forams_interpol <- ifelse(all_metadata$core == "LK21_ICst26", MB_proxy_interpol$forams_interpol[match(all_metadata$age_median, MB_proxy_interpol$Age_a)], all_metadata$forams_interpol)


#################################################################################

### PALAEOPROXIES NORTH GREENLAND

### Detlef et al. 2023

Ry19_12_detlef <- read_xlsx("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/proxies/oden-ryder-2019-sediment-detlef-lincoln-sea-1-2/biomarker-ryder19-12-GC1.xlsx")

### Jennings et al 2011: Ry19-12/Ry19-24
jennings2011_forams <- read_xlsx("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/proxies/jennings2011/jennings2011.xlsx", sheet = "Forams")
jennings2011_isotopes <- read_xlsx("C:/Users/vcl943/Documents/Analysis/R/palaeoenvironment/proxies/jennings2011/jennings2011.xlsx", sheet = "Cneoteretis_isotopes")

#################################################################################


### COMPOSITE PLOTS
### !adjust for each core!
#################################################################################

### ocean_temperature: combined plot
### !adjust for each core!
##################################################################################
ocean_temp <- ggplot() +
  geom_rect(data=rects, aes(ymin=-2, ymax=63, xmin=ymin, xmax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2) + proxy_theme + 
  theme(legend.position="none") + scale_x_reverse(limits = c(x_axis_max_age,-0.1), breaks =seq(12, -1, by = -1)) + 
  scale_y_continuous(name = "Atlantic \n foraminifera \n (ind. g⁻¹ sed)") + #,sec.axis = sec_axis(~./25, name="SST (°C)")) +
  ##St7
  #geom_point(data=DA17_st7_calc_form, aes(`Age [ka BP]`, `C. neoteretis [#/g]`+ `P. bulloides [#/g]`), color="#474F58", size = 1.5) +
  #geom_line(data = filter(DA17_st7_calc_form, is.na(`C. neoteretis [#/g]`+ `P. bulloides [#/g]`)==FALSE), aes(`Age [ka BP]`, `C. neoteretis [#/g]`+ `P. bulloides [#/g]`), color="#474F58", linewidth = 0.8, linetype = "dashed") 
  ##st26
  geom_point(data= st26_proxies, aes(median_age, `chilled_atl_forams/gram`), color="#474F58", size = 1.5) +
  geom_line(data = filter(st26_proxies, is.na(`chilled_atl_forams/gram`)==FALSE), aes(median_age, `chilled_atl_forams/gram`), color="#474F58", linewidth = 0.8, linetype = "dashed") 

nares_strait_Niridea <- ggplot() +
  proxy_theme + scale_x_reverse(limits = c(x_axis_max_age,-0.1)) +  
  geom_rect(data=rects, aes(ymin=-1, ymax=30, xmin=ymin, xmax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2) +
  geom_point(data= jennings2011_forams, aes((`CalAge yrBP`)/1000, `N. iridea`), color="grey70", size = 1.5) +
  geom_line(data = filter(jennings2011_forams, is.na(`N. iridea`)==FALSE), aes((`CalAge yrBP`)/1000, `N. iridea`), color="grey70", linewidth = 0.8, linetype = "dashed") +
  geom_point(data= jennings2011_isotopes, aes((`Cal Age, 3/31/11, Cneo`)/1000, (`C13Mean, C.neoteretis`+ 1)*15), color="#474F58", size = 1.5) +
  geom_line(data = filter(jennings2011_isotopes, is.na(`C13Mean, C.neoteretis`)==FALSE), aes((`Cal Age, 3/31/11, Cneo`)/1000, (`C13Mean, C.neoteretis`+ 1)*15), color="#474F58", linewidth = 0.8, linetype = "dashed") +
  scale_y_continuous(name = "N. iridea (%) \n δ13C C. neoteretis" ,sec.axis = sec_axis(~ (. /15) -1))

##################################################################################


### sea ice: combined plot 
### !adjust for each core!
################################################################################

#DA17st7
sea_ice <-  ggplot() +
  geom_rect(data=rects, aes(ymin=-0.005, ymax=1.5, xmin=ymin, xmax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2) +
  proxy_theme +scale_x_reverse(limits = c(x_axis_max_age,-0.1)) +
  geom_point(data=ice_data_PS93_025, aes(`Age [ka BP]`, `IP25/TOC [µg/g]`), color="#474F58", size = 1.5) +
  geom_line(data = filter(ice_data_PS93_025, is.na(`IP25/TOC [µg/g]`)==FALSE), aes(`Age [ka BP]`, `IP25/TOC [µg/g]`), color="#474F58", linewidth = 0.8, linetype = "dashed") +
  ylab("IP₂₅\n (µg g⁻¹ TOC)")

#####################################

#st26

sea_ice <- 
  ggplot(GeoB19927_noNA, aes(`Age [ka BP]`)) + proxy_theme + scale_x_reverse(limits = c(x_axis_max_age,-0.1)) +
  geom_rect(data=rects, aes(ymin=0, ymax=3.5, xmin=ymin, xmax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2) +
  geom_point(data= GeoB19927_noNA, aes(y=`IP25/TOC [µg/g]`), color="#474F58", size = 1.5) +
  geom_line(data = filter(GeoB19927_noNA, is.na(`IP25/TOC [µg/g]`)==FALSE), aes(y=`IP25/TOC [µg/g]`), color="#474F58", linewidth = 0.8, linetype = "dashed") +
  ylab("IP₂₅\n (µg g⁻¹ TOC)") 
  
################################################################################


### productivity
### !adjust for each core!
################################################################################


productivity <- GeoB19927_noNA %>%
  ggplot(aes(`Age [ka BP]`)) + proxy_theme +
  geom_rect(data=rects, aes(ymin=-1, ymax=38, xmin=ymin, xmax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2) + 
  geom_point(aes(y=`Brassicasterol/TOC [µg/g]`), color="grey70", size = 1.5) +
  geom_line(data = filter(GeoB19927_noNA, is.na(`Brassicasterol/TOC [µg/g]`)==FALSE), aes(y=`Brassicasterol/TOC [µg/g]`), color="grey70", linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(y=`Dinosterol/TOC [µg/g]`), color="#474F58", size = 1.5) +
  geom_line(data = filter(GeoB19927_noNA, is.na(`Dinosterol/TOC [µg/g]`)==FALSE), aes(y=`Dinosterol/TOC [µg/g]`), color="#474F58", linewidth = 0.8, linetype = "dashed") +
  scale_x_reverse(limits = c(x_axis_max_age,-0.1)) +
  scale_y_continuous(name = "Brassicasterol \n (µg g⁻¹ TOC) \n Dinosterol \n (µg g⁻¹ TOC)") 

################################################################################


Ry19_12_comb <- Ry19_12_detlef %>%
  ggplot(aes(`Age_(cal-yrs-BP)_dR=300±0`/1000)) + proxy_theme +
  geom_rect(data=rects, aes(ymin=-10, ymax=435, xmin=ymin, xmax=ymax), fill=rects$col, inherit.aes=FALSE, alpha=0.2) + 
  geom_point(aes(y=`Brassicasterol_(ug/g-TOC)`), color="grey70", size = 1.5) +
  geom_line(data = filter(Ry19_12_detlef, is.na(`Brassicasterol_(ug/g-TOC)`)==FALSE), aes(y=`Brassicasterol_(ug/g-TOC)`), color="grey70", linewidth = 0.8, linetype = "dashed") +
  scale_x_reverse(limits = c(x_axis_max_age,-0.1)) +
  geom_point(aes(y=(`IP25_(ug/g-TOC)`)*667), color="#474F58", size = 1.5) +
  geom_line(data = filter(Ry19_12_detlef, is.na(`IP25_(ug/g-TOC)`)==FALSE), aes(y=(`IP25_(ug/g-TOC)`)*667), color="#474F58", linewidth = 0.8, linetype = "dashed") +
  scale_y_continuous(name = "Brassicasterol  \n (µg g⁻¹ TOC) \n IP₂₅\n (µg g⁻¹ TOC) ", sec.axis = sec_axis(~ (. /667))) 


