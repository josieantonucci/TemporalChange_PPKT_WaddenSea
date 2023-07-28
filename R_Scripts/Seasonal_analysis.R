# Manuscript: Temporal change in phytoplankton diversity and FG composition
#review 1: calculating trends by season - Appendix

#load packages
library(vegan)
library(readr)
library(plyr)
library(lattice)
library(car)
library(ggplot2)
library(grid)
library(gridExtra)
library(calibrate)
library(maps)
library(mapdata)
library(Hmisc)
library(psych)
library(reshape2)
library(agricolae)
library(ggExtra)
library(metafor)
library(cowplot)
library(ggridges)
library(multcomp)
library(RColorBrewer)
library(weights)
library(lubridate)
library(tidyverse)
library(nycflights13)
library(hrbrthemes)
library(scales)
library(ggpubr)

library(ggpmisc)
library(effects)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library("dplyr")

# Set up the Working Directory
setwd("~.../TemporalChange_PPKT_WaddenSea")

setwd("~/Desktop/Job at AWI_HIFMB/nlwkn/my data/Interreg_Josie/InterReg-project/Manuscript/First manuscript_Temporal change/Repository_Git_TeporalChange/TemporalChange_PPKT_WaddenSea/Data")
# Data ----------------------
#download phytoplankton data for the Wadden Sea

PPKT_WS <- read_delim("PPKT_count_WaddenSea_1999_2018.csv", 
                      delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                      locale = locale(decimal_mark = ",", grouping_mark = "."), 
                      trim_ws = TRUE)
#How many stations?
length(unique(PPKT_WS$StationID)) #13 stations
#Station names:
unique(PPKT_WS$StationID)
# 1] "MARSDND"   "HUIBGOT"   "GROOTGND"  "DANTZGT"   "DOOVBWT"   "ROTTMPT3"  "BOOMKDP"   "TERSLG10"  "BOCHTVWTM" "Bork_W_1"  "JaBu_W_1" 
# [12] "Nney_W_2"  "WeMu_W_1" 


# Data Analysis ----------------------

# 1. Standing Diversity ----------------------

## 1.1. Gamma diversity: First calculate the median abundance of each species per station, per year, then calculate diversity -------
# ENSy = for each species, take mean biomass per year and calculate ENS then
# Sy = for each species, take mean biomass per year and calculate S then

### 1.1.1 Calculate the mean biomass per species, per year -----
# PPKT_WS.median.annual<-ddply(PPKT_WS,.(Country, StationID, Year, Species), 
#                              colwise(mean,  .(abundance_l)), na.rm=T)
###Seasonal
names(PPKT_WS)
str(PPKT_WS)

library(ggplot2)
theme_set(theme_classic())
library(hydroTSM)
# Install
#install.packages("wesanderson")
# Load
library(wesanderson)
library(ggthemes)
#reorder variables:
# Load the lubridate package
library(lubridate)
# Convert the "Date" column to a date format
PPKT_WS$Date <- as.Date(PPKT_WS$Date, format = "%d.%m.%Y")
# Extract the month from the "Date" column and create a new "Month" column
PPKT_WS$Month <- month(PPKT_WS$Date)
# convert the "Month" column to a character month name if desired
PPKT_WS$Month_name <- month(PPKT_WS$Date, label = TRUE)
# Define a function to map month to season
get_season <- function(month) {
  if (month %in% c(12, 1, 2)) {
    return("Winter")
  } else if (month %in% c(3, 4, 5)) {
    return("Spring")
  } else if (month %in% c(6, 7, 8)) {
    return("Summer")
  } else {
    return("Autumn")
  }
}

# Convert the Month column to seasons
PPKT_WS$Season <- sapply(PPKT_WS$Month, get_season)


### 1.1.1 Calculate the mean abundance per species, per year -----
PPKT_WS.median.seasonal<-ddply(PPKT_WS,.(Country, StationID, Year, Season, Species),
                               colwise(mean,  .(abundance_l)), na.rm=T)

### 1.1.2 Converting the data Between Wide And Long Forms -----
PPKT_WS.median.seasonal_wide<-dcast(PPKT_WS.median.seasonal, Country + StationID + Season + Year ~ Species, value.var="abundance_l")

### 1.1.3 Transforming all NAs in 0: in all rows -----
PPKT_WS.median.seasonal_wide[,5:443][is.na(PPKT_WS.median.seasonal_wide[,5:443])]<-0

data_1 <- PPKT_WS.median.seasonal_wide

# number of columns
lastspec2<-ncol(data_1)

### 1.1.4 Calculating Richness and ESN per year -----

##S
data_1$median.ENSy<-diversity(data_1[,5:lastspec2], "invsimpson") 
##ENS
data_1$median.Sy<-specnumber(data_1[,5:lastspec2])  

### 1.1.5 Select columns -----
data_1.seasonal<-data_1[,c("Country","StationID","Year","Season","median.Sy","median.ENSy")]

### 1.1.6 Statistical analysis ----- still wrong -------
#### a) General slope ----
#In order to test whether the change in median.Sy over years is significant within each season, you can include an interaction term between "Season" and "Year" in the model. However, to examine the significance of the change within each season, you need to include season-specific fixed effects in the model as well.
#S:
#model formula to test the change in median.Sy over years within each season
# Convert "Year" to a factor and include "Autumn" in the levels
#data_1.seasonal$Year <- factor(data_1.seasonal$Year, levels = unique(data_1.seasonal$Year))
# Fit the linear mixed-effects model with the interaction term and season-specific fixed effects
# MELM.Sy_1 <- lmer(median.Sy ~ Season * Year + (1|StationID) + (0 + Year | Season), data = data_1.seasonal)
# #or
# MELM.Sy_1 <-lmer(median.Sy ~ Year + (1|StationID) + (0 + Year | Season), data = data_1.seasonal)
# summary(MELM.Sy_1)
# tab_model(MELM.Sy_1) 
# fit.MELM.Sy_1 = as.data.frame(Effect(c("Year"),MELM.Sy_1))
# coef(summary(MELM.Sy_1))
# 
# #ENS
# MELM.ENSy_1 <- lmer(median.ENSy ~ Season * Year + (1|StationID) + (0 + Year | Season), data = data_1.seasonal)
# summary(MELM.ENSy_1)
# tab_model(MELM.ENSy_1)
# fit.MELM.ENSy_1 = as.data.frame(Effect(c("Year"),MELM.ENSy_1))
# coef(summary(MELM.ENSy_1))
# 
# tab_model(MELM.Sy_1, MELM.ENSy_1) 


#### b) Separate models per Country -----

### Statistical analysis: ----------
# # Fit the linear mixed-effects model
# model <- lmer(median.bm ~ Year * Season + (1 | StationID) + (Year * Season | Country), data = tot.biom.season)

###1) DE

## a) S
data_1.seasonal_DE <- data_1.seasonal %>% 
  filter(Country=="Germany")
S_seasonal_DE <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_DE) 
summary(S_seasonal_DE)
tab_model(S_seasonal_DE)

#spring:
data_1.seasonal_DE_spring <- data_1.seasonal_DE %>% 
  filter(Season =="Spring")
S_seasonal_DE_spring <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_DE_spring) 
summary(S_seasonal_DE_spring)
tab_model(S_seasonal_DE_spring)
#summer:
data_1.seasonal_DE_summer <- data_1.seasonal_DE %>% 
  filter(Season =="Summer")
S_seasonal_DE_summer<- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_DE_summer) 
summary(S_seasonal_DE_summer)
tab_model(S_seasonal_DE_summer)
#autumn:
data_1.seasonal_DE_aut <- data_1.seasonal_DE %>% 
  filter(Season =="Autumn")
S_seasonal_DE_aut <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_DE_aut) 
summary(S_seasonal_DE_aut)
tab_model(S_seasonal_DE_aut)
#winter:
data_1.seasonal_DE_win <- data_1.seasonal_DE %>% 
  filter(Season =="Winter")
S_seasonal_DE_win <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_DE_win) 
summary(S_seasonal_DE_win)
tab_model(S_seasonal_DE_win)

tab_model(S_seasonal_DE_spring,  S_seasonal_DE_summer,
          S_seasonal_DE_aut, S_seasonal_DE_win, 
          show.ci = FALSE,  digits = 3,
          dv.labels = c("S.Spring-DE","S.Summer-DE", "S.Autumn-DE", "S.Winter-DE"),
          wrap.labels = 6) 

## a) ENS
data_1.seasonal_DE <- data_1.seasonal %>% 
  filter(Country=="Germany")
ENS_seasonal_DE <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_DE) 
summary(ENS_seasonal_DE)
tab_model(ENS_seasonal_DE)

#spring:
data_1.seasonal_DE_spring <- data_1.seasonal_DE %>% 
  filter(Season =="Spring")
ENS_seasonal_DE_spring <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_DE_spring) 
summary(ENS_seasonal_DE_spring)
tab_model(ENS_seasonal_DE_spring)
#summer:
data_1.seasonal_DE_summer <- data_1.seasonal_DE %>% 
  filter(Season =="Summer")
ENS_seasonal_DE_summer<- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_DE_summer) 
summary(ENS_seasonal_DE_summer)
tab_model(ENS_seasonal_DE_summer)
#autumn:
data_1.seasonal_DE_aut <- data_1.seasonal_DE %>% 
  filter(Season =="Autumn")
ENS_seasonal_DE_aut <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_DE_aut) 
summary(ENS_seasonal_DE_aut)
tab_model(ENS_seasonal_DE_aut)
#winter:
data_1.seasonal_DE_win <- data_1.seasonal_DE %>% 
  filter(Season =="Winter")
ENS_seasonal_DE_win <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_DE_win) 
summary(ENS_seasonal_DE_win)
tab_model(ENS_seasonal_DE_win)

tab_model(ENS_seasonal_DE_spring,ENS_seasonal_DE_summer,
          ENS_seasonal_DE_aut, ENS_seasonal_DE_win, 
          show.ci = FALSE,  digits = 3,
          dv.labels = c("ENS.Spring-DE","ENS.Summer-DE", "ENS.Autumn-DE", "ENS.Winter-DE"),
          wrap.labels = 6) 


###2) NL

##a) S
data_1.seasonal_NL <- data_1.seasonal %>% 
  filter(Country=="Netherlands")
S_seasonal_NL <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_NL) 
summary(S_seasonal_NL)
tab_model(S_seasonal_NL)

#spring:
data_1.seasonal_NL_spring <- data_1.seasonal_NL %>% 
  filter(Season =="Spring")
S_seasonal_NL_spring <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_NL_spring) 
summary(S_seasonal_NL_spring)
tab_model(S_seasonal_NL_spring)
#summer:
data_1.seasonal_NL_summer <- data_1.seasonal_NL %>% 
  filter(Season =="Summer")
S_seasonal_NL_summer<- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_NL_summer) 
summary(S_seasonal_NL_summer)
tab_model(S_seasonal_NL_summer)
#autumn:
data_1.seasonal_NL_aut <- data_1.seasonal_NL %>% 
  filter(Season =="Autumn")
S_seasonal_NL_aut <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_NL_aut) 
summary(S_seasonal_NL_aut)
tab_model(S_seasonal_NL_aut)
#winter:
data_1.seasonal_NL_win <- data_1.seasonal_NL %>% 
  filter(Season =="Winter")
S_seasonal_NL_win <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.seasonal_NL_win) 
summary(S_seasonal_NL_win)
tab_model(S_seasonal_NL_win)

tab_model(S_seasonal_NL_spring,  S_seasonal_NL_summer,
          S_seasonal_NL_aut, S_seasonal_NL_win, 
          show.ci = FALSE,  digits = 3,
          dv.labels = c("S.Spring-NL","S.Summer-NL", "S.Autumn-NL", "S.Winter-NL"),
          wrap.labels = 6) 


## b) ENS
data_1.seasonal_NL <- data_1.seasonal %>% 
  filter(Country=="Netherlands")
ENS_seasonal_NL <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_NL) 
summary(ENS_seasonal_NL)
tab_model(ENS_seasonal_NL)

#spring:
data_1.seasonal_NL_spring <- data_1.seasonal_NL %>% 
  filter(Season =="Spring")
ENS_seasonal_NL_spring <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_NL_spring) 
summary(ENS_seasonal_NL_spring)
tab_model(ENS_seasonal_NL_spring)
#summer:
data_1.seasonal_NL_summer <- data_1.seasonal_NL %>% 
  filter(Season =="Summer")
ENS_seasonal_NL_summer<- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_NL_summer) 
summary(ENS_seasonal_NL_summer)
tab_model(ENS_seasonal_NL_summer)
#autumn:
data_1.seasonal_NL_aut <- data_1.seasonal_NL %>% 
  filter(Season =="Autumn")
ENS_seasonal_NL_aut <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_NL_aut) 
summary(ENS_seasonal_NL_aut)
tab_model(ENS_seasonal_NL_aut)
#winter:
data_1.seasonal_NL_win <- data_1.seasonal_NL %>% 
  filter(Season =="Winter")
ENS_seasonal_NL_win <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.seasonal_NL_win) 
summary(ENS_seasonal_NL_win)
tab_model(ENS_seasonal_NL_win)

tab_model(ENS_seasonal_NL_spring,ENS_seasonal_NL_summer,
          ENS_seasonal_NL_aut, ENS_seasonal_NL_win, 
          show.ci = FALSE,  digits = 3,
          dv.labels = c("ENS.Spring-NL","ENS.Summer-NL", "ENS.Autumn-NL", "ENS.Winter-NL"),
          wrap.labels = 6) 

### 1.1.7 Plots -----

# Change order of stations
data_1.seasonal$StationID <- factor(data_1.seasonal$StationID,
                                    c("MARSDND", "DOOVBWT","BOOMKDP", "TERSLG10", "DANTZGT", "ROTTMPT3", "HUIBGOT",  "BOCHTVWTM", "GROOTGND", 
                                      "Bork_W_1", "Nney_W_2", "JaBu_W_1", "WeMu_W_1"))
# Colors
Colors_coastal_NL_DE_2 <-  c("#fff5eb","#fee6ce","#fdd0a2", #"MARSDND", "DOOVBWT","BOOMKDP"
                                      "#fdae6b", #"TERSLG10"
                                      "#fd8d3c", # "DANTZGT"
                                      "#f16913",  #"ROTTMPT3"
                                      "#d94801", #"HUIBGOT"
                                      "#a63603", #"BOCHTVWTM"
                                      "#7f2704", #"GROOTGND",
                                      "#80cdc1", "#35978f", "#01665e", "#003c30")  # "Bork_W_1", "Nney_W_2", "JaBu_W_1", "WeMu_W_1"
                                      


# Richness (S)

# Plot with p-values added
p.S.seasonal <- ggplot(data_1.seasonal, aes(x = Year, y = median.Sy)) +
  geom_line(aes(col = StationID, linetype = Country), size = 0.8, alpha = 0.8, position = position_dodge(width = 0.3)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  #geom_line(data = fit.MELM.Sy_1, aes(y = fit), col = "black", size = 0.8) +
  #geom_ribbon(data = fit.MELM.Sy_1, aes(y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  scale_color_manual(values = c(Colors_coastal_NL_DE_2)) +
  ylab("Annual Richness") +
  xlab("Year") +
  theme_classic() +
  theme(text = element_text(size = 18), legend.position = "right") +
  theme(strip.text.x = element_text(size = 18)) +
  theme(legend.key.width = unit(1, "cm"), legend.box = "vertical", legend.margin = margin()) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  facet_wrap(~ Season, ncol = 2) 
# geom_text(aes(label = paste0("p-value = ", round(p_value, 3))), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 4, show.legend = FALSE)
p.S.seasonal 

setwd("~/Desktop/Job at AWI_HIFMB/nlwkn/my data/Interreg_Josie/InterReg-project/Manuscript/First manuscript_Temporal change/Repository_Git_TeporalChange/TemporalChange_PPKT_WaddenSea/Review1_R_script/Seasonal_plots")

png(filename="Richness_seasonal.png",
    type="cairo",
    units="in", 
    pointsize=14, 
    width=10, 
    height=8,  
    res=300)
print(p.S.seasonal)
dev.off()

# ENS

p.ENS.seasonal <- ggplot(data_1.seasonal, aes(x = Year, y = median.ENSy)) +
  geom_line(aes(col = StationID, linetype = Country), size = 0.8, alpha = 0.8, position = position_dodge(width = 0.3)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  #  geom_line(data = fit.MELM.ENSy_1, aes(y = fit), col = "black", size = 0.8) +
  #  geom_ribbon(data = fit.MELM.ENSy_1, aes(y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  scale_color_manual(values = c(Colors_coastal_NL_DE_2)) +
  ylab("Annual ENS")+
  xlab("Year") +
  theme_classic() +
  theme(text = element_text(size = 18), legend.position = "right") +
  theme(strip.text.x = element_text(size = 18)) +
  theme(legend.key.width = unit(1, "cm"), legend.box = "vertical", legend.margin = margin()) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  facet_wrap(~ Season, ncol = 2) 
#geom_text(aes(label = paste0("p-value = ", round(p_value, 3))), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 4, show.legend = FALSE)
p.ENS.seasonal 


setwd("~/Desktop/Job at AWI_HIFMB/nlwkn/my data/Interreg_Josie/InterReg-project/Manuscript/First manuscript_Temporal change/Repository_Git_TeporalChange/TemporalChange_PPKT_WaddenSea/Review1_R_script/Seasonal_plots")

png(filename="ENS_seasonal.png",
    type="cairo",
    units="in", 
    pointsize=14, 
    width=10, 
    height=8,  
    res=300)
print(p.ENS.seasonal)
dev.off()




# 3. Biomass and functional groups ----------------------

##3.1 Total biomass per sample -----
tot.biom<-ddply(PPKT_WS,.(Country, StationID, Date), 
                colwise(sum,  .(carbon_l)), na.rm=T)
str(tot.biom)

#add a year column 
tot.biom$Year<-year(tot.biom$Date)
#add a month  column 
tot.biom$Month<-month(tot.biom$Date)
names(tot.biom)

# Define a function to map month to season
get_season <- function(month) {
  if (month %in% c(12, 1, 2)) {
    return("Winter")
  } else if (month %in% c(3, 4, 5)) {
    return("Spring")
  } else if (month %in% c(6, 7, 8)) {
    return("Summer")
  } else {
    return("Autumn")
  }
}

# Convert the Month column to seasons
tot.biom$Season <- sapply(tot.biom$Month, get_season)


## 3.2 Calculate the mean and median biomass per year and season-----
tot.biom.season<- ddply(tot.biom, .(Country, StationID, Year, Season), 
                        summarise, 
                        mean.bm = mean(carbon_l), 
                        sd.bm = sd(carbon_l),
                        max.bm = max(carbon_l),
                        median.bm = median(carbon_l))


## 3.4 Plots -----

# Change order of stations
tot.biom.season$StationID <- factor(tot.biom.season$StationID,
                                    c("MARSDND", "DOOVBWT","BOOMKDP", "TERSLG10", "DANTZGT", "ROTTMPT3", "HUIBGOT",  "BOCHTVWTM", "GROOTGND", 
                                      "Bork_W_1", "Nney_W_2", "JaBu_W_1", "WeMu_W_1"))
names(tot.biom.season)

### Statistical analysis: ----------
# # Fit the linear mixed-effects model
# model <- lmer(median.bm ~ Year * Season + (1 | StationID) + (Year * Season | Country), data = tot.biom.season)

###DE
tot.biom.season_DE <- tot.biom.season %>% 
  filter(Country=="Germany")
tot.biom.DE <- lmer(log(median.bm+1) ~Season * Year + (1|StationID) + (0 + Year | Season), data=tot.biom.season_DE) 
summary(tot.biom.DE)
tab_model(tot.biom.DE)
#spring:
tot.biom.season_DE_spring <- tot.biom.season_DE %>% 
  filter(Season =="Spring")
tot.biom.DE_spring <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.season_DE_spring) 
summary(tot.biom.DE_spring)
tab_model(tot.biom.DE_spring) #-0.02	-0.07 – 0.03	0.403
#summer:
tot.biom.season_DE_summer <- tot.biom.season_DE %>% 
  filter(Season =="Summer")
tot.biom.DE_summer<- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.season_DE_summer) 
summary(tot.biom.DE_summer)
tab_model(tot.biom.DE_summer)#-0.01	-0.04 – 0.01	0.211
#autumn:
tot.biom.season_DE_aut <- tot.biom.season_DE %>% 
  filter(Season =="Autumn")
tot.biom.DE_aut <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.season_DE_aut) 
summary(tot.biom.DE_aut)
tab_model(tot.biom.DE_aut) #0.01	-0.03 – 0.04	0.692
#winter:
tot.biom.season_DE_win <- tot.biom.season_DE %>% 
  filter(Season =="Winter")
tot.biom.DE_win <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.season_DE_win) 
summary(tot.biom.DE_win)
tab_model(tot.biom.DE_win)#Year	-0.09	-0.12 – -0.07	<0.001

tab_model(tot.biom.DE_spring,  tot.biom.DE_summer,
          tot.biom.DE_aut, tot.biom.DE_win, 
          show.ci = FALSE,  digits = 3,
          dv.labels = c("BM.Spring-DE","BM.Summer-DE", "BM.Autumn-DE", "BM.Winter-DE"),
          wrap.labels = 6) 


###NL
tot.biom.season_NL <- tot.biom.season %>% 
  filter(Country=="Netherlands")
tot.biom.NL <- lmer(log(median.bm+1) ~Season + (1|StationID) + (0 + Year | Season), data=tot.biom.season_NL) 
summary(tot.biom.NL)
tab_model(tot.biom.NL)
#spring:
tot.biom.season_NL_spring <- tot.biom.season_NL %>% 
  filter(Season =="Spring")
tot.biom.NL_spring <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.season_NL) 
summary(tot.biom.NL_spring)
tab_model(tot.biom.NL_spring) #0.10	0.09 – 0.12	<0.001****
#summer:
tot.biom.season_NL_summer <- tot.biom.season_NL %>% 
  filter(Season =="Summer")
tot.biom.NL_summer<- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.season_NL_summer) 
summary(tot.biom.NL_summer)
tab_model(tot.biom.NL_summer)#0.07	0.05 – 0.09	<0.001***
#autumn:
tot.biom.season_NL_aut <- tot.biom.season_NL %>% 
  filter(Season =="Autumn")
tot.biom.NL_aut <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.season_NL_aut) 
summary(tot.biom.NL_aut)
tab_model(tot.biom.NL_aut) #0.11	0.08 – 0.13	<0.001***
#winter:
tot.biom.season_NL_win <- tot.biom.season_NL %>% 
  filter(Season =="Winter")
tot.biom.NL_win <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.season_NL_win) 
summary(tot.biom.NL_win)
tab_model(tot.biom.NL_win)#0.15	0.11 – 0.19	<0.001***

tab_model(tot.biom.NL_spring,  tot.biom.NL_summer,
          tot.biom.NL_aut, tot.biom.NL_win, 
          show.ci = FALSE,  digits = 3,
          dv.labels = c("BM.Spring-NL","BM.Summer-NL", "BM.Autumn-NL", "BM.Winter-NL"),
          wrap.labels = 6) 


### 3.4.1 Median Biomass -----
p.bm_log <- ggplot(tot.biom.season, aes(x=Year, y=log(median.bm+1))) +
  geom_line(aes(Year,log(median.bm+1), col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
  labs(y = expression(paste("Biomass ", "[LN ","µgC ", L^{-1}, "]")))+
  xlab("year")+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = "right") +
  theme(strip.text.x = element_text(size = 18))  +
  theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  facet_wrap(~ Season, ncol = 2) 
#geom_text(aes(label = paste0("p-value = ", round(p_value, 3))), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 4, show.legend = FALSE)
p.bm_log

setwd("~/Desktop/Job at AWI_HIFMB/nlwkn/my data/Interreg_Josie/InterReg-project/Manuscript/First manuscript_Temporal change/Repository_Git_TeporalChange/TemporalChange_PPKT_WaddenSea/Review1_R_script/Seasonal_plots")

png(filename="LNBiomass_seasonal.png",
    type="cairo",
    units="in", 
    pointsize=14, 
    width=10, 
    height=8,  
    res=300)
print(p.bm_log)
dev.off()

##3.3 Functional group: Calculate the mean and median biomass per sample per functional group -----
#sum per sample
tx.sample<-ddply(PPKT_WS,.(Country, StationID, Date, Functional_group), 
                 colwise(sum,  .(carbon_l)), na.rm=T)
# Cyanobateria considered as "other"
tx.sample <- tx.sample %>% 
  mutate(Functional_group = case_when(
    Functional_group == "Cyanobacteria" ~ "other", 
    Functional_group == "Diatoms" ~ "Diatoms", 
    Functional_group == "Dinoflagellates" ~ "Dinoflagellates", 
    Functional_group == "Flagellates" ~ "Flagellates", 
    Functional_group == "Phaeocystis" ~ "Phaeocystis", 
    Functional_group == "other" ~ "other"
  ))    

### 3.3.1 Calculate the median func. group biomass per year -----

#add a month  column 
tx.sample$Month<-month(tx.sample$Date)
names(tot.biom)

# Define a function to map month to season
get_season <- function(month) {
  if (month %in% c(12, 1, 2)) {
    return("Winter")
  } else if (month %in% c(3, 4, 5)) {
    return("Spring")
  } else if (month %in% c(6, 7, 8)) {
    return("Summer")
  } else {
    return("Autumn")
  }
}

# Convert the Month column to seasons
tx.sample$Season <- sapply(tx.sample$Month, get_season)

#add a year column 
tx.sample$Year<-year(tx.sample$Date)


# median per year
tx.seasonal.med<-ddply(tx.sample,.(Country, StationID, Year, Season, Functional_group), 
                       colwise(median,  .(carbon_l)), na.rm=T)

# Change order of stations
tx.seasonal.med$StationID <- factor(tx.seasonal.med$StationID,
                                    c("MARSDND", "DOOVBWT","BOOMKDP", "TERSLG10", "DANTZGT", "ROTTMPT3", "HUIBGOT",  "BOCHTVWTM", "GROOTGND", 
                                      "Bork_W_1", "Nney_W_2", "JaBu_W_1", "WeMu_W_1"))
tx.seasonal.med$Functional_group <- factor(tx.seasonal.med$Functional_group,
                                           c("Diatoms", "Dinoflagellates","Flagellates", "Phaeocystis", "other"))
# Colors
Colors_coastal_NL_DE_2 <-  c("#fff5eb","#fee6ce","#fdd0a2", #"MARSDND", "DOOVBWT","BOOMKDP"
                                      "#fdae6b", #"TERSLG10"
                                      "#fd8d3c", # "DANTZGT"
                                      "#f16913",  #"ROTTMPT3"
                                      "#d94801", #"HUIBGOT"
                                      "#a63603", #"BOCHTVWTM"
                                      "#7f2704", #"GROOTGND",
                                      "#80cdc1", "#35978f", "#01665e", "#003c30")  # "Bork_W_1", "Nney_W_2", "JaBu_W_1", "WeMu_W_1"
                                      

###Plot ------
tx.seasonal.med_1 = filter(tx.seasonal.med, Functional_group !="other")

p.fg <- ggplot(tx.seasonal.med, aes(x=Year, y=log(carbon_l+1))) +
  geom_line(aes(Year,log(carbon_l+1), col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
  labs(y = expression(paste("Biomass ", "[LN ","µgC ", L^{-1}, "]")))+
  xlab("year")+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = "right") +
  theme(strip.text.x = element_text(size = 18))  +
  theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  #facet_grid(Season ~ Functional_group, margins=TRUE)
  facet_grid(Functional_group ~ Season, scales='free')
p.fg


setwd("~/Desktop/Job at AWI_HIFMB/nlwkn/my data/Interreg_Josie/InterReg-project/Manuscript/First manuscript_Temporal change/Repository_Git_TeporalChange/TemporalChange_PPKT_WaddenSea/Review1_R_script/Seasonal_plots")

png(filename="LNBiomass_FG_seasonal.png",
    type="cairo",
    units="in", 
    pointsize=14, 
    width=12, 
    height=9,  
    res=300)
print(p.fg)
dev.off()

