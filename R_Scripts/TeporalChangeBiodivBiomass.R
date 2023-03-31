# Manuscript: Temporal change in phytoplankton diversity and FG composition


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


# Data ----------------------
#download phytoplankton data for the Wadden Sea

PPKT_WS <- read_delim("PPKT_count_WaddenSea_1999_2018.csv", 
                          delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                          locale = locale(decimal_mark = ",", grouping_mark = "."), 
                          trim_ws = TRUE)
#How many stations?
length(unique(PPKT_WS$StationID)) #13 stations


# Data Analysis ----------------------

# 1. Standing Diversity ----------------------


## 1.1. Gamma diversity: First calculate the median abundance of each species per station, per year, then calculate diversity -------
# ENSy = for each species, take mean biomass per year and calculate ENS then
# Sy = for each species, take mean biomass per year and calculate S then

### 1.1.1 Calculate the mean biomass per species, per year -----
PPKT_WS.median.annual<-ddply(PPKT_WS,.(Country, StationID, Year, Species), 
                        colwise(mean,  .(abundance_l)), na.rm=T)

### 1.1.2 Converting the data Between Wide And Long Forms -----
PPKT_WS.median.annual_wide<-dcast(PPKT_WS.median.annual, Country + StationID + Year ~ Species, value.var="abundance_l")

### 1.1.3 Transforming all NAs in 0: in all rows -----
PPKT_WS.median.annual_wide[,4:442][is.na(PPKT_WS.median.annual_wide[,4:442])]<-0

data_1 <- PPKT_WS.median.annual_wide

# number of columns
lastspec2<-ncol(data_1)

### 1.1.4 Calculating Richness and ESN per year -----

##S
data_1$median.ENSy<-diversity(data_1[,4:lastspec2], "invsimpson") 
##ENS
data_1$median.Sy<-specnumber(data_1[,4:lastspec2])  

### 1.1.5 Select columns -----
data_1.annual<-data_1[,c("Country","StationID","Year","median.Sy","median.ENSy")]

### 1.1.6 Statistical analysis -----
#### a) General slope ----
#S:
MELM.Sy_1 <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.annual) 
summary(MELM.Sy_1)
tab_model(MELM.Sy_1) #cond. R2 = 0.648 slope -0.86, p <0.001
fit.MELM.Sy_1 = as.data.frame(Effect(c("Year"),MELM.Sy_1))
coef(summary(MELM.Sy_1))

#ENS
MELM.ENSy_1 <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.annual) 
summary(MELM.ENSy_1)
tab_model(MELM.ENSy_1) #cond. R2 = 0.257 slope -0.13, p =0.009
fit.MELM.ENSy_1 = as.data.frame(Effect(c("Year"),MELM.ENSy_1))
coef(summary(MELM.ENSy_1))

#### b) Separate models per Country (In the Appendix) -----

#German Stations
  #S
data_1.annual_DE <- data_1.annual %>% 
  filter(Country=="Germany")
MELM.Sy_DE <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.annual_DE) 
summary(MELM.Sy_DE)
tab_model(MELM.Sy_DE)
  #ENS
MELM.ENSy_DE <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.annual_DE) 
summary(MELM.ENSy_DE)
tab_model(MELM.ENSy_DE)

#Dutch Stations
  #S
data_1.annual_NL <- data_1.annual %>% 
  filter(Country=="Netherlands")
MELM.Sy_NL <- lmer(median.Sy ~ Year + (1|StationID), data=data_1.annual_NL) 
summary(MELM.Sy_NL)
tab_model(MELM.Sy_NL)
  #ENS
MELM.ENSy_NL <- lmer(median.ENSy ~ Year + (1|StationID), data=data_1.annual_NL) 
summary(MELM.ENSy_NL)
tab_model(MELM.ENSy_NL)


### 1.1.7 Plots -----

# Change order of stations
data_1.annual$StationID <- factor(data_1.annual$StationID,
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
p.S <- ggplot(data_1.annual, 
         aes(x=Year, y=median.Sy)) +
  geom_line(aes(Year,median.Sy, col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  geom_line(data = fit.MELM.Sy_1, aes(y=fit), col="black", size=0.8) +
  geom_ribbon(data = fit.MELM.Sy_1, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
  scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
  ylab("Annual Richness")+
  xlab("year")+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = "right") +
  theme(strip.text.x = element_text(size = 18))  +
  theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
p.S

p.S_lm <- p.S + annotate("text", x = 2016, y = 125, label = "p<0.001",
           parse = TRUE, size = 4) 
p.S_lm

# ENS 
p.ENS <- ggplot(data_1.annual, 
              aes(x=Year, y=median.ENSy)) +
  geom_line(aes(Year,median.ENSy, col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  geom_line(data = fit.MELM.ENSy_1, aes(y=fit), col="black", size=0.8) +
  geom_ribbon(data = fit.MELM.ENSy_1, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
  scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
  ylab("Annual ENS")+
  xlab("year")+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = "right") +
  theme(strip.text.x = element_text(size = 18))  +
  theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
p.ENS

p.ENS_lm <- p.ENS + annotate("text", x = 2000, y = 18, label = "p==0.009",
                             parse = TRUE, size = 4) 
p.ENS_lm


# 2. Turnover ----------------------

## 2.1. mean of the biomass per station, year and species -----
data_mean.annual<-ddply(PPKT_WS,.(Country, StationID, Year, Species), 
                        colwise(mean,  .(carbon_l)), na.rm=T)

## 2.2. dcast: Convert Data Between Wide And Long Forms.
data_mean.annual_wideBM<-dcast(data_mean.annual,Country + StationID + Year ~ Species, value.var="carbon_l")

## 2.3. Transforming all NAs in 0: 
data_mean.annual_wideBM[,4:442][is.na(data_2_mean.annual_wideBM[,4:442])]<-0

data_2 <- data_mean.annual_wideBM

## 2.4 Building the function to calculate turnover -----
length(unique(data_2$StationID)) #13

max(data_2$Year)
#[1] 2018
min(data_2$Year)
#[1] 1999
max(data_2$Year)-min(data_2$Year)+1
#[1] 20
#number of rows
nrow(data_2)
#[1] 213
#number of columns with Species names
ncol(data_2)-3
#[1]439

antal<-ncol(data_2) #312L - data length
#selecting the names of the columns with species names.
nms<-names(data_2)[4:dim(data_2)[2]]
#unique function in R –unique(), eliminates duplicate elements/rows from a vector, data frame or array.
#selects the Stations 
st<-sort(unique(data_2$StationID))
#selects only the first station
d<-subset(data_2,StationID==st[1]) 
#selects only the first year of the first station 
d<-subset(d,Year==sort(unique(d$Year))[1]) 

## 2.5 Function to calculate turnover -----
f<-function(d){
  print(unique(d$StationID))
  #SELECT ALL DATA IN YEAR RANGE
  d2<-subset(data_2,StationID==unique(d$StationID) & Year>=unique(d$Year))
  
  #IF ENOUGH DATA, PROCEED
  if(dim(d2)[1]>1){
    
    #LOOP THROUGH ALL POSSIBLE YEAR COMBINATIONS
    yrs<-sort(unique(d2$Year))[-1]
    l<-list()
    for(i in 1:length(yrs)){
      d3<-rbind(subset(d2,Year==unique(d$Year)),subset(d2,Year==yrs[i]))
      d3<-d3[order(d3$Year,decreasing=TRUE),]
      
      if(dim(d3)[1]>1){
        #SUBTRACT VALUE IN SECOND YEAR FROM VALUE IN FIRST YEAR
        datout<-data.frame(t((delta=apply((d3)[,3:dim(d3)[2]],2,diff))/(apply((d3)[,3:dim(d3)[2]],2,sum)))) ##diff: Lagged Differences
        #ADD SPECIES NAMES
        names(datout)<-nms
        #ADD AVERAGE AND DELTA YEAR
        datout$avgyear<-(min(d2$Year)+yrs[i])/2
        datout$deltayear<-yrs[i]-min(d2$Year)
        ##calculate Simpson indexes and diversity before, after and mixed
        datout$SimpBefore <- rowSums(d3[2, 3:antal]^2)  #Simpson index and diversity Year 1
        datout$DivBefore  <- 1/datout$SimpBefore
        datout$SimpAfter  <- rowSums(d3[1, 3:antal]^2)  #Simpson index and diversity Year 2
        datout$DivAfter   <-  1/datout$SimpAfter
        #Simpson index Year 1-2 (mixed)
        datout$SimpMix    <- rowSums(d3[1, 3:antal] * d3[2, 3:antal], na.rm = FALSE, dims = 1)
        #Effective number of common species
        datout$DivCommon  <-  datout$SimpMix/(datout$SimpBefore * datout$SimpAfter)
        #Effective number of emigrating species
        datout$DivExt     <- datout$DivBefore - datout$DivCommon
        #Effective number of immigrants
        datout$DivImm     <- datout$DivAfter - datout$DivCommon
        #Effective total species number
        datout$DivTot     <- datout$DivBefore + datout$DivAfter - datout$DivCommon
        #Generalized SER index
        datout$TORpie     <- (datout$DivImm + datout$DivExt)/datout$DivTot
        l[[i]]<-datout
      } else NULL
    }
    return(data.frame(do.call('rbind',l)))
  } else NULL
}

dataout<-ddply(data_2,.(StationID,Year),.fun=f,.progress='text') 
names(dataout) 

## 2.5 getting the diversity estimates -----

# Use 0 or 1 here to calculate using richness (1) or Simpson (0)
# use richness 

# number of immigrating species, replace selection of columns
# note: it is -1 because we substract year1-year2!!!!
dataout$Imm <- rowSums(dataout[,3:antal] == '-1')
# number of emigrating species
dataout$Ext <- rowSums(dataout[,3:antal] == '1')
# total number of species
dataout$totS<- rowSums(!is.na(dataout[,3:antal]))
# number of persisting species
dataout$Spers<- dataout$totS-dataout$Imm-dataout$Ext
# species turnover ratio, between 0 and 1
dataout$TOR<- (dataout$Imm+dataout$Ext)/dataout$totS
# change in species richness
dataout$deltaS<-(dataout$Imm-dataout$Ext)

### 2.5.1 record the distance in years between samples -----
dataout$dist <- dataout$deltayear
names(dataout)

## Legend - columns in dataout --------
# abundance based indexes  ??? was not bases on biomass?
# SimpBefore, SimpAfter  simpsons index in year1 and year 2
# DivBefore,  DivAfter   effective species (ESN) number (1/Simppson) in year1 and year 2
# SimpMix      "mixed" simpson index = \sum_i p_i p'_i  
# DivCommon   ESN of common species                                
# DivExt      ESN, extinct
# DivImm      ESN Immigrated
# DivTot      ESN total                                   
# TORpie      SERa, abundance based turnover rate
# 
# richness based indexes
# Imm         number of immigrating species
# Ext         number of extinct species                                      
# totS        total number of species
# Spers       number of persisting (common) species
# TOR         SERr, richness based turnover rate
# deltaS      change in species richness
# Int.sp      change in species proportion of persisting species, between 0 and 1
# dist        the distance in years between samples 

# summarize the information into a new data frame
turnover1<-dataout[,c("StationID", "Year","Imm" ,"Ext" ,"totS","Spers", "TOR", "deltaS","dist")]
# change in species richness based on species abundances
dataout$DivdeltaS<-(dataout$DivImm-dataout$DivExt)

# record the distance in years between samples
dataout$dist <- dataout$deltayear
turnover2<-dataout[,c("StationID", "Year","TORpie", "dist")]

#merge files
data.turnover <-merge(turnover1, turnover2, by=c("StationID", "Year", "dist"))

### 2.5.2 Immediate turnover: select the rows with dist less than 2 -----
data.turnover_annual<-data.turnover[data.turnover$dist<2,]

#rename columns
data.turnover_annual <- data.turnover_annual %>% 
  rename(SERr = TOR, 
         SERa = TORpie)

# Turnover file: selecting columns 
data.turnover_annual<-data.turnover_annual[,c("StationID", "Year", "SERr","SERa","Imm","Ext")]

### 2.5.3 Cumulative turnover -----
#rename columns
data.turnover_cum <- data.turnover %>% 
  rename(SERr = TOR, 
         SERa = TORpie)

### add a column with the country names: 
#data.turnover_annual <- data.turnover_annual %>% 
data.turnover_cum <- data.turnover_cum %>% 
  mutate(Country = case_when(
    StationID == "MARSDND" ~ "Netherlands",
    StationID == "DOOVBWT" ~ "Netherlands",
    StationID == "BOOMKDP" ~ "Netherlands",
    StationID == "TERSLG10" ~ "Netherlands",
    StationID == "DANTZGT" ~ "Netherlands",
    StationID == "ROTTMPT3" ~ "Netherlands",
    StationID == "HUIBGOT" ~ "Netherlands",
    StationID == "BOCHTVWTM" ~ "Netherlands", 
    StationID == "GROOTGND" ~ "Netherlands",
    
    StationID == "Bork_W_1" ~ "Germany",
    StationID == "Nney_W_2" ~ "Germany",
    StationID == "JaBu_W_1" ~ "Germany",
    StationID == "WeMu_W_1" ~ "Germany"))

## 2.6 Statistical analysis -----
 
  ### 2.6.1 Immediate turnover -----
  #### SERa -----
    #### a) General slope ----
    MELM.SERa_1 <- lmer(SERa ~ Year + (1|StationID), data=data.turnover_annual) 
    summary(MELM.SERa_1)
    tab_model(MELM.SERa_1) #cond. R2 = 0.390 slope 0.00 p=0.665
    fit.MELM.SERa = as.data.frame(Effect(c("Year"),MELM.SERa_1))
    coef(summary(MELM.SERa_1))
    
    #### b) Separate models per Country (In the Appendix) ----
    
    #DE
    data.turnover_annual_DE <- data.turnover_annual %>% 
      filter(Country=="Germany")
    MELM.SERa_DE_4.1 <- lmer(SERa ~ Year + (1|StationID), data=data.turnover_annual_DE) 
    summary(MELM.SERa_DE_4.1)
    tab_model(MELM.SERa_DE_4.1)
    #NL
    data.turnover_annual_NL <- data.turnover_annual %>% 
      filter(Country=="Netherlands")
    MELM.SERa_NL_4.1 <- lmer(SERa ~ Year + (1|StationID), data=data.turnover_annual_NL) 
    summary(MELM.SERa_NL_4.1)
    tab_model(MELM.SERa_NL_4.1)
    
  #### SERr -----
    #### a) General slope ----
    MELM.SERr_1 <- lmer(SERr ~ Year + (1|StationID), data=data.turnover_annual) 
    summary(MELM.SERr_1)
    tab_model(MELM.SERr_1) #cond. R2 = 0.390 slope 0.00 p=0.665
    fit.MELM.SERr = as.data.frame(Effect(c("Year"),MELM.SERr_1))
    coef(summary(MELM.SERr_1))
    
    #### b) Separate models per Country (In the Appendix) ----
    
    #DE
    data.turnover_annual_DE <- data.turnover_annual %>% 
      filter(Country=="Germany")
    MELM.SERr_DE_4.1 <- lmer(SERr ~ Year + (1|StationID), data=data.turnover_annual_DE) 
    summary(MELM.SERr_DE_4.1)
    tab_model(MELM.SERr_DE_4.1)
    #NL
    data.turnover_annual_NL <- data.turnover_annual %>% 
      filter(Country=="Netherlands")
    MELM.SERr_NL_4.1 <- lmer(SERr ~ Year + (1|StationID), data=data.turnover_annual_NL) 
    summary(MELM.SERr_NL_4.1)
    tab_model(MELM.SERr_NL_4.1)
    
    
  ### 2.6.2 Cumulative turnover -----
    #### SERa -----
    #### a) General slope ----
    MELM.SERa_1 <- lmer(SERa ~ dist + (1|StationID), data=data.turnover_cum) 
    summary(MELM.SERa_1)
    tab_model(MELM.SERa_1) 
    fit.MELM.SERa = as.data.frame(Effect(c("dist"),MELM.SERa_1))
    coef(summary(MELM.SERa_1))
    
    #### b) Separate models per Country (In the Appendix) ----
    #DE
    data.turnover_cum_DE <- data.turnover_cum %>% 
      filter(Country=="Germany")
    MELM.SERa_DE_4.1 <- lmer(SERa ~ dist + (1|StationID), data=data.turnover_cum_DE) 
    summary(MELM.SERa_DE_4.1)
    tab_model(MELM.SERa_DE_4.1)
    #NL
    data.turnover_cum_NL <- data.turnover_cum %>% 
      filter(Country=="Netherlands")
    MELM.SERa_NL_4.1 <- lmer(SERa ~ dist + (1|StationID), data=data.turnover_cum_NL) 
    summary(MELM.SERa_NL_4.1)
    tab_model(MELM.SERa_NL_4.1)
    
    #### SERr -----
    #### a) General slope ----
    MELM.SERr_1 <- lmer(SERr ~ dist + (1|StationID), data=data.turnover_cum) 
    summary(MELM.SERr_1)
    tab_model(MELM.SERr_1) 
    fit.MELM.SERr = as.data.frame(Effect(c("dist"),MELM.SERr_1))
    coef(summary(MELM.SERr_1))
    
    #### b) Separate models per Country (In the Appendix) ----
    #DE
    data.turnover_cum_DE <- data.turnover_cum %>% 
      filter(Country=="Germany")
    MELM.SERr_DE_4.1 <- lmer(SERr ~ dist + (1|StationID), data=data.turnover_cum_DE) 
    summary(MELM.SERr_DE_4.1)
    tab_model(MELM.SERr_DE_4.1)
    #NL
    data.turnover_cum_NL <- data.turnover_cum %>% 
      filter(Country=="Netherlands")
    MELM.SERr_NL_4.1 <- lmer(SERr ~ dist + (1|StationID), data=data.turnover_cum_NL) 
    summary(MELM.SERr_NL_4.1)
    tab_model(MELM.SERr_NL_4.1)
    
    
## 2.7 Plots -----

### 2.7.1 Immediate turnover -----
# Change order of stations
data.turnover_annual$StationID <- factor(data.turnover_annual$StationID,
                                  c("MARSDND", "DOOVBWT","BOOMKDP", "TERSLG10", "DANTZGT", "ROTTMPT3", "HUIBGOT",  "BOCHTVWTM", "GROOTGND", 
                                    "Bork_W_1", "Nney_W_2", "JaBu_W_1", "WeMu_W_1"))

#SERa
    p.SERa <- 
      ggplot(data.turnover_annual, 
             aes(x=Year, y=SERa)) +
      geom_line(aes(Year,SERa, col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
      scale_linetype_manual(values=c("dashed", "solid"))+
      geom_line(data = fit.MELM.SERa, aes(y=fit), col="black", size=0.8) +
      geom_ribbon(data = fit.MELM.SERa, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
      # scale_color_manual(values=c(Colors_NL_DE_8))+  
      scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
      ylab("Annual SERa")+
      xlab("year")+
      theme_classic()+
      theme(text = element_text(size = 18), 
            legend.position = "right") +
      theme(strip.text.x = element_text(size = 18))  +
      theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
      theme(panel.background = element_rect(fill = "white", colour = "black"))
    p.SERa
    p.SERa_lm <- p.SERa + 
      annotate("text", x = 2001, y = 0.05, label = "p==0.66",
               parse = TRUE, size = 4) 
    p.SERa_lm

#SERr
    p.SERr <- 
      ggplot(data.turnover_annual, 
             aes(x=Year, y=SERr)) +
      geom_line(aes(Year,SERr, col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
      scale_linetype_manual(values=c("dashed", "solid"))+
      geom_line(data = fit.MELM.SERr, aes(y=fit), col="black", size=0.8) +
      geom_ribbon(data = fit.MELM.SERr, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
      # scale_color_manual(values=c(Colors_NL_DE_8))+  
      scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
      ylab("Annual SERr")+
      xlab("year")+
      theme_classic()+
      theme(text = element_text(size = 18), 
            legend.position = "right") +
      theme(strip.text.x = element_text(size = 18))  +
      theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
      theme(panel.background = element_rect(fill = "white", colour = "black"))
    p.SERr
    p.SERr_lm <- p.SERr +  annotate("text", x = 2001, y = 0.24, label = "p==0.41",
               parse = TRUE, size = 4) 
    p.SERr_lm

### 2.7.2 Cumulative turnover -----

#SERa
    p.SERa_cum <- ggplot(data.turnover_cum, aes(x=dist, y=SERa)) +
      geom_point(aes(x=dist, y=SERa, col=StationID, shape=Country, group=Country),size=1.5, alpha=.4, position="jitter", stroke=.8) +
      geom_line(data = fit.MELM.SERa, aes(y=fit), col="black", size=0.8) +
      geom_ribbon(data = fit.MELM.SERa, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
      scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
      ylab("Cumulative SERa")+
      xlab("Temporal distance [years]")+
      theme_classic()+
      theme(text = element_text(size = 18), 
            legend.position = "right") +
      theme(strip.text.x = element_text(size = 18))  +
      theme(panel.background = element_rect(fill = "white", colour = "black"))
    p.SERa_cum
    p.SERa_lm_cum <- p.SERa_cum_1 + 
      annotate("text", x = 17, y = 0.1, label = "p<0.001",
               parse = TRUE, size = 4) 
    p.SERa_lm_cum

#SERr
    
    p.SERr_cum <- ggplot(data.turnover_cum, aes(x=dist, y=SERr)) +
      geom_point(aes(x=dist, y=SERr, col=StationID, shape=Country, group=Country),size=1.5, alpha=.4, position="jitter", stroke=.8) +
      geom_line(data = fit.MELM.SERr, aes(y=fit), col="black", size=0.8) +
      geom_ribbon(data = fit.MELM.SERr, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
      
      scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
      ylab("Cumulative SERr")+
      xlab("Temporal distance [years]")+
      theme_classic()+
      theme(text = element_text(size = 18), 
            legend.position = "right") +
      theme(strip.text.x = element_text(size = 18))  +
      theme(panel.background = element_rect(fill = "white", colour = "black"))
    p.SERr_cum
    p.SERr_lm_cum <- p.SERr_cum_1 + 
      annotate("text", x = 17, y = 0.24, label = "p<0.001",
               parse = TRUE, size = 4) 
    p.SERr_lm_cum
    
# 3. Biomass and functional groups ----------------------

    ##3.1 Total biomass per sample -----
        tot.biom<-ddply(PPKT_WS,.(Country, StationID, Date), 
                    colwise(sum,  .(carbon_l)), na.rm=T)
    
    #Add a sample date column: 
    tot.biom$sampledate<-dmy(tot.biom$Date)
    #add a year column 
    tot.biom$Year<-year(tot.biom$sampledate)
    #add a month  column 
    tot.biom$month<-month(tot.biom$sampledate, label = TRUE, abbr = TRUE)
    names(tot.biom)
    
    ## 3.2 Calculate the mean and median biomass per year -----
    tot.biom.annual<- ddply(tot.biom, .(Country, StationID, Year), 
                            summarise, 
                            mean.bm = mean(carbon_l), 
                            sd.bm = sd(carbon_l),
                            max.bm = max(carbon_l),
                            median.bm = median(carbon_l))
    
        ### 3.2.1 Statistical analysis -----
        MELM.bio_log <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.annual) 
        summary(MELM.bio_log)
        tab_model(MELM.bio_log) #cond. R2 = 0.870, slope 0.06***, p<0.001
        coef(summary(MELM.bio_log))
        fit.MELM.bio_log = as.data.frame(Effect(c("Year"), MELM.bio_log)) ##worked!!!!
        
        #one model per country
        tot.biom.annual_NL <- tot.biom.annual %>% 
          filter(Country=="Netherlands")
        MELM.bio_log_NL <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.annual_NL) 
        summary(MELM.bio_log_NL)
        tab_model(MELM.bio_log_NL) 
        
        tot.biom.annual_DE <- tot.biom.annual %>% 
          filter(Country=="Germany")
        MELM.bio_log_DE <- lmer(log(median.bm+1) ~ Year + (1|StationID), data=tot.biom.annual_DE) 
        summary(MELM.bio_log_DE)
        tab_model(MELM.bio_log_DE) 
    
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
    
        #Add a sample date column: 
        tx.sample$sampledate<-dmy(tx.sample$Date)
        #add a year column 
        tx.sample$Year<-year(tx.sample$sampledate)
        #month
        tx.sample$month<-month(tx.sample$sampledate, label = TRUE, abbr = TRUE)
    
        # median per year
        tx.year.med<-ddply(tx.sample,.(Country, StationID, Year, Functional_group), 
                       colwise(median,  .(carbon_l)), na.rm=T)

        # Converting the data Between Wide And Long Forms
        #dcast: 
        tx.year.med_wideBM<-dcast(tx.year.med, Country + StationID + Year ~ Functional_group, value.var="carbon_l")
        
        # Transforming all NAs in 0: 
        tx.year.med_wideBM[,4:8][is.na(tx.year.med_wideBM[,4:8])]<-0
        
        #rename:
        data_median_fg <- tx.year.med_wideBM %>%
          rename(median.diatoms = Diatoms,
                 median.dino = Dinoflagellates,
                 median.fla = Flagellates, 
                 median.phae = Phaeocystis, 
                 median.other = other)  

        ### 3.2.2 Statistical analysis -----
        
          ####Diatoms -----
          ##### a) General slope -----
          MELM.Diatoms_log <- lmer(log(median.diatoms+1) ~ Year + (1|StationID), data=data_median_fg) 
          summary(MELM.Diatoms_log)
          tab_model(MELM.Diatoms_log) 
          coef(summary(MELM.Diatoms_log))
          fit.MELM.Diatoms_log = as.data.frame(Effect(c("Year"), MELM.Diatoms_log)) 
          
          ##### b) Separate models per Country (In the Appendix) ------
          data_median_fg_NL <- data_median_fg %>% 
            filter(Country=="Netherlands")
          MELM.Diatoms_log_NL <- lmer(log(median.diatoms+1) ~ Year + (1|StationID), data=data_median_fg_NL) 
          summary(MELM.Diatoms_log_NL)
          tab_model(MELM.Diatoms_log_NL) 
          data_median_fg_DE <- data_median_fg %>% 
            filter(Country=="Germany")
          MELM.Diatoms_log_DE <- lmer(log(median.diatoms+1) ~ Year + (1|StationID), data=data_median_fg_DE) 
          summary(MELM.Diatoms_log_DE)
          tab_model(MELM.Diatoms_log_DE)

          ####Dinoflagellates -----
          ##### a) General slope -----
          MELM.Dino_log <- lmer(log(median.dino+1) ~ Year + (1|StationID), data=data_median_fg) 
          summary(MELM.Dino_log)
          tab_model(MELM.Dino_log) 
          fit.MELM.Dino_log = as.data.frame(Effect(c("Year"), MELM.Dino_log)) 
          
          ##### b) Separate models per Country (In the Appendix) ------
          MELM.Dino_log_NL <- lmer(log(median.dino+1) ~ Year + (1|StationID), data=data_median_fg_NL) 
          summary(MELM.Dino_log_NL)
          tab_model(MELM.Dino_log_NL)
        
          MELM.Dino_log_DE <- lmer(log(median.dino+1) ~ Year + (1|StationID), data=data_median_fg_DE) 
          summary(MELM.Dino_log_DE)
          tab_model(MELM.Dino_log_DE) 
          
          ####Flagellates -----
          ##### a) General slope -----
          MELM.Flagellates_log <- lmer(log(median.fla+1) ~ Year + (1|StationID), data=data_median_fg) 
          summary(MELM.Flagellates_log)
          tab_model(MELM.Flagellates_log) 
          fit.MELM.Flagellates_log = as.data.frame(Effect(c("Year"), MELM.Flagellates_log)) 
          
          ##### b) Separate models per Country (In the Appendix) ------
          MELM.Flagellates_log_NL <- lmer(log(median.fla+1) ~ Year + (1|StationID), data=data_median_fg_NL) 
          summary(MELM.Flagellates_log_NL)
          tab_model(MELM.Flagellates_log_NL)
          
          MELM.Flagellates_log_DE <- lmer(log(median.fla+1) ~ Year + (1|StationID), data=data_median_fg_DE) 
          summary(MELM.Flagellates_log_DE)
          tab_model(MELM.Flagellates_log_DE) 
          
          ####Phaeocystis -----
          ##### a) General slope -----
          MELM.Phaeocystis_log <- lmer(log(median.phae+1) ~ Year + (1|StationID), data=data_median_fg) 
          summary(MELM.Phaeocystis_log)
          tab_model(MELM.Phaeocystis_log) 
          fit.MELM.Phaeocystis_log = as.data.frame(Effect(c("Year"), MELM.Phaeocystis_log)) 
          
          ##### b) Separate models per Country (In the Appendix) ------
          MELM.Phaeocystis_log_NL <- lmer(log(median.phae+1) ~ Year + (1|StationID), data=data_median_fg_NL) 
          summary(MELM.Phaeocystis_log_NL)
          tab_model(MELM.Phaeocystis_log_NL)
          
          MELM.Phaeocystis_log_DE <- lmer(log(median.phae+1) ~ Year + (1|StationID), data=data_median_fg_DE) 
          summary(MELM.Phaeocystis_log_DE)
          tab_model(MELM.Phaeocystis_log_DE) 
        
          ####Other -----
          ##### a) General slope -----
          MELM.Other_log <- lmer(log(median.other+1) ~ Year + (1|StationID), data=data_median_fg) 
          summary(MELM.Other_log)
          tab_model(MELM.Other_log) 
          fit.MELM.Other_log = as.data.frame(Effect(c("Year"), MELM.Other_log)) 
          
          ##### b) -Separate models per Country (In the Appendix) -----
          MELM.Other_log_NL <- lmer(log(median.other+1) ~ Year + (1|StationID), data=data_median_fg_NL) 
          summary(MELM.Other_log_NL)
          tab_model(MELM.Other_log_NL)
          
          MELM.Other_log_DE <- lmer(log(median.other+1) ~ Year + (1|StationID), data=data_median_fg_DE) 
          summary(MELM.Other_log_DE)
          tab_model(MELM.Other_log_DE) 

          

## 3.4 Plots -----

# Change order of stations
tot.biom.annual$StationID <- factor(tot.biom.annual$StationID,
c("MARSDND", "DOOVBWT","BOOMKDP", "TERSLG10", "DANTZGT", "ROTTMPT3", "HUIBGOT",  "BOCHTVWTM", "GROOTGND", 
      "Bork_W_1", "Nney_W_2", "JaBu_W_1", "WeMu_W_1"))
        
   ### 3.4.1 Median Biomass -----
        p.bm_log <- ggplot(tot.biom.annual, aes(x=Year, y=log(median.bm+1))) +
          geom_line(aes(Year,log(median.bm+1), col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
          scale_linetype_manual(values=c("dashed", "solid"))+
          geom_line(data = fit.MELM.bio_log, aes(y=fit), col="black", size=0.8) +
          geom_ribbon(data = fit.MELM.bio_log, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
          #scale_color_manual(values=c(Colors_NL_DE_8))+  
          scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
          labs(y = expression(paste("Biomass ", "[LN ","µgC ", L^{-1}, "]")))+
          xlab("year")+
          theme_classic()+
          theme(text = element_text(size = 18), 
                legend.position = "right") +
          theme(strip.text.x = element_text(size = 18))  +
          theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
          theme(panel.background = element_rect(fill = "white", colour = "black"))
        p.bm_log
        
        p.bm_lm_log <- p.bm_log + 
          annotate("text", x = 2001, y = 9, label = "p<0.001",
                   parse = TRUE, size = 4) 
        p.bm_lm_log
        
        
   ### 3.4.2 Median Biomass - Functional Groups -----
        # Diatoms
        p.Dia_log <- ggplot(data_median_fg, aes(x=Year, y=log(median.diatoms+1))) +
          geom_line(aes(Year,log(median.diatoms+1), col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
          scale_linetype_manual(values=c("dashed", "solid"))+
          geom_line(data = fit.MELM.Diatoms_log, aes(y=fit), col="black", size=0.8) +
          geom_ribbon(data = fit.MELM.Diatoms_log, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
          scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
          labs(y = expression(paste("Diatoms ", "[LN ","µgC ", L^{-1}, "]")))+
          xlab("year")+
          theme_classic()+
          theme(text = element_text(size = 18), 
                legend.position = "right") +
          theme(strip.text.x = element_text(size = 18))  +
          theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
          theme(panel.background = element_rect(fill = "white", colour = "black"))
        p.Dia_log
        p.Dia_lm_log <- p.Dia_log + 
          annotate("text", x = 2001, y = 9, label = "p<0.001",
                   parse = TRUE, size = 4) 
        p.Dia_lm_log
        
        #Dinoflagellates
        p.Dino_log <- ggplot(data_median_fg, aes(x=Year, y=log(median.dino+1))) +
          geom_line(aes(Year,log(median.dino+1), col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
          scale_linetype_manual(values=c("dashed", "solid"))+
          geom_line(data = fit.MELM.Dino_log, aes(y=fit), col="black", size=0.8) +
          geom_ribbon(data = fit.MELM.Dino_log, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
          scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
          labs(y = expression(paste("Dinoflagellates ", "[LN ","µgC ", L^{-1}, "]")))+
          xlab("year")+
          theme_classic()+
          theme(text = element_text(size = 18), 
                legend.position = "right") +
          theme(strip.text.x = element_text(size = 18))  +
          theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
          theme(panel.background = element_rect(fill = "white", colour = "black"))
        p.Dino_log
        p.Dino_lm_log <- p.Dino_log + 
          annotate("text", x = 2001, y = 3.5, label = "p==0.24",
                   parse = TRUE, size = 4) 
        p.Dino_lm_log

        #Flagellates
        p.Flag_log <- ggplot(data_median_fg, aes(x=Year, y=log(median.fla+1))) +
          geom_line(aes(Year,log(median.fla+1), col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
          scale_linetype_manual(values=c("dashed", "solid"))+
          geom_line(data = fit.MELM.Flagellates_log, aes(y=fit), col="black", size=0.8) +
          geom_ribbon(data = fit.MELM.Flagellates_log, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
          scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
          labs(y = expression(paste("Flagellates ", "[LN ","µgC ", L^{-1}, "]")))+
          xlab("year")+
          theme_classic()+
          theme(text = element_text(size = 18), 
                legend.position = "right") +
          theme(strip.text.x = element_text(size = 18))  +
          theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
          theme(panel.background = element_rect(fill = "white", colour = "black"))
        p.Flag_log
        p.Flag_lm_log <- p.Flag_log + 
          annotate("text", x = 2001, y = 5.5, label = "p<0.001",
                   parse = TRUE, size = 4) 
        p.Flag_lm_log        
        
        #Phaeocystis
        p.pha_log <- ggplot(data_median_fg, aes(x=Year, y=log(median.phae+1))) +
          geom_line(aes(Year,log(median.phae+1), col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
          scale_linetype_manual(values=c("dashed", "solid"))+
          geom_line(data = fit.MELM.Phaeocystis_log, aes(y=fit), col="black", size=0.8) +
          geom_ribbon(data = fit.MELM.Phaeocystis_log, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
          scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
          labs(y = expression(paste(italic("Phaeocystis "), "[LN ","µgC ", L^{-1}, "]")))+
          xlab("year")+
          theme_classic()+
          theme(text = element_text(size = 18), 
                legend.position = "right") +
          theme(strip.text.x = element_text(size = 18))  +
          theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
          theme(panel.background = element_rect(fill = "white", colour = "black"))
        p.pha_log
        p.pha_lm_log <- p.pha_log + 
          annotate("text", x = 2015, y = 5, label = "p==0.05",
                   parse = TRUE, size = 4) 
        p.pha_lm_log
        
        #Other
        p.other_log <- ggplot(data_median_fg, aes(x=Year, y=log(median_othercya+1))) +
          geom_line(aes(Year,log(median.other+1), col= StationID, linetype=Country), size=0.8,alpha=0.8,position=position_dodge(width=.3))+
          scale_linetype_manual(values=c("dashed", "solid"))+
          geom_line(data = fit.MELM.Other_log, aes(y=fit), col="black", size=0.8) +
          geom_ribbon(data = fit.MELM.Other_log, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2)+
          scale_color_manual(values=c(Colors_coastal_NL_DE_2))+
          labs(y = expression(paste("Other ", "[LN ","µgC ", L^{-1}, "]")))+
          xlab("year")+
          theme_classic()+
          theme(text = element_text(size = 18), 
                legend.position = "right") +
          theme(strip.text.x = element_text(size = 18))  +
          theme(legend.key.width = unit(1,"cm"), legend.box="vertical", legend.margin=margin())+
          theme(panel.background = element_rect(fill = "white", colour = "black"))
        p.other_log
        p.other_lm_log <- p.other_log + 
          annotate("text", x = 2000, y = 5.5, label = "p<0.001",
                   parse = TRUE, size = 4) 
        p.other_lm_log
        
        


# Sampling frequency (Appendix) ----------------------

#Add a sample date column: year - month - day
PPKT_WS$sampledate<-dmy(PPKT_WS$Date)
#add a year column 
PPKT_WS$year<-year(PPKT_WS$sampledate)
summary(PPKT_WS)
#add a month column
PPKT_WS$month<-month(PPKT_WS$sampledate, label = TRUE, abbr = TRUE)
#day of the month
PPKT_WS$day<-mday(PPKT_WS$sampledate)
length(unique(PPKT_WS$day)) #31

library(hydroTSM)
sampledate_vector <- PPKT_WS$sampledate
#default => "winter"= Dec, Jan, Feb; "spring"= Mar, Apr, May; "summer"=Jun, Jul, Aug; "autumn"= Sep, Oct, Nov
seasons <- time2season(sampledate_vector, out.fmt = "seasons", type="default")
#add another column called season where "winter"= Dec, Jan, Feb; "spring"= Mar, Apr, May; "summer"=Jun, Jul, Aug; "autumn"= Sep, Oct, Nov
PPKT_WS$season <- seasons
#correcting "autumn"
PPKT_WS$season <- gsub("autumm", "autumn", PPKT_WS$season)# sum of the abun per date, station

# sum of the abun per date, station
PPKT_WS_sum.ab <-ddply(PPKT_WS,.(Country, StationID, Date, month, season, year), 
                       colwise(sum,  .(abundance_l)), na.rm=T)

PPKT_WS_sum.ab$month <- factor(PPKT_WS_sum.ab$month, c("Dec", "Jan", "Feb","Mar","Apr","May", "Jun","Jul", "Aug","Sep", "Oct", "Nov"))
PPKT_WS_sum.ab$season <- factor(PPKT_WS_sum.ab$season, c("winter", "spring", "summer", "autumn"))
PPKT_WS_sum.ab$StationID <- factor(PPKT_WS_sum.ab$StationID,
                                   c("MARSDND", "DOOVBWT","BOOMKDP", "TERSLG10", "DANTZGT", "ROTTMPT3", "HUIBGOT",  "BOCHTVWTM", "GROOTGND", 
                                     "Bork_W_1", "Nney_W_2", "JaBu_W_1", "WeMu_W_1"))
library(ggthemes)
# Plot - sampling frequency
ggplot(PPKT_WS_sum.ab, aes(year, group=month)) +
  geom_histogram(aes(fill=season),
                 binwidth = 1,
                 col="black",
                 size=.1,
                 alpha=0.7) +  # change binwidth
  scale_fill_tableau()+
  labs(title="Phytoplankton sampling frequency",
       subtitle="Coastal stations",
       y="Count (samples)", x="Year") +
  facet_wrap(~StationID,scales="free_y", ncol=3)
