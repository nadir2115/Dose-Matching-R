#Matching problem on evaluating effects of dose- Nadir Nibras

# clear workspace variables
rm(list = ls()); 
# clear window (same as ctrl+L. )
cat("\014")   
# close all plots
graphics.off() 

library(tidyverse)
library(caret)
library(Matching)
library(dummy)
library(dummies)
library(optmatch)
library(MatchIt)

# set directory
setwd("C:/Users/nadir/Desktop/Matching/BKN 599 Advanced Data Analysis")

# Read files
evalData=read.csv("00042_evalDataset_notext.csv")
patChar= read.csv("PatientCharacteristics_withDose.csv")

# variables
threshold= 27
replace= 1
onehotencode= 1


names(patChar)[6]<-paste("pid")  # Changing column name for merging

# Remove sparse features/features with no use
evalData=evalData[-c(1:4, 7:119, 120:127, 129:203, 205:245, 248:259, 261:276, 278:364, 369:375)]  
patChar=patChar[-c(1:5, 7:8, 10:11, 13, 23:26, 28,31, 33)]

# Merge data
dataMerged= merge(evalData,patChar)
dataMerged$severity= ifelse(dataMerged$severity=="H", 1, 0)
names(dataMerged)[18]<-("Hispanic")  # Changing column name for merging



# Changing table to numeric
for (i in 2:ncol(dataMerged)){
  dataMerged[c(i)]= as.numeric(unlist(dataMerged[c(i)]))
}

# Extract visit datas
dataVisit1= subset(dataMerged, visitnum==1)                        
dataVisit2= subset(dataMerged, visitnum==2)                        

# Removing visit1 data with no available visit 2 features
temp= data.frame(dataVisit2$pid); 
names(temp)[1]<-paste("pid");
dataVisit1= merge(dataVisit1, temp)
rm(temp)

sum(is.na(dataVisit1))



# Impute missing values- is.na(dataVisit1[,i]) returns TRUE only for row index of column with NA
for(i in 2:ncol(dataVisit1)){
          dataVisit1[is.na(dataVisit1[,i]), i] <- median(dataVisit1[,i], na.rm = TRUE)
}  


if (onehotencode==1){
  dataVisit1$visitnum=NULL
  dataVisit1 = dummy.data.frame(dataVisit1, names=c("center"), sep="_") 
  dataVisit1 = dummy.data.frame(dataVisit1, names=c("DEMORace"), sep="_") 
  dataVisit1$Hispanic= ifelse(dataVisit1$Hispanic==1, 1, 0)
  dataVisit1$DEMORace_1=ifelse(dataVisit1$DEMORace_1==1|dataVisit1$DEMORace_2==1|dataVisit1$DEMORace_3==1,1,0)
  dataVisit1$DEMORace_2=NULL
  dataVisit1$DEMORace_3=NULL
  dataVisit1$DEMORace_7=NULL
  dataVisit1$stroke_wo_hemmo=ifelse(dataVisit1$CSstrktype==1,1,0)
  dataVisit1$CSstrktype=NULL
  dataVisit1 = dummy.data.frame(dataVisit1, names=c("CSstrkloc"), sep="_") 
  dataVisit1$CSstrkloc_3=NULL
  dataVisit1$finishedcollege=ifelse(dataVisit1$DEMOedu>4,1,0)
  dataVisit1$DEMOedu=NULL
}
if (onehotencode==0){
  dataVisit1=dataVisit1[-c(2,16:24,26,28)]
}

highdose= subset(dataVisit1, dose_hours>threshold)
lowdose= subset(dataVisit1, dose_hours<=threshold)

dataVisit1$treat=ifelse(dataVisit1$dose_hours<=threshold, 1, 0)

