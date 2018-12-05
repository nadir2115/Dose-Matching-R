#Matching problem on evaluating effects of dose- Nadir Nibras

# clear workspace variables
rm(list = ls()); 
# clear window (same as ctrl+L. )
cat("\014")   
# close all plots
graphics.off() 



library(Matching)

# set directory
setwd("C:/Users/nadir/Desktop/Matching/BKN 599 Advanced Data Analysis")

# Read files
evalData=read.csv("00042_evalDataset_notext.csv")
patChar= read.csv("PatientCharacteristics_withDose.csv")

# variables
thresh= 27
replace= 1
onehotencode= 0


names(patChar)[6]<-paste("pid")  # Changing column name for merging

# Remove sparse features/features with no use
evalData=evalData[-c(1:4, 7:119, 120:127, 129:203, 205:245, 248:259, 261:276, 278:364, 369:375)]  
patChar=patChar[-c(1:5, 7:8, 10:11, 13, 23:26, 28,31, 33)]

# Merge data
dataMerged= merge(evalData,patChar)

# Extract visit datas
dataVisit1= subset(dataMerged, visitnum==1)                        
dataVisit2= subset(dataMerged, visitnum==2)                        

# Removing visit1 data with no available visit 2 features
temp= data.frame(dataVisit2$pid); 
names(temp)[1]<-paste("pid");
dataVisit1= merge(dataVisit1, temp)
rm(temp)
