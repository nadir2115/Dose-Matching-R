9#Matching problem on evaluating effects of dose- Nadir Nibras

# clear workspace variables
rm(list = ls()); 
# clear window (same as ctrl+L. )
cat("\014")   
# close all plots
graphics.off() 

library(tidyverse)
library(Matching)
library(dummy)
library(dummies)
library(MatchIt)
library(RItools)
library(rbounds)

# set directory
setwd("C:/Users/nadir/Desktop/Matching/BKN 599 Advanced Data Analysis")

# variables
threshold= 27
replace= 1

# Read files
evalData=read.csv("00042_evalDataset_notext.csv")
patChar= read.csv("PatientCharacteristics_withDose.csv")

# Changing column name for merging
names(patChar)[6]<-paste("pid")  

# Remove sparse features/features with no use
evalData= evalData[-c(1:4, 7:203, 205:245, 247:276, 278:366, 369:376)]  
# Note: Leaving out FAS scores because many (18) missed values for visit 1 stats

patChar= patChar[c(6,9,12,18,19,27,29,30)]

# Merge data
dataMerged= merge(evalData,patChar)

# Feature manipulation
dataMerged$severity= ifelse(dataMerged$severity=="H", 1, 0)
dataMerged$finished_college=ifelse(dataMerged$DEMOedu>4,1,0)
dataMerged$DEMOedu=NULL

# Replace -99 with missing values
dataMerged[dataMerged == -99] <- NA


# Changing table to numeric
for (i in 2:ncol(dataMerged)){
  dataMerged[c(i)]= as.numeric(unlist(dataMerged[c(i)]))
}

# Extract visit datas
dataVisit1= subset(dataMerged, visitnum==1)                        
dataVisit2= subset(dataMerged, visitnum==2)                        

# Removing visit2 data with missing values for NIHtotal
dataVisit2= dataVisit2[complete.cases(dataVisit2[,5]),]

# Removing visit1 data with no available visit 2 features
temp= data.frame(dataVisit2$pid); 
names(temp)[1]<-paste("pid");   #renaming column name to "pid"
dataVisit1= merge(temp,dataVisit1)
rm(temp)

# Check number of missing values
sum(is.na(dataVisit1))
sum(is.na(dataVisit2))


# Impute missing values- is.na(dataVisit1[,i]) returns TRUE only for row index of column with NA
for(i in 2:ncol(dataVisit1)){
  dataVisit1[is.na(dataVisit1[,i]), i] <- median(dataVisit1[,i], na.rm = TRUE)
}  
for(i in 2:ncol(dataVisit2)){
  dataVisit2[is.na(dataVisit2[,i]), i] <- median(dataVisit2[,i], na.rm = TRUE)
}

# Check number of missing values AGAIN
sum(is.na(dataVisit1))
sum(is.na(dataVisit2))


# Creating new dataset for future calculations
data12= dataVisit1

# creating binary variable indicating whether treatment was performed
data12$treat=ifelse(dataVisit1$dose_hours<=threshold, 1, 0)

# # normalize by subtracting mean and dividing by std
# for(i in 2:(ncol(data12)-1)){
#   data12[,i] <- (data12[,i]-mean(data12[,i]))/(max(data12[,i])-min(data12[,i]))
# }

data12$sishand_change= dataVisit2$SIS_hand-dataVisit1$SIS_hand
data12$fugm_change= dataVisit2$ufugm-dataVisit1$ufugm
data12$logwmftMA_change= dataVisit2$log_mean_time_MA_PA-dataVisit1$log_mean_time_MA_PA
data12$logwmftLA_change= dataVisit2$log_mean_time_LA_PA-dataVisit1$log_mean_time_LA_PA
data12$FAS_change= dataVisit2$FAS_score_PA-dataVisit1$FAS_score_PA


highdosegroup= subset(data12, treat==0)
lowdosegroup= subset(data12, treat==1)

# highdose= subset(dataVisit1, treat==0)
# lowdose= subset(dataVisit1, treat==1)
# highdoseindex= which(dataVisit1$dose_hours>threshold)
# lowdoseindex= which(dataVisit1$dose_hours<=threshold)
# 
# Euclidean distance matching by Matlab rules ---------------------------------------------
# 
# euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
# lowdosetemp = lowdose[c(3:ncol(lowdose)-1)];
# highdosetemp = highdose[c(3:ncol(highdose)-1)];
# matchedindex_e=data.frame(matrix(0, nrow(lowdosetemp), 2))
# for (i in 1:nrow(lowdosetemp)){
# dist=1000;
# for (j in 1:nrow(highdosetemp)){
# if (euc.dist(lowdosetemp[i,],highdosetemp[j,])<dist){
# dist= euc.dist(lowdosetemp[i,],highdosetemp[j,])
# matchedindex_e[i,]= c(j,dist)
# }
# }
# }
# 
# matchedhighdose= highdose[matchedindex_e$X1,]


ggplot(data12,aes(fugm_change))+
  geom_histogram(fill="white", col="black", binwidth = 1)


# Youtube tutorial points -------------------------------------------------

# creates a marrix of columns/variables
X= cbind(data12$SIS_hand, data12$ufugm,data12$CSgender,data12$NIHtot, 
         data12$log_mean_time_MA_PA, data12$log_mean_time_LA_PA,
         data12$severity, data12$onset_to_rand,data12$age_at_rand,
         data12$concordance, data12$finished_college)

# Output variable that we will evaluate
y=data12$log_mean_time_MA_PA

# Generalized linear models are fit using the glm( ) function.

# Logistic Regression
# where treat is a binary factor and 
# x is a matrix of all predictor variables  
# Family objects provide a convenient way to specify the details of the models used by functions such as glm
# family= binomial() performs logistic regression
# treat ~ X means treat is the response of all other variables
model <- glm(treat ~ X, 
             family = binomial(), data = data12)

summary(model) #"Estimate" values indicate more likely to be in treated (LD) group if positive

# Match implements a variety of algorithms for multivariate matching including 
# propensity score, Mahalanobis and inverse variance matching. 
# The function is intended to be used in conjunction with the 
# MatchBalance function which determines the extent to which 
# Match has been able to achieve covariate balance


# average treated effect on the treated
rr1 =Match(Y=y, Tr= data12$treat, X= model$fitted.values, estimand = "ATT", 
           M=1, ties= TRUE, replace= TRUE)
summary(rr1) #Estimate is effect on Y if individuals are treated

rr2 =Match(Y=y, Tr= data12$treat, X= model$fitted.values, estimand = "ATE", 
           M=1, ties= TRUE, replace= TRUE)
summary(rr2)

# Assessing balance
MatchBalance(treat~ SIS_hand + ufugm + NIHtot + log_mean_time_MA_PA +
               log_mean_time_LA_PA +CSgender+ concordance + severity + 
               onset_to_rand + finished_college + age_at_rand,
             match.out= rr1, nboots=0,data= data12)
# MatchBalance(treat~X, match.out= rr1, nboots=0,data= data12)

qqplot(data12$SIS_hand[rr1$index.control],data12$SIS_hand[rr1$index.treated])
abline(coef=c(0,1),col=2)

qqplot(data12$wmft_mean_time_MA_PA[rr1$index.control],data12$wmft_mean_time_MA_PA[rr1$index.treated])
abline(coef=c(0,1),col=2)


