#Matching problem on evaluating effects of dose- Nadir Nibras

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
# library(optmatch)
library(MatchIt)
library(RItools)
library(rbounds)

# set directory
setwd("C:/Users/nadir/Desktop/Matching/BKN 599 Advanced Data Analysis")

# Read files
evalData=read.csv("00042_evalDataset_notext.csv")
patChar= read.csv("PatientCharacteristics_withDose.csv")

# variables
threshold= 27
replace= 1
onehotencode= 0

names(patChar)[6]<-paste("pid")  # Changing column name for merging

# Remove sparse features/features with no use
evalData=evalData[-c(1:4, 7:127, 129:203, 205:245, 247:259, 261:276, 278:364, 369:375)]  
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
dataVisit1= merge(temp,dataVisit1)
rm(temp)

# Replace -99 with missing values
dataVisit1[dataVisit1 == -99] <- NA
dataVisit2[dataVisit2 == -99] <- NA

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
  dataVisit1=dataVisit1[-c(2,15:23,25,27)]
  dataVisit2=dataVisit2[-c(2,15:23,25,27)]
}

data12= dataVisit1

data12$treat=ifelse(dataVisit1$dose_hours<=threshold, 1, 0)

# # normalize by subtracting mean and dividing by std
# for(i in 2:(ncol(data12)-1)){
#   data12[,i] <- (data12[,i]-mean(data12[,i]))/sd(data12[,i])
# }  

data12$eqindex_change= dataVisit2$EQ_index-dataVisit1$EQ_index
data12$sishand_change= dataVisit2$SIS_hand-dataVisit1$SIS_hand
data12$fugm_change= dataVisit2$ufugm-dataVisit1$ufugm
data12$rnli_change= dataVisit2$RNLIadj-dataVisit1$RNLIadj
data12$wmftMA_change= dataVisit2$wmft_mean_time_MA_PA-dataVisit1$wmft_mean_time_MA_PA
data12$wmftLA_change= dataVisit2$wmft_mean_time_LA_PA-dataVisit1$wmft_mean_time_LA_PA
data12$FAS_change= dataVisit2$FAS_score_PA-dataVisit1$FAS_score_PA

# highdose= subset(dataVisit1, treat==0)
# lowdose= subset(dataVisit1, treat==1)
# 
# highdoseindex= which(dataVisit1$dose_hours>threshold)
# lowdoseindex= which(dataVisit1$dose_hours<=threshold)
# 
# Euclidean distance matching by Matlab rules ---------------------------------------------
# 
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

X= cbind(data12$EQ_index,data12$SIS_hand,data12$ufugm,data12$RNLIadj,data12$NIHtot,
         data12$wmft_mean_time_MA_PA,data12$wmft_mean_time_LA_PA,data12$log_mean_time_MA_PA,
         data12$log_mean_time_LA_PA,data12$FAS_score_PA,data12$severity,data12$onset_to_rand, 
         data12$age_at_rand)

y=data12$sishand_change

model <- glm(treat ~ X, 
             family = binomial(), data = data12)

# EQ_index+SIS_hand+ufugm+RNLIadj+NIHtot+
#   wmft_mean_time_MA_PA+wmft_mean_time_LA_PA+log_mean_time_MA_PA+
#   log_mean_time_LA_PA+FAS_score_PA+severity+onset_to_rand +
#   age_at_rand

summary(model) #"estimate" value indicates more likely to be in treated (LD) group if positive

# average treated effect on the treated
rr1 =Match(Y=y, Tr= data12$treat, X= model$fitted.values, estimand = "ATT", 
           M=1, ties= TRUE, replace= TRUE)
summary(rr1) #Estimate is effect on Y if individuals are treated

rr2 =Match(Y=y, Tr= data12$treat, X= model$fitted.values, estimand = "ATE", 
           M=1, ties= TRUE, replace= TRUE)
summary(rr2)

# Assessing balance
MatchBalance(treat~ EQ_index + SIS_hand + ufugm + RNLIadj + NIHtot +
               wmft_mean_time_MA_PA + wmft_mean_time_LA_PA + log_mean_time_MA_PA +
               log_mean_time_LA_PA + FAS_score_PA + severity + onset_to_rand +
               age_at_rand,
             match.out= rr1, nboots=0,data= data12)
# MatchBalance(treat~X, match.out= rr1, nboots=0,data= data12)

qqplot(data12$SIS_hand[rr1$index.control],data12$SIS_hand[rr1$index.treated])
abline(coef=c(0,1),col=2)

qqplot(data12$wmft_mean_time_MA_PA[rr1$index.control],data12$wmft_mean_time_MA_PA[rr1$index.treated])
abline(coef=c(0,1),col=2)


# Genetic matching --------------------------------------------------------
gen1 = GenMatch(Tr= data12$treat, X=X, BalanceMatrix = X, pop.size=100)
mgen1=Match(Y=y, Tr=data12$treat, X=X, Weight.matrix = gen1)
summary(mgen1)

# Assessing balance
MatchBalance(treat~ EQ_index + SIS_hand + ufugm + RNLIadj + NIHtot +
               wmft_mean_time_MA_PA + wmft_mean_time_LA_PA + log_mean_time_MA_PA +
               log_mean_time_LA_PA + FAS_score_PA + severity + onset_to_rand +
               age_at_rand,
             data= data12, match.out= mgen1, nboots=0)

qqplot(data12$SIS_hand[mgen1$index.control],data12$SIS_hand[mgen1$index.treated])
abline(coef=c(0,1),col=2)

qqplot(data12$wmft_mean_time_MA_PA[mgen1$index.control],data12$wmft_mean_time_MA_PA[mgen1$index.treated])
abline(coef=c(0,1),col=2)

