#### Code to partition Argentina Observed Mortality into rates for never, current for 6 categories and former smokers
#### Based on Rosenberg et al, Risk Anal 2012 methodology
#### Version 1 - Rafael Meza, CISNET LWG, Oct 2024
#### Uses as inputs: Argentina mortality rates by birth-cohort, prevalence by cohort files, and smoothed RRs by CPD
#### CISNET LWG / CEDES collaboration

#In the CPD analysis, never- and former-smoker mortality rates are intermediate quantities used to obtain internally consistent CPD-specific current-smoker mortality rates. 
#Only current-smoker CPD categories are interpreted from this analysis.

rm(list = ls())
setwd("")

library(RColorBrewer)
library(dplyr)
library(tidyr)

#################
#################

#####
#####
### BC1. Mortality rates by cohort year

### Read mortality rates by cohort year

BCMortF=read.csv("../Fem_ArgentinaCohortMortality_1900_2100.csv",header = TRUE)
BCMortM=read.csv("../Male_ArgentinaCohortMortality_1900_2100.csv",header = TRUE)

names(BCMortF)=names(BCMortM)=1900:2100

#####
#####
### BC2. Smoking Data

### Read smoking prevalence by cohort data (from Argentina´s)
curr<-read.table("../../../../../APC Analysis/Arg APC analysis/2025/output_sex/new_SHG_current.txt",sep=',',skip=4,header=TRUE)
neve<-read.table("../../../../../APC Analysis/Arg APC analysis/2025/output_sex/new_SHG_never.txt",sep=',',skip=4,header=TRUE)
form<-read.table("../../../../../APC Analysis/Arg APC analysis/2025/output_sex/new_SHG_former.txt",sep=',',skip=4,header=TRUE)

#clean it up
M.smkcurr<-curr[curr$Sex==0,c(-1,-2,-3)]
F.smkcurr<-curr[curr$Sex==1,c(-1,-2,-3)]
M.smkneve<-neve[neve$Sex==0,c(-1,-2,-3)]
F.smkneve<-neve[neve$Sex==1,c(-1,-2,-3)]
M.smkform<-form[form$Sex==0,c(-1,-2,-3)]
F.smkform<-form[form$Sex==1,c(-1,-2,-3)]

names(M.smkcurr)<-names(F.smkcurr)<-c(paste(seq(1881,2100)))
names(M.smkneve)<-names(F.smkneve)<-c(paste(seq(1881,2100)))
names(M.smkform)<-names(F.smkform)<-c(paste(seq(1881,2100)))

rm(curr,neve,form)

BCNeveF=F.smkneve[,20:220] ##from 1900 to 2100
BCCurrF=F.smkcurr[,20:220]
BCFormF=F.smkform[,20:220]

BCNeveM=M.smkneve[,20:220]
BCCurrM=M.smkcurr[,20:220]
BCFormM=M.smkform[,20:220]

### Read prevalence of CPD categories for current smokers
A=read.table('../../../../../APC Analysis/Arg APC analysis/2025/cpd analysis/output_cpd/new_SHG_cpd_race0.txt',sep=',',skip=7,header=FALSE,stringsAsFactors = F,na.strings = c("."," .         ","NA"))

currCAT6<-as.data.frame(data.matrix(A))
names(currCAT6)<-c("Race","Sex","StartYOB","EndYOB","Age","Cat1","Cat2","Cat3","Cat4","Cat5","Cat6")
currCAT6[currCAT6$Age%in%(0:7),6:11]=0. # add 0s to ages 0:7

currCAT6M=currCAT6 %>% filter(Sex==0)
currCAT6F=currCAT6 %>% filter(Sex==1)

### get prevalence for each category
currC1M=currCAT6M[,c(3,5,6)]%>%   ##c(3,5,6)=YOB, Age, Cat1
  pivot_wider(names_from = StartYOB,values_from = Cat1)
currC1F=currCAT6F[,c(3,5,6)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat1)
currC2M=currCAT6M[,c(3,5,7)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat2)
currC2F=currCAT6F[,c(3,5,7)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat2)
currC3M=currCAT6M[,c(3,5,8)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat3)
currC3F=currCAT6F[,c(3,5,8)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat3)
currC4M=currCAT6M[,c(3,5,9)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat4)
currC4F=currCAT6F[,c(3,5,9)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat4)
currC5M=currCAT6M[,c(3,5,10)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat5)
currC5F=currCAT6F[,c(3,5,10)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat5)
currC6M=currCAT6M[,c(3,5,11)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat6)
currC6F=currCAT6F[,c(3,5,11)]%>% 
  pivot_wider(names_from = StartYOB,values_from = Cat6)

currC1M=currC1M[,-c(1)];currC1F=currC1F[,-c(1)] ##delete 1900 column (first column)
currC2M=currC2M[,-c(1)];currC2F=currC2F[,-c(1)]
currC3M=currC3M[,-c(1)];currC3F=currC3F[,-c(1)]
currC4M=currC4M[,-c(1)];currC4F=currC4F[,-c(1)]
currC5M=currC5M[,-c(1)];currC5F=currC5F[,-c(1)]
currC6M=currC6M[,-c(1)];currC6F=currC6F[,-c(1)]

# multiply by current smoking prevalence to get prevalence for each category at the population level (ie, not just among current smokers)
currC1M=BCCurrM*currC1M;currC1F=BCCurrF*currC1F ## current smoking prevalence * CPD cat prevalence
currC2M=BCCurrM*currC2M;currC2F=BCCurrF*currC2F
currC3M=BCCurrM*currC3M;currC3F=BCCurrF*currC3F
currC4M=BCCurrM*currC4M;currC4F=BCCurrF*currC4F
currC5M=BCCurrM*currC5M;currC5F=BCCurrF*currC5F
currC6M=BCCurrM*currC6M;currC6F=BCCurrF*currC6F

### check passed- should be 0 or almost 0 due to computing operations rounding errors
sum(BCCurrM-currC1M-currC2M-currC3M-currC4M-currC5M-currC6M) #-0.000000006315704
sum(BCCurrF-currC1F-currC2F-currC3F-currC4F-currC5F-currC6F) #0.00000001346094

#####
### BC3. Smoking Mortality Relative Risks

### Current smokers
### Smoothed Current smoking RRs by CPD, rows ages 32:87; cols CPDs 0:60 -- smoothing with LOESS, es smoothing de CPS-II
RRCPDCatM <- read.csv("../../SmoothedRRAllCigM_upto60cpd.csv", header = FALSE) 
RRCPDCatF <- read.csv("../../SmoothedRRAllCigF_upto60cpd.csv", header = FALSE)

### Keep only CPDs for the mean CPD in each category (Argentina data)
#CPD≤5 (3), 5<CPD≤15 (11), 15<CPD≤25 (20), 25<CPD≤35 (30), 35<CPD≤45 (40), and 45<CPD (65; but only up to 60)

RRCPDCatM=RRCPDCatM[,c(4,12,21,31,41,61)]
RRCPDCatF=RRCPDCatF[,c(4,12,21,31,41,61)]

### Keep only ages 40 onwards
RRCPDCatM=RRCPDCatM[9:56,] ##age 40 to 87
RRCPDCatF=RRCPDCatF[9:56,]

### complete to age 99, keep value at age 
for (age in 88:99){
  RRCPDCatM=rbind(RRCPDCatM,RRCPDCatM[48,]) ##valor de edad 87 (la última disponible)
  RRCPDCatF=rbind(RRCPDCatF,RRCPDCatF[48,])
}

## complete from ages 0-39, assuming RR=1
for (age in 0:39){
  RRCPDCatM=rbind(c(1,1,1,1,1,1),RRCPDCatM)
  RRCPDCatF=rbind(c(1,1,1,1,1,1),RRCPDCatF)
}

### Former smokers based on the US for now - Thun et al, NEJM 2013 - CPSII values
### Former smoker rates are not use by SHG. Instead, the current smoker rates converge to the never smokers as people quit 
RRFSFThun=1.33
RRFSMThun=1.42

## assume RR=1 before age 40, constant by age (ages 0 to 99)
RRFSF=RRFSM=matrix(1,100,201) # 1900 to 2100 birth cohorts
RRFSF[41:100,]=RRFSFThun ##row 41 is age 40
RRFSM[41:100,]=RRFSMThun

#####
### Calculation of mortality rates for never, current by category and former smokers
### Following the Rosenberg et al Risk Anal 2012  approach

BCMuNSF=as.matrix(BCMortF/(BCNeveF+currC1F*RRCPDCatF[,1]+currC2F*RRCPDCatF[,2]+currC3F*RRCPDCatF[,3]+currC4F*RRCPDCatF[,4]+currC5F*RRCPDCatF[,5]+currC6F*RRCPDCatF[,6]+BCFormF*RRFSF))
BCMuNSM=as.matrix(BCMortM/(BCNeveM+currC1M*RRCPDCatM[,1]+currC2M*RRCPDCatM[,2]+currC3M*RRCPDCatM[,3]+currC4M*RRCPDCatM[,4]+currC5M*RRCPDCatM[,5]+currC6M*RRCPDCatM[,6]+BCFormM*RRFSM))

BCMuCS1F=pmin(as.matrix(RRCPDCatF[,1]*BCMuNSF),0.999)
BCMuCS1M=pmin(as.matrix(RRCPDCatM[,1]*BCMuNSM),0.999)
BCMuCS2F=pmin(as.matrix(RRCPDCatF[,2]*BCMuNSF),0.999)
BCMuCS2M=pmin(as.matrix(RRCPDCatM[,2]*BCMuNSM),0.999)
BCMuCS3F=pmin(as.matrix(RRCPDCatF[,3]*BCMuNSF),0.999)
BCMuCS3M=pmin(as.matrix(RRCPDCatM[,3]*BCMuNSM),0.999)
BCMuCS4F=pmin(as.matrix(RRCPDCatF[,4]*BCMuNSF),0.999)
BCMuCS4M=pmin(as.matrix(RRCPDCatM[,4]*BCMuNSM),0.999)
BCMuCS5F=pmin(as.matrix(RRCPDCatF[,5]*BCMuNSF),0.999)
BCMuCS5M=pmin(as.matrix(RRCPDCatM[,5]*BCMuNSM),0.999)
BCMuCS6F=pmin(as.matrix(RRCPDCatF[,6]*BCMuNSF),0.999)
BCMuCS6M=pmin(as.matrix(RRCPDCatM[,6]*BCMuNSM),0.999)

BCMuFSF=pmin(as.matrix(RRFSF*BCMuNSF),0.999)
BCMuFSM=pmin(as.matrix(RRFSM*BCMuNSM),0.999)

### Export mortality files as csv
write.csv(BCMuNSF,"MuNSF_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS1F,"MuCS1F_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS2F,"MuCS2F_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS3F,"MuCS3F_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS4F,"MuCS4F_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS5F,"MuCS5F_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS6F,"MuCS6F_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuFSF,"MuFSF_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuNSM,"MuNSM_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS1M,"MuCS1M_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS2M,"MuCS2M_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS3M,"MuCS3M_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS4M,"MuCS4M_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS5M,"MuCS5M_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCS6M,"MuCS6M_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuFSM,"MuFSM_BirthYr_1900_2100.csv",row.names = FALSE)

## ============================================================
## Export long-format mortality table for CURRENT smokers only
## (6 CPD categories), with both rate and probability
## per (calendar year) = coh (birth cohort) + age
## ============================================================

ages    <- 0:99
cohorts <- 1900:2100

mat_to_long_cs <- function(mu_mat, sex_label, status_label, ages, cohorts) {
  mu_mat <- as.matrix(mu_mat)
  storage.mode(mu_mat) <- "numeric"
  
  # Ensure cohort column names exist (and are 1900:2100)
  if (is.null(colnames(mu_mat)) || length(colnames(mu_mat)) != length(cohorts)) {
    colnames(mu_mat) <- as.character(cohorts)
  }
  
  # Wide -> long
  df <- as.data.frame(mu_mat, check.names = FALSE)
  df$age <- ages
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = -age,
    names_to = "coh",
    values_to = "Mortality_rate"
  )
  
  # If any cohort names look like "X1900", strip the leading "X"
  df_long$coh <- as.integer(sub("^X", "", df_long$coh))
  
  # per = coh + age
  df_long$per <- df_long$coh + df_long$age
  
  # Convert rate <-> probability using 1-year time interval
  df_long$Mortality_probability <- 1 - exp(-df_long$Mortality_rate * 1)
  
  # Recompute rate from probability (per your required formula)
  df_long$Mortality_rate <- -log(1 - df_long$Mortality_probability) / 1
  
  df_long$SEX <- sex_label
  df_long$Smoking_Status <- status_label
  
  df_long %>%
    select(
      age,
      per,
      coh,
      SEX,
      Mortality_probability,
      Mortality_rate,
      Smoking_Status
    )
}

# Build the full dataset (Current smokers only: CS1–CS6)
df_mortality_cs <- bind_rows(
  mat_to_long_cs(BCMuCS1F, "Female", "CS1", ages, cohorts),
  mat_to_long_cs(BCMuCS2F, "Female", "CS2", ages, cohorts),
  mat_to_long_cs(BCMuCS3F, "Female", "CS3", ages, cohorts),
  mat_to_long_cs(BCMuCS4F, "Female", "CS4", ages, cohorts),
  mat_to_long_cs(BCMuCS5F, "Female", "CS5", ages, cohorts),
  mat_to_long_cs(BCMuCS6F, "Female", "CS6", ages, cohorts),
  
  mat_to_long_cs(BCMuCS1M, "Male", "CS1", ages, cohorts),
  mat_to_long_cs(BCMuCS2M, "Male", "CS2", ages, cohorts),
  mat_to_long_cs(BCMuCS3M, "Male", "CS3", ages, cohorts),
  mat_to_long_cs(BCMuCS4M, "Male", "CS4", ages, cohorts),
  mat_to_long_cs(BCMuCS5M, "Male", "CS5", ages, cohorts),
  mat_to_long_cs(BCMuCS6M, "Male", "CS6", ages, cohorts)
)

# Write CSV
write.csv(
  df_mortality_cs,
  file = "Argentina_Mortality_Age_Per_Coh_Sex_6CPD_Rate_Prob.csv",
  row.names = FALSE
)
