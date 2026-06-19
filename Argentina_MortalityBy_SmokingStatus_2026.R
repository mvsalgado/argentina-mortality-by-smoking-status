#### Code to partition Argentina Observed Mortality into rates for never, current and former smokers
#### Based on Rosenberg et al, Risk Anal 2012 methodology
#### Version 1 - Rafael Meza, CISNET LWG, Dec 2023
#### Includes calculation of mortality rates by birth-cohort
#### CISNET LWG / CEDES collaboration

library(readr)
library(readxl)
library(tidyverse)
library(Epi)
library(RColorBrewer)

rm(list = ls())

library(here)
here::i_am("ArgentinaMortalityBy_SmokingStatus_MVS_2026.R")

#####
### 1. Mortality rates by calendar year

### Read mortality rates by calendar year

MortF=read.csv("Fem_ArgentinaMortality_1994_2050.csv",header = TRUE)
MortM=read.csv("Male_ArgentinaMortality_1994_2050.csv",header = TRUE)

names(MortF)=names(MortM)=1994:2050

#####
### 2. Smoking Data

### Read smoking prevalence by cohort data (from Argentina´s SHG)

curr<-read.table("../../../../APC Analysis/Arg APC analysis/2025/output_sex/new_SHG_current.txt",sep=',',skip=4,header=TRUE)
neve<-read.table("../../../../APC Analysis/Arg APC analysis/2025/output_sex/new_SHG_never.txt",sep=',',skip=4,header=TRUE)
form<-read.table("../../../../APC Analysis/Arg APC analysis/2025/output_sex/new_SHG_former.txt",sep=',',skip=4,header=TRUE)

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

### Flip smoking prevalence to calendar year

CurrF=matrix(0,100,57) ##57 is the number of years between 1994 and 2050 (2050-1994+1)
FormF=matrix(0,100,57)
NeveF=matrix(0,100,57)
CurrM=matrix(0,100,57)
FormM=matrix(0,100,57)
NeveM=matrix(0,100,57)
colnames(CurrF)=colnames(FormF)=colnames(NeveF)=1994:2050
colnames(CurrM)=colnames(FormM)=colnames(NeveM)=1994:2050

for (year in 1994:2050){
  for (age in 0:99){
    i=age+1  ## age index
    j=year-1993 ## year index
    k=year-age-1880 ## cohort index
    CurrF[i,j]= F.smkcurr[i,k]## read age value for corresponding cohort prevalence. 
    FormF[i,j]= F.smkform[i,k]
    NeveF[i,j]= F.smkneve[i,k]
    CurrM[i,j]= M.smkcurr[i,k]
    FormM[i,j]= M.smkform[i,k]
    NeveM[i,j]= M.smkneve[i,k]
  }
}

#####
### 3. Smoking Mortality Relative Risks

### Based on the US - Thun et al, NEJM 2013 - CPSII values

RRCSFThun=2.08 ##RRCSF = RR current smoking female
RRFSFThun=1.33 ##RRFSF = RR former smoking female
RRCSMThun=2.33
RRFSMThun=1.42
  
RRCSF=RRFSF=RRCSM=RRFSM=matrix(1,100,57) # 57 years

## assume 1 before age 40, constant by age 
RRCSF[41:100,]=RRCSFThun ##row 41 is age 40 (row 1 is age 0)
RRFSF[41:100,]=RRFSFThun
RRCSM[41:100,]=RRCSMThun
RRFSM[41:100,]=RRFSMThun

#colnames(RRCSF)=colnames(RRFSF)=colnames(RRCSM)=colnames(RRFSM)=1994:2050 ##added by MVS

#####
### Calculation of mortality rates for never, current and former smokers
### Following the Rosenberg et al Risk Anal 2012  approach

MuNSF=as.matrix(MortF/(NeveF+CurrF*RRCSF+FormF*RRFSF))
MuNSM=as.matrix(MortM/(NeveM+CurrM*RRCSM+FormM*RRFSM))
MuCSF=pmin(as.matrix(RRCSF*MuNSF),0.999) ## pmin to ensure it is never bigger than 0.999
MuFSF=pmin(as.matrix(RRFSF*MuNSF),0.999)
MuCSM=pmin(as.matrix(RRCSM*MuNSM),0.999) 
MuFSM=pmin(as.matrix(RRFSM*MuNSM),0.999)


# Calculate the maximum mortality rate across all datasets
max_mortality <- max(c(MuCSF[,yr-1994+1], MuFSF[,yr-1994+1], MuNSF[,yr-1994+1], 
                       MuCSM[,yr-1994+1], MuFSM[,yr-1994+1], MuNSM[,yr-1994+1]), na.rm = TRUE)

### Export mortality files as csv
write.csv(MuNSF,"MuNSF_CalYr_1994_2050.csv",row.names = FALSE)
write.csv(MuCSF,"MuCSF_CalYr_1994_2050.csv",row.names = FALSE)
write.csv(MuFSF,"MuFSF_CalYr_1994_2050.csv",row.names = FALSE)
write.csv(MuNSM,"MuNSM_CalYr_1994_2050.csv",row.names = FALSE)
write.csv(MuCSM,"MuCSM_CalYr_1994_2050.csv",row.names = FALSE)
write.csv(MuFSM,"MuFSM_CalYr_1994_2050.csv",row.names = FALSE)

### Life expectancy

source("../LE_Calc.R")

LENF=LEcalc(MortF)
LENM=LEcalc(MortM)
LENSF=LEcalc(MuNSF)
LECSF=LEcalc(MuCSF)
LEFSF=LEcalc(MuFSF)
LENSM=LEcalc(MuNSM)
LECSM=LEcalc(MuCSM)
LEFSM=LEcalc(MuFSM)


LE40=matrix(0,2,8)
LE40[1,]=c(LENF[41,1],LENSF[41,1],LECSF[41,1],LEFSF[41,1],LENM[41,1],LENSM[41,1],LECSM[41,1],LEFSM[41,1])
LE40[2,]=c(LENF[41,25],LENSF[41,25],LECSF[41,25],LEFSF[41,25],LENM[41,25],LENSM[41,25],LECSM[41,25],LEFSM[41,25])

colnames(LE40)=c("Fem","NS Fem","CS Fem","FS Fem","Male","NS Male","CS Male","FS Male")
rownames(LE40)=c(1994,2018)

write.csv(LE40,"Arg_LE40_94_2018.csv")


#################
#################


#####
### BC1. Mortality rates by cohort year

### Read mortality rates by cohort year

BCMortF=read.csv("Fem_ArgentinaCohortMortality_1900_2100.csv",header = TRUE)
BCMortM=read.csv("Male_ArgentinaCohortMortality_1900_2100.csv",header = TRUE)

names(BCMortF)=names(BCMortM)=1900:2100

#####
### BC2. Smoking Data - already read
### we will use M.smkcurr, F.smkcurr, Msmkneve, Fsmkneve, Msmkform, Fsmkform

BCNeveF=F.smkneve[,20:220]
BCCurrF=F.smkcurr[,20:220]
BCFormF=F.smkform[,20:220]

BCNeveM=M.smkneve[,20:220]
BCCurrM=M.smkcurr[,20:220]
BCFormM=M.smkform[,20:220]

#####
### BC3. Smoking Mortality Relative Risks

RRCSFThun=2.08
RRFSFThun=1.33
RRCSMThun=2.33
RRFSMThun=1.42

RRCSF=RRFSF=RRCSM=RRFSM=matrix(1,100,201) # 57 years

## assume 1 before age 40, constant by age
RRCSF[41:100,]=RRCSFThun
RRFSF[41:100,]=RRFSFThun
RRCSM[41:100,]=RRCSMThun
RRFSM[41:100,]=RRFSMThun

#####
### Calculation of mortality rates for never, current and former smokers
### Following the Rosenberg et al Risk Anal 2012  approach

BCMuNSF=as.matrix(BCMortF/(BCNeveF+BCCurrF*RRCSF+BCFormF*RRFSF))
BCMuNSM=as.matrix(BCMortM/(BCNeveM+BCCurrM*RRCSM+BCFormM*RRFSM))
BCMuCSF=pmin(as.matrix(RRCSF*BCMuNSF),0.999)
BCMuFSF=pmin(as.matrix(RRFSF*BCMuNSF),0.999)
BCMuCSM=pmin(as.matrix(RRCSM*BCMuNSM),0.999)
BCMuFSM=pmin(as.matrix(RRFSM*BCMuNSM),0.999)

### Export mortality files as csv
write.csv(BCMuNSF,"MuNSF_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCSF,"MuCSF_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuFSF,"MuFSF_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuNSM,"MuNSM_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuCSM,"MuCSM_BirthYr_1900_2100.csv",row.names = FALSE)
write.csv(BCMuFSM,"MuFSM_BirthYr_1900_2100.csv",row.names = FALSE)

### Life expectancy

BCLENF=LEcalc(BCMortF)
BCLENM=LEcalc(BCMortM)
BCLENSF=LEcalc(BCMuNSF)
BCLECSF=LEcalc(BCMuCSF)
BCLEFSF=LEcalc(BCMuFSF)
BCLENSM=LEcalc(BCMuNSM)
BCLECSM=LEcalc(BCMuCSM)
BCLEFSM=LEcalc(BCMuFSM)



LE40_BC <- matrix(0, 8, 201)

for (c in 1:201) {  # cohorts 1900..2100
  LE40_BC[, c] <- c(
    BCLENF[41, c],  BCLENSF[41, c],  BCLECSF[41, c],  BCLEFSF[41, c],
    BCLENM[41, c],  BCLENSM[41, c],  BCLECSM[41, c],  BCLEFSM[41, c]
  )
}

rownames(LE40_BC) <- c("Fem","NS Fem","CS Fem","FS Fem","Male","NS Male","CS Male","FS Male")
colnames(LE40_BC) <- 1900:2100

write.csv(LE40_BC, "Arg_LE40_ByBirthCohort_1900_2100.csv", row.names = TRUE)

###############################################################
## =========================================================
## Export long-format mortality table (Age x Cohort matrices)
## Columns: age, per, coh, SEX, Mortality_probability, Mortality_rate, Smoking_Status
## Conversions (time = 1 year):
##   prob = 1 - exp(-rate * time)
##   rate = -log(1 - prob) / time
## =========================================================

# If column names got prefixed with "X" (e.g., X1900), remove it safely
strip_X <- function(x) as.integer(sub("^X", "", x))

# Convert an Age (rows) x Cohort (cols) mortality RATE matrix into long format
mat_agecoh_to_long <- function(mat_rate, sex_label, status_label, time = 1) {
  mat_rate <- as.matrix(mat_rate)
  
  if (is.null(colnames(mat_rate))) {
    stop("Input matrix must have cohort years as column names (e.g., 1900:2100).")
  }
  
  coh_years <- strip_X(colnames(mat_rate))
  ages <- 0:(nrow(mat_rate) - 1)
  
  df <- as.data.frame(mat_rate, check.names = FALSE)
  df$age <- ages
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = -age,
    names_to = "coh",
    values_to = "Mortality_rate"
  )
  
  df_long$coh <- strip_X(df_long$coh)
  df_long$per <- df_long$coh + df_long$age
  
  # Convert RATE -> PROBABILITY (time = 1 year)
  df_long$Mortality_probability <- 1 - exp(-df_long$Mortality_rate * time)
  
  # Numerical guard (optional)
  df_long$Mortality_probability <- pmin(pmax(df_long$Mortality_probability, 0), 1)
  
  df_long$SEX <- sex_label
  df_long$Smoking_Status <- status_label
  
  df_long <- df_long[, c(
    "age", "per", "coh", "SEX",
    "Mortality_probability", "Mortality_rate",
    "Smoking_Status"
  )]
  
  df_long
}

# ---- Build the full dataset from your cohort-based RATE matrices ----
df.Mortality <- dplyr::bind_rows(
  mat_agecoh_to_long(BCMuNSF, "Female", "Never"),
  mat_agecoh_to_long(BCMuFSF, "Female", "Former"),
  mat_agecoh_to_long(BCMuCSF, "Female", "Current"),
  mat_agecoh_to_long(BCMuNSM, "Male",   "Never"),
  mat_agecoh_to_long(BCMuFSM, "Male",   "Former"),
  mat_agecoh_to_long(BCMuCSM, "Male",   "Current")
)

# ---- Optional check: PROB -> RATE should reproduce the original rate (time = 1) ----
df.Mortality$Rate_check <- -log(1 - df.Mortality$Mortality_probability)  # / time, time=1
max_diff <- max(abs(df.Mortality$Mortality_rate - df.Mortality$Rate_check), na.rm = TRUE)
message("Max |Mortality_rate - (-log(1-Prob))| = ", signif(max_diff, 6))
df.Mortality$Rate_check <- NULL

# ---- Optional filter if you only want certain calendar years (uncomment as needed) ----
# df.Mortality <- dplyr::filter(df.Mortality, per >= 1994, per <= 2050)

# Write CSV
write.csv(
  df.Mortality,
  file = "Argentina_Mortality_Age_Per_Coh_Sex_SmokingStatus_Rate_Prob_ACM.csv",
  row.names = FALSE
)


