###Code to calculate mortality rates from Argentina using deaths and population files
### and then fill missing years and extrapolate to future years using age-period-cohort models
### V0 piloting. Rafael Meza, CISNET LWG / CEDES collaboration. Dec 2023 and may 2026

library(readr)
library(readxl)
library(tidyverse)
library(Epi)
library(dplyr)

rm(list = ls())

library(here)
here::i_am("Argentina_Mortality_Data_Projection_2026.R")

#####
#### Deaths file

# Read
def94a24 <- read_csv("../../DEIS/def94a24.csv")

# Eliminate those with no age
def94a24 = def94a24 %>% filter(!(uniedad%in%c(0,6,7,9)))

# replace those with age in hours, days or months as age 0 or those with more than 100 to 100
## still need to remove those with more than 12 months?
def94a24 = def94a24 %>% mutate(age=ifelse(uniedad%in%c(2,3,4,8),0,ifelse(uniedad==5,100,edad))) 
def94a24 <- def94a24[def94a24$ano != 1993, ] ##remove a few rows from 1993 (small number, there by mistake?)

## Separate Females and Males
def94a24F = def94a24[def94a24$sexo==2,]
def94a24M = def94a24[def94a24$sexo==1,]

### Count by year/age
deathsageyear = def94a24 %>% group_by(ano,age)%>%summarize(n=n())%>%ungroup()%>%complete(ano,age, fill = list(n = 0)) %>% filter(age>=0 & age<=99)
deathsageyearF = def94a24F %>% group_by(ano,age)%>%summarize(n=n())%>%ungroup()%>%complete(ano,age, fill = list(n = 0)) %>% filter(age>=0 & age<=99)
deathsageyearM = def94a24M %>% group_by(ano,age)%>%summarize(n=n())%>%ungroup()%>%complete(ano,age, fill = list(n = 0)) %>% filter(age>=0 & age<=99)

### Name as needed by Epi package (year, age, number of deaths)
colnames(deathsageyear)=c("P","A","D")
colnames(deathsageyearF)=c("P","A","D")
colnames(deathsageyearM)=c("P","A","D")

## delete death files from memmory
rm(def94a24,def94a24F,def94a24M)

#####
#### Pop  files 
#### Load Pop  files - UN projections
#2010-2100
#Females
pop1950_2100F=read.csv("../../../Argentina population/UN_Pop_1950-2100_BySingleAge_(2024_estimates)_Female.csv")
#Males
pop1950_2100M=read.csv("../../../Argentina population/UN_Pop_1950-2100_BySingleAge_(2024_estimates)_Male.csv")

# Keep only ages 0-99; delete first column (age)
pop1950_2100F = pop1950_2100F[1:100, 2:152]
pop1950_2100M = pop1950_2100M[1:100, 2:152]

# Remove the X before the year (in columns)
colnames(pop1950_2100F) <- sub("^X", "", colnames(pop1950_2100F))
colnames(pop1950_2100M) <- sub("^X", "", colnames(pop1950_2100M))

#####
#### Combine Deaths/Pop

deathsageyearF$Y=NA ## Y = population denominator
deathsageyearM$Y=NA

### Add population denominator Y to mortality datasets
### Mortality years: 1994-2024

for (yr in (1994:2024)) {
#for (yr in setdiff(1994:2024, 2021)) {  ##esto exluiría 2021 (COVID)
  # Find the population column corresponding to the current calendar year
  colyrF <- match(as.character(yr), colnames(pop1950_2100F))
  colyrM <- match(as.character(yr), colnames(pop1950_2100M))
  
  # Select rows in the mortality datasets for the current year and ages 0-99
  idxF <- deathsageyearF$P == yr & deathsageyearF$A >= 0 & deathsageyearF$A <= 99
  idxM <- deathsageyearM$P == yr & deathsageyearM$A >= 0 & deathsageyearM$A <= 99
  
  # Assign female population to Y
  deathsageyearF$Y[idxF] <- pop1950_2100F[
    deathsageyearF$A[idxF] + 1,
    colyrF
  ]
  
  # Assign male population to Y
  deathsageyearM$Y[idxM] <- pop1950_2100M[
    deathsageyearM$A[idxM] + 1,
    colyrM
  ]
}

deathsageyearF=na.omit(deathsageyearF)
deathsageyearM=na.omit(deathsageyearM)

deathsageyearF=deathsageyearF[,c(2,1,3,4)] ##reorder. now: A P D Y (age, year, age, pop)
deathsageyearM=deathsageyearM[,c(2,1,3,4)]


#####
#### FIT age-period-cohort models of mortality with Epi package
###WITHOUT SMOOTHING
#Females: 
ACPF=apc.fit(deathsageyearF,model='factor',dr.extr='Holford',parm='ACP',ref.c=1950)  # AC-P model - fit AC first - fix and fit F
APCF=apc.fit(deathsageyearF,model='factor',dr.extr='Holford',parm='APC',ref.p=2010)  # AC-P model - fit AC first - fix and fit F

#Males:
ACPM=apc.fit(deathsageyearM,model='factor',dr.extr='Holford',parm='ACP',ref.c=1950)  # AC-P model - fit AC first - fix and fit F
APCM=apc.fit(deathsageyearM,model='factor',dr.extr='Holford',parm='APC',ref.p=2010)  # AC-P model - fit AC first - fix and fit F

apc.plot(ACPF)
apc.lines(ACPM,col='blue')

par(mfrow=c(2,3))
plot(ACPF$Age[,1:2],main='Fem Age',log='y')
plot(ACPF$Per[,1:2],main='Fem Per')
plot(ACPF$Coh[,1:2],main='Fem Coh')
plot(ACPM$Age[,1:2],main='Male Age',log='y')
plot(ACPM$Per[,1:2],main='Male Per')
plot(ACPM$Coh[,1:2],main='Male Coh')

#####
### Project cohort effects until 2100. ### RM: updated May 8, 2024 
### Use last 20 years in y-log scale (linear)
### Project linearly, de-trending every 5 years (see Nordpred methodology)

x=2005:2024
yf=log(ACPF$Coh[(2005-1895+1):(2024-1895+1),2]) # log cohort effects
ym=log(ACPM$Coh[(2005-1895+1):(2024-1895+1),2])

CF=lm(yf~x)$coefficients ## fitted line constant and slope
CM=lm(ym~x)$coefficients

### Check fit - looks good
plot(x,yf)
lines(x,CF[1]+CF[2]*x)
plot(x,ym)
lines(x,CM[1]+CM[2]*x)

### Detrending consistent with Nordpred. See Chang J et al JGO 2018
detrend=(1-.08)^(1:76) # detrend 8% annually

for (i in 1:76){
  year=2024+i
  x=c(x,year)
  yf=c(yf,yf[19+i]+CF[2]*detrend[i]) # take previous value and increase with detrend slope
  ym=c(ym,ym[19+i]+CM[2]*detrend[i])
}

## check result in log and actual scale - looks good
par(mfrow=c(2,2))
plot(x,yf)
plot(x,exp(yf))
plot(x,ym)
plot(x,exp(ym))

### Cohort effects extrapolated -  right to 2100
BCEF=c(ACPF$Coh[,2],exp(yf[21:96]))
BCEM=c(ACPM$Coh[,2],exp(ym[21:96]))

# check - looks good
plot(1895:2100,BCEF)
lines(1895:2100,BCEM)

### Figure: projected cohort effects through 2100
log_BCEF <- log(BCEF)
log_BCEM <- log(BCEM)

years_cohort <- 1895:2100
keep_cohort <- years_cohort <= 2050

ylim_cohort <- range(
  c(log_BCEF[keep_cohort], log_BCEM[keep_cohort]),
  na.rm = TRUE
)

ylim_cohort <- ylim_cohort + c(-0.05, 0.05)

pdf("Projected_cohort_effects_until_2050.pdf")

plot(years_cohort[keep_cohort], log_BCEF[keep_cohort], type = "l",
     xlab = "Birth cohort",
     ylab = "Log cohort effect",
     ylim = ylim_cohort)

lines(years_cohort[keep_cohort], log_BCEM[keep_cohort], lty = 2)

legend("topright",
       legend = c("Females", "Males"),
       lty = c(1, 2),
       bty = "n")

dev.off()

PEF=cbind(rep(1,(2199-1895+1)))
PEF[2021-1895+1]=ACPF$Per[(2021-1994+1),2]
PEF[2022-1895+1]=ACPF$Per[(2022-1994+1),2]

PEM=cbind(rep(1,(2199-1895+1)))
PEM[2021-1895+1]=ACPM$Per[(2021-1994+1),2]
PEM[2022-1895+1]=ACPM$Per[(2022-1994+1),2]


#####
### Fill out mortality calendar year table from 1994-2050; period effects assumed to be 1
MortF=matrix(0,100,57)
MortM=matrix(0,100,57)

for (year in 1994:2050){
  for (age in 0:99){
    i=age+1  ## age index
    j=year-1993 ## year index
    k=year-age-1894 ## cohort index
    MortF[i,j]=ACPF$Age[i,2]*BCEF[k]*PEF[year-1895+1] ## apply ACP model Age*Cohort effects. Period effects are almost 1 (except PEF y PEM)
    MortM[i,j]=ACPM$Age[i,2]*BCEM[k]*PEM[year-1895+1]
  }
}

colnames(MortF)=1994:2050
colnames(MortM)=1994:2050


# Determine the maximum and minimum values of mortality rates across both MortF and MortM
min_mortality <- min(c(MortF, MortM), na.rm = TRUE)
max_mortality <- max(c(MortF, MortM), na.rm = TRUE)

### Export mortality file as csv
write.csv(MortF,"Fem_ArgentinaMortality_1994_2050.csv",row.names = FALSE)
write.csv(MortM,"Male_ArgentinaMortality_1994_2050.csv",row.names = FALSE)


#### Mortality by birth-cohort from cohorts born 1900 until 2100

BCMortF=matrix(0,100,201)
BCMortM=matrix(0,100,201)

for (byear in 1900:2100){
  k=byear-1894 ## cohort index
  j=byear-1899 ## column index for output matrix
  BCMortF[,j]=ACPF$Age[,2]*BCEF[k]*PEF[(byear-1895+1):((byear-1895+1+99))]
  BCMortM[,j]=ACPM$Age[,2]*BCEM[k]*PEM[(byear-1895+1):((byear-1895+1+99))]
}

colnames(BCMortF)=1900:2100
colnames(BCMortM)=1900:2100

### Export cohort mortality file as csv
write.csv(BCMortF,"Fem_ArgentinaCohortMortality_1900_2100.csv",row.names = FALSE)
write.csv(BCMortM,"Male_ArgentinaCohortMortality_1900_2100.csv",row.names = FALSE)


