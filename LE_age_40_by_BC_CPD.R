#Life expectancy tables by birth cohort and smoking status using 6 CPD categories for Current Smokers
#In the CPD analysis, never- and former-smoker mortality rates are intermediate quantities used to obtain internally consistent CPD-specific current-smoker mortality rates. 
#Only current-smoker CPD categories are interpreted from this analysis.

rm(list = ls())

setwd("")

library(tidyverse)
library(RColorBrewer)

source("../../LE_Calc.R")  

read_mu <- function(path){
  x <- read.csv(path, header = TRUE, check.names = FALSE)
  # if R added "X" to numeric column names, remove it
  colnames(x) <- sub("^X", "", colnames(x))
  as.matrix(x)
}

# --- Read all 16 mortality tables (age 0–99 x cohorts 1900–2100) ---
MuNSF  <- read_mu("MuNSF_BirthYr_1900_2100.csv")
MuCS1F <- read_mu("MuCS1F_BirthYr_1900_2100.csv")
MuCS2F <- read_mu("MuCS2F_BirthYr_1900_2100.csv")
MuCS3F <- read_mu("MuCS3F_BirthYr_1900_2100.csv")
MuCS4F <- read_mu("MuCS4F_BirthYr_1900_2100.csv")
MuCS5F <- read_mu("MuCS5F_BirthYr_1900_2100.csv")
MuCS6F <- read_mu("MuCS6F_BirthYr_1900_2100.csv")
MuFSF  <- read_mu("MuFSF_BirthYr_1900_2100.csv")

MuNSM  <- read_mu("MuNSM_BirthYr_1900_2100.csv")
MuCS1M <- read_mu("MuCS1M_BirthYr_1900_2100.csv")
MuCS2M <- read_mu("MuCS2M_BirthYr_1900_2100.csv")
MuCS3M <- read_mu("MuCS3M_BirthYr_1900_2100.csv")
MuCS4M <- read_mu("MuCS4M_BirthYr_1900_2100.csv")
MuCS5M <- read_mu("MuCS5M_BirthYr_1900_2100.csv")
MuCS6M <- read_mu("MuCS6M_BirthYr_1900_2100.csv")
MuFSM  <- read_mu("MuFSM_BirthYr_1900_2100.csv")

# --- Compute LE tables (same dimensions as inputs) ---
LEN_SF  <- LEcalc(MuNSF)
LECS1_F <- LEcalc(MuCS1F)
LECS2_F <- LEcalc(MuCS2F)
LECS3_F <- LEcalc(MuCS3F)
LECS4_F <- LEcalc(MuCS4F)
LECS5_F <- LEcalc(MuCS5F)
LECS6_F <- LEcalc(MuCS6F)
LEFS_F  <- LEcalc(MuFSF)

LEN_SM  <- LEcalc(MuNSM)
LECS1_M <- LEcalc(MuCS1M)
LECS2_M <- LEcalc(MuCS2M)
LECS3_M <- LEcalc(MuCS3M)
LECS4_M <- LEcalc(MuCS4M)
LECS5_M <- LEcalc(MuCS5M)
LECS6_M <- LEcalc(MuCS6M)
LEFS_M  <- LEcalc(MuFSM)

# --- Extract LE at age 40: row 41 (because age 40 -> 40+1) ---
cohorts <- as.integer(colnames(MuNSF))

LE40_by_cohort <- rbind(
  `NS_F`  = LEN_SF[41,],
  `CS1_F` = LECS1_F[41,],
  `CS2_F` = LECS2_F[41,],
  `CS3_F` = LECS3_F[41,],
  `CS4_F` = LECS4_F[41,],
  `CS5_F` = LECS5_F[41,],
  `CS6_F` = LECS6_F[41,],
  `FS_F`  = LEFS_F[41,],
  `NS_M`  = LEN_SM[41,],
  `CS1_M` = LECS1_M[41,],
  `CS2_M` = LECS2_M[41,],
  `CS3_M` = LECS3_M[41,],
  `CS4_M` = LECS4_M[41,],
  `CS5_M` = LECS5_M[41,],
  `CS6_M` = LECS6_M[41,],
  `FS_M`  = LEFS_M[41,]
)

colnames(LE40_by_cohort) <- cohorts  # 1900:2100
write.csv(LE40_by_cohort, "Arg_LE40_ByBirthCohort_CPD_1900_2100.csv", row.names = TRUE)
