# argentina-mortality-by-smoking-status
Argentina Mortality by Smoking Status

This repository contains R scripts used to estimate age-, sex-, and smoking status-specific mortality rates in Argentina and to generate life expectancy estimates by smoking status and smoking intensity.

## Associated manuscript

**Estimating mortality trends in Argentina by smoking status and smoking intensity**

M. Victoria Salgado, Jihyoun Jeon, Luis Zavala-Arciniega, Raul Mejia, Jamie Tam, Rafael Meza

American Journal of Preventive Medicine (AJPM)

## Repository contents

### 1. Argentina_Mortality_Data_Projection_2026.R

Estimation and projection of all-cause mortality rates in Argentina by age, sex, calendar year, and birth cohort.

### 2. Argentina_MortalityBy_SmokingStatus_2026.R

Estimation of mortality rates by smoking status (never, former, and current smokers).

### 3. Argentina_MortalityBy_SmokingStatus_6CPDcats_2026.R

Estimation of mortality rates among current smokers stratified by six categories of smoking intensity (cigarettes per day).

### 4. LE_age_40_by_BC_CPD.R

Generation of life expectancy at age 40 estimates by birth cohort and smoking intensity category.

### Supporting script

#### LE_Calc.R

Functions used for life table construction and life expectancy calculations. This script is sourced by other analytical workflows in the repository.

## Software

Analyses were performed in R.
