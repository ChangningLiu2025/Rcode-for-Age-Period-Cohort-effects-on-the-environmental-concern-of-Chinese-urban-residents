# Author：Changning Liu  2025
# e-mail:liuchn6@mail2.sysu.edu.cn 
# Age-Period-Cohort effects on the environmental concern of Chinese urban residents

library(cem)
library(dplyr)
library(lme4)
library(nlme)
library(lmerTest)
library(mice)
library(tidyr)
library(MuMIn)
library(ggplot2)
library(ggnewscale)
library(reshape2)
library(patchwork)
library(purrr)
library(ggeffects)



# 1. Data Processing ====
## 1.1 Data Imputation ====
replace_values_NEP <- function(data, col_name) {
  data[[col_name]] <- ifelse(data[[col_name]] == 1, 5,
                             ifelse(data[[col_name]] == 2, 4,
                                    ifelse(data[[col_name]] == 4, 2,
                                           ifelse(data[[col_name]] == 5, 1, data[[col_name]]))))
  return(data)
}

### 1.1.1 Data Processing 2003 ====
CGSS2003 <- read.csv('CGSS2003.csv')  # The original data downloaded from the CGSS website.
# Translate the data from 2003 including: ID, gender, year of birth, 15 questions 
# from the New Ecological Paradigm (NEP) Scale, subjective well-being, 
# years of education, local government environmental protection efforts, and annual household income.
# CGSS2003 consists entirely of urban samples.
CGSS2003_filter <- CGSS2003[,c('id','sex','birth','k1_1','k1_2','k1_3','k1_4',
                               'k1_5','k1_6','k1_7','k1_8','k1_9','k1_10','k1_11',
                               'k1_12', 'k1_13','k1_14','k1_15',
                               'lifefeel','educ','hhyrinc')]
# Rename column names
names(CGSS2003_filter)[names(CGSS2003_filter) == "lifefeel"] <- "happiness"
names(CGSS2003_filter)[names(CGSS2003_filter) == "educ"] <- "education"
names(CGSS2003_filter)[names(CGSS2003_filter) == "hhyrinc"] <- "household_income"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_1"] <- "NEP1"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_2"] <- "NEP2"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_3"] <- "NEP3"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_4"] <- "NEP4"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_5"] <- "NEP5"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_6"] <- "NEP6"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_7"] <- "NEP7"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_8"] <- "NEP8"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_9"] <- "NEP9"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_10"] <- "NEP10"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_11"] <- "NEP11"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_12"] <- "NEP12"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_13"] <- "NEP13"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_14"] <- "NEP14"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k1_15"] <- "NEP15"

# Remove samples under the age of 18 (To meet the requirements of the subsequent analysis, 
# samples under the age of 20 are directly excluded here).
CGSS2003_filter$birth <- as.numeric(CGSS2003_filter$birth)
CGSS2003_filter$age <- 2003 - CGSS2003_filter$birth
CGSS2003_filter <- CGSS2003_filter[CGSS2003_filter$age >= 20, ]

# Change samples with an NEP response of -2 (no answer) to NA.
CGSS2003_filter$NEP1[CGSS2003_filter$NEP1 %in% c(-2)] <- NA
CGSS2003_filter$NEP2[CGSS2003_filter$NEP2 %in% c(-2)] <- NA
CGSS2003_filter$NEP3[CGSS2003_filter$NEP3 %in% c(-2)] <- NA
CGSS2003_filter$NEP4[CGSS2003_filter$NEP4 %in% c(-2)] <- NA
CGSS2003_filter$NEP5[CGSS2003_filter$NEP5 %in% c(-2)] <- NA
CGSS2003_filter$NEP6[CGSS2003_filter$NEP6 %in% c(-2)] <- NA
CGSS2003_filter$NEP7[CGSS2003_filter$NEP7 %in% c(-2)] <- NA
CGSS2003_filter$NEP8[CGSS2003_filter$NEP8 %in% c(-2)] <- NA
CGSS2003_filter$NEP9[CGSS2003_filter$NEP9 %in% c(-2)] <- NA
CGSS2003_filter$NEP10[CGSS2003_filter$NEP10 %in% c(-2)] <- NA
CGSS2003_filter$NEP11[CGSS2003_filter$NEP11 %in% c(-2)] <- NA
CGSS2003_filter$NEP12[CGSS2003_filter$NEP12 %in% c(-2)] <- NA
CGSS2003_filter$NEP13[CGSS2003_filter$NEP13 %in% c(-2)] <- NA
CGSS2003_filter$NEP14[CGSS2003_filter$NEP14 %in% c(-2)] <- NA
CGSS2003_filter$NEP15[CGSS2003_filter$NEP15 %in% c(-2)] <- NA

# Change samples with missing data in birth year, gender, subjective well-being, 
# years of education, local government protection efforts, and household annual income to NA.
CGSS2003_filter$sex[CGSS2003_filter$sex %in% c(-3, -2, -1)] <- NA 
CGSS2003_filter$birth[CGSS2003_filter$birth %in% c(-3, -2, -1)] <- NA 
CGSS2003_filter$happiness[CGSS2003_filter$happiness %in% c(-3, -2, -1)] <- NA 
CGSS2003_filter$education[CGSS2003_filter$education %in% c(-3, -2, -1)] <- NA 
CGSS2003_filter$household_income[CGSS2003_filter$household_income %in% c(-3, -2, -1)] <- NA 

# Due to the settings of the 2003 data, adjust the NEP values to ensure 
# consistency with the 2010 and 2021 data format, where:  
# 1 = Strongly disagree  
# 2 = Somewhat disagree  
# 3 = Neutral  
# 4 = Somewhat agree  
# 5 = Strongly agree
col_names_list_2003 <- c("NEP1", "NEP2", "NEP3",'NEP4','NEP5','NEP6','NEP7','NEP8','NEP9','NEP10','NEP11','NEP12','NEP13','NEP14',
                         'NEP15')
for (col_name in col_names_list_2003) {
  CGSS2003_filter <- replace_values_NEP(CGSS2003_filter, col_name)
}

# Define years of education based on the data as follows:  
#   - No formal education = 0  
# - Primary school and private school = 6  
# - Junior high school = 9  
# - Secondary vocational school, vocational high school, technical school, general high school = 12  
# - Associate degree (regular and adult higher education) = 15  
# - Bachelor's degree (regular and adult higher education) = 16  
# - Graduate degree and above = 19          
CGSS2003_filter$education <- ifelse(CGSS2003_filter$education %in% c(0), 0, 
                                    ifelse(CGSS2003_filter$education %in% c(1,11), 6, 
                                           ifelse(CGSS2003_filter$education %in% c(2), 9,
                                                  ifelse(CGSS2003_filter$education %in% c(3,4,5), 12,
                                                         ifelse(CGSS2003_filter$education %in% c(6,7), 15,
                                                                ifelse(CGSS2003_filter$education %in% c(8,9), 16,
                                                                       ifelse(CGSS2003_filter$education %in% c(10), 19,CGSS2003_filter$education)))))))


# Removed samples that did not answer any Chinese NEP questions 
selected_cols2003 <- c("NEP1", "NEP3",'NEP5','NEP7',
                       'NEP8','NEP9','NEP10','NEP11','NEP13',
                       'NEP15')
missing_values2003 <- rowSums(is.na(CGSS2003_filter[selected_cols2003]))
total_selected_cols2003 <- length(selected_cols2003)
missing_values_ratio2003 <- missing_values2003 / total_selected_cols2003
CGSS2003_filter$missing_values_ratio2003 <- missing_values_ratio2003
CGSS2003_filter <- CGSS2003_filter[CGSS2003_filter$missing_values_ratio2003 < 1, ]

# Remove ‘id’ and ‘missing_values_ratio2003’ to avoid interference during multiple imputation.
CGSS2003_filter <- subset(CGSS2003_filter, select = -id)
CGSS2003_filter <- subset(CGSS2003_filter, select = -missing_values_ratio2003)

# # Define the variables as follows:  
# - **Continuous variables**: Age, household annual income, years of education  
# - **Unordered categorical variables**: Gender, housing type, marital status  
# - **Ordered categorical variables**: NEP, subjective well-being
CGSS2003_filter$age <- as.numeric(CGSS2003_filter$age)
CGSS2003_filter$household_income <- as.numeric(CGSS2003_filter$household_income)
CGSS2003_filter$household_income <- scale(CGSS2003_filter$household_income) # Standardize household income to ensure consistency in variable scales.
CGSS2003_filter$household_income <- as.numeric(CGSS2003_filter$household_income)
CGSS2003_filter$education <- as.numeric(CGSS2003_filter$education)
CGSS2003_filter <- CGSS2003_filter %>% mutate_at(vars(sex), as.factor)
NEP_order <- c("1", "2", "3",'4','5')
CGSS2003_filter <- CGSS2003_filter %>% mutate_at(vars(NEP1,NEP2,NEP3,NEP4,NEP5,NEP6,NEP7,NEP8,NEP9,NEP10,NEP11,NEP12,NEP13,NEP14,NEP15), 
                                                 list(~ factor(., levels = NEP_order, ordered = TRUE)))
happy_order <- c("1", "2", "3",'4','5')
CGSS2003_filter <- CGSS2003_filter %>% mutate(happiness = factor(happiness, levels = happy_order, ordered = TRUE))

# Sort the data according to the specified order.
custom_order <- c("NEP1", "NEP2", "NEP3",'NEP4','NEP5','NEP6','NEP7','NEP8','NEP9','NEP10','NEP11','NEP12','NEP13','NEP14',
                  'NEP15','sex','age','happiness','education','household_income','birth')# ,rental,marriage
CGSS2003_filter <- CGSS2003_filter %>% select(all_of(custom_order))
CGSS2003_filter$Year <- 2003

# Write the data before imputation (i.e., CGSS2003 in Supplementary Table 7).
# write.csv(CGSS2003_filter, file = 'CGSS2003(插补前).csv', row.names = FALSE)

# 排除出生年份和时间（因为其与age完全共线性）
CGSS2003_filter <- subset(CGSS2003_filter, select = -c(Year, birth))

# Data imputation
impmice=mice(CGSS2003_filter,m=5,seed=2003)

# Write the 5 five sets of imputed datasets 2003.
n_iterations <- 5
for (i in 1:n_iterations) {
  assign(paste0("imp", i), complete(impmice, i))
  write.csv(get(paste0("imp", i)), file = paste0("2003_imp", i, ".csv"), row.names = FALSE)
}


### 1.1.2 Data Processing 2010 ====
CGSS2010 <- read.csv('CGSS2010.csv') # The original data downloaded from the CGSS website.
# Translate the data from 2003 including: ID, gender, year of birth, 15 questions 
# from the New Ecological Paradigm (NEP) Scale, subjective well-being, sample type
# years of education, local government environmental protection efforts, and annual household income.
CGSS2010_filter <- CGSS2010[,c('id','a2','a3a','l2501','l2502','l2503','l2504','l2505',
                               'l2506','l2507','l2508','l2509','l2510','l2511','l2512',
                               'l2513','l2514','l2515','a36','a7a','a62','s5')]

# Rename column names
names(CGSS2010_filter)[names(CGSS2010_filter) == "a2"] <- "sex"
names(CGSS2010_filter)[names(CGSS2010_filter) == "a3a"] <- "birth"
names(CGSS2010_filter)[names(CGSS2010_filter) == "a36"] <- "happiness"
names(CGSS2010_filter)[names(CGSS2010_filter) == "a7a"] <- "education"
names(CGSS2010_filter)[names(CGSS2010_filter) == "a62"] <- "household_income"
names(CGSS2010_filter)[names(CGSS2010_filter) == "s5"] <- "type"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2501"] <- "NEP1"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2502"] <- "NEP2"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2503"] <- "NEP3"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2504"] <- "NEP4"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2505"] <- "NEP5"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2506"] <- "NEP6"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2507"] <- "NEP7"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2508"] <- "NEP8"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2509"] <- "NEP9"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2510"] <- "NEP10"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2511"] <- "NEP11"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2512"] <- "NEP12"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2513"] <- "NEP13"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2514"] <- "NEP14"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l2515"] <- "NEP15"

# Remove samples under the age of 18 (To meet the requirements of the subsequent analysis, 
# samples under the age of 20 are directly excluded here).
CGSS2010_filter$birth <- as.numeric(CGSS2010_filter$birth)
CGSS2010_filter$age <- 2010 - CGSS2010_filter$birth
CGSS2010_filter <- CGSS2010_filter[CGSS2010_filter$age >= 20, ]
# Select the samples with the type "urban."
CGSS2010_filter <- CGSS2010_filter[CGSS2010_filter$type == 1, ]

# Change samples with an NEP response of -3（refuse to answer）,
# -2（don't know）,8 (unable to choose) to NA.
CGSS2010_filter$NEP1[CGSS2010_filter$NEP1 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP2[CGSS2010_filter$NEP2 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP3[CGSS2010_filter$NEP3 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP4[CGSS2010_filter$NEP4 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP5[CGSS2010_filter$NEP5 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP6[CGSS2010_filter$NEP6 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP7[CGSS2010_filter$NEP7 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP8[CGSS2010_filter$NEP8 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP9[CGSS2010_filter$NEP9 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP10[CGSS2010_filter$NEP10 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP11[CGSS2010_filter$NEP11 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP12[CGSS2010_filter$NEP12 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP13[CGSS2010_filter$NEP13 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP14[CGSS2010_filter$NEP14 %in% c(-3, -2, 8)] <- NA
CGSS2010_filter$NEP15[CGSS2010_filter$NEP15 %in% c(-3, -2, 8)] <- NA

# Change samples with missing data in birth year, gender, subjective well-being, 
# years of education, local government protection efforts, and household annual income to NA.
CGSS2010_filter$sex[CGSS2010_filter$sex %in% c(-3, -2, -1)] <- NA 
CGSS2010_filter$birth[CGSS2010_filter$birth %in% c(-3)] <- NA 
CGSS2010_filter$happiness[CGSS2010_filter$happiness %in% c(-3, -2)] <- NA 
CGSS2010_filter$education[CGSS2010_filter$education %in% c(-3, 14)] <- NA 
CGSS2010_filter$household_income[CGSS2010_filter$household_income %in% c(9999997, 9999998, 9999999)] <- NA 

# Define years of education based on the data as follows:  
#   - No formal education = 0  
# - Primary school and private school = 6  
# - Junior high school = 9  
# - Secondary vocational school, vocational high school, technical school, general high school = 12  
# - Associate degree (regular and adult higher education) = 15  
# - Bachelor's degree (regular and adult higher education) = 16  
# - Graduate degree and above = 19                   
CGSS2010_filter$education <- ifelse(CGSS2010_filter$education %in% c(1), 0, 
                                    ifelse(CGSS2010_filter$education %in% c(2,3), 6, 
                                           ifelse(CGSS2010_filter$education %in% c(4), 9,
                                                  ifelse(CGSS2010_filter$education %in% c(5,6,7,8), 12,
                                                         ifelse(CGSS2010_filter$education %in% c(9,10), 15,
                                                                ifelse(CGSS2010_filter$education %in% c(11,12), 16,
                                                                       ifelse(CGSS2010_filter$education %in% c(13), 19,CGSS2010_filter$education)))))))



# Removed samples that did not answer any Chinese NEP questions 
selected_cols2010 <- c("NEP1", "NEP3",'NEP5','NEP7',
                       'NEP8','NEP9','NEP10','NEP11','NEP13',
                       'NEP15')
missing_values2010 <- rowSums(is.na(CGSS2010_filter[selected_cols2010]))
total_selected_cols2010 <- length(selected_cols2010)
missing_values_ratio2010 <- missing_values2010 / total_selected_cols2010
CGSS2010_filter$missing_values_ratio2010 <- missing_values_ratio2010
CGSS2010_filter <- CGSS2010_filter[CGSS2010_filter$missing_values_ratio2010 < 1, ]

# Remove ‘id’, 'type' and ‘missing_values_ratio2003’ to avoid interference during multiple imputation.
CGSS2010_filter <- subset(CGSS2010_filter, select = -id)
CGSS2010_filter <- subset(CGSS2010_filter, select = -missing_values_ratio2010)
CGSS2010_filter <- subset(CGSS2010_filter, select = -type)

# # Define the variables as follows:  
# - **Continuous variables**: Age, household annual income, years of education  
# - **Unordered categorical variables**: Gender, housing type, marital status  
# - **Ordered categorical variables**: NEP, subjective well-being
CGSS2010_filter$age <- as.numeric(CGSS2010_filter$age)
CGSS2010_filter$household_income <- as.numeric(CGSS2010_filter$household_income)
CGSS2010_filter$household_income <- scale(CGSS2010_filter$household_income) # Standardize household income to ensure consistency in variable scales.
CGSS2010_filter$household_income <- as.numeric(CGSS2010_filter$household_income)
CGSS2010_filter$education <- as.numeric(CGSS2010_filter$education)
CGSS2010_filter <- CGSS2010_filter %>% mutate_at(vars(sex), as.factor)
NEP_order <- c("1", "2", "3",'4','5')
CGSS2010_filter <- CGSS2010_filter %>% mutate_at(vars(NEP1,NEP2,NEP3,NEP4,NEP5,NEP6,NEP7,NEP8,NEP9,NEP10,NEP11,NEP12,NEP13,NEP14,NEP15), 
                                                 list(~ factor(., levels = NEP_order, ordered = TRUE)))
happy_order <- c("1", "2", "3",'4','5')
CGSS2010_filter <- CGSS2010_filter %>% mutate(happiness = factor(happiness, levels = happy_order, ordered = TRUE))

# Sort the data according to the specified order.
custom_order <- c("NEP1", "NEP2", "NEP3",'NEP4','NEP5','NEP6','NEP7','NEP8','NEP9','NEP10','NEP11','NEP12','NEP13','NEP14',
                  'NEP15','sex','age','happiness','education','household_income','birth') #'rental','marriage',
CGSS2010_filter <- CGSS2010_filter %>% select(all_of(custom_order))
CGSS2010_filter$Year <- 2010

# Write the data before imputation (i.e., CGSS2010 in Supplementary Table 7).
# write.csv(CGSS2010_filter, file = 'CGSS2010(插补前).csv', row.names = FALSE)

# 排除出生年份和时间（因为其与age完全共线性）
CGSS2010_filter <- subset(CGSS2010_filter, select = -c(Year, birth))

# Data imputation
impmice=mice(CGSS2010_filter,m=5,seed=2003)

# Write the 5 five sets of imputed datasets 2010.
n_iterations <- 5
for (i in 1:n_iterations) {
  assign(paste0("imp", i), complete(impmice, i))
  write.csv(get(paste0("imp", i)), file = paste0("2010_imp", i, ".csv"), row.names = FALSE)
}


### 1.1.3 Data Processing 2021 ====
CGSS2021 <- read.csv('CGSS2021.csv') # The original data downloaded from the CGSS website.
# Translate the data from 2003 including: ID, gender, year of birth, 15 questions 
# from the New Ecological Paradigm (NEP) Scale, subjective well-being, sample type
# years of education, local government environmental protection efforts, and annual household income.
CGSS2021_filter <- CGSS2021[,c('id','A2','A3_1','H12_1','H12_2','H12_3','H12_4','H12_5','H12_6',
                               'H12_7','H12_8','H12_9','H12_10','H12_11','H12_12','H12_13','H12_14',
                               'H12_15','A36','A7a','A62','type')]

# Rename column names
names(CGSS2021_filter)[names(CGSS2021_filter) == "A2"] <- "sex"
names(CGSS2021_filter)[names(CGSS2021_filter) == "A3_1"] <- "birth"
names(CGSS2021_filter)[names(CGSS2021_filter) == "A36"] <- "happiness"
names(CGSS2021_filter)[names(CGSS2021_filter) == "A7a"] <- "education"
names(CGSS2021_filter)[names(CGSS2021_filter) == "A62"] <- "household_income"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_1"] <- "NEP1"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_2"] <- "NEP2"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_3"] <- "NEP3"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_4"] <- "NEP4"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_5"] <- "NEP5"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_6"] <- "NEP6"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_7"] <- "NEP7"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_8"] <- "NEP8"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_9"] <- "NEP9"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_10"] <- "NEP10"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_11"] <- "NEP11"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_12"] <- "NEP12"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_13"] <- "NEP13"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_14"] <- "NEP14"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H12_15"] <- "NEP15"

# Remove samples under the age of 18 (To meet the requirements of the subsequent analysis, 
# samples under the age of 20 are directly excluded here).
CGSS2021_filter$birth <- as.numeric(CGSS2021_filter$birth)
CGSS2021_filter$age <- 2021 - CGSS2021_filter$birth
CGSS2021_filter <- CGSS2021_filter[CGSS2021_filter$age >= 20, ]
# Select the samples with the type "urban."
CGSS2021_filter <- CGSS2021_filter[CGSS2021_filter$type == 1, ]

# Change samples with an NEP response of 98(unable to choose) to NA.
CGSS2021_filter$NEP1[CGSS2021_filter$NEP1 %in% c(98)] <- NA
CGSS2021_filter$NEP2[CGSS2021_filter$NEP2 %in% c(98)] <- NA
CGSS2021_filter$NEP3[CGSS2021_filter$NEP3 %in% c(98)] <- NA
CGSS2021_filter$NEP4[CGSS2021_filter$NEP4 %in% c(98)] <- NA
CGSS2021_filter$NEP5[CGSS2021_filter$NEP5 %in% c(98)] <- NA
CGSS2021_filter$NEP6[CGSS2021_filter$NEP6 %in% c(98)] <- NA
CGSS2021_filter$NEP7[CGSS2021_filter$NEP7 %in% c(98)] <- NA
CGSS2021_filter$NEP8[CGSS2021_filter$NEP8 %in% c(98)] <- NA
CGSS2021_filter$NEP9[CGSS2021_filter$NEP9 %in% c(98)] <- NA
CGSS2021_filter$NEP10[CGSS2021_filter$NEP10 %in% c(98)] <- NA
CGSS2021_filter$NEP11[CGSS2021_filter$NEP11 %in% c(98)] <- NA
CGSS2021_filter$NEP12[CGSS2021_filter$NEP12 %in% c(98)] <- NA
CGSS2021_filter$NEP13[CGSS2021_filter$NEP13 %in% c(98)] <- NA
CGSS2021_filter$NEP14[CGSS2021_filter$NEP14 %in% c(98)] <- NA
CGSS2021_filter$NEP15[CGSS2021_filter$NEP15 %in% c(98)] <- NA

# Change samples with missing data in birth year, gender, subjective well-being, 
# years of education, local government protection efforts, and household annual income to NA.
CGSS2021_filter$happiness[CGSS2021_filter$happiness %in% c(98,99)] <- NA
CGSS2021_filter$education[CGSS2021_filter$education %in% c(14)] <- NA
CGSS2021_filter$household_income[CGSS2021_filter$household_income %in% c(9999997,9999998, 9999999)] <- NA
# The household annual income of 9999996 is higher than one million, representing the top 1% and is the maximum value. 
# Assign it a value of 1,000,000, and after standardization, it will still be the maximum value.
CGSS2021_filter$household_income[CGSS2021_filter$household_income == 9999996] <- 1000000

# Define years of education based on the data as follows:  
#   - No formal education = 0  
# - Primary school and private school = 6  
# - Junior high school = 9  
# - Secondary vocational school, vocational high school, technical school, general high school = 12  
# - Associate degree (regular and adult higher education) = 15  
# - Bachelor's degree (regular and adult higher education) = 16  
# - Graduate degree and above = 19           
CGSS2021_filter$education <- ifelse(CGSS2021_filter$education %in% c(1), 0, 
                                    ifelse(CGSS2021_filter$education %in% c(2,3), 6, 
                                           ifelse(CGSS2021_filter$education %in% c(4), 9,
                                                  ifelse(CGSS2021_filter$education %in% c(5,6,7,8), 12,
                                                         ifelse(CGSS2021_filter$education %in% c(9,10), 15,
                                                                ifelse(CGSS2021_filter$education %in% c(11,12), 16,
                                                                       ifelse(CGSS2021_filter$education %in% c(13), 19,CGSS2021_filter$education)))))))


# Removed samples that did not answer any Chinese NEP questions 
selected_cols2021 <- c("NEP1", "NEP3",'NEP5','NEP7',
                       'NEP8','NEP9','NEP10','NEP11','NEP13',
                       'NEP15')
missing_values2021 <- rowSums(is.na(CGSS2021_filter[selected_cols2021]))
total_selected_cols2021 <- length(selected_cols2021)
missing_values_ratio2021 <- missing_values2021 / total_selected_cols2021
CGSS2021_filter$missing_values_ratio2021 <- missing_values_ratio2021
CGSS2021_filter <- CGSS2021_filter[CGSS2021_filter$missing_values_ratio2021 < 1, ]

# Remove ‘id’,'type' and ‘missing_values_ratio2003’ to avoid interference during multiple imputation.
CGSS2021_filter <- subset(CGSS2021_filter, select = -id)
CGSS2021_filter <- subset(CGSS2021_filter, select = -missing_values_ratio2021)
CGSS2021_filter <- subset(CGSS2021_filter, select = -type)

# # Define the variables as follows:  
# - **Continuous variables**: Age, household annual income, years of education  
# - **Unordered categorical variables**: Gender, housing type, marital status  
# - **Ordered categorical variables**: NEP, subjective well-being
CGSS2021_filter$age <- as.numeric(CGSS2021_filter$age)
CGSS2021_filter$household_income <- as.numeric(CGSS2021_filter$household_income)
CGSS2021_filter$household_income <- scale(CGSS2021_filter$household_income) # Standardize household income to ensure consistency in variable scales.
CGSS2021_filter$household_income <- as.numeric(CGSS2021_filter$household_income)
CGSS2021_filter$education <- as.numeric(CGSS2021_filter$education)
CGSS2021_filter <- CGSS2021_filter %>% mutate_at(vars(sex), as.factor)
NEP_order <- c("1", "2", "3",'4','5')
CGSS2021_filter <- CGSS2021_filter %>% mutate_at(vars(NEP1,NEP2,NEP3,NEP4,NEP5,NEP6,NEP7,NEP8,NEP9,NEP10,NEP11,NEP12,NEP13,NEP14,NEP15), 
                                                 list(~ factor(., levels = NEP_order, ordered = TRUE)))
happy_order <- c("1", "2", "3",'4','5')
CGSS2021_filter <- CGSS2021_filter %>% mutate(happiness = factor(happiness, levels = happy_order, ordered = TRUE))

# Sort the data according to the specified order.
custom_order <- c("NEP1", "NEP2", "NEP3",'NEP4','NEP5','NEP6','NEP7','NEP8','NEP9','NEP10','NEP11','NEP12','NEP13','NEP14',
                  'NEP15','sex','age','happiness','education','household_income','birth')#'rental','marriage',
CGSS2021_filter <- CGSS2021_filter %>% select(all_of(custom_order))
CGSS2021_filter$Year <- 2021

# Write the data before imputation (i.e., CGSS2021 in Supplementary Table 7).
# write.csv(CGSS2021_filter, file = 'CGSS2021(插补前).csv', row.names = FALSE)

# 排除出生年份和时间（因为其与age完全共线性）
CGSS2021_filter <- subset(CGSS2021_filter, select = -c(Year, birth))

# Data imputation
impmice=mice(CGSS2021_filter,m=5,seed=2003)

# Write the 5 five sets of imputed datasets 2021.
n_iterations <- 5
for (i in 1:n_iterations) {
  assign(paste0("imp", i), complete(impmice, i))
  write.csv(get(paste0("imp", i)), file = paste0("2021_imp", i, ".csv"), row.names = FALSE)
}


### 1.1.4 Figure 1 ====
# setwd('C:\\') Your own path
CGSS2003 <- read.csv('CGSS2003(插补前).csv') # The data removed before imputation in section 1.2.1
CGSS2010 <- read.csv('CGSS2010(插补前).csv') # The data removed before imputation in section 1.2.2
CGSS2021 <- read.csv('CGSS2021(插补前).csv') # The data removed before imputation in section 1.2.3
CGSS <- rbind(CGSS2003,CGSS2010,CGSS2021) 
CGSS <- CGSS %>% mutate(cohort_cat = case_when(
  birth >= 1920 & birth <= 1929 ~ 1,
  birth >= 1930 & birth <= 1939 ~ 2,
  birth >= 1940 & birth <= 1949 ~ 3,
  birth >= 1950 & birth <= 1959 ~ 4,
  birth >= 1960 & birth <= 1969 ~ 5,
  birth >= 1970 & birth <= 1979 ~ 6,
  birth >= 1980 & birth <= 1989 ~ 7,
  birth >= 1990 & birth <= 1999 ~ 8,
  birth >= 2000 & birth <= 2009 ~ 9,
  TRUE ~ NA_integer_  
))
CGSS <- CGSS[!is.na(CGSS$cohort_cat), ] 
CGSS <- CGSS %>%
  mutate(cohort_cat = recode(cohort_cat,
                             `1` = "1920-1929",
                             `2` = "1930-1939",
                             `3` = "1940-1949",
                             `4` = "1950-1959",
                             `5` = "1960-1969",
                             `6` = "1970-1979",
                             `7` = "1980-1989",
                             `8` = "1990-1999",
                             `9` = "2000-2009"))
CGSS$Year <- as.factor(CGSS$Year)
CGSS$cohort_cat <- as.factor(CGSS$cohort_cat)
col_names_list_reset <- c( "NEP2", 'NEP4','NEP6','NEP8','NEP10','NEP12','NEP14')
replace_values_NEP <- function(data, col_name) {
  data[[col_name]] <- ifelse(data[[col_name]] == 1, 5,
                             ifelse(data[[col_name]] == 2, 4,
                                    ifelse(data[[col_name]] == 4, 2,
                                           ifelse(data[[col_name]] == 5, 1, data[[col_name]]))))
  return(data)
}
for (col_name in col_names_list_reset) {
  CGSS <- replace_values_NEP(CGSS, col_name)
}
CGSS <- CGSS %>% mutate(CNEP=NEP1+NEP3+NEP5+NEP7+NEP8+NEP9+NEP10+NEP11+NEP13+NEP15)

#### 1.1.4.1 Age-CNEP Period Plots ====
point_colors <- c("#66c2a5", "#fc8d62", "#8da0cb")  # 温和的Set2颜色
line_colors <- c("#4cd964", "#ff9500", "#8a6de9")   # 更亮的版本

p1 <- ggplot(CGSS, aes(x = age, y = CNEP)) +
  geom_point(aes(color = Year), 
             size = 3.5, shape = 16, alpha = 0.1, 
             show.legend = FALSE) +  
  scale_color_manual(values = setNames(point_colors, levels(CGSS$Year))) +
  new_scale_color() +
  geom_smooth(aes(color = Year, group = Year), 
              method = "loess", se = FALSE, size = 1.2, n = 10000) +
  scale_color_manual(
    name = "Year",  
    values = setNames(line_colors, levels(CGSS$Year))
  ) +
  labs(x = "Age", y = "CNEP") +
  theme_minimal()

ggsave("Age_CNEP_Period_Plots.png", 
       plot = p1,
       width = 6,
       height = 4,
       dpi = 1200,
       units = "in")    

#### 1.1.4.1 Age-CNEP Cohort Plots ====
point_colors_set2 <- c(
  "#66c2a5",  # 1. 浅绿色
  "#fc8d62",  # 2. 浅橙色
  "#8da0cb",  # 3. 浅蓝紫色
  "#e78ac3",  # 4. 浅粉色
  "#a6d854",  # 5. 浅黄绿色
  "#ffd92f",  # 6. 浅黄色
  "#e5c494",  # 7. 浅米色
  "#b3b3b3",  # 8. 浅灰色
  "#80b1d3"   # 9. 补充：浅蓝色（来自Set3）
)
line_colors_brighter <- c(
  "#2ecc71",  # 1. 亮绿色 (更亮)
  "#ff7b24",  # 2. 亮橙色 (更亮)
  "#5e9cf0",  # 3. 亮蓝紫色 (更亮)
  "#ff6b9d",  # 4. 亮粉色 (更亮)
  "#8cc63f",  # 5. 亮黄绿色 (更亮)
  "#ffde17",  # 6. 亮黄色 (更亮)
  "#f7b977",  # 7. 亮米色 (更亮)
  "#808080",  # 8. 中灰色 (稍微深一点，因为浅灰不好调亮)
  "#45b7d1"   # 9. 亮蓝色 (更亮)
)

p2 <- ggplot(CGSS, aes(x = age, y = CNEP)) +
  geom_point(aes(color = cohort_cat), 
             size = 3.5, shape = 16, alpha = 0.1, 
             show.legend = FALSE) +  
  scale_color_manual(values = setNames(point_colors_set2, levels(CGSS$cohort_cat))) +
  new_scale_color() +
  geom_smooth(aes(color = cohort_cat, group = cohort_cat), 
              method = "loess", se = FALSE, size = 1.2, n = 10000) +
  scale_color_manual(
    name = "Cohort",  
    values = setNames(line_colors_brighter, levels(CGSS$cohort_cat))
  ) +
  labs(x = "Age", y = "CNEP") +
  theme_minimal()

ggsave("Age_CNEP_Cohort_Plots.png", 
       plot = p2,
       width = 7,
       height = 4,
       dpi = 1200,
       units = "in") 



#### 1.1.4.3 Figure 1 ====
library(patchwork)
Figure1 <- (p1 | p2) 
Figure1 <- Figure1 + plot_annotation(
  tag_levels = 'a') 
ggsave("Figure1.png", 
       plot = Figure1,
       width = 10,
       height = 4,
       dpi = 1200,
       units = "in") 


## 1.2 Combined the Imputed Datasets ====
# files_2003 <- paste0("2003_imp", 1:5, ".csv")
# files_2010 <- paste0("2010_imp", 1:5, ".csv")
# files_2021 <- paste0("2021_imp", 1:5, ".csv")

# combined the imputed datasets to create 125 complete datasets
for (i in 1:5) {
  for (j in 1:5) {
    for (k in 1:5) {
      # Construct the file path.
      file_2003 <- files_2003[i]
      file_2010 <- files_2010[j]
      file_2021 <- files_2021[k]
      
      # Read the data.(Your own path)
      df_2003 <- read.csv(paste0("D:\\博士阶段\\各类材料\\论文\\2026-APC\\数据\\插补后数据\\", file_2003))
      df_2010 <- read.csv(paste0("D:\\博士阶段\\各类材料\\论文\\2026-APC\\数据\\插补后数据\\", file_2010))
      df_2021 <- read.csv(paste0("D:\\博士阶段\\各类材料\\论文\\2026-APC\\数据\\插补后数据\\", file_2021))
      
      # Combined the data.
      combined_data <- rbind(df_2003, df_2010, df_2021)
      
      # Construct the output file path.
      output_file <- paste0("D:\\博士阶段\\各类材料\\论文\\2026-APC\\数据\\插补后数据\\", (i - 1) * 25 + (j - 1) * 5 + k, ".csv")
      
      # Write to a new file.
      write.csv(combined_data, file = output_file, row.names = FALSE)
    }
  }
}



# 2. Main Statistical Analysis ====
## 2.1 Main Analysis ====
# The data analysis section of the main text.
# 125 complete datasets

### 2.1.1 Read the Data ====
# setwd('D:\\') Your own path
# read all 125 complete datasets
file_list <- list.files(pattern = "^.*\\.csv$")
data_list <- lapply(file_list, read.csv) 

# Reverse the scores for reverse-coded questions.
col_names_list_reset <- c( "NEP2", 'NEP4','NEP6','NEP8','NEP10','NEP12','NEP14')
replace_values_NEP <- function(data, col_name) {
  data[[col_name]] <- ifelse(data[[col_name]] == 1, 5,
                             ifelse(data[[col_name]] == 2, 4,
                                    ifelse(data[[col_name]] == 4, 2,
                                           ifelse(data[[col_name]] == 5, 1, data[[col_name]]))))
  return(data)
}

data_list_modified <- lapply(data_list, function(data) {
  for (col_name in col_names_list_reset) {
    data <- replace_values_NEP(data, col_name)
  }
  return(data)
})

# Calculate the age used in the analysis (AgeR:divide by ten and subtract 2), as well as its square, 
# classify the generations(cohort) and calculate CNEP.
data_list_modified <- lapply(data_list_modified, function(data) {
  data <- data %>%
    mutate(CNEP = NEP1 + NEP3 + NEP5 + NEP7 + NEP8 + NEP9 + NEP10 + NEP11 + NEP13 + NEP15)
  data$birth <- data$Year - data$age
  data$AgeR = data$age/10 - 2
  data$AgeR2 = data$AgeR^2
  data <- data %>% mutate(cohort_cat = case_when(
    birth >= 1920 & birth <= 1929 ~ 1,
    birth >= 1930 & birth <= 1939 ~ 2,
    birth >= 1940 & birth <= 1949 ~ 3,
    birth >= 1950 & birth <= 1959 ~ 4,
    birth >= 1960 & birth <= 1969 ~ 5,
    birth >= 1970 & birth <= 1979 ~ 6,
    birth >= 1980 & birth <= 1989 ~ 7,
    birth >= 1990 & birth <= 1999 ~ 8,
    birth >= 2000 & birth <= 2009 ~ 9,
    TRUE ~ NA_integer_  # Define the samples that do not belong to these generations as NA.
  ))
  data <- data[!is.na(data$cohort_cat), ] 
  data$cohort_cat <- as.factor(data$cohort_cat)
  data <- data %>% mutate_at(vars(sex,Year,happiness), as.factor)
  data$birth_cent <- data$birth - mean(data$birth)
  return(data)
})

### 2.1.2 Perform Coarsened Exact Matching ====
matstrict <- cem("Year", datalist=data_list_modified, drop=c('NEP1','NEP2','NEP3','NEP4','NEP5','NEP6','NEP7','NEP8','NEP9',
                                                             'NEP10','NEP11','NEP12','NEP13','NEP14','NEP15',
                                                             'age','birth','Year','CNEP',
                                                             'AgeR','AgeR2','cohort_cat','birth_cent'))
matstrict
matched_data_list <- list()
for (i in seq_along(data_list)) {
  match_col <- paste0("match", i)
  matched_data_list[[i]] <- data_list_modified[[i]][matstrict[[match_col]]$matched, ]
}

### 2.1.3 Construct the HAPC model for each matched dataset ====
result_CNEP <- data.frame()
for (i in 1:length(matched_data_list)) {
  data <- matched_data_list[[i]]
  result_tem <- data.frame() 
  
  fullmodel <- lmer(CNEP ~ AgeR + AgeR2 + birth_cent + sex + happiness + education + household_income
                    + (1 | Year) + (1 | cohort_cat), data = data)
  summary_fullmodel <- summary(fullmodel)
  
  # Model information
  result_tem[1,'Model'] <- paste0("data", i)
  result_tem[1,'samplesize'] <- nrow(data)
  
  # variance_cohort & variance_period & variance_Residual
  random_effects_var <- as.data.frame(VarCorr(fullmodel))
  result_tem[1,'variance_cohort'] <- random_effects_var[1,4]
  result_tem[1,'variance_period'] <- random_effects_var[2,4]
  result_tem[1,'variance_Residual'] <- random_effects_var[3,4]
  
  # ICC
  result_tem[1,'ICC_cohort'] <- random_effects_var[1,4] / (random_effects_var[1,4] + random_effects_var[2,4] + random_effects_var[3,4])
  result_tem[1,'ICC_period'] <- random_effects_var[2,4] / (random_effects_var[1,4] + random_effects_var[2,4] + random_effects_var[3,4])
  
  # cohort_p & period_p
  significant <- ranova(fullmodel)
  result_tem[1,'cohort_p'] <- significant['(1 | cohort_cat)','Pr(>Chisq)']
  result_tem[1,'period_p'] <- significant['(1 | Year)','Pr(>Chisq)']
  
  # Cohort_effect1-9
  re <- ranef(fullmodel, condVar = TRUE)  
  
  cohort_re <- re$cohort_cat              
  cohort_postvar <- attr(re$cohort_cat, "postVar") 
  cohort_se <- sqrt(cohort_postvar[1,1,])
  cohort_summary <- data.frame(
    cohort = rownames(cohort_re),
    intercept = cohort_re[,1],
    lower = cohort_re[,1] - 1.96*cohort_se,
    upper = cohort_re[,1] + 1.96*cohort_se
  )
  result_tem[1,'Cohort1920-1929_intercept'] <- cohort_summary[1,'intercept']
  result_tem[1,'Cohort1920-1929_lower'] <- cohort_summary[1,'lower']
  result_tem[1,'Cohort1920-1929_upper'] <- cohort_summary[1,'upper']
  
  result_tem[1,'Cohort1930-1939_intercept'] <- cohort_summary[2,'intercept']
  result_tem[1,'Cohort1930-1939_lower'] <- cohort_summary[2,'lower']
  result_tem[1,'Cohort1930-1939_upper'] <- cohort_summary[2,'upper']
  
  result_tem[1,'Cohort1940-1949_intercept'] <- cohort_summary[3,'intercept']
  result_tem[1,'Cohort1940-1949_lower'] <- cohort_summary[3,'lower']
  result_tem[1,'Cohort1940-1949_upper'] <- cohort_summary[3,'upper']
  
  result_tem[1,'Cohort1950-1959_intercept'] <- cohort_summary[4,'intercept']
  result_tem[1,'Cohort1950-1959_lower'] <- cohort_summary[4,'lower']
  result_tem[1,'Cohort1950-1959_upper'] <- cohort_summary[4,'upper']
  
  result_tem[1,'Cohort1960-1969_intercept'] <- cohort_summary[5,'intercept']
  result_tem[1,'Cohort1960-1969_lower'] <- cohort_summary[5,'lower']
  result_tem[1,'Cohort1960-1969_upper'] <- cohort_summary[5,'upper']
  
  result_tem[1,'Cohort1970-1979_intercept'] <- cohort_summary[6,'intercept']
  result_tem[1,'Cohort1970-1979_lower'] <- cohort_summary[6,'lower']
  result_tem[1,'Cohort1970-1979_upper'] <- cohort_summary[6,'upper']
  
  result_tem[1,'Cohort1980-1989_intercept'] <- cohort_summary[7,'intercept']
  result_tem[1,'Cohort1980-1989_lower'] <- cohort_summary[7,'lower']
  result_tem[1,'Cohort1980-1989_upper'] <- cohort_summary[7,'upper']
  
  result_tem[1,'Cohort1990-1999_intercept'] <- cohort_summary[8,'intercept']
  result_tem[1,'Cohort1990-1999_lower'] <- cohort_summary[8,'lower']
  result_tem[1,'Cohort1990-1999_upper'] <- cohort_summary[8,'upper']
  
  result_tem[1,'Cohort2000-2009_intercept'] <- cohort_summary[9,'intercept']
  result_tem[1,'Cohort2000-2009_lower'] <- cohort_summary[9,'lower']
  result_tem[1,'Cohort2000-2009_upper'] <- cohort_summary[9,'upper']
  
  # Period_effect_dis 2003\2010\2021
  period_re <- re$Year             
  period_postvar <- attr(re$Year, "postVar") 
  period_se <- sqrt(period_postvar[1,1,])
  period_summary <- data.frame(
    period = rownames(period_re),
    intercept = period_re[,1],
    lower = period_re[,1] - 1.96*period_se,
    upper = period_re[,1] + 1.96*period_se
  )
  result_tem[1,'Period2003_intercept'] <- period_summary[1,'intercept']
  result_tem[1,'Period2003_lower'] <- period_summary[1,'lower']
  result_tem[1,'Period2003_upper'] <- period_summary[1,'upper']
  
  result_tem[1,'Period2010_intercept'] <- period_summary[2,'intercept']
  result_tem[1,'Period2010_lower'] <- period_summary[2,'lower']
  result_tem[1,'Period2010_upper'] <- period_summary[2,'upper']
  
  result_tem[1,'Period2021_intercept'] <- period_summary[3,'intercept']
  result_tem[1,'Period2021_lower'] <- period_summary[3,'lower']
  result_tem[1,'Period2021_upper'] <- period_summary[3,'upper']
  
  
  
  # age_effect & age_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  result_tem[1,'age_effect'] <- fixed_effect["AgeR", "Estimate"]
  result_tem[1,'age_p'] <- fixed_effect["AgeR", "Pr(>|t|)"]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect["AgeR2", "Estimate"]
  result_tem[1,'age2_p'] <- fixed_effect["AgeR2", "Pr(>|t|)"]
  
  # Cohort_effect_con & p
  result_tem[1,'Cohort_con_effect'] <- fixed_effect["birth_cent", "Estimate"]
  result_tem[1,'Cohort_con_p'] <- fixed_effect["birth_cent", "Pr(>|t|)"]
  
  # Control Variables
  result_tem[1,'female_effect'] <- fixed_effect["sex2", "Estimate"]
  result_tem[1,'female_p'] <- fixed_effect["sex2", "Pr(>|t|)"]
  result_tem[1,'education_effect'] <- fixed_effect['education',"Estimate"]
  result_tem[1,'education_p'] <- fixed_effect['education',"Pr(>|t|)"]
  result_tem[1,'household_income_effect'] <- fixed_effect['household_income',"Estimate"]
  result_tem[1,'household_income_p'] <- fixed_effect['household_income',"Pr(>|t|)"]
  result_tem[1,'happiness.2_effect'] <- fixed_effect["happiness2", "Estimate"]
  result_tem[1,'happiness.2_p'] <- fixed_effect["happiness2","Pr(>|t|)"]
  result_tem[1,'happiness.3_effect'] <- fixed_effect["happiness3","Estimate"]
  result_tem[1,'happiness.3_p'] <- fixed_effect["happiness3","Pr(>|t|)"]
  result_tem[1,'happiness.4_effect'] <- fixed_effect["happiness4","Estimate"]
  result_tem[1,'happiness.4_p'] <- fixed_effect["happiness4","Pr(>|t|)"]
  result_tem[1,'happiness.5_effect'] <- fixed_effect["happiness5","Estimate"]
  result_tem[1,'happiness.5_p'] <- fixed_effect["happiness5","Pr(>|t|)"]
  
  result_CNEP <- rbind(result_CNEP, result_tem)
}
# write.csv(result_CNEP[,1:65], file = 'result_CNEP.csv', row.names = FALSE)
# write.csv(result_CNEP[,10:36], file = 'result_cohort_random.csv', row.names = FALSE)
# write.csv(result_CNEP[,37:45], file = 'result_period_random.csv', row.names = FALSE)
# These two combined form Supplementary Table 1.
# For convenience in plotting, store them in two separate files here, and do not consider R2 when plotting.

### 2.1.4 Table 1 ====
result_CNEP <- read.csv('result_CNEP.csv')
extract_and_summarize <- function(df) {
  # 指定要提取的列名
  selected_cols <- c(
    "variance_cohort", "variance_period", "variance_Residual",
    "ICC_cohort", "ICC_period",
    "cohort_p", "period_p",
    "age_effect", "age_p",
    "age2_effect", "age2_p",
    "Cohort_con_effect", "Cohort_con_p",
    "female_effect", "female_p",
    "education_effect", "education_p",
    "household_income_effect", "household_income_p",
    "happiness.2_effect", "happiness.2_p",
    "happiness.3_effect", "happiness.3_p",
    "happiness.4_effect", "happiness.4_p",
    "happiness.5_effect", "happiness.5_p"
  )
  # 检查哪些列实际存在于数据框中
  existing_cols <- selected_cols[selected_cols %in% names(df)]
  if(length(existing_cols) == 0) {
    stop("没有找到指定的列")
  }
  if(length(existing_cols) < length(selected_cols)) {
    missing_cols <- setdiff(selected_cols, existing_cols)
    warning(paste("以下列不存在：", paste(missing_cols, collapse = ", ")))
  }
  # 提取数据并计算统计量
  result <- data.frame(
    Variable = existing_cols,
    Mean = sapply(df[existing_cols], mean, na.rm = TRUE),
    Variance = sapply(df[existing_cols], var, na.rm = TRUE)
  )
  # 设置行名
  rownames(result) <- NULL
  return(result)
}

# 使用函数
summary_stats <- extract_and_summarize(result_CNEP)
print(summary_stats)

write.csv(summary_stats[,1:3], file = 'Table1.csv', row.names = FALSE)


### 2.1.5 Figure 2 ====
df_cohort_random <- read.csv('result_cohort_random.csv')

intercept_cols <- grep("_intercept$", names(df_cohort_random), value = TRUE)
df_long_cohort_random <- df_cohort_random %>%
  select(all_of(intercept_cols)) %>%
  pivot_longer(cols = everything(), 
               names_to = "cohort", 
               values_to = "intercept") %>%
  mutate(period = gsub("_intercept", "", cohort))
df_long_cohort_random <- df_long_cohort_random %>%
  mutate(cohort = recode(cohort,
                         `Cohort1920.1929_intercept` = "1920-1929",
                         `Cohort1930.1939_intercept` = "1930-1939",
                         `Cohort1940.1949_intercept` = "1940-1949",
                         `Cohort1950.1959_intercept` = "1950-1959",
                         `Cohort1960.1969_intercept` = "1960-1969",
                         `Cohort1970.1979_intercept` = "1970-1979",
                         `Cohort1980.1989_intercept` = "1980-1989",
                         `Cohort1990.1999_intercept` = "1990-1999",
                         `Cohort2000.2009_intercept` = "2000-2009"))

p_cohort_random <- df_long_cohort_random %>%
  ggplot( aes(x=cohort, y=intercept, fill=cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8, alpha = 0.7) +
  geom_violin(width=1.1,color='white', adjust = 5,trim = FALSE) +
  scale_fill_manual(values = c(  "#66c2a5",  # 1. 浅绿色
                                 "#fc8d62",  # 2. 浅橙色
                                 "#8da0cb",  # 3. 浅蓝紫色
                                 "#e78ac3",  # 4. 浅粉色
                                 "#a6d854",  # 5. 浅黄绿色
                                 "#ffd92f",  # 6. 浅黄色
                                 "#e5c494",  # 7. 浅米色
                                 "#b3b3b3",  # 8. 浅灰色
                                 "#80b1d3")) +
  geom_boxplot(width=0.35, color="black",fill='white', alpha=0.2, outlier.shape = NA) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),  
    axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 8)),
    axis.title.y = element_text(size = 14, color = "black"),  
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1, vjust = 1, margin = margin(t = 5)), 
    axis.text.y = element_text(size = 10, color = "black"), 
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "black"), 
    axis.ticks = element_line(color = "black")  
  ) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  ggtitle("Cohort Level Intercept") + 
  xlab("Cohort") + 
  ylab("Intercept") 
p_cohort_random


df_period_random <- read.csv('result_period_random.csv')

intercept_cols <- grep("_intercept$", names(df_period_random), value = TRUE)
df_long_period_random <- df_period_random %>%
  select(all_of(intercept_cols)) %>%
  pivot_longer(cols = everything(), 
               names_to = "period", 
               values_to = "intercept") %>%
  mutate(period = gsub("_intercept", "", period))
df_long_period_random <- df_long_period_random %>%
  mutate(period = recode(period,
                         `Period2003` = "2003",
                         `Period2010` = "2010",
                         `Period2021` = "2021"))

p_period_random <- df_long_period_random %>%
  ggplot( aes(x=period, y=intercept, fill=period)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8, alpha = 0.7) +
  geom_violin(width=1,color='white', adjust = 5,trim = FALSE) +
  scale_fill_manual(values = c(  "#66c2a5",  # 1. 浅绿色
                                 "#fc8d62",  # 2. 浅橙色
                                 "#8da0cb"  # 3. 浅蓝紫色
  )) +
  geom_boxplot(width=0.2, color="black",fill='white', alpha=0.2, outlier.shape = NA) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),  
    axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 8)),
    axis.title.y = element_text(size = 14, color = "black"),  
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1, vjust = 1, margin = margin(t = 5)), 
    axis.text.y = element_text(size = 10, color = "black"), 
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "black"), 
    axis.ticks = element_line(color = "black")  
  ) +
  coord_cartesian(ylim = c(-0.8, 0.8)) +
  ggtitle("Period Level Intercept") + 
  xlab("Period") + 
  ylab("Intercept") 
p_period_random


combined_plot <- p_period_random + p_cohort_random + 
  plot_layout(ncol = 2) +  # 两列并排
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

# 保存图像
print(combined_plot)
ggsave("Figure2.png", combined_plot, width = 10, height = 5, dpi = 600)


### 2.1.6 Supplementary Fig.1 ====
df_period_random <- read.csv('result_period_random.csv')
df_period_random <- df_period_random %>%
  mutate(model = row_number())
df_long <- df_period_random %>%
  pivot_longer(
    cols = -model,
    names_to = c("period","type"),
    names_pattern = "Period(.*)_(.*)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  )
df_long$period <- factor(
  df_long$period,
  levels = c(
    "2003","2010","2021"
  )
)
S1 <- ggplot(df_long,
            aes(x = period,
                y = intercept,
                group = model)) +
  geom_point(size = 1.2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15
  ) +
  geom_line() +
  facet_wrap(~model, ncol = 10) +
  theme_bw() +
  labs(
    x = "Period",
    y = "Random Intercept",
    title = "Random Period Effects Across 125 Models",
  ) +
  theme(
    axis.text.x = element_text(size = 6),
    strip.text = element_text(size = 6),
    plot.caption = element_text(size = 10, hjust = 0)
  )
ggsave(
  "Supplementary Fig.1.png",
  S1,
  width = 20,
  height = 15,
  dpi = 300
)

### 2.1.5 Supplementary Fig.2 ====
df_cohort_random <- read.csv('result_cohort_random.csv')
df_cohort_random <- df_cohort_random %>%
  mutate(model = row_number())
df_long <- df_cohort_random %>%
  pivot_longer(
    cols = -model,
    names_to = c("cohort","type"),
    names_pattern = "Cohort(.*)_(.*)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  )
df_long$cohort <- factor(
  df_long$cohort,
  levels = c(
    "1920.1929","1930.1939","1940.1949","1950.1959",
    "1960.1969","1970.1979","1980.1989","1990.1999","2000.2009"
  )
)
cohort_levels <- c(
  "1920.1929","1930.1939","1940.1949","1950.1959",
  "1960.1969","1970.1979","1980.1989","1990.1999","2000.2009"
)

df_long <- df_long %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    cohort_id = as.numeric(cohort)
  )
cohort_legend <- paste0(
  1:9, ": ",
  gsub("\\.", "-", cohort_levels),
  collapse = "   "
)
S2 <- ggplot(df_long,
            aes(x = cohort_id,
                y = intercept,
                group = model)) +
  
  geom_point(size = 1.2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15
  ) +
  geom_line() +
  facet_wrap(~model, ncol = 10) +
  scale_x_continuous(
    breaks = 1:9,
    labels = 1:9
  ) +
  theme_bw() +
  labs(
    x = "Cohort ID",
    y = "Random Intercept",
    title = "Random Cohort Effects Across 125 Models",
    caption = cohort_legend
  ) +
  theme(
    axis.text.x = element_text(size = 6),
    strip.text = element_text(size = 6),
    plot.caption = element_text(size = 10, hjust = 0)
  )
ggsave(
  "Supplementary Fig.2.png",
  S2,
  width = 20,
  height = 15,
  dpi = 300
)

### 2.1.6 Supplementary Fig.3 ====
file_list <- list.files(pattern = "^.*\\.csv$")
data_list <- lapply(file_list, read.csv) 

# Reverse the scores for reverse-coded questions.
col_names_list_reset <- c( "NEP2", 'NEP4','NEP6','NEP8','NEP10','NEP12','NEP14')
replace_values_NEP <- function(data, col_name) {
  data[[col_name]] <- ifelse(data[[col_name]] == 1, 5,
                             ifelse(data[[col_name]] == 2, 4,
                                    ifelse(data[[col_name]] == 4, 2,
                                           ifelse(data[[col_name]] == 5, 1, data[[col_name]]))))
  return(data)
}

data_list_modified <- lapply(data_list, function(data) {
  for (col_name in col_names_list_reset) {
    data <- replace_values_NEP(data, col_name)
  }
  return(data)
})

# Calculate the age used in the analysis (AgeR:divide by ten and subtract 2), as well as its square, 
# classify the generations(cohort) and calculate CNEP.
data_list_modified <- lapply(data_list_modified, function(data) {
  data <- data %>%
    mutate(CNEP = NEP1 + NEP3 + NEP5 + NEP7 + NEP8 + NEP9 + NEP10 + NEP11 + NEP13 + NEP15)
  data$birth <- data$Year - data$age
  data$AgeR = data$age/10 - 2
  data$AgeR2 = data$AgeR^2
  data <- data %>% mutate(cohort_cat = case_when(
    birth >= 1920 & birth <= 1929 ~ 1,
    birth >= 1930 & birth <= 1939 ~ 2,
    birth >= 1940 & birth <= 1949 ~ 3,
    birth >= 1950 & birth <= 1959 ~ 4,
    birth >= 1960 & birth <= 1969 ~ 5,
    birth >= 1970 & birth <= 1979 ~ 6,
    birth >= 1980 & birth <= 1989 ~ 7,
    birth >= 1990 & birth <= 1999 ~ 8,
    birth >= 2000 & birth <= 2009 ~ 9,
    TRUE ~ NA_integer_  # Define the samples that do not belong to these generations as NA.
  ))
  data <- data[!is.na(data$cohort_cat), ] 
  data$cohort_cat <- as.factor(data$cohort_cat)
  data <- data %>% mutate_at(vars(sex,Year,happiness), as.factor)
  data$birth_cent <- data$birth - mean(data$birth)
  return(data)
})

matstrict <- cem("Year", datalist=data_list_modified, drop=c('NEP1','NEP2','NEP3','NEP4','NEP5','NEP6','NEP7','NEP8','NEP9',
                                                             'NEP10','NEP11','NEP12','NEP13','NEP14','NEP15',
                                                             'age','birth','Year','CNEP',
                                                             'AgeR','AgeR2','cohort_cat','birth_cent'))
matstrict
matched_data_list <- list()
for (i in seq_along(data_list)) {
  match_col <- paste0("match", i)
  matched_data_list[[i]] <- data_list_modified[[i]][matstrict[[match_col]]$matched, ]
}

plot_model <- function(data, index) {
  model <- lmer(
    CNEP ~ AgeR + I(AgeR^2) + birth_cent + sex + happiness +
      education + household_income +
      (1 | Year) + (1 | cohort_cat),
    data = data
  )
  pred <- ggeffect(model, terms = "AgeR [all]")
  p <- ggplot(pred, aes(x = x, y = predicted)) +
    geom_line(linewidth = 1, color = "blue") +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.2,
      fill = "blue"
    ) +
    labs(
      x = "AgeR",
      y = "Predicted CNEP",
      title = paste("Dataset", index)
    ) +
    theme_minimal()
  return(p)
}

plots <- map2(
  matched_data_list,
  seq_along(matched_data_list),
  plot_model
)
p_S3 <- wrap_plots(plots) +
  plot_layout(ncol = 10)

ggsave(
  "Supplementary Fig.3.png",
  p_S3,
  width = 20,
  height = 25,
  units = "in",
  dpi = 300
)


# 3. Sensitivity Analysis ====
## 3.1 Supplementary data files 2  ====
file_list <- list.files(pattern = "^.*\\.csv$")
data_list <- lapply(file_list, read.csv) 

# Reverse the scores for reverse-coded questions.
col_names_list_reset <- c( "NEP2", 'NEP4','NEP6','NEP8','NEP10','NEP12','NEP14')
replace_values_NEP <- function(data, col_name) {
  data[[col_name]] <- ifelse(data[[col_name]] == 1, 5,
                             ifelse(data[[col_name]] == 2, 4,
                                    ifelse(data[[col_name]] == 4, 2,
                                           ifelse(data[[col_name]] == 5, 1, data[[col_name]]))))
  return(data)
}

data_list_modified <- lapply(data_list, function(data) {
  for (col_name in col_names_list_reset) {
    data <- replace_values_NEP(data, col_name)
  }
  return(data)
})

# Calculate the age used in the analysis (AgeR:divide by ten and subtract 2), as well as its square, 
# classify the generations(cohort) and calculate CNEP.
data_list_modified <- lapply(data_list_modified, function(data) {
  data <- data %>%
    mutate(CNEP = NEP1 + NEP3 + NEP5 + NEP7 + NEP8 + NEP9 + NEP10 + NEP11 + NEP13 + NEP15)
  data$birth <- data$Year - data$age
  data$AgeR = data$age/10 - 2
  data$AgeR2 = data$AgeR^2
  data <- data %>% mutate(cohort_cat = case_when(
    birth >= 1930 & birth <= 1939 ~ 1,
    birth >= 1940 & birth <= 1949 ~ 2,
    birth >= 1950 & birth <= 1959 ~ 3,
    birth >= 1960 & birth <= 1969 ~ 4,
    birth >= 1970 & birth <= 1979 ~ 5,
    birth >= 1980 & birth <= 1989 ~ 6,
    birth >= 1990 & birth <= 1999 ~ 7,
    TRUE ~ NA_integer_  # Define the samples that do not belong to these generations as NA.
  ))
  data <- data[!is.na(data$cohort_cat), ] 
  data$cohort_cat <- as.factor(data$cohort_cat)
  data <- data %>% mutate_at(vars(sex,Year,happiness), as.factor)
  data$birth_cent <- data$birth - mean(data$birth)
  return(data)
})

###  Perform Coarsened Exact Matching 
matstrict <- cem("Year", datalist=data_list_modified, drop=c('NEP1','NEP2','NEP3','NEP4','NEP5','NEP6','NEP7','NEP8','NEP9',
                                                             'NEP10','NEP11','NEP12','NEP13','NEP14','NEP15',
                                                             'age','birth','Year','CNEP',
                                                             'AgeR','AgeR2','cohort_cat','birth_cent'))
matstrict
matched_data_list <- list()
for (i in seq_along(data_list)) {
  match_col <- paste0("match", i)
  matched_data_list[[i]] <- data_list_modified[[i]][matstrict[[match_col]]$matched, ]
}

###  Construct the HAPC model for each matched dataset 
result_CNEP_S2 <- data.frame()
for (i in 1:length(matched_data_list)) {
  data <- matched_data_list[[i]]
  result_tem <- data.frame() 
  
  fullmodel <- lmer(CNEP ~ AgeR + AgeR2 + birth_cent + sex + happiness + education + household_income
                    + (1 | Year) + (1 | cohort_cat), data = data)
  summary_fullmodel <- summary(fullmodel)
  
  # Model information
  result_tem[1,'Model'] <- paste0("data", i)
  result_tem[1,'samplesize'] <- nrow(data)
  
  # variance_cohort & variance_period & variance_Residual
  random_effects_var <- as.data.frame(VarCorr(fullmodel))
  result_tem[1,'variance_cohort'] <- random_effects_var[1,4]
  result_tem[1,'variance_period'] <- random_effects_var[2,4]
  result_tem[1,'variance_Residual'] <- random_effects_var[3,4]
  
  # ICC
  result_tem[1,'ICC_cohort'] <- random_effects_var[1,4] / (random_effects_var[1,4] + random_effects_var[2,4] + random_effects_var[3,4])
  result_tem[1,'ICC_period'] <- random_effects_var[2,4] / (random_effects_var[1,4] + random_effects_var[2,4] + random_effects_var[3,4])
  
  # cohort_p & period_p
  significant <- ranova(fullmodel)
  result_tem[1,'cohort_p'] <- significant['(1 | cohort_cat)','Pr(>Chisq)']
  result_tem[1,'period_p'] <- significant['(1 | Year)','Pr(>Chisq)']
  
  # Cohort_effect1-7
  re <- ranef(fullmodel, condVar = TRUE)  
  
  cohort_re <- re$cohort_cat              
  cohort_postvar <- attr(re$cohort_cat, "postVar") 
  cohort_se <- sqrt(cohort_postvar[1,1,])
  cohort_summary <- data.frame(
    cohort = rownames(cohort_re),
    intercept = cohort_re[,1],
    lower = cohort_re[,1] - 1.96*cohort_se,
    upper = cohort_re[,1] + 1.96*cohort_se
  )
  
  result_tem[1,'Cohort1930-1939_intercept'] <- cohort_summary[1,'intercept']
  result_tem[1,'Cohort1930-1939_lower'] <- cohort_summary[1,'lower']
  result_tem[1,'Cohort1930-1939_upper'] <- cohort_summary[1,'upper']
  
  result_tem[1,'Cohort1940-1949_intercept'] <- cohort_summary[2,'intercept']
  result_tem[1,'Cohort1940-1949_lower'] <- cohort_summary[2,'lower']
  result_tem[1,'Cohort1940-1949_upper'] <- cohort_summary[2,'upper']
  
  result_tem[1,'Cohort1950-1959_intercept'] <- cohort_summary[3,'intercept']
  result_tem[1,'Cohort1950-1959_lower'] <- cohort_summary[3,'lower']
  result_tem[1,'Cohort1950-1959_upper'] <- cohort_summary[3,'upper']
  
  result_tem[1,'Cohort1960-1969_intercept'] <- cohort_summary[4,'intercept']
  result_tem[1,'Cohort1960-1969_lower'] <- cohort_summary[4,'lower']
  result_tem[1,'Cohort1960-1969_upper'] <- cohort_summary[4,'upper']
  
  result_tem[1,'Cohort1970-1979_intercept'] <- cohort_summary[5,'intercept']
  result_tem[1,'Cohort1970-1979_lower'] <- cohort_summary[5,'lower']
  result_tem[1,'Cohort1970-1979_upper'] <- cohort_summary[5,'upper']
  
  result_tem[1,'Cohort1980-1989_intercept'] <- cohort_summary[6,'intercept']
  result_tem[1,'Cohort1980-1989_lower'] <- cohort_summary[6,'lower']
  result_tem[1,'Cohort1980-1989_upper'] <- cohort_summary[6,'upper']
  
  result_tem[1,'Cohort1990-1999_intercept'] <- cohort_summary[7,'intercept']
  result_tem[1,'Cohort1990-1999_lower'] <- cohort_summary[7,'lower']
  result_tem[1,'Cohort1990-1999_upper'] <- cohort_summary[7,'upper']
  
  # Period_effect_dis 2003\2010\2021
  period_re <- re$Year             
  period_postvar <- attr(re$Year, "postVar") 
  period_se <- sqrt(period_postvar[1,1,])
  period_summary <- data.frame(
    period = rownames(period_re),
    intercept = period_re[,1],
    lower = period_re[,1] - 1.96*period_se,
    upper = period_re[,1] + 1.96*period_se
  )
  result_tem[1,'Period2003_intercept'] <- period_summary[1,'intercept']
  result_tem[1,'Period2003_lower'] <- period_summary[1,'lower']
  result_tem[1,'Period2003_upper'] <- period_summary[1,'upper']
  
  result_tem[1,'Period2010_intercept'] <- period_summary[2,'intercept']
  result_tem[1,'Period2010_lower'] <- period_summary[2,'lower']
  result_tem[1,'Period2010_upper'] <- period_summary[2,'upper']
  
  result_tem[1,'Period2021_intercept'] <- period_summary[3,'intercept']
  result_tem[1,'Period2021_lower'] <- period_summary[3,'lower']
  result_tem[1,'Period2021_upper'] <- period_summary[3,'upper']
  
  
  
  # age_effect & age_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  result_tem[1,'age_effect'] <- fixed_effect["AgeR", "Estimate"]
  result_tem[1,'age_p'] <- fixed_effect["AgeR", "Pr(>|t|)"]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect["AgeR2", "Estimate"]
  result_tem[1,'age2_p'] <- fixed_effect["AgeR2", "Pr(>|t|)"]
  
  # Cohort_effect_con & p
  result_tem[1,'Cohort_con_effect'] <- fixed_effect["birth_cent", "Estimate"]
  result_tem[1,'Cohort_con_p'] <- fixed_effect["birth_cent", "Pr(>|t|)"]
  
  # Control Variables
  result_tem[1,'female_effect'] <- fixed_effect["sex2", "Estimate"]
  result_tem[1,'female_p'] <- fixed_effect["sex2", "Pr(>|t|)"]
  result_tem[1,'education_effect'] <- fixed_effect['education',"Estimate"]
  result_tem[1,'education_p'] <- fixed_effect['education',"Pr(>|t|)"]
  result_tem[1,'household_income_effect'] <- fixed_effect['household_income',"Estimate"]
  result_tem[1,'household_income_p'] <- fixed_effect['household_income',"Pr(>|t|)"]
  result_tem[1,'happiness.2_effect'] <- fixed_effect["happiness2", "Estimate"]
  result_tem[1,'happiness.2_p'] <- fixed_effect["happiness2","Pr(>|t|)"]
  result_tem[1,'happiness.3_effect'] <- fixed_effect["happiness3","Estimate"]
  result_tem[1,'happiness.3_p'] <- fixed_effect["happiness3","Pr(>|t|)"]
  result_tem[1,'happiness.4_effect'] <- fixed_effect["happiness4","Estimate"]
  result_tem[1,'happiness.4_p'] <- fixed_effect["happiness4","Pr(>|t|)"]
  result_tem[1,'happiness.5_effect'] <- fixed_effect["happiness5","Estimate"]
  result_tem[1,'happiness.5_p'] <- fixed_effect["happiness5","Pr(>|t|)"]
  
  result_CNEP_S2 <- rbind(result_CNEP_S2, result_tem)
}

## Supplementary data files 2
write.csv(result_CNEP_S2[,1:59], file = 'Supplementary data files 2.csv', row.names = FALSE)


## 3.2 Supplementary data files 3  ====
file_list <- list.files(pattern = "^.*\\.csv$")
data_list <- lapply(file_list, read.csv) 

# Reverse the scores for reverse-coded questions.
col_names_list_reset <- c( "NEP2", 'NEP4','NEP6','NEP8','NEP10','NEP12','NEP14')
replace_values_NEP <- function(data, col_name) {
  data[[col_name]] <- ifelse(data[[col_name]] == 1, 5,
                             ifelse(data[[col_name]] == 2, 4,
                                    ifelse(data[[col_name]] == 4, 2,
                                           ifelse(data[[col_name]] == 5, 1, data[[col_name]]))))
  return(data)
}

data_list_modified <- lapply(data_list, function(data) {
  for (col_name in col_names_list_reset) {
    data <- replace_values_NEP(data, col_name)
  }
  return(data)
})

# Calculate the age used in the analysis (AgeR:divide by ten and subtract 2), as well as its square, 
# classify the generations(cohort) and calculate CNEP.
data_list_modified <- lapply(data_list_modified, function(data) {
  data <- data %>%
    mutate(CNEP = NEP1 + NEP3 + NEP5 + NEP7 + NEP8 + NEP9 + NEP10 + NEP11 + NEP13 + NEP15)
  data$birth <- data$Year - data$age
  data$AgeR = data$age/10 - 2
  data$AgeR2 = data$AgeR^2
  data <- data %>% mutate(cohort_cat = case_when(
    birth >= 1920 & birth <= 1929 ~ 1,
    birth >= 1930 & birth <= 1939 ~ 2,
    birth >= 1940 & birth <= 1949 ~ 3,
    birth >= 1950 & birth <= 1959 ~ 4,
    birth >= 1960 & birth <= 1969 ~ 5,
    birth >= 1970 & birth <= 1979 ~ 6,
    birth >= 1980 & birth <= 1989 ~ 7,
    birth >= 1990 & birth <= 1999 ~ 8,
    birth >= 2000 & birth <= 2009 ~ 9,
    TRUE ~ NA_integer_  # Define the samples that do not belong to these generations as NA.
  ))
  data <- data[!is.na(data$cohort_cat), ] 
  data$cohort_cat <- as.factor(data$cohort_cat)
  data <- data %>% mutate_at(vars(sex,Year,happiness), as.factor)
  data$birth_cent <- data$birth - mean(data$birth)
  return(data)
})

### Perform Coarsened Exact Matching 
matstrict <- cem("Year", datalist=data_list_modified, drop=c('NEP1','NEP2','NEP3','NEP4','NEP5','NEP6','NEP7','NEP8','NEP9',
                                                             'NEP10','NEP11','NEP12','NEP13','NEP14','NEP15',
                                                             'age','birth','Year','CNEP',
                                                             'AgeR','AgeR2','cohort_cat','birth_cent'))
matstrict
matched_data_list <- list()
for (i in seq_along(data_list)) {
  match_col <- paste0("match", i)
  matched_data_list[[i]] <- data_list_modified[[i]][matstrict[[match_col]]$matched, ]
}

### Construct the HAPC model for each matched dataset
result_CNEP_S3 <- data.frame()
for (i in 1:length(matched_data_list)) {
  data <- matched_data_list[[i]]
  result_tem <- data.frame() 
  
  fullmodel <- lmer(CNEP ~ AgeR + AgeR2 + birth_cent + 
                    + (1 | Year) + (1 | cohort_cat), data = data)
  summary_fullmodel <- summary(fullmodel)
  
  # Model information
  result_tem[1,'Model'] <- paste0("data", i)
  result_tem[1,'samplesize'] <- nrow(data)
  
  # variance_cohort & variance_period & variance_Residual
  random_effects_var <- as.data.frame(VarCorr(fullmodel))
  result_tem[1,'variance_cohort'] <- random_effects_var[1,4]
  result_tem[1,'variance_period'] <- random_effects_var[2,4]
  result_tem[1,'variance_Residual'] <- random_effects_var[3,4]
  
  # ICC
  result_tem[1,'ICC_cohort'] <- random_effects_var[1,4] / (random_effects_var[1,4] + random_effects_var[2,4] + random_effects_var[3,4])
  result_tem[1,'ICC_period'] <- random_effects_var[2,4] / (random_effects_var[1,4] + random_effects_var[2,4] + random_effects_var[3,4])
  
  # cohort_p & period_p
  significant <- ranova(fullmodel)
  result_tem[1,'cohort_p'] <- significant['(1 | cohort_cat)','Pr(>Chisq)']
  result_tem[1,'period_p'] <- significant['(1 | Year)','Pr(>Chisq)']
  
  # Cohort_effect1-9
  re <- ranef(fullmodel, condVar = TRUE)  
  
  cohort_re <- re$cohort_cat              
  cohort_postvar <- attr(re$cohort_cat, "postVar") 
  cohort_se <- sqrt(cohort_postvar[1,1,])
  cohort_summary <- data.frame(
    cohort = rownames(cohort_re),
    intercept = cohort_re[,1],
    lower = cohort_re[,1] - 1.96*cohort_se,
    upper = cohort_re[,1] + 1.96*cohort_se
  )
  result_tem[1,'Cohort1920-1929_intercept'] <- cohort_summary[1,'intercept']
  result_tem[1,'Cohort1920-1929_lower'] <- cohort_summary[1,'lower']
  result_tem[1,'Cohort1920-1929_upper'] <- cohort_summary[1,'upper']
  
  result_tem[1,'Cohort1930-1939_intercept'] <- cohort_summary[2,'intercept']
  result_tem[1,'Cohort1930-1939_lower'] <- cohort_summary[2,'lower']
  result_tem[1,'Cohort1930-1939_upper'] <- cohort_summary[2,'upper']
  
  result_tem[1,'Cohort1940-1949_intercept'] <- cohort_summary[3,'intercept']
  result_tem[1,'Cohort1940-1949_lower'] <- cohort_summary[3,'lower']
  result_tem[1,'Cohort1940-1949_upper'] <- cohort_summary[3,'upper']
  
  result_tem[1,'Cohort1950-1959_intercept'] <- cohort_summary[4,'intercept']
  result_tem[1,'Cohort1950-1959_lower'] <- cohort_summary[4,'lower']
  result_tem[1,'Cohort1950-1959_upper'] <- cohort_summary[4,'upper']
  
  result_tem[1,'Cohort1960-1969_intercept'] <- cohort_summary[5,'intercept']
  result_tem[1,'Cohort1960-1969_lower'] <- cohort_summary[5,'lower']
  result_tem[1,'Cohort1960-1969_upper'] <- cohort_summary[5,'upper']
  
  result_tem[1,'Cohort1970-1979_intercept'] <- cohort_summary[6,'intercept']
  result_tem[1,'Cohort1970-1979_lower'] <- cohort_summary[6,'lower']
  result_tem[1,'Cohort1970-1979_upper'] <- cohort_summary[6,'upper']
  
  result_tem[1,'Cohort1980-1989_intercept'] <- cohort_summary[7,'intercept']
  result_tem[1,'Cohort1980-1989_lower'] <- cohort_summary[7,'lower']
  result_tem[1,'Cohort1980-1989_upper'] <- cohort_summary[7,'upper']
  
  result_tem[1,'Cohort1990-1999_intercept'] <- cohort_summary[8,'intercept']
  result_tem[1,'Cohort1990-1999_lower'] <- cohort_summary[8,'lower']
  result_tem[1,'Cohort1990-1999_upper'] <- cohort_summary[8,'upper']
  
  result_tem[1,'Cohort2000-2009_intercept'] <- cohort_summary[9,'intercept']
  result_tem[1,'Cohort2000-2009_lower'] <- cohort_summary[9,'lower']
  result_tem[1,'Cohort2000-2009_upper'] <- cohort_summary[9,'upper']
  
  # Period_effect_dis 2003\2010\2021
  period_re <- re$Year             
  period_postvar <- attr(re$Year, "postVar") 
  period_se <- sqrt(period_postvar[1,1,])
  period_summary <- data.frame(
    period = rownames(period_re),
    intercept = period_re[,1],
    lower = period_re[,1] - 1.96*period_se,
    upper = period_re[,1] + 1.96*period_se
  )
  result_tem[1,'Period2003_intercept'] <- period_summary[1,'intercept']
  result_tem[1,'Period2003_lower'] <- period_summary[1,'lower']
  result_tem[1,'Period2003_upper'] <- period_summary[1,'upper']
  
  result_tem[1,'Period2010_intercept'] <- period_summary[2,'intercept']
  result_tem[1,'Period2010_lower'] <- period_summary[2,'lower']
  result_tem[1,'Period2010_upper'] <- period_summary[2,'upper']
  
  result_tem[1,'Period2021_intercept'] <- period_summary[3,'intercept']
  result_tem[1,'Period2021_lower'] <- period_summary[3,'lower']
  result_tem[1,'Period2021_upper'] <- period_summary[3,'upper']
  
  # age_effect & age_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  result_tem[1,'age_effect'] <- fixed_effect["AgeR", "Estimate"]
  result_tem[1,'age_p'] <- fixed_effect["AgeR", "Pr(>|t|)"]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect["AgeR2", "Estimate"]
  result_tem[1,'age2_p'] <- fixed_effect["AgeR2", "Pr(>|t|)"]
  
  # Cohort_effect_con & p
  result_tem[1,'Cohort_con_effect'] <- fixed_effect["birth_cent", "Estimate"]
  result_tem[1,'Cohort_con_p'] <- fixed_effect["birth_cent", "Pr(>|t|)"]
  
  
  result_CNEP_S3 <- rbind(result_CNEP_S3, result_tem)
}
## Supplementary data files 3
write.csv(result_CNEP_S3[,1:51], file = 'Supplementary data files 3.csv', row.names = FALSE)
