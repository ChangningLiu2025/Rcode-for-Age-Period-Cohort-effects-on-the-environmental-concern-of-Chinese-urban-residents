# Author：Changning Liu  2025
# e-mail:liuchn6@mail2.sysu.edu.cn 
# Age-Period-Cohort effects on the environmental concern of Chinese urban residents

library(cem)
library(dplyr)
library(lme4)
library(lmerTest)
library(mice)
library(tidyr)
library(MuMIn)
library(ggplot2)
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
                               'lifefeel','k14','educ','hhyrinc')]
# Rename column names
names(CGSS2003_filter)[names(CGSS2003_filter) == "lifefeel"] <- "happiness"
names(CGSS2003_filter)[names(CGSS2003_filter) == "k14"] <- "gov_satisfaction"
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
CGSS2003_filter$gov_satisfaction[CGSS2003_filter$gov_satisfaction %in% c(-3, -2, -1)] <- NA 
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

# Define government efforts based on the data as follows:  
#   - "Overemphasizing economic development while neglecting environmental protection" = 1  
# - "Insufficient attention and inadequate investment in environmental protection" = 1  
# - "Unclear" = 2  
# - "Made efforts but with unsatisfactory results" = 2  
# - "Made significant efforts with some achievements" = 3  
# - "Achieved great success" = 3
CGSS2003_filter$gov_satisfaction <- ifelse(CGSS2003_filter$gov_satisfaction %in% c(1,2), 1, 
                                           ifelse(CGSS2003_filter$gov_satisfaction %in% c(3,6), 2, 
                                                  ifelse(CGSS2003_filter$gov_satisfaction %in% c(4,5), 3,CGSS2003_filter$gov_satisfaction)))

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
# - **Ordered categorical variables**: NEP, government satisfaction, subjective well-being
CGSS2003_filter$age <- as.numeric(CGSS2003_filter$age)
CGSS2003_filter$household_income <- as.numeric(CGSS2003_filter$household_income)
CGSS2003_filter$household_income <- scale(CGSS2003_filter$household_income) # Standardize household income to ensure consistency in variable scales.
CGSS2003_filter$household_income <- as.numeric(CGSS2003_filter$household_income)
CGSS2003_filter$education <- as.numeric(CGSS2003_filter$education)
CGSS2003_filter <- CGSS2003_filter %>% mutate_at(vars(sex), as.factor)
NEP_order <- c("1", "2", "3",'4','5')
CGSS2003_filter <- CGSS2003_filter %>% mutate_at(vars(NEP1,NEP2,NEP3,NEP4,NEP5,NEP6,NEP7,NEP8,NEP9,NEP10,NEP11,NEP12,NEP13,NEP14,NEP15), 
                                                 list(~ factor(., levels = NEP_order, ordered = TRUE)))
gov_order <- c("1", "2", "3")
CGSS2003_filter <- CGSS2003_filter %>% mutate(gov_satisfaction = factor(gov_satisfaction, levels = gov_order, ordered = TRUE))
happy_order <- c("1", "2", "3",'4','5')
CGSS2003_filter <- CGSS2003_filter %>% mutate(happiness = factor(happiness, levels = happy_order, ordered = TRUE))

# Sort the data according to the specified order.
custom_order <- c("NEP1", "NEP2", "NEP3",'NEP4','NEP5','NEP6','NEP7','NEP8','NEP9','NEP10','NEP11','NEP12','NEP13','NEP14',
                  'NEP15','sex','age','happiness','gov_satisfaction','education','household_income','birth')# ,rental,marriage
CGSS2003_filter <- CGSS2003_filter %>% select(all_of(custom_order))
CGSS2003_filter$Year <- 2003

# Write the data before imputation (i.e., CGSS2003 in Supplementary Table 7).
# write.csv(CGSS2003_filter, file = 'CGSS2003(插补前).csv', row.names = FALSE)

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
                               'l2513','l2514','l2515','a36','l16c','a7a','a62','s5')]

# Rename column names
names(CGSS2010_filter)[names(CGSS2010_filter) == "a2"] <- "sex"
names(CGSS2010_filter)[names(CGSS2010_filter) == "a3a"] <- "birth"
names(CGSS2010_filter)[names(CGSS2010_filter) == "a36"] <- "happiness"
names(CGSS2010_filter)[names(CGSS2010_filter) == "l16c"] <- "gov_satisfaction"
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
CGSS2010_filter$gov_satisfaction[CGSS2010_filter$gov_satisfaction %in% c(-3, -2, -1)] <- NA 
CGSS2010_filter$household_income[CGSS2010_filter$household_income %in% c(9999997, 9999998, 9999999)] <- NA 

# Define government efforts based on the data as follows:  
#   - "Overemphasizing economic development while neglecting environmental protection" = 1  
# - "Insufficient attention and inadequate investment in environmental protection" = 1  
# - "Unclear" = 2  
# - "Made efforts but with unsatisfactory results" = 2  
# - "Made significant efforts with some achievements" = 3  
# - "Achieved great success" = 3
CGSS2010_filter$gov_satisfaction <- ifelse(CGSS2010_filter$gov_satisfaction %in% c(1,2), 1, 
                                           ifelse(CGSS2010_filter$gov_satisfaction %in% c(3,8), 2, 
                                                  ifelse(CGSS2010_filter$gov_satisfaction %in% c(4,5), 3,CGSS2010_filter$gov_satisfaction)))

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
# - **Ordered categorical variables**: NEP, government satisfaction, subjective well-being
CGSS2010_filter$age <- as.numeric(CGSS2010_filter$age)
CGSS2010_filter$household_income <- as.numeric(CGSS2010_filter$household_income)
CGSS2010_filter$household_income <- scale(CGSS2010_filter$household_income) # Standardize household income to ensure consistency in variable scales.
CGSS2010_filter$household_income <- as.numeric(CGSS2010_filter$household_income)
CGSS2010_filter$education <- as.numeric(CGSS2010_filter$education)
CGSS2010_filter <- CGSS2010_filter %>% mutate_at(vars(sex), as.factor)
NEP_order <- c("1", "2", "3",'4','5')
CGSS2010_filter <- CGSS2010_filter %>% mutate_at(vars(NEP1,NEP2,NEP3,NEP4,NEP5,NEP6,NEP7,NEP8,NEP9,NEP10,NEP11,NEP12,NEP13,NEP14,NEP15), 
                                                 list(~ factor(., levels = NEP_order, ordered = TRUE)))
gov_order <- c("1", "2", "3")
CGSS2010_filter <- CGSS2010_filter %>% mutate(gov_satisfaction = factor(gov_satisfaction, levels = gov_order, ordered = TRUE))
happy_order <- c("1", "2", "3",'4','5')
CGSS2010_filter <- CGSS2010_filter %>% mutate(happiness = factor(happiness, levels = happy_order, ordered = TRUE))

# Sort the data according to the specified order.
custom_order <- c("NEP1", "NEP2", "NEP3",'NEP4','NEP5','NEP6','NEP7','NEP8','NEP9','NEP10','NEP11','NEP12','NEP13','NEP14',
                  'NEP15','sex','age','happiness','gov_satisfaction','education','household_income','birth') #'rental','marriage',
CGSS2010_filter <- CGSS2010_filter %>% select(all_of(custom_order))
CGSS2010_filter$Year <- 2010

# Write the data before imputation (i.e., CGSS2010 in Supplementary Table 7).
# write.csv(CGSS2010_filter, file = 'CGSS2010(插补前).csv', row.names = FALSE)

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
                               'H12_15','A36','A7a','H8','A62','type')]

# Rename column names
names(CGSS2021_filter)[names(CGSS2021_filter) == "A2"] <- "sex"
names(CGSS2021_filter)[names(CGSS2021_filter) == "A3_1"] <- "birth"
names(CGSS2021_filter)[names(CGSS2021_filter) == "A36"] <- "happiness"
names(CGSS2021_filter)[names(CGSS2021_filter) == "H8"] <- "gov_satisfaction"
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


# Define government efforts based on the data as follows:  
#   - "Overemphasizing economic development while neglecting environmental protection" = 1  
# - "Insufficient attention and inadequate investment in environmental protection" = 1  
# - "Unclear" = 2  
# - "Made efforts but with unsatisfactory results" = 2  
# - "Made significant efforts with some achievements" = 3  
# - "Achieved great success" = 3
CGSS2021_filter$gov_satisfaction <- ifelse(CGSS2021_filter$gov_satisfaction %in% c(1,2), 1, 
                                           ifelse(CGSS2021_filter$gov_satisfaction %in% c(3,98), 2, 
                                                  ifelse(CGSS2021_filter$gov_satisfaction %in% c(4,5), 3,CGSS2021_filter$gov_satisfaction)))

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
# - **Ordered categorical variables**: NEP, government satisfaction, subjective well-being
CGSS2021_filter$age <- as.numeric(CGSS2021_filter$age)
CGSS2021_filter$household_income <- as.numeric(CGSS2021_filter$household_income)
CGSS2021_filter$household_income <- scale(CGSS2021_filter$household_income) # Standardize household income to ensure consistency in variable scales.
CGSS2021_filter$household_income <- as.numeric(CGSS2021_filter$household_income)
CGSS2021_filter$education <- as.numeric(CGSS2021_filter$education)
CGSS2021_filter <- CGSS2021_filter %>% mutate_at(vars(sex), as.factor)
NEP_order <- c("1", "2", "3",'4','5')
CGSS2021_filter <- CGSS2021_filter %>% mutate_at(vars(NEP1,NEP2,NEP3,NEP4,NEP5,NEP6,NEP7,NEP8,NEP9,NEP10,NEP11,NEP12,NEP13,NEP14,NEP15), 
                                                 list(~ factor(., levels = NEP_order, ordered = TRUE)))
gov_order <- c("1", "2", "3")
CGSS2021_filter <- CGSS2021_filter %>% mutate(gov_satisfaction = factor(gov_satisfaction, levels = gov_order, ordered = TRUE))
happy_order <- c("1", "2", "3",'4','5')
CGSS2021_filter <- CGSS2021_filter %>% mutate(happiness = factor(happiness, levels = happy_order, ordered = TRUE))

# Sort the data according to the specified order.
custom_order <- c("NEP1", "NEP2", "NEP3",'NEP4','NEP5','NEP6','NEP7','NEP8','NEP9','NEP10','NEP11','NEP12','NEP13','NEP14',
                  'NEP15','sex','age','happiness','gov_satisfaction','education','household_income','birth')#'rental','marriage',
CGSS2021_filter <- CGSS2021_filter %>% select(all_of(custom_order))
CGSS2021_filter$Year <- 2021

# Write the data before imputation (i.e., CGSS2021 in Supplementary Table 7).
# write.csv(CGSS2021_filter, file = 'CGSS2021(插补前).csv', row.names = FALSE)

# Data imputation
impmice=mice(CGSS2021_filter,m=5,seed=2003)

# Write the 5 five sets of imputed datasets 2021.
n_iterations <- 5
for (i in 1:n_iterations) {
  assign(paste0("imp", i), complete(impmice, i))
  write.csv(get(paste0("imp", i)), file = paste0("2021_imp", i, ".csv"), row.names = FALSE)
}


### 1.1.4 Figure 1 ====
#### 1.1.4.1 Violin Plots ====
# setwd('C:\\') Your own path
CGSS2003 <- read.csv('CGSS2003(插补前).csv') # The data removed before imputation in section 1.2.1
CGSS2010 <- read.csv('CGSS2010(插补前).csv') # The data removed before imputation in section 1.2.2
CGSS2021 <- read.csv('CGSS2021(插补前).csv') # The data removed before imputation in section 1.2.3
CGSS <- rbind(CGSS2003,CGSS2010,CGSS2021) 
CGSS <- CGSS %>% mutate(cohort_cat = case_when(
  birth >= 1921 & birth <= 1930 ~ 1,
  birth >= 1931 & birth <= 1940 ~ 2,
  birth >= 1941 & birth <= 1950 ~ 3,
  birth >= 1951 & birth <= 1960 ~ 4,
  birth >= 1961 & birth <= 1970 ~ 5,
  birth >= 1971 & birth <= 1980 ~ 6,
  birth >= 1981 & birth <= 1990 ~ 7,
  birth >= 1991 & birth <= 2000 ~ 8,
  birth >= 2001 & birth <= 2010 ~ 9,
  TRUE ~ NA_integer_  
))
CGSS <- CGSS[!is.na(CGSS$cohort_cat), ] 
CGSS <- CGSS %>%
  mutate(cohort_cat = recode(cohort_cat,
                             `1` = "Cohort 1921-1930",
                             `2` = "Cohort 1931-1940",
                             `3` = "Cohort 1941-1950",
                             `4` = "Cohort 1951-1960",
                             `5` = "Cohort 1961-1970",
                             `6` = "Cohort 1971-1980",
                             `7` = "Cohort 1981-1990",
                             `8` = "Cohort 1991-2000",
                             `9` = "Cohort 2001-2010"))
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

##### 1.1.4.1.1 Violin plot for different years ====
p_period <- CGSS %>%
  ggplot(aes(x = Year, y = CNEP, fill = Year)) +
  geom_violin(width = 0.5, color = 'white', trim = TRUE) +  
  scale_fill_manual(values = c("#7179AD", "#2F66AC", "#3FB2C4")) +
  geom_boxplot(width = 0.1, color = "black", fill = 'white', alpha = 0.2) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),  
    axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 8)), 
    axis.title.y = element_text(size = 14, color = "black"), 
    axis.text.x = element_text(size = 10, color = "black", angle = 0, hjust = 0.5),  
    axis.text.y = element_text(size = 10, color = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "black"), 
    axis.ticks = element_line(color = "black") 
  ) +
  coord_cartesian(ylim = c(15, 51)) +
  ggtitle(" ") + 
  xlab("Year") + 
  ylab(" ") 

##### 1.1.4.1.2 Violin plot for different cohort ====
CGSS$cohort_cat <- gsub("^Cohort ", "", CGSS$cohort_cat)
p_cohort <- CGSS %>%
  ggplot( aes(x=cohort_cat, y=CNEP, fill=cohort_cat)) +
  geom_violin(width=0.5,color='white', adjust = 1) +
  scale_fill_manual(values = c("#A0B8C1","#88A7B6", "#8C9A94","#9B8C8A", "#9C7A7A", "#7F6C8A","#7C7F8C","#6C7B7F","#6B6D75")) +
  geom_boxplot(width=0.1, color="black",fill='white', alpha=0.2) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),  
    axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 8)),
    axis.title.y = element_text(size = 14, color = "black"),  
    axis.text.x = element_text(size = 8, color = "black", angle = 0, hjust = 0.5), 
    axis.text.y = element_text(size = 10, color = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "black"), 
    axis.ticks = element_line(color = "black")  
  ) +
  coord_cartesian(ylim = c(15, 51)) +
  ggtitle(" ") + 
  xlab("Cohort") + 
  ylab("CNEP Score")  

#### 1.1.4.2 Age Dot plot ====
p_age <- ggplot(CGSS, aes(age, CNEP)) +
  geom_point(size = 3.5, shape = 16 ,alpha = 0.1,color = "#3892AA", alpha = 0.6) +
  geom_smooth(method = "lm",color = "black",se = T,fill = "grey") +
  facet_grid( scales = "free_y") +
  theme_classic() +
  coord_cartesian(ylim = c(15, 52),xlim = c(20,80)) +
  xlab("Age") +
  ylab("CNEP Score") +
  theme(
    axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 8)), 
    axis.title.y = element_text(size = 14, color = "black"), 
    axis.text.x = element_text(size = 10, color = "black"), 
    axis.text.y = element_text(size = 10, color = "black"),  
    plot.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5)  
  )

#### 1.1.4.3 Figure 1 ====
library(patchwork)
Figure1 <- (p_age | p_period) / p_cohort
Figure1 <- Figure1 + plot_annotation(
  tag_levels = 'a') 
Figure1


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
      df_2003 <- read.csv(paste0("file path\\", file_2003))
      df_2010 <- read.csv(paste0("file path\\", file_2010))
      df_2021 <- read.csv(paste0("file path\\", file_2021))
      
      # Combined the data.
      combined_data <- rbind(df_2003, df_2010, df_2021)
      
      # Construct the output file path.
      output_file <- paste0("file path\\", (i - 1) * 25 + (j - 1) * 5 + k, ".csv")
      
      # Write to a new file.
      write.csv(combined_data, file = output_file, row.names = FALSE)
    }
  }
}



# 2. Separate Data Analyses ====
# Separate data analyses for the five imputed datasets of each year
# in 1.1.1  1.1.2  and  1.1.3

## 2.1 CGSS2003 ====
### 2.1.1 Read the data ====
# setwd('D:\\') Your own path
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

# Calculate the age used in the analysis (AgeR:divide by ten and subtract 2), as well as its square, and calculate CNEP.
data_list_modified <- lapply(data_list_modified, function(data) {
  data <- data %>%
    mutate(CNEP = NEP1 + NEP3 + NEP5 + NEP7 + NEP8 + NEP9 + NEP10 + NEP11 + NEP13 + NEP15)
  data$AgeR = data$age/10 - 2
  data$AgeR2 = data$AgeR^2
  data <- data %>% mutate_at(vars(sex,Year,gov_satisfaction,happiness), as.factor)
  return(data)
})

### 2.1.2 Data Analysis ====
result_2003 <- data.frame()
for (data in data_list_modified) {
  result_tem <- data.frame()
  fullmodel <- glm(CNEP ~ AgeR + AgeR2 + sex + happiness + gov_satisfaction + education + household_income ,
                   data = data)
  summary_fullmodel <- summary(fullmodel)
  # age_effect & age_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  result_tem[1,'age_effect'] <- fixed_effect[2,1]
  result_tem[1,'age_p'] <- fixed_effect[2,4]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect[3,1]
  result_tem[1,'age2_p'] <- fixed_effect[3,4]
  
  # Control Variables
  result_tem[1,'female_effect'] <- fixed_effect[4,1]
  result_tem[1,'female_p'] <- fixed_effect[4,4]
  result_tem[1,'happiness.2_effect'] <- fixed_effect[5,1]
  result_tem[1,'happiness.2_p'] <- fixed_effect[5,4]
  result_tem[1,'happiness.3_effect'] <- fixed_effect[6,1]
  result_tem[1,'happiness.3_p'] <- fixed_effect[6,4]
  result_tem[1,'happiness.4_effect'] <- fixed_effect[7,1]
  result_tem[1,'happiness.4_p'] <- fixed_effect[7,4]
  result_tem[1,'happiness.5_effect'] <- fixed_effect[8,1]
  result_tem[1,'happiness.5_p'] <- fixed_effect[8,4]
  result_tem[1,'satisfaction.2_effect'] <- fixed_effect[9,1]
  result_tem[1,'satisfaction.2_p'] <- fixed_effect[9,4]
  result_tem[1,'satisfaction.3_effect'] <- fixed_effect[10,1]
  result_tem[1,'satisfaction.3_p'] <- fixed_effect[10,4]
  result_tem[1,'education_effect'] <- fixed_effect[11,1]
  result_tem[1,'education_p'] <- fixed_effect[11,4]
  result_tem[1,'household_income_effect'] <- fixed_effect[12,1]
  result_tem[1,'household_income_p'] <- fixed_effect[12,4]
  
  result_tem[1,'r2m'] <- r.squaredGLMM(fullmodel)[,2]
  

  result_2003  <- rbind(result_2003 , result_tem)
}
# Output the results, i.e., the 2003 data in Supplementary Table 2.
# write.csv(result_2003[,1:23], file = 'result2003.csv', row.names = FALSE)

## 2.2 CGSS2010 ====
### 2.2.1 Read the data ====
# setwd('D:\\') Your own path
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

# Calculate the age used in the analysis (AgeR:divide by ten and subtract 2), as well as its square, and calculate CNEP.
data_list_modified <- lapply(data_list_modified, function(data) {
  data <- data %>%
    mutate(CNEP = NEP1 + NEP3 + NEP5 + NEP7 + NEP8 + NEP9 + NEP10 + NEP11 + NEP13 + NEP15)
  data$AgeR = data$age/10 - 2
  data$AgeR2 = data$AgeR^2
  data <- data %>% mutate_at(vars(sex,Year,gov_satisfaction,happiness), as.factor)
  return(data)
})

### 2.2.2 Data Analysis ====
result_2010 <- data.frame()
for (data in data_list_modified) {
  result_tem <- data.frame() 
  fullmodel <- glm(CNEP ~ AgeR + AgeR2 + sex + happiness + gov_satisfaction + education + household_income ,
                   data = data)
  summary_fullmodel <- summary(fullmodel)
  
  # age_effect & age_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  result_tem[1,'age_effect'] <- fixed_effect[2,1]
  result_tem[1,'age_p'] <- fixed_effect[2,4]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect[3,1]
  result_tem[1,'age2_p'] <- fixed_effect[3,4]
  
  # Control Variables
  result_tem[1,'female_effect'] <- fixed_effect[4,1]
  result_tem[1,'female_p'] <- fixed_effect[4,4]
  result_tem[1,'happiness.2_effect'] <- fixed_effect[5,1]
  result_tem[1,'happiness.2_p'] <- fixed_effect[5,4]
  result_tem[1,'happiness.3_effect'] <- fixed_effect[6,1]
  result_tem[1,'happiness.3_p'] <- fixed_effect[6,4]
  result_tem[1,'happiness.4_effect'] <- fixed_effect[7,1]
  result_tem[1,'happiness.4_p'] <- fixed_effect[7,4]
  result_tem[1,'happiness.5_effect'] <- fixed_effect[8,1]
  result_tem[1,'happiness.5_p'] <- fixed_effect[8,4]
  result_tem[1,'satisfaction.2_effect'] <- fixed_effect[9,1]
  result_tem[1,'satisfaction.2_p'] <- fixed_effect[9,4]
  result_tem[1,'satisfaction.3_effect'] <- fixed_effect[10,1]
  result_tem[1,'satisfaction.3_p'] <- fixed_effect[10,4]
  result_tem[1,'education_effect'] <- fixed_effect[11,1]
  result_tem[1,'education_p'] <- fixed_effect[11,4]
  result_tem[1,'household_income_effect'] <- fixed_effect[12,1]
  result_tem[1,'household_income_p'] <- fixed_effect[12,4]
  
  result_tem[1,'r2m'] <- r.squaredGLMM(fullmodel)[,2]
  
  result_2010  <- rbind(result_2010 , result_tem)
}
# Output the results, i.e., the 2010 data in Supplementary Table 2.
# write.csv(result_2010[,1:23], file = 'result2010.csv', row.names = FALSE)

## 2.3.CGSS2010 ====
### 2.3.1 Read the data ====
# setwd('D:\\') Your own path
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

# Calculate the age used in the analysis (AgeR:divide by ten and subtract 2), as well as its square, and calculate CNEP.
data_list_modified <- lapply(data_list_modified, function(data) {
  data <- data %>%
    mutate(CNEP = NEP1 + NEP3 + NEP5 + NEP7 + NEP8 + NEP9 + NEP10 + NEP11 + NEP13 + NEP15)
  data$AgeR = data$age/10 - 2
  data$AgeR2 = data$AgeR^2
  data <- data %>% mutate_at(vars(sex,Year,gov_satisfaction,happiness), as.factor)
  return(data)
})

### 2.3.2 Data Analysis ====
result_2021 <- data.frame()
for (data in data_list_modified) {
  result_tem <- data.frame() 
  fullmodel <- glm(CNEP ~ AgeR + AgeR2 + sex + happiness + gov_satisfaction + education + household_income ,
                   data = data)
  summary_fullmodel <- summary(fullmodel)
  
  # age_effect & age_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  result_tem[1,'age_effect'] <- fixed_effect[2,1]
  result_tem[1,'age_p'] <- fixed_effect[2,4]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect[3,1]
  result_tem[1,'age2_p'] <- fixed_effect[3,4]
  
  # Control Variables
  result_tem[1,'female_effect'] <- fixed_effect[4,1]
  result_tem[1,'female_p'] <- fixed_effect[4,4]
  result_tem[1,'happiness.2_effect'] <- fixed_effect[5,1]
  result_tem[1,'happiness.2_p'] <- fixed_effect[5,4]
  result_tem[1,'happiness.3_effect'] <- fixed_effect[6,1]
  result_tem[1,'happiness.3_p'] <- fixed_effect[6,4]
  result_tem[1,'happiness.4_effect'] <- fixed_effect[7,1]
  result_tem[1,'happiness.4_p'] <- fixed_effect[7,4]
  result_tem[1,'happiness.5_effect'] <- fixed_effect[8,1]
  result_tem[1,'happiness.5_p'] <- fixed_effect[8,4]
  result_tem[1,'satisfaction.2_effect'] <- fixed_effect[9,1]
  result_tem[1,'satisfaction.2_p'] <- fixed_effect[9,4]
  result_tem[1,'satisfaction.3_effect'] <- fixed_effect[10,1]
  result_tem[1,'satisfaction.3_p'] <- fixed_effect[10,4]
  result_tem[1,'education_effect'] <- fixed_effect[11,1]
  result_tem[1,'education_p'] <- fixed_effect[11,4]
  result_tem[1,'household_income_effect'] <- fixed_effect[12,1]
  result_tem[1,'household_income_p'] <- fixed_effect[12,4]
  
  result_tem[1,'r2m'] <- r.squaredGLMM(fullmodel)[,2]
  
  result_2021  <- rbind(result_2021 , result_tem)
}
# Output the results, i.e., the 2021 data in Supplementary Table 2.
# write.csv(result_2021[,1:23], file = 'result2021.csv', row.names = FALSE)



# 3. Main Statistical Analysis ====
## 3.1 Main Analysis ====
# The data analysis section of the main text.
# 125 complete datasets

### 3.1.1 Read the Data ====
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
    birth >= 1921 & birth <= 1930 ~ 1,
    birth >= 1931 & birth <= 1940 ~ 2,
    birth >= 1941 & birth <= 1950 ~ 3,
    birth >= 1951 & birth <= 1960 ~ 4,
    birth >= 1961 & birth <= 1970 ~ 5,
    birth >= 1971 & birth <= 1980 ~ 6,
    birth >= 1981 & birth <= 1990 ~ 7,
    birth >= 1991 & birth <= 2000 ~ 8,
    birth >= 2001 & birth <= 2010 ~ 9,
    TRUE ~ NA_integer_  # Define the samples that do not belong to these generations as NA.
  ))
  data <- data[!is.na(data$cohort_cat), ] 
  data$cohort_cat <- as.factor(data$cohort_cat)
  data <- data %>% mutate_at(vars(sex,Year,gov_satisfaction,happiness), as.factor)
  return(data)
})

### 3.1.2 Perform Coarsened Exact Matching ====
matstrict <- cem("Year", datalist=data_list_modified, drop=c('NEP1','NEP2','NEP3','NEP4','NEP5','NEP6','NEP7','NEP8','NEP9',
                                                             'NEP10','NEP11','NEP12','NEP13','NEP14','NEP15',
                                                             'age','birth','Year','NEP','CNEP',
                                                             'AgeR','AgeR2','cohort_cat'))
matstrict
matched_data_list <- list()
for (i in seq_along(data_list)) {
  match_col <- paste0("match", i)
  matched_data_list[[i]] <- data_list_modified[[i]][matstrict[[match_col]]$matched, ]
}

### 3.1.3 Construct the HAPC model for each matched dataset ====
result_CNEP <- data.frame()
for (data in matched_data_list) {
  result_tem <- data.frame() 
  
  fullmodel <- lmer(CNEP ~ AgeR + AgeR2 + sex + happiness + gov_satisfaction + education + household_income +
                      + (1 | cohort_cat) + (1 | Year),data = data)
  nullmodel <- glm(CNEP ~ sex + happiness + gov_satisfaction + education + household_income ,
                   data = data)
  summary_fullmodel <- summary(fullmodel)
  
  # anova_p
  anova <- anova(fullmodel,nullmodel)
  result_tem[1,'anova_p'] <- anova$`Pr(>Chisq)`[2]
  
  # variance_cohort
  random_effects_var <- as.data.frame(VarCorr(fullmodel))
  result_tem[1,'variance_cohort'] <- random_effects_var[1,4]
  
  # cohort_p
  significant <- ranova(fullmodel)
  result_tem[1,'cohort_p'] <- significant$`Pr(>Chisq)`[2]
  
  # variance_period
  result_tem[1,'variance_period'] <- random_effects_var[2,4]
  
  # period_p
  result_tem[1,'period_p'] <- significant$`Pr(>Chisq)`[3]
  
  # Residual
  result_tem[1,'Residual'] <- random_effects_var[3,4]
  
  # Cohort_effect1-9
  Cohort_effect <- ranef(fullmodel)$cohort_cat
  result_tem[1,'Cohort_effect1921-1930'] <- Cohort_effect[1,1]
  result_tem[1,'Cohort_effect1931-1940'] <- Cohort_effect[2,1]
  result_tem[1,'Cohort_effect1941-1950'] <- Cohort_effect[3,1]
  result_tem[1,'Cohort_effect1951-1960'] <- Cohort_effect[4,1]
  result_tem[1,'Cohort_effect1961-1970'] <- Cohort_effect[5,1]
  result_tem[1,'Cohort_effect1971-1980'] <- Cohort_effect[6,1]
  result_tem[1,'Cohort_effect1981-1990'] <- Cohort_effect[7,1]
  result_tem[1,'Cohort_effect1991-2000'] <- Cohort_effect[8,1]
  result_tem[1,'Cohort_effect2001-2010'] <- Cohort_effect[9,1]
  
  # Period_effect2003-2021 
  Period_effect <- ranef(fullmodel)$Year
  result_tem[1,'Period_effect2003'] <- Period_effect[1,1]
  result_tem[1,'Period_effect2010'] <- Period_effect[2,1]
  result_tem[1,'Period_effect2021'] <- Period_effect[3,1]
  
  # age_effect & age_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  result_tem[1,'age_effect'] <- fixed_effect[2,1]
  result_tem[1,'age_p'] <- fixed_effect[2,5]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect[3,1]
  result_tem[1,'age2_p'] <- fixed_effect[3,5]
  
  # Control Variables
  result_tem[1,'female_effect'] <- fixed_effect[4,1]
  result_tem[1,'female_p'] <- fixed_effect[4,5]
  result_tem[1,'happiness.2_effect'] <- fixed_effect[5,1]
  result_tem[1,'happiness.2_p'] <- fixed_effect[5,5]
  result_tem[1,'happiness.3_effect'] <- fixed_effect[6,1]
  result_tem[1,'happiness.3_p'] <- fixed_effect[6,5]
  result_tem[1,'happiness.4_effect'] <- fixed_effect[7,1]
  result_tem[1,'happiness.4_p'] <- fixed_effect[7,5]
  result_tem[1,'happiness.5_effect'] <- fixed_effect[8,1]
  result_tem[1,'happiness.5_p'] <- fixed_effect[8,5]
  result_tem[1,'satisfaction.2_effect'] <- fixed_effect[9,1]
  result_tem[1,'satisfaction.2_p'] <- fixed_effect[9,5]
  result_tem[1,'satisfaction.3_effect'] <- fixed_effect[10,1]
  result_tem[1,'satisfaction.3_p'] <- fixed_effect[10,5]
  result_tem[1,'education_effect'] <- fixed_effect[11,1]
  result_tem[1,'education_p'] <- fixed_effect[11,5]
  result_tem[1,'household_income_effect'] <- fixed_effect[12,1]
  result_tem[1,'household_income_p'] <- fixed_effect[12,5]
  
  # R2
  r2 <- r.squaredGLMM(fullmodel)
  result_tem[1,'HAPC_r2m'] <- r2[,1]
  result_tem[1,'HAPC_r2c'] <- r2[,2]
  result_tem[1,'nullmodel_r2m'] <- r.squaredGLMM(nullmodel)[,2]
  

  result_CNEP <- rbind(result_CNEP, result_tem)
}
# write.csv(result_CNEP[,41:43], file = 'result_CNEP_R2.csv', row.names = FALSE)
# write.csv(result_CNEP[,1:40], file = 'result_CNEP.csv', row.names = FALSE)
# These two combined form Supplementary Table 1.
# For convenience in plotting, store them in two separate files here, and do not consider R2 when plotting.

### 3.1.4 Figure 2 ====
result_CNEP <- read.csv('result_CNEP.csv') # The result_CNEP.csv in section 3.1.3 (the version without R2)

# Test each p-value.
t.test(result_CNEP$cohort_p, mu = 0.05, alternative = "less") #p-value = 1 
t.test(result_CNEP$period_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$age_p, mu = 0.05, alternative = "less") #p-value < 2.2e-16 *
t.test(result_CNEP$age2_p, mu = 0.05, alternative = "less") #p-value < 2.2e-16 *
t.test(result_CNEP$female_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$happiness.2_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP$happiness.3_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP$happiness.4_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP$happiness.5_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP$satisfaction.2_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$satisfaction.3_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$education_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$household_income_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***

names(result_CNEP)[names(result_CNEP) == "age_effect"] <- "Age"
names(result_CNEP)[names(result_CNEP) == "age2_effect"] <- "Age-quadratic"
names(result_CNEP)[names(result_CNEP) == "female_effect"] <- "Female(Ref:male)"
names(result_CNEP)[names(result_CNEP) == "happiness.2_effect"] <- "unhappy(Ref:Extremely unhappy)"
names(result_CNEP)[names(result_CNEP) == "happiness.3_effect"] <- "Neutral(Ref:Extremely unhappy)"
names(result_CNEP)[names(result_CNEP) == "happiness.4_effect"] <- "happy(Ref:Extremely unhappy)"
names(result_CNEP)[names(result_CNEP) == "happiness.5_effect"] <- "Extremely happy(Ref:Extremely unhappy)"
names(result_CNEP)[names(result_CNEP) == "satisfaction.2_effect"] <- "Neutral(Ref:dissatisfied)"
names(result_CNEP)[names(result_CNEP) == "satisfaction.3_effect"] <- "Satisfied(Ref:dissatisfied)"
names(result_CNEP)[names(result_CNEP) == "education_effect"] <- "Education"
names(result_CNEP)[names(result_CNEP) == "household_income_effect"] <- "Household income"

names(result_CNEP)[names(result_CNEP) == "Cohort_effect1921.1930"] <- "Cohort 1921-1930"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1931.1940"] <- "Cohort 1931-1940"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1941.1950"] <- "Cohort 1941-1950"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1951.1960"] <- "Cohort 1951-1960"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1961.1970"] <- "Cohort 1961-1970"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1971.1980"] <- "Cohort 1971-1980"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1981.1990"] <- "Cohort 1981-1990"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1991.2000"] <- "Cohort 1991-2000"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect2001.2010"] <- "Cohort 2001-2010"
names(result_CNEP)[names(result_CNEP) == "Period_effect2003"] <- "Period 2003"
names(result_CNEP)[names(result_CNEP) == "Period_effect2010"] <- "Period 2010"
names(result_CNEP)[names(result_CNEP) == "Period_effect2021"] <- "Period 2021"

CNEP_order <- c("Age", "Age-quadratic","Period 2003","Period 2010","Period 2021",
                "Cohort 1921-1930","Cohort 1931-1940","Cohort 1941-1950","Cohort 1951-1960",
                "Cohort 1961-1970","Cohort 1971-1980","Cohort 1981-1990","Cohort 1991-2000",
                "Cohort 2001-2010")
result_CNEP <- result_CNEP %>% select(all_of(CNEP_order))
CNEP <- melt(result_CNEP)
colnames(CNEP) <- c('sample','value')
CNEP$sample <- as.character(CNEP$sample)
CNEP <- CNEP[!grepl("_p$", CNEP$sample), ]
CNEP <- CNEP[!grepl("^variance", CNEP$sample), ]
CNEP <- CNEP[!grepl("^Resid", CNEP$sample), ]
CNEP$group <- ifelse(CNEP$sample %in% c("Age", "Age-quadratic",'Female(Ref:male)','unhappy(Ref:Extremely unhappy)','Neutral(Ref:Extremely unhappy)',
                                        'happy(Ref:Extremely unhappy)','Extremely happy(Ref:Extremely unhappy)','Neutral(Ref:dissatisfied)','Satisfied(Ref:dissatisfied)',
                                        'Education','Household income'), "fix", "random")
CNEP <- CNEP[rev(seq_len(nrow(CNEP))), ]
unique_strings <- unique(unlist(CNEP$sample))
CNEP$sample <- factor(CNEP$sample,levels = unique_strings)
CNEP$group <- as.factor(CNEP$group)
CNEP$facet <- rep("APC effect",times=1750)

CNEP <- CNEP %>%
  mutate(facet = case_when(
    grepl("^Cohort", sample) ~ "Cohort effect",
    grepl("^Period", sample) ~ "Period effect ***",
    grepl("^Age", sample) ~ "Age effect",
    TRUE ~ facet 
  ))

CNEPage <- CNEP[1501:1750,]
CNEPcohort <- CNEP[1:1125,]
CNEPperiod <- CNEP[1126:1500,]

#### 3.1.4.1 Age Part ====
p <- ggplot(CNEPage,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
CNEPage %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apcage <- ggplot(CNEPage,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#459f81"))+
  scale_color_manual(values = c("#459f81"))+
  scale_y_continuous(limits = c(-0.6, 0.6))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Coefficient",x=NULL)+
  annotate("rect", xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  annotate('text', label = '*', x =2, y =0.5, angle=-90, size =8,color="black") + 
  annotate('text', label = '*', x =1, y =-0.12, angle=-90, size =8,color="black") + 
  facet_grid(~ facet)
p_apcage

#### 3.1.4.2 Cohort Part ====
CNEPcohort$sample <- gsub("^Cohort ", "", CNEPcohort$sample)
CNEPcohort <- CNEPcohort %>%
  mutate(sample = factor(sample, levels = c("2001-2010", "1991-2000",'1981-1990','1971-1980','1961-1970',
                                            '1951-1960','1941-1950','1931-1940','1921-1930')))

p <- ggplot(CNEPcohort,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPcohort %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apccohort <- ggplot(CNEPcohort,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#3769AE"))+
  scale_color_manual(values = c("#3769AE"))+
  scale_y_continuous(limits = c(-0.3, 0.3))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Random Intercept",x=NULL)+
  annotate("rect", xmin = 0, xmax = 9.6, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  facet_grid(~ facet)
p_apccohort


#### 3.1.4.3 Period Part ====
CNEPperiod$sample <- gsub("^Period ", "", CNEPperiod$sample)
CNEPperiod <- CNEPperiod %>%
  mutate(sample = factor(sample, levels = c("2021", "2010",'2003')))

p <- ggplot(CNEPperiod,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPperiod %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apcperiod <- ggplot(CNEPperiod,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#3769AE"))+
  scale_color_manual(values = c("#3769AE"))+
  scale_y_continuous(limits = c(-2, 2))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) + 
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Random Intercept",x=NULL)+
  annotate("rect", xmin = 0, xmax = 3.6, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  facet_grid(~ facet)
p_apcperiod

#### 3.1.4.4 Figure 2 ====
Figure2 <- p_apcage / p_apcperiod / p_apccohort + 
  plot_layout(heights = c(2.2, 3, 9))
Figure2

### 3.1.5 Extended Data Fig1 ====
plot_model <- function(data, index) {
  model <- lmer(CNEP ~ AgeR + I(AgeR^2) + sex + happiness + gov_satisfaction + 
                  education + household_income + (1 | cohort_cat) + (1 | Year), 
                data = data)
  pred <- ggeffect(model, terms = "AgeR [all]")
  p <- ggplot(pred, aes(x = x, y = predicted)) +
    geom_line(size = 1, color = "blue") + 
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "blue") +
    labs(x = "AgeR", y = "Predicted CNEP", title = paste("Dataset", index)) +
    theme_minimal()
  return(p)
}

plots <- map2(matched_data_list, 1:length(data_list), plot_model)

ExtendedDataFig1 <- wrap_plots(plots) + plot_layout(ncol = 10)  
ggsave("ExtendedDataFig1.png", ExtendedDataFig1, width = 20, height = 25, units = "in", dpi = 300)
print(ExtendedDataFig1)

### 3.1.6 Supplementary Fig. 1 ====
result_CNEP <- read.csv('result_CNEP.csv') # The result_CNEP.csv in section 3.1.3 (the version without R2)
t.test(result_CNEP$cohort_p, mu = 0.05, alternative = "less") #p-value = 1 
t.test(result_CNEP$period_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$age_p, mu = 0.05, alternative = "less") #p-value < 2.2e-16 *
t.test(result_CNEP$age2_p, mu = 0.05, alternative = "less") #p-value < 2.2e-16 *
t.test(result_CNEP$female_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$happiness.2_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP$happiness.3_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP$happiness.4_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP$happiness.5_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP$satisfaction.2_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$satisfaction.3_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$education_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP$household_income_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***

names(result_CNEP)[names(result_CNEP) == "age_effect"] <- "Age"
names(result_CNEP)[names(result_CNEP) == "age2_effect"] <- "Age-quadratic"
names(result_CNEP)[names(result_CNEP) == "female_effect"] <- "Female(Ref:male)"
names(result_CNEP)[names(result_CNEP) == "happiness.2_effect"] <- "unhappy(Ref:Extremely unhappy)"
names(result_CNEP)[names(result_CNEP) == "happiness.3_effect"] <- "Neutral(Ref:Extremely unhappy)"
names(result_CNEP)[names(result_CNEP) == "happiness.4_effect"] <- "happy(Ref:Extremely unhappy)"
names(result_CNEP)[names(result_CNEP) == "happiness.5_effect"] <- "Extremely happy(Ref:Extremely unhappy)"
names(result_CNEP)[names(result_CNEP) == "satisfaction.2_effect"] <- "Neutral(Ref:dissatisfied)"
names(result_CNEP)[names(result_CNEP) == "satisfaction.3_effect"] <- "Satisfied(Ref:dissatisfied)"
names(result_CNEP)[names(result_CNEP) == "education_effect"] <- "Education"
names(result_CNEP)[names(result_CNEP) == "household_income_effect"] <- "Household income"

names(result_CNEP)[names(result_CNEP) == "Cohort_effect1921.1930"] <- "Cohort 1921-1930"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1931.1940"] <- "Cohort 1931-1940"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1941.1950"] <- "Cohort 1941-1950"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1951.1960"] <- "Cohort 1951-1960"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1961.1970"] <- "Cohort 1961-1970"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1971.1980"] <- "Cohort 1971-1980"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1981.1990"] <- "Cohort 1981-1990"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect1991.2000"] <- "Cohort 1991-2000"
names(result_CNEP)[names(result_CNEP) == "Cohort_effect2001.2010"] <- "Cohort 2001-2010"
names(result_CNEP)[names(result_CNEP) == "Period_effect2003"] <- "Period 2003"
names(result_CNEP)[names(result_CNEP) == "Period_effect2010"] <- "Period 2010"
names(result_CNEP)[names(result_CNEP) == "Period_effect2021"] <- "Period 2021"

CNEP_order <- c("Age", "Age-quadratic","Female(Ref:male)","unhappy(Ref:Extremely unhappy)","Neutral(Ref:Extremely unhappy)",
                "happy(Ref:Extremely unhappy)","Extremely happy(Ref:Extremely unhappy)","Neutral(Ref:dissatisfied)",
                "Satisfied(Ref:dissatisfied)","Education","Household income",
                "Period 2003","Period 2010","Period 2021",
                "Cohort 1921-1930","Cohort 1931-1940","Cohort 1941-1950","Cohort 1951-1960",
                "Cohort 1961-1970","Cohort 1971-1980","Cohort 1981-1990","Cohort 1991-2000",
                "Cohort 2001-2010")
result_CNEP <- result_CNEP %>% select(all_of(CNEP_order))
CNEP <- melt(result_CNEP)
colnames(CNEP) <- c('sample','value')
CNEP$sample <- as.character(CNEP$sample)
CNEP <- CNEP[!grepl("_p$", CNEP$sample), ]
CNEP <- CNEP[!grepl("^variance", CNEP$sample), ]
CNEP <- CNEP[!grepl("^Resid", CNEP$sample), ]
CNEP$group <- ifelse(CNEP$sample %in% c("Age", "Age-quadratic",'Female(Ref:male)','unhappy(Ref:Extremely unhappy)','Neutral(Ref:Extremely unhappy)',
                                        'happy(Ref:Extremely unhappy)','Extremely happy(Ref:Extremely unhappy)','Neutral(Ref:dissatisfied)','Satisfied(Ref:dissatisfied)',
                                        'Education','Household income'), "fix", "random")
CNEP <- CNEP[rev(seq_len(nrow(CNEP))), ]
unique_strings <- unique(unlist(CNEP$sample))
CNEP$sample <- factor(CNEP$sample,levels = unique_strings)
CNEP$group <- as.factor(CNEP$group)
CNEP$facet <- rep("APC effect",times=2875)

CNEP <- CNEP %>%
  mutate(facet = case_when(
    grepl("^Cohort", sample) ~ "Cohort effect",
    grepl("^Period", sample) ~ "Period effect ***",
    TRUE ~  "Fixed effect"  
  ))

CNEPage <- CNEP[1501:2875,]
CNEPcohort <- CNEP[1:1125,]
CNEPperiod <- CNEP[1126:1500,]

#### 3.1.6.1 Age Part ====
p <- ggplot(CNEPage,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
CNEPage %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apcage <- ggplot(CNEPage,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#459f81"))+
  scale_color_manual(values = c("#459f81"))+
  scale_y_continuous(limits = c(-2, 2))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Coefficient",x=NULL)+
  annotate("rect", xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  annotate('text', label = '*', x =11, y =0.52, angle=-90, size =8,color="black") + 
  annotate('text', label = '*', x =10, y =-0.2, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =9, y =-0.55, angle=-90, size =8,color="black") + 
  annotate('text', label = '*', x =8, y =-3, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =4, y =-1.3, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =3, y =-0.87, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =2, y =0.44, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =1, y =1.25, angle=-90, size =8,color="black") + 
  facet_grid(~ facet)
p_apcage

#### 3.1.6.2 Cohort Part ====
CNEPcohort$sample <- gsub("^Cohort ", "", CNEPcohort$sample)
CNEPcohort <- CNEPcohort %>%
  mutate(sample = factor(sample, levels = c("2001-2010", "1991-2000",'1981-1990','1971-1980','1961-1970',
                                            '1951-1960','1941-1950','1931-1940','1921-1930')))

p <- ggplot(CNEPcohort,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPcohort %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apccohort <- ggplot(CNEPcohort,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#3769AE"))+
  scale_color_manual(values = c("#3769AE"))+
  scale_y_continuous(limits = c(-0.4, 0.4))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Random Intercept",x=NULL)+
  annotate("rect", xmin = 0, xmax = 9.6, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  facet_grid(~ facet)
p_apccohort

#### 3.1.6.3 Period Part ====
CNEPperiod$sample <- gsub("^Period ", "", CNEPperiod$sample)
CNEPperiod <- CNEPperiod %>%
  mutate(sample = factor(sample, levels = c("2021", "2010",'2003')))

p <- ggplot(CNEPperiod,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPperiod %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apcperiod <- ggplot(CNEPperiod,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#3769AE"))+
  scale_color_manual(values = c("#3769AE"))+
  scale_y_continuous(limits = c(-2, 2))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Random Intercept",x=NULL)+
  annotate("rect", xmin = 0, xmax = 3.6, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  facet_grid(~ facet)
p_apcperiod

#### 3.1.6.4 Supplementary Fig 1 ====
SupplementaryFig1 <- p_apcage / p_apcperiod / p_apccohort + 
  plot_layout(heights = c(11, 3, 9))
SupplementaryFig1

## 3.2 Sensitivity analysis ====
### 3.2.1 The period variable was treated as a fixed effect ====
#### 3.2.1.1 Read the Data ====
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
    birth >= 1921 & birth <= 1930 ~ 1,
    birth >= 1931 & birth <= 1940 ~ 2,
    birth >= 1941 & birth <= 1950 ~ 3,
    birth >= 1951 & birth <= 1960 ~ 4,
    birth >= 1961 & birth <= 1970 ~ 5,
    birth >= 1971 & birth <= 1980 ~ 6,
    birth >= 1981 & birth <= 1990 ~ 7,
    birth >= 1991 & birth <= 2000 ~ 8,
    birth >= 2001 & birth <= 2010 ~ 9,
    TRUE ~ NA_integer_  
  ))
  data <- data[!is.na(data$cohort_cat), ] 
  data$cohort_cat <- as.factor(data$cohort_cat)
  data <- data %>% mutate_at(vars(sex,Year,gov_satisfaction,happiness), as.factor)
  return(data)
})

#### 3.2.1.2 Perform Coarsened Exact Matching ====
matstrict <- cem("Year", datalist=data_list_modified, drop=c('NEP1','NEP2','NEP3','NEP4','NEP5','NEP6','NEP7','NEP8','NEP9',
                                                             'NEP10','NEP11','NEP12','NEP13','NEP14','NEP15',
                                                             'age','birth','Year','NEP','CNEP',
                                                             'AgeR','AgeR2','cohort_cat'))
matstrict
matched_data_list <- list()
for (i in seq_along(data_list)) {
  match_col <- paste0("match", i)
  matched_data_list[[i]] <- data_list_modified[[i]][matstrict[[match_col]]$matched, ]
}

#### 3.2.1.3 Construct the HAPC model for each matched dataset ====
result_CNEP_sensitive_fix <- data.frame()
for (data in matched_data_list) {型
  result_tem <- data.frame() 
  
  fullmodel <- lmer(CNEP ~ Year + AgeR + AgeR2 + sex + happiness + gov_satisfaction + education + household_income +
                      + (1 | cohort_cat) ,data = data) # The period variable was treated as a fixed effect
  nullmodel <- glm(CNEP ~ happiness + gov_satisfaction + education + household_income ,
                   data = data)
  summary_fullmodel <- summary(fullmodel)
  
  # anova_p
  anova <- anova(fullmodel,nullmodel)
  result_tem[1,'anova_p'] <- anova$`Pr(>Chisq)`[2]
  
  # variance_cohort
  random_effects_var <- as.data.frame(VarCorr(fullmodel))
  result_tem[1,'variance_cohort'] <- random_effects_var[1,4]
  
  # cohort_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  
  significant <- ranova(fullmodel)
  result_tem[1,'cohort_p'] <- significant$`Pr(>Chisq)`[2]
  
  # period_p
  result_tem[1,'period_p2010'] <- fixed_effect[2,5]
  result_tem[1,'period_p2021'] <- fixed_effect[3,5]
  
  # Residual
  result_tem[1,'Residual'] <- random_effects_var[2,4] 
  
  # Cohort_effect1-9
  Cohort_effect <- ranef(fullmodel)$cohort_cat
  result_tem[1,'Cohort_effect1921-1930'] <- Cohort_effect[1,1]
  result_tem[1,'Cohort_effect1931-1940'] <- Cohort_effect[2,1]
  result_tem[1,'Cohort_effect1941-1950'] <- Cohort_effect[3,1]
  result_tem[1,'Cohort_effect1951-1960'] <- Cohort_effect[4,1]
  result_tem[1,'Cohort_effect1961-1970'] <- Cohort_effect[5,1]
  result_tem[1,'Cohort_effect1971-1980'] <- Cohort_effect[6,1]
  result_tem[1,'Cohort_effect1981-1990'] <- Cohort_effect[7,1]
  result_tem[1,'Cohort_effect1991-2000'] <- Cohort_effect[8,1]
  result_tem[1,'Cohort_effect2001-2010'] <- Cohort_effect[9,1]
  
  # Period_effect2003 
  Period_effect <- ranef(fullmodel)$Year
  result_tem[1,'Period_effect2010'] <- fixed_effect[2,1]
  result_tem[1,'Period_effect2021'] <- fixed_effect[3,1]
  
  # age_effect & age_p
  result_tem[1,'age_effect'] <- fixed_effect[4,1]
  result_tem[1,'age_p'] <- fixed_effect[4,5]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect[5,1]
  result_tem[1,'age2_p'] <- fixed_effect[5,5]
  
  # 控制变量
  result_tem[1,'female_effect'] <- fixed_effect[6,1]
  result_tem[1,'female_p'] <- fixed_effect[6,5]
  result_tem[1,'happiness.2_effect'] <- fixed_effect[7,1]
  result_tem[1,'happiness.2_p'] <- fixed_effect[7,5]
  result_tem[1,'happiness.3_effect'] <- fixed_effect[8,1]
  result_tem[1,'happiness.3_p'] <- fixed_effect[8,5]
  result_tem[1,'happiness.4_effect'] <- fixed_effect[9,1]
  result_tem[1,'happiness.4_p'] <- fixed_effect[9,5]
  result_tem[1,'happiness.5_effect'] <- fixed_effect[10,1]
  result_tem[1,'happiness.5_p'] <- fixed_effect[10,5]
  result_tem[1,'satisfaction.2_effect'] <- fixed_effect[11,1]
  result_tem[1,'satisfaction.2_p'] <- fixed_effect[11,5]
  result_tem[1,'satisfaction.3_effect'] <- fixed_effect[12,1]
  result_tem[1,'satisfaction.3_p'] <- fixed_effect[12,5]
  result_tem[1,'education_effect'] <- fixed_effect[13,1]
  result_tem[1,'education_p'] <- fixed_effect[13,5]
  result_tem[1,'household_income_effect'] <- fixed_effect[14,1]
  result_tem[1,'household_income_p'] <- fixed_effect[14,5]
  
  #R2
  r2 <- r.squaredGLMM(fullmodel)
  result_tem[1,'r2m'] <- r2[,1]
  result_tem[1,'r2c'] <- r2[,2]
  
  result_CNEP_sensitive_fix <- rbind(result_CNEP_sensitive_fix, result_tem)
}
# write.csv(result_CNEP_sensitive_fix[,40:41], file = 'result_CNEP_sensitive_fix_R2.csv', row.names = FALSE)
# write.csv(result_CNEP_sensitive_fix[,1:39], file = 'result_CNEP_sensitive_fix.csv', row.names = FALSE)
# These two combined form Supplementary Table 8.
# For convenience in plotting, store them in two separate files here, and do not consider R2 when plotting.

#### 3.2.1.4 Supplementary Fig 2 ====
result_CNEP_sensitive_fix <- read.csv('result_CNEP_sensitive_fix.csv') 
# The result_CNEP_sensitive_fix.csv in section 3.2.1.3 (the version without R2)

t.test(result_CNEP_sensitive_fix$cohort_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive_fix$period_p2010, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive_fix$period_p2021, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive_fix$age_p, mu = 0.05, alternative = "less") #p-value < 2.2e-16 *
t.test(result_CNEP_sensitive_fix$age2_p, mu = 0.0, alternative = "less") #p-value < 2.2e-16 *
t.test(result_CNEP_sensitive_fix$female_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive_fix$happiness.2_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive_fix$happiness.3_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive_fix$happiness.4_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive_fix$happiness.5_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive_fix$satisfaction.2_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive_fix$satisfaction.3_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive_fix$education_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive_fix$household_income_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***

names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Period_effect2010"] <- "Period 2010"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Period_effect2021"] <- "Period 2021"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "age_effect"] <- "Age"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "age2_effect"] <- "Age-quadratic"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "female_effect"] <- "Female(Ref:male)"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "happiness.2_effect"] <- "unhappy(Ref:Extremely unhappy)"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "happiness.3_effect"] <- "Neutral(Ref:Extremely unhappy)"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "happiness.4_effect"] <- "happy(Ref:Extremely unhappy)"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "happiness.5_effect"] <- "Extremely happy(Ref:Extremely unhappy)"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "satisfaction.2_effect"] <- "Neutral(Ref:dissatisfied)"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "satisfaction.3_effect"] <- "Satisfied(Ref:dissatisfied)"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "education_effect"] <- "Education"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "household_income_effect"] <- "Household income"

names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect1921.1930"] <- "Cohort 1921-1930"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect1931.1940"] <- "Cohort 1931-1940"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect1941.1950"] <- "Cohort 1941-1950"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect1951.1960"] <- "Cohort 1951-1960"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect1961.1970"] <- "Cohort 1961-1970"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect1971.1980"] <- "Cohort 1971-1980"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect1981.1990"] <- "Cohort 1981-1990"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect1991.2000"] <- "Cohort 1991-2000"
names(result_CNEP_sensitive_fix)[names(result_CNEP_sensitive_fix) == "Cohort_effect2001.2010"] <- "Cohort 2001-2010"

CNEP_order <- c("Age", "Age-quadratic",'Female(Ref:male)','unhappy(Ref:Extremely unhappy)','Neutral(Ref:Extremely unhappy)',
                'happy(Ref:Extremely unhappy)','Extremely happy(Ref:Extremely unhappy)','Neutral(Ref:dissatisfied)','Satisfied(Ref:dissatisfied)',
                'Education','Household income',"Period 2010","Period 2021",
                "Cohort 1921-1930","Cohort 1931-1940","Cohort 1941-1950","Cohort 1951-1960",
                "Cohort 1961-1970","Cohort 1971-1980","Cohort 1981-1990","Cohort 1991-2000",
                "Cohort 2001-2010")
result_CNEP_sensitive_fix <- result_CNEP_sensitive_fix %>% select(all_of(CNEP_order))
CNEP <- melt(result_CNEP_sensitive_fix)
colnames(CNEP) <- c('sample','value')
CNEP$sample <- as.character(CNEP$sample)
CNEP <- CNEP[!grepl("_p$", CNEP$sample), ]
CNEP <- CNEP[!grepl("^variance", CNEP$sample), ]
CNEP <- CNEP[!grepl("^Resid", CNEP$sample), ]
CNEP$group <- ifelse(CNEP$sample %in% c("Age", "Age-quadratic",'Female(Ref:male)','unhappy(Ref:Extremely unhappy)','Neutral(Ref:Extremely unhappy)',
                                        'happy(Ref:Extremely unhappy)','Extremely happy(Ref:Extremely unhappy)','Neutral(Ref:dissatisfied)','Satisfied(Ref:dissatisfied)',
                                        'Education','Household income',"Period 2010","Period 2021"), "fix", "random")
CNEP <- CNEP[rev(seq_len(nrow(CNEP))), ]
unique_strings <- unique(unlist(CNEP$sample))
CNEP$sample <- factor(CNEP$sample,levels = unique_strings)
CNEP$group <- as.factor(CNEP$group)
CNEP$facet <- rep("APC effect",times=2750)

CNEP <- CNEP %>%
  mutate(facet = case_when(
    grepl("^Cohort", sample) ~ "Cohort effect",
    TRUE ~  "Fixed effect"  
  ))

CNEPage <- CNEP[1126:2750,]
CNEPcohort <- CNEP[1:1125,]

##### 3.2.1.4.1 Age & Period Part ====
p <- ggplot(CNEPage,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPage %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apcage <- ggplot(CNEPage,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#459f81"))+
  scale_color_manual(values = c("#459f81"))+
  scale_y_continuous(limits = c(-2.5, 2.5))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Coefficient",x=NULL)+
  annotate("rect", xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  annotate('text', label = '*', x =13, y =0.7, angle=-90, size =8,color="black") + 
  annotate('text', label = '*', x =12, y =-0.2, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =11, y =-0.45, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =6, y =-1.35, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =5, y =-0.85, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =4, y =0.5, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =3, y =1.45, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =1, y =-2.2, angle=-90, size =8,color="black") + 
  facet_grid(~ facet)
p_apcage

##### 3.2.1.4.2 Cohort Part ====
CNEPcohort$sample <- gsub("^Cohort ", "", CNEPcohort$sample)
CNEPcohort <- CNEPcohort %>%
  mutate(sample = factor(sample, levels = c("2001-2010", "1991-2000",'1981-1990','1971-1980','1961-1970',
                                            '1951-1960','1941-1950','1931-1940','1921-1930')))

p <- ggplot(CNEPcohort,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPcohort %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean
p_apccohort <- ggplot(CNEPcohort,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#3769AE"))+
  scale_color_manual(values = c("#3769AE"))+
  scale_y_continuous(limits = c(-0.4, 0.4))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Random Intercept",x=NULL)+
  annotate("rect", xmin = 0, xmax = 9.6, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  facet_grid(~ facet)
p_apccohort

##### 3.2.1.4.3 Supplementary Fig 2 ====
SupplementaryFig2 <- p_apcage /  p_apccohort + 
  plot_layout(heights = c(14, 9))
SupplementaryFig2


### 3.2.2 two cohorts with small sample sizes (i.e., 1921-1930 and 2001-2010) were excluded ====
#### 3.2.2.1 Read the Data ====
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

# Calculate the age used in the analysis (AgeR:divide by ten and subtract 2), as well as its square, and calculate CNEP.
data_list_modified <- lapply(data_list_modified, function(data) {
  data <- data %>%
    mutate(CNEP = NEP1 + NEP3 + NEP5 + NEP7 + NEP8 + NEP9 + NEP10 + NEP11 + NEP13 + NEP15)
  data$birth <- data$Year - data$age
  data$AgeR = data$age/10 - 2
  data$AgeR2 = data$AgeR^2
  data <- data %>% mutate(cohort_cat = case_when(
    birth >= 1921 & birth <= 1930 ~ 1,
    birth >= 1931 & birth <= 1940 ~ 2,
    birth >= 1941 & birth <= 1950 ~ 3,
    birth >= 1951 & birth <= 1960 ~ 4,
    birth >= 1961 & birth <= 1970 ~ 5,
    birth >= 1971 & birth <= 1980 ~ 6,
    birth >= 1981 & birth <= 1990 ~ 7,
    birth >= 1991 & birth <= 2000 ~ 8,
    birth >= 2001 & birth <= 2010 ~ 9,
    TRUE ~ NA_integer_  
  ))
  data <- data[!is.na(data$cohort_cat), ] 
  data$cohort_cat <- as.factor(data$cohort_cat)
  # two cohorts with small sample sizes (i.e., 1921-1930 and 2001-2010) were excluded
  data <- data %>%
    filter(cohort_cat %in% 2:8) %>%
    mutate(cohort_cat = factor(cohort_cat))
  data <- data %>% mutate_at(vars(sex,Year,gov_satisfaction,happiness), as.factor)
  return(data)
})

#### 3.2.2.2 Perform Coarsened Exact Matching ====
matstrict <- cem("Year", datalist=data_list_modified, drop=c('NEP1','NEP2','NEP3','NEP4','NEP5','NEP6','NEP7','NEP8','NEP9',
                                                             'NEP10','NEP11','NEP12','NEP13','NEP14','NEP15',
                                                             'age','birth','Year','NEP','CNEP',
                                                             'AgeR','AgeR2','cohort_cat'))
matstrict
matched_data_list <- list()
for (i in seq_along(data_list)) {
  match_col <- paste0("match", i)
  matched_data_list[[i]] <- data_list_modified[[i]][matstrict[[match_col]]$matched, ]
}

#### 3.2.2.3 Construct the HAPC model for each matched dataset ====
result_CNEP_sensitive <- data.frame()
for (data in matched_data_list) {
  result_tem <- data.frame() 
  
  fullmodel <- lmer(CNEP ~ AgeR + AgeR2 + sex + happiness + gov_satisfaction + education + household_income +
                      + (1 | cohort_cat) + (1 | Year),data = data)
  nullmodel <- glm(CNEP ~ sex + happiness + gov_satisfaction + education + household_income ,
                   data = data)
  summary_fullmodel <- summary(fullmodel)
  
  # anova_p
  anova <- anova(fullmodel,nullmodel)
  result_tem[1,'anova_p'] <- anova$`Pr(>Chisq)`[2]
  
  # variance_cohort
  random_effects_var <- as.data.frame(VarCorr(fullmodel))
  result_tem[1,'variance_cohort'] <- random_effects_var[1,4]
  
  # cohort_p
  significant <- ranova(fullmodel)
  result_tem[1,'cohort_p'] <- significant$`Pr(>Chisq)`[2]
  
  # variance_period
  result_tem[1,'variance_period'] <- random_effects_var[2,4]
  
  # period_p
  result_tem[1,'period_p'] <- significant$`Pr(>Chisq)`[3]
  
  # Residual
  result_tem[1,'Residual'] <- random_effects_var[3,4]
  
  # Cohort_effect1-9
  Cohort_effect <- ranef(fullmodel)$cohort_cat
  result_tem[1,'Cohort_effect1931-1940'] <- Cohort_effect[1,1]
  result_tem[1,'Cohort_effect1941-1950'] <- Cohort_effect[2,1]
  result_tem[1,'Cohort_effect1951-1960'] <- Cohort_effect[3,1]
  result_tem[1,'Cohort_effect1961-1970'] <- Cohort_effect[4,1]
  result_tem[1,'Cohort_effect1971-1980'] <- Cohort_effect[5,1]
  result_tem[1,'Cohort_effect1981-1990'] <- Cohort_effect[6,1]
  result_tem[1,'Cohort_effect1991-2000'] <- Cohort_effect[7,1]
  
  # Period_effect2003-2021 
  Period_effect <- ranef(fullmodel)$Year
  result_tem[1,'Period_effect2003'] <- Period_effect[1,1]
  result_tem[1,'Period_effect2010'] <- Period_effect[2,1]
  result_tem[1,'Period_effect2021'] <- Period_effect[3,1]
  
  # age_effect & age_p
  fixed_effect <- as.data.frame(summary_fullmodel$coefficients)
  result_tem[1,'age_effect'] <- fixed_effect[2,1]
  result_tem[1,'age_p'] <- fixed_effect[2,5]
  
  # age2_effect & age2_p
  result_tem[1,'age2_effect'] <- fixed_effect[3,1]
  result_tem[1,'age2_p'] <- fixed_effect[3,5]
  
  # 控制变量
  result_tem[1,'female_effect'] <- fixed_effect[4,1]
  result_tem[1,'female_p'] <- fixed_effect[4,5]
  result_tem[1,'happiness.2_effect'] <- fixed_effect[5,1]
  result_tem[1,'happiness.2_p'] <- fixed_effect[5,5]
  result_tem[1,'happiness.3_effect'] <- fixed_effect[6,1]
  result_tem[1,'happiness.3_p'] <- fixed_effect[6,5]
  result_tem[1,'happiness.4_effect'] <- fixed_effect[7,1]
  result_tem[1,'happiness.4_p'] <- fixed_effect[7,5]
  result_tem[1,'happiness.5_effect'] <- fixed_effect[8,1]
  result_tem[1,'happiness.5_p'] <- fixed_effect[8,5]
  result_tem[1,'satisfaction.2_effect'] <- fixed_effect[9,1]
  result_tem[1,'satisfaction.2_p'] <- fixed_effect[9,5]
  result_tem[1,'satisfaction.3_effect'] <- fixed_effect[10,1]
  result_tem[1,'satisfaction.3_p'] <- fixed_effect[10,5]
  result_tem[1,'education_effect'] <- fixed_effect[11,1]
  result_tem[1,'education_p'] <- fixed_effect[11,5]
  result_tem[1,'household_income_effect'] <- fixed_effect[12,1]
  result_tem[1,'household_income_p'] <- fixed_effect[12,5]
  
  #R2
  r2 <- r.squaredGLMM(fullmodel)
  result_tem[1,'r2m'] <- r2[,1]
  result_tem[1,'r2c'] <- r2[,2]
  
  result_CNEP_sensitive <- rbind(result_CNEP_sensitive, result_tem)
}
# write.csv(result_CNEP_sensitive[,39:40], file = 'result_CNEP_sensitive_remove_R2.csv', row.names = FALSE)
# write.csv(result_CNEP_sensitive[,1:38], file = 'result_CNEP_sensitive_remove.csv', row.names = FALSE)
# These two combined form Supplementary Table 9.
# For convenience in plotting, store them in two separate files here, and do not consider R2 when plotting.

#### 3.2.2.4 Supplementary Fig 3 ====
result_CNEP_sensitive <- read.csv('result_CNEP_sensitive_remove.csv')
# The result_CNEP_sensitive_remove.csv in section 3.2.2.3 (the version without R2)

t.test(result_CNEP_sensitive$cohort_p, mu = 0.05, alternative = "less") #p-value = 0.01633 *
t.test(result_CNEP_sensitive$period_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive$age_p, mu = 0.05, alternative = "less") #p-value < 2.2e-16 *
t.test(result_CNEP_sensitive$age2_p, mu = 0.05, alternative = "less") #p-value < 2.2e-16 *
t.test(result_CNEP_sensitive$female_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive$happiness.2_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive$happiness.3_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive$happiness.4_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive$happiness.5_p, mu = 0.05, alternative = "less") #p-value = 1
t.test(result_CNEP_sensitive$satisfaction.2_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive$satisfaction.3_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive$education_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***
t.test(result_CNEP_sensitive$household_income_p, mu = 0.001, alternative = "less") #p-value < 2.2e-16 ***

names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "age_effect"] <- "Age"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "age2_effect"] <- "Age-quadratic"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "female_effect"] <- "Female(Ref:male)"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "happiness.2_effect"] <- "unhappy(Ref:Extremely unhappy)"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "happiness.3_effect"] <- "Neutral(Ref:Extremely unhappy)"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "happiness.4_effect"] <- "happy(Ref:Extremely unhappy)"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "happiness.5_effect"] <- "Extremely happy(Ref:Extremely unhappy)"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "satisfaction.2_effect"] <- "Neutral(Ref:dissatisfied)"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "satisfaction.3_effect"] <- "Satisfied(Ref:dissatisfied)"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "education_effect"] <- "Education"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "household_income_effect"] <- "Household income"

names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Cohort_effect1931.1940"] <- "Cohort 1931-1940"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Cohort_effect1941.1950"] <- "Cohort 1941-1950"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Cohort_effect1951.1960"] <- "Cohort 1951-1960"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Cohort_effect1961.1970"] <- "Cohort 1961-1970"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Cohort_effect1971.1980"] <- "Cohort 1971-1980"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Cohort_effect1981.1990"] <- "Cohort 1981-1990"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Cohort_effect1991.2000"] <- "Cohort 1991-2000"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Period_effect2003"] <- "Period 2003"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Period_effect2010"] <- "Period 2010"
names(result_CNEP_sensitive)[names(result_CNEP_sensitive) == "Period_effect2021"] <- "Period 2021"

CNEP_order <- c("Age", "Age-quadratic",'Female(Ref:male)','unhappy(Ref:Extremely unhappy)','Neutral(Ref:Extremely unhappy)',
                'happy(Ref:Extremely unhappy)','Extremely happy(Ref:Extremely unhappy)','Neutral(Ref:dissatisfied)','Satisfied(Ref:dissatisfied)',
                'Education','Household income','Period 2003',"Period 2010","Period 2021",
                "Cohort 1931-1940","Cohort 1941-1950","Cohort 1951-1960",
                "Cohort 1961-1970","Cohort 1971-1980","Cohort 1981-1990","Cohort 1991-2000")
result_CNEP_sensitive <- result_CNEP_sensitive %>% select(all_of(CNEP_order))
CNEP <- melt(result_CNEP_sensitive)
colnames(CNEP) <- c('sample','value')
CNEP$sample <- as.character(CNEP$sample)
CNEP <- CNEP[!grepl("_p$", CNEP$sample), ]
CNEP <- CNEP[!grepl("^variance", CNEP$sample), ]
CNEP <- CNEP[!grepl("^Resid", CNEP$sample), ]
CNEP$group <- ifelse(CNEP$sample %in% c("Age", "Age-quadratic",'Female(Ref:male)','unhappy(Ref:Extremely unhappy)','Neutral(Ref:Extremely unhappy)',
                                        'happy(Ref:Extremely unhappy)','Extremely happy(Ref:Extremely unhappy)','Neutral(Ref:dissatisfied)','Satisfied(Ref:dissatisfied)',
                                        'Education','Household income',"Period 2003","Period 2010","Period 2021"), "fix", "random")
CNEP <- CNEP[rev(seq_len(nrow(CNEP))), ]
unique_strings <- unique(unlist(CNEP$sample))
CNEP$sample <- factor(CNEP$sample,levels = unique_strings)
CNEP$group <- as.factor(CNEP$group)
CNEP$facet <- rep("APC effect",times=2625)

CNEP <- CNEP %>%
  mutate(facet = case_when(
    grepl("^Cohort", sample) ~ "Cohort effect *",
    grepl("^Period", sample) ~ "Period effect ***",
    TRUE ~  "Fixed effect"  
  ))

CNEPage <- CNEP[1251:2625,]
CNEPcohort <- CNEP[1:875,]
CNEPperiod <- CNEP[876:1250,]

##### 3.2.2.4.1 Age Part ====
p <- ggplot(CNEPage,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPage %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apcage <- ggplot(CNEPage,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#459f81"))+
  scale_color_manual(values = c("#459f81"))+
  scale_y_continuous(limits = c(-2, 2))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Coefficient",x=NULL)+
  annotate("rect", xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  annotate('text', label = '*', x =11, y =0.65, angle=-90, size =8,color="black") + 
  annotate('text', label = '*', x =10, y =-0.2, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =9, y =-0.55, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =4, y =-1.25, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =3, y =-0.9, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =2, y =0.45, angle=-90, size =8,color="black") + 
  annotate('text', label = '***', x =1, y =1.3, angle=-90, size =8,color="black") +
  facet_grid(~ facet)
p_apcage

##### 3.2.2.4.2 Cohort Part ====
CNEPcohort$sample <- gsub("^Cohort ", "", CNEPcohort$sample)
CNEPcohort <- CNEPcohort %>%
  mutate(sample = factor(sample, levels = c( "1991-2000",'1981-1990','1971-1980','1961-1970',
                                             '1951-1960','1941-1950','1931-1940')))

p <- ggplot(CNEPcohort,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPcohort %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean

p_apccohort <- ggplot(CNEPcohort,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#3769AE"))+
  scale_color_manual(values = c("#3769AE"))+
  scale_y_continuous(limits = c(-0.4, 0.4))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Random Intercept",x=NULL)+
  annotate("rect", xmin = 0, xmax = 7.6, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  facet_grid(~ facet)
p_apccohort

##### 3.2.2.4.3 Period Part ====
CNEPperiod$sample <- gsub("^Period ", "", CNEPperiod$sample)
CNEPperiod <- CNEPperiod %>%
  mutate(sample = factor(sample, levels = c("2021", "2010",'2003')))

p <- ggplot(CNEPperiod,aes(sample,value))+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.5,linewidth=0.5)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.5)
p
CNEPperiod %>% 
  group_by(sample) %>% 
  summarise(mean_value=mean(value)) %>%
  cbind(ggplot_build(p)$data[[1]]) -> mean
p_apcperiod <- ggplot(CNEPperiod,aes(sample,value))+
  geom_hline(yintercept = 0, linetype = 2, color = "grey60",linewidth=0.8)+
  stat_boxplot(aes(color=group),geom = "errorbar", width=0.3,size=0.6)+
  geom_boxplot(aes(fill=group,color=group),outlier.shape = 18,size=0.6)+
  geom_segment(mean,
               mapping=aes(x=xmin-0.25,xend=xmax+0.25,y=mean_value,yend=mean_value),
               color="white",size=0.5)+
  coord_flip()+
  scale_fill_manual(values = c("#3769AE"))+
  scale_color_manual(values = c("#3769AE"))+
  #y轴范围设置
  scale_y_continuous(limits = c(-2, 2))+
  scale_x_discrete(labels = function(x) gsub("Age-quadratic", "Age²", x)) +  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=11),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=13))+
  labs(y="Random Intercept",x=NULL)+
  annotate("rect", xmin = 0, xmax = 3.6, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="white") +
  facet_grid(~ facet)
p_apcperiod

##### 3.2.2.4.4 Supplementary Fig 3 ====
SupplementaryFig3 <- p_apcage / p_apcperiod / p_apccohort + 
  plot_layout(heights = c(11, 3, 7))
SupplementaryFig3