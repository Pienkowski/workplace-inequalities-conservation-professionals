#####################################################
######### Data exploration and pre-analysis #########
#####################################################

# The following explores the data used in the manuscript: Psychological distress and workplace risk inequalities among conservation professionals

### This script contains the following steps ###
# 1) Set the environment and load data 
# 2) Explore each variable and instrument 
# 3) Examine patterns of missing data 
# 4) Perform multiple imputation 
# 5) Post imputation manipulation 
# 6) Create a dataset to present summary statistics of the dataset  

######### 1) Set the environment and load data #########

# Load the data
DF_combined <- readRDS("DF_combined.rds") ### NB this dataset has been adapted by removing the organisational dummy variable.
# drop <- c("Organisation")
# DF_combined_redac = DF_combined[,!(names(DF_combined) %in% drop)]
# save(DF_combined_redac, file = "DF_combined_redac.Rdata") # Start 

# File containing variable information 
variables_df <- read.csv("Variables.csv")

### Set seed ###
set.seed(123)

# Load packages 
library(ggplot2)
library(psych)
library(lavaan)
library(semTools)
library(naniar)
library(mice)
library(fastDummies)
library(dplyr)

### Variables into groups ###
# K10
K10_col_name <- c("K10_1","K10_2","K10_3","K10_4","K10_5","K10_6","K10_7","K10_8","K10_9","K10_10")
K10_col_name_num <- paste0(c("K10_1","K10_2","K10_3","K10_4","K10_5","K10_6","K10_7","K10_8","K10_9","K10_10"),"_num")
K6_col_name <- c("K10_2","K10_4", "K10_5", "K10_8", "K10_9","K10_10")
K6_col_name_num <- paste0(c("K10_2","K10_4", "K10_5", "K10_8", "K10_9","K10_10"),"_num")

# SO
SO_col_name <- c("SO_1","SO_2","SO_3","SO_4","SO_5","SO_6","SO_7","SO_8","SO_9","SO_10","SO_11")
SO_col_name_num <- paste0(c("SO_1","SO_2","SO_3","SO_4","SO_5","SO_6","SO_7","SO_8","SO_9","SO_10","SO_11"),"_num")

# LOTR
LOTR_col_name <- c("LOTR_1","LOTR_2","LOTR_3","LOTR_4","LOTR_5","LOTR_6")
LOTR_pos <- c("LOTR_1","LOTR_3" ,"LOTR_6")
LOTR_neg <- c("LOTR_2","LOTR_4", "LOTR_5")
LOTR_col_name_num <- paste0(c("LOTR_1","LOTR_2","LOTR_3","LOTR_4","LOTR_5","LOTR_6"),"_num")
LOTR_pos_num <- paste0(c("LOTR_1","LOTR_3" ,"LOTR_6"),"_num")
LOTR_neg_num <- paste0(c("LOTR_2","LOTR_4", "LOTR_5"),"_num")

# PS
PS_col_name <- c("PS_1", "PS_2","PS_3")
PS_col_name_num <- paste0(c("PS_1", "PS_2","PS_3"),"_num")

# Social support
SS_col_name <- c("SS1", "SS2","SS3")
SS_col_name_num <- paste0(c("SS1", "SS2","SS3"),"_num")

# ERI
ERI_col_name <- c("ERI_1","ERI_2","ERI_3","ERI_4","ERI_5", "ERI_6","ERI_7","ERI_8","ERI_9","ERI_10", "ADD1","ADD2","ADD3","ADD4","ADD5")
ERI_normal<- c("ERI_1","ERI_2","ERI_3", "ADD1","ADD2" , "ADD3", "ERI_4", "ERI_8","ERI_9", "ERI_10", "ADD4","ADD5")
ERI_reverse <- c("ERI_5", "ERI_6","ERI_7")
ERI_col_name_num <- paste0(c("ERI_1","ERI_2","ERI_3","ERI_4","ERI_5", "ERI_6","ERI_7","ERI_8","ERI_9","ERI_10", "ADD1","ADD2","ADD3","ADD4","ADD5"),"_num")
ERI_normal_num<- paste0(c("ERI_1","ERI_2","ERI_3", "ADD1","ADD2" , "ADD3", "ERI_4", "ERI_8","ERI_9", "ERI_10", "ADD4","ADD5"), "_num")
ERI_reverse_num <- paste0(c("ERI_5", "ERI_6","ERI_7"), "_num")


######### 2) Explore each variable and instrument ######### 

###### Kessler 10 ######

### Plotting raw reponses ###
# Subset to K10 items 
K10_DF <- DF_combined[K10_col_name]

# create a DF with the proportions for the first K10 variable 
K10_prop <- data.frame(Variable = K10_col_name[1], cbind(t(prop.table(table(K10_DF[1])))))

# For loop appending the proportions for the remaining variables 
for (i in seq_along(2:(length(K10_col_name)))){
  K10_prop <- rbind(K10_prop, data.frame(Variable = K10_col_name[i+1], cbind(t(prop.table(table(K10_DF[i+1]))))))
}

# Rename K10 variable code to question
K10_prop$Variable <- subset(variables_df, Code %in% K10_col_name, select=c("Question"))$Question

# rename response levels
names(K10_prop) <- c("Variable", "None of the time", "A little of the time", "Some of the time","Most of the time", "All of the time")

# Plot the responses 
HH::likert(Variable~.,K10_prop, positive.order=TRUE,as.percent = TRUE,
           main="K10 responses",
           xlab="Percentage", ylab="Variable")

### Density ###
# Overall
ggplot(DF_combined, aes(x=K10_total)) + 
  geom_density()

# Save plot - jpeg
ggsave(
  "Figure_K10.jpeg",
  width = 90,
  height = 90,
  units = c("mm"),
  dpi = 800
)

# 20–24 are likely to have a mild level of distress
# 25–29 moderate level of distress
# 30–50 severe

# Nominal distress
round(table(DF_combined$K10_total < 20)/nrow(DF_combined),2)

# No distress
round(table(DF_combined$K10_total >= 30)/nrow(DF_combined),2)

# Between organisations 
ggplot(DF_combined, aes(x=K10_total, color=Organisation)) +
  geom_density()

# Between roles 
ggplot(DF_combined, aes(x=K10_total, color=position_simple)) +
  geom_density()

### K6 total 
DF_combined$K6_total <- rowSums(DF_combined[K6_col_name_num])
DF_combined$K6_total_scale <- scale(DF_combined$K6_total, center = T, scale = T)
DF_combined$K6_total_scale <- as.numeric(DF_combined$K6_total_scale)

# Association between K10 and K6 
cor.test(x=DF_combined$K10_total_scale, y=DF_combined$K6_total_scale, method = 'spearman')

###### LOTR ######

### Density ###
ggplot(DF_combined, aes(x=LOTR_total)) + 
  geom_density()

# Between organisations 
ggplot(DF_combined, aes(x=LOTR_total, color=Organisation)) +
  geom_density()

# Number of NA's 
table(is.na(DF_combined[,LOTR_col_name]))

###### Gender ###### 

### Table ###
table(DF_combined$gender)

# Recode 
DF_combined$gender_simple <- ifelse(DF_combined$gender == "Male", "Male", ifelse(
  DF_combined$gender == "Female", "Female", ifelse(
    DF_combined$gender %in% c("Prefer not to say", "Unknown", NA, "NA"), "Unknown", "ERROR"
  )))

###### Age ###### 

### Density ###
ggplot(DF_combined, aes(x=age_year)) + 
  geom_density()

# Between organisations 
ggplot(DF_combined, aes(x=age_year, color=Organisation)) +
  geom_density()

# Number of NAs
table(is.na(DF_combined$age_year))


###### Years in conservation ###### 

### Density ###
ggplot(DF_combined, aes(x=years_cons)) + 
  geom_density()

# Between organisations 
ggplot(DF_combined, aes(x=years_cons, color=Organisation)) +
  geom_density()

# Number of NAs
table(is.na(DF_combined$years_cons))


# Check for covariation between age and years of experience - moderate 
cor.test(x=DF_combined$age_year_scaled, y=DF_combined$years_cons_scaled, method = 'spearman')


###### Health ###### 

### Table ###
table(DF_combined$health)

# Number of NAs
table(is.na(DF_combined$health))

# Scale health
DF_combined$health_num_scaled <- scale(DF_combined$health_num, scale = T, center = T)

###### Social support ###### 

### Tables ###
# Overall
table(DF_combined$SS1)
table(DF_combined$SS2)
table(DF_combined$SS3)

# Between roles 
round(prop.table(table(DF_combined$SS1, DF_combined$position_simple)),2)
round(prop.table(table(DF_combined$SS2, DF_combined$position_simple)),2)
round(prop.table(table(DF_combined$SS3, DF_combined$position_simple)),2)


### Create composite variable for social support ### 
DF_combined$SS_comp <- rowMeans(DF_combined[SS_col_name_num]) 
DF_combined$SS_comp_scaled <- scale(DF_combined$SS_comp, scale = T, center = T)
DF_combined$SS_comp_scaled <- as.numeric(DF_combined$SS_comp_scaled)

###### ERI score ###### 
### Density ###
# Overall
ggplot(DF_combined, aes(x=ERI_n_num)) + 
  geom_density()

# Between organisations 
ggplot(DF_combined, aes(x=ERI_n_num, color=Organisation)) +
  geom_density()

# Between roles 
ggplot(DF_combined, aes(x=ERI_n_num, color=position_simple)) +
  geom_density()

###### Not feeling safe ###### 

### Tables ###

# Overall 
table(DF_combined$PS_1)
table(DF_combined$PS_1)
table(DF_combined$PS_1)

# Between roles 
round(prop.table(table(DF_combined$PS_1, DF_combined$position_simple)),2)
round(prop.table(table(DF_combined$PS_2, DF_combined$position_simple)),2)
round(prop.table(table(DF_combined$PS_3, DF_combined$position_simple)),2)

###### Position ###### 

### Table ###
# Simple 
table(DF_combined$position)

# All positions 
table(DF_combined$position_simple)


######  3) Examine patterns of missing data ###### 
### Removing observation that do not respond to a sufficient number of effort or reward items (thus generating NaN in the ERI calculation) ###
table(is.nan(DF_combined$ERI_n_num))
DF_combined_2 <- DF_combined[!is.nan(DF_combined$ERI_n_num),]
table(is.nan(DF_combined_2$ERI_n_num))

### Check complete cases for K10 ### 
table(is.na(DF_combined_2[K10_col_name_num]))


### Subset to data for use in the analysis ###
keep <- c("K10_1_num","K10_2_num","K10_3_num","K10_4_num", "K10_5_num","K10_6_num","K10_7_num","K10_8_num",
          "K10_9_num", "K10_10_num","PS_2_num","SS1_num","SS2_num"  ,"SS3_num","health_num", 
          "position_simple" ,  "gender_simple",  "years_cons_scaled",
          "age_year_scaled",  "ERI_o_scaled" ,  "ERI_n_scaled", "Organisation" )
DF1_analy <-  DF_combined_2[keep]

### Look at patterns of missing data ###
vis_miss(DF1_analy) + coord_flip()
gg_miss_upset(DF1_analy)

# Check the patterns of missingness in the variables to be imputed 
vis_miss(DF1_analy) + coord_flip()

# Save plot - jpeg
ggsave(
  "Figure_missing.jpeg",
  width = 100,
  height = 120,
  units = c("mm"),
  dpi = 800
)

# Save
save(DF_combined_2, file = "DF_combined_2.Rdata")


###### 4) Perform multiple imputation ######
### Subset to MI dataset ###
# I.e., removing un-scaled values, individual ERI elements, etc. 
imptuted_variables <- c("Location",  "K10_1_num","K10_2_num","K10_3_num","K10_4_num", "K10_5_num","K10_6_num","K10_7_num","K10_8_num",
          "K10_9_num", "K10_10_num","LOTR_total_scaled", "PS_1_num","PS_2_num","PS_3_num","SS1_num","SS2_num"  ,"SS3_num","health_num_scaled",
          "position_simple" , "education_simple", "gender_simple",  "years_cons_scaled",
          "age_year_scaled",   "ERI_n_scaled", "WH_scaled","Organisation" )
DF1_analy_imp <-  DF_combined_2[imptuted_variables]


### Specify the correct variable type ### 
# Ordinal variables 
ordinal_vars <- c("K10_1_num","K10_2_num","K10_3_num","K10_4_num", "K10_5_num","K10_6_num","K10_7_num","K10_8_num",
                  "K10_9_num", "K10_10_num", "PS_1_num","PS_2_num","PS_3_num","SS1_num","SS2_num"  ,"SS3_num","health_num_scaled")

# Numeric variables 
numeric_vars <- c("age_year_scaled", "years_cons_scaled", "WH_scaled", "ERI_n_scaled", "LOTR_total_scaled" )

# Factors 
factor_vars <- c("position_simple" , "education_simple", "gender_simple", "Location", "Organisation" )

# Exogenous variables 
exo_vars <- c("PS_1_num","PS_2_num","PS_3_num","SS1_num","SS2_num"  ,"SS3_num", "health_num_scaled")

# Did I miss any variable? 
setdiff(imptuted_variables, c(ordinal_vars,numeric_vars,factor_vars))
setdiff(c(ordinal_vars,numeric_vars,factor_vars),imptuted_variables)

### Function that ensures the data is in the correct format,
MI_function_1 <-function(DF1) {
  
  ### Correct data type ###
  # As ordinal
  DF1[, ordinal_vars] <- lapply(DF1[, ordinal_vars], as.ordered)
  
  # As numeric
  DF1[, numeric_vars] <- lapply(DF1[, numeric_vars], as.numeric)
  
  # As factor 
  DF1[, factor_vars] <- lapply(DF1[, factor_vars], as.factor)
  
  ### Number of imputed DF ###
  # 10 imputed DF
  N.Imp = 10
  
  ### determine the imputation method from the data type ##
  # The data type for each variable 
  str_out <- data.frame(capture.output(str(DF1)))
  
  # Delete the first row
  str_out <- data.frame(str_output = matrix(str_out[2:(nrow(str_out)),]))
  
  # Create a column that contain the appropriate model for each variable - this only works if the variable is of the correct data type in the first place 
  str_out$type <- ifelse(grepl("Ord.factor", str_out$str_output, fixed = TRUE)==T, "polr", 
                         ifelse(grepl("num", str_out$str_output, fixed = TRUE)==T, "pmm", 
                                ifelse(grepl("Factor", str_out$str_output, fixed = TRUE)==T, "polyreg",
                                       ifelse(grepl("int", str_out$str_output, fixed = TRUE)==T, "logreg", "ERROR"))))
  

  
  # Conduct the MI - with the number of datasets specified by N.Imp, and the estimation type specified by str_out$type (derived from the above)
  DF1_imp <- mice(DF1, m = N.Imp, method = str_out$type )
  
  # Print the first 50 logged events, if they occur 
  print(head(DF1_imp$loggedEvents, 50))
  
  # Return the imputed data
  return(DF1_imp)
}

### Conduct the imputation ### 
# Conduct the imputation using DF.DIS.1 (excluding country code)
mice.imp.DEEP_int <- MI_function_1(DF = DF1_analy_imp)

# Save the imputed data
save(mice.imp.DEEP_int, file = "mice.imp.DEEP_int.Rdata")
load("mice.imp.DEEP_int.Rdata") 


###


###### 5) Post imputation manipulation ######
### Function to add relivant variables 
MI_function_2 <-function(DF1, DF_imp) {
  
  ### Number of imputed DF ###
  # 10 imputed DF
  N.Imp = 10
  
  ### Extract each imputed dataset and perform additional manipulation ###
  # Create a list to store the imputed datasets 
  mice.imp <- list()
  
  # For i in each dataset
  for(i in 1:N.Imp) {
    
    ### Extract the imputed data
    mice.imp[[i]] <- mice::complete(DF_imp, action= i, inc=FALSE)
    

    ### Turn factor into series of binary variables and remove special characters 
    # Create the dummy columns 
    mice.imp[[i]] <- dummy_cols(mice.imp[[i]], select_columns = c("position_simple" , "education_simple", "gender_simple", "Location", "Organisation" ))
    
    # Remove spaces and other characters on column names
    colnames( mice.imp[[i]]) <- gsub(" ", "", colnames( mice.imp[[i]]), fixed = TRUE)
    colnames( mice.imp[[i]]) <- gsub("/", "", colnames( mice.imp[[i]]), fixed = TRUE)
    colnames( mice.imp[[i]]) <- gsub("-", "_", colnames( mice.imp[[i]]), fixed = TRUE)
    colnames( mice.imp[[i]]) <- gsub("&", "_", colnames( mice.imp[[i]]), fixed = TRUE)
    
    # Add back the ERI variables 
    mice.imp[[i]] <- cbind(mice.imp[[i]], DF1[,ERI_col_name_num])

    # Add back the unscaled items 
    mice.imp[[i]]$years_cons <- DF1$years_cons
    mice.imp[[i]]$age_year <- DF1$age_year
    mice.imp[[i]]$WH <- DF1$WH
    mice.imp[[i]]$LOTR_totaL  <- DF1$LOTR_total
    mice.imp[[i]]$health_num <- DF1$health_num
    mice.imp[[i]]$ERI_o_scaled <- DF1$ERI_o_scaled
            
    # # Convert exogenous ordinal variables to numeric 
    if(i == 1) { print(str(mice.imp[[i]][exo_vars][1]))} # Check the 
    if(i == 1) { print(table(mice.imp[[i]][exo_vars][1]))} # Check the 
    mice.imp[[i]][exo_vars] <- apply(mice.imp[[i]][exo_vars], 2, as.character)
    mice.imp[[i]][exo_vars] <- apply(mice.imp[[i]][exo_vars], 2, as.numeric)
    if(i == 1) { print(str(mice.imp[[i]][exo_vars][1]))}
    if(i == 1) { print(table(mice.imp[[i]][exo_vars][1]))}
    
    
  }
  
  # Return the manipulated DF 
  return(mice.imp)
}

# Implement the function 
mice.imp.DEEP <- MI_function_2(DF1 = DF_combined_2, DF_imp = mice.imp.DEEP_int)

# Save the data 
save(mice.imp.DEEP, file = "mice.imp.DEEP.Rdata")

###


#########  6) Create a dataset to present summary statistics of the dataset ######### 
# Select variables to include in the summary statistics table 
# Manuscript variables 
man_vars <- c("K10_total", "K6_total", "LOTR_total", "gender_simple", "age_year", "years_cons", 
              "health", "SS_comp_scaled", "ERI_o_num", "ERI_n_num",
              "position_simple", "PS_2",  "Organisation")

# Subset 
man_DF <- DF_combined_2[man_vars]

### Removing observation that do not respond to a sufficient number of effort or reward items (thus generating NaN in the ERI calculation) ###
table(is.nan(man_DF$ERI_n_num))
man_DF <- man_DF[!is.nan(man_DF$ERI_n_num),]
table(is.nan(man_DF$ERI_n_num))

# Drop observations with no K10 scores (should drop none)
man_DF <- man_DF[!is.na(man_DF$K10_total),]

### Rename variables ### 
man_DF %>% 
  rename(
    "Kessler-10 score" = "K10_total",
    "Kessler-6 score" = "K6_total",
    "LOTR score" = "LOTR_total", 
    "Gender" = "gender_simple", 
    "Age" = "age_year", 
    "Years in conservation" = "years_cons",
    "Physical health" = "health",
    "Social support" = "SS_comp_scaled",
    "Position" = "position_simple",
    "ERI (original)" = "ERI_o_num", 
    "ERI (adapted)" = "ERI_n_num", 
    "Dangerous situations" = "PS_2",
    "Organisation" = "Organisation"
  ) -> man_DF


# Save
save(man_DF, file = "man_DF.Rdata") # Start 

