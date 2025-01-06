setwd("C:/Users/shaim/Downloads/R datasets")
library(readr)
churn = read_csv('Churn_Modelling.csv')
head(churn)

churn2 <- churn[, -c(1, 2)] #removes customer id, and surname

#EDA 

#check for missing values 
colSums(is.na(churn2)) #no missing values
#check for duplicates 
duplicates <- churn2[duplicated(churn2), ] #no duplicates 

#High-leverage data/outlier detection

# CreditScore 
##method 1 IQR

churn_adj <- churn2[, -c( 2,3,7,8,9,11)]

Q1 <- quantile(churn_adj$CreditScore, 0.25)
Q3 <- quantile(churn_adj$CreditScore, 0.75)

IQR <- Q3 - Q1
multiplier <- 1.5
lower_bound <- Q1 - (IQR * multiplier)
upper_bound <- Q3 + (IQR * multiplier)

# Create a dataframe for boxplot
C_score_df <- data.frame(value = churn_adj$CreditScore)

# Outlier values
cs_outliers_IQR <- churn_adj$CreditScore[churn_adj$CreditScore < lower_bound | churn_adj$CreditScore > upper_bound]

## Method 2- Z-score
cs_score <- scale(churn_adj$CreditScore)

outliers <- abs(cs_score)>3
  
cs_outliers_z <-churn_adj$CreditScore[outliers]

#Age

#method 1- IQR

Q1 <- quantile(churn_adj$Age, 0.25)
Q3 <- quantile(churn_adj$Age, 0.75)

IQR <- Q3 - Q1
multiplier <- 1.5
lower_bound <- Q1 - (IQR * multiplier)
upper_bound <- Q3 + (IQR * multiplier)

# Create a dataframe for boxplot
Age_df <- data.frame(value = churn_adj$Age)

# Outlier values
Age_outliers <- churn_adj$Age[churn_adj$Age < lower_bound | churn_adj$Age > upper_bound]

#Method 2 - Z scores

Age_score <- scale(churn_adj$Age)

outliers <- abs(Age_score)>3

outlier_values <-churn_adj$Age[outliers]

#Tenure

#method 1- IQR

Q1 <- quantile(churn_adj$Tenure, 0.25)
Q3 <- quantile(churn_adj$Tenure, 0.75)

IQR <- Q3 - Q1
multiplier <- 1.5
lower_bound <- Q1 - (IQR * multiplier)
upper_bound <- Q3 + (IQR * multiplier)

# Create a dataframe for boxplot
Tenure_df <- data.frame(value = churn_adj$Tenure)

# Outlier values
Tenure_outliers <- churn_adj$Tenure[churn_adj$Tenure < lower_bound | churn_adj$Tenure > upper_bound]

#method 2- Z scores
Tenure_score <- scale(churn_adj$Tenure)

outliers <- abs(Tenure_score)>3

outlier_values <-churn_adj$Tenure[outliers]

# Balance

#method 1- IQR

Q1 <- quantile(churn_adj$Balance, 0.25)
Q3 <- quantile(churn_adj$Balance, 0.75)

IQR <- Q3 - Q1
multiplier <- 1.5
lower_bound <- Q1 - (IQR * multiplier)
upper_bound <- Q3 + (IQR * multiplier)

# Create a data frame for box plot
Balance_df <- data.frame(value = churn_adj$Balance)

# Outlier values
Balance_outliers <- churn_adj$Balance[churn_adj$Balance < lower_bound | churn_adj$Balance > upper_bound]

#method 2- Z scores
Balance_score <- scale(churn_adj$Balance)

outliers <- abs(Balance_score)>3

outlier_values <-churn_adj$Balance[outliers]

#Estimated Salary

#method 1- IQR

Q1 <- quantile(churn_adj$EstimatedSalary, 0.25)
Q3 <- quantile(churn_adj$EstimatedSalary, 0.75)

IQR <- Q3 - Q1
multiplier <- 1.5
lower_bound <- Q1 - (IQR * multiplier)
upper_bound <- Q3 + (IQR * multiplier)

# Create a data frame for box plot
EstimatedSalary_df <- data.frame(value = churn_adj$EstimatedSalary)

# Outlier values
EstimatedSalary_outliers <- churn_adj$EstimatedSalary[churn_adj$EstimatedSalary < lower_bound | churn_adj$EstimatedSalary > upper_bound]

#method 2- Z scores
EstimatedSalary_score <- scale(churn_adj$EstimatedSalary)

outliers <- abs(EstimatedSalary_score)>3

outlier_values <-churn_adj$EstimatedSalary[outliers]

######box plots######

par(mfrow = c(1, 5))

# Credit score
bplot_cs <- boxplot(C_score_df,
        main = "Box Plot for CreditScore",
        xlab = "CreditScore",
        ylab = "Values",
        col = c("lightblue"),
        border = "black",
        notch = TRUE)

# Age
bplot_age <- boxplot(Age_df,
        main = "Box Plot for Age",
        xlab = "Age",
        ylab = "Values",
        col = c("lightblue"),
        border = "black",
        notch = TRUE)

#Tenure
bplot_Tenure <- boxplot(Tenure_df,
        main = "Box Plot for Tenure",
        xlab = "Tenure",
        ylab = "Values",
        col = c("lightblue"),
        border = "black",
        notch = TRUE)

#Balance
bplot_Balance <- boxplot(Balance_df,
        main = "Box Plot for Balance",
        xlab = "Balance",
        ylab = "Values",
        col = c("lightblue"),
        border = "black",
        notch = TRUE)


#EstimatedSalary
bplot_EstimatedSalary <- boxplot(EstimatedSalary_df,
        main = "Box Plot for EstimatedSalary",
        xlab = "EstimatedSalary",
        ylab = "Values",
        col = c("lightblue"),
        border = "black",
        notch = TRUE)

######################deciding to remove or keep outliers################

# Filter outlier and typical credit scores
outlier_scores <- churn_adj %>%
  filter(CreditScore >= 350 & CreditScore <= 382)

typical_scores <- churn_adj %>%
  filter(!(CreditScore >= 350 & CreditScore <= 382))

outlier_stats <- outlier_scores %>%
  summarise(
    Mean_Age = mean(Age, na.rm = TRUE),
    #Median_Age = median(Age, na.rm = TRUE),
    Mean_Tenure = mean(Tenure, na.rm = TRUE),
    #Median_Tenure = median(Tenure, na.rm = TRUE),
    Mean_Balance = mean(Balance, na.rm = TRUE),
    #Median_Balance = median(Balance, na.rm = TRUE),
    Mean_EstimatedSalary = mean(EstimatedSalary, na.rm = TRUE),
    #Median_EstimatedSalary = median(EstimatedSalary, na.rm = TRUE)
  )

# Summary statistics for typical scores
typical_stats <- typical_scores %>%
  summarise(
    Mean_Age = mean(Age, na.rm = TRUE),
    #Median_Age = median(Age, na.rm = TRUE),
    Mean_Tenure = mean(Tenure, na.rm = TRUE),
    #Median_Tenure = median(Tenure, na.rm = TRUE),
    Mean_Balance = mean(Balance, na.rm = TRUE),
    #Median_Balance = median(Balance, na.rm = TRUE),
    Mean_EstimatedSalary = mean(EstimatedSalary, na.rm = TRUE),
  )

#remove outliers in credit scores and update the table 
churn_updated <- churn2 %>%
  filter(!(CreditScore >= 350 & CreditScore <= 382))

###distribution charts####

par(mfrow = c(1, 2))
#credit score 
hist(churn_updated$CreditScore, xlab = 'Credit Score', main = 'Distribution of Credit Score')

#age 
hist(churn_updated$Age, ylim = c(0, 2500), xlab = 'Age', main = 'Distribution of Age')

#Geography 
library(ggplot2)

ggplot(churn_updated, aes(x = Geography)) +
  geom_bar(fill = "grey") +
  labs(title = "Frequency of Customers by Geography",
       x = "Geography",
       y = "Count") +
  theme_minimal()

#Gender 
ggplot(churn_updated, aes(x = Gender)) +
  geom_bar(fill = "grey") +
  labs(title = "Frequency of Customers by Gender",
       x = "Gender",
       y = "Count") +
  scale_y_continuous(breaks = seq(0, 6000, by = 1000)) +  
  theme_minimal ()

#balance 
hist(churn_updated$Balance, ylim = c(0,4000), xlab = 'Balance', main = 'Distribution of Balance')

#Tenure 
hist(churn_updated$Tenure, ylim = c(0, 1200), xlab = 'Tenure', main = 'Distribution of Tenure')

#Estimated Salary 

hist(churn_updated$EstimatedSalary, xlab = 'EstimatedSalary', main = 'Distribution of EstimatedSalary')

#NumOfProducts 
ggplot(churn_updated, aes(x = NumOfProducts)) +
  geom_bar(fill = "grey") +
  labs(title = "Frequency of Customers by NumOfProducts",
       x = "NumOfProducts",
       y = "Count") +
  theme_minimal()

#HasCrCard 

churn_updated$HasCrCard <- as.factor(churn_updated$HasCrCard)

ggplot(churn_updated, aes(x = HasCrCard)) +
  geom_bar(fill = "grey") +
  labs(title = "Frequency of Customers by Creditcard availability",
       x = "Creditcard",
       y = "Count") +
  theme_minimal()

#IsActiveMember

churn_updated$IsActiveMember <- as.factor(churn_updated$IsActiveMember)

ggplot(churn_updated, aes(x = IsActiveMember)) +
  geom_bar(fill = "grey") +
  labs(title = "Frequency of Customers by Member status",
       x = "IsActiveMember",
       y = "Count") +
  theme_minimal()

#Exited

churn_updated$Exited <- as.factor(churn_updated$Exited)

ggplot(churn_updated, aes(x = Exited)) +
  geom_bar(fill = "grey") +
  labs(title = "Frequency of Customers by Exited",
       x = "Exited",
       y = "Count") +
  theme_minimal()


#Exited and IsActiveMember

ggplot(churn_updated, aes(x = Exited, fill = factor(IsActiveMember))) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("lightblue", "lightgreen"), 
                    name = "Is Active Member",
                    labels = c("No", "Yes")) +
  labs(title = "Frequency of Customers by Exited & Membership Status",
       x = "Exited",
       y = "Count") +
  theme_minimal()


#Exited and NumOfProducts

ggplot(churn_updated, aes(x = Exited, fill = factor(NumOfProducts))) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2", 
                    name = "Number of Products") +
  labs(title = "Frequency of Customers by Exited & Number of Products",
       x = "Exited",
       y = "Count") +
  theme_minimal()

#Exited and HasCrCard
ggplot(churn_updated, aes(x = Exited, fill = factor(HasCrCard))) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set5", 
                    name = "HasCrCard") +
  labs(title = "Frequency of Customers by Exited & Number of Products",
       x = "Exited",
       y = "Count") +
  theme_minimal()

######descriptive stats 

#remove all non numeric columns from the data frame to calculate descriptive stats

churn_updated_2 = churn_updated[, -c(2,3,8,9,11)]

describe(churn_updated_2)

#####Correlation Analysis#####

#find correlation between binary independent and binary dependent variables
#using phi coefficient 

library(psych)

#correlation between HasCrCard column and Exited 

phi_result1 <- phi(table(churn_updated$HasCrCard, churn_updated$Exited))

# Print the result
print(phi_result1) # -0.01, results show almost no correlation 

#correlation between IsActiveMember column and Exited 

phi_result2 <- phi(table(churn_updated$IsActiveMember, churn_updated$Exited))

# Print the result
print(phi_result2) #results show little correlation (-0.16)

#convert gender into binary variable 

churn_updated$Gender<-ifelse(churn_updated$Gender=="Male",1,0) #male 1 and female 0

#calculate phi correlation for gender and exited

phi_result3 <- phi(table(churn_updated$Gender, churn_updated$Exited))

# Print the result
print(phi_result3) #results show little correlation (-0.11)

#convert Geography to binary variable. 

churn_updated$Spain <- ifelse(churn_updated$Geography == "Spain", 1, 0)
churn_updated$Germany <- ifelse(churn_updated$Geography == "Germany", 1, 0)

# Drop the original geography column

churn_new <- churn_updated[, -c( 2)]

#calculate correlation between Spain and Exited
phi_result4 <- phi(table(churn_new$Spain, churn_new$Exited))

print(phi_result4) #results show almost no correlation (-0.05)

phi_result5 <- phi(table(churn_new$Germany, churn_new$Exited))

print(phi_result5) #results show little correlation (0.17)

#Create correlation chart for the numerical variables with Exited

churn_num = churn_new[, -c(2, 7,8, 11,12)] #removes binary independent variables

#convert exited to numeric variable 
churn_num$Exited <- as.numeric(as.character(churn_num$Exited))

#use point-biserial correlation for credit score and exited 
churn_cor1= cor.test(churn_num$CreditScore, churn_num$Exited) #-0.01838391 

#use point-biserial correlation for tenure and exited 
churn_cor2= cor.test(churn_num$Tenure, churn_num$Exited) #-0.01308019 

#use point-biserial correlation for Age and exited 
churn_cor3= cor.test(churn_num$Age, churn_num$Exited) #0.2851342  

#use point-biserial correlation for Balance and exited 
churn_cor4= cor.test(churn_num$Balance, churn_num$Exited) #0.1185169 

#use point-biserial correlation for NumOfProducts and exited 
churn_cor5= cor.test(churn_num$NumOfProducts, churn_num$Exited) #-0.04807463  

#use point-biserial correlation for EstimatedSalary and exited 
churn_cor6= cor.test(churn_num$EstimatedSalary, churn_num$Exited) #0.01017127  

#create a data frame grouping the calculated correlations for each numeric 
#variables with exited

Cor_values <- c(-0.01838391, -0.01308019, 0.2851342, 0.1185169 , -0.04807463, 0.01017127)
Predictors <- c("CreditScore", "Tenure", "Age", "Balance", "NumOfProducts", "EstimatedSalary")

corr_df<- data.frame(Exited = Cor_values, row.names = Predictors)

#convert the dataframe to matrix 

corr_matrix <- as.matrix(corr_df)

corrplot(corr_matrix)

#create a separate correlation chart between binary and numeric predictors 

#create a df for binary pred

churn_bin = churn_new[, -c(1,3,4,5,6,9,10)]

churn_num_pred <- churn_new[, -c(2,7,8,9,10,11,12)] 
churn_pred <- cbind(churn_bin, churn_num_pred)

churn_pred$Gender<- as.numeric(as.character(churn_pred$Gender))
churn_pred$HasCrCard<- as.numeric(as.character(churn_pred$HasCrCard))
churn_pred$IsActiveMember<- as.numeric(as.character(churn_pred$IsActiveMember))
churn_pred$Spain<- as.numeric(as.character(churn_pred$Spain))
churn_pred$Germany<- as.numeric(as.character(churn_pred$Germany))


churn_pred_cor <- cor(churn_pred)

corrplot(churn_pred_cor)

#change each binary data to numeric then calculate its correlation 

#Credit card and Gender

churn_cor2= cor(churn_num2) 
churn_num$Exited <- as.numeric(as.character(churn_num$Exited))


corrplot(churn_cor2)

#correlation chart for binary variables with Exited 

#create a data frame grouping the calculated correlations for each binary
#variables with exited

Cor_values <- c(-0.01, -0.16, -0.11, -0.05, 0.17)
Predictors <- c("HasCrCard", "IsActiveMember", "Gender", "Spain", "Germany")

corr_df<- data.frame(Exited = Cor_values, row.names = Predictors)

#convert the dataframe to matrix 

corr_matrix <- as.matrix(corr_df)

corrplot(corr_matrix)


#export the updated dataframe to excel 
install.packages("writexl")
library(writexl)
write_xlsx(
  churn_new,
  "C:/Users/shaim/Downloads/R datasets\\churn_new_outlier.xlsx",
col_names = TRUE)

