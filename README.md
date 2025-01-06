# RStudio_Projects

## Project Details 

**1- Predicting Chickenpox cases**

**Objective:** This project cleaned, transformed, and applied various predictive models on the chickenpox dataset to predict an outbreak or simply future cases. Several candidate models were trained and tested, then two promising models were selected based on forecasting accuarcy metrics. Those two models, namely ARIMA and Holt Winters model, were tuned and tested. Finally, using forecasting accuracy metrics, ARIMA was selected as the best model

Candidate models: Holt-Winter model, ARIMA, TBATS, and STL + ETS 

**Goal:** These predictions allow hospitals and other healthcare institutions to devise effective strategies that can assist them in combating the impacts of the anticipated influx of cases, mitigating the significant clinical expenses caused by the disease.

[Presentation Video](https://www.youtube.com/watch?v=CSUos9Z-z34&ab_channel=SheymaAbdikebir)

**Skills:** · R studio · Predictive Analytics · Data Visualization · Data Cleaning

**2- Customer Loyalty Prediction**

This project strived to provide insights into the specific predictors that the organization should optimize to improve customer loyalty. The dataset included various marketing-related predictors such as customer satisfaction, and negative publicity. etc. which were cleaned and transformed before feeding it to the Multiple Linear Regression model. 

Predictive accuracy was measured using 'Adjusted R-square ' and 'Mean Squared Error'.
Finally, repeated K-fold cross-validation was conducted to determine the performance of the model on new data.

**Skills:** · Predictive Analytics · Data Visualization · Data Transformation · Data Cleaning

**3- Amazon Musical Instruments Reviews NLP project**

The purpose of this project is to do topic modeling on the reviewText variable of Amazon’s musical instruments review data. This technique can allow Amazon to understand the different topics that the customers discussed in the reviews. Once topics are identified, the store can discover reviews that contain community guidelines violating keywords which can be removed to enhance customer experience. 

Data preparation and transformation techniques: 
- 7 missing rows were removed, LDA algorithm cannot function well with missing rows 
- Each review tokenized to single words by removing commas, periods, and symbols, to facilitate analysis 
- Each token standardized by being converted to lowercase, prevents case differences from impacting the analysis 
- Removed stopwords from the reviews 
- Numbers removed since initial results showed random numbers as heavily weighted
- Applied lemmatization as initial word cloud results showed words like guitar, guitars, string, and strings as heavily weighted

Selected parameters: 
- Chuck size = 10, reviews are generally short 
- Training iterations = 300, shows optimal results and avoids overfitting 
- Topics = 3, shows distinct themes    

Conclusion: 
- Guitar is the most frequent topic and token
- Word cloud showed that customers believe the instruments are of good quality, highlighting positive sentiment in general
- Reviews include the presence of community violating key terms. Hence, contextual analysis is required before removal
- Major topics: Guitars, Instrument’s sound-related topics, and instrument’s quality- related topics

**Skills:** · Natural Language Processing · Data Visualization · Data Transformation · Data Cleaning 

**4- Bank Customer Churn Analysis**

--Descriptive Analytics
Performed exploratory data analysis and cleaned a dataset with 10,001 customers, checked missing values and assessed whether to retain or remove outliers. Conducted data transformations to prepare the dataset for classification and used various correlation analysis methods (point biserial, phi, and Pearson’s) to identify relationships between variables.

--Predictive Analytics
Handled class imbalance by using oversampling, built and evaluated predictive models, including Logistic Regression, SVM, Decision Trees, and Random Forests. Scaled the data for Logistic Regression and SVM and assessed model performance using recall, F1-score, and accuracy metrics. Random Forest had the highest accuracy of approximately 94%, but I selected Logistic Regression with an accuracy of approximately 71% for its balance between complexity and interpretability, enabling actionable insights into the factors influencing customer churn.

--Prescriptive Analytics
Utilized Excel’s Solver and Solver Table to perform sensitivity analysis, minimizing churn by optimizing controllable decision variables while holding constants at median values. Then I proposed targeted marketing strategies based on model predictions and optimization results to enhance retention.

**Skills:** · Predictive Analytics · Prescriptive Analytics · Descriptive Analytics · Data Transformation

[Presentation Video](https://youtu.be/xtSQlXvSjfw)


