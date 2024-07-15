install.packages("readr")
library(readr)
setwd('C:/Users/shaim/Downloads/R datasets')
Cpox <- read.csv("hungary_Chickenpox_data.csv")
head(Cpox)
tail(Cpox)
colSums(is.na(Cpox)) #no missing values
duplicates <- Cpox[duplicated(Cpox), ] 
summary(Cpox)

Cpox$Date <- as.POSIXct(Cpox$Date, format = "%d/%m/%Y")
#check data type
data_types <- sapply(Cpox, class) #shows data type for each variable
print(data_types)

# create a new total cases column 
Cpox$total_cases <- rowSums(Cpox[, -1])


#create time object 
Total.ts <- ts(Cpox$total_cases, start = c(2005, 1), end = c(2014, 52), freq = 52)
#plot the time series for total cases 
plot(Cpox$Date, Cpox$total_cases, type = "l", xlab = "Date", ylab = "Total Cases", main = "Time Series plot of Total Chickenpox Cases in Hungary")

#create histogram for budapest 

hist(Cpox$BUDAPEST, 
     main = "Distribution of Cases in Budapest",  # Title of the histogram
     xlab = "Values",                 # Label for the x-axis
     ylab = "Frequency",              # Label for the y-axis
     col = "skyblue",                 # Color of the bars
     border = "black") 

#create histogram for Total cases

hist(Cpox$total_cases, 
     main = "Distribution of Total Cases in Hungary",  # Title of the histogram
     xlab = "Values",                 # Label for the x-axis
     ylab = "Frequency",              # Label for the y-axis
     col = "skyblue",                 # Color of the bars
     border = "black")

#create time object 
Bud.ts <- ts(Cpox$BUDAPEST, start = c(2005, 1),end = c(2014, 52), freq = 52)
#plot the time series for budapest
plot(Cpox$Date, Cpox$BUDAPEST, type = "l", xlab = "Date", ylab = "Cases", main = "Time Series plot of Chickenpox Cases in BUDAPEST")


library(forecast)

# Decompose the time series and plot it 
decomposed <- decompose(Bud.ts, type= 'additive')
plot(decomposed)

#Visualize trend through centered moving average 

library(zoo)

ma.centered <- ma(Bud.ts, order = 52)
plot(Bud.ts, ylab = "cases", xlab = "Time", bty = "l")
lines(ma.centered, lwd = 2, col = "red")
legend("topright", legend = c("Original Data", "Centered Moving Average"), lty = 1, col = c("black", "red"), cex = 0.65)

#data transformation 

Cpox2 <- read.csv("Chickenpox_copy.csv")
Cpox2$Date <- as.POSIXct(Cpox2$Date, format = "%d/%m/%Y")
#check data type
data_types2 <- sapply(Cpox2, class) #shows data type for each variable
print(data_types2)

Bud2.ts <- ts(Cpox2$BUDAPEST, start = c(2005, 1),end = c(2014, 52), freq = 52)

#plot the time series for budapest
plot(Cpox2$Date, Cpox2$BUDAPEST, type = "l", xlab = "Date", ylab = "Cases", main = "Time Series plot of Chickenpox Cases in Budapest")

##################################################################################################

#######baseline models#####
Budna.ts <- ts(Cpox2$BUDAPEST, start = c(2005, 1),end = c(2014, 52), freq = 52)

nValidna <- 52
(nTrainna <- length(Budna.ts) - nValidna)

train.tsna <- window(Budna.ts, start = c(2005, 1), end = c(2005, nTrainna))  
valid.tsna <- window(Budna.ts, start = c(2005, nTrainna + 1), end = c(2005, nTrainna + nValidna))
library(ggplot2)
library(forecast)

# Set training data from 2005 to 2014
# Plot some forecasts
autoplot(train.tsna) +
  autolayer(meanf(train.tsna, h=nValidna),
            series="Mean", PI=FALSE) +
  autolayer(naive(train.tsna, h=nValidna),
            series="Naïve", PI=FALSE) +
  autolayer(snaive(train.tsna, h=nValidna),
            series="Seasonal naïve", PI=FALSE) +
  ggtitle("Baseline Forecasts for Weekly Cases in Budapest") +
  xlab("Year") + ylab("Weekly Cases") +
  guides(colour=guide_legend(title="Forecast"))


# Generate the forecasts
mean_forecast <- meanf(train.tsna, h = nValidna)
naive_forecast <- naive(train.tsna, h = nValidna)
snaive_forecast <- snaive(train.tsna, h = nValidna)

# Calculate accuracy measures
accuracy(mean_forecast, valid.tsna)
accuracy(naive_forecast, valid.tsna)
accuracy(snaive_forecast, valid.tsna)

#Calculate residuals
snaive_residuals <- valid.tsna - snaive_forecast$mean


# Plot ACF of residuals
acf(coredata(snaive_residuals), main = "ACF of Seasonal Naive Residuals")


######################Build Predictive models####################################################

############Holt-Winters Model############

#build time series object

Bud2.ts <- ts(Cpox2$BUDAPEST, start = c(2005, 1),end = c(2014, 52), freq = 52)

#Split data into train and valid set 

nValid <- 52
(nTrain <- length(Bud2.ts) - nValid)

train.tsp <- window(Bud2.ts, start = c(2005, 1), end = c(2005, nTrain))  
valid.tsp <- window(Bud2.ts, start = c(2005, nTrain + 1), end = c(2005, nTrain + nValid))

#Fit an additive holt winter model 

fit <- HoltWinters(train.tsp, seasonal = "additive")
forecast <- forecast(fit, h = nValid)
accuracy(forecast, valid.tsp)

#Tune holt-winter model through grid search

# Define range of parameter values to search
alpha_values <- seq(0.1, 0.9, by = 0.1)
beta_values <- seq(0.1, 0.9, by = 0.1)
gamma_values <- seq(0.1, 0.9, by = 0.1)

# Specify seasonal periods (weekly data with annual and biyearly seasonality)
seasonality_values <- c(52, 104)  

# Initialize variables to store best parameters with minimum RMSE
best_alpha <- NULL
best_beta <- NULL
best_gamma <- NULL
best_seasonality <- NULL
min_rmse <- Inf

# Perform grid search
for (alpha in alpha_values) {
  for (beta in beta_values) {
    for (gamma in gamma_values) {
      for (seasonality in seasonality_values) {
        # Fit Holt-Winters model with current parameter values
        fit <- HoltWinters(train.tsp, alpha = alpha, beta = beta, gamma = gamma,
                           seasonal = "additive")
        
        # Generate forecasts
        hw_forecast <- forecast(fit, h = length(valid.tsp))
        
        # Calculate RMSE
        rmse <- sqrt(mean((hw_forecast$mean - valid.tsp)^2))
        
        # Update best parameters if RMSE is minimized
        if (rmse < min_rmse) {
          min_rmse <- rmse
          best_alpha <- alpha
          best_beta <- beta
          best_gamma <- gamma
          best_seasonality <- seasonality
        }
      }
    }
  }
}

# Print best parameters and minimum RMSE
cat("Best alpha:", best_alpha, "\n")
cat("Best beta:", best_beta, "\n")
cat("Best gamma:", best_gamma, "\n")
cat("Best seasonality:", best_seasonality, "\n")
cat("Minimum RMSE:", min_rmse, "\n")

#build model with best parameters 

hw_model <- HoltWinters(train.tsp, seasonal = "additive", alpha = 0.2, beta = 0.1, gamma = 0.4)

forecast <- forecast(hw_model, h = nValid)
accuracy(forecast, valid.tsp)
checkresiduals(hw_model)

#Visualization using Holt-Winters Model

Budtrainp.ts <- ts(Cpox2$BUDAPEST[1:nTrain], frequency = 1)
Budvalidp.ts <- ts(Cpox2$BUDAPEST[(nTrain + 1):(nTrain + nValid)], frequency = 1)

# Plot the time series without x-axis ticks
plot(Budtrainp.ts, xaxt = "n", xlim = c(1, nTrain + nValid), ylab = "Weekly Cases", main = "Forecast using Tuned Holt-Winters approach")

# Define the positions and labels for the x-axis
axis_pos <- seq(1, nTrain + nValid, length.out = 6)
axis_dates <- seq(2005, 2015, by = 2)  # Adjusted to start from 2006
axis_labels <- axis_dates

# Adjust axis_pos and axis_labels to have 5 positions and labels
axis_pos <- axis_pos[-6]  # Remove the last position
axis_labels <- axis_labels[-6]  # Remove the last label

# Add the x-axis labels
axis(1, at = axis_pos, labels = axis_labels)

# Specify the line details 

lines(x = seq(1, 468), forecast$fitted, lwd = 2, col = "blue")
lines(x = seq(469, 520), forecast$mean, lwd = 2, col = "green", lty = 2)
lines(x = seq(469, 520), valid.tsp)

#############ARIMA MODEL#########

full.ts <- ts(Cpox2$BUDAPEST, start = c(2005, 1),end = c(2014, 52), freq = 52)

#Take double difference of the time series object

diff<-diff(diff(full.ts,
          lag = 52), lag = 1)

par(mfrow = c(1, 2))
acf(coredata(diff))
pacf(coredata(diff))


nValid2 <- 52
(nTrain2 <- length(diff) - nValid2)

train.tsa <- window(diff, start = c(2006, 2), end = c(2006, nTrain2))  
valid.tsa <- window(diff, start = c(2006, nTrain2 + 1), end = c(2006, nTrain2 + nValid2))

######use Auto Arima function

AutoArimaModel <- auto.arima(train.tsa)

pred <- forecast(AutoArimaModel, h = nValid2, level = 0)
accuracy(pred, valid.tsa)
checkresiduals(AutoArimaModel)

#plot the model

plot(pred, ylab = "Weekly Cases", xlab= 'Time')
lines(valid.tsa, col ='yellow' )

legend("topright", legend = c("ARIMA forecasts", "Validation set"), 
       col = c("black", "yellow"), lty = 1,  cex = 0.65)



###########Tune the ARIMA Model ##########

arima1=arima(train.tsa, order=c(3,0,1), seasonal = list(order = c(0, 0, 0)))
arima2=arima(train.tsa, order=c(3,0,1), seasonal = list(order = c(0, 0, 1)))
arima3=arima(train.tsa, order=c(3,0,1), seasonal = list(order = c(0, 0, 2)))
arima4=arima(train.tsa, order=c(3,0,1), seasonal = list(order = c(1, 0, 1)))

#arima 1

arima1.pred<- forecast(arima1, h = nValid2, level = 0)
accuracy(arima1.pred, valid.tsa)
checkresiduals(arima1)

plot(arima1.pred, ylab = "Weekly Cases", xlab= 'Time')
lines(valid.tsa, col ='yellow' )

legend("topright", legend = c("ARIMA forecasts", "Validation set"), 
       col = c("black", "yellow"), lty = 1,  cex = 0.65)

#arima 2
arima2.pred<- forecast(arima2, h = nValid2, level = 0)
accuracy(arima2.pred, valid.tsa)
checkresiduals(arima2)

plot(arima2.pred, ylab = "Weekly Cases", xlab= 'Time')
lines(valid.tsa, col ='yellow' )

legend("topright", legend = c("ARIMA forecasts", "Validation set"), 
       col = c("black", "yellow"), lty = 1,  cex = 0.65)

#arima 3
arima3.pred<- forecast(arima3, h = nValid2, level = 0)
accuracy(arima3.pred, valid.tsa)
checkresiduals(arima3)

plot(arima3.pred, ylab = "Weekly Cases", xlab= 'Time')
lines(valid.tsa, col ='yellow' )

legend("topright", legend = c("ARIMA forecasts", "Validation set"), 
       col = c("black", "yellow"), lty = 1,  cex = 0.65)

#arima 4
arima4.pred<- forecast(arima4, h = nValid2, level = 0)
accuracy(arima4.pred, valid.tsa)
checkresiduals(arima4)

plot(arima4.pred, ylab = "Weekly Cases", xlab= 'Time')
lines(valid.tsa, col ='yellow' )

legend("topright", legend = c("ARIMA forecasts", "Validation set"), 
       col = c("black", "yellow"), lty = 1,  cex = 0.65)


######TBATS and STL+ Ets Models#####

Bud3.ts <- msts(Cpox2$BUDAPEST, seasonal.periods = c(52, 104))

nValid3 <- 52
(nTrain3 <- length(Bud3.ts) - nValid3)

train.msts <- window(Bud3.ts, start = c(1, 1), end = c(1, nTrain3))  
valid.msts<- window(Bud3.ts, start = c(1, nTrain3 + 1), end = c(1, nTrain3 + nValid3))

#Build Tbats model 

train.tbats <- tbats(train.msts)
train.tbats.pred <- forecast(train.tbats, h = nValid3)

#build Stl+ Ets model 

train.stlm <- stlm(train.msts, s.window = "periodic", method = "ets")
train.stlm.pred <- forecast(train.stlm, h = nValid3)
par(mfrow = c(1, 2))
plot(train.tbats.pred, xlab = "Year", ylab = "Weekly cases",
     main = "TBATS")
plot(train.stlm.pred, xlab = "Year", ylab = "Weekly cases",
     main = "STL + ETS")
accuracy(train.tbats.pred,valid.msts)
accuracy(train.stlm.pred,valid.msts)
checkresiduals(train.tbats)
checkresiduals(train.stlm)

################ARIMA 2 is the best model##############








