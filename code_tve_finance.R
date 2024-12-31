#------------------------------------------
# Extreme Value Theory Project
# Financial Data

# Authors: Tom Brault & Alexandre Maghames
#------------------------------------------

# Load necessary libraries
library(quantmod)
library(evd)
library(tseries)
# Clear the environment
rm(list=ls())

# define seed 
set.seed(2002)

# Define the ticker symbol
start_date <- "2007-01-01"
end_date <- "2023-12-31"

# import data 
ticker <- "MC.PA"
getSymbols(ticker, src = "yahoo",from=start_date,to=end_date)

par(mfrow=c(3,2))
ts.plot(MC.PA$MC.PA.Open, main = "Open Prices", col = "blue")
ts.plot(MC.PA$MC.PA.High, main = "High Prices", col = "green")
ts.plot(MC.PA$MC.PA.Low, main = "Low Prices", col = "red")
ts.plot(MC.PA$MC.PA.Close, main = "Close Prices", col = "purple")
ts.plot(MC.PA$MC.PA.Volume, main = "Adjusted Close Prices", col = "orange")
ts.plot(MC.PA$MC.PA.Adjusted, main = "Volume", col = "brown")

par(mfrow=c(1,1))
plot(MC.PA$MC.PA.Close, ylab="Closing Prices",xlab="Time", col = "darkblue",main="")

data_df <- as.data.frame(MC.PA)

# We will focus on closing prices (Close)
serie <- data_df$MC.PA.Close
log_returns <- -diff(log(serie))
par(mfrow=c(1,1))
ts.plot(serie)

plot(seq_along(log_returns), log_returns, type = "l", col = "darkblue",
     xlab = "Time Index", ylab = "Negative Log Returns")
abline(h=0.03,col="darkred") ## Example of chosen threshold 

col_sup_seuil = ifelse(log_returns > 0.03, "darkred", "darkgray")

for (i in 1:(length(log_returns) - 1)) {
  lines(c(i, i + 1), 
        log_returns[i:(i + 1)], 
        col = col_sup_seuil[i])
}
sum(col_sup_seuil=="darkred") # 191 values above the threshold


par(mfrow=c(1,2))
acf(log_returns)
pacf(log_returns)
adf.test(log_returns) # Low p-value (we cannot accept H0, the series is stationary)


#######################################################################
# fitted by GPD 
#######################################################################
par(mfrow=c(1,1))
# mrlplot: Empirical Mean Residual Life Plot
mrlplot(log_returns, main="",col="darkblue",ylim=c(0,0.06),xlim=c(-0.05,0.085))
abline(v=0.02,col="darkred")
abline(v=0.06,col="darkred")
grid()

par(mfrow=c(1,1))
# tcplot: Threshold Choice Plot
tcplot(log_returns, c(0.02,0.06),which=1)
abline(v=0.02,col="darkred")
abline(v=0.035,col="darkred")
grid()
tcplot(log_returns, c(0.02,0.06),which=2)
abline(v=0.02,col="darkred")
abline(v=0.03,col="darkred")
grid()

res <-numeric(10)
seuil_grid <- seq(from=0.02,to=0.03,length=10)

for(s_index in seq_along(seuil_grid)){
  seuil <- seuil_grid[s_index]
  res[s_index] <- sum(log_returns> seuil)
}
res # Number of observed values greater than each chosen threshold in seuil_grid
sum(log_returns> 0.03)

npp <- 250

# GPD fit with a threshold of 0.03
# 250 (number of business days in a year)
fitted <- fpot(log_returns,0.03,npp=npp)
# npp = 250: reasoning in years, npp = 22: reasoning in months
# 250: number of business days in a year, 22: average number of business days in a month
fitted
confint(fitted)

par(mfrow=c(2,2))
plot(fitted)

par(mfrow=c(1,1))
plot(fitted,which=1)
plot(fitted,which=4,xlab="Return Period (in years)")

prof <- profile(fitted)
prof
par(mfrow=c(1,2))
plot(prof)
confint(prof)

# with shape = 0
fitted_0 <- fpot(log_returns,0.03,npp=npp,shape=0)
fitted_0
confint(fitted_0)
par(mfrow=c(2,2))
plot(fitted_0)

test_dev <- anova(fitted,fitted_0)
#we cannot reject the null hypothesis
# xi = 0 

S_stat <- test_dev$Deviance[2] - test_dev$Deviance[1]
qchisq(p=0.95, df=1) # loi à H_0
S_stat > qchisq(p = 0.95,df = 1)
# We get the test statistic. We cannot reject H_0, i.e. we take shape = 0
p_val <- 1 - pchisq(S_stat, 1) # We get the p-value


par(mfrow=c(1,1))
rl(fitted_0)

# Value of loss corresponding to a period of return of 20 years (=250*20=5000 working days approximately).
inverse_function_gpd = function(x,tau){
  return(-tau * log(1-x))}

return_level_gpd = function(T,tau){
  u = 0.03
  lambda = length(fitted_0$exceedances)/length(log_returns) # Proportion of exceedances of the threshold (denoted p(u) in the lectures)
  return(u+inverse_function_gpd(1-1/(T*lambda),tau))
}

# T: return period (in business days)
return_level_gpd(T = 5000, tau = as.numeric(fitted_0$estimate["scale"]))
# Return level corresponding to a return period of 20 years: 0.102
# with the GEV, we found a return level of 0.129 approximately

# Value of loss corresponding to a return period of 10 years (=250*10=2500 working days approximately).
return_level_gpd(T = 2500, tau = as.numeric(fitted_0$estimate["scale"]))
# Return level corresponding to a return period of 10 years: 0.093
# with the GEV, we found a return level of 0.111 approximately

# Indices of exceedances of the return level 0.093
for (i in 1:length(log_returns)) {
  if (log_returns[i] > 0.102) {
    cat("Index:", i, "Value:", log_returns[i], "\n")
  }
}
# 3 exceedances of this threshold in less than 20 years with a GPD. Note that 
# these exceedances are fairly close together (see the graph at the bottom of the code)

# Period of return corresponding to a certain return level
return_period_gpd <- function(y_p,tau){
  # y_p: return, tau: modified scale parameter
  u = 0.03 # u: threshold (fixed in the model)
  lambda = length(fitted_0$exceedances)/length(log_returns) # Proportion of exceedances of the threshold (denoted p(u) in the lectures)
  return(1/lambda * exp((y_p-u)/tau))
}

# Period of return corresponding to a return level of 0.111 for the GPD
# Value corresponding to a return period of 10 years in the case of a GEV
return_period_gpd(y_p = 0.111, tau = as.numeric(fitted_0$estimate["scale"]))
# The value is given in days here: 38241.62
# Value in years:
return_period_gpd(y_p = 0.111, tau = as.numeric(fitted_0$estimate["scale"]))/250
# Return level: 1.97 year (almost 2 years)



#######################################################################
# fitted by GEV 
#######################################################################

data <- get(ticker)
close_prices <- data[, "MC.PA.Close"]
df <- data.frame(Date = index(close_prices), Close = coredata(close_prices))
df$log_ret <- -c(0,diff(log(df$MC.PA.Close)))
df$Year <- format(df$Date, "%Y")
df$Year_Month <- format(df$Date,"%Y-%m")

to_do_stat_desc <- as.data.frame(table(df$Year_Month))
min(to_do_stat_desc$Freq)
to_do_stat_desc[which(to_do_stat_desc$Freq==max(to_do_stat_desc$Freq)),]
sqrt(var(to_do_stat_desc$Freq))
mean(to_do_stat_desc$Freq)

max_year <- (aggregate(log_ret~Year,data = df,FUN = max))       
max_year_month <- (aggregate(log_ret~Year_Month,data = df,FUN = max))       
par(mfrow=c(1,1))
plot(max_year$log_ret,type="l")
plot(max_year_month$log_ret,type="l")

fitted_gev_year <- fgev(max_year$log_ret)
fitted_gev_year
par(mfrow=c(2,2))
plot(fitted_gev_year)
confint(fitted_gev_year)
par(mfrow=c(1,3))
prof_gev1 <- profile(fitted_gev_year)
plot(prof_gev1)
confint(prof_gev1)

fitted_gev_year_month <- fgev(max_year_month$log_ret)
fitted_gev_year_month
par(mfrow=c(2,2))
plot(fitted_gev_year_month)

par(mfrow=c(1,1))
plot(fitted_gev_year_month,which = 4,xlab="Return Period (in months)")
confint(fitted_gev_year_month)
par(mfrow=c(1,3))
prof_gev2 <- profile(fitted_gev_year_month) 
plot(prof_gev2)
confint(prof_gev2)

fitted_gev_year_month$estimate
qgev(1-1/20,0.02730844,0.01145445,0.16412963)
qgev(1-1/120,0.02730844,0.01145445,0.16412963)
qgev(1-1/240,0.02730844,0.01145445,0.16412963)


#-------------------------------------------------------------------------------
# Comparison of GEV and GPD methods
#-------------------------------------------------------------------------------

# Represent of exceedances depending on a certain threshold
threshold = 0.093
colors_exceedances = ifelse(df$log_ret > threshold, "red", "darkgray")

par(mfrow=c(1,1))
# Create the plot without lines first
plot(df$log_ret, 
     type = "n", # Create an empty plot
     main = "Exceedances above the value at risk corresponding to a return period of 10 years (GPD)",
     xlab = "Time Index", 
     ylab = "Log Return")

# Add lines segment by segment with the respective colors
for (i in 1:(length(df$log_ret) - 1)) {
  lines(c(i, i + 1), 
        df$log_ret[i:(i + 1)], 
        col = colors_exceedances[i])
}

abline(h = threshold, col = "darkred", lty = 2)
text(x = length(df$log_ret) * 0.9,  # Adjust x position to place near the end
     y = threshold + 0.01,          # Adjust y position slightly above the line
     labels = paste0("Threshold: ", threshold),
     col = "darkred", cex = 0.8)    # cex adjusts text size

which(colors_exceedances=="red")
