#------------------------------------------
# Projet de théorie des valeurs extrêmes
# Données financières

# Auteurs : Tom Brault & Alexandre Maghames 
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

# On va s'interesser aux prix de fermeture (Close)
serie <- data_df$MC.PA.Close
log_returns <- -diff(log(serie))
par(mfrow=c(1,1))
ts.plot(serie)

plot(seq_along(log_returns), log_returns, type = "l", col = "darkblue",
     xlab = "Time Index", ylab = "Negative Log Returns")
abline(h=0.03,col="red") # seuil choisi


par(mfrow=c(1,2))
acf(log_returns)
pacf(log_returns)
adf.test(log_returns) # faible p-val (on ne peux pas accepter H0 la serie est stationnaire)


# fitted by GPD 
par(mfrow=c(1,1))
# mrlplot: Empirical Mean Residual Life Plot
mrlplot(log_returns, main="",col="darkblue",ylim=c(0,0.06),xlim=c(-0.05,0.085))
abline(v=0.02,col="red")
abline(v=0.06,col="red")
grid()

par(mfrow=c(1,1))
# tcplot: Threshold Choice Plot
tcplot(log_returns, c(0.02,0.06),which=1)
abline(v=0.02,col="red")
abline(v=0.035,col="red")
grid()
tcplot(log_returns, c(0.02,0.06),which=2)
abline(v=0.02,col="red")
abline(v=0.03,col="red")
grid()

res <-numeric(10)
seuil_grid <- seq(from=0.02,to=0.03,length=10)

for(s_index in seq_along(seuil_grid)){
  seuil <- seuil_grid[s_index]
  res[s_index] <- sum(log_returns> seuil)
}
res
sum(log_returns> 0.03)

# GPD fit with a threshold of 0.03
# 250 (nombres de jours ouvrés dans une année)
fitted <- fpot(log_returns,0.03,npp=250)
# npp = 250: on raisonne en années, npp = 22: on raisonne en mois
# 250: nb de jours ouvrés dans une année, 22: nb de jours ouvrés dans un mois
# npp = 22 : nb de jours ouvrés moyen par mois
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
fitted_0 <- fpot(log_returns,0.045,npp=250,shape=0)
fitted_0
confint(fitted_0)
par(mfrow=c(2,2))
plot(fitted_0)

test_dev <- anova(fitted,fitted_0)
# on ne peut pas rejetter l'hypothèse nulle. 
# xi = 0 



S_stat <- test_dev$Deviance[2] - test_dev$Deviance[1]
qchisq(p=0.95, df=1) # loi à H_0
S_stat > qchisq(p=0.95, df=1)
# False : on ne peut pas rejeter H_0, i.e. on prend shape = 0







#######################################################################
# fitted by GEV 
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



#-------------------------------------------------------------------------------
# Comparison of GEV and GPD methods
#-------------------------------------------------------------------------------

# Represent of exceedances depending on a certain threshold
threshold = 0.093
colors_exceedances = ifelse(df$log_ret > threshold, "red", "darkgray")

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

plot(df$log_ret,type="l")

which(colors_exceedances=="red")

