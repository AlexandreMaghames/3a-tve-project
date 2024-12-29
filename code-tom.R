#------------------------------------------
# Projet de théorie des valeurs extrêmes
# Données financières
#------------------------------------------

rm(list = ls())

### Packages à utiliser

library(dplyr)
library(tidyr)
library(forcats)
library(FactoMineR)
library(ggplot2)
library(factoextra)
library(class)
library(quantmod) # package pour CA d'entreprises

set.seed(2002)

# Test avec Yahoo
#ticker = "AAPL" # code - entreprise
#getSymbols(ticker, src = "yahoo", from = "2020-01-01", to = "2023-01-01")

# Set date range
start_date <- "2007-01-01"
end_date <- "2023-12-31"

code_entreprise = "MC.PA"
# LVMH

# Fetch data for each ticker
getSymbols(code_entreprise, src = "yahoo", from = start_date, to = end_date)
# données disponibles sur Yahoo Finance

par(mfrow = c(2,2))

lvmh_ts <- ts(MC.PA)
plot(lvmh_ts)

# Convert to data frame data
lvmh_data <- as.data.frame(MC.PA)

cor(lvmh_data)

# On standardise les données (si nécessaire pour la suite)
lvmh_scaled <- scale(lvmh_data)

# On choisit de se focaliser sur la variable associée aux prix de fermeture/cloture
lvmh_closing_price <- lvmh_data["MC.PA.Close"]

lvmh_closing_price_ts <- ts(lvmh_closing_price)

diff_log_lvmh_closing_price <- -diff(log(lvmh_closing_price_ts))

par(mfrow = c(1,1))

plot(diff_log_lvmh_closing_price, main = "Negative log Returns of Closing Prices", 
     xlab = "Time Index", ylab = "Log Return", col = "darkblue")

par(mfrow = c(1,2))
acf(diff_log_lvmh_closing_price)
pacf(diff_log_lvmh_closing_price)

#-------------------------------------------------------------------------------
# GPD : Dépassement de seuil
#-------------------------------------------------------------------------------

library(evd)

# mrlplot: Empirical Mean Residual Life Plot
par(mfrow = c(1,1))
mrlplot(diff_log_lvmh_closing_price, 
        main = "Mean Residual Life Plot",
        col = "darkblue")
abline(v = 0.02, col="darkred")
abline(v = 0.06, col = "darkred")
# seuil entre 0.03 et 0.06 d'après le graphique

# tcplot: Threshold Choice Plot
par(mfrow = c(1,1))
tcplot(diff_log_lvmh_closing_price, tlim = c(0.02,0.06), which=1) # tlim: interval of thresholds
grid()
abline(v = 0.02, col="darkred")
abline(v = 0.035, col = "darkred")

par(mfrow = c(1,1))
tcplot(diff_log_lvmh_closing_price, tlim = c(0.02,0.06), which=2) # tlim: interval of thresholds
grid()
abline(v = 0.02, col="darkred")
abline(v = 0.03, col = "darkred")
# seuil choisi : 0.03

# Number of exceedances depending on the threshold value
threshold_values = seq(0.02,0.03,length.out=11)
number_exceedances = numeric(length(threshold_values))

for (i in seq_along(threshold_values)) {
  number_exceedances[i] = sum(diff_log_lvmh_closing_price > threshold_values[i])
}

n_exceedances_threshold_values = data.frame(threshold_values,number_exceedances)
n_exceedances_threshold_values
# En particulier : 191 dépassements pour le seuil u=0.03

# GPD fit with a threshold of 0.03
lvmh_gpd_fit = fpot(diff_log_lvmh_closing_price, 0.03, npp = 250)
# npp = 250: on raisonne en années, npp = 22: on raisonne en mois
# 250: nb de jours ouvrés dans une année, 22: nb de jours ouvrés dans un mois
# npp = 22 : nb de jours ouvrés moyen par mois
par(mfrow = c(2,2))
plot(lvmh_gpd_fit)

par(mfrow = c(1,1))
plot(lvmh_gpd_fit, which = 1) # Probability Plot
plot(lvmh_gpd_fit, which = 4) # Return Level Plot

lvmh_gpd_fit

confint(lvmh_gpd_fit)
# modified_scale != 0 et shape = 0

# with shape = 0:
lvmh_gpd_fit_shape0 = fpot(diff_log_lvmh_closing_price, 0.03, npp = 250, shape=0)

par(mfrow = c(2,2))
plot(lvmh_gpd_fit_shape0)

par(mfrow = c(1,1))
plot(lvmh_gpd_fit_shape0, which = 1) # Probability Plot
plot(lvmh_gpd_fit_shape0, which = 4) # Return Level Plot

# With larger scale for the return period
plot(lvmh_gpd_fit_shape0, which = 4, xlim = c(0.2,110), xaxt = "n") # Return Level Plot
axis(1, at = c(0.2,0.5,1,2,5,10,20,50,100,200), labels = c("0.2", "0.5", "1", "2", "5", "10", "20", "50", "100", "200"), las = 1)

confint(lvmh_gpd_fit_shape0)
# paramètre scale semble significatif, càd différent de 0

# Test du rapport de vraisemblance (déviance)
test_deviance_shape0 = anova(lvmh_gpd_fit, lvmh_gpd_fit_shape0)
test_deviance_shape0
# On ne peut pas rejeter l'hypothèse nulle selon laquelle shape = 0.
# Donc on pose shape = 0

# Autre méthode: réécriture du test de déviance
stat_test_deviance = test_deviance_shape0$Deviance[2] - test_deviance_shape0$Deviance[1]
qchisq(p=0.95, df=1) # loi à H_0
stat_test_deviance > qchisq(p=0.95, df=1)
# False : on ne peut pas rejeter H_0, i.e. on prend shape = 0

# Interpretation with the return level plot
rl(lvmh_gpd_fit_shape0) # Return Level Plot

# Value of loss corresponding to a period of return of 20 years (=250*20=5000 working days approximately).
inverse_function_gpd_shape0 = function(x,tau){
  return(-tau * log(1-x))}

return_level_gpd_shape0 = function(T,tau){
  u = 0.03
  lambda = length(lvmh_gpd_fit_shape0$exceedances)/length(diff_log_lvmh_closing_price) # proportion de dépassements du seuil (p(u) dans le cours)
  return(u+inverse_function_gpd_shape0(1-1/(T*lambda),tau))
}
# T: période de retour (en jours ouvrés)
return_level_gpd_shape0(T = 5000, tau = as.numeric(lvmh_gpd_fit_shape0$estimate["scale"]))
# Return level corresponding to a period of return of 20 years: 0.102
# with the GEV, we found a return level of 0.129 approximately

# Value of loss corresponding to a period of return of 10 years (=250*10=2500 working days approximately).
return_level_gpd_shape0(T = 2500, tau = as.numeric(lvmh_gpd_fit_shape0$estimate["scale"]))
# Return level corresponding to a period of return of 10 years: 0.093
# with the GEV, we found a return level of 0.111 approximately

# Indices des dépassements du niveau de retour 0.093
for (i in 1:length(diff_log_lvmh_closing_price)) {
  if (diff_log_lvmh_closing_price[i] > 0.102) {
    cat("Index:", i, "Value:", diff_log_lvmh_closing_price[i], "\n")
  }
}
# 3 dépassements de ce seuil en moins de 20 ans avec une GPD. Notons que 
# ces dépassements sont assez resserrés entre eux (voir graphique tout en bas du code)

# Period of return corresponding to a certain return level

return_period_gpd_shape0 = function(y_p,tau){
  # y_p: return, tau: modified scale parameter
  u = 0.03 # u: threshold (fixed in the model)
  lambda = length(lvmh_gpd_fit_shape0$exceedances)/length(diff_log_lvmh_closing_price) # proportion de dépassements du seuil (p(u) dans le cours)
  return(1/lambda * exp((y_p-u)/tau))
}

# Period of return corresponding to a return level of 0.111 for the GPD
# Value corresponding to a return period of 10 years in the case of a GEV
return_period_gpd_shape0(y_p = 0.111, tau = as.numeric(lvmh_gpd_fit_shape0$estimate["scale"]))
# The value is given in days here: 38241.62
# Value in years:
return_period_gpd_shape0(y_p = 0.111, tau = as.numeric(lvmh_gpd_fit_shape0$estimate["scale"]))/250
# Return level: 1.97 year (almost 2 years)

#-------------------------------------------------------------------------------
# GEV
#-------------------------------------------------------------------------------

# Dataframe avec les maxima par mois
# Un mois compte environ 22 jours ouvrés

n_log_returns_date = data.frame(date=index(lvmh_closing_price), close = coredata(lvmh_closing_price))
n_log_returns_date$n_log_returns = c(0,diff(log(n_log_returns_date$MC.PA.Close)))

head(n_log_returns_date)

max_closing_prices = n_log_returns_date %>%
  mutate(subsections = (row_number()-1) %/% 22) %>% # %/% : dividende ou quotient
  group_by(subsections) %>%
  slice_max(n_log_returns, with_ties = FALSE)

lvmh_gev_fit = fgev(max_closing_prices$n_log_returns)
lvmh_gev_fit
plot(lvmh_gev_fit)

confint(lvmh_gev_fit)
lvmh_gev_fit

#-------------------------------------------------------------------------------
# Comparison of GEV and GPD methods
#-------------------------------------------------------------------------------

# Represent of exceedances depending on a certain threshold
threshold = 0.093
colors_exceedances = ifelse(diff_log_lvmh_closing_price > threshold, "red", "darkgray")

# Create the plot without lines first
plot(diff_log_lvmh_closing_price, 
     type = "n", # Create an empty plot
     main = "Exceedances above the value at risk corresponding to a return period of 10 years (GPD)",
     xlab = "Time Index", 
     ylab = "Log Return")

# Add lines segment by segment with the respective colors
for (i in 1:(length(diff_log_lvmh_closing_price) - 1)) {
  lines(c(i, i + 1), 
        diff_log_lvmh_closing_price[i:(i + 1)], 
        col = colors_exceedances[i])
}

abline(h = threshold, col = "darkred", lty = 2)

