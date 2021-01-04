
rm(list=ls())
setwd("~/Google Drive (nmilikic@nd.edu)/Spring 2020/ACMS 60855 Spatio-Temporal Statistics for Environmental Applications/Final")

library(CARBayes)
library(spdep)
library(grid)

set.seed(1)

load("lipscotland.RData")

######################
# Exploratory Analysis
######################

# Getting aspect ratio for plotting (converting degrees latitude/longitude to distance to correspond to how Scotland would look on a map)
# Distance of 1 degree latitude and average distance of 1 degree longitude in Scotland
long = cos(mean(lipscotland$latitude) * pi / 180)
aspect.ratio = 1 / long

# Observed counts
summary(lipscotland$observed)
lipscotland$observed.plot = ifelse(lipscotland$observed > 20, 21, lipscotland$observed)
# Cap observations at 20 for plotting; anything over 20 is treated as just > 20
breaks.obs = c(0, 4, 8, 12, 16, 20, 24)
spplot(lipscotland, "observed.plot", at = breaks.obs, aspect = aspect.ratio)
grid.text("Observed Male Lip Cancer Cases", x = unit(0.82, "npc"), y = unit(0.50, "npc"), rot = 90)

# Expected counts
summary(lipscotland$expected)
lipscotland$expected.plot = ifelse(lipscotland$expected > 20, 21, lipscotland$expected)
# Same as before; some information is lost as some of the expected counts are much higher
# than 20 (up to 89)
breaks.exp = c(0, 4, 8, 12, 16, 20, 24)
spplot(lipscotland, "expected.plot", at = breaks.exp, aspect = aspect.ratio)
grid.text("Expected Male Lip Cancer Cases", x = unit(0.82, "npc"), y = unit(0.50, "npc"), rot = 90)

# Percentage employed in agriculture, fishing, or forestry
summary(lipscotland$pcaff)
breaks.pcaff = c(0, 4, 8, 12, 16, 20, 24)
spplot(lipscotland, "pcaff", at = breaks.pcaff, aspect = aspect.ratio)
grid.text("Percentage Employed in Agriculture, Fishing, or Forestry", x = unit(0.82, "npc"), y = unit(0.50, "npc"), rot = 90)

hist(lipscotland$observed, xlab = "", main = "Observed male lip cancer count")
hist(lipscotland$expected, breaks = seq(0, 90, 5), xlab = "", main = "Expected male lip cancer count")
hist(lipscotland$pcaff, xlab = "", main = "Percentage employed in agriculture, fishing, or\nforestry")

mean(lipscotland$observed)
var(lipscotland$observed)
# The mean and variance are not close to one another; this may be problematic in assuming a Poisson distribution

###################
# Non-spatial model
###################

# Is the effect of pcaff on region risk linear?
plot(lipscotland$pcaff, lipscotland$observed / lipscotland$expected, xlab = "pcaff", ylab = "Region Risk")
# Approximately; it's not obviously nonlinear

form = observed ~ offset(log(expected)) + pcaff
lip.nsa = glm(formula = form, family = poisson, data = lipscotland)
summary(lip.nsa)

mod.resid = lip.nsa$residuals
moran.mc(x = mod.resid, listw = W.list, nsim = 1000)
geary.mc(x = mod.resid, listw = W.list, nsim = 1000)

##################
# Spatial analysis
##################

lip.spatial = S.CARleroux(formula = form, data = lipscotland@data, family = "poisson", W = W, burnin = 20000, n.sample = 120000, thin = 10)
print(lip.spatial)

# Plot fitted values
lipscotland$fitted.car = fitted(lip.spatial)
lipscotland$fitted.car.plot = ifelse(lipscotland$fitted.car > 20, 21, lipscotland$fitted.car)
# Same as before, capping fitted values at 20 and treat any greater values as just > 20
spplot(lipscotland, "fitted.car.plot", at = breaks.obs, aspect = aspect.ratio)
grid.text("Fitted Male Lip Cancer Cases", x = unit(0.82, "npc"), y = unit(0.50, "npc"), rot = 90)

# Compare fitted values to true values
cbind(lipscotland$observed, lipscotland$fitted.car)

summarise.samples(lip.spatial$samples$beta, quantiles = c(0.5, 0.025, 0.975))
# Intercept
plot(lip.spatial$samples$beta[,1], ylab = "beta_0", xlab = "Sample")
# pcaff effect
plot(lip.spatial$samples$beta[,2], ylab = "beta_1", xlab = "Sample")
# tau2
plot(lip.spatial$samples$tau, ylab = "tau2", xlab = "Sample")
# rho
plot(lip.spatial$samples$rho, ylab = "rho", xlab = "Sample")



















