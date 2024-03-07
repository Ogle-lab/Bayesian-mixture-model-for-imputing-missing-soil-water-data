# Script for fitting Bayesian dynamic mixture weight model (using pseudodata
# with simulated missing data gaps), to obtain estimates of the mixture weights
# and covariate effects in the mixture model; produces posterior samples of the 
# effects parameters, which are subsequently used to impute missing SWC in the 
# observed timeseries from each site.

# Script can be easily modified to implement the fixed weights model.

# This script is for the MPJ site; scripts for other sites are identical, just
# some data are updated for each site.

# Load necessary libraries and modules:
library(rjags)
library(mcmcplots)
load.module("dic")


### Load / get relevant data.
load("./AuxiliaryData/a_data.RData")

# IMPORTANT: These data are specific to each SITE!
load("./AuxiliaryData/mpj.RData")
load("./AuxiliaryData/inits_timesince.RData")
inits2 = inits.mpj


### General set-up for running models.

# Set variables for jags.model and coda.samples:
n.adapt = 5000
n.iter2 = 10000
n.iter.more = 2000
n.iter.w = 2000


# Set b.use to indicate which, if any, coefficients are set to zero (b.use = 0);
# IMPORTANT: Need to update b.use for each site.
b.use = matrix(data = c(1, 1, 1, 
                        1, 1, 1, 
                        1, 1, 1, 
                        0, 1, 1, 
                        1, 1, 1), 
               nrow = 5, ncol = nc, byrow = TRUE)

# Create data list for jags.model
# IMPORTANT: Need to update a for each site.
dat2 <- list(a = a.mpj, b.use = b.use,
             start = 1, end = length(gstart), Nlayers = nc, 
             SWC_obs = obs.swc,
             SWC_linear.interp = lin.int.swc, 
             SWC_soilwat = sim.swc,
             delta.SWC_soilwat = delta.sim.swc,
             gmiss=gmiss,
             depth.id = depth.id,
             pptday = pptday,
             gapday = gapday,
             tempday = tempday,
             ave.temp = ave.temp,
             sd.temp = sd.temp)

# Initialize jags model for dynamic weights model:
jm2 <- jags.model(file = "10_Model_dynamic_weights_FINAL.R", 
                  data = dat2, 
                  n.chains = 3, 
                  inits = inits2, n.adapt = n.adapt)


# Obtain coda (mcmc) samples and monitor quantities of interest:
coda2 <- coda.samples(jm2,variable.names=c("b","b.temp","sig","mu.obs","tau", 
                                           "deviance","w.b",
                                           "Dsum"),
                      n.iter=n.iter2)

# Run coda.samples again to monitor observation-level predicted quantities:
coda2.w <- coda.samples(jm2,variable.names=c("w"), n.iter=n.iter.w)
coda2.obs <- coda.samples(jm2,variable.names=c("SWC_obs.rep"),
                          n.iter=n.iter.more)

# Compute posterior statistics:
stats2 <- summary(coda2)
postout2 <- stats2$statistics[,1:2]
postout2 <- cbind(postout2,stats2$quantiles[,c(1,3,5)])
stats2.obs <- summary(coda2.obs)
postout2.obs <- stats2.obs$statistics[,1:2]
postout2.obs <- cbind(postout2.obs,stats2.obs$quantiles[,c(1,3,5)])
stats2.w <- summary(coda2.w)
postout2.w <- stats2.w$statistics[,1:2]
postout2.w <- cbind(postout2.w,stats2.w$quantiles[,c(1,3,5)])


# Save coda (mcmc samples) and posterior statistics:
# IMPORTANT: Rename file appropriate to site.
save(coda2, postout2, postout2.obs, postout2.w,
     file = "posterior_dynamicW_mpj.RData")