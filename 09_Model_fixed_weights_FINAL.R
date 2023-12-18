# Bayesian mixture model for estimating missing soil water content (SWC) during
# data gaps, applied to pseudo data of varying gap characteristics. This model 
# assumes time-invariant (fixed), depth-specific mixture weights. The missing 
# SWC values are modeled as a weighted combination of the linearly interpolated 
# values and the adjusted SOILWAT2 values.

## Some data preparation:
data{
  for(i in start:end){
    # Compute adjusted SOILWAT2 values given fixed values for coefficients (a):
    soilwat.pred[i] <- a[1,depth.id[i]] + a[2,depth.id[i]]*SWC_soilwat[i] + 
      a[3,depth.id[i]]*delta.SWC_soilwat[i] + 
      a[4,depth.id[i]]*SWC_soilwat[i]*delta.SWC_soilwat[i]
  }
}

## Start model specification:
model{
  for(i in start:end){
    # Likelihood of SWC pseudo-data; precision (variance) is allowed to vary
    # by indicator for missing data gap (gmiss) and soil depth ID:
    SWC_obs[i] ~ dnorm(mu[i], tau[gmiss[i],depth.id[i]])
    # For the mean (mu) model, if observation i is before/after gap (gmiss = 1), 
    # then use a simple, depth-specific mean (mu.obs); if during the missing 
    # data gap (gmiss = 2), then use the mixture model (w.part):
    mu[i] <- (gmiss[i]==1)*mu.obs[depth.id[i]] + (gmiss[i]==2)*w.part[i]
    w.part[i] <- w[i]*SWC_linear.interp[i] + (1-w[i])*soilwat.pred[i]
    
    # Assume constant mixture weight:
    w[i] <- w.const[depth.id[i]]
    
    # Replicated data for evaluating model fit:
    SWC_obs.rep[i] ~ dnorm(mu[i], tau[depth.id[i]])
    # Squared deviation:
    Sqdiff[i] <- pow(SWC_obs[i]-SWC_obs.rep[i],2)
  }
  # Sum of squared deviations:
  Dsum <- sum(Sqdiff[])
  
  for(v in 1:Nlayers){
    # Relatively non-informative, independent priors for weights:
    w.const[v] ~ dunif(0,1)
    for(k in 1:2){  
      # Relatively non-informative prior for precision in likelihood:
      tau[k,v] ~ dgamma(0.1, 0.1)
      # Compute standard deviation:
      sig[k,v] <- 1/sqrt(tau[k,v])
    }
  }
}