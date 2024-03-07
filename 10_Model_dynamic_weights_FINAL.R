# Bayesian mixture model for estimating missing soil water content (SWC) during
# data gaps, applied to pseudo data of varying gap characteristics. This model 
# allows the mixture weights to vary over time ("dynamic weights"), in response 
# to environmental covariates. The missing SWC values are modeled as a weighted 
# combination of the linearly interpolated values and the adjusted SOILWAT2 values.


## Some data preparation:
data{
  for(i in start:end){
    # Compute adjusted SOILWAT2 values given fixed values for coefficients (a):
    soilwat.pred[i] <- a[1,depth.id[i]] + a[2,depth.id[i]]*SWC_soilwat[i] + 
      a[3,depth.id[i]]*delta.SWC_soilwat[i] + 
      a[4,depth.id[i]]*SWC_soilwat[i]*delta.SWC_soilwat[i]
    
    # Some relevant covariates:
    # Days since start of gap:
    days.ssog[i] <- gapday[i] - 1
    # Cumulative precipitation since start of gap:
    cum.ppt[i] <- pptday[i]*(days.ssog[i]+2)
    
    # Create standardized covariates to use in mixture model;
    # only center temp on overall mean; scale all 3 by their sd:
    Z[1,i] <- cum.ppt[i]/sd(cum.ppt[])
    Z[2,i] <- (tempday[i] - mean(tempday[]))/sd(tempday[])
    Z[3,i] <- days.ssog[i]/sd(days.ssog[])
    # Cumulative precipitation x average daily temperature interaction:
    Z[4,i] <- Z[1,i]*Z[2,i]
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
    
    # Define logit for mixture weight; assume that weight given to simple linear 
    # interpolation vs adjusted SOILWAT2 can vary, depending potentially on
    # cumulative daily precipitation received since start of gap, average daily 
    # temperature since start of gap, time (days) since gap started, and the
    # interaction between cumulative precip and average daily temp:
    logit(w[i]) <- 
      b[1,depth.id[i]] + 
      b[2,depth.id[i]]*Z[1,i] + 
      b[3,depth.id[i]]*Z[2,i] +
      b[4,depth.id[i]]*Z[3,i] +
      b[5,depth.id[i]]*Z[4,i]
    
    # Replicated data for evaluating model fit:
    SWC_obs.rep[i] ~ dnorm(mu[i], tau[gmiss[i],depth.id[i]])
    # Squared deviation:
    Sqdiff[i] <- pow(SWC_obs[i]-SWC_obs.rep[i],2)
  }
  # Sum of squared deviations:
  Dsum <- sum(Sqdiff[])
  
  # Independent priors for coefficients in the logit model (b's):
  for(v in 1:Nlayers){ 
    mu.obs[v] ~ dnorm(0,0.0001)
    for(p in 1:5){
      # Relatively non-informative priors:
      b.temp[p,v] ~ dnorm(0, 0.0001)
      # If b.use = 0 then the associated coefficient is set to zero, otherwise,
      # estimate the coefficient (b.use = 1):
      b[p,v] <- b.temp[p,v]*b.use[p,v]
    }
    for(k in 1:2){  
      # Relatively non-informative prior for precision in likelihood:
      tau[k,v] ~ dgamma(0.1, 0.1)
      # Compute standard deviation:
      sig[k,v] <- 1/sqrt(tau[k,v])
    }
    # If fixed weight for layer v (i.e., all coefficients, except intercept
    # are set to 0), then look at w.b (the weight for layer v); w.b may also be 
    # interpreted as the mixture weight when cum.ppt = 0, at the mean daily temp, 
    # and at the beginning of a missing data gap.
    w.b[v] <- ilogit(b[1,v])
  }
}