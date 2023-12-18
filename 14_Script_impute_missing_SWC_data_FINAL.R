# Script to impute missing soil water content (SWC) data for each site's
# SWC timeseries, based on results from fitting the Bayesian dynamic mixture
# model to the pseudo-data.

##### !! Select site to evaluate:
site.names = c("mpj","seg","ses","vcm","vcp","vcs","wjs")
site = "mpj"

# Load data created by running "13_Script_to_create_observed_data_FINAL.R"
# Site-specific observed SWC (with missing data gaps), SOILWAT2 data, and
# environmental data:
load(paste0("./AuxiliaryData/",site,"_obs_data.RData"))
# Fixed parameters for computing adjusted SOILWAT2 values:
load("./AuxiliaryData/a_data.RData")
a = get(paste0("a.",site))
# Covariate statistics (e.g., means and/or standard deviations) used for
# standardizing covariates in the Bayesian mixture model:
load(paste0("./AuxiliaryData/envstats_",site,".RData"))

# Load posterior samples prduced by running the Bayesian mixture model
# with the pseudo-data via "11_ScriptDynamicWeightModel_mpj_FINAL.R"; coda
# for MPJ are provided as an example.
load(paste0("./AuxiliaryData/posterior_dynamicW_",site,".RData"))
# Set thinning, c.thin, so that have about 6000-9000 mcmc samples to impute 
# missing SWC values; based on the coda samples we obtained, we used:
# mpj = 10 (2000/chain), seg = 10 (2000/chain), ses = 3 (3000/chain), 
# vcm = 10 (2500/chain), vcp = 10 (2000/chain), vcs = 5 (2000/chain), 
# wjs = 10 (2000/chain)
c.thin = c(10, 10, 3, 10, 10, 5, 10)

# If coda lists hasn't already been converted to matrix form, do so here.
# The example coda for MPJ has already been converted to matrix form, so no
# need to run this code:
mr = dim(coda2[[1]])[1]
mc = dim(coda2[[1]])[2]
coda2.temp = matrix(data = NA, nrow = 3*mr, ncol = mc)
coda2.temp[1:mr,] = coda2[[1]]
coda2.temp[(1+mr):(2*mr),] = coda2[[2]]
coda2.temp[(1+2*mr):(3*mr),] = coda2[[1]]
colnames(coda2.temp) = colnames(coda2[[1]])
coda2 = coda2.temp

# Total number of mcmc samples, across all chains:
n.its = dim(coda2)[1]

# Thin the coda samples (create coda1 as the thinned sample):
coda1 = coda2[seq(from=1, to=n.its, by = c.thin[site.names==site]),]
# Number of mcmc samples after thinning:
Nmcmc = dim(coda1)[1]

# Various variables needed to compute SWC via the mixture model:
Nlayers = nc
start.data = 1
end.data = dim(SWC_obs)[1]
SWC_obs.data = SWC_obs
linear.data = linear.interp
soilwat.data = soilwat.interp
delta.soilwat.data = delta.soilwat.data


############################################################
# Calculate predicted SWC for actual, observed timeseries:
############################################################

# Create array of b coefficients, to be used in the mixture model:
b = array(data=NA,dim=c(5,Nlayers,Nmcmc))
for(p in 1:5){ # parameters
  for(v in 1:Nlayers){ # depths or layers
    c.get = match(paste0("b[",p,",",v,"]"),colnames(coda1))
    b[p,v,] = as.vector(coda1[,c.get])
  }
}

# Get / create environmental data:
days.ssog = array(data = NA, dim = c(end.data,Nlayers))
cum.ppt = array(data = NA, dim = c(end.data,Nlayers)) 
Z = array(data = NA, dim = c(4,end.data,Nlayers))
# Days since start of gap:
days.ssog = gapday - 1
# Cumulative precip since start of gap:
cum.ppt = pptday*(days.ssog+2)

# Create standardized covariates to use in mixture model:
Z[1,,] = cum.ppt/sd.cum.ppt
Z[2,,] = (tempday - mean.tempday)/sd.tempday
Z[3,,] = days.ssog/sd.days.ssog
# cumulative ppt x average daily temp interaction:
Z[4,,] <- Z[1,,]*Z[2,,]

# Compute the adjusted SOILWAT2 values:
soilwat.pred.data = array(data = NA, dim = c(end.data, Nlayers))
for(v in 1:Nlayers){
  # Adjusted SOILWAT2 values:
  soilwat.pred.data[,v] <- a[1,v] + a[2,v]*soilwat.data[,v] + 
    a[3,v]*delta.soilwat.data[,v] + 
    a[4,v]*soilwat.data[,v]*delta.soilwat.data[,v]
}

# Compute (impute) the predicted SWC for each observation point, and for
# each MCMC sample of the b coefficients in the mixture model:
u.data = array(data = NA, dim = c(end.data, Nlayers, Nmcmc))
w.data = array(data = NA, dim = c(end.data, Nlayers, Nmcmc))
w.part = array(data = NA, dim = c(end.data, Nlayers, Nmcmc))
pred.SWC = array(data = NA, dim = c(end.data, Nlayers, Nmcmc))
for(m in 1:Nmcmc){
  for(v in 1:Nlayers){
    for(i in start.data:end.data){
      # Compute the logit-scale mixture weight:
      u.data[i,v,m] = b[1,v,m] + b[2,v,m]*Z[1,i,v] + b[3,v,m]*Z[2,i,v] + 
        b[4,v,m]*Z[3,i,v] + b[5,v,m]*Z[4,i,v]
      # Compute mixture weight (inverse-logit):
      w.data[i,v,m] = exp(u.data[i,v,m])/(1+exp(u.data[i,v,m]))
      
      # Predicted SWC:
      # If during the missing data gap (gmiss = 2), then use the mixture model
      # to impute missing values; otherwise, return the observed value:
      w.part[i,v,m] = w.data[i,v,m]*linear.data[i,v] + (1-w.data[i,v,m])*soilwat.pred.data[i,v]
      if(is.na(gmiss[i,v])==TRUE | gmiss[i,v]==1){
        # If not during a "true" missing data gap:
        pred.SWC[i,v,m] = SWC_obs.data[i,v]
      }else{
        # If during a missing data gap:
        pred.SWC[i,v,m] = w.part[i,v,m]
      }
    }
  }
}

# Summarize the imputed SWC values via posterior mean, median, standard
# deviation, and 2.5th and 97.5th percentiles:
mean.SWC = apply(pred.SWC,c(1,2),mean, na.rm = TRUE)
median.SWC = apply(pred.SWC,c(1,2),median, na.rm = TRUE)
sd.SWC = apply(pred.SWC,c(1,2),sd, na.rm = TRUE)
lower.SWC = apply(pred.SWC,c(1,2), quantile,probs = c(0.025), na.rm = TRUE)
upper.SWC = apply(pred.SWC,c(1,2), quantile,probs = c(0.975), na.rm = TRUE)
stats.SWC = cbind(mean.SWC, median.SWC, sd.SWC, lower.SWC, upper.SWC)
colnames(stats.SWC) = c(paste0("mean_",1:Nlayers),paste0("median_",1:Nlayers), 
                        paste0("sd_",1:Nlayers), paste0("2.5th_",1:Nlayers),
                        paste0("97.5th_",1:Nlayers))



# Get missing data code (1 = SWC was reported; 2 = SWC was missing and imputed)
gmiss.data = gmiss
colnames(gmiss.data) = paste0("code_",1:dim(gmiss)[2])

# Get date information from original data:
# Read in the full data for each site
fulldata.ses <- read.csv("01_US-Ses_SWC.csv")
fulldata.seg <- read.csv("02_US-Seg_SWC.csv")
fulldata.wjs <- read.csv("03_US-Wjs_SWC.csv")
fulldata.mpj <- read.csv("04_US-Mpj_SWC.csv")
fulldata.vcp <- read.csv("05_US-Vcp_SWC.csv")
fulldata.vcm <- read.csv("06_US-Vcm_SWC.csv")
fulldata.vcs <- read.csv("07_US-Vcs_SWC.csv")

dat = get(paste0("fulldata.",site))

dstamp = cbind(dat$Date,dat$Year,dat$DOY,dat$Month,dat$Day)
colnames(dstamp) = c("Date","Year","DOY","Month","Day")

# Combine all data into one array:
new.data = cbind(dstamp,gmiss.data,stats.SWC)

# Save data created:
# E.g., the data generated for MPJ are contained in "19_imputed_SWC_mpj.csv"
# See "15_README_imputed_SWC_data.txt" for a description of the columns 
# (variables) the resulting file, "imputed_SWC_SITE.csv"
write.csv(new.data,file=paste0("imputed_SWC_",site,".csv"))
