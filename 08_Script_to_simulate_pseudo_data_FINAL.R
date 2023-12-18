# Loads observed timeseries data for each site, then samples from observed data
# to create pseudo data representing a range of missing data gap characteristics.


## DON'T NEED TO RUN BLOCK OF CODE FROM LINE 8 TO LINE 26; CODE WAS
## ALREADY RUN TO PRODUCE "stratified_sampled_gap_lengths_USE.RData".
## Code simply provided for documentation purposes.

# For sampling gap lengths, draw for discrete uniform distributions with
# min (MinL) and max (MaxL) gap lengths (days) defined by:
MinL = c(1, 3, 5, 7, 9, 11, 14, 17, 20, 24, 29, 34, 41, 51, 71, 101)
MaxL = c(2, 4, 6, 8, 10, 13, 16, 19, 23, 28, 33, 40, 50, 70, 100, 300)
cbind(MinL,MaxL,MaxL-MinL)

# Simulate collection of gap lengths:
Ngaps = 150
gsim = c()
Nj = length(MinL)
for(i in 1:Ngaps){
  j = (i%%Nj)+1
  gsim[i] = rdunif(n=1,min=MinL[j],max=MaxL[j])
}
## DON'T NEED TO RUN ABOVE CODE.


## START HERE FOR SIMULATING PSEUDO DATA

# Load relevant libraries:
library(imputeTS)
library(extraDistr)

# Load gap lengths for 150 gaps (gsim = simulated gap lengths)
load("./AuxiliaryData/stratified_sampled_gap_lengths_USE.RData")

# Read in the full datasets for individual sites:
fulldata.ses <- read.csv("01_US-Ses_SWC.csv")
fulldata.seg <- read.csv("02_US-Seg_SWC.csv")
fulldata.wjs <- read.csv("03_US-Wjs_SWC.csv")
fulldata.mpj <- read.csv("04_US-Mpj_SWC.csv")
fulldata.vcp <- read.csv("05_US-Vcp_SWC.csv")
fulldata.vcm <- read.csv("06_US-Vcm_SWC.csv")
fulldata.vcs <- read.csv("07_US-Vcs_SWC.csv")

# Now, given above simulated gap lengths, simulate psuedo-data for each site.
# Sites: ses, seg, wjs, mpj, vcp, vcm, vcs

site.names = c("mpj","seg","ses","vcm","vcp","vcs","wjs")

# Loop through all sites:
for(s in 1:7){
  
  site = site.names[s]
  
  # Get original data for focal site:
  fulldata = get(paste0("fulldata.",site))
  
  # Load the indices for grabbing each site's observed SWC, soilwat SWC, and env data:
  source("./AuxiliaryData/site_data_indices.R")
  obs.id <- get(paste0("obs.",site)) # Column indices for observed SWC
  sim.id <- get(paste0("sim.",site)) # Column indices for SOILWAT2 SWC
  env.id <- get(paste0("env.",site)) # Column indices for environmental variables
  
  # Variable for matching observed soil depths to corresponding SOILWAT2 layers
  source("./AuxiliaryData/soilwat_match_depths.R")
  match = get(paste0("match.",site))
  
  # Potential environmental predictor variables:
  SWC_x <- fulldata[,env.id]
  
  # Pick temperature variable to use (e.g., max temp is picked here):
  Tuse.temp <- SWC_x$Input_AirTemp_max_C
  # If any missing temp data (rare), fill-in via linear interpolation:
  Tuse <- na_interpolation(Tuse.temp, option="linear")
  
  # Extract / get precipitation data:
  if(site=="vcs") SWC_x$P = SWC_x$Input_PPT_mm
  ppt = SWC_x$P
  
  # Defined / extract observed SWC data:
  SWC_obs <- fulldata[,obs.id]
  nr = dim(SWC_obs)[1] # number of days
  nc = dim(SWC_obs)[2] # number of depths
  
  # Set number of missing data gaps to simulate for each layer:
  Nrep = rep(Ngaps,times=nc)

  # Create imputed missing observed values based on linear interpolation of missing
  # values (NAs):
  linear.interp = array(data=NA, dim = dim(SWC_obs))
  for(v in 1:nc){
    linear.interp[,v] <- na_interpolation(SWC_obs[,v], option = "linear")
  }
  # Find the row or observation ID associated with the first and last
  # (non-missing) SWC in the timeseries for each depth:
  obs.start = rep(1,nc)
  obs.end = rep(1,nc)
  for(v in 1:nc){
    obs.start[v] <- min(which(SWC_obs[,v]>0))
    obs.end[v] <- max(which(SWC_obs[,v]>0))
  }
  
  # Get/extract SOILWAT2 simulated SWC:
  SWC_soilwatX <- fulldata[,sim.id]
  SWC_soilwat <- SWC_soilwatX[,1:2]
  for(v in 1:nc){
    # Average SOILWAT2 values of the 2 layers above/below the observed layer
    # (in come cases, will average over the same SOILWAT2 layer)
    SWC_soilwat[,v] <- (SWC_soilwatX[,match[v,1]] + SWC_soilwatX[,match[v,2]])/2
  }
  # Find the row or observation ID associated with the first and last simulated 
  # (non-missing) SWC from SOILWAT2
  sim.start = rep(1,nc)
  sim.end = rep(1,nc)
  for(v in 1:nc){
    sim.start[v] <- min(which(SWC_soilwat[,v]>0))
    sim.end[v] <- max(which(SWC_soilwat[,v]>0))
  }
  # SOILWAT2 simulations may be missing at the beginning of the timeseries, so
  # for simplicity, "impute" those missing values so we can avoid providing NAs
  # to JAGS, but the imputed values will not be used in the model:
  soilwat.interp = array(data=NA, dim = dim(SWC_soilwat))
  delta.soilwat.data = array(data=NA, dim = dim(SWC_soilwat))
  for(v in 1:nc){
    soilwat.interp[,v] <- na_interpolation(SWC_soilwat[,v], option = "linear")
  }
  # Compute difference in simulated SOILWAT2 SWC between day t and day (t-1):
  delta.soilwat.data[1,] = 0
  delta.soilwat.data[2:nr,] = soilwat.interp[2:nr,]-soilwat.interp[1:(nr-1),]
  

  # Now, simulate pseudo data by randomly drawing time periods for the original
  # timeseries. Repeat Nrep times for each soil depth. Randomly pick a day to 
  # start the missing data gap, then the length of the missing data gap is 
  # determined by the gsim values.
  # Get original observed SWC at the start and end; if NAs, then try again.
  # Use the observed values at the start and end to generate linearly interpolated values.
  # Grab the SOILWAT2 values and covariate data for the entire "missing data" 
  # period (start to finish).
  # Grab the observed SWC for the entire period (it's okay if it has some NAs in
  # the middle).
  obs.swc = c() # SWC drawn from observed timeseries
  lin.int.swc = c() # Linear interpolated SWC from start to end of period
  temp.lin.int = c() # Linearly interpolated temperature (if there are missing)
  sim.swc = c() # SOILWAT2 for period
  delta.sim.swc = c() # Change in SOILWAT2 for period
  ppt.swc = c()
  temp.use.swc = c()
  glength = c()
  gstart = c()
  gnum = c()
  depth.id = c()
  # Loop through all soil depths:
  for(v in 1:nc){
    nn = 0
    # Set range of possible values for the start or end of missing data period:
    SWCstart = max(obs.start[v], sim.start[v]) + 1
    SWCend = min(obs.end[v],sim.end[v])
    # Loop through all simulated gap periods:
    for(j in 1:Nrep[v]){
      nn = 0
      # Gap length (days):
      gap.length = gsim[j]
      while(nn == 0){
        # Randomly pick a day to start the missing data gap period:
        gap.start = rdunif(1,SWCstart,max(SWCstart+2,(SWCend-gap.length-2)))
        # Day to end gap:
        gap.end = gap.start+gap.length+1
        # Does the first or last day for the selected period have an NA for SWC?
        testNA = sum(is.na(SWC_obs[c(gap.start,gap.end),v])) + 
          sum(is.na(SWC_soilwat[c(gap.start,gap.end),v]))
        if(testNA > 0 | gap.end > SWCend){
          # If so, start over:
          nn = 0
        }
        else {
          # Get all data associated with the gap period:
          nn = 1
          obs.swc = c(obs.swc,SWC_obs[gap.start:gap.end,v])
          sim.swc = c(sim.swc,soilwat.interp[gap.start:gap.end,v])
          delta.sim.swc = c(delta.sim.swc,soilwat.interp[gap.start:gap.end,v] - 
                              soilwat.interp[(gap.start-1):(gap.end-1),v])
          ppt.swc = c(ppt.swc,SWC_x$P[gap.start:gap.end])
          temp.use.swc = c(temp.use.swc,Tuse[gap.start:gap.end])
          temp.lin.int = SWC_obs[SWCstart,v]
          # Linear interpolation of SWC and temperature from start to end of gap period:
          for(t in 1:(gap.length+2)){
            m = (SWC_obs[gap.end,v] - SWC_obs[gap.start,v])/((gap.length+2) - 1)
            temp.lin.int[t] = SWC_obs[gap.start,v] + m*(t-1)
          }
          lin.int.swc = c(lin.int.swc,temp.lin.int)
          # Gap length:
          glength = c(glength,rep(gap.length,gap.length+2))
          # First day of gap period:
          gstart = c(gstart,rep(gap.start,gap.length+2))
          # Soil depthID
          depth.id = c(depth.id,rep(v,gap.length+2))
          # Gap ID number:
          gnum = c(gnum,rep(j,gap.length+2))
        }
      }
    }
  }
  
  # Create mean daily precip since start of gap (pptday),
  # mean daily temp since start of gap (tempday), and time since
  # start of gap (gapday)
  gapday =c()
  pptday = c()
  tempday = c()
  pptcheck = c()
  gapday[1] = 0
  pptday[1] = ppt[gstart[1]]
  tempday[1] = Tuse[gstart[1]] 
  pptcheck[1] = ppt[gstart[1]]
  for(j in 2:length(gnum)){
    if(gnum[j]!=gnum[j-1]){
      gapday[j] = 0
      pptday[j] = ppt[gstart[j]]
      tempday[j] = Tuse[gstart[j]]
      pptcheck[j] = ppt[gstart[j]]
    }
    if(gnum[j]==gnum[j-1]){ 
      gapday[j] = gapday[j-1]+1
      pptday[j] = mean(ppt[gstart[j]:(gstart[j]+gapday[j])])
      tempday[j] = mean(Tuse[gstart[j]:(gstart[j]+gapday[j])])
      pptcheck[j] = ppt[(gstart[j]+gapday[j])]
    }
  }
  
  ave.temp = mean(Tuse,na.rm=TRUE)
  sd.temp = sd(Tuse,na.rm=TRUE)
  
  # Missing data indicator:
  gmiss = c()
  N = length(gapday)
  # gmiss = 1: day before/after gap (not part of missing data gap)
  # gmiss = 2: during the missing data gap.
  gmiss[1] = 1
  gmiss[N] = 1
  for(i in 2:(N-1)){
    # 1 = day before/after gap (not part of missing data gap)
    if(gapday[i] == 0 | gapday[i+1] == 0) gmiss[i] = 1 
    else gmiss[i] = 2
  }

  save(nc, nr, Nrep, glength, gstart, gnum, gapday, pptday, tempday, 
       pptcheck, depth.id, ave.temp, sd.temp, gmiss,
       obs.swc, sim.swc, delta.sim.swc, ppt.swc,
       temp.use.swc, lin.int.swc, file = paste0(site,".RData"))
  
}
