# Create covariate data needed for imputing missing SWC, including:
# mean daily precip since start of gap (pptday), mean daily temp since start 
# of gap (tempday), and time since start of gap (gapday).

library(imputeTS)

# Read in the full data for each site
fulldata.ses <- read.csv("01_US-Ses_SWC.csv")
fulldata.seg <- read.csv("02_US-Seg_SWC.csv")
fulldata.wjs <- read.csv("03_US-Wjs_SWC.csv")
fulldata.mpj <- read.csv("04_US-Mpj_SWC.csv")
fulldata.vcp <- read.csv("05_US-Vcp_SWC.csv")
fulldata.vcs <- read.csv("06_US-Vcs_SWC.csv")

# Now, create observed data files needed for running script to impute
# missing soil water content.

site.names = c("mpj", "seg", "ses", "vcp", "vcs", "wjs")
# Loop through each site
for (s in 1:6) {
  
  site = site.names[s]
  
  ## STEP 1: Compute standard deviations and means use for covariate 
  # standardizing in the Bayesian model that was applied to pseudo-data, 
  # thus, compute sd's and means based on the pseudo-data that were simulated.
  sd.cum.ppt = c()
  mean.tempday = c()
  sd.tempday = c()
  sd.days.ssog = c()

  # Load pseudo-data for site; the pseudo-data is created by running
  # "08_Script_to_simulate_pseudo_data_FINAL.R"
  load(paste0("./AuxiliaryData/",site, ".RData"))
  # Compute days since start of gap:
  days.ssog = gapday - 1
  # Compute cumulative precipitation since start of gap:
  cum.ppt = pptday * (days.ssog + 2)
  
  # Compute sd's and mean for cumulative precipitation (sd), mean daily
  # temperature (sd and mean), and days since start of gap (sd) based on
  # psuedo-data;
  sd.cum.ppt = sd(cum.ppt, na.rm = TRUE)
  mean.tempday = mean(tempday, na.rm = TRUE)
  sd.tempday = sd(tempday, na.rm = TRUE)
  sd.days.ssog = sd(days.ssog, na.rm = TRUE)
  
  # Save above environmental stats (sd's and means)
  save(sd.cum.ppt, mean.tempday, sd.tempday, sd.days.ssog, 
       file = paste0("envstats_", site, ".RData"))
  
  ## STEP 2: Create observed data needed for missing data imputation script.
  ## Much of this code looks similar to the code for simulating pseudo-data
  ## in "08_Script_to_simulate_pseudo_data_FINAL.R"
  
  # Get original data for site:
  fulldata = get(paste0("fulldata.", site))
  
  # Load the indices for grabbing each site's observed SWC, SOILWAT2 SWC,
  # and environmental (precip and temp) data:
  source("./AuxiliaryData/site_data_indices.R")
  obs.id <- get(paste0("obs.", site))
  sim.id <- get(paste0("sim.", site))
  env.id <- get(paste0("env.", site))
  
  # Variable for matching observed soil depths to corresponding SOILWAT2 layers
  source("./AuxiliaryData/soilwat_match_depths.R")
  match = get(paste0("match.", site))
  
  # Get environmental variables:
  SWC_x <- fulldata[, env.id]
  
  # Pick temperature variable to use (same as used for simulating pseudo-data):
  Tuse.temp <- SWC_x$Input_AirTemp_max_C

  # Get observed SWC
  SWC_obs <- fulldata[, obs.id]
  nr = dim(SWC_obs)[1] # number of days
  nc = dim(SWC_obs)[2] # number of depths
  
  # Create linearly interpolated SWC values for missing data gaps, to be used
  # in the mixture model to impute missing data:
  linear.interp = array(data = NA, dim = dim(SWC_obs))
  for (v in 1:nc) {
    linear.interp[, v] <-
      na_interpolation(SWC_obs[, v], option = "linear")
  }
  # Find the row or observation ID associated with the first and last
  # (non-missing) SWC:
  obs.start = rep(1, nc)
  obs.end = rep(1, nc)
  for (v in 1:nc) {
    obs.start[v] <- min(which(SWC_obs[, v] > 0))
    obs.end[v] <- max(which(SWC_obs[, v] > 0))
  }
  
  # Get SOILWAT2 simulated SWC (if observated layer lies between 2 SOILWAT2 
  # layers, then get the mean of the 2 SOILWAT2 layers):
  SWC_soilwatX <- fulldata[, sim.id]
  SWC_soilwat <- SWC_soilwatX[, 1:2]
  for (v in 1:nc) {
    SWC_soilwat[, v] <-
      (SWC_soilwatX[, match[v, 1]] + SWC_soilwatX[, match[v, 2]]) / 2
  }
  # Find the row or observation ID associated with the first and last simulated
  # (non-missing) SWC from SOILWAT2
  sim.start = rep(1, nc)
  sim.end = rep(1, nc)
  for (v in 1:nc) {
    sim.start[v] <- min(which(SWC_soilwat[, v] > 0))
    sim.end[v] <- max(which(SWC_soilwat[, v] > 0))
  }
  # SOILWAT2 simulated SWC may be missing at the beginning of the timeseries, so
  # for simplicity, "impute" those missing values so we can avoid providing NAs
  # to the imputation script:
  soilwat.interp = array(data = NA, dim = dim(SWC_soilwat))
  delta.soilwat.data = array(data = NA, dim = dim(SWC_soilwat))
  for (v in 1:nc) {
    soilwat.interp[, v] <-
      na_interpolation(SWC_soilwat[, v], option = "linear")
  }
  # Compute difference in simulated SOILWAT2 between day t and day (t-1):
  delta.soilwat.data[1, ] = 0
  delta.soilwat.data[2:nr, ] = soilwat.interp[2:nr, ] - soilwat.interp[1:(nr -
                                                                            1), ]
  # Get / define precipitation data:
  if (site == "vcs")
    SWC_x$P = SWC_x$Input_PPT_mm
  ppt = SWC_x$P
  
  # Define temperature data, and impute missing values (rare) via linear interpolation:
  Tuse <- na_interpolation(Tuse.temp, option = "linear")
  
  # Create mean daily precip since start of gap (pptday), mean daily temperature
  # since start of gap (tempday), time since start of gap (gapday), and indicator
  # for whether or not a day has a missing SWC value (gmiss):
  gapday = array(data = NA, dim = dim(SWC_obs))
  gmiss = array(data = NA, dim = dim(SWC_obs))
  pptday = array(data = NA, dim = dim(SWC_obs))
  tempday = array(data = NA, dim = dim(SWC_obs))
  for (v in 1:nc) {
    gapday[obs.start[v], v] = 0
    pptday[obs.start[v], v] = ppt[obs.start[v]]
    tempday[obs.start[v], v] = Tuse[obs.start[v]]
    gmiss[obs.start[v], v] = 1
    
    for (j in (obs.start[v] + 1):obs.end[v]) {
      if (is.na(SWC_obs[j, v]) == FALSE) {
        # Not a missing data gap
        gapday[j, v] = 0
        gmiss[j, v] = 1
        pptday[j, v] = ppt[j]
        tempday[j, v] = Tuse[j]
      } else{
        # During a gap
        gmiss[j, v] = 2
        if (is.na(SWC_obs[j, v]) == TRUE &
            is.na(SWC_obs[j - 1, v]) == FALSE) {
          # First day of a missing data gap
          gapday[j, v] = 1
        } else{
          # Subsequent days in a missing data gap
          gapday[j, v] = gapday[j - 1, v] + 1
        }
        # Mean precipitation since start of gap:
        pptday[j, v] = mean(ppt[(j - gapday[j, v]):j])
        # Mean temperature since start of gap:
        tempday[j, v] = mean(Tuse[(j - gapday[j, v]):j])
      }
    }
  }
  
  # Save site-specific observed data:
  save(nc, nr, gapday, pptday, tempday, gmiss, obs.start, obs.end, SWC_obs,
    linear.interp, soilwat.interp, delta.soilwat.data,
    file = paste0(site, "_obs_data.RData"))
}
