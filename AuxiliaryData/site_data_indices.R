# Column indices for getting obs SWC
obs.ses = 2:6
obs.seg = 2:6
obs.wjs = 69:71
obs.mpj = 69:71
obs.vcp = 69:74
obs.vcs = 69:72

# Column indices for getting SOILWAT2 SWC
sim.ses = 57:65
sim.seg = 57:65
sim.wjs = 50:58
sim.mpj = 50:58
sim.vcp = 50:58
sim.vcs = 50:58

# Column indices for getting env covariate data (P, Input_AirTemp_max_C, 
# Input_AirTemp_min_C, Input_PPT_mm):

env.ses = c(7,12:14)
env.seg = c(7,12:14)
env.wjs = c(72,5:7)
env.mpj = c(72,5:7)
env.vcp = c(75,5:7)
#env.vcs = c(73,5:7)
# Since P data (index/column 73) has lots of missing data for vcs,
# use Input_PPT_mm instead, which takes on exactly the same values as P
# when P is available.
env.vcs = c(7,5:7)