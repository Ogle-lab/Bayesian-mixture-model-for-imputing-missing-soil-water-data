## Script for evaluating mixing and convergence of MCMC chains

library(mcmcplots)

##################################
# Checking observed vs. predicted soil water content for pseudo-data:

# Set-up arrays for storing R2 and fit statistics:
R2.vals = array(data=NA,dim=c(6,1))
site.names = c("mpj","seg","ses","vcp","vcs","wjs")
rownames(R2.vals) = site.names
colnames(R2.vals) = "R2"
coef.vals = array(data=NA,dim=c(2*6,4))
rownames(coef.vals) = c("mpj","mpj","seg","seg","ses","ses","vcp","vcp",
                        "vcs","vcs","wjs","wjs")

##### !! Select site to evaluate:
site = "mpj"

# Load coda samples and posterior statistics (coda produced by running 
# 11_ScriptDynamicWeightModel_mpj_FINAL.R:
load(paste0("posterior_dynamicW_",site,".RData"))
# Load data (pseudodata):
load(paste0("./AuxiliaryData/",site,".RData"))

# Graphically check observed vs predicted SWC for pseudo-data that belonged
# to a missing data gap (i.e., SWC data treated as missing, but we know
# the actual SWC values, if not missing in the original timeseries):
n = length(obs.swc)
pred.swc = postout2.obs[1:n,1]
plot(obs.swc[gmiss==2],pred.swc[gmiss==2])
abline(a=0,b=1,col="red")
# Linear regression of observed vs predicted to check model fit:
fit = lm(obs.swc[gmiss==2] ~ pred.swc[gmiss==2],na.action = na.omit)
fit.summary = summary(fit)
fit.summary$r.squared

# Collect model fit stats:
ss = max((1:6)*match(site.names,site),na.rm=TRUE)
R2.vals[ss] = fit.summary$r.squared
coef.vals[(2*ss-1):(2*ss),] = fit.summary$coefficients

# Save model fit stats after evaluating all sites:
colnames(coef.vals) = colnames(fit.summary$coefficients)
save(R2.vals,coef.vals,file="model_fit_stats.RData")


######################################
# Check convergence and mixing

## Set-up arrays for storing convergence stats:
GelRaft = array(data=NA,dim=c(6,4))
rownames(GelRaft) = site.names
colnames(GelRaft) = c("min_psrf", "max_psrf", "median_raftery", "max_raftery")

#### !! Set "burn-in" (start) for computing Raftery and posterior stats (if needed)
r.start = 6000

# Look at history plots, autocorrelation, etc:
mcmcplot(window(coda2,thin=15))
mcmcplot(window(coda2,thin=15,start=r.start))

# Evaluate convergence via gelman diagnostic (potential scale reduction factor)
g.out = gelman.diag(coda2, multivariate = FALSE)
min(g.out$psrf[,1], na.rm = TRUE); max(g.out$psrf[,1], na.rm = TRUE)
GelRaft[ss,1:2] = c(min(g.out$psrf[,1], na.rm = TRUE), max(g.out$psrf[,1], na.rm = TRUE))

# Use Raftery diagnostic to estimates number of iterations needed
raft1<-raftery.diag(window(coda2,start=r.start)); raft1
test1<-c()
median1<-c()
for(i in 1:3){
  # For each chain, grab the maximum number of iterations needed (N)
  # across all quantities monitored:
  test1[i]<-max(raft1[[i]][[2]][,2], na.rm = TRUE)
  median1[i] <- median(raft1[[i]][[2]][,2], na.rm = TRUE)
}
# Per chain estimates of number of iterations:
median1/3
# find max number of iterations required among the 3 chains, 
# then divided by 3 (# chains) to get iterations per chain:
max(test1)/3

GelRaft[ss,3] = floor(mean(median1/3))
GelRaft[ss,4] = floor(max(test1)/3)