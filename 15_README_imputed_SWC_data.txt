Readme created on 11/13/2023

The original data files are named 0#_US-SITE_SWC.csv.

Data files containing the observed volumetric soil water content (SWC) when reported 
and containing the imputed missing SWC (during missing data gaps) are named:
"##_imputed_SWC_SITE.csv"

Example data (observed and imputed SWC) files are provided to match the dates (rows) in the original example data files
(named "0#_US-SITE_SWC.csv"), such that 100 observations of SWC are reported for at least one of the soil layers.
Other layers may have less than 100 reported observations, potentially zero observations, in the example datasets. If
SWC values are missing in the example data, these imputed SWC files contain the imputed SWC data, producing complete
timeseries of SWC for the subset of example data provided, unless the original example data included zero observations
for a particular soil layer, in which case, imputed values are not available in the example file. These example
data are simply provided to illustrate the imputed data products produced by the Bayesian mixture model and associated
R code.

SITE = mpj, seg, ses, vcp, vcs, and wjs

The "imputed" SWC data files are based on the information provided in the "original
data files" (noted above).

Columns in imputed SWC data files:

* Date, Year, DOY, Month, Day: Original date information provided in the soil water 
	content data from Marcy's lab group.
* code_X: Code indicating if data and depth X associated with a missing data gaps 
	(code = 2 if during a missing data gap), or not associated with a missing data
	gap (code = 1 if observed SWC reported, or code = NA if prior to the start of data
	collection). X denotes the "depth id" (X = 1, 2, 3, ...); see below for definitions.
* mean_X: Mean SWC for "depth id" X. When SWC was missing from the original data,
	the missing value was imputed based on the Bayesian imputation mixture model (i.e.,
	mixture of SOILWAT2 adjusted value and linearly interpolated value). The posterior mean
	of all imputed values (across all mcmc samples) is reported as "mean_X". When the 
	observed SWC was NOT missing, mean_X = observed SWC value (original data).
* sd_X: Standard deviation of the imputed missing value for "depth ID" X (posterior std dev 
	across all mcmc samples); see mean_X above for more details. When the observed SWC value 
	was not missing (i.e., reported), sd_X = 0.
* 2.5th_X: 2.5th percentile of imputed SWC values for "depth ID" X (based on all mcmc 
	samples of the imputed missing value).
* 97.5th_X: 97.5th percentile of imputed SWC values for "depth ID" X (based on all mcmc 
	samples of the imputed missing value).

**If mean_X, sd_X, 2.5th_X, and 97.th_X = NA, then there was insufficient information to
impute the missing data (e.g., can't impute missing SWC values at the beginning of a time
series; can only impute missing values sandwiched between observed values).


Site-specific depth IDs:

SES (5 depths):
X = 1 = 2 cm
X = 2 = 12 cm
X = 3 = 22 cm
X = 4 = 37 cm
X = 5 = 52 cm

SEG (5 depths):
X = 1 = 2 cm
X = 2 = 12 cm
X = 3 = 22 cm
X = 4 = 37 cm
X = 5 = 52 cm

WJS (3 depths):
X = 1 = 5 cm
X = 2 = 10 cm
X = 3 = 30 cm

MPJ (3 depths):
X = 1 = 5 cm
X = 2 = 10 cm
X = 3 = 30 cm 

VCP (6 depths):
X = 1 = 5 cm
X = 2 = 10 cm
X = 3 = 20 cm
X = 4 = 30 cm
X = 5 = 50 cm
X = 6 = 60 cm

VCS (4 depths):
X = 1 = 5 cm
X = 2 = 10 cm
X = 3 = 30 cm
X = 4 = 60 cm