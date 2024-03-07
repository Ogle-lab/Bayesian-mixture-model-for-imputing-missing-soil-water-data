Information for all data files named "0#_US-SITE_SWC.csv"

Example data files are provided such that 100 observations of SWC are reported for at least one of the soil layers.
Other layers may have less than 100 observations, potentially zero observations, in the example datasets. The example
data are simply provided to illustrate how site-level SWC data are used in the associated R code and Bayesian model
code files.

- one csv spreadsheet per site, with additional columns added to the below formatted spreadsheet of the SOILWAT2 output
- note that US-Ses and US-Seg have columns in a different order because the flux tower data predates the SOILWAT2 output
    - Obs_VWC_i-jcm = observed volumetric water content in soil layers at i-j cm depth
    - P = flux tower observed precipitation


20220711_NMEG_ObservationalPeriod_SoilCalibrate-SW2v700-FXW-neuroFX2021-r6-mix1

  * csv-formatted spreadsheet (one per simulation run):
    - rows (days)
    - columns:
      - Time: (Gregorian calendar) Year; day of year "DOY"; Month; Day
      - Daily forcing input: Input_AirTemp_max_C; Input_AirTemp_min_C; Input_PPT_mm
      - Daily simulation output from SOILWAT2:
        - Sim_SWE_mm = snowpack [SWE mm] (snow water equivalents)
        - Sim_Hoh_MJm-2 = extraterrestrial horizontal solar irradiation [MJ/m2]
        - Sim_Hgt_MJm-2 = global tilted irradiation [MJ/m2] = total incoming radiation at a site accounting for atmosphere, topography (slope, aspect, and elevation), and geographic location
        - Sim_Infiltration_mm = infiltration of rain and snowmelt into soil surface [mm]
        - Sim_DiffuseRecharge_mm = water percolation below simulated soil profile [mm]
        - Sim_ET_mm = evapotranspiration [mm]
        - Sim_T_mm = transpiration [mm]
        - Sim_E_mm = evaporation [mm]
        - Sim_E_Snowloss_mm = sublimation from snowpack [mm]
        - Sim_E_Baresoil_mm = evaporation from bare soil [mm]
        - Sim_E_InterceptedCanopy_mm = evaporation from water intercepted by vegetation [mm]
        - Sim_E_SurfaceWater_mm = evaporation from water intercepted by litter or ponded at surface [mm]
        - Sim_SurfaceTemp_C_lvl = surface temperature [C] for lvl = {min, avg, max}
        - Sim_SoilTemp_C_lvl_j_cm = soil temperature [C] at j cm depth for lvl = {min, avg, max}
        - Sim_SWAat30bar_mm_itoj_cm = amount of available water [mm] in soil layers at i-j cm depth held at > -3.0 MPa
        - Sim_VWC_itoj_cm = volumetric water content [m3/m3] in soil layers at i-j cm depth
        - Sim_SWP_MPa_itoj_cm = soil water potential [MPa] in soil layers at i-j cm depth

