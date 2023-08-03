# Proof of Concept Satellite Ammonia Project

## Background
This repository contains a sample of code produced for the STFC Air Quality Network funded proof of concept project titled: 'Improving Satellite Observations of Ammonia by Integrating Chemical Transport Modelling'. This project was a collaboration between Imperial College London (ICL) and the Rutherford Appleton Laboratory (RAL). 

The project aimed to improve the spatial and temporal scale of satellite-derived ammonia (NH3) available to RAL and the wider research and policy community. The specific objectives were:
1.	Establish the current level of agreement between satellite-derived NH3 estimates with state-of-the-art regional chemical transport model (CMAQ).
2.	Develop satellite NH3 estimates at improved spatial and temporal scales through the integration of modelled vertical profiles.
3.	Improve the agreement between satellite observations and ground-level measurements for individual episodes and longer timescales with a focus on intensive agricultural activity.

The satellite data, provided by the remote sensing group at RAL Space, used their Infrared and Microwave Sounder (IMS) retrieval scheme applied to the Cross-track Infrared Sounder (CrIS) and Across-track Microwave Sounder (ATMS) on the Suomi-NPP satellite. Each individual sounding of NH3 column-average volume mixing ratio was accompanied by its corresponding averaging kernel (AK) describing vertical sensitivity along with estimated random error and auxiliary information on cloud, temperature profile, other co-retrieved variables and data quality indicators.

The satellite vertical sensitivity to NH3 was investigated by integrating with the chemical transport model CMAQ data. Each satellite sounding was co-located with the corresponding CMAQ grid cell. The satellite averaging kernels were then applied to the z-dimension of the CMAQ data at the equivalent heights and the column average mixing ratio (CAMR) was calculated (CMAQ with Satellite sensitivity CAMR) from the resulting profile. The CMAQ with Satellite sensitivity CAMR was compared to the satellite CAMR and CMAQ CAMR. 

The comparisons were conducted using CMAQ model data gridded at spatial resolutions of 10km and 50km, for each season in the years 2018-2019. Different satellite filters were applied to the CMAQ data, investigating the effect of using various levels of cloud fraction and thermal contrast between the surface and 1km altitude level. Cloud fractions of 5%, 20%, 50% and 100% were applied as well as a temperature contrast filter of greater than 10°K (DT1000). All the tests included a measure of the spectral fit precision (cost function) too (only values <1000 were permitted).

Below is an example of CMAQ with satellite sensitivity CAMR for summer (june, july and august) season in 2018 using 10km CMAQ domain and 20% cloud filter. 
![alt text](https://github.com/annalisasheehan/PoC_Ammonia_Sample_Code/blob/main/CMAQ_AK_camr_mean_av_filtered_2018JJA_NH3_10kmdomain_cf0.2.png?raw=true)


## Contents of repository

This repository contains the following files: 

```
- CMAQ_SatelliteAK_CAMR.py
- utils/camr_functions.py
- HPC_runscripts/HPC_CMAQakCAMR_2018_Domain10km_cf02.pbs
```

### CMAQ_SatelliteAK_CAMR.py

This file takes the satellite and CMAQ data, filters the satellite data based on a cloud filter and cost function, and calculates seasonal averages for the CMAQ with satellite sensitivty CAMR values. 

### utils/camr_functions.py

This file contains the functions used in CMAQ_SatelliteAK_CAMR.py

### HPC_runscripts/HPC_CMAQakCAMR_2018_Domain10km_cf02.pbs

This file contains the specific HPC job request variables along with the command to run CMAQ_SatelliteAK_CAMR.py specifying the year, cloud fraction and CMAQ model domain size parameters. 