import monet as m
import monetio as mio
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import xarray as xr
import geopandas as gpd
import os
import argparse
from utils.camr_functions import *
matplotlib.use('Agg')
plt.style.use('classic')

# set variables from command line
parser = argparse.ArgumentParser()
parser.add_argument('--year', type=int, required=True)
parser.add_argument('--domain_size', type=int, required=True)
parser.add_argument('--cloud_fraction', type=float, required=True)
args = parser.parse_args()

# Year variable (2018, 2019 or 2020)
year = args.year
# CMAQ model domain size either 10km or 50km
c_size = args.domain_size
# % cloud cover in pixel filter
cf = args.cloud_fraction
# limit the satellite data to hours surrounding satellite overpass across the UK.
valid_hours = [10,11,12,13,14,15]

######################################################
# Set folder directories
main_dir = '/rds/general/project/erg_modelling_storage/live/SAQN_Ammonia_PoC'
sat_dir = '{m_dir}/Data/Satellite/Europe/{y}'.format(m_dir = main_dir, y=year)

plot_folder_satCAMR = '{m_dir}/Plots/Sample_Code/1_Satellite_CAMR/{size}km_fullZ_cf{cf}'.format(m_dir=main_dir, cf=cf, size=c_size)
if not os.path.exists(plot_folder_satCAMR):
    os.makedirs(plot_folder_satCAMR)
plot_folder_cmaqAK = '{m_dir}/Plots/Sample_Code/3_CMAQ_AK/{size}km_fullZ_cf{cf}'.format(m_dir=main_dir, cf=cf, size=c_size)
if not os.path.exists(plot_folder_cmaqAK):
    os.makedirs(plot_folder_cmaqAK)

data_folder_satCAMR = '{m_dir}/Data/Sample_Code/1_Satellite_CAMR/{size}km_fullZ_cf{cf}'.format(m_dir=main_dir, cf=cf, size=c_size)
if not os.path.exists(data_folder_satCAMR):
    os.makedirs(data_folder_satCAMR)
data_folder_cmaqAK = '{m_dir}/Data/Sample_Code/3_CMAQ_AK/{size}km_fullZ_cf{cf}'.format(m_dir=main_dir, cf=cf, size=c_size)
if not os.path.exists(data_folder_cmaqAK):
    os.makedirs(data_folder_cmaqAK)

######################################################
# Set parameters based on input variables
d_size, domain, max_d, cf_name = set_parameters_sat_agg(c_size, cf)

# Read in the CMAQ model domain extents
model_boxgdf = gpd.read_file('{m_dir}/Data/CMAQ_domain/{y}/{y}_d{d_size}.shp'.format(m_dir=main_dir, d_size=d_size, y=year))
model_boxgdf.to_crs(epsg=4326,inplace=True)

# Read in CMAQ data
cmaq_path = '/rds/general/project/mrc_funded_storage/live/nosha/cmaq/postproc/s00.for_annie/out/cmaq_nh3_data/netcdf'
cmaq = mio.cmaq.open_dataset('{folder}/COMBINE_ACONC_3d_nh3_v531_intel_cmaqurban_{y}_{d}{size}.nc'.format(folder = cmaq_path, y =year, size=c_size, d=domain))

# Create list of seasons to iterate through. Ignoring DJF due to low temperature contrast.
seasons_list = ['MAM','JJA', 'SON']

# Iterate through seasons to produce seasonal average CAMR
for season in seasons_list:
    print(season)

    # Set the month numbers for the season
    outname, season_nums, first_month, month_range = set_seasons_sat_agg(year, season, c_size)

    # Subset CMAQ data to season and hours of satellite overpass
    # Fill the xarray datasets with zeros ready to be filled with results
    sat_cmaqGrid, sat_mgf_cmaqGrid, cmaq_time_df = read_cmaq_season_data(cmaq, season_nums, valid_hours)

    # Get a list of satellite files within the season
    sat_files = get_list_sat_season_files(season, sat_dir, month_range)
    print(len(sat_files))

    # Iterate through satellite files extracting CAMR and other satellite variables to calculate CMAQ CAMR
    for sat_file in sat_files:

        # Open satellite sounding NetCDF
        satellite = xr.open_dataset(sat_file)

        # Extract satellite sensing time and day
        sat_day, sat_hour = extract_satellite_time_day(satellite)

        # If the satellite sounding covers the UK based on overpass time
        if sat_hour in valid_hours:
            # Extract satellite variables and filter data
            # Identify corresponding CMAQ time point index (cmaq_t_i)
            filtered_sat_gdf, cmaq_t_i = filter_satellite_data(sat_day, sat_hour, cf, cmaq_time_df, satellite)

            # Subset the CMAQ data to surface level NH3 values for the single day
            cmaq_sub, CMAQ_gdf = extract_CMAQ_day_2D_data(cmaq, sat_day, sat_hour)

            # Join satellite and CMAQ data frames
            sat_and_cmaq_subset = gpd.sjoin_nearest(filtered_sat_gdf, CMAQ_gdf, lsuffix='sat', rsuffix='CMAQ', distance_col="distances", max_distance=max_d)
            n_soundings = len(sat_and_cmaq_subset)

            # Extract satellite CAMR and variables, and calculate CMAQ with satellite sensitivity CAMR
            if n_soundings > 0:
                # Add in a new column to fill in with the CMAQ CAMR values
                sat_and_cmaq_subset['cmaq_camr'] = np.nan
                sat_and_cmaq_subset['sat_camr'] = np.nan

                # Extract satellite sounding variables
                ak_p_t, ak_df, p, sat_df_mgf, sat_df_profmgf = extract_satellite_variables(satellite)

                # Iterate through each of the filtered soundings to calculate the CMAQ CAMR
                for i in range(n_soundings):

                    # Interpolate the CMAQ NH3 data to the same pressure grid as the satellite data
                    # Apply satellite Averaging Kernels (A.K's)
                    # Calculate CMAQ with A.K CAMR
                    nh3_ak_sf, npi_value, x_num, y_num = calc_CMAQ_with_AK_CAMR(sat_and_cmaq_subset, i, ak_df, ak_p_t, p, cmaq_sub)

                    # Insert satellite and CMAQ CAMR values into output files
                    sat_and_cmaq_subset, sat_cmaqGrid, sat_mgf_cmaqGrid = update_output_files(sat_and_cmaq_subset, sat_cmaqGrid, sat_mgf_cmaqGrid, sat_df_mgf, npi_value, cmaq_t_i, x_num, y_num, nh3_ak_sf)

    # Make maps of the output and save results as NetCDF

    ### CMAQ with satellite sensitivity (A.K's) CAMR
    # Calculate the seasonal mean
    sat_cmaqGrid_mean = sat_cmaqGrid.mean(dim=['time'])
    # Plot seasonal average
    sat_cmaqGrid_mean.monet.quick_map(robust=True, figsize=(10,7), map_kws=dict(linewidth=1), cbar_kwargs={"label": "NH3 (ppb)"})
    plt.title("CMAQ*A.K. in CMAQ Filtered {out} domain".format(out=outname))
    plt.savefig('{p_folder}/CMAQ_AK_camr_mean_av_filtered_{out}domain_cf{cf}.png'.format(p_folder=plot_folder_cmaqAK, cf=cf, out=outname, size=c_size))

    ### Satellite CAMR plots
    # Calculate the seasonal mean
    sat_mgf_cmaqGrid_mean = sat_mgf_cmaqGrid.mean(dim=['time'])
    # Plot seasonal average
    sat_mgf_cmaqGrid_mean.monet.quick_map(robust=True, figsize=(10,7), map_kws=dict(linewidth=1), cbar_kwargs={"label": "NH3 (ppb)"})
    plt.title("Satellite CAMR in CMAQ Filtered {out} domain".format(out=outname))
    plt.savefig('{p_folder}/Satellite_camr_mean_av_filtered_{out}domain_cf{cf}.png'.format(p_folder=plot_folder_satCAMR, cf=cf, out=outname, size=c_size))

    # Save results as NetCDF
    sat_cmaqGrid.to_netcdf("{d_folder}/CMAQ_AK_camr_filtered_{out}_cf{cf}.nc".format(d_folder=data_folder_cmaqAK, cf=cf, out=outname, size=c_size))
    sat_mgf_cmaqGrid.to_netcdf("{d_folder}/Satellite_camr_filtered_{out}_cf{cf}.nc".format(d_folder=data_folder_satCAMR, cf=cf, out=outname, size=c_size))
