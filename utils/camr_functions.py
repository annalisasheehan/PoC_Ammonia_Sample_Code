import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
from scipy import interpolate
import matplotlib
import xarray as xr
import geopandas as gpd
import math
import glob
matplotlib.use('Agg')
plt.style.use('classic')

def set_parameters(year, season, c_size, cf):
    if c_size == 50:
        d_size = '01'
        domain = 'eu'
        max_d = 0.7
    else: 
        d_size = '02'
        domain = 'uk'
        max_d = 0.1
    
    outname = '{y}{s}_NH3_{size}km'.format(y=year, s=season, size=c_size)

    if cf == 0.05:
        cf_name = '_cf005'
    elif cf == 0.2:
        cf_name = '_cf02'
    elif cf == 1:
        cf_name = '_cf1'
    elif cf == 0.5:
        cf_name = '_cf05'
    else: 
        print('Check cloud fraction variable')
        
    if season == 'MAM':
        season_nums = [3,4,5]
        first_month = season_nums[0]
        month_range = '0[3-5]'
    elif season == 'JJA':
        season_nums = [6,7,8]
        first_month = season_nums[0]
        month_range = '0[6-8]'
    elif season == 'SON':
        season_nums = [9,10,11]
        first_month = season_nums[0]
        month_range = '0[9-11]'
    else:
        print('Check season variable')
        season_nums = 'nan'
        first_month = 'nan'
        month_range = 'nan'
    return d_size, domain, max_d, outname, cf_name, season_nums, first_month, month_range

def set_parameters_sat_agg(c_size, cf):
    if c_size == int(50):
        d_size = '01'
        domain = 'eu'
        max_d = 0.7 # max distance between CMAQ and satellite points in sjoin nearest
    else:
        d_size = '02'
        domain = 'uk'
        max_d = 0.1 # max distance between CMAQ and satellite points in sjoin nearest

    # cloud fraction
    if cf == 0.05:
        cf_name = '_cf005'
    elif cf == 0.2:
        cf_name = '_cf02'
    elif cf == 1:
        cf_name = '_cf1'
    elif cf == 0.5:
        cf_name = '_cf05'
    return d_size, domain, max_d, cf_name

def set_seasons_sat_agg(year, season, c_size):
    outname = '{y}{s}_NH3_{size}km'.format(y=year, s=season, size=c_size)
    if season == 'MAM':
        season_nums = [3,4,5]
        first_month = season_nums[0]
        month_range = '0[3-5]'
    elif season == 'JJA':
        season_nums = [6,7,8]
        first_month = season_nums[0]
        month_range = '0[6-8]'
    elif season == 'SON':
        season_nums = [9,10,11]
        first_month = season_nums[0]
        month_range = '0[9-11]'
    else:
        print('Check...')
        season_nums = 'nan'
        first_month = 'nan'
        month_range = 'nan'
    return outname, season_nums, first_month, month_range


def read_cmaq_season_data(cmaq, season_nums, valid_hours):
    
    # subset to season 
    cmaq_season = cmaq.isel(time=(cmaq.time.dt.month.isin(season_nums)))
    # subset to hours for which the satellite data overlap the UK
    cmaq_season_hrs = cmaq_season.isel(time=(cmaq_season.time.dt.hour.isin(valid_hours)))
    cmaq_season_3d = cmaq_season_hrs.NH3.isel(z=0)
    # fill the grid with zeros 
    sat_cmaqGrid = xr.zeros_like(cmaq_season_3d)
    sat_cmaqGrid = sat_cmaqGrid.where(sat_cmaqGrid != 0)
    
    cmaq_time = sat_cmaqGrid.time
    # convert the xarray datasets to pandas df
    cmaq_time_df = cmaq_time.to_dataframe()
    cmaq_time_df = cmaq_time_df.rename(columns={'time': 'Col_1'})
    cmaq_time_df = cmaq_time_df.reset_index().reset_index()
    # Copy empty xarray dataset to input the satellite CAMR value
    sat_mgf_cmaqGrid = sat_cmaqGrid.copy()

    return sat_cmaqGrid, sat_mgf_cmaqGrid, cmaq_time_df

def get_list_sat_season_files(season, sat_dir, month_range):
    if season == 'MAM' or season == 'JJA':
        sat_files = glob.glob("{sdir}/{months}/*/*.nc".format(sdir=sat_dir, months= month_range))
    elif season == 'SON':
        sat_files = []
        sat_files_s = glob.glob("{sdir}/09/*/*.nc".format(sdir=sat_dir))
        sat_files_o = glob.glob("{sdir}/10/*/*.nc".format(sdir=sat_dir))
        sat_files_n = glob.glob("{sdir}/11/*/*.nc".format(sdir=sat_dir))
        sat_files = sat_files_s + sat_files_o + sat_files_n
    return sat_files

def extract_satellite_time_day(satellite):
    """
    Extract the date and time of the satellite overpass
    """
    # extract satellite sensing time and day 
    sat_day_time_coords = satellite.set_coords(('latitude', 'longitude','sensingtime_day', 'sensingtime_msec'))
    sat_day_time_coords_npi = sat_day_time_coords.npi.isel(npi=0)
    day = float(sat_day_time_coords_npi.sensingtime_day.values)
    start_date = '01/01/2000'
    start_date = datetime.strptime(start_date, "%d/%m/%Y")
    sat_day = start_date + timedelta(days=day)
    time = float(sat_day_time_coords_npi.sensingtime_msec.values)/1000
    time /= 60 *60
    sat_hour = int(math.modf(time)[1])
    return sat_day, sat_hour

def extract_satellite_npi_npres(satellite):
    """
    Extract the npi and npres dimensions which contain the variables:
    - Latitude
    - Longitude 
    - Filters such as jy and cfr.
    """
    sat_day_time_coords = satellite.set_coords(('latitude', 'longitude','sensingtime_day', 'sensingtime_msec', 'mgf', 'jy', 'btd_flag', 'cfr', 'sp', 'tsk'))
    variables = ['latitude', 'longitude','sensingtime_day', 'sensingtime_msec', 'mgf', 'jy', 'btd_flag', 'cfr', 'sp', 'tsk']
    satellite_coords_sub = sat_day_time_coords[variables]
    # subset data set based on the dimensions 
    satellite_coords_sub_npi = satellite_coords_sub.npi # sp in hPa
    satellite_coords_sub_npres = satellite_coords_sub.npres 
    # convert the xarray datasets to pandas df
    sat_df_npi = satellite_coords_sub_npi.to_dataframe() # sp in hPa
    sat_df_npres = satellite_coords_sub_npres.to_dataframe() 
    # join on the index values which correspond to the npi and npres
    sat_df = pd.concat([sat_df_npi, sat_df_npres], axis=1)
    return sat_df

def extract_CMAQ_day_2D_data(cmaq, sat_day, sat_hour):
    """
    Subset the CMAQ data frame to a single day, ammonia and the surface level.
    This data frame will be used to find the nearest satellite sounding. 
    """
    ### subset the CMAQ data 
    # select a single day to be able to find the nearest cmaq x,y cell to satellite sounding
    cmaq_sub = cmaq.sel(time='{y}-{m}-{d}T{h}:00:00.000000000'.format(y=sat_day.year, m=sat_day.month, d=sat_day.day, h=sat_hour,))
    # select ammonia 
    cmaq_nh3 = cmaq_sub.NH3 # ppm
    # select the lowest level 
    cmaq_nh3_bottom = cmaq_nh3.isel(z=0) #ppm
    # convert to df and reset index 
    CMAQ_df = cmaq_nh3_bottom.to_dataframe().reset_index()
    # convert to gdf using lat and long values for points
    CMAQ_gdf = gpd.GeoDataFrame(CMAQ_df, geometry=gpd.points_from_xy(CMAQ_df.longitude, CMAQ_df.latitude))
    return cmaq_sub, CMAQ_gdf

def extract_satellite_pressure(satellite):
    # extract satellite pressure and a.k.'s
    satellite_coords = satellite.set_coords(('latitude', 'longitude', 'p')) # p= hPa
    # set coordinates in ncdf 
    variables = ['latitude', 'longitude', 'p']
    satellite_coords_sub = satellite_coords[variables] 
    # subset data set to pressure dimension and convert to d.f.
    satellite_coords_sub_nz = satellite_coords_sub.nz
    sat_df_nz = satellite_coords_sub_nz.to_dataframe()
    # create a numpy array with all the 101 pressure levels
    p = sat_df_nz['p'].to_numpy() # hPa
    return p

def extract_satellite_ak(satellite):
    satellite_coords = satellite.set_coords(('latitude', 'longitude', 'ak_mgf'))
    variables = ['latitude', 'longitude', 'ak_mgf'] # assume ak in ppb
    satellite_coords_ak = satellite_coords[variables]
    #select out ammonia 
    nh3 = satellite_coords_ak.isel(nmgf=1)
    ak_mgf_sub = nh3.ak_mgf
    ak_df = ak_mgf_sub.to_dataframe()
    # remove the duplicat ak_mgf column 
    ak_df.columns.values[1] = "b"
    ak_df = ak_df.drop(columns=['b'])
    return ak_df

def extract_satellite_CAMR(satellite):
    satellite_camr = satellite.set_coords(('mgf', 'profiles_mgf'))
    variables = ['mgf', 'profiles_mgf']
    satellite_camr_sub = satellite_camr[variables]
    satellite_camr_sub_nh3 = satellite_camr_sub.sel(nmgf=1)
    satellite_camr_sub_nh3_npres = satellite_camr_sub_nh3.npres
    sat_df_mgf = satellite_camr_sub_nh3_npres.to_dataframe()
    # assumed profile of NH3 throughout atmosphere (e.g. 1's)
    sat_camr_sub_nh3_prof_nz = satellite_camr_sub_nh3.nz
    sat_df_profmgf = sat_camr_sub_nh3_prof_nz.to_dataframe()
    # ; ppmv stored in file, convert to ppbv
    sat_df_profmgf['profiles_mgf'] = sat_df_profmgf['profiles_mgf'] * 1000
    return sat_df_mgf, sat_df_profmgf

def extract_corresponding_index_values(sat_and_cmaq_subset, i):
    """
    Retrieving satellite npi value and CMAQ x and y locations
    """
    npi_value = sat_and_cmaq_subset['npi'].values[i]
    x_num = sat_and_cmaq_subset['x'].values[i]
    y_num = sat_and_cmaq_subset['y'].values[i]

    return npi_value, x_num, y_num

def interpolate_sat_AKs(ak_df, npi_value, ak_p_t, p):
    # select the # pixel using npi value 
    ak_df_i = ak_df.xs(npi_value, level='npiak_mgf').reset_index()
    # get a.k. onto same number of levels as pressure grid 
    # convert a.k. values to numpy array
    ak_v = ak_df_i['ak_mgf'].to_numpy()
    # interpolate A.K. to 101 levels
    f = interpolate.interp1d(ak_p_t, ak_v, fill_value="extrapolate")
    # apply the interpolation to get the a.k. values at the satellite pressure levels 
    ak_new = f(p) #hPa
    return ak_new

def interpolate_cmaq_nh3(cmaq_sub, x_num, y_num, p):
    subset = cmaq_sub.isel(x=x_num, y=y_num)
    # convert to a df
    cmaq_df = subset.to_dataframe().reset_index()
    # drop duplicates because of dimensions
    cmaq_df_nd = cmaq_df.drop_duplicates(subset=['ZH', 'ZF', 'NH3', 'PRES'])
    # inverse table based on pressure column 
    cmaq_df_nd = cmaq_df_nd.sort_values(by=['PRES']) # pa
    # convert pressure and nh3 columns to np array 
    cmaq_pres = cmaq_df_nd['PRES'].to_numpy()
    cmaq_nh3 = cmaq_df_nd['NH3'].to_numpy()
    # convert cmaq pressure from pa to hpa
    cmaq_pres = cmaq_pres / 100
    # convert cmaq ppm to ppb 
    cmaq_nh3 = cmaq_nh3 * 1000
    # interpolate cmaq nh3 values 
    f = interpolate.interp1d(cmaq_pres,cmaq_nh3, fill_value="extrapolate")
    nh3_new = f(p) #hPa and ppb
    return nh3_new

def calc_CAMR(p, model_ppbv):
    """ 
    Calculate column average mixing ratio 
    """
    # set constants 
    g=9.80665e0             # m/s2
    na=6.0221367e23         # /mol
    mmair=28.964001e-3      # kg/mol

    # calculate amount of pressure between top and bottom of each level
    p_1_ = p[1:]
    p_o_1 = p[:-1]
    dp = p_o_1 - p_1_
    
    # convert back to satellite 101 levels
    dp_end = np.append(dp, 0)
    dp_beg = np.insert(dp, 0, 0)

    # calcuate weight 
    w=(dp_beg + dp_end)/2e9*na/g/mmair/1e2
    # transform model_ppbv
    model_tc = np.dot(model_ppbv.T,w)
    # divide by the sum of the weights to convert back to ppb
    model_tc_camr = model_tc/ sum(w)

    return model_tc_camr


def filter_satellite_data(sat_day, sat_hour, cf, cmaq_time_df, satellite):

    """
    For each satellite sounding: 
    - Identify corresponding CMAQ time point index
    - Extract satellite variables to filter on 
    - Filter data using predetermined filter values for cloud fraction, cost function and air surface temperature contrast 
    Return a geodataframe with the satellite data for the filter variables
    """
    # find the corresponding cmaq time point index
    t = '{y}-{m}-{d} {h}:00:00'.format(y=sat_day.year, m=sat_day.month, d=sat_day.day, h=sat_hour)
    index_row = cmaq_time_df.loc[cmaq_time_df['time']==t]
    cmaq_t_i = index_row.index.values.tolist()[0]
        
    # extract satellite npi, npres, lat/long, filter variables
    sat_df = extract_satellite_npi_npres(satellite)
            
    # filter data on cost function and cloud fraction
    filtered_sat_df = sat_df.loc[sat_df['jy'] < 1000]
    filtered_sat_df = filtered_sat_df.loc[filtered_sat_df['cfr'] < cf] 

    # convert to gdf using lat and long values for points
    filtered_sat_gdf = gpd.GeoDataFrame(filtered_sat_df, geometry=gpd.points_from_xy(filtered_sat_df.longitude, filtered_sat_df.latitude))
    return filtered_sat_gdf, cmaq_t_i  

def extract_satellite_variables(satellite):
    #### Satellite pressure
    # create a numpy array with all the 101 pressure levels
    p = extract_satellite_pressure(satellite) # hPa
    # select out every other pressure value which will correspond to the A.K. 
    ak_p_t = p[0::2]

    # Extract satellite A.K.s from satellite sounding NetCDF
    ak_df = extract_satellite_ak(satellite)

    # Extract satellite CAMR from satellite sounding NetCDF (converting sat prof_mgf into ppb)
    sat_df_mgf, sat_df_profmgf = extract_satellite_CAMR(satellite)
    return ak_p_t, ak_df, p, sat_df_mgf, sat_df_profmgf


def calc_CMAQ_with_AK_CAMR(sat_and_cmaq_subset, i, ak_df, ak_p_t, p, cmaq_sub):
    # extract the corresponding index values
    npi_value, x_num, y_num = extract_corresponding_index_values(sat_and_cmaq_subset, i)
    
    ### extract the satellite surface pressure value
    sp_i = sat_and_cmaq_subset['sp'].loc[npi_value] # hPa

    ### Interpolate A.K.s to 101 pressure levels
    # select the # pixel using npi value 
    ak_new = interpolate_sat_AKs(ak_df, npi_value, ak_p_t, p) #hPa

    #### interpolate CMAQ pressure grid
    nh3_new = interpolate_cmaq_nh3(cmaq_sub, x_num, y_num, p)#hPa and ppb

    # multiply a.k. with corresponding CMAQ NH3 value 
    nh3_ak = nh3_new * ak_new
    # sum up the nh3 * a.k. to get CAMR 
    nh3_ak_sf = nh3_ak.sum()
    return  nh3_ak_sf, npi_value, x_num, y_num

def update_output_files(sat_and_cmaq_subset, sat_cmaqGrid, sat_mgf_cmaqGrid, sat_df_mgf, npi_value, cmaq_t_i, x_num, y_num, nh3_ak_sf):
    # # insert the camr value into the table
    sat_and_cmaq_subset.loc[npi_value, 'cmaq_camr']= nh3_ak_sf

    # insert in the corresponding satellite CAMR
    sat_camr = sat_df_mgf['mgf'].loc[npi_value]
    sat_and_cmaq_subset.loc[npi_value, 'sat_camr']= sat_camr
    
    # insert into cmaq grid
    sat_cmaqGrid[cmaq_t_i,y_num,x_num]= nh3_ak_sf
    sat_mgf_cmaqGrid[cmaq_t_i,y_num,x_num] = sat_camr 
    return sat_and_cmaq_subset, sat_cmaqGrid, sat_mgf_cmaqGrid