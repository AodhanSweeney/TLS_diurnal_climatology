import numpy as np
import pandas as pd
import dask.dataframe as dd

def to_bin_hour(h):
    hour_step = 3
    binned_hour =int((np.floor(h / hour_step) * hour_step)/hour_step)
    return(binned_hour)

def daily_mean_remover(era_5_df):
    """
    Description
    ===========
    Find the diurnal mean of the ERA-5 TLS temperatures in each 5x10 degree box and remove it
    
    Arguments
    =========
    
    era_5_df
    --------
    Pandas dataframe of ERA-5 TLS temperatures.
    
    returns
    =======
    daily_mean_removed_df
    ----------------
    Pandas dataframe of ERA-5 TLS temperatures after daily mean has been removed
    """
    daily_mean_removed = []
    for year in range(2006, 2021):
        print(year)
        era_5_df_year = era_5_df[era_5_df['Year'] == str(year)]
        for day in era_5_df_year.Day.unique():
            era_5_day = era_5_df_year[era_5_df_year['Day'] == day]
            for lat in np.arange(-90, 90, 5):
                era_5_lat = era_5_day[era_5_day['Lat'] == lat]
                for lon in np.arange(-180, 180, 10):
                    era_5_lon = era_5_lat[era_5_lat['Lon'] == lon]
                    era_5_lon.Temp = era_5_lon.Temp - era_5_lon.Temp.mean()
                    daily_mean_removed.append(era_5_lon)
    daily_mean_removed_df = pd.concat(daily_mean_removed)
    return (daily_mean_removed_df)


def df_organizer(daily_mean_removed_df):
    """
    Description
    ===========
    Sorts dataframe so that chunking logic is consistent with where means will be drawn. Chunks are 
    created by latitude band.
    
    Arguments
    =========
    
    daily_mean_removed_df
    --------
    Pandas dataframe of ERA-5 TLS temperatures after the daily mean has been removed.
    
    returns
    =======
    daily_mean_removed_dd
    ----------------
    Dask dataframe of ERA-5 TLS temperatures after daily mean has been removed
    """
    daily_mean_removed_df['Day'] = daily_mean_removed_df['Day'].astype(int)
    daily_mean_removed_df['Hour'] = daily_mean_removed_df['Hour'].astype(int)
    daily_mean_removed_df['Year'] = daily_mean_removed_df['Year'].astype(int)
    daily_mean_removed_df['Lat'] = daily_mean_removed_df['Lat'].astype(int)
    daily_mean_removed_df['Lon'] = daily_mean_removed_df['Lon'].astype(int)
    daily_mean_removed_df['Temp'] = daily_mean_removed_df['Temp'].astype(float)
    
    # Sort data frame by lat so chunks will align correctly
    sorted_by_lat_df = daily_mean_removed_df.sort_values(by='Lat')
    daily_mean_removed_df = sorted_by_lat_df.reset_index(drop=True)
    daily_mean_removed_dd = dd.from_pandas(daily_mean_removed_df, npartitions=36, sort=False)
    return(daily_mean_removed_dd)

def diurnal_binner(season_of_cosmic_data):
    """
    Description
    ===========
    A variation of the gpsro_tools.diurnal_binner() function. This version uses dask to compute diurnal cycles because the ERA-5
    binning can be tedious because of the large amount of data
    
    Arguments
    =========
    
    season_of_cosmic_data
    --------
    Dask dataframe of ERA-5 TLS temperatures after the daily mean has been removed.
    
    returns
    =======
    cycles_in_boxes_over_globe
    ----------------
    List where diurnal cycles are binned into 36 latitude and 36 longitude bins.
    """
    cycles_in_boxes_over_globe = []
    for lat in np.arange(-90,90,5):
        cycle_in_boxes_at_lat = []
        season_of_data_at_lat = season_of_cosmic_data.loc[season_of_cosmic_data['Lat'] == lat]
        
        for lon in np.arange(-180, 180, 10):
            season_of_data_in_box = season_of_data_at_lat.loc[season_of_data_at_lat['Lon'] == lon]
            cycle_in_box = []
            
            for hour_idx in np.arange(0,8,1):
                season_of_data_in_hourbin = season_of_data_in_box.loc[season_of_data_in_box['hourbin'] == hour_idx]
                mean_of_hourbin = season_of_data_in_hourbin.Temp.mean()
                mean_of_hourbin = mean_of_hourbin.compute()
                cycle_in_box.append(mean_of_hourbin)
                
            cycle_in_boxes_at_lat.append(cycle_in_box)

        cycles_in_boxes_over_globe.append(cycle_in_boxes_at_lat)
        print('Lat done: ', lat)
        
    return(cycles_in_boxes_over_globe)
