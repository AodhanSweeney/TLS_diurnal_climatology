import numpy as np
import pandas as pd

def to_bin_lat(y):
    lat_step = 5.
    binned_lat = np.floor(y / lat_step) * lat_step
    return(binned_lat)

def to_bin_lon(x):
    lon_step = 10.
    binned_lon = np.floor(x / lon_step) * lon_step
    return(binned_lon)

def monthly_mean_map(season_of_cosmic_data):
    """
    Description
    ===========
    Find the diurnal mean of the TLS temperatures in each 5x10 degree box
    
    Arguments
    =========
    
    season_of_cosmic_data
    --------
    Pandas dataframe of GPS-RO TLS temperatures after ERA-5 daily mean has been subtracted
    
    returns
    =======
    monthly_mean_map
    ----------------
    Diurnal mean for one month of TLS temperature in every 5x10 box
    """
    
    #Group cosmic data by lats and lons
    season_of_cosmic_data["latbin"] = season_of_cosmic_data.Lat.map(to_bin_lat)
    season_of_cosmic_data["lonbin"] = season_of_cosmic_data.Lon.map(to_bin_lon)
    grid_groups = season_of_cosmic_data.groupby(("latbin", 'lonbin'))
    
    # Set lats and lons
    lats = np.arange(-90, 90., 5.)
    lons = np.arange(-180, 180., 10.)
    
    monthly_mean_map = []
    for i in range(0, len(lats)):
        monthly_mean_lat = []
        for j in range(0, len(lons)):
            try:
                #Subset monthly data by lat lon and drop NaN values
                key = (lats[i], lons[j])
                cosmic_grid_cell = grid_groups.get_group(key)
                cosmic_grid_cell.dropna(subset=["Temp"], inplace=True)

                #Find monthly average temperatures
                #NOTE: Because we are looking for deviations from daily means
                #      we must carefully remove the true daily average by averaging
                #      over all 8 local time hour bins
                daily_average_bins = []
                for hour_bin in np.linspace(0, 21, 8):
                    hour_bin_1 = cosmic_grid_cell[cosmic_grid_cell['Hour'] == int(hour_bin)]
                    hour_bin_2 = cosmic_grid_cell[cosmic_grid_cell['Hour'] == (int(hour_bin)+1)]
                    hour_bin_3 = cosmic_grid_cell[cosmic_grid_cell['Hour'] == (int(hour_bin)+2)]
                    hour_bin_total = pd.concat([hour_bin_1, hour_bin_2, hour_bin_3])
                    cosmic_temp_hour_bin = hour_bin_total['Temp'].mean()
                    daily_average_bins.append(cosmic_temp_hour_bin)

                daily_mean_temp = np.nanmean(daily_average_bins)
                monthly_mean_map.append([lats[i], lons[j], daily_mean_temp])
            except:
                t_list = [lats[i], lons[j], np.NaN]
                monthly_mean_map.append(t_list)
        
    return monthly_mean_map

def mean_map_removal(month_of_data_df, mean_TLS_map):
    """
    Description
    ===========
    Radio Occultation data differs from ERA-5 data because of bias inherent to ERA-5 in lower stratosphere. 
    mean_map_removal attempts to remove any bias in the diurnal cycles after the daily ERA-5 mean has been
    removed so that this error is not aliased into the diurnal cycles.
    
    Arguments
    =========
    
    month_of_data_df
    --------
    Pandas dataframe of GPS-RO TLS temperatures after ERA-5 daily mean has been subtracted
    
    mean_TLS_map
    ------------
    a 36x36 map of diurnal means to be removed from month_of_data_df
    
    returns
    =======
    
    era_5_df_new
    ------------
    Pandas dataframe with correctly wrapped longitude
    """
    
    # Set lats and lons
    lats = np.arange(-90, 90., 5.)
    lons = np.arange(-180, 180., 10.)
    
    mean_TLS_list = np.reshape(mean_TLS_map, (1296, 3)) # our map is 5degx10deg so 36binsx36bins=1296 globally
    mean_TLS_df = pd.DataFrame(mean_TLS_list, columns=['Lat', 'Lon', 'Temp'])

    # Make sure data is properly binned
    month_of_data_df["latbin"] = month_of_data_df.Lat.map(to_bin_lat)
    month_of_data_df["lonbin"] = month_of_data_df.Lon.map(to_bin_lon)
    
    occultations_mean_removed_roladex = []
    for lat_idx in range(0, len(lats)):
        lat_key = lat_idx*5. - 90.
        occultations_at_lat = month_of_data_df[month_of_data_df['latbin'] == lat_key]
        monthly_mean_at_lat = mean_TLS_df[mean_TLS_df['Lat'] == lat_key]
        for lon_idx in range(0, len(lons)):
            lon_key = lon_idx*10. - 180.
            occultations_in_box = occultations_at_lat[occultations_at_lat['lonbin'] == lon_key]
            monthly_mean_in_box = monthly_mean_at_lat[monthly_mean_at_lat['Lon'] == lon_key]
            occultations_in_box['Temp'] = occultations_in_box['Temp'] - monthly_mean_in_box['Temp'].iloc[0]
            occultations_mean_removed_roladex.append(occultations_in_box)
      
    mean_removed_df = pd.concat(occultations_mean_removed_roladex)
    return mean_removed_df

def era5_df_switcher(era_5_np):
    """
    Description
    ===========
    DJF data was initially processed without recognizing the longitudinal wrapping was off by 180degrees. Instead
    of reprocessing all data, just realign the data if needed. 
    
    Arguments
    =========
    
    era_5_np
    --------
    Numpy array of post processed era 5 data with incorrect longitudinal wrapping. 
    
    returns
    =======
    
    era_5_df_new
    ------------
    Pandas dataframe with correctly wrapped longitude
    """
    
    era_5_df = pd.DataFrame(era_5_np, columns=['Day', 'Hour', 'Year', 'Lat', 'Lon', 'Temp'])
    era_5_df_pos_lons = era_5_df[era_5_df['Lon'] >= 0]
    era_5_df_neg_lons = era_5_df[era_5_df['Lon'] < 0]
    era_5_df_pos_lons['Lon'] = era_5_df_pos_lons['Lon'] - 180.
    era_5_df_neg_lons['Lon'] = era_5_df_neg_lons['Lon'] + 180
    era_5_df_new = pd.concat([era_5_df_pos_lons, era_5_df_neg_lons])
    return era_5_df_new

def background_and_bias_remover(month_of_occultation_df, month_of_era5_df):
    """
    Description
    ===========
    background_and_bias_remover is a function used to remove both high frequency synoptic variability and low 
    frequency bias from GPS-RO data in an effort to create a composite of TLS diurnal cycles. The function
    returns a pandas DataFrame with processed data to be used for the composite.
    
    
    Arguments
    =========
    
    month_of_occultation_df 
    -----------------------
    month_of_occultation_df is a pandas DataFrame filled with TLS data from GPS-RO occultations for one month,
    spanning 2006-2020. TLS data has information of year, day, local time hour, latitude, longitude, and the TLS
    temperature
    
    month_of_era5_df
    ----------------
    month_of_era5_df is a dataframe with TLS data derived from interpolated ERA-5 reanalysis data. Like 
    month_of_occultation_df, month_of_era5_df is just data from 2006-2020 for one month. Data has information
    of year, day, local time hour, latitude, longitude, and the TLS temperature
    
    returns
    =======
    TLS temperatures after both the ERA-5 daily mean has been removed, and the monthly mean of the diurnal
    cycles has been removed
    """
    # Set lats and lons
    lat_bins = np.arange(-90, 90., 5.)    
    lon_bins = np.arange(-180, 180., 10.)

    bias_and_ERA5_removed_data = []
    #ERA5_removed_data = []
    for year in range(2006, 2021):
        print(year)
        year_of_occultations = month_of_occultation_df[month_of_occultation_df['Year'] == year]
    
        year_of_occultations["latbin"] = year_of_occultations.Lat.map(to_bin_lat)
        year_of_occultations['lonbin'] = year_of_occultations.Lon.map(to_bin_lon)
        year_of_era5 = month_of_era5_df[month_of_era5_df['Year'] == str(year)]
        month_of_daily_era_5_removed = []
        if year_of_era5.Temp.size > 0:
            for day in year_of_occultations['Day'].unique():
                day_of_occultations = year_of_occultations[year_of_occultations['Day'] == day]
                day_of_era5 = year_of_era5[year_of_era5['Day'] == day]

                for lat in lat_bins:
                    occultations_at_lat = day_of_occultations[day_of_occultations['latbin'] == lat]
                    era5_at_lat = day_of_era5[day_of_era5['Lat'] == lat]
                    for lon in lon_bins:
                        
                        occultations_at_lon = occultations_at_lat[occultations_at_lat['lonbin'] == lon]
                        era5_at_lon = era5_at_lat[era5_at_lat['Lon'] == lon]
                        era5_daily_temps = era5_at_lon.Temp.to_list()
                        daily_mean = np.nanmean(era5_daily_temps)

                        if occultations_at_lon.Temp.size > 0:
                            occultations_at_lon.Temp = occultations_at_lon.Temp - daily_mean
                            month_of_daily_era_5_removed.append(occultations_at_lon)

                        else:
                            continue  

            month_of_daily_era_5_removed_df = pd.concat(month_of_daily_era_5_removed)
            monthly_mean_bias_by_map = monthly_mean_map(month_of_daily_era_5_removed_df)
            monthly_mean_bias_removed = mean_map_removal(month_of_daily_era_5_removed_df,
                                                         monthly_mean_bias_by_map)
            bias_and_ERA5_removed_data.append(monthly_mean_bias_removed)
        else:
            continue
        
    bias_and_era_5_removed_df = pd.concat(bias_and_ERA5_removed_data)
    return(bias_and_era_5_removed_df)

def box_mean_remover(season_of_cosmic_data):   
    """
    Description
    ===========
    removes diurnal mean from each 5x10 degree box after all data processing. Really just a final centering step 
    
    
    Arguments
    =========
    
    season_of_cosmic_data 
    -----------------------
    Data after all daily and monthly mean removal processing
    
    
    returns
    =======
    pandas data frame where mean in every 5x10 box is zero
    """
    lats = np.arange(-90, 90., 5)
    lons = np.arange(-180, 180, 10)
    hour_map = []
    for i in range(0, len(lats)):
        lat_band_df = []
        for j in range(0, len(lons)):
            lat_band = season_of_cosmic_data[season_of_cosmic_data['latbin'] == lats[i]]
            lat_lon_box = lat_band[lat_band['lonbin'] == lons[j]]
            lat_lon_box.dropna(subset = ["Temp"], inplace=True)
            lat_lon_box['Temp'] = lat_lon_box['Temp'] - lat_lon_box.Temp.mean()
            lat_band_df.append(lat_lon_box)
        
        lat_band_df_concat = pd.concat(lat_band_df)
        lat_band_df_concat['Temp'] = lat_band_df_concat['Temp'] - lat_band_df_concat.Temp.mean()
        hour_map.append(lat_band_df_concat)        
        
    final_hour_map = pd.concat(hour_map)
    return final_hour_map

def diurnal_binner(season_of_cosmic_data):
    """
    Description
    ===========
    Creates diurnal cycles in 36x36 degree map from data produced from box_mean_remover
    
    
    Arguments
    =========
    
    season_of_cosmic_data 
    -----------------------
    Data after all daily and monthly mean removal processing, and after centering has taken place
    
    
    returns
    =======
    
    cycles_in_lats_over_globe
    -----------------
    diurnal cycles for 36 5 degree latitude bands, starting at south pole
    
    cycles_in_boxes_over_globe
    --------------------------
    diurnal cycles binned into 36 latitude and 36 longitude bins
    """
    
    #hour step 
    hour_step = 3
    to_bin_hour = lambda x: (np.floor(x / hour_step) * hour_step)/hour_step
    hours = np.arange(0, 8, 1)
    
    season_of_cosmic_data["hourbin"] = season_of_cosmic_data.Hour.map(to_bin_hour)
    
    cycles_in_lats_over_globe = []
    cycles_in_boxes_over_globe = []
    for lat in np.arange(-90,90,5):
        cycle_in_boxes_at_lat = []
        season_of_data_at_lat = season_of_cosmic_data[season_of_cosmic_data['latbin'] == lat]
        for lon in np.arange(-180, 180, 10):
            season_of_data_in_box = season_of_data_at_lat[season_of_data_at_lat['lonbin'] == lon]
            cycle_in_box = []
            for hour_idx in np.arange(0,8,1):
                season_of_data_in_hourbin = season_of_data_in_box[season_of_data_in_box['hourbin'] == hour_idx]
                mean_of_hourbin = season_of_data_in_hourbin.Temp.mean()
                cycle_in_box.append(mean_of_hourbin)
            cycle_in_boxes_at_lat.append(cycle_in_box)
        cycles_in_boxes_over_globe.append(cycle_in_boxes_at_lat)
        cycle_at_lat = np.nanmean(cycle_in_boxes_at_lat, axis=0)
        cycles_in_lats_over_globe.append(cycle_at_lat)
                
    return np.array(cycles_in_lats_over_globe), np.array(cycles_in_boxes_over_globe)

