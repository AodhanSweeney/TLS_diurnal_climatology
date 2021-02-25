import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import cartopy.crs as ccrs
import cartopy.util as cutil
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib.patches as mpatches

def area_bins(diurnal_cycles_array):
    lat_centroids = np.linspace(-87.5, 87.5, 36)
    array_of_weights = np.cos(np.deg2rad(lat_centroids))
    
    weights = [array_of_weights[:6], array_of_weights[6:10], array_of_weights[10:14], array_of_weights[14:16],
               array_of_weights[16:18], array_of_weights[18:20], array_of_weights[20:22],  array_of_weights[22:26],
               array_of_weights[26:30], array_of_weights[30:]]
    
    pole_90S_60S = []
    mid_60S_45_S = []
    mid_45S_30_S = []
    trop_20S_10_S = []
    trop_10S_eq = []
    trop_eq_10N = []
    trop_10N_20N = []
    mid_30N_45N = []
    mid_45N_60N = []
    pole_60N_90N = []
    
    for lat_idx in range(len(lat_centroids)):
        if lat_idx <=5:
            pole_90S_60S.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>5) and (lat_idx<=9):
            mid_60S_45_S.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>9) and (lat_idx <=13):
            mid_45S_30_S.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>13) and (lat_idx <=15):
            trop_20S_10_S.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>15) and (lat_idx <=17):
            trop_10S_eq.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>17) and (lat_idx <=19):
            trop_eq_10N.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>19) and (lat_idx <=21):
            trop_10N_20N.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>21) and (lat_idx <=25):
            mid_30N_45N.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>25) and (lat_idx <=29):
            mid_45N_60N.append(diurnal_cycles_array[lat_idx])
        elif (lat_idx>29):
            pole_60N_90N.append(diurnal_cycles_array[lat_idx])
    
    
    lat_bins_arrays = [pole_90S_60S, mid_60S_45_S, mid_45S_30_S, trop_20S_10_S, trop_10S_eq, trop_eq_10N, 
                       trop_10N_20N, mid_30N_45N, mid_45N_60N, pole_60N_90N]
    
    
    weighted_average_area_cycles = []
    for j in range(len(lat_bins_arrays)):
        masked_lat_bin_values = np.ma.masked_array(lat_bins_arrays[j])
        weight_vals = weights[j]
        mean_of_lat_bin = np.ma.average(masked_lat_bin_values, axis=0, weights=weight_vals)
        weighted_average_area_cycles.append(mean_of_lat_bin)
    return np.array(weighted_average_area_cycles)