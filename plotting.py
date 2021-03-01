import numpy as np
import pandas as pd
import matplotlib
import numpy.ma as ma
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

def tropical_land_ocean(diurnal_cycles):
    land_boxes = np.ones((36, 36, 8))

    # Southern Hemisphere Land
    land_boxes[14:15, 10:14] = 0
    land_boxes[15:18, 9:14] = 0
    land_boxes[16, 9:15] = 0
    land_boxes[15:18, 18:22] = 0 
    land_boxes[14:16, 19:23] = 0
    land_boxes[14:18, 29:33] = 0
    land_boxes[16:18, 27:35] = 0
    land_boxes[14:18, 29:33] = 0

    # Northern Hemisphere Land
    land_boxes[18:20, 9:13] = 0 
    land_boxes[20:22, 7:11] = 0
    land_boxes[18:22, 16:22] = 0
    land_boxes[18, 22] = 0
    land_boxes[21, 23] = 0
    land_boxes[18:22, 27:31] = 0
    land_boxes[21, 26] = 0
    land_boxes[19:22, 25:26] = 0
    land_boxes[19:22, 22:23] = 0

    tropical_land_diurnal_cycles = ma.array(diurnal_cycles, mask=land_boxes)
    tropical_ocean_diurnal_cycles = ma.array(diurnal_cycles, mask = np.abs(land_boxes-1))
    
    tropical_land_sh = ma.mean(ma.reshape(tropical_land_diurnal_cycles[14:18], (4*36, 8)), axis=0)
    tropical_land_nh = ma.mean(ma.reshape(tropical_land_diurnal_cycles[18:22], (4*36, 8)), axis=0)
    
    tropical_ocean_sh = ma.mean(ma.reshape(tropical_ocean_diurnal_cycles[14:18], (4*36, 8)), axis=0)
    tropical_ocean_nh = ma.mean(ma.reshape(tropical_ocean_diurnal_cycles[18:22], (4*36, 8)), axis=0)
    
    return(tropical_land_sh, tropical_ocean_sh, tropical_land_nh, tropical_ocean_nh)



def season_cycle_plotter(gpsro_cycles, ERA_5_cycles, waccm6_cycels, ccsm_cycles, season_string):
    fig, axs = plt.subplots(2, 5, figsize=(18, 8))
    lt_hours = np.linspace(1.5, 22.5, 8)
    plt.suptitle(' ')
    ###########################################################################################################################
    axs[0,0].plot(lt_hours, gpsro_cycles[-1] - np.nanmean(gpsro_cycles[-1]), linewidth=3, color='black', label='GPS-RO')
    axs[0,0].plot(lt_hours, ERA_5_cycles[-1] - np.nanmean(ERA_5_cycles[-1]), linewidth=2, color='firebrick', label='ERA-5', marker="s")
    axs[0,0].plot(lt_hours, waccm6_cycels[-1] - np.nanmean(waccm6_cycels[-1]), linewidth=2, color='royalblue', label='WACCM6', marker="^")
    axs[0,0].plot(lt_hours, ccsm_cycles[-1] - np.nanmean(ccsm_cycles[-1]), linewidth=2, color='seagreen', label='CCM3', marker='o')
    fig.text(0.12, 0.97, '{season_string}'.format(season_string=season_string), fontsize=30, verticalalignment='top')

    fig.legend(bbox_to_anchor=(.2,0.87,.4,0.2), loc="lower left", mode="expand",
               borderaxespad=0, ncol=4, frameon=False, prop={'size': 15})
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,0].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                         left=True, right=True, labelsize=15)
    axs[0,0].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,0].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,0].set_ylim(-.2,.2)
    axs[0,0].set_xlim(1.5,22.5)
    axs[0,0].set_title('90$\degree$N-60$\degree$N', fontsize=14)
    axs[0,0].set_ylabel('Diurnal Anomaly ($\degree$K)', fontsize=15)

    ###########################################################################################################################
    axs[0,1].plot(lt_hours, gpsro_cycles[-2] - np.nanmean(gpsro_cycles[-2]), linewidth=3, color='black', label='GPS-RO')
    axs[0,1].plot(lt_hours, ERA_5_cycles[-2] - np.nanmean(ERA_5_cycles[-2]), linewidth=2, color='firebrick',  marker="s")
    axs[0,1].plot(lt_hours, waccm6_cycels[-2] - np.nanmean(waccm6_cycels[-2]), linewidth=2, color='royalblue',  marker="^")
    axs[0,1].plot(lt_hours, ccsm_cycles[-2] - np.nanmean(ccsm_cycles[-2]), linewidth=2, color='seagreen',  marker='o')

    axs[0,1].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,1].xaxis.set_minor_locator(MultipleLocator(3))

    axs[0,1].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                         left=True, right=True, labelsize=15)
    axs[0,1].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,1].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,1].set_ylim(-.2,.2)
    axs[0,1].set_xlim(1.5,22.5)
    axs[0,1].set_title('60$\degree$N-40$\degree$N', fontsize=14)

    ###########################################################################################################################
    axs[0,2].plot(lt_hours, gpsro_cycles[-3] - np.nanmean(gpsro_cycles[-3]), linewidth=3, color='black', label='GPS-RO')
    axs[0,2].plot(lt_hours, ERA_5_cycles[-3] - np.nanmean(ERA_5_cycles[-3]), linewidth=2, color='firebrick',  marker="s")
    axs[0,2].plot(lt_hours, waccm6_cycels[-3] - np.nanmean(waccm6_cycels[-3]), linewidth=2, color='royalblue',  marker="^")
    axs[0,2].plot(lt_hours, ccsm_cycles[-3] - np.nanmean(ccsm_cycles[-3]), linewidth=2, color='seagreen',  marker='o')

    axs[0,2].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,2].xaxis.set_minor_locator(MultipleLocator(3))

    axs[0,2].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                         left=True, right=True, labelsize=15)
    axs[0,2].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,2].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,2].set_ylim(-.2,.2)
    axs[0,2].set_xlim(1.5,22.5)
    axs[0,2].set_title('40$\degree$N-20$\degree$N', fontsize=14)

    ###########################################################################################################################
    axs[0,3].plot(lt_hours, gpsro_cycles[-4] - np.nanmean(gpsro_cycles[-4]), linewidth=3, color='black', label='GPS-RO')
    axs[0,3].plot(lt_hours, ERA_5_cycles[-4] - np.nanmean(ERA_5_cycles[-4]), linewidth=2, color='firebrick',  marker="s")
    axs[0,3].plot(lt_hours, waccm6_cycels[-4] - np.nanmean(waccm6_cycels[-4]), linewidth=2, color='royalblue',  marker="^")
    axs[0,3].plot(lt_hours, ccsm_cycles[-4] - np.nanmean(ccsm_cycles[-4]), linewidth=2, color='seagreen',  marker='o')

    axs[0,3].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,3].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,3].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,3].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,3].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,3].set_ylim(-.2,.2)
    axs[0,3].set_xlim(1.5,22.5)
    axs[0,3].set_title('20$\degree$N-10$\degree$N', fontsize=14)

    ###########################################################################################################################
    axs[0,4].plot(lt_hours, gpsro_cycles[-5] - np.nanmean(gpsro_cycles[-5]), linewidth=3, color='black', label='GPS-RO')
    axs[0,4].plot(lt_hours, ERA_5_cycles[-5] - np.nanmean(ERA_5_cycles[-5]), linewidth=2, color='firebrick',  marker="s")
    axs[0,4].plot(lt_hours, waccm6_cycels[-5] - np.nanmean(waccm6_cycels[-5]), linewidth=2, color='royalblue',  marker="^")
    axs[0,4].plot(lt_hours, ccsm_cycles[-5] - np.nanmean(ccsm_cycles[-5]), linewidth=2, color='seagreen',  marker='o')

    axs[0,4].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,4].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,4].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,4].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,4].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,4].set_ylim(-.2,.2)
    axs[0,4].set_xlim(1.5,22.5)
    axs[0,4].set_title('10$\degree$N-Equator', fontsize=14)

    ###########################################################################################################################
    axs[1,0].plot(lt_hours, gpsro_cycles[-6] - np.nanmean(gpsro_cycles[-6]), linewidth=3, color='black', label='GPS-RO')
    axs[1,0].plot(lt_hours, ERA_5_cycles[-6] - np.nanmean(ERA_5_cycles[-6]), linewidth=2, color='firebrick',  marker="s")
    axs[1,0].plot(lt_hours, waccm6_cycels[-6] - np.nanmean(waccm6_cycels[-6]), linewidth=2, color='royalblue',  marker="^")
    axs[1,0].plot(lt_hours, ccsm_cycles[-6] - np.nanmean(ccsm_cycles[-6]), linewidth=2, color='seagreen',  marker='o')

    axs[1,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,0].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,0].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,0].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,0].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,0].set_ylim(-.2,.2)
    axs[1,0].set_xlim(1.5,22.5)
    axs[1,0].set_title('10$\degree$S-Equator', fontsize=14)
    axs[1,0].set_ylabel('Diurnal Anomaly ($\degree$K)', fontsize=15)
    axs[1,0].set_xlabel('Local Time Hour', fontsize=15)



    ###########################################################################################################################
    axs[1,1].plot(lt_hours, gpsro_cycles[-7] - np.nanmean(gpsro_cycles[-7]), linewidth=3, color='black', label='GPS-RO')
    axs[1,1].plot(lt_hours, ERA_5_cycles[-7] - np.nanmean(ERA_5_cycles[-7]), linewidth=2, color='firebrick',  marker="s")
    axs[1,1].plot(lt_hours, waccm6_cycels[-7] - np.nanmean(waccm6_cycels[-7]), linewidth=2, color='royalblue',  marker="^")
    axs[1,1].plot(lt_hours, ccsm_cycles[-7] - np.nanmean(ccsm_cycles[-7]), linewidth=2, color='seagreen',  marker='o')

    axs[1,1].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,1].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,1].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,1].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,1].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,1].set_ylim(-.2,.2)
    axs[1,1].set_xlim(1.5,22.5)
    axs[1,1].set_title('10$\degree$S-20$\degree$S', fontsize=14)
    axs[1,1].set_xlabel('Local Time Hour', fontsize=15)


    ###########################################################################################################################
    axs[1,2].plot(lt_hours, gpsro_cycles[-8] - np.nanmean(gpsro_cycles[-8]), linewidth=3, color='black', label='GPS-RO')
    axs[1,2].plot(lt_hours, ERA_5_cycles[-8] - np.nanmean(ERA_5_cycles[-8]), linewidth=2, color='firebrick',  marker="s")
    axs[1,2].plot(lt_hours, waccm6_cycels[-8] - np.nanmean(waccm6_cycels[-8]), linewidth=2, color='royalblue',  marker="^")
    axs[1,2].plot(lt_hours, ccsm_cycles[-8] - np.nanmean(ccsm_cycles[-8]), linewidth=2, color='seagreen',  marker='o')

    axs[1,2].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,2].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,2].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,2].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,2].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,2].set_ylim(-.2,.2)
    axs[1,2].set_xlim(1.5,22.5)
    axs[1,2].set_title('20$\degree$S-40$\degree$S', fontsize=14)
    axs[1,2].set_xlabel('Local Time Hour', fontsize=15)


    ###########################################################################################################################
    axs[1,3].plot(lt_hours, gpsro_cycles[-9] - np.nanmean(gpsro_cycles[-9]), linewidth=3, color='black', label='GPS-RO')
    axs[1,3].plot(lt_hours, ERA_5_cycles[-9] - np.nanmean(ERA_5_cycles[-9]), linewidth=2, color='firebrick',  marker="s")
    axs[1,3].plot(lt_hours, waccm6_cycels[-9] - np.nanmean(waccm6_cycels[-9]), linewidth=2, color='royalblue',  marker="^")
    axs[1,3].plot(lt_hours, ccsm_cycles[-9] - np.nanmean(ccsm_cycles[-9]), linewidth=2, color='seagreen',  marker='o')

    axs[1,3].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,3].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,3].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,3].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,3].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,3].set_ylim(-.2,.2)
    axs[1,3].set_xlim(1.5,22.5)
    axs[1,3].set_title('40$\degree$S-60$\degree$S', fontsize=14)
    axs[1,3].set_xlabel('Local Time Hour', fontsize=15)


    ###########################################################################################################################
    axs[1,4].plot(lt_hours, gpsro_cycles[-10] - np.nanmean(gpsro_cycles[-10]), linewidth=3, color='black', label='GPS-RO')
    axs[1,4].plot(lt_hours, ERA_5_cycles[-10] - np.nanmean(ERA_5_cycles[-10]), linewidth=2, color='firebrick',  marker="s")
    axs[1,4].plot(lt_hours, waccm6_cycels[-10] - np.nanmean(waccm6_cycels[-10]), linewidth=2, color='royalblue',  marker="^")
    axs[1,4].plot(lt_hours, ccsm_cycles[-10] - np.nanmean(ccsm_cycles[-10]), linewidth=2, color='seagreen',  marker='o')

    axs[1,4].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,4].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,4].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,4].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,4].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,4].set_ylim(-.2,.2)
    axs[1,4].set_xlim(1.5,22.5)
    axs[1,4].set_title('60$\degree$S-90$\degree$S', fontsize=14)
    axs[1,4].set_xlabel('Local Time Hour', fontsize=15)

    return(fig)

def seasonal_colocation_plotter(ERA_5_cycles, ERA_5_colocation_cycles, ERA_5_colocations_no_removal,
                                season_string):
    fig, axs = plt.subplots(2, 5, figsize=(18, 8))
    lt_hours = np.linspace(1.5, 22.5, 8)

    plt.suptitle(' ')
    ###########################################################################################################################
    axs[0,0].plot(lt_hours, ERA_5_cycles[-1] - np.nanmean(ERA_5_cycles[-1]), linewidth=4, color='black', label='ERA-5')
    axs[0,0].plot(lt_hours, ERA_5_colocation_cycles[-1] - np.nanmean(ERA_5_colocation_cycles[-1]), linewidth=2, color='red', label='Colocations \nmean removal')
    axs[0,0].plot(lt_hours, ERA_5_colocations_no_removal[-1] - np.nanmean(ERA_5_colocations_no_removal[-1]), linewidth=2, color='blue', label='Colocations \nw/o mean removal')
    fig.text(0.12, 0.97, '{season_string}'.format(season_string=season_string), fontsize=30, verticalalignment='top')
    fig.legend(bbox_to_anchor=(.2,0.87,.4,0.2), loc="lower left", mode="expand",
               borderaxespad=0, ncol=4, frameon=False, prop={'size': 15})
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,0].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                         left=True, right=True, labelsize=15)
    axs[0,0].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,0].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,0].set_ylim(-.2,.2)
    axs[0,0].set_xlim(1.5,22.5)
    axs[0,0].set_title('90$\degree$N-60$\degree$N', fontsize=14)
    axs[0,0].set_ylabel('Diurnal Anomaly ($\degree$K)', fontsize=15)

    ###########################################################################################################################
    axs[0,1].plot(lt_hours, ERA_5_cycles[-2] - np.nanmean(ERA_5_cycles[-2]), linewidth=4, color='black')
    axs[0,1].plot(lt_hours, ERA_5_colocation_cycles[-2] - np.nanmean(ERA_5_colocation_cycles[-2]), linewidth=2, color='red')
    axs[0,1].plot(lt_hours, ERA_5_colocations_no_removal[-2] - np.nanmean(ERA_5_colocations_no_removal[-2]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[0,1].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,1].xaxis.set_minor_locator(MultipleLocator(3))

    axs[0,1].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                         left=True, right=True, labelsize=15)
    axs[0,1].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,1].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,1].set_ylim(-.2,.2)
    axs[0,1].set_xlim(1.5,22.5)
    axs[0,1].set_title('60$\degree$N-40$\degree$N', fontsize=14)

    ###########################################################################################################################
    axs[0,2].plot(lt_hours, ERA_5_cycles[-3] - np.nanmean(ERA_5_cycles[-3]), linewidth=4, color='black')
    axs[0,2].plot(lt_hours, ERA_5_colocation_cycles[-3] - np.nanmean(ERA_5_colocation_cycles[-3]), linewidth=2, color='red')
    axs[0,2].plot(lt_hours, ERA_5_colocations_no_removal[-3] - np.nanmean(ERA_5_colocations_no_removal[-3]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[0,2].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,2].xaxis.set_minor_locator(MultipleLocator(3))

    axs[0,2].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                         left=True, right=True, labelsize=15)
    axs[0,2].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,2].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,2].set_ylim(-.2,.2)
    axs[0,2].set_xlim(1.5,22.5)
    axs[0,2].set_title('40$\degree$N-20$\degree$N', fontsize=14)

    ###########################################################################################################################
    axs[0,3].plot(lt_hours, ERA_5_cycles[-4] - np.nanmean(ERA_5_cycles[-4]), linewidth=4, color='black')
    axs[0,3].plot(lt_hours, ERA_5_colocation_cycles[-4] - np.nanmean(ERA_5_colocation_cycles[-4]), linewidth=2, color='red')
    axs[0,3].plot(lt_hours, ERA_5_colocations_no_removal[-4] - np.nanmean(ERA_5_colocations_no_removal[-4]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[0,3].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,3].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,3].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,3].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,3].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,3].set_ylim(-.2,.2)
    axs[0,3].set_xlim(1.5,22.5)
    axs[0,3].set_title('20$\degree$N-10$\degree$N', fontsize=14)

    ###########################################################################################################################
    axs[0,4].plot(lt_hours, ERA_5_cycles[-5] - np.nanmean(ERA_5_cycles[-5]), linewidth=4, color='black')
    axs[0,4].plot(lt_hours, ERA_5_colocation_cycles[-5] - np.nanmean(ERA_5_colocation_cycles[-5]), linewidth=2, color='red')
    axs[0,4].plot(lt_hours, ERA_5_colocations_no_removal[-5] - np.nanmean(ERA_5_colocations_no_removal[-5]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[0,4].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,4].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,4].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,4].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,4].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,4].set_ylim(-.2,.2)
    axs[0,4].set_xlim(1.5,22.5)
    axs[0,4].set_title('10$\degree$N-Equator', fontsize=14)

    ###########################################################################################################################
    axs[1,0].plot(lt_hours, ERA_5_cycles[-6] - np.nanmean(ERA_5_cycles[-6]), linewidth=4, color='black')
    axs[1,0].plot(lt_hours, ERA_5_colocation_cycles[-6] - np.nanmean(ERA_5_colocation_cycles[-6]), linewidth=2, color='red')
    axs[1,0].plot(lt_hours, ERA_5_colocations_no_removal[-6] - np.nanmean(ERA_5_colocations_no_removal[-6]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[1,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,0].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,0].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,0].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,0].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,0].set_ylim(-.2,.2)
    axs[1,0].set_xlim(1.5,22.5)
    axs[1,0].set_title('Equator-10$\degree$S', fontsize=14)
    axs[1,0].set_ylabel('Diurnal Anomaly ($\degree$K)', fontsize=15)
    axs[1,0].set_xlabel('Local Time Hour', fontsize=15)



    ###########################################################################################################################
    axs[1,1].plot(lt_hours, ERA_5_cycles[-7] - np.nanmean(ERA_5_cycles[-7]), linewidth=4, color='black')
    axs[1,1].plot(lt_hours, ERA_5_colocation_cycles[-7] - np.nanmean(ERA_5_colocation_cycles[-7]), linewidth=2, color='red')
    axs[1,1].plot(lt_hours, ERA_5_colocations_no_removal[-7] - np.nanmean(ERA_5_colocations_no_removal[-7]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[1,1].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,1].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,1].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,1].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,1].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,1].set_ylim(-.2,.2)
    axs[1,1].set_xlim(1.5,22.5)
    axs[1,1].set_title('10$\degree$S-20$\degree$S', fontsize=14)
    axs[1,1].set_xlabel('Local Time Hour', fontsize=15)


    ###########################################################################################################################
    axs[1,2].plot(lt_hours, ERA_5_cycles[-8] - np.nanmean(ERA_5_cycles[-8]), linewidth=4, color='black')
    axs[1,2].plot(lt_hours, ERA_5_colocation_cycles[-8] - np.nanmean(ERA_5_colocation_cycles[-8]), linewidth=2, color='red')
    axs[1,2].plot(lt_hours, ERA_5_colocations_no_removal[-8] - np.nanmean(ERA_5_colocations_no_removal[-8]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[1,2].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,2].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,2].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,2].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,2].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,2].set_ylim(-.2,.2)
    axs[1,2].set_xlim(1.5,22.5)
    axs[1,2].set_title('20$\degree$S-40$\degree$S', fontsize=14)
    axs[1,2].set_xlabel('Local Time Hour', fontsize=15)


    ###########################################################################################################################
    axs[1,3].plot(lt_hours, ERA_5_cycles[-9] - np.nanmean(ERA_5_cycles[-9]), linewidth=4, color='black')
    axs[1,3].plot(lt_hours, ERA_5_colocation_cycles[-9] - np.nanmean(ERA_5_colocation_cycles[-9]), linewidth=2, color='red')
    axs[1,3].plot(lt_hours, ERA_5_colocations_no_removal[-9] - np.nanmean(ERA_5_colocations_no_removal[-9]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[1,3].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,3].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,3].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,3].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,3].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,3].set_ylim(-.2,.2)
    axs[1,3].set_xlim(1.5,22.5)
    axs[1,3].set_title('40$\degree$S-60$\degree$S', fontsize=14)
    axs[1,3].set_xlabel('Local Time Hour', fontsize=15)


    ###########################################################################################################################
    axs[1,4].plot(lt_hours, ERA_5_cycles[-10] - np.nanmean(ERA_5_cycles[-10]), linewidth=4, color='black')
    axs[1,4].plot(lt_hours, ERA_5_colocation_cycles[-10] - np.nanmean(ERA_5_colocation_cycles[-10]), linewidth=2, color='red')
    axs[1,4].plot(lt_hours, ERA_5_colocations_no_removal[-10] - np.nanmean(ERA_5_colocations_no_removal[-10]), linewidth=2, color='blue', label='Colocations no mean removal')
    axs[1,4].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,4].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,4].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,4].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,4].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,4].set_ylim(-.2,.2)
    axs[1,4].set_xlim(1.5,22.5)
    axs[1,4].set_title('60$\degree$S-90$\degree$S', fontsize=14)
    axs[1,4].set_xlabel('Local Time Hour', fontsize=15)
    
    return(fig)

def land_ocean_plotting(nh_land_gpsro_djf, nh_land_era5_djf, nh_land_waccm6_djf, nh_land_ccm3_djf, 
                        nh_ocean_gpsro_djf, nh_ocean_era5_djf, nh_ocean_waccm6_djf, nh_ocean_ccm3_djf,
                        nh_land_gpsro_mam, nh_land_era5_mam, nh_land_waccm6_mam, nh_land_ccm3_mam, 
                        nh_ocean_gpsro_mam, nh_ocean_era5_mam, nh_ocean_waccm6_mam, nh_ocean_ccm3_mam,
                        nh_land_gpsro_jja, nh_land_era5_jja, nh_land_waccm6_jja, nh_land_ccm3_jja, 
                        nh_ocean_gpsro_jja, nh_ocean_era5_jja, nh_ocean_waccm6_jja, nh_ocean_ccm3_jja,
                        nh_land_gpsro_son, nh_land_era5_son, nh_land_waccm6_son, nh_land_ccm3_son, 
                        nh_ocean_gpsro_son, nh_ocean_era5_son, nh_ocean_waccm6_son, nh_ocean_ccm3_son,
                        sh_land_gpsro_djf, sh_land_era5_djf, sh_land_waccm6_djf, sh_land_ccm3_djf, 
                        sh_ocean_gpsro_djf, sh_ocean_era5_djf, sh_ocean_waccm6_djf, sh_ocean_ccm3_djf,
                        sh_land_gpsro_mam, sh_land_era5_mam, sh_land_waccm6_mam, sh_land_ccm3_mam, 
                        sh_ocean_gpsro_mam, sh_ocean_era5_mam, sh_ocean_waccm6_mam, sh_ocean_ccm3_mam,
                        sh_land_gpsro_jja, sh_land_era5_jja, sh_land_waccm6_jja, sh_land_ccm3_jja, 
                        sh_ocean_gpsro_jja, sh_ocean_era5_jja, sh_ocean_waccm6_jja, sh_ocean_ccm3_jja,
                        sh_land_gpsro_son, sh_land_era5_son, sh_land_waccm6_son, sh_land_ccm3_son, 
                        sh_ocean_gpsro_son, sh_ocean_era5_son, sh_ocean_waccm6_son, sh_ocean_ccm3_son
                        ):
    fig, axs = plt.subplots(2, 8, figsize=(18, 7))
    lt_hours = np.linspace(1.5, 22.5, 8)
    
    axs[0,0].plot(lt_hours,nh_land_gpsro_djf, label='GPS-RO', color='black', linewidth=3)
    axs[0,0].plot(lt_hours,nh_land_era5_djf, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[0,0].plot(lt_hours,nh_land_waccm6_djf, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[0,0].plot(lt_hours,nh_land_ccm3_djf, label='CCM3', color='seagreen', marker='o', linewidth=2)
    fig.text(0.12, 0.97, 'Land Ocean Contrast', fontsize=30, verticalalignment='top')

    fig.legend(bbox_to_anchor=(.4,0.87,.4,0.2), loc="lower left", mode="expand",
               borderaxespad=0, ncol=4, frameon=False, prop={'size': 15})
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,0].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,0].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,0].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,0].set_ylim(-.3, .3)
    axs[0,0].set_xlim(1.5,22.5)
    axs[0,0].set_ylabel('Equator-20$\degree$N \nDiurnal Amonamly ($\degree$K)', fontsize=14)
    axs[0,0].set_title('DJF Land', fontsize=14)

    ###############################################################################################
    axs[0,1].plot(lt_hours,nh_ocean_gpsro_djf, label='GPS-RO', color='black', linewidth=3)
    axs[0,1].plot(lt_hours,nh_ocean_era5_djf, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[0,1].plot(lt_hours,nh_ocean_waccm6_djf, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[0,1].plot(lt_hours,nh_ocean_ccm3_djf, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[0,1].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,1].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[0,1].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,1].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,1].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,1].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,1].set_ylim(-.3, .3)
    axs[0,1].set_xlim(1.5,22.5)
    axs[0,1].set_yticklabels([])
    axs[0,1].set_title('DJF Ocean', fontsize=14)
    ##############################################################################################
    axs[0,2].plot(lt_hours,nh_land_gpsro_mam, label='GPS-RO', color='black', linewidth=3)
    axs[0,2].plot(lt_hours,nh_land_era5_mam, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[0,2].plot(lt_hours,nh_land_waccm6_mam, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[0,2].plot(lt_hours,nh_land_ccm3_mam, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[0,2].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[0,2].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,2].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,2].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,2].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,2].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,2].set_ylim(-.3, .3)
    axs[0,2].set_xlim(1.5,22.5)
    axs[0,2].set_yticklabels([])
    axs[0,2].set_title('MAM Land', fontsize=14)
    ###############################################################################################
    axs[0,3].plot(lt_hours,nh_ocean_gpsro_mam, label='GPS-RO', color='black', linewidth=3)
    axs[0,3].plot(lt_hours,nh_ocean_era5_mam, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[0,3].plot(lt_hours,nh_ocean_waccm6_mam, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[0,3].plot(lt_hours,nh_ocean_ccm3_mam, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[0,3].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,3].yaxis.set_major_locator(MultipleLocator(0.15))

    axs[0,3].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,3].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,3].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,3].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,3].set_ylim(-.3, .3)
    axs[0,3].set_xlim(1.5,22.5)
    axs[0,3].set_yticklabels([])
    axs[0,3].set_title('MAM Ocean', fontsize=14)
    ###############################################################################################
    axs[0,4].plot(lt_hours,nh_land_gpsro_jja, label='GPS-RO', color='black', linewidth=3)
    axs[0,4].plot(lt_hours,nh_land_era5_jja, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[0,4].plot(lt_hours,nh_land_waccm6_jja, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[0,4].plot(lt_hours,nh_land_ccm3_jja, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[0,4].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[0,4].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,4].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,4].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,4].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,4].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,4].set_ylim(-.3, .3)
    axs[0,4].set_xlim(1.5,22.5)
    axs[0,4].set_yticklabels([])
    #axs[0,4].set_xticklabels([])
    axs[0,4].set_title('JJA land', fontsize=14)
    ###############################################################################################
    axs[0,5].plot(lt_hours,nh_ocean_gpsro_jja, label='GPS-RO', color='black', linewidth=3)
    axs[0,5].plot(lt_hours,nh_ocean_era5_jja, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[0,5].plot(lt_hours,nh_ocean_waccm6_jja, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[0,5].plot(lt_hours,nh_ocean_ccm3_jja, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[0,5].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[0,5].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,5].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,5].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,5].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,5].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,5].set_ylim(-.3, .3)
    axs[0,5].set_xlim(1.5,22.5)
    axs[0,5].set_yticklabels([])
    #axs[0,5].set_xticklabels([])
    axs[0,5].set_title('JJA Ocean', fontsize=14)

    ###############################################################################################
    axs[0,6].plot(lt_hours,nh_land_gpsro_son, label='GPS-RO', color='black', linewidth=3)
    axs[0,6].plot(lt_hours,nh_land_era5_son, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[0,6].plot(lt_hours,nh_land_waccm6_son, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[0,6].plot(lt_hours,nh_land_ccm3_son, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[0,6].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[0,6].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,6].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,6].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,6].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,6].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,6].set_ylim(-.3, .3)
    axs[0,6].set_xlim(1.5,22.5)
    axs[0,6].set_yticklabels([])
    axs[0,6].set_title('SON Land', fontsize=14)
    
    ###############################################################################################
    axs[0,7].plot(lt_hours,nh_ocean_gpsro_son, label='GPS-RO', color='black', linewidth=3)
    axs[0,7].plot(lt_hours,nh_ocean_era5_son, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[0,7].plot(lt_hours,nh_ocean_waccm6_son, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[0,7].plot(lt_hours,nh_ocean_ccm3_son, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[0,7].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[0,7].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,7].xaxis.set_minor_locator(MultipleLocator(3))
    axs[0,7].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[0,7].set_xticks(ticks=[3, 9, 15, 21])
    axs[0,7].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[0,7].set_ylim(-.3, .3)
    axs[0,7].set_xlim(1.5,22.5)
    axs[0,7].set_yticklabels([])
    axs[0,7].set_title('SON Ocean', fontsize=14)

    ###############################################################################################
    axs[1,0].plot(lt_hours,sh_land_gpsro_djf, label='GPS-RO', color='black', linewidth=3)
    axs[1,0].plot(lt_hours,sh_land_era5_djf, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[1,0].plot(lt_hours,sh_land_waccm6_djf, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[1,0].plot(lt_hours,sh_land_ccm3_djf, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[1,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,0].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,0].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,0].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[1,0].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,0].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,0].set_ylim(-.3, .3)
    axs[1,0].set_xlim(1.5,22.5)
    axs[1,0].set_ylabel('20$\degree$S-Equator \nDiurnal Amonamly ($\degree$K)', fontsize=14)
    axs[1,0].set_xlabel('LTH', fontsize=15)

    ###############################################################################################
    axs[1,1].plot(lt_hours,sh_ocean_gpsro_djf, label='GPS-RO', color='black', linewidth=3)
    axs[1,1].plot(lt_hours,sh_ocean_era5_djf, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[1,1].plot(lt_hours,sh_ocean_waccm6_djf, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[1,1].plot(lt_hours,sh_ocean_ccm3_djf, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[1,1].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[1,1].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,1].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,1].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,1].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,1].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,1].set_ylim(-.3, .3)
    axs[1,1].set_xlim(1.5,22.5)
    axs[1,1].set_yticklabels([])
    axs[1,1].set_xlabel('LTH', fontsize=15)

    ###############################################################################################
    axs[1,2].plot(lt_hours,sh_land_gpsro_mam, label='GPS-RO', color='black', linewidth=3)
    axs[1,2].plot(lt_hours,sh_land_era5_mam, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[1,2].plot(lt_hours,sh_land_waccm6_mam, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[1,2].plot(lt_hours,sh_land_ccm3_mam, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[1,2].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,2].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,2].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[1,2].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,2].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,2].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,2].set_ylim(-.3, .3)
    axs[1,2].set_xlim(1.5,22.5)
    axs[1,2].set_yticklabels([])
    axs[1,2].set_xlabel('LTH', fontsize=15)

    ###############################################################################################
    axs[1,3].plot(lt_hours,sh_ocean_gpsro_mam, label='GPS-RO', color='black', linewidth=3)
    axs[1,3].plot(lt_hours,sh_ocean_era5_mam, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[1,3].plot(lt_hours,sh_ocean_waccm6_mam, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[1,3].plot(lt_hours,sh_ocean_ccm3_mam, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[1,3].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[1,3].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,3].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,3].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,3].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,3].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,3].set_ylim(-.3, .3)
    axs[1,3].set_xlim(1.5,22.5)
    axs[1,3].set_yticklabels([])
    axs[1,3].set_xlabel('LTH', fontsize=15)

    ###############################################################################################
    axs[1,4].plot(lt_hours,sh_land_gpsro_jja, label='GPS-RO', color='black', linewidth=3)
    axs[1,4].plot(lt_hours,sh_land_era5_jja, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[1,4].plot(lt_hours,sh_land_waccm6_jja, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[1,4].plot(lt_hours,sh_land_ccm3_jja, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[1,4].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,4].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[1,4].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,4].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,4].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,4].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,4].set_ylim(-.3, .3)
    axs[1,4].set_xlim(1.5,22.5)
    axs[1,4].set_yticklabels([])
    axs[1,4].set_xlabel('LTH', fontsize=15)

    ###############################################################################################
    axs[1,5].plot(lt_hours,sh_ocean_gpsro_jja, label='GPS-RO', color='black', linewidth=3)
    axs[1,5].plot(lt_hours,sh_ocean_era5_jja, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[1,5].plot(lt_hours,sh_ocean_waccm6_jja, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[1,5].plot(lt_hours,sh_ocean_ccm3_jja, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[1,5].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[1,5].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,5].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,5].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                      left=True, right=True, labelsize=15)
    axs[1,5].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,5].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,5].set_ylim(-.3, .3)
    axs[1,5].set_xlim(1.5,22.5)
    axs[1,5].set_yticklabels([])
    axs[1,5].set_xlabel('LTH', fontsize=15)
    
    ###############################################################################################
    axs[1,6].plot(lt_hours,sh_land_gpsro_son, label='GPS-RO', color='black', linewidth=3)
    axs[1,6].plot(lt_hours,sh_land_era5_son, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[1,6].plot(lt_hours,sh_land_waccm6_son, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[1,6].plot(lt_hours,sh_land_ccm3_son, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[1,6].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[1,6].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,6].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,6].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True, 
                         left=True, right=True, labelsize=15)
    axs[1,6].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,6].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,6].set_ylim(-.3, .3)
    axs[1,6].set_xlim(1.5,22.5)
    axs[1,6].set_yticklabels([])
    axs[1,6].set_xlabel('LTH', fontsize=15)
    
    ###############################################################################################
    axs[1,7].plot(lt_hours,sh_ocean_gpsro_son, label='GPS-RO', color='black', linewidth=3)
    axs[1,7].plot(lt_hours,sh_ocean_era5_son, label='ERA-5', color='firebrick', marker='s', linewidth=2)
    axs[1,7].plot(lt_hours,sh_ocean_waccm6_son, label='WACCM6', color='dodgerblue', marker='^', linewidth=2)
    axs[1,7].plot(lt_hours,sh_ocean_ccm3_son, label='CCM3', color='seagreen', marker='o', linewidth=2)
    axs[1,7].yaxis.set_major_locator(MultipleLocator(0.15))
    axs[1,7].yaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1,7].xaxis.set_minor_locator(MultipleLocator(3))
    axs[1,7].tick_params(direction='in', length=8, width=2, colors='black', bottom=True, top=True,
                         left=True, right=True, labelsize=15)
    axs[1,7].set_xticks(ticks=[3, 9, 15, 21])
    axs[1,7].tick_params(which='minor', direction='in', length=4, width=1.5, left=True, right=True)
    axs[1,7].set_ylim(-.3, .3)
    axs[1,7].set_xlim(1.5,22.5)
    axs[1,7].set_yticklabels([])
    axs[1,7].set_xlabel('LTH', fontsize=15)
