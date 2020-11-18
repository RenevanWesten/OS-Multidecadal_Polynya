#Program plots the vertical OHC and ice thickness/fraction over the Polynya Region

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory 	= '../../../Data/CESM/'

def TrendRemover(time, data, trend_type):
	"""Removes trend of choice"""
	
	rank = polyfit(time, data, trend_type)
	fitting = 0.0 
		
	for rank_i in range(len(rank)):
			
		fitting += rank[rank_i] * (time**(len(rank) - 1 - rank_i))

	data -= fitting
	
	return data

def MovingAverage(time, data, moving_average):
	"""Determines moving average of time series"""
	
	#Empty array's for the smoothed time series
	time_2     = ma.masked_all(len(time) - moving_average + 1)
	data_2 	   = ma.masked_all(len(time_2))
	
	#Determine the so-called middle index where the moving average is determined
	selection  = (moving_average - 1) / 2		

	for time_i in range(selection, selection + len(time_2)):
		#Take moving average
		data_2[time_i - selection]  	= np.mean(data[time_i-selection:time_i-selection + moving_average], axis = 0)
		time_2[time_i - selection] 	= time[time_i]
	
	return time_2, data_2

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------
trend_type 	= 2	#Remove trend (1 = linear, 2 = quadratic)
moving_average 	= 60 	#Running mean

#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/Polynya_OHC_depth_sum_0-100m.nc', 'r')
	
time	= HEAT_data.variables['time'][:]     	
OHC_1	= HEAT_data.variables['OHC'][:] 		

HEAT_data.close()

#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/Polynya_OHC_depth_sum_200-1000m.nc', 'r')
	
time	= HEAT_data.variables['time'][:]     	
OHC_2	= HEAT_data.variables['OHC'][:] 		

HEAT_data.close()

#-----------------------------------------------------------------------------------------

ICE_data 	= netcdf.Dataset(directory+'Ice/Polynya_ice_fraction_thickness.nc', 'r')

time	 	= ICE_data.variables['time'][:]     	
thickness	= ICE_data.variables['hi'][:] * 100.0
fraction	= ICE_data.variables['aice'][:]

ICE_data.close()

#-----------------------------------------------------------------------------------------

time_year		= np.zeros(len(time) / 12)
ice_thickness_year	= np.zeros(len(time) / 12)
ice_fraction_year	= np.zeros(len(time) / 12)

for time_i in range(len(time) / 12):
	#Take the September ice fraction and thickness
	time_year[time_i] 		= np.mean(time[time_i * 12: (time_i + 1) * 12])
	ice_thickness_year[time_i]	= thickness[time_i * 12 + 8]
	ice_fraction_year[time_i]	= fraction[time_i * 12 + 8]

if trend_type > 0:
	#Remove trend
	OHC_1 	= TrendRemover(time, OHC_1, trend_type)
	OHC_2	= TrendRemover(time, OHC_2, trend_type)

if moving_average > 1:
	#Apply moving average
	OHC_1				= MovingAverage(time, OHC_1, moving_average)[1]
	time, OHC_2			= MovingAverage(time, OHC_2, moving_average)

#-----------------------------------------------------------------------------------------

fig, ax1 	= plt.subplots(figsize = (8, 6.8))

ax1.fill_between(np.arange(datetime.datetime(158, 1, 1).toordinal(), datetime.datetime(160, 1, 1).toordinal()), -1, 1, alpha=0.25, edgecolor='orange', facecolor='orange')
ax1.fill_between(np.arange(datetime.datetime(178, 1, 1).toordinal(), datetime.datetime(183, 1, 1).toordinal()), -1, 1, alpha=0.25, edgecolor='orange', facecolor='orange')
ax1.fill_between(np.arange(datetime.datetime(205, 1, 1).toordinal(), datetime.datetime(210, 1, 1).toordinal()), -1, 1, alpha=0.25, edgecolor='orange', facecolor='orange')
ax1.fill_between(np.arange(datetime.datetime(231, 1, 1).toordinal(), datetime.datetime(238, 1, 1).toordinal()), -1, 1, alpha=0.25, edgecolor='orange', facecolor='orange')

OHC_1_graph	= ax1.plot_date(time, OHC_1 * 5.0, '-r', linewidth = 2.0, label = 'OHC 100 m (5x)')
OHC_2_graph	= ax1.plot_date(time, OHC_2, '-k', linewidth = 2.0, label = 'OHC 200 - 1000 m')

ax1.set_xlabel('Model year')
ax1.set_ylabel('Ocean heat content anomaly (ZJ)')
ax1.set_ylim(-0.3, 0.3)
ax1.set_xlim(datetime.datetime(152, 6, 15).toordinal(), datetime.datetime(248, 6, 15).toordinal())
ax1.grid()

ax12 = ax1.twinx()
#Ice_thickness_graph = ax12.plot_date(time_year, ice_thickness_year, '-b', linewidth = 2.0, label = 'Thickness')
Ice_fraction_graph = ax12.plot_date(time_year, ice_fraction_year * 100.0, '-b', linewidth = 2.0, label = 'Sea-ice fraction')
ax12.set_ylim(0, 100)
ax12.set_ylabel('Sea-ice fraction ($\%$)', color ='b')
ax12.set_xlim(datetime.datetime(152, 6, 15).toordinal(), datetime.datetime(248, 6, 15).toordinal())

for tl in ax12.get_yticklabels():
    tl.set_color('b')

graphs	      = OHC_1_graph + OHC_2_graph + Ice_fraction_graph
legend_labels = [l.get_label() for l in graphs]
legend_2      = legend(graphs, legend_labels, loc='lower left', #bbox_to_anchor=(0.73, 0.142),
		  ncol=1, fancybox=True, numpoints = 1)

ax1.set_xticks([datetime.datetime(160, 1, 1).toordinal(), 
		datetime.datetime(170, 1, 1).toordinal(),
		datetime.datetime(180, 1, 1).toordinal(),
		datetime.datetime(190, 1, 1).toordinal(), 
		datetime.datetime(200, 1, 1).toordinal(),
		datetime.datetime(210, 1, 1).toordinal(), 
		datetime.datetime(220, 1, 1).toordinal(),
		datetime.datetime(230, 1, 1).toordinal(),
		datetime.datetime(240, 1, 1).toordinal()])

ax1.set_title('a) OHC and sea-ice fraction')

show()
