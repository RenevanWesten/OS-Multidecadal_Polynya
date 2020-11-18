#Program plots the SAM index and the Weddell Gyre strength

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats

#Making pathway to folder with all data
directory 	= '../../../Data/POP/'

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
depth_min	= 0
depth_max	= 6000
lon_Weddell	= 330

trend_type 	= 2	#Remove trend (1 = linear, 2 = quadratic)
moving_average 	= 60 	#Running mean
#-----------------------------------------------------------------------------------------

TEMP_data 	= netcdf.Dataset(directory+'Atmosphere/SAM_index.nc', 'r')

time 		= TEMP_data.variables['time'][:]	
SAM		= TEMP_data.variables['SAM'][:]

TEMP_data.close()


HEAT_data 	= netcdf.Dataset(directory+'/Ocean/Weddell_Gyre_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'_m_lon_'+str(lon_Weddell)+'E.nc', 'r')

time		= HEAT_data.variables['time'][:]     
transport	= -HEAT_data.variables['Transport'][:] 

HEAT_data.close()


#-----------------------------------------------------------------------------------------

if trend_type > 0:
	
	SAM 		= TrendRemover(time, SAM, trend_type)
	transport 	= TrendRemover(time, transport, trend_type)

if moving_average > 1:

	time_2, SAM	= MovingAverage(time, SAM, moving_average)
	time, transport	= MovingAverage(time, transport, moving_average)

#-----------------------------------------------------------------------------------------
#----------------------------------------PLOTTING-----------------------------------------
#-----------------------------------------------------------------------------------------	

fig, ax1 = plt.subplots()
ax1.grid()

SAM_graph = ax1.plot_date(time, SAM, '-k', linewidth = 2.0, label = 'SAM index')
ax1.set_xlabel('Model year')
ax1.set_ylabel('Southern Annular Mode index', color ='k')
ax1.set_ylim(-0.35, 0.35)


ax2 		= ax1.twinx()

Transport_graph = ax2.plot_date(time, transport * 0.5, '-r', linewidth = 2.0, label = 'Transport (0.5x)')

ax2.set_ylim(-4.7, 4.7)
ax2.set_ylabel('Westward transport anomaly (Sv)', color ='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')


ax1.set_xticks([datetime.datetime(180, 1, 1).toordinal(), 
		datetime.datetime(190, 1, 1).toordinal(),
		datetime.datetime(200, 1, 1).toordinal(),
		datetime.datetime(210, 1, 1).toordinal(), 
		datetime.datetime(220, 1, 1).toordinal(),
		datetime.datetime(230, 1, 1).toordinal(), 
		datetime.datetime(240, 1, 1).toordinal(),
		datetime.datetime(250, 1, 1).toordinal(),
		datetime.datetime(260, 1, 1).toordinal(),
		datetime.datetime(270, 1, 1).toordinal()])

graphs	      = SAM_graph + Transport_graph
legend_labels = [l.get_label() for l in graphs]

ax2.legend(graphs, legend_labels, loc='lower left', ncol=1, fancybox=True, shadow=False, numpoints = 1)

ax1.set_title('e) POP, SAM index and volume transport')

show()


