#Program plots the heat advection at the eastern boundary of the Polynya Region

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory 		= '../../../Data/CESM/'

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

moving_average 	= 60 	#Running mean

#-----------------------------------------------------------------------------------------

HEAT_data 	= netcdf.Dataset(directory+'Ocean/Polynya_HEAT_advection.nc', 'r')

#Writing data to correct variable	
time		= HEAT_data.variables['time'][:]     	  	
depth		= HEAT_data.variables['depth'][:]     	  	
u_heat		= HEAT_data.variables['UHEAT_EAST'][:] 		

HEAT_data.close()


#-----------------------------------------------------------------------------------------
TEMP_data 	= netcdf.Dataset(directory+'Ocean/Polynya_mixed_layer_depth.nc', 'r')

time 		= TEMP_data.variables['time'][:]     	
mixed		= TEMP_data.variables['XMXL'][:] 

TEMP_data.close()
#-----------------------------------------------------------------------------------------

time_year	= np.zeros(len(time) / 12)
mixed_year	= np.zeros(len(time) / 12)

for time_i in range(len(time) / 12):
	#Take the maximum for each year
	time_year[time_i] 	= np.mean(time[time_i * 12: (time_i + 1) * 12])
	mixed_year[time_i]	= mixed[time_i * 12: (time_i + 1) * 12].max()


#-----------------------------------------------------------------------------------------

if moving_average > 1:
	
	time_2, u_heat_2	= MovingAverage(time, u_heat[:, 0], moving_average)
	u_heat_2		= ma.masked_all((len(time_2), len(depth)))

	for depth_i in range(len(depth)):
		#Take moving average over each depth level
		time_2, u_heat_2[:, depth_i]	= MovingAverage(time, u_heat[:, depth_i], moving_average)
	
	#Save the smoothed files
	time	= time_2	
	u_heat	= u_heat_2		

#-----------------------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 2], 'hspace': 0.0})

x, y 	= np.meshgrid(time, depth)

levels 	= np.arange(-7, 7.1, 0.5)
CS 	= ax1.contourf(x, y, -u_heat.transpose() / 10**12.0, levels, extend = 'both', cmap = 'RdBu_r')


fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.10, 0.025, 0.8])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-7, 7.1, 1))
cbar.set_ylabel('Advective heat (TW)')

ax1.set_xlabel('Model year')
#ax1.set_ylabel("Depth (m)")
ax1.xaxis_date()
ax1.set_ylim(200, 0)


ax1.plot_date(time_year, mixed_year, '-k', linewidth = 3.0)

ax1.set_xticks([datetime.datetime(160, 1, 1).toordinal(), 
		datetime.datetime(170, 1, 1).toordinal(),
		datetime.datetime(180, 1, 1).toordinal(),
		datetime.datetime(190, 1, 1).toordinal(), 
		datetime.datetime(200, 1, 1).toordinal(),
		datetime.datetime(210, 1, 1).toordinal(), 
		datetime.datetime(220, 1, 1).toordinal(),
		datetime.datetime(230, 1, 1).toordinal(),
		datetime.datetime(240, 1, 1).toordinal()])

ax1.set_xticklabels([])
ax1.set_xlim(datetime.datetime(152, 6, 15).toordinal(), datetime.datetime(248, 6, 15).toordinal())


x, y 	= np.meshgrid(time, depth)

CS 	= ax2.contourf(x, y, -u_heat.transpose() / 10**12.0, levels, extend = 'both', cmap = 'RdBu_r')

ax2.set_xlabel('Model year')
ax2.set_ylabel("Depth (m)")
ax2.xaxis_date()
ax2.set_ylim(2000, 200)


ax2.plot_date(time_year, mixed_year, '-k', linewidth = 3.0)

ax2.set_xticks([datetime.datetime(160, 1, 1).toordinal(), 
		datetime.datetime(170, 1, 1).toordinal(),
		datetime.datetime(180, 1, 1).toordinal(),
		datetime.datetime(190, 1, 1).toordinal(), 
		datetime.datetime(200, 1, 1).toordinal(),
		datetime.datetime(210, 1, 1).toordinal(), 
		datetime.datetime(220, 1, 1).toordinal(),
		datetime.datetime(230, 1, 1).toordinal(),
		datetime.datetime(240, 1, 1).toordinal()])

ax2.set_yticks([500, 1000, 1500, 2000])

ax2.set_xlim(datetime.datetime(152, 6, 15).toordinal(), datetime.datetime(248, 6, 15).toordinal())

ax1.set_title('a) Advective heat and mixed layer')

show()

