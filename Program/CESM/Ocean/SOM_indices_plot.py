#Program plots the PC index and the SOM indices

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
	
	rank 	= polyfit(time, data, trend_type)
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

TEMP_data 	= netcdf.Dataset(directory+'Ocean/SOM_index.nc', 'r')

time 		=	TEMP_data.variables['time'][:]     	
SOM		=	TEMP_data.variables['SOM'][:] 

TEMP_data.close()

TEMP_data 	= netcdf.Dataset(directory+'Ocean/SOM_index_star.nc', 'r')

time 		=	TEMP_data.variables['time'][:]     	
SOM_2		=	TEMP_data.variables['SOM'][:] 

TEMP_data.close()


TEMP_data = netcdf.Dataset(directory+'Ocean/PC_SOM.nc', 'r')

#Writing data to correct variable	    	  	
time		= TEMP_data.variables['time'][:]    
PC		= TEMP_data.variables['PC_SOM'][:] 						

TEMP_data.close()

#Normalise by maximum
PC		= PC / (fabs(PC).max())

#-----------------------------------------------------------------------------------------

if trend_type > 0:
	
	SOM 		= TrendRemover(time, SOM, trend_type)
	SOM_2		= TrendRemover(time, SOM_2, trend_type)
	PC		= TrendRemover(time, PC, trend_type)

if moving_average > 1:

	SOM		= MovingAverage(time, SOM, moving_average)[1]
	SOM_2		= MovingAverage(time, SOM_2, moving_average)[1]
	time, PC	= MovingAverage(time, PC, moving_average)


#-----------------------------------------------------------------------------------------
#----------------------------------------PLOTTING-----------------------------------------
#-----------------------------------------------------------------------------------------	

left, width = .01, .5
bottom, height = 0.01, .5
right = left + width
top = bottom + height

fig, ax1 = plt.subplots()

#ax1.fill_between(np.arange(datetime.datetime(158, 1, 1).toordinal(), datetime.datetime(160, 1, 1).toordinal()), -1, 1, alpha=0.25, edgecolor='orange', facecolor='orange')
#ax1.fill_between(np.arange(datetime.datetime(178, 1, 1).toordinal(), datetime.datetime(183, 1, 1).toordinal()), -1, 1, alpha=0.25, edgecolor='orange', facecolor='orange')
#ax1.fill_between(np.arange(datetime.datetime(205, 1, 1).toordinal(), datetime.datetime(209, 1, 1).toordinal()), -1, 1, alpha=0.25, edgecolor='orange', facecolor='orange')
#ax1.fill_between(np.arange(datetime.datetime(231, 1, 1).toordinal(), datetime.datetime(237, 1, 1).toordinal()), -1, 1, alpha=0.25, edgecolor='orange', facecolor='orange')


ax1.grid()

#SOM_graph 	= ax1.plot_date(time, SOM, '-k', linewidth = 2.0, label = 'SOM')
SOM_graph_2 	= ax1.plot_date(time, SOM_2, '-b', linewidth = 2.0, label = 'SOM$^{*}$')

ax1.set_xlabel('Model year')
ax1.set_ylim(-0.4, 0.4)
ax1.set_ylabel('Temperature anomaly ($^{\circ}$C)')

ax2 		= ax1.twinx()

#SOM graph is plotted in second axis (for legend)
SOM_graph 	= ax2.plot_date(time, SOM * 1.0 / 0.4 , '-k', linewidth = 2.0, label = 'SOM')
PC_graph 	= ax2.plot_date(time, PC, '-r', linewidth = 2.0, label = 'PC 1')


ax2.set_ylim(-1, 1)
ax2.set_ylabel('Amplitude', color ='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')

graphs	      = SOM_graph + SOM_graph_2 + PC_graph
legend_labels = [l.get_label() for l in graphs]

ax1.legend(graphs, legend_labels, loc='lower left', ncol=1, fancybox=True, shadow=False, numpoints = 1)

ax1.set_xticks([datetime.datetime(160, 1, 1).toordinal(), 
		datetime.datetime(170, 1, 1).toordinal(),
		datetime.datetime(180, 1, 1).toordinal(),
		datetime.datetime(190, 1, 1).toordinal(), 
		datetime.datetime(200, 1, 1).toordinal(),
		datetime.datetime(210, 1, 1).toordinal(), 
		datetime.datetime(220, 1, 1).toordinal(),
		datetime.datetime(230, 1, 1).toordinal(),
		datetime.datetime(240, 1, 1).toordinal()])

ax1.set_title('d) CESM, SOM indices')
show()
