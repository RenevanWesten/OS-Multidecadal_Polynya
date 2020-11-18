#Lag-correlation analysis between the SOM (star) index and OHC of the polynya region

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
directory 	= '../../../Data/CESM/'

def CorrelationLag(series_1, series_2, lag_length = 50, name_series_1 = 'Series 1', Plot = True):
	"""Returns the correlation coefficient while lagging both series to each other
	time can be selected and is by default 50 time units"""
	
	lag_1, lag_2 = [], []
	correlation_1, correlation_2 = [], []
	error_low_1, error_high_1, error_low_2, error_high_2 = [], [], [], []

	for lag_i in range(lag_length+1): #Lagging the time series in both directions
		correlation_1.append(np.corrcoef(series_1[:len(series_1)-lag_i],series_2[lag_i:])[0][1])
		lag_1.append(lag_i)
		correlation_2.append(np.corrcoef(series_2[:len(series_2)-lag_i],series_1[lag_i:])[0][1])
		lag_2.append(-lag_i)
		
		#Determine errors on correlation
		low, high = FisherTransformation(correlation_1[-1], len(series_1[:len(series_1)-lag_i]))
		error_low_1.append(correlation_1[-1] - low), error_high_1.append(high - correlation_1[-1])
		low, high = FisherTransformation(correlation_2[-1], len(series_2[:len(series_2)-lag_i]))
		error_low_2.append(correlation_2[-1] - low), error_high_2.append(high - correlation_2[-1])
	
	#Pasting the time series together where lag = 0 is removed from one series (double)
	lag, correlation, error_low, error_high =  lag_2[::-1][:-1] + lag_1, correlation_2[::-1][:-1] + correlation_1, error_low_2[::-1][:-1] + error_low_1, error_high_2[::-1][:-1] + error_high_1
	
	if Plot == False:
		return np.asarray(lag), np.asarray(correlation), np.asarray(error_low), np.asarray(error_high)

	print 'Correlation between the time series'
	print '		Lag 		=', lag[len(lag) - lag_length - 1]
	print '		Correlation 	=', correlation[len(lag) - lag_length - 1]
	print '		Lower error	=', error_low[len(lag) - lag_length - 1]
	print '		Upper error	=', error_high[len(lag) - lag_length - 1]

	print
	print 'Maximum correlation between the time series'
	print '		Lag 		=', lag[argmax(correlation)]
	print '		Correlation 	=', correlation[argmax(correlation)]
	print '		Lower error	=', error_low[argmax(correlation)]
	print '		Upper error	=', error_high[argmax(correlation)]

	print
	print 'Minimum correlation between the time series'
	print '		Lag 		=', lag[argmin(correlation)]
	print '		Correlation 	=', correlation[argmin(correlation)]
	print '		Lower error	=', error_low[argmin(correlation)]
	print '		Upper error	=', error_high[argmin(correlation)]

	bottom, height 	= -0.95, 2.0
	middle		= 0
	left		= -1.0/3.0 * max(lag)
	right		= 1.0/25.0 * max(lag)

	axvline(x=0,color='k',ls='dashed')
	axhline(y=0,color='k',ls='dashed')

	text(right, bottom, name_series_1+' leads', horizontalalignment='left',verticalalignment='bottom')
	text(left, bottom, name_series_1+' lags',horizontalalignment='left',verticalalignment='bottom')		

	plot(lag, correlation, '-r', linewidth = 2.0)
	fill_between(lag, np.asarray(correlation) - np.asarray(error_low), np.asarray(correlation) + np.asarray(error_high), alpha=0.4, edgecolor='r', facecolor='r')	

	xlabel('Lag (months)')
	ylabel("Correlation coefficient")
	ylim([-1, 1])

	show()

def FisherTransformation(correlation, size_length):
	"""Determines error on correlation coefficient"""
	
	z = 0.5 * log((1.0 + correlation) / (1.0 - correlation)) #Fisher transformation
	
	error = 1.0/(sqrt(size_length - 3.0)) #Standard error

	err_low, err_up = z - 1.96 * error, z + 1.96 * error #95 Confidence interval
	
	err_low = (e**(2.0 * err_low) - 1.0) / (e**(2.0 * err_low) + 1.0) #Fisher inverse transformation
	err_up = (e**(2.0 * err_up) - 1.0) / (e**(2.0 * err_up) + 1.0)
	
	return err_low, err_up

def SignificanceCorrelationLag(series_1, series_2, lag_length = 50, name_series_1 = 'Series 1'):
	#Computes the significance of the lag

	lag, corr, corr_low, corr_high = CorrelationLag(series_1, series_2, lag_length, name_series_1, Plot = False)
	
	#Make empty array's to keep track of your significant lag
	significance_level 	= ma.masked_all(len(lag))
	
	for lag_i in range(len(lag)):
		print lag_i, len(lag)

		#First determine the auto-correlation for each time series
		auto_lag	= np.arange(250)
		auto_corr_1	= np.zeros(250)
		auto_corr_2	= np.zeros(250)

		for lag_j in range(len(auto_lag)):

			if lag[lag_i] < 0:
				#If series 2 leads series 1
				auto_series_1	= series_1[-lag[lag_i]:] #Lag is negative, adjust with minus
				auto_series_2 	= series_2[:len(series_1) + lag[lag_i]]

			else:
				#If series 1 leads series 2
				auto_series_1	= series_1[:len(series_1) - lag[lag_i]]
				auto_series_2 	= series_2[lag[lag_i]:]

			#Determine the auto-correlation
			auto_corr_1[lag_j] = np.corrcoef(auto_series_1[0:len(auto_series_1)-lag_j], auto_series_1[lag_j:])[0][1]
			auto_corr_2[lag_j] = np.corrcoef(auto_series_2[0:len(auto_series_1)-lag_j], auto_series_2[lag_j:])[0][1]

		#Determine the e-folding time for each time series and keep the maximum
		e_1	= np.where(auto_corr_1 < 1.0/np.e)[0][0]
		e_2	= np.where(auto_corr_2 < 1.0/np.e)[0][0]
		e_max	= max(e_1, e_2)

		#Determine the degrees of freedom and corresponding critical value
		dof	= len(auto_series_1) / e_max
		t_crit 	= stats.t.ppf(1 - 0.05, dof - 2)

		#Determine the minimum correlation to exceed the critical value
		corr_min 			= t_crit * np.sqrt( 1.0 / (dof - 2.0 + t_crit**2.0))
		significance_level[lag_i]	= corr_min

	return lag, corr, corr_low, corr_high, significance_level


def TrendRemover(time, data, trend_type):
	"""Removes trend of choice"""
	
	rank = polyfit(time, data, trend_type)
	fitting = 0.0 
		
	for rank_i in range(len(rank)):
			
		fitting += rank[rank_i] * (time**(len(rank) - 1 - rank_i))

	data -= fitting
	
	return data

def MonthRemover(time, data):
	"""Removes monthly climatology mean"""

	for month_i in range(12):
		#Loop over each month

		month_index	= np.arange(month_i, len(time), 12)
	
		#Determine the climatology mean per month
		data_month	= np.mean(data[month_index])

		#Subtract the climatology monthly mean
		data[month_index]	= data[month_index] - data_month
		
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
season_remove	= 0
moving_average 	= 60 	#Running mean

#-----------------------------------------------------------------------------------------

TEMP_data 	= netcdf.Dataset(directory+'Ocean/SOM_index_star.nc', 'r')

time 		= TEMP_data.variables['time'][:]     	
SOM		= TEMP_data.variables['SOM'][:] 

TEMP_data.close()

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

if trend_type > 0:
	
	SOM 	= TrendRemover(time, SOM, trend_type)
	OHC_1 	= TrendRemover(time, OHC_1, trend_type)
	OHC_2	= TrendRemover(time, OHC_2, trend_type)

if season_remove == 1:

	SOM	= MonthRemover(time, SOM)
	OHC_1	= MonthRemover(time, OHC_1)
	OHC_2	= MonthRemover(time, OHC_2)

if moving_average > 1:

	print 'Moving average is determined'
	time_2, SOM	= MovingAverage(time, SOM, moving_average)
	time_2, OHC_1	= MovingAverage(time, OHC_1, moving_average)
	time, OHC_2	= MovingAverage(time, OHC_2, moving_average)

#-----------------------------------------------------------------------------------------
#----------------------------------------PLOTTING-----------------------------------------
#-----------------------------------------------------------------------------------------	

left, width = .01, .5
bottom, height = 0.01, .5
right = left + width
top = bottom + height

fig, ax1 = plt.subplots()
for tl in ax1.get_yticklabels():
    tl.set_color('k')

ax1.grid()

OHC_1_graph = ax1.plot_date(time, OHC_1 * 5.0, '-b', linewidth = 2.0, label = 'OHC 100 m (5x)')
ax1.set_ylabel('Ocean heat content anomaly (ZJ)', color = 'k')
ax1.set_xlabel('Model year')


OHC_2_graph = ax1.plot_date(time, OHC_2, '-k', linewidth = 2.0, label = 'OHC 1000 m')
ax1.set_xlabel('Model year')
ax1.set_ylim([-0.35, 0.35])

ax2 = ax1.twinx()
SOM_graph = ax2.plot_date(time, SOM, '-r', linewidth = 2.0, label = 'SOM')
ax2.set_ylabel('Temperature anomaly ($^{\circ}$C)', color ='r')
ax2.set_ylim([-0.35, 0.35])

for tl in ax2.get_yticklabels():
    tl.set_color('r')

graphs	      = OHC_1_graph + OHC_2_graph + SOM_graph
legend_labels = [l.get_label() for l in graphs]
ax1.legend(graphs, legend_labels, loc='upper center', bbox_to_anchor=(0.5, 1.02),
 		  ncol=4, fancybox=True, shadow=False, numpoints = 1)

show()
#-----------------------------------------------------------------------------------------


lag_1, corr_1, corr_low_1, corr_high_1, significance_level_1 = SignificanceCorrelationLag(SOM, OHC_1, 150, 'SOM')
lag_2, corr_2, corr_low_2, corr_high_2, significance_level_2 = SignificanceCorrelationLag(SOM, OHC_2, 150, 'SOM')


#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

bottom, height 	= -0.95, 2.0
middle		= 0
left		= -1.0
right		= 1.0

ax.grid()
ax.axvline(x=0,color='k')
ax.axhline(y=0,color='k')	

ax.text(15, -0.9, 'SOM$^{*}$ leads', horizontalalignment='left',verticalalignment='bottom')
ax.text(-15, -0.87, r'$\bar{H}_{100}$ leads',horizontalalignment='right',verticalalignment='bottom')	
ax.text(-15, -0.95, r'$\bar{H}_{1000}$ leads',horizontalalignment='right',verticalalignment='bottom')	
	

SOM_vs_100_graph = ax.plot(lag_1, corr_1, '-r', linewidth = 2.0, label = r'SOM$^{*}$  vs. $\bar{H}_{100}$')
ax.fill_between(lag_1, corr_1 - corr_low_1, corr_1 + corr_high_1, alpha=0.4, edgecolor='r', facecolor='r')	

SOM_vs_1000_graph = ax.plot(lag_2, corr_2, '-b', linewidth = 2.0, label = r'SOM$^{*}$  vs. $\bar{H}_{1000}$')
ax.fill_between(lag_2, corr_2 - corr_low_2, corr_2 + corr_high_2, alpha=0.4, edgecolor='b', facecolor='b')	

ax.plot(lag_1, significance_level_1, '--r')
ax.plot(lag_2, significance_level_2, '--b')
ax.plot(lag_1, -significance_level_1, '--r')
ax.plot(lag_2, -significance_level_2, '--b')

ax.set_xlabel('Lag (months)')
ax.set_ylabel("Lag-correlation coefficient")
ax.set_xlim(-150, 150)
ax.set_ylim(-1, 1)

graphs	      = SOM_vs_100_graph + SOM_vs_1000_graph

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)

ax.set_title('a) Lag-correlation analysis')

show()