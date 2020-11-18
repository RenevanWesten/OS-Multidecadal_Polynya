#Program performs a red-noise test on the SOM

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

def MonthConverter(X):

	V = (1.0/(X)/12.0)

	return ["%.0f" % z for z in V]

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

trend_type 	= 0
moving_average 	= 0

surrogate	= 10000

label_level	= 10**6.0
period_min	= 5
#-----------------------------------------------------------------------------------------

TEMP_data = netcdf.Dataset(directory+'Ocean/PC_SOM.nc', 'r')

#Writing data to correct variable	    	  	
time		= TEMP_data.variables['time'][:]    
SOM		= TEMP_data.variables['PC_SOM'][:] 						

TEMP_data.close()

#-----------------------------------------------------------------------------------------

if trend_type > 0:
	print 'Data is detrended'
	SOM 		= TrendRemover(time, SOM, trend_type)

if moving_average > 1:
	print 'Data is smoothed'
	time, SOM	= MovingAverage(time, SOM, moving_average)

#SOM	= MonthRemover(time, SOM)
SOM	= (SOM - np.mean(SOM)) / np.std(SOM)

time, data_series	= time, SOM


#-----------------------------------------------------------------------------------------

#First determine the auto-correlation for each time series
auto_lag	= np.arange(250)
auto_corr	= np.zeros(250)

for lag_i in range(len(auto_lag)):

	#Determine the auto-correlation
	auto_corr[lag_i] = np.corrcoef(data_series[0:len(data_series)-lag_i], data_series[lag_i:])[0][1]

#Determine the e-folding time for each time series and keep the maximum
e_1	=  np.where(auto_corr < 1.0/np.e)[0][0]

#Determine the first coefficient for the AR(1) process
a 	= -1.0/(e_1) + 1.0

#Determine the variance in the time series and the last coefficeint for the AR(1) process
var	= np.var(data_series)
b	= np.sqrt((1.0 - a**2.0) * var)

#Make the Fourier Spectrum
freq_series 		= fft(data_series) 					#Take fourier spectrum
freq_series 		= ((real(freq_series)**2.0) + (imag(freq_series)**2.0)) #Determine power law (absolute value)
freq			= fftfreq(len(data_series))
freq_series_original 	= freq_series[:freq.argmax()]
freq			= freq[:freq.argmax()]

#Make surrogate fourier spectrum
surrogate_fourier	= ma.masked_all((surrogate, len(freq)))
surrogate_freq_90	= ma.masked_all(len(freq))
surrogate_freq_95	= ma.masked_all(len(freq))
surrogate_freq_99	= ma.masked_all(len(freq))

for surrogate_i in range(surrogate):
	#Now generate surrogate data
	spin_up		= 300
	dummy_data	= np.zeros(len(time))
	signal		= 0.0
	white_noise	= np.random.normal(0, 1, spin_up + len(dummy_data))

	for time_i in range(spin_up + len(dummy_data)):
		#Generate Red-noise spectrum with a spin-up (to remove the initial values)
		signal	= a * signal + b * white_noise[time_i]

		if time_i >= spin_up:
			#After the spin up, save the data
			dummy_data[time_i - spin_up]	= signal

	#Take the Fourier Spectrum
	freq_series 			= fft(dummy_data) 					#Take fourier spectrum
	freq_series 			= ((real(freq_series)**2.0) + (imag(freq_series)**2.0)) #Determine power law (absolute value)
	freq_series 			= freq_series[:len(freq)]
	surrogate_fourier[surrogate_i]	= freq_series


for freq_i in range(len(freq)):
	#Per frequency, sort the power spectrum
	surrogate_sort	= sorted(surrogate_fourier[:, freq_i])
	
	#Take the 90, 95, 99% confidence levels based on all the surrogate data
	surrogate_freq_90[freq_i]	= surrogate_sort[int(90.0 / 100.0 * surrogate)]
	surrogate_freq_95[freq_i]	= surrogate_sort[int(95.0 / 100.0 * surrogate)]
	surrogate_freq_99[freq_i]	= surrogate_sort[int(99.0 / 100.0 * surrogate)]

fig, ax1 = plt.subplots(figsize = (8.5, 6.375))

ax1.set_xlim([1.0/720,5.0/120.0])
ax1.set_xlabel('Frequency (month$^{-1}$)')
ax2 = ax1.twiny()

new_tick_locations = np.array([1.0/600, 1.0/360.0, 1.0/240, 1.0/120, 1.0/60, 1.0/36.0, 5.0/120])
ax2.set_xlabel('Period (years)')
ax2.set_xscale('log')

ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(MonthConverter(new_tick_locations))
ax2.grid()

ax1.plot(freq, freq_series_original, '-k', linewidth = 2)
sig_90_graph	= ax1.plot(freq, surrogate_freq_90, '-c', linewidth = 2, label = '90$\%$-Cl')
sig_95_graph	= ax1.plot(freq, surrogate_freq_95, '-b', linewidth = 2, label = '95$\%$-Cl')
sig_99_graph	= ax1.plot(freq, surrogate_freq_99, '-r', linewidth = 2, label = '99$\%$-Cl')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(10, 10000000)
ax1.set_ylabel('Power')


graphs	      	= sig_99_graph + sig_95_graph + sig_90_graph
legend_labels 	= [l.get_label() for l in graphs]
legend      	= legend(graphs, legend_labels, loc='upper right',
		  ncol=1, fancybox=True, numpoints = 1)


#Plot the period of significance
for freq_i in range(3, len(freq)):
	#Indicate period of significance
	period	= int(round(1.0 / freq[freq_i] / 12.0, 0))

	if period < period_min:
		break

	if freq_series_original[freq_i] > surrogate_freq_99[freq_i]:
		#99% confidence level
		color_sig	= 'r'

	elif freq_series_original[freq_i] > surrogate_freq_95[freq_i]:
		#95% confidence level
		color_sig	= 'b'

	elif freq_series_original[freq_i] > surrogate_freq_90[freq_i]:
		#90% confidence level
		color_sig	= 'c'

	if freq_series_original[freq_i] > surrogate_freq_90[freq_i]:
		#Only plot above 90% confidence level
		ax1.text(freq[freq_i], label_level, str(period), horizontalalignment='center', verticalalignment='bottom', color = color_sig, fontsize=13)
		ax1.plot([freq[freq_i], freq[freq_i]], [freq_series_original[freq_i], label_level], ':'+color_sig, linewidth = 1.5)


left, width = .01, .5
bottom, height = 0.01, .5
right = left + width
top = bottom + height

ax1.text(0.01, 0.95, 'b) PC 1 SOM',
		horizontalalignment='left',
		verticalalignment='bottom',
		transform=ax1.transAxes)

show()
