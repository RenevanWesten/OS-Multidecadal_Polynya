#Program plots the heat and salt advection

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
directory 		= '../../../Data/CESM/'

def YearlyConverter(time, data, month_start = 1, month_end = 12):
	"""Determines yearly averaged, over different months of choice,
	default is set to January - December"""

	#Take twice the amount of years for the month day
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31., 31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days[month_start - 1:month_end]
	month_days	= month_days / np.sum(month_days)
	
	if month_end <= 12:
		#Normal average over a single year, for example, February 100 - December 100
		time_year		= np.zeros(len(time) / 12)

	else:
		#If you take the average, for example, over November 100 - May 101
		#Take year 101 as the average over this period
		#There is one year less compared to the period analysed
		time_year		= np.zeros(len(time) / 12 - 1)

	#-----------------------------------------------------------------------------------------
	data_year	= ma.masked_all(len(time_year))

	for year_i in range(len(time_year)):
		#Determine the SSH over the selected months

		#Determine the time mean
		time_year[year_i] 		= np.sum(time[year_i * 12 + month_start - 1: year_i * 12 + month_end] * month_days, axis = 0)

		#Determine the time mean over the months of choice
		data_year[year_i]		= np.sum(data[year_i * 12 + month_start - 1: year_i * 12 + month_end] * month_days, axis = 0)

	return time_year, data_year

def SignificanceCorrelation(series_1, series_2):
	#Computes the significance of the lag

	corr	= np.corrcoef(series_1, series_2)[0][1]

	#First determine the auto-correlation for each time series
	auto_lag	= np.arange(50)
	auto_corr_1	= np.zeros(50)
	auto_corr_2	= np.zeros(50)

	for lag_i in range(len(auto_lag)):

		#Determine the auto-correlation
		auto_corr_1[lag_i] = np.corrcoef(series_1[0:len(series_1)-lag_i], series_1[lag_i:])[0][1]
		auto_corr_2[lag_i] = np.corrcoef(series_2[0:len(series_2)-lag_i], series_2[lag_i:])[0][1]
	
	#Determine the e-folding time for each time series and keep the maximum
	e_1	= np.where(auto_corr_1 < 1.0/np.e)[0][0]
	e_2	= np.where(auto_corr_2 < 1.0/np.e)[0][0]
	e_max	= max(e_1, e_2)

	#Determine the degrees of freedom and corresponding critical value
	dof	= len(series_1) / e_max
	t_crit 	= stats.t.ppf(1 - 0.05, dof - 2)

	#Determine the minimum correlation to exceed the critical value
	corr_min 			= t_crit * np.sqrt( 1.0 / (dof - 2.0 + t_crit**2.0))

	if fabs(corr) > corr_min:
		print 'Significant correlation'
		print 'Correlation coefficient:', corr

	print

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

month_start	= 1
month_end	= 12

#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/Polynya_HEAT_advection.nc', 'r')

#Writing data to correct variable	
time		= HEAT_data.variables['time'][:]     	  	
depth		= HEAT_data.variables['depth'][:]     	  	
layer		= HEAT_data.variables['layer'][:]     	  	
v_heat_2	= HEAT_data.variables['VHEAT_NORTH'][:] 		
u_heat_2	= HEAT_data.variables['UHEAT_EAST'][:] 		
v_heat_1	= HEAT_data.variables['VHEAT_SOUTH'][:] 		
u_heat_1	= HEAT_data.variables['UHEAT_WEST'][:] 		

HEAT_data.close()

#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/Polynya_SALT_advection.nc', 'r')

#Writing data to correct variable	
time		= HEAT_data.variables['time'][:]     	  	
depth		= HEAT_data.variables['depth'][:]     	  	
layer		= HEAT_data.variables['layer'][:]     	  	
v_salt_2	= HEAT_data.variables['VSALT_NORTH'][:] 		
u_salt_2	= HEAT_data.variables['USALT_EAST'][:] 		
v_salt_1	= HEAT_data.variables['VSALT_SOUTH'][:] 		
u_salt_1	= HEAT_data.variables['USALT_WEST'][:] 		

HEAT_data.close()

#-----------------------------------------------------------------------------------------

ICE_data 	= netcdf.Dataset(directory+'Ice/Polynya_ice_fraction_thickness.nc', 'r')

time_ice 	= ICE_data.variables['time'][:]     	
ice		= ICE_data.variables['hi'][:] 

ICE_data.close()

#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/Polynya_TEMP_depth.nc', 'r')

#Writing data to correct variable	
time		= HEAT_data.variables['time'][:]     	  	
depth		= HEAT_data.variables['depth'][:]     	  	
temp_depth	= HEAT_data.variables['TEMP'][:] 			

HEAT_data.close()

#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/Polynya_SALT_depth.nc', 'r')

#Writing data to correct variable	
time		= HEAT_data.variables['time'][:]     	  	
depth		= HEAT_data.variables['depth'][:]     	  	
salt_depth	= HEAT_data.variables['SALT'][:] 			

HEAT_data.close()

#-----------------------------------------------------------------------------------------

ice_year		= np.zeros(len(time) / 12)

for time_i in range(len(time) / 12):
	#Take the maximum for each year
	ice_year[time_i]		= np.max(ice[time_i * 12: time_i * 12 + 12])

#-----------------------------------------------------------------------------------------

depth_min_index_0	= (fabs(depth - 0)).argmin()
depth_max_index_100	= (fabs(depth - 100)).argmin() + 1
depth_min_index_200	= (fabs(depth - 200)).argmin()
depth_max_index_1000	= (fabs(depth - 1000)).argmin() + 1

#Take the normalised area over the depth layers
temp_depth_100	= temp_depth[:, depth_min_index_0:depth_max_index_100]
temp_depth_1000	= temp_depth[:, depth_min_index_200:depth_max_index_1000]
salt_surface	= salt_depth[:, 0]
salt_depth_100	= salt_depth[:, depth_min_index_0:depth_max_index_100]
salt_depth_1000	= salt_depth[:, depth_min_index_200:depth_max_index_1000]

temp_depth_100	= np.sum(temp_depth_100  *  layer[depth_min_index_0:depth_max_index_100] / np.sum(layer[depth_min_index_0:depth_max_index_100]), axis = 1) 
temp_depth_1000	= np.sum(temp_depth_1000 *  layer[depth_min_index_200:depth_max_index_1000]  / np.sum(layer[depth_min_index_200:depth_max_index_1000]), axis = 1)
salt_depth_100	= np.sum(salt_depth_100  *  layer[depth_min_index_0:depth_max_index_100] / np.sum(layer[depth_min_index_0:depth_max_index_100]), axis = 1) 
salt_depth_1000	= np.sum(salt_depth_1000 *  layer[depth_min_index_200:depth_max_index_1000]  / np.sum(layer[depth_min_index_200:depth_max_index_1000]), axis = 1)

u_heat_all	= u_heat_1 - u_heat_2 + v_heat_1 - v_heat_2
u_heat_surface	= u_heat_all[:, 0]
u_heat_all_100	= np.sum(u_heat_all[:, depth_min_index_0:depth_max_index_100], axis = 1)
u_heat_all_1000	= np.sum(u_heat_all[:, depth_min_index_200:depth_max_index_1000], axis = 1)

u_salt_all	= u_salt_1 - u_salt_2 + v_salt_1 - v_salt_2
u_salt_surface	= u_salt_all[:, 0]
u_salt_all_100	= np.sum(u_salt_all[:, depth_min_index_0:depth_max_index_100], axis = 1)
u_salt_all_1000	= np.sum(u_salt_all[:, depth_min_index_200:depth_max_index_1000], axis = 1)

u_heat_1	= np.sum(u_heat_1[:, depth_min_index_200:depth_max_index_1000], axis = 1)
u_heat_2	= np.sum(u_heat_2[:, depth_min_index_200:depth_max_index_1000], axis = 1)
v_heat_1	= np.sum(v_heat_1[:, depth_min_index_200:depth_max_index_1000], axis = 1)
v_heat_2	= np.sum(v_heat_2[:, depth_min_index_200:depth_max_index_1000], axis = 1)

#Take the yearly averages
time_year, u_heat_all_100		= YearlyConverter(time, u_heat_all_100, month_start, month_end)
time_year, u_heat_all_1000		= YearlyConverter(time, u_heat_all_1000, month_start, month_end)
time_year, u_salt_surface		= YearlyConverter(time, u_salt_surface, month_start, month_end)
time_year, u_salt_all_100		= YearlyConverter(time, u_salt_all_100, month_start, month_end)
time_year, u_salt_all_1000		= YearlyConverter(time, u_salt_all_1000, month_start, month_end)

time_year, u_heat_1			= YearlyConverter(time, u_heat_1, month_start, month_end)
time_year, u_heat_2			= YearlyConverter(time, u_heat_2, month_start, month_end)
time_year, v_heat_1			= YearlyConverter(time, v_heat_1, month_start, month_end)
time_year, v_heat_2			= YearlyConverter(time, v_heat_2, month_start, month_end)

time_year, temp_depth_diff		= YearlyConverter(time, temp_depth_1000 - temp_depth_100, month_start, month_end)
time_year, salt_depth_diff		= YearlyConverter(time, salt_depth_1000 - salt_depth_100, month_start, month_end)
time_year, salt_surface			= YearlyConverter(time, salt_surface, month_start, month_end)
time_year, temp_depth_100		= YearlyConverter(time, temp_depth_100, month_start, month_end)
time_year, temp_depth_1000		= YearlyConverter(time, temp_depth_1000, month_start, month_end)
time_year, salt_depth_100		= YearlyConverter(time, salt_depth_100, month_start, month_end)
time_year, salt_depth_1000		= YearlyConverter(time, salt_depth_1000, month_start, month_end)

#-----------------------------------------------------------------------------------------
SignificanceCorrelation(u_heat_all_1000, temp_depth_diff)
SignificanceCorrelation(u_salt_all_1000, salt_depth_diff)
#-----------------------------------------------------------------------------------------

fig, ax1 	= plt.subplots()

ax1.fill_between(np.arange(datetime.datetime(231, 1, 1).toordinal(), datetime.datetime(237, 12, 31).toordinal()), -100, 100, alpha=0.25, edgecolor='orange', facecolor='orange')

heat_100_graph		= ax1.plot_date(time_year, u_heat_all_100 * 5.0 / (10**12.0), '-r', linewidth = 2.0, label = r'$\bar{F}_{100}^{\mathrm{Net}}$ (5x)')
heat_1000_graph		= ax1.plot_date(time_year, u_heat_all_1000 * 5.0 / (10**12.0), '-k', linewidth = 2.0, label = r'$\bar{F}_{1000}^{\mathrm{Net}}$ (5x)')

heat_north_graph	= ax1.plot_date(time_year, - v_heat_2 / (10**12.0), '--r', linewidth = 2.0, label = r'$\bar{F}_{1000}^{\mathrm{North}}$')
heat_east_graph		= ax1.plot_date(time_year, - u_heat_2 / (10**12.0), '-b', linewidth = 2.0, label = r'$\bar{F}_{1000}^{\mathrm{East}}$')
heat_south_graph	= ax1.plot_date(time_year, v_heat_1  / (10**12.0), '--b', linewidth = 2.0, label = r'$\bar{F}_{1000}^{\mathrm{South}}$')
heat_west_graph		= ax1.plot_date(time_year, u_heat_1 / (10**12.0), '--k', linewidth = 2.0, label = r'$\bar{F}_{1000}^{\mathrm{West}}$')

ax1.set_xlabel('Model year')
ax1.set_ylabel('Advective heat (TW)', color = 'k')
ax1.set_ylim(-60, 60)

ax1.grid()

ax2 = ax1.twinx()

temp_diff_graph		= ax2.plot_date(time_year, temp_depth_diff, '-c', linewidth = 2.0, label = r'$\bar{T}_{1000} - \bar{T}_{100}$')

ax2.set_ylim(0, 3.0)
ax2.set_ylabel('Temperature difference ($^{\circ}$C)', color ='k')

ax1.set_xticks([datetime.datetime(210, 1, 1).toordinal(), 
		datetime.datetime(215, 1, 1).toordinal(),
		datetime.datetime(220, 1, 1).toordinal(),
		datetime.datetime(225, 1, 1).toordinal(), 
		datetime.datetime(230, 1, 1).toordinal(),
		datetime.datetime(235, 1, 1).toordinal(), 
		datetime.datetime(240, 1, 1).toordinal()])

ax1.set_xlim(datetime.datetime(210, 1, 1).toordinal(),datetime.datetime(240, 1, 1).toordinal())

graphs	      = heat_100_graph + heat_1000_graph + temp_diff_graph + heat_east_graph + heat_west_graph + heat_north_graph + heat_south_graph

legend_labels = [l.get_label() for l in graphs]
legend_1      = ax2.legend(graphs, legend_labels, loc='lower left', 
 		  ncol=3, fancybox=True, numpoints = 1)

ax1.set_title('c) Advective heat and temperature difference')

#-----------------------------------------------------------------------------------------

fig, ax1 	= plt.subplots()

ax1.fill_between(np.arange(datetime.datetime(231, 1, 1).toordinal(), datetime.datetime(237, 12, 31).toordinal()), -100, 100, alpha=0.25, edgecolor='orange', facecolor='orange')

salt_100_graph		= ax1.plot_date(time_year, u_salt_all_100, '-r', linewidth = 2.0, label = r'$\bar{P}_{100}^{\mathrm{Net}}$')
salt_1000_graph		= ax1.plot_date(time_year, u_salt_all_1000, '-k', linewidth = 2.0, label = r'$\bar{P}_{1000}^{\mathrm{Net}}$')
salt_surface_graph	= ax1.plot_date(time_year, (salt_surface - salt_surface[60]) * 10.0, '-b', linewidth = 2.0, label = r'$S_{surface}$ (10x)')

ax1.set_xlabel('Model year')
ax1.set_ylabel('Advective salt (Sv Psu) and salinity anomaly (Psu)', color ='k')
ax1.grid()
ax1.set_ylim(-17.5, 17.5)

ax2 = ax1.twinx()

ice_graph		= ax2.plot_date(time_year, ice_year / 2, '--b', linewidth = 2.0, label = '$H_{ice}$ (0.5x)')
salt_diff_graph		= ax2.plot_date(time_year, salt_depth_diff, '-c', linewidth = 2.0, label = r'$\bar{S}_{1000} - \bar{S}_{100}$')

ax2.set_ylim(0, 0.65)
ax2.set_ylabel('Salinity difference (Psu) and sea-ice thickness (m)', color ='k')

ax1.set_xticks([datetime.datetime(210, 1, 1).toordinal(), 
		datetime.datetime(215, 1, 1).toordinal(),
		datetime.datetime(220, 1, 1).toordinal(),
		datetime.datetime(225, 1, 1).toordinal(), 
		datetime.datetime(230, 1, 1).toordinal(),
		datetime.datetime(235, 1, 1).toordinal(), 
		datetime.datetime(240, 1, 1).toordinal()])

ax1.set_xlim(datetime.datetime(210, 1, 1).toordinal(),datetime.datetime(240, 1, 1).toordinal())


graphs	      = salt_100_graph + salt_1000_graph + salt_surface_graph + ice_graph + salt_diff_graph

legend_labels = [l.get_label() for l in graphs]
legend_1      = ax1.legend(graphs, legend_labels, loc='upper center',
 		  ncol=3, fancybox=True, numpoints = 1)

ax1.set_title('d) Advective salt and salinity difference')

show()
