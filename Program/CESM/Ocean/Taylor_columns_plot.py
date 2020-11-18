#Program plots the Taylor Columns near the Polynya Region

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

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/Taylor_columns.nc', 'r')

#Writing data to correct variable	
time		= HEAT_data.variables['time'][:]     	  	
lon		= HEAT_data.variables['lon'][:] 
depth		= HEAT_data.variables['depth'][:]     	  	    	  	
mixed_max	= HEAT_data.variables['XMXL'][:] 			
dens		= HEAT_data.variables['PD'][:]			

HEAT_data.close()

#Take all September profiles
mixed_max	= mixed_max[np.arange(8, len(time), 12)]
dens		= dens[np.arange(8, len(time), 12)]

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between(lon, y1 = np.zeros(len(lon)) + depth[0], y2 = np.zeros(len(lon)) + depth[-1], color = 'gray', alpha = 0.50)

x, y	= np.meshgrid(lon, depth)
CS	= contourf(x, y, dens[75], levels = np.arange(0, 0.61, 0.025), extend = 'both', cmap = 'Spectral_r')
cbar	= colorbar(CS, ticks = np.arange(0, 0.61, 0.1))
cbar.set_label('Potential density with respect to surface (kg m$^{-3}$)')

ax.plot(lon, mixed_max[75], '-k', linewidth = 2.0)

ax.set_xlabel('Longitude')
ax.set_ylabel('Depth (m)')

ax.set_ylim(4000, 0)

ax.set_title('a) Potential density, September model year 225')

ax.set_xticks([-20, -10, 0, 10, 20, 30])
ax.set_xticklabels(['20$^{\circ}$W', '10$^{\circ}$W', '0$^{\circ}$', '10$^{\circ}$E', '20$^{\circ}$E', '30$^{\circ}$E'])

#-----------------------------------------------------------------------------------------

fig, ax		= subplots()
time		= np.arange(150, 251)

levels		= np.arange(0, 0.61, 0.025)
x, y		= np.meshgrid(lon, time)
CS		= contourf(x, y, dens[:, 21], levels, extend = 'both', cmap = 'Spectral_r')
cbar		= colorbar(CS, ticks = np.arange(0, 0.61, 0.1))
cbar.set_label('Potential density with respect to surface (kg m$^{-3}$)')

ax.set_xlabel('Longitude')
ax.set_ylabel('Model year')

ax.set_yticks(np.arange(150, 251, 10))

ax.set_xticks([-20, -10, 0, 10, 20, 30])
ax.set_xticklabels(['20$^{\circ}$W', '10$^{\circ}$W', '0$^{\circ}$', '10$^{\circ}$E', '20$^{\circ}$E', '30$^{\circ}$E'])

ax.set_title('b) September potential density at 900 m depth')

show()
#-----------------------------------------------------------------------------------------


