#Program determines the average position of the Maud Rise Polynyas

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

ICE_data = netcdf.Dataset(directory+'Ice/Polynya_fields.nc', 'r')

#Writing data to correct variable
time		= ICE_data.variables['time'][:]		
lon		= ICE_data.variables['lon'][:]     		
lat		= ICE_data.variables['lat'][:]     	
position	= ICE_data.variables['polynya'][:] 

ICE_data.close()

#-----------------------------------------------------------------------------------------
levels 		= [0.9, 1.0, 1.1]
x, y 		= np.meshgrid(lon, lat)
CS		= contourf(x, y, position[27], levels)
contours	= CS.allsegs

x_polynya	= contours[0][0][:,0]  	#Take the x-coordinate of contour
y_polynya	= contours[0][0][:,1]	#Take the y-coordinate of contour
close('all')

position_average	= np.sum(position, axis = 0) / len(time) 

fig, ax	= subplots()
levels 	= np.arange(0, 0.71, 0.05)
x, y 	= np.meshgrid(lon, -lat)
CS 	= contourf(x, y, position_average, levels, extend = 'max', cmap = 'Spectral_r')
cbar	= colorbar(CS, ticks = np.arange(0, 0.71, 0.1))
cbar.set_label('Polynya probability function')
ax.set_ylim(68, 62)
ax.set_xlim(-1.5, 20)
ax.set_xlabel('Longitude ($^{\circ}$E)')
ax.set_ylabel('Latitude ($^{\circ}$S)')
grid()

plot([2, 2, 11, 11, 2], [66.5, 63.5, 63.5, 66.5, 66.5], '--k', linewidth = 3.0)
ax.set_title('b) Polynya probability function')

#-----------------------------------------------------------------------------------------

polynya_10_october_2017 = np.loadtxt('../../../Data/SSMI-SSMR/Ice/Polynya_SSMI-SSMR_10_October_2017.txt')
polynya_10_october_2017_lon, polynya_10_october_2017_lat = polynya_10_october_2017[0], polynya_10_october_2017[1]

#-----------------------------------------------------------------------------------------

#Plot individual polynyas
Polynya_graph	= plot(x_polynya, -y_polynya, '-k', linewidth = 3.0, label = 'MRP')	
OBS_graph	= plot(polynya_10_october_2017_lon, -polynya_10_october_2017_lat, '-r', linewidth = 3.0, label = 'MRP 10 October 2017')

#-----------------------------------------------------------------------------------------

graphs	      = Polynya_graph + OBS_graph
legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower center',
 		  ncol=3, fancybox=True, shadow=False, numpoints = 1)

show()
