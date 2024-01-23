#Program plots the AMOC hysteresis

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker


#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

	fh.close()

	return time, transport

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

time, transport             = ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_hysteresis.nc')	
time_0600, transport_0600	= ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_0600.nc')	
time_1500, transport_1500	= ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_1500.nc')		
time_2900, transport_2900	= ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_2900.nc')	
time_3800, transport_3800	= ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_3800.nc')	
#-----------------------------------------------------------------------------------------

    
fig, ax	= subplots()

ax.fill_between([-100, 2500], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot([0.1/0.0003, 0.1/0.0003], [-5, 33], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.2/0.0003, 0.2/0.0003], [-5, 33], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.3/0.0003, 0.3/0.0003], [-5, 33], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.4/0.0003, 0.4/0.0003], [-5, 33], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.5/0.0003, 0.5/0.0003], [-5, 33], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.6/0.0003, 0.6/0.0003], [-5, 33], linestyle = '--', color = 'c', linewidth = 1)

ax.text(0.1/0.0003, 33, '0.1 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.2/0.0003, 33, '0.2 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.3/0.0003, 33, '0.3 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.4/0.0003, 33, '0.4 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.5/0.0003, 33, '0.5 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.6/0.0003, 33, '0.6 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)


ax.plot(4400 - time, transport, '-r', linewidth = 0.5)
ax.plot(time, transport, '-k', linewidth = 0.5)


y_error_0600		= np.zeros((2, 1)) + 1
y_error_0600[:, 0]	= np.mean(transport_0600[-50:]) - np.min(transport_0600[-50:]), np.max(transport_0600[-50:]) - np.mean(transport_0600[-50:])
ax.errorbar(600, np.mean(transport_0600[-50:]), color = 'b', marker = 's', markerfacecolor = 'dodgerblue', markeredgecolor = 'royalblue', yerr = y_error_0600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1500		= np.zeros((2, 1)) + 1
y_error_1500[:, 0]	= np.mean(transport_1500[-50:]) - np.min(transport_1500[-50:]), np.max(transport_1500[-50:]) - np.mean(transport_1500[-50:])
graph_3 		= ax.errorbar(1500, np.mean(transport_1500[-50:]), color = 'b', marker = 's', markerfacecolor = 'dodgerblue', markeredgecolor = 'royalblue', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')

y_error_2900		= np.zeros((2, 1)) + 1
y_error_2900[:, 0]	= np.mean(transport_2900[-50:]) - np.min(transport_2900[-50:]), np.max(transport_2900[-50:]) - np.mean(transport_2900[-50:])
ax.errorbar(4400 - 2900, np.mean(transport_2900[-50:]), color = 'b', marker = 's', markerfacecolor = 'dodgerblue', markeredgecolor = 'royalblue', yerr = y_error_2900, linewidth = 2.0, capsize=4)

y_error_3800		= np.zeros((2, 1)) + 1
y_error_3800[:, 0]	= np.mean(transport_3800[-50:]) - np.min(transport_3800[-50:]), np.max(transport_3800[-50:]) - np.mean(transport_3800[-50:])
ax.errorbar(4400 - 3800, np.mean(transport_3800[-50:]), color = 'b', marker = 's', markerfacecolor = 'dodgerblue', markeredgecolor = 'royalblue', yerr = y_error_3800, linewidth = 0.0, elinewidth = 2.0, capsize=4)

ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-2, 35)
ax.grid()

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = '1 - 2200')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = '2201 - 4400')

legend_1	= ax.legend(loc=(0.315, 0.19), ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

#-----------------------------------------------------------------------------------------

box1 = TextArea("1/", textprops=dict(color="k"))
box2 = TextArea("4400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.023, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("500/", textprops=dict(color="k"))
box2 = TextArea("3900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.165, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("1000/", textprops=dict(color="k"))
box2 = TextArea("3400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.379, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("1500/", textprops=dict(color="k"))
box2 = TextArea("2900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.605, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("2000/", textprops=dict(color="k"))
box2 = TextArea("2400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.835, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

ax.set_title('a) AMOC strength at 26$^{\circ}$N')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.63, 0.40, 0.35, 0.35], projection = ccrs.Orthographic(-30, 10))

ax2.coastlines(resolution='50m')
ax2.gridlines()
ax2.add_feature(cfeature.LAND, zorder=10)
ax2.set_global()


lon     = np.arange(-1, 360)
lat     = np.arange(-90, 91)
field   = np.ones((len(lat), len(lon))) * -0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())


lon     = np.arange(-100, -5)
lat     = np.arange(20, 43)
field   = np.ones((len(lat), len(lon))) * 0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

lon     = np.arange(-100, 3)
lat     = np.arange(42, 51)
field   = np.ones((len(lat), len(lon))) * 0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

ax2.text(320, 38, '$+F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=11, transform=ccrs.PlateCarree())
ax2.text(340, -10, '$-F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=11, transform=ccrs.PlateCarree())

x_1	= np.arange(-81, -9.99, 0.1)
y_1	= np.zeros(len(x_1)) + 26.0
y_2	= np.arange(24, 28.01, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(24, 28.01, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

x_1	= np.arange(-60, 20.01, 0.1)
y_1	= np.zeros(len(x_1)) - 34
y_2	= np.arange(-37, -30.99, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(-37, -30.99, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-100, 2500], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(4400 - time[2199:], transport[2199:], '-r', linewidth = 0.5, label = '1 - 2200', alpha = 0.3)
ax.plot(time[:2200], transport[:2200], '-k', linewidth = 0.5, label = '2201 - 4400', alpha = 0.3)

ax.plot(time_0600, transport_0600, '-b', linewidth = 0.5)
ax.plot(time_1500, transport_1500, '-b', linewidth = 0.5)
ax.plot(4400 - time_2900, transport_2900, '-', color = 'royalblue', linewidth = 0.5)
ax.plot(4400 - time_3800, transport_3800, '-', color = 'royalblue', linewidth = 0.5)

ax.quiver(600, 18, 3, 0, scale = 40, color = 'b', zorder = 10, width = 0.005)
ax.quiver(1500, 15, 3, 0, scale = 40, color = 'b', zorder = 10, width = 0.005)
ax.quiver(1500, 3, -3, 0, scale = 40, color = 'royalblue', zorder = 10, width = 0.005)
ax.quiver(600, 7, -3, 0, scale = 40, color = 'royalblue', zorder = 10, width = 0.005)

ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-2, 35)
ax.grid()

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = '1 - 2200')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = '2201 - 4400')
graph_3		= ax.plot([-100, -100], [-100, -100], '-b', linewidth = 1.5, label = 'Branched')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc=(0.315, 0.19), ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

ax.set_title('c) AMOC strength at 26$^{\circ}$N, branched simulations')

#-----------------------------------------------------------------------------------------

box1 = TextArea("1/", textprops=dict(color="k"))
box2 = TextArea("4400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.023, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("500/", textprops=dict(color="k"))
box2 = TextArea("3900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.165, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("1000/", textprops=dict(color="k"))
box2 = TextArea("3400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.379, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("1500/", textprops=dict(color="k"))
box2 = TextArea("2900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.605, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------

box1 = TextArea("2000/", textprops=dict(color="k"))
box2 = TextArea("2400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.835, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)

#-----------------------------------------------------------------------------------------
print()
print('Branch 0600')
print('Steady state:', np.mean(transport_0600[-50:]), 'Sv')
print('Hysteresis  :', np.mean(transport[574:624]), 'Sv')
print()

print('Branch 1500')
print('Steady state:', np.mean(transport_1500[-50:]), 'Sv')
print('Hysteresis  :', np.mean(transport[1474:1524]), 'Sv')
print()

print('Branch 2900')
print('Steady state:', np.mean(transport_2900[-50:]), 'Sv')
print('Hysteresis  :', np.mean(transport[2874:2924]), 'Sv')
print()

print('Branch 3800')
print('Steady state:', np.mean(transport_3800[-50:]), 'Sv')
print('Hysteresis  :', np.mean(transport[3774:3824]), 'Sv')
print()

#-----------------------------------------------------------------------------------------
time_var		= np.arange(50)

#Determine the detrended time series
trend, base 	= polyfit(time_var, transport_0600[-50:], 1)
print('AMOC variance branch 0600:', np.var(transport_0600[-50:] - ((trend * time_var) + base)), 'Sv^2')

trend, base 	= polyfit(time_var, transport_1500[-50:], 1)
print('AMOC variance branch 1500:', np.var(transport_1500[-50:] - ((trend * time_var) + base)), 'Sv^2')



show()