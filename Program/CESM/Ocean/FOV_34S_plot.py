#Program plots the FOV hysteresis

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_OV'][:]	#Fresh water

	fh.close()

	return time, FOV

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time, FOV_34S	        = ReadinData(directory+'Ocean/FOV_index_section_34S.nc')
time_0600, FOV_34S_0600	= ReadinData(directory+'Ocean/FOV_index_section_34S_branch_0600.nc')	
time_1500, FOV_34S_1500	= ReadinData(directory+'Ocean/FOV_index_section_34S_branch_1500.nc')	
time_2900, FOV_34S_2900	= ReadinData(directory+'Ocean/FOV_index_section_34S_branch_2900.nc')	
time_3800, FOV_34S_3800	= ReadinData(directory+'Ocean/FOV_index_section_34S_branch_3800.nc')	
#----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-100, 2500], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.plot([0.1/0.0003, 0.1/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.2/0.0003, 0.2/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.3/0.0003, 0.3/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.4/0.0003, 0.4/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.5/0.0003, 0.5/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.6/0.0003, 0.6/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)

ax.text(0.1/0.0003, 0.312, '0.1 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.2/0.0003, 0.312, '0.2 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.3/0.0003, 0.312, '0.3 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.4/0.0003, 0.312, '0.4 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.5/0.0003, 0.312, '0.5 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.6/0.0003, 0.312, '0.6 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)

ax.plot(4400 - time, FOV_34S, '-r', linewidth = 0.5)
ax.plot(time, FOV_34S, '-k', linewidth = 0.5)

y_error_0600		= np.zeros((2, 1)) + 1
y_error_0600[:, 0]	= np.mean(FOV_34S_0600[-50:]) - np.min(FOV_34S_0600[-50:]), np.max(FOV_34S_0600[-50:]) - np.mean(FOV_34S_0600[-50:])
ax.errorbar(600, np.mean(FOV_34S_0600[-50:]), color = 'b', marker = 's', markerfacecolor = 'dodgerblue', markeredgecolor = 'royalblue', yerr = y_error_0600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1500		= np.zeros((2, 1)) + 1
y_error_1500[:, 0]	= np.mean(FOV_34S_1500[-50:]) - np.min(FOV_34S_1500[-50:]), np.max(FOV_34S_1500[-50:]) - np.mean(FOV_34S_1500[-50:])
graph_3 		= ax.errorbar(1500, np.mean(FOV_34S_1500[-50:]), color = 'b', marker = 's', markerfacecolor = 'dodgerblue', markeredgecolor = 'royalblue', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')

y_error_2900		= np.zeros((2, 1)) + 1
y_error_2900[:, 0]	= np.mean(FOV_34S_2900[-50:]) - np.min(FOV_34S_2900[-50:]), np.max(FOV_34S_2900[-50:]) - np.mean(FOV_34S_2900[-50:])
ax.errorbar(4400 - 2900, np.mean(FOV_34S_2900[-50:]), color = 'b', marker = 's', markerfacecolor = 'dodgerblue', markeredgecolor = 'royalblue', yerr = y_error_2900, linewidth = 2.0, capsize=4)

y_error_3800		= np.zeros((2, 1)) + 1
y_error_3800[:, 0]	= np.mean(FOV_34S_3800[-50:]) - np.min(FOV_34S_3800[-50:]), np.max(FOV_34S_3800[-50:]) - np.mean(FOV_34S_3800[-50:])
ax.errorbar(4400 - 3800, np.mean(FOV_34S_3800[-50:]), color = 'b', marker = 's', markerfacecolor = 'dodgerblue', markeredgecolor = 'royalblue', yerr = y_error_3800, linewidth = 0.0, elinewidth = 2.0, capsize=4)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = '1 - 2200')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = '2201 - 4400')

legend_1	= ax.legend(loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.35, 0.35)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

#-----------------------------------------------------------------------------------------

box1 = TextArea("1/", textprops=dict(color="k"))
box2 = TextArea("4400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

#anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.049, -0.063), bbox_transform=ax.transAxes, borderpad=0.)
anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.023, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("500/", textprops=dict(color="k"))
box2 = TextArea("3900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

#anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.181, -0.063), bbox_transform=ax.transAxes, borderpad=0.)
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

ax.set_title('b) $F_{\mathrm{ovS}}$')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-100, 2500], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(4400 - time[2199:], FOV_34S[2199:], '-r', linewidth = 0.5, alpha = 0.3)
ax.plot(time[:2200], FOV_34S[:2200], '-k', linewidth = 0.5, alpha = 0.3)

ax.plot(4400 - time_2900, FOV_34S_2900, '-', color = 'royalblue', linewidth = 0.5)
ax.plot(4400 - time_3800, FOV_34S_3800, '-', color = 'royalblue',  linewidth = 0.5)
ax.plot(time_0600, FOV_34S_0600, '-b', linewidth = 0.5)
ax.plot(time_1500, FOV_34S_1500, '-b', linewidth = 0.5)

ax.quiver(600, 0.15, 3, 0, scale = 40, color = 'b', zorder = 10, width = 0.005)
ax.quiver(1500, -0.06, 3, 0, scale = 40, color = 'b', zorder = 10, width = 0.005)
ax.quiver(1500, 0.22, -3, 0, scale = 40, color = 'royalblue', zorder = 10, width = 0.005)
ax.quiver(600, -0.01, -3, 0, scale = 40, color = 'royalblue', zorder = 10, width = 0.005)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = '1 - 2200')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = '2201 - 4400')
graph_3		= ax.plot([-100, -100], [-100, -100], '-b', linewidth = 1.5, label = 'Branched')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.35, 0.35)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xticklabels(['1', '500/3900', '1000/3400', '1500/2900', '2000/2400'])
ax.tick_params(axis='x', colors='white')

#-----------------------------------------------------------------------------------------

box1 = TextArea("1/", textprops=dict(color="k"))
box2 = TextArea("4400 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

#anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.049, -0.063), bbox_transform=ax.transAxes, borderpad=0.)
anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(-0.023, -0.063), bbox_transform=ax.transAxes, borderpad=0.)

ax.add_artist(anchored_box)
#-----------------------------------------------------------------------------------------

box1 = TextArea("500/", textprops=dict(color="k"))
box2 = TextArea("3900 ", textprops=dict(color="r"))

box = HPacker(children=[box1, box2], align="center", pad=0, sep=0)

#anchored_box = AnchoredOffsetbox(loc=3, child=box, pad=0., frameon=False, bbox_to_anchor=(0.181, -0.063), bbox_transform=ax.transAxes, borderpad=0.)
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

ax.set_title('d) $F_{\mathrm{ovS}}$, branched simulations')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
print()
print('Branch 0600')
print('Steady state:', np.mean(FOV_34S_0600[-50:]), 'Sv')
print('Hysteresis  :', np.mean(FOV_34S[574:624]), 'Sv')
print()

print('Branch 1500')
print('Steady state:', np.mean(FOV_34S_1500[-50:]), 'Sv')
print('Hysteresis  :', np.mean(FOV_34S[1474:1524]), 'Sv')
print()

print('Branch 2900')
print('Steady state:', np.mean(FOV_34S_2900[-50:]), 'Sv')
print('Hysteresis  :', np.mean(FOV_34S[2874:2924]), 'Sv')
print()

print('Branch 3800')
print('Steady state:', np.mean(FOV_34S_3800[-50:]), 'Sv')
print('Hysteresis  :', np.mean(FOV_34S[3774:3824]), 'Sv')
print()

#-----------------------------------------------------------------------------------------
time_var		= np.arange(50)

#Determine the detrended time series
trend, base 	= polyfit(time_var, FOV_34S_0600[-50:], 1)
print('AMOC variance branch 0600:', np.var(FOV_34S_0600[-50:] - ((trend * time_var) + base)), 'Sv^2')

trend, base 	= polyfit(time_var, FOV_34S_1500[-50:], 1)
print('AMOC variance branch 1500:', np.var(FOV_34S_1500[-50:] - ((trend * time_var) + base)), 'Sv^2')

show()

