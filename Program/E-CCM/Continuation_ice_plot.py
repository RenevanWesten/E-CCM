#Plots the steady states for the E-CCM (temperature only)

from pylab import *
import numpy as np
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory 	= '../../Data/E-CCM/'
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Continuation/Continuation_MOC_box_model_ice_lambda_a_0_0000035_r_min_0_*.nc')
files.sort()

q_min_all   = np.zeros(len(files))
tipping_1   = np.zeros(len(files))
tipping_2   = np.zeros(len(files))
tipping_3   = ma.masked_all(len(files))

for file_i in range(len(files)):

	#Get the q min factor
	q_min_all[file_i]    = float('0.'+str(files[file_i][-6:-3]))

	fh = netcdf.Dataset(files[file_i], 'r')

	E_A	    = fh.variables['E_A'][:]     	  	
	q_N	    = fh.variables['q_n'][:]     	  	
	F_ov    = fh.variables['F_ovS'][:]   
	T_t     = fh.variables['T_t'][:]     	  
	T_n     = fh.variables['T_n'][:]     	  
	T_ts    = fh.variables['T_ts'][:]     	  
	T_s     = fh.variables['T_s'][:]     	

	fh.close()

	E_A_1, q_N_1, F_ov_1, T_n_1, T_ts_1	= E_A[0], q_N[0], F_ov[0], T_n[0], T_ts[0]	#AMOC on
	E_A_2, q_N_2, F_ov_2, T_n_2, T_ts_2	= E_A[1], q_N[1], F_ov[1], T_n[1], T_ts[1]	#AMOC weak
	E_A_8, q_N_8, F_ov_8, T_n_8, T_ts_8	= E_A[7], q_N[7], F_ov[7], T_n[7], T_ts[7]	#AMOC off

	tipping_1[file_i]       = np.max(E_A_1)
	tipping_2[file_i]       = np.min(E_A_8)

	if np.any(E_A_2.mask == False):
		#Third branch possible
		tipping_3[file_i]       = np.min(E_A_2)

#-----------------------------------------------------------------------------------------
fig, ax = subplots()

graph_tipping_1 = ax.plot(tipping_1, q_min_all, '-r', linewidth = 2.0, label = '$E_{A}^1$')
graph_tipping_2  = ax.plot(tipping_2, q_min_all, '-b', linewidth = 2.0, label = '$E_{A}^2$')
graph_tipping_3  = ax.plot(tipping_3, q_min_all, '-k', linewidth = 2.0, label = '$E_{A}^3$')

ax.set_xlim(0, 0.60)
ax.set_ylim(0.018, 0.102)
ax.grid()

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel(r'$r_{q}^{\mathrm{min}}$')

graphs	      = graph_tipping_1 + graph_tipping_2 + graph_tipping_3

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)


ax.quiver(0.167, 0.038, 0.75, 0, scale = 5, color = 'k', zorder = 10, width = 0.005)
ax.quiver(0.21, 0.054, 0.44, 0.25, scale = 5, color = 'k', zorder = 10, width = 0.005)
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Continuation/Continuation_MOC_box_model_ice_lambda_a_0_0000035_r_min_0_038.nc', 'r')

#Writing data to correct variable	
E_A	= fh.variables['E_A'][:]     	  	
q_N	= fh.variables['q_n'][:]  

fh.close()	   	

E_A_1, q_N_1	= E_A[1], q_N[1]	#Weak state 1
E_A_2, q_N_2	= E_A[2], q_N[2]	#Weak state 2
E_A_3, q_N_3	= E_A[4], q_N[4]	#Unstable state 1
E_A_4, q_N_4	= E_A[5], q_N[5]	#Unstable state 2
E_A_5, q_N_5	= E_A[7], q_N[7]	#Off state

ax2 	= fig.add_axes([0.47, 0.18, 0.25, 0.20])
ax2.plot(E_A_1, q_N_1, '-k', linewidth = 1.5)
ax2.plot(E_A_2, q_N_2, '-k', linewidth = 1.5)
ax2.plot(E_A_3, q_N_3, ':k', linewidth = 1.5)
ax2.plot(E_A_4, q_N_4, ':k', linewidth = 1.5)
ax2.plot(E_A_5, q_N_5, '-b', linewidth = 1.5)


ax2.set_xlim(0.135, 0.215)
ax2.set_ylim(-0.2, 3.5)
ax2.set_yticks([0, 3])
ax2.set_xticks([0.15, 0.20])
ax2.grid()
ax2.set_title('$q_N^{\mathrm{ice}}$, $r_q^{\mathrm{min}}$ = 0.038', fontsize = 10)

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Continuation/Continuation_MOC_box_model_ice_lambda_a_0_0000035_r_min_0_054.nc', 'r')

#Writing data to correct variable	
E_A	= fh.variables['E_A'][:]     	  	
q_N	= fh.variables['q_n'][:]  

fh.close()	   	

E_A_1, q_N_1	= E_A[1], q_N[1]	#Weak state 1
E_A_2, q_N_2	= E_A[2], q_N[2]	#Weak state 2
E_A_3, q_N_3	= E_A[4], q_N[4]	#Unstable state 1
E_A_4, q_N_4	= E_A[5], q_N[5]	#Unstable state 2
E_A_5, q_N_5	= E_A[7], q_N[7]	#Off state

ax3 	= fig.add_axes([0.47, 0.50, 0.25, 0.20])
ax3.plot(E_A_1, q_N_1, '-k', linewidth = 1.5)
ax3.plot(E_A_2, q_N_2, '-k', linewidth = 1.5)
ax3.plot(E_A_3, q_N_3, ':k', linewidth = 1.5)
ax3.plot(E_A_4, q_N_4, ':k', linewidth = 1.5)
ax3.plot(E_A_5, q_N_5, '-b', linewidth = 1.5)


ax3.set_xlim(0.135, 0.215)
ax3.set_ylim(-0.2, 3.5)
ax3.set_yticks([0, 3])
ax3.set_xticks([0.15, 0.20])
ax3.grid()
ax3.set_title('$q_N^{\mathrm{ice}}$, $r_q^{\mathrm{min}}$ = 0.054', fontsize = 10)

ax.set_title('d) Bifurcation points')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig1, ax1 = subplots()
fig2, ax2 = subplots()
fig3, ax3 = subplots()
fig4, ax4 = subplots()

fh = netcdf.Dataset(files[1], 'r')

E_A	    = fh.variables['E_A'][:]     	  	
q_N	    = fh.variables['q_n'][:]  
T_n     = fh.variables['T_n'][:]     	  
T_ts    = fh.variables['T_ts'][:]
F_ov    = fh.variables['F_ovS'][:]  
ice	    = fh.variables['Ice'][:] 
  	   	
fh.close()

E_A_1, q_N_1, F_ov_1, T_n_1, T_ts_1, ice_1	= E_A[0], q_N[0], F_ov[0], T_n[0], T_ts[0], ice[0]	#On state
E_A_2, q_N_2, F_ov_2, T_n_2, T_ts_2, ice_2	= E_A[7], q_N[7], F_ov[7], T_n[7], T_ts[7], ice[7]	#Off state
E_A_3, q_N_3, F_ov_3, T_n_3, T_ts_3, ice_3	= E_A[1], q_N[1], F_ov[1], T_n[1], T_ts[1], ice[1]	#Weak state
E_A_4, q_N_4, F_ov_4, T_n_4, T_ts_4, ice_4	= E_A[4], q_N[4], F_ov[4], T_n[4], T_ts[4], ice[4]	#Unstable state

ax1.plot(E_A_4, q_N_4, ':k', linewidth = 1.5)
ax1.plot(E_A_1, q_N_1, '-r', linewidth = 2.0)
ax1.plot(E_A_2, q_N_2, '-b', linewidth = 2.0)
ax1.plot(E_A_3, q_N_3, '-k', linewidth = 2.0)
ax1.scatter(E_A_1[np.argmax(E_A_1)], q_N_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)
ax1.scatter(E_A_2[np.argmin(E_A_2)], q_N_2[np.argmin(E_A_2)], s = 30, marker = 'D', color = 'b', zorder = 10)
ax1.scatter(E_A_3[np.argmin(E_A_3)], q_N_3[np.argmin(E_A_3)], s = 30, color = 'k', zorder = 10)

ax1.text(E_A_1[np.argmax(E_A_1)]+0.01, q_N_1[np.argmax(E_A_1)], '$E_A^1$', verticalalignment='center', horizontalalignment='left', color = 'red', fontsize = 10)
ax1.text(E_A_2[np.argmin(E_A_2)], q_N_2[np.argmin(E_A_2)]+0.3, '$E_A^2$', verticalalignment='bottom', horizontalalignment='center', color = 'b', fontsize = 10)
ax1.text(E_A_3[np.argmin(E_A_3)], q_N_3[np.argmin(E_A_3)]-0.3, '$E_A^3$', verticalalignment='top', horizontalalignment='center', color = 'k', fontsize = 10)

ax2.plot(E_A_4, T_n_4, ':k', linewidth = 1.5)
ax2.plot(E_A_1, T_n_1, '-r', linewidth = 2.0)
ax2.plot(E_A_2, T_n_2, '-b', linewidth = 2.0)
ax2.plot(E_A_3, T_n_3, '-k', linewidth = 2.0)
ax2.scatter(E_A_1[np.argmax(E_A_1)], T_n_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)
ax2.scatter(E_A_2[np.argmin(E_A_2)], T_n_2[np.argmin(E_A_2)], s = 30, marker = 'D', color = 'b', edgecolor = 'dodgerblue', zorder = 10)
ax2.scatter(E_A_3[np.argmin(E_A_3)], T_n_3[np.argmin(E_A_3)], s = 30, color = 'k', edgecolor = 'gray', zorder = 10)

ax2_2 	= ax2.twinx()
ax2_2.plot(E_A_4, ice_4, linestyle = 'dashdot', color = 'gray', linewidth = 1.5)
ax2_2.plot(E_A_1, ice_1, '--', color = 'firebrick', linewidth = 2.0)
ax2_2.plot(E_A_2, ice_2, '--', color = 'cyan', linewidth = 2.0)
ax2_2.plot(E_A_3, ice_3, '--', color = 'grey', linewidth = 2.0)
ax2_2.scatter(E_A_1[np.argmax(E_A_1)], ice_1[np.argmax(E_A_1)], s = 30, color = 'firebrick', edgecolor = 'orangered', zorder = 10)
ax2_2.scatter(E_A_2[np.argmin(E_A_2)], ice_2[np.argmin(E_A_2)], s = 30, marker = 'D', color = 'cyan', edgecolor = 'darkturquoise', zorder = 10)
ax2_2.scatter(E_A_3[np.argmin(E_A_3)], ice_3[np.argmin(E_A_3)], s = 30, color = 'grey', edgecolor = 'dimgrey', zorder = 10)

ax3.plot(E_A_4, T_ts_4, ':k', linewidth = 1.5)
ax3.plot(E_A_1, T_ts_1, '-r', linewidth = 2.0)
ax3.plot(E_A_2, T_ts_2, '-b', linewidth = 2.0)
ax3.plot(E_A_3, T_ts_3, '-k', linewidth = 2.0)
ax3.scatter(E_A_1[np.argmax(E_A_1)], T_ts_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)
ax3.scatter(E_A_2[np.argmin(E_A_2)], T_ts_2[np.argmin(E_A_2)], s = 30, marker = 'D', color = 'b', edgecolor = 'dodgerblue', zorder = 10)
ax3.scatter(E_A_3[np.argmin(E_A_3)], T_ts_3[np.argmin(E_A_3)], s = 30, color = 'k', edgecolor = 'gray', zorder = 10)

ax4.plot(E_A_4, F_ov_4, ':k', linewidth = 1.5)
ax4.plot(E_A_1, F_ov_1, '-r', linewidth = 2.0)
ax4.plot(E_A_2, F_ov_2, '-b', linewidth = 2.0)
ax4.plot(E_A_3, F_ov_3, '-k', linewidth = 2.0)
ax4.scatter(E_A_1[np.argmax(E_A_1)], F_ov_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)
ax4.scatter(E_A_2[np.argmin(E_A_2)], F_ov_2[np.argmin(E_A_2)], s = 30, marker = 'D', color = 'b', edgecolor = 'dodgerblue', zorder = 10)
ax4.scatter(E_A_3[np.argmin(E_A_3)], F_ov_3[np.argmin(E_A_3)], s = 30, color = 'k', edgecolor = 'gray', zorder = 10)
#-----------------------------------------------------------------------------------------

ax1.set_xlim(0, 0.60)
ax1.set_ylim(-0.5, 20)
ax1.set_xlabel('Freshwater forcing (Sv)')
ax1.set_ylabel('Volume transport (Sv)')
ax1.set_yticks([0, 5, 10, 15, 20])
ax1.grid()

graph_AMOC_on   = ax1.plot(E_A_1, q_N_1-100, '-r', linewidth = 2.0, label = 'AMOC on')
graph_AMOC_off  = ax1.plot(E_A_1, q_N_1-100, '-b', linewidth = 2.0, label = 'AMOC off')
graph_AMOC_weak = ax1.plot(E_A_1, q_N_1-100, '-k', linewidth = 2.0, label = 'AMOC weak')

graphs		= graph_AMOC_on + graph_AMOC_off + graph_AMOC_weak
legend_labels   = [l.get_label() for l in graphs]
legend_1	= ax1.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax1.set_title('a) AMOC strength ($q_N^{\mathrm{ice}}$)')

#-----------------------------------------------------------------------------------------

ax2.set_xlim(0, 0.60)
ax2.set_ylim(-1, 6)
ax2.set_xlabel('Freshwater forcing (Sv)')
ax2.set_ylabel('Temperature ($^{\circ}$C)')
ax2.grid()

ax2_2.set_ylim(-5, 105)
ax2_2.set_ylabel('Sea-ice fraction ($\%$)')

legend_1	= ax2.legend(graphs, legend_labels, loc='center left', ncol=1, framealpha = 1.0, numpoints = 1)

graph_temp = ax2.plot(E_A_1, T_n_1-100, '-k', linewidth = 2.0, label = 'Temperature, $T_\mathrm{n}$')
graph_ice  = ax2.plot(E_A_1, T_n_1-100, '--k', linewidth = 2.0, label = 'Sea-ice fraction, $f_\mathrm{n}$')

graphs_2	  = graph_temp + graph_ice
legend_labels_2   = [l.get_label() for l in graphs_2]

ax2.legend(graphs_2, legend_labels_2, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)
ax2.add_artist(legend_1)

ax2.set_title('b) Temperature ($T_\mathrm{n}$) and sea-ice fraction ($f_\mathrm{n}$) box n')

#-----------------------------------------------------------------------------------------

ax3.set_xlim(0, 0.60)
ax3.set_ylim(9, 12)
ax3.set_xlabel('Freshwater forcing (Sv)')
ax3.set_ylabel('Temperature ($^{\circ}$C)')
ax3.grid()

legend_1	= ax3.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax3.set_title('Temperature ts-box ($T_{ts}$)')

#-----------------------------------------------------------------------------------------

ax4.set_xlim(0, 0.60)
ax4.set_ylim(-0.4, 0.1)
ax4.set_xlabel('Freshwater forcing (Sv)')
ax4.set_ylabel('Freshwater transport (Sv)')
ax4.grid()

legend_1	= ax4.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax4.set_title('c) $F_{\mathrm{ovS}}$')

ax4_2 	= fig4.add_axes([0.22, 0.165, 0.25, 0.20])

ax4_2.plot(np.insert(E_A_4, 0, E_A_1[np.argmax(E_A_1)]), np.insert(F_ov_4, 0, F_ov_1[np.argmax(E_A_1)]), ':k', linewidth = 1.5)
ax4_2.plot(E_A_1, F_ov_1, '-r', linewidth = 2.0)
ax4_2.scatter(E_A_1[np.argmax(E_A_1)], F_ov_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)


ax4_2.set_xlim(0.47, 0.49)
ax4_2.set_ylim(-0.362, -0.348)
ax4_2.set_yticks([-0.36, -0.35])
ax4_2.grid()
ax4_2.set_title('Around $F_{\mathrm{ovS}}$ minimum', fontsize = 10)

show()



    


