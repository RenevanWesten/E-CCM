#Plots the steady states for the E-CCM (temperature only)

from pylab import *
import numpy as np
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

#Making pathway to folder with all data
directory 	= '../../Data/E-CCM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Continuation/Continuation_MOC_box_model_temp_*.nc')
files.sort()

lambda_all  = np.zeros(len(files))
T_tropic    = np.zeros(len(files))
T_polar     = np.zeros(len(files))
T_n_all     = np.zeros((len(files), 2))
T_ts_all    = np.zeros((len(files), 2))
T_t_all     = np.zeros((len(files), 2))
T_s_all     = np.zeros((len(files), 2))
tipping_1   = np.zeros(len(files))
tipping_2   = np.zeros(len(files))

for file_i in range(len(files)):
    
	#Get the lambda value
	lambda_all[file_i]   = float('0.'+str(files[file_i][-10:-3])) * 10**6.0

	#Get the steady state and use these as input parameters, but these are freely allowed to evolve
	T_ts_a, T_n_a, S_t, S_n, S_ts, S_s, T_t, T_n, T_ts, T_s, T_d, D = np.loadtxt(directory+'Initalisation/MOC_box_model_lambda_a_0_'+files[file_i][-10:-3]+'.txt')

	T_tropic[file_i]      	= T_ts_a
	T_polar[file_i]    	    = T_n_a
	T_n_all[file_i, 0]   	= T_n
	T_ts_all[file_i, 0]   	= T_ts
	T_s_all[file_i, 0]   	= T_s
	T_t_all[file_i, 0]   	= T_t

	fh = netcdf.Dataset(files[file_i], 'r')

	E_A	    = fh.variables['E_A'][:]     	  	
	q_N	    = fh.variables['q_n'][:]     	  	
	F_ov    = fh.variables['F_ovS'][:]   
	T_t     = fh.variables['T_t'][:]     	  
	T_n     = fh.variables['T_n'][:]     	  
	T_ts    = fh.variables['T_ts'][:]     	  
	T_s     = fh.variables['T_s'][:]     	

	fh.close()

	E_A_1, q_N_1, F_ov_1	= E_A[0], q_N[0], F_ov[0]
	E_A_2, q_N_2, F_ov_2	= E_A[1], q_N[1], F_ov[1]
	E_A_3, q_N_3, F_ov_3	= E_A[2], q_N[2], F_ov[2]

	tipping_1[file_i]       = np.max(E_A_1)
	tipping_2[file_i]       = np.min(E_A_2)

#-----------------------------------------------------------------------------------------

fig, ax = subplots()

graph_T_tropic = ax.plot(lambda_all, T_tropic, '--', color = 'firebrick', linewidth = 2.0, label = '$T_{\mathrm{ts}}^a$')
graph_T_polar  = ax.plot(lambda_all, T_polar, '--', color = 'gray', linewidth = 2.0, label = '$T_{\mathrm{n}}^a$')
graph_T_ts     = ax.plot(lambda_all, T_ts_all[:, 0], '-', color = 'k', linewidth = 2.0, label = '$T_{\mathrm{ts}}$')
graph_T_n      = ax.plot(lambda_all, T_n_all[:, 0], '-', color = 'b', linewidth = 2.0, label = '$T_{\mathrm{n}}$')
graph_T_t     = ax.plot(lambda_all, T_t_all[:, 0], '-', color = 'r', linewidth = 2.0, label = '$T_{\mathrm{t}}$')
graph_T_s     = ax.plot(lambda_all, T_s_all[:, 0], '-', color = 'c', linewidth = 2.0, label = '$T_{\mathrm{s}}$')
ax.set_xlim(2.7, 10.3)
ax.set_ylim(-5, 20)
ax.grid()

ax.set_xlabel(r'$\lambda^a$ ($\times$ $10^{-6}$ m s$^{-1}$)')
ax.set_ylabel('Temperature ($^{\circ}$C)')

graphs	      = graph_T_tropic + graph_T_polar + graph_T_t + graph_T_ts + graph_T_n + graph_T_s

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=2, framealpha = 1.0, numpoints = 1)

ax.set_title('b) Steady states for $E_A = 0.0644$ Sv')

#-----------------------------------------------------------------------------------------

#The tipping points for the CC-model
tipping_1_CC_model  = 0.3532
tipping_2_CC_model  = 0.0598

#-----------------------------------------------------------------------------------------
    
fig, ax = subplots()

ax.plot([tipping_1_CC_model, tipping_1_CC_model], [-1, 11], '--', color = 'firebrick', linewidth = 1.5)
ax.plot([tipping_2_CC_model, tipping_2_CC_model], [-1, 11], '--', color = 'royalblue', linewidth = 1.5)

graph_tipping_1 = ax.plot(tipping_1, lambda_all, '-r', linewidth = 2.0, label = '$E_{A}^1$')
graph_tipping_2  = ax.plot(tipping_2, lambda_all, '-b', linewidth = 2.0, label = '$E_{A}^2$')

ax.set_xlim(0, 0.60)
ax.set_ylim(2.7, 10.3)
ax.grid()

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel(r'$\lambda^a$ ($\times 10^{-6}$ m s$^{-1}$)')

graphs	      = graph_tipping_1 + graph_tipping_2

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) Bifurcation points')

#-----------------------------------------------------------------------------------------

fig1, ax1 = subplots()
fig2, ax2 = subplots()
fig3, ax3 = subplots()
fig4, ax4 = subplots()

fh = netcdf.Dataset(files[2], 'r')

E_A	    = fh.variables['E_A'][:]     	  	
q_N	    = fh.variables['q_n'][:]  
T_n     = fh.variables['T_n'][:]     	  
T_ts    = fh.variables['T_ts'][:]
F_ov    = fh.variables['F_ovS'][:]   
      	  	   	
fh.close()

E_A_1, q_N_1, F_ov_1, T_n_1, T_ts_1	= E_A[0], q_N[0], F_ov[0], T_n[0], T_ts[0]
E_A_2, q_N_2, F_ov_2, T_n_2, T_ts_2	= E_A[1], q_N[1], F_ov[1], T_n[1], T_ts[1]
E_A_3, q_N_3, F_ov_3, T_n_3, T_ts_3	= E_A[2], q_N[2], F_ov[2], T_n[2], T_ts[2]

ax1.plot(E_A_3, q_N_3, ':k', linewidth = 1.5)
ax1.plot(E_A_1, q_N_1, '-r', linewidth = 2.0)
ax1.plot(E_A_2, q_N_2, '-b', linewidth = 2.0)
ax1.scatter(E_A_1[np.argmax(E_A_1)], q_N_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)
ax1.scatter(E_A_2[np.argmin(E_A_2)], q_N_2[np.argmin(E_A_2)], s = 30, color = 'b', edgecolor = 'dodgerblue', zorder = 10)

ax1.text(E_A_1[np.argmax(E_A_1)]+0.01, q_N_1[np.argmax(E_A_1)], '$E_A^1$', verticalalignment='center', horizontalalignment='left', color = 'red', fontsize = 10)
ax1.text(E_A_2[np.argmin(E_A_2)], q_N_2[np.argmin(E_A_2)]+0.3, '$E_A^2$', verticalalignment='bottom', horizontalalignment='center', color = 'b', fontsize = 10)

ax2.plot(E_A_3, T_n_3, ':k', linewidth = 1.5)
ax2.plot(E_A_1, T_n_1, '-r', linewidth = 2.0)
ax2.plot(E_A_2, T_n_2, '-b', linewidth = 2.0)
ax2.scatter(E_A_1[np.argmax(E_A_1)], T_n_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)
ax2.scatter(E_A_2[np.argmin(E_A_2)], T_n_2[np.argmin(E_A_2)], s = 30, color = 'b', edgecolor = 'dodgerblue', zorder = 10)

ax3.plot(E_A_3, T_ts_3, ':k', linewidth = 1.5)
ax3.plot(E_A_1, T_ts_1, '-r', linewidth = 2.0)
ax3.plot(E_A_2, T_ts_2, '-b', linewidth = 2.0)
ax3.scatter(E_A_1[np.argmax(E_A_1)], T_ts_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)
ax3.scatter(E_A_2[np.argmin(E_A_2)], T_ts_2[np.argmin(E_A_2)], s = 30, color = 'b', edgecolor = 'dodgerblue', zorder = 10)

ax4.plot(E_A_3, F_ov_3, ':k', linewidth = 1.5)
ax4.plot(E_A_1, F_ov_1, '-r', linewidth = 2.0)
ax4.plot(E_A_2, F_ov_2, '-b', linewidth = 2.0)
ax4.scatter(E_A_1[np.argmax(E_A_1)], F_ov_1[np.argmax(E_A_1)], s = 30, color = 'r', edgecolor = 'darkred', zorder = 10)
ax4.scatter(E_A_2[np.argmin(E_A_2)], F_ov_2[np.argmin(E_A_2)], s = 30, color = 'b', edgecolor = 'dodgerblue', zorder = 10)
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(files[14], 'r')

E_A	    = fh.variables['E_A'][:]     	  	
q_N	    = fh.variables['q_n'][:]     	  	   	
T_n     = fh.variables['T_n'][:]     	  
T_ts    = fh.variables['T_ts'][:]
F_ov    = fh.variables['F_ovS'][:] 

fh.close()

E_A_1, q_N_1, F_ov_1, T_n_1, T_ts_1	= E_A[0], q_N[0], F_ov[0], T_n[0], T_ts[0]
E_A_2, q_N_2, F_ov_2, T_n_2, T_ts_2	= E_A[1], q_N[1], F_ov[1], T_n[1], T_ts[1]
E_A_3, q_N_3, F_ov_3, T_n_3, T_ts_3	= E_A[2], q_N[2], F_ov[2], T_n[2], T_ts[2]

ax1.plot(E_A_3, q_N_3, linestyle = 'dashdot', color = 'gray', linewidth = 1.5)
ax1.plot(E_A_1, q_N_1, '--', color = 'firebrick', linewidth = 2.0)
ax1.plot(E_A_2, q_N_2, '--', color = 'cyan', linewidth = 2.0)
ax1.scatter(E_A_1[np.argmax(E_A_1)], q_N_1[np.argmax(E_A_1)], s = 30, color = 'firebrick', edgecolor = 'orangered', zorder = 10)
ax1.scatter(E_A_2[np.argmin(E_A_2)], q_N_2[np.argmin(E_A_2)], s = 30, color = 'cyan', edgecolor = 'darkturquoise', zorder = 10)

ax1.text(E_A_1[np.argmax(E_A_1)]+0.01, q_N_1[np.argmax(E_A_1)], '$E_A^1$', verticalalignment='center', horizontalalignment='left', color = 'firebrick', fontsize = 10)
ax1.text(E_A_2[np.argmin(E_A_2)], q_N_2[np.argmin(E_A_2)]+0.3, '$E_A^2$', verticalalignment='bottom', horizontalalignment='center', color = 'cyan', fontsize = 10)

ax2.plot(E_A_3, T_n_3, linestyle = 'dashdot', color = 'gray', linewidth = 1.5)
ax2.plot(E_A_1, T_n_1, '--', color = 'firebrick', linewidth = 2.0)
ax2.plot(E_A_2, T_n_2, '--', color = 'cyan', linewidth = 2.0)
ax2.scatter(E_A_1[np.argmax(E_A_1)], T_n_1[np.argmax(E_A_1)], s = 30, color = 'firebrick', edgecolor = 'orangered', zorder = 10)
ax2.scatter(E_A_2[np.argmin(E_A_2)], T_n_2[np.argmin(E_A_2)], s = 30, color = 'cyan', edgecolor = 'darkturquoise', zorder = 10)

ax3.plot(E_A_3, T_ts_3, linestyle = 'dashdot', color = 'gray', linewidth = 1.5)
ax3.plot(E_A_1, T_ts_1, '--', color = 'firebrick', linewidth = 2.0)
ax3.plot(E_A_2, T_ts_2, '--', color = 'cyan', linewidth = 2.0)
ax3.scatter(E_A_1[np.argmax(E_A_1)], T_ts_1[np.argmax(E_A_1)], s = 30, color = 'firebrick', edgecolor = 'orangered', zorder = 10)
ax3.scatter(E_A_2[np.argmin(E_A_2)], T_ts_2[np.argmin(E_A_2)], s = 30, color = 'cyan', edgecolor = 'darkturquoise', zorder = 10)

ax4.plot(E_A_3, F_ov_3, linestyle = 'dashdot', color = 'gray', linewidth = 1.5)
ax4.plot(E_A_1, F_ov_1, '--', color = 'firebrick', linewidth = 2.0)
ax4.plot(E_A_2, F_ov_2, '--', color = 'cyan', linewidth = 2.0)
ax4.scatter(E_A_1[np.argmax(E_A_1)], F_ov_1[np.argmax(E_A_1)], s = 30, color = 'firebrick', edgecolor = 'orangered', zorder = 10)
ax4.scatter(E_A_2[np.argmin(E_A_2)], F_ov_2[np.argmin(E_A_2)], s = 30, color = 'cyan', edgecolor = 'darkturquoise', zorder = 10)

#-----------------------------------------------------------------------------------------

ax1.set_xlim(0, 0.60)
ax1.set_ylim(-0.5, 20)
ax1.set_xlabel('Freshwater forcing (Sv)')
ax1.set_ylabel('Volume transport (Sv)')
ax1.set_yticks([0, 5, 10, 15, 20])
ax1.grid()


graph_tipping_1 = ax1.plot(E_A_1, q_N_1-100, '-k', linewidth = 2.0, label = r'$\lambda^a = 3.5 \times 10^{-6}$')
graph_tipping_2 = ax1.plot(E_A_1, q_N_1-100, '--k', linewidth = 2.0, label = r'$\lambda^a = 9.5 \times 10^{-6}$')

graphs_1	= graph_tipping_1 + graph_tipping_2
legend_labels_1 = [l.get_label() for l in graphs_1]
legend_1	= ax1.legend(graphs_1, legend_labels_1, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

graph_AMOC_on  = ax1.plot(E_A_1, q_N_1-100, '-r', linewidth = 2.0, label = 'AMOC on')
graph_AMOC_off = ax1.plot(E_A_1, q_N_1-100, '-b', linewidth = 2.0, label = 'AMOC off')

graphs_2	= graph_AMOC_on + graph_AMOC_off

legend_labels_2 = [l.get_label() for l in graphs_2]
legend_2	= ax1.legend(graphs_2, legend_labels_2, loc = 'center left', ncol=1, framealpha = 1.0)
ax1.add_artist(legend_1)

ax1.set_title('a) AMOC strength ($q_N$)')

#-----------------------------------------------------------------------------------------

ax2.set_xlim(0, 0.60)
ax2.set_ylim(-1, 6)
ax2.set_xlabel('Freshwater forcing (Sv)')
ax2.set_ylabel('Temperature ($^{\circ}$C)')
ax2.grid()

legend_1	= ax2.legend(graphs_1, legend_labels_1, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)
ax2.legend(graphs_2, legend_labels_2, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)
ax2.add_artist(legend_1)

ax2.set_title('b) Temperature box n ($T_\mathrm{n}$)')

#-----------------------------------------------------------------------------------------

ax3.set_xlim(0, 0.60)
ax3.set_ylim(9, 12)
ax3.set_xlabel('Freshwater forcing (Sv)')
ax3.set_ylabel('Temperature ($^{\circ}$C)')
ax3.grid()

legend_1	= ax3.legend(graphs_1, legend_labels_1, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)
ax3.legend(graphs_2, legend_labels_2, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)
ax3.add_artist(legend_1)

ax3.set_title('Temperature ts-box ($T_{ts}$)')

#-----------------------------------------------------------------------------------------

ax4.set_xlim(0, 0.60)
ax4.set_ylim(-0.4, 0.1)
ax4.set_xlabel('Freshwater forcing (Sv)')
ax4.set_ylabel('Freshwater transport (Sv)')
ax4.grid()

legend_1	= ax4.legend(graphs_1, legend_labels_1, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)
ax4.legend(graphs_2, legend_labels_2, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)
ax4.add_artist(legend_1)

ax4.set_title('c) $F_{\mathrm{ovS}}$')

show()



    

