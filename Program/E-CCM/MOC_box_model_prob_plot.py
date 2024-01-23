#Plots the transition probabilties for the E-CCM

from pylab import *
import numpy as np
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory 	= '../../Data/E-CCM/'


def ReadinData(filename):
    
    fh = netcdf.Dataset(filename, 'r')
    
    E_A         =	fh.variables['E_A'][:] 
    f_noise     = 	fh.variables['f_noise'][:]      	
    prob        = 	fh.variables['prob'][:]     
        
    fh.close()
    
    return E_A, f_noise, prob

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------


#CCM
E_A, f_noise, prob_1        = ReadinData(directory+'TAMS/TAMS_start_on_nc_10_tmax_100.nc')
E_A, f_noise, prob_2        = ReadinData(directory+'TAMS/TAMS_start_off_nc_10_tmax_100.nc')
E_A, f_noise_cst, prob_3    = ReadinData(directory+'TAMS/TAMS_start_on_nc_10_tmax_100_noise_0_025.nc')
E_A, f_noise_cst, prob_4    = ReadinData(directory+'TAMS/TAMS_start_off_nc_10_tmax_100_noise_0_025.nc')

#E-CCM (temperature only)
E_A, f_noise, prob_5        = ReadinData(directory+'TAMS/TAMS_temp_start_on_nc_10_tmax_100.nc')
E_A, f_noise, prob_6        = ReadinData(directory+'TAMS/TAMS_temp_start_off_nc_10_tmax_100.nc')
E_A, f_noise_cst, prob_7    = ReadinData(directory+'TAMS/TAMS_temp_start_on_nc_10_tmax_100_noise_0_025.nc')
E_A, f_noise_cst, prob_8    = ReadinData(directory+'TAMS/TAMS_temp_start_off_nc_10_tmax_100_noise_0_025.nc')

#E-CCM
E_A, f_noise, prob_9        = ReadinData(directory+'TAMS/TAMS_temp_ice_start_on_nc_10_tmax_100.nc')
E_A, f_noise, prob_10       = ReadinData(directory+'TAMS/TAMS_temp_ice_start_off_nc_10_tmax_100.nc')
E_A, f_noise_cst, prob_11   = ReadinData(directory+'TAMS/TAMS_temp_ice_start_on_nc_10_tmax_100_noise_0_025.nc')
E_A, f_noise_cst, prob_12   = ReadinData(directory+'TAMS/TAMS_temp_ice_start_off_nc_10_tmax_100_noise_0_025.nc')
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-1, 1], y1 = np.zeros(2) - 1, y2 = np.zeros(2) + 1, color = 'gray', alpha = 0.25)

CS	= pcolor(E_A, f_noise, np.mean(prob_1, axis = 0), norm=colors.LogNorm(10**(-8), 10**0), cmap = 'Spectral_r', shading = 'auto')
cbar	= colorbar(CS, extend = 'min', ticks = 10**(np.arange(-8, 0.1, 1)))
cbar.set_label('Transition probability')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Noise amplitude ($f_{\sigma}$)')
ax.set_xlim(0, 0.5)
ax.set_ylim(0, 0.5)

ax.plot(E_A[1:], 0.025 / E_A[1:], ':k', linewidth = 1.0)

ax.set_title(r'a) On-to-off transitions (CCM)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-1, 1], y1 = np.zeros(2) - 1, y2 = np.zeros(2) + 1, color = 'gray', alpha = 0.25)

CS	= pcolor(E_A, f_noise, np.mean(prob_2, axis = 0), norm=colors.LogNorm(10**(-8), 10**0), cmap = 'Spectral_r', shading = 'auto')
cbar	= colorbar(CS, extend = 'min', ticks = 10**(np.arange(-8, 0.1, 1)))
cbar.set_label('Transition probability')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Noise amplitude ($f_{\sigma}$)')
ax.set_xlim(0, 0.5)
ax.set_ylim(0, 0.5)

ax.plot(E_A[1:], 0.025 / E_A[1:], ':k', linewidth = 1.0)

ax.set_title(r'b) Off-to-on transitions (CCM)')

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

ax.fill_between(E_A, y1 = np.min(prob_3, axis = 0), y2 = np.max(prob_3, axis = 0), color = 'red', alpha = 0.20)
ax.fill_between(E_A, y1 = np.min(prob_4, axis = 0), y2 = np.max(prob_4, axis = 0), color = 'blue', alpha = 0.20)

graph_1 = ax.plot(E_A, np.mean(prob_3, axis = 0), '-r', linewidth = 2.0, label = 'On-to-off transitions')
graph_2 = ax.plot(E_A, np.mean(prob_4, axis = 0), '-b', linewidth = 2.0, label = 'Off-to-on transitions')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Transition probability')
ax.set_xlim(0, 0.5)
ax.set_ylim(10**(-8.0), 10**(0.0))
ax.set_yscale('log')
ax.set_yticks(10**(np.arange(-8, 0.1, 1)))
ax.grid()

legend_1	= ax.legend(loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title(r'c) $f_{\sigma} E_A = 0.025$ Sv (CCM)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-1, 1], y1 = np.zeros(2) - 1, y2 = np.zeros(2) + 1, color = 'gray', alpha = 0.25)

CS	= pcolor(E_A, f_noise, np.mean(prob_5, axis = 0), norm=colors.LogNorm(10**(-8), 10**0), cmap = 'Spectral_r', shading = 'auto')
cbar	= colorbar(CS, extend = 'min', ticks = 10**(np.arange(-8, 0.1, 1)))
cbar.set_label('Transition probability')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Noise amplitude ($f_{\sigma}$)')
ax.set_xlim(0, 0.5)
ax.set_ylim(0, 0.5)

ax.plot(E_A[1:], 0.025 / E_A[1:], ':k', linewidth = 1.0)

ax.set_title(r'd) On-to-off transitions (E-CCM, temperature only)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-1, 1], y1 = np.zeros(2) - 1, y2 = np.zeros(2) + 1, color = 'gray', alpha = 0.25)

CS	= pcolor(E_A, f_noise, np.mean(prob_6, axis = 0), norm=colors.LogNorm(10**(-8), 10**0), cmap = 'Spectral_r', shading = 'auto')
cbar	= colorbar(CS, extend = 'min', ticks = 10**(np.arange(-8, 0.1, 1)))
cbar.set_label('Transition probability')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Noise amplitude ($f_{\sigma}$)')
ax.set_xlim(0, 0.5)
ax.set_ylim(0, 0.5)

ax.plot(E_A[1:], 0.025 / E_A[1:], ':k', linewidth = 1.0)

ax.set_title(r'e) Off-to-on transitions (E-CCM, temperature only)')

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

ax.fill_between(E_A, y1 = np.min(prob_7, axis = 0), y2 = np.max(prob_7, axis = 0), color = 'red', alpha = 0.20)
ax.fill_between(E_A, y1 = np.min(prob_8, axis = 0), y2 = np.max(prob_8, axis = 0), color = 'blue', alpha = 0.20)

graph_1 = ax.plot(E_A, np.mean(prob_7, axis = 0), '-r', linewidth = 2.0, label = 'On-to-off transitions')
graph_2 = ax.plot(E_A, np.mean(prob_8, axis = 0), '-b', linewidth = 2.0, label = 'Off-to-on transitions')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Transition probability')
ax.set_xlim(0, 0.5)
ax.set_ylim(10**(-8.0), 10**(0.0))
ax.set_yscale('log')
ax.set_yticks(10**(np.arange(-8, 0.1, 1)))
ax.grid()

legend_1	= ax.legend(loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title(r'f) $f_{\sigma} E_A = 0.025$ Sv (E-CCM, temperature only)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-1, 1], y1 = np.zeros(2) - 1, y2 = np.zeros(2) + 1, color = 'gray', alpha = 0.25)

CS	= pcolor(E_A, f_noise, np.mean(prob_9, axis = 0), norm=colors.LogNorm(10**(-8), 10**0), cmap = 'Spectral_r', shading = 'auto')
cbar	= colorbar(CS, extend = 'min', ticks = 10**(np.arange(-8, 0.1, 1)))
cbar.set_label('Transition probability')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Noise amplitude ($f_{\sigma}$)')
ax.set_xlim(0, 0.5)
ax.set_ylim(0, 0.5)

ax.plot(E_A[1:], 0.025 / E_A[1:], ':k', linewidth = 1.0)

ax.set_title(r'g) On-to-off transitions (E-CCM)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-1, 1], y1 = np.zeros(2) - 1, y2 = np.zeros(2) + 1, color = 'gray', alpha = 0.25)

CS	= pcolor(E_A, f_noise, np.mean(prob_10, axis = 0), norm=colors.LogNorm(10**(-8), 10**0), cmap = 'Spectral_r', shading = 'auto')
cbar	= colorbar(CS, extend = 'min', ticks = 10**(np.arange(-8, 0.1, 1)))
cbar.set_label('Transition probability')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Noise amplitude ($f_{\sigma}$)')
ax.set_xlim(0, 0.5)
ax.set_ylim(0, 0.5)

ax.plot(E_A[1:], 0.025 / E_A[1:], ':k', linewidth = 1.0)

ax.set_title(r'h) Off-to-on transitions (E-CCM)')

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

ax.fill_between(E_A, y1 = np.min(prob_11, axis = 0), y2 = np.max(prob_11, axis = 0), color = 'red', alpha = 0.20)
ax.fill_between(E_A, y1 = np.min(prob_12, axis = 0), y2 = np.max(prob_12, axis = 0), color = 'blue', alpha = 0.20)

graph_1 = ax.plot(E_A, np.mean(prob_11, axis = 0), '-r', linewidth = 2.0, label = 'On-to-off transitions')
graph_2 = ax.plot(E_A, np.mean(prob_12, axis = 0), '-b', linewidth = 2.0, label = 'Off-to-on transitions')

ax.set_xlabel('Freshwater forcing (Sv)')
ax.set_ylabel('Transition probability')
ax.set_xlim(0, 0.5)
ax.set_ylim(10**(-8.0), 10**(0.0))
ax.set_yscale('log')
ax.set_yticks(10**(np.arange(-8, 0.1, 1)))
ax.grid()

legend_1	= ax.legend(loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

#-----------------------------------------------------------------------------------------
ax2 	= fig.add_axes([0.50, 0.418, 0.42, 0.20])

prob_10	= ma.masked_where(prob_10 < 0, prob_10)

ax2.plot([0.16285, 0.16285], [10**(-3.0), 10], '--k', linewidth = 1.0)
ax2.plot(E_A, np.mean(prob_10[:, 60], axis = 0), '-b', linewidth = 2.0)	#Plot for f = 0.3

ax2.set_yscale('log')
ax2.set_ylim(10**(-2.0), 10**(0.0))
ax2.set_xlim(0, 0.5)
ax2.set_xlabel('Freshwater forcing (Sv)')
ax2.set_ylabel('$p$')
ax2.grid()

ax2.text(0.16285 - 0.01, 1.5 * 10**(-2.0), '$F_{\mathrm{ovS}} < 0$', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 10)
ax2.text(0.16285 + 0.01, 1.5 * 10**(-2.0), r'$F_{\mathrm{ovS}} \gtrapprox 0$', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize = 10)

ax2.set_title('Off/weak-to-on transitions, $f_{\sigma} = 0.3$', fontsize = 10)

#-----------------------------------------------------------------------------------------

ax.set_title(r'i) $f_{\sigma} E_A = 0.025$ Sv (E-CCM)')

show()