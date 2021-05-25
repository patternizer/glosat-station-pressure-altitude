#!/usr/bin/env python

#-----------------------------------------------------------------------
# PROGRAM: glosat-hpa-altitude-converter.py
#-----------------------------------------------------------------------
# Version 0.1
# 24 May, 2021
# Dr Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
# michael DOT a DOT taylor AT uea DOT ac DOT uk
#-----------------------------------------------------------------------

# Dataframe libraries:
import numpy as np
import pandas as pd
import xarray as xr

# Maths libraries:
import math
from scipy.interpolate import griddata
from scipy import spatial
import scipy
import scipy.stats as stats    

# Plotting libraries
import matplotlib.pyplot as plt
from matplotlib import image

#-----------------------------------------------------------------------------
# SETTINGS
#-----------------------------------------------------------------------------

g = 9.80665
R = 287.00
pressure_level_file = 'DATA/air_z_3hourly_BN.nc'
altitude_file = 'DATA/lut.csv'
image_file = 'IMAGES/ben-nevis-tranquility.jpg'
fontsize = 16

#-----------------------------------------------------------------------------
# METHODS
#-----------------------------------------------------------------------------

def altitude_to_isa(p0, t0, a, h0, h1):
    
	if a != 0:
		t1 = t0 + a * (h1 - h0)
		p1 = p0 * (t1 / t0) ** (-g / a / R)
	else:
		t1 = t0
		p1 = p0 * math.exp(-g / R / t0 * (h1 - h0))
	return t1, p1

def hpa_to_altitude(hpa):

    h = ( (1.0 - (hpa/1013.25)**0.190284) )*145366.45
    altitude = h*0.3048
    return altitude

def isa(altitude):
    
	a = [-0.0065, 0, 0.0010, 0.0028, 0, -0.0028, -0.0020]
	h = [11000, 20000, 32000, 47000, 51000, 71000, 85000]
	p0 = 101325
	t0 = 288.15
	prevh = 0
	if altitude < 0 or altitude > 85000:
		print("You are in space! Please choose a more Earthly value [0, 85000]")
		return
	for i in range(0, 7):
		if altitude <= h[i]:
			temperature, pressure = altitude_to_isa(p0, t0, a[i], prevh, altitude)
			break;
		else:
			t0, p0 = altitude_to_isa(p0, t0, a[i], prevh, h[i])
			prevh = h[i]
	density = pressure / (R * temperature)
	pressure = pressure/100.0
	temperature = temperature-273.15
	strformat = 'Temperature: {0:.3f} °C \nPressure: {1:.3f} hPa \nDensity: {2:.3f} kg/m³'
	print(strformat.format(temperature, pressure, density))    
	return pressure, temperature, density
	
def gamma_fit(y):
	
	mask = np.isfinite(y)
	z = y[mask] 
	n = len(z) 
	z_map = np.linspace(z.min(), z.max(), n)
	disttype = 'gamma' 
	dist = getattr(scipy.stats, disttype)
	param = dist.fit(z) 
	gamma_pdf = dist.pdf(z_map, param[0], param[1], param[2])	
	return z_map, gamma_pdf, param

#-----------------------------------------------------------------------------
# LOAD: 20CRv3 pressure level data
#-----------------------------------------------------------------------------

df = xr.open_dataset(pressure_level_file, decode_cf=True)
levels = np.array(df.level)

#-----------------------------------------------------------------------------
# LOAD: look-up table
#-----------------------------------------------------------------------------

lut = pd.read_csv(altitude_file, index_col=0)
station_altitudes = lut.station_elevation.values

#-----------------------------------------------------------------------------
# CREATE: hPa --> altitude [m] look-up table
#-----------------------------------------------------------------------------

altitudes = []
for j in range(len(levels)):
    h = hpa_to_altitude(levels[j])
    altitudes.append(h)                   

#-----------------------------------------------------------------------------
# CREATE: altitude [m] --> hPa look-up table
#-----------------------------------------------------------------------------

p_ben_nevis, t_ben_nevis, d_ben_nevis = isa(1345.0)
p_everest, t_everest, d_everest = isa(8848.86)

pressures = []
densities = []
temperatures = []
for j in range(len(altitudes)):
    p, t, d = isa(altitudes[j])
    pressures.append(p)                   
    densities.append(d)                   
    temperatures.append(t)                   

#-----------------------------------------------------------------------------
# CREATE: station altitude [m] --> hPa look-up table
#-----------------------------------------------------------------------------

station_pressures = []
station_densities = []
station_temperatures = []
for j in range(len(station_altitudes)):
    if station_altitudes[j] < 0:
        p, t, d = [np.nan, np.nan, np.nan]
    else:
        p, t, d = isa(station_altitudes[j])
    station_pressures.append(np.round(p,1))                   
    station_densities.append(d)                   
    station_temperatures.append(t)                   

#-----------------------------------------------------------------------------
# FIND: nearest 20CRv3 pressure level
#-----------------------------------------------------------------------------
    
station_20CRv3_level_idx = [ np.absolute(levels-station_pressures[i]).argmin() for i in range(len(lut)) ] 
station_20CRv3_level = [ int(levels[station_20CRv3_level_idx[i]]) for i in range(len(lut)) ] 

lut['hPa'] = station_pressures
lut['20CRv3 hPa'] = station_20CRv3_level
lut['20CRv3 hPa level'] = station_20CRv3_level_idx
lut['20CRv3 hPa'][lut.station_elevation==-999] = np.nan
lut['20CRv3 hPa level'][lut.station_elevation==-999] = np.nan

#-----------------------------------------------------------------------------
# FIT: gamma function to distribution of station hPa - 20CRv3 nearest level hPa differencea
#-----------------------------------------------------------------------------
    
x = lut['hPa']-lut['20CRv3 hPa']
bias = np.nanmean(x)       
mask = np.isfinite(x) 
y = x[mask]            
#z = (y-y.mean())/y.std()
z = y
n = len(z)
z_map = np.linspace(z.min(), z.max(), n)
z_maxsigma = np.ceil(np.max([np.abs(z_map.min()),np.abs(z_map.min())]))

disttype = 'gamma'
dist = getattr(scipy.stats, disttype)
param = dist.fit(z)
gamma_pdf = dist.pdf(z_map, param[0], param[1], param[2])	
gamma = dist.rvs(param[0], loc=param[1], scale=param[2], size=n)
gamma_sorted = np.sort(gamma)
gamma_bins = np.percentile(gamma_sorted, range(0,101))

#-----------------------------------------------------------------------------
# WRITE: look-up table to CSV
#-----------------------------------------------------------------------------

lut.to_csv( 'lut_hPa.csv')

#-----------------------------------------------------------------------------
# PLOT: altitude v pressure curves
#-----------------------------------------------------------------------------
    
img = image.imread(image_file)

figstr = '20crv3-pressure-curve-fit.png'
titlestr = '20CRv3 pressure altitude curve'

fig,ax = plt.subplots(figsize=(15,10))
plt.plot(levels, altitudes, 'o', color='purple', alpha=0.7, lw=3, ls='-', label='ISA conversion')
plt.plot(pressures, altitudes, 'o', color='cyan', alpha=0.7, lw=3,  ls='-', label='NIST conversion')
#altitude = 1345 # Ben Nevis
plt.axhline(y=1345.0, color='orange', ls='--', lw=1)
plt.axhline(y=8848.86, color='lime', ls='--', lw=1)
plt.axvline(x=p_ben_nevis, color='orange', ls='--', lw=1, label=r'Ben Nevis: 1345m, '+str(np.round(p_ben_nevis,1))+'hPa')
plt.axvline(x=p_everest, color='lime', ls='--', lw=1, label=r'Mt. Everest: 8848m, '+str(np.round(p_everest,1))+'hPa')
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
aspect = img.shape[0] / img.shape[1] * (xmax - xmin)/(ymax - ymin)
plt.imshow(img, zorder=0, extent=[xmin, xmax, ymin, ymax], aspect=aspect, alpha=0.8)
plt.tick_params(labelsize=fontsize)
plt.legend(loc='lower left', bbox_to_anchor=(0, -0.6), ncol=1, facecolor='lightgrey', framealpha=1, fontsize=fontsize)    
fig.subplots_adjust(left=None, bottom=0.4, right=None, top=None, wspace=None, hspace=None)             
plt.xlabel("Pressure, hPa", fontsize=fontsize)
plt.ylabel("Altitude, m", fontsize=fontsize)
plt.title(titlestr, fontsize=fontsize, fontweight='normal')
plt.savefig(figstr, dpi=300)
plt.close('all')

#-----------------------------------------------------------------------------
# PLOT: histogram of station hPa - 20CRv3 nearest level hPa differences
#-----------------------------------------------------------------------------
    
figstr = 'station-20crv3-pressure-histogram.png'
titlestr = 'Station - nearest 20CRv3 pressure level histogram'
                
fig,ax = plt.subplots(figsize=(15,10))
h1 = ax.hist(gamma, bins=100, color='cyan', alpha=0.2, density=True, label=r'$\Gamma$ fit (draws): n='+str(n))
#h2 = ax.hist(z, bins=100, color='grey', alpha=0.2, density=True, label='Difference distribution: n='+str(n))
h3 = plt.plot(z_map, gamma_pdf, color='teal', lw=5, label=r'$\Gamma$ fit:'+
         r' ($\alpha=$'+str(np.round(param[0],2))+','+
         r' loc='+str(np.round(param[1],2))+','+
         r' scale='+str(np.round(param[2],2))+')')
#ymin = np.min([h1[1].min(),h2[1].min()])
#ymax = np.max([h1[1].max(),h2[1].max()])
ymin = np.min(h1[1].min())
ymax = np.max(h1[1].max())
plt.axvline(gamma_bins[50], ymin, ymax, color='black', ls='--', label=r'$\Gamma$ fit: median bias='+str(np.round(gamma_bins[50],3))+'hPa')
plt.legend(loc='lower left', bbox_to_anchor=(0, -0.4), markerscale=1, ncol=1, facecolor='lightgrey', framealpha=1, fontsize=fontsize)    
fig.subplots_adjust(left=None, bottom=0.3, right=None, top=None, wspace=None, hspace=None)     
plt.xlabel(r'Station-20CRv3 pressure difference, $hPa$', fontsize=fontsize)
plt.ylabel(r'Kernel density estimate (KDE)', fontsize=fontsize)
ax.set_xlim(-z_maxsigma,z_maxsigma)
plt.tick_params(axis='both', which='major', labelsize=fontsize)
plt.title(titlestr, fontsize=fontsize)
plt.savefig(figstr)
plt.close('all')       

#-----------------------------------------------------------------------------
print('** END')





