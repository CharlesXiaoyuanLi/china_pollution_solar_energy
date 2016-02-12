## This program calculates and plots the monthly mean time series 

from netCDF4 import Dataset
import numpy as np
import scipy.ndimage
import types

### Read data files

# Open netCDF file
nc_file = '../macc/macc_aero_anthsfc_2003_2014_monthly.nc'
fh = Dataset(nc_file, 'r')

time = fh.variables['t'][:]
lonst = fh.variables['longitude'][:]
latst = fh.variables['latitude'][:]
anthsfct = fh.variables['anthsrf'][:][:][:]

ntime = time.size

nlat = latst.size
nlon = lonst.size

# Extract data from the following domain
lonllim = 70
lonrlim = 140
latblim = 10
latulim = 54

# Calculate the index boundaries of the domain
lonlind = np.where((lonst >= lonllim) & (lonst < lonllim+1))[0]
lonrind = np.where((lonst >lonrlim-1) & (lonst <= lonrlim))[0]
latbind = np.where((latst >= latblim) & (latst < latblim+1))[0]
latuind = np.where((latst > latulim-1) & (latst <= latulim))[0]

lons = lonst[lonlind:lonrind]
lats = latst[latuind:latbind]
anthsfc = anthsfct[:,latuind:latbind,lonlind:lonrind]

print 'Domain: longitude('+str(lonllim)+'~'+str(lonrlim)+'),'+'latitude('+str(latblim)+'~'+str(latulim)+')'
print 'nlon = '+str(lons.size)
print 'nlat = '+str(lats.size)

units = fh.variables['anthsrf'].units

fh.close()

### Calculate total power at each grid

# Calcualte the area weight array
R = 6378100  # radius of the earth in [m] http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html 
# aw = length in lat  *  length in lon
aw = 2. * np.pi * R * np.cos(np.pi / 180. * lats) / nlon    *    np.pi / 180. * (lats[0] - lats[1]) * R

# Convert W m-2 to kW at each grid point
anthsfc_power = np.zeros( anthsfc.shape )
for i in np.arange(ntime):
    temp = aw * anthsfc[i,:,:].T / 1000.
    anthsfc_power[i,:,:] = temp.T


### Create China Provincial Matrix

# Create a 2D array representing lat x lon
cn_pro_map_t = np.zeros( (nlat,nlon) )

# Read value into the 2D array from the file
f = open("../cn_province_"+str(nlat)+"x"+str(nlon)+".txt","r")
fstring = f.read().split()
istr = 0
for i in range(nlat):
    for j in range(nlon):
        cn_pro_map_t[i,j] = int(fstring[istr])
        istr += 1

# Only take the part for the required domain
cn_pro_map = cn_pro_map_t[latuind:latbind,lonlind:lonrind]

f.close()


### Read in China Province Namelist

f = open("../cn_province_nmlist.txt","r")
pro_nm = f.read().split("\n")
f.close()

print pro_nm

### Calculate total area and power of each province

pro_area = np.zeros(31)
pro_power = np.zeros((ntime,31))

drop = np.zeros(31)

for i in np.arange(cn_pro_map.shape[0]):
    for j in np.arange(cn_pro_map.shape[1]):
        
        if cn_pro_map[i,j] != 0:
            aind = cn_pro_map[i,j]

            pro_area[aind-1] += aw[i]  # total area
            pro_power[:,aind-1] += anthsfc_power[:,i,j]  # total power

            drop[aind-1] += 1

print drop
print pro_area
print pro_power.shape


### Calculate power/energy density for each province

# Create array for days in month
comm_yr = np.array([31,28,31,30,31,30,31,31,30,31,30,31]) 
leap_yr = np.array([31,29,31,30,31,30,31,31,30,31,30,31])

dinmon = comm_yr # 2003

for yr in np.arange(2004,2015,1): #2004-2014
    if yr%4 == 0:
        dinmon = np.append(dinmon,leap_yr)
    else:
        dinmon = np.append(dinmon,comm_yr)
# calculate hours in month
hrinmon = dinmon * 12.

# Calculate density
pro_eneden = np.zeros((ntime,31))
pro_powerden = np.zeros((ntime,31))
for i in np.arange(31):  #loop through provinces
    pro_powerden[:,i] = pro_power[:,i] / pro_area[i] *1000           # [W/m^2]
    pro_eneden[:,i] = pro_power[:,i] * hrinmon / pro_area[i]    # [kW h / m^2]

print pro_eneden

# Calculate 12-year (2003-2014) averaged monthly-mean provincial energy density loss
pro_eneden_avg = np.zeros( (12,31) )

for i in range(12):
    pro_eneden_avg[i,:] = np.mean(pro_eneden[i:pro_eneden.shape[0]:12,:],axis=0)


# Read Provincial monthly-mean cloud fraction
pro_cloud_frac = np.zeros( (12,31) )
print pro_cloud_frac.shape
f = open('../pro_cloud_frac.txt','r')
fstring = f.read().split()
istr=0
for i in range(31):
    for j in range(12):
        pro_cloud_frac[j][i] = float(fstring[istr])
        istr += 1

f.close()
del(fstring)

# Read Provincial monthly-mean cosine zenith angle
pro_cosz = np.zeros( (12,31) )

f = open('../pro_cosz.txt','r')
fstring = f.read().split()
istr=0
for i in range(31):
    for j in range(12):
        pro_cosz[j][i] = float(fstring[istr])
        istr += 1

f.close()
del(fstring)


# Calculate all-sky normal 12-year (2003-2014) averaged monthly-mean provincial
# energy density loss due to aerosol
pro_nml_allsky_eneden_loss_month = pro_eneden_avg / pro_cosz * pro_cloud_frac

# calcualte averaged annual total
pro_nml_allsky_eneden_loss_ann = np.sum(pro_nml_allsky_eneden_loss_month,axis=0) 

# Print to file
f = open('pro_nml_allsky_eneden_loss_ann.txt','w+')
for i in range(31):
    f.write( str(pro_nml_allsky_eneden_loss_ann[i]) + '\n')

f.close()

# Print Beijing Powerden

f = open('beijing_powerden_loss.txt','w+')
for i in range(ntime):
    f.write( str(pro_powerden[i,30]) + '\n' )

f.close()




exit(0)
# Plot

#units = '($kW \cdot h / m^2$)'
units = '($W / m^2$)'
#ylim = (-25,0)  # eneden
ylim = (-60,0)  # pwden
#varnm = 'eneden'
varnm = 'pwden'

import matplotlib.pyplot as plt
from pandas import *

# ts = DataFrame(pro_eneden[:,:], index=date_range('1/2003',periods=144))
for pro_n in np.arange(31):

    plt.figure()

    ts = Series(pro_powerden[:,pro_n], index=date_range('1/1/2003',periods=144,freq='M'))

    ts.plot()
    
    plt.ylim(ylim)
    plt.title(pro_nm[pro_n]+'\t\tMonthly-mean\t\t'+units)
    plt.savefig('./Figures/Power Density/monthly/'+pro_nm[pro_n]+'_'+varnm+'_monthly_ts',format='ps')
