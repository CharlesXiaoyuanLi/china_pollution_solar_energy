## This program calculates and plots the monthly mean time series 

import numpy as np
import scipy.ndimage
import types

### Create shortwave and cosz data Matrix
ntime = 12
nlat = 180
nlon = 360

time = np.arange(1,13,1)
latst = np.arange(89,-91,-1)
lonst = np.arange(0,360,1)

swvt = np.zeros((ntime, nlat, nlon))
coszt = np.zeros((ntime, nlat, nlon))
clrt = np.zeros((ntime, nlat, nlon))

### Read data files

# Open shortwave file
swv_file = './22yr_swv_dwn_last_ann'
fh = open(swv_file, 'r')

fstring = fh.read().split()
istr=0
for s in np.arange(len(fstring)/15):
    flat = fstring[istr]
    flon = fstring[istr+1]
    i = 89 - int(flat)
    if flon >= 0:
        j = int(flon)
    else:
        j = int(flon) + 360
    for tm in np.arange(12):
        swvt[tm,i,j] = float(fstring[istr+2+tm])
    istr += 15

fh.close()
del(fstring)

# Open clearsky file
clr_file = './22yr_clr_sky_last_ann'
fh = open(clr_file, 'r')

fstring = fh.read().split()
istr=0
for s in np.arange(len(fstring)/15):
    flat = fstring[istr]
    flon = fstring[istr+1]
    i = 89 - int(flat)
    if flon >= 0:
        j = int(flon)
    else:
        j = int(flon) + 360
    for tm in np.arange(12):
        clrt[tm,i,j] = float(fstring[istr+2+tm])
    istr += 15

fh.close()
del(fstring)


# Open cosine zenith angle file
zmt_file = './nasa_sse_6_cos_zmt.txt'
fh = open(zmt_file, 'r')

fstring = fh.read().split()
istr=0
for s in np.arange(len(fstring)/14):
    flat = fstring[istr]
    flon = fstring[istr+1]
    i = 89 - int(flat)
    if flon >= 0:
        j = int(flon)
    else:
        j = int(flon) + 360
    for tm in np.arange(12):
        if fstring[istr+2+tm] != 'na':
            coszt[tm,i,j] = float(fstring[istr+2+tm])
        else:
            coszt[tm,i,j] = 1
    istr += 14

fh.close()
del(fstring)

# Extract data from the following domain
lonllim = 70
lonrlim = 140
latblim = 10
latulim = 54

# Calculate the index boundaries of the domain
lonlind = np.where((lonst == lonllim))[0]
lonrind = np.where((lonst == lonrlim))[0]
latbind = np.where((latst == latblim))[0]
latuind = np.where((latst == latulim))[0]

lons = lonst[lonlind:lonrind]
lats = latst[latuind:latbind]
swv = swvt[:,latuind:latbind,lonlind:lonrind]   # units:kwh/m2/day
cosz = coszt[:,latuind:latbind,lonlind:lonrind]
clr = clrt[:,latuind:latbind,lonlind:lonrind]   # units:kwh/m2/day

print 'Domain: longitude('+str(lonllim)+'~'+str(lonrlim)+'),'+'latitude('+str(latblim)+'~'+str(latulim)+')'
print 'nlon = '+str(lons.size)
print 'nlat = '+str(lats.size)

### Calculate total power at each grid

# Calcualte the area weight array
R = 6378100  # radius of the earth in [m] http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html 
# aw = length in lat  *  length in lon
aw = 2. * np.pi * R * np.cos(np.pi / 180. * lats) / nlon    *    np.pi / 180. * (lats[0] - lats[1]) * R

# Convert kWh m-2 day-1 to kWh day-1 at each grid point
swv_hzl = np.zeros( swv.shape )
swv_nml = np.zeros( swv.shape )
clr_hzl = np.zeros( clr.shape )
for i in np.arange(ntime):
    temp = aw * swv[i,:,:].T
    swv_hzl[i,:,:] = temp.T
    swv_nml[i,:,:] = temp.T / cosz[i,:,:]
    temp_clr = aw * clr[i,:,:].T
    clr_hzl[i,:,:] = temp_clr.T

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

### Calculate total area and energy of each province

pro_area = np.zeros(31)
pro_swv_hzl = np.zeros((ntime,31))
pro_swv_nml = np.zeros((ntime,31))
pro_clr_hzl = np.zeros((ntime,31))

drop = np.zeros(31)

for i in np.arange(cn_pro_map.shape[0]):
    for j in np.arange(cn_pro_map.shape[1]):
        
        if cn_pro_map[i,j] != 0:
            aind = cn_pro_map[i,j]

            pro_area[aind-1] += aw[i]  # total area
            pro_swv_hzl[:,aind-1] += swv_hzl[:,i,j]  # total horizontal energy
            pro_swv_nml[:,aind-1] += swv_nml[:,i,j]  # total normal energy
            pro_clr_hzl[:,aind-1] += clr_hzl[:,i,j]  # total clear-sky energy

            drop[aind-1] += 1

print drop


### Calculate power/energy density for each province

# Create array for days in month
dinmon = np.array([31,28,31,30,31,30,31,31,30,31,30,31]) 

# Calculate density
pro_hzl_eneden = np.zeros((ntime,31))
pro_nml_eneden = np.zeros((ntime,31))
pro_clr_hzl_eneden = np.zeros((ntime,31))
for i in np.arange(31):  #loop through provinces
    pro_hzl_eneden[:,i] = pro_swv_hzl[:,i] * dinmon / pro_area[i]    # [kW h / m^2]
    pro_nml_eneden[:,i] = pro_swv_nml[:,i] * dinmon / pro_area[i]    # [kW h / m^2]
    pro_clr_hzl_eneden[:,i] = pro_clr_hzl[:,i] * dinmon / pro_area[i]


# Provincial monthly-mean cloud fraction
pro_cloud_frac = pro_hzl_eneden / pro_clr_hzl_eneden

f = open('pro_cloud_frac.txt','w+')
for i in range(31):
    for j in range(ntime):
        f.write(str(pro_cloud_frac[j][i]) + ' ')
    f.write('\n')
f.close()

# Provincial monthly-mean cosine zenith angle
pro_cosz = pro_hzl_eneden / pro_nml_eneden

f = open('pro_cosz.txt','w+')
for i in range(31):
    for j in range(ntime):
        f.write(str(pro_cosz[j][i]) + ' ')
    f.write('\n')
f.close()


# Provincial annual total normal energy Density
pro_nml_ann_eneden = np.sum(pro_nml_eneden, axis=0 )

print pro_nml_ann_eneden
#print np.sum(pro_hzl_eneden, axis=0)
f = open('pro_nml_ann_eneden.txt','w+')
for i in range(31):
    f.write(str(pro_nml_ann_eneden[i]) + '\n')
f.close()

# Print Beijing powerden

f = open('beijing_powerden.txt','w+')
for i in range(12):
    f.write( str(pro_swv_hzl[i,30]/pro_area[30]) + '\n' )
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

    ts = Series(pro_powerden[:,pro_n], index=date_range('1/1/1983',periods=144,freq='M'))

    ts.plot()
    
    plt.ylim(ylim)
    plt.title(pro_nm[pro_n]+'\t\tMonthly-mean\t\t'+units)
    plt.savefig('./Figures/Power Density/monthly/'+pro_nm[pro_n]+'_'+varnm+'_monthly_ts',format='ps')
