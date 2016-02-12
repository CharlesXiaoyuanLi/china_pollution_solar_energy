## This program calculates and plots the daily mean time series

from netCDF4 import Dataset
import numpy as np
import scipy.ndimage
import types


for year in np.arange(2003,2015,1):

### Read data files

# Open netCDF file
    nc_file = '../macc/macc_aero_anthsfc_'+str(year)+'_daily.nc'
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

# Calculate density
    pro_eneden = np.zeros((ntime,31))
    pro_powerden = np.zeros((ntime,31))
    for i in np.arange(31):  #loop through provinces
        pro_powerden[:,i] = pro_power[:,i] / pro_area[i] *1000      # [W/m^2]
        pro_eneden[:,i] = pro_power[:,i] * 12 / pro_area[i]    # [kW h / m^2]

    print pro_eneden

# Plot

    units = '($W / m^2$)'
    varnm = 'pwden'
    varlim = (-120,0) #pwden
#    varlim = (-1.4,0) #eneden

    import matplotlib.pyplot as plt
    from pandas import *

# ts = DataFrame(pro_eneden[:,:], index=date_range('1/2003',periods=144))
    for pro_n in np.arange(31):

        plt.figure()

        if year%4 == 0:
            peri = 366
        else:
            peri = 365

        ts = Series(pro_powerden[:,pro_n],index=date_range('1/1/'+str(year),periods=peri))

        ts.plot()
        
        plt.ylim(varlim)
        plt.title(pro_nm[pro_n]+'\t\tDaily-mean\t\t'+units)
        plt.savefig('./Figures/Power Density/daily/'+pro_nm[pro_n]+'/'+pro_nm[pro_n]+'_'+str(year)+'_daily_'+varnm+'_ts',format='ps')
