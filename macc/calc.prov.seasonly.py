### This program calculates seasonly mean data
### and save them to ncfile.

from netCDF4 import Dataset
import numpy as np
import scipy.ndimage
import types
import sys

nprov = 31

# Create 12-year daily mean time series

ttime = 4*12

prov_pwden = np.zeros((ttime,nprov),dtype=float)

dim_time = np.zeros(ttime)

pos_time = 0


# Get daily data from each year file
for year in np.arange(2003,2015,1):

### Read data files

# Open netCDF file
    nc_file = '../macc/macc_aero_'+str(year)+'_seasonly.nc'
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

    pro_area = np.zeros(nprov)
    pro_power = np.zeros((ntime,nprov))

    drop = np.zeros(nprov)

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
    pro_powerden = np.zeros((ntime,nprov))
    for i in np.arange(nprov):  #loop through provinces
        if pro_area[i] == 0:
            pro_powerden[:,i] = 0
        else:
            pro_powerden[:,i] = pro_power[:,i] / pro_area[i] *1000      # [W/m^2]

# write to large array
    prov_pwden[pos_time:pos_time+4,:] = pro_powerden
    dim_time[pos_time:pos_time+4] = np.arange(0,4) + pos_time
    pos_time += 4

### Write the data to ncfile
outfile = Dataset('prov.seasonly.powerden.nc','w',format='NETCDF4')

outfile.createDimension('province',nprov)
outfile.createDimension('t',ttime)

prov_dim = outfile.createVariable('province',np.int32,('province',))
time_dim   = outfile.createVariable('t',np.int32,('t',))


description = ''
rgid = 0
for prov in pro_nm:
    if rgid == 31:
        break
    description += prov + ': ' + str(rgid) + ',  '
    rgid += 1
prov_dim.description =description

prov_pwden_out = outfile.createVariable('prov_pwden',np.float32,('province','t'))
prov_area_out   = outfile.createVariable('prov_area',np.float32,('province',))

prov_dim[:] = np.arange(0,nprov)
time_dim[:]   = dim_time
time_dim.units = 'seasons since DJF 2003'
time_dim.description = 'order of seasons: DJF, MAM, JJA, SON'

prov_pwden_out[:] = np.transpose(prov_pwden)
prov_pwden_out.units = 'W/m^2'
prov_area_out[:] = pro_area

outfile.close()
