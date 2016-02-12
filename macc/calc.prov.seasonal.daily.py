### This program reorganizes daily mean data into four seasons,
### and save them to ncfile.

from netCDF4 import Dataset
import numpy as np
import scipy.ndimage
import types
import sys

nprov = 31

# Create 12-year daily mean time series for each season

ndjf = 12*31+12*31+3*29+9*28
nmam = 12*31+12*30+12*31
njja = nmam
nson = 12*30+12*31+12*30

prov_pwden_djf = np.zeros((ndjf,nprov),dtype=float)
prov_pwden_mam = np.zeros((nmam,nprov),dtype=float)
prov_pwden_jja = np.zeros((njja,nprov),dtype=float)
prov_pwden_son = np.zeros((nson,nprov),dtype=float)

dim_djf = np.zeros(ndjf)
dim_mam = np.zeros(nmam)
dim_jja = np.zeros(njja)
dim_son = np.zeros(nson)

djf_array_comm = np.append(np.append(1200+np.arange(1,32),100+np.arange(1,32)),200+np.arange(1,29))
djf_array_leap = np.append(np.append(1200+np.arange(1,32),100+np.arange(1,32)),200+np.arange(1,30))
mam_array = np.append(np.append(300+np.arange(1,32),400+np.arange(1,31)),500+np.arange(1,32))
jja_array = np.append(np.append(600+np.arange(1,31),700+np.arange(1,32)),800+np.arange(1,32))
son_array = np.append(np.append(900+np.arange(1,31),1000+np.arange(1,32)),1100+np.arange(1,31))

pos_djf = 0
pos_mam = 0
pos_jja = 0
pos_son = 0

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
    pro_eneden = np.zeros((ntime,nprov))
    pro_powerden = np.zeros((ntime,nprov))
    for i in np.arange(nprov):  #loop through provinces
        pro_powerden[:,i] = pro_power[:,i] / pro_area[i] *1000      # [W/m^2]
        pro_eneden[:,i] = pro_power[:,i] * 12 / pro_area[i]    # [kW h / m^2]

    print pro_eneden

# Plot

    units = '($kW / m^2$)'
    varnm = 'pwden'
    varlim = (-120,0) #pwden
#    varlim = (-1.4,0) #eneden
    
    pos_prov = 0
# write to large array
    # DJF
    prov_pwden_djf[pos_djf:pos_djf+31,:] = pro_powerden[-31:,:]
    pos_djf += 31
    if year%4 == 0: 
        prov_pwden_djf[pos_djf:pos_djf+31+29,:] = pro_powerden[pos_prov:31+29,:]
        dim_djf[pos_djf-31:pos_djf+31+29] = year*10000+djf_array_leap
        pos_djf += 31+29
        pos_prov += 31+29
    else:
        prov_pwden_djf[pos_djf:pos_djf+31+28,:] = pro_powerden[pos_prov:31+28,:]
        dim_djf[pos_djf-31:pos_djf+31+28] = year*10000+djf_array_comm
        pos_djf += 31+28
        pos_prov += 31+28
    # MAM
    prov_pwden_mam[pos_mam:pos_mam+92,:] = pro_powerden[pos_prov:pos_prov+92,:]
    dim_mam[pos_mam:pos_mam+92] = year*10000+mam_array
    pos_mam += 92
    pos_prov += 92
    # JJA
    prov_pwden_jja[pos_jja:pos_jja+92,:] = pro_powerden[pos_prov:pos_prov+92,:]
    dim_jja[pos_jja:pos_jja+92] = year*10000+jja_array
    pos_jja += 92
    pos_prov += 92
    # SON
    prov_pwden_son[pos_son:pos_son+91,:] = pro_powerden[pos_prov:pos_prov+91,:]
    dim_son[pos_son:pos_son+91] = year*10000+son_array
    pos_son += 91
    pos_prov += 91

    print pos_prov+31

### Write the data to ncfile
outfile = Dataset('prov.seasonal.daily.powerden.nc','w',format='NETCDF4')

outfile.createDimension('province',nprov)
outfile.createDimension('time_djf',ndjf)
outfile.createDimension('time_mam',nmam)
outfile.createDimension('time_jja',njja)
outfile.createDimension('time_son',nson)

prov_dim = outfile.createVariable('province',np.int32,('province',))
djf_dim   = outfile.createVariable('time_djf',np.int32,('time_djf',))
mam_dim   = outfile.createVariable('time_mam',np.int32,('time_mam',))
jja_dim   = outfile.createVariable('time_jja',np.int32,('time_jja',))
son_dim   = outfile.createVariable('time_son',np.int32,('time_son',))


description = ''
rgid = 0
for prov in pro_nm:
    if rgid == 31:
        break
    description += prov + ': ' + str(rgid) + ',  '
    rgid += 1
prov_dim.description =description

prov_pwden_djf_out = outfile.createVariable('prov_pwden_djf',np.float32,('province','time_djf'))
prov_pwden_mam_out = outfile.createVariable('prov_pwden_mam',np.float32,('province','time_mam'))
prov_pwden_jja_out = outfile.createVariable('prov_pwden_jja',np.float32,('province','time_jja'))
prov_pwden_son_out = outfile.createVariable('prov_pwden_son',np.float32,('province','time_son'))

prov_dim[:] = np.arange(0,nprov)
djf_dim[:]   = dim_djf
mam_dim[:]   = dim_mam
jja_dim[:]   = dim_jja
son_dim[:]   = dim_son

prov_pwden_djf_out[:] = np.transpose(prov_pwden_djf)
prov_pwden_djf_out.units = 'W/m^2'
prov_pwden_mam_out[:] = np.transpose(prov_pwden_mam)
prov_pwden_mam_out.units = 'W/m^2'
prov_pwden_jja_out[:] = np.transpose(prov_pwden_jja)
prov_pwden_jja_out.units = 'W/m^2'
prov_pwden_son_out[:] = np.transpose(prov_pwden_son)
prov_pwden_son_out.units = 'W/m^2'

outfile.close()
