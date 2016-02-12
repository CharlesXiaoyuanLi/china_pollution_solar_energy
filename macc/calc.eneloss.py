## This program plots the seasonal/annual mean time series

from netCDF4 import Dataset
import numpy as np
import sys

### Read data files

# Open netCDF file
nc_file = '../macc/prov.monthly.powerden.nc'
fh = Dataset(nc_file, 'r')

time = fh.variables['t'][:]
pid  = fh.variables['province'][:]
prov_pwden_monthly_t = fh.variables['prov_pwden'][:][:]
prov_area  = fh.variables['prov_area'][:]

fh.close()

ntime = time.size
nprov = pid.size

prov_pwden_monthly = np.transpose(np.array([np.average(prov_pwden_monthly_t[:,i::12],axis=1) for i in
    np.arange(12)]))

prov_cldfrac = np.zeros((nprov,12))
prov_cosz    = np.zeros((nprov,12))

# Read cloud fraction and cosz into the 2D array from the file
f = open("../pro_cloud_frac.txt","r")
fstring = f.read().split()
istr = 0
for i in range(nprov):
    for j in range(12):
        prov_cldfrac[i,j] = float(fstring[istr])
        istr += 1
f.close()

f = open("../pro_cosz.txt","r")
fstring = f.read().split()
istr = 0
for i in range(nprov):
    for j in range(12):
        prov_cosz[i,j] = float(fstring[istr])
        istr += 1
f.close()

prov_eneloss_monthly = prov_pwden_monthly * prov_cldfrac / prov_cosz * 24. / 1000 # units: [kWh/m^2/day]

prov_eneloss_annual = np.average(prov_eneloss_monthly,axis=1)

f = open("./prov_unit_eneloss_annualavg_daily_all_sky.txt","w+")
for i in range(nprov):
    if i != 0:
        f.write('\n')
    f.write(str(prov_eneloss_annual[i]))
f.close()

f = open("./prov_unit_eneloss_annual_all_sky.txt","w+")
for i in range(nprov):
    if i != 0:
        f.write('\n')
    f.write(str(prov_eneloss_annual[i]*365))
f.close()


# Read base radiation data
nc_file = '../macc/prov.monthly.climatology.powerden.nc'
fh = Dataset(nc_file, 'r')

prov_pwden_monthly = fh.variables['prov_all_sky_pwden_hzl'][:][:]

prov_enebase_monthly = prov_pwden_monthly / prov_cosz  # units: [kWh/m^2/day]

prov_enebase_annual  = np.average(prov_enebase_monthly,axis=1)

fh.close()

f = open("./prov_unit_enebase_annualavg_daily_all_sky.txt","w+")
for i in range(nprov):
    if i != 0:
        f.write('\n')
    f.write(str(prov_enebase_annual[i]))
f.close()

f = open("./prov_unit_enebase_annual_all_sky.txt","w+")
for i in range(nprov):
    if i != 0:
        f.write('\n')
    f.write(str(prov_enebase_annual[i]*365))
f.close()


