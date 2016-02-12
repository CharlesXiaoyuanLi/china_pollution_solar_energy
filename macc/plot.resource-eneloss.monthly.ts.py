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

# Read base radiation data
nc_file = '../macc/prov.monthly.climatology.powerden.nc'
fh = Dataset(nc_file, 'r')

prov_pwden_monthly = fh.variables['prov_all_sky_pwden_hzl'][:][:]

prov_enebase_monthly = prov_pwden_monthly / prov_cosz  # units: [kWh/m^2/day]

fh.close()

# Pick Beijing for plot
bj_eneloss_monthly = prov_eneloss_monthly[-1,:]
bj_enebase_monthly = prov_enebase_monthly[-1,:]

bj_monthly = np.zeros((12,2))

bj_monthly[:,0] = bj_enebase_monthly
bj_monthly[:,1] = bj_eneloss_monthly

print bj_monthly

### Plot the timeseries
import matplotlib.pyplot as plt
import pandas as pd

indexstr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

ts = pd.DataFrame(bj_monthly,index=indexstr,columns=['Base Resource Potential','Loss due to Air Pollution'])

ts.plot(secondary_y='Loss due to Air Pollution',color=['k','tomato'],lw=2.,fontsize=12)

plt.gca().invert_yaxis()

plt.title('Beijing',fontsize=25)
plt.title('$kW \cdot h / m^2$',loc='right',fontsize=18)
plt.tight_layout()
plt.savefig('./Figures/Beijing_resource-eneloss_monthly.ps',format='ps')
plt.show()
