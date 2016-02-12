## This program plots the seasonal/annual mean time series

from netCDF4 import Dataset
import numpy as np
import sys

### Read data files

# Open netCDF file
nc_file = '../macc/prov.seasonly.powerden.nc'
fh = Dataset(nc_file, 'r')

time = fh.variables['t'][:]
pid  = fh.variables['province'][:]
prov_pwden_seasonal = fh.variables['prov_pwden'][:][:]
prov_area  = fh.variables['prov_area'][:]

ntime = time.size
nprov = pid.size

prov_pwden_annual = np.transpose(np.array([np.average(prov_pwden_seasonal[:,i:i+4],axis=1) for i in
    np.arange(0,ntime,4)]))

nyr = prov_pwden_annual.shape[1]

prov_pwden = np.zeros((nprov,nyr,5))

prov_pwden[:,:,0] = prov_pwden_annual

for i in np.arange(4):
    prov_pwden[:,:,i+1] = prov_pwden_seasonal[:,i::4]

cn_pwden = np.average(prov_pwden,axis=0,weights=prov_area)

print cn_pwden[:,0]

### Read in China Province Namelist

f = open("../cn_province_nmlist.txt","r")
pro_nm = f.read().split("\n")
f.close()

print pro_nm

### Plot

units = '$W / m^2$'
varnm = 'pwden'
varlim = (-120,0) #pwden
#varlim = (-1.4,0) #eneden

import matplotlib.pyplot as plt
import pandas as pd

# ts = DataFrame(pro_eneden[:,:], index=date_range('1/2003',periods=144))
plt.figure()

ts = pd.DataFrame(cn_pwden,index=[str(i) for i in np.arange(2003,2015)],columns=['Annual','DJF','MAM','JJA','SON'])

ts.plot()

plt.ylabel(units)

#plt.ylim(varlim)
plt.title('China')

plt.savefig('./Figures/annual_ts/China_annual_seasonal_'+varnm+'_ts',format='ps')

plt.close()

for pro_n in np.arange(nprov):

    plt.figure()

    ts = pd.DataFrame(prov_pwden[pro_n,:,:],index=[str(i) for i in np.arange(2003,2015)],columns=['Annual','DJF','MAM','JJA','SON'])

    ts.plot()
   
    plt.ylabel(units)
#    plt.ylim(varlim)
    plt.title(pro_nm[pro_n]+'\t\t'+units)
    plt.savefig('./Figures/annual_ts/'+pro_nm[pro_n]+'_annual_seasonal_'+varnm+'_ts',format='ps')
    plt.close()
