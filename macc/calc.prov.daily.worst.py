### This program plots daily mean data by season.

from netCDF4 import Dataset
import numpy as np

# Read in data
ncfile = './prov.seasonal.daily.eneloss.nc'
fh = Dataset(ncfile, 'r')

prov_eneloss_djf = fh.variables['prov_eneloss_djf'][:][:]
prov_eneloss_mam = fh.variables['prov_eneloss_mam'][:][:]
prov_eneloss_jja = fh.variables['prov_eneloss_jja'][:][:]
prov_eneloss_son = fh.variables['prov_eneloss_son'][:][:]

prov = fh.variables['province'][:]

# Read in province name
f = open("../cn_province_nmlist.txt","r")
pro_nm_t = f.read().split("\n")
f.close()
pro_nm = pro_nm_t[:-1] # dump the last empty element

prov_eneloss_all = np.concatenate((prov_eneloss_djf,prov_eneloss_mam,prov_eneloss_jja,prov_eneloss_son),axis=1)

prov_worst = np.array( [ prov_eneloss_all[i,:].min() for i in range(31) ] )
print prov_worst

f = open("../plot/prov_worst_daily_eneloss.txt","w+")
for i in range(31):
    if i != 0:
        f.write('\n')
    f.write(str(prov_worst[i]))
    
f.close()

ind_worst = np.where( prov_eneloss_all == prov_eneloss_all.min() )

t_ind = ind_worst[1][0]

prov_worst_day = prov_eneloss_all[:,t_ind]

f = open("../plot/prov_worst_day_eneloss.txt","w+")
for i in range(31):
    if i != 0:
        f.write('\n')
    f.write(str(prov_worst_day[i]))
    
f.close()

print t_ind

