### This program plots daily mean data by season.

from netCDF4 import Dataset
import numpy as np

# Read in data
ncfile = './prov.seasonal.daily.powerden.ratio.nc'
fh = Dataset(ncfile, 'r')

prov_pwden_djf = fh.variables['prov_pwden_ratio_djf'][:][:]
prov_pwden_mam = fh.variables['prov_pwden_ratio_mam'][:][:]
prov_pwden_jja = fh.variables['prov_pwden_ratio_jja'][:][:]
prov_pwden_son = fh.variables['prov_pwden_ratio_son'][:][:]

prov = fh.variables['province'][:]

# Read in province name
f = open("../cn_province_nmlist.txt","r")
pro_nm_t = f.read().split("\n")
f.close()
pro_nm = pro_nm_t[:-1] # dump the last empty element

### Plot frequency in histogram
import matplotlib.pyplot as plt

bins_set = np.arange(-50,2,2)
xticks = np.arange(-1,-51,-2)
xticklim = np.arange(0,-52,-2)
xticklabels = []

for i in np.arange(len(xticks)):
    xticklabels.append( str(xticklim[i])+' ~ '+str(xticklim[i+1]) )
print xticklabels

yinterval = np.arange(0,1.05,0.05)

normed_value = 2

for province,iprov in zip(pro_nm,prov):
    print province
    if province == 'Shanghai': continue
    plt.figure()
    plt.suptitle(province,fontsize=20)
    plt.subplot(221)
    hist, bins = np.histogram(prov_pwden_djf[iprov,:],bins=bins_set,density=True)
    hist *= normed_value
    maxnum = np.amax(hist)
    print bins
    widths = np.diff(bins)*0.8
    plt.bar(bins[:-1]+5*0.2*0.5,hist,color='tomato',width=widths)
#    plt.hist(prov_pwden_djf[iprov,:], bins=bins, histtype='bar',
#            color='tomato',rwidth=0.8,normed=True)
    plt.xlim(0,-50)
    plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
    plt.yticks(fontsize=7)
    plt.tick_params(axis='x',bottom='off')
    plt.title('DJF',loc='left')

    plt.subplot(222)
    hist, bins = np.histogram(prov_pwden_mam[iprov,:],bins=bins_set,density=True)
    hist *= normed_value
    if maxnum < np.amax(hist):
        maxnum = np.amax(hist)
    widths = np.diff(bins)*0.8
    plt.bar(bins[:-1]+5*0.2*0.5,hist,color='tomato',width=widths)
#    plt.hist(prov_pwden_mam[iprov,:], bins=bins, histtype='bar',
#            color='tomato',rwidth=0.8)
    plt.xlim(0,-50)
    plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
    plt.yticks(fontsize=7)
    plt.tick_params(axis='x',bottom='off')
    plt.title('MAM',loc='left')

    plt.subplot(223)
    hist, bins = np.histogram(prov_pwden_jja[iprov,:],bins=bins_set,density=True)
    hist *= normed_value
    if maxnum < np.amax(hist):
        maxnum = np.amax(hist)
    widths = np.diff(bins)*0.8
    plt.bar(bins[:-1]+5*0.2*0.5,hist,color='tomato',width=widths)
#    plt.hist(prov_pwden_jja[iprov,:], bins=bins, histtype='bar',
#            color='tomato',rwidth=0.8)
    plt.xlim(0,-50)
    plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
    plt.yticks(fontsize=7)

    plt.tick_params(axis='x',bottom='off')
    plt.title('JJA',loc='left')

    plt.subplot(224)
    hist, bins = np.histogram(prov_pwden_son[iprov,:],bins=bins_set,density=True)
    hist *= normed_value
    if maxnum < np.amax(hist):
        maxnum = np.amax(hist)
    widths = np.diff(bins)*0.8
    plt.bar(bins[:-1]+5*0.2*0.5,hist,color='tomato',width=widths)
#    plt.hist(prov_pwden_son[iprov,:], bins=bins, histtype='bar',
#            color='tomato',rwidth=0.8)
    plt.xlim(0,-50)
    plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
    plt.yticks(fontsize=7)

    plt.tick_params(axis='x',bottom='off')
    plt.title('SON',loc='left')

    # Create uniform limits for y-axis
    print maxnum
    maxnum_ind = np.where( yinterval >= maxnum )[0][0]
    ylimu = yinterval[maxnum_ind]
    plt.subplot(221)
    plt.ylim(0,ylimu)
    plt.subplot(222)
    plt.ylim(0,ylimu)
    plt.subplot(223)
    plt.ylim(0,ylimu)
    plt.subplot(224)
    plt.ylim(0,ylimu)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85) 
    plt.savefig('./Figures/season_freq/1Davg_ratio/'+province+'_1Davg_ratio_Freq_4seasons_aero_pwden',format='ps')
    plt.close()

    del hist,bins

