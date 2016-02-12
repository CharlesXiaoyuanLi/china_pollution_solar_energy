### This program plots daily mean data by season.

from netCDF4 import Dataset
import numpy as np
import sys

# Read in data
ncfile = './prov.seasonal.daily.powerden.ratio.nc'
fh = Dataset(ncfile, 'r')

cn_pwden_djf = fh.variables['cn_pwden_ratio_djf'][:]
cn_pwden_mam = fh.variables['cn_pwden_ratio_mam'][:]
cn_pwden_jja = fh.variables['cn_pwden_ratio_jja'][:]
cn_pwden_son = fh.variables['cn_pwden_ratio_son'][:]

prov_pwden_djf = fh.variables['prov_pwden_ratio_djf'][:][:]
prov_pwden_mam = fh.variables['prov_pwden_ratio_mam'][:][:]
prov_pwden_jja = fh.variables['prov_pwden_ratio_jja'][:][:]
prov_pwden_son = fh.variables['prov_pwden_ratio_son'][:][:]

prov = fh.variables['province'][:]

cn_pwden_all = np.concatenate((cn_pwden_djf,cn_pwden_mam,cn_pwden_jja,cn_pwden_son))

ind_worst = np.where( cn_pwden_all == cn_pwden_all.min() )
t_worst = ind_worst[0][0]
cn_worst = cn_pwden_all[t_worst]
print cn_worst
print ind_worst

prov_pwden_all = np.concatenate((prov_pwden_djf,prov_pwden_mam,prov_pwden_jja,prov_pwden_son),axis=1)
print prov_pwden_all.shape

prov_worst = np.array([prov_pwden_all[i,:].min() for i in range(31)])

print prov_worst

f = open("../plot/prov_worst_daily_ratio.txt","w+")
for i in range(31):
    if i != 0:
        f.write('\n')
    f.write(str(prov_worst[i]))
f.close()


sys.exit()

# Read in province name
f = open("../cn_province_nmlist.txt","r")
pro_nm_t = f.read().split("\n")
f.close()
pro_nm = pro_nm_t[:-1] # dump the last empty element

### Plot frequency in histogram
import matplotlib.pyplot as plt

bins_set = np.arange(-100,5,5)
xticks = np.arange(-2.5,-102.5,-5)
xticklim = np.arange(0,-105,-5)
xticklabels = []

for i in np.arange(len(xticks)):
    xticklabels.append( str(xticklim[i])+' ~ '+str(xticklim[i+1]) )
print xticklabels

yinterval = np.arange(0,1.05,0.05)

normed_value = 5

plt.figure()
plt.suptitle('China',fontsize=20)
plt.subplot(221)
hist, bins = np.histogram(cn_pwden_djf[:],bins=bins_set,density=True)
hist *= normed_value
maxnum = np.amax(hist)
print bins
widths = np.diff(bins)*0.8
plt.bar(bins[:-1]+5*0.2*0.5,hist,color='tomato',width=widths)
plt.xlim(0,-100)
plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
plt.yticks(fontsize=7)
plt.tick_params(axis='x',bottom='off')
plt.title('DJF',loc='left')

plt.subplot(222)
hist, bins = np.histogram(cn_pwden_mam[:],bins=bins_set,density=True)
hist *= normed_value
if maxnum < np.amax(hist):
    maxnum = np.amax(hist)
widths = np.diff(bins)*0.8
plt.bar(bins[:-1]+5*0.2*0.5,hist,color='tomato',width=widths)
plt.xlim(0,-100)
plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
plt.yticks(fontsize=7)
plt.tick_params(axis='x',bottom='off')
plt.title('MAM',loc='left')

plt.subplot(223)
hist, bins = np.histogram(cn_pwden_jja[:],bins=bins_set,density=True)
hist *= normed_value
if maxnum < np.amax(hist):
    maxnum = np.amax(hist)
widths = np.diff(bins)*0.8
plt.bar(bins[:-1]+5*0.2*0.5,hist,color='tomato',width=widths)
plt.xlim(0,-100)
plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
plt.yticks(fontsize=7)

plt.tick_params(axis='x',bottom='off')
plt.title('JJA',loc='left')

plt.subplot(224)
hist, bins = np.histogram(cn_pwden_son[:],bins=bins_set,density=True)
hist *= normed_value
if maxnum < np.amax(hist):
    maxnum = np.amax(hist)
widths = np.diff(bins)*0.8
plt.bar(bins[:-1]+5*0.2*0.5,hist,color='tomato',width=widths)
plt.xlim(0,-100)
plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
plt.yticks(fontsize=7)

plt.tick_params(axis='x',bottom='off')
plt.title('SON',loc='left')

# Create uniform limits for y-axis
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
plt.savefig('./Figures/season_freq/1Davg_ratio/China_1Davg_ratio_Freq_4seasons_aero_pwden',format='ps')
plt.close()

del hist,bins

sys.exit()

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
    plt.xlim(0,-100)
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
    plt.xlim(0,-100)
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
    plt.xlim(0,-100)
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
    plt.xlim(0,-100)
    plt.xticks(xticks,xticklabels,rotation='vertical',fontsize=7)
    plt.yticks(fontsize=7)

    plt.tick_params(axis='x',bottom='off')
    plt.title('SON',loc='left')

    # Create uniform limits for y-axis
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

