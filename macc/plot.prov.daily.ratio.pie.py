### This program plots daily mean data by season.

from netCDF4 import Dataset
import numpy as np
import sys
import pandas as pd

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

bins_set = [-100,-40,-20,-5,0]
xticks = [-70,-30,-12.5,-2.5]
xticklim = bins_set[::-1]
lglabels = []

for i in np.arange(len(xticks)):
    lglabels.append( str(xticklim[i])+'% ~ '+str(xticklim[i+1])+'%' )
print lglabels

colors = ['lightskyblue', 'yellowgreen', 'lightcoral', 'gold']

normed_value = [60*100,20*100,15*100,5*100]

for province,iprov in zip(pro_nm,prov):
    print province
    if province == 'Shanghai': continue
    fig = plt.figure()
    plt.suptitle(province,fontsize=20)
    plt.subplot(221)
    hist, bins = np.histogram(prov_pwden_djf[iprov,:],bins=bins_set,density=True)
    pcts = hist * normed_value
    pctlabels = ['%.1f%%' % x for x in pcts[::-1]]
    hist, bins = np.histogram(prov_pwden_djf[iprov,:],bins=bins_set,density=False)
    numlabels = ['%dd' % x for x in hist[::-1]]
    labels = [y+'\n('+x+')' for x,y in zip(pctlabels,numlabels)]
    pplt = plt.pie(pcts[::-1], labels=labels, colors=colors,
            labeldistance=1.1)
    plt.axis('equal')
    ndays = len(prov_pwden_djf[iprov,:])
    #plt.title('DJF('+str(ndays)+'d)',loc='center')
    plt.title('DJF',loc='center')

    plt.subplot(222)
    hist, bins = np.histogram(prov_pwden_mam[iprov,:],bins=bins_set,density=True)
    pcts = hist * normed_value
    pctlabels = ['%.1f%%' % x for x in pcts[::-1]]
    hist, bins = np.histogram(prov_pwden_mam[iprov,:],bins=bins_set,density=False)
    numlabels = ['%dd' % x for x in hist[::-1]]
    labels = [y+'\n('+x+')' for x,y in zip(pctlabels,numlabels)]
    pplt = plt.pie(pcts[::-1], labels=labels, colors=colors,
            labeldistance=1.1)
    plt.axis('equal')
    ndays = len(prov_pwden_mam[iprov,:])
    #plt.title('MAM('+str(ndays)+'d)',loc='center')
    plt.title('MAM',loc='center')

    plt.subplot(223)
    hist, bins = np.histogram(prov_pwden_jja[iprov,:],bins=bins_set,density=True)
    pcts = hist * normed_value
    pctlabels = ['%.1f%%' % x for x in pcts[::-1]]
    hist, bins = np.histogram(prov_pwden_jja[iprov,:],bins=bins_set,density=False)
    numlabels = ['%dd' % x for x in hist[::-1]]
    labels = [y+'\n('+x+')' for x,y in zip(pctlabels,numlabels)]
    pplt = plt.pie(pcts[::-1], labels=labels, colors=colors,
            labeldistance=1.1)
    plt.axis('equal')
    ndays = len(prov_pwden_jja[iprov,:])
    #plt.title('JJA('+str(ndays)+'d)',loc='center')
    plt.title('JJA',loc='center')

    plt.subplot(224)
    hist, bins = np.histogram(prov_pwden_son[iprov,:],bins=bins_set,density=True)
    pcts = hist * normed_value
    pctlabels = ['%.1f%%' % x for x in pcts[::-1]]
    hist, bins = np.histogram(prov_pwden_son[iprov,:],bins=bins_set,density=False)
    numlabels = ['%dd' % x for x in hist[::-1]]
    labels = [y+'\n('+x+')' for x,y in zip(pctlabels,numlabels)]
    pplt = plt.pie(pcts[::-1], labels=labels, colors=colors,
            labeldistance=1.1)
    plt.axis('equal')  
    ndays = len(prov_pwden_son[iprov,:])
    #plt.title('SON('+str(ndays)+'d)',loc='center')
    plt.title('SON',loc='center')

    fig.legend(pplt[0], lglabels, loc="center", title="% Loss")

#    plt.tight_layout()
    plt.subplots_adjust(top=0.85) 
    plt.savefig('./Figures/season_freq/1Davg_pie/'+province+'_1Davg_ratio_pie_Freq_4seasons_aero_pwden.ps',format='ps')
    plt.close()

    del hist,bins

