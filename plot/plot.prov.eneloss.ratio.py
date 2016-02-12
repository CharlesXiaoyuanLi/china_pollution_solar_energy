from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
import numpy as np
from matplotlib import cm, colors
import matplotlib as mpl
import sys

fig, ax = plt.subplots(frameon=False)
#ax = fig.add_subplot(111)

map = Basemap(llcrnrlon=80.,llcrnrlat=12.5,urcrnrlon=140.,urcrnrlat=52.5,resolution='i',projection='laea',lat_0=35.,
        lon_0=110.)

#map.drawmapboundary(fill_color='aqua')
#map.fillcontinents(color='#ddaa66') #lake_color='aqua'
#map.drawcoastlines()

map.readshapefile('../../maps_masks_shapefiles/china_provinces/china_provinces00/china_provinces00','cnprovinces', drawbounds = True)

patches = {}
subpatches = []
count = 0
for info, shape in zip(map.cnprovinces_info, map.cnprovinces):
    if count > 31:
        break
    if info['RINGNUM'] == 1:
        if info['SHAPENUM'] != 1: 
            patches[pronm] = subpatches
            subpatches = []
        pronm = info['EPROV']
        count += 1
    subpatches.append( Polygon(np.array(shape), True) )

print count

### Read in China Province Namelist                                                                                   
f = open("../cn_province_nmlist.txt","r")
pro_nm = f.read().split("\n")
f.close()

color = {}

### Read in data
pro_eneloss_ann = np.loadtxt("./prov_unit_eneloss_annualavg_daily_all_sky.txt")
pro_enebase_ann = np.loadtxt("./prov_unit_enebase_annualavg_daily_all_sky.txt")

pro_eneloss_ratio_ann = pro_eneloss_ann / pro_enebase_ann * 100


for i in range(len(pro_nm)-1):
    color[pro_nm[i]] = pro_eneloss_ratio_ann[i]
    if pro_nm[i] == 'Shanghai':
        color[pro_nm[i]] = color['Jiangsu']

cmap=cm.RdYlBu

bounds = np.linspace(-14,0,15)
norm = colors.BoundaryNorm(bounds, cmap.N)

#bounds = np.linspace(0,120,21)
#norm = colors.BoundaryNorm(bounds, ncolors=256)
#print norm(100)

#sys.exit()

for key, value in patches.iteritems():
    print key+" "+str(color[key])
    ax.add_collection(PatchCollection(value,
        facecolor=cm.RdYlBu(norm(color[key])),
        cmap=cm.RdYlBu, edgecolor='k', linewidth=1., zorder=2))

#    ax.add_collection(PatchCollection(value, facecolors=color[key], cmap=cm.RdYlBu, edgecolor='k',
#        linewidth=1., zorder=2, norm=norm))


### Plot colorbar
#set room for colorbar
fig.subplots_adjust(right = 0.85)


axcb = fig.add_axes([0.87,0.15,0.03,0.7])

cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, boundaries=bounds,
        orientation='vertical', ticks=bounds)

cb.set_label('%')

plt.savefig('prov_eneloss_ratio.ps',format='ps')

plt.show()
