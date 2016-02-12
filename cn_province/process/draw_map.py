import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from decimal import Decimal

# Create a 2D array representing lat x lon
nlat = 161
nlon = 320
cn_province_map = [[0 for x in range(nlon)] for x in range(nlat)]


# Import province boundaries
from subprocess import check_output

fnm = check_output("ls ../*txt", shell=True)

pid = 1

flist = open("cn_province_id_list.txt","w+")

for province in fnm.split():

    print province
    print pid

    x,y = [],[]
    f = open(province,'r')
    count = 0
    for l in f:
        for i in l.split():
            if count >= 5:
                if count % 2 != 0:
                    x.append(Decimal(i))
                else:
                    y.append(Decimal(i))
            count = count + 1

    f.close()

    del x[-5:-1]
    del x[-1]
    del y[-5:-1]
    del y[-1]

# Create Path from boundaries
    verts = []
    print len(x)
    for ind in range(len(x)):
        verts.append((x[ind],y[ind]))
    path = Path(verts)

# Assign value to grid point inside the boundaries
    grid_num = 0
    for i in range(nlat):
        for j in range(nlon):

            latn = 90 - i*1.125
            lonn = 0 + j*1.125

            if path.contains_point((lonn,latn)):
                cn_province_map[i][j] = pid
                grid_num += 1
                if pid == 2:
                    print str(latn) + "\t" + str(lonn)
    
    if not (province == "../heb1_out.txt" or province == "../heb2_out.txt"):   
        
        # print province namelist
        print province.replace("../","").replace("_out.txt","")
        flist.write(str(pid) + "\t" + province.replace("../","").replace("_out.txt","") + "\t" + str(grid_num) + "\n") 
        
        pid += 1

flist.close()

# Write cn_province_map to the file
f = open("cn_province_"+str(nlat)+"x"+str(nlon)+".txt","w+")
for i in range(nlat):
    for j in range(nlon):
        f.write(str(cn_province_map[i][j]) + " ")
    f.write("\n")
f.close()

#    if province == 'guangx':
#        fig = plt.figure()
#        ax = fig.add_subplot(111)
#    patch = patches.PathPatch(path, facecolor='orange', lw=2)
#    ax.add_patch(patch)

#ax.set_xlim(lonl,lonr)
#ax.set_ylim(latl,latu)
#plt.savefig('path',format='ps')
