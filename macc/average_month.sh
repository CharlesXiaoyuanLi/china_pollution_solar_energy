#! /bin/bash 

for fn in $(ls |grep aerosol_rf)
do
    ncra ${fn} ${fn/.nc/_monthly.nc}
done
