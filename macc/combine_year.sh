#! /bin/bash 

for year in {2003..2014}
do
    ncrcat *${year}*.nc macc_aero_${year}_daily.nc
    ncks -v anthsrf,anthdfs macc_aero_${year}_daily.nc macc_aero_anthsfc_${year}_daily.nc
    rm -f macc_aero_${year}_daily.nc
done
