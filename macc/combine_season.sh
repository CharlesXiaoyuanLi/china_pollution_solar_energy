#! /bin/bash 

for year in {2003..2014}
do
    ncrcat *${year}12_monthly.nc *${year}01_monthly.nc *${year}02_monthly.nc macc_aero_${year}_DJF_monthly.nc
    ncra macc_aero_${year}_DJF_monthly.nc macc_aero_${year}_DJF.nc
    rm -f macc_aero_${year}_DJF_monthly.nc

    ncrcat *${year}03_monthly.nc *${year}04_monthly.nc *${year}05_monthly.nc macc_aero_${year}_MAM_monthly.nc
    ncra macc_aero_${year}_MAM_monthly.nc macc_aero_${year}_MAM.nc
    rm -f macc_aero_${year}_MAM_monthly.nc

    ncrcat *${year}06_monthly.nc *${year}07_monthly.nc *${year}08_monthly.nc macc_aero_${year}_JJA_monthly.nc
    ncra macc_aero_${year}_JJA_monthly.nc macc_aero_${year}_JJA.nc
    rm -f macc_aero_${year}_JJA_monthly.nc

    ncrcat *${year}09_monthly.nc *${year}10_monthly.nc *${year}11_monthly.nc macc_aero_${year}_SON_monthly.nc
    ncra macc_aero_${year}_SON_monthly.nc macc_aero_${year}_SON.nc
    rm -f macc_aero_${year}_SON_monthly.nc

    ncrcat macc_aero_${year}_DJF.nc macc_aero_${year}_MAM.nc macc_aero_${year}_JJA.nc macc_aero_${year}_SON.nc macc_aero_${year}_seasonly.nc
    
    rm -f macc_aero_${year}_DJF.nc macc_aero_${year}_MAM.nc macc_aero_${year}_JJA.nc macc_aero_${year}_SON.nc

done

ncrcat -v anthsrf,anthdfs *seasonly.nc macc_aero_anthsfc_2003_2014_seasonly.nc
