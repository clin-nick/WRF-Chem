# WRF-CHEM
![image](https://user-images.githubusercontent.com/61849582/229594261-5ca76c36-66ea-4f5f-b051-72c15008582e.png)

# Installation
## Environment setting for running WRF-Chem in Kratos (WRF-Chem v4.0)
1. refer to /network/rit/lab/lulab/share_lib/bash_profile_clin (WRF-Chem section)
2. Cp /network/rit/lab/lulab/share_lib/WRF/HLC.TBL to your WRF run directory (e.g. /WRF/test/em_real/) 
3. Configure (Linux x86_64 options )
- WRF configuring: 15 (./configure (-D (for debug)))
- WPS configuring: 19
4.	Compile WRF
- salloc -p kratos -N1 (-w kratos-10)
- srun ./compile  em_real > log.compile (2>&1 &)
- exit (Releasing Resources)
## Reference
1.	WRF downloads
https://www2.mmm.ucar.edu/wrf/users/downloads.html ,
https://github.com/wrf-model
2.	Latest WRF updates
https://github.com/wrf-model/WRF/releases
3.	WRF online tutorial (environment)
https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php#STEP4
4.	WRF online tutorial (install and practice)
https://www2.mmm.ucar.edu/wrf/OnLineTutorial/Compile/arw_compile3.php
5.	Best practice of namelist.input (WRF) https://www2.mmm.ucar.edu/wrf/users/namelist_best_prac_wrf.html
6.	Best practice of namelist.wps (WPS) https://www2.mmm.ucar.edu/wrf/users/namelist_best_prac_wps.html
7.	ASRC wiki 
https://wiki.albany.edu/display/asrc/%27Kratos%27+cluster
8.	ASRC IT support
http://support.asrc.albany.edu/
# Steps for WRF-Chem Run
## Run WPS
1. Prepare MET IC/BC Input:
- GFS:  https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs ,
- NAM:  https://rda.ucar.edu/datasets/ds609.0/
2. Run geogrid
- Edit the configuration file “namelist.wps”
- Check domain: 
ncl util/plotgrids_new.ncl 
- Run geogrid.exe
