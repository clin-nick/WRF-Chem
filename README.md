# WRF-CHEM
![image](https://user-images.githubusercontent.com/61849582/229594261-5ca76c36-66ea-4f5f-b051-72c15008582e.png)

# Installation
## environment setting for running WRF-Chem in Kratos (WRF-Chem v4.0)
1. refer to /network/rit/lab/lulab/share_lib/bash_profile_clin (WRF-Chem section)
2. Cp /network/rit/lab/lulab/share_lib/WRF/HLC.TBL to your WRF run directory (e.g. /WRF/test/em_real/) 
3. Configure (Linux x86_64 options )
 -WRF configuring: 15 (./configure (-D (for debug)))
 -WPS configuring: 19
4.	Compile WRF
 -salloc -p kratos -N1 (-w kratos-10)
 -srun ./compile  em_real > log.compile 2>&1 &
 -exit (Releasing Resources)
