      module module_data_mosaic_main








      use module_data_mosaic_kind, only:  r8
      use module_data_mosaic_constants, only:  &
          avogad, deg2rad, pi, piover4, piover6, third


      implicit none

      integer, parameter ::   &
                ngas_com = 40,   &
                ngas_urb = 19,   &
                ngas_bio =  7,   &
                ngas_mar = 11
      
      integer, parameter ::   &
                  naer_tot = 24 		      

      integer, save ::   &
                naerbin  = -999888777         
      




      integer, parameter ::   &
                ncld_tot = 13,		   &  
                ncldbin  =  4,		   &  
                ncld     = 22		

      integer, parameter :: ngas_max = ngas_com + ngas_urb + ngas_bio + ngas_mar
      
      integer, parameter :: ncld_max = ncld_tot*ncldbin
      
      integer, save :: naer_max = -999888777  

      integer, save :: ntot_max = -999888777  
      

      integer, save ::   &
      		naerbin_used=0,   &   
      		ncldbin_used=0,   &   
      		ntot_used=ngas_max    

      integer, save ::   &
      		ipmcmos = 0        
      		                   

      real(r8), parameter :: press0_pa = 1.01325d5  
      real(r8), parameter :: mw_air = 28.966d0      




      integer, save ::	m_partmc_mosaic  

      integer, save ::	mgas, maer, mcld

      integer, save ::	maeroptic, mshellcore

      integer, save ::	msolar, mphoto




      end module module_data_mosaic_main
