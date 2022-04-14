      module PARAMS_MOD

      implicit none

      integer, parameter :: dp = selected_real_kind(14,300)

      real, parameter :: m2km = .001                          
      real, parameter :: ppm2vmr = 1.e-6                      
      real, parameter :: o2vmr = .2095                        
      real, parameter :: km2cm = 1.e5                         
      real, parameter :: m2s   = 60.                          

      REAL, PARAMETER :: pi = 3.1415926535898
      REAL, PARAMETER :: radius = 6.371E+3                    
      REAL, PARAMETER :: hc = 6.626068E-34 * 2.99792458E8
      REAL, PARAMETER :: largest=1.E+36

      REAL, PARAMETER :: pzero = +10./largest
      REAL, PARAMETER :: nzero = -10./largest

      REAL, PARAMETER :: precis = 1.e-7

      real :: lambda_cutoff                        

      end module PARAMS_MOD
