      module module_params

      implicit none





      INTEGER, PARAMETER :: kout=53

      INTEGER, PARAMETER :: kin=12



      integer, PARAMETER :: kz=125

      integer, PARAMETER :: kw=1000

      integer, PARAMETER :: kt=100



      integer, PARAMETER :: ks=60

      integer, PARAMETER :: kj=150

      integer, PARAMETER :: kdom=200

      real, PARAMETER :: deltax = 1.E-5




      real, PARAMETER :: pi=3.1415926535898


      real, PARAMETER :: radius=6.371E+3


      real, PARAMETER :: hc = 6.626068E-34 * 2.99792458E8


      real, PARAMETER :: largest=1.E+36


      real, PARAMETER :: pzero = +10./largest
      real, PARAMETER :: nzero = -10./largest


      real, PARAMETER :: precis = 1.e-7

      end module module_params
