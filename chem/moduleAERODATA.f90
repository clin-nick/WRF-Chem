      MODULE AERODATA

      INTEGER, PARAMETER :: NASECT = 3  
                                        
                                        


      REAL, PARAMETER :: DPUP = 10.0


      REAL, PARAMETER :: DPLOW = 0.0215

      REAL VRAT    

      REAL VRLOW   

      REAL VRHI    

      REAL DPBINMIN
                   





      REAL, PARAMETER :: DENSP = 1.352


      REAL :: DPCTR( NASECT )


      REAL :: PMCONC( NASECT )



      REAL :: SURFP, VOL, AEROMA, XNUM

      REAL :: AREA( NASECT )    

      END MODULE AERODATA

