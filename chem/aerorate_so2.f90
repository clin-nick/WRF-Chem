
      SUBROUTINE AERORATE_SO2( BTEMP, BPRESS, RTDAT_AE ) 


      USE AERODATA
      USE HETAERO


      IMPLICIT NONE


      REAL BTEMP     
      REAL BPRESS    
      REAL RTDAT_AE ( NRXNAERO )  
                                  
                                  



      INTEGER :: INASEC, IRXN  

      REAL VSP                  
      REAL DG                   
                                
      REAL, SAVE :: GAMMA( NRXNAERO )    

      REAL RS_TOT               
                                
      REAL RS                   
                                
      REAL TOTMA                
      REAL RCL                  
      REAL*8      PI 
      PARAMETER ( PI = 3.14159265358979324 )

      LOGICAL, SAVE :: FIRSTIME = .TRUE.




      IF ( IAERORATE == 0 ) RETURN  
                                    
      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.








         IF ( NGAMMA == 1 ) THEN        
            GAMMA( ISO2 )  = 1.0E-4

         ELSE IF ( NGAMMA == 2 ) THEN   
            GAMMA( ISO2 )  = 1.0E-5

         ELSE IF ( NGAMMA == 3 ) THEN   
            GAMMA( ISO2 )  = 0.1
         END IF

      END IF  




      DO INASEC  =  1,  NASECT

         TOTMA = 0.0






         TOTMA = PMCONC( INASEC )


         SURFP  = 4.0 * PI * ( DPCTR( INASEC ) / 2.0 )**2.0  
         VOL    = ( 4.0 * PI / 3.0 ) * ( DPCTR( INASEC )  &
                             / 2.0 )**3.0  
         AEROMA = VOL * DENSP * 1.0E-6 

         IF ( AEROMA > 0.0 ) THEN
            XNUM  = TOTMA / AEROMA
            AREA( INASEC ) = SURFP * XNUM * 1.0E-12               
         ELSE
            AREA( INASEC ) = 0.0
         END IF

      END DO



      DO IRXN = 1, NRXNAERO
         RS_TOT = 0.0
         DG = ( DG0( IRXN ) * ( PRESS0 / BPRESS ) * &
               ( BTEMP / TEMP0 )**1.75 ) * 1.0E-4
         VSP = SQRT( 8.0 * RG * BTEMP / PI / XMOLWEI( IRXN ) )
         DO INASEC = 1, NASECT
            RCL = ( DPCTR( INASEC ) / 2.0) * 1.0E-6
            RS  = 1.0 / ( RCL / DG + 4.0 / VSP / &
                          GAMMA( IRXN ) ) * AREA( INASEC )

            RS_TOT = RS_TOT + RS
         END DO



         
            RTDAT_AE( IRXN ) = RS_TOT

      END DO

      RETURN
      END
