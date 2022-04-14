



















      SUBROUTINE ISRP1R (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)



      CALL INIT1 (WI, RHI, TEMPI)



      IF (RH.GE.DRNH42S4) THEN         
         SULRATW = GETASR(WAER(2), RHI)     
      ELSE
         SULRATW = 2.0D0                    
      ENDIF
      SULRAT  = WAER(3)/WAER(2)         





      IF (SULRATW.LE.SULRAT) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'K2'
         CALL CALCK2                 
      ELSE

         IF (RH.LT.DRNH42S4) THEN
            SCASE = 'K1'
            CALL CALCK1              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'K2'
            CALL CALCK2              
         ENDIF
      ENDIF



      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN
      W(2) = WAER(2)
      W(3) = WAER(3)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL CALCB4                 
         SCASE = 'L4'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'B1'
            CALL CALCB1              
            SCASE = 'L1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN
            SCASE = 'B2'
            CALL CALCB2              
            SCASE = 'L2'

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'B3'
            CALL CALCB3              
            SCASE = 'L3'

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'B4'
            CALL CALCB4              
            SCASE = 'L4'
         ENDIF
      ENDIF

      CALL CALCNH3P          



      ELSEIF (SULRAT.LT.1.0) THEN
      W(2) = WAER(2)
      W(3) = WAER(3)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC2                 
         SCASE = 'M2'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'C1'
            CALL CALCC1              
            SCASE = 'M1'

         ELSEIF (DRNH4HS4.LE.RH) THEN
            SCASE = 'C2'
            CALL CALCC2              
            SCASE = 'M2'
         ENDIF
      ENDIF

      CALL CALCNH3P

      ENDIF
      RETURN



      END SUBROUTINE ISRP1R                 
















      SUBROUTINE ISRP2R (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ



      TRYLIQ = .TRUE.             

10    CALL INIT2 (WI, RHI, TEMPI)



      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN 
         SULRATW = GETASR(WAER(2), RHI)     
      ELSE
         SULRATW = 2.0D0                    
      ENDIF
      SULRAT = WAER(3)/WAER(2)





      IF (SULRATW.LE.SULRAT) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'N3'
         CALL CALCN3                 
      ELSE

         IF (RH.LT.DRNH4NO3) THEN
            SCASE = 'N1'
            CALL CALCN1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'N2'
            CALL CALCN2              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'N3'
            CALL CALCN3              
         ENDIF
      ENDIF







      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN
      W(2) = WAER(2)
      W(3) = WAER(3)
      W(4) = WAER(4)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL CALCB4                 
         SCASE = 'O4'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'B1'
            CALL CALCB1              
            SCASE = 'O1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN
            SCASE = 'B2'
            CALL CALCB2              
            SCASE = 'O2'

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'B3'
            CALL CALCB3              
            SCASE = 'O3'

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'B4'
            CALL CALCB4              
            SCASE = 'O4'
         ENDIF
      ENDIF



      MOLAL(7) = WAER(4)             
      MOLAL(1) = MOLAL(1) + WAER(4)  
      CALL CALCNAP            
      CALL CALCNH3P







      ELSEIF (SULRAT.LT.1.0) THEN
      W(2) = WAER(2)
      W(3) = WAER(3)
      W(4) = WAER(4)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC2                 
         SCASE = 'P2'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'C1'
            CALL CALCC1              
            SCASE = 'P1'

         ELSEIF (DRNH4HS4.LE.RH) THEN
            SCASE = 'C2'
            CALL CALCC2              
            SCASE = 'P2'
         ENDIF
      ENDIF



      MOLAL(7) = WAER(4)             
      MOLAL(1) = MOLAL(1) + WAER(4)  

      CALL CALCNAP                   
      CALL CALCNH3P
      ENDIF



      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0   &
                                          .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF

      RETURN



      END SUBROUTINE ISRP2R                 















      SUBROUTINE ISRP3R (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ








      TRYLIQ = .TRUE.             

10    CALL ISOINIT3 (WI, RHI, TEMPI) 











      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN  
         FRSO4   = WAER(2) - WAER(1)/2.0D0     
         FRSO4   = MAX(FRSO4, TINY)
         SRI     = GETASR(FRSO4, RHI)          
         SULRATW = (WAER(1)+FRSO4*SRI)/WAER(2) 
         SULRATW = MIN (SULRATW, 2.0D0)
      ELSE
         SULRATW = 2.0D0                     
      ENDIF
      SULRAT = (WAER(1)+WAER(3))/WAER(2)
      SODRAT = WAER(1)/WAER(2)





      IF (SULRATW.LE.SULRAT .AND. SODRAT.LT.2.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'Q5'
         CALL CALCQ5                 
         SCASE = 'Q5'
      ELSE

         IF (RH.LT.DRNH4NO3) THEN
            SCASE = 'Q1'
            CALL CALCQ1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH4CL) THEN
            SCASE = 'Q2'
            CALL CALCQ2              

         ELSEIF (DRNH4CL.LE.RH  .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'Q3'
            CALL CALCQ3              

        ELSEIF (DRNH42S4.LE.RH  .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'Q4'
            CALL CALCQ4              
            SCASE = 'Q4'

         ELSEIF (DRNA2SO4.LE.RH) THEN
            SCASE = 'Q5'
            CALL CALCQ5              
            SCASE = 'Q5'
         ENDIF
      ENDIF



      ELSE IF (SULRAT.GE.SULRATW .AND. SODRAT.GE.2.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'R6'
         CALL CALCR6                 
         SCASE = 'R6'
      ELSE

         IF (RH.LT.DRNH4NO3) THEN
            SCASE = 'R1'
            CALL CALCR1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN
            SCASE = 'R2'
            CALL CALCR2              

         ELSEIF (DRNANO3.LE.RH  .AND. RH.LT.DRNACL) THEN
            SCASE = 'R3'
            CALL CALCR3              

         ELSEIF (DRNACL.LE.RH   .AND. RH.LT.DRNH4CL) THEN
            SCASE = 'R4'
            CALL CALCR4              

         ELSEIF (DRNH4CL.LE.RH .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'R5'
            CALL CALCR5              
            SCASE = 'R5'

         ELSEIF (DRNA2SO4.LE.RH) THEN
            SCASE = 'R6'
            CALL CALCR6              
            SCASE = 'R6'
         ENDIF
      ENDIF



      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN
      DO 100 I=1,NCOMP
         W(I) = WAER(I)
100   CONTINUE

      IF(METSTBL.EQ.1) THEN
         SCASE = 'I6'
         CALL CALCI6                 
         SCASE = 'S6'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'I1'
            CALL CALCI1              
            SCASE = 'S1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN
            SCASE = 'I2'
            CALL CALCI2              
            SCASE = 'S2'

         ELSEIF (DRNAHSO4.LE.RH .AND. RH.LT.DRLC) THEN
            SCASE = 'I3'
            CALL CALCI3              
            SCASE = 'S3'

         ELSEIF (DRLC.LE.RH     .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'I4'
            CALL CALCI4              
            SCASE = 'S4'

         ELSEIF (DRNH42S4.LE.RH .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'I5'
            CALL CALCI5              
            SCASE = 'S5'

         ELSEIF (DRNA2SO4.LE.RH) THEN
            SCASE = 'I6'
            CALL CALCI6              
            SCASE = 'S6'
         ENDIF
      ENDIF

      CALL CALCNHP                
      CALL CALCNH3P



      ELSEIF (SULRAT.LT.1.0) THEN
      DO 200 I=1,NCOMP
         W(I) = WAER(I)
200   CONTINUE

      IF(METSTBL.EQ.1) THEN
         SCASE = 'J3'
         CALL CALCJ3                 
         SCASE = 'T3'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'J1'
            CALL CALCJ1              
            SCASE = 'T1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN
            SCASE = 'J2'
            CALL CALCJ2              
            SCASE = 'T2'

         ELSEIF (DRNAHSO4.LE.RH) THEN
            SCASE = 'J3'
            CALL CALCJ3
            SCASE = 'T3'
         ENDIF
      ENDIF

      CALL CALCNHP                
      CALL CALCNH3P

      ENDIF




      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0   &
                            .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF

      RETURN



      END SUBROUTINE ISRP3R                 
















      SUBROUTINE CALCK2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION NH4I, NH3GI, NH3AQ



      CALAOU   =.TRUE.     
      FRST     =.TRUE.
      CALAIN   =.TRUE.



      MOLALR(4)= MIN(WAER(2), 0.5d0*WAER(3))
      WATER    = MOLALR(4)/M0(4)  



      DO 10 I=1,NSWEEP

         A2   = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
         AKW  = XKW *RH*WATER*WATER

         NH4I = WAER(3)
         SO4I = WAER(2)
         HSO4I= ZERO

         CALL CALCPH (2.D0*SO4I - NH4I, HI, OHI)    

         NH3AQ = ZERO                               
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ (NH4I, OHI, DEL)
            NH4I  = MAX (NH4I-DEL, ZERO)
            OHI   = MAX (OHI -DEL, TINY)
            NH3AQ = DEL
            HI    = AKW/OHI
         ENDIF

         CALL CALCHS4 (HI, SO4I, ZERO, DEL)         
         SO4I  = SO4I - DEL
         HI    = HI   - DEL
         HSO4I = DEL

         NH3GI = NH4I/HI/A2   



         MOLAL(1) = HI
         MOLAL(3) = NH4I
         MOLAL(5) = SO4I
         MOLAL(6) = HSO4I
         COH      = OHI
         GASAQ(1) = NH3AQ
         GNH3     = NH3GI



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE

20    RETURN



      END SUBROUTINE CALCK2





















      SUBROUTINE CALCK1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      CNH42S4 = MIN(WAER(2),0.5d0*WAER(3))  
      GNH3    = ZERO

      W(2)    = CNH42S4
      W(3)    = 2.D0*CNH42S4 + GNH3

      RETURN



      END SUBROUTINE CALCK1


















      SUBROUTINE CALCN3
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION NH4I, NO3I, NH3AQ, NO3AQ






      CALAOU =.TRUE.              
      FRST   =.TRUE.
      CALAIN =.TRUE.



      MOLALR(4) = MIN(WAER(2),0.5d0*WAER(3))       
      AML5      = MAX(WAER(3)-2.D0*MOLALR(4),ZERO) 
      MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO)     
      WATER     = MOLALR(4)/M0(4) + MOLALR(5)/M0(5)
      WATER     = MAX(WATER, TINY)



      DO 10 I=1,NSWEEP
         A2    = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.

         A3    = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4    = XK7*(WATER/GAMA(4))**3.0
         AKW   = XKW *RH*WATER*WATER



         NH4I  = WAER(3)
         NO3I  = WAER(4)
         SO4I  = WAER(2)
         HSO4I = ZERO

         CALL CALCPH (2.D0*SO4I + NO3I - NH4I, HI, OHI)



         NH3AQ = ZERO
         NO3AQ = ZERO
         GG    = 2.D0*SO4I + NO3I - NH4I
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
         ELSE
            HI    = ZERO
            CALL CALCNIAQ2 (GG, NO3I, HI, NO3AQ) 



            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
         ENDIF



         MOLAL (1) = HI
         MOLAL (3) = NH4I
         MOLAL (5) = SO4I
         MOLAL (6) = HSO4I
         MOLAL (7) = NO3I
         COH       = OHI

         CNH42S4   = ZERO
         CNH4NO3   = ZERO

         GASAQ(1)  = NH3AQ
         GASAQ(3)  = NO3AQ

         GHNO3     = HI*NO3I/A3
         GNH3      = NH4I/HI/A2   



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    RETURN



      END SUBROUTINE CALCN3

















      SUBROUTINE CALCN2
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CHI1   = MIN(WAER(2),0.5d0*WAER(3))     
      CHI2   = MAX(WAER(3) - 2.D0*CHI1, ZERO) 
      CHI3   = MAX(WAER(4) - CHI2, ZERO)      

      PSI2   = CHI2
      PSI3   = CHI3

      CALAOU = .TRUE.              
      PSI1LO = TINY                
      PSI1HI = CHI1                



      X1 = PSI1HI
      Y1 = FUNCN2 (X1)
      IF (Y1.LE.EPS) RETURN   
      YHI= Y1                 



      DX = (PSI1HI-PSI1LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, ZERO)
         Y2 = FUNCN2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (ABS(Y2) .LT. EPS) THEN   
         RETURN



      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = CHI4
         YY = FUNCN2(P4)
         GOTO 50



      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         P4 = TINY
         YY = FUNCN2(P4)
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCN2')    
         RETURN
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCN2 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCN2')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCN2 (X3)
50    CONTINUE
      RETURN



      END SUBROUTINE CALCN2













      DOUBLE PRECISION FUNCTION FUNCN2 (P1)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION NH4I, NO3I, NH3AQ, NO3AQ






      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI1   = P1



      DO 10 I=1,NSWEEP
         A2    = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.

         A3    = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4    = XK7*(WATER/GAMA(4))**3.0
         AKW   = XKW *RH*WATER*WATER



         NH4I  = 2.D0*PSI1 + PSI2
         NO3I  = PSI2 + PSI3
         SO4I  = PSI1
         HSO4I = ZERO

         CALL CALCPH (2.D0*SO4I + NO3I - NH4I, HI, OHI)



         NH3AQ = ZERO
         NO3AQ = ZERO
         GG    = 2.D0*SO4I + NO3I - NH4I
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
         ELSE
            HI    = ZERO
            CALL CALCNIAQ2 (GG, NO3I, HI, NO3AQ) 



            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
         ENDIF



         MOLAL (1) = HI
         MOLAL (3) = NH4I
         MOLAL (5) = SO4I
         MOLAL (6) = HSO4I
         MOLAL (7) = NO3I
         COH       = OHI

         CNH42S4   = CHI1 - PSI1
         CNH4NO3   = ZERO

         GASAQ(1)  = NH3AQ
         GASAQ(3)  = NO3AQ

         GHNO3     = HI*NO3I/A3
         GNH3      = NH4I/HI/A2   



         CALL CALCMR



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    FUNCN2= NH4I*NH4I*SO4I/A4 - ONE
      RETURN



      END FUNCTION FUNCN2     





















      SUBROUTINE CALCN1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCN1A, CALCN2



      IF (RH.LT.DRMASAN) THEN
         SCASE = 'N1 ; SUBCASE 1'
         CALL CALCN1A              
         SCASE = 'N1 ; SUBCASE 1'
      ELSE
         SCASE = 'N1 ; SUBCASE 2'
         CALL CALCMDRP (RH, DRMASAN, DRNH4NO3, CALCN1A, CALCN2)
         SCASE = 'N1 ; SUBCASE 2'
      ENDIF

      RETURN



      END SUBROUTINE CALCN1




















      SUBROUTINE CALCN1A
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









      PSI1    = WAER(4)




      PSI2    = MAX(MIN(WAER(2),0.5d0*(WAER(3)-PSI1)),TINY)

      CNH4NO3 = PSI1
      CNH42S4 = PSI2



      GNH3    = ZERO
      GHNO3   = ZERO

      W(2)    = PSI2
      W(3)    = GNH3  + PSI1 + 2.0*PSI2
      W(4)    = GHNO3 + PSI1

      RETURN



      END SUBROUTINE CALCN1A

















      SUBROUTINE CALCQ5
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ






      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.



      CALL CALCQ1A

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4

      CALL CALCMR           

      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP
      AKW = XKW*RH*WATER*WATER               



      NAI    = WAER(1)
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4

      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ

      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNACL   = ZERO
      CNANO3  = ZERO
      CNA2SO4 = ZERO

      RETURN



      END SUBROUTINE CALCQ5
















      SUBROUTINE CALCQ4
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      LOGICAL PSCONV1
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ






      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV1 =.TRUE.
      PSI1O   =-GREAT
      ROOT3   = ZERO



      CALL CALCQ1A

      CHI1   = CNA2SO4      

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4

      CALL CALCMR           

      NAI    = WAER(1)      
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP
      A5  = XK5 *(WATER/GAMA(2))**3.         
      AKW = XKW*RH*WATER*WATER               



      IF (NAI*NAI*SO4I .GT. A5) THEN
         BB =-(WAER(2) + WAER(1))
         CC = WAER(1)*WAER(2) + 0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*WAER(2) - A5)
         CALL POLY3(BB, CC, DD, ROOT3, ISLV)
         IF (ISLV.NE.0) ROOT3 = TINY
         ROOT3 = MIN (ROOT3, WAER(1)/2.0, WAER(2), CHI1)
         ROOT3 = MAX (ROOT3, ZERO)
         PSI1  = CHI1-ROOT3
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      NAI = WAER(1) - 2.D0*ROOT3
      SO4I= WAER(2) - ROOT3
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         IF (PSCONV1) GOTO 20
      ENDIF
10    CONTINUE




20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4

      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ

      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNACL   = ZERO
      CNANO3  = ZERO
      CNA2SO4 = CHI1 - PSI1

      RETURN



      END SUBROUTINE CALCQ4

















      SUBROUTINE CALCQ3
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL EXNO, EXCL
      EXTERNAL CALCQ1A, CALCQ4



      EXNO = WAER(4).GT.TINY
      EXCL = WAER(5).GT.TINY

      IF (EXNO .OR. EXCL) THEN             
         SCASE = 'Q3 ; SUBCASE 1'
         CALL CALCQ3A
         SCASE = 'Q3 ; SUBCASE 1'

      ELSE                                 
         IF (RH.LT.DRMG3) THEN
            SCASE = 'Q3 ; SUBCASE 2'
            CALL CALCQ1A             
            SCASE = 'Q3 ; SUBCASE 2'
         ELSE
            SCASE = 'Q3 ; SUBCASE 3' 
            CALL CALCMDRP (RH, DRMG3, DRNH42S4, CALCQ1A, CALCQ4)
            SCASE = 'Q3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCQ3



















      SUBROUTINE CALCQ3A
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      LOGICAL PSCONV1, PSCONV6
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ






      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV1 =.TRUE.
      PSCONV6 =.TRUE.

      PSI1O   =-GREAT
      PSI6O   =-GREAT

      ROOT1   = ZERO
      ROOT3   = ZERO



      CALL CALCQ1A

      CHI1   = CNA2SO4      
      CHI4   = CNH4CL
      CHI6   = CNH42S4

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4

      CALL CALCMR           

      NAI    = WAER(1)      
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP
      A5  = XK5 *(WATER/GAMA(2))**3.         
      A7  = XK7 *(WATER/GAMA(4))**3.         
      AKW = XKW*RH*WATER*WATER               



      IF (NAI*NAI*SO4I .GT. A5) THEN
         BB =-(WAER(2) + WAER(1) - ROOT1)
         CC = WAER(1)*(WAER(2) - ROOT1) + 0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*(WAER(2) - ROOT1) - A5)
         CALL POLY3(BB, CC, DD, ROOT3, ISLV)
         IF (ISLV.NE.0) ROOT3 = TINY
         ROOT3 = MIN (ROOT3, WAER(1)/2.0, WAER(2) - ROOT1, CHI1)
         ROOT3 = MAX (ROOT3, ZERO)
         PSI1  = CHI1-ROOT3
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*NH4I*SO4I .GT. A4) THEN
         BB =-(WAER(2)+WAER(3)-ROOT3)
         CC =  WAER(3)*(WAER(2)-ROOT3+0.5D0*WAER(3))
         DD =-((WAER(2)-ROOT3)*WAER(3)**2.D0 + A4)/4.D0
         CALL POLY3(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MIN(ROOT1, WAER(3), WAER(2)-ROOT3, CHI6)
         ROOT1 = MAX(ROOT1, ZERO)
         PSI6  = CHI6-ROOT1
      ENDIF
      PSCONV6 = ABS(PSI6-PSI6O) .LE. EPS*PSI6O
      PSI6O   = PSI6



      NAI = WAER(1) - 2.D0*ROOT3
      SO4I= WAER(2) - ROOT1 - ROOT3
      NH4I= WAER(3) - 2.D0*ROOT1
      NO3I= WAER(4)
      CLI = WAER(5)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         IF (PSCONV1 .AND. PSCONV6) GOTO 20
      ENDIF
10    CONTINUE




20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4

      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ

      CNH42S4 = CHI6 - PSI6
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNACL   = ZERO
      CNANO3  = ZERO
      CNA2SO4 = CHI1 - PSI1

      RETURN



      END SUBROUTINE CALCQ3A

















      SUBROUTINE CALCQ2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL EXNO, EXCL
      EXTERNAL CALCQ1A, CALCQ3A, CALCQ4



      EXNO = WAER(4).GT.TINY
      EXCL = WAER(5).GT.TINY

      IF (EXNO) THEN                       
         SCASE = 'Q2 ; SUBCASE 1'
         CALL CALCQ2A
         SCASE = 'Q2 ; SUBCASE 1'

      ELSEIF (.NOT.EXNO .AND. EXCL) THEN   
         IF (RH.LT.DRMG2) THEN
            SCASE = 'Q2 ; SUBCASE 2'
            CALL CALCQ1A             
            SCASE = 'Q2 ; SUBCASE 2'
         ELSE
            SCASE = 'Q2 ; SUBCASE 3' 
            CALL CALCMDRP (RH, DRMG2, DRNH4CL, CALCQ1A, CALCQ3A)
            SCASE = 'Q2 ; SUBCASE 3'
         ENDIF

      ELSE                                 
         IF (RH.LT.DRMG3) THEN
            SCASE = 'Q2 ; SUBCASE 2'
            CALL CALCQ1A             
            SCASE = 'Q2 ; SUBCASE 2'
         ELSE
            SCASE = 'Q2 ; SUBCASE 4' 
            CALL CALCMDRP (RH, DRMG3, DRNH42S4, CALCQ1A, CALCQ4)
            SCASE = 'Q2 ; SUBCASE 4'
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCQ2


















      SUBROUTINE CALCQ2A
      USE SOLUT
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      LOGICAL PSCONV1, PSCONV4, PSCONV6
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ






      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV1 =.TRUE.
      PSCONV4 =.TRUE.
      PSCONV6 =.TRUE.

      PSI1O   =-GREAT
      PSI4O   =-GREAT
      PSI6O   =-GREAT

      ROOT1   = ZERO
      ROOT2   = ZERO
      ROOT3   = ZERO



      CALL CALCQ1A

      CHI1   = CNA2SO4      
      CHI4   = CNH4CL
      CHI6   = CNH42S4

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4

      CALL CALCMR           

      NAI    = WAER(1)      
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP
      A5  = XK5 *(WATER/GAMA(2))**3.         
      A14 = XK14*(WATER/GAMA(6))**2.         
      A7  = XK7 *(WATER/GAMA(4))**3.         
      AKW = XKW*RH*WATER*WATER               



      IF (NH4I*CLI .GT. A14) THEN
         BB    =-(WAER(3) + WAER(5) - 2.D0*ROOT1)
         CC    = WAER(5)*(WAER(3) - 2.D0*ROOT1) - A14
         DD    = BB*BB - 4.D0*CC
         IF (DD.LT.ZERO) THEN
            ROOT2 = ZERO
         ELSE
            DD    = SQRT(DD)
            ROOT2A= 0.5D0*(-BB+DD)
            ROOT2B= 0.5D0*(-BB-DD)
            IF (ZERO.LE.ROOT2A) THEN
               ROOT2 = ROOT2A
            ELSE
               ROOT2 = ROOT2B
            ENDIF
            ROOT2 = MIN(ROOT2, WAER(5), WAER(3) - 2.D0*ROOT1, CHI4)
            ROOT2 = MAX(ROOT2, ZERO)
            PSI4  = CHI4 - ROOT2
         ENDIF
      ENDIF
      PSCONV4 = ABS(PSI4-PSI4O) .LE. EPS*PSI4O
      PSI4O   = PSI4



      IF (NAI*NAI*SO4I .GT. A5) THEN
         BB =-(WAER(2) + WAER(1) - ROOT1)
         CC = WAER(1)*(WAER(2) - ROOT1) + 0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*(WAER(2) - ROOT1) - A5)
         CALL POLY3(BB, CC, DD, ROOT3, ISLV)
         IF (ISLV.NE.0) ROOT3 = TINY
         ROOT3 = MIN (ROOT3, WAER(1)/2.0, WAER(2) - ROOT1, CHI1)
         ROOT3 = MAX (ROOT3, ZERO)
         PSI1  = CHI1-ROOT3
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*NH4I*SO4I .GT. A4) THEN
         BB =-(WAER(2)+WAER(3)-ROOT2-ROOT3)
         CC = (WAER(3)-ROOT2)*(WAER(2)-ROOT3+0.5D0*(WAER(3)-ROOT2))
         DD =-((WAER(2)-ROOT3)*(WAER(3)-ROOT2)**2.D0 + A4)/4.D0
         CALL POLY3(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MIN(ROOT1, WAER(3)-ROOT2, WAER(2)-ROOT3, CHI6)
         ROOT1 = MAX(ROOT1, ZERO)
         PSI6  = CHI6-ROOT1
      ENDIF
      PSCONV6 = ABS(PSI6-PSI6O) .LE. EPS*PSI6O
      PSI6O   = PSI6



      NAI = WAER(1) - 2.D0*ROOT3
      SO4I= WAER(2) - ROOT1 - ROOT3
      NH4I= WAER(3) - ROOT2 - 2.D0*ROOT1
      NO3I= WAER(4)
      CLI = WAER(5) - ROOT2



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         IF (PSCONV1 .AND. PSCONV4 .AND. PSCONV6) GOTO 20
      ENDIF
10    CONTINUE




20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4

      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ

      CNH42S4 = CHI6 - PSI6
      CNH4NO3 = ZERO
      CNH4CL  = CHI4 - PSI4
      CNACL   = ZERO
      CNANO3  = ZERO
      CNA2SO4 = CHI1 - PSI1

      RETURN



      END SUBROUTINE CALCQ2A

















      SUBROUTINE CALCQ1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL EXNO, EXCL
      EXTERNAL CALCQ1A, CALCQ2A, CALCQ3A, CALCQ4



      EXNO = WAER(4).GT.TINY
      EXCL = WAER(5).GT.TINY

      IF (EXNO .AND. EXCL) THEN           
         IF (RH.LT.DRMG1) THEN
            SCASE = 'Q1 ; SUBCASE 1'
            CALL CALCQ1A             
            SCASE = 'Q1 ; SUBCASE 1'
         ELSE
            SCASE = 'Q1 ; SUBCASE 2' 
            CALL CALCMDRP (RH, DRMG1, DRNH4NO3, CALCQ1A, CALCQ2A)
            SCASE = 'Q1 ; SUBCASE 2'
         ENDIF

      ELSE IF (EXNO .AND. .NOT.EXCL) THEN 
         IF (RH.LT.DRMQ1) THEN
            SCASE = 'Q1 ; SUBCASE 1'
            CALL CALCQ1A             
            SCASE = 'Q1 ; SUBCASE 1'
         ELSE
            SCASE = 'Q1 ; SUBCASE 3' 
            CALL CALCMDRP (RH, DRMQ1, DRNH4NO3, CALCQ1A, CALCQ2A)
            SCASE = 'Q1 ; SUBCASE 3'
         ENDIF

      ELSE IF (.NOT.EXNO .AND. EXCL) THEN 
         IF (RH.LT.DRMG2) THEN
            SCASE = 'Q1 ; SUBCASE 1'
            CALL CALCQ1A             
            SCASE = 'Q1 ; SUBCASE 1'
         ELSE
            SCASE = 'Q1 ; SUBCASE 4' 
            CALL CALCMDRP (RH, DRMG2, DRNH4CL, CALCQ1A, CALCQ3A)
            SCASE = 'Q1 ; SUBCASE 4'
         ENDIF

      ELSE                                
         IF (RH.LT.DRMG3) THEN
            SCASE = 'Q1 ; SUBCASE 1'
            CALL CALCQ1A             
            SCASE = 'Q1 ; SUBCASE 1'
         ELSE
            SCASE = 'Q1 ; SUBCASE 5' 
            CALL CALCMDRP (RH, DRMG3, DRNH42S4, CALCQ1A, CALCQ4)
            SCASE = 'Q1 ; SUBCASE 5'
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCQ1



















      SUBROUTINE CALCQ1A
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      CNA2SO4 = 0.5d0*WAER(1)
      FRSO4   = MAX (WAER(2)-CNA2SO4, ZERO)

      CNH42S4 = MAX (MIN(FRSO4,0.5d0*WAER(3)), TINY)
      FRNH3   = MAX (WAER(3)-2.D0*CNH42S4, ZERO)

      CNH4NO3 = MIN (FRNH3, WAER(4))

      FRNH3   = MAX (FRNH3-CNH4NO3, ZERO)

      CNH4CL  = MIN (FRNH3, WAER(5))

      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)



      WATER   = ZERO

      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      RETURN



      END SUBROUTINE CALCQ1A
















      SUBROUTINE CALCR6
      USE SOLUT
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ






      CALL CALCR1A

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      FRST   = .TRUE.
      CALAIN = .TRUE.
      CALAOU = .TRUE.



      CALL CALCMR



      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP
      AKW = XKW*RH*WATER*WATER                        

      NAI    = WAER(1)
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)



      GG  = 2.D0*WAER(2) + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
      A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
      A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

      GNH3     = NH4I/HI/A2
      GHNO3    = HI*NO3I/A3
      GHCL     = HI*CLI /A4

      GASAQ(1) = NH3AQ
      GASAQ(2) = CLAQ
      GASAQ(3) = NO3AQ

      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4CL   = ZERO
      CNACL    = ZERO
      CNANO3   = ZERO
      CNA2SO4  = ZERO

      RETURN



      END SUBROUTINE CALCR6
















      SUBROUTINE CALCR5
      USE SOLUT
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      LOGICAL PSCONV
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ




      LOGICAL  NEAN, NEAC, NESN, NESC



      CALL CALCR1A                             

      NEAN = CNH4NO3.LE.TINY    
      NEAC = CNH4CL .LE.TINY    
      NESN = CNANO3 .LE.TINY    
      NESC = CNACL  .LE.TINY    
      IF (NEAN .AND. NEAC .AND. NESN .AND. NESC) RETURN

      CHI1   = CNA2SO4

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      PSIO   =-GREAT



      CALL CALCMR

      FRST   = .TRUE.
      CALAIN = .TRUE.
      CALAOU = .TRUE.
      PSCONV = .FALSE.



      NAI    = WAER(1)
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP
      A5  = XK5*(WATER/GAMA(2))**3.                   
      AKW = XKW*RH*WATER*WATER                        



      ROOT = ZERO
      IF (NAI*NAI*SO4I .GT. A5) THEN
         BB =-3.D0*CHI1
         CC = 3.D0*CHI1**2.0
         DD =-CHI1**3.0 + 0.25D0*A5
         CALL POLY3(BB, CC, DD, ROOT, ISLV)
         IF (ISLV.NE.0) ROOT = TINY
         ROOT = MIN (MAX(ROOT,ZERO), CHI1)
         PSI1 = CHI1-ROOT
      ENDIF
      PSCONV = ABS(PSI1-PSIO) .LE. EPS*PSIO
      PSIO   = PSI1



      NAI  = WAER(1) - 2.D0*ROOT
      SO4I = WAER(2) - ROOT
      NH4I = WAER(3)
      NO3I = WAER(4)
      CLI  = WAER(5)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         IF (PSCONV) GOTO 20
      ENDIF
10    CONTINUE




20    A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 

      A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
      A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

      GNH3     = NH4I/HI/A2  
      GHNO3    = HI*NO3I/A3
      GHCL     = HI*CLI /A4

      GASAQ(1) = NH3AQ
      GASAQ(2) = CLAQ
      GASAQ(3) = NO3AQ

      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4CL   = ZERO
      CNACL    = ZERO
      CNANO3   = ZERO
      CNA2SO4  = CHI1 - PSI1

      RETURN



      END SUBROUTINE CALCR5

















      SUBROUTINE CALCR4
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCR1A, CALCR5



      SCASE = 'R4 ; SUBCASE 2'
      CALL CALCR1A              
      SCASE = 'R4 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN .OR. EXSN .OR. EXSC) THEN   
         IF (RH.GE.DRMH1) THEN
            SCASE = 'R4 ; SUBCASE 1'
            CALL CALCR4A
            SCASE = 'R4 ; SUBCASE 1'
         ENDIF

      ELSE IF (EXAC) THEN                  
         IF (RH.GE.DRMR5) THEN
            SCASE = 'R4 ; SUBCASE 3'
            CALL CALCMDRP (RH, DRMR5, DRNH4CL, CALCR1A, CALCR5)
            SCASE = 'R4 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCR4



















      SUBROUTINE CALCR4A
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      LOGICAL PSCONV1, PSCONV4
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ






      FRST    = .TRUE.
      CALAIN  = .TRUE.
      CALAOU  = .TRUE.
      PSCONV1 = .FALSE.
      PSCONV4 = .FALSE.
      PSIO1   =-GREAT
      PSIO4   =-GREAT



      CALL CALCR1A

      CHI1   = CNA2SO4      
      CHI4   = CNH4CL

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      CALL CALCMR           

      NAI    = WAER(1)      
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP
      A5  = XK5 *(WATER/GAMA(2))**3.                  
      A14 = XK14*(WATER/GAMA(6))**2.                  
      AKW = XKW*RH*WATER*WATER                        



      ROOT = ZERO
      IF (NAI*NAI*SO4I .GT. A5) THEN
         BB =-3.D0*CHI1
         CC = 3.D0*CHI1**2.0
         DD =-CHI1**3.0 + 0.25D0*A5
         CALL POLY3(BB, CC, DD, ROOT, ISLV)
         IF (ISLV.NE.0) ROOT = TINY
         ROOT = MIN (MAX(ROOT,ZERO), CHI1)
         PSI1 = CHI1-ROOT
         NAI  = WAER(1) - 2.D0*ROOT
         SO4I = WAER(2) - ROOT
      ENDIF
      PSCONV1 = ABS(PSI1-PSIO1) .LE. EPS*PSIO1
      PSIO1   = PSI1



      ROOT = ZERO
      IF (NH4I*CLI .GT. A14) THEN
         BB   =-(NH4I + CLI)
         CC   =-A14 + NH4I*CLI
         DD   = BB*BB - 4.D0*CC
         ROOT = 0.5D0*(-BB-SQRT(DD))
         IF (ROOT.GT.TINY) THEN
            ROOT    = MIN(ROOT, CHI4)
            PSI4    = CHI4 - ROOT
            NH4I    = WAER(3) - ROOT
            CLI     = WAER(5) - ROOT
         ENDIF
      ENDIF
      PSCONV4 = ABS(PSI4-PSIO4) .LE. EPS*PSIO4
      PSIO4   = PSI4

      NO3I   = WAER(4)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         IF (PSCONV1 .AND. PSCONV4) GOTO 20
      ENDIF
10    CONTINUE




20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4

      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ

      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = CHI4 - PSI4
      CNACL   = ZERO
      CNANO3  = ZERO
      CNA2SO4 = CHI1 - PSI1

      RETURN



      END SUBROUTINE CALCR4A

















      SUBROUTINE CALCR3
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCR1A, CALCR4A, CALCR5



      SCASE = 'R3 ; SUBCASE 2'
      CALL CALCR1A              
      SCASE = 'R3 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN .OR. EXSN) THEN                   
         IF (RH.GE.DRMH1) THEN
            SCASE = 'R3 ; SUBCASE 1'
            CALL CALCR3A
            SCASE = 'R3 ; SUBCASE 1'
         ENDIF

      ELSE IF (.NOT.EXAN .AND. .NOT.EXSN) THEN   
         IF      (     EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN
               SCASE = 'R3 ; SUBCASE 3'
               CALL CALCMDRP (RH, DRMR4, DRNACL, CALCR1A, CALCR4A)
               SCASE = 'R3 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN
               SCASE = 'R3 ; SUBCASE 4'
               CALL CALCMDRP (RH, DRMR2, DRNACL, CALCR1A, CALCR4A)
               SCASE = 'R3 ; SUBCASE 4'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN
               SCASE = 'R3 ; SUBCASE 5'
               CALL CALCMDRP (RH, DRMR5, DRNACL, CALCR1A, CALCR5)
               SCASE = 'R3 ; SUBCASE 5'
            ENDIF
         ENDIF

      ENDIF

      RETURN



      END SUBROUTINE CALCR3


















      SUBROUTINE CALCR3A
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      LOGICAL PSCONV1, PSCONV3, PSCONV4
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ






      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.
      PSCONV1 =.TRUE.
      PSCONV3 =.TRUE.
      PSCONV4 =.TRUE.
      PSI1O   =-GREAT
      PSI3O   =-GREAT
      PSI4O   =-GREAT
      ROOT1   = ZERO
      ROOT2   = ZERO
      ROOT3   = ZERO



      CALL CALCR1A

      CHI1   = CNA2SO4      
      CHI4   = CNH4CL
      CHI3   = CNACL

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      CALL CALCMR           

      NAI    = WAER(1)      
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO

      MOLAL(1) = ZERO
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I

      CALL CALCACT          



      DO 10 I=1,NSWEEP
      A5  = XK5 *(WATER/GAMA(2))**3.                  
      A8  = XK8 *(WATER/GAMA(1))**2.                  
      A14 = XK14*(WATER/GAMA(6))**2.                  
      AKW = XKW*RH*WATER*WATER                        



      IF (NH4I*CLI .GT. A14) THEN
         BB    =-(WAER(3) + WAER(5) - ROOT3)
         CC    =-A14 + NH4I*(WAER(5) - ROOT3)
         DD    = MAX(BB*BB - 4.D0*CC, ZERO)
         ROOT2A= 0.5D0*(-BB+SQRT(DD))
         ROOT2B= 0.5D0*(-BB-SQRT(DD))
         IF (ZERO.LE.ROOT2A) THEN
            ROOT2 = ROOT2A
         ELSE
            ROOT2 = ROOT2B
         ENDIF
         ROOT2 = MIN(MAX(ZERO, ROOT2), MAX(WAER(5)-ROOT3,ZERO),   &
                     CHI4, WAER(3))
         PSI4  = CHI4 - ROOT2
      ENDIF
      PSCONV4 = ABS(PSI4-PSI4O) .LE. EPS*PSI4O
      PSI4O   = PSI4



      IF (NAI*NAI*SO4I .GT. A5) THEN
         BB =-(CHI1 + WAER(1) - ROOT3)
         CC = 0.25D0*(WAER(1) - ROOT3)*(4.D0*CHI1+WAER(1)-ROOT3)
         DD =-0.25D0*(CHI1*(WAER(1)-ROOT3)**2.D0 - A5)
         CALL POLY3(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MIN (MAX(ROOT1,ZERO), MAX(WAER(1)-ROOT3,ZERO),   &
                      CHI1, WAER(2))
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      NAI = WAER(1) - (2.D0*ROOT1 + ROOT3)
      SO4I= WAER(2) - ROOT1
      NH4I= WAER(3) - ROOT2
      CLI = WAER(5) - (ROOT3 + ROOT2)
      NO3I= WAER(4)



      IF (NAI*CLI .GT. A8) THEN
         BB    =-((CHI1-2.D0*ROOT1) + (WAER(5) - ROOT2))
         CC    = (CHI1-2.D0*ROOT1)*(WAER(5) - ROOT2) - A8
         DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT3A= 0.5D0*(-BB-SQRT(DD))
         ROOT3B= 0.5D0*(-BB+SQRT(DD))
         IF (ZERO.LE.ROOT3A) THEN
            ROOT3 = ROOT3A
         ELSE
            ROOT3 = ROOT3B
         ENDIF
         ROOT3   = MIN(MAX(ROOT3, ZERO), CHI3)
         PSI3    = CHI3-ROOT3
      ENDIF
      PSCONV3 = ABS(PSI3-PSI3O) .LE. EPS*PSI3O
      PSI3O   = PSI3



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         IF (PSCONV1.AND.PSCONV3.AND.PSCONV4) GOTO 20
      ENDIF
10    CONTINUE




20    IF (CLI.LE.TINY .AND. WAER(5).GT.TINY) THEN 
         DO 30 I=1,NIONS
            MOLAL(I) = ZERO
30       CONTINUE
         DO 40 I=1,NGASAQ
            GASAQ(I) = ZERO
40       CONTINUE
         CALL CALCR1A
      ELSE
         A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
         A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
         A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

         GNH3    = NH4I/HI/A2
         GHNO3   = HI*NO3I/A3
         GHCL    = HI*CLI /A4

         GASAQ(1)= NH3AQ
         GASAQ(2)= CLAQ
         GASAQ(3)= NO3AQ

         CNH42S4 = ZERO
         CNH4NO3 = ZERO
         CNH4CL  = CHI4 - PSI4
         CNACL   = CHI3 - PSI3
         CNANO3  = ZERO
         CNA2SO4 = CHI1 - PSI1
      ENDIF

      RETURN



      END SUBROUTINE CALCR3A

















      SUBROUTINE CALCR2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCR1A, CALCR3A, CALCR4A, CALCR5



      SCASE = 'R2 ; SUBCASE 2'
      CALL CALCR1A              
      SCASE = 'R2 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN) THEN                             
         IF (RH.GE.DRMH1) THEN
            SCASE = 'R2 ; SUBCASE 1'
            CALL CALCR2A
            SCASE = 'R2 ; SUBCASE 1'
         ENDIF

      ELSE IF (.NOT.EXAN) THEN                   
         IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMH2) THEN
               SCASE = 'R2 ; SUBCASE 3'
               CALL CALCMDRP (RH, DRMH2, DRNANO3, CALCR1A, CALCR3A)
               SCASE = 'R2 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR1) THEN
               SCASE = 'R2 ; SUBCASE 4'
               CALL CALCMDRP (RH, DRMR1, DRNANO3, CALCR1A, CALCR3A)
               SCASE = 'R2 ; SUBCASE 4'
            ENDIF

         ELSE IF (.NOT.EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN
               SCASE = 'R2 ; SUBCASE 5'
               CALL CALCMDRP (RH, DRMR2, DRNACL, CALCR1A, CALCR4A)
               SCASE = 'R2 ; SUBCASE 5'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR3) THEN
               SCASE = 'R2 ; SUBCASE 6'
               CALL CALCMDRP (RH, DRMR3, DRNANO3, CALCR1A, CALCR3A)
               SCASE = 'R2 ; SUBCASE 6'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN
               SCASE = 'R2 ; SUBCASE 7'
               CALL CALCMDRP (RH, DRMR4, DRNACL, CALCR1A, CALCR4A)
               SCASE = 'R2 ; SUBCASE 7'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN
               SCASE = 'R2 ; SUBCASE 8'
               CALL CALCMDRP (RH, DRMR5, DRNH4CL, CALCR1A, CALCR5)
               SCASE = 'R2 ; SUBCASE 8'
            ENDIF

         ELSE IF (     EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR6) THEN
               SCASE = 'R2 ; SUBCASE 9'
               CALL CALCMDRP (RH, DRMR6, DRNANO3, CALCR1A, CALCR3A)
               SCASE = 'R2 ; SUBCASE 9'
            ENDIF
         ENDIF

      ENDIF

      RETURN



      END SUBROUTINE CALCR2


















      SUBROUTINE CALCR2A
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      LOGICAL PSCONV1, PSCONV2, PSCONV3, PSCONV4
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ






      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV1 =.TRUE.
      PSCONV2 =.TRUE.
      PSCONV3 =.TRUE.
      PSCONV4 =.TRUE.

      PSI1O   =-GREAT
      PSI2O   =-GREAT
      PSI3O   =-GREAT
      PSI4O   =-GREAT

      ROOT1   = ZERO
      ROOT2   = ZERO
      ROOT3   = ZERO
      ROOT4   = ZERO



      CALL CALCR1A

      CHI1   = CNA2SO4      
      CHI2   = CNANO3
      CHI3   = CNACL
      CHI4   = CNH4CL

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      CALL CALCMR           

      NAI    = WAER(1)      
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO

      MOLAL(1) = ZERO
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I

      CALL CALCACT          



      DO 10 I=1,NSWEEP
      A5  = XK5 *(WATER/GAMA(2))**3.                  
      A8  = XK8 *(WATER/GAMA(1))**2.                  
      A9  = XK9 *(WATER/GAMA(3))**2.                  
      A14 = XK14*(WATER/GAMA(6))**2.                  
      AKW = XKW*RH*WATER*WATER                        



      IF (NH4I*CLI .GT. A14) THEN
         BB    =-(WAER(3) + WAER(5) - ROOT3)
         CC    = NH4I*(WAER(5) - ROOT3) - A14
         DD    = MAX(BB*BB - 4.D0*CC, ZERO)
         DD    = SQRT(DD)
         ROOT2A= 0.5D0*(-BB+DD)
         ROOT2B= 0.5D0*(-BB-DD)
         IF (ZERO.LE.ROOT2A) THEN
            ROOT2 = ROOT2A
         ELSE
            ROOT2 = ROOT2B
         ENDIF
         ROOT2 = MIN(MAX(ROOT2, ZERO), CHI4)
         PSI4  = CHI4 - ROOT2
      ENDIF
      PSCONV4 = ABS(PSI4-PSI4O) .LE. EPS*PSI4O
      PSI4O   = PSI4



      IF (NAI*NAI*SO4I .GT. A5) THEN
         BB =-(WAER(2) + WAER(1) - ROOT3 - ROOT4)
         CC = WAER(1)*(2.D0*ROOT3 + 2.D0*ROOT4 - 4.D0*WAER(2) - ONE)   &
             -(ROOT3 + ROOT4)**2.0 + 4.D0*WAER(2)*(ROOT3 + ROOT4)
         CC =-0.25*CC
         DD = WAER(1)*WAER(2)*(ONE - 2.D0*ROOT3 - 2.D0*ROOT4) +   &
              WAER(2)*(ROOT3 + ROOT4)**2.0 - A5
         DD =-0.25*DD
         CALL POLY3(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MIN (MAX(ROOT1,ZERO), CHI1)
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NAI*NO3I .GT. A9) THEN
         BB    =-(WAER(4) + WAER(1) - 2.D0*ROOT1 - ROOT3)
         CC    = WAER(4)*(WAER(1) - 2.D0*ROOT1 - ROOT3) - A9
         DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT4A= 0.5D0*(-BB-DD)
         ROOT4B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT4A) THEN
            ROOT4 = ROOT4A
         ELSE
            ROOT4 = ROOT4B
         ENDIF
         ROOT4 = MIN(MAX(ROOT4, ZERO), CHI2)
         PSI2  = CHI2-ROOT4
      ENDIF
      PSCONV2 = ABS(PSI2-PSI2O) .LE. EPS*PSI2O
      PSI2O   = PSI2



      NAI = WAER(1) - (2.D0*ROOT1 + ROOT3 + ROOT4)
      SO4I= WAER(2) - ROOT1
      NH4I= WAER(3) - ROOT2
      NO3I= WAER(4) - ROOT4
      CLI = WAER(5) - (ROOT3 + ROOT2)



      IF (NAI*CLI .GT. A8) THEN
         BB    =-(WAER(1) - 2.D0*ROOT1 + WAER(5) - ROOT2 - ROOT4)
         CC    = (WAER(5) + ROOT2)*(WAER(1) - 2.D0*ROOT1 - ROOT4) - A8
         DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT3A= 0.5D0*(-BB-DD)
         ROOT3B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT3A) THEN
            ROOT3 = ROOT3A
         ELSE
            ROOT3 = ROOT3B
         ENDIF
         ROOT3   = MIN(MAX(ROOT3, ZERO), CHI3)
         PSI3    = CHI3-ROOT3
      ENDIF
      PSCONV3 = ABS(PSI3-PSI3O) .LE. EPS*PSI3O
      PSI3O   = PSI3



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF



      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         IF (PSCONV1.AND.PSCONV2.AND.PSCONV3.AND.PSCONV4) GOTO 20
      ENDIF
10    CONTINUE




20    IF (CLI.LE.TINY .AND. WAER(5).GT.TINY) THEN 
         DO 30 I=1,NIONS
            MOLAL(I) = ZERO
30       CONTINUE
         DO 40 I=1,NGASAQ
            GASAQ(I) = ZERO
40       CONTINUE
         CALL CALCR1A
      ELSE                                     
         A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
         A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
         A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

         GNH3    = NH4I/HI/A2
         GHNO3   = HI*NO3I/A3
         GHCL    = HI*CLI /A4

         GASAQ(1)= NH3AQ
         GASAQ(2)= CLAQ
         GASAQ(3)= NO3AQ

         CNH42S4 = ZERO
         CNH4NO3 = ZERO
         CNH4CL  = CHI4 - PSI4
         CNACL   = CHI3 - PSI3
         CNANO3  = CHI2 - PSI2
         CNA2SO4 = CHI1 - PSI1
      ENDIF

      RETURN



      END SUBROUTINE CALCR2A

















      SUBROUTINE CALCR1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCR1A, CALCR2A, CALCR3A, CALCR4A, CALCR5



      SCASE = 'R1 ; SUBCASE 1'
      CALL CALCR1A              
      SCASE = 'R1 ; SUBCASE 1'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN.AND.EXAC.AND.EXSC.AND.EXSN) THEN  
         IF (RH.GE.DRMH1) THEN
            SCASE = 'R1 ; SUBCASE 2'  
            CALL CALCMDRP (RH, DRMH1, DRNH4NO3, CALCR1A, CALCR2A)
            SCASE = 'R1 ; SUBCASE 2'
         ENDIF

      ELSE IF (.NOT.EXAN) THEN                   
         IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMH2) THEN
               SCASE = 'R1 ; SUBCASE 3'
               CALL CALCMDRP (RH, DRMH2, DRNANO3, CALCR1A, CALCR3A)
               SCASE = 'R1 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR1) THEN
               SCASE = 'R1 ; SUBCASE 4'
               CALL CALCMDRP (RH, DRMR1, DRNANO3, CALCR1A, CALCR3A)
               SCASE = 'R1 ; SUBCASE 4'
            ENDIF

         ELSE IF (.NOT.EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN
               SCASE = 'R1 ; SUBCASE 5'
               CALL CALCMDRP (RH, DRMR2, DRNACL, CALCR1A, CALCR3A) 
               SCASE = 'R1 ; SUBCASE 5'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR3) THEN
               SCASE = 'R1 ; SUBCASE 6'
               CALL CALCMDRP (RH, DRMR3, DRNANO3, CALCR1A, CALCR3A)
               SCASE = 'R1 ; SUBCASE 6'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN
               SCASE = 'R1 ; SUBCASE 7'
               CALL CALCMDRP (RH, DRMR4, DRNACL, CALCR1A, CALCR3A) 
               SCASE = 'R1 ; SUBCASE 7'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN
               SCASE = 'R1 ; SUBCASE 8'
               CALL CALCMDRP (RH, DRMR5, DRNH4CL, CALCR1A, CALCR3A) 
               SCASE = 'R1 ; SUBCASE 8'
            ENDIF

         ELSE IF (     EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR6) THEN
               SCASE = 'R1 ; SUBCASE 9'
               CALL CALCMDRP (RH, DRMR6, DRNANO3, CALCR1A, CALCR3A)
               SCASE = 'R1 ; SUBCASE 9'
            ENDIF
         ENDIF

      ELSE IF (.NOT.EXAC) THEN                   
         IF      (     EXAN .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR7) THEN
               SCASE = 'R1 ; SUBCASE 10'
               CALL CALCMDRP (RH, DRMR7, DRNH4NO3, CALCR1A, CALCR2A)
               SCASE = 'R1 ; SUBCASE 10'
            ENDIF

         ELSE IF (     EXAN .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR8) THEN
               SCASE = 'R1 ; SUBCASE 11'
               CALL CALCMDRP (RH, DRMR8, DRNH4NO3, CALCR1A, CALCR2A)
               SCASE = 'R1 ; SUBCASE 11'
            ENDIF

         ELSE IF (     EXAN .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR9) THEN
               SCASE = 'R1 ; SUBCASE 12'
               CALL CALCMDRP (RH, DRMR9, DRNH4NO3, CALCR1A, CALCR2A)
               SCASE = 'R1 ; SUBCASE 12'
            ENDIF

         ELSE IF (     EXAN .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR10) THEN
               SCASE = 'R1 ; SUBCASE 13'
               CALL CALCMDRP (RH, DRMR10, DRNH4NO3, CALCR1A, CALCR2A)
               SCASE = 'R1 ; SUBCASE 13'
            ENDIF
         ENDIF

      ELSE IF (.NOT.EXSN) THEN                  
         IF      (     EXAN .AND.      EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR11) THEN
               SCASE = 'R1 ; SUBCASE 14'
               CALL CALCMDRP (RH, DRMR11, DRNH4NO3, CALCR1A, CALCR2A)
               SCASE = 'R1 ; SUBCASE 14'
            ENDIF

         ELSE IF (     EXAN .AND.      EXAC .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR12) THEN
               SCASE = 'R1 ; SUBCASE 15'
               CALL CALCMDRP (RH, DRMR12, DRNH4NO3, CALCR1A, CALCR2A)
               SCASE = 'R1 ; SUBCASE 15'
            ENDIF
         ENDIF

      ELSE IF (.NOT.EXSC) THEN                  
         IF      (     EXAN .AND.      EXAC .AND.      EXSN) THEN
            IF (RH.GE.DRMR13) THEN
               SCASE = 'R1 ; SUBCASE 16'
               CALL CALCMDRP (RH, DRMR13, DRNH4NO3, CALCR1A, CALCR2A)
               SCASE = 'R1 ; SUBCASE 16'
            ENDIF
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCR1



















      SUBROUTINE CALCR1A
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      CNA2SO4 = WAER(2)
      FRNA    = MAX (WAER(1)-2*CNA2SO4, ZERO)

      CNH42S4 = ZERO

      CNANO3  = MIN (FRNA, WAER(4))
      FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
      FRNA    = MAX (FRNA-CNANO3, ZERO)

      CNACL   = MIN (FRNA, WAER(5))
      FRCL    = MAX (WAER(5)-CNACL, ZERO)
      FRNA    = MAX (FRNA-CNACL, ZERO)

      CNH4NO3 = MIN (FRNO3, WAER(3))
      FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
      FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)

      CNH4CL  = MIN (FRCL, FRNH3)
      FRCL    = MAX (FRCL-CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)



      WATER   = ZERO

      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      RETURN



      END SUBROUTINE CALCR1A
