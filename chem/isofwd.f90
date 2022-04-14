




















      SUBROUTINE ISRP1F (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)



      CALL INIT1 (WI, RHI, TEMPI)



      SULRAT = W(3)/W(2)





      IF (2.0.LE.SULRAT) THEN
      DC   = W(3) - 2.001D0*W(2)  
      W(3) = W(3) + MAX(-DC, ZERO)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'A2'
         CALL CALCA2                 
      ELSE

         IF (RH.LT.DRNH42S4) THEN
            SCASE = 'A1'
            CALL CALCA1              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'A2'
            CALL CALCA2              
         ENDIF
      ENDIF



      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL CALCB4                 
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'B1'
            CALL CALCB1              

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN
            SCASE = 'B2'
            CALL CALCB2              

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'B3'
            CALL CALCB3              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'B4'
            CALL CALCB4              
         ENDIF
      ENDIF
      CALL CALCNH3



      ELSEIF (SULRAT.LT.1.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC2                 
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'C1'
            CALL CALCC1              

         ELSEIF (DRNH4HS4.LE.RH) THEN
            SCASE = 'C2'
            CALL CALCC2              

         ENDIF
      ENDIF
      CALL CALCNH3
      ENDIF



      RETURN



      END SUBROUTINE ISRP1F                 
















      SUBROUTINE ISRP2F (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)



      CALL INIT2 (WI, RHI, TEMPI)



      SULRAT = W(3)/W(2)





      IF (2.0.LE.SULRAT) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'D3'
         CALL CALCD3                 
      ELSE

         IF (RH.LT.DRNH4NO3) THEN
            SCASE = 'D1'
            CALL CALCD1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'D2'
            CALL CALCD2              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'D3'
            CALL CALCD3              
         ENDIF
      ENDIF







      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL CALCB4                 
         SCASE = 'E4'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'B1'
            CALL CALCB1              
            SCASE = 'E1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN
            SCASE = 'B2'
            CALL CALCB2              
            SCASE = 'E2'

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'B3'
            CALL CALCB3              
            SCASE = 'E3'

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'B4'
            CALL CALCB4              
            SCASE = 'E4'
         ENDIF
      ENDIF

      CALL CALCNA                 







      ELSEIF (SULRAT.LT.1.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC2                 
         SCASE = 'F2'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'C1'
            CALL CALCC1              
            SCASE = 'F1'

         ELSEIF (DRNH4HS4.LE.RH) THEN
            SCASE = 'C2'
            CALL CALCC2              
            SCASE = 'F2'
         ENDIF
      ENDIF

      CALL CALCNA                 
      ENDIF



      RETURN



      END SUBROUTINE ISRP2F                 
















      SUBROUTINE ISRP3F (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)



      WI(3) = MAX (WI(3), 1.D-10)  
      WI(5) = MAX (WI(5), 1.D-10)  



      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
         WI(1) = 1.D-10  
         WI(2) = 1.D-10  
      ENDIF



      CALL ISOINIT3 (WI, RHI, TEMPI)



      REST = 2.D0*W(2) + W(4) + W(5)
      IF (W(1).GT.REST) THEN            
         W(1) = (ONE-1D-6)*REST         
         CALL PUSHERR (0050, 'ISRP3F')  
      ENDIF



      SULRAT = (W(1)+W(3))/W(2)
      SODRAT = W(1)/W(2)





      IF (2.0.LE.SULRAT .AND. SODRAT.LT.2.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'G5'
         CALL CALCG5                 
      ELSE

         IF (RH.LT.DRNH4NO3) THEN
            SCASE = 'G1'
            CALL CALCG1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH4CL) THEN
            SCASE = 'G2'
            CALL CALCG2              

         ELSEIF (DRNH4CL.LE.RH  .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'G3'
            CALL CALCG3              

        ELSEIF (DRNH42S4.LE.RH  .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'G4'
            CALL CALCG4              

         ELSEIF (DRNA2SO4.LE.RH) THEN
            SCASE = 'G5'
            CALL CALCG5              
         ENDIF
      ENDIF



      ELSE IF (SULRAT.GE.2.0 .AND. SODRAT.GE.2.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'H6'
         CALL CALCH6                 
      ELSE

         IF (RH.LT.DRNH4NO3) THEN
            SCASE = 'H1'
            CALL CALCH1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN
            SCASE = 'H2'
            CALL CALCH2              

         ELSEIF (DRNANO3.LE.RH  .AND. RH.LT.DRNACL) THEN
            SCASE = 'H3'
            CALL CALCH3              

         ELSEIF (DRNACL.LE.RH   .AND. RH.LT.DRNH4Cl) THEN
            SCASE = 'H4'
            CALL CALCH4              

         ELSEIF (DRNH4Cl.LE.RH .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'H5'
            CALL CALCH5              

         ELSEIF (DRNA2SO4.LE.RH) THEN
            SCASE = 'H6'
            CALL CALCH6              
         ENDIF
      ENDIF



      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'I6'
         CALL CALCI6                 
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'I1'
            CALL CALCI1              

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN
            SCASE = 'I2'
            CALL CALCI2              

         ELSEIF (DRNAHSO4.LE.RH .AND. RH.LT.DRLC) THEN
            SCASE = 'I3'
            CALL CALCI3              

         ELSEIF (DRLC.LE.RH     .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'I4'
            CALL CALCI4              

         ELSEIF (DRNH42S4.LE.RH .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'I5'
            CALL CALCI5              

         ELSEIF (DRNA2SO4.LE.RH) THEN
            SCASE = 'I6'
            CALL CALCI6              
         ENDIF
      ENDIF

      CALL CALCNHA                
      CALL CALCNH3                



      ELSEIF (SULRAT.LT.1.0) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'J3'
         CALL CALCJ3                 
      ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'J1'
            CALL CALCJ1              

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN
            SCASE = 'J2'
            CALL CALCJ2              

         ELSEIF (DRNAHSO4.LE.RH) THEN
            SCASE = 'J3'
            CALL CALCJ3
         ENDIF
      ENDIF

      CALL CALCNHA                
      CALL CALCNH3                
      ENDIF



      RETURN



      END SUBROUTINE ISRP3F                 























      SUBROUTINE CALCA2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      CALAOU    =.TRUE.       
      OMELO     = TINY        
      OMEHI     = 2.0D0*W(2)  



      MOLAL(5) = W(2)
      MOLAL(6) = ZERO
      CALL CALCMR



      X1 = OMEHI
      Y1 = FUNCA2 (X1)
      IF (ABS(Y1).LE.EPS) RETURN



      DX = (OMEHI-OMELO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, OMELO)
         Y2 = FUNCA2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE
      IF (ABS(Y2).LE.EPS) THEN
         RETURN
      ELSE
         CALL PUSHERR (0001, 'CALCA2')    
         RETURN
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCA2 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCA2')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCA2 (X3)
      RETURN



      END SUBROUTINE CALCA2













      DOUBLE PRECISION FUNCTION FUNCA2 (OMEGI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION LAMDA



      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI    = W(2)         



      DO 10 I=1,NSWEEP
         A1    = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         A2    = XK2*R*TEMP/XKW*(GAMA(8)/GAMA(9))**2.
         A3    = XKW*RH*WATER*WATER

         LAMDA = PSI/(A1/OMEGI+ONE)
         ZETA  = A3/OMEGI



         MOLAL (1) = OMEGI                        
         MOLAL (3) = W(3)/(ONE/A2/OMEGI + ONE)    
         MOLAL (5) = MAX(PSI-LAMDA,TINY)          
         MOLAL (6) = LAMDA                        
         GNH3      = MAX (W(3)-MOLAL(3), TINY)    
         COH       = ZETA                         



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    DENOM = (2.0*MOLAL(5)+MOLAL(6))
      FUNCA2= (MOLAL(3)/DENOM - ONE) + MOLAL(1)/DENOM
      RETURN



      END FUNCTION FUNCA2        






















      SUBROUTINE CALCA1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      CNH42S4 = W(2)
      GNH3    = MAX (W(3)-2.0*CNH42S4, ZERO)
      RETURN



      END SUBROUTINE CALCA1
























      SUBROUTINE CALCB4
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.



      CALL CALCB1A         
      MOLALR(13) = CLC
      MOLALR(9)  = CNH4HS4
      MOLALR(4)  = CNH42S4
      CLC        = ZERO
      CNH4HS4    = ZERO
      CNH42S4    = ZERO
      WATER      = MOLALR(13)/M0(13)+MOLALR(9)/M0(9)+MOLALR(4)/M0(4)

      MOLAL(3)   = W(3)   

      DO 20 I=1,NSWEEP
         AK1   = XK1*((GAMA(8)/GAMA(7))**2.)*(WATER/GAMA(7))
         BET   = W(2)
         GAM   = MOLAL(3)

         BB    = BET + AK1 - GAM
         CC    =-AK1*BET
         DD    = BB*BB - 4.D0*CC



         MOLAL (5) = MAX(TINY,MIN(0.5*(-BB + SQRT(DD)), W(2))) 
         MOLAL (6) = MAX(TINY,MIN(W(2)-MOLAL(5),W(2)))         
         MOLAL (1) = MAX(TINY,MIN(AK1*MOLAL(6)/MOLAL(5),W(2))) 
         CALL CALCMR                                           



         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT
20    CONTINUE

30    RETURN



      END SUBROUTINE CALCB4


















      SUBROUTINE CALCB3
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      X = MAX(2*W(2)-W(3), ZERO)   
      Y = MAX(W(3)  -W(2), ZERO)   



      IF (X.LT.Y) THEN             
         SCASE   = 'B3 ; SUBCASE 1'
         TLC     = X
         TNH42S4 = Y-X
         CALL CALCB3A (TLC,TNH42S4)      
      ELSE
         SCASE   = 'B3 ; SUBCASE 2'
         TLC     = Y
         TNH4HS4 = X-Y
         CALL CALCB3B (TLC,TNH4HS4)      
      ENDIF

      RETURN



      END SUBROUTINE CALCB3



























      SUBROUTINE CALCB3A (TLC, TNH42S4)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      CALAOU = .TRUE.         
      ZLO    = ZERO           
      ZHI    = TNH42S4        



      Z1 = ZLO
      Y1 = FUNCB3A (Z1, TLC, TNH42S4)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1



      DZ = (ZHI-ZLO)/REAL(NDIV)
      DO 10 I=1,NDIV
         Z2 = Z1+DZ
         Y2 = FUNCB3A (Z2, TLC, TNH42S4)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         Z1 = Z2
         Y1 = Y2
10    CONTINUE



      YHI= Y1                      
      IF (ABS(Y2) .LT. EPS) THEN   
         RETURN



      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         Z1 = ZHI
         Z2 = ZHI
         GOTO 40



      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Z1 = ZLO
         Z2 = ZLO
         GOTO 40
      ELSE
         CALL PUSHERR (0001, 'CALCB3A')    
         RETURN
      ENDIF



20    DO 30 I=1,MAXIT
         Z3 = 0.5*(Z1+Z2)
         Y3 = FUNCB3A (Z3, TLC, TNH42S4)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            Z2    = Z3
         ELSE
            Y1    = Y3
            Z1    = Z3
         ENDIF
         IF (ABS(Z2-Z1) .LE. EPS*Z1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCB3A')    



40    ZK = 0.5*(Z1+Z2)
      Y3 = FUNCB3A (ZK, TLC, TNH42S4)

      RETURN



      END SUBROUTINE CALCB3A               













      DOUBLE PRECISION FUNCTION FUNCB3A (ZK, Y, X)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION KK



      FRST   = .TRUE.
      CALAIN = .TRUE.
      DO 20 I=1,NSWEEP
         GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         DD    = SQRT( (ZK+GRAT1+Y)**2. + 4.0*Y*GRAT1)
         KK    = 0.5*(-(ZK+GRAT1+Y) + DD )



         MOLAL (1) = KK                
         MOLAL (5) = KK+ZK+Y           
         MOLAL (6) = MAX (Y-KK, TINY)  
         MOLAL (3) = 3.0*Y+2*ZK        
         CNH42S4   = X-ZK              
         CALL CALCMR                   



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE




30    FUNCB3A= MOLAL(5)*MOLAL(3)**2.0
      FUNCB3A= FUNCB3A/(XK7*(WATER/GAMA(4))**3.0) - ONE
      RETURN



      END FUNCTION FUNCB3A           






















      SUBROUTINE CALCB3B (Y, X)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION KK

      CALAOU = .FALSE.        
      FRST   = .FALSE.
      CALAIN = .TRUE.



      DO 20 I=1,NSWEEP
         GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         DD    = SQRT( (GRAT1+Y)**2. + 4.0*(X+Y)*GRAT1)
         KK    = 0.5*(-(GRAT1+Y) + DD )



         MOLAL (1) = KK                   
         MOLAL (5) = Y+KK                 
         MOLAL (6) = MAX (X+Y-KK, TINY)   
         MOLAL (3) = 3.0*Y+X              
         CALL CALCMR                      



         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT
20    CONTINUE

30    RETURN



      END SUBROUTINE CALCB3B       






















      SUBROUTINE CALCB2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      X = MAX(2*W(2)-W(3), TINY)   
      Y = MAX(W(3)  -W(2), TINY)   



      IF (X.LE.Y) THEN             
         SCASE = 'B2 ; SUBCASE 1'
         CALL CALCB2A (X,Y-X)      
      ELSE
         SCASE = 'B2 ; SUBCASE 2'
         CALL CALCB2B (Y,X-Y)      
      ENDIF

      RETURN



      END SUBROUTINE CALCB2





























      SUBROUTINE CALCB2A (TLC, TNH42S4)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (RH.LT.DRMLCAS) THEN
         SCASE   = 'B2 ; SUBCASE A1'    
         CLC     = TLC
         CNH42S4 = TNH42S4
         SCASE   = 'B2 ; SUBCASE A1'
      ELSE
         SCASE = 'B2 ; SUBCASE A2'
         CALL CALCB2A2 (TLC, TNH42S4)   
         SCASE = 'B2 ; SUBCASE A2'
      ENDIF

      RETURN



      END SUBROUTINE CALCB2A               


























      SUBROUTINE CALCB2A2 (TLC, TNH42S4)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (WFTYP.EQ.0) THEN
         WF = ZERO
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (DRLC-RH)/(DRLC-DRMLCAS)
      ENDIF
      ONEMWF  = ONE - WF



      CLCO     = TLC                     
      CNH42SO  = TNH42S4



      CLC     = ZERO
      CNH42S4 = ZERO
      CALL CALCB3                        



      MOLAL(1)= ONEMWF*MOLAL(1)                                   
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + 3.D0*(CLCO-CLC)) 
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               
      MOLAL(6)= ONEMWF*(CLCO-CLC)                                 

      WATER   = ONEMWF*WATER

      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4

      RETURN



      END SUBROUTINE CALCB2A2               




























      SUBROUTINE CALCB2B (TLC,TNH4HS4)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      CALAOU = .TRUE.       
      ZLO    = ZERO
      ZHI    = TLC          



      X1 = ZHI
      Y1 = FUNCB2B (X1,TNH4HS4,TLC)
      IF (ABS(Y1).LE.EPS) RETURN
      YHI= Y1                        



      DX = (ZHI-ZLO)/NDIV
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCB2B (X2,TNH4HS4,TLC)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (ABS(Y2) .LT. EPS) THEN   
         RETURN



      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         X1 = ZHI
         X2 = ZHI
         GOTO 40



      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         X1 = ZLO
         X2 = ZLO
         GOTO 40
      ELSE
         CALL PUSHERR (0001, 'CALCB2B')    
         RETURN
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCB2B (X3,TNH4HS4,TLC)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCB2B')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCB2B (X3,TNH4HS4,TLC)

      RETURN



      END SUBROUTINE CALCB2B              













      DOUBLE PRECISION FUNCTION FUNCB2B (X,TNH4HS4,TLC)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      FRST   = .TRUE.
      CALAIN = .TRUE.
      DO 20 I=1,NSWEEP
         GRAT2 = XK1*WATER*(GAMA(8)/GAMA(7))**2./GAMA(7)
         PARM  = X+GRAT2
         DELTA = PARM*PARM + 4.0*(X+TNH4HS4)*GRAT2 
         OMEGA = 0.5*(-PARM + SQRT(DELTA))         



         MOLAL (1) = OMEGA                         
         MOLAL (3) = 3.0*X+TNH4HS4                 
         MOLAL (5) = X+OMEGA                       
         MOLAL (6) = MAX (X+TNH4HS4-OMEGA, TINY)   
         CLC       = MAX(TLC-X,ZERO)               
         CALL CALCMR                               



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE




30    FUNCB2B= (MOLAL(3)**3.)*MOLAL(5)*MOLAL(6)
      FUNCB2B= FUNCB2B/(XK13*(WATER/GAMA(13))**5.) - ONE
      RETURN



      END FUNCTION FUNCB2B                
























      SUBROUTINE CALCB1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (RH.LT.DRMLCAB) THEN
         SCASE = 'B1 ; SUBCASE 1'
         CALL CALCB1A              
         SCASE = 'B1 ; SUBCASE 1'
      ELSE
         SCASE = 'B1 ; SUBCASE 2'
         CALL CALCB1B              
         SCASE = 'B1 ; SUBCASE 2'
      ENDIF

      RETURN



      END SUBROUTINE CALCB1



























      SUBROUTINE CALCB1A
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      X = 2*W(2)-W(3)       
      Y = W(3)-W(2)         



      IF (X.LE.Y) THEN      
         CLC     = X        
         CNH4HS4 = ZERO
         CNH42S4 = Y-X
      ELSE
         CLC     = Y        
         CNH4HS4 = X-Y
         CNH42S4 = ZERO
      ENDIF
      RETURN



      END SUBROUTINE CALCB1A




























      SUBROUTINE CALCB1B
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (WFTYP.EQ.0) THEN
         WF = ZERO
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (DRNH4HS4-RH)/(DRNH4HS4-DRMLCAB)
      ENDIF
      ONEMWF  = ONE - WF



      CALL CALCB1A
      CLCO     = CLC               
      CNH42SO  = CNH42S4
      CNH4HSO  = CNH4HS4



      CLC     = ZERO
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CALL CALCB2                  



      MOLAL(1)= ONEMWF*MOLAL(1)                                   
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + (CNH4HSO-CNH4HS4)   &
                      + 3.D0*(CLCO-CLC))                          
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               
      MOLAL(6)= ONEMWF*(CNH4HSO-CNH4HS4 + CLCO-CLC)               

      WATER   = ONEMWF*WATER

      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4

      RETURN



      END SUBROUTINE CALCB1B



















      SUBROUTINE CALCC2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION LAMDA, KAPA

      CALAOU =.TRUE.         
      FRST   =.TRUE.
      CALAIN =.TRUE.



      LAMDA  = W(3)           
      PSI    = W(2)-W(3)      
      DO 20 I=1,NSWEEP
         PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         BB    = PSI+PARM
         CC    =-PARM*(LAMDA+PSI)
         KAPA  = 0.5*(-BB+SQRT(BB*BB-4.0*CC))



         MOLAL(1) = PSI+KAPA                               
         MOLAL(3) = LAMDA                                  
         MOLAL(5) = KAPA                                   
         MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              
         CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3), ZERO)  
         CALL CALCMR                                       



         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT
20    CONTINUE

30    RETURN



      END SUBROUTINE CALCC2





















      SUBROUTINE CALCC1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION KLO, KHI

      CALAOU = .TRUE.    
      KLO    = TINY
      KHI    = W(3)



      X1 = KLO
      Y1 = FUNCC1 (X1)
      IF (ABS(Y1).LE.EPS) GOTO 50
      YLO= Y1



      DX = (KHI-KLO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCC1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO) GOTO 20 
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YHI= Y2                 
      IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50



      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         GOTO 50



      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         X1 = KLO
         X2 = KLO
         GOTO 40
      ELSE
         CALL PUSHERR (0001, 'CALCC1')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCC1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCC1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCC1 (X3)

50    RETURN



      END SUBROUTINE CALCC1


















      DOUBLE PRECISION FUNCTION FUNCC1 (KAPA)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION KAPA, LAMDA



      FRST   = .TRUE.
      CALAIN = .TRUE.

      PSI = W(2)-W(3)
      DO 20 I=1,NSWEEP
         PAR1  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
         PAR2  = XK12*(WATER/GAMA(9))**2.0
         BB    = PSI + PAR1
         CC    =-PAR1*(PSI+KAPA)
         LAMDA = 0.5*(-BB+SQRT(BB*BB-4*CC))



         MOLAL(1) = PSI+LAMDA                    
         MOLAL(3) = KAPA                         
         MOLAL(5) = LAMDA                        
         MOLAL(6) = MAX (ZERO, PSI+KAPA-LAMDA)   
         CNH4HS4  = MAX(W(3)-MOLAL(3), ZERO)     
         CH2SO4   = MAX(PSI, ZERO)               
         CALL CALCMR                             



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE




30    FUNCC1= (MOLAL(3)*MOLAL(6)/PAR2) - ONE
      RETURN



      END FUNCTION FUNCC1       


















      SUBROUTINE CALCD3
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CALL CALCD1A



      CHI1 = CNH4NO3               
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3

      PSI1 = CNH4NO3               
      PSI2 = CHI2
      PSI3 = ZERO
      PSI4 = ZERO

      MOLAL(5) = ZERO
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL CALCMR                  

      CALAOU = .TRUE.              
      PSI4LO = TINY                
      PSI4HI = CHI4                



60    X1 = PSI4LO
      Y1 = FUNCD3 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1                 



      DX = (PSI4HI-PSI4LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCD3 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YHI= Y1                      
      IF (ABS(Y2) .LT. EPS) THEN   
         RETURN






      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = TINY 
         YY = FUNCD3(P4)
         GOTO 50






      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) 
         IF (PSI4LO.LT.-(PSI1+PSI2)) THEN
            CALL PUSHERR (0001, 'CALCD3')  
            RETURN
         ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  
            GOTO 60                        
         ENDIF
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCD3 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCD3')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCD3 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF
      RETURN



      END SUBROUTINE CALCD3













      DOUBLE PRECISION FUNCTION FUNCD3 (P4)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4



      DO 10 I=1,NSWEEP
         A2   = XK7*(WATER/GAMA(4))**3.0
         A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7   = XKW *RH*WATER*WATER

         PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4)
         PSI3 = MIN(MAX(PSI3, ZERO), CHI3)

         BB   = PSI4 - PSI3


         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF (DENM.LE.TINY) THEN       
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.0*A7/ABB 
         ENDIF
         AHI = 2.0*A7/DENM



         MOLAL (1) = AHI                             
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2         
         MOLAL (5) = PSI2                            
         MOLAL (6) = ZERO                            
         MOLAL (7) = PSI3 + PSI1                     
         CNH42S4   = CHI2 - PSI2                     
         CNH4NO3   = ZERO                            
         GHNO3     = CHI3 - PSI3                     
         GNH3      = CHI4 - PSI4                     
         CALL CALCMR                                 



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    CONTINUE

      FUNCD3= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE
      RETURN



      END FUNCTION FUNCD3     


















      SUBROUTINE CALCD2
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CALL CALCD1A



      CHI1 = CNH4NO3               
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3

      PSI1 = CNH4NO3               
      PSI2 = CNH42S4
      PSI3 = ZERO
      PSI4 = ZERO

      MOLAL(5) = ZERO
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL CALCMR                  

      CALAOU = .TRUE.              
      PSI4LO = TINY                
      PSI4HI = CHI4                



60    X1 = PSI4LO
      Y1 = FUNCD2 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1                 



      DX   = (PSI4HI-PSI4LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCD2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) THEN



             IF (Y1 .LE. Y2) GOTO 20  
         ENDIF
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YHI= Y1                      
      IF (ABS(Y2) .LT. EPS) THEN   
         RETURN






      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = TINY 
         YY = FUNCD2(P4)
         GOTO 50






      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) 
         IF (PSI4LO.LT.-(PSI1+PSI2)) THEN
            CALL PUSHERR (0001, 'CALCD2')  
            RETURN
         ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  
            GOTO 60                        
         ENDIF
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCD2 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCD2')    



40    X3 = MIN(X1,X2)   
      Y3 = FUNCD2 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF
      RETURN



      END SUBROUTINE CALCD2













      DOUBLE PRECISION FUNCTION FUNCD2 (P4)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CALL RSTGAM       
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4
      PSI2   = CHI2



      DO 10 I=1,NSWEEP
         A2  = XK7*(WATER/GAMA(4))**3.0
         A3  = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7  = XKW *RH*WATER*WATER

         IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN
            PSI14 = PSI1+PSI4
            CALL POLY3 (PSI14,0.25*PSI14**2.,-A2/4.D0, PSI2, ISLV)  
            IF (ISLV.EQ.0) THEN
                PSI2 = MIN (PSI2, CHI2)
            ELSE
                PSI2 = ZERO
            ENDIF
         ENDIF

         PSI3  = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3  = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4)


         BB   = PSI4-PSI3 



         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF (DENM.LE.TINY) THEN       
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.d0*A7/ABB 
         ENDIF
         AHI = 2.d0*A7/DENM



         MOLAL (1) = AHI                              
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2          
         MOLAL (5) = PSI2                             
         MOLAL (6) = ZERO                             
         MOLAL (7) = PSI3 + PSI1                      
         CNH42S4   = CHI2 - PSI2                      
         CNH4NO3   = ZERO                             
         GHNO3     = CHI3 - PSI3                      
         GNH3      = CHI4 - PSI4                      
         CALL CALCMR                                  



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    CONTINUE

      FUNCD2= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE
      RETURN



      END FUNCTION FUNCD2     






















      SUBROUTINE CALCD1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCD1A, CALCD2



      IF (RH.LT.DRMASAN) THEN
         SCASE = 'D1 ; SUBCASE 1'   
         CALL CALCD1A
         SCASE = 'D1 ; SUBCASE 1'
      ELSE
         SCASE = 'D1 ; SUBCASE 2'   
         CALL CALCMDRH (RH, DRMASAN, DRNH4NO3, CALCD1A, CALCD2)
         SCASE = 'D1 ; SUBCASE 2'
      ENDIF

      RETURN



      END SUBROUTINE CALCD1

























      SUBROUTINE CALCD1A
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      PARM    = XK10/(R*TEMP)/(R*TEMP)



      CNH42S4 = W(2)
      X       = MAX(ZERO, MIN(W(3)-2.0*CNH42S4, W(4)))  
      PS      = MAX(W(3) - X - 2.0*CNH42S4, ZERO)
      OM      = MAX(W(4) - X, ZERO)

      OMPS    = OM+PS
      DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              
      ZE      = MIN(X, 0.5*(-OMPS + DIAK))              



      CNH4NO3 = X  - ZE    
      GNH3    = PS + ZE    
      GHNO3   = OM + ZE    

      RETURN



      END SUBROUTINE CALCD1A


















      SUBROUTINE CALCG5
      USE ISRPIA
      USE CASEG 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









      CALAOU = .TRUE.
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)

      PSI1   = CHI1
      PSI2   = CHI2
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    

      WATER  = CHI2/M0(4) + CHI1/M0(2)



      X1 = PSI6LO
      Y1 = FUNCG5A (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCG5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG5A (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCG5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCG5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCG5A (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    
         MOLAL(5) = MOLAL(5) - DELTA                    
         MOLAL(6) = DELTA                               
      ENDIF

      RETURN



      END SUBROUTINE CALCG5






















      DOUBLE PRECISION FUNCTION FUNCG5A (X)
      USE ISRPIA
      USE CASEG 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









      PSI6   = X
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      AKK = A4*A6



      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF


      IF(W(2).GT.TINY) THEN       
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF



      MOLAL (2) = 2.0D0*PSI1                          
      MOLAL (3) = 2.0*PSI2 + PSI4                     
      MOLAL (4) = PSI6                                
      MOLAL (5) = PSI2 + PSI1                         
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5                                

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)              
      GHNO3     = MAX(CHI5 - PSI5, TINY)              
      GHCL      = MAX(CHI6 - PSI6, TINY)              

      CNH42S4   = ZERO                                
      CNH4NO3   = ZERO                                
      CNH4CL    = ZERO                                

      CALL CALCMR                                     



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCG5A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END FUNCTION FUNCG5A    



















      SUBROUTINE CALCG4
      USE ISRPIA
      USE CASEG 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









      CALAOU = .TRUE.
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)

      PSI2   = CHI2
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    

      WATER  = CHI2/M0(4) + CHI1/M0(2)



      X1 = PSI6LO
      Y1 = FUNCG4A (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2  = X1+DX
         Y2  = FUNCG4A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1  = X2
         Y1  = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG4A (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCG4A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCG4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCG4A (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END SUBROUTINE CALCG4






















      DOUBLE PRECISION FUNCTION FUNCG4A (X)
      USE ISRPIA
      USE CASEG
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)



      DOUBLE PRECISION  NAI, NH4I, NO3I






      PSI6   = X
      PSI1   = CHI1
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0



      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF


      IF(W(2).GT.TINY) THEN       
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO) 
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF



      NH4I = 2.0*PSI2 + PSI4
      CLI  = PSI6
      SO4I = PSI2 + PSI1
      NO3I = PSI5
      NAI  = 2.0D0*PSI1

      CALL CALCPH(2.d0*SO4I+NO3I+CLI-NAI-NH4I, HI, OHI)



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ELSE
         PSI1 = ZERO
      ENDIF



      MOLAL (1) = HI
      MOLAL (2) = NAI
      MOLAL (3) = NH4I
      MOLAL (4) = CLI
      MOLAL (5) = SO4I
      MOLAL (6) = ZERO
      MOLAL (7) = NO3I



      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNH4CL    = ZERO
      CNA2SO4   = MAX(CHI1-PSI1,ZERO)



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCG4A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END FUNCTION FUNCG4A    



















      SUBROUTINE CALCG3
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCG1A, CALCG4



      IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN 
         SCASE = 'G3 ; SUBCASE 1'
         CALL CALCG3A
         SCASE = 'G3 ; SUBCASE 1'
      ELSE                                      
         SCASE = 'G1 ; SUBCASE 1'
         CALL CALCG1A
         SCASE = 'G1 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMG3) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCG1A
            SCASE = 'G3 ; SUBCASE 2'
            RETURN
         ELSE
            SCASE = 'G3 ; SUBCASE 3'  
            CALL CALCMDRH (RH, DRMG3, DRNH42S4, CALCG1A, CALCG4)
            SCASE = 'G3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCG3




















      SUBROUTINE CALCG3A
      USE ISRPIA
      USE CASEG
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









      CALAOU = .TRUE.
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    

      WATER  = TINY



      X1 = PSI6LO
      Y1 = FUNCG3A (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2  = X1+DX
         Y2  = FUNCG3A (X2)

         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1  = X2
         Y1  = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG3A (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCG3A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCG3A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCG3A (X3)



50    CONTINUE



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ELSE
         PSI1 = ZERO
      ENDIF
      MOLAL(2) = 2.0D0*PSI1               
      MOLAL(5) = MOLAL(5) + PSI1          
      CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   



      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END SUBROUTINE CALCG3A






















      DOUBLE PRECISION FUNCTION FUNCG3A (X)
      USE ISRPIA
      USE CASEG
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









      PSI6   = X
      PSI2   = CHI2
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0



      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF


      IF(W(2).GT.TINY) THEN       
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)  
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN
         CALL POLY3 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (PSI20, CHI2)
      ENDIF



      MOLAL (2) = ZERO                                
      MOLAL (3) = 2.0*PSI2 + PSI4                     
      MOLAL (4) = PSI6                                
      MOLAL (5) = PSI2                                
      MOLAL (6) = ZERO                                
      MOLAL (7) = PSI5                                

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)              
      GHNO3     = MAX(CHI5 - PSI5, TINY)              
      GHCL      = MAX(CHI6 - PSI6, TINY)              

      CNH42S4   = CHI2 - PSI2                         
      CNH4NO3   = ZERO                                
      CNH4CL    = ZERO                                

      CALL CALCMR                                     



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCG3A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END FUNCTION FUNCG3A    



















      SUBROUTINE CALCG2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCG1A, CALCG3A, CALCG4



      IF (W(4).GT.TINY) THEN        
         SCASE = 'G2 ; SUBCASE 1'
         CALL CALCG2A
         SCASE = 'G2 ; SUBCASE 1'
      ELSE                          
         SCASE = 'G1 ; SUBCASE 1'
         CALL CALCG1A
         SCASE = 'G1 ; SUBCASE 1'
      ENDIF



      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMG2) THEN             
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCG1A
            SCASE = 'G2 ; SUBCASE 2'
         ELSE
            IF (W(5).GT. TINY) THEN
               SCASE = 'G2 ; SUBCASE 3'    
               CALL CALCMDRH (RH, DRMG2, DRNH4CL, CALCG1A, CALCG3A)
               SCASE = 'G2 ; SUBCASE 3'
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMG3) THEN
               SCASE = 'G2 ; SUBCASE 4'    
               CALL CALCMDRH (RH, DRMG3, DRNH42S4, CALCG1A, CALCG4)
               SCASE = 'G2 ; SUBCASE 4'
            ELSE
               WATER = TINY
               DO 20 I=1,NIONS
                  MOLAL(I) = ZERO
20             CONTINUE
               CALL CALCG1A
               SCASE = 'G2 ; SUBCASE 2'
            ENDIF
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCG2




















      SUBROUTINE CALCG2A
      USE ISRPIA
      USE CASEG
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









      CALAOU = .TRUE.
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)

      PSI6LO = TINY
      PSI6HI = CHI6-TINY

      WATER  = TINY



      X1 = PSI6LO
      Y1 = FUNCG2A (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCG2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) WATER = TINY
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCG2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCG2A')    



40    X3 = 0.5*(X1+X2)
      IF (X3.LE.TINY2) THEN   
         WATER = TINY
      ELSE
         Y3 = FUNCG2A (X3)
      ENDIF



50    CONTINUE



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ELSE
         PSI1 = ZERO
      ENDIF
      MOLAL(2) = 2.0D0*PSI1               
      MOLAL(5) = MOLAL(5) + PSI1          
      CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   



      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA     
         MOLAL(5) = MOLAL(5) - DELTA     
         MOLAL(6) = DELTA                
      ENDIF

      RETURN



      END SUBROUTINE CALCG2A






















      DOUBLE PRECISION FUNCTION FUNCG2A (X)
      USE ISRPIA
      USE CASEG
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









      PSI6   = X
      PSI2   = CHI2
      PSI3   = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0

      DENO = MAX(CHI6-PSI6-PSI3, ZERO)
      PSI5 = CHI5/((A6/A5)*(DENO/PSI6) + ONE)

      PSI4 = MIN(PSI5+PSI6,CHI4)

      IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN
         CALL POLY3 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (PSI20, CHI2)
      ENDIF



      MOLAL (2) = ZERO                             
      MOLAL (3) = 2.0*PSI2 + PSI4                  
      MOLAL (4) = PSI6                             
      MOLAL (5) = PSI2                             
      MOLAL (6) = ZERO                             
      MOLAL (7) = PSI5                             


      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = MAX(CHI2 - PSI2, ZERO)
      CNH4NO3   = ZERO



      A3   = XK6 /(R*TEMP*R*TEMP)
      IF (GNH3*GHCL.GT.A3) THEN
         DELT = MIN(GNH3, GHCL)
         BB = -(GNH3+GHCL)
         CC = GNH3*GHCL-A3
         DD = BB*BB - 4.D0*CC
         PSI31 = 0.5D0*(-BB + SQRT(DD))
         PSI32 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
            PSI3 = PSI31
         ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
            PSI3 = PSI32
         ELSE
            PSI3 = ZERO
         ENDIF
      ELSE
         PSI3 = ZERO
      ENDIF



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3



      CALL CALCMR



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    IF (CHI4.LE.TINY) THEN
         FUNCG2A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
      ELSE
         FUNCG2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
      ENDIF

      RETURN



      END FUNCTION FUNCG2A    























      SUBROUTINE CALCG1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCG1A, CALCG2A



      IF (RH.LT.DRMG1) THEN
         SCASE = 'G1 ; SUBCASE 1'
         CALL CALCG1A              
         SCASE = 'G1 ; SUBCASE 1'
      ELSE
         SCASE = 'G1 ; SUBCASE 2'  
         CALL CALCMDRH (RH, DRMG1, DRNH4NO3, CALCG1A, CALCG2A)
         SCASE = 'G1 ; SUBCASE 2'
      ENDIF

      RETURN



      END SUBROUTINE CALCG1

























      SUBROUTINE CALCG1A
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2



      CNA2SO4 = 0.5*W(1)
      CNH42S4 = W(2) - CNA2SO4



      ALF     = W(3) - 2.0*CNH42S4
      BET     = W(5)
      GAM     = W(4)

      RTSQ    = R*TEMP*R*TEMP
      A1      = XK6/RTSQ
      A2      = XK10/RTSQ

      THETA1  = GAM - BET*(A2/A1)
      THETA2  = A2/A1



      BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
      CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
      DD      = BB*BB - 4.0D0*CC
      IF (DD.LT.ZERO) GOTO 100   



      SQDD    = SQRT(DD)
      KAPA1   = 0.5D0*(-BB+SQDD)
      KAPA2   = 0.5D0*(-BB-SQDD)
      LAMDA1  = THETA1 + THETA2*KAPA1
      LAMDA2  = THETA1 + THETA2*KAPA2

      IF (KAPA1.GE.ZERO .AND. LAMDA1.GE.ZERO) THEN
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND.   &
             BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND.   &
             BET-KAPA2.GE.ZERO .AND. GAM-LAMDA2.GE.ZERO) THEN
             KAPA = KAPA2
             LAMDA= LAMDA2
             GOTO 200
         ENDIF
      ENDIF



100   KAPA  = ZERO
      LAMDA = ZERO
      DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
      DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)



      IF (DD1.GE.ZERO) THEN
         SQDD1 = SQRT(DD1)
         KAPA1 = 0.5D0*(ALF+BET + SQDD1)
         KAPA2 = 0.5D0*(ALF+BET - SQDD1)

         IF (KAPA1.GE.ZERO .AND. KAPA1.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA1
         ELSE IF (KAPA2.GE.ZERO .AND. KAPA2.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA2
         ELSE
            KAPA = ZERO
         ENDIF
      ENDIF



      IF (DD2.GE.ZERO) THEN
         SQDD2 = SQRT(DD2)
         LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
         LAMDA2= 0.5D0*(ALF+GAM - SQDD2)

         IF (LAMDA1.GE.ZERO .AND. LAMDA1.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1
         ELSE IF (LAMDA2.GE.ZERO .AND. LAMDA2.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
         ELSE
            LAMDA = ZERO
         ENDIF
      ENDIF



      IF (KAPA.GT.ZERO .AND. LAMDA.GT.ZERO) THEN
         IF (BET .LT. LAMDA/THETA1) THEN
            KAPA = ZERO
         ELSE
            LAMDA= ZERO
         ENDIF
      ENDIF



200   CONTINUE
      CNH4NO3 = LAMDA
      CNH4CL  = KAPA

      GNH3    = MAX(ALF - KAPA - LAMDA, ZERO)
      GHNO3   = MAX(GAM - LAMDA, ZERO)
      GHCL    = MAX(BET - KAPA, ZERO)

      RETURN



      END SUBROUTINE CALCG1A


















      SUBROUTINE CALCH6
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      CALAOU = .TRUE.
      CHI1   = W(2)                                
      CHI2   = ZERO                                
      CHI3   = ZERO                                
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
      CHI8   = MIN (FRNA, W(4))                    
      CHI4   = W(3)                                
      CHI5   = MAX (W(4)-CHI8, ZERO)               
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    
      CHI6   = MAX (W(5)-CHI7, ZERO)               

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCH6A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCH6A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH6A (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH6A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCH6')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH6A (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END SUBROUTINE CALCH6






















      DOUBLE PRECISION FUNCTION FUNCH6A (X)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.



      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF



      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               
      MOLAL (3) = PSI4                                  
      MOLAL (4) = PSI6 + PSI7                           
      MOLAL (5) = PSI2 + PSI1                           
      MOLAL (6) = ZERO                                  
      MOLAL (7) = PSI5 + PSI8                           

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)

      CALL CALCMR                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH6A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      RETURN



      END FUNCTION FUNCH6A    



















      SUBROUTINE CALCH5
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      IF (W(4).LE.TINY .AND. W(5).LE.TINY) THEN
         SCASE = 'H5'
         CALL CALCH1A
         SCASE = 'H5'
         RETURN
      ENDIF



      CALAOU = .TRUE.
      CHI1   = W(2)                                
      CHI2   = ZERO                                
      CHI3   = ZERO                                
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
      CHI8   = MIN (FRNA, W(4))                    
      CHI4   = W(3)                                
      CHI5   = MAX (W(4)-CHI8, ZERO)               
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    
      CHI6   = MAX (W(5)-CHI7, ZERO)               

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCH5A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCH5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH5A (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCH5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH5A (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END SUBROUTINE CALCH5






















      DOUBLE PRECISION FUNCTION FUNCH5A (X)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.



      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     
         AA = PSI7+PSI8
         BB = AA*AA
         CC =-A1/4.D0
         CALL POLY3 (AA, BB, CC, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ENDIF



      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                
      MOLAL (3) = PSI4                                   
      MOLAL (4) = PSI6 + PSI7                            
      MOLAL (5) = PSI2 + PSI1                            
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5 + PSI8                            

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)

      CALL CALCMR                               



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      RETURN



      END FUNCTION FUNCH5A    



















      SUBROUTINE CALCH4
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      IF (W(4).LE.TINY .AND. W(5).LE.TINY) THEN
         SCASE = 'H4'
         CALL CALCH1A
         SCASE = 'H4'
         RETURN
      ENDIF



      CALAOU = .TRUE.
      CHI1   = W(2)                                
      CHI2   = ZERO                                
      CHI3   = ZERO                                
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
      CHI8   = MIN (FRNA, W(4))                    
      CHI4   = W(3)                                
      CHI5   = MAX (W(4)-CHI8, ZERO)               
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    
      CHI6   = MAX (W(5)-CHI7, ZERO)               

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCH4A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCH4A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH4A (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH4A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCH4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH4A (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                      
         MOLAL(5) = MOLAL(5) - DELTA                      
         MOLAL(6) = DELTA                                 
      ENDIF

      RETURN



      END SUBROUTINE CALCH4






















      DOUBLE PRECISION FUNCTION FUNCH4A (X)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.



      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     
         AA = PSI7+PSI8
         BB = AA*AA
         CC =-A1/4.D0
         CALL POLY3 (AA, BB, CC, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ENDIF



      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                
      MOLAL (3) = PSI4                                   
      MOLAL (4) = PSI6 + PSI7                            
      MOLAL (5) = PSI2 + PSI1                            
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5 + PSI8                            

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)



      A3   = XK6 /(R*TEMP*R*TEMP)
      DELT = MIN(GNH3, GHCL)
      BB = -(GNH3+GHCL)
      CC = GNH3*GHCL-A3
      DD = BB*BB - 4.D0*CC
      PSI31 = 0.5D0*(-BB + SQRT(DD))
      PSI32 = 0.5D0*(-BB - SQRT(DD))
      IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
         PSI3 = PSI31
      ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
         PSI3 = PSI32
      ELSE
         PSI3 = ZERO
      ENDIF



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR                           



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH4A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      RETURN



      END FUNCTION FUNCH4A    



















      SUBROUTINE CALCH3
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      IF (W(4).LE.TINY) THEN        
         SCASE = 'H3'
         CALL CALCH1A
         SCASE = 'H3'
         RETURN
      ENDIF



      CALAOU = .TRUE.
      CHI1   = W(2)                                
      CHI2   = ZERO                                
      CHI3   = ZERO                                
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
      CHI8   = MIN (FRNA, W(4))                    
      CHI4   = W(3)                                
      CHI5   = MAX (W(4)-CHI8, ZERO)               
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    
      CHI6   = MAX (W(5)-CHI7, ZERO)               

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCH3A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCH3A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH3A (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH3A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCH3')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH3A (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END SUBROUTINE CALCH3






















      DOUBLE PRECISION FUNCTION FUNCH3A (X)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.



      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
         PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF

      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     
         AA = PSI7+PSI8
         BB = AA*AA
         CC =-A1/4.D0
         CALL POLY3 (AA, BB, CC, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ENDIF



      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1             
      MOLAL (3) = PSI4                                
      MOLAL (4) = PSI6 + PSI7                         
      MOLAL (5) = PSI2 + PSI1                         
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5 + PSI8                         

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)



      A3   = XK6 /(R*TEMP*R*TEMP)
      DELT = MIN(GNH3, GHCL)
      BB = -(GNH3+GHCL)
      CC = GNH3*GHCL-A3
      DD = BB*BB - 4.D0*CC
      PSI31 = 0.5D0*(-BB + SQRT(DD))
      PSI32 = 0.5D0*(-BB - SQRT(DD))
      IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
         PSI3 = PSI31
      ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
         PSI3 = PSI32
      ELSE
         PSI3 = ZERO
      ENDIF



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR                                 



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH3A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      RETURN



      END FUNCTION FUNCH3A    



























      SUBROUTINE CALCH2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCH1A, CALCH3



      IF (W(4).GT.TINY) THEN        
         SCASE = 'H2 ; SUBCASE 1'
         CALL CALCH2A
         SCASE = 'H2 ; SUBCASE 1'
      ELSE                          
         SCASE = 'H2 ; SUBCASE 1'
         CALL CALCH1A
         SCASE = 'H2 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY .AND. RH.LT.DRMH2) THEN      
         SCASE = 'H2 ; SUBCASE 2'

      ELSEIF (WATER.LE.TINY .AND. RH.GE.DRMH2) THEN  
         SCASE = 'H2 ; SUBCASE 3'
         CALL CALCMDRH (RH, DRMH2, DRNANO3, CALCH1A, CALCH3)
         SCASE = 'H2 ; SUBCASE 3'
      ENDIF

      RETURN



      END SUBROUTINE CALCH2






















      SUBROUTINE CALCH2A
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      CALAOU = .TRUE.
      CHI1   = W(2)                                
      CHI2   = ZERO                                
      CHI3   = ZERO                                
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
      CHI8   = MIN (FRNA, W(4))                    
      CHI4   = W(3)                                
      CHI5   = MAX (W(4)-CHI8, ZERO)               
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    
      CHI6   = MAX (W(5)-CHI7, ZERO)               

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCH2A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCH2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2 .GT. EPS) Y2 = FUNCH2A (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCH2A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH2A (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    
         MOLAL(5) = MOLAL(5) - DELTA                    
         MOLAL(6) = DELTA                               
      ENDIF

      RETURN



      END SUBROUTINE CALCH2A






















      DOUBLE PRECISION FUNCTION FUNCH2A (X)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A64 = (XK3*XK2/XKW)*(GAMA(10)/GAMA(5)/GAMA(11))**2.0
      A64 = A64*(R*TEMP*WATER)**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.



      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
         PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF

      IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     
         DIAK = (PSI7-PSI5)**2.D0 + 4.D0*A8
         PSI8 = 0.5D0*( -(PSI7+PSI5) + SQRT(DIAK) )
         PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
      ENDIF

      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     
         AA = PSI7+PSI8
         BB = AA*AA
         CC =-A1/4.D0
         CALL POLY3 (AA, BB, CC, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ENDIF



      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                 
      MOLAL (3) = PSI4                                    
      MOLAL (4) = PSI6 + PSI7                             
      MOLAL (5) = PSI2 + PSI1                             
      MOLAL (6) = ZERO                                    
      MOLAL (7) = PSI5 + PSI8                             

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)



      A3   = XK6 /(R*TEMP*R*TEMP)
      DELT = MIN(GNH3, GHCL)
      BB = -(GNH3+GHCL)
      CC = GNH3*GHCL-A3
      DD = BB*BB - 4.D0*CC
      PSI31 = 0.5D0*(-BB + SQRT(DD))
      PSI32 = 0.5D0*(-BB - SQRT(DD))
      IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
         PSI3 = PSI31
      ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
         PSI3 = PSI32
      ELSE
         PSI3 = ZERO
      ENDIF



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR                        



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A64 - ONE

      RETURN



      END FUNCTION FUNCH2A    
























      SUBROUTINE CALCH1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCH1A, CALCH2A



      IF (RH.LT.DRMH1) THEN
         SCASE = 'H1 ; SUBCASE 1'
         CALL CALCH1A              
         SCASE = 'H1 ; SUBCASE 1'
      ELSE
         SCASE = 'H1 ; SUBCASE 2'  
         CALL CALCMDRH (RH, DRMH1, DRNH4NO3, CALCH1A, CALCH2A)
         SCASE = 'H1 ; SUBCASE 2'
      ENDIF

      RETURN



      END SUBROUTINE CALCH1




















      SUBROUTINE CALCH1A
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR,   &
                       NO3FR



      CNA2SO4 = W(2)
      CNH42S4 = ZERO
      NAFR    = MAX (W(1)-2*CNA2SO4, ZERO)
      CNANO3  = MIN (NAFR, W(4))
      NO3FR   = MAX (W(4)-CNANO3, ZERO)
      CNACL   = MIN (MAX(NAFR-CNANO3, ZERO), W(5))
      CLFR    = MAX (W(5)-CNACL, ZERO)



      ALF     = W(3)                     
      BET     = CLFR                     
      GAM     = NO3FR                    

      RTSQ    = R*TEMP*R*TEMP
      A1      = XK6/RTSQ
      A2      = XK10/RTSQ

      THETA1  = GAM - BET*(A2/A1)
      THETA2  = A2/A1



      BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
      CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
      DD      = BB*BB - 4.0D0*CC
      IF (DD.LT.ZERO) GOTO 100   



      SQDD    = SQRT(DD)
      KAPA1   = 0.5D0*(-BB+SQDD)
      KAPA2   = 0.5D0*(-BB-SQDD)
      LAMDA1  = THETA1 + THETA2*KAPA1
      LAMDA2  = THETA1 + THETA2*KAPA2

      IF (KAPA1.GE.ZERO .AND. LAMDA1.GE.ZERO) THEN
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND.   &
             BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND.   &
             BET-KAPA2.GE.ZERO .AND. GAM-LAMDA2.GE.ZERO) THEN
             KAPA = KAPA2
             LAMDA= LAMDA2
             GOTO 200
         ENDIF
      ENDIF



100   KAPA  = ZERO
      LAMDA = ZERO
      DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
      DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)



      IF (DD1.GE.ZERO) THEN
         SQDD1 = SQRT(DD1)
         KAPA1 = 0.5D0*(ALF+BET + SQDD1)
         KAPA2 = 0.5D0*(ALF+BET - SQDD1)

         IF (KAPA1.GE.ZERO .AND. KAPA1.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA1
         ELSE IF (KAPA2.GE.ZERO .AND. KAPA2.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA2
         ELSE
            KAPA = ZERO
         ENDIF
      ENDIF



      IF (DD2.GE.ZERO) THEN
         SQDD2 = SQRT(DD2)
         LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
         LAMDA2= 0.5D0*(ALF+GAM - SQDD2)

         IF (LAMDA1.GE.ZERO .AND. LAMDA1.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1
         ELSE IF (LAMDA2.GE.ZERO .AND. LAMDA2.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
         ELSE
            LAMDA = ZERO
         ENDIF
      ENDIF



      IF (KAPA.GT.ZERO .AND. LAMDA.GT.ZERO) THEN
         IF (BET .LT. LAMDA/THETA1) THEN
            KAPA = ZERO
         ELSE
            LAMDA= ZERO
         ENDIF
      ENDIF



200   CONTINUE
      CNH4NO3 = LAMDA
      CNH4CL  = KAPA

      GNH3    = ALF - KAPA - LAMDA
      GHNO3   = GAM - LAMDA
      GHCL    = BET - KAPA

      RETURN



      END SUBROUTINE CALCH1A


















      SUBROUTINE CALCI6
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CALL CALCI1A



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = CNA2SO4
      PSI5 = CNH42S4

      CALAOU = .TRUE.              
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.



      BB   = PSI2 + PSI4 + PSI5 + A6                    
      CC   =-A6*(PSI2 + PSI3 + PSI1)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))



      MOLAL (1) = PSI6                                    
      MOLAL (2) = 2.D0*PSI4 + PSI3                        
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1            
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6               
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6               
      CLC       = ZERO
      CNAHSO4   = ZERO
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = ZERO
      CNH4HS4   = ZERO
      CALL CALCMR                                         



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE

20    RETURN



      END SUBROUTINE CALCI6



















      SUBROUTINE CALCI5
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CALL CALCI1A



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = CNH42S4

      CALAOU =.TRUE.               
      PSI4LO = ZERO                
      PSI4HI = CHI4                



      IF (CHI4.LE.TINY) THEN
         Y1 = FUNCI5A (ZERO)
         GOTO 50
      ENDIF



      X1 = PSI4HI
      Y1 = FUNCI5A (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCI5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCI5A (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCI5')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCI5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI5A (X3)

50    RETURN



      END SUBROUTINE CALCI5






















      DOUBLE PRECISION FUNCTION FUNCI5A (P4)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      PSI4   = P4     
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4 = XK5 *(WATER/GAMA(2))**3.0
      A5 = XK7 *(WATER/GAMA(4))**3.0
      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.



      BB   = PSI2 + PSI4 + PSI5 + A6                    
      CC   =-A6*(PSI2 + PSI3 + PSI1)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))



      MOLAL (1) = PSI6                            
      MOLAL (2) = 2.D0*PSI4 + PSI3                
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1    
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6       
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6       
      CLC       = ZERO
      CNAHSO4   = ZERO
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = ZERO
      CNH4HS4   = ZERO
      CALL CALCMR                                 



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0
      FUNCI5A= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END FUNCTION FUNCI5A     


















      SUBROUTINE CALCI4
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CALL CALCI1A



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = ZERO

      CALAOU = .TRUE.              
      PSI4LO = ZERO                
      PSI4HI = CHI4                



      IF (CHI4.LE.TINY) THEN
         Y1 = FUNCI4A (ZERO)
         GOTO 50
      ENDIF



      X1 = PSI4HI
      Y1 = FUNCI4A (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCI4A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCI4A (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCI4')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI4A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCI4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI4A (X3)

50    RETURN



      END SUBROUTINE CALCI4






















      DOUBLE PRECISION FUNCTION FUNCI4A (P4)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      PSI4   = P4     
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4 = XK5 *(WATER/GAMA(2))**3.0
      A5 = XK7 *(WATER/GAMA(4))**3.0
      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      A7 = SQRT(A4/A5)



      BB   = PSI2 + PSI4 + PSI5 + A6                    
      CC   =-A6*(PSI2 + PSI3 + PSI1)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))

      PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7
      PSI5 = MIN (PSI5, CHI5)



      MOLAL (1) = PSI6                            
      MOLAL (2) = 2.D0*PSI4 + PSI3                
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1    
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6       
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6       
      CLC       = ZERO
      CNAHSO4   = ZERO
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = CHI5 - PSI5
      CNH4HS4   = ZERO
      CALL CALCMR                                 



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0
      FUNCI4A= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END FUNCTION FUNCI4A     


























      SUBROUTINE CALCI3
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCI1A, CALCI4



      CALL CALCI1A



      IF (CNH4HS4.GT.TINY .OR. CNAHSO4.GT.TINY) THEN
         SCASE = 'I3 ; SUBCASE 1'
         CALL CALCI3A                     
         SCASE = 'I3 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMI3) THEN         
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCI1A
            SCASE = 'I3 ; SUBCASE 2'

         ELSEIF (RH.GE.DRMI3) THEN     
            SCASE = 'I3 ; SUBCASE 3'
            CALL CALCMDRH (RH, DRMI3, DRLC, CALCI1A, CALCI4)
            SCASE = 'I3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCI3





















      SUBROUTINE CALCI3A
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CALL CALCI1A         



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               
      PSI2 = ZERO
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = ZERO

      CALAOU = .TRUE.              
      PSI2LO = ZERO                
      PSI2HI = CHI2                



      X1 = PSI2HI
      Y1 = FUNCI3A (X1)
      YHI= Y1                      



      IF (YHI.LT.EPS) GOTO 50



      DX = (PSI2HI-PSI2LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = FUNCI3A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCI3A (ZERO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI3A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCI3A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI3A (X3)

50    RETURN



      END SUBROUTINE CALCI3A



















      DOUBLE PRECISION FUNCTION FUNCI3A (P2)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      PSI2   = P2                  
      PSI4LO = ZERO                
      PSI4HI = CHI4                



      IF (CHI4.LE.TINY) THEN
         FUNCI3A = FUNCI3B (ZERO)
         GOTO 50
      ENDIF



      X1 = PSI4HI
      Y1 = FUNCI3B (X1)
      IF (ABS(Y1).LE.EPS) GOTO 50
      YHI= Y1                      



      IF (YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI4LO)
         Y2 = FUNCI3B (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCI3B (PSI4LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI3B (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0004, 'FUNCI3A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI3B (X3)



50    A2      = XK13*(WATER/GAMA(13))**5.0
      FUNCI3A = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.D0/A2 - ONE
      RETURN



      END FUNCTION FUNCI3A     


















      DOUBLE PRECISION FUNCTION FUNCI3B (P4)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      PSI4   = P4



      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4 = XK5*(WATER/GAMA(2))**3.0
      A5 = XK7*(WATER/GAMA(4))**3.0
      A6 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      A7 = SQRT(A4/A5)



      BB   = PSI2 + PSI4 + PSI5 + A6                    
      CC   =-A6*(PSI2 + PSI3 + PSI1)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))

      PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7
      PSI5 = MIN (PSI5, CHI5)



      MOLAL(1) = PSI6                                  
      MOLAL(2) = 2.D0*PSI4 + PSI3                      
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1          
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6             
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 - PSI6, TINY)  
      CLC      = MAX(CHI2 - PSI2, ZERO)
      CNAHSO4  = ZERO
      CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
      CNH42S4  = MAX(CHI5 - PSI5, ZERO)
      CNH4HS4  = ZERO
      CALL CALCMR                                       



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0
      FUNCI3B= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END FUNCTION FUNCI3B     


























      SUBROUTINE CALCI2
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCI1A, CALCI3A



      CALL CALCI1A



      IF (CNH4HS4.GT.TINY) THEN
         SCASE = 'I2 ; SUBCASE 1'
         CALL CALCI2A
         SCASE = 'I2 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMI2) THEN         
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCI1A
            SCASE = 'I2 ; SUBCASE 2'

         ELSEIF (RH.GE.DRMI2) THEN     
            SCASE = 'I2 ; SUBCASE 3'
            CALL CALCMDRH (RH, DRMI2, DRNAHSO4, CALCI1A, CALCI3A)
            SCASE = 'I2 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END SUBROUTINE CALCI2




















      SUBROUTINE CALCI2A
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      CALL CALCI1A    



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               
      PSI2 = ZERO
      PSI3 = ZERO
      PSI4 = ZERO
      PSI5 = ZERO

      CALAOU = .TRUE.              
      PSI2LO = ZERO                
      PSI2HI = CHI2                



      X1 = PSI2HI
      Y1 = FUNCI2A (X1)
      YHI= Y1                      



      IF (YHI.LT.EPS) GOTO 50



      DX = (PSI2HI-PSI2LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = FUNCI2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCI3A (ZERO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCI2A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI2A (X3)

50    RETURN



      END SUBROUTINE CALCI2A






















      DOUBLE PRECISION FUNCTION FUNCI2A (P2)
      USE ISRPIA
      USE SOLUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI2   = P2                  
      PSI3   = CHI3
      PSI4   = CHI4
      PSI5   = CHI5
      PSI6   = ZERO



      DO 10 I=1,NSWEEP

      A3 = XK11*(WATER/GAMA(12))**2.0
      A4 = XK5 *(WATER/GAMA(2))**3.0
      A5 = XK7 *(WATER/GAMA(4))**3.0
      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      A7 = SQRT(A4/A5)



      IF (CHI5.GT.TINY .AND. WATER.GT.TINY) THEN
         PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7
         PSI5 = MAX(MIN (PSI5, CHI5), TINY)
      ENDIF

      IF (CHI4.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = PSI2+PSI5+PSI6+PSI3
         BB   = PSI3*AA
         CC   = 0.25D0*(PSI3*PSI3*(PSI2+PSI5+PSI6)-A4)
         CALL POLY3 (AA, BB, CC, PSI4, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI4 = MIN (PSI4, CHI4)
         ELSE
            PSI4 = ZERO
         ENDIF
      ENDIF

      IF (CHI3.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = 2.D0*PSI4 + PSI2 + PSI1 - PSI6
         BB   = 2.D0*PSI4*(PSI2 + PSI1 - PSI6) - A3
         CC   = ZERO
         CALL POLY3 (AA, BB, CC, PSI3, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI3 = MIN (PSI3, CHI3)
         ELSE
            PSI3 = ZERO
         ENDIF
      ENDIF

      BB   = PSI2 + PSI4 + PSI5 + A6                    
      CC   =-A6*(PSI2 + PSI3 + PSI1)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))



      MOLAL (1) = PSI6                           
      MOLAL (2) = 2.D0*PSI4 + PSI3               
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1   
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6      
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6      
      CLC       = CHI2 - PSI2
      CNAHSO4   = CHI3 - PSI3
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = CHI5 - PSI5
      CNH4HS4   = ZERO
      CALL CALCMR                                



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A2      = XK13*(WATER/GAMA(13))**5.0
      FUNCI2A = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.D0/A2 - ONE
      RETURN



      END FUNCTION FUNCI2A     























      SUBROUTINE CALCI1
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL CALCI1A, CALCI2A



      IF (RH.LT.DRMI1) THEN
         SCASE = 'I1 ; SUBCASE 1'
         CALL CALCI1A              
         SCASE = 'I1 ; SUBCASE 1'
      ELSE
         SCASE = 'I1 ; SUBCASE 2'  
         CALL CALCMDRH (RH, DRMI1, DRNH4HS4, CALCI1A, CALCI2A)
         SCASE = 'I1 ; SUBCASE 2'
      ENDIF





      RETURN



      END SUBROUTINE CALCI1




















      SUBROUTINE CALCI1A
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      CNA2SO4 = 0.5D0*W(1)
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      FRSO4   = MAX(W(2)-CNA2SO4, ZERO)

      CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
      FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
      FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)

      IF (FRSO4.LE.TINY) THEN
         CLC     = MAX(CLC - FRNH4, ZERO)
         CNH42S4 = 2.D0*FRNH4

      ELSEIF (FRNH4.LE.TINY) THEN
         CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
         CLC     = MAX(CLC-FRSO4, ZERO)
         IF (CNA2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
      ENDIF



      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO

      RETURN



      END SUBROUTINE CALCI1A

















      SUBROUTINE CALCJ3
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      DOUBLE PRECISION LAMDA, KAPA



      CALAOU = .TRUE.              
      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  
      CHI1   = W(1)                           
      CHI2   = W(3)                           
      PSI1   = CHI1
      PSI2   = CHI2                           



      DO 10 I=1,NSWEEP

      A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0



      BB   = A3+LAMDA                        
      CC   =-A3*(LAMDA + PSI1 + PSI2)
      DD   = BB*BB-4.D0*CC
      KAPA = 0.5D0*(-BB+SQRT(DD))



      MOLAL (1) = LAMDA + KAPA                 
      MOLAL (2) = PSI1                         
      MOLAL (3) = PSI2                         
      MOLAL (4) = ZERO                         
      MOLAL (5) = KAPA                         
      MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA   
      MOLAL (7) = ZERO                         

      CNAHSO4   = ZERO
      CNH4HS4   = ZERO

      CALL CALCMR                              



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 50
      ENDIF
10    CONTINUE

50    RETURN



      END SUBROUTINE CALCJ3


















      SUBROUTINE CALCJ2
      USE ISRPIA
      USE CASEJ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      CALAOU = .TRUE.              
      CHI1   = W(1)                
      CHI2   = W(3)                
      PSI1LO = TINY                
      PSI1HI = CHI1                



      X1 = PSI1HI
      Y1 = FUNCJ2 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI1HI-PSI1LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCJ2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCJ2 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCJ2')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCJ2 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCJ2')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCJ2 (X3)

50    RETURN



      END SUBROUTINE CALCJ2






















      DOUBLE PRECISION FUNCTION FUNCJ2 (P1)
      USE CASEJ
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  
      PSI1   = P1
      PSI2   = CHI2                           



      DO 10 I=1,NSWEEP

      A1 = XK11 *(WATER/GAMA(12))**2.0
      A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0



      BB   = A3+LAMDA                        
      CC   =-A3*(LAMDA + PSI1 + PSI2)
      DD   = BB*BB-4.D0*CC
      KAPA = 0.5D0*(-BB+SQRT(DD))



      MOLAL (1) = LAMDA + KAPA                  
      MOLAL (2) = PSI1                          
      MOLAL (3) = PSI2                          
      MOLAL (4) = ZERO                          
      MOLAL (5) = KAPA                          
      MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA    
      MOLAL (7) = ZERO                          

      CNAHSO4   = MAX(CHI1-PSI1,ZERO)
      CNH4HS4   = ZERO

      CALL CALCMR                               



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCJ2 = MOLAL(2)*MOLAL(6)/A1 - ONE



      END FUNCTION FUNCJ2     



















      SUBROUTINE CALCJ1
      USE ISRPIA
      USE CASEJ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      CALAOU =.TRUE.               
      CHI1   = W(1)                
      CHI2   = W(3)                

      PSI1LO = TINY                
      PSI1HI = CHI1                



      X1 = PSI1HI
      Y1 = FUNCJ1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI1HI-PSI1LO)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCJ1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCJ1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCJ1')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCJ1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCJ1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCJ1 (X3)

50    RETURN



      END SUBROUTINE CALCJ1






















      DOUBLE PRECISION FUNCTION FUNCJ1 (P1)
      USE ISRPIA
      USE CASEJ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)







      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  
      PSI1   = P1



      DO 10 I=1,NSWEEP

      A1 = XK11 *(WATER/GAMA(12))**2.0
      A2 = XK12 *(WATER/GAMA(09))**2.0
      A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0

      PSI2 = 0.5*(-(LAMDA+PSI1) + SQRT((LAMDA+PSI1)**2.D0+4.D0*A2))  
      PSI2 = MIN (PSI2, CHI2)

      BB   = A3+LAMDA                        
      CC   =-A3*(LAMDA + PSI2 + PSI1)
      DD   = BB*BB-4.D0*CC
      KAPA = 0.5D0*(-BB+SQRT(DD))



      MOLAL (1) = LAMDA + KAPA                  
      MOLAL (2) = PSI1                          
      MOLAL (3) = PSI2                          
      MOLAL (4) = ZERO
      MOLAL (5) = KAPA                          
      MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA    
      MOLAL (7) = ZERO

      CNAHSO4   = MAX(CHI1-PSI1,ZERO)
      CNH4HS4   = MAX(CHI2-PSI2,ZERO)

      CALL CALCMR                               



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCJ1 = MOLAL(2)*MOLAL(6)/A1 - ONE



      END FUNCTION FUNCJ1     

