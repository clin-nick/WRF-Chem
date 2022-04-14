




























































































































      SUBROUTINE ISOROPIA (WI, RHI, TEMPI,  CNTRL,    &
                           WT, GAS, AERLIQ, AERSLD, SCASI, OTHER)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NCTRL=2,NOTHER=6)
      CHARACTER SCASI*15
      DIMENSION WI(NCOMP), WT(NCOMP),   GAS(NGASAQ),  AERSLD(NSLDS),    &
                AERLIQ(NIONS+NGASAQ+2), CNTRL(NCTRL), OTHER(NOTHER)



      IPROB   = NINT(CNTRL(1))



      METSTBL = NINT(CNTRL(2))



50    IF (IPROB.EQ.0) THEN
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN 
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        
            CALL ISRP1F (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              
            CALL ISRP2F (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3F (WI, RHI, TEMPI)
         ENDIF



      ELSE
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN 
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        
            CALL ISRP1R (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              
            CALL ISRP2R (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3R (WI, RHI, TEMPI)
         ENDIF
      ENDIF



      IF (NADJ.EQ.1) CALL ADJUST (WI)










      GAS(1) = GNH3                
      GAS(2) = GHNO3
      GAS(3) = GHCL

      DO 10 I=1,NIONS              
         AERLIQ(I) = MOLAL(I)
  10  CONTINUE
      DO 20 I=1,NGASAQ
         AERLIQ(NIONS+1+I) = GASAQ(I)
  20  CONTINUE
      AERLIQ(NIONS+1)        = WATER*1.0D3/18.0D0
      AERLIQ(NIONS+NGASAQ+2) = COH

      AERSLD(1) = CNANO3           
      AERSLD(2) = CNH4NO3
      AERSLD(3) = CNACL
      AERSLD(4) = CNH4CL
      AERSLD(5) = CNA2SO4
      AERSLD(6) = CNH42S4
      AERSLD(7) = CNAHSO4
      AERSLD(8) = CNH4HS4
      AERSLD(9) = CLC

      IF(WATER.LE.TINY) THEN       
        OTHER(1) = 1.d0
      ELSE
        OTHER(1) = 0.d0
      ENDIF

      OTHER(2) = SULRAT            
      OTHER(3) = SULRATW
      OTHER(4) = SODRAT
      OTHER(5) = IONIC
      OTHER(6) = ICLACT

      SCASI = SCASE

      WT(1) = WI(1)                
      WT(2) = WI(2)
      WT(3) = WI(3) 
      WT(4) = WI(4)
      WT(5) = WI(5)
      IF (IPROB.GT.0 .AND. WATER.GT.TINY) THEN 
         WT(3) = WT(3) + GNH3 
         WT(4) = WT(4) + GHNO3
         WT(5) = WT(5) + GHCL
      ENDIF

      RETURN



      END







































































      SUBROUTINE SETPARM (WFTYPI,  IACALCI, EPSI, MAXITI, NSWEEPI,    &
                          EPSACTI, NDIVI, NADJI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER  WFTYPI



      IF (WFTYPI .GE. 0)   WFTYP  = WFTYPI
      IF (IACALCI.GE. 0)   IACALC = IACALCI
      IF (EPSI   .GE.ZERO) EPS    = EPSI
      IF (MAXITI .GT. 0)   MAXIT  = MAXITI
      IF (NSWEEPI.GT. 0)   NSWEEP = NSWEEPI
      IF (EPSACTI.GE.ZERO) EPSACT = EPSACTI
      IF (NDIVI  .GT. 0)   NDIV   = NDIVI
      IF (NADJI  .GE. 0)   NADJ   = NADJI



      RETURN
      END






















      SUBROUTINE GETPARM (WFTYPI,  IACALCI, EPSI, MAXITI, NSWEEPI,    &
                          EPSACTI, NDIVI, NADJI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER  WFTYPI



      WFTYPI  = WFTYP
      IACALCI = IACALC
      EPSI    = EPS
      MAXITI  = MAXIT
      NSWEEPI = NSWEEP
      EPSACTI = EPSACT
      NDIVI   = NDIV
      NADJI   = NADJ



      RETURN
      END



















      BLOCK DATA BLKISO
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)








      END














      SUBROUTINE INIT1 (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)
      REAL      IC,GII,GI0,XX,LN10
      PARAMETER (LN10=2.3025851)



      IF (IPROB.EQ.0) THEN                 
         DO 10 I=1,NCOMP
            W(I) = MAX(WI(I), TINY)
10       CONTINUE
      ELSE
         DO 15 I=1,NCOMP                   
            WAER(I) = MAX(WI(I), TINY)
            W(I)    = ZERO
15       CONTINUE
      ENDIF
      RH      = RHI
      TEMP    = TEMPI



      XK1  = 1.015e-2  
      XK21 = 57.639    
      XK22 = 1.805e-5  
      XK7  = 1.817     
      XK12 = 1.382e2   
      XK13 = 29.268    
      XKW  = 1.010e-14 

      IF (INT(TEMP) .NE. 298) THEN   
         T0  = 298.15
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2 = XK21*XK22       



      DRH2SO4  = 0.0000D0
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRLC     = 0.6900D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15d0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF) 
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF) 
         DRLC     = DRLC    *EXP(186.*TCF) 
      ENDIF



      DRMLCAB = 0.3780D0              
      DRMLCAS = 0.6900D0              









      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY

      DO 20 I=1,NPAIR
         MOLALR(I)=ZERO
         GAMA(I)  =0.1
         GAMIN(I) =GREAT
         GAMOU(I) =GREAT
         M0(I)    =1d5
 20   CONTINUE

      DO 30 I=1,NPAIR
         GAMA(I) = 0.1d0
 30   CONTINUE

      DO 40 I=1,NIONS
         MOLAL(I)=ZERO
40    CONTINUE
      COH = ZERO

      DO 50 I=1,NGASAQ
         GASAQ(I)=ZERO
50    CONTINUE



      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO



      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO



      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  
      IRH    = MAX (IRH, 1)

























      M0(04) = AWAS(IRH)      























      M0(07) = AWSA(IRH)      







      M0(08) = AWSA(IRH)      








      M0(09) = AWAB(IRH)      















      M0(13) = AWLC(IRH)      











      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = '??'
      SULRATW = 2.D0
      SODRAT  = ZERO
      NOFER   = 0
      STKOFL  =.FALSE.
      DO 60 I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
   60 CONTINUE



      END

















      SUBROUTINE INIT2 (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)
      REAL      IC,GII,GI0,XX,LN10
      PARAMETER (LN10=2.3025851)



      IF (IPROB.EQ.0) THEN                 
         DO 10 I=1,NCOMP
            W(I) = MAX(WI(I), TINY)
10       CONTINUE
      ELSE
         DO 15 I=1,NCOMP                   
            WAER(I) = MAX(WI(I), TINY)
            W(I)    = ZERO
15       CONTINUE
      ENDIF
      RH      = RHI
      TEMP    = TEMPI



      XK1  = 1.015e-2  
      XK21 = 57.639    
      XK22 = 1.805e-5  
      XK4  = 2.511e6   

      XK41 = 2.100e5   
      XK7  = 1.817     
      XK10 = 5.746e-17 

      XK12 = 1.382e2   
      XK13 = 29.268    
      XKW  = 1.010e-14 

      IF (INT(TEMP) .NE. 298) THEN   
         T0  = 298.15D0
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) 

         XK41= XK41*EXP( 29.17*(T0T-1.0) + 16.830*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) 

         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2  = XK21*XK22       
      XK42 = XK4/XK41



      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRNH4NO3 = 0.6183D0
      DRLC     = 0.6900D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15D0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNH4NO3 = DRNH4NO3*EXP(852.*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF) 
         DRLC     = DRLC    *EXP(186.*TCF) 
         DRNH4NO3 = MIN (DRNH4NO3,DRNH42S4) 
      ENDIF



      DRMLCAB = 0.3780D0              
      DRMLCAS = 0.6900D0              
      DRMASAN = 0.6000D0              










      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY

      DO 20 I=1,NPAIR
         MOLALR(I)=ZERO
         GAMA(I)  =0.1
         GAMIN(I) =GREAT
         GAMOU(I) =GREAT
         M0(I)    =1d5
 20   CONTINUE

      DO 30 I=1,NPAIR
         GAMA(I) = 0.1d0
 30   CONTINUE

      DO 40 I=1,NIONS
         MOLAL(I)=ZERO
40    CONTINUE
      COH = ZERO

      DO 50 I=1,NGASAQ
         GASAQ(I)=ZERO
50    CONTINUE



      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO



      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO



      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  
      IRH    = MAX (IRH, 1)

























      M0(04) = AWAS(IRH)      







      M0(05) = AWAN(IRH)      















      M0(07) = AWSA(IRH)      







      M0(08) = AWSA(IRH)      








      M0(09) = AWAB(IRH)      















      M0(13) = AWLC(IRH)      











      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = '??'
      SULRATW = 2.D0
      SODRAT  = ZERO
      NOFER   = 0
      STKOFL  =.FALSE.
      DO 60 I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
   60 CONTINUE



      END




















      SUBROUTINE ISOINIT3 (WI, RHI, TEMPI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)
      REAL      IC,GII,GI0,XX,LN10
      PARAMETER (LN10=2.3025851)



      IF (IPROB.EQ.0) THEN                 
         DO 10 I=1,NCOMP
            W(I) = MAX(WI(I), TINY)
10       CONTINUE
      ELSE
         DO 15 I=1,NCOMP                   
            WAER(I) = MAX(WI(I), TINY)
            W(I)    = ZERO
15       CONTINUE
      ENDIF
      RH      = RHI
      TEMP    = TEMPI



      XK1  = 1.015D-2  
      XK21 = 57.639D0  
      XK22 = 1.805D-5  
      XK3  = 1.971D6   
      XK31 = 2.500e3   
      XK4  = 2.511e6   

      XK41 = 2.100e5   
      XK5  = 0.4799D0  
      XK6  = 1.086D-16 
      XK7  = 1.817D0   
      XK8  = 37.661D0  
      XK10 = 5.746D-17 

      XK11 = 2.413D4   
      XK12 = 1.382D2   
      XK13 = 29.268D0  
      XK14 = 22.05D0   
      XKW  = 1.010D-14 
      XK9  = 11.977D0  

      IF (INT(TEMP) .NE. 298) THEN   
         T0  = 298.15D0
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK3 = XK3 *EXP( 30.20*(T0T-1.0) + 19.910*COEF)
         XK31= XK31*EXP( 30.20*(T0T-1.0) + 19.910*COEF)
         XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) 

         XK41= XK41*EXP( 29.17*(T0T-1.0) + 16.830*COEF)
         XK5 = XK5 *EXP(  0.98*(T0T-1.0) + 39.500*COEF)
         XK6 = XK6 *EXP(-71.00*(T0T-1.0) +  2.400*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK8 = XK8 *EXP( -1.56*(T0T-1.0) + 16.900*COEF)
         XK9 = XK9 *EXP( -8.22*(T0T-1.0) + 16.010*COEF)
         XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) 

         XK11= XK11*EXP(  0.79*(T0T-1.0) + 14.746*COEF)
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XK14= XK14*EXP( 24.55*(T0T-1.0) + 16.900*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2  = XK21*XK22       
      XK42 = XK4/XK41
      XK32 = XK3/XK31



      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRLC     = 0.6900D0
      DRNACL   = 0.7528D0
      DRNANO3  = 0.7379D0
      DRNH4CL  = 0.7710D0
      DRNH4NO3 = 0.6183D0
      DRNA2SO4 = 0.9300D0
      DRNAHSO4 = 0.5200D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15D0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNACL   = DRNACL  *EXP( 25.*TCF)
         DRNANO3  = DRNANO3 *EXP(304.*TCF)
         DRNA2SO4 = DRNA2SO4*EXP( 80.*TCF)
         DRNH4NO3 = DRNH4NO3*EXP(852.*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF) 
         DRLC     = DRLC    *EXP(186.*TCF)
         DRNH4CL  = DRNH4Cl *EXP(239.*TCF)
         DRNAHSO4 = DRNAHSO4*EXP(-45.*TCF) 



         DRNH4NO3  = MIN (DRNH4NO3, DRNH4CL, DRNH42S4, DRNANO3, DRNACL)
         DRNANO3   = MIN (DRNANO3, DRNACL)
         DRNH4CL   = MIN (DRNH4Cl, DRNH42S4)

      ENDIF



      DRMLCAB = 0.378D0    
      DRMLCAS = 0.690D0    
      DRMASAN = 0.600D0    
      DRMG1   = 0.460D0    
      DRMG2   = 0.691D0    
      DRMG3   = 0.697D0    
      DRMH1   = 0.240D0    
      DRMH2   = 0.596D0    
      DRMI1   = 0.240D0    
      DRMI2   = 0.363D0    
      DRMI3   = 0.610D0    
      DRMQ1   = 0.494D0    
      DRMR1   = 0.663D0    
      DRMR2   = 0.735D0    
      DRMR3   = 0.673D0    
      DRMR4   = 0.694D0    
      DRMR5   = 0.731D0    
      DRMR6   = 0.596D0    
      DRMR7   = 0.380D0    
      DRMR8   = 0.380D0    
      DRMR9   = 0.494D0    
      DRMR10  = 0.476D0    
      DRMR11  = 0.340D0    
      DRMR12  = 0.460D0    
      DRMR13  = 0.438D0    
































      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY

      DO 20 I=1,NPAIR
         MOLALR(I)=ZERO
         GAMA(I)  =0.1
         GAMIN(I) =GREAT
         GAMOU(I) =GREAT
         M0(I)    =1d5
 20   CONTINUE

      DO 30 I=1,NPAIR
         GAMA(I) = 0.1d0
 30   CONTINUE

      DO 40 I=1,NIONS
         MOLAL(I)=ZERO
40    CONTINUE
      COH = ZERO

      DO 50 I=1,NGASAQ
         GASAQ(I)=ZERO
50    CONTINUE



      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO



      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO



      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  
      IRH    = MAX (IRH, 1)

      M0(01) = AWSC(IRH)      







      M0(02) = AWSS(IRH)      







      M0(03) = AWSN(IRH)      







      M0(04) = AWAS(IRH)      







      M0(05) = AWAN(IRH)      







      M0(06) = AWAC(IRH)      







      M0(07) = AWSA(IRH)      







      M0(08) = AWSA(IRH)      








      M0(09) = AWAB(IRH)      







      M0(12) = AWSB(IRH)      







      M0(13) = AWLC(IRH)      











      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = '??'
      SULRATW = 2.D0
      NOFER   = 0
      STKOFL  =.FALSE.
      DO 60 I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
   60 CONTINUE



      END
      
















      SUBROUTINE ADJUST (WI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION WI(*)



      IF (IPROB.EQ.0) THEN         
         EXNH4 = GNH3 + MOLAL(3) + CNH4CL + CNH4NO3 + CNH4HS4   &
                      + 2D0*CNH42S4       + 3D0*CLC   &
                -WI(3)
      ELSE
         EXNH4 = MOLAL(3) + CNH4CL + CNH4NO3 + CNH4HS4 + 2D0*CNH42S4   &
                          + 3D0*CLC   &
                -WI(3)

      ENDIF
      EXNH4 = MAX(EXNH4,ZERO)
      IF (EXNH4.LT.TINY) GOTO 20    

      IF (MOLAL(3).GT.EXNH4) THEN   
         MOLAL(3) = MOLAL(3) - EXNH4
         GOTO 20
      ELSE
         EXNH4    = EXNH4 - MOLAL(3)
         MOLAL(3) = ZERO
      ENDIF

      IF (CNH4CL.GT.EXNH4) THEN     
         CNH4CL   = CNH4CL - EXNH4  
         GHCL     = GHCL   + EXNH4  
         GOTO 20
      ELSE                          
         GHCL     = GHCL   + CNH4CL 
         EXNH4    = EXNH4  - CNH4CL 
         CNH4CL   = ZERO            
      ENDIF

      IF (CNH4NO3.GT.EXNH4) THEN    
         CNH4NO3  = CNH4NO3- EXNH4  
         GHNO3    = GHNO3  + EXNH4  
         GOTO 20
      ELSE                          
         GHNO3    = GHNO3  + CNH4NO3
         EXNH4    = EXNH4  - CNH4NO3
         CNH4NO3  = ZERO            
      ENDIF

      IF (CLC.GT.3d0*EXNH4) THEN    
         CLC      = CLC - EXNH4/3d0 
         GOTO 20
      ELSE                          
         EXNH4    = EXNH4 - 3d0*CLC 
         CLC      = ZERO            
      ENDIF

      IF (CNH4HS4.GT.EXNH4) THEN    
         CNH4HS4  = CNH4HS4- EXNH4  
         GOTO 20
      ELSE                          
         EXNH4    = EXNH4  - CNH4HS4
         CNH4HS4  = ZERO            
      ENDIF

      IF (CNH42S4.GT.EXNH4) THEN    
         CNH42S4  = CNH42S4- EXNH4  
         GOTO 20
      ELSE                          
         EXNH4    = EXNH4  - CNH42S4
         CNH42S4  = ZERO            
      ENDIF



 20   IF (IPROB.EQ.0) THEN         
         EXNO3 = GHNO3 + MOLAL(7) + CNH4NO3   &
                -WI(4)
      ELSE
         EXNO3 = MOLAL(7) + CNH4NO3   &
                -WI(4)
      ENDIF
      EXNO3 = MAX(EXNO3,ZERO)
      IF (EXNO3.LT.TINY) GOTO 30    

      IF (MOLAL(7).GT.EXNO3) THEN   
         MOLAL(7) = MOLAL(7) - EXNO3
         GOTO 30
      ELSE
         EXNO3    = EXNO3 - MOLAL(7)
         MOLAL(7) = ZERO
      ENDIF

      IF (CNH4NO3.GT.EXNO3) THEN    
         CNH4NO3  = CNH4NO3- EXNO3  
         GNH3     = GNH3   + EXNO3  
         GOTO 30
      ELSE                          
         GNH3     = GNH3   + CNH4NO3
         EXNO3    = EXNO3  - CNH4NO3
         CNH4NO3  = ZERO            
      ENDIF



 30   IF (IPROB.EQ.0) THEN         
         EXCl = GHCL + MOLAL(4) + CNH4CL   &
               -WI(5)
      ELSE
         EXCl = MOLAL(4) + CNH4CL   &
               -WI(5)
      ENDIF
      EXCl = MAX(EXCl,ZERO)
      IF (EXCl.LT.TINY) GOTO 40    

      IF (MOLAL(4).GT.EXCL) THEN   
         MOLAL(4) = MOLAL(4) - EXCL
         GOTO 40
      ELSE
         EXCL     = EXCL - MOLAL(4)
         MOLAL(4) = ZERO
      ENDIF

      IF (CNH4CL.GT.EXCL) THEN      
         CNH4CL   = CNH4CL - EXCL   
         GHCL     = GHCL   + EXCL   
         GOTO 40
      ELSE                          
         GHCL     = GHCL   + CNH4CL 
         EXCL     = EXCL   - CNH4CL 
         CNH4CL   = ZERO            
      ENDIF



 40   EXS4 = MOLAL(5) + MOLAL(6) + 2.d0*CLC + CNH42S4 + CNH4HS4 +   &
             CNA2SO4  + CNAHSO4 - WI(2)
      EXS4 = MAX(EXS4,ZERO)        
      IF (EXS4.LT.TINY) GOTO 50    

      IF (MOLAL(6).GT.EXS4) THEN   
         MOLAL(6) = MOLAL(6) - EXS4
         GOTO 50
      ELSE
         EXS4     = EXS4 - MOLAL(6)
         MOLAL(6) = ZERO
      ENDIF

      IF (MOLAL(5).GT.EXS4) THEN   
         MOLAL(5) = MOLAL(5) - EXS4
         GOTO 50
      ELSE
         EXS4     = EXS4 - MOLAL(5)
         MOLAL(5) = ZERO
      ENDIF

      IF (CLC.GT.2d0*EXS4) THEN     
         CLC      = CLC - EXS4/2d0  
         GNH3     = GNH3 +1.5d0*EXS4
         GOTO 50
      ELSE                          
         GNH3     = GNH3 + 1.5d0*CLC
         EXS4     = EXS4 - 2d0*CLC  
         CLC      = ZERO            
      ENDIF

      IF (CNH4HS4.GT.EXS4) THEN     
         CNH4HS4  = CNH4HS4 - EXS4  
         GNH3     = GNH3 + EXS4     
         GOTO 50
      ELSE                          
         GNH3     = GNH3 + CNH4HS4  
         EXS4     = EXS4  - CNH4HS4 
         CNH4HS4  = ZERO            
      ENDIF

      IF (CNH42S4.GT.EXS4) THEN     
         CNH42S4  = CNH42S4- EXS4   
         GNH3     = GNH3 + 2.d0*EXS4
         GOTO 50
      ELSE                          
         GNH3     = GNH3+2.d0*CNH42S4 
         EXS4     = EXS4  - CNH42S4 
         CNH42S4  = ZERO            
      ENDIF



 50   RETURN
      END
      














      DOUBLE PRECISION FUNCTION GETASR (SO4I, RHI)
      USE ASRC
      PARAMETER (NSO4S=14, NRHS=20, NASRD=NSO4S*NRHS)

      DOUBLE PRECISION SO4I, RHI











      RAT    = SO4I/1.E-9    
      A1     = INT(ALOG10(RAT))                   
      IA1    = INT(RAT/2.5/10.0**A1)

      INDS   = 4.0*A1 + MIN(IA1,4)
      INDS   = MIN(MAX(0, INDS), NSO4S-1) + 1     

      INDR   = INT(99.0-RHI*100.0) + 1
      INDR   = MIN(MAX(1, INDR), NRHS)            



      INDSL  = INDS
      INDSH  = MIN(INDSL+1, NSO4S)
      IPOSL  = (INDSL-1)*NRHS + INDR              
      IPOSH  = (INDSH-1)*NRHS + INDR              

      WF     = (SO4I-ASSO4(INDSL))/(ASSO4(INDSH)-ASSO4(INDSL) + 1e-7)
      WF     = MIN(MAX(WF, 0.0), 1.0)

      GETASR = WF*ASRAT(IPOSH) + (1.0-WF)*ASRAT(IPOSL)



      RETURN
      END

















      BLOCK DATA AERSR
      USE ASRC
      PARAMETER (NSO4S=14, NRHS=20, NASRD=NSO4S*NRHS)






































































       END
      
       



















      SUBROUTINE CALCHA
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION KAPA




      X    = W(5) 
      DELT = 0.0d0
      IF (WATER.GT.TINY) THEN
         KAPA = MOLAL(1)
         ALFA = XK3*R*TEMP*(WATER/GAMA(11))**2.0
         DIAK = SQRT( (KAPA+ALFA)**2.0 + 4.0*ALFA*X)
         DELT = 0.5*(-(KAPA+ALFA) + DIAK)




      ENDIF



      GHCL     = MAX(X-DELT, 0.0d0)  



      MOLAL(4) = DELT                
      MOLAL(1) = MOLAL(1) + DELT     

      RETURN



      END


























      SUBROUTINE CALCHAP
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (WATER.LE.TINY) RETURN



      CALL CALCCLAQ (MOLAL(4), MOLAL(1), DELT)
      ALFA     = XK3*R*TEMP*(WATER/GAMA(11))**2.0
      GASAQ(3) = DELT
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(4) = MOLAL(4) - DELT
      GHCL     = MOLAL(1)*MOLAL(4)/ALFA

      RETURN



      END




















      SUBROUTINE CALCNA
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION KAPA




      X    = W(4) 
      DELT = 0.0d0
      IF (WATER.GT.TINY) THEN
         KAPA = MOLAL(1)
         ALFA = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         DIAK = SQRT( (KAPA+ALFA)**2.0 + 4.0*ALFA*X)
         DELT = 0.5*(-(KAPA+ALFA) + DIAK)




      ENDIF



      GHNO3    = MAX(X-DELT, 0.0d0)  



      MOLAL(7) = DELT                
      MOLAL(1) = MOLAL(1) + DELT     

      RETURN



      END























      SUBROUTINE CALCNAP
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (WATER.LE.TINY) RETURN



      CALL CALCNIAQ (MOLAL(7), MOLAL(1), DELT)
      ALFA     = XK4*R*TEMP*(WATER/GAMA(10))**2.0
      GASAQ(3) = DELT
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(7) = MOLAL(7) - DELT
      GHNO3    = MOLAL(1)*MOLAL(7)/ALFA
      


      RETURN



      END





















      SUBROUTINE CALCNH3
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (WATER.LE.TINY) RETURN



      A1   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      CHI1 = MOLAL(3)
      CHI2 = MOLAL(1)

      BB   =(CHI2 + ONE/A1)          
      CC   =-CHI1/A1             
      DIAK = SQRT(BB*BB - 4.D0*CC)   
      PSI  = 0.5*(-BB + DIAK)        
      PSI  = MAX(TINY, MIN(PSI,CHI1))



      GNH3     = PSI                 



      MOLAL(3) = CHI1 - PSI          
      MOLAL(1) = CHI2 + PSI          

      RETURN



      END





















      SUBROUTINE CALCNH3P
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (WATER.LE.TINY) RETURN



      A1   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      GNH3 = MOLAL(3)/MOLAL(1)/A1

      RETURN



      END


















      SUBROUTINE CALCNHA
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION M1, M2, M3
      CHARACTER ERRINF*40     



      IF (WATER.LE.TINY) THEN
         GOTO 55



      ELSEIF (W(5).LE.TINY .AND. W(4).LE.TINY) THEN
         GOTO 60



      ELSE IF (W(5).LE.TINY) THEN
         CALL CALCNA              
         GOTO 60



      ELSE IF (W(4).LE.TINY) THEN
         CALL CALCHA              
         GOTO 60
      ENDIF



      A3 = XK4*R*TEMP*(WATER/GAMA(10))**2.0   
      A4 = XK3*R*TEMP*(WATER/GAMA(11))**2.0   



      DELCL = ZERO
      DELNO = ZERO

      OMEGA = MOLAL(1)       
      CHI3  = W(4)           
      CHI4  = W(5)           

      C1    = A3*CHI3
      C2    = A4*CHI4
      C3    = A3 - A4

      M1    = (C1 + C2 + (OMEGA+A4)*C3)/C3
      M2    = ((OMEGA+A4)*C2 - A4*C3*CHI4)/C3
      M3    =-A4*C2*CHI4/C3



      CALL POLY3 (M1, M2, M3, DELCL, ISLV) 
      IF (ISLV.NE.0) THEN
         DELCL = TINY       
         WRITE (ERRINF,'(1PE7.1)') TINY
         CALL PUSHERR (0022, ERRINF)    
      ENDIF
      DELCL = MIN(DELCL, CHI4)

      DELNO = C1*DELCL/(C2 + C3*DELCL)  
      DELNO = MIN(DELNO, CHI3)

      IF (DELCL.LT.ZERO .OR. DELNO.LT.ZERO .OR.   &
         DELCL.GT.CHI4 .OR. DELNO.GT.CHI3       ) THEN
         DELCL = TINY  
         DELNO = TINY
         WRITE (ERRINF,'(1PE7.1)') TINY
         CALL PUSHERR (0022, ERRINF)    
      ENDIF










50    MOLAL(1) = MOLAL(1) + (DELNO+DELCL)  
      MOLAL(4) = MOLAL(4) + DELCL          
      MOLAL(7) = MOLAL(7) + DELNO          



55    GHCL     = MAX(W(5) - MOLAL(4), TINY)
      GHNO3    = MAX(W(4) - MOLAL(7), TINY)

60    RETURN



      END





















      SUBROUTINE CALCNHP
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      IF (WATER.LE.TINY) RETURN



      A3       = XK3*R*TEMP*(WATER/GAMA(11))**2.0
      A4       = XK4*R*TEMP*(WATER/GAMA(10))**2.0
      MOLAL(1) = MOLAL(1) + WAER(4) + WAER(5)  




      CALL CALCNIAQ (WAER(4), MOLAL(1)+MOLAL(7)+MOLAL(4), DELT)
      MOLAL(1) = MOLAL(1) - DELT 
      MOLAL(7) = WAER(4)  - DELT  
      GASAQ(3) = DELT

      CALL CALCCLAQ (WAER(5), MOLAL(1)+MOLAL(7)+MOLAL(4), DELT)
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(4) = WAER(5)  - DELT  
      GASAQ(2) = DELT

      GHNO3    = MOLAL(1)*MOLAL(7)/A4
      GHCL     = MOLAL(1)*MOLAL(4)/A3

      RETURN



      END














      SUBROUTINE CALCAMAQ (NH4I, OHI, DELT)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION NH4I




      A22  = XK22/XKW/WATER*(GAMA(8)/GAMA(9))**2. 
      AKW  = XKW *RH*WATER*WATER



      OM1  = NH4I          
      OM2  = OHI
      BB   =-(OM1+OM2+A22*AKW)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)



      IF (DEL1.LT.ZERO) THEN                 
         IF (DEL2.GT.NH4I .OR. DEL2.GT.OHI) THEN
            DELT = ZERO
         ELSE
            DELT = DEL2
         ENDIF
      ELSE
         DELT = DEL1
      ENDIF








      RETURN



      END

















      SUBROUTINE CALCAMAQ2 (GGNH3, NH4I, OHI, NH3AQ)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION NH4I, NH3AQ



      A22  = XK22/XKW/WATER*(GAMA(8)/GAMA(9))**2. 
      AKW  = XKW *RH*WATER*WATER



      ALF1 = NH4I - GGNH3
      ALF2 = GGNH3
      BB   = ALF1 + A22*AKW
      CC   =-A22*AKW*ALF2
      DEL  = 0.5D0*(-BB + SQRT(BB*BB-4.D0*CC))



      NH4I  = ALF1 + DEL
      OHI   = DEL
      IF (OHI.LE.TINY) OHI = SQRT(AKW)   
      NH3AQ = ALF2 - DEL 

      RETURN



      END

















      SUBROUTINE CALCCLAQ (CLI, HI, DELT)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION CLI



      A32  = XK32*WATER/(GAMA(11))**2. 



      OM1  = CLI          
      OM2  = HI
      BB   =-(OM1+OM2+A32)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)



      IF (DEL1.LT.ZERO) THEN                 
         IF (DEL2.LT.ZERO .OR. DEL2.GT.CLI .OR. DEL2.GT.HI) THEN
            DELT = ZERO
         ELSE
            DELT = DEL2
         ENDIF
      ELSE
         DELT = DEL1
      ENDIF

      RETURN



      END

















      SUBROUTINE CALCCLAQ2 (GGCL, CLI, HI, CLAQ)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION CLI



      A32  = XK32*WATER/(GAMA(11))**2. 
      AKW  = XKW *RH*WATER*WATER



      ALF1  = CLI - GGCL
      ALF2  = GGCL
      COEF  = (ALF1+A32)
      DEL1  = 0.5*(-COEF + SQRT(COEF*COEF+4.D0*A32*ALF2))



      CLI  = ALF1 + DEL1
      HI   = DEL1
      IF (HI.LE.TINY) HI = SQRT(AKW)   
      CLAQ = ALF2 - DEL1

      RETURN



      END

















      SUBROUTINE CALCNIAQ (NO3I, HI, DELT)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION NO3I, HI, DELT



      A42  = XK42*WATER/(GAMA(10))**2. 



      OM1  = NO3I          
      OM2  = HI
      BB   =-(OM1+OM2+A42)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)



      IF (DEL1.LT.ZERO .OR. DEL1.GT.HI .OR. DEL1.GT.NO3I) THEN
         print *, DELT
         DELT = ZERO
      ELSE
         DELT = DEL1
         RETURN
      ENDIF

      IF (DEL2.LT.ZERO .OR. DEL2.GT.NO3I .OR. DEL2.GT.HI) THEN
         DELT = ZERO
      ELSE
         DELT = DEL2
      ENDIF

      RETURN



      END

















      SUBROUTINE CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION NO3I, NO3AQ



      A42  = XK42*WATER/(GAMA(10))**2. 
      AKW  = XKW *RH*WATER*WATER



      ALF1  = NO3I - GGNO3
      ALF2  = GGNO3
      ALF3  = HI

      BB    = ALF3 + ALF1 + A42
      CC    = ALF3*ALF1 - A42*ALF2
      DEL1  = 0.5*(-BB + SQRT(BB*BB-4.D0*CC))



      NO3I  = ALF1 + DEL1
      HI    = ALF3 + DEL1
      IF (HI.LE.TINY) HI = SQRT(AKW)   
      NO3AQ = ALF2 - DEL1

      RETURN



      END

















      SUBROUTINE CALCMR
      USE SOLUT
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)




      CHARACTER SC*1



      SC =SCASE(1:1)                   



      IF (SC.EQ.'A') THEN      
         MOLALR(4) = MOLAL(5)+MOLAL(6) 



      ELSE IF (SC.EQ.'B') THEN
         SO4I  = MOLAL(5)-MOLAL(1)     
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF (SO4I.LT.HSO4I) THEN                
            MOLALR(13) = SO4I                   
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)  
         ELSE                                   
            MOLALR(13) = HSO4I                  
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)  
         ENDIF



      ELSE IF (SC.EQ.'C') THEN
         MOLALR(9) = MOLAL(3)                     
         MOLALR(7) = MAX(W(2)-W(3), ZERO)         



      ELSE IF (SC.EQ.'D') THEN      
         MOLALR(4) = MOLAL(5) + MOLAL(6)          
         AML5      = MOLAL(3)-2.D0*MOLALR(4)      
         MOLALR(5) = MAX(MIN(AML5,MOLAL(7)), ZERO)



      ELSE IF (SC.EQ.'E') THEN      
         SO4I  = MAX(MOLAL(5)-MOLAL(1),ZERO)      
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF (SO4I.LT.HSO4I) THEN                
            MOLALR(13) = SO4I                     
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)    
         ELSE                                   
            MOLALR(13) = HSO4I                    
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)    
         ENDIF



      ELSE IF (SC.EQ.'F') THEN      
         MOLALR(9) = MOLAL(3)                              
         MOLALR(7) = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3),ZERO)  



      ELSE IF (SC.EQ.'G') THEN      
         MOLALR(2) = 0.5*MOLAL(2)                          
         TOTS4     = MOLAL(5)+MOLAL(6)                     
         MOLALR(4) = MAX(TOTS4 - MOLALR(2), ZERO)          
         FRNH4     = MAX(MOLAL(3) - 2.D0*MOLALR(4), ZERO)
         MOLALR(5) = MIN(MOLAL(7),FRNH4)                   
         FRNH4     = MAX(FRNH4 - MOLALR(5), ZERO)
         MOLALR(6) = MIN(MOLAL(4), FRNH4)                  




      ELSE IF (SC.EQ.'H') THEN      
         MOLALR(1) = PSI7                                  
         MOLALR(2) = PSI1                                  
         MOLALR(3) = PSI8                                  
         MOLALR(4) = ZERO                                  
         FRNO3     = MAX(MOLAL(7) - MOLALR(3), ZERO)       
         FRCL      = MAX(MOLAL(4) - MOLALR(1), ZERO)       
         MOLALR(5) = MIN(MOLAL(3),FRNO3)                   
         FRNH4     = MAX(MOLAL(3) - MOLALR(5), ZERO)       
         MOLALR(6) = MIN(FRCL, FRNH4)                      




      ELSE IF (SC.EQ.'I') THEN      
         MOLALR(04) = PSI5                                 
         MOLALR(02) = PSI4                                 
         MOLALR(09) = PSI1                                 
         MOLALR(12) = PSI3                                 
         MOLALR(13) = PSI2                                 



      ELSE IF (SC.EQ.'J') THEN      
         MOLALR(09) = MOLAL(3)                             
         MOLALR(12) = MOLAL(2)                             
         MOLALR(07) = MOLAL(5)+MOLAL(6)-MOLAL(3)-MOLAL(2)  
         MOLALR(07) = MAX(MOLALR(07),ZERO)





      ELSE IF (SC.EQ.'N') THEN      
         MOLALR(4) = MOLAL(5) + MOLAL(6)          
         AML5      = WAER(3)-2.D0*MOLALR(4)       
         MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO) 



      ELSE IF (SC.EQ.'Q') THEN      
         MOLALR(2) = PSI1                                  
         MOLALR(4) = PSI6                                  
         MOLALR(5) = PSI5                                  
         MOLALR(6) = PSI4                                  



      ELSE IF (SC.EQ.'R') THEN      
         MOLALR(1) = PSI3                                  
         MOLALR(2) = PSI1                                  
         MOLALR(3) = PSI2                                  
         MOLALR(4) = ZERO                                  
         MOLALR(5) = PSI5                                  
         MOLALR(6) = PSI4                                  



      ELSE
         CALL PUSHERR (1001, ' ') 
      ENDIF



      WATER = ZERO
      DO 10 I=1,NPAIR
         WATER = WATER + MOLALR(I)/M0(I)
10    CONTINUE
      WATER = MAX(WATER, TINY)

      RETURN



      END

















      SUBROUTINE CALCMDRH (RHI, RHDRY, RHLIQ, DRYCASE, LIQCASE)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL DRYCASE, LIQCASE



      IF (WFTYP.EQ.0) THEN
         WF = ONE
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
      ENDIF
      ONEMWF  = ONE - WF



      CALL DRYCASE
      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  

      CNH42SO = CNH42S4                  
      CNH4HSO = CNH4HS4
      CLCO    = CLC 
      CNH4N3O = CNH4NO3
      CNH4CLO = CNH4CL
      CNA2SO  = CNA2SO4
      CNAHSO  = CNAHSO4
      CNANO   = CNANO3
      CNACLO  = CNACL
      GNH3O   = GNH3
      GHNO3O  = GHNO3
      GHCLO   = GHCL



      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CLC     = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CNAHSO4 = ZERO
      CNANO3  = ZERO
      CNACL   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
      CALL LIQCASE                   



      IF (WATER.LE.TINY) THEN
         DO 100 I=1,NIONS
            MOLAL(I)= ZERO           
  100    CONTINUE
         WATER   = ZERO

         CNH42S4 = CNH42SO           
         CNA2SO4 = CNA2SO
         CNAHSO4 = CNAHSO
         CNH4HS4 = CNH4HSO
         CLC     = CLCO
         CNH4NO3 = CNH4N3O
         CNANO3  = CNANO
         CNACL   = CNACLO                                                  
         CNH4CL  = CNH4CLO 

         GNH3    = GNH3O             
         GHNO3   = GHNO3O
         GHCL    = GHCLO

         GOTO 200
      ENDIF



      DAMSUL  = CNH42SO - CNH42S4
      DSOSUL  = CNA2SO  - CNA2SO4
      DAMBIS  = CNH4HSO - CNH4HS4
      DSOBIS  = CNAHSO  - CNAHSO4
      DLC     = CLCO    - CLC
      DAMNIT  = CNH4N3O - CNH4NO3
      DAMCHL  = CNH4CLO - CNH4CL
      DSONIT  = CNANO   - CNANO3
      DSOCHL  = CNACLO  - CNACL



      DAMG    = GNH3O   - GNH3 
      DHAG    = GHCLO   - GHCL
      DNAG    = GHNO3O  - GHNO3





      MOLAL(1)= ONEMWF*MOLAL(1)                                 
      MOLAL(2)= ONEMWF*(2.D0*DSOSUL + DSOBIS + DSONIT + DSOCHL) 
      MOLAL(3)= ONEMWF*(2.D0*DAMSUL + DAMG   + DAMBIS + DAMCHL +   &
                        3.D0*DLC    + DAMNIT )                  
      MOLAL(4)= ONEMWF*(     DAMCHL + DSOCHL + DHAG)            
      MOLAL(5)= ONEMWF*(     DAMSUL + DSOSUL + DLC - MOLAL(6))  
      MOLAL(6)= ONEMWF*(   MOLAL(6) + DSOBIS + DAMBIS + DLC)    
      MOLAL(7)= ONEMWF*(     DAMNIT + DSONIT + DNAG)            
      WATER   = ONEMWF*WATER



      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
      CNACL   = WF*CNACLO  + ONEMWF*CNACL
      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL



      GNH3    = WF*GNH3O   + ONEMWF*GNH3
      GHNO3   = WF*GHNO3O  + ONEMWF*GHNO3
      GHCL    = WF*GHCLO   + ONEMWF*GHCL



200   RETURN



      END























      SUBROUTINE CALCMDRP (RHI, RHDRY, RHLIQ, DRYCASE, LIQCASE)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      EXTERNAL DRYCASE, LIQCASE



      IF (WFTYP.EQ.0) THEN
         WF = ONE
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
      ENDIF
      ONEMWF  = ONE - WF



      CALL DRYCASE
      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  

      CNH42SO = CNH42S4              
      CNH4HSO = CNH4HS4
      CLCO    = CLC 
      CNH4N3O = CNH4NO3
      CNH4CLO = CNH4CL
      CNA2SO  = CNA2SO4
      CNAHSO  = CNAHSO4
      CNANO   = CNANO3
      CNACLO  = CNACL



      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CLC     = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CNAHSO4 = ZERO
      CNANO3  = ZERO
      CNACL   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
      CALL LIQCASE                   



      IF (WATER.LE.TINY) THEN
         WATER = ZERO
         DO 100 I=1,NIONS
            MOLAL(I)= ZERO
 100     CONTINUE
         CALL DRYCASE
         GOTO 200
      ENDIF



      DAMBIS  = CNH4HSO - CNH4HS4
      DSOBIS  = CNAHSO  - CNAHSO4
      DLC     = CLCO    - CLC





      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
      CNACL   = WF*CNACLO  + ONEMWF*CNACL
      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL



      WATER   = ONEMWF*WATER

      MOLAL(2)= WAER(1) - 2.D0*CNA2SO4 - CNAHSO4 - CNANO3 -        &
                               CNACL                            
      MOLAL(3)= WAER(3) - 2.D0*CNH42S4 - CNH4HS4 - CNH4CL -    &
                          3.D0*CLC     - CNH4NO3                
      MOLAL(4)= WAER(5) - CNACL - CNH4CL                        
      MOLAL(7)= WAER(4) - CNANO3 - CNH4NO3                      
      MOLAL(6)= ONEMWF*(MOLAL(6) + DSOBIS + DAMBIS + DLC)       
      MOLAL(5)= WAER(2) - MOLAL(6) - CLC - CNH42S4 - CNA2SO4    

      A8      = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      IF (MOLAL(5).LE.TINY) THEN
         HIEQ = SQRT(XKW *RH*WATER*WATER)  
      ELSE
         HIEQ = A8*MOLAL(6)/MOLAL(5)          
      ENDIF
      HIEN    = MOLAL(4) + MOLAL(7) + MOLAL(6) + 2.D0*MOLAL(5) -   &
                MOLAL(2) - MOLAL(3)
      MOLAL(1)= MAX (HIEQ, HIEN)                                



      A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. 
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        

      GNH3    = MOLAL(3)/MAX(MOLAL(1),TINY)/A2
      GHNO3   = MOLAL(1)*MOLAL(7)/A3
      GHCL    = MOLAL(1)*MOLAL(4)/A4

200   RETURN



      END













      SUBROUTINE CALCHS4 (HI, SO4I, HSO4I, DELTA)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)





      IF (WATER.LE.1d1*TINY) THEN
         DELTA = ZERO 
         RETURN
      ENDIF



      A8 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.

      BB =-(HI + SO4I + A8)
      CC = HI*SO4I - HSO4I*A8
      DD = BB*BB - 4.D0*CC

      IF (DD.GE.ZERO) THEN
         SQDD   = SQRT(DD)
         DELTA1 = 0.5*(-BB + SQDD)
         DELTA2 = 0.5*(-BB - SQDD)
         IF (HSO4I.LE.TINY) THEN
            DELTA = DELTA2
         ELSEIF( HI*SO4I .GE. A8*HSO4I ) THEN
            DELTA = DELTA2
         ELSEIF( HI*SO4I .LT. A8*HSO4I ) THEN
            DELTA = DELTA1
         ELSE
            DELTA = ZERO
         ENDIF
      ELSE
         DELTA  = ZERO
      ENDIF











      RETURN



      END














      SUBROUTINE CALCPH (GG, HI, OHI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      AKW  = XKW *RH*WATER*WATER
      CN   = SQRT(AKW)



      IF (GG.GT.TINY) THEN                        
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = MAX(0.5D0*(-BB + SQRT(DD)),CN)
         OHI= AKW/HI
      ELSE                                        
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= MAX(0.5D0*(-BB + SQRT(DD)),CN)
         HI = AKW/OHI
      ENDIF

      RETURN



      END
















      SUBROUTINE CALCACT
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      REAL EX10, URF
      REAL G0(3,4),ZPL,ZMI,AGAMA,SION,H,CH,F1(3),F2(4)
      DOUBLE PRECISION MPL, XIJ, YJI

      PARAMETER (LN10=2.30258509299404568402D0)

      G(I,J)= (F1(I)/Z(I) + F2(J)/Z(J+3)) / (Z(I)+Z(J+3)) - H



      IF (FRST) THEN               
         DO 10 I=1,NPAIR
            GAMOU(I) = GAMA(I)
10       CONTINUE
      ENDIF

      DO 20 I=1,NPAIR              
         GAMIN(I) = GAMA(I)
20    CONTINUE



      IONIC=0.0
      DO 30 I=1,NIONS
         IONIC=IONIC + MOLAL(I)*Z(I)*Z(I)
30    CONTINUE
      IONIC = MAX(MIN(0.5*IONIC/WATER,500.d0), TINY)






      IF (IACALC.EQ.0) THEN              
         CALL KMFUL (IONIC, SNGL(TEMP),G0(2,1),G0(2,2),G0(2,4),   &
                     G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),   &
                     G0(1,4),G0(1,1),G0(2,3))
      ELSE                               
         CALL KMTAB (IONIC, SNGL(TEMP),G0(2,1),G0(2,2),G0(2,4),   &
                     G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),   &
                     G0(1,4),G0(1,1),G0(2,3))
      ENDIF



      AGAMA = 0.511*(298.0/TEMP)**1.5    
      SION  = SQRT(IONIC)
      H     = AGAMA*SION/(1+SION)

      DO 100 I=1,3
         F1(I)=0.0
         F2(I)=0.0
100   CONTINUE
      F2(4)=0.0

      DO 110 I=1,3
         ZPL = Z(I)
         MPL = MOLAL(I)/WATER
         DO 110 J=1,4
            ZMI   = Z(J+3)
            CH    = 0.25*(ZPL+ZMI)*(ZPL+ZMI)/IONIC
            XIJ   = CH*MPL
            YJI   = CH*MOLAL(J+3)/WATER
            F1(I) = F1(I) + SNGL(YJI*(G0(I,J) + ZPL*ZMI*H))
            F2(J) = F2(J) + SNGL(XIJ*(G0(I,J) + ZPL*ZMI*H))
110   CONTINUE



      GAMA(01) = G(2,1)*ZZ(01)                     
      GAMA(02) = G(2,2)*ZZ(02)                     
      GAMA(03) = G(2,4)*ZZ(03)                     
      GAMA(04) = G(3,2)*ZZ(04)                     
      GAMA(05) = G(3,4)*ZZ(05)                     
      GAMA(06) = G(3,1)*ZZ(06)                     
      GAMA(07) = G(1,2)*ZZ(07)                     
      GAMA(08) = G(1,3)*ZZ(08)                     
      GAMA(09) = G(3,3)*ZZ(09)                     
      GAMA(10) = G(1,4)*ZZ(10)                     
      GAMA(11) = G(1,1)*ZZ(11)                     
      GAMA(12) = G(2,3)*ZZ(12)                     
      GAMA(13) = 0.20*(3.0*GAMA(04)+2.0*GAMA(09))  





      DO 200 I=1,NPAIR
         GAMA(I)=MAX(-11.0d0, MIN(GAMA(I),11.0d0) ) 

         GAMA(I)=EXP(LN10*GAMA(I))


  200 CONTINUE





      IF (FRST) THEN          
         ERROU = ZERO                    
         DO 210 I=1,NPAIR
            ERROU=MAX(ERROU, ABS((GAMOU(I)-GAMA(I))/GAMOU(I)))
210      CONTINUE
         CALAOU = ERROU .GE. EPSACT      
         FRST   =.FALSE.
      ENDIF



      ERRIN = ZERO                       
      DO 220 I=1,NPAIR
         ERRIN = MAX (ERRIN, ABS((GAMIN(I)-GAMA(I))/GAMIN(I)))
220   CONTINUE
      CALAIN = ERRIN .GE. EPSACT

      ICLACT = ICLACT + 1                



      RETURN
      END















      SUBROUTINE RSTGAM
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)


      DO 10 I=1, NPAIR
         GAMA(I) = 0.1
10    CONTINUE



      RETURN
      END      













      SUBROUTINE KMFUL (IONIC,TEMP,G01,G02,G03,G04,G05,G06,G07,G08,G09,   &
                        G10,G11,G12)
      REAL Ionic, TEMP
      DATA Z01,Z02,Z03,Z04,Z05,Z06,Z07,Z08,Z10,Z11   &
          /1,  2,  1,  2,  1,  1,  2,  1,  1,  1/

      SION = SQRT(IONIC)



      CALL MKBI(2.230, IONIC, SION, Z01, G01)
      CALL MKBI(-0.19, IONIC, SION, Z02, G02)
      CALL MKBI(-0.39, IONIC, SION, Z03, G03)
      CALL MKBI(-0.25, IONIC, SION, Z04, G04)
      CALL MKBI(-1.15, IONIC, SION, Z05, G05)
      CALL MKBI(0.820, IONIC, SION, Z06, G06)
      CALL MKBI(-.100, IONIC, SION, Z07, G07)
      CALL MKBI(8.000, IONIC, SION, Z08, G08)
      CALL MKBI(2.600, IONIC, SION, Z10, G10)
      CALL MKBI(6.000, IONIC, SION, Z11, G11)



      TI  = TEMP-273.0
      TC  = TI-25.0
      IF (ABS(TC) .GT. 1.0) THEN
         CF1 = 1.125-0.005*TI
         CF2 = (0.125-0.005*TI)*(0.039*IONIC**0.92-0.41*SION/(1.+SION))
         G01 = CF1*G01 - CF2*Z01
         G02 = CF1*G02 - CF2*Z02
         G03 = CF1*G03 - CF2*Z03
         G04 = CF1*G04 - CF2*Z04
         G05 = CF1*G05 - CF2*Z05
         G06 = CF1*G06 - CF2*Z06
         G07 = CF1*G07 - CF2*Z07
         G08 = CF1*G08 - CF2*Z08
         G10 = CF1*G10 - CF2*Z10
         G11 = CF1*G11 - CF2*Z11
      ENDIF

      G09 = G06 + G08 - G11
      G12 = G01 + G08 - G11



      RETURN
      END















      SUBROUTINE MKBI(Q,IONIC,SION,ZIP,BI)

      REAL IONIC

      B=.75-.065*Q
      C= 1.0
      IF (IONIC.LT.6.0) C=1.+.055*Q*EXP(-.023*IONIC*IONIC*IONIC)
      XX=-0.5107*SION/(1.+C*SION)
      BI=(1.+B*(1.+.1*IONIC)**Q-B)
      BI=ZIP*ALOG10(BI) + ZIP*XX

      RETURN
      END
















      SUBROUTINE KMTAB (IN,TEMP,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,   &
                                G11,G12)
      REAL IN, Temp



      IND = NINT((TEMP-198.0)/25.0) + 1
      IND = MIN(MAX(IND,1),6)



      IF (IND.EQ.1) THEN
         CALL KM198 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSEIF (IND.EQ.2) THEN
         CALL KM223 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSEIF (IND.EQ.3) THEN
         CALL KM248 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSEIF (IND.EQ.4) THEN
         CALL KM273 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSEIF (IND.EQ.5) THEN
         CALL KM298 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSE
         CALL KM323 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ENDIF



      RETURN
      END


      INTEGER FUNCTION IBACPOS(IN)






      implicit none
      real IN
      IF (IN .LE. 0.300000E+02) THEN
         ibacpos = MIN(NINT( 0.200000E+02*IN) + 1, 600)
      ELSE
         ibacpos =   600+NINT( 0.200000E+01*IN- 0.600000E+02)
      ENDIF
      ibacpos = min(ibacpos, 741)
      return
      end



















      SUBROUTINE KM198 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,   &
                           G11,G12)
       USE KMC198








      REAL IN



      ipos = ibacpos(IN)



      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)



      RETURN
      END




















      SUBROUTINE KM223 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,   &
                           G11,G12)
      USE KMC223








      REAL IN



      ipos = ibacpos(IN)



      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)



      RETURN
      END





















      SUBROUTINE KM248 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,   &
                           G11,G12)
       USE KMC248








      REAL IN



      ipos = ibacpos(IN)



      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)



      RETURN
      END




















      SUBROUTINE KM273 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,   &
                           G11,G12)
      USE KMC273








      REAL IN



      ipos = ibacpos(IN)



      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)



      RETURN
      END




















      SUBROUTINE KM298 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,   &
                           G11,G12)
      USE KMC298








      REAL IN



      ipos = ibacpos(IN)



      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)



      RETURN
      END




















      SUBROUTINE KM323 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,   &
                           G11,G12)
      USE KMC323








      REAL IN



      ipos = ibacpos(IN)



      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)



      RETURN
      END





      BLOCK DATA KMCF198
       USE KMC198









      END



      BLOCK DATA KMCF223
      USE KMC223









      END



      BLOCK DATA KMCF248
       USE KMC248








      END



      BLOCK DATA KMCF273
      USE KMC273








      END



      BLOCK DATA KMCF298
      USE KMC298








      END



      BLOCK DATA KMCF323
      USE KMC323








      END

























      SUBROUTINE CHRBLN (STR, IBLK)


      CHARACTER*(*) STR

      IBLK = 1                       
      ILEN = LEN(STR)                
      DO 10 i=ILEN,1,-1
         IF (STR(i:i).NE.' ' .AND. STR(i:i).NE.CHAR(0)) THEN
            IBLK = i
            RETURN
         ENDIF
10    CONTINUE
      RETURN

      END























      SUBROUTINE SHFTRGHT (CHR)


      CHARACTER CHR*(*)

      I1  = LEN(CHR)             
      CALL CHRBLN(CHR,I2)        
      IF (I2.EQ.I1) RETURN

      DO 10 I=I2,1,-1            
         CHR(I1+I-I2:I1+I-I2) = CHR(I:I)
         CHR(I:I) = ' '
10    CONTINUE
      RETURN

      END































      SUBROUTINE RPLSTR (STRING, OLD, NEW, IERR)


      CHARACTER STRING*(*), OLD*(*), NEW*(*)



      ILO = LEN(OLD)



      IP = INDEX(NEW,OLD)
      IF (IP.NE.0) THEN
         IERR = 1
         RETURN
      ELSE
         IERR = 0
      ENDIF



10    IP = INDEX(STRING,OLD)      
      IF (IP.EQ.0) RETURN         
      STRING(IP:IP+ILO-1) = NEW   
      GOTO 10                     

      END
        

































      SUBROUTINE INPTD (VAR, DEF, PROMPT, PRFMT, IERR)


      CHARACTER PROMPT*(*), PRFMT*(*), BUFFER*128
      DOUBLE PRECISION DEF, VAR
      INTEGER IERR

      IERR = 0



      WRITE (BUFFER, FMT=PRFMT, ERR=10) DEF
      CALL CHRBLN (BUFFER, IEND)



      WRITE (*,*) PROMPT,' [',BUFFER(1:IEND),']: '
      READ  (*, '(A)', ERR=20, END=20) BUFFER
      CALL CHRBLN (BUFFER,IEND)



      IF (IEND.EQ.1 .AND. BUFFER(1:1).EQ.' ') THEN
         VAR = DEF
      ELSE
         READ (BUFFER, *, ERR=20, END=20) VAR
      ENDIF



30    RETURN



10    IERR = 1       
      GOTO 30

20    IERR = 2       
      GOTO 30

      END
























      SUBROUTINE Pushend (Iunit)



      LOGICAL OPNED



      INQUIRE (UNIT=Iunit, OPENED=OPNED)
      IF (.NOT.OPNED) GOTO 25



10    READ (Iunit,'()', ERR=20, END=20)
      GOTO 10



20    BACKSPACE (Iunit)
25    RETURN
      END






























      SUBROUTINE Appendext (Filename, Defext, Overwrite)


      CHARACTER*(*) Filename, Defext
      LOGICAL       Overwrite

      CALL CHRBLN (Filename, Iend)
      IF (Filename(1:1).EQ.' ' .AND. Iend.EQ.1) RETURN  
      Idot = INDEX (Filename, '.')                      
      IF (Idot.EQ.0) Filename = Filename(1:Iend)//Defext
      IF (Overwrite .AND. Idot.NE.0)   &
                    Filename = Filename(:Idot-1)//Defext
      RETURN
      END





































      SUBROUTINE POLY3 (A1, A2, A3, ROOT, ISLV)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (EXPON=1.D0/3.D0,     ZERO=0.D0, THET1=120.D0/180.D0,    &
                 THET2=240.D0/180.D0, PI=3.14159265358932, EPS=1D-50)
      DOUBLE PRECISION  X(3)



      IF (ABS(A3).LE.EPS) THEN 
         ISLV = 1
         IX   = 1
         X(1) = ZERO
         D    = A1*A1-4.D0*A2
         IF (D.GE.ZERO) THEN
            IX   = 3
            SQD  = SQRT(D)
            X(2) = 0.5*(-A1+SQD)
            X(3) = 0.5*(-A1-SQD)
         ENDIF
      ELSE





         ISLV= 1
         Q   = (3.D0*A2 - A1*A1)/9.D0
         R   = (9.D0*A1*A2 - 27.D0*A3 - 2.D0*A1*A1*A1)/54.D0
         D   = Q*Q*Q + R*R





         IF (D.LT.-EPS) THEN        
            IX   = 3
            THET = EXPON*ACOS(R/SQRT(-Q*Q*Q))
            COEF = 2.D0*SQRT(-Q)
            X(1) = COEF*COS(THET)            - EXPON*A1
            X(2) = COEF*COS(THET + THET1*PI) - EXPON*A1
            X(3) = COEF*COS(THET + THET2*PI) - EXPON*A1



         ELSE IF (D.LE.EPS) THEN    
            IX   = 2
            SSIG = SIGN (1.D0, R)
            S    = SSIG*(ABS(R))**EXPON
            X(1) = 2.D0*S  - EXPON*A1
            X(2) =     -S  - EXPON*A1



         ELSE                       
            IX   = 1
            SQD  = SQRT(D)
            SSIG = SIGN (1.D0, R+SQD)       
            TSIG = SIGN (1.D0, R-SQD)
            S    = SSIG*(ABS(R+SQD))**EXPON 
            T    = TSIG*(ABS(R-SQD))**EXPON
            X(1) = S + T - EXPON*A1
         ENDIF
      ENDIF



      ROOT = 1.D30
      DO 10 I=1,IX
         IF (X(I).GT.ZERO) THEN
            ROOT = MIN (ROOT, X(I))
            ISLV = 0
         ENDIF
10    CONTINUE



      RETURN
      END





















      SUBROUTINE POLY3B (A1, A2, A3, RTLW, RTHI, ROOT, ISLV)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ZERO=0.D0, EPS=1D-15, MAXIT=100, NDIV=5)

      FUNC(X) = X**3.d0 + A1*X**2.0 + A2*X + A3



      X1   = RTLW
      Y1   = FUNC(X1)
      IF (ABS(Y1).LE.EPS) THEN     
         ROOT = RTLW
         GOTO 50
      ENDIF



      DX = (RTHI-RTLW)/REAL(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNC (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO) GOTO 20 
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .LT. EPS) THEN   
         ROOT = X2
      ELSE
         ROOT = 1.d30
         ISLV = 1
      ENDIF
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNC (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE



40    X3   = 0.5*(X1+X2)
      Y3   = FUNC (X3)
      ROOT = X3
      ISLV = 0

50    RETURN



      END
      







































      FUNCTION EX10(X,K)
      USE EXPNC
      REAL    X, EX10, Y , K
      INTEGER K1, K2




      Y    = MAX(-K, MIN(X,K))   



      K1   = INT(Y)
      K2   = INT(100*(Y-K1))



      EX10 = AINT10(K1+10)*ADEC10(K2+100)



      RETURN
      END















      BLOCK DATA EXPON



      USE EXPNC





























































      END

























































































      SUBROUTINE ISOPLUS (WI,  RHI,    TEMPI,  IPROBI,    &
                          GAS, AERLIQ, AERSLD, DRYI   )
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP), GAS(NGASAQ), AERLIQ(NIONS+NGASAQ+1),   &
                AERSLD(NSLDS)
      LOGICAL   DRYI



      IPROB = IPROBI



      IF (IPROB.EQ.0) THEN
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN 
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        
            CALL ISRP1F (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              
            CALL ISRP2F (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3F (WI, RHI, TEMPI)
         ENDIF



      ELSE
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN 
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        
            CALL ISRP1R (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              
            CALL ISRP2R (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3R (WI, RHI, TEMPI)
         ENDIF
      ENDIF



      GAS(1) = GNH3
      GAS(2) = GHNO3
      GAS(3) = GHCL

      DO 10 I=1,NIONS
         AERLIQ(I) = MOLAL(I)
  10  CONTINUE
      AERLIQ(NIONS+1) = WATER*1.0D3/18.0D0
      DO 20 I=1,NGASAQ
         AERLIQ(NIONS+1+I) = GASAQ(I)
  20  CONTINUE

      AERSLD(1) = CNANO3
      AERSLD(2) = CNH4NO3
      AERSLD(3) = CNACL
      AERSLD(4) = CNH4CL
      AERSLD(5) = CNA2SO4
      AERSLD(6) = CNH42S4
      AERSLD(7) = CNAHSO4
      AERSLD(8) = CNH4HS4
      AERSLD(9) = CLC

      DRYI = WATER.LE.TINY

      RETURN



      END








































      SUBROUTINE ISRPIAA (WI, RHI, TEMPI, IPROBI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WI(NCOMP)



      IPROB = IPROBI



      IF (IPROB.EQ.0) THEN
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN 
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        
            CALL ISRP1F (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              
            CALL ISRP2F (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3F (WI, RHI, TEMPI)
         ENDIF



      ELSE
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN 
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        
            CALL ISRP1R (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              
            CALL ISRP2R (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3R (WI, RHI, TEMPI)
         ENDIF
      ENDIF



      DRYF = WATER.LE.TINY



      IF (IPROB.EQ.0) THEN
         CONTINUE
      ELSE
         W(1) = WAER(1)
         W(2) = WAER(2)
         W(3) = WAER(3) 
         W(4) = WAER(4)
         W(5) = WAER(5)

         IF (.NOT.DRYF) THEN
            W(3) = W(3) + GNH3 
            W(4) = W(4) + GHNO3
            W(5) = W(5) + GHCL
         ENDIF
      ENDIF

      RETURN



      END













      SUBROUTINE PUSHERR (IERR,ERRINF)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      CHARACTER ERRINF*(*) 



      IF (NOFER.LT.NERRMX) THEN   
         NOFER         = NOFER + 1 
         ERRSTK(NOFER) = IERR
         ERRMSG(NOFER) = ERRINF   
         STKOFL        =.FALSE.
      ELSE
         STKOFL        =.TRUE.      
      ENDIF



      END
      















      SUBROUTINE ISERRINF (ERRSTKI, ERRMSGI, NOFERI, STKOFLI)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      CHARACTER ERRMSGI*40
      INTEGER   ERRSTKI
      LOGICAL   STKOFLI
      DIMENSION ERRMSGI(NERRMX), ERRSTKI(NERRMX)



      DO 10 I=1,NOFER              
        ERRSTKI(I) = ERRSTK(I)
        ERRMSGI(I) = ERRMSG(I)
  10  CONTINUE

      STKOFLI = STKOFL
      NOFERI  = NOFER

      RETURN



      END
      















      SUBROUTINE ERRSTAT (IO,IERR,ERRINF)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      CHARACTER CER*4, NCIS*29, NCIF*27, NSIS*26, NSIF*24, ERRINF*(*)
      DATA NCIS /'NO CONVERGENCE IN SUBROUTINE '/,   &
           NCIF /'NO CONVERGENCE IN FUNCTION '  /,   &
           NSIS /'NO SOLUTION IN SUBROUTINE '   /,   &
           NSIF /'NO SOLUTION IN FUNCTION '     /



      WRITE (CER,'(I4)') IERR
      CALL RPLSTR (CER, ' ', '0',IOK)   
      CALL CHRBLN (ERRINF, IEND)        



      IF (IERR.EQ.0) THEN
         WRITE (IO,1000) 'NO ERRORS DETECTED '
         GOTO 10

      ELSE IF (IERR.LT.0) THEN
         WRITE (IO,1000) 'ERROR STACK EXHAUSTED '
         GOTO 10

      ELSE IF (IERR.GT.1000) THEN
         WRITE (IO,1100) 'FATAL',CER

      ELSE
         WRITE (IO,1100) 'WARNING',CER
      ENDIF





      IF (IERR.EQ.1001) THEN 
         CALL CHRBLN (SCASE, IEND)
         WRITE (IO,1000) 'CASE NOT SUPPORTED IN CALCMR ['//SCASE(1:IEND)   &
                         //']'

      ELSEIF (IERR.EQ.1002) THEN 
         CALL CHRBLN (SCASE, IEND)
         WRITE (IO,1000) 'CASE NOT SUPPORTED ['//SCASE(1:IEND)//']'



      ELSEIF (IERR.EQ.0001) THEN 
         WRITE (IO,1000) NSIS,ERRINF

      ELSEIF (IERR.EQ.0002) THEN 
         WRITE (IO,1000) NCIS,ERRINF

      ELSEIF (IERR.EQ.0003) THEN 
         WRITE (IO,1000) NSIF,ERRINF

      ELSEIF (IERR.EQ.0004) THEN 
         WRITE (IO,1000) NCIF,ERRINF

      ELSE IF (IERR.EQ.0019) THEN
         WRITE (IO,1000) 'HNO3(aq) AFFECTS H+, WHICH '//   &
                         'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE (IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF (IERR.EQ.0020) THEN
         IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN
            WRITE (IO,1000) 'HSO4-SO4 EQUILIBRIUM MIGHT AFFECT HNO3,'   &
                          //'HCL DISSOLUTION'
         ELSE
            WRITE (IO,1000) 'HSO4-SO4 EQUILIBRIUM MIGHT AFFECT NH3 '   &
                          //'DISSOLUTION'
         ENDIF
         WRITE (IO,1000) 'DIRECT DECREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF (IERR.EQ.0021) THEN
         WRITE (IO,1000) 'HNO3(aq),HCL(aq) AFFECT H+, WHICH '//   &
                         'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE (IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF (IERR.EQ.0022) THEN
         WRITE (IO,1000) 'HCL(g) EQUILIBRIUM YIELDS NONPHYSICAL '//   &
                         'DISSOLUTION'
         WRITE (IO,1000) 'A TINY AMOUNT [',ERRINF(1:IEND),'] IS '//   &
                         'ASSUMED TO BE DISSOLVED'

      ELSEIF (IERR.EQ.0033) THEN
         WRITE (IO,1000) 'HCL(aq) AFFECTS H+, WHICH '//   &
                         'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE (IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSEIF (IERR.EQ.0050) THEN
         WRITE (IO,1000) 'TOO MUCH SODIUM GIVEN AS INPUT.'
         WRITE (IO,1000) 'REDUCED TO COMPLETELY NEUTRALIZE SO4,Cl,NO3.'
         WRITE (IO,1000) 'EXCESS SODIUM IS IGNORED.'

      ELSE
         WRITE (IO,1000) 'NO DIAGNOSTIC MESSAGE AVAILABLE'
      ENDIF

10    RETURN



1000  FORMAT (1X,A:A:A:A:A)
1100  FORMAT (1X,A,' ERROR [',A4,']:')



      END



















































      SUBROUTINE ISORINF (VERSI, NCMP, NION, NAQGAS, NSOL, NERR, TIN,   &
                          GRT)
      USE ISRPIA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      CHARACTER VERSI*(*)



      VERSI  = VERSION
      NCMP   = NCOMP
      NION   = NIONS
      NAQGAS = NGASAQ
      NSOL   = NSLDS
      NERR   = NERRMX
      TIN    = TINY
      GRT    = GREAT

      RETURN



      END
