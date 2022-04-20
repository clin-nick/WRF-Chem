

























































































































MODULE module_mp_fer_hires








      REAL,PRIVATE,SAVE ::  ABFR, CBFR, CIACW, CIACR, C_N0r0,           &
     &  C_NR, CRAIN, &    
     &  ARAUT, BRAUT, CN0r0, CN0r_DMRmin, CN0r_DMRmax, CRACW, ESW0,     &
     &  RFmax, RQR_DRmin, RQR_DRmax, RR_DRmin, RR_DR1, RR_DR2,          &
     &  RR_DR3, RR_DR4, RR_DR5, RR_DRmax, BETA6, PI_E, RFmx1, ARcw,     & 
     &  RH_NgC, RH_NgT, RQhail, AVhail, BVhail, QAUT0                

      INTEGER,PRIVATE,PARAMETER :: INDEXRstrmax=500      
      INTEGER, PRIVATE,PARAMETER :: MY_T1=1, MY_T2=35
      REAL,PRIVATE,DIMENSION(MY_T1:MY_T2),SAVE :: MY_GROWTH

      REAL, PRIVATE,PARAMETER :: DMImin=.05e-3, DMImax=1.e-3,           &
     &      DelDMI=1.e-6,XMImin=1.e6*DMImin
      REAL, PUBLIC,PARAMETER :: XMImax=1.e6*DMImax, XMIexp=.0536
      INTEGER, PUBLIC,PARAMETER :: MDImin=XMImin, MDImax=XMImax
      REAL, PRIVATE,DIMENSION(MDImin:MDImax) ::                         &
     &      ACCRI,VSNOWI,VENTI1,VENTI2
      REAL, PUBLIC,DIMENSION(MDImin:MDImax) :: SDENS    

      REAL, PRIVATE,PARAMETER :: DMRmin=.05e-3, DMRmax=1.0e-3,           &
     &      DelDMR=1.e-6, XMRmin=1.e6*DMRmin, XMRmax=1.e6*DMRmax
      INTEGER, PUBLIC,PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax

      REAL, PRIVATE,DIMENSION(MDRmin:MDRmax)::                           &
     &      ACCRR,MASSR,RRATE,VRAIN,VENTR1,VENTR2

      INTEGER, PRIVATE,PARAMETER :: Nrime=40
      REAL, DIMENSION(2:9,0:Nrime),PRIVATE,SAVE :: VEL_RF

      INTEGER,PARAMETER :: NX=7501
      REAL, PARAMETER :: XMIN=180.0,XMAX=330.0
      REAL, DIMENSION(NX),PRIVATE,SAVE :: TBPVS,TBPVS0
      REAL, PRIVATE,SAVE :: C1XPVS0,C2XPVS0,C1XPVS,C2XPVS

      REAL, PRIVATE,PARAMETER ::                                        &

     &   CP=1004.6, EPSQ=1.E-12, GRAV=9.806, RHOL=1000., RD=287.04      &
     &  ,RV=461.5, T0C=273.15, XLV=2.5E6, XLF=3.3358e+5                 &
     &  ,pi=3.141592653589793                                           &

     &  ,EPS=RD/RV, EPS1=RV/RD-1., EPSQ1=1.001*EPSQ                     &
     &  ,RCP=1./CP, RCPRV=RCP/RV, RGRAV=1./GRAV, RRHOL=1./RHOL          &
     &  ,XLS=XLV+XLF, XLV1=XLV/CP, XLF1=XLF/CP, XLS1=XLS/CP             &
     &  ,XLV2=XLV*XLV/RV, XLS2=XLS*XLS/RV                               &
     &  ,XLV3=XLV*XLV*RCPRV, XLS3=XLS*XLS*RCPRV                         &


     &  ,CLIMIT=10.*EPSQ, CLIMIT1=-CLIMIT                               &
     &  ,C1=1./3.                                                       &
     &  ,DMR1=.1E-3, DMR2=.2E-3, DMR3=.32E-3, DMR4=0.45E-3              &
     &  ,DMR5=0.67E-3                                                   &
     &  ,XMR1=1.e6*DMR1, XMR2=1.e6*DMR2, XMR3=1.e6*DMR3                 &
     &  ,XMR4=1.e6*DMR4, XMR5=1.e6*DMR5, RQRmix=0.05E-3, RQSmix=1.E-3   & 
     &  ,Cdry=1.634e13, Cwet=1./.224       
      INTEGER, PARAMETER :: MDR1=XMR1, MDR2=XMR2, MDR3=XMR3, MDR4=XMR4  &
     &  , MDR5=XMR5


LOGICAL, SAVE :: WARN1=.TRUE.,WARN2=.TRUE.,WARN3=.TRUE.,WARN5=.TRUE.
REAL, SAVE :: Pwarn=75.E2, QTwarn=1.E-3
INTEGER, PARAMETER :: MAX_ITERATIONS=10























      REAL, PUBLIC,PARAMETER ::                                         &
     &  RHgrd_in=1.                                                     &  
     & ,RHgrd_out=0.975                                                 &  
     & ,P_RHgrd_out=850.E2                                              &  
     & ,T_ICE=-40.                                                      &
     & ,T_ICEK=T0C+T_ICE                                                &
     & ,T_ICE_init=-12.                                                 &
     & ,NSI_max=250.E3                                                  &
     & ,NLImin=1.E3                                                     &
     & ,N0r0=8.E6                                                       &
     & ,N0rmin=1.E4                                                     &


     & ,NCW=250.E6                                                         

      LOGICAL, PARAMETER :: PRINT_diag=.FALSE.  

      REAL, PUBLIC,DIMENSION(MDImin:MDImax) :: MASSI


      CONTAINS



      SUBROUTINE FER_HIRES (itimestep,DT,DX,DY,GID,RAINNC,RAINNCV,    & 
     &                      dz8w,rho_phy,p_phy,pi_phy,th_phy,qv,qt,   & 
     &                      LOWLYR,SR,                                &
     &                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY,         &
     &                      QC,QR,QI,                                 &
     &                      ids,ide, jds,jde, kds,kde,                &
     &                      ims,ime, jms,jme, kms,kme,                &
     &                      its,ite, jts,jte, kts,kte                 )










      IMPLICIT NONE


      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,ITIMESTEP,GID  

      REAL, INTENT(IN) 	    :: DT,DX,DY
      REAL, INTENT(IN),     DIMENSION(ims:ime, kms:kme, jms:jme)::      &
     &                      dz8w,p_phy,pi_phy,rho_phy
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme)::      &
     &                      th_phy,qv,qt,qc,qr,qi
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme ) ::    &
     &                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme)           ::     &
     &                                                   RAINNC,RAINNCV
      REAL, INTENT(OUT),    DIMENSION(ims:ime,jms:jme):: SR





      INTEGER, DIMENSION( ims:ime, jms:jme ),INTENT(INOUT) :: LOWLYR













      REAL,  DIMENSION( ims:ime, kms:kme, jms:jme ) ::                  &
     &       TLATGS_PHY,TRAIN_PHY
      REAL,  DIMENSION(ims:ime,jms:jme):: APREC,PREC,ACPREC
      REAL,  DIMENSION(its:ite, kts:kte, jts:jte):: t_phy

      INTEGER :: I,J,K,KFLIP
      REAL :: WC





















          CALL MY_GROWTH_RATES (DT)





        CIACW=DT*0.25*PI_E*0.5*(1.E5)**C1




        CIACR=PI_E*DT




        CRACW=DT*0.25*PI_E*1.0



        BRAUT=DT*1.1E10*BETA6/NCW











      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
        t_phy(i,k,j) = th_phy(i,k,j)*pi_phy(i,k,j)
        qv(i,k,j)=qv(i,k,j)/(1.+qv(i,k,j)) 
      ENDDO
      ENDDO
      ENDDO



      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
	 TLATGS_PHY (i,k,j)=0.
	 TRAIN_PHY  (i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO

      DO j = jts,jte
      DO i = its,ite
         ACPREC(i,j)=0.
         APREC (i,j)=0.
         PREC  (i,j)=0.
         SR    (i,j)=0.
      ENDDO
      ENDDO


      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
         QT(I,K,J)=QC(I,K,J)+QR(I,K,J)+QI(I,K,J)
         IF (QI(I,K,J) <= EPSQ) THEN
            F_ICE_PHY(I,K,J)=0.
            F_RIMEF_PHY(I,K,J)=1.
            IF (T_PHY(I,K,J) < T_ICEK) F_ICE_PHY(I,K,J)=1.
         ELSE
            F_ICE_PHY(I,K,J)=MAX( 0., MIN(1., QI(I,K,J)/QT(I,K,J) ) )
         ENDIF
         IF (QR(I,K,J) <= EPSQ) THEN
            F_RAIN_PHY(I,K,J)=0.
         ELSE
            F_RAIN_PHY(I,K,J)=QR(I,K,J)/(QR(I,K,J)+QC(I,K,J))
         ENDIF
      ENDDO
      ENDDO
      ENDDO



      CALL EGCP01DRV(GID,DT,LOWLYR,                                     &
     &               APREC,PREC,ACPREC,SR,                              &
     &               dz8w,rho_phy,qt,t_phy,qv,F_ICE_PHY,P_PHY,          &
     &               F_RAIN_PHY,F_RIMEF_PHY,TLATGS_PHY,TRAIN_PHY,       &
     &               ids,ide, jds,jde, kds,kde,		                &
     &               ims,ime, jms,jme, kms,kme,		                &
     &               its,ite, jts,jte, kts,kte		          )


     DO j = jts,jte
        DO k = kts,kte
	DO i = its,ite
  	   th_phy(i,k,j) = t_phy(i,k,j)/pi_phy(i,k,j)
           qv(i,k,j)=qv(i,k,j)/(1.-qv(i,k,j))  
           WC=qt(I,K,J)
           QI(I,K,J)=0.
           QR(I,K,J)=0.
           QC(I,K,J)=0.
           IF(F_ICE_PHY(I,K,J)>=1.)THEN
             QI(I,K,J)=WC
           ELSEIF(F_ICE_PHY(I,K,J)<=0.)THEN
             QC(I,K,J)=WC
           ELSE
             QI(I,K,J)=F_ICE_PHY(I,K,J)*WC
             QC(I,K,J)=WC-QI(I,K,J)
           ENDIF

           IF(QC(I,K,J)>0..AND.F_RAIN_PHY(I,K,J)>0.)THEN
             IF(F_RAIN_PHY(I,K,J).GE.1.)THEN
               QR(I,K,J)=QC(I,K,J)
               QC(I,K,J)=0.
             ELSE
               QR(I,K,J)=F_RAIN_PHY(I,K,J)*QC(I,K,J)
               QC(I,K,J)=QC(I,K,J)-QR(I,K,J)
             ENDIF
          endif
	ENDDO
        ENDDO
     ENDDO



       DO j=jts,jte
       DO i=its,ite
          RAINNC(i,j)=APREC(i,j)*1000.+RAINNC(i,j)
          RAINNCV(i,j)=APREC(i,j)*1000.
       ENDDO
       ENDDO
















  END SUBROUTINE FER_HIRES







      SUBROUTINE FER_HIRES_ADVECT (itimestep,DT,DX,DY,GID,RAINNC,RAINNCV,    & 
     &                      dz8w,rho_phy,p_phy,pi_phy,th_phy,qv,   & 
     &                      LOWLYR,SR,                                &
     &                      QC,QR,QI,QRIMEF,                          &
     &                      ids,ide, jds,jde, kds,kde,                &
     &                      ims,ime, jms,jme, kms,kme,                &
     &                      its,ite, jts,jte, kts,kte                 )










      IMPLICIT NONE


      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,ITIMESTEP,GID  

      REAL, INTENT(IN) 	    :: DT,DX,DY
      REAL, INTENT(IN),     DIMENSION(ims:ime, kms:kme, jms:jme)::      &
     &                      dz8w,p_phy,pi_phy,rho_phy
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme)::      &
     &                      th_phy,qv,qc,qr,qi,qrimef
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme)           ::     &
     &                                                   RAINNC,RAINNCV
      REAL, INTENT(OUT),    DIMENSION(ims:ime,jms:jme):: SR





      INTEGER, DIMENSION( ims:ime, jms:jme ),INTENT(INOUT) :: LOWLYR





      REAL,                 DIMENSION(ims:ime, kms:kme, jms:jme ) ::    &
     &                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY, QT








      REAL,  DIMENSION( ims:ime, kms:kme, jms:jme ) ::                  &
     &       TLATGS_PHY,TRAIN_PHY
      REAL,  DIMENSION(ims:ime,jms:jme):: APREC,PREC,ACPREC
      REAL,  DIMENSION(its:ite, kts:kte, jts:jte):: t_phy

      INTEGER :: I,J,K,KFLIP
      REAL :: WC





















          CALL MY_GROWTH_RATES (DT)





        CIACW=DT*0.25*PI_E*0.5*(1.E5)**C1




        CIACR=PI_E*DT




        CRACW=DT*0.25*PI_E*1.0



        BRAUT=DT*1.1E10*BETA6/NCW











      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
        t_phy(i,k,j) = th_phy(i,k,j)*pi_phy(i,k,j)
        qv(i,k,j)=qv(i,k,j)/(1.+qv(i,k,j)) 
      ENDDO
      ENDDO
      ENDDO



      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
	 TLATGS_PHY (i,k,j)=0.
	 TRAIN_PHY  (i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO

      DO j = jts,jte
      DO i = its,ite
         ACPREC(i,j)=0.
         APREC (i,j)=0.
         PREC  (i,j)=0.
         SR    (i,j)=0.
      ENDDO
      ENDDO



      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
         QT(I,K,J)=QC(I,K,J)+QR(I,K,J)+QI(I,K,J)
         IF (QI(I,K,J) <= EPSQ) THEN
            F_ICE_PHY(I,K,J)=0.
            F_RIMEF_PHY(I,K,J)=1.
            IF (T_PHY(I,K,J) < T_ICEK) F_ICE_PHY(I,K,J)=1.
         ELSE
            F_ICE_PHY(I,K,J)=MAX( 0., MIN(1., QI(I,K,J)/QT(I,K,J) ) )
            F_RIMEF_PHY(I,K,J)=QRIMEF(I,K,J)/QI(I,K,J)
         ENDIF
         IF (QR(I,K,J) <= EPSQ) THEN
            F_RAIN_PHY(I,K,J)=0.
         ELSE
            F_RAIN_PHY(I,K,J)=QR(I,K,J)/(QR(I,K,J)+QC(I,K,J))
         ENDIF
      ENDDO
      ENDDO
      ENDDO



      CALL EGCP01DRV(GID,DT,LOWLYR,                                     &
     &               APREC,PREC,ACPREC,SR,                              &
     &               dz8w,rho_phy,qt,t_phy,qv,F_ICE_PHY,P_PHY,          &
     &               F_RAIN_PHY,F_RIMEF_PHY,TLATGS_PHY,TRAIN_PHY,       &
     &               ids,ide, jds,jde, kds,kde,		                &
     &               ims,ime, jms,jme, kms,kme,		                &
     &               its,ite, jts,jte, kts,kte		          )


     DO j = jts,jte
        DO k = kts,kte
	DO i = its,ite
  	   th_phy(i,k,j) = t_phy(i,k,j)/pi_phy(i,k,j)
           qv(i,k,j)=qv(i,k,j)/(1.-qv(i,k,j))  
           WC=qt(I,K,J)
           QI(I,K,J)=0.
           QR(I,K,J)=0.
           QC(I,K,J)=0.
           IF(F_ICE_PHY(I,K,J)>=1.)THEN
             QI(I,K,J)=WC
           ELSEIF(F_ICE_PHY(I,K,J)<=0.)THEN
             QC(I,K,J)=WC
           ELSE
             QI(I,K,J)=F_ICE_PHY(I,K,J)*WC
             QC(I,K,J)=WC-QI(I,K,J)
           ENDIF

           IF(QC(I,K,J)>0..AND.F_RAIN_PHY(I,K,J)>0.)THEN
             IF(F_RAIN_PHY(I,K,J).GE.1.)THEN
               QR(I,K,J)=QC(I,K,J)
               QC(I,K,J)=0.
             ELSE
               QR(I,K,J)=F_RAIN_PHY(I,K,J)*QC(I,K,J)
               QC(I,K,J)=QC(I,K,J)-QR(I,K,J)
             ENDIF
          endif
          QRIMEF(I,K,J)=QI(I,K,J)*F_RIMEF_PHY(I,K,J)
	ENDDO
        ENDDO
     ENDDO



       DO j=jts,jte
       DO i=its,ite
          RAINNC(i,j)=APREC(i,j)*1000.+RAINNC(i,j)
          RAINNCV(i,j)=APREC(i,j)*1000.
       ENDDO
       ENDDO
















  END SUBROUTINE FER_HIRES_ADVECT




      SUBROUTINE EGCP01DRV(GID,                            & 
     &  DTPH,LOWLYR,APREC,PREC,ACPREC,SR,                  &
     &  dz8w,RHO_PHY,CWM_PHY,T_PHY,Q_PHY,F_ICE_PHY,P_PHY,  &
     &  F_RAIN_PHY,F_RIMEF_PHY,TLATGS_PHY,TRAIN_PHY,       &
     &  ids,ide, jds,jde, kds,kde,                         &
     &  ims,ime, jms,jme, kms,kme,                         &
     &  its,ite, jts,jte, kts,kte)
















      IMPLICIT NONE



      INTEGER,INTENT(IN ) :: ids,ide, jds,jde, kds,kde                  &
     &                      ,ims,ime, jms,jme, kms,kme                  &
     &                      ,its,ite, jts,jte, kts,kte
      INTEGER,INTENT(IN ) :: GID     
      REAL,INTENT(IN) :: DTPH
      INTEGER, DIMENSION( ims:ime, jms:jme ),INTENT(INOUT) :: LOWLYR
      REAL,DIMENSION(ims:ime,jms:jme),INTENT(INOUT) ::                  &
     &                         APREC,PREC,ACPREC,SR
      REAL,DIMENSION( its:ite, kts:kte, jts:jte ),INTENT(INOUT) :: t_phy
      REAL,DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN) ::         &
     &                                             dz8w,P_PHY,RHO_PHY
      REAL,DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(INOUT) ::      &
     &   CWM_PHY, F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY,TLATGS_PHY           &
     &   ,Q_PHY,TRAIN_PHY

















      INTEGER :: LSFC,I,J,I_index,J_index,L,K,KFLIP

      REAL,DIMENSION(its:ite,jts:jte,kts:kte) ::                        &
     &   CWM,T,Q,TRAIN,TLATGS,P
      REAL,DIMENSION(kts:kte,its:ite,jts:jte) :: F_ice,F_rain,F_RimeF       
      INTEGER,DIMENSION(its:ite,jts:jte) :: LMH
      REAL :: TC,WC,QI,QR,QW,Fice,Frain,DUM,ASNOW,ARAIN
      REAL,DIMENSION(kts:kte) :: P_col,Q_col,T_col,QV_col,WC_col,       &
     & RimeF_col,QI_col,QR_col,QW_col, THICK_col, RHC_col, DPCOL,       &
     & pcond1d,pidep1d,                                                 &
     & piacw1d,piacwi1d,piacwr1d,piacr1d,picnd1d,pievp1d,pimlt1d,       &
     & praut1d,pracw1d,prevp1d,pisub1d,pevap1d, DBZ_col,NR_col,NS_col,  &
     & vsnow1d,vrain11d,vrain21d,vci1d,NSmICE1d,INDEXS1d,INDEXR1d,      & 
     & RFlag1d   
      REAL,DIMENSION(2) :: PRECtot,PRECmax


        DO J=JTS,JTE    
        DO I=ITS,ITE  
           LMH(I,J) = KTE-LOWLYR(I,J)+1
        ENDDO
        ENDDO


        DO 98  J=JTS,JTE    
        DO 98  I=ITS,ITE  
           DO L=KTS,KTE
             KFLIP=KTE+1-L
             CWM(I,J,L)=CWM_PHY(I,KFLIP,J)
             T(I,J,L)=T_PHY(I,KFLIP,J)
             Q(I,J,L)=Q_PHY(I,KFLIP,J)
             P(I,J,L)=P_PHY(I,KFLIP,J)
             TLATGS(I,J,L)=TLATGS_PHY(I,KFLIP,J)
             TRAIN(I,J,L)=TRAIN_PHY(I,KFLIP,J)
             F_ice(L,I,J)=F_ice_PHY(I,KFLIP,J)
             F_rain(L,I,J)=F_rain_PHY(I,KFLIP,J)
             F_RimeF(L,I,J)=F_RimeF_PHY(I,KFLIP,J)
           ENDDO
98      CONTINUE
     
       DO 100 J=JTS,JTE    
        DO 100 I=ITS,ITE  
          LSFC=LMH(I,J)                      

          DO K=KTS,KTE
            KFLIP=KTE+1-K
            DPCOL(K)=RHO_PHY(I,KFLIP,J)*GRAV*dz8w(I,KFLIP,J)
          ENDDO

   
   
   
          IF (CWM(I,J,1) .LE. EPSQ) CWM(I,J,1)=EPSQ
          F_ice(1,I,J)=1.
          F_rain(1,I,J)=0.
          F_RimeF(1,I,J)=1.
          DO L=1,LSFC
      
      
      
            P_col(L)=P(I,J,L)
      
      
      
            THICK_col(L)=DPCOL(L)*RGRAV
            T_col(L)=T(I,J,L)
            TC=T_col(L)-T0C
            QV_col(L)=max(EPSQ, Q(I,J,L))
            IF (CWM(I,J,L) .LE. EPSQ1) THEN
              WC_col(L)=0.
              IF (TC .LT. T_ICE) THEN
                F_ice(L,I,J)=1.
              ELSE
                F_ice(L,I,J)=0.
              ENDIF
              F_rain(L,I,J)=0.
              F_RimeF(L,I,J)=1.
            ELSE
              WC_col(L)=CWM(I,J,L)


IF (WC_col(L)>QTwarn .AND. P_col(L)<Pwarn .AND. TC==TC) THEN
   WRITE(0,*) 'WARN4: >1 g/kg condensate in stratosphere; I,J,L,TC,P,QT=',   &
              I,J,L,TC,.01*P_col(L),1000.*WC_col(L)
   QTwarn=MAX(WC_col(L),10.*QTwarn)
   Pwarn=MIN(P_col(L),0.5*Pwarn)
ENDIF

IF (WARN5 .AND. TC/=TC) THEN
   WRITE(0,*) 'WARN5: NaN temperature; I,J,L,P=',I,J,L,.01*P_col(L)
   WARN5=.FALSE.
ENDIF

            ENDIF
            IF (T_ICE<=-100.) F_ice(L,I,J)=0.  
      
      
      
      
            WC=WC_col(L)
            QI=0.
            QR=0.
            QW=0.
            Fice=F_ice(L,I,J)
            Frain=F_rain(L,I,J)
            IF (Fice .GE. 1.) THEN
              QI=WC
            ELSE IF (Fice .LE. 0.) THEN
              QW=WC
            ELSE
              QI=Fice*WC
              QW=WC-QI
            ENDIF
            IF (QW.GT.0. .AND. Frain.GT.0.) THEN
              IF (Frain .GE. 1.) THEN
                QR=QW
                QW=0.
              ELSE
                QR=Frain*QW
                QW=QW-QR
              ENDIF
            ENDIF
            IF (QI .LE. 0.) F_RimeF(L,I,J)=1.
            RimeF_col(L)=F_RimeF(L,I,J)
            QI_col(L)=QI
            QR_col(L)=QR
            QW_col(L)=QW




            IF(GID .EQ. 1 .AND. P_col(L)<P_RHgrd_out) THEN  
              RHC_col(L)=RHgrd_out        
            ELSE
              RHC_col(L)=RHgrd_in       
            ENDIF

          ENDDO


   
   
   
          I_index=I
          J_index=J

       CALL EGCP01COLUMN ( ARAIN, ASNOW, DTPH, RHC_col,                 &
     & I_index, J_index, LSFC,                                          &
     & P_col, QI_col, QR_col, QV_col, QW_col, RimeF_col, T_col,         &
     & THICK_col, WC_col, KTS,KTE,              pcond1d,pidep1d,        &
     & piacw1d,piacwi1d,piacwr1d,piacr1d,picnd1d,pievp1d,pimlt1d,       &
     & praut1d,pracw1d,prevp1d,pisub1d,pevap1d, DBZ_col,NR_col,NS_col,  &
     & vsnow1d,vrain11d,vrain21d,vci1d,NSmICE1d,INDEXS1d,INDEXR1d,      & 
     & RFlag1d)   

   


   
   
   
          DO L=1,LSFC
            TRAIN(I,J,L)=(T_col(L)-T(I,J,L))/DTPH
            TLATGS(I,J,L)=T_col(L)-T(I,J,L)
            T(I,J,L)=T_col(L)
            Q(I,J,L)=QV_col(L)
            CWM(I,J,L)=WC_col(L)
      
      
      
            IF (QI_col(L) .LE. EPSQ) THEN
              F_ice(L,I,J)=0.
              IF (T_col(L) .LT. T_ICEK) F_ice(L,I,J)=1.
              F_RimeF(L,I,J)=1.
            ELSE
              F_ice(L,I,J)=MAX( 0., MIN(1., QI_col(L)/WC_col(L)) )
              F_RimeF(L,I,J)=MAX(1., RimeF_col(L))
            ENDIF
            IF (QR_col(L) .LE. EPSQ) THEN
              DUM=0
            ELSE
              DUM=QR_col(L)/(QR_col(L)+QW_col(L))
            ENDIF
            F_rain(L,I,J)=DUM
      
          ENDDO
   
   
   
   
   
   
        APREC(I,J)=(ARAIN+ASNOW)*RRHOL       
        PREC(I,J)=PREC(I,J)+APREC(I,J)
        ACPREC(I,J)=ACPREC(I,J)+APREC(I,J)
        IF(APREC(I,J) .LT. 1.E-8) THEN
          SR(I,J)=0.
        ELSE
          SR(I,J)=RRHOL*ASNOW/APREC(I,J)
        ENDIF














100   CONTINUE                          
        DO 101 J=JTS,JTE    
        DO 101 I=ITS,ITE  
           DO L=KTS,KTE
              KFLIP=KTE+1-L
             CWM_PHY(I,KFLIP,J)=CWM(I,J,L)
             T_PHY(I,KFLIP,J)=T(I,J,L)
             Q_PHY(I,KFLIP,J)=Q(I,J,L)
             TLATGS_PHY(I,KFLIP,J)=TLATGS(I,J,L)
             TRAIN_PHY(I,KFLIP,J)=TRAIN(I,J,L)
             F_ice_PHY(I,KFLIP,J)=F_ice(L,I,J)
             F_rain_PHY(I,KFLIP,J)=F_rain(L,I,J)
             F_RimeF_PHY(I,KFLIP,J)=F_RimeF(L,I,J)
           ENDDO
101     CONTINUE

      END SUBROUTINE EGCP01DRV













































      SUBROUTINE EGCP01COLUMN ( ARAIN, ASNOW, DTPH, RHC_col,             &
     & I_index, J_index, LSFC,                                           &
     & P_col, QI_col, QR_col, Q_col, QW_col, RimeF_col, T_col,           &
     & THICK_col, WC_col ,KTS,KTE,pcond1d,pidep1d,                       &
     & piacw1d,piacwi1d,piacwr1d,piacr1d,picnd1d,pievp1d,pimlt1d,        &
     & praut1d,pracw1d,prevp1d,pisub1d,pevap1d, DBZ_col,NR_col,NS_col,   &
     & vsnow1d,vrain11d,vrain21d,vci1d,NSmICE1d,INDEXS1d,INDEXR1d,       &  
     & RFlag1d)   













































































      IMPLICIT NONE

      INTEGER,INTENT(IN) :: KTS,KTE,I_index, J_index, LSFC
      REAL,INTENT(IN)    :: DTPH
      REAL,INTENT(INOUT) ::  ARAIN, ASNOW
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: P_col,QI_col,QR_col,RHC_col  &
     & ,Q_col ,QW_col, RimeF_col, T_col, THICK_col,WC_col,pcond1d       &
     & ,pidep1d,piacw1d,piacwi1d,piacwr1d,piacr1d,picnd1d,pievp1d       &
     & ,pimlt1d,praut1d,pracw1d,prevp1d,pisub1d,pevap1d,DBZ_col,NR_col  &
     & ,NS_col,vsnow1d,vrain11d,vrain21d,vci1d,NSmICE1d,INDEXS1d        &  
     & ,INDEXR1d,RFlag1d    











































      REAL, PARAMETER :: TOLER=5.E-7, C2=1./6., RHO0=1.194,             &
                         Xratio=.025, Zmin=0.01, DBZmin=-20.











      REAL, PARAMETER :: BLEND=1.



      LOGICAL, PARAMETER :: PRINT_diag=.TRUE.





      REAL :: EMAIRI, N0r, NLICE, NSmICE, NInuclei, Nrain, Nsnow, Nmix
      LOGICAL :: CLEAR, ICE_logical, DBG_logical, RAIN_logical,         &
                 STRAT, DRZL
      INTEGER :: INDEX_MY,INDEXR,INDEXR1,INDEXR2,INDEXS,IPASS,ITDX,IXRF,&
     &           IXS,LBEF,L,INDEXRmin,INDEXS0C,IDR        


      REAL :: ABI,ABW,AIEVP,ARAINnew,ASNOWnew,BLDTRH,BUDGET,            &
     &        CREVP,DELI,DELR,DELT,DELV,DELW,DENOMF,                    &
     &        DENOMI,DENOMW,DENOMWI,DIDEP,                              &
     &        DIEVP,DIFFUS,DLI,DTRHO,DUM,DUM1,DUM2,DUM3,                &
     &        DWV0,DWVI,DWVR,DYNVIS,ESI,ESW,FIR,FLIMASS,                &
     &        FWR,FWS,GAMMAR,GAMMAS,                                    &
     &        PCOND,PIACR,PIACW,PIACWI,PIACWR,PICND,PIDEP,PIDEP_max,    &
     &        PIEVP,PILOSS,PIMLT,PINIT,PP,PRACW,PRAUT,PREVP,PRLOSS,     &
     &        QI,QInew,QLICE,QR,QRnew,QSI,QSIgrd,QSInew,QSW,QSW0,       &
     &        QSWgrd,QSWnew,QT,QTICE,QTnew,QTRAIN,Q,QW,QWnew,Rcw,       &
     &        RFACTOR,RFmx,RFrime,RHO,RIMEF,RIMEF1,RQR,RR,RRHO,SFACTOR, &
     &        TC,TCC,TFACTOR,THERM_COND,THICK,TK,TK2,TNEW,              &
     &        TOT_ICE,TOT_ICEnew,TOT_RAIN,TOT_RAINnew,                  &
     &        VEL_INC,VENTR,VENTIL,VENTIS,VRAIN1,VRAIN2,VRIMEF,VSNOW,   &
     &        VSNOW1,WC,WCnew,WSgrd,WS,WSnew,WV,WVnew,                  &
     &        XLI,XLIMASS,XRF, RHgrd,                                   &
     &        NSImax,QRdum,QSmICE,QLgIce,RQLICE,VCI,TIMLT,              &
     &        RQSnew,RQRnew,Zrain,Zsnow,Ztot,RHOX0C,RFnew,PSDEP,DELS     
      REAL, SAVE :: Revised_LICE=1.e-3    






      ARAIN=0.                
      VRAIN1=0.               
      VSNOW1=0.               
      ASNOW=0.                
      INDEXS0C=MDImin         
      RHOX0C=22.5             
      TIMLT=0.                
      STRAT=.FALSE.           
      DRZL=.FALSE.            





big_loop: DO L=1,LSFC
        pcond1d(L)=0.
        pidep1d(L)=0.
        piacw1d(L)=0.
        piacwi1d(L)=0.
        piacwr1d(L)=0.
        piacr1d(L)=0.
        picnd1d(L)=0.
        pievp1d(L)=0.
        pimlt1d(L)=0.
        praut1d(L)=0.
        pracw1d(L)=0.
        prevp1d(L)=0.
        pisub1d(L)=0.
        pevap1d(L)=0.
        vsnow1d(L)=0.
        vrain11d(L)=0.
        vrain21d(L)=0.
        vci1d(L)=0.
        NSmICE1d(L)=0.
        DBZ_col(L)=DBZmin
        NR_col(L)=0.
        NS_col(L)=0.
        INDEXR1d(L)=0.
        INDEXS1d(L)=0.
        RFlag1d(L)=0.   






        IF (ARAIN .GT. CLIMIT) THEN
          CLEAR=.FALSE.
          VRAIN1=0.
        ELSE
          CLEAR=.TRUE.
          ARAIN=0.
        ENDIF






        IF (ASNOW .GT. CLIMIT) THEN
          CLEAR=.FALSE.
          VSNOW1=0.
        ELSE
          ASNOW=0.
        ENDIF





        TK=T_col(L)         
        TC=TK-T0C           
        PP=P_col(L)         
        Q=Q_col(L)          
        WV=Q/(1.-Q)         
        WC=WC_col(L)        
        RHgrd=RHC_col(L)    







        ESW=MIN(1000.*FPVS0(TK),0.99*PP) 
        QSW=EPS*ESW/(PP-ESW)             
        WS=QSW                           
        QSI=QSW                          
        IF (TC .LT. 0.) THEN
          ESI=MIN(1000.*FPVS(TK),0.99*PP)  
          QSI=EPS*ESI/(PP-ESI)           
          WS=QSI                         
        ENDIF



        QSWgrd=RHgrd*QSW
        QSIgrd=RHgrd*QSI
        WSgrd=RHgrd*WS



        IF (WV.GT.WSgrd .OR. WC.GT.EPSQ) CLEAR=.FALSE.





        IF (CLEAR) THEN
          STRAT=.FALSE.     
          DRZL=.FALSE.      
          INDEXRmin=MDRmin  
          TIMLT=0.          
          CYCLE big_loop
        ENDIF









          RHO=PP/(RD*TK*(1.+EPS1*Q)) 
          RRHO=1./RHO                
          DTRHO=DTPH*RHO             
          BLDTRH=BLEND*DTRHO         
          THICK=THICK_col(L)         

          ARAINnew=0.                
          ASNOWnew=0.                
          QI=QI_col(L)               
          QInew=0.                   
          QR=QR_col(L)               
          QRnew=0.                   
          QW=QW_col(L)               
          QWnew=0.                   

          PCOND=0.                   
          PIDEP=0.                   
          PINIT=0.                   
          PIACW=0.                   
          PIACWI=0.                  
          PIACWR=0.                  
          PIACR=0.                   
          PICND=0.                   
          PIEVP=0.                   
          PIMLT=0.                   
          PRAUT=0.                   
          PRACW=0.                   
          PREVP=0.                   
          NSmICE=0.                  
          Nrain=0.                   
          Nsnow=0.                   
          RQRnew=0.                  
          RQSnew=0.                  
          Zrain=0.                   
          Zsnow=0.                   
          Ztot=0.                    
          INDEXR=MDRmin              
          INDEXR1=INDEXR             
          INDEXR2=INDEXR             
          N0r=0.                     
          DUM1=MIN(0.,TC)
          DUM=XMImax*EXP(XMIexp*DUM1)
          INDEXS=MIN(MDImax, MAX(MDImin, INT(DUM) ) )   
          VCI=0.                     
          VSNOW=0.                   
          VRAIN2=0.                  
          RimeF1=1.                  























          TK2=1./(TK*TK)             






          TFACTOR=SQRT(TK*TK*TK)/(TK+120.)
          DYNVIS=1.496E-6*TFACTOR
          THERM_COND=2.116E-3*TFACTOR
          DIFFUS=8.794E-5*TK**1.81/PP




          GAMMAS=MIN(1.5, (1.E5/PP)**C1)    



          GAMMAR=(RHO0/RHO)**.4







          IF (TC.LT.0. .OR. QI.GT. EPSQ .OR. ASNOW.GT.CLIMIT) THEN
            ICE_logical=.TRUE.
          ELSE
            ICE_logical=.FALSE.
            QLICE=0.
            QTICE=0.
          ENDIF
          IF (T_ICE <= -100.) THEN
            ICE_logical=.FALSE.
            QLICE=0.
            QTICE=0.
          ENDIF



          RAIN_logical=.FALSE.
          IF (ARAIN.GT.CLIMIT .OR. QR.GT.EPSQ) RAIN_logical=.TRUE.

ice_test: IF (ICE_logical) THEN






























            NInuclei=0.
            NSmICE=0.
            QSmICE=0. 
            Rcw=0.
            IF (TC<0.) THEN




              NSImax=MAX(NSI_max, 0.1*RHO*QI/MASSI(MDImin) )  










              NInuclei=MIN(0.01*EXP(-0.6*TC), NSImax)         
              NInuclei=MIN(5.*EXP(-0.304*TC), NInuclei)       
              IF (QI>EPSQ) THEN
                DUM=RRHO*MASSI(MDImin)
                NSmICE=MIN(NInuclei, QI/DUM)
                QSmICE=NSmICE*DUM
              ENDIF         
            ENDIF            
  init_ice: IF (QI<=EPSQ .AND. ASNOW<=CLIMIT) THEN
              TOT_ICE=0.
              PILOSS=0.
              RimeF1=1.
              VrimeF=1.
              VEL_INC=GAMMAS
              VSNOW=0.
              VSNOW1=0.
              VCI=0.
              EMAIRI=THICK
              XLIMASS=RimeF1*MASSI(INDEXS)
              FLIMASS=1.
              QLICE=0.
              RQLICE=0.
              QTICE=0.
              NLICE=0.
            ELSE  init_ice
   
   
   
   
   



              TOT_ICE=THICK*QI+BLEND*ASNOW
              PILOSS=-TOT_ICE/THICK
              QLgICE=MAX(0., QI-QSmICE)         
              VCI=GAMMAS*VSNOWI(MDImin)



              LBEF=MAX(1,L-1)
              RimeF1=(RimeF_col(L)*THICK*QLgICE                         &
     &                +RimeF_col(LBEF)*BLEND*ASNOW)/TOT_ICE


              IF (RimeF1>2.) THEN
                DUM3=RH_NgC*(RHO*QLgICE)**C1         
                DUM2=RH_NgT*(RHO*QLgICE)**C1         
                IF (RimeF1>=10.) THEN
                  DUM=DUM3
                ELSE IF (RimeF1>=5.) THEN
                  DUM=0.2*(RimeF1-5.)                
                  DUM=DUM3*DUM+DUM2*(1.-DUM)
                ELSE
                  DUM1=REAL(INDEXS)                  
                  DUM=0.33333*(RimeF1-2.)            
                  DUM=DUM2*DUM+DUM1*(1.-DUM)
                ENDIF
                INDEXS=MIN(MDImax, MAX(MDImin, INT(DUM) ) )
              ENDIF

              EMAIRI=THICK+BLDTRH*VSNOW1
              QLICE=(THICK*QLgICE+BLEND*ASNOW)/EMAIRI   

              RQLICE=RHO*QLICE
              QTICE=QLICE+QSmICE
              FLIMASS=QLICE/QTICE





    two_pass: DO IPASS=1,2





                DUM=1.E-6*REAL(INDEXS)          
                RFmx=RFmx1*DUM*DUM*DUM/MASSI(INDEXS)
                RimeF1=MIN(RimeF1, RFmx)

    vel_rime:   IF (RimeF1<=1.) THEN
                  RimeF1=1.
                  VrimeF=1.
                ELSE   vel_rime
  
                  RimeF1=MIN(RimeF1, RFmax)
                  IXS=MAX(2, MIN(INDEXS/100, 9))
                  XRF=10.492*ALOG(RimeF1)
                  IXRF=MAX(0, MIN(INT(XRF), Nrime))
                  IF (IXRF .GE. Nrime) THEN
                    VrimeF=VEL_RF(IXS,Nrime)
                  ELSE
                    VrimeF=VEL_RF(IXS,IXRF)+(XRF-FLOAT(IXRF))*            &
       &                   (VEL_RF(IXS,IXRF+1)-VEL_RF(IXS,IXRF))
                  ENDIF
                  VrimeF=MAX(1., VrimeF)
                ENDIF   vel_rime


                VEL_INC=GAMMAS*VrimeF           
                VSNOW=VEL_INC*VSNOWI(INDEXS)
                IF (RimeF1>=5. .AND. INDEXS==MDImax .AND. RQLICE>RQhail) THEN

                  DUM=GAMMAS*AVhail*RQLICE**BVhail
                  IF (DUM>VSNOW) THEN
                    VSNOW=DUM
                    VEL_INC=VSNOW/VSNOWI(INDEXS)
                  ENDIF
                ENDIF
                XLIMASS=RimeF1*MASSI(INDEXS)
                NLICE=RQLICE/XLIMASS




                IF (IPASS>=2 .OR.                                       &
                    (NLICE>=NLImin .AND. NLICE<=NSI_max)) EXIT two_pass






                NLICE=MAX(NLImin, MIN(NSI_max, NLICE) )


                XLI=RQLICE/(NLICE*RimeF1)         
new_size:       IF (XLI<=MASSI(MDImin) ) THEN
                  INDEXS=MDImin
                ELSE IF (XLI<=MASSI(450) ) THEN   new_size
                  DLI=9.5885E5*XLI**.42066         
                  INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                ELSE IF (XLI<MASSI(MDImax) ) THEN   new_size
                  DLI=3.9751E6*XLI**.49870         
                  INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                ELSE  new_size
                  INDEXS=MDImax






                ENDIF    new_size

              ENDDO      two_pass

            ENDIF        init_ice
          ENDIF          ice_test


          IF(TC<=0.) THEN
            STRAT=.FALSE.
            INDEXRmin=MDRmin
            TIMLT=0.
            INDEXS0C=INDEXS
            RHOX0C=22.5*MAX(1.,MIN(RimeF1,40.))     
          ELSE   
            IF(.NOT.RAIN_logical) THEN
              STRAT=.FALSE.     
              INDEXRmin=MDRmin  
              IF(.NOT.ICE_logical) TIMLT=0.
            ELSE

              IF(TIMLT>EPSQ .AND. RHOX0C<=225.) THEN
                STRAT=.TRUE.
              ELSE
                STRAT=.FALSE.     
                INDEXRmin=MDRmin  
              ENDIF
              IF(STRAT .AND. INDEXRmin<=MDRmin) THEN
                INDEXRmin=INDEXS0C*(0.001*RHOX0C)**C1
                INDEXRmin=MAX(MDRmin, MIN(INDEXRmin, INDEXRstrmax) )
              ENDIF
            ENDIF
          ENDIF

          IF(STRAT .OR. TIMLT>EPSQ) THEN
            DRZL=.FALSE.
          ELSE

            DRZL=.TRUE.    
          ENDIF









          IF (QW.GT.EPSQ .AND. TC.GE.T_ICE) THEN

            DUM2=RHO*QW
            IF (DUM2>QAUT0) THEN

              DUM2=DUM2*DUM2
              DUM=BRAUT*DUM2*QW
              DUM1=ARAUT*DUM2
              PRAUT=MIN(QW, DUM*(1.-EXP(-DUM1*DUM1)) )
            ENDIF
            IF (QLICE .GT. EPSQ) THEN
      
      
      
      
              FWS=MIN(1., CIACW*VEL_INC*NLICE*ACCRI(INDEXS) ) 
              PIACW=FWS*QW
              IF (TC<0.) THEN
                PIACWI=PIACW            
                Rcw=ARcw*DUM2**C1       
              ENDIF
            ENDIF           
          ENDIF             





ice_only: IF (TC.LT.T_ICE .AND. (WV.GT.QSWgrd .OR. QW.GT.EPSQ)) THEN
   
   
   
   
            PIACW=QW
            PIACWI=PIACW
            Rcw=0.                             
                                               
            DUM1=TK+XLF1*PIACW                 
            DUM2=WV                            
            DUM=MIN(1000.*FPVS(DUM1),0.99*PP)  
            DUM=RHgrd*EPS*DUM/(PP-DUM)         


IF (WARN1 .AND. DUM1<XMIN) THEN
   WRITE(0,*) 'WARN1: Water saturation T<180K; I,J,L,TC,P=',   &
              I_index,J_index,L,DUM1-T0C,.01*PP
   WARN1=.FALSE.
ENDIF


            IF (DUM2>DUM)                                               &
     &          PIDEP=DEPOSIT(PP,DUM1,DUM2,RHgrd,I_index,J_index,L)

            DWVi=0.    
   
          ELSE IF (TC<0.) THEN  ice_only
   
   
   
   
            DENOMI=1.+XLS3*QSI*TK2

            DWVi=WV-QSIgrd          
            PIDEP_max=MAX(PILOSS, DWVi/DENOMI)
            DIDEP=0.     
            IF (QTICE .GT. 0.) THEN
      
      
      
      
      
      
      
      
              ABI=1./(RHO*XLS2*QSI*TK2/THERM_COND+1./DIFFUS)
      
      
      
      
      
      
      
              DUM=(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
              SFACTOR=SQRT(GAMMAS)*DUM              
              VENTIS=(VENTI1(MDImin)+SFACTOR*VENTI2(MDImin))*NSmICE
              SFACTOR=SQRT(VEL_INC)*DUM
              VENTIL=(VENTI1(INDEXS)+SFACTOR*VENTI2(INDEXS))*NLICE
              DIDEP=ABI*(VENTIL+VENTIS)*DTPH
            ENDIF   
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
            IF (WV>QSWgrd .AND. TC<T_ICE_init .AND. NSmICE<NInuclei) THEN
              INDEX_MY=MAX(MY_T1, MIN( INT(.5-TC), MY_T2 ) )
      
              DUM=MAX(0., NInuclei-NSmICE)
              
              PINIT=MAX(0., DUM*MY_GROWTH(INDEX_MY)*RRHO)
            ENDIF
      
      
      
            IF (DWVi>0.) THEN
              PIDEP=MIN(DWVI*DIDEP+PINIT, PIDEP_max)
            ELSE IF (DWVi<0.) THEN
              PIDEP=MAX(DWVi*DIDEP, PIDEP_max)
            ENDIF
   
          ENDIF  ice_only




          IF (TC>=T_ICE .AND. (QW>EPSQ .OR. WV>QSWgrd)) THEN

            DUM1=QW-PIACWI
            DUM2=TK+XLS1*PIDEP+XLF1*PIACWI
            DUM3=WV-PIDEP
            PCOND=CONDENSE(PP,DUM1,DUM2,DUM3,RHgrd,I_index,J_index,L)

          ENDIF




          TCC=TC+XLV1*PCOND+XLS1*PIDEP+XLF1*PIACWI
          IF (TCC>0.) THEN
            PIACWI=0.
            Rcw=0.              
            PIDEP=0.            
            TCC=TC+XLV1*PCOND   
          ENDIF

          QSW0=0.
          IF (TC.GT.0. .AND. TCC.GT.0. .AND. ICE_logical) THEN
   
   
   
   
   
   
            SFACTOR=SQRT(VEL_INC)*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
            VENTIL=NLICE*(VENTI1(INDEXS)+SFACTOR*VENTI2(INDEXS))
            AIEVP=VENTIL*DIFFUS*DTPH
            IF (AIEVP .LT. Xratio) THEN
              DIEVP=AIEVP
            ELSE
              DIEVP=1.-EXP(-AIEVP)
            ENDIF
            QSW0=EPS*ESW0/(PP-ESW0)
            DWV0=MIN(WV,QSW)-QSW0
            DUM=QW+PCOND
            IF (WV.LT.QSW .AND. DUM.LE.EPSQ) THEN
   
   
   
   
              DUM=DWV0*DIEVP
              PIEVP=MAX( MIN(0., DUM), PILOSS)
              PICND=MAX(0., DUM)
            ENDIF            
            PIMLT=THERM_COND*TCC*VENTIL*RRHO*DTPH/XLF
   
   
   
            DUM1=MAX( 0., (TCC+XLV1*PIEVP)/XLF1 )
            PIMLT=MIN(PIMLT, DUM1)
   
   
   
            DUM=PIEVP-PIMLT
            IF (DUM .LT. PILOSS) THEN
              DUM1=PILOSS/DUM
              PIMLT=PIMLT*DUM1
              PIEVP=PIEVP*DUM1
            ENDIF           
          ENDIF             















          TOT_RAIN=0.
          QTRAIN=0.
          PRLOSS=0.
          RQR=0.
do_rain:  IF (RAIN_logical) THEN
            TOT_RAIN=THICK*QR+BLEND*ARAIN   
            PRLOSS=-TOT_RAIN/THICK

            QTRAIN=TOT_RAIN/(THICK+BLDTRH*VRAIN1) 
            RQR=RHO*QTRAIN                        
            IF(STRAT .AND. RQR>1.E-3) STRAT=.FALSE.    
            IF(DRZL .AND. PIMLT>EPSQ) DRZL=.FALSE.     

DSD1:       IF (RQR<=RQR_DRmin) THEN

              INDEXR=MDRmin
              N0r=RQR/MASSR(INDEXR)
            ELSE IF (DRZL) THEN  DSD1

              DUM=MAX(1.0, 0.5E-3/RQR)
              N0r=MIN(1.E9, N0r0*DUM*SQRT(DUM) )   
              DUM1=RQR/(PI*RHOL*N0r)
              DUM=1.E6*SQRT(SQRT(DUM1))
              INDEXR=MAX(XMRmin, MIN(DUM, XMRmax) )
              N0r=RQR/MASSR(INDEXR)

            ELSE IF (RQR>=RQR_DRmax) THEN  DSD1

              INDEXR=MDRmax
              N0r=RQR/MASSR(INDEXR)
            ELSE  DSD1

              N0r=N0r0
              DUM=CN0r0*SQRT(SQRT(RQR))
              INDEXR=MAX(XMRmin, MIN(DUM, XMRmax) )
              IF (STRAT .AND. INDEXR<INDEXRmin) THEN

                INDEXR=INDEXRmin
                N0r=RQR/MASSR(INDEXR)
              ENDIF
            ENDIF  DSD1

            VRAIN2=GAMMAR*VRAIN(INDEXR)    


            IF (TC .LT. T_ICE) THEN
              PIACR=-PRLOSS
            ELSE
              DWVr=WV-PCOND-QSW
              DUM=QW+PCOND
              IF (DWVr.LT.0. .AND. DUM.LE.EPSQ) THEN
      
      
      
      
      
      
      
      
      
      
                RFACTOR=SQRT(GAMMAR)*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
                ABW=1./(RHO*XLV2*QSW*TK2/THERM_COND+1./DIFFUS)
      
      
      
      
                VENTR=N0r*(VENTR1(INDEXR)+RFACTOR*VENTR2(INDEXR))
                CREVP=ABW*VENTR*DTPH
                PREVP=MAX(DWVr*CREVP, PRLOSS)
              ELSE IF (QW .GT. EPSQ) THEN
                FWR=CRACW*GAMMAR*N0r*ACCRR(INDEXR)
                PRACW=MIN(1.,FWR)*QW
              ENDIF           
      
              IF (ICE_logical .AND. TC<0. .AND. TCC<0.) THEN
         
         
         
         
                DUM=.001*FLOAT(INDEXR)
                DUM=(EXP(ABFR*TC)-1.)*DUM*DUM*DUM*DUM*DUM*DUM*DUM
                PIACR=MIN(CBFR*N0r*RRHO*DUM, QTRAIN)

                FIR=0.

                 IF (QLICE .GT. EPSQ) THEN
             
             
             
                   DUM=GAMMAR*VRAIN(INDEXR)
                   DUM1=DUM-VSNOW
             
             
             
             
             
                   DUM2=SQRT(DUM1*DUM1+.04*DUM*VSNOW)
                   DUM1=5.E-12*INDEXR*INDEXR+2.E-12*INDEXR*INDEXS        &
      &                 +.5E-12*INDEXS*INDEXS
                   FIR=MIN(1., CIACR*NLICE*DUM1*DUM2)
             
             
             
                   PIACR=MIN(PIACR+FIR*QTRAIN, QTRAIN)
                 ENDIF        
                DUM=PREVP-PIACR
                If (DUM .LT. PRLOSS) THEN
                  DUM1=PRLOSS/DUM
                  PREVP=DUM1*PREVP
                  PIACR=DUM1*PIACR
                ENDIF        
              ENDIF          
            ENDIF            
          ENDIF  do_rain     

          INDEXR1=INDEXR












          DUM1=PIACW+PRAUT+PRACW-MIN(0.,PCOND)
          IF (DUM1 .GT. QW) THEN
            DUM=QW/DUM1
            PIACW=DUM*PIACW
            PIACWI=DUM*PIACWI
            PRAUT=DUM*PRAUT
            PRACW=DUM*PRACW
            IF (PCOND .LT. 0.) PCOND=DUM*PCOND
          ENDIF
          PIACWR=PIACW-PIACWI          



          DELW=PCOND-PIACW-PRAUT-PRACW
          QWnew=QW+DELW
          IF (QWnew .LE. EPSQ) QWnew=0.
          IF (QW.GT.0. .AND. QWnew.NE.0.) THEN
            DUM=QWnew/QW
            IF (DUM .LT. TOLER) QWnew=0.
          ENDIF



          DELT= XLV1*(PCOND+PIEVP+PICND+PREVP)                          &
     &         +XLS1*PIDEP+XLF1*(PIACWI+PIACR-PIMLT)
          Tnew=TK+DELT

          DELV=-PCOND-PIDEP-PIEVP-PICND-PREVP
          WVnew=WV+DELV























          DELI=0.
          RimeF=1.
          IF (ICE_logical) THEN
            DELI=PIDEP+PIEVP+PIACWI+PIACR-PIMLT
            TOT_ICEnew=TOT_ICE+THICK*DELI
            IF (TOT_ICE.GT.0. .AND. TOT_ICEnew.NE.0.) THEN
              DUM=TOT_ICEnew/TOT_ICE
              IF (DUM .LT. TOLER) TOT_ICEnew=0.
            ENDIF
            IF (TOT_ICEnew .LE. CLIMIT) THEN
              TOT_ICEnew=0.
              RimeF=1.
              QInew=0.
              ASNOWnew=0.
            ELSE




              IF (PINIT<=EPSQ) THEN
                PSDEP=MAX(0., PIDEP)
              ELSE

                PSDEP=0.
              ENDIF
              DELS=PIACWI+PIACR+PSDEP
              IF (DELS<=EPSQ .OR. Tnew>=T0C) THEN
                RimeF=RimeF1
              ELSE
















































                DUM=1.E-6*REAL(INDEXS)              
                DUM1=DUM*DUM*DUM/MASSI(INDEXS)      
                RFmx=RFmx1*DUM1                     
                IF (PIACWI>0. .AND. Rcw>0.) THEN
                  DUM=-Rcw*VSNOW/MIN(-0.5,TC)
                  DUM=MIN(12.14, MAX(0.275, DUM) )
                  DUM=300.*DUM**0.44                
                  DUM=MIN(900., MAX(170., DUM) )    
                  RFrime=PI*DUM*DUM1                
                ELSE
                  RFrime=1.                         
                                                    
                ENDIF

                RFnew=(PSDEP+RFrime*PIACWI+RFmx*PIACR)/DELS
                DUM2=TOT_ICE+THICK*DELS
                DUM3=RimeF1*TOT_ICE+RFnew*THICK*DELS
                IF (DUM2<=CLIMIT) THEN
                  RimeF=RimeF1
                ELSE
                  RimeF=MIN(RFmx, MAX(1., DUM3/DUM2) )
                ENDIF

              ENDIF       

              DUM1=BLDTRH*(FLIMASS*VSNOW+(1.-FLIMASS)*VCI)
              QInew=TOT_ICEnew/(THICK+DUM1)
              IF (QInew .LE. EPSQ) QInew=0.
              IF (QI.GT.0. .AND. QInew.NE.0.) THEN
                DUM=QInew/QI
                IF (DUM .LT. TOLER) QInew=0.
              ENDIF
              ASNOWnew=QInew*DUM1
              IF (ASNOW.GT.0. .AND. ASNOWnew.NE.0.) THEN
                DUM=ASNOWnew/ASNOW
                IF (DUM .LT. TOLER) ASNOWnew=0.
              ENDIF
              RQSnew=FLIMASS*RHO*QInew
              IF (RQSnew>0.) THEN    

                Nsnow=RQSnew/(RimeF*MASSI(INDEXS))
                IF (RQSnew>0.0025 .AND. RimeF>5.) THEN

                  IF (RQSnew<=0.005) THEN
                    DUM=MIN(1., 400.*RQSnew-1.)       
                    Nsnow=1.E3*DUM+Nsnow*(1.-DUM)     
                  ELSE




                    DUM=180.*(RQSnew-0.005)           
                    Nsnow=1.E3*(1.-MIN(DUM,0.8))      
                  ENDIF
                ENDIF
                Zsnow=Cdry*RQSnew*RQSnew/Nsnow
              ENDIF   
            ENDIF         
          ENDIF           














          TIMLT=TIMLT+PIMLT*THICK     
          DELR=PRAUT+PRACW+PIACWR-PIACR+PIMLT+PREVP+PICND
          TOT_RAINnew=TOT_RAIN+THICK*DELR
          IF (TOT_RAIN.GT.0. .AND. TOT_RAINnew.NE.0.) THEN
            DUM=TOT_RAINnew/TOT_RAIN
            IF (DUM .LT. TOLER) TOT_RAINnew=0.
          ENDIF
update_rn: IF (TOT_RAINnew .LE. CLIMIT) THEN
            TOT_RAINnew=0.
            VRAIN2=0.
            QRnew=0.
            ARAINnew=0.
          ELSE  update_rn

            QRnew=TOT_RAINnew/(THICK+BLDTRH*VRAIN2)
            RQRnew=RHO*QRnew


            IF(STRAT .AND. RQRnew>1.E-3) STRAT=.FALSE.  
            IF(DRZL .AND. TIMLT>EPSQ) DRZL=.FALSE.      



            INDEXR2=INDEXR1
rain_pass:  DO IPASS=1,3
DSD2:         IF (RQRnew<=RQR_DRmin) THEN

                INDEXR=MDRmin
                N0r=RQRnew/MASSR(INDEXR)
              ELSE IF (DRZL) THEN  DSD2

                DUM=MAX(1.0, 0.5E-3/RQRnew)
                N0r=MIN(1.E9, N0r0*DUM*SQRT(DUM) )   
                DUM1=RQRnew/(PI*RHOL*N0r)
                DUM=1.E6*SQRT(SQRT(DUM1))
                INDEXR=MAX(XMRmin, MIN(DUM, XMRmax) )
                N0r=RQRnew/MASSR(INDEXR)
              ELSE IF (RQRnew>=RQR_DRmax) THEN

                INDEXR=MDRmax
                N0r=RQRnew/MASSR(INDEXR)
              ELSE DSD2

                N0r=N0r0
                DUM=CN0r0*SQRT(SQRT(RQRnew))
                INDEXR=MAX(XMRmin, MIN(DUM, XMRmax) )
                IF (STRAT .AND. INDEXR<INDEXRmin) THEN

                  INDEXR=INDEXRmin
                  N0r=RQRnew/MASSR(INDEXR)
                ENDIF
              ENDIF  DSD2
              VRAIN2=GAMMAR*VRAIN(INDEXR)
              QRnew=TOT_RAINnew/(THICK+BLDTRH*VRAIN2)
              RQRnew=RHO*QRnew
              IDR=ABS(INDEXR-INDEXR2)
              INDEXR2=INDEXR   
              IF(IDR<=5) EXIT  rain_pass    

            ENDDO  rain_pass

            IF (QRnew .LE. EPSQ) QRnew=0.
            IF (QR.GT.0. .AND. QRnew.NE.0.) THEN
              DUM=QRnew/QR
              IF (DUM .LT. TOLER) QRnew=0.
            ENDIF
            ARAINnew=BLDTRH*VRAIN2*QRnew
            IF (ARAIN.GT.0. .AND. ARAINnew.NE.0.) THEN
              DUM=ARAINnew/ARAIN
              IF (DUM .LT. TOLER) ARAINnew=0.
            ENDIF
            RQRnew=RHO*QRnew

            IF (RQRnew>EPSQ) THEN

              Nrain=N0r*1.E-6*REAL(INDEXR2)
              Zrain=Crain*RQRnew*RQRnew/Nrain   
            ENDIF

          ENDIF  update_rn

          WCnew=QWnew+QRnew+QInew







          QT=THICK*(WV+WC)+ARAIN+ASNOW
          QTnew=THICK*(WVnew+WCnew)+ARAINnew+ASNOWnew
          BUDGET=QT-QTnew



          DBG_logical=.FALSE.
          DUM=ABS(BUDGET)
          IF (DUM .GT. TOLER) THEN
            DUM=DUM/MIN(QT, QTnew)
            IF (DUM .GT. TOLER) DBG_logical=.TRUE.
          ENDIF

          IF (PRINT_diag) THEN
            ESW=MIN(1000.*FPVS0(Tnew),0.99*PP)
            QSWnew=(RHgrd+0.001)*EPS*ESW/(PP-ESW)
            IF( (QWnew>EPSQ .OR. QRnew>EPSQ .OR. WVnew>QSWnew)          &
     &         .AND. TC<T_ICE) DBG_logical=.TRUE.
          ENDIF

          IF ((WVnew.LT.EPSQ .OR. DBG_logical) .AND. PRINT_diag) THEN
   
            WRITE(0,"(/2(a,i4),2(a,i2))") '{} i=',I_index,' j=',J_index,&
     &                                    ' L=',L,' LSFC=',LSFC
   
            ESW=MIN(1000.*FPVS0(Tnew),0.99*PP)
            QSWnew=EPS*ESW/(PP-ESW)
            IF (TC.LT.0. .OR. Tnew .LT. 0.) THEN
              ESI=MIN(1000.*FPVS(Tnew),0.99*PP)
              QSInew=EPS*ESI/(PP-ESI)
            ELSE
              QSI=QSW
              QSInew=QSWnew
            ENDIF
            WSnew=QSInew
            WRITE(0,"(4(a12,g11.4,1x))")                                   &
     & '{} TCold=',TC,'TCnew=',Tnew-T0C,'P=',.01*PP,'RHO=',RHO,            &
     & '{} THICK=',THICK,'RHold=',WV/WS,'RHnew=',WVnew/WSnew,              &
     &   'RHgrd=',RHgrd,                                                   &
     & '{} RHWold=',WV/QSW,'RHWnew=',WVnew/QSWnew,'RHIold=',WV/QSI,        &
     &   'RHInew=',WVnew/QSInew,                                           &
     & '{} QSWold=',QSW,'QSWnew=',QSWnew,'QSIold=',QSI,'QSInew=',QSInew,   &
     & '{} WSold=',WS,'WSnew=',WSnew,'WVold=',WV,'WVnew=',WVnew,           &
     & '{} WCold=',WC,'WCnew=',WCnew,'QWold=',QW,'QWnew=',QWnew,           &
     & '{} QIold=',QI,'QInew=',QInew,'QRold=',QR,'QRnew=',QRnew,           &
     & '{} ARAINold=',ARAIN,'ARAINnew=',ARAINnew,'ASNOWold=',ASNOW,        &
     &   'ASNOWnew=',ASNOWnew,                                             &
     & '{} TOT_RAIN=',TOT_RAIN,'TOT_RAINnew=',TOT_RAINnew,                 &
     &   'TOT_ICE=',TOT_ICE,'TOT_ICEnew=',TOT_ICEnew,                      &
     & '{} BUDGET=',BUDGET,'QTold=',QT,'QTnew=',QTnew                       
   
            WRITE(0,"(4(a12,g11.4,1x))")                                   &
     & '{} DELT=',DELT,'DELV=',DELV,'DELW=',DELW,'DELI=',DELI,             &
     & '{} DELR=',DELR,'PCOND=',PCOND,'PIDEP=',PIDEP,'PIEVP=',PIEVP,       &
     & '{} PICND=',PICND,'PREVP=',PREVP,'PRAUT=',PRAUT,'PRACW=',PRACW,     &
     & '{} PIACW=',PIACW,'PIACWI=',PIACWI,'PIACWR=',PIACWR,'PIMLT=',       &
     &    PIMLT,                                                           &
     & '{} PIACR=',PIACR                                                    
   
            WRITE(0,"(4(a13,L2))")                                         &
     & '{} ICE_logical=',ICE_logical,'RAIN_logical=',RAIN_logical,         &
     &    'STRAT=',STRAT,'DRZL=',DRZL
   
            IF (ICE_logical) WRITE(0,"(4(a12,g11.4,1x))")                  &
     & '{} RimeF1=',RimeF1,'GAMMAS=',GAMMAS,'VrimeF=',VrimeF,              &
     &   'VSNOW=',VSNOW,                                                   &
     & '{} QSmICE=',QSmICE,'XLIMASS=',XLIMASS,'RQLICE=',RQLICE,            &
     &   'QTICE=',QTICE,                                                   &
     & '{} NLICE=',NLICE,'NSmICE=',NSmICE,'PILOSS=',PILOSS,                &
     &   'EMAIRI=',EMAIRI,                                                 &
     & '{} RimeF=',RimeF,'VCI=',VCI,'INDEXS=',FLOAT(INDEXS),               &
     &    'FLIMASS=',FLIMASS,                                              &
     & '{} Nsnow=',Nsnow,'Zsnow=',Zsnow,'TIMLT=',TIMLT,'RHOX0C=',RHOX0C,   &
     & '{} INDEXS0C=',FLOAT(INDEXS0C)
   
            IF (RAIN_logical) WRITE(0,"(4(a12,g11.4,1x))")                 &
     & '{} INDEXR1=',FLOAT(INDEXR1),'INDEXR=',FLOAT(INDEXR),               &
     &   'GAMMAR=',GAMMAR,'N0r=',N0r,                                      &
     & '{} VRAIN1=',VRAIN1,'VRAIN2=',VRAIN2,'QTRAIN=',QTRAIN,'RQR=',RQR,   &
     & '{} PRLOSS=',PRLOSS,'VOLR1=',THICK+BLDTRH*VRAIN1,                   &
     &   'VOLR2=',THICK+BLDTRH*VRAIN2,'INDEXR2=',INDEXR2,                  &
     & '{} Nrain=',Nrain,'Zrain=',Zrain
   
            IF (PRACW .GT. 0.) WRITE(0,"(a12,g11.4,1x)") '{} FWR=',FWR
   
            IF (PIACR .GT. 0.) WRITE(0,"(a12,g11.4,1x)") '{} FIR=',FIR
   
            DUM=PIMLT+PICND-PREVP-PIEVP
            IF (DUM.GT.0. .or. DWVi.NE.0.)                                 &
     &        WRITE(0,"(4(a12,g11.4,1x))")                                 &
     & '{} TFACTOR=',TFACTOR,'DYNVIS=',DYNVIS,                             &
     &   'THERM_CON=',THERM_COND,'DIFFUS=',DIFFUS
   
            IF (PREVP .LT. 0.) WRITE(0,"(4(a12,g11.4,1x))")                &
     & '{} RFACTOR=',RFACTOR,'ABW=',ABW,'VENTR=',VENTR,'CREVP=',CREVP,     &
     & '{} DWVr=',DWVr,'DENOMW=',DENOMW
   
            IF (PIDEP.NE.0. .AND. DWVi.NE.0.)                              &
     &        WRITE(0,"(4(a12,g11.4,1x))")                                 &
     & '{} DWVi=',DWVi,'DENOMI=',DENOMI,'PIDEP_max=',PIDEP_max,            &
     &   'SFACTOR=',SFACTOR,                                               &
     & '{} ABI=',ABI,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),           &
     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                                &
     & '{} VENTIS=',VENTIS,'DIDEP=',DIDEP
   
            IF (PIDEP.GT.0. .AND. PCOND.NE.0.)                             &
     &        WRITE(0,"(4(a12,g11.4,1x))")                                 &
     & '{} DENOMW=',DENOMW,'DENOMWI=',DENOMWI,'DENOMF=',DENOMF,            &
     &    'DUM2=',PCOND-PIACW
   
            IF (FWS .GT. 0.) WRITE(0,"(4(a12,g11.4,1x))")                  &
     & '{} FWS=',FWS                     
   
            DUM=PIMLT+PICND-PIEVP
            IF (DUM.GT. 0.) WRITE(0,"(4(a12,g11.4,1x))")                   &
     & '{} SFACTOR=',SFACTOR,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),   &
     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                                &
     & '{} AIEVP=',AIEVP,'DIEVP=',DIEVP,'QSW0=',QSW0,'DWV0=',DWV0       
   
          ENDIF





          T_col(L)=Tnew                           
          Q_col(L)=max(EPSQ, WVnew/(1.+WVnew))    
          WC_col(L)=max(EPSQ, WCnew)              
          QI_col(L)=max(EPSQ, QInew)              
          QR_col(L)=max(EPSQ, QRnew)              
          QW_col(L)=max(EPSQ, QWnew)              
          RimeF_col(L)=RimeF                      
          NR_col(L)=Nrain                         
          NS_col(L)=Nsnow                         



          Ztot=MAX(Ztot, Zrain+Zsnow)
          Ztot=Zrain+Zsnow
          IF (Ztot>Zmin) DBZ_col(L)=10.*ALOG10(Ztot)     
          ASNOW=ASNOWnew                          
          VSNOW1=VSNOW                            
          ARAIN=ARAINnew                          
          VRAIN1=VRAIN2                           



          if(pcond.gt.0)then
            pcond1d(L)=pcond
          elseif(pcond.lt.0)then
            pevap1d(L)=pcond
          endif
          if(pidep.gt.0)then
            pidep1d(L)=pidep
          elseif(pidep.lt.0)then
            pisub1d(L)=pidep
          endif
          piacw1d(L)=piacw
          piacwi1d(L)=piacwi
          piacwr1d(L)=piacwr
          piacr1d(L)=piacr
          picnd1d(L)=picnd
          pievp1d(L)=pievp
          pimlt1d(L)=pimlt
          praut1d(L)=praut
          pracw1d(L)=pracw
          prevp1d(L)=prevp
          if (qinew>EPSQ) then
            vsnow1d(L)=vsnow

            if (FLIMASS<1.) then 
              vci1d(L)=vci
              NSmICE1d(L)=NSmICE
            endif

            INDEXS1d(L)=FLOAT(INDEXS)
          endif
          if (qrnew>EPSQ) then
            vrain11d(L)=vrain1
            vrain21d(L)=vrain2
            INDEXR1d(L)=FLOAT(INDEXR2)

            IDR=0
            IF(STRAT) IDR=1
            IF(DRZL) IDR=IDR+2
            RFlag1d(L)=IDR

          endif



      ENDDO  big_loop







        CONTAINS




      REAL FUNCTION CONDENSE (PP,QW,TK,WV,RHgrd,I,J,L)







      IMPLICIT NONE

      INTEGER, PARAMETER :: HIGH_PRES=Selected_Real_Kind(15)
      REAL (KIND=HIGH_PRES), PARAMETER ::                               &
     & RHLIMIT=.001, RHLIMIT1=-RHLIMIT
      REAL (KIND=HIGH_PRES) :: COND, SSAT, WCdum

      REAL,INTENT(IN) :: QW,PP,WV,TK,RHgrd
      REAL WVdum,Tdum,DWV,WS,ESW

integer,INTENT(IN) :: I,J,L     
integer :: niter
real :: DWVinp,SSATinp




      Tdum=TK
      WVdum=WV
      WCdum=QW
      ESW=MIN(1000.*FPVS0(Tdum),0.99*PP)        
      WS=RHgrd*EPS*ESW/(PP-ESW)                 
      DWV=WVdum-WS                              
      SSAT=DWV/WS                               
      CONDENSE=0.

DWVinp=DWV     
SSATinp=SSAT

      DO NITER=1,MAX_ITERATIONS
        COND=DWV/(1.+XLV3*WS/(Tdum*Tdum))       
        COND=MAX(COND, -WCdum)                  
        Tdum=Tdum+XLV1*COND                     
        WVdum=WVdum-COND                        
        WCdum=WCdum+COND                        
        CONDENSE=CONDENSE+COND                  
        ESW=MIN(1000.*FPVS0(Tdum),0.99*PP)      
        WS=RHgrd*EPS*ESW/(PP-ESW)               
        DWV=WVdum-WS                            
        SSAT=DWV/WS                             
        IF (SSAT>=RHLIMIT1 .AND. SSAT<=RHLIMIT) EXIT   
        IF (SSAT<RHLIMIT1 .AND. WCdum<=EPSQ) EXIT      
      ENDDO


IF (NITER>MAX_ITERATIONS) THEN

   IF (WARN3) THEN
      write(0,*) 'WARN3: Too many iterations in function CONDENSE; ', &
         ' I,J,L,TC,SSAT,QW,DWV=',I,J,L,TK-T0C,SSATinp,1000.*QW,DWVinp
      WARN3=.FALSE.
   ENDIF
   SSAT=0.
   CONDENSE=DWVinp
ENDIF


      END FUNCTION CONDENSE





      REAL FUNCTION DEPOSIT (PP,Tdum,WVdum,RHgrd,I,J,L)   




      IMPLICIT NONE      

      INTEGER, PARAMETER :: HIGH_PRES=Selected_Real_Kind(15)
      REAL (KIND=HIGH_PRES), PARAMETER :: RHLIMIT=.001,                 &
     & RHLIMIT1=-RHLIMIT
      REAL (KIND=HIGH_PRES) :: DEP, SSAT

      real,INTENT(IN) ::  PP,RHgrd
      real,INTENT(INOUT) ::  WVdum,Tdum
      real ESI,WS,DWV

integer,INTENT(IN) :: I,J,L   
integer :: niter
real :: Tinp,DWVinp,SSATinp




      ESI=MIN(1000.*FPVS(Tdum),0.99*PP)         
      WS=RHgrd*EPS*ESI/(PP-ESI)                 
      DWV=WVdum-WS                              
      SSAT=DWV/WS                               
      DEPOSIT=0.

Tinp=Tdum     
DWVinp=DWV
SSATinp=SSAT

      DO NITER=1,MAX_ITERATIONS
        DEP=DWV/(1.+XLS3*WS/(Tdum*Tdum))        
        Tdum=Tdum+XLS1*DEP                      
        WVdum=WVdum-DEP                         
        DEPOSIT=DEPOSIT+DEP                     
        ESI=MIN(1000.*FPVS(Tdum),0.99*PP)       
        WS=RHgrd*EPS*ESI/(PP-ESI)               
        DWV=WVdum-WS                            
        SSAT=DWV/WS                             
        IF (SSAT>=RHLIMIT1 .AND. SSAT<=RHLIMIT) EXIT   
      ENDDO


IF (NITER>MAX_ITERATIONS) THEN

   IF (WARN2) THEN
      write(0,*) 'WARN2: Too many iterations in function DEPOSIT; ', &
         ' I,J,L,TC,SSAT,DWV=',I,J,L,Tinp-T0C,SSATinp,DWVinp
      WARN2=.FALSE.
   ENDIF
   SSAT=0.
   DEPOSIT=DWVinp
ENDIF


      END FUNCTION DEPOSIT







      INTEGER FUNCTION GET_INDEXR(RR)
      IMPLICIT NONE
      REAL, INTENT(IN) :: RR
      IF (RR .LE. RR_DRmin) THEN




        GET_INDEXR=MDRmin
      ELSE IF (RR .LE. RR_DR1) THEN








        GET_INDEXR=INT( 1.123E3*RR**.1947 + .5 )
        GET_INDEXR=MAX( MDRmin, MIN(GET_INDEXR, MDR1) )
      ELSE IF (RR .LE. RR_DR2) THEN








        GET_INDEXR=INT( 1.225E3*RR**.2017 + .5 )
        GET_INDEXR=MAX( MDR1, MIN(GET_INDEXR, MDR2) )
      ELSE IF (RR .LE. RR_DR3) THEN








        GET_INDEXR=INT( 1.3006E3*RR**.2083 + .5 )
        GET_INDEXR=MAX( MDR2, MIN(GET_INDEXR, MDR3) )
      ELSE IF (RR .LE. RR_DR4) THEN








        GET_INDEXR=INT( 1.354E3*RR**.2143 + .5 )
        GET_INDEXR=MAX( MDR3, MIN(GET_INDEXR, MDR4) )
      ELSE IF (RR .LE. RR_DR5) THEN








        GET_INDEXR=INT( 1.404E3*RR**.2213 + .5 )
        GET_INDEXR=MAX( MDR4, MIN(GET_INDEXR, MDR5) )
      ELSE IF (RR .LE. RR_DRmax) THEN








        GET_INDEXR=INT( 1.4457E3*RR**.2303 + .5 )
        GET_INDEXR=MAX( MDR5, MIN(GET_INDEXR, MDRmax) )
      ELSE 




        GET_INDEXR=MDRmax
      ENDIF              

      END FUNCTION GET_INDEXR

      END SUBROUTINE EGCP01COLUMN









      SUBROUTINE fer_hires_init (GSMDT,DT,DELX,DELY,LOWLYR,restart,         &

     &   ALLOWED_TO_READ,                                               &
     &   IDS,IDE,JDS,JDE,KDS,KDE,                                       &
     &   IMS,IME,JMS,JME,KMS,KME,                                       &
     &   ITS,ITE,JTS,JTE,KTS,KTE,                                       &
     &   F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY)

























































      IMPLICIT NONE















      INTEGER, PARAMETER :: MDR1=XMR1, MDR2=XMR2, MDR3=XMR3


      integer,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     & 
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE       

       INTEGER, DIMENSION(ims:ime,jms:jme),INTENT(INOUT) :: LOWLYR

      real, INTENT(IN) ::  DELX,DELY


      real,DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(OUT),OPTIONAL :: &
     &  F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY
      real,INTENT(IN) :: DT,GSMDT
      LOGICAL,INTENT(IN) :: allowed_to_read,restart




      REAL :: BBFR,DTPH,DX,Thour_print,RDIS
      INTEGER :: I,IM,J,L,K,JTF,KTF,ITF
      INTEGER :: etampnew_unit1
      LOGICAL :: opened
      LOGICAL , EXTERNAL      :: wrf_dm_on_monitor
      CHARACTER*80 errmess



      JTF=MIN0(JTE,JDE-1)
      KTF=MIN0(KTE,KDE-1)
      ITF=MIN0(ITE,IDE-1)

      DO J=JTS,JTF
      DO I=ITS,ITF
        LOWLYR(I,J)=1
      ENDDO
      ENDDO

      IF(.NOT.RESTART .AND. ALLOWED_TO_READ .AND. present(F_ICE_PHY)) THEN    
        CALL wrf_debug(1,'WARNING: F_ICE_PHY,F_RAIN_PHY AND F_RIMEF_PHY IS REINITIALIZED')   
        DO J = jts,jte
        DO K = kts,kte
        DO I= its,ite
          F_ICE_PHY(i,k,j)=0.
          F_RAIN_PHY(i,k,j)=0.
          F_RIMEF_PHY(i,k,j)=1.
        ENDDO
        ENDDO
        ENDDO
      ENDIF


      IF(ALLOWED_TO_READ)THEN


        DX=SQRT((DELX)**2+(DELY)**2)/1000.    
        DX=MIN(100., MAX(5., DX) )






        DTPH=MAX(GSMDT*60.,DT)
        DTPH=NINT(DTPH/DT)*DT



        CALL GPVS



        IF ( wrf_dm_on_monitor() ) THEN
          DO i = 31,99
            INQUIRE ( i , OPENED = opened )
            IF ( .NOT. opened ) THEN
              etampnew_unit1 = i
              GOTO 2061
            ENDIF
          ENDDO
          etampnew_unit1 = -1
 2061     CONTINUE
        ENDIF

        CALL wrf_dm_bcast_bytes ( etampnew_unit1 , 4 )

        IF ( etampnew_unit1 < 0 ) THEN
          CALL wrf_error_fatal3("<stdin>",2961,&
'module_mp_fer_hires: fer_hires_init: Can not find unused fortran unit to read in lookup table.' )
        ENDIF

        IF ( wrf_dm_on_monitor() ) THEN



          print*,'open ETAMPNEW_DATA.expanded_rain, in fer_hires' 
          OPEN(UNIT=etampnew_unit1,FILE="ETAMPNEW_DATA.expanded_rain",  &
     &        FORM="UNFORMATTED",STATUS="OLD",ERR=9061)

          READ(etampnew_unit1) VENTR1
          READ(etampnew_unit1) VENTR2
          READ(etampnew_unit1) ACCRR
          READ(etampnew_unit1) MASSR
          READ(etampnew_unit1) VRAIN
          READ(etampnew_unit1) RRATE
          READ(etampnew_unit1) VENTI1
          READ(etampnew_unit1) VENTI2
          READ(etampnew_unit1) ACCRI
          READ(etampnew_unit1) MASSI
          READ(etampnew_unit1) VSNOWI
          READ(etampnew_unit1) VEL_RF

          CLOSE (etampnew_unit1)
        ENDIF

        CALL wrf_dm_bcast_bytes ( VENTR1 , size ( VENTR1 ) * 4 )
        CALL wrf_dm_bcast_bytes ( VENTR2 , size ( VENTR2 ) * 4 )
        CALL wrf_dm_bcast_bytes ( ACCRR , size ( ACCRR ) * 4 )
        CALL wrf_dm_bcast_bytes ( MASSR , size ( MASSR ) * 4 )
        CALL wrf_dm_bcast_bytes ( VRAIN , size ( VRAIN ) * 4 )
        CALL wrf_dm_bcast_bytes ( RRATE , size ( RRATE ) * 4 )
        CALL wrf_dm_bcast_bytes ( VENTI1 , size ( VENTI1 ) * 4 )
        CALL wrf_dm_bcast_bytes ( VENTI2 , size ( VENTI2 ) * 4 )
        CALL wrf_dm_bcast_bytes ( ACCRI , size ( ACCRI ) * 4 )
        CALL wrf_dm_bcast_bytes ( MASSI , size ( MASSI ) * 4 )
        CALL wrf_dm_bcast_bytes ( VSNOWI , size ( VSNOWI ) * 4 )
        CALL wrf_dm_bcast_bytes ( VEL_RF , size ( VEL_RF ) * 4 )




        CALL MY_GROWTH_RATES (DTPH)

        PI_E=ACOS(-1.)




        ABFR=-0.66
        BBFR=100.
        CBFR=20.*PI_E*PI_E*BBFR*RHOL*1.E-21*DTPH   





        CIACW=0.5*DTPH*0.25*PI_E  




        CIACR=PI_E*DTPH







        RR_DRmin=N0r0*RRATE(MDRmin)     
        RR_DR1=N0r0*RRATE(MDR1)         
        RR_DR2=N0r0*RRATE(MDR2)         
        RR_DR3=N0r0*RRATE(MDR3)         
        RR_DR4=N0r0*RRATE(MDR4)         
        RR_DR5=N0r0*RRATE(MDR5)         
        RR_DRmax=N0r0*RRATE(MDRmax)     

        RQR_DRmin=N0r0*MASSR(MDRmin)    
        RQR_DRmax=N0r0*MASSR(MDRmax)    
        C_NR=1./(PI*RHOL)   
        Crain=720.E18*C_NR*C_NR    
        C_N0r0=PI_E*RHOL*N0r0
        CN0r0=1.E6/SQRT(SQRT(C_N0r0))
        CN0r_DMRmin=1./(PI_E*RHOL*DMRmin*DMRmin*DMRmin*DMRmin)
        CN0r_DMRmax=1./(PI_E*RHOL*DMRmax*DMRmax*DMRmax*DMRmax)




        CRACW=DTPH*0.25*PI_E*1.0

        ESW0=1000.*FPVS0(T0C)     
        RFmax=1.1**Nrime          
        RFmx1=PI*900.             










        RDIS=0.5  
        BETA6=( (1.+3.*RDIS*RDIS)*(1.+4.*RDIS*RDIS)*(1.+5.*RDIS*RDIS)/  &
     &         ((1.+RDIS*RDIS)*(1.+2.*RDIS*RDIS) ) )





        ARAUT=1.03e19/(NCW*SQRT(NCW))
        BRAUT=DTPH*1.1E10*BETA6/NCW









        QAUT0=PI*RHOL*NCW*(20.E-6)**3/6.     



        ARcw=1.E6*(0.75/(PI*NCW*RHOL))**C1










        RH_NgC=PI*500.*1.E3                
        RQhail=RH_NgC*(1.E-3)**3           
        Bvhail=0.82*C1                     
        Avhail=1353.*(1./RH_NgC)**Bvhail   
        RH_NgC=1.E6*(1./RH_NgC)**C1        








        RH_NgT=PI*300.*5.E3                
        RH_NgT=1.E6*(1./RH_NgT)**C1        














        DO I=MDImin,MDImax
          SDENS(I)=PI_E*1.5E-15*FLOAT(I*I*I)/MASSI(I)
        ENDDO

        Thour_print=-DTPH/3600.


      ENDIF  

      RETURN



9061 CONTINUE
      WRITE( errmess , '(A,I4)' )                                        &
       'module_mp_hwrf: error opening ETAMPNEW_DATA on unit '          &
     &, etampnew_unit1
      CALL wrf_error_fatal3("<stdin>",3148,&
errmess)


      END SUBROUTINE fer_hires_init

      SUBROUTINE MY_GROWTH_RATES (DTPH)












      IMPLICIT NONE

      REAL,INTENT(IN) :: DTPH

      REAL  DT_ICE
      REAL,DIMENSION(35) :: MY_600



      DATA MY_600 /                                                     &
     & 5.5e-8, 1.4E-7, 2.8E-7, 6.E-7, 3.3E-6,                           & 
     & 2.E-6, 9.E-7, 8.8E-7, 8.2E-7, 9.4e-7,                            & 
     & 1.2E-6, 1.85E-6, 5.5E-6, 1.5E-5, 1.7E-5,                         & 
     & 1.5E-5, 1.E-5, 3.4E-6, 1.85E-6, 1.35E-6,                         & 
     & 1.05E-6, 1.E-6, 9.5E-7, 9.0E-7, 9.5E-7,                          & 
     & 9.5E-7, 9.E-7, 9.E-7, 9.E-7, 9.E-7,                              & 
     & 9.E-7, 9.E-7, 9.E-7, 9.E-7, 9.E-7 /        



      DT_ICE=(DTPH/600.)**1.5
      MY_GROWTH=DT_ICE*MY_600*1.E-3    



      END SUBROUTINE MY_GROWTH_RATES





      SUBROUTINE GPVS






























      IMPLICIT NONE
      real :: X,XINC,T
      integer :: JX

      XINC=(XMAX-XMIN)/(NX-1)
      C1XPVS=1.-XMIN/XINC
      C2XPVS=1./XINC
      C1XPVS0=1.-XMIN/XINC
      C2XPVS0=1./XINC

      DO JX=1,NX
        X=XMIN+(JX-1)*XINC
        T=X
        TBPVS(JX)=FPVSX(T)
        TBPVS0(JX)=FPVSX0(T)
      ENDDO

      END SUBROUTINE GPVS



                     REAL   FUNCTION FPVS(T)































      IMPLICIT NONE
      real,INTENT(IN) :: T
      real XJ
      integer :: JX

      IF (T>=XMIN .AND. T<=XMAX) THEN
        XJ=MIN(MAX(C1XPVS+C2XPVS*T,1.),FLOAT(NX))
        JX=MIN(XJ,NX-1.)
        FPVS=TBPVS(JX)+(XJ-JX)*(TBPVS(JX+1)-TBPVS(JX))
      ELSE IF (T>XMAX) THEN


        FPVS=0.61078*exp(17.2694*(T-273.16)/(T-35.86))
      ELSE 


        FPVS=0.61078*exp(21.8746*(T-273.16)/(T-7.66))
      ENDIF

      END FUNCTION FPVS


                     REAL FUNCTION FPVS0(T)

      IMPLICIT NONE
      real,INTENT(IN) :: T
      real :: XJ1
      integer :: JX1

      IF (T>=XMIN .AND. T<=XMAX) THEN
        XJ1=MIN(MAX(C1XPVS0+C2XPVS0*T,1.),FLOAT(NX))
        JX1=MIN(XJ1,NX-1.)
        FPVS0=TBPVS0(JX1)+(XJ1-JX1)*(TBPVS0(JX1+1)-TBPVS0(JX1))
      ELSE


        FPVS0=0.61078*exp(17.2694*(T-273.16)/(T-35.86))
      ENDIF

      END FUNCTION FPVS0



                    REAL FUNCTION FPVSX(T)



































      IMPLICIT NONE

       real, parameter :: TTP=2.7316E+2,HVAP=2.5000E+6,PSAT=6.1078E+2   &
      ,         CLIQ=4.1855E+3,CVAP= 1.8460E+3                          &
      ,         CICE=2.1060E+3,HSUB=2.8340E+6

      real, parameter :: PSATK=PSAT*1.E-3
      real, parameter :: DLDT=CVAP-CLIQ,XA=-DLDT/RV,XB=XA+HVAP/(RV*TTP)
      real, parameter :: DLDTI=CVAP-CICE                                &
      ,                  XAI=-DLDTI/RV,XBI=XAI+HSUB/(RV*TTP)
      real T,TR

      TR=TTP/T

      IF(T.GE.TTP)THEN
        FPVSX=PSATK*(TR**XA)*EXP(XB*(1.-TR))
      ELSE
        FPVSX=PSATK*(TR**XAI)*EXP(XBI*(1.-TR))
      ENDIF

      END FUNCTION FPVSX


                 REAL   FUNCTION FPVSX0(T)

      IMPLICIT NONE
      real,parameter :: TTP=2.7316E+2,HVAP=2.5000E+6,PSAT=6.1078E+2     &
      ,         CLIQ=4.1855E+3,CVAP=1.8460E+3                           &
      ,         CICE=2.1060E+3,HSUB=2.8340E+6
      real,PARAMETER :: PSATK=PSAT*1.E-3
      real,PARAMETER :: DLDT=CVAP-CLIQ,XA=-DLDT/RV,XB=XA+HVAP/(RV*TTP)
      real,PARAMETER :: DLDTI=CVAP-CICE                                 &
      ,                 XAI=-DLDT/RV,XBI=XA+HSUB/(RV*TTP)
      real :: T,TR

      TR=TTP/T
      FPVSX0=PSATK*(TR**XA)*EXP(XB*(1.-TR))

      END FUNCTION FPVSX0

      END MODULE module_mp_fer_hires
