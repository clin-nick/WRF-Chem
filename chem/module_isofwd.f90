
















      SUBROUTINE ISRP1F2p1 (WI, RHI, TEMPI)
      INCLUDE 'module_isrpia_inc.F'
      DIMENSION WI(NCOMP)







      CALL INIT12p1 (WI, RHI, TEMPI)



      SULRAT = W(3)/W(2)





      IF (2.0.LE.SULRAT) THEN 
      DC   = W(3) - 2.001D0*W(2)  
      W(3) = W(3) + MAX(-DC, ZERO)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'A2'
         CALL CALCA22p1                 
      ELSE

         IF (RH.LT.DRNH42S4) THEN    
            SCASE = 'A1'
            CALL CALCA12p1              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'A2'
            CALL CALCA22p1              
         ENDIF
      ENDIF



      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 

      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL CALCB42p1                 
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'B1'
            CALL CALCB12p1              

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'B2'
            CALL CALCB22p1              

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'B3'
            CALL CALCB32p1              

         ELSEIF (DRNH42S4.LE.RH) THEN         
            SCASE = 'B4'
            CALL CALCB42p1              
         ENDIF
      ENDIF
      CALL CALCNH32p1



      ELSEIF (SULRAT.LT.1.0) THEN             

      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC22p1                 
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'C1'
            CALL CALCC12p1              

         ELSEIF (DRNH4HS4.LE.RH) THEN         
            SCASE = 'C2'
            CALL CALCC22p1              

         ENDIF
      ENDIF
      CALL CALCNH32p1
      ENDIF



      RETURN



      END
















      SUBROUTINE ISRP2F2p1 (WI, RHI, TEMPI)
      INCLUDE 'module_isrpia_inc.F'
      DIMENSION WI(NCOMP)







      CALL INIT22p1 (WI, RHI, TEMPI)



      SULRAT = W(3)/W(2)





      IF (2.0.LE.SULRAT) THEN                

      IF(METSTBL.EQ.1) THEN
         SCASE = 'D3'
         CALL CALCD32p1                 
      ELSE

         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'D1'
            CALL CALCD12p1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'D2'
            CALL CALCD22p1              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'D3'
            CALL CALCD32p1              
         ENDIF
      ENDIF







      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 

      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL CALCB42p1                 
         SCASE = 'E4'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'B1'
            CALL CALCB12p1              
            SCASE = 'E1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'B2'
            CALL CALCB22p1              
            SCASE = 'E2'

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'B3'
            CALL CALCB32p1              
            SCASE = 'E3'

         ELSEIF (DRNH42S4.LE.RH) THEN         
            SCASE = 'B4'
            CALL CALCB42p1              
            SCASE = 'E4'
         ENDIF
      ENDIF

      CALL CALCNA2p1                 







      ELSEIF (SULRAT.LT.1.0) THEN             

      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC22p1                 
         SCASE = 'F2'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'C1'
            CALL CALCC12p1              
            SCASE = 'F1'

         ELSEIF (DRNH4HS4.LE.RH) THEN         
            SCASE = 'C2'
            CALL CALCC22p1              
            SCASE = 'F2'
         ENDIF
      ENDIF

      CALL CALCNA2p1                 
      ENDIF



      RETURN



      END
















      SUBROUTINE ISRP3F2p1 (WI, RHI, TEMPI)
      INCLUDE 'module_isrpia_inc.F'
      DIMENSION WI(NCOMP)






      WI(3) = MAX (WI(3), 1.D-10)  
      WI(5) = MAX (WI(5), 1.D-10)  



      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
         WI(1) = 1.D-10  
         WI(2) = 1.D-10  
      ENDIF



      CALL ISOINIT32p1 (WI, RHI, TEMPI)



      REST = 2.D0*W(2) + W(4) + W(5) 
      IF (W(1).GT.REST) THEN            
         W(1) = (ONE-1D-6)*REST         
         CALL PUSHERR2p1 (0050, 'ISRP3F')  
      ENDIF



      SULRAT = (W(1)+W(3))/W(2)
      SODRAT = W(1)/W(2)





      IF (2.0.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                

      IF(METSTBL.EQ.1) THEN
         SCASE = 'G5'
         CALL CALCG52p1                 
      ELSE

         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'G1'
            CALL CALCG12p1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH4CL) THEN         
            SCASE = 'G2'
            CALL CALCG22p1              

         ELSEIF (DRNH4CL.LE.RH  .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'G3'
            CALL CALCG32p1              

        ELSEIF (DRNH42S4.LE.RH  .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'G4'
            CALL CALCG42p1              

         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'G5'
            CALL CALCG52p1              
         ENDIF
      ENDIF



      ELSE IF (SULRAT.GE.2.0 .AND. SODRAT.GE.2.0) THEN                

      IF(METSTBL.EQ.1) THEN
         SCASE = 'H6'
         CALL CALCH62p1                 
      ELSE

         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'H1'
            CALL CALCH12p1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN         
            SCASE = 'H2'
            CALL CALCH22p1              

         ELSEIF (DRNANO3.LE.RH  .AND. RH.LT.DRNACL) THEN         
            SCASE = 'H3'
            CALL CALCH32p1              

         ELSEIF (DRNACL.LE.RH   .AND. RH.LT.DRNH4Cl) THEN         
            SCASE = 'H4'
            CALL CALCH42p1              

         ELSEIF (DRNH4Cl.LE.RH .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'H5'
            CALL CALCH52p1              

         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'H6'
            CALL CALCH62p1              
         ENDIF
      ENDIF



      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 

      IF(METSTBL.EQ.1) THEN
         SCASE = 'I6'
         CALL CALCI62p1                 
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'I1'
            CALL CALCI12p1              

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN         
            SCASE = 'I2'
            CALL CALCI22p1              

         ELSEIF (DRNAHSO4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'I3'
            CALL CALCI32p1              

         ELSEIF (DRLC.LE.RH     .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'I4'
            CALL CALCI42p1              

         ELSEIF (DRNH42S4.LE.RH .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'I5'
            CALL CALCI52p1              

         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'I6'
            CALL CALCI62p1              
         ENDIF
      ENDIF

      CALL CALCNHA2p1                
      CALL CALCNH32p1                



      ELSEIF (SULRAT.LT.1.0) THEN             

      IF(METSTBL.EQ.1) THEN
         SCASE = 'J3'
         CALL CALCJ32p1                 
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'J1'
            CALL CALCJ12p1              

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN         
            SCASE = 'J2'
            CALL CALCJ22p1              

         ELSEIF (DRNAHSO4.LE.RH) THEN         
            SCASE = 'J3'
            CALL CALCJ32p1              
         ENDIF
      ENDIF

      CALL CALCNHA2p1                
      CALL CALCNH32p1                
      ENDIF



      RETURN



      END

















      SUBROUTINE ISRP4F2p1 (WI, RHI, TEMPI)
      INCLUDE 'module_isrpia_inc.F'
      DIMENSION WI(NCOMP)
      DOUBLE PRECISION NAFRI, NO3FRI















      CALL INIT42p1 (WI, RHI, TEMPI)



      REST = 2.D0*W(2) + W(4) + W(5)

      IF (W(1)+W(6)+W(7)+W(8).GT.REST) THEN

      CCASO4I  = MIN (W(2),W(6))
      FRSO4I   = MAX (W(2) - CCASO4I, ZERO)
      CAFRI    = MAX (W(6) - CCASO4I, ZERO)
      CCANO32I = MIN (CAFRI, 0.5D0*W(4))
      CAFRI    = MAX (CAFRI - CCANO32I, ZERO)
      NO3FRI   = MAX (W(4) - 2.D0*CCANO32I, ZERO)
      CCACL2I  = MIN (CAFRI, 0.5D0*W(5))
      CLFRI    = MAX (W(5) - 2.D0*CCACL2I, ZERO)
      REST1    = 2.D0*FRSO4I + NO3FRI + CLFRI

      CNA2SO4I = MIN (FRSO4I, 0.5D0*W(1))
      FRSO4I   = MAX (FRSO4I - CNA2SO4I, ZERO)
      NAFRI    = MAX (W(1) - 2.D0*CNA2SO4I, ZERO)
      CNACLI   = MIN (NAFRI, CLFRI)
      NAFRI    = MAX (NAFRI - CNACLI, ZERO)
      CLFRI    = MAX (CLFRI - CNACLI, ZERO)
      CNANO3I  = MIN (NAFRI, NO3FRI)
      NO3FR    = MAX (NO3FRI - CNANO3I, ZERO)
      REST2    = 2.D0*FRSO4I + NO3FRI + CLFRI

      CMGSO4I  = MIN (FRSO4I, W(8))
      FRMGI    = MAX (W(8) - CMGSO4I, ZERO)
      FRSO4I   = MAX (FRSO4I - CMGSO4I, ZERO)
      CMGNO32I = MIN (FRMGI, 0.5D0*NO3FRI)
      FRMGI    = MAX (FRMGI - CMGNO32I, ZERO)
      NO3FRI   = MAX (NO3FRI - 2.D0*CMGNO32I, ZERO)
      CMGCL2I  = MIN (FRMGI, 0.5D0*CLFRI)
      CLFRI    = MAX (CLFRI - 2.D0*CMGCL2I, ZERO)
      REST3    = 2.D0*FRSO4I + NO3FRI + CLFRI

         IF (W(6).GT.REST) THEN                       
             W(6) = (ONE-1D-6)*REST              
             W(1)= ZERO                          
             W(7)= ZERO                          
             W(8)= ZERO                          
             CALL PUSHERR2p1 (0051, 'ISRP4F')       

         ELSE IF (W(1).GT.REST1) THEN                 
             W(1) = (ONE-1D-6)*REST1             
             W(7)= ZERO                          
             W(8)= ZERO                          
             CALL PUSHERR2p1 (0052, 'ISRP4F')       

         ELSE IF (W(8).GT.REST2) THEN                 
             W(8) = (ONE-1D-6)*REST2             
             W(7)= ZERO                          
             CALL PUSHERR2p1 (0053, 'ISRP4F')       

         ELSE IF (W(7).GT.REST3) THEN                 
             W(7) = (ONE-1D-6)*REST3             
             CALL PUSHERR2p1 (0054, 'ISRP4F')       
         ENDIF
      ENDIF



      SO4RAT  = (W(1)+W(3)+W(6)+W(7)+W(8))/W(2)
      CRNARAT = (W(1)+W(6)+W(7)+W(8))/W(2)
      CRRAT   = (W(6)+W(7)+W(8))/W(2)





      IF (2.0.LE.SO4RAT .AND. CRNARAT.LT.2.0) THEN

       IF(METSTBL.EQ.1) THEN
         SCASE = 'O7'
         CALL CALCO72p1                 
       ELSE

         IF (RH.LT.DRNH4NO3) THEN
            SCASE = 'O1'
            CALL CALCO12p1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH4CL) THEN
            SCASE = 'O2'
            CALL CALCO22p1              

         ELSEIF (DRNH4CL.LE.RH  .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'O3'
            CALL CALCO32p1              

         ELSEIF (DRNH42S4.LE.RH .AND. RH.LT.DRMGSO4) THEN
            SCASE = 'O4'
            CALL CALCO42p1              

         ELSEIF (DRMGSO4.LE.RH .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'O5'
            CALL CALCO52p1              

         ELSEIF (DRNA2SO4.LE.RH .AND. RH.LT.DRK2SO4) THEN
            SCASE = 'O6'
            CALL CALCO62p1              

         ELSEIF (DRK2SO4.LE.RH) THEN
            SCASE = 'O7'
            CALL CALCO72p1              
         ENDIF
       ENDIF



      ELSEIF (SO4RAT.GE.2.0 .AND. CRNARAT.GE.2.0) THEN

       IF (CRRAT.LE.2.0) THEN

        IF(METSTBL.EQ.1) THEN
         SCASE = 'M8'
         CALL CALCM82p1                 
        ELSE

           IF (RH.LT.DRNH4NO3) THEN
             SCASE = 'M1'
             CALL CALCM12p1            

           ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN
             SCASE = 'M2'
             CALL CALCM22p1            

           ELSEIF (DRNANO3.LE.RH  .AND. RH.LT.DRNACL) THEN
             SCASE = 'M3'
             CALL CALCM32p1            

           ELSEIF (DRNACL.LE.RH   .AND. RH.LT.DRNH4Cl) THEN
             SCASE = 'M4'
             CALL CALCM42p1            

           ELSEIF (DRNH4Cl.LE.RH .AND. RH.LT.DRMGSO4) THEN
             SCASE = 'M5'
             CALL CALCM52p1            

           ELSEIF (DRMGSO4.LE.RH .AND. RH.LT.DRNA2SO4) THEN
             SCASE = 'M6'
             CALL CALCM62p1            

           ELSEIF (DRNA2SO4.LE.RH .AND. RH.LT.DRK2SO4) THEN
             SCASE = 'M7'
             CALL CALCM72p1            

           ELSEIF (DRK2SO4.LE.RH) THEN
             SCASE = 'M8'
             CALL CALCM82p1            
           ENDIF
        ENDIF




       ELSEIF (CRRAT.GT.2.0) THEN

        IF(METSTBL.EQ.1) THEN
         SCASE = 'P13'
         CALL CALCP132p1                 
        ELSE

           IF (RH.LT.DRCACL2) THEN
             SCASE = 'P1'
             CALL CALCP12p1             


           ELSEIF (DRCACL2.LE.RH .AND. RH.LT.DRMGCL2) THEN
             SCASE = 'P2'
             CALL CALCP22p1            


           ELSEIF (DRMGCL2.LE.RH  .AND. RH.LT.DRCANO32) THEN
             SCASE = 'P3'
             CALL CALCP32p1            


           ELSEIF (DRCANO32.LE.RH   .AND. RH.LT.DRMGNO32) THEN
             SCASE = 'P4'
             CALL CALCP42p1            


           ELSEIF (DRMGNO32.LE.RH .AND. RH.LT.DRNH4NO3) THEN
             SCASE = 'P5'
             CALL CALCP52p1            


           ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN
             SCASE = 'P6'
             CALL CALCP62p1            

           ELSEIF (DRNANO3.LE.RH .AND. RH.LT.DRNACL) THEN
             SCASE = 'P7'
             CALL CALCP72p1            

           ELSEIF (DRNACL.LE.RH .AND. RH.LT.DRNH4CL) THEN
             SCASE = 'P8'
             CALL CALCP82p1            

           ELSEIF (DRNH4CL.LE.RH .AND. RH.LT.DRKCL) THEN
             SCASE = 'P9'
             CALL CALCP92p1            

           ELSEIF (DRKCL.LE.RH .AND. RH.LT.DRMGSO4) THEN
             SCASE = 'P10'
             CALL CALCP102p1            

           ELSEIF (DRMGSO4.LE.RH .AND. RH.LT.DRKNO3) THEN
             SCASE = 'P11'
             CALL CALCP112p1            

           ELSEIF (DRKNO3.LE.RH .AND. RH.LT.DRK2SO4) THEN
             SCASE = 'P12'
             CALL CALCP122p1            

           ELSEIF (DRK2SO4.LE.RH) THEN
             SCASE = 'P13'
             CALL CALCP132p1            
           ENDIF
         ENDIF

       ENDIF



      ELSEIF (1.0.LE.SO4RAT .AND. SO4RAT.LT.2.0) THEN

       IF(METSTBL.EQ.1) THEN
         SCASE = 'L9'
         CALL CALCL92p1                
       ELSE

         IF (RH.LT.DRNH4HS4) THEN
            SCASE = 'L1'
            CALL CALCL12p1            

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN
            SCASE = 'L2'
            CALL CALCL22p1            

         ELSEIF (DRNAHSO4.LE.RH .AND. RH.LT.DRLC) THEN
            SCASE = 'L3'
            CALL CALCL32p1            

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'L4'
            CALL CALCL42p1            

         ELSEIF (DRNH42S4.LE.RH .AND. RH.LT.DRKHSO4) THEN
            SCASE = 'L5'
            CALL CALCL52p1            

         ELSEIF (DRKHSO4.LE.RH .AND. RH.LT.DRMGSO4) THEN
            SCASE = 'L6'
            CALL CALCL62p1            

         ELSEIF (DRMGSO4.LE.RH .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'L7'
            CALL CALCL72p1            

         ELSEIF (DRNA2SO4.LE.RH .AND. RH.LT.DRK2SO4) THEN
            SCASE = 'L8'
            CALL CALCL82p1            

	 ELSEIF (DRK2SO4.LE.RH) THEN
            SCASE = 'L9'
            CALL CALCL92p1            
	 ENDIF
       ENDIF

      CALL CALCNHA2p1                
      CALL CALCNH32p1                



      ELSEIF (SO4RAT.LT.1.0) THEN

       IF(METSTBL.EQ.1) THEN
         SCASE = 'K4'
         CALL CALCK42p1                 
       ELSE

         IF (RH.LT.DRNH4HS4) THEN                   
            SCASE = 'K1'
            CALL CALCK12p1           

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN
            SCASE = 'K2'
            CALL CALCK22p1           

         ELSEIF (DRNAHSO4.LE.RH .AND. RH.LT.DRKHSO4) THEN
            SCASE = 'K3'
            CALL CALCK32p1           

         ELSEIF (DRKHSO4.LE.RH) THEN
            SCASE = 'K4'
            CALL CALCK42p1           
         ENDIF
       ENDIF

      CALL CALCNHA2p1                  
      CALL CALCNH32p1                  

      ENDIF

      RETURN
      END
























      SUBROUTINE CALCA22p1
      INCLUDE 'module_isrpia_inc.F'



      CALAOU    =.TRUE.       
      OMELO     = TINY        
      OMEHI     = 2.0D0*W(2)  



      MOLAL(5) = W(2)
      MOLAL(6) = ZERO
      CALL CALCMR2p1



      X1 = OMEHI
      Y1 = FUNCA22p1 (X1)
      IF (ABS(Y1).LE.EPS) RETURN



      DX = (OMEHI-OMELO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, OMELO)
         Y2 = FUNCA22p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE
      IF (ABS(Y2).LE.EPS) THEN
         RETURN
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCA2')    
         RETURN
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCA22p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCA2')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCA22p1 (X3)
      RETURN



      END













      DOUBLE PRECISION FUNCTION FUNCA22p1 (OMEGI)
      INCLUDE 'module_isrpia_inc.F'
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
         MOLAL (5) = MAX(PSI-LAMDA,TINY)                          
         MOLAL (3) = MAX(W(3)/(ONE/A2/OMEGI + ONE), 2.*MOLAL(5))  
         MOLAL (6) = LAMDA                                        
         GNH3      = MAX (W(3)-MOLAL(3), TINY)                    
         COH       = ZETA                                         



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT2p1     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    DENOM = (2.0*MOLAL(5)+MOLAL(6))
      FUNCA22p1= (MOLAL(3)/DENOM - ONE) + MOLAL(1)/DENOM
      RETURN



      END






















      SUBROUTINE CALCA12p1
      INCLUDE 'module_isrpia_inc.F'

      CNH42S4 = W(2)
      GNH3    = MAX (W(3)-2.0*CNH42S4, ZERO)
      RETURN



      END
























      SUBROUTINE CALCB42p1
      INCLUDE 'module_isrpia_inc.F'



      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.



      CALL CALCB1A2p1         
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
         CALL CALCMR2p1                                           



         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT2p1
20    CONTINUE

30    RETURN



      END


















      SUBROUTINE CALCB32p1
      INCLUDE 'module_isrpia_inc.F'



      X = MAX(2*W(2)-W(3), ZERO)   
      Y = MAX(W(3)  -W(2), ZERO)   



      IF (X.LT.Y) THEN             
         SCASE   = 'B3 ; SUBCASE 1'
         TLC     = X
         TNH42S4 = Y-X
         CALL CALCB3A2p1 (TLC,TNH42S4)      
      ELSE
         SCASE   = 'B3 ; SUBCASE 2'
         TLC     = Y
         TNH4HS4 = X-Y
         CALL CALCB3B2p1 (TLC,TNH4HS4)      
      ENDIF

      RETURN



      END



























      SUBROUTINE CALCB3A2p1 (TLC, TNH42S4)
      INCLUDE 'module_isrpia_inc.F'

      CALAOU = .TRUE.         
      ZLO    = ZERO           
      ZHI    = TNH42S4        



      Z1 = ZLO
      Y1 = FUNCB3A2p1 (Z1, TLC, TNH42S4)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1



      DZ = (ZHI-ZLO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         Z2 = Z1+DZ
         Y2 = FUNCB3A2p1 (Z2, TLC, TNH42S4)
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
         CALL PUSHERR2p1 (0001, 'CALCB3A')    
         RETURN
      ENDIF



20    DO 30 I=1,MAXIT
         Z3 = 0.5*(Z1+Z2)
         Y3 = FUNCB3A2p1 (Z3, TLC, TNH42S4)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            Z2    = Z3
         ELSE
            Y1    = Y3
            Z1    = Z3
         ENDIF
         IF (ABS(Z2-Z1) .LE. EPS*Z1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCB3A')    



40    ZK = 0.5*(Z1+Z2)
      Y3 = FUNCB3A2p1 (ZK, TLC, TNH42S4)

      RETURN



      END













      DOUBLE PRECISION FUNCTION FUNCB3A2p1 (ZK, Y, X)
      INCLUDE 'module_isrpia_inc.F'
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
         CALL CALCMR2p1                   



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT2p1     
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE




30    FUNCB3A2p1= MOLAL(5)*MOLAL(3)**2.0
      FUNCB3A2p1= FUNCB3A2p1/(XK7*(WATER/GAMA(4))**3.0) - ONE
      RETURN



      END






















      SUBROUTINE CALCB3B2p1 (Y, X)
      INCLUDE 'module_isrpia_inc.F'
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
         CALL CALCMR2p1                      



         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT2p1     
20    CONTINUE

30    RETURN



      END






















      SUBROUTINE CALCB22p1
      INCLUDE 'module_isrpia_inc.F'



      X = MAX(2*W(2)-W(3), TINY)   
      Y = MAX(W(3)  -W(2), TINY)   



      IF (X.LE.Y) THEN             
         SCASE = 'B2 ; SUBCASE 1'
         CALL CALCB2A2p1 (X,Y-X)      
      ELSE
         SCASE = 'B2 ; SUBCASE 2'
         CALL CALCB2B2p1 (Y,X-Y)      
      ENDIF

      RETURN



      END





























      SUBROUTINE CALCB2A2p1 (TLC, TNH42S4)
      INCLUDE 'module_isrpia_inc.F'



      IF (RH.LT.DRMLCAS) THEN    
         SCASE   = 'B2 ; SUBCASE A1'    
         CLC     = TLC
         CNH42S4 = TNH42S4
         SCASE   = 'B2 ; SUBCASE A1'
      ELSE
         SCASE = 'B2 ; SUBCASE A2'
         CALL CALCB2A22p1 (TLC, TNH42S4)   
         SCASE = 'B2 ; SUBCASE A2'
      ENDIF

      RETURN



      END


























      SUBROUTINE CALCB2A22p1 (TLC, TNH42S4)
      INCLUDE 'module_isrpia_inc.F'



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
      CALL CALCB32p1                        



      MOLAL(1)= ONEMWF*MOLAL(1)                                   
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + 3.D0*(CLCO-CLC)) 
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               
      MOLAL(6)= ONEMWF*(CLCO-CLC)                                 

      WATER   = ONEMWF*WATER

      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4

      RETURN



      END




























      SUBROUTINE CALCB2B2p1 (TLC,TNH4HS4)
      INCLUDE 'module_isrpia_inc.F'

      CALAOU = .TRUE.       
      ZLO    = ZERO
      ZHI    = TLC          



      X1 = ZHI
      Y1 = FUNCB2B2p1 (X1,TNH4HS4,TLC)
      IF (ABS(Y1).LE.EPS) RETURN
      YHI= Y1                        



      DX = (ZHI-ZLO)/NDIV
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCB2B2p1 (X2,TNH4HS4,TLC)
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
         CALL PUSHERR2p1 (0001, 'CALCB2B')    
         RETURN
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCB2B2p1 (X3,TNH4HS4,TLC)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCB2B')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCB2B2p1 (X3,TNH4HS4,TLC)

      RETURN



      END













      DOUBLE PRECISION FUNCTION FUNCB2B2p1 (X,TNH4HS4,TLC)
      INCLUDE 'module_isrpia_inc.F'



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
         CALL CALCMR2p1                               



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT2p1     
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE




30    FUNCB2B2p1= (MOLAL(3)**3.)*MOLAL(5)*MOLAL(6)
      FUNCB2B2p1= FUNCB2B2p1/(XK13*(WATER/GAMA(13))**5.) - ONE
      RETURN



      END
























      SUBROUTINE CALCB12p1
      INCLUDE 'module_isrpia_inc.F'



      IF (RH.LT.DRMLCAB) THEN    
         SCASE = 'B1 ; SUBCASE 1'  
         CALL CALCB1A2p1              
         SCASE = 'B1 ; SUBCASE 1'
      ELSE
         SCASE = 'B1 ; SUBCASE 2'
         CALL CALCB1B2p1              
         SCASE = 'B1 ; SUBCASE 2'
      ENDIF

      RETURN



      END



























      SUBROUTINE CALCB1A2p1
      INCLUDE 'module_isrpia_inc.F'



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



      END




























      SUBROUTINE CALCB1B2p1
      INCLUDE 'module_isrpia_inc.F'



      IF (WFTYP.EQ.0) THEN
         WF = ZERO
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (DRNH4HS4-RH)/(DRNH4HS4-DRMLCAB)
      ENDIF
      ONEMWF  = ONE - WF



      CALL CALCB1A2p1
      CLCO     = CLC               
      CNH42SO  = CNH42S4
      CNH4HSO  = CNH4HS4



      CLC     = ZERO
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CALL CALCB22p1                  



      MOLAL(1)= ONEMWF*MOLAL(1)                                   
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + (CNH4HSO-CNH4HS4) & 
                     + 3.D0*(CLCO-CLC))                          
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               
      MOLAL(6)= ONEMWF*(CNH4HSO-CNH4HS4 + CLCO-CLC)               

      WATER   = ONEMWF*WATER

      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4

      RETURN



      END



















      SUBROUTINE CALCC22p1
      INCLUDE 'module_isrpia_inc.F'
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
         CALL CALCMR2p1                                       



         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT2p1     
20    CONTINUE

30    RETURN



      END





















      SUBROUTINE CALCC12p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION KLO, KHI

      CALAOU = .TRUE.    
      KLO    = TINY    
      KHI    = W(3)



      X1 = KLO
      Y1 = FUNCC12p1 (X1)
      IF (ABS(Y1).LE.EPS) GOTO 50
      YLO= Y1



      DX = (KHI-KLO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCC12p1 (X2)
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
         CALL PUSHERR2p1 (0001, 'CALCC1')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCC12p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCC1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCC12p1 (X3)

50    RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCC12p1 (KAPA)
      INCLUDE 'module_isrpia_inc.F'
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
         CALL CALCMR2p1                             



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT2p1     
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE




30    FUNCC12p1= (MOLAL(3)*MOLAL(6)/PAR2) - ONE
      RETURN



      END


















      SUBROUTINE CALCD32p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,    &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,    &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,  &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,      &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,&
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCD1A2p1



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
      CALL CALCMR2p1                  

      CALAOU = .TRUE.              
      PSI4LO = TINY                
      PSI4HI = CHI4                



60    X1 = PSI4LO
      Y1 = FUNCD32p1 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1                 



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCD32p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YHI= Y1                      
      IF (ABS(Y2) .LT. EPS) THEN   
         RETURN






      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = TINY 
         YY = FUNCD32p1(P4)
         GOTO 50






      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) 
         IF (PSI4LO.LT.-(PSI1+PSI2)) THEN
            CALL PUSHERR2p1 (0001, 'CALCD3')  
            RETURN
         ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR2p1                  
            GOTO 60                        
         ENDIF
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCD32p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCD3')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCD32p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF
      RETURN



      END













      DOUBLE PRECISION FUNCTION FUNCD32p1 (P4)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,    &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,    &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,  &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,      &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,&
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
         CALL CALCMR2p1                                 



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT2p1     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    CONTINUE

      FUNCD32p1= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 
      RETURN



      END


















      SUBROUTINE CALCD22p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCD1A2p1



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
      CALL CALCMR2p1                  

      CALAOU = .TRUE.              
      PSI4LO = TINY                
      PSI4HI = CHI4                



60    X1 = PSI4LO
      Y1 = FUNCD22p1 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1                 



      DX   = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCD22p1 (X2)
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
         YY = FUNCD22p1(P4)
         GOTO 50






      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) 
         IF (PSI4LO.LT.-(PSI1+PSI2)) THEN
            CALL PUSHERR2p1 (0001, 'CALCD2')  
            RETURN
         ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR2p1                  
            GOTO 60                        
         ENDIF
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCD22p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCD2')    



40    X3 = MIN(X1,X2)   
      Y3 = FUNCD22p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF
      RETURN



      END













      DOUBLE PRECISION FUNCTION FUNCD22p1 (P4)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL RSTGAM2p1       
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
            CALL POLY32p1 (PSI14,0.25*PSI14**2.,-A2/4.D0, PSI2, ISLV)  
            IF (ISLV.EQ.0) THEN
                PSI2 = MIN (PSI2, CHI2)
            ELSE
                PSI2 = TINY
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
         CALL CALCMR2p1                                  



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT2p1     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    CONTINUE

      FUNCD22p1= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 
      RETURN



      END






















      SUBROUTINE CALCD12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCD1A2p1, CALCD22p1



      IF (RH.LT.DRMASAN) THEN    
         SCASE = 'D1 ; SUBCASE 1'   
         CALL CALCD1A2p1            
         SCASE = 'D1 ; SUBCASE 1'
      ELSE
         SCASE = 'D1 ; SUBCASE 2'   
         CALL CALCMDRH2p1 (RH, DRMASAN, DRNH4NO3, CALCD1A2p1, CALCD22p1)
         SCASE = 'D1 ; SUBCASE 2'
      ENDIF

      RETURN



      END

























      SUBROUTINE CALCD1A2p1
      INCLUDE 'module_isrpia_inc.F'



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



      END


















      SUBROUTINE CALCG52p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEG2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,  &
                     PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,   &
                     A1,   A2,   A3,   A4,   A5,   A6,   A7



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
      Y1 = FUNCG5A2p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50  





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCG5A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG5A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCG5A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCG5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCG5A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    
         MOLAL(5) = MOLAL(5) - DELTA                    
         MOLAL(6) = DELTA                               
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCG5A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEG2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,  &
                     PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,   &
                     A1,   A2,   A3,   A4,   A5,   A6,   A7



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
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)              
      GHNO3     = MAX(CHI5 - PSI5, TINY)              
      GHCL      = MAX(CHI6 - PSI6, TINY)              

      CNH42S4   = ZERO                                
      CNH4NO3   = ZERO                                
      CNH4CL    = ZERO                                

      CALL CALCMR2p1                                     



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCG5A2p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END



















      SUBROUTINE CALCG42p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEG2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,  &
                     PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,   &
                     A1,   A2,   A3,   A4,   A5,   A6,   A7



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
      Y1 = FUNCG4A2p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50  





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2  = X1+DX
         Y2  = FUNCG4A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1  = X2
         Y1  = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG4A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCG4A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCG4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCG4A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCG4A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA, NAI, NH4I, NO3I
      COMMON /CASEG2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,  &
                     PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,   &
                     A1,   A2,   A3,   A4,   A5,   A6,   A7



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

      CALL CALCPH2p1(2.d0*SO4I+NO3I+CLI-NAI-NH4I, HI, OHI)



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCG4A2p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END



















      SUBROUTINE CALCG32p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCG1A2p1, CALCG42p1



      IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN 
         SCASE = 'G3 ; SUBCASE 1'  
         CALL CALCG3A2p1
         SCASE = 'G3 ; SUBCASE 1' 
      ELSE                                      
         SCASE = 'G1 ; SUBCASE 1'  
         CALL CALCG1A2p1
         SCASE = 'G1 ; SUBCASE 1'  
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMG3) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCG1A2p1
            SCASE = 'G3 ; SUBCASE 2'  
            RETURN
         ELSE
            SCASE = 'G3 ; SUBCASE 3'  
            CALL CALCMDRH2p1 (RH, DRMG3, DRNH42S4, CALCG1A2p1, CALCG42p1)
            SCASE = 'G3 ; SUBCASE 3'  
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCG3A2p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEG2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
                     PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,  &
                     A1,   A2,   A3,   A4,   A5,   A6,   A7



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
      Y1 = FUNCG3A2p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50  





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2  = X1+DX 
         Y2  = FUNCG3A2p1 (X2)

         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1  = X2
         Y1  = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG3A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCG3A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCG3A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCG3A2p1 (X3)



50    CONTINUE



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
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
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCG3A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEG2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
                     PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,  &
                     A1,   A2,   A3,   A4,   A5,   A6,   A7



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
         CALL POLY32p1 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (PSI20, CHI2)
      ENDIF



      MOLAL (2) = ZERO                                
      MOLAL (3) = 2.0*PSI2 + PSI4                     
      MOLAL (4) = PSI6                                
      MOLAL (5) = PSI2                                
      MOLAL (6) = ZERO                                
      MOLAL (7) = PSI5                                

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)              
      GHNO3     = MAX(CHI5 - PSI5, TINY)              
      GHCL      = MAX(CHI6 - PSI6, TINY)              

      CNH42S4   = CHI2 - PSI2                         
      CNH4NO3   = ZERO                                
      CNH4CL    = ZERO                                

      CALL CALCMR2p1                                     



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCG3A2p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END



















      SUBROUTINE CALCG22p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCG1A2p1, CALCG3A2p1, CALCG42p1



      IF (W(4).GT.TINY) THEN        
         SCASE = 'G2 ; SUBCASE 1'  
         CALL CALCG2A2p1
         SCASE = 'G2 ; SUBCASE 1' 
      ELSE                          
         SCASE = 'G1 ; SUBCASE 1'  
         CALL CALCG1A2p1
         SCASE = 'G1 ; SUBCASE 1'  
      ENDIF



      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMG2) THEN             
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCG1A2p1
            SCASE = 'G2 ; SUBCASE 2'  
         ELSE
            IF (W(5).GT. TINY) THEN
               SCASE = 'G2 ; SUBCASE 3'    
               CALL CALCMDRH2p1 (RH, DRMG2, DRNH4CL, CALCG1A2p1, CALCG3A2p1)
               SCASE = 'G2 ; SUBCASE 3'  
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMG3) THEN
               SCASE = 'G2 ; SUBCASE 4'    
               CALL CALCMDRH2p1 (RH, DRMG3, DRNH42S4, CALCG1A2p1, CALCG42p1)
               SCASE = 'G2 ; SUBCASE 4'  
            ELSE
               WATER = TINY
               DO 20 I=1,NIONS
                  MOLAL(I) = ZERO
20             CONTINUE
               CALL CALCG1A2p1
               SCASE = 'G2 ; SUBCASE 2'  
            ENDIF
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCG2A2p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEG2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,  &
                     PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,   &
                     A1,   A2,   A3,   A4,   A5,   A6,   A7



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
      Y1 = FUNCG2A2p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50  





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCG2A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) WATER = TINY
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCG2A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCG2A')    



40    X3 = 0.5*(X1+X2)
      IF (X3.LE.TINY2) THEN   
         WATER = TINY
      ELSE
         Y3 = FUNCG2A2p1 (X3)
      ENDIF



50    CONTINUE



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
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
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA     
         MOLAL(5) = MOLAL(5) - DELTA     
         MOLAL(6) = DELTA                
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCG2A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEG2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,   &
                     PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,    &
                     A1,   A2,   A3,   A4,   A5,   A6,   A7



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
         CALL POLY32p1 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (PSI20, CHI2)
      ENDIF



      MOLAL (2) = ZERO                             
      MOLAL (3) = 2.0*PSI2 + PSI4                  
      MOLAL (4) = PSI6                             
      MOLAL (5) = PSI2                             
      MOLAL (6) = ZERO                             
      MOLAL (7) = PSI5                             


      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH2p1 (SMIN, HI, OHI)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    IF (CHI4.LE.TINY) THEN
         FUNCG2A2p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
      ELSE
         FUNCG2A2p1 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
      ENDIF

      RETURN



      END























      SUBROUTINE CALCG12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCG1A2p1, CALCG2A2p1



      IF (RH.LT.DRMG1) THEN    
         SCASE = 'G1 ; SUBCASE 1'  
         CALL CALCG1A2p1              
         SCASE = 'G1 ; SUBCASE 1'
      ELSE
         SCASE = 'G1 ; SUBCASE 2'  
         CALL CALCMDRH2p1 (RH, DRMG1, DRNH4NO3, CALCG1A2p1, CALCG2A2p1)
         SCASE = 'G1 ; SUBCASE 2'
      ENDIF

      RETURN



      END

























      SUBROUTINE CALCG1A2p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2



      CNA2SO4 = MIN (0.5*W(1), W(2))
      FRNA    = MAX(W(1) - 2.D0*CNA2SO4, ZERO)
      SO4FR   = MAX(W(2) - CNA2SO4, ZERO)

      CNH42S4 = MAX (SO4FR , ZERO)                  



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
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND. &
             BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND. &
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



      END


















      SUBROUTINE CALCH62p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
      Y1 = FUNCH6A2p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH6A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH6A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH6A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCH6')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH6A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCH6A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH6A2p1 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      RETURN



      END



















      SUBROUTINE CALCH52p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      IF (W(4).LE.TINY .AND. W(5).LE.TINY) THEN  
         SCASE = 'H5'  
         CALL CALCH1A2p1
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
      Y1 = FUNCH5A2p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH5A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH5A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH5A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCH5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH5A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCH5A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
         CALL POLY32p1 (AA, BB, CC, PSI1, ISLV)
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
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 

      CALL CALCMR2p1                               



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH5A2p1 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      RETURN



      END



















      SUBROUTINE CALCH42p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      IF (W(4).LE.TINY .AND. W(5).LE.TINY) THEN  
         SCASE = 'H4'  
         CALL CALCH1A2p1
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
      Y1 = FUNCH4A2p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH4A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH4A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH4A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCH4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH4A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                      
         MOLAL(5) = MOLAL(5) - DELTA                      
         MOLAL(6) = DELTA                                 
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCH4A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
         CALL POLY32p1 (AA, BB, CC, PSI1, ISLV)
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
      CALL CALCPH2p1 (SMIN, HI, OHI)
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

      CALL CALCMR2p1                           



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH4A2p1 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      RETURN



      END



















      SUBROUTINE CALCH32p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      IF (W(4).LE.TINY) THEN        
         SCASE = 'H3'  
         CALL CALCH1A2p1
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
      Y1 = FUNCH3A2p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH3A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH3A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH3A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCH3')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH3A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCH3A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
         CALL POLY32p1 (AA, BB, CC, PSI1, ISLV)
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
      CALL CALCPH2p1 (SMIN, HI, OHI)
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

      CALL CALCMR2p1                                 



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH3A2p1 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      RETURN



      END



























      SUBROUTINE CALCH22p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCH1A2p1, CALCH32p1



      IF (W(4).GT.TINY) THEN        
         SCASE = 'H2 ; SUBCASE 1'  
         CALL CALCH2A2p1                                   
         SCASE = 'H2 ; SUBCASE 1'  
      ELSE                          
         SCASE = 'H2 ; SUBCASE 1'  
         CALL CALCH1A2p1
         SCASE = 'H2 ; SUBCASE 1'  
      ENDIF

      IF (WATER.LE.TINY .AND. RH.LT.DRMH2) THEN      
         SCASE = 'H2 ; SUBCASE 2'  

      ELSEIF (WATER.LE.TINY .AND. RH.GE.DRMH2) THEN  
         SCASE = 'H2 ; SUBCASE 3'
         CALL CALCMDRH2p1 (RH, DRMH2, DRNANO3, CALCH1A2p1, CALCH32p1)
         SCASE = 'H2 ; SUBCASE 3'
      ENDIF

      RETURN



      END






















      SUBROUTINE CALCH2A2p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
      Y1 = FUNCH2A2p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH2A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2 .GT. EPS) Y2 = FUNCH2A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCH2A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCH2A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCH2A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    
         MOLAL(5) = MOLAL(5) - DELTA                    
         MOLAL(6) = DELTA                               
      ENDIF

      RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCH2A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
         CALL POLY32p1 (AA, BB, CC, PSI1, ISLV)
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
      CALL CALCPH2p1 (SMIN, HI, OHI)
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

      CALL CALCMR2p1                        



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCH2A2p1 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A64 - ONE

      RETURN



      END
























      SUBROUTINE CALCH12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCH1A2p1, CALCH2A2p1



      IF (RH.LT.DRMH1) THEN    
         SCASE = 'H1 ; SUBCASE 1'  
         CALL CALCH1A2p1              
         SCASE = 'H1 ; SUBCASE 1'
      ELSE
         SCASE = 'H1 ; SUBCASE 2'  
         CALL CALCMDRH2p1 (RH, DRMH1, DRNH4NO3, CALCH1A2p1, CALCH2A2p1)
         SCASE = 'H1 ; SUBCASE 2'
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCH1A2p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR, &
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
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND. &
             BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND.  &
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



      END


















      SUBROUTINE CALCI62p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCI1A2p1



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
      CALL CALCMR2p1                                         



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE

20    RETURN



      END



















      SUBROUTINE CALCI52p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCI1A2p1



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
         Y1 = FUNCI5A2p1 (ZERO)
         GOTO 50
      ENDIF



      X1 = PSI4HI
      Y1 = FUNCI5A2p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCI5A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCI5A2p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCI5')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI5A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCI5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI5A2p1 (X3)

50    RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCI5A2p1 (P4)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
      CALL CALCMR2p1                                 



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0    
      FUNCI5A2p1= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END


















      SUBROUTINE CALCI42p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCI1A2p1



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
         Y1 = FUNCI4A2p1 (ZERO)
         GOTO 50
      ENDIF



      X1 = PSI4HI
      Y1 = FUNCI4A2p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCI4A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCI4A2p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCI4')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI4A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCI4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI4A2p1 (X3)

50    RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCI4A2p1 (P4)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
      PSI5 = MAX (MIN (PSI5, CHI5), ZERO)



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
      CALL CALCMR2p1                                 



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0    
      FUNCI4A2p1= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END


























      SUBROUTINE CALCI32p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCI1A2p1, CALCI42p1



      CALL CALCI1A2p1



      IF (CNH4HS4.GT.TINY .OR. CNAHSO4.GT.TINY) THEN
         SCASE = 'I3 ; SUBCASE 1'  
         CALL CALCI3A2p1                     
         SCASE = 'I3 ; SUBCASE 1'  
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMI3) THEN         
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCI1A2p1
            SCASE = 'I3 ; SUBCASE 2'  

         ELSEIF (RH.GE.DRMI3) THEN     
            SCASE = 'I3 ; SUBCASE 3'
            CALL CALCMDRH2p1 (RH, DRMI3, DRLC, CALCI1A2p1, CALCI42p1)
            SCASE = 'I3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END





















      SUBROUTINE CALCI3A2p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCI1A2p1         



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
      Y1 = FUNCI3A2p1 (X1)
      YHI= Y1                      



      IF (YHI.LT.EPS) GOTO 50



      DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = FUNCI3A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCI3A2p1 (ZERO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI3A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCI3A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI3A2p1 (X3)

50    RETURN



      END



















      DOUBLE PRECISION FUNCTION FUNCI3A2p1 (P2)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI2   = P2                  
      PSI4LO = ZERO                
      PSI4HI = CHI4                



      IF (CHI4.LE.TINY) THEN
         FUNCI3A2p1 = FUNCI3B2p1 (ZERO)
         GOTO 50
      ENDIF



      X1 = PSI4HI
      Y1 = FUNCI3B2p1 (X1)
      IF (ABS(Y1).LE.EPS) GOTO 50
      YHI= Y1                      



      IF (YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI4LO)
         Y2 = FUNCI3B2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCI3B2p1 (PSI4LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI3B2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0004, 'FUNCI3A2p1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI3B2p1 (X3)



50    A2      = XK13*(WATER/GAMA(13))**5.0
      FUNCI3A2p1 = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.D0/A2 - ONE
      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCI3B2p1 (P4)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
      PSI5 = MAX (MIN (PSI5, CHI5), ZERO)



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
      CALL CALCMR2p1                                       



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0    
      FUNCI3B2p1= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END


























      SUBROUTINE CALCI22p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCI1A2p1, CALCI3A2p1



      CALL CALCI1A2p1



      IF (CNH4HS4.GT.TINY) THEN
         SCASE = 'I2 ; SUBCASE 1'  
         CALL CALCI2A2p1                       
         SCASE = 'I2 ; SUBCASE 1'  
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMI2) THEN         
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCI1A2p1
            SCASE = 'I2 ; SUBCASE 2'  

         ELSEIF (RH.GE.DRMI2) THEN     
            SCASE = 'I2 ; SUBCASE 3'
            CALL CALCMDRH2p1 (RH, DRMI2, DRNAHSO4, CALCI1A2p1, CALCI3A2p1)
            SCASE = 'I2 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCI2A2p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCI1A2p1    



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
      Y1 = FUNCI2A2p1 (X1)
      YHI= Y1                      



      IF (YHI.LT.EPS) GOTO 50



      DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = FUNCI2A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCI2A2p1 (ZERO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCI2A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCI2A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCI2A2p1 (X3)

50    RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCI2A2p1 (P2)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,    &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,    &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,  &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,      &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,&
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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
         CALL POLY32p1 (AA, BB, CC, PSI4, ISLV)
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
         CALL POLY32p1 (AA, BB, CC, PSI3, ISLV)
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
      CALL CALCMR2p1                                



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A2      = XK13*(WATER/GAMA(13))**5.0
      FUNCI2A2p1 = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.D0/A2 - ONE
      RETURN



      END























      SUBROUTINE CALCI12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCI1A2p1, CALCI2A2p1



      IF (RH.LT.DRMI1) THEN    
         SCASE = 'I1 ; SUBCASE 1'  
         CALL CALCI1A2p1              
         SCASE = 'I1 ; SUBCASE 1'
      ELSE
         SCASE = 'I1 ; SUBCASE 2'  
         CALL CALCMDRH2p1 (RH, DRMI1, DRNH4HS4, CALCI1A2p1, CALCI2A2p1)
         SCASE = 'I1 ; SUBCASE 2'
      ENDIF





      RETURN



      END




















      SUBROUTINE CALCI1A2p1
      INCLUDE 'module_isrpia_inc.F'



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



      END

















      SUBROUTINE CALCJ32p1
      INCLUDE 'module_isrpia_inc.F'

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

      CALL CALCMR2p1                              



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 50
      ENDIF
10    CONTINUE

50    RETURN



      END


















      SUBROUTINE CALCJ22p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEJ2p1/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, &
                     A1,   A2,   A3



      CALAOU = .TRUE.              
      CHI1   = W(1)                
      CHI2   = W(3)                
      PSI1LO = TINY                
      PSI1HI = CHI1                



      X1 = PSI1HI
      Y1 = FUNCJ22p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCJ22p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCJ22p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCJ2')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCJ22p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCJ2')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCJ22p1 (X3)

50    RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCJ22p1 (P1)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEJ2p1/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, &
                     A1,   A2,   A3



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

      CALL CALCMR2p1                               



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCJ22p1 = MOLAL(2)*MOLAL(6)/A1 - ONE



      END



















      SUBROUTINE CALCJ12p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEJ2p1/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, &
                     A1,   A2,   A3



      CALAOU =.TRUE.               
      CHI1   = W(1)                
      CHI2   = W(3)                

      PSI1LO = TINY                
      PSI1HI = CHI1                



      X1 = PSI1HI
      Y1 = FUNCJ12p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCJ12p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCJ12p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCJ1')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCJ12p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCJ1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCJ12p1 (X3)

50    RETURN



      END






















      DOUBLE PRECISION FUNCTION FUNCJ12p1 (P1)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEJ2p1/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, &
                     A1,   A2,   A3



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

      CALL CALCMR2p1                               



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCJ12p1 = MOLAL(2)*MOLAL(6)/A1 - ONE



      END


















      SUBROUTINE CALCO72p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,      &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,     &
                     A5, A6, A7, A8, A9



      CALAOU = .TRUE.
      CHI9   = MIN (W(6), W(2))                     
      SO4FR  = MAX (W(2)-CHI9, ZERO)
      CAFR   = MAX (W(6)-CHI9, ZERO)
      CHI7   = MIN (0.5D0*W(7), SO4FR)              
      FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
      SO4FR  = MAX (SO4FR - CHI7, ZERO)
      CHI1   = MIN (0.5D0*W(1), SO4FR)              
      NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
      SO4FR  = MAX (SO4FR - CHI1, ZERO)
      CHI8   = MIN (W(8), SO4FR)                    
      FRMG    = MAX(W(8) - CHI8, ZERO)
      SO4FR   = MAX(SO4FR - CHI8, ZERO)
      CHI3   = ZERO
      CHI5   = W(4)
      CHI6   = W(5)
      CHI2   = MAX (SO4FR, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)

      PSI1   = CHI1
      PSI2   = CHI2
      PSI3   = ZERO
      PSI4   = ZERO
      PSI5   = ZERO
      PSI6   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    

      WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
      WATER  = MAX (WATER , TINY)



      X1 = PSI6LO
      Y1 = FUNCO72p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCO72p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCO72p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCO72p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCO7')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCO72p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    
         MOLAL(5) = MOLAL(5) - DELTA                    
         MOLAL(6) = DELTA                               
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCO72p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,      &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,     &
                     A5, A6, A7, A8, A9



      PSI6   = X
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0


      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
         PSI5 = MIN (PSI5,CHI5)
      ELSE
         PSI5 = TINY
      ENDIF


      IF(W(2).GT.TINY) THEN       
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
      ELSE
         PSI4 = TINY
      ENDIF



      MOLAL (2) = 2.0D0*PSI1                       
      MOLAL (3) = 2.0D0*PSI2 + PSI4                
      MOLAL (4) = PSI6                             
      MOLAL (5) = PSI1+PSI2+PSI7+PSI8              
      MOLAL (6) = ZERO                             
      MOLAL (7) = PSI5                             
      MOLAL (8) = ZERO                             
      MOLAL (9) = 2.0D0*PSI7                       
      MOLAL (10)= PSI8                             




       SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                   -MOLAL(9)-2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNA2SO4  = ZERO
      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4Cl   = ZERO
      CK2SO4   = ZERO
      CMGSO4   = ZERO
      CCASO4   = CHI9



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCO72p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END


















      SUBROUTINE CALCO62p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,  &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,       &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,      &
                     A5, A6, A7, A8, A9



      CALAOU = .TRUE.
      CHI9   = MIN (W(6), W(2))                     
      SO4FR  = MAX (W(2)-CHI9, ZERO)
      CAFR   = MAX (W(6)-CHI9, ZERO)
      CHI7   = MIN (0.5D0*W(7), SO4FR)              
      FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
      SO4FR  = MAX (SO4FR - CHI7, ZERO)
      CHI1   = MIN (0.5D0*W(1), SO4FR)              
      NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
      SO4FR  = MAX (SO4FR - CHI1, ZERO)
      CHI8   = MIN (W(8), SO4FR)                    
      FRMG    = MAX(W(8) - CHI8, ZERO)
      SO4FR   = MAX(SO4FR - CHI8, ZERO)
      CHI3   = ZERO
      CHI5   = W(4)
      CHI6   = W(5)
      CHI2   = MAX (SO4FR, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)


      PSI1   = CHI1
      PSI2   = CHI2
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = CHI8
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    

      WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
      WATER  = MAX (WATER , TINY)



      X1 = PSI6LO
      Y1 = FUNCO62p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCO62p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCO62p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCO62p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCO6')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCO62p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    
         MOLAL(5) = MOLAL(5) - DELTA                    
         MOLAL(6) = DELTA                               
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCO62p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,  &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,       &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,      &
                     A5, A6, A7, A8, A9



      PSI6   = X
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK17 *(WATER/GAMA(17))**3.0


      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
         PSI5 = MIN (PSI5,CHI5)
      ELSE
         PSI5 = TINY
      ENDIF


      IF(W(2).GT.TINY) THEN       
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI1+PSI2+PSI8, ZERO, -A7/4.D0, PSI7, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
         ELSE
             PSI7 = ZERO
         ENDIF
      ELSE
         PSI7 = ZERO
      ENDIF




      MOLAL (2) = 2.0D0*PSI1                       
      MOLAL (3) = 2.0D0*PSI2 + PSI4                
      MOLAL (4) = PSI6                             
      MOLAL (5) = PSI1+PSI2+PSI7+PSI8              
      MOLAL (6) = ZERO                             
      MOLAL (7) = PSI5                             
      MOLAL (8) = ZERO                             
      MOLAL (9) = 2.0D0*PSI7                       
      MOLAL (10)= PSI8                             






       SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                   -MOLAL(9)-2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3     = MAX(CHI4 - PSI4, TINY)
      GHNO3    = MAX(CHI5 - PSI5, TINY)
      GHCL     = MAX(CHI6 - PSI6, TINY)

      CNA2SO4  = ZERO
      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4Cl   = ZERO
      CK2SO4   = MAX(CHI7 - PSI7, TINY)
      CMGSO4   = ZERO
      CCASO4   = CHI9



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCO62p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END



















      SUBROUTINE CALCO52p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,      &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,     &
                     A5, A6, A7, A8, A9



      CALAOU = .TRUE.
      CHI9   = MIN (W(6), W(2))                     
      SO4FR  = MAX (W(2)-CHI9, ZERO)
      CAFR   = MAX (W(6)-CHI9, ZERO)
      CHI7   = MIN (0.5D0*W(7), SO4FR)              
      FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
      SO4FR  = MAX (SO4FR - CHI7, ZERO)
      CHI1   = MIN (0.5D0*W(1), SO4FR)              
      NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
      SO4FR  = MAX (SO4FR - CHI1, ZERO)
      CHI8   = MIN (W(8), SO4FR)                    
      FRMG    = MAX(W(8) - CHI8, ZERO)
      SO4FR   = MAX(SO4FR - CHI8, ZERO)
      CHI3   = ZERO
      CHI5   = W(4)
      CHI6   = W(5)
      CHI2   = MAX (SO4FR, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)

      PSI1   = ZERO
      PSI2   = CHI2
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = CHI8
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    

      WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
      WATER  = MAX (WATER , TINY)



      X1 = PSI6LO
      Y1 = FUNCO52p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCO52p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCO52p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCO52p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCO5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCO52p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    
         MOLAL(5) = MOLAL(5) - DELTA                    
         MOLAL(6) = DELTA                               
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCO52p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,      &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,     &
                     A5, A6, A7, A8, A9



      PSI6   = X
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK17 *(WATER/GAMA(17))**3.0


      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
         PSI5 = MIN (PSI5,CHI5)
      ELSE
         PSI5 = TINY
      ENDIF


      IF(W(2).GT.TINY) THEN       
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 ((PSI2+PSI8)/(SQRT(A1/A7)+1.D0), ZERO, &
                      -A7/4.D0/(SQRT(A1/A7)+1.D0), PSI7, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
         ELSE
             PSI7 = ZERO
         ENDIF
      ELSE
         PSI7 = ZERO
      ENDIF

      IF (CHI1.GE.TINY) THEN                              
         PSI1   = SQRT(A1/A7)*PSI7
         PSI1   = MIN(PSI1,CHI1)
      ELSE
         PSI1 = ZERO
      ENDIF




      MOLAL (2) = 2.0D0*PSI1                       
      MOLAL (3) = 2.0D0*PSI2 + PSI4                
      MOLAL (4) = PSI6                             
      MOLAL (5) = PSI1+PSI2+PSI7+PSI8              
      MOLAL (6) = ZERO                             
      MOLAL (7) = PSI5                             
      MOLAL (8) = ZERO                             
      MOLAL (9) = 2.0D0*PSI7                       
      MOLAL (10)= PSI8                             






       SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                   -MOLAL(9)-2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3     = MAX(CHI4 - PSI4, TINY)
      GHNO3    = MAX(CHI5 - PSI5, TINY)
      GHCL     = MAX(CHI6 - PSI6, TINY)

      CNA2SO4  = MAX(CHI1 - PSI1, TINY)
      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4Cl   = ZERO
      CK2SO4   = MAX(CHI7 - PSI7, TINY)
      CMGSO4   = ZERO
      CCASO4   = CHI9



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCO52p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END



















      SUBROUTINE CALCO42p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,  &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,       &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,      &
                     A5, A6, A7, A8, A9



      CALAOU = .TRUE.
      CHI9   = MIN (W(6), W(2))                     
      SO4FR  = MAX (W(2)-CHI9, ZERO)
      CAFR   = MAX (W(6)-CHI9, ZERO)
      CHI7   = MIN (0.5D0*W(7), SO4FR)              
      FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
      SO4FR  = MAX (SO4FR - CHI7, ZERO)
      CHI1   = MIN (0.5D0*W(1), SO4FR)              
      NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
      SO4FR  = MAX (SO4FR - CHI1, ZERO)
      CHI8   = MIN (W(8), SO4FR)                    
      FRMG    = MAX(W(8) - CHI8, ZERO)
      SO4FR   = MAX(SO4FR - CHI8, ZERO)
      CHI3   = ZERO
      CHI5   = W(4)
      CHI6   = W(5)
      CHI2   = MAX (SO4FR, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)

      PSI2   = CHI2
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = CHI8
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    

      WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
      WATER  = MAX (WATER , TINY)



      X1 = PSI6LO
      Y1 = FUNCO42p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCO42p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCO42p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCO42p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCO42p1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCO42p1 (X3)



50    CONTINUE



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2+PSI7+PSI8, ZERO, -A1/4.D0, PSI1, ISLV)
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
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA     
         MOLAL(5) = MOLAL(5) - DELTA     
         MOLAL(6) = DELTA                
      ENDIF

      RETURN



      END



















      DOUBLE PRECISION FUNCTION FUNCO42p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,  &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,       &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,      &
                     A5, A6, A7, A8, A9



      PSI6   = X
      PSI2   = CHI2
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK17 *(WATER/GAMA(17))**3.0




      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
         PSI5 = MIN (PSI5,CHI5)
      ELSE
         PSI5 = TINY
      ENDIF


      IF(W(2).GT.TINY) THEN       
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2+PSI8, ZERO, -A7/4.D0, PSI7, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
         ELSE
             PSI7 = ZERO
         ENDIF
      ELSE
         PSI7 = ZERO
      ENDIF



      MOLAL (2) = ZERO                             
      MOLAL (3) = 2.0D0*PSI2 + PSI4                
      MOLAL (4) = PSI6                             
      MOLAL (5) = PSI2+PSI7+PSI8                   
      MOLAL (6) = ZERO                             
      MOLAL (7) = PSI5                             
      MOLAL (8) = ZERO                             
      MOLAL (9) = 2.0D0*PSI7                       
      MOLAL (10)= PSI8                             






       SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                   -MOLAL(9)-2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3     = MAX(CHI4 - PSI4, TINY)
      GHNO3    = MAX(CHI5 - PSI5, TINY)
      GHCL     = MAX(CHI6 - PSI6, TINY)

      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4Cl   = ZERO
      CK2SO4   = MAX(CHI7 - PSI7, TINY)
      CMGSO4   = ZERO
      CCASO4   = CHI9

      CALL CALCMR2p1                                     



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCO42p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END


















      SUBROUTINE CALCO32p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCO1A2p1, CALCO42p1



      IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN 
         SCASE = 'O3 ; SUBCASE 1'
         CALL CALCO3A2p1
         SCASE = 'O3 ; SUBCASE 1'
      ELSE                                      
         SCASE = 'O1 ; SUBCASE 1'
         CALL CALCO1A2p1
         SCASE = 'O1 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMO3) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCO1A2p1
            SCASE = 'O3 ; SUBCASE 2'
            RETURN
         ELSE
            SCASE = 'O3 ; SUBCASE 3'  
            CALL CALCMDRH22p1 (RH, DRMO3, DRNH42S4, CALCO1A2p1, CALCO42p1)
            SCASE = 'O3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCO3A2p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,   &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,        &
                     PSI6, PSI7, PSI8, PSI9, A1,  A2,  A3,  A4,        &
                     A5,  A6,  A7,  A8,  A9



      CALAOU = .TRUE.
      CHI9   = MIN (W(6), W(2))                     
      SO4FR  = MAX (W(2)-CHI9, ZERO)
      CAFR   = MAX (W(6)-CHI9, ZERO)
      CHI7   = MIN (0.5D0*W(7), SO4FR)              
      FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
      SO4FR  = MAX (SO4FR - CHI7, ZERO)
      CHI1   = MIN (0.5D0*W(1), SO4FR)              
      NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
      SO4FR  = MAX (SO4FR - CHI1, ZERO)
      CHI8   = MIN (W(8), SO4FR)                    
      FRMG    = MAX(W(8) - CHI8, ZERO)
      SO4FR   = MAX(SO4FR - CHI8, ZERO)
      CHI3   = ZERO
      CHI5   = W(4)
      CHI6   = W(5)
      CHI2   = MAX (SO4FR, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)

      PSI8   = CHI8
      PSI6LO = TINY
      PSI6HI = CHI6-TINY

      WATER  = TINY



      X1 = PSI6LO
      Y1 = FUNCO3A2p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCO3A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCO3A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCO3A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCO3A2p1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCO3A2p1 (X3)



50    CONTINUE



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2+PSI7+PSI8, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (max (PSI1, zero), CHI1)
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
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA     
         MOLAL(5) = MOLAL(5) - DELTA     
         MOLAL(6) = DELTA                
      ENDIF

      RETURN



      END



















      DOUBLE PRECISION FUNCTION FUNCO3A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,  &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,       &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,      &
                     A5, A6, A7, A8, A9



      PSI2   = CHI2
      PSI8   = CHI8
      PSI3   = ZERO
      PSI6   = X

      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0D0
      A2  = XK7 *(WATER/GAMA(4))**3.0D0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0D0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0D0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0D0
      A7  = XK17 *(WATER/GAMA(17))**3.0D0

      A65 = A6/A5



      DENO = MAX(CHI6-PSI6-PSI3, ZERO)
      PSI5 = PSI6*CHI5/(A6/A5*DENO + PSI6)
      PSI5 = MIN(MAX(PSI5,ZERO),CHI5)


      IF(W(2).GT.TINY) THEN       
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)  
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF
         PSI4 = MIN (MAX (PSI4,ZERO), CHI4)

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2+PSI8, ZERO, -A7/4.D0, PSI7, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
         ELSE
             PSI7 = ZERO
         ENDIF
      ELSE
         PSI7 = ZERO
      ENDIF

      IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN
         CALL POLY32p1 (PSI7+PSI8+PSI4, PSI4*(PSI7+PSI8)+              &
                     PSI4*PSI4/4.D0, (PSI4*PSI4*(PSI7+PSI8)-A2)     &
                     /4.D0,PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (MAX(PSI20,ZERO), CHI2)
      ENDIF






      MOLAL (2) = ZERO                             
      MOLAL (3) = 2.0D0*PSI2 + PSI4                
      MOLAL (4) = PSI6                             
      MOLAL (5) = PSI2+PSI7+PSI8                   
      MOLAL (6) = ZERO                             
      MOLAL (7) = PSI5                             
      MOLAL (8) = ZERO                             
      MOLAL (9) = 2.0D0*PSI7                       
      MOLAL (10)= PSI8                             


      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                   -MOLAL(9)-2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)


      CNH42S4  = MAX(CHI2 - PSI2, ZERO)
      CNH4NO3  = ZERO
      CNH4Cl   = ZERO
      CK2SO4   = MAX(CHI7 - PSI7, ZERO)
      CMGSO4   = ZERO
      CCASO4   = CHI9



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCO3A2p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


      RETURN



      END



















      SUBROUTINE CALCO22p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCO1A2p1, CALCO3A2p1, CALCO42p1



      IF (W(4).GT.TINY) THEN        
         SCASE = 'O2 ; SUBCASE 1'
         CALL CALCO2A2p1
         SCASE = 'O2 ; SUBCASE 1'
      ELSE                          
         SCASE = 'O1 ; SUBCASE 1'
         CALL CALCO1A2p1
         SCASE = 'O1 ; SUBCASE 1'
      ENDIF



      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMO2) THEN             
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCO1A2p1
            SCASE = 'O2 ; SUBCASE 2'
         ELSE
            IF (W(5).GT. TINY) THEN
               SCASE = 'O2 ; SUBCASE 3'    
               CALL CALCMDRH22p1 (RH, DRMO2, DRNH4CL, CALCO1A2p1, CALCO3A2p1)
               SCASE = 'O2 ; SUBCASE 3'
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMO3) THEN
               SCASE = 'O2 ; SUBCASE 4'    
               CALL CALCMDRH22p1 (RH, DRMO3, DRNH42S4, CALCO1A2p1, CALCO42p1)
               SCASE = 'O2 ; SUBCASE 4'
            ELSE
               WATER = TINY
               DO 20 I=1,NIONS
                  MOLAL(I) = ZERO
20             CONTINUE
               CALL CALCO1A2p1
               SCASE = 'O2 ; SUBCASE 2'
            ENDIF
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCO2A2p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,   &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,        &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,       &
                     A5, A6, A7, A8, A9



      CALAOU = .TRUE.
      CHI9   = MIN (W(6), W(2))                     
      SO4FR  = MAX (W(2)-CHI9, ZERO)
      CAFR   = MAX (W(6)-CHI9, ZERO)
      CHI7   = MIN (0.5D0*W(7), SO4FR)              
      FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
      SO4FR  = MAX (SO4FR - CHI7, ZERO)
      CHI1   = MIN (0.5D0*W(1), SO4FR)              
      NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
      SO4FR  = MAX (SO4FR - CHI1, ZERO)
      CHI8   = MIN (W(8), SO4FR)                    
      FRMG    = MAX(W(8) - CHI8, ZERO)
      SO4FR   = MAX(SO4FR - CHI8, ZERO)
      CHI3   = ZERO
      CHI5   = W(4)
      CHI6   = W(5)
      CHI2   = MAX (SO4FR, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)

      PSI8   = CHI8
      PSI6LO = TINY
      PSI6HI = CHI6-TINY

      WATER  = TINY



      X1 = PSI6LO
      Y1 = FUNCO2A2p1 (X1)
      IF (CHI6.LE.TINY) GOTO 50





      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCO2A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) WATER = TINY
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCO2A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCO2A')    



40    X3 = 0.5*(X1+X2)
      IF (X3.LE.TINY2) THEN   
         WATER = TINY
      ELSE
         Y3 = FUNCO2A2p1 (X3)
      ENDIF



50    CONTINUE



      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2+PSI7+PSI8, ZERO, -A1/4.D0, PSI1, ISLV)
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
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA     
         MOLAL(5) = MOLAL(5) - DELTA     
         MOLAL(6) = DELTA                
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCO2A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA
      COMMON /CASEO2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,          &
                     PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,         &
                     A5, A6, A7, A8, A9



      PSI6   = X
      PSI2   = CHI2
      PSI3   = ZERO

      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0D0
      A2  = XK7 *(WATER/GAMA(4))**3.0D0
      A3  = XK6 /(R*TEMP*R*TEMP)
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0D0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0D0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0D0
      A65 = A6/A5
      A7  = XK17 *(WATER/GAMA(17))**3.0D0


      DENO = MAX(CHI6-PSI6-PSI3, ZERO)
      PSI5 = PSI6*CHI5/(A6/A5*DENO + PSI6)
      PSI5 = MIN(PSI5,CHI5)

      PSI4 = MIN(PSI5+PSI6,CHI4)


      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN        
         CALL POLY32p1 (PSI2+PSI8, ZERO, -A7/4.D0, PSI7, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
         ELSE
             PSI7 = ZERO
         ENDIF
      ELSE
         PSI7 = ZERO
      ENDIF

      IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN
         CALL POLY32p1 (PSI7+PSI8+PSI4, PSI4*(PSI7+PSI8)+              &
                     PSI4*PSI4/4.D0, (PSI4*PSI4*(PSI7+PSI8)-A2)     &
                     /4.D0,PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (MAX(PSI20,ZERO), CHI2)
      ENDIF






      MOLAL (2) = ZERO                             
      MOLAL (3) = 2.0D0*PSI2 + PSI4                
      MOLAL (4) = PSI6                             
      MOLAL (5) = PSI2+PSI7+PSI8                   
      MOLAL (6) = ZERO                             
      MOLAL (7) = PSI5                             
      MOLAL (8) = ZERO                             
      MOLAL (9) = 2.0D0*PSI7                       
      MOLAL (10)= PSI8                             


      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                   -MOLAL(9)-2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI



      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)


      CNH42S4  = MAX(CHI2 - PSI2, ZERO)
      CNH4NO3  = ZERO
      CK2SO4   = MAX(CHI7 - PSI7, ZERO)
      CMGSO4   = ZERO
      CCASO4   = CHI9
      



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
         PSI3 = MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE







20         FUNCO2A2p1 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE


      RETURN



      END






















      SUBROUTINE CALCO12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCO1A2p1, CALCO2A2p1



      IF (RH.LT.DRMO1) THEN
         SCASE = 'O1 ; SUBCASE 1'
         CALL CALCO1A2p1              
         SCASE = 'O1 ; SUBCASE 1'
      ELSE
         SCASE = 'O1 ; SUBCASE 2'  
         CALL CALCMDRH22p1 (RH, DRMO1, DRNH4NO3, CALCO1A2p1, CALCO2A2p1)
         SCASE = 'O1 ; SUBCASE 2'
      ENDIF

      RETURN



      END






















      SUBROUTINE CALCO1A2p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2



      CCASO4  = MIN (W(6), W(2))                    
      SO4FR   = MAX(W(2) - CCASO4, ZERO)
      CAFR    = MAX(W(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*W(7), SO4FR)             
      FRK     = MAX(W(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX(SO4FR - CK2SO4, ZERO)
      CNA2SO4 = MIN (0.5D0*W(1), SO4FR)             
      FRNA    = MAX(W(1) - 2.D0*CNA2SO4, ZERO)
      SO4FR   = MAX(SO4FR - CNA2SO4, ZERO)
      CMGSO4  = MIN (W(8), SO4FR)                   
      FRMG    = MAX(W(8) - CMGSO4, ZERO)
      SO4FR   = MAX(SO4FR - CMGSO4, ZERO)

      CNH42S4 = MAX (SO4FR , ZERO)                  



      ALF     = W(3) - 2.0D0*CNH42S4
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
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND. &
             BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND. &
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



      END



















      SUBROUTINE CALCM82p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    
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
      Y1 = FUNCM82p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCM82p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCM82p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCM82p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCM8')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCM82p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCM82p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,       &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,       &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,     &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,         &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,   &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = CHI9
      PSI10  = CHI10
      PSI11  = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP


      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0






      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF



      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               
      MOLAL (3) = PSI4                                  
      MOLAL (4) = PSI6 + PSI7                           
      MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            
      MOLAL (6) = ZERO                                  
      MOLAL (7) = PSI5 + PSI8                           
      MOLAL (8) = PSI11                                 
      MOLAL (9) = 2.D0*PSI9                             
      MOLAL (10)= PSI10                                 

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CNA2SO4   = ZERO
      CK2SO4    = ZERO
      CMGSO4    = ZERO
      CCASO4    = CHI11

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCM82p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



















      SUBROUTINE CALCM72p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    
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
      Y1 = FUNCM72p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCM72p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCM72p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCM72p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCM7')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCM72p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END



















      DOUBLE PRECISION FUNCTION FUNCM72p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP


      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0






      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
      CALL POLY32p1 (PSI1+PSI10,ZERO,-A9/4.D0, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MAX (MIN (PSI9,CHI9), ZERO)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF



      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               
      MOLAL (3) = PSI4                                  
      MOLAL (4) = PSI6 + PSI7                           
      MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            
      MOLAL (6) = ZERO                                  
      MOLAL (7) = PSI5 + PSI8                           
      MOLAL (8) = PSI11                                 
      MOLAL (9) = 2.D0*PSI9                             
      MOLAL (10)= PSI10                                 

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CNA2SO4   = ZERO
      CK2SO4    = MAX(CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCM72p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END

















      SUBROUTINE CALCM62p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    
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
      Y1 = FUNCM62p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCM62p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCM62p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCM62p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCM6')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCM62p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCM62p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0






      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN   
      RIZ = SQRT(A9/A1)
      AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.D0+RIZ)*(PSI7+PSI8)) &
             /(1.D0+RIZ)
      BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.D0+RIZ))/(1.D0+RIZ)
      CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
             -A1/4.D0)/(1.D0+RIZ)




      CALL POLY32p1 (AA,BB,CC,PSI1,ISLV)
        IF (ISLV.EQ.0) THEN
            PSI1 = MIN (PSI1,CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
      ENDIF








      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN   
      CALL POLY32p1 (PSI1+PSI10,ZERO,-A9/4.D0, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (PSI9,CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF



      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               
      MOLAL (3) = PSI4                                  
      MOLAL (4) = PSI6 + PSI7                           
      MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            
      MOLAL (6) = ZERO                                  
      MOLAL (7) = PSI5 + PSI8                           
      MOLAL (8) = PSI11                                 
      MOLAL (9) = 2.D0*PSI9                             
      MOLAL (10)= PSI10                                 

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
      CK2SO4    = MAX(CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCM62p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END


















      SUBROUTINE CALCM52p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    
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
      Y1 = FUNCM52p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCM52p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCM52p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCM52p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCM5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCM52p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCM52p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,    &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,    &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,  &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,      &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,&
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0






      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN   
      RIZ = SQRT(A9/A1)
      AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.D0+RIZ)*(PSI7+PSI8)) &
             /(1.D0+RIZ)
      BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.D0+RIZ))/(1.D0+RIZ)
      CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
             -A1/4.D0)/(1.D0+RIZ)




      CALL POLY32p1 (AA,BB,CC,PSI1,ISLV)
        IF (ISLV.EQ.0) THEN
            PSI1 = MIN (PSI1,CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
      ENDIF

      IF (CHI9.GE.TINY .AND. WATER.GT.TINY) THEN
         PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
         PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
      ELSE
         PSI9  = ZERO
      ENDIF












      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               
      MOLAL (3) = PSI4                                  
      MOLAL (4) = PSI6 + PSI7                           
      MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            
      MOLAL (6) = ZERO                                  
      MOLAL (7) = PSI5 + PSI8                           
      MOLAL (8) = PSI11                                 
      MOLAL (9) = 2.D0*PSI9                             
      MOLAL (10)= PSI10                                 

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
      CK2SO4    = MAX(CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCM52p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END


















      SUBROUTINE CALCM42p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      IF (W(4).LE.TINY .AND. W(5).LE.TINY) THEN
         SCASE = 'M4 ; SUBCASE 1'
         CALL CALCM1A2p1
         SCASE = 'M4 ; SUBCASE 1'
         RETURN
      ENDIF



      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    
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
      Y1 = FUNCM42p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCM42p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCM42p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCM42p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCM4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCM42p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCM42p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A3  = XK6 /(R*TEMP*R*TEMP)
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0






      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,TINY),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN   
      RIZ = SQRT(A9/A1)
      AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.D0+RIZ)*(PSI7+PSI8)) &
             /(1.D0+RIZ)
      BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.D0+RIZ))/(1.D0+RIZ)
      CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
             -A1/4.D0)/(1.D0+RIZ)




      CALL POLY32p1 (AA,BB,CC,PSI1,ISLV)
        IF (ISLV.EQ.0) THEN
            PSI1 = MIN (PSI1,CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
      ENDIF

      IF (CHI9.GE.TINY .AND. WATER.GT.TINY) THEN
         PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
         PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
      ELSE
         PSI9  = ZERO
      ENDIF












      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               
      MOLAL (3) = PSI4                                  
      MOLAL (4) = PSI6 + PSI7                           
      MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            
      MOLAL (6) = ZERO                                  
      MOLAL (7) = PSI5 + PSI8                           
      MOLAL (8) = PSI11                                 
      MOLAL (9) = 2.D0*PSI9                             
      MOLAL (10)= PSI10                                 

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
      CK2SO4    = MAX(CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11



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
      PSI3 = MAX (MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6), ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCM42p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END


















      SUBROUTINE CALCM32p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      IF (W(4).LE.TINY) THEN        
         SCASE = 'M3 ; SUBCASE 1'
         CALL CALCM1A2p1
         SCASE = 'M3 ; SUBCASE 1'
         RETURN
      ENDIF



      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    
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
      Y1 = FUNCM32p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCM32p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCM32p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCM32p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCM3')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCM32p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCM32p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A3  = XK6 /(R*TEMP*R*TEMP)
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A10 = XK23 *(WATER/GAMA(21))**2.0





      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,TINY),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF









      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
         PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF


      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN   
      RIZ = SQRT(A9/A1)
      AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.D0+RIZ)*(PSI7+PSI8)) &
             /(1.D0+RIZ)
      BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.D0+RIZ))/(1.D0+RIZ)
      CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
             -A1/4.D0)/(1.D0+RIZ)




      CALL POLY32p1 (AA,BB,CC,PSI1,ISLV)
        IF (ISLV.EQ.0) THEN
            PSI1 = MIN (PSI1,CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
      ENDIF

      IF (CHI9.GE.TINY) THEN
         PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
         PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
      ELSE
         PSI9  = ZERO
      ENDIF












      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               
      MOLAL (3) = PSI4                                  
      MOLAL (4) = PSI6 + PSI7                           
      MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            
      MOLAL (6) = ZERO                                  
      MOLAL (7) = PSI5 + PSI8                           
      MOLAL (8) = PSI11                                 
      MOLAL (9) = 2.D0*PSI9                             
      MOLAL (10)= PSI10                                 

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = ZERO
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
      CK2SO4    = MAX(CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11



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
      PSI3 = MAX (MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6), ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCM32p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



























      SUBROUTINE CALCM22p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCM1A2p1, CALCM32p1



      CALL CALCM1A2p1

      IF (CNH4NO3.GT.TINY) THEN        
         SCASE = 'M2 ; SUBCASE 1'
         CALL CALCM2A2p1
         SCASE = 'M2 ; SUBCASE 1'
      ELSE                          
         SCASE = 'M2 ; SUBCASE 1'
         CALL CALCM1A2p1
         SCASE = 'M2 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY .AND. RH.LT.DRMM2) THEN      
         SCASE = 'M2 ; SUBCASE 2'

      ELSEIF (WATER.LE.TINY .AND. RH.GE.DRMM2) THEN  
         SCASE = 'M2 ; SUBCASE 3'
         CALL CALCMDRH22p1 (RH, DRMM2, DRNANO3, CALCM1A2p1, CALCM32p1)
         SCASE = 'M2 ; SUBCASE 3'
      ENDIF

      RETURN



      END


















      SUBROUTINE CALCM2A2p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    
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
      Y1 = FUNCM2A2p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCM2A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCM2A2p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCM2A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCM2A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCM2A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCM2A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.0
      A3  = XK6 /(R*TEMP*R*TEMP)
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A64 = (XK3*XK2/XKW)*(GAMA(10)/GAMA(5)/GAMA(11))**2.0
      A64 = A64*(R*TEMP*WATER)**2.0




      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,TINY),CHI4)
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
      RIZ = SQRT(A9/A1)
      AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.D0+RIZ)*(PSI7+PSI8)) &
             /(1.D0+RIZ)
      BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.D0+RIZ))/(1.D0+RIZ)
      CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
             -A1/4.D0)/(1.D0+RIZ)





      CALL POLY32p1 (AA,BB,CC,PSI1,ISLV)
        IF (ISLV.EQ.0) THEN
            PSI1 = MIN (PSI1,CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
      ENDIF

      IF (CHI9.GE.TINY .AND. WATER.GT.TINY) THEN

         PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
         PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
      ELSE
         PSI9  = ZERO
      ENDIF












      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               
      MOLAL (3) = PSI4                                  
      MOLAL (4) = PSI6 + PSI7                           
      MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            
      MOLAL (6) = ZERO                                  
      MOLAL (7) = PSI5 + PSI8                           
      MOLAL (8) = PSI11                                 
      MOLAL (9) = 2.D0*PSI9                             
      MOLAL (10)= PSI10                                 

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
      CK2SO4    = MAX(CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11



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
      PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6), ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCM2A2p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END






















      SUBROUTINE CALCM12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCM1A2p1, CALCM2A2p1



      IF (RH.LT.DRMM1) THEN
         SCASE = 'M1 ; SUBCASE 1'
         CALL CALCM1A2p1              
         SCASE = 'M1 ; SUBCASE 1'
      ELSE
         SCASE = 'M1 ; SUBCASE 2'  
         CALL CALCMDRH22p1 (RH, DRMM1, DRNH4NO3, CALCM1A2p1, CALCM2A2p1)
         SCASE = 'M1 ; SUBCASE 2'
      ENDIF

      RETURN



      END

















      SUBROUTINE CALCM1A2p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR, &
                       NO3FR



      CCASO4  = MIN (W(6), W(2))                    
      SO4FR   = MAX(W(2) - CCASO4, ZERO)
      CAFR    = MAX(W(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*W(7), SO4FR)             
      FRK     = MAX(W(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX(SO4FR - CK2SO4, ZERO)
      CMGSO4  = MIN (W(8), SO4FR)                   
      FRMG    = MAX(W(8) - CMGSO4, ZERO)
      SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
      CNA2SO4 = MAX (SO4FR,ZERO)                    
      NAFR    = MAX (W(1)-2.D0*CNA2SO4, ZERO)
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
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND. &
             BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND. &
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



      END




















      SUBROUTINE CALCP132p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP132p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP132p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP132p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP132p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP13')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP132p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP132p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI4   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = CHI9
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = CHI13
      PSI14  = CHI14
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0



      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
             A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF



      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH4NO3   = ZERO
      CNH4CL    = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CK2SO4    = ZERO
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = ZERO
      CKCL      = ZERO
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP132p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



















      SUBROUTINE CALCP122p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP122p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP122p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP122p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP122p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP12')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP122p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP122p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI4   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = CHI13
      PSI14  = CHI14
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0



      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
             A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH4NO3   = ZERO
      CNH4CL    = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = ZERO
      CKCL      = ZERO
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP122p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END




















      SUBROUTINE CALCP112p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP112p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP112p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP112p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP112p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP11')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP112p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP112p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = CHI14
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0



      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
             A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
        DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 =0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH4NO3   = ZERO
      CNH4CL    = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = ZERO
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP112p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END




















      SUBROUTINE CALCP102p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP102p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP102p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP102p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP102p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP10')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP102p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP102p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = CHI14
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0



      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
             A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 =0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH4NO3   = ZERO
      CNH4CL    = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = ZERO
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP102p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END




















      SUBROUTINE CALCP92p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP92p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP92p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP92p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP92p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP9')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP92p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP92p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = ZERO
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0
      A14 = XK20 *(WATER/GAMA(20))**2.0



      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
             A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI14.GT.TINY .AND. WATER.GT.TINY) THEN          
         PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
                 PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
         PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH4NO3   = ZERO
      CNH4CL    = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = MAX (CHI14 - PSI14, ZERO)
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP92p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



















      SUBROUTINE CALCP82p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP82p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP82p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP82p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP82p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP8')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP82p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP82p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = ZERO
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0
      A14 = XK20 *(WATER/GAMA(20))**2.0



      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
             A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI14.GT.TINY .AND. WATER.GT.TINY) THEN          
         PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
                 PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
         PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH4NO3   = ZERO

      CNACL     = ZERO
      CNANO3    = ZERO
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = MAX (CHI14 - PSI14, ZERO)
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO



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
      PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP82p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



















      SUBROUTINE CALCP72p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP72p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP72p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP72p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP72p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP7')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP72p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP72p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = CHI8
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = ZERO
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0
      A14 = XK20 *(WATER/GAMA(20))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0



      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
             A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI14.GT.TINY .AND. WATER.GT.TINY) THEN          
         PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
                 PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
         PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
         GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
         DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
         PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH4NO3   = ZERO

      CNACL     = MAX (CHI7 - PSI7, ZERO)
      CNANO3    = ZERO
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = MAX (CHI14 - PSI14, ZERO)
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO



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
      PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP72p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



















      SUBROUTINE CALCP62p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP62p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP62p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP62p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP62p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP6')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP62p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP62p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = ZERO
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = ZERO
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0
      A14 = XK20 *(WATER/GAMA(20))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0



      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
             A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI14.GT.TINY .AND. WATER.GT.TINY) THEN          
         PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
                 PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
         PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
         GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
         DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
         PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF

      IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     




          PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
                 PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
          PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)

      CNH4NO3   = ZERO

      CNACL     = MAX (CHI7 - PSI7, ZERO)
      CNANO3    = MAX (CHI8 - PSI8, ZERO)
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = MAX (CHI14 - PSI14, ZERO)
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO



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
      PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP62p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



















      SUBROUTINE CALCP52p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCP1A2p1, CALCP62p1



      IF (W(4).GT.TINY)   THEN 
         SCASE = 'P5 ; SUBCASE 1'
         CALL CALCP5A2p1
         SCASE = 'P5 ; SUBCASE 1'
      ELSE                                      
         SCASE = 'P1 ; SUBCASE 1'
         CALL CALCP1A2p1
         SCASE = 'P1 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMP5) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCP1A2p1
            SCASE = 'P5 ; SUBCASE 2'
            RETURN
         ELSE
            SCASE = 'P5 ; SUBCASE 3'  

            CALL CALCMDRH22p1 (RH, DRMP5, DRNH4NO3, CALCP1A2p1, CALCP62p1)
            SCASE = 'P5 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCP5A2p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP52p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP52p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP52p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP52p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP52p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP52p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = ZERO
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = ZERO
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0
      A14 = XK20 *(WATER/GAMA(20))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0



      PSI5 = (CHI5-PSI2)*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) &
             - A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI14.GT.TINY .AND. WATER.GT.TINY) THEN          
         PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
                 PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
         PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
         GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
         DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
         PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF

      IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     




          PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
                 PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
          PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)



      CNACL     = MAX (CHI7 - PSI7, ZERO)
      CNANO3    = MAX (CHI8 - PSI8, ZERO)
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = MAX (CHI14 - PSI14, ZERO)
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO



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
      PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3



      A2   = XK10 /(R*TEMP*R*TEMP)
      IF (GNH3*GHNO3.GT.A2) THEN
         DELT = MIN(GNH3, GHNO3)
         BB = -(GNH3+GHNO3)
         CC = GNH3*GHNO3-A2
         DD = BB*BB - 4.D0*CC
         PSI21 = 0.5D0*(-BB + SQRT(DD))
         PSI22 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI21.GT.ZERO .AND. PSI21.GT.ZERO) THEN
            PSI2 = PSI21
         ELSEIF (DELT-PSI22.GT.ZERO .AND. PSI22.GT.ZERO) THEN
            PSI2 = PSI22
         ELSE
            PSI2 = ZERO
         ENDIF
      ELSE
         PSI2 = ZERO
      ENDIF
      PSI2 = MAX(MIN(MIN(PSI2,CHI4-PSI4-PSI3),CHI5-PSI5), ZERO)



      GNH3    = MAX(GNH3 - PSI2, TINY)
      GHCL    = MAX(GHNO3 - PSI2, TINY)
      CNH4NO3 = PSI2

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP52p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



















      SUBROUTINE CALCP42p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCP1A2p1, CALCP5A2p1



      IF (W(4).GT.TINY)   THEN 
         SCASE = 'P4 ; SUBCASE 1'
         CALL CALCP4A2p1
         SCASE = 'P4 ; SUBCASE 1'
      ELSE                                      
         SCASE = 'P1 ; SUBCASE 1'
         CALL CALCP1A2p1
         SCASE = 'P1 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMP4) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCP1A2p1
            SCASE = 'P4 ; SUBCASE 2'
            RETURN
         ELSE
            SCASE = 'P4 ; SUBCASE 3'  

            CALL CALCMDRH22p1 (RH, DRMP4, DRMGNO32, CALCP1A2p1, CALCP5A2p1)
            SCASE = 'P4 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCP4A2p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP42p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP42p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP42p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP42p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP42p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END




















      DOUBLE PRECISION FUNCTION FUNCP42p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = ZERO
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = ZERO
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0
      A14 = XK20 *(WATER/GAMA(20))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0



      PSI5 = (CHI5-PSI2)*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) &
             - A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI14.GT.TINY .AND. WATER.GT.TINY) THEN          
         PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
                 PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
         PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
         GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
         DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
         PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF

      IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     




          PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
                 PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
          PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)



      CNACL     = MAX (CHI7 - PSI7, ZERO)
      CNANO3    = MAX (CHI8 - PSI8, ZERO)
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = MAX (CHI14 - PSI14, ZERO)
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO



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
      PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3



      A2   = XK10 /(R*TEMP*R*TEMP)
      IF (GNH3*GHNO3.GT.A2) THEN
         DELT = MIN(GNH3, GHNO3)
         BB = -(GNH3+GHNO3)
         CC = GNH3*GHNO3-A2
         DD = BB*BB - 4.D0*CC
         PSI21 = 0.5D0*(-BB + SQRT(DD))
         PSI22 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI21.GT.ZERO .AND. PSI21.GT.ZERO) THEN
            PSI2 = PSI21
         ELSEIF (DELT-PSI22.GT.ZERO .AND. PSI22.GT.ZERO) THEN
            PSI2 = PSI22
         ELSE
            PSI2 = ZERO
         ENDIF
      ELSE
         PSI2 = ZERO
      ENDIF
      PSI2 = MAX(MIN(MIN(PSI2,CHI4-PSI4-PSI3),CHI5-PSI5), ZERO)



      GNH3    = MAX(GNH3 - PSI2, TINY)
      GHCL    = MAX(GHNO3 - PSI2, TINY)
      CNH4NO3 = PSI2

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP42p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



















      SUBROUTINE CALCP32p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCP1A2p1, CALCP4A2p1



      IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN 
         SCASE = 'P3 ; SUBCASE 1'
         CALL CALCP3A2p1
         SCASE = 'P3 ; SUBCASE 1'
      ELSE                                      
         SCASE = 'P1 ; SUBCASE 1'
         CALL CALCP1A2p1
         SCASE = 'P1 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMP3) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCP1A2p1
            SCASE = 'P3 ; SUBCASE 2'
            RETURN
         ELSE
            SCASE = 'P3 ; SUBCASE 3'  

            CALL CALCMDRH22p1 (RH, DRMP3, DRCANO32, CALCP1A2p1, CALCP4A2p1)
            SCASE = 'P3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCP3A2p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP32p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP32p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP32p1 (PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP32p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP3')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP32p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP32p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = ZERO
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = ZERO
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0
      A14 = XK20 *(WATER/GAMA(20))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0



      PSI5 = (CHI5-PSI2)*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) &
             - A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI14.GT.TINY .AND. WATER.GT.TINY) THEN          
         PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
                 PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
         PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
         GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
         DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
         PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF

      IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     




          PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
                 PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
          PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)



      CNACL     = MAX (CHI7 - PSI7, ZERO)
      CNANO3    = MAX (CHI8 - PSI8, ZERO)
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = MAX (CHI14 - PSI14, ZERO)
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO



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
      PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6), ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3



      A2   = XK10 /(R*TEMP*R*TEMP)
      IF (GNH3*GHNO3.GT.A2) THEN
         DELT = MIN(GNH3, GHNO3)
         BB = -(GNH3+GHNO3)
         CC = GNH3*GHNO3-A2
         DD = BB*BB - 4.D0*CC
         PSI21 = 0.5D0*(-BB + SQRT(DD))
         PSI22 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI21.GT.ZERO .AND. PSI21.GT.ZERO) THEN
            PSI2 = PSI21
         ELSEIF (DELT-PSI22.GT.ZERO .AND. PSI22.GT.ZERO) THEN
            PSI2 = PSI22
         ELSE
            PSI2 = ZERO
         ENDIF
      ELSE
         PSI2 = ZERO
      ENDIF
      PSI2 = MAX(MIN(MIN(PSI2,CHI4-PSI4-PSI3),CHI5-PSI5),ZERO)



      GNH3    = MAX(GNH3 - PSI2, TINY)
      GHCL    = MAX(GHNO3 - PSI2, TINY)
      CNH4NO3 = PSI2

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP32p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END



























      SUBROUTINE CALCP22p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCP1A2p1, CALCP3A2p1, CALCP4A2p1, CALCP5A2p1, CALCP62p1



      CALL CALCP1A2p1



      IF (CCACL2.GT.TINY) THEN
         SCASE = 'P2 ; SUBCASE 1'
         CALL CALCP2A2p1
         SCASE = 'P2 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMP2) THEN             
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCP1A2p1
            SCASE = 'P2 ; SUBCASE 2'
         ELSE
            IF (CMGCL2.GT. TINY) THEN
               SCASE = 'P2 ; SUBCASE 3'    

               CALL CALCMDRH22p1 (RH, DRMP2, DRMGCL2, CALCP1A2p1, CALCP3A2p1)
               SCASE = 'P2 ; SUBCASE 3'
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMP3 .AND. RH.LT.DRMP4) THEN
               SCASE = 'P2 ; SUBCASE 4'    

               CALL CALCMDRH22p1 (RH, DRMP3, DRCANO32, CALCP1A2p1, CALCP4A2p1)
               SCASE = 'P2 ; SUBCASE 4'
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMP4 .AND. RH.LT.DRMP5) THEN
               SCASE = 'P2 ; SUBCASE 5'    

               CALL CALCMDRH22p1 (RH, DRMP4, DRMGNO32, CALCP1A2p1, CALCP5A2p1)
               SCASE = 'P2 ; SUBCASE 5'
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMP5) THEN
               SCASE = 'P2 ; SUBCASE 6'    

               CALL CALCMDRH22p1 (RH, DRMP5, DRNH4NO3, CALCP1A2p1, CALCP62p1)
               SCASE = 'P2 ; SUBCASE 6'
            ELSE
               WATER = TINY
               DO 20 I=1,NIONS
                  MOLAL(I) = ZERO
20             CONTINUE
               CALL CALCP1A2p1
               SCASE = 'P2 ; SUBCASE 2'
            ENDIF
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCP2A2p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               
      CHI6    = FRCL                                
      CHI4    = W(3)                                

      CHI3    = ZERO                                
      CHI1    = ZERO
      CHI2    = ZERO

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    



      X1 = PSI6LO
      Y1 = FUNCP2A2p1 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP2A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (ABS(Y2) .GT. EPS) Y2 = FUNCP2A2p1(PSI6LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCP2A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCP2A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCP2A2p1 (X3)



50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS42p1 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     
         MOLAL(5) = MOLAL(5) - DELTA                     
         MOLAL(6) = DELTA                                
      ENDIF

      RETURN



      END





















      DOUBLE PRECISION FUNCTION FUNCP2A2p1 (X)
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = ZERO
      PSI8   = ZERO
      PSI9   = ZERO
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = ZERO
      PSI14  = ZERO
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A9  = XK17 *(WATER/GAMA(17))**3.0
      A13 = XK19 *(WATER/GAMA(19))**2.0
      A14 = XK20 *(WATER/GAMA(20))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0



      PSI5 = (CHI5-PSI2)*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) &
             - A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
             2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN (MAX (PSI5, TINY) , CHI5)

      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI13.GT.TINY .AND. WATER.GT.TINY) THEN          
         VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
         GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
         DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
         PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
         PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
      ENDIF

      IF (CHI14.GT.TINY .AND. WATER.GT.TINY) THEN          
         PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
                 PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
         PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
      ENDIF

      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN          
         BBP = PSI10+PSI13+PSI14
         CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
         DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
      CALL POLY32p1 (BBP, CCP, DDP, PSI9, ISLV)
        IF (ISLV.EQ.0) THEN
            PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
        ELSE
            PSI9 = ZERO
        ENDIF
      ENDIF

      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     
         VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
         GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
         DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
         PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF

      IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     




          PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
                 PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
          PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
      ENDIF




      MOLAL (2) = PSI8 + PSI7                                     
      MOLAL (3) = PSI4                                            
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   
      MOLAL (5) = PSI9 + PSI10                                    
      MOLAL (6) = ZERO                                            
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   
      MOLAL (8) = PSI11 + PSI12 + PSI17                           
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       
      MOLAL (10)= PSI10 + PSI15 + PSI16                           


























      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                  - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH2p1 (SMIN, HI, OHI)
      MOLAL (1) = HI


      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)



      CNACL     = MAX (CHI7 - PSI7, ZERO)
      CNANO3    = MAX (CHI8 - PSI8, ZERO)
      CK2SO4    = MAX (CHI9 - PSI9, ZERO)
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = MAX (CHI13 - PSI13, ZERO)
      CKCL      = MAX (CHI14 - PSI14, ZERO)
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO



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
      PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)



      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3



      A2   = XK10 /(R*TEMP*R*TEMP)
      IF (GNH3*GHNO3.GT.A2) THEN
         DELT = MIN(GNH3, GHNO3)
         BB = -(GNH3+GHNO3)
         CC = GNH3*GHNO3-A2
         DD = BB*BB - 4.D0*CC
         PSI21 = 0.5D0*(-BB + SQRT(DD))
         PSI22 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI21.GT.ZERO .AND. PSI21.GT.ZERO) THEN
            PSI2 = PSI21
         ELSEIF (DELT-PSI22.GT.ZERO .AND. PSI22.GT.ZERO) THEN
            PSI2 = PSI22
         ELSE
            PSI2 = ZERO
         ENDIF
      ELSE
         PSI2 = ZERO
      ENDIF
      PSI2 = MAX(MIN(MIN(PSI2,CHI4-PSI4-PSI3),CHI5-PSI5),ZERO)



      GNH3    = MAX (GNH3 - PSI2, TINY)
      GHCL    = MAX (GHNO3 - PSI2, TINY)
      CNH4NO3 = PSI2

      CALL CALCMR2p1                                    



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE




20    FUNCP2A2p1 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      RETURN



      END
























      SUBROUTINE CALCP12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCP1A2p1, CALCP2A2p1



      IF (RH.LT.DRMP1) THEN
         SCASE = 'P1 ; SUBCASE 1'
         CALL CALCP1A2p1              
         SCASE = 'P1 ; SUBCASE 1'
      ELSE
         SCASE = 'P1 ; SUBCASE 2'  
         CALL CALCMDRH22p1 (RH, DRMP1, DRCACL2, CALCP1A2p1, CALCP2A2p1)
         SCASE = 'P1 ; SUBCASE 2'
      ENDIF


      RETURN



      END



















      SUBROUTINE CALCP1A2p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR, &
                       NO3FR



      CCASO4  = MIN (W(2), W(6))                    
      CAFR    = MAX (W(6) - CCASO4, ZERO)
      SO4FR   = MAX (W(2) - CCASO4, ZERO)
      CK2SO4  = MIN (SO4FR, 0.5D0*W(7))             
      FRK     = MAX (W(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
      CMGSO4  = SO4FR                               
      FRMG    = MAX (W(8) - CMGSO4, ZERO)
      CNACL   = MIN (W(1), W(5))                    
      NAFR    = MAX (W(1) - CNACL, ZERO)
      CLFR    = MAX (W(5) - CNACL, ZERO)
      CCANO32 = MIN (CAFR, 0.5D0*W(4))              
      CAFR    = MAX (CAFR - CCANO32, ZERO)
      NO3FR   = MAX (W(4) - 2.D0*CCANO32, ZERO)
      CCACL2  = MIN (CAFR, 0.5D0*CLFR)              
      CAFR    = MAX (CAFR - CCACL2, ZERO)
      CLFR    = MAX (CLFR - 2.D0*CCACL2, ZERO)
      CMGNO32 = MIN (FRMG, 0.5D0*NO3FR)             
      FRMG    = MAX (FRMG - CMGNO32, ZERO)
      NO3FR   = MAX (NO3FR - 2.D0*CMGNO32, ZERO)
      CMGCL2  = MIN (FRMG, 0.5D0*CLFR)              
      FRMG    = MAX (FRMG - CMGCL2, ZERO)
      CLFR    = MAX (CLFR - 2.D0*CMGCL2, ZERO)
      CNANO3  = MIN (NAFR, NO3FR)                   
      NAFR    = MAX (NAFR - CNANO3, ZERO)
      NO3FR   = MAX (NO3FR - CNANO3, ZERO)
      CKCL    = MIN (FRK, CLFR)                     
      FRK     = MAX (FRK - CKCL, ZERO)
      CLFR    = MAX (CLFR - CKCL, ZERO)
      CKNO3   = MIN (FRK, NO3FR)                    
      FRK     = MAX (FRK - CKNO3, ZERO)
      NO3FR   = MAX (NO3FR - CKNO3, ZERO)



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
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND. &
             BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND. &
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



      END



















      SUBROUTINE CALCL92p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,       &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,       &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,     &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,         &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,   &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCL1A2p1



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = CNA2SO4
      PSI5 = CNH42S4
      PSI6 = CK2SO4
      PSI7 = CMGSO4
      PSI8 = CKHSO4

      CALAOU = .TRUE.              
      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A9 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.



      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9              
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = MAX(BB*BB - 4.D0*CC, ZERO)
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)



      MOLAL(1) = LAMDA                                            
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         
      MOLAL(6) = PSI2 + PSI3 + PSI1 + PSI8 - LAMDA                
      MOLAL(9) = PSI8 + 2.0D0*PSI6                                
      MOLAL(10)= PSI7                                             

      CLC      = ZERO
      CNAHSO4  = ZERO
      CNA2SO4  = ZERO
      CNH42S4  = ZERO
      CNH4HS4  = ZERO
      CK2SO4   = ZERO
      CMGSO4   = ZERO
      CKHSO4   = ZERO

      CALL CALCMR2p1                                         




      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE

20    RETURN



      END


















      SUBROUTINE CALCL82p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCL1A2p1



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = CNA2SO4
      PSI5 = CNH42S4
      PSI6 = ZERO
      PSI7 = CMGSO4
      PSI8 = CKHSO4

      CALAOU = .TRUE.              
      PSI6LO = ZERO                
      PSI6HI = CHI6                



       IF (CHI6.LE.TINY) THEN
         Y1 = FUNCL82p1 (ZERO)
         GOTO 50
      ENDIF

      X1 = PSI6HI
      Y1 = FUNCL82p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCL82p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCL82p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCL8')    
         GOTO 50
      ENDIF


20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL82p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCL8')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL82p1 (X3)

50    RETURN



      END




















      DOUBLE PRECISION FUNCTION FUNCL82p1 (P6)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI6   = P6



      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)



      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9              
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = BB*BB - 4.D0*CC
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)



      MOLAL(1) = LAMDA                                            
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     
      MOLAL(9) = PSI8 + 2.0*PSI6                                  
      MOLAL(10)= PSI7                                             

      CLC      = ZERO
      CNAHSO4  = ZERO
      CNA2SO4  = ZERO
      CNH42S4  = ZERO
      CNH4HS4  = ZERO
      CK2SO4   = MAX(CHI6 - PSI6, ZERO)
      CMGSO4   = ZERO
      CKHSO4   = ZERO
      CALL CALCMR2p1                                       



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A6 = XK17*(WATER/GAMA(17))**3.0
      FUNCL82p1 = MOLAL(9)*MOLAL(9)*MOLAL(5)/A6 - ONE
      RETURN



      END



















      SUBROUTINE CALCL72p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,       &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,       &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,     &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,         &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,   &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCL1A2p1



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = CNH42S4
      PSI6 = ZERO
      PSI7 = CMGSO4
      PSI8 = CKHSO4

      CALAOU = .TRUE.              
      PSI4LO = ZERO                
      PSI4HI = CHI4                



       IF (CHI4.LE.TINY) THEN
         Y1 = FUNCL72p1 (ZERO)
         GOTO 50
      ENDIF

      X1 = PSI4HI
      Y1 = FUNCL72p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCL72p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCL72p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCL7')    
         GOTO 50
      ENDIF


20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL72p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCL7')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL72p1 (X3)

50    RETURN



      END




















      DOUBLE PRECISION FUNCTION FUNCL72p1 (P4)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI4   = P4



      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4 = XK5 *(WATER/GAMA(2))**3.0
      A6 = XK17*(WATER/GAMA(17))**3.0
      A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)






      IF (CHI6.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
         BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
         CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
         CALL POLY32p1 (AA, BB, CC, PSI6, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI6 = MIN (PSI6, CHI6)
         ELSE
            PSI6 = ZERO
         ENDIF
      ENDIF

      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9              
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = BB*BB - 4.D0*CC
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)




      MOLAL(1) = LAMDA                                            
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     
      MOLAL(9) = PSI8 + 2.0*PSI6                                  
      MOLAL(10)= PSI7                                             

      CLC      = ZERO
      CNAHSO4  = ZERO
      CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
      CNH42S4  = ZERO
      CNH4HS4  = ZERO
      CK2SO4   = MAX(CHI6 - PSI6, ZERO)
      CMGSO4   = ZERO
      CKHSO4   = ZERO
      CALL CALCMR2p1                                       



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0
      FUNCL72p1 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END




















      SUBROUTINE CALCL62p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCL1A2p1



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = CNH42S4
      PSI6 = ZERO
      PSI7 = ZERO
      PSI8 = CKHSO4

      CALAOU = .TRUE.              
      PSI4LO = ZERO                
      PSI4HI = CHI4                



       IF (CHI4.LE.TINY) THEN
         Y1 = FUNCL62p1 (ZERO)
         GOTO 50
      ENDIF

      X1 = PSI4HI
      Y1 = FUNCL62p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCL62p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCL62p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCL6')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL62p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCL6')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL62p1 (X3)

50    RETURN



      END



















      DOUBLE PRECISION FUNCTION FUNCL62p1 (P4)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI4   = P4



      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4 = XK5*(WATER/GAMA(2))**3.0
      A6 = XK17*(WATER/GAMA(17))**3.0
      A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)






      IF (CHI6.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
         BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
         CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
         CALL POLY32p1 (AA, BB, CC, PSI6, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI6 = MIN (PSI6, CHI6)
         ELSE
            PSI6 = ZERO
         ENDIF
      ENDIF

      PSI7 = CHI7

      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = BB*BB - 4.D0*CC
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)



      MOLAL(1) = LAMDA                                            
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     
      MOLAL(9) = PSI8 + 2.0*PSI6                                  
      MOLAL(10)= PSI7                                             

      CLC      = ZERO
      CNAHSO4  = ZERO
      CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
      CNH42S4  = ZERO
      CNH4HS4  = ZERO
      CK2SO4   = MAX(CHI6 - PSI6, ZERO)
      CMGSO4   = ZERO
      CKHSO4   = ZERO
      CALL CALCMR2p1                                       



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4 = XK5 *(WATER/GAMA(2))**3.0
      FUNCL62p1 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END



















      SUBROUTINE CALCL52p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCL1A2p1



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = CNH42S4
      PSI6 = ZERO
      PSI7 = ZERO
      PSI8 = ZERO

      CALAOU = .TRUE.              
      PSI4LO = ZERO                
      PSI4HI = CHI4                




      IF (CHI4.LE.TINY) THEN
         Y1 = FUNCL52p1 (ZERO)
         GOTO 50
      ENDIF

      X1 = PSI4HI
      Y1 = FUNCL52p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50




      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI4LO)
         Y2 = FUNCL52p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCL52p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCL5')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL52p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCL5')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL52p1 (X3)

50    RETURN



      END




















      DOUBLE PRECISION FUNCTION FUNCL52p1 (P4)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI4   = P4



      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4 = XK5*(WATER/GAMA(2))**3.0
      A6 = XK17*(WATER/GAMA(17))**3.0
      A8 = XK18*(WATER/GAMA(18))**2.0
      A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)






      IF (CHI6.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
         BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
         CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
         CALL POLY32p1 (AA, BB, CC, PSI6, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI6 = MIN (PSI6, CHI6)
         ELSE
            PSI6 = ZERO
         ENDIF
      ENDIF

      PSI7 = CHI7

      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = MAX(BB*BB - 4.D0*CC, ZERO)
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)

      BITA = PSI3 + PSI2 + PSI1 + 2.D0*PSI6 - LAMDA
      CAMA = 2.D0*PSI6*(PSI3 + PSI2 + PSI1 - LAMDA) - A8
      DELT  = MAX(BITA*BITA - 4.D0*CAMA, ZERO)
      PSI8 = 0.5D0*(-BITA + SQRT(DELT))
      PSI8 = MIN(MAX (PSI8, ZERO), CHI8)



      MOLAL(1) = LAMDA                                            
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     
      MOLAL(9) = PSI8 + 2.0D0*PSI6                                
      MOLAL(10)= PSI7                                             

      CLC      = ZERO
      CNAHSO4  = ZERO
      CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
      CNH42S4  = ZERO
      CNH4HS4  = ZERO
      CK2SO4   = MAX(CHI6 - PSI6, ZERO)
      CMGSO4   = ZERO
      CKHSO4   = MAX(CHI8 - PSI8, ZERO)

      CALL CALCMR2p1                                       




      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0
      FUNCL52p1 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE

      RETURN



      END



















      SUBROUTINE CALCL42p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCL1A2p1



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4

      PSI1 = CNH4HS4               
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = ZERO
      PSI6 = ZERO
      PSI7 = ZERO
      PSI8 = ZERO

      CALAOU = .TRUE.              
      PSI4LO = ZERO                
      PSI4HI = CHI4                

      IF (CHI4.LE.TINY) THEN
         Y1 = FUNCL42p1 (ZERO)
         GOTO 50
      ENDIF



      X1 = PSI4HI
      Y1 = FUNCL42p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCL42p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCL42p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCL4')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL42p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCL4')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL42p1 (X3)

50    RETURN



      END




















      DOUBLE PRECISION FUNCTION FUNCL42p1 (P4)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI4   = P4



      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4 = XK5*(WATER/GAMA(2))**3.0
      A5 = XK7*(WATER/GAMA(4))**3.0
      A6 = XK17*(WATER/GAMA(17))**3.0
      A8 = XK18*(WATER/GAMA(18))**2.0
      A9 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0



      PSI5 = (PSI3 + 2.D0*PSI4 - SQRT(A4/A5)*(3.D0*PSI2 + PSI1))  &
              /2.D0/SQRT(A4/A5)
      PSI5 = MAX (MIN (PSI5, CHI5), ZERO)

      PSI7 = CHI7

      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = MAX(BB*BB - 4.D0*CC, ZERO)
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)




      IF (CHI6.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
         BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
         CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
         CALL POLY32p1 (AA, BB, CC, PSI6, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI6 = MIN (PSI6, CHI6)
         ELSE
            PSI6 = ZERO
         ENDIF
      ENDIF

      BITA = PSI3 + PSI2 + PSI1 + 2.D0*PSI6 - LAMDA
      CAMA = 2.D0*PSI6*(PSI3 + PSI2 + PSI1 - LAMDA) - A8
      DELT  = MAX(BITA*BITA - 4.D0*CAMA, ZERO)
      PSI8 = 0.5D0*(-BITA + SQRT(DELT))
      PSI8 = MIN(MAX (PSI8, ZERO), CHI8)



      MOLAL(1) = LAMDA                                            
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     
      MOLAL(9) = PSI8 + 2.0D0*PSI6                                
      MOLAL(10)= PSI7                                             

      CLC      = ZERO
      CNAHSO4  = ZERO
      CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
      CNH42S4  = MAX(CHI5 - PSI5, ZERO)
      CNH4HS4  = ZERO
      CK2SO4   = MAX(CHI6 - PSI6, ZERO)
      CMGSO4   = ZERO
      CKHSO4   = MAX(CHI8 - PSI8, ZERO)
      CALL CALCMR2p1                                       



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0
      FUNCL42p1 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END

























      SUBROUTINE CALCL32p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCL1A2p1, CALCL42p1



      CALL CALCL1A2p1



      IF (CNH4HS4.GT.TINY .OR. CNAHSO4.GT.TINY) THEN
         SCASE = 'L3 ; SUBCASE 1'
         CALL CALCL3A2p1                     
         SCASE = 'L3 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRML3) THEN         
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCL1A2p1
            SCASE = 'L3 ; SUBCASE 2'

         ELSEIF (RH.GE.DRML3) THEN     
            SCASE = 'L3 ; SUBCASE 3'
            CALL CALCMDRH22p1 (RH, DRML3, DRLC, CALCL1A2p1, CALCL42p1)
            SCASE = 'L3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCL3A2p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCL1A2p1



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4

      PSI1 = CNH4HS4               
      PSI2 = ZERO
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = ZERO
      PSI6 = ZERO
      PSI7 = ZERO
      PSI8 = ZERO

      CALAOU = .TRUE.              
      PSI2LO = ZERO                
      PSI2HI = CHI2                



      X1 = PSI2HI
      Y1 = FUNCL3A2p1 (X1)
      YHI= Y1                      



      IF (YHI.LT.EPS) GOTO 50



      DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = FUNCL3A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCL3A2p1 (ZERO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL3A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCL3A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL3A2p1 (X3)

50    RETURN



      END



















      DOUBLE PRECISION FUNCTION FUNCL3A2p1 (P2)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17




      PSI2   = P2                  
      PSI4LO = ZERO                
      PSI4HI = CHI4                



      IF (CHI4.LE.TINY) THEN
         FUNCL3A2p1 = FUNCL3B2p1 (ZERO)
         GOTO 50
      ENDIF



      X1 = PSI4HI
      Y1 = FUNCL3B2p1 (X1)
      IF (ABS(Y1).LE.EPS) GOTO 50
      YHI= Y1                      



      IF (YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI4LO)
         Y2 = FUNCL3B2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCL3B2p1 (PSI4LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL3B2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0004, 'FUNCL3A2p1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL3B2p1 (X3)



50    A2      = XK13*(WATER/GAMA(13))**5.0
      FUNCL3A2p1 = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.0/A2 - ONE
      RETURN



      END




















      DOUBLE PRECISION FUNCTION FUNCL3B2p1 (P4)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI4   = P4

      FRST   = .TRUE.
      CALAIN = .TRUE.



      DO 10 I=1,NSWEEP

      A4 = XK5*(WATER/GAMA(2))**3.0
      A5 = XK7*(WATER/GAMA(4))**3.0
      A6 = XK17*(WATER/GAMA(17))**3.0
      A8 = XK18*(WATER/GAMA(18))**2.0
      A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)



      PSI5 = (PSI3 + 2.D0*PSI4 - SQRT(A4/A5)*(3.D0*PSI2 + PSI1)) & 
              /2.D0/SQRT(A4/A5)
      PSI5 = MAX (MIN (PSI5, CHI5), ZERO)

      PSI7 = CHI7

      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = MAX(BB*BB - 4.D0*CC, ZERO)
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)




      IF (CHI6.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
         BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
         CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
         CALL POLY32p1 (AA, BB, CC, PSI6, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI6 = MIN (PSI6, CHI6)
         ELSE
            PSI6 = ZERO
         ENDIF
      ENDIF

      BITA = PSI3 + PSI2 + PSI1 + 2.D0*PSI6 - LAMDA
      CAMA = 2.D0*PSI6*(PSI3 + PSI2 + PSI1 - LAMDA) - A8
      DELT  = MAX(BITA*BITA - 4.D0*CAMA, ZERO)
      PSI8 = 0.5D0*(-BITA + SQRT(DELT))
      PSI8 = MIN(MAX (PSI8, ZERO), CHI8)



      MOLAL(1) = LAMDA                                            
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     
      MOLAL(9) = PSI8 + 2.0D0*PSI6                                
      MOLAL(10)= PSI7                                             

      CLC      = MAX(CHI2 - PSI2, ZERO)
      CNAHSO4  = ZERO
      CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
      CNH42S4  = MAX(CHI5 - PSI5, ZERO)
      CNH4HS4  = ZERO
      CK2SO4   = MAX(CHI6 - PSI6, ZERO)
      CMGSO4   = MAX(CHI7 - PSI7, ZERO)
      CKHSO4   = MAX(CHI8 - PSI8, ZERO)
      CALL CALCMR2p1                                       



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0
      FUNCL3B2p1 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END


























      SUBROUTINE CALCL22p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCL1A2p1, CALCL3A2p1



      CALL CALCL1A2p1



      IF (CNH4HS4.GT.TINY) THEN
         SCASE = 'L2 ; SUBCASE 1'
         CALL CALCL2A2p1
         SCASE = 'L2 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRML2) THEN         
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCL1A2p1
            SCASE = 'L2 ; SUBCASE 2'

         ELSEIF (RH.GE.DRML2) THEN     
            SCASE = 'L2 ; SUBCASE 3'
            CALL CALCMDRH22p1 (RH, DRML2, DRNAHSO4, CALCL1A2p1, CALCL3A2p1)
            SCASE = 'L2 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCL2A2p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,       &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,       &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,     &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,         &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,   &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CHI1 = CNH4HS4               
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4


      PSI1 = CNH4HS4               
      PSI2 = ZERO
      PSI3 = ZERO
      PSI4 = ZERO
      PSI5 = ZERO
      PSI6 = ZERO
      PSI7 = ZERO
      PSI8 = ZERO

      CALAOU = .TRUE.              
      PSI2LO = ZERO                
      PSI2HI = CHI2                



      X1 = PSI2HI
      Y1 = FUNCL2A2p1 (X1)
      YHI= Y1                      



      IF (YHI.LT.EPS) GOTO 50



      DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = FUNCL2A2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCL2A2p1 (ZERO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL2A2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCL2A')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL2A2p1 (X3)

50    RETURN



      END



















      DOUBLE PRECISION FUNCTION FUNCL2A2p1 (P2)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17




      PSI2   = P2                  
      PSI4LO = ZERO                
      PSI4HI = CHI4                




      IF (CHI4.LE.TINY) THEN
         FUNCL2A2p1 = FUNCL2B2p1 (ZERO)
         GOTO 50
      ENDIF




      X1 = PSI4HI
      Y1 = FUNCL2B2p1 (X1)

      IF (ABS(Y1).LE.EPS) GOTO 50
      YHI= Y1                      



      IF (YHI.LT.ZERO) GOTO 50



      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI4LO)
         Y2 = FUNCL2B2p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      IF (Y2.GT.EPS) Y2 = FUNCL2B2p1 (PSI4LO)
      GOTO 50



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCL2B2p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0004, 'FUNCL2A2p1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCL2B2p1 (X3)



50    A2      = XK13*(WATER/GAMA(13))**5.0
      FUNCL2A2p1 = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.0/A2 - ONE
      RETURN



      END



















      DOUBLE PRECISION FUNCTION FUNCL2B2p1 (P4)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      PSI4   = P4                  



      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI3   = CHI3
      PSI5   = CHI5
      LAMDA  = ZERO
      PSI6   = CHI6
      PSI8   = CHI8



      DO 10 I=1,NSWEEP

      A3 = XK11*(WATER/GAMA(12))**2.0
      A4 = XK5*(WATER/GAMA(2))**3.0
      A5 = XK7*(WATER/GAMA(4))**3.0
      A6 = XK17*(WATER/GAMA(17))**3.0
      A8 = XK18*(WATER/GAMA(18))**2.0
      A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)



      PSI5 = (PSI3 + 2.D0*PSI4 - SQRT(A4/A5)*(3.D0*PSI2 + PSI1)) &
              /2.D0/SQRT(A4/A5)
      PSI5 = MAX (MIN (PSI5, CHI5), ZERO)

      IF (CHI3.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = 2.D0*PSI4 + PSI2 + PSI1 + PSI8 - LAMDA
         BB   = 2.D0*PSI4*(PSI2 + PSI1 + PSI8 - LAMDA) - A3
         CC   = ZERO
         CALL POLY32p1 (AA, BB, CC, PSI3, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI3 = MIN (PSI3, CHI3)
         ELSE
            PSI3 = ZERO
         ENDIF
      ENDIF

      PSI7 = CHI7

      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = MAX(BB*BB - 4.D0*CC, ZERO)
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)




      IF (CHI6.GT.TINY .AND. WATER.GT.TINY) THEN
         AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
         BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
         CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
         CALL POLY32p1 (AA, BB, CC, PSI6, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI6 = MIN (PSI6, CHI6)
         ELSE
            PSI6 = ZERO
         ENDIF
      ENDIF

      BITA = PSI3 + PSI2 + PSI1 + 2.D0*PSI6 - LAMDA              
      CAMA = 2.D0*PSI6*(PSI3 + PSI2 + PSI1 - LAMDA) - A8
      DELT  = MAX(BITA*BITA - 4.D0*CAMA, ZERO)
      PSI8 = 0.5D0*(-BITA + SQRT(DELT))
      PSI8 = MIN(MAX (PSI8, ZERO), CHI8)



      MOLAL(1) = LAMDA                                            
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     
      MOLAL(9) = PSI8 + 2.0D0*PSI6                                
      MOLAL(10)= PSI7                                             

      CLC      = MAX(CHI2 - PSI2, ZERO)
      CNAHSO4  = MAX(CHI3 - PSI3, ZERO)
      CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
      CNH42S4  = MAX(CHI5 - PSI5, ZERO)
      CNH4HS4  = ZERO
      CK2SO4   = MAX(CHI6 - PSI6, ZERO)
      CMGSO4   = MAX(CHI7 - PSI7, ZERO)
      CKHSO4   = MAX(CHI8 - PSI8, ZERO)
      CALL CALCMR2p1                                       



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    A4     = XK5 *(WATER/GAMA(2))**3.0
      FUNCL2B2p1 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN



      END






















      SUBROUTINE CALCL12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCL1A2p1, CALCL2A2p1



      IF (RH.LT.DRML1) THEN
         SCASE = 'L1 ; SUBCASE 1'
         CALL CALCL1A2p1              
         SCASE = 'L1 ; SUBCASE 1'
      ELSE
         SCASE = 'L1 ; SUBCASE 2'  
         CALL CALCMDRH22p1 (RH, DRML1, DRNH4HS4, CALCL1A2p1, CALCL2A2p1)
         SCASE = 'L1 ; SUBCASE 2'
      ENDIF





      RETURN



      END


















      SUBROUTINE CALCL1A2p1
      INCLUDE 'module_isrpia_inc.F'



      CCASO4  = MIN (W(6), W(2))                    
      FRSO4   = MAX(W(2) - CCASO4, ZERO)
      CAFR    = MAX(W(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*W(7), FRSO4)             
      FRK     = MAX(W(7) - 2.D0*CK2SO4, ZERO)
      FRSO4   = MAX(FRSO4 - CK2SO4, ZERO)
      CNA2SO4 = MIN (0.5D0*W(1), FRSO4)             
      FRNA    = MAX(W(1) - 2.D0*CNA2SO4, ZERO)
      FRSO4   = MAX(FRSO4 - CNA2SO4, ZERO)
      CMGSO4  = MIN (W(8), FRSO4)                   
      FRMG    = MAX(W(8) - CMGSO4, ZERO)
      FRSO4   = MAX(FRSO4 - CMGSO4, ZERO)

      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      CKHSO4  = ZERO

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
            FRSO4  = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
         IF (CK2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CKHSO4 = 2.D0*FRSO4
            CK2SO4 = MAX(CK2SO4-FRSO4, ZERO)
       ENDIF
      ENDIF



      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO

      RETURN



      END



















      SUBROUTINE CALCK42p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEK2p1/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
                     A1,   A2,   A3,   A4



      CALAOU =.TRUE.               
      FRST   = .TRUE.
      CALAIN = .TRUE.

      CHI1   = W(3)                
      CHI2   = W(1)                
      CHI3   = W(7)                
      CHI4   = W(8)                

      LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  
      PSI1   = CHI1                            
      PSI2   = CHI2                            
      PSI3   = CHI3                            
      PSI4   = CHI4                            



      DO 10 I=1,NSWEEP

      A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0

      BB   = A4+LAMDA+PSI4                               
      CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
      DD   = MAX(BB*BB-4.D0*CC, ZERO)
      KAPA = 0.5D0*(-BB+SQRT(DD))



      MOLAL (1) = MAX(LAMDA + KAPA, TINY)                         
      MOLAL (2) = PSI2                                            
      MOLAL (3) = PSI1                                            
      MOLAL (5) = MAX(KAPA + PSI4, ZERO)                          
      MOLAL (6) = MAX(LAMDA + PSI1 + PSI2 + PSI3 - KAPA, ZERO)    
      MOLAL (9) = PSI3                                            
      MOLAL (10)= PSI4                                            

      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CKHSO4  = ZERO
      CCASO4  = W(6)
      CMGSO4  = ZERO

      CALL CALCMR2p1                                      



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE

20    RETURN



      END


















      SUBROUTINE CALCK32p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEK2p1/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
                     A1,   A2,   A3,   A4



      CALAOU =.TRUE.               
      CHI1   = W(3)                
      CHI2   = W(1)                
      CHI3   = W(7)                
      CHI4   = W(8)                

      PSI3LO = TINY                
      PSI3HI = CHI3                



      X1 = PSI3HI
      Y1 = FUNCK32p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI3HI-PSI3LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCK32p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCK32p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCK3')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCK32p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCK3')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCK32p1 (X3)

50    RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCK32p1 (P1)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEK2p1/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
                     A1,   A2,   A3,   A4



      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  
      PSI3   = P1
      PSI1   = CHI1                             
      PSI2   = CHI2                             
      PSI4   = CHI4                             




      DO 10 I=1,NSWEEP

      A3 = XK18 *(WATER/GAMA(18))**2.0
      A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0


      BB   = A4+LAMDA+PSI4                             
      CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
      DD   = MAX(BB*BB-4.D0*CC, ZERO)
      KAPA = 0.5D0*(-BB+SQRT(DD))



      MOLAL (1) = MAX(LAMDA + KAPA, ZERO)                
      MOLAL (2) = PSI2                                   
      MOLAL (3) = PSI1                                   
      MOLAL (4) = ZERO
      MOLAL (5) = MAX(KAPA + PSI4, ZERO)                 
      MOLAL (6) = MAX(LAMDA+PSI1+PSI2+PSI3-KAPA,ZERO)    
      MOLAL (7) = ZERO
      MOLAL (8) = ZERO
      MOLAL (9) = PSI3                                   
      MOLAL (10)= PSI4

      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CKHSO4  = CHI3-PSI3
      CCASO4  = W(6)
      CMGSO4  = ZERO

      CALL CALCMR2p1                                      



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCK32p1 = MOLAL(9)*MOLAL(6)/A3 - ONE



      END


















      SUBROUTINE CALCK22p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEK2p1/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
                     A1,   A2,   A3,   A4



      CALAOU =.TRUE.               
      CHI1   = W(3)                
      CHI2   = W(1)                
      CHI3   = W(7)                
      CHI4   = W(8)                

      PSI3LO = TINY                
      PSI3HI = CHI3                



      X1 = PSI3HI
      Y1 = FUNCK22p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI3HI-PSI3LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCK22p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCK22p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCK2')    
         GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCK22p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCK2')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCK22p1 (X3)

50    RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCK22p1 (P1)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEK2p1/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
                     A1,   A2,   A3,   A4



      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  
      PSI3   = P1
      PSI1   = CHI1                              
      PSI4   = CHI4                              



      DO 10 I=1,NSWEEP

      A2 = XK11 *(WATER/GAMA(12))**2.0
      A3 = XK18 *(WATER/GAMA(18))**2.0
      A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0

      PSI2 = A2/A3*PSI3                                   
      PSI2 = MIN(MAX(PSI2, ZERO),CHI2)

      BB   = A4+LAMDA+PSI4                                
      CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
      DD   = MAX(BB*BB-4.D0*CC, ZERO)
      KAPA = 0.5D0*(-BB+SQRT(DD))



      MOLAL (1) = MAX(LAMDA + KAPA, ZERO)                
      MOLAL (2) = PSI2                                   
      MOLAL (3) = PSI1                                   
      MOLAL (4) = ZERO
      MOLAL (5) = MAX(KAPA + PSI4, ZERO)                 
      MOLAL (6) = MAX(LAMDA+PSI1+PSI2+PSI3-KAPA,ZERO)    
      MOLAL (7) = ZERO
      MOLAL (8) = ZERO
      MOLAL (9) = PSI3                                   
      MOLAL (10)= PSI4

      CNH4HS4 = ZERO
      CNAHSO4 = CHI2-PSI2
      CKHSO4  = CHI3-PSI3
      CCASO4  = W(6)
      CMGSO4  = ZERO

      CALL CALCMR2p1                                      



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCK22p1 = MOLAL(9)*MOLAL(6)/A3 - ONE



      END


















      SUBROUTINE CALCK12p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEK2p1/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
                     A1,   A2,   A3,   A4




      CALAOU =.TRUE.               
      CHI1   = W(3)                
      CHI2   = W(1)                
      CHI3   = W(7)                
      CHI4   = W(8)                

      PSI3LO = TINY                
      PSI3HI = CHI3                



      X1 = PSI3HI
      Y1 = FUNCK12p1 (X1)
      YHI= Y1                      



      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50



      DX = (PSI3HI-PSI3LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = FUNCK12p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = FUNCK12p1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN       
         GOTO 50
      ELSE
        CALL PUSHERR2p1 (0001, 'CALCK1')    
        GOTO 50
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCK12p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCK1')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCK12p1 (X3)

50    RETURN



      END


















      DOUBLE PRECISION FUNCTION FUNCK12p1 (P1)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEK2p1/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
                     A1,   A2,   A3,   A4



      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  
      PSI3   = P1
      PSI4   = CHI4                                    



      DO 10 I=1,NSWEEP

      A1 = XK12 *(WATER/GAMA(09))**2.0
      A2 = XK11 *(WATER/GAMA(12))**2.0
      A3 = XK18 *(WATER/GAMA(18))**2.0
      A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0

      PSI1 = A1/A3*PSI3                                   
      PSI1 = MIN(MAX(PSI1, ZERO),CHI1)

      PSI2 = A2/A3*PSI3                                   
      PSI2 = MIN(MAX(PSI2, ZERO),CHI2)

      BB   = A4+LAMDA+PSI4                                
      CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
      DD   = MAX(BB*BB-4.D0*CC, ZERO)
      KAPA = 0.5D0*(-BB+SQRT(DD))



      MOLAL (1) = MAX(LAMDA + KAPA, ZERO)              
      MOLAL (2) = PSI2                                 
      MOLAL (3) = PSI1                                 
      MOLAL (4) = ZERO                                 
      MOLAL (5) = MAX(KAPA + PSI4, ZERO)               
      MOLAL (6) = MAX(LAMDA+PSI1+PSI2+PSI3-KAPA,ZERO)  
      MOLAL (7) = ZERO                                 
      MOLAL (8) = ZERO                                 
      MOLAL (9) = PSI3                                 
      MOLAL (10)= PSI4                                 

      CNH4HS4 = CHI1-PSI1
      CNAHSO4 = CHI2-PSI2
      CKHSO4  = CHI3-PSI3
      CCASO4  = W(6)
      CMGSO4  = ZERO

      CALL CALCMR2p1                                      



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1

      ELSE
         GOTO 20
      ENDIF
10    CONTINUE



20    FUNCK12p1 = MOLAL(9)*MOLAL(6)/A3 - ONE



      END


