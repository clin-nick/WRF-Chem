















      SUBROUTINE ISRP1R2p1 (WI, RHI, TEMPI)
      INCLUDE 'module_isrpia_inc.F'
      DIMENSION WI(NCOMP)



      CALL INIT12p1 (WI, RHI, TEMPI)



      IF (RH.GE.DRNH42S4) THEN         
         SULRATW = GETASR2p1(WAER(2), RHI)     
      ELSE
         SULRATW = 2.0D0                    
      ENDIF
      SULRAT  = WAER(3)/WAER(2)         





      IF (SULRATW.LE.SULRAT) THEN

      IF(METSTBL.EQ.1) THEN
         SCASE = 'S2'
         CALL CALCS22p1                 
      ELSE

         IF (RH.LT.DRNH42S4) THEN    
            SCASE = 'S1'
            CALL CALCS12p1              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'S2'
            CALL CALCS22p1              
         ENDIF
      ENDIF



      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN
      W(2) = WAER(2)
      W(3) = WAER(3)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL CALCB42p1                 
         SCASE = 'B4'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'B1'
            CALL CALCB12p1              
            SCASE = 'B1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'B2'
            CALL CALCB22p1              
            SCASE = 'B2'

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'B3'
            CALL CALCB32p1              
            SCASE = 'B3'

         ELSEIF (DRNH42S4.LE.RH) THEN         
            SCASE = 'B4'
            CALL CALCB42p1              
            SCASE = 'B4'
         ENDIF
      ENDIF

      CALL CALCNH3P2p1          



      ELSEIF (SULRAT.LT.1.0) THEN             
      W(2) = WAER(2)
      W(3) = WAER(3)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC22p1                 
         SCASE = 'C2'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'C1'
            CALL CALCC12p1              
            SCASE = 'C1'

         ELSEIF (DRNH4HS4.LE.RH) THEN         
            SCASE = 'C2'
            CALL CALCC22p1              
            SCASE = 'C2'
         ENDIF
      ENDIF

      CALL CALCNH3P2p1

      ENDIF
      RETURN



      END
















      SUBROUTINE ISRP2R2p1 (WI, RHI, TEMPI)
      INCLUDE 'module_isrpia_inc.F'
      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ



      TRYLIQ = .TRUE.             

10    CALL INIT22p1 (WI, RHI, TEMPI)



      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN 
         SULRATW = GETASR2p1(WAER(2), RHI)     
      ELSE
         SULRATW = 2.0D0                    
      ENDIF
      SULRAT = WAER(3)/WAER(2)





      IF (SULRATW.LE.SULRAT) THEN                

      IF(METSTBL.EQ.1) THEN
         SCASE = 'N3'
         CALL CALCN32p1                 
      ELSE

         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'N1'
            CALL CALCN12p1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'N2'
            CALL CALCN22p1              

         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'N3'
            CALL CALCN32p1              
         ENDIF
      ENDIF







      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
      W(2) = WAER(2)
      W(3) = WAER(3)
      W(4) = WAER(4)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL CALCB42p1                 
         SCASE = 'B4'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'B1'
            CALL CALCB12p1              
            SCASE = 'B1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'B2'
            CALL CALCB22p1              
            SCASE = 'B2'

         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'B3'
            CALL CALCB32p1              
            SCASE = 'B3'

         ELSEIF (DRNH42S4.LE.RH) THEN         
            SCASE = 'B4'
            CALL CALCB42p1              
            SCASE = 'B4'
         ENDIF
      ENDIF



      MOLAL(7) = WAER(4)             
      MOLAL(1) = MOLAL(1) + WAER(4)  
      CALL CALCNAP2p1            
      CALL CALCNH3P2p1







      ELSEIF (SULRAT.LT.1.0) THEN             
      W(2) = WAER(2)
      W(3) = WAER(3)
      W(4) = WAER(4)

      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC22p1                 
         SCASE = 'C2'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'C1'
            CALL CALCC12p1              
            SCASE = 'C1'

         ELSEIF (DRNH4HS4.LE.RH) THEN         
            SCASE = 'C2'
            CALL CALCC22p1              
            SCASE = 'C2'
         ENDIF
      ENDIF



      MOLAL(7) = WAER(4)             
      MOLAL(1) = MOLAL(1) + WAER(4)  

      CALL CALCNAP2p1                   
      CALL CALCNH3P2p1
      ENDIF



      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0  &
                                          .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF

      RETURN



      END















      SUBROUTINE ISRP3R2p1 (WI, RHI, TEMPI)
      INCLUDE 'module_isrpia_inc.F'
      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ








      TRYLIQ = .TRUE.             

10    CALL ISOINIT32p1 (WI, RHI, TEMPI) 











      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN  
         FRSO4   = WAER(2) - WAER(1)/2.0D0     
         FRSO4   = MAX(FRSO4, TINY)
         SRI     = GETASR2p1(FRSO4, RHI)          
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
         CALL CALCQ52p1                 
         SCASE = 'Q5'
      ELSE

         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'Q1'
            CALL CALCQ12p1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH4CL) THEN         
            SCASE = 'Q2'
            CALL CALCQ22p1              

         ELSEIF (DRNH4CL.LE.RH  .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'Q3'
            CALL CALCQ32p1              

        ELSEIF (DRNH42S4.LE.RH  .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'Q4'
            CALL CALCQ42p1              
            SCASE = 'Q4'

         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'Q5'
            CALL CALCQ52p1              
            SCASE = 'Q5'
         ENDIF
      ENDIF



      ELSE IF (SULRAT.GE.SULRATW .AND. SODRAT.GE.2.0) THEN                

      IF(METSTBL.EQ.1) THEN
         SCASE = 'R6'
         CALL CALCR62p1                 
         SCASE = 'R6'
      ELSE

         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'R1'
            CALL CALCR12p1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN         
            SCASE = 'R2'
            CALL CALCR22p1              

         ELSEIF (DRNANO3.LE.RH  .AND. RH.LT.DRNACL) THEN         
            SCASE = 'R3'
            CALL CALCR32p1              

         ELSEIF (DRNACL.LE.RH   .AND. RH.LT.DRNH4CL) THEN         
            SCASE = 'R4'
            CALL CALCR42p1              

         ELSEIF (DRNH4CL.LE.RH .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'R5'
            CALL CALCR52p1              
            SCASE = 'R5'

         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'R6'
            CALL CALCR62p1              
            SCASE = 'R6'
         ENDIF
      ENDIF



      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
      DO 100 I=1,NCOMP
         W(I) = WAER(I)
100   CONTINUE

      IF(METSTBL.EQ.1) THEN
         SCASE = 'I6'
         CALL CALCI62p1                 
         SCASE = 'I6'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'I1'
            CALL CALCI12p1              
            SCASE = 'I1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN         
            SCASE = 'I2'
            CALL CALCI22p1              
            SCASE = 'I2'

         ELSEIF (DRNAHSO4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'I3'
            CALL CALCI32p1              
            SCASE = 'I3'

         ELSEIF (DRLC.LE.RH     .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'I4'
            CALL CALCI42p1              
            SCASE = 'I4'

         ELSEIF (DRNH42S4.LE.RH .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'I5'
            CALL CALCI52p1              
            SCASE = 'I5'

         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'I6'
            CALL CALCI62p1              
            SCASE = 'I6'
         ENDIF
      ENDIF

      CALL CALCNHP2p1                
      CALL CALCNH3P2p1



      ELSEIF (SULRAT.LT.1.0) THEN             
      DO 200 I=1,NCOMP
         W(I) = WAER(I)
200   CONTINUE

      IF(METSTBL.EQ.1) THEN
         SCASE = 'J3'
         CALL CALCJ32p1                 
         SCASE = 'J3'
      ELSE

         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'J1'
            CALL CALCJ12p1              
            SCASE = 'J1'

         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN         
            SCASE = 'J2'
            CALL CALCJ22p1              
            SCASE = 'J2'

         ELSEIF (DRNAHSO4.LE.RH) THEN         
            SCASE = 'J3'
            CALL CALCJ32p1              
            SCASE = 'J3'
         ENDIF
      ENDIF

      CALL CALCNHP2p1                
      CALL CALCNH3P2p1

      ENDIF




      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0  &
                            .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF

      RETURN



      END
















      SUBROUTINE ISRP4R2p1 (WI, RHI, TEMPI)
      INCLUDE 'module_isrpia_inc.F'
      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ








      TRYLIQ  = .TRUE.             
      IPROB   = 1            


10    CALL INIT42p1 (WI, RHI, TEMPI) 











      IF (TRYLIQ) THEN                               
         FRSO4   = WAER(2) - WAER(1)/2.0D0 &
                 - WAER(6) - WAER(7)/2.0D0 - WAER(8) 
         FRSO4   = MAX(FRSO4, TINY)
         SRI     = GETASR2p1(FRSO4, RHI)                
         SULRATW = (WAER(1)+FRSO4*SRI+WAER(6) &
                    +WAER(7)+WAER(8))/WAER(2)       
         SULRATW = MIN (SULRATW, 2.0D0)
      ELSE
         SULRATW = 2.0D0                     
      ENDIF
      SO4RAT = (WAER(1)+WAER(3)+WAER(6)+WAER(7)+WAER(8))/WAER(2)
      CRNARAT = (WAER(1)+WAER(6)+WAER(7)+WAER(8))/WAER(2)
      CRRAT  = (WAER(6)+WAER(7)+WAER(8))/WAER(2)





      IF (SULRATW.LE.SO4RAT .AND. CRNARAT.LT.2.0) THEN

       IF(METSTBL.EQ.1) THEN
         SCASE = 'V7'
         CALL CALCV72p1                 
       ELSE

         IF (RH.LT.DRNH4NO3) THEN
            SCASE = 'V1'
            CALL CALCV12p1              

         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH4CL) THEN
            SCASE = 'V2'
            CALL CALCV22p1              

         ELSEIF (DRNH4CL.LE.RH  .AND. RH.LT.DRNH42S4) THEN
            SCASE = 'V3'
            CALL CALCV32p1              

         ELSEIF (DRNH42S4.LE.RH  .AND. RH.LT.DRMGSO4) THEN
            SCASE = 'V4'
            CALL CALCV42p1              

         ELSEIF (DRMGSO4.LE.RH .AND. RH.LT.DRNA2SO4) THEN
            SCASE = 'V5'
            CALL CALCV52p1              

         ELSEIF (DRNA2SO4.LE.RH .AND. RH.LT.DRK2SO4) THEN
            SCASE = 'V6'
            CALL CALCV62p1              

         ELSEIF (DRK2SO4.LE.RH) THEN
            SCASE = 'V7'
            CALL CALCV72p1              
         ENDIF
       ENDIF



      ELSEIF (SO4RAT.GE.SULRATW .AND. CRNARAT.GE.2.0) THEN

       IF (CRRAT.LE.2.0) THEN

        IF(METSTBL.EQ.1) THEN
         SCASE = 'U8'
         CALL CALCU82p1                 
        ELSE

           IF (RH.LT.DRNH4NO3) THEN
             SCASE = 'U1'
             CALL CALCU12p1             

           ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN
             SCASE = 'U2'
             CALL CALCU22p1            

           ELSEIF (DRNANO3.LE.RH  .AND. RH.LT.DRNACL) THEN
             SCASE = 'U3'
             CALL CALCU32p1            

           ELSEIF (DRNACL.LE.RH   .AND. RH.LT.DRNH4Cl) THEN
             SCASE = 'U4'
             CALL CALCU42p1            

           ELSEIF (DRNH4Cl.LE.RH .AND. RH.LT.DRMGSO4) THEN
             SCASE = 'U5'
             CALL CALCU52p1            

           ELSEIF (DRMGSO4.LE.RH .AND. RH.LT.DRNA2SO4) THEN
             SCASE = 'U6'
             CALL CALCU62p1            

           ELSEIF (DRNA2SO4.LE.RH .AND. RH.LT.DRK2SO4) THEN
             SCASE = 'U7'
             CALL CALCU72p1            

           ELSEIF (DRK2SO4.LE.RH) THEN
             SCASE = 'U8'
             CALL CALCU82p1            
           ENDIF
        ENDIF



       ELSEIF (CRRAT.GT.2.0) THEN

        IF(METSTBL.EQ.1) THEN
         SCASE = 'W13'
         CALL CALCW132p1                 
        ELSE

           IF (RH.LT.DRCACL2) THEN
             SCASE = 'W1'
             CALL CALCW12p1             


           ELSEIF (DRCACL2.LE.RH .AND. RH.LT.DRMGCL2) THEN
             SCASE = 'W2'
             CALL CALCW22p1            


           ELSEIF (DRMGCL2.LE.RH  .AND. RH.LT.DRCANO32) THEN
             SCASE = 'W3'
             CALL CALCW32p1            


           ELSEIF (DRCANO32.LE.RH   .AND. RH.LT.DRMGNO32) THEN
             SCASE = 'W4'
             CALL CALCW42p1            


           ELSEIF (DRMGNO32.LE.RH .AND. RH.LT.DRNH4NO3) THEN
             SCASE = 'W5'
             CALL CALCW52p1            


           ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN
             SCASE = 'W6'
             CALL CALCW62p1            

           ELSEIF (DRNANO3.LE.RH .AND. RH.LT.DRNACL) THEN
             SCASE = 'W7'
             CALL CALCW72p1            

           ELSEIF (DRNACL.LE.RH .AND. RH.LT.DRNH4CL) THEN
             SCASE = 'W8'
             CALL CALCW82p1            

           ELSEIF (DRNH4CL.LE.RH .AND. RH.LT.DRKCL) THEN
             SCASE = 'W9'
             CALL CALCW92p1            

           ELSEIF (DRKCL.LE.RH .AND. RH.LT.DRMGSO4) THEN
             SCASE = 'W10'
             CALL CALCW102p1            

           ELSEIF (DRMGSO4.LE.RH .AND. RH.LT.DRKNO3) THEN
             SCASE = 'W11'
             CALL CALCW112p1            

           ELSEIF (DRKNO3.LE.RH .AND. RH.LT.DRK2SO4) THEN
             SCASE = 'W12'
             CALL CALCW122p1            

           ELSEIF (DRK2SO4.LE.RH) THEN
             SCASE = 'W13'
             CALL CALCW132p1            
           ENDIF
         ENDIF

       ENDIF



      ELSEIF (1.0.LE.SO4RAT .AND. SO4RAT.LT.SULRATW) THEN
      DO 800 I=1,NCOMP
         W(I) = WAER(I)
 800  CONTINUE

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

      CALL CALCNHP2p1                
      CALL CALCNH3P2p1               



      ELSEIF (SO4RAT.LT.1.0) THEN
      DO 900 I=1,NCOMP
         W(I) = WAER(I)
 900  CONTINUE

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

      CALL CALCNHP2p1                  
      CALL CALCNH3P2p1                 

      ENDIF




      IF (SULRATW.LE.SO4RAT .AND. SO4RAT.LT.2.0 &
                            .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF

      RETURN



      END
















      SUBROUTINE CALCS22p1
      INCLUDE 'module_isrpia_inc.F'
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

         CALL CALCPH2p1 (2.D0*SO4I - NH4I, HI, OHI)    

         NH3AQ = ZERO                               
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ2p1 (NH4I, OHI, DEL)
            NH4I  = MAX (NH4I-DEL, ZERO) 
            OHI   = MAX (OHI -DEL, TINY)
            NH3AQ = DEL
            HI    = AKW/OHI
         ENDIF

         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)         
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
            CALL CALCACT2p1     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE

20    RETURN



      END





















      SUBROUTINE CALCS12p1
      INCLUDE 'module_isrpia_inc.F'

      CNH42S4 = MIN(WAER(2),0.5d0*WAER(3))  
      GNH3    = ZERO

      W(2)    = CNH42S4
      W(3)    = 2.D0*CNH42S4 + GNH3

      RETURN



      END


















      SUBROUTINE CALCN32p1
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION NH4I, NO3I, NH3AQ, NO3AQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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

         CALL CALCPH2p1 (2.D0*SO4I + NO3I - NH4I, HI, OHI)



         NH3AQ = ZERO
         NO3AQ = ZERO
         GG    = 2.D0*SO4I + NO3I - NH4I
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
         ELSE
            HI    = ZERO
            CALL CALCNIAQ22p1 (GG, NO3I, HI, NO3AQ) 



            CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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
            CALL CALCACT2p1     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    RETURN



      END

















      SUBROUTINE CALCN22p1
      INCLUDE 'module_isrpia_inc.F'

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CHI1   = MIN(WAER(2),0.5d0*WAER(3))     
      CHI2   = MAX(WAER(3) - 2.D0*CHI1, ZERO) 
      CHI3   = MAX(WAER(4) - CHI2, ZERO)      

      PSI2   = CHI2
      PSI3   = CHI3

      CALAOU = .TRUE.              
      PSI1LO = TINY                
      PSI1HI = CHI1                



      X1 = PSI1HI
      Y1 = FUNCN22p1 (X1)
      IF (Y1.LE.EPS) RETURN   
      YHI= Y1                 



      DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, ZERO)
         Y2 = FUNCN22p1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  
         X1 = X2
         Y1 = Y2
10    CONTINUE



      YLO= Y1                      
      IF (ABS(Y2) .LT. EPS) THEN   
         RETURN



      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = CHI4
         YY = FUNCN22p1(P4)
         GOTO 50



      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         P4 = TINY
         YY = FUNCN22p1(P4)
         GOTO 50
      ELSE
         CALL PUSHERR2p1 (0001, 'CALCN2')    
         RETURN
      ENDIF



20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCN22p1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR2p1 (0002, 'CALCN2')    



40    X3 = 0.5*(X1+X2)
      Y3 = FUNCN22p1 (X3)
50    CONTINUE
      RETURN



      END













      DOUBLE PRECISION FUNCTION FUNCN22p1 (P1)
      INCLUDE 'module_isrpia_inc.F'
      DOUBLE PRECISION NH4I, NO3I, NH3AQ, NO3AQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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

         CALL CALCPH2p1 (2.D0*SO4I + NO3I - NH4I, HI, OHI)



         NH3AQ = ZERO
         NO3AQ = ZERO
         GG    = 2.D0*SO4I + NO3I - NH4I
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
         ELSE
            HI    = ZERO
            CALL CALCNIAQ22p1 (GG, NO3I, HI, NO3AQ) 



            CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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



         CALL CALCMR2p1



         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT2p1     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE



20    FUNCN22p1= NH4I*NH4I*SO4I/A4 - ONE 
      RETURN



      END





















      SUBROUTINE CALCN12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCN1A2p1, CALCN22p1



      IF (RH.LT.DRMASAN) THEN    
         SCASE = 'N1 ; SUBCASE 1'  
         CALL CALCN1A2p1              
         SCASE = 'N1 ; SUBCASE 1'
      ELSE
         SCASE = 'N1 ; SUBCASE 2'  
         CALL CALCMDRP2p1 (RH, DRMASAN, DRNH4NO3, CALCN1A2p1, CALCN22p1)
         SCASE = 'N1 ; SUBCASE 2'
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCN1A2p1
      INCLUDE 'module_isrpia_inc.F'








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



      END

















      SUBROUTINE CALCQ52p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE. 
      CALAOU  =.TRUE.



      CALL CALCQ1A2p1

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4

      CALL CALCMR2p1           

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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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
         CALL CALCACT2p1
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



      END

















      SUBROUTINE CALCQ42p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV1
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE. 
      CALAOU  =.TRUE.

      PSCONV1 =.TRUE.
      PSI1O   =-GREAT
      ROOT3   = ZERO



      CALL CALCQ1A2p1

      CHI1   = CNA2SO4      

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4

      CALL CALCMR2p1           

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
         CALL POLY32p1(BB, CC, DD, ROOT3, ISLV)
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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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



      END

















      SUBROUTINE CALCQ32p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL EXNO, EXCL
      EXTERNAL CALCQ1A2p1, CALCQ42p1



      EXNO = WAER(4).GT.TINY   
      EXCL = WAER(5).GT.TINY   

      IF (EXNO .OR. EXCL) THEN             
         SCASE = 'Q3 ; SUBCASE 1'  
         CALL CALCQ3A2p1                                   
         SCASE = 'Q3 ; SUBCASE 1' 

      ELSE                                 
         IF (RH.LT.DRMG3) THEN    
            SCASE = 'Q3 ; SUBCASE 2'  
            CALL CALCQ1A2p1             
            SCASE = 'Q3 ; SUBCASE 2'
         ELSE
            SCASE = 'Q3 ; SUBCASE 3' 
            CALL CALCMDRP2p1 (RH, DRMG3, DRNH42S4, CALCQ1A2p1, CALCQ42p1)
            SCASE = 'Q3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCQ3A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV1, PSCONV6
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE. 
      CALAOU  =.TRUE.

      PSCONV1 =.TRUE.
      PSCONV6 =.TRUE.

      PSI1O   =-GREAT
      PSI6O   =-GREAT

      ROOT1   = ZERO
      ROOT3   = ZERO



      CALL CALCQ1A2p1

      CHI1   = CNA2SO4      
      CHI4   = CNH4CL
      CHI6   = CNH42S4

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4

      CALL CALCMR2p1           

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
         CALL POLY32p1(BB, CC, DD, ROOT3, ISLV)
         IF (ISLV.NE.0) ROOT3 = TINY
         ROOT3 = MIN (ROOT3, WAER(1)/2.0, WAER(2) - ROOT1, CHI1)
         ROOT3 = MAX (ROOT3, ZERO)
         PSI1  = CHI1-ROOT3
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*NH4I*SO4I .GT. A7) THEN
         BB =-(WAER(2)+WAER(3)-ROOT3)
         CC =  WAER(3)*(WAER(2)-ROOT3+0.5D0*WAER(3))
         DD =-((WAER(2)-ROOT3)*WAER(3)**2.D0 + A7)/4.D0
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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



      END

















      SUBROUTINE CALCQ22p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL EXNO, EXCL
      EXTERNAL CALCQ1A2p1, CALCQ3A2p1, CALCQ42p1



      EXNO = WAER(4).GT.TINY   
      EXCL = WAER(5).GT.TINY   

      IF (EXNO) THEN                       
         SCASE = 'Q2 ; SUBCASE 1'  
         CALL CALCQ2A2p1                                   
         SCASE = 'Q2 ; SUBCASE 1' 

      ELSEIF (.NOT.EXNO .AND. EXCL) THEN   
         IF (RH.LT.DRMG2) THEN    
            SCASE = 'Q2 ; SUBCASE 2'  
            CALL CALCQ1A2p1             
            SCASE = 'Q2 ; SUBCASE 2'
         ELSE
            SCASE = 'Q2 ; SUBCASE 3' 
            CALL CALCMDRP2p1 (RH, DRMG2, DRNH4CL, CALCQ1A2p1, CALCQ3A2p1)
            SCASE = 'Q2 ; SUBCASE 3'
         ENDIF

      ELSE                                 
         IF (RH.LT.DRMG3) THEN    
            SCASE = 'Q2 ; SUBCASE 2'  
            CALL CALCQ1A2p1             
            SCASE = 'Q2 ; SUBCASE 2'
         ELSE
            SCASE = 'Q2 ; SUBCASE 4' 
            CALL CALCMDRP2p1 (RH, DRMG3, DRNH42S4, CALCQ1A2p1, CALCQ42p1)
            SCASE = 'Q2 ; SUBCASE 4'
         ENDIF
      ENDIF

      RETURN



      END


















      SUBROUTINE CALCQ2A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV1, PSCONV4, PSCONV6
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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



      CALL CALCQ1A2p1

      CHI1   = CNA2SO4      
      CHI4   = CNH4CL
      CHI6   = CNH42S4

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4

      CALL CALCMR2p1           

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
         CALL POLY32p1(BB, CC, DD, ROOT3, ISLV)
         IF (ISLV.NE.0) ROOT3 = TINY
         ROOT3 = MIN (ROOT3, WAER(1)/2.0, WAER(2) - ROOT1, CHI1)
         ROOT3 = MAX (ROOT3, ZERO)
         PSI1  = CHI1-ROOT3
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*NH4I*SO4I .GT. A7) THEN
         BB =-(WAER(2)+WAER(3)-ROOT2-ROOT3)
         CC = (WAER(3)-ROOT2)*(WAER(2)-ROOT3+0.5D0*(WAER(3)-ROOT2))
         DD =-((WAER(2)-ROOT3)*(WAER(3)-ROOT2)**2.D0 + A7)/4.D0
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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



      END

















      SUBROUTINE CALCQ12p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL EXNO, EXCL
      EXTERNAL CALCQ1A2p1, CALCQ2A2p1, CALCQ3A2p1, CALCQ42p1



      EXNO = WAER(4).GT.TINY   
      EXCL = WAER(5).GT.TINY   

      IF (EXNO .AND. EXCL) THEN           
         IF (RH.LT.DRMG1) THEN    
            SCASE = 'Q1 ; SUBCASE 1'  
            CALL CALCQ1A2p1             
            SCASE = 'Q1 ; SUBCASE 1'
         ELSE
            SCASE = 'Q1 ; SUBCASE 2' 
            CALL CALCMDRP2p1 (RH, DRMG1, DRNH4NO3, CALCQ1A2p1, CALCQ2A2p1)
            SCASE = 'Q1 ; SUBCASE 2'
         ENDIF

      ELSE IF (EXNO .AND. .NOT.EXCL) THEN 
         IF (RH.LT.DRMQ1) THEN    
            SCASE = 'Q1 ; SUBCASE 1'  
            CALL CALCQ1A2p1             
            SCASE = 'Q1 ; SUBCASE 1'
         ELSE
            SCASE = 'Q1 ; SUBCASE 3' 
            CALL CALCMDRP2p1 (RH, DRMQ1, DRNH4NO3, CALCQ1A2p1, CALCQ2A2p1)
            SCASE = 'Q1 ; SUBCASE 3'
         ENDIF

      ELSE IF (.NOT.EXNO .AND. EXCL) THEN 
         IF (RH.LT.DRMG2) THEN    
            SCASE = 'Q1 ; SUBCASE 1'  
            CALL CALCQ1A2p1             
            SCASE = 'Q1 ; SUBCASE 1'
         ELSE
            SCASE = 'Q1 ; SUBCASE 4' 
            CALL CALCMDRP2p1 (RH, DRMG2, DRNH4CL, CALCQ1A2p1, CALCQ3A2p1)
            SCASE = 'Q1 ; SUBCASE 4'
         ENDIF

      ELSE                                
         IF (RH.LT.DRMG3) THEN    
            SCASE = 'Q1 ; SUBCASE 1'  
            CALL CALCQ1A2p1             
            SCASE = 'Q1 ; SUBCASE 1'
         ELSE
            SCASE = 'Q1 ; SUBCASE 5' 
            CALL CALCMDRP2p1 (RH, DRMG3, DRNH42S4, CALCQ1A2p1, CALCQ42p1)
            SCASE = 'Q1 ; SUBCASE 5'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCQ1A2p1
      INCLUDE 'module_isrpia_inc.F'



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



      END
















      SUBROUTINE CALCR62p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCR1A2p1

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      FRST   = .TRUE.
      CALAIN = .TRUE. 
      CALAOU = .TRUE. 



      CALL CALCMR2p1



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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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
         CALL CALCACT2p1
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



      END
















      SUBROUTINE CALCR52p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

      LOGICAL  NEAN, NEAC, NESN, NESC



      CALL CALCR1A2p1                             

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



      CALL CALCMR2p1

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
         CALL POLY32p1(BB, CC, DD, ROOT, ISLV)
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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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



      END

















      SUBROUTINE CALCR42p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCR1A2p1, CALCR52p1



      SCASE = 'R4 ; SUBCASE 2'  
      CALL CALCR1A2p1              
      SCASE = 'R4 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN .OR. EXSN .OR. EXSC) THEN   
         IF (RH.GE.DRMH1) THEN    
            SCASE = 'R4 ; SUBCASE 1' 
            CALL CALCR4A2p1
            SCASE = 'R4 ; SUBCASE 1'
         ENDIF

      ELSE IF (EXAC) THEN                  
         IF (RH.GE.DRMR5) THEN    
            SCASE = 'R4 ; SUBCASE 3'  
            CALL CALCMDRP2p1 (RH, DRMR5, DRNH4CL, CALCR1A2p1, CALCR52p1)
            SCASE = 'R4 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCR4A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV1, PSCONV4
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    = .TRUE.
      CALAIN  = .TRUE. 
      CALAOU  = .TRUE. 
      PSCONV1 = .FALSE.
      PSCONV4 = .FALSE.
      PSIO1   =-GREAT
      PSIO4   =-GREAT



      CALL CALCR1A2p1

      CHI1   = CNA2SO4      
      CHI4   = CNH4CL

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      CALL CALCMR2p1           

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
         CALL POLY32p1(BB, CC, DD, ROOT, ISLV)
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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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



      END

















      SUBROUTINE CALCR32p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCR1A2p1, CALCR4A2p1, CALCR52p1



      SCASE = 'R3 ; SUBCASE 2'  
      CALL CALCR1A2p1              
      SCASE = 'R3 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN .OR. EXSN) THEN                   
         IF (RH.GE.DRMH1) THEN    
            SCASE = 'R3 ; SUBCASE 1' 
            CALL CALCR3A2p1
            SCASE = 'R3 ; SUBCASE 1'
         ENDIF

      ELSE IF (.NOT.EXAN .AND. .NOT.EXSN) THEN   
         IF      (     EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN    
               SCASE = 'R3 ; SUBCASE 3'  
               CALL CALCMDRP2p1 (RH, DRMR4, DRNACL, CALCR1A2p1, CALCR4A2p1)
               SCASE = 'R3 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN    
               SCASE = 'R3 ; SUBCASE 4'  
               CALL CALCMDRP2p1 (RH, DRMR2, DRNACL, CALCR1A2p1, CALCR4A2p1)
               SCASE = 'R3 ; SUBCASE 4'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN    
               SCASE = 'R3 ; SUBCASE 5'  
               CALL CALCMDRP2p1 (RH, DRMR5, DRNACL, CALCR1A2p1, CALCR52p1)
               SCASE = 'R3 ; SUBCASE 5'
            ENDIF
         ENDIF

      ENDIF

      RETURN



      END


















      SUBROUTINE CALCR3A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV1, PSCONV3, PSCONV4
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,    &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,    &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,  &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,      &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,&
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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



      CALL CALCR1A2p1

      CHI1   = CNA2SO4      
      CHI4   = CNH4CL
      CHI3   = CNACL

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      CALL CALCMR2p1           

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

      CALL CALCACT2p1          



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
         ROOT2 = MIN(MAX(ZERO, ROOT2), MAX(WAER(5)-ROOT3,ZERO),  &
                     CHI4, WAER(3))
         PSI4  = CHI4 - ROOT2
      ENDIF
      PSCONV4 = ABS(PSI4-PSI4O) .LE. EPS*PSI4O
      PSI4O   = PSI4



      IF (NAI*NAI*SO4I .GT. A5) THEN
         BB =-(CHI1 + WAER(1) - ROOT3)
         CC = 0.25D0*(WAER(1) - ROOT3)*(4.D0*CHI1+WAER(1)-ROOT3)
         DD =-0.25D0*(CHI1*(WAER(1)-ROOT3)**2.D0 - A5) 
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MIN (MAX(ROOT1,ZERO), MAX(WAER(1)-ROOT3,ZERO),  &
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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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
         CALL CALCR1A2p1
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



      END

















      SUBROUTINE CALCR22p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCR1A2p1, CALCR3A2p1, CALCR4A2p1, CALCR52p1



      SCASE = 'R2 ; SUBCASE 2'  
      CALL CALCR1A2p1              
      SCASE = 'R2 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN) THEN                             
         IF (RH.GE.DRMH1) THEN    
            SCASE = 'R2 ; SUBCASE 1' 
            CALL CALCR2A2p1
            SCASE = 'R2 ; SUBCASE 1'
         ENDIF

      ELSE IF (.NOT.EXAN) THEN                   
         IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMH2) THEN    
               SCASE = 'R2 ; SUBCASE 3'  
               CALL CALCMDRP2p1 (RH, DRMH2, DRNANO3, CALCR1A2p1, CALCR3A2p1)
               SCASE = 'R2 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR1) THEN    
               SCASE = 'R2 ; SUBCASE 4'  
               CALL CALCMDRP2p1 (RH, DRMR1, DRNANO3, CALCR1A2p1, CALCR3A2p1)
               SCASE = 'R2 ; SUBCASE 4'
            ENDIF

         ELSE IF (.NOT.EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN    
               SCASE = 'R2 ; SUBCASE 5'  
               CALL CALCMDRP2p1 (RH, DRMR2, DRNACL, CALCR1A2p1, CALCR4A2p1)
               SCASE = 'R2 ; SUBCASE 5'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR3) THEN    
               SCASE = 'R2 ; SUBCASE 6'  
               CALL CALCMDRP2p1 (RH, DRMR3, DRNANO3, CALCR1A2p1, CALCR3A2p1)
               SCASE = 'R2 ; SUBCASE 6'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN    
               SCASE = 'R2 ; SUBCASE 7'  
               CALL CALCMDRP2p1 (RH, DRMR4, DRNACL, CALCR1A2p1, CALCR4A2p1)
               SCASE = 'R2 ; SUBCASE 7'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN    
               SCASE = 'R2 ; SUBCASE 8'  
               CALL CALCMDRP2p1 (RH, DRMR5, DRNH4CL, CALCR1A2p1, CALCR52p1)
               SCASE = 'R2 ; SUBCASE 8'
            ENDIF

         ELSE IF (     EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR6) THEN    
               SCASE = 'R2 ; SUBCASE 9'  
               CALL CALCMDRP2p1 (RH, DRMR6, DRNANO3, CALCR1A2p1, CALCR3A2p1)
               SCASE = 'R2 ; SUBCASE 9'
            ENDIF
         ENDIF

      ENDIF

      RETURN



      END


















      SUBROUTINE CALCR2A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV1, PSCONV2, PSCONV3, PSCONV4
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



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



      CALL CALCR1A2p1

      CHI1   = CNA2SO4      
      CHI2   = CNANO3
      CHI3   = CNACL
      CHI4   = CNH4CL

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3

      CALL CALCMR2p1           

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

      CALL CALCACT2p1          



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
         CC = WAER(1)*(2.D0*ROOT3 + 2.D0*ROOT4 - 4.D0*WAER(2) - ONE) &
             -(ROOT3 + ROOT4)**2.0 + 4.D0*WAER(2)*(ROOT3 + ROOT4)
         CC =-0.25*CC
         DD = WAER(1)*WAER(2)*(ONE - 2.D0*ROOT3 - 2.D0*ROOT4) + &
              WAER(2)*(ROOT3 + ROOT4)**2.0 - A5
         DD =-0.25*DD
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
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
         CALL CALCAMAQ22p1 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ22p1 (GGCL, CLI, HI, CLAQ) 
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ22p1 (GGNO3, NO3I, HI, NO3AQ)              
         ENDIF



         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
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



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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
         CALL CALCR1A2p1
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



      END

















      SUBROUTINE CALCR12p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCR1A2p1, CALCR2A2p1, CALCR3A2p1, CALCR4A2p1, CALCR52p1



      SCASE = 'R1 ; SUBCASE 1'  
      CALL CALCR1A2p1              
      SCASE = 'R1 ; SUBCASE 1'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN.AND.EXAC.AND.EXSC.AND.EXSN) THEN  
         IF (RH.GE.DRMH1) THEN    
            SCASE = 'R1 ; SUBCASE 2'  
            CALL CALCMDRP2p1 (RH, DRMH1, DRNH4NO3, CALCR1A2p1, CALCR2A2p1)
            SCASE = 'R1 ; SUBCASE 2'
         ENDIF

      ELSE IF (.NOT.EXAN) THEN                   
         IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMH2) THEN    
               SCASE = 'R1 ; SUBCASE 3'  
               CALL CALCMDRP2p1 (RH, DRMH2, DRNANO3, CALCR1A2p1, CALCR3A2p1)
               SCASE = 'R1 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR1) THEN    
               SCASE = 'R1 ; SUBCASE 4'  
               CALL CALCMDRP2p1 (RH, DRMR1, DRNANO3, CALCR1A2p1, CALCR3A2p1)
               SCASE = 'R1 ; SUBCASE 4'
            ENDIF

         ELSE IF (.NOT.EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN    
               SCASE = 'R1 ; SUBCASE 5'  
               CALL CALCMDRP2p1 (RH, DRMR2, DRNACL, CALCR1A2p1, CALCR3A2p1) 
               SCASE = 'R1 ; SUBCASE 5'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR3) THEN    
               SCASE = 'R1 ; SUBCASE 6'  
               CALL CALCMDRP2p1 (RH, DRMR3, DRNANO3, CALCR1A2p1, CALCR3A2p1)
               SCASE = 'R1 ; SUBCASE 6'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN    
               SCASE = 'R1 ; SUBCASE 7'  
               CALL CALCMDRP2p1 (RH, DRMR4, DRNACL, CALCR1A2p1, CALCR3A2p1) 
               SCASE = 'R1 ; SUBCASE 7'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN    
               SCASE = 'R1 ; SUBCASE 8'  
               CALL CALCMDRP2p1 (RH, DRMR5, DRNH4CL, CALCR1A2p1, CALCR3A2p1) 
               SCASE = 'R1 ; SUBCASE 8'
            ENDIF

         ELSE IF (     EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR6) THEN    
               SCASE = 'R1 ; SUBCASE 9'  
               CALL CALCMDRP2p1 (RH, DRMR6, DRNANO3, CALCR1A2p1, CALCR3A2p1)
               SCASE = 'R1 ; SUBCASE 9'
            ENDIF
         ENDIF

      ELSE IF (.NOT.EXAC) THEN                   
         IF      (     EXAN .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR7) THEN    
               SCASE = 'R1 ; SUBCASE 10'  
               CALL CALCMDRP2p1 (RH, DRMR7, DRNH4NO3, CALCR1A2p1, CALCR2A2p1)
               SCASE = 'R1 ; SUBCASE 10'
            ENDIF

         ELSE IF (     EXAN .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR8) THEN    
               SCASE = 'R1 ; SUBCASE 11'  
               CALL CALCMDRP2p1 (RH, DRMR8, DRNH4NO3, CALCR1A2p1, CALCR2A2p1)
               SCASE = 'R1 ; SUBCASE 11'
            ENDIF

         ELSE IF (     EXAN .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR9) THEN    
               SCASE = 'R1 ; SUBCASE 12'  
               CALL CALCMDRP2p1 (RH, DRMR9, DRNH4NO3, CALCR1A2p1, CALCR2A2p1)
               SCASE = 'R1 ; SUBCASE 12'
            ENDIF

         ELSE IF (     EXAN .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR10) THEN    
               SCASE = 'R1 ; SUBCASE 13'  
               CALL CALCMDRP2p1 (RH, DRMR10, DRNH4NO3, CALCR1A2p1, CALCR2A2p1)
               SCASE = 'R1 ; SUBCASE 13'
            ENDIF
         ENDIF

      ELSE IF (.NOT.EXSN) THEN                  
         IF      (     EXAN .AND.      EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR11) THEN    
               SCASE = 'R1 ; SUBCASE 14'  
               CALL CALCMDRP2p1 (RH, DRMR11, DRNH4NO3, CALCR1A2p1, CALCR2A2p1)
               SCASE = 'R1 ; SUBCASE 14'
            ENDIF

         ELSE IF (     EXAN .AND.      EXAC .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR12) THEN    
               SCASE = 'R1 ; SUBCASE 15'  
               CALL CALCMDRP2p1 (RH, DRMR12, DRNH4NO3, CALCR1A2p1, CALCR2A2p1)
               SCASE = 'R1 ; SUBCASE 15'
            ENDIF
         ENDIF

      ELSE IF (.NOT.EXSC) THEN                  
         IF      (     EXAN .AND.      EXAC .AND.      EXSN) THEN
            IF (RH.GE.DRMR13) THEN    
               SCASE = 'R1 ; SUBCASE 16'  
               CALL CALCMDRP2p1 (RH, DRMR13, DRNH4NO3, CALCR1A2p1, CALCR2A2p1)
               SCASE = 'R1 ; SUBCASE 16'
            ENDIF
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCR1A2p1
      INCLUDE 'module_isrpia_inc.F'



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



      END















      SUBROUTINE CALCV72p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,       &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,       &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,     &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,         &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,   &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.



      CALL CALCV1A2p1

      CHI9   = CCASO4

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      AKW = XKW*RH*WATER*WATER               



      NAI    = WAER(1)
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      KI     = WAER(7)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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




      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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
      CNA2SO4 = ZERO
      CMGSO4  = ZERO
      CK2SO4  = ZERO
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END



















      SUBROUTINE CALCV62p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.TRUE.
      PSI70   =-GREAT                                 
      ROOT7   = ZERO



      CALL CALCV1A2p1

      CHI9   = CCASO4
      CHI7   = CK2SO4       

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7))
         CC = WAER(7)*(WAER(2)-WAER(6)) + 0.25D0*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*WAER(2) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MIN (ROOT7,WAER(7)/2.0,MAX(WAER(2)-WAER(6),ZERO),CHI7)
         ROOT7 = MAX (ROOT7, ZERO)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT7, ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      NAI    = WAER(1)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7) GOTO 20
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
      CNA2SO4 = ZERO
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END


















      SUBROUTINE CALCV52p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.TRUE.
      PSCONV1 =.TRUE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO



      CALL CALCV1A2p1

      CHI9   = CCASO4
      CHI7   = CK2SO4       
      CHI1   = CNA2SO4

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
         CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX(WAER(2)-WAER(6) - ROOT1, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
         CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MAX (ROOT1, ZERO)
         ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
                 MAX ((WAER(2)-WAER(6)) - ROOT7, ZERO), CHI1)
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
      SO4I   = MAX ((WAER(2)-WAER(6)) - ROOT7 - ROOT1, ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1) GOTO 20
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
      CNA2SO4 = CHI1 - PSI1
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END



















      SUBROUTINE CALCV42p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.TRUE.
      PSCONV1 =.TRUE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO



      CALL CALCV1A2p1

      CHI9   = CCASO4
      CHI7   = CK2SO4       
      CHI1   = CNA2SO4
      CHI8   = CMGSO4

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
         CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX((WAER(2)-WAER(6)) - ROOT1, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
         CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MAX (ROOT1, ZERO)
         ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
                 MAX ((WAER(2)-WAER(6)) - ROOT7, ZERO), CHI1)
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
      SO4I   = MAX ((WAER(2)-WAER(6)) - ROOT7 - ROOT1, ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1) GOTO 20
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
      CNA2SO4 = CHI1 - PSI1
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END


















      SUBROUTINE CALCV32p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL EXNO, EXCL
      EXTERNAL CALCV1A2p1, CALCV42p1



      EXNO = WAER(4).GT.TINY
      EXCL = WAER(5).GT.TINY

      IF (EXNO .OR. EXCL) THEN             
         SCASE = 'V3 ; SUBCASE 1'
         CALL CALCV3A2p1
         SCASE = 'V3 ; SUBCASE 1'

      ELSE                                 
         IF (RH.LT.DRMO3) THEN
            SCASE = 'V3 ; SUBCASE 2'
            CALL CALCV1A2p1             
            SCASE = 'V3 ; SUBCASE 2'
         ELSE
            SCASE = 'V3 ; SUBCASE 3' 
            CALL CALCMDRPII2p1 (RH, DRMO3, DRNH42S4, CALCV1A2p1, CALCV42p1)
            SCASE = 'V3 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCV3A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1, PSCONV6
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.TRUE.
      PSCONV1 =.TRUE.
      PSCONV6 =.TRUE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT
      PSI60   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO
      ROOT6   = ZERO



      CALL CALCV1A2p1

      CHI9   = CCASO4
      CHI7   = CK2SO4       
      CHI1   = CNA2SO4
      CHI8   = CMGSO4
      CHI6   = CNH42S4

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      A6  = XK7 *(WATER/GAMA(4))**3.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1 - ROOT6)
         CC = WAER(7)*((WAER(2) - WAER(6)) - ROOT1 - ROOT6) + &
              0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6))-ROOT1-ROOT6)-A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX (WAER(2)-WAER(6)-ROOT1-ROOT6, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7 - ROOT6)
         CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7 - ROOT6) + &
              0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6))-ROOT7-ROOT6)-A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MAX (ROOT1, ZERO)
         ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
                      MAX (WAER(2)-WAER(6)-ROOT7-ROOT6, ZERO), CHI1)
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*NH4I*SO4I .GT. A6) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(3) - ROOT7 - ROOT1)
         CC = WAER(3)*((WAER(2)-WAER(6)) - ROOT7 - ROOT1) + &
              0.25*WAER(3)*WAER(3)
         DD =-0.25*(WAER(3)*WAER(3)*((WAER(2)-WAER(6))-ROOT7-ROOT1)-A6)
         CALL POLY32p1(BB, CC, DD, ROOT6, ISLV)
         IF (ISLV.NE.0) ROOT6 = TINY
         ROOT6 = MAX (ROOT6, ZERO)
         ROOT6 = MIN (ROOT6, WAER(3)/2.0, &
                      MAX (WAER(2)-WAER(6)-ROOT7-ROOT1, ZERO), CHI6)
         PSI6  = CHI6-ROOT6
      ENDIF
      PSCONV6 = ABS(PSI6-PSI60) .LE. EPS*PSI60
      PSI60   = PSI6


      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1 - ROOT6, ZERO)
      NH4I   = MAX (WAER(3) - 2.D0*ROOT6, ZERO)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV6) GOTO 20
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
      CNA2SO4 = CHI1 - PSI1
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END


















      SUBROUTINE CALCV22p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL EXNO, EXCL
      EXTERNAL CALCV1A2p1, CALCV3A2p1, CALCV42p1



      EXNO = WAER(4).GT.TINY
      EXCL = WAER(5).GT.TINY

      IF (EXNO) THEN                       
         SCASE = 'V2 ; SUBCASE 1'
         CALL CALCV2A2p1
         SCASE = 'V2 ; SUBCASE 1'

      ELSEIF (.NOT.EXNO .AND. EXCL) THEN   
         IF (RH.LT.DRMO2) THEN
            SCASE = 'V2 ; SUBCASE 2'
            CALL CALCV1A2p1             
            SCASE = 'V2 ; SUBCASE 2'
         ELSE
            SCASE = 'V2 ; SUBCASE 3' 
            CALL CALCMDRPII2p1 (RH, DRMO2, DRNH4CL, CALCV1A2p1, CALCV3A2p1)
            SCASE = 'V2 ; SUBCASE 3'
         ENDIF

      ELSE                                 
         IF (RH.LT.DRMO3) THEN
            SCASE = 'V2 ; SUBCASE 2'
            CALL CALCV1A2p1             
            SCASE = 'V2 ; SUBCASE 2'
         ELSE
            SCASE = 'V2 ; SUBCASE 4' 
            CALL CALCMDRPII2p1 (RH, DRMO3, DRNH42S4, CALCV1A2p1, CALCV42p1)
            SCASE = 'V2 ; SUBCASE 4'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCV2A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1, PSCONV6, PSCONV4
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.TRUE.
      PSCONV1 =.TRUE.
      PSCONV6 =.TRUE.
      PSCONV4 =.TRUE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT
      PSI60   =-GREAT
      PSI40   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO
      ROOT6   = ZERO
      ROOT4   = ZERO



      CALL CALCV1A2p1

      CHI9   = CCASO4
      CHI7   = CK2SO4       
      CHI1   = CNA2SO4
      CHI8   = CMGSO4
      CHI6   = CNH42S4
      CHI4   = CNH4CL

      PSI1   = CNA2SO4      
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      A6  = XK7 *(WATER/GAMA(4))**3.0        
      A14 = XK14*(WATER/GAMA(6))**2.         
      AKW = XKW*RH*WATER*WATER               



      IF (NH4I*CLI .GT. A14) THEN
         BB    =-(WAER(3) + WAER(5) - 2.D0*ROOT6)
         CC    = WAER(5)*(WAER(3) - 2.D0*ROOT6) - A14
         DD    = BB*BB - 4.D0*CC
         IF (DD.LT.ZERO) THEN
            ROOT4 = ZERO
         ELSE
            DD    = SQRT(DD)
            ROOT4A= 0.5D0*(-BB+DD)
            ROOT4B= 0.5D0*(-BB-DD)
            IF (ZERO.LE.ROOT4A) THEN
               ROOT4 = ROOT4A
            ELSE
               ROOT4 = ROOT4B
            ENDIF
            ROOT4 = MAX(ROOT4, ZERO)
            ROOT4 = MIN(ROOT4, WAER(5), &
                        MAX (WAER(3) - 2.D0*ROOT6, ZERO), CHI4)
            PSI4  = CHI4 - ROOT4
         ENDIF
      ENDIF
      PSCONV4 = ABS(PSI4-PSI40) .LE. EPS*PSI40
      PSI40   = PSI4



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2) - WAER(6)) + WAER(7) - ROOT1 - ROOT6)
         CC = WAER(7)*((WAER(2) - WAER(6)) - ROOT1 - ROOT6) &
              + 0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6))-ROOT1-ROOT6)-A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX (WAER(2)-WAER(6)-ROOT1-ROOT6, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-((WAER(2) - WAER(6)) + WAER(1) - ROOT7 - ROOT6)
         CC = WAER(1)*((WAER(2) - WAER(6)) - ROOT7 - ROOT6) + &
              0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6))-ROOT7-ROOT6)-A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MAX (ROOT1, ZERO)
         ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
                      MAX (WAER(2)-WAER(6)-ROOT7-ROOT6, ZERO), CHI1)
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*NH4I*SO4I .GT. A6) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(3) - ROOT7 - ROOT1 - ROOT4)
         CC = WAER(3)*((WAER(2)-WAER(6)) - ROOT7 - ROOT1) + 0.25* &
            (WAER(3)-ROOT4)**2.0 + ROOT4*(ROOT1+ROOT7-(WAER(2)-WAER(6)))
         DD =-0.25*((WAER(3)-ROOT4)**2.0 * &
                    ((WAER(2)-WAER(6))-ROOT7-ROOT1) - A6)
         CALL POLY32p1(BB, CC, DD, ROOT6, ISLV)
         IF (ISLV.NE.0) ROOT6 = TINY
         ROOT6 = MAX (ROOT6, ZERO)
         ROOT6 = MIN (ROOT6, WAER(3)/2.0, &
                      MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1, ZERO), CHI6)
         PSI6  = CHI6-ROOT6
      ENDIF
      PSCONV6 = ABS(PSI6-PSI60) .LE. EPS*PSI60
      PSI60   = PSI6



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1 - ROOT6, ZERO)
      NH4I   = MAX (WAER(3) - 2.D0*ROOT6, ZERO)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV6 .AND. PSCONV4) GOTO 20
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
      CNA2SO4 = CHI1 - PSI1
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END


















      SUBROUTINE CALCV12p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL EXNO, EXCL
      EXTERNAL CALCV1A2p1, CALCV2A2p1, CALCV3A2p1, CALCV42p1



      EXNO = WAER(4).GT.TINY
      EXCL = WAER(5).GT.TINY

      IF (EXNO .AND. EXCL) THEN           
         IF (RH.LT.DRMO1) THEN
            SCASE = 'V1 ; SUBCASE 1'
            CALL CALCV1A2p1             
            SCASE = 'V1 ; SUBCASE 1'
         ELSE
            SCASE = 'V1 ; SUBCASE 2' 
            CALL CALCMDRPII2p1 (RH, DRMO1, DRNH4NO3, CALCV1A2p1, CALCV2A2p1)
            SCASE = 'V1 ; SUBCASE 2'
         ENDIF

      ELSE IF (EXNO .AND. .NOT.EXCL) THEN 
         IF (RH.LT.DRMV1) THEN
            SCASE = 'V1 ; SUBCASE 1'
            CALL CALCV1A2p1             
            SCASE = 'V1 ; SUBCASE 1'
         ELSE
            SCASE = 'V1 ; SUBCASE 3' 
            CALL CALCMDRPII2p1 (RH, DRMV1, DRNH4NO3, CALCV1A2p1, CALCV2A2p1)
            SCASE = 'V1 ; SUBCASE 3'
         ENDIF

      ELSE IF (.NOT.EXNO .AND. EXCL) THEN 
         IF (RH.LT.DRMO2) THEN
            SCASE = 'V1 ; SUBCASE 1'
            CALL CALCV1A2p1             
            SCASE = 'V1 ; SUBCASE 1'
         ELSE
            SCASE = 'V1 ; SUBCASE 4' 
            CALL CALCMDRPII2p1 (RH, DRMO2, DRNH4CL, CALCV1A2p1, CALCV3A2p1)
            SCASE = 'V1 ; SUBCASE 4'
         ENDIF

      ELSE                                
         IF (RH.LT.DRMO3) THEN
            SCASE = 'V1 ; SUBCASE 1'
            CALL CALCV1A2p1             
            SCASE = 'V1 ; SUBCASE 1'
         ELSE
            SCASE = 'V1 ; SUBCASE 5' 
            CALL CALCMDRPII2p1 (RH, DRMO3, DRNH42S4, CALCV1A2p1, CALCV42p1)
            SCASE = 'V1 ; SUBCASE 5'
         ENDIF
      ENDIF

      RETURN















      END


















      SUBROUTINE CALCV1A2p1
      INCLUDE 'module_isrpia_inc.F'



      CCASO4  = MIN (WAER(6), WAER(2))                     
      SO4FR   = MAX (WAER(2) - CCASO4, ZERO)
      CAFR    = MAX (WAER(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*WAER(7), SO4FR)                 
      FRK     = MAX (WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
      CNA2SO4 = MIN (0.5D0*WAER(1), SO4FR)                 
      NAFR    = MAX (WAER(1) - 2.D0*CNA2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CNA2SO4, ZERO)
      CMGSO4  = MIN (WAER(8), SO4FR)                       
      FRMG    = MAX(WAER(8) - CMGSO4, ZERO)
      SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
      CNH42S4 = MAX (MIN (SO4FR , 0.5d0*WAER(3)) , TINY)
      FRNH3   = MAX (WAER(3) - 2.D0*CNH42S4, ZERO)

      CNH4NO3 = MIN (FRNH3, WAER(4))

      FRNH3   = MAX (FRNH3 - CNH4NO3, ZERO)

      CNH4CL  = MIN (FRNH3, WAER(5))

      FRNH3   = MAX (FRNH3 - CNH4CL, ZERO)



      WATER   = ZERO

      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      RETURN



      END
















      SUBROUTINE CALCU82p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      CALL CALCU1A2p1

      CHI9   = CCASO4        

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      FRST   = .TRUE.
      CALAIN = .TRUE.
      CALAOU = .TRUE.



      CALL CALCMR2p1



      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      AKW = XKW*RH*WATER*WATER                        

      NAI    = WAER(1)
      SO4I   = MAX(WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      KI     = WAER(7)
      MGI    = WAER(8)




      GG  = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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
      IF (HI.LE.TINY) HI = SQRT(AKW)




      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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
      CMGSO4   = ZERO
      CK2SO4   = ZERO
      CCASO4   = MIN (WAER(6), WAER(2))

      RETURN



      END



















      SUBROUTINE CALCU72p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.TRUE.
      PSI70   =-GREAT                                 
      ROOT7   = ZERO



      CALL CALCU1A2p1

      CHI7   = CK2SO4       
      CHI9   = CCASO4

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4


      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

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
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI

      CALL CALCACT2p1          



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7))
         CC = WAER(7)*(WAER(2)-WAER(6)) + 0.25D0*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*(WAER(2)-WAER(6)) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7,WAER(7)/2.0,MAX(WAER(2)-WAER(6),ZERO),CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      SO4I   = MAX (WAER(2) - WAER(6) - ROOT7, ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      NAI    = WAER(1)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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





      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7) GOTO 20
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
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END


















      SUBROUTINE CALCU62p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.TRUE.
      PSCONV1 =.TRUE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO



      CALL CALCU1A2p1

      CHI1   = CNA2SO4            
      CHI7   = CK2SO4
      CHI9   = CCASO4

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

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
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI

      CALL CALCACT2p1          



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
         CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX((WAER(2)-WAER(6)) - ROOT1,ZERO), CHI7)
         PSI7  = CHI7-ROOT7

      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
         CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MAX (ROOT1, ZERO)
         ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
                      MAX((WAER(2)-WAER(6)) - ROOT7, ZERO) ,CHI1)
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
      SO4I   = MAX (WAER(2) - WAER(6) - ROOT7 - ROOT1, ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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





      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1) GOTO 20
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
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END


















      SUBROUTINE CALCU52p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.TRUE.
      PSCONV1 =.TRUE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO



      CALL CALCU1A2p1

      CHI1   = CNA2SO4            
      CHI7   = CK2SO4
      CHI8   = CMGSO4
      CHI9   = CCASO4

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

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
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI

      CALL CALCACT2p1          



      DO 10 I=1,NSWEEP
      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
         CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX(WAER(2)-WAER(6)-ROOT1, ZERO),CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
         CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MAX (ROOT1, ZERO)
         ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
                      MAX(WAER(2)-WAER(6)-ROOT7, ZERO),CHI1)
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1, ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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





      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1) GOTO 20
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
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END


















      SUBROUTINE CALCU42p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCU1A2p1, CALCU52p1



      SCASE = 'U4 ; SUBCASE 2'
      CALL CALCU1A2p1              
      SCASE = 'U4 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN .OR. EXSN .OR. EXSC) THEN   
         IF (RH.GE.DRMM1) THEN
            SCASE = 'U4 ; SUBCASE 1'
            CALL CALCU4A2p1
            SCASE = 'U4 ; SUBCASE 1'
         ENDIF

      ELSE IF (EXAC) THEN                  
         IF (RH.GE.DRMR5) THEN
            SCASE = 'U4 ; SUBCASE 3'
            CALL CALCMDRPII2p1 (RH, DRMR5, DRNH4CL, CALCU1A2p1, CALCU52p1)
            SCASE = 'U4 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCU4A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1, PSCONV4
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.FALSE.
      PSCONV1 =.FALSE.
      PSCONV4 =.FALSE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT
      PSI40   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO
      ROOT4   = ZERO



      CALL CALCU1A2p1

      CHI1   = CNA2SO4            
      CHI4   = CNH4CL
      CHI7   = CK2SO4
      CHI8   = CMGSO4
      CHI9   = CCASO4

      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

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
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI

      CALL CALCACT2p1          



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      A14 = XK14*(WATER/GAMA(6))**2.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
         CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX(WAER(2)-WAER(6)-ROOT1, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
         CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
         DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MAX (ROOT1, ZERO)
         ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
                      MAX (WAER(2)-WAER(6)-ROOT7, ZERO), CHI1)
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*CLI .GT. A14) THEN
         BB   =-(NH4I + CLI)
         CC   =-A14 + NH4I*CLI
         DD   = BB*BB - 4.D0*CC
         ROOT4 = 0.5D0*(-BB-SQRT(DD))
         IF (ROOT4.GT.TINY) THEN
            ROOT4    = MIN(MAX (ROOT4, ZERO), CHI4)
            PSI4    = CHI4 - ROOT4
         ENDIF
      ENDIF
      PSCONV4 = ABS(PSI4-PSI40) .LE. EPS*PSI40
      PSI40   = PSI4



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
      SO4I   = MAX (WAER(2) - WAER(6) - ROOT7 - ROOT1, ZERO)
      NH4I   = MAX (WAER(3) - ROOT4, ZERO)
      NO3I   = WAER(4)
      CLI    = MAX (WAER(5) - ROOT4, ZERO)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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





      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV4) GOTO 20
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
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))

      RETURN



      END


















      SUBROUTINE CALCU32p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCU1A2p1, CALCU4A2p1, CALCU52p1



      SCASE = 'U3 ; SUBCASE 2'
      CALL CALCU1A2p1              
      SCASE = 'U3 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN .OR. EXSN) THEN                   
         IF (RH.GE.DRMM1) THEN
            SCASE = 'U3 ; SUBCASE 1'
            CALL CALCU3A2p1
            SCASE = 'U3 ; SUBCASE 1'
         ENDIF

      ELSE IF (.NOT.EXAN .AND. .NOT.EXSN) THEN   
         IF      (     EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN
               SCASE = 'U3 ; SUBCASE 3'
               CALL CALCMDRPII2p1 (RH, DRMR4, DRNACL, CALCU1A2p1, CALCU4A2p1)
               SCASE = 'U3 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN
               SCASE = 'U3 ; SUBCASE 4'
               CALL CALCMDRPII2p1 (RH, DRMR2, DRNACL, CALCU1A2p1, CALCU4A2p1)
               SCASE = 'U3 ; SUBCASE 4'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN
               SCASE = 'U3 ; SUBCASE 5'
               CALL CALCMDRPII2p1 (RH, DRMR5, DRNACL, CALCU1A2p1, CALCU52p1)
               SCASE = 'U3 ; SUBCASE 5'
            ENDIF
         ENDIF

      ENDIF

      RETURN



      END



















      SUBROUTINE CALCU3A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1, PSCONV4, PSCONV3
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.FALSE.
      PSCONV1 =.FALSE.
      PSCONV4 =.FALSE.
      PSCONV3 =.FALSE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT
      PSI40   =-GREAT
      PSI30   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO
      ROOT4   = ZERO
      ROOT3   = ZERO



      CALL CALCU1A2p1

      CHI1   = CNA2SO4            
      CHI3   = CNACL
      CHI4   = CNH4CL
      CHI7   = CK2SO4
      CHI8   = CMGSO4
      CHI9   = CCASO4

      PSI1   = CNA2SO4      
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

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
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI

      CALL CALCACT2p1          



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      A14 = XK14*(WATER/GAMA(6))**2.0        
      AKW = XKW*RH*WATER*WATER               
      A8  = XK8 *(WATER/GAMA(1))**2.0        



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
         CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX(WAER(2)-WAER(6)-ROOT1, ZERO),CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-(((WAER(2)-WAER(6))-ROOT7)*(WAER(1) - ROOT3))
         CC = ((WAER(2) - WAER(6)) - ROOT7)*(WAER(1) - ROOT3) + &
                0.25D0*(WAER(1) - ROOT3)**2.
         DD =-0.25D0*(((WAER(2) - WAER(6)) - ROOT7)* &
                        (WAER(1) - ROOT3)**2.D0 - A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MIN (MAX(ROOT1, ZERO), MAX(WAER(1) - ROOT3, ZERO), &
                      CHI1, MAX(WAER(2)-WAER(6), ZERO))
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*CLI .GT. A14) THEN
         BB    =-(WAER(3) + WAER(5) - ROOT4)
         CC    =-A14 + NH4I*(WAER(5) - ROOT4)
         DD    = MAX(BB*BB - 4.D0*CC, ZERO)
         ROOT4A= 0.5D0*(-BB+SQRT(DD))
         ROOT4B= 0.5D0*(-BB-SQRT(DD))
         IF (ZERO.LE.ROOT4A) THEN
            ROOT4 = ROOT4A
         ELSE
            ROOT4 = ROOT4B
         ENDIF
         ROOT4 = MIN(MAX(ZERO, ROOT4), MAX(WAER(5)-ROOT3,ZERO), &
                     CHI4, WAER(3))
         PSI4  = CHI4 - ROOT4
      ENDIF
      PSCONV4 = ABS(PSI4-PSI40) .LE. EPS*PSI40
      PSI40   = PSI4



      IF (NAI*CLI .GT. A8) THEN
         BB    =-((CHI1-2.D0*ROOT1) + (WAER(5) - ROOT4))
         CC    = (CHI1-2.D0*ROOT1)*(WAER(5) - ROOT4) - A8
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
      PSCONV3 = ABS(PSI3-PSI30) .LE. EPS*PSI30
      PSI30   = PSI3



      KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.D0*ROOT1 - ROOT3, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1, ZERO)
      NH4I   = MAX (WAER(3) - ROOT4, ZERO)
      NO3I   = WAER(4)
      CLI    = MAX (WAER(5) - ROOT4 - ROOT3, ZERO)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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





      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV4 .AND. PSCONV3) GOTO 20
      ENDIF
10    CONTINUE




20    IF (CLI.LE.TINY .AND. WAER(5).GT.TINY) THEN 
         DO 30 I=1,NIONS
            MOLAL(I) = ZERO
30       CONTINUE
         DO 40 I=1,NGASAQ
            GASAQ(I) = ZERO
40       CONTINUE
         CALL CALCU1A2p1
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
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))
      ENDIF

      RETURN



      END


















      SUBROUTINE CALCU22p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCU1A2p1, CALCU3A2p1, CALCU4A2p1, CALCU52p1



      SCASE = 'U2 ; SUBCASE 2'
      CALL CALCU1A2p1              
      SCASE = 'U2 ; SUBCASE 2'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN) THEN                             
        IF (RH.GE.DRMM1) THEN
           SCASE = 'U2 ; SUBCASE 1'
           CALL CALCU2A2p1
           SCASE = 'U2 ; SUBCASE 1'
        ENDIF

      ELSE IF (.NOT.EXAN) THEN                   
        IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
           IF (RH.GE.DRMM2) THEN
              SCASE = 'U2 ; SUBCASE 3'
               CALL CALCMDRPII2p1 (RH, DRMM2, DRNANO3, CALCU1A2p1, CALCU3A2p1)
               SCASE = 'U2 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR1) THEN
               SCASE = 'U2 ; SUBCASE 4'
               CALL CALCMDRPII2p1 (RH, DRMR1, DRNANO3, CALCU1A2p1, CALCU3A2p1)
               SCASE = 'U2 ; SUBCASE 4'
            ENDIF

         ELSE IF (.NOT.EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN
               SCASE = 'U2 ; SUBCASE 5'
               CALL CALCMDRPII2p1 (RH, DRMR2, DRNACL, CALCU1A2p1, CALCU4A2p1)
               SCASE = 'U2 ; SUBCASE 5'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR3) THEN
               SCASE = 'U2 ; SUBCASE 6'
               CALL CALCMDRPII2p1 (RH, DRMR3, DRNANO3, CALCU1A2p1, CALCU3A2p1)
               SCASE = 'U2 ; SUBCASE 6'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN
               SCASE = 'U2 ; SUBCASE 7'
               CALL CALCMDRPII2p1 (RH, DRMR4, DRNACL, CALCU1A2p1, CALCU4A2p1)
               SCASE = 'U2 ; SUBCASE 7'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN
               SCASE = 'U2 ; SUBCASE 8'
               CALL CALCMDRPII2p1 (RH, DRMR5, DRNH4CL, CALCU1A2p1, CALCU52p1)
               SCASE = 'U2 ; SUBCASE 8'
            ENDIF

         ELSE IF (     EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
           IF (RH.GE.DRMR6) THEN
              SCASE = 'U2 ; SUBCASE 9'
               CALL CALCMDRPII2p1 (RH, DRMR6, DRNANO3, CALCU1A2p1, CALCU3A2p1)
               SCASE = 'U2 ; SUBCASE 9'
            ENDIF
         ENDIF

      ENDIF

      RETURN


























      END



















      SUBROUTINE CALCU2A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV7, PSCONV1, PSCONV4, PSCONV3, PSCONV5
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV7 =.FALSE.
      PSCONV1 =.FALSE.
      PSCONV4 =.FALSE.
      PSCONV3 =.FALSE.
      PSCONV5 =.FALSE.

      PSI70   =-GREAT                                 
      PSI1O   =-GREAT
      PSI40   =-GREAT
      PSI30   =-GREAT
      PSI50   =-GREAT

      ROOT7   = ZERO
      ROOT1   = ZERO
      ROOT4   = ZERO
      ROOT3   = ZERO
      ROOT5   = ZERO



      CALL CALCU1A2p1

      CHI1   = CNA2SO4            
      CHI2   = CNANO3
      CHI3   = CNACL
      CHI4   = CNH4CL
      CHI7   = CK2SO4
      CHI8   = CMGSO4
      CHI9   = CCASO4

      PSI1   = CNA2SO4      
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

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
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI

      CALL CALCACT2p1          



      DO 10 I=1,NSWEEP

      A7  = XK17 *(WATER/GAMA(17))**3.0      
      A1  = XK5 *(WATER/GAMA(2))**3.0        
      A14 = XK14*(WATER/GAMA(6))**2.0        
      A8  = XK8 *(WATER/GAMA(1))**2.0        
      A9  = XK9 *(WATER/GAMA(3))**2.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A7) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
         CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
         CALL POLY32p1(BB, CC, DD, ROOT7, ISLV)
         IF (ISLV.NE.0) ROOT7 = TINY
         ROOT7 = MAX (ROOT7, ZERO)
         ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
                      MAX(WAER(2)-WAER(6)-ROOT1, ZERO),CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI70) .LE. EPS*PSI70
      PSI70   = PSI7



      IF (NAI*NAI*SO4I .GT. A1) THEN
         BB =-(((WAER(2)-WAER(6))-ROOT7)*(WAER(1) - ROOT3 - ROOT5))
         CC = ((WAER(2)-WAER(6)) - ROOT7)*(WAER(1) - ROOT3 - ROOT5) + &
               0.25D0*(WAER(1) - ROOT3 - ROOT5)**2.0
         DD =-0.25D0*(((WAER(2) - WAER(6)) - ROOT7)* &
                      (WAER(1) - ROOT3 - ROOT5)**2.D0 - A1)
         CALL POLY32p1(BB, CC, DD, ROOT1, ISLV)
         IF (ISLV.NE.0) ROOT1 = TINY
         ROOT1 = MIN (MAX(ROOT1,ZERO), MAX(WAER(1)-ROOT3-ROOT5,ZERO), &
                      CHI1, MAX(WAER(2)-WAER(6),ZERO))
         PSI1  = CHI1-ROOT1
      ENDIF
      PSCONV1 = ABS(PSI1-PSI1O) .LE. EPS*PSI1O
      PSI1O   = PSI1



      IF (NH4I*CLI .GT. A14) THEN
         BB    =-(WAER(3) + WAER(5) - ROOT4)
         CC    =-A14 + NH4I*(WAER(5) - ROOT4)
         DD    = MAX(BB*BB - 4.D0*CC, ZERO)
         ROOT4A= 0.5D0*(-BB+SQRT(DD))
         ROOT4B= 0.5D0*(-BB-SQRT(DD))
         IF (ZERO.LE.ROOT4A) THEN
            ROOT4 = ROOT4A
         ELSE
            ROOT4 = ROOT4B
         ENDIF
         ROOT4 = MIN(MAX(ZERO, ROOT4), MAX(WAER(5)-ROOT3,ZERO), &
                     CHI4, WAER(3))
         PSI4  = CHI4 - ROOT4
      ENDIF
      PSCONV4 = ABS(PSI4-PSI40) .LE. EPS*PSI40
      PSI40   = PSI4



      IF (NAI*CLI .GT. A8) THEN
         BB    =-((CHI1-2.D0*ROOT1-ROOT5) + (WAER(5) - ROOT4))
         CC    = (CHI1-2.D0*ROOT1-ROOT5)*(WAER(5) - ROOT4) - A8
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
      PSCONV3 = ABS(PSI3-PSI30) .LE. EPS*PSI30
      PSI30   = PSI3



      IF (NAI*NO3I .GT. A9) THEN
         BB    =-(WAER(4) + WAER(1) - 2.D0*ROOT1 - ROOT3)
         CC    = WAER(4)*(WAER(1) - 2.D0*ROOT1 - ROOT3) - A9
         DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT5A= 0.5D0*(-BB-DD)
         ROOT5B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT5A) THEN
            ROOT5 = ROOT5A
         ELSE
            ROOT5 = ROOT5B
         ENDIF
         ROOT5 = MIN(MAX(ROOT5, ZERO), CHI2)
         PSI2  = CHI2-ROOT5
      ENDIF

      PSCONV5 = ABS(PSI2-PSI20) .LE. EPS*PSI20
      PSI20   = PSI2



      KI     = MAX (WAER(7) - 2.0D0*ROOT7, ZERO)
      NAI    = MAX (WAER(1) - 2.0D0*ROOT1 - ROOT3 - ROOT5, ZERO)
      SO4I   = MAX (WAER(2) - WAER(6) - ROOT7 - ROOT1, ZERO)
      NH4I   = MAX (WAER(3) - ROOT4, ZERO)
      NO3I   = MAX (WAER(4) - ROOT5, ZERO)
      CLI    = MAX (WAER(5) - ROOT4 - ROOT3, ZERO)
      CAI    = ZERO
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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





      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV4 .AND. PSCONV3 &
              .AND. PSCONV5) GOTO 20
      ENDIF
10    CONTINUE




20    IF (CLI.LE.TINY .AND. WAER(5).GT.TINY) THEN 
         DO 30 I=1,NIONS
            MOLAL(I) = ZERO
30       CONTINUE
         DO 40 I=1,NGASAQ
            GASAQ(I) = ZERO
40       CONTINUE
         CALL CALCU1A2p1
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
      CMGSO4  = ZERO
      CK2SO4  = CHI7 - PSI7
      CCASO4  = MIN (WAER(6), WAER(2))
      ENDIF

      RETURN



      END


















      SUBROUTINE CALCU12p1
      INCLUDE 'module_isrpia_inc.F'
      LOGICAL  EXAN, EXAC, EXSN, EXSC
      EXTERNAL CALCU1A2p1, CALCU2A2p1, CALCU3A2p1, CALCU4A2p1, CALCU52p1



      SCASE = 'U1 ; SUBCASE 1'
      CALL CALCU1A2p1              
      SCASE = 'U1 ; SUBCASE 1'

      EXAN = CNH4NO3.GT.TINY    
      EXAC = CNH4CL .GT.TINY    
      EXSN = CNANO3 .GT.TINY    
      EXSC = CNACL  .GT.TINY    



      IF (EXAN.OR.EXAC.OR.EXSC.OR.EXSN) THEN  
         IF (RH.GE.DRMM1) THEN
            SCASE = 'U1 ; SUBCASE 2'  
            CALL CALCMDRPII2p1 (RH, DRMM1, DRNH4NO3, CALCU1A2p1, CALCU2A2p1)
            SCASE = 'U1 ; SUBCASE 2'
         ENDIF

      ELSE IF (.NOT.EXAN) THEN                   
         IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMM2) THEN
               SCASE = 'U1 ; SUBCASE 3'
               CALL CALCMDRPII2p1 (RH, DRMM2, DRNANO3, CALCU1A2p1, CALCU3A2p1)
               SCASE = 'U1 ; SUBCASE 3'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR1) THEN
              SCASE = 'U1 ; SUBCASE 4'
               CALL CALCMDRPII2p1 (RH, DRMR1, DRNANO3, CALCU1A2p1, CALCU3A2p1)
               SCASE = 'U1 ; SUBCASE 4'
            ENDIF

         ELSE IF (.NOT.EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR2) THEN
               SCASE = 'U1 ; SUBCASE 5'
               CALL CALCMDRPII2p1 (RH, DRMR2, DRNACL, CALCU1A2p1, CALCU3A2p1) 
               SCASE = 'U1 ; SUBCASE 5'
            ENDIF

         ELSE IF (.NOT.EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR3) THEN
               SCASE = 'U1 ; SUBCASE 6'
               CALL CALCMDRPII2p1 (RH, DRMR3, DRNANO3, CALCU1A2p1, CALCU3A2p1)
               SCASE = 'U1 ; SUBCASE 6'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR4) THEN
               SCASE = 'U1 ; SUBCASE 7'
               CALL CALCMDRPII2p1 (RH, DRMR4, DRNACL, CALCU1A2p1, CALCU3A2p1) 
               SCASE = 'U1 ; SUBCASE 7'
            ENDIF

         ELSE IF (     EXAC .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR5) THEN
               SCASE = 'U1 ; SUBCASE 8'
               CALL CALCMDRPII2p1 (RH, DRMR5, DRNH4CL, CALCU1A2p1, CALCU3A2p1) 
               SCASE = 'U1 ; SUBCASE 8'
            ENDIF

         ELSE IF (     EXAC .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR6) THEN
               SCASE = 'U1 ; SUBCASE 9'
               CALL CALCMDRPII2p1 (RH, DRMR6, DRNANO3, CALCU1A2p1, CALCU3A2p1)
               SCASE = 'U1 ; SUBCASE 9'
            ENDIF
         ENDIF

      ELSE IF (.NOT.EXAC) THEN                   
         IF      (     EXAN .AND.      EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR7) THEN
               SCASE = 'U1 ; SUBCASE 10'
               CALL CALCMDRPII2p1 (RH, DRMR7, DRNH4NO3, CALCU1A2p1, CALCU2A2p1)
               SCASE = 'U1 ; SUBCASE 10'
            ENDIF

         ELSE IF (     EXAN .AND. .NOT.EXSN .AND.      EXSC) THEN
            IF (RH.GE.DRMR8) THEN
              SCASE = 'U1 ; SUBCASE 11'
              CALL CALCMDRPII2p1 (RH, DRMR8, DRNH4NO3, CALCU1A2p1, CALCU2A2p1)
               SCASE = 'U1 ; SUBCASE 11'
            ENDIF

         ELSE IF (     EXAN .AND. .NOT.EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR9) THEN
               SCASE = 'U1 ; SUBCASE 12'
               CALL CALCMDRPII2p1 (RH, DRMR9, DRNH4NO3, CALCU1A2p1, CALCU2A2p1)
               SCASE = 'U1 ; SUBCASE 12'
            ENDIF

         ELSE IF (     EXAN .AND.      EXSN .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR10) THEN
               SCASE = 'U1 ; SUBCASE 13'
               CALL CALCMDRPII2p1 (RH, DRMR10, DRNH4NO3, CALCU1A2p1, CALCU2A2p1)
               SCASE = 'U1 ; SUBCASE 13'
           ENDIF
        ENDIF

      ELSE IF (.NOT.EXSN) THEN                  
         IF      (     EXAN .AND.      EXAC .AND.      EXSC) THEN
            IF (RH.GE.DRMR11) THEN
               SCASE = 'U1 ; SUBCASE 14'
               CALL CALCMDRPII2p1 (RH, DRMR11, DRNH4NO3, CALCU1A2p1, CALCU2A2p1)
              SCASE = 'U1 ; SUBCASE 14'
           ENDIF

        ELSE IF (     EXAN .AND.      EXAC .AND. .NOT.EXSC) THEN
            IF (RH.GE.DRMR12) THEN
               SCASE = 'U1 ; SUBCASE 15'
               CALL CALCMDRPII2p1 (RH, DRMR12, DRNH4NO3, CALCU1A2p1, CALCU2A2p1)
               SCASE = 'U1 ; SUBCASE 15'
            ENDIF
         ENDIF

      ELSE IF (.NOT.EXSC) THEN                  
         IF      (     EXAN .AND.      EXAC .AND.      EXSN) THEN
            IF (RH.GE.DRMR13) THEN
               SCASE = 'U1 ; SUBCASE 16'
               CALL CALCMDRPII2p1 (RH, DRMR13, DRNH4NO3, CALCU1A2p1, CALCU2A2p1)
               SCASE = 'U1 ; SUBCASE 16'
            ENDIF
         ENDIF
      ENDIF

      RETURN
















      END

















      SUBROUTINE CALCU1A2p1
      INCLUDE 'module_isrpia_inc.F'



      CCASO4  = MIN (WAER(6), WAER(2))                 
      SO4FR   = MAX(WAER(2) - CCASO4, ZERO)
      CAFR    = MAX(WAER(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*WAER(7), SO4FR)             
      FRK     = MAX(WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX(SO4FR - CK2SO4, ZERO)
      CMGSO4  = MIN (WAER(8), SO4FR)                   
      FRMG    = MAX(WAER(8) - CMGSO4, ZERO)
      SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
      CNA2SO4 = MAX (SO4FR, ZERO)                      
      FRNA    = MAX (WAER(1) - 2.D0*CNA2SO4, ZERO)

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



      END




















      SUBROUTINE CALCW132p1
      INCLUDE 'module_isrpia_inc.F'

      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.



      CALL CALCW1A2p1

      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      AKW = XKW*RH*WATER*WATER               



      NAI    = WAER(1)
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      KI     = WAER(7)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
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
      CMGSO4  = ZERO
      CK2SO4  = ZERO
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = ZERO
      KCL     = ZERO
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



















      SUBROUTINE CALCW122p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSI9O   =-GREAT                                 
      ROOT9   = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7))
         CC = WAER(7)*(WAER(2)-WAER(6)) + 0.25D0*WAER(7)*WAER(7)
         DD =-0.25*(WAER(7)*WAER(7)*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0, (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      KI     = MAX (WAER(7) - 2.D0*ROOT9, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      NAI    = WAER(1)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9) GOTO 20
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
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = ZERO
      KCL     = ZERO
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



















      SUBROUTINE CALCW112p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT                                

      ROOT9   = ZERO
      ROOT13  = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13)
         CC = (WAER(7)-ROOT13)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13)**2.0
         DD =-0.25*((WAER(7)-ROOT13)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9,WAER(7)/2.0-ROOT13,(WAER(2)-WAER(6)),CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9)
         CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = WAER(3)
      NO3I   = MAX (WAER(4) - ROOT13, ZERO)
      CLI    = WAER(5)
      CAI    = ZERO
      NAI    = WAER(1)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13) GOTO 20
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
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = ZERO
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



















      SUBROUTINE CALCW102p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT                                

      ROOT9   = ZERO
      ROOT13  = ZERO



      CALL CALCW1A2p1


      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13)
         CC = (WAER(7)-ROOT13)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13)**2.0
         DD =-0.25*((WAER(7)-ROOT13)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9,WAER(7)/2.0-ROOT13,(WAER(2)-WAER(6)),CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9)
         CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = WAER(3)
      NO3I   = MAX (WAER(4) - ROOT13, ZERO)
      CLI    = WAER(5)
      CAI    = ZERO
      NAI    = WAER(1)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13) GOTO 20
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
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = ZERO
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



















      SUBROUTINE CALCW92p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13, PSCONV14
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.
      PSCONV14=.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT
      PSI14O  =-GREAT                              

      ROOT9   = ZERO
      ROOT13  = ZERO
      ROOT14  = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI14  = CKCL
      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      A14 = XK20 *(WATER/GAMA(20))**2.0      
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
         CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
         DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
                      (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9 - ROOT14)
         CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      IF (KI*CLI .GT. A14) THEN
         BB     =-(WAER(5) + WAER(7) - 2.D0*ROOT9 - ROOT13)
         CC     = WAER(5)*(WAER(7) - 2.D0*ROOT9 - ROOT13) - A14
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT14A= 0.5D0*(-BB-DD)
         ROOT14B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT14A) THEN
            ROOT14 = ROOT14A
         ELSE
            ROOT14 = ROOT14B
         ENDIF
         ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
         PSI14  = CHI14-ROOT14
      ENDIF
      PSCONV14 = ABS(PSI14-PSI14O) .LE. EPS*PSI14O
      PSI14O   = PSI14



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = WAER(3)
      NO3I   = MAX (WAER(4) - ROOT13, ZERO)
      CLI    = MAX (WAER(5) - ROOT14, ZERO)
      CAI    = ZERO
      NAI    = WAER(1)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14) GOTO 20
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
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = CHI14 - PSI14
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



















      SUBROUTINE CALCW82p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13, PSCONV14, PSCONV5
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.
      PSCONV14=.TRUE.
      PSCONV5 =.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT
      PSI14O  =-GREAT
      PSI5O   =-GREAT                              

      ROOT9   = ZERO
      ROOT13  = ZERO
      ROOT14  = ZERO
      ROOT5   = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI14  = CKCL
      CHI5   = CNH4CL
      CHI11  = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      A14 = XK20 *(WATER/GAMA(20))**2.0      
      A5  = XK14*(WATER/GAMA(6))**2.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
         CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
         DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
                      (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9 - ROOT14)
         CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      IF (KI*CLI .GT. A14) THEN
         BB     =-(WAER(5) - ROOT5 + WAER(7) - 2.D0*ROOT9 - ROOT13)
         CC     = (WAER(5)-ROOT5)*(WAER(7) - 2.D0*ROOT9 - ROOT13) - A14
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT14A= 0.5D0*(-BB-DD)
         ROOT14B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT14A) THEN
            ROOT14 = ROOT14A
         ELSE
            ROOT14 = ROOT14B
         ENDIF
         ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
         PSI14  = CHI14-ROOT14
      ENDIF
      PSCONV14 = ABS(PSI14-PSI14O) .LE. EPS*PSI14O
      PSI14O   = PSI14



      IF (NH4I*CLI .GT. A5) THEN
         BB     =-(WAER(5) + WAER(3) - ROOT14)
         CC     = (WAER(5)-ROOT14)*WAER(3) - A5
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT5A = 0.5D0*(-BB-DD)
         ROOT5B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT5A) THEN
            ROOT5 = ROOT5A
         ELSE
            ROOT5 = ROOT5B
         ENDIF
         ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
         PSI5  = CHI5-ROOT5
      ENDIF
      PSCONV5 = ABS(PSI5-PSI5O) .LE. EPS*PSI5O
      PSI5O   = PSI5



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = MAX (WAER(3) - ROOT5, ZERO)
      NO3I   = MAX (WAER(4) - ROOT13, ZERO)
      CLI    = MAX (WAER(5) - ROOT14 - ROOT5, ZERO)
      CAI    = ZERO
      NAI    = WAER(1)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND.PSCONV5) GOTO 20
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
      CNH4CL  = CHI5 - PSI5
      CNACL   = ZERO
      CNANO3  = ZERO
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = CHI14 - PSI14
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



















      SUBROUTINE CALCW72p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.
      PSCONV14=.TRUE.
      PSCONV5 =.TRUE.
      PSCONV7 =.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT
      PSI14O  =-GREAT
      PSI5O  =-GREAT
      PSI7O  =-GREAT                            

      ROOT9   = ZERO
      ROOT13  = ZERO
      ROOT14  = ZERO
      ROOT5   = ZERO
      ROOT7   = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI14  = CKCL
      CHI5   = CNH4CL
      CHI7   = CNACL
      CHI11  = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      A14 = XK20 *(WATER/GAMA(20))**2.0      
      A5  = XK14*(WATER/GAMA(6))**2.0        
      A7  = XK8 *(WATER/GAMA(1))**2.0        
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
         CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
         DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
                      (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9 - ROOT14)
         CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      IF (KI*CLI .GT. A14) THEN
         BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
         CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT14A= 0.5D0*(-BB-DD)
         ROOT14B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT14A) THEN
            ROOT14 = ROOT14A
         ELSE
            ROOT14 = ROOT14B
         ENDIF
         ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
         PSI14  = CHI14-ROOT14
      ENDIF
      PSCONV14 = ABS(PSI14-PSI14O) .LE. EPS*PSI14O
      PSI14O   = PSI14



      IF (NH4I*CLI .GT. A5) THEN
         BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
         CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT5A = 0.5D0*(-BB-DD)
         ROOT5B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT5A) THEN
            ROOT5 = ROOT5A
         ELSE
            ROOT5 = ROOT5B
         ENDIF
         ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
         PSI5  = CHI5-ROOT5
      ENDIF
      PSCONV5 = ABS(PSI5-PSI5O) .LE. EPS*PSI5O
      PSI5O   = PSI5



      IF (NAI*CLI .GT. A7) THEN
         BB     =-(WAER(5) + WAER(1) - ROOT14 - ROOT5)
         CC     = (WAER(5) - ROOT14 - ROOT5)*WAER(1) - A7
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT7A = 0.5D0*(-BB-DD)
         ROOT7B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT7A) THEN
            ROOT7 = ROOT7A
         ELSE
            ROOT7 = ROOT7B
         ENDIF
         ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI7O) .LE. EPS*PSI7O
      PSI7O   = PSI7



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = MAX (WAER(3) - ROOT5, ZERO)
      NO3I   = MAX (WAER(4) - ROOT13, ZERO)
      CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
      CAI    = ZERO
      NAI    = MAX (WAER(1) - ROOT7, ZERO)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7) GOTO 20
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
      CNH4CL  = CHI5 - PSI5
      CNACL   = CHI7 - PSI7
      CNANO3  = ZERO
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = CHI14 - PSI14
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END




















      SUBROUTINE CALCW62p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.
      PSCONV14=.TRUE.
      PSCONV5 =.TRUE.
      PSCONV7 =.TRUE.
      PSCONV8 =.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT
      PSI14O  =-GREAT
      PSI5O   =-GREAT
      PSI7O   =-GREAT
      PSI8O   =-GREAT                     

      ROOT9   = ZERO
      ROOT13  = ZERO
      ROOT14  = ZERO
      ROOT5   = ZERO
      ROOT7   = ZERO
      ROOT8   = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI14  = CKCL
      CHI5   = CNH4CL
      CHI7   = CNACL
      CHI8   = CNANO3
      CHI11  = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      A14 = XK20 *(WATER/GAMA(20))**2.0      
      A5  = XK14*(WATER/GAMA(6))**2.0        
      A7  = XK8 *(WATER/GAMA(1))**2.0        
      A8  = XK9 *(WATER/GAMA(3))**2.         
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
         CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
         DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
                     (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
         CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      IF (KI*CLI .GT. A14) THEN
         BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
         CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT14A= 0.5D0*(-BB-DD)
         ROOT14B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT14A) THEN
            ROOT14 = ROOT14A
         ELSE
            ROOT14 = ROOT14B
         ENDIF
         ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
         PSI14  = CHI14-ROOT14
      ENDIF
      PSCONV14 = ABS(PSI14-PSI14O) .LE. EPS*PSI14O
      PSI14O   = PSI14



      IF (NH4I*CLI .GT. A5) THEN
         BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
         CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT5A = 0.5D0*(-BB-DD)
         ROOT5B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT5A) THEN
            ROOT5 = ROOT5A
         ELSE
            ROOT5 = ROOT5B
         ENDIF
         ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
         PSI5  = CHI5-ROOT5
      ENDIF
      PSCONV5 = ABS(PSI5-PSI5O) .LE. EPS*PSI5O
      PSI5O   = PSI5



      IF (NAI*CLI .GT. A7) THEN
         BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
         CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT7A = 0.5D0*(-BB-DD)
         ROOT7B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT7A) THEN
            ROOT7 = ROOT7A
         ELSE
            ROOT7 = ROOT7B
         ENDIF
         ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI7O) .LE. EPS*PSI7O
      PSI7O   = PSI7



      IF (NAI*NO3I .GT. A8) THEN
         BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
         CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT8A = 0.5D0*(-BB-DD)
         ROOT8B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT8A) THEN
            ROOT8 = ROOT8A
         ELSE
            ROOT8 = ROOT8B
         ENDIF
         ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
         PSI8  = CHI8-ROOT8
      ENDIF
      PSCONV8 = ABS(PSI8-PSI8O) .LE. EPS*PSI8O
      PSI8O   = PSI8



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = MAX (WAER(3) - ROOT5, ZERO)
      NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
      CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
      CAI    = ZERO
      NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
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
      CNH4CL  = CHI5 - PSI5
      CNACL   = CHI7 - PSI7
      CNANO3  = CHI8 - PSI8
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = CHI14 - PSI14
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END




















      SUBROUTINE CALCW52p1
      INCLUDE 'module_isrpia_inc.F'

      EXTERNAL CALCW1A2p1, CALCW62p1



      IF (WAER(4).GT.TINY)   THEN 
         SCASE = 'W5 ; SUBCASE 1'
         CALL CALCW5A2p1
         SCASE = 'W5 ; SUBCASE 1'
      ELSE                                      
         SCASE = 'W1 ; SUBCASE 1'
         CALL CALCW1A2p1
         SCASE = 'W1 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMP5) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCW1A2p1
            SCASE = 'W5 ; SUBCASE 2'
            RETURN
         ELSE
            SCASE = 'W5 ; SUBCASE 3'  

            CALL CALCMDRPII2p1 (RH, DRMP5, DRNH4NO3, CALCW1A2p1, CALCW62p1)
            SCASE = 'W5 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCW5A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.
      PSCONV14=.TRUE.
      PSCONV5 =.TRUE.
      PSCONV7 =.TRUE.
      PSCONV8 =.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT
      PSI14O  =-GREAT
      PSI5O   =-GREAT
      PSI7O   =-GREAT
      PSI8O   =-GREAT                

      ROOT9   = ZERO
      ROOT13  = ZERO
      ROOT14  = ZERO
      ROOT5   = ZERO
      ROOT7   = ZERO
      ROOT8   = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI14  = CKCL
      CHI5   = CNH4CL
      CHI7   = CNACL
      CHI8   = CNANO3
      CHI6   = CNH4NO3
      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      A14 = XK20 *(WATER/GAMA(20))**2.0      
      A5  = XK14*(WATER/GAMA(6))**2.0        
      A7  = XK8 *(WATER/GAMA(1))**2.0        
      A8  = XK9 *(WATER/GAMA(3))**2.         
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
         CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
         DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
                      (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
         CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      IF (KI*CLI .GT. A14) THEN
         BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
         CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT14A= 0.5D0*(-BB-DD)
         ROOT14B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT14A) THEN
            ROOT14 = ROOT14A
         ELSE
            ROOT14 = ROOT14B
         ENDIF
         ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
         PSI14  = CHI14-ROOT14
      ENDIF
      PSCONV14 = ABS(PSI14-PSI14O) .LE. EPS*PSI14O
      PSI14O   = PSI14



      IF (NH4I*CLI .GT. A5) THEN
         BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
         CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT5A = 0.5D0*(-BB-DD)
         ROOT5B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT5A) THEN
            ROOT5 = ROOT5A
         ELSE
            ROOT5 = ROOT5B
         ENDIF
         ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
         PSI5  = CHI5-ROOT5
      ENDIF
      PSCONV5 = ABS(PSI5-PSI5O) .LE. EPS*PSI5O
      PSI5O   = PSI5



      IF (NAI*CLI .GT. A7) THEN
         BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
         CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT7A = 0.5D0*(-BB-DD)
         ROOT7B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT7A) THEN
            ROOT7 = ROOT7A
         ELSE
            ROOT7 = ROOT7B
         ENDIF
         ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI7O) .LE. EPS*PSI7O
      PSI7O   = PSI7



      IF (NAI*NO3I .GT. A8) THEN
         BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
         CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT8A = 0.5D0*(-BB-DD)
         ROOT8B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT8A) THEN
            ROOT8 = ROOT8A
         ELSE
            ROOT8 = ROOT8B
         ENDIF
         ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
         PSI8  = CHI8-ROOT8
      ENDIF
      PSCONV8 = ABS(PSI8-PSI8O) .LE. EPS*PSI8O
      PSI8O   = PSI8



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = MAX (WAER(3) - ROOT5, ZERO)
      NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
      CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
      CAI    = ZERO
      NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
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
      CNH4CL  = CHI5 - PSI5
      CNACL   = CHI7 - PSI7
      CNANO3  = CHI8 - PSI8
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = CHI14 - PSI14
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



















      SUBROUTINE CALCW42p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCW1A2p1, CALCW5A2p1



      IF (WAER(4).GT.TINY)   THEN 
         SCASE = 'W4 ; SUBCASE 1'
         CALL CALCW4A2p1
         SCASE = 'W4 ; SUBCASE 1'
      ELSE                                      
         SCASE = 'W1 ; SUBCASE 1'
         CALL CALCW1A2p1
         SCASE = 'W1 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMP4) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCW1A2p1
            SCASE = 'W4 ; SUBCASE 2'
            RETURN
         ELSE
            SCASE = 'W4 ; SUBCASE 3'  

            CALL CALCMDRPII2p1 (RH, DRMP4, DRMGNO32, CALCW1A2p1, CALCW5A2p1)
            SCASE = 'W4 ; SUBCASE 3'
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCW4A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.
      PSCONV14=.TRUE.
      PSCONV5 =.TRUE.
      PSCONV7 =.TRUE.
      PSCONV8 =.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT
      PSI14O  =-GREAT
      PSI5O   =-GREAT
      PSI7O   =-GREAT
      PSI8O   =-GREAT                

      ROOT9   = ZERO
      ROOT13  = ZERO
      ROOT14  = ZERO
      ROOT5   = ZERO
      ROOT7   = ZERO
      ROOT8   = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI14  = CKCL
      CHI5   = CNH4CL
      CHI7   = CNACL
      CHI8   = CNANO3
      CHI6   = CNH4NO3
      CHI15  = CMGNO32
      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      A14 = XK20 *(WATER/GAMA(20))**2.0      
      A5  = XK14*(WATER/GAMA(6))**2.0        
      A7  = XK8 *(WATER/GAMA(1))**2.0        
      A8  = XK9 *(WATER/GAMA(3))**2.         
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
         CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
         DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
                      (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
         CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      IF (KI*CLI .GT. A14) THEN
         BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
         CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT14A= 0.5D0*(-BB-DD)
         ROOT14B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT14A) THEN
            ROOT14 = ROOT14A
         ELSE
            ROOT14 = ROOT14B
         ENDIF
         ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
         PSI14  = CHI14-ROOT14
      ENDIF
      PSCONV14 = ABS(PSI14-PSI14O) .LE. EPS*PSI14O
      PSI14O   = PSI14



      IF (NH4I*CLI .GT. A5) THEN
         BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
         CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT5A = 0.5D0*(-BB-DD)
         ROOT5B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT5A) THEN
            ROOT5 = ROOT5A
         ELSE
            ROOT5 = ROOT5B
         ENDIF
         ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
         PSI5  = CHI5-ROOT5
      ENDIF
      PSCONV5 = ABS(PSI5-PSI5O) .LE. EPS*PSI5O
      PSI5O   = PSI5



      IF (NAI*CLI .GT. A7) THEN
         BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
         CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT7A = 0.5D0*(-BB-DD)
         ROOT7B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT7A) THEN
            ROOT7 = ROOT7A
         ELSE
            ROOT7 = ROOT7B
         ENDIF
         ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI7O) .LE. EPS*PSI7O
      PSI7O   = PSI7



      IF (NAI*NO3I .GT. A8) THEN
         BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
         CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT8A = 0.5D0*(-BB-DD)
         ROOT8B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT8A) THEN
            ROOT8 = ROOT8A
         ELSE
            ROOT8 = ROOT8B
         ENDIF
         ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
         PSI8  = CHI8-ROOT8
      ENDIF
      PSCONV8 = ABS(PSI8-PSI8O) .LE. EPS*PSI8O
      PSI8O   = PSI8



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = MAX (WAER(3) - ROOT5, ZERO)
      NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
      CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
      CAI    = ZERO
      NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
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
      CNH4CL  = CHI5 - PSI5
      CNACL   = CHI7 - PSI7
      CNANO3  = CHI8 - PSI8
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = CHI14 - PSI14
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



















      SUBROUTINE CALCW32p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCW1A2p1, CALCW4A2p1













      CALL CALCW1A2p1
      
      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMP3) THEN        
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCW1A2p1
            SCASE = 'W3 ; SUBCASE 2'
            RETURN
         ELSE
            SCASE = 'W3 ; SUBCASE 3'  

            CALL CALCMDRPII2p1 (RH, DRMP3, DRCANO32, CALCW1A2p1, CALCW4A2p1)
            SCASE = 'W3 ; SUBCASE 3'
         ENDIF
      ELSE                                      
         SCASE = 'W3 ; SUBCASE 1'
         CALL CALCW3A2p1
         SCASE = 'W3 ; SUBCASE 1'
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCW3A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.
      PSCONV14=.TRUE.
      PSCONV5 =.TRUE.
      PSCONV7 =.TRUE.
      PSCONV8 =.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT
      PSI14O  =-GREAT
      PSI5O   =-GREAT
      PSI7O   =-GREAT
      PSI8O   =-GREAT                

      ROOT9   = ZERO
      ROOT13  = ZERO
      ROOT14  = ZERO
      ROOT5   = ZERO
      ROOT7   = ZERO
      ROOT8   = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI14  = CKCL
      CHI5   = CNH4CL
      CHI7   = CNACL
      CHI8   = CNANO3
      CHI6   = CNH4NO3
      CHI15  = CMGNO32
      CHI12  = CCANO32
      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      A14 = XK20 *(WATER/GAMA(20))**2.0      
      A5  = XK14*(WATER/GAMA(6))**2.0        
      A7  = XK8 *(WATER/GAMA(1))**2.0        
      A8  = XK9 *(WATER/GAMA(3))**2.         
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
         CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
         DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
                      (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
         CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      IF (KI*CLI .GT. A14) THEN
         BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
         CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT14A= 0.5D0*(-BB-DD)
         ROOT14B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT14A) THEN
            ROOT14 = ROOT14A
         ELSE
            ROOT14 = ROOT14B
         ENDIF
         ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
         PSI14  = CHI14-ROOT14
      ENDIF
      PSCONV14 = ABS(PSI14-PSI14O) .LE. EPS*PSI14O
      PSI14O   = PSI14



      IF (NH4I*CLI .GT. A5) THEN
         BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
         CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT5A = 0.5D0*(-BB-DD)
         ROOT5B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT5A) THEN
            ROOT5 = ROOT5A
         ELSE
            ROOT5 = ROOT5B
         ENDIF
         ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
         PSI5  = CHI5-ROOT5
      ENDIF
      PSCONV5 = ABS(PSI5-PSI5O) .LE. EPS*PSI5O
      PSI5O   = PSI5



      IF (NAI*CLI .GT. A7) THEN
         BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
         CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT7A = 0.5D0*(-BB-DD)
         ROOT7B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT7A) THEN
            ROOT7 = ROOT7A
         ELSE
            ROOT7 = ROOT7B
         ENDIF
         ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI7O) .LE. EPS*PSI7O
      PSI7O   = PSI7



      IF (NAI*NO3I .GT. A8) THEN
         BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
         CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT8A = 0.5D0*(-BB-DD)
         ROOT8B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT8A) THEN
            ROOT8 = ROOT8A
         ELSE
            ROOT8 = ROOT8B
         ENDIF
         ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
         PSI8  = CHI8-ROOT8
      ENDIF
      PSCONV8 = ABS(PSI8-PSI8O) .LE. EPS*PSI8O
      PSI8O   = PSI8



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = MAX (WAER(3) - ROOT5, ZERO)
      NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
      CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
      CAI    = ZERO
      NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
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
      CNH4CL  = CHI5 - PSI5
      CNACL   = CHI7 - PSI7
      CNANO3  = CHI8 - PSI8
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = CHI14 - PSI14
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END



























      SUBROUTINE CALCW22p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCW1A2p1, CALCW3A2p1, CALCW4A2p1, CALCW5A2p1, CALCW62p1



      CALL CALCW1A2p1



      IF (CCACL2.GT.TINY) THEN
         SCASE = 'W2 ; SUBCASE 1'
         CALL CALCW2A2p1
         SCASE = 'W2 ; SUBCASE 1'
      ENDIF

      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMP2) THEN             
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL CALCW1A2p1
            SCASE = 'W2 ; SUBCASE 2'
         ELSE
            IF (CMGCL2.GT. TINY) THEN
               SCASE = 'W2 ; SUBCASE 3'    

               CALL CALCMDRPII2p1 (RH, DRMP2, DRMGCL2, CALCW1A2p1, CALCW3A2p1)
               SCASE = 'W2 ; SUBCASE 3'
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMP3 .AND. RH.LT.DRMP4) THEN
               SCASE = 'W2 ; SUBCASE 4'    

               CALL CALCMDRPII2p1 (RH, DRMP3, DRCANO32, CALCW1A2p1, CALCW4A2p1)
               SCASE = 'W2 ; SUBCASE 4'
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMP4 .AND. RH.LT.DRMP5) THEN
               SCASE = 'W2 ; SUBCASE 5'    

               CALL CALCMDRPII2p1 (RH, DRMP4, DRMGNO32, CALCW1A2p1, CALCW5A2p1)
               SCASE = 'W2 ; SUBCASE 5'
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMP5) THEN
               SCASE = 'W2 ; SUBCASE 6'    

               CALL CALCMDRPII2p1 (RH, DRMP5, DRNH4NO3, CALCW1A2p1, CALCW62p1)
               SCASE = 'W2 ; SUBCASE 6'
            ELSE
               WATER = TINY
               DO 20 I=1,NIONS
                  MOLAL(I) = ZERO
20             CONTINUE
               CALL CALCW1A2p1
               SCASE = 'W2 ; SUBCASE 2'
            ENDIF
         ENDIF
      ENDIF

      RETURN



      END




















      SUBROUTINE CALCW2A2p1
      INCLUDE 'module_isrpia_inc.F'

      LOGICAL PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

      COMMON /SOLUT2p1/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                     CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                     CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                     PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                     PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                     A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17



      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.

      PSCONV9 =.TRUE.
      PSCONV13=.TRUE.
      PSCONV14=.TRUE.
      PSCONV5 =.TRUE.
      PSCONV7 =.TRUE.
      PSCONV8 =.TRUE.

      PSI9O   =-GREAT
      PSI13O  =-GREAT
      PSI14O  =-GREAT
      PSI5O   =-GREAT
      PSI7O   =-GREAT
      PSI8O   =-GREAT                

      ROOT9   = ZERO
      ROOT13  = ZERO
      ROOT14  = ZERO
      ROOT5   = ZERO
      ROOT7   = ZERO
      ROOT8   = ZERO



      CALL CALCW1A2p1

      CHI9   = CK2SO4       
      CHI13  = CKNO3
      CHI10  = CMGSO4
      CHI14  = CKCL
      CHI5   = CNH4CL
      CHI7   = CNACL
      CHI8   = CNANO3
      CHI6   = CNH4NO3
      CHI15  = CMGNO32
      CHI12  = CCANO32
      CHI16  = CMGCL2
      CHI11   = CCASO4

      PSI1   = CNA2SO4      
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2

      CALL CALCMR2p1           

      NAI    = WAER(1)      
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = WAER(6)
      KI     = WAER(7)
      MGI    = WAER(8)

      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO



      DO 10 I=1,NSWEEP

      A9  = XK17 *(WATER/GAMA(17))**3.0      
      A13 = XK19 *(WATER/GAMA(19))**2.0      
      A14 = XK20 *(WATER/GAMA(20))**2.0      
      A5  = XK14*(WATER/GAMA(6))**2.0        
      A7  = XK8 *(WATER/GAMA(1))**2.0        
      A8  = XK9 *(WATER/GAMA(3))**2.         
      AKW = XKW*RH*WATER*WATER               



      IF (KI*KI*SO4I .GT. A9) THEN
         BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
         CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
               0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
         DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
         CALL POLY32p1(BB, CC, DD, ROOT9, ISLV)
         IF (ISLV.NE.0) ROOT9 = TINY
         ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
                      (WAER(2)-WAER(6)), CHI9)
         ROOT9 = MAX (ROOT9, ZERO)
         PSI9  = CHI9 - ROOT9
      ENDIF
      PSCONV9 = ABS(PSI9-PSI9O) .LE. EPS*PSI9O
      PSI9O   = PSI9



      IF (KI*NO3I .GT. A13) THEN
         BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
         CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT13A= 0.5D0*(-BB-DD)
         ROOT13B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT13A) THEN
            ROOT13 = ROOT13A
         ELSE
            ROOT13 = ROOT13B
         ENDIF
         ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
         PSI13  = CHI13-ROOT13
      ENDIF
      PSCONV13 = ABS(PSI13-PSI13O) .LE. EPS*PSI13O
      PSI13O   = PSI13



      IF (KI*CLI .GT. A14) THEN
         BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
         CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT14A= 0.5D0*(-BB-DD)
         ROOT14B= 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT14A) THEN
            ROOT14 = ROOT14A
         ELSE
            ROOT14 = ROOT14B
         ENDIF
         ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
         PSI14  = CHI14-ROOT14
      ENDIF
      PSCONV14 = ABS(PSI14-PSI14O) .LE. EPS*PSI14O
      PSI14O   = PSI14



      IF (NH4I*CLI .GT. A5) THEN
         BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
         CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT5A = 0.5D0*(-BB-DD)
         ROOT5B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT5A) THEN
            ROOT5 = ROOT5A
         ELSE
            ROOT5 = ROOT5B
         ENDIF
         ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
         PSI5  = CHI5-ROOT5
      ENDIF
      PSCONV5 = ABS(PSI5-PSI5O) .LE. EPS*PSI5O
      PSI5O   = PSI5



      IF (NAI*CLI .GT. A7) THEN
         BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
         CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT7A = 0.5D0*(-BB-DD)
         ROOT7B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT7A) THEN
            ROOT7 = ROOT7A
         ELSE
            ROOT7 = ROOT7B
         ENDIF
         ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
         PSI7  = CHI7-ROOT7
      ENDIF
      PSCONV7 = ABS(PSI7-PSI7O) .LE. EPS*PSI7O
      PSI7O   = PSI7



      IF (NAI*NO3I .GT. A8) THEN
         BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
         CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
         DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
         ROOT8A = 0.5D0*(-BB-DD)
         ROOT8B = 0.5D0*(-BB+DD)
         IF (ZERO.LE.ROOT8A) THEN
            ROOT8 = ROOT8A
         ELSE
            ROOT8 = ROOT8B
         ENDIF
         ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
         PSI8  = CHI8-ROOT8
      ENDIF
      PSCONV8 = ABS(PSI8-PSI8O) .LE. EPS*PSI8O
      PSI8O   = PSI8



      KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
      SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
      NH4I   = MAX (WAER(3) - ROOT5, ZERO)
      NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
      CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
      CAI    = ZERO
      NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
      MGI    = WAER(8)



      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
             - 2.D0*CAI - KI - 2.D0*MGI
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



      IF (HI.GT.OHI) THEN















         CALL CALCHS42p1 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL

      OHI   = AKW/HI

      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF



      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI



      CALL CALCMR2p1



      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT2p1
      ELSE
         IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
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
      CNH4CL  = CHI5 - PSI5
      CNACL   = CHI7 - PSI7
      CNANO3  = CHI8 - PSI8
      CMGSO4  = ZERO
      CK2SO4  = CHI9 - PSI9
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = CHI13 - PSI13
      KCL     = CHI14 - PSI14
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      RETURN



      END























      SUBROUTINE CALCW12p1
      INCLUDE 'module_isrpia_inc.F'
      EXTERNAL CALCW1A2p1, CALCW2A2p1



      IF (RH.LT.DRMP1) THEN
         SCASE = 'W1 ; SUBCASE 1'
         CALL CALCW1A2p1              
         SCASE = 'W1 ; SUBCASE 1'
      ELSE
         SCASE = 'W1 ; SUBCASE 2'  
         CALL CALCMDRPII2p1 (RH, DRMP1, DRCACL2, CALCW1A2p1, CALCW2A2p1)
         SCASE = 'W1 ; SUBCASE 2'
      ENDIF

      RETURN



      END



















      SUBROUTINE CALCW1A2p1
      INCLUDE 'module_isrpia_inc.F'



      CCASO4  = MIN (WAER(2), WAER(6))              
      CAFR    = MAX (WAER(6) - CCASO4, ZERO)
      SO4FR   = MAX (WAER(2) - CCASO4, ZERO)
      CK2SO4  = MIN (SO4FR, 0.5D0*WAER(7))          
      FRK     = MAX (WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
      CMGSO4  = SO4FR                               
      FRMG    = MAX (WAER(8) - CMGSO4, ZERO)
      CNACL   = MIN (WAER(1), WAER(5))              
      FRNA    = MAX (WAER(1) - CNACL, ZERO)
      CLFR    = MAX (WAER(5) - CNACL, ZERO)
      CCACL2  = MIN (CAFR, 0.5D0*CLFR)              
      CAFR    = MAX (CAFR - CCACL2, ZERO)
      CLFR    = MAX (WAER(5) - 2.D0*CCACL2, ZERO)
      CCANO32 = MIN (CAFR, 0.5D0*WAER(4))           
      CAFR    = MAX (CAFR - CCANO32, ZERO)
      FRNO3   = MAX (WAER(4) - 2.D0*CCANO32, ZERO)
      CMGCL2  = MIN (FRMG, 0.5D0*CLFR)              
      FRMG    = MAX (FRMG - CMGCL2, ZERO)
      CLFR    = MAX (CLFR - 2.D0*CMGCL2, ZERO)
      CMGNO32 = MIN (FRMG, 0.5D0*FRNO3)             
      FRMG    = MAX (FRMG - CMGNO32, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CMGNO32, ZERO)
      CNANO3  = MIN (FRNA, FRNO3)                   
      FRNA    = MAX (FRNA - CNANO3, ZERO)
      FRNO3   = MAX (FRNO3 - CNANO3, ZERO)
      CKCL    = MIN (FRK, CLFR)                     
      FRK     = MAX (FRK - CKCL, ZERO)
      CLFR    = MAX (CLFR - CKCL, ZERO)
      CKNO3   = MIN (FRK, FRNO3)                    
      FRK     = MAX (FRK - CKNO3, ZERO)
      FRNO3   = MAX (FRNO3 - CKNO3, ZERO)



      WATER   = ZERO

      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      RETURN



      END


