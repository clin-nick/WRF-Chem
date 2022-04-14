MODULE module_bioemi_beis314



























    CONTAINS
      SUBROUTINE bio_emissions_beis314(id,config_flags,ktau,curr_secs, &
               dtstep,julday,gmt,xlat,xlong,t_phy,p_phy,gsw,           &
               sebio_iso,sebio_oli,sebio_api,sebio_lim,sebio_xyl,      &
               sebio_hc3,sebio_ete,sebio_olt,sebio_ket,sebio_ald,      &
               sebio_hcho,sebio_eth,sebio_ora2,sebio_co,sebio_nr,      &
               sebio_sesq,sebio_mbo,                                   & 
               noag_grow,noag_nongrow,nononag,slai,                    &
               ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,           &
               ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,           &
               ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,   &
               ebio_sesq,ebio_mbo,                                     & 
               ids,ide, jds,jde, kds,kde,                              &
               ims,ime, jms,jme, kms,kme,                              &
               its,ite, jts,jte, kts,kte                               )

  USE module_configure
  USE module_state_description

      IMPLICIT NONE


      TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags


      INTEGER,   INTENT(IN   ) :: id,                                  &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte

      INTEGER, INTENT (IN)  ::    ktau,                                & 
                                  julday                                 

      REAL(KIND=8), INTENT(IN) :: curr_secs                              

      REAL, INTENT (IN)   ::      gmt,dtstep

      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
          INTENT(IN   ) ::                                             &
                                  t_phy,                               & 
                                  p_phy                                  

      REAL,  DIMENSION( ims:ime , jms:jme ),                           &
          INTENT(IN   ) ::                                             &
                                  xlat,                                & 
                                  xlong,                               & 
                                  gsw                                    


      REAL,  DIMENSION( ims:ime , jms:jme ),                           &
          INTENT(IN   ) ::                                             &
               sebio_iso,sebio_oli,sebio_api,sebio_lim,sebio_xyl,      &
               sebio_hc3,sebio_ete,sebio_olt,sebio_ket,sebio_ald,      &
               sebio_hcho,sebio_eth,sebio_ora2,sebio_co,sebio_nr,      &
               sebio_sesq,sebio_mbo,                                   &
               noag_grow,noag_nongrow,nononag


      REAL,  DIMENSION( ims:ime , jms:jme ),                           &
          INTENT(IN   ) ::        slai 


      REAL,  DIMENSION( ims:ime , jms:jme ),                           &
          INTENT(INOUT  ) ::                                           &
               ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,           &
               ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,           &
               ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,   &
               ebio_sesq,ebio_mbo   



      INTEGER :: i,j



      REAL  ::  tair      
      REAL  ::  tsolar    
      REAL  ::  pres      
      REAL  ::  ylat      
      REAL  ::  ylong     

      REAL  :: se_iso,se_oli,se_api,se_lim,se_xyl,      &
               se_hc3,se_ete,se_olt,se_ket,se_ald,      &
               se_hcho,se_eth,se_ora2,se_co,se_nr,      &
               se_mbo,se_sesq,                          & 
               growagno,ngrowagno,nonagno

      REAL  ::  tlai  

      REAL  :: e_no



      REAL  ::  ct, dt       
      REAL  ::  cfno         
      REAL  ::  cfovoc       
      REAL  ::  par          
      REAL  ::  csubl        
      REAL  ::  zen          
      REAL  ::  coszen       
      REAL  ::  pardb        
      REAL  ::  pardif       
      REAL :: gmtp           


      INTEGER , PARAMETER ::  ldev = 6    
      CHARACTER*256   ::   mesg
















                         



      gmtp = curr_secs/3600._8

      gmtp=mod(gmt+gmtp,24.)
      write(mesg,*) 'calculate beis314 emissions at gmtp = ',gmtp
      call wrf_debug(15,mesg)
      DO 100 j=jts,jte  
      DO 100 i=its,ite  

           tair = t_phy(i,kts,j)
           pres = .01*p_phy(i,kts,j)
           ylat = xlat(i,j)
           ylong = xlong(i,j)
           tsolar = gsw(i,j)
           tlai = slai(i,j)
           se_iso  = sebio_iso(i,j)
           se_oli  = sebio_oli(i,j)
           se_api  = sebio_api(i,j)
           se_lim  = sebio_lim(i,j)
           se_xyl  = sebio_xyl(i,j)
           se_hc3  = sebio_hc3(i,j)
           se_ete  = sebio_ete(i,j)
           se_olt  = sebio_olt(i,j)
           se_ket  = sebio_ket(i,j)
           se_ald  = sebio_ald(i,j)
           se_hcho  = sebio_hcho(i,j)
           se_eth  = sebio_eth(i,j)
           se_ora2  = sebio_ora2(i,j)
           se_co  = sebio_co(i,j)
           se_nr  = sebio_nr(i,j)

           
           se_sesq  = sebio_sesq(i,j)
           se_mbo  = sebio_mbo(i,j)

           growagno = noag_grow(i,j)
           ngrowagno = noag_nongrow(i,j) 
           nonagno = nononag(i,j)



           IF (tair .LT. 200.0) THEN



               WRITE( ldev, * ) mesg
           END IF

           IF (tair .GT. 315.0 ) THEN




               tair = 315.0

           ENDIF



           dt   = 28668.514 / tair
           ct   = EXP( 37.711 - 0.398570815 * dt ) /    &
                         (1.0 + EXP( 91.301 - dt ) )



           CALL calc_zenithb(ylat,-ylong,julday,gmtp,zen)
           coszen = COS(zen)


           CALL getpar( tsolar, pres, zen, pardb, pardif )
           par = pardb + pardif



           IF ( par .LT. 0.00 .OR. par .GT. 2600.0 ) THEN




           ENDIF


           IF ( tlai .GT. 10.0 ) THEN




           ENDIF


           csubl = 0.0


           IF ( par .LE. 0.01 .OR. coszen .LE. 0.02079483 ) THEN
               ebio_iso(i,j) = 0.0

           ELSE


               IF ( tlai .GT. 0.1 ) THEN
                     csubl = clnew( zen, pardb, pardif, tlai )


               ELSE  
                     csubl  = cguen( par )

               ENDIF

               ebio_iso(i,j) = se_iso * ct * csubl

           ENDIF





           cfovoc = EXP( 0.09 * ( tair - 303.0 ) )

           ebio_oli(i,j)   =  se_oli * cfovoc
           ebio_api(i,j)   =  se_api * cfovoc
           ebio_lim(i,j)   =  se_lim * cfovoc
           ebio_xyl(i,j)   =  se_xyl * cfovoc
           ebio_hc3(i,j)   =  se_hc3 * cfovoc
           ebio_ete(i,j)   =  se_ete * cfovoc
           ebio_olt(i,j)   =  se_olt * cfovoc
           ebio_ket(i,j)   =  se_ket * cfovoc
           ebio_ald(i,j)   =  se_ald * cfovoc
           ebio_hcho(i,j)  =  se_hcho * cfovoc
           ebio_eth(i,j)   =  se_eth * cfovoc
           ebio_ora2(i,j)  =  se_ora2 * cfovoc
           ebio_co(i,j)    =  se_co * cfovoc
           ebio_nr(i,j)    =  se_nr * cfovoc

           
           ebio_sesq(i,j)  =  se_sesq * cfovoc
           ebio_mbo(i,j)   =  se_mbo * cfovoc



           CALL hrno( julday, growagno, ngrowagno, nonagno, tair, e_no)

           ebio_no(i,j) = e_no

 100  CONTINUE

      RETURN










94010   FORMAT( A, F10.2, 1X, A, I4, A1, I4)
94020   FORMAT( A, F10.2, 1X, A, I4, A1, I4, 1X, A)



            CONTAINS

            REAL FUNCTION clnew( zen, pardb, pardif, tlai )






            IMPLICIT NONE

            REAL, INTENT (IN) :: pardb    
            REAL, INTENT (IN) :: pardif   
            REAL, INTENT (IN) :: zen      
            REAL, INTENT (IN) :: tlai     

            REAL kbe                
            REAL ALPHA              
            REAL KD                 
            REAL SQALPHA            
            REAL canparscat         
            REAL canpardif          
            REAL parshade           
            REAL parsun             
            REAL laisun             
            REAL fracsun            
            REAL fracshade          

            ALPHA = 0.8
            SQALPHA = SQRT(0.8)
            KD = 0.68



            kbe = 0.5 * SQRT(1. + TAN( zen ) * TAN( zen ))
 


            canparscat = 0.5 * pardb * (EXP(-1.*sqalpha*kbe*tlai) -    &
                     EXP(-1.* kbe * tlai))



            canpardif  = pardif * (1.-EXP(-1.*sqalpha*kd*tlai)) /      &
                     (sqalpha*kd*tlai)



            parshade   = canpardif + canparscat
	    parsun     = kbe * pardb  + parshade
	    laisun     = (1. - EXP( -1. * kbe * tlai))/kbe
	    fracsun    = laisun/tlai
	    fracshade  = 1. - fracsun







            clnew =fracsun * cguen(parsun) + fracshade * cguen(parshade)

            RETURN 
            END FUNCTION clnew

            REAL FUNCTION cguen( partmp ) 






            IMPLICIT NONE
            REAL, INTENT (IN) :: partmp
            REAL, PARAMETER :: alpha = 0.001
            REAL, PARAMETER :: cl = 1.42

            IF ( partmp .LE. 0.01) THEN
               cguen = 0.0
            ELSE
               cguen = (alpha *cl * partmp) /                         &
                   SQRT(1. + alpha * alpha * partmp * partmp)
            ENDIF
 
            RETURN
            END FUNCTION cguen

      END SUBROUTINE bio_emissions_beis314



      SUBROUTINE calc_zenithb(lat,long,ijd,gmt,zenith)
        
        
        
        
        
        
        
        
        
        
        
        
        
        CHARACTER*256   ::   mesg
        REAL :: gmt, lat, long, zenith
        INTEGER :: ijd
        
        REAL :: csz, cw, d, decl, dr, ec, epsi, eqt, eyt, feqt, feqt1, &
          feqt2, feqt3, feqt4, feqt5, feqt6, feqt7, lbgmt, lzgmt, ml, pepsi, &
          pi, ra, rdecl, reqt, rlt, rml, rphi, rra, ssw, sw, tab, w, wr, &
          yt, zpt, zr
        INTEGER :: jd
        
        INTRINSIC acos, atan, cos, min, sin, tan
        
        pi = 3.1415926535590
        dr = pi/180.
        rlt = lat*dr
        rphi = long*dr

        

        jd = ijd

        d = jd + gmt/24.0
        
        ml = 279.2801988 + .9856473354*d + 2.267E-13*d*d
        rml = ml*dr

        
        
        
        
        w = 282.4932328 + 4.70684E-5*d + 3.39E-13*d*d
        wr = w*dr
        ec = 1.6720041E-2 - 1.1444E-9*d - 9.4E-17*d*d
        epsi = 23.44266511 - 3.5626E-7*d - 1.23E-15*d*d
        pepsi = epsi*dr
        yt = (tan(pepsi/2.0))**2
        cw = cos(wr)
        sw = sin(wr)
        ssw = sin(2.0*wr)
        eyt = 2.*ec*yt
        feqt1 = sin(rml)*(-eyt*cw-2.*ec*cw)
        feqt2 = cos(rml)*(2.*ec*sw-eyt*sw)
        feqt3 = sin(2.*rml)*(yt-(5.*ec**2/4.)*(cw**2-sw**2))
        feqt4 = cos(2.*rml)*(5.*ec**2*ssw/4.)
        feqt5 = sin(3.*rml)*(eyt*cw)
        feqt6 = cos(3.*rml)*(-eyt*sw)
        feqt7 = -sin(4.*rml)*(.5*yt**2)
        feqt = feqt1 + feqt2 + feqt3 + feqt4 + feqt5 + feqt6 + feqt7
        eqt = feqt*13751.0

        
        reqt = eqt/240.
        
        ra = ml - reqt
        rra = ra*dr
        
        tab = 0.43360*sin(rra)
        rdecl = atan(tab)
        decl = rdecl/dr
        
        lbgmt = 12.0 - eqt/3600. + long*24./360.
        lzgmt = 15.0*(gmt-lbgmt)
        zpt = lzgmt*dr
        csz = sin(rlt)*sin(rdecl) + cos(rlt)*cos(rdecl)*cos(zpt)
        if(csz.gt.1) then
           write(mesg,*) 'calczen,csz ',csz
           call wrf_debug(15,mesg)
        endif
        csz = min(1.,csz)
        zr = acos(csz)


        zenith = zr 

        RETURN

      END SUBROUTINE calc_zenithb




        SUBROUTINE getpar( tsolar, pres, zen, pardb, pardif )








































      IMPLICIT NONE



        REAL , INTENT  (IN) :: tsolar   
        REAL , INTENT  (IN) :: pres     
        REAL , INTENT  (IN) :: zen      
 


        REAL, INTENT (OUT) :: pardb     
        REAL, INTENT (OUT) :: pardif    



        REAL, PARAMETER :: watt2umol = 4.6  


        REAL ratio		
        REAL ot                 
        REAL rdvis              
        REAL rfvis              
        REAL wa                 
        REAL rdir               
        REAL rfir               
        REAL rvt                
        REAL rirt               
        REAL fvis               
        REAL fvb                
        REAL fvd                







        IF (zen .GE. 1.51844 .OR. tsolar .LE. 0.) THEN
             pardb  = 0.
             pardif = 0.
             RETURN
        ENDIF
	   


        ot    = pres / 1013.25 / COS(zen)              
        rdvis = 600. * EXP(-.185 * ot) * COS(zen)      
        rfvis = 0.42 * (600 - rdvis) * COS(zen)        
        wa    = 1320 * .077 * (2. * ot)**0.3           
        rdir  = (720. * EXP(-0.06 * ot)-wa) * COS(zen) 
        rfir  = 0.65 * (720. - wa - rdir) * COS(zen)   

        rvt   = rdvis + rfvis                    
        rirt  = rdir + rfir                      
        fvis  = rvt/(rirt + rvt)                 
        ratio = tsolar /(rirt + rvt)             



        IF (ratio .GE. 0.89) THEN
           fvb = rdvis/rvt * 0.941124
        ELSE IF (ratio .LE. 0.21) THEN
           fvb = rdvis/rvt * 9.55E-3
        ELSE
           fvb = rdvis/rvt * (1.-((0.9 - ratio)/0.7)**0.666667)
        ENDIF
        fvd = 1. - fvb



        pardb  = tsolar * fvis * fvb * watt2umol	
        pardif = tsolar * fvis * fvd * watt2umol      


        RETURN 








        END SUBROUTINE getpar

      SUBROUTINE hrno( julday, growagno, ngrowagno, nonagno, tairin, e_no)


















































      IMPLICIT NONE



        INTEGER, INTENT (IN)  :: julday   


        REAL, INTENT (IN)  ::  tairin     
        REAL, INTENT (IN)  ::  growagno     
        REAL, INTENT (IN)  ::  ngrowagno    
        REAL, INTENT (IN)  ::  nonagno      

        REAL, INTENT (OUT)  ::  e_no      



        REAL            cfno         
        REAL            cfnograss    
        REAL            tsoi         
        REAL            tair         

        REAL          :: cfnowet, cfnodry

        INTEGER gday


        tair = tairin



        gday = growseason(julday)

        
        IF (gday .eq. 0) THEN         
           IF ( tair .GT. 303.00 ) tair = 303.00

           IF ( tair .GT. 268.8690 ) THEN  
              cfno = EXP( 0.04686 * tair  -  14.30579 ) 
           ELSE
              cfno = 0.0
           ENDIF

           e_no =                   &
                 ngrowagno * cfno   &     
                 +  nonagno * cfno   

        ELSE 

           tsoi = 0.72*tair+82.28
           IF (tsoi .LE. 273.16) tsoi = 273.16
           IF (tsoi .GE. 303.16) tsoi = 303.16

           cfnodry = (1./3.)*(1./30.)*(tsoi-273.16)  
           IF (tsoi .LE. 283.16) THEN       
              cfnowet = (tsoi-273.16)*EXP(-0.103*30.0)*0.28 
           ELSE                             
              cfnowet =  EXP(0.103*(tsoi-273.16))   &
                         *EXP(-0.103*30.0)
           ENDIF
           cfno = 0.5*cfnowet + 0.5*cfnodry

           IF ( tair .GT. 303.00 ) tair = 303.00

           IF ( tair .GT. 268.8690 ) THEN  
              cfnograss = EXP( 0.04686 * tair  -  14.30579 ) 
           ELSE
              cfnograss = 0.0
           ENDIF

           e_no =  growagno * cfno *fertilizer_adj(julday)*veg_adj(julday)   &
                  +  nonagno * cfnograss                   

        ENDIF

        RETURN


        CONTAINS

        REAL FUNCTION fertilizer_adj(julday)











        implicit none
        integer julday



       integer gday



     gday = growseason(julday)

     
      
      IF (gday .EQ. 0) THEN
          fertilizer_adj = 0.0
      ELSEIF ((gday .GE. 1) .AND. (gday .LT. 30)) THEN 
          fertilizer_adj = 1.0
      ELSEIF (gday .GE. 30)   THEN
          fertilizer_adj = 1.0+30.0/184.0-float(gday)/184.0
      ELSE
          write (*,*) 'ERROR: invalid Julian day'
	  write (*,*) 'julday = ', julday
	  write (*,*) 'growing season day = ',gday
	  CALL wrf_error_fatal3("<stdin>",793,&
'INVALID GROWING SEASON DAY')
      ENDIF
	
      RETURN

      END FUNCTION fertilizer_adj


      REAL FUNCTION veg_adj(julday)











      implicit none
  
       integer julday





      integer gday




      gday = growseason(julday)

      
      
      IF (gday .LE. 30) THEN
          veg_adj = 1.0
      ELSEIF ((gday .GT. 30) .AND. (gday .LT. 60)) THEN 
          veg_adj = 1.5-(float(gday)/60.0)
      ELSEIF (gday .GE. 60) THEN 
          veg_adj = 0.5
      ELSE
          write (*,*) 'ERROR: invalid Julian day'
	  write (*,*) 'julday = ', julday
	  write (*,*) 'growing season day = ',gday
	  CALL wrf_error_fatal3("<stdin>",841,&
'veg_adj: INVALID GROWING SEASON DAY' )
      ENDIF


      RETURN


      END FUNCTION veg_adj      

     END SUBROUTINE hrno 

      INTEGER FUNCTION growseason(julday)








      implicit none       
      integer julday










      integer gsjulian_start
      integer gsjulian_end

      data gsjulian_start /91/ 
      data gsjulian_end /304/ 
	 
      IF      ((julday .GE. gsjulian_start)       &
         .AND. (julday .LE. gsjulian_end)) THEN   
       
         growseason = julday-gsjulian_start+1

	  
      ELSEIF  ((julday .GE. 1)     &       
         .AND. (julday .LE. 366)) THEN      
     
         growseason = 0
	 
      ELSE
          write (*,*) 'ERROR: Invalid julday '
	  write (*,*) 'julday = ',julday
	  CALL wrf_error_fatal3("<stdin>",894,&
'growseason: INVALID JULIAN DAY')
      ENDIF


      RETURN
      END FUNCTION growseason


END MODULE module_bioemi_beis314
