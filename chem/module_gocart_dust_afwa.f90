MODULE GOCART_DUST_AFWA





  USE module_data_gocart_dust

  IMPLICIT NONE  

  INTRINSIC max, min

CONTAINS
  SUBROUTINE gocart_dust_afwa_driver(ktau,dt,config_flags,julday,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,erod_dri,dustin,snowh,zs,   &
         ivgtyp,isltyp,vegfra,lai_vegmask,xland,xlat,xlong,gsw,dx,g,emis_dust,      &
         ust,znt,clay_wrf,sand_wrf,clay_nga,sand_nga,afwa_dustloft,                 &
         tot_dust,tot_edust,vis_dust,alpha,gamma,smtune,ustune,            &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_configure
  USE module_state_description
  USE module_data_sorgam, ONLY: factnuma,factnumc,soilfac

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: julday, ktau,                            &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                              moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                           chem
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL,           &
         INTENT(INOUT ) ::                                                 &
         emis_dust
   REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,     &
         INTENT(IN   ) ::                            smois
   REAL, DIMENSION( config_flags%num_soil_layers ) ,                       &
         INTENT(IN   ) ::                            zs
   REAL, DIMENSION( ims:ime , jms:jme, ndcls )             ,               &
         INTENT(IN   ) ::                            erod,erod_dri
   REAL, DIMENSION( ims:ime , jms:jme, 5 )                 ,               &
         INTENT(INOUT) ::                            dustin
   REAL, DIMENSION( ims:ime , jms:jme )                    ,               &
         INTENT(IN   ) ::                                                  &
                                                     u10,                  &
                                                     v10,                  &
                                                     gsw,                  &
                                                     vegfra,               &
                                                     lai_vegmask,          &
                                                     xland,                &
                                                     xlat,                 &
                                                     xlong,                &
                                                     ust,                  &
                                                     znt,                  &
                                                     clay_wrf,             &
                                                     sand_wrf,             &
                                                     clay_nga,             &
                                                     sand_nga,             &
                                                     snowh
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                         &
         INTENT(IN   ) ::                                                  &
                                                     alt,                  &
                                                     t_phy,                &
                                                     dz8w,p8w,             &
                                                     u_phy,v_phy,rho_phy
  REAL, DIMENSION( ims:ime , jms:jme )                    ,                &
         INTENT(  OUT) ::                            afwa_dustloft,        &
                                                     tot_edust
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                          &
         INTENT(  OUT) ::                            tot_dust,             &
                                                     vis_dust
  REAL, INTENT(IN   ) :: dt,dx,g

  

  INTEGER :: nmx,smx,i,j,k,imx,jmx,lmx,lhave
  INTEGER,DIMENSION (1,1) :: ilwi
  REAL*8, DIMENSION (1,1) :: erodtot
  REAL*8, DIMENSION (1,1) :: vegmask
  REAL*8, DIMENSION (1,1) :: gravsm
  REAL, DIMENSION( ims:ime , jms:jme ) :: clay,sand
  REAL*8, DIMENSION (1,1) :: drylimit
  REAL*8, DIMENSION (5)   :: tc,bems
  REAL*8, DIMENSION (1,1) :: airden,airmas,ustar
  REAL*8, DIMENSION (1) :: dxy
  REAL*8, DIMENSION (3) :: massfrac
  REAL*8                :: volsm
  REAL :: conver,converi
  REAL :: psi,ustart,w10
  REAL*8 :: zwant
  REAL, INTENT(IN   ) :: alpha, gamma, smtune, ustune
  INTEGER :: smois_opt

  conver=1.e-9
  converi=1.e9

  

  imx=1
  jmx=1
  lmx=1
  nmx=ndust
  smx=ngsalt

  k=kts
  DO j=jts,jte
  DO i=its,ite

    

    afwa_dustloft(i,j)=-99.



    IF (xland(i,j) .lt. 1.5) THEN
      ilwi(1,1)=1

      
      

      IF(config_flags%chem_opt == CB05_SORG_AQ_KPP .or. &
         config_flags%chem_opt == CB05_SORG_VBS_AQ_KPP ) then
        tc(1)=0.0
        tc(2)=0.0
        tc(3)=0.0
        tc(4)=0.0
        tc(5)=0.0
      ELSE
        tc(1)=chem(i,kts,j,p_dust_1)*conver
        tc(2)=chem(i,kts,j,p_dust_2)*conver
        tc(3)=chem(i,kts,j,p_dust_3)*conver
        tc(4)=chem(i,kts,j,p_dust_4)*conver
        tc(5)=chem(i,kts,j,p_dust_5)*conver
      END IF

      

      airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g
      airden(1,1)=rho_phy(i,kts,j)
      ustar(1,1)=ust(i,j)
      dxy(1)=dx*dx

      
      
      
      
      ustar(1,1)=ustar(1,1) * ustune
 

      
      
      
      
      

      IF (config_flags%dust_dsr .eq. 1) then  
          IF (erod_dri(i,j,1).ge.0) THEN  
             erodtot(1,1) = SUM(erod_dri(i,j,:))
          ELSE 
             erodtot(1,1)=SUM(erod(i,j,:))
          ENDIF
      ELSE 
          erodtot(1,1)=SUM(erod(i,j,:))
      ENDIF

      

      IF (config_flags%dust_veg .eq. 1) then

         
         
         

         IF (vegfra(i,j) .ge. 5.) then
             vegmask(1,1) = 0.0
         ELSE
             vegmask(1,1) = 1.0
         ENDIF
      ELSE IF (config_flags%dust_veg .eq. 2) then

         
         

         vegmask(1,1) = lai_vegmask(i,j)
      ELSE

         

         IF (erod(i,j,1) .eq. 0) THEN
            vegmask(1,1) = 0.0
         ELSE
            vegmask(1,1) = 1.0
         ENDIF
      ENDIF

      

      erodtot(1,1) = erodtot(1,1) * vegmask(1,1)

      
      
      
      
      
      
      

      IF (config_flags%dust_soils .eq. 1) then
         IF (clay_nga(i,j) .ge.0) then
             clay(i,j) = clay_nga(i,j)
             sand(i,j) = sand_nga(i,j)
         ELSE
             clay(i,j) =clay_wrf(i,j)
             sand(i,j) =sand_wrf(i,j)
         ENDIF
      ELSE
         clay(i,j) =clay_wrf(i,j)
         sand(i,j) =sand_wrf(i,j)
      ENDIF

      

      massfrac(1)=clay(i,j)
      massfrac(2)=1-(clay(i,j)+sand(i,j))
      massfrac(3)=sand(i,j)


      
      
      
      
      

      IF (znt(i,j) .gt. 0.2) then
        ilwi(1,1)=0
      ENDIF

      

      IF (isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
          isltyp(i,j) .eq. 18) then
        ilwi(1,1)=0
      ENDIF

      


      IF (snowh(i,j) .gt. 0.01) then 
        ilwi(1,1)=0
      ENDIF

      
      

      volsm=max(smois(i,1,j)*smtune,0.)

      

      gravsm(1,1)=100*volsm/((1.-porosity(isltyp(i,j)))*(2.65*(1-clay(i,j))+2.50*clay(i,j)))

      
      
      
      

      smois_opt = 0
      IF (config_flags%dust_smois == 1) then
         sfc_select: SELECT CASE(config_flags%sf_surface_physics)
            CASE ( RUCLSMSCHEME, PXLSMSCHEME )
               drylimit(1,1) =0.035*(13.52*clay(i,j)+3.53)**2.68
               smois_opt = 1
            CASE ( LSMSCHEME )
               drylimit(1,1) =0.0756*(15.127*clay(i,j)+3.09)**2.3211
               smois_opt = 1
            CASE DEFAULT


               

               drylimit(1,1)=14.0*clay(i,j)*clay(i,j)+17.0*clay(i,j)
         END SELECT sfc_select
      ELSE

         

         drylimit(1,1)=14.0*clay(i,j)*clay(i,j)+17.0*clay(i,j)
      END IF
 
      

      call source_dust(imx, jmx, lmx, nmx, smx, dt, tc, ustar, massfrac, &
                       erodtot, ilwi, dxy, gravsm, volsm, airden, airmas, &
                       bems, ustart, g, drylimit, alpha, gamma, smois_opt)

      IF(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
        dustin(i,j,1:5)=tc(1:5)*converi
      ELSE IF(config_flags%chem_opt == CB05_SORG_AQ_KPP .or.  &
              config_flags%chem_opt == CB05_SORG_VBS_AQ_KPP) then









         chem(i,kts,j,p_p25j)=chem(i,kts,j,p_p25j) + 0.03*sum(tc(1:5))*converi
         chem(i,kts,j,p_soila)=chem(i,kts,j,p_soila) + 0.97*1.02*sum(tc(1:5))*converi

         chem(i,kts,j,p_ac0) =  chem(i,kts,j,p_ac0) +  0.03*sum(tc(1:5))*converi*factnuma*soilfac
         chem(i,kts,j,p_corn) =  chem(i,kts,j,p_corn) + 0.97*1.02*sum(tc(1:5))*converi*factnumc*soilfac
      ELSE
        chem(i,kts,j,p_dust_1)=tc(1)*converi
        chem(i,kts,j,p_dust_2)=tc(2)*converi
        chem(i,kts,j,p_dust_3)=tc(3)*converi
        chem(i,kts,j,p_dust_4)=tc(4)*converi
        chem(i,kts,j,p_dust_5)=tc(5)*converi
      ENDIF

      

      psi=0.
      w10=(u10(i,j)**2.+v10(i,j)**2.)**0.5
      IF (ustar(1,1) .ne. 0. .and. znt(i,j) .ne. 0.) THEN
        psi=0.4*w10/ustar(1,1)-LOG(10.0/znt(i,j))
      ENDIF
      IF (erodtot(1,1) .gt. 0.) then
        afwa_dustloft(i,j)=ustune*w10-ustart*(LOG(10.0/znt(i,j))+psi)/0.4
      ENDIF

      

      emis_dust(i,1,j,p_edust1)=bems(1)
      emis_dust(i,1,j,p_edust2)=bems(2)
      emis_dust(i,1,j,p_edust3)=bems(3)
      emis_dust(i,1,j,p_edust4)=bems(4)
      emis_dust(i,1,j,p_edust5)=bems(5)

      

      tot_edust(i,j)=(bems(1)+bems(2)+bems(3)+bems(4)+bems(5))

    ENDIF

    
    
    
    

    DO k=kts,kte
      tot_dust(i,k,j)=(chem(i,k,j,p_dust_1)+chem(i,k,j,p_dust_2)+       &
                       chem(i,k,j,p_dust_3)+chem(i,k,j,p_dust_4)+       &
                       chem(i,k,j,p_dust_5))*rho_phy(i,k,j)
      IF ( tot_dust(i,k,j) .gt. 0. ) THEN
        vis_dust(i,k,j)=MIN(3.912/                                      &
                          ((1.470E-6*chem(i,k,j,p_dust_1)+              &
                            7.877E-7*chem(i,k,j,p_dust_2)+              &
                            4.623E-7*chem(i,k,j,p_dust_3)+              &
                            2.429E-7*chem(i,k,j,p_dust_4)+              &
                            1.387E-7*chem(i,k,j,p_dust_5))*             &
                           rho_phy(i,k,j)),999999.)
      ELSE
        vis_dust(i,k,j)=999999.
      ENDIF
    ENDDO

  ENDDO
  ENDDO

  END SUBROUTINE gocart_dust_afwa_driver

  SUBROUTINE source_dust(imx, jmx, lmx, nmx, smx, dt1, tc, ustar, massfrac,&
                         erod, ilwi, dxy, gravsm, volsm, airden, airmas, &
                         bems, ustart, g0, drylimit, alpha, gamma, smois_opt)

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  INTEGER, INTENT(IN)   :: nmx,imx,jmx,lmx,smx
  INTEGER, INTENT(IN)   :: ilwi(imx,jmx)
  REAL*8, INTENT(IN)    :: erod(imx,jmx)
  REAL*8, INTENT(IN)    :: ustar(imx,jmx)
  REAL*8, INTENT(IN)    :: gravsm(imx,jmx)
  REAL*8, INTENT(IN)    :: drylimit(imx,jmx) 
  REAL*8, INTENT(IN)    :: dxy(jmx)
  REAL*8, INTENT(IN)    :: airden(imx,jmx,lmx), airmas(imx,jmx,lmx)
  REAL*8, INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8, INTENT(OUT)   :: bems(imx,jmx,nmx) 
  REAL, INTENT(IN)    :: g0,dt1
  REAL, INTENT(OUT)   :: ustart
  INTEGER, INTENT(IN) :: smois_opt
  REAL*8, INTENT(IN)    :: volsm

  REAL*8    :: den(smx), diam(smx)
  REAL*8    :: dvol(nmx), dlndp(nmx)

  REAL*8    :: dsurface(smx), ds_rel(smx)
  REAL*8    :: massfrac(3)
  REAL*8    :: u_ts0, u_ts, dsrc, srce, dmass, dvol_tot
  REAL*8    :: emit, emit_vol
  REAL      :: rhoa, g
  REAL*8    :: salt, stotal
  INTEGER   :: i, j, m, n, s




  REAL, INTENT(IN)  :: alpha
  REAL, PARAMETER :: betamax=5.25E-4
  REAL*8 :: beta






  
  REAL, INTENT(IN) :: gamma





  REAL, PARAMETER :: cmb=1.0    

















  DO n=1,smx
    dmass=massfrac(spoint(n))*frac_salt(n)
    dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))  
  ENDDO
  







  stotal=SUM(dsurface(:))
  DO n=1,smx
    ds_rel(n)=dsurface(n)/stotal
  ENDDO







  g = g0*1.0E2                          
  emit=0.0

  DO n = 1, smx                         
    den(n) = den_salt(n)*1.0D-3         
    diam(n) = 2.0*reff_salt(n)*1.0D2    
    DO i = 1,imx
      DO j = 1,jmx
        rhoa = airden(i,j,1)*1.0D-3     
 
        
        

        u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 

        
        
        
        
        
        

        
        
        
        

        IF (smois_opt .EQ. 1) THEN
          IF (100.*volsm > drylimit(i,j)) THEN
            u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*((100.*volsm)-drylimit(i,j))**0.68)))
          ELSE
            u_ts = u_ts0
          ENDIF
        ELSE
          IF (gravsm(i,j) > drylimit(i,j)) THEN
            u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*(gravsm(i,j)-drylimit(i,j))**0.68)))
          ELSE
            u_ts = u_ts0
          END IF 
        END IF

        
        

        IF (n .eq. 7) THEN  
          ustart = u_ts
        ENDIF
 
        
        
 
        IF (ustar(i,j) .gt. u_ts .and. erod(i,j) .gt. 0.0 .and. ilwi(i,j) == 1) THEN
          salt = cmb*ds_rel(n)*(airden(i,j,1)/g0)*(ustar(i,j)**3)* &
                 (1. + u_ts/ustar(i,j))*(1. - (u_ts**2)/(ustar(i,j)**2))
        ELSE 
          salt = 0.D0
        ENDIF
 
        
        
        
        
        
 
        beta=10**(13.6*massfrac(1)-6.0)  
        IF (beta .gt. betamax) THEN
          beta=betamax
        ENDIF
        emit=emit+salt*(erod(i,j)**gamma)*alpha*beta    
      END DO
    END DO
  END DO                                 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  













 
  
  

  DO n=1,nmx  
    DO i=1,imx
      DO j=1,jmx

        

        dsrc = emit*distr_dust(n)*dxy(j)*dt1  
        IF (dsrc < 0.0) dsrc = 0.0

        

        tc(i,j,1,n) = tc(i,j,1,n) + dsrc / airmas(i,j,1) 
        bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) 
      END DO
    END DO
  END DO

  END SUBROUTINE source_dust


END MODULE GOCART_DUST_AFWA
