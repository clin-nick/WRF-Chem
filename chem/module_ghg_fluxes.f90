MODULE module_ghg_fluxes
   USE module_configure
   USE module_state_description

   IMPLICIT NONE


































CONTAINS

SUBROUTINE add_ghg_fluxes (  ids,ide, jds,jde, kds,kde,             &     
                             ims,ime, jms,jme, kms,kme,             &     
                             its,ite, jts,jte, kts,kte,             &     
                                                                    
                             dtstep, dz8w, config_flags, rho_phy,   &
                             chem, emis_ant,eghg_bio,ebio_co2oce    )





TYPE(grid_config_rec_type), INTENT(IN   )    :: config_flags

INTEGER, INTENT(IN   ) ::   ids,ide, jds,jde, kds,kde,        &
                            ims,ime, jms,jme, kms,kme,        &
                            its,ite, jts,jte, kts,kte

REAL, INTENT(IN   )    ::   dtstep

REAL, DIMENSION( ims:ime,jms:jme ), INTENT(IN ) :: ebio_co2oce
REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT ) ::      chem
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) ::   rho_phy, dz8w
REAL, DIMENSION( ims:ime, kms:config_flags%kemit, jms:jme, num_emis_ant ), INTENT(IN ) :: emis_ant
REAL, DIMENSION( ims:ime, 1,jms:jme, num_eghg_bio ), INTENT(IN ) :: eghg_bio

INTEGER :: i,j,k
REAL    :: conv_rho

call wrf_debug(15,'add_ghg_fluxes')


DO j=jts,jte
   DO i=its,ite
      
      DO k=kts,min(config_flags%kemit,kte)
         conv_rho=8.0461e-6/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)  
         chem(i,k,j,p_co2_ant)= chem(i,k,j,p_co2_ant) + conv_rho* emis_ant(i,k,j,p_e_co2)
         chem(i,k,j,p_co2_tst)= chem(i,k,j,p_co2_tst) + conv_rho* emis_ant(i,k,j,p_e_co2tst)
         chem(i,k,j,p_co_ant) = chem(i,k,j,p_co_ant)  + conv_rho* emis_ant(i,k,j,p_e_co)

      
      if (k==1) then
         chem(i,1,j,p_co2_bio)= chem(i,1,j,p_co2_bio) + conv_rho* (eghg_bio(i,1,j,p_ebio_gee) + eghg_bio(i,1,j,p_ebio_res))   
         chem(i,1,j,p_co2_oce)= chem(i,1,j,p_co2_oce) + conv_rho* ebio_co2oce(i,j)
      end if

      ENDDO
   ENDDO
ENDDO


IF(config_flags%chem_opt==GHG_TRACER) THEN
    DO j=jts,jte
       DO i=its,ite
          
          DO k=kts,min(config_flags%kemit,kte)
             conv_rho=8.0461e-6/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)  
             chem(i,k,j,p_ch4_ant)= chem(i,k,j,p_ch4_ant) + conv_rho* emis_ant(i,k,j,p_e_ch4)
             chem(i,k,j,p_ch4_tst)= chem(i,k,j,p_ch4_tst) + conv_rho* emis_ant(i,k,j,p_e_ch4tst)
             chem(i,k,j,p_co_tst) = chem(i,k,j,p_co_tst)  + conv_rho* emis_ant(i,k,j,p_e_cotst)

          
          if (k==1) then
             chem(i,1,j,p_ch4_bio)= chem(i,1,j,p_ch4_bio) + conv_rho* (eghg_bio(i,1,j,p_ebio_ch4wet) + eghg_bio(i,1,j,p_ebio_ch4soil) &
                                                                       + eghg_bio(i,1,j,p_ebio_ch4term))
          end if

          ENDDO
       ENDDO
    ENDDO

END IF

END SUBROUTINE add_ghg_fluxes

SUBROUTINE VPRM                 (  ids,ide, jds,jde,                 &
                                   ims,ime, jms,jme,                 &
                                   its,ite, jts,jte,                 &

                                   vprm_in,rad0, lambda,            &
                                   alpha, RESP0,                    &
                                   T2,RAD, eghg_bio                  )


INTEGER,  INTENT(IN   )   ::  ids,ide, jds,jde, &
                              ims,ime, jms,jme, &
                              its,ite, jts,jte
INTEGER :: i,j,m

REAL, PARAMETER :: const= 3.6e+3 

REAL, DIMENSION (ims:ime, jms:jme), INTENT(IN)    ::  T2,  RAD

REAL, DIMENSION( ims:ime, 8, jms:jme, num_vprm_in), INTENT(IN)   ::  vprm_in
REAL, DIMENSION( ims:ime, 1,jms:jme, num_eghg_bio ), INTENT(INOUT ) ::  eghg_bio

REAL, DIMENSION(num_vprm_in) :: rad0, lambda, alpha, RESP0
REAL  ::  a1,a2,a3,Tair,Tscale,Wscale,Pscale,GEE_frac,RESP_frac,evithresh,RADscale










REAL, DIMENSION(8) ::  Tmax,Tmin,Topt

DATA Tmin   /0.,0.,0.,2.,2.,5.,2.,0./
DATA Topt   /20.,20.,20.,20.,20.,22.,18.,0./
DATA Tmax   /8*40./

DO j=jts,min(jte,jde-1)
DO i=its,min(ite,ide-1)

   GEE_frac= 0.
   RESP_frac= 0.

   Tair= T2(i,j)-273.15
   veg_frac_loop: DO m=1,7

        if (vprm_in(i,m,j,p_vegfra_vprm)<1.e-8) CYCLE  

        a1= Tair-Tmin(m)
        a2= Tair-Tmax(m)
        a3= Tair-Topt(m)


        if (a1<0. .OR. a2>0.) then
            Tscale= 0.
        else
            Tscale=a1*a2/(a1*a2 - a3**2)
        end if

        if (Tscale<0.) then
            Tscale=0.
        end if

       
        if (m==4 .OR. m==7) then  
            if (vprm_in(i,m,j,p_lswi_max)<1e-7) then  
                Wscale= 0.
            else
                Wscale= (vprm_in(i,m,j,p_lswi)-vprm_in(i,m,j,p_lswi_min))/(vprm_in(i,m,j,p_lswi_max)-vprm_in(i,m,j,p_lswi_min))
            end if
        else
            Wscale= (1.+vprm_in(i,m,j,p_lswi))/(1.+vprm_in(i,m,j,p_lswi_max))
        end if

        
        if (m==1) then  
            Pscale= 1.
        else if (m==5 .OR. m==7) then  
            Pscale= (1.+vprm_in(i,m,j,p_lswi))/2.
        else                           
            evithresh= vprm_in(i,m,j,p_evi_min) + 0.55*(vprm_in(i,m,j,p_evi_max)-vprm_in(i,m,j,p_evi_min))
            if (vprm_in(i,m,j,p_evi)>=evithresh) then  
               Pscale= 1.
            else
               Pscale=(1.+vprm_in(i,m,j,p_lswi))/2.  
            end if
        end if

        RADscale= 1./(1. + RAD(i,j)/rad0(m))
        GEE_frac= lambda(m)*Tscale*Pscale*Wscale*RADscale* vprm_in(i,m,j,p_evi)* RAD(i,j)*vprm_in(i,m,j,p_vegfra_vprm) + GEE_frac

        RESP_frac= (alpha(m)*Tair + RESP0(m))*vprm_in(i,m,j,p_vegfra_vprm) + RESP_frac


   ENDDO veg_frac_loop
   
   eghg_bio(i,1,j,p_ebio_gee)= min(0.0,const*GEE_frac)
   eghg_bio(i,1,j,p_ebio_res)= max(0.0,const*RESP_frac)

  ENDDO
  ENDDO

END SUBROUTINE VPRM






SUBROUTINE KAPLAN             ( ids,ide, jds,jde,                        &
                                ims,ime, jms,jme,                        &
                                its,ite, jts,jte,                        &

                                dt, T_soil, Soil_M, wet_in,              &
                                soil_type, T_skin, eghg_bio,             &
                                num_soil_layers,E_f,M_s                  )



INTEGER,  INTENT(IN   )   ::  ids,ide, jds,jde, &
                              ims,ime, jms,jme, &
                              its,ite, jts,jte, &
                              num_soil_layers
INTEGER :: i,j
INTEGER, DIMENSION (ims:ime, jms:jme), INTENT(IN)  :: soil_type

REAL, PARAMETER :: const= 3.6e9/(0.012+4.*0.001) 
REAL, DIMENSION (ims:ime,1,jms:jme, num_wet_in), INTENT(IN)    :: wet_in
REAL, DIMENSION (ims:ime, jms:jme), INTENT(IN)    :: T_skin  
REAL, DIMENSION (ims:ime, num_soil_layers, jms:jme), INTENT(IN) :: Soil_M, T_soil
REAL, DIMENSION (ims:ime, 1,jms:jme, num_eghg_bio), INTENT(INOUT) ::  eghg_bio
REAL, INTENT(IN) :: dt,M_s,E_f

REAL, DIMENSION( ims:ime, jms:jme) :: sm_01, sm_02, sm_03, sm_04, g_T, f_SM, &
                                      k_r, Carbon, HR, T_flood, sm_tot, T_2s, &
                                      T_peat, P_l, soil_sat

REAL  :: tau_10 
































  tau_10 = 2.86  

 DO j=jts,min(jte,jde-1)
  DO i=its,min(ite,ide-1)



    eghg_bio(i,1,j,p_ebio_ch4wet)=0.
    HR(i,j)=0.
    k_r(i,j)=0.
    g_T(i,j)=0.
    T_flood(i,j)=0.
    T_peat(i,j)=0.
    P_l(i,j)=0.
    f_SM(i,j)=0.
    sm_01(i,j)=0.
    sm_02(i,j)=0.
    sm_03(i,j)=0.
    sm_04(i,j)=0.    
    Carbon(i,j)=0.
    sm_tot(i,j)=0.
    T_2s(i,j)=0.



  IF(soil_type(i,j) == 1)then
   soil_sat(i,j) = 0.339
  ELSE IF(soil_type(i,j) == 2) then
   soil_sat(i,j) = 0.421
  ELSE IF(soil_type(i,j) == 3) then
   soil_sat(i,j) = 0.434
  ELSE IF(soil_type(i,j) == 4) then
   soil_sat(i,j) = 0.476
  ELSE IF(soil_type(i,j) == 5) then
   soil_sat(i,j) = 0.476
  ELSE IF(soil_type(i,j) == 6) then
   soil_sat(i,j) = 0.439
  ELSE IF(soil_type(i,j) == 7) then
   soil_sat(i,j) = 0.404
  ELSE IF(soil_type(i,j) == 8) then
   soil_sat(i,j) = 0.464
  ELSE IF(soil_type(i,j) == 9) then
   soil_sat(i,j) = 0.465
  ELSE IF(soil_type(i,j) == 10) then
   soil_sat(i,j) = 0.406
  ELSE IF(soil_type(i,j) == 11) then
   soil_sat(i,j) = 0.468
  ELSE IF(soil_type(i,j) == 12) then
   soil_sat(i,j) = 0.468
  ELSE IF(soil_type(i,j) == 13) then
   soil_sat(i,j) = 0.439
  ELSE IF(soil_type(i,j) == 14) then
   soil_sat(i,j) = 1.0
  ELSE IF(soil_type(i,j) == 15) then
   soil_sat(i,j) = 0.20
  ELSE IF(soil_type(i,j) == 16) then
   soil_sat(i,j) = 0.421
  ELSE IF(soil_type(i,j) == 17) then
   soil_sat(i,j) = 0.468
  ELSE IF(soil_type(i,j) == 18) then
   soil_sat(i,j) = 0.200
  ELSE IF(soil_type(i,j) == 19) then
   soil_sat(i,j) = 0.339
  ELSE
   soil_sat(i, j) = 0.4
  END IF   




  sm_01(i,j)  = Soil_M(i,1,j)
  sm_02(i,j)  = Soil_M(i,2,j)
  sm_03(i,j)  = Soil_M(i,3,j)
  sm_04(i,j)  = Soil_M(i,4,j)
 
  sm_tot(i,j) = 0.5*(sm_01(i,j)+sm_02(i,j))

  IF(sm_tot(i,j) < 0.0) then
     sm_tot(i,j) = 0.0
  END IF

  IF ((T_Soil(i,1,j) == 273.16).OR.(T_Soil(i,1,j) == 0.0)) then
     T_2s(i, j) = T_skin(i, j) - 273.15
  ELSE
     T_2s(i,j)= T_Soil(i,1,j)-273.15
  END IF



    g_T(i,j) = EXP(308.56*(1.0/56.02 - 1.0/(T_2s(i,j) + 46.02)))




    f_SM(i,j) = 0.25 + 0.75 * (sm_tot(i,j)/soil_sat(i,j))



    k_r(i,j) = ((1/tau_10)*g_T(i,j)*f_SM(i,j))/(12.0*24.0*30.0)



    Carbon(i,j) = wet_in(i,1,j,p_Cpool)*(1.0 - EXP(-k_r(i,j)*1.0))



    Carbon(i,j) = Carbon(i,j)/(3600.0)



    HR(i,j) =  0.7*Carbon(i,j)



  T_flood(i,j) = HR(i,j)*M_s*wet_in(i,1,j,p_wetmap)



  T_peat(i,j) = HR(i,j)*E_f*wet_in(i,1,j,p_wetmap)



   P_l(i,j) = EXP((wet_in(i,1,j,p_T_ann) - 303.0)/8.0)
   eghg_bio(i,1,j,p_ebio_ch4wet) = P_l(i,j)*T_flood(i,j)+(1.0 - P_l(i,j))*T_peat(i,j)

  
  eghg_bio(i,1,j,p_ebio_ch4wet) = eghg_bio(i,1,j,p_ebio_ch4wet)*const

  END DO

  END DO

END SUBROUTINE KAPLAN 

  SUBROUTINE TERMITE (                    ids,ide, jds,jde,                            &
                                          ims,ime, jms,jme,                            &
                                          its,ite, jts,jte,                            &

                                          dt,eghg_bio, vtype_ter,                      &
                                          biomt_par, emit_par                          )



INTEGER,  INTENT(IN   )   ::  ids,ide, jds,jde, &
                              ims,ime, jms,jme, &
                              its,ite, jts,jte
INTEGER :: i,j,vtype

INTEGER, DIMENSION (ims:ime, jms:jme), INTENT(IN)    ::  vtype_ter
REAL, DIMENSION (ims:ime, 1,jms:jme, num_eghg_bio), INTENT(INOUT) ::  eghg_bio
REAL, INTENT(IN) :: dt

REAL, PARAMETER :: const= 3.6e9/(0.012+4.*0.001) 
REAL, DIMENSION(14), INTENT(IN) :: biomt_par
REAL, DIMENSION(14), INTENT(IN) :: emit_par

CHARACTER*256 :: mminlu_loc

 DO j=jts,min(jte,jde-1)
  DO i=its,min(ite,ide-1)

   CALL nl_get_mminlu(1 , mminlu_loc)

    IF(mminlu_loc .EQ. 'USGS')THEN

     vtype = INT(vtype_ter(i,j))

     SELECT CASE (vtype)

       CASE (1)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (2)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(10)*emit_par(10)
       CASE (3)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(10)*emit_par(10)
       CASE (4)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(11)*emit_par(11)
       CASE (5)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(11)*emit_par(11)
       CASE (6)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(11)*emit_par(11)
       CASE (7)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(7)*emit_par(7)
       CASE (8)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(12)*emit_par(12)
       CASE (9)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(12)*emit_par(12)
       CASE (10)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(5)*emit_par(5)
       CASE (11)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(4)*emit_par(4)
       CASE (12)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(4)*emit_par(4)
       CASE (13)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(1)*emit_par(1)
       CASE (14)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(2)*emit_par(2)
       CASE (15)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(4)*emit_par(4)
       CASE (16)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (17)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (18)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (19)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(13)*emit_par(13)
       CASE (20)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (21)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (22)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (23)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (24)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE default
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       END SELECT


    ELSE IF(mminlu_loc .EQ. 'MODIFIED_IGBP_MODIS_NOAH')THEN

     vtype = INT(vtype_ter(i,j))

     SELECT CASE (vtype)

       CASE (1)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(2)*emit_par(2)
       CASE (2)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(1)*emit_par(1)
       CASE (3)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(4)*emit_par(4)
       CASE (4)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(4)*emit_par(4)
       CASE (5)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(4)*emit_par(4)
       CASE (6)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(12)*emit_par(12)
       CASE (7)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(12)*emit_par(12)
       CASE (8)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(6)*emit_par(6)
       CASE (9)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(5)*emit_par(5)
       CASE (10)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(7)*emit_par(7)
       CASE (11)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (12)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(10)*emit_par(10)
       CASE (13)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (14)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(11)*emit_par(11)
       CASE (15)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (16)
          eghg_bio(i,1,j,p_ebio_ch4term) = (1.0/3.6)*1.0e-12*biomt_par(13)*emit_par(13)
       CASE (17)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (18)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (19)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE (20)
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       CASE default
          eghg_bio(i,1,j,p_ebio_ch4term) = 0.0
       END SELECT

    END IF

    eghg_bio(i,1,j,p_ebio_ch4term) = eghg_bio(i,1,j,p_ebio_ch4term)*const

  END DO
 END DO


END SUBROUTINE TERMITE





  SUBROUTINE SOILUPTAKE            ( ids,ide, jds,jde,                      &
                                     ims,ime, jms,jme,                      &
                                     its,ite, jts,jte,                      &

                                     Soil_M, soil_typ, eghg_bio,            &
                                     rain_1, rain_2,                        &
                                     potevap, sfevap, landuse, T2, dt,      &
                                     num_soil_layers, wet_in                )

  

  INTEGER,  INTENT(IN   )   ::  ids,ide, jds,jde, &
                                ims,ime, jms,jme, &
                                its,ite, jts,jte, &
                                num_soil_layers
  INTEGER :: i,j

  INTEGER, DIMENSION (ims:ime, jms:jme), INTENT(IN)    :: soil_typ
    

  REAL, PARAMETER :: const= 3.6e9/(0.012+4.*0.001) 
  REAL, DIMENSION (ims:ime, jms:jme), INTENT(IN)  ::  rain_1, T2,       &
                                                      rain_2, landuse,  &
                                                      potevap, sfevap

  REAL, DIMENSION (ims:ime, num_soil_layers, jms:jme), INTENT(IN) :: Soil_M
  REAL, INTENT(IN) :: dt

  REAL, DIMENSION (ims:ime, 1,jms:jme, num_eghg_bio), INTENT(INOUT) ::  eghg_bio
  REAL, DIMENSION (ims:ime, 1,jms:jme, num_wet_in), INTENT(IN) ::  wet_in

  REAL, DIMENSION(ims:ime, jms:jme) :: sm_01, sm_02, sm_03, sm_04,    &
                                       i_cult, r_n, r_sm, vero, k_d, &
                                       G_soil, G_t, sm_res, dch4,     &
                                       phi_soil, eps, b_f, sand_c,    &
                                       clay_c, vero1

  REAL :: dch4_0, k_0, z_dsoil, conv_Fk, conv_ch4, methane_0

  CHARACTER*256 :: mminlu_loc



   dch4_0    =       0.196         
   k_0       =       0.00087       
                                   
   z_dsoil   =       6.0           
   conv_Fk   =       616.9         
   conv_ch4  =       (28.97/26.0)*1e6  

   methane_0 =       1.77     

 DO j=jts,min(jte,jde-1)
  DO i=its,min(ite,ide-1)



  sm_01(i, j)  = 0.0
  sm_02(i, j)  = 0.0
  sm_03(i, j)  = 0.0
  sm_04(i, j)  = 0.0
  i_cult(i, j) = 0.0
  r_n(i, j)    = 0.0
  r_sm(i, j)   = 0.0
  vero(i, j)   = 0.0
  k_d(i, j)    = 0.0
  G_soil(i, j) = 0.0
  G_t(i, j)    = 0.0
  sm_res(i, j) = 0.0
  dch4(i, j)   = 0.0
  phi_soil(i, j) = 0.0
  eps(i, j)    = 0.0
  b_f(i, j)    = 0.0
  sand_c(i, j) = 0.0
  clay_c(i, j) = 0.0
  vero1(i, j)  = 0.0



  sm_01(i,j)  = Soil_M(i,1,j)
  sm_02(i,j)  = Soil_M(i,2,j)
  sm_03(i,j)  = Soil_M(i,3,j)
  sm_04(i,j)  = Soil_M(i,4,j)





     G_t(i, j) = 1.0 + 0.0055 * (T2(i, j)-273.15)










  IF (soil_typ(i, j) == 1) then   
         phi_soil(i, j) = 0.339
         eps(i, j)      = 0.103
         sand_c(i,j)    = 0.92
         clay_c(i,j)    = 0.03
         b_f(i, j)      = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                          *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 2) then  
         phi_soil(i, j) = 0.421
         eps(i, j)      = 0.038
         sand_c(i,j)    = 0.82
         clay_c(i,j)    = 0.06
         b_f(i, j)      = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                           *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 3) then 
         phi_soil(i, j) = 0.434
         eps(i, j) = 0.051
         sand_c(i,j) = 0.58
         clay_c(i,j) = 0.10
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 4) then 
         phi_soil(i, j) = 0.476
         eps(i, j) = 0.116
         sand_c(i,j) = 0.17
         clay_c(i,j) = 0.13
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 5) then 
         phi_soil(i, j) = 0.476
         eps(i, j) = 0.093
         sand_c(i,j) = 0.10
         clay_c(i,j) = 0.10  
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 6) then  
         phi_soil(i, j) = 0.439
         eps(i, j) = 0.11
         sand_c(i,j) = 0.43
         clay_c(i,j) = 0.18
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 7) then  
         phi_soil(i, j) = 0.464
         eps(i, j) = 0.09
         sand_c(i,j) = 0.58
         clay_c(i,j) = 0.27
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 8) then  
         phi_soil(i, j) = 0.404
         eps(i, j) = 0.077
         sand_c(i,j) = 0.10
         clay_c(i,j) = 0.34
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 9) then  
         phi_soil(i, j) = 0.465
         eps(i, j) = 0.083
         sand_c(i,j) = 0.32
         clay_c(i,j) = 0.34
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 10) then  
         phi_soil(i, j) = 0.406
         eps(i, j) = 0.068
         sand_c(i,j) = 0.52
         clay_c(i,j) = 0.42
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 11) then 
         phi_soil(i, j) = 0.468
         eps(i, j) = 0.064
         sand_c(i,j) = 0.06
         clay_c(i,j) = 0.47
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 12) then  
         phi_soil(i, j) = 0.468
         eps(i, j) = 0.056
         sand_c(i,j) = 0.22
         clay_c(i,j) = 0.58
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 13) then 
         phi_soil(i, j) = 0.439
         eps(i, j) = 0.11
         sand_c(i,j) = 0.05
         clay_c(i,j) = 0.05 
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 14) then  
         phi_soil(i, j) = 1.0
         eps(i, j) = 0.0
         sand_c(i,j) = 0.60
         clay_c(i,j) = 0.40  
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 15) then  
         phi_soil(i, j) = 0.2
         eps(i, j) = 0.03
         sand_c(i,j) = 0.07
         clay_c(i,j) = 0.06
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 16) then  
         phi_soil(i, j) = 0.421
         eps(i, j) = 0.138
         sand_c(i,j) = 0.25
         clay_c(i,j) = 0.25
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 17) then  
         phi_soil(i, j) = 0.468
         eps(i, j) = 0.014
         sand_c(i,j) = 0.60
         clay_c(i,j) = 0.40
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
 ELSE IF (soil_typ(i, j) == 18) then  
         phi_soil(i, j) = 0.468
         eps(i, j) = 0.03
         sand_c(i,j) = 0.52
         clay_c(i,j) = 0.48
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE IF (soil_typ(i, j) == 19) then  
         phi_soil(i, j) = 0.339
         eps(i, j) = 0.103
         sand_c(i,j) = 0.92
         clay_c(i,j) = 0.03
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
  ELSE
         phi_soil(i, j) = 0.339
         eps(i, j) = 0.0
         sand_c(i,j) = 0.92
         clay_c(i,j) = 0.03
         b_f(i, j) = -3.140 - 0.00222* clay_c(i, j)**2 - 3.484*10e-5  &
                     *sand_c(i, j)**2 *clay_c(i, j)
   END IF
                   
   G_soil(i, j) = (phi_soil(i, j)**(4.0/3.0))* &
                  ((eps(i, j)/phi_soil(i, j))**(1.5 + 3.0/b_f(i, j)))



   dch4(i, j) = G_soil(i, j)*G_t(i, j)*dch4_0
                                                            







   IF(T2(i,j)<= 0.0) then

   vero(i,j) = 0.0

   ELSE IF(T2(i,j) > 0.0) then

   vero1(i, j) = (T2(i, j)-273.15)**4.0

   vero(i,j) = EXP((0.0693*(T2(i,j)-273.15))-(8.56*10e-7*vero1(i, j)))


   END IF



 CALL nl_get_mminlu(1, mminlu_loc)

    IF(mminlu_loc .EQ. 'USGS') THEN 

       IF(landuse(i, j) <= 5.0) THEN

            i_cult(i, j) = 1.0

       ELSE IF(landuse(i, j) == 6.0) THEN

            i_cult(i, j) = 0.5

       ELSE

            i_cult(i, j) = 0.0

       END IF

   ELSE IF(mminlu_loc .EQ. 'MODIFIED_IGBP_MODIS_NOAH')THEN

       IF(landuse(i, j) == 12.0 .OR. landuse(i,j) == 32.0 .OR. landuse(i,j) == 33.0 .OR. landuse(i,j) == 13.0) THEN

             i_cult(i, j) = 1.0

       ELSE IF(landuse(i, j) == 14.0 .OR. landuse(i,j) == 31.0) THEN

             i_cult(i, j) = 0.5

       ELSE

             i_cult(i, j) = 0.0

       END IF

  END IF 

  r_n(i, j) = 1.0 - (0.75*i_cult(i, j))






   IF(sfevap(i, j) == 0.0) THEN

    sm_res(i, j) = 2.0

  ELSE

     sm_res(i, j) = ((rain_1(i, j) + rain_2(i, j))*30.0*24.*3600./dt*0.001 &
                     + sm_02(i, j))/sfevap(i, j)

  END IF

  IF(sm_res(i, j) > 1.0)THEN

    r_sm(i, j) = 1.0

  ELSE IF(sm_res(i, j) <= 1.0)THEN

    r_sm(i, j) = sm_res(i, j)

  END IF

  k_d(i, j) = r_n(i, j)*r_sm(i, j)*vero(i, j)*k_0







  IF(dch4(i,j) /= 0.0 .AND. wet_in(i,1,j,p_wetmap) < 0.10 .AND. soil_typ(i,j) /= 14 .AND. soil_typ(i,j) /=1) then

  eghg_bio(i,1,j,p_ebio_ch4soil) = ((methane_0*dch4(i, j))/z_dsoil)* &
                                 (1.0 - dch4(i, j)/(dch4(i, j)+k_d(i, j)*z_dsoil))*conv_Fk


  eghg_bio(i,1,j,p_ebio_ch4soil) = -eghg_bio(i,1,j,p_ebio_ch4soil) * 1.0e-6/(24.0*3600.0)*const 


  ELSE

  eghg_bio(i,1,j,p_ebio_ch4soil) = 0.0

  END IF

  IF(eghg_bio(i,1,j,p_ebio_ch4soil) > 0.0) then

   eghg_bio(i,1,j,p_ebio_ch4soil) = 0.0

  END IF

  END DO
  END DO

  END SUBROUTINE SOILUPTAKE

END MODULE module_ghg_fluxes
