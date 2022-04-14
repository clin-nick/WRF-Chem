







      SUBROUTINE aerosols_driver (id,curr_secs,ktau,dtstep,ktauc,          &
               config_flags,dtstepc,dx,                                    &
               alt,t_phy,moist,aerwrf,p8w,t8w,p_phy,chem,rho_phy,dz8w,rh,  &
               z,z_at_w,pbl_h,cldfra,cldfra_mp_all,vbs_nbin,               &
               gamn2o5,cn2o5,kn2o5,yclno2,snu,sac,          &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,  &
               cvaro2,cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,vcsulf_old,&
               vdrog3, vdrog3_vbs,brch_ratio,dgnum,dgnumwet,wetdens_ap,    &
               del_h2so4_gasprod,dvmrdt_sv13d,dvmrcwdt_sv13d,              &
               is_CAMMGMP_used,                                            &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )


   USE module_configure
   USE module_state_description
   USE module_model_constants



   USE module_aerosols_sorgam
   USE module_gocart_aerosols
   USE module_data_sorgam
   USE module_mosaic_driver, only:  mosaic_aerchem_driver
   USE module_aerosols_soa_vbs, only: soa_vbs_driver
   USE module_aerosols_sorgam_vbs, only: sorgam_vbs_driver
   USE module_data_soa_vbs, only: ldrog_vbs
   USE module_cam_mam_aerchem_driver, only:  cam_mam_aerchem_driver
   USE modal_aero_data, only:  ntot_amode_cam_mam => ntot_amode
   USE module_cam_support, only: gas_pcnst => gas_pcnst_modal_aero, &
        gas_pcnst_pos => gas_pcnst_modal_aero_pos

   
   
   
   
   


   IMPLICIT NONE



































































   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
   LOGICAL,      INTENT(IN)       :: is_CAMMGMP_used
   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id,ktau,ktauc,vbs_nbin(1) 
   REAL(KIND=8), INTENT(IN   ) :: curr_secs
   REAL,         INTENT(IN   ) :: dtstep,dtstepc,dx



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),         &
         INTENT(IN ) ::                                   moist



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                    &
         INTENT(INOUT ) ::                                          &
           gamn2o5,cn2o5,kn2o5,yclno2,snu,sac,                   &
           h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
           cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,brch_ratio

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, ntot_amode_cam_mam ), &
         INTENT(INOUT ) ::                                   &
           dgnum, dgnumwet, wetdens_ap


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), &
         INTENT(IN ) ::                                   &
           del_h2so4_gasprod



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                    &
         INTENT(INOUT ) ::                                          &
               aerwrf



   REAL,  DIMENSION(ims:ime,kms:kme-0,jms:jme,ldrog),               &
               INTENT(IN   ) ::                                     &
                                                     VDROG3

   REAL,  DIMENSION(ims:ime,kms:kme-0,jms:jme,ldrog_vbs),           &
               INTENT(IN   ) ::                                     &
                                                     VDROG3_VBS


   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
          INTENT(IN   ) ::                                          &
                                                        alt,        &
                                                      t_phy,        &
                                                      p_phy,        &
                                                      dz8w,         &
                                                      rh,           & 
                                                      z    ,        &
                                            t8w,p8w,z_at_w ,        &
                                                    rho_phy,        &
                                                     cldfra,        &
                                                     cldfra_mp_all    

   REAL,  DIMENSION( ims:ime , jms:jme )                   ,        &
          INTENT(IN   ) ::                                          &
                                                      pbl_h



     REAL, dimension (ims:ime,kms:kme-0,jms:jme),                   &
               INTENT(INOUT) ::                                     &
                               vcsulf_old




   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, gas_pcnst_pos ),          &
        INTENT(IN ) ::                                  dvmrdt_sv13d,dvmrcwdt_sv13d 


     integer :: ii,jj,kk








   cps_select: SELECT CASE(config_flags%chem_opt)

   CASE (GOCART_SIMPLE,GOCARTRACM_KPP,GOCARTRADM2,MOZCART_KPP,T1_MOZCART_KPP)
      call gocart_aerosols_driver(ktauc,dtstepc,config_flags,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,dx,g,         &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
   CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_AQCHEM,RADM2SORG_KPP,CBMZSORG,CBMZSORG_AQ, &
          CB05_SORG_AQ_KPP) 
       CALL wrf_debug(15,'aerosols_driver calling sorgam_driver')
       do ii=its,ite
          do kk=kts,kte
             do jj=jts,jte
                if(chem(ii,kk,jj,p_nu0).lt.1.e07)then
                   chem(ii,kk,jj,p_nu0)=1.e7
                endif
             enddo
          enddo
       enddo
       call sorgam_driver (id,ktauc,dtstepc,t_phy,moist,aerwrf,p8w,t8w, &
               alt,p_phy,chem,rho_phy,dz8w,z,z_at_w,                    &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,vcsulf_old,    &
               vdrog3,                                                  &
               config_flags%kemit,                                      &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )
   CASE (CB05_SORG_VBS_AQ_KPP)
       CALL wrf_debug(15,'aerosols_driver calling sorgam_vbs_driver')
       do ii=its,ite
          do kk=kts,kte
             do jj=jts,jte
                if(chem(ii,kk,jj,p_nu0).lt.1.e07)then
                   chem(ii,kk,jj,p_nu0)=1.e7
                endif
             enddo
          enddo
       enddo
       call sorgam_vbs_driver (id,ktauc,dtstepc,t_phy,moist,aerwrf,p8w,t8w, &
            alt,p_phy,chem,rho_phy,dz8w,rh,z,z_at_w,                 &
            h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,               &
            vcsulf_old,vdrog3_vbs,                                   &
            config_flags%kemit,brch_ratio,                           &
            ids,ide, jds,jde, kds,kde,                               &
            ims,ime, jms,jme, kms,kme,                               &
            its,ite, jts,jte, kts,kte                                )

   CASE (RACMSORG_AQ,RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP,RACMSORG_KPP,RACM_ESRLSORG_KPP)

       CALL wrf_debug(15,'aerosols_driver calling sorgam_driver')
       do ii=its,ite
          do kk=kts,kte
             do jj=jts,jte
                if(chem(ii,kk,jj,p_nu0).lt.1.e07)then
                   chem(ii,kk,jj,p_nu0)=1.e7
                endif
             enddo
          enddo
       enddo
       call sorgam_driver (id,ktauc,dtstepc,t_phy,moist,aerwrf,p8w,t8w, &
               alt,p_phy,chem,rho_phy,dz8w,z,z_at_w,                    &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,vcsulf_old,    &
               vdrog3,                                                  &
               config_flags%kemit,                                      &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )

   CASE (CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_KPP, CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_4BIN_AQ, &
         CBMZ_MOSAIC_8BIN_AQ, SAPRC99_MOSAIC_4BIN_VBS2_KPP, &
         MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP, &
         SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_KPP, & 
         CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ)
       CALL wrf_debug(15,'aerosols_driver calling mosaic_aerchem_driver')
       CALL mosaic_aerchem_driver(                                      &
            id, curr_secs, ktau, dtstep, ktauc, dtstepc, config_flags,  &
            t_phy, rho_phy, p_phy,                                      &
            moist, chem,vbs_nbin,                                       &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )
   CASE ( RACM_SOA_VBS_KPP, RACM_SOA_VBS_AQCHEM_KPP, RACM_SOA_VBS_HET_KPP )
       CALL wrf_debug(15,'aerosols_driver calling soa_vbs_driver')
       do ii=its,ite
          do kk=kts,kte
             do jj=jts,jte
                if(chem(ii,kk,jj,p_nu0).lt.1.e07)then
                   chem(ii,kk,jj,p_nu0)=1.e7
                endif
             enddo
          enddo
       enddo
       call soa_vbs_driver ( id,ktauc,dtstepc,t_phy,moist,aerwrf,p8w,t8w, &
            alt,p_phy,chem,rho_phy,dz8w,rh,z,z_at_w,                 &
            gamn2o5,cn2o5,kn2o5,yclno2,snu,sac,                    &
            h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,               &
            vcsulf_old,vdrog3_vbs,                                   &
            config_flags%kemit,brch_ratio,                           &
            ids,ide, jds,jde, kds,kde,                               &
            ims,ime, jms,jme, kms,kme,                               &
            its,ite, jts,jte, kts,kte                                ) 

   CASE (CBMZ_CAM_MAM3_NOAQ, CBMZ_CAM_MAM3_AQ, CBMZ_CAM_MAM7_NOAQ, CBMZ_CAM_MAM7_AQ)
       CALL wrf_debug(15,'aerosols_driver calling cam_mam_aerchem_driver')
       CALL cam_mam_aerchem_driver(                                     &
            id, curr_secs, ktau, dtstep, ktauc, dtstepc, config_flags,  &
            t_phy, rho_phy, p_phy, p8w, alt, z, z_at_w, pbl_h, cldfra,  &
            cldfra_mp_all, moist, chem,                                 &
            dgnum, dgnumwet, wetdens_ap, del_h2so4_gasprod,             &
            dvmrdt_sv13d,dvmrcwdt_sv13d,                                & 
            is_CAMMGMP_used,                                            &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )


   CASE DEFAULT 

   END SELECT cps_select

   END SUBROUTINE aerosols_driver







   SUBROUTINE sum_pm_driver ( config_flags,                            &
             alt, chem, h2oaj, h2oai,                                  &
             pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,               &
             tsoa,asoa,bsoa,                                           &
             hoa_a01,hoa_a02,hoa_a03,hoa_a04,                          &
             hoa_a05,hoa_a06,hoa_a07,hoa_a08,                          & 
             bboa_a01,bboa_a02,bboa_a03,bboa_a04,                      &
             bboa_a05,bboa_a06,bboa_a07,bboa_a08,                      &
             soa_a01,soa_a02,soa_a03,soa_a04,                          &
             soa_a05,soa_a06,soa_a07,soa_a08,                          &
             bbsoa_a01,bbsoa_a02,bbsoa_a03,bbsoa_a04,                  &
             bbsoa_a05,bbsoa_a06,bbsoa_a07,bbsoa_a08,                  &
             hsoa_a01,hsoa_a02,hsoa_a03,hsoa_a04,                      &
             hsoa_a05,hsoa_a06,hsoa_a07,hsoa_a08,                      &
             biog_a01,biog_a02,biog_a03,biog_a04,                      &
             biog_a05,biog_a06,biog_a07,biog_a08,                      &
             asmpsoa_a01,asmpsoa_a02,asmpsoa_a03,asmpsoa_a04,          &
             arosoa_a01,arosoa_a02,arosoa_a03,arosoa_a04,              &
             arosoa_a05,arosoa_a06,arosoa_a07,arosoa_a08,              &
             totoa_a01,totoa_a02,totoa_a03,totoa_a04,                  &
             totoa_a05,totoa_a06,totoa_a07,totoa_a08,                  &
             hsoa_c,hsoa_o,bbsoa_c,bbsoa_o,                            &
             biog_v1,biog_v2,biog_v3,biog_v4,                          &
             ant_v1,ant_v2,ant_v3,ant_v4,                              &
             smpa_v1,smpbb_v1,                                         &
             
             hoa_cw01,    hoa_cw02,    hoa_cw03,    hoa_cw04,          &
             hoa_cw05,    hoa_cw06,    hoa_cw07,    hoa_cw08,          &
             bboa_cw01,   bboa_cw02,   bboa_cw03,   bboa_cw04,         &
             bboa_cw05,   bboa_cw06,   bboa_cw07,   bboa_cw08,         &
             soa_cw01,    soa_cw02,    soa_cw03,    soa_cw04,          &
             soa_cw05,    soa_cw06,    soa_cw07,    soa_cw08,          &
             bbsoa_cw01,  bbsoa_cw02,  bbsoa_cw03,  bbsoa_cw04,        &
             bbsoa_cw05,  bbsoa_cw06,  bbsoa_cw07,  bbsoa_cw08,        &
             biog_cw01,   biog_cw02,   biog_cw03,   biog_cw04,         &
             biog_cw05,   biog_cw06,   biog_cw07,   biog_cw08,         &
             hsoa_cw01,   hsoa_cw02,   hsoa_cw03,   hsoa_cw04,         &
             hsoa_cw05,   hsoa_cw06,   hsoa_cw07,   hsoa_cw08,         &
             arosoa_cw01, arosoa_cw02, arosoa_cw03, arosoa_cw04,       &
             arosoa_cw05, arosoa_cw06, arosoa_cw07, arosoa_cw08,       &
             totoa_cw01,  totoa_cw02,  totoa_cw03,  totoa_cw04,        &
             totoa_cw05,  totoa_cw06,  totoa_cw07,  totoa_cw08,        &        
             hsoa_cw_c,   hsoa_cw_o,   bbsoa_cw_c,  bbsoa_cw_o,        &
             biog_cw_v1,                                               &
             ant_cw_v1,                                                &
             
             ids,ide, jds,jde, kds,kde,                                &
             ims,ime, jms,jme, kms,kme,                                &
             its,ite, jts,jte, kts,kte                                  )


   USE module_configure
   USE module_aerosols_sorgam, only: sum_pm_sorgam
   USE module_mosaic_driver, only: sum_pm_mosaic,sum_pm_mosaic_vbs2,sum_pm_mosaic_vbs0,sum_pm_mosaic_vbs4,&
                                   sum_vbs9,sum_vbs2,sum_vbs0,sum_vbs4,sum_aq_vbs2
   USE module_gocart_aerosols, only: sum_pm_gocart
   USE module_aerosols_soa_vbs, only: sum_pm_soa_vbs
   USE module_aerosols_sorgam_vbs, only: sum_pm_sorgam_vbs

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::                                   &
                                      ids,ide, jds,jde, kds,kde,       &
                                      ims,ime, jms,jme, kms,kme,       &
                                      its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(IN ) :: chem

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(IN ) :: alt
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         OPTIONAL,                                                     &
         INTENT(IN ) :: h2oaj,h2oai

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         OPTIONAL,                                                     &
         INTENT(OUT) :: pm2_5_dry,pm2_5_water,pm2_5_dry_ec,pm10,       &
             tsoa,asoa,bsoa,                                           &
             hoa_a01,hoa_a02,hoa_a03,hoa_a04,                          &
             hoa_a05,hoa_a06,hoa_a07,hoa_a08,                          &
             bboa_a01,bboa_a02,bboa_a03,bboa_a04,                      &
             bboa_a05,bboa_a06,bboa_a07,bboa_a08,                      & 
             soa_a01,soa_a02,soa_a03,soa_a04,                          &
             soa_a05,soa_a06,soa_a07,soa_a08,                          & 
             bbsoa_a01,bbsoa_a02,bbsoa_a03,bbsoa_a04,                  &
             bbsoa_a05,bbsoa_a06,bbsoa_a07,bbsoa_a08,                  &       
             hsoa_a01,hsoa_a02,hsoa_a03,hsoa_a04,                      &
             hsoa_a05,hsoa_a06,hsoa_a07,hsoa_a08,                      &  
             biog_a01,biog_a02,biog_a03,biog_a04,                      &
             biog_a05,biog_a06,biog_a07,biog_a08,                      & 
             arosoa_a01,arosoa_a02,arosoa_a03,arosoa_a04,              &
             arosoa_a05,arosoa_a06,arosoa_a07,arosoa_a08,              &   
             totoa_a01,totoa_a02,totoa_a03,totoa_a04,                  &
             totoa_a05,totoa_a06,totoa_a07,totoa_a08,                  &    
             hsoa_c,hsoa_o,bbsoa_c,bbsoa_o,                            &
             biog_v1,biog_v2,biog_v3,biog_v4,                          &
             ant_v1,ant_v2,ant_v3,ant_v4,                              &
             smpa_v1,                                                  &
             smpbb_v1,                                                 &
             asmpsoa_a01,asmpsoa_a02,asmpsoa_a03,asmpsoa_a04,       &
             
             hoa_cw01,    hoa_cw02,    hoa_cw03,    hoa_cw04,          &
             hoa_cw05,    hoa_cw06,    hoa_cw07,    hoa_cw08,          &
             bboa_cw01,   bboa_cw02,   bboa_cw03,   bboa_cw04,         &
             bboa_cw05,   bboa_cw06,   bboa_cw07,   bboa_cw08,         &
             soa_cw01,    soa_cw02,    soa_cw03,    soa_cw04,          &
             soa_cw05,    soa_cw06,    soa_cw07,    soa_cw08,          &
             bbsoa_cw01,  bbsoa_cw02,  bbsoa_cw03,  bbsoa_cw04,        &
             bbsoa_cw05,  bbsoa_cw06,  bbsoa_cw07,  bbsoa_cw08,        &
             biog_cw01,   biog_cw02,   biog_cw03,   biog_cw04,         &
             biog_cw05,   biog_cw06,   biog_cw07,   biog_cw08,         &
             hsoa_cw01,   hsoa_cw02,   hsoa_cw03,   hsoa_cw04,         &
             hsoa_cw05,   hsoa_cw06,   hsoa_cw07,   hsoa_cw08,         &
             arosoa_cw01, arosoa_cw02, arosoa_cw03, arosoa_cw04,       &
             arosoa_cw05, arosoa_cw06, arosoa_cw07, arosoa_cw08,       &
             totoa_cw01,  totoa_cw02,  totoa_cw03,  totoa_cw04,        &
             totoa_cw05,  totoa_cw06,  totoa_cw07,  totoa_cw08,        &        
             hsoa_cw_c,   hsoa_cw_o,   bbsoa_cw_c,  bbsoa_cw_o,        &
             biog_cw_v1,                                               &
             ant_cw_v1                                                 
             
 

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags




   sum_pm_select: SELECT CASE(config_flags%chem_opt)

   CASE (GOCART_SIMPLE,GOCARTRACM_KPP,GOCARTRADM2,MOZCART_KPP,T1_MOZCART_KPP)
       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_gocart')
       CALL sum_pm_gocart (                                            &
            alt, chem,pm2_5_dry, pm2_5_dry_ec, pm10,                   &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )
   CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_AQCHEM,RACMSORG_AQ,RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP,RADM2SORG_KPP, &
         RACMSORG_KPP,RACM_ESRLSORG_KPP,CBMZSORG,CBMZSORG_AQ,CB05_SORG_AQ_KPP)
       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_sorgam')
       CALL sum_pm_sorgam (                                            &
            alt, chem, h2oaj, h2oai,                                   &
            pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                &
            config_flags%dust_opt,ids,ide, jds,jde, kds,kde,           &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )
   CASE (CB05_SORG_VBS_AQ_KPP)
       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_sorgam_vbs')
       CALL sum_pm_sorgam_vbs (                                        &
            alt, chem, h2oaj, h2oai,                                   &
            pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                &
            tsoa,asoa,bsoa,                                            &
            config_flags%dust_opt,ids,ide, jds,jde, kds,kde,           &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )
   CASE (RACM_SOA_VBS_KPP,RACM_SOA_VBS_AQCHEM_KPP,RACM_SOA_VBS_HET_KPP)
       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_soa_vbs')
       CALL sum_pm_soa_vbs (                                           &
            alt, chem, h2oaj, h2oai,                                   &
            pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                &
            config_flags%dust_opt,ids,ide, jds,jde, kds,kde,           &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

   CASE (CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_KPP, CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
         CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
         CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP)
       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_mosaic')
       call sum_pm_mosaic (                                            &
            alt, chem,                                                 &
            pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

   CASE (SAPRC99_MOSAIC_4BIN_VBS2_KPP)

       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_mosaic_vbs2')
       call sum_pm_mosaic_vbs2 (                                       &
            alt, chem,                                                 &
            pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

       CALL wrf_debug(15,'sum_pm_driver: calling sum_vbs2')
       call sum_vbs2 ( config_flags%aero_diag_opt,                     &
             alt, chem,                                                &
             hoa_a01,hoa_a02,hoa_a03,hoa_a04,                          &
             hoa_a05,hoa_a06,hoa_a07,hoa_a08,                          &
             bboa_a01,bboa_a02,bboa_a03,bboa_a04,                      &
             bboa_a05,bboa_a06,bboa_a07,bboa_a08,                      &
             soa_a01,soa_a02,soa_a03,soa_a04,                          &
             soa_a05,soa_a06,soa_a07,soa_a08,                          &
             bbsoa_a01,bbsoa_a02,bbsoa_a03,bbsoa_a04,                  &
             bbsoa_a05,bbsoa_a06,bbsoa_a07,bbsoa_a08,                  &
             hsoa_a01,hsoa_a02,hsoa_a03,hsoa_a04,                      &
             hsoa_a05,hsoa_a06,hsoa_a07,hsoa_a08,                      &
             biog_a01,biog_a02,biog_a03,biog_a04,                      &
             biog_a05,biog_a06,biog_a07,biog_a08,                      &
             arosoa_a01,arosoa_a02,arosoa_a03,arosoa_a04,              &
             arosoa_a05,arosoa_a06,arosoa_a07,arosoa_a08,              &
             totoa_a01,totoa_a02,totoa_a03,totoa_a04,                  &
             totoa_a05,totoa_a06,totoa_a07,totoa_a08,                  &
             hsoa_c,hsoa_o,bbsoa_c,bbsoa_o,                            &
             biog_v1,biog_v2,biog_v3,biog_v4,                          &
             ant_v1,ant_v2,ant_v3,ant_v4,                              &
             ids,ide, jds,jde, kds,kde,                                &
             ims,ime, jms,jme, kms,kme,                                &
             its,ite, jts,jte, kts,kte                                  )

   CASE (MOZART_MOSAIC_4BIN_KPP)

       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_mosaic_vbs0')
       call sum_pm_mosaic_vbs0 (                                       &
            alt, chem,                                                 &
            pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

       CALL wrf_debug(15,'sum_pm_driver: calling sum_vbs0')
       call sum_vbs0 ( config_flags%aero_diag_opt,                     &
             alt, chem,                                                &
             hoa_a01,hoa_a02,hoa_a03,hoa_a04,                          &
             bboa_a01,bboa_a02,bboa_a03,bboa_a04,                      &
             soa_a01,soa_a02,soa_a03,soa_a04,                          &
             bbsoa_a01,bbsoa_a02,bbsoa_a03,bbsoa_a04,                  &
             biog_a01,biog_a02,biog_a03,biog_a04,                      &
             asmpsoa_a01,asmpsoa_a02,asmpsoa_a03,asmpsoa_a04,              &
             arosoa_a01,arosoa_a02,arosoa_a03,arosoa_a04,              &
             totoa_a01,totoa_a02,totoa_a03,totoa_a04,                  &
             biog_v1,biog_v2,biog_v3,biog_v4,                          &
             ant_v1,ant_v2,ant_v3,ant_v4,                              &
             smpa_v1,smpbb_v1,                              &
             ids,ide, jds,jde, kds,kde,                                &
             ims,ime, jms,jme, kms,kme,                                &
             its,ite, jts,jte, kts,kte                                  )

   CASE (MOZART_MOSAIC_4BIN_AQ_KPP)

       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_mosaic_vbs4')
       call sum_pm_mosaic_vbs4 (                                       &
            alt, chem,                                                 &
            pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

       CALL wrf_debug(15,'sum_pm_driver: calling sum_vbs4')
       call sum_vbs4 ( config_flags%aero_diag_opt,                     &
             alt, chem,                                                &
             hoa_a01,hoa_a02,hoa_a03,hoa_a04,                          &
             soa_a01,soa_a02,soa_a03,soa_a04,                          &
             biog_a01,biog_a02,biog_a03,biog_a04,                      &
             totoa_a01,totoa_a02,totoa_a03,totoa_a04,                  &
             biog_v1,biog_v2,biog_v3,biog_v4,                          &
             ant_v1,ant_v2,ant_v3,ant_v4,                              &
             ids,ide, jds,jde, kds,kde,                                &
             ims,ime, jms,jme, kms,kme,                                &
             its,ite, jts,jte, kts,kte                                 )
       
    CASE (SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_8BIN_VBS2_KPP)
       
       CALL wrf_debug(15,'sum_pm_driver: calling sum_pm_mosaic_vbs2')
       call sum_pm_mosaic_vbs2 (                                       &
            alt, chem,                                                 &
            pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

       CALL wrf_debug(15,'sum_pm_driver: calling sum_vbs2')
       call sum_vbs2 ( config_flags%aero_diag_opt,                     &
             alt, chem,                                                &
             hoa_a01,hoa_a02,hoa_a03,hoa_a04,                          &
             hoa_a05,hoa_a06,hoa_a07,hoa_a08,                          &
             bboa_a01,bboa_a02,bboa_a03,bboa_a04,                      &
             bboa_a05,bboa_a06,bboa_a07,bboa_a08,                      &
             soa_a01,soa_a02,soa_a03,soa_a04,                          &
             soa_a05,soa_a06,soa_a07,soa_a08,                          &
             bbsoa_a01,bbsoa_a02,bbsoa_a03,bbsoa_a04,                  &
             bbsoa_a05,bbsoa_a06,bbsoa_a07,bbsoa_a08,                  &
             hsoa_a01,hsoa_a02,hsoa_a03,hsoa_a04,                      &
             hsoa_a05,hsoa_a06,hsoa_a07,hsoa_a08,                      &
             biog_a01,biog_a02,biog_a03,biog_a04,                      &
             biog_a05,biog_a06,biog_a07,biog_a08,                      &
             arosoa_a01,arosoa_a02,arosoa_a03,arosoa_a04,              &
             arosoa_a05,arosoa_a06,arosoa_a07,arosoa_a08,              &
             totoa_a01,totoa_a02,totoa_a03,totoa_a04,                  &
             totoa_a05,totoa_a06,totoa_a07,totoa_a08,                  &
             hsoa_c,hsoa_o,bbsoa_c,bbsoa_o,                            &
             biog_v1,biog_v2,biog_v3,biog_v4,                          &
             ant_v1,ant_v2,ant_v3,ant_v4,                              &
             ids,ide, jds,jde, kds,kde,                                &
             ims,ime, jms,jme, kms,kme,                                &
             its,ite, jts,jte, kts,kte                                  )



       IF( config_flags%aero_cw_diag_opt == diag_cw_aero ) THEN
       CALL wrf_debug(15,'sum_pm_driver: calling sum_aq_vbs2')
       call sum_aq_vbs2 (                                                  &
             alt, chem,                                                    &
             hoa_cw01,hoa_cw02,hoa_cw03,hoa_cw04,hoa_cw05,hoa_cw06,hoa_cw07,hoa_cw08,                          &
             bboa_cw01,bboa_cw02,bboa_cw03,bboa_cw04,bboa_cw05,bboa_cw06,bboa_cw07,bboa_cw08,                  &
             soa_cw01,soa_cw02,soa_cw03,soa_cw04,soa_cw05,soa_cw06,soa_cw07,soa_cw08,                          &
             bbsoa_cw01,bbsoa_cw02,bbsoa_cw03,bbsoa_cw04,bbsoa_cw05,bbsoa_cw06,bbsoa_cw07,bbsoa_cw08,          &
             hsoa_cw01,hsoa_cw02,hsoa_cw03,hsoa_cw04,hsoa_cw05,hsoa_cw06,hsoa_cw07,hsoa_cw08,                  &
             biog_cw01,biog_cw02,biog_cw03,biog_cw04,biog_cw05,biog_cw06,biog_cw07,biog_cw08,                  &
             arosoa_cw01,arosoa_cw02,arosoa_cw03,arosoa_cw04,arosoa_cw05,arosoa_cw06,arosoa_cw07,arosoa_cw08,  &
             totoa_cw01,totoa_cw02,totoa_cw03,totoa_cw04,totoa_cw05,totoa_cw06,totoa_cw07,totoa_cw08,          &
             hsoa_cw_c,hsoa_cw_o,bbsoa_cw_c,bbsoa_cw_o,                    &
             biog_cw_v1,                                                   &
             ant_cw_v1,                                                    &
             ids,ide, jds,jde, kds,kde,                                &
             ims,ime, jms,jme, kms,kme,                                &
             its,ite, jts,jte, kts,kte                                  )
       ENDIF
       

   CASE DEFAULT 

   END SELECT sum_pm_select

   END SUBROUTINE sum_pm_driver
