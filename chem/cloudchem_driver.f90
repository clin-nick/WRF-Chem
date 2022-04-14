











      SUBROUTINE cloudchem_driver(                                   &
               id, ktau, ktauc, dtstep, dtstepc, config_flags,       &
               t_phy, p_phy, rho_phy, alt, dz8w,                           &
               p8w, prain3d,scalar,dvmrdt_sv13d,dvmrcwdt_sv13d,      & 
               f_ice_phy, f_rain_phy, cldfrai, cldfral,              &
	       moist, cldfra, cldfra_mp_all, ph_no2,                 &
	       chem, gas_aqfrac, numgas_aqfrac,                      &
               is_CAMMGMP_used,                                      &
               ids,ide, jds,jde, kds,kde,                            &
               ims,ime, jms,jme, kms,kme,                            &
               its,ite, jts,jte, kts,kte                             )







   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE module_cam_support, only: gas_pcnst => gas_pcnst_modal_aero, &
        gas_pcnst_pos => gas_pcnst_modal_aero_pos
   USE module_mosaic_cloudchem,  only: mosaic_cloudchem_driver
   USE module_sorgam_cloudchem,  only: sorgam_cloudchem_driver
   USE module_sorgam_vbs_cloudchem, only: sorgam_vbs_cloudchem_driver
   USE module_cam_mam_cloudchem, only: cam_mam_cloudchem_driver
   USE module_sorgam_aqchem, only: sorgam_aqchem_driver
   USE module_sorgam_vbs_aqchem, only: sorgam_vbs_aqchem_driver

   
   
   
   


   IMPLICIT NONE









































































   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
   LOGICAL,      INTENT(IN)       :: is_CAMMGMP_used
   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc,              &
                                      numgas_aqfrac
      REAL,      INTENT(IN   ) :: dtstep, dtstepc



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),         &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_scalar ),         &
         INTENT(IN ) ::                                   scalar    




   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, gas_pcnst_pos ),          &
        INTENT(OUT ) ::                                  dvmrdt_sv13d,dvmrcwdt_sv13d 


   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
          INTENT(IN   ) ::                                          &
                                t_phy,                              &
                                p_phy,                              &
                                rho_phy,                            &
                                alt,                                &
                                dz8w,                               &
                                cldfra,                             &
                                ph_no2,                             &
                                p8w,                                & 
                                prain3d,                            &
                                F_ICE_PHY,                          &
                                F_RAIN_PHY,                         &
                                cldfrai,                            &
                                cldfral,                            &
                                cldfra_mp_all



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(INOUT ) ::                                gas_aqfrac




     integer :: ii,jj,kk













   cps_select: SELECT CASE(config_flags%chem_opt)

   CASE ( CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
        CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, CBMZ_MOSAIC_DMS_4BIN_AQ,            &
        CBMZ_MOSAIC_DMS_8BIN_AQ, CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP, &
        MOZART_MOSAIC_4BIN_AQ_KPP,        &
        SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP)

       call wrf_debug(15, &
       'cloudchem_driver calling mosaic_cloudchem_driver')
       call mosaic_cloudchem_driver(                  &
            id, ktau, ktauc, dtstepc, config_flags,   &
            p_phy, t_phy, rho_phy, alt,               &
            cldfra, ph_no2,                           &
            moist, chem,                              &
            gas_aqfrac, numgas_aqfrac,                &
            ids,ide, jds,jde, kds,kde,                &
            ims,ime, jms,jme, kms,kme,                &
            its,ite, jts,jte, kts,kte )

   CASE ( RADM2SORG_AQ, RACMSORG_AQ, CBMZSORG_AQ )

       call wrf_debug(15, &
       'cloudchem_driver calling sorgam_cloudchem_driver')
       call sorgam_cloudchem_driver(                  &
            id, ktau, ktauc, dtstepc, config_flags,   &
            p_phy, t_phy, rho_phy, alt,               &
            cldfra, ph_no2,                           &
            moist, chem,                              &
            gas_aqfrac, numgas_aqfrac,                &
            ids,ide, jds,jde, kds,kde,                &
            ims,ime, jms,jme, kms,kme,                &
            its,ite, jts,jte, kts,kte )
    CASE (CBMZ_CAM_MAM3_NOAQ,CBMZ_CAM_MAM3_AQ,CBMZ_CAM_MAM7_NOAQ,CBMZ_CAM_MAM7_AQ)       
       CALL wrf_debug(15,'cloudchem_driver calling mam_cloudchem_driver')       
       call cam_mam_cloudchem_driver (                &
            
            dvmrdt_sv13d,dvmrcwdt_sv13d,              & 
            
            chem,                                     &
            
            moist, scalar, p8w, prain3d, p_phy,       &
            t_phy, dtstepc, ktau,alt, f_ice_phy,      &
            f_rain_phy, cldfra, cldfra_mp_all,        &
            cldfrai, cldfral, is_CAMMGMP_used,        & 
            ids,ide, jds,jde, kds,kde,                &
            ims,ime, jms,jme, kms,kme,                &
            its,ite, jts,jte, kts,kte                 )
   CASE ( CB05_SORG_VBS_AQ_KPP )

       call wrf_debug(15, &
       'cloudchem_driver calling sorgam_vbs_aqchem_driver')
       call sorgam_vbs_aqchem_driver(                 &
            id, ktau, ktauc, dtstepc, config_flags,   &
            p_phy, t_phy, rho_phy, alt, dz8w,         &
            moist, chem,                              &
            gas_aqfrac, numgas_aqfrac,                &
            ids,ide, jds,jde, kds,kde,                &
            ims,ime, jms,jme, kms,kme,                &
            its,ite, jts,jte, kts,kte )


   CASE ( RADM2SORG_AQCHEM, RACMSORG_AQCHEM_KPP, RACM_ESRLSORG_AQCHEM_KPP, &
          CB05_SORG_AQ_KPP,RACM_SOA_VBS_AQCHEM_KPP )

       call wrf_debug(15, &
       'cloudchem_driver calling sorgam_aqchem_driver')
       call sorgam_aqchem_driver(                  &
            id, ktau, ktauc, dtstepc, config_flags,   &
            p_phy, t_phy, rho_phy, alt, dz8w,         &
            moist, chem,                              &
            gas_aqfrac, numgas_aqfrac,                &
            ids,ide, jds,jde, kds,kde,                &
            ims,ime, jms,jme, kms,kme,                &
            its,ite, jts,jte, kts,kte )

   CASE DEFAULT

   END SELECT cps_select

   END SUBROUTINE cloudchem_driver

