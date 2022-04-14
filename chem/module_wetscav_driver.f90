








MODULE module_wetscav_driver

  REAL, PARAMETER :: mwso4 = 96.00   
  REAL, PARAMETER :: mwno3 = 62.0    

CONTAINS




      subroutine wetscav_driver( id, ktau, dtstep, ktauc, config_flags,    &
               dtstepc, alt, t_phy, moist, p8w,                            &
               t8w, dx, dy, p_phy, chem,                                   &
               rho_phy, cldfra, cldfra2, rainprod, evapprod,               &
               hno3_col_mdel, qlsink, precr, preci, precs, precg,          &
               wdflx, &
	       gas_aqfrac, numgas_aqfrac, dz8w,                            &
	       h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,           &
               cvaro2,cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,           &
               wd_no3,wd_so4,wd_nh4,wd_oa,                                 &
               wd_so2, wd_sulf, wd_hno3, wd_nh3,                           &
               wd_cvasoa, wd_cvbsoa, wd_asoa, wd_bsoa,                     &
               qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp,                         &


               p_hyd,scalar,dgnum4d,dgnumwet4d,dlf3d,dlf2_3d,qme3d,prain3d,&
               nevapr3d,rate1ord_cw2pr_st3d,shfrc3d,cmfmc,cmfmc2,evapcsh,  &
               icwmrsh,rprdsh,evapcdp3d,icwmrdp3d,rprddp3d,fracis3d,       &
               f_ice_phy,f_rain_phy,cldfrai,cldfral,cldfra_mp_all,         &
               is_CAMMGMP_used,                                            &

               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
















   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE modal_aero_data,        only:  ntot_amode               
   USE module_mozcart_wetscav, only:  wetscav_mozcart
   USE module_mosaic_wetscav, only:  wetscav_cbmz_mosaic, wetscav_mozart_mosaic

   
   USE module_aerosols_sorgam_vbs, only:  wetscav_sorgam_vbs_driver
   USE module_cam_mam_wetscav, only:  wetscav_cam_mam_driver
   USE module_cam_support,     only:  pcnst =>pcnst_runtime

   


   USE module_prep_wetscav_sorgam, only: wetscav_sorgam_driver, wetscav_soa_vbs_driver
   USE module_data_mosaic_asect, only: mw_so4_aer, mw_no3_aer, mw_nh4_aer, &
                                       mw_smpa_aer, mw_smpbb_aer, mw_oc_aer, &
                                       mw_glysoa_r1_aer, mw_glysoa_r2_aer, &
                                       mw_glysoa_sfc_aer, mw_glysoa_nh4_aer, &
                                       mw_glysoa_oh_aer,&
                         mw_asoaX_aer, mw_asoa1_aer, mw_asoa2_aer, mw_asoa3_aer, mw_asoa4_aer, &
                         mw_bsoaX_aer, mw_bsoa1_aer, mw_bsoa2_aer, mw_bsoa3_aer, mw_bsoa4_aer, &
                                       mw_biog1_c_aer, mw_biog1_o_aer

   IMPLICIT NONE



























































































   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
   LOGICAL,      INTENT(IN)       :: is_CAMMGMP_used
   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc, numgas_aqfrac

      REAL,      INTENT(IN   ) :: dtstep,dtstepc
      REAL,      INTENT(IN   ) :: dx, dy



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),         &
         INTENT(INOUT ) ::                                   moist



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_scalar ),         &
         INTENT(INOUT) ::                                   scalar




   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(IN ) ::                                   gas_aqfrac




   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                    &
         INTENT(INOUT ) ::                                          &
           h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
           cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2



   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: wd_no3,wd_so4, &
                                                          wd_nh4, wd_oa, &
                                                          wd_so2, wd_sulf, &
                                                          wd_hno3, wd_nh3
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: &
                   wd_cvasoa, wd_cvbsoa, wd_asoa, wd_bsoa


   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
          INTENT(IN   ) ::                                          &
                                                      alt,          &
                                                      t_phy,        &
                                                      p_phy,        &
                                                      t8w,          &
                                                      p8w,          &
                                                      dz8w,         &
	                            qlsink,precr,preci,precs,precg, &
                                                    rho_phy,cldfra, &

                                     cldfrai,cldfral,cldfra_mp_all, &
                                               p_hyd,dlf3d,dlf2_3d, & 
                                            qme3d,prain3d,nevapr3d, & 
                                       rate1ord_cw2pr_st3d,shfrc3d, &
                                    cmfmc,cmfmc2,evapcsh,icwmrsh,   &
                                        rprdsh,evapcdp3d,icwmrdp3d, &
                                     rprddp3d,f_ice_phy,f_rain_phy 
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
         INTENT(INOUT ) ::                             cldfra2,     &
                                                       rainprod,    &
                                                       evapprod
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
         INTENT(IN) ::                                 qc_b4mp,     &
                                                       qv_b4mp,     &
                                                       qi_b4mp,     &
                                                       qs_b4mp
   REAL,  DIMENSION( ims:ime , jms:jme )         ,                  &
         INTENT(INOUT ) ::                             hno3_col_mdel

   REAL, DIMENSION( ims:ime, jms:jme, num_chem )            ,       &
          INTENT(INOUT) ::                                          &
                                                       wdflx 


REAL,  DIMENSION( ims:ime , kms:kme , jms:jme , ntot_amode ),       &
          INTENT(IN   ) ::                       dgnum4d,dgnumwet4d   



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme, pcnst )         ,        &
          INTENT(OUT  ) ::                                          &
                                                        fracis3d


   integer :: ii,jj,kk
   REAL, DIMENSION( ims:ime, jms:jme, num_chem ) :: qsrflx, delta_mass_col 
   REAL ::  tmp_minval = 1.0e7
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) :: rainrate, evaprate




  REAL, DIMENSION( ims:ime , jms:jme ) :: wdi_no3,wdi_so4,wdi_nh4,wdi_oa, &
                                          wdi_so2, wdi_sulf, wdi_hno3, wdi_nh3
  REAL, DIMENSION( ims:ime , jms:jme ) :: wdi_cvasoa,wdi_cvbsoa,wdi_asoa,wdi_bsoa












   cps_select: SELECT CASE(config_flags%chem_opt)

   CASE ( RADM2SORG_AQ, RADM2SORG_AQCHEM, RACMSORG_AQ, RACMSORG_AQCHEM_KPP, RACM_ESRLSORG_AQCHEM_KPP, CBMZSORG_AQ, &
          CB05_SORG_AQ_KPP)
       CALL wrf_debug(15,'wetscav_driver calling sorgam_wetscav_driver' )
       call wetscav_sorgam_driver (id,ktau,dtstep,ktauc,config_flags,      &
               dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,        &
               qlsink,precr,preci,precs,precg, qsrflx,                     &
               gas_aqfrac, numgas_aqfrac,                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )

       tmp_minval = 1.0e7
       do jj=jts,jte
       do kk=kts,kte
       do ii=its,ite
          if (chem(ii,kk,jj,p_nu0) .lt. tmp_minval) then
             chem(ii,kk,jj,p_nu0) = tmp_minval
          endif
       enddo
       enddo
       enddo


     CASE ( RACM_SOA_VBS_AQCHEM_KPP )
       CALL wrf_debug(15,'wetscav_driver calling soa_vbs_wetscav_driver' )
       call wetscav_soa_vbs_driver (id,ktau,dtstep,ktauc,config_flags,     &
               dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,        &
               qlsink,precr,preci,precs,precg, qsrflx,                     &
               gas_aqfrac, numgas_aqfrac,                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )

       tmp_minval = 1.0e7
       do jj=jts,jte
       do kk=kts,kte
       do ii=its,ite
          if (chem(ii,kk,jj,p_nu0) .lt. tmp_minval) then
             chem(ii,kk,jj,p_nu0) = tmp_minval
          endif
       enddo
       enddo
       enddo

   CASE (CB05_SORG_VBS_AQ_KPP )
       CALL wrf_debug(15,'wetscav_driver calling sorgam_wetscav_driver' )
       call wetscav_sorgam_vbs_driver (id,ktau,dtstep,ktauc,config_flags,      &
               dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,        &
               qlsink,precr,preci,precs,precg, qsrflx,                     &
               gas_aqfrac, numgas_aqfrac,                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )

       tmp_minval = 1.0e7
       do jj=jts,jte
       do kk=kts,kte
       do ii=its,ite
          if (chem(ii,kk,jj,p_nu0) .lt. tmp_minval) then
             chem(ii,kk,jj,p_nu0) = tmp_minval
          endif
       enddo
       enddo
       enddo

   CASE (CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN)
       CALL wrf_error_fatal3("<stdin>",362,&
'Wet scavenging is currently not possible with MOSAIC unless aqueous aerosols are turned on.')

   CASE (CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
         CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP,SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP )
       CALL wrf_debug(15,'wetscav_driver calling mosaic_wetscav_driver')
       call wetscav_cbmz_mosaic (id,ktau,dtstep,ktauc,config_flags,        &
               dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,        &
               qlsink,precr,preci,precs,precg, qsrflx,                     &
               gas_aqfrac, numgas_aqfrac,                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
   CASE (MOZART_KPP,MOZCART_KPP,T1_MOZCART_KPP,MOZART_MOSAIC_4BIN_KPP)
       CALL wrf_debug(15,'wetscav_driver calling wetscav_mozcart')
       if( config_flags%mp_physics == THOMPSON ) then
         rainrate(:,:,:) = rainprod(:,:,:)
         evaprate(:,:,:) = evapprod(:,:,:)
       elseif( config_flags%mp_physics == WSM6SCHEME ) then
         rainrate(:,:,:) = rainprod(:,:,:)
         evaprate(:,:,:) = evapprod(:,:,:)
       elseif( config_flags%mp_physics == CAMMGMPSCHEME ) then
         rainrate(:,:,:) = prain3d(:,:,:)
         evaprate(:,:,:) = nevapr3d(:,:,:)
       elseif( config_flags%mp_physics == MORR_TWO_MOMENT ) then
         rainrate(:,:,:) = rainprod(:,:,:)
         evaprate(:,:,:) = evapprod(:,:,:)
       else
         rainrate(:,:,:) = 0.
         evaprate(:,:,:) = 0.
       endif
       call wetscav_mozcart( id, ktau, dtstep, ktauc, config_flags,             &
                             dtstepc, t_phy, p8w, t8w, p_phy,                   &
                             chem, rho_phy, cldfra2, rainrate, evaprate,        &
                             qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp,                &
                             gas_aqfrac, numgas_aqfrac, dz8w, dx, dy,           &
                             moist(ims,kms,jms,p_qv), moist(ims,kms,jms,p_qc),  &
                             moist(ims,kms,jms,p_qi), moist(ims,kms,jms,p_qs),  &

                             delta_mass_col,                                    &
                             ids,ide, jds,jde, kds,kde,                         &
                             ims,ime, jms,jme, kms,kme,                         &
                             its,ite, jts,jte, kts,kte                          )

        hno3_col_mdel = delta_mass_col(:,:,p_hno3)
   
   
   CASE (MOZART_MOSAIC_4BIN_AQ_KPP)
       CALL wrf_debug(15,'wetscav_driver calling mosaic_wetscav_driver')
       if( config_flags%mp_physics == THOMPSON ) then
         rainrate(:,:,:) = rainprod(:,:,:)
         evaprate(:,:,:) = evapprod(:,:,:)
       elseif( config_flags%mp_physics == WSM6SCHEME ) then
         rainrate(:,:,:) = rainprod(:,:,:)
         evaprate(:,:,:) = evapprod(:,:,:)
       elseif( config_flags%mp_physics == CAMMGMPSCHEME ) then
         rainrate(:,:,:) = prain3d(:,:,:)
         evaprate(:,:,:) = nevapr3d(:,:,:)
       elseif( config_flags%mp_physics == MORR_TWO_MOMENT ) then
         rainrate(:,:,:) = rainprod(:,:,:)
         evaprate(:,:,:) = evapprod(:,:,:)
       else
         rainrate(:,:,:) = 0.
         evaprate(:,:,:) = 0.
       endif
       call wetscav_mozart_mosaic (id,ktau,dtstep,ktauc,config_flags,        &
               dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,        &
               qlsink,precr,preci,precs,precg, qsrflx,                     &
               gas_aqfrac, numgas_aqfrac,                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
       CALL wrf_debug(15,'wetscav_driver calling wetscav_mozcart')
       call wetscav_mozcart( id, ktau, dtstep, ktauc, config_flags,             &
                             dtstepc, t_phy, p8w, t8w, p_phy,                   &
                             chem, rho_phy, cldfra2, rainrate, evaprate,        &
                             qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp,                &
                             gas_aqfrac, numgas_aqfrac, dz8w, dx, dy,           &
                             moist(ims,kms,jms,p_qv), moist(ims,kms,jms,p_qc),  &
                             moist(ims,kms,jms,p_qi), moist(ims,kms,jms,p_qs),  &

                             delta_mass_col,                                    &
                             ids,ide, jds,jde, kds,kde,                         &
                             ims,ime, jms,jme, kms,kme,                         &
                             its,ite, jts,jte, kts,kte                          )
        hno3_col_mdel = delta_mass_col(:,:,p_hno3)

CASE (CBMZ_CAM_MAM3_NOAQ,CBMZ_CAM_MAM3_AQ,CBMZ_CAM_MAM7_NOAQ,CBMZ_CAM_MAM7_AQ)       
       CALL wrf_debug(15,'wetscav_driver calling wetscav_cam_mam_driver')
       call wetscav_cam_mam_driver (ktau,p_hyd,p8w,t_phy,dgnum4d,          &
            dgnumwet4d,dlf3d,dlf2_3d,dtstep,qme3d,prain3d,nevapr3d,        &
            rate1ord_cw2pr_st3d,shfrc3d,cmfmc,cmfmc2,evapcsh,icwmrsh,      &
            rprdsh,evapcdp3d,icwmrdp3d,rprddp3d,moist(ims,kms,jms,P_QS),   &
            f_ice_phy,f_rain_phy,config_flags,cldfra_mp_all,cldfrai,       &
            cldfral,cldfra,is_CAMMGMP_used,                                &                             
            ids,ide, jds,jde, kds,kde,                                     &
            ims,ime, jms,jme, kms,kme,                                     &
            its,ite, jts,jte, kts,kte,                                     &
            
            moist(ims,kms,jms,P_QV),moist(ims,kms,jms,P_QC),               &
            moist(ims,kms,jms,P_QI),scalar(ims,kms,jms,P_QNI),             &
            scalar(ims,kms,jms,P_QNC),chem,                                &
            
            fracis3D                                                       )
   CASE DEFAULT

   END SELECT cps_select

   
   
   
    SELECT CASE(config_flags%chem_opt)
      
      CASE ( RADM2SORG_AQCHEM,RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP )
      
      do jj=jts,jte
        do ii=its,ite
          
          
          wdi_no3(ii,jj) = - 0.001*qsrflx(ii,jj,p_no3cwj)/mwno3 &
                           - 0.001*qsrflx(ii,jj,p_no3cwi)/mwno3  
          
          
          wd_no3(ii,jj) = wd_no3(ii,jj) + wdi_no3(ii,jj)         
          
          
          wdi_so4(ii,jj) = - 0.001*qsrflx(ii,jj,p_so4cwj)/mwso4 &
                           - 0.001*qsrflx(ii,jj,p_so4cwi)/mwso4 &
                           - qsrflx(ii,jj,p_sulf) &
                           - qsrflx(ii,jj,p_so2)                 
          
          
          wd_so4(ii,jj) = wd_so4(ii,jj) + wdi_so4(ii,jj)         
          
        enddo
      enddo
    
      CASE ( MOZART_MOSAIC_4BIN_AQ_KPP)

      do jj=jts,jte
        do ii=its,ite

          
          wdi_no3(ii,jj) = - 0.001*qsrflx(ii,jj,p_no3_cw01)/mw_no3_aer &
                           - 0.001*qsrflx(ii,jj,p_no3_cw02)/mw_no3_aer &
                           - 0.001*qsrflx(ii,jj,p_no3_cw03)/mw_no3_aer &
                           - 0.001*qsrflx(ii,jj,p_no3_cw04)/mw_no3_aer 

          
          wd_no3(ii,jj) = wd_no3(ii,jj) + wdi_no3(ii,jj)         

          
          wdi_so4(ii,jj) = - 0.001*qsrflx(ii,jj,p_so4_cw01)/mw_so4_aer &
                           - 0.001*qsrflx(ii,jj,p_so4_cw02)/mw_so4_aer &
                           - 0.001*qsrflx(ii,jj,p_so4_cw03)/mw_so4_aer &
                           - 0.001*qsrflx(ii,jj,p_so4_cw04)/mw_so4_aer 

          
          wd_so4(ii,jj) = wd_so4(ii,jj) + wdi_so4(ii,jj)         

          
          wdi_nh4(ii,jj) = - 0.001*qsrflx(ii,jj,p_nh4_cw01)/mw_nh4_aer &
                           - 0.001*qsrflx(ii,jj,p_nh4_cw02)/mw_nh4_aer &
                           - 0.001*qsrflx(ii,jj,p_nh4_cw03)/mw_nh4_aer &
                           - 0.001*qsrflx(ii,jj,p_nh4_cw04)/mw_nh4_aer 

          
          wd_nh4(ii,jj) = wd_nh4(ii,jj) + wdi_nh4(ii,jj)         

          
          wdi_oa(ii,jj)  = - 0.001*qsrflx(ii,jj,p_oc_cw01)/mw_oc_aer &
                           - 0.001*qsrflx(ii,jj,p_oc_cw02)/mw_oc_aer &
                           - 0.001*qsrflx(ii,jj,p_oc_cw03)/mw_oc_aer &
                           - 0.001*qsrflx(ii,jj,p_oc_cw04)/mw_oc_aer &
                           - 0.001*qsrflx(ii,jj,p_asoaX_cw01)/mw_asoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_asoaX_cw02)/mw_asoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_asoaX_cw03)/mw_asoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_asoaX_cw04)/mw_asoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa1_cw01)/mw_asoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa1_cw02)/mw_asoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa1_cw03)/mw_asoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa1_cw04)/mw_asoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa2_cw01)/mw_asoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa2_cw02)/mw_asoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa2_cw03)/mw_asoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa2_cw04)/mw_asoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa3_cw01)/mw_asoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa3_cw02)/mw_asoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa3_cw03)/mw_asoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa3_cw04)/mw_asoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa4_cw01)/mw_asoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa4_cw02)/mw_asoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa4_cw03)/mw_asoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa4_cw04)/mw_asoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoaX_cw01)/mw_bsoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoaX_cw02)/mw_bsoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoaX_cw03)/mw_bsoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoaX_cw04)/mw_bsoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa1_cw01)/mw_bsoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa1_cw02)/mw_bsoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa1_cw03)/mw_bsoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa1_cw04)/mw_bsoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa2_cw01)/mw_bsoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa2_cw02)/mw_bsoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa2_cw03)/mw_bsoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa2_cw04)/mw_bsoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa3_cw01)/mw_bsoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa3_cw02)/mw_bsoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa3_cw03)/mw_bsoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa3_cw04)/mw_bsoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa4_cw01)/mw_bsoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa4_cw02)/mw_bsoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa4_cw03)/mw_bsoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa4_cw04)/mw_bsoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_r1_cw01)/mw_glysoa_r1_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_r1_cw02)/mw_glysoa_r1_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_r1_cw03)/mw_glysoa_r1_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_r1_cw04)/mw_glysoa_r1_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_r2_cw01)/mw_glysoa_r2_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_r2_cw02)/mw_glysoa_r2_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_r2_cw03)/mw_glysoa_r2_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_r2_cw04)/mw_glysoa_r2_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_sfc_cw01)/mw_glysoa_sfc_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_sfc_cw02)/mw_glysoa_sfc_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_sfc_cw03)/mw_glysoa_sfc_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_sfc_cw04)/mw_glysoa_sfc_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_nh4_cw01)/mw_glysoa_nh4_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_nh4_cw02)/mw_glysoa_nh4_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_nh4_cw03)/mw_glysoa_nh4_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_nh4_cw04)/mw_glysoa_nh4_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_oh_cw01)/mw_glysoa_oh_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_oh_cw02)/mw_glysoa_oh_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_oh_cw03)/mw_glysoa_oh_aer &
                           - 0.001*qsrflx(ii,jj,p_glysoa_oh_cw04)/mw_glysoa_oh_aer  

          
          wd_oa(ii,jj) = wd_oa(ii,jj) + wdi_oa(ii,jj)         

          wdi_asoa(ii,jj)= - 0.001*qsrflx(ii,jj,p_asoaX_cw01)/mw_asoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_asoaX_cw02)/mw_asoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_asoaX_cw03)/mw_asoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_asoaX_cw04)/mw_asoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa1_cw01)/mw_asoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa1_cw02)/mw_asoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa1_cw03)/mw_asoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa1_cw04)/mw_asoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa2_cw01)/mw_asoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa2_cw02)/mw_asoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa2_cw03)/mw_asoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa2_cw04)/mw_asoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa3_cw01)/mw_asoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa3_cw02)/mw_asoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa3_cw03)/mw_asoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa3_cw04)/mw_asoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa4_cw01)/mw_asoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa4_cw02)/mw_asoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa4_cw03)/mw_asoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_asoa4_cw04)/mw_asoa4_aer 

          
          wd_asoa(ii,jj) = wd_asoa(ii,jj) + wdi_asoa(ii,jj)         

          wdi_bsoa(ii,jj)= - 0.001*qsrflx(ii,jj,p_bsoaX_cw01)/mw_bsoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoaX_cw02)/mw_bsoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoaX_cw03)/mw_bsoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoaX_cw04)/mw_bsoaX_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa1_cw01)/mw_bsoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa1_cw02)/mw_bsoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa1_cw03)/mw_bsoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa1_cw04)/mw_bsoa1_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa2_cw01)/mw_bsoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa2_cw02)/mw_bsoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa2_cw03)/mw_bsoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa2_cw04)/mw_bsoa2_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa3_cw01)/mw_bsoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa3_cw02)/mw_bsoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa3_cw03)/mw_bsoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa3_cw04)/mw_bsoa3_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa4_cw01)/mw_bsoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa4_cw02)/mw_bsoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa4_cw03)/mw_bsoa4_aer &
                           - 0.001*qsrflx(ii,jj,p_bsoa4_cw04)/mw_bsoa4_aer 

          
          wd_bsoa(ii,jj) = wd_bsoa(ii,jj) + wdi_bsoa(ii,jj)         

          wdi_cvasoa(ii,jj) = - delta_mass_col(ii,jj,p_cvasoaX) &
                              - delta_mass_col(ii,jj,p_cvasoa1) &
                              - delta_mass_col(ii,jj,p_cvasoa2) &
                              - delta_mass_col(ii,jj,p_cvasoa3) &
                              - delta_mass_col(ii,jj,p_cvasoa4) 

          wdi_cvasoa(ii,jj) = wdi_cvasoa(ii,jj) / (dx * dy) / (150.0 * 1e-3) * 1e3 

          
          wd_cvasoa(ii,jj)  = wd_cvasoa(ii,jj) + wdi_cvasoa(ii,jj)         

          wdi_cvbsoa(ii,jj) = - delta_mass_col(ii,jj,p_cvbsoaX) &
                              - delta_mass_col(ii,jj,p_cvbsoa1) &
                              - delta_mass_col(ii,jj,p_cvbsoa2) &
                              - delta_mass_col(ii,jj,p_cvbsoa3) &
                              - delta_mass_col(ii,jj,p_cvbsoa4) 

          wdi_cvbsoa(ii,jj) = wdi_cvbsoa(ii,jj) / (dx * dy) / (180.0 * 1e-3) * 1e3 

          
          wd_cvbsoa(ii,jj)  = wd_cvbsoa(ii,jj) + wdi_cvbsoa(ii,jj)         

          

          wdi_hno3(ii,jj)    = - delta_mass_col(ii,jj,p_hno3) / (63.0123405 * 1e-3) 
          wdi_hno3(ii,jj)    = wdi_hno3(ii,jj) / (dx * dy) * 1e3 
          wd_hno3(ii,jj)     = wd_hno3(ii,jj) + wdi_hno3(ii,jj)               

          wdi_so2(ii,jj)    = - delta_mass_col(ii,jj,p_so2) / (63.961901 * 1e-3) 
          wdi_so2(ii,jj)    = wdi_so2(ii,jj) / (dx * dy) * 1e3 
          wd_so2(ii,jj)     = wd_so2(ii,jj) + wdi_so2(ii,jj)               

          wdi_sulf(ii,jj)    = - delta_mass_col(ii,jj,p_sulf) / (98.078 * 1e-3) 
          wdi_sulf(ii,jj)    = wdi_sulf(ii,jj) / (dx * dy) * 1e3 
          wd_sulf(ii,jj)     = wd_sulf(ii,jj) + wdi_sulf(ii,jj)               

          wdi_nh3(ii,jj)    = - delta_mass_col(ii,jj,p_nh3) / (17.0289402 * 1e-3) 
          wdi_nh3(ii,jj)    = wdi_nh3(ii,jj) / (dx * dy) * 1e3 
          wd_nh3(ii,jj)     = wd_nh3(ii,jj) + wdi_nh3(ii,jj)               

        enddo
      enddo
    

    CASE DEFAULT
      
    END SELECT

    if( config_flags%diagnostic_dep == 1 ) then
       wdflx(its:ite,jts:jte,:)=wdflx(its:ite,jts:jte,:)-qsrflx(its:ite,jts:jte,:)
    end if

   end subroutine wetscav_driver


END MODULE module_wetscav_driver


