













      SUBROUTINE optical_driver(id,curr_secs,dtstep,config_flags,haveaer,&
               chem,dz8w,alt,relhum,                                     &
               h2oai,h2oaj,                                              &
               tauaer1,tauaer2,tauaer3,tauaer4,                          &
               
               extaer1,extaer2,extaer3,extaer4,                          &
               gaer1,gaer2,gaer3,gaer4,                                  &
               waer1,waer2,waer3,waer4,                                  &
               bscoef1,bscoef2,bscoef3,bscoef4,                          &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                      &
               totoa_a01,totoa_a02,totoa_a03,totoa_a04,                  &
               totoa_a05,totoa_a06,totoa_a07,totoa_a08,                  & 
               extaerlw1,extaerlw2,extaerlw3,extaerlw4,extaerlw5,extaerlw6, &
               extaerlw7,extaerlw8,extaerlw9,extaerlw10,extaerlw11,extaerlw12, &
               extaerlw13,extaerlw14,extaerlw15,extaerlw16,  &
               tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,tauaerlw5,tauaerlw6, & 
               tauaerlw7,tauaerlw8,tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12, & 
               tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16,  & 
               ids,ide, jds,jde, kds,kde,                                &
               ims,ime, jms,jme, kms,kme,                                &
               its,ite, jts,jte, kts,kte                                 )


   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE module_optical_averaging
   USE module_data_mosaic_therm, only: nbin_a
   USE module_data_rrtmgaeropt, only: nswbands,nlwbands 
   USE module_peg_util, only:  peg_error_fatal, peg_message
   use infnan,                 only: inf
   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) :: id,                                  &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
   REAL(KIND=8), INTENT(IN   ) :: curr_secs
   REAL,         INTENT(IN   ) :: dtstep



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(IN ) ::  chem

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(IN ) ::  relhum, dz8w, alt, h2oai, h2oaj,              &
                         totoa_a01, totoa_a02, totoa_a03, totoa_a04,   &
                         totoa_a05,totoa_a06,totoa_a07,totoa_a08



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(INOUT ) ::                                             &
           tauaer1, tauaer2, tauaer3, tauaer4,                         &
           
           extaer1, extaer2, extaer3, extaer4,                         &
           gaer1, gaer2, gaer3, gaer4,                                 &
           waer1, waer2, waer3, waer4,                                 &
           bscoef1, bscoef2, bscoef3, bscoef4                              
   
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme),                &
         INTENT(INOUT ) :: extaerlw1,extaerlw2,extaerlw3,extaerlw4,extaerlw5, & 
                           extaerlw6,extaerlw7,extaerlw8,extaerlw9,extaerlw10, & 
                           extaerlw11,extaerlw12,extaerlw13,extaerlw14,extaerlw15, &
                           extaerlw16,   & 
                           tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,tauaerlw5, &
                           tauaerlw6,tauaerlw7,tauaerlw8,tauaerlw9,tauaerlw10, & 
                           tauaerlw11,tauaerlw12,tauaerlw13,tauaerlw14,tauaerlw15, &
                           tauaerlw16
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:4) ::  & 
         tauaersw,extaersw,gaersw,waersw,bscoefsw 
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:4 ),                  &  
         INTENT(INOUT ) ::                                             &
           l2aer, l3aer, l4aer, l5aer, l6aer, l7aer

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:nlwbands) ::  &
         extaerlw,tauaerlw 

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   character*100 msg
   integer lunerr 

   LOGICAL, INTENT(IN) :: haveaer



      logical processingAerosols
      integer nbin_o
      integer option_method, option_mie








   select case (config_flags%chem_opt)
   case ( RADM2SORG,           RADM2SORG_KPP,      RADM2SORG_AQ, RADM2SORG_AQCHEM, &
          GOCART_SIMPLE,       RACMSORG_KPP,       RACMSORG_AQ,  RACMSORG_AQCHEM_KPP, &
          RACM_ESRLSORG_AQCHEM_KPP, RACM_SOA_VBS_KPP, RACM_SOA_VBS_AQCHEM_KPP,        &
          RACM_SOA_VBS_HET_KPP,        &
          GOCARTRACM_KPP,      GOCARTRADM2,  &
          RACM_ESRLSORG_KPP,   MOZCART_KPP,        T1_MOZCART_KPP,  &
          CBMZ_MOSAIC_4BIN,    CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_KPP,   &
          CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, CBMZSORG, CBMZSORG_AQ, &
          CBMZ_MOSAIC_DMS_4BIN,    CBMZ_MOSAIC_DMS_8BIN,   &
          CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &    
          SAPRC99_MOSAIC_4BIN_VBS2_KPP, &
          MOZART_MOSAIC_4BIN_KPP , MOZART_MOSAIC_4BIN_AQ_KPP, &
          CBMZ_CAM_MAM3_NOAQ,CBMZ_CAM_MAM7_NOAQ,  CBMZ_CAM_MAM3_AQ,  &
          CBMZ_CAM_MAM7_AQ, CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP, &
          SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_KPP,  &
          CB05_SORG_AQ_KPP, CB05_SORG_VBS_AQ_KPP )
      processingAerosols = .true.
      call wrf_debug(15,'optical driver: process aerosols true')
   case default
      processingAerosols = .false.
      call wrf_debug(15,'optical driver: process aerosols false')
   end select

  if( processingAerosols ) then






   select case (config_flags%chem_opt)
   case ( RADM2SORG, RACM_ESRLSORG_KPP, RADM2SORG_KPP, RADM2SORG_AQ, RADM2SORG_AQCHEM, &
          GOCARTRACM_KPP,      GOCARTRADM2,      &
          GOCART_SIMPLE,       RACMSORG_KPP,       RACMSORG_AQ,      RACMSORG_AQCHEM_KPP, &
          RACM_ESRLSORG_AQCHEM_KPP, RACM_SOA_VBS_KPP, RACM_SOA_VBS_AQCHEM_KPP,            &
          RACM_SOA_VBS_HET_KPP, CBMZSORG, CBMZSORG_AQ, MOZCART_KPP, T1_MOZCART_KPP,     &
          CBMZ_CAM_MAM3_NOAQ,  CBMZ_CAM_MAM7_NOAQ,  CBMZ_CAM_MAM3_AQ,  CBMZ_CAM_MAM7_AQ, &
          CB05_SORG_AQ_KPP, CB05_SORG_VBS_AQ_KPP )
     nbin_o = 8
   case (CBMZ_MOSAIC_4BIN,    CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_KPP,  &
         CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
         CBMZ_MOSAIC_DMS_4BIN,    CBMZ_MOSAIC_DMS_8BIN,   &
         CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &    
         SAPRC99_MOSAIC_4BIN_VBS2_KPP, &
         MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP, &
         CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP, &
         SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_KPP )
     nbin_o = nbin_a
   end select

     call wrf_debug(15,'optical averaging')
     aer_op_opt_select: SELECT CASE(config_flags%aer_op_opt)
     CASE (VOLUME_APPROX)
       option_method=1
       option_mie=1
     CASE (MAXWELL_APPROX)
       option_method=2
       option_mie=1
     CASE (VOLUME_EXACT)
       option_method=1
       option_mie=2
     CASE (MAXWELL_EXACT)
       option_method=2
       option_mie=2
     CASE (SHELL_EXACT)
       option_method=3
       option_mie=2
     CASE DEFAULT
        if( config_flags%aer_op_opt > 0 ) then
           call wrf_message('WARNING: Invalid aer_op_opt. Defaulting to VOLUME_APPROX.')
           option_method=1
           option_mie=1
        end if
     END SELECT aer_op_opt_select

     if( config_flags%aer_op_opt > 0 ) then
        call wrf_debug(15,'optical driver: call optical averaging')











        
        tauaersw(:,:,:,:) = inf 
        extaersw(:,:,:,:) = inf
        gaersw(:,:,:,:)   = inf
        waersw(:,:,:,:)   = inf
        bscoefsw(:,:,:,:) = inf
        
        extaerlw(:,:,:,:) = inf
        tauaerlw(:,:,:,:) = inf

        call optical_averaging(id,curr_secs,dtstep,config_flags,     &
             nbin_o,haveaer,option_method,option_mie,chem,dz8w,alt,  &
             relhum,h2oai,h2oaj,                                     &




             tauaersw,extaersw,gaersw,waersw,bscoefsw,               &
             l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                    &
             totoa_a01,totoa_a02,totoa_a03,totoa_a04,                &
             totoa_a05,totoa_a06,totoa_a07,totoa_a08,                &
             tauaerlw,extaerlw,                                      &
             ids,ide, jds,jde, kds,kde,                              &
             ims,ime, jms,jme, kms,kme,                              &
             its,ite, jts,jte, kts,kte                               )
             
             tauaer1=tauaersw(:,:,:,1)
             tauaer2=tauaersw(:,:,:,2)
             tauaer3=tauaersw(:,:,:,3)
             tauaer4=tauaersw(:,:,:,4)
             extaer1=extaersw(:,:,:,1)
             extaer2=extaersw(:,:,:,2)
             extaer3=extaersw(:,:,:,3)
             extaer4=extaersw(:,:,:,4)
             gaer1=gaersw(:,:,:,1)
             gaer2=gaersw(:,:,:,2)
             gaer3=gaersw(:,:,:,3)
             gaer4=gaersw(:,:,:,4)
             waer1=waersw(:,:,:,1)
             waer2=waersw(:,:,:,2)
             waer3=waersw(:,:,:,3)
             waer4=waersw(:,:,:,4)
             bscoef1=bscoefsw(:,:,:,1)
             bscoef2=bscoefsw(:,:,:,2)
             bscoef3=bscoefsw(:,:,:,3)
             bscoef4=bscoefsw(:,:,:,4)
             
             extaerlw1=extaerlw(:,:,:,1)
             extaerlw2=extaerlw(:,:,:,2)
             extaerlw3=extaerlw(:,:,:,3)
             extaerlw4=extaerlw(:,:,:,4)
             extaerlw5=extaerlw(:,:,:,5)
             extaerlw6=extaerlw(:,:,:,6)
             extaerlw7=extaerlw(:,:,:,7)
             extaerlw8=extaerlw(:,:,:,8)
             extaerlw9=extaerlw(:,:,:,9)
             extaerlw10=extaerlw(:,:,:,10)
             extaerlw11=extaerlw(:,:,:,11)
             extaerlw12=extaerlw(:,:,:,12)
             extaerlw13=extaerlw(:,:,:,13)
             extaerlw14=extaerlw(:,:,:,14)
             extaerlw15=extaerlw(:,:,:,15)
             extaerlw16=extaerlw(:,:,:,16)
             tauaerlw1=tauaerlw(:,:,:,1)
             tauaerlw2=tauaerlw(:,:,:,2)
             tauaerlw3=tauaerlw(:,:,:,3)
             tauaerlw4=tauaerlw(:,:,:,4)
             tauaerlw5=tauaerlw(:,:,:,5)
             tauaerlw6=tauaerlw(:,:,:,6)
             tauaerlw7=tauaerlw(:,:,:,7)
             tauaerlw8=tauaerlw(:,:,:,8)
             tauaerlw9=tauaerlw(:,:,:,9)
             tauaerlw10=tauaerlw(:,:,:,10)
             tauaerlw11=tauaerlw(:,:,:,11)
             tauaerlw12=tauaerlw(:,:,:,12)
             tauaerlw13=tauaerlw(:,:,:,13)
             tauaerlw14=tauaerlw(:,:,:,14)
             tauaerlw15=tauaerlw(:,:,:,15)
             tauaerlw16=tauaerlw(:,:,:,16)
        call wrf_debug(15,'optical driver: after call optical averaging')
     else
        
        
        
     end if

   endif
   return

END SUBROUTINE optical_driver
