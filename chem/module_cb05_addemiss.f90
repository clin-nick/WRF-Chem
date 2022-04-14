




























   MODULE module_cb05_addemiss




























CONTAINS





   subroutine cb05_addemiss_anthro( id, dtstep, dz8w, config_flags,          &
                rho_phy, chem, emis_ant,                                     &
                ids,ide, jds,jde, kds, kde,                                  &
                ims,ime, jms,jme, kms, kme,                                  &
                its,ite, jts,jte, kts, kte )


  USE module_configure
  USE module_state_description
  USE module_data_radm2

   IMPLICIT NONE


   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   REAL,      INTENT(IN   ) ::                                             &
                             dtstep


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem




   REAL, DIMENSION( ims:ime, kms:config_flags%kemit, jms:jme,num_emis_ant),&
         INTENT(IN ) ::                                                    &
                               emis_ant


   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   dz8w, rho_phy


    integer i,j,k
    real, parameter :: efact1 = 1.0/60.0
    real :: conv
   double precision :: chem_sum(num_chem)





         call wrf_debug(15,'cb05_addemiss_anthro')



      do 100 j=jts,jte  
      do 100 i=its,ite  

      DO k=kts,min(config_flags%kemit,kte)


        conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)

        chem(i,k,j,p_no2)  = chem(i,k,j,p_no2)                        &
                           + emis_ant(i,k,j,p_e_no2)*conv
        chem(i,k,j,p_xyl)  = chem(i,k,j,p_xyl)                        &
                           + emis_ant(i,k,j,p_e_xyl)*conv
        chem(i,k,j,p_tol)  = chem(i,k,j,p_tol)                        &
                           + emis_ant(i,k,j,p_e_tol)*conv
        chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                        &
                           + emis_ant(i,k,j,p_e_so2)*conv
        chem(i,k,j,p_no)   = chem(i,k,j,p_no)                         &
                           + emis_ant(i,k,j,p_e_no)*conv
        chem(i,k,j,p_nh3)  = chem(i,k,j,p_nh3)                        &
                           + emis_ant(i,k,j,p_e_nh3)*conv
        chem(i,k,j,p_hcl)  = chem(i,k,j,p_hcl)                        &
                           + emis_ant(i,k,j,p_e_hcl)*conv
        chem(i,k,j,p_co)   = chem(i,k,j,p_co)                         &
                           + emis_ant(i,k,j,p_e_co)*conv
        chem(i,k,j,p_aldx) = chem(i,k,j,p_aldx)                       &
                           + emis_ant(i,k,j,p_e_aldx)*conv

        if (config_flags%bio_emiss_opt == 0) then
            chem(i,k,j,p_terp) = chem(i,k,j,p_terp)                       &
                           + emis_ant(i,k,j,p_e_terp)*conv
        end if


        if ( (config_flags%emiss_opt == 14) ) then
            chem(i,k,j,p_par) = chem(i,k,j,p_par)             &
                + conv*                                       &
                  ( 2.9*emis_ant(i,k,j,p_e_hc3)       &
                  + 4.8*emis_ant(i,k,j,p_e_hc5) + 7.9*emis_ant(i,k,j,p_e_hc8)       &
                  + 0.9*emis_ant(i,k,j,p_e_ket) )
            chem(i,k,j,p_aacd) = chem(i,k,j,p_aacd)                       &
                           + emis_ant(i,k,j,p_e_ora2)*conv
            chem(i,k,j,p_ole)  = chem(i,k,j,p_ole)                        &
                           + emis_ant(i,k,j,p_e_olt)*conv
            chem(i,k,j,p_iole) = chem(i,k,j,p_iole)                       &
                           + emis_ant(i,k,j,p_e_oli)*conv
            chem(i,k,j,p_eth)  = chem(i,k,j,p_eth)                        &
                           + emis_ant(i,k,j,p_e_ol2)*conv
            chem(i,k,j,p_form) = chem(i,k,j,p_form)                       &
                           + emis_ant(i,k,j,p_e_hcho)*conv
            chem(i,k,j,p_etha) = chem(i,k,j,p_etha)                       &
                           + emis_ant(i,k,j,p_e_eth)*conv
            chem(i,k,j,p_cres) = chem(i,k,j,p_cres)                       &
                           + emis_ant(i,k,j,p_e_csl)*conv
            chem(i,k,j,p_meoh) = chem(i,k,j,p_meoh)                       &
                           + emis_ant(i,k,j,p_e_ch3oh)*conv
            chem(i,k,j,p_etoh) = chem(i,k,j,p_etoh)                       &
                           + emis_ant(i,k,j,p_e_c2h5oh)*conv
            chem(i,k,j,p_ald2) = chem(i,k,j,p_ald2)                       &
                           + emis_ant(i,k,j,p_e_ald)*conv

            if (config_flags%bio_emiss_opt == 0) then
             chem(i,k,j,p_isop) = chem(i,k,j,p_isop)                       &
                           + emis_ant(i,k,j,p_e_iso)*conv
            end if


        else
            chem(i,k,j,p_par) = chem(i,k,j,p_par)             &
                + conv*emis_ant(i,k,j,p_e_par)
            chem(i,k,j,p_ole)  = chem(i,k,j,p_ole)                        &
                           + emis_ant(i,k,j,p_e_ole)*conv
            chem(i,k,j,p_iole) = chem(i,k,j,p_iole)                       &
                           + emis_ant(i,k,j,p_e_iole)*conv
            chem(i,k,j,p_eth)  = chem(i,k,j,p_eth)                        &
                           + emis_ant(i,k,j,p_e_eth)*conv
            chem(i,k,j,p_form) = chem(i,k,j,p_form)                       &
                           + emis_ant(i,k,j,p_e_form)*conv
            chem(i,k,j,p_etha) = chem(i,k,j,p_etha)                       &
                           + emis_ant(i,k,j,p_e_etha)*conv
            chem(i,k,j,p_cres) = chem(i,k,j,p_cres)                       &
                           + emis_ant(i,k,j,p_e_cres)*conv                &
                           + emis_ant(i,k,j,p_e_phen)*conv 
            chem(i,k,j,p_meoh) = chem(i,k,j,p_meoh)                       &
                           + emis_ant(i,k,j,p_e_meoh)*conv
            chem(i,k,j,p_etoh) = chem(i,k,j,p_etoh)                       &
                           + emis_ant(i,k,j,p_e_etoh)*conv
            chem(i,k,j,p_ald2) = chem(i,k,j,p_ald2)                       &
                           + emis_ant(i,k,j,p_e_ald2)*conv
            chem(i,k,j,p_meo2) = chem(i,k,j,p_meo2)                       &
                           + emis_ant(i,k,j,p_e_meo2)*conv
            chem(i,k,j,p_sulf) = chem(i,k,j,p_sulf)                       &
                           + emis_ant(i,k,j,p_e_psulf)*conv
            chem(i,k,j,p_mgly) = chem(i,k,j,p_mgly)                       &
                           + emis_ant(i,k,j,p_e_mgly)*conv
            chem(i,k,j,p_facd) = chem(i,k,j,p_facd)                       &
                           + emis_ant(i,k,j,p_e_hcooh)*conv
            chem(i,k,j,p_aacd) = chem(i,k,j,p_aacd)                       &
                           + emis_ant(i,k,j,p_e_ccooh)*conv
            chem(i,k,j,p_ispd) = chem(i,k,j,p_ispd)                       &
                           + emis_ant(i,k,j,p_e_iprod)*conv



            if (config_flags%bio_emiss_opt == 0) then
             chem(i,k,j,p_isop) = chem(i,k,j,p_isop)                       &
                           + emis_ant(i,k,j,p_e_isop)*conv
            end if
        end if

      END DO                                                          
 100  continue

    END subroutine cb05_addemiss_anthro





  subroutine cb05_addemiss_bio( id, dtstep, dz8w, config_flags,       &
        rho_phy, chem, e_bio, ne_area, e_iso,                         &
        ids,ide, jds,jde, kds,kde,                                    &
        ims,ime, jms,jme, kms,kme,                                    &
        its,ite, jts,jte, kts,kte                                     )

  USE module_configure
  USE module_state_description
  USE module_data_radm2
  USE module_aerosols_sorgam

  IMPLICIT NONE


   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id, ne_area,                             &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL,      INTENT(IN   ) ::    dtstep

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::        chem

   REAL, DIMENSION( ims:ime, jms:jme,ne_area ),                            &
         INTENT(IN ) ::           e_bio
         
   REAL, DIMENSION( ims:ime, kms:config_flags%kemit, jms:jme ),            &
         INTENT(IN ) ::           e_iso

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::        dz8w, rho_phy            



   integer i,j,k,n
   real, parameter :: efact1 = 1.0/60.0
   double precision :: chem_sum(num_chem)






   if (config_flags%bio_emiss_opt == GUNTHER1) then

      do j=jts,jte  
      do i=its,ite  
        chem(i,kts,j,p_isop) = chem(i,kts,j,p_isop)    &
                          + e_bio(i,j,liso)/(dz8w(i,kts,j)*60.)*dtstep

        chem(i,kts,j,p_terp) = chem(i,kts,j,p_terp)    &
                          + e_bio(i,j,ltpan)/(dz8w(i,kts,j)*60.)*dtstep
      end do
      end do

   end if


   END subroutine cb05_addemiss_bio


END MODULE module_cb05_addemiss
