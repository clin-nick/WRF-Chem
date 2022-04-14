

MODULE module_cbm4_addemiss


   integer, parameter :: cbm4_addemiss_masscheck = -1
                       



CONTAINS




   subroutine cbm4_addemiss_anthro( id, dtstep, dz8w, config_flags,       &
               rho_phy, chem,emis_ant,                                    &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )






  USE module_configure
  USE module_state_description
  USE module_data_radm2

  IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem




   REAL, DIMENSION( ims:ime, kms:config_flags%kemit, jms:jme,num_emis_ant),&
         INTENT(IN ) ::                                                    &
                         emis_ant

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   dz8w, rho_phy


   integer :: i,j,k
   real, parameter :: efact1 = 1.0/60.0
   real :: conv
   double precision :: chem_sum(num_chem)



























      do 100 j=jts,jte  
      do 100 i=its,ite  

      DO k=kts,min(config_flags%kemit,kte)


        conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)




























        chem(i,k,j,p_cres)  =  chem(i,k,j,p_cres)                       &
                         +emis_ant(i,k,j,p_e_csl)  * conv 
        chem(i,k,j,p_no)    = chem(i,k,j,p_no)                          &
                         +emis_ant(i,k,j,p_e_no)   * conv
        chem(i,k,j,p_ald2)  = chem(i,k,j,p_ald2)                        &
                         +emis_ant(i,k,j,p_e_ald)  * conv
        chem(i,k,j,p_hcho)  = chem(i,k,j,p_hcho)                        &
                         +emis_ant(i,k,j,p_e_hcho) * conv
        chem(i,k,j,p_eth)   = chem(i,k,j,p_eth)                         &
                         +emis_ant(i,k,j,p_e_ol2)  * conv
        chem(i,k,j,p_co)    = chem(i,k,j,p_co)                          &
                         +emis_ant(i,k,j,p_e_co)   * conv
        chem(i,k,j,p_tol)   = chem(i,k,j,p_tol)                         &
                         +emis_ant(i,k,j,p_e_tol)  * conv
        chem(i,k,j,p_xyl)   = chem(i,k,j,p_xyl)                         &
                         +emis_ant(i,k,j,p_e_xyl)  * conv       



        if ( (config_flags%emiss_inpt_opt == EMISS_INPT_DEFAULT) .or.   &
             (config_flags%emiss_inpt_opt == EMISS_INPT_PNNL_RS) ) then
            chem(i,k,j,p_par) = chem(i,k,j,p_par)             &
                + conv*                                       &
                  ( 0.4*emis_ant(i,k,j,p_e_ald) + 2.9*emis_ant(i,k,j,p_e_hc3)       &
                  + 4.8*emis_ant(i,k,j,p_e_hc5) + 7.9*emis_ant(i,k,j,p_e_hc8)       &
                  + 0.9*emis_ant(i,k,j,p_e_ket) + 2.8*emis_ant(i,k,j,p_e_oli)       &
                  + 1.8*emis_ant(i,k,j,p_e_olt) + 1.0*emis_ant(i,k,j,p_e_ora2) )
            chem(i,k,j,p_ole)  = chem(i,k,j,p_ole)                         &
             +(emis_ant(i,k,j,p_e_oli)+emis_ant(i,k,j,p_e_olt))*conv




        elseif(config_flags%emiss_inpt_opt == EMISS_INPT_CB4) then
            chem(i,k,j,p_par)  = chem(i,k,j,p_par)             &
                + conv*emis_ant(i,k,j,p_e_hc5)
            chem(i,k,j,p_no2)  = chem(i,k,j,p_no2)             &
                + conv*emis_ant(i,k,j,p_e_no2)
            chem(i,k,j,p_ole)  = chem(i,k,j,p_ole)             &
                 +emis_ant(i,k,j,p_e_oli)*conv
        end if

      END DO                                                          
 100  continue























   END subroutine cbm4_addemiss_anthro




  subroutine cbm4_addemiss_bio( id, dtstep, dz8w, config_flags,       &
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
        chem(i,kts,j,p_so2) = chem(i,kts,j,p_so2)    &
                          + e_bio(i,j,lso2)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_sulf) = chem(i,kts,j,p_sulf)    &
                          + e_bio(i,j,lsulf)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_no2) = chem(i,kts,j,p_no2)    &
                          + e_bio(i,j,lno2)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_no) = chem(i,kts,j,p_no)    &
                          + e_bio(i,j,lno)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_o3) = chem(i,kts,j,p_o3)    &
                          + e_bio(i,j,lo3)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_hno3) = chem(i,kts,j,p_hno3)    &
                          + e_bio(i,j,lhno3)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_h2o2) = chem(i,kts,j,p_h2o2)    &
                          + e_bio(i,j,lh2o2)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_ald2) = chem(i,kts,j,p_ald2)    &
                          + e_bio(i,j,lald)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_hcho) = chem(i,kts,j,p_hcho)    &
                          + e_bio(i,j,lhcho)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_no2) = chem(i,kts,j,p_no2)    &
                          + e_bio(i,j,lno2)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_pan) = chem(i,kts,j,p_pan)    &
                          + e_bio(i,j,lpan)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_co) = chem(i,kts,j,p_co)    &
                          + e_bio(i,j,lco)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_eth) = chem(i,kts,j,p_eth)    &
                          + e_bio(i,j,lol2)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_ole) = chem(i,kts,j,p_ole)    &
                          + (e_bio(i,j,lolt)+ e_bio(i,j,loli))  &
                          /(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_tol) = chem(i,kts,j,p_tol)    &
                          + e_bio(i,j,ltol)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_xyl) = chem(i,kts,j,p_xyl)    &
                          + e_bio(i,j,lxyl)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_onit) = chem(i,kts,j,p_onit)    &
                          + e_bio(i,j,lonit)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_csl) = chem(i,kts,j,p_csl)    &
                          + e_bio(i,j,lcsl)/(dz8w(i,kts,j)*60.)*dtstep
        chem(i,kts,j,p_iso) = chem(i,kts,j,p_iso)    &
                          + e_bio(i,j,liso)/(dz8w(i,kts,j)*60.)*dtstep
      end do
      end do



      do j = jts, jte
      do i = its, ite
         chem(i,kts,j,p_par)  =  chem(i,kts,j,p_par)               &
             + (dtstep/(dz8w(i,kts,j)*60.))*                       &
               ( 0.4*e_bio(i,j,lald) + 2.9*e_bio(i,j,lhc3)         &
               + 4.8*e_bio(i,j,lhc5) + 7.9*e_bio(i,j,lhc8)         &
               + 0.9*e_bio(i,j,lket) + 2.8*e_bio(i,j,loli)         &
               + 1.8*e_bio(i,j,lolt) + 1.0*e_bio(i,j,lora2)        )
      end do
      end do




















      end if





   if (config_flags%bio_emiss_opt /= GUNTHER1) then












      do j = jts, jte
      do k = kts, min(config_flags%kemit,kte)
      do i = its, ite
         chem(i,k,j,p_iso) = chem(i,k,j,p_iso) + e_iso(i,k,j)              &
              *4.828e-4/rho_phy(i,k,j)*(dtstep/(dz8w(i,k,j)*60.))
      end do
      end do
      end do












   end if


   END subroutine cbm4_addemiss_bio


END MODULE module_cbm4_addemiss



