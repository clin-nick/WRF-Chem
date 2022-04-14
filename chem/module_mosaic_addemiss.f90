







MODULE module_mosaic_addemiss









   integer, parameter :: mosaic_addemiss_active = 1
                       
                       

   integer, parameter :: mosaic_addemiss_masscheck = -1
                       

CONTAINS




   subroutine mosaic_addemiss( id, dtstep, u10, v10, alt, dz8w, xland,     &
               config_flags, chem, slai, ust, smois, ivgtyp, isltyp,       &
               emis_ant,ebu,biom_active,dust_opt,                          &
               
               ktau,p8w, u_phy,v_phy,rho_phy,g,dx,erod,  &  
               dust_emiss_active, seasalt_emiss_active,                    &
               dust_flux, seas_flux,                                       &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )






   USE module_configure, only:  grid_config_rec_type
   USE module_state_description


   USE module_data_mosaic_asect

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   INTEGER, INTENT(IN) ::    dust_emiss_active, seasalt_emiss_active, biom_active, dust_opt

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: dust_flux, seas_flux

   
   INTEGER, INTENT(IN) :: ktau
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )  ,               &
          INTENT(IN) ::  p8w, u_phy,v_phy,rho_phy
   REAL, INTENT(IN   ) ::  dx, g
   REAL, DIMENSION( ims:ime, jms:jme,3),&
         INTENT(IN ) ::     erod


   REAL, INTENT(IN   ) ::    dtstep


   REAL,  DIMENSION( ims:ime , jms:jme )         ,                         &
          INTENT(IN   ) ::   u10, v10, xland, slai, ust
   INTEGER,  DIMENSION( ims:ime , jms:jme )      ,                         &
          INTENT(IN   ) ::   ivgtyp, isltyp


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem




   REAL, DIMENSION( ims:ime,   1:config_flags%kemit, jms:jme,num_emis_ant ),            &
         INTENT(IN ) ::                                                    &
                         emis_ant

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,num_ebu ),             &
         INTENT(IN    ) ::                                             &
         ebu


   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w

   REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,      &
          INTENT(INOUT) ::   smois


   integer i, j, k, l, n
   integer iphase, itype
   integer p1st

   real, parameter :: efact1 = 1.0
   real :: aem_so4, aem_no3, aem_cl, aem_msa, aem_co3, aem_nh4,   &
           aem_na, aem_ca, aem_oin, aem_oc, aem_bc, aem_num
   real dum, fact







   real, save :: fr8b_aem_sorgam_i(8) =   & 
                 (/ 0.965,  0.035,  0.000,  0.000,   &
                    0.000,  0.000,  0.000,  0.000 /)



   real, save :: fr8b_aem_sorgam_j(8) =   &
                 (/ 0.026,  0.147,  0.350,  0.332,   &
                    0.125,  0.019,  0.001,  0.000/)




   real, save :: fr8b_aem_sorgam_c(8) =   &
                 (/ 0.000,  0.000,  0.000,  0.002,   &
                    0.021,  0.110,  0.275,  0.592 /)









   real, save :: fr8b_aem_mosaic_f(8) =   &
       (/ 0.060, 0.045, 0.245, 0.400, 0.100, 0.150, 0., 0./) 







   real, save :: fr8b_aem_mosaic_c(8) =   &
       (/ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.300, 0.700 /) 







real, save :: fr8b_bburn_mosaic(8) = &
      (/ 0.0092, 0.1385, 0.4548, 0.3388, 0.0567, 0.002, 0.0, 0.0/) 





   real, save :: fr8b_tno_bc1(8) =   &
       (/ 0.0494,  0.3795,  0.4714,  0.0967,  0.003,  0.0,  0.0, 0.0 /)



   real, save :: fr8b_tno_ec25(8) =   &
       (/ 0.0,  0.0,  0.0,  0.0, 0.40, 0.60,  0.0,  0.0 /)



   real, save :: fr8b_tno_oc_dom(8) =   &
       (/ 0.0358, 0.1325, 0.2704, 0.3502, 0.1904, 0.0657, 0.0, 0.0 /)



   real, save :: fr8b_tno_oc_tra(8) =   &
       (/ 0.0063, 0.0877, 0.3496, 0.4054, 0.1376, 0.0134, 0.0, 0.0 /)



   real, save :: fr8b_tno_mosaic_f(8) =   &
       (/ 0.060, 0.045, 0.245, 0.400, 0.100, 0.150, 0.0, 0.0/)



   real, save :: fr8b_tno_mosaic_c(8) =   &
       (/ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.300, 0.700 /)









   real :: fr_aem_sorgam_i(8)
   real :: fr_aem_sorgam_j(8)
   real :: fr_aem_sorgam_c(8)
   real :: fr_aem_mosaic_f(8)
   real :: fr_aem_mosaic_c(8)
   real :: fr_tno_bc1(8)
   real :: fr_tno_ec25(8)
   real :: fr_tno_oc_dom(8)
   real :: fr_tno_oc_tra(8)
   real :: fr_tno_mosaic_f(8)
   real :: fr_tno_mosaic_c(8)

   real :: fr_aem_gc2mosaic_f(8) 
   real :: fr_aem_gc2mosaic_c(8)
   real :: bburn_mosaic_f(8)       
   real :: bburn_mosaic_c(8)       

   double precision :: chem_sum(num_chem)

   character(len=80) :: msg



   itype = 1
   iphase = ai_phase

























        if ((nsize_aer(itype) .ne. 4) .and. (nsize_aer(itype) .ne. 8)) then
          write(msg,'(a,i5)')   &
           'subr mosaic_addemiss - nsize_aer(itype) must be ' //   &
           '4 or 8 but = ',  nsize_aer(itype)
          call wrf_error_fatal3("<stdin>",261,&
msg )
        end if

        fr_aem_sorgam_i(:) = 0.0
        fr_aem_sorgam_j(:) = 0.0
        fr_aem_sorgam_c(:) = 0.0
        fr_aem_mosaic_f(:) = 0.0
        fr_aem_mosaic_c(:) = 0.0
        fr_tno_bc1(:) = 0.0
        fr_tno_ec25(:) = 0.0
        fr_tno_oc_dom(:) = 0.0
        fr_tno_oc_tra(:) = 0.0
        fr_tno_mosaic_f(:) = 0.0
        fr_tno_mosaic_c(:) = 0.0
        fr_aem_gc2mosaic_f(:) = 0.0
        fr_aem_gc2mosaic_c(:) = 0.0

	emiss_inpt_select_1: SELECT CASE( config_flags%emiss_inpt_opt )

	  CASE( emiss_inpt_default, emiss_inpt_pnnl_rs )
	    if (nsize_aer(itype) .eq. 8) then
		fr_aem_sorgam_i(:) = fr8b_aem_sorgam_i(:)
		fr_aem_sorgam_j(:) = fr8b_aem_sorgam_j(:)
		fr_aem_sorgam_c(:) = fr8b_aem_sorgam_c(:)
	    else if (nsize_aer(itype) .eq. 4) then
		do n = 1, nsize_aer(itype)
		    fr_aem_sorgam_i(n) = fr8b_aem_sorgam_i(2*n-1)   &
		                       + fr8b_aem_sorgam_i(2*n)
		    fr_aem_sorgam_j(n) = fr8b_aem_sorgam_j(2*n-1)   &
		                       + fr8b_aem_sorgam_j(2*n)
		    fr_aem_sorgam_c(n) = fr8b_aem_sorgam_c(2*n-1)   &
		                       + fr8b_aem_sorgam_c(2*n)
		end do
	    end if








            if( config_flags%emiss_opt == 5 .or. config_flags%emiss_opt == 6 ) then
              CALL wrf_debug(15,'mosaic_addemiss: emiss_opt = eptec')
              CALL wrf_debug(15,'mosaic_addemiss: gocart speciation being mapped to mosaic')
  
              fr_aem_sorgam_i(:) = 0.0
              fr_aem_sorgam_j(:) = 0.0
              fr_aem_sorgam_c(:) = 0.0


              if (nsize_aer(itype) .eq. 8) then
                fr_aem_gc2mosaic_f(:) = fr8b_aem_mosaic_f(:)
                fr_aem_gc2mosaic_c(:) = fr8b_aem_mosaic_c(:)
              elseif (nsize_aer(itype) .eq. 4) then
                do n = 1, nsize_aer(itype)
                  fr_aem_gc2mosaic_f(n) = fr8b_aem_mosaic_f(2*n-1) + fr8b_aem_mosaic_f(2*n)
                  fr_aem_gc2mosaic_c(n) = fr8b_aem_mosaic_c(2*n-1) + fr8b_aem_mosaic_c(2*n)
                end do
              endif
            endif

          CASE( emiss_inpt_pnnl_cm )
	    if (nsize_aer(itype) .eq. 8) then
		fr_aem_mosaic_f(:) = fr8b_aem_mosaic_f(:)
		fr_aem_mosaic_c(:) = fr8b_aem_mosaic_c(:)
	    else if (nsize_aer(itype) .eq. 4) then
		do n = 1, nsize_aer(itype)
		    fr_aem_mosaic_f(n) = fr8b_aem_mosaic_f(2*n-1)   &
		                       + fr8b_aem_mosaic_f(2*n)
		    fr_aem_mosaic_c(n) = fr8b_aem_mosaic_c(2*n-1)   &
		                       + fr8b_aem_mosaic_c(2*n)
		end do
	    end if

	  CASE( emiss_inpt_tno )		
	    if (nsize_aer(itype) .eq. 8) then
		fr_tno_bc1(:) = fr8b_tno_bc1(:)
		fr_tno_ec25(:) = fr8b_tno_ec25(:)
		fr_tno_oc_dom(:) = fr8b_tno_oc_dom(:)
		fr_tno_oc_tra(:) = fr8b_tno_oc_tra(:)
        fr_tno_mosaic_c(:) = fr8b_tno_mosaic_c(:)
        fr_tno_mosaic_f(:) = fr8b_tno_mosaic_f(:)
	    else if (nsize_aer(itype) .eq. 4) then
		do n = 1, nsize_aer(itype)
		    fr_tno_bc1(n) = fr8b_tno_bc1(2*n-1)   &
		                  + fr8b_tno_bc1(2*n)
		    fr_tno_ec25(n) = fr8b_tno_ec25(2*n-1)   &
		                   + fr8b_tno_ec25(2*n)
		    fr_tno_oc_dom(n) = fr8b_tno_oc_dom(2*n-1)   &
		                     + fr8b_tno_oc_dom(2*n)
		    fr_tno_oc_tra(n) = fr8b_tno_oc_tra(2*n-1)   &
		                     + fr8b_tno_oc_tra(2*n)
		    fr_tno_mosaic_c(n) = fr8b_tno_mosaic_c(2*n-1)   &
		                       + fr8b_tno_mosaic_c(2*n)
		    fr_tno_mosaic_f(n) = fr8b_tno_mosaic_f(2*n-1)   &
		                       + fr8b_tno_mosaic_f(2*n)
		end do
	    end if

	  CASE DEFAULT
	    return

	END SELECT emiss_inpt_select_1



        fire_inpt_select: SELECT CASE (config_flags%biomass_burn_opt)
          CASE (BIOMASSB,BIOMASSB_MOZC)
            if (nsize_aer(itype) .eq. 8) then
              bburn_mosaic_f(:) = fr8b_bburn_mosaic(:)

              bburn_mosaic_c(:) = fr8b_aem_mosaic_c(:)

            else if (nsize_aer(itype) .eq. 4) then
              do n = 1, nsize_aer(itype)
                bburn_mosaic_f(n) = fr8b_bburn_mosaic(2*n-1) + fr8b_bburn_mosaic(2*n)

                bburn_mosaic_c(n) = fr8b_aem_mosaic_c(2*n-1) + fr8b_aem_mosaic_c(2*n)

              end do
          end if
        END SELECT fire_inpt_select



        if (mosaic_addemiss_active <= 0 .and. biom_active <= 0) then
	    fr_aem_sorgam_i(:) = 0.0
	    fr_aem_sorgam_j(:) = 0.0
	    fr_aem_sorgam_c(:) = 0.0
	    fr_aem_mosaic_f(:) = 0.0
	    fr_aem_mosaic_c(:) = 0.0
	    fr_tno_bc1(:) = 0.0
	    fr_tno_ec25(:) = 0.0
	    fr_tno_oc_dom(:) = 0.0
	    fr_tno_oc_tra(:) = 0.0
	    fr_tno_mosaic_f(:) = 0.0
	    fr_tno_mosaic_c(:) = 0.0
	end if



	if (mosaic_addemiss_masscheck > 0) call addemiss_masscheck(            &
               id, config_flags, 1, 'mosaic_ademiss',                      &
               dtstep, efact1, dz8w, chem, chem_sum,                       &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte,                                  &
               14,                                                         &
               emis_ant(ims,  1,jms,p_e_pm_10),emis_ant(ims,  1,jms,p_e_pm_25), &
               emis_ant(ims,  1,jms,p_e_pm25i),emis_ant(ims,  1,jms,p_e_pm25j), &
               emis_ant(ims,  1,jms,p_e_eci),emis_ant(ims,  1,jms,p_e_ecj),     &
               emis_ant(ims,  1,jms,p_e_orgi),emis_ant(ims,  1,jms,p_e_orgj),   &
               emis_ant(ims,  1,jms,p_e_so4j),emis_ant(ims,  1,jms,p_e_so4c),   &
               emis_ant(ims,  1,jms,p_e_no3j),emis_ant(ims,  1,jms,p_e_no3c),   &
               emis_ant(ims,  1,jms,p_e_orgc),emis_ant(ims,  1,jms,p_e_ecc),    &
               emis_ant(ims,  1,jms,p_e_ecc),emis_ant(ims,  1,jms,p_e_ecc),     &
               emis_ant(ims,  1,jms,p_e_ecc),emis_ant(ims,  1,jms,p_e_ecc),     &
               emis_ant(ims,  1,jms,p_e_ecc),emis_ant(ims,  1,jms,p_e_ecc),     &
               emis_ant(ims,  1,jms,p_e_ecc))


	p1st = param_first_scalar




	do 1900 n = 1, nsize_aer(itype)

	do 1830 j = jts, jte
	do 1820 k =   1, min(config_flags%kemit,kte)
	do 1810 i = its, ite



	aem_so4 = fr_aem_mosaic_f(n)*emis_ant(i,k,j,p_e_so4j)  &
	        + fr_aem_mosaic_c(n)*emis_ant(i,k,j,p_e_so4c)  &
	        + fr_aem_sorgam_i(n)*emis_ant(i,k,j,p_e_so4i)  &
	        + fr_aem_sorgam_j(n)*emis_ant(i,k,j,p_e_so4j)  

	aem_no3 = fr_aem_mosaic_f(n)*emis_ant(i,k,j,p_e_no3j)   &
	        + fr_aem_mosaic_c(n)*emis_ant(i,k,j,p_e_no3c)   &
	        + fr_aem_sorgam_i(n)*emis_ant(i,k,j,p_e_no3i)   &
	        + fr_aem_sorgam_j(n)*emis_ant(i,k,j,p_e_no3j)  

	aem_oc  = fr_aem_mosaic_f(n)*emis_ant(i,k,j,p_e_orgj)   &
	        + fr_aem_mosaic_c(n)*emis_ant(i,k,j,p_e_orgc)   &
	        + fr_aem_sorgam_i(n)*emis_ant(i,k,j,p_e_orgi)   &
	        + fr_aem_sorgam_j(n)*emis_ant(i,k,j,p_e_orgj)   &
	        + fr_tno_oc_dom(n)*emis_ant(i,k,j,p_e_oc_dom)   &
	        + fr_tno_oc_tra(n)*emis_ant(i,k,j,p_e_oc_tra)   &
	        + fr_tno_mosaic_c(n)*emis_ant(i,k,j,p_e_oc_25_10) &  

	        + fr_aem_gc2mosaic_f(n)*emis_ant(i,k,j,p_e_oc) 


         chem_select_1 : SELECT CASE( config_flags%chem_opt )
         CASE(SAPRC99_MOSAIC_4BIN_VBS2_KPP, SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, &
              SAPRC99_MOSAIC_8BIN_VBS2_KPP)
             aem_oc = 0.0
         END SELECT  chem_select_1 

	aem_bc  = fr_aem_mosaic_f(n)*emis_ant(i,k,j,p_e_ecj)  &
	        + fr_aem_mosaic_c(n)*emis_ant(i,k,j,p_e_ecc)  &
	        + fr_aem_sorgam_i(n)*emis_ant(i,k,j,p_e_eci)  &
	        + fr_aem_sorgam_j(n)*emis_ant(i,k,j,p_e_ecj)  &
	        + fr_tno_bc1(n)*emis_ant(i,k,j,p_e_bc_1)      &
	        + fr_tno_ec25(n)*emis_ant(i,k,j,p_e_ec_1_25)  &
	        + fr_tno_mosaic_c(n)*emis_ant(i,k,j,p_e_ec_25_10) &

	        + fr_aem_gc2mosaic_f(n)*emis_ant(i,k,j,p_e_bc) 


	aem_oin = fr_aem_mosaic_f(n)*emis_ant(i,k,j,p_e_pm25j)   &
	        + fr_aem_mosaic_c(n)*emis_ant(i,k,j,p_e_pm_10)   &
	        + fr_aem_sorgam_i(n)*emis_ant(i,k,j,p_e_pm25i)   &
	        + fr_aem_sorgam_j(n)*emis_ant(i,k,j,p_e_pm25j)   &

	        + fr_aem_sorgam_c(n)*emis_ant(i,k,j,p_e_pm_10)   &

	        + fr_tno_mosaic_f(n)*emis_ant(i,k,j,p_e_oin_25)  &
                + fr_tno_mosaic_c(n)*emis_ant(i,k,j,p_e_oin_10)  &

	        + fr_aem_gc2mosaic_f(n)*emis_ant(i,k,j,p_e_pm25) &
	        + fr_aem_gc2mosaic_f(n)*emis_ant(i,k,j,p_e_pm10) 

  


	aem_nh4 = 0.0
	aem_na  = 0.0
	aem_cl  = 0.0
	aem_ca  = 0.0
	aem_co3 = 0.0
	aem_msa = 0.0

  chem_select_2 : SELECT CASE( config_flags%chem_opt )
  CASE(MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP)
    aem_nh4 = fr_aem_mosaic_f(n)*emis_ant(i,k,j,p_e_nh4j)   &
            + fr_aem_mosaic_c(n)*emis_ant(i,k,j,p_e_nh4c)   &
            + fr_aem_sorgam_i(n)*emis_ant(i,k,j,p_e_nh4i)   &
            + fr_aem_sorgam_j(n)*emis_ant(i,k,j,p_e_nh4j)

    aem_na = fr_aem_mosaic_f(n)*emis_ant(i,k,j,p_e_naj)   &
            + fr_aem_mosaic_c(n)*emis_ant(i,k,j,p_e_nac)   &
            + fr_aem_sorgam_i(n)*emis_ant(i,k,j,p_e_nai)   &
            + fr_aem_sorgam_j(n)*emis_ant(i,k,j,p_e_naj)

    aem_cl = fr_aem_mosaic_f(n)*emis_ant(i,k,j,p_e_clj)   &
            + fr_aem_mosaic_c(n)*emis_ant(i,k,j,p_e_clc)   &
            + fr_aem_sorgam_i(n)*emis_ant(i,k,j,p_e_cli)   &
            + fr_aem_sorgam_j(n)*emis_ant(i,k,j,p_e_clj)
   END SELECT chem_select_2



	aem_num =   &
	(aem_so4/dens_so4_aer) + (aem_no3/dens_no3_aer) +   &
	(aem_cl /dens_cl_aer ) + (aem_msa/dens_msa_aer) +   &
	(aem_co3/dens_co3_aer) + (aem_nh4/dens_nh4_aer) +   &
	(aem_na /dens_na_aer ) + (aem_ca /dens_ca_aer ) +   &
	(aem_oin/dens_oin_aer) + (aem_oc /dens_oc_aer ) +   &
	(aem_bc /dens_bc_aer )



	aem_num = aem_num*1.0e-6/volumcen_sect(n,itype)



	fact = (dtstep/dz8w(i,k,j))*alt(i,k,j)


	l = lptr_so4_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_so4*fact

	l = lptr_no3_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_no3*fact

	l = lptr_cl_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_cl*fact

	l = lptr_msa_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_msa*fact

	l = lptr_co3_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_co3*fact

	l = lptr_nh4_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_nh4*fact

	l = lptr_na_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_na*fact

	l = lptr_ca_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_ca*fact

	l = lptr_oin_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_oin*fact

	l = lptr_oc_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_oc*fact

	l = lptr_bc_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_bc*fact

	l = numptr_aer(n,itype,iphase)
	if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_num*fact

1810	continue
1820	continue
1830	continue

1900	continue

  if (num_ebu <= 0 ) then
    return
  endif





  aem_so4 = 0.
  aem_no3 = 0.
  aem_oc  = 0.
  aem_bc  = 0.
  aem_oin = 0.
  aem_num = 0.

 fire_select:  SELECT CASE(config_flags%biomass_burn_opt)
   CASE (BIOMASSB,BIOMASSB_MOZC) 
      CALL wrf_debug(15,'mosaic_addemiss: adding fire emissions to MOSAIC')



size_loop: &
        do n = 1, nsize_aer(itype)
          do j = jts, jte
            do k = kts, kte 
              do i = its, ite




                aem_so4 = bburn_mosaic_f(n)*ebu(i,k,j,p_ebu_sulf)*0 
                aem_oc  = bburn_mosaic_f(n)*ebu(i,k,j,p_ebu_oc)    
                aem_bc  = bburn_mosaic_f(n)*ebu(i,k,j,p_ebu_bc) 


                aem_oin = bburn_mosaic_f(n)*ebu(i,k,j,p_ebu_pm25)  &
                        + bburn_mosaic_c(n)*ebu(i,k,j,p_ebu_pm10) 

         chem_select_3 : SELECT CASE( config_flags%chem_opt )
            CASE(SAPRC99_MOSAIC_4BIN_VBS2_KPP) 
              aem_oc = 0.0
            END SELECT  chem_select_3 
 

          aem_nh4 = 0.0
          aem_na  = 0.0
          aem_cl  = 0.0
          aem_ca  = 0.0
          aem_co3 = 0.0
          aem_msa = 0.0

	
	
		aem_num =   &
		(aem_so4/dens_so4_aer) + (aem_no3/dens_no3_aer) +   &
		(aem_cl /dens_cl_aer ) + (aem_msa/dens_msa_aer) +   &
		(aem_co3/dens_co3_aer) + (aem_nh4/dens_nh4_aer) +   &
		(aem_na /dens_na_aer ) + (aem_ca /dens_ca_aer ) +   &
		(aem_oin/dens_oin_aer) + (aem_oc /dens_oc_aer ) +   &
		(aem_bc /dens_bc_aer )

	
	
		aem_num = aem_num*1.0e-6/volumcen_sect(n,itype)

	
	
		fact = (dtstep/dz8w(i,k,j))*alt(i,k,j)

	
		l = lptr_so4_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_so4*fact

		l = lptr_no3_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_no3*fact

		l = lptr_cl_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_cl*fact

		l = lptr_msa_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_msa*fact

		l = lptr_co3_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_co3*fact

		l = lptr_nh4_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_nh4*fact

		l = lptr_na_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_na*fact

		l = lptr_ca_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_ca*fact

		l = lptr_oin_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_oin*fact

		l = lptr_oc_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_oc*fact

		l = lptr_bc_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_bc*fact

		l = numptr_aer(n,itype,iphase)
		if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) + aem_num*fact

	      end do
	    end do
	  end do
	end do size_loop
   END SELECT fire_select




	if (mosaic_addemiss_masscheck > 0) call addemiss_masscheck(        &
               id, config_flags, 2, 'mosaic_ademiss',                      &
               dtstep, efact1, dz8w, chem, chem_sum,                       &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte,                                  &
               14,                                                         &
               emis_ant(ims,  1,jms,p_e_pm_10),emis_ant(ims,  1,jms,p_e_pm_25), &
               emis_ant(ims,  1,jms,p_e_pm25i),emis_ant(ims,  1,jms,p_e_pm25j), &
               emis_ant(ims,  1,jms,p_e_eci),emis_ant(ims,  1,jms,p_e_ecj),     &
               emis_ant(ims,  1,jms,p_e_orgi),emis_ant(ims,  1,jms,p_e_orgj),   &
               emis_ant(ims,  1,jms,p_e_so4j),emis_ant(ims,  1,jms,p_e_so4c),   &
               emis_ant(ims,  1,jms,p_e_no3j),emis_ant(ims,  1,jms,p_e_no3c),   &
               emis_ant(ims,  1,jms,p_e_orgc),emis_ant(ims,  1,jms,p_e_ecc),    &
               emis_ant(ims,  1,jms,p_e_ecc),emis_ant(ims,  1,jms,p_e_ecc),     &
               emis_ant(ims,  1,jms,p_e_ecc),emis_ant(ims,  1,jms,p_e_ecc),     &
               emis_ant(ims,  1,jms,p_e_ecc),emis_ant(ims,  1,jms,p_e_ecc),     &
               emis_ant(ims,  1,jms,p_e_ecc))



	if (seasalt_emiss_active > 0)   &
	    call mosaic_seasalt_emiss(                                     &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               seas_flux,                                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte, seasalt_emiss_active             )










        if (dust_opt == 2) then 
            
            call wrf_message("WARNING: You are calling DUSTRAN dust emission scheme with MOSAIC, which is highly experimental and not recommended for use. Please use dust_opt==13")
            
            call mosaic_dust_emiss( slai, ust, smois, ivgtyp, isltyp,      &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
        endif
        
        if (dust_opt == 13) then 
        
         call mosaic_dust_gocartemis (ktau,dtstep,config_flags%num_soil_layers,alt,u_phy,  &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,                   &
         ivgtyp,isltyp,xland,dx,g,                                         &
         dust_flux,                                                        &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte             ) 
        endif


	return


   END subroutine mosaic_addemiss




   subroutine mosaic_seasalt_emiss(                                        &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               seas_flux,                                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte, seasalt_emiss_active             )






   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem, param_first_scalar
   USE module_data_mosaic_asect

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep
   INTEGER, INTENT(IN) ::    seasalt_emiss_active


   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: seas_flux


	integer i, j, k, l, l_na, l_cl, n, l_oc
	integer iphase, itype
	integer p1st

	real dum, dumdlo, dumdhi, dumoceanfrac, dumspd10
	real factaa, factbb, fracna, fraccl, fracorg

	real :: ssemfact_numb( maxd_asize, maxd_atype )
	real :: ssemfact_mass( maxd_asize, maxd_atype )


	
	
	
    real, save :: alpha_f1(4) =   &
                 (/ 12.328, 38.077, 102.31, 281.65 /)
    real, save :: alpha_f2(4) =   &
                 (/ 2.2958, 8.0935, 25.251, 46.80 /)
    real, save :: alpha_f3(4) =   &
                 (/ 0.00452, 0.00705, 0.00080, 0.000761 /)
    real, save :: beta_f1(4) =   &
                 (/ 0.0311, -0.031633, 0.013154, -0.0017762 /)
    real, save :: beta_f2(4) =   &
                 (/ -13.916, 35.73, -9.7651, 1.1665 /)
    real, save :: beta_f3(4) =   &
                 (/ 4747.8, 12920.0, 7313.4, 6610.0 /)
	real :: nti(4), dp0gi(4)
	
	
	real, save :: oc02um(2) = (/ 0.0, 280.0 /)		
	
	
	
	real, save :: org_frac_low_activity(8) =   &
       (/ 0.05,  0.05,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 /)
	real, save :: org_frac_high_activity(8) =   &
       (/ 0.2,  0.2,  0.1,  0.01,  0.01,  0.01,  0.01,  0.01 /)
	

    p1st = PARAM_FIRST_SCALAR


	itype = 1
	iphase = ai_phase

	
	if(seasalt_emiss_active .eq. 3 .or. seasalt_emiss_active .eq. 4)then
		do i=1,4
			nti(i)   = beta_f1(i) * oc02um(seasalt_emiss_active-2)**2.0 + beta_f2(i) &
							* oc02um(seasalt_emiss_active-2) + beta_f3(i)
			dp0gi(i) = alpha_f1(i) + alpha_f2(i) &
							* exp(-alpha_f3(i)*oc02um(seasalt_emiss_active-2))
		end do
	end if



	do n = 1, nsize_aer(itype)

	    
	    if(seasalt_emiss_active == 1)then



	    	dumdlo = max( dlo_sect(n,itype), 0.02e-4 )
	    	dumdhi = max( dhi_sect(n,itype), 0.02e-4 )
	    	call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
			ssemfact_numb(n,itype), dum, ssemfact_mass(n,itype) )
		elseif(seasalt_emiss_active .eq. 3 .or. seasalt_emiss_active .eq. 4)then
		    call seasalt_emit_feuntes_1bin( dlo_sect(n,itype), dhi_sect(n,itype),  &
		        ssemfact_numb(n,itype), dum, ssemfact_mass(n,itype), nti, dp0gi, oc02um(seasalt_emiss_active-2) )
		endif


	    ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*1.0e6
	end do

  seas_flux(:,:) = 0.0


	k = kts
	do j = jts, jte
	do i = its, ite
	
		
		
		
		if( xland(i,j) < 1.5 ) cycle
	
		
		
		dumoceanfrac = 1. 
		dumspd10 = dumoceanfrac* &
			 ( (u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5*3.41) )
	
		
		
		
		factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
	
		factbb = factaa * dumspd10
	
		if(seasalt_emiss_active == 1)then
		
		seas_flux(i,j) = dumspd10*SUM(ssemfact_mass(1:nsize_aer(itype),itype))
		
			
			fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
			fraccl = 1.0 - fracna
		
			do n = 1, nsize_aer(itype)
		
				
				l_na = lptr_na_aer(n,itype,iphase)
				l_cl = lptr_cl_aer(n,itype,iphase)
				if ((l_na >= p1st) .and. (l_cl >= p1st)) then
			
					chem(i,k,j,l_na) = chem(i,k,j,l_na) +   &
						factbb * ssemfact_mass(n,itype) * fracna
			
					chem(i,k,j,l_cl) = chem(i,k,j,l_cl) +   &
						factbb * ssemfact_mass(n,itype) * fraccl
			
					l = numptr_aer(n,itype,iphase)
					if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
										factbb * ssemfact_numb(n,itype)
		
				end if
				
			end do 
			
			
		elseif(seasalt_emiss_active.eq.3 .or. seasalt_emiss_active.eq.4)then
			do n = 1, nsize_aer(itype)

				
				
				if(seasalt_emiss_active.eq.3)then
					fracorg = org_frac_low_activity(n)
				elseif(seasalt_emiss_active.eq.4)then
					fracorg = org_frac_high_activity(n)
				endif
				fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
				fraccl = 1.0 - fracna
				fracna = (1.0-fracorg)*fracna
				fraccl = (1.0-fracorg)*fraccl
		
				
				l_na = lptr_na_aer(n,itype,iphase)
				l_cl = lptr_cl_aer(n,itype,iphase)
				l_oc = lptr_oc_aer(n,itype,iphase)
				if ((l_na >= p1st) .and. (l_cl >= p1st) .and. (l_oc >= p1st)) then
			
					chem(i,k,j,l_na) = chem(i,k,j,l_na) +   &
						factbb * ssemfact_mass(n,itype) * fracna
			
					chem(i,k,j,l_cl) = chem(i,k,j,l_cl) +   &
						factbb * ssemfact_mass(n,itype) * fraccl

					chem(i,k,j,l_oc) = chem(i,k,j,l_oc) +   &
						factbb * ssemfact_mass(n,itype) * fracorg
			
					l = numptr_aer(n,itype,iphase)
					if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
										factbb * ssemfact_numb(n,itype)
		
				end if
				
			end do 
		endif

	end do 
	end do 

	return

   END subroutine mosaic_seasalt_emiss







	subroutine seasalt_emitfactors_1bin( ireduce_smallr_emit,	&
      		dpdrylo_cm, dpdryhi_cm,	  &
                emitfact_numb, emitfact_surf, emitfact_mass )























	implicit none


	integer ireduce_smallr_emit
	real dpdrylo_cm, dpdryhi_cm,				&
                emitfact_numb, emitfact_surf, emitfact_mass


	integer isub_bin, nsub_bin

	real alnrdrylo
	real drydens, drydens_43pi_em12, x_4pi_em8
	real dum, dumadjust, dumb, dumexpb
	real dumsum_na, dumsum_ma, dumsum_sa
	real drwet, dlnrdry
	real df0drwet, df0dlnrdry, df0dlnrdry_star
	real relhum
	real rdry, rdrylo, rdryhi, rdryaa, rdrybb
	real rdrylowermost, rdryuppermost, rdry_star
	real rwet, rwetaa, rwetbb
	real rdry_cm, rwet_cm
	real sigmag_star
	real xmdry, xsdry

	real pi
	parameter (pi = 3.1415936536)



	real c1, c2, c3, c4, onethird
	parameter (c1 = 0.7674)
	parameter (c2 = 3.079)
	parameter (c3 = 2.573e-11)
	parameter (c4 = -1.424)
	parameter (onethird = 1.0/3.0)



	drydens = 2.165

	drydens_43pi_em12 = drydens*(4.0/3.0)*pi*1.0e-12

	x_4pi_em8 = 4.0*pi*1.0e-8

	relhum = 0.80



	rdry_star = 0.1
	if (ireduce_smallr_emit .le. 0) rdry_star = -1.0e20


	sigmag_star = 1.9


	dumsum_na = 0.0
	dumsum_sa = 0.0
	dumsum_ma = 0.0



        rdrylowermost = dpdrylo_cm*0.5e4
        rdryuppermost = dpdryhi_cm*0.5e4








	if (rdryuppermost .le. rdry_star) goto 2000



        rdrylo = max( rdrylowermost, rdry_star )
        rdryhi = rdryuppermost

	nsub_bin = 1000

	alnrdrylo = log( rdrylo )
	dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin


	rdrybb = exp( alnrdrylo )
	rdry_cm = rdrybb*1.0e-4
	rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      		( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
	rwetbb = rwet_cm*1.0e4

	do 1900 isub_bin = 1, nsub_bin



	rdryaa = rdrybb
	rwetaa = rwetbb


	dum = alnrdrylo + isub_bin*dlnrdry
	rdrybb = exp( dum )

	rdry_cm = rdrybb*1.0e-4
	rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      		( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
	rwetbb = rwet_cm*1.0e4


	rdry = sqrt(rdryaa * rdrybb)
	rwet = sqrt(rwetaa * rwetbb)
	drwet = rwetbb - rwetaa


	xmdry = drydens_43pi_em12 * (rdry**3.0)


	xsdry = x_4pi_em8 * (rdry**2.0)



	dumb = ( 0.380 - log10(rwet) ) / 0.650
	dumexpb = exp( -dumb*dumb)
	df0drwet = 1.373 * (rwet**(-3.0)) * 			&
      		(1.0 + 0.057*(rwet**1.05)) * 			&
      		(10.0**(1.19*dumexpb))

	dumsum_na = dumsum_na + drwet*df0drwet
	dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
	dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

1900	continue













2000	if (rdrylowermost .ge. rdry_star) goto 3000


	rdryaa = 0.99*rdry_star
	rdry_cm = rdryaa*1.0e-4
	rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      		( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
	rwetaa = rwet_cm*1.0e4

	rdrybb = 1.01*rdry_star
	rdry_cm = rdrybb*1.0e-4
	rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      		( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
	rwetbb = rwet_cm*1.0e4

	rwet = 0.5*(rwetaa + rwetbb)
	dumb = ( 0.380 - log10(rwet) ) / 0.650
	dumexpb = exp( -dumb*dumb)
	df0drwet = 1.373 * (rwet**(-3.0)) * 			&
      		(1.0 + 0.057*(rwet**1.05)) * 			&
      		(10.0**(1.19*dumexpb))

	drwet = rwetbb - rwetaa
	dlnrdry = log( rdrybb/rdryaa )
	df0dlnrdry_star = df0drwet * (drwet/dlnrdry)




        rdrylo = rdrylowermost
        rdryhi = min( rdryuppermost, rdry_star )

	nsub_bin = 1000

	alnrdrylo = log( rdrylo )
	dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

	do 2900 isub_bin = 1, nsub_bin


	dum = alnrdrylo + (isub_bin-0.5)*dlnrdry
	rdry = exp( dum )


	xmdry = drydens_43pi_em12 * (rdry**3.0)


	xsdry = x_4pi_em8 * (rdry**2.0)


	dum = log( rdry/rdry_star ) / log( sigmag_star )
	dumadjust = exp( -0.5*dum*dum )

	df0dlnrdry = df0dlnrdry_star * dumadjust

	dumsum_na = dumsum_na + dlnrdry*df0dlnrdry
	dumsum_ma = dumsum_ma + dlnrdry*df0dlnrdry*xmdry
	dumsum_sa = dumsum_sa + dlnrdry*df0dlnrdry*xsdry

2900	continue





3000	emitfact_numb = dumsum_na
	emitfact_mass = dumsum_ma
	emitfact_surf = dumsum_sa

	return
	end subroutine seasalt_emitfactors_1bin







	subroutine seasalt_emit_feuntes_1bin( &
      		dpdrylo_cm, dpdryhi_cm,	  &
                emitfact_numb, emitfact_surf, emitfact_mass, nti, dp0gi, oc02um  )






















	implicit none

	
	
	real dpdrylo_cm, dpdryhi_cm,				&
                emitfact_numb, emitfact_surf, emitfact_mass
    real, intent(in) :: nti(4), dp0gi(4), oc02um

	
	integer isub_bin, nsub_bin, jd, nsub_lower_bin, nsub_upper_bin

	real alnrdrylo
	real drydens, drydens_43pi_em12, x_4pi_em8, drydens_16pi_em12, x_pi_em8
	real dum, dumadjust, dumb, dumexpb
	real dumsum_na, dumsum_ma, dumsum_sa
	real drwet, dlnrdry, ddwet, ddry, dwet
	real df0drwet, df0dlnrdry, df0dlnrdry_star
	real relhum
	real rdry, rdrylo, rdryhi, rdryaa, rdrybb
	real rdrylowermost, rdryuppermost, rdry_star
	real rwet, rwetaa, rwetbb
	real rdry_cm, rwet_cm
	real sigmag_star
	real xmdry, xsdry
	
	real ddrylo, ddryhi, alogddrylo
	real ddrybb, ddry_cm, dwet_cm, dwetbb
	real ddryaa, dwetaa, dlogddry, df0dlogddry

	real pi
	parameter (pi = 3.1415936536)

	
	
	real c1, c2, c3, c4, onethird
	parameter (c1 = 0.7674)
	parameter (c2 = 3.079)
	parameter (c3 = 2.573e-11)
	parameter (c4 = -1.424)
	parameter (onethird = 1.0/3.0)


	
	
	real, save :: width_ssinorg_f(4) = &
				 (/ 1.55, 1.7, 1.5, 1.7 /)
	real, save :: width_ssorg_f(4) = &
				 (/ 1.55, 1.9, 1.5, 1.7 /)
	real :: width_ss_f(4), frac
	
	
	
	
	
	
	
	real, parameter :: scale_factor = 58.3/(0.0146)*3.84e-6 


	
	if(oc02um.gt.0e0)then
		width_ss_f = width_ssorg_f
	else
		width_ss_f = width_ssinorg_f
	end if


	
	drydens = 2.165
	
	drydens_43pi_em12 = drydens*(4.0/3.0)*pi*1.0e-12
	
	drydens_16pi_em12 = drydens*(1.0/6.0)*pi*1.0e-12
	
	x_4pi_em8 = 4.0*pi*1.0e-8
	
	x_pi_em8 = pi*1.0e-8
	
	relhum = 0.80

	
	
	
	rdry_star = 0.45 / 2.0
	
	
	
	sigmag_star = 1.9

	
	dumsum_na = 0.0
	dumsum_sa = 0.0
	dumsum_ma = 0.0

	
	
    rdrylowermost = dpdrylo_cm*0.5e4
    rdryuppermost = dpdryhi_cm*0.5e4


	
	
	
	
	if (rdrylowermost .ge. rdry_star) then

		
		
        rdrylo = rdrylowermost
        rdryhi = rdryuppermost

		nsub_bin = 1000

		alnrdrylo = log( rdrylo )
		dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

		
		rdrybb = exp( alnrdrylo )
		rdry_cm = rdrybb*1.0e-4
		rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      			( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
		rwetbb = rwet_cm*1.0e4

		do isub_bin = 1, nsub_bin

			
			
			rdryaa = rdrybb
			rwetaa = rwetbb

			
			dum = alnrdrylo + isub_bin*dlnrdry
			rdrybb = exp( dum )

			rdry_cm = rdrybb*1.0e-4
			rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      				( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
			rwetbb = rwet_cm*1.0e4

			
			rdry = sqrt(rdryaa * rdrybb)
			rwet = sqrt(rwetaa * rwetbb)
			drwet = rwetbb - rwetaa

			
			xmdry = drydens_43pi_em12 * (rdry**3.0)

			
			xsdry = x_4pi_em8 * (rdry**2.0)

			
			
			dumb = ( 0.380 - log10(rwet) ) / 0.650
			dumexpb = exp( -dumb*dumb)
			df0drwet = 1.373 * (rwet**(-3.0)) * 			&
      				(1.0 + 0.057*(rwet**1.05)) * 			&
      				(10.0**(1.19*dumexpb))

			dumsum_na = dumsum_na + drwet*df0drwet
			dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
			dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

		end do


	
	
	
	
	
	
	elseif (rdryuppermost .gt. rdry_star) then

		
		frac = (log(rdry_star)-log(rdrylowermost)) / (log(rdryuppermost)-log(rdrylowermost))
		nsub_lower_bin = floor(frac*1000.0)	 
		nsub_upper_bin = 1000-nsub_lower_bin 


	
		
		
        ddrylo = rdrylowermost*2.0
        ddryhi = rdry_star*2.0

		nsub_bin = nsub_lower_bin

		alogddrylo = log10( ddrylo )
		dlogddry = (log10( ddryhi ) - alogddrylo)/nsub_bin

		
		ddrybb = 10.0**( alogddrylo )
		ddry_cm = ddrybb*1.0e-4
		dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/		&
      			( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
		dwetbb = dwet_cm*1.0e4

		do isub_bin = 1, nsub_bin

			
			
			ddryaa = ddrybb
			

			
			dum = alogddrylo + isub_bin*dlogddry
			ddrybb = 10.0**( dum )

			ddry_cm = ddrybb*1.0e-4
			
      		
			

			
			ddry = sqrt(ddryaa * ddrybb)
			

			
			xmdry = drydens_16pi_em12 * (ddry**3.0)

			
			xsdry = x_pi_em8 * (ddry**2.0)

			
			
			df0dlogddry = 0.0
			do jd = 1, 4
				df0dlogddry = df0dlogddry + nti(jd)/( (2.0*pi)**0.5 * log10(width_ss_f(jd)) ) &
							* exp(-1.0/2.0 * (log10(ddry*1e3/dp0gi(jd))/log10(width_ss_f(jd)))**2.0 )
			end do

			dumsum_na = dumsum_na + dlogddry*df0dlogddry*scale_factor
			dumsum_ma = dumsum_ma + dlogddry*df0dlogddry*scale_factor*xmdry
			dumsum_sa = dumsum_sa + dlogddry*df0dlogddry*scale_factor*xsdry

		end do

	
	
		
		
        rdrylo = rdry_star
        rdryhi = rdryuppermost

		nsub_bin = nsub_upper_bin

		alnrdrylo = log( rdrylo )
		dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

		
		rdrybb = exp( alnrdrylo )
		rdry_cm = rdrybb*1.0e-4
		rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      			( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
		rwetbb = rwet_cm*1.0e4

		do isub_bin = 1, nsub_bin

			
			
			rdryaa = rdrybb
			rwetaa = rwetbb

			
			dum = alnrdrylo + isub_bin*dlnrdry
			rdrybb = exp( dum )

			rdry_cm = rdrybb*1.0e-4
			rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      				( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
			rwetbb = rwet_cm*1.0e4

			
			rdry = sqrt(rdryaa * rdrybb)
			rwet = sqrt(rwetaa * rwetbb)
			drwet = rwetbb - rwetaa

			
			xmdry = drydens_43pi_em12 * (rdry**3.0)

			
			xsdry = x_4pi_em8 * (rdry**2.0)

			
			
			dumb = ( 0.380 - log10(rwet) ) / 0.650
			dumexpb = exp( -dumb*dumb)
			df0drwet = 1.373 * (rwet**(-3.0)) * 			&
      				(1.0 + 0.057*(rwet**1.05)) * 			&
      				(10.0**(1.19*dumexpb))

			dumsum_na = dumsum_na + drwet*df0drwet
			dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
			dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

		end do



	
	
	
	
	else
	
		
		
        ddrylo = rdrylowermost*2.0
        ddryhi = rdryuppermost*2.0

		nsub_bin = 1000

		alogddrylo = log10( ddrylo )
		dlogddry = (log10( ddryhi ) - alogddrylo)/nsub_bin

		
		ddrybb = 10.0**( alogddrylo )
		ddry_cm = ddrybb*1.0e-4
		dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/		&
      			( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
		dwetbb = dwet_cm*1.0e4

		do isub_bin = 1, nsub_bin

			
			
			ddryaa = ddrybb
			

			
			dum = alogddrylo + isub_bin*dlogddry
			ddrybb = 10.0**( dum )

			ddry_cm = ddrybb*1.0e-4
			
      		
			

			
			ddry = sqrt(ddryaa * ddrybb)
			

			
			xmdry = drydens_16pi_em12 * (ddry**3.0)

			
			xsdry = x_pi_em8 * (ddry**2.0)

			
			
			df0dlogddry = 0.0
			do jd = 1, 4
				df0dlogddry = df0dlogddry + nti(jd)/( (2.0*pi)**0.5 * log10(width_ss_f(jd)) ) &
							* exp(-1.0/2.0 * (log10(ddry*1e3/dp0gi(jd))/log10(width_ss_f(jd)))**2.0 )
			end do

			dumsum_na = dumsum_na + dlogddry*df0dlogddry*scale_factor
			dumsum_ma = dumsum_ma + dlogddry*df0dlogddry*scale_factor*xmdry
			dumsum_sa = dumsum_sa + dlogddry*df0dlogddry*scale_factor*xsdry

		end do
	


	end if
	
	
	
	emitfact_numb = dumsum_na
	emitfact_mass = dumsum_ma
	emitfact_surf = dumsum_sa


	return
	end subroutine seasalt_emit_feuntes_1bin






END MODULE module_mosaic_addemiss



   subroutine mosaic_dust_emiss(  slai,ust, smois, ivgtyp, isltyp,         &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )


























   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem, param_first_scalar    
   USE module_data_mosaic_asect

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep


   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland, slai, ust
   INTEGER,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   ivgtyp, isltyp


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w

   REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,     &
          INTENT(INOUT) ::   smois


        integer i, j, k, l, l_oin, l_ca, l_co3, n, ii
        integer iphase, itype, izob
        integer p1st

        real dum, dumdlo, dumdhi, dumlandfrac, dumspd10
        real factaa, factbb, fracoin, fracca, fracco3, fractot
        real ustart, ustar1, ustart0
        real alphamask, f8, f50, f51, f52, wetfactor, sumdelta, ftot
        real smois_grav, wp, pclay
        real :: beta(4,7)
        real :: gamma(4), delta(4)
        real :: sz(8)
        real :: dustflux, densdust, mass1part
        real :: c_const

        p1st = param_first_scalar








        beta(1,1)=0.12
        beta(2,1)=0.04
        beta(3,1)=0.04
        beta(4,1)=0.80
        beta(1,2)=0.34
        beta(2,2)=0.28
        beta(3,2)=0.28
        beta(4,2)=0.10
        beta(1,3)=0.45
        beta(2,3)=0.15
        beta(3,3)=0.15
        beta(4,3)=0.25
        beta(1,4)=0.12
        beta(2,4)=0.09
        beta(3,4)=0.09
        beta(4,4)=0.70
        beta(1,5)=0.40
        beta(2,5)=0.05
        beta(3,5)=0.05
        beta(4,5)=0.50
        beta(1,6)=0.34
        beta(2,6)=0.18
        beta(3,6)=0.18
        beta(4,6)=0.30
        beta(1,7)=0.22
        beta(2,7)=0.09
        beta(3,7)=0.09
        beta(4,7)=0.60
        gamma(1)=0.08
        gamma(2)=1.00
        gamma(3)=1.00
        gamma(4)=0.12

















        itype=1
	if (nsize_aer(itype) .eq. 8) then
          
          
          
          
          
          
          
          
          
          
          
          sz(1)=0
          sz(2)=1.78751e-06
          sz(3)=0.000273786
          sz(4)=0.00847978
          sz(5)=0.056055
          sz(6)=0.0951896
          sz(7)=0.17
          sz(8)=0.67
	else if (nsize_aer(itype) .eq. 4) then
          
          
          
          
          
          
          sz(1)=1.78751e-06
          sz(2)=0.00875357
          sz(3)=0.1512446
          sz(4)=0.84
          sz(5)=0.0
          sz(6)=0.0
          sz(7)=0.0
          sz(8)=0.0
        endif


        itype = 1
        iphase = ai_phase


        k = kts
        do 1830 j = jts, jte
        do 1820 i = its, ite

    if( xland(i,j) > 1.5 ) cycle



        dumlandfrac = 1.
        dumspd10=(u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5)
        if(dumspd10 >= 5.0) then
           dumspd10 = dumlandfrac* &
         ( dumspd10*dumspd10*(dumspd10-5.0))
         else
            dumspd10=0.
         endif












         alphamask=0.001
         if (ivgtyp(i,j) .eq. 7) then
           f8=0.005
           f50=0.00
           f51=0.10
           f51=0.066
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i,j) .eq. 8) then
           f8=0.010
           f50=0.00
           f51=0.00
           f52=0.15
           f52=0.10
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i,j) .eq. 10) then
           f8=0.00
           f50=0.00
           f51=0.01
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif


















         izob=0
         if(isltyp(i,j).eq.1) izob=1
         if(isltyp(i,j).eq.2) izob=1
         if(isltyp(i,j).eq.3) izob=4
         if(isltyp(i,j).eq.4) izob=2
         if(isltyp(i,j).eq.5) izob=2
         if(isltyp(i,j).eq.6) izob=2
         if(isltyp(i,j).eq.7) izob=7
         if(isltyp(i,j).eq.8) izob=2
         if(isltyp(i,j).eq.9) izob=6
         if(isltyp(i,j).eq.10) izob=5
         if(isltyp(i,j).eq.11) izob=2
         if(isltyp(i,j).eq.12) izob=3
         if(isltyp(i,j).ge.13) izob=0
         if(izob.eq.0) goto 1840



         do ii=1,4
           delta(ii)=0.0
         enddo
         sumdelta=0.0
         do ii=1,4
           delta(ii)=beta(ii,izob)*gamma(ii)
           if(ii.lt.4) then
             sumdelta=sumdelta+delta(ii)
           endif
         enddo
         do ii=1,4
           delta(ii)=delta(ii)/sumdelta
         enddo









         pclay=beta(1,izob)*100.
         wp=0.0014*pclay*pclay+0.17*pclay
         smois_grav=(smois(i,1,j)/2.6)*100.
         if(smois_grav.gt.wp) then
           wetfactor=sqrt(1.0+1.21*(smois_grav-wp)**0.68)
         else
           wetfactor=1.0
         endif





         c_const=1.e-14  

         ustar1=ust(i,j)*100.0
         if(ustar1.gt.100.0) ustar1=100.0
         ustart0=20.0
         ustart=ustart0*wetfactor
         if(ustar1.le.ustart) then
           dustflux=0.0
         else
           dustflux=c_const*(ustar1**4)*(1.0-(ustart/ustar1))
         endif
         dustflux=dustflux*10.0

         ftot=0.0
         do ii=1,2
           ftot=ftot+dustflux*alphamask*delta(ii)
         enddo

         ftot=ftot*1.0e+09


         factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
         factbb = factaa * ftot
         fracoin = 0.97
         fracca = 0.03*0.4
         fracco3 = 0.03*0.6
         fractot = fracoin + fracca + fracco3

         do 1810 n = 1, nsize_aer(itype)
            l_oin = lptr_oin_aer(n,itype,iphase)
            l_ca = lptr_ca_aer(n,itype,iphase)
            l_co3 = lptr_co3_aer(n,itype,iphase)
            if (l_oin >= p1st) then
               chem(i,k,j,l_oin) = chem(i,k,j,l_oin) +      &
               factbb * sz(n) * fracoin
               if (l_ca >= p1st)                            &
                  chem(i,k,j,l_ca) = chem(i,k,j,l_ca) +     &
                  factbb * sz(n) * fracca
               if (l_co3 >= p1st)                           &
                  chem(i,k,j,l_co3) = chem(i,k,j,l_co3) +   &
                  factbb * sz(n) * fracco3

               densdust=2.5
               mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
               l = numptr_aer(n,itype,iphase)
               if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
                        factbb * sz(n) * fractot / mass1part
            end if
1810    continue

1840    continue

1820    continue
1830    continue

        return

   END subroutine mosaic_dust_emiss




  subroutine mosaic_dust_gocartemis (ktau,dt,num_soil_layers,alt,u_phy,    &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,                   &
         ivgtyp,isltyp,xland,dx,g,                                         &
         dust_flux,                                                        &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_data_gocart_dust
  USE module_configure
  USE module_state_description
  USE module_model_constants, ONLY: mwdry
  USE module_data_mosaic_asect
  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ktau, num_soil_layers,           &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(INOUT) ::                               smois
   REAL,  DIMENSION( ims:ime , jms:jme, 3 )                   ,               &
          INTENT(IN   ) ::    erod
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     xland
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                     dz8w,p8w,             &
                                              u_phy,v_phy,rho_phy

  REAL, INTENT(IN   ) :: dt,dx,g
  REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: dust_flux



  integer :: nmx,i,j,k,ndt,imx,jmx,lmx
  integer ilwi, start_month
  real*8, DIMENSION (3) :: erodin
  real*8, DIMENSION (5) :: bems
  real*8  w10m,gwet,airden,airmas
  real*8  cdustemis,jdustemis,cdustcon,jdustcon
  real*8  cdustdens,jdustdens,mass1part,jdustdiam,cdustdiam,dp_meanvol_tmp
  real factaa, factbb, fracoin, fracca, fracco3, fractot
  real*8  dxy
  real*8  conver,converi
  real dttt
  real soilfacj,rhosoilj,rhosoilc
  integer p1st,l_oin,l,n
  real :: densdust
  real sz(8),ftot,frac_10um
  real totalemis,accfrac,corfrac,rscale
  integer iphase, itype

  rscale=1.01  
  accfrac=0.07              
  corfrac=0.93              
  accfrac=accfrac*rscale
  corfrac=corfrac*rscale

  p1st = param_first_scalar


















        sz(1)=1.0e-8
        sz(2)=1.0e-6
        sz(3)=3.0e-4
        sz(4)=3.5e-3
        sz(5)=0.018
        sz(6)=0.070
        sz(7)=0.2595
        sz(8)=0.42




        itype = 1
        iphase = ai_phase

        
        IF (nsize_aer(itype) == 4) THEN
          sz(1) = sz(1) + sz(2)
          sz(2) = sz(3) + sz(4)
          sz(3) = sz(5) + sz(6)
          sz(4) = sz(7) + sz(8)
        ENDIF

  conver=1.e-9
  converi=1.e9

  dust_flux(:,:) = 0.0



  nmx=5
  k=kts
  do j=jts,jte
  do i=its,ite


     if(xland(i,j).lt.1.5)then

     ilwi=1
     start_month = 3   
     w10m=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
     airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g   


     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))

    
    
     erodin(1)=erod(i,j,1)
     erodin(2)=erod(i,j,2)
     erodin(3)=erod(i,j,3)


     gwet=smois(i,1,j)/porosity(isltyp(i,j))
     ndt=ifix(dt)
     airden=rho_phy(i,kts,j)
     dxy=dx*dx

    call mosaic_source_du( nmx, dt, &
                     erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,start_month,g)


    
    
    totalemis=(sum(bems(1:5))/dt)*converi/dxy

    dust_flux(i,j) = totalemis

    
                              
    jdustemis = totalemis*accfrac   
    cdustemis = totalemis*corfrac   

    ftot=jdustemis+cdustemis  



         factaa = (dt/dz8w(i,k,j))*alt(i,k,j)
         factbb = factaa * ftot
         fracoin = 0.97
         fracca = 0.03*0.4
         fracco3 = 0.03*0.6
         fractot = fracoin

         do 1810 n = 1, nsize_aer(itype)
            l_oin = lptr_oin_aer(n,itype,iphase)
            if (l_oin >= p1st) then
               chem(i,k,j,l_oin) = chem(i,k,j,l_oin) +      &
               factbb * sz(n) * fracoin

               if (n <= 5) densdust=2.5
               if (n > 5 ) densdust=2.65
               
               if (nsize_aer(itype) == 4) then
                 if (n <= 2) densdust=2.5
                 if (n > 2 ) densdust=2.65
               endif
               mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
               l = numptr_aer(n,itype,iphase)
               if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
                        factbb * sz(n) * fractot / mass1part
            end if
1810    continue

     endif 
  enddo 
  enddo 


end subroutine mosaic_dust_gocartemis

  SUBROUTINE mosaic_source_du( nmx, dt1, &
                     erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,month,g0)























 USE module_data_gocart_dust

  INTEGER, INTENT(IN)    :: nmx
  REAL*8,    INTENT(IN)  :: erod(ndcls)
  INTEGER, INTENT(IN)    :: ilwi,month

  REAL*8,    INTENT(IN)    :: w10m, gwet
  REAL*8,    INTENT(IN)    :: dxy
  REAL*8,    INTENT(IN)    :: airden, airmas
  REAL*8,    INTENT(OUT)   :: bems(nmx)

  REAL*8    :: den(nmx), diam(nmx)
  REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
  REAL, intent(in)    :: g0
  REAL    :: rhoa, g,dt1
  INTEGER :: i, j, n, m, k

  
  
   ch_dust(:,:)=1.0D-9  
  

  
  DO n = 1, nmx
     
     den(n) = den_dust(n)*1.0D-3
     diam(n) = 2.0*reff_dust(n)*1.0D2
     g = g0*1.0E2
     
     m = ipoint(n)
     tsrc = 0.0
              rhoa = airden*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)

              
             IF (gwet < 0.5) THEN  

                 u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet))))
              ELSE
                 
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(m)*dxy  
              IF (ilwi == 1 ) THEN
                 dsrc = ch_dust(n,month)*srce*w10m**2 &
                      * (w10m - u_ts)*dt1  
              ELSE
                 dsrc = 0.0
              END IF
              IF (dsrc < 0.0) dsrc = 0.0

              
              
              bems(n) = dsrc     

  ENDDO

END SUBROUTINE mosaic_source_du




   subroutine addemiss_masscheck( id, config_flags, iflagaa, fromwhere,    &
               dtstep, efact1, dz8w, chem, chem_sum,                       &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte,                                  &
               nemit,                                                      &
               e01, e02, e03, e04, e05, e06, e07, e08, e09, e10,           &
               e11, e12, e13, e14, e15, e16, e17, e18, e19, e20, e21       )










   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id, iflagaa,                             &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte,               &
                                  nemit

   REAL, INTENT(IN   ) ::    dtstep, efact1


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(IN    ) ::   chem


   DOUBLE PRECISION, DIMENSION( num_chem ),                                &
         INTENT(INOUT ) ::   chem_sum


   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   dz8w



  REAL,  DIMENSION( ims:ime , kms:config_flags%kemit , jms:jme ),          &
          INTENT(IN   ) ::                                                 &
               e01, e02, e03, e04, e05, e06, e07, e08, e09, e10,           &
               e11, e12, e13, e14, e15, e16, e17, e18, e19, e20, e21

   character(len=*), intent(in) :: fromwhere


   integer, parameter :: nemit_maxd = 21
   integer :: i, j, k, l
   double precision :: chem_sum_prev
   real :: fact
   real :: emit_sum(nemit_maxd)





	do 1900 l = 1, num_chem

	chem_sum_prev = chem_sum(l)
	chem_sum(l) = 0.0

	do j = jts, jte
	do k = kts, kte
	do i = its, ite
	    chem_sum(l) = chem_sum(l) + dble( chem(i,k,j,l)*dz8w(i,k,j) )
	end do
	end do
	end do

	if (iflagaa == 2) chem_sum(l) =  (chem_sum(l) - chem_sum_prev)

1900	continue

	if (iflagaa /= 2) return



	emit_sum(:) = 0.0

	do 2900 l = 1, min(nemit,nemit_maxd)
	do j = jts, jte
	do k = kts, min(config_flags%kemit,kte)
	do i = its, ite
	    if (l== 1) emit_sum(l) = emit_sum(l) + e01(i,k,j)
	    if (l== 2) emit_sum(l) = emit_sum(l) + e02(i,k,j)
	    if (l== 3) emit_sum(l) = emit_sum(l) + e03(i,k,j)
	    if (l== 4) emit_sum(l) = emit_sum(l) + e04(i,k,j)
	    if (l== 5) emit_sum(l) = emit_sum(l) + e05(i,k,j)
	    if (l== 6) emit_sum(l) = emit_sum(l) + e06(i,k,j)
	    if (l== 7) emit_sum(l) = emit_sum(l) + e07(i,k,j)
	    if (l== 8) emit_sum(l) = emit_sum(l) + e08(i,k,j)
	    if (l== 9) emit_sum(l) = emit_sum(l) + e09(i,k,j)
	    if (l==10) emit_sum(l) = emit_sum(l) + e10(i,k,j)

	    if (l==11) emit_sum(l) = emit_sum(l) + e11(i,k,j)
	    if (l==12) emit_sum(l) = emit_sum(l) + e12(i,k,j)
	    if (l==13) emit_sum(l) = emit_sum(l) + e13(i,k,j)
	    if (l==14) emit_sum(l) = emit_sum(l) + e14(i,k,j)
	    if (l==15) emit_sum(l) = emit_sum(l) + e15(i,k,j)
	    if (l==16) emit_sum(l) = emit_sum(l) + e16(i,k,j)
	    if (l==17) emit_sum(l) = emit_sum(l) + e17(i,k,j)
	    if (l==18) emit_sum(l) = emit_sum(l) + e18(i,k,j)
	    if (l==19) emit_sum(l) = emit_sum(l) + e19(i,k,j)
	    if (l==20) emit_sum(l) = emit_sum(l) + e20(i,k,j)

	    if (l==21) emit_sum(l) = emit_sum(l) + e21(i,k,j)
	end do
	end do
	end do
2900	continue


	print 9110, fromwhere, its, ite, jts, jte
	print 9100, 'chem_sum'
	fact = 1.0/(dtstep*efact1)
	print 9120, (l, fact*chem_sum(l), l=1,num_chem)
	print 9100, 'emit_sum'
	print 9120, (l, emit_sum(l), l=1,min(nemit,nemit_maxd))

9100	format( a )
9110	format( / 'addemiss_masscheck output, fromwhere = ', a /   &
	'its, ite, jts, jte =', 4i5  )
9120	format( 5( i5, 1pe11.3 ) )


	return
   END subroutine addemiss_masscheck

