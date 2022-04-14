MODULE module_cam_mam_addemiss
















   private
   public :: cam_mam_addemiss


CONTAINS




   subroutine cam_mam_addemiss( id, dtstep, u10, v10, alt, dz8w, xland,    &
               config_flags, chem, slai, ust, smois, ivgtyp, isltyp,       &
               emis_ant,ebio_iso,ebio_olt,ebio_oli,rho_phy,                &
               dust_emiss_active, seasalt_emiss_active,                    &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )

















   USE module_configure, only:  grid_config_rec_type
   USE module_state_description



   USE module_data_cam_mam_asect, only:  ai_phase,   &
         dens_so4_aer, dens_nh4_aer, dens_pom_aer,   &
         dens_bc_aer, dens_dust_aer, dens_seas_aer,   &
         lptr_so4_aer, lptr_nh4_aer, lptr_pom_aer,   &
         lptr_bc_aer, lptr_dust_aer, lptr_seas_aer,   &
         maxd_asize, maxd_atype, nsize_aer, ntype_aer, numptr_aer

   USE modal_aero_data, only:  &
         modeptr_accum, modeptr_aitken, modeptr_coarse,  &
         modeptr_coardust, modeptr_finedust, modeptr_pcarbon, &
         modeptr_fineseas,  modeptr_coarseas

   USE module_cam_mam_init, only:  &
         pom_emit_1p4_factor, so4_emit_1p2_factor

   USE module_data_sorgam, only:  &
         dgvem_i, dgvem_j, dgvem_c, sgem_i, sgem_j, sgem_c

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   INTEGER, INTENT(IN) ::    dust_emiss_active, seasalt_emiss_active

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

   REAL,  DIMENSION( ims:ime , jms:jme )      ,                            &
          INTENT(IN   ) ::   ebio_iso, ebio_olt, ebio_oli
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::   rho_phy



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w

   REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,      &
          INTENT(IN   ) ::   smois





   integer, parameter :: anth_emiss_size_opt = 1   


   integer, parameter :: dust_emiss_size_opt = 1   


   integer, parameter :: seas_emiss_size_opt = 1   

   integer :: i, iphase, isize, itype
   integer :: itype_accum, itype_aitken, itype_pcarbon, &
              itype_finedust, itype_coardust
   integer :: j, k
   integer :: l, lemit, lemit_type, lnumb
   integer :: m, n
   integer :: p1st

   real, parameter :: pi = 3.1415926536
   real :: dp_anth_emiss_aitken, dp_anth_emiss_accum, dp_anth_emiss_coarse
   real :: dlo_seas_emiss( maxd_asize, maxd_atype ), &
           dhi_seas_emiss( maxd_asize, maxd_atype )
   real :: fact_mass, factbb_mass, factbb_numb
   real :: tmp_dens, tmp_dp, tmp_vol1p
   real :: total_soag, fact_soag


   character(len=100) :: msg


	iphase = ai_phase










	emiss_inpt_select_1: SELECT CASE( config_flags%emiss_inpt_opt )


	    CASE( emiss_inpt_default, emiss_inpt_pnnl_rs, emiss_inpt_pnnl_mam )
		write(*,*) 'cam_mam_addemiss DOING emiss - emiss_inpt_opt =', &
			config_flags%emiss_inpt_opt

	    CASE DEFAULT
		write(*,*) 'cam_mam_addemiss SKIPPING emiss - emiss_inpt_opt =', &
			config_flags%emiss_inpt_opt
		return

	END SELECT emiss_inpt_select_1


	p1st = param_first_scalar

	itype_accum    = modeptr_accum
	itype_aitken   = modeptr_aitken
	itype_pcarbon  = modeptr_accum
	itype_finedust = modeptr_accum
	itype_coardust = modeptr_coarse


	if (anth_emiss_size_opt == 1) then
	    
	    dp_anth_emiss_aitken = 0.0504e-4
	    dp_anth_emiss_accum  = 0.134e-4
	    dp_anth_emiss_coarse = 2.06e-4
	else if (anth_emiss_size_opt == 2) then
	    
	    dp_anth_emiss_aitken = 1.0e2 * dgvem_i / exp( 1.5 * (log(sgem_i)**2) )
	    dp_anth_emiss_accum  = 1.0e2 * dgvem_j / exp( 1.5 * (log(sgem_j)**2) )
	    dp_anth_emiss_coarse = 1.0e2 * dgvem_c / exp( 1.5 * (log(sgem_c)**2) )
	else
	    write(msg,'(2a,i7)') 'subr cam_mam_addemiss', &
		' - illegal anth_emiss_size_opt = ', anth_emiss_size_opt
	    call wrf_error_fatal3("<stdin>",197,&
msg )
	end if



	emiss_inpt_select_2: SELECT CASE( config_flags%emiss_inpt_opt )

	CASE( emiss_inpt_pnnl_rs )

	do 1900 lemit_type = 1, 11

	iphase = 1
	isize = 1











	if (lemit_type == 1) then
	    lemit = p_e_so4i
	    itype = itype_aitken
	    l = lptr_so4_aer(isize,itype,iphase)
	    tmp_dens = dens_so4_aer
	    factbb_mass = so4_emit_1p2_factor
	    tmp_dp = dp_anth_emiss_aitken
	else if (lemit_type == 2) then
	    lemit = p_e_so4j
	    itype = itype_accum
	    l = lptr_so4_aer(isize,itype,iphase)
	    tmp_dens = dens_so4_aer


	    factbb_mass = so4_emit_1p2_factor
	    tmp_dp = dp_anth_emiss_accum

	else if (lemit_type == 3) then
	    lemit = p_e_no3i
	    cycle   
	else if (lemit_type == 4) then
	    lemit = p_e_no3j
	    cycle   

	else if (lemit_type == 5) then
	    lemit = p_e_orgi
	    itype = itype_pcarbon
	    l = lptr_pom_aer(isize,itype,iphase)
	    tmp_dens = dens_pom_aer


	    factbb_mass = pom_emit_1p4_factor
	    tmp_dp = dp_anth_emiss_accum
	else if (lemit_type == 6) then
	    lemit = p_e_orgj
	    itype = itype_pcarbon
	    l = lptr_pom_aer(isize,itype,iphase)
	    tmp_dens = dens_pom_aer
	    factbb_mass = pom_emit_1p4_factor
	    tmp_dp = dp_anth_emiss_accum

	else if (lemit_type == 7) then
	    lemit = p_e_eci
	    itype = itype_pcarbon
	    l = lptr_bc_aer(isize,itype,iphase)
	    tmp_dens = dens_bc_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_accum
	else if (lemit_type == 8) then
	    lemit = p_e_ecj
	    itype = itype_pcarbon
	    l = lptr_bc_aer(isize,itype,iphase)
	    tmp_dens = dens_bc_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_accum

	else if (lemit_type == 9) then
	    lemit = p_e_pm25i
	    itype = itype_finedust
	    l = lptr_dust_aer(isize,itype,iphase)
	    tmp_dens = dens_dust_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_accum
	else if (lemit_type == 10) then
	    lemit = p_e_pm25j
	    itype = itype_finedust
	    l = lptr_dust_aer(isize,itype,iphase)
	    tmp_dens = dens_dust_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_accum
	else if (lemit_type == 11) then
	    lemit = p_e_pm_10
	    itype = itype_coardust
	    l = lptr_dust_aer(isize,itype,iphase)
	    tmp_dens = dens_dust_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_coarse
	else
	    cycle
	end if

	if ((l < p1st) .or. (l > num_chem)) cycle

	lnumb = numptr_aer(isize,itype,iphase)
	if ((lnumb < p1st) .or. (lnumb > num_chem)) lnumb = -999888777

	tmp_vol1p = (pi/6.0)*(tmp_dp**3)   

	factbb_numb = 1.0e-6/(tmp_dens*tmp_vol1p)

	do j = jts, jte
	do k =   1, min(config_flags%kemit,kte)
	do i = its, ite



	fact_mass = (dtstep/dz8w(i,k,j))*alt(i,k,j)*factbb_mass

	chem(i,k,j,l) = chem(i,k,j,l) + emis_ant(i,k,j,lemit)*fact_mass

	if (lnumb > 0) then
	chem(i,k,j,lnumb) = chem(i,k,j,lnumb) + &
		emis_ant(i,k,j,lemit)*fact_mass*factbb_numb
	end if

	end do 
	end do 
	end do 

1900	continue








        do j = jts, jte
        do k =   1, 1
        do i = its, ite
          fact_soag = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
          total_soag = (emis_ant(i,k,j,p_e_tol) * (92.14*0.15/12.0))              + &
                       (emis_ant(i,k,j,p_e_xyl) * (106.16*0.15/12.0))             + &
                       (emis_ant(i,k,j,p_e_hc5) * (77.6*0.05/12.0))               + &
                       (emis_ant(i,k,j,p_e_olt) * (72.34*0.05/12.0))              + &
                       (emis_ant(i,k,j,p_e_oli) * (75.78*0.05/12.0))              + &
                       (ebio_iso(i,j) * (68.11*0.04/12.0))







          chem(i,k,j,p_soag) = chem(i,k,j,p_soag) + total_soag*fact_soag
        end do 
        end do 
        end do 
        do j = jts, jte
        do k =   2, min(config_flags%kemit,kte)
        do i = its, ite
          fact_soag = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
          total_soag = (emis_ant(i,k,j,p_e_tol) * (92.14*0.15/12.0))              + &
                       (emis_ant(i,k,j,p_e_xyl) * (106.16*0.15/12.0))             + &
                       (emis_ant(i,k,j,p_e_hc5) * (77.6*0.05/12.0))               + &
                       (emis_ant(i,k,j,p_e_olt) * (72.34*0.05/12.0))              + &
                       (emis_ant(i,k,j,p_e_oli) * (75.78*0.05/12.0))





          chem(i,k,j,p_soag) = chem(i,k,j,p_soag) + total_soag*fact_soag
        end do 
        end do 
        end do 

	CASE( emiss_inpt_pnnl_mam )

	do 1910 lemit_type = 1, 14

	iphase = 1
        isize = 1
	if (lemit_type == 1) then
	    lemit = p_e_so4i
	    itype = itype_aitken
	    l = lptr_so4_aer(isize,itype,iphase)
	    factbb_mass = so4_emit_1p2_factor
	else if (lemit_type == 2) then
	    lemit = p_e_so4j
	    itype = itype_accum
	    l = lptr_so4_aer(isize,itype,iphase)


	    factbb_mass = so4_emit_1p2_factor
	else if (lemit_type == 3) then
	    lemit = p_e_orgj
	    itype = itype_pcarbon
	    l = lptr_pom_aer(isize,itype,iphase)
	    factbb_mass = pom_emit_1p4_factor
	else if (lemit_type == 4) then
	    lemit = p_e_ecj
	    itype = itype_pcarbon
	    l = lptr_bc_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 5) then
	    lemit = p_e_dust_a1
	    itype = itype_finedust
	    l = lptr_dust_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 6) then
	    lemit = p_e_dust_a3
	    itype = itype_coardust
	    l = lptr_dust_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 7) then
	    lemit = p_e_ncl_a1
	    itype = itype_accum
	    l = lptr_seas_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 8) then
	    lemit = p_e_ncl_a2
	    itype = itype_aitken
	    l = lptr_seas_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 9) then
	    lemit = p_e_ncl_a3
	    itype = itype_coardust
	    l = lptr_seas_aer(isize,itype,iphase)
	    factbb_mass = 1.0






        else if (lemit_type == 10) then
            lemit = p_e_so4i_num
            itype = itype_aitken
            l = numptr_aer(isize,itype,iphase)
        else if (lemit_type == 11) then
            lemit = p_e_so4j_num
            itype = itype_accum
            l = numptr_aer(isize,itype,iphase)
        else if (lemit_type == 12) then
            lemit = p_e_orgj_num
            itype = itype_accum
            l = numptr_aer(isize,itype,iphase)
        else if (lemit_type == 13) then
            lemit = p_e_ecj_num
            itype = itype_accum
            l = numptr_aer(isize,itype,iphase)
        else if (lemit_type == 14) then
            lemit = p_e_num_a3
            itype = itype_coardust
            l = numptr_aer(isize,itype,iphase)
        else
	    cycle
	end if

	if ((l < p1st) .or. (l > num_chem)) cycle

	do j = jts, jte
	do k =   1, min(config_flags%kemit,kte)
	do i = its, ite



	fact_mass = (dtstep/dz8w(i,k,j))*alt(i,k,j)*factbb_mass

        factbb_numb = (dtstep/dz8w(i,k,j))*alt(i,k,j)
        if (lemit_type < 10) then 
	chem(i,k,j,l) = chem(i,k,j,l) + emis_ant(i,k,j,lemit)*fact_mass
        else
        chem(i,k,j,l) = chem(i,k,j,l) + emis_ant(i,k,j,lemit)*factbb_numb
        end if

	end do 
	end do 
	end do 

1910	continue



	do j = jts, jte
	do k =   1, min(config_flags%kemit,kte)
	do i = its, ite
	  fact_soag = (dtstep/dz8w(i,k,j))*alt(i,k,j)/1000.0
          total_soag = emis_ant(i,k,j,p_e_soag_bigalk)   + &
                       emis_ant(i,k,j,p_e_soag_bigene)   + &
                       emis_ant(i,k,j,p_e_soag_isoprene) + &
                       emis_ant(i,k,j,p_e_soag_terpene)  + &
                       emis_ant(i,k,j,p_e_soag_toluene)
	  chem(i,k,j,p_soag) = chem(i,k,j,p_soag) + total_soag*fact_soag
	end do 
	end do 
	end do 


	END SELECT emiss_inpt_select_2







	dlo_seas_emiss(:,:) = 0.0
	dhi_seas_emiss(:,:) = 0.0

	if (seas_emiss_size_opt == 1) then

	    dlo_seas_emiss(1,modeptr_aitken  ) = 0.02e-4
	    dhi_seas_emiss(1,modeptr_aitken  ) = 0.08e-4
	    dlo_seas_emiss(1,modeptr_accum   ) = 0.08e-4
	    dhi_seas_emiss(1,modeptr_accum   ) = 1.00e-4
	    dlo_seas_emiss(1,modeptr_coarse  ) = 1.00e-4
	    dhi_seas_emiss(1,modeptr_coarse  ) = 10.0e-4

	else if (seas_emiss_size_opt == 2) then

	    dlo_seas_emiss(1,modeptr_accum   ) = 0.10e-4
	    dhi_seas_emiss(1,modeptr_accum   ) = 1.25e-4
	    dlo_seas_emiss(1,modeptr_coarse  ) = 1.25e-4
	    dhi_seas_emiss(1,modeptr_coarse  ) = 10.0e-4

	else
	    write(msg,'(2a,i7)') 'subr cam_mam_addemiss', &
		' - illegal seas_emiss_size_opt = ', seas_emiss_size_opt
	    call wrf_error_fatal3("<stdin>",533,&
msg )
	end if


	if (config_flags%seas_opt == 2) then
	    call cam_mam_mosaic_seasalt_emiss(                             &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dlo_seas_emiss, dhi_seas_emiss,                             &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
	end if







	if (config_flags%dust_opt == 2) then
            call cam_mam_mosaic_dust_emiss(                                &
               slai, ust, smois, ivgtyp, isltyp,                           &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dust_emiss_size_opt,                                        &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
	end if


	return


   END subroutine cam_mam_addemiss




   subroutine cam_mam_mosaic_seasalt_emiss(                                &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dlo_seas_emiss, dhi_seas_emiss,                             &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )






   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem, param_first_scalar
   USE module_data_cam_mam_asect, only:  ai_phase,   &
         lptr_seas_aer,   &
         maxd_asize, maxd_atype, nsize_aer, ntype_aer, numptr_aer


   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep


   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w


   REAL,  DIMENSION( maxd_asize, maxd_atype ), &
          INTENT(IN   ) ::  dlo_seas_emiss, dhi_seas_emiss


	integer i, j, k, l, l_seas, n
	integer iphase, itype
	integer p1st

	real dum, dumdlo, dumdhi, dumoceanfrac, dumspd10
	real factaa, factbb

	real :: ssemfact_numb( maxd_asize, maxd_atype )
	real :: ssemfact_mass( maxd_asize, maxd_atype )

    write(*,*) 'in subr cam_mam_mosaic_seasalt_emiss'
    p1st = PARAM_FIRST_SCALAR


	iphase = ai_phase



	ssemfact_mass(:,:) = 0.0
	ssemfact_numb(:,:) = 0.0
	do itype = 1, ntype_aer
	do n = 1, nsize_aer(itype)
	    dumdlo = max( dlo_seas_emiss(n,itype), 0.1e-4 )
	    dumdhi = max( dhi_seas_emiss(n,itype), 0.1e-4 )
	    if (dumdlo >= dumdhi) cycle
	    call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
		ssemfact_numb(n,itype), dum, ssemfact_mass(n,itype) )


	    ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*1.0e6
	end do
	end do



	k = kts
	do 1830 j = jts, jte
	do 1820 i = its, ite

    
    
    
	if( xland(i,j) < 1.5 ) cycle

    
    
	dumoceanfrac = 1. 
	dumspd10 = dumoceanfrac* &
         ( (u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5*3.41) )




	factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)

	factbb = factaa * dumspd10

	do 1815 itype = 1, ntype_aer
	do 1810 n = 1, nsize_aer(itype)
	    if (ssemfact_mass(n,itype) <= 0.0) cycle


	    l_seas = lptr_seas_aer(n,itype,iphase)
	    if (l_seas < p1st) cycle

	    chem(i,k,j,l_seas) = chem(i,k,j,l_seas) +   &
			factbb * ssemfact_mass(n,itype)

	    l = numptr_aer(n,itype,iphase)
	    if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
			factbb * ssemfact_numb(n,itype)

1810	continue
1815	continue

1820	continue
1830	continue

	return

   END subroutine cam_mam_mosaic_seasalt_emiss







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




   subroutine cam_mam_mosaic_dust_emiss(  slai,ust, smois, ivgtyp, isltyp, &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dust_emiss_size_opt,                                        &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )


























   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem, param_first_scalar

   USE module_data_cam_mam_asect, only:  ai_phase,   &
         lptr_dust_aer, ntype_aer, numptr_aer

   USE modal_aero_data, only:  &
         modeptr_accum, modeptr_coarse, modeptr_coardust, modeptr_finedust

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  dust_emiss_size_opt,                     &
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
          INTENT(IN   ) ::   smois


        integer, parameter :: max_dust_mode = 2
        integer :: i, ii, j, k, l, l_dust, l_numb, n
        integer :: imode, iphase, isize, itype, izob
        integer :: p1st
        integer :: isize_dust_mode(max_dust_mode)
        integer :: itype_dust_mode(max_dust_mode)
        integer :: n_dust_mode

        real :: dum, dumdlo, dumdhi, dumlandfrac, dumspd10
        real :: factaa, factbb, factcc, fracoin, fracca, fracco3, fractot
        real :: ustart, ustar1, ustart0
        real :: alphamask, f8, f50, f51, f52, wetfactor, sumdelta, ftot
        real :: smois_grav, wp, pclay
        real :: beta(4,7)
        real :: gamma(4), delta(4)
        real :: sz(8)
        real :: dustflux, densdust
        real :: mass1part_mos8bin(8)
        real :: sz_wght_dust_mode(8,max_dust_mode)
        real :: dcen_mos8bin(8)

        character(len=100) :: msg

        write(*,*) 'in subr cam_mam_mosaic_dust_emiss'
        p1st = param_first_scalar
















        n_dust_mode = 2
        itype_dust_mode(1) = modeptr_accum
        itype_dust_mode(2) = modeptr_coarse
        isize_dust_mode(:) = 1
        if (dust_emiss_size_opt == 1) then


           sz_wght_dust_mode(1:8,1) = (/ 1.,  1.,  1.,  1.,  0.6, 0.,  0.,  0.  /)
           sz_wght_dust_mode(1:8,2) = (/ 0.,  0.,  0.,  0.,  0.4, 1.,  1.,  1.  /)
        else if (dust_emiss_size_opt == 2) then


           sz_wght_dust_mode(1:8,1) = (/ 1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.  /)
           sz_wght_dust_mode(1:8,2) = (/ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.  /)
        else
            write(msg,'(2a,i7)') 'subr cam_mam_mosaic_dust_emiss', &
                ' - illegal dust_emiss_size_opt = ', dust_emiss_size_opt
            call wrf_error_fatal3("<stdin>",1082,&
msg )
        end if


        dcen_mos8bin(8) = 1.0e-4*(10.0/sqrt(2.0))
        do n = 7, 1, -1
           dcen_mos8bin(n) = dcen_mos8bin(n+1)*0.5
        end do


        densdust=2.5
        do n = 1, 8
           mass1part_mos8bin(n)=0.523598*(dcen_mos8bin(n)**3)*densdust*1.0e06
        end do










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





















        sz(1)=0.0
        sz(2)=0.0
        sz(3)=0.0005
        sz(4)=0.0095
        sz(5)=0.01
        sz(6)=0.06
        sz(7)=0.20
        sz(8)=0.72

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
         if(izob.eq.0) goto 1820



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






         ustar1=ust(i,j)*100.0
         if(ustar1.gt.100.0) ustar1=100.0
         ustart0=20.0
         ustart=ustart0*wetfactor
         if(ustar1.le.ustart) then
           dustflux=0.0
         else
           dustflux=1.0e-14*(ustar1**4)*(1.0-(ustart/ustar1))
         endif
         dustflux=dustflux*10.0

         ftot=0.0
         do ii=1,2
           ftot=ftot+dustflux*alphamask*delta(ii)
         enddo

         ftot=ftot*1.0e+09


         factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
         factbb = factaa * ftot
         fractot = 1.0

         do imode = 1, n_dust_mode

         itype = itype_dust_mode(imode)
         isize = isize_dust_mode(imode)

         l_dust = lptr_dust_aer(isize,itype,iphase)
         if ((l_dust < p1st) .or. (l_dust > num_chem)) cycle
         l_numb = numptr_aer(isize,itype,iphase)
         if ((l_numb < p1st) .or. (l_numb > num_chem)) l_numb = -1



         do n = 1, 8

            if (sz_wght_dust_mode(n,imode) <= 0.0) cycle
            factcc = factbb * sz(n) * fractot * sz_wght_dust_mode(n,imode)

            chem(i,k,j,l_dust) = chem(i,k,j,l_dust) + factcc

            if (l_numb >= p1st) &
            chem(i,k,j,l_numb) = chem(i,k,j,l_numb) + factcc/mass1part_mos8bin(n)

        end do 

        end do 


1820    continue
1830    continue

        return

   END subroutine cam_mam_mosaic_dust_emiss



END MODULE module_cam_mam_addemiss
