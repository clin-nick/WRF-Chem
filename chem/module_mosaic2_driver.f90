





























































































	module module_mosaic2_driver








	implicit none

	contains

















	subroutine mosaic2_aerchem_driver(                        &
		id, curr_secs, ktau, dtstep, ktauc, dtstepc,      &
		config_flags,                                     &
		t_phy, rho_phy, p_phy,                            &
		moist, chem,                                      &
		ids,ide, jds,jde, kds,kde,                        &
		ims,ime, jms,jme, kms,kme,                        &
		its,ite, jts,jte, kts,kte                         )


	use module_configure, only:  grid_config_rec_type, &
	    p_qv, p_so2, p_ho2, p_so4aj, p_corn, p_hcl, p_mtf, &
	    p_so4_a01, p_water_a01, p_num_a01, &
	    p_so4_a04, p_water_a04, p_num_a04
	use module_state_description, only:  num_moist, num_chem, param_first_scalar

	use module_data_mosaic_kind, only:  r8

	use module_data_mosaic_aero, only:  nbin_a_max, msoa_flag1, msoa_vbs_info
	use module_data_mosaic_main, only:  ntot_used

	use module_mosaic_aerdynam_intr, only:  aerosoldynamics




	use module_peg_util, only:  peg_error_fatal, peg_message

	implicit none





























	integer, intent(in) ::              &
		id, ktau, ktauc,                &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte












    real(kind=8), intent(in) :: curr_secs
	real, intent(in) :: dtstep, dtstepc



	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		t_phy, rho_phy, p_phy




	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
		moist


 
	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem



	type(grid_config_rec_type), intent(in) :: config_flags




	integer, parameter :: debug_level=0
	integer :: i, itmpa, istat, it, j, jt, jtmpa, k, kt, l, n
	integer :: levdbg_err, levdbg_info
	integer :: i_force_dump, mode_force_dump
	integer :: idiagaa, idiagbb, ijkcount
	integer :: iter_mesa(nbin_a_max), jaerosolstate(nbin_a_max)
	
        real(r8) :: dtchem
        real(r8) :: tempbox, presbox, airdenbox, relhumbox, swdownbox
        real(r8) :: rbox(ntot_used), rbox0(ntot_used)
        real(r8) :: dp_dry_a(nbin_a_max), dp_wet_a(nbin_a_max)
	real, dimension( ims:ime, jms:jme ) :: swdown

	character*100 msg


        if ( config_flags%aerchem_onoff <= 0 ) return   



	swdown = 0.0  
	              



	mode_force_dump = 0
	levdbg_err = 0
	levdbg_info = 15
	ijkcount = 0
	if ((jde-jds > 2) .or. (ide-ids > 2)) then
	    idiagaa = 0     
	    idiagbb = 0
	else
	    idiagaa = 100   
	    idiagbb = 100
	end if

	if (idiagaa > 0) print 93010, 'entered mosaic2_aerchem_driver - ktau =', ktau



        if (debug_level .ge. 15 .or. idiagaa >= 100) then
        if (ktauc .le. 2) then
        print 93010, ' '
        print 93010, 'rcetestc diagnostics from mosaic2_aerchem_driver'
        print '(a,5i9)', 'msoa_flag1, msoa_vbs_info(1:3)', msoa_flag1, msoa_vbs_info(1:3)
        print 93010, 'id, chem_opt, ktau, ktauc    ',   &
             id, config_flags%chem_opt, ktau, ktauc
        print 93020, 'dtstep, dtstepc                 ',   &
             dtstep, dtstepc
        print 93010, 'ids/e, j, k', ids, ide, jds, jde, kds, kde
        print 93010, 'ims/e, j, k', ims, ime, jms, jme, kms, kme
        print 93010, 'its/e, j, k', its, ite, jts, jte, kts, kte
        print 93010, 'num_chem, param_first_scalar ', num_chem, param_first_scalar
        print 93010, 'p_so2, p_ho2                 ', p_so2, p_ho2
        print 93010, 'p_so4aj, p_corn, p_hcl, p_mtf', p_so4aj, p_corn, p_hcl, p_mtf
        print 93010, 'p_so4_a01, p_water, p_num_a01', p_so4_a01, p_water_a01, p_num_a01
        print 93010, 'p_so4_a04, p_water, p_num_a04', p_so4_a04, p_water_a04, p_num_a04

        k = kts
        print 93020, 't, p, rho, qv at its/kts /jts', t_phy(its,k,jts),   &
	    p_phy(its,k,jts), rho_phy(its,k,jts), moist(its,k,jts,p_qv)
        k = (kts + kte)/2
        print 93020, 't, p, rho, qv at its/ktmi/jts', t_phy(its,k,jts),   &
	    p_phy(its,k,jts), rho_phy(its,k,jts), moist(its,k,jts,p_qv)
        k = kte
        print 93020, 't, p, rho, qv at its/kte /jts', t_phy(its,k,jts),   &
	    p_phy(its,k,jts), rho_phy(its,k,jts), moist(its,k,jts,p_qv)
93010	format( a, 8(1x,i6) )
93020	format( a, 8(1p,e14.6) )
        end if
        end if











main_jt_loop: &
	do jt = jts, jte
main_kt_loop: &
	do kt = kts, kte
main_it_loop: &
	do it = its, ite

	ijkcount = ijkcount + 1
	dtchem = dtstepc




	i_force_dump = 0










	if (idiagaa >= 100) print 93010, &
	    'calling mosaic2_map 0   - ktau,ktauc,id,i,k,j =', ktau, ktauc, id, it, kt, jt
	call  mosaic2_map_tofrom_host( 0,                     &
		ids,ide, jds,jde, kds,kde,                    &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      kt,                         &
		t_phy, p_phy, rho_phy, swdown, moist, chem,   &
		rbox, tempbox, presbox, airdenbox,            &
		relhumbox, swdownbox                          )


	rbox0(1:nbin_a_max) = rbox(1:nbin_a_max)


	if (idiagaa >= 100) print 93010, &
	    'calling aerosoldynamics - ktau,ktauc,id,i,k,j =', ktau, ktauc, id, it, kt, jt

        call aerosoldynamics(                             & 
             idiagbb,                                     &
             id, it, jt, kt,                              &
             ktau, ktauc, dtchem,                         &
             tempbox, presbox, airdenbox, relhumbox,      &
             swdownbox,                                   &
             rbox, iter_mesa, jaerosolstate,              & 
             dp_dry_a, dp_wet_a                           )

        call wrf_debug(300,"mosaic_aerchem_driver: back from aerchemistry")



	if (idiagaa >= 100) print 93010, &
	    'calling mosaic2_map 1   - ktau,ktauc,id,i,k,j =', ktau, ktauc, id, it, kt, jt
	call  mosaic2_map_tofrom_host( 1,                     &
		ids,ide, jds,jde, kds,kde,                    &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      kt,                         &
		t_phy, p_phy, rho_phy, swdown, moist, chem,   &
		rbox, tempbox, presbox, airdenbox,            &
		relhumbox, swdownbox                          )


        end do main_it_loop
        end do main_kt_loop
        end do main_jt_loop



	if (idiagaa > 0) print 93010, 'leaving mosaic2_aerchem_driver - ktau =', ktau

	return
	end subroutine mosaic2_aerchem_driver




	subroutine mosaic2_map_tofrom_host( imap,             &
		ids,ide, jds,jde, kds,kde,                    &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      kt,                         &
		t_phy, p_phy, rho_phy, swdown, moist, chem,   &
		rbox, tempbox, presbox, airdenbox,            &
		relhumbox, swdownbox                          )

	use module_data_mosaic_kind, only:  r8
	use module_data_mosaic_main, only:  ntot_used
	use module_mosaic_csuesat, only:  esat_gchm
	use module_configure, only:  p_qv, p_qc
	use module_state_description, only:  num_moist, num_chem, param_first_scalar
	use module_peg_util, only:  peg_error_fatal, peg_message

	implicit none




	integer, intent(in) :: imap

	integer, intent(in) :: ids, ide, jds, jde, kds, kde
	integer, intent(in) :: ims, ime, jms, jme, kms, kme
	integer, intent(in) :: its, ite, jts, jte, kts, kte

	integer, intent(in) :: it, jt, kt
   
	real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: &
		t_phy, rho_phy, p_phy

        real, intent(in), dimension( ims:ime, jms:jme ) :: swdown

	real, intent(in), &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
		moist
 
	real, intent(inout), &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem

        real(r8), intent(inout), dimension( ntot_used ) :: rbox

        real(r8), intent(inout) :: tempbox, presbox, airdenbox, relhumbox, swdownbox



        integer :: l
	real, parameter :: eps=0.622
	real :: tmp_esat, tmp_pres, tmp_rh, tmp_temp, tmp_vap, tmp_vapsat


        if (imap >= 1) go to 20000

        do l = 1, ntot_used
            if (l >= param_first_scalar .and. l <= num_chem) then
                rbox(l) = max( chem(it,kt,jt,l), 0.0_r8 )
            else
                rbox(l) = 0.0_r8
            end if
        end do


	tempbox   = t_phy(it,kt,jt)
	presbox   = p_phy(it,kt,jt)
	airdenbox = rho_phy(it,kt,jt)
	swdownbox = swdown(it,jt)


	tmp_temp   = t_phy(it,kt,jt)             
	tmp_pres   = p_phy(it,kt,jt)*10.0        
	tmp_esat   = esat_gchm( tmp_temp )       
	tmp_vapsat = tmp_esat / (tmp_pres - (1.0-eps)*tmp_esat)  
	tmp_vap    = moist(it,kt,jt,p_qv) / eps  
	tmp_rh     = tmp_vap / max( tmp_vapsat, 1.0e-20 )

	tmp_rh     = max( 0.0, min( 0.98, tmp_rh ) )  
	relhumbox  = tmp_rh


	return


20000   continue

        do l = 1, ntot_used
            if (l >= param_first_scalar .and. l <= num_chem) then
                chem(it,kt,jt,l) = max( rbox(l), 0.0_r8 )
            end if
        end do

	return
	end subroutine mosaic2_map_tofrom_host




	subroutine init_data_mosaic2_asect( id, config_flags, is_aerosol, vbs_nbin)




	use module_configure, only:   grid_config_rec_type
	use module_state_description, only:  num_chem
	use module_peg_util, only:  peg_error_fatal, peg_message

	use module_data_mosaic_aero, only: &
            mcoag_flag1, mmovesect_flag1, mnewnuc_flag1, &
            nbin_a, nbin_a_max
	use module_data_mosaic_asecthp
	use module_data_mosaic_main, only: &
            naerbin, naerbin_used, naer_max, naer_tot, &
            ntot_max, ntot_used
	use module_data_mosaic_boxmod, only:   name_rbox

	implicit none

	integer, intent(in) :: id
	integer, intent(in) :: vbs_nbin(*)
	type(grid_config_rec_type), intent(in) :: config_flags
	logical, intent(out) :: is_aerosol(num_chem)


	integer :: i, itype, itmpa, j, jtmpa, l, n
	integer :: vbs_uq_aqsoa, vbs_uq_par
	character(len=200) :: msg



        ntot_max  = num_chem
        ntot_used = num_chem


        vbs_uq_aqsoa = 0
        vbs_uq_par = 0


        if ( .not. allocated( name_rbox ) ) allocate( name_rbox(ntot_used) ) 
	do l = 1, ntot_used
	    write( name_rbox(l), '(a,i4.4,15x)' ) 'r', l
	end do





 	call mosaic2_wrfchem_init( 1, vbs_nbin, vbs_uq_aqsoa, vbs_uq_par )



	call peg_message( lunerr, 'call mosaic2_set_mastercomp' )
	call mosaic2_set_mastercomp( )



 	call peg_message( lunerr, 'call mosaic2_set_ntype' )
 	call mosaic2_set_ntype( config_flags%chem_opt, lunerr )



 	call peg_message( lunerr, 'call mosaic2_set_nphase' )
 	call mosaic2_set_nphase( config_flags%chem_opt, lunerr )



 	call peg_message( lunerr, 'call mosaic2_set_nsize' )
 	call mosaic2_set_nsize( config_flags%chem_opt, lunerr )





	n = 0
	do itype = 1, ntype_aer
	    n = n + nsize_aer(itype)
	end do
	nbin_a = n
	if (nbin_a .gt. nbin_a_max) then
	    call peg_error_fatal( lunerr,   &
		'init_data_mosaic2_asect - nbin_a > nbin_a_maxd' )
	end if
        naerbin      = nbin_a
        naerbin_used = nbin_a
        naer_max = naer_tot*naerbin


	call peg_message( lunerr, 'mosaic2_set_3dbin_1dbin_ptrs' )


	call mosaic2_set_3dbin_1dbin_ptrs( config_flags%chem_opt,                           0, lunerr )






	call peg_message( lunerr, 'call mosaic2_set_bin_sizes' )
	call mosaic2_set_bin_sizes( config_flags%chem_opt, lunerr )





	call init_data_mosaic2_ptr( id, config_flags%chem_opt, is_aerosol )





 	call mosaic2_wrfchem_init( 2, vbs_nbin, vbs_uq_aqsoa, vbs_uq_par )





        call set_rbox_gas_ptrs
        call set_rbox_aer_ptrs





        itmpa = max( config_flags%mosaic_aerchem_optaa, 0 )
        if (itmpa >= 10000) then
            jtmpa = mod(itmpa,100)/10  
            if (jtmpa == 0) mmovesect_flag1 = 0   
            if (jtmpa == 1) mmovesect_flag1 = 10  
            if (jtmpa == 2) mmovesect_flag1 = 20  
            if (jtmpa == 5) mmovesect_flag1 = 50  
            if (jtmpa == 6) mmovesect_flag1 = 60  

            jtmpa = mod(itmpa,1000)/100  
            if (jtmpa == 0) mnewnuc_flag1 = 0  
            if (jtmpa == 1) mnewnuc_flag1 = 1  
            if (jtmpa == 2) mnewnuc_flag1 = 2  
            if (jtmpa == 3) mnewnuc_flag1 = 3  
            if (jtmpa == 5) mnewnuc_flag1 = 11 
            if (jtmpa == 6) mnewnuc_flag1 = 12 

            jtmpa = mod(itmpa,10000)/1000  
            if (jtmpa == 0) mcoag_flag1 = 0  
            if (jtmpa == 1) mcoag_flag1 = 1  
            if (jtmpa == 2) mcoag_flag1 = 10 
            if (jtmpa == 6) mcoag_flag1 = 60 

            if (ntype_aer > 1) then
                
                if (mmovesect_flag1 == 10) mmovesect_flag1 = 50
                if (mmovesect_flag1 == 20) mmovesect_flag1 = 60
                if (mcoag_flag1 ==  1) mcoag_flag1 = 60
                if (mcoag_flag1 == 10) mcoag_flag1 = 60
            end if
        end if 
	write(msg,'(a,4i10)') &
           'mosaic_aerchem_optaa, mmovesect_flag1, mnewnuc_flag1, mcoag_flag1 =', &
            config_flags%mosaic_aerchem_optaa, mmovesect_flag1, mnewnuc_flag1, mcoag_flag1
	call peg_message( lunerr, 'call mosaic2_set_bin_sizes' )


























    

	end subroutine init_data_mosaic2_asect



	subroutine init_data_mosaic2_ptr( id, chem_opt, is_aerosol )




	use module_configure
	use module_state_description, only:  num_chem, p1st => param_first_scalar

	use module_data_mosaic_asecthp
	use module_data_mosaic_main, only:  ntot_used
	use module_peg_util, only:  peg_error_fatal, peg_message


	implicit none


        integer, intent(in)  :: id, chem_opt
        logical, intent(out) :: is_aerosol(num_chem)


	integer l, ll, n
	integer icomp, isize, itype, iphase

	character*200 msg





	call peg_message( lunerr, 'call mosaic2_set_all_lnw_ptr' )
	call mosaic2_set_all_lnw_ptr( chem_opt, ntot_used, lunerr )




	call peg_message( lunerr, 'call mosaic2_set_ncomp' )
	call mosaic2_set_ncomp( chem_opt, lunerr )




	call peg_message( lunerr, 'call mosaic2_set_massptr' )
	call mosaic2_set_massptr( chem_opt, ntot_used, lunerr )
	call peg_message( lunerr, 'done mosaic2_set_massptr' )





	call mosaic2_output_wrfchem_pointers( id )













	if ( 1 == 0 ) then
	end if 





	call mosaic2_set_gas_kptrs_names






	call mosaic2_set_otheraa( chem_opt, is_aerosol )


	return
	end subroutine init_data_mosaic2_ptr

  

	subroutine mosaic2_wrfchem_init( ipass, vbs_nbin, vbs_uq_aqsoa, vbs_uq_par )











	use module_data_mosaic_kind, only: r8
	use module_data_mosaic_aero,      only: &
	    alpha_ASTEM, ptol_mol_ASTEM, rtol_eqb_ASTEM, &
	    mGAS_AER_XFER, mDYNAMIC_SOLVER, &
	    mhyst_method, mhyst_uporlo_waterhyst, &
	    msize_framework, msectional, &
	    method_bcfrac, method_kappa, &
	    maersize_init_flag1, mcoag_flag1, ifreq_coag, &
	    mmovesect_flag1, mnewnuc_flag1, msectional_flag1, msectional_flag2, &
	    msoa_flag1, msoa_vbs_info, &
	    use_cam5mam_soa_params, use_cam5mam_accom_coefs

	use module_data_mosaic_main,      only: &
	    ipmcmos, m_partmc_mosaic, maer, mcld, mgas, &
	    maeroptic, mphoto, mshellcore, msolar

	use module_data_mosaic_asecthp,   only: ntype_aer, ntype_md1_aer, ntype_md2_aer

	use module_data_mosaic_constants, only: pi, piover4, piover6, deg2rad, third, avogad


	use module_mosaic_init_aerpar,    only: mosaic_init_aer_params
	   
	integer, intent(in) :: ipass
	integer, intent(in) :: vbs_nbin(*), vbs_uq_aqsoa, vbs_uq_par



        if (ipass == 1) then
           msoa_vbs_info(:) = -99
           msoa_vbs_info(1) = vbs_nbin(1)

           if (vbs_nbin(1) > 0) then  
               msoa_flag1 = 1000
               msoa_vbs_info(2) = vbs_uq_aqsoa
               msoa_vbs_info(3) = vbs_uq_par
           else
               msoa_flag1 = 1
               msoa_flag1 = 0  
           end if

           print '(/a,5i9)', 'msoa_flag1, msoa_vbs_info(1:3)', msoa_flag1, msoa_vbs_info(1:3)
           return
        end if



















	use_cam5mam_soa_params  = 0  
	use_cam5mam_accom_coefs = 0  




	mhyst_method    = mhyst_uporlo_waterhyst  
	mGAS_AER_XFER   = 1    
	mDYNAMIC_SOLVER = 1    

	msize_framework = msectional  
	alpha_ASTEM     = 0.05 
	rtol_eqb_ASTEM  = 0.01 
	ptol_mol_ASTEM  = 0.01 
	ipmcmos         = 0    
	m_partmc_mosaic = 0    









	method_bcfrac       = 1      
	method_kappa        = 11     
	maersize_init_flag1 = 0      

	ifreq_coag          = 1      
	if (ntype_aer <= 1) then
	    mcoag_flag1      = 10     
	    mmovesect_flag1  = 20     
	else
	    mcoag_flag1      = 60     
	    mmovesect_flag1  = 60     
	end if
	mnewnuc_flag1       = 3      
	                             
	msectional_flag1    = 1      

	msectional_flag2 = 0
	if (ntype_aer > 1) msectional_flag2 = 1

	mgas                = 0      
	maer                = 1      
	mcld                = 0      
	maeroptic           = 0      
	mshellcore          = 0      
	msolar              = 0      
	mphoto              = 0      

	call mosaic_init_aer_params

	end subroutine mosaic2_wrfchem_init



        subroutine mosaic2_output_wrfchem_pointers( id )

        use module_data_mosaic_asecthp
        use module_peg_util, only:  peg_error_fatal, peg_message
        use module_scalar_tables, only:  chem_dname_table

        implicit none


        integer, intent(in) :: id


        integer l, ll, lu, n, ns
        integer icomp, isize, itype, iphase, jt, jp

        character*200 msg




9350    format( a, 32(1x,i4) )
	msg = ' '
	call peg_message( lunout, msg )
	msg = 'output from subr mosaic2_output_wrfchem_pointers'
	call peg_message( lunout, trim(msg) )
	write(msg,9350) 'nphase_aer =     ', nphase_aer
	call peg_message( lunout, trim(msg) )

	do iphase=1,nphase_aer

	write(msg,9350) 'iphase =     ', iphase
	call peg_message( lunout, trim(msg) )
	write(msg,9350) 'ntype_aer =     ', ntype_aer
	call peg_message( lunout, trim(msg) )
        write(msg,9350) 'ncomp_aer =     ', ncomp_aer
        call peg_message( lunout, trim(msg) )
      
	do itype=1,ntype_aer

	write(msg,9350) 'itype =     ', itype
	call peg_message( lunout, trim(msg) )
	write(msg,9350) 'nsize_aer = ', nsize_aer(itype)
	call peg_message( lunout, trim(msg) )

        jt = itype
        jp = iphase
        ns = nsize_aer(itype)
        lu = lunout

        call mosaic2_out1ptraa( "numptr_aer         ", numptr_aer(1:ns,jt,jp),         ns, lu )
        call mosaic2_out1ptraa( "hyswptr_aer        ", hyswptr_aer(1:ns,jt),           ns, lu )
        call mosaic2_out1ptraa( "waterptr_aer       ", waterptr_aer(1:ns,jt),          ns, lu )

        call mosaic2_out1ptraa( "lptr_so4_aer       ", lptr_so4_aer(1:ns,jt,jp),       ns, lu )
        call mosaic2_out1ptraa( "lptr_no3_aer       ", lptr_no3_aer(1:ns,jt,jp),       ns, lu )
        call mosaic2_out1ptraa( "lptr_cl_aer        ", lptr_cl_aer(1:ns,jt,jp),        ns, lu )
        call mosaic2_out1ptraa( "lptr_msa_aer       ", lptr_msa_aer(1:ns,jt,jp),       ns, lu )
        call mosaic2_out1ptraa( "lptr_co3_aer       ", lptr_co3_aer(1:ns,jt,jp),       ns, lu )
        call mosaic2_out1ptraa( "lptr_nh4_aer       ", lptr_nh4_aer(1:ns,jt,jp),       ns, lu )
        call mosaic2_out1ptraa( "lptr_na_aer        ", lptr_na_aer(1:ns,jt,jp),        ns, lu )
        call mosaic2_out1ptraa( "lptr_ca_aer        ", lptr_ca_aer(1:ns,jt,jp),        ns, lu )
        call mosaic2_out1ptraa( "lptr_oin_aer       ", lptr_oin_aer(1:ns,jt,jp),       ns, lu )
        call mosaic2_out1ptraa( "lptr_aro1_aer      ", lptr_aro1_aer(1:ns,jt,jp),      ns, lu )
        call mosaic2_out1ptraa( "lptr_aro2_aer      ", lptr_aro2_aer(1:ns,jt,jp),      ns, lu )
        call mosaic2_out1ptraa( "lptr_alk1_aer      ", lptr_alk1_aer(1:ns,jt,jp),      ns, lu )
        call mosaic2_out1ptraa( "lptr_ole1_aer      ", lptr_ole1_aer(1:ns,jt,jp),      ns, lu )
        call mosaic2_out1ptraa( "lptr_api1_aer      ", lptr_api1_aer(1:ns,jt,jp),      ns, lu )
        call mosaic2_out1ptraa( "lptr_api2_aer      ", lptr_api2_aer(1:ns,jt,jp),      ns, lu )
        call mosaic2_out1ptraa( "lptr_lim1_aer      ", lptr_lim1_aer(1:ns,jt,jp),      ns, lu )
        call mosaic2_out1ptraa( "lptr_lim2_aer      ", lptr_lim2_aer(1:ns,jt,jp),      ns, lu )


        call mosaic2_out1ptraa( "lptr_oc_aer        ", lptr_oc_aer(1:ns,jt,jp),        ns, lu )
        call mosaic2_out1ptraa( "lptr_bc_aer        ", lptr_bc_aer(1:ns,jt,jp),        ns, lu )





        call mosaic2_out1ptraa( "lptr_pcg1_b_c_aer  ", lptr_pcg1_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg2_b_c_aer  ", lptr_pcg2_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg3_b_c_aer  ", lptr_pcg3_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg4_b_c_aer  ", lptr_pcg4_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg5_b_c_aer  ", lptr_pcg5_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg6_b_c_aer  ", lptr_pcg6_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg7_b_c_aer  ", lptr_pcg7_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg8_b_c_aer  ", lptr_pcg8_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg9_b_c_aer  ", lptr_pcg9_b_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg1_b_o_aer  ", lptr_pcg1_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg2_b_o_aer  ", lptr_pcg2_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg3_b_o_aer  ", lptr_pcg3_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg4_b_o_aer  ", lptr_pcg4_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg5_b_o_aer  ", lptr_pcg5_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg6_b_o_aer  ", lptr_pcg6_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg7_b_o_aer  ", lptr_pcg7_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg8_b_o_aer  ", lptr_pcg8_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg9_b_o_aer  ", lptr_pcg9_b_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg1_b_c_aer ", lptr_opcg1_b_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg2_b_c_aer ", lptr_opcg2_b_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg3_b_c_aer ", lptr_opcg3_b_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg4_b_c_aer ", lptr_opcg4_b_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg5_b_c_aer ", lptr_opcg5_b_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg6_b_c_aer ", lptr_opcg6_b_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg7_b_c_aer ", lptr_opcg7_b_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg8_b_c_aer ", lptr_opcg8_b_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg1_b_o_aer ", lptr_opcg1_b_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg2_b_o_aer ", lptr_opcg2_b_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg3_b_o_aer ", lptr_opcg3_b_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg4_b_o_aer ", lptr_opcg4_b_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg5_b_o_aer ", lptr_opcg5_b_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg6_b_o_aer ", lptr_opcg6_b_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg7_b_o_aer ", lptr_opcg7_b_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg8_b_o_aer ", lptr_opcg8_b_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg1_f_c_aer  ", lptr_pcg1_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg2_f_c_aer  ", lptr_pcg2_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg3_f_c_aer  ", lptr_pcg3_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg4_f_c_aer  ", lptr_pcg4_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg5_f_c_aer  ", lptr_pcg5_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg6_f_c_aer  ", lptr_pcg6_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg7_f_c_aer  ", lptr_pcg7_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg8_f_c_aer  ", lptr_pcg8_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg9_f_c_aer  ", lptr_pcg9_f_c_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg1_f_o_aer  ", lptr_pcg1_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg2_f_o_aer  ", lptr_pcg2_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg3_f_o_aer  ", lptr_pcg3_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg4_f_o_aer  ", lptr_pcg4_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg5_f_o_aer  ", lptr_pcg5_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg6_f_o_aer  ", lptr_pcg6_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg7_f_o_aer  ", lptr_pcg7_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg8_f_o_aer  ", lptr_pcg8_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_pcg9_f_o_aer  ", lptr_pcg9_f_o_aer(1:ns,jt,jp),  ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg1_f_c_aer ", lptr_opcg1_f_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg2_f_c_aer ", lptr_opcg2_f_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg3_f_c_aer ", lptr_opcg3_f_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg4_f_c_aer ", lptr_opcg4_f_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg5_f_c_aer ", lptr_opcg5_f_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg6_f_c_aer ", lptr_opcg6_f_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg7_f_c_aer ", lptr_opcg7_f_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg8_f_c_aer ", lptr_opcg8_f_c_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg1_f_o_aer ", lptr_opcg1_f_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg2_f_o_aer ", lptr_opcg2_f_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg3_f_o_aer ", lptr_opcg3_f_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg4_f_o_aer ", lptr_opcg4_f_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg5_f_o_aer ", lptr_opcg5_f_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg6_f_o_aer ", lptr_opcg6_f_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg7_f_o_aer ", lptr_opcg7_f_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "lptr_opcg8_f_o_aer ", lptr_opcg8_f_o_aer(1:ns,jt,jp), ns, lu )
        call mosaic2_out1ptraa( "ant1_c_aer         ", lptr_ant1_c_aer(1:ns,jt,jp),    ns, lu )
        call mosaic2_out1ptraa( "ant2_c_aer         ", lptr_ant2_c_aer(1:ns,jt,jp),    ns, lu )
        call mosaic2_out1ptraa( "ant3_c_aer         ", lptr_ant3_c_aer(1:ns,jt,jp),    ns, lu )
        call mosaic2_out1ptraa( "ant4_c_aer         ", lptr_ant4_c_aer(1:ns,jt,jp),    ns, lu )
        call mosaic2_out1ptraa( "ant1_o_aer         ", lptr_ant1_o_aer(1:ns,jt,jp),    ns, lu )
        call mosaic2_out1ptraa( "ant2_o_aer         ", lptr_ant2_o_aer(1:ns,jt,jp),    ns, lu )
        call mosaic2_out1ptraa( "ant3_o_aer         ", lptr_ant3_o_aer(1:ns,jt,jp),    ns, lu )
        call mosaic2_out1ptraa( "ant4_o_aer         ", lptr_ant4_o_aer(1:ns,jt,jp),    ns, lu )
        call mosaic2_out1ptraa( "biog1_c_aer        ", lptr_biog1_c_aer(1:ns,jt,jp),   ns, lu )
        call mosaic2_out1ptraa( "biog2_c_aer        ", lptr_biog2_c_aer(1:ns,jt,jp),   ns, lu )
        call mosaic2_out1ptraa( "biog3_c_aer        ", lptr_biog3_c_aer(1:ns,jt,jp),   ns, lu )
        call mosaic2_out1ptraa( "biog4_c_aer        ", lptr_biog4_c_aer(1:ns,jt,jp),   ns, lu )
        call mosaic2_out1ptraa( "biog1_o_aer        ", lptr_biog1_o_aer(1:ns,jt,jp),   ns, lu )
        call mosaic2_out1ptraa( "biog2_o_aer        ", lptr_biog2_o_aer(1:ns,jt,jp),   ns, lu )
        call mosaic2_out1ptraa( "biog3_o_aer        ", lptr_biog3_o_aer(1:ns,jt,jp),   ns, lu )
        call mosaic2_out1ptraa( "biog4_o_aer        ", lptr_biog4_o_aer(1:ns,jt,jp),   ns, lu )





9352    format( 5a )
	do ll = 1, ncomp_plustracer_aer(itype)
	    write(msg,9350) 'massptr_aer(), ll',    &
		(massptr_aer(ll,n,itype,iphase), n=1,nsize_aer(itype)), ll
	    call peg_message( lunout, trim(msg) )

            write(msg,9352) 'chem_dname_table =  ', &
                (chem_dname_table(id,(massptr_aer(ll,n,itype,iphase)))(1:14),n=1,4)
	    call peg_message( lunout, trim(msg) )
	end do 

	end do 
	end do 

        return
        end subroutine mosaic2_output_wrfchem_pointers



	subroutine mosaic2_out1ptraa( lptrname, lptr, nsz, lunout )

	use module_peg_util, only:  peg_message

        integer, intent(in) :: nsz, lunout, lptr(1:nsz)
        character(len=*), intent(in) :: lptrname
        
	character*200 msg

9350    format( a, 32(1x,i4) )
        write(msg,9350) lptrname, lptr(1:nsz)
        call peg_message( lunout, trim(msg) )

        return
        end subroutine mosaic2_out1ptraa



	subroutine mosaic2_set_gas_kptrs_names



	use module_state_description, only:  num_chem, p1st => param_first_scalar
        use module_configure, only:   p_h2so4, p_sulf
	use module_peg_util, only:  peg_error_fatal, peg_message

	use module_data_mosaic_main, only:  ntot_used
	use module_data_mosaic_asecthp, only:   &
            lunout, lunerr
	use module_data_mosaic_boxmod, only:   name_rbox, &
            kh2so4, khno3,  khcl,   knh3,   kmsa,   ko3, &
            kso2,   kh2o2,  khcho,  koh,    kho2,        &
            kno3,   kno,    kno2,   khono,  kpan,        &
            kn2o5,  kclno2,                              &
            karo1,     karo2,     kalk1,     kole1,                &
            kapi1,     kapi2,     klim1,     klim2,                &
            kpcg1_b_c, kpcg2_b_c, kpcg3_b_c, kpcg4_b_c, kpcg5_b_c, &
            kpcg6_b_c, kpcg7_b_c, kpcg8_b_c, kpcg9_b_c,            &
            kpcg1_b_o, kpcg2_b_o, kpcg3_b_o, kpcg4_b_o, kpcg5_b_o, &
            kpcg6_b_o, kpcg7_b_o, kpcg8_b_o, kpcg9_b_o,            &
            kpcg1_f_c, kpcg2_f_c, kpcg3_f_c, kpcg4_f_c, kpcg5_f_c, &
            kpcg6_f_c, kpcg7_f_c, kpcg8_f_c, kpcg9_f_c,            &
            kpcg1_f_o, kpcg2_f_o, kpcg3_f_o, kpcg4_f_o, kpcg5_f_o, &
            kpcg6_f_o, kpcg7_f_o, kpcg8_f_o, kpcg9_f_o,            &
            kopcg1_b_c, kopcg2_b_c, kopcg3_b_c, kopcg4_b_c, &
            kopcg5_b_c, kopcg6_b_c, kopcg7_b_c, kopcg8_b_c, &
            kopcg1_b_o, kopcg2_b_o, kopcg3_b_o, kopcg4_b_o, &
            kopcg5_b_o, kopcg6_b_o, kopcg7_b_o, kopcg8_b_o, &
            kopcg1_f_c, kopcg2_f_c, kopcg3_f_c, kopcg4_f_c, &
            kopcg5_f_c, kopcg6_f_c, kopcg7_f_c, kopcg8_f_c, &
            kopcg1_f_o, kopcg2_f_o, kopcg3_f_o, kopcg4_f_o, &
            kopcg5_f_o, kopcg6_f_o, kopcg7_f_o, kopcg8_f_o, &
            kant1_c,  kant2_c,  kant3_c,  kant4_c, &
            kant1_o,  kant2_o,  kant3_o,  kant4_o, &
            kbiog1_c, kbiog2_c, kbiog3_c, kbiog4_c, &
            kbiog1_o, kbiog2_o, kbiog3_o, kbiog4_o, &
            ksmpa, ksmpbb, &

            khcooh, kch3o2, kch3oh, kch3ooh

	implicit none




	integer l, ll, n

	character*200 msg

	msg = ' '
	call peg_message( lunout, msg )
	msg = 'output from subr mosaic2_set_gas_kptrs_names'
	call peg_message( lunout, trim(msg) )




        if (p_h2so4 .ge. p1st .and. p_h2so4 .le. ntot_used) then 
            kh2so4 = p_h2so4
       
        else if (p_sulf .ge. p1st .and. p_sulf .le. ntot_used) then 
            kh2so4 = p_sulf
	else
	    msg = '*** subr mosaic2_set_gas_kptrs_names - ptr error for h2so4'
	    call peg_error_fatal( lunerr, msg )
	end if

        if ((kh2so4 .ge. p1st) .and. (kh2so4 .le. ntot_used)) then
             name_rbox(kh2so4 ) = 'h2so4'
        endif






        call mosaic2_set_one_gas_kptr_name( "ho", l )
        if ((l .ge. p1st) .and. (l .le. ntot_used)) then
            koh = l
            name_rbox(koh ) = 'oh'
        endif

        call mosaic2_set_one_gas_kptr_name( "hno3",      khno3      )
        call mosaic2_set_one_gas_kptr_name( "nh3",       knh3       )
        call mosaic2_set_one_gas_kptr_name( "hcl",       khcl       )
        call mosaic2_set_one_gas_kptr_name( "msa",       kmsa       )
        call mosaic2_set_one_gas_kptr_name( "aro1",      karo1      )
        call mosaic2_set_one_gas_kptr_name( "aro2",      karo2      )
        call mosaic2_set_one_gas_kptr_name( "alk1",      kalk1      )
        call mosaic2_set_one_gas_kptr_name( "ole1",      kole1      )
        call mosaic2_set_one_gas_kptr_name( "api1",      kapi1      )
        call mosaic2_set_one_gas_kptr_name( "api2",      kapi2      )
        call mosaic2_set_one_gas_kptr_name( "lim1",      klim1      )
        call mosaic2_set_one_gas_kptr_name( "lim2",      klim2      )
        call mosaic2_set_one_gas_kptr_name( "n2o5",      kn2o5      )
        call mosaic2_set_one_gas_kptr_name( "clno2",     kclno2     )
        call mosaic2_set_one_gas_kptr_name( "o3",        ko3        )
        call mosaic2_set_one_gas_kptr_name( "pcg1_b_c",  kpcg1_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg2_b_c",  kpcg2_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg3_b_c",  kpcg3_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg4_b_c",  kpcg4_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg5_b_c",  kpcg5_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg6_b_c",  kpcg6_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg7_b_c",  kpcg7_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg8_b_c",  kpcg8_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg9_b_c",  kpcg9_b_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg1_b_o",  kpcg1_b_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg2_b_o",  kpcg2_b_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg3_b_o",  kpcg3_b_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg4_b_o",  kpcg4_b_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg5_b_o",  kpcg5_b_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg6_b_o",  kpcg6_b_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg7_b_o",  kpcg7_b_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg8_b_o",  kpcg8_b_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg9_b_o",  kpcg9_b_o  )
        call mosaic2_set_one_gas_kptr_name( "opcg1_b_c", kopcg1_b_c )
        call mosaic2_set_one_gas_kptr_name( "opcg2_b_c", kopcg2_b_c )
        call mosaic2_set_one_gas_kptr_name( "opcg3_b_c", kopcg3_b_c )
        call mosaic2_set_one_gas_kptr_name( "opcg4_b_c", kopcg4_b_c )
        call mosaic2_set_one_gas_kptr_name( "opcg5_b_c", kopcg5_b_c )
        call mosaic2_set_one_gas_kptr_name( "opcg6_b_c", kopcg6_b_c )
        call mosaic2_set_one_gas_kptr_name( "opcg7_b_c", kopcg7_b_c )
        call mosaic2_set_one_gas_kptr_name( "opcg8_b_c", kopcg8_b_c )
        call mosaic2_set_one_gas_kptr_name( "opcg1_b_o", kopcg1_b_o )
        call mosaic2_set_one_gas_kptr_name( "opcg2_b_o", kopcg2_b_o )
        call mosaic2_set_one_gas_kptr_name( "opcg3_b_o", kopcg3_b_o )
        call mosaic2_set_one_gas_kptr_name( "opcg4_b_o", kopcg4_b_o )
        call mosaic2_set_one_gas_kptr_name( "opcg5_b_o", kopcg5_b_o )
        call mosaic2_set_one_gas_kptr_name( "opcg6_b_o", kopcg6_b_o )
        call mosaic2_set_one_gas_kptr_name( "opcg7_b_o", kopcg7_b_o )
        call mosaic2_set_one_gas_kptr_name( "opcg8_b_o", kopcg8_b_o )
        call mosaic2_set_one_gas_kptr_name( "pcg1_f_c",  kpcg1_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg2_f_c",  kpcg2_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg3_f_c",  kpcg3_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg4_f_c",  kpcg4_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg5_f_c",  kpcg5_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg6_f_c",  kpcg6_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg7_f_c",  kpcg7_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg8_f_c",  kpcg8_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg9_f_c",  kpcg9_f_c  )
        call mosaic2_set_one_gas_kptr_name( "pcg1_f_o",  kpcg1_f_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg2_f_o",  kpcg2_f_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg3_f_o",  kpcg3_f_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg4_f_o",  kpcg4_f_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg5_f_o",  kpcg5_f_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg6_f_o",  kpcg6_f_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg7_f_o",  kpcg7_f_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg8_f_o",  kpcg8_f_o  )
        call mosaic2_set_one_gas_kptr_name( "pcg9_f_o",  kpcg9_f_o  )
        call mosaic2_set_one_gas_kptr_name( "opcg1_f_c", kopcg1_f_c )
        call mosaic2_set_one_gas_kptr_name( "opcg2_f_c", kopcg2_f_c )
        call mosaic2_set_one_gas_kptr_name( "opcg3_f_c", kopcg3_f_c )
        call mosaic2_set_one_gas_kptr_name( "opcg4_f_c", kopcg4_f_c )
        call mosaic2_set_one_gas_kptr_name( "opcg5_f_c", kopcg5_f_c )
        call mosaic2_set_one_gas_kptr_name( "opcg6_f_c", kopcg6_f_c )
        call mosaic2_set_one_gas_kptr_name( "opcg7_f_c", kopcg7_f_c )
        call mosaic2_set_one_gas_kptr_name( "opcg8_f_c", kopcg8_f_c )
        call mosaic2_set_one_gas_kptr_name( "opcg1_f_o", kopcg1_f_o )
        call mosaic2_set_one_gas_kptr_name( "opcg2_f_o", kopcg2_f_o )
        call mosaic2_set_one_gas_kptr_name( "opcg3_f_o", kopcg3_f_o )
        call mosaic2_set_one_gas_kptr_name( "opcg4_f_o", kopcg4_f_o )
        call mosaic2_set_one_gas_kptr_name( "opcg5_f_o", kopcg5_f_o )
        call mosaic2_set_one_gas_kptr_name( "opcg6_f_o", kopcg6_f_o )
        call mosaic2_set_one_gas_kptr_name( "opcg7_f_o", kopcg7_f_o )
        call mosaic2_set_one_gas_kptr_name( "opcg8_f_o", kopcg8_f_o )
        call mosaic2_set_one_gas_kptr_name( "smpa",      ksmpa      )
        call mosaic2_set_one_gas_kptr_name( "smpbb",     ksmpbb     )


        call mosaic2_set_one_gas_kptr_name( "ant1_c",    kant1_c    )
        call mosaic2_set_one_gas_kptr_name( "ant2_c",    kant2_c    )
        call mosaic2_set_one_gas_kptr_name( "ant3_c",    kant3_c    )
        call mosaic2_set_one_gas_kptr_name( "ant4_c",    kant4_c    )
        call mosaic2_set_one_gas_kptr_name( "ant1_o",    kant1_o    )
        call mosaic2_set_one_gas_kptr_name( "ant2_o",    kant2_o    )
        call mosaic2_set_one_gas_kptr_name( "ant3_o",    kant3_o    )
        call mosaic2_set_one_gas_kptr_name( "ant4_o",    kant4_o    )
        call mosaic2_set_one_gas_kptr_name( "biog1_c",   kbiog1_c   )
        call mosaic2_set_one_gas_kptr_name( "biog2_c",   kbiog2_c   )
        call mosaic2_set_one_gas_kptr_name( "biog3_c",   kbiog3_c   )
        call mosaic2_set_one_gas_kptr_name( "biog4_c",   kbiog4_c   )
        call mosaic2_set_one_gas_kptr_name( "biog1_o",   kbiog1_o   )
        call mosaic2_set_one_gas_kptr_name( "biog2_o",   kbiog2_o   )
        call mosaic2_set_one_gas_kptr_name( "biog3_o",   kbiog3_o   )
        call mosaic2_set_one_gas_kptr_name( "biog4_o",   kbiog4_o   )
        call mosaic2_set_one_gas_kptr_name( "so2",       kso2       )
        call mosaic2_set_one_gas_kptr_name( "h2o2",      kh2o2      )
        call mosaic2_set_one_gas_kptr_name( "hcho",      khcho      )
        call mosaic2_set_one_gas_kptr_name( "hcooh",     khcooh     )
        call mosaic2_set_one_gas_kptr_name( "ho2",       kho2       )
        call mosaic2_set_one_gas_kptr_name( "no3",       kno3       )
        call mosaic2_set_one_gas_kptr_name( "no",        kno        )
        call mosaic2_set_one_gas_kptr_name( "no2",       kno2       )
        call mosaic2_set_one_gas_kptr_name( "hono",      khono      )
        call mosaic2_set_one_gas_kptr_name( "pan",       kpan       )
        call mosaic2_set_one_gas_kptr_name( "ch3o2",     kch3o2     )
        call mosaic2_set_one_gas_kptr_name( "ch3oh",     kch3oh     )
        call mosaic2_set_one_gas_kptr_name( "ch3ooh",    kch3ooh    )

	return
	end subroutine mosaic2_set_gas_kptrs_names



        subroutine mosaic2_set_one_gas_kptr_name( gasname, kptr )

	use module_data_mosaic_asecthp, only:  lunerr, lunout
	use module_data_mosaic_boxmod, only:  name_rbox
        use module_data_mosaic_main, only:  ntot_used
        use module_state_description, only:  p1st => param_first_scalar
	use module_peg_util, only:  peg_error_fatal, peg_message

        integer, intent(inout) :: kptr
        character(len=*), intent(in) :: gasname

        integer :: l
	character(len=40)  :: chemname
	character(len=200) :: msg

        chemname = gasname
	call mosaic2_find_chemname_in_table( chemname, l )

        if (l > ntot_used) then
	    write(msg,'(2a,2(1x,i5),2x,a)') 'mosaic2_set_one_gas_kptr_name error', &
                ' - l, ntot_used, gasname = ', l, ntot_used, trim(gasname)
	    call peg_error_fatal( lunerr, msg )
        end if

        if (l >= p1st) then
            kptr = l
            name_rbox(l) = gasname
        end if

	write(msg,'(2i12,2a)') l, kptr, ' = k', trim(chemname)
	call peg_message( lunout, trim(msg) )

        return
        end subroutine mosaic2_set_one_gas_kptr_name



        subroutine set_rbox_gas_ptrs




        use module_state_description, only:  p1st => param_first_scalar

        use module_data_mosaic_aero, only: &
            gas_name, ngas_aerchtot, &
            ih2so4_g,     ihno3_g,      ihcl_g,      inh3_g,        &
            imsa_g,                                                 &
            iaro1_g,      iaro2_g,      ialk1_g,     iole1_g,       &
            iapi1_g,      iapi2_g,      ilim1_g,     ilim2_g,       &
            in2o5_g,      iclno2_g,                                 &
            ipcg1_b_c_g,  ipcg2_b_c_g,  ipcg3_b_c_g,  ipcg4_b_c_g,  &
            ipcg5_b_c_g,  ipcg6_b_c_g,  ipcg7_b_c_g,  ipcg8_b_c_g,  &
            ipcg9_b_c_g,                                            &
            ipcg1_b_o_g,  ipcg2_b_o_g,  ipcg3_b_o_g,  ipcg4_b_o_g,  &
            ipcg5_b_o_g,  ipcg6_b_o_g,  ipcg7_b_o_g,  ipcg8_b_o_g,  &
            ipcg9_b_o_g,                                            &
            iopcg1_b_c_g, iopcg2_b_c_g, iopcg3_b_c_g, iopcg4_b_c_g, &
            iopcg5_b_c_g, iopcg6_b_c_g, iopcg7_b_c_g, iopcg8_b_c_g, &
            iopcg1_b_o_g, iopcg2_b_o_g, iopcg3_b_o_g, iopcg4_b_o_g, &
            iopcg5_b_o_g, iopcg6_b_o_g, iopcg7_b_o_g, iopcg8_b_o_g, &
            ipcg1_f_c_g,  ipcg2_f_c_g,  ipcg3_f_c_g,  ipcg4_f_c_g,  &
            ipcg5_f_c_g,  ipcg6_f_c_g,  ipcg7_f_c_g,  ipcg8_f_c_g,  &
            ipcg9_f_c_g,                                            &
            ipcg1_f_o_g,  ipcg2_f_o_g,  ipcg3_f_o_g,  ipcg4_f_o_g,  &
            ipcg5_f_o_g,  ipcg6_f_o_g,  ipcg7_f_o_g,  ipcg8_f_o_g,  &
            ipcg9_f_o_g,                                            &
            iopcg1_f_c_g, iopcg2_f_c_g, iopcg3_f_c_g, iopcg4_f_c_g, &
            iopcg5_f_c_g, iopcg6_f_c_g, iopcg7_f_c_g, iopcg8_f_c_g, &
            iopcg1_f_o_g, iopcg2_f_o_g, iopcg3_f_o_g, iopcg4_f_o_g, &
            iopcg5_f_o_g, iopcg6_f_o_g, iopcg7_f_o_g, iopcg8_f_o_g, &
            iant1_c_g,    iant2_c_g,    iant3_c_g,    iant4_c_g,    &
            iant1_o_g,    iant2_o_g,    iant3_o_g,    iant4_o_g,    &
            ibiog1_c_g,   ibiog2_c_g,   ibiog3_c_g,   ibiog4_c_g,   &
            ibiog1_o_g,   ibiog2_o_g,   ibiog3_o_g,   ibiog4_o_g,   &
            ismpa_g,      ismpbb_g



        use module_data_mosaic_asecthp, only:  rbox_gas_ptr

        use module_data_mosaic_boxmod, only: &
            name_rbox, &
            kh2so4, khno3,  khcl,   knh3,   kmsa,   ko3, &
            kso2,   kh2o2,  khcho,  koh,    kho2,        &
            kno3,   kno,    kno2,   khono,  kpan,        &
            kn2o5,  kclno2,                              &
            karo1,     karo2,     kalk1,     kole1,                &
            kapi1,     kapi2,     klim1,     klim2,                &
            kpcg1_b_c, kpcg2_b_c, kpcg3_b_c, kpcg4_b_c, kpcg5_b_c, &
            kpcg6_b_c, kpcg7_b_c, kpcg8_b_c, kpcg9_b_c,            &
            kpcg1_b_o, kpcg2_b_o, kpcg3_b_o, kpcg4_b_o, kpcg5_b_o, &
            kpcg6_b_o, kpcg7_b_o, kpcg8_b_o, kpcg9_b_o,            &
            kpcg1_f_c, kpcg2_f_c, kpcg3_f_c, kpcg4_f_c, kpcg5_f_c, &
            kpcg6_f_c, kpcg7_f_c, kpcg8_f_c, kpcg9_f_c,            &
            kpcg1_f_o, kpcg2_f_o, kpcg3_f_o, kpcg4_f_o, kpcg5_f_o, &
            kpcg6_f_o, kpcg7_f_o, kpcg8_f_o, kpcg9_f_o,            &
            kopcg1_b_c, kopcg2_b_c, kopcg3_b_c, kopcg4_b_c, &
            kopcg5_b_c, kopcg6_b_c, kopcg7_b_c, kopcg8_b_c, &
            kopcg1_b_o, kopcg2_b_o, kopcg3_b_o, kopcg4_b_o, &
            kopcg5_b_o, kopcg6_b_o, kopcg7_b_o, kopcg8_b_o, &
            kopcg1_f_c, kopcg2_f_c, kopcg3_f_c, kopcg4_f_c, &
            kopcg5_f_c, kopcg6_f_c, kopcg7_f_c, kopcg8_f_c, &
            kopcg1_f_o, kopcg2_f_o, kopcg3_f_o, kopcg4_f_o, &
            kopcg5_f_o, kopcg6_f_o, kopcg7_f_o, kopcg8_f_o, &
            kant1_c,  kant2_c,  kant3_c,  kant4_c, &
            kant1_o,  kant2_o,  kant3_o,  kant4_o, &
            kbiog1_c, kbiog2_c, kbiog3_c, kbiog4_c, &
            kbiog1_o, kbiog2_o, kbiog3_o, kbiog4_o, &
            ksmpa, ksmpbb, &

            khcooh, kch3o2, kch3oh, kch3ooh

        use module_data_mosaic_main, only: &
            ntot_used
        
        integer :: igas, l


        rbox_gas_ptr(1:ngas_aerchtot) = -1

        call set_1_rbox_gas_ptr( ih2so4_g,     kh2so4     )
        call set_1_rbox_gas_ptr( ihno3_g,      khno3      )
        call set_1_rbox_gas_ptr( ihcl_g,       khcl       )
        call set_1_rbox_gas_ptr( inh3_g,       knh3       )
        call set_1_rbox_gas_ptr( imsa_g,       kmsa       )
        call set_1_rbox_gas_ptr( iaro1_g,      karo1      )
        call set_1_rbox_gas_ptr( iaro2_g,      karo2      )
        call set_1_rbox_gas_ptr( ialk1_g,      kalk1      )
        call set_1_rbox_gas_ptr( iole1_g,      kole1      )
        call set_1_rbox_gas_ptr( iapi1_g,      kapi1      )
        call set_1_rbox_gas_ptr( iapi2_g,      kapi2      )
        call set_1_rbox_gas_ptr( ilim1_g,      klim1      )
        call set_1_rbox_gas_ptr( ilim2_g,      klim2      )
        call set_1_rbox_gas_ptr( in2o5_g,      kn2o5      )
        call set_1_rbox_gas_ptr( iclno2_g,     kclno2     )
        call set_1_rbox_gas_ptr( ipcg1_b_c_g,  kpcg1_b_c  )
        call set_1_rbox_gas_ptr( ipcg2_b_c_g,  kpcg2_b_c  )
        call set_1_rbox_gas_ptr( ipcg3_b_c_g,  kpcg3_b_c  )
        call set_1_rbox_gas_ptr( ipcg4_b_c_g,  kpcg4_b_c  )
        call set_1_rbox_gas_ptr( ipcg5_b_c_g,  kpcg5_b_c  )
        call set_1_rbox_gas_ptr( ipcg6_b_c_g,  kpcg6_b_c  )
        call set_1_rbox_gas_ptr( ipcg7_b_c_g,  kpcg7_b_c  )
        call set_1_rbox_gas_ptr( ipcg8_b_c_g,  kpcg8_b_c  )
        call set_1_rbox_gas_ptr( ipcg9_b_c_g,  kpcg9_b_c  )
        call set_1_rbox_gas_ptr( ipcg1_b_o_g,  kpcg1_b_o  )
        call set_1_rbox_gas_ptr( ipcg2_b_o_g,  kpcg2_b_o  )
        call set_1_rbox_gas_ptr( ipcg3_b_o_g,  kpcg3_b_o  )
        call set_1_rbox_gas_ptr( ipcg4_b_o_g,  kpcg4_b_o  )
        call set_1_rbox_gas_ptr( ipcg5_b_o_g,  kpcg5_b_o  )
        call set_1_rbox_gas_ptr( ipcg6_b_o_g,  kpcg6_b_o  )
        call set_1_rbox_gas_ptr( ipcg7_b_o_g,  kpcg7_b_o  )
        call set_1_rbox_gas_ptr( ipcg8_b_o_g,  kpcg8_b_o  )
        call set_1_rbox_gas_ptr( ipcg9_b_o_g,  kpcg9_b_o  )
        call set_1_rbox_gas_ptr( iopcg1_b_c_g, kopcg1_b_c )
        call set_1_rbox_gas_ptr( iopcg2_b_c_g, kopcg2_b_c )
        call set_1_rbox_gas_ptr( iopcg3_b_c_g, kopcg3_b_c )
        call set_1_rbox_gas_ptr( iopcg4_b_c_g, kopcg4_b_c )
        call set_1_rbox_gas_ptr( iopcg5_b_c_g, kopcg5_b_c )
        call set_1_rbox_gas_ptr( iopcg6_b_c_g, kopcg6_b_c )
        call set_1_rbox_gas_ptr( iopcg7_b_c_g, kopcg7_b_c )
        call set_1_rbox_gas_ptr( iopcg8_b_c_g, kopcg8_b_c )
        call set_1_rbox_gas_ptr( iopcg1_b_o_g, kopcg1_b_o )
        call set_1_rbox_gas_ptr( iopcg2_b_o_g, kopcg2_b_o )
        call set_1_rbox_gas_ptr( iopcg3_b_o_g, kopcg3_b_o )
        call set_1_rbox_gas_ptr( iopcg4_b_o_g, kopcg4_b_o )
        call set_1_rbox_gas_ptr( iopcg5_b_o_g, kopcg5_b_o )
        call set_1_rbox_gas_ptr( iopcg6_b_o_g, kopcg6_b_o )
        call set_1_rbox_gas_ptr( iopcg7_b_o_g, kopcg7_b_o )
        call set_1_rbox_gas_ptr( iopcg8_b_o_g, kopcg8_b_o )
        call set_1_rbox_gas_ptr( ipcg1_f_c_g,  kpcg1_f_c  )
        call set_1_rbox_gas_ptr( ipcg2_f_c_g,  kpcg2_f_c  )
        call set_1_rbox_gas_ptr( ipcg3_f_c_g,  kpcg3_f_c  )
        call set_1_rbox_gas_ptr( ipcg4_f_c_g,  kpcg4_f_c  )
        call set_1_rbox_gas_ptr( ipcg5_f_c_g,  kpcg5_f_c  )
        call set_1_rbox_gas_ptr( ipcg6_f_c_g,  kpcg6_f_c  )
        call set_1_rbox_gas_ptr( ipcg7_f_c_g,  kpcg7_f_c  )
        call set_1_rbox_gas_ptr( ipcg8_f_c_g,  kpcg8_f_c  )
        call set_1_rbox_gas_ptr( ipcg9_f_c_g,  kpcg9_f_c  )
        call set_1_rbox_gas_ptr( ipcg1_f_o_g,  kpcg1_f_o  )
        call set_1_rbox_gas_ptr( ipcg2_f_o_g,  kpcg2_f_o  )
        call set_1_rbox_gas_ptr( ipcg3_f_o_g,  kpcg3_f_o  )
        call set_1_rbox_gas_ptr( ipcg4_f_o_g,  kpcg4_f_o  )
        call set_1_rbox_gas_ptr( ipcg5_f_o_g,  kpcg5_f_o  )
        call set_1_rbox_gas_ptr( ipcg6_f_o_g,  kpcg6_f_o  )
        call set_1_rbox_gas_ptr( ipcg7_f_o_g,  kpcg7_f_o  )
        call set_1_rbox_gas_ptr( ipcg8_f_o_g,  kpcg8_f_o  )
        call set_1_rbox_gas_ptr( ipcg9_f_o_g,  kpcg9_f_o  )
        call set_1_rbox_gas_ptr( iopcg1_f_c_g, kopcg1_f_c )
        call set_1_rbox_gas_ptr( iopcg2_f_c_g, kopcg2_f_c )
        call set_1_rbox_gas_ptr( iopcg3_f_c_g, kopcg3_f_c )
        call set_1_rbox_gas_ptr( iopcg4_f_c_g, kopcg4_f_c )
        call set_1_rbox_gas_ptr( iopcg5_f_c_g, kopcg5_f_c )
        call set_1_rbox_gas_ptr( iopcg6_f_c_g, kopcg6_f_c )
        call set_1_rbox_gas_ptr( iopcg7_f_c_g, kopcg7_f_c )
        call set_1_rbox_gas_ptr( iopcg8_f_c_g, kopcg8_f_c )
        call set_1_rbox_gas_ptr( iopcg1_f_o_g, kopcg1_f_o )
        call set_1_rbox_gas_ptr( iopcg2_f_o_g, kopcg2_f_o )
        call set_1_rbox_gas_ptr( iopcg3_f_o_g, kopcg3_f_o )
        call set_1_rbox_gas_ptr( iopcg4_f_o_g, kopcg4_f_o )
        call set_1_rbox_gas_ptr( iopcg5_f_o_g, kopcg5_f_o )
        call set_1_rbox_gas_ptr( iopcg6_f_o_g, kopcg6_f_o )
        call set_1_rbox_gas_ptr( iopcg7_f_o_g, kopcg7_f_o )
        call set_1_rbox_gas_ptr( iopcg8_f_o_g, kopcg8_f_o )
        call set_1_rbox_gas_ptr( iant1_c_g,    kant1_c    )
        call set_1_rbox_gas_ptr( iant2_c_g,    kant2_c    )
        call set_1_rbox_gas_ptr( iant3_c_g,    kant3_c    )
        call set_1_rbox_gas_ptr( iant4_c_g,    kant4_c    )
        call set_1_rbox_gas_ptr( iant1_o_g,    kant1_o    )
        call set_1_rbox_gas_ptr( iant2_o_g,    kant2_o    )
        call set_1_rbox_gas_ptr( iant3_o_g,    kant3_o    )
        call set_1_rbox_gas_ptr( iant4_o_g,    kant4_o    )
        call set_1_rbox_gas_ptr( ibiog1_c_g,   kbiog1_c   )
        call set_1_rbox_gas_ptr( ibiog2_c_g,   kbiog2_c   )
        call set_1_rbox_gas_ptr( ibiog3_c_g,   kbiog3_c   )
        call set_1_rbox_gas_ptr( ibiog4_c_g,   kbiog4_c   )
        call set_1_rbox_gas_ptr( ibiog1_o_g,   kbiog1_o   )
        call set_1_rbox_gas_ptr( ibiog2_o_g,   kbiog2_o   )
        call set_1_rbox_gas_ptr( ibiog3_o_g,   kbiog3_o   )
        call set_1_rbox_gas_ptr( ibiog4_o_g,   kbiog4_o   )
        call set_1_rbox_gas_ptr( ismpa_g,      ksmpa      )
        call set_1_rbox_gas_ptr( ismpbb_g,     ksmpbb     )



        write(*,'(/a)') 'set_rbox_gas_ptrs'
        do igas = 1, ngas_aerchtot
           l = rbox_gas_ptr(igas)
           if ( (l < p1st) .or. (l > ntot_used) ) then
              rbox_gas_ptr(igas) = -1
              write(*,'(2a,2i7)')       'rbox_gas_ptr  ', gas_name(igas), igas, -999
           else
              write(*,'(2a,2i7,2x,a7)') 'rbox_gas_ptr  ', gas_name(igas), igas, l, name_rbox(l)
           end if
        end do

        return
        end subroutine set_rbox_gas_ptrs



        subroutine set_1_rbox_gas_ptr( igas, lgas )

        use module_data_mosaic_aero, only:  ngas_aerchtot

        use module_data_mosaic_asecthp, only:  rbox_gas_ptr

        use module_data_mosaic_main, only:  ntot_used
        
        use module_state_description, only:  p1st => param_first_scalar


        integer, intent(in)  :: igas, lgas

        if ( igas < 1    .or. igas > ngas_aerchtot) return
        if ( lgas < p1st .or. lgas > ntot_used    ) return

        rbox_gas_ptr(igas) = lgas

        return
        end subroutine set_1_rbox_gas_ptr
        


	subroutine mosaic2_set_otheraa( chem_opt, is_aerosol )



	use module_configure
	use module_state_description, only:  num_chem, p1st => param_first_scalar
	use module_scalar_tables, only:  chem_dname_table

	use module_data_mosaic_asecthp
	use module_data_mosaic_main, only:  ntot_used
	use module_data_mosaic_boxmod, only:   name_rbox
	use module_peg_util, only:  peg_error_fatal, peg_message


        integer, intent(in)  :: chem_opt
        logical, intent(inout) :: is_aerosol(num_chem)


	integer i, j, l, ll, n
	integer icomp, isize, itype, iphase

	character*200 msg
	character*20  tmpch20


	msg = ' '
	call peg_message( lunout, msg )
	msg = 'output from subr mosaic2_set_otheraa'
	call peg_message( lunout, trim(msg) )

	is_aerosol(1:num_chem) = .false.

	do iphase = 1, nphase_aer
	do itype  = 1, ntype_aer
	do isize  = 1, nsize_aer(itype)

	do icomp = 1, ncomp_aer(itype)
	    l = massptr_aer(icomp,isize,itype,iphase)
	    if (l >= p1st .and. l <= num_chem) then
		is_aerosol(l) = .true.
		name_rbox(l) = chem_dname_table(1,l)
	    end if
	end do 

	    l = numptr_aer(isize,itype,iphase)
	    if (l >= p1st .and. l <= num_chem) then
		is_aerosol(l) = .true.
		name_rbox(l) = chem_dname_table(1,l)
	    end if

	    if (iphase == ai_phase) then
	    l = waterptr_aer(isize,itype)
	    if (l >= p1st .and. l <= num_chem) then
		is_aerosol(l) = .true.
		name_rbox(l) = chem_dname_table(1,l)
	    end if
	    l = hyswptr_aer(isize,itype)
	    if (l >= p1st .and. l <= num_chem) then
		is_aerosol(l) = .true.
		name_rbox(l) = chem_dname_table(1,l)
	    end if
	    end if 

	end do 
	end do 
	end do 


       do l = 1, ntot_used
            tmpch20 = chem_dname_table(1,l)(1:20)
            do i = 1, 20
               j = ichar(tmpch20(i:i))
               if (j < 32 .or. j > 126) tmpch20(i:i) = '?'
            end do
	    msg = ' '
	    if (l > num_chem) then
		write(msg,'(a,i5,2(2x,a20))') 'name      ', l, name_rbox(l)(1:20)
	    else if (name_rbox(l) == chem_dname_table(1,l)) then
		write(msg,'(a,i5,2(2x,a20))') 'chem_dname', l, tmpch20
	    else
		write(msg,'(a,i5,2(2x,a20))') 'name_dname', l, name_rbox(l)(1:20), tmpch20
	    end if
	    call peg_message( lunout, trim(msg) )
        end do


	return
	end subroutine mosaic2_set_otheraa



	subroutine mosaic2_set_all_lnw_ptr( chem_opt, lmax, lunerr )

	use module_data_mosaic_asecthp, only:  &
	ai_phase, &
	maxd_asize, maxd_atype, maxd_aphase, &
	numptr_aer, hyswptr_aer, waterptr_aer, &
        lptr_so4_aer,       lptr_no3_aer,       lptr_cl_aer,        lptr_msa_aer,        &
        lptr_co3_aer,       lptr_nh4_aer,       lptr_na_aer,        lptr_ca_aer,         &
        lptr_oc_aer,        lptr_bc_aer,        lptr_oin_aer,                            &
        lptr_aro1_aer,      lptr_aro2_aer,      lptr_alk1_aer,      lptr_ole1_aer,       &
        lptr_api1_aer,      lptr_api2_aer,      lptr_lim1_aer,      lptr_lim2_aer,       &




                            lptr_pcg1_b_c_aer,  lptr_pcg2_b_c_aer,  lptr_pcg3_b_c_aer,   &
        lptr_pcg4_b_c_aer,  lptr_pcg5_b_c_aer,  lptr_pcg6_b_c_aer,  lptr_pcg7_b_c_aer,   &
        lptr_pcg8_b_c_aer,  lptr_pcg9_b_c_aer,  lptr_pcg1_b_o_aer,  lptr_pcg2_b_o_aer,   &
        lptr_pcg3_b_o_aer,  lptr_pcg4_b_o_aer,  lptr_pcg5_b_o_aer,  lptr_pcg6_b_o_aer,   &
        lptr_pcg7_b_o_aer,  lptr_pcg8_b_o_aer,  lptr_pcg9_b_o_aer,  lptr_opcg1_b_c_aer,  &
        lptr_opcg2_b_c_aer, lptr_opcg3_b_c_aer, lptr_opcg4_b_c_aer, lptr_opcg5_b_c_aer,  &
        lptr_opcg6_b_c_aer, lptr_opcg7_b_c_aer, lptr_opcg8_b_c_aer, lptr_opcg1_b_o_aer,  &
        lptr_opcg2_b_o_aer, lptr_opcg3_b_o_aer, lptr_opcg4_b_o_aer, lptr_opcg5_b_o_aer,  &
        lptr_opcg6_b_o_aer, lptr_opcg7_b_o_aer, lptr_opcg8_b_o_aer, lptr_pcg1_f_c_aer,   &
        lptr_pcg2_f_c_aer,  lptr_pcg3_f_c_aer,  lptr_pcg4_f_c_aer,  lptr_pcg5_f_c_aer,   &
        lptr_pcg6_f_c_aer,  lptr_pcg7_f_c_aer,  lptr_pcg8_f_c_aer,  lptr_pcg9_f_c_aer,   &
        lptr_pcg1_f_o_aer,  lptr_pcg2_f_o_aer,  lptr_pcg3_f_o_aer,  lptr_pcg4_f_o_aer,   &
        lptr_pcg5_f_o_aer,  lptr_pcg6_f_o_aer,  lptr_pcg7_f_o_aer,  lptr_pcg8_f_o_aer,   &
        lptr_pcg9_f_o_aer,  lptr_opcg1_f_c_aer, lptr_opcg2_f_c_aer, lptr_opcg3_f_c_aer,  &
        lptr_opcg4_f_c_aer, lptr_opcg5_f_c_aer, lptr_opcg6_f_c_aer, lptr_opcg7_f_c_aer,  &
        lptr_opcg8_f_c_aer, lptr_opcg1_f_o_aer, lptr_opcg2_f_o_aer, lptr_opcg3_f_o_aer,  &
        lptr_opcg4_f_o_aer, lptr_opcg5_f_o_aer, lptr_opcg6_f_o_aer, lptr_opcg7_f_o_aer,  &
        lptr_opcg8_f_o_aer, lptr_smpa_aer,      lptr_smpbb_aer,                          &


                                                                    lptr_ant1_c_aer,     &
        lptr_ant2_c_aer,    lptr_ant3_c_aer,    lptr_ant4_c_aer,    lptr_ant1_o_aer,     &
        lptr_ant2_o_aer,    lptr_ant3_o_aer,    lptr_ant4_o_aer,    lptr_biog1_c_aer,    &
        lptr_biog2_c_aer,   lptr_biog3_c_aer,   lptr_biog4_c_aer,   lptr_biog1_o_aer,    &
        lptr_biog2_o_aer,   lptr_biog3_o_aer,   lptr_biog4_o_aer

	integer, intent(in) :: chem_opt, lmax, lunerr

	integer :: lptr_tmp(maxd_asize,maxd_atype,maxd_aphase)



	call mosaic2_set_one_lnw_ptr( numptr_aer,   'num',   1, 0, lmax, lunerr )


	call mosaic2_set_one_lnw_ptr( lptr_tmp,     'water', 1, 1, lmax, lunerr )
	waterptr_aer(1:maxd_asize,1:maxd_atype) = lptr_tmp(1:maxd_asize,1:maxd_atype,ai_phase)
	call mosaic2_set_one_lnw_ptr( lptr_tmp,     'hysw',  1, 1, lmax, lunerr )
	hyswptr_aer(1:maxd_asize,1:maxd_atype)  = lptr_tmp(1:maxd_asize,1:maxd_atype,ai_phase)


	call mosaic2_set_one_lnw_ptr( lptr_so4_aer,       "so4",        1, 0, lmax, lunerr )
 

        call mosaic2_set_one_lnw_ptr( lptr_no3_aer,       "no3",        0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_cl_aer,        "cl",         0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_msa_aer,       "msa",        0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_co3_aer,       "co3",        0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_nh4_aer,       "nh4",        0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_na_aer,        "na",         0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ca_aer,        "ca",         0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_oc_aer,        "oc",         0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_bc_aer,        "bc",         0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_oin_aer,       "oin",        0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_aro1_aer,      "aro1",       0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_aro2_aer,      "aro2",       0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_alk1_aer,      "alk1",       0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ole1_aer,      "ole1",       0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_api1_aer,      "api1",       0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_api2_aer,      "api2",       0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_lim1_aer,      "lim1",       0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_lim2_aer,      "lim2",       0, 0, lmax, lunerr )










        call mosaic2_set_one_lnw_ptr( lptr_pcg1_b_c_aer,  "pcg1_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg2_b_c_aer,  "pcg2_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg3_b_c_aer,  "pcg3_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg4_b_c_aer,  "pcg4_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg5_b_c_aer,  "pcg5_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg6_b_c_aer,  "pcg6_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg7_b_c_aer,  "pcg7_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg8_b_c_aer,  "pcg8_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg9_b_c_aer,  "pcg9_b_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg1_b_o_aer,  "pcg1_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg2_b_o_aer,  "pcg2_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg3_b_o_aer,  "pcg3_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg4_b_o_aer,  "pcg4_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg5_b_o_aer,  "pcg5_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg6_b_o_aer,  "pcg6_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg7_b_o_aer,  "pcg7_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg8_b_o_aer,  "pcg8_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg9_b_o_aer,  "pcg9_b_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg1_b_c_aer, "opcg1_b_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg2_b_c_aer, "opcg2_b_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg3_b_c_aer, "opcg3_b_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg4_b_c_aer, "opcg4_b_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg5_b_c_aer, "opcg5_b_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg6_b_c_aer, "opcg6_b_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg7_b_c_aer, "opcg7_b_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg8_b_c_aer, "opcg8_b_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg1_b_o_aer, "opcg1_b_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg2_b_o_aer, "opcg2_b_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg3_b_o_aer, "opcg3_b_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg4_b_o_aer, "opcg4_b_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg5_b_o_aer, "opcg5_b_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg6_b_o_aer, "opcg6_b_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg7_b_o_aer, "opcg7_b_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg8_b_o_aer, "opcg8_b_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg1_f_c_aer,  "pcg1_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg2_f_c_aer,  "pcg2_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg3_f_c_aer,  "pcg3_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg4_f_c_aer,  "pcg4_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg5_f_c_aer,  "pcg5_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg6_f_c_aer,  "pcg6_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg7_f_c_aer,  "pcg7_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg8_f_c_aer,  "pcg8_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg9_f_c_aer,  "pcg9_f_c",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg1_f_o_aer,  "pcg1_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg2_f_o_aer,  "pcg2_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg3_f_o_aer,  "pcg3_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg4_f_o_aer,  "pcg4_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg5_f_o_aer,  "pcg5_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg6_f_o_aer,  "pcg6_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg7_f_o_aer,  "pcg7_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg8_f_o_aer,  "pcg8_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_pcg9_f_o_aer,  "pcg9_f_o",   0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg1_f_c_aer, "opcg1_f_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg2_f_c_aer, "opcg2_f_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg3_f_c_aer, "opcg3_f_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg4_f_c_aer, "opcg4_f_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg5_f_c_aer, "opcg5_f_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg6_f_c_aer, "opcg6_f_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg7_f_c_aer, "opcg7_f_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg8_f_c_aer, "opcg8_f_c",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg1_f_o_aer, "opcg1_f_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg2_f_o_aer, "opcg2_f_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg3_f_o_aer, "opcg3_f_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg4_f_o_aer, "opcg4_f_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg5_f_o_aer, "opcg5_f_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg6_f_o_aer, "opcg6_f_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg7_f_o_aer, "opcg7_f_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_opcg8_f_o_aer, "opcg8_f_o",  0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_smpa_aer,      "smpa",       0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_smpbb_aer,     "smpbb",      0, 0, lmax, lunerr )




        call mosaic2_set_one_lnw_ptr( lptr_ant1_c_aer,    "ant1_c",     0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ant2_c_aer,    "ant2_c",     0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ant3_c_aer,    "ant3_c",     0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ant4_c_aer,    "ant4_c",     0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ant1_o_aer,    "ant1_o",     0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ant2_o_aer,    "ant2_o",     0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ant3_o_aer,    "ant3_o",     0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_ant4_o_aer,    "ant4_o",     0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_biog1_c_aer,   "biog1_c",    0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_biog2_c_aer,   "biog2_c",    0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_biog3_c_aer,   "biog3_c",    0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_biog4_c_aer,   "biog4_c",    0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_biog1_o_aer,   "biog1_o",    0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_biog2_o_aer,   "biog2_o",    0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_biog3_o_aer,   "biog3_o",    0, 0, lmax, lunerr )
        call mosaic2_set_one_lnw_ptr( lptr_biog4_o_aer,   "biog4_o",    0, 0, lmax, lunerr )


	return
	end subroutine mosaic2_set_all_lnw_ptr



	subroutine mosaic2_set_one_lnw_ptr( lptr_aer, sname, required, ai_phase_only, lmax, lunerr )

	use module_data_mosaic_asecthp, only:  &
		ai_phase, &
		identical_comps_optaa, &
		maxd_asize, maxd_atype, maxd_aphase, &
		nsize_aer, ntype_aer, nphase_aer
	use module_peg_util, only:  peg_error_fatal

	integer, intent(out) :: lptr_aer(maxd_asize,maxd_atype,maxd_aphase)
	integer, intent(in)  :: required        
	integer, intent(in)  :: ai_phase_only   
	integer, intent(in)  :: lmax, lunerr

	character(len=*), intent(in)  :: sname

	integer :: ierr, iphase, isize, itype, l, ntot, ntotsv
	character(len=40)  :: chemname
	character(len=200) :: msg


	if (ai_phase /= 1) then
	    write(msg,*) 'mosaic2_set_one_lnw_ptr error 0 - bad ai_phase = ', ai_phase
	    call peg_error_fatal( lunerr, msg )
	end if


	lptr_aer(1:maxd_asize,1:maxd_atype,1:maxd_aphase) = 1

	do itype = 1, ntype_aer

	do iphase = 1, nphase_aer
	if ((ai_phase_only > 0) .and. (iphase /= 1)) cycle

	ntot = 0

	do isize = 1, nsize_aer(itype)
	    call mosaic2_form_chemname( sname, isize, itype, iphase, lunerr, chemname )
	    call mosaic2_find_chemname_in_table( chemname, l )

	    if (l <= 0) cycle

	    if (l > lmax) then
		write(msg,*) 'mosaic2_set_one_lnw_ptr error 1 - chemname, l, lmax = ', &
		    chemname, l, lmax
		call peg_error_fatal( lunerr, msg )
	    end if
	    lptr_aer(isize,itype,iphase) = l
	    ntot = ntot + 1
	end do

	if (itype == 1 .and. iphase == 1) ntotsv = ntot


	ierr = 0
	if (ntotsv == nsize_aer(itype)) ierr = 1
	if ( (required <= 0) .and. (ntotsv == 0) ) ierr = 1
	if ( ierr <= 0 ) then
	    write(msg,*) 'mosaic2_set_one_lnw_ptr error 2 - sname, itype, iphase, ntotsv, nsize = ', &
		sname, itype, iphase, ntotsv, nsize_aer(itype)
	    call peg_error_fatal( lunerr, msg )
	end if

	if (identical_comps_optaa > 0) then

	if ( ntotsv /= ntot ) then
	    write(msg,*) 'mosaic2_set_one_lnw_ptr error 3 - sname, itype, iphase, ntotsv, ntot = ', &
		sname, itype, iphase, ntotsv, ntot
	    call peg_error_fatal( lunerr, msg )
	end if
	end if

	end do 

	end do 

	return
	end subroutine mosaic2_set_one_lnw_ptr



        subroutine set_rbox_aer_ptrs




        use module_state_description, only:  p1st => param_first_scalar

        use module_data_mosaic_aero, only: &
            aer_name, naer, nbin_a, nbin_a_max, &
            iso4_a,     ino3_a,     icl_a,     inh4_a,     ico3_a,   &
            imsa_a,     ina_a,      ica_a,     ioc_a,      ibc_a,    &
            ioin_a,     iaro1_a,    iaro2_a,   ialk1_a,    iole1_a,  &
            iapi1_a,    iapi2_a,    ilim1_a,   ilim2_a,              &


            ipcg1_b_c_a,  ipcg2_b_c_a,  ipcg3_b_c_a,  ipcg4_b_c_a,  &
            ipcg5_b_c_a,  ipcg6_b_c_a,  ipcg7_b_c_a,  ipcg8_b_c_a,  &
            ipcg9_b_c_a,                                            &
            ipcg1_b_o_a,  ipcg2_b_o_a,  ipcg3_b_o_a,  ipcg4_b_o_a,  &
            ipcg5_b_o_a,  ipcg6_b_o_a,  ipcg7_b_o_a,  ipcg8_b_o_a,  &
            ipcg9_b_o_a,                                            &
            iopcg1_b_c_a, iopcg2_b_c_a, iopcg3_b_c_a, iopcg4_b_c_a, &
            iopcg5_b_c_a, iopcg6_b_c_a, iopcg7_b_c_a, iopcg8_b_c_a, &
            iopcg1_b_o_a, iopcg2_b_o_a, iopcg3_b_o_a, iopcg4_b_o_a, &
            iopcg5_b_o_a, iopcg6_b_o_a, iopcg7_b_o_a, iopcg8_b_o_a, &
            ipcg1_f_c_a,  ipcg2_f_c_a,  ipcg3_f_c_a,  ipcg4_f_c_a,  &
            ipcg5_f_c_a,  ipcg6_f_c_a,  ipcg7_f_c_a,  ipcg8_f_c_a,  &
            ipcg9_f_c_a,                                            &
            ipcg1_f_o_a,  ipcg2_f_o_a,  ipcg3_f_o_a,  ipcg4_f_o_a,  &
            ipcg5_f_o_a,  ipcg6_f_o_a,  ipcg7_f_o_a,  ipcg8_f_o_a,  &
            ipcg9_f_o_a,                                            &
            iopcg1_f_c_a, iopcg2_f_c_a, iopcg3_f_c_a, iopcg4_f_c_a, &
            iopcg5_f_c_a, iopcg6_f_c_a, iopcg7_f_c_a, iopcg8_f_c_a, &
            iopcg1_f_o_a, iopcg2_f_o_a, iopcg3_f_o_a, iopcg4_f_o_a, &
            iopcg5_f_o_a, iopcg6_f_o_a, iopcg7_f_o_a, iopcg8_f_o_a, &
            iant1_c_a,    iant2_c_a,    iant3_c_a,    iant4_c_a,    &
            iant1_o_a,    iant2_o_a,    iant3_o_a,    iant4_o_a,    &
            ibiog1_c_a,   ibiog2_c_a,   ibiog3_c_a,   ibiog4_c_a,   &
            ibiog1_o_a,   ibiog2_o_a,   ibiog3_o_a,   ibiog4_o_a,   &
            ismpa_a,      ismpbb_a





        use module_data_mosaic_asecthp, only: &
            ai_phase, &
            isize_of_ibin, itype_of_ibin, &
            rbox_aer_ptr, &
            numptr_aer, hyswptr_aer, waterptr_aer, &
            lptr_so4_aer,       lptr_no3_aer,       lptr_cl_aer,        lptr_msa_aer,        &
            lptr_co3_aer,       lptr_nh4_aer,       lptr_na_aer,        lptr_ca_aer,         &
            lptr_oc_aer,        lptr_bc_aer,        lptr_oin_aer,                            &
            lptr_aro1_aer,      lptr_aro2_aer,      lptr_alk1_aer,      lptr_ole1_aer,       &
            lptr_api1_aer,      lptr_api2_aer,      lptr_lim1_aer,      lptr_lim2_aer,       &




                                lptr_pcg1_b_c_aer,  lptr_pcg2_b_c_aer,  lptr_pcg3_b_c_aer,   &
            lptr_pcg4_b_c_aer,  lptr_pcg5_b_c_aer,  lptr_pcg6_b_c_aer,  lptr_pcg7_b_c_aer,   &
            lptr_pcg8_b_c_aer,  lptr_pcg9_b_c_aer,  lptr_pcg1_b_o_aer,  lptr_pcg2_b_o_aer,   &
            lptr_pcg3_b_o_aer,  lptr_pcg4_b_o_aer,  lptr_pcg5_b_o_aer,  lptr_pcg6_b_o_aer,   &
            lptr_pcg7_b_o_aer,  lptr_pcg8_b_o_aer,  lptr_pcg9_b_o_aer,  lptr_opcg1_b_c_aer,  &
            lptr_opcg2_b_c_aer, lptr_opcg3_b_c_aer, lptr_opcg4_b_c_aer, lptr_opcg5_b_c_aer,  &
            lptr_opcg6_b_c_aer, lptr_opcg7_b_c_aer, lptr_opcg8_b_c_aer, lptr_opcg1_b_o_aer,  &
            lptr_opcg2_b_o_aer, lptr_opcg3_b_o_aer, lptr_opcg4_b_o_aer, lptr_opcg5_b_o_aer,  &
            lptr_opcg6_b_o_aer, lptr_opcg7_b_o_aer, lptr_opcg8_b_o_aer, lptr_pcg1_f_c_aer,   &
            lptr_pcg2_f_c_aer,  lptr_pcg3_f_c_aer,  lptr_pcg4_f_c_aer,  lptr_pcg5_f_c_aer,   &
            lptr_pcg6_f_c_aer,  lptr_pcg7_f_c_aer,  lptr_pcg8_f_c_aer,  lptr_pcg9_f_c_aer,   &
            lptr_pcg1_f_o_aer,  lptr_pcg2_f_o_aer,  lptr_pcg3_f_o_aer,  lptr_pcg4_f_o_aer,   &
            lptr_pcg5_f_o_aer,  lptr_pcg6_f_o_aer,  lptr_pcg7_f_o_aer,  lptr_pcg8_f_o_aer,   &
            lptr_pcg9_f_o_aer,  lptr_opcg1_f_c_aer, lptr_opcg2_f_c_aer, lptr_opcg3_f_c_aer,  &
            lptr_opcg4_f_c_aer, lptr_opcg5_f_c_aer, lptr_opcg6_f_c_aer, lptr_opcg7_f_c_aer,  &
            lptr_opcg8_f_c_aer, lptr_opcg1_f_o_aer, lptr_opcg2_f_o_aer, lptr_opcg3_f_o_aer,  &
            lptr_opcg4_f_o_aer, lptr_opcg5_f_o_aer, lptr_opcg6_f_o_aer, lptr_opcg7_f_o_aer,  &
            lptr_opcg8_f_o_aer, lptr_smpa_aer,      lptr_smpbb_aer,                          &


                                                                        lptr_ant1_c_aer,     &
            lptr_ant2_c_aer,    lptr_ant3_c_aer,    lptr_ant4_c_aer,    lptr_ant1_o_aer,     &
            lptr_ant2_o_aer,    lptr_ant3_o_aer,    lptr_ant4_o_aer,    lptr_biog1_c_aer,    &
            lptr_biog2_c_aer,   lptr_biog3_c_aer,   lptr_biog4_c_aer,   lptr_biog1_o_aer,    &
            lptr_biog2_o_aer,   lptr_biog3_o_aer,   lptr_biog4_o_aer

        use module_data_mosaic_boxmod, only:  name_rbox

        use module_data_mosaic_main, only:  ntot_used


        integer :: iaer, ibin, iph, isz, ity, l

        character(len=16) :: tmp_name


        rbox_aer_ptr(-3:naer,1:nbin_a_max) = -1

        iph = ai_phase

        do ibin = 1, nbin_a

        isz = isize_of_ibin(ibin)
        ity = itype_of_ibin(ibin)

        rbox_aer_ptr( -1, ibin ) = numptr_aer(isz,ity,iph)
        rbox_aer_ptr( -2, ibin ) = waterptr_aer(isz,ity)
        rbox_aer_ptr( -3, ibin ) = hyswptr_aer(isz,ity)

        call set_1_rbox_aer_ptr( ibin, iso4_a,         lptr_so4_aer(      isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ino3_a,         lptr_no3_aer(      isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, icl_a,          lptr_cl_aer(       isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, inh4_a,         lptr_nh4_aer(      isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ico3_a,         lptr_co3_aer(      isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, imsa_a,         lptr_msa_aer(      isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ina_a,          lptr_na_aer(       isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ica_a,          lptr_ca_aer(       isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ioc_a,          lptr_oc_aer(       isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibc_a,          lptr_bc_aer(       isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ioin_a,         lptr_oin_aer(      isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iaro1_a,        lptr_aro1_aer(     isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iaro2_a,        lptr_aro2_aer(     isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ialk1_a,        lptr_alk1_aer(     isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iole1_a,        lptr_ole1_aer(     isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iapi1_a,        lptr_api1_aer(     isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iapi2_a,        lptr_api2_aer(     isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ilim1_a,        lptr_lim1_aer(     isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ilim2_a,        lptr_lim2_aer(     isz,ity,iph) )





        call set_1_rbox_aer_ptr( ibin, ipcg1_b_c_a,    lptr_pcg1_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg2_b_c_a,    lptr_pcg2_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg3_b_c_a,    lptr_pcg3_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg4_b_c_a,    lptr_pcg4_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg5_b_c_a,    lptr_pcg5_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg6_b_c_a,    lptr_pcg6_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg7_b_c_a,    lptr_pcg7_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg8_b_c_a,    lptr_pcg8_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg9_b_c_a,    lptr_pcg9_b_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg1_b_o_a,    lptr_pcg1_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg2_b_o_a,    lptr_pcg2_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg3_b_o_a,    lptr_pcg3_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg4_b_o_a,    lptr_pcg4_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg5_b_o_a,    lptr_pcg5_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg6_b_o_a,    lptr_pcg6_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg7_b_o_a,    lptr_pcg7_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg8_b_o_a,    lptr_pcg8_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg9_b_o_a,    lptr_pcg9_b_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg1_b_c_a,   lptr_opcg1_b_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg2_b_c_a,   lptr_opcg2_b_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg3_b_c_a,   lptr_opcg3_b_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg4_b_c_a,   lptr_opcg4_b_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg5_b_c_a,   lptr_opcg5_b_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg6_b_c_a,   lptr_opcg6_b_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg7_b_c_a,   lptr_opcg7_b_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg8_b_c_a,   lptr_opcg8_b_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg1_b_o_a,   lptr_opcg1_b_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg2_b_o_a,   lptr_opcg2_b_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg3_b_o_a,   lptr_opcg3_b_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg4_b_o_a,   lptr_opcg4_b_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg5_b_o_a,   lptr_opcg5_b_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg6_b_o_a,   lptr_opcg6_b_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg7_b_o_a,   lptr_opcg7_b_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg8_b_o_a,   lptr_opcg8_b_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg1_f_c_a,    lptr_pcg1_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg2_f_c_a,    lptr_pcg2_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg3_f_c_a,    lptr_pcg3_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg4_f_c_a,    lptr_pcg4_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg5_f_c_a,    lptr_pcg5_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg6_f_c_a,    lptr_pcg6_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg7_f_c_a,    lptr_pcg7_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg8_f_c_a,    lptr_pcg8_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg9_f_c_a,    lptr_pcg9_f_c_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg1_f_o_a,    lptr_pcg1_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg2_f_o_a,    lptr_pcg2_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg3_f_o_a,    lptr_pcg3_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg4_f_o_a,    lptr_pcg4_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg5_f_o_a,    lptr_pcg5_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg6_f_o_a,    lptr_pcg6_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg7_f_o_a,    lptr_pcg7_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg8_f_o_a,    lptr_pcg8_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ipcg9_f_o_a,    lptr_pcg9_f_o_aer( isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg1_f_c_a,   lptr_opcg1_f_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg2_f_c_a,   lptr_opcg2_f_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg3_f_c_a,   lptr_opcg3_f_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg4_f_c_a,   lptr_opcg4_f_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg5_f_c_a,   lptr_opcg5_f_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg6_f_c_a,   lptr_opcg6_f_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg7_f_c_a,   lptr_opcg7_f_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg8_f_c_a,   lptr_opcg8_f_c_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg1_f_o_a,   lptr_opcg1_f_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg2_f_o_a,   lptr_opcg2_f_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg3_f_o_a,   lptr_opcg3_f_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg4_f_o_a,   lptr_opcg4_f_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg5_f_o_a,   lptr_opcg5_f_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg6_f_o_a,   lptr_opcg6_f_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg7_f_o_a,   lptr_opcg7_f_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iopcg8_f_o_a,   lptr_opcg8_f_o_aer(isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iant1_c_a,      lptr_ant1_c_aer(   isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iant2_c_a,      lptr_ant2_c_aer(   isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iant3_c_a,      lptr_ant3_c_aer(   isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iant4_c_a,      lptr_ant4_c_aer(   isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iant1_o_a,      lptr_ant1_o_aer(   isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iant2_o_a,      lptr_ant2_o_aer(   isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iant3_o_a,      lptr_ant3_o_aer(   isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, iant4_o_a,      lptr_ant4_o_aer(   isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibiog1_c_a,     lptr_biog1_c_aer(  isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibiog2_c_a,     lptr_biog2_c_aer(  isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibiog3_c_a,     lptr_biog3_c_aer(  isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibiog4_c_a,     lptr_biog4_c_aer(  isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibiog1_o_a,     lptr_biog1_o_aer(  isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibiog2_o_a,     lptr_biog2_o_aer(  isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibiog3_o_a,     lptr_biog3_o_aer(  isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ibiog4_o_a,     lptr_biog4_o_aer(  isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ismpa_a,        lptr_smpa_aer(     isz,ity,iph) )
        call set_1_rbox_aer_ptr( ibin, ismpbb_a,       lptr_smpbb_aer(    isz,ity,iph) )










        end do 

        write(*,'(/a)') 'set_rbox_aer_ptrs'
        do ibin = 1, nbin_a
        write(*,'(a)')
        do iaer = -3, naer
           if      (iaer == -3) then
              tmp_name = 'hyswtr'
           else if (iaer == -2) then
              tmp_name = 'water'
           else if (iaer == -1) then
              tmp_name = 'num'
           else if (iaer ==  0) then
              cycle
           else
              tmp_name = aer_name(iaer)
           end if
           l = rbox_aer_ptr(iaer,ibin)
           if ( (l < p1st) .or.  (l > ntot_used) ) then
              rbox_aer_ptr(iaer,ibin) = -1
              if ( ibin <= 1 ) then
                 write(*,'(2a,i4,i5,i7,2x,a)') 'rbox_aer_ptr  ', tmp_name, ibin, iaer, -999
              else
                 if ( rbox_aer_ptr(iaer,ibin-1) > 0 ) &
                 write(*,'(2a,i4,i5,i7,2x,a)') 'rbox_aer_ptr  ', tmp_name, ibin, iaer, -999
              end if
           else
              write(*,'(2a,i4,i5,i7,2x,a)') 'rbox_aer_ptr  ', tmp_name, ibin, iaer, l, name_rbox(l)
           end if
        end do 
        end do 

        return
        end subroutine set_rbox_aer_ptrs



        subroutine set_1_rbox_aer_ptr( ibin, iaer, laer )

        use module_data_mosaic_aero, only:  naer, nbin_a_max

        use module_data_mosaic_asecthp, only:  rbox_aer_ptr

        use module_data_mosaic_main, only:  ntot_used
        
        use module_state_description, only:  p1st => param_first_scalar


        integer, intent(in)  :: ibin, iaer, laer

        if ( ibin < 1    .or. ibin > nbin_a_max) return
        if ( iaer < 1    .or. iaer > naer      ) return
        if ( laer < p1st .or. laer > ntot_used ) return

        rbox_aer_ptr(iaer,ibin) = laer

        return
        end subroutine set_1_rbox_aer_ptr
        


	subroutine mosaic2_set_massptr( chem_opt, lmax, lunerr )


	use module_data_mosaic_asecthp, only:  &
		identical_comps_optaa, &
		massptr_aer, &
                maxd_acomp, maxd_asize, maxd_atype, maxd_aphase, &
		ncomp_aer, nsize_aer, ntype_aer, nphase_aer, &
		sname_aer
	use module_peg_util, only:  peg_error_fatal

	integer, intent(in)  :: chem_opt, lmax, lunerr

	integer :: icomp, iphase, isize, itype, l, ntot, ntotsvaa, ntotsvbb
	character(len=40)  :: chemname
	character(len=200) :: msg



	massptr_aer(1:maxd_acomp,1:maxd_asize,1:maxd_atype,1:maxd_aphase) = -999888777
	ntotsvaa = -999888777 ; ntotsvbb = -999888777

	do itype = 1, ntype_aer

	do iphase = 1, nphase_aer

	do isize = 1, nsize_aer(itype)

	ntot = 0
	do icomp = 1, ncomp_aer(itype)
	    call mosaic2_form_chemname( sname_aer(icomp,itype), isize, itype, iphase, lunerr, chemname )
	    call mosaic2_find_chemname_in_table( chemname, l )
	    if (l <= 0) cycle

	    if (l > lmax) then
		write(msg,*) 'mosaic2_set_massptr error 1 - chemname, l, lmax = ', &
		    chemname, l, lmax
		call peg_error_fatal( lunerr, msg )
	    end if
	    ntot = ntot + 1
	    massptr_aer(icomp,isize,itype,iphase) = l
	end do 

	if ( itype == 1 .and. iphase == 1 .and. isize == 1 ) ntotsvaa = ntot
	if (                  iphase == 1 .and. isize == 1 ) ntotsvbb = ntot


	if (ntotsvbb /= ntot) then
	    write(msg,*) 'mosaic2_set_massptr error 2 - itype, iphase, isize, ntotsvbb, ntot = ', &
		itype, iphase, isize, ntotsvbb, ntot
	    call peg_error_fatal( lunerr, msg )
	end if

	if (identical_comps_optaa > 0) then

	if (ntotsvaa /= ntotsvbb) then
	    write(msg,*) 'mosaic2_set_massptr error 3 - itype, iphase, isize, ntotsvaa, ntotsvbb = ', &
		itype, iphase, isize, ntotsvaa, ntotsvbb
	    call peg_error_fatal( lunerr, msg )
	end if
	end if

	end do 

	end do 

	end do 

	return
	end subroutine mosaic2_set_massptr



	subroutine mosaic2_set_ncomp( chem_opt, lunerr )

	use module_data_mosaic_asecthp, only:  &
		dens_aer, dens_mastercomp_aer, &
		hygro_aer, hygro_mastercomp_aer, &
		identical_comps_optaa, is_tracer_mastercomp_aer, &
		mastercompptr_aer, maxd_acomp, maxd_atype, &
		mw_aer, mw_mastercomp_aer, &
		name_aer, name_mastercomp_aer, &
		ncomp_aer, ncomp_plustracer_aer, &
		ntot_mastercomp_aer, ntype_aer, &
		sname_aer, sname_mastercomp_aer, &
		type_chars_aer

	use module_peg_util, only:  peg_error_fatal

	integer, intent(in) :: chem_opt, lunerr

	integer :: icomp, icompmc, ipass, itype, l
	character(len=40)  :: chemname
	character(len=200) :: msg






	mastercompptr_aer(1:maxd_acomp,1:maxd_atype) = -999888777
	ncomp_aer(1:maxd_atype) = 0
	ncomp_plustracer_aer(1:maxd_atype) = 0

	name_aer( 1:maxd_acomp,1:maxd_atype) = 'empty'
	sname_aer(1:maxd_acomp,1:maxd_atype) = 'empty'

	dens_aer( 1:maxd_acomp,1:maxd_atype) = 1.0
	mw_aer(   1:maxd_acomp,1:maxd_atype) = 1.0
	hygro_aer(1:maxd_acomp,1:maxd_atype) = 1.0


	do itype = 1, ntype_aer

	icomp = 0

	do ipass = 1, 2

	do icompmc = 1, ntot_mastercomp_aer

	    if (ipass == 1) then
		if ( is_tracer_mastercomp_aer(icompmc) .eqv. .true. ) cycle
	    else
		if ( is_tracer_mastercomp_aer(icompmc) .eqv. .false. ) cycle
	    end if

	    chemname = trim(sname_mastercomp_aer(icompmc)) // '_a01' // type_chars_aer(itype)
	    call mosaic2_find_chemname_in_table( chemname, l )
	    if (l <= 0) cycle

	    icomp = icomp + 1
	    name_aer(icomp,itype) = name_mastercomp_aer(icompmc)
	    sname_aer(icomp,itype) = sname_mastercomp_aer(icompmc)
	    dens_aer(icomp,itype) = dens_mastercomp_aer(icompmc)
	    mw_aer(icomp,itype) = mw_mastercomp_aer(icompmc)
	    hygro_aer(icomp,itype) = hygro_mastercomp_aer(icompmc)

	end do 

	if (ipass == 1) then
	    ncomp_aer(itype) = icomp
	else
	    ncomp_plustracer_aer(itype) = icomp
	end if

	end do 

	end do 

	if (identical_comps_optaa > 0) then

	    do itype = 2, ntype_aer
	    do icomp = 1, max( ncomp_plustracer_aer(1), ncomp_plustracer_aer(itype) )
		if (mastercompptr_aer(icomp,1) /= mastercompptr_aer(icomp,itype)) then
		    write(msg,*) 'mosaic2_set_ncomp error 1 - itype, icomp, mcompptr1, mcompptr2ntot = ', &
			itype, icomp, mastercompptr_aer(icomp,1), mastercompptr_aer(icomp,itype)
		    call peg_error_fatal( lunerr, msg )
		end if
	    end do
	    end do
	end if

	return
	end subroutine mosaic2_set_ncomp



	subroutine mosaic2_set_bin_sizes( chem_opt, lunerr )

	use module_data_mosaic_kind, only:  r8
	use module_data_mosaic_asecthp, only:  &
		dcen_sect, dcut_sect, dhi_sect, dlo_sect, &
		maxd_asize, maxd_atype, &
		nsize_aer, ntype_aer, &
		sigmag_aer, &
		volumcen_sect, volumcut_sect, volumhi_sect, volumlo_sect
	use module_data_mosaic_constants, only: pi
	use module_peg_util, only:  peg_error_fatal

	integer, intent(in) :: chem_opt, lunerr

	integer :: itype, n, nhi
	real :: tmpa

	dlo_sect( 1:maxd_asize,1:maxd_atype) = 0.0
	dhi_sect( 1:maxd_asize,1:maxd_atype) = 0.0
	dcen_sect(1:maxd_asize,1:maxd_atype) = 0.0
	dcut_sect(0:maxd_asize,1:maxd_atype) = 0.0
	volumlo_sect( 1:maxd_asize,1:maxd_atype) = 0.0
	volumhi_sect( 1:maxd_asize,1:maxd_atype) = 0.0
	volumcen_sect(1:maxd_asize,1:maxd_atype) = 0.0
	volumcut_sect(0:maxd_asize,1:maxd_atype) = 0.0
	sigmag_aer(1:maxd_asize,1:maxd_atype) = 1.0

        do itype = 1, ntype_aer
	    nhi = nsize_aer(itype)
	    dlo_sect(1,itype) = 3.90625e-6_r8
	    dhi_sect(nhi,itype) = 10.0e-4_r8

	    tmpa = log( dhi_sect(nhi,itype)/dlo_sect(1,itype) ) / nhi
	    do n = 2, nhi
		dlo_sect(n,itype) = dlo_sect(1,itype) * exp( (n-1)*tmpa )
		dhi_sect(n-1,itype) = dlo_sect(n,itype)
	    end do
	    do n = 1, nhi
		dcen_sect(n,itype) = sqrt( dlo_sect(n,itype)*dhi_sect(n,itype) )
		volumlo_sect(n,itype) = (pi/6.) * (dlo_sect(n,itype)**3)
		volumhi_sect(n,itype) = (pi/6.) * (dhi_sect(n,itype)**3)
		volumcen_sect(n,itype) = (pi/6.) * (dcen_sect(n,itype)**3)
		sigmag_aer(n,itype) = (dhi_sect(n,itype)/dlo_sect(n,itype))**0.289
		dcut_sect(n,itype) = dhi_sect(n,itype)
		volumcut_sect(n,itype) = volumhi_sect(n,itype)
	    end do
	    dcut_sect(0,itype) = dlo_sect(1,itype)
	    volumcut_sect(0,itype) = volumlo_sect(1,itype)
	end do

	return
	end subroutine mosaic2_set_bin_sizes



	subroutine mosaic2_set_nsize( chem_opt, lunerr )

	use module_data_mosaic_asecthp, only:  &
		identical_sizes_optaa, &
		maxd_asize, nphase_aer, nsize_aer, ntype_aer
	use module_peg_util, only:  peg_error_fatal

	integer, intent(in) :: chem_opt, lunerr

	integer :: iphase, isize, itype, l, nsizetmp
	character(len=40)  :: chemname
	character(len=200) :: msg




	do itype = 1, ntype_aer

	do iphase = 1, nphase_aer

	nsizetmp = 0
	do isize = 1, maxd_asize+10
	    call mosaic2_form_chemname( 'so4', isize, itype, iphase, lunerr, chemname )
	    call mosaic2_find_chemname_in_table( chemname, l )
	    if (l <= 0) exit
	    nsizetmp = nsizetmp + 1
	end do

	if (iphase == 1) then
	    nsize_aer(itype) = nsizetmp
	else if (nsize_aer(itype) /= nsizetmp) then
	    write(msg,*) 'mosaic2_set_nsize error 1 - itype, iphase, nsize_aer, nsizetmp =', &
		itype, iphase, nsize_aer(itype), nsizetmp
	    call peg_error_fatal( lunerr, msg )
	end if

	end do 

	if ((nsize_aer(itype) < 1) .or. (nsize_aer(itype) > maxd_asize)) then
	    write(msg,*) 'mosaic2_set_nsize error 2 - itype, nsize_aer, maxd_asize =', itype, nsize_aer(itype), maxd_asize
	    call peg_error_fatal( lunerr, msg )
	end if

	if (identical_sizes_optaa > 0) then
	if (nsize_aer(itype) /= nsize_aer(1)) then
	    write(msg,*) 'mosaic2_set_nsize error 3 - itype, nsize, nsize =', itype, nsize_aer(itype), nsize_aer(1)
	    call peg_error_fatal( lunerr, msg )
	end if
	end if

	end do 

	return
	end subroutine mosaic2_set_nsize



	subroutine mosaic2_set_nphase( chem_opt, lunerr )

	use module_data_mosaic_asecthp, only:  &
		ai_phase, cw_phase, ci_phase, &
		identical_phases_optaa, &
		maxd_aphase, nphase_aer, ntype_aer, &
		phase_chars_aer, type_chars_aer
	use module_peg_util, only:  peg_error_fatal

	integer, intent(in) :: chem_opt, lunerr

	integer :: iphase, itype, l, nphasetmp
	character(len=40)  :: chemname
	character(len=200) :: msg
	character(len=3)   :: phase_chars_tmp







	nphase_aer = 0
	phase_chars_aer(1:maxd_aphase) = '_empty'

	do itype = 1, ntype_aer
	nphasetmp = 0

	do iphase = 1, maxd_aphase
	    if (iphase == 1) then
		phase_chars_tmp = '_a'
	    else if (iphase == 2) then
		phase_chars_tmp = '_cw'
	    else if (iphase == 3) then
		phase_chars_tmp = '_ci'
	    else 
		exit
	    end if
	    chemname = 'so4' // trim(phase_chars_tmp) // '01' // type_chars_aer(1)
	    call mosaic2_find_chemname_in_table( chemname, l )
	    if (l <= 0) cycle

	    nphasetmp = nphasetmp + 1


	    if (nphasetmp /= iphase) then
		write(msg,*) 'mosaic2_set_nphase error - nphasetmp, iphase =', nphasetmp, iphase
		call peg_error_fatal( lunerr, msg )
	    end if

	    if (iphase == 1) then
		ai_phase = nphasetmp
	    else if (iphase == 2) then
		cw_phase = nphasetmp

	    else if (iphase == 3) then
		ci_phase = nphasetmp
	    else
		cycle
	    end if
	    phase_chars_aer(iphase) = phase_chars_tmp
	end do 

	if (itype == 1) then
	    nphase_aer = nphasetmp
	    if ((nphase_aer < 1) .or. (nphase_aer > maxd_aphase)) then
		write(msg,*) 'mosaic2_set_nphase error 1 - nphase_aer, maxd_aphase =', nphase_aer, maxd_aphase
		call peg_error_fatal( lunerr, msg )
	    end if
	else
	    if (identical_phases_optaa > 0) then
	    if (nphase_aer /= nphasetmp) then
		write(msg,*) 'mosaic2_set_nphase error 2 - iphase, nphase_aer, nphasetmp =', iphase, nphase_aer, nphasetmp
		call peg_error_fatal( lunerr, msg )
	    end if
	    end if
	end if

	end do 

	return
	end subroutine mosaic2_set_nphase



	subroutine mosaic2_set_ntype( chem_opt, lunerr )

	use module_data_mosaic_asecthp, only:  maxd_atype, ntype_aer, type_chars_aer
	use module_peg_util, only:  peg_error_fatal

	integer, intent(in) :: chem_opt, lunerr

        integer :: itype, l
	character(len=40)  :: chemname
	character(len=200) :: msg




        ntype_aer = 0
        do itype = 1, 99
            if (itype == 1) then
                chemname = 'so4_a01'
            else
	        write( chemname, '(a,i2.2)' ) 'so4_a01_t', itype
            end if
	    call mosaic2_find_chemname_in_table( chemname, l )
	    if (l <= 0) exit
            ntype_aer = itype
        end do

	if (ntype_aer < 1 .or. ntype_aer > maxd_atype) then
	    write(msg,*) 'mosaic2_set_ntype error 1 - ' // &
                'ntype_aer, maxd_atype =', ntype_aer, maxd_atype
	    call peg_error_fatal( lunerr, msg )
	end if


	type_chars_aer(1:maxd_atype) = 'empty'
        do itype = 1, maxd_atype
	    type_chars_aer(itype) = ' '
            if (itype == 1) cycle
            write(type_chars_aer(itype), '(a,i2.2)' ) '_t', itype
        end do

	return
	end subroutine mosaic2_set_ntype



	subroutine mosaic2_set_3dbin_1dbin_ptrs( chem_opt, aer_extmix_opt, lunerr )



	use module_data_mosaic_kind, only:  r8
	use module_data_mosaic_asecthp, only:  &
		ibin_of_isize_itype, isize_of_ibin, itype_of_ibin, &
		itype_of_itype_md1md2, itype_md1_of_itype, itype_md2_of_itype, &
	        maxd_asize, maxd_atype, maxd_atype_md1, maxd_atype_md2, &
		nbin_a_max, nsize_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer, &
                xcut_atype_md1, xcut_atype_md2

	use module_peg_util, only:  peg_message, peg_error_fatal

	integer, intent(in)  :: chem_opt, aer_extmix_opt, lunerr

	integer :: iok, isize, itype, it1, it2, n
	character(len=230) :: msg



        iok = 0
        if      (ntype_aer == 1) then
            if ( aer_extmix_opt ==    0 .or. aer_extmix_opt ==    1 ) iok = 1
        else if (ntype_aer == 2) then
            if ( aer_extmix_opt == 2100 .or. aer_extmix_opt == 2101 .or. &
                 aer_extmix_opt == 2200 .or. aer_extmix_opt == 2201 ) iok = 1
        else if (ntype_aer == 3) then
            if ( aer_extmix_opt == 3100 .or. aer_extmix_opt == 3101 .or. &
                 aer_extmix_opt == 3200 .or. aer_extmix_opt == 3201 ) iok = 1
        else if (ntype_aer == 4) then
            if ( aer_extmix_opt == 4300 .or. aer_extmix_opt == 4301 ) iok = 1
        end if
        if (iok /= 1) then
	    write(msg,*) 'mosaic2_set_3dbin_1dbin_ptrs error - ' // &
                'ntype_aer, aer_extmix_opt =', ntype_aer, aer_extmix_opt
	    call peg_error_fatal( lunerr, msg )
        end if



        ntype_md1_aer = 1
        ntype_md2_aer = 1
        xcut_atype_md1(0:1) = (/ -0.10_r8, 1.10_r8 /)
        xcut_atype_md2(0:1) = (/ -0.10_r8, 5.10_r8 /)

        if ( ((aer_extmix_opt ==    0) .and. (ntype_aer == 1)) .or. &
             ((aer_extmix_opt ==    1) .and. (ntype_aer == 1)) ) then
            continue

        else if ( ((aer_extmix_opt == 2100) .and. (ntype_aer == 2)) .or. &
                  ((aer_extmix_opt == 2101) .and. (ntype_aer == 2)) ) then 
            ntype_md1_aer = 2
            xcut_atype_md1(0:2) = (/ -0.10_r8, 0.26_r8, 1.10_r8 /)

        else if ( ((aer_extmix_opt == 2200) .and. (ntype_aer == 2)) .or. &
                  ((aer_extmix_opt == 2201) .and. (ntype_aer == 2)) ) then 
            ntype_md2_aer = 2
            xcut_atype_md2(0:2) = (/ -0.10_r8, 0.20_r8, 5.10_r8 /)

        else if ( ((aer_extmix_opt == 3100) .and. (ntype_aer == 3)) .or. &
                  ((aer_extmix_opt == 3101) .and. (ntype_aer == 3)) ) then 
            ntype_md1_aer = 3
            xcut_atype_md1(0:3) = (/ -0.10_r8, 0.10_r8, 0.36_r8, 1.10_r8 /)

        else if ( ((aer_extmix_opt == 3200) .and. (ntype_aer == 3)) .or. &
                  ((aer_extmix_opt == 3201) .and. (ntype_aer == 3)) ) then 
            ntype_md2_aer = 3
            xcut_atype_md2(0:3) = (/ -0.10_r8, 0.06_r8, 0.20_r8, 5.10_r8 /)

        else if ( ((aer_extmix_opt == 4300) .and. (ntype_aer == 4)) .or. &
                  ((aer_extmix_opt == 4301) .and. (ntype_aer == 4)) ) then 
            ntype_md1_aer = 2
            ntype_md2_aer = 2
            xcut_atype_md1(0:2) = (/ -0.10_r8, 0.30_r8, 1.10_r8 /)
            xcut_atype_md2(0:2) = (/ -0.10_r8, 0.10_r8, 5.10_r8 /)

        end if

        if ( ntype_md1_aer < 1 .or. ntype_md1_aer > maxd_atype_md1 ) then
	    write(msg,*) 'mosaic2_set_3dbin_1dbin_ptrs error - ' // &
                'ntype_md1_aer, maxd_atype_md1 =', &
                ntype_md1_aer, maxd_atype_md1
	    call peg_error_fatal( lunerr, msg )
        end if
        if ( ntype_md2_aer < 1 .or. ntype_md2_aer > maxd_atype_md2 ) then
	    write(msg,*) 'mosaic2_set_3dbin_1dbin_ptrs error - ' // &
                'ntype_md2_aer, maxd_atype_md2 =', &
                ntype_md2_aer, maxd_atype_md2
	    call peg_error_fatal( lunerr, msg )
        end if
        if ( ntype_md1_aer*ntype_md2_aer /= ntype_aer ) then
	    write(msg,*) 'mosaic2_set_3dbin_1dbin_ptrs error - ' // &
            'ntype_md1_aer, ntype_md2_aer, ntype_aer =', ntype_md1_aer, ntype_md2_aer, ntype_aer
	    call peg_error_fatal( lunerr, msg )
        end if


        if ((aer_extmix_opt >= 2000) .and. (mod(aer_extmix_opt,100) == 0)) then
            xcut_atype_md1(0:4) = (/ -0.10_r8, 1.10_r8, 2.10_r8, 3.10_r8, 4.10_r8 /)
            xcut_atype_md2(0:4) = (/ -0.10_r8, 5.10_r8, 6.10_r8, 7.10_r8, 8.10_r8 /)
        end if


        do it1 = ntype_md1_aer+1, maxd_atype_md1
            xcut_atype_md1(it1) = xcut_atype_md1(it1-1) + 1.0_r8
        end do
        do it2 = ntype_md2_aer+1, maxd_atype_md2
            xcut_atype_md2(it2) = xcut_atype_md2(it2-1) + 1.0_r8
        end do


	ibin_of_isize_itype(1:maxd_asize,1:maxd_atype) = -999888777
	isize_of_ibin(1:nbin_a_max) = -999888777
	itype_of_ibin(1:nbin_a_max) = -999888777
	n = 0

	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    n = n + 1
	    ibin_of_isize_itype(isize,itype) = n
	    isize_of_ibin(n) = isize
	    itype_of_ibin(n) = itype
	end do 
	end do 

	itype_of_itype_md1md2(1:maxd_atype_md1,1:maxd_atype_md2) = -999888777
	itype_md1_of_itype(1:maxd_atype) = -999888777
	itype_md2_of_itype(1:maxd_atype) = -999888777
        itype = 0

        do it2 = 1, ntype_md2_aer
        do it1 = 1, ntype_md1_aer
            itype = itype + 1
	    itype_of_itype_md1md2(it1,it2) = itype
	    itype_md1_of_itype(itype) = it1
	    itype_md2_of_itype(itype) = it2
        end do
        end do

	write(msg,'(a,3i10)') 'ntype_aer, ntype_md1_aer, ntype_md2_aer', &
                               ntype_aer, ntype_md1_aer, ntype_md2_aer
	call peg_message( lunerr, msg )
	
	write(msg,'(a,1p,21e10.2)') 'xcut_atype_md1', xcut_atype_md1(0:min(maxd_atype_md1,20))
	call peg_message( lunerr, msg )
	write(msg,'(a,1p,21e10.2)') 'xcut_atype_md2', xcut_atype_md2(0:min(maxd_atype_md2,20))
	call peg_message( lunerr, msg )

	return
	end subroutine mosaic2_set_3dbin_1dbin_ptrs



	subroutine mosaic2_form_chemname( sname, isize, itype, iphase, lunerr, chemname )

	use module_data_mosaic_asecthp, only:  phase_chars_aer, type_chars_aer
	use module_peg_util, only:  peg_error_fatal

	integer, intent(in) :: isize, itype, iphase, lunerr
	character(len=*), intent(in) :: sname
	character(len=*), intent(inout) :: chemname

	integer :: ierr
	character(len=3)   :: size_chars
	character(len=200) :: msg

	size_chars = ' '
	ierr = 0

	if ( isize < 1 .or. isize > 999 ) then
	    ierr = 1
	else if (isize <= 99) then
	    write(size_chars(1:2),'(i2.2)') isize
	else
	    write(size_chars(1:3),'(i3.3)') isize
	end if

	if (ierr > 0) then
	    write(msg,*) 'mosaic2_form_chemname error 1 - sname, isize = ', trim(sname), isize
	    call peg_error_fatal( lunerr, msg )
	end if

	chemname = trim(sname) // trim(phase_chars_aer(iphase)) // trim(size_chars) // type_chars_aer(itype)

	return
	end subroutine mosaic2_form_chemname



	subroutine mosaic2_find_chemname_in_table( chemname, chemindx )

	use module_state_description, only:  num_chem,  param_first_scalar
	use module_scalar_tables, only:  chem_dname_table

	integer, intent(out) :: chemindx
	character(len=*), intent(in) :: chemname

	integer :: l

	chemindx = 0
	do l = param_first_scalar, num_chem
	    if (chem_dname_table(1,l) == chemname) then
		chemindx = l
		return
	    end if
	end do

	return
	end subroutine mosaic2_find_chemname_in_table



        subroutine mosaic2_set_mastercomp

        use module_data_mosaic_asecthp

        integer :: l




        name_mastercomp_aer( 1:maxd_acomp ) = 'empty'
        sname_mastercomp_aer(1:maxd_acomp ) = 'empty'
        dens_mastercomp_aer( 1:maxd_acomp ) = 1.0
        mw_mastercomp_aer(   1:maxd_acomp ) = 1.0
        hygro_mastercomp_aer(1:maxd_acomp ) = 0.0

        is_tracer_mastercomp_aer(1:maxd_acomp) = .false.


	
	ntot_mastercomp_aer = 109 



	l = 1
	mastercompindx_so4_aer = l
	name_mastercomp_aer( l ) = 'sulfate'
	sname_mastercomp_aer( l ) = 'so4'
	dens_mastercomp_aer( l ) =  dens_so4_aer
	mw_mastercomp_aer(   l ) =    mw_so4_aer
	hygro_mastercomp_aer(l ) = hygro_so4_aer

	l = 2
	mastercompindx_no3_aer = l
	name_mastercomp_aer( l ) = 'nitrate'
	sname_mastercomp_aer( l ) = 'no3'
	dens_mastercomp_aer( l ) =  dens_no3_aer
	mw_mastercomp_aer(   l ) =    mw_no3_aer
	hygro_mastercomp_aer(l ) = hygro_no3_aer

	l = 3
	mastercompindx_cl_aer = l
	name_mastercomp_aer( l ) = 'chloride'
	sname_mastercomp_aer( l ) = 'cl'
	dens_mastercomp_aer( l ) =  dens_cl_aer
	mw_mastercomp_aer(   l ) =    mw_cl_aer
	hygro_mastercomp_aer(l ) = hygro_cl_aer

        l = 4
        mastercompindx_msa_aer = l
        name_mastercomp_aer( l ) = 'msa'
        dens_mastercomp_aer( l ) =  dens_msa_aer
        mw_mastercomp_aer(   l ) =    mw_msa_aer
        hygro_mastercomp_aer(l ) = hygro_msa_aer

	l = 5
	mastercompindx_co3_aer = l
	name_mastercomp_aer( l ) = 'carbonate'
	sname_mastercomp_aer( l ) = 'co3'
	dens_mastercomp_aer( l ) =  dens_co3_aer
	mw_mastercomp_aer(   l ) =    mw_co3_aer
	hygro_mastercomp_aer(l ) = hygro_co3_aer

	l = 6
	mastercompindx_nh4_aer = l
	name_mastercomp_aer( l ) = 'ammonium'
	sname_mastercomp_aer( l ) = 'nh4'
	dens_mastercomp_aer( l ) =  dens_nh4_aer
	mw_mastercomp_aer(   l ) =    mw_nh4_aer
	hygro_mastercomp_aer(l ) = hygro_nh4_aer

	l = 7
	mastercompindx_na_aer = l
	name_mastercomp_aer( l ) = 'sodium'
	sname_mastercomp_aer( l ) = 'na'
	dens_mastercomp_aer( l ) =  dens_na_aer
	mw_mastercomp_aer(   l ) =    mw_na_aer
	hygro_mastercomp_aer(l ) = hygro_na_aer

	l = 8
	mastercompindx_ca_aer = l
	name_mastercomp_aer( l ) = 'calcium'
	sname_mastercomp_aer( l ) = 'ca'
	dens_mastercomp_aer( l ) =  dens_ca_aer
	mw_mastercomp_aer(   l ) =    mw_ca_aer
	hygro_mastercomp_aer(l ) = hygro_ca_aer

	l = 9
	mastercompindx_oin_aer = l
	name_mastercomp_aer( l ) = 'otherinorg'
	sname_mastercomp_aer( l ) = 'oin'
	dens_mastercomp_aer( l ) =  dens_oin_aer
	mw_mastercomp_aer(   l ) =    mw_oin_aer
	hygro_mastercomp_aer(l ) = hygro_oin_aer








	l = 11
	mastercompindx_oc_aer = l
	name_mastercomp_aer( l ) = 'organic-c'
	sname_mastercomp_aer( l ) = 'oc'
	dens_mastercomp_aer( l ) =  dens_oc_aer
	mw_mastercomp_aer(   l ) =    mw_oc_aer
	hygro_mastercomp_aer(l ) = hygro_oc_aer

	l = 12
	mastercompindx_bc_aer = l
	name_mastercomp_aer( l ) = 'black-c'
	sname_mastercomp_aer( l ) = 'bc'
	dens_mastercomp_aer( l ) =  dens_bc_aer
	mw_mastercomp_aer(   l ) =    mw_bc_aer
	hygro_mastercomp_aer(l ) = hygro_bc_aer

        l = 13
        mastercompindx_pcg1_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg1_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg1_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg1_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg1_b_c_aer

        l = 14
        mastercompindx_pcg2_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg2_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg2_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg2_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg2_b_c_aer

        l = 15
        mastercompindx_pcg3_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg3_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg3_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg3_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg3_b_c_aer

        l = 16
        mastercompindx_pcg4_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg4_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg4_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg4_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg4_b_c_aer

        l = 17
        mastercompindx_pcg5_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg5_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg5_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg5_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg5_b_c_aer

        l = 18
        mastercompindx_pcg6_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg6_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg6_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg6_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg6_b_c_aer

        l = 19
        mastercompindx_pcg7_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg7_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg7_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg7_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg7_b_c_aer

        l = 20
        mastercompindx_pcg8_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg8_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg8_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg8_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg8_b_c_aer

        l = 21
        mastercompindx_pcg9_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg9_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg9_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg9_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg9_b_c_aer

        l = 22
        mastercompindx_pcg1_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg1_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg1_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg1_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg1_b_o_aer

        l = 23
        mastercompindx_pcg2_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg2_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg2_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg2_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg2_b_o_aer

        l = 24
        mastercompindx_pcg3_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg3_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg3_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg3_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg3_b_o_aer

        l = 25
        mastercompindx_pcg4_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg4_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg4_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg4_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg4_b_o_aer

        l = 26
        mastercompindx_pcg5_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg5_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg5_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg5_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg5_b_o_aer

        l = 27
        mastercompindx_pcg6_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg6_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg6_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg6_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg6_b_o_aer

        l = 28
        mastercompindx_pcg7_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg7_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg7_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg7_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg7_b_o_aer

        l = 29
        mastercompindx_pcg8_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg8_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg8_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg8_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg8_b_o_aer

        l = 30
        mastercompindx_pcg9_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg9_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg9_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg9_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg9_b_o_aer

        l = 31
        mastercompindx_opcg1_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg1_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg1_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg1_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg1_b_c_aer

        l = 32
        mastercompindx_opcg2_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg2_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg2_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg2_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg2_b_c_aer

        l = 33
        mastercompindx_opcg3_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg3_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg3_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg3_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg3_b_c_aer

        l = 34
        mastercompindx_opcg4_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg4_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg4_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg4_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg4_b_c_aer

        l = 35
        mastercompindx_opcg5_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg5_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg5_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg5_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg5_b_c_aer

        l = 36
        mastercompindx_opcg6_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg6_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg6_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg6_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg6_b_c_aer

        l = 37
        mastercompindx_opcg7_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg7_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg7_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg7_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg7_b_c_aer

        l = 38
        mastercompindx_opcg8_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg8_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg8_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg8_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg8_b_c_aer

        l = 39
        mastercompindx_opcg1_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg1_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg1_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg1_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg1_b_o_aer

        l = 40
        mastercompindx_opcg2_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg2_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg2_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg2_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg2_b_o_aer

        l = 41
        mastercompindx_opcg3_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg3_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg3_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg3_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg3_b_o_aer

        l = 42
        mastercompindx_opcg4_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg4_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg4_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg4_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg4_b_o_aer

        l = 43
        mastercompindx_opcg5_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg5_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg5_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg5_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg5_b_o_aer

        l = 44
        mastercompindx_opcg6_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg6_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg6_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg6_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg6_b_o_aer

        l = 45
        mastercompindx_opcg7_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg7_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg7_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg7_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg7_b_o_aer

        l = 46
        mastercompindx_opcg8_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg8_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg8_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg8_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg8_b_o_aer

        l = 47
        mastercompindx_pcg1_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg1_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg1_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg1_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg1_f_c_aer

        l = 48
        mastercompindx_pcg2_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg2_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg2_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg2_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg2_f_c_aer

        l = 49
        mastercompindx_pcg3_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg3_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg3_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg3_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg3_f_c_aer

        l = 50
        mastercompindx_pcg4_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg4_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg4_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg4_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg4_f_c_aer

        l = 51
        mastercompindx_pcg5_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg5_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg5_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg5_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg5_f_c_aer

        l = 52
        mastercompindx_pcg6_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg6_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg6_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg6_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg6_f_c_aer

        l = 53
        mastercompindx_pcg7_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg7_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg7_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg7_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg7_f_c_aer

        l = 54
        mastercompindx_pcg8_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg8_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg8_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg8_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg8_f_c_aer

        l = 55
        mastercompindx_pcg9_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg9_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg9_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg9_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg9_f_c_aer

        l = 56
        mastercompindx_pcg1_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg1_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg1_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg1_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg1_f_o_aer

        l = 57
        mastercompindx_pcg2_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg2_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg2_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg2_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg2_f_o_aer

        l = 58
        mastercompindx_pcg3_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg3_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg3_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg3_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg3_f_o_aer

        l = 59
        mastercompindx_pcg4_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg4_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg4_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg4_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg4_f_o_aer

        l = 60
        mastercompindx_pcg5_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg5_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg5_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg5_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg5_f_o_aer

        l = 61
        mastercompindx_pcg6_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg6_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg6_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg6_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg6_f_o_aer

        l = 62
        mastercompindx_pcg7_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg7_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg7_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg7_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg7_f_o_aer

        l = 63
        mastercompindx_pcg8_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg8_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg8_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg8_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg8_f_o_aer

        l = 64
        mastercompindx_pcg9_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg9_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg9_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg9_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg9_f_o_aer

        l = 65
        mastercompindx_opcg1_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg1_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg1_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg1_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg1_f_c_aer

        l = 66
        mastercompindx_opcg2_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg2_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg2_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg2_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg2_f_c_aer

        l = 67
        mastercompindx_opcg3_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg3_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg3_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg3_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg3_f_c_aer

        l = 68
        mastercompindx_opcg4_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg4_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg4_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg4_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg4_f_c_aer

        l = 69
        mastercompindx_opcg5_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg5_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg5_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg5_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg5_f_c_aer

        l = 70
        mastercompindx_opcg6_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg6_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg6_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg6_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg6_f_c_aer

        l = 71
        mastercompindx_opcg7_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg7_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg7_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg7_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg7_f_c_aer

        l = 72
        mastercompindx_opcg8_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg8_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg8_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg8_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg8_f_c_aer

        l = 73
        mastercompindx_opcg1_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg1_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg1_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg1_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg1_f_o_aer

        l = 74
        mastercompindx_opcg2_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg2_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg2_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg2_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg2_f_o_aer

        l = 75
        mastercompindx_opcg3_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg3_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg3_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg3_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg3_f_o_aer

        l = 76
        mastercompindx_opcg4_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg4_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg4_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg4_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg4_f_o_aer

        l = 77
        mastercompindx_opcg5_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg5_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg5_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg5_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg5_f_o_aer

        l = 78
        mastercompindx_opcg6_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg6_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg6_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg6_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg6_f_o_aer

        l = 79
        mastercompindx_opcg7_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg7_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg7_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg7_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg7_f_o_aer

        l = 80
        mastercompindx_opcg8_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg8_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg8_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg8_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg8_f_o_aer

        l = 81
        mastercompindx_ant1_c_aer = l
        name_mastercomp_aer( l ) = 'ant1_c'
        dens_mastercomp_aer( l ) =  dens_ant1_c_aer
        mw_mastercomp_aer(   l ) =    mw_ant1_c_aer
        hygro_mastercomp_aer(l ) = hygro_ant1_c_aer

        l = 82
        mastercompindx_ant2_c_aer = l
        name_mastercomp_aer( l ) = 'ant2_c'
        dens_mastercomp_aer( l ) =  dens_ant2_c_aer
        mw_mastercomp_aer(   l ) =    mw_ant2_c_aer
        hygro_mastercomp_aer(l ) = hygro_ant2_c_aer

        l = 83
        mastercompindx_ant3_c_aer = l
        name_mastercomp_aer( l ) = 'ant3_c'
        dens_mastercomp_aer( l ) =  dens_ant3_c_aer
        mw_mastercomp_aer(   l ) =    mw_ant3_c_aer
        hygro_mastercomp_aer(l ) = hygro_ant3_c_aer

        l = 84
        mastercompindx_ant4_c_aer = l
        name_mastercomp_aer( l ) = 'ant4_c'
        dens_mastercomp_aer( l ) =  dens_ant4_c_aer
        mw_mastercomp_aer(   l ) =    mw_ant4_c_aer
        hygro_mastercomp_aer(l ) = hygro_ant4_c_aer

        l = 85
        mastercompindx_ant1_o_aer = l
        name_mastercomp_aer( l ) = 'ant1_o'
        dens_mastercomp_aer( l ) =  dens_ant1_o_aer
        mw_mastercomp_aer(   l ) =    mw_ant1_o_aer
        hygro_mastercomp_aer(l ) = hygro_ant1_o_aer

        l = 86
        mastercompindx_ant2_o_aer = l
        name_mastercomp_aer( l ) = 'ant2_o'
        dens_mastercomp_aer( l ) =  dens_ant2_o_aer
        mw_mastercomp_aer(   l ) =    mw_ant2_o_aer
        hygro_mastercomp_aer(l ) = hygro_ant2_o_aer

        l = 87
        mastercompindx_ant3_o_aer = l
        name_mastercomp_aer( l ) = 'ant3_o'
        dens_mastercomp_aer( l ) =  dens_ant3_o_aer
        mw_mastercomp_aer(   l ) =    mw_ant3_o_aer
        hygro_mastercomp_aer(l ) = hygro_ant3_o_aer

        l = 88
        mastercompindx_ant4_o_aer = l
        name_mastercomp_aer( l ) = 'ant4_o'
        dens_mastercomp_aer( l ) =  dens_ant4_o_aer
        mw_mastercomp_aer(   l ) =    mw_ant4_o_aer
        hygro_mastercomp_aer(l ) = hygro_ant4_o_aer


        l = 89
        mastercompindx_biog1_c_aer = l
        name_mastercomp_aer( l ) = 'biog1_c'
        dens_mastercomp_aer( l ) =  dens_biog1_c_aer
        mw_mastercomp_aer(   l ) =    mw_biog1_c_aer
        hygro_mastercomp_aer(l ) = hygro_biog1_c_aer

        l = 90
        mastercompindx_biog2_c_aer = l
        name_mastercomp_aer( l ) = 'biog2_c'
        dens_mastercomp_aer( l ) =  dens_biog2_c_aer
        mw_mastercomp_aer(   l ) =    mw_biog2_c_aer
        hygro_mastercomp_aer(l ) = hygro_biog2_c_aer

        l = 91
        mastercompindx_biog3_c_aer = l
        name_mastercomp_aer( l ) = 'biog3_c'
        dens_mastercomp_aer( l ) =  dens_biog3_c_aer
        mw_mastercomp_aer(   l ) =    mw_biog3_c_aer
        hygro_mastercomp_aer(l ) = hygro_biog3_c_aer

        l = 92
        mastercompindx_biog4_c_aer = l
        name_mastercomp_aer( l ) = 'biog4_c'
        dens_mastercomp_aer( l ) =  dens_biog4_c_aer
        mw_mastercomp_aer(   l ) =    mw_biog4_c_aer
        hygro_mastercomp_aer(l ) = hygro_biog4_c_aer

        l = 93
        mastercompindx_biog1_o_aer = l
        name_mastercomp_aer( l ) = 'biog1_o'
        dens_mastercomp_aer( l ) =  dens_biog1_o_aer
        mw_mastercomp_aer(   l ) =    mw_biog1_o_aer
        hygro_mastercomp_aer(l ) = hygro_biog1_o_aer

        l = 94
        mastercompindx_biog2_o_aer = l
        name_mastercomp_aer( l ) = 'biog2_o'
        dens_mastercomp_aer( l ) =  dens_biog2_o_aer
        mw_mastercomp_aer(   l ) =    mw_biog2_o_aer
        hygro_mastercomp_aer(l ) = hygro_biog2_o_aer

        l = 95
        mastercompindx_biog3_o_aer = l
        name_mastercomp_aer( l ) = 'biog3_o'
        dens_mastercomp_aer( l ) =  dens_biog3_o_aer
        mw_mastercomp_aer(   l ) =    mw_biog3_o_aer
        hygro_mastercomp_aer(l ) = hygro_biog3_o_aer

        l = 96
        mastercompindx_biog4_o_aer = l
        name_mastercomp_aer( l ) = 'biog4_o'
        dens_mastercomp_aer( l ) =  dens_biog4_o_aer
        mw_mastercomp_aer(   l ) =    mw_biog4_o_aer
        hygro_mastercomp_aer(l ) = hygro_biog4_o_aer

        l = 95
        mastercompindx_smpa_aer = l
        name_mastercomp_aer( l ) = 'smpa'
        dens_mastercomp_aer( l ) =  dens_smpa_aer
        mw_mastercomp_aer(   l ) =    mw_smpa_aer
        hygro_mastercomp_aer(l ) = hygro_smpa_aer

        l = 96
        mastercompindx_smpbb_aer = l
        name_mastercomp_aer( l ) = 'smpbb'
        dens_mastercomp_aer( l ) =  dens_smpbb_aer
        mw_mastercomp_aer(   l ) =    mw_smpbb_aer
        hygro_mastercomp_aer(l ) = hygro_smpbb_aer































































































        do l = 1, ntot_mastercomp_aer
            if ( sname_mastercomp_aer( l ) == 'empty' ) &
                 sname_mastercomp_aer( l ) = name_mastercomp_aer( l )
        end do

	return
	end subroutine mosaic2_set_mastercomp




	end module module_mosaic2_driver

