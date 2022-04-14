








	module module_sorgam_vbs_cloudchem



	integer, parameter :: l_so4_aqyy = 1
	integer, parameter :: l_no3_aqyy = 2
	integer, parameter :: l_cl_aqyy  = 3
	integer, parameter :: l_nh4_aqyy = 4
	integer, parameter :: l_na_aqyy  = 5
	integer, parameter :: l_oin_aqyy = 6
	integer, parameter :: l_bc_aqyy  = 7
	integer, parameter :: l_oc_aqyy  = 8

	integer, parameter :: nyyy = 8


	real, parameter :: smallvolaa = 0.5e-18



	contains




	subroutine sorgam_vbs_cloudchem_driver(   &
	    id, ktau, ktauc, dtstepc, config_flags,   &
	    p_phy, t_phy, rho_phy, alt,   &
	    cldfra, ph_no2,   &
	    moist, chem,   &
	    gas_aqfrac, numgas_aqfrac,   &
	    ids,ide, jds,jde, kds,kde,   &
	    ims,ime, jms,jme, kms,kme,   &
	    its,ite, jts,jte, kts,kte )

	use module_state_description, only:   &
		num_moist, num_chem, p_qc

	use module_configure, only:  grid_config_rec_type

	use module_data_sorgam_vbs, only:  cw_phase, nphase_aer


	implicit none


	integer, intent(in) ::   &
		id, ktau, ktauc,   &
		numgas_aqfrac,   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte













	type(grid_config_rec_type), intent(in) :: config_flags


	real, intent(in) ::   &
	    dtstepc


        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                p_phy, t_phy, rho_phy, alt, cldfra, ph_no2







        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
                moist



        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
                chem



        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme, numgas_aqfrac ) :: &
                gas_aqfrac




	integer :: it, jt, kt, l
	integer :: icase
	integer :: igaschem_onoff, iphotol_onoff, iradical_onoff

	real :: rbox(num_chem)
	real :: gas_aqfrac_box(numgas_aqfrac)
	real :: ph_aq_box
	real, parameter :: qcldwtr_cutoff = 1.0e-6
	real :: qcldwtr



        if ((cw_phase .le. 0) .or. (cw_phase .gt. nphase_aer)) then
            print *, '*** sorgam_vbs_cloudchem_driver - cw_phase not active'
            return
        end if

	print 93010, 'entering sorgam_vbs_cloudchem_driver - ktau =', ktau

	icase = 0


	iphotol_onoff = 0
	if (config_flags%phot_opt .gt. 0) iphotol_onoff = 1

	igaschem_onoff = 0
	if (config_flags%gaschem_onoff .gt. 0) igaschem_onoff = 1



	if ((igaschem_onoff .le. 0) .or. (iphotol_onoff  .le. 0)) then
	    iradical_onoff = 0
	else
	    iradical_onoff = 1
	end if

	iradical_onoff = 0


	do 3920 jt = jts, jte
	do 3910 it = its, ite

	do 3800 kt = kts, kte

	qcldwtr = moist(it,kt,jt,p_qc)
	if (qcldwtr .le. qcldwtr_cutoff) goto 3800


	icase = icase + 1


	if (ktau .eq. -13579) then


	  call sorgam_cloudchem_dumpaa(   &
	    id, ktau, ktauc, dtstepc, config_flags,   &
	    p_phy, t_phy, rho_phy, alt,   &
	    cldfra, ph_no2,   &
	    moist, chem,   &
	    gas_aqfrac, numgas_aqfrac,   &
	    ids,ide, jds,jde, kds,kde,   &
	    ims,ime, jms,jme, kms,kme,   &
	    its,ite, jts,jte, kts,kte,   &
	    qcldwtr_cutoff,   &
	    it, jt, kt )
	end if


	rbox(1:num_chem) = chem(it,kt,jt,1:num_chem)



	call sorgam_cloudchem_1box(   &
	    id, ktau, ktauc, dtstepc,   &
	    iphotol_onoff, iradical_onoff,   &
	    ph_no2(it,kt,jt),   &
	    ph_aq_box, gas_aqfrac_box,   &
	    numgas_aqfrac, it, jt, kt, icase,   &
	    rbox, qcldwtr,   &
	    t_phy(it,kt,jt), p_phy(it,kt,jt), rho_phy(it,kt,jt),&
            config_flags)


	chem(it,kt,jt,1:num_chem) = rbox(1:num_chem)
	gas_aqfrac(it,kt,jt,:) = gas_aqfrac_box(:) 


3800	continue

3910	continue
3920	continue

	print 93010, 'leaving  sorgam_vbs_cloudchem_driver - ktau =', ktau, icase
93010	format( a, 8(1x,i6) )

	return
	end subroutine sorgam_vbs_cloudchem_driver




	subroutine sorgam_cloudchem_1box(   &
	    id, ktau, ktauc, dtstepc,   &
	    iphotol_onoff, iradical_onoff,   &
	    photol_no2_box,   &
	    ph_aq_box, gas_aqfrac_box,   &
	    numgas_aqfrac, it, jt, kt, icase,   &
	    rbox, qcw_box, temp_box, pres_box, rho_box, &
            config_flags )

        use module_configure, only: grid_config_rec_type

	use module_state_description, only:   &
		num_moist, num_chem

	use module_data_sorgam_vbs, only:   &
		msectional, maxd_asize, maxd_atype,   &
		cw_phase, nsize_aer, ntype_aer, do_cloudchem_aer,   &
		lptr_so4_aer, lptr_no3_aer, lptr_nh4_aer,   &
		lptr_orgpa_aer, lptr_ec_aer, lptr_p25_aer, &
                lptr_cl_aer, lptr_na_aer

	use module_data_cmu_bulkaqchem, only:   &
	        meqn1max


	implicit none



        type(grid_config_rec_type), intent(in) :: config_flags

	integer, intent(in) ::   &
		id, ktau, ktauc,   &
		numgas_aqfrac, it, jt, kt,   &
		icase, iphotol_onoff, iradical_onoff

	real, intent(in) ::   &
	    dtstepc, photol_no2_box,   &
	    qcw_box,   &		
	    temp_box,   &		
	    pres_box,   &		
	    rho_box			

	real, intent(inout) :: ph_aq_box

        real, intent(inout), dimension( num_chem ) :: rbox


        real, intent(inout), dimension( numgas_aqfrac ) :: gas_aqfrac_box


	integer :: iphase
	integer :: icase_in, idecomp_hmsa_hso5,   &
		iradical_in, istat_aqop

	integer :: lptr_yyy_cwaer(maxd_asize,maxd_atype,nyyy)

	real :: co2_mixrat_in
	real :: ph_cmuaq_cur
	real :: photol_no2_in
	real :: xprescribe_ph

        real :: yaq_beg(meqn1max), yaq_end(meqn1max)
	real :: rbox_sv1(num_chem)
	real :: rbulk_cwaer(nyyy,2)

        real, dimension( maxd_asize, maxd_atype ) :: fr_partit_cw
        real, dimension( 2, 3 ) :: xvol_old





	iphase = cw_phase
	lptr_yyy_cwaer(:,:,l_so4_aqyy) = lptr_so4_aer(:,:,iphase)
	lptr_yyy_cwaer(:,:,l_no3_aqyy) = lptr_no3_aer(:,:,iphase)
	lptr_yyy_cwaer(:,:,l_nh4_aqyy) = lptr_nh4_aer(:,:,iphase)
	lptr_yyy_cwaer(:,:,l_oin_aqyy) = lptr_p25_aer(:,:,iphase)
	lptr_yyy_cwaer(:,:,l_bc_aqyy ) = lptr_ec_aer( :,:,iphase)
	lptr_yyy_cwaer(:,:,l_oc_aqyy ) = lptr_orgpa_aer( :,:,iphase)

        lptr_yyy_cwaer(:,:,l_cl_aqyy ) = lptr_cl_aer(:,:,iphase)
        lptr_yyy_cwaer(:,:,l_na_aqyy ) = lptr_na_aer(:,:,iphase)








	icase_in = icase
	iradical_in = 1
	idecomp_hmsa_hso5 = 1

	co2_mixrat_in = 350.0

	photol_no2_in = photol_no2_box


	xprescribe_ph = -1.0e31



	if ((iphotol_onoff .le. 0) .or. (iradical_onoff .le. 0)) then
	    photol_no2_in = 0.0
	    iradical_in = 0
	end if


	rbox_sv1(:) = rbox(:)
	gas_aqfrac_box(:) = 0.0



	call sorgam_interface_to_aqoperator1(   &
	    istat_aqop,   &
	    dtstepc,   &
	    rbox, gas_aqfrac_box,   &
	    qcw_box, temp_box, pres_box, rho_box,   &
	    rbulk_cwaer, lptr_yyy_cwaer,   &
	    co2_mixrat_in, photol_no2_in, xprescribe_ph,   &
	    iradical_in, idecomp_hmsa_hso5,   &
	    yaq_beg, yaq_end, ph_cmuaq_cur,   &
	    numgas_aqfrac, id, it, jt, kt, ktau, icase_in, &
            config_flags )

	ph_aq_box = ph_cmuaq_cur









	call sorgam_partition_cldwtr(   &
	    rbox, fr_partit_cw, xvol_old,   &
	    id, it, jt, kt, icase_in )






	call sorgam_distribute_bulk_changes(   &
	    rbox, rbox_sv1, fr_partit_cw,   &
	    rbulk_cwaer, lptr_yyy_cwaer,   &
	    id, it, jt, kt, icase_in )





	if (msectional .lt. 1000000000) then
	    call sorgam_cloudchem_apply_mode_transfer(   &
		rbox, rbox_sv1, xvol_old,   &
		id, it, jt, kt, icase_in )
	end if



	return
	end subroutine sorgam_cloudchem_1box




	subroutine sorgam_interface_to_aqoperator1(   &
	    istat_aqop,   &
	    dtstepc,   &
	    rbox, gas_aqfrac_box,   &
	    qcw_box, temp_box, pres_box, rho_box,   &
	    rbulk_cwaer, lptr_yyy_cwaer,   &
	    co2_mixrat_in, photol_no2_in, xprescribe_ph,   &
	    iradical_in, idecomp_hmsa_hso5,   &
	    yaq_beg, yaq_end, ph_cmuaq_cur,   &
	    numgas_aqfrac, id, it, jt, kt, ktau, icase, &
            config_flags )

        use module_configure, only: grid_config_rec_type

	use module_state_description, only:   &
		num_chem, param_first_scalar, p_qc,   &
		p_nh3, p_hno3, p_hcl, p_sulf, p_hcho,   &
		p_ora1, p_so2, p_h2o2, p_o3, p_ho,   &
		p_ho2, p_no3, p_no, p_no2, p_hono,   &
		p_pan, p_ch3o2, p_ch3oh, p_op1, &
                p_form, p_facd, p_oh, p_meo2, p_meoh, p_mepx, &
                CB05_SORG_VBS_AQ_KPP

	use module_data_cmu_bulkaqchem, only:   &
	        meqn1max, naers, ngas,   &
	        na4, naa, nac, nae, nah, nahmsa, nahso5,   &
	        nan, nao, nar, nas, naw,   &
	        ng4, nga, ngc, ngch3co3h, ngch3o2, ngch3o2h, ngch3oh,   &
	        ngh2o2, nghcho, nghcooh, nghno2, ngho2,   &
	        ngn, ngno, ngno2, ngno3, ngo3, ngoh, ngpan, ngso2

	use module_cmu_bulkaqchem, only:  aqoperator1

	use module_data_sorgam_vbs, only:   &
        	maxd_asize, maxd_atype,   &
		cw_phase, nsize_aer, ntype_aer, do_cloudchem_aer,   &
		lptr_so4_aer, lptr_no3_aer, lptr_nh4_aer,   &
		lptr_orgpa_aer, lptr_ec_aer, lptr_p25_aer, &
                lptr_cl_aer, lptr_na_aer


	implicit none



        type(grid_config_rec_type), intent(in) :: config_flags

	integer, intent(in) ::   &
		iradical_in, idecomp_hmsa_hso5,   &
	        numgas_aqfrac, id, it, jt, kt, ktau, icase
	integer, intent(inout) ::   &
		istat_aqop

	integer, intent(in) :: lptr_yyy_cwaer(maxd_asize,maxd_atype,nyyy)

	real, intent(in) ::   &
	    dtstepc, co2_mixrat_in,   &
	    photol_no2_in, xprescribe_ph,   &
	    qcw_box, temp_box, pres_box, rho_box

	real, intent(inout) :: ph_cmuaq_cur

        real, intent(inout), dimension( num_chem ) :: rbox	

        real, intent(inout), dimension( numgas_aqfrac ) :: gas_aqfrac_box

        real, intent(inout), dimension( nyyy, 2 ) :: rbulk_cwaer

        real, intent(inout), dimension( meqn1max ) :: yaq_beg, yaq_end



	integer :: i, iphase, isize, itype
	integer :: iaq, istat_fatal, istat_warn
	integer :: l, lyyy
	integer :: p1st

	real, parameter :: eps=0.622   


	real :: cair_moleperm3
	real :: dum, dumb
	real :: factgas, factlwc, factpatm, factphoto
	real :: factaerbc, factaercl, factaerna, factaernh4,   &
		factaerno3, factaeroc, factaeroin, factaerso4
	real :: lwc
	real :: p_atm, photo_in
	real :: rh
	real :: temp, tstep_beg_sec, tstep_end_sec
	real :: totsulf_beg, totsulf_end
	real :: gas(ngas), aerosol(naers)
	real :: gas_aqfrac_cmu(ngas)

	double precision tstep_beg_sec_dp, tstep_end_sec_dp,   &
	  temp_dp, p_atm_dp, lwc_dp, rh_dp,   &
	  co2_mixrat_in_dp, photo_in_dp, ph_cmuaq_cur_dp,   &
	  xprescribe_ph_dp
	double precision gas_dp(ngas), gas_aqfrac_cmu_dp(ngas),   &
	  aerosol_dp(naers), yaq_beg_dp(meqn1max), yaq_end_dp(meqn1max)



	p1st = param_first_scalar






	factpatm = 1.0/1.01325e5

	factlwc = 1.0e3*rho_box

	factphoto = 1.6


	factgas = 1.0


	dum = rho_box
	factaerso4   = dum
	factaerno3   = dum
	factaercl    = dum
	factaernh4   = dum
	factaerna    = dum
	factaeroin   = dum
	factaeroc    = dum
	factaerbc    = dum





	temp = temp_box

	lwc = qcw_box * factlwc
	p_atm = pres_box * factpatm





	p_atm = (rho_box/28.966)*0.082058e0*temp

	photo_in = photol_no2_in * factphoto

	rh = 1.0
	iaq = 1

	tstep_beg_sec = 0.0
	tstep_end_sec = dtstepc


	gas(:) = 0.0

	if (p_nh3   >= p1st) gas(nga     ) = rbox(p_nh3  )*factgas
	if (p_hno3  >= p1st) gas(ngn     ) = rbox(p_hno3 )*factgas
	if (p_hcl   >= p1st) gas(ngc     ) = rbox(p_hcl  )*factgas
	if (p_sulf  >= p1st) gas(ng4     ) = rbox(p_sulf )*factgas



	if (p_so2   >= p1st) gas(ngso2   ) = rbox(p_so2  )*factgas
	if (p_h2o2  >= p1st) gas(ngh2o2  ) = rbox(p_h2o2 )*factgas
	if (p_o3    >= p1st) gas(ngo3    ) = rbox(p_o3   )*factgas

	if (p_ho2   >= p1st) gas(ngho2   ) = rbox(p_ho2  )*factgas
	if (p_no3   >= p1st) gas(ngno3   ) = rbox(p_no3  )*factgas

	if (p_no    >= p1st) gas(ngno    ) = rbox(p_no   )*factgas
	if (p_no2   >= p1st) gas(ngno2   ) = rbox(p_no2  )*factgas
	if (p_hono  >= p1st) gas(nghno2  ) = rbox(p_hono )*factgas
	if (p_pan   >= p1st) gas(ngpan   ) = rbox(p_pan  )*factgas





        if ((config_flags%chem_opt == CB05_SORG_VBS_AQ_KPP)) then

          if (p_form  >= p1st) gas(nghcho  ) = rbox(p_form )*factgas
          if (p_facd  >= p1st) gas(nghcooh ) = rbox(p_facd )*factgas
          if (p_oh    >= p1st) gas(ngoh    ) = rbox(p_oh   )*factgas
          if (p_meo2  >= p1st) gas(ngch3o2 ) = rbox(p_meo2 )*factgas
          if (p_meoh  >= p1st) gas(ngch3oh ) = rbox(p_meoh )*factgas
          if (p_mepx  >= p1st) gas(ngch3o2h) = rbox(p_mepx )*factgas

        else

          if (p_hcho  >= p1st) gas(nghcho  ) = rbox(p_hcho )*factgas
          if (p_ora1  >= p1st) gas(nghcooh ) = rbox(p_ora1 )*factgas
          if (p_ho    >= p1st) gas(ngoh    ) = rbox(p_ho   )*factgas
          if (p_ch3o2 >= p1st) gas(ngch3o2 ) = rbox(p_ch3o2)*factgas
          if (p_ch3oh >= p1st) gas(ngch3oh ) = rbox(p_ch3oh)*factgas
          if (p_op1   >= p1st) gas(ngch3o2h) = rbox(p_op1  )*factgas

        endif



	aerosol(:) = 0.0
	rbulk_cwaer(:,:) = 0.0

	iphase = cw_phase
	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	if ( .not. do_cloudchem_aer(isize,itype) ) cycle

	do lyyy = 1, nyyy

	l = lptr_yyy_cwaer(isize,itype,lyyy)
	if (l .ge. p1st) rbulk_cwaer(lyyy,1) = rbulk_cwaer(lyyy,1) + rbox(l)

	end do

	end do
	end do


	aerosol(na4) = rbulk_cwaer(l_so4_aqyy,1) * factaerso4
	aerosol(nan) = rbulk_cwaer(l_no3_aqyy,1) * factaerno3
	aerosol(nac) = rbulk_cwaer(l_cl_aqyy, 1) * factaercl
	aerosol(naa) = rbulk_cwaer(l_nh4_aqyy,1) * factaernh4
	aerosol(nas) = rbulk_cwaer(l_na_aqyy, 1) * factaerna
	aerosol(nar) = rbulk_cwaer(l_oin_aqyy,1) * factaeroin
	aerosol(nae) = rbulk_cwaer(l_bc_aqyy, 1) * factaerbc
	aerosol(nao) = rbulk_cwaer(l_oc_aqyy, 1) * factaeroc









	tstep_beg_sec_dp = 0.0d0
	if (tstep_beg_sec .ne. 0.0) tstep_beg_sec_dp = tstep_beg_sec
	tstep_end_sec_dp = 0.0d0
	if (tstep_end_sec .ne. 0.0) tstep_end_sec_dp = tstep_end_sec
	temp_dp = 0.0d0
	if (temp .ne. 0.0) temp_dp = temp
	p_atm_dp = 0.0d0
	if (p_atm .ne. 0.0) p_atm_dp = p_atm
	lwc_dp = 0.0d0
	if (lwc .ne. 0.0) lwc_dp = lwc
	rh_dp = 0.0d0
	if (rh .ne. 0.0) rh_dp = rh
	co2_mixrat_in_dp = 0.0d0
	if (co2_mixrat_in .ne. 0.0) co2_mixrat_in_dp = co2_mixrat_in
	photo_in_dp = 0.0d0
	if (photo_in .ne. 0.0) photo_in_dp = photo_in
	xprescribe_ph_dp = 0.0d0
	if (xprescribe_ph .ne. 0.0) xprescribe_ph_dp = xprescribe_ph
	ph_cmuaq_cur_dp = 0.0d0
	if (ph_cmuaq_cur .ne. 0.0) ph_cmuaq_cur_dp = ph_cmuaq_cur

	do i = 1, ngas
	    gas_dp(i) = 0.0d0
	    if (gas(i) .ne. 0.0) gas_dp(i) = gas(i)
	end do
	do i = 1, naers
	    aerosol_dp(i) = 0.0d0
	    if (aerosol(i) .ne. 0.0) aerosol_dp(i) = aerosol(i)
	end do
	do i = 1, ngas
	    gas_aqfrac_cmu_dp(i) = 0.0d0
	    if (gas_aqfrac_cmu(i) .ne. 0.0) gas_aqfrac_cmu_dp(i) = gas_aqfrac_cmu(i)
	end do
	do i = 1, meqn1max
	    yaq_beg_dp(i) = 0.0d0
	    if (yaq_beg(i) .ne. 0.0) yaq_beg_dp(i) = yaq_beg(i)
	end do
	do i = 1, meqn1max
	    yaq_end_dp(i) = 0.0d0
	    if (yaq_end(i) .ne. 0.0) yaq_end_dp(i) = yaq_end(i)
	end do



	cair_moleperm3 = 1.0e3*p_atm_dp/(0.082058e0*temp_dp)
	totsulf_beg = ( aerosol_dp(na4)/96.   &
	              + aerosol_dp(nahso5)/113. + aerosol_dp(nahmsa)/111.   &
	              + (gas_dp(ngso2) + gas_dp(ng4))*cair_moleperm3 )*96.0










	call aqoperator1(   &
	    istat_fatal, istat_warn,   &
	    tstep_beg_sec_dp, tstep_end_sec_dp,   &
	    gas_dp, aerosol_dp, gas_aqfrac_cmu_dp,   &
	    temp_dp, p_atm_dp, lwc_dp, rh_dp,   &
	    co2_mixrat_in_dp, photo_in_dp, xprescribe_ph_dp,   &
	    iradical_in, idecomp_hmsa_hso5, iaq,   &
	    yaq_beg_dp, yaq_end_dp, ph_cmuaq_cur_dp )

	totsulf_end = ( aerosol_dp(na4)/96.   &
	              + aerosol_dp(nahso5)/113. + aerosol_dp(nahmsa)/111.   &
	              + (gas_dp(ngso2) + gas_dp(ng4))*cair_moleperm3 )*96.0



	tstep_beg_sec = tstep_beg_sec_dp
	tstep_end_sec = tstep_end_sec_dp
	temp = temp_dp
	p_atm = p_atm_dp
	lwc = lwc_dp
	rh = rh_dp



	ph_cmuaq_cur = ph_cmuaq_cur_dp

	do i = 1, ngas
	    gas(i) = gas_dp(i)
	end do
	do i = 1, naers
	    aerosol(i) = aerosol_dp(i)
	end do
	do i = 1, ngas
	    gas_aqfrac_cmu(i) = gas_aqfrac_cmu_dp(i)
	end do
	do i = 1, meqn1max
	    yaq_beg(i) = yaq_beg_dp(i)
	end do
	do i = 1, meqn1max
	    yaq_end(i) = yaq_end_dp(i)
	end do





	istat_aqop = 0
	if (istat_fatal .ne. 0) then
	    write(6,*)   &
		'*** sorgam_cloudchem_driver, subr interface_to_aqoperator1'
	    write(6,'(a,4i5,2i10)')   &
		'    id,it,jt,kt, istat_fatal, warn =',   &
		id, it, jt, kt, istat_fatal, istat_warn
	    istat_aqop = -10
	end if





	dum = totsulf_end - totsulf_beg
	dumb = max( totsulf_beg, totsulf_end )
	if (abs(dum) .gt. max(1.0e-3,1.0e-3*dumb)) then
	    write(6,*)   &
		'*** sorgam_cloudchem_driver, sulfur balance warning'
	    write(6,'(a,4i5,1p,3e12.4)')   &
		'    id,it,jt,kt, total_sulfur_beg, _end, _error =',   &
		id, it, jt, kt, totsulf_beg, totsulf_end, dum
	end if




	gas_aqfrac_box(:) = 0.0

	if (p_nh3   >= p1st) then
	    rbox(p_nh3  ) = gas(nga     )/factgas
	    if (p_nh3   <= numgas_aqfrac)   &
		gas_aqfrac_box(p_nh3   ) = gas_aqfrac_cmu(nga     )
	end if
	if (p_hno3  >= p1st) then
	    rbox(p_hno3 ) = gas(ngn     )/factgas
	    if (p_hno3  <= numgas_aqfrac)   &
		gas_aqfrac_box(p_hno3  ) = gas_aqfrac_cmu(ngn     )
	end if
	if (p_hcl   >= p1st) then
	    rbox(p_hcl  ) = gas(ngc     )/factgas
	    if (p_hcl   <= numgas_aqfrac)   &
		gas_aqfrac_box(p_hcl   ) = gas_aqfrac_cmu(ngc     )
	end if
	if (p_sulf  >= p1st) then
	    rbox(p_sulf ) = gas(ng4     )/factgas
	    if (p_sulf  <= numgas_aqfrac)   &
		gas_aqfrac_box(p_sulf  ) = gas_aqfrac_cmu(ng4     )
	end if











	if (p_so2   >= p1st) then
	    rbox(p_so2  ) = gas(ngso2   )/factgas
	    if (p_so2   <= numgas_aqfrac)   &
		gas_aqfrac_box(p_so2   ) = gas_aqfrac_cmu(ngso2   )
	end if
	if (p_h2o2  >= p1st) then
	    rbox(p_h2o2 ) = gas(ngh2o2  )/factgas
	    if (p_h2o2  <= numgas_aqfrac)   &
		gas_aqfrac_box(p_h2o2  ) = gas_aqfrac_cmu(ngh2o2  )
	end if
	if (p_o3    >= p1st) then
	    rbox(p_o3   ) = gas(ngo3    )/factgas
	    if (p_o3    <= numgas_aqfrac)   &
		gas_aqfrac_box(p_o3    ) = gas_aqfrac_cmu(ngo3    )
	end if





	if (p_ho2   >= p1st) then
	    rbox(p_ho2  ) = gas(ngho2   )/factgas
	    if (p_ho2   <= numgas_aqfrac)   &
		gas_aqfrac_box(p_ho2   ) = gas_aqfrac_cmu(ngho2   )
	end if
	if (p_no3   >= p1st) then
	    rbox(p_no3  ) = gas(ngno3   )/factgas
	    if (p_no3   <= numgas_aqfrac)   &
		gas_aqfrac_box(p_no3   ) = gas_aqfrac_cmu(ngno3   )
	end if

	if (p_no    >= p1st) then
	    rbox(p_no   ) = gas(ngno    )/factgas
	    if (p_no    <= numgas_aqfrac)   &
		gas_aqfrac_box(p_no    ) = gas_aqfrac_cmu(ngno    )
	end if
	if (p_no2   >= p1st) then
	    rbox(p_no2  ) = gas(ngno2   )/factgas
	    if (p_no2   <= numgas_aqfrac)   &
		gas_aqfrac_box(p_no2   ) = gas_aqfrac_cmu(ngno2   )
	end if
	if (p_hono  >= p1st) then
	    rbox(p_hono ) = gas(nghno2  )/factgas
	    if (p_hono  <= numgas_aqfrac)   &
		gas_aqfrac_box(p_hono  ) = gas_aqfrac_cmu(nghno2  )
	end if
	if (p_pan   >= p1st) then
	    rbox(p_pan  ) = gas(ngpan   )/factgas
	    if (p_pan   <= numgas_aqfrac)   &
		gas_aqfrac_box(p_pan   ) = gas_aqfrac_cmu(ngpan   )
	end if

















       if ( (config_flags%chem_opt == CB05_SORG_VBS_AQ_KPP) ) then

          if (p_form  >= p1st) then
              rbox(p_form ) = gas(nghcho  )/factgas
              if (p_form  <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_form  ) = gas_aqfrac_cmu(nghcho  )
          end if
          if (p_facd  >= p1st) then
              rbox(p_facd ) = gas(nghcooh )/factgas
              if (p_facd  <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_facd  ) = gas_aqfrac_cmu(nghcooh )
          end if
          if (p_oh    >= p1st) then
              rbox(p_oh   ) = gas(ngoh    )/factgas
              if (p_oh    <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_oh    ) = gas_aqfrac_cmu(ngoh    )
          end if
          if (p_meo2  >= p1st) then
              rbox(p_meo2 ) = gas(ngch3o2 )/factgas
              if (p_meo2  <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_meo2  ) = gas_aqfrac_cmu(ngch3o2 )
          end if

          if (p_meoh  >= p1st) then
              rbox(p_meoh ) = gas(ngch3oh )/factgas
              if (p_meoh  <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_meoh  ) = gas_aqfrac_cmu(ngch3oh )
          end if
          if (p_mepx  >= p1st) then
              rbox(p_mepx ) = gas(ngch3o2h)/factgas
              if (p_mepx  <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_mepx  ) = gas_aqfrac_cmu(ngch3o2h)
          end if

       else
          if (p_hcho  >= p1st) then
              rbox(p_hcho ) = gas(nghcho  )/factgas
              if (p_hcho  <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_hcho  ) = gas_aqfrac_cmu(nghcho  )
          end if
          if (p_ora1  >= p1st) then
              rbox(p_ora1 ) = gas(nghcooh )/factgas
              if (p_ora1  <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_ora1  ) = gas_aqfrac_cmu(nghcooh )
          end if
          if (p_ho    >= p1st) then
              rbox(p_ho   ) = gas(ngoh    )/factgas
              if (p_ho    <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_ho    ) = gas_aqfrac_cmu(ngoh    )
          end if
          if (p_ch3o2 >= p1st) then
              rbox(p_ch3o2) = gas(ngch3o2 )/factgas
              if (p_ch3o2 <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_ch3o2 ) = gas_aqfrac_cmu(ngch3o2 )
          end if
          if (p_ch3oh >= p1st) then
              rbox(p_ch3oh) = gas(ngch3oh )/factgas
              if (p_ch3oh <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_ch3oh ) = gas_aqfrac_cmu(ngch3oh )
          end if
          if (p_op1   >= p1st) then
              rbox(p_op1  ) = gas(ngch3o2h)/factgas
              if (p_op1   <= numgas_aqfrac)   &
                  gas_aqfrac_box(p_op1   ) = gas_aqfrac_cmu(ngch3o2h)
          end if

       end if



	rbulk_cwaer(l_so4_aqyy,2) = aerosol(na4)/factaerso4
	rbulk_cwaer(l_no3_aqyy,2) = aerosol(nan)/factaerno3
	rbulk_cwaer(l_cl_aqyy, 2) = aerosol(nac)/factaercl
	rbulk_cwaer(l_nh4_aqyy,2) = aerosol(naa)/factaernh4
	rbulk_cwaer(l_na_aqyy, 2) = aerosol(nas)/factaerna
	rbulk_cwaer(l_oin_aqyy,2) = aerosol(nar)/factaeroin
	rbulk_cwaer(l_bc_aqyy, 2) = aerosol(nae)/factaerbc
	rbulk_cwaer(l_oc_aqyy, 2) = aerosol(nao)/factaeroc




	return
	end subroutine sorgam_interface_to_aqoperator1




	subroutine sorgam_partition_cldwtr(   &
	    rbox, fr_partit_cw, xvol_old,   &
	    id, it, jt, kt, icase )

	use module_state_description, only:   &
		param_first_scalar, num_chem

	use module_data_sorgam_vbs, only:   &
        	maxd_asize, maxd_atype,   &
		ai_phase, cw_phase, nsize_aer, ntype_aer, ncomp_aer,   &
		do_cloudchem_aer, massptr_aer, numptr_aer,   &
		dens_aer, sigmag_aer,   &
		dcen_sect, dlo_sect, dhi_sect,   &
		volumcen_sect, volumlo_sect, volumhi_sect


	implicit none


	integer, intent(in) :: id, it, jt, kt, icase

        real, intent(inout), dimension( 1:num_chem ) :: rbox

        real, intent(inout), dimension( maxd_asize, maxd_atype ) ::   &
		fr_partit_cw

        real, intent(inout), dimension( 2, 3 ) :: xvol_old


	integer :: isize, itype
	integer :: jdone_mass, jdone_numb, jpos, jpos_mass, jpos_numb
	integer :: l, ll
	integer :: p1st

	real, parameter :: partit_wght_mass = 0.5

	real :: tmpa, tmpb, tmpc
	real :: tmp_cwvolfrac, tmp_lnsg
	real :: tmass, tnumb, umass, unumb, wmass, wnumb
	real :: xmass_c, xmass_a, xmass_t, xvolu_c, xvolu_a, xvolu_t
	real :: xnumb_c1, xnumb_a1, xnumb_t1, xnumb_c2, xnumb_a2, xnumb_t2
	real, dimension( maxd_asize, maxd_atype ) :: fmass, fnumb, xmass, xnumb, xnumbsv


	p1st = PARAM_FIRST_SCALAR

	tmass = 0.0
	tnumb = 0.0
	umass = 0.0
	unumb = 0.0









	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    if ( .not. do_cloudchem_aer(isize,itype) ) cycle
	    xmass_c = 0.0
	    xvolu_c = 0.0
	    xvolu_a = 0.0
	    do ll = 1, ncomp_aer(itype)
		l = massptr_aer(ll,isize,itype,cw_phase)
		if (l .ge. p1st) then
		    tmpa = max( 0.0, rbox(l) )*1.0e-6
		    xmass_c = xmass_c + tmpa
		    xvolu_c = xvolu_c + tmpa/dens_aer(ll,itype)
		end if
		l = massptr_aer(ll,isize,itype,ai_phase)
		if (l .ge. p1st) then
		    tmpa = max( 0.0, rbox(l) )*1.0e-6
		    xvolu_a = xvolu_a + tmpa/dens_aer(ll,itype)
		end if
	    end do

	    xnumb_c1 = max( 0.0, rbox(numptr_aer(isize,itype,cw_phase)) )
	    xnumb_a1 = max( 0.0, rbox(numptr_aer(isize,itype,ai_phase)) )
	    xnumbsv(isize,itype) = xnumb_c1
	    xnumb_t1 = xnumb_a1 + xnumb_c1
	    xvolu_t = xvolu_a + xvolu_c



	    if (xvolu_t < smallvolaa) then
		xnumb_t2 = xvolu_t/volumcen_sect(isize,itype)
	    else if (xnumb_t1 < xvolu_t/volumhi_sect(isize,itype)) then
		xnumb_t2 = xvolu_t/volumhi_sect(isize,itype)
	    else if (xnumb_t1 > xvolu_t/volumlo_sect(isize,itype)) then
		xnumb_t2 = xvolu_t/volumlo_sect(isize,itype)
	    else
		xnumb_t2 = xnumb_t1
	    end if



	    tmp_cwvolfrac = xvolu_c/max(xvolu_t,1.e-30)
	    tmp_lnsg = log(sigmag_aer(isize,itype))
	    if ((xvolu_c < smallvolaa) .or. (tmp_cwvolfrac < 1.0e-10)) then


		xnumb_c2 = xnumb_t2*tmp_cwvolfrac
		tmpa = -7.0 ; tmpb = -7.0 ; tmpc = -7.0
	    else


		tmpa = norm01_uptail_inv( tmp_cwvolfrac )

		tmpb = tmpa + 3.0*tmp_lnsg

		tmpc = norm01_uptail( tmpb )



		xnumb_c2 = max( xnumb_c1,  xnumb_t2*tmpc )



		xnumb_c2 = min( xnumb_c2,  xnumb_t2*tmp_cwvolfrac )
	    end if
	    xnumb_a2 = xnumb_t2 - xnumb_c2
	    rbox(numptr_aer(isize,itype,cw_phase)) = xnumb_c2
	    rbox(numptr_aer(isize,itype,ai_phase)) = xnumb_a2


	    if (xmass_c .lt. 1.0e-37) xmass_c = 0.0
	    xmass(isize,itype) = xmass_c
	    if (xnumb_c2 .lt. 1.0e-37) xnumb_c2 = 0.0
	    xnumb(isize,itype) = xnumb_c2
	    xnumbsv(isize,itype) = xnumb_c1

	    tmass = tmass + xmass(isize,itype)
	    tnumb = tnumb + xnumb(isize,itype)
	    umass = max( umass, xmass(isize,itype) )
	    unumb = max( unumb, xnumb(isize,itype) )

	    if ((itype == 1) .and. (isize <= 2)) then
		xvol_old(isize,1) = xvolu_c
		xvol_old(isize,2) = xvolu_a
		xvol_old(isize,3) = xvolu_t
	    end if
	end do
	end do





	jdone_mass = 0
	jdone_numb = 0
	jpos_mass = 0
	jpos_numb = 0
	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    if ( .not. do_cloudchem_aer(isize,itype) ) cycle
	    fmass(isize,itype) = 0.0
	    if (tmass .ge. 1.0e-35) then
		fmass(isize,itype) = xmass(isize,itype)/tmass
	    else if (umass .gt. 0.0) then
		if ( (jdone_mass .eq. 0) .and.   &
		     (xmass(isize,itype) .eq. umass) ) then
		    jdone_mass = 1
		    fmass(isize,itype) = 1.0
		end if
	    end if
	    if (fmass(isize,itype) .gt. 0) jpos_mass = jpos_mass + 1

	    fnumb(isize,itype) = 0.0
	    if (tnumb .ge. 1.0e-35) then
		fnumb(isize,itype) = xnumb(isize,itype)/tnumb
	    else if (unumb .gt. 0.0) then
		if ( (jdone_numb .eq. 0) .and.   &
		     (xnumb(isize,itype) .eq. unumb) ) then
		    jdone_numb = 1
		    fnumb(isize,itype) = 1.0
		end if
	    end if
	    if (fnumb(isize,itype) .gt. 0) jpos_numb = jpos_numb + 1
	end do
	end do


	if ((jpos_mass .eq. 1) .or. (jpos_numb .eq. 1)) then
	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    if ( .not. do_cloudchem_aer(isize,itype) ) cycle
	    if (jpos_mass .eq. 1) then
		if (fmass(isize,itype) .gt. 0) fmass(isize,itype) = 1.0
	    end if
	    if (jpos_numb .eq. 1) then
		if (fnumb(isize,itype) .gt. 0) fnumb(isize,itype) = 1.0
	    end if
	end do
	end do
	end if








	fr_partit_cw(:,:) = 0.0
	if ((jpos_mass .eq. 0) .and. (jpos_numb .eq. 0)) then
	    itype = 1
	    isize = (nsize_aer(itype)+1)/2
	    fr_partit_cw(isize,itype) = 1.0

	else if (jpos_mass .eq. 0) then
	    fr_partit_cw(:,:) = fnumb(:,:)

	else if (jpos_numb .eq. 0) then
	    fr_partit_cw(:,:) = fmass(:,:)

	else
	    wmass = max( 0.0, min( 1.0, partit_wght_mass ) )
	    wnumb = 1.0 - wmass
	    fr_partit_cw(:,:) =  wmass*fmass(:,:) + wnumb*fnumb(:,:)

	    jpos = 0
	    do itype = 1, ntype_aer
	    do isize = 1, nsize_aer(itype)
		if ( .not. do_cloudchem_aer(isize,itype) ) cycle
		if (fr_partit_cw(isize,itype) .gt. 0.0) jpos = jpos + 1
	    end do
	    end do


	    if (jpos .eq. 1) then
	    do itype = 1, ntype_aer
	    do isize = 1, nsize_aer(itype)
		if ( .not. do_cloudchem_aer(isize,itype) ) cycle
		if (fr_partit_cw(isize,itype) .gt. 0.0)   &
			fr_partit_cw(isize,itype) = 1.0
	    end do
	    end do
	    end if
	end if




	return
	end subroutine sorgam_partition_cldwtr




	subroutine sorgam_distribute_bulk_changes(   &
		rbox, rbox_sv1, fr_partit_cw,   &
		rbulk_cwaer, lptr_yyy_cwaer,   &
		id, it, jt, kt, icase )

	use module_state_description, only:   &
        	param_first_scalar, num_chem

	use module_scalar_tables, only:  chem_dname_table

	use module_data_sorgam_vbs, only:   &
		maxd_asize, maxd_atype,   &
		cw_phase, nsize_aer, ntype_aer, do_cloudchem_aer,   &
		lptr_so4_aer, lptr_no3_aer, lptr_nh4_aer,   &
		lptr_orgpa_aer, lptr_ec_aer, lptr_p25_aer


	implicit none


	integer, intent(in) :: id, it, jt, kt, icase

	integer, intent(in) :: lptr_yyy_cwaer(maxd_asize,maxd_atype,nyyy)

        real, intent(inout), dimension( 1:num_chem ) :: rbox, rbox_sv1

        real, intent(in), dimension( maxd_asize, maxd_atype ) ::   &
		fr_partit_cw

        real, intent(in), dimension( nyyy, 2 ) :: rbulk_cwaer



	integer :: iphase, isize, itype
	integer :: idone, icount, ncount
	integer :: jpos, jpos_sv
	integer :: l, lunxxaa, lunxxbb, lyyy
	integer :: p1st

	real :: duma, dumb, dumc
	real :: fr, frsum_cur
        real :: fr_cur(maxd_asize,maxd_atype)
	real :: del_r_current, del_r_remain
	real :: del_rbulk_cwaer(nyyy)


	p1st = param_first_scalar

	do lyyy = 1, nyyy
	    del_rbulk_cwaer(lyyy) = rbulk_cwaer(lyyy,2) - rbulk_cwaer(lyyy,1)
	end do

	iphase = cw_phase


	jpos = 0
	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    if ( .not. do_cloudchem_aer(isize,itype) ) cycle
	    if (fr_partit_cw(isize,itype) .gt. 0) jpos = jpos + 1
	end do
	end do
	jpos_sv = jpos




	if (jpos_sv .eq. 1) then
	    do lyyy = 1, nyyy

	    do itype = 1, ntype_aer
	    do isize = 1, nsize_aer(itype)
		if ( .not. do_cloudchem_aer(isize,itype) ) cycle
		fr = fr_partit_cw(isize,itype)
		if (fr .eq. 1.0) then
		    l = lptr_yyy_cwaer(isize,itype,lyyy)
		    if (l .ge. p1st) rbox(l) = rbulk_cwaer(lyyy,2)
		end if
	    end do
	    end do

	    end do
	    goto 7900
	end if


	do 3900 lyyy = 1, nyyy




	if (del_rbulk_cwaer(lyyy) .eq. 0.0) then
	    goto 3900
	else if (del_rbulk_cwaer(lyyy) .gt. 0.0) then
	    do itype = 1, ntype_aer
	    do isize = 1, nsize_aer(itype)
		if ( .not. do_cloudchem_aer(isize,itype) ) cycle
		fr = fr_partit_cw(isize,itype)
		if (fr .gt. 0.0) then
		    l = lptr_yyy_cwaer(isize,itype,lyyy)
		    if (l .ge. p1st) then
			rbox(l) = rbox(l) + fr*del_rbulk_cwaer(lyyy)
		    end if
		end if
	    end do
	    end do

	    goto 3900
	end if





	del_r_remain = del_rbulk_cwaer(lyyy)
	fr_cur(:,:) = fr_partit_cw(:,:)

	ncount = max( 1, jpos_sv*2 )
	icount = 0


	do while (icount .le. ncount)

	icount = icount + 1
	del_r_current = del_r_remain
	jpos = 0
	frsum_cur = 0.0

	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	if ( .not. do_cloudchem_aer(isize,itype) ) cycle

	fr = fr_cur(isize,itype)

	if (fr .gt. 0.0) then
	    l = lptr_yyy_cwaer(isize,itype,lyyy)
	    if (l .ge. p1st) then
		duma = fr*del_r_current
		dumb = rbox(l) + duma
		if (dumb .gt. 0.0) then
		    jpos = jpos + 1
		else if (dumb .eq. 0.0) then
		    fr_cur(isize,itype) = 0.0
		else
		    duma = -rbox(l)
		    dumb = 0.0
		    fr_cur(isize,itype) = 0.0
		end if
		del_r_remain = del_r_remain - duma
		rbox(l) = dumb
		frsum_cur = frsum_cur + fr_cur(isize,itype)
	    else
		fr_cur(isize,itype) = 0.0
	    end if
	end if

	end do	
	end do	


	if (jpos .eq. jpos_sv) then
	    idone = 1

	else if (del_r_remain .ge. 0.0) then
	    idone = 2

	else if (abs(del_r_remain) .le. 1.0e-7*abs(del_rbulk_cwaer(lyyy))) then
	    idone = 3

	else if (frsum_cur .le. 0.0) then
	    idone = 4

	else if (jpos .le. 0) then
	    idone = 5
	else
	    idone = 0
	end if


	if (idone .gt. 0) then
	    lunxxaa = 6
	    if ((lunxxaa .gt. 0) .and. (icount .gt. (1+jpos_sv)/2)) then
		write(lunxxaa,9800)   &
		'distribute_bulk_changes - icount>jpos_sv/2 - i,j,k'
		write(lunxxaa,9810)  it, jt, kt
		write(lunxxaa,9800) 'icase, lyyy, idone, icount, jpos, jpos_sv'
		write(lunxxaa,9810)  icase, lyyy, idone, icount, jpos, jpos_sv
	    end if
	    goto 3900	
	end if


	fr_cur(:,:) = fr_cur(:,:)/frsum_cur

	end do	



	lunxxbb = 6
	if (lunxxbb .gt. 0) then
	    write(lunxxbb,9800)
	    write(lunxxbb,9800)   &
		'distribute_bulk_changes - icount>ncount - i,j,k'
	    write(lunxxbb,9810)  it, jt, kt
	    write(lunxxbb,9800) 'icase, lyyy, icount, ncount, jpos_sv, jpos'
	    write(lunxxbb,9810)  icase, lyyy, icount, ncount, jpos_sv, jpos
	    write(lunxxbb,9800) 'rbulk_cwaer(1), del_rbulk_cwaer, del_r_remain, frsum_cur, (frsum_cur-1.0)'
	    write(lunxxbb,9820)  rbulk_cwaer(lyyy,1), del_rbulk_cwaer(lyyy),   &
		del_r_remain, frsum_cur, (frsum_cur-1.0)
	end if
9800	format( a )
9801	format( 3a )
9810	format( 7i10 )
9820	format( 7(1pe10.2) )
9840	format( 2i3, 5(1pe14.6) )


3900	continue

7900	continue




	return
	end subroutine sorgam_distribute_bulk_changes




	subroutine sorgam_cloudchem_apply_mode_transfer(   &
		rbox, rbox_sv1, xvol_old,   &
		id, it, jt, kt, icase )

	use module_state_description, only:   &
        	param_first_scalar, num_chem

	use module_scalar_tables, only:  chem_dname_table

	use module_data_sorgam_vbs, only:   &
		pirs,   &
		msectional,   &
		maxd_asize, maxd_atype,   &
		ai_phase, cw_phase, nsize_aer, ntype_aer, ncomp_aer,   &
		do_cloudchem_aer, massptr_aer, numptr_aer, dens_aer,   &
		sigmag_aer, dcen_sect, dlo_sect, dhi_sect,   &
		volumcen_sect, volumlo_sect, volumhi_sect,   &
		lptr_so4_aer, lptr_nh4_aer, lptr_p25_aer

	use module_aerosols_sorgam_vbs, only:  getaf


	implicit none


	integer, intent(in) :: id, it, jt, kt, icase

        real, intent(inout), dimension( 1:num_chem ) :: rbox, rbox_sv1

        real, intent(in), dimension( 2, 3 ) :: xvol_old



	integer :: idum_msect
	integer :: ii, isize, isize_ait, isize_acc, itype
	integer :: jj
	integer :: l, lfrm, ltoo, ll
	integer :: lptr_dum(maxd_asize,maxd_atype)
	integer :: p1st

	logical :: skip_xfer

	real :: delvol(2)
	real :: fracrem_num, fracrem_vol, fracxfr_num, fracxfr_vol
	real :: rbox_sv2(1:num_chem)
	real :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf
	real :: tmp_cwnumfrac, tmp_cwvolfrac
	real :: tmp_dpmeanvol, tmp_lnsg
	real :: xcut_num, xcut_vol
	real :: xdgnum_aaa(2), xlnsg(2)
	real :: xnum_aaa(2,3)
	real :: xvol_aaa(2,3)



	p1st = param_first_scalar




	skip_xfer = .false.

	do ii = 1, 2
	itype = 1
	isize = ii


	tmpa = 0.0
	do ll = 1, ncomp_aer(itype)
	    l = massptr_aer(ll,isize,itype,cw_phase)
	    if (l >= p1st) tmpa = tmpa + rbox(l)/dens_aer(ll,itype)
	end do
	xvol_aaa(ii,1) = tmpa*1.0e-6
	xvol_aaa(ii,2) = xvol_old(ii,2)
	xnum_aaa(ii,1) = rbox(numptr_aer(isize,itype,cw_phase))
	xnum_aaa(ii,2) = rbox(numptr_aer(isize,itype,ai_phase))

	xvol_aaa(ii,3) = xvol_aaa(ii,1) + xvol_aaa(ii,2)
	xnum_aaa(ii,3) = xnum_aaa(ii,1) + xnum_aaa(ii,2)
	delvol(ii) = xvol_aaa(ii,1) - xvol_old(ii,1)


	if (ii == 1) then
	    if (xvol_aaa(ii,3) < smallvolaa) then
		skip_xfer = .true.
		exit
	    end if
	end if


	if (xvol_aaa(ii,3) < smallvolaa) then
	    tmp_dpmeanvol = dcen_sect(isize,itype)
	else
	    tmp_dpmeanvol = xvol_aaa(ii,3)/xnum_aaa(ii,3)
	    tmp_dpmeanvol = (tmp_dpmeanvol*6.0/pirs)**0.33333333
	end if
	xlnsg(ii) = log(sigmag_aer(isize,itype))
	xdgnum_aaa(ii) = tmp_dpmeanvol*exp(-1.5*xlnsg(ii)*xlnsg(ii))

	tmp_cwvolfrac = xvol_aaa(ii,1)/max(xvol_aaa(ii,3),1.e-30)


	end do  



	if ( skip_xfer ) return
	if (delvol(1) > delvol(2)) then
	    continue
	else if ( (xdgnum_aaa(1) > 0.03e-4) .and. (xnum_aaa(1,3) > xnum_aaa(2,3)) ) then
	    continue
	else
	    return
	end if











	tmpa = sqrt(2.0)

	xcut_num = tmpa * getaf( xnum_aaa(1,3), xnum_aaa(2,3), &
		xdgnum_aaa(1), xdgnum_aaa(2), xlnsg(1), xlnsg(2), tmpa )



	tmpd = xcut_num
	tmpc = 3.0*xlnsg(1)
	xcut_vol = max( xcut_num-tmpc, 0.0 )
	xcut_num = xcut_vol + tmpc
	fracxfr_vol = norm01_uptail( xcut_vol )
	fracxfr_num = norm01_uptail( xcut_num )
	tmpe = fracxfr_vol ; tmpf = fracxfr_num

	tmp_cwvolfrac = xvol_aaa(1,1)/max(xvol_aaa(1,3),1.e-30)
	tmp_cwnumfrac = xnum_aaa(1,1)/max(xnum_aaa(1,3),1.e-30)
	if ( (fracxfr_vol >= tmp_cwvolfrac) .or. &
	     (fracxfr_num >= tmp_cwnumfrac) ) then


	    fracxfr_num = 1.0
	    fracxfr_vol = 1.0
	else



	    fracxfr_vol = fracxfr_vol/max(1.0e-10,tmp_cwvolfrac)
	    fracxfr_num = fracxfr_num/max(1.0e-10,tmp_cwnumfrac)

	    fracxfr_num = min( fracxfr_num, fracxfr_vol )
	end if
	fracrem_vol = 1.0 - fracxfr_vol
	fracrem_num = 1.0 - fracxfr_num
	if ( skip_xfer ) return


	rbox_sv2(:) = rbox(:)
	itype = 1
	isize_ait = 1
	isize_acc = 2

	lfrm = numptr_aer(isize_ait,itype,cw_phase)
	ltoo = numptr_aer(isize_acc,itype,cw_phase)
	rbox(ltoo) = rbox(ltoo) + rbox(lfrm)*fracxfr_num
	rbox(lfrm) = rbox(lfrm)*fracrem_num

	do ll = 1, ncomp_aer(itype)
	    lfrm = massptr_aer(ll,isize_ait,itype,cw_phase)
	    ltoo = massptr_aer(ll,isize_acc,itype,cw_phase)
	    if (lfrm >= p1st) then
		if (ltoo >= p1st) rbox(ltoo) = rbox(ltoo) &
		                             + rbox(lfrm)*fracxfr_vol
		rbox(lfrm) = rbox(lfrm)*fracrem_vol
	    end if
	end do




	return
	end subroutine sorgam_cloudchem_apply_mode_transfer




	real function norm01_uptail( x )






	implicit none
	real x, xabs
	real*8 erfc_approx, tmpa, t, z

	xabs = abs(x)
	if (xabs >= 12.962359) then
	    if (x > 0.0) then
		norm01_uptail = 0.0
	    else
		norm01_uptail = 1.0
	    end if
	    return
	end if

	z = xabs / sqrt(2.0_8)
	t = 1.0_8/(1.0_8 + 0.5_8*z)







	tmpa =  ( -z*z - 1.26551223_8 + t*(1.00002368_8 + t*(0.37409196_8 +   &
      	  t*(0.09678418_8 + t*(-0.18628806_8 + t*(0.27886807_8 +   &
      	                                   t*(-1.13520398_8 +   &
      	  t*(1.48851587_8 + t*(-0.82215223_8 + t*0.17087277_8 )))))))))

	erfc_approx = t * exp(tmpa)
	if (x .lt. 0.0) erfc_approx = 2.0_8 - erfc_approx

	norm01_uptail = 0.5_8 * erfc_approx

	return
	end function norm01_uptail


	real function norm01_uptail_inv( x )






	implicit none


	real x


	integer niter
	real dfdyinv, f, pi, sqrt2pi, tmpa, y, ynew

	parameter (pi = 3.1415926535897932384626434)

	if (x .le. 1.0e-38) then
	    norm01_uptail_inv = 12.962359
	    return
	else if (x .ge. 1.0) then
	    norm01_uptail_inv = -12.962359
	    return
	end if

	sqrt2pi = sqrt( 2.0*pi )





	tmpa = x
	tmpa = max( 0.0, min( 1.0, tmpa ) )
	tmpa = 4.0*tmpa*(1.0 - tmpa)
	tmpa = max( 1.0e-38, min( 1.0, tmpa ) )
	y = sqrt( -(pi/2.0)*log(tmpa) )
	if (x > 0.5) y = -y

	f = norm01_uptail(y) - x
	do niter = 1, 100

	    dfdyinv = -sqrt2pi * exp( 0.5*y*y )
	    ynew = y - f*dfdyinv
	    f = norm01_uptail(ynew) - x

	    if ( (ynew == y) .or.   &
		 (abs(f) <= abs(x)*1.0e-5) ) then
		exit
	    end if
	    y = ynew
	end do
9100	format( 'niter/x/f/y/ynew', i5, 4(1pe16.8) )

	norm01_uptail_inv = ynew
	return
	end function norm01_uptail_inv




	subroutine sorgam_cloudchem_dumpaa(   &
	    id, ktau, ktauc, dtstepc, config_flags,   &
	    p_phy, t_phy, rho_phy, alt,   &
	    cldfra, ph_no2,   &
	    moist, chem,   &
	    gas_aqfrac, numgas_aqfrac,   &
	    ids,ide, jds,jde, kds,kde,   &
	    ims,ime, jms,jme, kms,kme,   &
	    its,ite, jts,jte, kts,kte,   &
	    qcldwtr_cutoff,   &
	    itcur, jtcur, ktcur )

	use module_state_description, only:   &
		num_moist, num_chem, p_qc
	use module_scalar_tables, only:  chem_dname_table
	use module_configure, only:  grid_config_rec_type
	use module_data_sorgam_vbs


	implicit none


	integer, intent(in) ::   &
		id, ktau, ktauc,   &
		numgas_aqfrac,   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte,   &
		itcur, jtcur, ktcur













	type(grid_config_rec_type), intent(in) :: config_flags


	real, intent(in) ::   &
	    dtstepc, qcldwtr_cutoff


        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                p_phy, t_phy, rho_phy, alt, cldfra, ph_no2







        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
                moist



        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
                chem



        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme, numgas_aqfrac ) :: &
                gas_aqfrac




	integer :: it, jt, kt, l, ll, n
	integer :: isize, itype

	real :: dumai, dumcw
	real :: qcldwtr


	it = itcur
	jt = jtcur
	kt = ktcur

	write(*,*)
	write(*,*)
	write(*,*)
	write(*,9100)
	write(*,9102) ktau, it, jt, kt
9100	format( 7('----------') )
9102	format(   &
	'sorgam_cloudchem_dumpaa - ktau, i, j, k =', 4i5 )

	do 2900 itype = 1, ntype_aer
	do 2900 isize = 1, nsize_aer(itype)
	if ( .not. do_cloudchem_aer(isize,itype) ) goto 2900

	write(*,9110) isize
9110	format( / 'isize, itype =', 2i3 /   &
	'  k    cldwtr      mass-ai   numb-ai      mass-cw   numb-cw' )

	do 2800 kt = kte, kts, -1

	dumai = 0.0
	dumcw = 0.0
	do ll = 1, ncomp_aer(itype)
	    l = massptr_aer(ll,isize,itype,1)
	    dumai = dumai + chem(it,kt,jt,l)
	    l = massptr_aer(ll,isize,itype,2)
	    dumcw = dumcw + chem(it,kt,jt,l)
	end do
	write(*,9120) kt,   &
		moist(it,kt,jt,p_qc),   &
		dumai, chem(it,kt,jt,numptr_aer(isize,itype,1)),   &
		dumcw, chem(it,kt,jt,numptr_aer(isize,itype,2))
9120	format( i3, 1p, e10.2, 2(3x, 2e10.2) )

2800	continue
2900	continue

	write(*,*)
	write(*,9100)
	write(*,*)


	kt = ktcur
	if ((ktau .eq. 30) .and. (it .eq. 23) .and.   &
	      (jt .eq.  1) .and. (kt .eq. 11)) then
	    qcldwtr = moist(it,kt,jt,p_qc)
	    write(*,*)
	    write(*,*)
	    write(*,9102) ktau, it, jt, kt
	    write(*,*)
	    write( *, '(3(1pe10.2,3x,a))' )   &
		(chem(it,kt,jt,l), chem_dname_table(id,l)(1:12), l=1,num_chem)
	    write(*,*)
	    write( *, '(3(1pe10.2,3x,a))' )   &
		  p_phy(it,kt,jt), 'p_phy     ',   &
		  t_phy(it,kt,jt), 't_phy     ',   &
		rho_phy(it,kt,jt), 'rho_phy   ',   &
		    alt(it,kt,jt), 'alt       ',   &
		          qcldwtr, 'qcldwtr   ',   &
		   qcldwtr_cutoff, 'qcldwtrcut'
	    write(*,*)
	    write(*,9100)
	    write(*,*)
	end if


	return
	end subroutine sorgam_cloudchem_dumpaa




	end module module_sorgam_vbs_cloudchem
