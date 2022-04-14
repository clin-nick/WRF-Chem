








	module module_mosaic_cloudchem



	integer, parameter :: l_so4_aqyy = 1
	integer, parameter :: l_no3_aqyy = 2
	integer, parameter :: l_cl_aqyy  = 3
	integer, parameter :: l_nh4_aqyy = 4
	integer, parameter :: l_na_aqyy  = 5
	integer, parameter :: l_oin_aqyy = 6
	integer, parameter :: l_bc_aqyy  = 7
	integer, parameter :: l_oc_aqyy  = 8

	integer, parameter :: nyyy = 8



	contains




	subroutine mosaic_cloudchem_driver(   &
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

	use module_data_mosaic_asect, only:  cw_phase, nphase_aer

	use module_data_mosaic_other, only:  k_pegbegin, name

	use module_mosaic_driver, only:  mapaer_tofrom_host


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




	integer :: it, jt, kt, kpeg, k_pegshift, l, mpeg
	integer :: icase
	integer :: igaschem_onoff, iphotol_onoff, iradical_onoff

	real :: gas_aqfrac_box(numgas_aqfrac)
	real :: ph_aq_box
	real, parameter :: qcldwtr_cutoff = 1.0e-6
	real :: qcldwtr



        if ((cw_phase .le. 0) .or. (cw_phase .gt. nphase_aer)) then
            print *, '*** mosaic_cloudchem_driver - cw_phase not active'
            return
        end if

	print 93010, 'entering mosaic_cloudchem_driver - ktau =', ktau

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


	k_pegshift = k_pegbegin - kts
	kpeg = kt + k_pegshift
	mpeg = 1
	icase = icase + 1


	if (ktau .eq. -13579) then


	  call mosaic_cloudchem_dumpaa(   &
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


	call mapaer_tofrom_host( 0,                           &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      kt, kt,                     &
		num_moist, num_chem, moist, chem,             &
		t_phy, p_phy, rho_phy                         )



	call mosaic_cloudchem_1box(   &
	    id, ktau, ktauc, dtstepc,   &
	    iphotol_onoff, iradical_onoff,   &
	    ph_no2(it,kt,jt),   &
	    ph_aq_box, gas_aqfrac_box,   &
	    numgas_aqfrac, it, jt, kt, kpeg, mpeg, icase )


	call mapaer_tofrom_host( 1,                           &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      kt, kt,                     &
		num_moist, num_chem, moist, chem,             &
		t_phy, p_phy, rho_phy                         )

	gas_aqfrac(it,kt,jt,:) = gas_aqfrac_box(:) 


3800	continue

3910	continue
3920	continue

	print 93010, 'leaving  mosaic_cloudchem_driver - ktau =', ktau, icase
93010	format( a, 8(1x,i6) )

	return
	end subroutine mosaic_cloudchem_driver




	subroutine mosaic_cloudchem_1box(   &
	    id, ktau, ktauc, dtstepc,   &
	    iphotol_onoff, iradical_onoff,   &
	    photol_no2_box,   &
	    ph_aq_box, gas_aqfrac_box,   &
	    numgas_aqfrac, it, jt, kt, kpeg, mpeg, icase )

	use module_state_description, only:   &
		num_moist, num_chem

	use module_data_mosaic_asect, only:   &
		msectional,   &
        	maxd_asize, maxd_atype,   &
		cw_phase, nsize_aer, ntype_aer,   &
		lptr_so4_aer, lptr_no3_aer, lptr_cl_aer, lptr_co3_aer,   &
		lptr_msa_aer, lptr_nh4_aer, lptr_na_aer, lptr_ca_aer,   &
		lptr_oin_aer, lptr_bc_aer, lptr_oc_aer

	use module_data_mosaic_other, only:   &
		l2maxd, ltot2, rsub

	use module_data_cmu_bulkaqchem, only:   &
	        meqn1max


	implicit none


	integer, intent(in) ::   &
		id, ktau, ktauc,   &
		numgas_aqfrac, it, jt, kt, kpeg, mpeg,   &
		icase, iphotol_onoff, iradical_onoff

	real, intent(in) ::   &
	    dtstepc, photol_no2_box

	real, intent(inout) :: ph_aq_box

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
	real :: rbox(l2maxd), rbox_sv1(l2maxd)
	real :: rbulk_cwaer(nyyy,2)

        real, dimension( maxd_asize, maxd_atype ) :: fr_partit_cw





	iphase = cw_phase
	lptr_yyy_cwaer(:,:,l_so4_aqyy) = lptr_so4_aer(:,:,iphase)
	lptr_yyy_cwaer(:,:,l_no3_aqyy) = lptr_no3_aer(:,:,iphase)
	lptr_yyy_cwaer(:,:,l_cl_aqyy ) = lptr_cl_aer( :,:,iphase)
	lptr_yyy_cwaer(:,:,l_nh4_aqyy) = lptr_nh4_aer(:,:,iphase)
	lptr_yyy_cwaer(:,:,l_na_aqyy ) = lptr_na_aer( :,:,iphase)
	lptr_yyy_cwaer(:,:,l_oin_aqyy) = lptr_oin_aer(:,:,iphase)
	lptr_yyy_cwaer(:,:,l_bc_aqyy ) = lptr_bc_aer( :,:,iphase)
	lptr_yyy_cwaer(:,:,l_oc_aqyy ) = lptr_oc_aer( :,:,iphase)




	rbox(1:ltot2) = max( 0.0, rsub(1:ltot2,kpeg,mpeg) )
	rbox_sv1(1:ltot2) = rbox(1:ltot2)






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


	gas_aqfrac_box(:) = 0.0



	call interface_to_aqoperator1(   &
	    istat_aqop,   &
	    dtstepc,   &
	    rbox, gas_aqfrac_box,   &
	    rbulk_cwaer, lptr_yyy_cwaer,   &
	    co2_mixrat_in, photol_no2_in, xprescribe_ph,   &
	    iradical_in, idecomp_hmsa_hso5,   &
	    yaq_beg, yaq_end, ph_cmuaq_cur,   &
	    numgas_aqfrac, id, it, jt, kt, kpeg, mpeg, ktau, icase_in )

	ph_aq_box = ph_cmuaq_cur









	call partition_cldwtr(   &
	    rbox, fr_partit_cw,   &
	    it, jt, kt, kpeg, mpeg, icase_in )






	call distribute_bulk_changes(   &
	    rbox, rbox_sv1, fr_partit_cw,   &
	    rbulk_cwaer, lptr_yyy_cwaer,   &
	    it, jt, kt, kpeg, mpeg, icase_in )





	rsub(1:ltot2,kpeg,mpeg) = max( 0.0, rbox(1:ltot2) )





	if (msectional .lt. 1000000000) then
	    call cloudchem_apply_move_sections(   &
		rbox, rbox_sv1,   &
		it, jt, kt, kpeg, mpeg, icase_in )
	end if



	return
	end subroutine mosaic_cloudchem_1box




	subroutine interface_to_aqoperator1(   &
	    istat_aqop,   &
	    dtstepc,   &
	    rbox, gas_aqfrac_box,   &
	    rbulk_cwaer, lptr_yyy_cwaer,   &
	    co2_mixrat_in, photol_no2_in, xprescribe_ph,   &
	    iradical_in, idecomp_hmsa_hso5,   &
	    yaq_beg, yaq_end, ph_cmuaq_cur,   &
	    numgas_aqfrac, id, it, jt, kt, kpeg, mpeg, ktau, icase )

	use module_state_description, only:   &
		num_chem, param_first_scalar, p_qc,   &
		p_nh3, p_hno3, p_hcl, p_sulf, p_h2so4, p_hcho,   &
		p_ora1, p_so2, p_h2o2, p_o3, p_ho,   &
		p_ho2, p_no3, p_no, p_no2, p_hono,   &

		p_pan, p_ch3o2, p_ch3oh, p_op1,      &
		p_hcooh, p_ch3ooh

	use module_data_cmu_bulkaqchem, only:   &
	        meqn1max, naers, ngas,   &
	        na4, naa, nac, nae, nah, nahmsa, nahso5,   &
	        nan, nao, nar, nas, naw,   &
	        ng4, nga, ngc, ngch3co3h, ngch3o2, ngch3o2h, ngch3oh,   &
	        ngh2o2, nghcho, nghcooh, nghno2, ngho2,   &
	        ngn, ngno, ngno2, ngno3, ngo3, ngoh, ngpan, ngso2

	use module_cmu_bulkaqchem, only:  aqoperator1

	use module_data_mosaic_asect, only:   &
        	maxd_asize, maxd_atype,   &
		cw_phase, nsize_aer, ntype_aer,   &
		lptr_so4_aer, lptr_no3_aer, lptr_cl_aer, lptr_co3_aer,   &
		lptr_msa_aer, lptr_nh4_aer, lptr_na_aer, lptr_ca_aer,   &
		lptr_oin_aer, lptr_bc_aer, lptr_oc_aer,   &
		mw_cl_aer, mw_na_aer, mw_nh4_aer, mw_no3_aer, mw_so4_aer

	use module_data_mosaic_other, only:   &
		aboxtest_units_convert, cairclm,   &
		ktemp, l2maxd, ptotclm, rcldwtr_sub


	implicit none


	integer, intent(in) ::   &
		iradical_in, idecomp_hmsa_hso5,   &
	        numgas_aqfrac, id, it, jt, kt, kpeg, mpeg, ktau, icase
	integer, intent(inout) ::   &
		istat_aqop

	integer, intent(in) :: lptr_yyy_cwaer(maxd_asize,maxd_atype,nyyy)

	real, intent(in) ::   &
	    dtstepc, co2_mixrat_in,   &
	    photol_no2_in, xprescribe_ph

	real, intent(inout) :: ph_cmuaq_cur

        real, intent(inout), dimension( 1:l2maxd ) :: rbox

        real, intent(inout), dimension( nyyy, 2 ) :: rbulk_cwaer

        real, intent(inout), dimension( 1:numgas_aqfrac ) :: gas_aqfrac_box

        real, intent(inout), dimension( meqn1max ) :: yaq_beg, yaq_end



	integer :: i, iphase, isize, itype
	integer :: iaq, istat_fatal, istat_warn
	integer :: l, lunxx, lyyy
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






	factpatm = 1.0/1.01325e6

	factlwc = 28.966*eps*1.0e6*cairclm(kpeg)

	factphoto = 1.6


	factgas = 1.0e6


	dum = cairclm(kpeg)*1.0e12
	factaerso4   = dum*mw_so4_aer
	factaerno3   = dum*mw_no3_aer
	factaercl    = dum*mw_cl_aer
	factaernh4   = dum*mw_nh4_aer
	factaerna    = dum*mw_na_aer
	factaeroin   = dum
	factaeroc    = dum
	factaerbc    = dum










	if (aboxtest_units_convert .eq. 10) then
	    factpatm = 1.0
	    factlwc = 1.0
	    factphoto = 1.0
	    factgas = 1.0
	    factaerso4   = 1.0
	    factaerno3   = 1.0
	    factaercl    = 1.0
	    factaernh4   = 1.0
	    factaerna    = 1.0
	    factaeroin   = 1.0
	    factaeroc    = 1.0
	    factaerbc    = 1.0
	end if




	temp = rbox(ktemp)

	lwc = rcldwtr_sub(kpeg,mpeg) * factlwc
	p_atm = ptotclm(kpeg) * factpatm


	p_atm = cairclm(kpeg)*1.0e3*0.082058e0*temp

	photo_in = photol_no2_in * factphoto

	rh = 1.0
	iaq = 1

	tstep_beg_sec = 0.0
	tstep_end_sec = dtstepc


	gas(:) = 0.0

	gas(nga     ) = rbox(p_nh3   ) * factgas
	gas(ngn     ) = rbox(p_hno3  ) * factgas
	if(p_hcl > param_first_scalar )   gas(ngc ) = rbox(p_hcl   ) * factgas
	if(p_sulf > param_first_scalar )  gas(ng4 ) = rbox(p_sulf  ) * factgas
	if(p_h2so4 > param_first_scalar ) gas(ng4 ) = rbox(p_h2so4 ) * factgas


	gas(nghcho  ) = rbox(p_hcho  ) * factgas
	if(p_ora1 > param_first_scalar )  gas(nghcooh ) = rbox(p_ora1  ) * factgas
        if(p_hcooh > param_first_scalar ) gas(nghcooh ) = rbox(p_hcooh ) * factgas
	gas(ngso2   ) = rbox(p_so2   ) * factgas
	gas(ngh2o2  ) = rbox(p_h2o2  ) * factgas
	gas(ngo3    ) = rbox(p_o3    ) * factgas
	gas(ngoh    ) = rbox(p_ho    ) * factgas
	gas(ngho2   ) = rbox(p_ho2   ) * factgas
	gas(ngno3   ) = rbox(p_no3   ) * factgas

	gas(ngno    ) = rbox(p_no    ) * factgas
	gas(ngno2   ) = rbox(p_no2   ) * factgas
	gas(nghno2  ) = rbox(p_hono  ) * factgas
	gas(ngpan   ) = rbox(p_pan   ) * factgas
	gas(ngch3o2 ) = rbox(p_ch3o2 ) * factgas
	gas(ngch3oh ) = rbox(p_ch3oh ) * factgas
        if(p_op1 > param_first_scalar )    gas(ngch3o2h) = rbox(p_op1   ) * factgas
        if(p_ch3ooh > param_first_scalar ) gas(ngch3o2h) = rbox(p_ch3ooh) * factgas


	aerosol(:) = 0.0
	rbulk_cwaer(:,:) = 0.0

	iphase = cw_phase
	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)

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
		'*** mosaic_cloudchem_driver, subr interface_to_aqoperator1'
	    write(6,'(a,4i5,2i10)')   &
		'    id,it,jt,kt, istat_fatal, warn =',   &
		id, it, jt, kt, istat_fatal, istat_warn
	    istat_aqop = -10
	end if





	dum = totsulf_end - totsulf_beg
	dumb = max( totsulf_beg, totsulf_end )
	if (abs(dum) .gt. max(1.0e-3,1.0e-3*dumb)) then
	    write(6,*)   &
		'*** mosaic_cloudchem_driver, sulfur balance warning'
	    write(6,'(a,4i5,1p,3e12.4)')   &
		'    id,it,jt,kt, total_sulfur_beg, _end, _error =',   &
		id, it, jt, kt, totsulf_beg, totsulf_end, dum
	end if




	rbox(p_nh3   ) = gas(nga     ) / factgas
	rbox(p_hno3  ) = gas(ngn     ) / factgas
	if(p_hcl .gt. param_first_scalar)   rbox(p_hcl   ) = gas(ngc     ) / factgas
	if(p_sulf .gt. param_first_scalar)  rbox(p_sulf  ) = gas(ng4) / factgas
	if(p_h2so4 .gt. param_first_scalar) rbox(p_h2so4 ) = gas(ng4) / factgas


	rbox(p_hcho  ) = gas(nghcho  ) / factgas
	if(p_ora1 > param_first_scalar ) rbox(p_ora1  ) = gas(nghcooh ) / factgas
        if(p_hcooh > param_first_scalar ) rbox(p_hcooh) = gas(nghcooh ) / factgas
	rbox(p_so2   ) = gas(ngso2   ) / factgas
	rbox(p_h2o2  ) = gas(ngh2o2  ) / factgas
	rbox(p_o3    ) = gas(ngo3    ) / factgas
	rbox(p_ho    ) = gas(ngoh    ) / factgas
	rbox(p_ho2   ) = gas(ngho2   ) / factgas
	rbox(p_no3   ) = gas(ngno3   ) / factgas

	rbox(p_no    ) = gas(ngno    ) / factgas
	rbox(p_no2   ) = gas(ngno2   ) / factgas
	rbox(p_hono  ) = gas(nghno2  ) / factgas
	rbox(p_pan   ) = gas(ngpan   ) / factgas
	rbox(p_ch3o2 ) = gas(ngch3o2 ) / factgas
	rbox(p_ch3oh ) = gas(ngch3oh ) / factgas
	if(p_op1 > param_first_scalar )    rbox(p_op1   ) = gas(ngch3o2h) / factgas
        if(p_ch3ooh > param_first_scalar ) rbox(p_ch3ooh) = gas(ngch3o2h) / factgas

	gas_aqfrac_box(:) = 0.0

	if (p_nh3   .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_nh3   ) = gas_aqfrac_cmu(nga     )
	if (p_hno3  .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_hno3  ) = gas_aqfrac_cmu(ngn     )
	if (p_hcl   .le. numgas_aqfrac .and. p_hcl .gt. param_first_scalar)   &
		gas_aqfrac_box(p_hcl   ) = gas_aqfrac_cmu(ngc     )
	if (p_sulf  .le. numgas_aqfrac .and. p_sulf .gt. param_first_scalar)   &
		gas_aqfrac_box(p_sulf  ) = gas_aqfrac_cmu(ng4     )
	if (p_h2so4  .le. numgas_aqfrac .and. p_h2so4 .gt. param_first_scalar) &
		gas_aqfrac_box(p_h2so4 ) = gas_aqfrac_cmu(ng4     )

	if (p_hcho  .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_hcho  ) = gas_aqfrac_cmu(nghcho  )
	if (p_ora1  .le. numgas_aqfrac .and. p_ora1 .gt. param_first_scalar)   &
		gas_aqfrac_box(p_ora1  ) = gas_aqfrac_cmu(nghcooh )
        if (p_hcooh .le. numgas_aqfrac .and. p_hcooh .gt. param_first_scalar)   &
                gas_aqfrac_box(p_hcooh ) = gas_aqfrac_cmu(nghcooh )
	if (p_so2   .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_so2   ) = gas_aqfrac_cmu(ngso2   )
	if (p_h2o2  .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_h2o2  ) = gas_aqfrac_cmu(ngh2o2  )
	if (p_o3    .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_o3    ) = gas_aqfrac_cmu(ngo3    )
	if (p_ho    .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_ho    ) = gas_aqfrac_cmu(ngoh    )
	if (p_ho2   .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_ho2   ) = gas_aqfrac_cmu(ngho2   )
	if (p_no3   .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_no3   ) = gas_aqfrac_cmu(ngno3   )

	if (p_no    .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_no    ) = gas_aqfrac_cmu(ngno    )
	if (p_no2   .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_no2   ) = gas_aqfrac_cmu(ngno2   )
	if (p_hono  .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_hono  ) = gas_aqfrac_cmu(nghno2  )
	if (p_pan   .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_pan   ) = gas_aqfrac_cmu(ngpan   )
	if (p_ch3o2 .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_ch3o2 ) = gas_aqfrac_cmu(ngch3o2 )
	if (p_ch3oh .le. numgas_aqfrac)   &
		gas_aqfrac_box(p_ch3oh ) = gas_aqfrac_cmu(ngch3oh )
	if (p_op1   .le. numgas_aqfrac .and. p_op1 .gt. param_first_scalar)   &
		gas_aqfrac_box(p_op1   ) = gas_aqfrac_cmu(ngch3o2h)
        if (p_ch3ooh.le. numgas_aqfrac .and. p_ch3ooh .gt. param_first_scalar)   &
                gas_aqfrac_box(p_ch3ooh) = gas_aqfrac_cmu(ngch3o2h)

	rbulk_cwaer(l_so4_aqyy,2) = aerosol(na4) / factaerso4
	rbulk_cwaer(l_no3_aqyy,2) = aerosol(nan) / factaerno3
	rbulk_cwaer(l_cl_aqyy, 2) = aerosol(nac) / factaercl
	rbulk_cwaer(l_nh4_aqyy,2) = aerosol(naa) / factaernh4
	rbulk_cwaer(l_na_aqyy, 2) = aerosol(nas) / factaerna
	rbulk_cwaer(l_oin_aqyy,2) = aerosol(nar) / factaeroin
	rbulk_cwaer(l_bc_aqyy, 2) = aerosol(nae) / factaerbc
	rbulk_cwaer(l_oc_aqyy, 2) = aerosol(nao) / factaeroc





	return
	end subroutine interface_to_aqoperator1




	subroutine partition_cldwtr(   &
	    rbox, fr_partit_cw,   &
	    it, jt, kt, kpeg, mpeg, icase )

	use module_state_description, only:   &
        param_first_scalar

	use module_data_mosaic_asect, only:   &
        	maxd_asize, maxd_atype,   &
		cw_phase, nsize_aer, ntype_aer, ncomp_aer,   &
		massptr_aer, numptr_aer,   &
		dens_aer, mw_aer, volumlo_sect, volumhi_sect

	use module_data_mosaic_other, only:   &
		aboxtest_units_convert, cairclm,   &
		ktemp, l2maxd, ptotclm, rcldwtr_sub


	implicit none


	integer, intent(in) :: it, jt, kt, kpeg, mpeg, icase

        real, intent(inout), dimension( 1:l2maxd ) :: rbox

        real, intent(inout), dimension( maxd_asize, maxd_atype ) ::   &
		fr_partit_cw


	integer :: iphase, isize, itype
	integer :: jdone_mass, jdone_numb, jpos, jpos_mass, jpos_numb
	integer :: l, ll, lunxx
	integer :: p1st

	real, parameter :: partit_wght_mass = 0.5

	real :: dum, duma, dumb, dumc, dummass, dumnumb, dumvolu
	real :: tmass, tnumb, umass, unumb, wmass, wnumb
    real, dimension( maxd_asize, maxd_atype ) :: fmass, fnumb, xmass, xnumb

    p1st = PARAM_FIRST_SCALAR

	iphase = cw_phase
	tmass = 0.0
	tnumb = 0.0
	umass = 0.0
	unumb = 0.0








	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    dummass = 0.0
	    dumvolu = 0.0
	    do ll = 1, ncomp_aer(itype)
		l = massptr_aer(ll,isize,itype,iphase)
		if (l .ge. p1st) then
		    dum = max( 0.0, rbox(l) )*mw_aer(ll,itype)
		    dummass = dummass + dum
		    dumvolu = dumvolu + dum/dens_aer(ll,itype)
		end if
	    end do

	    l = numptr_aer(isize,itype,iphase)
	    dumnumb = max( 0.0, rbox(l) )
	    if (dumnumb .gt. dumvolu/volumlo_sect(isize,itype)) then
	        dumnumb = dumvolu/volumlo_sect(isize,itype)
		rbox(l) = dumnumb
	    else if (dumnumb .lt. dumvolu/volumhi_sect(isize,itype)) then
	        dumnumb = dumvolu/volumhi_sect(isize,itype)
		rbox(l) = dumnumb
	    end if

	    if (dummass .lt. 1.0e-37) dummass = 0.0
	    xmass(isize,itype) = dummass
	    if (dumnumb .lt. 1.0e-37) dumnumb = 0.0
	    xnumb(isize,itype) = dumnumb

	    tmass = tmass + xmass(isize,itype)
	    tnumb = tnumb + xnumb(isize,itype)
	    umass = max( umass, xmass(isize,itype) )
	    unumb = max( unumb, xnumb(isize,itype) )
	end do
	end do





	jdone_mass = 0
	jdone_numb = 0
	jpos_mass = 0
	jpos_numb = 0
	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
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
		if (fr_partit_cw(isize,itype) .gt. 0.0) jpos = jpos + 1
	    end do
	    end do


	    if (jpos .eq. 1) then
	    do itype = 1, ntype_aer
	    do isize = 1, nsize_aer(itype)
		if (fr_partit_cw(isize,itype) .gt. 0.0)   &
			fr_partit_cw(isize,itype) = 1.0
	    end do
	    end do
	    end if
	end if





	return
	end subroutine partition_cldwtr




	subroutine distribute_bulk_changes(   &
		rbox, rbox_sv1, fr_partit_cw,   &
		rbulk_cwaer, lptr_yyy_cwaer,   &
		it, jt, kt, kpeg, mpeg, icase )

	use module_state_description, only:   &
		param_first_scalar

	use module_data_mosaic_asect, only:   &
		maxd_asize, maxd_atype,   &
		cw_phase, nsize_aer, ntype_aer,   &
		lptr_so4_aer, lptr_no3_aer, lptr_cl_aer, lptr_co3_aer,   &
		lptr_msa_aer, lptr_nh4_aer, lptr_na_aer, lptr_ca_aer,   &
		lptr_oin_aer, lptr_bc_aer, lptr_oc_aer

	use module_data_mosaic_other, only:  l2maxd, lunout, name


	implicit none


	integer, intent(in) :: it, jt, kt, kpeg, mpeg, icase

	integer, intent(in) :: lptr_yyy_cwaer(maxd_asize,maxd_atype,nyyy)

        real, intent(inout), dimension( 1:l2maxd ) :: rbox, rbox_sv1

        real, intent(in), dimension( maxd_asize, maxd_atype ) ::   &
		fr_partit_cw

        real, intent(in), dimension( nyyy, 2 ) :: rbulk_cwaer



	integer :: iphase, isize, itype
	integer :: idone, icount, ncount
	integer :: jpos, jpos_sv
	integer :: l, lunxx, lunxxaa, lunxxbb, lyyy
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
	    if (fr_partit_cw(isize,itype) .gt. 0) jpos = jpos + 1
	end do
	end do
	jpos_sv = jpos




	if (jpos_sv .eq. 1) then
	    do lyyy = 1, nyyy

	    do itype = 1, ntype_aer
	    do isize = 1, nsize_aer(itype)
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
	end subroutine distribute_bulk_changes




	subroutine cloudchem_apply_move_sections(   &
		rbox, rbox_sv1,   &
		it, jt, kt, kpeg, mpeg, icase )

	use module_state_description, only:   &
		param_first_scalar

	use module_data_mosaic_asect, only:   &
		msectional,   &
		maxd_asize, maxd_atype,   &
		cw_phase, nsize_aer, ntype_aer, ncomp_aer,   &
		massptr_aer, numptr_aer, mw_aer, dens_aer,   &
		lptr_so4_aer, lptr_no3_aer, lptr_cl_aer, lptr_co3_aer,   &
		lptr_msa_aer, lptr_nh4_aer, lptr_na_aer, lptr_ca_aer,   &
		lptr_oin_aer, lptr_bc_aer, lptr_oc_aer,   &
		drymass_aftgrow, drymass_pregrow,   &
		drydens_aftgrow, drydens_pregrow

	use module_data_mosaic_other, only:  l2maxd, name, rsub

	use module_mosaic_movesect, only:  move_sections


	implicit none


	integer, intent(in) :: it, jt, kt, kpeg, mpeg, icase

        real, intent(inout), dimension( 1:l2maxd ) :: rbox, rbox_sv1



	integer :: idum_msect
	integer :: iphase, isize, itype
	integer :: l, ll, lunxx
	integer :: p1st
	integer :: lptr_dum(maxd_asize,maxd_atype)

	real :: densdefault
	real :: dmaft, dmpre, dvaft, dvpre
	real :: duma, dumb, dumc
	real :: smallmassbb


	p1st = param_first_scalar
	iphase = cw_phase





	densdefault = 2.0
	smallmassbb = 1.0e-30

	do 1800 itype = 1, ntype_aer
	do 1800 isize = 1, nsize_aer(itype)
	    dmaft = 0.0
	    dmpre = 0.0
	    dvaft = 0.0
	    dvpre = 0.0

	    do ll = 1, ncomp_aer(itype)
		l = massptr_aer(ll,isize,itype,iphase)
		if (l .ge. p1st) then
		    duma = mw_aer(ll,itype)
		    dmaft = dmaft + duma*rbox(l)
		    dmpre = dmpre + duma*rbox_sv1(l)

		    duma = duma/dens_aer(ll,itype)
		    dvaft = dvaft + duma*rbox(l)
		    dvpre = dvpre + duma*rbox_sv1(l)
		end if
	    end do

	    drymass_aftgrow(isize,itype) = dmaft
	    drymass_pregrow(isize,itype) = dmpre

	    if (min(dmaft,dvaft) .le. smallmassbb) then
		drydens_aftgrow(isize,itype) = densdefault
	    else
		drydens_aftgrow(isize,itype) = dmaft/dvaft
	    end if
	    if (min(dmpre,dvpre) .le. smallmassbb) then
		drydens_pregrow(isize,itype) = densdefault
	    else
		drydens_pregrow(isize,itype) = dmpre/dvpre
	    end if

1800	continue




	idum_msect = msectional



	call move_sections( 2, it, jt, kpeg, mpeg )

	msectional = idum_msect





	return
	end subroutine cloudchem_apply_move_sections




	subroutine mosaic_cloudchem_dumpaa(   &
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

	use module_configure, only:  grid_config_rec_type

	use module_data_mosaic_asect

	use module_data_mosaic_other, only:  k_pegbegin, name

	use module_mosaic_driver, only:  mapaer_tofrom_host


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
	'mosaic_cloudchem_dumpaa - ktau, i, j, k =', 4i5 )

	itype = 1
	do 2900 isize = 1, nsize_aer(itype)

	write(*,9110) isize
9110	format( / 'isize =', i3 /   &
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
		(chem(it,kt,jt,l), name(l)(1:10), l=1,num_chem)
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
	end subroutine mosaic_cloudchem_dumpaa




	end module module_mosaic_cloudchem
