








MODULE module_aer_drydep

IMPLICIT NONE


CONTAINS



	subroutine aer_drydep_driver(                                      &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		gmt, julday,                                               &
		t_phy, rho_phy, p_phy,                                     &
		alt, p8w, t8w, dz8w, z, z_at_w,                            &
		ust, aer_res, ivgtyp, vegfra, pbl, rmol, znt,              &
		moist, chem, ddvel,                                        &
		h2oai, h2oaj, numgas,                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )


































	use module_configure, only:  grid_config_rec_type, num_moist, num_chem
	use module_state_description, only:  param_first_scalar

	use module_data_mosaic_asect
	use module_data_mosaic_other

        use modal_aero_data, only:  ntot_amode_cam_mam => ntot_amode
	implicit none


	integer, intent(in) ::   &
		id, ktau, julday, aer_mech_id,   &
		numgas

	integer, intent(in) ::   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(in) :: dtstep, gmt

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		t_phy, rho_phy, p_phy, &
		alt, p8w, t8w, dz8w, z, z_at_w

	integer, intent(in),   &
		dimension( ims:ime, jms:jme ) :: &
		ivgtyp

	real, intent(in),   &
		dimension( ims:ime, jms:jme ) :: &
		ust, vegfra, pbl, rmol, znt

	real, intent(in),   &
		dimension( its:ite, jts:jte ) :: &
		aer_res

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
		moist
 
	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem

	real, intent(inout),   &
		dimension( its:ite, jts:jte, 1:num_chem ) :: &
		ddvel

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		h2oai, h2oaj

	type(grid_config_rec_type), intent(in) :: config_flags



	integer :: i, ibin, iflag_coarse, iflag_nac, &
	           iok_aer_drydep_opt, ioptaa, iphase, iseason, isurftype, itype, iwetsurf
	integer :: itmpa
	integer :: j
	integer :: n, nbin_add, numaer

	real :: airdensbox, airkinvisc
	real :: depresist_a, depresist_d0, depresist_d3
	real :: depresist_unstabpblfact
	real :: freepath
	real :: presbox
	real :: tempbox
	real :: tmpa, tmp_lnsg, tmp_num
	real :: ustar
	real :: vsettl_0, vsettl_3
	real :: wetdgnum, wetdens, wstar
	real :: wetdgnum_si, wetdens_si

	real, allocatable :: alnsg_add(:)
	real, allocatable :: depvel0_add(:), depvel3_add(:)
	real, allocatable :: wetdgnum_add(:), wetdens_add(:)
	real, allocatable :: yextra_add(:,:)

	character(len=160) :: txtaa


	write(*,'(a,2(1x,i10))') &
		'*** aer_drydep_driver - aer_mech_id, aer_drydep_opt = ', &
		aer_mech_id, config_flags%aer_drydep_opt


	if ((aer_mech_id == 1) .or.   &
	    (aer_mech_id == 2)) then
	    nbin_add = 3

	else if (aer_mech_id == 3) then
	    nbin_add = sum( nsize_aer(1:ntype_aer) )
        else if (aer_mech_id == 4) then
            nbin_add = ntot_amode_cam_mam
	else
	    write(txtaa,'(a,2(1x,i10))') &
		'*** aer_drydep_driver - bad aer_mech_id, aer_drydep_opt = ', &
		aer_mech_id, config_flags%aer_drydep_opt
	    write(*,'(a)') trim(txtaa)
	    call wrf_error_fatal3("<stdin>",166,&
txtaa )
	end if

	allocate( alnsg_add(    nbin_add ) )
	allocate( depvel0_add(  nbin_add ) )
	allocate( depvel3_add(  nbin_add ) )
	allocate( wetdgnum_add( nbin_add ) )
	allocate( wetdens_add(  nbin_add ) )
	allocate( yextra_add(   nbin_add, 10 ) )


main_j_loop:   &
	do j = jts, jte
main_i_loop:   &
	do i = its, ite



	if ((aer_mech_id == 1) .or.   &
	    (aer_mech_id == 2)) then
	    call sorgam_aer_drydep_prep(                                   &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		gmt, julday,                                               &
		t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
		moist, chem,                                               &
		h2oai, h2oaj,                                              &
		nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
		yextra_add,                                                &
		i, j,                                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )

	else if (aer_mech_id == 3) then
	    call mosaic_aer_drydep_prep(                                   &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		gmt, julday,                                               &
		t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
		moist, chem,                                               &
		nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
		yextra_add,                                                &
		i, j,                                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )
       else if (aer_mech_id == 4) then
            call cam_mam_aer_drydep_prep(                                  &
                id, ktau, dtstep, config_flags, aer_mech_id,               &
                gmt, julday,                                               &
                t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
                moist, chem,                                               &
                nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
                yextra_add,                                                &
                i, j,                                                      &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte                                  )

	else
	    txtaa = '*** aer_drydep_driver - fatal error 200'
	    call wrf_error_fatal3("<stdin>",227,&
txtaa )
	end if



	do ibin = 1, nbin_add

	ustar = ust(i,j)
	depresist_a = aer_res(i,j)
	tempbox = t_phy(i,kts,j)
	presbox = p_phy(i,kts,j)

	wetdgnum = wetdgnum_add(ibin)
	wetdens  = wetdens_add(ibin)
	tmp_lnsg = alnsg_add(ibin)

	iok_aer_drydep_opt = 0

	if ((config_flags%aer_drydep_opt >= 100) .and.   &
	    (config_flags%aer_drydep_opt <= 199)) then
	    airkinvisc = 0.0
	    freepath = 0.0
	    ustar = max( 1.e-1, ust(i,j) )   
	    tmpa = .001*p_phy(i,kts,j)
	    presbox = 1.e3*tmpa

	    airdensbox = 1.0/alt(i,kts,j)

	    iflag_coarse = 0
	    if (aer_mech_id <= 2) then
		if (ibin == 3) iflag_coarse = 1
	    else
		
		if (yextra_add(ibin,9) >= 1.0e-6) iflag_coarse = 1
	    end if

	    if (config_flags%aer_drydep_opt == 101) then
		iok_aer_drydep_opt = 1

		call sorgam_aer_drydepvel_2(   &
		    iflag_coarse,   &
		    wetdgnum, tmp_lnsg, wetdens,   &
		    tempbox, presbox, airdensbox,   &
		    ustar, depresist_a,   &
		    pbl(i,j), znt(i,j), rmol(i,j),   &
		    airkinvisc, freepath,   &
		    depvel0_add(ibin), depvel3_add(ibin),   &
		    depresist_d0, depresist_d3,   &
		    vsettl_0, vsettl_3 )

	    else if (config_flags%aer_drydep_opt == 111) then
		iok_aer_drydep_opt = 1
		wstar = 0.0

		call sorgam_aer_drydepvel_1(   &
		    iflag_coarse,   &
		    wetdgnum, tmp_lnsg, wetdens,   &
		    tempbox, presbox, airdensbox,   &
		    ustar, wstar, depresist_a,   &
		    airkinvisc, freepath,   &
		    depvel0_add(ibin), depvel3_add(ibin),   &
		    depresist_d0, depresist_d3,   &
		    vsettl_0, vsettl_3 )

	    end if
	
	    depresist_unstabpblfact = 0.0



	else if ((config_flags%aer_drydep_opt >= 200) .and.   &
	         (config_flags%aer_drydep_opt <= 299)) then
	    airkinvisc = 0.0
	    freepath = 0.0
	    depresist_unstabpblfact = 0.0

	    if ((config_flags%aer_drydep_opt == 201) .or.   &
	        (config_flags%aer_drydep_opt == 211)) then
		iok_aer_drydep_opt = 1
		airdensbox = rho_phy(i,kts,j)
		ioptaa = 1
		if (config_flags%aer_drydep_opt == 211) ioptaa = 2

		call mosaic_aer_drydepvel_1(   &
		    ioptaa, wetdgnum, tmp_lnsg, wetdens,   &
		    tempbox, airdensbox, ustar, depresist_a,   &
		    airkinvisc, freepath, depresist_unstabpblfact,   &
		    depvel0_add(ibin), depvel3_add(ibin),   &
		    depresist_d0, depresist_d3,   &
		    vsettl_0, vsettl_3 )

	    end if
	
	else if ((config_flags%aer_drydep_opt >= 300) .and.   &
	         (config_flags%aer_drydep_opt <= 399)) then
	    airkinvisc = 0.0
	    freepath = 0.0
	    depresist_unstabpblfact = 0.0

	    if ((config_flags%aer_drydep_opt == 301) .or. &
	        (config_flags%aer_drydep_opt == 302) .or. &
	        (config_flags%aer_drydep_opt == 311) .or. &
	        (config_flags%aer_drydep_opt == 312)) then
		iok_aer_drydep_opt = 1
		airdensbox = rho_phy(i,kts,j)
		isurftype = ivgtyp(i,j)
		
		iwetsurf = 0
		
		
		iseason = 1
		if (julday < 90 .or. julday > 270) iseason = 2
		if ((config_flags%aer_drydep_opt == 302) .or. &
		    (config_flags%aer_drydep_opt == 312)) iseason = 1
		ioptaa = 1
		if ((config_flags%aer_drydep_opt == 311) .or. &
		    (config_flags%aer_drydep_opt == 312)) ioptaa = 2

		call zhang2001_aer_drydepvel_1(   &
		    wetdgnum, tmp_lnsg, wetdens,   &
		    tempbox, airdensbox, ustar, depresist_a,   &
		    ioptaa, iseason, isurftype, 2, iwetsurf,   &
		    airkinvisc, freepath,   &
		    1, 1,   &
		    depvel0_add(ibin), depvel3_add(ibin),   &
		    depresist_d0, depresist_d3,   &
		    vsettl_0, vsettl_3 )

	    end if
	
	end if

	if (iok_aer_drydep_opt <= 0) then
	    write(txtaa,'(a,2(1x,i10))') &
		'*** aer_drydep_driver - bad aer_mech_id, aer_drydep_opt = ', &
		aer_mech_id, config_flags%aer_drydep_opt
	    write(*,'(a)') trim(txtaa)
	    call wrf_error_fatal3("<stdin>",365,&
txtaa )
	end if


	end do   



	if ((aer_mech_id == 1) .or.   &
	    (aer_mech_id == 2)) then
	    call sorgam_aer_drydep_load_ddvel(                             &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		numgas, ddvel,                                             &
		nbin_add, depvel0_add, depvel3_add,                        &
		i, j,                                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )

	else if (aer_mech_id == 3) then
	    call mosaic_aer_drydep_load_ddvel(                             &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		ddvel,                                                     &
		nbin_add, depvel0_add, depvel3_add,                        &
		i, j,                                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )
        else if (aer_mech_id == 4) then
            call cam_mam_aer_drydep_load_ddvel(                            &
                id, ktau, dtstep, config_flags, aer_mech_id,               &
                ddvel,                                                     &
                nbin_add, depvel0_add, depvel3_add,                        &
                i, j,                                                      &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte                                  )

	else
	    txtaa = '*** aer_drydep_driver - fatal error 400'
	    call wrf_error_fatal3("<stdin>",406,&
txtaa )
	end if


	end do main_i_loop
	end do main_j_loop



	deallocate( alnsg_add    )
	deallocate( depvel0_add  )
	deallocate( depvel3_add  )
	deallocate( wetdgnum_add )
	deallocate( wetdens_add  )
	deallocate( yextra_add   )


80000	continue
	return
	end subroutine aer_drydep_driver




	subroutine mosaic_aer_drydep_prep(                                 &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		gmt, julday,                                               &
		t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
		moist, chem,                                               &
		nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
		yextra_add,                                                &
		i, j,                                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )

	use module_configure, only:  num_moist, num_chem, &
		grid_config_rec_type
	use module_state_description, only:  param_first_scalar

	use module_data_mosaic_asect
	use module_data_mosaic_other, only:  lunerr, lunout, pi

	implicit none


	integer, intent(in) ::   &
		id, ktau, julday, aer_mech_id, nbin_add,   &
		i, j

	integer, intent(in) ::   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(in) :: dtstep, gmt

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		t_phy, rho_phy, p_phy, alt, p8w, t8w

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
		moist

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem

	real, intent(inout),   &
		dimension( 1:nbin_add ) :: &
		alnsg_add, wetdgnum_add, wetdens_add

	real, intent(inout),   &
		dimension( 1:nbin_add, 10 ) :: &
		yextra_add

	type(grid_config_rec_type), intent(in) :: config_flags



	integer :: ibin, ijcount, iok, iphase, itype
	integer :: itmpa
	integer :: k
	integer :: l, ll, l1, n
	integer :: p1st

	real, parameter :: densdefault = 2.0
	real, parameter :: smallmassaa = 1.0e-20
	real, parameter :: smallmassbb = 1.0e-30
	real, parameter :: piover6 = pi/6.0
	real, parameter :: onethird = 1.0/3.0

	real :: drydens, drydp, drymass, dryvol
	real :: rnum
	real :: tmpa, tmp_lnsg
	real :: wetdgnum, wetdens, wetdp, wetmass, wetvol
	real :: wetdgnum_si, wetdens_si

	character(len=160) :: txtaa


	n = sum( nsize_aer(1:ntype_aer) )
	if (n /= nbin_add) then
	    write(txtaa,'(a,2(1x,i10))') &
		'*** mosaic_aer_drydep_prep - bad nbin_add, sum(nsize)', &
		nbin_add, n
	    write(*,'(a)') trim(txtaa)
	    call wrf_error_fatal3("<stdin>",515,&
txtaa )
	    stop
	end if



	p1st = param_first_scalar
	lunerr = -1
	lunout = -1

	ijcount = 0

	k = kts

	ijcount = ijcount + 1


	iphase = ai_phase


	ibin = 0
	do 2900 itype = 1, ntype_aer
	do 2800 n = 1, nsize_aer(itype)
	ibin = ibin + 1


	dryvol = 0.0
	drymass = 0.0
	do ll = 1, ncomp_aer(itype)
	    l1 = massptr_aer(ll,n,itype,iphase)
	    tmpa = max( chem(i,k,j,l1), 0.0 )
	    drymass = drymass + tmpa
	    dryvol = dryvol + tmpa/dens_aer(ll,itype)
	end do

	l1 = waterptr_aer(n,itype)
	tmpa = max( chem(i,k,j,l1), 0.0 )
	wetmass = drymass + tmpa
	wetvol = dryvol + tmpa/dens_water_aer

	l1 = numptr_aer(n,itype,iphase)
	rnum = max( chem(i,k,j,l1), 0.0 )

	drymass = drymass*28.966e-9 
	wetmass = wetmass*28.966e-9
	dryvol  = dryvol *28.966e-9 
	wetvol  = wetvol *28.966e-9
	rnum    = rnum*28.966e-3    

        if (drymass <= smallmassbb) then
            drydp = dcen_sect(n,itype)
            drydens = densdefault
            wetdp = drydp
            wetdens = drydens

        else
            if (drymass <= smallmassaa) then
                wetmass = drymass
                wetvol = dryvol
            end if
            drydens = drymass/dryvol
            wetdens = wetmass/wetvol

            if (rnum >= dryvol/volumlo_sect(n,itype)) then
                drydp = dlo_sect(n,itype)
            else if (rnum <= dryvol/volumhi_sect(n,itype)) then
                drydp = dhi_sect(n,itype)
            else
                drydp = (dryvol/(piover6*rnum))**onethird
            end if

            if (abs(wetvol) > (1000.*abs(dryvol))) then
              tmpa=10.0
            else
              tmpa=abs(wetvol/dryvol)**onethird
              tmpa=max(1.0,min(tmpa,10.0))
            endif
            wetdp = drydp*tmpa

        end if




	tmp_lnsg = log( 1.0 )
	wetdgnum = wetdp * exp( -1.5*tmp_lnsg*tmp_lnsg )

	wetdgnum_add(ibin) = wetdgnum * 1.0e-2 
	wetdens_add(ibin) = wetdens * 1.0e3    
	alnsg_add(ibin) = tmp_lnsg


	yextra_add(ibin,1) = drymass/28.966e-9 
	yextra_add(ibin,5) = wetmass/28.966e-9
	yextra_add(ibin,2) = dryvol/28.966e-6  
	yextra_add(ibin,6) = wetvol/28.966e-6
	yextra_add(ibin,3) = drydens*1.0e3     
	yextra_add(ibin,7) = wetdens*1.0e3
	yextra_add(ibin,4) = drydp*1.0e-2      
	yextra_add(ibin,8) = wetdp*1.0e-2

	yextra_add(ibin,9) = dcen_sect(n,itype)*1.0e-2   
	yextra_add(ibin,10) = rnum/28.966e-3    


2800	continue
2900	continue


	return
	end subroutine mosaic_aer_drydep_prep




	subroutine mosaic_aer_drydep_load_ddvel(                           &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		ddvel,                                                     &
		nbin_add, depvel0_add, depvel3_add,                        &
		i, j,                                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )

	use module_configure, only:  num_moist, num_chem, &
		grid_config_rec_type
	use module_state_description, only:  param_first_scalar

	use module_data_mosaic_asect

	implicit none


	integer, intent(in) ::   &
		id, ktau, aer_mech_id, nbin_add,   &
		i, j

	integer, intent(in) ::   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(in) :: dtstep

	real, intent(in),   &
		dimension( 1:nbin_add ) :: &
		depvel0_add, depvel3_add

	real, intent(out),   &
		dimension( its:ite, jts:jte, 1:num_chem ) :: &
		ddvel

	type(grid_config_rec_type), intent(in) :: config_flags



	integer :: ibin, iphase, itype
	integer :: l, ll, n

	character(len=160) :: txtaa


	n = sum( nsize_aer(1:ntype_aer) )
	if (n /= nbin_add) then
	    write(txtaa,'(a,2(1x,i10))') &
		'*** mosaic_aer_drydep_load_ddvel - bad nbin_add, sum(nsize)', &
		nbin_add, n
	    write(*,'(a)') trim(txtaa)
	    call wrf_error_fatal3("<stdin>",684,&
txtaa )
	    stop
	end if



	do iphase = 1, nphase_aer
	ibin = 0
	do itype = 1, ntype_aer
	do n = 1, nsize_aer(itype)
	    ibin = ibin + 1
	    do ll = -2, ncomp_plustracer_aer(itype)
		if (ll == -2) then
		    l = numptr_aer(n,itype,iphase)
		else if (ll == -1) then
		    l = -1
		    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
		else if (ll == 0) then
		    l = -1
		    if (iphase == ai_phase) l = hyswptr_aer(n,itype)
		else
		    l = massptr_aer(ll,n,itype,iphase)
		end if
		if (iphase > 1) l = -1

		if ((l >= param_first_scalar) .and. (l <= num_chem)) then
		    ddvel(i,j,l) = depvel3_add(ibin)
		end if
	    end do
	end do
	end do
	end do


	return
	end subroutine mosaic_aer_drydep_load_ddvel




	subroutine mosaic_aer_drydepvel_1(   &
	    ioptaa, dgnum, alnsg, aerodens,   &
	    temp, airdens, ustar, depresist_a,   &
	    airkinvisc, freepath, depresist_unstabpblfact,   &
	    depvel_0, depvel_3,   &
	    depresist_d0, depresist_d3,   &
	    vsettl_0, vsettl_3 )






























	implicit none

	integer, intent(in) ::   &
	    ioptaa
	real, intent(in) ::   &
	    dgnum, alnsg, aerodens,   &
	    temp, airdens, ustar, depresist_a
	real, intent(inout) ::   &
	    airkinvisc, freepath, depresist_unstabpblfact
	real, intent(out) ::   &
	    depvel_0, depvel_3,   &
	    depresist_d0, depresist_d3,   &
	    vsettl_0, vsettl_3

	real :: aerodiffus_0, schmidt_0, stokes_0, facdepresist_d0
	real :: aerodiffus_3, schmidt_3, stokes_3, facdepresist_d3
	common / aerosol_depvel_cmn01 /   &
      		aerodiffus_0, schmidt_0, stokes_0, facdepresist_d0,   &
      		aerodiffus_3, schmidt_3, stokes_3, facdepresist_d3

	real :: aerodiffus_dgnum, alnsg2,   &
	    stickfrac,   &
	    tmpa, tmpb, tmp_dvm,   &
	    vsettl_dgnum, xknudsen, xknudsenfact

	real, parameter :: pi = 3.1415926536

	real, parameter :: gravity = 9.80616

	real, parameter :: boltzmann = 1.3807e-23


	if (airkinvisc <= 0.0) then

	    airkinvisc = ( 1.8325e-5 * (416.16/(temp+120.0)) *   &
				((temp/296.16)**1.5) ) / airdens
	end if
	if (freepath <= 0.0) then

	    freepath = 7.39758e-2 * airkinvisc / sqrt(temp)
	end if
	if (depresist_unstabpblfact <= 0) then

	    depresist_unstabpblfact = 1.0
	end if

	xknudsen = 2.*freepath/dgnum
	xknudsenfact = xknudsen*1.246
	alnsg2 = alnsg*alnsg

	vsettl_dgnum = (gravity*aerodens*dgnum*dgnum)/   &
      		(18.*airkinvisc*airdens)
	vsettl_0 = vsettl_dgnum *   &
      		( exp(2.*alnsg2) + xknudsenfact*exp(0.5*alnsg2) )
	vsettl_3 = vsettl_dgnum *   &
      		( exp(8.*alnsg2) + xknudsenfact*exp(3.5*alnsg2) )

	aerodiffus_dgnum = (boltzmann*temp)/   &
      		(3.*pi*airkinvisc*airdens*dgnum)
	aerodiffus_0 = aerodiffus_dgnum *   &
      		( exp(+0.5*alnsg2) + xknudsenfact*exp(+2.*alnsg2) )
	aerodiffus_3 = aerodiffus_dgnum *   &
      		( exp(-2.5*alnsg2) + xknudsenfact*exp(-4.*alnsg2) )

	schmidt_0 = airkinvisc/aerodiffus_0
	schmidt_3 = airkinvisc/aerodiffus_3

	stokes_0 = ustar*ustar*vsettl_0/(gravity*airkinvisc)
	stokes_3 = ustar*ustar*vsettl_3/(gravity*airkinvisc)
	
	tmp_dvm = dgnum * exp(1.5*alnsg2) 

	tmpa = (schmidt_0**(-0.66666666)) +   &
      		(10.**(-3./max(0.03,stokes_0)))

	tmpb = tmpa*ustar*depresist_unstabpblfact
	if (ioptaa == 2) then
	    if (tmp_dvm >= 5.0e-6) then
		stickfrac = exp( -2.0*sqrt(max(0.0,stokes_0)) )
		tmpb = tmpb*stickfrac
	    end if
	end if
	depresist_d0 = 1./max( tmpb, 1.e-22 )
	facdepresist_d0 = tmpa

	tmpa = (schmidt_3**(-0.66666666)) +   &
      		(10.**(-3./max(0.03,stokes_3)))

	tmpb = tmpa*ustar*depresist_unstabpblfact
	if (ioptaa == 2) then
	    if (tmp_dvm >= 5.0e-6) then
		stickfrac = exp( -2.0*sqrt(max(0.0,stokes_3)) )
		tmpb = tmpb*stickfrac
	    end if
	end if
	depresist_d3 = 1./max( tmpb, 1.e-22 )
	facdepresist_d3 = tmpa


	tmpa = depresist_a + depresist_d3 +   &
      		depresist_a*depresist_d3*vsettl_3
	depvel_3 = vsettl_3 + (1./tmpa)

	tmpa = depresist_a + depresist_d0 +   &
      		depresist_a*depresist_d0*vsettl_0
	depvel_0 = vsettl_0 + (1./tmpa)

	return
	end subroutine mosaic_aer_drydepvel_1




	subroutine sorgam_aer_drydep_prep(                                 &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		gmt, julday,                                               &
		t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
		moist, chem,                                               &
		h2oai, h2oaj,                                              &
		nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
		yextra_add,                                                &
		i, j,                                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )

	use module_configure, only:  num_moist, num_chem, &
		grid_config_rec_type
	use module_state_description

	use module_data_sorgam
	use module_aerosols_sorgam, only:  modpar

	implicit none


	integer, intent(in) ::   &
		id, ktau, julday, aer_mech_id,   &
		nbin_add, i, j

	integer, intent(in) ::   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(in) :: dtstep, gmt

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		t_phy, rho_phy, p_phy, alt, p8w, t8w

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
		moist

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		h2oai, h2oaj

	real, intent(inout),   &
		dimension( 1:nbin_add ) :: &
		alnsg_add, wetdgnum_add, wetdens_add

	real, intent(inout),   &
		dimension( 1:nbin_add, 10 ) :: &
		yextra_add

	type(grid_config_rec_type), intent(in) :: config_flags







      integer, parameter :: blksize=1



      integer, parameter :: nspcsda=l1ae  

      integer, parameter :: numcells=1
      
      


      real cblk(blksize,nspcsda) 
                                   

      real blkprs(blksize)         
      real blkta(blksize)          
      real blkdens(blksize)        



      real xlm( blksize )           
      real amu( blksize )           
      

      real vsed( blksize, naspcssed) 
      real vdep( blksize, naspcsdep) 


      real dgnuc( blksize )         
      real dgacc( blksize )         
      real dgcor( blksize )         
      


      real pmassn( blksize )        
      real pmassa( blksize )        
      real pmassc( blksize )        


      real pdensn( blksize )        
      real pdensa( blksize )        
      real pdensc( blksize )        


      real knnuc ( blksize )        
      real knacc ( blksize )        
      real kncor ( blksize )        

      integer :: k,l

      real :: convfac2
      real, dimension( ims:ime, jms:jme ) :: aer_res, ust   
      real, dimension( kts:kte ) :: p   

      character(len=160) :: txtaa


      if (3 /= nbin_add) then
          write(txtaa,'(a,2(1x,i10))') &
            '*** sorgam_aer_drydep_prep - bad nbin_add', &
            nbin_add
          write(*,'(a)') trim(txtaa)
          call wrf_error_fatal3("<stdin>",1003,&
txtaa )
          stop
      end if




      cblk=epsilc









      k=kts

      p(k)  = .001*p_phy(i,k,j)
      convfac2 = 1./alt(i,k,j)
      blkdens(blksize) = convfac2
      blkta(blksize)   = t_phy(i,k,j)
      blkprs(blksize)  = 1.e3*p(k)


      cblk(1,vso4aj   ) =   max( epsilc, chem(i,k,j,p_so4aj)*convfac2 )
      cblk(1,vso4ai   ) =   max( epsilc, chem(i,k,j,p_so4ai)*convfac2 )
      cblk(1,vnh4aj   ) =   max( epsilc, chem(i,k,j,p_nh4aj)*convfac2 )
      cblk(1,vnh4ai   ) =   max( epsilc, chem(i,k,j,p_nh4ai)*convfac2 )
      cblk(1,vno3aj   ) =   max( epsilc, chem(i,k,j,p_no3aj)*convfac2 )
      cblk(1,vno3ai   ) =   max( epsilc, chem(i,k,j,p_no3ai)*convfac2 )
      if (p_naaj >= param_first_scalar) &
         cblk(1,vnaaj ) =   max( epsilc, chem(i,k,j,p_naaj)*convfac2 )
      if (p_naai >= param_first_scalar) &
         cblk(1,vnaai ) =   max( epsilc, chem(i,k,j,p_naai)*convfac2 )
      if (p_claj >= param_first_scalar) &
         cblk(1,vclaj ) =   max( epsilc, chem(i,k,j,p_claj)*convfac2 )
      if (p_clai >= param_first_scalar) &
         cblk(1,vclai ) =   max( epsilc, chem(i,k,j,p_clai)*convfac2 )
      cblk(1,vorgaro1j) =   max( epsilc, chem(i,k,j,p_orgaro1j)*convfac2 )
      cblk(1,vorgaro1i) =   max( epsilc, chem(i,k,j,p_orgaro1i)*convfac2 )
      cblk(1,vorgaro2j) =   max( epsilc, chem(i,k,j,p_orgaro2j)*convfac2 )
      cblk(1,vorgaro2i) =   max( epsilc, chem(i,k,j,p_orgaro2i)*convfac2 )
      cblk(1,vorgalk1j) =   max( epsilc, chem(i,k,j,p_orgalk1j)*convfac2 )
      cblk(1,vorgalk1i) =   max( epsilc, chem(i,k,j,p_orgalk1i)*convfac2 )
      cblk(1,vorgole1j) =   max( epsilc, chem(i,k,j,p_orgole1j)*convfac2 )
      cblk(1,vorgole1i) =   max( epsilc, chem(i,k,j,p_orgole1i)*convfac2 )
      cblk(1,vorgba1j ) =   max( epsilc, chem(i,k,j,p_orgba1j)*convfac2 )
      cblk(1,vorgba1i ) =   max( epsilc, chem(i,k,j,p_orgba1i)*convfac2 )
      cblk(1,vorgba2j ) =   max( epsilc, chem(i,k,j,p_orgba2j)*convfac2 )
      cblk(1,vorgba2i ) =   max( epsilc, chem(i,k,j,p_orgba2i)*convfac2 )
      cblk(1,vorgba3j ) =   max( epsilc, chem(i,k,j,p_orgba3j)*convfac2 )
      cblk(1,vorgba3i ) =   max( epsilc, chem(i,k,j,p_orgba3i)*convfac2 )
      cblk(1,vorgba4j ) =   max( epsilc, chem(i,k,j,p_orgba4j)*convfac2 )
      cblk(1,vorgba4i ) =   max( epsilc, chem(i,k,j,p_orgba4i)*convfac2 )
      cblk(1,vorgpaj  ) =   max( epsilc, chem(i,k,j,p_orgpaj)*convfac2 )
      cblk(1,vorgpai  ) =   max( epsilc, chem(i,k,j,p_orgpai)*convfac2 )
      cblk(1,vecj     ) =   max( epsilc, chem(i,k,j,p_ecj)*convfac2 )
      cblk(1,veci     ) =   max( epsilc, chem(i,k,j,p_eci)*convfac2 )
      cblk(1,vp25aj   ) =   max( epsilc, chem(i,k,j,p_p25j)*convfac2 )
      cblk(1,vp25ai   ) =   max( epsilc, chem(i,k,j,p_p25i)*convfac2 )
      cblk(1,vantha   ) =   max( epsilc, chem(i,k,j,p_antha)*convfac2 )
      cblk(1,vseas    ) =   max( epsilc, chem(i,k,j,p_seas)*convfac2 )
      cblk(1,vsoila   ) =   max( epsilc, chem(i,k,j,p_soila)*convfac2 )
      cblk(1,vnu0     ) =   max( epsilc, chem(i,k,j,p_nu0)*convfac2 )
      cblk(1,vac0     ) =   max( epsilc, chem(i,k,j,p_ac0)*convfac2 )
      cblk(1,vcorn    ) =   max( epsilc, chem(i,k,j,p_corn)*convfac2 )

      cblk(1,vh2oaj   ) =   h2oaj(i,k,j)
      cblk(1,vh2oai   ) =   h2oai(i,k,j)


      call modpar(  blksize, nspcsda, numcells, &
           cblk,                                &
           blkta, blkprs,                       &
           pmassn, pmassa, pmassc,              &
           pdensn, pdensa, pdensc,              &
           xlm, amu,                            &
           dgnuc, dgacc, dgcor,                 &
           knnuc, knacc, kncor                  )                                   

      wetdgnum_add(1) = dgnuc(1)
      wetdgnum_add(2) = dgacc(1)
      wetdgnum_add(3) = dgcor(1)

      wetdens_add(1) = pdensn(1)
      wetdens_add(2) = pdensa(1)
      wetdens_add(3) = pdensc(1)

      alnsg_add(1) = xxlsgn
      alnsg_add(2) = xxlsga
      alnsg_add(3) = xxlsgc

      yextra_add(1:3,1:8) = 0.0

      yextra_add(1,5) = pmassn(1)/blkdens(1)   
      yextra_add(2,5) = pmassa(1)/blkdens(1)
      yextra_add(3,5) = pmassc(1)/blkdens(1)

      yextra_add(1,6) = cblk(1,vnu3 )*(pirs/6.0)*1.0e9/blkdens(1)   
      yextra_add(2,6) = cblk(1,vac3 )*(pirs/6.0)*1.0e9/blkdens(1)
      yextra_add(3,6) = cblk(1,vcor3)*(pirs/6.0)*1.0e9/blkdens(1)

      yextra_add(1,7) = pdensn(1)   
      yextra_add(2,7) = pdensa(1)
      yextra_add(3,7) = pdensc(1)

      yextra_add(1,8) = dgnuc(1)*exp(1.5*(alnsg_add(1)**2))   
      yextra_add(2,8) = dgacc(1)*exp(1.5*(alnsg_add(2)**2))
      yextra_add(3,8) = dgcor(1)*exp(1.5*(alnsg_add(3)**2))

      yextra_add(1,9) = dginin
      yextra_add(2,9) = dginia
      yextra_add(3,9) = dginic

      yextra_add(1,10) = cblk(1,vnu0 )/blkdens(1)   
      yextra_add(2,10) = cblk(1,vac0 )/blkdens(1)   
      yextra_add(3,10) = cblk(1,vcorn)/blkdens(1)   

	return
	end subroutine sorgam_aer_drydep_prep




	subroutine sorgam_aer_drydep_load_ddvel(                           &
		id, ktau, dtstep, config_flags, aer_mech_id,               &
		numgas, ddvel,                                             &
		nbin_add, depvel0_add, depvel3_add,                        &
		i, j,                                                      &
		ids,ide, jds,jde, kds,kde,                                 &
		ims,ime, jms,jme, kms,kme,                                 &
		its,ite, jts,jte, kts,kte                                  )

	use module_configure, only:  num_moist, num_chem, &
		grid_config_rec_type
	use module_state_description, only:  param_first_scalar, &
		p_naai, p_naaj, p_clai, p_claj

	use module_data_sorgam

	implicit none


	integer, intent(in) ::   &
		id, ktau, aer_mech_id, &
		numgas, nbin_add, i, j

	integer, intent(in) ::   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(in) :: dtstep

	real, intent(in),   &
		dimension( 1:nbin_add ) :: &
		depvel0_add, depvel3_add

	real, intent(out),   &
		dimension( its:ite, jts:jte, 1:num_chem ) :: &
		ddvel

	type(grid_config_rec_type), intent(in) :: config_flags



	integer :: l

	real :: vgsa(num_chem)

	character(len=160) :: txtaa


	if (3 /= nbin_add) then
	    write(txtaa,'(a,2(1x,i10))') &
		'*** sorgam_aer_drydep_load_ddvel - bad nbin_add', &
		nbin_add
	    write(*,'(a)') trim(txtaa)
	    call wrf_error_fatal3("<stdin>",1184,&
txtaa )
	    stop
	end if


        vgsa( : )  =  0.0



        vgsa( vso4aj   )  =  depvel3_add( 2 )  
        vgsa( vso4ai   )  =  depvel3_add( 1 )  

        vgsa( vnh4aj   )  =  vgsa( vso4aj )
        vgsa( vnh4ai   )  =  vgsa( vso4ai )
        vgsa( vno3aj   )  =  vgsa( vso4aj )
        vgsa( vno3ai   )  =  vgsa( vso4ai )
        if (p_naaj >= param_first_scalar) vgsa( vnaaj   )  =  vgsa( vso4aj )
        if (p_naai >= param_first_scalar) vgsa( vnaai   )  =  vgsa( vso4ai )
        if (p_claj >= param_first_scalar) vgsa( vclaj   )  =  vgsa( vso4aj )
        if (p_clai >= param_first_scalar) vgsa( vclai   )  =  vgsa( vso4ai )
        vgsa( vorgaro1j)  =  vgsa( vso4aj )
        vgsa( vorgaro1i)  =  vgsa( vso4ai )
        vgsa( vorgaro2j)  =  vgsa( vso4aj )
        vgsa( vorgaro2i)  =  vgsa( vso4ai )
        vgsa( vorgalk1j)  =  vgsa( vso4aj )
        vgsa( vorgalk1i)  =  vgsa( vso4ai )
        vgsa( vorgole1j)  =  vgsa( vso4aj )
        vgsa( vorgole1i)  =  vgsa( vso4ai )
        vgsa( vorgba1j )  =  vgsa( vso4aj )
        vgsa( vorgba1i )  =  vgsa( vso4ai )
        vgsa( vorgba2j )  =  vgsa( vso4aj )
        vgsa( vorgba2i )  =  vgsa( vso4ai )
        vgsa( vorgba3j )  =  vgsa( vso4aj )
        vgsa( vorgba3i )  =  vgsa( vso4ai )
        vgsa( vorgba4j )  =  vgsa( vso4aj )
        vgsa( vorgba4i )  =  vgsa( vso4ai )
        vgsa( vorgpaj  )  =  vgsa( vso4aj )
        vgsa( vorgpai  )  =  vgsa( vso4ai )
        vgsa( vecj     )  =  vgsa( vso4aj )
        vgsa( veci     )  =  vgsa( vso4ai )
        vgsa( vp25aj   )  =  vgsa( vso4aj )
        vgsa( vp25ai   )  =  vgsa( vso4ai )


        vgsa( vantha   )  =  depvel3_add( 3 )
        vgsa( vseas    )  =  vgsa( vantha )
        vgsa( vsoila   )  =  vgsa( vantha )




        vgsa( vnu0     )  =  depvel0_add( 1 )
        vgsa( vac0     )  =  depvel0_add( 2 )
        vgsa( vcorn    )  =  depvel0_add( 3 )

        do l = 1, (num_chem - numgas )
           ddvel(i, j, l+numgas ) =  vgsa( l )
        end do


	return
	end subroutine sorgam_aer_drydep_load_ddvel




	subroutine sorgam_aer_drydepvel_1(   &
	    iflag_coarse,   &
	    dgnum, alnsg, aerodens,   &
	    temp, pres, airdens, ustar, wstar, depresist_a,   &
	    airkinvisc, freepath,   &
	    depvel_0, depvel_3,   &
	    depresist_d0, depresist_d3,   &
	    vsettl_0, vsettl_3 )




























	implicit none

	integer, intent(in) ::   &
	    iflag_coarse
	real, intent(in) ::   &
	    dgnum, alnsg, aerodens,   &
	    temp, pres, airdens, ustar, wstar, depresist_a
	real, intent(inout) ::   &
	    airkinvisc, freepath
	real, intent(out) ::   &
	    depvel_0, depvel_3,   &
	    depresist_d0, depresist_d3,   &
	    vsettl_0, vsettl_3




        real, parameter :: avo=6.0221367e23

        real, parameter :: rgasuniv=8.314510

        real, parameter :: boltz=rgasuniv/avo

        real, parameter :: grav=9.80622

        real, parameter :: pss0=101325.0

        real, parameter :: tss0=288.15

        real, parameter :: two3=2.0/3.0

        real*8, parameter :: pirs=3.14159265358979324

        real, parameter :: threepi=3.0*pirs


	real, parameter :: bhat = 1.246   
	real :: alnsgy2, amu
	real :: blkta, blkprs, blkdens
	real :: dconst1, dconst1y, dconst2, dconst3y
	real :: dchat0y, dchat3y
	real :: dgyyy
	real :: ey1, esy04, esy08, esy16, esy20, esy28, esy32, esy64
	real :: esym20, esym32
	real :: knyyy
	real :: nu
	real :: pdensy
	real :: ra, rd0y, rd3y
	real :: sc0y, sc3y, st0y, st3y
	real :: ustfac, utscale
	real :: vghat0y, vghat3y
	real :: xlm
	real :: yextra_vdvg(100,3)


	blkta = temp
	blkprs = pres
	blkdens = airdens


	if (airkinvisc <= 0.0) then
	    amu = 1.458e-6*blkta*sqrt(blkta)/(blkta+110.4)
	    airkinvisc = amu/airdens
	else
	    amu = airkinvisc*airdens
	end if


	if (freepath <= 0.0) then
	    freepath = 6.6328e-8*pss0*blkta/(tss0*blkprs)
	end if
	xlm = freepath

	ra = depresist_a

	dgyyy = dgnum
	pdensy = aerodens

        alnsgy2 = alnsg**2
        ey1 = exp(0.125*alnsgy2)
        esy04 = ey1**4
        esy08 = esy04*esy04

        esy16 = esy08*esy08
        esy20 = esy16*esy04

        esy28 = esy20*esy08
        esy32 = esy16*esy16

        esy64 = esy32*esy32
        esym20 = 1.0/esy20
        esym32 = 1.0/esy32

	knyyy = 2.0*xlm/dgyyy

	dconst1 = boltz * blkta / ( threepi * amu )
	dconst1y = dconst1 / dgyyy

	dconst2 = grav / ( 18.0 * amu )
	dconst3y = dconst2 * pdensy * dgyyy**2

	dchat0y = dconst1y * ( esy04  + bhat * knyyy * esy16 )
	dchat3y = dconst1y * ( esym20 + bhat * knyyy * esym32 )
	vghat0y = dconst3y * ( esy16  + bhat * knyyy * esy04 )
	vghat3y = dconst3y * ( esy64  + bhat * knyyy * esy28 )

	nu = amu / blkdens 
	ustfac = ustar * ustar / ( grav * nu)
	utscale = ustar + 0.24 * wstar * wstar / ustar

	sc0y = nu / dchat0y      
	if (iflag_coarse > 0) then
	    rd0y = 1.0 / ( utscale * ( sc0y**(-two3) ) ) 
	else
	    st0y = max( vghat0y * ustfac , 0.01)
	    rd0y = 1.0 / ( utscale * ( sc0y**(-two3) + 10.0**(-3.0 / st0y) ) ) 
	end if
      
	depvel_0 = vghat0y + 1.0 / ( ra + rd0y + rd0y * ra * vghat0y )
	vsettl_0 = vghat0y 
	depresist_d0 = rd0y


	sc3y = nu / dchat3y      
	if (iflag_coarse > 0) then
	    rd3y = 1.0 / ( utscale * ( sc3y**(-two3) ) ) 
	else
	    st3y = max( vghat3y * ustfac , 0.01)
	    rd3y = 1.0 / ( utscale * ( sc3y**(-two3) + 10.0**(-3.0 / st3y) ) ) 
	end if
      
	depvel_3 = vghat3y + 1.0 / ( ra + rd3y + rd3y * ra * vghat3y )
	vsettl_3 = vghat3y 
	depresist_d3 = rd3y


	return
	end subroutine sorgam_aer_drydepvel_1




	subroutine sorgam_aer_drydepvel_2(   &
	    iflag_coarse,   &
	    dgnum, alnsg, aerodens,   &
	    temp, pres, airdens, ustar, depresist_a,   &
	    pblh, zntt, rmolm,   &
	    airkinvisc, freepath,   &
	    depvel_0, depvel_3,   &
	    depresist_d0, depresist_d3,   &
	    vsettl_0, vsettl_3 )





























	implicit none

	integer, intent(in) ::   &
	    iflag_coarse
	real, intent(in) ::   &
	    dgnum, alnsg, aerodens,   &
	    temp, pres, airdens, ustar, depresist_a, pblh, zntt, rmolm
	real, intent(inout) ::   &
	    airkinvisc, freepath
	real, intent(out) ::   &
	    depvel_0, depvel_3,   &
	    depresist_d0, depresist_d3,   &
	    vsettl_0, vsettl_3


	
	integer, parameter :: idowescor = 0



        real, parameter :: avo=6.0221367e23

        real, parameter :: rgasuniv=8.314510

        real, parameter :: boltz=rgasuniv/avo

        real, parameter :: grav=9.80622

        real, parameter :: pss0=101325.0

        real, parameter :: tss0=288.15
        real*8, parameter :: pirs=3.14159265358979324
	real, parameter :: sqrt2=1.4142135623731
	real, parameter :: sqrtpi=1.7724539
        real, parameter :: two3=2.0/3.0
        real, parameter :: threepi=3.0*pirs


	real, parameter :: bhat = 1.246   
	real, parameter :: colctr_bigd=2.e-3,colctr_smald=20.e-6  
	                 

	integer :: n

	real :: alnsgy2, amu
	real :: blkta, blkprs, blkdens
	real :: cunq, czh 
	real :: dconst1, dconst2, dconst3, dconst3y
	real :: dgyyy, dq
	real :: eff_dif, eff_imp, eff_int
	real :: ey1, esy04, esy08, esy16, esy20, esy28, esy32, esy64
	real :: knyyy, knq
	real :: nu
	real :: pdensy
	real :: ra, rbcor, rsurfq
	real :: scq, stq, sum0, sum0b, sum3, sum3b
	real :: tmpa
	real :: utscale
	real :: vdplim, vsedq
	real :: xlm
	real :: yextra_vdvg(100,3)

	integer, parameter :: ngausdv = 7
	real, save :: y_gq(ngausdv), wgaus(ngausdv)


	
	
	y_gq = (/ -2.651961356835233, &
	-1.673551628767471, -0.816287882858965, -0.0, &
	0.816287882858965, 1.673551628767471, 2.651961356835233 /)
	wgaus = (/ 0.0009717812450995, &
	0.05451558281913, 0.4256072526101, 0.8102646175568, &
	0.4256072526101, 0.05451558281913, 0.0009717812450995 /)

	blkta = temp
	blkprs = pres
	blkdens = airdens


	if (airkinvisc <= 0.0) then
	    amu = 1.458e-6*blkta*sqrt(blkta)/(blkta+110.4)
	    airkinvisc = amu/airdens
	else
	    amu = airkinvisc*airdens
	end if


	if (freepath <= 0.0) then
	    freepath = 6.6328e-8*pss0*blkta/(tss0*blkprs)
	end if
	xlm = freepath

	ra = depresist_a

	dgyyy = dgnum
	pdensy = aerodens

        alnsgy2 = alnsg**2
        ey1 = exp(0.125*alnsgy2)
        esy04 = ey1**4
        esy08 = esy04*esy04
        esy16 = esy08*esy08
        esy20 = esy16*esy04
        esy28 = esy20*esy08
        esy32 = esy16*esy16
        esy64 = esy32*esy32

	dconst1 = boltz * blkta / ( threepi * amu )
	dconst2 = grav / ( 18.0 * amu )
	dconst3 =  ustar/(9.*amu*colctr_bigd)

	nu = amu / blkdens
                                                                                                                                 
	utscale =  1.
	if (idowescor.eq.1) then
	
	    if (rmolm.lt.0.) then
		czh = -1.*pblh*rmolm
		if (czh.gt.30.0) then
		    utscale=0.45*(czh**0.6667)
		else
		    utscale=1.+((-300.*rmolm)**0.6667)
		endif
	    endif
	endif   
	utscale = ustar*utscale

        sum0 = 0. ; sum0b = 0.
        sum3 = 0. ; sum3b = 0.

        do n = 1,ngausdv
        dq = dgyyy*exp(y_gq(n)*sqrt2*alnsg)  
        knq = 2.*xlm/dq  
        cunq = 1.+knq*(1.257+.4*exp(-1.1/knq))  
        vsedq = pdensy*dconst2*cunq*dq*dq  
        scq = nu*dq/dconst1/cunq  
        eff_dif = scq**(-two3)    
        stq = dconst3*pdensy*dq**2  
        eff_imp = (stq/(0.8+stq))**2   

        eff_int = (0.00116+0.0061*zntt)*dq/1.414e-7 
        if (iflag_coarse > 0) then
           eff_int = min(1.,eff_int)
           rbcor=exp(-2.0*(stq**0.5)) 
        else
           rbcor = 1. 
        endif
        vdplim = utscale*(eff_dif+eff_imp+eff_int)*rbcor

        vdplim = min(vdplim,.02)

        vdplim = max(vdplim,1.e-20)
        rsurfq = ra+1./vdplim
        sum0 = sum0 + wgaus(n)*(vsedq + 1./rsurfq)  
        sum3 = sum3 + wgaus(n)*(vsedq + 1./rsurfq)*dq**3  
        sum0b = sum0b + wgaus(n)*(1./rsurfq)  
        sum3b = sum3b + wgaus(n)*(1./rsurfq)*dq**3  
        end do 

	
        depvel_0 = sum0/sqrtpi  
	
        depvel_3 = sum3/(sqrtpi*exp((1.5*sqrt2*alnsg)**2)*dgyyy**3) 

        sum0b = sum0b/sqrtpi  
        sum3b = sum3b/(sqrtpi*exp((1.5*sqrt2*alnsg)**2)*dgyyy**3) 

	tmpa = 1.0/max( sum0b, 1.0e-10 )
	depresist_d0 = max( (tmpa - ra), 1.0e-10 )
	tmpa = 1.0/max( sum3b, 1.0e-10 )
	depresist_d3 = max( (tmpa - ra), 1.0e-10 )









	vsettl_0 = max( (depvel_0 - sum0b), 0.0 )
	vsettl_3 = max( (depvel_3 - sum3b), 0.0 )


	return
	end subroutine sorgam_aer_drydepvel_2




     subroutine zhang2001_aer_drydepvel_1( &
         dg_num_in, lnsig_in, dens_part_in, &
         airtemp_in, airdens_in, ustar_in, resist_a_in, &
         zhang_optaa, iseason_in, &
         isurftype_in, modesurftype, iwetsurf, &
         airkinvisc, freepath,   &
         ido_mom0, ido_mom3,   &
         depvel_0, depvel_3,   &
         resist_s0, resist_s3,   &
         vsettl_0, vsettl_3 )







      implicit none


      real     :: dg_num_in    
      real     :: lnsig_in     
      real     :: dens_part_in 
      real     :: airdens_in   
      real     :: airtemp_in   
      integer  :: zhang_optaa  
      integer  :: iseason_in   
      integer  :: isurftype_in 
      integer  :: modesurftype 
                               
      integer  :: iwetsurf     
      real     :: resist_a_in  
      real     :: ustar_in     
      real     :: airkinvisc   
      real     :: freepath     
      integer  :: ido_mom0, ido_mom3 


      real     :: depvel_0, depvel_3    
                                        
      real     :: resist_s0, resist_s3  
      real     :: vsettl_0, vsettl_3    


      integer  :: iseason, isurftype
      integer  :: moment     

      real, parameter :: gravit = 9.806
      real, parameter :: boltzmann=1.381e-23 
      real, parameter :: pi = 3.1415926536

      real     :: airdens
      real     :: airtemp
      real     :: beta
      real     :: cuncorrect     
      real     :: dens_part
      real     :: depvel         
      real     :: dg_num
      real     :: diffus_part    
      real     :: dispersion     
                                
      real     :: dp_mom         
      real     :: dp_vm          
      real     :: e_brown        
      real     :: e_impact       
      real     :: e_intcept      
      real     :: e0_brown, e0_impact, e0_intcept   
      real     :: lnsig
      real     :: resist_a
      real     :: resist_sfc     
      real     :: schmidt        
      real     :: stickfrac      
      real     :: stokes         
      real     :: tmpa
      real     :: ustar
      real     :: vsettl         
      real     :: viscosity      

      real, save :: gammax(15)     
      data gammax/ 0.56, 0.58, 0.56, 0.56, 0.56, &
                   0.54, 0.54, 0.54, 0.54, 0.54, &
                   0.54, 0.54, 0.50, 0.50, 0.56/

      real, save :: alpha(15)      
      data alpha/   1.0,   0.6,   1.1,   0.8,   0.8, &
                    1.2,   1.2,  50.0,  50.0,   1.3, &
                    2.0,  50.0, 100.0, 100.0,   1.5/

      real, save :: rad_collector(5,15) 
      data rad_collector/ &
           0.002, 0.002, 0.002, 0.002, 0.002, &  
           0.005, 0.005, 0.005, 0.005, 0.005, &  
           0.002, 0.002, 0.005, 0.005, 0.002, &  
           0.005, 0.005, 0.010, 0.010, 0.005, &  
           0.005, 0.005, 0.005, 0.005, 0.005, &  
           0.002, 0.002, 0.005, 0.005, 0.002, &  
           0.002, 0.002, 0.005, 0.005, 0.002, &  
           9.999, 9.999, 9.999, 9.999, 9.999, &  
           9.999, 9.999, 9.999, 9.999, 9.999, &  
           0.010, 0.010, 0.010, 0.010, 0.010, &  
           0.010, 0.010, 0.010, 0.010, 0.010, &  
           9.999, 9.999, 9.999, 9.999, 9.999, &  
           9.999, 9.999, 9.999, 9.999, 9.999, &  
           9.999, 9.999, 9.999, 9.999, 9.999, &  
           0.010, 0.010, 0.010, 0.010, 0.010  /  

      integer, save  :: isurftype_zhang_from_usgs(25) 
      data isurftype_zhang_from_usgs / &  
               
        15, &  
         7, &  
         7, &  
         7, &  
         7, &  
        10, &  
         6, &  
         6, &  
         6, &  
         6, &  
         4, &  
         3, &  
         2, &  
         1, &  
         5, &  
        14, &  
        11, &  
        10, &  
         8, &  
         9, &  
        10, &  
         9, &  
         8, &  
        12, &  
         8  /  
               
               
               
               
               
           
 

      airtemp   = max( 50.0,   min( 350.0,  airtemp_in ) )
      airdens   = max( 0.01,   min( 5.0,    airdens_in ) )
      dg_num    = max( 1.0e-9, min( 1.0e-2, dg_num_in ) )
      dens_part = max( 0.1e3,  min( 50.0e3, dens_part_in ) )
      lnsig     = max( 0.0,    min( log(5.0), lnsig_in ) )
      resist_a  = max( 1.0e-2, min( 1.0e6,  resist_a_in ) )
      ustar     = max( 1.0e-6, min( 1.0e2,  ustar_in ) )

      isurftype = 8
      if (modesurftype <= 1) then
         if ((isurftype_in >= 1) .and. (isurftype_in <= 15)) &
            isurftype = isurftype_in
      else
         
         if ((isurftype_in >= 1) .and. (isurftype_in <= 25)) &
            isurftype = isurftype_zhang_from_usgs(isurftype_in)
      end if

      iseason = 1
      if ((iseason_in >= 1) .and. (iseason_in <= 5)) &
         iseason = iseason_in













      if (zhang_optaa == 2) then
         e0_brown   = 1.0
         e0_impact  = 1.0
         e0_intcept = 4.0
      else
         e0_brown   = 3.0
         e0_impact  = 3.0
         e0_intcept = 3.0
      end if





      if (freepath <= 0.0) then
         freepath = 6.6328e-8*(1.2255/airdens)
      end if

      if (airkinvisc <= 0.0) then


         viscosity = 1.458e-6*airtemp*sqrt(airtemp)/(airtemp+110.4)
         airkinvisc = viscosity/airdens
      else
         viscosity = airkinvisc*airdens
      end if

      dp_vm = dg_num*exp(1.5*lnsig*lnsig)

      do moment = 0, 3, 3

      if (moment == 0) then
         if (ido_mom0 <= 0) then
            depvel_0 = 0.0
            vsettl_0 = 0.0
            resist_s0 = 1.0e30
            cycle
         end if
      else if (moment == 3) then
         if (ido_mom3 <= 0) then
            depvel_3 = 0.0
            vsettl_3 = 0.0
            resist_s3 = 1.0e30
            cycle
         end if
      else
         cycle
      end if


      dp_mom = dg_num*exp(float(moment)*lnsig*lnsig)
      dp_mom = min( 100.0e-6, dp_mom )

      tmpa = 2.0*freepath/dp_mom
      cuncorrect = 1. + tmpa*(1.257 + 0.4*exp(-1.10/tmpa))
      dispersion = exp(2.*lnsig*lnsig)
      vsettl = ( dens_part*(dp_mom**2)*gravit*   &
                        cuncorrect / (18.0*viscosity) ) * dispersion

      diffus_part = boltzmann*airtemp*cuncorrect/   &
                         (3.0*pi*viscosity*dp_mom)
      schmidt=airkinvisc/diffus_part
      e_brown = schmidt**(-gammax(isurftype))

      if (isurftype==8  .or. isurftype==9  .or. &
          isurftype==12 .or. isurftype==13 .or. isurftype==14) then

         stokes = vsettl*ustar*ustar/(gravit*airkinvisc)
         e_intcept = 0.
      else

         stokes = vsettl*ustar/(gravit*rad_collector(iseason,isurftype))
         if (zhang_optaa == 2) then
            tmpa = 0.5*dp_mom/rad_collector(iseason,isurftype)
            e_intcept = tmpa
         else
            tmpa = dp_mom/rad_collector(iseason,isurftype)
            e_intcept = 0.5 * tmpa**2
         endif
      endif

      beta = 2.0
      e_impact = ( stokes / (alpha(isurftype)+stokes) )**beta


      if ((iwetsurf > 0) .or. &
          (isurftype == 13) .or. (isurftype == 14) .or. &
          (dp_vm < 5.0e-6)) then
         stickfrac = 1.0
      else
         stickfrac = exp(-sqrt(stokes))
      end if



      tmpa = ustar*(e0_brown*e_brown + e0_intcept*e_intcept + e0_impact*e_impact)*stickfrac
      resist_sfc = 1./max(tmpa,1.0e-30)



      depvel = vsettl + 1./(resist_a + resist_sfc + resist_a*resist_sfc*vsettl)

      if (moment == 0) then
         depvel_0 = depvel
         vsettl_0 = vsettl
         resist_s0 = resist_sfc
      else
         depvel_3 = depvel
         vsettl_3 = vsettl
         resist_s3 = resist_sfc
      end if

      end do 

      return
      end subroutine zhang2001_aer_drydepvel_1


        subroutine cam_mam_aer_drydep_prep(                                 &
                id, ktau, dtstep, config_flags, aer_mech_id,               &
                gmt, julday,                                               &
                t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
                moist, chem,                                               &
                nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
                yextra_add,                                                &
                i, j,                                                      &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte                                  )

        use module_configure, only:  num_moist, num_chem, &
                grid_config_rec_type
        use module_state_description, only:  param_first_scalar

        use module_data_cam_mam_asect

        implicit none


        integer, intent(in) ::   &
                id, ktau, julday, aer_mech_id, nbin_add,   &
                i, j

        integer, intent(in) ::   &
                ids, ide, jds, jde, kds, kde,   &
                ims, ime, jms, jme, kms, kme,   &
                its, ite, jts, jte, kts, kte

        real, intent(in) :: dtstep, gmt

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                t_phy, rho_phy, p_phy, alt, p8w, t8w

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
                moist

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
                chem

        real, intent(inout),   &
                dimension( 1:nbin_add ) :: &
                alnsg_add, wetdgnum_add, wetdens_add

        real, intent(inout),   &
                dimension( 1:nbin_add, 10 ) :: &
                yextra_add

        type(grid_config_rec_type), intent(in) :: config_flags



        integer :: ibin, iok, iphase, itype
        integer :: itmpa
        integer :: k
        integer :: l, ll, l1, n
        integer :: p1st

        real, parameter :: densdefault = 2.0
        real, parameter :: smallmassaa = 1.0e-20
        real, parameter :: smallmassbb = 1.0e-30
        real, parameter :: pi = 3.1415926536
        real, parameter :: piover6 = pi/6.0
        real, parameter :: onethird = 1.0/3.0

        real :: drydens, drydp, drymass, dryvol
        real :: rnum
        real :: tmpa, tmp_lnsg
        real :: wetdgnum, wetdens, wetdp, wetmass, wetvol
        real :: wetdgnum_si, wetdens_si

        character(len=160) :: txtaa


        n = sum( nsize_aer(1:ntype_aer) )
        if (n /= nbin_add) then
            write(txtaa,'(a,2(1x,i10))') &
                '*** cam_mam_aer_drydep_prep - bad nbin_add, sum(nsize)', &
                nbin_add, n
            write(*,'(a)') trim(txtaa)
            call wrf_error_fatal3("<stdin>",2036,&
txtaa )
            stop
        end if



        p1st = param_first_scalar

        k = kts


        iphase = ai_phase


        ibin = 0
        do 2900 itype = 1, ntype_aer
        do 2800 n = 1, nsize_aer(itype)
        ibin = ibin + 1


        dryvol = 0.0
        drymass = 0.0
        do ll = 1, ncomp_aer(itype)
            l1 = massptr_aer(ll,n,itype,iphase)
            tmpa = max( chem(i,k,j,l1), 0.0 )
            drymass = drymass + tmpa
            dryvol = dryvol + tmpa/dens_aer(ll,itype)
        end do

        l1 = waterptr_aer(n,itype)
        tmpa = max( chem(i,k,j,l1), 0.0 )
        wetmass = drymass + tmpa
        wetvol = dryvol + tmpa/dens_water_aer

        l1 = numptr_aer(n,itype,iphase)
        rnum = max( chem(i,k,j,l1), 0.0 )

        drymass = drymass*28.966e-9 
        wetmass = wetmass*28.966e-9
        dryvol  = dryvol *28.966e-9 
        wetvol  = wetvol *28.966e-9
        rnum    = rnum*28.966e-3    

        if (drymass <= smallmassbb) then
            drydp = dcen_sect(n,itype)
            drydens = densdefault
            wetdp = drydp
            wetdens = drydens

        else
            if (drymass <= smallmassaa) then
                wetmass = drymass
                wetvol = dryvol
            end if
            drydens = drymass/dryvol
            wetdens = wetmass/wetvol

            if (rnum >= dryvol/volumlo_sect(n,itype)) then
                drydp = dlo_sect(n,itype)
            else if (rnum <= dryvol/volumhi_sect(n,itype)) then
                drydp = dhi_sect(n,itype)
            else
                drydp = (dryvol/(piover6*rnum))**onethird
            end if

            if (abs(wetvol) > (1000.*abs(dryvol))) then
              tmpa=10.0
            else
              tmpa=abs(wetvol/dryvol)**onethird
              tmpa=max(1.0,min(tmpa,10.0))
            endif
            wetdp = drydp*tmpa

        end if




        tmp_lnsg = log( sigmag_aer(n,itype) )
        wetdgnum = wetdp * exp( -1.5*tmp_lnsg*tmp_lnsg )

        wetdgnum_add(ibin) = wetdgnum * 1.0e-2 
        wetdens_add(ibin) = wetdens * 1.0e3    
        alnsg_add(ibin) = tmp_lnsg


        yextra_add(ibin,1) = drymass/28.966e-9 
        yextra_add(ibin,5) = wetmass/28.966e-9
        yextra_add(ibin,2) = dryvol/28.966e-6  
        yextra_add(ibin,6) = wetvol/28.966e-6
        yextra_add(ibin,3) = drydens*1.0e3     
        yextra_add(ibin,7) = wetdens*1.0e3
        yextra_add(ibin,4) = drydp*1.0e-2      
        yextra_add(ibin,8) = wetdp*1.0e-2

        yextra_add(ibin,9) = dcen_sect(n,itype)*1.0e-2   
        yextra_add(ibin,10) = rnum/28.966e-3    


2800    continue
2900    continue


        return
        end subroutine cam_mam_aer_drydep_prep




        subroutine cam_mam_aer_drydep_load_ddvel(                           &
                id, ktau, dtstep, config_flags, aer_mech_id,               &
                ddvel,                                                     &
                nbin_add, depvel0_add, depvel3_add,                        &
                i, j,                                                      &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte                                  )

        use module_configure, only:  num_moist, num_chem, &
                grid_config_rec_type
        use module_state_description, only:  param_first_scalar

        use module_data_cam_mam_asect

        implicit none


        integer, intent(in) ::   &
                id, ktau, aer_mech_id, nbin_add,   &
                i, j

        integer, intent(in) ::   &
                ids, ide, jds, jde, kds, kde,   &
                ims, ime, jms, jme, kms, kme,   &
                its, ite, jts, jte, kts, kte

        real, intent(in) :: dtstep

        real, intent(in),   &
                dimension( 1:nbin_add ) :: &
                depvel0_add, depvel3_add

        real, intent(out),   &
                dimension( its:ite, jts:jte, 1:num_chem ) :: &
                ddvel

        type(grid_config_rec_type), intent(in) :: config_flags



        integer :: ibin, iphase, itype
        integer :: l, ll, n

        character(len=160) :: txtaa


        n = sum( nsize_aer(1:ntype_aer) )
        if (n /= nbin_add) then
            write(txtaa,'(a,2(1x,i10))') &
                '*** cam_mam_aer_drydep_load_ddvel - bad nbin_add, sum(nsize)', &
                nbin_add, n
            write(*,'(a)') trim(txtaa)
            call wrf_error_fatal3("<stdin>",2199,&
txtaa )
            stop
        end if



        do iphase = 1, nphase_aer
        ibin = 0
        do itype = 1, ntype_aer
        do n = 1, nsize_aer(itype)
            ibin = ibin + 1
            do ll = -2, ncomp_plustracer_aer(itype)
                if (iphase /= ai_phase) cycle

                if (ll == -2) then
                    l = numptr_aer(n,itype,iphase)
                else if (ll == -1) then
                    l = -1
                    if (iphase == ai_phase) l = waterptr_aer(n,itype)
                else if (ll == 0) then
                    l = -1
                    if (iphase == ai_phase) l = hyswptr_aer(n,itype)
                else
                    l = massptr_aer(ll,n,itype,iphase)
                end if

                if ((l >= param_first_scalar) .and. (l <= num_chem)) then
                    if (l == numptr_aer(n,itype,iphase)) then
                        ddvel(i,j,l) = depvel0_add(ibin)
                    else
                        ddvel(i,j,l) = depvel3_add(ibin)
                    end if
                end if
            end do
        end do
        end do
        end do


        return
        end subroutine cam_mam_aer_drydep_load_ddvel




END MODULE module_aer_drydep
