







	module module_mosaic_drydep









	contains



	subroutine mosaic_drydep_driver(                                   &
		id, curr_secs, ktau, dtstep, config_flags,                     &
		gmt, julday,                                                   &
		t_phy, rho_phy, p_phy,                                         &
		ust, aer_res,                                                  &
		moist, chem, ddvel,                                            &
		ids,ide, jds,jde, kds,kde,                                     &
		ims,ime, jms,jme, kms,kme,                                     &
		its,ite, jts,jte, kts,kte                                      )

	use module_configure, only:  grid_config_rec_type, num_moist, num_chem
	use module_state_description, only:  param_first_scalar

	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_mosaic_driver, only:  mapaer_tofrom_host
	use module_peg_util, only:  peg_error_fatal

	implicit none


	integer, intent(in) ::   &
		id, ktau, julday

	integer, intent(in) ::   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

    real(kind=8), intent(in) :: curr_secs

	real, intent(in) :: dtstep, gmt

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		t_phy, rho_phy, p_phy

	real, intent(in),   &
		dimension( ims:ime, jms:jme ) :: &
		ust

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


	type(grid_config_rec_type), intent(in) :: config_flags



	integer idum, jdum
	integer it, jt, kt
	integer iphase, itype
	integer ktmaps, ktmape
	integer ll, l1, n
	integer levdbg_err, levdbg_info

	integer idiagaa_dum, ijcount_dum

	real dum
	real vdep_aer(maxd_asize,maxd_atype,maxd_aphase)

	character*100 msg






























	ktmaps = kts
	ktmape = kts
	ktot = 1



	lunerr = -1
	lunout = -1
	levdbg_err = 0
        levdbg_info = 15

	iymdcur = 20030822
	ihmscur = 0
	ihmscur = nint( mod( real(gmt*3600.,8)+curr_secs, 86400.0_8 ) )

	t = curr_secs
	ncorecnt = ktau - 1












	itot = ite
	jtot = jte
	nsubareas = 1

	ijcount_dum = 0


	do 2920 jt = jts, jte
	do 2910 it = its, ite

	ijcount_dum = ijcount_dum + 1


	call mapaer_tofrom_host( 0,                       &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      ktmaps,ktmape,              &
		num_moist, num_chem, moist, chem,             &
		t_phy, p_phy, rho_phy                         )



	idiagaa_dum = 1
	idiagaa_dum = 0
	if ((jt.ne.jts) .and. (jt.ne.jte) .and.   &
			(jt.ne.(jts+jte)/2)) idiagaa_dum = 0
	if ((it.ne.its) .and. (it.ne.ite) .and.   &
			(it.ne.(its+ite)/2)) idiagaa_dum = 0

	call mosaic_drydep_1clm( idiagaa_dum, it, jt,   &
		ust(it,jt), aer_res(it,jt), vdep_aer )



	do iphase = 1, nphase_aer
	do itype = 1, ntype_aer
	do n = 1, nsize_aer(itype)
	    do ll = -2, ncomp_plustracer_aer(itype)
		if (ll .eq. -2) then
        	    l1 = numptr_aer(n,itype,iphase)
		else if (ll .eq. -1) then
		    l1 = -1
		    if (iphase .eq. ai_phase) l1 = waterptr_aer(n,itype)
		else if (ll .eq. 0) then
		    l1 = -1
		    if (iphase .eq. ai_phase) l1 = hyswptr_aer(n,itype)
		else
		    l1 = massptr_aer(ll,n,itype,iphase)
		end if
		if (l1 .ge. param_first_scalar) then
		    ddvel(it,jt,l1) = vdep_aer(n,itype,iphase)
		end if
	    end do
	end do
	end do
	end do


2910	continue
2920	continue



	return
	end subroutine mosaic_drydep_driver



	subroutine mosaic_drydep_1clm( idiagaa, it, jt,   &
		ustar_in, depresist_a_in, vdep_aer )

	use module_configure, only:  grid_config_rec_type

	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_mosaic_driver, only:  mapaer_tofrom_host
	use module_peg_util, only:  peg_error_fatal

	implicit none


	integer, intent(in) :: idiagaa, it, jt


	real, intent(in) :: ustar_in

	real, intent(in) :: depresist_a_in


	real, intent(inout) :: vdep_aer(maxd_asize,maxd_atype,maxd_aphase)


	real, parameter :: densdefault = 2.0
	real, parameter :: smallmassaa = 1.0e-20
	real, parameter :: smallmassbb = 1.0e-30
	real, parameter :: piover6 = pi/6.0
	real, parameter :: onethird = 1.0/3.0

	integer iphase, itype, k, ll, l1, m, n

	real airdens, airkinvisc
	real depresist_a, depresist_unstabpblfact
	real depresist_d0, depresist_d3
	real depvel_a0, depvel_a3
	real drydens, drydp, drymass, dryvol
	real dum, dumalnsg, dumfact, dummass
	real freepath
	real rnum
	real temp
	real ustar
	real vsettl_0, vsettl_3
	real wetdgnum, wetdens, wetdp, wetmass, wetvol



	k = 1
	m = 1


	temp = rsub(ktemp,k,m)


	airdens = cairclm(1)*28.966

	airkinvisc = ( 1.8325e-4 * (416.16/(temp+120.0)) *   &
      				((temp/296.16)**1.5) ) / airdens

	freepath = 7.39758e-4 * airkinvisc / sqrt(temp)

	ustar = ustar_in * 100.0

	depresist_a = depresist_a_in * 0.01


	depresist_unstabpblfact = 1.0



	vdep_aer(:,:,:) = 0.0


	iphase = ai_phase


	do 2900 itype = 1, ntype_aer
	do 2800 n = 1, nsize_aer(itype)

	dryvol = 0.0
	drymass = 0.0
	do ll = 1, ncomp_aer(itype)
	    l1 = massptr_aer(ll,n,itype,iphase)
	    dummass = rsub(l1,k,m)*mw_aer(ll,itype)
	    drymass = drymass + dummass
	    dryvol = dryvol + dummass/dens_aer(ll,itype)
	end do

	l1 = waterptr_aer(n,itype)
	dummass = rsub(l1,k,m)*mw_water_aer
	wetmass = drymass + dummass
	wetvol = dryvol + dummass/dens_water_aer

	l1 = numptr_aer(n,itype,iphase)
	rnum = rsub(l1,k,m)

	if (drymass .le. smallmassbb) then
	    drydp = dcen_sect(n,itype)
	    drydens = densdefault
	    wetdp = drydp
	    wetdens = drydens
	    goto 1900
	end if


	if (drymass .le. smallmassaa) then
	    wetmass = drymass
	    wetvol = dryvol
	end if
	drydens = drymass/dryvol
	wetdens = wetmass/wetvol


	if (rnum .ge. dryvol/volumlo_sect(n,itype)) then
	    drydp = dlo_sect(n,itype)
	else if (rnum .le. dryvol/volumhi_sect(n,itype)) then
	    drydp = dhi_sect(n,itype)
	else
	    drydp = (dryvol/(piover6*rnum))**onethird
	end if



        if(abs(wetvol).gt.(1000.*abs(dryvol))) then
          dumfact=10.0
        else
          dumfact=abs(wetvol/dryvol)**onethird
          dumfact=max(1.0,min(dumfact,10.0))
        endif

	wetdp = drydp*dumfact

1900	continue






	dumalnsg = log( 1.0 )
	wetdgnum = wetdp * exp( -1.5*dumalnsg*dumalnsg )
	call aerosol_depvel_2(   &
      	    wetdgnum, dumalnsg, wetdens,   &
      	    temp, airdens, airkinvisc, freepath,   &
      	    ustar, depresist_unstabpblfact,   &
      	    depresist_d0, vsettl_0,   &
      	    depresist_d3, vsettl_3 )




	dum = depresist_a + depresist_d3 +   &
      		    depresist_a*depresist_d3*vsettl_3
	depvel_a3 = vsettl_3 + (1./dum)

	dum = depresist_a + depresist_d0 +   &
      		    depresist_a*depresist_d0*vsettl_0
	depvel_a0 = vsettl_0 + (1./dum)


	vdep_aer(n,itype,iphase) = depvel_a3 * 0.01





	if (idiagaa>0) print 9310, it, jt, n, itype, iphase,   &
		dcen_sect(n,itype), drydp, wetdp,   &
		drydens, wetdens, vdep_aer(n,itype,iphase),   &
		vsettl_3, depresist_d3, depresist_a
9310	format( 'aerdep', 2i4, 3i3, 1p, 3e10.2,   &
		2x, 0p, 2f5.2, 2x, 1p, 4e10.2 )


2800	continue
2900	continue


	return
	end subroutine mosaic_drydep_1clm




	subroutine aerosol_depvel_2(   &
      	    dgnum, alnsg, aerodens,   &
      	    temp, airdens, airkinvisc, freepath,   &
      	    ustar, depresist_unstabpblfact,   &
      	    depresist_d0, vsettl_0,   &
      	    depresist_d3, vsettl_3 )



























	implicit none

	real dgnum, alnsg, aerodens,   &
      	    temp, airdens, airkinvisc, freepath,   &
      	    ustar, depresist_unstabpblfact,   &
      	    depresist_d0, vsettl_0,   &
      	    depresist_d3, vsettl_3

	real aerodiffus_0, schmidt_0, stokes_0, facdepresist_d0
	real aerodiffus_3, schmidt_3, stokes_3, facdepresist_d3
	common / aerosol_depvel_cmn01 /   &
      		aerodiffus_0, schmidt_0, stokes_0, facdepresist_d0,   &
      		aerodiffus_3, schmidt_3, stokes_3, facdepresist_d3

	real xknudsen, xknudsenfact, alnsg2, duma, dumb,   &
      	vsettl_dgnum, aerodiffus_dgnum

	real pi
	parameter (pi = 3.1415926536)

	real gravity
	parameter (gravity = 980.616)

	real boltzmann
	parameter (boltzmann = 1.3807e-16)

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
	
	duma = (schmidt_0**(-0.66666666)) +   &
      		(10.**(-3./max(0.03,stokes_0)))

	dumb = duma*ustar*depresist_unstabpblfact
	depresist_d0 = 1./max( dumb, 1.e-20 )
	facdepresist_d0 = duma

	duma = (schmidt_3**(-0.66666666)) +   &
      		(10.**(-3./max(0.03,stokes_3)))

	dumb = duma*ustar*depresist_unstabpblfact
	depresist_d3 = 1./max( dumb, 1.e-20 )
	facdepresist_d3 = duma

	return
	end subroutine aerosol_depvel_2 



	end module module_mosaic_drydep
