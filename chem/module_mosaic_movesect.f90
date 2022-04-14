







	module module_mosaic_movesect


	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_peg_util









	integer, parameter :: apply_n1_inflow = 0

	contains





























	subroutine move_sections( iflag, iclm, jclm, k, m)














































	implicit none








	integer iflag, iclm, jclm, k, m


	integer idiag_movesect, iphase, itype,   &
	  l, ll, llhysw, llwater, lnew, lold, l3,   &
      	  method_movesect, n, nnew, nold
	integer nnewsave(2,maxd_asize)

	real densdefault, densh2o, smallmassaa, smallmassbb
	real delta_water_conform1, delta_numb_conform1

	real drydenspp(maxd_asize), drydensxx0(maxd_asize),   &
      	     drydensxx(maxd_asize), drydensyy(maxd_asize)
	real drymasspp(maxd_asize), drymassxx0(maxd_asize),   &
      	     drymassxx(maxd_asize), drymassyy(maxd_asize)
	real dryvolxx(maxd_asize), dryvolyy(maxd_asize)
	real rmassxx(maxd_acomp+2,maxd_asize),   &
      	     rmassyy(maxd_acomp+2,maxd_asize)
	real rnumbpp(maxd_asize), rnumbxx0(maxd_asize),   &
      	     rnumbxx(maxd_asize), rnumbyy(maxd_asize)
	real specdensxx(maxd_acomp), specmwxx(maxd_acomp)
	real xferfracvol(2,maxd_asize), xferfracnum(2,maxd_asize)
	real wetvolxx(maxd_asize), wetvolyy(maxd_asize)
	real wetmassxx(maxd_asize), wetmassyy(maxd_asize)

	character*160 msg





	if (msectional .le. 0) return
	if (ntype_aer .le. 0) return
	if (nphase_aer .le. 0) return

 










	method_movesect = mod( msectional, 100 )
	if (method_movesect .le. 10) method_movesect = 10


	idiag_movesect  = mod( msectional, 10000 )/100

	if      ((method_movesect .eq. 10) .or.   &
      		 (method_movesect .eq. 11) .or.   &
      		 (method_movesect .eq. 20) .or.   &
      		 (method_movesect .eq. 21) .or.   &
      		 (method_movesect .eq. 30) .or.   &
      		 (method_movesect .eq. 31)) then
	    continue
	else if ((method_movesect .eq. 19) .or.   &
      		 (method_movesect .eq. 29) .or.   &
      		 (method_movesect .eq. 39)) then
	    return
	else
	    msg = '*** subr move_sections error - ' //   &
		'illegal value for msectional'
	    call peg_error_fatal( lunerr, msg )
	end if

	if (iabs(iflag) .eq. 1) then
	    iphase = ai_phase
	else if (iabs(iflag) .eq. 2) then
	    iphase = cw_phase
	    if (nphase_aer .lt. 2) then
		msg = '*** subr move_sections error - ' //   &
		    'iflag=2 (after aqueous chemistry) but nphase_aer < 2'
		call peg_error_fatal( lunerr, msg )
	    else if (cw_phase .ne. 2) then
		msg = '*** subr move_sections error - ' //   &
		    'iflag=2 (after aqueous chemistry) but cw_phase .ne. 2'
		call peg_error_fatal( lunerr, msg )
	    end if
	else
	    msg = '*** subr move_sections error - ' //   &
		'iabs(iflag) must be 1 or 2'
	    call peg_error_fatal( lunerr, msg )
	end if




	if (iflag .le. 0) then
	    write(msg,9040) 'method', method_movesect
	    call peg_message( lunout, msg )
	    write(msg,9040) 'idiag ', idiag_movesect
	    call peg_message( lunout, msg )
	    call move_sections_checkptrs( iflag, iclm, jclm, k, m )
	    return
	end if
9040	format( '*** subr move_sections - ', a, ' =', i6 )


	if (idiag_movesect .eq. 70) then
	    msg = ' '
	    call peg_message( lunout, msg )
	    write(msg,9060) iclm, jclm, k, m, msectional
	    call peg_message( lunout, msg )
	end if
9060	format( '*** subr move_sections diagnostics i,j,k,m,msect =', 4i4, i6 )


	densdefault = 2.0
	densh2o = 1.0



	smallmassaa = 1.0e-20


	smallmassbb = 1.0e-30



	do 1900 itype = 1, ntype_aer

	if (nsize_aer(itype) .le. 0) goto 1900

	call move_sections_initial_conform(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drydensxx0, drymassxx0, rnumbxx0,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  specmwxx, specdensxx )

	if (method_movesect .le. 19) then
	call move_sections_calc_movingcenter(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  xferfracvol, xferfracnum )
	else
	call move_sections_calc_masnumadvect(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drydensxx0, drymassxx0, rnumbxx0,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  xferfracvol, xferfracnum )
	end if

	call move_sections_apply_moves(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy,   &
      	  xferfracvol, xferfracnum )



	if ((apply_n1_inflow .eq. 1) .and.   &
	    (iphase .eq. ai_phase)) then
	call move_sections_apply_n1_inflow(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy,   &
      	  xferfracvol, xferfracnum,   &
      	  specmwxx, specdensxx )
	end if




1900	continue

	return
	end subroutine move_sections                          



	subroutine move_sections_initial_conform(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drydensxx0, drymassxx0, rnumbxx0,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  specmwxx, specdensxx )











	implicit none








	integer iflag, iclm, jclm, iphase, itype, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	real densdefault, densh2o, smallmassaa, smallmassbb
	real delta_water_conform1, delta_numb_conform1
	real drydenspp(maxd_asize), drydensxx0(maxd_asize),   &
      	     drydensxx(maxd_asize)
	real drymasspp(maxd_asize), drymassxx0(maxd_asize),   &
      	     drymassxx(maxd_asize)
	real dryvolxx(maxd_asize)
	real rmassxx(maxd_acomp+2,maxd_asize)
	real rnumbpp(maxd_asize), rnumbxx0(maxd_asize),   &
      	     rnumbxx(maxd_asize)
	real specdensxx(maxd_acomp), specmwxx(maxd_acomp)
	real wetvolxx(maxd_asize)
	real wetmassxx(maxd_asize)



	integer l, ll, lnew, lold, l3, n, nnew, nold

	real dummass, dumnum, dumnum_at_dhi, dumnum_at_dlo, dumr,   &
      	  dumvol, dumvol1p, dumwatrmass



	do l = 1, ltot2
	    rsub(l,k,m) = max( 0., rsub(l,k,m) )
	end do


	llhysw = ncomp_plustracer_aer(itype) + 1
	llwater = ncomp_plustracer_aer(itype) + 2
	do n = 1, nsize_aer(itype)
	    do ll = 1, ncomp_plustracer_aer(itype)
		l = massptr_aer(ll,n,itype,iphase)
		rmassxx(ll,n) = rsub(l,k,m)
	    end do
	    rmassxx(llhysw,n) = 0.
	    l = 0
	    if (iphase .eq. ai_phase) l = hyswptr_aer(n,itype)
	    if (l .gt. 0) rmassxx(llhysw,n) = rsub(l,k,m)
	    rmassxx(llwater,n) = 0.
	    l = 0
	    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
	    if (l .gt. 0) rmassxx(llwater,n) = rsub(l,k,m)

	    rnumbxx(n)  = rsub(numptr_aer(n,itype,iphase),k,m)
	    rnumbxx0(n) = rnumbxx(n)
	    rnumbpp(n)  = rnumbxx(n)

	    drydenspp(n) = drydens_pregrow(n,itype)
	    drymasspp(n) = drymass_pregrow(n,itype)

	    drydensxx(n) = drydens_aftgrow(n,itype)
	    drymassxx(n) = drymass_aftgrow(n,itype)
	    drydensxx0(n) = drydensxx(n)
	    drymassxx0(n) = drymassxx(n)
	end do


	do ll = 1, ncomp_plustracer_aer(itype)
	    specmwxx(ll) = mw_aer(ll,itype)
	    specdensxx(ll) = dens_aer(ll,itype)
	end do

	delta_water_conform1 = 0.0
	delta_numb_conform1 = 0.0


	do 1390 n = 1, nsize_aer(itype)








	if ( (drydensxx(n) .lt.  0.1) .or.   &
	     (drydensxx(n) .gt. 20.0) .or.   &
      	     (drymassxx(n) .le. smallmassaa) ) then
	    dummass = 0.
	    dumvol = 0.
	    do ll = 1, ncomp_aer(itype)
		dumr = rmassxx(ll,n)*specmwxx(ll)
		dummass = dummass + dumr
		dumvol  = dumvol  + dumr/specdensxx(ll)
	    end do
	    drymassxx(n) = dummass
	    if (min(dummass,dumvol) .le. smallmassbb) then
		drydensxx(n) = densdefault
		dumvol = dummass/densdefault
		dumnum = dummass/(volumcen_sect(n,itype)*densdefault)
	    else
		drydensxx(n) = dummass/dumvol
		dumnum = rnumbxx(n)
		dumnum_at_dhi = dumvol/volumhi_sect(n,itype)
		dumnum_at_dlo = dumvol/volumlo_sect(n,itype)
		dumnum = max( dumnum_at_dhi, min( dumnum_at_dlo, dumnum ) )
	    end if
	    delta_numb_conform1 = delta_numb_conform1 + dumnum - rnumbxx(n)
	    rnumbxx(n) = dumnum
	    rnumbpp(n) = rnumbxx(n)
	    delta_water_conform1 = delta_water_conform1 - rmassxx(llwater,n) 
	    rmassxx(llwater,n) = 0.
	end if



	dryvolxx(n) = drymassxx(n)/drydensxx(n)
	dumwatrmass = rmassxx(llwater,n)*mw_water_aer
	wetmassxx(n) = drymassxx(n) + dumwatrmass
	wetvolxx(n) = dryvolxx(n) + dumwatrmass/densh2o

1390	continue

	return
	end subroutine move_sections_initial_conform                          



	subroutine move_sections_calc_movingcenter(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  xferfracvol, xferfracnum )








	implicit none








	integer iflag, iclm, jclm, iphase, itype, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize)
	real densdefault, densh2o, smallmassaa, smallmassbb
	real drydensxx(maxd_asize)
	real drymassxx(maxd_asize)
	real dryvolxx(maxd_asize)
	real rmassxx(maxd_acomp+2,maxd_asize)
	real rnumbxx(maxd_asize)
	real xferfracvol(2,maxd_asize), xferfracnum(2,maxd_asize)
	real wetmassxx(maxd_asize)
	real wetvolxx(maxd_asize)


	integer ll, n, ndum, nnew, nold
	real dumnum, dumvol, dumvol1p, sixoverpi, third
	character*160 msg


	sixoverpi = 6.0/pi
	third = 1.0/3.0





	do 1390 n = 1, nsize_aer(itype)

	nnew = n


	if (drymassxx(n) .le. smallmassaa) goto 1290

	dumvol = dryvolxx(n)
	dumnum = rnumbxx(n)



	if (dumnum .le. dumvol/volumhi_sect(nsize_aer(itype),itype)) then
	    nnew = nsize_aer(itype)
	    goto 1290

	else if (dumnum .ge. dumvol/volumlo_sect(1,itype)) then
	    nnew = 1
	    goto 1290
	end if


	dumvol1p = dumvol/dumnum
	if (dumvol1p .gt. volumhi_sect(n,itype)) then
	    do while ( (nnew .lt. nsize_aer(itype)) .and.   &
      		       (dumvol1p .gt. volumhi_sect(nnew,itype)) )
		nnew = nnew + 1
	    end do

	else if (dumvol1p .lt. volumlo_sect(n,itype)) then
	    do while ( (nnew .gt. 1) .and.   &
      		       (dumvol1p .lt. volumlo_sect(nnew,itype)) )
		nnew = nnew - 1
	    end do

	end if

1290	nnewsave(1,n) = nnew
	nnewsave(2,n) = 0

	xferfracvol(1,n) = 1.0
	xferfracvol(2,n) = 0.0
	xferfracnum(1,n) = 1.0
	xferfracnum(2,n) = 0.0

1390	continue




	if (idiag_movesect .eq. 70) then
	    ndum = 0
	    do n = 1, nsize_aer(itype)
		if (nnewsave(1,n) .ne. n) ndum = ndum + 1
	    end do
	    if (ndum .gt. 0) then
		write(msg,97751) 'YES', iclm, jclm, k, m,   &
      		ndum, (nnewsave(1,n), n=1,nsize_aer(itype))
		call peg_message( lunout, msg )
	    else
		write(msg,97751) 'NO ', iclm, jclm, k, m,   &
      		ndum, (nnewsave(1,n), n=1,nsize_aer(itype))
		call peg_message( lunout, msg )
	    end if
	end if
97751	format( 'movesect', a, 4i3, 3x, i3, 3x, 14i3 )

	return
	end subroutine move_sections_calc_movingcenter                          



	subroutine move_sections_calc_masnumadvect(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drydensxx0, drymassxx0, rnumbxx0,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  xferfracvol, xferfracnum )













	implicit none








	integer iflag, iclm, jclm, iphase, itype, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize)

	real densdefault, densh2o, smallmassaa, smallmassbb
	real drydenspp(maxd_asize), drydensxx0(maxd_asize),   &
      	     drydensxx(maxd_asize)
	real drymasspp(maxd_asize), drymassxx0(maxd_asize),   &
      	     drymassxx(maxd_asize)
	real dryvolxx(maxd_asize)
	real rmassxx(maxd_acomp+2,maxd_asize)
	real rnumbpp(maxd_asize), rnumbxx0(maxd_asize),   &
      	     rnumbxx(maxd_asize)
	real xferfracvol(2,maxd_asize), xferfracnum(2,maxd_asize)
	real wetvolxx(maxd_asize)
	real wetmassxx(maxd_asize)


	integer ierr, n, nnew, nnew2
	integer iforce_movecenter(maxd_asize)

	real dum1, dum2, dum3
	real dumaa, dumbb, dumgamma, dumratio
	real dumfracnum, dumfracvol
	real dumntot
	real dumv
	real dumvbar_aft, dumvbar_pre
	real dumvcutlo_nnew_pre, dumvcuthi_nnew_pre
	real dumvlo_pre, dumvhi_pre, dumvdel_pre
	real dumvtot_aft, dumvtot_pre
	real dumzlo, dumzhi
	real sixoverpi, third

	character*4 dumch4
	character*1 dumch1
	character*160 msg


	sixoverpi = 6.0/pi
	third = 1.0/3.0











	do 3900 n = 1, nsize_aer(itype)

	nnew = n
	iforce_movecenter(n) = 0

	xferfracvol(1,n) = 1.0
	xferfracvol(2,n) = 0.0
	xferfracnum(1,n) = 1.0
	xferfracnum(2,n) = 0.0

	dumvtot_aft = -1.0
	dumvtot_pre = -1.0
	dumvbar_aft = -1.0
	dumvbar_pre = -1.0
	dumvlo_pre = -1.0
	dumvhi_pre = -1.0
	dumgamma = -1.0
	dumratio = -1.0
	dumvcutlo_nnew_pre = volumlo_sect(nnew,itype)*(dumvbar_pre/dumvbar_aft)
	dumvcuthi_nnew_pre = volumhi_sect(nnew,itype)*(dumvbar_pre/dumvbar_aft)
	dumfracvol = -1.0
	dumfracnum = -1.0


	if (drymassxx(n) .le. smallmassaa) then
	    iforce_movecenter(n) = 1
	    goto 1290
	end if

	dumvtot_aft = dryvolxx(n)
	dumntot = rnumbxx(n)



	if (dumntot .le. dumvtot_aft/volumhi_sect(nsize_aer(itype),itype)) then
	    nnew = nsize_aer(itype)
	    iforce_movecenter(n) = 2
	    goto 1290
	else if (dumntot .ge. dumvtot_aft/volumlo_sect(1,itype)) then
	    nnew = 1
	    iforce_movecenter(n) = 3
	    goto 1290
	end if



	dumvbar_aft = dumvtot_aft/dumntot
	if (dumvbar_aft .gt. volumhi_sect(n,itype)) then
	    do while ( (nnew .lt. nsize_aer(itype)) .and.   &
      		       (dumvbar_aft .gt. volumhi_sect(nnew,itype)) )
		nnew = nnew + 1
	    end do

	else if (dumvbar_aft .lt. volumlo_sect(n,itype)) then
	    do while ( (nnew .gt. 1) .and.   &
      		       (dumvbar_aft .lt. volumlo_sect(nnew,itype)) )
		nnew = nnew - 1
	    end do

	end if

1290	nnewsave(1,n) = nnew
	nnewsave(2,n) = 0

	if (iforce_movecenter(n) .gt. 0) goto 3700








	if ( (drydenspp(n) .lt.  0.1) .or.   &
	     (drydenspp(n) .gt. 20.0) .or.   &
      	     (drymasspp(n) .le. smallmassaa) ) then
	    iforce_movecenter(n) = 11
	    goto 3700
	end if

	dumvtot_pre = drymasspp(n)/drydenspp(n)

	dumvlo_pre = volumlo_sect(n,itype)
	dumvhi_pre = volumhi_sect(n,itype)
	dumvdel_pre = dumvhi_pre - dumvlo_pre



	dumv = dumvhi_pre - 0.01*dumvdel_pre
	if (dumntot .le. dumvtot_pre/dumv) then
	    iforce_movecenter(n) = 12
	    goto 3700
	end if
	dumv = dumvlo_pre + 0.01*dumvdel_pre
	if (dumntot .ge. dumvtot_pre/dumv) then
	    iforce_movecenter(n) = 13
	    goto 3700
	end if


	dumvbar_pre = dumvtot_pre/dumntot
	dumgamma = (dumvhi_pre/dumvlo_pre) - 1.0
	dumratio = dumvbar_pre/dumvlo_pre

	if (dumratio .le. (1.0001 + dumgamma/3.0)) then
	    dumv = dumvlo_pre + 3.0*(dumvbar_pre-dumvlo_pre)
	    dumvhi_pre = min( dumvhi_pre, dumv )
	    dumvdel_pre = dumvhi_pre - dumvlo_pre
	    dumgamma = (dumvhi_pre/dumvlo_pre) - 1.0
	    dumratio = dumvbar_pre/dumvlo_pre
	else if (dumratio .ge. (0.9999 + dumgamma*2.0/3.0)) then
	    dumv = dumvhi_pre + 3.0*(dumvbar_pre-dumvhi_pre)
	    dumvlo_pre = max( dumvlo_pre, dumv )
	    dumvdel_pre = dumvhi_pre - dumvlo_pre
	    dumgamma = (dumvhi_pre/dumvlo_pre) - 1.0
	    dumratio = dumvbar_pre/dumvlo_pre
	end if

	dumbb = (dumratio - 1.0 - 0.5*dumgamma)*12.0/dumgamma
	dumaa = 1.0 - 0.5*dumbb



	dumvcutlo_nnew_pre = volumlo_sect(nnew,itype)*(dumvbar_pre/dumvbar_aft)
	dumvcuthi_nnew_pre = volumhi_sect(nnew,itype)*(dumvbar_pre/dumvbar_aft)




	if (nnew .eq. 1) then
	    if (dumvhi_pre .le. dumvcuthi_nnew_pre) then
		iforce_movecenter(n) = 21
	    else
		nnew2 = nnew + 1
	    end if
	else if (nnew .eq. nsize_aer(itype)) then
	    if (dumvlo_pre .ge. dumvcutlo_nnew_pre) then
		iforce_movecenter(n) = 22
	    else
		nnew2 = nnew - 1
	    end if
	else
	    if ((dumvlo_pre .ge. dumvcutlo_nnew_pre) .and.   &
      		(dumvhi_pre .le. dumvcuthi_nnew_pre)) then
		iforce_movecenter(n) = 23
	    else if (dumvlo_pre .lt. dumvcutlo_nnew_pre) then
		nnew2 = nnew - 1
	    else
		nnew2 = nnew + 1
	    end if
	end if
	if (iforce_movecenter(n) .gt. 0) goto 3700



	dumzlo = (dumvcutlo_nnew_pre - dumvlo_pre)/dumvdel_pre
	dumzhi = (dumvcuthi_nnew_pre - dumvlo_pre)/dumvdel_pre
	dumzlo = max( dumzlo, 0.0 )
	dumzhi = min( dumzhi, 1.0 )
	dum1 =  dumzhi    - dumzlo
	dum2 = (dumzhi**2 - dumzlo**2)*0.5
	dum3 = (dumzhi**3 - dumzlo**3)/3.0
	dumfracnum = dumaa*dum1 + dumbb*dum2
	dumfracvol = (dumvlo_pre/dumvbar_pre) * (dumaa*dum1 +   &
      		(dumaa*dumgamma + dumbb)*dum2 + (dumbb*dumgamma)*dum3)

	if ((dumfracnum .le. 0.0) .or. (dumfracvol .le. 0.0)) then
	    iforce_movecenter(n) = 31
	    nnewsave(1,n) = nnew2
	else if ((dumfracnum .ge. 1.0) .or. (dumfracvol .ge. 1.0)) then
	    iforce_movecenter(n) = 32
	end if
	if (iforce_movecenter(n) .gt. 0) goto 3700

	nnewsave(2,n) = nnew2

	xferfracvol(1,n) = dumfracvol
	xferfracvol(2,n) = 1.0 - dumfracvol
	xferfracnum(1,n) = dumfracnum
	xferfracnum(2,n) = 1.0 - dumfracnum

3700	continue


	if (idiag_movesect .ne. 70) goto 3800

	if (nnewsave(2,n) .eq. 0) then
	    if (nnewsave(1,n) .eq. 0) then
		dumch4 = 'NO X'
	    else if (nnewsave(1,n) .eq. n) then
		dumch4 = 'NO A'
	    else
		dumch4 = 'YESA'
	    end if
	else if (nnewsave(1,n) .eq. 0) then
	    if (nnewsave(2,n) .eq. n) then
		dumch4 = 'NO B'
	    else
		dumch4 = 'YESB'
	    end if
	else if (nnewsave(2,n) .eq. n) then
	    if (nnewsave(1,n) .eq. n) then
		dumch4 = 'NO Y'
	    else
		dumch4 = 'YESC'
	    end if
	else if (nnewsave(1,n) .eq. n) then
	    dumch4 = 'YESD'
	else
	    dumch4 = 'YESE'
	end if

	dumch1 = '+'
	if (drymasspp(n) .gt. drymassxx(n)) dumch1 = '-'
		
	msg = ' '
	call peg_message( lunout, msg )
	write(msg,97010) dumch1, dumch4, iclm, jclm, k, m,   &
      		n, nnewsave(1,n), nnewsave(2,n), iforce_movecenter(n)
	call peg_message( lunout, msg )
	write(msg,97020) 'pre mass, dens      ',   &
      		drymasspp(n), drydenspp(n)
	call peg_message( lunout, msg )
	write(msg,97020) 'aft mass, dens, numb',   &
      		drymassxx(n), drydensxx(n), rnumbxx(n)
	call peg_message( lunout, msg )
	if ((drydensxx(n) .ne. drydensxx0(n)) .or.   &
      	    (drymassxx(n) .ne. drymassxx0(n)) .or.   &
      	    (rnumbxx(n)   .ne. rnumbxx0(n)  )) then
      	    write(msg,97020) 'aft0 mas, dens, numb',   &
      		drymassxx0(n), drydensxx0(n), rnumbxx0(n)
	    call peg_message( lunout, msg )
	end if
	write(msg,97020) 'vlop0, vbarp,  vhip0',   &
      		volumlo_sect(n,itype), dumvbar_pre, volumhi_sect(n,itype)
	call peg_message( lunout, msg )
	write(msg,97020) 'vlop , vbarp,  vhip ',   &
      		dumvlo_pre, dumvbar_pre, dumvhi_pre
	call peg_message( lunout, msg )
	write(msg,97020) 'vloax, vbarax, vhiax',   &
      		dumvcutlo_nnew_pre, dumvbar_pre, dumvcuthi_nnew_pre
	call peg_message( lunout, msg )
	write(msg,97020) 'vloa0, vbara,  vhia0',   &
      		volumlo_sect(nnew,itype), dumvbar_aft, volumhi_sect(nnew,itype)
	call peg_message( lunout, msg )
	write(msg,97020) 'dumfrvol, num, ratio',   &
      		dumfracvol, dumfracnum, dumratio
	call peg_message( lunout, msg )
	write(msg,97020) 'frvol,num1; vol,num2',   &
      		xferfracvol(1,n), xferfracnum(1,n),   &
      		xferfracvol(2,n), xferfracnum(2,n)
	call peg_message( lunout, msg )

97010	format( 'movesect', 2a, 7x, 4i3, 4x,   &
      		'n,nnews', 3i3, 4x, 'iforce', i3.2 )
97020	format( a, 1p, 4e13.4 )

3800	continue






	ierr = 0
	if (nnewsave(1,n) .eq. nnewsave(2,n)) then
	    ierr = 1
	else if (nnewsave(1,n)*nnewsave(2,n) .ne. 0) then
	    if (iabs(nnewsave(1,n)-nnewsave(2,n)) .ne. 1) ierr = 1
	end if
	if (ierr .gt. 0) then
	    write(msg,97010) 'E', 'RROR', iclm, jclm, k, m,   &
      		n, nnewsave(1,n), nnewsave(2,n), iforce_movecenter(n)
	    call peg_message( lunout, msg )
	end if




	if ((method_movesect .ge. 30) .and. (method_movesect .le. 39)) then
	    nnewsave(1,n) = nnew
	    nnewsave(2,n) = 0
	    xferfracvol(1,n) = 1.0
	    xferfracvol(2,n) = 0.0
	    xferfracnum(1,n) = 1.0
	    xferfracnum(2,n) = 0.0
	end if

3900	continue

	return
	end subroutine move_sections_calc_masnumadvect                          



	subroutine move_sections_apply_moves(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy,   &
      	  xferfracvol, xferfracnum )




	implicit none








	integer iflag, iclm, jclm, iphase, itype, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize)

	real densdefault, densh2o, smallmassaa, smallmassbb
	real delta_water_conform1, delta_numb_conform1
	real drydenspp(maxd_asize)
	real drymasspp(maxd_asize)
	real drymassxx(maxd_asize), drymassyy(maxd_asize)
	real dryvolxx(maxd_asize), dryvolyy(maxd_asize)
	real rmassxx(maxd_acomp+2,maxd_asize),   &
      	     rmassyy(maxd_acomp+2,maxd_asize)
	real rnumbpp(maxd_asize)
	real rnumbxx(maxd_asize), rnumbyy(maxd_asize)
	real xferfracvol(2,maxd_asize), xferfracnum(2,maxd_asize)
	real wetvolxx(maxd_asize), wetvolyy(maxd_asize)
	real wetmassxx(maxd_asize), wetmassyy(maxd_asize)


	integer jj, l, ll, n, ndum, nnew, nold
	integer jja, jjb, jjc

	real delta_numb_conform2,   &
	  dumbot, dumnum, dumnum_at_dhi, dumnum_at_dlo,   &
      	  dumvol, dumvol1p, dumxfvol, dumxfnum, sixoverpi, third
	real dumpp(maxd_asize), dumxx(maxd_asize), dumyy(maxd_asize),   &
	  dumout(maxd_asize)

	character*160 msg
	character*8 dumch8
	character*4 dumch4


	sixoverpi = 6.0/pi
	third = 1.0/3.0






	do 1900 n = 1, nsize_aer(itype)

	if ( (nnewsave(1,n) .eq. n) .and.   &
      	     (nnewsave(2,n) .eq. 0) ) then




	    drymassyy(n) = drymassxx(n)
	    dryvolyy(n) = dryvolxx(n)
	    wetmassyy(n) = wetmassxx(n)
	    wetvolyy(n) = wetvolxx(n)
	    rnumbyy(n) = rnumbxx(n)
	    do ll = 1, ncomp_plustracer_aer(itype) + 2
		rmassyy(ll,n) = rmassxx(ll,n)
	    end do

	else




	    drymassyy(n) = 0.0
	    dryvolyy(n) = 0.0
	    wetmassyy(n) = 0.0
	    wetvolyy(n) = 0.0
	    rnumbyy(n) = 0.0
	    do ll = 1, ncomp_plustracer_aer(itype) + 2
		rmassyy(ll,n) = 0.0
	    end do

	end if

1900	continue




	do 2900 nold = 1, nsize_aer(itype)

	if ( (nnewsave(1,nold) .eq. nold) .and.   &
      	     (nnewsave(2,nold) .eq. 0   ) ) goto 2900

	do 2800 jj = 1, 2

	nnew = nnewsave(jj,nold)
	if (nnew .le. 0) goto 2800

	dumxfvol = xferfracvol(jj,nold)
	dumxfnum = xferfracnum(jj,nold)

	do ll = 1, ncomp_plustracer_aer(itype) + 2
	    rmassyy(ll,nnew) = rmassyy(ll,nnew) + rmassxx(ll,nold)*dumxfvol
	end do
	rnumbyy(nnew) = rnumbyy(nnew) + rnumbxx(nold)*dumxfnum

	drymassyy(nnew) = drymassyy(nnew) + drymassxx(nold)*dumxfvol
	dryvolyy(nnew)  = dryvolyy(nnew)  + dryvolxx(nold)*dumxfvol
	wetmassyy(nnew) = wetmassyy(nnew) + wetmassxx(nold)*dumxfvol
	wetvolyy(nnew)  = wetvolyy(nnew)  + wetvolxx(nold)*dumxfvol

2800	continue

2900	continue









	call move_sections_conserve_check(   &
	  1, iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy )

	delta_numb_conform2 = 0.0

	do 3900 n = 1, nsize_aer(itype)

	dumvol = dryvolyy(n)
	if (min(drymassyy(n),dumvol) .le. smallmassbb) then
	    dumvol = drymassyy(n)/densdefault
	    dumnum = drymassyy(n)/(volumlo_sect(n,itype)*densdefault)
	    delta_numb_conform2 = delta_numb_conform2 + dumnum - rnumbyy(n)
	    rnumbyy(n) = dumnum
	    adrydens_sub(n,itype,k,m) = densdefault
	    awetdens_sub(n,itype,k,m) = densdefault
	    admeandry_sub(n,itype,k,m) = dcen_sect(n,itype)
	    admeanwet_sub(n,itype,k,m) = dcen_sect(n,itype)
	else
	    dumnum = rnumbyy(n)
	    dumnum_at_dhi = dumvol/volumhi_sect(n,itype)
	    dumnum_at_dlo = dumvol/volumlo_sect(n,itype)
	    dumnum = max( dumnum_at_dhi, min( dumnum_at_dlo, dumnum ) )
	    delta_numb_conform2 = delta_numb_conform2 + dumnum - rnumbyy(n)
	    rnumbyy(n) = dumnum
	    adrydens_sub(n,itype,k,m) = drymassyy(n)/dumvol
	    dumvol1p = dumvol/dumnum
	    admeandry_sub(n,itype,k,m) = (dumvol1p*sixoverpi)**third
	    awetdens_sub(n,itype,k,m) = wetmassyy(n)/wetvolyy(n)
	    dumvol1p = wetvolyy(n)/dumnum
	    admeanwet_sub(n,itype,k,m) = min( 100.*dcen_sect(n,itype),   &
      			(dumvol1p*sixoverpi)**third )
	end if
	aqvoldry_sub(n,itype,k,m) = dumvol
	aqmassdry_sub(n,itype,k,m) = drymassyy(n)






	do ll = 1, ncomp_plustracer_aer(itype)
	    l = massptr_aer(ll,n,itype,iphase)
	    rsub(l,k,m) = rmassyy(ll,n)
	end do
	l = 0
	if (iphase .eq. ai_phase) then
	    l = waterptr_aer(n,itype)
	    if (l .gt. 0) rsub(l,k,m) = rmassyy(llwater,n)
	    l = hyswptr_aer(n,itype)
	    if (l .gt. 0) rsub(l,k,m) = rmassyy(llhysw,n)
	end if
	rsub(numptr_aer(n,itype,iphase),k,m) = rnumbyy(n)

3900	continue

	delta_numb_conform1 = delta_numb_conform1 + delta_numb_conform2

	call move_sections_conserve_check(   &
	  2, iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy )




	if (idiag_movesect .ne. 70) goto 4900

	ndum = 0
	do n = 1, nsize_aer(itype)
	    if (nnewsave(1,n)+nnewsave(2,n) .ne. n) ndum = ndum + 1
	end do










	dumch4 = 'NONE'
	if (ndum .gt. 0) dumch4 = 'SOME'
	msg = ' '
	call peg_message( lunout, msg )
	write(msg,97010) dumch4, iclm, jclm, k, m, ndum
	call peg_message( lunout, msg )
	do jjb = 1, nsize_aer(itype), 10
	    jjc = min( jjb+9, nsize_aer(itype) )
	    write(msg,97011) (nnewsave(1,n), nnewsave(2,n), n=jjb,jjc)
	    call peg_message( lunout, msg )
	end do











	dumbot = log( volumhi_sect(1,itype)/volumlo_sect(1,itype) )
	do n = 1, nsize_aer(itype)
	    dumpp(n) = -9.99
	    dumxx(n) = -9.99
	    dumyy(n) = -9.99
	    if ( (drydenspp(n) .gt. 0.5) .and.   &
      	         (drymasspp(n) .gt. smallmassaa) ) then
      		dumvol = drymasspp(n)/drydenspp(n)
		if ((rnumbpp(n) .ge. 1.0e-35) .and.   &
      		    (dumvol .ge. 1.0e-35)) then
		    dumvol1p = dumvol/rnumbpp(n)
		    dumpp(n) = 1.0 + log(dumvol1p/volumlo_sect(1,itype))/dumbot
		end if
	    end if
	    if ((rnumbxx(n) .ge. 1.0e-35) .and.   &
      		(dryvolxx(n) .ge. 1.0e-35)) then
		dumvol1p = dryvolxx(n)/rnumbxx(n)
		dumxx(n) = 1.0 + log(dumvol1p/volumlo_sect(1,itype))/dumbot
	    end if
	    if ((rnumbyy(n) .ge. 1.0e-35) .and.   &
      		(dryvolyy(n) .ge. 1.0e-35)) then
		dumvol1p = dryvolyy(n)/rnumbyy(n)
		dumyy(n) = 1.0 + log(dumvol1p/volumlo_sect(1,itype))/dumbot
	    end if
	end do






	do jja = 1, 7
	    if      (jja .eq. 1) then
		dumch8 = 'rnumbold'
		dumout(:) = rnumbxx(:)
	    else if (jja .eq. 2) then
		dumch8 = 'rnumbnew'
		dumout(:) = rnumbyy(:)
	    else if (jja .eq. 3) then
		dumch8 = 'drvolold'
		dumout(:) = dryvolxx(:)
	    else if (jja .eq. 4) then
		dumch8 = 'drvolnew'
		dumout(:) = dryvolyy(:)
	    else if (jja .eq. 5) then
		dumch8 = 'lnvolold'
		dumout(:) = dumxx(:)
	    else if (jja .eq. 6) then
		dumch8 = 'lnvolnew'
		dumout(:) = dumyy(:)
	    else if (jja .eq. 7) then
		dumch8 = 'lnvolpre'
		dumout(:) = dumpp(:)
	    end if
	    do jjb = 1, nsize_aer(itype), 10
		jjc = min( jjb+9, nsize_aer(itype) )
		if (jja .le. 4) then
		    write(msg,97020) dumch8, (dumout(n), n=jjb,jjc)
		else
		    write(msg,97030) dumch8, (dumout(n), n=jjb,jjc)
		end if
		call peg_message( lunout, msg )
		dumch8 = ' '
	    end do
	end do




97010	format( 'movesectapply', a, 4i3, 3x, i3 )
97011	format( 5x, 10(3x,2i3) )
97020	format( a, 1p, 10e9.1 )
97030	format( a,     10f9.3 )

4900	continue
	return
	end subroutine move_sections_apply_moves                          





	subroutine move_sections_apply_n1_inflow(   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy,   &
      	  xferfracvol, xferfracnum,   &
      	  specmwxx, specdensxx )









	implicit none








	integer iflag, iclm, jclm, iphase, itype, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize)

	real densdefault, densh2o, smallmassaa, smallmassbb
	real delta_water_conform1, delta_numb_conform1
	real drymassxx(maxd_asize), drymassyy(maxd_asize)
	real dryvolxx(maxd_asize), dryvolyy(maxd_asize)
	real rmassxx(maxd_acomp+2,maxd_asize),   &
      	     rmassyy(maxd_acomp+2,maxd_asize)
	real rnumbxx(maxd_asize), rnumbyy(maxd_asize)
	real xferfracvol(2,maxd_asize), xferfracnum(2,maxd_asize)
	real wetvolxx(maxd_asize), wetvolyy(maxd_asize)
	real wetmassxx(maxd_asize), wetmassyy(maxd_asize)
	real specdensxx(maxd_acomp), specmwxx(maxd_acomp)



	integer jj, l, ll, n, nnew, nold

	real deltanum, deltavol, dumvol1p, dumxfvol, dumxfnum





	nold = 1
	n = nold
	if ( (nnewsave(1,nold) .eq. nold) .and.   &
      	     (nnewsave(2,nold) .eq. 0   ) ) goto 3900

	dumxfnum = 0.0
	do 2800 jj = 1, 2
	    nnew = nnewsave(jj,nold)
	    if (nnew .le. 0) goto 2800
	    if (nnew .eq. nold) goto 2800
	    dumxfnum = dumxfnum + xferfracnum(jj,nold)
2800	continue






	dumvol1p = 0.95*volumlo_sect(n,itype) + 0.05*volumhi_sect(n,itype)
	deltanum = dumxfnum*rnumbxx(n)
	deltavol = deltanum*dumvol1p
	if (dumxfnum .le. 0.0) goto 3900
	if (deltanum .le. 0.0) goto 3900
	if (deltavol .le. 0.0) goto 3900







	if (dryvolxx(n) .gt. smallmassbb) then
	    dumxfvol = deltavol/dryvolxx(n)
	    do ll = 1, ncomp_plustracer_aer(itype) + 2
		rmassyy(ll,n) = rmassyy(ll,n) + dumxfvol*rmassxx(ll,n)
	    end do
	else
	    ll = 1
	    rmassyy(ll,n) = rmassyy(ll,n) + deltavol*specdensxx(ll)/specmwxx(ll)
	end if
	rnumbyy(n) = rnumbyy(n) + deltanum




	do ll = 1, ncomp_plustracer_aer(itype)
	    l = massptr_aer(ll,n,itype,iphase)
	    rsub(l,k,m) = rmassyy(ll,n)
	end do
	if (iphase .eq. ai_phase) then
	    l = waterptr_aer(n,itype)
	    if (l .gt. 0) rsub(l,k,m) = rmassyy(llwater,n)
	    l = hyswptr_aer(n,itype)
	    if (l .gt. 0) rsub(l,k,m) = rmassyy(llhysw,n)
	end if
	rsub(numptr_aer(n,itype,iphase),k,m) = rnumbyy(n)


3900	continue
	return
	end subroutine move_sections_apply_n1_inflow



	subroutine move_sections_conserve_check( ipass,   &
	  iflag, iclm, jclm, k, m, iphase, itype,   &
      	  method_movesect, idiag_movesect, llhysw, llwater, nnewsave,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy )
















	implicit none








	integer ipass, iflag, iclm, jclm, iphase, itype, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize)
	real densdefault, densh2o, smallmassaa, smallmassbb
	real delta_water_conform1, delta_numb_conform1
	real drymassxx(maxd_asize), drymassyy(maxd_asize)
	real dryvolxx(maxd_asize), dryvolyy(maxd_asize)
	real rmassxx(maxd_acomp+2,maxd_asize),   &
      	     rmassyy(maxd_acomp+2,maxd_asize)
	real rnumbxx(maxd_asize), rnumbyy(maxd_asize)
	real wetvolxx(maxd_asize), wetvolyy(maxd_asize)
	real wetmassxx(maxd_asize), wetmassyy(maxd_asize)


	integer jj, l, ll, llworst, llworstb, n
	integer nerr, nerrmax
	save nerr, nerrmax
	data nerr, nerrmax / 0, 999 /

	real dumbot, dumtop, dumtoler, dumerr, dumworst, dumworstb
	real duma, dumb, dumc, dume
	real thesum(4,maxd_acomp+7)
	save thesum

	character*8 dumname(maxd_acomp+7)
	character*160 msg


	if (ipass .eq. 2) goto 2000

	do ll = 1, ncomp_plustracer_aer(itype)+7
	do jj = 1, 4
	    thesum(jj,ll) = 0.0
	end do
	end do

	do n = 1, nsize_aer(itype)
	    do ll = 1, ncomp_plustracer_aer(itype)+2
		thesum(1,ll) = thesum(1,ll) + rmassxx(ll,n)
		thesum(2,ll) = thesum(2,ll) + rmassyy(ll,n)
	    end do
	    ll = ncomp_plustracer_aer(itype)+3
	    thesum(1,ll) = thesum(1,ll) + rnumbxx(n)
	    thesum(2,ll) = thesum(2,ll) + rnumbyy(n)
	    ll = ncomp_plustracer_aer(itype)+4
	    thesum(1,ll) = thesum(1,ll) + drymassxx(n)
	    thesum(2,ll) = thesum(2,ll) + drymassyy(n)
	    ll = ncomp_plustracer_aer(itype)+5
	    thesum(1,ll) = thesum(1,ll) + dryvolxx(n)
	    thesum(2,ll) = thesum(2,ll) + dryvolyy(n)
	    ll = ncomp_plustracer_aer(itype)+6
	    thesum(1,ll) = thesum(1,ll) + wetmassxx(n)
	    thesum(2,ll) = thesum(2,ll) + wetmassyy(n)
	    ll = ncomp_plustracer_aer(itype)+7
	    thesum(1,ll) = thesum(1,ll) + wetvolxx(n)
	    thesum(2,ll) = thesum(2,ll) + wetvolyy(n)
	end do


2000	continue




	do n = 1, nsize_aer(itype)
	    do ll = 1, ncomp_plustracer_aer(itype)+3
		if (ll .le. ncomp_plustracer_aer(itype)) then
		    l = massptr_aer(ll,n,itype,iphase)
		else if (ll .eq. ncomp_plustracer_aer(itype)+1) then
		    l = 0
		    if (iphase .eq. ai_phase) l = hyswptr_aer(n,itype)
		else if (ll .eq. ncomp_plustracer_aer(itype)+2) then
		    l = 0
		    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
		else
		    l = numptr_aer(n,itype,iphase)
		end if
		if (l .gt. 0)   &
		    thesum(ipass+2,ll) = thesum(ipass+2,ll) + rsub(l,k,m)
	    end do
	end do
	if (ipass .eq. 2) then
	    ll = ncomp_plustracer_aer(itype)+2
	    thesum(3,ll) = thesum(3,ll) + delta_water_conform1
	    ll = ncomp_plustracer_aer(itype)+3
	    thesum(3,ll) = thesum(3,ll) + delta_numb_conform1
	end if






	do ll = 1, ncomp_plustracer_aer(itype)+7
	    dumname(ll) = ' '
	    write(dumname(ll),'(i4.4)') ll
	    if (ll .le. ncomp_plustracer_aer(itype)) dumname(ll) =   &
      	    		name(massptr_aer(ll,1,itype,iphase))(1:4)
	    if (ll .eq. ncomp_plustracer_aer(itype)+1) dumname(ll) = 'hysw'
	    if (ll .eq. ncomp_plustracer_aer(itype)+2) dumname(ll) = 'watr'
	    if (ll .eq. ncomp_plustracer_aer(itype)+3) dumname(ll) = 'numb'
	    if (ll .eq. ncomp_plustracer_aer(itype)+4) dumname(ll) = 'drymass'
	    if (ll .eq. ncomp_plustracer_aer(itype)+5) dumname(ll) = 'dryvol'
	    if (ll .eq. ncomp_plustracer_aer(itype)+6) dumname(ll) = 'wetmass'
	    if (ll .eq. ncomp_plustracer_aer(itype)+7) dumname(ll) = 'wetvol'
	end do

	jj = 2*ipass - 1
	dumworst = 0.0
	dumworstb = 0.0
	llworst = 0
	llworstb = 0
	do ll = 1, ncomp_plustracer_aer(itype)+7
	    dumtop = thesum(jj+1,ll) - thesum(jj,ll)


	    dumbot = max( abs(thesum(jj,ll)), abs(thesum(jj+1,ll)), 1.0e-30 )
	    dumerr = dumtop/dumbot







	    if (ipass .eq. 2) then
		dumc = 1.0
		if (ll .eq. ncomp_plustracer_aer(itype)+2) then
		    dumc = delta_water_conform1
		else if (ll .eq. ncomp_plustracer_aer(itype)+3) then
		    dumc = delta_numb_conform1
		end if
		if (dumc .lt. 0.0) then
		    duma = thesum(3,ll) - dumc
		    dumb = thesum(4,ll) - dumc
		    dumtop = dumb - duma

		    dumbot = max( abs(duma), abs(dumb), 1.0e-30 )
		    dume = dumtop/dumbot
		    if (abs(dume) .lt. abs(dumerr)) dumerr = dume
		end if
	    end if

	    if (abs(dumerr) .gt. abs(dumworst)) then
		llworstb = llworst
		dumworstb = dumworst
		llworst = ll
		dumworst = dumerr
	    end if
	end do

	dumtoler = 1.0e-6
	if (abs(dumworst) .gt. dumtoler) then
	    nerr = nerr + 1
	    if (nerr .le. nerrmax) then
		msg = ' '
		call peg_message( lunout, msg )
		write(msg,97110) iclm, jclm, k, m, ipass, llworst
		call peg_message( lunout, msg )
		write(msg,97120) '    nnew(1,n)',   &
      			(nnewsave(1,n), n=1,nsize_aer(itype))
		call peg_message( lunout, msg )
		write(msg,97120) '    nnew(2,n)',   &
      			(nnewsave(2,n), n=1,nsize_aer(itype))
		call peg_message( lunout, msg )

		ll = llworst
		if (ll .eq. 0) ll = ncomp_plustracer_aer(itype)+2
		write(msg,97130) 'name/relerr/thesum', jj, '/thesum', jj+1,   &
      			dumname(ll), dumworst, thesum(jj,ll), thesum(jj+1,ll)
		call peg_message( lunout, msg )

		if ( (ll .eq. ncomp_plustracer_aer(itype)+3) .and.   &
		     (abs(dumworstb) .gt. dumtoler) ) then
		    ll = max( 1, llworstb )
		    dumtop = thesum(jj+1,ll) - thesum(jj,ll)

		    dumbot = max( abs(thesum(jj,ll)), abs(thesum(jj+1,ll)), 1.0e-30 )
		    dumerr = dumtop/dumbot
		    write(msg,97130) 'name/relerr/thesum', jj, '/thesum', jj+1,   &
      			dumname(ll), dumerr, thesum(jj,ll), thesum(jj+1,ll)
		    call peg_message( lunout, msg )
		end if
	    end if
	end if

97110	format( 'movesect conserve ERROR - i/j/k/m/pass/llworst',   &
		4i3, 2x, 2i3 )
97120	format( a, 64i3 )
97130	format( a, i1, a, i1, 2x, a, 1p, 3e16.7 )

	return
	end subroutine move_sections_conserve_check        



	subroutine test_move_sections( iflag_test, iclm, jclm, k, m )




	implicit none








	integer iflag_test, iclm, jclm, k, m


	integer idiag_movesect, iflag, ii, iphase, itype, jj,   &
	  l, ll, n, nn
	integer ientryno
	save ientryno
	data ientryno / 0 /

	real dumnumb, dumvolpre, dumvolaft
	real dumsv_rsub(l2maxd)
	real dumsv_drymass_pregrow(maxd_asize,maxd_atype)
	real dumsv_drymass_aftgrow(maxd_asize,maxd_atype)
	real dumsv_drydens_pregrow(maxd_asize,maxd_atype)
	real dumsv_drydens_aftgrow(maxd_asize,maxd_atype)

	character*160 msg

	integer maxvolfactpre, maxvolfactaft
	parameter (maxvolfactpre=15, maxvolfactaft=23)

	real dumvolfactpre(maxvolfactpre)
	data dumvolfactpre /   &
      	2.0, 0.0, 1.0e-20, 0.5, 0.9,   &
      	1.0, 1.01, 1.1, 2.0, 4.0, 7.9, 7.99, 8.0,   &
      	8.1, 16.0 /

	real dumvolfactaft(maxvolfactaft)
	data dumvolfactaft /   &
      	4.0, 0.0, 1.0e-20, 0.01, 0.02, 0.05, 0.1, 0.5, 0.9,   &
      	1.0, 1.01, 1.1, 2.0, 4.0, 7.9, 7.99, 8.0,   &
      	8.1, 16.0, 32., 64., 128., 256. /






	if (msectional .le. 0) return
	if (nsize_aer(1) .le. 0) return

	idiag_movesect = mod( msectional, 10000 )/100
	if (idiag_movesect .ne. 70) return
	
	ientryno = ientryno + 1
	if (ientryno .gt. 2) return




	do l = 1, ltot2
	    dumsv_rsub(l) = rsub(l,k,m)
	end do
	do itype = 1, ntype_aer
	    do n = 1, nsize_aer(itype)
		dumsv_drymass_pregrow(n,itype) = drymass_pregrow(n,itype)
		dumsv_drymass_aftgrow(n,itype) = drymass_aftgrow(n,itype)
		dumsv_drydens_pregrow(n,itype) = drydens_pregrow(n,itype)
		dumsv_drydens_aftgrow(n,itype) = drydens_aftgrow(n,itype)
	    end do
	end do




	do 3900 iflag = 1, 1

	iphase = ai_phase
	if (iabs(iflag) .eq. 2) iphase = cw_phase

	do 3800 itype = 1, ntype_aer

	do 2900 nn = 1, nsize_aer(itype)

	do 2800 ii = 1, maxvolfactpre

	do 2700 jj = 1, maxvolfactaft


	do n = 1, nsize_aer(itype)
	    do ll = 1, ncomp_plustracer_aer(itype)
		rsub(massptr_aer(ll,n,itype,iphase),k,m) = 0.0
	    end do
	    l = 0
	    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
	    if (l .gt. 0) rsub(l,k,m) = 0.0
	    l = numptr_aer(n,itype,iphase)
	    if (l .gt. 0) rsub(l,k,m) = 0.0
	    drymass_pregrow(n,itype) = 0.0
	    drymass_aftgrow(n,itype) = 0.0
	    drydens_pregrow(n,itype) = -1.0
	    drydens_aftgrow(n,itype) = -1.0
	end do


	n = nn
	dumnumb = 1.0e7
	rsub(numptr_aer(n,itype,iphase),k,m) = dumnumb
	ll = 1
	l = massptr_aer(ll,n,itype,iphase)

	dumvolpre = volumlo_sect(n,itype)*dumvolfactpre(ii)*dumnumb
	drydens_pregrow(n,itype) = dens_aer(ll,itype)
	drymass_pregrow(n,itype) = dumvolpre*drydens_pregrow(n,itype)
	if (ii .eq. 1) drydens_pregrow(n,itype) = -1.0
	
	dumvolaft = volumlo_sect(n,itype)*dumvolfactaft(jj)*dumnumb
	drydens_aftgrow(n,itype) = dens_aer(ll,itype)
	drymass_aftgrow(n,itype) = dumvolaft*drydens_aftgrow(n,itype)
	if (jj .eq. 1) drydens_aftgrow(n,itype) = -1.0
	
	rsub(l,k,m) = drymass_aftgrow(n,itype)/mw_aer(ll,itype)
	
	msg = ' '
	call peg_message( lunout, msg )
	write(msg,98010) nn, ii, jj
	call peg_message( lunout, msg )
	write(msg,98011) dumvolfactpre(ii), dumvolfactaft(jj)
	call peg_message( lunout, msg )


	call move_sections( iflag, iclm, jclm, k, m )

	msg = ' '
	call peg_message( lunout, msg )
	write(msg,98010) nn, ii, jj
	call peg_message( lunout, msg )
	write(msg,98011) dumvolfactpre(ii), dumvolfactaft(jj)
	call peg_message( lunout, msg )

2700	continue
2800	continue
2900	continue

3800	continue
3900	continue

98010	format( 'test_move_sections output - nn, ii, jj =', 3i3 )
98011	format( 'volfactpre, volfactaft =', 1p, 2e12.4 )




	do l = 1, ltot2
	    rsub(l,k,m) = dumsv_rsub(l)
	end do
	do itype = 1, ntype_aer
	do n = 1, nsize_aer(itype)
	    drymass_pregrow(n,itype) = dumsv_drymass_pregrow(n,itype)
	    drymass_aftgrow(n,itype) = dumsv_drymass_aftgrow(n,itype)
	    drydens_pregrow(n,itype) = dumsv_drydens_pregrow(n,itype)
	    drydens_aftgrow(n,itype) = dumsv_drydens_aftgrow(n,itype)
	end do
	end do

	return
	end subroutine test_move_sections                           



	subroutine move_sections_checkptrs( iflag, iclm, jclm, k, m )



	implicit none








	integer iflag, iclm, jclm, k, m


	integer l, itype, iphase, n, ndum
	character*160 msg

	do 1900 itype = 1, ntype_aer
	do 1800 iphase = 1, nphase_aer

	ndum = 0
	do n = 1, nsize_aer(itype)
	    l = numptr_aer(n,itype,iphase)
	    if ((l .le. 0) .or. (l .gt. ltot2)) then
		msg = '*** subr move_sections error - ' //   &
			'numptr_amode not defined'
		call peg_message( lunerr, msg )
		write(msg,9030) 'mode, numptr =', n, l
		call peg_message( lunerr, msg )
		write(msg,9030) 'iphase, itype =', iphase, itype
		call peg_message( lunerr, msg )
		call peg_error_fatal( lunerr, msg )
	    end if


	    l = 0
	    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
	    if ((l .gt. 0) .and. (l .le. ltot2)) ndum = ndum + 1
	end do
	if ((ndum .ne. 0) .and. (ndum .ne. nsize_aer(itype))) then
	    msg = '*** subr move_sections error - ' //   &
      		'waterptr_aer must be on/off for all modes'
	    call peg_message( lunerr, msg )
	    write(msg,9030) 'iphase, itype =', iphase, itype
	    call peg_message( lunerr, msg )
	    call peg_error_fatal( lunerr, msg )
	end if
9030	format( a, 2(1x,i6) )

1800	continue
1900	continue

	return
	end subroutine move_sections_checkptrs                           



	end module module_mosaic_movesect
