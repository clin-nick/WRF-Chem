
















	module module_mosaic_initmixrats







	USE module_state_description

	integer, parameter :: mosaic_init_wrf_mixrats_flagaa = 1
                            

	contains



	subroutine mosaic_init_wrf_mixrats(  &
		iflagaa, config_flags,       &
		chem, alt, z_at_w, g,        &
		ids,ide, jds,jde, kds,kde,   &
		ims,ime, jms,jme, kms,kme,   &
		its,ite, jts,jte, kts,kte    )







	use module_configure, only:  grid_config_rec_type
	use module_state_description, only:  num_chem, param_first_scalar,   &
			aer_ic_pnnl
	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_peg_util, only:  peg_message, peg_error_fatal

	implicit none



	type(grid_config_rec_type), intent(in) :: config_flags

	integer, intent(in) ::   &
		iflagaa,   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem

	real, intent(in),      &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
        alt, z_at_w
    real, intent(in) :: g


        integer :: i, ic, j, k


	return
	end subroutine mosaic_init_wrf_mixrats



	subroutine mosaic_init_wrf_mixrats_opt1(   &
		iflagaa, config_flags,   &
		chem,   &
		ids,ide, jds,jde, kds,kde,   &
		ims,ime, jms,jme, kms,kme,   &
		its,ite, jts,jte, kts,kte    )









	use module_configure, only:  grid_config_rec_type
	use module_state_description, only:  num_chem, param_first_scalar
	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_peg_util, only:  peg_message, peg_error_fatal

	implicit none


	type(grid_config_rec_type), intent(in) :: config_flags

	integer, intent(in) ::   &
		iflagaa,   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem


	integer i, j, k, l, ll, l1, l3, l1dum, m, mdum, nsm
	integer it, jt, kt

	real dum, dumdp, dumrsfc, dumvol,   &
      	  xlo, xhi,   &
	  dumvol1p,  &
      	  pdummb, zdumkm, zscalekm, zfactor

	real vtot_nsm_ofmode(maxd_asize)
	real dumarr(maxd_acomp+5)

	real erfc


	real fracvol, tlo, thi

	integer nmaxd_nsm
	parameter (nmaxd_nsm = 4)

	integer iphase, itype, ntot_nsm
	integer iiprof_nsm(nmaxd_nsm)
	integer lldum_so4, lldum_nh4, lldum_oc, lldum_bc,   &
		lldum_oin, lldum_na, lldum_cl, lldum_hysw

	real sx_nsm(nmaxd_nsm), sxr2_nsm(nmaxd_nsm),   &
      	  x0_nsm(nmaxd_nsm), x3_nsm(nmaxd_nsm),   &
      	  rtot_nsm(maxd_acomp,nmaxd_nsm),   &
      	  vtot_nsm(nmaxd_nsm), xntot_nsm(nmaxd_nsm)

	real dgnum_nsm(nmaxd_nsm), sigmag_nsm(nmaxd_nsm)
	real aaprof_nsm(maxd_acomp+1,nmaxd_nsm)
	real bbprof_nsm(nmaxd_nsm)

	character*80 msg
	character*10 dumname



	if (mosaic_init_wrf_mixrats_flagaa .le. 0) return



	itype = 1
	iphase = ai_phase
	m = 1


	iiprof_nsm(:) = 1
	aaprof_nsm(:,:) = 0.0
	bbprof_nsm(:) = 0.0

	ntot_nsm = 4
	ntot_nsm = min( ntot_nsm, nsize_aer(itype) )

	lldum_so4 = 0
	lldum_nh4 = 0
	lldum_oc  = 0
	lldum_bc  = 0
	lldum_oin = 0
	lldum_na  = 0
	lldum_cl  = 0
	lldum_hysw = 0

	do ll = 1, ncomp_plustracer_aer(itype)
	    if (massptr_aer(ll,m,itype,iphase) .eq.   &
		lptr_so4_aer(m,itype,iphase)) lldum_so4 = ll
	    if (massptr_aer(ll,m,itype,iphase) .eq.   &
		lptr_nh4_aer(m,itype,iphase)) lldum_nh4 = ll
	    if (massptr_aer(ll,m,itype,iphase) .eq.   &
		lptr_oc_aer(m,itype,iphase))  lldum_oc  = ll
	    if (massptr_aer(ll,m,itype,iphase) .eq.   &
		lptr_bc_aer(m,itype,iphase))  lldum_bc  = ll
	    if (massptr_aer(ll,m,itype,iphase) .eq.   &
		lptr_oin_aer(m,itype,iphase)) lldum_oin = ll
	    if (massptr_aer(ll,m,itype,iphase) .eq.   &
		lptr_na_aer(m,itype,iphase))  lldum_na  = ll
	    if (massptr_aer(ll,m,itype,iphase) .eq.   &
		lptr_cl_aer(m,itype,iphase))  lldum_cl  = ll
	end do
	if (hyswptr_aer(m,itype) .gt. 0)   &
		lldum_hysw = ncomp_plustracer_aer(itype) + 1

	msg = ' '
	if (lldum_so4 .le. 0)   &
		msg = '*** subr mosaic_init_wrf_mixrats - lldum_so4 = 0'
	if (lldum_nh4 .le. 0)   &
		msg = '*** subr mosaic_init_wrf_mixrats - lldum_nh4 = 0'
	if (lldum_oc .le. 0)   &
		msg = '*** subr mosaic_init_wrf_mixrats - lldum_oc = 0'
	if (lldum_bc .le. 0)   &
		msg = '*** subr mosaic_init_wrf_mixrats - lldum_bc = 0'
	if (lldum_oin .le. 0)   &
		msg = '*** subr mosaic_init_wrf_mixrats - lldum_oin = 0'
	if (lldum_na .le. 0)   &
		msg = '*** subr mosaic_init_wrf_mixrats - lldum_na = 0'
	if (lldum_cl .le. 0)   &
		msg = '*** subr mosaic_init_wrf_mixrats - lldum_cl = 0'
	if (lldum_hysw .le. 0)   &
		msg = '*** subr mosaic_init_wrf_mixrats - lldum_hysw = 0'
	if (msg .ne. ' ') call peg_error_fatal( lunerr, msg )


	do nsm = 1, ntot_nsm

	if (nsm .eq. 1) then

	    dgnum_nsm( nsm) = 0.15e-4
	    sigmag_nsm(nsm) = 2.0
	    aaprof_nsm(lldum_so4,nsm) = 0.50
	    aaprof_nsm(lldum_nh4,nsm) = aaprof_nsm(lldum_so4,nsm)   &
					* (mw_nh4_aer/mw_so4_aer)
	    aaprof_nsm(lldum_oc,nsm) = 0.25
	    aaprof_nsm(lldum_bc,nsm) = 0.05
	    aaprof_nsm(lldum_hysw,nsm) = aaprof_nsm(lldum_so4,nsm) * 0.2

	else if (nsm .eq. 2) then

	    dgnum_nsm( nsm) = 0.03e-4
	    sigmag_nsm(nsm) = 2.0
	    aaprof_nsm(lldum_so4,nsm) = 0.50 * 0.020
	    aaprof_nsm(lldum_nh4,nsm) = aaprof_nsm(lldum_so4,nsm)   &
					* (mw_nh4_aer/mw_so4_aer)
	    aaprof_nsm(lldum_oc,nsm) = 0.25 * 0.020
	    aaprof_nsm(lldum_bc,nsm) = 0.05 * 0.020
	    aaprof_nsm(lldum_hysw,nsm) = aaprof_nsm(lldum_so4,nsm) * 0.2

	else if (nsm .eq. 3) then

	    dgnum_nsm( nsm) = 1.0e-4
	    sigmag_nsm(nsm) = 2.0
	    aaprof_nsm(lldum_oin,nsm) = 0.5
	    aaprof_nsm(lldum_hysw,nsm) = aaprof_nsm( 9,nsm) * 1.0e-3

	else if (nsm .eq. 4) then

	    dgnum_nsm( nsm) = 2.0e-4
	    sigmag_nsm(nsm) = 2.0
	    aaprof_nsm(lldum_cl,nsm) = 0.1
	    aaprof_nsm(lldum_na,nsm) = aaprof_nsm(lldum_cl,nsm)   &
					* (mw_na_aer/mw_cl_aer)
	    aaprof_nsm(lldum_hysw,nsm) = aaprof_nsm(lldum_cl,nsm) * 0.2

	end if

	end do


	if (iflagaa .gt. 0) then
	    do nsm = 1, ntot_nsm
		if (nsm .ne. iflagaa) aaprof_nsm(:,nsm) = 0.0
	    end do
	end if








	do nsm = 1, ntot_nsm
	    sx_nsm(nsm) = alog( sigmag_nsm(nsm) )
	    sxr2_nsm(nsm) = sx_nsm(nsm) * sqrt(2.0)
	    x0_nsm(nsm) = alog( dgnum_nsm(nsm) )
	    x3_nsm(nsm) = x0_nsm(nsm) + 3.0*sx_nsm(nsm)*sx_nsm(nsm)
	end do


	rclm(:,:) = 0.






	do 12900 k = 1, 1



	zdumkm = 0.0














	do nsm = 1, ntot_nsm

	if (iiprof_nsm(nsm) .eq. 2) then
	    zscalekm = bbprof_nsm(nsm)
	    zfactor = exp( -zdumkm/zscalekm )
	else if (iiprof_nsm(nsm) .eq. 3) then
	    zscalekm = bbprof_nsm(nsm)
	    zfactor = max( 0., (1. - zdumkm/zscalekm) )
	else
	    zfactor = 1.0
	end if

	    do ll = 1, ncomp_plustracer_aer(itype) + 1
		rtot_nsm(ll,nsm) = aaprof_nsm(ll,nsm) * zfactor
	    end do

	end do



	do nsm = 1, ntot_nsm
	    dumvol = 0.
	    do ll = 1, ncomp_aer(itype)
		dum = 1.0e-6*rtot_nsm(ll,nsm)/dens_aer(ll,itype)
		dumvol = dumvol + max( 0., dum )
	    end do
	    vtot_nsm(nsm) = dumvol
	end do


	do 12700 m = 1, nsize_aer(itype)

	vtot_nsm_ofmode(m) = 0.0

	do 12500 nsm = 1, ntot_nsm



	xlo = alog( dlo_sect(m,itype) )
	xhi = alog( dhi_sect(m,itype) )

	tlo = (xlo - x3_nsm(nsm))/sxr2_nsm(nsm)
	thi = (xhi - x3_nsm(nsm))/sxr2_nsm(nsm)
	if (tlo .ge. 0.) then

	    fracvol = 0.5*( erfc_num_recipes(tlo) - erfc_num_recipes(thi) )
	else

	    fracvol = 0.5*( erfc_num_recipes(-thi) - erfc_num_recipes(-tlo) )
	end if
	fracvol = max( fracvol, 0.0 )

	vtot_nsm_ofmode(m) = vtot_nsm_ofmode(m)  + vtot_nsm(nsm)*fracvol



	do ll = 1, ncomp_plustracer_aer(itype)
	    rclm( k, massptr_aer(ll,m,itype,iphase) ) =   &
      	    rclm( k, massptr_aer(ll,m,itype,iphase) ) +   &
      		fracvol*rtot_nsm(ll,nsm)
	end do

	if ((iphase .eq. ai_phase) .and.   &
	    (lldum_hysw .gt. 0) .and.   &
	    (hyswptr_aer(m,itype) .gt. 0)) then

	    rclm( k, hyswptr_aer(m,itype) ) =   &
      	    rclm( k, hyswptr_aer(m,itype) ) +   &
      		fracvol*rtot_nsm(lldum_hysw,nsm)
	end if

12500	continue


	dum = sqrt( dlo_sect(m,itype)*dhi_sect(m,itype) )
	dumvol1p = (pi/6.0)*(dum**3)
	rclm( k, numptr_aer(m,itype,iphase) ) = vtot_nsm_ofmode(m)/dumvol1p


	if ((iphase .eq. ai_phase) .and.   &
	    (lldum_hysw .gt. 0) .and.   &
	    (hyswptr_aer(m,itype) .gt. 0) .and.   &
	    (waterptr_aer(m,itype) .gt. 0)) then

      	    rclm( k, waterptr_aer(m,itype) ) =    &
	    rclm( k, hyswptr_aer(m,itype) )
	end if

12700	continue

12900	continue









	dumarr(:) = 0.0
	msg = ' '
	call peg_message( lunout, msg )
	msg = '*** subr mosaic_init_wrf_mixrats_opt1 results'
	call peg_message( lunout, msg )
	msg = '    mass in ug/m3     number in #/m3     volume in cm3/m3'
	call peg_message( lunout, msg )
	msg = ' '
	call peg_message( lunout, msg )
	msg = ' mode  l  l1  species      conc'
	call peg_message( lunout, msg )

        do 14390 mdum = 1, nsize_aer(itype)+1
	m = min( mdum, nsize_aer(itype) )
	msg = ' '
	call peg_message( lunout, msg )
        do 14350 l = 1, ncomp_plustracer_aer(itype)+4

            if (l .le. ncomp_plustracer_aer(itype)) then
                l1 = massptr_aer(l,m,itype,iphase)
		dumname = name_aer(l,itype)
		dum = rclm(1,l1)
            else if (l .eq. ncomp_plustracer_aer(itype)+1) then
                l1 = hyswptr_aer(m,itype)
		dumname = 'hystwatr'
		dum = rclm(1,l1)
            else if (l .eq. ncomp_plustracer_aer(itype)+2) then
                l1 = waterptr_aer(m,itype)
		dumname = 'water'
		dum = rclm(1,l1)
            else if (l .eq. ncomp_plustracer_aer(itype)+3) then
                l1 = numptr_aer(m,itype,iphase)
		dumname = 'number'
		dum = rclm(1,l1)
            else if (l .eq. ncomp_plustracer_aer(itype)+4) then
                l1 = 0
		dumname = 'volume'
		dum = vtot_nsm_ofmode(m)
	    else
		dumname = '=BADBAD='
		l1 = -1
		dum = -1.0
            end if

	    l1dum = l1
	    if (aboxtest_testmode .gt. 0) l1dum = max( l1-1, 0 )

	    if (mdum .le. nsize_aer(itype)) then
		dumarr(l) = dumarr(l) + dum
		write(msg,9620) m, l, l1dum, dumname, dum
	    else
		write(msg,9625) l, dumname, dumarr(l)
	    end if
	    call peg_message( lunout, msg )

14350   continue
14390   continue

9620    format( 3i4, 2x, a, 3(1pe12.3) )
9625    format( ' sum', i4, '   -', 2x, a, 3(1pe12.3) )





        do 16390 m = 1, nsize_aer(itype)
        do 16350 l = 1, 15

	    if (l .eq. 1) then
		l1 = lptr_so4_aer(m,itype,iphase)
	    else if (l .eq. 2) then
		l1 = lptr_no3_aer(m,itype,iphase)
	    else if (l .eq. 3) then
		l1 = lptr_cl_aer(m,itype,iphase)
	    else if (l .eq. 4) then
		l1 = lptr_msa_aer(m,itype,iphase)
	    else if (l .eq. 5) then
		l1 = lptr_co3_aer(m,itype,iphase)
	    else if (l .eq. 6) then
		l1 = lptr_nh4_aer(m,itype,iphase)
	    else if (l .eq. 7) then
		l1 = lptr_na_aer(m,itype,iphase)
	    else if (l .eq. 8) then
		l1 = lptr_ca_aer(m,itype,iphase)
	    else if (l .eq. 9) then
		l1 = lptr_oin_aer(m,itype,iphase)
	    else if (l .eq. 10) then
		l1 = lptr_oc_aer(m,itype,iphase)
	    else if (l .eq. 11) then
		l1 = lptr_bc_aer(m,itype,iphase)
	    else if (l .eq. 12) then
		l1 = hyswptr_aer(m,itype)
	    else if (l .eq. 13) then
		l1 = waterptr_aer(m,itype)
	    else if (l .eq. 14) then
		l1 = numptr_aer(m,itype,iphase)
	    else
		goto 16350
	    end if
	    l3 = l1

	    if ((l1 .gt. 0) .and. (l1 .le. lmaxd) .and.   &
	        (l3 .ge. param_first_scalar)) then
		do it = its, ite

		do kt = kts, kte-1
		do jt = jts, jte
		    chem(it,kt,jt,l3) = rclm(1,l1)
		end do
		end do
		end do
	    end if

16350   continue
16390   continue



	return
	end subroutine mosaic_init_wrf_mixrats_opt1




	subroutine mosaic_init_wrf_mixrats_opt2(   &
		iflagaa, config_flags,             &
		chem, z_at_w, g,                   &
		ids,ide, jds,jde, kds,kde,         &
		ims,ime, jms,jme, kms,kme,         &
		its,ite, jts,jte, kts,kte          )










	use module_configure, only:  grid_config_rec_type
	use module_state_description, only:  num_chem, param_first_scalar
	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_peg_util, only:  peg_message, peg_error_fatal

	implicit none


	type(grid_config_rec_type), intent(in) :: config_flags

	integer, intent(in) ::   &
		iflagaa,   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(inout),     &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem

        real, intent(in),        &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
                z_at_w
        real :: g


	integer l, l1, l3, m, mdum
	integer iphase, itype
	integer it, jt, kt


	real dum, dumvol1p, mult
	real qcoar, qfine, qval

	real :: vtot_ofmode(maxd_asize)
	real :: dumarr(maxd_acomp+5)
	real :: fr_coar(8), fr_fine(8)






	real, save :: fr8b_coar(8) =   &
        (/ 0.0,  0.0,  0.0,  0.019,  0.0212,  0.1101,  0.2751, 0.7018 /) 




	real, save :: fr8b_fine(8) =   &
        (/ 0.0245, 0.1472, 0.3501, 0.3317, 0.1251, 0.0186, 0.0011, 0.0/) 





	real :: qfine_so4, qfine_no3, qfine_cl, qfine_msa,   &
		qfine_co3, qfine_nh4, qfine_na, qfine_ca, qfine_oin,   &
		qfine_oc, qfine_bc, qfine_hysw, qfine_watr, qfine_vol
	real :: qcoar_so4, qcoar_no3, qcoar_cl, qcoar_msa,   &
		qcoar_co3, qcoar_nh4, qcoar_na, qcoar_ca, qcoar_oin,   &
		qcoar_oc, qcoar_bc, qcoar_hysw, qcoar_watr, qcoar_vol


        real, dimension( ims:ime, kms:kme, jms:jme ) :: z

	character*80 msg
	character*10 dumname



	itype = 1
	iphase = ai_phase
	m = 1



    do jt = jts, min(jte,jde-1)
       do kt = kts, kte-1
          do it = its, min(ite,ide-1)
             z(it,kt,jt) = (z_at_w(it,kt,jt)+z_at_w(it,kt+1,jt))*0.5
          end do
       end do
    end do



	if (nsize_aer(itype) == 8) then
	    fr_coar(:) = fr8b_coar(:)
	    fr_fine(:) = fr8b_fine(:)
	else if (nsize_aer(itype) == 4) then
	    do m = 1, nsize_aer(itype)
		fr_coar(m) = fr8b_coar(2*m) + fr8b_coar(2*m-1)
		fr_fine(m) = fr8b_fine(2*m) + fr8b_fine(2*m-1)
	    end do
	else
	    write(msg,'(a,i5)')   &
		'subr mosaic_init_wrf_mixrats_opt2' //   &
		' - nsize_aer(itype) must be 4 or 8 but = ', nsize_aer(itype)
	    call peg_error_fatal( lunout, msg )
	end if






	rclm(:,:) = 0.0
	vtot_ofmode(:) = 0.0







	qfine_so4 = 2.14
	qcoar_so4 = 0.242
	qfine_no3 = 0.11
	qcoar_no3 = 0.03

	qfine_cl  = 0.14     
	qcoar_cl  = 0.139
	qfine_msa = 0.0
	qcoar_msa = 0.0
	qfine_co3 = 0.0
	qcoar_co3 = 0.0
	qfine_nh4 = 0.83
	qcoar_nh4 = 0.10

	qfine_na  = 0.1      
	qcoar_na  = 0.09
	qfine_ca  = 0.0
	qcoar_ca  = 0.0

	qfine_oin = 3.48     
	qcoar_oin = 0.35
	qfine_oc  = 1.00
	qcoar_oc  = 1.50
	qfine_bc  = 0.2
	qcoar_bc  = 0.075
	qfine_hysw = 0.0
	qcoar_hysw = 0.0
	qfine_watr = 0.0
	qcoar_watr = 0.0
































	qfine_vol = 1.0e-6 * (   &
	    (qfine_so4/dens_so4_aer) + (qfine_no3/dens_no3_aer) +   &
	    (qfine_cl /dens_cl_aer ) + (qfine_msa/dens_msa_aer) +   &
	    (qfine_co3/dens_co3_aer) + (qfine_nh4/dens_nh4_aer) +   &
	    (qfine_na /dens_na_aer ) + (qfine_ca /dens_ca_aer ) +   &
	    (qfine_oin/dens_oin_aer) + (qfine_oc /dens_oc_aer ) +   &
	    (qfine_bc /dens_bc_aer ) )
	qcoar_vol =  1.0e-6 * (  &
	    (qcoar_so4/dens_so4_aer) + (qcoar_no3/dens_no3_aer) +   &
	    (qcoar_cl /dens_cl_aer ) + (qcoar_msa/dens_msa_aer) +   &
	    (qcoar_co3/dens_co3_aer) + (qcoar_nh4/dens_nh4_aer) +   &
	    (qcoar_na /dens_na_aer ) + (qcoar_ca /dens_ca_aer ) +   &
	    (qcoar_oin/dens_oin_aer) + (qcoar_oc /dens_oc_aer ) +   &
	    (qcoar_bc /dens_bc_aer ) )

        do 2900 m = 1, nsize_aer(itype)
        do 2800 l = 1, 15

	    if (l .eq. 1) then
		l1 = lptr_so4_aer(m,itype,iphase)
		qfine = qfine_so4
		qcoar = qcoar_so4
	    else if (l .eq. 2) then
		l1 = lptr_no3_aer(m,itype,iphase)
		qfine = qfine_no3
		qcoar = qcoar_no3
	    else if (l .eq. 3) then
		l1 = lptr_cl_aer(m,itype,iphase)
		qfine = qfine_cl
		qcoar = qcoar_cl
	    else if (l .eq. 4) then
		l1 = lptr_msa_aer(m,itype,iphase)
		qfine = qfine_msa
		qcoar = qcoar_msa
	    else if (l .eq. 5) then
		l1 = lptr_co3_aer(m,itype,iphase)
		qfine = qfine_co3
		qcoar = qcoar_co3
	    else if (l .eq. 6) then
		l1 = lptr_nh4_aer(m,itype,iphase)
		qfine = qfine_nh4
		qcoar = qcoar_nh4
	    else if (l .eq. 7) then
		l1 = lptr_na_aer(m,itype,iphase)
		qfine = qfine_na
		qcoar = qcoar_na
	    else if (l .eq. 8) then
		l1 = lptr_ca_aer(m,itype,iphase)
		qfine = qfine_ca
		qcoar = qcoar_ca
	    else if (l .eq. 9) then
		l1 = lptr_oin_aer(m,itype,iphase)
		qfine = qfine_oin
		qcoar = qcoar_oin
	    else if (l .eq. 10) then
		l1 = lptr_oc_aer(m,itype,iphase)
		qfine = qfine_oc
		qcoar = qcoar_oc
	    else if (l .eq. 11) then
		l1 = lptr_bc_aer(m,itype,iphase)
		qfine = qfine_bc
		qcoar = qcoar_bc
	    else if (l .eq. 12) then
		l1 = hyswptr_aer(m,itype)
		qfine = qfine_hysw
		qcoar = qcoar_hysw
	    else if (l .eq. 13) then
		l1 = waterptr_aer(m,itype)
		qfine = qfine_watr
		qcoar = qcoar_watr
	    else if (l .eq. 14) then
		l1 = numptr_aer(m,itype,iphase)
		dumvol1p = sqrt(volumlo_sect(m,itype)*volumhi_sect(m,itype))
		qfine = qfine_vol/dumvol1p
		qcoar = qcoar_vol/dumvol1p
		vtot_ofmode(m) =   &
			qfine_vol*fr_fine(m) + qcoar_vol*fr_coar(m)
	    else
		goto 2800
	    end if
	    l3 = l1

	    if ((l1 .gt. 0) .and. (l1 .le. lmaxd) .and.   &
	        (l3 .ge. param_first_scalar)) then
		qval = qfine*fr_fine(m) + qcoar*fr_coar(m)
		rclm(1,l1) = qval

		do jt = jts, min(jte,jde-1)
		do kt = kts, kte-1
		do it = its, min(ite,ide-1)







































		   if( z(it,kt,jt) <= 2000. ) then
		      mult = 1.0
		   elseif( z(it,kt,jt) > 2000. &
                        .and. z(it,kt,jt) <= 3000. ) then
		      mult = 1.0 - 0.00075*(z(it,kt,jt)-2000.)
		   elseif( z(it,kt,jt) > 3000. &
                        .and. z(it,kt,jt) <= 5000. ) then
		      mult = 0.25 - 4.166666667e-5*(z(it,kt,jt)-3000.)
		   else
		      mult = 0.125
           end if

		    chem(it,kt,jt,l3) = mult*rclm(1,l1)


		end do
		end do
		end do
	    end if


2800   continue
2900   continue




	dumarr(:) = 0.0
	msg = ' '
	call peg_message( lunout, msg )
	msg = '*** subr mosaic_init_wrf_mixrats_opt2 results'
	call peg_message( lunout, msg )
	msg = '    mass in ug/m3     number in #/m3     volume in cm3/m3'
	call peg_message( lunout, msg )
	msg = ' '
	call peg_message( lunout, msg )
	msg = ' mode  l  l1  species      conc'
	call peg_message( lunout, msg )

        do 3190 mdum = 1, nsize_aer(itype)+1
	m = min( mdum, nsize_aer(itype) )
	msg = ' '
	call peg_message( lunout, msg )
        do 3150 l = 1, ncomp_plustracer_aer(itype)+4

            if (l .le. ncomp_plustracer_aer(itype)) then
                l1 = massptr_aer(l,m,itype,iphase)
		dumname = name_aer(l,itype)
		dum = rclm(1,l1)
            else if (l .eq. ncomp_plustracer_aer(itype)+1) then
                l1 = hyswptr_aer(m,itype)
		dumname = 'hystwatr'
		dum = rclm(1,l1)
            else if (l .eq. ncomp_plustracer_aer(itype)+2) then
                l1 = waterptr_aer(m,itype)
		dumname = 'water'
		dum = rclm(1,l1)
            else if (l .eq. ncomp_plustracer_aer(itype)+3) then
                l1 = numptr_aer(m,itype,iphase)
		dumname = 'number'
		dum = rclm(1,l1)
            else if (l .eq. ncomp_plustracer_aer(itype)+4) then
                l1 = 0
		dumname = 'volume'
		dum = vtot_ofmode(m)
	    else
		dumname = '=BADBAD='
		l1 = -1
		dum = -1.0
            end if

	    if (mdum .le. nsize_aer(itype)) then
		dumarr(l) = dumarr(l) + dum
		write(msg,9620) m, l, l1, dumname, dum
	    else
		write(msg,9625) l, dumname, dumarr(l)
	    end if
	    call peg_message( lunout, msg )

3150   continue
3190   continue

9620    format( 3i4, 2x, a, 3(1pe12.3) )
9625    format( ' sum', i4, '   -', 2x, a, 3(1pe12.3) )



	return
	end subroutine mosaic_init_wrf_mixrats_opt2



	real function erfc_num_recipes( x )



	implicit none
	real x
	double precision erfc_dbl, dum, t, zz

	zz = abs(x)
	t = 1.0/(1.0 + 0.5*zz)







	dum =  ( -zz*zz - 1.26551223 + t*(1.00002368 + t*(0.37409196 +   &
      	  t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +   &
      	                                   t*(-1.13520398 +   &
      	  t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))

	erfc_dbl = t * exp(dum)
	if (x .lt. 0.0) erfc_dbl = 2.0d0 - erfc_dbl

	erfc_num_recipes = erfc_dbl

	return
	end function erfc_num_recipes     



	end module module_mosaic_initmixrats





 	subroutine bdy_chem_value_mosaic ( chem_bv, alt, zz, nch, config_flags )









	use module_configure, only:  grid_config_rec_type
	use module_state_description, only:  param_first_scalar,   &
			aer_bc_pnnl
	use module_data_mosaic_asect
	use module_data_mosaic_other
	implicit none


	REAL,    intent(OUT)  :: chem_bv    
	REAL,    intent(IN)   :: alt        
	REAL,    intent(IN)   :: zz         
	INTEGER, intent(IN)   :: nch        
	TYPE (grid_config_rec_type), intent(in) :: config_flags


	integer :: iphase, itype, m
    logical :: foundit

	real, parameter :: chem_bv_def = 1.0e-20

	real :: dumvol1p, mult
        integer :: iflg
	real :: qcoar, qfine, qval

	real :: fr_coar(8), fr_fine(8)





	real, save :: fr8b_coar(8) =   &
        (/ 0.0,  0.0,  0.0,  0.019,  0.0212,  0.1101,  0.2751, 0.7018 /) 



	real, save :: fr8b_fine(8) =   &
        (/ 0.0245, 0.1472, 0.3501, 0.3317, 0.1251, 0.0186, 0.0011, 0.0/) 





	real :: qfine_so4, qfine_no3, qfine_cl, qfine_msa,   &
		qfine_co3, qfine_nh4, qfine_na, qfine_ca, qfine_oin,   &
		qfine_oc, qfine_bc, qfine_hysw, qfine_watr, qfine_vol
	real :: qcoar_so4, qcoar_no3, qcoar_cl, qcoar_msa,   &
		qcoar_co3, qcoar_nh4, qcoar_na, qcoar_ca, qcoar_oin,   &
		qcoar_oc, qcoar_bc, qcoar_hysw, qcoar_watr, qcoar_vol

	character*80 msg





	chem_bv = chem_bv_def
	if (config_flags%aer_bc_opt /= aer_bc_pnnl) return
	if (nch < param_first_scalar) return



	itype = 1
	iphase = ai_phase
	m = 1 







	if (nsize_aer(itype) == 8) then
	    fr_coar(:) = fr8b_coar(:)
	    fr_fine(:) = fr8b_fine(:)
	else if (nsize_aer(itype) == 4) then
	    do m = 1, nsize_aer(itype)
		fr_coar(m) = fr8b_coar(2*m) + fr8b_coar(2*m-1)
		fr_fine(m) = fr8b_fine(2*m) + fr8b_fine(2*m-1)
	    end do
	else
	    write(msg,'(a,i5)')   &
		'subr bdy_chem_value_mosaic' //   &
		' - nsize_aer(itype) must be 4 or 8 but = ', nsize_aer(itype)
	    call wrf_error_fatal3("<stdin>",1086,&
msg )
	end if






































        if( zz <= 2000. ) then
           mult = 1.0
        elseif( zz > 2000. &
             .and. zz <= 3000. ) then
           mult = 1.0 - 0.00075*(zz-2000.)
        elseif( zz > 3000. &
             .and. zz <= 5000. ) then
           mult = 0.25 - 4.166666667e-5*(zz-3000.)
        else
           mult = 0.125
        end if





	qfine_so4 = mult*2.14
	qcoar_so4 = mult*0.242
	qfine_no3 = mult*0.11
	qcoar_no3 = mult*0.03

	qfine_cl  = mult*0.14     
	qcoar_cl  = mult*0.139
	qfine_msa = mult*0.0
	qcoar_msa = mult*0.0
	qfine_co3 = mult*0.0
	qcoar_co3 = mult*0.0
	qfine_nh4 = mult*0.83
	qcoar_nh4 = mult*0.10

	qfine_na  = mult*0.1      
	qcoar_na  = mult*0.09
	qfine_ca  = mult*0.0
	qcoar_ca  = mult*0.0


	qfine_oin = mult*3.48     
	qcoar_oin = mult*0.35
	qfine_oc  = mult*1.00
	qcoar_oc  = mult*1.50
	qfine_bc  = mult*0.2
	qcoar_bc  = mult*0.075
	qfine_hysw = mult*0.0
	qcoar_hysw = mult*0.0
	qfine_watr = mult*0.0
	qcoar_watr = mult*0.0































	qfine_vol = 1.0e-6 * (   &
	    (qfine_so4/dens_so4_aer) + (qfine_no3/dens_no3_aer) +   &
	    (qfine_cl /dens_cl_aer ) + (qfine_msa/dens_msa_aer) +   &
	    (qfine_co3/dens_co3_aer) + (qfine_nh4/dens_nh4_aer) +   &
	    (qfine_na /dens_na_aer ) + (qfine_ca /dens_ca_aer ) +   &
	    (qfine_oin/dens_oin_aer) + (qfine_oc /dens_oc_aer ) +   &
	    (qfine_bc /dens_bc_aer ) )
	qcoar_vol = 1.0e-6 * (   &
	    (qcoar_so4/dens_so4_aer) + (qcoar_no3/dens_no3_aer) +   &
	    (qcoar_cl /dens_cl_aer ) + (qcoar_msa/dens_msa_aer) +   &
	    (qcoar_co3/dens_co3_aer) + (qcoar_nh4/dens_nh4_aer) +   &
	    (qcoar_na /dens_na_aer ) + (qcoar_ca /dens_ca_aer ) +   &
	    (qcoar_oin/dens_oin_aer) + (qcoar_oc /dens_oc_aer ) +   &
	    (qcoar_bc /dens_bc_aer ) )

	qfine = -1.0e30
	qcoar = -1.0e30


        do 2900 m = 1, nsize_aer(itype)

	    if (nch .eq. lptr_so4_aer(m,itype,iphase)) then
		qfine = qfine_so4
		qcoar = qcoar_so4
	    else if (nch .eq. lptr_no3_aer(m,itype,iphase)) then
		qfine = qfine_no3
		qcoar = qcoar_no3
	    else if (nch .eq. lptr_cl_aer(m,itype,iphase)) then
		qfine = qfine_cl
		qcoar = qcoar_cl
	    else if (nch .eq. lptr_msa_aer(m,itype,iphase)) then
		qfine = qfine_msa
		qcoar = qcoar_msa
	    else if (nch .eq. lptr_co3_aer(m,itype,iphase)) then
		qfine = qfine_co3
		qcoar = qcoar_co3
	    else if (nch .eq. lptr_nh4_aer(m,itype,iphase)) then
		qfine = qfine_nh4
		qcoar = qcoar_nh4
	    else if (nch .eq. lptr_na_aer(m,itype,iphase)) then
		qfine = qfine_na
		qcoar = qcoar_na
	    else if (nch .eq. lptr_ca_aer(m,itype,iphase)) then
		qfine = qfine_ca
		qcoar = qcoar_ca
	    else if (nch .eq. lptr_oin_aer(m,itype,iphase)) then
		qfine = qfine_oin
		qcoar = qcoar_oin
	    else if (nch .eq. lptr_oc_aer(m,itype,iphase)) then
		qfine = qfine_oc
		qcoar = qcoar_oc
	    else if (nch .eq. lptr_bc_aer(m,itype,iphase)) then
		qfine = qfine_bc
		qcoar = qcoar_bc
	    else if (nch .eq. hyswptr_aer(m,itype)) then
		qfine = qfine_hysw
		qcoar = qcoar_hysw
	    else if (nch .eq. waterptr_aer(m,itype)) then
		qfine = qfine_watr
		qcoar = qcoar_watr
	    else if (nch .eq. numptr_aer(m,itype,iphase)) then
		dumvol1p = sqrt(volumlo_sect(m,itype)*volumhi_sect(m,itype))
		qfine = qfine_vol/dumvol1p
		qcoar = qcoar_vol/dumvol1p
	    end if

	    if ((qfine >= 0.0) .and. (qcoar >= 0.0)) then
		qval = qfine*fr_fine(m) + qcoar*fr_coar(m)
		chem_bv = qval*alt
		goto 2910
	    end if

2900   continue
2910   continue

	return
 	end subroutine bdy_chem_value_mosaic


