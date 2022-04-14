


















  module modal_aero_rename


  use shr_kind_mod, only : r8 => shr_kind_r8
  use modal_aero_data, only:  maxd_aspectype
  use module_cam_support, only: &
        pcnstxx => gas_pcnst_modal_aero, &
        endrun

  implicit none
  private
  save


  public modal_aero_rename_sub, modal_aero_rename_init


  integer, parameter, public :: maxpair_renamexf = 1
  integer, parameter, public :: maxspec_renamexf = maxd_aspectype

  integer, public :: npair_renamexf = -123456789
  integer, public :: modefrm_renamexf(maxpair_renamexf)
  integer, public :: modetoo_renamexf(maxpair_renamexf)
  integer, public :: nspecfrm_renamexf(maxpair_renamexf)
  integer, public :: lspecfrmc_renamexf(maxspec_renamexf,maxpair_renamexf)
  integer, public :: lspecfrma_renamexf(maxspec_renamexf,maxpair_renamexf)
  integer, public :: lspectooc_renamexf(maxspec_renamexf,maxpair_renamexf)
  integer, public :: lspectooa_renamexf(maxspec_renamexf,maxpair_renamexf)

















  contains
                                                                                                                                            





	subroutine modal_aero_rename_sub(                       &
                        fromwhere,         lchnk,               &
                        ncol,              nstep,               &
                        loffset,           deltat,              &
                        latndx,            lonndx,              &
                        pdel,                                   &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx            )


   use shr_kind_mod, only : r4 => shr_kind_r4
   use modal_aero_data
   use physconst, only: gravit, mwdry
   use module_cam_support, only: &
        pcnst => pcnst_runtime, &
        pcols, pver
   use constituents, only: cnst_name
   use physconst, only: pi

   implicit none



   character(len=*), intent(in) :: fromwhere    
                                                
   integer,  intent(in)    :: lchnk                
   integer,  intent(in)    :: ncol                 
   integer,  intent(in)    :: nstep                
   integer,  intent(in)    :: loffset              
   real(r8), intent(in)    :: deltat               
   integer,  intent(in)    :: latndx(pcols), lonndx(pcols)

   real(r8), intent(in)    :: pdel(pcols,pver)     
   real(r8), intent(in)    :: q(ncol,pver,pcnstxx) 
                                                   
                                                   
   real(r8), intent(in)    :: qqcw(ncol,pver,pcnstxx) 

   real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)  
                              
                              
                              
                              
   real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
   real(r8), intent(in)    :: dqdt_other(ncol,pver,pcnstxx)  
                              
                              
                              
                              
                              
   real(r8), intent(in)    :: dqqcwdt_other(ncol,pver,pcnstxx)  
   logical,  intent(inout) :: dotendrn(pcnstxx) 
                              
   logical,  intent(inout) :: dotendqqcwrn(pcnstxx)

   logical,  intent(in)    :: is_dorename_atik          
   logical,  intent(in)    :: dorename_atik(ncol,pver) 
                                                        
   integer,  intent(in)    :: jsrflx_rename        
   integer,  intent(in)    :: nsrflx               

   real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
                              
   real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)

















   integer, parameter :: ldiag1=-1
   integer :: i, icol_diag, ipair, iq, j, k, l, l1, l2, la, lc, lunout
   integer :: lsfrma, lsfrmc, lstooa, lstooc
   integer :: mfrm, mtoo, n, n1, n2, ntot_msa_a
   integer :: idomode(ntot_amode)
   integer, save :: lun = -1  
                              


   real (r8) :: deldryvol_a(ncol,pver,ntot_amode)
   real (r8) :: deldryvol_c(ncol,pver,ntot_amode)
   real (r8) :: deltatinv
   real (r8) :: dp_belowcut(maxpair_renamexf)
   real (r8) :: dp_cut(maxpair_renamexf)
   real (r8) :: dgn_aftr, dgn_xfer
   real (r8) :: dgn_t_new, dgn_t_old
   real (r8) :: dryvol_t_del, dryvol_t_new
   real (r8) :: dryvol_t_old, dryvol_t_oldbnd
   real (r8) :: dryvol_a(ncol,pver,ntot_amode)
   real (r8) :: dryvol_c(ncol,pver,ntot_amode)
   real (r8) :: dryvol_smallest(ntot_amode)
   real (r8) :: dum
   real (r8) :: dum3alnsg2(maxpair_renamexf)
   real (r8) :: dum_m2v, dum_m2vdt
   real (r8) :: factoraa(ntot_amode)
   real (r8) :: factoryy(ntot_amode)
   real (r8) :: frelax
   real (r8) :: lndp_cut(maxpair_renamexf)
   real (r8) :: lndgn_new, lndgn_old
   real (r8) :: lndgv_new, lndgv_old
   real (r8) :: num_t_old, num_t_oldbnd
   real (r8) :: onethird
   real (r8) :: pdel_fac
   real (r8) :: tailfr_volnew, tailfr_volold
   real (r8) :: tailfr_numnew, tailfr_numold
   real (r8) :: v2nhirlx(ntot_amode), v2nlorlx(ntot_amode)
   real (r8) :: xfercoef, xfertend
   real (r8) :: xferfrac_vol, xferfrac_num, xferfrac_max

   real (r4) :: yn_tail, yv_tail
   real (r4) :: erfc



	lunout = 6


 	lun = 6














	if (npair_renamexf .le. 0) return








	deltatinv = 1.0/(deltat*(1.0d0 + 1.0d-15))
	onethird = 1.0/3.0
	frelax = 27.0
	xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   

	do n = 1, ntot_amode
	    idomode(n) = 0
	end do

	do ipair = 1, npair_renamexf
	    if (ipair .gt. 1) goto 8100
	    idomode(modefrm_renamexf(ipair)) = 1

	    mfrm = modefrm_renamexf(ipair)
	    mtoo = modetoo_renamexf(ipair)
	    factoraa(mfrm) = (pi/6.)*exp(4.5*(alnsg_amode(mfrm)**2))
	    factoraa(mtoo) = (pi/6.)*exp(4.5*(alnsg_amode(mtoo)**2))
	    factoryy(mfrm) = sqrt( 0.5 )/alnsg_amode(mfrm)



	    dryvol_smallest(mfrm) = 1.0e-25
	    v2nlorlx(mfrm) = voltonumblo_amode(mfrm)*frelax
	    v2nhirlx(mfrm) = voltonumbhi_amode(mfrm)/frelax

	    dum3alnsg2(ipair) = 3.0 * (alnsg_amode(mfrm)**2)
	    dp_cut(ipair) = sqrt(   &
		dgnum_amode(mfrm)*exp(1.5*(alnsg_amode(mfrm)**2)) *   &
		dgnum_amode(mtoo)*exp(1.5*(alnsg_amode(mtoo)**2)) )
	    lndp_cut(ipair) = log( dp_cut(ipair) )
	    dp_belowcut(ipair) = 0.99*dp_cut(ipair)
	end do

	do n = 1, ntot_amode
	    if (idomode(n) .gt. 0) then
		dryvol_a(1:ncol,:,n) = 0.0_r8
		dryvol_c(1:ncol,:,n) = 0.0_r8
		deldryvol_a(1:ncol,:,n) = 0.0_r8
		deldryvol_c(1:ncol,:,n) = 0.0_r8
		do l1 = 1, nspec_amode(n)
		    l2 = lspectype_amode(l1,n)


		    dum_m2v = specmw_amode(l2) / specdens_amode(l2)
		    dum_m2vdt = dum_m2v*deltat
		    la = lmassptr_amode(l1,n)-loffset
		    if (la > 0) then
		    dryvol_a(1:ncol,:,n) = dryvol_a(1:ncol,:,n)    &
			+ dum_m2v*max( 0.0_r8,   &
                          q(1:ncol,:,la)-deltat*dqdt_other(1:ncol,:,la) )
		    deldryvol_a(1:ncol,:,n) = deldryvol_a(1:ncol,:,n)    &
			+ (dqdt_other(1:ncol,:,la) + dqdt(1:ncol,:,la))*dum_m2vdt
		    end if

		    lc = lmassptrcw_amode(l1,n)-loffset
		    if (lc > 0) then
		    dryvol_c(1:ncol,:,n) = dryvol_c(1:ncol,:,n)    &
			+ dum_m2v*max( 0.0_r8,   &
                          qqcw(1:ncol,:,lc)-deltat*dqqcwdt_other(1:ncol,:,lc) )
		    deldryvol_c(1:ncol,:,n) = deldryvol_c(1:ncol,:,n)    &
			+ (dqqcwdt_other(1:ncol,:,lc) +   &
			         dqqcwdt(1:ncol,:,lc))*dum_m2vdt
		    end if
		end do
	    end if
	end do






mainloop1_k:  do k = 1, pver
mainloop1_i:  do i = 1, ncol


	if (is_dorename_atik) then
	    if (.not. dorename_atik(i,k)) cycle mainloop1_i
	end if
	pdel_fac = pdel(i,k)/gravit




mainloop1_ipair:  do ipair = 1, npair_renamexf

	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)





	dryvol_t_old = dryvol_a(i,k,mfrm) + dryvol_c(i,k,mfrm)
	dryvol_t_del = deldryvol_a(i,k,mfrm) + deldryvol_c(i,k,mfrm)
	dryvol_t_new = dryvol_t_old + dryvol_t_del
	dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )


	if (dryvol_t_new .le. dryvol_smallest(mfrm)) cycle mainloop1_ipair
	if (dryvol_t_del .le. 1.0e-6*dryvol_t_oldbnd) cycle mainloop1_ipair


	num_t_old = q(i,k,numptr_amode(mfrm)-loffset)
	num_t_old = num_t_old + qqcw(i,k,numptrcw_amode(mfrm)-loffset)
	num_t_old = max( 0.0_r8, num_t_old )
	dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
	num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
	num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )


	dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
	if (dgn_t_new .le. dgnum_amode(mfrm)) cycle mainloop1_ipair


	lndgn_new = log( dgn_t_new )
	lndgv_new = lndgn_new + dum3alnsg2(ipair)
	yn_tail = (lndp_cut(ipair) - lndgn_new)*factoryy(mfrm)
	yv_tail = (lndp_cut(ipair) - lndgv_new)*factoryy(mfrm)
	tailfr_numnew = 0.5*erfc( yn_tail )
	tailfr_volnew = 0.5*erfc( yv_tail )


	dgn_t_old =   &
		(dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird


	if (dgn_t_new .ge. dp_cut(ipair)) then
	    dgn_t_old = min( dgn_t_old, dp_belowcut(ipair) )
	end if
	lndgn_old = log( dgn_t_old )
	lndgv_old = lndgn_old + dum3alnsg2(ipair)
	yn_tail = (lndp_cut(ipair) - lndgn_old)*factoryy(mfrm)
	yv_tail = (lndp_cut(ipair) - lndgv_old)*factoryy(mfrm)
	tailfr_numold = 0.5*erfc( yn_tail )
	tailfr_volold = 0.5*erfc( yv_tail )



	dum = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_old
	if (dum .le. 0.0_r8) cycle mainloop1_ipair

	xferfrac_vol = min( dum, dryvol_t_new )/dryvol_t_new
	xferfrac_vol = min( xferfrac_vol, xferfrac_max ) 
	xferfrac_num = tailfr_numnew - tailfr_numold
	xferfrac_num = max( 0.0_r8, min( xferfrac_num, xferfrac_vol ) )


 	if (ldiag1 > 0) then
 	icol_diag = -1
 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) icol_diag = i
 	if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
 
 	write(lun,97010) fromwhere, nstep, latndx(i), lonndx(i), k, ipair
 	write(lun,97020) 'drv old/oldbnd/new/del     ',   &
 		dryvol_t_old, dryvol_t_oldbnd, dryvol_t_new, dryvol_t_del
 	write(lun,97020) 'num old/oldbnd, dgnold/new ',   &
 		num_t_old, num_t_oldbnd, dgn_t_old, dgn_t_new
 	write(lun,97020) 'tailfr v_old/new, n_old/new',   &
 		tailfr_volold, tailfr_volnew, tailfr_numold, tailfr_numnew
 	dum = max(1.0d-10,xferfrac_vol) / max(1.0d-10,xferfrac_num)
 	dgn_xfer = dgn_t_new * dum**onethird
 	dum = max(1.0d-10,(1.0d0-xferfrac_vol)) /   &
               max(1.0d-10,(1.0d0-xferfrac_num))
 	dgn_aftr = dgn_t_new * dum**onethird
 	write(lun,97020) 'xferfrac_v/n; dgn_xfer/aftr',   &
 		xferfrac_vol, xferfrac_num, dgn_xfer, dgn_aftr

 97010	format( / 'RENAME ', a, '  nx,lat,lon,k,ip', i8, 4i4 )
 97020	format( a, 6(1pe15.7) )
 	end if
 	end if






	j = jsrflx_rename
	do iq = 1, nspecfrm_renamexf(ipair)
	    xfercoef = xferfrac_vol*deltatinv
	    if (iq .eq. 1) xfercoef = xferfrac_num*deltatinv

	    lsfrma = lspecfrma_renamexf(iq,ipair)-loffset
	    lsfrmc = lspecfrmc_renamexf(iq,ipair)-loffset
	    lstooa = lspectooa_renamexf(iq,ipair)-loffset
	    lstooc = lspectooc_renamexf(iq,ipair)-loffset

	    if (lsfrma .gt. 0) then
		xfertend = xfercoef*max( 0.0_r8,   &
			    (q(i,k,lsfrma)+dqdt(i,k,lsfrma)*deltat) )


                if (ldiag1 > 0) then
                if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
                  if (lstooa .gt. 0) then
                    write(*,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend),   &
                        deltat*dqdt(i,k,lstooa), deltat*(dqdt(i,k,lstooa) + xfertend)
                  else
                    write(*,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend)
                  end if
                end if
                end if



		dqdt(i,k,lsfrma) = dqdt(i,k,lsfrma) - xfertend
		qsrflx(i,lsfrma,j) = qsrflx(i,lsfrma,j) - xfertend*pdel_fac
		if (lstooa .gt. 0) then
		    dqdt(i,k,lstooa) = dqdt(i,k,lstooa) + xfertend
		    qsrflx(i,lstooa,j) = qsrflx(i,lstooa,j) + xfertend*pdel_fac
		end if
	    end if

	    if (lsfrmc .gt. 0) then
		xfertend = xfercoef*max( 0.0_r8,   &
			    (qqcw(i,k,lsfrmc)+dqqcwdt(i,k,lsfrmc)*deltat) )
		dqqcwdt(i,k,lsfrmc) = dqqcwdt(i,k,lsfrmc) - xfertend
		qqcwsrflx(i,lsfrmc,j) = qqcwsrflx(i,lsfrmc,j) - xfertend*pdel_fac
		if (lstooc .gt. 0) then
		    dqqcwdt(i,k,lstooc) = dqqcwdt(i,k,lstooc) + xfertend
		    qqcwsrflx(i,lstooc,j) = qqcwsrflx(i,lstooc,j) + xfertend*pdel_fac
		end if
	    end if

	end do   


	end do mainloop1_ipair


	end do mainloop1_i
	end do mainloop1_k




	dotendrn(:) = .false.
	dotendqqcwrn(:) = .false.
	do ipair = 1, npair_renamexf
	do iq = 1, nspecfrm_renamexf(ipair)
	    lsfrma = lspecfrma_renamexf(iq,ipair) - loffset
	    lsfrmc = lspecfrmc_renamexf(iq,ipair) - loffset
	    lstooa = lspectooa_renamexf(iq,ipair) - loffset
	    lstooc = lspectooc_renamexf(iq,ipair) - loffset
	    if (lsfrma .gt. 0) then
		dotendrn(lsfrma) = .true.
		if (lstooa .gt. 0) dotendrn(lstooa) = .true.
	    end if
	    if (lsfrmc .gt. 0) then
		dotendqqcwrn(lsfrmc) = .true.
		if (lstooc .gt. 0) dotendqqcwrn(lstooc) = .true.
	    end if
	end do
	end do


	return





8100	write(lunout,9050) ipair
	call endrun( 'modal_aero_rename_sub error' )
9050	format( / '*** subr. modal_aero_rename_sub ***' /   &
      	    4x, 'aerosol renaming not implemented for ipair =', i5 )


	end subroutine modal_aero_rename_sub




	subroutine modal_aero_rename_init






	use modal_aero_data
        use module_cam_support, only: &
		pcnst => pcnst_runtime, &
		masterproc
        use constituents, only: cnst_name

	implicit none


	integer ipair, iq, iqfrm, iqfrm_aa, iqtoo, iqtoo_aa,   &
      	  lsfrma, lsfrmc, lstooa, lstooc, lunout,   &
      	  mfrm, mtoo, n1, n2, nsamefrm, nsametoo, nspec


	lunout = 6




	n1 = modeptr_accum
	n2 = modeptr_aitken
	if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
	    npair_renamexf = 1
	    modefrm_renamexf(1) = n2
	    modetoo_renamexf(1) = n1
	else
	    npair_renamexf = 0
	    return
	end if





	do 1900 ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)

	nspec = 0
	do 1490 iqfrm = -1, nspec_amode(mfrm)
	    iqtoo = iqfrm
	    if (iqfrm .eq. -1) then
		lsfrma = numptr_amode(mfrm)
		lstooa = numptr_amode(mtoo)
		lsfrmc = numptrcw_amode(mfrm)
		lstooc = numptrcw_amode(mtoo)
	    else if (iqfrm .eq. 0) then

                goto 1490




	    else
		lsfrma = lmassptr_amode(iqfrm,mfrm)
		lsfrmc = lmassptrcw_amode(iqfrm,mfrm)
		lstooa = 0
		lstooc = 0
	    end if

	    if ((lsfrma .lt. 1) .or. (lsfrma .gt. pcnst)) then
		write(lunout,9100) mfrm, iqfrm, lsfrma
		call endrun( 'modal_aero_rename_init error' )
	    end if
	    if (iqfrm .le. 0) goto 1430

	    if ((lsfrmc .lt. 1) .or. (lsfrmc .gt. pcnst)) then
		write(lunout,9102) mfrm, iqfrm, lsfrmc
		call endrun( 'modal_aero_rename_init error' )
	    end if




	    iqfrm_aa = 1
	    iqtoo_aa = 1
	    if (iqfrm .gt. nspec_amode(mfrm)) then
		iqfrm_aa = nspec_amode(mfrm) + 1
		iqtoo_aa = nspec_amode(mtoo) + 1
	    end if
	    nsamefrm = 0
	    do iq = iqfrm_aa, iqfrm
		if ( lspectype_amode(iq   ,mfrm) .eq.   &
      		     lspectype_amode(iqfrm,mfrm) ) then
		    nsamefrm = nsamefrm + 1
		end if
	    end do
	    nsametoo = 0
	    do iqtoo = iqtoo_aa, nspec_amode(mtoo)
		if ( lspectype_amode(iqtoo,mtoo) .eq.   &
      		     lspectype_amode(iqfrm,mfrm) ) then
		    nsametoo = nsametoo + 1
		    if (nsametoo .eq. nsamefrm) then
			lstooc = lmassptrcw_amode(iqtoo,mtoo)
			lstooa = lmassptr_amode(iqtoo,mtoo)
			goto 1430
		    end if
		end if
	    end do

1430	    nspec = nspec + 1
	    if ((lstooc .lt. 1) .or. (lstooc .gt. pcnst)) lstooc = 0
	    if ((lstooa .lt. 1) .or. (lstooa .gt. pcnst)) lstooa = 0
	    if (lstooa .eq. 0) then
		write(lunout,9104) mfrm, iqfrm, lsfrma, iqtoo, lstooa
		call endrun( 'modal_aero_rename_init error' )
	    end if
	    if ((lstooc .eq. 0) .and. (iqfrm .ne. 0)) then
		write(lunout,9104) mfrm, iqfrm, lsfrmc, iqtoo, lstooc
		call endrun( 'modal_aero_rename_init error' )
	    end if
	    lspecfrma_renamexf(nspec,ipair) = lsfrma
	    lspectooa_renamexf(nspec,ipair) = lstooa
	    lspecfrmc_renamexf(nspec,ipair) = lsfrmc
	    lspectooc_renamexf(nspec,ipair) = lstooc
1490	continue

	nspecfrm_renamexf(ipair) = nspec
1900	continue

9100	format( / '*** subr. modal_aero_rename_init' /   &
      	'lspecfrma out of range' /   &
      	'modefrm, ispecfrm, lspecfrma =', 3i6 / )
9102	format( / '*** subr. modal_aero_rename_init' /   &
      	'lspecfrmc out of range' /   &
      	'modefrm, ispecfrm, lspecfrmc =', 3i6 / )
9104	format( / '*** subr. modal_aero_rename_init' /   &
      	'lspectooa out of range' /   &
      	'modefrm, ispecfrm, lspecfrma, ispectoo, lspectooa =', 5i6 / )
9106	format( / '*** subr. modal_aero_rename_init' /   &
      	'lspectooc out of range' /   &
      	'modefrm, ispecfrm, lspecfrmc, ispectoo, lspectooc =', 5i6 / )




	if ( masterproc ) then

	write(lunout,9310)

	do 2900 ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)
	write(lunout,9320) ipair, mfrm, mtoo

	do iq = 1, nspecfrm_renamexf(ipair)
	    lsfrma = lspecfrma_renamexf(iq,ipair)
	    lstooa = lspectooa_renamexf(iq,ipair)
	    lsfrmc = lspecfrmc_renamexf(iq,ipair)
	    lstooc = lspectooc_renamexf(iq,ipair)
	    if (lstooa .gt. 0) then
		write(lunout,9330) lsfrma, cnst_name(lsfrma),   &
				   lstooa, cnst_name(lstooa)
	    else
		write(lunout,9340) lsfrma, cnst_name(lsfrma)
	    end if
	    if (lstooc .gt. 0) then
		write(lunout,9330) lsfrmc, cnst_name_cw(lsfrmc),   &
				   lstooc, cnst_name_cw(lstooc)
	    else if (lsfrmc .gt. 0) then
		write(lunout,9340) lsfrmc, cnst_name_cw(lsfrmc)
	    else
		write(lunout,9350)
	    end if
	end do

2900	continue
	write(lunout,*)

	end if 

9310	format( / 'subr. modal_aero_rename_init' )
9320	format( 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )
9350	format( 5x, 'no corresponding activated species' )

	return
	end subroutine modal_aero_rename_init



   end module modal_aero_rename
