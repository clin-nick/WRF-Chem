	module module_mosaic_coag3d



	use module_data_mosaic_kind
	use module_peg_util



	implicit none



	contains




	subroutine mosaic_coag_3d_1box( istat_coag,   &
	    idiagbb_in, itstep,   &
	    dtcoag_in,   &
	    temp_box, patm_box, rhoair_box, rbox,   &
	    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
	    adrydens_box, awetdens_box,   &
	    adrydpav_box, awetdpav_box, adryqmas_box, mcoag_flag1 )






	use module_data_mosaic_main, only:  ntot_used, pi, piover6, third
	
        use module_data_mosaic_aero, only:  ifreq_coag 
	use module_data_mosaic_asecthp, only: &
	    ai_phase, hyswptr_aer, &
	    lunerr, lunout, &
	    maxd_acomp, maxd_asize, maxd_atype, &
	    ncomp_aer, ncomp_plustracer_aer, nsize_aer, ntype_aer, &
	    massptr_aer, numptr_aer, waterptr_aer, &
	    dens_aer, dens_water_aer, &
	    dcen_sect, dhi_sect, dlo_sect, &
	    smallmassbb, &
	    volumcen_sect, volumcut_sect, volumhi_sect, volumlo_sect



	integer, intent(inout)  :: istat_coag  
	integer, intent(in)     :: idiagbb_in  
	integer, intent(in)     :: itstep      
        integer, intent(in)     :: mcoag_flag1
	real(r8), intent(in)    :: dtcoag_in   
	real(r8), intent(in)    :: temp_box    
	real(r8), intent(in)    :: patm_box    
	real(r8), intent(in)    :: rhoair_box    

	real(r8), intent(inout) :: rbox(ntot_used)




	real(r8), intent(in) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam








	real(r8), intent(inout) :: adrydens_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: awetdens_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: adrydpav_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: awetdpav_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: adryqmas_box(maxd_asize,maxd_atype)









	integer, parameter :: coag_method = 1
	integer, parameter :: ncomp_coag_maxd = maxd_acomp + 3

	integer :: l, ll, lla, llb, llx
	integer :: isize, itype, iphase
	integer :: iconform_numb
	integer :: idiagbb, idiag_coag_onoff
	integer :: lunout_coag
	integer :: ncomp_coag, nsize, nsubstep_coag, ntau_coag
	integer,save :: ncountaa(10), ncountbb(0:ncomp_coag_maxd)

	real(r8), parameter :: factsafe = 1.00001









	real(r8) :: densdefault
	real(r8) :: dtcoag
	real(r8) :: fact_apdia3, fact_apvolumr
	real(r8) :: press_cgs   
	real(r8) :: tmpa, tmpb
	real(r8) :: wwdens, wwdpav, wwmass, wwvolu
	real(r8) :: xxdens, xxdpav, xxmass, xxnumb, xxvolu
	real(r8) :: xxmasswet, xxvoluwet

	real(r8) :: adryqvol_box(maxd_asize,maxd_atype) 
	real(r8) :: dpdry(maxd_asize,maxd_atype), dpwet(maxd_asize,maxd_atype)
	real(r8) :: densdry(maxd_asize,maxd_atype), denswet(maxd_asize,maxd_atype)
	real(r8) :: num_distrib(maxd_asize,maxd_atype)
	real(r8) :: sumnew(ncomp_coag_maxd), sumold(ncomp_coag_maxd)
	real(r8) :: sum_num(2), sum_vol(0:ncomp_coag_maxd,2)
	real(r8) :: vol_distrib(maxd_asize,maxd_atype,ncomp_coag_maxd)
	real(r8) :: vpcut(0:maxd_asize), vpdry(maxd_asize,maxd_atype)
	real(r8) :: xxvolubb(maxd_asize,maxd_atype)

	character(len=120) :: msg



	if (mcoag_flag1 <= 0) return
	if (ifreq_coag > 1) then
           if (mod(itstep,ifreq_coag) /= 0) return
           dtcoag = dtcoag_in*ifreq_coag
        else
           dtcoag = dtcoag_in
        end if

	istat_coag = 0

	lunout_coag = 6
	if (lunout .gt. 0) lunout_coag = lunout
	write(lunout_coag,'(/a,i10)') 'mosaic_coag_3d_1box - itstep ', itstep



	idiagbb = idiagbb_in
	idiag_coag_onoff = +1
	if (idiagbb <= 0) idiag_coag_onoff = 0

	iphase = ai_phase
	fact_apdia3 = fact_apdiam*fact_apdiam*fact_apdiam
	fact_apvolumr = fact_apmassmr/fact_apdens

	nsubstep_coag = 1
        press_cgs = patm_box*1.01325e6_r8


	ncountaa(:) = 0   
	ncountbb(:) = 0



	itype = 1
	nsize = nsize_aer(itype)
	ncomp_coag = ncomp_plustracer_aer(itype) + 3
	densdefault = dens_aer(1,itype)*fact_apdens
	vpcut(0:nsize) = volumcut_sect(0:nsize,itype)*fact_apdia3

	ncountaa(1) = ncountaa(1) + 1



	vol_distrib(:,:,:) = 0.0
	sumold(:) = 0.0

itype_loop_aa: &
	do itype = 1, ntype_aer

isize_loop_aa: &
	do isize = 1, nsize
	    l = numptr_aer(isize,itype,iphase)
	    xxnumb = rbox(l)*fact_apnumbmr
	    xxmass = adryqmas_box(isize,itype)*fact_apmassmr
	    xxdens = adrydens_box(isize,itype)*fact_apdens
	    iconform_numb = 1

	    if ( (xxdens .lt. 0.1) .or. (xxdens .gt. 20.0) ) then


		ncountaa(2) = ncountaa(2) + 1
		xxvolu = 0.0
		if (idiagbb .ge. 200)   &
		    write(lunout_coag,'(a,i10,2i4,1p,4e10.2)') 'coagaa-2a',   &
		    itstep, itype, isize, xxmass, xxvolu, xxdens
		xxmass = 0.0
		do ll = 1, ncomp_aer(itype)
		    l = massptr_aer(ll,isize,itype,iphase)
		    if (l .ge. 1) then
			tmpa = max( 0.0_r8, rbox(l) )
			xxmass = xxmass + tmpa
			xxvolu = xxvolu + tmpa/dens_aer(ll,itype)
		    end if
		end do
		xxmass = xxmass*fact_apmassmr
		xxvolu = xxvolu*fact_apvolumr
		if (xxmass .gt. 0.99*smallmassbb) then
		    xxdens = xxmass/xxvolu
		    xxvolu = xxmass/xxdens
		    if (idiagbb .ge. 200)   &
			write(lunout_coag,'(a,i10,2i4,1p,4e10.2)') 'coagaa-2b',   &
			itstep, itype, isize, xxmass, xxvolu, xxdens
		end if
	    else
		xxvolu = xxmass/xxdens
	    end if

	    if (xxmass .le. 1.01*smallmassbb) then


		ncountaa(3) = ncountaa(3) + 1
		xxdens = densdefault
		xxvolu = xxmass/xxdens
		xxnumb = xxmass/(volumcen_sect(isize,itype)*fact_apdia3*xxdens)
		l = hyswptr_aer(isize,itype)
		if (l .ge. 1) rbox(l) = 0.0
		l = waterptr_aer(isize,itype)
		if (l .ge. 1) rbox(l) = 0.0
		iconform_numb = 0
		if (idiagbb .ge. 210)   &
		    write(lunout_coag,'(a,i10,2i4,1p,4e10.2)') 'coagaa-2c',   &
		    itstep, itype, isize, xxmass, xxvolu, xxdens
	    else
		xxvolu = xxmass/xxdens
	    end if


	    if (iconform_numb .gt. 0) then
		if (xxnumb .gt.   &
			xxvolu/(factsafe*volumlo_sect(isize,itype)*fact_apdia3)) then
		    ncountaa(4) = ncountaa(4) + 1
		    tmpa = xxnumb
		    xxnumb = xxvolu/(factsafe*volumlo_sect(isize,itype)*fact_apdia3)
		    if (idiagbb .ge. 200)   &
			write(lunout_coag,'(a,i10,2i4,1p,3e12.4)') 'coagcc-4 ',   &
			itstep, itype, isize, xxvolu, tmpa, xxnumb
		else if (xxnumb .lt.   &
			xxvolu*factsafe/(volumhi_sect(isize,itype)*fact_apdia3)) then
		    ncountaa(5) = ncountaa(5) + 1
		    tmpa = xxnumb
		    xxnumb = xxvolu*factsafe/(volumhi_sect(isize,itype)*fact_apdia3)
		    if (idiagbb .ge. 200)   &
			write(lunout_coag,'(a,i10,2i4,1p,3e12.4)') 'coagcc-5 ',   &
			itstep, itype, isize, xxvolu, tmpa, xxnumb
		    if (idiagbb .ge. 200)   &
			write(lunout_coag,'(a,10x, 8x,1p,4e12.4)') 'coagcc-5 ',   &
			dlo_sect(isize,itype), dlo_sect(isize,itype)*fact_apdiam, &
			volumlo_sect(isize,itype), volumlo_sect(isize,itype)*fact_apdia3
		    if (idiagbb .ge. 200)   &
			write(lunout_coag,'(a,10x, 8x,1p,4e12.4)') 'coagcc-5 ',   &
			dhi_sect(isize,itype), dhi_sect(isize,itype)*fact_apdiam, &
			volumhi_sect(isize,itype), volumhi_sect(isize,itype)*fact_apdia3
		end if
	    end if


	    l = numptr_aer(isize,itype,iphase)
	    rbox(l) = xxnumb/fact_apnumbmr
	    adrydens_box(isize,itype) = xxdens/fact_apdens
	    adryqmas_box(isize,itype) = xxmass/fact_apmassmr
	    adryqvol_box(isize,itype) = xxvolu/fact_apvolumr










	    num_distrib(isize,itype) = xxnumb*rhoair_box

	    do ll = 1, ncomp_plustracer_aer(itype)
		l = massptr_aer(ll,isize,itype,iphase)
		tmpa = 0.0
		if (l .ge. 1) tmpa = max( 0.0_r8, rbox(l) )*fact_apmassmr
		vol_distrib(isize,itype,ll)   = tmpa
		sumold(ll) = sumold(ll) + tmpa
	    end do

	    do llx = 1, 3
		ll = (ncomp_coag-3) + llx
		tmpa = 0.0
		if (llx .eq. 1) then
		    l = hyswptr_aer(isize,itype)
		    if (l .ge. 1) tmpa = max( 0.0_r8, rbox(l) )*fact_apmassmr
		else if (llx .eq. 2) then
		    l = waterptr_aer(isize,itype)
		    if (l .ge. 1) tmpa = max( 0.0_r8, rbox(l) )*fact_apmassmr
		else
		    tmpa = max( 0.0_r8, adryqvol_box(isize,itype) )*fact_apvolumr
		end if
		vol_distrib(isize,itype,ll)   = tmpa
		sumold(ll) = sumold(ll) + tmpa
	    end do


	    if (xxmass .le. 1.01*smallmassbb) then
		dpdry(isize,itype) = dcen_sect(isize,itype)*fact_apdiam
		dpwet(isize,itype) = dpdry(isize,itype)
		denswet(isize,itype) = xxdens
	    else
		dpdry(isize,itype) = (xxvolu/(xxnumb*piover6))**third
		dpdry(isize,itype) = max( dpdry(isize,itype), dlo_sect(isize,itype)*fact_apdiam )
		l = waterptr_aer(isize,itype)
		if (l .ge. 1) then
		    tmpa = max( 0.0_r8, rbox(l) )*fact_apmassmr
		    xxmasswet = xxmass + tmpa
		    xxvoluwet = xxvolu + tmpa/(dens_water_aer*fact_apdens)
		    dpwet(isize,itype) = (xxvoluwet/(xxnumb*piover6))**third
		    dpwet(isize,itype) = min( dpwet(isize,itype), 30.0_r8*dhi_sect(isize,itype)*fact_apdiam )
		    denswet(isize,itype) = (xxmasswet/xxvoluwet)
		else
		    dpwet(isize,itype) = dpdry(isize,itype)
		    denswet(isize,itype) = xxdens
		end if
	    end if
	    densdry(isize,itype) = xxdens

	    vpdry(isize,itype) = piover6 * (dpdry(isize,itype)**3)

	end do isize_loop_aa

	end do itype_loop_aa





	if (idiag_coag_onoff > 0) call coag_3d_conserve_check(   &
	    nsize, maxd_asize, ntype_aer, maxd_atype, ncomp_coag, ncomp_coag_maxd, 1, lunout_coag,   &
	    dpdry, num_distrib, vol_distrib, sum_num, sum_vol )

	write(*,'(/a/)') 'mosaic_coag_3d_1box calling movingcenter_coag_3d'
	call movingcenter_coag_3d(   &
		idiagbb, lunout_coag, nsize, maxd_asize, ntype_aer, maxd_atype, iphase, &
		ncomp_coag, ncomp_coag_maxd,   &
		temp_box, press_cgs, dtcoag, nsubstep_coag,   &
		vpcut, vpdry, dpdry, dpwet, densdry, denswet, num_distrib, vol_distrib )
	write(*,'(/a/)') 'mosaic_coag_3d_1box backfrm movingcenter_coag_3d'

	if (idiag_coag_onoff > 0) call coag_3d_conserve_check(   &
	    nsize, maxd_asize, ntype_aer, maxd_atype, ncomp_coag, ncomp_coag_maxd, 2, lunout_coag,   &
	    dpdry, num_distrib, vol_distrib, sum_num, sum_vol )





	sumnew(:) = 0.0

itype_loop_yy: &
	do itype = 1, ntype_aer

isize_loop_xx: &
	do isize = 1, nsize
	    do ll = 1, ncomp_coag
		sumnew(ll) = sumnew(ll) + max( 0.0_r8, vol_distrib(isize,itype,ll) )
	    end do

	    l = numptr_aer(isize,itype,iphase)
	    rbox(l) = max( 0.0_r8, num_distrib(isize,itype)/(rhoair_box*fact_apnumbmr) )


	    xxmass = 0.0
	    xxvolu = 0.0
	    do ll = 1, ncomp_aer(itype)
		l = massptr_aer(ll,isize,itype,iphase)
		if (l .ge. 1) then
		    tmpa = max( 0.0_r8, vol_distrib(isize,itype,ll) )
		    rbox(l) = tmpa/fact_apmassmr
		    xxmass = xxmass + tmpa
		    xxvolu = xxvolu + tmpa/(dens_aer(ll,itype)*fact_apdens)
		end if
	    end do
	    adryqmas_box(isize,itype) = xxmass/fact_apmassmr
	    xxvolubb(isize,itype) = xxvolu

	    ll = (ncomp_coag-3) + 1
	    l = hyswptr_aer(isize,itype)
	    if (l .ge. 1) rbox(l) = max( 0.0_r8, vol_distrib(isize,itype,ll) )/fact_apmassmr

	    ll = (ncomp_coag-3) + 2
	    l = waterptr_aer(isize,itype)
	    if (l .ge. 1) rbox(l) = max( 0.0_r8, vol_distrib(isize,itype,ll) )/fact_apmassmr

	    ll = (ncomp_coag-3) + 3
	    adryqvol_box(isize,itype) = max( 0.0_r8, vol_distrib(isize,itype,ll) )/fact_apvolumr
	end do isize_loop_xx






isize_loop_yy: &
	do isize = 1, nsize

	    xxmass = adryqmas_box(isize,itype)*fact_apmassmr
	    xxvolu = adryqvol_box(isize,itype)*fact_apvolumr
	    l = numptr_aer(isize,itype,iphase)
	    xxnumb = rbox(l)*fact_apnumbmr
	    iconform_numb = 1


	    if (xxmass .le. smallmassbb) then
		ncountaa(8) = ncountaa(8) + 1
		xxdens = densdefault
		xxvolu = xxmass/xxdens
		xxnumb = xxmass/(volumcen_sect(isize,itype)*fact_apdia3*xxdens)
		xxdpav = dcen_sect(isize,itype)*fact_apdiam
		wwdens = xxdens
		wwdpav = xxdpav
		l = hyswptr_aer(isize,itype)
		if (l .ge. 1) rbox(l) = 0.0
		l = waterptr_aer(isize,itype)
		if (l .ge. 1) rbox(l) = 0.0
		iconform_numb = 0
		if (idiagbb .ge. 210)   &
		    write(lunout_coag,'(a,i10,2i4,1p,4e10.2)') 'coagaa-7a',   &
		    itstep, itype, isize, xxmass, xxvolu, xxdens
	    else if (xxmass .gt. 1000.0*xxvolu) then


		xxdens = 1000.0
	    else 
		xxdens = xxmass/xxvolu
	    end if

	    if ((xxdens .lt. 0.1) .or. (xxdens .gt. 20.0)) then


		ncountaa(7) = ncountaa(7) + 1
		if (idiagbb .ge. 200)   &
		    write(lunout_coag,'(a,i10,2i4,1p,4e10.2)') 'coagaa-7b',   &
		    itstep, itype, isize, xxmass, xxvolu, xxdens
		xxvolu = xxvolubb(isize,itype)
		xxdens = xxmass/xxvolu
		if (idiagbb .ge. 200)   &
		    write(lunout_coag,'(a,18x,1p,4e10.2)') 'coagaa-7c',   &
		    xxmass, xxvolu, xxdens
	    end if


	    if (iconform_numb .gt. 0) then
		if (xxnumb .gt. xxvolu/(volumlo_sect(isize,itype)*fact_apdia3)) then
		    ncountaa(9) = ncountaa(9) + 1
		    tmpa = xxnumb
		    xxnumb = xxvolu/(volumlo_sect(isize,itype)*fact_apdia3)
		    xxdpav = dlo_sect(isize,itype)*fact_apdiam
		    if (idiagbb .ge. 200)   &
			write(lunout_coag,'(a,i10,2i4,1p,3e12.4)') 'coagcc-9 ',   &
			itstep, itype, isize, xxvolu, tmpa, xxnumb
		else if (xxnumb .lt. xxvolu/(volumhi_sect(isize,itype)*fact_apdia3)) then
		    ncountaa(10) = ncountaa(10) + 1
		    tmpa = xxnumb
		    xxnumb = xxvolu/(volumhi_sect(isize,itype)*fact_apdia3)
		    xxdpav = dhi_sect(isize,itype)*fact_apdiam
		    if (idiagbb .ge. 200)   &
			write(lunout_coag,'(a,i10,2i4,1p,3e12.4)') 'coagcc-10',   &
			itstep, itype, isize, xxvolu, tmpa, xxnumb
		else
		    xxdpav = (xxvolu/(xxnumb*piover6))**third
		end if

		tmpb = 0.0
		l = waterptr_aer(isize,itype)
		if (l .ge. 1) tmpb = max(0.0_r8,rbox(l))*fact_apmassmr
		wwmass = xxmass + tmpb
		wwvolu = xxvolu + tmpb/(dens_water_aer*fact_apdens)
		wwdens = wwmass/wwvolu
		wwdpav = xxdpav*((wwvolu/xxvolu)**third)
	    end if


	    l = numptr_aer(isize,itype,iphase)
	    rbox(l) = xxnumb/fact_apnumbmr
	    adrydens_box(isize,itype) = xxdens/fact_apdens
	    adrydpav_box(isize,itype) = xxdpav/fact_apdiam
	    adryqmas_box(isize,itype) = xxmass/fact_apmassmr
	    adryqvol_box(isize,itype) = xxvolu/fact_apvolumr
	    awetdens_box(isize,itype) = wwdens/fact_apdens
	    awetdpav_box(isize,itype) = wwdpav/fact_apdiam

	    dpdry(isize,itype) = xxdpav

	end do isize_loop_yy

	end do itype_loop_yy


	if (idiag_coag_onoff > 0) call coag_3d_conserve_check(   &
	    nsize, maxd_asize, ntype_aer, maxd_atype, ncomp_coag, ncomp_coag_maxd, 3, lunout_coag,   &
	    dpdry, num_distrib, vol_distrib, sum_num, sum_vol )



	do ll = 1, ncomp_coag
	    
	    tmpa = max( sumold(ll), sumnew(ll), 1.0e-35_r8 )
	    if (abs(sumold(ll)-sumnew(ll)) .gt. 1.0e-6*tmpa) then
		ncountbb(ll) = ncountbb(ll) + 1
		ncountbb(0) = ncountbb(0) + 1
	    end if
	end do

	if (idiagbb .ge. 100) then
	    write(msg,93030) 'coagbb ncntaa ', ncountaa(1:10)
	    call peg_message( lunerr, msg )
	    if (ncountbb(0) .gt. 0) then
		do llx = 1, (ncomp_coag+9)/10
		    llb = llx*10
		    lla = llb - 9
		    llb = min( llb, ncomp_coag)
		    write(msg,93032) 'coagbb ncntbb',   &
			mod(llx,10), ncountbb(lla:llb)
		    call peg_message( lunerr, msg )
		end do
	    end if
	end if
93020	format( a, 1p, 10e10.2 )
93030	format( a, 1p, 10i10 )
93032	format( a, 1p, i1, 10i10 )


	return
	end subroutine mosaic_coag_3d_1box



      subroutine coag_3d_conserve_check(   &
        nbin, nbin_maxd, ntype, ntype_maxd, ncomp, ncomp_maxd, iflagaa, lunout,   &
        dpdry, num_distrib, vol_distrib, sum_num, sum_vol )

      implicit none



      integer, intent(in) :: nbin            
      integer, intent(in) :: nbin_maxd       
      integer, intent(in) :: ntype           
      integer, intent(in) :: ntype_maxd      
      integer, intent(in) :: ncomp           
      integer, intent(in) :: ncomp_maxd      
      integer, intent(in) :: lunout          
      integer, intent(in) :: iflagaa         

      real(r8), intent(in) :: dpdry(nbin_maxd,ntype_maxd)
      real(r8), intent(in) :: num_distrib(nbin_maxd,ntype_maxd)
      real(r8), intent(in) :: vol_distrib(nbin_maxd,ntype_maxd,ncomp_maxd)
      real(r8), intent(inout) :: sum_num(2)
      real(r8), intent(inout) :: sum_vol(0:ncomp_maxd,2)


      integer :: i, itype, l
      real(r8), parameter :: pi = 3.14159265358979323846_r8
      real(r8) :: voldry(nbin,ntype)


      if (iflagaa == 1) then
         i = 1
      else if (iflagaa == 2) then
         i = 2
      else if (iflagaa == 3) then
         i = 2
      else
         return
      end if



      voldry(1:nbin,1:ntype) = (pi/6.0)*(dpdry(1:nbin,1:ntype)**3)
      sum_num(i) = sum( num_distrib(1:nbin,1:ntype) )
      sum_vol(0,i) = sum( num_distrib(1:nbin,1:ntype)*voldry(1:nbin,1:ntype) )
      do l = 1, ncomp
         sum_vol(l,i) = sum( vol_distrib(1:nbin,1:ntype,l) )
      end do


      if ((iflagaa /= 2) .and. (iflagaa /= 3)) return
      write(lunout,'(/a,i3)') 'coag_3d_conserve_check -- iflagaa = ', iflagaa

      write(lunout,'(2a,1p,2e14.6)')   &
               'number total ', '    initial, final   =',   &
               sum_num(1:2)
      write(lunout,'(2a,1p,2e14.6)')   &
               'number total ', '   (final/initial)   =',   &
               sum_num(2)/sum_num(1)







      if (iflagaa == 3) then
      write(lunout,'(2a,1p,2e14.6)')   &
               'volume total ', '    initial, final   =',   &
               sum_vol(0,1:2)
      write(lunout,'(2a,1p,2e14.6)')   &
               'volume total ', '   (initial/final)-1 =',   &
               (sum_vol(0,1)/sum_vol(0,2))-1.0
      end if

      do l = 1, ncomp
         if (abs(sum_vol(l,2)) > 0.0) then
            write(lunout,'(a,i4,a,1p,2e14.6)')   &
               'volume l=',  l, '   (initial/final)-1 =',   &
               (sum_vol(l,1)/sum_vol(l,2))-1.0
         else
            write(lunout,'(a,i4,a,1p,2e14.6)')   &
               'volume l=',  l, '    initial, final   =',   &
               sum_vol(l,1), sum_vol(l,2)
         end if
      end do

      return
      end subroutine coag_3d_conserve_check




      subroutine movingcenter_coag_3d(   &
        idiagbb, lunout, nbin, nbin_maxd, ntype, ntype_maxd, iphase, &
        ncomp, ncomp_maxd,   &
        tempk, press, deltat, nsubstep_in,   &
        vpcut, vpdry_in, dpdry, dpwet, densdry, denswet,   &
        num_distrib, vol_distrib )

      use module_data_mosaic_aero, only:  &
         method_bcfrac, method_kappa, msectional_flag2
      use module_data_mosaic_asecthp, only:  &
         dens_aer, hygro_aer, &
         itype_md1_of_itype, itype_md2_of_itype, itype_of_itype_md1md2, &
         lptr_bc_aer, massptr_aer, &
         ncomp_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer, &
         smallmassbb, xcut_atype_md1, xcut_atype_md2

      implicit none













      integer, intent(in) :: idiagbb
      integer, intent(in) :: lunout
      integer, intent(in) :: nbin            
      integer, intent(in) :: nbin_maxd       
      integer, intent(in) :: ntype           
      integer, intent(in) :: ntype_maxd      
      integer, intent(in) :: iphase
      integer, intent(in) :: ncomp           
      integer, intent(in) :: ncomp_maxd      
      integer, intent(in) :: nsubstep_in     

      real(r8), intent(in) :: tempk               
      real(r8), intent(in) :: press               
      real(r8), intent(in) :: deltat              
      real(r8), intent(in) :: dpdry(nbin_maxd,ntype_maxd)    
      real(r8), intent(in) :: dpwet(nbin_maxd,ntype_maxd)    
      real(r8), intent(in) :: densdry(nbin_maxd,ntype_maxd)  
      real(r8), intent(in) :: denswet(nbin_maxd,ntype_maxd)  
      real(r8), intent(in) :: vpdry_in(nbin_maxd,ntype_maxd) 
      real(r8), intent(in) :: vpcut(0:nbin_maxd)  


      real(r8), intent(inout) :: num_distrib(nbin_maxd,ntype_maxd)

      real(r8), intent(inout) :: vol_distrib(nbin_maxd,ntype_maxd,ncomp_maxd)





















      character(len=256) :: errmsg
      integer :: i, isubstep, itype, itype1, itype2, iwrite193
      integer :: j, jtype, jtype1, jtype2
      integer :: k, ktype, ktype1, ktype2, ktmp1, ktmp2
      integer :: l, method_ktype, nsubstep

      real(r8), parameter :: frac_loss_limit = 0.8_r8
      real(r8), parameter :: pi = 3.14159265358979323846_r8

      real(r8) :: del_num, del_voli, del_volj, deltat_sub
      real(r8) :: tmp_bcfrac, tmp_beta
      real(r8) :: tmp_densbc
      real(r8) :: tmp_kappa, tmp_kappaxvolu
      real(r8) :: tmp_massbc, tmp_massdry
      real(r8) :: tmp_volu
      real(r8) :: tmpa, tmpb, tmpc, tmpi, tmpj
      real(r8) :: vpdry_ipj

      real(r8) ::   &
             betaxh(nbin,nbin),   &
             bcfrac_bin(nbin,ntype),   &
             cnum(nbin,ntype), cnum_old(nbin,ntype),   &
             cvol(nbin,ntype,ncomp), cvol_old(nbin,ntype,ncomp),   &
             kappa_bin(nbin,ntype),   &
             tmpvecb(nbin),   &
             vpdry(nbin,ntype)



      method_ktype = 0
      if ((ntype > 1) .and. (msectional_flag2 > 0)) then
         if (ncomp /= ncomp_aer(1) + 3) then
            write(errmsg,'(/2a,2(1x,i5)/)') '*** movingcenter_coag_3d - ', &
               'mismatch of ncomp & ncomp_aer(1)', ncomp, ncomp_aer(1)
            call wrf_error_fatal3("<stdin>",742,&
trim(adjustl(errmsg)))
         end if
         if ((method_bcfrac /= 1) .and. (method_bcfrac /= 11)) then
            write(errmsg,'(/2a,2(1x,i5)/)') '*** movingcenter_coag_3d - ', &
               'method_bcfrac must be 1 or 11 but is', method_bcfrac
            call wrf_error_fatal3("<stdin>",748,&
trim(adjustl(errmsg)))
         end if
         if ((method_kappa /= 11) .and. (method_kappa /= 12)) then
            write(errmsg,'(/2a,2(1x,i5)/)') '*** movingcenter_coag_3d - ', &
               'method_kappa must be 11 or 12 but is', method_kappa
            call wrf_error_fatal3("<stdin>",754,&
trim(adjustl(errmsg)))
         end if
         method_ktype = 1
      end if

      iwrite193 = 0


      tmp_densbc = -2.0
      do itype = 1, ntype
      do l = 1, ncomp_aer(itype)
         if (massptr_aer(l,1,itype,iphase) == lptr_bc_aer(1,itype,iphase)) then
            tmp_densbc = dens_aer(l,itype)
         end if
      end do
      end do
      if (tmp_densbc < -1.0) then
         call wrf_error_fatal3("<stdin>",772,&
'*** movingcenter_coag_3d - cannot find bc')        
      end if





      nsubstep = nsubstep_in
      deltat_sub = deltat/nsubstep   
      tmpa = 0.0
      do jtype = 1, ntype
         tmpvecb(:) = 0.0

         do itype = 1, ntype

         call brownian_kernel_3dsub(   &
            nbin, nbin_maxd, ntype, ntype_maxd, itype, jtype, &
            tempk, press, dpwet, denswet, betaxh )

         do j = 1, nbin
         do i = 1, nbin
            if ((i /= j) .or. (itype /= jtype)) then
               tmpvecb(j) = tmpvecb(j) + betaxh(i,j)*num_distrib(i,itype)
            else
               tmpvecb(j) = tmpvecb(j) + betaxh(i,j)*num_distrib(i,itype)*0.5
            end if
         end do   
         end do   
         end do   

         do j = 1, nbin
            tmpa = max( tmpa, tmpvecb(j) )
         end do   
      end do   

      if (idiagbb > 0) write(lunout,'(/a,i10,1p,4e12.4)') 'coag aa nsub, dtsub, fracloss', &
            nsubstep, deltat_sub, tmpa*deltat_sub, frac_loss_limit
      if (tmpa*deltat_sub > frac_loss_limit) then
         tmpb = tmpa*deltat/frac_loss_limit
         isubstep = int(tmpb) + 1
         tmpb = deltat/isubstep
         if (idiagbb > 0) write(lunout,'( a,i10,1p,4e12.4)') 'coag bb nsub, dtsub, fracloss', &
            isubstep, tmpb, tmpa*tmpb
         nsubstep = isubstep
         deltat_sub = deltat/nsubstep
      end if





      cnum(1:nbin,1:ntype) = num_distrib(1:nbin,1:ntype)
      cvol(1:nbin,1:ntype,1:ncomp) = vol_distrib(1:nbin,1:ntype,1:ncomp)





solve_isubstep_loop:   &
      do isubstep = 1, nsubstep




         cnum_old(1:nbin,1:ntype) = cnum(1:nbin,1:ntype) 
         cvol_old(1:nbin,1:ntype,1:ncomp) = cvol(1:nbin,1:ntype,1:ncomp)


         if (isubstep == 1) then
             vpdry(1:nbin,1:ntype) = vpdry_in(1:nbin,1:ntype)







         end if

         if (method_ktype == 1) then

            do itype = 1, ntype
            do i = 1, nbin
               tmp_massbc  = 0.0
               tmp_massdry = 0.0
               tmp_kappaxvolu  = 0.0
               tmp_volu = 0.0
               do l = 1, ncomp_aer(itype)
                  tmpa = vol_distrib(i,itype,l)
                  tmpc = tmpa/dens_aer(l,itype)
                  if (method_bcfrac == 11) then
                     
                     tmpb = tmpc
                  else
                     
                     tmpb = tmpa
                  end if
                  tmp_massdry = tmp_massdry + tmpb


                  if (massptr_aer(l,i,itype,iphase) == lptr_bc_aer(i,itype,iphase)) then
                     tmp_massbc = tmpb
                     if (method_kappa == 12) cycle
                     
                  end if
                  tmp_volu = tmp_volu + tmpc
                  tmp_kappaxvolu = tmp_kappaxvolu + tmpc*hygro_aer(l,itype)
               end do   

               if (tmp_massdry <= smallmassbb) then
                  itype1 = itype_md1_of_itype( itype )
                  bcfrac_bin(i,itype) = 0.5*( xcut_atype_md1(itype1-1) + xcut_atype_md1(itype1) )
               else
                  bcfrac_bin(i,itype) = tmp_massbc/tmp_massdry
               end if
               bcfrac_bin(i,itype) = max( 0.0_r8, xcut_atype_md1(0), &
                                          bcfrac_bin(i,itype) )
               bcfrac_bin(i,itype) = min( 1.0_r8, xcut_atype_md1(ntype_md1_aer), &
                                          bcfrac_bin(i,itype) )

               if (tmp_volu <= smallmassbb) then
                  itype2 = itype_md2_of_itype( itype )
                  kappa_bin(i,itype)  = sqrt( max( 0.0_r8, &
                                              xcut_atype_md2(itype2-1) * xcut_atype_md2(itype2) ) )
               else
                  kappa_bin(i,itype) = tmp_kappaxvolu/tmp_volu
               end if
               kappa_bin(i,itype) = max( 1.0e-10_r8, xcut_atype_md2(0), &
                                          kappa_bin(i,itype) )
               kappa_bin(i,itype) = min( 2.0_r8,     xcut_atype_md2(ntype_md2_aer), &
                                          kappa_bin(i,itype) )
            end do   
            end do   
         end if   



         do jtype = 1, ntype
         jtype1 = itype_md1_of_itype( jtype )
         jtype2 = itype_md2_of_itype( jtype )

         do itype = 1, jtype   
         itype1 = itype_md1_of_itype( itype )
         itype2 = itype_md2_of_itype( itype )



         call brownian_kernel_3dsub(   &
            nbin, nbin_maxd, ntype, ntype_maxd, itype, jtype, &
            tempk, press, dpwet, denswet, betaxh )

         if (iwrite193 > 0) &
            write(193,'(/a,2(4x,2i4))') 'i/jtype1 then ...2', itype1, jtype1, itype2, jtype2
         do j = 1, nbin
         do i = 1, nbin
            if (itype == jtype) then
               if (i > j) cycle   
            end if

            
            vpdry_ipj = vpdry(i,itype) + vpdry(j,jtype)
            k = max( i, j )
            do while (k > -99)
               if (vpdry_ipj > vpcut(k)) then
                  k = k+1
                  if (k >= nbin) exit   
               else if (vpdry_ipj < vpcut(k-1)) then
                  k = k-1
                  if (k <= 1) exit   
               else
                  exit   
               end if
            end do   
            k = max( 1, min( nbin, k ) )

            if (method_ktype == 1) then
               if (method_bcfrac == 11) then
                  tmpi = vpdry(i,itype)
                  tmpj = vpdry(j,jtype)
               else
                  tmpi = vpdry(i,itype)*densdry(i,itype)
                  tmpj = vpdry(j,jtype)*densdry(j,jtype)
               end if
               tmpa = tmpi/max( (tmpi + tmpj), 1.0e-35_r8 )
               tmpa = max( 0.0_r8, min( 1.0_r8, tmpa ) )
               tmp_bcfrac = bcfrac_bin(i,itype)*tmpa + bcfrac_bin(j,jtype)*(1.0_r8 - tmpa)

               ktype1 = ntype_md1_aer
               do ktmp1 = 1, ntype_md1_aer - 1
                   if (tmp_bcfrac <= xcut_atype_md1(ktmp1)) then
                       ktype1 = ktmp1
                       exit
                   end if
               end do

               if (method_kappa == 11) then
                  tmpi = vpdry(i,itype)
                  tmpj = vpdry(j,jtype)
               else if (method_bcfrac == 11) then
                  tmpa = 1.0_r8 - bcfrac_bin(i,itype)
                  tmpi = vpdry(i,itype)*max(tmpa,1.0e-10_r8)
                  tmpa = 1.0_r8 - bcfrac_bin(j,jtype)
                  tmpj = vpdry(j,jtype)*max(tmpa,1.0e-10_r8)
               else
                  tmpa = 1.0_r8 - bcfrac_bin(i,itype)*(densdry(i,itype)/tmp_densbc)
                  tmpi = vpdry(i,itype)*max(tmpa,1.0e-10_r8)
                  tmpa = 1.0_r8 - bcfrac_bin(j,jtype)*(densdry(j,jtype)/tmp_densbc)
                  tmpj = vpdry(j,jtype)*max(tmpa,1.0e-10_r8)
               end if
               tmpa = tmpi/max( (tmpi + tmpj), 1.0e-35_r8 )
               tmpa = max( 0.0_r8, min( 1.0_r8, tmpa ) )
               tmp_kappa = kappa_bin(i,itype)*tmpa + kappa_bin(j,jtype)*(1.0_r8 - tmpa)

               ktype2 = ntype_md2_aer
               do ktmp2 = 1, ntype_md2_aer - 1
                   if (tmp_kappa <= xcut_atype_md2(ktmp2)) then
                       ktype2 = ktmp2
                       exit
                   end if
               end do

               ktype = itype_of_itype_md1md2(ktype1,ktype2)

               if (iwrite193 > 0) then




















               if ( ( j==6 .or. j==12 .or. j==18 .or. j==24 ) .and. &
                    ( abs(i-j)<=4 .or. abs(i-j)==8 ) .and. &
                    ( jtype1==1 .or. jtype1==5 .or. jtype1==10 .or. jtype1==15 .or. &
                      itype1==1 .or. itype1==5 .or. itype1==10 .or. itype1==15 ) .and. &
                    ( abs(jtype1-itype1)<=1 .or. abs(jtype1-itype1)==5 .or. abs(jtype1-itype1)==10 ) .and. &
                    ( jtype2==1 .or. jtype2==5 .or. jtype2==9 .or. &
                      itype2==1 .or. itype2==5 .or. itype2==9 ) .and. &
                    ( abs(jtype2-itype2)<=1 .or. abs(jtype2-itype2)==4 .or. abs(jtype2-itype2)==8 ) ) then
                  write(193,'(a,3i3,1p,5e9.2,0p,2f6.3,3i3,3f6.3,3i3,3f7.4)') 'ijk&types', &
                     i, j, k, dpdry(i,itype), dpdry(j,jtype), vpdry(i,itype), vpdry(j,jtype), vpdry_ipj, &
                     densdry(i,itype), densdry(j,jtype), &
                     itype1, jtype1, ktype1, bcfrac_bin(i,itype), bcfrac_bin(j,jtype), tmp_bcfrac, &
                     itype2, jtype2, ktype2, kappa_bin(i,itype), kappa_bin(j,jtype), tmp_kappa
               end if

               end if  

            else
               if (itype == jtype) then
                  ktype = itype
               else
                  ktype = ntype
               end if
            end if



            tmp_beta = betaxh(i,j) * deltat_sub

            
            if ((i == j) .and. (itype == jtype)) then
               del_num = 0.5*tmp_beta*cnum_old(i,itype)*cnum_old(i,itype)
               if ((i == k) .and. (itype == ktype)) then
                  
                  
                  cnum(i,itype) = cnum(i,itype) - del_num
                  cycle
               else
                  
                  
                  
                  cnum(i,itype) = cnum(i,itype) - del_num*2.0
                  cnum(k,ktype) = cnum(k,ktype) + del_num
                  do l = 1, ncomp
                     del_voli = tmp_beta*cnum_old(i,itype)*cvol_old(i,itype,l)
                     cvol(i,itype,l) = cvol(i,itype,l) - del_voli
                     cvol(k,ktype,l) = cvol(k,ktype,l) + del_voli
                  end do
               end if

            else   
               del_num = tmp_beta*cnum_old(i,itype)*cnum_old(j,jtype)
               cnum(i,itype) = cnum(i,itype) - del_num
               if ((j == k) .and. (jtype == ktype)) then
                  
                  
                  
                  do l = 1, ncomp
                     del_voli = tmp_beta*cnum_old(j,jtype)*cvol_old(i,itype,l)
                     cvol(i,itype,l) = cvol(i,itype,l) - del_voli
                     cvol(k,ktype,l) = cvol(k,ktype,l) + del_voli
                  end do
               else
                  
                  
                  cnum(j,jtype) = cnum(j,jtype) - del_num
                  cnum(k,ktype) = cnum(k,ktype) + del_num
                  do l = 1, ncomp
                     del_voli = tmp_beta*cnum_old(j,jtype)*cvol_old(i,itype,l)
                     del_volj = tmp_beta*cnum_old(i,itype)*cvol_old(j,jtype,l)
                     cvol(i,itype,l) = cvol(i,itype,l) - del_voli
                     cvol(j,jtype,l) = cvol(j,jtype,l) - del_volj
                     cvol(k,ktype,l) = cvol(k,ktype,l) + del_voli + del_volj
                  end do
               end if

            end if

         end do   
         end do   

         end do   
         end do   

      end do solve_isubstep_loop



      num_distrib(1:nbin,1:ntype) = cnum(1:nbin,1:ntype)
      do l = 1, ncomp
         vol_distrib(1:nbin,1:ntype,l) = cvol(1:nbin,1:ntype,l)
      end do



      return
      end subroutine movingcenter_coag_3d




      subroutine brownian_kernel_3dsub(   &
         nbin, nbin_maxd, ntype, ntype_maxd, itype, jtype, &
         tempk, press, dpwet, denswet, bckernel )






      implicit none


      integer, intent(in) :: nbin            
      integer, intent(in) :: nbin_maxd       
      integer, intent(in) :: ntype           
      integer, intent(in) :: ntype_maxd      
      integer, intent(in) :: itype, jtype    

      real(r8), intent(in)  :: tempk            
      real(r8), intent(in)  :: press            
      real(r8), intent(in)  :: dpwet(nbin_maxd,ntype_maxd)      
      real(r8), intent(in)  :: denswet(nbin_maxd,ntype_maxd)    
      real(r8), intent(out) :: bckernel(nbin,nbin)  


      integer i, j

      real(r8), parameter :: pi = 3.14159265358979323846_r8

      real(r8)   &
             apfreepath, avogad,   &
             boltz,   &
             cunning,   &
             deltasq_i, deltasq_j,   &
             diffus_i, diffus_j, diffus_ipj,   &
             dpwet_i, dpwet_j, dpwet_ipj,   &
             gasfreepathx2, gasspeed,   &
             knud,   &
             mass_i, mass_j, mwair,   &
             pi6,   &
             rgas, rhoair,   &
             speedsq_i, speedsq_j,   &
             tmp1,   &
             viscosd, viscosk










      pi6 = pi/6.0

      boltz = 1.38065e-16
      avogad = 6.02214e+23
      mwair = 28.966
      rgas = 8.31447e7

      rhoair = press*mwair/(rgas*tempk)

      viscosd = 1.8325e-04*(416.16/(tempk+120.0)) * (tempk/296.16)**1.5
      viscosk = viscosd/rhoair
      gasspeed = sqrt(8.0*boltz*tempk*avogad/(pi*mwair))
      gasfreepathx2 = 4.0*viscosk/gasspeed









       bckernel(:,:) = -2.0_r8

       do i = 1, nbin

          dpwet_i = dpwet(i,itype)                    
          mass_i = pi6*denswet(i,itype)*(dpwet_i**3)  
          knud = gasfreepathx2/dpwet_i  
          cunning = 1.0 + knud*(1.249 + 0.42*exp(-0.87/knud))
          diffus_i = boltz*tempk*cunning/(3.0*pi*dpwet_i*viscosd)
          speedsq_i = 8.0*boltz*tempk/(pi*mass_i)
          apfreepath = 8.0*diffus_i/(pi*sqrt(speedsq_i))
          tmp1 = (dpwet_i + apfreepath)**3   &
               - (dpwet_i*dpwet_i + apfreepath*apfreepath)**1.5
          deltasq_i = ( tmp1/(3.0*dpwet_i*apfreepath) - dpwet_i )**2

          do j = 1, nbin

           if ( (itype == jtype) .and.   &
                (bckernel(j,i) > -1.0_r8) ) then
             
             bckernel(i,j) = bckernel(j,i) 

           else

             dpwet_j = dpwet(j,jtype)                    
             mass_j = pi6*denswet(j,jtype)*(dpwet_j**3)  
             knud = gasfreepathx2/dpwet_j  
             cunning = 1.0 + knud*(1.249 + 0.42*exp(-0.87/knud))
             diffus_j = boltz*tempk*cunning/(3.0*pi*dpwet_j*viscosd)
             speedsq_j = 8.0*boltz*tempk/(pi*mass_j)
             apfreepath = 8.0*diffus_j/(pi*sqrt(speedsq_j))
             tmp1 = (dpwet_j + apfreepath)**3   &
                  - (dpwet_j*dpwet_j + apfreepath*apfreepath)**1.5
             deltasq_j = ( tmp1/(3.0*dpwet_j*apfreepath) - dpwet_j )**2

             dpwet_ipj = dpwet_i + dpwet_j
             diffus_ipj = diffus_i + diffus_j 
             tmp1 = (dpwet_ipj/(dpwet_ipj + 2.0*sqrt(deltasq_i + deltasq_j)))   &
                  + (8.0*diffus_ipj/(dpwet_ipj*sqrt(speedsq_i + speedsq_j)))
             bckernel(i,j) = max( 0.0_r8, 2.0*pi*dpwet_ipj*diffus_ipj/tmp1 )

           end if

          end do   
       end do   

       return
       end subroutine brownian_kernel_3dsub







	end module module_mosaic_coag3d
