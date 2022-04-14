	module module_mosaic_newnucb














	use module_data_mosaic_kind, only:  r8
	use module_data_mosaic_main, only:  pi 
	use module_peg_util



	implicit none
        
        
	real(r8), parameter, public :: adjust_factor_dnaitdt = 1.0_r8  
	real(r8), parameter, public :: adjust_factor_bin_tern_ratenucl = 1.0_r8  
	real(r8), parameter, public :: adjust_factor_pbl_ratenucl = 1.0_r8  
        

	contains




	subroutine mosaic_newnuc_1box(   &
		istat_newnuc, idiagbb_in,  &
		itstep, itype_newnuc, dtchem,   &
		temp_box, patm_box, cair_box, rh_box,   &
		dens_nh4so4a_in, mw_so4a, mw_nh4a, mw_h2o, mw_air,   &
		fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr,   &
		rbox0, rbox,   &
		adrydens_box, awetdens_box,   &
		adrydpav_box, awetdpav_box, adryqmas_box )














	use module_data_mosaic_main, only:  ntot_used, piover6, third
	use module_data_mosaic_boxmod, only:  kh2so4, knh3, kso2
	use module_data_mosaic_aero, only:  mnewnuc_flag1
	use module_data_mosaic_asecthp, only:   &
		ai_phase, dens_aer, dens_water_aer,   &
		dcen_sect, dlo_sect, dhi_sect,   &
		hyswptr_aer, lptr_nh4_aer, lptr_so4_aer, lunerr,   &
		massptr_aer, maxd_asize, maxd_atype,   &
		ncomp_aer, nsize_aer, numptr_aer, smallmassbb,   &
		volumcen_sect, volumlo_sect, volumhi_sect, waterptr_aer



	integer,  intent(inout) :: istat_newnuc   
	integer,  intent(in)    :: idiagbb_in     
	integer,  intent(in)    :: itstep         
	integer,  intent(in)    :: itype_newnuc   
	real(r8), intent(in)    :: dtchem         
	real(r8), intent(in)    :: temp_box       
	real(r8), intent(in)    :: patm_box       
	real(r8), intent(in)    :: cair_box       
	real(r8), intent(in)    :: rh_box         

	real(r8), intent(in) :: dens_nh4so4a_in   
	real(r8), intent(in) :: mw_so4a, mw_nh4a, mw_h2o, mw_air   

	real(r8), intent(in) :: fact_apmassmr, fact_apnumbmr, &
                                fact_apdens, fact_apdiam, fact_gasmr















	real(r8), intent(in)    :: rbox0(ntot_used)
	real(r8), intent(inout) :: rbox(ntot_used)



	real(r8), intent(inout) :: adrydens_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: awetdens_box(maxd_asize,maxd_atype)


	real(r8), intent(inout) :: adrydpav_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: awetdpav_box(maxd_asize,maxd_atype)

	real(r8), intent(inout) :: adryqmas_box(maxd_asize,maxd_atype)



	integer, parameter :: p1st = 1

	integer :: it=1, jt=1, k=1
	integer :: l, ll
	integer :: isize, itype, iphase
	integer :: iconform_numb
	integer :: idiagbb
	integer :: ldiagaa, lundiagaa
	integer :: nsize
	integer, save :: ncount(10)


	real(r8), parameter :: qh2so4_cutoff = 4.0e-16_r8


        real(r8), parameter :: a_zsr_xx1 =  1.15510
        real(r8), parameter :: a_zsr_xx2 = -3.20815
        real(r8), parameter :: a_zsr_xx3 =  2.71141
        real(r8), parameter :: a_zsr_xx4 =  2.01155
        real(r8), parameter :: a_zsr_xx5 = -4.71014
        real(r8), parameter :: a_zsr_xx6 =  2.04616
        real(r8), parameter :: b_zsr_xx  = 29.4779

	real(r8) :: aw
	real(r8) :: densdefault
	real(r8) :: dens_nh4so4a, dens_nh4so4a_kgm3
	real(r8) :: dtnuc
	real(r8) :: duma, dumb, dumc
	real(r8) :: fact_apvolu
	real(r8) :: h2so4_uptkrate
	real(r8) :: press_pa
	real(r8) :: qh2so4_avg, qh2so4_cur, qh2so4_del  
	real(r8) :: qnh3_avg, qnh3_cur, qnh3_del        
	real(r8) :: qso4a_del, qnh4a_del                
	real(r8) :: qnuma_del                           
	real(r8) :: tmpa, tmpb, tmpc, tmp_q2, tmp_q3
	real(r8) :: wwdens, wwdpav, wwmass, wwvolu
	real(r8) :: xxdens, xxdpav, xxmass, xxnumb, xxvolu

	real(r8),save :: tmpveca(10), tmpvecb(10), tmpvecc(10), tmpvecd(10), tmpvece(10), tmpvecf(10) 
	real(r8) :: dplom_nuc(maxd_asize), dphim_nuc(maxd_asize)
	real(r8) :: volumlo_nuc(maxd_asize), volumhi_nuc(maxd_asize)

	character(len=100) :: msg


        if ( kh2so4 < 1 .or. kh2so4 > ntot_used .or. &
             knh3   < 1 .or. knh3   > ntot_used .or. &
             kso2   < 1 .or. kso2   > ntot_used ) then
           write(msg,'(a,3i15)') &
              'module_mosaic_newnucb bad kh2so4/knh3/kso2', & 
              kh2so4, knh3, kso2
           call peg_message( lunerr, msg )
           call peg_error_fatal( lunerr, msg )
           call wrf_error_fatal3("<stdin>",171,&
trim(adjustl(msg)))
        end if


	istat_newnuc = 0
	if (mnewnuc_flag1 <= 0) then 
           
           return
	else if ( mnewnuc_flag1 ==  1 .or. &
             mnewnuc_flag1 ==  2 .or. &
             mnewnuc_flag1 ==  3 .or. &
             mnewnuc_flag1 == 11 .or. &
             mnewnuc_flag1 == 12 ) then
           continue
	else
           

		call peg_message( lunerr,   &
		'*** mosaic_newnuc_1box -- illegal mnewnuc_flag1' )
	    istat_newnuc = -1
	    return

	end if



	dtnuc = dtchem
	idiagbb = idiagbb_in
        lundiagaa = 93

	itype = itype_newnuc
	iphase = ai_phase
	nsize = nsize_aer(itype)
	fact_apvolu = fact_apdiam*fact_apdiam*fact_apdiam
	densdefault = dens_aer(1,itype)*fact_apdens
	volumlo_nuc(1:nsize) = volumlo_sect(1:nsize,itype)*fact_apvolu
	volumhi_nuc(1:nsize) = volumhi_sect(1:nsize,itype)*fact_apvolu
	dplom_nuc(1:nsize) = dlo_sect(1:nsize,itype)*0.01*fact_apdiam   
	dphim_nuc(1:nsize) = dhi_sect(1:nsize,itype)*0.01*fact_apdiam   










	    tmpveca(:) = 0.0         
	    tmpvecb(:) = +1.0e35     
	    tmpvecc(:) = -1.0e35     
	    tmpvecd(:) = 0.0         
	    tmpvece(:) = 0.0         
	    ncount(:) = 0



	ncount(1) = ncount(1) + 1

	qh2so4_cur = max(0.0_r8,rbox(kh2so4))*fact_gasmr

	if (qh2so4_cur <= qh2so4_cutoff) goto 2700

	qnh3_cur   = max(0.0_r8,rbox(knh3))*fact_gasmr
	qh2so4_avg = 0.5*( qh2so4_cur + max(0.0_r8,rbox0(kh2so4))*fact_gasmr )
	qnh3_avg   = 0.5*( qnh3_cur   + max(0.0_r8,rbox0(knh3))*fact_gasmr )




	tmp_q3 = qh2so4_cur

	tmp_q2 = max(0.0_r8,rbox0(kh2so4))*fact_gasmr


	if (tmp_q2 > tmp_q3) then
	   tmpc = tmp_q2 * exp( -20.0_r8 )
	   if (tmp_q3 <= tmpc) then
	      tmpb = 20.0_r8
	   else
	      tmpb = log( tmp_q2/tmp_q3 )
	   end if
	else
	   tmpb = 0.0
	end if
	h2so4_uptkrate = tmpb/dtchem


	qh2so4_del = 0.0
	qnh3_del = 0.0
	qnuma_del = 0.0
	qso4a_del = 0.0
	qnh4a_del = 0.0



	dens_nh4so4a = dens_nh4so4a_in*fact_apdens

	isize = 0


	if (mnewnuc_flag1 /= 3) then


            ldiagaa = idiagbb
            press_pa = patm_box*1.01325e5   
            dens_nh4so4a_kgm3 = dens_nh4so4a*1.0e3

            call mer07_veh02_nuc_mosaic_1box(   &
               mnewnuc_flag1, dtnuc, temp_box, rh_box, press_pa,   &
               qh2so4_cur, qh2so4_avg, qnh3_cur, qnh3_avg, h2so4_uptkrate,   &
               nsize, maxd_asize, dplom_nuc, dphim_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del,   &
               qh2so4_del, qnh3_del, dens_nh4so4a_kgm3, &
               itstep, ldiagaa, lundiagaa )








            dens_nh4so4a = dens_nh4so4a_kgm3/1.0e3

	else

            call wexler_nuc_mosaic_1box(   &
               dtnuc, temp_box, rh_box, cair_box,   &
               qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
               nsize, maxd_asize, volumlo_nuc, volumhi_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del,   &
               qh2so4_del, qnh3_del, dens_nh4so4a )

	
	end if



	tmpveca(1) = temp_box
	tmpveca(2) = rh_box
	tmpveca(3) = rbox(kso2)*fact_gasmr
	tmpveca(4) = qh2so4_avg
	tmpveca(5) = qnh3_avg
	tmpveca(6) = qnuma_del
	do l = 1, 6
	    tmpvecb(l) = min( tmpvecb(l), tmpveca(l) )
	    tmpvecc(l) = max( tmpvecc(l), tmpveca(l) )
	    tmpvecd(l) = tmpvecd(l) + tmpveca(l)
	    if (qnuma_del .gt. tmpvece(6)) tmpvece(l) = tmpveca(l)
	end do



	if (qnuma_del .le. 0.0) goto 2700


	if (isize .ne. 1) ncount(3) = ncount(3) + 1
	if ((isize .lt. 1) .or. (isize .gt. nsize)) then
	    write(msg,93010) 'newnucxx bad isize_nuc' , it, jt, k,   &
		isize, nsize
	    call peg_message( lunerr, msg )
	    goto 2700
	end if
93010	format( a, 3i3, 1p, 9e10.2 )


	ncount(2) = ncount(2) + 1


	rbox(kh2so4) = max( 0.0_r8, rbox(kh2so4) + qh2so4_del/fact_gasmr )
	rbox(knh3  ) = max( 0.0_r8, rbox(knh3  ) + qnh3_del/fact_gasmr )

	l = lptr_so4_aer(isize,itype,iphase)
	if (l .ge. p1st) then

	    tmpa = qso4a_del*mw_so4a/mw_air   
            tmpvecf(1) = rbox(l)
	    rbox(l) = rbox(l) + tmpa/fact_apmassmr
            
            tmpvecf(2) = rbox(l)
	    
	    tmpvecf(1:2) = tmpvecf(1:2) * fact_apmassmr  * 1.0e6              * (1.0e6_r8*cair_box*mw_air)
	    tmpvecf(3) = tmpvecf(2) - tmpvecf(1)
	    if (idiagbb >= 10) &
                 write(lundiagaa,'(a,1p,3e12.4)') 'so4 (ug/m3) o/n/d', tmpvecf(1:3)
            
	end if
	l = lptr_nh4_aer(isize,itype,iphase)
	if (l .ge. p1st) then

	    tmpa = qnh4a_del*mw_nh4a/mw_air   
            tmpvecf(1) = rbox(l)
	    rbox(l) = rbox(l) + tmpa/fact_apmassmr
            
            tmpvecf(2) = rbox(l)
	    
	    tmpvecf(1:2) = tmpvecf(1:2) * fact_apmassmr  * 1.0e6              * (1.0e6_r8*cair_box*mw_air)
	    tmpvecf(3) = tmpvecf(2) - tmpvecf(1)
	    if (idiagbb >= 10) &
	    write(lundiagaa,'(a,1p,3e12.4)') 'nh4 (ug/m3) o/n/d', tmpvecf(1:3)
            
	end if
	l = numptr_aer(isize,itype,iphase)

	tmpa = qnuma_del/mw_air   
        tmpvecf(1) = rbox(l)
	rbox(l) = rbox(l) + tmpa/fact_apnumbmr
        
        tmpvecf(2) = rbox(l)
	
	tmpvecf(1:2) = tmpvecf(1:2) * fact_apnumbmr  * (cair_box*mw_air)
	tmpvecf(3) = tmpvecf(2) - tmpvecf(1)
	if (idiagbb >= 10) &
	write(lundiagaa,'(a,1p,3e12.4)') 'num (#/cm3) o/n/d', tmpvecf(1:3)
        
	xxnumb = rbox(l)*fact_apnumbmr





	l = waterptr_aer(isize,itype)
	if ((rh_box .gt. 0.10) .and. (l .ge. p1st)) then
	    aw = min( rh_box, 0.98_r8 )
	    if (aw .lt. 0.97) then
		duma =       a_zsr_xx1 +   &
		        aw*( a_zsr_xx2 +   &
		        aw*( a_zsr_xx3 +   &
		        aw*( a_zsr_xx4 +   &
		        aw*( a_zsr_xx5 +   &
		        aw*  a_zsr_xx6 ))))
	    else
		dumb = -b_zsr_xx*log(aw)
	        dumb = max( dumb, 0.5_r8 )
		duma = 1.0/(1.0 + 55.509/dumb)
	    end if
	    duma = max( duma, 0.01_r8 )
	    dumc = (1.0 - duma)/duma

	    tmpa = qso4a_del*dumc*mw_h2o/mw_air
	    rbox(l) = rbox(l) + tmpa/fact_apmassmr
	end if






	xxmass = adryqmas_box(isize,itype)*fact_apmassmr
	xxdens = adrydens_box( isize,itype)*fact_apdens
	iconform_numb = 1

	if ((xxdens .lt. 0.1) .or. (xxdens .gt. 20.0)) then

	    continue
	else


	    xxvolu = xxmass/xxdens


	    duma = (qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mw_air   
	    xxmass = xxmass + duma
	    xxvolu  = xxvolu  + duma/dens_nh4so4a
	    if (xxmass .le. smallmassbb) then

		xxdens = 0.001
	    else if (xxmass .gt. 1000.0*xxvolu) then


		xxdens = 1000.0
	    else 
		xxdens = xxmass/xxvolu
	    end if
	end if

	if ((xxdens .lt. 0.1) .or. (xxdens .gt. 20.0)) then


	    ncount(4) = ncount(4) + 1
	    xxmass = 0.0
	    xxvolu  = 0.0
	    do ll = 1, ncomp_aer(itype)
		l = massptr_aer(ll,isize,itype,iphase)
		if (l .ge. p1st) then



		    duma = max( 0.0_r8, rbox(l) )*fact_apmassmr
		    xxmass = xxmass + duma
		    xxvolu = xxvolu + duma/(dens_aer(ll,itype)*fact_apdens)
		end if
	    end do
	end if

	if (xxmass .le. smallmassbb) then


	    ncount(5) = ncount(5) + 1
	    xxdens = densdefault
	    xxvolu = xxmass/xxdens
	    xxnumb = xxmass/(volumcen_sect(isize,itype)*fact_apvolu*xxdens)
	    xxdpav = dcen_sect(isize,itype)*fact_apdiam
	    wwdens = xxdens
	    wwdpav = xxdpav
	    iconform_numb = 0
	    l = waterptr_aer(isize,itype)
	    if (l .ge. p1st) rbox(l) = 0.0
	    l = hyswptr_aer(isize,itype)
	    if (l .ge. p1st) rbox(l) = 0.0
	else
	    xxdens = xxmass/xxvolu
	end if

	if (iconform_numb .gt. 0) then

	    if ( xxnumb .gt. xxvolu/volumlo_nuc(isize) ) then
		ncount(6) = ncount(6) + 1
		xxnumb = xxvolu/volumlo_nuc(isize)
		xxdpav = dlo_sect(isize,itype)*fact_apdiam
	    else if ( xxnumb .lt. xxvolu/volumhi_nuc(isize) ) then
		ncount(7) = ncount(7) + 1
		xxnumb = xxvolu/volumhi_nuc(isize)
		xxdpav = dhi_sect(isize,itype)*fact_apdiam
	    else
		xxdpav = (xxvolu/(xxnumb*piover6))**third
	    end if

	    tmpb = 0.0
	    l = waterptr_aer(isize,itype)

	    if (l .ge. p1st) tmpb = max(0.0_r8,rbox(l))*fact_apmassmr
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
	awetdens_box(isize,itype) = wwdens/fact_apdens
	awetdpav_box(isize,itype) = wwdpav/fact_apdiam


2700	continue


	if (idiagbb .ge. 100) then



            
	    write(lundiagaa,'(2a)')   &
		'newnucbb names  temp      rh        ',   &
		'so2       h2so4_avg nh3_avg   numa_del'
	    if (idiagbb .ge. 110) then
	      write(lundiagaa,93020) 'newnucbb mins ', tmpvecb(1:6)
	      write(lundiagaa,93020) 'newnucbb maxs ', tmpvecc(1:6)
	      duma = max( 1, ncount(1) ) 
	      write(lundiagaa,93020) 'newnucbb avgs ', tmpvecd(1:6)/duma
	      write(lundiagaa,93020) 'newnucbb hinuc', tmpvece(1:6)
	      write(lundiagaa,93020) 'newnucbb dtnuc', dtnuc
	    end if
	    write(lundiagaa,93030) 'newnucbb ncnt ', ncount(1:7)
	end if
93020	format( a, 1p, 10e10.2 )
93030	format( a, 1p, 10i10 )






	return
	end subroutine mosaic_newnuc_1box





        subroutine mer07_veh02_nuc_mosaic_1box(   &
           newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in,   &
           qh2so4_cur, qh2so4_avg, qnh3_cur, qnh3_avg, h2so4_uptkrate,   &
           nsize, maxd_asize, dplom_sect, dphim_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a, &
           itstep, ldiagaa, lundiagaa ) 











































      implicit none


        real(r8), intent(in) :: dtnuc             
        real(r8), intent(in) :: temp_in           
        real(r8), intent(in) :: rh_in             
        real(r8), intent(in) :: press_in          

        real(r8), intent(in) :: qh2so4_cur, qh2so4_avg
        real(r8), intent(in) :: qnh3_cur, qnh3_avg
             
             
             
        real(r8), intent(in) :: h2so4_uptkrate    

        integer,  intent(in) :: newnuc_method_flagaa    
                                                        
        integer,  intent(in) :: nsize                   
        integer,  intent(in) :: maxd_asize              
        real(r8), intent(in) :: dplom_sect(maxd_asize)  
        real(r8), intent(in) :: dphim_sect(maxd_asize)  
        integer,  intent(in) :: itstep 
        integer,  intent(in) :: ldiagaa, lundiagaa


        integer,  intent(out) :: isize_nuc        
        real(r8), intent(out) :: qnuma_del        
        real(r8), intent(out) :: qso4a_del        
        real(r8), intent(out) :: qnh4a_del        
        real(r8), intent(out) :: qh2so4_del       
        real(r8), intent(out) :: qnh3_del         
                                                  
        real(r8), intent(inout) :: dens_nh4so4a   





        real(r8) :: ratenuclt        
        real(r8) :: rateloge         
        real(r8) :: cnum_h2so4       
        real(r8) :: cnum_nh3         
        real(r8) :: cnum_tot         
        real(r8) :: radius_cluster   



        integer :: i
        integer :: igrow
        integer, save :: icase = 0, icase_reldiffmax = 0

        integer :: lun
        integer :: newnuc_method_flagaa2

        real(r8), parameter :: onethird = 1.0/3.0
        real(r8), parameter :: avogad = 6.022e23   
        real(r8), parameter :: mw_air = 28.966     

        real(r8), parameter :: accom_coef_h2so4 = 0.65   







        real(r8), parameter :: dens_ammsulf   = 1.770e3
        real(r8), parameter :: dens_ammbisulf = 1.770e3
        real(r8), parameter :: dens_sulfacid  = 1.770e3
        real(r8), parameter :: dens_water     = 1.0e3




        real(r8), parameter :: mw_ammsulf   = 132.0
        real(r8), parameter :: mw_ammbisulf = 114.0
        real(r8), parameter :: mw_sulfacid  =  96.0

        real(r8), parameter :: mw_so4a      =  96.0
        real(r8), parameter :: mw_nh4a      =  18.0
        real(r8), parameter :: mw_water     =  18.0

        real(r8), save :: reldiffmax = 0.0

        real(r8) cair                     
        real(r8) cs_prime_kk              
        real(r8) cs_kk                    
        real(r8) dens_part                
        real(r8) dfin_kk, dnuc_kk         
        real(r8) dpdry_clus               
        real(r8) dpdry_part               
        real(r8) tmpa, tmpb, tmpc, tmpe, tmpq
        real(r8) tmpa1, tmpb1
        real(r8) tmp_m1, tmp_m2, tmp_m3, tmp_n1, tmp_n2, tmp_n3
        real(r8) tmp_spd                  
        real(r8) factor_kk
        real(r8) fogas, foso4a, fonh4a, fonuma
        real(r8) freduce                  
                                          
        real(r8) freducea, freduceb
        real(r8) gamma_kk                 
        real(r8) gr_kk                    
        real(r8) kgaero_per_moleso4a      
        real(r8) mass_part                
        real(r8) molenh4a_per_moleso4a    
        real(r8) nh3ppt, nh3ppt_bb        
        real(r8) nu_kk                    
        real(r8) qmolnh4a_del_max         
        real(r8) qmolso4a_del_max         
        real(r8) ratenuclt_bb             
        real(r8) ratenuclt_kk             
        real(r8) rh_bb                    
        real(r8) so4vol_in                
        real(r8) so4vol_bb                
        real(r8) temp_bb                  
        real(r8) voldry_clus              
        real(r8) voldry_part              
        real(r8) wetvol_dryvol            
        real(r8) wet_volfrac_so4a         







        isize_nuc = 1
        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0









        cair = press_in/(temp_in*8.3144)
        so4vol_in = qh2so4_avg * cair * avogad * 1.0e-6
        nh3ppt    = qnh3_avg * 1.0e12
        newnuc_method_flagaa2 = 0

        if ((newnuc_method_flagaa == 1) .and. (nh3ppt >= 0.1)) then
            if (so4vol_in < 5.0e4) return


            temp_bb = max( 235.0_r8, min( 295.0_r8, temp_in ) )
            rh_bb = max( 0.05_r8, min( 0.95_r8, rh_in ) )
            so4vol_bb = max( 5.0e4_r8, min( 1.0e9_r8, so4vol_in ) )
            nh3ppt_bb = max( 0.1_r8, min( 1.0e3_r8, nh3ppt ) )
            call ternary_nuc_merik2007(   &
                temp_bb, rh_bb, so4vol_bb, nh3ppt_bb,   &
                rateloge,   &
                cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )
            newnuc_method_flagaa2 = 1

        else if ((newnuc_method_flagaa == 1) .or. &
                 (newnuc_method_flagaa == 2)) then

            if (so4vol_in < 1.0e4) return

            temp_bb = max( 230.15_r8, min( 305.15_r8, temp_in ) )
            rh_bb = max( 1.0e-4_r8, min( 1.0_r8, rh_in ) )
            so4vol_bb = max( 1.0e4_r8, min( 1.0e11_r8, so4vol_in ) )
            call binary_nuc_vehk2002(   &
                temp_bb, rh_bb, so4vol_bb,   &
                ratenuclt, rateloge,   &
                cnum_h2so4, cnum_tot, radius_cluster )
            cnum_nh3 = 0.0
            newnuc_method_flagaa2 = 2
            
         else if ((newnuc_method_flagaa == 11) .or. &
              (newnuc_method_flagaa == 12)) then
            call pbl_nuc_wang2008( so4vol_in,   &
                 newnuc_method_flagaa, newnuc_method_flagaa2,   &
                 ratenuclt, rateloge,   &
                 cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )
            
        else
            return
        end if
        if (ldiagaa > 0) &
        write(lundiagaa,'(a,i8,1p,e10.2)') 'newnuc - itstep, ratenucl =', itstep, exp( rateloge )


        if (rateloge  .le. -13.82_r8) return

        ratenuclt = exp( rateloge )
        ratenuclt_bb = ratenuclt*1.0e6_r8



        tmpa = max( 0.10_r8, min( 0.95_r8, rh_in ) )
        wetvol_dryvol = 1.0 - 0.56/log(tmpa)




        voldry_clus = ( max(cnum_h2so4,1.0_r8)*mw_so4a + cnum_nh3*mw_nh4a ) /   &
                      (1.0e3*dens_sulfacid*avogad)
        dpdry_clus = (voldry_clus*6.0/pi)**onethird

        isize_nuc = 1
        dpdry_part = dplom_sect(1)
        if (dpdry_clus <= dplom_sect(1)) then
           igrow = 1   
        else if (dpdry_clus >= dphim_sect(nsize)) then
           igrow = 0
           isize_nuc = nsize
           dpdry_part = dphim_sect(nsize)
        else
           igrow = 0
           do i = 1, nsize
              if (dpdry_clus < dphim_sect(i)) then
                 isize_nuc = i
                 dpdry_part = dpdry_clus
                 dpdry_part = min( dpdry_part, dphim_sect(i) )
                 dpdry_part = max( dpdry_part, dplom_sect(i) )
                 exit
              end if
           end do
        end if
        voldry_part = (pi/6.0)*(dpdry_part**3)
        if (ldiagaa > 0) &
             write(lundiagaa,'(a,i8,i4)') 'newnuc - itstep, igrow    =', itstep, igrow










        if (igrow .le. 0) then

           tmp_n1 = 0.0
           tmp_n2 = 0.0
           tmp_n3 = 1.0
        else if (qnh3_cur .ge. qh2so4_cur) then


           tmp_n1 = (qnh3_cur/qh2so4_cur) - 1.0
           tmp_n1 = max( 0.0_r8, min( 1.0_r8, tmp_n1 ) )
           tmp_n2 = 1.0 - tmp_n1
           tmp_n3 = 0.0
        else


           tmp_n1 = 0.0
           tmp_n2 = (qnh3_cur/qh2so4_cur)
           tmp_n2 = max( 0.0_r8, min( 1.0_r8, tmp_n2 ) )
           tmp_n3 = 1.0 - tmp_n2
	end if

        tmp_m1 = tmp_n1*mw_ammsulf
        tmp_m2 = tmp_n2*mw_ammbisulf
        tmp_m3 = tmp_n3*mw_sulfacid
        dens_part = (tmp_m1 + tmp_m2 + tmp_m3)/   &
           ((tmp_m1/dens_ammsulf) + (tmp_m2/dens_ammbisulf)   &
                                  + (tmp_m3/dens_sulfacid))

	if (abs(dens_nh4so4a-1750.0) .le. 250.0) then
	    dens_part = dens_nh4so4a
	else
            dens_nh4so4a = dens_part
	end if

        mass_part  = voldry_part*dens_part 

        molenh4a_per_moleso4a = 2.0*tmp_n1 + tmp_n2  

        kgaero_per_moleso4a = 1.0e-3*(tmp_m1 + tmp_m2 + tmp_m3)  


        tmpb = 1.0 + molenh4a_per_moleso4a*17.0/98.0
        wet_volfrac_so4a = 1.0 / ( wetvol_dryvol * tmpb )





        if (igrow <=  0) then
            factor_kk = 1.0

        else


            tmp_spd = 14.7*sqrt(temp_in)   
            gr_kk = 3.0e-9*tmp_spd*mw_sulfacid*so4vol_in/   &
                    (dens_part*wet_volfrac_so4a)





            dfin_kk = 1.0e9 * dpdry_part * (wetvol_dryvol**onethird)

            dnuc_kk = 2.0*radius_cluster
            dnuc_kk = max( dnuc_kk, 1.0_r8 )


            gamma_kk = 0.23 * (dnuc_kk)**0.2   &
                     * (dfin_kk/3.0)**0.075   &
                     * (dens_part*1.0e-3)**(-0.33)   &
                     * (temp_in/293.0)**(-0.75)











            tmpa = h2so4_uptkrate * 3600.0
            tmpa1 = tmpa
            tmpa = max( tmpa, 0.0_r8 )

            tmpb = 6.7037e-6 * (temp_in**0.75) / cair
            tmpb1 = tmpb         
            tmpb = tmpb*3600.0   
            cs_prime_kk = tmpa/(4.0*pi*tmpb*accom_coef_h2so4)
            cs_kk = cs_prime_kk*4.0*pi*tmpb1


            nu_kk = gamma_kk*cs_prime_kk/gr_kk

            factor_kk = exp( (nu_kk/dfin_kk) - (nu_kk/dnuc_kk) )

        end if
        ratenuclt_kk = ratenuclt_bb*factor_kk



        tmpa = max( 0.0_r8, (ratenuclt_kk*dtnuc*mass_part) )

        tmpe = tmpa/(kgaero_per_moleso4a*cair)


        qmolso4a_del_max = tmpe


        freducea = 1.0
        if (qmolso4a_del_max .gt. qh2so4_cur) then
           freducea = qh2so4_cur/qmolso4a_del_max
        end if


        freduceb = 1.0
        if (molenh4a_per_moleso4a .ge. 1.0e-10) then

           qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
           if (qmolnh4a_del_max .gt. qnh3_cur) then
              freduceb = qnh3_cur/qmolnh4a_del_max
           end if
        end if
        freduce = min( freducea, freduceb )



        if (freduce*ratenuclt_kk .le. 1.0e-12) return















        tmpa = 0.9999
        qh2so4_del = min( tmpa*qh2so4_cur, freduce*qmolso4a_del_max )
        qnh3_del   = min( tmpa*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del


        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del

        qnuma_del = 1.0e-3*(qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part



        tmpa = max( 0.0_r8, (ratenuclt_kk*dtnuc/cair) )

        tmpb = tmpa*freduce

        tmpc = (tmpb - qnuma_del)/max(tmpb, qnuma_del, 1.0e-35_r8)






        if (ldiagaa <= 99) return
        lun = lundiagaa

        icase = icase + 1
        if (abs(tmpc) .gt. abs(reldiffmax)) then
           reldiffmax = tmpc
           icase_reldiffmax = icase
        end if



           write(lun,'(a,2i9,1p,e10.2)')   &
               'vehkam bin-nuc icase, icase_rdmax =',   &
               icase, icase_reldiffmax, reldiffmax
           if (freduceb .lt. freducea) then
              if (abs(freducea-freduceb) .gt.   &
                   3.0e-7*max(freduceb,freducea)) write(lun,'(a,1p,2e15.7)')   &
                 'freducea, b =', freducea, freduceb
           end if







        fogas  = 1.0
        foso4a = 1.0
        fonh4a = 1.0
        fonuma = 1.0




        write(lun,'(a,2i5)') 'newnuc_method_flagaa/aa2',   &
           newnuc_method_flagaa, newnuc_method_flagaa2

        write(lun,9210)
        write(lun,9201) temp_in, rh_in,   &
           ratenuclt, 2.0*radius_cluster*1.0e-7, dpdry_part*1.0e2,   &
           voldry_part*1.0e6, float(igrow)
        write(lun,9215)
        write(lun,9201)   &
           qh2so4_avg*fogas, qnh3_avg*fogas,  &
           qh2so4_cur*fogas, qnh3_cur*fogas,  &
           qh2so4_del*fogas, qnh3_del*fogas,  &
           qso4a_del*foso4a, qnh4a_del*fonh4a

        write(lun,9220)
        write(lun,9201)   &
           dtnuc, dens_nh4so4a*1.0e-3,   &
           (qnh3_cur/qh2so4_cur), molenh4a_per_moleso4a,   &
           qnuma_del*fonuma, tmpb*fonuma, tmpc, freduce





        write(lun,9230)
        write(lun,9201)   &
           press_in, cair*1.0e-6, so4vol_in,   &
           wet_volfrac_so4a, wetvol_dryvol, dens_part*1.0e-3

        if (igrow > 0) then
        write(lun,9240)
        write(lun,9201)   &
           tmp_spd, gr_kk, dnuc_kk, dfin_kk,   &
           gamma_kk, tmpa1, tmpb1, cs_kk

        write(lun,9250)
        write(lun,9201)   &
           cs_prime_kk, nu_kk, factor_kk, ratenuclt,   &
           ratenuclt_kk*1.0e-6
        end if

9201    format ( 1p, 40e10.2  )
9210    format (   &
        '      temp        rh',   &
        '   ratenuc  dia_clus ddry_part',   &
        ' vdry_part     igrow' )
9215    format (   &
        '  h2so4avg  h2so4pre',   &
        '  h2so4cur   nh3_cur',   &
        '  h2so4del   nh3_del',   &
        '  so4a_del  nh4a_del' )
9220    format (    &
        '     dtnuc    dens_a   nh/so g   nh/so a',   &
        '  numa_del  numa_dl2   reldiff   freduce' )
9230    format (   &
        '  press_in      cair so4_volin',   &
        ' wet_volfr wetv_dryv dens_part' )
9240    format (   &
        '   tmp_spd     gr_kk   dnuc_kk   dfin_kk',   &
        '  gamma_kk     tmpa1     tmpb1     cs_kk' )
9250    format (   &
        ' cs_pri_kk     nu_kk factor_kk ratenuclt',   &
        ' ratenu_kk' )


        return
        end subroutine mer07_veh02_nuc_mosaic_1box






        subroutine pbl_nuc_wang2008( so4vol,   &
            newnuc_method_flagaa, newnuc_method_flagaa2,   &
            ratenucl, rateloge,   &
            cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )







        implicit none


        real(r8), intent(in) :: so4vol            
        integer, intent(in)  :: newnuc_method_flagaa  
                                


        integer, intent(inout)  :: newnuc_method_flagaa2
        real(r8), intent(inout) :: ratenucl         
        real(r8), intent(inout) :: rateloge         

        real(r8), intent(inout) :: cnum_tot         
                                                    
        real(r8), intent(inout) :: cnum_h2so4       
        real(r8), intent(inout) :: cnum_nh3         
        real(r8), intent(inout) :: radius_cluster   



        real(r8) :: tmp_diam, tmp_mass, tmp_volu
        real(r8) :: tmp_rateloge, tmp_ratenucl





        if (newnuc_method_flagaa == 11) then
           tmp_ratenucl = 1.0e-6_r8 * so4vol
        else if (newnuc_method_flagaa == 12) then
           tmp_ratenucl = 1.0e-12_r8 * (so4vol**2)
        else
           return
        end if
        tmp_ratenucl = tmp_ratenucl * adjust_factor_pbl_ratenucl
        tmp_rateloge = log( tmp_ratenucl )




        rateloge = tmp_rateloge
        ratenucl = tmp_ratenucl
        newnuc_method_flagaa2 = newnuc_method_flagaa



        radius_cluster = 0.5_r8




        tmp_diam = radius_cluster * 2.0e-7_r8   
        tmp_volu = (tmp_diam**3) * (pi/6.0_r8)  
        tmp_mass = tmp_volu * 1.8_r8            
        cnum_h2so4 = (tmp_mass / 98.0_r8) * 6.023e23_r8   
        cnum_tot = cnum_h2so4
        cnum_nh3 = 0.0_r8


        return
        end subroutine pbl_nuc_wang2008






        subroutine binary_nuc_vehk2002( temp, rh, so4vol,   &
            ratenucl, rateloge,   &
            cnum_h2so4, cnum_tot, radius_cluster )









        implicit none


        real(r8), intent(in) :: temp              
        real(r8), intent(in) :: rh                
        real(r8), intent(in) :: so4vol            


        real(r8), intent(out) :: ratenucl         
        real(r8), intent(out) :: rateloge         

        real(r8), intent(out) :: cnum_h2so4       
                                                  
        real(r8), intent(out) :: cnum_tot         
                                                  
        real(r8), intent(out) :: radius_cluster   



        real(r8) :: crit_x
        real(r8) :: acoe, bcoe, ccoe, dcoe, ecoe, fcoe, gcoe, hcoe, icoe, jcoe
        real(r8) :: tmpa, tmpb





        crit_x = 0.740997 - 0.00266379 * temp   &
               - 0.00349998 * log (so4vol)   &
               + 0.0000504022 * temp * log (so4vol)   &
               + 0.00201048 * log (rh)   &
               - 0.000183289 * temp * log (rh)   &
               + 0.00157407 * (log (rh)) ** 2.0   &
               - 0.0000179059 * temp * (log (rh)) ** 2.0   &
               + 0.000184403 * (log (rh)) ** 3.0   &
               - 1.50345e-6 * temp * (log (rh)) ** 3.0



        acoe    = 0.14309+2.21956*temp   &
                - 0.0273911 * temp**2.0   &
                + 0.0000722811 * temp**3.0 + 5.91822/crit_x

        bcoe    = 0.117489 + 0.462532 *temp   &
                - 0.0118059 * temp**2.0   &
                + 0.0000404196 * temp**3.0 + 15.7963/crit_x

        ccoe    = -0.215554-0.0810269 * temp   &
                + 0.00143581 * temp**2.0   &
                - 4.7758e-6 * temp**3.0   &
                - 2.91297/crit_x

        dcoe    = -3.58856+0.049508 * temp   &
                - 0.00021382 * temp**2.0   &
                + 3.10801e-7 * temp**3.0   &
                - 0.0293333/crit_x

        ecoe    = 1.14598 - 0.600796 * temp   &
                + 0.00864245 * temp**2.0   &
                - 0.0000228947 * temp**3.0   &
                - 8.44985/crit_x

        fcoe    = 2.15855 + 0.0808121 * temp   &
                -0.000407382 * temp**2.0   &
                -4.01957e-7 * temp**3.0   &
                + 0.721326/crit_x

        gcoe    = 1.6241 - 0.0160106 * temp   &
                + 0.0000377124 * temp**2.0   &
                + 3.21794e-8 * temp**3.0   &
                - 0.0113255/crit_x

        hcoe    = 9.71682 - 0.115048 * temp   &
                + 0.000157098 * temp**2.0   &
                + 4.00914e-7 * temp**3.0   &
                + 0.71186/crit_x

        icoe    = -1.05611 + 0.00903378 * temp   &
                - 0.0000198417 * temp**2.0   &
                + 2.46048e-8  * temp**3.0   &
                - 0.0579087/crit_x

        jcoe    = -0.148712 + 0.00283508 * temp   &
                - 9.24619e-6  * temp**2.0   &
                + 5.00427e-9 * temp**3.0   &
                - 0.0127081/crit_x

        tmpa     =     (   &
                  acoe   &
                + bcoe * log (rh)   &
                + ccoe * ( log (rh))**2.0   &
                + dcoe * ( log (rh))**3.0   &
                + ecoe * log (so4vol)   &
                + fcoe * (log (rh)) * (log (so4vol))   &
                + gcoe * ((log (rh) ) **2.0)   &
                       * (log (so4vol))   &
                + hcoe * (log (so4vol)) **2.0   &
                + icoe * log (rh)   &
                       * ((log (so4vol)) **2.0)   &
                + jcoe * (log (so4vol)) **3.0   &
                )
        rateloge = tmpa
        tmpa = min( tmpa, log(1.0e38_r8) )
        ratenucl = exp ( tmpa )





        acoe    = -0.00295413 - 0.0976834*temp   &
                + 0.00102485 * temp**2.0   &
                - 2.18646e-6 * temp**3.0 - 0.101717/crit_x

        bcoe    = -0.00205064 - 0.00758504*temp   &
                + 0.000192654 * temp**2.0   &
                - 6.7043e-7 * temp**3.0 - 0.255774/crit_x

        ccoe    = +0.00322308 + 0.000852637 * temp   &
                - 0.0000154757 * temp**2.0   &
                + 5.66661e-8 * temp**3.0   &
                + 0.0338444/crit_x

        dcoe    = +0.0474323 - 0.000625104 * temp   &
                + 2.65066e-6 * temp**2.0   &
                - 3.67471e-9 * temp**3.0   &
                - 0.000267251/crit_x

        ecoe    = -0.0125211 + 0.00580655 * temp   &
                - 0.000101674 * temp**2.0   &
                + 2.88195e-7 * temp**3.0   &
                + 0.0942243/crit_x

        fcoe    = -0.038546 - 0.000672316 * temp   &
                + 2.60288e-6 * temp**2.0   &
                + 1.19416e-8 * temp**3.0   &
                - 0.00851515/crit_x

        gcoe    = -0.0183749 + 0.000172072 * temp   &
                - 3.71766e-7 * temp**2.0   &
                - 5.14875e-10 * temp**3.0   &
                + 0.00026866/crit_x

        hcoe    = -0.0619974 + 0.000906958 * temp   &
                - 9.11728e-7 * temp**2.0   &
                - 5.36796e-9 * temp**3.0   &
                - 0.00774234/crit_x

        icoe    = +0.0121827 - 0.00010665 * temp   &
                + 2.5346e-7 * temp**2.0   &
                - 3.63519e-10 * temp**3.0   &
                + 0.000610065/crit_x

        jcoe    = +0.000320184 - 0.0000174762 * temp   &
                + 6.06504e-8 * temp**2.0   &
                - 1.4177e-11 * temp**3.0   &
                + 0.000135751/crit_x

        cnum_tot = exp (   &
                  acoe   &
                + bcoe * log (rh)   &
                + ccoe * ( log (rh))**2.0   &
                + dcoe * ( log (rh))**3.0   &
                + ecoe * log (so4vol)   &
                + fcoe * (log (rh)) * (log (so4vol))   &
                + gcoe * ((log (rh) ) **2.0)   &
                       * (log (so4vol))   &
                + hcoe * (log (so4vol)) **2.0   &
                + icoe * log (rh)   &
                       * ((log (so4vol)) **2.0)   &
                + jcoe * (log (so4vol)) **3.0   &
                )

        cnum_h2so4 = cnum_tot * crit_x


        radius_cluster = exp( -1.6524245 + 0.42316402*crit_x   &
                              + 0.3346648*log(cnum_tot) )
      

      return
      end subroutine binary_nuc_vehk2002





subroutine ternary_nuc_merik2007( t, rh, c2, c3, j_log, ntot, nacid, namm, r )























implicit none

real(r8), intent(in) :: t, rh, c2, c3
real(r8), intent(out) :: j_log, ntot, nacid, namm, r
real(r8) :: j, t_onset

t_onset=143.6002929064716 + 1.0178856665693992*rh + &
   10.196398812974294*log(c2) - &
   0.1849879416839113*log(c2)**2 - 17.161783213150173*log(c3) + &
   (109.92469248546053*log(c3))/log(c2) + &
   0.7734119613144357*log(c2)*log(c3) - 0.15576469879527022*log(c3)**2

if(t_onset.gt.t) then 

   j_log=-12.861848898625231 + 4.905527742256349*c3 - 358.2337705052991*rh -& 
   0.05463019231872484*c3*t + 4.8630382337426985*rh*t + &
   0.00020258394697064567*c3*t**2 - 0.02175548069741675*rh*t**2 - &
   2.502406532869512e-7*c3*t**3 + 0.00003212869941055865*rh*t**3 - &
   4.39129415725234e6/log(c2)**2 + (56383.93843154586*t)/log(c2)**2 -& 
   (239.835990963361*t**2)/log(c2)**2 + &
   (0.33765136625580167*t**3)/log(c2)**2 - &
   (629.7882041830943*rh)/(c3**3*log(c2)) + &
   (7.772806552631709*rh*t)/(c3**3*log(c2)) - &
   (0.031974053936299256*rh*t**2)/(c3**3*log(c2)) + &
   (0.00004383764128775082*rh*t**3)/(c3**3*log(c2)) + &
   1200.472096232311*log(c2) - 17.37107890065621*t*log(c2) + &
   0.08170681335921742*t**2*log(c2) - 0.00012534476159729881*t**3*log(c2) - &
   14.833042158178936*log(c2)**2 + 0.2932631303555295*t*log(c2)**2 - &
   0.0016497524241142845*t**2*log(c2)**2 + &
   2.844074805239367e-6*t**3*log(c2)**2 - 231375.56676032578*log(c3) - &
   100.21645273730675*rh*log(c3) + 2919.2852552424706*t*log(c3) + &
   0.977886555834732*rh*t*log(c3) - 12.286497122264588*t**2*log(c3) - &
   0.0030511783284506377*rh*t**2*log(c3) + &
   0.017249301826661612*t**3*log(c3) + 2.967320346100855e-6*rh*t**3*log(c3) + &
   (2.360931724951942e6*log(c3))/log(c2) - &
   (29752.130254319443*t*log(c3))/log(c2) + &
   (125.04965118142027*t**2*log(c3))/log(c2) - &
   (0.1752996881934318*t**3*log(c3))/log(c2) + &
   5599.912337254629*log(c2)*log(c3) - 70.70896612937771*t*log(c2)*log(c3) + &
   0.2978801613269466*t**2*log(c2)*log(c3) - &
   0.00041866525019504*t**3*log(c2)*log(c3) + 75061.15281456841*log(c3)**2 - &
   931.8802278173565*t*log(c3)**2 + 3.863266220840964*t**2*log(c3)**2 - &
   0.005349472062284983*t**3*log(c3)**2 - &
   (732006.8180571689*log(c3)**2)/log(c2) + &
   (9100.06398573816*t*log(c3)**2)/log(c2) - &
   (37.771091915932004*t**2*log(c3)**2)/log(c2) + &
   (0.05235455395566905*t**3*log(c3)**2)/log(c2) - &
   1911.0303773001353*log(c2)*log(c3)**2 + &
   23.6903969622286*t*log(c2)*log(c3)**2 - &
   0.09807872005428583*t**2*log(c2)*log(c3)**2 + &
   0.00013564560238552576*t**3*log(c2)*log(c3)**2 - &
   3180.5610833308*log(c3)**3 + 39.08268568672095*t*log(c3)**3 - &
   0.16048521066690752*t**2*log(c3)**3 + &
   0.00022031380023793877*t**3*log(c3)**3 + &
   (40751.075322248245*log(c3)**3)/log(c2) - &
   (501.66977622013934*t*log(c3)**3)/log(c2) + &
   (2.063469732254135*t**2*log(c3)**3)/log(c2) - &
   (0.002836873785758324*t**3*log(c3)**3)/log(c2) + &
   2.792313345723013*log(c2)**2*log(c3)**3 - &
   0.03422552111802899*t*log(c2)**2*log(c3)**3 + &
   0.00014019195277521142*t**2*log(c2)**2*log(c3)**3 - &
   1.9201227328396297e-7*t**3*log(c2)**2*log(c3)**3 - &
   980.923146020468*log(rh) + 10.054155220444462*t*log(rh) - &
   0.03306644502023841*t**2*log(rh) + 0.000034274041225891804*t**3*log(rh) + &
   (16597.75554295064*log(rh))/log(c2) - &
   (175.2365504237746*t*log(rh))/log(c2) + &
   (0.6033215603167458*t**2*log(rh))/log(c2) - &
   (0.0006731787599587544*t**3*log(rh))/log(c2) - &
   89.38961120336789*log(c3)*log(rh) + 1.153344219304926*t*log(c3)*log(rh) - &
   0.004954549700267233*t**2*log(c3)*log(rh) + &
   7.096309866238719e-6*t**3*log(c3)*log(rh) + &
   3.1712136610383244*log(c3)**3*log(rh) - &
   0.037822330602328806*t*log(c3)**3*log(rh) + &
   0.0001500555743561457*t**2*log(c3)**3*log(rh) - &
   1.9828365865570703e-7*t**3*log(c3)**3*log(rh)

   j=exp(j_log)

   ntot=57.40091052369212 - 0.2996341884645408*t + &
   0.0007395477768531926*t**2 - &
   5.090604835032423*log(c2) + 0.011016634044531128*t*log(c2) + &
   0.06750032251225707*log(c2)**2 - 0.8102831333223962*log(c3) + &
   0.015905081275952426*t*log(c3) - 0.2044174683159531*log(c2)*log(c3) + &
   0.08918159167625832*log(c3)**2 - 0.0004969033586666147*t*log(c3)**2 + &
   0.005704394549007816*log(c3)**3 + 3.4098703903474368*log(j) - &
   0.014916956508210809*t*log(j) + 0.08459090011666293*log(c3)*log(j) - &
   0.00014800625143907616*t*log(c3)*log(j) + 0.00503804694656905*log(j)**2
 
   r=3.2888553966535506e-10 - 3.374171768439839e-12*t + &
   1.8347359507774313e-14*t**2 + 2.5419844298881856e-12*log(c2) - &
   9.498107643050827e-14*t*log(c2) + 7.446266520834559e-13*log(c2)**2 + &
   2.4303397746137294e-11*log(c3) + 1.589324325956633e-14*t*log(c3) - &
   2.034596219775266e-12*log(c2)*log(c3) - 5.59303954457172e-13*log(c3)**2 - &
   4.889507104645867e-16*t*log(c3)**2 + 1.3847024107506764e-13*log(c3)**3 + &
   4.141077193427042e-15*log(j) - 2.6813110884009767e-14*t*log(j) + &
   1.2879071621313094e-12*log(c3)*log(j) - &
   3.80352446061867e-15*t*log(c3)*log(j) - 1.8790172502456827e-14*log(j)**2
 
   nacid=-4.7154180661803595 + 0.13436423483953885*t - & 
   0.00047184686478816176*t**2 - & 
   2.564010713640308*log(c2) + 0.011353312899114723*t*log(c2) + &
   0.0010801941974317014*log(c2)**2 + 0.5171368624197119*log(c3) - &
   0.0027882479896204665*t*log(c3) + 0.8066971907026886*log(c3)**2 - & 
   0.0031849094214409335*t*log(c3)**2 - 0.09951184152927882*log(c3)**3 + &
   0.00040072788891745513*t*log(c3)**3 + 1.3276469271073974*log(j) - &
   0.006167654171986281*t*log(j) - 0.11061390967822708*log(c3)*log(j) + &
   0.0004367575329273496*t*log(c3)*log(j) + 0.000916366357266258*log(j)**2
 
   namm=71.20073903979772 - 0.8409600103431923*t + &
   0.0024803006590334922*t**2 + &
   2.7798606841602607*log(c2) - 0.01475023348171676*t*log(c2) + &
   0.012264508212031405*log(c2)**2 - 2.009926050440182*log(c3) + &
   0.008689123511431527*t*log(c3) - 0.009141180198955415*log(c2)*log(c3) + &
   0.1374122553905617*log(c3)**2 - 0.0006253227821679215*t*log(c3)**2 + &
   0.00009377332742098946*log(c3)**3 + 0.5202974341687757*log(j) - &
   0.002419872323052805*t*log(j) + 0.07916392322884074*log(c3)*log(j) - &
   0.0003021586030317366*t*log(c3)*log(j) + 0.0046977006608603395*log(j)**2

else

   j_log=-300.
end if

return

end  subroutine ternary_nuc_merik2007





        subroutine wexler_nuc_mosaic_1box(   &
           dtnuc, temp_in, rh_in, cair,   &
           qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
           nsize, maxd_asize, volumlo_sect, volumhi_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a )






















        implicit none


        real(r8), intent(in) :: dtnuc             
        real(r8), intent(in) :: temp_in           
        real(r8), intent(in) :: rh_in             
        real(r8), intent(in) :: cair              

        real(r8), intent(in) :: qh2so4_avg, qh2so4_cur   
        real(r8), intent(in) :: qnh3_avg, qnh3_cur       
             
             

        integer, intent(in) :: nsize                    
        integer, intent(in) :: maxd_asize               
        real(r8), intent(in) :: volumlo_sect(maxd_asize)    
        real(r8), intent(in) :: volumhi_sect(maxd_asize)    


        integer, intent(out) :: isize_nuc     
        real(r8), intent(out) :: qnuma_del        
        real(r8), intent(out) :: qso4a_del        
        real(r8), intent(out) :: qnh4a_del        
        real(r8), intent(out) :: qh2so4_del       
        real(r8), intent(out) :: qnh3_del         
                                              


        real(r8), intent(inout) :: dens_nh4so4a   



        integer i
        integer, save :: icase = 0, icase_reldiffmax = 0

        real(r8), parameter :: pi = 3.1415926536
        real(r8), parameter :: avogad = 6.022e23   
        real(r8), parameter :: mw_air = 28.966     



        real(r8), parameter :: dens_ammsulf   = 1.769
        real(r8), parameter :: dens_ammbisulf = 1.78
        real(r8), parameter :: dens_sulfacid  = 1.841




        real(r8), parameter :: mw_ammsulf   = 132.0
        real(r8), parameter :: mw_ammbisulf = 114.0
        real(r8), parameter :: mw_sulfacid  =  96.0

        real(r8), parameter :: mw_so4a      =  96.0
        real(r8), parameter :: mw_nh4a      =  18.0

        real(r8), save :: reldiffmax = 0.0

        real(r8) ch2so4_crit              
        real(r8) dens_part                
        real(r8) duma, dumb, dumc, dume
        real(r8) dum_m1, dum_m2, dum_m3, dum_n1, dum_n2, dum_n3
        real(r8) fogas, foso4a, fonh4a, fonuma
        real(r8) mass_part                
        real(r8) molenh4a_per_moleso4a    
        real(r8) qh2so4_crit              
        real(r8) qh2so4_avail             
        real(r8) vol_part                 





        isize_nuc = 1
        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0




        ch2so4_crit = 0.16 * exp( 0.1*temp_in - 3.5*rh_in - 27.7 )


        qh2so4_crit = (ch2so4_crit*1.0e-12/98.0)/cair
        qh2so4_avail = qh2so4_cur - qh2so4_crit



        if (qh2so4_avail .le. 4.0e-18) then
           return
        end if


        isize_nuc = 1
        vol_part = volumlo_sect(1)









        if (qnh3_cur .ge. qh2so4_avail) then


           dum_n1 = (qnh3_cur/qh2so4_avail) - 1.0
           dum_n1 = max( 0.0_r8, min( 1.0_r8, dum_n1 ) )
           dum_n2 = 1.0 - dum_n1
           dum_n3 = 0.0
        else


           dum_n1 = 0.0
           dum_n2 = (qnh3_cur/qh2so4_avail)
           dum_n2 = max( 0.0_r8, min( 1.0_r8, dum_n2 ) )
           dum_n3 = 1.0 - dum_n2
	end if

        dum_m1 = dum_n1*mw_ammsulf
        dum_m2 = dum_n2*mw_ammbisulf
        dum_m3 = dum_n3*mw_sulfacid
        dens_part = (dum_m1 + dum_m2 + dum_m3)/   &
           ((dum_m1/dens_ammsulf) + (dum_m2/dens_ammbisulf)   &
                                  + (dum_m3/dens_sulfacid))

	if (abs(dens_nh4so4a-1.75) .le. 0.25) then
	    dens_part = dens_nh4so4a
	else
            dens_nh4so4a = dens_part
	end if
        mass_part  = vol_part*dens_part 
        molenh4a_per_moleso4a = 2.0*dum_n1 + dum_n2



        duma = 0.9999
        qh2so4_del = min( duma*qh2so4_cur, qh2so4_avail )
        qnh3_del   = min( duma*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del


        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del

        qnuma_del = (qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part


        return
        end subroutine wexler_nuc_mosaic_1box








	end module module_mosaic_newnucb
