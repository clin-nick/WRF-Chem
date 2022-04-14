







	module module_mosaic_newnuc



	use module_peg_util



	implicit none



	contains




	subroutine mosaic_newnuc_1clm( istat_newnuc,   &
		it, jt, kclm_calcbgn, kclm_calcend,   &
		idiagbb_in,   &
		dtchem, dtnuc_in, rsub0,   &
		id, ktau, ktauc, its, ite, jts, jte, kts, kte )












	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_state_description, only:  param_first_scalar


	integer, intent(inout) :: istat_newnuc    
	integer, intent(in) ::   &
		it, jt, kclm_calcbgn, kclm_calcend,   &
		idiagbb_in,   &
		id, ktau, ktauc, its, ite, jts, jte, kts, kte
	real, intent(in) :: dtchem, dtnuc_in
	real, intent(in) :: rsub0(l2maxd,kmaxd,nsubareamaxd)














	integer, parameter :: newnuc_method = 2

	integer :: k, l, ll, m
	integer :: isize, itype, iphase
	integer :: iconform_numb
	integer :: idiagbb
	integer :: nsize, ntau_nuc
	integer, save :: ncount(10)
	integer :: p1st

	real, parameter :: densdefault = 2.0
	real, parameter :: smallmassbb = 1.0e-30


        real, parameter :: a_zsr_xx1 =  1.15510
        real, parameter :: a_zsr_xx2 = -3.20815
        real, parameter :: a_zsr_xx3 =  2.71141
        real, parameter :: a_zsr_xx4 =  2.01155
        real, parameter :: a_zsr_xx5 = -4.71014
        real, parameter :: a_zsr_xx6 =  2.04616
        real, parameter :: b_zsr_xx  = 29.4779

	real :: aw
	real :: cair_box
	real :: dens_nh4so4a, dtnuc
	real :: duma, dumb, dumc
	real :: rh_box
	real :: qh2so4_avg, qh2so4_cur, qh2so4_del 
	real :: qnh3_avg, qnh3_cur, qnh3_del 
	real :: qnuma_del, qso4a_del, qnh4a_del
	real :: temp_box
	real :: xxdens, xxmass, xxnumb, xxvolu

	real,save :: dumveca(10), dumvecb(10), dumvecc(10), dumvecd(10), dumvece(10)
	real :: volumlo_nuc(maxd_asize), volumhi_nuc(maxd_asize)

	character(len=100) :: msg

    p1st = PARAM_FIRST_SCALAR


	istat_newnuc = 0
	if (newnuc_method .ne. 2) then
	    if ((it .eq. its) .and. (jt .eq. jts))   &
		call peg_message( lunerr,   &
		'*** mosaic_newnuc_1clm -- illegal newnuc_method' )
	    istat_newnuc = -1
	    return
	end if




	ntau_nuc = nint( dtnuc_in/dtchem )
	ntau_nuc = max( 1, ntau_nuc )
	if (mod(ktau,ntau_nuc) .ne. 0) return
	dtnuc = dtchem*ntau_nuc



	idiagbb = idiagbb_in

	itype = 1
	iphase = ai_phase
	nsize = nsize_aer(itype)
	volumlo_nuc(1:nsize) = volumlo_sect(1:nsize,itype)
	volumhi_nuc(1:nsize) = volumhi_sect(1:nsize,itype)



	do 2900 m = 1, nsubareas

	do 2800 k = kclm_calcbgn, kclm_calcend



	if ((it .eq. its) .and.   &
	    (jt .eq. jts) .and. (k .eq. kclm_calcbgn)) then
	    dumveca(:) = 0.0         
	    dumvecb(:) = +1.0e35     
	    dumvecc(:) = -1.0e35     
	    dumvecd(:) = 0.0         
	    dumvece(:) = 0.0         
	    ncount(:) = 0
	end if


	ncount(1) = ncount(1) + 1
	if (afracsubarea(k,m) .lt. 1.e-4) goto 2700

	cair_box = cairclm(k)
	temp_box = rsub(ktemp,k,m)
	rh_box = relhumclm(k)

	qh2so4_cur = max(0.0,rsub(kh2so4,k,m))
	qnh3_cur   = max(0.0,rsub(knh3,k,m))
	qh2so4_avg = 0.5*( qh2so4_cur + max(0.0,rsub0(kh2so4,k,m)) )
	qnh3_avg   = 0.5*( qnh3_cur   + max(0.0,rsub0(knh3,k,m)) )

	qh2so4_del = 0.0
	qnh3_del = 0.0
	qnuma_del = 0.0
	qso4a_del = 0.0
	qnh4a_del = 0.0

	dens_nh4so4a = dens_so4_aer

	isize = 0


	if (newnuc_method .eq. 1) then
            call ternary_nuc_mosaic_1box(   &
               dtnuc, temp_box, rh_box, cair_box,   &
               qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
               nsize, maxd_asize, volumlo_nuc, volumhi_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del,   &
               qh2so4_del, qnh3_del, dens_nh4so4a )
	else if (newnuc_method .eq. 2) then
            call wexler_nuc_mosaic_1box(   &
               dtnuc, temp_box, rh_box, cair_box,   &
               qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
               nsize, maxd_asize, volumlo_nuc, volumhi_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del,   &
               qh2so4_del, qnh3_del, dens_nh4so4a )
	else
	    istat_newnuc = -1
	    return
	end if



	dumveca(1) = temp_box
	dumveca(2) = rh_box
	dumveca(3) = rsub(kso2,k,m)
	dumveca(4) = qh2so4_avg
	dumveca(5) = qnh3_avg
	dumveca(6) = qnuma_del
	do l = 1, 6
	    dumvecb(l) = min( dumvecb(l), dumveca(l) )
	    dumvecc(l) = max( dumvecc(l), dumveca(l) )
	    dumvecd(l) = dumvecd(l) + dumveca(l)
	    if (qnuma_del .gt. dumvece(6)) dumvece(l) = dumveca(l)
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


	rsub(kh2so4,k,m) = max( 0.0, rsub(kh2so4,k,m) + qh2so4_del )
	rsub(knh3,  k,m) = max( 0.0, rsub(knh3,  k,m) + qnh3_del )

	l = lptr_so4_aer(isize,itype,iphase)
	if (l .ge. p1st) then
	    rsub(l,k,m) = rsub(l,k,m) + qso4a_del
	end if
	l = lptr_nh4_aer(isize,itype,iphase)
	if (l .ge. p1st) then
	    rsub(l,k,m) = rsub(l,k,m) + qnh4a_del
	end if
	l = numptr_aer(isize,itype,iphase)
	rsub(l,k,m) = rsub(l,k,m) + qnuma_del
	xxnumb = rsub(l,k,m)





	l = waterptr_aer(isize,itype)
	if ((rh_box .gt. 0.10) .and. (l .ge. p1st)) then
	    aw = min( rh_box, 0.98 )
	    if (aw .lt. 0.97) then
		duma =       a_zsr_xx1 +   &
		        aw*( a_zsr_xx2 +   &
		        aw*( a_zsr_xx3 +   &
		        aw*( a_zsr_xx4 +   &
		        aw*( a_zsr_xx5 +   &
		        aw*  a_zsr_xx6 ))))
	    else
		dumb = -b_zsr_xx*log(aw)
	        dumb = max( dumb, 0.5 )
		duma = 1.0/(1.0 + 55.509/dumb)
	    end if
	    duma = max( duma, 0.01 )
	    dumc = (1.0 - duma)/duma
	    rsub(l,k,m) = rsub(l,k,m) + qso4a_del*dumc
	end if






	xxmass = aqmassdry_sub(isize,itype,k,m)
	xxdens = adrydens_sub( isize,itype,k,m)
	iconform_numb = 1

	if ((xxdens .lt. 0.1) .or. (xxdens .gt. 20.0)) then

	    continue
	else


	    xxvolu = xxmass/xxdens
	    duma = qso4a_del*mw_so4_aer + qnh4a_del*mw_nh4_aer
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
		    duma = max( 0.0, rsub(l,k,m) )*mw_aer(ll,itype)
		    xxmass = xxmass + duma
		    xxvolu = xxvolu + duma/dens_aer(ll,itype)
		end if
	    end do
	end if

	if (xxmass .le. smallmassbb) then


	    ncount(5) = ncount(5) + 1
	    xxdens = densdefault
	    xxvolu = xxmass/xxdens
	    xxnumb = xxmass/(volumcen_sect(isize,itype)*xxdens)
	    iconform_numb = 0
	    l = waterptr_aer(isize,itype)
	    if (l .ge. p1st) rsub(l,k,m) = 0.0
	    l = hyswptr_aer(isize,itype)
	    if (l .ge. p1st) rsub(l,k,m) = 0.0
	else
	    xxdens = xxmass/xxvolu
	end if

	if (iconform_numb .gt. 0) then

	    if (xxnumb .gt. xxvolu/volumlo_sect(isize,itype)) then
		ncount(6) = ncount(6) + 1
		xxnumb = xxvolu/volumlo_sect(isize,itype)
	    else if (xxnumb .lt. xxvolu/volumhi_sect(isize,itype)) then
		ncount(7) = ncount(7) + 1
		xxnumb = xxvolu/volumhi_sect(isize,itype)
	    end if
	end if


	l = numptr_aer(isize,itype,iphase)
	rsub(l,k,m) = xxnumb
	adrydens_sub( isize,itype,k,m) = xxdens
	aqmassdry_sub(isize,itype,k,m) = xxmass
	aqvoldry_sub( isize,itype,k,m) = xxvolu


2700	continue


	if ((idiagbb .ge. 100) .and.   &
	    (it .eq. ite) .and.    &
	    (jt .eq. jte) .and. (k .eq. kclm_calcend)) then
	    if (idiagbb .ge. 110) then
	      write(msg,93020) 'newnucbb mins ', dumvecb(1:6)
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb maxs ', dumvecc(1:6)
	      call peg_message( lunerr, msg )
	      duma = max( 1, ncount(1) ) 
	      write(msg,93020) 'newnucbb avgs ', dumvecd(1:6)/duma
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb hinuc', dumvece(1:6)
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb dtnuc', dtnuc
	      call peg_message( lunerr, msg )
	    end if
	    write(msg,93030) 'newnucbb ncnt ', ncount(1:7)
	    call peg_message( lunerr, msg )
	end if
93020	format( a, 1p, 10e10.2 )
93030	format( a, 1p, 10i10 )


2800	continue	

2900	continue	


	return
	end subroutine mosaic_newnuc_1clm






        subroutine ternary_nuc_mosaic_1box(   &
           dtnuc, temp_in, rh_in, cair,   &
           qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
           nsize, maxd_asize, volumlo_sect, volumhi_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a )
























        implicit none


        real, intent(in) :: dtnuc             
        real, intent(in) :: temp_in           
        real, intent(in) :: rh_in             
        real, intent(in) :: cair              

        real, intent(in) :: qh2so4_avg, qh2so4_cur   
        real, intent(in) :: qnh3_avg, qnh3_cur       
             
             

        integer, intent(in) :: nsize                    
        integer, intent(in) :: maxd_asize               
        real, intent(in) :: volumlo_sect(maxd_asize)    
        real, intent(in) :: volumhi_sect(maxd_asize)    


        integer, intent(out) :: isize_nuc     
        real, intent(out) :: qnuma_del        
        real, intent(out) :: qso4a_del        
        real, intent(out) :: qnh4a_del        
        real, intent(out) :: qh2so4_del       
        real, intent(out) :: qnh3_del         
                                              


        real, intent(inout) :: dens_nh4so4a   





        real :: ratenuclt        
        real :: rateloge         
        real :: cnum_h2so4       
        real :: cnum_nh3         
        real :: cnum_tot         
        real :: radius_cluster   







        integer i
        integer, save :: icase = 0, icase_reldiffmax = 0

        real, parameter :: pi = 3.1415926536
        real, parameter :: avogad = 6.022e23   
        real, parameter :: mw_air = 28.966     



        real, parameter :: dens_ammsulf   = 1.769
        real, parameter :: dens_ammbisulf = 1.78
        real, parameter :: dens_sulfacid  = 1.841




        real, parameter :: mw_ammsulf   = 132.0
        real, parameter :: mw_ammbisulf = 114.0
        real, parameter :: mw_sulfacid  =  96.0

        real, parameter :: mw_so4a      =  96.0
        real, parameter :: mw_nh4a      =  18.0

        real, save :: reldiffmax = 0.0

        real dens_part                
        real duma, dumb, dumc, dume
        real dum_m1, dum_m2, dum_m3, dum_n1, dum_n2, dum_n3
        real fogas, foso4a, fonh4a, fonuma
        real freduce                  
                                      
        real freducea, freduceb
        real gramaero_per_moleso4a    
        real mass_part                
        real molenh4a_per_moleso4a    
        real nh3conc_in               
        real so4vol_in                
        real qmolnh4a_del_max         
        real qmolso4a_del_max         
        real vol_cluster              
        real vol_part                 






        isize_nuc = 1
        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0
        if (qh2so4_avg .le. 4.0e-16) return
        if (qh2so4_cur .le. 4.0e-16) return








        nh3conc_in = qnh3_avg/1.0e-12



        so4vol_in  = (qh2so4_avg) * cair * avogad

        call ternary_nuc_napari(   &
            temp_in, rh_in, nh3conc_in, so4vol_in,   &
            ratenuclt, rateloge,   &
            cnum_h2so4, cnum_nh3, cnum_tot, radius_cluster )



        if (ratenuclt .le. 1.0e-6) return




        vol_cluster = (pi*4.0/3.0)* (radius_cluster**3) * 1.0e-21
        isize_nuc = 1
        vol_part = volumlo_sect(1)
        if (vol_cluster .le. volumlo_sect(1)) then
           continue
        else if (vol_cluster .ge. volumhi_sect(nsize)) then
           isize_nuc = nsize
           vol_part = volumhi_sect(nsize)
        else
           do i = 1, nsize
              if (vol_cluster .lt. volumhi_sect(i)) then
                 isize_nuc = i
                 vol_part = vol_cluster
                 vol_part = min( vol_part, volumhi_sect(i) )
                 vol_part = max( vol_part, volumlo_sect(i) )
                 exit
              end if
           end do
        end if










        if (qnh3_cur .ge. qh2so4_cur) then


           dum_n1 = (qnh3_cur/qh2so4_cur) - 1.0
           dum_n1 = max( 0.0, min( 1.0, dum_n1 ) )
           dum_n2 = 1.0 - dum_n1
           dum_n3 = 0.0
        else


           dum_n1 = 0.0
           dum_n2 = (qnh3_cur/qh2so4_cur)
           dum_n2 = max( 0.0, min( 1.0, dum_n2 ) )
           dum_n3 = 1.0 - dum_n2
	end if

        dum_m1 = dum_n1*mw_ammsulf
        dum_m2 = dum_n2*mw_ammbisulf
        dum_m3 = dum_n3*mw_sulfacid
        dens_part = (dum_m1 + dum_m2 + dum_m3)/   &
           ((dum_m1/dens_ammsulf) + (dum_m2/dens_ammbisulf)   &
                                  + (dum_m3/dens_sulfacid))

	if (abs(dens_nh4so4a-1.8) .le. 0.2) then
	    dens_part = dens_nh4so4a
	else
            dens_nh4so4a = dens_part
	end if
        mass_part  = vol_part*dens_part 
        molenh4a_per_moleso4a = 2.0*dum_n1 + dum_n2  
        gramaero_per_moleso4a = dum_m1 + dum_m2 + dum_m3  



        duma = max( 0.0, (ratenuclt*dtnuc*mass_part) )

        dumc = duma/gramaero_per_moleso4a

        dume = dumc/cair


        qmolso4a_del_max = dume


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



        if (freduce*ratenuclt .le. 1.0e-6) return














        duma = 0.9999
        qh2so4_del = min( duma*qh2so4_cur, freduce*qmolso4a_del_max )
        qnh3_del   = min( duma*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del


        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del

        qnuma_del = (qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part


        return
        end subroutine ternary_nuc_mosaic_1box




        subroutine ternary_nuc_napari(   &
           temp_in, rh_in, nh3conc_in, so4vol_in,   &
           ratenuclt, rateloge,   &
           cnum_h2so4, cnum_nh3, cnum_tot, radius_cluster )












































      implicit none


        real, intent(in) :: temp_in           
        real, intent(in) :: rh_in             
        real, intent(in) :: nh3conc_in        
        real, intent(in) :: so4vol_in         


        real, intent(out) :: ratenuclt        
        real, intent(out) :: rateloge         

        real, intent(out) :: cnum_h2so4       
                                              
        real, intent(out) :: cnum_nh3         
                                              
        real, intent(out) :: cnum_tot         
                                              
        real, intent(out) :: radius_cluster   


        integer ncoeff
        parameter ( ncoeff = 4 )      
                                      
        integer npoly
        parameter ( npoly = 20 )      

        integer n                     


        real f ( npoly )


        real a  ( ncoeff, npoly )

        real temp, rh, nh3conc, so4vol     
        real log_rh, log_nh3conc, log_so4vol


        data a  / -0.355297,  -33.8449,      0.34536,    -8.24007e-4,   &
                   3.13735,    -0.772861,    5.61204e-3, -9.74576e-6,   &
                  19.0359,     -0.170957,    4.79808e-4, -4.14699e-7,   &
                   1.07605,     1.48932,    -7.96052e-3,  7.61229e-6,   &
                   6.0916,     -1.25378,     9.39836e-3, -1.74927e-5,   &
                   0.31176,     1.64009,    -3.43852e-3, -1.09753e-5,   &
                  -2.00738e-2, -0.752115,    5.25813e-3, -8.98038e-6,   &
                   0.165536,    3.26623,    -4.89703e-2,  1.46967e-4,   &
                   6.52645,    -0.258002,    1.43456e-3, -2.02036e-6,   &
                   3.68024,    -0.204098,    1.06259e-3, -1.2656e-6 ,   &
                  -6.6514e-2,  -7.82382,     1.22938e-2,  6.18554e-5,   &
                   0.65874,     0.190542,   -1.65718e-3,  3.41744e-6,   &
                   5.99321e-2,  5.96475,    -3.62432e-2,  4.93337e-5,   &
                  -0.732731,   -1.84179e-2,  1.47186e-4, -2.37711e-7,   &
                   0.728429,    3.64736,    -2.7422e-2,   4.93478e-5,   &
                  41.3016,     -0.35752,     9.04383e-4, -5.73788e-7,   &
                  -0.160336,    8.89881e-3, -5.39514e-5,  8.39522e-8,   &
                   8.57868,    -0.112358,    4.72626e-4, -6.48365e-7,   &
                   5.30167e-2, -1.98815,     1.57827e-2, -2.93564e-5,   &
                  -2.32736,     2.34646e-2, -7.6519e-5,   8.0459e-8 /





        temp    = max( 240.15, min (300.15, temp_in ) )
        rh      = max( 0.05,   min (0.95,   rh_in ) )
        so4vol  = max( 1.0e4,  min (1.0e9,  so4vol_in ) )
        nh3conc = max( 0.1,    min (100.0,  nh3conc_in ) )





        do n = 1, npoly
            f ( n )   = a ( 1, n ) + a ( 2, n ) * temp   &
                      + a ( 3, n ) * ( temp ) ** 2.0   &
                      + a ( 4, n ) * ( temp ) ** 3.0
        end do




        log_rh = log ( rh )
        log_nh3conc = log ( nh3conc )
        log_so4vol = log ( so4vol )
        rateloge = -84.7551   &
                 + f ( 1 ) / log_so4vol   &
                 + f ( 2 )  * ( log_so4vol )   &
                 + f ( 3 )  * ( log_so4vol ) **2.0   &
                 + f ( 4 )  * ( log_nh3conc )   &
                 + f ( 5 )  * ( log_nh3conc ) **2.0   &
                 + f ( 6 )  * rh   &
                 + f ( 7 )  * ( log_rh )   &
                 + f ( 8 )  * ( log_nh3conc /   &
                              log_so4vol )   &
                 + f ( 9 )  * ( log_nh3conc *   &
                              log_so4vol )   &
                 + f ( 10 ) * rh  *   &
                              ( log_so4vol )   &
                 + f ( 11 ) * ( rh /   &
                              log_so4vol )   &
                 + f ( 12 ) * ( rh *   &
                              log_nh3conc )   &
                 + f ( 13 ) * ( log_rh /   &
                              log_so4vol )   &
                 + f ( 14 ) * ( log_rh *   &
                              log_nh3conc )   &
                 + f ( 15 ) * (( log_nh3conc ) ** 2.0   &
                              / log_so4vol )   &
                 + f ( 16 ) * ( log_so4vol *   &
                              ( log_nh3conc ) ** 2.0 )   &
                 + f ( 17 ) * (( log_so4vol ) ** 2.0 *   &
                              log_nh3conc )   &
                 + f ( 18 ) * ( rh  *   &
                              ( log_nh3conc ) ** 2.0 )   &
                 + f ( 19 ) * ( rh  *  log_nh3conc   &
                              / log_so4vol )   &
                 + f ( 20 ) * (( log_so4vol ) ** 2.0 *   &
                              ( log_nh3conc ) ** 2.0 )

        ratenuclt = exp ( rateloge )



        cnum_h2so4 = 38.1645 + 0.774106 * rateloge   &
                   + 0.00298879 * ( rateloge ) ** 2.0   &
                   - 0.357605 * temp   &
                   - 0.00366358 * temp * rateloge   &
                   + 0.0008553 * ( temp ) ** 2.0

        cnum_nh3   = 26.8982 + 0.682905 * rateloge   &
                   + 0.00357521 * ( rateloge ) ** 2.0   &
                   - 0.265748 * temp   &
                   - 0.00341895 * temp * rateloge   &
                   + 0.000673454 * ( temp ) ** 2.0

        cnum_tot   = 79.3484 + 1.7384 * rateloge   &
                   + 0.00711403 * ( rateloge ) ** 2.0   &
                   - 0.744993 * temp   &
                   - 0.00820608 * temp * rateloge   &
                   + 0.0017855 * ( temp ) ** 2.0

        radius_cluster = 0.141027 - 0.00122625 * rateloge   &
                   - 7.82211e-6 * ( rateloge ) ** 2.0   &
                   - 0.00156727 * temp   &
                   - 0.00003076 * temp * rateloge   &
                   + 0.0000108375 * ( temp ) ** 2.0

        return
        end subroutine ternary_nuc_napari




        subroutine wexler_nuc_mosaic_1box(   &
           dtnuc, temp_in, rh_in, cair,   &
           qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
           nsize, maxd_asize, volumlo_sect, volumhi_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a )






















        implicit none


        real, intent(in) :: dtnuc             
        real, intent(in) :: temp_in           
        real, intent(in) :: rh_in             
        real, intent(in) :: cair              

        real, intent(in) :: qh2so4_avg, qh2so4_cur   
        real, intent(in) :: qnh3_avg, qnh3_cur       
             
             

        integer, intent(in) :: nsize                    
        integer, intent(in) :: maxd_asize               
        real, intent(in) :: volumlo_sect(maxd_asize)    
        real, intent(in) :: volumhi_sect(maxd_asize)    


        integer, intent(out) :: isize_nuc     
        real, intent(out) :: qnuma_del        
        real, intent(out) :: qso4a_del        
        real, intent(out) :: qnh4a_del        
        real, intent(out) :: qh2so4_del       
        real, intent(out) :: qnh3_del         
                                              


        real, intent(inout) :: dens_nh4so4a   



        integer i
        integer, save :: icase = 0, icase_reldiffmax = 0

        real, parameter :: pi = 3.1415926536
        real, parameter :: avogad = 6.022e23   
        real, parameter :: mw_air = 28.966     



        real, parameter :: dens_ammsulf   = 1.769
        real, parameter :: dens_ammbisulf = 1.78
        real, parameter :: dens_sulfacid  = 1.841




        real, parameter :: mw_ammsulf   = 132.0
        real, parameter :: mw_ammbisulf = 114.0
        real, parameter :: mw_sulfacid  =  96.0

        real, parameter :: mw_so4a      =  96.0
        real, parameter :: mw_nh4a      =  18.0

        real, save :: reldiffmax = 0.0

        real ch2so4_crit              
        real dens_part                
        real duma, dumb, dumc, dume
        real dum_m1, dum_m2, dum_m3, dum_n1, dum_n2, dum_n3
        real fogas, foso4a, fonh4a, fonuma
        real mass_part                
        real molenh4a_per_moleso4a    
        real qh2so4_crit              
        real qh2so4_avail             
        real vol_part                 





        isize_nuc = 1
        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0




        ch2so4_crit = 0.16 * exp( 0.1*temp_in - 3.5*rh_in - 27.7 )


        qh2so4_crit = (ch2so4_crit*1.0e-12/98.0)/cair
        qh2so4_avail = (qh2so4_cur - qh2so4_crit)*dtnuc



        if (qh2so4_avail .le. 4.0e-18) then
           return
        end if


        isize_nuc = 1
        vol_part = volumlo_sect(1)









        if (qnh3_cur .ge. qh2so4_avail) then


           dum_n1 = (qnh3_cur/qh2so4_avail) - 1.0
           dum_n1 = max( 0.0, min( 1.0, dum_n1 ) )
           dum_n2 = 1.0 - dum_n1
           dum_n3 = 0.0
        else


           dum_n1 = 0.0
           dum_n2 = (qnh3_cur/qh2so4_avail)
           dum_n2 = max( 0.0, min( 1.0, dum_n2 ) )
           dum_n3 = 1.0 - dum_n2
	end if

        dum_m1 = dum_n1*mw_ammsulf
        dum_m2 = dum_n2*mw_ammbisulf
        dum_m3 = dum_n3*mw_sulfacid
        dens_part = (dum_m1 + dum_m2 + dum_m3)/   &
           ((dum_m1/dens_ammsulf) + (dum_m2/dens_ammbisulf)   &
                                  + (dum_m3/dens_sulfacid))

	if (abs(dens_nh4so4a-1.8) .le. 0.2) then
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








	end module module_mosaic_newnuc
