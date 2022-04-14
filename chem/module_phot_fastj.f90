
































































	module module_phot_fastj
	integer, parameter :: lunerr = -1

	contains

       subroutine fastj_driver(id,curr_secs,dtstep,config_flags,       &
               gmt,julday,t_phy,moist,p8w,p_phy,                       &
               chem,rho_phy,dz8w,xlat,xlong,z_at_w,                    &
               ph_o2,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,  &
               ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,    &
               ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,         &
               ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,&
               ph_n2o5,                                                &
               tauaer1,tauaer2,tauaer3,tauaer4,                        &
               gaer1,gaer2,gaer3,gaer4,                                &
               waer1,waer2,waer3,waer4,                                &
               bscoef1,bscoef2,bscoef3,bscoef4,                        &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                    &
               ids,ide, jds,jde, kds,kde,                              &
               ims,ime, jms,jme, kms,kme,                              &
               its,ite, jts,jte, kts,kte                               )


   USE module_configure
   USE module_state_description
   USE module_data_mosaic_therm, only: nbin_a, nbin_a_maxd

   USE module_data_mosaic_other, only: kmaxd, nsubareas
   USE module_fastj_mie


   USE module_fastj_data, only:  nb, nc

   implicit none


   integer, parameter :: single = 4        
   integer, parameter :: double = 8        
   integer,parameter :: ipar_fastj=1,jpar=1
   integer,parameter :: jppj=14        
   logical,parameter :: ldeg45=.false. 
   integer,save :: lpar           
   integer,save :: jpnl           
   real(kind=double), dimension(ipar_fastj) :: xgrd 
   real(kind=double), dimension(jpar) :: ygrd           
   real(kind=double), dimension(jpar) :: ydgrd          
   real(kind=double), dimension(kmaxd+1) :: etaa    
   real(kind=double), dimension(kmaxd+1) :: etab    
   real(kind=double) ::  tau_fastj          
   integer month_fastj        
   integer iday_fastj         
   integer nspint           
   parameter ( nspint = 4 ) 
   real, dimension (nspint),save :: wavmid 
   real, dimension (nspint, kmaxd+1),save :: sizeaer,extaer,waer,gaer,tauaer,bscoef
   real, dimension (nspint, kmaxd+1),save :: l2,l3,l4,l5,l6,l7
   data wavmid     &
       / 0.30e-4, 0.40e-4, 0.60e-4 ,0.999e-04 /
   integer nphoto_fastj
   parameter (nphoto_fastj = 14)
   integer   &
   lfastj_no2,   lfastj_o3a,   lfastj_o3b,    lfastj_h2o2,   &
   lfastj_hchoa, lfastj_hchob, lfastj_ch3ooh, lfastj_no3x,   &
   lfastj_no3l,  lfastj_hono,  lfastj_n2o5,   lfastj_hno3,   &
   lfastj_hno4     
   parameter( lfastj_no2   = 1 )
   parameter( lfastj_o3a   = 2 )
   parameter( lfastj_o3b   = 3 )
   parameter( lfastj_h2o2  = 4 )
   parameter( lfastj_hchoa = 5 )
   parameter( lfastj_hchob = 6 )
   parameter( lfastj_ch3ooh= 7 )
   parameter( lfastj_no3x  = 8 )
   parameter( lfastj_no3l  = 9 )
   parameter( lfastj_hono  = 10 )
   parameter( lfastj_n2o5  = 11 )
   parameter( lfastj_hno3  = 12 )
   parameter( lfastj_hno4  = 13 )

   INTEGER,      INTENT(IN   ) :: id,julday,                    &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte
   REAL(KIND=8), INTENT(IN   ) :: curr_secs
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                    &
         INTENT(INOUT ) ::                                          &
           ph_o2,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,   &
           ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,     &
           ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,          &
           ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob, &
           ph_n2o5
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                    &
         INTENT(IN ) ::                                             &
           tauaer1,tauaer2,tauaer3,tauaer4,                         &
           gaer1,gaer2,gaer3,gaer4,                                 &
           waer1,waer2,waer3,waer4,                                 &
           bscoef1,bscoef2,bscoef3,bscoef4                          
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:4 ),               &
         INTENT(IN ) ::                                             &
           l2aer,l3aer,l4aer,l5aer,l6aer,l7aer                      

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),         &
         INTENT(INOUT ) ::                                chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                      t_phy,    &
                                                      p_phy,    &
                                                       dz8w,    &
                                                        p8w,    &
                                                    rho_phy,    &
                                                      z_at_w
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,    &     
          INTENT(IN   ) ::                             xlat,    &
                                                      xlong
   REAL,      INTENT(IN   ) ::                                  &
                                                 dtstep,gmt

   TYPE(grid_config_rec_type), INTENT(IN  ) :: config_flags




      integer iclm, jclm

      real number_bin(nbin_a_maxd,kmaxd) 
      real radius_wet(nbin_a_maxd,kmaxd) 
      complex refindx(nbin_a_maxd,kmaxd)

      integer(kind=8) :: ixhour
      integer i,j, k, nsub
	  real(kind=8) :: xtime, xhour
      real xmin, gmtp, tmidh
      real sla, slo
      real psfc

      real cos_sza
      real, dimension(kts:kte) :: temp, ozone, dz 
      real, dimension(0:kte) :: pbnd
      real, dimension(kts:kte) :: cloudmr, airdensity, relhum
      real, dimension(kts:kte+1) :: zatw
 	     
      real valuej(kte,nphoto_fastj)

      logical processingAerosols




    lpar = kte
    jpnl = kte
    nb   = lpar + 1 
    nc   = 2*nb     
	nsubareas = 1
	if ((kte > kmaxd) .or. (lpar <= 0)) then
	    write( wrf_err_message, '(a,4i5)' )   &
		'*** subr fastj_driver -- ' //   &
		'lpar, kmaxd, kts, kte', lpar, kmaxd, kts, kte
        call wrf_message( trim(wrf_err_message) )
	    wrf_err_message = '*** subr fastj_driver -- ' //   &
             'kte>kmaxd OR lpar<=0'
	    call wrf_error_fatal3("<stdin>",228,&
wrf_err_message )
	end if




    select case (config_flags%chem_opt)
    case ( CBMZ_MOSAIC_4BIN,    CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_KPP,  &
           CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
           CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, &
           CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
           MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP, &
           CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, &
           SAPRC99_MOSAIC_8BIN_VBS2_KPP )
       processingAerosols = .true.
    case default
       processingAerosols = .false.
    end select



















    xtime    = curr_secs/60._8 + real(dtstep,8)/120._8
    ixhour   = int(gmt + 0.01,8) + int(xtime/60._8,8)
    xhour    = real(ixhour,8)	
    xmin     = 60.*gmt + real(xtime-xhour*60_8,8)
    gmtp     = mod(xhour,24._8)
    tmidh    = gmtp + xmin/60.


    do nsub = 1, nsubareas
	do jclm = jts, jte
	do iclm = its, ite

       do k = kts, lpar
          dz(k) = dz8w(iclm, k, jclm)	
       end do

       if( processingAerosols ) then
         do k = kts, lpar
           l2(1,k)=l2aer(iclm,k,jclm,1)
           l2(2,k)=l2aer(iclm,k,jclm,2)
           l2(3,k)=l2aer(iclm,k,jclm,3)
           l2(4,k)=l2aer(iclm,k,jclm,4)
           l3(1,k)=l3aer(iclm,k,jclm,1)
           l3(2,k)=l3aer(iclm,k,jclm,2)
           l3(3,k)=l3aer(iclm,k,jclm,3)
           l3(4,k)=l3aer(iclm,k,jclm,4)
           l4(1,k)=l4aer(iclm,k,jclm,1)
           l4(2,k)=l4aer(iclm,k,jclm,2)
           l4(3,k)=l4aer(iclm,k,jclm,3)
           l4(4,k)=l4aer(iclm,k,jclm,4)
           l5(1,k)=l5aer(iclm,k,jclm,1)
           l5(2,k)=l5aer(iclm,k,jclm,2)
           l5(3,k)=l5aer(iclm,k,jclm,3)
           l5(4,k)=l5aer(iclm,k,jclm,4)
           l6(1,k)=l6aer(iclm,k,jclm,1)
           l6(2,k)=l6aer(iclm,k,jclm,2)
           l6(3,k)=l6aer(iclm,k,jclm,3)
           l6(4,k)=l6aer(iclm,k,jclm,4)
           l7(1,k)=l7aer(iclm,k,jclm,1)
           l7(2,k)=l7aer(iclm,k,jclm,2)
           l7(3,k)=l7aer(iclm,k,jclm,3)
           l7(4,k)=l7aer(iclm,k,jclm,4)
         enddo


















       end if


	  sla = xlat(iclm,jclm)
	  slo = xlong(iclm,jclm)

	  psfc = p8w(iclm,1,jclm) * 10. 
	  do k = kts, lpar
	    pbnd(k) = p8w(iclm,k+1,jclm) *10.  
	    temp(k) = t_phy(iclm,k,jclm)
	    ozone(k) = chem(iclm,k,jclm,p_o3) / 1.0e6	
        cloudmr(k) = moist(iclm,k,jclm,p_qc)/0.622
        airdensity(k) = rho_phy(iclm,k,jclm)/28.966e3
        relhum(k) = MIN( .95, moist(iclm,k,jclm,p_qv) / &
                       (3.80*exp(17.27*(t_phy(iclm,k,jclm)-273.)/ &
                       (t_phy(iclm,k,jclm)-36.))/(.01*p_phy(iclm,k,jclm))))
        relhum(k) = MAX(.001,relhum(k))
        zatw(k)=z_at_w(iclm,k,jclm)
	  end do
      zatw(lpar+1)=z_at_w(iclm,lpar+1,jclm)

	  CALL wrf_debug(250,'fastj_driver: calling interface_fastj')
	  call interface_fastj(tmidh,sla,slo,julday,           &
           pbnd, psfc, temp, ozone,                        &
           dz, cloudmr, airdensity, relhum, zatw,          &
           iclm, jclm, lpar, jpnl,                         &
           curr_secs, valuej, cos_sza, processingAerosols, &
           sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)

      CALL wrf_debug(250,'fastj_driver: calling mapJrates_tofrom_host')
      call mapJrates_tofrom_host( 0,                         	    &
           ims,ime, jms,jme, kms,kme,                    		    &
           its,ite, jts,jte, kts,kte,                    		    &
           iclm, jclm, kts,lpar,                 		            &
           valuej,  						                        &
           ph_o2,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,   &
           ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,   	&
           ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,        	&
           ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob, &
           ph_n2o5                                                  )





















    end do
    end do
    end do
    
    return
    end subroutine  fastj_driver	     


	subroutine mapJrates_tofrom_host( iflag,                        &
		ims,ime, jms,jme, kms,kme,                    		        &
		its,ite, jts,jte, kts,kte,                    		        &
        iclm, jclm, ktmaps,ktmape,              		            &
        valuej,  						                            &
        ph_o2,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,      &
        ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,   	    &
        ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,        	    &
        ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,    &
        ph_n2o5                                                     )

	USE module_data_cbmz

   implicit none

   integer nphoto_fastj
   parameter (nphoto_fastj = 14)
   integer   &
   lfastj_no2,   lfastj_o3a,   lfastj_o3b,    lfastj_h2o2,   &
   lfastj_hchoa, lfastj_hchob, lfastj_ch3ooh, lfastj_no3x,   &
   lfastj_no3l,  lfastj_hono,  lfastj_n2o5,   lfastj_hno3,   &
   lfastj_hno4     
   parameter( lfastj_no2   = 1 )
   parameter( lfastj_o3a   = 2 )
   parameter( lfastj_o3b   = 3 )
   parameter( lfastj_h2o2  = 4 )
   parameter( lfastj_hchoa = 5 )
   parameter( lfastj_hchob = 6 )
   parameter( lfastj_ch3ooh= 7 )
   parameter( lfastj_no3x  = 8 )
   parameter( lfastj_no3l  = 9 )
   parameter( lfastj_hono  = 10 )
   parameter( lfastj_n2o5  = 11 )
   parameter( lfastj_hno3  = 12 )
   parameter( lfastj_hno4  = 13 )

   INTEGER,      INTENT(IN   ) :: iflag,                        &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte, 	&
                                  iclm, jclm, ktmaps, ktmape
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                   &
         INTENT(INOUT ) ::                                         &
           ph_o2,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,  &
           ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,    &
           ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,         &
           ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,&
           ph_n2o5

   REAL, DIMENSION( kte,nphoto_fastj ), INTENT(INOUT) :: valuej


	real ft
	integer kt

	ft = 60.

    if (iflag .gt. 0) go to 2000

	do kt = ktmaps, ktmape
	  ph_no2(iclm,kt,jclm)       = valuej(kt,lfastj_no2) * ft
	  ph_no3o(iclm,kt,jclm)      = valuej(kt,lfastj_no3x) * ft
	  ph_no3o2(iclm,kt,jclm)     = valuej(kt,lfastj_no3l) * ft
	  ph_o33p(iclm,kt,jclm)      = valuej(kt,lfastj_o3a) * ft
	  ph_o31d(iclm,kt,jclm)      = valuej(kt,lfastj_o3b) * ft
	  ph_hno2(iclm,kt,jclm)      = valuej(kt,lfastj_hono) * ft
	  ph_hno3(iclm,kt,jclm)      = valuej(kt,lfastj_hno3) * ft
	  ph_hno4(iclm,kt,jclm)      = valuej(kt,lfastj_hno4) * ft
	  ph_h2o2(iclm,kt,jclm)      = valuej(kt,lfastj_h2o2) * ft
	  ph_ch3o2h(iclm,kt,jclm)    = valuej(kt,lfastj_ch3ooh) * ft
	  ph_ch2or(iclm,kt,jclm)     = valuej(kt,lfastj_hchoa) * ft
	  ph_ch2om(iclm,kt,jclm)     = valuej(kt,lfastj_hchob) * ft
	  ph_n2o5(iclm,kt,jclm)      = valuej(kt,lfastj_n2o5) * ft

	  ph_o2(iclm,kt,jclm)        = 0.0
	  ph_ch3cho(iclm,kt,jclm)    = 0.0
	  ph_ch3coch3(iclm,kt,jclm)  = 0.0
	  ph_ch3coc2h5(iclm,kt,jclm) = 0.0
	  ph_hcocho(iclm,kt,jclm)    = 0.0
	  ph_ch3cocho(iclm,kt,jclm)  = 0.0
	  ph_hcochest(iclm,kt,jclm)  = 0.0
	  ph_ch3coo2h(iclm,kt,jclm)  = 0.0
	  ph_ch3ono2(iclm,kt,jclm)   = 0.0
	  ph_hcochob(iclm,kt,jclm)   = 0.0

	end do
	return 	

2000	continue

 	do kt = ktmaps, ktmape
	  valuej(kt,lfastj_no2)    = ph_no2(iclm,kt,jclm) /  ft
	  valuej(kt,lfastj_no3x)   = ph_no3o(iclm,kt,jclm) / ft
	  valuej(kt,lfastj_no3l)   = ph_no3o2(iclm,kt,jclm)/ ft
	  valuej(kt,lfastj_o3a)    = ph_o33p(iclm,kt,jclm) / ft
	  valuej(kt,lfastj_o3b)    = ph_o31d(iclm,kt,jclm) / ft
	  valuej(kt,lfastj_hono)   = ph_hno2(iclm,kt,jclm) / ft
	  valuej(kt,lfastj_hno3)   = ph_hno3(iclm,kt,jclm) / ft
	  valuej(kt,lfastj_hno4)   = ph_hno4(iclm,kt,jclm) / ft
	  valuej(kt,lfastj_h2o2)   = ph_h2o2(iclm,kt,jclm) / ft
	  valuej(kt,lfastj_ch3ooh) = ph_ch3o2h(iclm,kt,jclm)/ft
	  valuej(kt,lfastj_hchoa)  = ph_ch2or(iclm,kt,jclm)/ ft
	  valuej(kt,lfastj_hchob)  = ph_ch2om(iclm,kt,jclm)/ ft
	  valuej(kt,lfastj_n2o5)   = ph_n2o5(iclm,kt,jclm) / ft
	end do

    return	

    end subroutine mapJrates_tofrom_host



      	     	

        subroutine interface_fastj(tmidh,sla,slo,julian_day,   &
             pbnd, psfc, temp, ozone,                          &
             dz, cloudmr, airdensity, relhum, zatw,            &
             isvode, jsvode, lpar, jpnl,                       &
      	     curr_secs, valuej, cos_sza, processingAerosols,   &
             sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)














































        USE module_data_mosaic_other, only : kmaxd
	USE module_peg_util, only:  peg_message, peg_error_fatal
	
 	IMPLICIT NONE


   integer, parameter :: iprint = 0
   integer, parameter :: single = 4        
   integer, parameter :: double = 8        
   integer,parameter :: ipar_fastj=1,jpar=1
   integer,parameter :: jppj=14        
   logical,parameter :: ldeg45=.false. 
   integer lpar           
   integer jpnl           
   real(kind=double), dimension(ipar_fastj) :: xgrd 
   real(kind=double), dimension(jpar) :: ygrd           
   real(kind=double), dimension(jpar) :: ydgrd          
   real(kind=double), dimension(kmaxd+1) :: etaa    
   real(kind=double), dimension(kmaxd+1) :: etab    
   real(kind=double) ::  tau_fastj          
   integer month_fastj        
   integer iday_fastj         
   integer nphoto_fastj
   parameter (nphoto_fastj = 14)
   integer   &
   lfastj_no2,   lfastj_o3a,   lfastj_o3b,    lfastj_h2o2,   &
   lfastj_hchoa, lfastj_hchob, lfastj_ch3ooh, lfastj_no3x,   &
   lfastj_no3l,  lfastj_hono,  lfastj_n2o5,   lfastj_hno3,   &
   lfastj_hno4
   parameter( lfastj_no2   = 1 )
   parameter( lfastj_o3a   = 2 )
   parameter( lfastj_o3b   = 3 )
   parameter( lfastj_h2o2  = 4 )
   parameter( lfastj_hchoa = 5 )
   parameter( lfastj_hchob = 6 )
   parameter( lfastj_ch3ooh= 7 )
   parameter( lfastj_no3x  = 8 )
   parameter( lfastj_no3l  = 9 )
   parameter( lfastj_hono  = 10 )
   parameter( lfastj_n2o5  = 11 )
   parameter( lfastj_hno3  = 12 )
   parameter( lfastj_hno4  = 13 )
   integer nspint           
   parameter ( nspint = 4 ) 
   real, dimension (nspint),save :: wavmid 
   real, dimension (nspint, kmaxd+1) :: sizeaer,extaer,waer,gaer,tauaer
   real, dimension (nspint, kmaxd+1) :: l2,l3,l4,l5,l6,l7
   data wavmid     &
       / 0.30e-4, 0.40e-4, 0.60e-4 ,0.999e-04 /

	real pbnd(0:lpar), psfc
	real temp(lpar), ozone(lpar), surface_albedo
	real dz(lpar), cloudmr(lpar), airdensity(lpar), relhum(lpar), zatw(lpar+1)
    real(kind=8) :: curr_secs
	integer isvode, jsvode
	
	real cos_sza
    integer,parameter :: lunout=41
		
	real valuej(lpar,nphoto_fastj)
	
    real hl,rhl,factor1,part1,part2,cfrac,rhfrac
    real emziohl(lpar+1),clwp(lpar)

 	real valuej_no3rate(lpar)

	real*8 lat,lon
    real*8 jvalue(lpar,nphoto_fastj)
    real sza
	real tau1
	real tmidh, sla, slo	

	integer julian_day,iozone1
    integer,parameter :: nfastj_rxns = 14
	integer k, l
 		
	real surface_pressure_mb, tauaer_550,   &
         col_press_mb,col_temp_K,col_ozone,col_optical_depth
	dimension col_press_mb(lpar+2),col_temp_K(lpar+1),   &
        	col_ozone(lpar+1),col_optical_depth(lpar+1)
    character*80 msg	









      logical processingAerosols



	lat  = sla
	lon  = slo










    hl=1080.+2000.0*cos(lat*0.017454329)
    rhl=1.0/hl
	do k =1, lpar+1
       emziohl(k)=exp(-zatw(k)*rhl)
    enddo
	do k =1, lpar
       clwp(k)=0.18*hl*(emziohl(k)-emziohl(k+1))
    enddo



    factor1=1500.0
	do k =1, lpar
	   col_optical_depth(k) = 0.0 		
       cfrac=0.0
       cloudmr(k)=0.0
       if(cloudmr(k).gt.0.0) cfrac=1.0

       part1=cloudmr(k)*cfrac*18.0*airdensity(k)*dz(k)*100.0
       if(relhum(k).lt.0.8) then
          rhfrac=0.0
       elseif(relhum(k).le.1.0.and.relhum(k).ge.0.8) then

          rhfrac=(relhum(k)-0.8)/0.2
       else
          rhfrac=1.0
       endif
       if(rhfrac.ge.0.01) then

          part2=rhfrac*clwp(k)/1.e4
       else
          part2=0.0
       endif
       if(cfrac.gt.0) part2=0.0
       col_optical_depth(k) = factor1*(part1+part2)




	end do
    col_optical_depth(lpar+1) = 0.0 		
	if (.not.processingAerosols) then

	   call set_common_mie(lpar, &
           sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)
	end if      			 





	surface_pressure_mb = psfc * 0.001
	tau1 = tmidh
	col_press_mb(1) = psfc * 0.001	
	iozone1 = lpar
	do k =1, lpar
       col_press_mb(k+1) = pbnd(k) * 0.001
	   col_temp_K(k) = temp(k)
	   col_ozone(k) = ozone(k)
	end do

 	surface_albedo=0.055


        if (processingAerosols) then
            tauaer_550 = 0.0 	
              			
              			
         else
            tauaer_550 = 0.05 	
         end if
	
	  CALL wrf_debug(250,'interface_fastj: calling fastj')
      call fastj(isvode,jsvode,lat,lon,surface_pressure_mb,surface_albedo,   &
           julian_day,  tau1,   &
          col_press_mb, col_temp_K, col_optical_depth, col_ozone,   &
          iozone1,tauaer_550,jvalue,sza,lpar,jpnl, &
          sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)
     	
	
	cos_sza = cos(sza*3.141592653/180.)
	




	   do k = 1, lpar
	     valuej(k, lfastj_no2) = jvalue(k,lfastj_no2)
	     valuej(k, lfastj_o3a) = jvalue(k,lfastj_o3a)
	     valuej(k, lfastj_o3b) = jvalue(k,lfastj_o3b)
	     valuej(k, lfastj_h2o2) = jvalue(k,lfastj_h2o2)
	     valuej(k, lfastj_hchoa) = jvalue(k,lfastj_hchoa)
	     valuej(k, lfastj_hchob) = jvalue(k,lfastj_hchob)
	     valuej(k, lfastj_ch3ooh) = jvalue(k,lfastj_ch3ooh)
	     valuej(k, lfastj_no3x) = jvalue(k,lfastj_no3x)
	     valuej(k, lfastj_no3l) = jvalue(k,lfastj_no3l)
	     valuej(k, lfastj_hono) = jvalue(k,lfastj_hono)
	     valuej(k, lfastj_n2o5) = jvalue(k,lfastj_n2o5)
	     valuej(k, lfastj_hno3) = jvalue(k,lfastj_hno3)
	     valuej(k, lfastj_hno4) = jvalue(k,lfastj_hno4)
	  end do

	  do k = 1, lpar
         valuej(k,nphoto_fastj)=0.0
         do l = 1, nphoto_fastj-1
            if (valuej(k,l) .lt. 0) then
               write( msg, '(a,f14.2,4i4,1x,e11.4)' )   &
                    'FASTJ negative Jrate ' //   &
                    'tsec i j k l J(k,l)', curr_secs,isvode,jsvode,k,l,valuej(k,l)
               call peg_message( lunerr, msg )
               valuej(k,l) = 0.0




            end if
         end do
      end do







	end subroutine interface_fastj                          


	subroutine set_common_mie(lpar, &
          sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)

















        USE module_data_mosaic_other, only : kmaxd
	
	IMPLICIT NONE

        integer lpar             
        integer nspint           
        parameter ( nspint = 4 ) 
        real, dimension (nspint),save :: wavmid 
        real, dimension (nspint, kmaxd+1) :: sizeaer,extaer,waer,gaer,tauaer
        real, dimension (nspint, kmaxd+1) :: l2,l3,l4,l5,l6,l7
        data wavmid     &
            / 0.30e-4, 0.40e-4, 0.60e-4 ,0.999e-04 /



        integer klevel   
        integer ns       



        do 1000 ns=1,nspint
        do 1000 klevel = 1, lpar
          tauaer(ns,klevel)=0.
          waer(ns,klevel)=0.
          gaer(ns,klevel)=0.
	  sizeaer(ns,klevel)=0.0
	  extaer(ns,klevel)=0.0
	  l2(ns,klevel)=0.0
	  l3(ns,klevel)=0.0
  	  l4(ns,klevel)=0.0
	  l5(ns,klevel)=0.0
	  l6(ns,klevel)=0.0
	  l7(ns,klevel)=0.0
1000    continue

	return
	end subroutine set_common_mie      

    subroutine fastj(isvode,jsvode,lat,lon,surface_pressure,surface_albedo,   &
      julian_day,tau1, pressure, temperature, optical_depth, my_ozone1,   &
      iozone1,tauaer_550_1,jvalue,SZA_dum,lpar,jpnl, &
      sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)





































        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	
	IMPLICIT NONE


       integer, parameter :: iprint = 0
       integer, parameter :: single = 4        



       integer,parameter :: ipar_fastj=1,jpar=1

       logical,parameter :: ldeg45=.false. 

       integer lpar           
       integer jpnl           


       real(kind=double), dimension(ipar_fastj) :: xgrd 
       real(kind=double), dimension(jpar) :: ygrd           
       real(kind=double), dimension(jpar) :: ydgrd          
       real(kind=double), dimension(kmaxd+1) :: etaa    
       real(kind=double), dimension(kmaxd+1) :: etab    
       real(kind=double) ::  tau_fastj          
       integer month_fastj        
       integer iday_fastj         
       real(kind=double), dimension(ipar_fastj,jpar) ::   P , SA
       real(kind=double), dimension(ipar_fastj,jpar,kmaxd+1) :: T, OD
       real(kind=double) my_ozone(kmaxd)       
       real tauaer_550
       integer iozone
       integer nslat               
       integer nslon               
       save :: nslat, nslon
       integer nspint           
       parameter ( nspint = 4 ) 
       real, dimension (nspint),save :: wavmid 
       real, dimension (nspint, kmaxd+1) :: sizeaer,extaer,waer,gaer,tauaer
       real, dimension (nspint, kmaxd+1) :: l2,l3,l4,l5,l6,l7
       data wavmid     &
           / 0.30e-4, 0.40e-4, 0.60e-4 ,0.999e-04 /

	
    integer julian_day
    real surface_pressure,surface_albedo,pressure(lpar+2),   &
         temperature(lpar+1)
	real optical_depth(lpar+1)
	real tau1
    real*8 pi_fastj,lat,lon,timej,jvalue(lpar,jppj)
    integer isvode,jsvode

	integer iozone1,i
	real my_ozone1(lpar+1)

	real tauaer_550_1
	real sza_dum
	
	integer ientryno_fastj
	save ientryno_fastj
        data ientryno_fastj / 0 /

	


      nslat = 1
      nslon = 1
      pi_fastj=3.141592653589793D0







        do i=1,lpar
           pj(i)=pressure(i)
           T(nslon,nslat,i)=temperature(i)
           OD(nslon,nslat,i)=optical_depth(i)
        enddo
        pj(lpar+1) = pressure(i)
     

	SA(nslon,nslat)=surface_albedo

	iozone=iozone1
	do i=1,iozone1
       my_ozone(i)=my_ozone1(i)
	enddo

	tau_fastj=tau1 
	iday_fastj=julian_day	

	tauaer_550=tauaer_550_1

      month_fastj=int(dble(iday_fastj)*12.d0/365.d0)+1    
      xgrd(nslon)=lon*pi_fastj/180.d0
      ygrd(nslat)=lat*pi_fastj/180.d0
      ydgrd(nslat)=lat


	if (ientryno_fastj .eq. 0) then
            call inphot2
            ientryno_fastj = 1
        end if    


	timej=0.0 
    call photoj(isvode,jsvode,jvalue,timej,nslat,nslon,iozone,tauaer_550, &
        my_ozone,p,t,od,sa,lpar,jpnl, &
        xgrd,ygrd,tau_fastj,month_fastj,iday_fastj,ydgrd, &
        sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)
	sza_dum=SZA

	return
      end subroutine fastj






      subroutine inphot2

















        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data

	IMPLICIT NONE


        integer, parameter :: iprint = 0
        integer, parameter :: single = 4        

       integer,parameter :: ipar_fastj=1,jpar=1

       logical,parameter :: ldeg45=.false. 
       integer lpar           
       integer jpnl           
       real(kind=double), dimension(ipar_fastj) :: xgrd 
       real(kind=double), dimension(jpar) :: ygrd           
       real(kind=double), dimension(jpar) :: ydgrd          
       real(kind=double), dimension(kmaxd+1) :: etaa    
       real(kind=double), dimension(kmaxd+1) :: etab    
       real(kind=double) ::  tau_fastj          
       integer month_fastj        
       integer iday_fastj         








	call rd_tjpl2








      return
      end subroutine inphot2


      subroutine photoj(isvode,jsvode,zpj,timej,nslat,nslon,iozone,tauaer_550_1, &
        my_ozone,p,t,od,sa,lpar,jpnl,xgrd,ygrd,tau_fastj,month_fastj,iday_fastj, &
        ydgrd,sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)















        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	
	IMPLICIT NONE


       integer, parameter :: iprint = 0
       integer, parameter :: single = 4        

       integer,parameter :: ipar_fastj=1,jpar=1

       logical,parameter :: ldeg45=.false. 
       integer lpar           
       integer jpnl           
       real(kind=double), dimension(ipar_fastj) :: xgrd 
       real(kind=double), dimension(jpar) :: ygrd           
       real(kind=double), dimension(jpar) :: ydgrd          
       real(kind=double), dimension(kmaxd+1) :: etaa    
       real(kind=double), dimension(kmaxd+1) :: etab    
       real(kind=double) ::  tau_fastj          
       integer month_fastj        
       integer iday_fastj         
       real(kind=double), dimension(ipar_fastj,jpar) ::   P , SA
       real(kind=double), dimension(ipar_fastj,jpar,kmaxd+1) :: T, OD
       real(kind=double) my_ozone(kmaxd)       
       real tauaer_550_1
       integer iozone
       integer nslat               
       integer nslon               
       integer nspint           
       parameter ( nspint = 4 ) 
       real, dimension (nspint),save :: wavmid 
       real, dimension (nspint, kmaxd+1) :: sizeaer,extaer,waer,gaer,tauaer
       real, dimension (nspint, kmaxd+1) :: l2,l3,l4,l5,l6,l7
       data wavmid     &
           / 0.30e-4, 0.40e-4, 0.60e-4 ,0.999e-04 /

      real*8 zpj(lpar,jppj),timej,solf
	real*8 pi_fastj

	integer i,j
    integer isvode,jsvode



      do i=1,jpnl
        do j=1,jppj
          zj(i,j)=0.D0
          zpj(i,j)=0.D0
        enddo
      enddo


      CALL SOLAR2(timej,nslat,nslon, &
      xgrd,ygrd,tau_fastj,month_fastj,iday_fastj)
      if(SZA.gt.szamax) go to 10


      CALL SET_PROF(isvode,jsvode,nslat,nslon,iozone,tauaer_550_1, &
        my_ozone,p,t,od,sa,lpar,jpnl,month_fastj,ydgrd)


       if(iprint.ne.0)CALL PRTATM(3,nslat,nslon,tau_fastj,month_fastj,ydgrd)  



      CALL JVALUE(isvode,jsvode,lpar,jpnl, &
        ydgrd,sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)









	pi_fastj=3.1415926536d0
       solf=1.d0-(0.034d0*cos(dble(iday_fastj-186)*2.d0   &
             *pi_fastj/365.d0))
	if(iprint.ne.0)then

          write(*,'('' solf = '', f10.5)')solf
	endif


      CALL JRATET(solf,nslat,nslon,p,t,od,sa,lpar,jpnl)




      do i=1,jpnl
        do j=1,jppj
          zpj(i,j)= zj(i,j)
        enddo
      enddo



  10  if((.not.ldeg45.and.nslon.eq.37.and.nslat.eq.36).or.   &
               (ldeg45.and.nslon.eq.19.and.nslat.eq.18)) then
        i=min(jppj,8)

      endif

      return


      end subroutine photoj           


      subroutine set_prof(isvode,jsvode,nslat,nslon,iozone,tauaer_550, &
        my_ozone,p,t,od,sa,lpar,jpnl,month_fastj,ydgrd)






















      
      USE module_data_mosaic_other, only : kmaxd
      USE module_fastj_data
      
      IMPLICIT NONE


      integer, parameter :: iprint = 0
      integer, parameter :: single = 4        

      integer,parameter :: ipar_fastj=1,jpar=1

      logical,parameter :: ldeg45=.false. 
      integer lpar           
      integer jpnl           
      real(kind=double), dimension(ipar_fastj) :: xgrd 
      real(kind=double), dimension(jpar) :: ygrd           
      real(kind=double), dimension(jpar) :: ydgrd          
      real(kind=double), dimension(kmaxd+1) :: etaa    
      real(kind=double), dimension(kmaxd+1) :: etab    
      real(kind=double) ::  tau_fastj          
      integer month_fastj        
      integer iday_fastj         
      real(kind=double), dimension(ipar_fastj,jpar) ::   P , SA
      real(kind=double), dimension(ipar_fastj,jpar,kmaxd+1) :: T, OD
      real(kind=double) my_ozone(kmaxd)       
      real tauaer_550
      integer iozone
      integer nslat               
      integer nslon               

      real*8  pstd(52),oref2(51),tref2(51),bref2(51)
      real*8  odcol(lpar),dlogp,f0,t0,b0,pb,pc,xc,masfac,scaleh
      real vis, aerd1, aerd2

      integer  i, k, l, m
      integer isvode,jsvode

       pj(NB+1) = 0.d0 


      call CLDSRF(isvode,jsvode,odcol,nslat,nslon,p,t,od,sa,lpar,jpnl)


      masfac=100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)






      pstd(1) = max(pj(1),1000.d0)
      pstd(2) = 1000.d0*10.d0**(-1.d0/16.d0)
      dlogp = 10.d0**(-2.d0/16.d0)
      do i=3,51
        pstd(i) = pstd(i-1)*dlogp
      enddo
      pstd(52) = 0.d0


      m = max(1,min(12,month_fastj))
      l = max(1,min(18,(int(ydgrd(nslat))+99)/10))


      do i=1,51
        oref2(i)=oref(i,l,m)
        tref2(i)=tref(i,l,m)
        bref2(i)=bref(i)
      enddo





      do i = 1,NB
        F0 = 0.d0
        T0 = 0.d0
        B0 = 0.d0
        do k = 1,51
          PC = min(pj(i),pstd(k))
          PB = max(pj(i+1),pstd(k+1))
          if(PC.gt.PB) then
            XC = (PC-PB)/(pj(i)-pj(i+1))
            F0 = F0 + oref2(k)*XC
            T0 = T0 + tref2(k)*XC
            B0 = B0 + bref2(k)*XC
          endif
        enddo
        TJ(i) = T0
        DO3(i)= F0*1.d-6
        DBC(i)= B0
      enddo





        do i=1,lpar 
        if(i.le.iozone)DO3(i) = my_ozone(i) 




	TJ(i)=T(nslon,nslat,i)
        enddo
	if(lpar+1.le.iozone)then
        DO3(lpar+1) = my_ozone(lpar+1)  
	endif








      z(1) = 0.d0
      do i=1,lpar
        scaleh=1.3806d-19*masfac*TJ(i)
        z(i+1) = z(i)-(log(pj(i+1)/pj(i))*scaleh)
      enddo





      do i=1,lpar


	vis=23.0
	call aeroden(z(i)/100000.,vis,aerd1) 
	call aeroden(z(i+1)/100000.,vis,aerd2)

	AER(1,i)=(z(i+1)-z(i))/100000.*(aerd1+aerd2)/2./4287.55*tauaer_550


        if(T(nslon,nslat,I).gt.233.d0) then
          AER(2,i) = odcol(i)
          AER(3,i) = 0.d0
        else
          AER(2,i) = 0.d0
          AER(3,i) = odcol(i)
        endif
      enddo
      do k=1,MX
        AER(k,lpar+1) = 0.d0 
      enddo

	AER(1,lpar+1)=2.0*AER(1,lpar) 


      do i=1,NB
        DM(i)  = (PJ(i)-PJ(i+1))*masfac
        DO3(i) = DO3(i)*DM(i)
      enddo

      end subroutine set_prof



      SUBROUTINE CLDSRF(isvode,jsvode,odcol,nslat,nslon,p,t,od,sa, &
         lpar,jpnl)










        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	
	IMPLICIT NONE


      integer, parameter :: iprint = 0
      integer, parameter :: single = 4        

      integer,parameter :: ipar_fastj=1,jpar=1

      logical,parameter :: ldeg45=.false. 
      integer lpar           
      integer jpnl           
      real(kind=double), dimension(ipar_fastj) :: xgrd 
      real(kind=double), dimension(jpar) :: ygrd           
      real(kind=double), dimension(jpar) :: ydgrd          
      real(kind=double), dimension(kmaxd+1) :: etaa    
      real(kind=double), dimension(kmaxd+1) :: etab    
      real(kind=double) ::  tau_fastj          
      integer month_fastj        
      integer iday_fastj         
      integer nslat               
      integer nslon               
      real(kind=double), dimension(ipar_fastj,jpar) ::   P , SA
      real(kind=double), dimension(ipar_fastj,jpar,kmaxd+1) :: T, OD

      integer i, j, k
      integer isvode, jsvode
      real*8  odcol(lpar), odsum, odmax, odtot


      nlbatm = 1


      RFLECT = dble(SA(nslon,nslat))
      RFLECT = max(0.d0,min(1.d0,RFLECT))


      do k=1,MX
        do i=1,NB
          AER(k,i) = 0.d0
        enddo
      enddo


      odmax = 200.d0
      odsum =   0.d0
      do i=1,lpar
        odcol(i) = dble(OD(nslon,nslat,i))
        odsum    = odsum + odcol(i)

      enddo
      if(odsum.gt.odmax) then
        odsum = odmax/odsum
        do i=1,lpar
          odcol(i) = odcol(i)*odsum
        enddo
        odsum = odmax
      endif

      odtot=0.d0
      jadsub(nb)=0
      jadsub(nb-1)=0
      do i=nb-1,1,-1
        k=2*i
        jadsub(k)=0
        jadsub(k-1)=0
        odtot=odtot+odcol(i)
        if(odcol(i).gt.0.d0.and.dtausub.gt.0.d0) then
          if(odtot.le.dtausub) then
            jadsub(k)=1
            jadsub(k-1)=1

          elseif(odtot.gt.dtausub) then
            jadsub(k)=1
            jadsub(k-1)=0
            do j=1,2*(i-1)
              jadsub(j)=0
            enddo

            go to 20
          endif
        endif
      enddo
 20   continue

      return
      end SUBROUTINE CLDSRF       


      subroutine solar2(timej,nslat,nslon, &
      xgrd,ygrd,tau_fastj,month_fastj,iday_fastj)







        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	IMPLICIT NONE


      integer, parameter :: iprint = 0
      integer, parameter :: single = 4        

      integer,parameter :: ipar_fastj=1,jpar=1

      logical,parameter :: ldeg45=.false. 
      integer lpar           
      integer jpnl           
      real(kind=double), dimension(ipar_fastj) :: xgrd 
      real(kind=double), dimension(jpar) :: ygrd           
      real(kind=double), dimension(jpar) :: ydgrd          
      real(kind=double), dimension(kmaxd+1) :: etaa    
      real(kind=double), dimension(kmaxd+1) :: etab    
      real(kind=double) ::  tau_fastj          
      integer month_fastj        
      integer iday_fastj         
      integer nslat               
      integer nslon               

      real*8 pi_fastj, pi180, loct, timej
      real*8 sindec, soldek, cosdec, sinlat, sollat, coslat, cosz

      pi_fastj=3.141592653589793D0
      pi180=pi_fastj/180.d0
      sindec=0.3978d0*sin(0.9863d0*(dble(iday_fastj)-80.d0)*pi180)
      soldek=asin(sindec)
      cosdec=cos(soldek)
      sinlat=sin(ygrd(nslat))
      sollat=asin(sinlat)
      coslat=cos(sollat)

      loct = (((tau_fastj+timej)*15.d0)-180.d0)*pi180 + xgrd(nslon)
      cosz = cosdec*coslat*cos(loct) + sindec*sinlat
      sza  = acos(cosz)/pi180
      U0 = cos(SZA*pi180)

      return
      end subroutine solar2       




      SUBROUTINE JRATET(SOLF,nslat,nslon,p,t,od,sa,lpar,jpnl)














        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	
	IMPLICIT NONE


      integer, parameter :: iprint = 0
      integer, parameter :: single = 4        

      integer,parameter :: ipar_fastj=1,jpar=1

      logical,parameter :: ldeg45=.false. 
      integer lpar           
      integer jpnl           
      real(kind=double), dimension(ipar_fastj) :: xgrd 
      real(kind=double), dimension(jpar) :: ygrd           
      real(kind=double), dimension(jpar) :: ydgrd          
      real(kind=double), dimension(kmaxd+1) :: etaa    
      real(kind=double), dimension(kmaxd+1) :: etab    
      real(kind=double) ::  tau_fastj          
      integer month_fastj        
      integer iday_fastj         
      integer nslat               
      integer nslon               
      real(kind=double), dimension(ipar_fastj,jpar) ::   P , SA
      real(kind=double), dimension(ipar_fastj,jpar,kmaxd+1) :: T, OD

      integer i, j, k
      real*8 qo2tot, qo3tot, qo31d, qo33p, qqqt

      real*8  solf, tfact

      do I=1,jpnl
       VALJ(1) = 0.d0
       VALJ(2) = 0.d0
       VALJ(3) = 0.d0
       do K=NW1,NW2                       
         QO2TOT= xseco2(K,dble(T(nslon,nslat,I)))
         VALJ(1) = VALJ(1) + QO2TOT*FFF(K,I)
         QO3TOT= xseco3(K,dble(T(nslon,nslat,I)))
         QO31D = xsec1d(K,dble(T(nslon,nslat,I)))*QO3TOT
         QO33P = QO3TOT - QO31D
         VALJ(2) = VALJ(2) + QO33P*FFF(K,I)
         VALJ(3) = VALJ(3) + QO31D*FFF(K,I)
       enddo

       do J=4,NJVAL
         VALJ(J) = 0.d0
         TFACT = 0.d0
         if(TQQ(2,J).gt.TQQ(1,J)) TFACT = max(0.d0,min(1.d0,   &
              (T(nslon,nslat,I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))
         do K=NW1,NW2
           QQQT = QQQ(K,1,J-3) + (QQQ(K,2,J-3) - QQQ(K,1,J-3))*TFACT
           VALJ(J) = VALJ(J) + QQQT*FFF(K,I)





         enddo
       enddo
       do j=1,jppj
         zj(i,j)=VALJ(jind(j))*jfacta(j)*SOLF
       enddo

       do j=1,nhz
         zj(i,hzind(j))=hztoa(j)*fhz(i)*SOLF
       enddo
      enddo
      return
      end SUBROUTINE JRATET




      SUBROUTINE PRTATM(N,nslat,nslon,tau_fastj,month_fastj,ydgrd)







        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	


      integer, parameter :: iprint = 0
      integer, parameter :: single = 4        

      integer,parameter :: ipar_fastj=1,jpar=1

      logical,parameter :: ldeg45=.false. 
      integer lpar           
      integer jpnl           
      real(kind=double), dimension(ipar_fastj) :: xgrd 
      real(kind=double), dimension(jpar) :: ygrd           
      real(kind=double), dimension(jpar) :: ydgrd          
      real(kind=double), dimension(kmaxd+1) :: etaa    
      real(kind=double), dimension(kmaxd+1) :: etab    
      real(kind=double) ::  tau_fastj          
      integer month_fastj        
      integer iday_fastj         
      integer nslat               
      integer nslon               

      integer n, i, k, l, m
      real*8 COLO3(NB),COLO2(NB),COLAX(MX,NB),ZKM,ZSTAR,PJC
      real*8 climat(9),masfac,dlogp
      if(N.eq.0) return

      COLO3(NB) = DO3(NB)
      COLO2(NB) = DM(NB)*0.20948d0
      do K=1,MX
        COLAX(K,NB) = AER(K,NB)
      enddo
      do I=NB-1,1,-1
        COLO3(i) = COLO3(i+1)+DO3(i)
        COLO2(i) = COLO2(i+1)+DM(i)*0.20948d0
        do K=1,MX
          COLAX(k,i) = COLAX(k,i+1)+AER(k,i)
        enddo
      enddo
      write(*,1200) ' Tau=',tau_fastj,'  SZA=',sza
      write(*,1200) ' O3-column(DU)=',COLO3(1)/2.687d16,   &
                    '  column aerosol @1000nm=',(COLAX(K,1),K=1,MX)

      if(N.gt.1) then
        write(*,1000) (' AER-X ','col-AER',k=1,mx)
        do I=NB,1,-1
          PJC = PJ(I)
          ZKM =1.d-5*Z(I)
          ZSTAR = 16.d0*DLOG10(1013.d0/PJC)
          write(*,1100) I,ZKM,ZSTAR,DM(I),DO3(I),1.d6*DO3(I)/DM(I),   &
               TJ(I),PJC,COLO3(I),COLO2(I),(AER(K,I),COLAX(K,I),K=1,MX)
        enddo
      endif


      if(N.gt.2) then
        do i=1,9
          climat(i)=0.d0
        enddo
        m = max(1,min(12,month_fastj))
        l = max(1,min(18,(int(ydgrd(nslat))+99)/10))
        masfac=100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)
        write(*,*) 'Specified Climatology'
        write(*,1000)
        do i=51,1,-1
          dlogp=10.d0**(-1.d0/16.d0)
          PJC = 1000.d0*dlogp**(2*i-2)
          climat(1) = 16.d0*DLOG10(1000.D0/PJC)
          climat(2) = climat(1)
          climat(3) = PJC*(1.d0/dlogp-dlogp)*masfac
          if(i.eq.1) climat(3)=PJC*(1.d0-dlogp)*masfac
          climat(4)=climat(3)*oref(i,l,m)*1.d-6
          climat(5)=oref(i,l,m)
          climat(6)=tref(i,l,m)
          climat(7)=PJC
          climat(8)=climat(8)+climat(4)
          climat(9)=climat(9)+climat(3)*0.20948d0
          write(*,1100) I,(climat(k),k=1,9)
        enddo
        write(*,1200) ' O3-column(DU)=',climat(8)/2.687d16
      endif
      return
 1000 format(5X,'Zkm',3X,'Z*',8X,'M',8X,'O3',6X,'f-O3',5X,'T',7X,'P',6x,   &
          'col-O3',3X,'col-O2',2X,10(a7,2x))
 1100 format(1X,I2,0P,2F6.2,1P,2E10.3,0P,F7.3,F8.2,F10.4,1P,10E9.2)
 1200 format(A,F8.1,A,10(1pE10.3))
      end SUBROUTINE PRTATM   

      SUBROUTINE JVALUE(isvode,jsvode,lpar,jpnl, &
        ydgrd,sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)















        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	USE module_peg_util, only:  peg_message


        integer, parameter :: iprint = 0
        integer, parameter :: single = 4        

       integer,parameter :: ipar_fastj=1,jpar=1

       logical,parameter :: ldeg45=.false. 
       integer lpar           
       integer jpnl           
      real(kind=double), dimension(ipar_fastj) :: xgrd 
      real(kind=double), dimension(jpar) :: ygrd           
      real(kind=double), dimension(jpar) :: ydgrd          
      real(kind=double), dimension(kmaxd+1) :: etaa    
      real(kind=double), dimension(kmaxd+1) :: etab    
      real(kind=double) ::  tau_fastj          
      integer month_fastj        
      integer iday_fastj         
      real(kind=double), dimension(ipar_fastj,jpar) ::   P , SA
      real(kind=double), dimension(ipar_fastj,jpar,kmaxd+1) :: T, OD
      integer nspint           
      parameter ( nspint = 4 ) 
      real, dimension (nspint),save :: wavmid 
      real, dimension (nspint, kmaxd+1) :: sizeaer,extaer,waer,gaer,tauaer
      real, dimension (nspint, kmaxd+1) :: l2,l3,l4,l5,l6,l7
      data wavmid     &
        / 0.30e-4, 0.40e-4, 0.60e-4 ,0.999e-04 /

      integer j, k

      real*8  wave
      real*8  AVGF(lpar),XQO3(NB),XQO2(NB)



    integer isvode,jsvode
	character*80 msg

      do J=1,jpnl
        do K=NW1,NW2
          FFF(K,J) = 0.d0
        enddo
        FHZ(J) = 0.d0
      enddo


      if(SZA.gt.szamax) GOTO 99


      CALL SPHERE(lpar)


      do K=NW1,NW2
        WAVE = WL(K)
        do J=1,NB
          XQO3(J) = xseco3(K,dble(TJ(J)))
        enddo
        do J=1,NB
          XQO2(J) = xseco2(K,dble(TJ(J)))
        enddo

        CALL OPMIE(isvode,jsvode,K,WAVE,XQO2,XQO3,AVGF,lpar,jpnl, &
          ydgrd,sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)

        do J=1,jpnl
          FFF(K,J) = FFF(K,J) + FL(K)*AVGF(J)

          if ( FFF(K,J) .lt. 0) then         
	     write( msg, '(a,2i4,e14.6)' )   &
                  'FASTJ neg actinic flux ' //   &
                  'k j FFF(K,J) ', k, j, fff(k,j)
             call peg_message( lunerr, msg )         
          end if

        enddo
      enddo


      if(NHZ.gt.0) then
        K=NW2+1
        WAVE = 204.d0
        do J=1,NB
          XQO3(J) = HZO3
          XQO2(J) = HZO2
        enddo
        CALL OPMIE(isvode,jsvode,K,WAVE,XQO2,XQO3,AVGF,lpar,jpnl, &
          ydgrd,sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)
        do J=1,jpnl
          if(z(j).gt.1.d6) FHZ(J)=AVGF(J)
        enddo
      endif

   99 continue
 1000 format('  SZA=',f6.1,' Reflectvty=',f6.3,' OD=',10(1pe10.3))

      return
      end SUBROUTINE JVALUE

      FUNCTION xseco3(K,TTT)




	USE module_fastj_data
	
      integer k

      real*8 ttt, xseco3
      xseco3  =   &
        flint(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))
      return
      end FUNCTION xseco3       

      FUNCTION xsec1d(K,TTT)




	USE module_fastj_data
	
      integer k

      real*8 ttt,  xsec1d
      xsec1d =   &
        flint(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))
      return
      end FUNCTION xsec1d       

      FUNCTION xseco2(K,TTT)




	USE module_fastj_data
	
      integer k

      real*8 ttt,  xseco2
      xseco2 =   &
        flint(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))
      return
      end FUNCTION xseco2       

      REAL*8 FUNCTION flint (TINT,T1,T2,T3,F1,F2,F3)



      real*8 TINT,T1,T2,T3,F1,F2,F3
      IF (TINT .LE. T2)  THEN
        IF (TINT .LE. T1)  THEN
          flint  = F1
        ELSE
          flint = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT .GE. T3)  THEN
          flint  = F3
        ELSE
          flint = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      return
      end FUNCTION flint

      SUBROUTINE SPHERE(lpar)
















	USE module_fastj_data     
      

      integer lpar

      integer i, j, k, ii
      real*8 airmas, gmu, xmu1, xmu2, xl, diff
      REAL*8 Ux,H,RZ(NB),RQ(NB),ZBYR


      AIRMAS(Ux,H) = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0-   &
               0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*H))/   &
                                                   (1.0d0+0.625d0*H)))

      GMU = U0
      RZ(1)=RAD+Z(1)
      ZBYR = ZZHT/RAD
      DO 2 II=2,NB
        RZ(II) = RAD + Z(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2
    2 CONTINUE
      IF (GMU.LT.0.0D0) THEN
        TANHT = RZ(nlbatm)/DSQRT(1.0D0-GMU**2)
      ELSE
        TANHT = RZ(nlbatm)
      ENDIF



      DO 16 J=1,NB
        DO K=1,NB
          AMF(K,J)=0.D0
        ENDDO


        IF (RZ(J).LT.TANHT) GOTO 16

        XMU1=ABS(GMU)
        DO 12 I=J,lpar
          XMU2=DSQRT(1.0D0-RQ(I)*(1.0D0-XMU1**2))
          XL=RZ(I+1)*XMU2-RZ(I)*XMU1
          AMF(I,J)=XL/(RZ(I+1)-RZ(I))
          XMU1=XMU2
   12   CONTINUE

        AMF(NB,J)=AIRMAS(XMU1,ZBYR)


        IF (GMU.GE.0.0D0) GOTO 16
        XMU1=ABS(GMU)

        DO 14 II=J-1,1,-1
          DIFF=RZ(II+1)*DSQRT(1.0D0-XMU1**2)-RZ(II)
          if(II.eq.1) DIFF=max(DIFF,0.d0)   

          IF (DIFF.LT.0.0D0) THEN
            XMU2=DSQRT(1.0D0-(1.0D0-XMU1**2)/RQ(II))
            XL=ABS(RZ(II+1)*XMU1-RZ(II)*XMU2)
            AMF(II,J)=2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1=XMU2

          ELSE
            XL=RZ(II+1)*XMU1*2.0D0


            AMF(II,J)=XL/(RZ(II+1)-RZ(II))
            GOTO 16
          ENDIF
   14   CONTINUE

   16 CONTINUE
      RETURN
      END SUBROUTINE SPHERE


      SUBROUTINE OPMIE(isvode,jsvode,KW,WAVEL,XQO2,XQO3,FMEAN,lpar,jpnl, &
        ydgrd,sizeaer,extaer,waer,gaer,tauaer,l2,l3,l4,l5,l6,l7)









































      
        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	USE module_peg_util, only:  peg_message, peg_error_fatal    


      integer, parameter :: iprint = 0
      integer, parameter :: single = 4        

      integer,parameter :: ipar_fastj=1,jpar=1

      logical,parameter :: ldeg45=.false. 
      integer lpar           
      integer jpnl           
      real(kind=double), dimension(ipar_fastj) :: xgrd 
      real(kind=double), dimension(jpar) :: ygrd           
      real(kind=double), dimension(jpar) :: ydgrd          
      real(kind=double), dimension(kmaxd+1) :: etaa    
      real(kind=double), dimension(kmaxd+1) :: etab    
      real(kind=double) ::  tau_fastj          
      integer month_fastj        
      integer iday_fastj         
      integer nspint           
      parameter ( nspint = 4 ) 
      real, dimension (nspint),save :: wavmid 
      real, dimension (nspint, kmaxd+1) :: sizeaer,extaer,waer,gaer,tauaer
      real, dimension (nspint, kmaxd+1) :: l2,l3,l4,l5,l6,l7
      data wavmid     &
          / 0.30e-4, 0.40e-4, 0.60e-4 ,0.999e-04 /
      INTEGER    NL, N__, M__
      PARAMETER (NL=500, N__=2*NL, M__=4)  
      REAL(kind=double), dimension(M__) :: A,C1,H,V1,WT,EMU
      REAL(kind=double), dimension(M__,M__) :: B,AA,CC,S,W,U1
      REAL(kind=double), dimension(M__,2*M__) :: PM
      REAL(kind=double), dimension(2*M__) :: PM0
      REAL(kind=double), dimension(2*M__,N__) :: POMEGA
      REAL(kind=double), dimension(N__) :: ZTAU, FZ, FJ
      REAL(kind=double), dimension(M__,M__,N__) :: DD
      REAL(kind=double), dimension(M__,N__) :: RR
      REAL(kind=double) :: ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER ND,N,M,MFIT

      integer jndlev(lpar),jaddlv(nc),jaddto(nc+1)
      integer KW,km,i,j,k,l,ix,j1
      integer isvode,jsvode
      real*8 QXMIE(MX),XLAER(MX),SSALB(MX)
      real*8 xlo2,xlo3,xlray,xltau,zk,taudn,tauup,zk2
      real*8 WAVEL,XQO2(NB),XQO3(NB),FMEAN(lpar),POMEGAJ(2*M__,NC+1)
      real*8 DTAUX(NB),PIRAY(NB),PIAER(MX,NB),TTAU(NC+1),FTAU(NC+1)
      real*8 ftaulog,dttau,dpomega(2*M__)
      real*8 ftaulog2,dttau2,dpomega2(2*M__)

	real*8 PIAER_MX1(NB)

	 character*80 msg


                              KM=1
      if( WAVEL .gt. 355.d0 ) KM=2
      if( WAVEL .gt. 500.d0 ) KM=3





	ang=log(QAA(1,MIEDX(1))/QAA(4,MIEDX(1)))/log(300./999.)
       do I=1,MX

         QXMIE(I) = QAA(KM,MIEDX(I))/QAA(4,MIEDX(I))




         if(I.eq.1) QXMIE(I) = (WAVEL/550.0)**ang
         SSALB(I) = SSA(KM,MIEDX(I)) 
      enddo



      do j=1,nc+1
        ttau(j)=0.d0
        ftau(j)=0.d0
      enddo


      J1 = NLBATM
      do J=J1,NB
        XLO3=DO3(J)*XQO3(J)
        XLO2=DM(J)*XQO2(J)*0.20948d0
        XLRAY=DM(J)*QRAYL(KW)


        do I=1,MX



          XLAER(I)=AER(I,J)*QXMIE(I)
        enddo

        DTAUX(J)=XLO3+XLO2+XLRAY
        do I=1,MX
          DTAUX(J)=DTAUX(J)+XLAER(I)


        enddo




	dtaux(j)=dtaux(j)+tauaer(km,j)




        PIRAY(J)=XLRAY/DTAUX(J)
        do I=1,MX
          PIAER(I,J)=SSALB(I)*XLAER(I)/DTAUX(J)
        enddo


        PIAER_MX1(J)=waer(km,j)*tauaer(km,j)/DTAUX(J)

      enddo 



      N = M__
      MFIT = 2*M__
      do j=j1,NB  
        do i=1,MFIT
          pomegaj(i,j) = PIRAY(J)*PAA(I,KM,1) 
          do k=1,MX 
            pomegaj(i,j) = pomegaj(i,j) + PIAER(K,J)*PAA(I,KM,MIEDX(K))
          enddo
        enddo




           pomegaj(1,j) = pomegaj(1,j) + PIAER_MX1(J)*1.0 
           pomegaj(2,j) = pomegaj(2,j) + PIAER_MX1(J)*gaer(KM,j)*3.0 
           pomegaj(3,j) = pomegaj(3,j) + PIAER_MX1(J)*l2(KM,j)
           pomegaj(4,j) = pomegaj(4,j) + PIAER_MX1(J)*l3(KM,j)
           pomegaj(5,j) = pomegaj(5,j) + PIAER_MX1(J)*l4(KM,j)
           pomegaj(6,j) = pomegaj(6,j) + PIAER_MX1(J)*l5(KM,j)
           pomegaj(7,j) = pomegaj(7,j) + PIAER_MX1(J)*l6(KM,j)
           pomegaj(8,j) = pomegaj(7,j) + PIAER_MX1(J)*l7(KM,j)


      enddo


      do J=J1,NB
        if(AMF(J,J).gt.0.0D0) then
          XLTAU=0.0D0
          do I=1,NB
            XLTAU=XLTAU + DTAUX(I)*AMF(I,J)
          enddo
          if(XLTAU.gt.450.d0) then   
            FTAU(j)=0.d0
          else
            FTAU(J)=DEXP(-XLTAU)
          endif
        else
          FTAU(J)=0.0D0
        endif
      enddo
      if(U0.gt.0.D0) then
        ZFLUX = U0*FTAU(J1)*RFLECT/(1.d0+RFLECT)
      else
        ZFLUX = 0.d0
      endif









      J1=2*J1-1
      do j=1,lpar
        jndlev(j)=2*j
      enddo


      TTAU(NC+1)=0.0D0
      do J=NC,J1,-1
        I=(J+1)/2
        TTAU(J)=TTAU(J+1) + 0.5d0*DTAUX(I)
        jaddlv(j)=int(0.5d0*DTAUX(I)/dtaumax)

        if(jadsub(j).gt.0) then
          jadsub(j)=min(jaddlv(j)+1,nint(dtausub))*(nint(dsubdiv)-1)
          jaddlv(j)=jaddlv(j)+jadsub(j)
        endif


      enddo


      FTAU(NC+1)=1.0D0
      do J=NC-1,J1,-2
        I=(J+1)/2
        FTAU(J)=FTAU(I)
      enddo
      do J=NC,J1,-2
        FTAU(J)=sqrt(FTAU(J+1)*FTAU(J-1))
      enddo



      do j=NC,J1,-2
        do i=1,MFIT
          pomegaj(i,j) = pomegaj(i,j/2)
        enddo
      enddo
      do j=J1+2,nc,2
        taudn = ttau(j-1)-ttau(j)
        tauup = ttau(j)-ttau(j+1)
        do i=1,MFIT
          pomegaj(i,j) = (pomegaj(i,j-1)*taudn +   &
                          pomegaj(i,j+1)*tauup) / (taudn+tauup)
        enddo
      enddo

      do i=1,MFIT
        pomegaj(i,J1)   = pomegaj(i,J1+1)
        pomegaj(i,nc+1) = pomegaj(i,nc)
      enddo












      do j=1,nc+1
        jaddto(j)=0
      enddo

      jaddto(J1)=jaddlv(J1)
      do j=J1+1,nc
        jaddto(j)=jaddto(j-1)+jaddlv(j)
      enddo
      if((jaddto(nc)+nc).gt.nl) then

         write ( msg, '(a, 2i6)' )  &
           'FASTJ  Max NL exceeded '  //  &
           'jaddto(nc)+nc NL', jaddto(nc)+nc,NL
         call peg_message( lunerr, msg )
         msg = 'FASTJ subr OPMIE error. Max NL exceeded'
         call peg_error_fatal( lunerr, msg )      


      endif
      do i=1,lpar
        jndlev(i)=jndlev(i)+jaddto(jndlev(i)-1)
      enddo
      jaddto(nc)=jaddlv(nc)
      do j=nc-1,J1,-1
        jaddto(j)=jaddto(j+1)+jaddlv(j)
      enddo
















      do k=1,N__
        do i=1,MFIT
          pomega(i,k) = 0.d0
        enddo
        ztau(k) = 0.d0
        fz(k)   = 0.d0
      enddo


      do j=J1,nc+1
        k = 2*(nc+1-j)+2*jaddto(j)+1
        ztau(k)= ttau(j)
        fz(k)  = ftau(j)
        do i=1,MFIT
          pomega(i,k) = pomegaj(i,j)
        enddo
      enddo



















      do j=nc,J1,-1
          zk = 0.5d0/(1.d0+dble(jaddlv(j)-jadsub(j)))
          dttau = (ttau(j)-ttau(j+1))*zk
          do i=1,MFIT
            dpomega(i) = (pomegaj(i,j)-pomegaj(i,j+1))*zk
          enddo

          if(ftau(j+1).eq.0.d0) then
            ftaulog=0.d0
          else
            ftaulog = ftau(j)/ftau(j+1)
            if(ftaulog.lt.1.d-150) then
              ftaulog=1.0d-05
            else
              ftaulog=exp(log(ftaulog)*zk)
            endif
          endif
          k = 2*(nc-j+jaddto(j)-jaddlv(j))+1   
          l = 0

          if(jadsub(j).ne.0) then
            l=jadsub(j)/nint(dsubdiv-1)
            zk2=1.d0/dsubdiv
            dttau2=dttau*zk2
            ftaulog2=ftaulog**zk2
            do i=1,MFIT
              dpomega2(i)=dpomega(i)*zk2
            enddo
            do ix=1,2*(jadsub(j)+l)
              ztau(k+1) = ztau(k) + dttau2
              fz(k+1) = fz(k)*ftaulog2
              do i=1,MFIT
                pomega(i,k+1) = pomega(i,k) + dpomega2(i)
              enddo
              k = k+1
            enddo
          endif
          l = 2*(jaddlv(j)-jadsub(j)-l)+1


          do ix=1,l
            ztau(k+1) = ztau(k) + dttau
            fz(k+1) = fz(k)*ftaulog
            do i=1,MFIT
              pomega(i,k+1) = pomega(i,k) + dpomega(i)
            enddo
            k = k+1
          enddo
















      enddo


      ND = 2*(NC+jaddto(J1)-J1)  + 3
      if(nd.gt.N__) then      
         write ( msg, '(a, 2i6)' )  &
           'FASTJ  Max N__ exceeded '  //  &
           'ND N__', ND, N__
         call peg_message( lunerr, msg )
         msg = 'FASTJ subr OPMIE error. Max N__ exceeded'
         call peg_error_fatal( lunerr, msg )


      endif



      ZTAU(ND+1) = ZTAU(ND)*1.000005d0
      ZTAU(ND+2) = ZTAU(ND)*1.000010d0
      zk=max(abs(U0),0.01d0)
      zk=dexp(-ZTAU(ND)*5.d-6/zk)
      FZ(ND+1) = FZ(ND)*zk
      FZ(ND+2) = FZ(ND+1)*zk
      do I=1,MFIT
        POMEGA(I,ND+1)   = POMEGA(I,ND)
        POMEGA(I,ND+2)   = POMEGA(I,ND)
      enddo
      ND = ND+2

      ZU0 = U0
      ZREFL = RFLECT


      CALL MIESCT(ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU,ZFLUX, &
        ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)


      l=2*(NC+jaddto(J1))+3
      do j=1,lpar
        k=l-(2*jndlev(j))
        if(k.gt.ND-2) then
          FMEAN(j) = 0.d0
        else
          FMEAN(j) = FJ(k)
        endif
      enddo

      return
 1000 format(1x,i3,3(2x,1pe10.4),1x,i3)
 1300 format(1x,50(i3))
 1500 format(' Too many levels in photolysis code: need ',i3,' but ',a,   &
             ' dimensioned as ',i3)
      END SUBROUTINE OPMIE                          


      subroutine EFOLD (F0, F1, N, F)
































      implicit none
      real*8 F0,F1,F(250)  
      integer N
      integer I
      real*8 A,DX,D,DSQ,DDP1, B(101),R(101)

      if(F0.eq.0.d0) then
        do I=1,N
          F(I)=0.d0
        enddo
        return
      elseif(F1.eq.0.d0) then
        A = DLOG(F0/1.d-250)
      else
        A = DLOG(F0/F1)
      endif

      DX = float(N)/A
      D = DX*DX
      DSQ = D*D
      DDP1 = D+D+1.d0

      B(2) = DDP1
      R(2) = +D*F0
      do I=3,N
        B(I) = DDP1 - DSQ/B(I-1)
        R(I) = +D*R(I-1)/B(I-1)
      enddo
      F(N+1) = F1
      do I=N,2,-1
        F(I) = (R(I) + D*F(I+1))/B(I)
      enddo
      F(1) = F0
      return
      end subroutine EFOLD               


      subroutine CH_PROF(ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU, &
        ZFLUX,ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)




      USE module_peg_util, only:  peg_message
      
      implicit none

      integer, parameter :: single = 4      
      integer, parameter :: double = 8      
      INTEGER    NL, N__, M__
      PARAMETER (NL=350, N__=2*NL, M__=4)
      REAL(kind=double), dimension(M__) :: A,C1,H,V1,WT,EMU
      REAL(kind=double), dimension(M__,M__) :: B,AA,CC,S,W,U1
      REAL(kind=double), dimension(M__,2*M__) :: PM
      REAL(kind=double), dimension(2*M__) :: PM0
      REAL(kind=double), dimension(2*M__,N__) :: POMEGA
      REAL(kind=double), dimension(N__) :: ZTAU, FZ, FJ
      REAL(kind=double), dimension(M__,M__,N__) :: DD
      REAL(kind=double), dimension(M__,N__) :: RR
      REAL(kind=double) :: ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER ND,N,M,MFIT

      integer i,j
      character*80 msg

      do i=1,ND
        if(ztau(i).ne.0.d0) then
          write ( msg, '(a, i3, 2(1x,1pe9.3))' ) 	&
              'FASTJ subr CH_PROF ztau ne 0. check pomega. ' //  &
              'k ztau(k) fz(k) ', i,ztau(i),fz(i)
          call peg_message( lunerr, msg )

        endif
      enddo
      return
 1100 format(1x,a3,4(a9,2x))
 1200 format(1x,i3,11(1x,1pe9.3))
      end subroutine CH_PROF


      SUBROUTINE MIESCT(ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU, &
        ZFLUX,ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)





























 

      
       implicit none    

      integer, parameter :: single = 4      
      integer, parameter :: double = 8      
      INTEGER    NL, N__, M__
      PARAMETER (NL=350, N__=2*NL, M__=4)
      REAL(kind=double), dimension(M__) :: A,C1,H,V1,WT,EMU
      REAL(kind=double), dimension(M__,M__) :: B,AA,CC,S,W,U1
      REAL(kind=double), dimension(M__,2*M__) :: PM
      REAL(kind=double), dimension(2*M__) :: PM0
      REAL(kind=double), dimension(2*M__,N__) :: POMEGA
      REAL(kind=double), dimension(N__) :: ZTAU, FZ, FJ
      REAL(kind=double), dimension(M__,M__,N__) :: DD
      REAL(kind=double), dimension(M__,N__) :: RR
      REAL(kind=double) :: ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER ND,N,M,MFIT

      integer i, id, im
      real*8  cmeq1


      CALL GAUSSP (N,EMU,WT)


      ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M=1
      DO I=1,N
        CALL LEGND0 (EMU(I),PM0,MFIT)
        DO IM=M,MFIT
          PM(I,IM) = PM0(IM)
        ENDDO
      ENDDO

      CMEQ1 = 0.25D0
      CALL LEGND0 (-ZU0,PM0,MFIT)
      DO IM=M,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      ENDDO


      CALL BLKSLV(ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU,ZFLUX, &
        ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)

      DO ID=1,ND,2
        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
      ENDDO

      RETURN
      END SUBROUTINE MIESCT

      SUBROUTINE BLKSLV(ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU, &
        ZFLUX,ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)






	implicit none

      integer, parameter :: single = 4      
      integer, parameter :: double = 8      
      INTEGER    NL, N__, M__
      PARAMETER (NL=350, N__=2*NL, M__=4)
      REAL(kind=double), dimension(M__) :: A,C1,H,V1,WT,EMU
      REAL(kind=double), dimension(M__,M__) :: B,AA,CC,S,W,U1
      REAL(kind=double), dimension(M__,2*M__) :: PM
      REAL(kind=double), dimension(2*M__) :: PM0
      REAL(kind=double), dimension(2*M__,N__) :: POMEGA
      REAL(kind=double), dimension(N__) :: ZTAU, FZ, FJ
      REAL(kind=double), dimension(M__,M__,N__) :: DD
      REAL(kind=double), dimension(M__,N__) :: RR
      REAL(kind=double) :: ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER ND,N,M,MFIT

      integer i, j, k, id
      real*8  thesum

      CALL GEN(1,ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU,ZFLUX, &
        ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)
      CALL MATIN4 (B)
      DO I=1,N
         RR(I,1) = 0.0d0
        DO J=1,N
          THESUM = 0.0d0
         DO K=1,N
          THESUM = THESUM - B(I,K)*CC(K,J)
         ENDDO
         DD(I,J,1) = THESUM
         RR(I,1) = RR(I,1) + B(I,J)*H(J)
        ENDDO
      ENDDO

      DO ID=2,ND-1
        CALL GEN(ID,ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU,ZFLUX, &
          ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)
        DO I=1,N
          DO J=1,N
          B(I,J) = B(I,J) + A(I)*DD(I,J,ID-1)
          ENDDO
          H(I) = H(I) - A(I)*RR(I,ID-1)
        ENDDO
        CALL MATIN4 (B)
        DO I=1,N
          RR(I,ID) = 0.0d0
          DO J=1,N
          RR(I,ID) = RR(I,ID) + B(I,J)*H(J)
          DD(I,J,ID) = - B(I,J)*C1(J)
          ENDDO
        ENDDO
      ENDDO

      CALL GEN(ND,ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU,ZFLUX, &
        ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)
      DO I=1,N
        DO J=1,N
          THESUM = 0.0d0
          DO K=1,N
          THESUM = THESUM + AA(I,K)*DD(K,J,ND-1)
          ENDDO
        B(I,J) = B(I,J) + THESUM
        H(I) = H(I) - AA(I,J)*RR(J,ND-1)
        ENDDO
      ENDDO
      CALL MATIN4 (B)
      DO I=1,N
        RR(I,ND) = 0.0d0
        DO J=1,N
        RR(I,ND) = RR(I,ND) + B(I,J)*H(J)
        ENDDO
      ENDDO

      DO ID=ND-1,1,-1
       DO I=1,N
        DO J=1,N
         RR(I,ID) = RR(I,ID) + DD(I,J,ID)*RR(J,ID+1)
        ENDDO
       ENDDO
      ENDDO

      DO ID=1,ND,2
        FJ(ID) = 0.0d0
       DO I=1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)
       ENDDO
      ENDDO
      DO ID=2,ND,2
        FJ(ID) = 0.0d0
       DO I=1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)*EMU(I)
       ENDDO
      ENDDO





      RETURN
      END SUBROUTINE BLKSLV


      SUBROUTINE CH_FLUX(ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU, &
        ZFLUX,ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)







      	IMPLICIT NONE

      integer, parameter :: single = 4      
      integer, parameter :: double = 8      
      INTEGER    NL, N__, M__
      PARAMETER (NL=350, N__=2*NL, M__=4)
      REAL(kind=double), dimension(M__) :: A,C1,H,V1,WT,EMU
      REAL(kind=double), dimension(M__,M__) :: B,AA,CC,S,W,U1
      REAL(kind=double), dimension(M__,2*M__) :: PM
      REAL(kind=double), dimension(2*M__) :: PM0
      REAL(kind=double), dimension(2*M__,N__) :: POMEGA
      REAL(kind=double), dimension(N__) :: ZTAU, FZ, FJ
      REAL(kind=double), dimension(M__,M__,N__) :: DD
      REAL(kind=double), dimension(M__,N__) :: RR
      REAL(kind=double) :: ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER ND,N,M,MFIT

      integer I,ID
      real*8 FJCHEK(N__),FZMEAN


      DO ID=1,ND,2
        FJCHEK(ID) = 0.0d0
       DO I=1,N
        FJCHEK(ID) = FJCHEK(ID) + RR(I,ID)*WT(I)*EMU(i)
       ENDDO
      ENDDO


      DO ID=2,ND,2
       DO I=1,N
        FJCHEK(ID) = FJ(ID)
       ENDDO
      ENDDO



      WRITE(34,1200)
      DO ID=2,ND,2
        FZMEAN=sqrt(FZ(ID)*FZ(ID-1))

        WRITE(34,1000) ID, ZU0*FZMEAN-2.0*(FJCHEK(id)-FJCHEK(id-1)),   &
                                     2.0*(FJCHEK(id)+FJCHEK(id-1)),   &
                                     2.0*(FJCHEK(id)+FJCHEK(id-1))/   &
                         (ZU0*FZMEAN-2.0*(FJCHEK(id)-FJCHEK(id-1)))
      ENDDO
      RETURN
 1000 FORMAT(1x,i3,1p,2E12.4,1x,0p,f9.4)
 1200 FORMAT(1x,'Lev',3x,'Downward',4x,'Upward',7x,'Ratio')
      END SUBROUTINE CH_FLUX

      SUBROUTINE NOABS(XLO3,XLO2,XLRAY,BCAER,RFLECT)





      IMPLICIT NONE
      real*8 XLO3,XLO2,XLRAY,BCAER,RFLECT
      XLO3=0.d0
      XLO2=0.d0
      XLRAY=XLRAY*1.d-10
      BCAER=0.d0
      RFLECT=1.d0
      RETURN
      END SUBROUTINE NOABS                              

      SUBROUTINE GEN(ID,ND,N,M,MFIT,POMEGA,PM,PM0,FZ,WT,EMU,ZTAU, &
        ZFLUX,ZREFL,FJ,A,C1,H,V1,B,AA,CC,S,W,U1,DD,RR,RADIUS,ZU0)






	IMPLICIT NONE

      integer, parameter :: single = 4      
      integer, parameter :: double = 8      
      INTEGER    NL, N__, M__
      PARAMETER (NL=350, N__=2*NL, M__=4)
      REAL(kind=double), dimension(M__) :: A,C1,H,V1,WT,EMU
      REAL(kind=double), dimension(M__,M__) :: B,AA,CC,S,W,U1
      REAL(kind=double), dimension(M__,2*M__) :: PM
      REAL(kind=double), dimension(2*M__) :: PM0
      REAL(kind=double), dimension(2*M__,N__) :: POMEGA
      REAL(kind=double), dimension(N__) :: ZTAU, FZ, FJ
      REAL(kind=double), dimension(M__,M__,N__) :: DD
      REAL(kind=double), dimension(M__,N__) :: RR
      REAL(kind=double) :: ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER ND,N,M,MFIT

      integer id, id0, id1, im, i, j, k, mstart
      real*8  sum0, sum1, sum2, sum3
      real*8  deltau, d1, d2, surfac

      IF(ID.EQ.1 .OR. ID.EQ.ND) THEN

       ID0 = ID
       ID1 = ID+1
       IF(ID.GE.ND) ID1 = ID-1
       DO 10 I=1,N
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
        DO IM=M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
        DO IM=M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
         H(I) = 0.5d0*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         A(I) = 0.5d0*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
        DO J=1,I
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
         DO IM=M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         DO IM=M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         W(I,J) = - SUM1*WT(J)
         W(J,I) = - SUM1*WT(I)
         U1(I,J) = - SUM3*WT(J)
         U1(J,I) = - SUM3*WT(I)
          SUM0 = 0.5d0*(SUM0 + SUM2)
         B(I,J) = - SUM0*WT(J)
         B(J,I) = - SUM0*WT(I)
        ENDDO
         S(I,I) = S(I,I) + 1.0d0
         W(I,I) = W(I,I) + 1.0d0
         U1(I,I) = U1(I,I) + 1.0d0
         B(I,I) = B(I,I) + 1.0d0
   10  CONTINUE
       DO I=1,N
         SUM0 = 0.0d0
        DO J=1,N
         SUM0 = SUM0 + S(I,J)*A(J)/EMU(J)
        ENDDO
        C1(I) = SUM0
       ENDDO
       DO I=1,N
        DO J=1,N
          SUM0 = 0.0d0
          SUM2 = 0.0d0
         DO K=1,N
          SUM0 = SUM0 + S(J,K)*W(K,I)/EMU(K)
          SUM2 = SUM2 + S(J,K)*U1(K,I)/EMU(K)
         ENDDO
         A(J) = SUM0
         V1(J) = SUM2
        ENDDO
        DO J=1,N
         W(J,I) = A(J)
         U1(J,I) = V1(J)
        ENDDO
       ENDDO
       IF (ID.EQ.1) THEN

        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25d0*DELTAU
        DO I=1,N
          D1 = EMU(I)/DELTAU
          DO J=1,N
           B(I,J) = B(I,J) + D2*W(I,J)
           CC(I,J) = D2*U1(I,J)
          ENDDO
          B(I,I) = B(I,I) + D1
          CC(I,I) = CC(I,I) - D1

          H(I) = H(I) + 2.0d0*D2*C1(I)
          A(I) = 0.0d0
        ENDDO
       ELSE

        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25d0*DELTAU
        SURFAC = 4.0d0*ZREFL/(1.0d0 + ZREFL)
        DO I=1,N
          D1 = EMU(I)/DELTAU
          H(I) = H(I) - 2.0d0*D2*C1(I)
           SUM0 = 0.0d0
          DO J=1,N
           SUM0 = SUM0 + W(I,J)
          ENDDO
           SUM0 = D1 + D2*SUM0
           SUM1 = SURFAC*SUM0
          DO J=1,N
           B(I,J) = B(I,J) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
          ENDDO
          B(I,I) = B(I,I) + D1
          H(I) = H(I) + SUM0*ZFLUX
          DO J=1,N
           AA(I,J) = - D2*U1(I,J)
          ENDDO
           AA(I,I) = AA(I,I) + D1
           C1(I) = 0.0d0
        ENDDO
       ENDIF

      ELSE
        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = M + MOD(ID+1,2)
        DO I=1,N
          A(I) = EMU(I)/DELTAU
          C1(I) = -A(I)
           SUM0 = 0.0d0
          DO IM=MSTART,MFIT,2
           SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM0(IM)
          ENDDO
          H(I) = SUM0*FZ(ID)
          DO J=1,I
            SUM0 = 0.0d0
           DO IM=MSTART,MFIT,2
            SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM(J,IM)
           ENDDO
            B(I,J) =  - SUM0*WT(J)
            B(J,I) =  - SUM0*WT(I)
          ENDDO
          B(I,I) = B(I,I) + 1.0d0
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE GEN    

      SUBROUTINE LEGND0 (X,PL,N)


      IMPLICIT NONE
      INTEGER N,I
      REAL*8 X,PL(N),DEN

        PL(1) = 1.D0
        PL(2) = X
        DO I=3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.D0/DEN) - PL(I-2)*(1.d0-1.D0/DEN)
        ENDDO
      RETURN
      END SUBROUTINE LEGND0         

      SUBROUTINE MATIN4 (A)



      IMPLICIT NONE
      REAL*8 A(4,4)

      A(2,1) = A(2,1)/A(1,1)
      A(2,2) = A(2,2)-A(2,1)*A(1,2)
      A(2,3) = A(2,3)-A(2,1)*A(1,3)
      A(2,4) = A(2,4)-A(2,1)*A(1,4)
      A(3,1) = A(3,1)/A(1,1)
      A(3,2) = (A(3,2)-A(3,1)*A(1,2))/A(2,2)
      A(3,3) = A(3,3)-A(3,1)*A(1,3)-A(3,2)*A(2,3)
      A(3,4) = A(3,4)-A(3,1)*A(1,4)-A(3,2)*A(2,4)
      A(4,1) = A(4,1)/A(1,1)
      A(4,2) = (A(4,2)-A(4,1)*A(1,2))/A(2,2)
      A(4,3) = (A(4,3)-A(4,1)*A(1,3)-A(4,2)*A(2,3))/A(3,3)
      A(4,4) = A(4,4)-A(4,1)*A(1,4)-A(4,2)*A(2,4)-A(4,3)*A(3,4)

      A(4,3) = -A(4,3)
      A(4,2) = -A(4,2)-A(4,3)*A(3,2)
      A(4,1) = -A(4,1)-A(4,2)*A(2,1)-A(4,3)*A(3,1)
      A(3,2) = -A(3,2)
      A(3,1) = -A(3,1)-A(3,2)*A(2,1)
      A(2,1) = -A(2,1)

      A(4,4) = 1.D0/A(4,4)
      A(3,4) = -A(3,4)*A(4,4)/A(3,3)
      A(3,3) = 1.D0/A(3,3)
      A(2,4) = -(A(2,3)*A(3,4)+A(2,4)*A(4,4))/A(2,2)
      A(2,3) = -A(2,3)*A(3,3)/A(2,2)
      A(2,2) = 1.D0/A(2,2)
      A(1,4) = -(A(1,2)*A(2,4)+A(1,3)*A(3,4)+A(1,4)*A(4,4))/A(1,1)
      A(1,3) = -(A(1,2)*A(2,3)+A(1,3)*A(3,3))/A(1,1)
      A(1,2) = -A(1,2)*A(2,2)/A(1,1)
      A(1,1) = 1.D0/A(1,1)

      A(1,1) = A(1,1)+A(1,2)*A(2,1)+A(1,3)*A(3,1)+A(1,4)*A(4,1)
      A(1,2) = A(1,2)+A(1,3)*A(3,2)+A(1,4)*A(4,2)
      A(1,3) = A(1,3)+A(1,4)*A(4,3)
      A(2,1) = A(2,2)*A(2,1)+A(2,3)*A(3,1)+A(2,4)*A(4,1)
      A(2,2) = A(2,2)+A(2,3)*A(3,2)+A(2,4)*A(4,2)
      A(2,3) = A(2,3)+A(2,4)*A(4,3)
      A(3,1) = A(3,3)*A(3,1)+A(3,4)*A(4,1)
      A(3,2) = A(3,3)*A(3,2)+A(3,4)*A(4,2)
      A(3,3) = A(3,3)+A(3,4)*A(4,3)
      A(4,1) = A(4,4)*A(4,1)
      A(4,2) = A(4,4)*A(4,2)
      A(4,3) = A(4,4)*A(4,3)
      RETURN
      END SUBROUTINE MATIN4    

      SUBROUTINE GAUSSP (N,XPT,XWT)



      IMPLICIT NONE
      INTEGER N,I
      REAL*8  XPT(N),XWT(N)
      REAL*8 GPT4(4),GWT4(4)
      DATA GPT4/.06943184420297D0,.33000947820757D0,.66999052179243D0,   &
                .93056815579703D0/
      DATA GWT4/.17392742256873D0,.32607257743127D0,.32607257743127D0,   &
                .17392742256873D0/
      N = 4
      DO I=1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      ENDDO
      RETURN
      END SUBROUTINE GAUSSP            

      subroutine aeroden(zz,v,aerd)

















	implicit none
	integer mz, nz
        parameter (mz=33)
        real v,aerd
        real*8 zz 
        real alt, aden05, aden23, aer05,aer23      
      dimension alt(mz),aden05(mz),aden23(mz)

	
	real z, f, wth
	integer k, kp
      save alt,aden05,aden23,nz
      

      data nz/mz/

      data alt/   &
              0.0,        1.0,        2.0,        3.0,        4.0,   &
              5.0,        6.0,        7.0,        8.0,        9.0,   &
             10.0,       11.0,       12.0,       13.0,       14.0,   &
             15.0,       16.0,       17.0,       18.0,       19.0,   &
             20.0,       21.0,       22.0,       23.0,       24.0,   &
             25.0,       30.0,       35.0,       40.0,       45.0,   &
             50.0,       70.0,      100.0/
      data aden05/   &
        1.378E+04,  5.030E+03,  1.844E+03,  6.731E+02,  2.453E+02,   &
        8.987E+01,  6.337E+01,  5.890E+01,  6.069E+01,  5.818E+01,   &
        5.675E+01,  5.317E+01,  5.585E+01,  5.156E+01,  5.048E+01,   &
        4.744E+01,  4.511E+01,  4.458E+01,  4.314E+01,  3.634E+01,   &
        2.667E+01,  1.933E+01,  1.455E+01,  1.113E+01,  8.826E+00,   &
        7.429E+00,  2.238E+00,  5.890E-01,  1.550E-01,  4.082E-02,   &
        1.078E-02,  5.550E-05,  1.969E-08/
      data aden23/   &
        2.828E+03,  1.244E+03,  5.371E+02,  2.256E+02,  1.192E+02,   &
        8.987E+01,  6.337E+01,  5.890E+01,  6.069E+01,  5.818E+01,   &
        5.675E+01,  5.317E+01,  5.585E+01,  5.156E+01,  5.048E+01,   &
        4.744E+01,  4.511E+01,  4.458E+01,  4.314E+01,  3.634E+01,   &
        2.667E+01,  1.933E+01,  1.455E+01,  1.113E+01,  8.826E+00,   &
        7.429E+00,  2.238E+00,  5.890E-01,  1.550E-01,  4.082E-02,   &
        1.078E-02,  5.550E-05,  1.969E-08/

      z=max(0.,min(100.,real(zz)))
      aerd=0.
      if(z.gt.alt(nz)) return

      call locate(alt,nz,z,k)
      kp=k+1
      f=(z-alt(k))/(alt(kp)-alt(k))

      if(min(aden05(k),aden05(kp),aden23(k),aden23(kp)).le.0.) then
        aer05=aden05(k)*(1.-f)+aden05(kp)*f
        aer23=aden23(k)*(1.-f)+aden23(kp)*f
      else
        aer05=aden05(k)*(aden05(kp)/aden05(k))**f
        aer23=aden23(k)*(aden23(kp)/aden23(k))**f
      endif

      wth=(1./v-1/5.)/(1./23.-1./5.)
      wth=max(0.,min(1.,wth))

      aerd=(1.-wth)*aer05+wth*aer23




      return
	end subroutine aeroden           

      subroutine locate(xx,n,x,j)

























      implicit none
      integer j,n
      real x,xx(n)
      integer jl,jm,ju



      if(x.eq.xx(1)) then
        j=1
        return
      endif
      if(x.eq.xx(n)) then
        j=n-1
        return
      endif
      jl=1
      ju=n
10    if(ju-jl.gt.1) then
        jm=(ju+jl)/2
         if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      end subroutine locate          

      subroutine rd_tjpl2

































        USE module_data_mosaic_other, only : kmaxd
	USE module_fastj_data
	USE module_peg_util, only:  peg_message, peg_error_fatal
	
	IMPLICIT NONE


        integer, parameter :: iprint = 0
        integer, parameter :: single = 4        

        integer,parameter :: ipar_fastj=1,jpar=1

        logical,parameter :: ldeg45=.false. 
        integer lpar           
        integer jpnl           
        real(kind=double), dimension(ipar_fastj) :: xgrd 
        real(kind=double), dimension(jpar) :: ygrd           
        real(kind=double), dimension(jpar) :: ydgrd          
        real(kind=double), dimension(kmaxd+1) :: etaa    
        real(kind=double), dimension(kmaxd+1) :: etab    
        real(kind=double) ::  tau_fastj          
        integer month_fastj        
        integer iday_fastj         

         integer i, j, k
	 character*7  lpdep(3)
	 character*80 msg

         if(NJVAL.gt.NS) then

            write ( msg, '(a, 2i6)' )  &
              'FASTJ  # xsect supplied > max allowed  '  //  &
              'NJVAL NS ', NJVAL, NS
            call peg_message( lunerr, msg )
            msg = 		&
             'FASTJ Setup Error: # xsect supplied > max allowed. Increase NS'
            call peg_error_fatal( lunerr, msg )              


         endif

         if(NAA.gt.NP) then
            write ( msg, '(a, 2i6)' )  &
              'FASTJ  # aerosol/cloud types > NP  '  //  &
              'NAA NP ', NAA ,NP
            call peg_message( lunerr, msg )
            msg = 		&
              'FASTJ Setup Error: Too many phase functions supplied. Increase NP'
            call peg_error_fatal( lunerr, msg )                       


         endif

      do j=1,jppj
        jind(j)=0
      enddo
      do j=1,NJVAL
        jpdep(j)=0
      enddo
      do j=1,nh
        hzind(j)=0
      enddo


      do j=1,NJVAL
        do k=1,jppj
          if(jlabel(k).eq.titlej(1,j)) jind(k)=j


        enddo
        do k=1,npdep
          if(lpdep(k).eq.titlej(1,j)) jpdep(j)=k
        enddo
      enddo
      do k=1,jppj
        if(jfacta(k).eq.0.d0) then 

           write ( msg, '(a, i6)' )  &
             'FASTJ  Not using photolysis reaction ' , k 
           call peg_message( lunerr, msg ) 
        end if           
        if(jind(k).eq.0) then
          if(jfacta(k).eq.0.d0) then
            jind(k)=1
          else
              write ( msg, '(a, i6)' )  &
               'FASTJ  Which J-rate for photolysis reaction ' , k 
              call peg_message( lunerr, msg ) 


              msg = 'FASTJ subr rd_tjpl2 Unknown Jrate. Incorrect FASTJ setup'
              call peg_error_fatal( lunerr, msg )      
          endif
        endif
      enddo

      i=0
      do j=1,nhz
        do k=1,jppj
          if(jlabel(k).eq.hzlab(j)) then
            i=i+1
            hzind(i)=k
            hztoa(i)=hztmp(j)*jfacta(k)
          endif
        enddo
      enddo
      nhz=i
      if(nhz.eq.0) then
        if(iprint.ne.0) then
           write ( msg, '(a)' )  &
            'FASTJ  Not using Herzberg bin '    
           call peg_message( lunerr, msg )        

         end if
      else
        if(iprint.ne.0) then
           write ( msg, '(a)' )  & 
              'FASTJ Using Herzberg bin for: ' 
           call peg_message( lunerr, msg )
           write( msg, '(a,10a7)' )	&
              'FASTJ ', (jlabel(hzind(i)),i=1,nhz)     

        end if 
      endif











         return
         end subroutine rd_tjpl2



end module module_phot_fastj
