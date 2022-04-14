      MODULE tuv_subs

      use params_mod, only : dp

      IMPLICIT none

      private
      public :: tuv_radfld, sundis, calc_zenith

      integer :: nstr = 1        

      CONTAINS

      SUBROUTINE tuv_radfld( nlambda_start, cld_od_opt, cldfrac, nlyr, nwave, &
                             zenith, z, albedo, &
                             aircol, o2col, o3col, so2col, no2col, &
                             tauaer300, tauaer400, tauaer600, tauaer999, &
                             waer300, waer400, waer600, waer999, &
                             gaer300, gaer400, gaer600, gaer999, &
                             dtaer, omaer, gaer, dtcld, omcld, gcld, &
                             has_aer_ra_feedback, &
                             qll, dobsi, o3_xs, no2_xs, o2_xs, &
                             so2_xs, wmin, wc, tlev, srb_o2_xs, radfld, efld, &
                             e_dir, e_dn, e_up, &
                             dir_fld, dwn_fld, up_fld, dt_cld )



  
      use srb,       only : la_srb
      use rad_trans, only : rtlink




      integer, intent(in)  :: nlambda_start
      integer, intent(in)  :: nlyr
      integer, intent(in)  :: nwave
      integer, intent(in)  :: cld_od_opt
      real, intent(in)  :: zenith
      real, intent(in)  :: dobsi
      real, intent(in)  :: wmin
      real, intent(in)  :: z(:)
      real, intent(in)  :: albedo(:)
      real, intent(in)  :: aircol(:)
      real, intent(in)  :: o2col(:)
      real, intent(in)  :: o3col(:)
      real, intent(in)  :: so2col(:)
      real, intent(in)  :: no2col(:)
      real, intent(in)  :: tauaer300(:)
      real, intent(in)  :: tauaer400(:)
      real, intent(in)  :: tauaer600(:)
      real, intent(in)  :: tauaer999(:)
      real, intent(in)  :: waer300(:)
      real, intent(in)  :: waer400(:)
      real, intent(in)  :: waer600(:)
      real, intent(in)  :: waer999(:)
      real, intent(in)  :: gaer300(:)
      real, intent(in)  :: gaer400(:)
      real, intent(in)  :: gaer600(:)
      real, intent(in)  :: gaer999(:)
      real, intent(in)  :: qll(:)
      real, intent(in)  :: wc(:)
      real, intent(in)  :: tlev(:)
      real, intent(in)  :: cldfrac(:)
      real, intent(in)  :: o2_xs(:)
      real, intent(in)  :: so2_xs(:)
      real, intent(in)  :: o3_xs(:,:)
      real, intent(in)  :: no2_xs(:,:)
      real, intent(out) :: srb_o2_xs(:,:)
      real, intent(out) :: radfld(:,:)
      real, intent(out) :: efld(:,:)
      real, intent(inout)  :: dir_fld(:,:), dwn_fld(:,:), up_fld(:,:)
      real, intent(inout)  :: e_dir(:,:), e_dn(:,:), e_up(:,:)
      real, intent(inout)  :: dt_cld(:)
      real, intent(inout)  :: dtaer(:,:), omaer(:,:), gaer(:,:)
      real, intent(inout)  :: dtcld(:,:), omcld(:,:), gcld(:,:)
      logical, intent(in)  :: has_aer_ra_feedback




      integer :: k, n
      integer :: wn
      integer :: n_radlev, n_radlevp1
      integer :: nid(0:nlyr)
      real :: dtrl(nlyr,nwave)
      real :: dto2(nlyr,nwave)
      real :: dto3(nlyr,nwave)
      real :: dtso2(nlyr,nwave)
      real :: dtno2(nlyr,nwave)


      real :: dtsnw(nlyr,nwave)





      real :: omsnw(nlyr,nwave)
      real :: gsnw(nlyr,nwave)

      real :: edir(nlyr+1)
      real :: edn(nlyr+1)
      real :: eup(nlyr+1)
      real :: fdir(nlyr+1)
      real :: fdn(nlyr+1)
      real :: fup(nlyr+1)

      real :: vcol(nlyr)
      real :: scol(nlyr)

      real :: dsdh(0:nlyr,nlyr)

      n_radlev = size( radfld,dim=2 )
      n_radlevp1 = n_radlev + 1

      do wn = 1,nwave
        omcld(:,wn) = 0.
        omaer(:,wn) = 0.
        omsnw(:,wn) = 0.
        gcld(:,wn)  = 0.
        gaer(:,wn)  = 0.
        gsnw(:,wn)  = 0.
        dtcld(:,wn) = 0.
        dtaer(:,wn) = 0.
        dtsnw(:,wn) = 0.
      end do

      call odrl( wc, aircol, dtrl )
      call seto2( o2col, o2_xs, dto2 )
      call odo3( o3col, o3_xs, dto3, dobsi )
      call setso2( so2col, so2_xs, dtso2 )
      call setno2( no2col, no2_xs, dtno2 )



      if( has_aer_ra_feedback ) then
        call setaer( nlambda_start, wc, tauaer300, tauaer400, &
                     tauaer600, tauaer999, waer300, &
                     waer400, waer600, waer999,     &
                     gaer300, gaer400, gaer600,     &
                     gaer999, dtaer, omaer, gaer )
      endif



      call setcld( nlambda_start, cld_od_opt, z, qll, cldfrac, &
                   dtcld, omcld, gcld )
      dt_cld(:n_radlev) = dtcld(2:n_radlevp1,1)

      call sphers( nlyr, z, zenith, dsdh, nid )


      call airmas( nlyr, dsdh, nid, aircol, vcol, scol )
      call la_srb( nlyr, z, tlev, wmin, &
                   vcol, scol, o2_xs, dto2, srb_o2_xs )

      do wn = nlambda_start,nwave
        call rtlink( &
           nstr, nlyr+1, nlyr, nwave, &
           wn, albedo(wn), zenith, &
           dsdh, nid, &
           dtrl,  &
           dto3,  &
           dto2, &
           dtso2, &
           dtno2,  &
           dtcld, omcld, gcld, &
           dtaer, omaer, gaer, &
           dtsnw, omsnw, gsnw, &
           edir, edn, eup, fdir, fdn, fup )








        radfld(wn,1:n_radlev) = fdir(2:n_radlevp1) + fdn(2:n_radlevp1) + fup(2:n_radlevp1)
        efld(1:n_radlev,wn)    = edir(2:n_radlevp1) + edn(2:n_radlevp1) + eup(2:n_radlevp1)
        dir_fld(1:n_radlev,wn) = fdir(2:n_radlevp1)
        dwn_fld(1:n_radlev,wn) = fdn(2:n_radlevp1)
        up_fld(1:n_radlev,wn)  = fup(2:n_radlevp1)
        e_dir(1:n_radlev,wn)   = edir(2:n_radlevp1)
        e_dn(1:n_radlev,wn)    = edn(2:n_radlevp1)
        e_up(1:n_radlev,wn)    = eup(2:n_radlevp1)
      end do

      END SUBROUTINE tuv_radfld

      SUBROUTINE odrl( wc, aircol, dtrl )














      REAL,    intent(in)  :: aircol(:)
      REAL,    intent(in)  :: wc(:)
      REAL,    intent(out) :: dtrl(:,:)




      INTEGER :: nwave, nlyr
      INTEGER :: wn
      REAL    :: srayl, wmicrn, xx 
      
      nwave = size( wc )
      nlyr  = size( aircol )



      DO wn = 1,nwave






        wmicrn =  wc(wn)*1.E-3
        IF( wmicrn <= 0.55 ) THEN
          xx = 3.6772 + 0.389*wmicrn + 0.09426/wmicrn
        ELSE
          xx = 4.04
        ENDIF
        srayl = 4.02e-28/(wmicrn)**xx




        dtrl(:nlyr,wn) = aircol(:nlyr)*srayl
      END DO

      END SUBROUTINE odrl

      SUBROUTINE odo3( o3col, o3xs, dto3, dobsi )

















      REAL, intent(in)    :: dobsi
      REAL, intent(in)    :: o3col(:)
      REAL, intent(in)    :: o3xs(:,:)
      REAL, intent(inout) :: dto3(:,:)

      INTEGER :: nlyr, nwave
      INTEGER :: wn
      REAL    :: dob_at_grnd, scale_fac

      nwave = size(o3xs,dim=1)
      nlyr  = size(o3col)

      if( dobsi == 0. ) then



        DO wn = 1,nwave
          dto3(:nlyr,wn) = o3col(:nlyr) * o3xs(wn,:nlyr)
        END DO
      else



        dob_at_grnd = sum( o3col(:nlyr) )/2.687e16
        scale_fac   = dobsi/dob_at_grnd
        DO wn = 1,nwave
          dto3(:nlyr,wn) = scale_fac * o3col(:nlyr) * o3xs(wn,:nlyr)
        END DO
      endif

      END SUBROUTINE odo3

      SUBROUTINE seto2( o2col, o2xs1, dto2 )














      REAL, intent(in)  :: o2col(:)
      REAL, intent(in)  :: o2xs1(:)
      REAL, intent(out) :: dto2(:,:)




      INTEGER :: nlyr, nwave
      INTEGER :: wn

      nwave = size(o2xs1)
      nlyr  = size(o2col)





      DO wn = 1,nwave
        dto2(:nlyr,wn) = o2col(:nlyr) * o2xs1(wn)
      ENDDO

      END SUBROUTINE seto2

      SUBROUTINE setso2( colso2, so2_xs, dtso2 )

















      REAL,    intent(in)  :: colso2(:)
      REAL,    intent(in)  :: so2_xs(:)
      REAL,    intent(out) :: dtso2(:,:)




      integer :: nwave, nlyr
      integer :: wn

      nwave = size( so2_xs )
      nlyr  = size( colso2 )

      DO wn = 1,nwave
        dtso2(:nlyr,wn) = colso2(:nlyr)*so2_xs(wn)
      END DO

      END SUBROUTINE setso2

      SUBROUTINE setno2( colno2, no2_xs, dtno2 )

















      REAL, intent(in)    :: colno2(:)
      REAL, intent(in)    :: no2_xs(:,:)
      REAL, intent(inout) :: dtno2(:,:)

      INTEGER :: nlyr, nwave
      INTEGER :: wn

      nwave = size(no2_xs,dim=1)
      nlyr  = size(colno2)

      DO wn = 1,nwave
        dtno2(:nlyr,wn) = colno2(:nlyr) * no2_xs(wn,:nlyr)
      END DO

      END SUBROUTINE setno2

      subroutine setaer( nlambda_start, wc, tauaer300, tauaer400, &
                         tauaer600, tauaer999,               &
                         waer300, waer400, waer600, waer999, &
                         gaer300, gaer400, gaer600, gaer999, &
                         dtaer, omaer, gaer )



















      integer, intent(in) :: nlambda_start
      real, intent(in)  :: wc(:)
      real, intent(in)  :: tauaer300(:), tauaer400(:),    &
                           tauaer600(:), tauaer999(:)
      real, intent(in)  :: waer300(:), waer400(:),        &
                           waer600(:), waer999(:)
      real, intent(in)  :: gaer300(:), gaer400(:),        &
                           gaer600(:), gaer999(:)
      real, intent(out) :: dtaer(:,:), omaer(:,:), gaer(:,:)




      real, parameter :: thresh = 1.e-9
      integer     :: k, wn, nlyr, nwave
      real        :: ang, slope, wfac

      nlyr =  size(dtaer,dim=1)
      nwave = size(dtaer,dim=2)

wave_loop: &
      do wn = nlambda_start,nwave
        wfac = wc(wn)*1.e-3 - .6
        do k = 1,nlyr-1



          if( tauaer300(k) > thresh .and. tauaer999(k) > thresh ) then
            ang = log(tauaer300(k)/tauaer999(k))/log(0.999/0.3)
            dtaer(k,wn) = tauaer400(k)*(0.4/(wc(wn)*1.e-3))**ang



            slope = 5.*(waer600(k) - waer400(k))
            omaer(k,wn) = slope*wfac + waer600(k)
            omaer(k,wn) = max( .4,min( 1.,omaer(k,wn) ) )



            slope = 5.*(gaer600(k) - gaer400(k))
            gaer(k,wn) = slope*wfac + gaer600(k)
            gaer(k,wn) = max( .5,min( 1.,gaer(k,wn) ) )
          endif
        end do
      end do wave_loop

      end subroutine setaer

      subroutine setcld( nlambda_start, cld_od_opt, z, xlwc, cldfrac, &
                         dtcld, omcld, gcld )























      integer, intent(in) :: nlambda_start
      integer, intent(in) :: cld_od_opt
      real, intent(in)  :: z(:)
      real, intent(in)  :: xlwc(:)
      real, intent(in)  :: cldfrac(:)
      real, intent(inout) :: dtcld(:,:)
      real, intent(inout) :: omcld(:,:)
      real, intent(inout) :: gcld(:,:)



      real, parameter :: km2m = 1.e3        
      real, parameter :: wden = 1.e6        
      real, parameter :: re = 10.0 * 1.e-6  
      real, parameter :: fac = 1./(wden*re)

      integer  :: astat
      integer  :: wn
      integer  :: nlyr, nwave
      real, allocatable :: wrk(:), layer_cldfrac(:)

      nlyr  = size(dtcld,dim=1)
      nwave = size(dtcld,dim=2)

      allocate( wrk(nlyr),layer_cldfrac(nlyr),stat=astat )
      if( astat /= 0 ) then
        call wrf_error_fatal3("<stdin>",529,&
'setcld: failed to allocate wrk' )
      endif




      wrk(1:nlyr-1) = (z(2:nlyr) - z(1:nlyr-1))*km2m   
      wrk(1:nlyr-1) = 1.5 * .5*(xlwc(1:nlyr-1) + xlwc(2:nlyr))*wrk(1:nlyr-1)*fac
      wrk(1:nlyr-1) = max( wrk(1:nlyr-1),0. )
      if( cld_od_opt == 2 ) then
        layer_cldfrac(1:nlyr-1) = .5*(cldfrac(1:nlyr-1) + cldfrac(2:nlyr))
        wrk(1:nlyr-1) = wrk(1:nlyr-1)*layer_cldfrac(1:nlyr-1)*sqrt( layer_cldfrac(1:nlyr-1) )
      endif




      if( any( wrk(1:nlyr-1) > 0. ) ) then
        do wn = nlambda_start,nwave
          dtcld(1:nlyr-1,wn) = wrk(1:nlyr-1)
          omcld(1:nlyr-1,wn) = .9999
          gcld (1:nlyr-1,wn) = .85
        end do
      endif

      if( allocated( wrk ) ) then
        deallocate( wrk )
      endif
      if( allocated( layer_cldfrac ) ) then
        deallocate( layer_cldfrac )
      endif

      end subroutine setcld

      REAL FUNCTION sundis( julday )












      implicit none




      integer, intent(in) :: julday




      real(dp), parameter :: pi = 3.1415926_dp

      real(dp) :: dayn, thet0
      real(dp) :: sinth, costh, sin2th, cos2th




      dayn = real(julday - 1,kind=dp) + .5_dp



      thet0 = 2._dp*pi*dayn/365._dp







      sinth  = sin( thet0 )
      costh  = cos( thet0 )
      sin2th = 2._dp*sinth*costh
      cos2th = costh*costh - sinth*sinth
      sundis  = real( 1.000110_dp + .034221_dp*costh + .001280_dp*sinth &
                                  + .000719_dp*cos2th + .000077_dp*sin2th )

      END FUNCTION sundis

      subroutine calc_zenith( lat, long, julday, gmt, zenith, &
                              its, ite, jts, jte, &
                              ims, ime, jms, jme )













        integer, intent(in)  :: julday
        integer, intent(in)  :: its,ite
        integer, intent(in)  :: jts,jte
        integer, intent(in)  :: ims,ime
        integer, intent(in)  :: jms,jme
        real(dp), intent(in) :: gmt
        real, intent(in)     :: lat(ims:ime,jms:jme)
        real, intent(in)     :: long(ims:ime,jms:jme)
        real, intent(out)    :: zenith(ims:ime,jms:jme)




      real(dp), parameter :: d2r = 3.1415926_dp/180.0_dp
      real(dp), parameter :: r2d = 1.0_dp/d2r

      integer  :: i, j
      real(dp) :: caz, csz, cw, d, ec, epsi, eqt, eyt, feqt, feqt1, &
          feqt2, feqt3, feqt4, feqt5, feqt6, feqt7, lbgmt, lzgmt, ml, pepsi, &
          pi, ra, raz, rdecl, reqt, rlt, rml, rra, ssw, sw, tab, w, wr, &
          yt, zpt, zr

      d = real(julday,dp) + gmt/24.0_dp



      ml  = 279.2801988_dp + d*(.9856473354_dp + 2.267E-13_dp*d)
      rml = ml*d2r






      w = 282.4932328_dp + d*(4.70684E-5_dp + 3.39E-13_dp*d)
      wr = w*d2r
      ec   = 1.6720041E-2_dp - d*(1.1444E-9_dp + 9.4E-17_dp*d)
      epsi = 23.44266511_dp - d*(3.5626E-7_dp + 1.23E-15_dp*d)
      pepsi = epsi*d2r
      yt = (tan(pepsi/2.0_dp))**2
      cw = cos(wr)
      sw = sin(wr)
      ssw = sin(2.0_dp*wr)
      eyt = 2._dp*ec*yt
      feqt1 = -sin(rml)*cw*(eyt + 2._dp*ec)
      feqt2 = cos(rml)*sw*(2._dp*ec - eyt)
      feqt3 = sin(2._dp*rml)*(yt - (5._dp*ec**2/4._dp)*(cw**2 - sw**2))
      feqt4 = cos(2._dp*rml)*(5._dp*ec**2*ssw/4._dp)
      feqt5 = sin(3._dp*rml)*(eyt*cw)
      feqt6 = -cos(3._dp*rml)*(eyt*sw)
      feqt7 = -sin(4._dp*rml)*(.5_dp*yt**2)
      feqt  = feqt1 + feqt2 + feqt3 + feqt4 + feqt5 + feqt6 + feqt7
      eqt   = feqt*13751.0_dp




      reqt = eqt/240._dp



      ra  = ml - reqt
      rra = ra*d2r



      tab = 0.43360_dp*sin(rra)
      rdecl = atan(tab)
      do j = jts,jte
        do i = its,ite



          lbgmt = 12.0_dp - eqt/3600._dp + real(long(i,j),dp)*24._dp/360._dp
          lzgmt = 15.0_dp*(gmt - lbgmt)
          zpt   = lzgmt*d2r
          rlt   = real(lat(i,j),dp)*d2r
          csz   = sin(rlt)*sin(rdecl) + cos(rlt)*cos(rdecl)*cos(zpt)
          csz   = min( 1._dp,csz )
          zr    = acos(csz)
          zenith(i,j) = real( zr/d2r,4 )
        end do
      end do

      end subroutine calc_zenith

      SUBROUTINE sphers( nlyr, z, zen, dsdh, nid )



































      use params_mod, only : pi, radius

      IMPLICIT NONE




      INTEGER, intent(in) :: nlyr
      REAL,    intent(in) :: zen
      REAL,    intent(in) :: z(:)

      INTEGER, intent(inout) :: nid(0:nlyr)
      REAL,    intent(inout) :: dsdh(0:nlyr,nlyr)





      REAL, parameter ::  dr = pi/180.

      INTEGER :: j, jm1, k
      INTEGER :: id
      REAL    :: re
      REAL    :: zd(0:nlyr)
      REAL    :: ds_dh(1:nlyr)
      REAL(8) :: zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm

      zenrad = REAL( zen*dr,8 )



      re = radius + z(1)





      zd(0:nlyr) = z(nlyr+1:1:-1) - z(1)




      nid(0:nlyr) = 0




layer_loop : &
      DO k = 0, nlyr
        ds_dh(:) = 0.
        rpsinz = real(re + zd(k),8) * SIN(zenrad)

        IF( zen <= 90.0 .or. rpsinz >= real(re,8) ) THEN



          id = k 
          IF( zen > 90.0 ) THEN
            DO j = 1,nlyr
              IF( rpsinz < real(zd(j-1) + re,8) .AND. &
                  rpsinz >= real(zd(j) + re,8) ) then
                id = j
              ENDIF
            END DO
          END IF
 
          DO j = 1, id
            jm1 = j - 1

            IF( j /= id .or. k /= id .or. zen <= 90.0 ) then
              sm = 1.0_8
            ELSE
              sm = -1.0_8
            ENDIF
            rj   = real(re + zd(jm1),8)
            rjp1 = real(re + zd(j),8)
            dhj  = zd(jm1) - zd(j)
 
            ga = max( rj*rj - rpsinz*rpsinz,0.0_8 )
            gb = max( rjp1*rjp1 - rpsinz*rpsinz,0.0_8 )

            IF( id > k .AND. j == id ) THEN
              dsj = SQRT( ga )
            ELSE
              dsj = SQRT( ga ) - sm*SQRT( gb )
            END IF
            ds_dh(j) = real( dsj/dhj,4 )
          END DO
          nid(k) = id
        ELSE
          nid(k) = -1
        ENDIF
        dsdh(k,:) = ds_dh(:)

      END DO layer_loop

      END SUBROUTINE sphers

      SUBROUTINE airmas( nlyr, dsdh, nid, cz, vcol, scol )



















      use params_mod, only : largest

      IMPLICIT NONE




      INTEGER, intent(in) :: nlyr
      INTEGER, intent(in) :: nid(0:nlyr)
      REAL,    intent(in) :: dsdh(0:nlyr,nlyr)
      REAL,    intent(in) :: cz(nlyr)

      REAL, intent(inout) :: vcol(nlyr), scol(nlyr)




      INTEGER :: lyr, j, nlev, nlevi
      REAL    :: sum, vsum




      nlev = nlyr + 1
      vsum = 0.
      DO lyr = 1, nlyr
        nlevi = nlev - lyr
        vsum = vsum + cz(nlevi)
        vcol(nlevi) = vsum
        sum = 0.
        IF( nid(lyr) < 0 ) THEN
          sum = largest
        ELSE



          DO j = 1, MIN(nid(lyr), lyr)
            sum = sum + cz(nlev-j)*dsdh(lyr,j)
          END DO



           DO j = MIN(nid(lyr),lyr)+1, nid(lyr)
             sum = sum + 2.*cz(nlev-j)*dsdh(lyr,j)
           END DO
        ENDIF
        scol(nlevi) = sum
      END DO

      END SUBROUTINE airmas

      END MODULE tuv_subs
