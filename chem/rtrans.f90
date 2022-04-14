







      MODULE rad_trans

      IMPLICIT NONE

      private
      public :: rtlink

      CONTAINS

      SUBROUTINE rtlink( &
           nstr, nlev, nlyr, nwave, &
           iw, albedo, zen, &
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

      use params_mod, only : largest, pi




      INTEGER, intent(in) :: nstr
      INTEGER, intent(in) :: nlev, nlyr
      INTEGER, intent(in) :: nwave, iw
      REAL, intent(in)    :: albedo
      REAL, intent(in)    :: zen
      REAL, intent(in)    :: dsdh(0:nlyr,nlyr)
      INTEGER, intent(in) :: nid(0:nlyr)

      REAL, intent(in)    :: dtrl(nlyr,nwave)
      REAL, intent(in)    :: dto3(nlyr,nwave), dto2(nlyr,nwave)
      REAL, intent(in)    :: dtso2(nlyr,nwave), dtno2(nlyr,nwave)
      REAL, intent(in)    :: dtcld(nlyr,nwave), omcld(nlyr,nwave), gcld(nlyr,nwave)
      REAL, intent(in)    :: dtaer(nlyr,nwave), omaer(nlyr,nwave), gaer(nlyr,nwave)
      REAL, intent(in)    :: dtsnw(nlyr,nwave), omsnw(nlyr,nwave), gsnw(nlyr,nwave)

      REAL, intent(out)   :: edir(nlev), edn(nlev), eup(nlev)
      REAL, intent(out)   :: fdir(nlev), fdn(nlev), fup(nlev)






      REAL, PARAMETER :: dr = pi/180.

      INTEGER :: k, kk
      REAL    :: dtabs,dtsct,dscld,dsaer,dssnw,dacld,daaer,dasnw
      REAL    :: dt(nlyr), om(nlyr), g(nlyr)





      LOGICAL, parameter :: delta = .true.
      REAL ediri(nlev), edni(nlev), eupi(nlev)
      REAL fdiri(nlev), fdni(nlev), fupi(nlev)




      fdir(1:nlev) = 0.
      fup(1:nlev)  = 0.
      fdn(1:nlev)  = 0.
      edir(1:nlev) = 0.
      eup(1:nlev)  = 0.
      edn(1:nlev)  = 0.

      DO k = 1, nlyr
        dscld = dtcld(k,iw)*omcld(k,iw)
        dacld = dtcld(k,iw)*(1.-omcld(k,iw))

        dsaer = dtaer(k,iw)*omaer(k,iw)
        daaer = dtaer(k,iw)*(1.-omaer(k,iw))

        dssnw = dtsnw(k,iw)*omsnw(k,iw)
        dasnw = dtsnw(k,iw)*(1.-omsnw(k,iw))

        dtsct = dtrl(k,iw) + dscld + dsaer + dssnw
        dtabs = dtso2(k,iw) + dto2(k,iw) + dto3(k,iw) & 
              + dtno2(k,iw) + dacld + daaer + dasnw

        dtabs = max( dtabs,1./largest )
        dtsct = max( dtsct,1./largest )




        kk = nlyr - k + 1
        dt(kk) = dtsct + dtabs
        g(kk)  = (gcld(k,iw)*dscld + gsnw(k,iw)*dssnw + gaer(k,iw)*dsaer)/dtsct
        IF( dtsct /= 1./largest ) then
          om(kk) = dtsct/(dtsct + dtabs)
        ELSE
          om(kk) = 1./largest
        ENDIF
      END DO   


      CALL ps2str( nlyr, nlev, zen, albedo, &
                   dt, om, g, &
                   dsdh, nid, delta, &
                   fdiri, fupi, fdni, ediri, eupi, edni)




      fdir(1:nlev) = fdiri(nlev:1:-1)
      fup(1:nlev)  = fupi(nlev:1:-1)
      fdn(1:nlev)  = fdni(nlev:1:-1)
      edir(1:nlev) = ediri(nlev:1:-1)
      eup(1:nlev)  = eupi(nlev:1:-1)
      edn(1:nlev)  = edni(nlev:1:-1)

      END SUBROUTINE rtlink

      SUBROUTINE ps2str( &
           nlyr, nlev, zen, rsfc, &
           tauu, omu, gu, &
           dsdh, nid, delta, &
           fdr, fup, fdn, edr, eup, edn)






































      use params_mod, only : pi, precis, largest




      INTEGER, intent(in) :: nlyr, nlev
      REAL, intent(in)    :: zen, rsfc
      REAL, intent(in)    :: tauu(nlyr), omu(nlyr), gu(nlyr)
      REAL, intent(in)    :: dsdh(0:nlyr,nlyr)
      INTEGER, intent(in) :: nid(0:nlyr)
      LOGICAL, intent(in) :: delta

      REAL, intent(out) :: fup(nlev), fdn(nlev), fdr(nlev)
      REAL, intent(out) :: eup(nlev), edn(nlev), edr(nlev)




      REAL, PARAMETER    :: eps = 1.E-3



      REAL, PARAMETER    :: pifs = 1.
      REAL, PARAMETER    :: fdn0 = 0.

      REAL :: mu, sum
      REAL :: tausla(0:nlyr)
      REAL :: tauc(0:nlyr)
      REAL :: mu2(0:nlyr)




      INTEGER :: row, nlyrm1
      REAL :: lam(nlyr), taun(nlyr), bgam(nlyr)
      REAL :: e1(nlyr), e2(nlyr), e3(nlyr), e4(nlyr)
      REAL :: cup(nlyr), cdn(nlyr), cuptn(nlyr), cdntn(nlyr)
      REAL :: mu1(nlyr)
      REAL :: a(2*nlyr), b(2*nlyr), d(2*nlyr), e(2*nlyr), y(2*nlyr)

      REAL :: f, g, om, tmpg
      REAL :: gam1, gam2, gam3, gam4
      REAL :: gi(nlyr), omi(nlyr)










      INTEGER :: mrows, mrowsm1, mrowsm2
      REAL    :: expon, expon0, expon1, divisr, tmp, up, dn
      REAL    :: ssfc
      REAL    :: wrk, wrk1

      INTEGER :: i, im1, j, k












      mu = COS( zen*pi/180. )









      tauc(0:nlyr)   = 0.
      tausla(0:nlyr) = 0.
      mu2(0:nlyr)    = 1./SQRT(largest)

      IF( .NOT. delta ) THEN
        gi(1:nlyr)   = gu(1:nlyr)
        omi(1:nlyr)  = omu(1:nlyr)
        taun(1:nlyr) = tauu(1:nlyr)
      ELSE 





        DO k = 1, nlyr
          f       = gu(k)*gu(k)
          wrk     = 1. - f
          wrk1    = 1. - omu(k)*f
          gi(k)   = (gu(k) - f)/wrk
          omi(k)  = wrk*omu(k)/wrk1
          taun(k) = wrk1*tauu(k)
        ENDDO
      END IF






      IF( zen > 90.0 ) THEN
        IF(nid(0) < 0) THEN
          tausla(0) = largest
        ELSE
          sum = 0.0
          DO j = 1, nid(0)
            sum = sum + 2.*taun(j)*dsdh(0,j)
          END DO
          tausla(0) = sum 
        END IF
      END IF
  
layer_loop : &
      DO i = 1, nlyr
        im1 = i - 1
        g  = gi(i)
        om = omi(i)
        tauc(i) = tauc(im1) + taun(i)




        tmpg = MIN( abs(g),1. - precis )
        g    = SIGN( tmpg,g )
        om   = MIN( om,1. - precis )




        IF(nid(i) < 0) THEN
          tausla(i) = largest
        ELSE
          sum = 0.0
          DO j = 1, MIN(nid(i),i)
             sum = sum + taun(j)*dsdh(i,j)
          ENDDO
          DO j = MIN(nid(i),i)+1,nid(i)
             sum = sum + 2.*taun(j)*dsdh(i,j)
          ENDDO
          tausla(i) = sum 
          IF(tausla(i) == tausla(im1)) THEN
            mu2(i) = SQRT(largest)
          ELSE
            mu2(i) = (tauc(i) - tauc(im1))/(tausla(i) - tausla(im1))
            mu2(i) = SIGN( MAX( ABS(mu2(i)),1./SQRT(largest) ), mu2(i) )
          END IF
        END IF






        gam1 =  .25*(7. - om*(4. + 3.*g))
        gam2 = -.25*(1. - om*(4. - 3.*g))
        gam3 = .25*(2. - 3.*g*mu)
        gam4 = 1. - gam3
        mu1(i) = 0.5







        lam(i) = sqrt(gam1*gam1 - gam2*gam2)

        IF( gam2 /= 0.) THEN
          bgam(i) = (gam1 - lam(i))/gam2
        ELSE
          bgam(i) = 0.
        ENDIF

        expon = EXP( -lam(i)*taun(i) )




        e1(i) = 1. + bgam(i)*expon
        e2(i) = 1. - bgam(i)*expon
        e3(i) = bgam(i) + expon
        e4(i) = bgam(i) - expon







        expon0 = EXP( -tausla(im1) )
        expon1 = EXP( -tausla(i) )
          
        divisr = lam(i)*lam(i) - 1./(mu2(i)*mu2(i))
        tmp    = MAX( eps,abs(divisr) )
        divisr = SIGN( tmp,divisr )

        up = om*pifs*((gam1 - 1./mu2(i))*gam3 + gam4*gam2)/divisr
        dn = om*pifs*((gam1 + 1./mu2(i))*gam4 + gam2*gam3)/divisr
         




        cup(i) = up*expon0
        cdn(i) = dn*expon0
        cuptn(i) = up*expon1
        cdntn(i) = dn*expon1
      end do layer_loop





      ssfc = rsfc*mu*EXP( -tausla(nlyr) )*pifs




      mrows   = 2*nlyr     
      mrowsm1 = mrows - 1
      mrowsm2 = mrows - 2
      nlyrm1  = nlyr - 1
      




      a(1) = 0.
      b(1) = e1(1)
      d(1) = -e2(1)
      e(1) = fdn0 - cdn(1)




      a(3:mrowsm1:2) = e2(1:nlyrm1)*e3(1:nlyrm1) - e4(1:nlyrm1)*e1(1:nlyrm1)
      b(3:mrowsm1:2) = e1(1:nlyrm1)*e1(2:nlyr) - e3(1:nlyrm1)*e3(2:nlyr)
      d(3:mrowsm1:2) = e3(1:nlyrm1)*e4(2:nlyr) - e1(1:nlyrm1)*e2(2:nlyr)
      e(3:mrowsm1:2) = e3(1:nlyrm1)*(cup(2:nlyr) - cuptn(1:nlyrm1))  &
                     + e1(1:nlyrm1)*(cdntn(1:nlyrm1) - cdn(2:nlyr))




      a(2:mrowsm2:2) = e2(2:nlyr)*e1(1:nlyrm1) - e3(1:nlyrm1)*e4(2:nlyr)
      b(2:mrowsm2:2) = e2(1:nlyrm1)*e2(2:nlyr) - e4(1:nlyrm1)*e4(2:nlyr)
      d(2:mrowsm2:2) = e1(2:nlyr)*e4(2:nlyr) - e2(2:nlyr)*e3(2:nlyr)
      e(2:mrowsm2:2) = (cup(2:nlyr) - cuptn(1:nlyrm1))*e2(2:nlyr) & 
                     - (cdn(2:nlyr) - cdntn(1:nlyrm1))*e4(2:nlyr)




      a(mrows) = e1(nlyr) - rsfc*e3(nlyr)
      b(mrows) = e2(nlyr) - rsfc*e4(nlyr)
      d(mrows) = 0.
      e(mrows) = ssfc - cuptn(nlyr) + rsfc*cdntn(nlyr)




      CALL tridag( a, b, d, e, y, mrows )






      fdr(1) = EXP( -tausla(0) )
      edr(1) = mu * fdr(1)
      edn(1) = fdn0
      eup(1) = y(1)*e3(1) - y(2)*e4(1) + cup(1)
      fdn(1) = edn(1)/mu1(1)
      fup(1) = eup(1)/mu1(1)

      fdr(2:nlev) = EXP( -tausla(1:nlyr) )
      edr(2:nlev) = mu *fdr(2:nlev)
      edn(2:nlev) = y(1:mrowsm1:2)*e3(1:nlyr) + y(2:mrows:2)*e4(1:nlyr) + cdntn(1:nlyr)
      eup(2:nlev) = y(1:mrowsm1:2)*e1(1:nlyr) + y(2:mrows:2)*e2(1:nlyr) + cuptn(1:nlyr)
      fdn(2:nlev) = edn(2:nlev)/mu1(1:nlyr)
      fup(2:nlev) = eup(2:nlev)/mu1(1:nlyr)

      END SUBROUTINE ps2str

      SUBROUTINE tridag( a, b, c, r, u, n)







      INTEGER, intent(in) :: n
      REAL,    intent(in) :: a(n), b(n), c(n), r(n)
      REAL, intent(out)   :: u(n)




      INTEGER :: j, jm1, jp1
      REAL    :: bet
      REAL    :: gam(n)
      CHARACTER(len=64) :: err_msg

      IF (b(1) == 0.) then
        call wrf_error_fatal3("<stdin>",490,&
'tridag: zero pivot @ n == 1' )
      ENDIF
      bet  = b(1)
      u(1) = r(1)/bet
      DO j = 2, n   
         jm1 = j - 1
         gam(j) = c(jm1)/bet
         bet    = b(j) - a(j)*gam(j)
         IF (bet == 0.) then
           write(err_msg,'(''tridag: zero pivot @ n = '',i4)') j
           call wrf_error_fatal3("<stdin>",501,&
trim(err_msg) )
         ENDIF
         u(j) = (r(j) - a(j)*u(jm1))/bet
      END DO

      DO j = n - 1, 1, -1  
         jp1 = j + 1
         u(j) = u(j) - gam(jp1)*u(jp1)
      END DO

      END SUBROUTINE tridag

      END MODULE rad_trans
