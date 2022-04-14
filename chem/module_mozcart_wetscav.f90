
MODULE module_mozcart_wetscav

   USE module_configure
   USE module_state_description

IMPLICIT NONE

private
public :: wetscav_mozcart_init
public :: wetscav_mozcart

save

   real, parameter :: zero = 0.
   real, parameter :: one  = 1.
   real, parameter :: four = 4.
   real(8), parameter :: oner8  = 1._8
   real(8), parameter :: fourr8 = 4._8
   real, parameter :: adj_factor = one + 10.*epsilon( one )
   real, parameter :: TICE = 273.
   real, parameter :: TMIX = 258.
   integer, parameter :: idiag = 0
   integer, parameter :: jdiag = 0
   integer, parameter :: kdiag = 18

   integer :: hetcnt
   integer :: hno3_ndx = 0
   integer, allocatable :: wrf2tab(:)
   REAL, allocatable    :: mol_wght(:)
   logical, allocatable :: ice_uptake(:)

   type wet_scav
     character(len=12) :: name
     integer :: p_ndx
     real    :: heff(6)
     real    :: molecw
     logical :: ice_uptake
     real    :: reteff
   end type wet_scav

   type(wet_scav), allocatable :: wet_scav_tab(:)

   LOGICAL, EXTERNAL  :: wrf_dm_on_monitor

CONTAINS

   subroutine wetscav_mozcart_init( id, numgas, config_flags )



   use module_scalar_tables, only : chem_dname_table
   use module_chem_utilities, only : upcase
   use module_HLawConst, only : nHLC, HLC




   integer, intent(in)                           :: id
   integer, intent(in)                           :: numgas
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags




   integer :: m, m1, m2
   integer :: astat
   character(len=12)  :: wrf_spc_name
   character(len=64)  :: HL_tbl_name
   character(len=256) :: message

is_allocated : &
   if( .not. allocated( wet_scav_tab ) ) then
     call wrf_message( ' ' )
     call wrf_message( 'wetscav_mozcart_init: species names' )
     do m = param_first_scalar,numgas,6
       write(message,'(6a12)') chem_dname_table(id,m:min(m+5,numgas))
       call wrf_message( trim(message) )
     enddo
     call wrf_message( ' ' )




     hetcnt = 0
     do m = 1,nHLC
       HL_tbl_name = HLC(m)%name
       call upcase( HL_tbl_name )
       do m1 = param_first_scalar,numgas
         wrf_spc_name = chem_dname_table(id,m1)
         call upcase( wrf_spc_name )
         if( trim(HL_tbl_name) == trim(wrf_spc_name) ) then
           hetcnt = hetcnt + 1 
           exit
         endif
       end do
     end do

     if( hetcnt > 0 ) then
       allocate( wet_scav_tab(hetcnt),stat=astat )
       if( astat /= 0 ) then
         call wrf_error_fatal3("<stdin>",102,&
"mozcart_wetscav_init: failed to allocate wet_scav_tab")
       endif
     else
       call wrf_message( ' ' )
       call wrf_message( 'wetscav_mozcart_init: There are no wet scavenging mozcart species' )
       call wrf_message( ' ' )
       return
     endif




     m2 = 0
     do m = 1,nHLC
       HL_tbl_name = HLC(m)%name
       call upcase( HL_tbl_name )
       do m1 = param_first_scalar,numgas
         wrf_spc_name = chem_dname_table(id,m1)
         call upcase( wrf_spc_name )
         if( trim(HL_tbl_name) == trim(wrf_spc_name) ) then
           m2 = m2 + 1
           wet_scav_tab(m2)%name  = chem_dname_table(id,m1)
           wet_scav_tab(m2)%p_ndx = m1
           wet_scav_tab(m2)%molecw = HLC(m)%mw
           wet_scav_tab(m2)%heff(1:6) = HLC(m)%hcnst(1:6)
           exit
         endif
       end do
     end do

     wet_scav_tab(:)%ice_uptake = .false.
     wet_scav_tab(:)%reteff = 0.
     do m = 1,hetcnt
       if( wet_scav_tab(m)%p_ndx == p_h2o2 ) then
         wet_scav_tab(m)%ice_uptake = .true.
         wet_scav_tab(m)%reteff     = .64
       elseif( wet_scav_tab(m)%p_ndx == p_hno3 ) then
         wet_scav_tab(m)%ice_uptake = .true.
         wet_scav_tab(m)%reteff     = 1.
         hno3_ndx = m
       elseif( wet_scav_tab(m)%p_ndx == p_hcooh ) then
         wet_scav_tab(m)%ice_uptake = .true.
         wet_scav_tab(m)%reteff     = .64
       elseif( wet_scav_tab(m)%p_ndx == p_ch3ooh ) then
         wet_scav_tab(m)%ice_uptake = .true.
         wet_scav_tab(m)%reteff     = .02
       elseif( wet_scav_tab(m)%p_ndx == p_so2 ) then
         wet_scav_tab(m)%ice_uptake = .true.
         wet_scav_tab(m)%reteff     = .02
       elseif( wet_scav_tab(m)%p_ndx == p_hcooh ) then
         wet_scav_tab(m)%ice_uptake = .true.
         wet_scav_tab(m)%reteff     = .68
       endif
     end do

     allocate( wrf2tab(hetcnt),mol_wght(hetcnt),ice_uptake(hetcnt),stat=astat )
     if( astat /= 0 ) then
       call wrf_error_fatal3("<stdin>",160,&
"mozcart_wetscav_init: failed to allocate wrf2tab,mol_wght,ice_uptake")
     endif
     do m = 1,hetcnt
       wrf2tab(m)    = m
       mol_wght(m)   = wet_scav_tab(m)%molecw
       ice_uptake(m) = wet_scav_tab(m)%ice_uptake
     end do

     call wrf_message( 'wetscav_mozcart_init: Wet scavenging mozcart species' )
     do m = 1,hetcnt
       write(message,*) m,wrf2tab(m),trim(wet_scav_tab(wrf2tab(m))%name),mol_wght(m),ice_uptake(m)
       call wrf_message( trim(message) )
     end do

     if( wrf_dm_on_monitor() ) then
       call wrf_debug( 100,' ' )
       write(message,*) 'wetscav_mozcart_init: hetcnt = ',hetcnt
       call wrf_debug( 100, trim(message) )
       write(message,*) 'wetscav_mozcart_init: hno3_ndx = ',hno3_ndx
       call wrf_debug( 100, trim(message) )
       call wrf_debug( 100,' ' )
     endif
   endif is_allocated

   end subroutine wetscav_mozcart_init

   subroutine wetscav_mozcart( id, ktau, dtstep, ktauc, config_flags,                      &
                               dtstepc, t_phy, p8w, t8w, p_phy,                            &
                               chem, rho_phy, cldfra, rainprod, evapprod,                  &
                               qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp,                         &
                               gas_aqfrac, numgas_aqfrac, dz8w, dx, dy,                    &
                               qv, qc, qi, qs,                                             &


                               delta_mass_col,                                             &

                               ids,ide, jds,jde, kds,kde,                                  &
                               ims,ime, jms,jme, kms,kme,                                  &
                               its,ite, jts,jte, kts,kte                                   )

   USE module_configure, only: grid_config_rec_type
   USE module_state_description
   USE module_model_constants, only: g, mwdry




   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc, numgas_aqfrac
   REAL, INTENT(IN   ) :: dtstep, dtstepc
   REAL, INTENT(IN   ) :: dx, dy



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(IN ) ::                                   gas_aqfrac



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
         INTENT(IN   ) ::                                           &
                                                      t_phy,        &
                                                      p_phy,        &
                                                      t8w,          &
                                                      p8w,          &
                                                      dz8w,         &
                                                    rho_phy

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 QV, &
                                                                 QI, &
                                                                 QC, &
                                                                 QS
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 rainprod, &
                                                                 evapprod
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 qv_b4mp, &
                                                                 qc_b4mp, &
                                                                 qi_b4mp, &
                                                                 qs_b4mp

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
         INTENT(INOUT ) ::                                cldfra





   REAL,  DIMENSION( ims:ime , jms:jme, num_chem )         ,           &
         INTENT(OUT ) ::                                delta_mass_col





   real, parameter    :: t0 = 298.


   real, parameter    :: hion = 1.e-5
   real, parameter    :: hion_inv = 1./hion
   real, parameter    :: henry_thres = 1.e4

   integer :: i, j, k, ktem1, m, m1
   integer :: pndx
   integer :: cld_col_cnt
   integer :: precip_col_cnt
   integer :: max_ndx(3)
   integer :: wrk_ndx(3)
   REAL :: area
   REAL :: e298, dhr
   REAL :: percent_cld
   REAL :: percent_precip
   REAL :: max_rls
   REAL :: layer_mass(kts:kte)
   REAL :: delp(kts:kte)
   REAL :: p(kts:kte)
   REAL :: t(kts:kte)
   REAL :: rls(kts:kte)
   REAL :: evaprate(kts:kte)
   REAL :: totprec(kts:kte)
   REAL :: totevap(kts:kte)
   REAL :: wrk(kts:kte)
   REAL :: tfac(kts:kte)
   REAL :: dk1s(kts:kte)
   REAL :: dk2s(kts:kte)
   REAL :: kh(kts:kte)
   REAL :: diff(its:ite,kts:kte,jts:jte)
   REAL :: hno3(kts:kte)
   REAL :: wrk_mass(hetcnt)
   REAL :: trc_mass(kts:kte,hetcnt)
   REAL :: heff(kts:kte,hetcnt)

   logical :: is_hno3
   logical :: tckaqb(hetcnt)


   REAL :: reteff(hetcnt)


   character(len=128) :: message

has_wet_scav : &
   if( hetcnt > 0 ) then



     CALL cal_cldfra3( CLDFRA, qc_b4mp, qi_b4mp, qs_b4mp,         &
                      ids,ide, jds,jde, kds,kde,                  &
                      ims,ime, jms,jme, kms,kme,                  &
                      its,ite, jts,jte, kts,kte                   )




     ktem1 = kte - 1
     area = dx * dy
     diff(:,:,:) = 0.
     cld_col_cnt = 0
     precip_col_cnt = 0


     delta_mass_col(:,:,:) = 0.

     max_rls       = 0.
jloop : &
     do j = jts,jte
iloop : &
       do i = its,ite
           t(kts:kte)      = t_phy(i,kts:kte,j)
           tfac(kts:ktem1) = (t0 - t(kts:ktem1))/(t0*t(kts:ktem1))
           p(kts:kte)      = p_phy(i,kts:kte,j)*.01
           delp(kts:ktem1) = p8w(i,kts:ktem1,j) - p8w(i,kts+1:kte,j)
           layer_mass(kts:ktem1) = area*delp(kts:ktem1)/g
           totprec(kts:ktem1) = rainprod(i,kts:ktem1,j)*layer_mass(kts:ktem1)
           totevap(kts:ktem1) = evapprod(i,kts:ktem1,j)*layer_mass(kts:ktem1)
           rls(kte)      = 0.
           evaprate(kte) = 0.
           do k = ktem1,kts,-1
             rls(k) = max( 0.,totprec(k)-totevap(k)+rls(k+1) )
             evaprate(k) = min( 1.,totevap(k)/(rls(k+1)+1.e-20) )
           end do
column_has_precip : &
         if( any( rls(kts:ktem1) > 0. ) ) then
           if( maxval(rls(kts:ktem1)) >= max_rls ) then
             max_rls = max( max_rls,maxval(rls(kts:ktem1)) )
             max_ndx(3:3) = maxloc(rls(kts:ktem1))
             max_ndx(1:2) = (/ i,j /)
             if( max_ndx(3) /= kts ) then
               write(message,'(''wetscav: max rls not at srf; time,i,j,k = '',4i6)') ktau,max_ndx(:)
               call wrf_debug( 100,trim(message) )
             endif
           endif
           precip_col_cnt = precip_col_cnt + 1
species_loop : &
           do m = 1,hetcnt
             m1 = wrf2tab(m)
             pndx = wet_scav_tab(m1)%p_ndx
             if( pndx == p_hno3 ) then
               hno3(kts:kte)   = chem(i,kts:kte,j,p_hno3)
             endif
             wrk(kts:ktem1)        = 1.e-6*mol_wght(m)*chem(i,kts:ktem1,j,pndx)/mwdry
             trc_mass(kts:ktem1,m) = wrk(kts:ktem1)*layer_mass(kts:ktem1)
             wrk_mass(m)    = sum( trc_mass(kts:ktem1,m) )
             e298 = wet_scav_tab(m1)%heff(1)
             dhr  = wet_scav_tab(m1)%heff(2)
             kh(kts:ktem1) = e298 * exp( dhr*tfac(kts:ktem1) )
             e298 = wet_scav_tab(m1)%heff(3)
             dhr  = wet_scav_tab(m1)%heff(4)
             if( e298 /= 0. ) then
               dk1s(kts:ktem1) = e298*exp( dhr*tfac(kts:ktem1) )
             else
               dk1s(kts:ktem1) = 0.
             endif
             e298 = wet_scav_tab(m1)%heff(5)
             dhr  = wet_scav_tab(m1)%heff(6)
             if( e298 /= 0. ) then
               dk2s(kts:ktem1) = e298*exp( dhr*tfac(kts:ktem1) )
             else
               dk2s(kts:ktem1) = 0.
             endif
             if( pndx /= p_nh3 ) then
               heff(kts:ktem1,m) = kh(kts:ktem1)*(1. + dk1s(kts:ktem1)*hion_inv*(1. + dk2s(kts:ktem1)*hion_inv))
             else



               heff(kts:ktem1,m) = kh(kts:ktem1)*(1. + dk1s(kts:ktem1)*hion/dk2s(kts:ktem1))
             endif
             tckaqb(m) = any( heff(kts:ktem1,m) > henry_thres )

             reteff(m) = wet_scav_tab(m1)%reteff

           end do species_loop




           CALL washout( kte-kts+1, hetcnt, dtstep, trc_mass, layer_mass, &
                         p, dz8w(i,kts:kte,j), rls, qc_b4mp(i,kts:kte,j), qi_b4mp(i,kts:kte,j), &
                         cldfra(i,kts:kte,j), t, evaprate, area, heff, &


                         mol_wght, tckaqb, ice_uptake, i, j, reteff )


species_loop1 : &
           do m = 1,hetcnt
             m1 = wrf2tab(m)
             pndx = wet_scav_tab(m1)%p_ndx
             is_hno3 = pndx == p_hno3




             delta_mass_col(i,j,pndx) = sum( trc_mass(kts:ktem1,m) ) - wrk_mass(m)


             wrk(kts:ktem1) = 1.e6*mwdry*trc_mass(kts:ktem1,m)/mol_wght(m)
             chem(i,kts:ktem1,j,pndx) = wrk(kts:ktem1)/layer_mass(kts:ktem1)
             if( is_hno3 ) then
               diff(i,kts:ktem1,j) = 100.*(chem(i,kts:ktem1,j,p_hno3) - hno3(kts:ktem1))/hno3(kts:ktem1)
             endif
           end do species_loop1
         endif column_has_precip
       end do iloop
     end do jloop

     write(message,'(''washout: max rls @ (i,j,k) '',3i4,'' = '',1pg15.7)') max_ndx(:),max_rls
     call wrf_debug( 100,trim(message) )

     percent_precip = 100.*real(precip_col_cnt)/real((ite-its+1)*(jte-jts+1)) 
     write(*,*) 'wetscav_mozcart: percent columns with precip = ',percent_precip,'%'
   endif has_wet_scav

   end subroutine wetscav_mozcart

   SUBROUTINE cal_cldfra( CLDFRA,QC,QI,                               &
                          ids,ide, jds,jde, kds,kde,                  &
                          ims,ime, jms,jme, kms,kme,                  &
                          its,ite, jts,jte, kts,kte                   )

















   IMPLICIT NONE

   INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                             CLDFRA

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 QI, &
                                                                 QC




   REAL, parameter :: thresh = 1.e-9

   INTEGER :: j,k

   DO j = jts,jte
     DO k = kts,kte-1
       where( (qc(its:ite,k,j) + qi(its:ite,k,j)) > thresh )
         cldfra(its:ite,k,j) = one
       elsewhere
         cldfra(its:ite,k,j) = zero
       endwhere
     ENDDO
   ENDDO

   END SUBROUTINE cal_cldfra

   SUBROUTINE cal_cldfra3( CLDFRA,QC,QI, QS,                          &
                          ids,ide, jds,jde, kds,kde,                  &
                          ims,ime, jms,jme, kms,kme,                  &
                          its,ite, jts,jte, kts,kte                   )

















   IMPLICIT NONE

   INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                             CLDFRA

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 QI, &
                                                                 QC, &
                                                                 QS



   REAL, parameter :: thresh = 1.e-9

   INTEGER :: j

   DO j = jts,jte
     where( (qc(its:ite,kts:kte,j) + qi(its:ite,kts:kte,j) + qs(its:ite,kts:kte,j)) > thresh )
       cldfra(its:ite,kts:kte,j) = one
     elsewhere
       cldfra(its:ite,kts:kte,j) = zero
     endwhere
   ENDDO

   END SUBROUTINE cal_cldfra3

   SUBROUTINE cal_cldfra2( CLDFRA, QV, QC, QI, QS,                     &
                           t_phy, p_phy,                               &
                           ids,ide, jds,jde, kds,kde,                  &
                           ims,ime, jms,jme, kms,kme,                  &
                           its,ite, jts,jte, kts,kte                   )




   INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                             CLDFRA

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 QV, &
                                                                 QI, &
                                                                 QC, &
                                                                 QS, &
                                                              t_phy, &
                                                              p_phy





   integer, parameter  :: idbg = 144, kdbg = 15, jdbg = 147
   REAL    , PARAMETER :: ALPHA0 = 100.
   REAL    , PARAMETER :: GAMMA = 0.49
   REAL    , PARAMETER :: QCLDMIN = 1.E-12
   REAL    , PARAMETER :: PEXP = 0.25
   REAL    , PARAMETER :: RHGRID =1.0
   REAL    , PARAMETER :: SVP1 = 0.61078
   REAL    , PARAMETER :: SVP2 = 17.2693882
   REAL    , PARAMETER :: SVPI2 = 21.8745584
   REAL    , PARAMETER :: SVP3 = 35.86
   REAL    , PARAMETER :: SVPI3 = 7.66
   REAL    , PARAMETER :: SVPT0 = 273.15
   REAL    , PARAMETER :: r_d = 287.
   REAL    , PARAMETER :: r_v = 461.6
   REAL    , PARAMETER :: ep_2 = r_d/r_v

   INTEGER :: i,j,k
   INTEGER :: imax, jmax, kmax
   REAL    :: RHUM, tc, esw, esi, weight, qvsw, qvsi, qvs_weight
   REAL    :: QCLD, DENOM, ARG, SUBSAT, wrk
   REAL    :: relhum_max, wrk_max












































    CLDFRA(its:ite,kts:kte,jts:jte) = 0.
    relhum_max = -100.
    wrk_max    = -10000.
    imax = 0; kmax = 0; jmax = 0

    DO j = jts,jte
      DO k = kts,kte
        DO i = its,ite



          QCLD = QI(i,k,j) + QC(i,k,j) + QS(i,k,j)
has_cloud : &
          IF( QCLD >= QCLDMIN ) THEN
            tc   = t_phy(i,k,j) - SVPT0
            esw  = 1000.0 * SVP1 * EXP( SVP2  * tc / ( t_phy(i,k,j) - SVP3  ) )
            esi  = 1000.0 * SVP1 * EXP( SVPI2 * tc / ( t_phy(i,k,j) - SVPI3 ) )
            QVSW = EP_2 * esw / ( p_phy(i,k,j) - esw )
            QVSI = EP_2 * esi / ( p_phy(i,k,j) - esi )

            weight     = (QI(i,k,j) + QS(i,k,j)) / QCLD
            QVS_WEIGHT = (1. - weight)*QVSW + weight*QVSI
            RHUM       = QV(i,k,j)/QVS_WEIGHT              



            IF( RHUM >= RHGRID )THEN




              CLDFRA(i,k,j) = 1.
            ELSE




              SUBSAT = MAX( 1.E-10,RHGRID*QVS_WEIGHT - QV(i,k,j) )
              DENOM  = SUBSAT**GAMMA
              ARG    = MAX( -6.9,-ALPHA0*QCLD/DENOM )    



              RHUM = MAX( 1.E-10, RHUM )
              wrk  = (RHUM/RHGRID)**PEXP*(1. - EXP( ARG ))
              if( rhum >= relhum_max ) then
                relhum_max = rhum
                imax = i
                kmax = k
                jmax = j
              endif
              IF( wrk >= .01 ) then
                CLDFRA(i,k,j) = wrk
                if( wrk >= wrk_max ) then
                  wrk_max = wrk
                endif
              ENDIF
            ENDIF
          ENDIF has_cloud
        END DO
      END DO
    END DO

   END SUBROUTINE cal_cldfra2

   subroutine WASHOUT( LPAR, NTRACE, DTSCAV, QTTJFL, QM, &
                       POFL, DELZ, RLS, CLWC, CIWC, &
                       CFR, TEM, EVAPRATE, GAREA, HSTAR, &


                       TCMASS, TCKAQB, TCNION, ii, jj, RETEFF )










      



      integer, intent(in)  :: LPAR
      integer, intent(in)  :: NTRACE
      integer, intent(in)  :: ii, jj
      real,  intent(in)    :: DTSCAV
      real,  intent(in)    :: GAREA
      real,  intent(in)    :: QM(LPAR)
      real,  intent(in)    :: POFL(LPAR)
      real,  intent(in)    :: DELZ(LPAR)
      real,  intent(in)    :: RLS(LPAR)
      real,  intent(in)    :: CLWC(LPAR)
      real,  intent(in)    :: CIWC(LPAR)
      real,  intent(in)    :: CFR(LPAR)
      real,  intent(in)    :: TEM(LPAR)
      real,  intent(in)    :: EVAPRATE(LPAR)
      real,  intent(in)    :: TCMASS(NTRACE)
      real,  intent(in)    :: HSTAR(LPAR,NTRACE)
      real,  intent(inout) :: QTTJFL(LPAR,NTRACE)
      logical, intent(in)  :: TCKAQB(NTRACE)
      logical, intent(in)  :: TCNION(NTRACE)

      real,  intent(in)    :: RETEFF(NTRACE)





      integer  :: I, J, L, LE, LM1, N
      integer  :: LWASHTYP, LICETYP

      real :: WRK, RNEW_TST
      real :: CLWX
      real :: RNEW, RPRECIP, DELTARIMEMASS, DELTARIME, RAMPCT
      real :: MASSLOSS
      real :: DOR, DNEW, DEMP, COLEFFSNOW, RHOSNOW
      real :: WEMP, REMP, RRAIN, RWASH
      real :: QTPRECIP, QTRAIN, QTCXA, QTAX

      real :: FAMA, RAMA, DAMA, FCA, RCA, DCA
      real :: FAX, RAX, DAX, FCXA, RCXA, DCXA, FCXB, RCXB, DCXB
      real :: RAXADJ, FAXADJ, RAXADJF
      real :: QTDISCF, QTDISRIME, QTDISCXA
      real :: QTEVAPAXP, QTEVAPAXW, QTEVAPAX
      real :: QTWASHAX
      real :: QTEVAPCXAP, QTEVAPCXAW, QTEVAPCXA
      real :: QTWASHCXA, QTRIMECXA
      real :: QTRAINCXA, QTRAINCXB
      real :: QTTOPCA, QTTOPAA, QTTOPCAX, QTTOPAAX

      real :: AMPCT, AMCLPCT, CLNEWPCT, CLNEWAMPCT, CLOLDPCT, CLOLDAMPCT
      real :: RAXLOC, RCXALOC, RCXBLOC, RCALOC, RAMALOC, RCXPCT

      real :: QTNETLCXA, QTNETLCXB, QTNETLAX, QTNETL
      real :: QTDISSTAR

      real :: CFXX(lpar)
      real :: QTT(lpar)
      real :: QTTNEW(lpar)
      real :: rls_wrk(lpar)
      real :: rnew_wrk(lpar)
      real :: rca_wrk(lpar)
      real :: fca_wrk(lpar)
      real :: rcxa_wrk(lpar)
      real :: fcxa_wrk(lpar)
      real :: rcxb_wrk(lpar)
      real :: fcxb_wrk(lpar)
      real :: rax_wrk(lpar,2)
      real :: fax_wrk(lpar,2)
      real :: rama_wrk(lpar)
      real :: fama_wrk(lpar)
      real :: deltarime_wrk(lpar)
      real :: clwx_wrk(lpar)
      real :: frc(lpar,3)
      real :: rlsog(lpar)

      logical :: is_hno3
      logical :: rls_flag(lpar)
      logical :: rnew_flag(lpar)
      logical :: cf_trigger(lpar)
      logical :: freezing(lpar)

      character(len=132) :: message
      

      real, parameter  :: TFROZ      = 240.
      real, parameter  :: CFMIN      = 1.0
      real, parameter  :: CWMIN      = 1.0e-9
      real, parameter  :: DMIN       = 1.0e-1    
      real, parameter  :: VOLPOW     = 1./3.
      real, parameter  :: RHORAIN    = 1.0e3     
      real, parameter  :: RHOSNOWFIX = 1.0e2     
      real, parameter  :: COLEFFRAIN = 0.7
      real, parameter  :: COLEFFAER  = 0.05




      LE = LPAR - 1   
      rls_flag(1:le) = rls(1:le) > zero
      freezing(1:le) = tem(1:le) < tice
      rlsog(1:le) = rls(1:le)/garea

species_loop : &
     do N = 1,NTRACE
       QTT(:lpar)    = QTTJFL(:lpar,N)
       QTTNEW(:lpar) = QTTJFL(:lpar,N)
       is_hno3 = n == hno3_ndx
       if( is_hno3 ) then
         rca_wrk(:lpar) = zero
         fca_wrk(:lpar) = zero
         rcxa_wrk(:lpar) = zero
         fcxa_wrk(:lpar) = zero
         rcxb_wrk(:lpar) = zero
         fcxb_wrk(:lpar) = zero
         rls_wrk(:lpar) = zero
         rnew_wrk(:lpar) = zero
         cf_trigger(:lpar) = .false.
         clwx_wrk(:lpar) = -9999.
         deltarime_wrk(:lpar) = -9999.
         rax_wrk(:lpar,:) = zero
         fax_wrk(:lpar,:) = zero
       endif




       if( TCKAQB(N) ) then
         LWASHTYP = 1
       else
         LWASHTYP = 2
       end if



       if( TCNION(N) ) then
         LICETYP = 1
       else
         LICETYP = 2
       end if




       QTTOPAA = zero
       QTTOPCA = zero

       RCA  = zero
       FCA  = zero
       DCA  = zero
       RAMA = zero
       FAMA = zero
       DAMA = zero

       AMPCT      = zero
       AMCLPCT    = zero
       CLNEWPCT   = zero
       CLNEWAMPCT = zero
       CLOLDPCT   = zero
       CLOLDAMPCT = zero



       if( RLS(LE) > zero ) then
         CFXX(LE) = max( CFMIN,CFR(LE) )
       else
         CFXX(LE) = CFR(LE)
       endif

       rnew_flag(1:le) = .false.

level_loop : &
       do L = LE,1,-1
         LM1  = L - 1
         FAX  = zero
         RAX  = zero
         DAX  = zero
         FCXA = zero
         FCXB = zero
         DCXA = zero
         DCXB = zero
         RCXA = zero
         RCXB = zero

         QTDISCF   = zero
         QTDISRIME = zero
         QTDISCXA  = zero

         QTEVAPAXP = zero
         QTEVAPAXW = zero
         QTEVAPAX  = zero
         QTWASHAX  = zero

         QTEVAPCXAP = zero
         QTEVAPCXAW = zero
         QTEVAPCXA  = zero
         QTRIMECXA  = zero
         QTWASHCXA  = zero
         QTRAINCXA  = zero
         QTRAINCXB  = zero
         
         RAMPCT = zero
         RCXPCT = zero

         RCXALOC = zero
         RCXBLOC = zero
         RAXLOC  = zero
         RAMALOC = zero
         RCALOC  = zero

         RPRECIP       = zero
         DELTARIMEMASS = zero
         DELTARIME     = zero
         DOR           = zero
         DNEW          = zero

         QTTOPAAX = zero
         QTTOPCAX = zero

has_rls : &
         if( rls_flag(l) ) then








           FAX = max( zero,FAMA*(one - evaprate(l)) )
           RAX = RAMA								     
           if( FAMA > zero ) then
             if( freezing(l) ) then
               DAX = DAMA      
             else
               DAX = four    
             endif
           else
             DAX = zero
           endif

           if( RAMA > zero ) then
             QTEVAPAXP = min( QTTOPAA,EVAPRATE(L)*QTTOPAA )
           else
             QTEVAPAXP = zero
           endif
           if( is_hno3 ) then
             rax_wrk(l,1) = rax
             fax_wrk(l,1) = fax
           endif





           WRK = RAX*FAX + RCA*FCA
           if( WRK > 0. ) then
             RNEW_TST = RLS(L)/(GAREA * WRK)
           else
             RNEW_TST = 10.
           endif
           RNEW = RLSOG(L) - (RAX*FAX + RCA*FCA)     
           rnew_wrk(l) = rnew_tst



has_rnew:  if( rlsog(l) > adj_factor*(rax*fax + rca*fca) ) then





             if( cfxx(l) == zero ) then
               write(*,*) 'offline inputs'
               write(*,*) qttjfl(:,n)
               write(*,*) qm(:)
               write(*,*) pofl(:)
               write(*,*) delz(:)
               write(*,*) rls(:)
               write(*,*) clwc(:)
               write(*,*) ciwc(:)
               write(*,*) cfr(:)
               write(*,*) tem(:)
               write(*,*) evaprate(:)
               write(*,*) hstar(:,n)
               write(message,'('' washout: cloud fraction == 0 @ i,j,l,n = '',4i4)') ii,jj,l,n
               call wrf_debug( 15, trim(message) )
               QTTJFL(:lpar,N) = QTT(:lpar)
               cycle species_loop
             endif
             rnew_flag(l) = .true.
             CLWX = max( CLWC(L)+CIWC(L),CWMIN*CFXX(L) )
             if( is_hno3 ) then
               clwx_wrk(l) = clwx
             endif



             FCXA = FCA
             FCXB = max( zero,CFXX(L)-FCXA )







is_freezing : &
             if( freezing(l) ) then
               COLEFFSNOW = exp( 2.5e-2*(TEM(L) - TICE) )
               if( TEM(L) <= TFROZ ) then
                 RHOSNOW = RHOSNOWFIX
               else
                 RHOSNOW = 0.303*(TEM(L) - TFROZ)*RHOSNOWFIX
               endif
               if( FCXA > zero ) then
                 if( DCA > zero ) then
                   DELTARIMEMASS = CLWX*QM(L)*(FCXA/CFXX(L))* &
                     (one - exp( (-COLEFFSNOW/(DCA*1.e-3))*((RCA)/(2.*RHOSNOW))*DTSCAV ))   
                 else
                   DELTARIMEMASS = zero
                 endif
               else
                 DELTARIMEMASS = zero
               endif




               if( FCXA > zero ) then
                 DELTARIME = min( RNEW/FCXA,DELTARIMEMASS/(FCXA*GAREA*DTSCAV) ) 
               else
                 DELTARIME = zero
               endif
               if( is_hno3 ) then
                 deltarime_wrk(l) = deltarime
               endif



               if( RCA > zero ) then
                 DOR = max( DMIN,(((RCA+DELTARIME)/RCA)**VOLPOW)*DCA )
               else
                 DOR = zero
               endif






               RPRECIP = (RNEW-(DELTARIME*FCXA))/CFXX(L) 



               RCXA = RCA + DELTARIME + RPRECIP          
               RCXB = RPRECIP                            










               if( RPRECIP > zero ) then
                 WEMP = (CLWX*QM(L))/(GAREA*CFXX(L)*DELZ(L)) 
                 REMP = RPRECIP/((RHORAIN/1.e3))             
                 DNEW = DEMPIRICAL( WEMP, REMP )
                 DNEW = max( DMIN,DNEW )
                 if( FCXB > zero ) then
                   DCXB = DNEW
                 else
                   DCXB = zero
                 endif
               else
                 DCXB = zero
               endif

               if( FCXA > zero ) then
                 WEMP = (CLWX*QM(L)*(FCXA/CFXX(L)))/(GAREA*FCXA*DELZ(L)) 
                 REMP = RCXA/((RHORAIN/1.e3))                         
                 DEMP = DEMPIRICAL( WEMP, REMP )
                 DCXA = ((RCA+DELTARIME)/RCXA)*DOR + (RPRECIP/RCXA)*DNEW
                 DCXA = max( DEMP,DCXA )
                 DCXA = max( DMIN,DCXA )
               else
                 DCXA = zero
               endif

               if( QTT(L) > zero ) then   








                 if( RPRECIP > zero ) then
                   if( LICETYP == 1 ) then
                     RRAIN = RPRECIP*GAREA                                  
                     call DISGAS( CLWX, CFXX(L), TCMASS(N), HSTAR(L,N), &
                                  TEM(L),POFL(L),QM(L),                 &


                                  QTT(L)*CFXX(L),QTDISCF, is_hno3, RETEFF(N) )

                     call RAINGAS( RRAIN, DTSCAV, CLWX, CFXX(L),        &
                                   QM(L), QTT(L), QTDISCF, QTRAIN )
                     WRK       = QTRAIN/CFXX(L)
                     QTRAINCXA = FCXA*WRK
                     QTRAINCXB = FCXB*WRK
                   elseif( LICETYP == 2 ) then
                     QTRAINCXA = zero
                     QTRAINCXB = zero
                   endif
                 endif







                 if( DELTARIME > zero ) then
                   if( LICETYP == 1 ) then
                     if( TEM(L) <= TFROZ ) then
                       RHOSNOW = RHOSNOWFIX
                     else
                       RHOSNOW = 0.303*(TEM(L) - TFROZ)*RHOSNOWFIX
                     endif
                     QTCXA = QTT(L)*FCXA
                     call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),   &
                                  HSTAR(L,N), TEM(L), POFL(L),            &


                                  QM(L), QTCXA, QTDISRIME, is_hno3, RETEFF(N) )       

                     QTDISSTAR = (QTDISRIME*QTCXA)/(QTDISRIME + QTCXA)
                     QTRIMECXA = QTCXA*                             &
                        (one - exp((-COLEFFSNOW/(DCA*1.e-3))*       &
                        (RCA/(2.*RHOSNOW))*                         &  
                        (QTDISSTAR/QTCXA)*DTSCAV))
                     QTRIMECXA = min( QTRIMECXA, &               
                        ((RNEW*GAREA*DTSCAV)/(CLWX*QM(L)*(FCXA/CFXX(L))))*QTDISSTAR)
                   elseif( LICETYP == 2 ) then
                     QTRIMECXA = zero
                   endif
                 endif
               else
                 QTRAINCXA = zero
                 QTRAINCXB = zero
                 QTRIMECXA = zero
               endif



               QTWASHCXA = zero
               QTEVAPCXA = zero






             else is_freezing
               if( FCXA > zero ) then
                 DELTARIMEMASS = (CLWX*QM(L))*(FCXA/CFXX(L))*           &
                   (one - exp( -0.24*COLEFFRAIN*((RCA)**0.75)*DTSCAV ))  
               else
                 DELTARIMEMASS = zero
               endif




               if( FCXA > zero ) then
                 DELTARIME = min( RNEW/FCXA,DELTARIMEMASS/(FCXA*GAREA*DTSCAV) ) 
               else
                 DELTARIME = zero
               endif



               RPRECIP = (RNEW-(DELTARIME*FCXA))/CFXX(L)       

               RCXA = RCA + DELTARIME + RPRECIP            
               RCXB = RPRECIP                              
               DCXA = FOUR  
               if( FCXB > zero ) then
                 DCXB = FOUR
               else
                 DCXB = zero
               endif




               if( QTT(L) > zero ) then
                 if( RPRECIP > zero ) then
                   RRAIN = (RPRECIP*GAREA) 
                   call DISGAS( CLWX, CFXX(L), TCMASS(N), HSTAR(L,N), &
                                TEM(L), POFL(L), QM(L),               &


                                QTT(L)*CFXX(L), QTDISCF, is_hno3, RETEFF(N) )

                   call RAINGAS( RRAIN, DTSCAV, CLWX, CFXX(L),        &
                                 QM(L), QTT(L), QTDISCF, QTRAIN )
                   WRK       = QTRAIN/CFXX(L)
                   QTRAINCXA = FCXA*WRK
                   QTRAINCXB = FCXB*WRK
                 endif






                 if( DELTARIME > zero ) then
                   QTCXA = QTT(L)*FCXA
                   call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),    &
                                HSTAR(L,N), TEM(L), POFL(L),             &


                                QM(L), QTCXA, QTDISRIME, is_hno3, RETEFF(N) )

                   QTDISSTAR = (QTDISRIME*QTCXA)/(QTDISRIME + QTCXA)
                   QTRIMECXA = QTCXA*                              &
                      (one - exp(-0.24*COLEFFRAIN*                 &
                      ((RCA)**0.75)*                               & 
                      (QTDISSTAR/QTCXA)*DTSCAV))               
                   QTRIMECXA = min( QTRIMECXA, &
                      ((RNEW*GAREA*DTSCAV)/(CLWX*QM(L)*(FCXA/CFXX(L))))*QTDISSTAR)
                 else
                   QTRIMECXA = zero
                 endif
               else
                 QTRAINCXA = zero
                 QTRAINCXB = zero
                 QTRIMECXA = zero
               endif





               if( RCA > zero ) then
                 QTPRECIP = FCXA*QTT(L) - QTDISRIME
                 if( LWASHTYP == 1 ) then
                   if( QTPRECIP > zero ) then
                     QTWASHCXA = QTPRECIP*(one - exp( -0.24*COLEFFAER*((RCA)**0.75)*DTSCAV ))   
                   else
                     QTWASHCXA = zero
                   endif
                   QTEVAPCXA = zero
                 elseif( LWASHTYP == 2 ) then
                   RWASH = RCA*GAREA                                
                   if( QTPRECIP > zero ) then
                     call WASHGAS( RWASH, FCA, DTSCAV, QTTOPCA+QTRIMECXA, &
                                   HSTAR(L,N), TEM(L), POFL(L),           &
                                   QM(L), QTPRECIP, QTWASHCXA, QTEVAPCXA )
                   else
                     QTWASHCXA = zero
                     QTEVAPCXA = zero
                   endif
                 endif
               endif
             endif is_freezing





           else has_rnew
             CLWX = CLWC(L) + CIWC(L)
             if( is_hno3 ) then
               clwx_wrk(l) = clwx
             endif
             FCXA = FCA
             FCXB = max( zero,CFXX(L)-FCXA )
             RCXB = zero
             DCXB = zero
             QTRAINCXA = zero
             QTRAINCXB = zero
             QTRIMECXA = zero







             if( FCXA > zero ) then
               RCXA = min( RCA,RLS(L)/(GAREA*FCXA) )     
               if( FAX > zero .and. ((RCXA+1.e-12) < RLS(L)/(GAREA*FCXA)) ) then
                 RAXADJF = RLS(L)/GAREA - RCXA*FCXA
                 RAMPCT = RAXADJF/(RAX*FAX)
                 FAXADJ = RAMPCT*FAX
                 if( FAXADJ > zero ) then
                   RAXADJ = RAXADJF/FAXADJ
                 else
                   RAXADJ = zero
                 endif
               else
                 RAXADJ = zero
                 RAMPCT = zero
                 FAXADJ = zero
               endif
             else
               RCXA = zero
               if( FAX > zero ) then
                 RAXADJF = RLS(L)/GAREA
                 RAMPCT = RAXADJF/(RAX*FAX)
                 FAXADJ = RAMPCT*FAX
                 if( FAXADJ > zero ) then
                   RAXADJ = RAXADJF/FAXADJ
                 else
                   RAXADJ = zero
                 endif              
               else
                 RAXADJ = zero
                 RAMPCT = zero
                 FAXADJ = zero
               endif
             endif
  
             QTEVAPAXP = min( QTTOPAA,QTTOPAA - (RAMPCT*(QTTOPAA-QTEVAPAXP)) )
             FAX = FAXADJ
             RAX = RAXADJ






             if( RCXA <= zero ) then
               QTEVAPCXA = QTTOPCA
               RCXA = zero
               DCXA = zero
             else















is_freezing_a : &
               if( freezing(l) ) then
                 QTWASHCXA = zero
                 DCXA = ((RCXA/RCA)**VOLPOW)*DCA
                 if( LICETYP == 1 ) then
                   if( TEM(L) <= TMIX ) then
                     MASSLOSS = (RCA-RCXA)*FCXA*GAREA*DTSCAV



                     call DISGAS( (MASSLOSS/QM(L)), FCXA, TCMASS(N),   &
                                   HSTAR(L,N), TEM(L), POFL(L),        &


                                   QM(L), QTT(L), QTEVAPCXA, is_hno3, RETEFF(N) )

                     QTEVAPCXA = min( QTTOPCA,QTEVAPCXA )
                   else
                     QTEVAPCXA = zero
                   endif
                 elseif( LICETYP == 2 ) then   
                   QTEVAPCXA = zero
                 endif
               else is_freezing_a
                 QTEVAPCXAP = (RCA - RCXA)/RCA*QTTOPCA
                 DCXA = FOUR
                 QTCXA = FCXA*QTT(L)
                 if( LWASHTYP == 1 ) then
                   if( QTT(L) > zero ) then
                     call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),   &
                                  HSTAR(L,N), TEM(L), POFL(L),            &


                                  QM(L), QTCXA, QTDISCXA, is_hno3, RETEFF(N) )

                     if( QTCXA > QTDISCXA ) then
                       QTWASHCXA = (QTCXA - QTDISCXA)*(one - exp( -0.24*COLEFFAER*((RCXA)**0.75)*DTSCAV )) 
                     else
                       QTWASHCXA = zero
                     endif
                     QTEVAPCXAW = zero
                   else
                     QTWASHCXA  = zero
                     QTEVAPCXAW = zero
                   endif
                 elseif (LWASHTYP == 2 ) then
                   RWASH = RCXA*GAREA                         
                   call WASHGAS( RWASH, FCXA, DTSCAV, QTTOPCA, HSTAR(L,N), &
                                 TEM(L), POFL(L), QM(L),                   &
                                 QTCXA-QTDISCXA, QTWASHCXA, QTEVAPCXAW )
                 endif
                 QTEVAPCXA = QTEVAPCXAP + QTEVAPCXAW
               endif is_freezing_a
             endif
           endif has_rnew






           if( RAX > zero ) then

             if( .not. freezing(l) ) then
               QTAX = FAX*QTT(L)
               if( LWASHTYP == 1 ) then
                 QTWASHAX = QTAX*                        &
                    (one - exp(-0.24*COLEFFAER*       &
                   ((RAX)**0.75)*DTSCAV))  
                 QTEVAPAXW = zero
               elseif( LWASHTYP == 2 ) then
                 RWASH = RAX*GAREA   
                 call WASHGAS( RWASH, FAX, DTSCAV, QTTOPAA, HSTAR(L,N), &
                               TEM(L), POFL(L), QM(L), QTAX,            &
                               QTWASHAX, QTEVAPAXW )
               endif
             else
               QTEVAPAXW = zero
               QTWASHAX  = zero
             endif
           else
             QTEVAPAXW = zero
             QTWASHAX  = zero
           endif
           QTEVAPAX = QTEVAPAXP + QTEVAPAXW






           if( is_hno3 ) then
             rls_wrk(l) = rls(l)/garea
             rca_wrk(l) = rca
             fca_wrk(l) = fca
             rcxa_wrk(l) = rcxa
             fcxa_wrk(l) = fcxa
             rcxb_wrk(l) = rcxb
             fcxb_wrk(l) = fcxb
             rax_wrk(l,2) = rax
             fax_wrk(l,2) = fax
           endif
upper_level : &
           if( L > 1 ) then
             FAMA = max( FCXA + FCXB + FAX - CFR(LM1),zero )
             if( FAX > zero ) then
               RAXLOC = RAX/FAX
             else
               RAXLOC = zero
             endif
             if( FCXA > zero ) then
               RCXALOC = RCXA/FCXA
             else
               RCXALOC = zero
             endif
             if( FCXB > zero ) then
               RCXBLOC = RCXB/FCXB
             else
               RCXBLOC = zero
             endif

             if( CFR(LM1) >= CFMIN ) then
               CFXX(LM1) = CFR(LM1)
             else
               if( adj_factor*RLSOG(LM1) >= (RCXA*FCXA + RCXB*FCXB + RAX*FAX)*(one - EVAPRATE(LM1)) ) then
                 CFXX(LM1) = CFMIN
                 cf_trigger(lm1) = .true.
               else
                 CFXX(LM1) = CFR(LM1)
               endif
             endif




             if( FAX > zero ) then
               RAXLOC = RAX/FAX
               AMPCT = max( zero,min( one,(CFXX(L) + FAX - CFXX(LM1))/FAX ) )
               AMCLPCT = one - AMPCT
             else
               RAXLOC  = zero
               AMPCT   = zero
               AMCLPCT = zero
             endif
             if( FCXB > zero ) then
               RCXBLOC = RCXB/FCXB
               CLNEWPCT = max( zero,min( (CFXX(LM1) - FCXA)/FCXB,one ) )
               CLNEWAMPCT = one - CLNEWPCT
             else
               RCXBLOC    = zero
               CLNEWPCT   = zero
               CLNEWAMPCT = zero
             endif
             if( FCXA > zero ) then
               RCXALOC = RCXA/FCXA
               CLOLDPCT = max( zero,min( CFXX(LM1)/FCXA,one ) )
               CLOLDAMPCT = one - CLOLDPCT
             else
               RCXALOC    = zero
               CLOLDPCT   = zero
               CLOLDAMPCT = zero
             endif



             FCA = min( CFXX(LM1),FCXA*CLOLDPCT + CLNEWPCT*FCXB + AMCLPCT*FAX )
             if( FCA > zero ) then



               RCA = (RCXA*FCXA*CLOLDPCT + RCXB*FCXB*CLNEWPCT + RAX*FAX*AMCLPCT)/FCA
               if( RCA > zero ) then
                 DCA = (RCXA*FCXA*CLOLDPCT)/(RCA*FCA)*DCXA + & 
                       (RCXB*FCXB*CLNEWPCT)/(RCA*FCA)*DCXB + &
                       (RAX*FAX*AMCLPCT)/(RCA*FCA)*DAX
               else
                 DCA = zero
               endif
             else
               FCA = zero
               DCA = zero
               RCA = zero
             endif

             FAMA = FCXA + FCXB + FAX - CFXX(LM1)
             if( FAMA > zero ) then
               RAMA = (RCXA*FCXA*CLOLDAMPCT + RCXB*FCXB*CLNEWAMPCT + RAX*FAX*AMPCT)/FAMA
	       if( RAMA > zero ) then
                 DAMA = (RCXA*FCXA*CLOLDAMPCT)/(RAMA*FAMA)*DCXA +  &
                        (RCXB*FCXB*CLNEWAMPCT)/(RAMA*FAMA)*DCXB +  &
                        (RAX*FAX*AMPCT)/(RAMA*FAMA)*DAX
	       else
		  FAMA = zero
                  DAMA = zero
	       endif
             else
               FAMA = zero
               DAMA = zero
               RAMA = zero
             endif
           else upper_level
             AMPCT      = zero
             AMCLPCT    = zero
             CLNEWPCT   = zero
             CLNEWAMPCT = zero
             CLOLDPCT   = zero
             CLOLDAMPCT = zero
           endif upper_level
         else has_rls
	   RNEW = zero
           QTEVAPCXA = QTTOPCA
           QTEVAPAX = QTTOPAA
           if( L > 1 ) then
             if( RLS(LM1) > zero ) then
               CFXX(LM1) = max( CFMIN,CFR(LM1) )





             else
               CFXX(LM1) = CFR(LM1)
             endif
           endif
           AMPCT      = zero
           AMCLPCT    = zero
           CLNEWPCT   = zero
           CLNEWAMPCT = zero
           CLOLDPCT   = zero
           CLOLDAMPCT = zero
           RCA        = zero
           RAMA       = zero
           FCA        = zero
           FAMA       = zero
           DCA        = zero
           DAMA       = zero
         endif has_rls

         if( is_hno3 ) then
           fama_wrk(l) = fama
           rama_wrk(l) = rama
         endif



         QTNETLCXA = QTRAINCXA + QTRIMECXA + QTWASHCXA - QTEVAPCXA
         QTNETLCXA = min( QTT(L)*FCXA,QTNETLCXA )
   
         QTNETLCXB =QTRAINCXB
         QTNETLCXB = min( QTT(L)*FCXB,QTNETLCXB )

         QTNETLAX = QTWASHAX - QTEVAPAX
         QTNETLAX = min( QTT(L)*FAX,QTNETLAX )
              
         QTTNEW(L) = QTT(L) - (QTNETLCXA + QTNETLCXB + QTNETLAX)


         QTTOPCAX = (QTTOPCA + QTNETLCXA)*CLOLDPCT + QTNETLCXB*CLNEWPCT + (QTTOPAA + QTNETLAX)*AMCLPCT
         QTTOPAAX = (QTTOPCA + QTNETLCXA)*CLOLDAMPCT + QTNETLCXB*CLNEWAMPCT + (QTTOPAA + QTNETLAX)*AMPCT
         QTTOPCA = QTTOPCAX
         QTTOPAA = QTTOPAAX
       end do level_loop




       QTTJFL(:le,N) = QTTNEW(:le)
     end do species_loop

end subroutine WASHOUT

subroutine DISGAS( CLWX, CFX, MOLMASS, HSTAR, &
                   TM, PR, QM, &


                   QT, QTDIS, is_hno3, RETEFF )





      real, intent(in) :: CLWX,CFX    
      real, intent(in) :: MOLMASS     
      real, intent(in) :: HSTAR       
      real, intent(in) :: TM          
      real, intent(in) :: PR          
      real, intent(in) :: QM          
      real, intent(in) :: QT          
      real, intent(out) :: QTDIS      


      logical, intent(in) :: is_hno3
      real, intent(in) :: RETEFF      

 









      real  :: MUEMP











      if( TM >= TICE ) then
        QTDIS = (HSTAR*(QT/(QM*CFX))*0.029*(PR/1.0e3))*(CLWX*QM)


      elseif( TM <= TMIX .and. is_hno3 ) then

        MUEMP = exp( -14.2252 + TM*(1.55704e-1 - 7.1929e-4*TM) )
        QTDIS = MUEMP*(MOLMASS/18.)*(CLWX*QM)
      else
        QTDIS = RETEFF*((HSTAR*(QT/(QM*CFX))*0.029*(PR/1.0e3))*(CLWX*QM))
      endif

end subroutine DISGAS

subroutine RAINGAS( RRAIN, DTSCAV, CLWX, CFX, QM, QT, QTDIS, QTRAIN )


















      real, intent(in) :: RRAIN       
      real, intent(in) :: DTSCAV      
      real, intent(in) :: CLWX,CFX    
      real, intent(in) :: QM          
      real, intent(in) :: QT          
      real, intent(in) :: QTDIS       
      real, intent(out) :: QTRAIN     




      real ::  QTLF, QTDISSTAR

      QTDISSTAR = (QTDIS*(QT*CFX))/(QTDIS+(QT*CFX))
 



      QTLF = (RRAIN*QTDISSTAR)/(CLWX*QM*QT*CFX)





      QTRAIN = QT*CFX*(one - exp( -DTSCAV*QTLF ))

end subroutine RAINGAS

subroutine WASHGAS( RWASH, BOXF, DTSCAV, QTRTOP, HSTAR, TM, PR, QM, &
                    QT, QTWASH, QTEVAP )















      real, intent(in)  :: RWASH   
      real, intent(in)  :: BOXF    
      real, intent(in)  :: DTSCAV  
      real, intent(in)  :: QTRTOP  
      real, intent(in)  :: HSTAR   
      real, intent(in)  :: TM      
      real, intent(in)  :: PR      
      real, intent(in)  :: QT      
      real, intent(in)  :: QM      
      real, intent(out) :: QTWASH  
      real, intent(out) :: QTEVAP  
      



      real            :: FWASH, QTMAX, QTDIF










      FWASH = (RWASH*HSTAR*29.e-6*PR)/(QM*BOXF)



      QTMAX = QT*FWASH*DTSCAV
      if( QTMAX > QTRTOP ) then



         QTDIF = min( QT,QTMAX-QTRTOP )
         QTWASH = QTDIF * (one - exp( -DTSCAV*FWASH ))
         QTEVAP = zero
      else



         QTWASH = zero
         QTEVAP = QTRTOP - QTMAX
      endif
     
end subroutine WASHGAS

function DEMPIRICAL( CWATER, RRATE )




      real, intent(in)  :: CWATER   
      real, intent(in)  :: RRATE




      real(8), parameter :: big_diameter = 100._8            
      real(8), parameter :: const0       = .638_8
      real(8), parameter :: const1       = oner8 + const0

      real(8) :: RRATEX, WX, THETA, PHI, ETA, BETA, ALPHA, BEE
      real(8) :: GAMTHETA, GAMBETA
      real(8) :: numer, denom
      real(8) :: diameter                      

      real :: DEMPIRICAL

      if( cwater > 0._8 ) then
        RRATEX = real(RRATE,kind=8)*3600._8       
        WX     = real(CWATER,kind=8)*1.0e3_8      

        if( RRATEX > 0.04_8 ) then
           THETA = exp( -1.43_8*log10( 7._8*RRATEX ) ) + 2.8_8
        else
           THETA = 5._8
        endif

        PHI = RRATEX/(3600._8*10._8) 
        ETA = exp( 3.01_8*THETA - 10.5_8 )
        BETA  = THETA/const1
        ALPHA = exp( FOURR8*(BETA - 3.5_8) )
        BEE   = const0*BETA - ONER8
        GAMTHETA = real( GAMMA( real(THETA,kind=4) ),kind=8 )
        GAMBETA  = real( GAMMA( real(BETA + ONER8,kind=4) ),kind=8 )

        numer = WX*ETA*GAMTHETA
        denom = 1.0e6_8*ALPHA*PHI*GAMBETA
        diameter = ((numer/denom)**(-oner8/BEE))*10._8
        diameter = min( big_diameter,diameter )
        DEMPIRICAL = real( diameter )

      else
        DEMPIRICAL = 0.
      endif

end function DEMPIRICAL

function GAMMA( X )










        real, intent(in)  :: X




        real, parameter :: PI = 3.141592653589793e0

        integer :: k, M, M1
        real    :: GR, R, Z 
        real    :: G(26)




        real :: GAMMA

        DATA G/1.0e0,0.5772156649015329,                     &
            -0.6558780715202538e0, -0.420026350340952e-1,    &
            0.1665386113822915e0,-.421977345555443e-1,       &
            -.96219715278770e-2, .72189432466630e-2,         &
            -.11651675918591e-2, -.2152416741149e-3,         &
            .1280502823882e-3, -.201348547807e-4,            &
            -.12504934821e-5, .11330272320e-5,               &
            -.2056338417e-6, .61160950e-8,                   &
            .50020075e-8, -.11812746e-8,                     &
            .1043427e-9, .77823e-11,                         &
            -.36968e-11, .51e-12,                            &
            -.206e-13, -.54e-14, .14e-14, .1e-15/

is_integer : &
        IF( x == real( int(x) ) ) then
          IF( X > zero ) THEN
            GAMMA = ONE
            M1 = INT(X) - 1
            DO K = 2,M1
              GAMMA = GAMMA*real(K)
            END DO
          ELSE
            GAMMA = 1.0e36
          ENDIF
        ELSE is_integer
          IF( ABS(X) > ONE ) THEN
            Z = ABS(X)
            M = INT(Z)
            R = ONE
            DO K = 1,M
              R = R*(Z - real(k))
            END DO
            Z = Z - real(M)
          ELSE
            Z = X
          ENDIF
          GR = G(26)
          DO K = 25,1,-1
            GR = GR*Z + G(K)
          end DO
          GAMMA = ONE/(GR*Z)
          IF( ABS(X) > ONE ) THEN
            GAMMA = GAMMA*R
            IF( X < zero ) then
              GAMMA = -PI/(X*GAMMA*SIN( PI*X ))
            ENDIF
          ENDIF
        ENDIF is_integer

END function GAMMA

END MODULE module_mozcart_wetscav
