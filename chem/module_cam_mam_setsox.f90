

      module MO_SETSOX

      use shr_kind_mod, only : r8 => shr_kind_r8
      use module_cam_support, only: iulog
      use modal_aero_data

      private
      public :: sox_inti, setsox
      public :: has_sox

      save

      integer, target    ::  spc_ids(8)
      integer, pointer   ::  id_so2, id_so4_1, id_so4_2, id_so4_3, &
                             id_h2o2, id_ox, id_ho2, id_h2so4
      integer            ::  id_nh3
      logical            ::  inv_o3
      integer            ::  id_hno3, id_msa
      logical            ::  has_sox = .true.
      logical :: inv_so2, inv_nh3, inv_hno3, inv_h2o2, inv_ox, inv_nh4no3, inv_ho2
      contains

      subroutine sox_inti




      use mo_chem_utls, only : get_spc_ndx, get_inv_ndx
      use module_cam_support, only: masterproc
      
      implicit none

      id_so2    => spc_ids(1)
      id_so4_1  => spc_ids(2)
      id_so4_2  => spc_ids(3)
      id_so4_3  => spc_ids(4)

      id_h2o2   => spc_ids(5)
      id_ox     => spc_ids(6)
      id_ho2    => spc_ids(7)
      id_h2so4  => spc_ids(8)



      id_so4_1  = lptr_so4_cw_amode(1)
      id_so4_2  = lptr_so4_cw_amode(2)
      id_so4_3  = lptr_so4_cw_amode(3)

      id_h2so4  = get_spc_ndx( 'H2SO4' )
      id_msa    = get_spc_ndx( 'MSA' )

      inv_so2 = .false.
      id_so2    = get_inv_ndx( 'SO2' )
      inv_so2 = id_so2 > 0
      if ( .not. inv_so2 ) then
         id_so2    = get_spc_ndx( 'SO2' )
      endif

      inv_NH3 = .false.
      id_NH3    = get_inv_ndx( 'NH3' )
      inv_NH3 = id_NH3 > 0
      if ( .not. inv_NH3 ) then
         id_NH3    = get_spc_ndx( 'NH3' )
      endif


      inv_HNO3 = .false.
      id_HNO3    = get_inv_ndx( 'HNO3' )
      inv_HNO3 = id_hno3 > 0
      if ( .not. inv_HNO3 ) then
         id_HNO3    = get_spc_ndx( 'HNO3' )
      endif

      inv_H2O2 = .false.
      id_H2O2    = get_inv_ndx( 'H2O2' )
      inv_H2O2 = id_H2O2 > 0
      if ( .not. inv_H2O2 ) then
         id_H2O2    = get_spc_ndx( 'H2O2' )
      endif

      inv_HO2 = .false.
      id_HO2    = get_inv_ndx( 'HO2' )
      inv_HO2 = id_HO2 > 0
      if ( .not. inv_HO2 ) then
         id_HO2    = get_spc_ndx( 'HO2' )
      endif

      inv_o3    = get_inv_ndx( 'O3' ) > 0
      if (inv_o3) then
         id_ox  = get_inv_ndx( 'O3' )
      else
         id_ox  = get_spc_ndx( 'O3' )
      endif
      inv_ho2   = get_inv_ndx( 'HO2' ) > 0
      if (inv_ho2) then
         id_ho2 = get_inv_ndx( 'HO2' )
      else
         id_ho2 = get_spc_ndx( 'HO2' )
      endif

      has_sox = all( spc_ids(1:8) > 0 )

      if (masterproc) then
         write(iulog,*) 'sox_inti: has_sox = ',has_sox
         call wrf_message(iulog)
      endif

      if( has_sox ) then
         if (masterproc) then
            write(iulog,*) '-----------------------------------------'
            call wrf_message(iulog)
            write(iulog,*) 'mozart will do sox aerosols'
            call wrf_message(iulog)
            write(iulog,*) '-----------------------------------------'
            call wrf_message(iulog)
         endif
      end if

      end subroutine sox_inti

      subroutine SETSOX( ncol,   &
           press,  &
           dtime,  &
           tfld,   &
           qfld,   &
           lwc,    &
           lchnk,  &
           pdel,   &
           mbar,   &
           clwlrat,&
           cldfrc, &
           cldnum, &
           qcw,    &
           loffset,&
           xhnm,   &
           qin, invariants )

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        use module_cam_support, only: pcols, pver, gas_pcnst => gas_pcnst_modal_aero, &
             nfs
        use module_data_cam_mam_asect, only:  adv_mass => mw_q_mo_array
        use physconst,    only : mwdry, gravit
        use physconst,    only : pi
        use module_cam_support, only: outfld
        
        implicit none
        
        
        
        
        integer, intent(in)                                  :: ncol
        real(r8), intent(in)                                     :: dtime                     
        real(r8), dimension(ncol ,pver,gas_pcnst), intent(inout) :: qin                       
        real(r8), dimension(ncol ,pver), intent(in)              :: xhnm                      
        real(r8), dimension(pcols,pver), intent(in)              :: tfld, &                   
             qfld, &                   
             press                     
        real(r8), dimension(ncol ,pver), intent(in)              :: lwc  
        real(r8), intent(in)    ::  invariants(ncol,pver,nfs)
        integer,  intent(in)                                     :: lchnk                     
        real(r8), dimension(pcols,pver), intent(in)              :: pdel                      
        real(r8), dimension(ncol ,pver), intent(in)              :: mbar                      
        real(r8), dimension(ncol ,pver), intent(in)              :: cldfrc, &                 
                                                                    clwlrat,&                 
                                                                    cldnum                    
        real(r8), intent(inout) ::  qcw(ncol,pver,gas_pcnst)       
        integer,  intent(in)    ::  loffset                        
        
        
        
        
        
        integer, parameter :: itermax = 20
        real(r8), parameter ::  ph0 = 5.0_r8  
        real(r8), parameter ::  const0 = 1.e3_r8/6.023e23_r8
        real(r8), parameter ::  xa0 = 11._r8,   &
             xb0 = -.1_r8,   &
             xa1 = 1.053_r8, &
             xb1 = -4.368_r8,&
             xa2 = 1.016_r8, &
             xb2 = -2.54_r8, &
             xa3 = .816e-32_r8, &
             xb3 = .259_r8

        real(r8), parameter ::  kh0 = 9.e3_r8, &           
             kh1 = 2.05e-5_r8, &        
             kh2 = 8.6e5_r8,   &        
             kh3 = 1.e8_r8,    &        
             Ra = 8314._r8/101325._r8, &   
             xkw = 1.e-14_r8            
        real(r8), parameter :: small_value = 1.e-20_r8
        
        integer  :: l, n, m
        integer  :: ntot_msa_c

        integer  :: id_so4_1a, id_so4_2a, id_so4_3a, id_so4_4a, id_so4_5a, id_so4_6a, &
                    id_nh4_1a, id_nh4_2a, id_nh4_3a, id_nh4_4a, id_nh4_5a, id_nh4_6a

        real(r8) :: delso4_hprxn, delso4_o3rxn,                  &
                    dso4dt_aqrxn, dso4dt_hprxn,                  &
                    dso4dt_gasuptk, dmsadt_gasuptk,              &
                    dmsadt_gasuptk_tomsa, dmsadt_gasuptk_toso4,  &
                    dqdt_aq, dqdt_wr, dqdt,                      &
                    rad_cd, radxnum_cd, num_cd,                  &
                    gasdiffus, gasspeed, knudsen, uptkrate,      &
                    fuchs_sutugin, volx34pi_cd,                  &
                    fwetrem, sumf,                               &
                    frso2_g, frso2_c, frh2o2_g, frh2o2_c
        real(r8) :: delnh3, delnh4

        real(r8) :: faqgain_msa(ntot_amode), faqgain_so4(ntot_amode), &
                    qnum_c(ntot_amode)
        real(r8) :: dqdt_aqso4(ncol,pver,gas_pcnst), &
                    dqdt_aqh2so4(ncol,pver,gas_pcnst), &
                    dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver), &
                    xphlwc(ncol,pver), &
                    sflx(1:ncol)
        integer  :: k, i, iter, file
        real(r8) :: wrk, delta
        real(r8) :: xph0, aden, xk, xe, x2
        real(r8) :: tz, xl, px, qz, pz, es, qs, patm
        real(r8) :: Eso2, Eso4, Ehno3, Eco2, Eh2o, Enh3
        real(r8) :: hno3g, nh3g, so2g, h2o2g, co2g, o3g
        real(r8) :: hno3a, nh3a, so2a, h2o2a, co2a, o3a
        real(r8) :: rah2o2, rao3, pso4, ccc
        real(r8) :: cnh3, chno3, com, com1, com2, xra
        
        
        
        
        
        real(r8) :: kh4    
        real(r8) :: xam    
        real(r8) :: ho2s   
        real(r8) :: r1h2o2 
        real(r8) :: r2h2o2 

        real(r8), dimension(ncol,pver)  ::             &
             xhno3, xh2o2, xso2, xso4,&
             xnh3, xnh4, xo3,         &
             xlwc, cfact, &
             xph, xho2,         &
             xh2so4, xmsa, xso4_init, xso4c, xnh4c, &
             hehno3, &            
             heh2o2, &            
             heso2,  &            
             henh3,  &            
             heo3              
        real(r8), dimension(ncol)  :: work1
        logical :: converged

        id_so4_1a = id_so4_1 - loffset
        id_so4_2a = id_so4_2 - loffset
        id_so4_3a = id_so4_3 - loffset


        dqdt_aqso4(:,:,:) = 0.0_r8
        dqdt_aqh2so4(:,:,:) = 0.0_r8
        dqdt_aqhprxn(:,:) = 0.0_r8
        dqdt_aqo3rxn(:,:) = 0.0_r8
        
        
        
        
        
        
        
        
        
        
        xph0 = 10._r8**(-ph0)                      

        do k = 1,pver
           cfact(:,k) = xhnm(:,k)     &          
                * 1.e6_r8             &          
                * 1.38e-23_r8/287._r8 &          
                * 1.e-3_r8                       
        end do

        do k = 1,pver
           xph(:,k) = xph0                                
           xlwc(:,k) = lwc(:,k) *cfact(:,k)               

           if ( inv_so2 ) then
              xso2 (:,k) = invariants(:,k,id_so2)                   
           else
              xso2 (:,k) = qin(:,k,id_so2)                   
           endif

           if (id_hno3 > 0) then
             xhno3(:,k) = qin(:,k,id_hno3)
           else
             xhno3(:,k) = 0.0_r8
           endif

           if ( inv_h2o2 ) then
              xh2o2 (:,k) = invariants(:,k,id_h2o2)                   
           else
              xh2o2 (:,k) = qin(:,k,id_h2o2)                   
           endif

           if (id_nh3  > 0) then
             xnh3 (:,k) = qin(:,k,id_nh3)
           else
             xnh3 (:,k) = 0.0_r8
           endif


           if ( inv_o3 ) then
             xo3  (:,k) = invariants(:,k,id_ox)/xhnm(:,k) 
           else
             xo3  (:,k) = qin(:,k,id_ox)                  
           endif
           if ( inv_ho2 ) then
             xho2 (:,k) = invariants(:,k,id_ho2)/xhnm(:,k)
           else
             xho2 (:,k) = qin(:,k,id_ho2)                 
           endif

           xso4c(:,k) = qcw(:,k,id_so4_1a)   &
                      + qcw(:,k,id_so4_2a)   &
                      + qcw(:,k,id_so4_3a)
           xh2so4(:,k)= qin(:,k,id_h2so4)
           if (id_msa > 0) xmsa (:,k) = qin(:,k,id_msa)

           xnh4c(:,k) = xso4c(:,k) * 1.0_r8     
        end do

        
        
        
        do k = 1,pver                                             
           do i = 1,ncol
           if (cldfrc(i,k) >=  1.0e-5_r8) then
              xso4(i,k) = xso4c(i,k) / cldfrc(i,k)
              xnh4(i,k) = xnh4c(i,k) / cldfrc(i,k)
              xl = xlwc(i,k) / cldfrc(i,k)     
              if( xl >= 1.e-8_r8 ) then
                 work1(i) = 1._r8 / tfld(i,k) - 1._r8 / 298._r8
                 
                 
                 
                 do iter = 1,itermax
                    xk = 2.1e5_r8 *EXP( 8700._r8*work1(i) )
                    xe = 15.4_r8
                    hehno3(i,k)  = xk*(1._r8 + xe/xph(i,k))
                    
                    
                    
                    xk = 7.4e4_r8   *EXP( 6621._r8*work1(i) )
                    xe = 2.2e-12_r8 *EXP(-3730._r8*work1(i) )
                    heh2o2(i,k)  = xk*(1._r8 + xe/xph(i,k))
                    
                    
                    
                    xk = 1.23_r8  *EXP( 3120._r8*work1(i) )
                    xe = 1.7e-2_r8*EXP( 2090._r8*work1(i) )
                    x2 = 6.0e-8_r8*EXP( 1120._r8*work1(i) )
                    wrk = xe/xph(i,k)
                    heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))
                    
                    
                    
                    xk = 58._r8   *EXP( 4085._r8*work1(i) )
                    xe = 1.7e-5_r8*EXP(-4325._r8*work1(i) )
                    henh3(i,k)  = xk*(1._r8 + xe*xph(i,k)/xkw)
                    
                    
                    
                    pz = .01_r8*press(i,k)       
                    tz = tfld(i,k)
                    patm = pz/1013._r8
                    xam  = press(i,k)/(1.38e-23_r8*tz)  
                    
                    
                    
                    px = hehno3(i,k) * Ra * tz * xl
                    hno3g = xhno3(i,k)/(1._r8 + px)
                    xk = 2.1e5_r8 *EXP( 8700._r8*work1(i) )
                    xe = 15.4_r8
                    Ehno3 = xk*xe*hno3g *patm
                    
                    
                    
                    px = heso2(i,k) * Ra * tz * xl
                    so2g =  xso2(i,k)/(1._r8+ px)
                    xk = 1.23_r8  *EXP( 3120._r8*work1(i) )
                    xe = 1.7e-2_r8*EXP( 2090._r8*work1(i) )
                    Eso2 = xk*xe*so2g *patm
                    
                    
                    
                    px = henh3(i,k) * Ra * tz * xl
                    nh3g = xnh4(i,k)/px
                    xk = 58._r8   *EXP( 4085._r8*work1(i) )
                    xe = 1.7e-5_r8*EXP( -4325._r8*work1(i) )
                    Enh3 = xk*xe*nh3g/xkw *patm
                    
                    
                    
                    Eh2o = xkw
                    
                    
                    
                    co2g = 330.e-6_r8                            
                    xk = 3.1e-2_r8*EXP( 2423._r8*work1(i) )
                    xe = 4.3e-7_r8*EXP(-913._r8 *work1(i) )
                    Eco2 = xk*xe*co2g  *patm
                    
                    
                    
                    com2 = (Eh2o + Ehno3 + Eso2 + Eco2)  &
                         / (1._r8 + Enh3 )
                    com2 = MAX( com2,1.e-20_r8 )
                    xph(i,k) = SQRT( com2 )
                    
                    
                    
                    Eso4 = xso4(i,k)*xhnm(i,k)   &         
                         *const0/xl
                    xph(i,k) =  MIN( 1.e-2_r8,MAX( 1.e-7_r8,xph(i,k) + 2._r8*Eso4 ) )
                    if( iter > 1 ) then
                       delta = ABS( (xph(i,k) - delta)/delta )
                       converged = delta < .01_r8
                       if( converged ) then
                          exit
                       else
                          delta = xph(i,k)
                       end if
                    else
                       delta = xph(i,k)
                    end if
                 end do
                 if( .not. converged ) then
                    write(iulog,*) 'SETSOX: pH failed to converge @ (',i,',',k,'), % change=', &
                         100._r8*delta
                    call wrf_message(iulog)
                 end if
              else
                 xph(i,k) =  1.e-7_r8
              end if
           end if   
           end do
        end do  

        
        
        
        do k = 1,pver
           do i = 1,ncol
              work1(i) = 1._r8 / tfld(i,k) - 1._r8 / 298._r8
              tz = tfld(i,k)
            if (cldfrc(i,k) >=  1.0e-5_r8) then
              xl = xlwc(i,k) / cldfrc(i,k)
              patm = press(i,k)/101300._r8        
              xam  = press(i,k)/(1.38e-23_r8*tz)  

              
              
              
              xk = 7.4e4_r8   *EXP( 6621._r8*work1(i) )
              xe = 2.2e-12_r8 *EXP(-3730._r8*work1(i) )
              heh2o2(i,k)  = xk*(1._r8 + xe/xph(i,k))

              
              
              
              xk = 1.23_r8  *EXP( 3120._r8*work1(i) )
              xe = 1.7e-2_r8*EXP( 2090._r8*work1(i) )
              x2 = 6.0e-8_r8*EXP( 1120._r8*work1(i) )

              wrk = xe/xph(i,k)
              heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))

              
              
              
              xk = 58._r8   *EXP( 4085._r8*work1(i) )
              xe = 1.7e-5_r8*EXP(-4325._r8*work1(i) )
              henh3(i,k)  = xk*(1._r8 + xe*xph(i,k)/xkw)

              
              
              
              xk = 1.15e-2_r8 *EXP( 2560._r8*work1(i) )
              heo3(i,k) = xk

              
              
              
              
              kh4 = (kh2 + kh3*kh1/xph(i,k)) / ((1._r8 + kh1/xph(i,k))**2)
              ho2s = kh0*xho2(i,k)*patm*(1._r8 + kh1/xph(i,k))  
              r1h2o2 = kh4*ho2s*ho2s                         
              r2h2o2 = r1h2o2*xl               &             
                             /const0*1.e+6_r8  &             
                             /xam
              r2h2o2 = 0.0_r8
              xh2o2(i,k) = xh2o2(i,k) + r2h2o2*dtime         

              
              
              
              
              
              
              px = heh2o2(i,k) * Ra * tz * xl
              h2o2g =  xh2o2(i,k)/(1._r8+ px)

              
              
              
              px = heso2(i,k) * Ra * tz * xl
              so2g =  xso2(i,k)/(1._r8+ px)

              
              
              
              px = heo3(i,k) * Ra * tz * xl
              o3g =  xo3(i,k)/(1._r8+ px)

              
              
              
              px = henh3(i,k) * Ra * tz * xl
              if(xl .ge. 1.e-8_r8) then
                nh3g = xnh4(i,k)/px
              endif
              
              
              
              
              

              
              
              
              rah2o2 = 8.e4_r8 * EXP( -3650._r8*work1(i) )  &
                   / (.1_r8 + xph(i,k))

              
              
              
              rao3   = 4.39e11_r8 * EXP(-4131._r8/tz)  &
                   + 2.56e3_r8  * EXP(-996._r8 /tz) /xph(i,k)

              
              
              
              
              
              
              
              
              
              
              
              
              
              
              

              IF (XL .ge. 1.e-8_r8) THEN    

                 pso4 = rah2o2 * 7.4e4_r8*EXP(6621._r8*work1(i)) * h2o2g * patm &
                               * 1.23_r8 *EXP(3120._r8*work1(i)) * so2g  * patm
                 pso4 = pso4       &                          
                      * xl         &                          
                      / const0     &                          
                      / xhnm(i,k)

                 ccc = pso4*dtime
                 ccc = max(ccc, 1.e-30_r8)

                 xso4_init(i,k)=xso4(i,k)
                 IF (xh2o2(i,k) .gt. xso2(i,k)) THEN
                    if (ccc .gt. xso2(i,k)) then
                       xso4(i,k)=xso4(i,k)+xso2(i,k)
                       xh2o2(i,k)=xh2o2(i,k)-xso2(i,k)
                       xso2(i,k)=1.e-20_r8
                    else
                       xso4(i,k)  = xso4(i,k)  + ccc
                       xh2o2(i,k) = xh2o2(i,k) - ccc
                       xso2(i,k)  = xso2(i,k)  - ccc
                    end if

                 ELSE
                    if (ccc  .gt. xh2o2(i,k)) then
                       xso4(i,k)=xso4(i,k)+xh2o2(i,k)
                       xso2(i,k)=xso2(i,k)-xh2o2(i,k)
                       xh2o2(i,k)=1.e-20_r8
                    else
                       xso4(i,k)  = xso4(i,k)  + ccc
                       xh2o2(i,k) = xh2o2(i,k) - ccc
                       xso2(i,k)  = xso2(i,k)  - ccc
                    end if
                 END IF

                 delso4_hprxn = xso4(i,k) - xso4_init(i,k)
                 
                 
                 

                 pso4 = rao3 * heo3(i,k)*o3g*patm * heso2(i,k)*so2g*patm   
                 pso4 = pso4        &                                
                      * xl          &                                
                      / const0      &                                
                      / xhnm(i,k)                                    
                 ccc = pso4*dtime
                 ccc = max(ccc, 1.e-30_r8)

                 xso4_init(i,k)=xso4(i,k)

                 if (ccc .gt. xso2(i,k)) then
                    xso4(i,k)=xso4(i,k)+xso2(i,k)
                    xso2(i,k)=1.e-20_r8
                 else
                    xso4(i,k)  = xso4(i,k)  + ccc
                    xso2(i,k)  = xso2(i,k)  - ccc
                 end if

                 delso4_o3rxn = xso4(i,k) - xso4_init(i,k)












        do n = 1, ntot_amode
            qnum_c(n) = 0.0
            l = numptrcw_amode(n) - loffset
            if (l > 0) qnum_c(n) = max( 0.0_r8, qcw(i,k,l) )
        end do


        n = modeptr_accum
        if (n <= 0) n = 1
        qnum_c(n) = max( 1.0e-10_r8, qnum_c(n) )



        sumf = 0.0
        do n = 1, ntot_amode
            faqgain_so4(n) = 0.0
            if (lptr_so4_cw_amode(n) > 0) then
                faqgain_so4(n) = qnum_c(n)
                sumf = sumf + faqgain_so4(n)
            end if
        end do

        if (sumf > 0.0) then
            do n = 1, ntot_amode
                faqgain_so4(n) = faqgain_so4(n) / sumf
            end do
        end if



        ntot_msa_c = 0
        sumf = 0.0
        do n = 1, ntot_amode
            faqgain_msa(n) = 0.0
            if (lptr_msa_cw_amode(n) > 0) then
                faqgain_msa(n) = qnum_c(n)
                ntot_msa_c = ntot_msa_c + 1
            end if
            sumf = sumf + faqgain_msa(n)
        end do

        if (sumf > 0.0) then
            do n = 1, ntot_amode
                faqgain_msa(n) = faqgain_msa(n) / sumf
            end do
        end if










        num_cd = 1.0e-3_r8*cldnum(i,k)*cfact(i,k)/cldfrc(i,k)
        num_cd = max( num_cd, 0.0_r8 )






        volx34pi_cd = xl*0.75_r8/pi


        radxnum_cd = (volx34pi_cd*num_cd*num_cd)**0.3333333_r8


        if (radxnum_cd .le. volx34pi_cd*4.0e4_r8) then
            radxnum_cd = volx34pi_cd*4.0e4_r8
            rad_cd = 50.0e-4_r8
        else if (radxnum_cd .ge. volx34pi_cd*4.0e8_r8) then
            radxnum_cd = volx34pi_cd*4.0e8_r8
            rad_cd = 0.5e-4_r8
        else
            rad_cd = radxnum_cd/num_cd
        end if



        gasdiffus = 0.557_r8 * (tfld(i,k)**1.75_r8) / press(i,k)


        gasspeed  = 1.455e4_r8 * sqrt(tfld(i,k)/98.0_r8)


        knudsen = 3.0_r8*gasdiffus/(gasspeed*rad_cd)





        fuchs_sutugin = (0.4875_r8*(1. + knudsen)) /   &
                        (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)


        uptkrate = 12.56637_r8*radxnum_cd*gasdiffus*fuchs_sutugin


        uptkrate = (1.0_r8 - exp(-min(100._r8,dtime*uptkrate))) / dtime



        dso4dt_gasuptk = xh2so4(i,k) * uptkrate
        if (id_msa > 0) then
           dmsadt_gasuptk = xmsa(i,k) * uptkrate
        else
           dmsadt_gasuptk = 0.0_r8
        end if


        dmsadt_gasuptk_toso4 = 0.0_r8
        dmsadt_gasuptk_tomsa = dmsadt_gasuptk
        if (ntot_msa_c == 0) then
            dmsadt_gasuptk_tomsa = 0.0_r8
            dmsadt_gasuptk_toso4 = dmsadt_gasuptk
        end if







        dso4dt_aqrxn = (delso4_o3rxn + delso4_hprxn) / dtime
        dso4dt_hprxn = delso4_hprxn / dtime



        fwetrem = 0.0    


        do n = 1, ntot_amode
            l = lptr_so4_cw_amode(n) - loffset
            if (l > 0) then
                dqdt_aqso4(i,k,l) = faqgain_so4(n)*dso4dt_aqrxn*cldfrc(i,k)
                dqdt_aqh2so4(i,k,l) = faqgain_so4(n)*  &
                         (dso4dt_gasuptk + dmsadt_gasuptk_toso4)*cldfrc(i,k)
                dqdt_aq = dqdt_aqso4(i,k,l) + dqdt_aqh2so4(i,k,l)
                dqdt_wr = -fwetrem*dqdt_aq
                dqdt= dqdt_aq + dqdt_wr
                qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime
            end if

            l = lptr_msa_cw_amode(n) - loffset
            if (l > 0) then
                dqdt_aq = faqgain_msa(n)*dmsadt_gasuptk_tomsa*cldfrc(i,k)
                dqdt_wr = -fwetrem*dqdt_aq
                dqdt = dqdt_aq + dqdt_wr
                qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime
            end if

            l = lptr_nh4_cw_amode(n) - loffset
            if (l > 0) then
                if (delnh4 > 0.0_r8) then
                   dqdt_aq = faqgain_so4(n)*delnh4/dtime*cldfrc(i,k)
                   dqdt = dqdt_aq
                   qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime
                else
                   dqdt = (qcw(i,k,l)/max(xnh4c(i,k),1.0e-35_r8)) &
                           *delnh4/dtime*cldfrc(i,k)
                   qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime
                endif
            end if
        end do









        qin(i,k,id_h2so4) = qin(i,k,id_h2so4) - dso4dt_gasuptk * dtime * cldfrc(i,k)
        if (id_msa > 0) qin(i,k,id_msa) = qin(i,k,id_msa) - dmsadt_gasuptk * dtime * cldfrc(i,k)


        frso2_g = 1.0_r8 / ( 1.0_r8 + heso2(i,k) * Ra * tz * xl )
        frh2o2_g = 1.0_r8 / ( 1.0_r8 + heh2o2(i,k) * Ra * tz * xl )


        frso2_c = max( 0.0_r8, (1.0_r8-frso2_g) )
        frh2o2_c = max( 0.0_r8, (1.0_r8-frh2o2_g) )



        fwetrem = 0.0   

        dqdt_wr = -fwetrem*xso2(i,k)/dtime*cldfrc(i,k)
        dqdt_aq = -dso4dt_aqrxn*cldfrc(i,k)
        dqdt = dqdt_aq + dqdt_wr
        qin(i,k,id_so2) = qin(i,k,id_so2) + dqdt * dtime



        fwetrem = 0.0   

        dqdt_wr = -fwetrem*xh2o2(i,k)/dtime*cldfrc(i,k)
        dqdt_aq = -dso4dt_hprxn*cldfrc(i,k)
        dqdt = dqdt_aq + dqdt_wr
        qin(i,k,id_h2o2) = qin(i,k,id_h2o2) + dqdt * dtime




        dqdt_aqhprxn(i,k) = dso4dt_hprxn*cldfrc(i,k)
        dqdt_aqo3rxn(i,k) = (dso4dt_aqrxn - dso4dt_hprxn)*cldfrc(i,k)

              END IF 

              end if 

           end do
        end do

        
        
        
        do k = 1,pver
           qin(:,k,id_so2)   =  MAX( qin(:,k,id_so2), small_value )
           qcw(:,k,id_so4_1a) =  MAX( qcw(:,k,id_so4_1a), small_value )
           qcw(:,k,id_so4_2a) =  MAX( qcw(:,k,id_so4_2a), small_value )
           qcw(:,k,id_so4_3a) =  MAX( qcw(:,k,id_so4_3a), small_value )

        end do

        do n = 1, ntot_amode
            m = lptr_so4_cw_amode(n)
            l = m - loffset
            if (l > 0) then
              sflx(:)=0._r8
              do k=1,pver
                 do i=1,ncol
                    sflx(i)=sflx(i)+dqdt_aqso4(i,k,l)*adv_mass(l)/mbar(i,k)  &
                                   *pdel(i,k)/gravit  
                 enddo
              enddo
              call outfld( trim(cnst_name_cw(m))//'AQSO4', sflx(:ncol), ncol, lchnk)

              sflx(:)=0._r8
              do k=1,pver
                 do i=1,ncol
                    sflx(i)=sflx(i)+dqdt_aqh2so4(i,k,l)*adv_mass(l)/mbar(i,k)  &
                                   *pdel(i,k)/gravit  
                 enddo
              enddo
              call outfld( trim(cnst_name_cw(m))//'AQH2SO4', sflx(:ncol), ncol, lchnk)
            endif
        end do

        sflx(:)=0._r8
        do k=1,pver
           do i=1,ncol
              sflx(i)=sflx(i)+dqdt_aqhprxn(i,k)*specmw_so4_amode/mbar(i,k)  &
                             *pdel(i,k)/gravit  
           enddo
        enddo
        call outfld( 'AQSO4_H2O2', sflx(:ncol), ncol, lchnk)
        sflx(:)=0._r8
        do k=1,pver
           do i=1,ncol
              sflx(i)=sflx(i)+dqdt_aqo3rxn(i,k)*specmw_so4_amode/mbar(i,k)  &
                             *pdel(i,k)/gravit  
           enddo
        enddo
        call outfld( 'AQSO4_O3', sflx(:ncol), ncol, lchnk)

        xphlwc(:,:) = 0._r8
        do k = 1,pver
           do i = 1,ncol
              if (cldfrc(i,k) >=  1.0e-5_r8) then
              if( lwc(i,k) >= 1.0e-8_r8 ) then
                xphlwc(i,k) = -1._r8*log10(xph(i,k)) * lwc(i,k)
              endif
              endif
           enddo
        enddo
        call outfld( 'XPH_LWC', xphlwc(:ncol,:), ncol ,lchnk )

      end subroutine SETSOX

      end module MO_SETSOX
