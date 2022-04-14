      
      
      
      
      
      
      
      

      MODULE module_mosaic_gly

      IMPLICIT NONE

      INTEGER, PARAMETER :: nspecs   = 13,                             &
                            igly_g   =  1,                             &
                            igly_r1  =  2,                             &
                            igly_r2  =  3,                             &
                            igly_nh4 =  4,                             &
                            igly_sfc =  5,                             &
                            igly_oh  =  6,                             &
                            ic_as    =  7,                             &
                            ic_an    =  8,                             &
                            ia_nh4   =  9,                             &
                            ioh_g    = 10,                             &
                            iph      = 11,                             &
                            iwater   = 12,                             &
                            iarea    = 13

      
      LOGICAL, PARAMETER :: lfast_tau1 = .FALSE.
      

      CONTAINS

      
      
      SUBROUTINE rk4(y, dydx, n, x, h, yout, derivs)

      INTEGER n
      REAL(kind=8) :: h, x, dydx(n), y(n), yout(n)
      EXTERNAL derivs
      INTEGER i
      REAL(kind=8) :: h6, hh, xh, dym(nspecs), dyt(nspecs), yt(nspecs)

      hh=h*0.5
      h6=h/6.
      xh=x+hh

      DO i=1, n
        yt(i) = y(i) + hh * dydx(i)
      ENDDO

      CALL derivs(xh, yt, h, dyt)

      DO i=1, n
        yt(i) = y(i) + hh * dyt(i)
      ENDDO

      CALL derivs(xh, yt, h, dym)

      DO i=1, n
        yt(i) = y(i) + h * dym(i)
        dym(i) = dyt(i) + dym(i)
      ENDDO

      CALL derivs(x+h, yt, h, dyt)

      DO i=1, n
        yout(i) = y(i) + h6 * ( dydx(i) + dyt(i) + 2. * dym(i))
      ENDDO

      RETURN

      END SUBROUTINE rk4

      
      
      SUBROUTINE glysoa_simple(dtchem)

      USE module_data_mosaic_therm, ONLY: t_k, area_wet_a, gas, aer,   &
                                          jtotal, igly, iglysoa_sfc_a, nbin_a

      IMPLICIT NONE

      REAL(kind=8), INTENT(IN)  :: dtchem
      REAL(kind=8)              :: omega, gamma_gly, A, delta_gly, frac_A
      INTEGER                   :: ibin

      
      omega = 1.455e4 * sqrt(t_k / 58.0_8)

      
      
      
      
      
      
      
      gamma_gly = 3.3E-3

      
      A = 0.0
      DO ibin = 1, nbin_a
        A = A + area_wet_a(ibin)
      ENDDO

      
      IF (A > 0.0) THEN
        
        
        delta_gly = 0.25 * gamma_gly * A * omega * gas(igly) * dtchem

        
        delta_gly = MIN(gas(igly), delta_gly)

        
        gas(igly) = gas(igly) - delta_gly

        
        DO ibin = 1, nbin_a
          frac_A = area_wet_a(ibin) / A
          
          aer(iglysoa_sfc_a, jtotal, ibin) = aer(iglysoa_sfc_a, jtotal, ibin) &
                                       + frac_A * delta_gly
        ENDDO
      ENDIF

      END SUBROUTINE glysoa_simple

      SUBROUTINE glysoa_complex_derivs(x, y, dt, dydx)

      USE module_data_mosaic_therm, ONLY: conv1a,  & 
                                          p_atm,   & 
                                          t_k        

      REAL(kind=8), INTENT(IN)  :: x, y(nspecs), dt
      REAL(kind=8), INTENT(OUT) :: dydx(nspecs)

      REAL(kind=8), PARAMETER   :: eps      = 1.e-16 , & 
                                   Kh_water = 4.19e5 , & 
                                   Kh_oh    = 25.0,    & 
                                   k_oh     = 1.1e9      

      REAL(kind=8)              :: gly_g_atm,  & 
                                   f_A1,       & 
                                   tau1,       & 
                                   tau2,       & 
                                   oh_g_atm,   & 
                                   oh_a,       & 
                                   c_tot,      & 
                                   Kh_eq,      & 
                                   gly_ptot_eq,& 
                                   gly_r1_eq,  & 
                                   gly_r2_eq,  & 
                                   anh4,       & 
                                   kII, kI,    & 
                                   omega         

      
      REAL(kind=8)              :: dg_r1,      & 
                                   dr1_r2,     & 
                                   dr1_nh4,    & 
                                   dr1_oh,     & 
                                   dg_sfc        

      REAL(kind=8)              :: accloss,    & 
                                   scaling       

      dg_r1   = 0.0
      dr1_r2  = 0.0
      dr1_nh4 = 0.0
      dr1_oh  = 0.0
      dg_sfc  = 0.0

      
      

      gly_g_atm = y(igly_g) / conv1a 
      gly_g_atm = gly_g_atm * p_atm  

      
      f_A1 = 0.5
      tau1 = 2.5e2 
      tau2 = 5.5e3 
      IF ( y(ic_as) + y(ic_an) .GT. 12.0 ) THEN
        
        f_A1 = 0.6667
        tau1 = 4.4e4 
        tau2 = 4.7e4 
      ENDIF

      
      
      c_tot = MIN( 12.0, y(ic_as) + y(ic_an) )


      
      
      Kh_eq = Kh_water / 10**(-0.24D0 * c_tot) 

      
      
      gly_ptot_eq = gly_g_atm * Kh_eq * y(iwater) * 1e9 

      gly_r1_eq   = gly_ptot_eq *        f_A1  
      gly_r2_eq   = gly_ptot_eq * (1.0 - f_A1) 

      
      
      
      IF (.NOT. lfast_tau1) THEN
      
        dg_r1       = (1.0/tau1) * (gly_r1_eq - y(igly_r1))
      
      ENDIF
      
      
      dr1_r2      = (1.0/tau2) * (gly_r2_eq - y(igly_r2))

      
      

      oh_g_atm = y(ioh_g) / conv1a                  
      oh_g_atm = oh_g_atm * p_atm                   
      oh_a     = oh_g_atm * Kh_oh * y(iwater) * 1e9 

      dr1_oh   = k_oh * y(igly_r1) * oh_a           

      
      

      
      anh4 = MAX( 0.0, MIN( 4.0, y(ia_nh4)) ) 
      kII   = 2.e-10 * exp(1.5 * anh4) * exp(2.5  * y(iph))

      
      kI = kII * y(igly_r1) / y(iwater) * 1e-9 

      dr1_nh4 = kI * y(igly_r1) 

      
      

      
      omega = 1.455e4 * sqrt(t_k / 58.0D0)

      
      
      
      dg_sfc = 0.25D0 * 1.e-3 * y(iarea) * omega * y(igly_g)

      
      

      

      IF ( y(igly_g) < eps ) THEN
        dg_r1  = 0.0
        dg_sfc = 0.0
      ELSE
        accloss = (dg_r1 + dg_sfc) * dt
        IF ( ( y(igly_g) - accloss ) < eps ) THEN
          scaling = y(igly_g) / (accloss + eps)
          dg_r1  = dg_r1  * scaling
          dg_sfc = dg_sfc  * scaling
        ENDIF
      ENDIF

      IF ( y(igly_r1) < eps ) THEN
        dr1_r2  = 0.0
        dr1_nh4 = 0.0
        dr1_oh  = 0.0
      ELSE
        accloss = (-dg_r1 + dr1_r2 + dr1_nh4 + dr1_oh) * dt
        IF ( ( y(igly_r1) - accloss ) < eps) THEN
          scaling = y(igly_r1) / (accloss + eps)
          dr1_r2  = dr1_r2  * scaling
          dr1_nh4 = dr1_nh4 * scaling
          dr1_oh  = dr1_oh  * scaling
        ENDIF
      ENDIF

      IF ( y(igly_r2) < eps ) THEN
        dr1_r2  = MAX( 0.0, dr1_r2 )
      ELSE
        accloss = -dr1_r2 * dt
        IF ( ( y(igly_r2) - accloss ) < eps ) THEN
          scaling = y(igly_r2) / (accloss + eps)
          dr1_r2  = dr1_r2  * scaling
        ENDIF
      ENDIF

      
      dydx(igly_g)   = -dg_r1                          -dg_sfc
      dydx(igly_r1)  =  dg_r1 -dr1_r2 -dr1_nh4 -dr1_oh
      dydx(igly_r2)  =         dr1_r2
      dydx(igly_nh4) =                +dr1_nh4
      dydx(igly_oh)  =                          dr1_oh
      dydx(igly_sfc) =                                  dg_sfc

      END SUBROUTINE glysoa_complex_derivs

      SUBROUTINE glysoa_complex(dtchem)

      USE module_data_mosaic_therm, ONLY : jaerosolstate, all_liquid, mixed, &
                                           jtotal, jliquid, nbin_a, &
                                           area_wet_a, gas, water_a, aer, mc, &
                                           ph, a_nh4, c_as, c_an, &
                                           igly, iho, &
                                           iglysoa_r1_a, iglysoa_r2_a, &
                                           iglysoa_oh_a, &
                                           iglysoa_nh4_a, iglysoa_sfc_a, &
                                           iso4_a, ino3_a, inh4_a, jc_h

      
      USE module_data_mosaic_therm, ONLY: conv1a,  & 
                                          p_atm      
      

      REAL(kind=8), INTENT(IN)  :: dtchem

      REAL(kind=8)              :: A, conv, y(nspecs), yout(nspecs), &
                                   dydx(nspecs), gly_g

      INTEGER                   :: i, ii, nbin_proc, bin_proc(nbin_a)

      
      REAL(kind=8), PARAMETER   :: Kh_water = 4.19e5     

      REAL(kind=8)              :: gly_g_atm,  & 
                                   f_A1,       & 
                                   c_tot,      & 
                                   Kh_eq,      & 
                                   gly_ptot_eq,& 
                                   gly_r1_eq,  & 
                                   deltagly      
      


      
      
      nbin_proc   = 0
      bin_proc(:) = -1 
                       
      DO i = 1, nbin_a
        IF (jaerosolstate(i) == all_liquid .OR. &
            jaerosolstate(i) == mixed) THEN
          nbin_proc = nbin_proc + 1
          bin_proc(nbin_proc) = i
        ENDIF
      ENDDO
      IF (nbin_proc == 0) RETURN

      
      

      A     =     0.0 
      DO i = 1, nbin_proc
        ii = bin_proc(i)
        A = A + area_wet_a(ii)
      ENDDO
      IF (A <= 0) RETURN

      
      

      ph(:) = -9999.0 
      a_nh4(:) = 0.0  
      c_as(:) = 0.0   
      c_an(:) = 0.0   

      
      

      
      gly_g     = gas(igly)         

      DO i = 1, nbin_proc

        ii = bin_proc(i)

        
        

        conv        = 1.e-9 / water_a(ii) 

        y(:)        = 0.0

        
        y(igly_g)   = gly_g
        y(igly_r1)  = aer(iglysoa_r1_a,jtotal,ii)
        y(igly_r2)  = aer(iglysoa_r2_a,jtotal,ii)
        y(igly_nh4) = aer(iglysoa_nh4_a,jtotal,ii)
        y(igly_sfc) = aer(iglysoa_sfc_a,jtotal,ii)
        y(igly_oh)  = aer(iglysoa_oh_a,jtotal,ii)

        y(ic_as)    = aer(iso4_a,jliquid,ii) * conv 
        y(ic_an)    = aer(ino3_a,jliquid,ii) * conv 
        y(ia_nh4)   = aer(inh4_a,jliquid,ii) * conv 

        y(iph)      = MIN(14.0D0, MAX(0.0D0, -log10(mc(jc_h,ii)) ))

        y(iwater)   = water_a(ii) 

        y(ioh_g)    = gas(iho) 

        y(iarea)    = area_wet_a(ii) 

        
        IF (lfast_tau1) THEN

          
          
          gly_g_atm = y(igly_g) / conv1a 
          gly_g_atm = gly_g_atm * p_atm  

          
          f_A1 = 0.5
          IF ( y(ic_as) + y(ic_an) .GT. 12.0 ) THEN
            f_A1 = 0.6667
          ENDIF

          
          
          c_tot = MIN( 12.0, y(ic_as) + y(ic_an) )


          
          
          Kh_eq = Kh_water / 10**(-0.24D0 * c_tot) 

          
          
          gly_ptot_eq = gly_g_atm * Kh_eq * y(iwater) * 1e9 

          gly_r1_eq   = gly_ptot_eq *        f_A1  

          deltagly    = gly_r1_eq - y(igly_r1)

          y(igly_g)   = y(igly_g) - deltagly
          y(igly_r1)  = y(igly_r1) + deltagly

          IF (y(igly_g) < 0.0 .OR. y(igly_r1) < 0.0) THEN
            WRITE(*,*) "THIS IS NOT RIGHT: ",y(igly_g), y(igly_r1),deltagly
          ENDIF
        ENDIF
        

        
        
        CALL glysoa_complex_derivs(1._8, y, dtchem, dydx)
        CALL rk4(y, dydx, nspecs, 1._8, dtchem, yout, glysoa_complex_derivs)

        
        aer(iglysoa_r1_a,jtotal,ii)  = yout(igly_r1)
        aer(iglysoa_r2_a,jtotal,ii)  = yout(igly_r2)
        aer(iglysoa_nh4_a,jtotal,ii) = yout(igly_nh4)
        aer(iglysoa_sfc_a,jtotal,ii) = yout(igly_sfc)
        aer(iglysoa_oh_a,jtotal,ii)  = yout(igly_oh)

        
        
        gly_g                        = yout(igly_g)

        
        c_as(ii)  = yout(ic_as)
        c_an(ii)  = yout(ic_an)
        ph(ii)    = yout(iph)
        a_nh4(ii) = yout(ia_nh4)

      ENDDO

      
      gas(igly) = gly_g

      END SUBROUTINE glysoa_complex


      END MODULE module_mosaic_gly

