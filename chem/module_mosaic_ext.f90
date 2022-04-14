module module_mosaic_ext
contains


      
      subroutine dumpxx( dumpch, dtchem, t_k, p_atm, ah2o, &
         jaerosolstate, jphase, jhyst_leg, &
         aer, gas, num_a, water_a, water_a_hyst, dp_dry_a, &
         mosaic_vars_aa )

      use module_data_mosaic_kind,  only: r8

      use module_data_mosaic_aero,  only: &
         iso4_a,       ino3_a,       icl_a,        inh4_a,       &
         ina_a,        ioin_a,       ioc_a,        ibc_a,        &
         ipcg2_b_c_a,  ipcg2_b_o_a,  ipcg1_b_c_a,  ipcg1_b_o_a,  &
         iopcg1_b_c_a, iopcg1_b_o_a, ipcg2_f_c_a,  ipcg2_f_o_a,  &
         ipcg1_f_c_a,  ipcg1_f_o_a,  iopcg1_f_c_a, iopcg1_f_o_a, &
         iant1_c_a,    iant1_o_a,    ibiog1_c_a,   ibiog1_o_a,   &
         aer_name, gas_name, &
         jtotal, naer, nbin_a, nbin_a_max, ngas_aerchtot, ngas_volatile, &
         mosaic_vars_aa_type

      implicit none

      integer, intent(in), dimension(nbin_a_max) :: jaerosolstate, jphase, jhyst_leg

      real(r8), intent(in) :: dtchem, t_k, p_atm, ah2o
      real(r8), intent(in), dimension(naer,3,nbin_a_max) :: aer
      real(r8), intent(in), dimension(ngas_volatile) :: gas
      real(r8), intent(in), dimension(nbin_a_max) :: num_a, water_a, water_a_hyst, dp_dry_a

      character(len=*), intent(in) :: dumpch

      type (mosaic_vars_aa_type), intent(in) :: mosaic_vars_aa

      integer, parameter :: lunaa = 131
      integer, parameter :: naer_dump = 24
      integer :: ii, iv, ix, jj, kk, ktau
      integer :: iaer_dump(naer_dump) = 24
      character(len=80)  :: fmtaa
      character(len=2)   :: tmpch2
      character(len=256) :: tmpchaa
      character(len= 32) :: tmpchbb




      if ( mosaic_vars_aa%idiagbb_host < 100 ) return

      ii = mosaic_vars_aa%hostgridinfo(2)
      jj = mosaic_vars_aa%hostgridinfo(3)
      kk = mosaic_vars_aa%hostgridinfo(4)
      if ( ii*jj*kk /= 1 ) return

      ktau = mosaic_vars_aa%it_host
      iaer_dump = &
         (/ iso4_a,       ino3_a,       icl_a,        inh4_a,      &
            ina_a,        ioin_a,       ioc_a,        ibc_a,       &
            ipcg2_b_c_a,  ipcg2_b_o_a,  ipcg1_b_c_a,  ipcg1_b_o_a, &
            iopcg1_b_c_a, iopcg1_b_o_a, ipcg2_f_c_a,  ipcg2_f_o_a, &
            ipcg1_f_c_a,  ipcg1_f_o_a,  iopcg1_f_c_a, iopcg1_f_o_a, &
            iant1_c_a,    iant1_o_a,    ibiog1_c_a,   ibiog1_o_a /)



      tmpch2 = dumpch
      write(lunaa,'(/2a,i5,3i3)') tmpch2, 'dump', &
         ktau, ii, jj, kk
      write(lunaa,'(a,f11.2,f11.4,f11.4      )') 't p a ', t_k, p_atm, ah2o





      fmtaa = '(a,8i11)'
      if (nbin_a > 8) fmtaa = '(a,8i11/(10x,8i11))'
      write(lunaa,fmtaa) 'jstate    ', jaerosolstate(1:nbin_a)
      write(lunaa,fmtaa) 'jphase    ', jphase(1:nbin_a)
      write(lunaa,fmtaa) 'jhyst     ', jhyst_leg(1:nbin_a)

      fmtaa = '(a,1p,8e11.3)'
      if (nbin_a > 8) fmtaa = '(a,1p,8e11.3/(10x,1p,8e11.3))'
      write(lunaa,fmtaa) 'num       ', num_a(1:nbin_a)
      write(lunaa,fmtaa) 'dpdry     ', dp_dry_a(1:nbin_a)
      write(lunaa,fmtaa) 'water     ', water_a(1:nbin_a)
      write(lunaa,fmtaa) 'hyswtr    ', water_a_hyst(1:nbin_a)












       


      do ix = 1, naer_dump
         iv = iaer_dump(ix)
         if (iv < 1 .or. iv > naer) cycle
         fmtaa = '(a,1p,8e11.3)'
         write(tmpchaa,fmtaa) aer_name(iv)(1:10), aer(iv,jtotal,1:min(nbin_a,8))
         if (iv <= ngas_aerchtot) then
            write(tmpchbb,fmtaa) gas_name(iv)(1:10), gas(iv)
            tmpchaa = trim(tmpchaa) // '  ' // trim(tmpchbb)
         end if
         write(lunaa,'(a)') trim(tmpchaa)
         fmtaa = '(10x,1p,8e11.3)'
         if (nbin_a > 8) write(lunaa,fmtaa) aer(iv,jtotal,9:nbin_a)
      end do 


      return
      end subroutine dumpxx





  
  
  
  
  
  
  subroutine aerosol_phase_state( ibin, jaerosolstate, jphase, aer,  &
       jhyst_leg, electrolyte, epercent, kel, activity, mc, num_a, mass_wet_a, mass_dry_a,    &
       mass_soluble_a, vol_dry_a, vol_wet_a, water_a, water_a_hyst, water_a_up, aH2O_a,     &
       aH2O, ma, gam, log_gamZ, zc, za, gam_ratio,     &
       xeq_a, na_Ma, nc_Mc, xeq_c,       mw_electrolyte, partial_molar_vol, sigma_soln, T_K, & 
       RH_pc, mw_aer_mac, dens_aer_mac, sigma_water, Keq_ll, Keq_sl, MW_a, MW_c,             &
       growth_factor, MDRH, MDRH_T, molality0, rtol_mesa, jsalt_present, jsalt_index,       &
       jsulf_poor, jsulf_rich, phi_salt_old,                                     &
       kappa_nonelectro, mosaic_vars_aa )

    use module_data_mosaic_aero,  only: r8, nbin_a_max,                                   &
         ngas_aerchtot, ngas_volatile, nelectrolyte,                                      &
         Ncation, naer, jtotal, all_solid, jhyst_up, all_liquid, Nanion, nrxn_aer_ll,     &
         nrxn_aer_sl, nsalt, MDRH_T_NUM, jsulf_poor_NUM, jsulf_rich_NUM,                  &
         inh4_a, ina_a, ica_a, ico3_a, imsa_a, icl_a, ino3_a, iso4_a,                     & 
         a_zsr,  b_zsr,  aw_min,                                                          &
         mosaic_vars_aa_type


    implicit none
    

    integer, intent(in):: ibin
    integer, intent(in), dimension(nsalt) :: jsalt_index
    integer, intent(in), dimension(jsulf_poor_NUM) :: jsulf_poor
    integer, intent(in), dimension(jsulf_rich_NUM) :: jsulf_rich

    real(r8), intent(in) :: aH2O,T_K,RH_pc,rtol_mesa
    real(r8), intent(in), dimension(naer) :: mw_aer_mac,dens_aer_mac
    real(r8), intent(in), dimension(Ncation) :: zc,MW_c
    real(r8), intent(in), dimension(Nanion)  :: za,MW_a
    real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
    real(r8), intent(in), dimension(ngas_aerchtot) :: partial_molar_vol
    real(r8), intent(in), dimension(naer) :: kappa_nonelectro

    
    integer, intent(inout), dimension(nsalt) :: jsalt_present
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase,jhyst_leg

    real(r8), intent(inout) :: sigma_water
    real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
    real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a,mass_wet_a,mass_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_soluble_a,gam_ratio
    real(r8), intent(inout), dimension(nbin_a_max) :: vol_dry_a,vol_wet_a,water_a
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a_hyst,water_a_up,aH2O_a
    real(r8), intent(inout), dimension(nbin_a_max) :: sigma_soln,growth_factor,MDRH
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0            
    real(r8), intent(inout), dimension (nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension (nrxn_aer_sl) :: Keq_sl
    real(r8), intent(inout), dimension(MDRH_T_NUM) :: MDRH_T
    real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kel
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
    real(r8), intent(inout), dimension(nsalt) :: phi_salt_old

    type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

    
    character(len=256) :: errmsg
    integer, parameter :: aer_pha_sta_diagaa = -1 
    integer, parameter :: iter_kelvin_method =  3 
    
    
    
    integer, parameter :: iter_kelvin_meth1_max = 10
    integer, parameter :: iter_kelvin_meth2_max = 100
    integer :: iaer, iv, itmpa
    integer :: iter_kelvin, iter_kelvin_meth1, iter_kelvin_state
    integer :: js, je

    real(r8) :: aer_H
    real(r8):: aH2O_range_bisect_toler
    real(r8) :: aH2O_a_new, aH2O_a_old, aH2O_a_oldn, aH2O_a_oldp, aH2O_a_del_state3
    real(r8), dimension(nbin_a_max) :: DpmV
    real(r8), dimension(nbin_a_max) :: kelvin
    real(r8) :: kelvin_old, kelvin_oldn, kelvin_oldp
    real(r8) :: kelvin_toler
    real(r8) :: rel_err, rel_err_old, rel_err_old2, rel_err_oldn, rel_err_oldp
    real(r8) :: term, tmpa
    real(r8) :: water_a_old, water_a_oldn, water_a_oldp


    if (aer_pha_sta_diagaa >= 3) &
    write(*,'(/a,5i5,2f12.8,1p,2e11.3)') 'aer_pha_sta_a', ibin, jhyst_leg(ibin), jaerosolstate(ibin), -1, 0, aH2O, aH2O_a(ibin)
    
    aH2O_a(ibin) = aH2O
    kelvin(ibin) = 1.0
    do iv = 1, ngas_aerchtot
       kel(iv,ibin) = 1.0
    enddo







    kelvin_toler = 1.e-6_r8 * max( 1.0_r8-aH2O, 1.0e-4_r8 )
    aH2O_range_bisect_toler = 1.e-6_r8 * max( 1.0_r8-aH2O, 1.0e-4_r8 )


    
    mass_dry_a(ibin) = 0.0          
    vol_dry_a(ibin)  = 0.0          

    aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
         aer(ino3_a,jtotal,ibin) +   &
         aer(icl_a,jtotal,ibin)  +   &
         aer(imsa_a,jtotal,ibin) +   &
         2.*aer(ico3_a,jtotal,ibin))-   &
         (2.*aer(ica_a,jtotal,ibin)  +   &
         aer(ina_a,jtotal,ibin)  +   &
         aer(inh4_a,jtotal,ibin))
    aer_H = max(aer_H, 0.0d0)

    do iaer = 1, naer
       mass_dry_a(ibin) = mass_dry_a(ibin) +   &
            aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)  
       vol_dry_a(ibin)  = vol_dry_a(ibin) +   &
            aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)       
    enddo
    mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
    vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

    mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15                      
    vol_dry_a(ibin)  = vol_dry_a(ibin)*1.e-15                               

    
    mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3               
    vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3                


    water_a_up(ibin) = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)   

    iter_kelvin = 0
    iter_kelvin_meth1 = 0

    iter_kelvin_state = 0
    if (iter_kelvin_method == 2) iter_kelvin_state = 2

    aH2O_a_old = aH2O
    kelvin_old = 1.0_r8
    rel_err_old = 1.0e30_r8
    rel_err_old2 = 1.0e30_r8
    water_a_old = 0.0_r8

    aH2O_a_del_state3 = 1.0e-3_r8
    aH2O_a_oldn = aH2O
    aH2O_a_oldp = aH2O
    kelvin_oldp = 1.0_r8
    kelvin_oldn = 1.0_r8
    rel_err_oldn = 1.0e30_r8
    rel_err_oldp = 1.0e30_r8
    water_a_oldp = 0.0_r8
    water_a_oldn = 0.0_r8
    aH2O_a_new = aH2O    


10  iter_kelvin = iter_kelvin + 1
    aH2O_a(ibin) = aH2O_a_new


      do je = 1, nelectrolyte
        molality0(je,ibin) = bin_molality(je,ibin,aH2O_a,b_zsr,a_zsr,aw_min)  
      enddo
    call MESA( ibin, jaerosolstate, jphase, aer, jhyst_leg,         &
         electrolyte, epercent, activity, mc, num_a, mass_wet_a, mass_dry_a,              &
         mass_soluble_a, vol_dry_a, vol_wet_a, water_a, water_a_hyst, water_a_up, aH2O_a, &
         aH2O, ma, gam, log_gamZ, zc, za, gam_ratio, &
         xeq_a, na_Ma, nc_Mc, xeq_c, mw_electrolyte, mw_aer_mac, dens_aer_mac, Keq_ll,     &
         Keq_sl, MW_c, MW_a, growth_factor, MDRH, MDRH_T, molality0, rtol_mesa,            &
         jsalt_present, jsalt_index, jsulf_poor, jsulf_rich, phi_salt_old,    &
         kappa_nonelectro, mosaic_vars_aa )

    if(jaerosolstate(ibin) .eq. all_solid)then
       if (aer_pha_sta_diagaa >= 2) &
       write(*,'(a,5i5,2f12.8,1p,2e11.3)') 'aer_pha_sta_b', ibin, jhyst_leg(ibin), jaerosolstate(ibin), &
          iter_kelvin_state, iter_kelvin, aH2O, aH2O_a(ibin)
       return
    endif
    
    mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3               
    vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3                
 
    call calculate_kelvin(ibin,num_a,vol_wet_a,aH2O_a,DpmV,kelvin,sigma_soln,T_K,  &
         sigma_water)
    
    kelvin(ibin) = max( kelvin(ibin), 1.0_r8 )
    if (water_a(ibin) <= 0.0_r8) kelvin(ibin) = 1.0_r8

    aH2O_a_new = aH2O/kelvin(ibin)













    rel_err = (aH2O_a(ibin)*kelvin(ibin) - aH2O) / max( aH2O, 0.01_r8 )

    if (aer_pha_sta_diagaa >= 10) &
    write(*,'(a,2i5, 1p,e10.2, 0p,f14.10, 2x,2f14.10, 2x,1p,2e18.10)') &
       'iter_kelvin', iter_kelvin_state, iter_kelvin, rel_err, kelvin(ibin), &
       aH2O_a(ibin), aH2O_a_new, water_a_old, water_a(ibin)

    if (abs(rel_err) <= kelvin_toler) then
       iter_kelvin_state = iter_kelvin_state + 100
       goto 90
    end if

    if (iter_kelvin_state <= 0) then
       
       itmpa = 0
       if (iter_kelvin >= iter_kelvin_meth1_max) then
          itmpa = 1
       else if (iter_kelvin >= iter_kelvin_meth1_max) then
          tmpa = min( rel_err_old, rel_err_old2 )
          if (tmpa < 0.0_r8 .and. rel_err <= tmpa) itmpa = 1
          tmpa = max( rel_err_old, rel_err_old2 )
          if (tmpa > 0.0_r8 .and. rel_err >= tmpa) itmpa = 1
       end if

       if (itmpa > 0) then
          if (iter_kelvin_method <= 1) then
             
             
             
             aH2O_a(ibin) = aH2O_a_new   
             if (aer_pha_sta_diagaa >= 1) &
             write(*,'(a,5i5,2f12.8,1p,3e11.3)') 'iter_kelv_err', ibin, jhyst_leg(ibin), jaerosolstate(ibin), &
                iter_kelvin_state, iter_kelvin, aH2O, aH2O_a(ibin), rel_err, kelvin_toler
             iter_kelvin_state = 100
             goto 90
          else
             
             iter_kelvin_state = 1
             iter_kelvin_meth1 = iter_kelvin
          end if
       else
          
          aH2O_a_old = aH2O_a(ibin)
          kelvin_old = kelvin(ibin)
          rel_err_old2 = rel_err_old
          rel_err_old = rel_err
          water_a_old  = water_a(ibin)
       
       
          goto 10
       end if
    endif

    if (iter_kelvin_state == 1) then
       
       iter_kelvin_state = 2
       if (rel_err < 0.0_r8) then
          
          aH2O_a_new = aH2O
          goto 10
       else
          
          
          continue
       end if
    end if

    if (iter_kelvin_state == 2) then
       
       
       if (rel_err < 0.0_r8) then
          
          if (aer_pha_sta_diagaa >= 1) &
             write(*,'(a,5i5,2f12.8,1p,3e11.3)') 'iter_kelv_er2', ibin, jhyst_leg(ibin), jaerosolstate(ibin), &
                iter_kelvin_state, iter_kelvin, aH2O, aH2O_a(ibin), rel_err, kelvin_toler
          iter_kelvin_state = 100
          goto 90
       end if
       
       
       aH2O_a_oldp = aH2O_a(ibin)
       kelvin_oldp = kelvin(ibin)
       rel_err_oldp = rel_err
       water_a_oldp  = water_a(ibin)
       aH2O_a_new = min( aH2O/kelvin(ibin), 0.999999_r8 )   
       iter_kelvin_state = 3
       goto 10
    end if

    if (iter_kelvin_state == 3) then
       
       
       if (rel_err < 0.0_r8) then
          
          
          aH2O_a_oldn = aH2O_a(ibin)
          kelvin_oldn = kelvin(ibin)
          rel_err_oldn = rel_err
          water_a_oldn  = water_a(ibin)
          aH2O_a_new = 0.5_r8*(aH2O_a_oldn + aH2O_a_oldp)
          iter_kelvin_state = 4
          goto 10
       else
          
          if ( (rel_err >= rel_err_oldp) .or. &
               (aH2O_a_del_state3 >= 0.999_r8) ) then
             
             if (aer_pha_sta_diagaa >= 1) &
                write(*,'(a,5i5,2f12.8,1p,3e11.3)') 'iter_kelv_er3', ibin, jhyst_leg(ibin), jaerosolstate(ibin), &
                   iter_kelvin_state, iter_kelvin, aH2O, aH2O_a(ibin), rel_err, kelvin_toler
             iter_kelvin_state = 200
             goto 90
          else
             
             
             
             aH2O_a_oldp = aH2O_a(ibin)
             kelvin_oldp = kelvin(ibin)
             rel_err_oldp = rel_err
             water_a_oldp  = water_a(ibin)
             aH2O_a_new = aH2O_a(ibin) - aH2O_a_del_state3
             aH2O_a_del_state3 = aH2O_a_del_state3*1.5_r8
             if (aH2O_a_new .le. 0.01_r8) then
                aH2O_a_new = 0.01_r8
                aH2O_a_del_state3 = 1.0_r8
             end if
             goto 10
          end if
       end if
    end if

    if (iter_kelvin_state == 4) then
       
       if ( iter_kelvin >= iter_kelvin_meth2_max + iter_kelvin_meth1 ) then
          
          if (aer_pha_sta_diagaa >= 1) &
             write(*,'(a,5i5,2f12.8,1p,3e11.3)') 'iter_kelv_er4', ibin, jhyst_leg(ibin), jaerosolstate(ibin), &
                iter_kelvin_state, iter_kelvin, aH2O, aH2O_a(ibin), rel_err, kelvin_toler
          iter_kelvin_state = 301
          goto 90
       else if ( abs(aH2O_a_oldp - aH2O_a_oldn) <= aH2O_range_bisect_toler ) then
          




          iter_kelvin_state = 302
          goto 90
       end if
       
       
       if (rel_err >= 0.0_r8) then
          if (rel_err >= rel_err_oldp) then
             
             
             if (aer_pha_sta_diagaa >= 1) &
                write(*,'(a,5i5,2f12.8,1p,3e11.3)') 'iter_kelv_er6', ibin, jhyst_leg(ibin), jaerosolstate(ibin), &
                   iter_kelvin_state, iter_kelvin, aH2O, aH2O_a(ibin), rel_err, kelvin_toler
             iter_kelvin_state = 303
             goto 90
          else
             
             aH2O_a_oldp = aH2O_a(ibin)
             kelvin_oldp = kelvin(ibin)
             rel_err_oldp = rel_err
             water_a_oldp  = water_a(ibin)
          end if
       else
          if (rel_err <= rel_err_oldn) then
             
             
             if (aer_pha_sta_diagaa >= 1) &
                write(*,'(a,5i5,2f12.8,1p,3e11.3)') 'iter_kelv_er7', ibin, jhyst_leg(ibin), jaerosolstate(ibin), &
                   iter_kelvin_state, iter_kelvin, aH2O, aH2O_a(ibin), rel_err, kelvin_toler
             iter_kelvin_state = 304
             goto 90
          else
             
             aH2O_a_oldn = aH2O_a(ibin)
             kelvin_oldn = kelvin(ibin)
             rel_err_oldn = rel_err
             water_a_oldn  = water_a(ibin)
          end if
       end if
       aH2O_a_new = 0.5_r8*(aH2O_a_oldn + aH2O_a_oldp)
       goto 10
    end if

    write(errmsg,'(a,4i5)') 'iter_kelv fatal err 1', ibin, iter_kelvin, iter_kelvin_state
    call wrf_error_fatal3("<stdin>",530,&
trim(adjustl(errmsg)))


    
90  if (iter_kelvin_state == 200) then
       
       if (abs(rel_err_oldp) < abs(rel_err)) then
          aH2O_a(ibin) = aH2O_a_oldp
          rel_err = rel_err_oldp
       end if
    else if (iter_kelvin_state >= 300 .and. iter_kelvin_state <= 304) then
       
       tmpa = min( abs(rel_err_oldn), abs(rel_err_oldp), abs(rel_err) )
       if (abs(rel_err_oldp) == tmpa) then
          aH2O_a(ibin) = aH2O_a_oldp
          rel_err = rel_err_oldp
       else if (abs(rel_err_oldn) == tmpa) then
          aH2O_a(ibin) = aH2O_a_oldn
          rel_err = rel_err_oldn
       end if
    end if

    if(jaerosolstate(ibin) .eq. all_liquid)jhyst_leg(ibin) = jhyst_up

    
    do iv = 1,  ngas_aerchtot
       term = 4.*sigma_soln(ibin)*partial_molar_vol(iv)/   &
            (8.3144e7*T_K*DpmV(ibin))
       kel(iv,ibin) = 1. + term*(1. + 0.5*term*(1. + term/3.))
    enddo

    if (aer_pha_sta_diagaa >= 2) &
    write(*,'(a,5i5,2f12.8,1p,e11.3,e14.5)') 'aer_pha_sta_c', ibin, jhyst_leg(ibin), jaerosolstate(ibin), &
       iter_kelvin_state, iter_kelvin, aH2O, aH2O_a(ibin), rel_err, water_a(ibin)
    return
  end subroutine aerosol_phase_state



  
  
  
  
  
  
  
  
  
  
  
  
  subroutine MESA( ibin, jaerosolstate, jphase, aer, jhyst_leg,      &
       electrolyte, epercent, activity, mc, num_a, mass_wet_a, mass_dry_a, mass_soluble_a,    &
       vol_dry_a, vol_wet_a, water_a, water_a_hyst, water_a_up, aH2O_a, aH2O,                &
       ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a,     &
       na_Ma, nc_Mc, xeq_c, mw_electrolyte, mw_aer_mac, dens_aer_mac, Keq_ll, Keq_sl, MW_c,    &
       MW_a, growth_factor, MDRH, MDRH_T, molality0, rtol_mesa, jsalt_present, jsalt_index,   &
       jsulf_poor, jsulf_rich, phi_salt_old,                                       &
       kappa_nonelectro, mosaic_vars_aa )

    use module_data_mosaic_aero,  only: r8, nbin_a_max, nelectrolyte, Ncation, naer,        &
         jtotal, all_solid, jsolid, all_liquid, jliquid, jhyst_lo, mhyst_uporlo_jhyst,       &
         jhyst_up, mhyst_uporlo_waterhyst, nsoluble, nsalt, Nanion, nrxn_aer_sl,            &
         nrxn_aer_ll, MDRH_T_NUM, jsulf_poor_NUM, jsulf_rich_NUM,                         &
         ptol_mol_astem,  mhyst_force_lo,  mhyst_force_up,                               &
         jcacl2, jcano3, mhyst_method, ioin_a, ibc_a, jcaco3, jcaso4,                     & 
         mosaic_vars_aa_type



    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(inout), dimension(nbin_a_max)  :: jhyst_leg
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase
    integer, intent(in), dimension(nsalt) :: jsalt_index
    integer, intent(inout), dimension(nsalt) :: jsalt_present
    integer, intent(in), dimension(jsulf_poor_NUM) :: jsulf_poor
    integer, intent(in), dimension(jsulf_rich_NUM) :: jsulf_rich

    real(r8), intent(in) :: aH2O,rtol_mesa
    real(r8), intent(in), dimension(naer) :: mw_aer_mac,dens_aer_mac
    real(r8), intent(in), dimension(Ncation) :: zc,MW_c
    real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
    real(r8), intent(in), dimension(Nanion)  :: za,MW_a
    real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a,mass_wet_a,mass_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_soluble_a,vol_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: vol_wet_a,gam_ratio
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a,water_a_hyst,water_a_up
    real(r8), intent(inout), dimension(nbin_a_max) :: aH2O_a,growth_factor,MDRH
    real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nrxn_aer_sl) :: Keq_sl
    real(r8), intent(inout), dimension(MDRH_T_NUM) :: MDRH_T
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
    real(r8), intent(inout), dimension(nsalt) :: phi_salt_old
    real(r8), intent(in), dimension(naer) :: kappa_nonelectro

    type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

    
    integer :: idissolved, j_index, jsalt_dum, jdum, js, je 
    real(r8) :: CRH, solids, sum_soluble, sum_insoluble, XT 
    
    
    real(r8) :: H_ion, sum_dum


    
    
    sum_dum = 0.0
    do je = 1, nelectrolyte
       sum_dum = sum_dum + electrolyte(je,jtotal,ibin)
    enddo

    if(sum_dum .eq. 0.)sum_dum = 1.0

    do je = 1, nelectrolyte
       epercent(je,jtotal,ibin) = 100.*electrolyte(je,jtotal,ibin)/sum_dum
    enddo


    call calculate_XT(ibin,jtotal,XT,aer)




    jsalt_dum = 0 
    do js = 1, nsalt
       jsalt_present(js) = 0                        

       if(epercent(js,jtotal,ibin) .gt. ptol_mol_astem)then
          jsalt_present(js) = 1                     
          jsalt_dum = jsalt_dum + jsalt_index(js)   
       endif
    enddo


    if( (epercent(jcano3,jtotal,ibin) .gt. ptol_mol_astem) .or. &
        (epercent(jcacl2,jtotal,ibin) .gt. ptol_mol_astem) )then
      CRH = 0.0  
    else
      CRH = 0.35 
    endif

  
    if(jsalt_dum .eq. 0)then			

       CRH = 0.0
       MDRH(ibin) = 0.0
       jaerosolstate(ibin) = all_solid
       jphase(ibin)    = jsolid
       jhyst_leg(ibin) = jhyst_lo
       call adjust_solid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent,water_a)
       water_a(ibin) = aerosol_water(jtotal,ibin,jaerosolstate,jphase,jhyst_leg,   &	
        electrolyte,aer,kappa_nonelectro,num_a,mass_dry_a,mass_soluble_a,aH2O,molality0)
       return

    elseif(XT .lt. 1. .and. XT .gt. 0.0)then  	
       MDRH(ibin) = 0.0
    elseif(XT .ge. 2.0 .or. XT .lt. 0.0)then  	
       j_index = jsulf_poor(jsalt_dum)          
       MDRH(ibin) = MDRH_T(j_index)
    else					
       j_index = jsulf_rich(jsalt_dum)          
       MDRH(ibin) = MDRH_T(j_index)
    endif

    CRH = min(CRH, MDRH(ibin)/100.0)		




    
    
    if( aH2O_a(ibin).lt.CRH )then
       jaerosolstate(ibin) = all_solid
       jphase(ibin)    = jsolid
       jhyst_leg(ibin) = jhyst_lo
       call adjust_solid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent,water_a)
       water_a(ibin) = aerosol_water(jtotal,ibin,jaerosolstate,jphase,jhyst_leg,   &	
        electrolyte,aer,kappa_nonelectro,num_a,mass_dry_a,mass_soluble_a,aH2O,molality0)
       return
    endif


    
    jdum = 0
    if (mhyst_method == mhyst_uporlo_jhyst) then         
       if (jhyst_leg(ibin) == jhyst_up) jdum = 1
    elseif (mhyst_method == mhyst_uporlo_waterhyst) then 
       if (water_a_hyst(ibin) > 0.5*water_a_up(ibin)) jdum = 1
       
    elseif (mhyst_method == mhyst_force_lo) then
       jdum = 0
    elseif (mhyst_method == mhyst_force_up) then
       jdum = 1
       
    else
       call wrf_error_fatal3("<stdin>",739,&
'*** MESA - bad mhyst_method')
    endif
    if (jdum == 1) then 
       call do_full_deliquescence(ibin,aer,electrolyte)

       
       
       
































          jaerosolstate(ibin) = all_liquid
          jhyst_leg(ibin) = jhyst_up
          jphase(ibin) = jliquid
          water_a(ibin) = aerosol_water(jtotal,ibin,jaerosolstate,jphase,jhyst_leg,   &
               electrolyte,aer,kappa_nonelectro,num_a,mass_dry_a,mass_soluble_a,aH2O,molality0)
          if(water_a(ibin) .le. 0.0)then     
             jdum = 0
             jaerosolstate(ibin) = all_solid 
             jphase(ibin)    = jsolid
             jhyst_leg(ibin) = jhyst_lo
             call adjust_solid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent,water_a)
          else
             call adjust_liquid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent)
             call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,      &
                  electrolyte,activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a, &
                  aH2O,ma,gam,log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)
          endif

          
          mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3 
          vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3  
          growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)   

          return




    endif 


    
    if(aH2O_a(ibin)*100. .lt. MDRH(ibin)) then
       jaerosolstate(ibin) = all_solid
       jphase(ibin) = jsolid
       jhyst_leg(ibin) = jhyst_lo
       call adjust_solid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent,water_a)
       return
    endif


    
10  call do_full_deliquescence(ibin,aer,electrolyte)
    call MESA_PTC( ibin, jaerosolstate, jphase, aer, jhyst_leg,                   &
         electrolyte, epercent, activity, mc, num_a, mass_dry_a, mass_wet_a,      &
         mass_soluble_a, vol_dry_a, vol_wet_a, water_a, aH2O,                     &
         ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c,        &
         mw_electrolyte, mw_aer_mac, dens_aer_mac, Keq_sl, MW_c, MW_a, Keq_ll,    &
         growth_factor, molality0, rtol_mesa, jsalt_present, phi_salt_old,        &
         kappa_nonelectro, mosaic_vars_aa                                    )     
    return
  end subroutine MESA



  
  
  
  
  
  
  subroutine calculate_kelvin(ibin,num_a,vol_wet_a,aH2O_a,DpmV,kelvin,sigma_soln,  &
       T_K,sigma_water)
    use module_data_mosaic_constants, only:  pi
    use module_data_mosaic_aero, only: r8,nbin_a_max                                   

    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(in) :: T_K,sigma_water
    real(r8), intent(in), dimension(nbin_a_max) :: num_a
    real(r8), intent(inout), dimension(nbin_a_max) :: sigma_soln
    real(r8), intent(inout), dimension(nbin_a_max) ::vol_wet_a,aH2O_a,DpmV,kelvin
    
    integer je
    real(r8) :: term, sum_dum
    real(r8), dimension(nbin_a_max) :: volume_a

    volume_a(ibin) = vol_wet_a(ibin)                                
    DpmV(ibin)=(6.*volume_a(ibin)/(num_a(ibin)*pi))**(1./3.)        


    
    
    
    
    
    
    


    
    sigma_soln(ibin) = sigma_water + 49.0*(1. - aH2O_a(ibin))       



    term = 72.*sigma_soln(ibin)/(8.3144e7*T_K*DpmV(ibin))           

    kelvin(ibin) = 1. + term*(1. + 0.5*term*(1. + term/3.))


    return
  end subroutine calculate_kelvin



  
  
  
  
  
  
  subroutine calculate_XT(ibin,jp,XT,aer)
    use module_data_mosaic_aero, only: r8,naer,nbin_a_max,                         &
         imsa_a,iso4_a,ica_a,ina_a,inh4_a

    implicit none

    
    integer, intent(in) :: ibin, jp
    real(r8), intent(inout) :: XT
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer


    if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
       XT   = ( aer(inh4_a,jp,ibin) +   &
            aer(ina_a,jp,ibin)  +   &
            2.*aer(ica_a,jp,ibin) )/   &
            (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
    else
       XT   = -1.0
    endif


    return
  end subroutine calculate_XT



  
  
  
  
  
  
  subroutine adjust_solid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent,  &
       water_a)

    use module_data_mosaic_aero, only: r8,nbin_a_max,naer,nelectrolyte,jsolid,     &
         jhyst_lo,jtotal,jliquid,                                                  &
         inh4_a,ino3_a,icl_a                                                        

    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(inout), dimension(nbin_a_max) :: jphase,jhyst_leg
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte,epercent
    
    integer iaer, je


    jphase(ibin)    = jsolid
    jhyst_leg(ibin) = jhyst_lo   
    water_a(ibin)   = 0.0

    
    do iaer = 1, naer
       aer(iaer, jsolid, ibin) = aer(iaer,jtotal,ibin)
       aer(iaer, jliquid,ibin) = 0.0
    enddo

    
    do je = 1, nelectrolyte
       electrolyte(je,jliquid,ibin) = 0.0
       epercent(je,jliquid,ibin)    = 0.0
       electrolyte(je,jsolid,ibin)  = electrolyte(je,jtotal,ibin)
       epercent(je,jsolid,ibin)     = epercent(je,jtotal,ibin)
    enddo

    
    aer(inh4_a,jtotal,ibin) = aer(inh4_a,jsolid,ibin)
    aer(ino3_a,jtotal,ibin) = aer(ino3_a,jsolid,ibin)
    aer(icl_a,jtotal,ibin)  = aer(icl_a,jsolid,ibin)


    return
  end subroutine adjust_solid_aerosol



  
  
  
  
  
  
  subroutine adjust_liquid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent) 

    use module_data_mosaic_aero, only: r8,nbin_a_max,naer,nelectrolyte,jliquid,    &
         jhyst_up,jsolid,jtotal,                                                   &
         jcaco3,jcaso4,                                                            &
         inh4_a,ina_a,ica_a,imsa_a,icl_a,ino3_a,iso4_a

    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(inout), dimension(nbin_a_max) :: jphase,jhyst_leg

    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
    
    integer :: iaer, je

    jphase(ibin)    = jliquid
    jhyst_leg(ibin) = jhyst_up   

    
    do je = 1, nelectrolyte
       electrolyte(je,jsolid,ibin)  = 0.0
       epercent(je,jsolid,ibin)     = 0.0
       electrolyte(je,jliquid,ibin) = electrolyte(je,jtotal,ibin)
       epercent(je,jliquid,ibin)    = epercent(je,jtotal,ibin)
    enddo
    
    electrolyte(jcaco3,jsolid,ibin) = electrolyte(jcaco3,jtotal,ibin)
    electrolyte(jcaso4,jsolid,ibin) = electrolyte(jcaso4,jtotal,ibin)
    epercent(jcaco3,jsolid,ibin)    = epercent(jcaco3,jtotal,ibin)
    epercent(jcaso4,jsolid,ibin)    = epercent(jcaso4,jtotal,ibin)
    electrolyte(jcaco3,jliquid,ibin)= 0.0
    electrolyte(jcaso4,jliquid,ibin)= 0.0
    epercent(jcaco3,jliquid,ibin)   = 0.0
    epercent(jcaso4,jliquid,ibin)   = 0.0


    
    
    do iaer = 1, naer
    aer(iaer,jsolid,ibin)  = aer(iaer,jtotal,ibin)
    end do
    aer(iso4_a,jsolid,ibin) = electrolyte(jcaso4,jsolid,ibin)
    aer(ino3_a,jsolid,ibin) = 0.0
    aer(icl_a,jsolid,ibin)  = 0.0
    aer(inh4_a,jsolid,ibin) = 0.0
    aer(imsa_a,jsolid,ibin) = 0.0
    aer(ina_a,jsolid,ibin)  = 0.0
    aer(ica_a,jsolid,ibin)  = electrolyte(jcaco3,jsolid,ibin) +   &
                              electrolyte(jcaso4,jsolid,ibin)













    
    do iaer = 1, naer
    aer(iaer,jliquid,ibin)  = 0.0
    end do
    aer(iso4_a,jliquid,ibin) = aer(iso4_a,jtotal,ibin) -   &
                               aer(iso4_a,jsolid,ibin)
    aer(iso4_a,jliquid,ibin) = max(0.d0, aer(iso4_a,jliquid,ibin)) 
    aer(ino3_a,jliquid,ibin) = aer(ino3_a,jtotal,ibin)
    aer(icl_a,jliquid,ibin)  = aer(icl_a,jtotal,ibin)
    aer(inh4_a,jliquid,ibin) = aer(inh4_a,jtotal,ibin)
    aer(imsa_a,jliquid,ibin) = aer(imsa_a,jtotal,ibin)
    aer(ina_a,jliquid,ibin)  = aer(ina_a,jtotal,ibin)
    aer(ica_a,jliquid,ibin)  = aer(ica_a,jtotal,ibin) -   &
                               aer(ica_a,jsolid,ibin)
    aer(ica_a,jliquid,ibin)  = max(0.d0, aer(ica_a,jliquid,ibin)) 













    return
  end subroutine adjust_liquid_aerosol



  
  
  
  
  
  
  
  
  
  
  subroutine do_full_deliquescence(ibin,aer,electrolyte)    
    use module_data_mosaic_aero, only: r8,naer,nbin_a_max,nelectrolyte,jtotal,jsolid, &
         jliquid,                                                                     &
         jcacl2,jcano3,ioin_a,jcaco3,jcaso4,                                          &
         inh4_a,ina_a,ica_a,imsa_a,icl_a,ino3_a,iso4_a



    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer ::  iaer, js

    
    do js = 1, nelectrolyte
       electrolyte(js,jsolid,ibin)  = 0.0
       electrolyte(js,jliquid,ibin) = electrolyte(js,jtotal,ibin)
    enddo
    
    
    electrolyte(jcaco3,jsolid,ibin) = electrolyte(jcaco3,jtotal,ibin)
    electrolyte(jcaso4,jsolid,ibin) = electrolyte(jcaso4,jtotal,ibin)
    electrolyte(jcaco3,jliquid,ibin)= 0.0
    electrolyte(jcaso4,jliquid,ibin)= 0.0


    
    
    do iaer = 1, naer
    aer(iaer,jsolid,ibin)= aer(iaer,jtotal,ibin)
    end do
    aer(iso4_a,jsolid,ibin) = electrolyte(jcaso4,jsolid,ibin)
    aer(ino3_a,jsolid,ibin) = 0.0
    aer(icl_a, jsolid,ibin) = 0.0
    aer(inh4_a,jsolid,ibin) = 0.0
    aer(imsa_a,jsolid,ibin) = 0.0
    aer(ina_a, jsolid,ibin) = 0.0
    aer(ica_a, jsolid,ibin) = electrolyte(jcaco3,jsolid,ibin) +   &
                              electrolyte(jcaso4,jsolid,ibin)













    
    do iaer = 1, naer
    aer(iaer,jliquid,ibin) = 0.0
    end do
    aer(iso4_a,jliquid,ibin) = max(0.0_r8, aer(iso4_a,jtotal,ibin) -   &
                               electrolyte(jcaso4,jsolid,ibin))      
    aer(ino3_a,jliquid,ibin) = aer(ino3_a,jtotal,ibin)
    aer(icl_a, jliquid,ibin) = aer(icl_a,jtotal,ibin)
    aer(inh4_a,jliquid,ibin) = aer(inh4_a,jtotal,ibin)
    aer(imsa_a,jliquid,ibin) = aer(imsa_a,jtotal,ibin)
    aer(ina_a, jliquid,ibin) = aer(ina_a,jtotal,ibin)
    aer(ica_a, jliquid,ibin) = electrolyte(jcano3,jtotal,ibin) +   &
                               electrolyte(jcacl2,jtotal,ibin)













    return
  end subroutine do_full_deliquescence
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine MESA_PTC(ibin, jaerosolstate, jphase, aer, jhyst_leg,  &
       electrolyte, epercent, activity, mc, num_a, mass_dry_a, mass_wet_a, mass_soluble_a,    &
       vol_dry_a, vol_wet_a, water_a, aH2O, ma, gam,   &
       log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c, mw_electrolyte, mw_aer_mac,     &
       dens_aer_mac, Keq_sl, MW_c, MW_a, Keq_ll, growth_factor, molality0, rtol_mesa,         &
       jsalt_present, phi_salt_old,                                                &
       kappa_nonelectro, mosaic_vars_aa )                

    use module_data_mosaic_aero,  only: r8, nbin_a_max, nelectrolyte, Ncation, naer, nsalt,  &
         jhyst_lo, mixed, all_liquid, jsolid, jliquid, jtotal, mYES,                         &
         all_solid, Nanion, nrxn_aer_sl, nrxn_aer_ll,                                     &
         ino3_a, iso4_a, ioc_a, ilim1_a, ilim2_a, inh4_a, ina_a, ica_a, ico3_a, imsa_a, icl_a, &
         mosaic_vars_aa_type

    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(inout), dimension(nsalt) :: jsalt_present
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase,jhyst_leg

    real(r8), intent(in) :: aH2O,rtol_mesa
    real(r8), intent(in), dimension(naer) :: mw_aer_mac,dens_aer_mac
    real(r8), intent(in), dimension(Ncation) :: zc,MW_c
    real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
    real(r8), intent(in), dimension(Nanion)  :: za,MW_a
    real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a,mass_dry_a,mass_wet_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_soluble_a,vol_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: growth_factor
    real(r8), intent(inout), dimension(nbin_a_max) :: vol_wet_a,water_a,gam_ratio
    real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
    real(r8), intent(inout), dimension(nrxn_aer_sl) :: Keq_sl
    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
    real(r8), intent(inout), dimension(nsalt) :: phi_salt_old
    real(r8), intent(in), dimension(naer) :: kappa_nonelectro

    type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

    
    integer iaer, iconverge, iconverge_flux, iconverge_mass,   &
         idissolved, itdum, js, je, jp

    real(r8) :: tau_p(nsalt), tau_d(nsalt)
    real(r8) :: frac_solid, sumflux, hsalt_min, alpha, XT, dumdum,   &
         H_ion
    real(r8) :: phi_prod, alpha_fac, sum_dum
    real(r8) :: aer_H,hsalt_max
    real(r8), dimension(nelectrolyte) :: eleliquid
    real(r8), dimension(nbin_a_max) :: mass_dry_salt
    real(r8), dimension(nsalt) :: phi_salt,flux_sl,phi_bar,alpha_salt
    real(r8), dimension(nsalt) :: sat_ratio,hsalt
  
    
    

    
    itdum = 0               
    hsalt_max = 1.e25



    do js = 1, nsalt
       hsalt(js)     = 0.0
       sat_ratio(js) = 0.0
       phi_salt(js)  = 0.0
       flux_sl(js)   = 0.0
    enddo



    
    sum_dum = 0.0
    do je = 1, nelectrolyte
       sum_dum = sum_dum + electrolyte(je,jtotal,ibin)
    enddo

    if(sum_dum .eq. 0.)sum_dum = 1.0

    do je = 1, nelectrolyte
       epercent(je,jtotal,ibin) = 100.*electrolyte(je,jtotal,ibin)/sum_dum
    enddo
    



    do js = 1, nsalt
       jsalt_present(js) = 0                        
       if(epercent(js,jtotal,ibin) .gt. 1.0)then
          jsalt_present(js) = 1                     
       endif
    enddo


    mass_dry_a(ibin) = 0.0

    aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
         aer(ino3_a,jtotal,ibin) +   &
         aer(icl_a,jtotal,ibin)  +   &
         aer(imsa_a,jtotal,ibin) +   &
         2.*aer(ico3_a,jtotal,ibin))-   &
         (2.*aer(ica_a,jtotal,ibin)  +   &
         aer(ina_a,jtotal,ibin)  +   &
         aer(inh4_a,jtotal,ibin))
    aer_H = max(aer_H, 0.0d0)

    do iaer = 1, naer
       mass_dry_a(ibin) = mass_dry_a(ibin) +   &
            aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)  
       vol_dry_a(ibin)  = vol_dry_a(ibin) +   &
            aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)       
    enddo
    mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
    vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

    mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15                      
    vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15                                

    mass_dry_salt(ibin) = 0.0               
    do je = 1, nsalt
       mass_dry_salt(ibin) = mass_dry_salt(ibin) +   &
            electrolyte(je,jtotal,ibin)*mw_electrolyte(je)*1.e-15   
    enddo

    mosaic_vars_aa%jMESA_call = mosaic_vars_aa%jMESA_call + 1
    
    

    do 500 itdum = 1, mosaic_vars_aa%Nmax_MESA
       
       
       
       call MESA_flux_salt(ibin,jaerosolstate,jphase, aer,jhyst_leg,electrolyte, &
            epercent,activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,ma,&
            gam,log_gamZ,zc,za,gam_ratio,xeq_a,na_Ma,nc_Mc,xeq_c,mw_electrolyte, &
            Keq_sl,MW_c,MW_a,Keq_ll,eleliquid,flux_sl,phi_salt,sat_ratio,        &
            molality0,jsalt_present,kappa_nonelectro)

       
       
       call MESA_convergence_criterion(ibin,iconverge_mass,iconverge_flux,idissolved, &
            aer,electrolyte,mass_dry_a,mass_dry_salt,mw_electrolyte,mw_aer_mac, &
            flux_sl,phi_salt,rtol_mesa)
       
       if(iconverge_mass .eq. mYES)then
          mosaic_vars_aa%iter_MESA(ibin) = mosaic_vars_aa%iter_MESA(ibin) + itdum
          mosaic_vars_aa%niter_MESA = mosaic_vars_aa%niter_MESA + float(itdum)
          mosaic_vars_aa%niter_MESA_max = max( mosaic_vars_aa%niter_MESA_max, itdum)
          jaerosolstate(ibin) = all_solid
          call adjust_solid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent,   &
               water_a)
          jhyst_leg(ibin) = jhyst_lo
          growth_factor(ibin) = 1.0
          return
       elseif(iconverge_flux .eq. mYES)then
          mosaic_vars_aa%iter_MESA(ibin) = mosaic_vars_aa%iter_MESA(ibin) + itdum
          mosaic_vars_aa%niter_MESA = mosaic_vars_aa%niter_MESA + float(itdum)
          mosaic_vars_aa%niter_MESA_max = max( mosaic_vars_aa%niter_MESA_max, itdum)
          jaerosolstate(ibin) = mixed
          vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3          
          growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)           
          
          if(idissolved .eq. myes)then
             jaerosolstate(ibin) = all_liquid
             
          else
             jaerosolstate(ibin) = mixed
             jhyst_leg(ibin) = jhyst_lo
          endif
             
          
          
          
          
          
          
          
          
          
          
          
          
          
          return
       endif
       
       
       hsalt_min = 1.e25
      
       do js = 1, nsalt
          
          phi_prod = phi_salt(js) * phi_salt_old(js)

          if(itdum .gt. 1 .and. phi_prod .gt. 0.0)then
             phi_bar(js) = (abs(phi_salt(js))-abs(phi_salt_old(js)))/   &
                  alpha_salt(js)
          else
             phi_bar(js) = 0.0                      
          endif

          if(phi_bar(js) .lt. 0.0)then              
             phi_bar(js) = max(phi_bar(js), -10.0d0)
             alpha_fac = 3.0*exp(phi_bar(js))
             alpha_salt(js) = min(alpha_fac*abs(phi_salt(js)), 0.9d0)
          elseif(phi_bar(js) .gt. 0.0)then  
             alpha_salt(js) = min(abs(phi_salt(js)), 0.5d0)
          else                                      
             alpha_salt(js) = min(abs(phi_salt(js))/3.0d0, 0.5d0)
          endif
          
          
          
          phi_salt_old(js) = phi_salt(js)           
          

          if(flux_sl(js) .gt. 0.)then
             
             tau_p(js) = eleliquid(js)/flux_sl(js)  
             if(tau_p(js) .eq. 0.0)then
                hsalt(js) = 1.e25
                flux_sl(js) = 0.0
                phi_salt(js)= 0.0
             else
                hsalt(js) = alpha_salt(js)*tau_p(js)
             endif
             
          elseif(flux_sl(js) .lt. 0.)then
             
             tau_p(js) = -eleliquid(js)/flux_sl(js) 
             tau_d(js) = -electrolyte(js,jsolid,ibin)/flux_sl(js) 
             if(tau_p(js) .eq. 0.0)then
                hsalt(js) = alpha_salt(js)*tau_d(js)
             else
                hsalt(js) = alpha_salt(js)*min(tau_p(js),tau_d(js))
             endif
             
          else
             
             hsalt(js) = 1.e25
             
          endif
          
          hsalt_min = min(hsalt(js), hsalt_min)
          
       enddo

       
       
       
       do js = 1, nsalt
          electrolyte(js,jsolid,ibin) = (   &
               (electrolyte(js,jsolid,ibin))  +   &
               (hsalt(js)) * (flux_sl(js)) )
       enddo
       
       
       
       call electrolytes_to_ions(jsolid,ibin,aer,electrolyte)
       
       
       
       do iaer = 1, naer
          aer(iaer,jliquid,ibin) = ( (aer(iaer,jtotal,ibin)) -   &
               (aer(iaer,jsolid,ibin)) )
       enddo
       
       
       

       
500 continue     
    
    mosaic_vars_aa%jMESA_fail = mosaic_vars_aa%jMESA_fail + 1
    mosaic_vars_aa%iter_MESA(ibin) = mosaic_vars_aa%iter_MESA(ibin) + itdum
    mosaic_vars_aa%niter_MESA = mosaic_vars_aa%niter_MESA + float(itdum)
    jaerosolstate(ibin) = mixed
    jhyst_leg(ibin) = jhyst_lo
    mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3    
    vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3     
    growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)      
   
    return
  end subroutine MESA_PTC



  
  
  
  
  
  
  subroutine MESA_flux_salt(ibin, jaerosolstate,jphase,aer,jhyst_leg,electrolyte,  &
       epercent,activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,ma,gam,   &
       log_gamZ,zc,za,gam_ratio,xeq_a,na_Ma,nc_Mc,xeq_c,mw_electrolyte,Keq_sl,MW_c,&
       MW_a,Keq_ll,eleliquid,flux_sl,phi_salt,sat_ratio,molality0,jsalt_present,   &
       kappa_nonelectro                                                            )      

    use module_data_mosaic_aero, only: r8,nbin_a_max,nelectrolyte,Ncation,naer,    &
         jliquid,nsalt,jsolid,Nanion,nrxn_aer_sl,nrxn_aer_ll,nrxn_aer_sl,          &
         jna3hso4,ica_a,jcano3,jcacl2                                               

    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(inout), dimension(nsalt) :: jsalt_present
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase,jhyst_leg

    real(r8), intent(in) :: aH2O
    real(r8), intent(inout), dimension(nsalt) :: flux_sl,phi_salt,sat_ratio
    real(r8), intent(in), dimension(Ncation) :: zc,MW_c
    real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
    real(r8), intent(in), dimension(Nanion)  :: za,MW_a
    real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a,mass_dry_a,gam_ratio
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_soluble_a,water_a
    real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
    real(r8), intent(inout), dimension(nelectrolyte) :: eleliquid
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nrxn_aer_sl) :: Keq_sl
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte,epercent
    real(r8), intent(in), dimension(naer) :: kappa_nonelectro

    
    integer js, je
    real(r8) :: XT, calcium, sum_salt, sum_dum 
    real(r8), dimension(nsalt) :: frac_salt_liq,frac_salt_solid


    
    call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,   &
         nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)
    call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,electrolyte,   &
         activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,ma,gam,log_gamZ, &
         gam_ratio,Keq_ll,molality0,kappa_nonelectro)
    activity(jna3hso4,ibin)   = 0.0

    if(water_a(ibin) .le. 0.0)then
       do js = 1, nsalt
          flux_sl(js) = 0.0
       enddo
       return
    endif


    call MESA_estimate_eleliquid(ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,nc_Mc,  &
         xeq_c,mw_electrolyte,MW_c,MW_a,eleliquid)

    calcium = aer(ica_a,jliquid,ibin)



    
    sum_dum = 0.0
    do je = 1, nelectrolyte
       sum_dum = sum_dum + electrolyte(je,jliquid,ibin)
    enddo

    if(sum_dum .eq. 0.)sum_dum = 1.0

    do je = 1, nelectrolyte
       epercent(je,jliquid,ibin) = 100.*electrolyte(je,jliquid,ibin)/sum_dum
    enddo
    



    
    sum_salt = 0.0
    do js = 1, nsalt
       sum_salt = sum_salt + electrolyte(js,jsolid,ibin)
    enddo

    if(sum_salt .eq. 0.0)sum_salt = 1.0
    do js = 1, nsalt
       frac_salt_solid(js) = electrolyte(js,jsolid,ibin)/sum_salt
       frac_salt_liq(js)   = epercent(js,jliquid,ibin)/100.
    enddo

    
    do js = 1, nsalt             

       
       sat_ratio(js) = activity(js,ibin)/Keq_sl(js)
       
       phi_salt(js)  = (sat_ratio(js) - 1.0)/max(sat_ratio(js),1.0d0)

       
       if(sat_ratio(js)       .lt. 1.00 .and.   &
            frac_salt_solid(js) .lt. 0.01 .and.   &
            frac_salt_solid(js) .gt. 0.0)then
          call MESA_dissolve_small_salt(ibin,js,aer,electrolyte)
          call MESA_estimate_eleliquid(ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,  &
               nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a,eleliquid)
          sat_ratio(js) = activity(js,ibin)/Keq_sl(js)
       endif

       
       flux_sl(js) = sat_ratio(js) - 1.0

       
       if( (sat_ratio(js)               .lt. 1.0 .and.   &
            electrolyte(js,jsolid,ibin) .eq. 0.0) .or.   &
            (calcium .gt. 0.0 .and. frac_salt_liq(js).lt.0.01).or.   &
            (calcium .gt. 0.0 .and. jsalt_present(js).eq.0) )then
          flux_sl(js) = 0.0
          phi_salt(js)= 0.0
       endif

    enddo


    
    sat_ratio(jcano3) = 1.0
    phi_salt(jcano3)  = 0.0
    flux_sl(jcano3)   = 0.0

    sat_ratio(jcacl2) = 1.0
    phi_salt(jcacl2)  = 0.0
    flux_sl(jcacl2)   = 0.0


    return
  end subroutine MESA_flux_salt

 
  
  
  
  
  
  subroutine compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,           &
       electrolyte,activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,ma,gam,&
       log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)

    use module_data_mosaic_aero, only: r8,nbin_a_max,nelectrolyte,Ncation,naer,    &
         jliquid,Nanion,nrxn_aer_ll,                                               &
         iso4_a,ja_so4,ja_hso4,ino3_a,ja_no3,icl_a,ja_cl,imsa_a,ja_msa,ica_a,jc_ca,&
         inh4_a,jc_nh4,ina_a,jc_na,jc_h,jhcl,jhno3,jcacl2,jcano3,jnacl,jnano3,     &
         jna2so4,jnh4so4,jnh4cl,jnh4no3,jlvcite,jnh4hso4,jnh4msa,jna3hso4,jnahso4, &
         jnamsa,jcamsa2,jh2so4,jhhso4,jmsa                                          

    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase,jhyst_leg

    real(r8), intent(in) :: aH2O
    real(r8), intent(in), dimension(nbin_a_max) :: num_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a,mass_soluble_a
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a,gam_ratio
    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(in), dimension(naer) :: kappa_nonelectro

    
    real(r8), dimension(nelectrolyte) :: log_gam
    integer jp, jA
    real(r8) :: XT, xmol(Nelectrolyte), sum_elec, dumK, c_bal, a_c 
    real(r8) :: quad, aq, bq, cq, xq, dum, mSULF
    


    water_a(ibin) = aerosol_water(jliquid,ibin,jaerosolstate,jphase, &
         jhyst_leg,electrolyte,aer,kappa_nonelectro,num_a,mass_dry_a,mass_soluble_a,aH2O, &
         molality0)      
    if(water_a(ibin) .eq. 0.0)return


    call calculate_XT(ibin,jliquid,XT,aer)


    if(XT.ge.2.0 .or. XT.lt.0.)then   
       


       
       ma(ja_so4,ibin)  = 1.e-9*aer(iso4_a,jliquid,ibin)/water_a(ibin)
       ma(ja_hso4,ibin) = 0.0
       ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
       ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)
       ma(ja_msa,ibin)  = 1.e-9*aer(imsa_a,jliquid,ibin)/water_a(ibin)

       
       mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
       mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
       mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)
       a_c              = (   &
            (2.*ma(ja_so4,ibin)+   &
            ma(ja_no3,ibin)+   &
            ma(ja_cl,ibin) +   &
            ma(ja_msa,ibin)) -   &
            (2.*mc(jc_ca,ibin) +   &
            mc(jc_nh4,ibin)+   &
            mc(jc_na,ibin)) )

       mc(jc_h,ibin) = 0.5*( (a_c) +   &
            (sqrt(a_c**2 + 4.*Keq_ll(3))) )

       if(mc(jc_h,ibin) .le. 0.0)then   
          mc(jc_h,ibin) = 1.e-10
       endif


       jp = jliquid


       sum_elec = 2.*electrolyte(jnh4no3,jp,ibin) +   &
            2.*electrolyte(jnh4cl,jp,ibin)  +   &
            3.*electrolyte(jnh4so4,jp,ibin) +   &
            3.*electrolyte(jna2so4,jp,ibin) +   &
            2.*electrolyte(jnano3,jp,ibin)  +   &
            2.*electrolyte(jnacl,jp,ibin)   +   &
            3.*electrolyte(jcano3,jp,ibin)  +   &
            3.*electrolyte(jcacl2,jp,ibin)  +   &
            2.*electrolyte(jhno3,jp,ibin)   +   &
            2.*electrolyte(jhcl,jp,ibin)

       if(sum_elec .eq. 0.0)then
          do jA = 1, nelectrolyte
             gam(jA,ibin) = 1.0
          enddo
          goto 10
       endif


       
       xmol(jnh4no3) = 2.*electrolyte(jnh4no3,jp,ibin)/sum_elec
       xmol(jnh4cl)  = 2.*electrolyte(jnh4cl,jp,ibin) /sum_elec
       xmol(jnh4so4) = 3.*electrolyte(jnh4so4,jp,ibin)/sum_elec
       xmol(jna2so4) = 3.*electrolyte(jna2so4,jp,ibin)/sum_elec
       xmol(jnano3)  = 2.*electrolyte(jnano3,jp,ibin) /sum_elec
       xmol(jnacl)   = 2.*electrolyte(jnacl,jp,ibin)  /sum_elec
       xmol(jcano3)  = 3.*electrolyte(jcano3,jp,ibin) /sum_elec
       xmol(jcacl2)  = 3.*electrolyte(jcacl2,jp,ibin) /sum_elec
       xmol(jhno3)   = 2.*electrolyte(jhno3,jp,ibin)  /sum_elec
       xmol(jhcl)    = 2.*electrolyte(jhcl,jp,ibin)   /sum_elec


       jA = jnh4so4
       if(xmol(jA).gt.0.0)then
          log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
               xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
               xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
               xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
               xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
               xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
               xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
               xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
               xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
               xmol(jhcl)   *log_gamZ(jA,jhcl)
          gam(jA,ibin) = 10.**log_gam(jA)
          activity(jnh4so4,ibin) = mc(jc_nh4,ibin)**2 * ma(ja_so4,ibin) *   &
               gam(jnh4so4,ibin)**3
       endif





       jA = jnh4no3

          log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
               xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
               xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
               xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
               xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
               xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
               xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
               xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
               xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
               xmol(jhcl)   *log_gamZ(jA,jhcl)
          gam(jA,ibin) = 10.**log_gam(jA)
          activity(jnh4no3,ibin) = mc(jc_nh4,ibin) * ma(ja_no3,ibin) *   &
               gam(jnh4no3,ibin)**2



       jA = jnh4cl
       if(xmol(jA).gt.0.0)then
          log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
               xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
               xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
               xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
               xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
               xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
               xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
               xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
               xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
               xmol(jhcl)   *log_gamZ(jA,jhcl)
          gam(jA,ibin) = 10.**log_gam(jA)
          activity(jnh4cl,ibin)  = mc(jc_nh4,ibin) * ma(ja_cl,ibin) *   &
               gam(jnh4cl,ibin)**2
       endif


       jA = jna2so4
       if(xmol(jA).gt.0.0)then
          log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
               xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
               xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
               xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
               xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
               xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
               xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
               xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
               xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
               xmol(jhcl)   *log_gamZ(jA,jhcl)
          gam(jA,ibin) = 10.**log_gam(jA)
          activity(jna2so4,ibin) = mc(jc_na,ibin)**2 * ma(ja_so4,ibin) *   &
               gam(jna2so4,ibin)**3
       endif


       jA = jnano3
       if(xmol(jA).gt.0.0)then
          log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
               xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
               xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
               xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
               xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
               xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
               xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
               xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
               xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
               xmol(jhcl)   *log_gamZ(jA,jhcl)
          gam(jA,ibin) = 10.**log_gam(jA)
          activity(jnano3,ibin)  = mc(jc_na,ibin) * ma(ja_no3,ibin) *   &
               gam(jnano3,ibin)**2
       endif



       jA = jnacl
       if(xmol(jA).gt.0.0)then
          log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
               xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
               xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
               xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
               xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
               xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
               xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
               xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
               xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
               xmol(jhcl)   *log_gamZ(jA,jhcl)
          gam(jA,ibin) = 10.**log_gam(jA)
          activity(jnacl,ibin)   = mc(jc_na,ibin) * ma(ja_cl,ibin) *   &
               gam(jnacl,ibin)**2
       endif



       
       
       
       
       



       
       
       
       
       

       jA = jcano3
       if(xmol(jA).gt.0.0)then
          log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
               xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
               xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
               xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
               xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
               xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
               xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
               xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
               xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
               xmol(jhcl)   *log_gamZ(jA,jhcl)
          gam(jA,ibin) = 10.**log_gam(jA)
          activity(jcano3,ibin)  = mc(jc_ca,ibin) * ma(ja_no3,ibin)**2 *   &
               gam(jcano3,ibin)**3
       endif



       jA = jcacl2
       if(xmol(jA).gt.0.0)then
          log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
               xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
               xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
               xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
               xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
               xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
               xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
               xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
               xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
               xmol(jhcl)   *log_gamZ(jA,jhcl)
          gam(jA,ibin) = 10.**log_gam(jA)
          activity(jcacl2,ibin)  = mc(jc_ca,ibin) * ma(ja_cl,ibin)**2 *   &
               gam(jcacl2,ibin)**3
       endif


       jA = jhno3
       log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
            xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
            xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
            xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
            xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
            xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
            xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
            xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
            xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)   *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)
       activity(jhno3,ibin)   = mc(jc_h,ibin) * ma(ja_no3,ibin) *   &
            gam(jhno3,ibin)**2


       jA = jhcl
       log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
            xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
            xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
            xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
            xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
            xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
            xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
            xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
            xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)   *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)
       activity(jhcl,ibin)    = mc(jc_h,ibin) * ma(ja_cl,ibin) *   &
            gam(jhcl,ibin)**2

       
10     gam(jlvcite,ibin) = 1.0

       gam(jnh4hso4,ibin)= 1.0

       gam(jnh4msa,ibin) = 1.0

       gam(jna3hso4,ibin) = 1.0

       gam(jnahso4,ibin) = 1.0

       gam(jnamsa,ibin)  = 1.0

       gam(jcamsa2,ibin) = 1.0

       activity(jlvcite,ibin) = 0.0

       activity(jnh4hso4,ibin)= 0.0

       activity(jnh4msa,ibin) = mc(jc_nh4,ibin) * ma(ja_msa,ibin) *   &
            gam(jnh4msa,ibin)**2

       activity(jna3hso4,ibin)= 0.0

       activity(jnahso4,ibin) = 0.0

       activity(jnamsa,ibin) = mc(jc_na,ibin) * ma(ja_msa,ibin) *   &
            gam(jnamsa,ibin)**2

       activity(jcamsa2,ibin) = mc(jc_ca,ibin) * ma(ja_msa,ibin)**2 *   &
            gam(jcamsa2,ibin)**3

       gam_ratio(ibin) = gam(jnh4no3,ibin)**2/gam(jhno3,ibin)**2


    else
       

       jp = jliquid

       sum_elec = 3.*electrolyte(jh2so4,jp,ibin)    +   &
            2.*electrolyte(jnh4hso4,jp,ibin)  +   &
            5.*electrolyte(jlvcite,jp,ibin)   +   &
            3.*electrolyte(jnh4so4,jp,ibin)   +   &
            2.*electrolyte(jnahso4,jp,ibin)   +   &
            5.*electrolyte(jna3hso4,jp,ibin)  +   &
            3.*electrolyte(jna2so4,jp,ibin)   +   &
            2.*electrolyte(jhno3,jp,ibin)     +   &
            2.*electrolyte(jhcl,jp,ibin)


       if(sum_elec .eq. 0.0)then
          do jA = 1, nelectrolyte
             gam(jA,ibin) = 1.0
          enddo
          goto 20
       endif


       xmol(jh2so4)  = 3.*electrolyte(jh2so4,jp,ibin)/sum_elec
       xmol(jnh4hso4)= 2.*electrolyte(jnh4hso4,jp,ibin)/sum_elec
       xmol(jlvcite) = 5.*electrolyte(jlvcite,jp,ibin)/sum_elec
       xmol(jnh4so4) = 3.*electrolyte(jnh4so4,jp,ibin)/sum_elec
       xmol(jnahso4) = 2.*electrolyte(jnahso4,jp,ibin)/sum_elec
       xmol(jna3hso4)= 5.*electrolyte(jna3hso4,jp,ibin)/sum_elec
       xmol(jna2so4) = 3.*electrolyte(jna2so4,jp,ibin)/sum_elec
       xmol(jhno3)   = 2.*electrolyte(jhno3,jp,ibin)/sum_elec
       xmol(jhcl)    = 2.*electrolyte(jhcl,jp,ibin)/sum_elec


       
       jA = jh2so4
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


       
       jA = jhhso4
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


       
       jA = jnh4hso4
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


       
       jA = jlvcite
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


       
       jA = jnh4so4
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


       
       jA = jnahso4
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


       
       jA = jna3hso4
       
       
       
       
       
       
       
       
       
       
       gam(jA,ibin) = 1.0


       
       jA = jna2so4
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


       
       jA = jhno3
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


       
       jA = jhcl
       log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
            xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
            xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
            xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
            xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
            xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
            xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
            xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
            xmol(jhcl)    *log_gamZ(jA,jhcl)
       gam(jA,ibin) = 10.**log_gam(jA)


20     gam(jnh4no3,ibin) = 1.0
       gam(jnh4cl,ibin)  = 1.0
       gam(jnano3,ibin)  = 1.0
       gam(jnacl,ibin)   = 1.0
       gam(jcano3,ibin)  = 1.0
       gam(jcacl2,ibin)  = 1.0

       gam(jnh4msa,ibin) = 1.0
       gam(jnamsa,ibin)  = 1.0
       gam(jcamsa2,ibin) = 1.0



       
       
       mc(jc_ca,ibin)   = 1.e-9*aer(ica_a,jliquid,ibin)/water_a(ibin)
       mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
       mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)

       
       mSULF            = 1.e-9*aer(iso4_a,jliquid,ibin)/water_a(ibin)
       ma(ja_hso4,ibin) = 0.0
       ma(ja_so4,ibin)  = 0.0
       ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
       ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)
       ma(ja_msa,ibin)  = 1.e-9*aer(imsa_a,jliquid,ibin)/water_a(ibin)

       gam_ratio(ibin)  = gam(jnh4hso4,ibin)**2/gam(jhhso4,ibin)**2
       dumK = Keq_ll(1)*gam(jhhso4,ibin)**2/gam(jh2so4,ibin)**3

       c_bal =  mc(jc_nh4,ibin) + mc(jc_na,ibin) + 2.*mc(jc_ca,ibin)   &
            - ma(ja_no3,ibin) - ma(ja_cl,ibin) - mSULF - ma(ja_msa,ibin)

       aq = 1.0
       bq = dumK + c_bal
       cq = dumK*(c_bal - mSULF)


       
       if(bq .ne. 0.0)then
          xq = 4.*(1./bq)*(cq/bq)
       else
          xq = 1.e+6
       endif

       if(abs(xq) .lt. 1.e-6)then
          dum = xq*(0.5 + xq*(0.125 + xq*0.0625))
          quad = (-0.5*bq/aq)*dum
          if(quad .lt. 0.)then
             quad = -bq/aq - quad
          endif
       else
          quad = 0.5*(-bq+sqrt(bq*bq - 4.*cq))
       endif
       

       mc(jc_h,ibin) = max(quad, 1.d-7)
       ma(ja_so4,ibin) = mSULF*dumK/(mc(jc_h,ibin) + dumK)
       ma(ja_hso4,ibin)= mSULF - ma(ja_so4,ibin)

       activity(jcamsa2,ibin) = mc(jc_ca,ibin) * ma(ja_msa,ibin)**2 *   &
            gam(jcamsa2,ibin)**3

       activity(jnh4so4,ibin) = mc(jc_nh4,ibin)**2 * ma(ja_so4,ibin) *   &
            gam(jnh4so4,ibin)**3

       activity(jlvcite,ibin) = mc(jc_nh4,ibin)**3 * ma(ja_hso4,ibin) *   &
            ma(ja_so4,ibin) * gam(jlvcite,ibin)**5

       activity(jnh4hso4,ibin)= mc(jc_nh4,ibin) * ma(ja_hso4,ibin) *   &
            gam(jnh4hso4,ibin)**2

       activity(jnh4msa,ibin) = mc(jc_nh4,ibin) * ma(ja_msa,ibin) *   &
            gam(jnh4msa,ibin)**2

       activity(jna2so4,ibin) = mc(jc_na,ibin)**2 * ma(ja_so4,ibin) *   &
            gam(jna2so4,ibin)**3

       activity(jnahso4,ibin) = mc(jc_na,ibin) * ma(ja_hso4,ibin) *   &
            gam(jnahso4,ibin)**2

       activity(jnamsa,ibin)  = mc(jc_na,ibin) * ma(ja_msa,ibin) *   &
            gam(jnamsa,ibin)**2

       
       

       activity(jna3hso4,ibin)= 0.0

       activity(jhno3,ibin)   = mc(jc_h,ibin) * ma(ja_no3,ibin) *   &
            gam(jhno3,ibin)**2

       activity(jhcl,ibin)    = mc(jc_h,ibin) * ma(ja_cl,ibin) *   &
            gam(jhcl,ibin)**2

       activity(jmsa,ibin)    = mc(jc_h,ibin) * ma(ja_msa,ibin) *   &
            gam(jmsa,ibin)**2


       
       activity(jnh4no3,ibin) = 0.0

       activity(jnh4cl,ibin)  = 0.0

       activity(jnano3,ibin)  = 0.0

       activity(jnacl,ibin)   = 0.0

       activity(jcano3,ibin)  = 0.0

       activity(jcacl2,ibin)  = 0.0


    endif
    return
  end subroutine compute_activities



  
  
  
  
  
  
  subroutine MESA_convergence_criterion(ibin,iconverge_mass,iconverge_flux,        &
       idissolved,aer,electrolyte,mass_dry_a,mass_dry_salt,                        &
       mw_electrolyte,mw_aer_mac,flux_sl,phi_salt,rtol_mesa)  

    use module_data_mosaic_aero, only: r8,nbin_a_max,naer,nelectrolyte,nsalt,      &
         jsolid, mYES, xhyst_up_crustal_thresh,                                    &
         mno,ioin_a,jcaso4,jcaco3         

    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(inout) :: iconverge_mass, iconverge_flux, idissolved
    real(r8), intent(in) :: rtol_mesa
    real(r8), intent(inout), dimension(nsalt) :: flux_sl,phi_salt
    real(r8), intent(in),    dimension(nbin_a_max) :: mass_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_salt
    real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
    real(r8), intent(in), dimension(naer) :: mw_aer_mac
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer je, js, iaer
    real(r8) :: mass_solid, mass_solid_salt,frac_solid, XT, H_ion,   &
         crustal_solids, sumflux


    idissolved = mno             

    
    iconverge_mass = mNO 

    
    
    
    
    
    

    mass_solid_salt = 0.0
    do je = 1, nsalt
       mass_solid_salt = mass_solid_salt +   &
            electrolyte(je,jsolid,ibin)*mw_electrolyte(je)*1.e-15        
    enddo



    


    if(mass_dry_salt(ibin) .le. 0.0)then
      frac_solid = 0.0
    else
      frac_solid = mass_solid_salt/mass_dry_salt(ibin)
    endif


    if(frac_solid .ge. 0.98)then
       iconverge_mass = mYES
       return
    endif



    
    iconverge_flux = mYES
    do js = 1, nsalt
       if(abs(phi_salt(js)).gt. rtol_mesa)then
          iconverge_flux = mNO
          return
       endif
    enddo



    

    sumflux = 0.0
    do js = 1, nsalt
       sumflux = sumflux + abs(flux_sl(js))
    enddo




    crustal_solids = electrolyte(jcaco3,jsolid,ibin)*mw_electrolyte(jcaco3) +  &
                     electrolyte(jcaso4,jsolid,ibin)*mw_electrolyte(jcaso4) +  &
                     aer(ioin_a,jsolid,ibin)*mw_aer_mac(ioin_a)


    if ( sumflux .eq. 0.0 .and. &
         crustal_solids .le. xhyst_up_crustal_thresh*(mass_dry_a(ibin)*1.0e15) ) then
       
       idissolved = myes
    endif



    return
  end subroutine MESA_convergence_criterion



  
  
  
  
  
  
  subroutine electrolytes_to_ions(jp,ibin,aer,electrolyte)

    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         jh2so4,jnh4hso4,jlvcite,jnh4so4,jnahso4,jna3hso4,jna2so4,jcaso4,iso4_a,   &
         jhno3,jnh4no3,jcano3,jnano3,ino3_a,jhcl,jnh4cl,jcacl2,jnacl,icl_a,jmsa,   &
         jcamsa2,jnamsa,jnh4msa,imsa_a,jcaco3,ico3_a,ica_a,ina_a,inh4_a             
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    real(r8) :: sum_dum


    aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
         electrolyte(jna2so4,jp,ibin) +   &
         2.*electrolyte(jna3hso4,jp,ibin)+   &
         electrolyte(jnahso4,jp,ibin) +   &
         electrolyte(jnh4so4,jp,ibin) +   &
         2.*electrolyte(jlvcite,jp,ibin) +   &
         electrolyte(jnh4hso4,jp,ibin)+   &
         electrolyte(jh2so4,jp,ibin)

    aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
         2.*electrolyte(jcano3,jp,ibin)  +   &
         electrolyte(jnh4no3,jp,ibin) +   &
         electrolyte(jhno3,jp,ibin)

    aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
         2.*electrolyte(jcacl2,jp,ibin)  +   &
         electrolyte(jnh4cl,jp,ibin)  +   &
         electrolyte(jhcl,jp,ibin)

    aer(imsa_a,jp,ibin) = electrolyte(jnh4msa,jp,ibin) +   &
         electrolyte(jnamsa,jp,ibin)  +   &
         2.*electrolyte(jcamsa2,jp,ibin) +   &
         electrolyte(jmsa,jp,ibin)

    aer(ico3_a,jp,ibin) = electrolyte(jcaco3,jp,ibin)

    aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
         electrolyte(jcano3,jp,ibin)  +   &
         electrolyte(jcacl2,jp,ibin)  +   &
         electrolyte(jcaco3,jp,ibin)  +   &
         electrolyte(jcamsa2,jp,ibin)

    aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
         electrolyte(jnacl,jp,ibin)   +   &
         2.*electrolyte(jna2so4,jp,ibin) +   &
         3.*electrolyte(jna3hso4,jp,ibin)+   &
         electrolyte(jnahso4,jp,ibin) +   &
         electrolyte(jnamsa,jp,ibin)

    aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
         electrolyte(jnh4cl,jp,ibin)  +   &
         2.*electrolyte(jnh4so4,jp,ibin) +   &
         3.*electrolyte(jlvcite,jp,ibin) +   &
         electrolyte(jnh4hso4,jp,ibin)+   &
         electrolyte(jnh4msa,jp,ibin)


    return
  end subroutine electrolytes_to_ions



  
  
  
  
  
  
  
  
  
  
  subroutine ions_to_electrolytes(jp,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,    &
       nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)

    use module_data_mosaic_aero, only: r8,naer,nbin_a_max,nelectrolyte,ncation,    &
         nanion,jliquid,jsolid,                                                    &
         ica_a,iso4_a,jcaso4,imsa_a,ina_a,inh4_a,ja_hso4,ja_so4,ino3_a,ja_no3,     &
         icl_a,ja_cl,ja_msa,jc_ca,jc_na,jc_nh4,jc_h,jna2so4,jnahso4,jnamsa,jnano3, &
         jnacl,jnh4so4,jnh4hso4,jnh4msa,jnh4no3,jnh4cl,jcano3,jcacl2,jcamsa2,      &
         jh2so4,jhno3,jhcl,jmsa,jlvcite,jna3hso4                                    
    implicit none

    
    integer, intent(in) :: ibin, jp
    real(r8), intent(inout) :: XT
    real(r8), intent(in), dimension(Ncation) :: zc,MW_c
    real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
    real(r8), intent(in), dimension(Nanion) :: za,MW_a
    real(r8), intent(inout), dimension(Nanion) :: xeq_a,na_Ma
    real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    character(len=256) :: errmsg
    integer iaer, je, jc, ja, icase
    real(r8) :: store(naer), sum_dum, sum_naza, sum_nczc, sum_na_nh4,   &
         f_nh4, f_na, xh, xb, xl, xs, cat_net, rem_nh4, rem_na
    real(r8) :: nc(ncation), na(nanion)




    if(jp .ne. jliquid)then
       call wrf_message(' jp must be jliquid')
       call wrf_message(' in ions_to_electrolytes sub')
       write(errmsg,*)' wrong jp = ', jp
       call wrf_error_fatal3("<stdin>",2499,&
trim(adjustl(errmsg)))
       
    endif

    
    
    
    


    
    store(ica_a)  = aer(ica_a, jp,ibin)
    store(iso4_a) = aer(iso4_a,jp,ibin)

    call form_caso4(store,jp,ibin,electrolyte)

    if(jp .eq. jliquid)then 
       aer(ica_a,jliquid,ibin) = aer(ica_a,jliquid,ibin) -   &
            electrolyte(jcaso4,jliquid,ibin)

       aer(iso4_a,jliquid,ibin)= aer(iso4_a,jliquid,ibin)-   &
            electrolyte(jcaso4,jliquid,ibin)

       aer(ica_a,jsolid,ibin)  = aer(ica_a,jsolid,ibin) +   &
            electrolyte(jcaso4,jliquid,ibin)

       aer(iso4_a,jsolid,ibin) = aer(iso4_a,jsolid,ibin) +   &
            electrolyte(jcaso4,jliquid,ibin)

       electrolyte(jcaso4,jsolid,ibin)=electrolyte(jcaso4,jsolid,ibin)   &
            +electrolyte(jcaso4,jliquid,ibin)
       electrolyte(jcaso4,jliquid,ibin)= 0.0
    endif


    
    

    if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
       XT   = ( aer(inh4_a,jp,ibin) +   &
            aer(ina_a,jp,ibin)  +   &
            2.*aer(ica_a,jp,ibin) )/   &
            (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
    else
       XT   = -1.0
    endif





    if(XT.ge.2.0 .or. XT.lt.0.)then     
       icase = 1  
    else
       icase = 2  
    endif


    
    do je = 1, nelectrolyte
       electrolyte(je,jp,ibin) = 0.0
    enddo

    
    
    

    if(icase.eq.1)then 

       na(ja_hso4)= 0.0
       na(ja_so4) = aer(iso4_a,jp,ibin)
       na(ja_no3) = aer(ino3_a,jp,ibin)
       na(ja_cl)  = aer(icl_a, jp,ibin)
       na(ja_msa) = aer(imsa_a,jp,ibin)

       nc(jc_ca)  = aer(ica_a, jp,ibin)
       nc(jc_na)  = aer(ina_a, jp,ibin)
       nc(jc_nh4) = aer(inh4_a,jp,ibin)

       cat_net = (   &
            (2.*na(ja_so4)+na(ja_no3)+na(ja_cl)+na(ja_msa)) -   &
            (2.*nc(jc_ca) +nc(jc_nh4)+nc(jc_na)) )

       if(cat_net .lt. 0.0)then

          nc(jc_h) = 0.0

       else  

          nc(jc_h) = cat_net

       endif


       
       sum_naza = 0.0
       do ja = 1, nanion
          sum_naza = sum_naza + na(ja)*za(ja)
       enddo

       sum_nczc = 0.0
       do jc = 1, ncation
          sum_nczc = sum_nczc + nc(jc)*zc(jc)
       enddo

       if(sum_naza .eq. 0. .or. sum_nczc .eq. 0.)then 



          return
       endif

       do ja = 1, nanion
          xeq_a(ja) = na(ja)*za(ja)/sum_naza
       enddo

       do jc = 1, ncation
          xeq_c(jc) = nc(jc)*zc(jc)/sum_nczc
       enddo

       na_Ma(ja_so4) = na(ja_so4) *MW_a(ja_so4)
       na_Ma(ja_no3) = na(ja_no3) *MW_a(ja_no3)
       na_Ma(ja_cl)  = na(ja_cl)  *MW_a(ja_cl)
       na_Ma(ja_msa) = na(ja_msa) *MW_a(ja_msa)
       na_Ma(ja_hso4)= na(ja_hso4)*MW_a(ja_hso4)

       nc_Mc(jc_ca)  = nc(jc_ca) *MW_c(jc_ca)
       nc_Mc(jc_na)  = nc(jc_na) *MW_c(jc_na)
       nc_Mc(jc_nh4) = nc(jc_nh4)*MW_c(jc_nh4)
       nc_Mc(jc_h)   = nc(jc_h)  *MW_c(jc_h)


       
       if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_so4) .gt. 0.)then
          electrolyte(jna2so4,jp,ibin) = (xeq_c(jc_na) *na_Ma(ja_so4) +   &
               xeq_a(ja_so4)*nc_Mc(jc_na))/   &
               mw_electrolyte(jna2so4)
       endif

       electrolyte(jnahso4,jp,ibin) = 0.0

       if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
          electrolyte(jnamsa,jp,ibin)  = (xeq_c(jc_na) *na_Ma(ja_msa) +   &
               xeq_a(ja_msa)*nc_Mc(jc_na))/   &
               mw_electrolyte(jnamsa)
       endif

       if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
          electrolyte(jnano3,jp,ibin)  = (xeq_c(jc_na) *na_Ma(ja_no3) +   &
               xeq_a(ja_no3)*nc_Mc(jc_na))/   &
               mw_electrolyte(jnano3)
       endif

       if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
          electrolyte(jnacl,jp,ibin)   = (xeq_c(jc_na) *na_Ma(ja_cl) +   &
               xeq_a(ja_cl) *nc_Mc(jc_na))/   &
               mw_electrolyte(jnacl)
       endif

       if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_so4) .gt. 0.)then
          electrolyte(jnh4so4,jp,ibin) = (xeq_c(jc_nh4)*na_Ma(ja_so4) +   &
               xeq_a(ja_so4)*nc_Mc(jc_nh4))/   &
               mw_electrolyte(jnh4so4)
       endif

       electrolyte(jnh4hso4,jp,ibin)= 0.0

       if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
          electrolyte(jnh4msa,jp,ibin) = (xeq_c(jc_nh4)*na_Ma(ja_msa) +   &
               xeq_a(ja_msa)*nc_Mc(jc_nh4))/   &
               mw_electrolyte(jnh4msa)
       endif

       if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
          electrolyte(jnh4no3,jp,ibin) = (xeq_c(jc_nh4)*na_Ma(ja_no3) +   &
               xeq_a(ja_no3)*nc_Mc(jc_nh4))/   &
               mw_electrolyte(jnh4no3)
       endif

       if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
          electrolyte(jnh4cl,jp,ibin)  = (xeq_c(jc_nh4)*na_Ma(ja_cl) +   &
               xeq_a(ja_cl) *nc_Mc(jc_nh4))/   &
               mw_electrolyte(jnh4cl)
       endif

       if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.0)then
          electrolyte(jcano3, jp,ibin) = (xeq_c(jc_ca) *na_Ma(ja_no3) +   &
               xeq_a(ja_no3)*nc_Mc(jc_ca))/   &
               mw_electrolyte(jcano3)
       endif

       if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
          electrolyte(jcacl2,jp,ibin)  = (xeq_c(jc_ca) *na_Ma(ja_cl) +   &
               xeq_a(ja_cl) *nc_Mc(jc_ca))/   &
               mw_electrolyte(jcacl2)
       endif

       if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
          electrolyte(jcamsa2,jp,ibin) = (xeq_c(jc_ca) *na_Ma(ja_msa) +   &
               xeq_a(ja_msa) *nc_Mc(jc_ca))/   &
               mw_electrolyte(jcamsa2)
       endif

       electrolyte(jh2so4, jp,ibin) = 0.0

       if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
          electrolyte(jhno3,jp,ibin)     = (xeq_c(jc_h)  *na_Ma(ja_no3) +   &
               xeq_a(ja_no3)*nc_Mc(jc_h))/   &
               mw_electrolyte(jhno3)
       endif

       if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
          electrolyte(jhcl,jp,ibin)    = (xeq_c(jc_h) *na_Ma(ja_cl) +   &
               xeq_a(ja_cl)*nc_Mc(jc_h))/   &
               mw_electrolyte(jhcl)
       endif

       if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
          electrolyte(jmsa,jp,ibin)    = (xeq_c(jc_h) *na_Ma(ja_msa) +   &
               xeq_a(ja_msa)*nc_Mc(jc_h))/   &
               mw_electrolyte(jmsa)
       endif

       

    elseif(icase.eq.2)then 

       store(imsa_a) = aer(imsa_a,jp,ibin)
       store(ica_a)  = aer(ica_a, jp,ibin)

       call form_camsa2(store,jp,ibin,electrolyte)

       sum_na_nh4 = aer(ina_a,jp,ibin) + aer(inh4_a,jp,ibin)

       if(sum_na_nh4 .gt. 0.0)then
          f_na  = aer(ina_a,jp,ibin)/sum_na_nh4
          f_nh4 = aer(inh4_a,jp,ibin)/sum_na_nh4
       else
          f_na  = 0.0
          f_nh4 = 0.0
       endif

       
       if(sum_na_nh4 .gt. store(imsa_a))then
          electrolyte(jnamsa,jp,ibin)  = f_na *store(imsa_a)
          electrolyte(jnh4msa,jp,ibin) = f_nh4*store(imsa_a)
          rem_na = max(0.0_r8, aer(ina_a,jp,ibin) - electrolyte(jnamsa,jp,ibin))  
          rem_nh4= max(0.0_r8, aer(inh4_a,jp,ibin)- electrolyte(jnh4msa,jp,ibin)) 
       else
          electrolyte(jnamsa,jp,ibin)  = aer(ina_a,jp,ibin)
          electrolyte(jnh4msa,jp,ibin) = aer(inh4_a,jp,ibin)
          electrolyte(jmsa,jp,ibin)    = max(0.0_r8, store(imsa_a) - sum_na_nh4) 
          rem_nh4 = 0.0  
          rem_na  = 0.0  
       endif


       
       if(aer(iso4_a,jp,ibin).gt.0.0)then
          XT = (rem_nh4 + rem_na)/aer(iso4_a,jp,ibin)
       else
          goto 10
       endif

       if(XT .le. 1.0)then            
          xh = max(0.0_r8, (1.0_r8 - XT))   
          xb = XT
          electrolyte(jh2so4,jp,ibin)   = xh*aer(iso4_a,jp,ibin)
          electrolyte(jnh4hso4,jp,ibin) = xb*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jnahso4,jp,ibin)  = xb*f_na *aer(iso4_a,jp,ibin)
       elseif(XT .le. 1.5)then    
          xb = max(0.0_r8, 3.0_r8 - 2.0_r8*XT) 
          xl = max(0.0_r8, XT - 1.0_r8)     
          electrolyte(jnh4hso4,jp,ibin) = xb*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jnahso4,jp,ibin)  = xb*f_na *aer(iso4_a,jp,ibin)
          electrolyte(jlvcite,jp,ibin)  = xl*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna3hso4,jp,ibin) = xl*f_na *aer(iso4_a,jp,ibin)
       else                       
          xl = max(0.0_r8, 2.0_r8 - XT)     
          xs = max(0.0_r8, 2.0_r8*XT - 3.0_r8) 
          electrolyte(jlvcite,jp,ibin)  = xl*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna3hso4,jp,ibin) = xl*f_na *aer(iso4_a,jp,ibin)
          electrolyte(jnh4so4,jp,ibin)  = xs*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna2so4,jp,ibin)  = xs*f_na *aer(iso4_a,jp,ibin)
       endif

       electrolyte(jhno3,jp,ibin) = aer(ino3_a,jp,ibin)
       electrolyte(jhcl,jp,ibin)  = aer(icl_a,jp,ibin)

    endif
    
    
    
10  sum_dum = 0.0
    
    
    
    
    
    
    
    
    
    
    

    return
  end subroutine ions_to_electrolytes



  
  
  
  
  
  
  
  
  
  
  subroutine MESA_estimate_eleliquid(ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,    &
       nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a,eleliquid)    
    use module_data_mosaic_aero, only: r8,naer,nbin_a_max,nelectrolyte,ncation,    &
         nanion,jliquid,                                                           &
         jh2so4,jhno3,jhcl,jmsa,jlvcite,jnh4no3,jnh4cl,jcamsa2,jcano3,jcacl2,      &
         jnano3,jnacl,jnh4so4,jnh4hso4,jnh4msa,jna2so4,jnahso4,jnamsa,iso4_a,      &
         ja_so4,ja_no3,ja_cl,imsa_a,ja_msa,jc_ca,ina_a,jc_na,inh4_a,jc_nh4,jc_h,   &
         ica_a,ino3_a,icl_a,ja_hso4

    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout) :: XT
    real(r8), intent(in), dimension(Ncation) :: zc,MW_c
    real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
    real(r8), intent(in), dimension(Nanion) :: za,MW_a
    real(r8), intent(inout), dimension(Nanion) :: xeq_a,na_Ma
    real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
    real(r8), intent(inout), dimension(nelectrolyte) :: eleliquid
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer iaer, je, jc, ja, icase, jp
    real(r8) :: store(naer), sum_dum, sum_naza, sum_nczc, sum_na_nh4,   &
         f_nh4, f_na, xh, xb, xl, xs, XT_d, XNa_d, XNH4_d,   &
         xdum, dum, cat_net
    real(r8) :: nc(ncation), na(nanion)
    real(r8) :: dum_ca, dum_no3, dum_cl, cano3, cacl2

    

    
    do iaer =  1, naer
       aer(iaer,jliquid,ibin) = max(0.0d0, aer(iaer,jliquid,ibin))
    enddo


    
    call calculate_XT(ibin,jliquid,XT,aer)

    if(XT .ge. 2.0 .or. XT.lt.0.)then
       icase = 1 
    else
       icase = 2 
    endif


    
    do je = 1, nelectrolyte
       eleliquid(je) = 0.0
    enddo

    
    
    

    jp = jliquid

    if(icase.eq.1)then 

       dum_ca  = aer(ica_a,jp,ibin)
       dum_no3 = aer(ino3_a,jp,ibin)
       dum_cl  = aer(icl_a,jp,ibin)

       cano3   = min(dum_ca, 0.5*dum_no3)
       dum_ca  = max(0.d0, dum_ca - cano3)
       dum_no3 = max(0.d0, dum_no3 - 2.*cano3)

       cacl2   = min(dum_ca, 0.5*dum_cl)
       dum_ca  = max(0.d0, dum_ca - cacl2)
       dum_cl  = max(0.d0, dum_cl - 2.*cacl2)

       na(ja_hso4)= 0.0
       na(ja_so4) = aer(iso4_a,jp,ibin)
       na(ja_no3) = aer(ino3_a,jp,ibin)
       na(ja_cl)  = aer(icl_a, jp,ibin)
       na(ja_msa) = aer(imsa_a,jp,ibin)

       nc(jc_ca)  = aer(ica_a, jp,ibin)
       nc(jc_na)  = aer(ina_a, jp,ibin)
       nc(jc_nh4) = aer(inh4_a,jp,ibin)

       cat_net = (   &
            (2.*na(ja_so4)+na(ja_no3)+na(ja_cl)+na(ja_msa)) -   &
            (2.*nc(jc_ca) +nc(jc_nh4)+nc(jc_na)) )   

       if(cat_net .lt. 0.0)then

          nc(jc_h) = 0.0

       else  

          nc(jc_h) = cat_net

       endif


       
       sum_naza = 0.0
       do ja = 1, nanion
          sum_naza = sum_naza + na(ja)*za(ja)
       enddo

       sum_nczc = 0.0
       do jc = 1, ncation
          sum_nczc = sum_nczc + nc(jc)*zc(jc)
       enddo

       if(sum_naza .eq. 0. .or. sum_nczc .eq. 0.)then
          write(6,*)'ionic concentrations are zero in ibin', ibin
          write(6,*)'sum_naza = ', sum_naza
          write(6,*)'sum_nczc = ', sum_nczc
          return
       endif

       do ja = 1, nanion
          xeq_a(ja) = na(ja)*za(ja)/sum_naza
       enddo

       do jc = 1, ncation
          xeq_c(jc) = nc(jc)*zc(jc)/sum_nczc
       enddo

       na_Ma(ja_so4) = na(ja_so4) *MW_a(ja_so4)
       na_Ma(ja_no3) = na(ja_no3) *MW_a(ja_no3)
       na_Ma(ja_cl)  = na(ja_cl)  *MW_a(ja_cl)
       na_Ma(ja_hso4)= na(ja_hso4)*MW_a(ja_hso4)
       na_Ma(ja_msa) = na(ja_msa) *MW_a(ja_msa)

       nc_Mc(jc_ca)  = nc(jc_ca) *MW_c(jc_ca)
       nc_Mc(jc_na)  = nc(jc_na) *MW_c(jc_na)
       nc_Mc(jc_nh4) = nc(jc_nh4)*MW_c(jc_nh4)
       nc_Mc(jc_h)   = nc(jc_h)  *MW_c(jc_h)


       
       eleliquid(jna2so4) = (xeq_c(jc_na) *na_Ma(ja_so4) +   &
            xeq_a(ja_so4)*nc_Mc(jc_na))/   &
            mw_electrolyte(jna2so4)

       eleliquid(jnahso4) = (xeq_c(jc_na) *na_Ma(ja_hso4) +   &
            xeq_a(ja_hso4)*nc_Mc(jc_na))/   &
            mw_electrolyte(jnahso4)

       eleliquid(jnamsa)  = (xeq_c(jc_na) *na_Ma(ja_msa) +   &
            xeq_a(ja_msa)*nc_Mc(jc_na))/   &
            mw_electrolyte(jnamsa)

       eleliquid(jnano3)  = (xeq_c(jc_na) *na_Ma(ja_no3) +   &
            xeq_a(ja_no3)*nc_Mc(jc_na))/   &
            mw_electrolyte(jnano3)

       eleliquid(jnacl)   = (xeq_c(jc_na) *na_Ma(ja_cl) +   &
            xeq_a(ja_cl) *nc_Mc(jc_na))/   &
            mw_electrolyte(jnacl)

       eleliquid(jnh4so4) = (xeq_c(jc_nh4)*na_Ma(ja_so4) +   &
            xeq_a(ja_so4)*nc_Mc(jc_nh4))/   &
            mw_electrolyte(jnh4so4)

       eleliquid(jnh4hso4)= (xeq_c(jc_nh4)*na_Ma(ja_hso4) +   &
            xeq_a(ja_hso4)*nc_Mc(jc_nh4))/   &
            mw_electrolyte(jnh4hso4)

       eleliquid(jnh4msa) = (xeq_c(jc_nh4) *na_Ma(ja_msa) +   &
            xeq_a(ja_msa)*nc_Mc(jc_nh4))/   &
            mw_electrolyte(jnh4msa)

       eleliquid(jnh4no3) = (xeq_c(jc_nh4)*na_Ma(ja_no3) +   &
            xeq_a(ja_no3)*nc_Mc(jc_nh4))/   &
            mw_electrolyte(jnh4no3)

       eleliquid(jnh4cl)  = (xeq_c(jc_nh4)*na_Ma(ja_cl) +   &
            xeq_a(ja_cl) *nc_Mc(jc_nh4))/   &
            mw_electrolyte(jnh4cl)

       eleliquid(jcamsa2) = (xeq_c(jc_ca) *na_Ma(ja_msa) +   &
            xeq_a(ja_msa)*nc_Mc(jc_ca))/   &
            mw_electrolyte(jcamsa2)

       eleliquid(jcano3)  = (xeq_c(jc_ca) *na_Ma(ja_no3) +   &
            xeq_a(ja_no3)*nc_Mc(jc_ca))/   &
            mw_electrolyte(jcano3)

       eleliquid(jcacl2)  = (xeq_c(jc_ca) *na_Ma(ja_cl) +   &
            xeq_a(ja_cl) *nc_Mc(jc_ca))/   &
            mw_electrolyte(jcacl2)

       eleliquid(jh2so4)  = (xeq_c(jc_h)   *na_Ma(ja_hso4) +   &
            xeq_a(ja_hso4)*nc_Mc(jc_h))/   &
            mw_electrolyte(jh2so4)

       eleliquid(jhno3)   = (xeq_c(jc_h)  *na_Ma(ja_no3) +   &
            xeq_a(ja_no3)*nc_Mc(jc_h))/   &
            mw_electrolyte(jhno3)

       eleliquid(jhcl)    = (xeq_c(jc_h) *na_Ma(ja_cl) +   &
            xeq_a(ja_cl)*nc_Mc(jc_h))/   &
            mw_electrolyte(jhcl)

       eleliquid(jmsa)    = (xeq_c(jc_h)  *na_Ma(ja_msa) +   &
            xeq_a(ja_msa)*nc_Mc(jc_h))/   &
            mw_electrolyte(jmsa)

       

    elseif(icase.eq.2)then 

       jp = jliquid

       store(iso4_a) = aer(iso4_a,jp,ibin)
       store(imsa_a) = aer(imsa_a,jp,ibin)
       store(inh4_a) = aer(inh4_a,jp,ibin)
       store(ina_a)  = aer(ina_a, jp,ibin)
       store(ica_a)  = aer(ica_a, jp,ibin)

       call form_camsa2(store,jp,ibin,electrolyte)

       sum_na_nh4 = store(ina_a) + store(inh4_a)
       if(sum_na_nh4 .gt. 0.0)then
          f_nh4 = store(inh4_a)/sum_na_nh4
          f_na  = store(ina_a)/sum_na_nh4
       else
          f_nh4 = 0.0
          f_na  = 0.0
       endif

       
       if(sum_na_nh4 .gt. store(imsa_a))then
          eleliquid(jnh4msa) = f_nh4*store(imsa_a)
          eleliquid(jnamsa)  = f_na *store(imsa_a)
          store(inh4_a)= store(inh4_a)-eleliquid(jnh4msa) 
          store(ina_a) = store(ina_a) -eleliquid(jnamsa)  
       else
          eleliquid(jnh4msa) = store(inh4_a)
          eleliquid(jnamsa)  = store(ina_a)
          eleliquid(jmsa)    = store(imsa_a) - sum_na_nh4
          store(inh4_a)= 0.0  
          store(ina_a) = 0.0  
       endif

       if(store(iso4_a).eq.0.0)goto 10

       XT_d  = XT
       XNa_d = 1. + 0.5*store(ina_a)/store(iso4_a)
       xdum = store(iso4_a) - store(inh4_a)

       dum = ( (2.*store(iso4_a)) -   &
            (store(ina_a)) )
       if(store(inh4_a) .gt. 0.0 .and. dum .gt. 0.0)then
          XNH4_d = 2.*store(inh4_a)/   &
               (2.*store(iso4_a) - store(ina_a))
       else
          XNH4_d = 0.0
       endif


       IF(store(inh4_a) .gt. 0.0)THEN
          if(XT_d .ge. XNa_d)then
             eleliquid(jna2so4) = 0.5*store(ina_a)

             if(XNH4_d .ge. 5./3.)then
                eleliquid(jnh4so4) = 1.5*store(ina_a)   &
                     - 3.*xdum - store(inh4_a)
                eleliquid(jlvcite) = 2.*xdum + store(inh4_a)   &
                     - store(ina_a)
             elseif(XNH4_d .ge. 1.5)then
                eleliquid(jnh4so4) = store(inh4_a)/5.
                eleliquid(jlvcite) = store(inh4_a)/5.
             elseif(XNH4_d .ge. 1.0)then
                eleliquid(jnh4so4) = store(inh4_a)/6.
                eleliquid(jlvcite) = store(inh4_a)/6.
                eleliquid(jnh4hso4)= store(inh4_a)/6.
             endif

          elseif(XT_d .gt. 1.0)then
             eleliquid(jnh4so4)  = store(inh4_a)/6.
             eleliquid(jlvcite)  = store(inh4_a)/6.
             eleliquid(jnh4hso4) = store(inh4_a)/6.
             eleliquid(jna2so4)  = store(ina_a)/3.
             eleliquid(jnahso4)  = store(ina_a)/3.
          elseif(XT_d .le. 1.0)then
             eleliquid(jna2so4)  = store(ina_a)/4.
             eleliquid(jnahso4)  = store(ina_a)/2.
             eleliquid(jlvcite)  = store(inh4_a)/6.
             eleliquid(jnh4hso4) = store(inh4_a)/2.
          endif

       ELSE

          if(XT_d .gt. 1.0)then
             eleliquid(jna2so4) = store(ina_a) - store(iso4_a)
             eleliquid(jnahso4) = 2.*store(iso4_a) -   &
                  store(ina_a)
          else
             eleliquid(jna2so4) = store(ina_a)/4.
             eleliquid(jnahso4) = store(ina_a)/2.
          endif


       ENDIF



    endif
    


10  return
  end subroutine MESA_estimate_eleliquid



  
  
  
  
  
  
  subroutine MESA_dissolve_small_salt(ibin,js,aer,electrolyte)

    use module_data_mosaic_aero, only:r8,naer,nbin_a_max,nelectrolyte,jsolid,      &
         jliquid,                                                                  &
         jh2so4,jhno3,jhcl,jlvcite,jnh4no3,jnh4cl,jcamsa2,jcano3,jcacl2,jnano3,    &
         jnacl,jnh4so4,jnh4hso4,jnh4msa,jna2so4,jnahso4,jnamsa,iso4_a,ina_a,       &
         inh4_a,jna3hso4,jcaso4,jcaco3,ica_a,ino3_a,icl_a

    implicit none

    
    integer, intent(in) :: ibin, js
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer :: jp

    jp = jsolid


    if(js .eq. jnh4so4)then
       aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
            2.*electrolyte(js,jsolid,ibin)
       aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
            electrolyte(jnh4cl,jp,ibin)  +   &
            2.*electrolyte(jnh4so4,jp,ibin) +   &
            3.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jnh4msa,jp,ibin)

       aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
            electrolyte(jna2so4,jp,ibin) +   &
            2.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnh4so4,jp,ibin) +   &
            2.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jh2so4,jp,ibin)
       return
    endif


    if(js .eq. jlvcite)then
       aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
            3.*electrolyte(js,jsolid,ibin)
       aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
            2.*electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
            electrolyte(jnh4cl,jp,ibin)  +   &
            2.*electrolyte(jnh4so4,jp,ibin) +   &
            3.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jnh4msa,jp,ibin)

       aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
            electrolyte(jna2so4,jp,ibin) +   &
            2.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnh4so4,jp,ibin) +   &
            2.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jh2so4,jp,ibin)
       return
    endif


    if(js .eq. jnh4hso4)then
       aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)
       aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
            electrolyte(jnh4cl,jp,ibin)  +   &
            2.*electrolyte(jnh4so4,jp,ibin) +   &
            3.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jnh4msa,jp,ibin)

       aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
            electrolyte(jna2so4,jp,ibin) +   &
            2.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnh4so4,jp,ibin) +   &
            2.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jh2so4,jp,ibin)
       return
    endif


    if(js .eq. jna2so4)then
       aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
            2.*electrolyte(js,jsolid,ibin)
       aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
            electrolyte(jnacl,jp,ibin)   +   &
            2.*electrolyte(jna2so4,jp,ibin) +   &
            3.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnamsa,jp,ibin)

       aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
            electrolyte(jna2so4,jp,ibin) +   &
            2.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnh4so4,jp,ibin) +   &
            2.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jh2so4,jp,ibin)
       return
    endif


    if(js .eq. jna3hso4)then
       aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
            3.*electrolyte(js,jsolid,ibin)
       aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
            2.*electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
            electrolyte(jnacl,jp,ibin)   +   &
            2.*electrolyte(jna2so4,jp,ibin) +   &
            3.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnamsa,jp,ibin)

       aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
            electrolyte(jna2so4,jp,ibin) +   &
            2.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnh4so4,jp,ibin) +   &
            2.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jh2so4,jp,ibin)
       return
    endif


    if(js .eq. jnahso4)then
       aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)
       aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
            electrolyte(jnacl,jp,ibin)   +   &
            2.*electrolyte(jna2so4,jp,ibin) +   &
            3.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnamsa,jp,ibin)

       aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
            electrolyte(jna2so4,jp,ibin) +   &
            2.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnh4so4,jp,ibin) +   &
            2.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jh2so4,jp,ibin)
       return
    endif


    if(js .eq. jnh4no3)then
       aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)
       aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
            electrolyte(jnh4cl,jp,ibin)  +   &
            2.*electrolyte(jnh4so4,jp,ibin) +   &
            3.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jnh4msa,jp,ibin)

       aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
            2.*electrolyte(jcano3,jp,ibin)  +   &
            electrolyte(jnh4no3,jp,ibin) +   &
            electrolyte(jhno3,jp,ibin)
       return
    endif


    if(js .eq. jnh4cl)then
       aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)
       aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
            electrolyte(jnh4cl,jp,ibin)  +   &
            2.*electrolyte(jnh4so4,jp,ibin) +   &
            3.*electrolyte(jlvcite,jp,ibin) +   &
            electrolyte(jnh4hso4,jp,ibin)+   &
            electrolyte(jnh4msa,jp,ibin)

       aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
            2.*electrolyte(jcacl2,jp,ibin)  +   &
            electrolyte(jnh4cl,jp,ibin)  +   &
            electrolyte(jhcl,jp,ibin)
       return
    endif


    if(js .eq. jnano3)then
       aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)
       aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
            electrolyte(jnacl,jp,ibin)   +   &
            2.*electrolyte(jna2so4,jp,ibin) +   &
            3.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnamsa,jp,ibin)

       aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
            2.*electrolyte(jcano3,jp,ibin)  +   &
            electrolyte(jnh4no3,jp,ibin) +   &
            electrolyte(jhno3,jp,ibin)
       return
    endif


    if(js .eq. jnacl)then
       aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)
       aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
            electrolyte(jnacl,jp,ibin)   +   &
            2.*electrolyte(jna2so4,jp,ibin) +   &
            3.*electrolyte(jna3hso4,jp,ibin)+   &
            electrolyte(jnahso4,jp,ibin) +   &
            electrolyte(jnamsa,jp,ibin)

       aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
            2.*electrolyte(jcacl2,jp,ibin)  +   &
            electrolyte(jnh4cl,jp,ibin)  +   &
            electrolyte(jhcl,jp,ibin)
       return
    endif


    if(js .eq. jcano3)then
       aer(ica_a,jliquid,ibin)  = aer(ica_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)
       aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
            2.*electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
            electrolyte(jcano3,jp,ibin)  +   &
            electrolyte(jcacl2,jp,ibin)  +   &
            electrolyte(jcaco3,jp,ibin)  +   &
            electrolyte(jcamsa2,jp,ibin)

       aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
            2.*electrolyte(jcano3,jp,ibin)  +   &
            electrolyte(jnh4no3,jp,ibin) +   &
            electrolyte(jhno3,jp,ibin)
       return
    endif


    if(js .eq. jcacl2)then
       aer(ica_a,jliquid,ibin) = aer(ica_a,jliquid,ibin) +   &
            electrolyte(js,jsolid,ibin)
       aer(icl_a,jliquid,ibin) = aer(icl_a,jliquid,ibin) +   &
            2.*electrolyte(js,jsolid,ibin)

       electrolyte(js,jsolid,ibin) = 0.0

       aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
            electrolyte(jcano3,jp,ibin)  +   &
            electrolyte(jcacl2,jp,ibin)  +   &
            electrolyte(jcaco3,jp,ibin)  +   &
            electrolyte(jcamsa2,jp,ibin)

       aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
            2.*electrolyte(jcacl2,jp,ibin)  +   &
            electrolyte(jnh4cl,jp,ibin)  +   &
            electrolyte(jhcl,jp,ibin)
       return
    endif

    return
  end subroutine MESA_dissolve_small_salt



  
  
  
  
  
  
  subroutine form_caso4(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,ica_a,jcaso4

    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jcaso4,jp,ibin) = min(store(ica_a),store(iso4_a))
    store(ica_a)  = ( (store(ica_a)) -   &
         (electrolyte(jcaso4,jp,ibin)) )
    store(iso4_a) = ( (store(iso4_a)) -   &
         (electrolyte(jcaso4,jp,ibin)) )
    store(ica_a)  = max(0.d0, store(ica_a))
    store(iso4_a) = max(0.d0, store(iso4_a))

    return
  end subroutine form_caso4



  subroutine form_camsa2(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         imsa_a,ica_a,jcamsa2
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jcamsa2,jp,ibin) = min(store(ica_a),0.5*store(imsa_a))
    store(ica_a)  = ( (store(ica_a)) -   &
         (electrolyte(jcamsa2,jp,ibin)) )
    store(imsa_a) = ( (store(imsa_a)) -   &
         (2.*electrolyte(jcamsa2,jp,ibin)) )
    store(ica_a)  = max(0.d0, store(ica_a))
    store(imsa_a) = max(0.d0, store(imsa_a))

    return
  end subroutine form_camsa2



  
  
  
  
  
  
  
  subroutine aerosolmtc( jaerosolstate, num_a, Dp_wet_a, sigmag_a, P_atm, T_K, kg )
    
    use module_data_mosaic_aero,  only: r8, nbin_a_max, nbin_a, naer, naercomp,             &
         ngas_aerchtot, ngas_volatile, nelectrolyte, ngas_ioa,                              &
         mMODAL, no_aerosol, mUNSTRUCTURED, mSECTIONAL, mSIZE_FRAMEWORK,                    &
         isoa_first, mw_gas, v_molar_gas,                                                   &
         i_gas2bin_uptk_flag, m_gas2bin_uptk_flag,                                          &
         use_cam5mam_accom_coefs, mosaic_vars_aa_type


    implicit none
    
    
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate

    real(r8), intent(in) :: P_atm,T_K
    real(r8), intent(in), dimension(nbin_a_max) :: num_a
    real(r8), intent(in), dimension(nbin_a_max) :: sigmag_a
    real(r8), intent(inout), dimension(nbin_a_max) :: Dp_wet_a
    real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg

    
    integer nghq
    parameter (nghq = 2)         
    integer ibin, iq, iv
    real(r8) :: tworootpi, root2, beta
    parameter (tworootpi = 3.5449077, root2 = 1.4142135, beta = 2.0)
    real(r8) :: cdum, Dp, Dp_avg, Fkn, Kn, lnsg, lnDpgn, lnDp, speed,   &
                sumghq, tmpa
    real(r8) :: xghq(nghq), wghq(nghq)                           
    real(r8) :: accom(ngas_aerchtot)                             
    real(r8) :: freepath(ngas_aerchtot)                          
    real(r8) :: Dg(ngas_aerchtot)                                
    
    
    

    
    tmpa = 0.1
    if ( use_cam5mam_accom_coefs > 0 ) tmpa = 0.65
    accom(1:ngas_aerchtot) = tmpa  














    
    xghq(1) =  0.70710678
    xghq(2) = -0.70710678
    wghq(1) =  0.88622693
    wghq(2) =  0.88622693


    
    
    do iv = 1, ngas_ioa
       speed  = mean_molecular_speed(T_K,mw_gas(iv))     
       Dg(iv) = gas_diffusivity(T_K,P_atm,mw_gas(iv),v_molar_gas(iv)) 
       freepath(iv) = 3.*Dg(iv)/speed                    
    enddo

    
    do iv = isoa_first, ngas_volatile
       speed = mean_molecular_speed(T_K,mw_gas(iv))      
       Dg(iv) = 0.1                                      
       freepath(iv) = 3.*Dg(iv)/speed
    enddo


    do iv = (ngas_volatile+1), ngas_aerchtot
       speed = mean_molecular_speed(t_k,mw_gas(iv))    
       dg(iv) = gas_diffusivity(t_k,p_atm,mw_gas(iv),v_molar_gas(iv)) 
       freepath(iv) = 3.*dg(iv)/speed                  
    enddo


    

    if (mSIZE_FRAMEWORK .eq. mMODAL) then

       
       do 10 ibin = 1, nbin_a

          if(jaerosolstate(ibin) .eq. no_aerosol)goto 10

          lnsg   = log(sigmag_a(ibin))

          
          
          
          
          
          lnDpgn = log(Dp_wet_a(ibin)) - 1.5*lnsg*lnsg

          cdum   = tworootpi*num_a(ibin)*   &
               exp(beta*lnDpgn + 0.5*(beta*lnsg)**2)

          do 20 iv = 1, ngas_aerchtot

             sumghq = 0.0_r8
             do 30 iq = 1, nghq  
                lnDp = lnDpgn + beta*lnsg**2 + root2*lnsg*xghq(iq)
                Dp = exp(lnDp)
                Kn = 2.*freepath(iv)/Dp
                Fkn = fuchs_sutugin(Kn,accom(iv))
                sumghq = sumghq + wghq(iq)*Dp*Fkn/(Dp**beta)
30              continue

                kg(iv,ibin) = cdum*Dg(iv)*sumghq         

20              continue
10     continue
                
    elseif ((mSIZE_FRAMEWORK .eq. mSECTIONAL   ) .or. &
         (mSIZE_FRAMEWORK .eq. mUNSTRUCTURED)) then
       
       
       do 11 ibin = 1, nbin_a
          
          if(jaerosolstate(ibin) .eq. no_aerosol)goto 11
          
          cdum  = 6.283185*Dp_wet_a(ibin)*num_a(ibin)
          
          do 21 iv = 1, ngas_aerchtot
             Kn = 2.*freepath(iv)/Dp_wet_a(ibin)
             Fkn = fuchs_sutugin(Kn,accom(iv))
             kg(iv,ibin) = cdum*Dg(iv)*Fkn              
21           continue
             
11     continue
            
    else
       call wrf_message('Error in the choice of mSIZE_FRAMEWORK')
       call wrf_error_fatal3("<stdin>",3665,&
'Stopping in subr. aerosolmtc')       
    endif


    if (m_gas2bin_uptk_flag <= 0) then
       
       do ibin = 1, nbin_a
          do iv = 1, ngas_aerchtot
             if (i_gas2bin_uptk_flag(iv,ibin) <= 0) kg(iv,ibin) = 0.0
          end do
       end do
    end if

    return
  end subroutine aerosolmtc



  
  
  
  
  
  
  subroutine calc_dry_n_wet_aerosol_props(                          &
     ibin, jaerosolstate, aer, electrolyte, water_a, num_a,         &  
     dens_comp_a, mw_comp_a, dens_aer_mac, mw_aer_mac, ref_index_a, &  
     Dp_dry_a, Dp_wet_a, dp_core_a,                                 &  
     area_dry_a, area_wet_a, mass_dry_a, mass_wet_a,                &  
     vol_dry_a, vol_wet_a, dens_dry_a, dens_wet_a,                  &  
     ri_shell_a, ri_core_a, ri_avg_a                                )  


    use module_data_mosaic_constants,  only: piover4,piover6,third
    use module_data_mosaic_aero,  only: r8,nbin_a_max,naer,nelectrolyte,naercomp,  &
         ngas_soa,no_aerosol,msectional,                                           &
         maeroptic_aero,msize_framework,                                           &
         inh4_a,ina_a,ica_a,ico3_a,imsa_a,icl_a,ino3_a,jtotal,iso4_a,ioc_a,joc,    &
         ibc_a,jbc,ioin_a,join,jh2o, isoa_first, jsoa_first

    use module_data_mosaic_asecthp, only: dcen_sect,isize_of_ibin,itype_of_ibin

    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(in), dimension(nbin_a_max) :: jaerosolstate

    real(r8), intent(in), dimension(nbin_a_max) :: num_a
    real(r8), intent(in), dimension(naer) :: mw_aer_mac,dens_aer_mac
    real(r8), intent(in), dimension(naercomp) :: dens_comp_a,mw_comp_a
    real(r8), intent(inout), dimension(nbin_a_max) :: Dp_dry_a,Dp_wet_a
    real(r8), intent(inout), dimension(nbin_a_max) :: dp_core_a,vol_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: vol_wet_a,dens_wet_a,water_a
    real(r8), intent(inout), dimension(nbin_a_max) :: area_dry_a,area_wet_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a,mass_wet_a
    real(r8), intent(inout), dimension(nbin_a_max) :: dens_dry_a
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    complex, intent(in), dimension(naercomp) :: ref_index_a
    complex, intent(inout), dimension(nbin_a_max) :: ri_avg_a,ri_core_a,ri_shell_a
    
    integer i, iaer, isize, itype, j, jc, je, k
    real(r8) :: aer_H, duma, vol_core, vol_shell, vol_dum
    real(r8),dimension(naercomp) :: comp_a
    complex rixvol_tot, rixvol_core, rixvol_shell


    
    mass_dry_a(ibin) = 0.0                
    vol_dry_a(ibin)  = 0.0                
    area_dry_a(ibin) = 0.0                

    if(jaerosolstate(ibin) .ne. no_aerosol)then

       aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
            aer(ino3_a,jtotal,ibin) +   &
            aer(icl_a,jtotal,ibin)  +   &
            aer(imsa_a,jtotal,ibin) +   &
            2.*aer(ico3_a,jtotal,ibin))-   &
            (2.*aer(ica_a,jtotal,ibin)  +   &
            aer(ina_a,jtotal,ibin)  +   &
            aer(inh4_a,jtotal,ibin))
       aer_H = max(aer_H, 0.0d0)

       do iaer = 1, naer
          mass_dry_a(ibin) = mass_dry_a(ibin) +   &
               aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)     
          vol_dry_a(ibin) = vol_dry_a(ibin) +   &
               aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)          
       enddo
       mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
       vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

       mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15                 
       vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15                   

       
       mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3  
       vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3   

       
       dens_dry_a(ibin) = mass_dry_a(ibin)/vol_dry_a(ibin)                
       dens_wet_a(ibin) = mass_wet_a(ibin)/vol_wet_a(ibin)                

       
       Dp_dry_a(ibin)=(vol_dry_a(ibin)/(piover6*num_a(ibin)))**third      
       Dp_wet_a(ibin)=(vol_wet_a(ibin)/(piover6*num_a(ibin)))**third      

       
       area_dry_a(ibin)= piover4*num_a(ibin)*Dp_dry_a(ibin)**2    
       area_wet_a(ibin)= piover4*num_a(ibin)*Dp_wet_a(ibin)**2    

       
       

       
       
       if (maeroptic_aero <= 0) goto 100

       do je = 1, nelectrolyte
          comp_a(je)=electrolyte(je,jtotal,ibin)*mw_comp_a(je)*1.e-15     
       enddo
       comp_a(joc)  = aer(ioc_a,  jtotal,ibin)*mw_comp_a(joc  )*1.e-15    
       comp_a(jbc)  = aer(ibc_a,  jtotal,ibin)*mw_comp_a(jbc  )*1.e-15    
       comp_a(join) = aer(ioin_a, jtotal,ibin)*mw_comp_a(join )*1.e-15    









       do k = 1, ngas_soa
       j = jsoa_first + k - 1
       i = isoa_first + k - 1
       comp_a(j) = aer(i,jtotal,ibin)*mw_comp_a(j)*1.e-15    
       end do

       comp_a(jh2o) = water_a(ibin)*1.e-3                         

       rixvol_tot   = (0.0,0.0)
       do jc = 1, naercomp
          comp_a(jc) = max( 0.0d0, comp_a(jc) )
          rixvol_tot = rixvol_tot   &
               + ref_index_a(jc)*comp_a(jc)/dens_comp_a(jc)
       enddo
       ri_avg_a(ibin) = rixvol_tot/vol_wet_a(ibin)

       
       
       
       ri_shell_a(ibin) = ri_avg_a(ibin)
       ri_core_a(ibin)  = (0.0,0.0)
       Dp_core_a(ibin)  = 0.0

       
       
       jc = jbc
       rixvol_core  = ref_index_a(jc)*comp_a(jc)/dens_comp_a(jc)
       vol_core = comp_a(jc)/dens_comp_a(jc)
       vol_core = max( 0.0d0, min( vol_core, vol_wet_a(ibin) ) )

       
       
       
       vol_dum = max( 1.0d-22, 1.0d-9*vol_wet_a(ibin) )
       vol_shell = vol_wet_a(ibin) - vol_core
       if (vol_core >= vol_dum) then
          if (vol_shell < vol_dum) then
             ri_shell_a(ibin)  = (0.0,0.0)
             ri_core_a(ibin) = ri_avg_a(ibin)
             Dp_core_a(ibin) = Dp_wet_a(ibin)
          else
             ri_core_a(ibin) = rixvol_core/vol_core
             Dp_core_a(ibin) = Dp_wet_a(ibin)   &
                  * (vol_core/vol_wet_a(ibin))**third

             if (vol_shell >= vol_dum) then
                rixvol_shell = rixvol_tot - rixvol_core
                ri_shell_a(ibin) = rixvol_shell/vol_shell
             else
                ri_shell_a(ibin) = (0.0,0.0)
             endif
          endif
       endif

    else
       

       dens_dry_a(ibin) = 1.0      
       dens_wet_a(ibin) = 1.0      
       
       
       if (msize_framework == msectional) then
          isize = isize_of_ibin(ibin)
          itype = itype_of_ibin(ibin)
          Dp_dry_a(ibin) = dcen_sect(isize,itype)
          Dp_wet_a(ibin) = Dp_dry_a(ibin)
       end if

       ri_avg_a(ibin) = (1.5,0.0)
       ri_shell_a(ibin) = (1.5,0.0)
       ri_core_a(ibin)  = (0.0,0.0)
       Dp_core_a(ibin)  = 0.0

    endif   


100 continue

    return
  end subroutine calc_dry_n_wet_aerosol_props



  
  
  
  
  
  
  subroutine form_electrolytes(jp,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)
    use module_data_mosaic_aero, only: r8,naer,nbin_a_max,           &
         ngas_aerchtot,ngas_volatile,nelectrolyte,jsolid,            &
         imsa_a,iso4_a,ica_a,ina_a,inh4_a,ino3_a,icl_a,ico3_a

    implicit none

    
    integer, intent(in) :: ibin, jp
    real(r8), intent(inout) :: XT
    real(r8), intent(inout) :: tot_cl_in
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(ngas_volatile) :: total_species
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer i, iXT_case, j, je
    real(r8) :: sum_dum, XNa_prime, XNH4_prime, XT_prime
    real(r8) :: store(naer)

    
    
    
    


    

    if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
       XT   = ( aer(inh4_a,jp,ibin) +   &
            aer(ina_a,jp,ibin)  +   &
            2.*aer(ica_a,jp,ibin) )/   &
            (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
    else
       XT   = -1.0
    endif





    if(XT .ge. 2.0 .or. XT.lt.0.)then         
       iXT_case = 1       
    else
       iXT_case = 2       
    endif

    
    
    
    store(iso4_a) = aer(iso4_a,jp,ibin)
    store(ino3_a) = aer(ino3_a,jp,ibin)
    store(icl_a)  = aer(icl_a, jp,ibin)
    store(imsa_a) = aer(imsa_a,jp,ibin)
    store(ico3_a) = aer(ico3_a,jp,ibin)
    store(inh4_a) = aer(inh4_a,jp,ibin)
    store(ina_a)  = aer(ina_a, jp,ibin)
    store(ica_a)  = aer(ica_a, jp,ibin)

    do j=1,nelectrolyte
       electrolyte(j,jp,ibin) = 0.0
    enddo

    
    
    
    if(iXT_case.eq.1)then

       
       call form_caso4(store,jp,ibin,electrolyte)
       call form_camsa2(store,jp,ibin,electrolyte)
       call form_na2so4(store,jp,ibin,electrolyte)
       call form_namsa(store,jp,ibin,electrolyte)
       call form_cano3(store,jp,ibin,electrolyte)
       call form_nano3(store,jp,ibin,electrolyte)
       call form_nacl(store,jp,ibin,aer,gas,electrolyte,total_species,tot_cl_in)
       call form_cacl2(store,jp,ibin,electrolyte)
       call form_caco3(store,jp,ibin,aer,electrolyte)
       call form_nh4so4(store,jp,ibin,electrolyte)
       call form_nh4msa(store,jp,ibin,electrolyte)
       call form_nh4no3(store,jp,ibin,electrolyte)
       call form_nh4cl(store,jp,ibin,electrolyte)
       call form_msa(store,jp,ibin,electrolyte)

       if(jp .eq. jsolid)then
          call degas_hno3(store,jp,ibin,aer,gas,electrolyte)
          call degas_hcl(store,jp,ibin,aer,gas,electrolyte)
          call degas_nh3(store,jp,ibin,aer,gas)
       else
          call form_hno3(store,jp,ibin,electrolyte)
          call form_hcl(store,jp,ibin,electrolyte)
          call degas_nh3(store,jp,ibin,aer,gas)
       endif



    elseif(iXT_case.eq.2)then

       

       call form_caso4(store,jp,ibin,electrolyte)
       call form_camsa2(store,jp,ibin,electrolyte)
       call form_namsa(store,jp,ibin,electrolyte)
       call form_nh4msa(store,jp,ibin,electrolyte)
       call form_msa(store,jp,ibin,electrolyte)

       if(store(iso4_a).eq.0.0)goto 10


       XT_prime =(store(ina_a)+store(inh4_a))/   &
            store(iso4_a)
       XNa_prime=0.5*store(ina_a)/store(iso4_a) + 1.

       if(XT_prime.ge.XNa_prime)then
          call form_na2so4(store,jp,ibin,electrolyte)
          XNH4_prime = 0.0
          if(store(iso4_a).gt.1.e-15)then
             XNH4_prime = store(inh4_a)/store(iso4_a)
          endif

          if(XNH4_prime .ge. 1.5)then
             call form_nh4so4_lvcite(store,jp,ibin,electrolyte)
          else
             call form_lvcite_nh4hso4(store,jp,ibin,electrolyte)
          endif

       elseif(XT_prime.ge.1.)then
          call form_nh4hso4(store,jp,ibin,electrolyte)
          call form_na2so4_nahso4(store,jp,ibin,electrolyte)
       elseif(XT_prime.lt.1.)then
          call form_nahso4(store,jp,ibin,electrolyte)
          call form_nh4hso4(store,jp,ibin,electrolyte)
          call form_h2so4(store,jp,ibin,electrolyte)
       endif

10     if(jp .eq. jsolid)then
          call degas_hno3(store,jp,ibin,aer,gas,electrolyte)
          call degas_hcl(store,jp,ibin,aer,gas,electrolyte)
          call degas_nh3(store,jp,ibin,aer,gas)
       else
          call form_hno3(store,jp,ibin,electrolyte)
          call form_hcl(store,jp,ibin,electrolyte)
          call degas_nh3(store,jp,ibin,aer,gas)
       endif

    endif 


    
    call electrolytes_to_ions(jp, ibin,aer,electrolyte)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    return
  end subroutine form_electrolytes



   subroutine form_na2so4(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,ina_a,jna2so4
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store(naer)
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    electrolyte(jna2so4,jp,ibin) = min(.5*store(ina_a),   &
         store(iso4_a))
    store(ina_a) =( (store(ina_a)) -   &
         (2.*electrolyte(jna2so4,jp,ibin)) )
    store(iso4_a)=( (store(iso4_a)) -   &
         (electrolyte(jna2so4,jp,ibin)) )
    store(ina_a) =max(0.d0, store(ina_a))
    store(iso4_a)=max(0.d0, store(iso4_a))

    return
  end subroutine form_na2so4



  subroutine form_nahso4(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,ina_a,jnahso4

    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnahso4,jp,ibin) = min(store(ina_a),   &
         store(iso4_a))
    store(ina_a)  = ( (store(ina_a)) -   &
         (electrolyte(jnahso4,jp,ibin)) )
    store(iso4_a) = ( (store(iso4_a)) -   &
         (electrolyte(jnahso4,jp,ibin)) )
    store(ina_a)  = max(0.d0, store(ina_a))
    store(iso4_a) = max(0.d0, store(iso4_a))

    return
  end subroutine form_nahso4



  subroutine form_namsa(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         imsa_a,ina_a,jnamsa
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer)  :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnamsa,jp,ibin) = min(store(ina_a),   &
         store(imsa_a))
    store(ina_a)  = ( (store(ina_a)) -   &
         (electrolyte(jnamsa,jp,ibin)) )
    store(imsa_a) = ( (store(imsa_a)) -   &
         (electrolyte(jnamsa,jp,ibin)) )
    store(ina_a)  = max(0.d0, store(ina_a))
    store(imsa_a) = max(0.d0, store(imsa_a))

    return
  end subroutine form_namsa



  subroutine form_nano3(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ino3_a,ina_a,jnano3
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnano3,jp,ibin)=min(store(ina_a),store(ino3_a))
    store(ina_a)  = ( (store(ina_a)) -   &
         (electrolyte(jnano3,jp,ibin)) )
    store(ino3_a) = ( (store(ino3_a)) -   &
         (electrolyte(jnano3,jp,ibin)) )
    store(ina_a)  = max(0.d0, store(ina_a))
    store(ino3_a) = max(0.d0, store(ino3_a))

    return
  end subroutine form_nano3



  subroutine form_cano3(store,jp,ibin,electrolyte)        
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ino3_a,ica_a,jcano3
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jcano3,jp,ibin) = min(store(ica_a),0.5*store(ino3_a))

    store(ica_a)  = ( (store(ica_a)) -   &
         (electrolyte(jcano3,jp,ibin)) )
    store(ino3_a) = ( (store(ino3_a)) -   &
         (2.*electrolyte(jcano3,jp,ibin)) )
    store(ica_a)  = max(0.d0, store(ica_a))
    store(ino3_a) = max(0.d0, store(ino3_a))

    return
  end subroutine form_cano3



  subroutine form_cacl2(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         icl_a,ica_a,jcacl2

    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jcacl2,jp,ibin) = min(store(ica_a),0.5*store(icl_a))

    store(ica_a)  = ( (store(ica_a)) -   &
         (electrolyte(jcacl2,jp,ibin)) )
    store(icl_a)  = ( (store(icl_a)) -   &
         (2.*electrolyte(jcacl2,jp,ibin)) )
    store(ica_a)  = max(0.d0, store(ica_a))
    store(icl_a)  = max(0.d0, store(icl_a))

    return
  end subroutine form_cacl2

  
  subroutine form_caco3(store,jp,ibin,aer,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,jsolid,     &
         jtotal,                                                                   &
         ica_a,jcaco3,ico3_a

    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    if(jp.eq.jtotal .or. jp.eq.jsolid)then
       electrolyte(jcaco3,jp,ibin) = store(ica_a)

       aer(ico3_a,jp,ibin)= electrolyte(jcaco3,jp,ibin)   

       store(ica_a) = 0.0
       store(ico3_a)= 0.0
    endif

    return
  end subroutine form_caco3
  
  
  
  subroutine form_nacl(store,jp,ibin,aer,gas,electrolyte,total_species,tot_cl_in)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,ngas_volatile,jtotal,jsolid,jliquid,                        &
         ina_a,jnacl,icl_a,ihcl_g

    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(ngas_volatile) :: total_species
    real(r8), intent(inout) :: tot_cl_in
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    

    electrolyte(jnacl,jp,ibin) = store(ina_a)

    store(ina_a) = 0.0
    store(icl_a) = ( (store(icl_a)) -   &
         (electrolyte(jnacl,jp,ibin)) )

    if(store(icl_a) .lt. 0.)then                          
       aer(icl_a,jp,ibin)= aer(icl_a,jp,ibin)- store(icl_a)       

       if(jp .ne. jtotal)then
          aer(icl_a,jtotal,ibin)= aer(icl_a,jliquid,ibin)+ &      
               aer(icl_a,jsolid,ibin)
       endif

       gas(ihcl_g) = gas(ihcl_g) + store(icl_a)                   

       if(gas(ihcl_g) .lt. 0.0)then
          total_species(ihcl_g) = total_species(ihcl_g) - gas(ihcl_g)     
          tot_cl_in = tot_cl_in - gas(ihcl_g)                             
       endif

       gas(ihcl_g) = max(0.d0, gas(ihcl_g))                               
       store(icl_a) = 0.                                          

    endif

    store(icl_a) = max(0.d0, store(icl_a))

    return
  end subroutine form_nacl



  subroutine form_nh4so4(store,jp,ibin,electrolyte)       
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,inh4_a,jnh4so4
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnh4so4,jp,ibin)= min(.5*store(inh4_a),   &
         store(iso4_a))
    store(inh4_a)= ( (store(inh4_a)) -   &
         (2.*electrolyte(jnh4so4,jp,ibin)) )
    store(iso4_a)= ( (store(iso4_a)) -   &
         (electrolyte(jnh4so4,jp,ibin)) )
    store(inh4_a) = max(0.d0, store(inh4_a))
    store(iso4_a) = max(0.d0, store(iso4_a))

    return
  end subroutine form_nh4so4



  subroutine form_nh4hso4(store,jp,ibin,electrolyte)      
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,inh4_a,jnh4hso4
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnh4hso4,jp,ibin) = min(store(inh4_a),   &
         store(iso4_a))
    store(inh4_a)= ( (store(inh4_a)) -   &
         (electrolyte(jnh4hso4,jp,ibin)) )
    store(iso4_a)= ( (store(iso4_a)) -   &
         (electrolyte(jnh4hso4,jp,ibin)) )
    store(inh4_a) = max(0.d0, store(inh4_a))
    store(iso4_a) = max(0.d0, store(iso4_a))

    return
  end subroutine form_nh4hso4



  subroutine form_nh4msa(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         imsa_a,inh4_a,jnh4msa
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnh4msa,jp,ibin) = min(store(inh4_a),   &
         store(imsa_a))
    store(inh4_a) = ( (store(inh4_a)) -   &
         (electrolyte(jnh4msa,jp,ibin)) )
    store(imsa_a) = ( (store(imsa_a)) -   &
         (electrolyte(jnh4msa,jp,ibin)) )
    store(inh4_a) = max(0.d0, store(inh4_a))
    store(imsa_a) = max(0.d0, store(imsa_a))

    return
  end subroutine form_nh4msa



  subroutine form_nh4cl(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         icl_a,inh4_a,jnh4cl
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnh4cl,jp,ibin) = min(store(inh4_a),   &
         store(icl_a))
    store(inh4_a) = ( (store(inh4_a)) -   &
         (electrolyte(jnh4cl,jp,ibin)) )
    store(icl_a)  = ( (store(icl_a)) -   &
         (electrolyte(jnh4cl,jp,ibin)) )
    store(inh4_a) = max(0.d0, store(inh4_a))
    store(icl_a)  = max(0.d0, store(icl_a))

    return
  end subroutine form_nh4cl



  subroutine form_nh4no3(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ino3_a,inh4_a,jnh4no3
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnh4no3,jp,ibin) = min(store(inh4_a),   &
         store(ino3_a))
    store(inh4_a) = ( (store(inh4_a)) -   &
         (electrolyte(jnh4no3,jp,ibin)) )
    store(ino3_a) = ( (store(ino3_a)) -   &
         (electrolyte(jnh4no3,jp,ibin)) )
    store(inh4_a) = max(0.d0, store(inh4_a))
    store(ino3_a) = max(0.d0, store(ino3_a))

    return
  end subroutine form_nh4no3



  subroutine form_nh4so4_lvcite(store,jp,ibin,electrolyte) 
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,inh4_a,jnh4so4,jlvcite
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jnh4so4,jp,ibin)= ( (2.*store(inh4_a)) -   &
         (3.*store(iso4_a)) )
    electrolyte(jlvcite,jp,ibin)= ( (2.*store(iso4_a)) -   &
         (store(inh4_a)) )
    electrolyte(jnh4so4,jp,ibin)= max(0.d0,   &
         electrolyte(jnh4so4,jp,ibin))
    electrolyte(jlvcite,jp,ibin)= max(0.d0,   &
         electrolyte(jlvcite,jp,ibin))
    store(inh4_a) = 0.
    store(iso4_a) = 0.

    return
  end subroutine form_nh4so4_lvcite



  subroutine form_lvcite_nh4hso4(store,jp,ibin,electrolyte) 
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,inh4_a,jlvcite,jnh4hso4
    implicit none

    
    integer, intent(in) ::  jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jlvcite,jp,ibin) = ( (store(inh4_a)) -   &
         (store(iso4_a)) )
    electrolyte(jnh4hso4,jp,ibin)= ( (3.*store(iso4_a)) -   &
         (2.*store(inh4_a)) )
    electrolyte(jlvcite,jp,ibin) = max(0.d0,   &
         electrolyte(jlvcite,jp,ibin))
    electrolyte(jnh4hso4,jp,ibin)= max(0.d0,   &
         electrolyte(jnh4hso4,jp,ibin))
    store(inh4_a) = 0.
    store(iso4_a) = 0.

    return
  end subroutine form_lvcite_nh4hso4



  subroutine form_na2so4_nahso4(store,jp,ibin,electrolyte) 
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,ina_a,jna2so4,jnahso4
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jna2so4,jp,ibin)= ( (store(ina_a)) -   &
         (store(iso4_a)) )
    electrolyte(jnahso4,jp,ibin)= ( (2.*store(iso4_a))-   &
         (store(ina_a)) )
    electrolyte(jna2so4,jp,ibin)= max(0.d0,   &
         electrolyte(jna2so4,jp,ibin))
    electrolyte(jnahso4,jp,ibin)= max(0.d0,   &
         electrolyte(jnahso4,jp,ibin))
    store(ina_a)  = 0.
    store(iso4_a) = 0.



    return
  end subroutine form_na2so4_nahso4



  subroutine form_h2so4(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         iso4_a,jh2so4
    implicit none

    
    integer, intent(in) ::  jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jh2so4,jp,ibin) = max(0.0d0, store(iso4_a))
    store(iso4_a) = 0.0

    return
  end subroutine form_h2so4



  subroutine form_msa(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         imsa_a,jmsa
    implicit none

    
    integer, intent(in) ::  jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jmsa,jp,ibin) = max(0.0d0, store(imsa_a))
    store(imsa_a) = 0.0

    return
  end subroutine form_msa



  subroutine form_hno3(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ino3_a,jhno3
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jhno3,jp,ibin) = max(0.0d0, store(ino3_a))
    store(ino3_a) = 0.0

    return
  end subroutine form_hno3



  subroutine form_hcl(store,jp,ibin,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         icl_a,jhcl
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    electrolyte(jhcl,jp,ibin) = max(0.0d0, store(icl_a))
    store(icl_a) = 0.0

    return
  end subroutine form_hcl



 subroutine degas_hno3(store,jp,ibin,aer,gas,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jtotal,jliquid,jsolid,                                      &
         ino3_a,ihno3_g,jhno3

    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout), dimension(naer) :: store
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    store(ino3_a) = max(0.0d0, store(ino3_a))
    gas(ihno3_g) = gas(ihno3_g) + store(ino3_a)
    aer(ino3_a,jp,ibin) = ( (aer(ino3_a,jp,ibin)) -   &
         (store(ino3_a)) )
    aer(ino3_a,jp,ibin) = max(0.0d0,aer(ino3_a,jp,ibin))

    
    if(jp .ne. jtotal)then
       aer(ino3_a,jtotal,ibin) = aer(ino3_a,jsolid, ibin) +   &
            aer(ino3_a,jliquid,ibin)
    endif

    electrolyte(jhno3,jp,ibin) = 0.0
    store(ino3_a) = 0.0

    return
  end subroutine degas_hno3



  subroutine degas_hcl(store,jp,ibin,aer,gas,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jtotal,jliquid,jsolid,                                      &
         icl_a,ihcl_g,jhcl
    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout) :: store(naer)
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    store(icl_a) = max(0.0d0, store(icl_a))
    gas(ihcl_g) = gas(ihcl_g) + store(icl_a)
    aer(icl_a,jp,ibin) = ( (aer(icl_a,jp,ibin)) -   &
         (store(icl_a)) )
    aer(icl_a,jp,ibin) = max(0.0d0,aer(icl_a,jp,ibin))

    
    if(jp .ne. jtotal)then
       aer(icl_a,jtotal,ibin) = aer(icl_a,jsolid, ibin) +   &
            aer(icl_a,jliquid,ibin)
    endif

    electrolyte(jhcl,jp,ibin) = 0.0
    store(icl_a) = 0.0

    return
  end subroutine degas_hcl



  subroutine degas_nh3(store,jp,ibin,aer,gas)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jtotal,jliquid,jsolid,                                      &
         inh3_g,inh4_a

    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(inout) :: store(naer)
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer

    store(inh4_a) = max(0.0d0, store(inh4_a))
    gas(inh3_g) = gas(inh3_g) + store(inh4_a)
    aer(inh4_a,jp,ibin) = ( (aer(inh4_a,jp,ibin)) -   &
         (store(inh4_a)) )
    aer(inh4_a,jp,ibin) = max(0.0d0,aer(inh4_a,jp,ibin))

    
    if(jp .ne. jtotal)then
       aer(inh4_a,jtotal,ibin)= aer(inh4_a,jsolid, ibin) +   &
            aer(inh4_a,jliquid,ibin)
    endif

    store(inh4_a) = 0.0

    return
  end subroutine degas_nh3



  
  
  
  
  
  
  
  
  subroutine absorb_tiny_nh4no3(ibin,aer,gas,electrolyte,delta_nh3_max,            &
       delta_hno3_max,electrolyte_sum)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jtotal,jliquid,jsolid,                                      &
         inh4_a,ino3_a,inh3_g,ihno3_g
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(nbin_a_max) :: delta_nh3_max,delta_hno3_max
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(3,nbin_a_max) :: electrolyte_sum
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer je
    real(r8) :: small_aer, small_gas, small_amt



    
    electrolyte_sum(jtotal,ibin) = 0.0
    do je = 1, nelectrolyte
       electrolyte_sum(jtotal,ibin) = electrolyte_sum(jtotal,ibin) + &
            electrolyte(je,jtotal,ibin)
    enddo
    


    small_gas = 0.01 * min(delta_nh3_max(ibin),delta_hno3_max(ibin))
    small_aer = 0.01 * electrolyte_sum(jtotal,ibin)
    if(small_aer .eq. 0.0)small_aer = small_gas

    small_amt = min(small_gas, small_aer)

    aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) + small_amt
    aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) + small_amt

    
    aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
         aer(inh4_a,jliquid,ibin)
    aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
         aer(ino3_a,jliquid,ibin)

    
    gas(inh3_g)  = ((gas(inh3_g)) - (small_amt))
    gas(ihno3_g) = ((gas(ihno3_g)) - (small_amt))

    return
  end subroutine absorb_tiny_nh4no3



  
  
  subroutine absorb_tiny_nh4cl(ibin,aer,gas,electrolyte,delta_nh3_max,             &
       delta_hcl_max,electrolyte_sum)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jtotal,jliquid,jsolid,                                      &
         inh4_a,icl_a,inh3_g,ihcl_g
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(nbin_a_max) :: delta_nh3_max,delta_hcl_max
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(3,nbin_a_max) :: electrolyte_sum
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer je
    real(r8) :: small_aer, small_gas, small_amt


    
    electrolyte_sum(jtotal,ibin) = 0.0
    do je = 1, nelectrolyte
       electrolyte_sum(jtotal,ibin) = electrolyte_sum(jtotal,ibin) + &
            electrolyte(je,jtotal,ibin)
    enddo
    



    small_gas = 0.01 * min(delta_nh3_max(ibin), delta_hcl_max(ibin))
    small_aer = 0.01 * electrolyte_sum(jtotal,ibin)
    if(small_aer .eq. 0.0)small_aer = small_gas

    small_amt = min(small_gas, small_aer)

    aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) + small_amt
    aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin)  + small_amt

    
    aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
         aer(inh4_a,jliquid,ibin)
    aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
         aer(icl_a,jliquid,ibin)

    
    gas(inh3_g) = ((gas(inh3_g)) - (small_amt))
    gas(ihcl_g) = ((gas(ihcl_g)) - (small_amt))

    return
  end subroutine absorb_tiny_nh4cl



  
  
  subroutine absorb_tiny_hno3(ibin,aer,gas,delta_hno3_max)        
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jliquid,jsolid,jtotal,                                      &
         icl_a,ino3_a,ihno3_g,ihcl_g
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(nbin_a_max) :: delta_hno3_max
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    
    real(r8) :: small_aer, small_amt, small_gas

    small_gas = 0.01 * delta_hno3_max(ibin)
    small_aer = 0.01 * aer(icl_a,jliquid,ibin)

    small_amt = min(small_gas, small_aer)

    
    aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) + small_amt
    aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
         aer(ino3_a,jliquid,ibin)
    gas(ihno3_g) = ((gas(ihno3_g))-(small_amt))

    
    aer(icl_a,jliquid,ibin)  = ((aer(icl_a,jliquid,ibin))-   &
         (small_amt))
    aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin) +   &
         aer(icl_a,jliquid,ibin)

    
    gas(ihcl_g) = gas(ihcl_g) + small_amt

    return
  end subroutine absorb_tiny_hno3



  
  
  subroutine absorb_tiny_hcl(ibin,aer,gas,delta_hcl_max)  
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jliquid,jtotal,jsolid,                                      &
         ino3_a,icl_a,ihcl_g,ihno3_g
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(nbin_a_max) :: delta_hcl_max
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    
    real(r8) :: small_aer, small_amt, small_gas

    small_gas = 0.01 * delta_hcl_max(ibin)
    small_aer = 0.01 * aer(ino3_a,jliquid,ibin)

    small_amt = min(small_gas, small_aer)

    
    aer(icl_a,jliquid,ibin)= aer(icl_a,jliquid,ibin) + small_amt
    aer(icl_a,jtotal,ibin) = aer(icl_a,jsolid,ibin) +   &
         aer(icl_a,jliquid,ibin)
    gas(ihcl_g) = ((gas(ihcl_g))-(small_amt))

    
    aer(ino3_a,jliquid,ibin) = ((aer(ino3_a,jliquid,ibin))-   &
         (small_amt))
    aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
         aer(ino3_a,jliquid,ibin)

    
    gas(ihno3_g) = gas(ihno3_g) + small_amt

    return
  end subroutine absorb_tiny_hcl
  


  
  
  subroutine degas_tiny_nh4no3(ibin,aer,gas,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jliquid,jsolid,jtotal,                                      &
         jnh4no3,inh4_a,ino3_a,inh3_g,ihno3_g
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    real(r8) :: small_amt

    small_amt = 0.01 * electrolyte(jnh4no3,jliquid,ibin)

    aer(inh4_a,jliquid,ibin) = ((aer(inh4_a,jliquid,ibin))-   &
         (small_amt))
    aer(ino3_a,jliquid,ibin) = ((aer(ino3_a,jliquid,ibin))-   &
         (small_amt))

    
    aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
         aer(inh4_a,jliquid,ibin)
    aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
         aer(ino3_a,jliquid,ibin)

    
    gas(inh3_g)  = gas(inh3_g)  + small_amt
    gas(ihno3_g) = gas(ihno3_g) + small_amt

    return
  end subroutine degas_tiny_nh4no3




  
  
  subroutine degas_tiny_nh4cl(ibin,aer,gas,electrolyte)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jliquid,jsolid,jtotal,                                      &
         jnh4cl,inh4_a,icl_a,inh3_g,ihcl_g
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    real(r8) :: small_amt


    small_amt = 0.01 * electrolyte(jnh4cl,jliquid,ibin)

    aer(inh4_a,jliquid,ibin) = ((aer(inh4_a,jliquid,ibin))-   &
         (small_amt))
    aer(icl_a,jliquid,ibin)  = ((aer(icl_a,jliquid,ibin))-   &
         (small_amt))

    
    aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
         aer(inh4_a,jliquid,ibin)
    aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
         aer(icl_a,jliquid,ibin)

    
    gas(inh3_g) = gas(inh3_g) + small_amt
    gas(ihcl_g) = gas(ihcl_g) + small_amt

    return
  end subroutine degas_tiny_nh4cl



  
  
  
  
  
  
  subroutine equilibrate_acids(ibin,aer,gas,electrolyte,activity,mc,water_a,       &
       total_species,tot_cl_in,ma,gam,Keq_ll,Keq_gl)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,ngas_volatile,Ncation,Nanion,nrxn_aer_gl,nrxn_aer_ll,       &
         ihno3_g,ihcl_g
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(ngas_volatile) :: total_species
    real(r8), intent(inout) :: tot_cl_in
    real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    



    if(gas(ihcl_g)*gas(ihno3_g) .gt. 0.)then
       call equilibrate_hcl_and_hno3(ibin,aer,gas,electrolyte,activity,mc,water_a, &
            total_species,tot_cl_in,ma,gam,Keq_ll,Keq_gl)
    elseif(gas(ihcl_g) .gt. 0.)then
       call equilibrate_hcl(ibin,aer,gas,electrolyte,activity,mc,water_a,          &
            total_species,tot_cl_in,ma,gam,Keq_ll,Keq_gl)
    elseif(gas(ihno3_g) .gt. 0.)then
       call equilibrate_hno3(ibin,aer,gas,electrolyte,activity,mc,water_a,         &
            total_species,tot_cl_in,ma,gam,Keq_ll,Keq_gl)
    endif


    return
  end subroutine equilibrate_acids



  
  subroutine equilibrate_hcl(ibin,aer,gas,electrolyte,activity,mc,water_a,         &
       total_species,tot_cl_in,ma,gam,Keq_ll,Keq_gl)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,ngas_volatile,Ncation,jliquid,jsolid,jtotal,Nanion,         &
         nrxn_aer_gl,nrxn_aer_ll,                                                  &
         ja_so4,ja_hso4,ihcl_g,icl_a,jhcl,ino3_a,ica_a,inh4_a,ina_a,jc_h,jc_ca,    &
         jc_nh4,jc_na,ja_cl,ja_no3,jhno3,jnh4cl
    
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(ngas_volatile) :: total_species
    real(r8), intent(inout) :: tot_cl_in
    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nrxn_aer_gl) ::Keq_gl
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    real(r8) :: a, aerH, aerHSO4, aerSO4, b, c, dum, Kdash_hcl, mH, Tcl,   &
         W, XT, Z
    

    aerSO4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
    aerHSO4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

    Tcl = aer(icl_a,jliquid,ibin) + gas(ihcl_g)   
    Kdash_hcl = Keq_gl(4)*1.e+18/gam(jhcl,ibin)**2        
    Z = (   aer(ina_a, jliquid,ibin) +               &  
         aer(inh4_a,jliquid,ibin) +   &
         2.*aer(ica_a, jliquid,ibin) ) -   &
         (2.*aerSO4  +   &
         aerHSO4 +   &
         aer(ino3_a,jliquid,ibin) )


    W     = water_a(ibin)                         

    Kdash_hcl = Keq_gl(4)*1.e+18/gam(jhcl,ibin)**2        
    a = 1.0
    b = ((Kdash_hcl*W) + (Z/W))*1.e-9
    c = Kdash_hcl*(Z - Tcl)*1.e-18


    dum = ((b*b)-(4.*a*c))
    if (dum .lt. 0.) return               


    if(c .lt. 0.)then
       mH = quadratic(a,b,c)      
       aerH = mH*W*1.e+9
       aer(icl_a,jliquid,ibin) = ((aerH) + (Z))
    else
       mH = sqrt(Keq_ll(3))
    endif

    call form_electrolytes(jliquid,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)

    
    gas(ihcl_g) = ( (Tcl)  - (aer(icl_a,jliquid,ibin))  )


    
    ma(ja_so4,ibin)  = 1.e-9*aerSO4/water_a(ibin)
    ma(ja_hso4,ibin) = 1.e-9*aerHSO4/water_a(ibin)
    ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
    ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

    mc(jc_h,ibin)    = mH
    mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
    mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
    mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)


    
    activity(jhcl,ibin)    = mc(jc_h,ibin)  *ma(ja_cl,ibin)  *   &
         gam(jhcl,ibin)**2

    activity(jhno3,ibin)   = mc(jc_h,ibin)  *ma(ja_no3,ibin) *   &
         gam(jhno3,ibin)**2

    activity(jnh4cl,ibin)  = mc(jc_nh4,ibin)*ma(ja_cl,ibin) *   &
         gam(jnh4cl,ibin)**2


    
    aer(icl_a,jtotal,ibin) = aer(icl_a,jliquid,ibin) +   &
         aer(icl_a,jsolid,ibin)

    electrolyte(jhcl,jtotal,ibin) = electrolyte(jhcl,jliquid,ibin)

    return
  end subroutine equilibrate_hcl



  
  subroutine equilibrate_hno3(ibin,aer,gas,electrolyte,activity,mc,water_a,        &
       total_species,tot_cl_in,ma,gam,Keq_ll,Keq_gl)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,ngas_volatile,Ncation,jliquid,jsolid,jtotal,Nanion,         &
         nrxn_aer_gl,nrxn_aer_ll,                                                  &
         ja_so4,ja_hso4,ihno3_g,ino3_a,jhno3,icl_a,ica_a,inh4_a,ina_a,jc_h,jc_ca,  &
         jc_nh4,jc_na,ja_cl,jhcl,ja_no3,jnh4no3
    
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(ngas_volatile) :: total_species
    real(r8), intent(inout) :: tot_cl_in
    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    real(r8) :: a, aerH, aerHSO4, aerSO4, b, c, dum, Kdash_hno3, mH,   &
         Tno3, W, XT, Z
    

    aerSO4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
    aerHSO4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

    Tno3 = aer(ino3_a,jliquid,ibin) + gas(ihno3_g)        
    Kdash_hno3 = Keq_gl(3)*1.e+18/gam(jhno3,ibin)**2      
    Z = (   aer(ina_a, jliquid,ibin) +               &  
         aer(inh4_a,jliquid,ibin) +   &
         2.*aer(ica_a, jliquid,ibin) ) -   &
         (2.*aerSO4  +   &
         aerHSO4 +   &
         aer(icl_a,jliquid,ibin) )


    W     = water_a(ibin)                         

    Kdash_hno3 = Keq_gl(3)*1.e+18/gam(jhno3,ibin)**2      
    a = 1.0
    b = ((Kdash_hno3*W) + (Z/W))*1.e-9
    c = Kdash_hno3*(Z - Tno3)*1.e-18

    dum = ((b*b)-(4.*a*c))
    if (dum .lt. 0.) return               



    if(c .lt. 0.)then
       mH = quadratic(a,b,c)      
       aerH = mH*W*1.e+9
       aer(ino3_a,jliquid,ibin) = ((aerH) + (Z))
    else
       mH = sqrt(Keq_ll(3))
    endif

    call form_electrolytes(jliquid,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)

    
    gas(ihno3_g)= ( (Tno3) - (aer(ino3_a,jliquid,ibin)) )


    
    ma(ja_so4,ibin)  = 1.e-9*aerSO4/water_a(ibin)
    ma(ja_hso4,ibin) = 1.e-9*aerHSO4/water_a(ibin)
    ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
    ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

    mc(jc_h,ibin)    = mH
    mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
    mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
    mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)


    
    activity(jhcl,ibin)    = mc(jc_h,ibin)  *ma(ja_cl,ibin)  *   &
         gam(jhcl,ibin)**2

    activity(jhno3,ibin)   = mc(jc_h,ibin)  *ma(ja_no3,ibin) *   &
         gam(jhno3,ibin)**2

    activity(jnh4no3,ibin) = mc(jc_nh4,ibin)*ma(ja_no3,ibin) *   &
         gam(jnh4no3,ibin)**2


    
    aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
         aer(ino3_a,jsolid,ibin)

    electrolyte(jhno3,jtotal,ibin) = electrolyte(jhno3,jliquid,ibin)

    return
  end subroutine equilibrate_hno3



  
  subroutine equilibrate_hcl_and_hno3(ibin,aer,gas,electrolyte,activity,mc,        &
       water_a,total_species,tot_cl_in,ma,gam,Keq_ll,Keq_gl)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,ngas_volatile,Ncation,jliquid,jsolid,jtotal,Nanion,         &
         nrxn_aer_gl,nrxn_aer_ll,                                                  &
         ja_so4,ja_hso4,ihcl_g,icl_a,ihno3_g,ino3_a,jhcl,jhno3,             &
         ica_a,inh4_a,ina_a,jc_h,jc_ca,jc_nh4,jc_na,ja_cl,ja_no3,jnh4no3,   &
         jnh4cl
    
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(ngas_volatile) :: total_species
    real(r8), intent(inout) :: tot_cl_in
    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    real(r8) :: aerH, aerHSO4, aerSO4, Kdash_hcl, Kdash_hno3,   &
         mH, p, q, r, Tcl, Tno3, W, XT, Z
    


    aerSO4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
    aerHSO4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

    Tcl  = aer(icl_a,jliquid,ibin)  + gas(ihcl_g) 
    Tno3 = aer(ino3_a,jliquid,ibin) + gas(ihno3_g)        

    Kdash_hcl  = Keq_gl(4)*1.e+18/gam(jhcl,ibin)**2       
    Kdash_hno3 = Keq_gl(3)*1.e+18/gam(jhno3,ibin)**2      

    Z = (   aer(ina_a, jliquid,ibin) +               &  
         aer(inh4_a,jliquid,ibin) +   &
         2.*aer(ica_a, jliquid,ibin) ) -   &
         (2.*aerSO4 + aerHSO4 )


    W = water_a(ibin)

    Kdash_hcl  = Keq_gl(4)*1.e+18/gam(jhcl,ibin)**2       
    Kdash_hno3 = Keq_gl(3)*1.e+18/gam(jhno3,ibin)**2      

    p = (Z/W + W*(Kdash_hcl + Kdash_hno3))*1.e-9

    q = 1.e-18*Kdash_hcl*Kdash_hno3*W**2  +   &
         1.e-18*Z*(Kdash_hcl + Kdash_hno3) -   &
         1.e-18*Kdash_hcl*Tcl -   &
         1.e-18*Kdash_hno3*Tno3

    r = 1.e-18*Kdash_hcl*Kdash_hno3*W*(Z - Tcl - Tno3)*1.e-9

    mH = cubic(p,q,r)

    if(mH .gt. 0.0)then
       aerH = mH*W*1.e+9
       aer(ino3_a,jliquid,ibin) = Kdash_hno3*W*W*Tno3/   &
            (aerH + Kdash_hno3*W*W)
       aer(icl_a, jliquid,ibin) = Kdash_hcl*W*W*Tcl/   &
            (aerH + Kdash_hcl*W*W)
    else
       mH = sqrt(Keq_ll(3))
    endif

    call form_electrolytes(jliquid,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)

    
    gas(ihno3_g)= ( (Tno3) - (aer(ino3_a,jliquid,ibin)) )
    gas(ihcl_g) = ( (Tcl)  - (aer(icl_a,jliquid,ibin))  )


    
    ma(ja_so4,ibin)  = 1.e-9*aerSO4/water_a(ibin)
    ma(ja_hso4,ibin) = 1.e-9*aerHSO4/water_a(ibin)
    ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
    ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

    mc(jc_h,ibin)    = mH
    mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
    mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
    mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)


    
    activity(jhcl,ibin)    = mc(jc_h,ibin)*ma(ja_cl,ibin)   *   &
         gam(jhcl,ibin)**2

    activity(jhno3,ibin)   = mc(jc_h,ibin)*ma(ja_no3,ibin)  *   &
         gam(jhno3,ibin)**2

    activity(jnh4no3,ibin) = mc(jc_nh4,ibin)*ma(ja_no3,ibin)*   &
         gam(jnh4no3,ibin)**2

    activity(jnh4cl,ibin)  = mc(jc_nh4,ibin)*ma(ja_cl,ibin) *   &
         gam(jnh4cl,ibin)**2


    
    aer(icl_a,jtotal,ibin)  = aer(icl_a,jliquid,ibin) +   &
         aer(icl_a,jsolid,ibin)

    aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
         aer(ino3_a,jsolid,ibin)

    electrolyte(jhno3,jtotal,ibin) = electrolyte(jhno3,jliquid,ibin)
    electrolyte(jhcl, jtotal,ibin) = electrolyte(jhcl, jliquid,ibin)

    return
  end subroutine equilibrate_hcl_and_hno3



  
  
  
  
  
  
  
  
  subroutine degas_solid_nh4no3(ibin,aer,gas,electrolyte,Keq_sg)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jsolid,jliquid,jtotal,nrxn_aer_sg,                          &
         ihno3_g,inh3_g,jnh4no3,inh4_a,ino3_a

    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer jp
    real(r8) :: a, b, c, xgas, XT
    


    jp = jsolid

    a = 1.0
    b = gas(inh3_g) + gas(ihno3_g)
    c = gas(inh3_g)*gas(ihno3_g) - Keq_sg(1)
    xgas = quadratic(a,b,c)

    if(xgas .ge. electrolyte(jnh4no3,jp,ibin))then 

       gas(inh3_g) = gas(inh3_g)  + electrolyte(jnh4no3,jp,ibin)
       gas(ihno3_g)= gas(ihno3_g) + electrolyte(jnh4no3,jp,ibin)
       aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) -   &
            electrolyte(jnh4no3,jp,ibin)
       aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) -   &
            electrolyte(jnh4no3,jp,ibin)

    else  

       gas(inh3_g) = gas(inh3_g)  + xgas
       gas(ihno3_g)= gas(ihno3_g) + xgas
       aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) - xgas
       aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) - xgas
    endif


    
    aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
         aer(inh4_a,jliquid,ibin)
    aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
         aer(ino3_a,jliquid,ibin)

    return
  end subroutine degas_solid_nh4no3



  
  subroutine degas_solid_nh4cl(ibin,aer,gas,electrolyte,Keq_sg)
    use module_data_mosaic_aero, only: r8,naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jsolid,jliquid,jtotal,nrxn_aer_sg,                          &
         ihcl_g,inh3_g,jnh4cl,inh4_a,icl_a
    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer jp
    real(r8) :: a, b, c, xgas, XT
    


    jp = jsolid

    a = 1.0
    b = gas(inh3_g) + gas(ihcl_g)
    c = gas(inh3_g)*gas(ihcl_g) - Keq_sg(2)
    xgas = quadratic(a,b,c)

    if(xgas .ge. electrolyte(jnh4cl,jp,ibin))then 

       gas(inh3_g) = gas(inh3_g) + electrolyte(jnh4cl,jp,ibin)
       gas(ihcl_g) = gas(ihcl_g) + electrolyte(jnh4cl,jp,ibin)
       aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) -   &
            electrolyte(jnh4cl,jp,ibin)
       aer(icl_a,jp,ibin)  = aer(icl_a,jp,ibin) -   &
            electrolyte(jnh4cl,jp,ibin)

    else  

       gas(inh3_g) = gas(inh3_g) + xgas
       gas(ihcl_g) = gas(ihcl_g) + xgas
       aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) - xgas
       aer(icl_a,jp,ibin)  = aer(icl_a,jp,ibin)  - xgas

    endif


    
    aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
         aer(inh4_a,jliquid,ibin)
    aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
         aer(icl_a,jliquid,ibin)

    return
  end subroutine degas_solid_nh4cl



  
  
  
  
  
  
  subroutine conform_electrolytes(jp,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)

    use module_data_mosaic_aero, only: r8,naer,nbin_a_max,           &
         ngas_aerchtot, ngas_volatile, nelectrolyte,                 &
         imsa_a,iso4_a,ica_a,ina_a,inh4_a,ino3_a,icl_a,ico3_a

    implicit none

    
    integer, intent(in) :: ibin, jp
    real(r8), intent(inout) :: XT
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(ngas_volatile) :: total_species
    real(r8), intent(inout) :: tot_cl_in
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    integer i, iXT_case, je
    real(r8) :: sum_dum, XNa_prime, XNH4_prime, XT_prime
    real(r8) :: store(naer)

    
    
    
    


    

    if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
       XT   = ( aer(inh4_a,jp,ibin) +   &
            aer(ina_a,jp,ibin)  +   &
            2.*aer(ica_a,jp,ibin) )/   &
            (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
    else
       XT   = -1.0
    endif



    if(XT .ge. 2.0 .or. XT.lt.0.)then 
       iXT_case = 1       
    else
       iXT_case = 2       
    endif

    
    
    
    store(iso4_a) = aer(iso4_a,jp,ibin)
    store(ino3_a) = aer(ino3_a,jp,ibin)
    store(icl_a)  = aer(icl_a, jp,ibin)
    store(imsa_a) = aer(imsa_a,jp,ibin)
    store(ico3_a) = aer(ico3_a,jp,ibin)
    store(inh4_a) = aer(inh4_a,jp,ibin)
    store(ina_a)  = aer(ina_a, jp,ibin)
    store(ica_a)  = aer(ica_a, jp,ibin)

    do je=1,nelectrolyte
       electrolyte(je,jp,ibin) = 0.0
    enddo

    
    
    
    if(iXT_case.eq.1)then

       

       call form_caso4(store,jp,ibin,electrolyte)
       call form_camsa2(store,jp,ibin,electrolyte)
       call form_na2so4(store,jp,ibin,electrolyte)
       call form_namsa(store,jp,ibin,electrolyte)
       call form_cano3(store,jp,ibin,electrolyte)
       call form_nano3(store,jp,ibin,electrolyte)
       call form_nacl(store,jp,ibin,aer,gas,electrolyte,total_species,tot_cl_in)
       call form_cacl2(store,jp,ibin,electrolyte)
       call form_caco3(store,jp,ibin,aer,electrolyte)
       call form_nh4so4(store,jp,ibin,electrolyte)
       call form_nh4msa(store,jp,ibin,electrolyte)
       call form_nh4no3(store,jp,ibin,electrolyte)
       call form_nh4cl(store,jp,ibin,electrolyte)
       call form_msa(store,jp,ibin,electrolyte)
       call degas_hno3(store,jp,ibin,aer,gas,electrolyte)
       call degas_hcl(store,jp,ibin,aer,gas,electrolyte)
       call degas_nh3(store,jp,ibin,aer,gas)

    elseif(iXT_case.eq.2)then

       

       call form_caso4(store,jp,ibin,electrolyte)
       call form_camsa2(store,jp,ibin,electrolyte)
       call form_namsa(store,jp,ibin,electrolyte)
       call form_nh4msa(store,jp,ibin,electrolyte)
       call form_msa(store,jp,ibin,electrolyte)

       if(store(iso4_a).eq.0.0)goto 10


       XT_prime =(store(ina_a)+store(inh4_a))/   &
            store(iso4_a)
       XNa_prime=0.5*store(ina_a)/store(iso4_a) + 1.

       if(XT_prime.ge.XNa_prime)then
          call form_na2so4(store,jp,ibin,electrolyte)
          XNH4_prime = 0.0
          if(store(iso4_a).gt.1.e-15)then
             XNH4_prime = store(inh4_a)/store(iso4_a)
          endif

          if(XNH4_prime .ge. 1.5)then
             call form_nh4so4_lvcite(store,jp,ibin,electrolyte)
          else
             call form_lvcite_nh4hso4(store,jp,ibin,electrolyte)
          endif

       elseif(XT_prime.ge.1.)then
          call form_nh4hso4(store,jp,ibin,electrolyte)
          call form_na2so4_nahso4(store,jp,ibin,electrolyte)
       elseif(XT_prime.lt.1.)then
          call form_nahso4(store,jp,ibin,electrolyte)
          call form_nh4hso4(store,jp,ibin,electrolyte)
          call form_h2so4(store,jp,ibin,electrolyte)
       endif

10     call degas_hno3(store,jp,ibin,aer,gas,electrolyte)
       call degas_hcl(store,jp,ibin,aer,gas,electrolyte)
       call degas_nh3(store,jp,ibin,aer,gas)

    endif 


    
    call electrolytes_to_ions(jp, ibin,aer,electrolyte)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    return
  end subroutine conform_electrolytes



  
  
  
  function cubic( psngl, qsngl, rsngl )
    use module_data_mosaic_kind, only:  r8
    implicit none
    real(r8) :: cubic
    
    real(r8) :: psngl, qsngl, rsngl
    
    real(r8) :: p, q, r, A, B, D, M, N, third, y
    real(r8) :: k, phi, thesign, x(3), duma
    integer icase, kk

    third = 1.d0/3.d0

    q = (qsngl)
    p = (psngl)
    r = (rsngl)

    A = (1.d0/3.d0)*((3.d0*q) - (p*p))
    B = (1.d0/27.d0)*((2.d0*p*p*p) - (9.d0*p*q) + (27.d0*r))

    D = ( ((A*A*A)/27.d0) + ((B*B)/4.d0) )

    if(D .gt. 0.)then     
       icase = 1
    elseif(D .eq. 0.)then 
       icase = 2
    else  
       icase = 3
    endif


    goto (1,2,3), icase

    
1   thesign = 1.
    if(B .gt. 0.)then
       B = -B
       thesign = -1.
    endif

    M = thesign*((-B/2.d0) + (sqrt(D)))**(third)
    N = thesign*((-B/2.d0) - (sqrt(D)))**(third)

    cubic = ( (M) + (N) - (p/3.d0) )
    return

    
2   thesign = 1.
    if(B .gt. 0.)then
       B = -B
       thesign = -1.
    endif

    M = thesign*(-B/2.d0)**third
    N = M

    x(1) = ( (M) + (N) - (p/3.d0) )
    x(2) = ( (-M/2.d0) + (-N/2.d0) - (p/3.d0) )
    x(2) = ( (-M/2.d0) + (-N/2.d0) - (p/3.d0) )

    cubic = 0.
    do kk = 1, 3
       if(x(kk).gt.cubic) cubic = x(kk)
    enddo
    return

    
3   if(B.gt.0.)then
       thesign = -1.
    elseif(B.lt.0.)then
       thesign = 1.
    endif

    
    
    duma = thesign*sqrt( (B*B/4.d0)/(-A*A*A/27.d0) )
    duma = min( duma, +1.0d0 )
    duma = max( duma, -1.0d0 )
    phi  = acos( duma )   


    cubic = 0.
    do kk = 1, 3
       k = kk-1
       y = 2.*Sqrt(-A/3.)*cos(phi + 120.*k*0.017453293)
       x(kk) = ((y) - (p/3.d0))
       if(x(kk).gt.cubic) cubic = x(kk)
    enddo
    return

  end function cubic
   


  
  function quadratic(a,b,c)
    use module_data_mosaic_kind, only:  r8
    implicit none
    real(r8) :: quadratic
    
    real(r8) :: a, b, c
    
    real(r8) :: x, dum, quad1, quad2


    if(b .ne. 0.0)then
       x = 4.*(a/b)*(c/b)
    else
       x = 1.e+6
    endif

    if(abs(x) .lt. 1.e-6)then
       dum = ( (0.5*x) +   &
            (0.125*x**2) +   &
            (0.0625*x**3) )

       quadratic = (-0.5*b/a)*dum

       if(quadratic .lt. 0.)then
          quadratic = -b/a - quadratic
       endif

    else
       quad1 = ((-b)+sqrt((b*b)-(4.*a*c)))/   &
            (2.*a)
       quad2 = ((-b)-sqrt((b*b)-(4.*a*c)))/   &
            (2.*a)

       quadratic = max(quad1, quad2)
    endif

    return
  end function quadratic
  


  
  function mean_molecular_speed(T, MW)    
    use module_data_mosaic_kind, only:  r8
    implicit none
    real(r8) :: mean_molecular_speed
    
    real(r8) :: T, MW     

    mean_molecular_speed = 1.455e4 * sqrt(T/MW)

    return
  end function mean_molecular_speed
  

  
  function gas_diffusivity(T, P, MW, Vm)  
    use module_data_mosaic_kind, only:  r8
    use module_data_mosaic_constants, only:  third
    implicit none
    real(r8) :: gas_diffusivity
    
    real(r8) :: MW, Vm, T, P      


    gas_diffusivity = (1.0e-3 * T**1.75 * sqrt(1./MW + 0.035))/   &
         (P * (Vm**third + 2.7189)**2)


    return
  end function gas_diffusivity
  


  
  function fuchs_sutugin(rkn,a)
    use module_data_mosaic_kind, only:  r8
    implicit none
    real(r8) :: fuchs_sutugin
    
    real(r8) :: rkn, a
    
    real(r8) :: rnum, denom


    rnum  = 0.75*a*(1. + rkn)
    denom = rkn**2 + rkn + 0.283*rkn*a + 0.75*a
    fuchs_sutugin = rnum/denom

    return
  end function fuchs_sutugin
  


  
  
  
  function aerosol_water_up( ibin, electrolyte, aer, kappa_nonelectro, a_zsr ) 

    use module_data_mosaic_aero, only: r8, nelectrolyte, naer, nbin_a_max, jtotal, &
        nsalt, ioc_a, ibc_a, ilim2_a, ioin_a, dens_aer_mac, mw_aer_mac 

    use module_data_mosaic_asecthp, only: dens_water_aer

    implicit none

    real(r8) :: aerosol_water_up
    
    integer, intent(in) :: ibin
    real(r8), intent(in), dimension (6,nelectrolyte) :: a_zsr
    real(r8), intent(in), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), intent(in), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(in), dimension(naer) :: kappa_nonelectro

    
    integer :: iaer, jp, je
    real(r8) :: tmpa, tmpb, aH2O_60  
    
    


    aH2O_60 = 0.6

    jp = jtotal
    tmpa = 0.0_r8

    do je = 1, (nsalt+4)  
       tmpa = tmpa + electrolyte(je,jp,ibin)/bin_molality_60(je,a_zsr)
    enddo







    tmpb = 0.0_r8
    do iaer = 1, naer
       if (kappa_nonelectro(iaer) > 0.0_r8) then
          tmpb = tmpb + (aer(iaer,jp,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer))*kappa_nonelectro(iaer)
       end if
    end do
    tmpa = tmpa + tmpb * dens_water_aer * 1.0e-3 * aH2O_60/(1.0-aH2O_60)   

    aerosol_water_up = tmpa*1.e-9  

    return
  end function aerosol_water_up
  


  
  function bin_molality_60(je,a_zsr)            
    use module_data_mosaic_aero, only: r8,nelectrolyte

    implicit none

    real(r8) :: bin_molality_60
    
    integer, intent(in) ::  je
    real(r8), intent(in), dimension (6,nelectrolyte) :: a_zsr
    
    real(r8) :: aw, xm


    aw = 0.6_r8

    xm =     a_zsr(1,je) +   &
         aw*(a_zsr(2,je) +   &
         aw*(a_zsr(3,je) +   &
         aw*(a_zsr(4,je) +   &
         aw*(a_zsr(5,je) +   &
         aw* a_zsr(6,je) ))))

    bin_molality_60 = 55.509_r8*xm/(1. - xm)

    return
  end function bin_molality_60
  

  
  
  
  function aerosol_water( jp, ibin, jaerosolstate, jphase, jhyst_leg, electrolyte, aer,   &
           kappa_nonelectro, num_a, mass_dry_a, mass_soluble_a, aH2O, molality0 ) 

    use module_data_mosaic_aero, only: r8, nbin_a_max, nelectrolyte, nsoluble, naer,   &
         all_solid, jsolid, jhyst_lo, ioc_a, ibc_a, ilim2_a, ioin_a,   &   
         jtotal, ah2o_max, dens_aer_mac, mw_aer_mac, ename

    use module_data_mosaic_asecthp, only: dens_water_aer

    implicit none

    real(r8) :: aerosol_water
    
    integer, intent(in) :: jp, ibin
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase,jhyst_leg

    real(r8), intent(in) :: aH2O
    real(r8), intent(in), dimension(nbin_a_max) :: num_a,mass_dry_a,mass_soluble_a
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(in), dimension(naer) :: kappa_nonelectro

    
    integer :: iaer, iclm_aer, jclm_aer, je
    real(r8) :: tmpa, tmpb
    
    real(r8) :: bin_molality



    tmpa = 0.0_r8
    if (jaerosolstate(ibin) .ne. all_solid) then                     


    do je = 1, 19   
       tmpa = tmpa + electrolyte(je,jp,ibin)/molality0(je,ibin)      
    enddo
    endif










    tmpb = 0.0_r8
    do iaer = 1, naer
       if (kappa_nonelectro(iaer) > 0.0_r8) then
          tmpb = tmpb + (aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer))*kappa_nonelectro(iaer)
       end if
    end do
    tmpa = tmpa + tmpb * dens_water_aer * 1.0e-3 * ah2o/(1.0-min(ah2o,ah2o_max))  

    aerosol_water = tmpa*1.e-9  

                 
    if (aerosol_water .le. 0.0) then 
       iclm_aer = 0 
       jclm_aer = 0 

       





       
       
       





       

       jaerosolstate(ibin) = all_solid
       jphase(ibin)    = jsolid
       jhyst_leg(ibin) = jhyst_lo

    endif

44  format(a7, 2x, e11.3)


    return
  end function aerosol_water





  
  function bin_molality(je,ibin,aH2O_a,b_zsr,a_zsr,aw_min)
    use module_data_mosaic_aero, only:r8,  nbin_a_max, nelectrolyte

    implicit none

    real(r8) :: bin_molality
    
    integer, intent(in) :: je, ibin
    real(r8), intent(in), dimension(nbin_a_max) :: aH2O_a
    real(r8), intent(in), dimension(nelectrolyte) :: b_zsr,aw_min
    real(r8), intent(in), dimension (6,nelectrolyte) :: a_zsr
    
    real(r8) :: aw, xm


    aw = max(aH2O_a(ibin), aw_min(je))
    aw = min(aw, 0.999999_r8)


    if(aw .lt. 0.97_r8)then

       xm =     a_zsr(1,je) +   &
            aw*(a_zsr(2,je) +   &
            aw*(a_zsr(3,je) +   &
            aw*(a_zsr(4,je) +   &
            aw*(a_zsr(5,je) +   &
            aw* a_zsr(6,je) ))))

       bin_molality = 55.509_r8*xm/(1. - xm)

    else

       bin_molality = -b_zsr(je)*log(aw)

    endif


    return
  end function bin_molality
  







end module module_mosaic_ext
