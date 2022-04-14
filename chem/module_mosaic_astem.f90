  module module_mosaic_astem


  use module_mosaic_support, only: mosaic_warn_mess, mosaic_err_mess
  use module_data_mosaic_kind, only: r8

  implicit none


  contains










  
  subroutine ASTEM(   mcall_print_aer,                                               &
       dtchem,        sigmag_a,  aH2O,   T_K,          RH_pc,     P_atm,             &
       kappa_nonelectro,                                                             &
       jaerosolstate, flux_s,            flux_l,       volatile_s,iprint_input,      &
       phi_volatile_s,phi_volatile_l,    jphase,       aer,       kg,       gas,     &
       gas_avg,       gas_netprod_otrproc,                                           &
       jhyst_leg,     electrolyte,       epercent,     activity,  mc,       sat_soa, &
       num_a,         Dp_dry_a,          Dp_wet_a,     dp_core_a, mass_dry_a,        &
       mass_soluble_a,vol_dry_a,         dens_dry_a,   water_a,   water_a_hyst,      &
       water_a_up,    aH2O_a,            total_species,tot_cl_in, ma,       gam,     &
       log_gamZ,      gam_ratio,         Keq_ll,       Keq_gl,    Keq_sg,   Kp_nh4cl,&
       Kp_nh4no3,     sigma_water,       Keq_sl,       MDRH_T,    molality0,         &
       uptkrate_h2so4,                   mosaic_vars_aa,                             &
       area_dry_a,    area_wet_a,        mass_wet_a,vol_wet_a,                       &
       dens_wet_a,    ri_shell_a,        ri_avg_a,     ri_core_a                     )
  
  use module_data_mosaic_aero, only: r8, nbin_a_max,                               &
       ngas_aerchtot, ngas_volatile, nelectrolyte,                                 &
       Ncation, naer, mYES, no_aerosol, Nanion, nrxn_aer_gl, nrxn_aer_ll,          &
       nrxn_aer_sg, nrxn_aer_sl, naercomp, nsalt, MDRH_T_NUM, jsulf_poor_NUM,      &
       jsulf_rich_NUM, nbin_a, zc, za,                                             &
       a_zsr, b_zsr, mw_electrolyte, partial_molar_vol, mw_aer_mac, dens_aer_mac,  &
       MW_a, MW_c, dens_comp_a, mw_comp_a, ref_index_a, rtol_mesa, jsalt_index,    &
       jsulf_poor, jsulf_rich, ih2so4_g,                                           &
       iso4_a, jtotal, msoa_flag1,                                                 & 
       mosaic_vars_aa_type
  
  use module_mosaic_ext, only: aerosol_phase_state,calc_dry_n_wet_aerosol_props,   &
       aerosolmtc, dumpxx
  
  use module_mosaic_soa_vbs, only: mosaic_soa_vbs_intr
  


  
  
  
  
  integer, intent(in) :: mcall_print_aer

  real(r8), intent(in) :: dtchem
  real(r8), intent(in) :: aH2O
  real(r8), intent(in) :: T_K, RH_pc, P_atm
  real(r8), intent(in), dimension(nbin_a_max) :: sigmag_a
  real(r8), intent(in), dimension(ngas_aerchtot) :: gas_netprod_otrproc
  real(r8), intent(in), dimension(naer) :: kappa_nonelectro
  
  
  integer, intent(inout) :: iprint_input

  integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate, jphase
  integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg

  real(r8), intent(inout) :: tot_cl_in
  real(r8), intent(inout) :: Kp_nh4cl, Kp_nh4no3
  real(r8), intent(inout) :: sigma_water
     
  real(r8), intent(inout), dimension(nbin_a_max)     :: num_a, Dp_dry_a, Dp_wet_a
  real(r8), intent(inout), dimension(nbin_a_max)     :: dp_core_a
  real(r8), intent(inout), dimension(nbin_a_max)     :: mass_dry_a, mass_soluble_a
  real(r8), intent(inout), dimension(nbin_a_max)     :: vol_dry_a, dens_dry_a
  real(r8), intent(inout), dimension(nbin_a_max)     :: water_a
  real(r8), intent(inout), dimension(nbin_a_max)     :: water_a_hyst,water_a_up
  real(r8), intent(inout), dimension(nbin_a_max)     :: aH2O_a
  real(r8), intent(inout), dimension(nbin_a_max)     :: gam_ratio
  real(r8), intent(inout), dimension(ngas_aerchtot)  :: gas
  real(r8), intent(inout), dimension(ngas_aerchtot)  :: gas_avg  
  real(r8), intent(inout)                            :: uptkrate_h2so4  

  type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

  
  
  
  
  
  real(r8), intent(inout), dimension(ngas_volatile)  :: sat_soa
  real(r8), intent(inout), dimension(ngas_volatile)  :: total_species
  real(r8), intent(inout), dimension(nrxn_aer_ll)    :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl)    :: Keq_gl
  real(r8), intent(inout), dimension(nrxn_aer_sg)    :: Keq_sg
  real(r8), intent(inout), dimension(nrxn_aer_sl)    :: Keq_sl
  real(r8), intent(inout), dimension(MDRH_T_NUM)     :: MDRH_T
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max)   :: molality0 

  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max)  :: flux_s,flux_l
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max)  :: volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max)  :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max)  :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max)  :: kg
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max)   :: activity
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max)   :: gam
  real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
  real(r8), intent(inout), dimension(Ncation,nbin_a_max)        :: mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max)         :: ma

  real(r8), intent(inout), dimension(naer,3,nbin_a_max)         :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent

  
  
  real(r8), intent(out), dimension(nbin_a_max) :: area_dry_a
  real(r8), intent(out), dimension(nbin_a_max) :: area_wet_a,mass_wet_a  
  real(r8), intent(out), dimension(nbin_a_max) :: vol_wet_a
  real(r8), intent(out), dimension(nbin_a_max) :: dens_wet_a

  complex,  intent(out), dimension(nbin_a_max) :: ri_shell_a,ri_avg_a,ri_core_a

  
  integer :: ibin, iv, itmpa

  integer, dimension(nsalt) :: jsalt_present
  integer, dimension(ngas_volatile,3,nbin_a_max) :: integrate


  real(r8) :: Keq_nh4cl
  real(r8) :: swdown_cell

  real(r8), dimension(nsalt) :: phi_salt_old
  real(r8), dimension(Ncation) :: nc_Mc,xeq_c
  real(r8), dimension(Nanion)  :: xeq_a,na_Ma    
  real(r8), dimension(nbin_a_max) :: sigma_soln
  real(r8), dimension(nbin_a_max) :: delta_nh3_max,delta_hno3_max
  real(r8), dimension(nbin_a_max) :: delta_hcl_max
  real(r8), dimension(nbin_a_max) :: growth_factor
  real(r8), dimension(nbin_a_max) :: MDRH

  real(r8), dimension(ngas_volatile) ::sfc_a
  real(r8), dimension(ngas_volatile,nbin_a_max) :: Heff
  real(r8), dimension(ngas_aerchtot,nbin_a_max) :: kel


  call dumpxx( 'cc', dtchem, t_k, p_atm, ah2o, &
         jaerosolstate, jphase, jhyst_leg, &
         aer, gas, num_a, water_a, water_a_hyst, dp_dry_a, &
         mosaic_vars_aa )

  mosaic_vars_aa%niter_MESA     = 0.0_r8
  mosaic_vars_aa%niter_MESA_max = 0
  mosaic_vars_aa%jMESA_fail     = 0
  mosaic_vars_aa%jMESA_call     = 0
  mosaic_vars_aa%iter_MESA(:)   = 0

  phi_salt_old(:)  = 0.0_r8
  integrate(:,:,:) = 0.0_r8 
  heff(:,:)        = 0.0_r8 
  
  gas_avg(:) = gas(:)  

  
  mosaic_vars_aa%jASTEM_call  = mosaic_vars_aa%jASTEM_call + 1
  
  
  iprint_input = mYES

  
  do ibin = 1, nbin_a
      area_dry_a(ibin) = 0.0_r8 
      area_wet_a(ibin) = 0.0_r8 
      mass_wet_a(ibin) = 0.0_r8 
     if(jaerosolstate(ibin) .ne. no_aerosol)then
        call aerosol_phase_state( ibin, jaerosolstate, jphase,  &
             aer, jhyst_leg, electrolyte, epercent, kel, activity, mc, num_a, mass_wet_a, &
             mass_dry_a, mass_soluble_a, vol_dry_a, vol_wet_a, water_a, water_a_hyst,  &
             water_a_up, aH2O_a, aH2O, ma, gam, & 
             log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c,              & 
             mw_electrolyte, partial_molar_vol, sigma_soln, T_K, RH_pc, mw_aer_mac,    &
             dens_aer_mac, sigma_water, Keq_ll, Keq_sl, MW_a, MW_c, growth_factor, MDRH, &
             MDRH_T, molality0, rtol_mesa, jsalt_present, jsalt_index, jsulf_poor,     &
             jsulf_rich, phi_salt_old,                                      &
             kappa_nonelectro, mosaic_vars_aa )

        call calc_dry_n_wet_aerosol_props(                          &
           ibin, jaerosolstate, aer, electrolyte, water_a, num_a,         &  
           dens_comp_a, mw_comp_a, dens_aer_mac, mw_aer_mac, ref_index_a, &  
           Dp_dry_a, Dp_wet_a, dp_core_a,                                 &  
           area_dry_a, area_wet_a, mass_dry_a, mass_wet_a,                &  
           vol_dry_a, vol_wet_a, dens_dry_a, dens_wet_a,                  &  
           ri_shell_a, ri_core_a, ri_avg_a                                )  
     endif
  enddo
  call check_astem_negative( 1, mosaic_vars_aa%xnerr_astem_negative, &
                                mosaic_vars_aa%fix_astem_negative, aer, gas )

  
  if (mcall_print_aer == 2) then
     
     
     
  endif			
  

  
  call aerosolmtc( jaerosolstate, num_a, Dp_wet_a, sigmag_a, P_atm, T_K, kg )

  uptkrate_h2so4 = sum( kg(ih2so4_g,1:nbin_a) )

  call dumpxx( 'ee', dtchem, t_k, p_atm, ah2o, &
         jaerosolstate, jphase, jhyst_leg, &
         aer, gas, num_a, water_a, water_a_hyst, dp_dry_a, &
         mosaic_vars_aa )

  
  call ASTEM_non_volatiles( dtchem, jaerosolstate, jphase, aer,  &
       kg, gas, gas_avg, gas_netprod_otrproc,                                        &
       jhyst_leg, electrolyte, epercent, kel, activity, mc, delta_nh3_max,              &
       delta_hno3_max, delta_hcl_max, num_a, mass_wet_a, mass_dry_a, mass_soluble_a,   &
       vol_dry_a, vol_wet_a, water_a, water_a_hyst, water_a_up, aH2O_a, total_species,  &
       tot_cl_in,                                                                 &
       aH2O, ma, gam, log_gamZ, zc, za,          &
       gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c, a_zsr, mw_electrolyte, partial_molar_vol,  &
       sigma_soln, T_K, RH_pc, mw_aer_mac, dens_aer_mac, sigma_water, Keq_ll, Keq_sl,    &
       MW_a, MW_c, growth_factor, MDRH, MDRH_T, molality0, rtol_mesa, jsalt_present,     &
       jsalt_index, jsulf_poor, jsulf_rich, phi_salt_old,                     &
       kappa_nonelectro, mosaic_vars_aa )	

  call dumpxx( 'gg', dtchem, t_k, p_atm, ah2o, &
         jaerosolstate, jphase, jhyst_leg, &
         aer, gas, num_a, water_a, water_a_hyst, dp_dry_a, &
         mosaic_vars_aa )

  call check_astem_negative( 2, mosaic_vars_aa%xnerr_astem_negative, &
                                mosaic_vars_aa%fix_astem_negative, aer, gas )

  
  call ASTEM_semi_volatiles( iprint_input, dtchem, jaerosolstate,                     &
       sfc_a, flux_s, flux_l, Heff, volatile_s, phi_volatile_s,   &
       jphase, aer, kg, gas, jhyst_leg, electrolyte, epercent, kel, activity, mc,          &
       delta_nh3_max, delta_hno3_max, delta_hcl_max,                                 &
       num_a, mass_dry_a, mass_wet_a, mass_soluble_a,   &
       vol_dry_a, vol_wet_a, water_a, total_species, tot_cl_in,                       &
       aH2O_a, aH2O, ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c, mw_electrolyte,    &
       Keq_ll, Keq_gl, Keq_sg, Kp_nh4cl, Kp_nh4no3, Keq_nh4cl, MW_c, MW_a, mw_aer_mac,    &
       dens_aer_mac, Keq_sl, growth_factor, MDRH, MDRH_T, molality0, rtol_mesa,         &
       jsalt_present, jsalt_index, jsulf_poor, jsulf_rich, phi_salt_old,    &
       integrate, phi_volatile_l,                                                      &
       kappa_nonelectro, mosaic_vars_aa )	

  call dumpxx( 'ii', dtchem, t_k, p_atm, ah2o, &
         jaerosolstate, jphase, jhyst_leg, &
         aer, gas, num_a, water_a, water_a_hyst, dp_dry_a, &
         mosaic_vars_aa )

  if (mosaic_vars_aa%f_mos_fail > 0 ) then
     return
  endif
  call check_astem_negative( 3, mosaic_vars_aa%xnerr_astem_negative, &
                                mosaic_vars_aa%fix_astem_negative, aer, gas )

  

  if (msoa_flag1 == 1) then
     itmpa = 1
     call ASTEM_secondary_organics(dtchem,jaerosolstate,sfc_a,Heff,phi_volatile_l,  &
        integrate,aer,kg,gas,sat_soa,total_species) 

  else if (msoa_flag1 == 1000) then
     itmpa = 2
     swdown_cell = mosaic_vars_aa%swdown
     call mosaic_soa_vbs_intr( &
        dtchem, p_atm, t_k, swdown_cell, &
        jaerosolstate, &
        aer, gas, water_a, area_wet_a, dp_wet_a, &
        kg, sat_soa, total_species, &
        ma, mc, mosaic_vars_aa )

  else
     itmpa = 0
  end if

  if (itmpa > 0) then
  call check_astem_negative( 4, mosaic_vars_aa%xnerr_astem_negative, &
                                mosaic_vars_aa%fix_astem_negative, aer, gas )
  end if

  do iv = 1, ngas_aerchtot
     if (iv == ih2so4_g) cycle
     
     gas_avg(iv) = 0.5_r8*(gas_avg(iv) + gas(iv))
  end do

  call dumpxx( 'kk', dtchem, t_k, p_atm, ah2o, &
         jaerosolstate, jphase, jhyst_leg, &
         aer, gas, num_a, water_a, water_a_hyst, dp_dry_a, &
         mosaic_vars_aa )

  return
end subroutine ASTEM




  subroutine check_astem_negative( n, xnerr_astem_negative, &
                                      fix_astem_negative, aer, gas )






  use module_data_mosaic_aero, only: naer, nbin_a, nbin_a_max, ngas_aerchtot

  integer,  intent(in)    :: n, fix_astem_negative
  real(r8), intent(inout) :: xnerr_astem_negative(5,4)
  real(r8), intent(inout), dimension(ngas_aerchtot)     :: gas
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer

  character(len = 100) :: tmp_str
  integer  :: iaer, ibin, igas, j, m
  real(r8) :: tmpa

  if ( n<1 .or. n>4 ) then
     write(tmp_str,'(/a,i10/)') '*** check_astem_negative fatal error, n =', n
     call mosaic_err_mess(tmp_str)
  end if

  do igas = 1, ngas_aerchtot
     tmpa = gas(igas)
     if (tmpa >= 0.0_r8) then
        cycle
     else if (tmpa <= -1.0e-5_r8 ) then
        m = 1
     else if (tmpa <= -1.0e-10_r8) then
        m = 2
     else if (tmpa <= -1.0e-20_r8) then
        m = 3
     else if (tmpa <= -1.0e-30_r8) then
        m = 4
     else
        m = 5
     end if
     xnerr_astem_negative(m,n) = xnerr_astem_negative(m,n) + 1.0_r8
     if (fix_astem_negative > 0) gas(igas) = 0.0_r8
  end do

  do ibin = 1, nbin_a
     do j = 1, 3
        do iaer = 1, naer
           tmpa = aer(iaer,j,ibin)
           if (tmpa >= 0.0_r8) then
              cycle
           else if (tmpa <= -1.0e-5_r8 ) then
              m = 1
           else if (tmpa <= -1.0e-10_r8) then
              m = 2
           else if (tmpa <= -1.0e-20_r8) then
              m = 3
           else if (tmpa <= -1.0e-30_r8) then
              m = 4
           else
              m = 5
           end if
           xnerr_astem_negative(m,n) = xnerr_astem_negative(m,n) + 1.0_r8
           if (fix_astem_negative > 0) aer(iaer,j,ibin) = 0.0_r8
        end do
     enddo
  enddo

  end subroutine check_astem_negative









subroutine ASTEM_semi_volatiles( iprint_input,  dtchem, jaerosolstate,                &
     sfc_a, flux_s, flux_l, Heff, volatile_s, phi_volatile_s,     &
     jphase, aer, kg, gas, jhyst_leg, electrolyte, epercent, kel, activity, mc,            &
     delta_nh3_max, delta_hno3_max, delta_hcl_max,                                &
     num_a, mass_dry_a, mass_wet_a, mass_soluble_a,     &
     vol_dry_a, vol_wet_a, water_a, total_species, tot_cl_in,                         &
     aH2O_a, aH2O, ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c, mw_electrolyte, Keq_ll,  &
     Keq_gl, Keq_sg, Kp_nh4cl, Kp_nh4no3, Keq_nh4cl, MW_c, MW_a, mw_aer_mac,             &
     dens_aer_mac, Keq_sl, growth_factor, MDRH, MDRH_T, molality0, rtol_mesa,           &
     jsalt_present, jsalt_index, jsulf_poor, jsulf_rich, phi_salt_old,      &
     integrate, phi_volatile_l,                                                        &
     kappa_nonelectro, mosaic_vars_aa )

  use module_data_mosaic_aero, only: nbin_a_max, nbin_a,                         &
       ngas_aerchtot, ngas_volatile, nelectrolyte,                               &
       Ncation, naer, mYES, mNO, ngas_ioa, jsolid, jliquid, all_solid, all_liquid,        &
       mixed, no_aerosol, jtotal, jhyst_lo, Nanion, nrxn_aer_gl, nrxn_aer_ll,           &
       nrxn_aer_sg, nrxn_aer_sl, nsalt, MDRH_T_NUM, jsulf_poor_NUM, jsulf_rich_NUM,    &
       jnh4cl, jnh4no3,                                                            &
       iso4_a, inh3_g, ihno3_g, ihcl_g,                                            & 
       mosaic_vars_aa_type

  use module_mosaic_ext, only: do_full_deliquescence,form_electrolytes
  

  
  
  integer, intent(inout) :: iprint_input
  integer, intent(in), dimension(nsalt) :: jsalt_index
  integer, intent(inout), dimension(nsalt) :: jsalt_present
  integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase
  integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg
  integer, intent(in), dimension(jsulf_poor_NUM) :: jsulf_poor
  integer, intent(in), dimension(jsulf_rich_NUM) :: jsulf_rich
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(in)  :: dtchem
  real(r8), intent(in) :: aH2O,rtol_mesa
  real(r8), intent(inout) :: Kp_nh4cl,Kp_nh4no3,Keq_nh4cl
  real(r8), intent(in), dimension(naer) :: mw_aer_mac,dens_aer_mac
  real(r8), intent(in), dimension(Ncation) :: zc,MW_c
  real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
  real(r8), intent(in), dimension(Nanion)  :: za,MW_a
  real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
  real(r8), intent(inout), dimension(nbin_a_max) :: delta_nh3_max,delta_hno3_max
  real(r8), intent(inout), dimension(nbin_a_max) :: delta_hcl_max,water_a,num_a
  real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a,mass_wet_a,MDRH
  real(r8), intent(inout), dimension(nbin_a_max) :: mass_soluble_a,vol_dry_a
  real(r8), intent(inout), dimension(nbin_a_max) :: aH2O_a,vol_wet_a,gam_ratio,growth_factor
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a, total_species
  real(r8), intent(inout) :: tot_cl_in
  real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(nrxn_aer_sl) :: Keq_sl
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(MDRH_T_NUM) :: MDRH_T
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: flux_s,flux_l
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: Heff,volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg,kel
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


  
  
  character(len=500) :: tmp_str
  integer ibin, iv, jp, ieqblm_ASTEM, islow_intermassxfer		
  integer, dimension(nbin_a_max) :: idry_case3a
  
  real(r8) :: dtmax, t_new, t_old, t_out, XT,kelvin_nh4no3
  real(r8) :: sum1, sum2, sum3, sum4, sum4a, sum4b, h_flux_s
  real(r8) :: phi_nh4no3_s, phi_nh4cl_s,kelvin_nh4cl,Keq_nh4no3
  real(r8) :: sumkg_nh3,sumkg_hno3,sumkg_hcl 				
  real(r8), dimension(nbin_a_max) :: kgfrac_nh3,kgfrac_hno3,kgfrac_hcl	
  real(r8), dimension(ngas_volatile) :: sum_phi_volatile_s, sum_phi_volatile_l, sum_phi_volatile
  real(r8), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,df_gas_l
  real(r8), dimension(ngas_volatile,nbin_a_max) :: h_s_i_m
  real(r8), dimension(3,nbin_a_max) :: electrolyte_sum


  
  t_old = 0.0
  t_out = dtchem
  
  
  mosaic_vars_aa%isteps_ASTEM = 0
  do ibin = 1, nbin_a
     mosaic_vars_aa%iter_MESA(ibin) = 0
  enddo


  sumkg_nh3   = 0.0
  sumkg_hno3  = 0.0
  sumkg_hcl   = 0.0
  do ibin = 1, nbin_a
     sumkg_nh3   = sumkg_nh3   + kg(inh3_g,ibin)
     sumkg_hno3  = sumkg_hno3  + kg(ihno3_g,ibin)
     sumkg_hcl   = sumkg_hcl   + kg(ihcl_g,ibin)
  enddo
  do ibin = 1, nbin_a
     kgfrac_nh3(ibin)  = kg(inh3_g,ibin)/sumkg_nh3
     kgfrac_hno3(ibin) = kg(ihno3_g,ibin)/sumkg_hno3
     kgfrac_hcl(ibin)  = kg(ihcl_g,ibin)/sumkg_hcl
  enddo


  
  
  
  
10 mosaic_vars_aa%isteps_ASTEM = mosaic_vars_aa%isteps_ASTEM + 1
  
  
  phi_nh4no3_s = 0.0
  phi_nh4cl_s  = 0.0
  ieqblm_ASTEM = mYES			
  
  do 501 ibin = 1, nbin_a
     
     idry_case3a(ibin) = mNO			
     
     do iv = 1, ngas_ioa
        sfc_a(iv)                  = gas(iv)
        df_gas_s(iv,ibin)          = 0.0
        df_gas_l(iv,ibin)          = 0.0
        flux_s(iv,ibin)            = 0.0
        flux_l(iv,ibin)            = 0.0
        Heff(iv,ibin)              = 0.0
        volatile_s(iv,ibin)        = 0.0
        phi_volatile_s(iv,ibin)    = 0.0
        phi_volatile_l(iv,ibin)    = 0.0
        integrate(iv,jsolid,ibin)  = mNO	
        integrate(iv,jliquid,ibin) = mNO	
     enddo




     if(jaerosolstate(ibin) .ne. no_aerosol)then
        delta_nh3_max(ibin) = 0.1*gas(inh3_g)*kgfrac_nh3(ibin)
        delta_hno3_max(ibin)= 0.1*gas(ihno3_g)*kgfrac_hno3(ibin)
        delta_hcl_max(ibin) = 0.1*gas(ihcl_g)*kgfrac_hcl(ibin)
     endif



     if(jaerosolstate(ibin) .eq. all_solid)then
        jphase(ibin) = jsolid
        call ASTEM_flux_dry(ibin, phi_nh4no3_s, phi_nh4cl_s, ieqblm_ASTEM,           &
             idry_case3a, sfc_a, df_gas_s, flux_s, phi_volatile_s, integrate, aer, kg,   &
             gas, electrolyte, epercent, Keq_sg)

     elseif(jaerosolstate(ibin) .eq. all_liquid)then
        jphase(ibin) = jliquid
        call ASTEM_flux_wet(ibin, ieqblm_ASTEM, sfc_a, df_gas_s, df_gas_l,            &
             jaerosolstate, flux_s, Heff, phi_volatile_s, phi_volatile_l, integrate,   &
             jphase, aer, kg, gas, jhyst_leg, electrolyte, kel, activity, mc,             &
             delta_nh3_max, delta_hno3_max, delta_hcl_max, Keq_nh4cl, Keq_nh4no3,     &
             num_a, electrolyte_sum, mass_dry_a, mass_soluble_a, water_a, aH2O,        &
             kelvin_nh4no3, kelvin_nh4cl, ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a,    &
             na_Ma, nc_Mc, xeq_c, mw_electrolyte, Kp_nh4cl, Kp_nh4no3, Keq_gl, Keq_ll,   &
             MW_c, MW_a, total_species, tot_cl_in, molality0,                          &
             kappa_nonelectro, mosaic_vars_aa                                      )

     elseif(jaerosolstate(ibin) .eq. mixed)then
        call ASTEM_flux_mix(ibin, phi_nh4no3_s, phi_nh4cl_s, ieqblm_ASTEM,           &
             idry_case3a, sfc_a, df_gas_s, df_gas_l, jaerosolstate, flux_s, Heff,       &
             phi_volatile_s, phi_volatile_l, integrate, jphase, aer, kg, gas, jhyst_leg, &
             electrolyte, epercent, kel, activity, mc, delta_nh3_max, delta_hno3_max,   &
             delta_hcl_max, Keq_nh4cl, Keq_nh4no3, num_a, electrolyte_sum, mass_dry_a, &
             mass_soluble_a, water_a, aH2O, kelvin_nh4no3, kelvin_nh4cl, ma, gam,       &
             log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c, mw_electrolyte,     &
             Kp_nh4cl, Kp_nh4no3, Keq_ll, Keq_gl, Keq_sg, MW_c, MW_a, total_species,     &
             tot_cl_in, molality0, kappa_nonelectro, mosaic_vars_aa )	

     endif
     
501 continue

  if(ieqblm_ASTEM .eq. mYES)goto 30	




  islow_intermassxfer = mNO 
  if(mosaic_vars_aa%isteps_ASTEM .gt. 20)then
    islow_intermassxfer = mYes 

    do iv = 2, 4 
      sum_phi_volatile_s(iv) = sum(abs(phi_volatile_s(iv,1:nbin_a)))
      sum_phi_volatile_l(iv) = sum(abs(phi_volatile_l(iv,1:nbin_a)))
      sum_phi_volatile(iv) = sum_phi_volatile_s(iv) + sum_phi_volatile_l(iv)

      if(gas(iv) .gt. 0.01 .and. sum_phi_volatile(iv) .gt. 0.01)islow_intermassxfer = mNO

    enddo
  endif

  if(islow_intermassxfer .eq. mYES)goto 30 



  
  
11 call ASTEM_calculate_dtmax( dtchem,  dtmax, jaerosolstate, idry_case3a, df_gas_s,   &
        flux_s, volatile_s, phi_volatile_l, integrate, aer, kg, gas, electrolyte,        &
        h_s_i_m, mosaic_vars_aa )
  t_new = t_old + dtmax	
  if(t_new .gt. t_out)then	
     dtmax = t_out - t_old
     t_new = t_out*1.01
  endif
  
  
  
  

  do 20 iv = 2, 4
     
     sum1 = 0.0
     sum2 = 0.0
     sum3 = 0.0
     sum4 = 0.0
     sum4a= 0.0
     sum4b= 0.0
     
     do 21 ibin = 1, nbin_a
        if(jaerosolstate(ibin) .eq. no_aerosol)goto 21
        
        jp = jliquid
        sum1 = sum1 + aer(iv,jp,ibin)/   &
             (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin)*integrate(iv,jp,ibin))
        
        sum2 = sum2 + kg(iv,ibin)*integrate(iv,jp,ibin)/   &
             (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin)*integrate(iv,jp,ibin))
        
        jp = jsolid
        sum3 = sum3 + aer(iv,jp,ibin)
        
        if(flux_s(iv,ibin) .gt. 0.)then
           h_flux_s = dtmax*flux_s(iv,ibin)
           sum4a = sum4a + h_flux_s
           aer(iv,jp,ibin) = aer(iv,jp,ibin) + h_flux_s
        elseif(flux_s(iv,ibin) .lt. 0.)then
           h_flux_s = min(h_s_i_m(iv,ibin),dtmax)*flux_s(iv,ibin)
           sum4b = sum4b + h_flux_s
           aer(iv,jp,ibin) = aer(iv,jp,ibin) + h_flux_s
           aer(iv,jp,ibin) = max(aer(iv,jp,ibin), 0.0d0)
        endif
        
21   continue
        
     sum4 = sum4a + sum4b
     
     
     
     gas(iv) = (total_species(iv) - (sum1 + sum3 + sum4) )/   &
          (1. + dtmax*sum2)
     gas(iv) = max(gas(iv), 0.0d0)
     
     
     
     
     do 22 ibin = 1, nbin_a
        
        if(integrate(iv,jliquid,ibin) .eq. mYES)then
           aer(iv,jliquid,ibin) =   &
                (aer(iv,jliquid,ibin) + dtmax*kg(iv,ibin)*gas(iv))/   &
                (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin))
        endif
        
22   continue
        
        
20 continue
  
  
        
        
  
  
  
  
  do 40 ibin = 1, nbin_a
     if(jaerosolstate(ibin) .eq. no_aerosol)goto 40
     
     if(jphase(ibin) .eq. jsolid)then
        call form_electrolytes(jsolid,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)  
     elseif(jphase(ibin) .eq. jliquid)then
        call form_electrolytes(jliquid,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in) 
     elseif(jphase(ibin) .eq. jtotal)then
        call form_electrolytes(jsolid,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)  
        call form_electrolytes(jliquid,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in) 
     endif
     
     
     
     do iv = 2, ngas_ioa
        aer(iv,jtotal,ibin)=aer(iv,jsolid,ibin)+aer(iv,jliquid,ibin)
     enddo
     
     
     
     call form_electrolytes(jtotal,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)	
     
     
     
     
     if(jhyst_leg(ibin) .eq. jhyst_lo)then
        call ASTEM_update_phase_eqblm(ibin, jaerosolstate,     &
             jphase, aer, jhyst_leg, electrolyte, epercent, activity, mc, num_a,         &
             mass_dry_a, mass_wet_a, mass_soluble_a, vol_dry_a, vol_wet_a, water_a,    &
             aH2O_a, aH2O, ma, gam, log_gamZ, zc, za,    &
             gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c, mw_electrolyte, mw_aer_mac,         &
             dens_aer_mac, Keq_sl, MW_c, MW_a, Keq_ll, growth_factor, MDRH, MDRH_T,      &
             molality0, rtol_mesa, jsalt_present, jsalt_index, jsulf_poor, jsulf_rich, &
             phi_salt_old,                                                   &
             kappa_nonelectro, mosaic_vars_aa )
        if (mosaic_vars_aa%f_mos_fail > 0) then
           return
        endif
     else
        call do_full_deliquescence(ibin,aer,electrolyte)		
     endif
     
     
40 continue
  
     
  
  t_old = t_new
  
  if(mosaic_vars_aa%isteps_ASTEM .ge. mosaic_vars_aa%nmax_ASTEM)then     
     mosaic_vars_aa%jASTEM_fail = mosaic_vars_aa%jASTEM_fail + 1
     write(tmp_str,*)'ASTEM internal steps exceeded', mosaic_vars_aa%nmax_ASTEM
     call mosaic_warn_mess(trim(adjustl(tmp_str)))

     write(tmp_str,*)'ibin =', ibin
     call mosaic_warn_mess(trim(adjustl(tmp_str)))

     if(iprint_input .eq. mYES)then
        
        iprint_input = mNO
     endif
     goto 30
  elseif(t_new .lt. t_out)then
     goto 10
  endif
  
  
  
  if(t_new .lt. 0.9999*t_out) goto 10
  
30 mosaic_vars_aa%cumul_steps_ASTEM = mosaic_vars_aa%cumul_steps_ASTEM + mosaic_vars_aa%isteps_ASTEM
   mosaic_vars_aa%isteps_ASTEM_max = max( mosaic_vars_aa%isteps_ASTEM_max, mosaic_vars_aa%isteps_ASTEM )
  
  
  
  
  
  
  
  do ibin = 1, nbin_a
     
     if(jaerosolstate(ibin) .eq. mixed)then
        if( electrolyte(jnh4no3,jsolid,ibin).gt. 0.0 .or.   &
             electrolyte(jnh4cl, jsolid,ibin).gt. 0.0 )then
           call ASTEM_flux_mix(ibin, phi_nh4no3_s, phi_nh4cl_s, ieqblm_ASTEM,        &
                idry_case3a, sfc_a, df_gas_s, df_gas_l, jaerosolstate, flux_s, Heff,    &
                phi_volatile_s, phi_volatile_l, integrate, jphase, aer, kg, gas,        &
                jhyst_leg, electrolyte, epercent, kel, activity, mc, delta_nh3_max,     &
                delta_hno3_max, delta_hcl_max, Keq_nh4cl, Keq_nh4no3, num_a,          &
                electrolyte_sum, mass_dry_a, mass_soluble_a, water_a, aH2O,           &
                kelvin_nh4no3, kelvin_nh4cl, ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, &
                na_Ma, nc_Mc, xeq_c, mw_electrolyte, Kp_nh4cl, Kp_nh4no3, Keq_ll,       &
                Keq_gl, Keq_sg, MW_c, MW_a, total_species, tot_cl_in, molality0,        &
                kappa_nonelectro, mosaic_vars_aa )		
        else
           jphase(ibin) = jliquid
           call ASTEM_flux_wet(ibin, ieqblm_ASTEM, sfc_a, df_gas_s, df_gas_l,         &
                jaerosolstate, flux_s, Heff, phi_volatile_s, phi_volatile_l,          &
                integrate, jphase, aer, kg, gas, jhyst_leg, electrolyte, kel, activity,   &
                mc, delta_nh3_max, delta_hno3_max, delta_hcl_max, Keq_nh4cl,          &
                Keq_nh4no3, num_a, electrolyte_sum, mass_dry_a, mass_soluble_a,       &
                water_a, aH2O, kelvin_nh4no3, kelvin_nh4cl, ma, gam, log_gamZ, zc, za,    &
                gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c, mw_electrolyte, Kp_nh4cl,        &
                Kp_nh4no3, Keq_gl, Keq_ll, MW_c, MW_a, total_species, tot_cl_in,        &
                molality0, kappa_nonelectro, mosaic_vars_aa )
        endif
     endif
     
  enddo
  
  
  return
end subroutine ASTEM_semi_volatiles









subroutine ASTEM_calculate_dtmax( dtchem, dtmax, jaerosolstate, idry_case3a,       &
     df_gas_s, flux_s, volatile_s, phi_volatile_l, integrate, aer, kg, gas, electrolyte, &
     h_s_i_m, mosaic_vars_aa )

  use module_data_mosaic_aero, only: r8, nbin_a_max,                       &
       ngas_aerchtot, ngas_volatile, nelectrolyte,                         &
       naer, ngas_ioa, mYES, jliquid, jsolid, no_aerosol,                  &
       nbin_a, alpha_astem,                                                &
       jnh4no3, ino3_a, jnh4cl, inh4_a, icl_a,                             &
       mosaic_vars_aa_type


  
  
  integer, intent(in), dimension(nbin_a_max) :: jaerosolstate,idry_case3a 
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(in)  :: dtchem
  real(r8), intent(out) :: dtmax
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: flux_s,volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: h_s_i_m
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

  type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

  
  character(len=500) :: tmp_str
  integer ::  ibin, iv
  real(r8) :: alpha, h_gas, h_sub_max,h_gas_i(ngas_ioa), h_gas_l, h_gas_s
  real(r8) :: sum_kg_phi, sum_kg_phi_pos, sum_kg_phi_neg, sumflux_s	
  real(r8), dimension(ngas_volatile) :: sum_bin_s,sum_vdf_s,sum_vol_s
  real(r8), dimension(ngas_volatile) :: avg_df_gas_s  
  
  h_sub_max = dtchem/5.0	
  
  
  
  
  
  
  h_gas_s = 2.e16
  
  do 5 iv = 2, ngas_ioa
     h_gas_i(iv) = 1.e16
     sumflux_s = 0.0
     do ibin = 1, nbin_a
        if(flux_s(iv,ibin) .gt. 0.0)then
           sumflux_s = sumflux_s + flux_s(iv,ibin)
        endif
     enddo
     
     if(sumflux_s .gt. 0.0)then
        h_gas_i(iv) = 0.1*gas(iv)/sumflux_s
        h_gas_s     = min(h_gas_s, h_gas_i(iv))
     endif
     
5 continue
     
     
  
  
     
  h_gas_l = 2.e16
  do 6 iv = 2, ngas_ioa
     h_gas_i(iv) = 1.e16
     sum_kg_phi = 0.0
     sum_kg_phi_pos = 0.0
     sum_kg_phi_neg = 0.0
     do ibin = 1, nbin_a
        if(integrate(iv,jliquid,ibin) .eq. mYES)then





          if(phi_volatile_l(iv,ibin) .gt. 0.0)then
           sum_kg_phi_pos = sum_kg_phi_pos + abs(phi_volatile_l(iv,ibin))*kg(iv,ibin)
          else
           sum_kg_phi_neg = sum_kg_phi_neg + abs(phi_volatile_l(iv,ibin))*kg(iv,ibin)
          endif


        endif
     enddo
     
     sum_kg_phi = max(sum_kg_phi_pos, sum_kg_phi_neg) 

     if(sum_kg_phi .gt. 0.0)then
        h_gas_i(iv) = alpha_astem/sum_kg_phi
        h_gas_l     = min(h_gas_l, h_gas_i(iv))
     endif

6 continue

  h_gas = min(h_gas_s, h_gas_l)
  h_gas = min(h_gas, h_sub_max)
  
  
  
  
  
  
  
  do ibin = 1, nbin_a
     
     volatile_s(ino3_a,ibin) = electrolyte(jnh4no3,jsolid,ibin)
     volatile_s(inh4_a,ibin) = electrolyte(jnh4cl,jsolid,ibin) +   &
          electrolyte(jnh4no3,jsolid,ibin)
     
     if(idry_case3a(ibin) .eq. mYES)then
        volatile_s(icl_a,ibin)  = aer(icl_a,jsolid,ibin)
     else
        volatile_s(icl_a,ibin)  = electrolyte(jnh4cl,jsolid,ibin)
     endif
     
  enddo
  
  
  
  do iv = 2, ngas_ioa
     
     sum_bin_s(iv) = 0.0
     sum_vdf_s(iv) = 0.0
     sum_vol_s(iv) = 0.0
     
     do ibin = 1, nbin_a
        if(flux_s(iv,ibin) .lt. 0.)then	
           sum_bin_s(iv) = sum_bin_s(iv) + 1.0
           sum_vdf_s(iv) = sum_vdf_s(iv) +   &
                volatile_s(iv,ibin)*df_gas_s(iv,ibin)
           sum_vol_s(iv) = sum_vol_s(iv) + volatile_s(iv,ibin)
        endif
     enddo
     
     if(sum_vol_s(iv) .gt. 0.0)then
        avg_df_gas_s(iv) = sum_vdf_s(iv)/sum_vol_s(iv)
     else
        avg_df_gas_s(iv) = 1.0 
     endif
     
  enddo
  
  
  
  
  
  do 20 ibin = 1, nbin_a
     
     if(jaerosolstate(ibin) .eq. no_aerosol) goto 20
     
     do 10 iv = 2, ngas_ioa
        
        if(flux_s(iv,ibin) .lt. 0.)then				
           
           alpha = abs(avg_df_gas_s(iv))/   &
                (volatile_s(iv,ibin)*sum_bin_s(iv))
           alpha = min(alpha, 1.0d0)
           
           if(idry_case3a(ibin) .eq. mYES)alpha = 1.0
           
           h_s_i_m(iv,ibin) =   &
                -alpha*volatile_s(iv,ibin)/flux_s(iv,ibin)
           
        endif
        
10   continue
        
        
20 continue
        
        
  dtmax = min(dtchem, h_gas)



  
  if(dtmax .eq. 0.0)then
     write(tmp_str,*)' dtmax = ', dtmax
     call mosaic_warn_mess(trim(adjustl(tmp_str)))
  endif
  
  return
end subroutine ASTEM_calculate_dtmax














subroutine ASTEM_update_phase_eqblm(ibin, jaerosolstate,       &
     jphase, aer, jhyst_leg, electrolyte, epercent, activity, mc, num_a, mass_dry_a,    &
     mass_wet_a, mass_soluble_a, vol_dry_a, vol_wet_a, water_a, aH2O_a, aH2O,           &
     ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc,                         &
     xeq_c, mw_electrolyte, mw_aer_mac, dens_aer_mac, Keq_sl, MW_c, MW_a, Keq_ll,       &
     growth_factor, MDRH, MDRH_T, molality0, rtol_mesa, jsalt_present, jsalt_index,     &
     jsulf_poor, jsulf_rich, phi_salt_old,                                              &
     kappa_nonelectro, mosaic_vars_aa )

  use module_data_mosaic_aero,  only: r8, nbin_a_max, nelectrolyte, Ncation, naer,      &
       jtotal, nsalt, all_solid, jsolid, all_liquid, jliquid, jhyst_lo, jhyst_up,       &
       Nanion, nrxn_aer_ll, nrxn_aer_sl, MDRH_T_NUM,  &
       jsulf_poor_NUM, jsulf_rich_NUM,                                                  &
       ptol_mol_astem, mhyst_force_lo,  mhyst_force_up,                                 &
       jcacl2, jcano3, mhyst_method,                                                    &
       mosaic_vars_aa_type

  use module_mosaic_ext,  only: do_full_deliquescence, adjust_solid_aerosol,            &
       MESA_PTC, calculate_XT, aerosol_water, adjust_liquid_aerosol,                    &
       compute_activities
  
  
  
  integer, intent(in):: ibin
  integer, intent(in), dimension(nsalt) :: jsalt_index
  integer, intent(inout), dimension(nsalt) :: jsalt_present
  integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase,jhyst_leg
  integer, intent(in), dimension(jsulf_poor_NUM) :: jsulf_poor
  integer, intent(in), dimension(jsulf_rich_NUM) :: jsulf_rich
  
  real(r8), intent(in) :: aH2O,rtol_mesa
  real(r8), intent(in), dimension(naer) :: mw_aer_mac,dens_aer_mac
  real(r8), intent(in), dimension(Ncation) :: zc,MW_c
  real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
  real(r8), intent(in), dimension(Nanion)  :: za,MW_a
  real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
  real(r8), intent(inout), dimension(nbin_a_max) :: num_a,mass_dry_a,mass_wet_a
  real(r8), intent(inout), dimension(nbin_a_max) :: mass_soluble_a,vol_dry_a
  real(r8), intent(inout), dimension(nbin_a_max) :: vol_wet_a,water_a,gam_ratio
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

  
  integer jsalt_dum, js, j_index, je
  real(r8) :: CRH, XT, sum_dum
  
  
  
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

    
    if(mhyst_method == mhyst_force_up .or. jhyst_leg(ibin) == jhyst_up) then 
       call do_full_deliquescence(ibin,aer,electrolyte) 
       jaerosolstate(ibin) = all_liquid
       jhyst_leg(ibin) = jhyst_up
       jphase(ibin) = jliquid
       water_a(ibin) = aerosol_water(jtotal,ibin,jaerosolstate,jphase,jhyst_leg,   &
          electrolyte,aer,kappa_nonelectro,num_a,mass_dry_a,mass_soluble_a,aH2O,molality0)

       if(water_a(ibin) .le. 0.0)then     
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

       return
    endif


    
    if(aH2O*100. .lt. MDRH(ibin)) then
       jaerosolstate(ibin) = all_solid
       jphase(ibin) = jsolid
       call adjust_solid_aerosol(ibin,jphase,aer,jhyst_leg,electrolyte,epercent,water_a)
       return
    endif
  
  
    
10  if(jphase(ibin) .eq. jsolid)then
      call do_full_deliquescence(ibin,aer,electrolyte)
      call MESA_PTC( ibin, jaerosolstate, jphase, aer,          &
          jhyst_leg, electrolyte, epercent, activity, mc, num_a, mass_dry_a, mass_wet_a, &
          mass_soluble_a, vol_dry_a, vol_wet_a, water_a, aH2O,          &
          ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma,     &
          nc_Mc, xeq_c, mw_electrolyte, mw_aer_mac, dens_aer_mac, Keq_sl, MW_c, MW_a,    &
          Keq_ll, growth_factor, molality0, rtol_mesa, jsalt_present,                  &
          phi_salt_old, kappa_nonelectro, mosaic_vars_aa )
    else
      call MESA_PTC( ibin, jaerosolstate, jphase, aer,          &
          jhyst_leg, electrolyte, epercent, activity, mc, num_a, mass_dry_a, mass_wet_a, &
          mass_soluble_a, vol_dry_a, vol_wet_a, water_a, aH2O,          &
          ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma,     &
          nc_Mc, xeq_c, mw_electrolyte, mw_aer_mac, dens_aer_mac, Keq_sl, MW_c, MW_a,    &
          Keq_ll, growth_factor, molality0, rtol_mesa, jsalt_present,                  &
          phi_salt_old, kappa_nonelectro, mosaic_vars_aa )
    endif  
    return

  end subroutine ASTEM_update_phase_eqblm













subroutine ASTEM_flux_wet(ibin, ieqblm_ASTEM, sfc_a, df_gas_s, df_gas_l,              &
     jaerosolstate, flux_s, Heff, phi_volatile_s, phi_volatile_l, integrate, jphase,    &
     aer, kg, gas, jhyst_leg, electrolyte, kel, activity, mc, delta_nh3_max,              &
     delta_hno3_max, delta_hcl_max, Keq_nh4cl, Keq_nh4no3, num_a, electrolyte_sum,     &
     mass_dry_a, mass_soluble_a, water_a, aH2O, kelvin_nh4no3, kelvin_nh4cl, ma, gam,    &
     log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc, xeq_c, mw_electrolyte, Kp_nh4cl,    &
     Kp_nh4no3, Keq_gl, Keq_ll, MW_c, MW_a, total_species, tot_cl_in, molality0,         &
     kappa_nonelectro, mosaic_vars_aa                                            )

  use module_data_mosaic_aero, only: r8, nbin_a_max,                             &
       ngas_aerchtot, ngas_volatile, nelectrolyte,                               &
       Ncation, naer, jliquid, jsolid, mNO, mYES, Nanion, nrxn_aer_gl, nrxn_aer_ll,   &
       jcaco3, inh4_a, inh3_g, ihno3_g, ino3_a, ihcl_g, icl_a, jnh4no3, jnh4cl,       &
       mosaic_vars_aa_type

  use module_mosaic_ext,  only: compute_activities, ions_to_electrolytes,           &
       absorb_tiny_nh4no3, absorb_tiny_nh4cl, absorb_tiny_hno3, absorb_tiny_hcl
  
  
  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase,jhyst_leg
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(in) :: aH2O
  real(r8), intent(inout) :: Keq_nh4cl,Keq_nh4no3,kelvin_nh4no3,kelvin_nh4cl
  real(r8), intent(inout) :: Kp_nh4cl,Kp_nh4no3
  real(r8), intent(in), dimension(Ncation) :: zc,MW_c
  real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
  real(r8), intent(in), dimension(Nanion)  :: za,MW_a
  real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
  real(r8), intent(inout), dimension(nbin_a_max) :: delta_nh3_max,delta_hno3_max
  real(r8), intent(inout), dimension(nbin_a_max) :: delta_hcl_max
  real(r8), intent(inout), dimension(nbin_a_max) :: num_a,mass_dry_a,gam_ratio
  real(r8), intent(inout), dimension(nbin_a_max) :: mass_soluble_a,water_a
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a, total_species
  real(r8), intent(inout) :: tot_cl_in
  real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_l
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: flux_s,Heff
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg,kel
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
  real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  real(r8), intent(inout), dimension(3,nbin_a_max) :: electrolyte_sum
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(in), dimension(naer) :: kappa_nonelectro

  type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

  
  character(len=500) :: tmp_str
  integer iv, iadjust, iadjust_intermed
  real(r8) :: XT, g_nh3_hno3, g_nh3_hcl, a_nh4_no3, a_nh4_cl

  call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,   &
       nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)  	
  call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg, &
       electrolyte,activity,mc,num_a,mass_dry_a,mass_soluble_a,    &
       water_a,aH2O,ma,gam,log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)
  
  if(water_a(ibin) .eq. 0.0)then
     write(tmp_str,*)'Water is zero in liquid phase'
     call mosaic_warn_mess(trim(adjustl(tmp_str)))
     write(tmp_str,*)'Stopping in ASTEM_flux_wet'     
     call mosaic_warn_mess(trim(adjustl(tmp_str)))
     mosaic_vars_aa%zero_water_flag = .true.
  endif
  
  
  
  
  if(electrolyte(jcaco3,jsolid,ibin) .gt. 0.0)then
     call ASTEM_flux_wet_case1(ibin,ieqblm_ASTEM,sfc_a,df_gas_s,flux_s,          &
          phi_volatile_s,integrate,jphase,kg,gas,mc,Keq_ll)
     return
  endif
  
  
  
  

  if(XT.lt.2.0 .and. XT.ge.0.)then  
     call ASTEM_flux_wet_case2(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,            &
          phi_volatile_l,integrate,gas,kel,mc,water_a,ma,gam,gam_ratio,Keq_ll,   &
          Keq_gl)
     return
  endif
  
  
  
  if( (gas(inh3_g)+aer(inh4_a,jliquid,ibin)) .lt. 1.e-25)goto 10  
  
  
  
  
  
  iadjust = mNO		
  iadjust_intermed = mNO	
  
  
  g_nh3_hno3 = gas(inh3_g)*gas(ihno3_g)
  a_nh4_no3  = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)
  
  if(g_nh3_hno3 .gt. 0. .and. a_nh4_no3 .eq. 0.)then
     call absorb_tiny_nh4no3(ibin,aer,gas,electrolyte,delta_nh3_max,             &
          delta_hno3_max,electrolyte_sum)
     iadjust = mYES
     iadjust_intermed = mYES
  endif
  
  if(iadjust_intermed .eq. mYES)then
     call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,&
          nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)  	
     iadjust_intermed = mNO	
  endif
  
  
  g_nh3_hcl = gas(inh3_g)*gas(ihcl_g)
  a_nh4_cl  = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)
  
  if(g_nh3_hcl .gt. 0. .and. a_nh4_cl .eq. 0.)then
     call absorb_tiny_nh4cl(ibin,aer,gas,electrolyte,delta_nh3_max,delta_hcl_max,&
          electrolyte_sum)
     iadjust = mYES
     iadjust_intermed = mYES
  endif
  
  if(iadjust_intermed .eq. mYES)then
     call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,&
          nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)  	
  endif
  
  if(iadjust .eq. mYES)then
     call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,electrolyte,&
          activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,ma,gam,       &
          log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)			
  endif
  
  
  
  
  
  kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
  Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	
  
  kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
  Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl	
  
  call ASTEM_flux_wet_case3(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,phi_volatile_l,&
       integrate,kg,gas,kel,mc,Keq_nh4cl,Keq_nh4no3,water_a,ma,gam,gam_ratio,    &
       Keq_ll,Keq_gl,aer,total_species,tot_cl_in,activity,electrolyte)
  
  return
  
  
  
  
  
  
10 iadjust = mNO		
  iadjust_intermed = mNO	
  
  
  if(gas(ihno3_g).gt.0. .and. aer(ino3_a,jliquid,ibin).eq.0. .and.   &
       aer(icl_a,jliquid,ibin) .gt. 0.0)then
     call absorb_tiny_hno3(ibin,aer,gas,delta_hno3_max)	
     iadjust = mYES
     iadjust_intermed = mYES
  endif
  
  if(iadjust_intermed .eq. mYES)then
     call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,&
          nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)  	
     iadjust_intermed = mNO	
  endif
  
  
  if(gas(ihcl_g).gt.0. .and. aer(icl_a,jliquid,ibin) .eq. 0. .and.   &
       aer(ino3_a,jliquid,ibin) .gt. 0.0)then
     call absorb_tiny_hcl(ibin,aer,gas,delta_hcl_max)	
     iadjust = mYES
     iadjust_intermed = mYES
  endif
  
  if(iadjust_intermed .eq. mYES)then
     call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,&
          nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)  	
  endif
  
  if(iadjust .eq. mYES)then
     call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,electrolyte,&
          activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,ma,gam,       &
          log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)			
  endif
  
  
  
  call ASTEM_flux_wet_case4(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,phi_volatile_l,&
       integrate,kg,gas,kel,mc,water_a,ma,gam,Keq_ll,Keq_gl)

  
  return
end subroutine ASTEM_flux_wet












subroutine ASTEM_flux_wet_case1(ibin,ieqblm_ASTEM,sfc_a,df_gas_s,flux_s,         &
     phi_volatile_s,integrate,jphase,kg,gas,mc,Keq_ll)

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       Ncation,mYES, jsolid,mNO,nrxn_aer_ll,                                       &
       jc_h,ihno3_g,ihcl_g
  
  
  
  integer, intent(in):: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(nbin_a_max) :: jphase
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
  
  
  integer iv
  
  mc(jc_h,ibin) = sqrt(Keq_ll(3))
  
  
  if(gas(ihno3_g) .gt. 1.e-6)then
     sfc_a(ihno3_g) = 0.0
     df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
     phi_volatile_s(ihno3_g,ibin) = 1.0
     flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
     integrate(ihno3_g,jsolid,ibin) = mYES
     jphase(ibin) = jsolid
     ieqblm_ASTEM = mNO
  endif
  
  if(gas(ihcl_g) .gt. 1.e-6)then
     sfc_a(ihcl_g)  = 0.0
     df_gas_s(ihcl_g,ibin) = gas(ihcl_g)
     phi_volatile_s(ihcl_g,ibin) = 1.0
     flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
     integrate(ihcl_g,jsolid,ibin)  = mYES
     jphase(ibin) = jsolid
     ieqblm_ASTEM = mNO
  endif
  
  return
end subroutine ASTEM_flux_wet_case1






subroutine ASTEM_flux_wet_case2(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,           &
     phi_volatile_l,integrate,gas,kel,mc,water_a,ma,gam,gam_ratio,Keq_ll,Keq_gl)

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       Ncation,mYES, jliquid,mNO,Nanion,nelectrolyte,nrxn_aer_gl,nrxn_aer_ll,      &
       jc_h,jc_nh4,inh3_g,jhno3,ja_no3,ihno3_g,jhcl,ja_cl,ihcl_g

  
  
  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a,gam_ratio
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_l,Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: gam
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  
  real(r8) :: dum_hno3, dum_hcl, dum_nh3
  
  
  sfc_a(inh3_g)  = kel(inh3_g,ibin)*   &
       gam_ratio(ibin)*mc(jc_nh4,ibin)*Keq_ll(3)/   &
       (mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
  
  sfc_a(ihno3_g) = kel(ihno3_g,ibin)*   &
       mc(jc_h,ibin)*ma(ja_no3,ibin)*gam(jhno3,ibin)**2/   &
       Keq_gl(3)
  
  sfc_a(ihcl_g)  = kel(ihcl_g,ibin)*   &
       mc(jc_h,ibin)*ma(ja_cl,ibin)*gam(jhcl,ibin)**2/   &
       Keq_gl(4)
  
  dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
  dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
  dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))
  
  
  
  if(dum_hno3 .gt. 0.0)then
     df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
     phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
  else
     phi_volatile_l(ihno3_g,ibin)= 0.0
  endif
  
  if(dum_hcl .gt. 0.0)then
     df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
     phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
  else
     phi_volatile_l(ihcl_g,ibin) = 0.0
  endif
  
  if(dum_nh3 .gt. 0.0)then
     df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
     phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
  else
     phi_volatile_l(inh3_g,ibin) = 0.0
  endif
  
  
  
  
  
  
  
  
  
  
  
  
  if(dum_hno3 .gt. 0.0)then
     Heff(ihno3_g,ibin)=   &
          kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
          (water_a(ibin)*Keq_gl(3))
     integrate(ihno3_g,jliquid,ibin)= mYES
     ieqblm_ASTEM = mNO
  endif
  
  if(dum_hcl .gt. 0.0)then
     Heff(ihcl_g,ibin)=   &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
          (water_a(ibin)*Keq_gl(4))
     integrate(ihcl_g,jliquid,ibin) = mYES
     ieqblm_ASTEM = mNO
  endif
  
  if(dum_nh3 .gt. 0.0)then
     Heff(inh3_g,ibin) =   &
          kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/   &
          (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
     integrate(inh3_g,jliquid,ibin) = mYES
     ieqblm_ASTEM = mNO
  endif
  
  
  return
end subroutine ASTEM_flux_wet_case2






subroutine ASTEM_flux_wet_case3(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,           &
     phi_volatile_l,integrate,kg,gas,kel,mc,Keq_nh4cl,Keq_nh4no3,water_a,ma,gam, &
     gam_ratio,Keq_ll,Keq_gl,aer,total_species,tot_cl_in,activity,electrolyte)

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       Ncation,mYES, jliquid,mNO,Nanion,nelectrolyte,nrxn_aer_gl,nrxn_aer_ll,naer, &
       inh3_g,ihcl_g,ihno3_g,ja_no3,jhno3,jc_h,ja_cl,jhcl,jc_nh4

  use module_mosaic_ext, only: quadratic,equilibrate_acids
  
  
  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(inout) :: Keq_nh4cl,Keq_nh4no3
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a,gam_ratio
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a, total_species
  real(r8), intent(inout) :: tot_cl_in
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_l,Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg,kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: gam,activity
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  
  real(r8) :: a, b, c, dum_hno3, dum_hcl, dum_nh3
  
  
  
  a =   kg(inh3_g,ibin)
  b = - kg(inh3_g,ibin)*gas(inh3_g)   &
       + kg(ihno3_g,ibin)*gas(ihno3_g)   &
       + kg(ihcl_g,ibin)*gas(ihcl_g)
  c = -(kg(ihno3_g,ibin)*Keq_nh4no3 + kg(ihcl_g,ibin)*Keq_nh4cl)
  
  sfc_a(inh3_g)  = quadratic(a,b,c)
  sfc_a(ihno3_g) = Keq_nh4no3/max(sfc_a(inh3_g),1.d-20)
  sfc_a(ihcl_g)  = Keq_nh4cl/max(sfc_a(inh3_g),1.d-20)
  
  
  
  if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
     mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
  elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
     mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
  else
     call equilibrate_acids(ibin,aer,gas,electrolyte,activity,mc,water_a,       &
       total_species,tot_cl_in,ma,gam,Keq_ll,Keq_gl)	
     mc(jc_h,ibin)  = max(mc(jc_h,ibin), sqrt(Keq_ll(3)))
     
     sfc_a(inh3_g)  = kel(inh3_g,ibin)*   &
          gam_ratio(ibin)*mc(jc_nh4,ibin)*Keq_ll(3)/   &
          (mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
     
     sfc_a(ihno3_g) = kel(ihno3_g,ibin)*   &
          mc(jc_h,ibin)*ma(ja_no3,ibin)*gam(jhno3,ibin)**2/   &
          Keq_gl(3)
     sfc_a(ihcl_g)  = kel(ihcl_g,ibin)*   &
          mc(jc_h,ibin)*ma(ja_cl,ibin)*gam(jhcl,ibin)**2/   &
          Keq_gl(4)
  endif
  
  dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
  dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
  dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))
  
  
  if(dum_hno3 .gt. 0.0)then
     df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
     phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
  else
     phi_volatile_l(ihno3_g,ibin)= 0.0
  endif
  
  if(dum_hcl .gt. 0.0)then
     df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
     phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
  else
     phi_volatile_l(ihcl_g,ibin) = 0.0
  endif
  
  if(dum_nh3 .gt. 0.0)then
     df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
     phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
  else
     phi_volatile_l(inh3_g,ibin) = 0.0
  endif
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(dum_hno3 .gt. 0.0)then
     Heff(ihno3_g,ibin)=   &
          kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
          (water_a(ibin)*Keq_gl(3))
     integrate(ihno3_g,jliquid,ibin)= mYES
     ieqblm_ASTEM = mNO
  endif
  
  if(dum_hcl .gt. 0.0)then
     Heff(ihcl_g,ibin)=   &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
          (water_a(ibin)*Keq_gl(4))
     integrate(ihcl_g,jliquid,ibin) = mYES
     ieqblm_ASTEM = mNO
  endif
  
  if(dum_nh3 .gt. 0.0)then
     Heff(inh3_g,ibin) =   &
          kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/   &
          (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
     integrate(inh3_g,jliquid,ibin) = mYES
     ieqblm_ASTEM = mNO
  endif    
  
  return
end subroutine ASTEM_flux_wet_case3
      





subroutine ASTEM_flux_wet_case3a(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,          &
     phi_volatile_l,integrate,kg,gas,kel,mc,Keq_nh4no3,water_a,ma,gam,gam_ratio, &
     Keq_ll,Keq_gl)    

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       Ncation,mYES, jliquid,mNO,Nanion,nelectrolyte,nrxn_aer_gl,nrxn_aer_ll,      &
       inh3_g,ihno3_g,ja_no3,jhno3,jc_h
  
  use module_mosaic_ext, only: quadratic


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(inout) :: Keq_nh4no3
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a,gam_ratio
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_l,Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg,kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: gam
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  
  real(r8) :: a, b, c, dum_hno3, dum_nh3
  
  


  a =   kg(inh3_g,ibin)
  b = - kg(inh3_g,ibin)*gas(inh3_g)   &
       + kg(ihno3_g,ibin)*gas(ihno3_g)
  c = -(kg(ihno3_g,ibin)*Keq_nh4no3)

  sfc_a(inh3_g)  = quadratic(a,b,c)
  sfc_a(ihno3_g) = Keq_nh4no3/sfc_a(inh3_g)


  
  if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
     mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
  else
     mc(jc_h,ibin) = sqrt(Keq_ll(3))
  endif


  
  dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
  dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))

  
  if(dum_hno3 .gt. 0.0)then
     df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
     phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
  else
     phi_volatile_l(ihno3_g,ibin)= 0.0
  endif

  if(dum_nh3 .gt. 0.0)then
     df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
     phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
  else
     phi_volatile_l(inh3_g,ibin) = 0.0
  endif


  
  
  
  
  
  


  
  Heff(ihno3_g,ibin)=   &
       kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
       (water_a(ibin)*Keq_gl(3))
  integrate(ihno3_g,jliquid,ibin)= mYES


  Heff(inh3_g,ibin) =   &
       kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/   &
       (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
  integrate(inh3_g,jliquid,ibin) = mYES


  ieqblm_ASTEM = mNO


  return
end subroutine ASTEM_flux_wet_case3a






subroutine ASTEM_flux_wet_case3b(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,          &
     phi_volatile_l,integrate,kg,gas,kel,mc,Keq_nh4cl,water_a,ma,gam,gam_ratio,  &
     Keq_ll,Keq_gl)     

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       Ncation,mYES, jliquid,mNO,Nanion,nelectrolyte,nrxn_aer_gl,nrxn_aer_ll,      &
       inh3_g,ihcl_g,ja_cl,jhcl,jc_h

  use module_mosaic_ext, only: quadratic


  
  integer,intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(inout) :: Keq_nh4cl
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a,gam_ratio
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_l,Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg,kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: gam
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  
  real(r8) :: a, b, c, dum_hcl, dum_nh3
  
  


  a =   kg(inh3_g,ibin)
  b = - kg(inh3_g,ibin)*gas(inh3_g)   &
       + kg(ihcl_g,ibin)*gas(ihcl_g)
  c = -(kg(ihcl_g,ibin)*Keq_nh4cl)

  sfc_a(inh3_g)  = quadratic(a,b,c)
  sfc_a(ihcl_g)  = Keq_nh4cl /sfc_a(inh3_g)


  
  if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
     mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
  else
     mc(jc_h,ibin) = sqrt(Keq_ll(3))
  endif


  
  dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
  dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))


  
  if(dum_hcl .gt. 0.0)then
     df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
     phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
  else
     phi_volatile_l(ihcl_g,ibin) = 0.0
  endif

  if(dum_nh3 .gt. 0.0)then
     df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
     phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
  else
     phi_volatile_l(inh3_g,ibin) = 0.0
  endif



  
  
  
  
  
  



  
  Heff(ihcl_g,ibin)=   &
       kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
       (water_a(ibin)*Keq_gl(4))
  integrate(ihcl_g,jliquid,ibin) = mYES


  Heff(inh3_g,ibin) =   &
       kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/   &
       (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
  integrate(inh3_g,jliquid,ibin) = mYES


  ieqblm_ASTEM = mNO



  return
end subroutine ASTEM_flux_wet_case3b






subroutine ASTEM_flux_wet_case4(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,           &
     phi_volatile_l,integrate,kg,gas,kel,mc,water_a,ma,gam,Keq_ll,Keq_gl)

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       Ncation,mYES, jliquid,mNO,Nanion,nelectrolyte,nrxn_aer_gl,nrxn_aer_ll,      &
       jhno3,ja_no3,ihno3_g,jhcl,ja_cl,ihcl_g,jc_h


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(inout), dimension(nbin_a_max) :: water_a
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_l,Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg,kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max):: gam
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) ::mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  
  real(r8) :: dum_numer, dum_denom, gas_eqb_ratio, dum_hno3, dum_hcl


  dum_numer = kel(ihno3_g,ibin)*Keq_gl(4)*ma(ja_no3,ibin)*   &
       gam(jhno3,ibin)**2
  dum_denom = kel(ihcl_g,ibin)*Keq_gl(3)*ma(ja_cl ,ibin)*   &
       gam(jhcl,ibin)**2


  if(dum_denom .eq. 0.0 .or. dum_numer .eq. 0.0)then
     mc(jc_h,ibin) = sqrt(Keq_ll(3))
     return
  endif

  gas_eqb_ratio = dum_numer/dum_denom   


  
  sfc_a(ihcl_g) =   &
       ( kg(ihno3_g,ibin)*gas(ihno3_g) + kg(ihcl_g,ibin)*gas(ihcl_g) )/   &
       ( kg(ihcl_g,ibin) + gas_eqb_ratio*kg(ihno3_g,ibin) )
  sfc_a(ihno3_g)= gas_eqb_ratio*sfc_a(ihcl_g)


  
  if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
     mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
  elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
     mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
  else
     mc(jc_h,ibin) = sqrt(Keq_ll(3))
  endif


  
  dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
  dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))

  
  if(dum_hno3 .gt. 0.0)then
     df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
     phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
  else
     phi_volatile_l(ihno3_g,ibin)= 0.0
  endif

  if(dum_hcl .gt. 0.0)then
     df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
     phi_volatile_l(ihcl_g,ibin)= df_gas_l(ihcl_g,ibin)/dum_hcl
  else
     phi_volatile_l(ihcl_g,ibin)= 0.0
  endif


  
  
  
  
  
  



  
  Heff(ihno3_g,ibin)=   &
       kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
       (water_a(ibin)*Keq_gl(3))
  integrate(ihno3_g,jliquid,ibin)= mYES


  Heff(ihcl_g,ibin)=   &
       kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
       (water_a(ibin)*Keq_gl(4))
  integrate(ihcl_g,jliquid,ibin) = mYES


  ieqblm_ASTEM = mNO



  return
end subroutine ASTEM_flux_wet_case4














subroutine ASTEM_flux_dry(ibin, phi_nh4no3_s,phi_nh4cl_s,ieqblm_ASTEM,           &
     idry_case3a,sfc_a,df_gas_s,flux_s,phi_volatile_s,integrate,aer,kg,gas,      &
     electrolyte,epercent,Keq_sg)
  
  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       naer,jsolid, nrxn_aer_sg,nelectrolyte,                                      &
       jcaco3,jcacl2,jnacl,ihno3_g,jnh4cl,ihcl_g,inh3_g,jnh4no3

  use module_mosaic_ext, only: calculate_XT


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(nbin_a_max) :: idry_case3a
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(out) :: phi_nh4no3_s,phi_nh4cl_s
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
  
  integer iv
  real(r8) :: XT, prod_nh4no3, prod_nh4cl, volatile_cl



  call calculate_XT(ibin,jsolid,XT,aer)

  
  

  if(electrolyte(jcaco3,jsolid,ibin) .gt. 0.0)then
     call ASTEM_flux_dry_case1(ibin,ieqblm_ASTEM,sfc_a,df_gas_s,flux_s,          &
          phi_volatile_s,integrate,kg,gas)

     return
  endif

  
  


  if(XT.lt.2.0 .and. XT.ge.0.)then   
     call ASTEM_flux_dry_case2(ibin,ieqblm_ASTEM,sfc_a,df_gas_s,flux_s,          &
          phi_volatile_s,integrate,kg,gas)

     return
  endif

  
  

  volatile_cl  = electrolyte(jnacl,jsolid,ibin) +   &
       electrolyte(jcacl2,jsolid,ibin)


  if(volatile_cl .gt. 0.0 .and. gas(ihno3_g).gt. 0.0 )then

     call ASTEM_flux_dry_case3a(ibin,ieqblm_ASTEM,idry_case3a,sfc_a,df_gas_s,    &
          flux_s,phi_volatile_s,integrate,aer,kg,gas)

     prod_nh4cl = max( (gas(inh3_g)*gas(ihcl_g)-Keq_sg(2)), 0.0d0) +   &
          electrolyte(jnh4cl, jsolid,ibin)

     if(prod_nh4cl .gt. 0.0)then
        call ASTEM_flux_dry_case3b(ibin,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s, &
             flux_s,phi_volatile_s,integrate,aer,kg,gas,electrolyte,epercent,    &
             Keq_sg)
     endif

     return
  endif

  
  

  prod_nh4no3 = max( (gas(inh3_g)*gas(ihno3_g)-Keq_sg(1)), 0.0d0) +   &
       electrolyte(jnh4no3,jsolid,ibin)
  prod_nh4cl  = max( (gas(inh3_g)*gas(ihcl_g) -Keq_sg(2)), 0.0d0) +   &
       electrolyte(jnh4cl, jsolid,ibin)

  if(prod_nh4no3 .gt. 0.0 .or. prod_nh4cl .gt. 0.0)then
     call ASTEM_flux_dry_case4(ibin,phi_nh4no3_s,phi_nh4cl_s,ieqblm_ASTEM,sfc_a, &
          df_gas_s,flux_s,phi_volatile_s,integrate,kg,gas,electrolyte,epercent,  &
          Keq_sg,aer)
     return
  endif

  

  return
end subroutine ASTEM_flux_dry














subroutine ASTEM_flux_dry_case1(ibin,ieqblm_ASTEM,sfc_a,df_gas_s,flux_s,         &
     phi_volatile_s,integrate,kg,gas)
  
  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       mYES,jsolid,mNO, ihno3_g,ihcl_g


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg


  if(gas(ihno3_g) .gt. 1.e-6)then
     sfc_a(ihno3_g) = 0.0
     df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
     phi_volatile_s(ihno3_g,ibin) = 1.0
     flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
     integrate(ihno3_g,jsolid,ibin) = mYES
     ieqblm_ASTEM = mNO
  endif

  if(gas(ihcl_g) .gt. 1.e-6)then
     sfc_a(ihcl_g)  = 0.0
     df_gas_s(ihcl_g,ibin) = gas(ihcl_g)
     phi_volatile_s(ihcl_g,ibin) = 1.0
     flux_s(ihcl_g,ibin)  = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
     integrate(ihcl_g,jsolid,ibin)  = mYES
     ieqblm_ASTEM = mNO
  endif


  return
end subroutine ASTEM_flux_dry_case1






subroutine ASTEM_flux_dry_case2(ibin,ieqblm_ASTEM,sfc_a,df_gas_s,flux_s,         &
     phi_volatile_s,integrate,kg,gas) 

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       mYES,jsolid,mNO, inh3_g



  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg


  if(gas(inh3_g).gt.1.e-6)then
     sfc_a(inh3_g) = 0.0
     df_gas_s(inh3_g,ibin) = gas(inh3_g)
     phi_volatile_s(inh3_g,ibin)  = 1.0
     flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*gas(inh3_g)
     integrate(inh3_g,jsolid,ibin) = mYES
     ieqblm_ASTEM = mNO
  endif


  return
end subroutine ASTEM_flux_dry_case2






subroutine ASTEM_flux_dry_case3a(ibin,ieqblm_ASTEM,idry_case3a,sfc_a,df_gas_s,   &
     flux_s,phi_volatile_s,integrate,aer,kg,gas)

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       naer,jsolid,mYES,mNO,                                                       &
       ihno3_g,icl_a,ihcl_g


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(nbin_a_max) :: idry_case3a  
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer


  if(gas(ihno3_g) .gt. 1.e-6)then
     sfc_a(ihno3_g) = 0.0
     sfc_a(ihcl_g)  = gas(ihcl_g) + aer(icl_a,jsolid,ibin)

     df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
     df_gas_s(ihcl_g,ibin)  = -aer(icl_a,jsolid,ibin)

     flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*gas(ihno3_g)
     flux_s(ihcl_g,ibin)  = -flux_s(ihno3_g,ibin)

     phi_volatile_s(ihno3_g,ibin) = 1.0
     phi_volatile_s(ihcl_g,ibin)=df_gas_s(ihcl_g,ibin)/sfc_a(ihcl_g)

     integrate(ihno3_g,jsolid,ibin) = mYES
     integrate(ihcl_g,jsolid,ibin)  = mYES

     idry_case3a(ibin) = mYES
     ieqblm_ASTEM = mNO
  endif

  return
end subroutine ASTEM_flux_dry_case3a







subroutine ASTEM_flux_dry_case3b(ibin, phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s,  &
     flux_s,phi_volatile_s,integrate,aer,kg,gas,electrolyte,epercent,Keq_sg)      

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       naer,nsalt, jsolid,mYES,mNO,nrxn_aer_sg,                                    &
       rtol_eqb_ASTEM,ptol_mol_ASTEM,                                              &
       nelectrolyte,jnh4cl,ihcl_g,inh3_g,icl_a

  use module_mosaic_ext, only: quadratic,degas_solid_nh4cl
  

  
  integer, intent(in) :: ibin
  integer, intent(inout) ::ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(out) :: phi_nh4cl_s
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
  
  integer iactive_nh4cl, js
  real(r8) :: a, b, c
  real(r8) :: sum_dum
  
  



  
  sum_dum = 0.0
  do js = 1, nsalt
     sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
  enddo

  if(sum_dum .eq. 0.)sum_dum = 1.0

  epercent(jnh4cl,jsolid,ibin) = 100.*electrolyte(jnh4cl,jsolid,ibin)/sum_dum
  



  
  
  iactive_nh4cl  = 1


  
  phi_nh4cl_s = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/   &
       max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))


  
  
  
  if( abs(phi_nh4cl_s) .lt. rtol_eqb_ASTEM )then
     iactive_nh4cl = 0
  elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and.   &
       epercent(jnh4cl, jsolid,ibin) .le. ptol_mol_ASTEM)then
     iactive_nh4cl = 0
     if(epercent(jnh4cl, jsolid,ibin) .gt. 0.0)then
        call degas_solid_nh4cl(ibin,aer,gas,electrolyte,Keq_sg)
     endif
  endif


  
  if(iactive_nh4cl .eq. 0)return


  
  


  a =   kg(inh3_g,ibin)
  b = - kg(inh3_g,ibin)*gas(inh3_g)   &
       + kg(ihcl_g,ibin)*gas(ihcl_g)
  c = -(kg(ihcl_g,ibin)*Keq_sg(2))

  sfc_a(inh3_g) = quadratic(a,b,c)
  sfc_a(ihcl_g) = Keq_sg(2)/sfc_a(inh3_g)

  df_gas_s(ihcl_g,ibin) = gas(ihcl_g) - sfc_a(ihcl_g)
  df_gas_s(inh3_g,ibin) = gas(inh3_g) - sfc_a(inh3_g)

  flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*df_gas_s(inh3_g,ibin)
  flux_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin) + flux_s(inh3_g,ibin)

  phi_volatile_s(inh3_g,ibin) = phi_nh4cl_s

  if(flux_s(ihcl_g,ibin) .gt. 0.0)then
     df_gas_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin)/kg(ihcl_g,ibin)        
     phi_volatile_s(ihcl_g,ibin) = phi_nh4cl_s
  else
     sfc_a(ihcl_g)  = gas(ihcl_g) + aer(icl_a,jsolid,ibin)
     df_gas_s(ihcl_g,ibin) = -aer(icl_a,jsolid,ibin)
     phi_volatile_s(ihcl_g,ibin)=df_gas_s(ihcl_g,ibin)/sfc_a(ihcl_g)  
  endif

  integrate(inh3_g,jsolid,ibin) = mYES
  integrate(ihcl_g,jsolid,ibin) = mYES  

  ieqblm_ASTEM = mNO

  return
end subroutine ASTEM_flux_dry_case3b







subroutine ASTEM_flux_dry_case4(ibin, phi_nh4no3_s,phi_nh4cl_s,ieqblm_ASTEM,     &
     sfc_a,df_gas_s,flux_s,phi_volatile_s,integrate,kg,gas,electrolyte,epercent, &
     Keq_sg,aer)       

  use module_data_mosaic_aero, only: r8, nbin_a_max,                             &
       ngas_aerchtot, ngas_volatile, nelectrolyte,                               &
       nsalt,jsolid,nrxn_aer_sg,naer,                                            &
       rtol_eqb_ASTEM,ptol_mol_ASTEM,                                            &
       jnh4no3,jnh4cl,ihno3_g,inh3_g,ihcl_g

  use module_mosaic_ext, only: quadratic,degas_solid_nh4no3,degas_solid_nh4cl
  

  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(out) :: phi_nh4no3_s,phi_nh4cl_s
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
  
  integer iactive_nh4no3, iactive_nh4cl, iactive, js
  real(r8) :: a, b, c
  real(r8) :: sum_dum
  
  



  
  sum_dum = 0.0
  do js = 1, nsalt
     sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
  enddo

  if(sum_dum .eq. 0.)sum_dum = 1.0

  epercent(jnh4no3,jsolid,ibin) = 100.*electrolyte(jnh4no3,jsolid,ibin)/sum_dum
  epercent(jnh4cl, jsolid,ibin) = 100.*electrolyte(jnh4cl, jsolid,ibin)/sum_dum
  


  
  
  iactive_nh4no3 = 1
  iactive_nh4cl  = 2


  
  phi_nh4no3_s = (gas(inh3_g)*gas(ihno3_g) - Keq_sg(1))/   &
       max(gas(inh3_g)*gas(ihno3_g),Keq_sg(1))
  phi_nh4cl_s  = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/   &
       max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))


  
  

  
  if( abs(phi_nh4no3_s) .lt. rtol_eqb_ASTEM )then
     iactive_nh4no3 = 0
  elseif(gas(inh3_g)*gas(ihno3_g) .lt. Keq_sg(1) .and.   &
       epercent(jnh4no3,jsolid,ibin) .le. ptol_mol_ASTEM)then
     iactive_nh4no3 = 0
     if(epercent(jnh4no3,jsolid,ibin) .gt. 0.0)then
        call degas_solid_nh4no3(ibin,aer,gas,electrolyte,Keq_sg)
     endif
  endif

  
  if( abs(phi_nh4cl_s) .lt. rtol_eqb_ASTEM )then
     iactive_nh4cl = 0
  elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and.   &
       epercent(jnh4cl, jsolid,ibin) .le. ptol_mol_ASTEM)then
     iactive_nh4cl = 0
     if(epercent(jnh4cl, jsolid,ibin) .gt. 0.0)then
        call degas_solid_nh4cl(ibin,aer,gas,electrolyte,Keq_sg)
     endif
  endif


  iactive = iactive_nh4no3 + iactive_nh4cl

  
  if(iactive .eq. 0)return


  goto (1,2,3),iactive

  
  
1 call ASTEM_flux_dry_case4a(ibin,phi_nh4no3_s,ieqblm_ASTEM,sfc_a,df_gas_s,      &
       flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg)
  return


  
  
2 call ASTEM_flux_dry_case4b(ibin,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s,flux_s,&
       phi_volatile_s,integrate,kg,gas,Keq_sg)
  return


  
  
3 call ASTEM_flux_dry_case4ab(ibin,phi_nh4no3_s,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,  &
       df_gas_s,flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg)


  return
end subroutine ASTEM_flux_dry_case4






subroutine ASTEM_flux_dry_case4a(ibin, phi_nh4no3_s,ieqblm_ASTEM,sfc_a,df_gas_s, &
     flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg) 

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       jsolid,mYES,mNO, nrxn_aer_sg,                                               &
       ihno3_g,inh3_g

  use module_mosaic_ext, only: quadratic


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(in) :: phi_nh4no3_s
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  
  real(r8) :: a, b, c
  
  



  a =   kg(inh3_g,ibin)
  b = - kg(inh3_g,ibin)*gas(inh3_g)   &
       + kg(ihno3_g,ibin)*gas(ihno3_g)
  c = -(kg(ihno3_g,ibin)*Keq_sg(1))

  sfc_a(inh3_g)  = quadratic(a,b,c)
  sfc_a(ihno3_g) = Keq_sg(1)/sfc_a(inh3_g)

  integrate(ihno3_g,jsolid,ibin) = mYES
  integrate(inh3_g,jsolid,ibin)  = mYES

  df_gas_s(ihno3_g,ibin)=gas(ihno3_g)-sfc_a(ihno3_g)
  df_gas_s(inh3_g,ibin) =gas(inh3_g) -sfc_a(inh3_g)

  phi_volatile_s(ihno3_g,ibin)= phi_nh4no3_s
  phi_volatile_s(inh3_g,ibin) = phi_nh4no3_s

  flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
  flux_s(inh3_g,ibin)  = flux_s(ihno3_g,ibin)

  ieqblm_ASTEM = mNO

  return
end subroutine ASTEM_flux_dry_case4a






subroutine ASTEM_flux_dry_case4b(ibin,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s,  &
     flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg) 

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       mYES,jsolid, mNO,nrxn_aer_sg,                                               &
       inh3_g,ihcl_g

  use module_mosaic_ext, only: quadratic
  

  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(in) :: phi_nh4cl_s
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  
  real(r8) :: a, b, c
  
  


  a =   kg(inh3_g,ibin)
  b = - kg(inh3_g,ibin)*gas(inh3_g)   &
       + kg(ihcl_g,ibin)*gas(ihcl_g)
  c = -(kg(ihcl_g,ibin)*Keq_sg(2))

  sfc_a(inh3_g) = quadratic(a,b,c)
  sfc_a(ihcl_g) = Keq_sg(2) /sfc_a(inh3_g)

  integrate(ihcl_g,jsolid,ibin) = mYES
  integrate(inh3_g,jsolid,ibin) = mYES

  df_gas_s(ihcl_g,ibin) = gas(ihcl_g)-sfc_a(ihcl_g)
  df_gas_s(inh3_g,ibin) = gas(inh3_g)-sfc_a(inh3_g)

  phi_volatile_s(ihcl_g,ibin) = phi_nh4cl_s
  phi_volatile_s(inh3_g,ibin) = phi_nh4cl_s

  flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
  flux_s(inh3_g,ibin) = flux_s(ihcl_g,ibin)

  ieqblm_ASTEM = mNO

  return
end subroutine ASTEM_flux_dry_case4b







subroutine ASTEM_flux_dry_case4ab(ibin, phi_nh4no3_s, phi_nh4cl_s,ieqblm_ASTEM,  &
     sfc_a,df_gas_s,flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg)    

  use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       mNO,nrxn_aer_sg,&
       ihcl_g,ihno3_g,inh3_g

  use module_mosaic_ext, only: quadratic


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(in) :: phi_nh4no3_s,phi_nh4cl_s
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,flux_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  
  real(r8) :: a,b,c,flux_nh3_est, flux_nh3_max, ratio_flux
  
  

  call ASTEM_flux_dry_case4a(ibin,phi_nh4no3_s,ieqblm_ASTEM,sfc_a,df_gas_s,      &
       flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg)
  call ASTEM_flux_dry_case4b(ibin,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s,       &
       flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg)


  

  flux_nh3_est = flux_s(ihno3_g,ibin)+flux_s(ihcl_g,ibin)
  flux_nh3_max = kg(inh3_g,ibin)*gas(inh3_g)


  if(flux_nh3_est .le. flux_nh3_max)then

     flux_s(inh3_g,ibin) = flux_nh3_est                 
     sfc_a(inh3_g)       = gas(inh3_g) -                           &  
          flux_s(inh3_g,ibin)/kg(inh3_g,ibin)
     phi_volatile_s(inh3_g,ibin) = max(abs(phi_nh4no3_s),   &
          abs(phi_nh4cl_s))

  else                  

     ratio_flux          = flux_nh3_max/flux_nh3_est
     flux_s(inh3_g,ibin) = flux_nh3_max
     flux_s(ihno3_g,ibin)= flux_s(ihno3_g,ibin)*ratio_flux
     flux_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin) *ratio_flux

     sfc_a(inh3_g) = 0.0
     sfc_a(ihno3_g)= gas(ihno3_g) -                        &  
          flux_s(ihno3_g,ibin)/kg(ihno3_g,ibin)
     sfc_a(ihcl_g) = gas(ihcl_g)  -                        &  
          flux_s(ihcl_g,ibin)/kg(ihcl_g,ibin)

     df_gas_s(inh3_g,ibin) =gas(inh3_g) -sfc_a(inh3_g)
     df_gas_s(ihno3_g,ibin)=gas(ihno3_g)-sfc_a(ihno3_g)
     df_gas_s(ihcl_g,ibin) =gas(ihcl_g) -sfc_a(ihcl_g)

     phi_volatile_s(inh3_g,ibin) = max(abs(phi_nh4no3_s),   &
          abs(phi_nh4cl_s))


  endif

  ieqblm_ASTEM = mNO

  return
end subroutine ASTEM_flux_dry_case4ab














subroutine ASTEM_flux_mix(ibin, phi_nh4no3_s, phi_nh4cl_s, ieqblm_ASTEM, idry_case3a, &
     sfc_a, df_gas_s, df_gas_l, jaerosolstate, flux_s, Heff, phi_volatile_s,            &
     phi_volatile_l, integrate, jphase, aer, kg, gas, jhyst_leg, electrolyte, epercent,   &
     kel, activity, mc, delta_nh3_max, delta_hno3_max, delta_hcl_max, Keq_nh4cl,        &
     Keq_nh4no3, num_a, electrolyte_sum, mass_dry_a, mass_soluble_a, water_a, aH2O,     &
     kelvin_nh4no3, kelvin_nh4cl, ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma,      &
     nc_Mc, xeq_c, mw_electrolyte, Kp_nh4cl, Kp_nh4no3, Keq_ll, Keq_gl, Keq_sg, MW_c,     &
     MW_a, total_species, tot_cl_in, molality0, kappa_nonelectro, mosaic_vars_aa )

  use module_data_mosaic_aero, only: r8, nbin_a_max,                             &
       ngas_aerchtot, ngas_volatile, nelectrolyte,                               &
       Ncation, naer, jliquid, nsalt, jsolid, mNO, mYES, jtotal, Nanion, nrxn_aer_gl,      &
       nrxn_aer_ll, nrxn_aer_sg,                                                   &
       jcaco3, jcacl2, jnacl, ihno3_g, jnh4cl, ihcl_g, inh3_g, jnh4no3, ja_no3, jc_h,      &
       ja_cl, jhcl, icl_a, inh4_a, ino3_a, jhno3, mosaic_vars_aa_type

  use module_mosaic_ext,  only: compute_activities, ions_to_electrolytes,           &
       absorb_tiny_nh4cl, degas_tiny_nh4cl, absorb_tiny_nh4no3, degas_tiny_nh4no3,   &
       absorb_tiny_hno3, absorb_tiny_hcl


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(nbin_a_max) :: idry_case3a,jaerosolstate      
  integer, intent(inout), dimension(nbin_a_max) :: jphase,jhyst_leg
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(out) :: phi_nh4no3_s, phi_nh4cl_s
  real(r8), intent(in) :: aH2O
  real(r8), intent(inout) :: Keq_nh4cl,Keq_nh4no3,kelvin_nh4no3,Kp_nh4cl
  real(r8), intent(inout) :: kelvin_nh4cl,Kp_nh4no3
  real(r8), intent(in), dimension(Ncation) :: zc,MW_c
  real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
  real(r8), intent(in), dimension(Nanion)  :: za,MW_a
  real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
  real(r8), intent(inout), dimension(nbin_a_max) :: delta_nh3_max,delta_hno3_max
  real(r8), intent(inout), dimension(nbin_a_max) :: delta_hcl_max,mass_soluble_a
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a,num_a,mass_dry_a
  real(r8), intent(inout), dimension(nbin_a_max) :: gam_ratio
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a, total_species
  real(r8), intent(inout) :: tot_cl_in
  real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_l
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: flux_s,Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg, kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max)  :: activity,gam
  real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) ::mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  real(r8), intent(inout), dimension(3,nbin_a_max) :: electrolyte_sum
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
  real(r8), intent(in), dimension(naer) :: kappa_nonelectro

  type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

  
  character(len=500) :: tmp_str
  integer iv, iadjust, iadjust_intermed, js
  real(r8) :: XT,g_nh3_hno3,g_nh3_hcl,a_nh4_no3,a_nh4_cl,a_no3,a_cl,prod_nh4no3
  real(r8) :: volatile_cl,sum_dum,prod_nh4cl


  call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,   &
       nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)    
  call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,electrolyte,   &
       activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,ma,gam,log_gamZ, &
       gam_ratio,Keq_ll,molality0,kappa_nonelectro)

  if(water_a(ibin) .eq. 0.0)then
     write(tmp_str,*)'Water is zero in liquid phase'
     call mosaic_warn_mess(trim(adjustl(tmp_str)))    
     write(tmp_str,*)'Stopping in ASTEM_flux_wet'
     call mosaic_warn_mess(trim(adjustl(tmp_str)))    
     mosaic_vars_aa%zero_water_flag = .true.
  endif


  
  sum_dum = 0.0
  do js = 1, nsalt
     sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
  enddo

  if(sum_dum .eq. 0.)sum_dum = 1.0

  epercent(jcaco3,jsolid,ibin) = 100.*electrolyte(jcaco3,jsolid,ibin)/sum_dum
  



  
  

  if(epercent(jcaco3,jsolid,ibin) .gt. 0.0)then
     jphase(ibin) = jliquid
     call ASTEM_flux_wet_case1(ibin,ieqblm_ASTEM,sfc_a,df_gas_s,flux_s,          &
          phi_volatile_s,integrate,jphase,kg,gas,mc,Keq_ll)
     return
  endif

  
  


  if(XT.lt.2.0 .and. XT.ge.0.)then   
     jphase(ibin) = jliquid
     call ASTEM_flux_wet_case2(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,            &
          phi_volatile_l,integrate,gas,kel,mc,water_a,ma,gam,gam_ratio,Keq_ll,   &
          Keq_gl)
     return
  endif

  
  

  volatile_cl  = electrolyte(jnacl,jsolid,ibin) +   &
       electrolyte(jcacl2,jsolid,ibin)


  if(volatile_cl .gt. 0.0 .and. gas(ihno3_g).gt. 0.0 )then

     call ASTEM_flux_dry_case3a(ibin,ieqblm_ASTEM,idry_case3a,sfc_a,df_gas_s,    &
          flux_s,phi_volatile_s,integrate,aer,kg,gas)

     prod_nh4cl = max( (gas(inh3_g)*gas(ihcl_g)-Keq_sg(2)), 0.0d0) +   &
          electrolyte(jnh4cl, jsolid,ibin)

     if(prod_nh4cl .gt. 0.0)then
        call ASTEM_flux_dry_case3b(ibin,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s, &
             flux_s,phi_volatile_s,integrate,aer,kg,gas,electrolyte,epercent,    &
             Keq_sg)
     endif

     jphase(ibin) = jsolid

     return
  endif

  
  

  if( electrolyte(jnh4no3,jsolid,ibin).gt.0. .and.   &
       electrolyte(jnh4cl,jsolid,ibin) .gt.0. )then
     jphase(ibin) = jsolid
     call ASTEM_flux_dry_case4(ibin,phi_nh4no3_s,phi_nh4cl_s,ieqblm_ASTEM,sfc_a, &
          df_gas_s,flux_s,phi_volatile_s,integrate,kg,gas,electrolyte,epercent,  &
          Keq_sg,aer)

     if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
             (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
     elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
             (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
     else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
     endif

     return

  elseif( electrolyte(jnh4no3,jsolid,ibin).gt.0. )then
     
     g_nh3_hcl= gas(inh3_g)*gas(ihcl_g)
     a_nh4_cl = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)

     iadjust = mNO              
     if(g_nh3_hcl .gt. 0.0 .and. a_nh4_cl .eq. 0.0)then
        call absorb_tiny_nh4cl(ibin,aer,gas,electrolyte,delta_nh3_max,           &
             delta_hcl_max,electrolyte_sum)
        iadjust = mYES
     elseif(g_nh3_hcl .eq. 0.0 .and. a_nh4_cl .gt. 0.0)then
        call degas_tiny_nh4cl(ibin,aer,gas,electrolyte)
        iadjust = mYES
     endif

     if(iadjust .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,   &
             na_Ma,nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)      
        call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,         &
             electrolyte,activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,    &
             aH2O,ma,gam,log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)                      
     endif

     call ASTEM_flux_mix_case4a(ibin,phi_nh4no3_s,ieqblm_ASTEM,sfc_a,df_gas_s,   &
          df_gas_l,flux_s,Heff,phi_volatile_s,phi_volatile_l,integrate,jphase,kg,&
          gas,electrolyte,epercent,kel,activity,mc,Keq_nh4cl,water_a,            &
          kelvin_nh4cl,ma,gam,gam_ratio,Kp_nh4cl,Keq_ll,Keq_gl,Keq_sg,aer)      
     jphase(ibin) = jtotal
     return

  elseif( electrolyte(jnh4cl,jsolid,ibin).gt.0.)then
     
     g_nh3_hno3= gas(inh3_g)*gas(ihno3_g)
     a_nh4_no3 = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)

     iadjust = mNO              
     if(g_nh3_hno3 .gt. 0.0 .and. a_nh4_no3 .eq. 0.0)then
        call absorb_tiny_nh4no3(ibin,aer,gas,electrolyte,delta_nh3_max,          &
             delta_hno3_max,electrolyte_sum)
        iadjust = mYES
     elseif(g_nh3_hno3 .eq. 0.0 .and. a_nh4_no3 .gt. 0.0)then
        call degas_tiny_nh4no3(ibin,aer,gas,electrolyte)
        iadjust = mYES
     endif

     if(iadjust .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,   &
             na_Ma,nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)      
        call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,         &
             electrolyte,activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,    &
             aH2O,ma,gam,log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)                      
     endif

     kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
     Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3        

     call ASTEM_flux_mix_case4b(ibin,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s,    &
          df_gas_l,flux_s,Heff,phi_volatile_s,phi_volatile_l,integrate,jphase,kg,&
          gas,electrolyte,epercent,kel,activity,mc,Keq_nh4no3,water_a,           &
          kelvin_nh4no3,ma,gam,gam_ratio,Keq_ll,Keq_gl,Kp_nh4no3,Keq_sg,aer)    
     jphase(ibin) = jtotal
     return
  endif


  

  if( (gas(inh3_g)+aer(inh4_a,jliquid,ibin)) .lt. 1.e-25)goto 10  

  
  
  

  iadjust = mNO         
  iadjust_intermed = mNO        

  
  g_nh3_hno3 = gas(inh3_g)*gas(ihno3_g)
  a_nh4_no3  = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)

  if(g_nh3_hno3 .gt. 0. .and. a_nh4_no3 .eq. 0.)then
     call absorb_tiny_nh4no3(ibin,aer,gas,electrolyte,delta_nh3_max,             &
          delta_hno3_max,electrolyte_sum)
     iadjust = mYES
     iadjust_intermed = mYES
  endif

  if(iadjust_intermed .eq. mYES)then
     call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,&
          nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)         
     iadjust_intermed = mNO     
  endif

  
  g_nh3_hcl = gas(inh3_g)*gas(ihcl_g)
  a_nh4_cl  = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)

  if(g_nh3_hcl .gt. 0. .and. a_nh4_cl .eq. 0.)then
     call absorb_tiny_nh4cl(ibin,aer,gas,electrolyte,delta_nh3_max,delta_hcl_max,&
          electrolyte_sum)
     iadjust = mYES
     iadjust_intermed = mYES
  endif

  if(iadjust_intermed .eq. mYES)then
     call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,&
          nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)         
  endif

  if(iadjust .eq. mYES)then
     call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,electrolyte,&
          activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,ma,gam,       &
          log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)                 
  endif


  

  
  kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
  Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3   

  kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
  Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl       

  call ASTEM_flux_wet_case3(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,phi_volatile_l,&
       integrate,kg,gas,kel,mc,Keq_nh4cl,Keq_nh4no3,water_a,ma,gam,gam_ratio,    &
       Keq_ll,Keq_gl,aer,total_species,tot_cl_in,activity,electrolyte)
  jphase(ibin) = jliquid

  return


  
  
  

10 iadjust = mNO                
  iadjust_intermed = mNO        

  
  if(gas(ihno3_g).gt.0. .and. aer(ino3_a,jliquid,ibin).eq.0. .and.   &
       aer(icl_a,jliquid,ibin) .gt. 0.0)then
     call absorb_tiny_hno3(ibin,aer,gas,delta_hno3_max)        
     iadjust = mYES
     iadjust_intermed = mYES
  endif

  if(iadjust_intermed .eq. mYES)then
     call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,&
          nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)         
     iadjust_intermed = mNO     
  endif

  
  if(gas(ihcl_g).gt.0. .and. aer(icl_a,jliquid,ibin) .eq. 0. .and.   &
       aer(ino3_a,jliquid,ibin) .gt. 0.0)then
     call absorb_tiny_hcl(ibin,aer,gas,delta_hcl_max)                 
     iadjust = mYES
     iadjust_intermed = mYES
  endif

  if(iadjust_intermed .eq. mYES)then
     call ions_to_electrolytes(jliquid,ibin,XT,aer,electrolyte,zc,za,xeq_a,na_Ma,&
          nc_Mc,xeq_c,mw_electrolyte,MW_c,MW_a)         
  endif

  if(iadjust .eq. mYES)then
     call compute_activities(ibin,jaerosolstate,jphase,aer,jhyst_leg,            &
          electrolyte,activity,mc,num_a,mass_dry_a,mass_soluble_a,water_a,aH2O,  &
          ma,gam,log_gamZ,gam_ratio,Keq_ll,molality0,kappa_nonelectro)                 
  endif

  

  call ASTEM_flux_wet_case4(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,phi_volatile_l,&
       integrate,kg,gas,kel,mc,water_a,ma,gam,Keq_ll,Keq_gl)
  jphase(ibin) = jliquid

  return
end subroutine ASTEM_flux_mix








subroutine ASTEM_flux_mix_case4a(ibin, phi_nh4no3_s,ieqblm_ASTEM,sfc_a,df_gas_s, &
     df_gas_l,flux_s,Heff,phi_volatile_s,phi_volatile_l,integrate,jphase,kg,gas, &
     electrolyte,epercent,kel,activity,mc,Keq_nh4cl,water_a,kelvin_nh4cl,ma,gam, &
     gam_ratio,Kp_nh4cl,Keq_ll,Keq_gl,Keq_sg,aer)   

  use module_data_mosaic_aero, only: r8, nbin_a_max,                             &
       ngas_aerchtot, ngas_volatile, nelectrolyte,                               &
       Ncation,mYES,jsolid,mNO,jliquid,jtotal,Nanion,nrxn_aer_gl,nrxn_aer_ll,    &
       nrxn_aer_sg,naer,                                                         &
       rtol_eqb_ASTEM,ptol_mol_ASTEM,                                            &
       jnh4no3,ihno3_g,inh3_g,ihcl_g,jnh4cl,ja_no3,jhno3,jc_h,ja_cl,jhcl

  use module_mosaic_ext, only: degas_solid_nh4no3


  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(nbin_a_max) :: jphase
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(inout) :: Keq_nh4cl,kelvin_nh4cl,Kp_nh4cl
  real(r8), intent(out) :: phi_nh4no3_s
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a,gam_ratio
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s,df_gas_l
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: flux_s,Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg, kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max)  :: activity,gam
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) ::mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
  
  integer iactive_nh4no3, iactive_nh4cl, js
  real(r8) :: sum_dum


  
  iactive_nh4no3 = mYES
  iactive_nh4cl  = mYES


  
  sum_dum = 0.0
  do js = 1, nelectrolyte
     sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
  enddo

  if(sum_dum .eq. 0.)sum_dum = 1.0

  epercent(jnh4no3,jsolid,ibin) = 100.*electrolyte(jnh4no3,jsolid,ibin)/sum_dum
  



  
  phi_nh4no3_s = (gas(inh3_g)*gas(ihno3_g) - Keq_sg(1))/   &
       max(gas(inh3_g)*gas(ihno3_g),Keq_sg(1))

  
  kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
  Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl       


  
  
  
  if( abs(phi_nh4no3_s) .le. rtol_eqb_ASTEM )then
     iactive_nh4no3 = mNO
  elseif(gas(inh3_g)*gas(ihno3_g) .lt. Keq_sg(1) .and.   &
       epercent(jnh4no3,jsolid,ibin) .le. ptol_mol_ASTEM)then
     iactive_nh4no3 = mNO
     if(epercent(jnh4no3,jsolid,ibin) .gt. 0.0)then
        call degas_solid_nh4no3(ibin,aer,gas,electrolyte,Keq_sg)
     endif
  endif

  
  if( gas(inh3_g)*gas(ihcl_g).eq.0. .or. Keq_nh4cl.eq.0. )then
     iactive_nh4cl = mNO
  endif


  
  if(iactive_nh4no3 .eq. mYES)then

     jphase(ibin) = jsolid
     call ASTEM_flux_dry_case4a(ibin,phi_nh4no3_s,ieqblm_ASTEM,sfc_a,df_gas_s,   &
          flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg)      

     if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
             (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
     elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
             (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
     else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
     endif

  endif


  if(iactive_nh4cl .eq. mYES)then

     jphase(ibin) = jliquid
     call ASTEM_flux_wet_case3b(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,           &
          phi_volatile_l,integrate,kg,gas,kel,mc,Keq_nh4cl,water_a,ma,gam,       &
          gam_ratio,Keq_ll,Keq_gl)        

     if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
             (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
     else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
     endif

  endif


  if(iactive_nh4cl .eq. mYES .and. iactive_nh4no3 .eq. mYES)then
     jphase(ibin) = jtotal
  endif



  return
end subroutine ASTEM_flux_mix_case4a






subroutine ASTEM_flux_mix_case4b(ibin,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s,   &
     df_gas_l,flux_s,Heff,phi_volatile_s,phi_volatile_l,integrate,jphase,kg,gas, &
     electrolyte,epercent,kel,activity,mc,Keq_nh4no3,water_a,kelvin_nh4no3,ma,   &
     gam,gam_ratio,Keq_ll,Keq_gl,Kp_nh4no3,Keq_sg,aer) 

  use module_data_mosaic_aero, only: r8, nbin_a_max,                             &
       ngas_aerchtot, ngas_volatile, nelectrolyte,                               &
       Ncation,mYES,nsalt,jsolid,mNO,jliquid,jtotal,Nanion,nrxn_aer_gl,naer,     &
       nrxn_aer_ll,nrxn_aer_sg,                                                  &
       rtol_eqb_ASTEM,ptol_mol_ASTEM,                                            &
       jnh4cl,ihcl_g,inh3_g,ihno3_g,jnh4no3,ja_cl,jhcl,jc_h,ja_no3,jhno3

  use module_mosaic_ext, only: degas_solid_nh4cl

  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_ASTEM
  integer, intent(inout), dimension(nbin_a_max) :: jphase
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate

  real(r8), intent(inout) :: Keq_nh4no3,kelvin_nh4no3,Kp_nh4no3
  real(r8), intent(out) :: phi_nh4cl_s
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a,gam_ratio
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
  real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: df_gas_l
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: flux_s,Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg,kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max)  :: activity,gam
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
  
  integer iactive_nh4no3, iactive_nh4cl, js
  real(r8) :: sum_dum


  
  iactive_nh4cl  = mYES
  iactive_nh4no3 = mYES


  
  sum_dum = 0.0
  do js = 1, nsalt
     sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
  enddo

  if(sum_dum .eq. 0.)sum_dum = 1.0

  epercent(jnh4cl,jsolid,ibin) = 100.*electrolyte(jnh4cl,jsolid,ibin)/sum_dum
  


  
  phi_nh4cl_s  = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/   &
       max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))

  
  kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
  Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3   


  
  
  
  if( abs(phi_nh4cl_s) .le. rtol_eqb_ASTEM )then
     iactive_nh4cl = mNO
  elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and.   &
       epercent(jnh4cl,jsolid,ibin) .le. ptol_mol_ASTEM)then
     iactive_nh4cl = mNO
     if(epercent(jnh4cl,jsolid,ibin) .gt. 0.0)then
        call degas_solid_nh4cl(ibin,aer,gas,electrolyte,Keq_sg)
     endif
  endif

  
  if( gas(inh3_g)*gas(ihno3_g).eq.0. .or. Keq_nh4no3.eq.0. )then
     iactive_nh4no3 = mNO
  endif


  
  if(iactive_nh4cl .eq. mYES)then

     jphase(ibin) = jsolid
     call ASTEM_flux_dry_case4b(ibin,phi_nh4cl_s,ieqblm_ASTEM,sfc_a,df_gas_s,    &
          flux_s,phi_volatile_s,integrate,kg,gas,Keq_sg)    

     if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
             (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
     elseif(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
             (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
     else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
     endif

  endif


  if(iactive_nh4no3 .eq. mYES)then

     jphase(ibin) = jliquid
     call ASTEM_flux_wet_case3a(ibin,ieqblm_ASTEM,sfc_a,df_gas_l,Heff,           &
          phi_volatile_l,integrate,kg,gas,kel,mc,Keq_nh4no3,water_a,ma,gam,      &
          gam_ratio,Keq_ll,Keq_gl)    

     if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
             (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
     else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
     endif

  endif


  if(iactive_nh4cl .eq. mYES .and. iactive_nh4no3 .eq. mYES)then
     jphase(ibin) = jtotal
  endif



  return
end subroutine ASTEM_flux_mix_case4b










subroutine ASTEM_non_volatiles( dtchem,  jaerosolstate, jphase, &
     aer, kg, gas, gas_avg, gas_netprod_otrproc,                                      &
     jhyst_leg, electrolyte, epercent, kel, activity, mc, delta_nh3_max,                &
     delta_hno3_max, delta_hcl_max, num_a, mass_wet_a, mass_dry_a, mass_soluble_a,     &
     vol_dry_a, vol_wet_a, water_a, water_a_hyst, water_a_up, aH2O_a, total_species,    &
     tot_cl_in,                                                                   &
     aH2O, ma, gam, log_gamZ, zc, za, gam_ratio,  &
     xeq_a, na_Ma, nc_Mc, xeq_c, a_zsr, mw_electrolyte, partial_molar_vol, sigma_soln,   &
     T_K, RH_pc, mw_aer_mac, dens_aer_mac, sigma_water, Keq_ll, Keq_sl, MW_a, MW_c,       &
     growth_factor, MDRH, MDRH_T, molality0, rtol_mesa, jsalt_present, jsalt_index,     &
     jsulf_poor, jsulf_rich, phi_salt_old,                                   &
     kappa_nonelectro, mosaic_vars_aa ) 
  
  use module_data_mosaic_aero,  only: r8, nbin_a_max, nbin_a,   &
       ngas_aerchtot, ngas_volatile, nelectrolyte,    &
       Ncation, naer, no_aerosol, jtotal, mNO, mYES, Nanion, nrxn_aer_ll, nrxn_aer_sl,    &
       nsalt, MDRH_T_NUM, jsulf_poor_NUM, jsulf_rich_NUM,                            &
       ih2so4_g, imsa_g, inh3_g, ihno3_g, ihcl_g, iso4_a, imsa_a, jcaco3, jcano3, jnano3,  &
       jcacl2, jnacl, inh4_a, mosaic_vars_aa_type

  use module_mosaic_ext, only: aerosol_phase_state,conform_electrolytes


  
  integer, intent(in), dimension(nsalt) :: jsalt_index
  integer, intent(in), dimension(jsulf_poor_NUM) :: jsulf_poor
  integer, intent(in), dimension(jsulf_rich_NUM) :: jsulf_rich

  real(r8), intent(in) :: dtchem
  real(r8), intent(in) :: aH2O,T_K,RH_pc,rtol_mesa
  real(r8), intent(in), dimension(naer) :: mw_aer_mac,dens_aer_mac
  real(r8), intent(in), dimension(Ncation) :: zc,MW_c
  real(r8), intent(in), dimension(Nanion)  :: za,MW_a
  real(r8), intent(in),    dimension(ngas_aerchtot) :: gas_netprod_otrproc 
  real(r8), intent(in), dimension(ngas_aerchtot) :: partial_molar_vol
  real(r8), intent(in), dimension(nelectrolyte) :: mw_electrolyte
  real(r8), intent(in), dimension (6,nelectrolyte) :: a_zsr
  real(r8), intent(in), dimension(naer) :: kappa_nonelectro

  
  integer, intent(inout), dimension(nsalt) :: jsalt_present
  integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase,jhyst_leg

  real(r8), intent(inout) :: sigma_water
  real(r8), intent(inout), dimension(Ncation) :: nc_Mc,xeq_c
  real(r8), intent(inout), dimension(Nanion)  :: xeq_a,na_Ma
  real(r8), intent(inout), dimension(nbin_a_max) :: num_a,mass_soluble_a,vol_dry_a
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a,delta_nh3_max
  real(r8), intent(inout), dimension(nbin_a_max) :: delta_hno3_max,delta_hcl_max
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a_hyst,water_a_up,aH2O_a
  real(r8), intent(inout), dimension(nbin_a_max) :: mass_wet_a,mass_dry_a
  real(r8), intent(inout), dimension(nbin_a_max) :: growth_factor,MDRH
  real(r8), intent(inout), dimension(nbin_a_max) :: vol_wet_a,gam_ratio,sigma_soln
  real(r8), intent(inout), dimension(ngas_volatile) :: total_species
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas_avg  

  type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

  
  
  real(r8), intent(inout) :: tot_cl_in
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 
  real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
  real(r8), intent(inout), dimension(nrxn_aer_sl) :: Keq_sl
  real(r8), intent(inout), dimension(MDRH_T_NUM) :: MDRH_T
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg,kel
  real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
  real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
  real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
  real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: epercent
  real(r8), intent(inout), dimension(nsalt) :: phi_salt_old

  
  integer ibin,iupdate_phase_state
  real(r8) :: decay_h2so4,decay_msa,delta_h2so4,delta_tmsa,delta_nh3,delta_hno3                
  real(r8) :: delta_hcl,XT,sumkg_h2so4,sumkg_msa,sumkg_nh3,sumkg_hno3,sumkg_hcl
  real(r8) :: tmp_kxt, tmp_kxt2, tmp_pok, tmp_pxt, tmp_q1, tmp_q3, tmp_q4
  real(r8), dimension(nbin_a) :: delta_so4,delta_msa,delta_nh4
  real(r8), dimension(nbin_a) :: new_so4a, old_so4a 

  sumkg_h2so4 = 0.0
  sumkg_msa   = 0.0
  sumkg_nh3   = 0.0
  sumkg_hno3  = 0.0
  sumkg_hcl   = 0.0
  do ibin = 1, nbin_a
     sumkg_h2so4 = sumkg_h2so4 + kg(ih2so4_g,ibin)
     sumkg_msa   = sumkg_msa   + kg(imsa_g,ibin)
     sumkg_nh3   = sumkg_nh3   + kg(inh3_g,ibin)
     sumkg_hno3  = sumkg_hno3  + kg(ihno3_g,ibin)
     sumkg_hcl   = sumkg_hcl   + kg(ihcl_g,ibin)
  enddo



  
  
  tmp_q1 = gas(ih2so4_g)
  tmp_pxt = max( gas_netprod_otrproc(ih2so4_g)*dtchem, 0.0_r8 )
  tmp_kxt = sumkg_h2so4*dtchem
  old_so4a(1:nbin_a) = aer(iso4_a,jtotal,1:nbin_a)   
  if ( (tmp_q1+tmp_pxt > 1.e-14_r8) .and. &
       (tmp_kxt >= 1.0e-20_r8) ) then






     
     
     
     
     if (tmp_kxt > 0.001_r8) then
        
        tmp_pok = tmp_pxt/tmp_kxt
        tmp_q3 = (tmp_q1 - tmp_pok)*exp(-tmp_kxt) + tmp_pok
        tmp_q4 = (tmp_q1 - tmp_pok)*(1.0_r8 - exp(-tmp_kxt))/tmp_kxt + tmp_pok
     else
        
        tmp_kxt2 = tmp_kxt*tmp_kxt
        tmp_q3 = tmp_q1 *(1.0_r8 - tmp_kxt        + tmp_kxt2*0.5_r8) &
               + tmp_pxt*(1.0_r8 - tmp_kxt*0.5_r8 + tmp_kxt2/6.0_r8)
        tmp_q4 = tmp_q1 *(1.0_r8 - tmp_kxt*0.5_r8 + tmp_kxt2/6.0_r8) &
               + tmp_pxt*(0.5_r8 - tmp_kxt/6.0_r8 + tmp_kxt2/24.0_r8)
     end if
     gas(ih2so4_g) = tmp_q3
     gas_avg(ih2so4_g) = tmp_q4
     delta_h2so4 = (tmp_q1 + tmp_pxt) - tmp_q3   


     
     do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .ne. no_aerosol)then
           delta_so4(ibin) = delta_h2so4*kg(ih2so4_g,ibin)/sumkg_h2so4
           aer(iso4_a,jtotal,ibin) = aer(iso4_a,jtotal,ibin) +   &
                delta_so4(ibin)
        endif
     enddo

  else
     
     
     
     gas(ih2so4_g) = tmp_q1 + tmp_pxt
     gas_avg(ih2so4_g) = tmp_q1 + tmp_pxt*0.5_r8
     delta_h2so4 = 0.0
     do ibin = 1, nbin_a
        delta_so4(ibin) = 0.0
     enddo

  endif

  new_so4a(1:nbin_a) = aer(iso4_a,jtotal,1:nbin_a)   
  
  



  
  if(gas(imsa_g) .gt. 1.e-14)then

     
     decay_msa   = exp(-sumkg_msa*dtchem)
     delta_tmsa  = gas(imsa_g)*(1.0 - decay_msa)
     gas(imsa_g) = gas(imsa_g)*decay_msa

     
     do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .ne. no_aerosol)then
           delta_msa(ibin) = delta_tmsa*kg(imsa_g,ibin)/sumkg_msa
           aer(imsa_a,jtotal,ibin) = aer(imsa_a,jtotal,ibin) +   &
                delta_msa(ibin)
        endif
     enddo

  else

     delta_tmsa = 0.0
     do ibin = 1, nbin_a
        delta_msa(ibin) = 0.0
     enddo

  endif
  
  



  
  delta_nh3 = gas(inh3_g) *(1.0 - exp(-sumkg_nh3*dtchem))
  delta_hno3= gas(ihno3_g)*(1.0 - exp(-sumkg_hno3*dtchem))
  delta_hcl = gas(ihcl_g) *(1.0 - exp(-sumkg_hcl*dtchem))

  
  do ibin = 1, nbin_a
     if(jaerosolstate(ibin) .ne. no_aerosol)then
        delta_nh3_max(ibin) = delta_nh3*kg(inh3_g,ibin)/sumkg_nh3
        delta_hno3_max(ibin)= delta_hno3*kg(ihno3_g,ibin)/sumkg_hno3
        delta_hcl_max(ibin) = delta_hcl*kg(ihcl_g,ibin)/sumkg_hcl
     endif
  enddo


  if(delta_h2so4 .eq. 0.0 .and. delta_tmsa .eq. 0.0)then
     iupdate_phase_state = mNO
     goto 100
  endif


  
  do ibin = 1, nbin_a

     if(electrolyte(jnacl,jtotal,ibin)  .eq. 0.0 .and.   &
          electrolyte(jcacl2,jtotal,ibin) .eq. 0.0 .and.   &
          electrolyte(jnano3,jtotal,ibin) .eq. 0.0 .and.   &
          electrolyte(jcano3,jtotal,ibin) .eq. 0.0 .and.   &
          electrolyte(jcaco3,jtotal,ibin) .eq. 0.0 .and.   &
          jaerosolstate(ibin) .ne. no_aerosol)then

        delta_nh4(ibin) = min( (2.*delta_so4(ibin)+delta_msa(ibin)),   &
             delta_nh3_max(ibin) )

        aer(inh4_a,jtotal,ibin) = aer(inh4_a,jtotal,ibin) +        &  
             delta_nh4(ibin)

        gas(inh3_g) = gas(inh3_g) - delta_nh4(ibin)             

     else

        delta_nh4(ibin)     = 0.0

     endif

  enddo

  iupdate_phase_state = mYES


  
100 if(iupdate_phase_state .eq. mYES)then
     do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .ne. no_aerosol)then
           call conform_electrolytes(jtotal, ibin, XT, aer, gas, electrolyte,          &
                total_species, tot_cl_in)
           call aerosol_phase_state( ibin, jaerosolstate,                     &
                jphase, aer, jhyst_leg, electrolyte, epercent, kel, activity, mc, num_a, &
                mass_wet_a, mass_dry_a, mass_soluble_a, vol_dry_a, vol_wet_a, water_a,   &
                water_a_hyst, water_a_up, aH2O_a, aH2O,                  &
                ma, gam, log_gamZ, zc, za, gam_ratio, xeq_a, na_Ma, nc_Mc,   &
                xeq_c,        mw_electrolyte, partial_molar_vol, sigma_soln, T_K,        & 
                RH_pc, mw_aer_mac, dens_aer_mac, sigma_water, Keq_ll, Keq_sl, MW_a,      &
                MW_c, growth_factor, MDRH, MDRH_T, molality0, rtol_mesa, jsalt_present,  &
                jsalt_index, jsulf_poor, jsulf_rich, phi_salt_old,            &
                kappa_nonelectro, mosaic_vars_aa )
        endif
     enddo
  endif

  return
end subroutine ASTEM_non_volatiles













subroutine ASTEM_secondary_organics(dtchem, jaerosolstate,sfc_a,Heff,            &
     phi_volatile_l,integrate,aer,kg,gas,sat_soa,total_species)
  
  use module_data_mosaic_aero, only: nbin_a_max, nbin_a, naer, no_aerosol,   &
       ngas_aerchtot, ngas_volatile, jtotal,mYES,                            &
       isoa_first
  
  
  
  integer, intent(in), dimension(nbin_a_max) :: jaerosolstate  
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(in) :: dtchem
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a, sat_soa, total_species
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: Heff
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  
  integer ibin, iv, jp,ieqblm, nsteps_max,ieqblm_soa,isteps_SOA
  parameter(nsteps_max = 400)
  real(r8) :: dtmax, t_new, t_old, t_out
  real(r8) :: sum1, sum2
  
 
  
  t_old = 0.0
  t_out = dtchem
  isteps_SOA = 0

  
  do iv = isoa_first, ngas_volatile
     total_species(iv) = gas(iv)
     do ibin = 1, nbin_a
        if (jaerosolstate(ibin) .eq. no_aerosol) cycle
        total_species(iv) = total_species(iv) + aer(iv,jtotal,ibin)
     enddo
  enddo
  
  
  
  
10 isteps_SOA = isteps_SOA + 1
  
  
  ieqblm_soa = mYES			
  
  do 501 ibin = 1, nbin_a
     if (jaerosolstate(ibin) .eq. no_aerosol) goto 501
     
     call ASTEM_flux_soa(ibin,sfc_a,Heff,integrate,aer,gas,sat_soa,ieqblm_soa)
     
501 continue
  if(ieqblm_soa .eq. mYES)goto 30 
  
  
  
  

11 call ASTEM_dtmax_soa(dtchem, dtmax, phi_volatile_l,integrate,kg)
  t_new = t_old + dtmax	
  if(t_new .gt. t_out)then	
     dtmax = t_out - t_old
     t_new = t_out*1.01
  endif
  
  
  
  
  
  
  
  jp = jtotal
  
  do 20 iv = isoa_first, ngas_volatile
     
     sum1 = 0.0
     sum2 = 0.0
     
     do 21 ibin = 1, nbin_a
        if(jaerosolstate(ibin) .eq. no_aerosol)goto 21
        
        sum1 = sum1 + aer(iv,jp,ibin)/   &
             (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin)*integrate(iv,jp,ibin))
        sum2 = sum2 + kg(iv,ibin)*integrate(iv,jp,ibin)/   &
             (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin)*integrate(iv,jp,ibin))
        
21   continue
        
     
     gas(iv) = (total_species(iv) - sum1)/   &
          (1. + dtmax*sum2)
     
     
     do 22 ibin = 1, nbin_a
        if (jaerosolstate(ibin) .eq. no_aerosol) goto 22 
        
        if(integrate(iv,jp,ibin) .eq. mYES)then
           aer(iv,jp,ibin) =   &
                (aer(iv,jp,ibin) + dtmax*kg(iv,ibin)*gas(iv))/   &
                (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin))
        endif
        
22   continue
        
20 continue
  
  

        
  
  
  
  


  
  t_old = t_new
  
  if(t_new .lt. 0.9999*t_out) goto 10
  
  
  
30 continue
  
  
  return
end subroutine ASTEM_secondary_organics









subroutine ASTEM_flux_soa(ibin,sfc_a,Heff,integrate,aer,gas,sat_soa,ieqblm_soa)		

  use module_data_mosaic_aero, only: r8, ngas_aerchtot, ngas_volatile, &
       nbin_a_max,naer,mNO,jtotal, rtol_eqb_ASTEM,                     &
       ioc_a,isoa_first
  
  
  
  integer, intent(in) :: ibin
  integer, intent(inout) :: ieqblm_soa
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
  real(r8), intent(inout), dimension(ngas_volatile) :: sfc_a, sat_soa
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: Heff
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  
  integer iv, jp
  real(r8) :: dum, sum_dum, sum_soa, small_oc
  real(r8), dimension(ngas_volatile,nbin_a_max) :: df_gas_o,flux_o,phi_volatile_o

  small_oc  = 1.e-15		
  
  
  
  do iv = isoa_first, ngas_volatile
     sfc_a(iv)               = gas(iv)
     df_gas_o(iv,ibin)       = 0.0
     flux_o(iv,ibin)         = 0.0
     phi_volatile_o(iv,ibin) = 0.0
  enddo
  
  
  jp = jtotal
  
  
  sum_soa = 0.0
  do iv = isoa_first, ngas_volatile
     sum_soa = sum_soa + aer(iv,jp,ibin)
  enddo
  sum_soa = sum_soa + aer(ioc_a,jp,ibin)/200.  
  
  
  
  if(aer(ioc_a,jp,ibin) .eq. 0.0)then
     sum_dum = 0.0
     do iv = isoa_first, ngas_volatile
        sum_dum = sum_dum + (gas(iv)+aer(iv,jp,ibin))/sat_soa(iv)
     enddo
     
     if(sum_dum .le. 1.0)then	
        do iv = isoa_first, ngas_volatile
           gas(iv)         = gas(iv) + aer(iv,jp,ibin)
           aer(iv,jp,ibin) = 0.0
           integrate(iv,jp,ibin) = 0.0
        enddo
        return
     endif
     
     sum_soa = max(sum_soa, 1.d-10)
     
  endif
  
  
  
  
  
  do iv = isoa_first, ngas_volatile
     
     Heff(iv,ibin) = sat_soa(iv)/sum_soa
     sfc_a(iv) = aer(iv,jp,ibin)*Heff(iv,ibin)		
     df_gas_o(iv,ibin) = gas(iv) - sfc_a(iv)
     
     dum = max(sfc_a(iv),gas(iv))
     if(dum .gt. 0.0)then
        phi_volatile_o(iv,ibin) = df_gas_o(iv,ibin)/dum
     else
        phi_volatile_o(iv,ibin) = 0.0
     endif
     
     
     if(abs(phi_volatile_o(iv,ibin)) .le. rtol_eqb_ASTEM)then
        integrate(iv,jp,ibin) = 0.0
     else
        integrate(iv,jp,ibin) = 1.0
        ieqblm_soa = mNO
     endif
     
  enddo
  
  
  return
end subroutine ASTEM_flux_soa









subroutine ASTEM_dtmax_soa(dtchem, dtmax, phi_volatile_l,integrate,kg)		

  use module_data_mosaic_aero, only: r8, ngas_aerchtot, ngas_volatile, &
       nbin_a_max,jtotal,mYES, alpha_astem,nbin_a,                     &
       isoa_first
  
  
  
  integer, intent(inout), dimension(ngas_volatile,3,nbin_a_max) :: integrate
  
  real(r8), intent(in)  :: dtchem
  real(r8), intent(out) :: dtmax
  real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
  real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg      
  
  
  character(len=500) :: tmp_str
  integer ibin, iv, jp
  real(r8) :: h_gas, h_gas_i(ngas_volatile), h_sub_max,   &
       sum_kg_phi
  
  
  h_sub_max = dtchem/6.	
  
  jp = jtotal
  
  
  

  h_gas = 2.e16
  
  do 6 iv = isoa_first, ngas_volatile
     
     h_gas_i(iv) = 1.e16
     sum_kg_phi = 0.0
     
     do ibin = 1, nbin_a
        if(integrate(iv,jtotal,ibin) .eq. mYES)then
           sum_kg_phi = sum_kg_phi +   &
                abs(phi_volatile_l(iv,ibin))*kg(iv,ibin)
        endif
     enddo
     
     if(sum_kg_phi .gt. 0.0)then
        h_gas_i(iv) = alpha_astem/sum_kg_phi
        h_gas       = min(h_gas, h_gas_i(iv))
     endif
     
6 continue


  dtmax = min(h_gas, h_sub_max)
  
  
  if(dtmax .le. 1.0e-10)then
     write(tmp_str,*)' SOA dtmax = ', dtmax
     call mosaic_warn_mess(trim(adjustl(tmp_str))) 
  endif
  
  
  return
end subroutine ASTEM_dtmax_soa



end module module_mosaic_astem
