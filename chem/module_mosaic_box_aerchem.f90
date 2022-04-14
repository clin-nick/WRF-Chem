module module_mosaic_box_aerchem

use module_data_mosaic_kind, only: r8

implicit none

contains
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  

  
  
  
  
  
  

  subroutine mosaic_box_aerchemistry(        aH2O,               T_K,            &
       P_atm,                  RH_pc,        dtchem,                             &
       mcall_load_mosaic_parameters,         mcall_print_aer_in, sigmag_a,       &
       kappa_nonelectro,                                                         &
       jaerosolstate,          aer,                                              &
       num_a,                  water_a,      gas,                                &
       gas_avg,                gas_netprod_otrproc,              Dp_dry_a,       &
       dp_wet_a,               jhyst_leg,                                        &
       mosaic_vars_aa,                                                           &
       mass_dry_a_bgn,         mass_dry_a,                                       &
       dens_dry_a_bgn,         dens_dry_a,   water_a_hyst,       aH2O_a,         &
       uptkrate_h2so4,         gam_ratio,    jaerosolstate_bgn                   )

    use module_data_mosaic_aero, only:                                             &
         nbin_a_max, ngas_aerchtot, ngas_volatile, naer, nsalt,                    &
         Nanion, Ncation, nrxn_aer_sl, nrxn_aer_ll, nrxn_aer_gl, nrxn_aer_sg,      &
         MDRH_T_NUM, nelectrolyte,                                                 &
         jsalt_index, jsulf_poor, jsulf_rich, rtol_mesa, dens_aer_mac,             &
         mw_aer_mac, zc, MW_c, za, MW_a, mw_comp_a, dens_comp_a, b_zsr,aw_min,     &
         mw_electrolyte, partial_molar_vol, a_zsr, d_mdrh, b_mtem, ref_index_a,    &
         Nmax_mesa, nmax_ASTEM, mosaic_vars_aa_type
         
    implicit none

    
    integer, intent(in) :: mcall_load_mosaic_parameters, mcall_print_aer_in

    real(r8), intent(in) :: aH2O
    real(r8), intent(in) :: T_K, P_atm, RH_pc
    real(r8), intent(in) :: dtchem

    real(r8), intent(in), dimension(nbin_a_max)        :: sigmag_a
    real(r8), intent(in), dimension(naer)              :: kappa_nonelectro
                
    
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate
    integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg

    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nbin_a_max)        :: num_a
    real(r8), intent(inout), dimension(nbin_a_max)        :: water_a
    real(r8), intent(inout), dimension(ngas_aerchtot)     :: gas
    real(r8), intent(inout), dimension(ngas_aerchtot)     :: gas_avg  
    real(r8), intent(in),    dimension(ngas_aerchtot)     :: gas_netprod_otrproc
              
              
              
              
              
              
    real(r8), intent(inout), dimension(nbin_a_max)        :: Dp_dry_a, dp_wet_a

    
    
    type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

    
    integer, intent(out), dimension(nbin_a_max) :: jaerosolstate_bgn

    real(r8), intent(out), dimension(nbin_a_max) :: mass_dry_a_bgn
    real(r8), intent(out), dimension(nbin_a_max) :: mass_dry_a
    real(r8), intent(out), dimension(nbin_a_max) :: dens_dry_a_bgn
    real(r8), intent(out), dimension(nbin_a_max) :: dens_dry_a
    real(r8), intent(out), dimension(nbin_a_max) :: water_a_hyst
    real(r8), intent(out), dimension(nbin_a_max) :: aH2O_a
    real(r8), intent(out), dimension(nbin_a_max) :: gam_ratio
    real(r8), intent(out)                        :: uptkrate_h2so4  

    
    integer :: iprint_input, irepeat_mosaic
    integer :: mcall_print_aer
    integer, dimension(nbin_a_max) :: jphase

    integer :: iaer 

    real(r8) :: sigma_water,Kp_nh4cl
    real(r8) :: Kp_nh4no3,Kp_nh3
    real(r8) :: tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in, tot_ca_in

    real(r8), dimension(nbin_a_max) :: mass_soluble_a
    real(r8), dimension(ngas_volatile) :: sat_soa,total_species
    real(r8), dimension(nrxn_aer_sl) :: Keq_sl
    real(r8), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), dimension(nrxn_aer_gl) :: Keq_gl
    real(r8), dimension(nrxn_aer_sg) :: Keq_sg
    real(r8), dimension(MDRH_T_NUM) :: MDRH_T
    real(r8), dimension(nelectrolyte,nbin_a_max) :: molality0 
    real(r8), dimension(ngas_volatile,nbin_a_max) :: flux_s,flux_l
    real(r8), dimension(ngas_volatile,nbin_a_max) :: volatile_s
    real(r8), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
    real(r8), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
    real(r8), dimension(ngas_aerchtot,nbin_a_max) :: kg
    real(r8), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
    real(r8), dimension(Ncation,nbin_a_max) :: mc
    real(r8), dimension(Nanion,nbin_a_max) :: ma
    real(r8), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    
    call update_thermodynamic_constants(  aH2O,    T_K,                                & 
         sat_soa,    aH2O_a,   log_gamZ,  Keq_sl,  sigma_water,  Kp_nh4cl,             & 
         Kp_nh4no3,  Kp_nh3,   Keq_ll,    Keq_gl,  Keq_sg,       MDRH_T,               &
         molality0                                                                     )





    do irepeat_mosaic = 1, 1
       mcall_print_aer = mcall_print_aer_in
       if (irepeat_mosaic > 1) mcall_print_aer = 0

       call initialize_mosaic_variables(                                                & 
            jaerosolstate, flux_s, flux_l, volatile_s, phi_volatile_s, phi_volatile_l,  & 
            jphase, kg, electrolyte, activity, mc, mass_dry_a, mass_soluble_a,          &
            dens_dry_a, ma, gam, gam_ratio                                              )

       mosaic_vars_aa%isteps_astem = 0
       mosaic_vars_aa%isteps_astem_max = 0
       mosaic_vars_aa%jastem_call = 0
       mosaic_vars_aa%jmesa_call = 0
       mosaic_vars_aa%jmesa_fail = 0
       mosaic_vars_aa%niter_mesa_max = 0
       mosaic_vars_aa%nmax_astem = nmax_astem
       mosaic_vars_aa%nmax_mesa = nmax_mesa
       mosaic_vars_aa%cumul_steps_astem = 0.0_r8
       mosaic_vars_aa%niter_mesa = 0.0_r8
       uptkrate_h2so4 = 0.0_r8

       call overall_massbal_in( aer, gas, gas_netprod_otrproc, dtchem,                  & 
            total_species, tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in,    & 
            tot_ca_in )

       call MOSAIC_dynamic_solver(      mcall_print_aer,     dtchem,                    & 
            aH2O,           T_K,        RH_pc,               P_atm,                     &
            irepeat_mosaic, tot_cl_in,  sigmag_a,            kappa_nonelectro,          &
            jaerosolstate,  flux_s,     flux_l,              volatile_s,                & 
            phi_volatile_s, phi_volatile_l,                  jphase,           aer,     &
            kg,             gas,        gas_avg,             gas_netprod_otrproc,       &
            jhyst_leg,      electrolyte,                     activity,                  &
            mc,             sat_soa,    num_a,               Dp_dry_a,         Dp_wet_a,&
            mass_dry_a,     mass_soluble_a,                  dens_dry_a,       water_a, &
            gam,            log_gamZ,   gam_ratio,           Keq_ll,           Keq_gl,  &
            Keq_sg,         Keq_sl,     Kp_nh4cl,            Kp_nh4no3,        ma,      &
            sigma_water,    MDRH_T,     molality0,                                      &
            total_species,  aH2O_a,     uptkrate_h2so4,                                 &
            mosaic_vars_aa,                                                             &
            iprint_input,                                                               & 
            mass_dry_a_bgn, dens_dry_a_bgn,                                             &
            water_a_hyst,   jaerosolstate_bgn                                           )

       if (mosaic_vars_aa%f_mos_fail > 0) then
          return
       endif
       
       
       call overall_massbal_out( iprint_input, 0, mosaic_vars_aa%isteps_ASTEM, aer, gas, &
          tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in, tot_ca_in )

    enddo




    return
  end subroutine mosaic_box_aerchemistry



  
  
  
  
  
  
  subroutine MOSAIC_aerosol_water_only(      aH2O,         T_K,                  &
       P_atm,                  RH_pc,        dtchem,                             &
       kappa_nonelectro,                                                         &
       jaerosolstate,          jhyst_leg,                                        &
       aer,                    num_a,        water_a,      gas,                  &
       Dp_dry_a,               Dp_wet_a,                                         &
       mosaic_vars_aa,                                                           &
       mass_dry_a,             dens_dry_a                                        )

    use module_data_mosaic_aero,  only:                                                 &
         nbin_a_max, ngas_aerchtot, ngas_volatile, nelectrolyte,                        &
         Nanion, Ncation, naer, nbin_a, no_aerosol, jtotal,                             &
         jsalt_index, jsulf_poor, jsulf_rich,                                           &
         jhyst_lo, jhyst_up, mhyst_method, mhyst_uporlo_jhyst,                          &
         mhyst_uporlo_waterhyst, mhyst_force_lo, mhyst_force_up,                        &
         mSECTIONAL, mSIZE_FRAMEWORK, MDRH_T_NUM,                                       &
         nrxn_aer_gl, nrxn_aer_ll, nrxn_aer_sg, nrxn_aer_sl, nsalt,                     &
         zc, za, a_zsr, b_zsr, aw_min,                                                  &
         mw_electrolyte, partial_molar_vol,                                             &
         dens_aer_mac, mw_aer_mac, dens_comp_a, mw_comp_a, ref_index_a, MW_a, MW_c,     &
         density_max_allow, density_min_allow,                                          &
         rtol_mesa, nmax_astem, nmax_mesa, mosaic_vars_aa_type

    use module_data_mosaic_asecthp, only: isize_of_ibin, itype_of_ibin,dcen_sect       
    
    use module_mosaic_ext,        only: aerosol_water_up, aerosol_phase_state, &
                                        calc_dry_n_wet_aerosol_props, conform_electrolytes
    
    implicit none
    
    
    real(r8), intent(in) :: dtchem
    real(r8), intent(in) :: aH2O, T_K, RH_pc, P_atm

    real(r8), intent(in), dimension(naer)       :: kappa_nonelectro

    
    integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate

    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nbin_a_max) :: Dp_dry_a, Dp_wet_a
    real(r8), intent(inout), dimension(nbin_a_max) :: dens_dry_a
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a

    type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa
   
    

        

    
    character(len=256) :: errmsg
    integer :: ibin, isize, itype, iv
    integer :: irepeat_mosaic
    integer, dimension(nbin_a_max) :: jphase
    integer, dimension(nbin_a_max) :: jaerosolstate_bgn
    integer, dimension(ngas_volatile,3,nbin_a_max) :: integrate
    integer, dimension(nsalt) :: jsalt_present

    real(r8) :: Keq_nh4cl
    real(r8) :: Kp_nh3, Kp_nh4cl, Kp_nh4no3
    real(r8) :: sigma_water
    real(r8) :: tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in, tot_ca_in
    real(r8) :: XT

    real(r8), dimension(MDRH_T_NUM) :: MDRH_T
    real(r8), dimension(Nanion ) :: na_Ma, xeq_a
    real(r8), dimension(Ncation) :: nc_Mc, xeq_c
    real(r8), dimension(Nanion, nbin_a_max) :: ma
    real(r8), dimension(Ncation,nbin_a_max) :: mc
    real(r8), dimension(nsalt) :: phi_salt_old

    real(r8), dimension(nbin_a_max) :: aH2O_a
    real(r8), dimension(nbin_a_max) :: area_dry_a, area_wet_a
    real(r8), dimension(nbin_a_max) :: delta_hcl_max, delta_nh3_max, delta_hno3_max
    real(r8), dimension(nbin_a_max) :: dens_dry_a_bgn, dens_wet_a
    real(r8), dimension(nbin_a_max) :: dp_core_a
    real(r8), dimension(nbin_a_max) :: gam_ratio
    real(r8), dimension(nbin_a_max) :: growth_factor
    real(r8), dimension(nbin_a_max) :: mass_dry_a_bgn, mass_soluble_a, mass_wet_a
    real(r8), dimension(nbin_a_max) :: MDRH
    real(r8), dimension(nbin_a_max) :: sigma_soln
    real(r8), dimension(nbin_a_max) :: vol_dry_a, vol_wet_a
    real(r8), dimension(nbin_a_max) :: water_a_hyst, water_a_up

    real(r8), dimension(ngas_aerchtot) :: gas_netprod_otrproc
    real(r8), dimension(ngas_volatile) :: sfc_a
    real(r8), dimension(ngas_volatile) :: sat_soa, total_species
    real(r8), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l, phi_volatile_s
    real(r8), dimension(ngas_volatile,nbin_a_max) :: volatile_s
    real(r8), dimension(ngas_volatile,nbin_a_max) :: flux_s, flux_l
    real(r8), dimension(ngas_aerchtot,nbin_a_max) :: kel
    real(r8), dimension(ngas_aerchtot,nbin_a_max) :: kg

    real(r8), dimension(nelectrolyte,nbin_a_max) :: activity, gam
    real(r8), dimension(nelectrolyte,nbin_a_max) :: molality0 
    real(r8), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), dimension(nelectrolyte,3,nbin_a_max) :: epercent
    real(r8), dimension(nelectrolyte,nelectrolyte) :: log_gamZ

    real(r8), dimension(nrxn_aer_sl) :: Keq_sl
    real(r8), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), dimension(nrxn_aer_gl) :: Keq_gl
    real(r8), dimension(nrxn_aer_sg) :: Keq_sg

    complex, dimension(nbin_a_max) :: ri_shell_a,ri_avg_a,ri_core_a


    gas_netprod_otrproc(1:ngas_volatile) = 0.0_r8

    call update_thermodynamic_constants(  aH2O,    T_K,                                & 
         sat_soa,    aH2O_a,   log_gamZ,  Keq_sl,  sigma_water,  Kp_nh4cl,             & 
         Kp_nh4no3,  Kp_nh3,   Keq_ll,    Keq_gl,  Keq_sg,       MDRH_T,               &
         molality0                                                                     )

    irepeat_mosaic = 1

    call initialize_mosaic_variables(                                                & 
         jaerosolstate, flux_s, flux_l, volatile_s, phi_volatile_s, phi_volatile_l,  & 
         jphase, kg, electrolyte, activity, mc, mass_dry_a, mass_soluble_a,          &
         dens_dry_a, ma, gam, gam_ratio                                              )

    mosaic_vars_aa%isteps_astem = 0
    mosaic_vars_aa%isteps_astem_max = 0
    mosaic_vars_aa%jastem_call = 0
    mosaic_vars_aa%jmesa_call = 0
    mosaic_vars_aa%jmesa_fail = 0
    mosaic_vars_aa%niter_mesa_max = 0
    mosaic_vars_aa%nmax_astem = nmax_astem
    mosaic_vars_aa%nmax_mesa = nmax_mesa
    mosaic_vars_aa%cumul_steps_astem = 0.0_r8
    mosaic_vars_aa%niter_mesa = 0.0_r8

    call overall_massbal_in( aer, gas, gas_netprod_otrproc, dtchem,                  & 
         total_species, tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in,    & 
         tot_ca_in )


    vol_dry_a = 0.0_r8

    do ibin = 1, nbin_a
       call check_aerosol_mass(ibin, jaerosolstate,jphase,aer,num_a, mass_dry_a)
       jaerosolstate_bgn(ibin) = jaerosolstate(ibin)
       
       if(jaerosolstate(ibin) .ne. no_aerosol) then
          
          call conform_electrolytes(jtotal,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)        
          call check_aerosol_mass(ibin,jaerosolstate,jphase,aer,num_a, mass_dry_a) 
          
          jaerosolstate_bgn(ibin) = jaerosolstate(ibin)
          if(jaerosolstate(ibin) .ne. no_aerosol)then 
             
             call conform_aerosol_number(ibin,jaerosolstate,aer,num_a,vol_dry_a,Dp_dry_a) 
             
             
             
             
             if (mosaic_vars_aa%it_mosaic == 1) then
                if (mhyst_method == mhyst_uporlo_waterhyst) then
                   if(jhyst_leg(ibin) == jhyst_lo)then
                      water_a_hyst(ibin) = 0.0
                   else
                      water_a_up(ibin)   = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)	
                      water_a_hyst(ibin) = water_a_up(ibin)
                   endif
                else if (mhyst_method == mhyst_force_lo) then
                   jhyst_leg(ibin) = jhyst_lo
                   water_a_hyst(ibin) = 0.0
                else if (mhyst_method == mhyst_force_up) then
                   jhyst_leg(ibin)    = jhyst_up
                   water_a_up(ibin)   = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)	
                   water_a_hyst(ibin) = water_a_up(ibin)
                end if
             end if
             
          endif
       endif
       if (irepeat_mosaic == 1) then
          mass_dry_a_bgn(ibin) = mass_dry_a(ibin)
          if ( (jaerosolstate(ibin) .eq. no_aerosol) .or.   &
               (min(mass_dry_a(ibin),vol_dry_a(ibin)) .le. 1.0e-35) ) then
             call calc_aerosol_dry_density( ibin,aer,dens_dry_a)
             dens_dry_a_bgn(ibin) = dens_dry_a(ibin)
          else
             dens_dry_a_bgn(ibin) = mass_dry_a(ibin)/vol_dry_a(ibin)
          end if
          dens_dry_a_bgn(ibin) = max( density_min_allow, &
               min( density_max_allow, dens_dry_a_bgn(ibin) ) )
       end if
       
       if (jaerosolstate(ibin) .eq. no_aerosol) then
          if (msize_framework == msectional) then
             isize = isize_of_ibin(ibin)
             itype = itype_of_ibin(ibin)
             Dp_dry_a(ibin) = dcen_sect(isize,itype)
             Dp_wet_a(ibin) = Dp_dry_a(ibin)
          end if
       end if
       
    enddo
    


  
  do ibin = 1, nbin_a

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

        call calc_dry_n_wet_aerosol_props(                                &
           ibin, jaerosolstate, aer, electrolyte, water_a, num_a,         &  
           dens_comp_a, mw_comp_a, dens_aer_mac, mw_aer_mac, ref_index_a, &  
           Dp_dry_a, Dp_wet_a, dp_core_a,                                 &  
           area_dry_a, area_wet_a, mass_dry_a, mass_wet_a,                &  
           vol_dry_a, vol_wet_a, dens_dry_a, dens_wet_a,                  &  
           ri_shell_a, ri_core_a, ri_avg_a                                )  
     endif
  enddo
    
    
    do ibin = 1, nbin_a
       if(jaerosolstate(ibin).ne.no_aerosol) then 
          

          if (mhyst_method == mhyst_uporlo_jhyst .or. & 
              mhyst_method == mhyst_uporlo_waterhyst) then
             
             if(jhyst_leg(ibin) == jhyst_lo)then
                water_a_hyst(ibin) = 0.0
             else
                water_a_up(ibin)   = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)   
                water_a_hyst(ibin) = water_a_up(ibin)
             endif










          else if (mhyst_method == mhyst_force_lo) then
             jhyst_leg(ibin) = jhyst_lo
             water_a_hyst(ibin) = 0.0
          else if (mhyst_method == mhyst_force_up) then
             jhyst_leg(ibin) = jhyst_up
             water_a_up(ibin)   = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)   
             water_a_hyst(ibin) = water_a_up(ibin)
             
          else
             write(errmsg,*) '*** MOSAIC_aerosol_water - bad mhyst_method =', mhyst_method
             call wrf_error_fatal3("<stdin>",510,&
trim(adjustl(errmsg)))
          endif
          
       endif
       if ( (jaerosolstate(ibin) .eq. no_aerosol) .or.   &
            (min(mass_dry_a(ibin),vol_dry_a(ibin)) .le. 1.0e-35) ) then
          call calc_aerosol_dry_density( ibin,aer,dens_dry_a)
       end if
       dens_dry_a(ibin) = max( density_min_allow, &
            min( density_max_allow, dens_dry_a(ibin) ) )
       
    enddo

    return
  end subroutine MOSAIC_aerosol_water_only



  
  
  
  
  
  
 
  subroutine MOSAIC_dynamic_solver( mcall_print_aer,    dtchem,                    & 
       aH2O,           T_K,        RH_pc,               P_atm,                     &
       irepeat_mosaic, tot_cl_in,  sigmag_a,                                       &
       kappa_nonelectro,                                                           &
       jaerosolstate,  flux_s,     flux_l,              volatile_s,                & 
       phi_volatile_s, phi_volatile_l,                  jphase,           aer,     &
       kg,             gas,        gas_avg,             gas_netprod_otrproc,       &
       jhyst_leg,      electrolyte,                     activity,                  &
       mc,             sat_soa,    num_a,               Dp_dry_a,         Dp_wet_a,&
       mass_dry_a,     mass_soluble_a,                  dens_dry_a,       water_a, &
       gam,            log_gamZ,   gam_ratio,           Keq_ll,           Keq_gl,  &
       Keq_sg,         Keq_sl,     Kp_nh4cl,            Kp_nh4no3,        ma,      &
       sigma_water,    MDRH_T,     molality0,                                      &
       total_species,  aH2O_a,     uptkrate_h2so4,                                 &
       mosaic_vars_aa,                                                             &
       iprint_input,                                                               & 
       mass_dry_a_bgn, dens_dry_a_bgn,                                             &
       water_a_hyst,   jaerosolstate_bgn                                           )
       
    use module_data_mosaic_aero,  only: nbin_a_max, ngas_aerchtot, ngas_volatile,        &
         nelectrolyte,                                                                   &
         Ncation, naer, no_aerosol, jtotal, mhyst_uporlo_waterhyst, jhyst_lo,            &
         density_max_allow, density_min_allow, mSECTIONAL, mON, mASTEM, mLSODE,          &
         mhyst_uporlo_jhyst, jhyst_up, Nanion, nrxn_aer_gl, nrxn_aer_ll,                &
         nrxn_aer_sg, nrxn_aer_sl, nsalt, MDRH_T_NUM,  mhyst_force_lo,  mhyst_force_up,  &
         nbin_a, mSIZE_FRAMEWORK, mGAS_AER_XFER, mDYNAMIC_SOLVER, mhyst_method,         &
         zc, za, a_zsr, mw_electrolyte, partial_molar_vol, dens_aer_mac,      &
         mw_aer_mac,  dens_comp_a, mw_comp_a, ref_index_a, MW_a, MW_c, rtol_mesa,         &
         jsalt_index, jsulf_poor, jsulf_rich,                                          &
         iso4_a,                                                                       & 
         mosaic_vars_aa_type

    use module_data_mosaic_asecthp, only: isize_of_ibin,itype_of_ibin,dcen_sect       
    
    use module_mosaic_astem,      only: ASTEM
    
    use module_mosaic_ext,        only: aerosol_water_up,calc_dry_n_wet_aerosol_props,&
                                        conform_electrolytes, dumpxx

    use module_mosaic_lsode,      only: mosaic_lsode
    
    implicit none
    
    
    integer, intent(in) :: mcall_print_aer
    integer, intent(in) :: irepeat_mosaic
    
    real(r8), intent(in) :: dtchem
    real(r8), intent(in) :: aH2O, T_K, RH_pc, P_atm

    real(r8), intent(in), dimension(nbin_a_max) :: sigmag_a
    real(r8), intent(in), dimension(naer)       :: kappa_nonelectro

    
    real(r8), intent(inout) :: Kp_nh4cl
    real(r8), intent(inout) :: Kp_nh4no3,sigma_water
    real(r8), intent(inout) :: tot_cl_in

    real(r8), intent(inout), dimension(MDRH_T_NUM) :: MDRH_T

    integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase

    real(r8), intent(inout), dimension(nbin_a_max) :: num_a, Dp_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: Dp_wet_a, gam_ratio
    real(r8), intent(inout), dimension(nbin_a_max) :: aH2O_a
    
    real(r8), intent(inout), dimension(nbin_a_max) :: dens_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_soluble_a
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(ngas_volatile) :: sat_soa, total_species
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas_avg  
    real(r8), intent(in),    dimension(ngas_aerchtot) :: gas_netprod_otrproc
              
              
              
              
              
              

    real(r8), intent(inout), dimension(nrxn_aer_ll) :: Keq_ll
    real(r8), intent(inout), dimension(nrxn_aer_gl) :: Keq_gl
    real(r8), intent(inout), dimension(nrxn_aer_sg) :: Keq_sg
    real(r8), intent(inout), dimension(nrxn_aer_sl) :: Keq_sl

    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: molality0 

    real(r8), intent(inout), dimension(Ncation,nbin_a_max) :: mc
    real(r8), intent(inout), dimension(Nanion,nbin_a_max) :: ma
    
    real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: flux_s,flux_l
    real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: volatile_s
    real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
    real(r8), intent(inout), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
    real(r8), intent(inout), dimension(ngas_aerchtot,nbin_a_max) :: kg
    real(r8), intent(inout), dimension(nelectrolyte,nbin_a_max) :: activity,gam
    real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ

    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    real(r8), intent(inout) :: uptkrate_h2so4  

    type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa
   
    
    integer, intent(out) :: iprint_input
    integer, intent(out), dimension(nbin_a_max) :: jaerosolstate_bgn

    real(r8), intent(out), dimension(nbin_a_max) :: water_a_hyst
    real(r8), intent(out), dimension(nbin_a_max) :: mass_dry_a_bgn,dens_dry_a_bgn
        
    
    character(len=256) :: errmsg
    integer ibin, isize, itype, iv

    real(r8) :: XT
    

    real(r8), dimension(nbin_a_max) :: area_dry_a,water_a_up
    real(r8), dimension(nbin_a_max) :: area_wet_a,mass_wet_a,vol_wet_a,dens_wet_a
    real(r8), dimension(nbin_a_max) :: vol_dry_a
    real(r8), dimension(nbin_a_max) :: dp_core_a

    real(r8), dimension(nelectrolyte,3,nbin_a_max) :: epercent

    complex, dimension(nbin_a_max) :: ri_shell_a,ri_avg_a,ri_core_a
    

    vol_dry_a = 0.0_r8

    
    mosaic_vars_aa%jASTEM_fail = 0
    mosaic_vars_aa%jASTEM_call       = 0
    mosaic_vars_aa%isteps_ASTEM      = 0
    mosaic_vars_aa%isteps_ASTEM_max  = 0
    mosaic_vars_aa%niter_MESA        = 0.0_r8
    mosaic_vars_aa%cumul_steps_ASTEM = 0.0_r8


    call dumpxx( 'aa', dtchem, t_k, p_atm, ah2o, &
         jaerosolstate, jphase, jhyst_leg, &
         aer, gas, num_a, water_a, water_a_hyst, dp_dry_a, &
         mosaic_vars_aa )


    do ibin = 1, nbin_a
       call check_aerosol_mass(ibin, jaerosolstate,jphase,aer,num_a, mass_dry_a)
       jaerosolstate_bgn(ibin) = jaerosolstate(ibin)
       
       if(jaerosolstate(ibin) .ne. no_aerosol) then
          
          
          
          call conform_electrolytes(jtotal,ibin,XT,aer,gas,electrolyte,total_species,tot_cl_in)        
          call check_aerosol_mass(ibin,jaerosolstate,jphase,aer,num_a, mass_dry_a) 
          
          jaerosolstate_bgn(ibin) = jaerosolstate(ibin)
          if(jaerosolstate(ibin) .ne. no_aerosol)then 
             
             
             call conform_aerosol_number(ibin,jaerosolstate,aer,num_a,vol_dry_a,Dp_dry_a) 
             
             
             
             
             if (mosaic_vars_aa%it_mosaic == 1) then
                if (mhyst_method == mhyst_uporlo_waterhyst) then
                   if(jhyst_leg(ibin) == jhyst_lo)then
                      water_a_hyst(ibin) = 0.0
                   else
                      water_a_up(ibin)   = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)	
                      water_a_hyst(ibin) = water_a_up(ibin)
                   endif
                else if (mhyst_method == mhyst_force_lo) then
                   jhyst_leg(ibin) = jhyst_lo
                   water_a_hyst(ibin) = 0.0
                else if (mhyst_method == mhyst_force_up) then
                   jhyst_leg(ibin)    = jhyst_up
                   water_a_up(ibin)   = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)	
                   water_a_hyst(ibin) = water_a_up(ibin)
                end if
             end if
             
          endif
       endif
       if (irepeat_mosaic == 1) then
          mass_dry_a_bgn(ibin) = mass_dry_a(ibin)
          if ( (jaerosolstate(ibin) .eq. no_aerosol) .or.   &
               (min(mass_dry_a(ibin),vol_dry_a(ibin)) .le. 1.0e-35) ) then
             call calc_aerosol_dry_density( ibin,aer,dens_dry_a)
             dens_dry_a_bgn(ibin) = dens_dry_a(ibin)
          else
             dens_dry_a_bgn(ibin) = mass_dry_a(ibin)/vol_dry_a(ibin)
          end if
          dens_dry_a_bgn(ibin) = max( density_min_allow, &
               min( density_max_allow, dens_dry_a_bgn(ibin) ) )
       end if
       
       if (jaerosolstate(ibin) .eq. no_aerosol) then
          if (msize_framework == msectional) then
             isize = isize_of_ibin(ibin)
             itype = itype_of_ibin(ibin)
             Dp_dry_a(ibin) = dcen_sect(isize,itype)
             Dp_wet_a(ibin) = Dp_dry_a(ibin)
          end if
       end if
       
    enddo
    
    
    

    
    
    
    
    if(mGAS_AER_XFER .eq. mON)then
       
       
       if(mDYNAMIC_SOLVER .eq. mASTEM)then
          call ASTEM( mcall_print_aer,          dtchem,           &
               sigmag_a,  aH2O,     T_K,         RH_pc,        P_atm,                        &
               kappa_nonelectro,                                                             &
               jaerosolstate, flux_s,            flux_l,       volatile_s, iprint_input,     &
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

          if (mosaic_vars_aa%f_mos_fail > 0) then
             return
          endif

          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
       elseif(mDYNAMIC_SOLVER .eq. mLSODE)then
          
          call MOSAIC_LSODE(dtchem)
          
       endif
       
    endif
    
    
    
    
    
    do ibin = 1, nbin_a
       if(jaerosolstate(ibin) .ne. no_aerosol)then
          call conform_aerosol_size( ibin,jaerosolstate,aer,num_a,       Dp_dry_a,  &
               vol_dry_a,mw_aer_mac,dens_aer_mac, mosaic_vars_aa )    
          if (mosaic_vars_aa%f_mos_fail > 0) then
             return
          endif
       endif
    enddo
    
    
    do ibin = 1, nbin_a
       if(jaerosolstate(ibin).ne.no_aerosol) then 
          

          if (mhyst_method == mhyst_uporlo_jhyst .or. & 
              mhyst_method == mhyst_uporlo_waterhyst) then
             
             if(jhyst_leg(ibin) == jhyst_lo)then
                water_a_hyst(ibin) = 0.0
             else
                water_a_up(ibin)   = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)   
                water_a_hyst(ibin) = water_a_up(ibin)
             endif










          else if (mhyst_method == mhyst_force_lo) then
             jhyst_leg(ibin) = jhyst_lo
             water_a_hyst(ibin) = 0.0
          else if (mhyst_method == mhyst_force_up) then
             jhyst_leg(ibin) = jhyst_up
             water_a_up(ibin)   = aerosol_water_up(ibin,electrolyte,aer,kappa_nonelectro,a_zsr)   
             water_a_hyst(ibin) = water_a_up(ibin)
             
          else
             write(errmsg,*) '*** MOSAIC_dynamic_solver - bad mhyst_method =', mhyst_method
             call wrf_error_fatal3("<stdin>",852,&
trim(adjustl(errmsg)))
          endif
          
          
          call calc_dry_n_wet_aerosol_props(                                &
               ibin, jaerosolstate, aer, electrolyte, water_a, num_a,         &  
               dens_comp_a, mw_comp_a, dens_aer_mac, mw_aer_mac, ref_index_a, &  
               Dp_dry_a, Dp_wet_a, dp_core_a,                                 &  
               area_dry_a, area_wet_a, mass_dry_a, mass_wet_a,                &  
               vol_dry_a, vol_wet_a, dens_dry_a, dens_wet_a,                  &  
               ri_shell_a, ri_core_a, ri_avg_a                                )  
          
       endif
       if ( (jaerosolstate(ibin) .eq. no_aerosol) .or.   &
            (min(mass_dry_a(ibin),vol_dry_a(ibin)) .le. 1.0e-35) ) then
          call calc_aerosol_dry_density( ibin,aer,dens_dry_a)
       end if
       dens_dry_a(ibin) = max( density_min_allow, &
            min( density_max_allow, dens_dry_a(ibin) ) )
       
    enddo
    
    if (mcall_print_aer == 1 .or. mcall_print_aer == 2) then
       
       
       
    end if


    call dumpxx( 'zz', dtchem, t_k, p_atm, ah2o, &
         jaerosolstate, jphase, jhyst_leg, &
         aer, gas, num_a, water_a, water_a_hyst, dp_dry_a, &
         mosaic_vars_aa )


    return
  end subroutine MOSAIC_dynamic_solver



  
  
  
  
  
  
  subroutine wall_loss(dtchem,aer,num_a)
    use module_data_mosaic_aero, only: nbin_a_max,naer,jtotal,jsolid,jliquid,   & 
         nbin_a 

    implicit none
    
    real(r8), intent(in) :: dtchem
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    
    integer  :: iaer, ibin
    real(r8) :: kwall


    kwall =  5.55e-5  

    do ibin = 1, nbin_a

       do iaer = 1, naer
          aer(iaer,jtotal,ibin)  = aer(iaer,jtotal,ibin)*exp(-kwall*dtchem)
          aer(iaer,jsolid,ibin)  = aer(iaer,jsolid,ibin)*exp(-kwall*dtchem)
          aer(iaer,jliquid,ibin) = aer(iaer,jliquid,ibin)*exp(-kwall*dtchem)
       enddo

       num_a(ibin) = num_a(ibin)*exp(-kwall*dtchem)

    enddo


    return
  end subroutine wall_loss



  
  
  
  
  
  
  subroutine initialize_mosaic_variables(                                          & 
       jaerosolstate, flux_s, flux_l, volatile_s, phi_volatile_s, phi_volatile_l,  & 
       jphase, kg, electrolyte, activity, mc, mass_dry_a, mass_soluble_a,          &
       dens_dry_a, ma, gam, gam_ratio                                              )

    use module_data_mosaic_aero, only: nbin_a_max,nbin_a,naer,                     &
         ngas_aerchtot, ngas_volatile,                                             &
         nelectrolyte,Ncation,ngas_ioa,jtotal,jsolid,jliquid,nanion


    implicit none

    
    integer, intent(out), dimension(nbin_a_max) :: jaerosolstate,jphase

    real(r8), intent(out), dimension(nbin_a_max) :: mass_dry_a,gam_ratio
    real(r8), intent(out), dimension(nbin_a_max) :: mass_soluble_a,dens_dry_a

    real(r8), intent(out), dimension(ngas_aerchtot,nbin_a_max) :: kg
    real(r8), intent(out), dimension(ngas_volatile,nbin_a_max) :: flux_s
    real(r8), intent(out), dimension(ngas_volatile,nbin_a_max) :: flux_l
    real(r8), intent(out), dimension(ngas_volatile,nbin_a_max) :: volatile_s
    real(r8), intent(out), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_s
    real(r8), intent(out), dimension(ngas_volatile,nbin_a_max) :: phi_volatile_l
    real(r8), intent(out), dimension(nelectrolyte,nbin_a_max)  :: activity,gam
    real(r8), intent(out), dimension(Ncation,nbin_a_max)       :: mc
    real(r8), intent(out), dimension(Nanion,nbin_a_max)        :: ma

    real(r8), intent(out), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

    
    integer iaer, ibin, iv, ja, jc, je

    phi_volatile_l(:,:) = 0.0_r8 

    
    do ibin = 1, nbin_a

       mass_dry_a(ibin)     = 0.0
       mass_soluble_a(ibin) = 0.0
       dens_dry_a(ibin)     =-1.0

       do je = 1, nelectrolyte
          electrolyte(je,jtotal,ibin)  = 0.0
          electrolyte(je,jsolid,ibin)  = 0.0
          electrolyte(je,jliquid,ibin) = 0.0
          activity(je,ibin)            = 0.0
          gam(je,ibin)                 = 0.0
       enddo

       gam_ratio(ibin)   = 0.0

       do iv = 1, ngas_ioa
          flux_s(iv,ibin)   = 0.0
          flux_l(iv,ibin)   = 0.0
          kg(iv,ibin)       = 0.0
          phi_volatile_s(iv,ibin) = 0.0
          phi_volatile_l(iv,ibin) = 0.0
          volatile_s(iv,ibin) = 0.0
       enddo


       jaerosolstate(ibin) = -1     
       jphase(ibin) = 0

       do jc = 1, ncation
          mc(jc,ibin) = 0.0
       enddo

       do ja = 1, nanion
          ma(ja,ibin) = 0.0
       enddo

    enddo   



    return
  end subroutine initialize_mosaic_variables



  subroutine overall_massbal_in( aer, gas, gas_netprod_otrproc, dtchem,            & 
       total_species, tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in,    & 
       tot_ca_in )


    use module_data_mosaic_aero, only: ngas_aerchtot, ngas_volatile,               &
         naer, nbin_a_max, jtotal, nbin_a,                                         &
         ih2so4_g,ihno3_g,ihcl_g,inh3_g,iso4_a,ino3_a,icl_a,inh4_a,ina_a,ica_a      
    
    implicit none

    
    real(r8), intent(in), dimension(ngas_aerchtot) :: gas, gas_netprod_otrproc
    real(r8), intent(in), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(in) :: dtchem

    real(r8), intent(out) :: tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in, tot_ca_in
    real(r8), intent(out), dimension(ngas_volatile) ::total_species


    
    integer ibin

    tot_so4_in = gas(ih2so4_g)
    tot_no3_in = gas(ihno3_g)
    tot_cl_in  = gas(ihcl_g)
    tot_nh4_in = gas(inh3_g)
    tot_na_in  = 0.0
    tot_ca_in  = 0.0

    tot_so4_in = gas(ih2so4_g) + max( gas_netprod_otrproc(ih2so4_g)*dtchem, 0.0_r8 )


    do ibin = 1, nbin_a
       tot_so4_in = tot_so4_in + aer(iso4_a,jtotal,ibin)
       tot_no3_in = tot_no3_in + aer(ino3_a,jtotal,ibin)
       tot_cl_in  = tot_cl_in  + aer(icl_a, jtotal,ibin)
       tot_nh4_in = tot_nh4_in + aer(inh4_a,jtotal,ibin)
       tot_na_in  = tot_na_in  + aer(ina_a,jtotal,ibin)
       tot_ca_in  = tot_ca_in  + aer(ica_a,jtotal,ibin)
    enddo


    total_species(inh3_g) = tot_nh4_in
    total_species(ihno3_g)= tot_no3_in
    total_species(ihcl_g) = tot_cl_in


    return
  end subroutine overall_massbal_in



  subroutine overall_massbal_out( iprint_input, mbin, isteps_ASTEM, aer, gas, &
    tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in, tot_ca_in )




    use module_data_mosaic_aero, only: ngas_aerchtot,naer,nbin_a_max,jtotal,    &
         mYES,mNO,                                                                 &
         nbin_a,                                                                   &
         ih2so4_g,ihno3_g,ihcl_g,inh3_g,iso4_a,ino3_a,icl_a,inh4_a,ina_a,ica_a      

    implicit none

    

    integer, intent(in)    ::  mbin, isteps_ASTEM
    integer, intent(inout) ::  iprint_input

    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout) :: tot_so4_in, tot_no3_in, tot_cl_in, tot_nh4_in, tot_na_in, tot_ca_in

    
    integer ibin
    real(r8) :: tot_so4_out, tot_no3_out, tot_cl_out, tot_nh4_out, tot_na_out, tot_ca_out
    real(r8) :: diff_so4, diff_no3, diff_cl, diff_nh4, diff_na, diff_ca
    real(r8) :: reldiff_so4, reldiff_no3, reldiff_cl, reldiff_nh4, reldiff_na, reldiff_ca



    tot_so4_out = gas(ih2so4_g)
    tot_no3_out = gas(ihno3_g)
    tot_cl_out  = gas(ihcl_g)
    tot_nh4_out = gas(inh3_g)
    tot_na_out  = 0.0
    tot_ca_out  = 0.0

    do ibin = 1, nbin_a
       tot_so4_out = tot_so4_out + aer(iso4_a,jtotal,ibin)
       tot_no3_out = tot_no3_out + aer(ino3_a,jtotal,ibin)
       tot_cl_out  = tot_cl_out  + aer(icl_a,jtotal,ibin)
       tot_nh4_out = tot_nh4_out + aer(inh4_a,jtotal,ibin)
       tot_na_out  = tot_na_out  + aer(ina_a,jtotal,ibin)
       tot_ca_out  = tot_ca_out  + aer(ica_a,jtotal,ibin)
    enddo

    diff_so4 = tot_so4_out - tot_so4_in
    diff_no3 = tot_no3_out - tot_no3_in
    diff_cl  = tot_cl_out  - tot_cl_in
    diff_nh4 = tot_nh4_out - tot_nh4_in
    diff_na  = tot_na_out  - tot_na_in
    diff_ca  = tot_ca_out  - tot_ca_in


    reldiff_so4 = 0.0
    if(tot_so4_in .gt. 1.e-25 .or. tot_so4_out .gt. 1.e-25)then
       reldiff_so4 = diff_so4/max(tot_so4_in, tot_so4_out)
    endif

    reldiff_no3 = 0.0
    if(tot_no3_in .gt. 1.e-25 .or. tot_no3_out .gt. 1.e-25)then
       reldiff_no3 = diff_no3/max(tot_no3_in, tot_no3_out)
    endif

    reldiff_cl = 0.0
    if(tot_cl_in .gt. 1.e-25 .or. tot_cl_out .gt. 1.e-25)then
       reldiff_cl = diff_cl/max(tot_cl_in, tot_cl_out)
    endif

    reldiff_nh4 = 0.0
    if(tot_nh4_in .gt. 1.e-25 .or. tot_nh4_out .gt. 1.e-25)then
       reldiff_nh4 = diff_nh4/max(tot_nh4_in, tot_nh4_out)
    endif

    reldiff_na = 0.0
    if(tot_na_in .gt. 1.e-25 .or. tot_na_out .gt. 1.e-25)then
       reldiff_na = diff_na/max(tot_na_in, tot_na_out)
    endif

    reldiff_ca = 0.0
    if(tot_ca_in .gt. 1.e-25 .or. tot_ca_out .gt. 1.e-25)then
       reldiff_ca = diff_ca/max(tot_ca_in, tot_ca_out)
    endif

    if( abs(reldiff_so4) .gt. 1.e-4 .or.   &
         abs(reldiff_no3) .gt. 1.e-4 .or.   &
         abs(reldiff_cl)  .gt. 1.e-4 .or.   &
         abs(reldiff_nh4) .gt. 1.e-4 .or.   &
         abs(reldiff_na)  .gt. 1.e-4 .or.   &
         abs(reldiff_ca)  .gt. 1.e-4)then


       if(iprint_input .eq. mYES)then
          write(6,*)'*** mbin = ', mbin, '  isteps = ', isteps_ASTEM
          write(6,*)'reldiff_so4 = ', reldiff_so4
          write(6,*)'reldiff_no3 = ', reldiff_no3
          write(6,*)'reldiff_cl  = ', reldiff_cl
          write(6,*)'reldiff_nh4 = ', reldiff_nh4
          write(6,*)'reldiff_na  = ', reldiff_na
          write(6,*)'reldiff_ca  = ', reldiff_ca
          
          iprint_input = mNO
       endif

       

    endif


    return
  end subroutine overall_massbal_out



  
  
  
  
  
  
  
  subroutine check_aerosol_mass(ibin, jaerosolstate,jphase,aer,num_a, mass_dry_a )

    use module_data_mosaic_aero, only: nbin_a_max,naer,jtotal,no_aerosol,       &
         mass_cutoff, mw_aer_mac,                                               &
         iso4_a,ino3_a,icl_a,imsa_a,ico3_a,ica_a,ina_a,inh4_a                       

    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(in), dimension(naer,3,nbin_a_max) :: aer

    
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate,jphase
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a, mass_dry_a

    
    integer iaer
    real(r8) :: drymass, aer_H

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
    enddo
    mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H

    drymass = mass_dry_a(ibin)                      
    mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15      

    if(drymass .lt. mass_cutoff)then                        
       jaerosolstate(ibin) = no_aerosol
       jphase(ibin) = 0
       if(drymass .eq. 0.)num_a(ibin) = 0.0
    endif

    return
  end subroutine check_aerosol_mass



  
  
  
  
  
  
  subroutine conform_aerosol_number(ibin,jaerosolstate,aer,num_a,vol_dry_a, Dp_dry_a)
    
    use module_data_mosaic_constants,  only: pi
    use module_data_mosaic_aero,  only: nbin_a_max,naer,mSECTIONAL,no_aerosol,  &
         jtotal, mw_aer_mac,dens_aer_mac,                                       &
         msize_framework,                                                       &
         iso4_a,ino3_a,icl_a,imsa_a,ico3_a,ica_a,ina_a,inh4_a                    

    use module_data_mosaic_asecthp, only:isize_of_ibin,itype_of_ibin,volumlo_sect,&
         volumhi_sect                                                               

    implicit none

    
    integer, intent(in) :: ibin
    integer, intent(in), dimension(nbin_a_max) :: jaerosolstate

    real(r8), intent(in), dimension(nbin_a_max) :: Dp_dry_a
    real(r8), intent(in), dimension(naer,3,nbin_a_max) :: aer

    
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a,vol_dry_a    


    
    integer :: iaer, isize, itype
    real(r8) :: num_at_dlo, num_at_dhi, numold
    real(r8) :: aer_H
    logical, parameter :: nonsect_set_number_always = .false.

    
    
    
    
    
    
    
    if (msize_framework /= msectional) then
       if (num_a(ibin) > 0.0) return
    end if

    vol_dry_a(ibin)  = 0.0          

    if(jaerosolstate(ibin) .eq. no_aerosol) return


    
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
       vol_dry_a(ibin) = vol_dry_a(ibin) +   &
            aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  
    enddo
    vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H
    vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15                           


    if (msize_framework /= msectional) then
       
       if (num_a(ibin) <= 0.0) then
          num_a(ibin) = vol_dry_a(ibin)/((pi/6.0_r8)*Dp_dry_a(ibin)**3)       
       end if
    else
       
       if (num_a(ibin) <= 0.0) then
          
          num_a(ibin) = vol_dry_a(ibin)/((pi/6.0_r8)*Dp_dry_a(ibin)**3)       
       else
          
          isize = isize_of_ibin( ibin )
          itype = itype_of_ibin( ibin )
          num_at_dlo = vol_dry_a(ibin)/volumlo_sect(isize,itype)
          num_at_dhi = vol_dry_a(ibin)/volumhi_sect(isize,itype)
          numold = num_a(ibin)
          num_a(ibin) = min( num_a(ibin), num_at_dlo )
          num_a(ibin) = max( num_a(ibin), num_at_dhi )
       end if
    end if


    return
  end subroutine conform_aerosol_number



  
  
  
  
  
  
  subroutine calc_aerosol_dry_density(ibin,aer,dens_dry_a)


    use module_data_mosaic_aero,  only: nbin_a_max,naer,jtotal,                 &
         inh4_a,ina_a,ica_a,ico3_a,imsa_a,icl_a,ino3_a,iso4_a,                  & 
         mw_aer_mac, dens_aer_mac

    

    implicit none

    
    integer, intent(in) :: ibin
    real(r8), intent(in), dimension(naer,3,nbin_a_max) :: aer

    
    real(r8), intent(inout), dimension(nbin_a_max) :: dens_dry_a

    
    integer :: iaer
    real(r8) :: aer_H
    real(r8) :: tmpa, tmp_volu, tmp_mass


    
    aer_H = ( 2.*max( 0.0_r8, aer(iso4_a,jtotal,ibin) ) +   &
                 max( 0.0_r8, aer(ino3_a,jtotal,ibin) ) +   &
                 max( 0.0_r8, aer(icl_a,jtotal,ibin) )  +   &
                 max( 0.0_r8, aer(imsa_a,jtotal,ibin) ) +   &
              2.*max( 0.0_r8, aer(ico3_a,jtotal,ibin) ) )   &
          - ( 2.*max( 0.0_r8, aer(ica_a,jtotal,ibin) )  +   &
                 max( 0.0_r8, aer(ina_a,jtotal,ibin) )  +   &
                 max( 0.0_r8, aer(inh4_a,jtotal,ibin) ) )
    aer_H = max( aer_H, 0.0_r8 )

    tmp_mass = aer_H
    tmp_volu = aer_H   

    do iaer = 1, naer
       tmpa = max( 0.0_r8, aer(iaer,jtotal,ibin) ) * mw_aer_mac(iaer)
       tmp_mass = tmp_mass + tmpa                     
       tmp_volu = tmp_volu + tmpa/dens_aer_mac(iaer)  
    enddo

    
    
    if (min(tmp_mass,tmp_volu) >= 1.0e-20) then
       dens_dry_a(ibin) = tmp_mass/tmp_volu   
    else
       dens_dry_a(ibin) = 1.0
    end if

    return
  end subroutine calc_aerosol_dry_density



  
  
  
  
  
  
  subroutine conform_aerosol_size( ibin, jaerosolstate, aer, num_a, Dp_dry_a,     &
       vol_dry_a, mw_aer_mac, dens_aer_mac, mosaic_vars_aa )   


    use module_data_mosaic_constants,  only : piover6,  third
    use module_data_mosaic_aero,  only : nbin_a_max, naer, no_aerosol, jtotal,       &
         inh4_a, ina_a, ica_a, ico3_a, imsa_a, icl_a, ino3_a, iso4_a,                &
         mosaic_vars_aa_type                                                     

    implicit none

    
    integer, intent(in):: ibin
    integer, intent(in), dimension(nbin_a_max) :: jaerosolstate

    real(r8), intent(in), dimension(nbin_a_max) :: num_a
    real(r8), intent(in), dimension(naer) :: mw_aer_mac,dens_aer_mac
    real(r8), intent(inout), dimension(nbin_a_max) ::        Dp_dry_a,vol_dry_a
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer

    type (mosaic_vars_aa_type), intent(inout) :: mosaic_vars_aa

    
    integer iaer
    real(r8) :: num_at_dlo, num_at_dhi
    real(r8) :: aer_H


    vol_dry_a(ibin)  = 0.0          

    if(jaerosolstate(ibin) .eq. no_aerosol) return

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
       vol_dry_a(ibin) = vol_dry_a(ibin) +   &
            aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)       
    enddo
    vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H
    vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15                                


    
    
    
    
    mosaic_vars_aa%f_mos_fail = -1
    if(vol_dry_a(ibin)<0.0_r8) then
       write(202,*)'EXITING due to negative vol_dry_a(',ibin,')=', &
          vol_dry_a(ibin), mosaic_vars_aa%it_mosaic, mosaic_vars_aa%hostgridinfo(1:3)
       mosaic_vars_aa%f_mos_fail = 1
       return
    endif
    Dp_dry_a(ibin) = (vol_dry_a(ibin)/(piover6*num_a(ibin)))**third

    return
  end subroutine conform_aerosol_size



  
  
  
  
  
  
  
  
  
  
  subroutine MTEM_compute_log_gamZ(aH2O,log_gamZ,b_mtem,aw_min)
    use module_data_mosaic_aero, only: nelectrolyte,                            &
         jhno3,jnh4so4,jnh4no3,jnh4cl,jna2so4,jnano3,jnacl,jcano3,jcacl2,jhcl,  &
         jh2so4,jnh4hso4,jlvcite,jnahso4,jna3hso4,jhhso4                            

    implicit none

    
    real(r8), intent(in) :: aH2O
    real(r8), intent(in), dimension(nelectrolyte) :: aw_min
    real(r8), intent(inout), dimension(nelectrolyte,nelectrolyte) :: log_gamZ
    real(r8), intent(in), dimension(6,nelectrolyte,nelectrolyte) :: b_mtem
    
    integer jA
    
    


    
    jA = jhno3
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)


    jA = jhcl
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)


    jA = jnh4so4
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)


    jA = jnh4no3
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jnh4cl
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jna2so4
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)


    jA = jnano3
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jnacl
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jcano3
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jcacl2
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    
    jA = jh2so4
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jhhso4
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jnh4hso4
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jlvcite
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jnahso4
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)


    jA = jna3hso4
    log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3,aH2O,b_mtem,aw_min)
    log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl,aH2O,b_mtem,aw_min)

    return
  end subroutine MTEM_compute_log_gamZ



  subroutine degas_acids(jp,ibin,XT,aer,gas,electrolyte)
    use module_data_mosaic_aero, only: naer,nelectrolyte,nbin_a_max,            &
         ngas_aerchtot,jliquid,jsolid,jtotal,                                      &
         jhno3,jhcl,ihno3_g,ihcl_g,ino3_a,icl_a

    implicit none

    
    integer, intent(in) :: jp, ibin
    real(r8), intent(in) :: XT
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
    
    real(r8) :: ehno3, ehcl



    if(jp .ne. jliquid)then
       write(6,*)'Error in degas_acids'
       write(6,*)'wrong jp'
    endif

    ehno3 = electrolyte(jhno3,jp,ibin)
    ehcl  = electrolyte(jhcl,jp,ibin)

    
    gas(ihno3_g) = gas(ihno3_g) + ehno3
    gas(ihcl_g)  = gas(ihcl_g)  + ehcl

    
    aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) - ehno3
    aer(icl_a, jp,ibin) = aer(icl_a, jp,ibin) - ehcl

    
    aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
         aer(ino3_a,jsolid, ibin)

    aer(icl_a,jtotal,ibin)  = aer(icl_a,jliquid,ibin) +   &
         aer(icl_a,jsolid, ibin)

    electrolyte(jhno3,jp,ibin) = 0.0
    electrolyte(jhcl,jp,ibin)  = 0.0

    return
  end subroutine degas_acids





  
  
  
  
  
  
  subroutine update_thermodynamic_constants(     aH2O,         T_K_in,               & 
       sat_soa,    aH2O_a,   log_gamZ,  Keq_sl,  sigma_water,  Kp_nh4cl,             & 
       Kp_nh4no3,  Kp_nh3,   Keq_ll,    Keq_gl,  Keq_sg,       MDRH_T,               &
       molality0                                                                     )

    use module_data_mosaic_aero, only: r8, nbin_a_max, ngas_volatile,  nelectrolyte,         &
         nrxn_aer_sg, nrxn_aer_gl, nrxn_aer_sl, nrxn_aer_ll, MDRH_T_NUM, d_mdrh_DIM2,        &
         nbin_a, b_mtem, b_zsr, a_zsr, aw_min, d_mdrh,                                       &
         jnh4so4, jlvcite, jnh4hso4, jnh4msa, jnh4no3, jnh4cl, jna2so4, jnahso4, jna3hso4,   &
         jnamsa, jnano3, jnacl, jcacl2, jcano3, jcamsa2, iaro1_g, iaro2_g, ialk1_g, iole1_g, &
         iapi1_g, iapi2_g, ilim1_g, ilim2_g, isoa_first, msoa_flag1,                         &
         use_cam5mam_soa_params

    use module_mosaic_ext, only: bin_molality

    use module_mosaic_soa_vbs, only: soa_vbs_update_thermcons

    implicit none

    
    real(r8), intent(in) :: aH2O, T_K_in

    real(r8), intent(out) :: sigma_water,Kp_nh4cl,Kp_nh4no3,Kp_nh3
    real(r8), intent(out), dimension(nbin_a_max)    :: aH2O_a
    real(r8), intent(out), dimension(ngas_volatile) :: sat_soa
    real(r8), intent(out), dimension(nrxn_aer_ll)   :: Keq_ll
    real(r8), intent(out), dimension(nrxn_aer_sl)   :: Keq_sl
    real(r8), intent(out), dimension(nrxn_aer_gl)   :: Keq_gl
    real(r8), intent(out), dimension(nrxn_aer_sg)   :: Keq_sg
    real(r8), intent(out), dimension(MDRH_T_NUM)    :: MDRH_T
    real(r8), intent(out), dimension(nelectrolyte,nbin_a_max)  :: molality0 
    real(r8), intent(out), dimension(nelectrolyte,nelectrolyte) :: log_gamZ

    
    integer :: ibin, iv, j_index, je
    logical :: use_sorgam_soa_species


    real(r8) :: Po_soa(ngas_volatile)
    real(r8) :: rt
    real(r8) :: sat_factor
    real(r8) :: T_K		
    real(r8) :: tr, term
    
    



    T_K = max( 220.0_r8, T_K_in )
    T_K = min( 330.0_r8, T_K )

    tr = 298.15                   
    rt = 82.056*T_K/(1.e9*1.e6)   

    
    Keq_gl(1)= 1.0                                        
    Keq_gl(2)= fn_Keq(57.64d0, 13.79d0, -5.39d0,T_K)*rt     
    Keq_gl(3)= fn_Keq(2.63d6,  29.17d0, 16.83d0,T_K)*rt     
    Keq_gl(4)= fn_Keq(2.00d6,  30.20d0, 19.91d0,T_K)*rt     

    
    Keq_ll(1)= fn_Keq(1.0502d-2, 8.85d0, 25.14d0,T_K)      
    Keq_ll(2)= fn_Keq(1.805d-5, -1.50d0, 26.92d0,T_K)      
    Keq_ll(3)= fn_Keq(1.01d-14,-22.52d0, 26.92d0,T_K)      


    Kp_nh3   = Keq_ll(3)/(Keq_ll(2)*Keq_gl(2))
    Kp_nh4no3= Kp_nh3/Keq_gl(3)
    Kp_nh4cl = Kp_nh3/Keq_gl(4)


    
    Keq_sg(1)= fn_Keq(4.72d-17,-74.38d0,6.12d0,T_K)/rt**2  
    Keq_sg(2)= fn_Keq(8.43d-17,-71.00d0,2.40d0,T_K)/rt**2  


    
    Keq_sl(jnh4so4) = fn_Keq(1.040d0,-2.65d0, 38.57d0, T_K)  
    Keq_sl(jlvcite) = fn_Keq(11.8d0, -5.19d0, 54.40d0, T_K)  
    Keq_sl(jnh4hso4)= fn_Keq(117.0d0,-2.87d0, 15.83d0, T_K)  
    Keq_sl(jnh4msa) = 1.e15                                      
    Keq_sl(jnh4no3) = fn_Keq(12.21d0,-10.4d0, 17.56d0, T_K)  
    Keq_sl(jnh4cl)  = fn_Keq(17.37d0,-6.03d0, 16.92d0, T_K)  
    Keq_sl(jna2so4) = fn_Keq(0.491d0, 0.98d0, 39.75d0, T_K)  
    Keq_sl(jnahso4) = fn_Keq(313.0d0, 0.8d0,  14.79d0, T_K)  
    Keq_sl(jna3hso4)= 1.e15                                      
    Keq_sl(jnamsa)  = 1.e15                                      
    Keq_sl(jnano3)  = fn_Keq(11.95d0,-8.22d0, 16.01d0, T_K)  
    Keq_sl(jnacl)   = fn_Keq(38.28d0,-1.52d0, 16.89d0, T_K)  
    Keq_sl(jcacl2)  = fn_Keq(8.0d11,  32.84d0,44.79d0, T_K)  
    Keq_sl(jcano3)  = fn_Keq(4.31d5,   7.83d0,42.01d0, T_K)  
    Keq_sl(jcamsa2) = 1.e15                                



    
    if (msoa_flag1 == 1) then
       use_sorgam_soa_species = .true.
    else
       use_sorgam_soa_species = .false.
    end if

    
    Po_soa(:) = 1.0e5_r8  
    if ( use_sorgam_soa_species ) then
    Po_soa(iaro1_g) = fn_Po(5.7d-5, 156.0d0, T_K) 
    Po_soa(iaro2_g) = fn_Po(1.6d-3, 156.0d0, T_K) 
    Po_soa(ialk1_g) = fn_Po(5.0d-6, 156.0d0, T_K) 
    Po_soa(iole1_g) = fn_Po(5.0d-6, 156.0d0, T_K) 
    Po_soa(iapi1_g) = fn_Po(4.0d-6, 156.0d0, T_K) 
    Po_soa(iapi2_g) = fn_Po(1.7d-4, 156.0d0, T_K) 
    Po_soa(ilim1_g) = fn_Po(2.5d-5, 156.0d0, T_K) 
    Po_soa(ilim2_g) = fn_Po(1.2d-4, 156.0d0, T_K) 
    end if

    sat_soa(:) = 0.0_r8
    sat_factor = 0.5  
    do iv = isoa_first, ngas_volatile
       
       sat_soa(iv) = sat_factor * 1.e9*Po_soa(iv)/(8.314*T_K)  
    enddo

    if ( ( use_cam5mam_soa_params > 0 ) .and. &
         ( 1 <= ilim2_g .and. ilim2_g <= ngas_volatile ) ) then 
       Po_soa(ilim2_g) = fn_Po(1.0d-10, 156.0d0, T_K) 
       sat_soa(ilim2_g) = 1.e9*Po_soa(ilim2_g)/(8.314*T_K)  
    end if

    
    
    

    if (msoa_flag1 >= 1000) call soa_vbs_update_thermcons( t_k, po_soa, sat_soa )



    
    term = (647.15 - T_K)/647.15
    sigma_water = 0.2358*term**1.256 * (1. - 0.625*term) 

    
    do j_index = 1, 63
       MDRH_T(j_index) = drh_mutual(j_index,T_K)
    enddo



    
    do ibin = 1, nbin_a
       aH2O_a(ibin) = aH2O                        

      do je = 1, nelectrolyte
        molality0(je,ibin) = bin_molality(je,ibin,aH2O_a,b_zsr,a_zsr,aw_min)  
      enddo

    enddo

    call MTEM_compute_log_gamZ(aH2O,log_gamZ,b_mtem,aw_min)              


    return
  end subroutine update_thermodynamic_constants





  
  
  
  
  
  



  
  function fn_Keq(Keq_298, a, b, T)
    implicit none
    real(r8) :: fn_Keq
    
    real(r8) :: Keq_298, a, b, T
    
    real(r8) :: tt


    tt = 298.15/T
    fn_Keq = Keq_298*exp(a*(tt-1.)+b*(1.+log(tt)-tt))

    return
  end function fn_Keq
  



  
  function fn_Po(Po_298, DH, T)   
    implicit none
    real(r8) :: fn_Po
    
    real(r8) :: Po_298, DH, T
    

    fn_Po = Po_298*exp(-(DH/8.314e-3)*(1./T - 3.354016435e-3))

    return
  end function fn_Po
  



  
  function drh_mutual(j_index,T_K)            
    use module_data_mosaic_aero, only: d_mdrh

    implicit none


    
    integer,  intent(in) ::  j_index
    real(r8), intent(in) :: T_K

    
    integer j
    real(r8) :: drh_mutual

    j = j_index

    if(j_index .eq. 7 .or. j_index .eq. 8 .or.   &
         (j_index.ge. 34 .and. j_index .le. 51))then

       drh_mutual = 10.0  

    else

       drh_mutual =  d_mdrh(j,1) + T_K*   &
            (d_mdrh(j,2) + T_K*   &
            (d_mdrh(j,3) + T_K*   &
            d_mdrh(j,4) )) + 1.0


       drh_mutual = max(   0.0_r8, drh_mutual )
       drh_mutual = min( 100.0_r8, drh_mutual )

    endif


    return
  end function drh_mutual
  







  
  function fnlog_gamZ(jA,jE,aH2O,b_mtem,aw_min) 
    use module_data_mosaic_aero, only: nelectrolyte

    implicit none

    real(r8) :: fnlog_gamZ
    
    integer, intent(in) :: jA, jE
    real(r8), intent(in) :: aH2O
    real(r8), intent(in),dimension(nelectrolyte) :: aw_min
    real(r8), intent(in), dimension(6,nelectrolyte,nelectrolyte) :: b_mtem
    
    real(r8) :: aw


    aw = max(aH2O, aw_min(jE))

    fnlog_gamZ = b_mtem(1,jA,jE) + aw*   &
         (b_mtem(2,jA,jE) + aw*   &
         (b_mtem(3,jA,jE) + aw*   &
         (b_mtem(4,jA,jE) + aw*   &
         (b_mtem(5,jA,jE) + aw*   &
         b_mtem(6,jA,jE) ))))

    return
  end function fnlog_gamZ
  



  
  
  
  
  
  subroutine quadratix(a,b,c, qx1,qx2)
    implicit none
    
    real(r8) :: a, b, c, qx1, qx2
    
    real(r8) :: x, dum


    if(b .ne. 0.0)then
       x = 4.*(a/b)*(c/b)
    else
       x = 1.e+6
    endif

    if(abs(x) .lt. 1.e-6)then
       dum = ( (0.5*x) +   &
            (0.125*x**2) +   &
            (0.0625*x**3) )

       qx1 = (-0.5*b/a)*dum
       qx2 = -b/a - qx1

    else

       qx1 = ((-b)+sqrt((b*b)-(4.*a*c)))/   &
            (2.*a)
       qx2 = ((-b)-sqrt((b*b)-(4.*a*c)))/   &
            (2.*a)

    endif

    return
  end subroutine quadratix



  
  
  
  
  
  
  subroutine aerosol_optical_properties(                                 &
          gas, aer, num_a, water_a,                                      & 
          dens_comp_a, mw_comp_a, dens_aer_mac, mw_aer_mac, ref_index_a, & 
          Dp_dry_a, Dp_wet_a, dp_core_a,                                 & 
          ri_shell_a, ri_core_a, ri_avg_a, jaerosolstate, jphase,        &
          tot_cl_in, tot_nh4_in, tot_no3_in, XT, area_dry_a, area_wet_a, &
          dens_dry_a, dens_wet_a, mass_dry_a, mass_wet_a, vol_dry_a,     &
          vol_wet_a, total_species, electrolyte     ) 

    use module_data_mosaic_aero, only: &
       icl_a, inh4_a, ino3_a, ihcl_g, inh3_g, ihno3_g, jtotal, &
       naer, naercomp, nbin_a, nbin_a_max, nelectrolyte,       &
       ngas_aerchtot, ngas_volatile, no_aerosol
    use module_mosaic_ext, only: calc_dry_n_wet_aerosol_props, &
         conform_electrolytes

    implicit none

    
    integer,  intent(inout), dimension(nbin_a_max) :: jaerosolstate, jphase 

    real(r8), intent(in), dimension(naer)       :: dens_aer_mac, mw_aer_mac
    real(r8), intent(in), dimension(naercomp)   :: dens_comp_a,mw_comp_a

    real(r8), intent(inout) :: tot_cl_in, tot_nh4_in, tot_no3_in, XT
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
    real(r8), intent(inout), dimension(nbin_a_max) :: water_a
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a
    real(r8), intent(inout), dimension(nbin_a_max) :: Dp_dry_a
    real(r8), intent(inout), dimension(nbin_a_max) :: Dp_wet_a, dp_core_a
    real(r8), intent(inout), dimension(nbin_a_max) :: area_dry_a, area_wet_a
    real(r8), intent(inout), dimension(nbin_a_max) :: dens_dry_a, dens_wet_a
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a, mass_wet_a
    real(r8), intent(inout), dimension(nbin_a_max) :: vol_dry_a, vol_wet_a
    real(r8), intent(inout), dimension(ngas_volatile) :: total_species    
    real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte


    complex,  intent(in), dimension(naercomp)   :: ref_index_a
    complex,  intent(inout), dimension(nbin_a_max) :: ri_shell_a, ri_avg_a, ri_core_a
    
    
    integer iaer, ibin, je, k

    
    do ibin = 1, nbin_a
       do je = 1, nelectrolyte
          electrolyte(je,jtotal,ibin)  = 0.0
       enddo
       jaerosolstate(ibin) = -1   
    enddo

    
    total_species(:) = 0.0_r8
    tot_no3_in = gas(ihno3_g)
    tot_cl_in  = gas(ihcl_g)
    tot_nh4_in = gas(inh3_g)
    do ibin = 1, nbin_a
       tot_no3_in = tot_no3_in + aer(ino3_a,jtotal,ibin)
       tot_cl_in  = tot_cl_in  + aer(icl_a, jtotal,ibin)
       tot_nh4_in = tot_nh4_in + aer(inh4_a,jtotal,ibin)
    enddo
    total_species(inh3_g) = tot_nh4_in
    total_species(ihno3_g)= tot_no3_in
    total_species(ihcl_g) = tot_cl_in


    
    do  ibin = 1, nbin_a
       
       call check_aerosol_mass( ibin, jaerosolstate, jphase, aer, num_a, mass_dry_a )
       
       if(jaerosolstate(ibin) .ne. no_aerosol) then
          
          
          call conform_electrolytes( jtotal, ibin, XT, aer, gas, electrolyte, total_species, tot_cl_in )
          
          
          call check_aerosol_mass( ibin, jaerosolstate, jphase, aer, num_a, mass_dry_a )
          
          if(jaerosolstate(ibin) .ne. no_aerosol) then
             
             call conform_aerosol_number( ibin, jaerosolstate, aer, num_a, vol_dry_a, Dp_dry_a)
             
             
             call calc_dry_n_wet_aerosol_props(                                &
                  ibin, jaerosolstate, aer, electrolyte, water_a, num_a,         &  
                  dens_comp_a, mw_comp_a, dens_aer_mac, mw_aer_mac, ref_index_a, &  
                  Dp_dry_a, Dp_wet_a, dp_core_a,                                 &  
                  area_dry_a, area_wet_a, mass_dry_a, mass_wet_a,                &  
                  vol_dry_a, vol_wet_a, dens_dry_a, dens_wet_a,                  &  
                  ri_shell_a, ri_core_a, ri_avg_a                                )  
          endif
       endif
       
    enddo
    
    return
  end subroutine aerosol_optical_properties


end module module_mosaic_box_aerchem
 
