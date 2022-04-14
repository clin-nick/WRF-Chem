  module module_mosaic_aerchem_intr


  implicit none


  contains


  
  subroutine aerchemistry(                                    &
     idiagbb_host,                                            & 
     hostgridinfo, it_host, it_mosaic, dtchem_in,             & 
     pr_atm, rh, te, cair_mol_m3, cair_mol_cc, swdownbox,     &
     jaerosolstate, jaerosolstate_bgn, jhyst_leg,             & 
     rbox, dp_dry_a, dp_wet_a, sigmag_a,                      & 
     gas_avg, gas_netprod_otrproc,                            & 
     mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a,  &
     aH2O_a, gam_ratio, iter_mesa_out                         ) 



  use module_data_mosaic_kind, only: r8
  use module_data_mosaic_main, only: &
       m_partmc_mosaic, ntot_max, ntot_used
  use module_data_mosaic_aero, only : &
       mosaic_vars_aa_type, &
       dens_aer_mac, &
       mw_aer_mac, mw_comp_a, msectional, msize_framework, &
       naer, nbin_a, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       nmax_astem, nmax_mesa, nsalt, &
       use_cam5mam_soa_params, use_cam5mam_accom_coefs
  use module_mosaic_box_aerchem, only: mosaic_box_aerchemistry


  
  integer,  intent(in)  :: idiagbb_host
  integer,  intent(in)  :: hostgridinfo(6), it_host, it_mosaic
  real(r8), intent(in)  :: dtchem_in, pr_atm, rh, te, cair_mol_m3, cair_mol_cc, swdownbox

  integer, intent(inout),  dimension(nbin_a_max) :: jaerosolstate, jaerosolstate_bgn
  integer, intent(inout),  dimension(nbin_a_max) :: jhyst_leg

  real(r8), intent(inout), dimension(ntot_used)     :: rbox
  real(r8), intent(inout), dimension(nbin_a_max)    :: dp_dry_a, dp_wet_a, sigmag_a
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas_avg  
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas_netprod_otrproc
            
            
            
            
  real(r8), intent(inout), dimension(nbin_a_max)    :: mass_dry_a_bgn, mass_dry_a
  real(r8), intent(inout), dimension(nbin_a_max)    :: dens_dry_a_bgn, dens_dry_a

  integer,  intent(out),   dimension(nbin_a_max)    :: iter_mesa_out
  real(r8), intent(out),   dimension(nbin_a_max)    :: aH2O_a, gam_ratio

  
  character(len=250) :: infile, tmp_str

  logical :: debug_mosaic = .false.

  integer :: ierr, ibin, igas, iaer, istate, iaer_in, istate_in, ibin_in
  integer :: ierror_grp1, ierror_grp2, istop_mosaic_error_grp1, istop_mosaic_error_grp2
  integer :: mcall_load_mosaic_parameters, mcall_print_aer_in
  integer :: n
  integer :: unitn

  real(r8) :: dtchem, RH_pc, aH2O, P_atm, T_K, aer_tmp
  real(r8), dimension(naer,3,nbin_a_max) :: aer
  real(r8), dimension(ngas_aerchtot)     :: gas
  real(r8), dimension(nbin_a_max)        :: num_a, water_a, water_a_hyst
  real(r8), dimension(naer)              :: kappa_nonelectro
  real(r8)                               :: uptkrate_h2so4  

  real(r8)                               :: xsv_misc(5)
  real(r8), dimension(naer,3,nbin_a_max) :: xsv_aer
  real(r8), dimension(ngas_aerchtot)     :: xsv_gas, xsv_gasavg, xsv_gasprod
  real(r8), dimension(nbin_a_max)        :: xsv_num, xsv_water, xsv_dpdry, xsv_dpwet, xsv_sigmag
  integer,  dimension(nbin_a_max)        :: jsv_jhyst, jsv_jstate

  type (mosaic_vars_aa_type) :: mosaic_vars_aa


  dtchem      = dtchem_in
  RH_pc       = RH                                    
  aH2O        = 0.01_r8*RH_pc                         
  P_atm       = pr_atm                                
  T_K         = te                                    

  
  
  
  
  
  
  if (it_mosaic <= 1) then
     mcall_load_mosaic_parameters = 1
     mcall_print_aer_in = 2
  else
     mcall_load_mosaic_parameters = 0
     mcall_print_aer_in = 1
  end if

  call set_kappa_nonelectro( kappa_nonelectro )

  
  call map_mosaic_species_aerchem_box( 0, jaerosolstate,  &          
       rbox, aer, gas, jhyst_leg, num_a, Dp_dry_a,        &
       sigmag_a, water_a, water_a_hyst, cair_mol_m3       )

  
  
  xsv_misc(1) = ah2o
  xsv_misc(2) = t_k
  xsv_misc(3) = p_atm
  xsv_misc(4) = rh_pc
  xsv_misc(5) = dtchem
  xsv_num(1:nbin_a_max)    = num_a(1:nbin_a_max)
  xsv_water(1:nbin_a_max)  = water_a(1:nbin_a_max)
  xsv_dpdry(1:nbin_a_max)  = dp_dry_a(1:nbin_a_max)
  xsv_dpwet(1:nbin_a_max)  = dp_wet_a(1:nbin_a_max)
  xsv_sigmag(1:nbin_a_max) = sigmag_a(1:nbin_a_max)
  jsv_jhyst(1:nbin_a_max)  = jhyst_leg(1:nbin_a_max)
  jsv_jstate(1:nbin_a_max) = jaerosolstate(1:nbin_a_max)
  xsv_gas(1:ngas_aerchtot)     = gas(1:ngas_aerchtot)
  xsv_gasavg(1:ngas_aerchtot)  = gas_avg(1:ngas_aerchtot)
  xsv_gasprod(1:ngas_aerchtot) = gas_netprod_otrproc(1:ngas_aerchtot)
  xsv_aer(1:naer,1:3,1:nbin_a_max) = aer(1:naer,1:3,1:nbin_a_max)


  
  
  
  
  if(debug_mosaic) then
     call wrf_error_fatal3("<stdin>",137,&
'module_mosaic_aerchem_intr - debug_mosaic must be false' )

     
     use_cam5mam_soa_params  = 1
     use_cam5mam_accom_coefs = 1

     
     
     
     unitn = 101
     infile = 'mosaic_error_48.bin'
     open( unitn, file=trim(infile), status='old', form='unformatted', CONVERT = 'BIG_ENDIAN' )
     
     read(unitn)aH2O
     read(unitn)T_K
     read(unitn)P_atm
     read(unitn)RH_pc
     read(unitn)dtchem
     
     do ibin = 1, nbin_a_max
        read(unitn)num_a(ibin),water_a(ibin),Dp_dry_a(ibin),        &
             sigmag_a(ibin),dp_wet_a(ibin),jhyst_leg(ibin),          &
             jaerosolstate(ibin)
     end do
     
     do igas = 1, ngas_aerchtot
        if (igas <= ngas_volatile) then
           read(unitn) gas(igas), gas_avg(igas), gas_netprod_otrproc(igas)
        else
           gas(igas) = 0.0 ; gas_avg(igas) = 0.0 ; gas_netprod_otrproc(igas) = 0.0
        end if
     enddo
     
     do ibin = 1, nbin_a_max
        do istate = 1, 3
           do iaer = 1 , naer
              read(unitn)iaer_in,istate_in,ibin_in, aer_tmp
              aer(iaer_in,istate_in,ibin_in) = aer_tmp                    
           end do
        end do
     end do
     close(unitn)

  endif
  
  


  
  

 
 


  allocate( mosaic_vars_aa%iter_mesa(nbin_a_max), stat=ierr )
  if (ierr /= 0) then
     call wrf_error_fatal3("<stdin>",195,&
'*** subr aerchemistry - allocate error for mosaic_vars_aa%iter_mesa')
  end if
  mosaic_vars_aa%it_host = it_host
  mosaic_vars_aa%it_mosaic = it_mosaic
  mosaic_vars_aa%hostgridinfo(1:6) = hostgridinfo(1:6)
  mosaic_vars_aa%idiagbb_host = idiagbb_host
  mosaic_vars_aa%f_mos_fail = -1
  mosaic_vars_aa%isteps_astem = 0
  mosaic_vars_aa%isteps_astem_max = 0
  mosaic_vars_aa%jastem_call = 0
  mosaic_vars_aa%jastem_fail = -1
  mosaic_vars_aa%jmesa_call = 0
  mosaic_vars_aa%jmesa_fail = 0
  mosaic_vars_aa%niter_mesa_max = 0
  mosaic_vars_aa%nmax_astem = nmax_astem
  mosaic_vars_aa%nmax_mesa = nmax_mesa
  mosaic_vars_aa%fix_astem_negative = 0
  mosaic_vars_aa%fix_astem_negative = 1
  
  mosaic_vars_aa%flag_itr_kel = .false.
  
  mosaic_vars_aa%zero_water_flag = .false.
  mosaic_vars_aa%cumul_steps_astem = 0.0_r8
  mosaic_vars_aa%niter_mesa = 0.0_r8
  mosaic_vars_aa%xnerr_astem_negative(:,:) = 0.0_r8
  mosaic_vars_aa%iter_mesa(1:nbin_a_max) = 0
  mosaic_vars_aa%swdown = swdownbox


  call mosaic_box_aerchemistry(              aH2O,               T_K,            &
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


     istop_mosaic_error_grp1 = 0
     istop_mosaic_error_grp2 = 0
     ierror_grp1 = 0
     ierror_grp2 = 0


  if (mosaic_vars_aa%jastem_fail > 0 .or. mosaic_vars_aa%zero_water_flag .or. mosaic_vars_aa%f_mos_fail > 0) then
     
     print '(2a/8i12)', '*** subr aerchemistry - ', &
        'astem_fail or zero_water error', &
        it_host, it_mosaic, hostgridinfo(1:6)
     write(tmp_str,*) 'Error in Mosaic,jASTEM_fail= ', mosaic_vars_aa%jASTEM_fail, &
        ' zero_water_flag: ', mosaic_vars_aa%zero_water_flag, &
        '  f_mos_fail:', mosaic_vars_aa%f_mos_fail
     print*, trim(adjustl(tmp_str))
     ierror_grp1 = 1
  endif

  iter_mesa_out(1:nbin_a_max) = mosaic_vars_aa%iter_mesa(1:nbin_a_max)

  if ( sum(mosaic_vars_aa%xnerr_astem_negative(:,:)) > 0.0_r8 ) then
     print '(2a/8i12)', '*** subr aerchemistry - ', &
        'astem_negative error', &
        it_host, it_mosaic, hostgridinfo(1:6)
     do n = 1, 4
        print '(i2,1p,5e10.2)', n, mosaic_vars_aa%xnerr_astem_negative(:,n)
     end do
     ierror_grp2 = 1
  end if

  if ( it_mosaic <= 3       .and.  hostgridinfo(2) == 1 .and. &
       hostgridinfo(3) == 1 .and.  hostgridinfo(4) == 1 ) then
     print '(2a/8i12)', '*** subr aerchemistry - ', &
        'artificial error', &
        it_host, it_mosaic, hostgridinfo(1:6)
     do n = 1, 4
        print '(i2,1p,5e10.2)', n, mosaic_vars_aa%xnerr_astem_negative(:,n)
     end do
     ierror_grp2 = 1
  end if

  if (ierror_grp1 > 0 .or. ierror_grp2 > 0) then
     
     print '(a)', 'naer, nbin, ngas'
     print '(8i14      )', naer, nbin_a, ngas_aerchtot
     print '(a)', 'jhyst'
     print '(8i14      )', jsv_jhyst(1:nbin_a)
     print '(8i14      )', jhyst_leg(1:nbin_a)
     print '(a)', 'jstate'
     print '(8i14      )', jsv_jstate(1:nbin_a)
     print '(8i14      )', jaerosolstate(1:nbin_a)
     print '(a)', 'misc'
     print '(1p,4e28.20)', xsv_misc(1:5)
     print '(1p,4e28.20)', ah2o, t_k, p_atm, rh_pc, dtchem
     print '(a)', 'num'
     print '(1p,4e28.20)', xsv_num(1:nbin_a)
     print '(1p,4e28.20)', num_a(1:nbin_a)
     print '(a)', 'water'
     print '(1p,4e28.20)', xsv_water(1:nbin_a)
     print '(1p,4e28.20)', water_a(1:nbin_a)
     print '(a)', 'dpdry'
     print '(1p,4e28.20)', xsv_dpdry(1:nbin_a)
     print '(1p,4e28.20)', dp_dry_a(1:nbin_a)
     print '(a)', 'dpwet'
     print '(1p,4e28.20)', xsv_dpwet(1:nbin_a)
     print '(1p,4e28.20)', dp_wet_a(1:nbin_a)
     print '(a)', 'sigmag'
     print '(1p,4e28.20)', xsv_sigmag(1:nbin_a)
     print '(1p,4e28.20)', sigmag_a(1:nbin_a)
     print '(a)', 'gas'
     print '(1p,4e28.20)', xsv_gas(1:ngas_aerchtot)
     print '(1p,4e28.20)', gas(1:ngas_aerchtot)
     print '(a)', 'gasavg'
     print '(1p,4e28.20)', xsv_gasavg(1:ngas_aerchtot)
     print '(1p,4e28.20)', gas_avg(1:ngas_aerchtot)
     print '(a)', 'gasprod'
     print '(1p,4e28.20)', xsv_gasprod(1:ngas_aerchtot)
     print '(1p,4e28.20)', gas_netprod_otrproc(1:ngas_aerchtot)
     print '(a)', 'aer'
     print '(1p,4e28.20)', xsv_aer(1:naer,1:3,1:nbin_a)
     print '(1p,4e28.20)', aer(1:naer,1:3,1:nbin_a)
  end if

  if ( (istop_mosaic_error_grp1 > 0 .and. ierror_grp1 > 0) .or. &
       (istop_mosaic_error_grp2 > 0 .and. ierror_grp1 > 0) ) then
     
     call wrf_error_fatal3("<stdin>",325,&
'Fortran Stop in subr aerchemistry')
     
  end if

  deallocate( mosaic_vars_aa%iter_mesa, stat=ierr )
  if (ierr /= 0) then
     print '(2a/8i12)', '*** subr aerchemistry - ', &
        'deallocate error for mosaic_vars_aa%iter_mesa', &
        it_host, it_mosaic, hostgridinfo(1:6)
     call wrf_error_fatal3("<stdin>",335,&
'Fortran Stop in subr aerchemistry')
     
  end if


  
  call map_mosaic_species_aerchem_box( 1, jaerosolstate, &          
       rbox, aer, gas, jhyst_leg, num_a, Dp_dry_a,       &
       sigmag_a, water_a, water_a_hyst, cair_mol_m3      )


  return
  end subroutine aerchemistry




      
      subroutine set_kappa_nonelectro( kappa_nonelectro )

      use module_data_mosaic_kind, only : r8

      use module_data_mosaic_aero, only : &
           ibc_a, ioc_a, ilim2_a, ioin_a, &
           ipcg1_b_c_a,  ipcg2_b_c_a,  ipcg3_b_c_a,  ipcg4_b_c_a, &
           ipcg5_b_c_a,  ipcg6_b_c_a,  ipcg7_b_c_a,  ipcg8_b_c_a,  ipcg9_b_c_a, &
           ipcg1_b_o_a,  ipcg2_b_o_a,  ipcg3_b_o_a,  ipcg4_b_o_a, &
           ipcg5_b_o_a,  ipcg6_b_o_a,  ipcg7_b_o_a,  ipcg8_b_o_a,  ipcg9_b_o_a, &
           iopcg1_b_c_a, iopcg2_b_c_a, iopcg3_b_c_a, iopcg4_b_c_a, &
           iopcg5_b_c_a, iopcg6_b_c_a, iopcg7_b_c_a, iopcg8_b_c_a, &
           iopcg1_b_o_a, iopcg2_b_o_a, iopcg3_b_o_a, iopcg4_b_o_a, &
           iopcg5_b_o_a, iopcg6_b_o_a, iopcg7_b_o_a, iopcg8_b_o_a, &
           ipcg1_f_c_a,  ipcg2_f_c_a,  ipcg3_f_c_a,  ipcg4_f_c_a, &
           ipcg5_f_c_a,  ipcg6_f_c_a,  ipcg7_f_c_a,  ipcg8_f_c_a,  ipcg9_f_c_a, &
           ipcg1_f_o_a,  ipcg2_f_o_a,  ipcg3_f_o_a,  ipcg4_f_o_a, &
           ipcg5_f_o_a,  ipcg6_f_o_a,  ipcg7_f_o_a,  ipcg8_f_o_a,  ipcg9_f_o_a, &
           iopcg1_f_c_a, iopcg2_f_c_a, iopcg3_f_c_a, iopcg4_f_c_a, &
           iopcg5_f_c_a, iopcg6_f_c_a, iopcg7_f_c_a, iopcg8_f_c_a, &
           iopcg1_f_o_a, iopcg2_f_o_a, iopcg3_f_o_a, iopcg4_f_o_a, &
           iopcg5_f_o_a, iopcg6_f_o_a, iopcg7_f_o_a, iopcg8_f_o_a, &
           iant1_c_a,  iant2_c_a,  iant3_c_a,  iant4_c_a, &
           iant1_o_a,  iant2_o_a,  iant3_o_a,  iant4_o_a, &
           ibiog1_c_a, ibiog2_c_a, ibiog3_c_a, ibiog4_c_a, &
           ibiog1_o_a, ibiog2_o_a, ibiog3_o_a, ibiog4_o_a, &
           ismpa_a, ismpbb_a, &
           msoa_flag1, naer

      use module_data_mosaic_asect, only: &
         hygro_oin_aer, hygro_oc_aer, hygro_bc_aer,  &
         hygro_pcg1_b_c_aer,  hygro_pcg2_b_c_aer,  hygro_pcg3_b_c_aer,  &
         hygro_pcg4_b_c_aer,  hygro_pcg5_b_c_aer,  hygro_pcg6_b_c_aer,  &
         hygro_pcg7_b_c_aer,  hygro_pcg8_b_c_aer,  hygro_pcg9_b_c_aer,  &
         hygro_pcg1_b_o_aer,  hygro_pcg2_b_o_aer,  hygro_pcg3_b_o_aer,  &
         hygro_pcg4_b_o_aer,  hygro_pcg5_b_o_aer,  hygro_pcg6_b_o_aer,  &
         hygro_pcg7_b_o_aer,  hygro_pcg8_b_o_aer,  hygro_pcg9_b_o_aer,  &
         hygro_opcg1_b_c_aer, hygro_opcg2_b_c_aer, hygro_opcg3_b_c_aer,  &
         hygro_opcg4_b_c_aer, hygro_opcg5_b_c_aer, hygro_opcg6_b_c_aer,  &
         hygro_opcg7_b_c_aer, hygro_opcg8_b_c_aer,  &
         hygro_opcg1_b_o_aer, hygro_opcg2_b_o_aer, hygro_opcg3_b_o_aer,  &
         hygro_opcg4_b_o_aer, hygro_opcg5_b_o_aer, hygro_opcg6_b_o_aer,  &
         hygro_opcg7_b_o_aer, hygro_opcg8_b_o_aer,  &
         hygro_pcg1_f_c_aer,  hygro_pcg2_f_c_aer,  hygro_pcg3_f_c_aer,  &
         hygro_pcg4_f_c_aer,  hygro_pcg5_f_c_aer,  hygro_pcg6_f_c_aer,  &
         hygro_pcg7_f_c_aer,  hygro_pcg8_f_c_aer,  hygro_pcg9_f_c_aer,  &
         hygro_pcg1_f_o_aer,  hygro_pcg2_f_o_aer,  hygro_pcg3_f_o_aer,  &
         hygro_pcg4_f_o_aer,  hygro_pcg5_f_o_aer,  hygro_pcg6_f_o_aer,  &
         hygro_pcg7_f_o_aer,  hygro_pcg8_f_o_aer,  hygro_pcg9_f_o_aer,  &
         hygro_opcg1_f_c_aer, hygro_opcg2_f_c_aer, hygro_opcg3_f_c_aer,  &
         hygro_opcg4_f_c_aer, hygro_opcg5_f_c_aer, hygro_opcg6_f_c_aer,  &
         hygro_opcg7_f_c_aer, hygro_opcg8_f_c_aer,  &
         hygro_opcg1_f_o_aer, hygro_opcg2_f_o_aer, hygro_opcg3_f_o_aer,  &
         hygro_opcg4_f_o_aer, hygro_opcg5_f_o_aer, hygro_opcg6_f_o_aer,  &
         hygro_opcg7_f_o_aer, hygro_opcg8_f_o_aer,  &
         hygro_ant1_c_aer,  hygro_ant2_c_aer,  hygro_ant3_c_aer,  hygro_ant4_c_aer,  &
         hygro_ant1_o_aer,  hygro_ant2_o_aer,  hygro_ant3_o_aer,  hygro_ant4_o_aer,  &
         hygro_biog1_c_aer, hygro_biog2_c_aer, hygro_biog3_c_aer, hygro_biog4_c_aer,  &
         hygro_biog1_o_aer, hygro_biog2_o_aer, hygro_biog3_o_aer, hygro_biog4_o_aer,  &
         hygro_smpa_aer, hygro_smpbb_aer

      real(r8), dimension(naer), intent(inout) :: kappa_nonelectro

      real(r8) :: kappa_soa

      kappa_nonelectro(1:naer) = 0.0_r8

      if (msoa_flag1 < 1000) then
      
      kappa_nonelectro(ibc_a  ) = 0.0001  
      kappa_nonelectro(ioc_a  ) = 0.0001  
      if (1 <= ilim2_a .and. ilim2_a <= naer) &
      kappa_nonelectro(ilim2_a) = 0.1     
      kappa_nonelectro(ioin_a ) = 0.06    

      else
      



      if (ioin_a        > 0) kappa_nonelectro(ioin_a       ) = hygro_oin_aer
      if (ioc_a         > 0) kappa_nonelectro(ioc_a        ) = hygro_oc_aer
      if (ibc_a         > 0) kappa_nonelectro(ibc_a        ) = hygro_bc_aer

      if (ipcg1_b_c_a   > 0) kappa_nonelectro(ipcg1_b_c_a  ) = hygro_pcg1_b_c_aer
      if (ipcg2_b_c_a   > 0) kappa_nonelectro(ipcg2_b_c_a  ) = hygro_pcg2_b_c_aer
      if (ipcg3_b_c_a   > 0) kappa_nonelectro(ipcg3_b_c_a  ) = hygro_pcg3_b_c_aer
      if (ipcg4_b_c_a   > 0) kappa_nonelectro(ipcg4_b_c_a  ) = hygro_pcg4_b_c_aer
      if (ipcg5_b_c_a   > 0) kappa_nonelectro(ipcg5_b_c_a  ) = hygro_pcg5_b_c_aer
      if (ipcg6_b_c_a   > 0) kappa_nonelectro(ipcg6_b_c_a  ) = hygro_pcg6_b_c_aer
      if (ipcg7_b_c_a   > 0) kappa_nonelectro(ipcg7_b_c_a  ) = hygro_pcg7_b_c_aer
      if (ipcg8_b_c_a   > 0) kappa_nonelectro(ipcg8_b_c_a  ) = hygro_pcg8_b_c_aer
      if (ipcg9_b_c_a   > 0) kappa_nonelectro(ipcg9_b_c_a  ) = hygro_pcg9_b_c_aer
      if (ipcg1_b_o_a   > 0) kappa_nonelectro(ipcg1_b_o_a  ) = hygro_pcg1_b_o_aer
      if (ipcg2_b_o_a   > 0) kappa_nonelectro(ipcg2_b_o_a  ) = hygro_pcg2_b_o_aer
      if (ipcg3_b_o_a   > 0) kappa_nonelectro(ipcg3_b_o_a  ) = hygro_pcg3_b_o_aer
      if (ipcg4_b_o_a   > 0) kappa_nonelectro(ipcg4_b_o_a  ) = hygro_pcg4_b_o_aer
      if (ipcg5_b_o_a   > 0) kappa_nonelectro(ipcg5_b_o_a  ) = hygro_pcg5_b_o_aer
      if (ipcg6_b_o_a   > 0) kappa_nonelectro(ipcg6_b_o_a  ) = hygro_pcg6_b_o_aer
      if (ipcg7_b_o_a   > 0) kappa_nonelectro(ipcg7_b_o_a  ) = hygro_pcg7_b_o_aer
      if (ipcg8_b_o_a   > 0) kappa_nonelectro(ipcg8_b_o_a  ) = hygro_pcg8_b_o_aer
      if (ipcg9_b_o_a   > 0) kappa_nonelectro(ipcg9_b_o_a  ) = hygro_pcg9_b_o_aer
      if (iopcg1_b_c_a  > 0) kappa_nonelectro(iopcg1_b_c_a ) = hygro_opcg1_b_c_aer
      if (iopcg2_b_c_a  > 0) kappa_nonelectro(iopcg2_b_c_a ) = hygro_opcg2_b_c_aer
      if (iopcg3_b_c_a  > 0) kappa_nonelectro(iopcg3_b_c_a ) = hygro_opcg3_b_c_aer
      if (iopcg4_b_c_a  > 0) kappa_nonelectro(iopcg4_b_c_a ) = hygro_opcg4_b_c_aer
      if (iopcg5_b_c_a  > 0) kappa_nonelectro(iopcg5_b_c_a ) = hygro_opcg5_b_c_aer
      if (iopcg6_b_c_a  > 0) kappa_nonelectro(iopcg6_b_c_a ) = hygro_opcg6_b_c_aer
      if (iopcg7_b_c_a  > 0) kappa_nonelectro(iopcg7_b_c_a ) = hygro_opcg7_b_c_aer
      if (iopcg8_b_c_a  > 0) kappa_nonelectro(iopcg8_b_c_a ) = hygro_opcg8_b_c_aer
      if (iopcg1_b_o_a  > 0) kappa_nonelectro(iopcg1_b_o_a ) = hygro_opcg1_b_o_aer
      if (iopcg2_b_o_a  > 0) kappa_nonelectro(iopcg2_b_o_a ) = hygro_opcg2_b_o_aer
      if (iopcg3_b_o_a  > 0) kappa_nonelectro(iopcg3_b_o_a ) = hygro_opcg3_b_o_aer
      if (iopcg4_b_o_a  > 0) kappa_nonelectro(iopcg4_b_o_a ) = hygro_opcg4_b_o_aer
      if (iopcg5_b_o_a  > 0) kappa_nonelectro(iopcg5_b_o_a ) = hygro_opcg5_b_o_aer
      if (iopcg6_b_o_a  > 0) kappa_nonelectro(iopcg6_b_o_a ) = hygro_opcg6_b_o_aer
      if (iopcg7_b_o_a  > 0) kappa_nonelectro(iopcg7_b_o_a ) = hygro_opcg7_b_o_aer
      if (iopcg8_b_o_a  > 0) kappa_nonelectro(iopcg8_b_o_a ) = hygro_opcg8_b_o_aer
      if (ipcg1_f_c_a   > 0) kappa_nonelectro(ipcg1_f_c_a  ) = hygro_pcg1_f_c_aer
      if (ipcg2_f_c_a   > 0) kappa_nonelectro(ipcg2_f_c_a  ) = hygro_pcg2_f_c_aer
      if (ipcg3_f_c_a   > 0) kappa_nonelectro(ipcg3_f_c_a  ) = hygro_pcg3_f_c_aer
      if (ipcg4_f_c_a   > 0) kappa_nonelectro(ipcg4_f_c_a  ) = hygro_pcg4_f_c_aer
      if (ipcg5_f_c_a   > 0) kappa_nonelectro(ipcg5_f_c_a  ) = hygro_pcg5_f_c_aer
      if (ipcg6_f_c_a   > 0) kappa_nonelectro(ipcg6_f_c_a  ) = hygro_pcg6_f_c_aer
      if (ipcg7_f_c_a   > 0) kappa_nonelectro(ipcg7_f_c_a  ) = hygro_pcg7_f_c_aer
      if (ipcg8_f_c_a   > 0) kappa_nonelectro(ipcg8_f_c_a  ) = hygro_pcg8_f_c_aer
      if (ipcg9_f_c_a   > 0) kappa_nonelectro(ipcg9_f_c_a  ) = hygro_pcg9_f_c_aer
      if (ipcg1_f_o_a   > 0) kappa_nonelectro(ipcg1_f_o_a  ) = hygro_pcg1_f_o_aer
      if (ipcg2_f_o_a   > 0) kappa_nonelectro(ipcg2_f_o_a  ) = hygro_pcg2_f_o_aer
      if (ipcg3_f_o_a   > 0) kappa_nonelectro(ipcg3_f_o_a  ) = hygro_pcg3_f_o_aer
      if (ipcg4_f_o_a   > 0) kappa_nonelectro(ipcg4_f_o_a  ) = hygro_pcg4_f_o_aer
      if (ipcg5_f_o_a   > 0) kappa_nonelectro(ipcg5_f_o_a  ) = hygro_pcg5_f_o_aer
      if (ipcg6_f_o_a   > 0) kappa_nonelectro(ipcg6_f_o_a  ) = hygro_pcg6_f_o_aer
      if (ipcg7_f_o_a   > 0) kappa_nonelectro(ipcg7_f_o_a  ) = hygro_pcg7_f_o_aer
      if (ipcg8_f_o_a   > 0) kappa_nonelectro(ipcg8_f_o_a  ) = hygro_pcg8_f_o_aer
      if (ipcg9_f_o_a   > 0) kappa_nonelectro(ipcg9_f_o_a  ) = hygro_pcg9_f_o_aer
      if (iopcg1_f_c_a  > 0) kappa_nonelectro(iopcg1_f_c_a ) = hygro_opcg1_f_c_aer
      if (iopcg2_f_c_a  > 0) kappa_nonelectro(iopcg2_f_c_a ) = hygro_opcg2_f_c_aer
      if (iopcg3_f_c_a  > 0) kappa_nonelectro(iopcg3_f_c_a ) = hygro_opcg3_f_c_aer
      if (iopcg4_f_c_a  > 0) kappa_nonelectro(iopcg4_f_c_a ) = hygro_opcg4_f_c_aer
      if (iopcg5_f_c_a  > 0) kappa_nonelectro(iopcg5_f_c_a ) = hygro_opcg5_f_c_aer
      if (iopcg6_f_c_a  > 0) kappa_nonelectro(iopcg6_f_c_a ) = hygro_opcg6_f_c_aer
      if (iopcg7_f_c_a  > 0) kappa_nonelectro(iopcg7_f_c_a ) = hygro_opcg7_f_c_aer
      if (iopcg8_f_c_a  > 0) kappa_nonelectro(iopcg8_f_c_a ) = hygro_opcg8_f_c_aer
      if (iopcg1_f_o_a  > 0) kappa_nonelectro(iopcg1_f_o_a ) = hygro_opcg1_f_o_aer
      if (iopcg2_f_o_a  > 0) kappa_nonelectro(iopcg2_f_o_a ) = hygro_opcg2_f_o_aer
      if (iopcg3_f_o_a  > 0) kappa_nonelectro(iopcg3_f_o_a ) = hygro_opcg3_f_o_aer
      if (iopcg4_f_o_a  > 0) kappa_nonelectro(iopcg4_f_o_a ) = hygro_opcg4_f_o_aer
      if (iopcg5_f_o_a  > 0) kappa_nonelectro(iopcg5_f_o_a ) = hygro_opcg5_f_o_aer
      if (iopcg6_f_o_a  > 0) kappa_nonelectro(iopcg6_f_o_a ) = hygro_opcg6_f_o_aer
      if (iopcg7_f_o_a  > 0) kappa_nonelectro(iopcg7_f_o_a ) = hygro_opcg7_f_o_aer
      if (iopcg8_f_o_a  > 0) kappa_nonelectro(iopcg8_f_o_a ) = hygro_opcg8_f_o_aer

      if (iant1_c_a     > 0) kappa_nonelectro(iant1_c_a    ) = hygro_ant1_c_aer
      if (iant2_c_a     > 0) kappa_nonelectro(iant2_c_a    ) = hygro_ant2_c_aer
      if (iant3_c_a     > 0) kappa_nonelectro(iant3_c_a    ) = hygro_ant3_c_aer
      if (iant4_c_a     > 0) kappa_nonelectro(iant4_c_a    ) = hygro_ant4_c_aer
      if (iant1_o_a     > 0) kappa_nonelectro(iant1_o_a    ) = hygro_ant1_o_aer
      if (iant2_o_a     > 0) kappa_nonelectro(iant2_o_a    ) = hygro_ant2_o_aer
      if (iant3_o_a     > 0) kappa_nonelectro(iant3_o_a    ) = hygro_ant3_o_aer
      if (iant4_o_a     > 0) kappa_nonelectro(iant4_o_a    ) = hygro_ant4_o_aer
      if (ibiog1_c_a    > 0) kappa_nonelectro(ibiog1_c_a   ) = hygro_biog1_c_aer
      if (ibiog2_c_a    > 0) kappa_nonelectro(ibiog2_c_a   ) = hygro_biog2_c_aer
      if (ibiog3_c_a    > 0) kappa_nonelectro(ibiog3_c_a   ) = hygro_biog3_c_aer
      if (ibiog4_c_a    > 0) kappa_nonelectro(ibiog4_c_a   ) = hygro_biog4_c_aer
      if (ibiog1_o_a    > 0) kappa_nonelectro(ibiog1_o_a   ) = hygro_biog1_o_aer
      if (ibiog2_o_a    > 0) kappa_nonelectro(ibiog2_o_a   ) = hygro_biog2_o_aer
      if (ibiog3_o_a    > 0) kappa_nonelectro(ibiog3_o_a   ) = hygro_biog3_o_aer
      if (ibiog4_o_a    > 0) kappa_nonelectro(ibiog4_o_a   ) = hygro_biog4_o_aer

      if (ismpa_a       > 0) kappa_nonelectro(ismpa_a      ) = hygro_smpa_aer
      if (ismpbb_a      > 0) kappa_nonelectro(ismpbb_a     ) = hygro_smpbb_aer

      end if 

      return
      end subroutine set_kappa_nonelectro


  
  
  
  
  
  
  
  
  subroutine map_mosaic_species_aerchem_box( imap, jaerosolstate,  &
       rbox, aer, gas, jhyst_leg, num_a, Dp_dry_a,                 &
       sigmag_a, water_a, water_a_hyst, cair_mol_m3                )

    use module_data_mosaic_kind, only: r8
    use module_data_mosaic_aero, only: &
         nbin_a_max, naer, ngas_aerchtot, jtotal,       & 
         nbin_a,                                        & 
         jhyst_lo, jhyst_up, jhyst_undefined,           &
         mhyst_method, mhyst_uporlo_waterhyst,          &
         mw_aer_mac,                                    &
         all_solid, all_liquid, mixed, no_aerosol
    use module_data_mosaic_asecthp, only: &
         rbox_aer_ptr, rbox_gas_ptr
    use module_data_mosaic_main, only: &
         m_partmc_mosaic, naer_tot, ngas_max, ntot_used, &
         avogad, mw_air, piover6


    
    integer, intent(in) :: imap
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate, jhyst_leg

    real(r8), intent(in) :: cair_mol_m3
    real(r8), intent(inout), dimension(ntot_used) :: rbox
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a, Dp_dry_a, sigmag_a, water_a, water_a_hyst
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer

    
    character(len=256) :: errmsg
    integer :: ibin, iaer, igas, l, noffset
    real(r8) :: conv_aer, conv_aerinv
    real(r8) :: conv_gas, conv_gasinv
    real(r8) :: conv_num, conv_numinv
    real(r8) :: conv_wat, conv_watinv
    real(r8) :: tmpa


    if ((imap < 0) .or. (imap > 1)) then
       write(errmsg,*)'*** map_mosaic_species_BOX fatal error - bad imap =', imap
       call wrf_error_fatal3("<stdin>",582,&
trim(adjustl(errmsg)))
    end if


    
    
    
    conv_gas = 1.e3_r8*cair_mol_m3
    conv_gasinv = 1.0_r8/conv_gas

    
    conv_aer = mw_air*cair_mol_m3
    conv_aerinv = 1.0_r8/conv_aer

    
    conv_wat = 1.e-12_r8*mw_air*cair_mol_m3
    conv_watinv = 1.0_r8/conv_wat

    
    conv_num = 1.e-9_r8*mw_air*cair_mol_m3
    conv_numinv = 1.0_r8/conv_num


    if (imap == 0) then    
       
       
       gas(1:ngas_aerchtot) = 0.0_r8
       aer(1:naer,3,1:nbin_a_max) = 0.0_r8
       num_a(1:nbin_a_max) = 0.0_r8
       water_a(1:nbin_a_max) = 0.0_r8
       water_a_hyst(1:nbin_a_max) = 0.0_r8

       
       do igas = 1, ngas_aerchtot
          l = rbox_gas_ptr( igas )
          if (l > 0) gas(igas) = rbox(l)*conv_gas
       end do


       
       
       if (m_partmc_mosaic <= 0) then
        
          
          do ibin = 1, nbin_a

             l = rbox_aer_ptr( -1, ibin )
             if (l > 0) &
             num_a(ibin)      = rbox(l)*conv_num    


             l = rbox_aer_ptr( -2, ibin )
             if (l > 0) &
             water_a(ibin)    = rbox(l)*conv_wat  


             if (mhyst_method == mhyst_uporlo_waterhyst) then
                
                l = rbox_aer_ptr( -3, ibin )
                if (l > 0) &
                water_a_hyst(ibin) = rbox(l)*conv_wat 

                
                jhyst_leg(ibin) = jhyst_undefined
             else
                
                
                water_a_hyst(ibin) = 0.0_r8
             end if

             do iaer = 1, naer
                
                
                l = rbox_aer_ptr( iaer, ibin )
                if (l > 0) &
                aer(iaer,jtotal,ibin) = rbox(l)*conv_aer/mw_aer_mac(iaer)

             enddo
             
          enddo
        endif

    else if (imap == 1) then
       
       

       do igas = 1, ngas_aerchtot
          l = rbox_gas_ptr( igas )
          if (l > 0) rbox(l) = gas(igas)*conv_gasinv
       end do

       
       
       if (m_partmc_mosaic <= 0) then
          
          
          do ibin = 1, nbin_a
             

             l = rbox_aer_ptr( -1, ibin )
             if (l > 0) &
             rbox(l)                   = num_a(ibin)*conv_numinv


             l = rbox_aer_ptr( -2, ibin )
             if (l > 0) &
             rbox(l)                   = water_a(ibin)*conv_watinv


             if (mhyst_method == mhyst_uporlo_waterhyst) then
                
                l = rbox_aer_ptr( -3, ibin )
                if ( jaerosolstate(ibin) == all_solid  .or. &
                     jaerosolstate(ibin) == all_liquid .or. &
                     jaerosolstate(ibin) == mixed      ) then
                   if (l > 0) &
                   rbox(l) = water_a_hyst(ibin)*conv_watinv

                else
                   if (l > 0) &
                   rbox(l) = 0.0_r8

                end if
                
             else
                
                
                if ( jaerosolstate(ibin) == all_solid  .or. &
                     jaerosolstate(ibin) == all_liquid .or. &
                     jaerosolstate(ibin) == mixed      ) then
                   continue

                else
                   jhyst_leg(ibin) = jhyst_undefined
                end if
             end if
             
             do iaer = 1, naer
                l = rbox_aer_ptr( iaer, ibin )
                if (l > 0) &
                rbox(l)                     = aer(iaer,jtotal,ibin)*conv_aerinv*mw_aer_mac(iaer)

             enddo
             
          enddo
       endif

    endif

    return
  end subroutine map_mosaic_species_aerchem_box



  end module module_mosaic_aerchem_intr
