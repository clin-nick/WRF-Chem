  module module_mosaic_aerdynam_intr


  implicit none


  contains


  



  subroutine aerosoldynamics(                     & 
     idiagbb_host,                                &
     id, ii, jj, kk,                              &
     ktau, ktauc, dtchem,                         &
     tempbox, presbox, airdenbox, relhumbox,      &
     swdownbox,                                   &
     rbox, iter_mesa, jaerosolstate,              & 
     dp_dry_a, dp_wet_a                           ) 






  use module_data_mosaic_kind, only: r8
  use module_data_mosaic_main, only: &
       mw_air, m_partmc_mosaic, ntot_max, ntot_used
  use module_data_mosaic_aero, only : &
       dens_aer_mac, mw_aer_mac, mw_comp_a, &
       jhyst_undefined, msectional, msize_framework, &
       naer, nbin_a, nbin_a_max, ngas_aerchtot, no_aerosol, nsalt

  use module_mosaic_aerchem_intr, only: aerchemistry
  use module_mosaic_sect_intr, only: sectional_interface_1


  
  integer,  intent(in)  :: idiagbb_host
  integer,  intent(in)  :: id, ii, jj, kk
  integer,  intent(in)  :: ktau, ktauc

  integer,  intent(inout), dimension(nbin_a_max) :: jaerosolstate

  integer,  intent(inout), dimension(nbin_a_max) :: iter_mesa

  real(r8), intent(in) :: dtchem       


  real(r8), intent(in) :: tempbox      
  real(r8), intent(in) :: presbox      
  real(r8), intent(in) :: airdenbox    
  real(r8), intent(in) :: relhumbox    
  real(r8), intent(in) :: swdownbox    


  real(r8), intent(inout), dimension(ntot_used)  :: rbox

  real(r8), intent(inout), dimension(nbin_a_max) :: dp_dry_a
  real(r8), intent(inout), dimension(nbin_a_max) :: dp_wet_a


  
  integer :: it_host, it_mosaic
  integer, dimension(6)          :: hostgridinfox
  integer, dimension(nbin_a_max) :: jaerosolstate_bgn
  integer, dimension(nbin_a_max) :: jhyst_leg, jhyst_leg_sv1

  logical :: ldiag1, ldiag2

  real(r8) :: cair_mol_cc, cair_mol_m3
  real(r8) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr
  real(r8) :: pr_atm
  real(r8) :: rh
  real(r8) :: te

  real(r8), dimension(ngas_aerchtot)     :: gas_avg  
  real(r8), dimension(ngas_aerchtot)     :: gas_netprod_otrproc
            
            
            
            
  real(r8), dimension(nbin_a_max)        :: aH2O_a
  real(r8), dimension(nbin_a_max)        :: gam_ratio
  real(r8), dimension(nbin_a_max)        :: dens_dry_a_bgn, dens_dry_a

  real(r8), dimension(nbin_a_max)        :: dp_dry_a_sv1
  real(r8), dimension(nbin_a_max)        :: mass_dry_a_bgn, mass_dry_a
  real(r8), dimension(nbin_a_max)        :: sigmag_a, sigmag_a_sv1

  real(r8), allocatable, dimension(:)    :: rbox_sv1


  if (idiagbb_host >= 100) then
     ldiag1 = .true. ; ldiag2 = .true.
  else
     ldiag1 = .false. ; ldiag2 = .false.
  end if

  if ( ldiag2 ) &
  call dump_aerdy( 'dya',                         & 
     id, ii, jj, kk,                              &
     ktau, ktauc, dtchem,                         &
     tempbox, presbox, airdenbox, relhumbox,      &
     rbox, iter_mesa, jaerosolstate,              & 
     dp_dry_a, dp_wet_a                           ) 

  if ( ldiag1 ) write(*,*) 'aerosoldynamics start'
  if ( (m_partmc_mosaic <= 0) .and. &
       (msize_framework == msectional) ) then
     if ( .not. allocated( rbox_sv1 ) ) allocate( rbox_sv1(1:ntot_used) )

  end if



  it_host     = ktau
  it_mosaic   = ktauc
  hostgridinfox(1) = id
  hostgridinfox(2) = ii
  hostgridinfox(3) = jj
  hostgridinfox(4) = kk
  hostgridinfox(5:6) = 0

  te = tempbox
  rh = relhumbox*100.0_r8  

  pr_atm = presbox / 1.01325e5_r8          
  pr_atm = presbox / 1.032e5_r8            

  cair_mol_m3 = airdenbox*1.0e3_r8/mw_air  
  cair_mol_cc = cair_mol_m3*1.e-6_r8       

  iter_mesa(1:nbin_a_max) = 0
  jaerosolstate(1:nbin_a_max) = no_aerosol
  jaerosolstate_bgn(1:nbin_a_max) = no_aerosol
  jhyst_leg(1:nbin_a_max) = jhyst_undefined
  jhyst_leg_sv1(1:nbin_a_max) = jhyst_undefined

  gas_avg(1:ngas_aerchtot) = 0.0_r8
  gas_netprod_otrproc(1:ngas_aerchtot) = 0.0_r8

  ah2o_a(1:nbin_a_max) = 0.0_r8
  gam_ratio(1:nbin_a_max) = 0.0_r8
  dens_dry_a(1:nbin_a_max) = 0.0_r8
  dens_dry_a_bgn(1:nbin_a_max) = 0.0_r8
  dp_wet_a(1:nbin_a_max) = 0.0_r8
  dp_dry_a(1:nbin_a_max) = 0.0_r8
  dp_dry_a_sv1(1:nbin_a_max) = 0.0_r8
  mass_dry_a(1:nbin_a_max) = 0.0_r8
  mass_dry_a_bgn(1:nbin_a_max) = 0.0_r8
  sigmag_a(1:nbin_a_max) = 0.0_r8
  sigmag_a_sv1(1:nbin_a_max) = 0.0_r8


  
  if ( ldiag1 ) write(*,*) 'aerosoldynamics call map 0'
  call aerodynam_map_mosaic_species( 0, &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       sigmag_a )





  if ( (m_partmc_mosaic <= 0) .and. &
       (msize_framework == msectional) ) then
     rbox_sv1(1:ntot_used) = rbox(1:ntot_used)
     dp_dry_a_sv1(1:nbin_a) = dp_dry_a(1:nbin_a)
     sigmag_a_sv1(1:nbin_a) = sigmag_a(1:nbin_a)
     jhyst_leg_sv1(1:nbin_a) = jhyst_leg(1:nbin_a)
  end if


  if ( ldiag2 ) &
  call dump_aerdy( 'dyc',                         & 
     id, ii, jj, kk,                              &
     ktau, ktauc, dtchem,                         &
     tempbox, presbox, airdenbox, relhumbox,      &
     rbox, iter_mesa, jaerosolstate,              & 
     dp_dry_a, dp_wet_a                           ) 

  if ( ldiag1 ) write(*,*) 'aerosoldynamics call aerchem'
  call aerchemistry(                                          &
     idiagbb_host,                                            &
     hostgridinfox, it_host, it_mosaic, dtchem,               & 
     pr_atm, rh, te, cair_mol_m3, cair_mol_cc, swdownbox,     &
     jaerosolstate, jaerosolstate_bgn, jhyst_leg,             & 
     rbox, dp_dry_a, dp_wet_a, sigmag_a,                      & 
     gas_avg, gas_netprod_otrproc,                            & 
     mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a,  &
     aH2O_a, gam_ratio, iter_mesa                             ) 


  if ( ldiag2 ) &
  call dump_aerdy( 'dye',                         & 
     id, ii, jj, kk,                              &
     ktau, ktauc, dtchem,                         &
     tempbox, presbox, airdenbox, relhumbox,      &
     rbox, iter_mesa, jaerosolstate,              & 
     dp_dry_a, dp_wet_a                           ) 

  
  
  
  
  if ( (m_partmc_mosaic <= 0) .and. &
       (msize_framework == msectional) ) then

      if ( ldiag1 ) write(*,*) 'aerosoldynamics call sectiface'
      call sectional_interface_1( &
          idiagbb_host, &
          id, ii, jj, kk, ktau, ktauc, dtchem, &
          rbox_sv1, rbox, &
          it_mosaic, jaerosolstate_bgn, jaerosolstate, &
          jhyst_leg, jhyst_leg_sv1, &
          dp_dry_a, dp_dry_a_sv1, &
          sigmag_a, sigmag_a_sv1, &
          mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
          fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
          pr_atm, rh, te, cair_mol_m3, cair_mol_cc )



  
  
  end if


  
  






  if ( ldiag1 ) write(*,*) 'aerosoldynamics end'


  if ( ldiag2 ) &
  call dump_aerdy( 'dyz',                         & 
     id, ii, jj, kk,                              &
     ktau, ktauc, dtchem,                         &
     tempbox, presbox, airdenbox, relhumbox,      &
     rbox, iter_mesa, jaerosolstate,              & 
     dp_dry_a, dp_wet_a                           ) 

  return
  end subroutine aerosoldynamics



  
  subroutine aerodynam_map_mosaic_species( imap, &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       sigmag_a )




    
    
    
    
    

    use module_data_mosaic_kind




    use module_data_mosaic_aero, only: &
       nbin_a, nbin_a_max







    use module_data_mosaic_asecthp, only: &
       dcen_sect, sigmag_aer, isize_of_ibin, itype_of_ibin







    
    integer, intent(in) :: imap




    real(r8), intent(inout) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr


    real(r8), intent(inout), dimension(nbin_a_max) :: sigmag_a


    
    integer :: ibin, iphase, isize, itype



    character(len=256) :: errmsg

    if ((imap < 0) .or. (imap > 1)) then
       write(errmsg,*) &
       '*** aerchem_map_mosaic_cnn_rbox fatal error - bad imap =', imap
       call wrf_error_fatal3("<stdin>",315,&
trim(adjustl(errmsg)))
    end if


    
    
    
    
    
    
    
    
    fact_gasmr = 1.0e-6_r8

    
    fact_apmassmr = 1.0e-9_r8

    
    fact_apnumbmr = 1.0e-3_r8

    
    
    fact_apdens = 1.0_r8

    
    
    fact_apdiam = 1.0_r8


    

    
    

    
    

    
    

    
    

    
    



    if (imap == 1) goto 20000


    
    
    






    do ibin = 1, nbin_a
       isize = isize_of_ibin(ibin)
       itype = itype_of_ibin(ibin)


       
       
       


       sigmag_a(ibin)          = sigmag_aer(isize,itype)









          
          
          



          
          
          
          











    end do

    return


20000 continue
    
    
    




















          
          
          


          
          
          
          





    return
  end subroutine aerodynam_map_mosaic_species



  
  subroutine dump_aerdy( dumpchaa,                & 
     id, ii, jj, kk,                              &
     ktau, ktauc, dtchem,                         &
     tempbox, presbox, airdenbox, relhumbox,      &
     rbox, iter_mesa, jaerosolstate,              & 
     dp_dry_a, dp_wet_a                           ) 

  use module_data_mosaic_kind, only: r8
  use module_data_mosaic_main, only: &
       ntot_used
  use module_data_mosaic_aero, only : &
       nbin_a_max
  use module_data_mosaic_asecthp, only : &
       nsize_aer, ntype_aer, numptr_aer, lptr_bc_aer, lptr_oc_aer, lptr_oin_aer

  use module_mosaic_aerchem_intr, only: aerchemistry
  use module_mosaic_sect_intr, only: sectional_interface_1


  
  character(len=*), intent(in) :: dumpchaa

  integer,  intent(in)  :: id, ii, jj, kk
  integer,  intent(in)  :: ktau, ktauc
  integer,  intent(inout), dimension(nbin_a_max) :: jaerosolstate
  integer,  intent(inout), dimension(nbin_a_max) :: iter_mesa

  real(r8), intent(in) :: dtchem       
  real(r8), intent(in) :: tempbox      
  real(r8), intent(in) :: presbox      
  real(r8), intent(in) :: airdenbox    
  real(r8), intent(in) :: relhumbox    

  real(r8), intent(inout), dimension(ntot_used)  :: rbox

  real(r8), intent(inout), dimension(nbin_a_max) :: dp_dry_a
  real(r8), intent(inout), dimension(nbin_a_max) :: dp_wet_a

  
  integer :: isize, itype
  real(r8) :: tmpa
  character(len=80) :: fmtaa

  if (ii*jj*kk /= 1) return

  fmtaa = '(2a,1p,e10.3,7e11.3)'
  if (nsize_aer(1)*ntype_aer > 8) fmtaa = '(2a,1p,e10.3,7e11.3/(7x,1p,e10.3,7e11.3))'


  tmpa = airdenbox*1.0e-6_r8
  write(131,fmtaa) dumpchaa, ' num', &
     ( ( rbox(numptr_aer(isize,itype,1))*tmpa, &
         isize=1,nsize_aer(itype) ), itype=1,ntype_aer )



  tmpa = airdenbox*1.0e3_r8
  write(131,fmtaa) dumpchaa, ' oin', &
     ( ( rbox(lptr_oin_aer(isize,itype,1))*tmpa, &
         isize=1,nsize_aer(itype) ), itype=1,ntype_aer )
  write(131,fmtaa) dumpchaa, ' oc ', &
     ( ( rbox(lptr_oc_aer(isize,itype,1))*tmpa, &
         isize=1,nsize_aer(itype) ), itype=1,ntype_aer )
  write(131,fmtaa) dumpchaa, ' bc ', &
     ( ( rbox(lptr_bc_aer(isize,itype,1))*tmpa, &
         isize=1,nsize_aer(itype) ), itype=1,ntype_aer )

  return
  end subroutine dump_aerdy



  end module module_mosaic_aerdynam_intr
