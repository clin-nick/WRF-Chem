



















module module_cam_mam_init
  
  private
  public :: cam_mam_init
  
  
  
  
  
  
  
  
  
  integer, parameter :: pom_emit_1p4_factor_flag = 0
  integer, parameter :: pom_icbc_1p4_factor_flag = 0
  real, public :: pom_emit_1p4_factor
  real, public :: pom_icbc_1p4_factor
  
  
  
  
  
  
  
  
  
  integer, parameter :: so4_emit_1p2_factor_flag = 1
  integer, parameter :: so4_icbc_1p2_factor_flag = 1
  integer, parameter :: so4_mwgt_1p2_factor_flag = 1
  real, public :: so4_emit_1p2_factor
  real, public :: so4_icbc_1p2_factor
  
  
  

  
  character(len=16),allocatable, public :: cnst_name_loc(:)     

  integer, parameter :: init_val = -999888777
  LOGICAL :: CAM_INITIALIZED_CHEM = .FALSE.
contains
  
  
  
  subroutine cam_mam_init(           &
       id, numgas, config_flags,       &
       ids,ide, jds,jde, kds,kde,      &
       ims,ime, jms,jme, kms,kme,      &
       its,ite, jts,jte, kts,kte       )
    
    
    
    
    
    
    
    
    use module_state_description, only: num_chem,CBMZ_CAM_MAM3_NOAQ,CBMZ_CAM_MAM3_AQ,CBMZ_CAM_MAM7_NOAQ,CBMZ_CAM_MAM7_AQ
    use module_configure, only:  grid_config_rec_type
    
    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst, only: epsilo, latvap, latice, rh2o, cpair, tmelt
    use module_cam_esinti, only: esinti
    
    use module_cam_support, only: pver, pverp, pcols, &
         pcnst => pcnst_runtime, &
         endrun, masterproc
    
    use modal_aero_data, only:  cnst_name_cw, ntot_aspectype, &
         qneg3_worst_thresh_amode, species_class, &
         specmw_amode, specname_amode, &
         specdens_so4_amode, specmw_so4_amode, &
         specdens_nh4_amode, specmw_nh4_amode, &
         specdens_no3_amode, specmw_no3_amode, &
         specdens_pom_amode, specmw_pom_amode, &
         specdens_soa_amode, specmw_soa_amode, &
         specdens_bc_amode, specmw_bc_amode, &
         specdens_dust_amode, specmw_dust_amode, &
         specdens_seasalt_amode, specmw_seasalt_amode
    
    use modal_aero_initialize_data, only: modal_aero_initialize, modal_aero_register, decouple_mam_mp
    use ndrop , only: activate_init
    USE module_cam_mam_cloudchem, only: cam_mam_cloudchem_inti
    USE module_cam_mam_gas_wetdep_driver, only: cam_mam_gas_wetdep_inti
    
    implicit none
    
    
    
    
    type(grid_config_rec_type), intent(in) :: config_flags
    
    integer, intent(in) ::   &
         id, numgas,   &
         ids, ide, jds, jde, kds, kde,   &
         ims, ime, jms, jme, kms, kme,   &
         its, ite, jts, jte, kts, kte
    
    
    
    
    
    integer :: ierr, l, m
    character(len=16)  :: tmpname
    character(len=160) :: msg
    
    
    
    if ( (config_flags%chem_opt /= CBMZ_CAM_MAM3_NOAQ) .and. &
         (config_flags%chem_opt /= CBMZ_CAM_MAM3_AQ  ) ) then
       call wrf_error_fatal3("<stdin>",131,&
'cam_mam_init - MODAL_AERO_3MODE is defined but chem_opt is not a CAM_MAM3 package' )
    end if
    
    
    
    if((config_flags%chem_opt == CBMZ_CAM_MAM3_NOAQ .OR. config_flags%chem_opt == CBMZ_CAM_MAM3_AQ) .AND.  config_flags%cam_mam_mode .NE. 3)then
       call wrf_error_fatal3("<stdin>",138,&
'CAM_MAM_INIT - For MODAL_AERO_3MODE (chem_opt - 503 CAM_MAM3 package), cam_mam_mode in namelist should be set to 3' )
    elseif((config_flags%chem_opt == CBMZ_CAM_MAM7_NOAQ .OR. config_flags%chem_opt == CBMZ_CAM_MAM7_AQ) .AND.  config_flags%cam_mam_mode .NE. 7)then
       call wrf_error_fatal3("<stdin>",141,&
'CAM_MAM_INIT - For MODAL_AERO_7MODE (chem_opt - 504 CAM_MAM7 package), cam_mam_mode in namelist should be set to 7')
    endif

    
    
    if((config_flags%chem_opt == CBMZ_CAM_MAM3_NOAQ .OR. config_flags%chem_opt == CBMZ_CAM_MAM3_AQ) .AND.  config_flags%cam_mam_nspec .NE. 85)then
       call wrf_error_fatal3("<stdin>",148,&
'CAM_MAM_INIT - For MODAL_AERO_3MODE (chem_opt - 503 CAM_MAM3 package), cam_mam_nspec in namelist should be set to 85' )
    elseif((config_flags%chem_opt == CBMZ_CAM_MAM7_NOAQ .OR. config_flags%chem_opt == CBMZ_CAM_MAM7_AQ) .AND.  config_flags%cam_mam_nspec .NE. 90)then
       
       call wrf_error_fatal3("<stdin>",152,&
'CAM_MAM_INIT - For MODAL_AERO_7MODE (chem_opt - 504 CAM_MAM7 package), cam_mam_nspec in namelist should be set to 90')
    endif
    
    
    
    
    write(*,'(/a)') 'cam_mam_init'
    write(*,*) 'id, num_chem, pcnst =', id, num_chem, pcnst
    
    write(*,'(a,3(4x,2i5))') 'ids/e, j..., k... ', ids,ide, jds,jde, kds,kde
    write(*,'(a,3(4x,2i5))') 'ims/e, j..., k... ', ims,ime, jms,jme, kms,kme
    write(*,'(a,3(4x,2i5))') 'its/e, j..., k... ', its,ite, jts,jte, kts,kte
    write(*,'(a,3(   i14))') 'pver, pverp, pcols', pver, pverp, pcols
    
    
    
    call esinti(epsilo, latvap, latice, rh2o, cpair, tmelt)
    
    
    
    if ( (max(pver,pverp) > 0) .and. &
         (pver /= kte-kts+1) .and. &
         (pverp /= pver+1) ) then
       write( msg, '(2a,3i15)' ) &
            'cam_mam_init fatal error ', &
            '- bad pver - id, pver, pverp = ', id, pver, pverp
       call wrf_error_fatal3("<stdin>",179,&
msg )
    end if
    if (pver <= 0) then
       pver = kde - kds
       pverp = pver + 1
    end if
    write(*,'(a,3(   i14))') 'pver, pverp, pcols', pver, pverp, pcols
    
    
    
    write(*,'(/a)') &
         'cam_mam_init calling cam_mam_init_set_cnst'
    call cam_mam_init_set_cnst( id, numgas, config_flags )
    


    write(*,'(/a)') &
         'cam_mam_init calling modal_aero_register'
    call modal_aero_register
    
    
    write(*,'(/a)') &
         'cam_mam_init calling modal_aero_initialize'
    call modal_aero_initialize

    

    call decouple_mam_mp(config_flags%CAM_MP_MAM_cpled)

    
    call activate_init
    
    if (so4_mwgt_1p2_factor_flag <= 0) then
       
       
       do m = 1, ntot_aspectype
          if      (specname_amode(m).eq.'sulfate   ') then

             specmw_so4_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'ammonium  ') then

             specmw_nh4_amode = specmw_amode(m)
          end if
       end do
    end if
    
    write(*,'(a,2f12.4)') 'so4 dens, mw', specdens_so4_amode, specmw_so4_amode
    write(*,'(a,2f12.4)') 'nh4 dens, mw', specdens_nh4_amode, specmw_nh4_amode
    write(*,'(a,2f12.4)') 'no3 dens, mw', specdens_no3_amode, specmw_no3_amode
    write(*,'(a,2f12.4)') 'pom dens, mw', specdens_pom_amode, specmw_pom_amode
    write(*,'(a,2f12.4)') 'soa dens, mw', specdens_soa_amode, specmw_soa_amode
    write(*,'(a,2f12.4)') 'bc  dens, mw', specdens_bc_amode, specmw_bc_amode
    write(*,'(a,2f12.4)') 'dst dens, mw', specdens_dust_amode, specmw_dust_amode
    write(*,'(a,2f12.4)') 'ncl dens, mw', specdens_seasalt_amode, specmw_seasalt_amode
    
    
    
    
    pom_emit_1p4_factor = 1.0
    pom_icbc_1p4_factor = 1.0
    if (pom_emit_1p4_factor_flag > 0) pom_emit_1p4_factor = 1.4
    if (pom_icbc_1p4_factor_flag > 0) pom_icbc_1p4_factor = 1.4
    
    so4_emit_1p2_factor = 1.0
    so4_icbc_1p2_factor = 1.0
    if (so4_emit_1p2_factor_flag > 0) so4_emit_1p2_factor = 115.0/96.0
    if (so4_icbc_1p2_factor_flag > 0) so4_icbc_1p2_factor = 115.0/96.0
    
    write(*,'(/a,2f10.4)') 'pom_emit_1p4_factor & _init_', &
         pom_emit_1p4_factor, pom_icbc_1p4_factor
    write(*,'( a,2f10.4)') 'so4_emit_1p2_factor & _init_', &
         so4_emit_1p2_factor, so4_icbc_1p2_factor
    
    
    
    write(*,'(/a)') &
         'cam_mam_init calling cam_mam_init_asect'
    call cam_mam_init_asect( id, config_flags )
    
    
    
    
    write(*,'(/a)') &
         'cam_mam_init calling cam_mam_init_other'
    call cam_mam_init_other( id, numgas, config_flags )

    
    call cam_mam_cloudchem_inti()

    
    call cam_mam_gas_wetdep_inti

    deallocate(cnst_name_loc)  
    
    write(*,'(/a)') &
         'cam_mam_init done'
    
    
    return
  end subroutine cam_mam_init
  
  
  
  subroutine cam_mam_init_set_cnst( id, numgas, config_flags )
    
    
    
    
    

    
    
    
    
    
    
    
    
    

    
    
    

    
    
    use module_configure, only:  grid_config_rec_type
    use module_state_description, only:  num_chem, param_first_scalar
    use module_scalar_tables, only:  chem_dname_table
    
    use module_cam_support,    only: pcnst => pcnst_runtime, &
         pcnst_non_chem => pcnst_non_chem_modal_aero, &
         gas_pcnst => gas_pcnst_modal_aero, &
         endrun, masterproc
    use modal_aero_data
    use constituents,             only: cnst_add
    use physconst,                only: cpair
    
    implicit none
    
    
    
    
    type(grid_config_rec_type), intent(in) :: config_flags
    integer, intent(in) :: id   
    integer, intent(in) :: numgas
    
    
    
    
    integer :: ierr, itmpa, dumind
    integer :: l, l2
    integer :: p1st
    character(len=360) :: msg
    
    
    
    
    write(*,*) 'cam_mam_init_set_cnst'
    write(*,*) 'id, num_chem         =', id, num_chem
    write(*,*) 'pcnst, gas_pcnst old =', pcnst, gas_pcnst
    
    
    p1st = param_first_scalar
    pcnst_non_chem = 5
    
    
    itmpa = pcnst_non_chem
    
    do l = p1st, numgas
       itmpa = itmpa + 1
    end do
    
    
    
    
    do l = numgas+1, num_chem
       if ( cam_mam_is_q_aerosol_species( l, numgas ) ) then
          itmpa = itmpa + 1
       end if
    end do
    
    if (pcnst == itmpa) then 
       
       gas_pcnst = pcnst - pcnst_non_chem
       write(*,*) 'pcnst, gas_pcnst new =', pcnst, gas_pcnst
    else if (pcnst /= itmpa) then
       write( msg, * ) &
            'CAM_MAM_INIT_SET_CNST fatal error: The value of PCNST should be set to:',itmpa, &
            ' in module_physics_init.F where pcnst is mentioned as', pcnst,'. ID is',id
       call wrf_error_fatal3("<stdin>",370,&
msg )
    end if
    
    
    if ( .not. allocated(cnst_name_loc) ) then
       allocate( cnst_name_loc(pcnst) )
    end if
    
    
    
    do l = 1, pcnst
       write( cnst_name_loc(l), '(a,i4.4)' ) 'empty_cnst_', l
    end do
    
    
    if (pcnst_non_chem >= 1) cnst_name_loc(1) = 'Q'
    if (pcnst_non_chem >= 2) cnst_name_loc(2) = 'CLDLIQ'
    if (pcnst_non_chem >= 3) cnst_name_loc(3) = 'CLDICE'
    if (pcnst_non_chem >= 4) cnst_name_loc(4) = 'NUMLIQ'
    if (pcnst_non_chem >= 5) cnst_name_loc(5) = 'NUMICE'
    
    l2 = pcnst_non_chem
    do l = p1st, numgas
       l2 = l2 + 1
       IF(.NOT.CAM_INITIALIZED_CHEM) call cnst_add(trim(adjustl(chem_dname_table(1,l))), 1.0_r8, cpair, 0._r8, dumind)   
       cnst_name_loc(l2) = chem_dname_table(1,l)
    end do
    
    do l = numgas+1, num_chem
       if ( cam_mam_is_q_aerosol_species( l, numgas ) ) then
          l2 = l2 + 1
          IF(.NOT.CAM_INITIALIZED_CHEM) call cnst_add(trim(adjustl(chem_dname_table(1,l))), 1.0_r8, cpair, 0._r8, dumind)      
          cnst_name_loc(l2) = chem_dname_table(1,l)
       end if
    end do
    if (l2 /= pcnst) then
       write( msg, '(2a,3i15)' ) &
            'cam_mam_init_set_cnst fatal error 101', &
            'for cnst_name, id, l2, pcnst = ', id, l2, pcnst
       call wrf_error_fatal3("<stdin>",410,&
msg )
    end if
    
    
    return
  end subroutine cam_mam_init_set_cnst
  
  
  
  function cam_mam_is_q_aerosol_species( lspec, numgas )
    
    
    
    
    
    
    
    
    use module_state_description, only:  num_chem, param_first_scalar
    use module_scalar_tables, only:  chem_dname_table
    use modal_aero_data, only:  ntot_amode
    
    implicit none
    
    logical ::  cam_mam_is_q_aerosol_species
    integer, intent(in) :: lspec, numgas
    
    integer :: i, n
    integer, parameter :: upper_to_lower = iachar('a')-iachar('A')
    character(len=32) :: tmpname
    character(len=1)  :: tmpch
    
    cam_mam_is_q_aerosol_species = .false.
    if ( (lspec <  param_first_scalar) .or. &
         (lspec <= numgas) .or. &
         (lspec >  num_chem) ) return
    
    call lower_case( chem_dname_table(1,lspec), tmpname )
    
    
    if (tmpname(1:5) == 'wtr_a') return
    
    n = len( trim(tmpname) )
    n = max( n, 4 )
    if (tmpname(n-2:n) == '_a1') cam_mam_is_q_aerosol_species = .true.
    if (tmpname(n-2:n) == '_a2') cam_mam_is_q_aerosol_species = .true.
    if (tmpname(n-2:n) == '_a3') cam_mam_is_q_aerosol_species = .true.
    if (ntot_amode == 3) return
    
    if (tmpname(n-2:n) == '_a4') cam_mam_is_q_aerosol_species = .true.
    if (tmpname(n-2:n) == '_a5') cam_mam_is_q_aerosol_species = .true.
    if (tmpname(n-2:n) == '_a6') cam_mam_is_q_aerosol_species = .true.
    if (tmpname(n-2:n) == '_a7') cam_mam_is_q_aerosol_species = .true.
    if (ntot_amode == 7) return
    
    return
  end function cam_mam_is_q_aerosol_species
  
  
  
  subroutine lower_case( txt_in, txt_lc )
    
    
    
    implicit none
    
    character(len=*), intent(in)  :: txt_in
    character(len=*), intent(out) :: txt_lc
    
    integer :: i, j
    integer, parameter :: iachar_lowera = iachar('a')
    integer, parameter :: iachar_uppera = iachar('A')
    integer, parameter :: iachar_upperz = iachar('Z')
    
    txt_lc = txt_in
    do i = 1, len( trim(txt_lc) )
       j = iachar( txt_lc(i:i) )
       if (j < iachar_uppera) cycle
       if (j > iachar_upperz) cycle
       txt_lc(i:i) = achar( j + iachar_lowera - iachar_uppera )
    end do
    
    return
  end subroutine lower_case
  
  
  
  subroutine cam_mam_init_asect( id, config_flags )
    
    
    
    
    
    
    
    
    
    
    
    
    use module_state_description, only: num_chem,param_first_scalar,CBMZ_CAM_MAM3_AQ,CBMZ_CAM_MAM7_AQ

    use module_configure, only:  grid_config_rec_type
    use module_scalar_tables, only:     chem_dname_table
    
    use modal_aero_data, only:  &
         dgnum_amode, lspectype_amode, nspec_amode, ntot_amode, ntot_aspectype, &
         sigmag_amode, specdens_amode, spechygro, specmw_amode, specname_amode
    
    use module_data_cam_mam_asect
    
    
    implicit none
    
    
    integer, intent(in) :: id
    type(grid_config_rec_type), intent(in) :: config_flags
    
    
    integer :: iphase, isize, itype
    integer :: l, l2, l3, l4, la, lc
    integer :: p1st
    
    real, parameter :: pi=3.1415926536
    real :: dp_meanvol_tmp
    
    character(len=160) :: msg
    character(len=32)  :: tmpname, tmpnamec, tmptxtaa
    
    
    p1st = param_first_scalar
    
    
    ntot_mastercomp_aer = ntot_aspectype
    do l = 1, ntot_mastercomp_aer
       name_mastercomp_aer(l) = specname_amode(l)
       
       
       dens_mastercomp_aer(l) = specdens_amode(l)*1.0e-3
       mw_mastercomp_aer(l) = specmw_amode(l)
       hygro_mastercomp_aer(l) = spechygro(l)
       if (name_mastercomp_aer(l) == 'sulfate') then
          mastercompindx_so4_aer = l
          namebb_mastercomp_aer(l) = 'so4'
          mw_so4_aer = mw_mastercomp_aer(l)
          dens_so4_aer = dens_mastercomp_aer(l)
       else if (name_mastercomp_aer(l) == 'ammonium') then
          mastercompindx_nh4_aer = l
          namebb_mastercomp_aer(l) = 'nh4'
          mw_nh4_aer = mw_mastercomp_aer(l)
          dens_nh4_aer = dens_mastercomp_aer(l)
       else if (name_mastercomp_aer(l) == 'nitrate') then
          mastercompindx_no3_aer = l
          namebb_mastercomp_aer(l) = 'no3'
          mw_no3_aer = mw_mastercomp_aer(l)
          dens_no3_aer = dens_mastercomp_aer(l)
       else if (name_mastercomp_aer(l) == 'p-organic') then
          mastercompindx_pom_aer = l
          namebb_mastercomp_aer(l) = 'pom'
          mw_pom_aer = mw_mastercomp_aer(l)
          dens_pom_aer = dens_mastercomp_aer(l)
       else if (name_mastercomp_aer(l) == 's-organic') then
          mastercompindx_soa_aer = l
          namebb_mastercomp_aer(l) = 'soa'
          mw_soa_aer = mw_mastercomp_aer(l)
          dens_soa_aer = dens_mastercomp_aer(l)
       else if (name_mastercomp_aer(l) == 'black-c') then
          mastercompindx_bc_aer = l
          namebb_mastercomp_aer(l) = 'bc'
          mw_bc_aer = mw_mastercomp_aer(l)
          dens_bc_aer = dens_mastercomp_aer(l)
       else if (name_mastercomp_aer(l) == 'seasalt') then
          mastercompindx_seas_aer = l
          namebb_mastercomp_aer(l) = 'ncl'
          mw_seas_aer = mw_mastercomp_aer(l)
          dens_seas_aer = dens_mastercomp_aer(l)
       else if (name_mastercomp_aer(l) == 'dust') then
          mastercompindx_dust_aer = l
          namebb_mastercomp_aer(l) = 'dst'
          mw_dust_aer = mw_mastercomp_aer(l)
          dens_dust_aer = dens_mastercomp_aer(l)
       else
          msg = '*** cam_mam_init_asect error 100 - mastercompindx'
          call wrf_message( msg )
          write( msg, '(a,i4,2x,a)' ) 'l, specname_amode = ', &
               l, specname_amode(l)
          call wrf_error_fatal3("<stdin>",597,&
msg )
       end if
    end do 
    
    
    
    nphase_aer = 1
    ai_phase = 1

          if ((config_flags%chem_opt == CBMZ_CAM_MAM3_AQ) .or. &
              (config_flags%chem_opt == CBMZ_CAM_MAM7_AQ  )) then
               nphase_aer = 2
               cw_phase = 2
         end if
    
    if (nphase_aer > 2) then
       msg = '*** cam_mam_init_asect error 120 - nphase_aer > 2'
       call wrf_error_fatal3("<stdin>",615,&
msg )
    end if
    
    
    nsize_aer(:) = 0
    ncomp_aer(:) = 0
    ncomp_plustracer_aer(:) = 0
    mastercompptr_aer(:,:)  = init_val
    massptr_aer(:,:,:,:)    = init_val
    waterptr_aer(:,:)       = init_val
    hyswptr_aer(:,:)        = init_val
    numptr_aer(:,:,:)       = init_val
    mprognum_aer(:,:,:)     = 0
    
    lptr_so4_aer(:,:,:)  = init_val
    lptr_nh4_aer(:,:,:)  = init_val
    lptr_no3_aer(:,:,:)  = init_val
    lptr_pom_aer(:,:,:)  = init_val
    lptr_soa_aer(:,:,:)  = init_val
    lptr_bc_aer(:,:,:)   = init_val
    lptr_dust_aer(:,:,:) = init_val
    lptr_seas_aer(:,:,:) = init_val
    
    volumcen_sect(:,:) = 0.0
    volumlo_sect(:,:) = 0.0
    volumhi_sect(:,:) = 0.0
    dcen_sect(:,:) = 0.0
    dlo_sect(:,:) = 0.0
    dhi_sect(:,:) = 0.0
    sigmag_aer(:,:) = 1.0
    
    
    
    
    
    
    
    
    ntype_aer = ntot_amode
    
    do itype = 1, ntype_aer
       
       nsize_aer(itype) = 1
       ncomp_aer(itype) = nspec_amode(itype)
       ncomp_plustracer_aer(itype) = ncomp_aer(itype)
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       isize = 1
       sigmag_aer(isize,itype) = sigmag_amode(itype)
       dp_meanvol_tmp = 1.0e2*dgnum_amode(itype) &   
            *exp( 1.5 * log(sigmag_aer(isize,itype))**2 )
       dcen_sect(isize,itype) = dp_meanvol_tmp
       dhi_sect( isize,itype) = dp_meanvol_tmp*4.0
       dlo_sect( isize,itype) = dp_meanvol_tmp/4.0
       
       do isize = 1, nsize_aer(itype)
          volumcen_sect(isize,itype) = (pi/6.0)*(dcen_sect(isize,itype)**3)
          volumlo_sect( isize,itype) = (pi/6.0)*(dlo_sect( isize,itype)**3)
          volumhi_sect( isize,itype) = (pi/6.0)*(dhi_sect( isize,itype)**3)
       end do
       write(*,'(a,i3,1p,5e11.3)') 'type, sg, dg, dp, vol', itype, &
            sigmag_aer(1,itype), dgnum_amode(itype), dp_meanvol_tmp, volumcen_sect(1,itype)
       
       
       do l = 1, ncomp_aer(itype)
          l2 = lspectype_amode(l,itype)
          if ((l2 > 0) .and. (l2 <= ntot_mastercomp_aer)) then
             mastercompptr_aer(l,itype) = l2
          else
             msg = '*** cam_mam_init_asect error 200 - mastercompptr'
             call wrf_message( msg )
             write( msg, '(a,4(1x,i10))' ) &
                  'itype, l, l2, ntot_mastercomp_aer = ', &
                  itype, l, l2, ntot_mastercomp_aer
             call wrf_error_fatal3("<stdin>",703,&
msg )
          end if
          dens_aer(l,itype) = dens_mastercomp_aer(l2)
          mw_aer(l,itype) = mw_mastercomp_aer(l2)
          hygro_aer(l,itype) = hygro_mastercomp_aer(l2)
          name_aer(l,itype) = name_mastercomp_aer(l2)
       end do
       
       isize = 1
       do l = -1, ncomp_aer(itype)
          do iphase = 1, nphase_aer
             if (l == -1) then
                tmpname = 'num'
             else if (l == 0) then
                tmpname = 'wtr'
                if (iphase > 1) cycle
             else
                l2 = lspectype_amode(l,itype)
                tmpname = namebb_mastercomp_aer(l2)
             end if
             
             if (iphase == 1) then
                tmpname = trim(tmpname) // '_a'
             else
                tmpname = trim(tmpname) // '_c'
             end if
             write( tmptxtaa, '(i1)' ) itype
             tmpname = trim(tmpname) // tmptxtaa(1:1)
             
             l3 = 0
             do l4 = p1st, num_chem
                if (chem_dname_table(1,l4) == tmpname) then
                   l3 = l4
                   exit
                end if
             end do
             if (l3 <= 0) then
                msg = '*** cam_mam_init_asect error 300' // &
                     ' - finding species - ' // tmpname
                call wrf_error_fatal3("<stdin>",743,&
msg )
             end if
             
             if (l == -1) then
                numptr_aer(isize,itype,iphase) = l3
                mprognum_aer(isize,itype,iphase) = 1
             else if (l == 0) then
                waterptr_aer(isize,itype) = l3
             else
                massptr_aer(l,isize,itype,iphase) = l3
                mastercompptr_aer(l,itype) = l2
                if (l2 == mastercompindx_so4_aer) then
                   lptr_so4_aer(isize,itype,iphase) = l3
                else if (l2 == mastercompindx_nh4_aer) then
                   lptr_nh4_aer(isize,itype,iphase) = l3
                   
                   
                else if (l2 == mastercompindx_pom_aer) then
                   lptr_pom_aer(isize,itype,iphase) = l3
                else if (l2 == mastercompindx_soa_aer) then
                   lptr_soa_aer(isize,itype,iphase) = l3
                else if (l2 == mastercompindx_bc_aer) then
                   lptr_bc_aer(isize,itype,iphase) = l3
                else if (l2 == mastercompindx_dust_aer) then
                   lptr_dust_aer(isize,itype,iphase) = l3
                else if (l2 == mastercompindx_seas_aer) then
                   lptr_seas_aer(isize,itype,iphase) = l3
                else
                   msg = '*** cam_mam_init_asect error 400' // &
                        ' - finding species type - ' // tmpname
                   call wrf_error_fatal3("<stdin>",774,&
msg )
                end if
             end if
             
          end do 
          
       end do 
       
    end do 
    
    
    
    write(*,'(/a)') 'cam_mam_init_asect diagnostics'
    
    write(*,'(/a,i5)') 'ntot_mastercomp_aer', ntot_mastercomp_aer
    write(*,'(a)') 'mastercomp name, l, mw, dens, hygro'
    do l = 1, ntot_mastercomp_aer
       write(*,'(a,i5,1p,3e12.4)') name_mastercomp_aer(l), l, &
            mw_mastercomp_aer(l), dens_mastercomp_aer(l), hygro_mastercomp_aer(l) 
    end do
    write(*,'(a)') &
         'mastercompindx_so4_aer, nh4, no3, pom, soa, bc, seas, dust'
    write(*,'(4i12)') &
         mastercompindx_so4_aer, mastercompindx_nh4_aer, mastercompindx_no3_aer, &
         mastercompindx_pom_aer, mastercompindx_soa_aer, mastercompindx_bc_aer, &
         mastercompindx_seas_aer, mastercompindx_dust_aer
    write(*,'(a)') '........... mw_so4_aer, nh4, no3, pom, soa, bc, seas, dust'
    write(*,'(1p,4e12.4)') &
         mw_so4_aer, mw_nh4_aer, mw_no3_aer, &
         mw_pom_aer, mw_soa_aer, mw_bc_aer, &
         mw_seas_aer, mw_dust_aer
    write(*,'(a)') '......... dens_so4_aer, nh4, no3, pom, soa, bc, seas, dust'
    write(*,'(1p,4e12.4)') &
         dens_so4_aer, dens_nh4_aer, dens_no3_aer, &
         dens_pom_aer, dens_soa_aer, dens_bc_aer, &
         dens_seas_aer, dens_dust_aer
    
    write(*,'(/a/6i12)') 'nphase_aer, ai_phase, cw_phase', &
         nphase_aer, ai_phase, cw_phase
    
    do itype = 1, ntype_aer
       do isize = 1, nsize_aer(itype)
          write(*,'(/a,2i5,a)') 'info for itype, isize = ', itype, isize, &
               'species;  id, name for ai & cw;  mw, dens, hygro'
          
          la = numptr_aer(isize,itype,1)
          lc = numptr_aer(isize,itype,2)
          tmpname = '---'
          if ((la >= p1st) .and. (la <= num_chem)) tmpname = chem_dname_table(1,la)
          tmpnamec = '---'
          if ((lc >= p1st) .and. (lc <= num_chem)) tmpnamec = chem_dname_table(1,la)
          write(*,'(a,i12,1x,a,i12,1x,a,1p,3e12.4)') 'number    ', &
               la, tmpname(1:10), lc, tmpnamec(1:10)
          
          do l = 1, ncomp_aer(itype)
             la = massptr_aer(l,isize,itype,1)
             lc = massptr_aer(l,isize,itype,2)
             tmpname = '---'
             if ((la >= p1st) .and. (la <= num_chem)) tmpname = chem_dname_table(1,la)
             tmpnamec = '---'
             if ((lc >= p1st) .and. (lc <= num_chem)) tmpnamec = chem_dname_table(1,la)
             write(*,'(a,i12,1x,a,i12,1x,a,1p,3e12.4)') name_aer(l,itype), &
                  la, tmpname(1:10), lc, tmpnamec(1:10), &
                  mw_aer(l,itype), dens_aer(l,itype), hygro_aer(l,itype) 
          end do 
          
          la = waterptr_aer(isize,itype)
          tmpname = '---'
          if ((la >= p1st) .and. (la <= num_chem)) tmpname = chem_dname_table(1,la)
          write(*,'(a,i12,1x,a,23x,1p,3e12.4)') 'water     ', &
               la, tmpname(1:10)
          
          la = hyswptr_aer(isize,itype)
          tmpname = '---'
          if ((la >= p1st) .and. (la <= num_chem)) tmpname = chem_dname_table(1,la)
          write(*,'(a,i12,1x,a,23x,1p,3e12.4)') 'hys-water ', &
               la, tmpname(1:10)
          
       end do 
    end do 
    
    
    return
  end subroutine cam_mam_init_asect
  
  
  
  
  subroutine cam_mam_init_other( id, numgas, config_flags )
    
    
    
    
    
    
    use shr_kind_mod,             only: r8 => shr_kind_r8
    use module_state_description, only: num_chem, param_first_scalar
    use module_configure,         only: grid_config_rec_type, &
         p_h2o2, p_hno3, p_nh3, p_o3, p_so2, p_soag, p_sulf
    use module_scalar_tables,     only: chem_dname_table
    
    use module_cam_support,       only: pcnst => pcnst_runtime, &
         pcnst_non_chem => pcnst_non_chem_modal_aero, &
         gas_pcnst => gas_pcnst_modal_aero
    use modal_aero_data,          only: cnst_name_cw
    use module_data_cam_mam_asect
    use constituents,             only: cnst_mw, cnst_rgas, cnst_cv, cnst_cp, &
         qmin,qmincg
    use physconst,                only: r_universal
    use infnan,                   only: nan


    

    
    
    
    implicit none
    
    
    integer, intent(in) :: id, numgas
    type(grid_config_rec_type), intent(in) :: config_flags
    
    
    integer :: iphase, isize, itype, dumind
    integer :: l, ll, l2, l3, l4
    integer :: p1st
    
    character(len=160) :: msg
    character(len=16)  :: tmpname, tmpname2, tmpname3
    real(r8)           :: qmin_gas, qmin_aer,qmin_num

    IF(.NOT.CAM_INITIALIZED_CHEM) THEN
       qmin_gas = 1.0E-17_r8 
       qmin_aer = 1.0E-14_r8 
       qmin_num = 1.0E+01_r8 
    ENDIF
    
    p1st = param_first_scalar
    
    
    if ( .not. allocated(lptr_chem_to_q) ) then
       allocate( lptr_chem_to_q(num_chem) )
    end if
    
    if ( .not. allocated(lptr_chem_to_qqcw) ) then
       allocate( lptr_chem_to_qqcw(num_chem) )
    end if
    
    if ( .not. allocated(factconv_chem_to_q) ) then
       allocate( factconv_chem_to_q(num_chem) )
    end if
    
    if ( .not. allocated(factconv_chem_to_qqcw) ) then
       allocate( factconv_chem_to_qqcw(num_chem) )
    end if
    
    if ( .not. allocated(mw_chem_array) ) then
       allocate( mw_chem_array(num_chem) )
    end if
    
    if ( .not. allocated(mw_q_array) ) then
       allocate( mw_q_array(pcnst) )
    end if
    
    if ( .not. allocated(mw_q_mo_array) ) then
       allocate( mw_q_mo_array(gas_pcnst) )
    end if
    
    
    
      lptr_chem_to_q(:) = init_val
      lptr_chem_to_qqcw(:) = init_val

      do l = p1st, num_chem
         tmpname = chem_dname_table(1,l)
         do l2 = 1, pcnst
            if (tmpname == cnst_name_loc(l2)) then
               lptr_chem_to_q(l) = l2
               exit
            end if
            if (tmpname == cnst_name_cw(l2)) then
               lptr_chem_to_qqcw(l) = l2
               exit
            end if
         end do
      end do



      factconv_chem_to_q(:) = 1.0
      factconv_chem_to_qqcw(:) = 1.0
      mw_chem_array(:) = 1.0
      mw_q_array(:) = 1.0
      mw_q_mo_array(:) = 1.0

      l2 = pcnst_non_chem
      do l = p1st, numgas
         l2 = lptr_chem_to_q(l)
         if ((l2 < 1) .or. (l2 > pcnst)) cycle


         if (l == p_sulf) mw_chem_array(l) = 96.0
         if (l == p_so2 ) mw_chem_array(l) = 64.0
         if (l == p_nh3 ) mw_chem_array(l) = 17.0
         if (l == p_hno3) mw_chem_array(l) = 63.0
         if (l == p_soag) mw_chem_array(l) = 12.0
         if (l == p_h2o2) mw_chem_array(l) = 34.0
         if (l == p_o3  ) mw_chem_array(l) = 48.0
         mw_q_array(l2) = mw_chem_array(l)
         IF(.NOT.CAM_INITIALIZED_CHEM) THEN
            cnst_mw(l2)    = mw_q_array(l2)
            cnst_rgas(l2) = r_universal * cnst_mw(l2)
            cnst_cp  (l2) = nan
            cnst_cv  (l2) = nan
            qmin     (l2) = qmin_gas
            qmincg   (l2) = qmin(l2)
         ENDIF
         
         factconv_chem_to_q(l) = 1.0e-6*mw_chem_array(l)/28.966
      end do

      iphase = ai_phase
      do itype = 1, ntype_aer
      do isize = 1, nsize_aer(itype)
         l = numptr_aer(isize,itype,iphase)
         mw_chem_array(l) = 1.0
         l2 = lptr_chem_to_q(l)
         if ((l2 >= 1) .and. (l2 <= pcnst)) then 
            mw_q_array(l2) = mw_chem_array(l)
            IF(.NOT.CAM_INITIALIZED_CHEM) THEN
               cnst_mw(l2)    = mw_q_array(l2)
               cnst_rgas(l2) = r_universal * cnst_mw(l2)
               cnst_cp  (l2) = nan
               cnst_cv  (l2) = nan
               qmin     (l2) = qmin_num
               qmincg   (l2) = qmin(l2)
            ENDIF
            
            factconv_chem_to_q(l) = 1.0
         end if

         l = waterptr_aer(isize,itype)
         if ((l >= p1st) .and. (l <= num_chem)) then 
            mw_chem_array(l) = 18.0
            
            factconv_chem_to_q(l) = 1.0e-9
         end if

         do ll = 1, ncomp_aer(itype)
            l = massptr_aer(ll,isize,itype,iphase)
            mw_chem_array(l) = mw_aer(ll,itype)
            l2 = lptr_chem_to_q(l)
            if ((l2 >= 1) .and. (l2 <= pcnst)) then 
               mw_q_array(l2) = mw_chem_array(l)
               IF(.NOT.CAM_INITIALIZED_CHEM) THEN
                  cnst_mw(l2)    = mw_q_array(l2)
                  cnst_rgas(l2) = r_universal * cnst_mw(l2)
                  cnst_cp  (l2) = nan
                  cnst_cv  (l2) = nan
                  qmin     (l2) = qmin_aer
                  qmincg   (l2) = qmin(l2)
               ENDIF
               
               factconv_chem_to_q(l) = 1.0e-9
            end if
         end do
      end do 
      end do 



      if (nphase_aer > 1) then
      iphase = cw_phase
      do itype = 1, ntype_aer
      do isize = 1, nsize_aer(itype)
         l = numptr_aer(isize,itype,iphase)
         if (l < p1st .or. l > num_chem) then
            write(*,'(a,10i10)') '*** cw_phase numb error', iphase, itype, isize, l
         else
         mw_chem_array(l) = 1.0
         l2 = lptr_chem_to_qqcw(l)
         if ((l2 >= 1) .and. (l2 <= pcnst)) then

            
            factconv_chem_to_q(l) = 1.0
         end if
         end if
         
         do ll = 1, ncomp_aer(itype)
            l = massptr_aer(ll,isize,itype,iphase)
            if (l < p1st .or. l > num_chem) then
               write(*,'(a,10i10)') '*** cw_phase mass error', iphase, itype, isize, ll, l
            else
            mw_chem_array(l) = mw_aer(ll,itype)
            l2 = lptr_chem_to_qqcw(l)
            if ((l2 >= 1) .and. (l2 <= pcnst)) then

               
               factconv_chem_to_q(l) = 1.0e-9
            end if
            end if
         end do
      end do 
      end do 
      end if



      mw_q_mo_array(1:gas_pcnst) = mw_q_array(pcnst_non_chem+1:pcnst)


      write( *, '(/2a)' ) &
         'l, cnst_name, chem_name, l3, chem_name2, ', &
         'lptr_chem_to_q, factconv_..., mw_...'
      do l = 1, max( pcnst, num_chem )
         tmpname  = ' '
         tmpname2 = ' '
         tmpname3 = ' '
         if (l <= pcnst) then
            tmpname = cnst_name_loc(l)
            do l2 = p1st, num_chem
               if (lptr_chem_to_q(l2) == l) tmpname2 = chem_dname_table(1,l2)
            end do
         end if
         l3 = l
         l4 = init_val
         if ((l3 >= p1st) .and. (l3 <= num_chem)) then
            tmpname3 = chem_dname_table(1,l3)
            l4 = lptr_chem_to_q(l3)
         end if
         if (l3 <= num_chem) then
            write( *, '(i4,2(2x,a),i6,2x,a,i12,1p,2e10.2)' ) &
               l, tmpname, tmpname2, l3, tmpname3, l4, &
               factconv_chem_to_q(l3), mw_chem_array(l3)
         else
            write( *, '(i4,2(2x,a),i6,2x,a,i12,1p,e10.2)' ) &
               l, tmpname, tmpname2
         end if
      end do

      IF(.NOT.CAM_INITIALIZED_CHEM) CAM_INITIALIZED_CHEM = .TRUE.
      return
      end subroutine cam_mam_init_other




      end module module_cam_mam_init

