











module modal_aero_initialize_data
  use module_cam_support,    only: iulog, endrun, masterproc, iam, &
       pcols, pver
  use modal_aero_data

  implicit none
  private
  public :: modal_aero_register
  public :: modal_aero_initialize
  public :: decouple_mam_mp

contains

  subroutine modal_aero_register
    use module_cam_support,    only: pcnst => pcnst_runtime
    use constituents,          only: cnst_name

    character(len=8)  :: &
         xname_massptr(maxd_aspectype,ntot_amode), &
         xname_massptrcw(maxd_aspectype,ntot_amode)
    character(len=10) :: xname_spectype(maxd_aspectype,ntot_amode)


    
    character(len=*), parameter ::     xname_numptr(ntot_amode)   = (/ 'num_a1  ', 'num_a2  ', &
         'num_a3  ' /)
    character(len=*), parameter ::     xname_numptrcw(ntot_amode) = (/ 'num_c1  ', 'num_c2  ', &
         'num_c3  ' /)



    integer :: m, l, iptr
    real pi
    character(len=3) :: trnum       

    pi = 4.*atan(1._r8)    

       


       
       
       
       
       xname_massptr(:nspec_amode(1),1)   = (/ 'so4_a1  ', &
            'pom_a1  ', 'soa_a1  ', 'bc_a1   ', &
            'dst_a1  ', 'ncl_a1  ' /)
       xname_massptrcw(:nspec_amode(1),1) = (/ 'so4_c1  ', &
            'pom_c1  ', 'soa_c1  ', 'bc_c1   ', &
            'dst_c1  ', 'ncl_c1  ' /)
       xname_spectype(:nspec_amode(1),1)  = (/ 'sulfate   ', &
            'p-organic ', 's-organic ', 'black-c   ', &
            'dust      ', 'seasalt   ' /)

       
       xname_massptr(:nspec_amode(2),2)   = (/ 'so4_a2  ', &
            'soa_a2  ', 'ncl_a2  ' /)
       xname_massptrcw(:nspec_amode(2),2) = (/ 'so4_c2  ', &
            'soa_c2  ', 'ncl_c2  ' /)
       xname_spectype(:nspec_amode(2),2)  = (/ 'sulfate   ', &
            's-organic ', 'seasalt   ' /)

       
       xname_massptr(:nspec_amode(3),3)   = (/ 'dst_a3  ', 'ncl_a3  ', 'so4_a3  ' /)
       xname_massptrcw(:nspec_amode(3),3) = (/ 'dst_c3  ', 'ncl_c3  ', 'so4_c3  ' /)
       xname_spectype(:nspec_amode(3),3)  = (/ 'dust      ', 'seasalt   ', 'sulfate   ' /)



    do m = 1, ntot_amode

       if (masterproc) then
          write(iulog,9231) m, modename_amode(m)
          write(iulog,9232)                                          &
               'nspec                       ',                         &
               nspec_amode(m)
          write(iulog,9232)                                          &
               'mprognum, mdiagnum, mprogsfc',                         &
               mprognum_amode(m), mdiagnum_amode(m), mprogsfc_amode(m)
          write(iulog,9232)                                          &
               'mcalcwater                  ',                         &
               mcalcwater_amode(m)
       endif

       
       
       alnsg_amode(m) = log( sigmag_amode(m) )

       voltonumb_amode(m) = 1. / ( (pi/6.)*                            &
            (dgnum_amode(m)**3.)*exp(4.5*alnsg_amode(m)**2.) )
       voltonumblo_amode(m) = 1. / ( (pi/6.)*                          &
            (dgnumlo_amode(m)**3.)*exp(4.5*alnsg_amode(m)**2.) )
       voltonumbhi_amode(m) = 1. / ( (pi/6.)*                          &
            (dgnumhi_amode(m)**3.)*exp(4.5*alnsg_amode(m)**2.) )

       alnv2n_amode(m)   = log( voltonumb_amode(m) )
       alnv2nlo_amode(m) = log( voltonumblo_amode(m) )
       alnv2nhi_amode(m) = log( voltonumbhi_amode(m) )

       
       call search_list_of_names(                                      &
            xname_numptr(m), numptr_amode(m), cnst_name, pcnst )
       if (numptr_amode(m) .le. 0) then
          write(iulog,9061) 'xname_numptr', xname_numptr(m), m
          call endrun()
       end if
       if (numptr_amode(m) .gt. pcnst) then
          write(iulog,9061) 'numptr_amode', numptr_amode(m), m
          write(iulog,9061) 'xname_numptr', xname_numptr(m), m
          call endrun()
       end if

       species_class(numptr_amode(m)) = spec_class_aerosol


       numptrcw_amode(m) = numptr_amode(m)  
       if (numptrcw_amode(m) .le. 0) then
          write(iulog,9061) 'xname_numptrcw', xname_numptrcw(m), m
          call endrun()
       end if
       if (numptrcw_amode(m) .gt. pcnst) then
          write(iulog,9061) 'numptrcw_amode', numptrcw_amode(m), m
          write(iulog,9061) 'xname_numptrcw', xname_numptrcw(m), m
          call endrun()
       end if
       species_class(numptrcw_amode(m)) = spec_class_aerosol

       
       if ( masterproc ) then
          write(iulog,9233) 'numptr         ',                           &
               numptr_amode(m), xname_numptr(m)
          write(iulog,9233) 'numptrcw       ',                           &
               numptrcw_amode(m), xname_numptrcw(m)
       end if


       
       do l = 1, nspec_amode(m)

          call search_list_of_names(                                  &
               xname_spectype(l,m), lspectype_amode(l,m),              &
               specname_amode, ntot_aspectype )
          if (lspectype_amode(l,m) .le. 0) then
             write(iulog,9062) 'xname_spectype', xname_spectype(l,m), l, m
             call endrun()
          end if

          call search_list_of_names(                                  &
               xname_massptr(l,m), lmassptr_amode(l,m), cnst_name, pcnst )
          if (lmassptr_amode(l,m) .le. 0) then
             write(iulog,9062) 'xname_massptr', xname_massptr(l,m), l, m
             call endrun()
          end if
          species_class(lmassptr_amode(l,m)) = spec_class_aerosol

          lmassptrcw_amode(l,m) = lmassptr_amode(l,m)  
          if (lmassptrcw_amode(l,m) .le. 0) then
             write(iulog,9062) 'xname_massptrcw', xname_massptrcw(l,m), l, m
             call endrun()
          end if
          species_class(lmassptrcw_amode(l,m)) = spec_class_aerosol

          if ( masterproc ) then
             write(iulog,9236) 'spec, spectype ', l,                    &
                  lspectype_amode(l,m), xname_spectype(l,m)
             write(iulog,9236) 'spec, massptr  ', l,                    &
                  lmassptr_amode(l,m), xname_massptr(l,m)
             write(iulog,9236) 'spec, massptrcw', l,                    &
                  lmassptrcw_amode(l,m), xname_massptrcw(l,m)
          end if

       enddo

       if ( masterproc ) write(iulog,*)


       
       write(unit=trnum,fmt='(i3)') m+100
       aodvisname(m) = 'AODVIS'//trnum(2:3)
       aodvislongname(m) = 'Aerosol optical depth for mode '//trnum(2:3)
       ssavisname(m) = 'SSAVIS'//trnum(2:3)
       ssavislongname(m) = 'Single-scatter albedo for mode '//trnum(2:3)
       fnactname(m) = 'FNACT'//trnum(2:3)
       fnactlongname(m) = 'Number faction activated for mode '//trnum(2:3)
       fmactname(m) = 'FMACT'//trnum(2:3)
       fmactlongname(m) = 'Fraction mass activated for mode'//trnum(2:3)
    end do
9230   format( // '*** init_aer_modes mode definitions' )
9231   format( 'mode = ', i4, ' = "', a, '"' )
9232   format( 4x, a, 4(1x, i5 ) )
9233   format( 4x, a15, 4x, i7, '="', a, '"' )
9236   format( 4x, a15, i4, i7, '="', a, '"' )
9061   format( '*** subr init_aer_modes - bad ', a /                   &
            5x, 'name, m =  ', a, 5x, i5 )
9062   format( '*** subr init_aer_modesaeromodeinit - bad ', a /                       &
            5x, 'name, l, m =  ', a, 5x, 2i5 )


  end subroutine modal_aero_register


  
  subroutine modal_aero_initialize
       use module_cam_support,    only: pcnst => pcnst_runtime, &
            addfld, add_default, phys_decomp
       use physconst,             only: rhoh2o, mwh2o
       use modal_aero_calcsize,   only: modal_aero_calcsize_init
       use modal_aero_coag,       only: modal_aero_coag_init
       use modal_aero_gasaerexch, only: modal_aero_gasaerexch_init
       use modal_aero_newnuc,     only: modal_aero_newnuc_init
       use modal_aero_rename,     only: modal_aero_rename_init
       use mz_aerosols_intr,      only: modal_aero_bcscavcoef_init

       implicit none

       
       
       
       integer l, m, i


       character(len=3) :: trnum       
       integer :: iaerosol, ibulk
       integer  :: numaerosols     
       character(len=20) :: bulkname
       complex, pointer  :: refindex_aer_sw(:), &
            refindex_aer_lw(:)
       real(r8) :: hygro_aer
       logical  :: history_aerosol      

       
       history_aerosol = .false.



       
       
       
       
       

       
       
       
       
       
       
       
       
       
       

       
       
       
       
       spechygro(:ntot_aspectype) = (/ 0.507, 0.507, 0.507, &
            0.100, 0.140, 1.0e-10, &
            1.160, 0.068 /)
       specrefndxsw(1,:ntot_aspectype)     = (/ (1.53,  0.01),   (1.53,  0.01),  (1.53,  0.01), &
            (1.55,  0.01),   (1.55,  0.01),  (1.90, 0.60), &
            (1.50, 1.0e-8), (1.50, 0.005) /)
       specrefndxlw(1,:ntot_aspectype)   = (/ (2.0, 0.5),   (2.0, 0.5), (2.0, 0.5), &
            (1.7, 0.5),   (1.7, 0.5), (2.22, 0.73), &
            (1.50, 0.02), (2.6, 0.6) /)
       do l = 1, ntot_aspectype
          specrefndxsw(:,l)=specrefndxsw(1,l)
          specrefndxlw(:,l)=specrefndxlw(1,l)
       end do
       
       
       specdens_amode(:ntot_aspectype) = (/1770.0,1770.0,1770.0, 1000.0, 1000.0, 1700.0,1900.0,2600.0 /)
       
9210   format( // '*** init_aer_modes aerosol species-types' )
9211   format( 'spectype =', i4)
9212   format( 4x, a, 3x, '"', a, '"' )
9213   format( 4x, a, 5(1pe14.5) )





       do i = 1, pcnst
          species_class(i) = spec_class_undefined
       end do



       
       call initaermodes_set_cnstnamecw()


       
       
       
       call initaermodes_setspecptrs

       if ( masterproc ) write(iulog,*)



       
       
       
       do m = 1, ntot_amode
          write( trnum, '(i3.3)' ) m
          
          call addfld( &
               'dgnd_a'//trnum(2:3), 'm', pver, 'A', &
               'dry dgnum, interstitial, mode '//trnum(2:3), phys_decomp )
          call addfld( &
               'dgnw_a'//trnum(2:3), 'm', pver, 'A', &
               'wet dgnum, interstitial, mode '//trnum(2:3), phys_decomp )
          call addfld( &
               'wat_a'//trnum(3:3), 'm', pver, 'A', &
               'aerosol water, interstitial, mode '//trnum(2:3), phys_decomp )
          if ( history_aerosol ) then    
            call add_default( 'dgnd_a'//trnum(2:3), 1, ' ' )
            call add_default( 'dgnw_a'//trnum(2:3), 1, ' ' )
            call add_default( 'wat_a'//trnum(3:3),  1, ' ' )     
          endif

          l = lptr_so4_cw_amode(m)
          if (l > 0) then
             call addfld (&
                  trim(cnst_name_cw(l))//'AQSO4','kg/m2/s ',1,  'A', &
                  trim(cnst_name_cw(l))//' aqueous phase chemistry',phys_decomp)
             call addfld (&
                  trim(cnst_name_cw(l))//'AQH2SO4','kg/m2/s ',1,  'A', &
                  trim(cnst_name_cw(l))//' aqueous phase chemistry',phys_decomp)
             if ( history_aerosol ) then 
                call add_default (trim(cnst_name_cw(l))//'AQSO4', 1, ' ')
                call add_default (trim(cnst_name_cw(l))//'AQH2SO4', 1, ' ')
             endif
          end if

       end do

       call addfld ('AQSO4_H2O2','kg/m2/s ',1,  'A', &
            'SO4 aqueous phase chemistry due to H2O2',phys_decomp)
       call addfld ('AQSO4_O3','kg/m2/s ',1,  'A', &
            'SO4 aqueous phase chemistry due to O3',phys_decomp)
       call addfld( 'XPH_LWC','kg/kg   ',pver, 'A', &
            'pH value multiplied by lwc', phys_decomp)

       if ( history_aerosol ) then    
          call add_default ('AQSO4_H2O2', 1, ' ')
          call add_default ('AQSO4_O3', 1, ' ')    
          call add_default ('XPH_LWC', 1, ' ')
       endif



       
       
       
       
       
       
       
       
       
       qneg3_worst_thresh_amode(:) = 0.0_r8
       do m = 1, ntot_amode
          l = numptr_amode(m)
          if ((l <= 0) .or. (l > pcnst)) cycle

          if      (m == modeptr_accum) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_aitken) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_pcarbon) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_ufine) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8

          else if (m == modeptr_fineseas) then
             qneg3_worst_thresh_amode(l) = 3.0e1_r8
          else if (m == modeptr_finedust) then
             qneg3_worst_thresh_amode(l) = 3.0e1_r8

          else
             qneg3_worst_thresh_amode(l) = 1.0e0_r8
          end if

          if ( masterproc ) write(iulog,'(i3,2x,a,1p,e12.3)') &
               m, modename_amode(m), qneg3_worst_thresh_amode(l)
       end do


       
       
       
       call modal_aero_rename_init
       
       call modal_aero_calcsize_init
       call modal_aero_gasaerexch_init
       
       call modal_aero_coag_init
       call modal_aero_newnuc_init
       call modal_aero_bcscavcoef_init

       return
     end subroutine modal_aero_initialize


     
     subroutine search_list_of_names(                                &
          name_to_find, name_id, list_of_names, list_length )
       
       
       
       
       
       
       
       
       
       character(len=*), intent(in):: name_to_find, list_of_names(:)
       integer, intent(in) :: list_length
       integer, intent(out) :: name_id
       
       integer :: i
       name_id = -999888777
       if (name_to_find .ne. ' ') then
          do i = 1, list_length
             if (name_to_find .eq. list_of_names(i)) then
                name_id = i
                exit
             end if
          end do
       end if
     end subroutine search_list_of_names


     
     subroutine initaermodes_setspecptrs
       
       
       
       
       
       
       
       
       implicit none

       
       integer l, l2, m
       character*8 dumname
       integer, parameter :: init_val=-999888777

       

       modeptr_accum = init_val
       modeptr_aitken = init_val
       modeptr_ufine = init_val
       modeptr_coarse = init_val
       modeptr_pcarbon = init_val
       modeptr_fineseas = init_val
       modeptr_finedust = init_val
       modeptr_coarseas = init_val
       modeptr_coardust = init_val
       do m = 1, ntot_amode
          if (modename_amode(m) .eq. 'accum') then
             modeptr_accum = m
          else if (modename_amode(m) .eq. 'aitken') then
             modeptr_aitken = m
          else if (modename_amode(m) .eq. 'ufine') then
             modeptr_ufine = m
          else if (modename_amode(m) .eq. 'coarse') then
             modeptr_coarse = m
          else if (modename_amode(m) .eq. 'primary carbon') then
             modeptr_pcarbon = m
          else if (modename_amode(m) .eq. 'fine seasalt') then
             modeptr_fineseas = m
          else if (modename_amode(m) .eq. 'fine dust') then
             modeptr_finedust = m
          else if (modename_amode(m) .eq. 'coarse seasalt') then
             modeptr_coarseas = m
          else if (modename_amode(m) .eq. 'coarse dust') then
             modeptr_coardust = m
          end if
       end do

       do m = 1, ntot_amode
          lptr_so4_a_amode(m)   = init_val
          lptr_so4_cw_amode(m)  = init_val
          lptr_msa_a_amode(m)   = init_val
          lptr_msa_cw_amode(m)  = init_val
          lptr_nh4_a_amode(m)   = init_val
          lptr_nh4_cw_amode(m)  = init_val
          lptr_no3_a_amode(m)   = init_val
          lptr_no3_cw_amode(m)  = init_val
          lptr_pom_a_amode(m)   = init_val
          lptr_pom_cw_amode(m)  = init_val
          lptr_soa_a_amode(m)   = init_val
          lptr_soa_cw_amode(m)  = init_val
          lptr_bc_a_amode(m)    = init_val
          lptr_bc_cw_amode(m)   = init_val
          lptr_nacl_a_amode(m)  = init_val
          lptr_nacl_cw_amode(m) = init_val
          lptr_dust_a_amode(m)  = init_val
          lptr_dust_cw_amode(m) = init_val
          do l = 1, nspec_amode(m)
             l2 = lspectype_amode(l,m)
             if ( (specname_amode(l2) .eq. 'sulfate') .and.  &
                  (lptr_so4_a_amode(m) .le. 0) ) then
                lptr_so4_a_amode(m)  = lmassptr_amode(l,m)
                lptr_so4_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'msa') .and.      &
                  (lptr_msa_a_amode(m) .le. 0) ) then
                lptr_msa_a_amode(m)  = lmassptr_amode(l,m)
                lptr_msa_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'ammonium') .and.  &
                  (lptr_nh4_a_amode(m) .le. 0) ) then
                lptr_nh4_a_amode(m)  = lmassptr_amode(l,m)
                lptr_nh4_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'nitrate') .and.  &
                  (lptr_no3_a_amode(m) .le. 0) ) then
                lptr_no3_a_amode(m)  = lmassptr_amode(l,m)
                lptr_no3_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'p-organic') .and.   &
                  (lptr_pom_a_amode(m) .le. 0) ) then
                lptr_pom_a_amode(m)  = lmassptr_amode(l,m)
                lptr_pom_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 's-organic') .and.   &
                  (lptr_soa_a_amode(m) .le. 0) ) then
                lptr_soa_a_amode(m)  = lmassptr_amode(l,m)
                lptr_soa_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'black-c') .and.  &
                  (lptr_bc_a_amode(m) .le. 0) ) then
                lptr_bc_a_amode(m)  = lmassptr_amode(l,m)
                lptr_bc_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'seasalt') .and.  &
                  (lptr_nacl_a_amode(m) .le. 0) ) then
                lptr_nacl_a_amode(m)  = lmassptr_amode(l,m)
                lptr_nacl_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'dust') .and.     &
                  (lptr_dust_a_amode(m) .le. 0) ) then
                lptr_dust_a_amode(m)  = lmassptr_amode(l,m)
                lptr_dust_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
          end do
       end do

       
       specdens_so4_amode = 2.0
       specdens_nh4_amode = 2.0
       specdens_no3_amode = 2.0
       specdens_pom_amode = 2.0
       specdens_soa_amode = 2.0
       specdens_bc_amode = 2.0
       specdens_dust_amode = 2.0
       specdens_seasalt_amode = 2.0
       specmw_so4_amode = 1.0
       specmw_nh4_amode = 1.0
       specmw_no3_amode = 1.0
       specmw_pom_amode = 1.0
       specmw_soa_amode = 1.0
       specmw_bc_amode = 1.0
       specmw_dust_amode = 1.0
       specmw_seasalt_amode = 1.0
       do m = 1, ntot_aspectype
          if      (specname_amode(m).eq.'sulfate   ') then
             specdens_so4_amode = specdens_amode(m)
             specmw_so4_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'ammonium  ') then
             specdens_nh4_amode = specdens_amode(m)
             specmw_nh4_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'nitrate   ') then
             specdens_no3_amode = specdens_amode(m)
             specmw_no3_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'p-organic ') then
             specdens_pom_amode = specdens_amode(m)
             specmw_pom_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'s-organic ') then
             specdens_soa_amode = specdens_amode(m)
             specmw_soa_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'black-c   ') then
             specdens_bc_amode = specdens_amode(m)
             specmw_bc_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'dust      ') then
             specdens_dust_amode = specdens_amode(m)
             specmw_dust_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'seasalt   ') then
             specdens_seasalt_amode = specdens_amode(m)
             specmw_seasalt_amode = specmw_amode(m)
          end if
       enddo

       
       if ( .not. ( masterproc ) ) return
       write(iulog,*) 'modeptr_accum    =', modeptr_accum
       write(iulog,*) 'modeptr_aitken   =', modeptr_aitken
       write(iulog,*) 'modeptr_ufine    =', modeptr_ufine
       write(iulog,*) 'modeptr_coarse   =', modeptr_coarse
       write(iulog,*) 'modeptr_pcarbon  =', modeptr_pcarbon
       write(iulog,*) 'modeptr_fineseas =', modeptr_fineseas
       write(iulog,*) 'modeptr_finedust =', modeptr_finedust
       write(iulog,*) 'modeptr_coarseas =', modeptr_coarseas
       write(iulog,*) 'modeptr_coardust =', modeptr_coardust

       dumname = 'none'
       write(iulog,9000) 'sulfate    '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_so4_a_amode(m), lptr_so4_cw_amode(m),  'so4' )
       end do

       write(iulog,9000) 'msa        '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_msa_a_amode(m), lptr_msa_cw_amode(m),  'msa' )
       end do

       write(iulog,9000) 'ammonium   '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_nh4_a_amode(m), lptr_nh4_cw_amode(m),  'nh4' )
       end do

       write(iulog,9000) 'nitrate    '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_no3_a_amode(m), lptr_no3_cw_amode(m),  'no3' )
       end do

       write(iulog,9000) 'p-organic  '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_pom_a_amode(m), lptr_pom_cw_amode(m),  'pom' )
       end do

       write(iulog,9000) 's-organic  '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_soa_a_amode(m), lptr_soa_cw_amode(m),  'soa' )
       end do

       write(iulog,9000) 'black-c    '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_bc_a_amode(m), lptr_bc_cw_amode(m),  'bc' )
       end do

       write(iulog,9000) 'seasalt   '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_nacl_a_amode(m), lptr_nacl_cw_amode(m),  'nacl' )
       end do

       write(iulog,9000) 'dust       '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_dust_a_amode(m), lptr_dust_cw_amode(m),  'dust' )
       end do

9000   format( a )
9230   format(                                                         &
            / 'mode-pointer output from subr initaermodes_setspecptrs' )
9240   format(                                                         &
            / 'species-pointer output from subr initaermodes_setspecptrs' / &
            'mode', 12x, 'id  name_a  ', 12x, 'id  name_cw' )

       return
     end subroutine initaermodes_setspecptrs


     
     subroutine initaermodes_setspecptrs_write2(                     &
          m, laptr, lcptr, txtdum )
       
       
       use module_cam_support,    only: pcnst => pcnst_runtime
       use constituents, only: cnst_name

       implicit none

       
       integer m, laptr, lcptr
       character*(*) txtdum

       
       character*8 dumnamea, dumnamec

       dumnamea = 'none'
       dumnamec = 'none'
       if (laptr .gt. 0) dumnamea = cnst_name(laptr)
       if (lcptr .gt. 0) dumnamec = cnst_name(lcptr)
       write(iulog,9241) m, laptr, dumnamea, lcptr, dumnamec, txtdum

9241   format( i4, 2( 2x, i12, 2x, a ),                                &
            4x, 'lptr_', a, '_a/cw_amode' )

       return
     end subroutine initaermodes_setspecptrs_write2


     
     subroutine initaermodes_set_cnstnamecw
       
       
       
       use module_cam_support,    only: pcnst => pcnst_runtime
       use constituents, only: cnst_name
       implicit none

       

       
       integer j, l, la, lc, ll, m

       
       cnst_name_cw = ' '
       do m = 1, ntot_amode
          do ll = 0, nspec_amode(m)
             if (ll == 0) then
                la = numptr_amode(m)
                lc = numptrcw_amode(m)
             else
                la = lmassptr_amode(ll,m)
                lc = lmassptrcw_amode(ll,m)
             end if
             if ((la < 1) .or. (la > pcnst) .or.   &
                  (lc < 1) .or. (lc > pcnst)) then
                write(*,'(/2a/a,5(1x,i10))')   &
                     '*** initaermodes_set_cnstnamecw error',   &
                     ' -- bad la or lc',   &
                     '    m, ll, la, lc, pcnst =', m, ll, la, lc, pcnst
                call endrun( '*** initaermodes_set_cnstnamecw error' )
             end if
             do j = 2, len( cnst_name(la) ) - 1
                if (cnst_name(la)(j:j+1) == '_a') then
                   cnst_name_cw(lc) = cnst_name(la)
                   cnst_name_cw(lc)(j:j+1) = '_c'
                   exit
                else if (cnst_name(la)(j:j+1) == '_A') then
                   cnst_name_cw(lc) = cnst_name(la)
                   cnst_name_cw(lc)(j:j+1) = '_C'
                   exit
                end if
             end do
             if (cnst_name_cw(lc) == ' ') then
                write(*,'(/2a/a,3(1x,i10),2x,a)')   &
                     '*** initaermodes_set_cnstnamecw error',   &
                     ' -- bad cnst_name(la)',   &
                     '    m, ll, la, cnst_name(la) =',   &
                     m, ll, la, cnst_name(la)
                call endrun( '*** initaermodes_set_cnstnamecw error' )
             end if
          end do   
       end do   

       if ( masterproc ) then
          write(*,'(/a)') 'l, cnst_name(l), cnst_name_cw(l)'
          do l = 1, pcnst
             write(*,'(i4,2(2x,a))') l, cnst_name(l), cnst_name_cw(l)
          end do
       end if

       return
     end subroutine initaermodes_set_cnstnamecw


     subroutine decouple_mam_mp(CAM_MP_MAM_cpled)











       implicit none
       logical, intent(in):: CAM_MP_MAM_cpled
       integer, parameter :: init_val=-999888777
       integer  :: n
       real(r8) :: pi, tmpsg_mp(ntot_amode)
       
       modeptr_accum_mp         = modeptr_accum
       modeptr_coarse_mp        = modeptr_coarse
       modeptr_coardust_mp      = modeptr_coardust 
       modeptr_aitken_mp        = modeptr_aitken

       numptrcw_amode_mp(:)     = numptrcw_amode(:)
       nspec_amode_mp(:)        = nspec_amode(:) 
       numptr_amode_mp(:)       = numptr_amode(:) 
       
       cnst_name_cw_mp(:)       = cnst_name_cw(:)
       sigmag_amode_mp(:)       = sigmag_amode(:)
       dgnum_amode_mp(:)        = dgnum_amode(:)
       
       voltonumbhi_amode_mp(:)  = voltonumbhi_amode(:)
       voltonumblo_amode_mp (:) = voltonumblo_amode(:)
       
       lptr_dust_a_amode_mp(:)  = lptr_dust_a_amode(:)
       lptr_nacl_a_amode_mp(:)  = lptr_nacl_a_amode(:)
       specdens_amode_mp(:)     = specdens_amode(:)
       alnsg_amode_mp(:)        = alnsg_amode(:)
       specmw_amode_mp(:)       = specmw_amode(:)
       spechygro_mp(:)          = spechygro(:)

       lmassptrcw_amode_mp(:,:) = lmassptrcw_amode(:,:) 
       lmassptr_amode_mp(:,:)   = lmassptr_amode(:,:) 
       lspectype_amode_mp(:,:)  = lspectype_amode(:,:)  

       if(.NOT.CAM_MP_MAM_cpled) then
          
          pi = 4.*atan(1._r8)
          
          
          
          
          
          modeptr_accum_mp  = 1
          modeptr_aitken_mp = 2
          modeptr_coarse_mp = 3
          
          
          
          nspec_amode_mp(:) = init_val
          lspectype_amode_mp(:,:) = init_val
          lmassptr_amode_mp(:,:) = init_val
          numptr_amode_mp(:) = init_val
          lptr_dust_a_amode_mp(:) = init_val
          lptr_nacl_a_amode_mp(:) = init_val
          
          n = modeptr_accum_mp
          nspec_amode_mp(n) = 1
          lspectype_amode_mp(1,n) = 1  
          lmassptr_amode_mp(1,n) = 6   
          numptr_amode_mp(n) = 7   
          
          n = modeptr_aitken_mp
          nspec_amode_mp(n) = 1
          lspectype_amode_mp(1,n) = 1  
          lmassptr_amode_mp(1,n) = 8   
          numptr_amode_mp(n) = 9   
          
          n = modeptr_coarse_mp
          nspec_amode_mp(n) = 2
          lspectype_amode_mp(1,n) = 2  
          lspectype_amode_mp(2,n) = 3  
          lmassptr_amode_mp(1,n) = 10  
          lmassptr_amode_mp(2,n) = 11  
          numptr_amode_mp(n) = 12  
          lptr_dust_a_amode_mp(n) = lmassptr_amode_mp(1,n)
          lptr_nacl_a_amode_mp(n) = lmassptr_amode_mp(2,n)
          
          lmassptrcw_amode_mp = lmassptr_amode_mp
          numptrcw_amode_mp = numptr_amode_mp
          
          msectional_mp = 0
          alnsg_amode_mp(:) = log( sigmag_amode_mp(:) )
          tmpsg_mp = exp( 4.5 * (alnsg_amode_mp(:)**2) )
          
          voltonumb_amode_mp(  :) = 1.0/( (pi/6.0) * (dgnum_amode_mp(  :)**3) * tmpsg_mp )
          voltonumblo_amode_mp(:) = 1.0/( (pi/6.0) * (dgnumlo_amode_mp(:)**3) * tmpsg_mp )
          voltonumbhi_amode_mp(:) = 1.0/( (pi/6.0) * (dgnumhi_amode_mp(:)**3) * tmpsg_mp )
          
          specdens_amode_mp(:) = 1.0e3   
          specmw_amode_mp(:) = 132.0     
          spechygro_mp(:) = 0.5          
       endif
       

     end subroutine decouple_mam_mp


     
   end module modal_aero_initialize_data

