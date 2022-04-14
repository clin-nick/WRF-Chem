module mo_chem_utls

  private
  public :: get_spc_ndx, get_het_ndx, get_inv_ndx

  save

contains

  integer function get_spc_ndx( spc_name )
    
    
    

    use constituents,        only : cnst_get_ind
    use module_cam_support, only : pcnst_non_chem => pcnst_non_chem_modal_aero

    implicit none

    
    
    
    character(len=*), intent(in) :: spc_name

    
    
    
    integer :: m

    get_spc_ndx = -1
    call cnst_get_ind(trim(adjustl(spc_name)),m,.false.)
    get_spc_ndx = m - pcnst_non_chem 

  end function get_spc_ndx

  integer function get_inv_ndx( invariant )
    
    
    


    implicit none

    
    
    
    character(len=*), intent(in) :: invariant

    
    
    
    integer :: m

    get_inv_ndx = -1
  end function get_inv_ndx


  integer function get_het_ndx( het_name )
    
    
    
    use module_cam_support, only: gas_wetdep_method, gas_wetdep_list, gas_wetdep_cnt

    implicit none

    
    
    
    character(len=*), intent(in) :: het_name

    
    
    
    integer :: m

    get_het_ndx=-1

    do m=1,gas_wetdep_cnt

       if( trim( het_name ) == trim( gas_wetdep_list(m) ) ) then
          get_het_ndx = get_spc_ndx( gas_wetdep_list(m) )
          return
       endif
  
    enddo

  end function get_het_ndx
end module mo_chem_utls
