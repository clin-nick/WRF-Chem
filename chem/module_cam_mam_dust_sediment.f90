









module dust_sediment_mod










  use shr_kind_mod,  only: r8=>shr_kind_r8




  private
  public :: dust_sediment_vel, dust_sediment_tend


  real (r8), parameter :: vland  = 2.8_r8            
  real (r8), parameter :: vocean = 1.5_r8            
  real (r8), parameter :: mxsedfac   = 0.99_r8       

contains


  subroutine dust_sediment_vel ( ncol )
  implicit none
  integer :: ncol
  return
end subroutine dust_sediment_vel 



  subroutine dust_sediment_tend ( ncol )
  implicit none
  integer :: ncol
  return
end subroutine dust_sediment_tend


end module dust_sediment_mod
