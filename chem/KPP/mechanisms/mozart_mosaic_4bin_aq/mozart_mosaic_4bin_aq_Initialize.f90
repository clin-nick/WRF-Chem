! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialization File
! 
! Generated by KPP-2.1 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : mozart_mosaic_4bin_aq_Initialize.f90
! Time                 : Thu Apr 14 13:13:45 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/mozart_mosaic_4bin_aq
! Equation file        : mozart_mosaic_4bin_aq.kpp
! Output root filename : mozart_mosaic_4bin_aq
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mozart_mosaic_4bin_aq_Initialize

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialize - function to initialize concentrations
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Initialize ( )


  USE mozart_mosaic_4bin_aq_Global

  INTEGER :: i
  REAL(kind=dp) :: x

  CFACTOR = 1.000000e+00_dp

  x = (0.)*CFACTOR
  DO i = 1, NVAR
    VAR(i) = x
  END DO

  x = (0.)*CFACTOR
  DO i = 1, NFIX
    FIX(i) = x
  END DO

! constant rate coefficients
  RCONST(334) = 4.5e-13
! END constant rate coefficients

! INLINED initializations

! End INLINED initializations

      
END SUBROUTINE Initialize

! End of Initialize function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE mozart_mosaic_4bin_aq_Initialize

