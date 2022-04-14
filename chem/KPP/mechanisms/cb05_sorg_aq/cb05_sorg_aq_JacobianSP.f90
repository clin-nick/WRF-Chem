! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Jacobian Data Structures File
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
! File                 : cb05_sorg_aq_JacobianSP.f90
! Time                 : Tue Apr 12 23:42:40 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/cb05_sorg_aq
! Equation file        : cb05_sorg_aq.kpp
! Output root filename : cb05_sorg_aq
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cb05_sorg_aq_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  3,  3,  3,  3,  4,  4, &
       4,  5,  5,  5,  6,  6,  6,  6,  6,  6,  7,  7, &
       7,  7,  7,  7,  8,  8,  8,  9,  9,  9, 10, 10, &
      10, 10, 10, 10, 11, 11, 11, 12, 12, 13, 13, 13, &
      14, 14, 14, 15, 15, 16, 16, 16, 17, 17, 17, 18, &
      18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, &
      22, 23, 23, 23, 24, 24, 24, 25, 25, 25, 26, 26, &
      26, 27, 27, 27, 28, 28, 28, 29, 29, 29, 29, 30, &
      30, 30, 31, 31, 31, 32, 32, 33, 33, 33, 34, 34, &
      35, 35, 35, 36, 36, 36, 37, 37, 38, 38, 38, 38, &
      38, 39, 39, 40, 40, 41, 41, 41, 42, 42, 42, 43, &
      43, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46, 47, &
      47, 48, 48, 48, 49, 49, 49, 49, 50, 50, 50, 50, &
      51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 53, &
      53, 53, 53, 54, 54, 54, 54, 54, 55, 55, 55, 55, &
      55, 55, 55, 56, 56, 56, 56, 57, 57, 57, 58, 58, &
      58, 58, 59, 59, 59, 60, 60, 60, 60, 60, 60, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 62, 62, 62, 62, 63, 63, 63, 63, 63, 63, 64, &
      64, 64, 64, 64, 65, 65, 65, 65, 65, 66, 66, 66, &
      66, 66, 66, 67, 67, 67, 67, 67, 67, 67, 68, 68, &
      68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 71, &
      71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 74, 74, 74, 74, &
      74, 74, 75, 75, 75, 75, 75, 76, 76, 76, 76, 76, &
      76, 77, 77, 77, 77, 77, 77, 77, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      79, 79, 79, 80, 80, 80, 80, 80, 80, 80, 81, 81, &
      81, 81, 81, 81, 81, 81, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 95, 95 /)
  INTEGER, PARAMETER, DIMENSION(99) :: LU_IROW_2 = (/ &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97 /)
  INTEGER, PARAMETER, DIMENSION(819) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 43, 96,  2, 43, 96,  3, 73, 96, 97,  4, 47, &
      96,  5, 47, 96,  6, 80, 86, 93, 96, 97,  7, 80, &
      86, 93, 96, 97,  8, 40, 96,  9, 40, 96, 10, 75, &
      86, 93, 96, 97, 11, 12, 96, 12, 96, 13, 15, 96, &
      14, 15, 96, 15, 96, 16, 18, 96, 17, 18, 96, 18, &
      96, 19, 23, 96, 20, 23, 96, 21, 23, 86, 22, 23, &
      86, 23, 86, 96, 24, 29, 96, 25, 29, 96, 26, 29, &
      86, 27, 29, 86, 28, 29, 97, 29, 86, 96, 97, 30, &
      32, 96, 31, 32, 96, 32, 96, 33, 34, 96, 34, 96, &
      35, 37, 96, 36, 37, 96, 37, 96, 38, 49, 60, 86, &
      96, 39, 64, 40, 96, 41, 87, 91, 42, 64, 90, 43, &
      96, 44, 45, 96, 44, 45, 86, 96, 46, 87, 97, 47, &
      96, 48, 88, 96, 49, 60, 86, 96, 50, 87, 95, 96, &
      51, 62, 74, 86, 95, 96, 43, 47, 52, 95, 96, 53, &
      87, 90, 96, 54, 90, 91, 94, 96, 55, 89, 90, 91, &
      92, 94, 96, 56, 87, 94, 96, 57, 88, 96, 58, 88, &
      92, 96, 59, 88, 96, 49, 60, 86, 90, 93, 96, 48, &
      57, 58, 59, 61, 76, 79, 80, 82, 83, 85, 88, 92, &
      96, 62, 85, 90, 95, 63, 73, 87, 90, 96, 97, 64, &
      86, 88, 90, 95, 65, 78, 89, 90, 96, 47, 66, 71, &
      81, 86, 96, 67, 74, 76, 77, 80, 88, 96, 62, 68, &
      85, 90, 92, 95, 96, 46, 69, 73, 81, 82, 83, 84, &
      85, 87, 90, 96, 97, 66, 67, 70, 71, 74, 75, 76, &
      77, 80, 81, 82, 83, 85, 86, 88, 93, 96, 97, 52, &
      71, 73, 86, 95, 96, 72, 79, 87, 88, 96, 43, 47, &
      52, 63, 73, 87, 90, 95, 96, 97, 74, 86, 88, 93, &
      96, 97, 75, 86, 93, 96, 97, 76, 86, 88, 93, 96, &
      97, 76, 77, 86, 88, 93, 96, 97, 57, 72, 75, 77, &
      78, 79, 80, 86, 87, 88, 89, 90, 93, 95, 96, 97, &
      47, 72, 75, 76, 77, 79, 80, 81, 84, 86, 87, 88 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      93, 96, 97, 80, 86, 87, 88, 93, 96, 97, 80, 81, &
      86, 87, 88, 93, 96, 97, 59, 65, 71, 72, 73, 74, &
      75, 76, 77, 78, 79, 80, 81, 82, 84, 86, 87, 88, &
      89, 90, 93, 95, 96, 97, 56, 57, 59, 65, 72, 76, &
      77, 78, 79, 80, 81, 83, 84, 86, 87, 88, 89, 90, &
      91, 92, 93, 94, 95, 96, 97, 52, 63, 72, 73, 75, &
      78, 79, 80, 81, 84, 86, 87, 88, 89, 90, 93, 95, &
      96, 97, 58, 59, 62, 68, 71, 73, 74, 75, 76, 77, &
      80, 81, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 49, 60, 71, 73, 74, 75, 76, 77, &
      80, 81, 86, 87, 88, 90, 91, 93, 94, 95, 96, 97, &
      41, 46, 50, 52, 53, 56, 62, 63, 64, 69, 72, 73, &
      74, 75, 76, 77, 79, 80, 81, 82, 83, 84, 85, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 39, &
      42, 48, 57, 58, 59, 61, 64, 67, 74, 76, 77, 79, &
      80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, &
      92, 93, 94, 95, 96, 97, 43, 47, 57, 59, 65, 66, &
      68, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, &
      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 40, 43, 44, 45, 47, 51, 52, 53, 57, 58, &
      59, 60, 62, 63, 64, 65, 66, 67, 68, 70, 71, 72, &
      73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, &
      97, 41, 54, 66, 71, 73, 81, 83, 84, 86, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 48, 54, 55, &
      68, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 45, 60, 74, 75, 76, 77, 80, &
      82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 56, 75, 80, 81, 82, 84, 86, 87, &
      88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 50, 52 /)
  INTEGER, PARAMETER, DIMENSION(99) :: LU_ICOL_2 = (/ &
      62, 64, 78, 79, 80, 81, 84, 85, 86, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 40, 42, 43, 44, &
      45, 47, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, &
      59, 60, 61, 62, 64, 65, 66, 67, 68, 69, 70, 71, &
      73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, &
      97, 46, 53, 69, 73, 74, 75, 76, 77, 80, 81, 82, &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97 /)
  INTEGER, PARAMETER, DIMENSION(819) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2 /)

  INTEGER, PARAMETER, DIMENSION(98) :: LU_CROW = (/ &
       1,  4,  7, 11, 14, 17, 23, 29, 32, 35, 41, 44, &
      46, 49, 52, 54, 57, 60, 62, 65, 68, 71, 74, 77, &
      80, 83, 86, 89, 92, 96, 99,102,104,107,109,112, &
     115,117,122,124,126,129,132,134,137,141,144,146, &
     149,153,157,163,168,172,177,184,188,191,195,198, &
     204,218,222,228,233,238,244,251,258,270,288,294, &
     299,309,315,320,326,333,349,364,371,379,403,428, &
     447,473,493,528,559,591,638,658,678,701,719,741, &
     794,820 /)

  INTEGER, PARAMETER, DIMENSION(98) :: LU_DIAG = (/ &
       1,  4,  7, 11, 14, 17, 23, 29, 32, 35, 41, 44, &
      46, 49, 52, 54, 57, 60, 62, 65, 68, 71, 74, 77, &
      80, 83, 86, 89, 92, 96, 99,102,104,107,109,112, &
     115,117,122,124,126,129,132,134,138,141,144,146, &
     149,153,157,165,168,172,177,184,188,191,195,199, &
     208,218,222,228,233,239,244,252,259,272,289,294, &
     303,309,315,320,327,337,354,364,372,392,414,437, &
     460,483,517,549,582,630,651,672,696,715,738,792, &
     819,820 /)


END MODULE cb05_sorg_aq_JacobianSP
