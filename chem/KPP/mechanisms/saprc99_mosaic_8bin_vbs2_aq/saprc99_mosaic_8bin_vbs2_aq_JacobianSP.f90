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
! File                 : saprc99_mosaic_8bin_vbs2_aq_JacobianSP.f90
! Time                 : Tue Apr 12 23:43:48 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/saprc99_mosaic_8bin_vbs2_aq
! Equation file        : saprc99_mosaic_8bin_vbs2_aq.kpp
! Output root filename : saprc99_mosaic_8bin_vbs2_aq
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE saprc99_mosaic_8bin_vbs2_aq_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  2,  3,  4,  5,  6,  7, &
       7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, &
       7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8, &
       8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, &
       9,  9,  9,  9, 10, 10, 10, 11, 11, 11, 11, 11, &
      12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, &
      13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, &
      15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 17, &
      17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, &
      20, 21, 21, 22, 22, 22, 22, 23, 23, 24, 24, 25, &
      25, 26, 26, 26, 27, 27, 28, 28, 28, 29, 29, 30, &
      30, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34, &
      35, 35, 35, 36, 36, 36, 37, 37, 38, 38, 38, 38, &
      38, 38, 39, 39, 40, 40, 40, 41, 41, 41, 42, 42, &
      43, 43, 44, 44, 44, 44, 45, 45, 45, 45, 46, 46, &
      46, 46, 47, 47, 47, 47, 47, 48, 48, 49, 49, 49, &
      49, 50, 50, 50, 50, 50, 51, 51, 52, 52, 52, 52, &
      53, 53, 53, 53, 54, 54, 54, 54, 54, 55, 55, 56, &
      56, 56, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, &
      59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, &
      60, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, &
      62, 62, 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 67, 67, 67, 67, 67, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 70, &
      70, 70, 70, 70, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 72, 72, &
      72, 72, 72, 72, 72, 73, 73, 73, 73, 73, 74, 74, &
      74, 74, 74, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 76, 76, 76, 76, 76, 76, 76, 77, 77, &
      77, 77, 77, 78, 78, 78, 78, 78, 78, 78, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 88, 88, 88, 88, 88, 88, 88, 88 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97 /)
  INTEGER, PARAMETER, DIMENSION(41) :: LU_IROW_3 = (/ &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98 /)
  INTEGER, PARAMETER, DIMENSION(1121) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 30, 95,  2, 56, 67, 86,  3,  4,  5,  6,  7, &
      46, 56, 65, 67, 69, 70, 72, 73, 74, 76, 77, 78, &
      86, 88, 95,  8, 67, 69, 77, 86, 89, 91, 94, 96, &
      97,  9, 69, 70, 73, 74, 76, 77, 86, 87, 90, 91, &
      94, 96, 97, 98, 10, 89, 96, 11, 87, 90, 96, 98, &
      12, 47, 67, 92, 93, 13, 47, 60, 67, 82, 86, 92, &
      93, 95, 14, 71, 88, 91, 93, 94, 97, 15, 71, 91, &
      94, 96, 97, 16, 42, 48, 51, 55, 69, 77, 95, 17, &
      70, 73, 74, 86, 93, 95, 18, 29, 30, 31, 37, 39, &
      41, 42, 43, 45, 48, 50, 51, 52, 53, 54, 55, 56, &
      57, 58, 60, 61, 62, 63, 64, 65, 66, 67, 69, 70, &
      72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 83, 84, &
      85, 86, 88, 92, 93, 95, 96,  3,  4,  5,  6, 19, &
      20, 21, 22, 23, 25, 26, 27, 28, 95, 20, 21, 25, &
      95, 21, 95, 22, 23, 27, 95, 23, 95, 24, 86, 25, &
      95, 25, 26, 95, 27, 95, 27, 28, 95, 29, 95, 30, &
      95, 31, 95, 32, 89, 92, 33, 92, 98, 34, 90, 92, &
      35, 87, 92, 36, 95, 96, 37, 95, 38, 48, 73, 74, &
      86, 95, 39, 95, 40, 92, 93, 41, 88, 95, 42, 95, &
      43, 95, 43, 44, 92, 95, 45, 94, 95, 96, 46, 80, &
      88, 96, 47, 59, 92, 93, 96, 48, 95, 49, 92, 95, &
      96, 50, 91, 94, 95, 97, 51, 95, 48, 51, 52, 95, &
      48, 51, 53, 95, 48, 51, 54, 93, 95, 55, 95, 56, &
      86, 95, 48, 51, 57, 86, 95, 58, 91, 95, 96, 97, &
      47, 59, 68, 92, 93, 96, 48, 51, 60, 77, 86, 93, &
      95, 51, 61, 68, 93, 95, 96, 48, 51, 52, 53, 54, &
      62, 72, 76, 78, 86, 93, 95, 52, 53, 55, 56, 57, &
      62, 63, 65, 67, 69, 70, 72, 73, 74, 75, 76, 77, &
      78, 79, 80, 82, 83, 86, 93, 95, 40, 54, 59, 60, &
      61, 62, 64, 68, 72, 75, 76, 77, 78, 79, 80, 83 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      86, 92, 93, 95, 96, 65, 82, 86, 93, 95, 37, 42, &
      43, 44, 55, 66, 69, 73, 74, 77, 81, 86, 92, 93, &
      95, 67, 82, 86, 93, 95, 54, 61, 68, 87, 88, 89, &
      90, 92, 93, 95, 96, 98, 69, 82, 86, 93, 95, 70, &
      82, 86, 93, 95, 42, 43, 52, 53, 55, 66, 69, 70, &
      71, 73, 74, 77, 78, 81, 82, 84, 85, 86, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 70, 72, &
      77, 82, 86, 93, 95, 73, 82, 86, 93, 95, 74, 82, &
      86, 93, 95, 48, 51, 52, 53, 56, 57, 61, 65, 68, &
      73, 74, 75, 76, 82, 86, 87, 88, 89, 90, 92, 93, &
      95, 96, 98, 70, 76, 77, 82, 86, 93, 95, 77, 82, &
      86, 93, 95, 70, 77, 78, 82, 86, 93, 95, 31, 39, &
      42, 43, 55, 65, 67, 69, 76, 77, 79, 81, 82, 83, &
      84, 85, 86, 87, 88, 89, 90, 93, 95, 98, 39, 42, &
      43, 45, 46, 50, 55, 56, 65, 66, 67, 69, 70, 72, &
      73, 74, 75, 76, 77, 78, 80, 81, 82, 84, 85, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      44, 69, 73, 74, 76, 77, 81, 82, 86, 88, 92, 93, &
      95, 97, 24, 65, 67, 69, 70, 72, 73, 74, 77, 78, &
      82, 86, 88, 92, 93, 95, 37, 42, 43, 52, 53, 55, &
      57, 58, 65, 67, 69, 72, 73, 74, 76, 77, 78, 81, &
      82, 83, 84, 85, 86, 88, 91, 92, 93, 95, 96, 97, &
      42, 43, 55, 67, 69, 72, 76, 77, 78, 81, 82, 84, &
      85, 86, 88, 91, 92, 93, 94, 95, 97, 42, 51, 55, &
      69, 70, 73, 74, 76, 77, 78, 81, 82, 85, 86, 88, &
      89, 90, 91, 92, 93, 94, 95, 97, 98, 56, 57, 65, &
      67, 69, 70, 72, 73, 74, 76, 77, 78, 82, 86, 87, &
      88, 89, 90, 92, 93, 95, 96, 98, 35, 70, 72, 76, &
      77, 78, 82, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 41, 46, 71, 73, 74, 77, 78, 80 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
      81, 82, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 32, 38, 48, 52, 53, 55, 62, &
      66, 69, 72, 73, 74, 76, 77, 78, 79, 81, 82, 83, &
      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 98, 34, 60, 77, 82, 86, 87, 88, 89, 90, &
      91, 92, 93, 94, 95, 96, 97, 98, 31, 37, 39, 42, &
      43, 48, 51, 52, 53, 54, 55, 56, 57, 58, 61, 65, &
      67, 68, 69, 70, 72, 73, 74, 76, 77, 78, 81, 82, &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 32, 33, 34, 35, 40, 41, 44, 46, &
      47, 49, 59, 64, 68, 70, 71, 72, 73, 74, 75, 76, &
      77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 40, 49, &
      54, 59, 60, 61, 62, 64, 65, 67, 68, 69, 70, 71, &
      72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, &
      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 98, 29, 44, 45, 55, 65, 66, 67, 69, 70, &
      73, 74, 77, 78, 79, 81, 82, 83, 84, 85, 86, 87, &
      88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98,  3, &
       5, 24, 25, 26, 27, 28, 29, 30, 31, 36, 37, 39, &
      41, 42, 43, 45, 48, 49, 50, 51, 52, 53, 54, 55, &
      56, 57, 58, 60, 61, 62, 63, 64, 65, 66, 67, 68, &
      69, 70, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, &
      82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 30, 36, 39, 41, 45, 46, 47, &
      48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 62, 63, &
      65, 67, 68, 69, 70, 72, 73, 74, 75, 76, 77, 78, &
      79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
      91, 92, 93, 94, 95, 96, 97, 98, 37, 42, 43, 48, &
      51, 55, 67, 69, 70, 73, 74, 76, 77, 78, 81, 82 /)
  INTEGER, PARAMETER, DIMENSION(41) :: LU_ICOL_3 = (/ &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 33, 72, 73, 74, 75, 76, 77, 78, &
      82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98 /)
  INTEGER, PARAMETER, DIMENSION(1121) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3 /)

  INTEGER, PARAMETER, DIMENSION(99) :: LU_CROW = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 96,103,152,166,170,172,176,178, &
     180,182,185,187,190,192,194,196,199,202,205,208, &
     211,213,219,221,224,227,229,231,235,239,243,248, &
     250,254,259,261,265,269,274,276,279,284,289,295, &
     302,308,320,345,366,371,386,391,403,408,413,443, &
     450,455,460,484,491,496,503,527,565,579,595,625, &
     646,670,693,713,738,772,789,833,875,916,948,1014, &
     1065,1097,1122 /)

  INTEGER, PARAMETER, DIMENSION(99) :: LU_DIAG = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 96,103,156,166,170,172,176,178, &
     180,183,185,188,190,192,194,196,199,202,205,208, &
     211,213,219,221,224,227,229,232,235,239,243,248, &
     250,254,259,263,267,271,274,276,281,284,290,297, &
     303,313,326,351,366,376,386,393,403,408,421,444, &
     450,455,471,485,491,498,513,547,571,589,614,636, &
     658,683,701,727,762,780,825,868,910,943,1010,1062, &
     1095,1121,1122 /)


END MODULE saprc99_mosaic_8bin_vbs2_aq_JacobianSP

