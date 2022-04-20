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
! File                 : racm_soa_vbs_het_JacobianSP.f90
! Time                 : Thu Apr 14 13:14:13 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/racm_soa_vbs_het
! Equation file        : racm_soa_vbs_het.kpp
! Output root filename : racm_soa_vbs_het
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE racm_soa_vbs_het_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  3,  3,  3,  3,  3,  3, &
       3,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4, &
       4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, &
       4,  4,  4,  4,  4,  4,  5,  5,  5,  6,  6,  6, &
       7,  7,  7,  8,  8,  9,  9,  9, 10, 10, 10, 11, &
      11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, &
      16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20, &
      20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 23, &
      24, 24, 24, 25, 25, 25, 25, 26, 26, 26, 26, 27, &
      27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 29, 29, &
      29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, &
      31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, &
      33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 34, 34, &
      34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, &
      35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, &
      36, 36, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, &
      38, 38, 39, 39, 39, 39, 39, 39, 40, 40, 40, 40, &
      40, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42, &
      43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, &
      44, 44, 44, 45, 45, 45, 45, 45, 45, 46, 46, 46, &
      46, 46, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, &
      48, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, &
      50, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, &
      53, 54, 54, 54, 54, 55, 55, 55, 55, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, &
      59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 61, 61, 61 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, 62, &
      62, 62, 62, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 82, 82, 82 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92 /)
  INTEGER, PARAMETER, DIMENSION(227) :: LU_IROW_3 = (/ &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97 /)
  INTEGER, PARAMETER, DIMENSION(1307) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 15, 91,  2, 52, 91,  3, 43, 49, 51, 53, 55, &
      58, 59, 64, 79, 82, 85, 91,  4, 43, 51, 54, 61, &
      62, 67, 70, 74, 76, 77, 79, 81, 82, 83, 84, 85, &
      86, 87, 89, 90, 93, 97,  5,  6, 91,  6,  7, 91, &
       7,  8, 91,  8, 91,  9, 10, 91, 10, 11, 91, 11, &
      12, 91, 12, 91, 13, 32, 14, 85, 15, 91, 16, 67, &
      89, 91, 17, 95, 96, 18, 88, 95, 19, 32, 89, 20, &
      56, 89, 91, 21, 79, 91, 22, 42, 91, 23, 88, 91, &
      24, 88, 91, 25, 89, 90, 91, 26, 56, 91, 95, 27, &
      89, 91, 95, 26, 28, 42, 56, 91, 94, 95, 29, 38, &
      39, 40, 91, 94, 95, 14, 30, 55, 85, 94, 95, 96, &
      31, 57, 89, 91, 95, 96, 32, 85, 88, 89, 94, 33, &
      47, 51, 54, 55, 59, 81, 82, 85, 89, 91, 23, 24, &
      34, 35, 41, 46, 49, 50, 59, 73, 78, 81, 88, 91, &
      35, 54, 81, 82, 85, 88, 91, 36, 53, 59, 81, 82, &
      88, 91, 37, 46, 58, 87, 90, 91, 23, 38, 85, 88, &
      91, 95, 24, 39, 85, 88, 91, 95, 40, 57, 85, 91, &
      95, 41, 81, 82, 88, 91, 42, 59, 67, 91, 94, 96, &
      43, 85, 91, 96, 44, 57, 60, 69, 73, 78, 79, 89, &
      91, 95, 96, 45, 79, 85, 89, 91, 93, 46, 81, 82, &
      88, 91, 47, 85, 91, 96, 48, 58, 85, 91, 93, 95, &
      96, 49, 81, 82, 88, 91, 50, 54, 81, 82, 85, 88, &
      91, 51, 85, 91, 96, 22, 36, 42, 47, 49, 51, 52, &
      53, 54, 55, 56, 58, 59, 60, 64, 67, 69, 73, 78, &
      79, 81, 82, 85, 88, 91, 94, 96, 53, 85, 88, 91, &
      96, 54, 85, 91, 96, 55, 85, 91, 96, 20, 26, 56, &
      59, 64, 85, 89, 91, 94, 95, 31, 38, 39, 40, 57, &
      85, 88, 89, 91, 95, 96, 58, 74, 85, 91, 95, 96, &
      59, 85, 88, 91, 96, 49, 60, 68, 71, 72, 79, 81, &
      82, 85, 86, 88, 90, 91, 93, 94, 96, 61, 82, 89 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      90, 91, 93, 94, 96, 53, 62, 85, 88, 89, 90, 91, &
      93, 94, 96, 51, 63, 85, 89, 90, 91, 93, 94, 96, &
      16, 43, 51, 55, 59, 63, 64, 67, 85, 88, 89, 90, &
      91, 93, 94, 96, 47, 65, 85, 89, 90, 91, 93, 94, &
      96, 21, 41, 43, 47, 54, 61, 65, 66, 70, 75, 76, &
      77, 79, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
      91, 92, 93, 94, 96, 55, 59, 67, 85, 88, 89, 90, &
      91, 93, 94, 96, 40, 57, 68, 85, 88, 89, 90, 91, &
      93, 94, 95, 96, 28, 37, 42, 46, 56, 58, 59, 64, &
      67, 68, 69, 71, 72, 74, 79, 81, 82, 85, 86, 87, &
      88, 89, 90, 91, 93, 94, 95, 96, 43, 54, 70, 81, &
      85, 89, 90, 91, 93, 94, 96, 38, 71, 85, 88, 89, &
      90, 91, 93, 94, 95, 96, 39, 72, 85, 88, 89, 90, &
      91, 93, 94, 95, 96, 22, 25, 37, 42, 43, 45, 46, &
      48, 49, 51, 53, 54, 55, 56, 58, 59, 60, 61, 62, &
      63, 64, 65, 67, 68, 70, 71, 72, 73, 74, 76, 77, &
      79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
      91, 93, 94, 95, 96, 97, 58, 74, 79, 85, 89, 90, &
      91, 93, 94, 95, 96, 31, 57, 58, 63, 65, 71, 72, &
      74, 75, 76, 77, 79, 83, 84, 85, 86, 88, 89, 90, &
      91, 93, 94, 95, 96, 41, 76, 81, 82, 88, 89, 90, &
      91, 93, 94, 96, 46, 77, 81, 82, 88, 89, 90, 91, &
      93, 94, 96, 21, 41, 43, 46, 47, 49, 50, 54, 61, &
      62, 65, 70, 75, 76, 77, 78, 79, 81, 82, 83, 84, &
      85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, &
      97, 71, 72, 79, 85, 88, 89, 90, 91, 93, 94, 95, &
      96, 23, 24, 41, 45, 46, 48, 49, 50, 53, 54, 55, &
      57, 58, 59, 74, 76, 77, 79, 80, 81, 82, 85, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 63, 67, &
      81, 85, 88, 89, 90, 91, 93, 94, 96, 51, 55, 67 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
      81, 82, 85, 88, 89, 90, 91, 93, 94, 96, 43, 47, &
      51, 53, 54, 55, 81, 82, 83, 84, 85, 88, 89, 90, &
      91, 93, 94, 96, 43, 47, 51, 53, 54, 55, 81, 82, &
      83, 84, 85, 88, 89, 90, 91, 93, 94, 96, 30, 38, &
      39, 40, 43, 47, 51, 53, 54, 55, 57, 58, 59, 64, &
      67, 74, 79, 81, 82, 85, 88, 89, 90, 91, 93, 94, &
      95, 96, 49, 75, 76, 77, 79, 81, 82, 83, 84, 85, &
      86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 47, 51, &
      54, 55, 66, 70, 75, 76, 77, 79, 81, 82, 83, 84, &
      85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, &
      13, 18, 19, 23, 24, 32, 34, 35, 36, 41, 46, 49, &
      50, 53, 54, 59, 73, 74, 76, 77, 78, 79, 80, 81, &
      82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 15, 21, 23, 24, 25, 27, 28, 31, &
      32, 33, 36, 37, 38, 39, 40, 41, 42, 45, 46, 47, &
      49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
      61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 72, 73, &
      74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, &
      86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      25, 35, 45, 54, 55, 59, 61, 62, 63, 65, 67, 68, &
      70, 71, 72, 74, 76, 77, 78, 79, 80, 81, 82, 83, &
      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 14, 15, 19, 20, 21, 22, 23, 24, 25, 26, &
      27, 28, 29, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
      41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, &
      53, 54, 55, 56, 57, 58, 59, 60, 64, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
      81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 61, 62, 63, 65, 68, 70, 71, &
      72, 74, 76, 77, 79, 80, 81, 82, 85, 86, 87, 88 /)
  INTEGER, PARAMETER, DIMENSION(227) :: LU_ICOL_3 = (/ &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 37, 45, 46, &
      48, 55, 56, 58, 59, 61, 62, 63, 64, 65, 66, 67, &
      68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 79, 80, &
      81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 29, 30, 32, 38, 39, 40, 55, &
      56, 57, 59, 61, 62, 63, 64, 65, 67, 68, 70, 71, &
      72, 74, 76, 77, 79, 80, 81, 82, 83, 84, 85, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 17, &
      18, 22, 26, 27, 29, 30, 31, 32, 38, 39, 40, 42, &
      44, 48, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, &
      65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 17, 27, 30, 43, &
      44, 47, 48, 51, 53, 54, 55, 57, 58, 59, 60, 61, &
      62, 63, 65, 68, 69, 70, 71, 72, 73, 74, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 47, 50, 51, 54, &
      66, 70, 75, 76, 77, 79, 81, 82, 83, 84, 85, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97 /)
  INTEGER, PARAMETER, DIMENSION(1307) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3 /)

  INTEGER, PARAMETER, DIMENSION(98) :: LU_CROW = (/ &
       1,  4,  7, 20, 43, 46, 49, 52, 54, 57, 60, 63, &
      65, 67, 69, 71, 75, 78, 81, 84, 88, 91, 94, 97, &
     100,104,108,112,119,126,133,139,144,155,169,176, &
     183,189,195,201,206,211,217,221,232,238,243,247, &
     254,259,266,270,297,302,306,310,320,331,337,342, &
     358,366,376,385,401,410,438,449,461,489,500,511, &
     522,571,582,606,617,628,662,674,707,718,731,749, &
     767,795,815,841,881,949,987,1062,1090,1134,1176,1233, &
     1281,1308 /)

  INTEGER, PARAMETER, DIMENSION(98) :: LU_DIAG = (/ &
       1,  4,  7, 20, 43, 46, 49, 52, 54, 57, 60, 63, &
      65, 67, 69, 71, 75, 78, 81, 84, 88, 91, 94, 97, &
     100,104,108,113,119,127,133,139,144,157,169,176, &
     183,190,196,201,206,211,217,221,232,238,243,247, &
     254,259,266,276,297,302,306,312,324,331,337,343, &
     358,367,377,391,402,417,440,451,471,491,501,512, &
     549,572,590,607,618,643,664,692,709,722,739,758, &
     786,805,831,871,940,979,1055,1084,1129,1172,1230,1279, &
     1307,1308 /)


END MODULE racm_soa_vbs_het_JacobianSP

