























MODULE nmhc9_JacobianSP

  PUBLIC
  SAVE





  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  2,  2,  3,  3,  3,  4,  4,  5,  5,  5, &
       6,  6,  6,  7,  7,  7,  8,  8,  8,  8,  8,  9, &
       9,  9,  9, 10, 10, 10, 10, 10, 10, 11, 11, 11, &
      11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, &
      14, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, &
      17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 20, &
      20, 20, 20, 21, 21, 21, 22, 22, 22, 22, 23, 23, &
      23, 23, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, &
      25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 28, 28, &
      28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 30, 30, &
      30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, &
      31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, &
      32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 33, &
      33, 33, 33, 33, 33, 34, 34, 34, 34, 34, 34, 34, &
      35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, &
      35, 35, 36, 36, 36, 36, 36, 36, 36, 36, 37, 37, &
      37, 37, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, &
      39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 40, &
      40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, &
      41, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, &
      42, 42, 42, 42, 43, 43, 43, 43, 43, 43, 43, 44, &
      44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 45, 45, &
      45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, &
      45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, &
      45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46, &
      46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, &
      47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, &
      47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, &
      49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50 /)
  INTEGER, PARAMETER, DIMENSION(220) :: LU_IROW_1 = (/ &
      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, &
      50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55, &
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57 /)
  INTEGER, PARAMETER, DIMENSION(580) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 55,  2, 55,  3, 20, 55,  4, 55,  5, 52, 54, &
       6, 50, 54,  7, 33, 55,  8, 38, 53, 55, 57,  9, &
      51, 55, 57, 10, 47, 50, 51, 55, 57, 11, 48, 55, &
      57, 12, 50, 55, 57, 13, 51, 54, 55, 14, 49, 54, &
      55, 15, 19, 37, 38, 41, 53, 55, 16, 47, 55, 57, &
      17, 34, 55, 57, 18, 34, 55, 56, 19, 53, 55, 20, &
      42, 55, 57, 21, 24, 53, 22, 54, 55, 57, 23, 36, &
      55, 57, 21, 24, 37, 53, 55,  5, 25, 45, 46, 52, &
      54, 55, 26, 40, 55, 57, 27, 43, 55, 57, 17, 18, &
      28, 34, 50, 55, 56, 57, 29, 49, 55, 57, 27, 30, &
      43, 44, 50, 55, 56, 57,  7, 19, 29, 31, 33, 35, &
      37, 38, 41, 45, 46, 49, 50, 52, 53, 55, 56, 57, &
      32, 34, 36, 37, 43, 47, 48, 49, 50, 53, 55, 33, &
      38, 48, 52, 55, 56,  4, 17, 34, 50, 55, 56, 57, &
      23, 29, 35, 36, 37, 39, 41, 48, 49, 50, 53, 55, &
      56, 57, 23, 28, 34, 36, 50, 55, 56, 57, 37, 52, &
      53, 55, 38, 52, 53, 55, 14, 26, 29, 33, 36, 38, &
      39, 40, 48, 49, 50, 52, 53, 54, 55, 56, 57, 19, &
      26, 37, 40, 52, 53, 55, 56, 57, 11, 33, 38, 41, &
      48, 50, 52, 53, 55, 56, 57, 20, 30, 42, 43, 44, &
      50, 55, 56, 57,  2, 27, 43, 50, 55, 56, 57, 37, &
      40, 42, 43, 44, 50, 52, 53, 55, 56, 57,  6,  7, &
      11, 12, 13, 19, 21, 24, 26, 29, 32, 33, 34, 36, &
      37, 38, 39, 40, 41, 43, 45, 47, 48, 49, 50, 51, &
      52, 53, 54, 55, 56, 57, 16, 20, 26, 27, 37, 40, &
      42, 43, 44, 46, 47, 50, 51, 52, 53, 55, 56, 57, &
       1,  4, 16, 27, 30, 37, 43, 44, 47, 50, 51, 52, &
      53, 55, 56, 57, 38, 48, 50, 52, 53, 55, 56, 57, &
      14, 29, 38, 41, 48, 49, 50, 52, 53, 54, 55, 56, &
      57,  6,  9, 10, 12, 21, 24, 28, 34, 36, 37, 38 /)
  INTEGER, PARAMETER, DIMENSION(220) :: LU_ICOL_1 = (/ &
      43, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, &
      57,  3,  9, 13, 20, 23, 28, 29, 30, 34, 35, 36, &
      37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, &
      50, 51, 52, 53, 54, 55, 56, 57,  5,  6, 22, 25, &
      37, 38, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, &
      55, 56, 57, 19, 21, 24, 37, 38, 41, 48, 50, 51, &
      52, 53, 54, 55, 56, 57,  5,  6,  7, 13, 14, 18, &
      22, 25, 33, 34, 36, 38, 40, 42, 43, 44, 45, 46, &
      47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,  1, &
       2,  4,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, &
      17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, &
      29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
      41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, &
      53, 54, 55, 56, 57, 34, 36, 40, 42, 43, 44, 47, &
      48, 49, 50, 51, 52, 53, 54, 55, 56, 57,  8, 11, &
      12, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, &
      29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, &
      42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, &
      54, 55, 56, 57 /)
  INTEGER, PARAMETER, DIMENSION(580) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1 /)

  INTEGER, PARAMETER, DIMENSION(58) :: LU_CROW = (/ &
       1,  3,  5,  8, 10, 13, 16, 19, 24, 28, 34, 38, &
      42, 46, 50, 57, 61, 65, 69, 72, 76, 79, 83, 87, &
      92, 99,103,107,115,119,127,145,156,162,169,183, &
     191,195,199,216,225,236,245,252,263,295,313,329, &
     337,350,374,405,424,439,468,522,539,581 /)

  INTEGER, PARAMETER, DIMENSION(58) :: LU_DIAG = (/ &
       1,  3,  5,  8, 10, 13, 16, 19, 24, 28, 34, 38, &
      42, 46, 50, 57, 61, 65, 69, 72, 76, 79, 83, 88, &
      93, 99,103,109,115,120,130,145,156,164,171,186, &
     191,195,205,219,228,238,247,256,283,304,321,330, &
     342,366,398,418,434,464,519,537,580,581 /)


END MODULE nmhc9_JacobianSP

