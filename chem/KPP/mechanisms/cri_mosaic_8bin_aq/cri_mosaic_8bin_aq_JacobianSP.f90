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
! File                 : cri_mosaic_8bin_aq_JacobianSP.f90
! Time                 : Tue Apr 12 23:42:58 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/cri_mosaic_8bin_aq
! Equation file        : cri_mosaic_8bin_aq.kpp
! Output root filename : cri_mosaic_8bin_aq
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cri_mosaic_8bin_aq_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  2,  3,  4,  4,  4,  5,  5,  5,  5,  6,  6, &
       6,  6,  6,  7,  7,  7,  8,  8,  9,  9, 10, 10, &
      10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, &
      16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, &
      22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, &
      28, 29, 29, 30, 30, 31, 31, 32, 32, 32, 33, 33, &
      34, 34, 34, 35, 35, 35, 36, 36, 37, 37, 37, 38, &
      38, 38, 39, 39, 39, 39, 40, 40, 40, 40, 41, 41, &
      41, 42, 42, 42, 43, 43, 43, 44, 44, 44, 45, 45, &
      45, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 48, &
      48, 48, 48, 48, 49, 49, 49, 49, 50, 50, 50, 50, &
      50, 51, 51, 51, 52, 52, 52, 52, 52, 53, 53, 53, &
      53, 53, 54, 54, 54, 54, 55, 55, 55, 55, 56, 56, &
      56, 56, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59, &
      59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 62, 62, &
      62, 62, 63, 63, 63, 63, 64, 64, 64, 64, 65, 65, &
      65, 65, 66, 66, 66, 66, 67, 67, 67, 67, 68, 68, &
      68, 68, 69, 69, 69, 69, 70, 70, 70, 70, 71, 71, &
      71, 71, 72, 72, 72, 72, 73, 73, 73, 73, 74, 74, &
      74, 74, 75, 75, 75, 75, 76, 76, 76, 76, 77, 77, &
      77, 77, 77, 77, 78, 78, 78, 79, 79, 79, 79, 80, &
      80, 80, 80, 80, 81, 81, 81, 81, 82, 82, 82, 82, &
      83, 83, 83, 83, 84, 84, 84, 84, 85, 85, 85, 85, &
      86, 86, 86, 86, 87, 87, 87, 87, 88, 88, 88, 88, &
      89, 89, 89, 89, 90, 90, 90, 90, 91, 91, 91, 91, &
      92, 92, 92, 92, 93, 93, 93, 93, 94, 94, 94, 94, &
      95, 95, 95, 95, 96, 96, 96, 96, 97, 97, 97, 97, &
      98, 98, 98, 98, 99, 99, 99, 99,100,100,100,100, &
     101,101,101,101,102,102,102,102,103,103,103,103, &
     104,104,104,104,104,105,105,105,105,105,106,106 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
     106,106,107,107,107,107,107,107,107,107,108,108, &
     108,108,108,108,108,108,108,108,109,109,109,110, &
     110,110,110,110,111,111,111,111,111,112,112,112, &
     112,113,113,113,113,114,114,114,114,115,115,115, &
     115,116,116,116,116,117,117,117,117,118,118,118, &
     118,119,119,119,119,120,120,120,120,120,121,121, &
     121,121,122,122,122,122,122,123,123,123,123,124, &
     124,124,124,125,125,125,125,126,126,126,126,127, &
     127,127,127,127,127,128,128,128,128,129,129,129, &
     129,130,130,130,130,130,131,131,131,131,131,132, &
     132,132,132,132,132,132,133,133,133,133,133,133, &
     133,133,133,133,133,133,134,134,134,134,134,135, &
     135,135,135,136,136,136,136,136,137,137,137,137, &
     137,138,138,138,138,138,138,138,138,139,139,139, &
     139,139,140,140,140,140,140,140,140,141,141,141, &
     141,141,141,141,142,142,142,142,142,143,143,143, &
     143,143,143,143,144,144,144,144,144,145,145,145, &
     145,145,146,146,146,146,146,147,147,147,147,147, &
     147,147,147,147,147,147,148,148,148,148,148,149, &
     149,149,149,149,150,150,150,150,151,151,151,151, &
     151,151,151,151,151,151,151,151,151,151,152,152, &
     152,152,152,152,152,152,153,153,153,153,153,153, &
     153,153,154,154,154,154,155,155,155,155,155,155, &
     155,155,155,155,156,156,156,156,156,156,157,157, &
     157,157,157,157,157,157,157,157,157,158,158,158, &
     158,158,158,158,158,159,159,159,159,160,160,160, &
     160,160,160,160,161,161,161,161,162,162,162,162, &
     163,163,163,163,163,163,163,163,163,164,164,164, &
     164,164,164,165,165,165,165,166,166,166,166,166, &
     166,166,166,166,166,166,166,167,167,167,167,167 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
     167,167,168,168,168,168,168,168,168,168,168,168, &
     168,169,169,169,169,169,169,169,169,170,170,170, &
     170,170,170,170,170,170,170,170,171,171,171,171, &
     171,171,171,172,172,172,172,172,172,172,173,173, &
     173,173,173,173,173,173,173,173,174,174,174,174, &
     174,174,174,175,175,175,175,175,175,175,175,176, &
     176,176,176,176,176,177,177,177,177,178,178,178, &
     178,178,179,179,179,179,179,179,179,179,179,179, &
     179,179,180,180,180,180,180,180,180,181,181,181, &
     181,181,181,181,181,181,181,181,181,181,181,182, &
     182,182,182,182,182,183,183,183,183,183,183,183, &
     183,183,183,183,183,183,184,184,184,184,184,184, &
     184,184,184,184,184,185,185,185,185,185,185,185, &
     185,186,186,186,186,186,186,186,186,187,187,187, &
     187,187,187,187,187,187,187,187,187,187,187,187, &
     187,187,187,187,187,187,187,187,187,187,187,187, &
     187,187,187,187,187,187,187,187,187,187,187,187, &
     187,187,187,187,187,188,188,188,188,188,188,188, &
     188,188,188,188,188,189,189,189,189,189,189,190, &
     190,190,190,190,190,190,190,190,190,190,190,190, &
     190,190,190,190,190,190,190,190,190,190,190,190, &
     190,190,190,190,190,190,191,191,191,191,191,191, &
     191,191,191,192,192,192,192,192,192,192,192,192, &
     192,192,192,192,192,192,192,192,192,192,192,192, &
     192,192,192,192,192,193,193,193,193,193,193,193, &
     194,194,194,194,194,194,194,194,195,195,195,195, &
     195,195,195,195,195,195,195,195,195,195,195,195, &
     195,195,195,195,195,195,195,195,195,196,196,196, &
     196,196,196,197,197,197,197,197,197,197,198,198, &
     198,198,198,198,198,198,198,198,198,198,198,198 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_3 = (/ &
     198,198,198,198,198,198,198,198,198,198,199,199, &
     199,199,199,199,199,199,199,199,199,199,199,199, &
     199,200,200,200,200,200,200,200,200,200,201,201, &
     201,201,201,201,201,201,201,201,202,202,202,202, &
     202,202,202,202,202,203,203,203,203,203,203,203, &
     203,203,203,203,203,204,204,204,204,204,204,204, &
     204,204,204,205,205,205,205,205,205,205,205,205, &
     206,206,206,206,206,206,206,206,207,207,207,207, &
     207,207,207,207,207,207,207,207,207,207,207,207, &
     207,207,207,207,207,207,207,207,207,207,207,207, &
     207,207,207,207,207,207,207,207,207,207,207,207, &
     207,207,207,207,207,207,207,207,207,207,208,208, &
     208,208,208,208,208,208,208,208,209,209,209,209, &
     209,209,209,209,209,209,209,210,210,210,210,210, &
     210,210,210,210,210,210,210,210,210,210,210,210, &
     210,210,210,210,210,210,210,210,210,210,210,210, &
     210,210,210,210,210,211,211,211,211,211,211,211, &
     211,211,212,212,212,212,212,212,212,212,212,212, &
     212,213,213,213,213,213,213,213,213,213,213,213, &
     213,213,213,213,213,214,214,214,214,214,214,214, &
     214,214,214,214,214,214,214,214,214,214,214,214, &
     214,214,214,214,214,214,214,214,215,215,215,215, &
     215,215,215,215,215,215,216,216,216,216,216,216, &
     216,216,216,216,216,216,216,216,216,216,216,216, &
     216,216,216,216,216,217,217,217,217,217,217,217, &
     217,217,218,218,218,218,218,218,218,218,218,218, &
     218,218,218,218,218,218,218,218,218,218,218,218, &
     218,218,218,218,219,219,219,219,219,219,219,219, &
     219,219,219,219,219,219,219,219,219,219,219,219, &
     219,219,219,219,219,219,219,219,219,219,219,219 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_4 = (/ &
     219,219,219,219,219,219,219,219,219,219,219,219, &
     219,219,219,220,220,220,220,220,220,220,220,221, &
     221,221,221,221,221,221,221,221,221,221,221,221, &
     221,221,221,221,221,221,221,221,221,221,221,221, &
     221,221,221,221,221,222,222,222,222,222,222,222, &
     222,222,222,222,223,223,223,223,223,223,223,223, &
     223,223,223,223,223,223,223,223,223,223,224,224, &
     224,224,224,224,224,224,224,224,224,224,224,224, &
     224,225,225,225,225,225,225,225,225,225,225,225, &
     225,225,225,225,226,226,226,226,226,226,226,226, &
     226,226,226,226,226,226,226,226,226,226,226,226, &
     226,226,227,227,227,227,227,227,227,227,227,227, &
     227,227,227,227,227,227,227,227,227,227,227,227, &
     227,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,228,228,228,228,228,228,228, &
     228,228,228,228,228,229,229,229,229,229,229,229, &
     229,229,229,229,229,229,229,229,229,229,229,229, &
     229,229,229,229,229,229,229,229,229,229,229,229, &
     229,229,229,229,229,229,229,229,229,229,229,229, &
     229,229,229,229,229,229,229,229,229,229,229,229, &
     229,229,229,229,229,229,229,229,229,229,229,229 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_5 = (/ &
     229,229,229,230,230,230,230,230,230,230,230,230, &
     230,230,230,230,230,230,230,230,230,230,230,230, &
     230,230,230,230,230,230,230,230,230,230,230,230, &
     230,230,230,230,230,230,230,230,230,230,230,230, &
     230,230,230,230,230,230,230,230,230,230,230,230, &
     230,230,230,230,230,230,230,230,230,230,230,230, &
     230,230,230,230,230,230,230,230,230,230,230,230, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,231,231,231,231,231,231,231,231,231,231,231, &
     231,232,232,232,232,232,232,232,232,232,232,232, &
     232,232,232,232,232,232,232,232,232,232,232,232, &
     232,232,232,232,232,232,232,232,232,232,232,232, &
     232,232,232,232,232,232,232,232,232,232,232,232, &
     232,232,232,232,232,232,232,232,232,232,232,232 /)
  INTEGER, PARAMETER, DIMENSION(94) :: LU_IROW_6 = (/ &
     232,232,232,232,232,232,232,232,232,232,232,232, &
     232,232,232,232,232,232,232,232,232,232,232,232, &
     232,232,232,232,232,232,232,232,232,232,232,232, &
     232,232,232,232,232,232,232,232,232,232,232,232, &
     232,232,232,232,232,232,232,232,233,233,233,233, &
     233,233,233,233,233,233,233,233,233,233,233,233, &
     233,233,233,233,233,233,233,233,233,233,233,233, &
     233,233,233,233,233,233,233,233,233,233 /)
  INTEGER, PARAMETER, DIMENSION(2254) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3, LU_IROW_4, &
    LU_IROW_5, LU_IROW_6 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1,  2,  3,  4,  5, 53,  5, 10,127,173,  6, 53, &
     110,228,231,  7, 42,231,  8,223,  9,231, 10,127, &
     231, 11,231, 12,231, 13,231, 14,231, 15,231, 16, &
     231, 17,231, 18,231, 19,231, 20,231, 21,231, 22, &
     231, 23,231, 24,231, 25,231, 26,231, 27,231, 28, &
     231, 29,231, 30,231, 31,231, 32,230,232, 33,231, &
      34, 97,231, 35,210,232, 36,231, 37,210,231, 38, &
      81,231, 39,102,103,231, 40,165,223,231, 41,229, &
     231, 42,109,231, 43,214,231, 44,156,231, 45,116, &
     231, 36, 46,154,177,223,231, 30, 47,230,231, 47, &
      48,230,231,232, 31, 49,230,231, 49, 50,230,231, &
     232, 51,216,231, 52,150,161,223,231, 53,171,223, &
     228,232, 54,218,229,231, 55,213,229,231, 56,216, &
     228,231, 57,206,228,231, 58,221,228,231, 59,208, &
     228,231, 60,148,228,231, 61,145,228,231, 62,204, &
     228,231, 63,228,231,232, 64,180,229,231, 65,144, &
     228,231, 66,198,229,231, 67,193,229,231, 68,193, &
     228,231, 69,222,228,231, 70,146,228,231, 71,212, &
     228,231, 72,158,228,231, 73,214,228,231, 74,156, &
     228,231, 75,215,228,231, 76,219,228,231, 77,153, &
     155,170,183,231, 78,162,231, 79,224,229,231, 80, &
     178,224,228,231, 81,220,228,231, 82,213,228,231, &
      83,167,228,231, 84,156,229,231, 85,188,228,231, &
      86,200,228,231, 87,201,231,232, 88,186,228,231, &
      89,164,228,231, 90,210,229,231, 91,214,229,231, &
      92,182,228,231, 93,210,228,231, 94,216,229,231, &
      95,205,228,231, 96,184,228,231, 97,212,229,231, &
      98,206,229,231, 99,185,229,231,100,215,229,231, &
     101,198,228,231,102,211,229,231,103,211,228,231, &
     104,209,225,229,231,105,174,202,229,231,106,176 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
     228,231,107,159,165,191,203,223,228,231, 75,100, &
     108,159,215,223,228,229,230,231,109,230,231, 42, &
     109,110,230,231,109,111,229,230,231,112,204,231, &
     232,113,188,231,232,114,180,228,231,115,201,228, &
     231,116,217,228,231,117,200,231,232,118,218,228, &
     231,119,164,229,231, 48,120,230,231,232,121,182, &
     229,231, 50,122,230,231,232,123,222,231,232,124, &
     194,228,231,125,185,228,231,126,169,228,231,127, &
     171,172,173,231,232,128,176,229,231,129,219,231, &
     232,130,174,202,228,231,131,189,196,228,231, 57, &
      98,132,206,228,229,231, 99,125,126,133,159,169, &
     185,223,228,229,230,231,134,157,197,228,231,135, &
     226,229,231, 78,136,162,230,231,137,209,225,228, &
     231, 79, 80,138,178,224,228,229,231,139,189,196, &
     229,231, 55, 82,140,213,228,229,231,135,141,175, &
     226,229,230,231,142,157,197,229,231,111,143,223, &
     229,230,231,232,144,177,228,229,230,145,161,228, &
     229,230,146,165,228,229,230, 59,140,147,196,197, &
     208,213,228,229,230,231,148,154,228,229,230,149, &
     215,228,229,230,150,223,230,231, 36, 89, 92,119, &
     121,151,164,182,189,196,228,229,230,231, 65,144, &
     152,177,228,229,230,231,106,128,153,176,228,229, &
     230,231,154,223,230,231, 89,119,155,164,182,196, &
     228,229,230,231, 26,156,228,229,230,231, 21, 22, &
      23, 24, 28, 29,157,228,229,230,231, 71, 72,158, &
     212,228,229,230,231,159,223,230,231,124,160,194, &
     228,229,230,231,161,223,230,231,162,205,229,231, &
      88,152,163,177,186,228,229,230,231, 30,164,228, &
     229,230,231,165,223,230,231, 67, 68, 70,146,165, &
     166,193,223,228,229,230,231,150,167,223,228,229 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
     230,231,104,120,137,168,209,225,228,229,230,231, &
     232,126,159,169,223,228,229,230,231,134,142,157, &
     170,189,196,197,228,229,230,231,110,171,172,223, &
     230,231,232,143,172,223,229,230,231,232,  8,127, &
     171,172,173,223,229,230,231,232, 98,174,206,228, &
     229,230,231,149,175,215,226,228,229,230,231, 25, &
     176,228,229,230,231,177,223,230,231,178,220,228, &
     229,230, 54,117,118,179,200,201,218,228,229,230, &
     231,232,177,180,223,228,229,230,231,105,122,123, &
     130,174,181,202,206,222,228,229,230,231,232, 31, &
     182,228,229,230,231, 92,121,131,139,182,183,189, &
     196,197,228,229,230,231, 96,136,160,162,184,194, &
     205,228,229,230,231,125,159,185,223,228,229,230, &
     231,152,177,186,223,228,229,230,231, 36, 88, 96, &
     112,113,117,129,136,150,151,152,154,159,160,161, &
     162,164,177,182,184,186,187,188,189,190,191,194, &
     196,199,200,201,203,204,205,207,219,223,227,228, &
     229,230,231,232,233, 78, 85,113,162,188,199,205, &
     228,229,230,231,232, 27,189,228,229,230,231, 58, &
     106,115,128,131,134,139,142,157,163,176,177,179, &
     182,186,189,190,191,196,197,200,201,203,218,221, &
     223,228,229,230,231,232, 64,114,180,191,223,228, &
     229,230,231, 47, 49,109,110,120,122,136,160,162, &
     166,191,192,193,194,199,203,205,207,223,227,228, &
     229,230,231,232,233,165,193,223,228,229,230,231, &
      95,124,194,205,228,229,230,231, 44, 74, 84, 95, &
      97,102,103,123,156,162,165,185,193,195,205,211, &
     212,215,222,223,228,229,230,231,232, 27,196,228, &
     229,230,231, 28, 29,197,228,229,230,231, 16, 17, &
      18, 19, 43, 77,153,154,155,164,170,176,182,183 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_3 = (/ &
     189,196,197,198,214,223,228,229,230,231, 66, 86, &
     101,115,191,198,199,200,201,214,223,228,229,230, &
     231, 86,117,200,203,228,229,230,231,232, 87,115, &
     191,201,223,228,229,230,231,232, 57,165,202,206, &
     223,228,229,230,231, 87,114,177,180,201,203,223, &
     228,229,230,231,232, 62,112,204,217,227,228,229, &
     230,231,232, 72, 95,158,205,212,228,229,230,231, &
     185,193,206,223,228,229,230,231, 37, 60, 83, 85, &
      90, 93, 99,101,111,113,114,118,124,125,126,129, &
     148,150,154,159,163,167,169,177,179,180,185,186, &
     188,194,198,199,200,201,203,205,207,210,212,214, &
     215,218,219,221,223,228,229,230,231,232, 33,141, &
     175,208,215,226,228,229,230,231, 38, 79, 81,178, &
     209,220,224,228,229,230,231, 12, 13, 14, 20, 35, &
      42, 52, 53, 76, 93,109,150,161,171,172,183,189, &
     195,196,197,205,210,211,212,215,219,222,223,228, &
     229,230,231,232,233,108,159,211,215,223,228,229, &
     230,231, 40, 69, 71,165,212,222,223,228,229,230, &
     231, 39, 45,102,103,116,149,161,211,213,215,217, &
     223,228,229,230,231,  9, 15, 33, 58, 62,135,141, &
     155,164,170,175,182,189,196,197,204,214,215,217, &
     221,226,227,228,229,230,231,232,133,159,169,185, &
     215,223,228,229,230,231, 26, 59,132,138,153,176, &
     178,183,189,196,197,206,208,215,216,220,223,224, &
     226,228,229,230,231,138,178,217,220,224,228,229, &
     230,231, 44, 51, 78,147,150,156,162,196,197,205, &
     208,212,213,215,216,217,218,220,222,223,224,226, &
     228,229,230,231, 33, 45, 76, 86,116,129,138,140, &
     141,147,160,163,175,177,178,179,186,190,191,194, &
     195,196,197,200,201,203,205,208,211,212,213,215 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_4 = (/ &
     217,218,219,220,221,222,223,224,226,228,229,230, &
     231,232,233,132,206,220,223,228,229,230,231, 38, &
      45, 81, 96,116,132,168,181,184,194,195,202,205, &
     206,209,211,212,215,217,220,221,222,223,224,225, &
     228,229,230,231,232, 69,123,166,193,222,223,228, &
     229,230,231,232,143,150,154,159,161,165,171,172, &
     173,177,191,203,223,228,229,230,231,232, 81,153, &
     176,181,202,206,220,222,223,224,228,229,230,231, &
     232, 34, 80, 97,178,212,220,222,223,224,225,228, &
     229,230,231,232, 11,103,116,168,170,189,196,197, &
     209,211,215,217,220,223,224,225,226,228,229,230, &
     231,232, 51, 56, 94,130,137,174,181,202,206,209, &
     216,220,222,223,224,225,226,227,228,229,230,231, &
     232, 10, 30, 31, 36, 37, 43, 44, 46, 51, 53, 56, &
      63, 65, 68, 73, 74, 75, 78, 82, 84, 85, 88, 89, &
      90, 91, 92, 93, 94,100,101,106,107,109,110,114, &
     115,118,119,121,124,125,127,128,130,131,134,135, &
     136,137,139,140,142,144,145,146,148,149,151,152, &
     153,154,155,156,157,158,159,161,162,164,165,166, &
     167,168,169,170,171,172,173,174,175,176,177,178, &
     179,180,181,182,183,184,185,186,187,188,189,190, &
     191,193,194,196,197,198,199,200,201,202,203,204, &
     205,206,207,208,209,210,211,212,213,214,215,216, &
     217,218,219,220,221,222,223,224,225,226,227,228, &
     229,230,231,232,233, 41,111,143,144,145,146,148, &
     149,154,156,157,158,161,164,165,167,169,171,172, &
     173,174,176,177,178,180,182,184,185,186,188,189, &
     193,194,196,197,198,199,200,201,202,203,204,205, &
     206,208,209,210,211,212,213,214,215,216,217,218, &
     219,220,221,222,223,224,225,226,227,228,229,230 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_5 = (/ &
     231,232,233, 32, 47, 49,109,110,120,122,136,144, &
     145,146,148,149,150,154,156,157,158,159,160,161, &
     162,164,165,166,167,169,173,174,176,177,178,180, &
     182,184,185,186,188,189,191,192,193,194,196,197, &
     198,199,200,201,202,203,204,205,206,207,208,209, &
     210,211,212,213,214,215,216,217,218,219,220,221, &
     222,223,224,225,226,227,228,229,230,231,232,233, &
       8,  9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
      21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, &
      34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, &
      47, 49, 51, 52, 54, 55, 56, 57, 58, 59, 60, 61, &
      62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, &
      74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, &
      86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      98, 99,100,101,102,103,104,105,106,107,108,109, &
     110,112,113,114,115,116,117,118,119,120,121,122, &
     123,124,125,126,127,128,129,130,131,132,133,134, &
     135,136,137,138,139,140,141,142,144,145,146,147, &
     148,150,151,152,153,154,155,156,157,158,159,160, &
     161,162,163,164,165,166,167,168,169,170,171,172, &
     173,174,175,176,177,178,179,180,181,182,183,184, &
     185,186,187,188,189,190,191,192,193,194,195,196, &
     197,198,199,200,201,202,203,204,205,206,207,208, &
     209,210,211,212,213,214,215,216,217,218,219,220, &
     221,222,223,224,225,226,227,228,229,230,231,232, &
     233, 32, 35, 41, 48, 50, 54, 55, 60, 61, 63, 64, &
      66, 67, 70, 79, 83, 84, 87, 90, 91, 94, 97, 98, &
      99,100,102,104,105,111,112,113,117,119,120,121, &
     122,123,126,128,129,135,139,142,143,144,145,146, &
     148,149,154,156,157,158,161,162,163,164,165,167 /)
  INTEGER, PARAMETER, DIMENSION(94) :: LU_ICOL_6 = (/ &
     169,171,172,173,174,176,177,178,180,182,184,185, &
     186,188,189,192,193,194,196,197,198,199,200,201, &
     202,203,204,205,206,207,208,209,210,211,212,213, &
     214,215,216,217,218,219,220,221,222,223,224,225, &
     226,227,228,229,230,231,232,233, 43, 61, 73, 82, &
      83, 91,112,118,135,137,140,145,161,167,168,175, &
     204,208,209,213,214,215,217,218,220,221,222,223, &
     224,225,226,227,228,229,230,231,232,233 /)
  INTEGER, PARAMETER, DIMENSION(2254) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3, LU_ICOL_4, &
    LU_ICOL_5, LU_ICOL_6 /)

  INTEGER, PARAMETER, DIMENSION(234) :: LU_CROW = (/ &
       1,  2,  3,  4,  7, 11, 16, 19, 21, 23, 26, 28, &
      30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, &
      54, 56, 58, 60, 62, 64, 66, 68, 71, 73, 76, 79, &
      81, 84, 87, 91, 95, 98,101,104,107,110,116,120, &
     125,129,134,137,142,147,151,155,159,163,167,171, &
     175,179,183,187,191,195,199,203,207,211,215,219, &
     223,227,231,235,239,245,248,252,257,261,265,269, &
     273,277,281,285,289,293,297,301,305,309,313,317, &
     321,325,329,333,337,341,345,349,354,359,363,371, &
     381,384,389,394,398,402,406,410,414,418,422,426, &
     431,435,440,444,448,452,456,462,466,470,475,480, &
     487,499,504,508,513,518,526,531,538,545,550,557, &
     562,567,572,583,588,593,597,611,619,627,631,641, &
     647,658,666,670,677,681,685,694,700,704,716,723, &
     734,742,753,760,767,777,784,792,798,802,807,819, &
     826,840,846,859,870,878,886,930,942,948,979,988, &
     1014,1021,1029,1054,1060,1067,1091,1106,1115,1125,1134,1146, &
     1156,1165,1173,1223,1233,1244,1278,1287,1298,1314,1341,1351, &
     1374,1383,1409,1456,1464,1494,1505,1523,1538,1553,1575,1598, &
     1734,1804,1885,2102,2217,2255 /)

  INTEGER, PARAMETER, DIMENSION(234) :: LU_DIAG = (/ &
       1,  2,  3,  4,  7, 11, 16, 19, 21, 23, 26, 28, &
      30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, &
      54, 56, 58, 60, 62, 64, 66, 68, 71, 73, 76, 79, &
      81, 84, 87, 91, 95, 98,101,104,107,111,117,121, &
     126,130,134,137,142,147,151,155,159,163,167,171, &
     175,179,183,187,191,195,199,203,207,211,215,219, &
     223,227,231,235,239,245,248,252,257,261,265,269, &
     273,277,281,285,289,293,297,301,305,309,313,317, &
     321,325,329,333,337,341,345,349,354,359,363,373, &
     381,386,390,394,398,402,406,410,414,418,422,427, &
     431,436,440,444,448,452,456,462,466,470,475,482, &
     490,499,504,509,513,520,526,533,539,545,551,557, &
     562,567,574,583,588,593,602,613,621,627,633,642, &
     653,660,666,671,677,681,687,695,700,709,717,726, &
     736,745,754,761,771,778,785,793,798,802,810,820, &
     831,841,851,863,872,880,907,934,943,964,982,999, &
     1015,1023,1042,1055,1062,1084,1097,1108,1118,1127,1139,1148, &
     1159,1167,1209,1226,1237,1265,1280,1291,1306,1330,1345,1365, &
     1376,1399,1443,1458,1484,1498,1517,1532,1547,1569,1592,1728, &
     1799,1881,2099,2215,2254,2255 /)


END MODULE cri_mosaic_8bin_aq_JacobianSP

