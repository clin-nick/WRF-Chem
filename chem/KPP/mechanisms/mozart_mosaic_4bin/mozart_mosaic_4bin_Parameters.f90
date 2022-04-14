! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
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
! File                 : mozart_mosaic_4bin_Parameters.f90
! Time                 : Tue Apr 12 23:43:06 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/mozart_mosaic_4bin
! Equation file        : mozart_mosaic_4bin.kpp
! Output root filename : mozart_mosaic_4bin
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mozart_mosaic_4bin_Parameters

  USE mozart_mosaic_4bin_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 138 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 136 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 129 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 2 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 292 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 137 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 1303 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 1469 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 137 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 0 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 
! PI - Value of pi
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_HONO = 1 
  INTEGER, PARAMETER :: ind_SO4 = 2 
  INTEGER, PARAMETER :: ind_NUME = 3 
  INTEGER, PARAMETER :: ind_DEN = 4 
  INTEGER, PARAMETER :: ind_BIOG1_c = 5 
  INTEGER, PARAMETER :: ind_BIOG1_o = 6 
  INTEGER, PARAMETER :: ind_SMPA = 7 
  INTEGER, PARAMETER :: ind_VOCA = 8 
  INTEGER, PARAMETER :: ind_SMPBB = 9 
  INTEGER, PARAMETER :: ind_VOCBB = 10 
  INTEGER, PARAMETER :: ind_NH3 = 11 
  INTEGER, PARAMETER :: ind_C2H6 = 12 
  INTEGER, PARAMETER :: ind_C3H8 = 13 
  INTEGER, PARAMETER :: ind_BIGENE = 14 
  INTEGER, PARAMETER :: ind_N2O = 15 
  INTEGER, PARAMETER :: ind_SO2 = 16 
  INTEGER, PARAMETER :: ind_PBZNIT = 17 
  INTEGER, PARAMETER :: ind_H2O2 = 18 
  INTEGER, PARAMETER :: ind_C2H2 = 19 
  INTEGER, PARAMETER :: ind_BENZENE = 20 
  INTEGER, PARAMETER :: ind_BEPOMUC = 21 
  INTEGER, PARAMETER :: ind_PHENOL = 22 
  INTEGER, PARAMETER :: ind_EO = 23 
  INTEGER, PARAMETER :: ind_TOLUENE = 24 
  INTEGER, PARAMETER :: ind_CRESOL = 25 
  INTEGER, PARAMETER :: ind_XOOH = 26 
  INTEGER, PARAMETER :: ind_PHENOOH = 27 
  INTEGER, PARAMETER :: ind_C6H5OOH = 28 
  INTEGER, PARAMETER :: ind_BENZOOH = 29 
  INTEGER, PARAMETER :: ind_BIGALD2 = 30 
  INTEGER, PARAMETER :: ind_TEPOMUC = 31 
  INTEGER, PARAMETER :: ind_BZOOH = 32 
  INTEGER, PARAMETER :: ind_BZALD = 33 
  INTEGER, PARAMETER :: ind_XYLOLOOH = 34 
  INTEGER, PARAMETER :: ind_XYLENOOH = 35 
  INTEGER, PARAMETER :: ind_N2O5 = 36 
  INTEGER, PARAMETER :: ind_XYLENES = 37 
  INTEGER, PARAMETER :: ind_XYLOL = 38 
  INTEGER, PARAMETER :: ind_DMS = 39 
  INTEGER, PARAMETER :: ind_BIGALD4 = 40 
  INTEGER, PARAMETER :: ind_C2H5OH = 41 
  INTEGER, PARAMETER :: ind_BIGALD3 = 42 
  INTEGER, PARAMETER :: ind_BIGALD1 = 43 
  INTEGER, PARAMETER :: ind_H2 = 44 
  INTEGER, PARAMETER :: ind_C2H5OOH = 45 
  INTEGER, PARAMETER :: ind_C3H7OOH = 46 
  INTEGER, PARAMETER :: ind_ROOH = 47 
  INTEGER, PARAMETER :: ind_ENEO2 = 48 
  INTEGER, PARAMETER :: ind_MEKOOH = 49 
  INTEGER, PARAMETER :: ind_CH3OOH = 50 
  INTEGER, PARAMETER :: ind_MACROOH = 51 
  INTEGER, PARAMETER :: ind_HO2NO2 = 52 
  INTEGER, PARAMETER :: ind_HYDRALD = 53 
  INTEGER, PARAMETER :: ind_BIGALK = 54 
  INTEGER, PARAMETER :: ind_C2H4 = 55 
  INTEGER, PARAMETER :: ind_HOCH2OO = 56 
  INTEGER, PARAMETER :: ind_TOLOOH = 57 
  INTEGER, PARAMETER :: ind_PHENO2 = 58 
  INTEGER, PARAMETER :: ind_PHENO = 59 
  INTEGER, PARAMETER :: ind_EO2 = 60 
  INTEGER, PARAMETER :: ind_BZOO = 61 
  INTEGER, PARAMETER :: ind_MEK = 62 
  INTEGER, PARAMETER :: ind_CH3COOH = 63 
  INTEGER, PARAMETER :: ind_CH3COOOH = 64 
  INTEGER, PARAMETER :: ind_TERPOOH = 65 
  INTEGER, PARAMETER :: ind_BENZO2 = 66 
  INTEGER, PARAMETER :: ind_ONIT = 67 
  INTEGER, PARAMETER :: ind_PAN = 68 
  INTEGER, PARAMETER :: ind_XYLOLO2 = 69 
  INTEGER, PARAMETER :: ind_POOH = 70 
  INTEGER, PARAMETER :: ind_HNO3 = 71 
  INTEGER, PARAMETER :: ind_CH4 = 72 
  INTEGER, PARAMETER :: ind_HMPROPO2 = 73 
  INTEGER, PARAMETER :: ind_ACBZO2 = 74 
  INTEGER, PARAMETER :: ind_ISOPOOH = 75 
  INTEGER, PARAMETER :: ind_MBOOOH = 76 
  INTEGER, PARAMETER :: ind_MPAN = 77 
  INTEGER, PARAMETER :: ind_HCOOH = 78 
  INTEGER, PARAMETER :: ind_TERP2OOH = 79 
  INTEGER, PARAMETER :: ind_O1D_CB4 = 80 
  INTEGER, PARAMETER :: ind_C6H5O2 = 81 
  INTEGER, PARAMETER :: ind_APIN = 82 
  INTEGER, PARAMETER :: ind_BPIN = 83 
  INTEGER, PARAMETER :: ind_ALKOOH = 84 
  INTEGER, PARAMETER :: ind_MEKO2 = 85 
  INTEGER, PARAMETER :: ind_PO2 = 86 
  INTEGER, PARAMETER :: ind_MALO2 = 87 
  INTEGER, PARAMETER :: ind_TOLO2 = 88 
  INTEGER, PARAMETER :: ind_MBO = 89 
  INTEGER, PARAMETER :: ind_XYLENO2 = 90 
  INTEGER, PARAMETER :: ind_DICARBO2 = 91 
  INTEGER, PARAMETER :: ind_O = 92 
  INTEGER, PARAMETER :: ind_C3H7O2 = 93 
  INTEGER, PARAMETER :: ind_MDIALO2 = 94 
  INTEGER, PARAMETER :: ind_CH3OH = 95 
  INTEGER, PARAMETER :: ind_GLYOXAL = 96 
  INTEGER, PARAMETER :: ind_C2H5O2 = 97 
  INTEGER, PARAMETER :: ind_BIGALD = 98 
  INTEGER, PARAMETER :: ind_ISOPNO3 = 99 
  INTEGER, PARAMETER :: ind_BCARY = 100 
  INTEGER, PARAMETER :: ind_LIMON = 101 
  INTEGER, PARAMETER :: ind_MBONO3O2 = 102 
  INTEGER, PARAMETER :: ind_HMPROP = 103 
  INTEGER, PARAMETER :: ind_MYRC = 104 
  INTEGER, PARAMETER :: ind_CH3COCH3 = 105 
  INTEGER, PARAMETER :: ind_ISOP = 106 
  INTEGER, PARAMETER :: ind_CO = 107 
  INTEGER, PARAMETER :: ind_GLYALD = 108 
  INTEGER, PARAMETER :: ind_C3H6 = 109 
  INTEGER, PARAMETER :: ind_TERPROD1 = 110 
  INTEGER, PARAMETER :: ind_MBOO2 = 111 
  INTEGER, PARAMETER :: ind_CH3CHO = 112 
  INTEGER, PARAMETER :: ind_ALKO2 = 113 
  INTEGER, PARAMETER :: ind_TERPROD2 = 114 
  INTEGER, PARAMETER :: ind_MACR = 115 
  INTEGER, PARAMETER :: ind_ONITR = 116 
  INTEGER, PARAMETER :: ind_CH3COCHO = 117 
  INTEGER, PARAMETER :: ind_TERPO2 = 118 
  INTEGER, PARAMETER :: ind_NTERPO2 = 119 
  INTEGER, PARAMETER :: ind_TERP2O2 = 120 
  INTEGER, PARAMETER :: ind_RO2 = 121 
  INTEGER, PARAMETER :: ind_HYAC = 122 
  INTEGER, PARAMETER :: ind_XO2 = 123 
  INTEGER, PARAMETER :: ind_MVK = 124 
  INTEGER, PARAMETER :: ind_MACRO2 = 125 
  INTEGER, PARAMETER :: ind_CH2O = 126 
  INTEGER, PARAMETER :: ind_ISOPO2 = 127 
  INTEGER, PARAMETER :: ind_CH3CO3 = 128 
  INTEGER, PARAMETER :: ind_CH3O2 = 129 
  INTEGER, PARAMETER :: ind_MCO3 = 130 
  INTEGER, PARAMETER :: ind_HO2 = 131 
  INTEGER, PARAMETER :: ind_O3 = 132 
  INTEGER, PARAMETER :: ind_NO = 133 
  INTEGER, PARAMETER :: ind_NO2 = 134 
  INTEGER, PARAMETER :: ind_NO3 = 135 
  INTEGER, PARAMETER :: ind_OH = 136 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_H2O = 137 
  INTEGER, PARAMETER :: ind_M = 138 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_M = 2 

END MODULE mozart_mosaic_4bin_Parameters

