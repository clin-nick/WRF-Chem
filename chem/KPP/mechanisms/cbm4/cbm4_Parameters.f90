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
! File                 : cbm4_Parameters.f90
! Time                 : Tue Apr 12 23:42:45 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/cbm4
! Equation file        : cbm4.kpp
! Output root filename : cbm4
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm4_Parameters

  USE cbm4_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 33 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 32 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 32 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 1 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 81 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 33 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 276 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 300 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 33 
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

  INTEGER, PARAMETER :: ind_O1D_CB4 = 1 
  INTEGER, PARAMETER :: ind_H2O2 = 2 
  INTEGER, PARAMETER :: ind_PAN = 3 
  INTEGER, PARAMETER :: ind_CRO = 4 
  INTEGER, PARAMETER :: ind_TOL = 5 
  INTEGER, PARAMETER :: ind_N2O5 = 6 
  INTEGER, PARAMETER :: ind_XYL = 7 
  INTEGER, PARAMETER :: ind_XO2N = 8 
  INTEGER, PARAMETER :: ind_HONO = 9 
  INTEGER, PARAMETER :: ind_PNA = 10 
  INTEGER, PARAMETER :: ind_TO2 = 11 
  INTEGER, PARAMETER :: ind_HNO3 = 12 
  INTEGER, PARAMETER :: ind_ROR = 13 
  INTEGER, PARAMETER :: ind_CRES = 14 
  INTEGER, PARAMETER :: ind_MGLY = 15 
  INTEGER, PARAMETER :: ind_CO = 16 
  INTEGER, PARAMETER :: ind_ETH = 17 
  INTEGER, PARAMETER :: ind_XO2 = 18 
  INTEGER, PARAMETER :: ind_OPEN = 19 
  INTEGER, PARAMETER :: ind_PAR = 20 
  INTEGER, PARAMETER :: ind_HCHO = 21 
  INTEGER, PARAMETER :: ind_ISOP = 22 
  INTEGER, PARAMETER :: ind_OLE = 23 
  INTEGER, PARAMETER :: ind_ALD2 = 24 
  INTEGER, PARAMETER :: ind_O3 = 25 
  INTEGER, PARAMETER :: ind_NO2 = 26 
  INTEGER, PARAMETER :: ind_HO = 27 
  INTEGER, PARAMETER :: ind_HO2 = 28 
  INTEGER, PARAMETER :: ind_O = 29 
  INTEGER, PARAMETER :: ind_NO3 = 30 
  INTEGER, PARAMETER :: ind_NO = 31 
  INTEGER, PARAMETER :: ind_C2O3 = 32 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_H2O = 33 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_H2O = 1 

END MODULE cbm4_Parameters

