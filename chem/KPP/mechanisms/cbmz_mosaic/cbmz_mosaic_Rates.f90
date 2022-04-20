! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The Reaction Rates File
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
! File                 : cbmz_mosaic_Rates.f90
! Time                 : Thu Apr 14 13:13:32 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/cbmz_mosaic
! Equation file        : cbmz_mosaic.kpp
! Output root filename : cbmz_mosaic
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbmz_mosaic_Rates

  USE cbmz_mosaic_Parameters
  USE cbmz_mosaic_Global
  IMPLICIT NONE

CONTAINS



! Begin Rate Law Functions from KPP_HOME/util/UserRateLaws

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  User-defined Rate Law functions
!  Note: the default argument type for rate laws, as read from the equations file, is single precision
!        but all the internal calculations are performed in double precision
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~>  Arrhenius
   REAL(kind=dp) FUNCTION ARR( A0,B0,C0 )
      REAL A0,B0,C0      
      ARR =  DBLE(A0) * EXP(-DBLE(B0)/TEMP) * (TEMP/300.0_dp)**DBLE(C0)
   END FUNCTION ARR        

!~~~> Simplified Arrhenius, with two arguments
!~~~> Note: The argument B0 has a changed sign when compared to ARR
   REAL(kind=dp) FUNCTION ARR2( A0,B0 )
      REAL A0,B0           
      ARR2 =  DBLE(A0) * EXP( DBLE(B0)/TEMP )              
   END FUNCTION ARR2          

   REAL(kind=dp) FUNCTION EP2(A0,C0,A2,C2,A3,C3)
      REAL A0,C0,A2,C2,A3,C3
      REAL(dp) K0,K2,K3            
      K0 = DBLE(A0) * EXP(-DBLE(C0)/TEMP)
      K2 = DBLE(A2) * EXP(-DBLE(C2)/TEMP)
      K3 = DBLE(A3) * EXP(-DBLE(C3)/TEMP)
      K3 = K3*CFACTOR*1.0E6_dp
      EP2 = K0 + K3/(1.0_dp+K3/K2 )
   END FUNCTION EP2

   REAL(kind=dp) FUNCTION EP3(A1,C1,A2,C2) 
      REAL A1, C1, A2, C2
      REAL(dp) K1, K2      
      K1 = DBLE(A1) * EXP(-DBLE(C1)/TEMP)
      K2 = DBLE(A2) * EXP(-DBLE(C2)/TEMP)
      EP3 = K1 + K2*(1.0E6_dp*CFACTOR)
   END FUNCTION EP3 

   REAL(kind=dp) FUNCTION FALL ( A0,B0,C0,A1,B1,C1,CF)
      REAL A0,B0,C0,A1,B1,C1,CF
      REAL(dp) K0, K1     
      K0 = DBLE(A0) * EXP(-DBLE(B0)/TEMP)* (TEMP/300.0_dp)**DBLE(C0)
      K1 = DBLE(A1) * EXP(-DBLE(B1)/TEMP)* (TEMP/300.0_dp)**DBLE(C1)
      K0 = K0*CFACTOR*1.0E6_dp
      K1 = K0/K1
      FALL = (K0/(1.0_dp+K1))*   &
           DBLE(CF)**(1.0_dp/(1.0_dp+(LOG10(K1))**2))
   END FUNCTION FALL

  !---------------------------------------------------------------------------

  ELEMENTAL REAL(dp) FUNCTION k_3rd(temp,cair,k0_300K,n,kinf_300K,m,fc)

    INTRINSIC LOG10

    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL,     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL,     INTENT(IN) :: n         ! exponent for low pressure limit
    REAL,     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL,     INTENT(IN) :: m         ! exponent for high pressure limit
    REAL,     INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
    REAL                 :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    k_3rd   = k0_T/(1._dp+k_ratio)*fc**(1._dp/(1._dp+LOG10(k_ratio)**2))

  END FUNCTION k_3rd

  !---------------------------------------------------------------------------

  ELEMENTAL REAL(dp) FUNCTION k_arr (k_298,tdep,temp)
    ! Arrhenius function

    REAL,     INTENT(IN) :: k_298 ! k at T = 298.15K
    REAL,     INTENT(IN) :: tdep  ! temperature dependence
    REAL(dp), INTENT(IN) :: temp  ! temperature

    INTRINSIC EXP

    k_arr = k_298 * EXP(tdep*(1._dp/temp-3.3540E-3_dp)) ! 1/298.15=3.3540e-3

  END FUNCTION k_arr

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  End of User-defined Rate Law functions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! End Rate Law Functions from KPP_HOME/util/UserRateLaws


! Begin INLINED Rate Law Functions


!__________________________________________________
 
    REAL(KIND=dp) FUNCTION ARR3( A0,B0, TEMP )
    REAL(KIND=dp),INTENT(IN) :: TEMP
    REAL(KIND=dp),INTENT(IN):: A0,B0
    ARR3 = A0 * EXP(- B0 /TEMP )
    END FUNCTION ARR3
!__________________________________________________

    REAL(KIND=dp) FUNCTION ARR3MS( A0,B0,TEMP,C_M )
    REAL(KIND=dp), INTENT(IN) :: A0,B0      
    REAL(KIND=dp), INTENT(IN) :: TEMP,C_M  

    ARR3MS = C_M*A0 *(TEMP/300._dp)**(-B0)
    END FUNCTION ARR3MS
!__________________________________________________

    REAL(KIND=dp) FUNCTION TROEMS(k0_300K,n,kinf_300K,m,TEMP,C_M)

    INTRINSIC LOG10

    REAL(KIND=dp), INTENT(IN) :: TEMP      ! TEMPerature [K]
    REAL(KIND=dp), INTENT(IN) :: C_M      ! air concentration [molecules/cm3]
    REAL(KIND=dp), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(KIND=dp), INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(KIND=dp), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(KIND=dp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(KIND=dp)             :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = TEMP/300._dp
    k0_T    = k0_300K   * zt_help**(n) * C_M ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    TROEMS   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

    END FUNCTION TROEMS
!__________________________________________________

    REAL(KIND=dp) FUNCTION TROEEMS(A, B, k0_300K,n,kinf_300K,m,TEMP,C_M)

    INTRINSIC LOG10

    REAL(KIND=dp), INTENT(IN) :: TEMP      ! TEMPerature [K]
    REAL(KIND=dp), INTENT(IN) :: C_M      ! air concentration [molecules/cm3]
    REAL(KIND=dp), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(KIND=dp), INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(KIND=dp), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(KIND=dp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(KIND=dp), INTENT(IN) :: A, B
    REAL(KIND=dp)             :: zt_help, k0_T, kinf_T, k_ratio, troe


    zt_help = TEMP/300._dp
    k0_T    = k0_300K   * zt_help**(n) * C_M ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    troe   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

    TROEEMS = A * EXP( - B / TEMP) * troe
    END FUNCTION TROEEMS

!__________________________________________________


   REAL(KIND=dp) FUNCTION k46( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   REAL(KIND=dp) :: k0, k2, k3 
   k0=7.2E-15_dp * EXP(785._dp/TEMP)
   k2=4.1E-16_dp * EXP(1440._dp/TEMP)
   k3=1.9E-33_dp * EXP(725._dp/TEMP)
   k46=k0+k3/(1+k3/k2)
   END FUNCTION k46
!__________________________________________________
   REAL(KIND=dp) FUNCTION RK_HO_HNO3( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   REAL(KIND=dp) :: k1, k2, k3 
   k1=7.2E-15_dp * EXP(785._dp/TEMP)
   k2=1.9E-33_dp * EXP(725._dp/TEMP)
   k3=4.1E-16_dp * EXP(1440._dp/TEMP)
   RK_HO_HNO3=k1+(C_M*k2)/(1+(C_M*k2)/k3)
   END FUNCTION RK_HO_HNO3
!__________________________________________________
   REAL(KIND=dp) FUNCTION RK_2HO2( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   REAL(KIND=dp) :: k1, k2, k3 
   k1=2.3E-13_dp * EXP(600._dp/TEMP)
   k2=1.7E-33_dp * EXP(1000._dp/TEMP)
   RK_2HO2=k1+(C_M*k2)
   END FUNCTION RK_2HO2
!__________________________________________________
   REAL(KIND=dp) FUNCTION RK_2HO2_H2O( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   REAL(KIND=dp) :: k1, k2, k3 
   k1=2.3E-13_dp * EXP(600._dp/TEMP)
   k2=1.7E-33_dp * EXP(1000._dp/TEMP) * C_M
   k3=1.4E-21_dp * EXP(2200._dp/TEMP)
   RK_2HO2_H2O=(k1+k2)*k3
   END FUNCTION RK_2HO2_H2O
!__________________________________________________
   REAL(KIND=dp) FUNCTION RK_CO_HO( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   RK_CO_HO =1.5e-13 * (1.0 + 8.18e-23 * TEMP * C_M)
   END FUNCTION RK_CO_HO
!__________________________________________________
   REAL(KIND=dp) FUNCTION peroxy(K,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10, &
     TEMP,C_M)
   REAL(KIND=dp), INTENT(IN) :: X1,X2,X3,X4,X5,X6,X7,X8,X9,X10
   REAL(KIND=dp), INTENT(IN) :: TEMP,C_M
   INTEGER :: nperox, I, J, K
   PARAMETER(nperox=10)
   REAL(KIND=dp) :: Aperox(nperox,nperox),Bperox(nperox,nperox)
   REAL(KIND=dp) :: RK_PEROX(nperox,nperox)
   REAL(KIND=dp) :: RK_PARAM(nperox),SPEROX(nperox)
!
   SPEROX(1)=X1
   SPEROX(2)=X2
   SPEROX(3)=X3
   SPEROX(4)=X4
   SPEROX(5)=X5
   SPEROX(6)=X6
   SPEROX(7)=X7
   SPEROX(8)=X8
   SPEROX(9)=X9
   SPEROX(10)=X10
!
   Aperox(1,1)=2.5e-13
   Aperox(2,2)=6.8e-14
   Aperox(3,3)=2.9e-12
   Aperox(4,4)=8.0e-12
   Aperox(5,5)=1.0e-12
   Aperox(6,6)=5.3e-16
   Aperox(7,7)=3.1e-14
   Aperox(8,8)=3.1e-14
   Aperox(9,9)=3.1e-14
   Aperox(10,10)=3.1e-14
   Bperox(1,1)=190.
   Bperox(2,2)=0.0
   Bperox(3,3)=500.
   Bperox(4,4)=0.0
   Bperox(5,5)=0.0
   Bperox(6,6)=1980.
   Bperox(7,7)=1000.
   Bperox(8,8)=1000.
   Bperox(9,9)=1000.
   Bperox(10,10)=1000.
   DO I=1,nperox
   DO J=1,nperox
     IF(I.NE.J) THEN
       Aperox(I,J)=2.0*SQRT(Aperox(I,I)*Aperox(J,J))
       Bperox(I,J)=0.5*(Bperox(I,I)+Bperox(J,J))
     ENDIF
   ENDDO
   ENDDO
   Aperox(3,1)=1.3e-12
   Aperox(1,3)=1.3e-12
   Bperox(3,1)=640.0
   Bperox(1,3)=640.0
!
   DO I=1,nperox
     RK_PARAM(I)=0.0
   ENDDO
   DO I=1,nperox
   DO J=1,nperox
     RK_PEROX(I,J)=ARR3(Aperox(I,J),-Bperox(I,J),TEMP)
     RK_PARAM(I)=RK_PARAM(I)+RK_PEROX(I,J)*SPEROX(J)
   ENDDO
   ENDDO
   peroxy=RK_PARAM(K)
!
   END FUNCTION peroxy
!__________________________________________________

! End INLINED Rate Law Functions

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_SUN - update SUN light using TIME
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE Update_SUN()
      !USE cbmz_mosaic_Parameters
      !USE cbmz_mosaic_Global

    IMPLICIT NONE

    REAL(kind=dp) SunRise, SunSet
    REAL(kind=dp) Thour, Tlocal, Ttmp 
   
    SunRise = 4.5_dp 
    SunSet  = 19.5_dp 
    Thour = TIME/3600.0_dp 
    Tlocal = Thour - (INT(Thour)/24)*24

    IF ((Tlocal>=SunRise).AND.(Tlocal<=SunSet)) THEN
       Ttmp = (2.0*Tlocal-SunRise-SunSet)/(SunSet-SunRise)
       IF (Ttmp.GT.0) THEN
          Ttmp =  Ttmp*Ttmp
       ELSE
          Ttmp = -Ttmp*Ttmp
       END IF
       SUN = ( 1.0_dp + COS(PI*Ttmp) )/2.0_dp 
    ELSE
       SUN = 0.0_dp 
    END IF

 END SUBROUTINE Update_SUN

! End of Update_SUN function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_RCONST - function to update rate constants
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Update_RCONST ( )




! Begin INLINED RCONST


! End INLINED RCONST

  RCONST(1) = (j(Pj_no2))
  RCONST(2) = (j(Pj_no3o))
  RCONST(3) = (j(Pj_hno2))
  RCONST(4) = (j(Pj_hno3))
  RCONST(5) = (j(Pj_hno4))
  RCONST(6) = (j(Pj_n2o5))
  RCONST(7) = (j(Pj_o33p))
  RCONST(8) = (j(Pj_o31d))
  RCONST(9) = (j(Pj_h2o2))
  RCONST(10) = (j(Pj_ch2or))
  RCONST(11) = (j(Pj_ch2om))
  RCONST(12) = (0.7*j(Pj_h2o2))
  RCONST(13) = (0.7*j(Pj_h2o2))
  RCONST(14) = (4.6D-4*j(Pj_no2))
  RCONST(15) = (7.8D-5*j(Pj_no2))
  RCONST(16) = (9.64*j(Pj_ch2or))
  RCONST(17) = (9.04*j(Pj_ch2or))
  RCONST(18) = (0.7*j(Pj_h2o2))
  RCONST(19) = (1.0D-4*j(Pj_no2))
  RCONST(20) = (0.025*j(Pj_ch2om))
  RCONST(21) = (.79*ARR3(1.8D-11,-110.0_dp,TEMP)+.21*ARR3(3.2D-11,-70.0_dp,TEMP))
  RCONST(22) = (2.2D-10)
  RCONST(23) = (.21*ARR3MS(6.0D-34,2.3_dp,TEMP,C_M))
  RCONST(24) = (ARR3(8.0D-12,2060._dp,TEMP))
  RCONST(25) = (ARR3(6.5D-12,-120._dp,TEMP))
  RCONST(26) = (TROEMS(9.0D-32,-2.0_dp,2.2D-11,0.0_dp,TEMP,C_M))
  RCONST(27) = (TROEMS(9.0D-32,-1.5_dp,3.0D-11,0.0_dp,TEMP,C_M))
  RCONST(28) = (ARR3(2.0D-12,1400._dp,TEMP))
  RCONST(29) = (ARR3(1.2D-13,2450._dp,TEMP))
  RCONST(30) = (ARR3(1.6D-12,940._dp,TEMP))
  RCONST(31) = (ARR3(1.1D-14,500._dp,TEMP))
  RCONST(32) = (5.8D-7*ARR3(5.5D-12,2000._dp,TEMP))
  RCONST(33) = (TROEMS(7.0D-31,-2.6_dp,3.6D-11,-0.1_dp,TEMP,C_M))
  RCONST(34) = (TROEMS(2.5D-30,-4.4_dp,1.6D-11,-1.7_dp,TEMP,C_M))
  RCONST(35) = (2.2D-11)
  RCONST(36) = (ARR3(1.8D-11,390._dp,TEMP))
  RCONST(37) = (RK_HO_HNO3(TEMP,C_M))
  RCONST(38) = (ARR3(1.3D-12,-380._dp,TEMP))
  RCONST(39) = (ARR3(4.8D-11,-250._dp,TEMP))
  RCONST(40) = (ARR3(2.9D-12,160._dp,TEMP))
  RCONST(41) = (RK_2HO2(TEMP,C_M))
  RCONST(42) = (RK_2HO2_H2O(TEMP,C_M))
  RCONST(43) = (ARR3(3.5D-12,-250._dp,TEMP))
  RCONST(44) = (TROEMS(1.8D-31,-3.2_dp,4.7D-12,-1.4_dp,TEMP,C_M))
  RCONST(45) = (5.0D-16)
  RCONST(46) = (TROEEMS(4.8D+26,10900._dp,1.8D-31,-3.2_dp,4.7D-12,-1.4_dp,TEMP,C_M))
  RCONST(47) = (ARR3(1.5D-11,-170._dp,TEMP))
  RCONST(48) = (ARR3(4.5D-14,1260._dp,TEMP))
  RCONST(49) = (TROEMS(2.2D-30,-3.9_dp,1.5D-12,-0.7_dp,TEMP,C_M))
  RCONST(50) = (ARR3(8.5D-13,2450._dp,TEMP))
  RCONST(51) = (3.5D-12)
  RCONST(52) = (2.0D-21)
  RCONST(53) = (TROEEMS(3.7D+26,11000._dp,2.2D-30,-3.9_dp,1.5D-12,-0.7_dp,TEMP,C_M))
  RCONST(54) = (RK_CO_HO(TEMP,C_M))
  RCONST(55) = (TROEMS(3.0D-31,-3.3_dp,1.5D-12,0.0_dp,TEMP,C_M))
  RCONST(56) = (TEMP**0.667*ARR3(2.8D-14,1575._dp,TEMP))
  RCONST(57) = (TEMP**2*ARR3(1.5D-17,492._dp,TEMP))
  RCONST(58) = (8.1D-13)
  RCONST(59) = (1.0D-11)
  RCONST(60) = (ARR3(3.4D-13,1900._dp,TEMP))
  RCONST(61) = (ARR3(5.6D-12,-270._dp,TEMP))
  RCONST(62) = (ARR3(1.4D-12,1900._dp,TEMP))
  RCONST(63) = (TEMP**2*ARR3(5.3D-18,230._dp,TEMP))
  RCONST(64) = (1.7D-11)
  RCONST(65) = (ARR3(1.4D-12,1900._dp,TEMP))
  RCONST(66) = (ARR3(1.2D-14,2630._dp,TEMP))
  RCONST(67) = (TROEMS(1.0D-28,-0.8_dp,8.8D-12,0.0_dp,TEMP,C_M))
  RCONST(68) = (ARR3(4.2D-15,1800._dp,TEMP))
  RCONST(69) = (ARR3(8.9D-16,392._dp,TEMP))
  RCONST(70) = (ARR3(5.8D-12,-478._dp,TEMP))
  RCONST(71) = (ARR3(2.9D-11,-255._dp,TEMP))
  RCONST(72) = (ARR3(3.1D-13,1010._dp,TEMP))
  RCONST(73) = (2.5D-12)
  RCONST(74) = (ARR3(2.1D-12,-322._dp,TEMP))
  RCONST(75) = (ARR3(1.7D-11,-116._dp,TEMP))
  RCONST(76) = (8.1D-12)
  RCONST(77) = (4.1D-11)
  RCONST(78) = (2.2D-11)
  RCONST(79) = (1.4D-11)
  RCONST(80) = (3.0D-11)
  RCONST(81) = (ARR3(5.4D-17,500._dp,TEMP))
  RCONST(82) = (ARR3(2.6D-11,-409._dp,TEMP))
  RCONST(83) = (ARR3(1.2D-14,2013._dp,TEMP))
  RCONST(84) = (ARR3(3.0D-12,446._dp,TEMP))
  RCONST(85) = (3.3D-11)
  RCONST(86) = (7.0D-18)
  RCONST(87) = (1.0D-15)
  RCONST(88) = (4.0D-12)
  RCONST(89) = (4.0D-12)
  RCONST(90) = (4.0D-12)
  RCONST(91) = (ARR3(1.7D-13,-1300._dp,TEMP))
  RCONST(92) = (ARR3(1.7D-13,-1300._dp,TEMP))
  RCONST(93) = (ARR3(1.7D-13,-1300._dp,TEMP))
  RCONST(94) = (ARR3(3.8D-12,-200._dp,TEMP))
  RCONST(95) = (ARR3(3.8D-12,-200._dp,TEMP))
  RCONST(96) = (ARR3(3.8D-12,-200._dp,TEMP))
  RCONST(97) = (ARR3(1.6D-11,540._dp,TEMP))
  RCONST(98) = (TROEMS(9.7D-29,-5.6_dp,9.3D-12,-1.5_dp,TEMP,C_M))
  RCONST(99) = (TROEEMS(1.1D+28,14000._dp,9.7D-29,-5.6_dp,9.3D-12,-1.5_dp,TEMP,C_M))
  RCONST(100) = (ARR3(6.7D-12,600._dp,TEMP))
  RCONST(101) = (ARR3(3.0D-12,-280._dp,TEMP))
  RCONST(102) = (ARR3(2.6D-12,-365._dp,TEMP))
  RCONST(103) = (4.0D-12)
  RCONST(104) = (ARR3(5.3D-12,-360._dp,TEMP))
  RCONST(105) = (4.0D-12)
  RCONST(106) = (4.0D-12)
  RCONST(107) = (4.0D-12)
  RCONST(108) = (1.1D-12)
  RCONST(109) = (2.5D-12)
  RCONST(110) = (2.5D-12)
  RCONST(111) = (4.0D-12)
  RCONST(112) = (1.2D-12)
  RCONST(113) = (4.0D-12)
  RCONST(114) = (2.5D-12)
  RCONST(115) = (ARR3(3.8D-13,-800._dp,TEMP))
  RCONST(116) = (ARR3(7.5D-13,-700._dp,TEMP))
  RCONST(117) = (ARR3(1.7D-13,-1300._dp,TEMP))
  RCONST(118) = (ARR3(4.5D-13,-1000._dp,TEMP))
  RCONST(119) = (ARR3(1.2D-13,-1300._dp,TEMP))
  RCONST(120) = (ARR3(1.7D-13,-1300._dp,TEMP))
  RCONST(121) = (ARR3(1.7D-13,-1300._dp,TEMP))
  RCONST(122) = (ARR3(7.0D-12,235._dp,TEMP))
  RCONST(123) = (1.0D-11)
! RCONST(124) = constant rate coefficient
! RCONST(125) = constant rate coefficient
! RCONST(126) = constant rate coefficient
  RCONST(127) = (peroxy(1,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(128) = (peroxy(2,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(129) = (peroxy(3,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(130) = (peroxy(4,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(131) = (peroxy(5,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(132) = (peroxy(6,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(133) = (peroxy(7,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(134) = (peroxy(8,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(135) = (peroxy(9,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(136) = (peroxy(10,CH3O2,ETHP,RO2,C2O3,ANO2,NAP,ISOPP,ISOPN,ISOPO2,XO2,TEMP,C_M))
  RCONST(137) = (ARR3(12.1D-12,-444._dp,TEMP))
  RCONST(138) = (ARR3(1.01D-15,732._dp,TEMP))
  RCONST(139) = (ARR3(1.19D-12,-490._dp,TEMP))
  RCONST(140) = (1.71D-10)
  RCONST(141) = (2.0D-16)
  RCONST(142) = (1.22D-11)
      
END SUBROUTINE Update_RCONST

! End of Update_RCONST function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_PHOTO - function to update photolytical rate constants
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Update_PHOTO ( )


   USE cbmz_mosaic_Global

  RCONST(1) = (j(Pj_no2))
  RCONST(2) = (j(Pj_no3o))
  RCONST(3) = (j(Pj_hno2))
  RCONST(4) = (j(Pj_hno3))
  RCONST(5) = (j(Pj_hno4))
  RCONST(6) = (j(Pj_n2o5))
  RCONST(7) = (j(Pj_o33p))
  RCONST(8) = (j(Pj_o31d))
  RCONST(9) = (j(Pj_h2o2))
  RCONST(10) = (j(Pj_ch2or))
  RCONST(11) = (j(Pj_ch2om))
  RCONST(12) = (0.7*j(Pj_h2o2))
  RCONST(13) = (0.7*j(Pj_h2o2))
  RCONST(14) = (4.6D-4*j(Pj_no2))
  RCONST(15) = (7.8D-5*j(Pj_no2))
  RCONST(16) = (9.64*j(Pj_ch2or))
  RCONST(17) = (9.04*j(Pj_ch2or))
  RCONST(18) = (0.7*j(Pj_h2o2))
  RCONST(19) = (1.0D-4*j(Pj_no2))
  RCONST(20) = (0.025*j(Pj_ch2om))
      
END SUBROUTINE Update_PHOTO

! End of Update_PHOTO function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE cbmz_mosaic_Rates

