























MODULE cbmz_bb_UpdateRconstWRF

  USE cbmz_bb_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE cbmz_bb_Update_RCONST(  &

ch3o2,ethp,ro2,c2o3,ano2,nap,isopp,isopn,isopo2,xo2, &


                j, nj,   &
                RCONST, &
             Pj_o31d, Pj_o33p, Pj_no2, Pj_no3o2, Pj_no3o, & 
             Pj_hno2, Pj_hno3, Pj_hno4, Pj_h2o2, Pj_ch2or, & 
             Pj_ch2om, Pj_ch3cho, Pj_ch3coch3, Pj_ch3coc2h5, Pj_hcocho, & 
             Pj_ch3cocho, Pj_hcochest, Pj_ch3o2h, Pj_ch3coo2h, Pj_ch3ono2, & 
             Pj_hcochob, Pj_macr, Pj_n2o5, Pj_o2, Pj_pan, & 
             Pj_acet, Pj_mglo, Pj_hno4_2, Pj_clno2, Pj_n2o, & 
             Pj_pooh, Pj_mpan, Pj_mvk, Pj_etooh, Pj_prooh, & 
             Pj_onitr, Pj_acetol, Pj_glyald, Pj_hyac, Pj_mek, & 
             Pj_open, Pj_gly, Pj_acetp, Pj_xooh, Pj_isooh, & 
             Pj_alkooh, Pj_mekooh, Pj_tolooh, Pj_terpooh, Pj_cl2, & 
             Pj_hocl, Pj_fmcl,  & 
                C_M, C_H2O, TEMP & 

)




   IMPLICIT NONE

  INTEGER, INTENT (IN ) :: nj 

    REAL(KIND=dp), DIMENSION(nj), INTENT(IN)  :: j


    REAL(KIND=dp), DIMENSION(NREACT), INTENT(OUT)  :: RCONST


    REAL(KIND=dp), INTENT(IN)  :: C_M, C_H2O,&
                                  TEMP





        INTEGER, INTENT(IN ) ::  & 
               Pj_o31d, Pj_o33p, Pj_no2, Pj_no3o2, Pj_no3o, & 
               Pj_hno2, Pj_hno3, Pj_hno4, Pj_h2o2, Pj_ch2or, & 
               Pj_ch2om, Pj_ch3cho, Pj_ch3coch3, Pj_ch3coc2h5, Pj_hcocho, & 
               Pj_ch3cocho, Pj_hcochest, Pj_ch3o2h, Pj_ch3coo2h, Pj_ch3ono2, & 
               Pj_hcochob, Pj_macr, Pj_n2o5, Pj_o2, Pj_pan, & 
               Pj_acet, Pj_mglo, Pj_hno4_2, Pj_clno2, Pj_n2o, & 
               Pj_pooh, Pj_mpan, Pj_mvk, Pj_etooh, Pj_prooh, & 
               Pj_onitr, Pj_acetol, Pj_glyald, Pj_hyac, Pj_mek, & 
               Pj_open, Pj_gly, Pj_acetp, Pj_xooh, Pj_isooh, & 
               Pj_alkooh, Pj_mekooh, Pj_tolooh, Pj_terpooh, Pj_cl2, & 
               Pj_hocl, Pj_fmcl

REAL(KIND=dp), INTENT(IN) ::  ch3o2,ethp,ro2,c2o3,ano2,nap,isopp,isopn,isopo2,xo2









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
  RCONST(124) = 1
  RCONST(125) = 1
  RCONST(126) = 1
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
END SUBROUTINE cbmz_bb_Update_RCONST













    REAL(kind=dp) FUNCTION ARR( A0,B0,C0, TEMP )
      REAL(kind=dp) :: TEMP
      REAL(kind=dp) A0,B0,C0
      ARR =  A0 * EXP( -B0 /TEMP ) * (TEMP/300._dp)**C0
    END FUNCTION ARR


   REAL(kind=dp) FUNCTION ARR2( A0,B0, TEMP )
      REAL(kind=dp) :: TEMP 
      REAL(kind=dp) A0,B0           
      ARR2 = A0 * EXP( -B0 /TEMP )              
   END FUNCTION ARR2          


   REAL(kind=dp) FUNCTION EP2( A0,C0,A2,C2,A3,C3,TEMP,cair)
      REAL(kind=dp) :: TEMP
      REAL(kind=dp) :: cair
      REAL(kind=dp) A0,C0,A2,C2,A3,C3
      REAL(kind=dp) K0,K2,K3

      K0 = A0 * EXP(-C0 /TEMP)
      K2 = A2 * EXP(-C2 /TEMP)
      K3 = A3 * EXP(-C3 /TEMP)

      K3 = K3 * cair
      EP2 = K0 + K3/(1._dp+K3/K2 )
   END FUNCTION EP2


    REAL(kind=dp) FUNCTION EP3(A1,C1,A2,C2,TEMP,cair)
      REAL(kind=dp) :: TEMP
      REAL(kind=dp) :: cair
      REAL(kind=dp) A1, C1, A2, C2
      REAL(kind=dp) K1, K2
 
      K1 = A1 * EXP(-C1 /TEMP)
      K2 = A2 * EXP(-C2 /TEMP)

      EP3 = K1 + K2*cair
    END FUNCTION EP3


    REAL(kind=dp) FUNCTION FALL( A0,B0,C0,A1,B1,C1,CF,TEMP,cair)
 
      INTRINSIC LOG10  

      REAL(kind=dp) :: TEMP
      REAL(kind=dp) :: cair
      REAL(kind=dp) A0,B0,C0,A1,B1,C1,CF
      REAL(kind=dp) K0, K1

      K0 = A0 * EXP(-B0 /TEMP)* (TEMP/300._dp)**C0
      K1 = A1 * EXP(-B1 /TEMP)* (TEMP/300._dp)**C1

      K0 = K0 * cair
      K1 = K0/K1
      FALL = (K0/(1._dp+K1))*CF**(1._dp/(1._dp+(LOG10(K1))**2))

    END FUNCTION FALL


    REAL(kind=dp) FUNCTION F2( A0,B0,C0,A1,B1,C1,CF,CN,TEMP,cair)

      INTRINSIC LOG10

      REAL(kind=dp) :: TEMP
      REAL(kind=dp) :: cair
      REAL(kind=dp) A0,B0,C0,A1,B1,C1,CF,CN
      REAL(kind=dp) K0, K1

      K0 = A0 * EXP(-B0 /TEMP)* (TEMP/300._dp)**C0
      K1 = A1 * EXP(-B1 /TEMP)* (TEMP/300._dp)**C1

      K0 = K0 * cair
      K1 = K0/K1
      F2 = (K0/(1._dp+K1))*CF**(1._dp/(1._dp+(LOG10(K1)/CN)**2))

    END FUNCTION F2
                                                                   



    REAL(kind=dp) FUNCTION TROE(k0_300K,n,kinf_300K,m,temp,cair)

    INTRINSIC LOG10

    REAL(kind=dp), INTENT(IN) :: temp      
    REAL(kind=dp), INTENT(IN) :: cair      
    REAL(kind=dp),          INTENT(IN) :: k0_300K   
    REAL(kind=dp),          INTENT(IN) :: n         
    REAL(kind=dp),          INTENT(IN) :: kinf_300K 
    REAL(kind=dp),          INTENT(IN) :: m         
    REAL(kind=dp)             :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair 
    kinf_T  = kinf_300K * zt_help**(m)        
    k_ratio = k0_T/kinf_T
    TROE   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

   END FUNCTION TROE






    REAL(kind=dp) FUNCTION TROEE(A, B, k0_300K,n,kinf_300K,m,temp,cair)

    INTRINSIC LOG10

    REAL(kind=dp), INTENT(IN) :: temp      
    REAL(kind=dp), INTENT(IN) :: cair      
    REAL(kind=dp),     INTENT(IN) :: k0_300K   
    REAL(kind=dp),     INTENT(IN) :: n         
    REAL(kind=dp),     INTENT(IN) :: kinf_300K 
    REAL(kind=dp),     INTENT(IN) :: m         
    REAL(kind=dp),     INTENT(IN) :: A, B 
    REAL(kind=dp)             :: zt_help, k0_T, kinf_T, k_ratio, troe
    

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair 
    kinf_T  = kinf_300K * zt_help**(m)        
    k_ratio = k0_T/kinf_T
    troe   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

    TROEE = A * EXP( - B / temp) * troe
    
    

  END FUNCTION TROEE




   REAL(kind=dp) FUNCTION THERMAL_T2(c, d ,temp)
    REAL(kind=dp), INTENT(IN) :: temp      
    REAL(kind=dp),     INTENT(IN) :: c, d


     THERMAL_T2= temp**2._dp * c * EXP(- d  / temp)

   END FUNCTION THERMAL_T2












 
    REAL(KIND=dp) FUNCTION ARR3( A0,B0, TEMP )
    REAL(KIND=dp),INTENT(IN) :: TEMP
    REAL(KIND=dp),INTENT(IN):: A0,B0
    ARR3 = A0 * EXP(- B0 /TEMP )
    END FUNCTION ARR3


    REAL(KIND=dp) FUNCTION ARR3MS( A0,B0,TEMP,C_M )
    REAL(KIND=dp), INTENT(IN) :: A0,B0      
    REAL(KIND=dp), INTENT(IN) :: TEMP,C_M  

    ARR3MS = C_M*A0 *(TEMP/300._dp)**(-B0)
    END FUNCTION ARR3MS


    REAL(KIND=dp) FUNCTION TROEMS(k0_300K,n,kinf_300K,m,TEMP,C_M)

    INTRINSIC LOG10

    REAL(KIND=dp), INTENT(IN) :: TEMP      
    REAL(KIND=dp), INTENT(IN) :: C_M      
    REAL(KIND=dp), INTENT(IN) :: k0_300K   
    REAL(KIND=dp), INTENT(IN) :: n         
    REAL(KIND=dp), INTENT(IN) :: kinf_300K 
    REAL(KIND=dp), INTENT(IN) :: m         
    REAL(KIND=dp)             :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = TEMP/300._dp
    k0_T    = k0_300K   * zt_help**(n) * C_M 
    kinf_T  = kinf_300K * zt_help**(m)        
    k_ratio = k0_T/kinf_T
    TROEMS   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

    END FUNCTION TROEMS


    REAL(KIND=dp) FUNCTION TROEEMS(A, B, k0_300K,n,kinf_300K,m,TEMP,C_M)

    INTRINSIC LOG10

    REAL(KIND=dp), INTENT(IN) :: TEMP      
    REAL(KIND=dp), INTENT(IN) :: C_M      
    REAL(KIND=dp), INTENT(IN) :: k0_300K   
    REAL(KIND=dp), INTENT(IN) :: n         
    REAL(KIND=dp), INTENT(IN) :: kinf_300K 
    REAL(KIND=dp), INTENT(IN) :: m         
    REAL(KIND=dp), INTENT(IN) :: A, B
    REAL(KIND=dp)             :: zt_help, k0_T, kinf_T, k_ratio, troe


    zt_help = TEMP/300._dp
    k0_T    = k0_300K   * zt_help**(n) * C_M 
    kinf_T  = kinf_300K * zt_help**(m)        
    k_ratio = k0_T/kinf_T
    troe   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

    TROEEMS = A * EXP( - B / TEMP) * troe
    END FUNCTION TROEEMS




   REAL(KIND=dp) FUNCTION k46( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   REAL(KIND=dp) :: k0, k2, k3 
   k0=7.2E-15_dp * EXP(785._dp/TEMP)
   k2=4.1E-16_dp * EXP(1440._dp/TEMP)
   k3=1.9E-33_dp * EXP(725._dp/TEMP)
   k46=k0+k3/(1+k3/k2)
   END FUNCTION k46

   REAL(KIND=dp) FUNCTION RK_HO_HNO3( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   REAL(KIND=dp) :: k1, k2, k3 
   k1=7.2E-15_dp * EXP(785._dp/TEMP)
   k2=1.9E-33_dp * EXP(725._dp/TEMP)
   k3=4.1E-16_dp * EXP(1440._dp/TEMP)
   RK_HO_HNO3=k1+(C_M*k2)/(1+(C_M*k2)/k3)
   END FUNCTION RK_HO_HNO3

   REAL(KIND=dp) FUNCTION RK_2HO2( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   REAL(KIND=dp) :: k1, k2, k3 
   k1=2.3E-13_dp * EXP(600._dp/TEMP)
   k2=1.7E-33_dp * EXP(1000._dp/TEMP)
   RK_2HO2=k1+(C_M*k2)
   END FUNCTION RK_2HO2

   REAL(KIND=dp) FUNCTION RK_2HO2_H2O( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   REAL(KIND=dp) :: k1, k2, k3 
   k1=2.3E-13_dp * EXP(600._dp/TEMP)
   k2=1.7E-33_dp * EXP(1000._dp/TEMP) * C_M
   k3=1.4E-21_dp * EXP(2200._dp/TEMP)
   RK_2HO2_H2O=(k1+k2)*k3
   END FUNCTION RK_2HO2_H2O

   REAL(KIND=dp) FUNCTION RK_CO_HO( TEMP, C_M )
   REAL(KIND=dp), INTENT(IN) :: TEMP, C_M 
   RK_CO_HO =1.5e-13 * (1.0 + 8.18e-23 * TEMP * C_M)
   END FUNCTION RK_CO_HO

   REAL(KIND=dp) FUNCTION peroxy(K,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10, &
     TEMP,C_M)
   REAL(KIND=dp), INTENT(IN) :: X1,X2,X3,X4,X5,X6,X7,X8,X9,X10
   REAL(KIND=dp), INTENT(IN) :: TEMP,C_M
   INTEGER :: nperox, I, J, K
   PARAMETER(nperox=10)
   REAL(KIND=dp) :: Aperox(nperox,nperox),Bperox(nperox,nperox)
   REAL(KIND=dp) :: RK_PEROX(nperox,nperox)
   REAL(KIND=dp) :: RK_PARAM(nperox),SPEROX(nperox)

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

   END FUNCTION peroxy





END MODULE cbmz_bb_UpdateRconstWRF

