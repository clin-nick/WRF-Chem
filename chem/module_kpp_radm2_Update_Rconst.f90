























MODULE radm2_UpdateRconstWRF

  USE radm2_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE radm2_Update_RCONST(  &


rc_n2o5, &


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


 REAL(KIND=dp) :: rc_n2o5








  RCONST(1) = (j(Pj_no2))
  RCONST(2) = (j(Pj_o31d))
  RCONST(3) = (j(Pj_o33p))
  RCONST(4) = (j(Pj_hno2))
  RCONST(5) = (j(Pj_hno3))
  RCONST(6) = (j(Pj_hno4))
  RCONST(7) = (j(Pj_no3o2))
  RCONST(8) = (j(Pj_no3o))
  RCONST(9) = (j(Pj_h2o2))
  RCONST(10) = (j(Pj_ch2om))
  RCONST(11) = (j(Pj_ch2or))
  RCONST(12) = (j(Pj_ch3cho))
  RCONST(13) = (j(Pj_ch3o2h))
  RCONST(14) = (j(Pj_ch3coch3))
  RCONST(15) = (j(Pj_ch3coo2h))
  RCONST(16) = (j(Pj_ch3coc2h5))
  RCONST(17) = (j(Pj_hcocho))
  RCONST(18) = (j(Pj_hcochob))
  RCONST(19) = (j(Pj_ch3cocho))
  RCONST(20) = (j(Pj_hcochest))
  RCONST(21) = (j(Pj_ch3ono2))
  RCONST(22) = (.20946e0*(C_M*6.00D-34*(TEMP/300.0)**(-2.3)))
  RCONST(23) = (ARR2(6.5D-12,-120.0_dp,TEMP))
  RCONST(24) = (.78084*ARR2(1.8D-11,-110.0_dp,TEMP)+.20946e0*ARR2(3.2D-11,-70.0_dp,TEMP))
  RCONST(25) = (2.2D-10)
  RCONST(26) = (ARR2(2.0D-12,1400.0_dp,TEMP))
  RCONST(27) = (ARR2(1.6D-12,940.0_dp,TEMP))
  RCONST(28) = (ARR2(1.1D-14,500.0_dp,TEMP))
  RCONST(29) = (ARR2(3.7D-12,-240.0_dp,TEMP))
  RCONST(30) = (TROE(1.80D-31,3.2_dp,4.70D-12,1.4_dp,TEMP,C_M))
  RCONST(31) = (TROEE(4.76D26,10900.0_dp,1.80D-31,3.2_dp,4.70D-12,1.4_dp,TEMP,C_M))
  RCONST(32) = ((2.2D-13*EXP(600./TEMP)+1.9D-33*C_M*EXP(980._dp/TEMP)))
  RCONST(33) = ((3.08D-34*EXP(2800._dp/TEMP)+2.66D-54*C_M*EXP(3180._dp/TEMP)))
  RCONST(34) = (ARR2(3.3D-12,200.0_dp,TEMP))
  RCONST(35) = (TROE(7.00D-31,2.6_dp,1.50D-11,0.5_dp,TEMP,C_M))
  RCONST(36) = (.20946e0*ARR2(3.3D-39,-530.0_dp,TEMP))
  RCONST(37) = (ARR2(1.4D-13,2500.0_dp,TEMP))
  RCONST(38) = (ARR2(1.7D-11,-150.0_dp,TEMP))
  RCONST(39) = (ARR2(2.5D-14,1230.0_dp,TEMP))
  RCONST(40) = (2.5D-12)
  RCONST(41) = (TROE(2.20D-30,4.3_dp,1.50D-12,0.5_dp,TEMP,C_M))
  RCONST(42) = (TROEE(9.09D26,11200.0_dp,2.20D-30,4.3_dp,1.50D-12,0.5_dp,TEMP,C_M))
  RCONST(43) = (rc_n2o5)
  RCONST(44) = (TROE(2.60D-30,3.2_dp,2.40D-11,1.3_dp,TEMP,C_M))
  RCONST(45) = (k46(TEMP,C_M))
  RCONST(46) = (ARR2(1.3D-12,-380.0_dp,TEMP))
  RCONST(47) = (ARR2(4.6D-11,-230.0_dp,TEMP))
  RCONST(48) = (TROE(3.00D-31,3.3_dp,1.50D-12,0.0_dp,TEMP,C_M))
  RCONST(49) = ((1.5D-13*(1._dp+2.439D-20*C_M)))
  RCONST(50) = (THERMAL_T2(6.95D-18,1280.0_dp,TEMP))
  RCONST(51) = (THERMAL_T2(1.37D-17,444.0_dp,TEMP))
  RCONST(52) = (ARR2(1.59D-11,540.0_dp,TEMP))
  RCONST(53) = (ARR2(1.73D-11,380.0_dp,TEMP))
  RCONST(54) = (ARR2(3.64D-11,380.0_dp,TEMP))
  RCONST(55) = (ARR2(2.15D-12,-411.0_dp,TEMP))
  RCONST(56) = (ARR2(5.32D-12,-504.0_dp,TEMP))
  RCONST(57) = (ARR2(1.07D-11,-549.0_dp,TEMP))
  RCONST(58) = (ARR2(2.1D-12,-322.0_dp,TEMP))
  RCONST(59) = (ARR2(1.89D-11,-116.0_dp,TEMP))
  RCONST(60) = (4.0D-11)
  RCONST(61) = (9.0D-12)
  RCONST(62) = (ARR2(6.87D-12,-256.0_dp,TEMP))
  RCONST(63) = (ARR2(1.2D-11,745.0_dp,TEMP))
  RCONST(64) = (1.15D-11)
  RCONST(65) = (1.7D-11)
  RCONST(66) = (2.8D-11)
  RCONST(67) = (1.0D-11)
  RCONST(68) = (1.0D-11)
  RCONST(69) = (1.0D-11)
  RCONST(70) = (THERMAL_T2(6.85D-18,444.0_dp,TEMP))
  RCONST(71) = (ARR2(1.55D-11,540.0_dp,TEMP))
  RCONST(72) = (ARR2(2.55D-11,-409.0_dp,TEMP))
  RCONST(73) = (ARR2(2.8D-12,-181.0_dp,TEMP))
  RCONST(74) = (ARR2(1.95D+16,13543.0_dp,TEMP))
  RCONST(75) = (4.7D-12)
  RCONST(76) = (ARR2(1.95D+16,13543.0_dp,TEMP))
  RCONST(77) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(78) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(79) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(80) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(81) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(82) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(83) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(84) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(85) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(86) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(87) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(88) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(89) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(90) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(91) = (ARR2(6.0D-13,2058.0_dp,TEMP))
  RCONST(92) = (ARR2(1.4D-12,1900.0_dp,TEMP))
  RCONST(93) = (ARR2(6.0D-13,2058.0_dp,TEMP))
  RCONST(94) = (ARR2(1.4D-12,1900.0_dp,TEMP))
  RCONST(95) = (ARR2(1.4D-12,1900.0_dp,TEMP))
  RCONST(96) = (2.2D-11)
  RCONST(97) = (ARR2(2.0D-12,2923.0_dp,TEMP))
  RCONST(98) = (ARR2(1.0D-11,1895.0_dp,TEMP))
  RCONST(99) = (ARR2(3.23D-11,975.0_dp,TEMP))
  RCONST(100) = (5.81D-13)
  RCONST(101) = (ARR2(1.2D-14,2633.0_dp,TEMP))
  RCONST(102) = (ARR2(1.32D-14,2105.0_dp,TEMP))
  RCONST(103) = (ARR2(7.29D-15,1136.0_dp,TEMP))
  RCONST(104) = (ARR2(1.23D-14,2013.0_dp,TEMP))
  RCONST(105) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(106) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(107) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(108) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(109) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(110) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(111) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(112) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(113) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(114) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(115) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(116) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(117) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(118) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(119) = (ARR2(1.9D-13,-220.0_dp,TEMP))
  RCONST(120) = (ARR2(1.4D-13,-220.0_dp,TEMP))
  RCONST(121) = (ARR2(4.2D-14,-220.0_dp,TEMP))
  RCONST(122) = (ARR2(3.4D-14,-220.0_dp,TEMP))
  RCONST(123) = (ARR2(2.9D-14,-220.0_dp,TEMP))
  RCONST(124) = (ARR2(1.4D-13,-220.0_dp,TEMP))
  RCONST(125) = (ARR2(1.4D-13,-220.0_dp,TEMP))
  RCONST(126) = (ARR2(1.7D-14,-220.0_dp,TEMP))
  RCONST(127) = (ARR2(1.7D-14,-220.0_dp,TEMP))
  RCONST(128) = (ARR2(9.6D-13,-220.0_dp,TEMP))
  RCONST(129) = (ARR2(1.7D-14,-220.0_dp,TEMP))
  RCONST(130) = (ARR2(1.7D-14,-220.0_dp,TEMP))
  RCONST(131) = (ARR2(9.6D-13,-220.0_dp,TEMP))
  RCONST(132) = (ARR2(3.4D-13,-220.0_dp,TEMP))
  RCONST(133) = (ARR2(1.0D-13,-220.0_dp,TEMP))
  RCONST(134) = (ARR2(8.4D-14,-220.0_dp,TEMP))
  RCONST(135) = (ARR2(7.2D-14,-220.0_dp,TEMP))
  RCONST(136) = (ARR2(3.4D-13,-220.0_dp,TEMP))
  RCONST(137) = (ARR2(3.4D-13,-220.0_dp,TEMP))
  RCONST(138) = (ARR2(4.2D-14,-220.0_dp,TEMP))
  RCONST(139) = (ARR2(4.2D-14,-220.0_dp,TEMP))
  RCONST(140) = (ARR2(1.19D-12,-220.0_dp,TEMP))
  RCONST(141) = (ARR2(4.2D-14,-220.0_dp,TEMP))
  RCONST(142) = (ARR2(4.2D-14,-220.0_dp,TEMP))
  RCONST(143) = (ARR2(1.19D-12,-220.0_dp,TEMP))
  RCONST(144) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(145) = (ARR2(1.7D-14,-220.0_dp,TEMP))
  RCONST(146) = (ARR2(4.2D-14,-220.0_dp,TEMP))
  RCONST(147) = (ARR2(3.6D-16,-220.0_dp,TEMP))
  RCONST(148) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(149) = (ARR2(4.2D-12,-180.0_dp,TEMP))
  RCONST(150) = (ARR2(7.7D-14,-1300.0_dp,TEMP))
  RCONST(151) = (ARR2(1.7D-14,-220.0_dp,TEMP))
  RCONST(152) = (ARR2(4.2D-14,-220.0_dp,TEMP))
  RCONST(153) = (ARR2(3.6D-16,-220.0_dp,TEMP))
  RCONST(154) = (ARR2(1.7D-14,-220.0_dp,TEMP))
  RCONST(155) = (ARR2(4.2D-14,-220.0_dp,TEMP))
  RCONST(156) = (ARR2(3.6D-16,-220.0_dp,TEMP))
END SUBROUTINE radm2_Update_RCONST













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











REAL(KIND=dp) FUNCTION k46( TEMP, C_M )
    REAL(KIND=dp), INTENT(IN) :: temp, c_m
    REAL(KIND=dp) :: k0, k2, k3 

   k0=7.2E-15_dp * EXP(785._dp/TEMP)
   k2=4.1E-16_dp * EXP(1440._dp/TEMP)
   k3=1.9E-33_dp * EXP(725._dp/TEMP)  * C_M

   k46=k0+k3/(1+k3/k2)


END FUNCTION k46







END MODULE radm2_UpdateRconstWRF

