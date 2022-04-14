























MODULE cb05_sorg_aq_UpdateRconstWRF

  USE cb05_sorg_aq_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE cb05_sorg_aq_Update_RCONST(  &



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










  RCONST(1) = (j(Pj_no2))
  RCONST(2) = (0.21*C_M*C_M*ARR(6.0D-34,0.0_dp,-2.4_dp,TEMP))
  RCONST(3) = (ARR2(3.0D-12,1500.0_dp,TEMP))
  RCONST(4) = (ARR2(5.6D-12,-180.0_dp,TEMP))
  RCONST(5) = (TROE(2.5D-31,1.8_dp,2.2D-11,0.7_dp,TEMP,C_M))
  RCONST(6) = (TROE(9.0D-32,1.5_dp,3.0D-11,0.0_dp,TEMP,C_M))
  RCONST(7) = (ARR2(1.2D-13,2450.0_dp,TEMP))
  RCONST(8) = (j(Pj_o33p))
  RCONST(9) = (j(Pj_o31d))
  RCONST(10) = (ARR2(2.1D-11,-102.0_dp,TEMP))
  RCONST(11) = (2.2D-10)
  RCONST(12) = (ARR2(1.7D-12,940.0_dp,TEMP))
  RCONST(13) = (ARR2(1.0D-14,490.0_dp,TEMP))
  RCONST(14) = (j(Pj_no3o))
  RCONST(15) = (j(Pj_no3o2))
  RCONST(16) = (ARR2(1.5D-11,-170.0_dp,TEMP))
  RCONST(17) = (ARR2(4.5D-14,1260.0_dp,TEMP))
  RCONST(18) = (TROE(2.0D-30,4.4_dp,1.4D-12,0.7_dp,TEMP,C_M))
  RCONST(19) = (2.5D-22)
  RCONST(20) = (1.8D-39)
  RCONST(21) = (FALL(1.0D-03,11000.0_dp,-3.5_dp,9.7D14,11080.0_dp,0.1_dp,0.45_dp,TEMP,C_M))
  RCONST(22) = (ARR2(3.3D-39,-530.0_dp,TEMP)*0.21*C_M)
  RCONST(23) = (5.0D-40)
  RCONST(24) = (TROE(7.0D-31,2.6_dp,3.6D-11,0.1_dp,TEMP,C_M))
  RCONST(25) = (j(Pj_hno2))
  RCONST(26) = (ARR2(1.8D-11,390.0_dp,TEMP))
  RCONST(27) = (1.0D-20)
  RCONST(28) = (TROE(2.0D-30,3.0_dp,2.5D-11,0.0_dp,TEMP,C_M))
  RCONST(29) = (EP2(2.4D-14,-460.0_dp,2.7D-17,-2199.0_dp,6.5D-34,-1335.0_dp,TEMP,C_M))
  RCONST(30) = (ARR2(3.5D-12,-250.0_dp,TEMP))
  RCONST(31) = (TROE(1.8D-31,3.2_dp,4.7D-12,0.0_dp,TEMP,C_M))
  RCONST(32) = (FALL(4.1D-05,10650.0_dp,0.0_dp,4.8D15,11170.0_dp,0.0_dp,0.6_dp,TEMP,C_M))
  RCONST(33) = (ARR2(1.3D-12,-380.0_dp,TEMP))
  RCONST(34) = (EP3(2.3D-13,-600.0_dp,1.7D-33,-1000.0_dp,TEMP,C_M))
  RCONST(35) = (EP3(3.22D-34,-2800.0_dp,2.38D-54,-3200.0_dp,TEMP,C_M))
  RCONST(36) = (j(Pj_h2o2))
  RCONST(37) = (ARR2(2.9D-12,160.0_dp,TEMP))
  RCONST(38) = (1.1D-10)
  RCONST(39) = (ARR2(5.5D-12,2000.0_dp,TEMP))
  RCONST(40) = (ARR2(2.2D-11,-120.0_dp,TEMP))
  RCONST(41) = (ARR2(4.2D-12,240.0_dp,TEMP))
  RCONST(42) = (TROE(6.9D-31,1.0_dp,2.6D-11,0.0_dp,TEMP,C_M))
  RCONST(43) = (ARR2(4.8D-11,-250.0_dp,TEMP))
  RCONST(44) = (ARR2(3.0D-11,-200.0_dp,TEMP))
  RCONST(45) = (ARR2(1.4D-12,2000.0_dp,TEMP))
  RCONST(46) = (1.0D-11)
  RCONST(47) = (2.2D-11)
  RCONST(48) = (3.5D-12)
  RCONST(49) = (1.0D-17)
  RCONST(50) = (ARR2(8.5D-13,2450.0_dp,TEMP))
  RCONST(51) = (j(Pj_hno4))
  RCONST(52) = (j(Pj_hno3))
  RCONST(53) = (j(Pj_n2o5))
  RCONST(54) = (ARR2(2.6D-12,-365.0_dp,TEMP))
  RCONST(55) = (ARR2(2.6D-12,-365.0_dp,TEMP))
  RCONST(56) = (ARR2(7.5D-13,-700.0_dp,TEMP))
  RCONST(57) = (ARR2(7.5D-13,-700.0_dp,TEMP))
  RCONST(58) = (6.8D-14)
  RCONST(59) = (6.8D-14)
  RCONST(60) = (6.8D-14)
  RCONST(61) = (ARR2(5.9D-13,360.0_dp,TEMP))
  RCONST(62) = (1.0D-4*j(Pj_no2))
  RCONST(63) = (ARR2(3.01D-12,-190.0_dp,TEMP))
  RCONST(64) = (0.7*j(Pj_h2o2))
  RCONST(65) = (EP3(1.44D-13,0.0_dp,3.43D-33,0.0_dp,TEMP,C_M))
  RCONST(66) = (ARR2(2.45D-12,1775.0_dp,TEMP))
  RCONST(67) = (ARR2(2.8D-12,-300.0_dp,TEMP))
  RCONST(68) = (ARR2(4.1D-13,-750.0_dp,TEMP))
  RCONST(69) = (ARR2(9.5D-14,-390.0_dp,TEMP))
  RCONST(70) = (ARR2(3.8D-12,-200.0_dp,TEMP))
  RCONST(71) = (0.7*j(Pj_h2o2))
  RCONST(72) = (ARR2(7.3D-12,620.0_dp,TEMP))
  RCONST(73) = (9.0D-12)
  RCONST(74) = (j(Pj_ch2or))
  RCONST(75) = (j(Pj_ch2om))
  RCONST(76) = (ARR2(3.4D-11,1600.0_dp,TEMP))
  RCONST(77) = (5.8D-16)
  RCONST(78) = (ARR2(9.7D-15,-625.0_dp,TEMP))
  RCONST(79) = (ARR2(2.4D+12,7000.0_dp,TEMP))
  RCONST(80) = (5.6D-12)
  RCONST(81) = (ARR2(5.6D-15,-2300.0_dp,TEMP))
  RCONST(82) = (4.0D-13)
  RCONST(83) = (ARR2(1.8D-11,1100.0_dp,TEMP))
  RCONST(84) = (ARR2(5.6D-12,-270.0_dp,TEMP))
  RCONST(85) = (ARR2(1.4D-12,1900.0_dp,TEMP))
  RCONST(86) = (4.6D-4*j(Pj_no2))
  RCONST(87) = (ARR2(8.1D-12,-270.0_dp,TEMP))
  RCONST(88) = (FALL(2.7D-28,0.0_dp,-7.1_dp,1.2D-11,0.0_dp,-0.9_dp,0.3_dp,TEMP,C_M))
  RCONST(89) = (FALL(4.9D-03,12100.0_dp,0.0_dp,5.4D16,13830.0_dp,0.0_dp,0.3_dp,TEMP,C_M))
  RCONST(90) = (j(Pj_pan))
  RCONST(91) = (ARR2(4.3D-13,-1040.0_dp,TEMP))
  RCONST(92) = (ARR2(2.0D-12,-500.0_dp,TEMP))
  RCONST(93) = (ARR2(4.4D-13,-1070.0_dp,TEMP))
  RCONST(94) = (ARR2(2.9D-12,-500.0_dp,TEMP))
  RCONST(95) = (ARR2(4.0D-13,-200.0_dp,TEMP))
  RCONST(96) = (0.0*0.7*j(Pj_h2o2))
  RCONST(97) = (ARR2(4.0D-13,-200.0_dp,TEMP))
  RCONST(98) = (ARR2(1.3D-11,870.0_dp,TEMP))
  RCONST(99) = (ARR2(5.1D-12,-405.0_dp,TEMP))
  RCONST(100) = (6.5D-15)
  RCONST(101) = (j(Pj_ch3cho))
  RCONST(102) = (ARR2(6.7D-12,-340.0_dp,TEMP))
  RCONST(103) = (FALL(2.7D-28,0.0_dp,-7.1_dp,1.2D-11,0.0_dp,-0.9_dp,0.3_dp,TEMP,C_M))
  RCONST(104) = (FALL(4.9D-03,12100.0_dp,0.0_dp,5.4D16,13830.0_dp,0.0_dp,0.3_dp,TEMP,C_M))
  RCONST(105) = (j(Pj_pan))
  RCONST(106) = (3.0D-13)
  RCONST(107) = (ARR2(4.3D-13,-1040.0_dp,TEMP))
  RCONST(108) = (ARR2(2.0D-12,-500.0_dp,TEMP))
  RCONST(109) = (ARR2(4.4D-13,-1070.0_dp,TEMP))
  RCONST(110) = (ARR2(2.9D-12,-500.0_dp,TEMP))
  RCONST(111) = (ARR2(2.9D-12,-500.0_dp,TEMP))
  RCONST(112) = (8.1D-13)
  RCONST(113) = (ARR2(1.0D15,8000.0_dp,TEMP))
  RCONST(114) = 1600
  RCONST(115) = (1.5D-11)
  RCONST(116) = (ARR2(1.0D-11,280.0_dp,TEMP))
  RCONST(117) = (3.2D-11)
  RCONST(118) = (ARR2(6.5D-15,1900.0_dp,TEMP))
  RCONST(119) = (ARR2(7.0D-13,2160.0_dp,TEMP))
  RCONST(120) = (ARR2(1.04D-11,792.0_dp,TEMP))
  RCONST(121) = (TROE(1.0D-28,0.8_dp,8.8D-12,0.0_dp,TEMP,C_M))
  RCONST(122) = (ARR2(1.2D-14,2630.0_dp,TEMP))
  RCONST(123) = (ARR2(3.3D-12,2880.0_dp,TEMP))
  RCONST(124) = (2.3D-11)
  RCONST(125) = (ARR2(1.0D-11,-550.0_dp,TEMP))
  RCONST(126) = (ARR2(8.4D-15,1100.0_dp,TEMP))
  RCONST(127) = (ARR2(9.6D-13,270.0_dp,TEMP))
  RCONST(128) = (ARR2(1.8D-12,-355.0_dp,TEMP))
  RCONST(129) = (8.1D-12)
  RCONST(130) = 4.2
  RCONST(131) = (4.1D-11)
  RCONST(132) = (2.2D-11)
  RCONST(133) = (1.4D-11)
  RCONST(134) = (5.5D-12)
  RCONST(135) = (9.0*j(Pj_ch2or))
  RCONST(136) = (3.0D-11)
  RCONST(137) = (ARR2(5.4D-17,500.0_dp,TEMP))
  RCONST(138) = (ARR2(1.7D-11,-116.0_dp,TEMP))
  RCONST(139) = (1.7D-11)
  RCONST(140) = (9.64*j(Pj_ch2or))
  RCONST(141) = (3.6D-11)
  RCONST(142) = (ARR2(2.54D-11,-407.6_dp,TEMP))
  RCONST(143) = (ARR2(7.86D-15,1912.0_dp,TEMP))
  RCONST(144) = (ARR2(3.03D-12,448.0_dp,TEMP))
  RCONST(145) = (3.36D-11)
  RCONST(146) = (7.1D-18)
  RCONST(147) = (1.0D-15)
  RCONST(148) = (0.0036*0.025*j(Pj_ch2om))
  RCONST(149) = (3.6D-11)
  RCONST(150) = (ARR2(1.5D-11,-449.0_dp,TEMP))
  RCONST(151) = (ARR2(1.2D-15,821.0_dp,TEMP))
  RCONST(152) = (ARR2(3.7D-12,-175.0_dp,TEMP))
  RCONST(153) = (TROE(3.3D-31,4.3_dp,1.6D-12,0.0_dp,TEMP,C_M))
  RCONST(154) = (ARR2(6.9D-12,230.0_dp,TEMP))
  RCONST(155) = (ARR2(8.7D-12,1070.0_dp,TEMP))
  RCONST(156) = (1.5D-19)
  RCONST(157) = (j(Pj_cl2))
  RCONST(158) = (j(Pj_hocl))
  RCONST(159) = (ARR2(2.3D-11,200.0_dp,TEMP))
  RCONST(160) = (1.63D-14)
  RCONST(161) = (ARR2(6.4D-12,-290.0_dp,TEMP))
  RCONST(162) = (ARR2(2.7D-12,-220.0_dp,TEMP))
  RCONST(163) = (5.0D-13)
  RCONST(164) = (j(Pj_fmcl))
  RCONST(165) = (ARR2(6.6D-12,1240.0_dp,TEMP))
  RCONST(166) = (5.0D-11)
  RCONST(167) = (ARR2(8.3D-11,100.0_dp,TEMP))
  RCONST(168) = (1.07D-10)
  RCONST(169) = (2.5D-10)
  RCONST(170) = (3.5D-10)
  RCONST(171) = (4.3D-10)
  RCONST(172) = (ARR2(8.2D-11,34.0_dp,TEMP))
  RCONST(173) = (7.9D-11)
  RCONST(174) = (1.3D-10)
  RCONST(175) = (5.5D-11)
  RCONST(176) = (ARR2(8.2D-11,-45.0_dp,TEMP))
  RCONST(177) = (ARR(6.58D-13,-58.0_dp,1.16_dp,TEMP))
  RCONST(178) = (3.0D-20)
  RCONST(179) = (8.7D-14)
  RCONST(180) = (8.5D-19)
  RCONST(181) = (2.93D-10)
  RCONST(182) = (1.71D-10)
  RCONST(183) = (2.52D-10)
  RCONST(184) = (5.37D-11)
  RCONST(185) = (8.66D-17)
  RCONST(186) = (7.89D-11)
  RCONST(187) = (1.36D-17)
  RCONST(188) = (2.31D-12)
  RCONST(189) = (2.7D-10)
  RCONST(190) = (1.97D-11)
  RCONST(191) = (7.7D-11)
END SUBROUTINE cb05_sorg_aq_Update_RCONST













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
   k3=1.9E-33_dp * EXP(725._dp/TEMP)

   k46=k0+k3/(1+k3/k2)


END FUNCTION k46




END MODULE cb05_sorg_aq_UpdateRconstWRF

