























MODULE racm_soa_vbs_UpdateRconstWRF

  USE racm_soa_vbs_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE racm_soa_vbs_Update_RCONST(  &



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
  RCONST(22) = (j(Pj_macr))
  RCONST(23) = (j(Pj_ch3coc2h5))
  RCONST(24) = ((C_M*6.00D-34*(TEMP/300.0)**(-2.4)))
  RCONST(25) = (ARR2(8.00D-12,2060.0_dp,TEMP))
  RCONST(26) = (.78084*ARR2(2.15D-11,-110.0_dp,TEMP)+.20946*ARR2(3.30D-11,-55.0_dp,TEMP))
  RCONST(27) = (ARR2(1.63D-10,-60.0_dp,TEMP))
  RCONST(28) = (ARR2(1.70D-12,940.0_dp,TEMP))
  RCONST(29) = (ARR2(1.0D-14,490.0_dp,TEMP))
  RCONST(30) = (ARR2(4.80D-11,-250.0_dp,TEMP))
  RCONST(31) = (1.8D-12)
  RCONST(32) = ((3.5D-13*EXP(430./TEMP)+1.7D-33*C_M*EXP(1000./TEMP)))
  RCONST(33) = ((4.9D-34*EXP(2630./TEMP)+2.38D-54*C_M*EXP(3200./TEMP)))
  RCONST(34) = (TROE(9.00D-32,1.5_dp,3.00D-11,0.0_dp,TEMP,C_M))
  RCONST(35) = (ARR2(5.1D-12,-210.0_dp,TEMP))
  RCONST(36) = (TROE(2.5D-31,1.8_dp,2.20D-11,0.7_dp,TEMP,C_M))
  RCONST(37) = (TROE(7.00D-31,2.6_dp,3.6D-11,0.1_dp,TEMP,C_M))
  RCONST(38) = (TROE(1.8D-30,3.0_dp,2.8D-11,0.0_dp,TEMP,C_M))
  RCONST(39) = (2.20D-11)
  RCONST(40) = (ARR2(3.50D-12,-250.0_dp,TEMP))
  RCONST(41) = (TROE(2.0D-31,3.4_dp,2.9D-12,1.1_dp,TEMP,C_M))
  RCONST(42) = (TROEE(4.76D26,10900.0_dp,2.0D-31,3.4_dp,2.9D-12,1.1_dp,TEMP,C_M))
  RCONST(43) = (3.50D-12)
  RCONST(44) = (ARR2(1.80D-11,390.0_dp,TEMP))
  RCONST(45) = (k45(TEMP,C_M))
  RCONST(46) = (ARR2(1.30D-12,-380.0_dp,TEMP))
  RCONST(47) = (ARR2(3.0D-12,1500.0_dp,TEMP))
  RCONST(48) = (ARR2(1.20D-13,2450.0_dp,TEMP))
  RCONST(49) = ((.20946D0*ARR2(3.30D-39,-530.0_dp,TEMP)))
  RCONST(50) = (ARR2(1.50D-11,-170.0_dp,TEMP))
  RCONST(51) = (ARR2(4.50D-14,1260.0_dp,TEMP))
  RCONST(52) = (TROE(2.0D-30,4.4_dp,1.4D-12,0.7_dp,TEMP,C_M))
  RCONST(53) = (TROEE(3.70D26,11000.0_dp,2.0D-30,4.4_dp,1.4D-12,0.7_dp,TEMP,C_M))
  RCONST(54) = (ARR2(8.50D-13,2450.0_dp,TEMP))
  RCONST(55) = ((5.31D-7*ARR2(2.8D-12,1800.0_dp,TEMP)))
  RCONST(56) = (TROE(3.3D-31,4.3_dp,1.6D-12,0.0_dp,TEMP,C_M))
  RCONST(57) = (k57(TEMP,C_M))
  RCONST(58) = (ARR2(5.60D-12,-270.0_dp,TEMP))
  RCONST(59) = (3.00D-12)
  RCONST(60) = (ARR2(2.45D-12,1775.0_dp,TEMP))
  RCONST(61) = (ARR2(8.7D-12,1070.0_dp,TEMP))
  RCONST(62) = (ARR2(5.26D-12,260.0_dp,TEMP))
  RCONST(63) = (ARR2(8.02D-12,155.0_dp,TEMP))
  RCONST(64) = (ARR2(1.64D-11,125.0_dp,TEMP))
  RCONST(65) = (TROE(1.0D-28,4.5_dp,8.8D-12,0.85_dp,TEMP,C_M))
  RCONST(66) = (ARR2(5.72D-12,-500.0_dp,TEMP))
  RCONST(67) = (ARR2(1.33D-11,-500.0_dp,TEMP))
  RCONST(68) = (ARR2(1.48D-11,-448.0_dp,TEMP))
  RCONST(69) = (ARR2(2.54D-11,-410.0_dp,TEMP))
  RCONST(70) = (ARR2(1.21D-11,-444.0_dp,TEMP))
  RCONST(71) = (1.71D-10)
  RCONST(72) = (ARR2(1.81D-12,-338.0_dp,TEMP))
  RCONST(73) = (ARR2(7.30D-12,-355.0_dp,TEMP))
  RCONST(74) = (6.8D-11)
  RCONST(75) = (ARR2(5.5D-12,-125.0_dp,TEMP))
  RCONST(76) = (ARR2(5.6D-12,-270.0_dp,TEMP))
  RCONST(77) = ((THERMAL_T2(5.68D-18,-92.0_dp,TEMP)))
  RCONST(78) = (3.00D-12)
  RCONST(79) = (1.15D-11)
  RCONST(80) = (1.72D-11)
  RCONST(81) = (.5*(4.13D-12*EXP(425./TEMP)+1.86D-11*EXP(175./TEMP)))
  RCONST(82) = (ARR2(2.80D-11,-175.0_dp,TEMP))
  RCONST(83) = (2.70D-10)
  RCONST(84) = (ARR2(3.8D-12,-200.0_dp,TEMP))
  RCONST(85) = (ARR2(3.40D-12,-190.0_dp,TEMP))
  RCONST(86) = (ARR2(3.8D-12,-200.0_dp,TEMP))
  RCONST(87) = (4.00D-14)
  RCONST(88) = (ARR2(3.25D-13,-500.0_dp,TEMP))
  RCONST(89) = (ARR2(5.31D-12,260.0_dp,TEMP))
  RCONST(90) = (ARR2(3.40D-13,1900.0_dp,TEMP))
  RCONST(91) = (ARR2(1.40D-12,1900.0_dp,TEMP))
  RCONST(92) = (ARR2(2.90D-12,1900.0_dp,TEMP))
  RCONST(93) = (ARR2(1.40D-12,1900.0_dp,TEMP))
  RCONST(94) = (3.00D-11)
  RCONST(95) = (ARR2(2.87D-13,1000.0_dp,TEMP))
  RCONST(96) = (2.20D-11)
  RCONST(97) = ((THERMAL_T2(4.88D-18,2282.0_dp,TEMP)))
  RCONST(98) = (ARR2(1.79D-13,450.0_dp,TEMP))
  RCONST(99) = (ARR2(8.64D-13,-450.0_dp,TEMP))
  RCONST(100) = (1.00D-13)
  RCONST(101) = (ARR2(3.03D-12,446.0_dp,TEMP))
  RCONST(102) = (ARR2(1.19D-12,-490.0_dp,TEMP))
  RCONST(103) = (1.22D-11)
  RCONST(104) = (ARR2(2.20D-14,500.0_dp,TEMP))
  RCONST(105) = (ARR2(1.2D-14,2630.0_dp,TEMP))
  RCONST(106) = (ARR2(4.33D-15,1800.0_dp,TEMP))
  RCONST(107) = (ARR2(4.40D-15,845.0_dp,TEMP))
  RCONST(108) = (ARR2(1.34D-14,2283.0_dp,TEMP))
  RCONST(109) = (ARR2(7.86D-15,1913.0_dp,TEMP))
  RCONST(110) = (ARR2(1.01D-15,732.0_dp,TEMP))
  RCONST(111) = (2.00D-16)
  RCONST(112) = (.5*(1.36D-15*EXP(-2112./TEMP)+7.51D-16*EXP(-1521./TEMP)))
  RCONST(113) = (2.00D-18)
  RCONST(114) = (ARR2(2.46D-15,1700.0_dp,TEMP))
  RCONST(115) = (2.00D-11)
  RCONST(116) = (1.00D-11)
  RCONST(117) = (3.60D-11)
  RCONST(118) = ((.20946D0*ARR2(1.66D-17,-1044.0_dp,TEMP)))
  RCONST(119) = (5.00D-11)
  RCONST(120) = (3.60D-11)
  RCONST(121) = ((.20946D0*ARR2(1.66D-17,-1044.0_dp,TEMP)))
  RCONST(122) = (1.00D-11)
  RCONST(123) = (3.60D-11)
  RCONST(124) = ((.20946D0*ARR2(1.66D-17,-1044.0_dp,TEMP)))
  RCONST(125) = (5.00D-11)
  RCONST(126) = (TROE(9.70D-29,5.6_dp,9.30D-12,1.5_dp,TEMP,C_M))
  RCONST(127) = (TROEE(1.11D28,14000.0_dp,9.70D-29,5.6_dp,9.30D-12,1.5_dp,TEMP,C_M))
  RCONST(128) = (TROE(9.70D-29,5.6_dp,9.30D-12,1.5_dp,TEMP,C_M))
  RCONST(129) = (TROEE(1.11D28,14000.0_dp,9.70D-29,5.6_dp,9.30D-12,1.5_dp,TEMP,C_M))
  RCONST(130) = (ARR2(2.8D-12,-300.0_dp,TEMP))
  RCONST(131) = (ARR2(2.6D-12,-365.0_dp,TEMP))
  RCONST(132) = (4.00D-12)
  RCONST(133) = (4.00D-12)
  RCONST(134) = (4.00D-12)
  RCONST(135) = (9.00D-12)
  RCONST(136) = (4.00D-12)
  RCONST(137) = (4.00D-12)
  RCONST(138) = (ARR2(2.43D-12,-360.0_dp,TEMP))
  RCONST(139) = (4.00D-12)
  RCONST(140) = (4.00D-12)
  RCONST(141) = (4.00D-12)
  RCONST(142) = (4.00D-12)
  RCONST(143) = (4.00D-12)
  RCONST(144) = (ARR2(8.1D-12,-270.0_dp,TEMP))
  RCONST(145) = (ARR2(8.1D-12,-270.0_dp,TEMP))
  RCONST(146) = (4.00D-12)
  RCONST(147) = (4.00D-12)
  RCONST(148) = (4.00D-12)
  RCONST(149) = (ARR2(4.1D-13,-750.0_dp,TEMP))
  RCONST(150) = (ARR2(7.4D-13,-700.0_dp,TEMP))
  RCONST(151) = (ARR2(1.66D-13,-1300.0_dp,TEMP))
  RCONST(152) = (ARR2(1.66D-13,-1300.0_dp,TEMP))
  RCONST(153) = (ARR2(1.66D-13,-1300.0_dp,TEMP))
  RCONST(154) = (ARR2(1.90D-13,-1300.0_dp,TEMP))
  RCONST(155) = (ARR2(1.66D-13,-1300.0_dp,TEMP))
  RCONST(156) = (ARR2(1.66D-13,-1300.0_dp,TEMP))
  RCONST(157) = (ARR2(2.05D-13,-1300.0_dp,TEMP))
  RCONST(158) = (1.50D-11)
  RCONST(159) = (1.50D-11)
  RCONST(160) = (ARR2(3.75D-13,-980.0_dp,TEMP))
  RCONST(161) = (ARR2(3.75D-13,-980.0_dp,TEMP))
  RCONST(162) = (ARR2(3.75D-13,-980.0_dp,TEMP))
  RCONST(163) = (4.3D-13*EXP(1040./TEMP)/(1.+0.027*EXP(660./TEMP)))
  RCONST(164) = (4.3D-13*EXP(1040./TEMP)/(1.+37.*EXP(-660./TEMP)))
  RCONST(165) = (4.3D-13*EXP(1040./TEMP)/(1.+0.027*EXP(660./TEMP)))
  RCONST(166) = (4.3D-13*EXP(1040./TEMP)/(1.+37.*EXP(-660./TEMP)))
  RCONST(167) = (ARR2(1.15D-13,-1300.0_dp,TEMP))
  RCONST(168) = (ARR2(1.66D-13,-1300.0_dp,TEMP))
  RCONST(169) = (ARR2(1.66D-13,-1300.0_dp,TEMP))
  RCONST(170) = (ARR2(9.5D-14,-390.0_dp,TEMP))
  RCONST(171) = (ARR2(1.18D-13,-158.0_dp,TEMP))
  RCONST(172) = (ARR2(9.46D-14,-431.0_dp,TEMP))
  RCONST(173) = (ARR2(1.00D-13,-467.0_dp,TEMP))
  RCONST(174) = (ARR2(4.34D-14,-633.0_dp,TEMP))
  RCONST(175) = (ARR2(1.71D-13,-708.0_dp,TEMP))
  RCONST(176) = (ARR2(1.46D-13,-708.0_dp,TEMP))
  RCONST(177) = (ARR2(9.18D-14,-708.0_dp,TEMP))
  RCONST(178) = (ARR2(1.36D-13,-708.0_dp,TEMP))
  RCONST(179) = (ARR2(3.56D-14,-708.0_dp,TEMP))
  RCONST(180) = (ARR2(3.56D-14,-708.0_dp,TEMP))
  RCONST(181) = (ARR2(3.56D-14,-708.0_dp,TEMP))
  RCONST(182) = (ARR2(3.56D-14,-708.0_dp,TEMP))
  RCONST(183) = (ARR2(3.56D-14,-708.0_dp,TEMP))
  RCONST(184) = (ARR2(1.8D-12,-500.0_dp,TEMP))
  RCONST(185) = (ARR2(2.0D-13,-500.0_dp,TEMP))
  RCONST(186) = (ARR2(1.8D-12,-500.0_dp,TEMP))
  RCONST(187) = (ARR2(2.0D-13,-500.0_dp,TEMP))
  RCONST(188) = (ARR2(6.91D-13,-508.0_dp,TEMP))
  RCONST(189) = (ARR2(1.60D-13,-708.0_dp,TEMP))
  RCONST(190) = (ARR2(9.68D-14,-708.0_dp,TEMP))
  RCONST(191) = (ARR2(1.03D-12,-211.0_dp,TEMP))
  RCONST(192) = (ARR2(6.90D-13,-460.0_dp,TEMP))
  RCONST(193) = (ARR2(5.59D-13,-522.0_dp,TEMP))
  RCONST(194) = (ARR2(2.47D-13,-683.0_dp,TEMP))
  RCONST(195) = (ARR2(9.48D-13,-765.0_dp,TEMP))
  RCONST(196) = (ARR2(8.11D-13,-765.0_dp,TEMP))
  RCONST(197) = (ARR2(5.09D-13,-765.0_dp,TEMP))
  RCONST(198) = (ARR2(7.60D-13,-765.0_dp,TEMP))
  RCONST(199) = (ARR2(7.40D-13,-765.0_dp,TEMP))
  RCONST(200) = (ARR2(7.40D-13,-765.0_dp,TEMP))
  RCONST(201) = (ARR2(7.40D-13,-765.0_dp,TEMP))
  RCONST(202) = (ARR2(7.40D-13,-765.0_dp,TEMP))
  RCONST(203) = (ARR2(7.40D-13,-765.0_dp,TEMP))
  RCONST(204) = (ARR2(2.5D-12,-500.0_dp,TEMP))
  RCONST(205) = (ARR2(2.5D-12,-500.0_dp,TEMP))
  RCONST(206) = (ARR2(7.51D-13,-565.0_dp,TEMP))
  RCONST(207) = (ARR2(8.85D-13,-765.0_dp,TEMP))
  RCONST(208) = (ARR2(5.37D-13,-765.0_dp,TEMP))
  RCONST(209) = (ARR2(7.00D-14,-1000.0_dp,TEMP))
  RCONST(210) = (ARR2(4.25D-14,-1000.0_dp,TEMP))
  RCONST(211) = (ARR2(2.96D-14,-1000.0_dp,TEMP))
  RCONST(212) = (1.20D-12)
  RCONST(213) = (1.20D-12)
  RCONST(214) = (1.20D-12)
  RCONST(215) = (1.20D-12)
  RCONST(216) = (1.20D-12)
  RCONST(217) = (1.20D-12)
  RCONST(218) = (1.20D-12)
  RCONST(219) = (1.20D-12)
  RCONST(220) = (3.2D-11)
  RCONST(221) = (1.20D-12)
  RCONST(222) = (1.20D-12)
  RCONST(223) = (1.20D-12)
  RCONST(224) = (1.20D-12)
  RCONST(225) = (1.20D-12)
  RCONST(226) = (4.00D-12)
  RCONST(227) = (4.00D-12)
  RCONST(228) = (1.20D-12)
  RCONST(229) = (1.20D-12)
  RCONST(230) = (1.20D-12)
  RCONST(231) = (ARR2(1.66D-13,-1300.0_dp,TEMP))
  RCONST(232) = (ARR2(5.99D-15,-1510.0_dp,TEMP))
  RCONST(233) = (ARR2(3.40D-14,-1560.0_dp,TEMP))
  RCONST(234) = (ARR2(7.13D-17,-2950.0_dp,TEMP))
  RCONST(235) = (4.00D-12)
  RCONST(236) = (1.20D-12)
  RCONST(237) = (2.00D-12)
  RCONST(238) = (1.00D-10)
  RCONST(239) = (1.30D-11)
  RCONST(240) = (ARR2(2.54D-12,-360.0_dp,TEMP))
  RCONST(241) = (ARR2(1.82D-13,-1300.0_dp,TEMP))
  RCONST(242) = (2.00D-12)
  RCONST(243) = (TROE(9.70D-29,5.6_dp,9.30D-12,1.5_dp,TEMP,C_M))
  RCONST(244) = (TROEE(1.11D28,14000.0_dp,9.70D-29,5.6_dp,9.30D-12,1.5_dp,TEMP,C_M))
  RCONST(245) = (2.52D-10)
  RCONST(246) = (5.60D-16)
  RCONST(247) = (2.20D-11)
  RCONST(248) = (ARR2(1.33D-11,-500.0_dp,TEMP))
  RCONST(249) = (ARR2(8.64D-13,-450.0_dp,TEMP))
  RCONST(250) = (ARR2(4.40D-15,845.0_dp,TEMP))
  RCONST(251) = (ARR2(1.0D-11,0.0_dp,TEMP))
  RCONST(252) = (ARR2(1.0D-11,0.0_dp,TEMP))
  RCONST(253) = (ARR2(1.0D-11,0.0_dp,TEMP))
  RCONST(254) = (ARR2(1.0D-11,0.0_dp,TEMP))
  RCONST(255) = (ARR2(1.0D-11,0.0_dp,TEMP))
  RCONST(256) = (ARR2(1.0D-11,0.0_dp,TEMP))
END SUBROUTINE racm_soa_vbs_Update_RCONST













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











REAL(KIND=dp) FUNCTION k45( TEMP, C_M )
    REAL(KIND=dp), INTENT(IN) :: temp, c_m
    REAL(KIND=dp) :: k0, k2, k3 

   k0=2.4E-14_dp * EXP(460._dp/TEMP)
   k2=2.7E-17_dp * EXP(2199._dp/TEMP)
   k3=6.5E-34_dp * EXP(1335._dp/TEMP) * c_m

   k45=k0+k3/(1+k3/k2)

END FUNCTION k45

REAL(kind=dp) FUNCTION k57( TEMP, C_M )

    INTRINSIC LOG10

    REAL(KIND=dp), INTENT(IN) :: temp      
    REAL(KIND=dp), INTENT(IN) :: c_m       
    REAL(KIND=dp) :: k0_300Kn   
    REAL(KIND=dp) :: nn         
    REAL(KIND=dp) :: kinf_300Kn 
    REAL(KIND=dp) :: mn         
    REAL(KIND=dp) :: zt_help, k0_T, kinf_T, k_ratio
    REAL(KIND=dp) :: k57troe, k57cact

    k0_300Kn = 5.9e-33_dp
    nn = 1.4_dp
    kinf_300Kn = 1.1e-12_dp
    mn = -1.3_dp

    zt_help = 300._dp/temp
    k0_T    = k0_300Kn   * zt_help**(nn) * c_m 
    kinf_T  = kinf_300Kn * zt_help**(mn)       
    k_ratio = k0_T/kinf_T
    k57troe   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

    k0_300Kn = 1.5e-13_dp
    nn = -0.6_dp
    kinf_300Kn = 2.9e9_dp
    mn = -6.1_dp

    k0_T    = k0_300Kn   * zt_help**(nn)
    kinf_T  = kinf_300Kn * zt_help**(mn) / c_m 
    k_ratio = k0_T/kinf_T
    k57cact = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

    k57 = k57troe + k57cact 

END FUNCTION k57





END MODULE racm_soa_vbs_UpdateRconstWRF

