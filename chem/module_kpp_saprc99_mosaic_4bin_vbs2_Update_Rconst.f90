























MODULE saprc99_mosaic_4bin_vbs2_UpdateRconstWRF

  USE saprc99_mosaic_4bin_vbs2_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE saprc99_mosaic_4bin_vbs2_Update_RCONST(  &

nume,den, &


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

REAL(KIND=dp), INTENT(IN) ::  nume,den









  RCONST(1) = (j(Pj_no2))
  RCONST(2) = (0.20946e0*C_M*ARR(5.68D-34,0.0_dp,-2.80_dp,TEMP))
  RCONST(3) = (ARR(8.00D-12,2060.0_dp,0.0_dp,TEMP))
  RCONST(4) = (ARR(1.00D-31,0.0_dp,-1.60_dp,TEMP))
  RCONST(5) = (ARR(6.50D-12,-120.0_dp,0.0_dp,TEMP))
  RCONST(6) = (FALL(9.00D-32,0.0_dp,-2.00_dp,2.20D-11,0.0_dp,0.0_dp,0.80_dp,TEMP,C_M))
  RCONST(7) = (ARR(1.80D-12,1370.0_dp,0.0_dp,TEMP))
  RCONST(8) = (ARR(1.40D-13,2470.0_dp,0.0_dp,TEMP))
  RCONST(9) = (ARR(1.80D-11,-110.0_dp,0.0_dp,TEMP))
  RCONST(10) = (0.20946e0*ARR(3.30D-39,-530.0_dp,0.0_dp,TEMP))
  RCONST(11) = (FALL(2.80D-30,0.0_dp,-3.50_dp,2.00D-12,0.0_dp,0.20_dp,0.45_dp,TEMP,C_M))
  RCONST(12) = (FALL(1.D-3,11000.0_dp,-3.5_dp,9.7D14,11080.0_dp,0.1_dp,0.45_dp,TEMP,C_M))
  RCONST(13) = (2.60D-22)
  RCONST(14) = (ARR(4.50D-14,1260.0_dp,0.0_dp,TEMP))
  RCONST(15) = (j(Pj_no3o2))
  RCONST(16) = (j(Pj_no3o))
  RCONST(17) = (j(Pj_o33p))
  RCONST(18) = (j(Pj_o31d))
  RCONST(19) = (2.20D-10)
  RCONST(20) = (ARR(2.09D-11,-95.0_dp,0.0_dp,TEMP))
  RCONST(21) = (FALL(7.00D-31,0.0_dp,-2.60_dp,3.60D-11,0.0_dp,-0.10_dp,0.60_dp,TEMP,C_M))
  RCONST(22) = (0.9000_dp*j(Pj_hno2))
  RCONST(23) = (0.1000_dp*j(Pj_hno2))
  RCONST(24) = (ARR(2.70D-12,-260.0_dp,0.0_dp,TEMP))
  RCONST(25) = (FALL(2.43D-30,0.0_dp,-3.10_dp,1.67D-11,0.0_dp,-2.10_dp,0.60_dp,TEMP,C_M))
  RCONST(26) = (2.00D-11)
  RCONST(27) = (EP2(7.20D-15,-785.0_dp,4.10D-16,-1440.0_dp,1.90D-33,-725.0_dp,TEMP,C_M))
  RCONST(28) = (j(Pj_hno3))
  RCONST(29) = (EP3(1.30D-13,0.0_dp,3.19D-33,0.0_dp,TEMP,C_M))
  RCONST(30) = (ARR(1.90D-12,1000.0_dp,0.0_dp,TEMP))
  RCONST(31) = (ARR(3.40D-12,-270.0_dp,0.0_dp,TEMP))
  RCONST(32) = (FALL(1.80D-31,0.0_dp,-3.20_dp,4.70D-12,0.0_dp,0.0_dp,0.60_dp,TEMP,C_M))
  RCONST(33) = (FALL(4.10D-05,10650.0_dp,0.0_dp,5.7D15,11170.0_dp,0.0_dp,0.5_dp,TEMP,C_M))
  RCONST(34) = (j(Pj_hno4))
  RCONST(35) = (ARR(1.50D-12,-360.0_dp,0.0_dp,TEMP))
  RCONST(36) = (ARR(1.40D-14,600.0_dp,0.0_dp,TEMP))
  RCONST(37) = (EP3(2.20D-13,-600.0_dp,1.85D-33,-980.0_dp,TEMP,C_M))
  RCONST(38) = (EP3(3.08D-34,-2800.0_dp,2.59D-54,-3180.0_dp,TEMP,C_M))
  RCONST(39) = (4.00D-12)
  RCONST(40) = (ARR(8.50D-13,2450.0_dp,0.0_dp,TEMP))
  RCONST(41) = (j(Pj_h2o2))
  RCONST(42) = (ARR(2.90D-12,160.0_dp,0.0_dp,TEMP))
  RCONST(43) = (ARR(4.80D-11,-250.0_dp,0.0_dp,TEMP))
  RCONST(44) = (FALL(4.00D-31,0.0_dp,-3.30_dp,2.00D-12,0.0_dp,0.0_dp,0.45_dp,TEMP,C_M))
  RCONST(45) = (5.40D-7*ARR(7.70D-12,2100.0_dp,0.0_dp,TEMP))
  RCONST(46) = (ARR(2.80D-12,-285.0_dp,0.0_dp,TEMP))
  RCONST(47) = (ARR(3.80D-13,-780.0_dp,0.0_dp,TEMP))
  RCONST(48) = (1.30D-12)
  RCONST(49) = (ARR(2.45D-14,-710.0_dp,0.0_dp,TEMP))
  RCONST(50) = (ARR(5.90D-13,509.0_dp,0.0_dp,TEMP))
  RCONST(51) = (ARR(2.70D-12,-360.0_dp,0.0_dp,TEMP))
  RCONST(52) = (ARR(1.90D-13,-1300.0_dp,0.0_dp,TEMP))
  RCONST(53) = (2.30D-12)
  RCONST(54) = (2.00D-13)
  RCONST(55) = (3.50D-14)
  RCONST(56) = (1.0_dp*ARR(2.70D-12,-360.0_dp,0.0_dp,TEMP))
  RCONST(57) = (1.0_dp*ARR(1.90D-13,-1300.0_dp,0.0_dp,TEMP))
  RCONST(58) = (1.0_dp*2.30D-12)
  RCONST(59) = (1.0_dp*2.00D-13)
  RCONST(60) = (1.0_dp*3.50D-14)
  RCONST(61) = (0.0_dp)
  RCONST(62) = (1.0_dp*ARR(2.70D-12,-360.0_dp,0.0_dp,TEMP))
  RCONST(63) = (1.0_dp*ARR(1.90D-13,-1300.0_dp,0.0_dp,TEMP))
  RCONST(64) = (1.0_dp*2.00D-13)
  RCONST(65) = (1.0_dp*2.30D-12)
  RCONST(66) = (1.0_dp*3.50D-14)
  RCONST(67) = (1.0_dp*3.50D-14)
  RCONST(68) = (1.0_dp*3.50D-14)
  RCONST(69) = (FALL(2.70D-28,0.0_dp,-7.10_dp,1.20D-11,0.0_dp,-0.90_dp,0.30_dp,TEMP,C_M))
  RCONST(70) = (FALL(4.90D-3,12100.0_dp,0.0_dp,4.0D16,13600.0_dp,0._dp,0.3_dp,TEMP,C_M))
  RCONST(71) = (ARR(7.80D-12,-300.0_dp,0.0_dp,TEMP))
  RCONST(72) = (ARR(4.30D-13,-1040.0_dp,0.0_dp,TEMP))
  RCONST(73) = (4.00D-12)
  RCONST(74) = (ARR(1.80D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(75) = (7.50D-12)
  RCONST(76) = (1.0_dp*7.50D-12)
  RCONST(77) = (1.0_dp*7.50D-12)
  RCONST(78) = (ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(79) = (ARR(1.20D-11,0.0_dp,-0.90_dp,TEMP))
  RCONST(80) = (ARR(2.00D15,12800.0_dp,0.0_dp,TEMP))
  RCONST(81) = (ARR(1.25D-11,-240.0_dp,0.0_dp,TEMP))
  RCONST(82) = (1.0_dp*ARR(4.30D-13,-1040.0_dp,0.0_dp,TEMP))
  RCONST(83) = (1.0_dp*4.00D-12)
  RCONST(84) = (1.0_dp*ARR(1.80D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(85) = (1.0_dp*7.50D-12)
  RCONST(86) = (1.0_dp*7.50D-12)
  RCONST(87) = (1.0_dp*7.50D-12)
  RCONST(88) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(89) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(90) = (1.37D-11)
  RCONST(91) = (ARR(7.90D16,14000.0_dp,0.0_dp,TEMP))
  RCONST(92) = (1.0_dp*ARR(1.25D-11,-240.0_dp,0.0_dp,TEMP))
  RCONST(93) = (1.0_dp*ARR(4.30D-13,-1040.0_dp,0.0_dp,TEMP))
  RCONST(94) = (1.0_dp*4.00D-12)
  RCONST(95) = (1.0_dp*ARR(1.80D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(96) = (1.0_dp*7.50D-12)
  RCONST(97) = (1.0_dp*7.50D-12)
  RCONST(98) = (1.0_dp*7.50D-12)
  RCONST(99) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(100) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(101) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(102) = (1.0_dp*ARR(1.20D-11,0.0_dp,-0.90_dp,TEMP))
  RCONST(103) = (ARR(1.60D16,13486.0_dp,0.0_dp,TEMP))
  RCONST(104) = (1.0_dp*ARR(1.25D-11,-240.0_dp,0.0_dp,TEMP))
  RCONST(105) = (1.0_dp*ARR(4.30D-13,-1040.0_dp,0.0_dp,TEMP))
  RCONST(106) = (1.0_dp*4.00D-12)
  RCONST(107) = (1.0_dp*ARR(1.80D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(108) = (1.0_dp*7.50D-12)
  RCONST(109) = (1.0_dp*7.50D-12)
  RCONST(110) = (1.0_dp*7.50D-12)
  RCONST(111) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(112) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(113) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(114) = (1.0_dp*ARR(2.90D-12,-500.0_dp,0.0_dp,TEMP))
  RCONST(115) = (2.40D-11)
  RCONST(116) = (ARR(7.50D14,8152.0_dp,0.0_dp,TEMP))
  RCONST(117) = (ARR(2.30D-11,-150.0_dp,0.0_dp,TEMP))
  RCONST(118) = (1.0_dp*ARR(1.90D-13,-1300.0_dp,0.0_dp,TEMP))
  RCONST(119) = (1.00D-03)
  RCONST(120) = (1.0_dp*ARR(7.50D14,8152.0_dp,0.0_dp,TEMP))
  RCONST(121) = (1.0_dp*ARR(2.30D-11,-150.0_dp,0.0_dp,TEMP))
  RCONST(122) = (1.0_dp*1.0_dp*ARR(1.90D-13,-1300.0_dp,0.0_dp,TEMP))
  RCONST(123) = (j(Pj_ch2or))
  RCONST(124) = (j(Pj_ch2om))
  RCONST(125) = (ARR(8.60D-12,-20.0_dp,0.0_dp,TEMP))
  RCONST(126) = (ARR(9.70D-15,-625.0_dp,0.0_dp,TEMP))
  RCONST(127) = (ARR(2.40D12,7000.0_dp,0.0_dp,TEMP))
  RCONST(128) = (1.0_dp*ARR(2.80D-12,-285.0_dp,0.0_dp,TEMP))
  RCONST(129) = (ARR(2.00D-12,2431.0_dp,0.0_dp,TEMP))
  RCONST(130) = (ARR(5.60D-12,-310.0_dp,0.0_dp,TEMP))
  RCONST(131) = (0.1900_dp*j(Pj_ch2or))
  RCONST(132) = (ARR(1.40D-12,1860.0_dp,0.0_dp,TEMP))
  RCONST(133) = (2.00D-11)
  RCONST(134) = (0.6500_dp*j(Pj_ch2or))
  RCONST(135) = (ARR(1.40D-12,1771.0_dp,0.0_dp,TEMP))
  RCONST(136) = (ARR(1.10D-12,520.0_dp,0.0_dp,TEMP))
  RCONST(137) = (0.0230_dp*j(Pj_ch2or))
  RCONST(138) = (ARR(1.30D-12,25.0_dp,2.0_dp,TEMP))
  RCONST(139) = (0.0650_dp*j(Pj_ch2or))
  RCONST(140) = (ARR(3.10D-12,360.0_dp,2.0_dp,TEMP))
  RCONST(141) = (ARR(0.0_dp,0.0_dp,1.0_dp,TEMP))
  RCONST(142) = (ARR(2.90D-12,-190.0_dp,0.0_dp,TEMP))
  RCONST(143) = (j(Pj_ch3o2h))
  RCONST(144) = (1.10D-11)
  RCONST(145) = (0.2100_dp*j(Pj_ch2or))
  RCONST(146) = (j(Pj_hcochob))
  RCONST(147) = (0.2000_dp*j(Pj_hcocho))
  RCONST(148) = (1.10D-11)
  RCONST(149) = (ARR(2.80D-12,2376.0_dp,0.0_dp,TEMP))
  RCONST(150) = (1.3000_dp*j(Pj_ch3cocho))
  RCONST(151) = (1.50D-11)
  RCONST(152) = (ARR(1.40D-12,1895.0_dp,0.0_dp,TEMP))
  RCONST(153) = (2.3000_dp*j(Pj_ch3cocho))
  RCONST(154) = (2.63D-11)
  RCONST(155) = (3.78D-12)
  RCONST(156) = (4.20D-11)
  RCONST(157) = (1.37D-11)
  RCONST(158) = (1.0_dp*2.63D-11)
  RCONST(159) = (1.29D-11)
  RCONST(160) = (1.7000_dp*j(Pj_ch2or))
  RCONST(161) = (ARR(1.40D-12,1872.0_dp,0.0_dp,TEMP))
  RCONST(162) = (ARR(1.86D-11,-176.0_dp,0.0_dp,TEMP))
  RCONST(163) = (ARR(1.36D-15,2114.0_dp,0.0_dp,TEMP))
  RCONST(164) = (ARR(1.50D-12,1726.0_dp,0.0_dp,TEMP))
  RCONST(165) = (6.34D-12)
  RCONST(166) = (0.0470_dp*j(Pj_ch2om))
  RCONST(167) = (ARR(4.14D-12,-453.0_dp,0.0_dp,TEMP))
  RCONST(168) = (ARR(7.51D-16,1520.0_dp,0.0_dp,TEMP))
  RCONST(169) = (4.32D-12)
  RCONST(170) = (0.6300_dp*j(Pj_macr))
  RCONST(171) = (6.19D-11)
  RCONST(172) = (4.18D-18)
  RCONST(173) = (1.00D-13)
  RCONST(174) = (0.0038_dp*j(Pj_hcochest))
  RCONST(175) = (1.50D-11)
  RCONST(176) = (0.3000_dp*j(Pj_ch3coc2h5))
  RCONST(177) = (7.80D-12)
  RCONST(178) = (1.2000_dp*j(Pj_ch3ono2))
  RCONST(179) = (5.00D-11)
  RCONST(180) = (2.00D-18)
  RCONST(181) = (5.00D-11)
  RCONST(182) = (2.0000_dp*j(Pj_hcochest))
  RCONST(183) = (5.00D-11)
  RCONST(184) = (6.8000_dp*j(Pj_hcochest))
  RCONST(185) = (ARR(2.15D-12,1735.0_dp,0.0_dp,TEMP))
  RCONST(186) = (ARR(1.96D-12,-438.0_dp,0.0_dp,TEMP))
  RCONST(187) = (ARR(9.14D-15,2580.0_dp,0.0_dp,TEMP))
  RCONST(188) = (ARR(4.39D-13,2282.0_dp,2.0_dp,TEMP))
  RCONST(189) = (ARR(1.04D-11,792.0_dp,0.0_dp,TEMP))
  RCONST(190) = (ARR(2.50D-11,-408.0_dp,0.0_dp,TEMP))
  RCONST(191) = (ARR(7.86D-15,1912.0_dp,0.0_dp,TEMP))
  RCONST(192) = (ARR(3.03D-12,448.0_dp,0.0_dp,TEMP))
  RCONST(193) = (3.60D-11)
  RCONST(194) = (ARR(1.83D-11,-449.0_dp,0.0_dp,TEMP))
  RCONST(195) = (ARR(1.08D-15,821.0_dp,0.0_dp,TEMP))
  RCONST(196) = (ARR(3.66D-12,-175.0_dp,0.0_dp,TEMP))
  RCONST(197) = (3.27D-11)
  RCONST(198) = (ARR(1.83D-11,-449.0_dp,0.0_dp,TEMP))
  RCONST(199) = (ARR(1.08D-15,821.0_dp,0.0_dp,TEMP))
  RCONST(200) = (ARR(3.66D-12,-175.0_dp,0.0_dp,TEMP))
  RCONST(201) = (3.27D-11)
  RCONST(202) = (ARR(1.37D-12,498.0_dp,2.0_dp,TEMP))
  RCONST(203) = (ARR(0.0_dp,0.0_dp,1.0_dp,TEMP))
  RCONST(204) = (ARR(9.87D-12,671.0_dp,0.0_dp,TEMP))
  RCONST(205) = (ARR(1.019D-11,434.0_dp,0.0_dp,TEMP))
  RCONST(206) = (ARR(5.946D-12,91.0_dp,0.0_dp,TEMP))
  RCONST(207) = (ARR(1.112D-11,52.0_dp,0.0_dp,TEMP))
  RCONST(208) = (ARR(1.81D-12,-355.0_dp,0.0_dp,TEMP))
  RCONST(209) = (2.640D-11)
  RCONST(210) = (ARR(7.095D-12,-451.0_dp,0.0_dp,TEMP))
  RCONST(211) = (ARR(2.617D-15,1640.0_dp,0.0_dp,TEMP))
  RCONST(212) = (ARR(4.453D-14,376.0_dp,0.0_dp,TEMP))
  RCONST(213) = (ARR(1.074D-11,234.0_dp,0.0_dp,TEMP))
  RCONST(214) = (ARR(1.743D-11,-384.0_dp,0.0_dp,TEMP))
  RCONST(215) = (ARR(5.022D-16,461.0_dp,0.0_dp,TEMP))
  RCONST(216) = (7.265D-13)
  RCONST(217) = (2.085D-11)
  RCONST(218) = (2.20D-10)
  RCONST(219) = (2.20D-10)
  RCONST(220) = (2.20D-10)
  RCONST(221) = (2.20D-10)
  RCONST(222) = (2.20D-10)
  RCONST(223) = (2.20D-10)
  RCONST(224) = (7.0D-7)
  RCONST(225) = (2.20D-10)
  RCONST(226) = (2.20D-10)
  RCONST(227) = (2.20D-10)
  RCONST(228) = (2.20D-10)
  RCONST(229) = (7.0D-7)
  RCONST(230) = (Keff(5.946D-12,91.0_dp,0.0_dp,TEMP,nume,den,0.011_dp,0.022_dp))
  RCONST(231) = (Keff(1.112D-11,52.0_dp,0.0_dp,TEMP,nume,den,0.064_dp,0.128_dp))
  RCONST(232) = (Keff(7.095D-12,-451.0_dp,0.0_dp,TEMP,nume,den,0.0002_dp,0.0012_dp))
  RCONST(233) = (Keff(1.743D-11,-384.0_dp,0.0_dp,TEMP,nume,den,0.0009_dp,0.0073_dp))
  RCONST(234) = (Keff(1.81D-12,-355.0_dp,0.0_dp,TEMP,nume,den,0.0322_dp,0.1222_dp))
  RCONST(235) = (Keff2(2.640D-11,nume,den,0.016_dp,0.1672_dp))
  RCONST(236) = (Keff(2.50D-11,-408.0_dp,0.0_dp,TEMP,nume,den,0.0104_dp,0.0027_dp))
  RCONST(237) = (Keff(1.83D-11,-449.0_dp,0.0_dp,TEMP,nume,den,0.036_dp,0.2065_dp))
  RCONST(238) = (Keff(1.83D-11,-449.0_dp,0.0_dp,TEMP,nume,den,0.6912_dp,0.3403_dp))
  RCONST(239) = (0.0D0)
  RCONST(240) = (0.0D0)
  RCONST(241) = (0.0D0)
  RCONST(242) = (0.0D0)
  RCONST(243) = (0.0D0)
  RCONST(244) = (0.0D0)
  RCONST(245) = (0.0D0)
  RCONST(246) = (0.0D0)
  RCONST(247) = (0.5714D-11)
  RCONST(248) = (0.5714D-11)
  RCONST(249) = (0.5714D-11)
  RCONST(250) = (0.5714D-11)
END SUBROUTINE saprc99_mosaic_4bin_vbs2_Update_RCONST













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













    REAL(KIND=dp) FUNCTION Keff ( A0,B0,C0, TEMP,X1,X2,y1,y2 )
    REAL(KIND=dp),INTENT(IN) :: X1,X2,y1,y2
    REAL(KIND=dp),INTENT(IN) :: TEMP
    REAL(KIND=dp),INTENT(IN):: A0,B0,C0
    Keff = A0 * EXP(- B0 /TEMP ) &
      *(TEMP/300._dp)**C0*(y1*X1/(X1 + X2 + 1.0e-35) &
       +y2*(1-X1/(X1 + X2 + 1.0e-35)))
    END FUNCTION Keff


    REAL(KIND=dp) FUNCTION Keff2 ( C0,X1,X2,y1,y2 )
    REAL(KIND=dp),INTENT(IN) :: X1,X2,y1,y2
    REAL(KIND=dp),INTENT(IN):: C0
    Keff2 = C0*(y1*X1/(X1 + X2 + 1.0e-35) &
       +y2*(1-X1/(X1 + X2 + 1.0e-35 )))
    END FUNCTION Keff2





END MODULE saprc99_mosaic_4bin_vbs2_UpdateRconstWRF

