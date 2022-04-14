























MODULE mozart_mosaic_4bin_UpdateRconstWRF

  USE mozart_mosaic_4bin_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE mozart_mosaic_4bin_Update_RCONST(  &


   rh, &
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


   real(kind=dp), intent(in) :: rh
   REAL(KIND=dp), INTENT(IN) ::  nume,den








  RCONST(1) = (.20946_dp*j(Pj_o2))
  RCONST(2) = (j(Pj_o31d))
  RCONST(3) = (j(Pj_o33p))
  RCONST(4) = (j(Pj_n2o))
  RCONST(5) = (j(Pj_no2))
  RCONST(6) = (j(Pj_n2o5))
  RCONST(7) = (j(Pj_hno3))
  RCONST(8) = (j(Pj_no3o))
  RCONST(9) = (j(Pj_hno4))
  RCONST(10) = (j(Pj_ch3o2h))
  RCONST(11) = (j(Pj_ch2or))
  RCONST(12) = (j(Pj_ch2om))
  RCONST(13) = (j(Pj_h2o2))
  RCONST(14) = (j(Pj_ch3cho))
  RCONST(15) = (j(Pj_pooh))
  RCONST(16) = (.28_dp*j(Pj_h2o2))
  RCONST(17) = (j(Pj_pan))
  RCONST(18) = (j(Pj_pan))
  RCONST(19) = (j(Pj_macr))
  RCONST(20) = (j(Pj_mvk))
  RCONST(21) = (j(Pj_ch3o2h))
  RCONST(22) = (j(Pj_ch3o2h))
  RCONST(23) = (j(Pj_ch3o2h))
  RCONST(24) = (j(Pj_ch3coch3))
  RCONST(25) = (j(Pj_ch3cocho))
  RCONST(26) = (j(Pj_ch3o2h))
  RCONST(27) = (j(Pj_ch3cho))
  RCONST(28) = (j(Pj_ch3o2h))
  RCONST(29) = (j(Pj_hyac))
  RCONST(30) = (j(Pj_glyald))
  RCONST(31) = (j(Pj_mek))
  RCONST(32) = (.2_dp*j(Pj_no2))
  RCONST(33) = (j(Pj_gly))
  RCONST(34) = (j(Pj_ch3o2h))
  RCONST(35) = (j(Pj_ch3o2h))
  RCONST(36) = (j(Pj_ch3o2h))
  RCONST(37) = (0.140*j(Pj_no2))
  RCONST(38) = (0.100*j(Pj_no2))
  RCONST(39) = (0.100*j(Pj_no2))
  RCONST(40) = (0.200*j(Pj_no2))
  RCONST(41) = (0.200*j(Pj_no2))
  RCONST(42) = (0.006*j(Pj_no2))
  RCONST(43) = (j(Pj_ch3o2h))
  RCONST(44) = (j(Pj_glyald))
  RCONST(45) = (j(Pj_ch3cho))
  RCONST(46) = (j(Pj_ch3cho))
  RCONST(47) = (j(Pj_ch3o2h))
  RCONST(48) = (j(Pj_ch3o2h))
  RCONST(49) = (j(Pj_hno2))
  RCONST(50) = (0.20946*(C_M*6.00e-34_dp*(TEMP/300._dp)**(-2.3_dp)))
  RCONST(51) = (ARR2(8.0e-12_dp,2060.0_dp,TEMP))
  RCONST(52) = (.79_dp*ARR2(2.1e-11_dp,-115.0_dp,TEMP)+.21_dp*ARR2(3.2e-11_dp,-70.0_dp,TEMP))
  RCONST(53) = (2.2e-10_dp)
  RCONST(54) = (1.1e-10_dp)
  RCONST(55) = (ARR2(5.5e-12_dp,2000.0_dp,TEMP))
  RCONST(56) = (ARR2(2.2e-11_dp,-120.0_dp,TEMP))
  RCONST(57) = (ARR2(3.0e-11_dp,-200.0_dp,TEMP))
  RCONST(58) = (ARR2(1.7e-12_dp,940.0_dp,TEMP))
  RCONST(59) = (ARR2(1.0e-14_dp,490.0_dp,TEMP))
  RCONST(60) = (usr9(temp,c_m,c_h2o))
  RCONST(61) = (ARR2(2.9e-12_dp,160.0_dp,TEMP))
  RCONST(62) = (ARR2(4.8e-11_dp,-250.0_dp,TEMP))
  RCONST(63) = (ARR2(4.2e-12_dp,240.0_dp,TEMP))
  RCONST(64) = (TROE(6.90e-31_dp,1.0_dp,2.60e-11_dp,0.0_dp,TEMP,C_M))
  RCONST(65) = (6.7e-11_dp)
  RCONST(66) = (4.9e-11_dp)
  RCONST(67) = (ARR2(3.5e-12_dp,-250.0_dp,TEMP))
  RCONST(68) = (ARR2(3.0e-12_dp,1500.0_dp,TEMP))
  RCONST(69) = (ARR2(5.6e-12_dp,-180.0_dp,TEMP))
  RCONST(70) = (ARR2(1.2e-13_dp,2450.0_dp,TEMP))
  RCONST(71) = (ARR2(2.3e-12_dp,-170.0_dp,TEMP))
  RCONST(72) = (TROE(2.e-30_dp,4.4_dp,1.4e-12_dp,.7_dp,TEMP,C_M))
  RCONST(73) = (TROEE(3.333e26_dp,10900._dp,2.2e-30_dp,4.4_dp,1.4e-12_dp,.7_dp,TEMP,C_M))
  RCONST(74) = (TROE(2.e-30_dp,3._dp,2.5e-11_dp,0._dp,TEMP,C_M))
  RCONST(75) = (usr5(TEMP,C_M))
  RCONST(76) = (ARR2(1.5e-11_dp,-170._dp,TEMP))
  RCONST(77) = (TROE(1.8e-31_dp,3.2_dp,4.7e-12_dp,1.4_dp,TEMP,C_M))
  RCONST(78) = (ARR2(1.3e-12_dp,-380._dp,TEMP))
  RCONST(79) = (TROEE(4.76e26_dp,10900._dp,1.8e-31_dp,3.2_dp,4.7e-12_dp,1.4_dp,TEMP,C_M))
  RCONST(80) = (usr16(rh,temp))
  RCONST(81) = (usr17(rh,temp))
  RCONST(82) = (usr17a(rh,temp))
  RCONST(83) = (ARR2(2.45e-12_dp,1775.0_dp,TEMP))
  RCONST(84) = (1.5e-10_dp)
  RCONST(85) = (ARR2(2.8e-12_dp,-300._dp,TEMP))
  RCONST(86) = (ARR2(5.e-13_dp,424._dp,TEMP))
  RCONST(87) = (ARR2(1.9e-14_dp,-706._dp,TEMP))
  RCONST(88) = (ARR2(4.1e-13_dp,-750._dp,TEMP))
  RCONST(89) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(90) = (ARR2(6.e-13_dp,2058._dp,TEMP))
  RCONST(91) = (9.e-12_dp)
  RCONST(92) = (usr8(temp,c_m))
  RCONST(93) = (TROE(1.e-28_dp,.8_dp,8.8e-12_dp,0._dp,TEMP,C_M))
  RCONST(94) = (ARR2(1.2e-14_dp,2630._dp,TEMP))
  RCONST(95) = (usr23(TEMP,C_M))
  RCONST(96) = (1.1e-11_dp)
  RCONST(97) = (ARR2(4.2e-12_dp,-180._dp,TEMP))
  RCONST(98) = (1e-14_dp)
  RCONST(99) = (ARR2(1.6e11_dp,4150._dp,TEMP))
  RCONST(100) = (ARR2(8.7e-12_dp,1070._dp,TEMP))
  RCONST(101) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(102) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(103) = (2.e-13_dp)
  RCONST(104) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(105) = (JPL_TROE(8.0e-27_dp,3.5_dp,3.e-11_dp,0._dp,.5_dp,TEMP,C_M))
  RCONST(106) = (ARR2(6.5e-15_dp,1900._dp,TEMP))
  RCONST(107) = (ARR2(4.6e-13_dp,1156._dp,TEMP))
  RCONST(108) = (ARR2(4.2e-12_dp,-180._dp,TEMP))
  RCONST(109) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(110) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(111) = (ARR2(5.6e-12_dp,-270._dp,TEMP))
  RCONST(112) = (ARR2(1.4e-12_dp,1900._dp,TEMP))
  RCONST(113) = (ARR2(8.1e-12_dp,-270._dp,TEMP))
  RCONST(114) = (TROE(8.5e-29_dp,6.5_dp,1.1e-11_dp,1._dp,TEMP,C_M))
  RCONST(115) = (ARR2(4.3e-13_dp,-1040._dp,TEMP))
  RCONST(116) = (ARR2(2.e-12_dp,-500._dp,TEMP))
  RCONST(117) = (1.e-12_dp)
  RCONST(118) = (TROEE(1.111e28_dp,14000._dp,8.5e-29_dp,6.5_dp,1.1e-11_dp,0._dp,TEMP,C_M))
  RCONST(119) = (ARR2(2.5e-12_dp,-500._dp,TEMP))
  RCONST(120) = (ARR2(1.e-11_dp,660._dp,TEMP))
  RCONST(121) = (ARR2(4.e-12_dp,-180._dp,TEMP))
  RCONST(122) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(123) = (ARR2(3.75e-13_dp,40._dp,TEMP))
  RCONST(124) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(125) = (3.82e-11_dp*exp(-2000._dp/temp)+1.33e-13_dp)
  RCONST(126) = (ARR2(2.9e-12_dp,-300._dp,TEMP))
  RCONST(127) = (ARR2(8.6e-13_dp,-700._dp,TEMP))
  RCONST(128) = (ARR2(2.e-12_dp,-500._dp,TEMP))
  RCONST(129) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(130) = (5.4e-11_dp)
  RCONST(131) = (ARR2(4.2e-12_dp,-180._dp,TEMP))
  RCONST(132) = (3.5e-12_dp)
  RCONST(133) = (ARR2(4.2e-12_dp,-180._dp,TEMP))
  RCONST(134) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(135) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(136) = (6.8e-13_dp)
  RCONST(137) = (ARR2(2.3e-12_dp,170._dp,TEMP))
  RCONST(138) = (ARR2(4.2e-12_dp,-180._dp,TEMP))
  RCONST(139) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(140) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(141) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(142) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(143) = (ARR2(2.54e-11_dp,-410._dp,TEMP))
  RCONST(144) = (ARR2(1.05e-14_dp,2000._dp,TEMP))
  RCONST(145) = (ARR2(4.4e-12_dp,-180._dp,TEMP))
  RCONST(146) = (2.4e-12_dp)
  RCONST(147) = (ARR2(8.e-13_dp,-700._dp,TEMP))
  RCONST(148) = (ARR2(1.52e-11_dp,-200._dp,TEMP))
  RCONST(149) = (ARR2(5.e-13_dp,-400._dp,TEMP))
  RCONST(150) = (1.4e-11_dp)
  RCONST(151) = (ARR2(4.13e-12_dp,-452._dp,TEMP))
  RCONST(152) = (ARR2(7.52e-16_dp,1521._dp,TEMP))
  RCONST(153) = (ARR2(1.86e-11_dp,-175._dp,TEMP))
  RCONST(154) = (ARR2(4.4e-15_dp,2500._dp,TEMP))
  RCONST(155) = (ARR2(2.7e-12_dp,-360._dp,TEMP))
  RCONST(156) = (ARR2(1.3e-13_dp,-360._dp,TEMP))
  RCONST(157) = (2.4e-12_dp)
  RCONST(158) = (ARR2(8.e-13_dp,-700._dp,TEMP))
  RCONST(159) = (ARR2(5.e-13_dp,-400._dp,TEMP))
  RCONST(160) = (1.4e-11_dp)
  RCONST(161) = (ARR2(2.3e-11_dp,-200._dp,TEMP))
  RCONST(162) = (ARR2(5.3e-12_dp,-360._dp,TEMP))
  RCONST(163) = (5.e-12_dp)
  RCONST(164) = (ARR2(4.3e-13_dp,-1040._dp,TEMP))
  RCONST(165) = (ARR2(2.e-12_dp,-500._dp,TEMP))
  RCONST(166) = (ARR2(4.6e-12_dp,-530._dp,TEMP))
  RCONST(167) = (ARR2(2.3e-12_dp,-530._dp,TEMP))
  RCONST(168) = (1.1e-11_dp*300._dp/(temp*c_m))
  RCONST(169) = (1.2221e17_dp*300._dp*exp(-14000._dp/temp)/(temp*c_m))
  RCONST(170) = (ARR2(2.3e-12_dp,193._dp,TEMP))
  RCONST(171) = (ARR2(4.7e-13_dp,-1220._dp,TEMP))
  RCONST(172) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(173) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(174) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(175) = (2.1e-12_dp)
  RCONST(176) = (2.8e-13_dp)
  RCONST(177) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(178) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(179) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(180) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(181) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(182) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(183) = (TROE(8.5e-29_dp,6.5_dp,1.1e-11_dp,1._dp,TEMP,C_M))
  RCONST(184) = (ARR2(7.5e-12_dp,-290._dp,TEMP))
  RCONST(185) = (ARR2(4.3e-13_dp,-1040._dp,TEMP))
  RCONST(186) = (ARR2(1.7e-12_dp,-352._dp,TEMP))
  RCONST(187) = (4.7e-11_dp)
  RCONST(188) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(189) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(190) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(191) = (ARR2(5.9e-12_dp,-225._dp,TEMP))
  RCONST(192) = (TROE(8.5e-29_dp,6.5_dp,1.1e-11_dp,1._dp,TEMP,C_M))
  RCONST(193) = (TROEE(1.111e28_dp,14000._dp,8.5e-29_dp,6.5_dp,1.1e-11_dp,0._dp,TEMP,C_M))
  RCONST(194) = (ARR2(7.5e-12_dp,-290._dp,TEMP))
  RCONST(195) = (ARR2(4.3e-13_dp,-1040._dp,TEMP))
  RCONST(196) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(197) = (ARR2(4.3e-13_dp,-1040._dp,TEMP))
  RCONST(198) = (ARR2(7.5e-12_dp,-290._dp,TEMP))
  RCONST(199) = (ARR2(4.3e-13_dp,-1040._dp,TEMP))
  RCONST(200) = (ARR2(7.5e-12_dp,-290._dp,TEMP))
  RCONST(201) = (TROE(8.5e-29_dp,6.5_dp,1.1e-11_dp,1._dp,TEMP,C_M))
  RCONST(202) = (TROE(8.5e-29_dp,6.5_dp,1.1e-11_dp,1._dp,TEMP,C_M))
  RCONST(203) = (1.7e-11_dp)
  RCONST(204) = (8.4e-11_dp)
  RCONST(205) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(206) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(207) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(208) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(209) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(210) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(211) = (ARR2(1.2e-11_dp,-440._dp,TEMP))
  RCONST(212) = (ARR2(1.6e-11_dp,-470._dp,TEMP))
  RCONST(213) = (ARR2(4.2e-11_dp,-400._dp,TEMP))
  RCONST(214) = (2.1e-10_dp)
  RCONST(215) = (2.0e-10_dp)
  RCONST(216) = (ARR2(6.3e-16_dp,580._dp,TEMP))
  RCONST(217) = (ARR2(1.7e-15_dp,1300._dp,TEMP))
  RCONST(218) = (ARR2(3.0e-15_dp,780._dp,TEMP))
  RCONST(219) = (4.7e-16_dp)
  RCONST(220) = (1.2e-14_dp)
  RCONST(221) = (1.1e-11_dp)
  RCONST(222) = (1.2e-11_dp)
  RCONST(223) = (1.9e-11_dp)
  RCONST(224) = (ARR2(4.2e-12_dp,-180._dp,TEMP))
  RCONST(225) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(226) = (ARR2(2.0e-12_dp,-500._dp,TEMP))
  RCONST(227) = (3.3e-11_dp)
  RCONST(228) = (2.3e-11_dp)
  RCONST(229) = (5.7e-11_dp)
  RCONST(230) = (1.0e-12_dp)
  RCONST(231) = (ARR2(4.2e-12_dp,-180._dp,TEMP))
  RCONST(232) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(233) = (ARR2(2.0e-12_dp,-500._dp,TEMP))
  RCONST(234) = (3.4e-11_dp)
  RCONST(235) = (ARR2(4.2e-12_dp,-180._dp,TEMP))
  RCONST(236) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(237) = (ARR2(2.0e-12_dp,-500._dp,TEMP))
  RCONST(238) = (2.4e-12_dp)
  RCONST(239) = (7.e-13_dp)
  RCONST(240) = (ARR2(3.03e-12_dp,446._dp,TEMP))
  RCONST(241) = (ARR2(2.7e-12_dp,-360._dp,TEMP))
  RCONST(242) = (2.4e-12_dp)
  RCONST(243) = (ARR2(8.e-13_dp,-700._dp,TEMP))
  RCONST(244) = (ARR2(8.4e-13_dp,-830._dp,TEMP))
  RCONST(245) = (ARR2(1.4e-12_dp,1860._dp,TEMP))
  RCONST(246) = (4.5e-11_dp)
  RCONST(247) = (ARR2(1.4e-12_dp,1860._dp,TEMP))
  RCONST(248) = (ARR2(1.86e-11_dp,-175._dp,TEMP))
  RCONST(249) = (ARR2(2.7e-12_dp,-360._dp,TEMP))
  RCONST(250) = (2.4e-12_dp)
  RCONST(251) = (ARR2(8.e-13_dp,-700._dp,TEMP))
  RCONST(252) = (ARR2(5.e-13_dp,-400._dp,TEMP))
  RCONST(253) = (ARR2(1.3e-12_dp,-640._dp,TEMP))
  RCONST(254) = (ARR2(1.9e-12_dp,-190._dp,TEMP))
  RCONST(255) = (THERMAL_T2(7.69e-17_dp,-253._dp,TEMP))
  RCONST(256) = (ARR2(7.3e-12_dp,620._dp,TEMP))
  RCONST(257) = (ARR2(6.9e-12_dp,230._dp,TEMP))
  RCONST(258) = (JPL_TROE(8.e-27_dp,3.5_dp,3.e-11_dp,0.0_dp,.5_dp,TEMP,C_M))
  RCONST(259) = (4.e-14_dp)
  RCONST(260) = (3.e-12_dp)
  RCONST(261) = (1.e-11_dp)
  RCONST(262) = (ARR2(9.6e-12_dp,234._dp,TEMP))
  RCONST(263) = (usr24(temp,c_m))
  RCONST(264) = (ARR2(1.9e-13_dp,-520._dp,TEMP))
  RCONST(265) = (ARR2(1.7e-12_dp,710._dp,TEMP))
  RCONST(266) = (usr26(rh,temp))
  RCONST(267) = (6.8e-14_dp)
  RCONST(268) = (ARR2(8.1e-12_dp,-610._dp,TEMP))
  RCONST(269) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(270) = (ARR2(3.75e-13_dp,40._dp,TEMP))
  RCONST(271) = (1.4e-11_dp)
  RCONST(272) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(273) = (ARR2(4.3e-13_dp,-1040._dp,TEMP))
  RCONST(274) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(275) = (ARR2(3.8e-12_dp,-200._dp,TEMP))
  RCONST(276) = (1e-17_dp)
  RCONST(277) = (ARR2(4.6e-14_dp,400._dp,TEMP))
  RCONST(278) = (ARR2(4.3e-13_dp,-1040._dp,TEMP))
  RCONST(279) = (ARR2(2.6e-12_dp,-365._dp,TEMP))
  RCONST(280) = (2.4e-12_dp)
  RCONST(281) = (TROE(5.5e-30_dp,0._dp,8.3e-13_dp,-2._dp,TEMP,C_M))
  RCONST(282) = 4.5e-13
  RCONST(283) = (ARR2(9.7e-15_dp,-625._dp,TEMP))
  RCONST(284) = (ARR2(2.4e12_dp,7000._dp,TEMP))
  RCONST(285) = (ARR2(2.6e-12_dp,-265._dp,TEMP))
  RCONST(286) = (ARR2(7.5e-13_dp,-700._dp,TEMP))
  RCONST(287) = (1.25D-11)
  RCONST(288) = (1.25D-11)
  RCONST(289) = (Keff(2.50D-11,-408.0_dp,0.0_dp,TEMP,nume,den,0.0104_dp,0.0078_dp))
  RCONST(290) = (Keff(1.2D-11,-440.0_dp,0.0_dp,TEMP,nume,den,0.036_dp,0.2065_dp))
  RCONST(291) = (Keff(1.6D-11,-470.0_dp,0.0_dp,TEMP,nume,den,0.036_dp,0.2065_dp))
  RCONST(292) = (Keff(4.2D-11,-400.0_dp,0.0_dp,TEMP,nume,den,0.036_dp,0.2065_dp))
END SUBROUTINE mozart_mosaic_4bin_Update_RCONST













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












REAL(kind=dp) FUNCTION JPL_TROE( k0_300K, n, kinf_300K, m, base, temp, cair )




    REAL(kind=dp), INTENT(IN) :: base      
    REAL(kind=dp), INTENT(IN) :: temp      
    REAL(kind=dp), INTENT(IN) :: cair      
    REAL(kind=dp), INTENT(IN) :: k0_300K   
    REAL(kind=dp), INTENT(IN) :: n         
    REAL(kind=dp), INTENT(IN) :: kinf_300K 
    REAL(kind=dp), INTENT(IN) :: m         




    REAL(kind=dp)  :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair 
    kinf_T  = kinf_300K * zt_help**(m)        
    k_ratio = k0_T/kinf_T

    JPL_TROE = k0_T/(1._dp + k_ratio)*base**(1._dp/(1._dp + LOG10(k_ratio)**2))

END FUNCTION JPL_TROE

REAL(KIND=dp) FUNCTION usr5( temp, c_m )

    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    REAL(KIND=dp) :: k0, k2

   k0 = c_m * 6.5e-34_dp * exp( 1335._dp/temp )
   k2 = exp( 2199._dp/temp )
   k0 = k0 /(1.0_dp + k0/(2.7e-17_dp*k2))
   k2 = exp( 460._dp/temp )

   usr5 = k0 + 2.4e-14_dp * k2

END FUNCTION usr5

REAL(KIND=dp) FUNCTION usr8( temp, c_m )

    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    REAL(KIND=dp), parameter :: boltz = 1.38044e-16_dp

    usr8 = 1.5e-13_dp * (1._dp + 6.e-7_dp*boltz*c_m*temp)

END FUNCTION usr8

REAL(KIND=dp) FUNCTION usr9( temp, c_m, c_h2o )

    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m
    REAL(KIND=dp), INTENT(IN) :: c_h2o

    REAL(KIND=dp) :: ko, kinf, fc

    if( c_h2o > 0._dp ) then
       ko   = 2.3e-13_dp * exp( 600._dp/temp )
       kinf = 1.7e-33_dp * c_m * exp( 1000._dp/temp )
       fc   = 1._dp/c_h2o + 1.4e-21_dp * exp( 2200._dp/temp )
       usr9 = (ko + kinf) * fc
    else
       usr9 = 0._dp
    end if

END FUNCTION usr9

REAL(KIND=dp) FUNCTION usr16( rh, temp )

    REAL(KIND=dp), INTENT(IN) :: rh                       
    REAL(KIND=dp), INTENT(IN) :: temp                     


    usr16 = 0._dp

END FUNCTION usr16

REAL(KIND=dp) FUNCTION usr17( rh, temp )

    REAL(KIND=dp), INTENT(IN) :: rh                       
    REAL(KIND=dp), INTENT(IN) :: temp                     

    usr17 = 0._dp

END FUNCTION usr17

REAL(KIND=dp) FUNCTION usr17a( rh, temp )

    REAL(KIND=dp), INTENT(IN) :: rh                       
    REAL(KIND=dp), INTENT(IN) :: temp                     

    usr17a = 0._dp

END FUNCTION usr17a

REAL(KIND=dp) FUNCTION usr23( temp, c_m )

    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    REAL(KIND=dp) :: fc, k0
    REAL(KIND=dp) :: wrk

    fc    = 3.e-11_dp * (300._dp/temp) ** 3.3_dp
    wrk   = fc * c_m
    k0    = wrk / (1._dp + wrk/1.5e-12_dp)
    usr23 = k0 * .6_dp ** (1._dp/(1._dp + (log10( wrk/1.5e-12_dp ))**2._dp))

END FUNCTION usr23

REAL(KIND=dp) FUNCTION usr24( temp, c_m )

    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    REAL(KIND=dp) :: ko, wrk

    wrk   = .21_dp*c_m
    ko    = 1._dp + 5.5e-31_dp*exp( 7460._dp/temp )*wrk
    usr24 = 1.7e-42_dp*exp( 7810._dp/temp )*wrk/ko

END FUNCTION usr24

REAL(KIND=dp) FUNCTION usr26( rh, temp )

    REAL(KIND=dp), INTENT(IN) :: rh                       
    REAL(KIND=dp), INTENT(IN) :: temp                     

    usr26 = 0._dp

END FUNCTION usr26





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








END MODULE mozart_mosaic_4bin_UpdateRconstWRF

