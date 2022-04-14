! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The KPP-WRF conform Reaction Rates File
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
! File                 : crimech_Update_Rconst.f90
! Time                 : Tue Apr 12 23:42:53 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/crimech
! Equation file        : crimech.kpp
! Output root filename : crimech
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE crimech_UpdateRconstWRF

  USE crimech_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE crimech_Update_RCONST(  &
!
#include <extra_args_update_rconst_crimech.inc>
!
#include <args_update_rconst.inc>
!
)


!
#include <decls_update_rconst.inc>
!
#include <extra_decls_update_rconst_crimech.inc>
!


! Begin INLINED RCONST


! End INLINED RCONST

  RCONST(1) = (j(Pj_o31d))
  RCONST(2) = (j(Pj_o33p))
  RCONST(3) = (j(Pj_h2o2))
  RCONST(4) = (j(Pj_no2))
  RCONST(5) = (j(Pj_no3o2))
  RCONST(6) = (j(Pj_no3o))
  RCONST(7) = (j(Pj_hno2))
  RCONST(8) = (j(Pj_hno3))
  RCONST(9) = (j(Pj_ch2or))
  RCONST(10) = (j(Pj_ch2om))
  RCONST(11) = (4.6D-4*j(Pj_no2))
  RCONST(12) = (4.19*4.6D-4*j(Pj_no2))
  RCONST(13) = (7.8D-5*j(Pj_no2))
  RCONST(14) = (7.047*7.8D-5*j(Pj_no2))
  RCONST(15) = (4.74*7.047*7.8D-5*j(Pj_no2))
  RCONST(16) = (1.33*7.047*7.8D-5*j(Pj_no2))
  RCONST(17) = (7.047*7.8D-5*j(Pj_no2))
  RCONST(18) = (7.047*7.8D-5*j(Pj_no2))
  RCONST(19) = (7.047*7.8D-5*j(Pj_no2))
  RCONST(20) = (3.00*7.047*7.8D-5*j(Pj_no2))
  RCONST(21) = (3.35*7.047*7.8D-5*j(Pj_no2))
  RCONST(22) = (9.64_dp*j(pj_ch2or))
  RCONST(23) = (2.0*1.9997*4.6D-4*j(Pj_no2))
  RCONST(24) = (9.64_dp*j(pj_ch2or))
  RCONST(25) = (32.6088*4.6D-4*j(Pj_no2))
  RCONST(26) = (2.149*32.6088*4.6D-4*j(Pj_no2))
  RCONST(27) = (2.149*32.6088*4.6D-4*j(Pj_no2))
  RCONST(28) = (2.149*32.6088*4.6D-4*j(Pj_no2))
  RCONST(29) = (2.0*1.9997*4.6D-4*j(Pj_no2))
  RCONST(30) = (1.9997*4.6D-4*j(Pj_no2))
  RCONST(31) = (1.155*4.6D-4*j(Pj_no2)+0.4933*4.6D-4*j(Pj_no2))
  RCONST(32) = (0.02*j(Pj_no2))
  RCONST(33) = (0.02*j(Pj_no2))
  RCONST(34) = (0.02*j(Pj_no2))
  RCONST(35) = (9.64_dp*j(pj_ch2or))
  RCONST(36) = (0.5*2.149*32.6088*4.6D-4*j(Pj_no2))
  RCONST(37) = (1.0D-4*j(Pj_no2))
  RCONST(38) = (2.3248*7.8D-5*j(Pj_no2))
  RCONST(39) = (3.079*7.8D-5*j(Pj_no2))
  RCONST(40) = (5.228*7.8D-5*j(Pj_no2))
  RCONST(41) = (0.398*3.079*7.8D-5*j(Pj_no2))
  RCONST(42) = (0.602*3.079*7.8D-5*j(Pj_no2))
  RCONST(43) = (3.079*7.8D-5*j(Pj_no2))
  RCONST(44) = (3.079*7.8D-5*j(Pj_no2))
  RCONST(45) = (5.228*7.8D-5*j(Pj_no2))
  RCONST(46) = (5.228*7.8D-5*j(Pj_no2))
  RCONST(47) = (5.228*7.8D-5*j(Pj_no2))
  RCONST(48) = (5.228*7.8D-5*j(Pj_no2))
  RCONST(49) = (0.7*j(Pj_h2o2))
  RCONST(50) = (0.7*j(Pj_h2o2))
  RCONST(51) = (0.7*j(Pj_h2o2))
  RCONST(52) = (0.7*j(Pj_h2o2))
  RCONST(53) = (0.398*0.7*j(Pj_h2o2))
  RCONST(54) = (0.602*0.7*j(Pj_h2o2))
  RCONST(55) = (0.7*j(Pj_h2o2))
  RCONST(56) = (0.7*j(Pj_h2o2))
  RCONST(57) = (0.7*j(Pj_h2o2))
  RCONST(58) = (0.7*j(Pj_h2o2))
  RCONST(59) = (0.7*j(Pj_h2o2))
  RCONST(60) = (0.7*j(Pj_h2o2))
  RCONST(61) = (0.7*j(Pj_h2o2))
  RCONST(62) = (0.7*j(Pj_h2o2))
  RCONST(63) = (0.7*j(Pj_h2o2))
  RCONST(64) = (0.252*0.7*j(Pj_h2o2))
  RCONST(65) = (0.748*0.7*j(Pj_h2o2))
  RCONST(66) = (0.7*j(Pj_h2o2))
  RCONST(67) = (0.7*j(Pj_h2o2))
  RCONST(68) = (0.7*j(Pj_h2o2))
  RCONST(69) = (0.7*j(Pj_h2o2))
  RCONST(70) = (0.7*j(Pj_h2o2))
  RCONST(71) = (0.7*j(Pj_h2o2))
  RCONST(72) = (0.7*j(Pj_h2o2))
  RCONST(73) = (0.7*j(Pj_h2o2))
  RCONST(74) = (0.7*j(Pj_h2o2))
  RCONST(75) = (0.7*j(Pj_h2o2))
  RCONST(76) = (0.7*j(Pj_h2o2))
  RCONST(77) = (0.7*j(Pj_h2o2))
  RCONST(78) = (0.7*j(Pj_h2o2))
  RCONST(79) = (0.7*j(Pj_h2o2))
  RCONST(80) = (0.7*j(Pj_h2o2))
  RCONST(81) = (0.7*j(Pj_h2o2))
  RCONST(82) = (0.7*j(Pj_h2o2))
  RCONST(83) = (0.7*j(Pj_h2o2))
  RCONST(84) = (0.7*j(Pj_h2o2))
  RCONST(85) = (0.7*j(Pj_h2o2))
  RCONST(86) = (0.7*j(Pj_h2o2))
  RCONST(87) = (0.7*j(Pj_h2o2))
  RCONST(88) = (0.7*j(Pj_h2o2))
  RCONST(89) = (0.7*j(Pj_h2o2))
  RCONST(90) = (0.7*j(Pj_h2o2))
  RCONST(91) = (0.7*j(Pj_h2o2))
  RCONST(92) = (0.7*j(Pj_h2o2))
  RCONST(93) = (0.02*0.36*j(Pj_no2))
  RCONST(94) = (0.02*0.45*j(Pj_no2))
  RCONST(95) = (0.02*0.45*j(Pj_no2))
  RCONST(96) = (.20946e0*(C_M*6.00D-34*(TEMP/300)**(-2.6)))
  RCONST(97) = (RJPL(1.3D-30,4.0_dp,7.5D-12,2.0_dp,C_M,TEMP)/(1.3D-28*exp(11200._dp/TEMP)))
  RCONST(98) = (ARR2(8.00D-12,2060.0_dp,TEMP))
  RCONST(99) = (TROE(9.00D-32,1.5_dp,3.00D-11,0.0_dp,TEMP,C_M))
  RCONST(100) = (ARR2(5.50D-12,-188.0_dp,TEMP))
  RCONST(101) = (TROE(9.00D-32,2.0_dp,2.20D-11,0.0_dp,TEMP,C_M))
  RCONST(102) = (.20946e0*ARR2(3.20D-11,-70.0_dp,TEMP)+.78084*ARR2(1.80D-11,-110.0_dp,TEMP))
  RCONST(103) = (ARR2(1.40D-12,1310.0_dp,TEMP))
  RCONST(104) = (ARR2(1.40D-13,2470.0_dp,TEMP))
  RCONST(105) = (.20946e0*ARR2(3.30D-39,-530.0_dp,TEMP))
  RCONST(106) = (ARR2(1.80D-11,-110.0_dp,TEMP))
  RCONST(107) = (ARR2(4.50D-14,1260.0_dp,TEMP))
  RCONST(108) = (TROE(2.20D-30,3.9_dp,1.50D-12,0.7_dp,TEMP,C_M))
  RCONST(109) = (TROEE(3.70D26,11000.0_dp,2.20D-30,3.9_dp,1.50D-12,0.7_dp,TEMP,C_M))
  RCONST(110) = (2.20D-10)
  RCONST(111) = (ARR2(1.70D-12,940.0_dp,TEMP))
  RCONST(112) = (5.31D-7*ARR2(7.70D-12,2100.0_dp,TEMP))
  RCONST(113) = (1.20D-13*(1.0+((0.6*C_M)/(2.652d+19*(273.0/temp)))))
  RCONST(114) = (ARR2(2.90D-12,160.0_dp,TEMP))
  RCONST(115) = (2.03D-16*((TEMP/300)**4.57)*EXP(693/TEMP))
  RCONST(116) = (ARR2(4.80D-11,-250.0_dp,TEMP))
  RCONST(117) = ((2.2D-13*EXP(600./TEMP)+1.9D-33*C_M*EXP(980._dp/TEMP)))
  RCONST(118) = ((3.08D-34*EXP(2800._dp/TEMP)+2.66D-54*C_M*EXP(3180._dp/TEMP)))
  RCONST(119) = (KMT_OH_NO(TEMP,C_M))
  RCONST(120) = (TROE(2.60D-30,3.2_dp,2.40D-11,1.3_dp,TEMP,C_M))
  RCONST(121) = (2.00D-11)
  RCONST(122) = (ARR2(3.60D-12,-270.0_dp,TEMP))
  RCONST(123) = (TROE(1.80D-31,3.2_dp,4.70D-12,1.4_dp,TEMP,C_M))
  RCONST(124) = (TROEE(4.76D26,10900.0_dp,1.80D-31,3.2_dp,4.70D-12,1.4_dp,TEMP,C_M))
  RCONST(125) = (ARR2(1.90D-12,-270.0_dp,TEMP))
  RCONST(126) = (4.00D-12)
  RCONST(127) = (ARR2(2.50D-12,-260.0_dp,TEMP))
  RCONST(128) = (k46(TEMP,C_M))
  RCONST(129) = (C_M*ARR2(4.00D-32,1000.0_dp,TEMP))
  RCONST(130) = (K47(TEMP,C_M))
  RCONST(131) = (.20946e0*ARR2(1.30D-12,-330.0_dp,TEMP))
  RCONST(132) = (ARR2(3.9d-41,-6830.6_dp,TEMP))
  RCONST(133) = (9.65D-20*TEMP**2.58*EXP(-1082/TEMP))
  RCONST(134) = (1.52D-17*TEMP**2*EXP(-498/TEMP))
  RCONST(135) = (1.55D-17*TEMP**2*EXP(-61/TEMP)*0.736)
  RCONST(136) = (1.55D-17*TEMP**2*EXP(-61/TEMP)*0.264)
  RCONST(137) = (1.69D-17*TEMP**2*EXP(145/TEMP))
  RCONST(138) = (KMT_IUPAC(8.6D-29,3.1_dp,9.0D-12,0.85_dp,0.48_dp,TEMP,C_M))
  RCONST(139) = (KMT_IUPAC(8.0D-27,3.5_dp,3.0D-11,1.0_dp,0.5_dp,TEMP,C_M))
  RCONST(140) = (ARR2(1.01D-11,-550.0_dp,TEMP))
  RCONST(141) = (2.10D-16)
  RCONST(142) = (9.40D-15)
  RCONST(143) = (3.90D-13)
  RCONST(144) = (0.13*ARR2(9.14D-15,2580.0_dp,TEMP))
  RCONST(145) = (0.87*ARR2(9.14D-15,2580.0_dp,TEMP))
  RCONST(146) = (0.36*ARR2(5.51D-15,1878.0_dp,TEMP))
  RCONST(147) = (0.64*ARR2(5.51D-15,1878.0_dp,TEMP))
  RCONST(148) = (0.69*ARR2(6.64D-15,1059.0_dp,TEMP))
  RCONST(149) = (0.31*ARR2(6.64D-15,1059.0_dp,TEMP))
  RCONST(150) = (ARR2(2.54D-11,-410.0_dp,TEMP))
  RCONST(151) = (ARR2(3.03D-12,446.0_dp,TEMP))
  RCONST(152) = (0.27*ARR2(7.86D-15,1913.0_dp,TEMP))
  RCONST(153) = (0.73*ARR2(7.86D-15,1913.0_dp,TEMP))
  RCONST(154) = (ARR2(1.20D-11,-444.0_dp,TEMP))
  RCONST(155) = (ARR2(1.19D-12,-490.0_dp,TEMP))
  RCONST(156) = (0.80*ARR2(1.01D-15,732.0_dp,TEMP))
  RCONST(157) = (0.075*ARR2(1.01D-15,732.0_dp,TEMP))
  RCONST(158) = (0.125*ARR2(1.01D-15,732.0_dp,TEMP))
  RCONST(159) = (ARR2(2.38D-11,-357.0_dp,TEMP))
  RCONST(160) = (2.51D-12)
  RCONST(161) = (1.50D-17*0.35)
  RCONST(162) = (1.50D-17*0.20)
  RCONST(163) = (1.50D-17*0.25)
  RCONST(164) = (1.50D-17*0.20)
  RCONST(165) = (0.364*KMT_IUPAC(5.0D-30,1.5_dp,1.0D-12,0.0_dp,0.37_dp,TEMP,C_M))
  RCONST(166) = (0.636*KMT_IUPAC(5.0D-30,1.5_dp,1.0D-12,0.0_dp,0.37_dp,TEMP,C_M))
  RCONST(167) = (0.47*ARR2(2.33D-12,193.0_dp,TEMP))
  RCONST(168) = (0.53*ARR2(2.33D-12,193.0_dp,TEMP))
  RCONST(169) = (0.82*ARR2(1.81D-12,-338.0_dp,TEMP))
  RCONST(170) = (0.18*ARR2(1.81D-12,-338.0_dp,TEMP))
  RCONST(171) = (1.36D-11*0.70)
  RCONST(172) = (1.36D-11*0.30)
  RCONST(173) = (1.20D-14*TEMP*EXP(287/TEMP))
  RCONST(174) = (ARR2(5.55D-12,-311.0_dp,TEMP))
  RCONST(175) = (1.96D-11)
  RCONST(176) = (5.80D-16)
  RCONST(177) = (1.44d-12*EXP(-1862.0/temp))
  RCONST(178) = (1.44d-12*EXP(-1862.0/temp)*2.4)
  RCONST(179) = (5.34D-18*TEMP**2*EXP(-230/TEMP))
  RCONST(180) = (3.24D-18*TEMP**2*EXP(414/TEMP))
  RCONST(181) = (6.01D-18*TEMP**2*EXP(170/TEMP))
  RCONST(182) = (6.18D-18*TEMP**2*EXP(532/TEMP)*0.887)
  RCONST(183) = (6.18D-18*TEMP**2*EXP(532/TEMP)*0.113)
  RCONST(184) = (5.53D-12*0.49)
  RCONST(185) = (5.53D-12*0.51)
  RCONST(186) = (4.06D-18*TEMP**2*EXP(788/TEMP)*0.86)
  RCONST(187) = (4.06D-18*TEMP**2*EXP(788/TEMP)*0.14)
  RCONST(188) = (4.50D-13)
  RCONST(189) = (8.00D-13)
  RCONST(190) = (0.999*ARR2(3.00D-12,-280.0_dp,TEMP))
  RCONST(191) = (0.991*ARR2(2.60D-12,-365.0_dp,TEMP))
  RCONST(192) = (0.980*ARR2(2.80D-12,-360.0_dp,TEMP))
  RCONST(193) = (0.958*ARR2(2.70D-12,-360.0_dp,TEMP))
  RCONST(194) = (2.40d-12*EXP(360.0/temp)*0.917*0.398)
  RCONST(195) = (2.40d-12*EXP(360.0/temp)*0.917*0.602)
  RCONST(196) = (2.40d-12*EXP(360.0/temp)*0.877)
  RCONST(197) = (2.40d-12*EXP(360.0/temp)*0.788)
  RCONST(198) = (2.40d-12*EXP(360.0/temp))
  RCONST(199) = (2.40d-12*EXP(360.0/temp))
  RCONST(200) = (2.40d-12*EXP(360.0/temp)*0.918)
  RCONST(201) = (2.40d-12*EXP(360.0/temp)*0.889*0.7)
  RCONST(202) = (2.40d-12*EXP(360.0/temp)*0.889*0.3)
  RCONST(203) = (2.40d-12*EXP(360.0/temp)*0.862)
  RCONST(204) = (2.40d-12*EXP(360.0/temp)*0.862)
  RCONST(205) = (2.40d-12*EXP(360.0/temp)*0.995*0.776)
  RCONST(206) = (2.40d-12*EXP(360.0/temp)*0.995*0.224)
  RCONST(207) = (2.40d-12*EXP(360.0/temp)*0.979)
  RCONST(208) = (2.40d-12*EXP(360.0/temp)*0.959)
  RCONST(209) = (2.40d-12*EXP(360.0/temp)*0.936)
  RCONST(210) = (2.40d-12*EXP(360.0/temp)*0.903)
  RCONST(211) = (2.40d-12*EXP(360.0/temp)*0.975)
  RCONST(212) = (2.40d-12*EXP(360.0/temp)*0.946)
  RCONST(213) = (8.10d-12*EXP(270.0/temp))
  RCONST(214) = (8.10d-12*EXP(270.0/temp))
  RCONST(215) = (8.10d-12*EXP(270.0/temp))
  RCONST(216) = (2.40d-12*EXP(360.0/temp))
  RCONST(217) = (2.40d-12*EXP(360.0/temp))
  RCONST(218) = (2.40d-12*EXP(360.0/temp))
  RCONST(219) = (2.40d-12*EXP(360.0/temp))
  RCONST(220) = (2.40d-12*EXP(360.0/temp)*0.900*0.252)
  RCONST(221) = (2.40d-12*EXP(360.0/temp)*0.900*0.748)
  RCONST(222) = (2.40d-12*EXP(360.0/temp)*0.7)
  RCONST(223) = (2.40d-12*EXP(360.0/temp)*0.3)
  RCONST(224) = (2.40d-12*EXP(360.0/temp)*0.5)
  RCONST(225) = (2.40d-12*EXP(360.0/temp)*0.3)
  RCONST(226) = (2.40d-12*EXP(360.0/temp)*0.2)
  RCONST(227) = (2.40d-12*EXP(360.0/temp))
  RCONST(228) = (2.40d-12*EXP(360.0/temp))
  RCONST(229) = (2.40d-12*EXP(360.0/temp))
  RCONST(230) = (2.40d-12*EXP(360.0/temp))
  RCONST(231) = (2.40d-12*EXP(360.0/temp))
  RCONST(232) = (2.40d-12*EXP(360.0/temp)*0.767*0.915)
  RCONST(233) = (2.40d-12*EXP(360.0/temp)*0.767*0.085)
  RCONST(234) = (2.40d-12*EXP(360.0/temp))
  RCONST(235) = (2.40d-12*EXP(360.0/temp))
  RCONST(236) = (2.40d-12*EXP(360.0/temp)*0.840)
  RCONST(237) = (2.40d-12*EXP(360.0/temp))
  RCONST(238) = (2.40d-12*EXP(360.0/temp))
  RCONST(239) = (2.40d-12*EXP(360.0/temp))
  RCONST(240) = (2.40d-12*EXP(360.0/temp))
  RCONST(241) = (2.40d-12*EXP(360.0/temp)*0.767*0.915)
  RCONST(242) = (2.40d-12*EXP(360.0/temp)*0.767*0.085)
  RCONST(243) = (2.40d-12*EXP(360.0/temp))
  RCONST(244) = (2.40d-12*EXP(360.0/temp)*0.843*0.6)
  RCONST(245) = (2.40d-12*EXP(360.0/temp)*0.843*0.4)
  RCONST(246) = (2.40d-12*EXP(360.0/temp)*0.700)
  RCONST(247) = (0.001*ARR2(3.00D-12,-280.0_dp,TEMP))
  RCONST(248) = (0.009*ARR2(2.60D-12,-365.0_dp,TEMP))
  RCONST(249) = (0.020*ARR2(2.80D-12,-360.0_dp,TEMP))
  RCONST(250) = (0.042*ARR2(2.70D-12,-360.0_dp,TEMP))
  RCONST(251) = (2.40d-12*EXP(360.0/temp)*0.083)
  RCONST(252) = (2.40d-12*EXP(360.0/temp)*0.123)
  RCONST(253) = (2.40d-12*EXP(360.0/temp)*0.212)
  RCONST(254) = (2.40d-12*EXP(360.0/temp)*0.005)
  RCONST(255) = (2.40d-12*EXP(360.0/temp)*0.021)
  RCONST(256) = (2.40d-12*EXP(360.0/temp)*0.041)
  RCONST(257) = (2.40d-12*EXP(360.0/temp)*0.064)
  RCONST(258) = (2.40d-12*EXP(360.0/temp)*0.097)
  RCONST(259) = (2.40d-12*EXP(360.0/temp)*0.025)
  RCONST(260) = (2.40d-12*EXP(360.0/temp)*0.054)
  RCONST(261) = (2.40d-12*EXP(360.0/temp)*0.100)
  RCONST(262) = (2.40d-12*EXP(360.0/temp)*0.082)
  RCONST(263) = (2.40d-12*EXP(360.0/temp)*0.111)
  RCONST(264) = (2.40d-12*EXP(360.0/temp)*0.138)
  RCONST(265) = (2.40d-12*EXP(360.0/temp)*0.138)
  RCONST(266) = (2.40d-12*EXP(360.0/temp)*0.233)
  RCONST(267) = (2.40d-12*EXP(360.0/temp)*0.160)
  RCONST(268) = (2.40d-12*EXP(360.0/temp)*0.233)
  RCONST(269) = (2.40d-12*EXP(360.0/temp)*0.157)
  RCONST(270) = (2.40d-12*EXP(360.0/temp)*0.300)
  RCONST(271) = (2.50d-12*0.40)
  RCONST(272) = (2.50d-12)
  RCONST(273) = (2.50d-12)
  RCONST(274) = (2.50d-12)
  RCONST(275) = (2.50d-12*0.398)
  RCONST(276) = (2.50d-12*0.602)
  RCONST(277) = (2.50d-12)
  RCONST(278) = (2.50d-12)
  RCONST(279) = (2.50d-12)
  RCONST(280) = (2.50d-12)
  RCONST(281) = (2.50d-12)
  RCONST(282) = (2.50d-12*0.7)
  RCONST(283) = (2.50d-12*0.3)
  RCONST(284) = (2.50d-12)
  RCONST(285) = (2.50d-12)
  RCONST(286) = (2.50d-12*0.776)
  RCONST(287) = (2.50d-12*0.224)
  RCONST(288) = (2.50d-12)
  RCONST(289) = (2.50d-12)
  RCONST(290) = (2.50d-12)
  RCONST(291) = (2.50d-12)
  RCONST(292) = (2.50d-12)
  RCONST(293) = (2.50d-12)
  RCONST(294) = (2.50d-12*1.60)
  RCONST(295) = (2.50d-12*1.60)
  RCONST(296) = (2.50d-12*1.60)
  RCONST(297) = (2.50d-12)
  RCONST(298) = (2.50d-12)
  RCONST(299) = (2.50d-12)
  RCONST(300) = (2.50d-12)
  RCONST(301) = (2.50d-12*0.252)
  RCONST(302) = (2.50d-12*0.748)
  RCONST(303) = (2.50d-12*0.7)
  RCONST(304) = (2.50d-12*0.3)
  RCONST(305) = (2.50d-12*0.5)
  RCONST(306) = (2.50d-12*0.3)
  RCONST(307) = (2.50d-12*0.2)
  RCONST(308) = (2.50d-12)
  RCONST(309) = (2.50d-12)
  RCONST(310) = (2.50d-12)
  RCONST(311) = (2.50d-12)
  RCONST(312) = (2.50d-12)
  RCONST(313) = (2.50d-12)
  RCONST(314) = (2.50d-12)
  RCONST(315) = (2.50d-12)
  RCONST(316) = (2.50d-12)
  RCONST(317) = (2.50d-12)
  RCONST(318) = (2.50d-12)
  RCONST(319) = (2.50d-12)
  RCONST(320) = (2.50d-12)
  RCONST(321) = (2.50d-12)
  RCONST(322) = (2.50d-12)
  RCONST(323) = (2.50d-12)
  RCONST(324) = (2.50d-12)
  RCONST(325) = (ARR2(3.80D-13,-780.0_dp,TEMP))
  RCONST(326) = (ARR2(7.50D-13,-700.0_dp,TEMP))
  RCONST(327) = (0.520*2.91d-13*EXP(1300.0/temp))
  RCONST(328) = (0.520*2.91d-13*EXP(1300.0/temp))
  RCONST(329) = (0.625*2.91d-13*EXP(1300.0/temp))
  RCONST(330) = (0.706*2.91d-13*EXP(1300.0/temp))
  RCONST(331) = (0.770*2.91d-13*EXP(1300.0/temp))
  RCONST(332) = (0.625*2.91d-13*EXP(1300.0/temp))
  RCONST(333) = (0.706*2.91d-13*EXP(1300.0/temp))
  RCONST(334) = (0.770*2.91d-13*EXP(1300.0/temp))
  RCONST(335) = (0.820*2.91d-13*EXP(1300.0/temp))
  RCONST(336) = (0.859*2.91d-13*EXP(1300.0/temp))
  RCONST(337) = (0.859*2.91d-13*EXP(1300.0/temp))
  RCONST(338) = (ARR2(2.03D-13,-1250.0_dp,TEMP))
  RCONST(339) = (0.520*2.91d-13*EXP(1300.0/temp))
  RCONST(340) = (0.625*2.91d-13*EXP(1300.0/temp))
  RCONST(341) = (0.706*2.91d-13*EXP(1300.0/temp))
  RCONST(342) = (0.770*2.91d-13*EXP(1300.0/temp))
  RCONST(343) = (0.706*2.91d-13*EXP(1300.0/temp))
  RCONST(344) = (0.770*2.91d-13*EXP(1300.0/temp))
  RCONST(345) = (4.30d-13*EXP(1040.0/temp))
  RCONST(346) = (4.30d-13*EXP(1040.0/temp))
  RCONST(347) = (4.30d-13*EXP(1040.0/temp))
  RCONST(348) = (0.520*2.91d-13*EXP(1300.0/temp))
  RCONST(349) = (0.625*2.91d-13*EXP(1300.0/temp))
  RCONST(350) = (0.706*2.91d-13*EXP(1300.0/temp))
  RCONST(351) = (0.770*2.91d-13*EXP(1300.0/temp))
  RCONST(352) = (0.770*2.91d-13*EXP(1300.0/temp))
  RCONST(353) = (0.706*2.91d-13*EXP(1300.0/temp))
  RCONST(354) = (0.625*2.91d-13*EXP(1300.0/temp))
  RCONST(355) = (0.387*2.91d-13*EXP(1300.0/temp))
  RCONST(356) = (0.520*2.91d-13*EXP(1300.0/temp))
  RCONST(357) = (0.625*2.91d-13*EXP(1300.0/temp))
  RCONST(358) = (0.770*2.91d-13*EXP(1300.0/temp))
  RCONST(359) = (0.625*2.91d-13*EXP(1300.0/temp))
  RCONST(360) = (0.914*2.91d-13*EXP(1300.0/temp))
  RCONST(361) = (0.914*2.91d-13*EXP(1300.0/temp))
  RCONST(362) = (0.914*2.91d-13*EXP(1300.0/temp))
  RCONST(363) = (0.890*2.91d-13*EXP(1300.0/temp))
  RCONST(364) = (0.890*2.91d-13*EXP(1300.0/temp))
  RCONST(365) = (0.890*2.91d-13*EXP(1300.0/temp))
  RCONST(366) = (0.770*2.91d-13*EXP(1300.0/temp))
  RCONST(367) = (0.706*2.91d-13*EXP(1300.0/temp))
  RCONST(368) = (0.914*2.91d-13*EXP(1300.0/temp))
  RCONST(369) = (0.890*2.91d-13*EXP(1300.0/temp))
  RCONST(370) = (0.890*2.91d-13*EXP(1300.0/temp))
  RCONST(371) = (0.914*2.91d-13*EXP(1300.0/temp))
  RCONST(372) = (0.33*RO2*ARR2(1.82D-13,-416.0_dp,TEMP))
  RCONST(373) = (0.335*RO2*ARR2(1.82D-13,-416.0_dp,TEMP))
  RCONST(374) = (0.335*RO2*ARR2(1.82D-13,-416.0_dp,TEMP))
  RCONST(375) = (3.10D-13*0.6*RO2)
  RCONST(376) = (3.10D-13*0.2*RO2)
  RCONST(377) = (3.10D-13*0.2*RO2)
  RCONST(378) = (6.00D-13*0.6*RO2)
  RCONST(379) = (6.00D-13*0.2*RO2)
  RCONST(380) = (6.00D-13*0.2*RO2)
  RCONST(381) = (4.00D-14*0.6*RO2)
  RCONST(382) = (4.00D-14*0.2*RO2)
  RCONST(383) = (4.00D-14*0.2*RO2)
  RCONST(384) = (2.50D-13*RO2*0.398)
  RCONST(385) = (2.50D-13*RO2*0.602)
  RCONST(386) = (8.80D-13*RO2)
  RCONST(387) = (8.80D-13*RO2)
  RCONST(388) = (8.80D-13*RO2)
  RCONST(389) = (8.80D-13*RO2*0.7)
  RCONST(390) = (8.80D-13*RO2*0.3)
  RCONST(391) = (8.80D-13*RO2)
  RCONST(392) = (8.80D-13*RO2)
  RCONST(393) = (2.50D-13*RO2)
  RCONST(394) = (2.50D-13*RO2)
  RCONST(395) = (2.00D-12*RO2*0.776)
  RCONST(396) = (2.00D-12*RO2*0.224)
  RCONST(397) = (8.80D-13*RO2)
  RCONST(398) = (8.80D-13*RO2)
  RCONST(399) = (8.80D-13*RO2)
  RCONST(400) = (8.80D-13*RO2)
  RCONST(401) = (8.80D-13*RO2)
  RCONST(402) = (8.80D-13*RO2)
  RCONST(403) = (1.00D-11*RO2)
  RCONST(404) = (1.00D-11*RO2)
  RCONST(405) = (1.00D-11*RO2)
  RCONST(406) = (1.40D-12*RO2)
  RCONST(407) = (1.40D-12*RO2)
  RCONST(408) = (1.40D-12*RO2)
  RCONST(409) = (1.40D-12*RO2)
  RCONST(410) = (1.71D-12*RO2*0.252)
  RCONST(411) = (1.71D-12*RO2*0.748)
  RCONST(412) = (2.00D-12*RO2*0.7)
  RCONST(413) = (2.00D-12*RO2*0.3)
  RCONST(414) = (2.00D-12*RO2*0.5)
  RCONST(415) = (2.00D-12*RO2*0.3)
  RCONST(416) = (2.00D-12*RO2*0.2)
  RCONST(417) = (6.00D-13*RO2)
  RCONST(418) = (2.30D-13*RO2)
  RCONST(419) = (2.50D-13*RO2)
  RCONST(420) = (1.30D-12*RO2)
  RCONST(421) = (9.60D-13*RO2)
  RCONST(422) = (2.85D-13*RO2)
  RCONST(423) = (1.00D-13*RO2)
  RCONST(424) = (2.00D-12*RO2)
  RCONST(425) = (1.30D-12*RO2)
  RCONST(426) = (6.70D-15*RO2)
  RCONST(427) = (6.70D-15*RO2)
  RCONST(428) = (8.80D-13*RO2)
  RCONST(429) = (2.00D-12*RO2)
  RCONST(430) = (2.00D-12*RO2)
  RCONST(431) = (2.50D-13*RO2)
  RCONST(432) = (2.50D-13*RO2)
  RCONST(433) = (9.20D-14*RO2)
  RCONST(434) = (1.87D-11)
  RCONST(435) = (4.36D-12)
  RCONST(436) = (3.24D-18*TEMP**2*EXP(414/TEMP))
  RCONST(437) = (3.00D-12)
  RCONST(438) = (5.86D-12)
  RCONST(439) = (1.65D-11)
  RCONST(440) = (1.25D-11)
  RCONST(441) = (2.50D-11)
  RCONST(442) = (1.44d-12*EXP(-1862.0/temp))
  RCONST(443) = (2.85D-18*0.59)
  RCONST(444) = (2.85D-18*0.41)
  RCONST(445) = (1.00D-11)
  RCONST(446) = (1.44d-12*EXP(-1862.0/temp))
  RCONST(447) = (1.14D-11)
  RCONST(448) = (1.72D-11)
  RCONST(449) = (2.40D-13)
  RCONST(450) = (1.38D-12)
  RCONST(451) = (4.81D-12)
  RCONST(452) = (4.79D-12)
  RCONST(453) = (4.52D-11)
  RCONST(454) = (1.44d-12*EXP(-1862.0/temp)*4.25)
  RCONST(455) = (2.40D-17*0.89)
  RCONST(456) = (2.40D-17*0.11)
  RCONST(457) = (4.16D-11)
  RCONST(458) = (1.30D-13)
  RCONST(459) = (5.20D-11*0.5)
  RCONST(460) = (5.58D-11*0.55)
  RCONST(461) = (7.00D-11*0.55)
  RCONST(462) = (4.20D-11)
  RCONST(463) = (1.00D-12)
  RCONST(464) = (1.00D-10)
  RCONST(465) = (3.80D-14)
  RCONST(466) = (1.44d-12*EXP(-1862.0/temp)*5.5)
  RCONST(467) = (6.65D-12)
  RCONST(468) = (1.55D-11)
  RCONST(469) = (4.55D-12)
  RCONST(470) = (ARR2(1.00D-14,-1060.0_dp,TEMP))
  RCONST(471) = (ARR2(4.40D-14,-720.0_dp,TEMP))
  RCONST(472) = (7.30D-13)
  RCONST(473) = (4.90D-13)
  RCONST(474) = (9.20D-13)
  RCONST(475) = (1.85D-12)
  RCONST(476) = (3.02D-12)
  RCONST(477) = (1.09D-12)
  RCONST(478) = (1.31D-12)
  RCONST(479) = (1.79D-12)
  RCONST(480) = (1.03D-11)
  RCONST(481) = (1.34D-11)
  RCONST(482) = (5.55D-11)
  RCONST(483) = (7.30D-11)
  RCONST(484) = (7.16D-11)
  RCONST(485) = (8.31D-11)
  RCONST(486) = (4.35D-12)
  RCONST(487) = (2.88D-12)
  RCONST(488) = (3.53D-12)
  RCONST(489) = (6.48D-12)
  RCONST(490) = (4.74D-12)
  RCONST(491) = (2.63D-11)
  RCONST(492) = (3.78D-12)
  RCONST(493) = (2.08D-12)
  RCONST(494) = (9.00D-13)
  RCONST(495) = (9.00D-14)
  RCONST(496) = (4.65D-11)
  RCONST(497) = (1.25D-11)
  RCONST(498) = (2.08D-12)
  RCONST(499) = (1.53D-12)
  RCONST(500) = (3.13D-13)
  RCONST(501) = (ARR2(1.90D-12,-190.0_dp,TEMP))
  RCONST(502) = (ARR2(1.00D-12,-190.0_dp,TEMP))
  RCONST(503) = (1.36D-11)
  RCONST(504) = (1.89D-11)
  RCONST(505) = (2.78D-11)
  RCONST(506) = (3.57D-11)
  RCONST(507) = (4.21D-11)
  RCONST(508) = (4.71D-11)
  RCONST(509) = (3.70D-12)
  RCONST(510) = (4.42D-12)
  RCONST(511) = (6.19D-12)
  RCONST(512) = (4.42D-12)
  RCONST(513) = (2.50D-11)
  RCONST(514) = (3.20D-11)
  RCONST(515) = (3.35D-11)
  RCONST(516) = (7.51D-11)
  RCONST(517) = (3.00D-11)
  RCONST(518) = (3.00D-11)
  RCONST(519) = (1.03D-10)
  RCONST(520) = (2.65D-11)
  RCONST(521) = (2.13D-11)
  RCONST(522) = (2.50D-11)
  RCONST(523) = (3.25D-11)
  RCONST(524) = (3.74D-11)
  RCONST(525) = (3.83D-11)
  RCONST(526) = (5.22D-12)
  RCONST(527) = (6.50D-12)
  RCONST(528) = (7.15D-12)
  RCONST(529) = (9.77D-11)
  RCONST(530) = (9.64D-11)
  RCONST(531) = (1.12D-10)
  RCONST(532) = (2.38D-11)
  RCONST(533) = (1.20D-11)
  RCONST(534) = (9.50D-12)
  RCONST(535) = (1.66D-11)
  RCONST(536) = (1.05D-11)
  RCONST(537) = (2.05D-11)
  RCONST(538) = (8.69D-11)
  RCONST(539) = (4.23D-12)
  RCONST(540) = (2.00D-11)
  RCONST(541) = (8.59D-11)
  RCONST(542) = (7.50D-11)
  RCONST(543) = (9.58D-12)
  RCONST(544) = (TROE(8.5d-29,6.5_dp,1.1d-11,1._dp,TEMP,C_M))
  RCONST(545) = (TROEE(1.111d28,14000._dp,8.5d-29,6.5_dp,1.1d-11,0._dp,TEMP,C_M))
  RCONST(546) = (TROE(8.5d-29,6.5_dp,1.1d-11,1._dp,TEMP,C_M))
  RCONST(547) = (TROEE(1.111d28,14000._dp,8.5d-29,6.5_dp,1.1d-11,0._dp,TEMP,C_M))
  RCONST(548) = (TROE(8.5d-29,6.5_dp,1.1d-11,1._dp,TEMP,C_M))
  RCONST(549) = (TROEE(1.111d28,14000._dp,8.5d-29,6.5_dp,1.1d-11,0._dp,TEMP,C_M))
  RCONST(550) = (ARR2(9.50D-13,650.0_dp,TEMP))
  RCONST(551) = (1.27D-12)
  RCONST(552) = (1.12D-12)
  RCONST(553) = (0.061*TROE(8.5d-29,6.5_dp,1.1d-11,1._dp,TEMP,C_M))
  RCONST(554) = (TROEE(1.111d28,14000._dp,8.5d-29,6.5_dp,1.1d-11,0._dp,TEMP,C_M))
  RCONST(555) = (0.041*TROE(8.5d-29,6.5_dp,1.1d-11,1._dp,TEMP,C_M))
  RCONST(556) = (TROEE(1.111d28,14000._dp,8.5d-29,6.5_dp,1.1d-11,0._dp,TEMP,C_M))
  RCONST(557) = (3.60D-12)
  RCONST(558) = (2.52D-11)
  RCONST(559) = (0.722*TROE(8.5d-29,6.5_dp,1.1d-11,1._dp,TEMP,C_M))
  RCONST(560) = (TROEE(1.111d28,14000._dp,8.5d-29,6.5_dp,1.1d-11,0._dp,TEMP,C_M))
  RCONST(561) = (3.66D-12)
  RCONST(562) = (1.50D-12)
  RCONST(563) = (5.20D-11*0.50)
  RCONST(564) = (5.58D-11*0.45)
  RCONST(565) = (7.00D-11*0.45)
  RCONST(566) = (7.33D-18*EXP(-809/TEMP)*TEMP**2)
  RCONST(567) = (6.14D-18*EXP(-389/TEMP)*TEMP**2)
  RCONST(568) = (1.80D-18*EXP(-129/TEMP)*TEMP**2)
  RCONST(569) = (2.25D-18*EXP(-910/TEMP)*TEMP**2)
  RCONST(570) = (ARR2(9.64D-12,1209.0_dp,TEMP))
  RCONST(571) = (ARR2(5.63D-13,-427.0_dp,TEMP))
  RCONST(572) = (ARR2(1.94D-12,-90.0_dp,TEMP))
  RCONST(573) = (ARR2(1.01D-12,-250.0_dp,TEMP))
  RCONST(574) = (RJPL(1.3D-30,4.0_dp,7.5D-12,2.0_dp,C_M,TEMP))
  RCONST(575) = (7.047*7.8D-5*j(Pj_no2))
  RCONST(576) = (32.6088*4.6D-4*j(Pj_no2))
  RCONST(577) = (5.228*7.8D-5*j(Pj_no2))
  RCONST(578) = (0.7*j(Pj_h2o2))
  RCONST(579) = (5.228*7.8D-5*j(Pj_no2))
  RCONST(580) = (0.7*j(Pj_h2o2))
  RCONST(581) = (0.02*0.45*j(Pj_no2))
  RCONST(582) = (0.02*0.55*j(Pj_no2))
  RCONST(583) = (2.40d-12*EXP(360.0/temp)*0.118)
  RCONST(584) = (5.37D-12)
  RCONST(585) = (3.22D-12)
  RCONST(586) = (1.33D-11)
  RCONST(587) = (1.44d-12*EXP(-1862.0/temp)*5.5)
  RCONST(588) = (3.27D-11*0.50)
  RCONST(589) = (3.27D-11*0.50)
  RCONST(590) = (3.25D-11*0.50)
  RCONST(591) = (3.25D-11*0.50)
  RCONST(592) = (5.67D-11)
  RCONST(593) = (1.19D-11)
  RCONST(594) = (1.86D-11)
  RCONST(595) = (1.18D-11)
  RCONST(596) = (2.40d-12*EXP(360.0/temp)*0.843)
  RCONST(597) = (2.40d-12*EXP(360.0/temp)*0.843)
  RCONST(598) = (2.40d-12*EXP(360.0/temp)*0.157)
  RCONST(599) = (2.40d-12*EXP(360.0/temp)*0.157)
  RCONST(600) = (2.50d-12)
  RCONST(601) = (2.50d-12)
  RCONST(602) = (0.520*2.91d-13*EXP(1300.0/temp)*0.890)
  RCONST(603) = (0.520*2.91d-13*EXP(1300.0/temp)*0.890)
  RCONST(604) = (8.80D-13*RO2)
  RCONST(605) = (8.80D-13*RO2)
  RCONST(606) = (9.45D-11)
  RCONST(607) = (1.28D-10)
  RCONST(608) = (5.67D-11)
  RCONST(609) = (2.40d-12*EXP(360.0/temp)*0.833)
  RCONST(610) = (2.40d-12*EXP(360.0/temp)*0.167)
  RCONST(611) = (2.50d-12)
  RCONST(612) = (0.520*2.91d-13*EXP(1300.0/temp)*0.914)
  RCONST(613) = (8.80D-13*RO2)
  RCONST(614) = (9.57D-11)
  RCONST(615) = (1.28D-10)
  RCONST(616) = (7.00D-11*0.55)
  RCONST(617) = (7.00D-11*0.45)
  RCONST(618) = (1.0_dp)
  RCONST(619) = (1.0_dp)
  RCONST(620) = (ARR2(1.12d-11,250._dp,temp))
  RCONST(621) = (iupac_ch3sch3(9.5d-39,5270._dp,7.5d-29,5610._dp,C_M*0.2_dp,temp))
  RCONST(622) = (ARR2(1.9d-13,-520._dp,temp))
  RCONST(623) = (ARR2(4.9d-12,-263._dp,temp))
  RCONST(624) = (1.0d-11)
  RCONST(625) = (ARR2(1.15d-12,-432._dp,temp))
  RCONST(626) = (ARR2(3.0d-11,-210._dp,temp))
  RCONST(627) = (1.2d-11)
  RCONST(628) = (6.0d-13)
  RCONST(629) = (ARR2(5.0d13,9673._dp,temp))
  RCONST(630) = (2.2d-12)
  RCONST(631) = (3.0d-13)
  RCONST(632) = (5.0d-11)
  RCONST(633) = (ARR2(1.36d14,11071._dp,temp))
  RCONST(634) = (8.7d-11)
  RCONST(635) = (9.d-11)
  RCONST(636) = (1.0d-13)
END SUBROUTINE crimech_Update_RCONST

! Begin Rate Law Functions from KPP_HOME/util/UserRateLaws

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  User-defined Rate Law functions

! !! FROM chem/KPP/kpp/kpp-2.1/util/WRFconform/WRFUserRateLaws !!!!
!
!  Note: the default argument type for rate laws, as read from the equations file, is single precision
!        but all the internal calculations are performed in double precision
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Arrhenius (added by psaide 15-07-2009)
    REAL(kind=dp) FUNCTION ARR( A0,B0,C0, TEMP )
      REAL(kind=dp) :: TEMP
      REAL(kind=dp) A0,B0,C0
      ARR =  A0 * EXP( -B0 /TEMP ) * (TEMP/300._dp)**C0
    END FUNCTION ARR

!~~~> Simplified Arrhenius, with two arguments
   REAL(kind=dp) FUNCTION ARR2( A0,B0, TEMP )
      REAL(kind=dp) :: TEMP 
      REAL(kind=dp) A0,B0           
      ARR2 = A0 * EXP( -B0 /TEMP )              
   END FUNCTION ARR2          

! EP2 (added by psaide 22-07-2009)
   REAL(kind=dp) FUNCTION EP2( A0,C0,A2,C2,A3,C3,TEMP,cair)
      REAL(kind=dp) :: TEMP
      REAL(kind=dp) :: cair
      REAL(kind=dp) A0,C0,A2,C2,A3,C3
      REAL(kind=dp) K0,K2,K3

      K0 = A0 * EXP(-C0 /TEMP)
      K2 = A2 * EXP(-C2 /TEMP)
      K3 = A3 * EXP(-C3 /TEMP)
!      K3 = K3 * 2.45E13_dp * 1.0E6_dp
      K3 = K3 * cair
      EP2 = K0 + K3/(1._dp+K3/K2 )
   END FUNCTION EP2

! EP3 (added by psaide 15-07-2009)
    REAL(kind=dp) FUNCTION EP3(A1,C1,A2,C2,TEMP,cair)
      REAL(kind=dp) :: TEMP
      REAL(kind=dp) :: cair
      REAL(kind=dp) A1, C1, A2, C2
      REAL(kind=dp) K1, K2
 
      K1 = A1 * EXP(-C1 /TEMP)
      K2 = A2 * EXP(-C2 /TEMP)
!      EP3 = K1 + K2*(1.0E6_dp * 2.45E13_dp)
      EP3 = K1 + K2*cair
    END FUNCTION EP3

! FALL (added by psaide 15-07-2009)
    REAL(kind=dp) FUNCTION FALL( A0,B0,C0,A1,B1,C1,CF,TEMP,cair)
 
      INTRINSIC LOG10  

      REAL(kind=dp) :: TEMP
      REAL(kind=dp) :: cair
      REAL(kind=dp) A0,B0,C0,A1,B1,C1,CF
      REAL(kind=dp) K0, K1

      K0 = A0 * EXP(-B0 /TEMP)* (TEMP/300._dp)**C0
      K1 = A1 * EXP(-B1 /TEMP)* (TEMP/300._dp)**C1
!      K0 = K0 * 2.45E13_dp * 1.0E6_dp
      K0 = K0 * cair
      K1 = K0/K1
      FALL = (K0/(1._dp+K1))*CF**(1._dp/(1._dp+(LOG10(K1))**2))

    END FUNCTION FALL

! FALL2 (added by psaide 28-10-2009)
    REAL(kind=dp) FUNCTION F2( A0,B0,C0,A1,B1,C1,CF,CN,TEMP,cair)

      INTRINSIC LOG10

      REAL(kind=dp) :: TEMP
      REAL(kind=dp) :: cair
      REAL(kind=dp) A0,B0,C0,A1,B1,C1,CF,CN
      REAL(kind=dp) K0, K1

      K0 = A0 * EXP(-B0 /TEMP)* (TEMP/300._dp)**C0
      K1 = A1 * EXP(-B1 /TEMP)* (TEMP/300._dp)**C1
!      K0 = K0 * 2.45E13_dp * 1.0E6_dp
      K0 = K0 * cair
      K1 = K0/K1
      F2 = (K0/(1._dp+K1))*CF**(1._dp/(1._dp+(LOG10(K1)/CN)**2))

    END FUNCTION F2
                                                                   
!------------------------------------
! Troe reactions (as in Stockwell et al, 1997)

    REAL(kind=dp) FUNCTION TROE(k0_300K,n,kinf_300K,m,temp,cair)

    INTRINSIC LOG10

    REAL(kind=dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(kind=dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(kind=dp),          INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(kind=dp),          INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(kind=dp),          INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(kind=dp),          INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(kind=dp)             :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    TROE   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

   END FUNCTION TROE



!-------------------------------------------
! Troe equilibrium reactions (as in Stockwell et al, 1997)

    REAL(kind=dp) FUNCTION TROEE(A, B, k0_300K,n,kinf_300K,m,temp,cair)

    INTRINSIC LOG10

    REAL(kind=dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(kind=dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(kind=dp),     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(kind=dp),     INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(kind=dp),     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(kind=dp),     INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(kind=dp),     INTENT(IN) :: A, B 
    REAL(kind=dp)             :: zt_help, k0_T, kinf_T, k_ratio, troe
    

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    troe   = k0_T/(1._dp+k_ratio)*0.6_dp**(1._dp/(1._dp+LOG10(k_ratio)**2))

    TROEE = A * EXP( - B / temp) * troe
    
    

  END FUNCTION TROEE

!------------------------
! k=T^2 C exp (-D/T) reactions

   REAL(kind=dp) FUNCTION THERMAL_T2(c, d ,temp)
    REAL(kind=dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(kind=dp),     INTENT(IN) :: c, d


     THERMAL_T2= temp**2._dp * c * EXP(- d  / temp)

   END FUNCTION THERMAL_T2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  End of User-defined Rate Law functions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! End Rate Law Functions from KPP_HOME/util/UserRateLaws


! Begin INLINED Rate Law Functions


!************** SPECIAL RATE FUNCTIONS **********************

REAL(KIND=dp) FUNCTION k46( TEMP, C_M )
    REAL(KIND=dp), INTENT(IN) :: temp, c_m
    REAL(KIND=dp) :: k0, k2, k3

   k0=7.2E-15_dp * EXP(785._dp/TEMP)
   k2=4.1E-16_dp * EXP(1440._dp/TEMP)
   k3=1.9E-33_dp * EXP(725._dp/TEMP)  * C_M

   k46=k0+k3/(1+k3/k2)

!   print*,'k46=',k46
END FUNCTION k46

REAL(KIND=dp) FUNCTION k47( TEMP, C_M )
    REAL(KIND=dp), INTENT(IN) :: temp, c_m
    REAL(KIND=dp) :: k0, ki, fc, x, ssign, f12
    k0 = 3.00d-31*((temp/300.0)**(-3.3))*C_M
    ki = 1.50d-12
    fc = 0.6
    x = 1.0d+0
    ssign = dsign(x,(k0-ki))
    f12 =10**(dlog10(fc)/(1.0+(ssign*(ABS(dlog10(k0/ki)))**(2.0))))
    k47=(k0*ki*f12)/(k0+ki)
!    print*,'k47=',k47
END FUNCTION k47

REAL(KIND=dp) FUNCTION k48( TEMP, C_M )
    REAL(KIND=dp), INTENT(IN) :: temp, c_m
    REAL(KIND=dp) :: k0, ki, fc, x, ssign, f17
    k0 = 5.00d-30*((temp/298.0)**(-1.5))*C_M
    ki = 9.40d-12*EXP(-700.0/temp)
    fc = (EXP(-temp/580.0) + EXP(-2320.0/temp))
    x = 1.0d+0
    ssign = dsign(x,(k0-ki))
    f17=10**(dlog10(fc)/(1.0+(ssign*(ABS(dlog10(k0/ki)))**(2.0))))
    k48=(k0*ki*f17)/(k0+ki)
!    print*,'k48=',k48
END FUNCTION k48

      REAL(KIND=dp) FUNCTION RJPL( K0300, Q, KU300, R, M, T )
      REAL(KIND=dp) :: k0300,q,ku300,r,m,t
      REAL(KIND=dp) :: tt,k0,ku,k0m,kk,lgkk,e,f
! JPL standard three body reaction rate format extended
      TT= T / 3.D2
      K0= K0300 * exp(-1._dp*Q*log(TT))
      KU= KU300 * exp(-1._dp*R*log(TT))
      K0M= K0 * M
      KK= K0M / KU
      LGKK=0.43429448190324926_dp * LOG(KK) ! = log10(KK)
      E=1.D0 / ( 1.D0 + LGKK*LGKK )
      F=exp(-0.5108256237659887_dp*E)       ! -0.51=log(0.6)
      RJPL = F * K0M / ( 1.D0 + KK )
!      print*,'RJPL=',RJPL
      END FUNCTION
!---------------------------------------------------------------------



      REAL(KIND=dp) FUNCTION RALKE( K0300, Q, KU, Fc, M, T )
      REAL(KIND=dp) :: k0300,q,m,t,Fc
      real(KIND=dp) :: tt,k0,ku,k0m,kk,lgkk,e,f
! special function for alkene+OH reactions
      TT= T / 3.D2
      K0= K0300 * exp(-1._dp*Q*log(TT))
      K0M= K0 * M
      KK= K0M / KU
      LGKK=0.43429448190324926_dp * LOG(KK) ! = log10(KK)
      E=1.D0 / ( 1.D0 + LGKK*LGKK )
      F=exp(log(Fc)*E)
      RALKE = F * K0M / ( 1.D0 + KK )
!      print*,'RALKE=',RALKE
      END FUNCTION


	real(kind=dp) function iupac_ch3sch3(a2,b2,a3,b3,cin_o2,temp)
		!rate calculation for CH3SCH3 + OH = CH3SCH3OO + H2O
		! from IUPAC report (www.iupac-kinetic.ch.cam.ac.uk)
		real(kind=dp) :: cin_o2, tr, temp
		real(kind=dp) :: a2, b2, a3, b3
	
		tr = 1._dp + ARR2(a3,b3,temp)*cin_o2
		iupac_ch3sch3 = ARR2(a2,b2,temp)*cin_o2/tr

	end function iupac_ch3sch3
    
!---------------------------------------------------------------------

!- SAN: adding standard 3-body reaction using convention of MCM & IUPAC recommenHNO3
! - Explicit form of TROE reactions
!	Based on Atkinson et. al. 2004

REAL(KIND=dp) FUNCTION KMT_IUPAC(k0_300K,n,kinf_300K,m,Fc,temp,cair)

    INTRINSIC LOG10

   	REAL(KIND=dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(KIND=dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(KIND=dp), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(KIND=dp), INTENT(IN) :: n         ! exponent for low pressure limit
    	!!! n.b. - remember to flip sign of exponents from IUPAC data sheets !!!
    REAL(KIND=dp), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(KIND=dp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(KIND=dp), INTENT(IN) :: Fc        ! Approximate broadening factor
    
    REAL(KIND=dp) :: zt_help, k0_T, kinf_T, k_ratio, Nint, F_exp

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    Nint = 0.75_dp - 1.27_dp*LOG10(Fc)
    ! Calculate explicit broadening factor:
    F_exp = Fc ** (1._dp / (1._dp + ( LOG10(k_ratio) / Nint )**2._dp ) )
        
    KMT_IUPAC   = k0_T/(1._dp+k_ratio) * F_exp

END FUNCTION KMT_IUPAC

!---------------------------------------------------------------------

!- SAN: Function for calculating NO + OH [+ M] 3-body reaction 
!- Explicit form of TROE reaction with temperature dependent Fc

REAL(KIND=dp) FUNCTION KMT_OH_NO(temp,cair)

    INTRINSIC LOG10

   	REAL(KIND=dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(KIND=dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    
    REAL(KIND=dp) :: k0_300K, n, kinf_300K, m, zt_help 
    REAL(KIND=dp) :: k0_T, kinf_T, k_ratio, Nint, Fc, F_exp
    
    k0_300K = 7.4D-31						! low pressure limit at 300 K
    n		= 2.4_dp						! exponent for low pressure limit
    kinf_300K = 3.3D-11						! high pressure limit at 300 K
    m		= 0.3_dp						! exponent for high pressure limit

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    
    ! OH + NO [+ M] uses temperature dependent Fc:
    Fc = exp(-temp / 1420._dp)
    
    Nint = 0.75_dp - 1.27_dp*LOG10(Fc)
    
    ! Calculate explicit broadening factor:
    F_exp = Fc ** (1._dp / (1._dp + ( LOG10(k_ratio) / Nint )**2._dp ) )
    
    KMT_OH_NO = k0_T/(1._dp+k_ratio) * F_exp

END FUNCTION KMT_OH_NO





! End INLINED Rate Law Functions


END MODULE crimech_UpdateRconstWRF

