























MODULE t1_mozcart_UpdateRconstWRF

  USE t1_mozcart_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE t1_mozcart_Update_RCONST(  &


   rh, aer_so4, aer_oc2, aer_bc2, &
   aero_srf_area_diag, &
   sulf_srf_area, oc_srf_area, bc_srf_area, &


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


   integer, intent(in)       :: aero_srf_area_diag
   real(kind=dp), intent(in) :: rh
   real(kind=dp), intent(in) :: aer_so4
   real(kind=dp), intent(in) :: aer_oc2
   real(kind=dp), intent(in) :: aer_bc2
   real, intent(inout) :: sulf_srf_area
   real, intent(inout) :: oc_srf_area
   real, intent(inout) :: bc_srf_area






   real(dp) :: aer_srf_area(3)
   real(dp) :: aer_diam(3)

   call aero_surfarea( aer_srf_area, aer_diam, rh, temp, &
                       aer_so4, aer_oc2, aer_bc2 )

   if( aero_srf_area_diag > 0 ) then
     sulf_srf_area = aer_srf_area(1)
     oc_srf_area   = aer_srf_area(2)
     bc_srf_area   = aer_srf_area(3)
   endif



  RCONST(1) = (j(Pj_h2o2))
  RCONST(2) = (j(Pj_o2))
  RCONST(3) = (j(Pj_o33p))
  RCONST(4) = (j(Pj_o31d))
  RCONST(5) = (j(Pj_hno3))
  RCONST(6) = (j(Pj_hno4))
  RCONST(7) = (j(Pj_n2o))
  RCONST(8) = (j(Pj_n2o5))
  RCONST(9) = (j(Pj_no2))
  RCONST(10) = (j(Pj_no3o))
  RCONST(11) = (j(Pj_ch3o2h))
  RCONST(12) = (j(Pj_ch3o2h))
  RCONST(13) = (j(Pj_ch3o2h))
  RCONST(14) = (0.1_dp*j(Pj_no2))
  RCONST(15) = (0.2_dp*j(Pj_no2))
  RCONST(16) = (0.14_dp*j(Pj_no2))
  RCONST(17) = (0.2_dp*j(Pj_no2))
  RCONST(18) = (0.2_dp*j(Pj_no2))
  RCONST(19) = (0.006_dp*j(Pj_no2))
  RCONST(20) = (j(Pj_ch3o2h))
  RCONST(21) = (j(Pj_ch3o2h))
  RCONST(22) = (j(Pj_ch3o2h))
  RCONST(23) = (j(Pj_ch3o2h))
  RCONST(24) = (j(Pj_ch2or))
  RCONST(25) = (j(Pj_ch2om))
  RCONST(26) = (j(Pj_ch3cho))
  RCONST(27) = (j(Pj_ch3coch3))
  RCONST(28) = (j(Pj_ch3cocho))
  RCONST(29) = (0.28_dp*j(Pj_h2o2))
  RCONST(30) = (j(Pj_ch3o2h))
  RCONST(31) = (j(Pj_ch3o2h))
  RCONST(32) = (j(Pj_glyald))
  RCONST(33) = (j(Pj_gly))
  RCONST(34) = (j(Pj_glyald))
  RCONST(35) = (j(Pj_ch2or))
  RCONST(36) = (0.006_dp*j(Pj_no2))
  RCONST(37) = (j(Pj_hyac))
  RCONST(38) = (j(Pj_ch3o2h))
  RCONST(39) = (j(Pj_ch3o2h))
  RCONST(40) = (j(Pj_macr))
  RCONST(41) = (j(Pj_ch3o2h))
  RCONST(42) = (j(Pj_ch3coch3))
  RCONST(43) = (j(Pj_ch3o2h))
  RCONST(44) = (j(Pj_pan))
  RCONST(45) = (j(Pj_mvk))
  RCONST(46) = (j(Pj_ch2or))
  RCONST(47) = (j(Pj_ch2or))
  RCONST(48) = (j(Pj_ch3o2h))
  RCONST(49) = (j(Pj_ch3cho))
  RCONST(50) = (j(Pj_pan))
  RCONST(51) = (j(Pj_ch3o2h))
  RCONST(52) = (j(Pj_ch3o2h))
  RCONST(53) = (j(Pj_ch3o2h))
  RCONST(54) = (0.1_dp*j(Pj_no2))
  RCONST(55) = (j(Pj_ch3o2h))
  RCONST(56) = (j(Pj_ch3o2h))
  RCONST(57) = (j(Pj_ch3o2h))
  RCONST(58) = (j(Pj_ch3cho))
  RCONST(59) = (j(Pj_ch3cho))
  RCONST(60) = (j(Pj_ch3o2h))
  RCONST(61) = (j(Pj_ch3o2h))
  RCONST(62) = (j(Pj_ch3o2h))
  RCONST(63) = (j(Pj_ch3o2h))
  RCONST(64) = (1.500000e-10_dp)
  RCONST(65) = (1.200000e-10_dp)
  RCONST(66) = (ARR2(1.630000e-10_dp,-60.00_dp,TEMP))
  RCONST(67) = (ARR2(2.150000e-11_dp,-110.00_dp,TEMP))
  RCONST(68) = (ARR2(3.300000e-11_dp,-55.00_dp,TEMP))
  RCONST(69) = (1.200000e-10_dp)
  RCONST(70) = (ARR2(8.000000e-12_dp,2060.00_dp,TEMP))
  RCONST(71) = (usr_O_O_M(temp))
  RCONST(72) = (usr_O_O2(temp))
  RCONST(73) = (ARR2(1.600000e-11_dp,4570.00_dp,TEMP))
  RCONST(74) = (ARR2(1.400000e-12_dp,2000.00_dp,TEMP))
  RCONST(75) = (ARR2(3.000000e-11_dp,-200.00_dp,TEMP))
  RCONST(76) = (ARR2(1.000000e-14_dp,490.00_dp,TEMP))
  RCONST(77) = (ARR2(2.800000e-12_dp,1800.00_dp,TEMP))
  RCONST(78) = (1.800000e-12_dp)
  RCONST(79) = (ARR2(4.800000e-11_dp,-250.00_dp,TEMP))
  RCONST(80) = (ARR2(1.800000e-11_dp,-180.00_dp,TEMP))
  RCONST(81) = (ARR2(1.700000e-12_dp,940.00_dp,TEMP))
  RCONST(82) = (1.800000e-12_dp)
  RCONST(83) = (JPL_TROE(6.900000e-31_dp,1.00_dp,2.600000e-11_dp,0.00_dp,0.60_dp,TEMP,C_M))
  RCONST(84) = (usr_HO2_HO2(temp,c_m,c_h2o))
  RCONST(85) = (ARR2(1.300000e-12_dp,-380.00_dp,TEMP))
  RCONST(86) = (ARR2(5.100000e-12_dp,-210.00_dp,TEMP))
  RCONST(87) = (ARR2(1.200000e-13_dp,2450.00_dp,TEMP))
  RCONST(88) = (JPL_TROE(2.500000e-31_dp,1.80_dp,2.200000e-11_dp,0.70_dp,0.60_dp,TEMP,C_M))
  RCONST(89) = (3.500000e-12_dp)
  RCONST(90) = (ARR2(1.500000e-11_dp,-170.00_dp,TEMP))
  RCONST(91) = (1.000000e-11_dp)
  RCONST(92) = (2.200000e-11_dp)
  RCONST(93) = (ARR2(3.300000e-12_dp,-270.00_dp,TEMP))
  RCONST(94) = (ARR2(3.000000e-12_dp,1500.00_dp,TEMP))
  RCONST(95) = (JPL_TROE(9.000000e-32_dp,1.50_dp,3.000000e-11_dp,0.00_dp,0.60_dp,TEMP,C_M))
  RCONST(96) = (ARR2(7.260000e-11_dp,-20.00_dp,TEMP))
  RCONST(97) = (ARR2(4.640000e-11_dp,-20.00_dp,TEMP))
  RCONST(98) = (JPL_TROE(1.900000e-31_dp,3.40_dp,4.000000e-12_dp,0.30_dp,0.60_dp,TEMP,C_M))
  RCONST(99) = (JPL_TROE(2.400000e-30_dp,3.00_dp,1.600000e-12_dp,-0.10_dp,0.60_dp,TEMP,C_M))
  RCONST(100) = (JPL_TROE(1.800000e-30_dp,3.00_dp,2.800000e-11_dp,0.00_dp,0.60_dp,TEMP,C_M))
  RCONST(101) = (usr_HNO3_OH(temp,c_m))
  RCONST(102) = (usr_HO2NO2_M(temp,c_m))
  RCONST(103) = (usr_N2O5_M(temp,c_m))
  RCONST(104) = (ARR2(9.700000e-15_dp,-625.00_dp,TEMP))
  RCONST(105) = (ARR2(6.000000e-13_dp,2058.00_dp,TEMP))
  RCONST(106) = (ARR2(3.400000e-11_dp,1600.00_dp,TEMP))
  RCONST(107) = (ARR2(5.500000e-12_dp,-125.00_dp,TEMP))
  RCONST(108) = (ARR2(5.000000e-13_dp,424.00_dp,TEMP))
  RCONST(109) = (ARR2(1.900000e-14_dp,-706.00_dp,TEMP))
  RCONST(110) = (ARR2(4.100000e-13_dp,-750.00_dp,TEMP))
  RCONST(111) = (ARR2(2.800000e-12_dp,-300.00_dp,TEMP))
  RCONST(112) = (ARR2(2.900000e-12_dp,345.00_dp,TEMP))
  RCONST(113) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(114) = (ARR2(2.450000e-12_dp,1775.00_dp,TEMP))
  RCONST(115) = (JPL_TROE(4.280000e-33_dp,0.00_dp,9.300000e-15_dp,-4.42_dp,0.80_dp,TEMP,C_M))
  RCONST(116) = (4.000000e-13_dp)
  RCONST(117) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(118) = (ARR2(2.400000e+12_dp,7000.00_dp,TEMP))
  RCONST(119) = (ARR2(2.600000e-12_dp,-265.00_dp,TEMP))
  RCONST(120) = (ARR2(1.080000e-10_dp,-105.00_dp,TEMP))
  RCONST(121) = (usr_CO_OH_a(temp,c_m))
  RCONST(122) = (JPL_TROE(5.500000e-30_dp,0.00_dp,8.300000e-13_dp,-2.00_dp,0.60_dp,TEMP,C_M))
  RCONST(123) = (ARR2(1.200000e-14_dp,2630.00_dp,TEMP))
  RCONST(124) = (6.800000e-14_dp)
  RCONST(125) = (2.000000e-13_dp)
  RCONST(126) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(127) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(128) = (ARR2(6.900000e-12_dp,230.00_dp,TEMP))
  RCONST(129) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(130) = (ARR2(7.660000e-12_dp,1020.00_dp,TEMP))
  RCONST(131) = (ARR2(1.400000e-12_dp,1900.00_dp,TEMP))
  RCONST(132) = (ARR2(4.630000e-12_dp,-350.00_dp,TEMP))
  RCONST(133) = (ARR2(7.800000e-13_dp,1050.00_dp,TEMP))
  RCONST(134) = (ARR2(2.900000e-12_dp,-500.00_dp,TEMP))
  RCONST(135) = (ARR2(2.000000e-12_dp,-500.00_dp,TEMP))
  RCONST(136) = (ARR2(4.300000e-13_dp,-1040.00_dp,TEMP))
  RCONST(137) = (ARR2(8.100000e-12_dp,-270.00_dp,TEMP))
  RCONST(138) = (7.000000e-13_dp)
  RCONST(139) = (1.000000e-12_dp)
  RCONST(140) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(141) = (ARR2(4.200000e-12_dp,-180.00_dp,TEMP))
  RCONST(142) = (ARR2(1.600000e+11_dp,4150.00_dp,TEMP))
  RCONST(143) = (1.000000e-14_dp)
  RCONST(144) = (1.000000e-11_dp)
  RCONST(145) = (1.150000e-11_dp)
  RCONST(146) = (4.000000e-14_dp)
  RCONST(147) = (JPL_TROE(8.600000e-29_dp,3.10_dp,9.000000e-12_dp,0.85_dp,0.48_dp,TEMP,C_M))
  RCONST(148) = (JPL_TROE(9.700000e-29_dp,5.60_dp,9.300000e-12_dp,1.50_dp,0.60_dp,TEMP,C_M))
  RCONST(149) = (usr_PAN_M(temp,c_m))
  RCONST(150) = (ARR2(4.600000e-13_dp,1156.00_dp,TEMP))
  RCONST(151) = (ARR2(6.500000e-15_dp,1900.00_dp,TEMP))
  RCONST(152) = (ARR2(3.750000e-13_dp,40.00_dp,TEMP))
  RCONST(153) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(154) = (ARR2(4.200000e-12_dp,-180.00_dp,TEMP))
  RCONST(155) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(156) = (ARR2(8.700000e-12_dp,615.00_dp,TEMP))
  RCONST(157) = (ARR2(1.400000e-12_dp,1860.00_dp,TEMP))
  RCONST(158) = (ARR2(8.400000e-13_dp,-830.00_dp,TEMP))
  RCONST(159) = (3.000000e-12_dp)
  RCONST(160) = (6.700000e-13_dp)
  RCONST(161) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(162) = (ARR2(4.200000e-12_dp,-180.00_dp,TEMP))
  RCONST(163) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(164) = (ARR2(7.100000e-13_dp,-500.00_dp,TEMP))
  RCONST(165) = (ARR2(8.600000e-13_dp,-700.00_dp,TEMP))
  RCONST(166) = (ARR2(2.900000e-12_dp,-300.00_dp,TEMP))
  RCONST(167) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(168) = (JPL_TROE(8.000000e-27_dp,3.50_dp,3.000000e-11_dp,0.00_dp,0.50_dp,TEMP,C_M))
  RCONST(169) = (usr_CH3COCH3_OH(temp))
  RCONST(170) = (3.500000e-13_dp)
  RCONST(171) = (5.400000e-11_dp)
  RCONST(172) = (ARR2(4.800000e-12_dp,-120.00_dp,TEMP))
  RCONST(173) = (ARR2(5.100000e-14_dp,-693.00_dp,TEMP))
  RCONST(174) = (2.000000e-12_dp)
  RCONST(175) = (1.400000e-11_dp)
  RCONST(176) = (ARR2(5.000000e-13_dp,-400.00_dp,TEMP))
  RCONST(177) = (ARR2(8.000000e-13_dp,-700.00_dp,TEMP))
  RCONST(178) = (2.400000e-12_dp)
  RCONST(179) = (ARR2(2.700000e-12_dp,-360.00_dp,TEMP))
  RCONST(180) = (ARR2(1.300000e-13_dp,-360.00_dp,TEMP))
  RCONST(181) = (ARR2(1.500000e-15_dp,2100.00_dp,TEMP))
  RCONST(182) = (ARR2(9.600000e-12_dp,-360.00_dp,TEMP))
  RCONST(183) = (ARR2(2.300000e-11_dp,-200.00_dp,TEMP))
  RCONST(184) = (ARR2(4.600000e-12_dp,-530.00_dp,TEMP))
  RCONST(185) = (ARR2(2.000000e-12_dp,-500.00_dp,TEMP))
  RCONST(186) = (ARR2(4.300000e-13_dp,-1040.00_dp,TEMP))
  RCONST(187) = (ARR2(2.300000e-12_dp,-530.00_dp,TEMP))
  RCONST(188) = (ARR2(5.300000e-12_dp,-360.00_dp,TEMP))
  RCONST(189) = (5.000000e-12_dp)
  RCONST(190) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(191) = (ARR2(4.200000e-12_dp,-180.00_dp,TEMP))
  RCONST(192) = (ARR2(2.300000e-12_dp,170.00_dp,TEMP))
  RCONST(193) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(194) = (JPL_TROE(8.000000e-27_dp,3.50_dp,3.000000e-11_dp,0.00_dp,0.50_dp,TEMP,C_M))
  RCONST(195) = (ARR2(8.500000e-16_dp,1520.00_dp,TEMP))
  RCONST(196) = (ARR2(4.130000e-12_dp,-452.00_dp,TEMP))
  RCONST(197) = (usr_MCO3_NO2(temp,c_m))
  RCONST(198) = (usr_MPAN_M(temp,c_m))
  RCONST(199) = (1.600000e-12_dp)
  RCONST(200) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(201) = (6.700000e-12_dp)
  RCONST(202) = (ARR2(5.400000e-14_dp,-870.00_dp,TEMP))
  RCONST(203) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(204) = (3.500000e-12_dp)
  RCONST(205) = (ARR2(1.860000e-11_dp,-175.00_dp,TEMP))
  RCONST(206) = (ARR2(1.860000e-11_dp,-175.00_dp,TEMP))
  RCONST(207) = (1.300000e-11_dp)
  RCONST(208) = (1.400000e-11_dp)
  RCONST(209) = (ARR2(5.000000e-13_dp,-400.00_dp,TEMP))
  RCONST(210) = (ARR2(8.000000e-13_dp,-700.00_dp,TEMP))
  RCONST(211) = (ARR2(4.400000e-12_dp,-180.00_dp,TEMP))
  RCONST(212) = (2.400000e-12_dp)
  RCONST(213) = (1.400000e-11_dp)
  RCONST(214) = (ARR2(5.000000e-13_dp,-400.00_dp,TEMP))
  RCONST(215) = (ARR2(8.000000e-13_dp,-700.00_dp,TEMP))
  RCONST(216) = (ARR2(1.600000e+9_dp,8300.00_dp,TEMP))
  RCONST(217) = (ARR2(4.400000e-12_dp,-180.00_dp,TEMP))
  RCONST(218) = (2.400000e-12_dp)
  RCONST(219) = (4.000000e-11_dp)
  RCONST(220) = (4.000000e-11_dp)
  RCONST(221) = (ARR2(3.030000e-12_dp,446.00_dp,TEMP))
  RCONST(222) = (1.400000e-11_dp)
  RCONST(223) = (ARR2(5.000000e-13_dp,-400.00_dp,TEMP))
  RCONST(224) = (ARR2(8.000000e-13_dp,-700.00_dp,TEMP))
  RCONST(225) = (ARR2(2.700000e-12_dp,-360.00_dp,TEMP))
  RCONST(226) = (2.400000e-12_dp)
  RCONST(227) = (4.000000e-11_dp)
  RCONST(228) = (ARR2(1.050000e-14_dp,2000.00_dp,TEMP))
  RCONST(229) = (ARR2(2.540000e-11_dp,-410.00_dp,TEMP))
  RCONST(230) = (ARR2(1.520000e-11_dp,-200.00_dp,TEMP))
  RCONST(231) = (7.000000e-11_dp)
  RCONST(232) = (1.000000e-10_dp)
  RCONST(233) = (ARR2(1.300000e-12_dp,-640.00_dp,TEMP))
  RCONST(234) = (ARR2(5.000000e-13_dp,-400.00_dp,TEMP))
  RCONST(235) = (ARR2(8.000000e-13_dp,-700.00_dp,TEMP))
  RCONST(236) = (ARR2(2.700000e-12_dp,-360.00_dp,TEMP))
  RCONST(237) = (2.400000e-12_dp)
  RCONST(238) = (ARR2(1.520000e-12_dp,-200.00_dp,TEMP))
  RCONST(239) = (ARR2(4.300000e-13_dp,-1040.00_dp,TEMP))
  RCONST(240) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(241) = (1.400000e-11_dp)
  RCONST(242) = (ARR2(4.600000e-14_dp,400.00_dp,TEMP))
  RCONST(243) = (ARR2(4.300000e-13_dp,-1040.00_dp,TEMP))
  RCONST(244) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(245) = (2.400000e-12_dp)
  RCONST(246) = (ARR2(3.750000e-13_dp,40.00_dp,TEMP))
  RCONST(247) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(248) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(249) = (1.000000e-17_dp)
  RCONST(250) = (ARR2(8.100000e-12_dp,-610.00_dp,TEMP))
  RCONST(251) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(252) = (ARR2(4.300000e-13_dp,-1040.00_dp,TEMP))
  RCONST(253) = (ARR2(7.500000e-12_dp,-290.00_dp,TEMP))
  RCONST(254) = (ARR2(2.300000e-12_dp,193.00_dp,TEMP))
  RCONST(255) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(256) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(257) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(258) = (ARR2(5.900000e-12_dp,-225.00_dp,TEMP))
  RCONST(259) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(260) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(261) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(262) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(263) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(264) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(265) = (4.700000e-11_dp)
  RCONST(266) = (ARR2(4.300000e-13_dp,-1040.00_dp,TEMP))
  RCONST(267) = (ARR2(7.500000e-12_dp,-290.00_dp,TEMP))
  RCONST(268) = (JPL_TROE(9.700000e-29_dp,5.60_dp,9.300000e-12_dp,1.50_dp,0.60_dp,TEMP,C_M))
  RCONST(269) = (ARR2(4.300000e-13_dp,-1040.00_dp,TEMP))
  RCONST(270) = (ARR2(7.500000e-12_dp,-290.00_dp,TEMP))
  RCONST(271) = (JPL_TROE(9.700000e-29_dp,5.60_dp,9.300000e-12_dp,1.50_dp,0.60_dp,TEMP,C_M))
  RCONST(272) = (ARR2(4.300000e-13_dp,-1040.00_dp,TEMP))
  RCONST(273) = (ARR2(7.500000e-12_dp,-290.00_dp,TEMP))
  RCONST(274) = (JPL_TROE(9.700000e-29_dp,5.60_dp,9.300000e-12_dp,1.50_dp,0.60_dp,TEMP,C_M))
  RCONST(275) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(276) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(277) = (ARR2(4.700000e-13_dp,-1220.00_dp,TEMP))
  RCONST(278) = (2.100000e-12_dp)
  RCONST(279) = (2.800000e-13_dp)
  RCONST(280) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(281) = (JPL_TROE(9.700000e-29_dp,5.60_dp,9.300000e-12_dp,1.50_dp,0.60_dp,TEMP,C_M))
  RCONST(282) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(283) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(284) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(285) = (ARR2(1.700000e-12_dp,-352.00_dp,TEMP))
  RCONST(286) = (usr_PBZNIT_M(temp,c_m))
  RCONST(287) = (1.700000e-11_dp)
  RCONST(288) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(289) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(290) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(291) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(292) = (ARR2(2.600000e-12_dp,-365.00_dp,TEMP))
  RCONST(293) = (8.400000e-11_dp)
  RCONST(294) = (ARR2(3.800000e-12_dp,-200.00_dp,TEMP))
  RCONST(295) = (ARR2(1.200000e-12_dp,-490.00_dp,TEMP))
  RCONST(296) = (ARR2(6.300000e-16_dp,580.00_dp,TEMP))
  RCONST(297) = (ARR2(1.200000e-11_dp,-440.00_dp,TEMP))
  RCONST(298) = (1.900000e-11_dp)
  RCONST(299) = (1.200000e-14_dp)
  RCONST(300) = (2.000000e-10_dp)
  RCONST(301) = (2.500000e-12_dp)
  RCONST(302) = (ARR2(1.700000e-15_dp,1300.00_dp,TEMP))
  RCONST(303) = (ARR2(1.600000e-11_dp,-470.00_dp,TEMP))
  RCONST(304) = (1.100000e-11_dp)
  RCONST(305) = (ARR2(3.000000e-15_dp,780.00_dp,TEMP))
  RCONST(306) = (ARR2(4.200000e-11_dp,-400.00_dp,TEMP))
  RCONST(307) = (1.200000e-11_dp)
  RCONST(308) = (4.700000e-16_dp)
  RCONST(309) = (2.100000e-10_dp)
  RCONST(310) = (ARR2(2.000000e-12_dp,-500.00_dp,TEMP))
  RCONST(311) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(312) = (ARR2(4.200000e-12_dp,-180.00_dp,TEMP))
  RCONST(313) = (2.400000e-12_dp)
  RCONST(314) = (2.000000e-11_dp)
  RCONST(315) = (ARR2(2.000000e-12_dp,-500.00_dp,TEMP))
  RCONST(316) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(317) = (ARR2(4.200000e-12_dp,-180.00_dp,TEMP))
  RCONST(318) = (2.300000e-11_dp)
  RCONST(319) = (2.000000e-11_dp)
  RCONST(320) = (ARR2(2.000000e-12_dp,-500.00_dp,TEMP))
  RCONST(321) = (ARR2(7.500000e-13_dp,-700.00_dp,TEMP))
  RCONST(322) = (ARR2(4.200000e-12_dp,-180.00_dp,TEMP))
  RCONST(323) = (3.300000e-11_dp)
  RCONST(324) = (1.000000e-12_dp)
  RCONST(325) = (5.700000e-11_dp)
  RCONST(326) = (3.400000e-11_dp)
  RCONST(327) = (ARR2(1.900000e-13_dp,-520.00_dp,TEMP))
  RCONST(328) = (ARR2(9.600000e-12_dp,234.00_dp,TEMP))
  RCONST(329) = (ARR2(1.700000e-12_dp,710.00_dp,TEMP))
  RCONST(330) = (usr_DMS_OH(temp,c_m))
  RCONST(331) = (usr_HO2_aer(aer_srf_area,aer_diam,temp))
  RCONST(332) = (usr_N2O5_aer(aer_srf_area,aer_diam,temp))
  RCONST(333) = (usr_NO2_aer(aer_srf_area,aer_diam,temp))
  RCONST(334) = (usr_NO3_aer(aer_srf_area,aer_diam,temp))
  RCONST(335) = (usr_HONITR_aer(aer_srf_area,aer_diam,temp))
  RCONST(336) = (usr_ISOPNIT_aer(aer_srf_area,aer_diam,temp))
  RCONST(337) = (usr_ISOPNIT_aer(aer_srf_area,aer_diam,temp))
  RCONST(338) = (usr_NC4CH2OH_aer(aer_srf_area,aer_diam,temp))
  RCONST(339) = (usr_NC4CHO_aer(aer_srf_area,aer_diam,temp))
  RCONST(340) = (usr_NTERPOOH_aer(aer_srf_area,aer_diam,temp))
  RCONST(341) = (usr_ONITR_aer(aer_srf_area,aer_diam,temp))
  RCONST(342) = (usr_TERPNIT_aer(aer_srf_area,aer_diam,temp))
  RCONST(343) = (usr_GLYOXAL_aer(aer_srf_area,temp))
  RCONST(344) = (usr_SO2_OH(temp,c_m))
END SUBROUTINE t1_mozcart_Update_RCONST













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

REAL(KIND=dp) FUNCTION usr_O_O2( temp )



    REAL(KIND=dp), INTENT(IN) :: temp

    usr_O_O2 = 6.00e-34_dp*(temp/300._dp)**(-2.4_dp)

END FUNCTION usr_O_O2

REAL(KIND=dp) FUNCTION usr_O_O_M( temp )



    REAL(KIND=dp), INTENT(IN) :: temp

    usr_O_O_M = 2.76e-34_dp * exp( 720.0_dp/temp)

END FUNCTION usr_O_O_M

REAL(KIND=dp) FUNCTION usr_HO2_HO2( temp, c_m, c_h2o )




    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m
    REAL(KIND=dp), INTENT(IN) :: c_h2o

    REAL(KIND=dp) :: ko, kinf, fc

    if( c_h2o > 0._dp ) then
       ko   = 2.3e-13_dp * exp( 600._dp/temp )
       kinf = 1.7e-33_dp * c_m * exp( 1000._dp/temp )
       fc = 1._dp + 1.4e-21_dp *c_h2o* exp( 2200._dp/temp )
       usr_HO2_HO2 = (ko + kinf) * fc
    else
       usr_HO2_HO2 = 0._dp
    end if
END FUNCTION usr_HO2_HO2

REAL(KIND=dp) FUNCTION usr_HNO3_OH( temp, c_m )



    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    REAL(KIND=dp) :: k0, k2

   k0 = c_m * 6.5e-34_dp * exp( 1335._dp/temp )
   k2 = exp( 2199._dp/temp )
   k0 = k0 /(1.0_dp + k0/(2.7e-17_dp*k2))
   k2 = exp( 460._dp/temp )

   usr_HNO3_OH = k0 + 2.4e-14_dp * k2

END FUNCTION usr_HNO3_OH


REAL(KIND=dp) FUNCTION usr_HO2NO2_M( temp, c_m )





    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    usr_HO2NO2_M = TROEE( 4.76e26_dp,10900._dp, 1.8e-31_dp , 3.2_dp , &
                          4.7e-12_dp , 1.4_dp , TEMP, C_M ) /c_m

END FUNCTION usr_HO2NO2_M

REAL(KIND=dp) FUNCTION usr_PAN_M( temp, c_m )




    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    usr_PAN_M = TROEE(1.111e28_dp, 14000._dp, 8.5e-29_dp, 6.5_dp, &
                       1.1e-11_dp, 0._dp, TEMP, C_M) /c_m

END FUNCTION usr_PAN_M

REAL(KIND=dp) FUNCTION usr_CH3COCH3_OH( temp )



    REAL(KIND=dp), INTENT(IN) :: temp

    usr_CH3COCH3_OH = 3.82e-11_dp*exp( -2000._dp/temp ) + 1.33e-13_dp

END FUNCTION usr_CH3COCH3_OH

REAL(KIND=dp) FUNCTION usr_MCO3_NO2( temp, c_m )



    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    usr_MCO3_NO2 = 1.1e-11_dp*300._dp/(temp*c_m)

END FUNCTION usr_MCO3_NO2


REAL(KIND=dp) FUNCTION usr_MPAN_M( temp, c_m )



    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    usr_MPAN_M = 1.2221e17_dp*300._dp*exp( -14000._dp/temp )/(temp*c_m)

END FUNCTION usr_MPAN_M

REAL(KIND=dp) FUNCTION usr_N2O5_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_n2o5 = .1_dp
    REAL(KIND=dp) :: c_n2o5, term

    n = size( aero_srf_area )

    c_n2o5 = 1.40e3_dp * sqrt( temp )
    term = 4._dp/(c_n2o5*gamma_n2o5)
    
    usr_N2O5_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_N2O5_aer

REAL(KIND=dp) FUNCTION usr_HONITR_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_honitr = .005_dp
    REAL(KIND=dp) :: c_honitr, term

    n = size( aero_srf_area )

    c_honitr = 1.26e3_dp * sqrt( temp )
    term = 4._dp/(c_honitr*gamma_honitr)
    
    usr_HONITR_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_HONITR_aer

REAL(KIND=dp) FUNCTION usr_ONITR_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_onitr = .005_dp
    REAL(KIND=dp) :: c_onitr, term

    n = size( aero_srf_area )

    c_onitr = 1.20e3_dp * sqrt( temp )
    term = 4._dp/(c_onitr*gamma_onitr)
    
    usr_ONITR_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_ONITR_aer

REAL(KIND=dp) FUNCTION usr_ISOPNIT_aer( aero_srf_area, aero_diam, temp )



    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_isopnit = .005_dp
    REAL(KIND=dp) :: c_isopnit, term

    n = size( aero_srf_area )

    c_isopnit = 1.20e3_dp * sqrt( temp )
    term = 4._dp/(c_isopnit*gamma_isopnit)
    
    usr_ISOPNIT_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_ISOPNIT_aer

REAL(KIND=dp) FUNCTION usr_TERPNIT_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_terpnit = .01_dp
    REAL(KIND=dp) :: c_terpnit, term

    n = size( aero_srf_area )

    c_terpnit = .992e3_dp * sqrt( temp )
    term = 4._dp/(c_terpnit*gamma_terpnit)
    
    usr_TERPNIT_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_TERPNIT_aer

REAL(KIND=dp) FUNCTION usr_NC4CH2OH_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_nc4ch2oh = .005_dp
    REAL(KIND=dp) :: c_nc4ch2oh, term

    n = size( aero_srf_area )

    c_nc4ch2oh = 1.20e3_dp * sqrt( temp )
    term = 4._dp/(c_nc4ch2oh*gamma_nc4ch2oh)
    
    usr_NC4CH2OH_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_NC4CH2OH_aer

REAL(KIND=dp) FUNCTION usr_NC4CHO_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_nc4cho = .005_dp
    REAL(KIND=dp) :: c_nc4cho, term

    n = size( aero_srf_area )

    c_nc4cho = 1.21e3_dp * sqrt( temp )
    term = 4._dp/(c_nc4cho*gamma_nc4cho)
    
    usr_NC4CHO_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_NC4CHO_aer

REAL(KIND=dp) FUNCTION usr_NTERPOOH_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_nterpooh = .01_dp
    REAL(KIND=dp) :: c_nterpooh, term

    n = size( aero_srf_area )

    c_nterpooh = .957e3_dp * sqrt( temp )
    term = 4._dp/(c_nterpooh*gamma_nterpooh)
    
    usr_NTERPOOH_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_NTERPOOH_aer

REAL(KIND=dp) FUNCTION usr_GLYOXAL_aer( aero_srf_area, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: gamma_glyoxal = .0002_dp
    REAL(KIND=dp) :: c_glyoxal, term

    n = size( aero_srf_area )

    c_glyoxal = 1.455e4_dp * sqrt( temp/58._dp )
    term = .25_dp * c_glyoxal * gamma_glyoxal
    
    usr_GLYOXAL_aer = sum( aero_srf_area(1:n) ) * term

END FUNCTION usr_GLYOXAL_aer

REAL(KIND=dp) FUNCTION usr_NO3_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_no3 = 1.e-3_dp
    REAL(KIND=dp) :: c_no3, term

    n = size( aero_srf_area )

    c_no3 = 1.85e3_dp * sqrt( temp )
    term = 4._dp/(c_no3*gamma_no3)
    
    usr_NO3_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_NO3_aer

REAL(KIND=dp) FUNCTION usr_NO2_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_no2 = 1.e-4_dp
    REAL(KIND=dp) :: c_no2, term

    n = size( aero_srf_area )

    c_no2 = 2.15e3_dp * sqrt( temp )
    term = 4._dp/(c_no2*gamma_no2)
    
    usr_NO2_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_NO2_aer

REAL(KIND=dp) FUNCTION usr_DMS_OH( temp, c_m )



    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    REAL(KIND=dp) :: ko, wrk

    wrk   = .21_dp*c_m
    ko    = 1._dp + 5.5e-31_dp*exp( 7460._dp/temp )*wrk
    usr_DMS_OH = 1.7e-42_dp*exp( 7810._dp/temp )*wrk/ko

END FUNCTION usr_DMS_OH


REAL(KIND=dp) FUNCTION usr_HO2_aer( aero_srf_area, aero_diam, temp )


    REAL(KIND=dp), INTENT(IN) :: aero_srf_area(:)         
    REAL(KIND=dp), INTENT(IN) :: aero_diam(:)             
    REAL(KIND=dp), INTENT(IN) :: temp                     

    INTEGER :: n
    REAL(KIND=dp), parameter :: dg = .1_dp
    REAL(KIND=dp), parameter :: gamma_ho2 = .2_dp
    REAL(KIND=dp) :: c_ho2, term

    n = size( aero_srf_area )

    c_ho2 = 2.53e3_dp * sqrt( temp )
    term = 4._dp/(c_ho2*gamma_ho2)
    
    usr_HO2_aer = &
     sum( aero_srf_area(1:n)/(.5_dp*aero_diam(1:n)/dg + term) )

END FUNCTION usr_HO2_aer

REAL(KIND=dp) FUNCTION usr_PBZNIT_M( temp, c_m )




    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    usr_PBZNIT_M = TROEE( 1.111e28_dp,14000._dp, 9.7e-29_dp , 5.6_dp , &
                          9.3e-12_dp , 0._dp , TEMP, C_M) /c_m

END FUNCTION usr_PBZNIT_M

REAL(KIND=dp) FUNCTION usr_N2O5_M( temp, c_m )





    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    usr_N2O5_M = TROEE(3.333e26_dp, 10900._dp, 2.2e-30_dp, 4.4_dp, &
                       1.4e-12_dp , .7_dp , TEMP, C_M ) /c_m

END FUNCTION usr_N2O5_M

REAL(KIND=dp) FUNCTION usr_SO2_OH( temp, c_m )



    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    REAL(KIND=dp) :: fc, k0
    REAL(KIND=dp) :: wrk

    fc    = 3.e-11_dp * (300._dp/temp) ** 3.3_dp
    wrk   = fc * c_m
    k0    = wrk / (1._dp + wrk/1.5e-12_dp)
    usr_SO2_OH = k0 * .6_dp ** (1._dp/(1._dp + &
                 (log10( wrk/1.5e-12_dp ))**2._dp))

END FUNCTION usr_SO2_OH

REAL(KIND=dp) FUNCTION usr_CO_OH_a( temp, c_m )



    REAL(KIND=dp), INTENT(IN) :: temp
    REAL(KIND=dp), INTENT(IN) :: c_m

    REAL(KIND=dp), parameter :: boltz = 1.38044e-16_dp

    usr_CO_OH_a = 1.5e-13_dp * (1._dp + 6.e-7_dp*boltz*c_m*temp)

END FUNCTION usr_CO_OH_a




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




    REAL(KIND=dp) FUNCTION vbs_yield ( nume, den, voc_idx, bin_idx )
    REAL(KIND=dp), INTENT(IN) :: nume, den
    INTEGER, INTENT(IN)       :: voc_idx, bin_idx
    INTEGER, PARAMETER        :: vbs_nbin = 4, vbs_nspec = 9

    
    REAL(KIND=dp)             :: vbs_alphlowN(vbs_nbin,vbs_nspec)
    REAL(KIND=dp)             :: vbs_alphhiN(vbs_nbin,vbs_nspec)
    REAL(KIND=dp)             :: vbs_mw_prec(vbs_nspec)
    
    REAL(KIND=dp), PARAMETER  :: dens_aer = 1.5
    
    REAL(KIND=dp), PARAMETER  :: mw_aer   = 250.0

    
    
    

    
    DATA vbs_alphlowN /   &
    0.0000, 0.0750, 0.0000, 0.0000,   & 
    0.0000, 0.3000, 0.0000, 0.0000,   & 
    0.0045, 0.0090, 0.0600, 0.2250,   & 
    0.0225, 0.0435, 0.1290, 0.3750,   & 
    0.0750, 0.2250, 0.3750, 0.5250,   & 
    0.0750, 0.3000, 0.3750, 0.5250,   & 
    0.0090, 0.0300, 0.0150, 0.0000,   & 
    0.0750, 0.1500, 0.7500, 0.9000,   & 
    0.1073, 0.0918, 0.3587, 0.6075/     

    
    DATA vbs_alphhiN /    &
    0.0000, 0.0375, 0.0000, 0.0000,   & 
    0.0000, 0.1500, 0.0000, 0.0000,   & 
    0.0008, 0.0045, 0.0375, 0.1500,   & 
    0.0030, 0.0255, 0.0825, 0.2700,   & 
    0.0030, 0.1650, 0.3000, 0.4350,   & 
    0.0015, 0.1950, 0.3000, 0.4350,   & 
    0.0003, 0.0225, 0.0150, 0.0000,   & 
    0.0750, 0.1500, 0.7500, 0.9000,   & 
    0.0120, 0.1215, 0.2010, 0.5070/     

    DATA vbs_mw_prec /    &
    120.0, & 
    150.0, & 
    120.0, & 
    120.0, & 
    150.0, & 
    150.0, & 
    136.0, & 
    250.0, & 
    180.0/   

    REAL(KIND=dp), PARAMETER  :: yields_dens_aer = 1.5 

    
    

    REAL(KIND=dp)             :: B, mw_ratio, dens_ratio

    
    
    
    B = nume / (nume + den + 1.0e-35_dp)

    
    mw_ratio = vbs_mw_prec(voc_idx)/mw_aer

    
    dens_ratio = dens_aer / yields_dens_aer

    vbs_yield = (vbs_alphhiN(bin_idx,voc_idx)  * B +             &
                 vbs_alphlowN(bin_idx,voc_idx) * (1.0_dp - B)) * &
                 dens_ratio * mw_ratio

    END FUNCTION vbs_yield

    SUBROUTINE aero_surfarea( aero_srf_area, aero_diam, rh, temp, &
                              aer_so4, aer_oc2, aer_bc2 )

    IMPLICIT NONE

    
    
    
    REAL(kind=dp), intent(in)  :: rh
    REAL(kind=dp), intent(in)  :: temp
    REAL(kind=dp), intent(in)  :: aer_so4, aer_oc2, aer_bc2
    REAL(kind=dp), intent(out) :: aero_srf_area(3)
    REAL(kind=dp), intent(out) :: aero_diam(3)

    
    
    
    
    real(dp), parameter :: rm_sulf  = 6.95e-6_dp
    real(dp), parameter :: dm_sulf  = 2._dp*rm_sulf
    real(dp), parameter :: sd_sulf  = 2.03_dp

    
    real(dp), parameter :: rm_orgc  = 2.12e-6_dp
    real(dp), parameter :: dm_orgc  = 2._dp*rm_orgc
    real(dp), parameter :: sd_orgc  = 2.20_dp

    
    real(dp), parameter :: rm_bc    = 1.18e-6_dp
    real(dp), parameter :: dm_bc    = 2._dp*rm_bc
    real(dp), parameter :: sd_bc    = 2.00_dp

    real(dp), parameter :: pi       = 3.1415926535897932384626433_dp

    integer  :: irh, rh_l, rh_u
    real(dp) :: log_sd_sulf, log_sd_orgc, log_sd_bc
    real(dp) :: dm_sulf_wet, dm_orgc_wet, dm_bc_wet
    real(dp) :: rfac_sulf, rfac_oc, rfac_bc
    real(dp) :: n, n_exp, factor, s_exp
    
    
    
    
    real(dp), dimension(7) :: table_rh, table_rfac_sulf
    real(dp), dimension(7) :: table_rfac_bc, table_rfac_oc

    data table_rh(1:7) &
        / 0.0_dp, 0.5_dp, 0.7_dp, 0.8_dp, 0.9_dp, 0.95_dp, 0.99_dp /
    data table_rfac_sulf(1:7) &
        / 1.0_dp, 1.4_dp, 1.5_dp, 1.6_dp, 1.8_dp, 1.9_dp,  2.2_dp /
    data table_rfac_oc(1:7) &
        / 1.0_dp, 1.2_dp, 1.4_dp, 1.5_dp, 1.6_dp, 1.8_dp,  2.2_dp /
    data table_rfac_bc(1:7) &
        / 1.0_dp, 1.0_dp, 1.0_dp, 1.2_dp, 1.4_dp, 1.5_dp,  1.9_dp /

    log_sd_sulf = log( sd_sulf )
    log_sd_orgc = log( sd_orgc )
    log_sd_bc   = log( sd_bc )

    
    
    
    n_exp = exp( -4.5_dp*log(sd_sulf)*log(sd_sulf) )
    
    
    
    if (rh >= table_rh(7)) then
      rfac_sulf = table_rfac_sulf(7)
      rfac_oc = table_rfac_oc(7)
      rfac_bc = table_rfac_bc(7)
    else
      do irh = 2,7
        if (rh <= table_rh(irh)) then
          exit
        end if
      end do
      rh_l = irh-1
      rh_u = irh

      factor = (rh - table_rh(rh_l))/(table_rh(rh_u) - table_rh(rh_l))

      rfac_sulf = table_rfac_sulf(rh_l) &
                + factor*(table_rfac_sulf(rh_u) - table_rfac_sulf(rh_l))
      rfac_oc = table_rfac_oc(rh_u) &
              + factor*(table_rfac_oc(rh_u) - table_rfac_oc(rh_l))
      rfac_bc = table_rfac_bc(rh_u) &
              + factor*(table_rfac_bc(rh_u) - table_rfac_bc(rh_l))
    end if

    dm_sulf_wet = dm_sulf * rfac_sulf
    dm_orgc_wet = dm_orgc * rfac_oc
    dm_bc_wet = dm_bc * rfac_bc

    
    dm_bc_wet   = min(dm_bc_wet  ,50.e-6_dp)
    dm_orgc_wet = min(dm_orgc_wet,50.e-6_dp)
    
    aero_diam(:) = (/ dm_sulf_wet, dm_orgc_wet, dm_bc_wet /)

    n = aer_so4 * (6._dp/pi)*(1._dp/(dm_sulf**3))*n_exp
    s_exp = exp( 2._dp*log_sd_sulf*log_sd_sulf )
    aero_srf_area(1) = n * pi * (dm_sulf_wet*dm_sulf_wet) * s_exp

    n = aer_oc2 * (6._dp/pi)*(1._dp/(dm_orgc**3))*n_exp
    s_exp = exp( 2._dp*log_sd_orgc*log_sd_orgc )
    aero_srf_area(2) = n * pi * (dm_orgc_wet*dm_orgc_wet) * s_exp

    n = aer_bc2 * (6._dp/pi)*(1._dp/(dm_bc**3))*n_exp
    s_exp = exp( 2._dp*log_sd_bc*log_sd_bc )
    aero_srf_area(3) = n * pi * (dm_bc_wet*dm_bc_wet) * s_exp

    END SUBROUTINE aero_surfarea






END MODULE t1_mozcart_UpdateRconstWRF

