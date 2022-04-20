! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Linear Algebra Data and Routines File
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
! File                 : t1_mozcart_LinearAlgebra.f90
! Time                 : Thu Apr 14 13:14:36 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/t1_mozcart
! Equation file        : t1_mozcart.kpp
! Output root filename : t1_mozcart
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE t1_mozcart_LinearAlgebra

  USE t1_mozcart_Parameters
  USE t1_mozcart_JacobianSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! KppSolveTR - sparse, transposed back substitution
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      X         - Vector for variables
!      XX        - Vector for output variables
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE KppSolveTR ( JVS, X, XX )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! X - Vector for variables
  REAL(kind=dp) :: X(NVAR)
! XX - Vector for output variables
  REAL(kind=dp) :: XX(NVAR)

  XX(1) = X(1)/JVS(1)
  XX(2) = X(2)/JVS(4)
  XX(3) = X(3)/JVS(6)
  XX(4) = X(4)/JVS(8)
  XX(5) = X(5)/JVS(10)
  XX(6) = X(6)/JVS(12)
  XX(7) = (X(7)-JVS(2)*XX(1))/(JVS(14))
  XX(8) = X(8)/JVS(18)
  XX(9) = X(9)/JVS(21)
  XX(10) = X(10)/JVS(24)
  XX(11) = X(11)/JVS(27)
  XX(12) = X(12)/JVS(29)
  XX(13) = X(13)/JVS(32)
  XX(14) = X(14)/JVS(35)
  XX(15) = X(15)/JVS(37)
  XX(16) = X(16)/JVS(40)
  XX(17) = X(17)/JVS(43)
  XX(18) = X(18)/JVS(46)
  XX(19) = X(19)/JVS(49)
  XX(20) = X(20)/JVS(52)
  XX(21) = (X(21)-JVS(50)*XX(19))/(JVS(56))
  XX(22) = X(22)/JVS(59)
  XX(23) = (X(23)-JVS(15)*XX(7))/(JVS(61))
  XX(24) = X(24)/JVS(64)
  XX(25) = X(25)/JVS(67)
  XX(26) = X(26)/JVS(72)
  XX(27) = X(27)/JVS(77)
  XX(28) = X(28)/JVS(81)
  XX(29) = X(29)/JVS(85)
  XX(30) = X(30)/JVS(90)
  XX(31) = X(31)/JVS(94)
  XX(32) = X(32)/JVS(100)
  XX(33) = X(33)/JVS(104)
  XX(34) = X(34)/JVS(108)
  XX(35) = X(35)/JVS(112)
  XX(36) = X(36)/JVS(119)
  XX(37) = X(37)/JVS(123)
  XX(38) = X(38)/JVS(127)
  XX(39) = X(39)/JVS(131)
  XX(40) = (X(40)-JVS(68)*XX(25))/(JVS(135))
  XX(41) = X(41)/JVS(139)
  XX(42) = X(42)/JVS(143)
  XX(43) = X(43)/JVS(147)
  XX(44) = X(44)/JVS(155)
  XX(45) = X(45)/JVS(158)
  XX(46) = X(46)/JVS(166)
  XX(47) = X(47)/JVS(169)
  XX(48) = X(48)/JVS(173)
  XX(49) = X(49)/JVS(177)
  XX(50) = X(50)/JVS(181)
  XX(51) = X(51)/JVS(185)
  XX(52) = (X(52)-JVS(69)*XX(25)-JVS(136)*XX(40))/(JVS(192))
  XX(53) = X(53)/JVS(196)
  XX(54) = (X(54)-JVS(128)*XX(38))/(JVS(203))
  XX(55) = X(55)/JVS(209)
  XX(56) = X(56)/JVS(218)
  XX(57) = X(57)/JVS(223)
  XX(58) = X(58)/JVS(231)
  XX(59) = X(59)/JVS(238)
  XX(60) = X(60)/JVS(242)
  XX(61) = (X(61)-JVS(174)*XX(48)-JVS(210)*XX(55))/(JVS(248))
  XX(62) = (X(62)-JVS(186)*XX(51))/(JVS(254))
  XX(63) = X(63)/JVS(258)
  XX(64) = X(64)/JVS(264)
  XX(65) = (X(65)-JVS(25)*XX(10))/(JVS(273))
  XX(66) = X(66)/JVS(278)
  XX(67) = X(67)/JVS(284)
  XX(68) = X(68)/JVS(288)
  XX(69) = X(69)/JVS(292)
  XX(70) = X(70)/JVS(297)
  XX(71) = (X(71)-JVS(279)*XX(66))/(JVS(301))
  XX(72) = X(72)/JVS(306)
  XX(73) = X(73)/JVS(311)
  XX(74) = X(74)/JVS(320)
  XX(75) = (X(75)-JVS(19)*XX(8)-JVS(91)*XX(30))/(JVS(325))
  XX(76) = X(76)/JVS(330)
  XX(77) = (X(77)-JVS(82)*XX(28)-JVS(232)*XX(58))/(JVS(339))
  XX(78) = (X(78)-JVS(259)*XX(63))/(JVS(345))
  XX(79) = (X(79)-JVS(95)*XX(31)-JVS(211)*XX(55)-JVS(312)*XX(73))/(JVS(349))
  XX(80) = (X(80)-JVS(260)*XX(63))/(JVS(353))
  XX(81) = (X(81)-JVS(124)*XX(37))/(JVS(361))
  XX(82) = (X(82)-JVS(239)*XX(59))/(JVS(367))
  XX(83) = X(83)/JVS(375)
  XX(84) = X(84)/JVS(391)
  XX(85) = (X(85)-JVS(96)*XX(31)-JVS(140)*XX(41)-JVS(212)*XX(55)-JVS(313)*XX(73))/(JVS(398))
  XX(86) = X(86)/JVS(406)
  XX(87) = (X(87)-JVS(97)*XX(31)-JVS(213)*XX(55)-JVS(314)*XX(73)-JVS(350)*XX(79)-JVS(407)*XX(86))/(JVS(415))
  XX(88) = (X(88)-JVS(22)*XX(9))/(JVS(420))
  XX(89) = (X(89)-JVS(376)*XX(83))/(JVS(426))
  XX(90) = (X(90)-JVS(98)*XX(31)-JVS(141)*XX(41)-JVS(214)*XX(55)-JVS(315)*XX(73)-JVS(399)*XX(85)-JVS(408)*XX(86))&
             &/(JVS(432))
  XX(91) = X(91)/JVS(441)
  XX(92) = (X(92)-JVS(109)*XX(34))/(JVS(448))
  XX(93) = (X(93)-JVS(13)*XX(6)-JVS(38)*XX(15)-JVS(280)*XX(66)-JVS(302)*XX(71))/(JVS(457))
  XX(94) = X(94)/JVS(463)
  XX(95) = (X(95)-JVS(86)*XX(29)-JVS(144)*XX(42)-JVS(464)*XX(94))/(JVS(482))
  XX(96) = X(96)/JVS(489)
  XX(97) = (X(97)-JVS(265)*XX(64))/(JVS(498))
  XX(98) = (X(98)-JVS(377)*XX(83))/(JVS(506))
  XX(99) = X(99)/JVS(513)
  XX(100) = (X(100)-JVS(159)*XX(45)-JVS(378)*XX(83)-JVS(392)*XX(84))/(JVS(523))
  XX(101) = (X(101)-JVS(160)*XX(45)-JVS(379)*XX(83)-JVS(393)*XX(84))/(JVS(527))
  XX(102) = X(102)/JVS(546)
  XX(103) = (X(103)-JVS(547)*XX(102))/(JVS(568))
  XX(104) = (X(104)-JVS(161)*XX(45)-JVS(380)*XX(83)-JVS(394)*XX(84)-JVS(548)*XX(102))/(JVS(575))
  XX(105) = (X(105)-JVS(162)*XX(45)-JVS(381)*XX(83)-JVS(395)*XX(84)-JVS(549)*XX(102))/(JVS(579))
  XX(106) = (X(106)-JVS(163)*XX(45)-JVS(382)*XX(83)-JVS(396)*XX(84)-JVS(550)*XX(102))/(JVS(583))
  XX(107) = (X(107)-JVS(261)*XX(63)-JVS(346)*XX(78)-JVS(354)*XX(80)-JVS(362)*XX(81)-JVS(483)*XX(95))/(JVS(594))
  XX(108) = X(108)/JVS(619)
  XX(109) = X(109)/JVS(656)
  XX(110) = (X(110)-JVS(551)*XX(102)-JVS(657)*XX(109))/(JVS(689))
  XX(111) = (X(111)-JVS(285)*XX(67)-JVS(465)*XX(94)-JVS(499)*XX(97)-JVS(620)*XX(108)-JVS(658)*XX(109)-JVS(690)*XX(110))&
              &/(JVS(707))
  XX(112) = (X(112)-JVS(621)*XX(108)-JVS(659)*XX(109))/(JVS(720))
  XX(113) = X(113)/JVS(737)
  XX(114) = (X(114)-JVS(148)*XX(43)-JVS(303)*XX(71)-JVS(368)*XX(82)-JVS(383)*XX(83)-JVS(458)*XX(93)-JVS(514)*XX(99)&
              &-JVS(660)*XX(109))/(JVS(756))
  XX(115) = (X(115)-JVS(661)*XX(109)-JVS(738)*XX(113))/(JVS(775))
  XX(116) = (X(116)-JVS(73)*XX(26)-JVS(321)*XX(74)-JVS(466)*XX(94)-JVS(622)*XX(108)-JVS(662)*XX(109)-JVS(691)*XX(110)&
              &-JVS(721)*XX(112)-JVS(739)*XX(113))/(JVS(785))
  XX(117) = (X(117)-JVS(316)*XX(73)-JVS(515)*XX(99)-JVS(552)*XX(102)-JVS(569)*XX(103)-JVS(663)*XX(109)-JVS(740)*XX(113))&
              &/(JVS(792))
  XX(118) = (X(118)-JVS(182)*XX(50)-JVS(331)*XX(76)-JVS(467)*XX(94)-JVS(741)*XX(113))/(JVS(805))
  XX(119) = (X(119)-JVS(243)*XX(60)-JVS(332)*XX(76)-JVS(468)*XX(94)-JVS(623)*XX(108)-JVS(742)*XX(113))/(JVS(819))
  XX(120) = (X(120)-JVS(664)*XX(109)-JVS(722)*XX(112)-JVS(786)*XX(116)-JVS(806)*XX(118))/(JVS(836))
  XX(121) = (X(121)-JVS(132)*XX(39)-JVS(469)*XX(94))/(JVS(855))
  XX(122) = (X(122)-JVS(78)*XX(27)-JVS(281)*XX(66)-JVS(459)*XX(93)-JVS(665)*XX(109))/(JVS(868))
  XX(123) = X(123)/JVS(884)
  XX(124) = (X(124)-JVS(74)*XX(26)-JVS(624)*XX(108)-JVS(666)*XX(109)-JVS(692)*XX(110)-JVS(743)*XX(113)-JVS(776)*XX(115)&
              &-JVS(856)*XX(121)-JVS(885)*XX(123))/(JVS(904))
  XX(125) = (X(125)-JVS(120)*XX(36)-JVS(470)*XX(94)-JVS(667)*XX(109)-JVS(693)*XX(110)-JVS(886)*XX(123)-JVS(905)*XX(124))&
              &/(JVS(914))
  XX(126) = (X(126)-JVS(668)*XX(109)-JVS(744)*XX(113))/(JVS(945))
  XX(127) = (X(127)-JVS(53)*XX(20)-JVS(471)*XX(94)-JVS(553)*XX(102)-JVS(669)*XX(109)-JVS(694)*XX(110)-JVS(887)*XX(123)&
              &-JVS(946)*XX(126))/(JVS(966))
  XX(128) = (X(128)-JVS(307)*XX(72)-JVS(421)*XX(88)-JVS(472)*XX(94)-JVS(695)*XX(110)-JVS(745)*XX(113)-JVS(888)*XX(123)&
              &-JVS(906)*XX(124)-JVS(947)*XX(126)-JVS(967)*XX(127))/(JVS(977))
  XX(129) = (X(129)-JVS(170)*XX(47)-JVS(178)*XX(49)-JVS(422)*XX(88)-JVS(473)*XX(94)-JVS(516)*XX(99)-JVS(696)*XX(110)&
              &-JVS(746)*XX(113)-JVS(793)*XX(117)-JVS(948)*XX(126)-JVS(968)*XX(127))/(JVS(986))
  XX(130) = (X(130)-JVS(65)*XX(24)-JVS(113)*XX(35)-JVS(293)*XX(69)-JVS(317)*XX(73)-JVS(423)*XX(88)-JVS(474)*XX(94)&
              &-JVS(517)*XX(99)-JVS(554)*XX(102)-JVS(570)*XX(103)-JVS(670)*XX(109)-JVS(697)*XX(110)-JVS(747)*XX(113)&
              &-JVS(889)*XX(123)-JVS(907)*XX(124)-JVS(949)*XX(126)-JVS(969)*XX(127))/(JVS(995))
  XX(131) = (X(131)-JVS(384)*XX(83)-JVS(671)*XX(109)-JVS(915)*XX(125)-JVS(950)*XX(126))/(JVS(1008))
  XX(132) = (X(132)-JVS(197)*XX(53)-JVS(282)*XX(66)-JVS(385)*XX(83)-JVS(460)*XX(93)-JVS(672)*XX(109)-JVS(748)*XX(113)&
              &-JVS(869)*XX(122))/(JVS(1070))
  XX(133) = (X(133)-JVS(386)*XX(83)-JVS(673)*XX(109)-JVS(757)*XX(114)-JVS(777)*XX(115)-JVS(870)*XX(122)-JVS(890)*XX(123)&
              &-JVS(916)*XX(125)-JVS(951)*XX(126)-JVS(1071)*XX(132))/(JVS(1086))
  XX(134) = (X(134)-JVS(16)*XX(7)-JVS(41)*XX(16)-JVS(62)*XX(23)-JVS(114)*XX(35)-JVS(156)*XX(44)-JVS(298)*XX(70)-JVS(427)&
              &*XX(89)-JVS(490)*XX(96)-JVS(500)*XX(97)-JVS(507)*XX(98)-JVS(518)*XX(99)-JVS(524)*XX(100)-JVS(528)*XX(101)&
              &-JVS(555)*XX(102)-JVS(576)*XX(104)-JVS(580)*XX(105)-JVS(584)*XX(106)-JVS(595)*XX(107)-JVS(625)*XX(108)&
              &-JVS(674)*XX(109)-JVS(698)*XX(110)-JVS(708)*XX(111)-JVS(723)*XX(112)-JVS(749)*XX(113)-JVS(758)*XX(114)&
              &-JVS(778)*XX(115)-JVS(787)*XX(116)-JVS(794)*XX(117)-JVS(807)*XX(118)-JVS(820)*XX(119)-JVS(837)*XX(120)&
              &-JVS(857)*XX(121)-JVS(871)*XX(122)-JVS(891)*XX(123)-JVS(908)*XX(124)-JVS(917)*XX(125)-JVS(952)*XX(126)&
              &-JVS(970)*XX(127)-JVS(978)*XX(128)-JVS(987)*XX(129)-JVS(996)*XX(130)-JVS(1009)*XX(131)-JVS(1072)*XX(132)&
              &-JVS(1087)*XX(133))/(JVS(1127))
  XX(135) = (X(135)-JVS(149)*XX(43)-JVS(164)*XX(45)-JVS(167)*XX(46)-JVS(233)*XX(58)-JVS(304)*XX(71)-JVS(326)*XX(75)&
              &-JVS(340)*XX(77)-JVS(387)*XX(83)-JVS(397)*XX(84)-JVS(428)*XX(89)-JVS(461)*XX(93)-JVS(491)*XX(96)-JVS(501)&
              &*XX(97)-JVS(508)*XX(98)-JVS(525)*XX(100)-JVS(529)*XX(101)-JVS(556)*XX(102)-JVS(577)*XX(104)-JVS(581)*XX(105)&
              &-JVS(585)*XX(106)-JVS(596)*XX(107)-JVS(626)*XX(108)-JVS(675)*XX(109)-JVS(699)*XX(110)-JVS(709)*XX(111)&
              &-JVS(724)*XX(112)-JVS(759)*XX(114)-JVS(779)*XX(115)-JVS(808)*XX(118)-JVS(821)*XX(119)-JVS(838)*XX(120)&
              &-JVS(858)*XX(121)-JVS(872)*XX(122)-JVS(892)*XX(123)-JVS(909)*XX(124)-JVS(953)*XX(126)-JVS(979)*XX(128)&
              &-JVS(988)*XX(129)-JVS(997)*XX(130)-JVS(1010)*XX(131)-JVS(1073)*XX(132)-JVS(1088)*XX(133)-JVS(1128)*XX(134))&
              &/(JVS(1153))
  XX(136) = (X(136)-JVS(115)*XX(35)-JVS(150)*XX(43)-JVS(219)*XX(56)-JVS(224)*XX(57)-JVS(557)*XX(102)-JVS(676)*XX(109)&
              &-JVS(700)*XX(110)-JVS(795)*XX(117)-JVS(893)*XX(123)-JVS(918)*XX(125)-JVS(954)*XX(126)-JVS(971)*XX(127)&
              &-JVS(980)*XX(128)-JVS(989)*XX(129)-JVS(998)*XX(130)-JVS(1011)*XX(131)-JVS(1074)*XX(132)-JVS(1089)*XX(133)&
              &-JVS(1129)*XX(134)-JVS(1154)*XX(135))/(JVS(1205))
  XX(137) = (X(137)-JVS(26)*XX(10)-JVS(42)*XX(16)-JVS(105)*XX(33)-JVS(225)*XX(57)-JVS(234)*XX(58)-JVS(274)*XX(65)&
              &-JVS(289)*XX(68)-JVS(341)*XX(77)-JVS(409)*XX(86)-JVS(442)*XX(91)-JVS(558)*XX(102)-JVS(571)*XX(103)-JVS(677)&
              &*XX(109)-JVS(750)*XX(113)-JVS(873)*XX(122)-JVS(894)*XX(123)-JVS(955)*XX(126)-JVS(1075)*XX(132)-JVS(1130)&
              &*XX(134)-JVS(1155)*XX(135)-JVS(1206)*XX(136))/(JVS(1272))
  XX(138) = (X(138)-JVS(3)*XX(1)-JVS(5)*XX(2)-JVS(7)*XX(3)-JVS(9)*XX(4)-JVS(11)*XX(5)-JVS(17)*XX(7)-JVS(23)*XX(9)&
              &-JVS(28)*XX(11)-JVS(30)*XX(12)-JVS(33)*XX(13)-JVS(36)*XX(14)-JVS(39)*XX(15)-JVS(44)*XX(17)-JVS(47)*XX(18)&
              &-JVS(51)*XX(19)-JVS(54)*XX(20)-JVS(57)*XX(21)-JVS(60)*XX(22)-JVS(63)*XX(23)-JVS(66)*XX(24)-JVS(70)*XX(25)&
              &-JVS(75)*XX(26)-JVS(79)*XX(27)-JVS(83)*XX(28)-JVS(87)*XX(29)-JVS(101)*XX(32)-JVS(106)*XX(33)-JVS(110)*XX(34)&
              &-JVS(116)*XX(35)-JVS(121)*XX(36)-JVS(125)*XX(37)-JVS(129)*XX(38)-JVS(133)*XX(39)-JVS(137)*XX(40)-JVS(145)&
              &*XX(42)-JVS(151)*XX(43)-JVS(157)*XX(44)-JVS(165)*XX(45)-JVS(168)*XX(46)-JVS(171)*XX(47)-JVS(175)*XX(48)&
              &-JVS(179)*XX(49)-JVS(183)*XX(50)-JVS(187)*XX(51)-JVS(193)*XX(52)-JVS(204)*XX(54)-JVS(215)*XX(55)-JVS(220)&
              &*XX(56)-JVS(226)*XX(57)-JVS(235)*XX(58)-JVS(240)*XX(59)-JVS(244)*XX(60)-JVS(249)*XX(61)-JVS(255)*XX(62)&
              &-JVS(262)*XX(63)-JVS(266)*XX(64)-JVS(275)*XX(65)-JVS(283)*XX(66)-JVS(286)*XX(67)-JVS(290)*XX(68)-JVS(294)&
              &*XX(69)-JVS(299)*XX(70)-JVS(305)*XX(71)-JVS(308)*XX(72)-JVS(318)*XX(73)-JVS(322)*XX(74)-JVS(327)*XX(75)&
              &-JVS(333)*XX(76)-JVS(342)*XX(77)-JVS(347)*XX(78)-JVS(351)*XX(79)-JVS(355)*XX(80)-JVS(363)*XX(81)-JVS(369)&
              &*XX(82)-JVS(388)*XX(83)-JVS(400)*XX(85)-JVS(410)*XX(86)-JVS(416)*XX(87)-JVS(424)*XX(88)-JVS(429)*XX(89)&
              &-JVS(433)*XX(90)-JVS(443)*XX(91)-JVS(449)*XX(92)-JVS(462)*XX(93)-JVS(475)*XX(94)-JVS(484)*XX(95)-JVS(492)&
              &*XX(96)-JVS(502)*XX(97)-JVS(509)*XX(98)-JVS(519)*XX(99)-JVS(526)*XX(100)-JVS(530)*XX(101)-JVS(559)*XX(102)&
              &-JVS(572)*XX(103)-JVS(578)*XX(104)-JVS(582)*XX(105)-JVS(586)*XX(106)-JVS(597)*XX(107)-JVS(627)*XX(108)&
              &-JVS(678)*XX(109)-JVS(701)*XX(110)-JVS(710)*XX(111)-JVS(725)*XX(112)-JVS(751)*XX(113)-JVS(760)*XX(114)&
              &-JVS(780)*XX(115)-JVS(788)*XX(116)-JVS(796)*XX(117)-JVS(809)*XX(118)-JVS(822)*XX(119)-JVS(839)*XX(120)&
              &-JVS(859)*XX(121)-JVS(874)*XX(122)-JVS(895)*XX(123)-JVS(910)*XX(124)-JVS(919)*XX(125)-JVS(956)*XX(126)&
              &-JVS(972)*XX(127)-JVS(981)*XX(128)-JVS(990)*XX(129)-JVS(999)*XX(130)-JVS(1012)*XX(131)-JVS(1076)*XX(132)&
              &-JVS(1090)*XX(133)-JVS(1131)*XX(134)-JVS(1156)*XX(135)-JVS(1207)*XX(136)-JVS(1273)*XX(137))/(JVS(1397))
  XX(139) = (X(139)-JVS(152)*XX(43)-JVS(221)*XX(56)-JVS(291)*XX(68)-JVS(896)*XX(123)-JVS(957)*XX(126)-JVS(1077)*XX(132)&
              &-JVS(1132)*XX(134)-JVS(1157)*XX(135)-JVS(1208)*XX(136)-JVS(1274)*XX(137)-JVS(1398)*XX(138))/(JVS(1412))
  XX(140) = (X(140)-JVS(88)*XX(29)-JVS(102)*XX(32)-JVS(117)*XX(35)-JVS(153)*XX(43)-JVS(180)*XX(49)-JVS(334)*XX(76)&
              &-JVS(450)*XX(92)-JVS(476)*XX(94)-JVS(485)*XX(95)-JVS(503)*XX(97)-JVS(520)*XX(99)-JVS(560)*XX(102)-JVS(628)&
              &*XX(108)-JVS(679)*XX(109)-JVS(702)*XX(110)-JVS(711)*XX(111)-JVS(726)*XX(112)-JVS(752)*XX(113)-JVS(781)&
              &*XX(115)-JVS(789)*XX(116)-JVS(797)*XX(117)-JVS(810)*XX(118)-JVS(823)*XX(119)-JVS(840)*XX(120)-JVS(860)&
              &*XX(121)-JVS(897)*XX(123)-JVS(920)*XX(125)-JVS(958)*XX(126)-JVS(973)*XX(127)-JVS(982)*XX(128)-JVS(991)&
              &*XX(129)-JVS(1000)*XX(130)-JVS(1013)*XX(131)-JVS(1078)*XX(132)-JVS(1091)*XX(133)-JVS(1133)*XX(134)-JVS(1158)&
              &*XX(135)-JVS(1209)*XX(136)-JVS(1275)*XX(137)-JVS(1399)*XX(138)-JVS(1413)*XX(139))/(JVS(1455))
  XX(141) = (X(141)-JVS(20)*XX(8)-JVS(55)*XX(20)-JVS(80)*XX(27)-JVS(84)*XX(28)-JVS(92)*XX(30)-JVS(103)*XX(32)-JVS(107)&
              &*XX(33)-JVS(111)*XX(34)-JVS(122)*XX(36)-JVS(126)*XX(37)-JVS(130)*XX(38)-JVS(134)*XX(39)-JVS(138)*XX(40)&
              &-JVS(146)*XX(42)-JVS(154)*XX(43)-JVS(172)*XX(47)-JVS(176)*XX(48)-JVS(184)*XX(50)-JVS(188)*XX(51)-JVS(194)&
              &*XX(52)-JVS(198)*XX(53)-JVS(205)*XX(54)-JVS(216)*XX(55)-JVS(222)*XX(56)-JVS(236)*XX(58)-JVS(241)*XX(59)&
              &-JVS(245)*XX(60)-JVS(250)*XX(61)-JVS(256)*XX(62)-JVS(267)*XX(64)-JVS(276)*XX(65)-JVS(287)*XX(67)-JVS(323)&
              &*XX(74)-JVS(328)*XX(75)-JVS(343)*XX(77)-JVS(348)*XX(78)-JVS(352)*XX(79)-JVS(364)*XX(81)-JVS(370)*XX(82)&
              &-JVS(389)*XX(83)-JVS(401)*XX(85)-JVS(411)*XX(86)-JVS(417)*XX(87)-JVS(425)*XX(88)-JVS(434)*XX(90)-JVS(444)&
              &*XX(91)-JVS(451)*XX(92)-JVS(486)*XX(95)-JVS(493)*XX(96)-JVS(504)*XX(97)-JVS(521)*XX(99)-JVS(561)*XX(102)&
              &-JVS(573)*XX(103)-JVS(598)*XX(107)-JVS(629)*XX(108)-JVS(680)*XX(109)-JVS(703)*XX(110)-JVS(712)*XX(111)&
              &-JVS(727)*XX(112)-JVS(753)*XX(113)-JVS(782)*XX(115)-JVS(790)*XX(116)-JVS(811)*XX(118)-JVS(824)*XX(119)&
              &-JVS(841)*XX(120)-JVS(861)*XX(121)-JVS(875)*XX(122)-JVS(898)*XX(123)-JVS(911)*XX(124)-JVS(921)*XX(125)&
              &-JVS(959)*XX(126)-JVS(974)*XX(127)-JVS(983)*XX(128)-JVS(992)*XX(129)-JVS(1001)*XX(130)-JVS(1014)*XX(131)&
              &-JVS(1079)*XX(132)-JVS(1092)*XX(133)-JVS(1134)*XX(134)-JVS(1159)*XX(135)-JVS(1210)*XX(136)-JVS(1276)*XX(137)&
              &-JVS(1400)*XX(138)-JVS(1414)*XX(139)-JVS(1456)*XX(140))/(JVS(1570))
  XX(142) = (X(142)-JVS(71)*XX(25)-JVS(76)*XX(26)-JVS(93)*XX(30)-JVS(99)*XX(31)-JVS(118)*XX(35)-JVS(142)*XX(41)-JVS(195)&
              &*XX(52)-JVS(199)*XX(53)-JVS(206)*XX(54)-JVS(217)*XX(55)-JVS(237)*XX(58)-JVS(251)*XX(61)-JVS(257)*XX(62)&
              &-JVS(263)*XX(63)-JVS(268)*XX(64)-JVS(277)*XX(65)-JVS(295)*XX(69)-JVS(300)*XX(70)-JVS(309)*XX(72)-JVS(319)&
              &*XX(73)-JVS(329)*XX(75)-JVS(335)*XX(76)-JVS(344)*XX(77)-JVS(356)*XX(80)-JVS(365)*XX(81)-JVS(371)*XX(82)&
              &-JVS(390)*XX(83)-JVS(412)*XX(86)-JVS(418)*XX(87)-JVS(435)*XX(90)-JVS(445)*XX(91)-JVS(452)*XX(92)-JVS(487)&
              &*XX(95)-JVS(494)*XX(96)-JVS(505)*XX(97)-JVS(522)*XX(99)-JVS(562)*XX(102)-JVS(574)*XX(103)-JVS(599)*XX(107)&
              &-JVS(630)*XX(108)-JVS(681)*XX(109)-JVS(704)*XX(110)-JVS(713)*XX(111)-JVS(728)*XX(112)-JVS(754)*XX(113)&
              &-JVS(783)*XX(115)-JVS(791)*XX(116)-JVS(798)*XX(117)-JVS(812)*XX(118)-JVS(825)*XX(119)-JVS(842)*XX(120)&
              &-JVS(862)*XX(121)-JVS(876)*XX(122)-JVS(899)*XX(123)-JVS(912)*XX(124)-JVS(922)*XX(125)-JVS(960)*XX(126)&
              &-JVS(975)*XX(127)-JVS(984)*XX(128)-JVS(993)*XX(129)-JVS(1002)*XX(130)-JVS(1015)*XX(131)-JVS(1080)*XX(132)&
              &-JVS(1093)*XX(133)-JVS(1135)*XX(134)-JVS(1160)*XX(135)-JVS(1211)*XX(136)-JVS(1277)*XX(137)-JVS(1401)*XX(138)&
              &-JVS(1415)*XX(139)-JVS(1457)*XX(140)-JVS(1571)*XX(141))/(JVS(1622))
  XX(142) = XX(142)
  XX(141) = XX(141)-JVS(1621)*XX(142)
  XX(140) = XX(140)-JVS(1569)*XX(141)-JVS(1620)*XX(142)
  XX(139) = XX(139)-JVS(1454)*XX(140)-JVS(1568)*XX(141)-JVS(1619)*XX(142)
  XX(138) = XX(138)-JVS(1411)*XX(139)-JVS(1453)*XX(140)-JVS(1567)*XX(141)-JVS(1618)*XX(142)
  XX(137) = XX(137)-JVS(1396)*XX(138)-JVS(1410)*XX(139)-JVS(1452)*XX(140)-JVS(1566)*XX(141)-JVS(1617)*XX(142)
  XX(136) = XX(136)-JVS(1271)*XX(137)-JVS(1395)*XX(138)-JVS(1409)*XX(139)-JVS(1451)*XX(140)-JVS(1565)*XX(141)-JVS(1616)&
              &*XX(142)
  XX(135) = XX(135)-JVS(1204)*XX(136)-JVS(1270)*XX(137)-JVS(1394)*XX(138)-JVS(1408)*XX(139)-JVS(1450)*XX(140)-JVS(1564)&
              &*XX(141)-JVS(1615)*XX(142)
  XX(134) = XX(134)-JVS(1152)*XX(135)-JVS(1203)*XX(136)-JVS(1269)*XX(137)-JVS(1393)*XX(138)-JVS(1407)*XX(139)-JVS(1449)&
              &*XX(140)-JVS(1563)*XX(141)-JVS(1614)*XX(142)
  XX(133) = XX(133)-JVS(1126)*XX(134)-JVS(1151)*XX(135)-JVS(1202)*XX(136)-JVS(1268)*XX(137)-JVS(1392)*XX(138)-JVS(1406)&
              &*XX(139)-JVS(1448)*XX(140)-JVS(1562)*XX(141)-JVS(1613)*XX(142)
  XX(132) = XX(132)-JVS(1125)*XX(134)-JVS(1150)*XX(135)-JVS(1267)*XX(137)-JVS(1391)*XX(138)-JVS(1447)*XX(140)-JVS(1561)&
              &*XX(141)-JVS(1612)*XX(142)
  XX(131) = XX(131)-JVS(1069)*XX(132)-JVS(1124)*XX(134)-JVS(1149)*XX(135)-JVS(1201)*XX(136)-JVS(1266)*XX(137)-JVS(1390)&
              &*XX(138)-JVS(1405)*XX(139)-JVS(1446)*XX(140)-JVS(1560)*XX(141)-JVS(1611)*XX(142)
  XX(130) = XX(130)-JVS(1007)*XX(131)-JVS(1068)*XX(132)-JVS(1085)*XX(133)-JVS(1123)*XX(134)-JVS(1200)*XX(136)-JVS(1265)&
              &*XX(137)-JVS(1389)*XX(138)-JVS(1445)*XX(140)-JVS(1559)*XX(141)-JVS(1610)*XX(142)
  XX(129) = XX(129)-JVS(1006)*XX(131)-JVS(1067)*XX(132)-JVS(1084)*XX(133)-JVS(1122)*XX(134)-JVS(1199)*XX(136)-JVS(1264)&
              &*XX(137)-JVS(1388)*XX(138)-JVS(1444)*XX(140)-JVS(1558)*XX(141)-JVS(1609)*XX(142)
  XX(128) = XX(128)-JVS(1005)*XX(131)-JVS(1066)*XX(132)-JVS(1083)*XX(133)-JVS(1121)*XX(134)-JVS(1198)*XX(136)-JVS(1263)&
              &*XX(137)-JVS(1387)*XX(138)-JVS(1443)*XX(140)-JVS(1557)*XX(141)-JVS(1608)*XX(142)
  XX(127) = XX(127)-JVS(1065)*XX(132)-JVS(1120)*XX(134)-JVS(1197)*XX(136)-JVS(1262)*XX(137)-JVS(1386)*XX(138)-JVS(1442)&
              &*XX(140)-JVS(1556)*XX(141)-JVS(1607)*XX(142)
  XX(126) = XX(126)-JVS(1119)*XX(134)-JVS(1196)*XX(136)-JVS(1261)*XX(137)-JVS(1385)*XX(138)-JVS(1555)*XX(141)
  XX(125) = XX(125)-JVS(944)*XX(126)-JVS(1064)*XX(132)-JVS(1118)*XX(134)-JVS(1195)*XX(136)-JVS(1260)*XX(137)-JVS(1384)&
              &*XX(138)-JVS(1404)*XX(139)-JVS(1441)*XX(140)-JVS(1554)*XX(141)-JVS(1606)*XX(142)
  XX(124) = XX(124)-JVS(943)*XX(126)-JVS(1063)*XX(132)-JVS(1117)*XX(134)-JVS(1194)*XX(136)-JVS(1259)*XX(137)-JVS(1383)&
              &*XX(138)-JVS(1440)*XX(140)-JVS(1553)*XX(141)-JVS(1605)*XX(142)
  XX(123) = XX(123)-JVS(942)*XX(126)-JVS(1062)*XX(132)-JVS(1193)*XX(136)-JVS(1382)*XX(138)-JVS(1552)*XX(141)
  XX(122) = XX(122)-JVS(1061)*XX(132)-JVS(1116)*XX(134)-JVS(1148)*XX(135)-JVS(1258)*XX(137)-JVS(1381)*XX(138)-JVS(1439)&
              &*XX(140)-JVS(1551)*XX(141)-JVS(1604)*XX(142)
  XX(121) = XX(121)-JVS(883)*XX(123)-JVS(941)*XX(126)-JVS(1060)*XX(132)-JVS(1192)*XX(136)-JVS(1257)*XX(137)-JVS(1380)&
              &*XX(138)-JVS(1438)*XX(140)-JVS(1550)*XX(141)-JVS(1603)*XX(142)
  XX(120) = XX(120)-JVS(854)*XX(121)-JVS(1059)*XX(132)-JVS(1115)*XX(134)-JVS(1191)*XX(136)-JVS(1256)*XX(137)-JVS(1379)&
              &*XX(138)-JVS(1437)*XX(140)-JVS(1549)*XX(141)-JVS(1602)*XX(142)
  XX(119) = XX(119)-JVS(835)*XX(120)-JVS(853)*XX(121)-JVS(1058)*XX(132)-JVS(1114)*XX(134)-JVS(1190)*XX(136)-JVS(1255)&
              &*XX(137)-JVS(1378)*XX(138)-JVS(1436)*XX(140)-JVS(1548)*XX(141)-JVS(1601)*XX(142)
  XX(118) = XX(118)-JVS(834)*XX(120)-JVS(1057)*XX(132)-JVS(1113)*XX(134)-JVS(1254)*XX(137)-JVS(1377)*XX(138)-JVS(1435)&
              &*XX(140)-JVS(1547)*XX(141)-JVS(1600)*XX(142)
  XX(117) = XX(117)-JVS(940)*XX(126)-JVS(1056)*XX(132)-JVS(1112)*XX(134)-JVS(1189)*XX(136)-JVS(1253)*XX(137)-JVS(1376)&
              &*XX(138)-JVS(1434)*XX(140)-JVS(1546)*XX(141)-JVS(1599)*XX(142)
  XX(116) = XX(116)-JVS(852)*XX(121)-JVS(1055)*XX(132)-JVS(1111)*XX(134)-JVS(1188)*XX(136)-JVS(1252)*XX(137)-JVS(1375)&
              &*XX(138)-JVS(1433)*XX(140)-JVS(1545)*XX(141)-JVS(1598)*XX(142)
  XX(115) = XX(115)-JVS(1110)*XX(134)-JVS(1187)*XX(136)-JVS(1251)*XX(137)-JVS(1374)*XX(138)-JVS(1432)*XX(140)-JVS(1544)&
              &*XX(141)
  XX(114) = XX(114)-JVS(774)*XX(115)-JVS(867)*XX(122)-JVS(882)*XX(123)-JVS(939)*XX(126)-JVS(1054)*XX(132)-JVS(1109)&
              &*XX(134)-JVS(1147)*XX(135)-JVS(1186)*XX(136)-JVS(1250)*XX(137)-JVS(1373)*XX(138)-JVS(1431)*XX(140)-JVS(1543)&
              &*XX(141)-JVS(1597)*XX(142)
  XX(113) = XX(113)-JVS(1108)*XX(134)-JVS(1249)*XX(137)-JVS(1372)*XX(138)
  XX(112) = XX(112)-JVS(851)*XX(121)-JVS(1053)*XX(132)-JVS(1185)*XX(136)-JVS(1371)*XX(138)-JVS(1430)*XX(140)-JVS(1542)&
              &*XX(141)
  XX(111) = XX(111)-JVS(850)*XX(121)-JVS(1052)*XX(132)-JVS(1184)*XX(136)-JVS(1248)*XX(137)-JVS(1370)*XX(138)-JVS(1429)&
              &*XX(140)-JVS(1541)*XX(141)-JVS(1596)*XX(142)
  XX(110) = XX(110)-JVS(1051)*XX(132)-JVS(1369)*XX(138)-JVS(1540)*XX(141)
  XX(109) = XX(109)-JVS(1368)*XX(138)-JVS(1539)*XX(141)
  XX(108) = XX(108)-JVS(849)*XX(121)-JVS(1183)*XX(136)-JVS(1367)*XX(138)-JVS(1428)*XX(140)
  XX(107) = XX(107)-JVS(618)*XX(108)-JVS(773)*XX(115)-JVS(1050)*XX(132)-JVS(1182)*XX(136)-JVS(1247)*XX(137)-JVS(1366)&
              &*XX(138)-JVS(1427)*XX(140)-JVS(1538)*XX(141)-JVS(1595)*XX(142)
  XX(106) = XX(106)-JVS(593)*XX(107)-JVS(617)*XX(108)-JVS(655)*XX(109)-JVS(719)*XX(112)-JVS(804)*XX(118)-JVS(818)&
              &*XX(119)-JVS(833)*XX(120)-JVS(848)*XX(121)-JVS(938)*XX(126)-JVS(1049)*XX(132)-JVS(1107)*XX(134)-JVS(1146)&
              &*XX(135)-JVS(1181)*XX(136)-JVS(1365)*XX(138)-JVS(1537)*XX(141)
  XX(105) = XX(105)-JVS(592)*XX(107)-JVS(616)*XX(108)-JVS(654)*XX(109)-JVS(718)*XX(112)-JVS(803)*XX(118)-JVS(817)&
              &*XX(119)-JVS(832)*XX(120)-JVS(847)*XX(121)-JVS(937)*XX(126)-JVS(1048)*XX(132)-JVS(1106)*XX(134)-JVS(1145)&
              &*XX(135)-JVS(1180)*XX(136)-JVS(1364)*XX(138)-JVS(1536)*XX(141)
  XX(104) = XX(104)-JVS(591)*XX(107)-JVS(615)*XX(108)-JVS(653)*XX(109)-JVS(717)*XX(112)-JVS(802)*XX(118)-JVS(816)&
              &*XX(119)-JVS(831)*XX(120)-JVS(846)*XX(121)-JVS(936)*XX(126)-JVS(1047)*XX(132)-JVS(1105)*XX(134)-JVS(1144)&
              &*XX(135)-JVS(1179)*XX(136)-JVS(1363)*XX(138)-JVS(1535)*XX(141)
  XX(103) = XX(103)-JVS(652)*XX(109)-JVS(935)*XX(126)-JVS(1246)*XX(137)-JVS(1362)*XX(138)-JVS(1426)*XX(140)-JVS(1534)&
              &*XX(141)-JVS(1594)*XX(142)
  XX(102) = XX(102)-JVS(651)*XX(109)-JVS(1361)*XX(138)-JVS(1533)*XX(141)
  XX(101) = XX(101)-JVS(545)*XX(102)-JVS(590)*XX(107)-JVS(614)*XX(108)-JVS(650)*XX(109)-JVS(716)*XX(112)-JVS(801)&
              &*XX(118)-JVS(815)*XX(119)-JVS(830)*XX(120)-JVS(845)*XX(121)-JVS(934)*XX(126)-JVS(1046)*XX(132)-JVS(1104)&
              &*XX(134)-JVS(1143)*XX(135)-JVS(1178)*XX(136)-JVS(1360)*XX(138)-JVS(1532)*XX(141)
  XX(100) = XX(100)-JVS(544)*XX(102)-JVS(589)*XX(107)-JVS(613)*XX(108)-JVS(649)*XX(109)-JVS(715)*XX(112)-JVS(800)&
              &*XX(118)-JVS(814)*XX(119)-JVS(829)*XX(120)-JVS(844)*XX(121)-JVS(933)*XX(126)-JVS(1045)*XX(132)-JVS(1103)&
              &*XX(134)-JVS(1142)*XX(135)-JVS(1177)*XX(136)-JVS(1359)*XX(138)-JVS(1531)*XX(141)
  XX(99) = XX(99)-JVS(932)*XX(126)-JVS(1044)*XX(132)-JVS(1176)*XX(136)-JVS(1245)*XX(137)-JVS(1358)*XX(138)
  XX(98) = XX(98)-JVS(648)*XX(109)-JVS(755)*XX(114)-JVS(976)*XX(128)-JVS(985)*XX(129)-JVS(994)*XX(130)-JVS(1004)*XX(131)&
             &-JVS(1043)*XX(132)-JVS(1082)*XX(133)-JVS(1102)*XX(134)-JVS(1141)*XX(135)-JVS(1175)*XX(136)-JVS(1357)*XX(138)&
             &-JVS(1425)*XX(140)-JVS(1530)*XX(141)
  XX(97) = XX(97)-JVS(612)*XX(108)-JVS(647)*XX(109)-JVS(1244)*XX(137)-JVS(1356)*XX(138)-JVS(1529)*XX(141)-JVS(1593)&
             &*XX(142)
  XX(96) = XX(96)-JVS(497)*XX(97)-JVS(611)*XX(108)-JVS(903)*XX(124)-JVS(1042)*XX(132)-JVS(1101)*XX(134)-JVS(1243)&
             &*XX(137)-JVS(1528)*XX(141)-JVS(1592)*XX(142)
  XX(95) = XX(95)-JVS(772)*XX(115)-JVS(1041)*XX(132)-JVS(1242)*XX(137)-JVS(1355)*XX(138)-JVS(1424)*XX(140)-JVS(1527)&
             &*XX(141)-JVS(1591)*XX(142)
  XX(94) = XX(94)-JVS(1040)*XX(132)-JVS(1354)*XX(138)-JVS(1526)*XX(141)
  XX(93) = XX(93)-JVS(866)*XX(122)-JVS(1039)*XX(132)-JVS(1140)*XX(135)-JVS(1353)*XX(138)-JVS(1423)*XX(140)-JVS(1525)&
             &*XX(141)-JVS(1590)*XX(142)
  XX(92) = XX(92)-JVS(610)*XX(108)-JVS(771)*XX(115)-JVS(1038)*XX(132)-JVS(1241)*XX(137)-JVS(1352)*XX(138)-JVS(1422)&
             &*XX(140)-JVS(1524)*XX(141)-JVS(1589)*XX(142)
  XX(91) = XX(91)-JVS(646)*XX(109)-JVS(931)*XX(126)-JVS(1240)*XX(137)-JVS(1351)*XX(138)-JVS(1421)*XX(140)-JVS(1523)&
             &*XX(141)-JVS(1588)*XX(142)
  XX(90) = XX(90)-JVS(440)*XX(91)-JVS(543)*XX(102)-JVS(567)*XX(103)-JVS(645)*XX(109)-JVS(930)*XX(126)-JVS(1174)*XX(136)&
             &-JVS(1239)*XX(137)-JVS(1350)*XX(138)-JVS(1522)*XX(141)-JVS(1587)*XX(142)
  XX(89) = XX(89)-JVS(488)*XX(96)-JVS(496)*XX(97)-JVS(609)*XX(108)-JVS(644)*XX(109)-JVS(706)*XX(111)-JVS(1037)*XX(132)&
             &-JVS(1100)*XX(134)-JVS(1139)*XX(135)-JVS(1349)*XX(138)-JVS(1521)*XX(141)
  XX(88) = XX(88)-JVS(965)*XX(127)-JVS(1003)*XX(131)-JVS(1036)*XX(132)-JVS(1081)*XX(133)-JVS(1348)*XX(138)-JVS(1520)&
             &*XX(141)
  XX(87) = XX(87)-JVS(439)*XX(91)-JVS(542)*XX(102)-JVS(566)*XX(103)-JVS(643)*XX(109)-JVS(929)*XX(126)-JVS(1238)*XX(137)&
             &-JVS(1347)*XX(138)-JVS(1519)*XX(141)-JVS(1586)*XX(142)
  XX(86) = XX(86)-JVS(541)*XX(102)-JVS(642)*XX(109)-JVS(1237)*XX(137)-JVS(1518)*XX(141)-JVS(1585)*XX(142)
  XX(85) = XX(85)-JVS(405)*XX(86)-JVS(431)*XX(90)-JVS(438)*XX(91)-JVS(540)*XX(102)-JVS(565)*XX(103)-JVS(641)*XX(109)&
             &-JVS(928)*XX(126)-JVS(1173)*XX(136)-JVS(1346)*XX(138)-JVS(1517)*XX(141)
  XX(84) = XX(84)-JVS(539)*XX(102)-JVS(640)*XX(109)-JVS(927)*XX(126)-JVS(1172)*XX(136)-JVS(1516)*XX(141)
  XX(83) = XX(83)-JVS(1345)*XX(138)-JVS(1515)*XX(141)
  XX(82) = XX(82)-JVS(770)*XX(115)-JVS(881)*XX(123)-JVS(1035)*XX(132)-JVS(1236)*XX(137)-JVS(1344)*XX(138)-JVS(1514)&
             &*XX(141)-JVS(1584)*XX(142)
  XX(81) = XX(81)-JVS(769)*XX(115)-JVS(1171)*XX(136)-JVS(1235)*XX(137)-JVS(1343)*XX(138)-JVS(1513)*XX(141)-JVS(1583)&
             &*XX(142)
  XX(80) = XX(80)-JVS(360)*XX(81)-JVS(481)*XX(95)-JVS(608)*XX(108)-JVS(768)*XX(115)-JVS(1034)*XX(132)-JVS(1170)*XX(136)&
             &-JVS(1234)*XX(137)-JVS(1342)*XX(138)-JVS(1512)*XX(141)
  XX(79) = XX(79)-JVS(404)*XX(86)-JVS(414)*XX(87)-JVS(437)*XX(91)-JVS(538)*XX(102)-JVS(564)*XX(103)-JVS(639)*XX(109)&
             &-JVS(926)*XX(126)-JVS(1341)*XX(138)-JVS(1511)*XX(141)
  XX(78) = XX(78)-JVS(359)*XX(81)-JVS(480)*XX(95)-JVS(588)*XX(107)-JVS(607)*XX(108)-JVS(767)*XX(115)-JVS(1033)*XX(132)&
             &-JVS(1169)*XX(136)-JVS(1340)*XX(138)-JVS(1510)*XX(141)
  XX(77) = XX(77)-JVS(1138)*XX(135)-JVS(1233)*XX(137)-JVS(1339)*XX(138)-JVS(1509)*XX(141)-JVS(1582)*XX(142)
  XX(76) = XX(76)-JVS(736)*XX(113)-JVS(828)*XX(120)-JVS(1232)*XX(137)-JVS(1338)*XX(138)-JVS(1508)*XX(141)
  XX(75) = XX(75)-JVS(688)*XX(110)-JVS(1032)*XX(132)-JVS(1231)*XX(137)-JVS(1337)*XX(138)-JVS(1507)*XX(141)-JVS(1581)&
             &*XX(142)
  XX(74) = XX(74)-JVS(606)*XX(108)-JVS(638)*XX(109)-JVS(687)*XX(110)-JVS(714)*XX(112)-JVS(784)*XX(116)-JVS(1031)*XX(132)&
             &-JVS(1336)*XX(138)-JVS(1506)*XX(141)
  XX(73) = XX(73)-JVS(563)*XX(103)-JVS(637)*XX(109)-JVS(1505)*XX(141)
  XX(72) = XX(72)-JVS(686)*XX(110)-JVS(735)*XX(113)-JVS(880)*XX(123)-JVS(902)*XX(124)-JVS(1030)*XX(132)-JVS(1230)&
             &*XX(137)-JVS(1335)*XX(138)-JVS(1504)*XX(141)
  XX(71) = XX(71)-JVS(456)*XX(93)-JVS(865)*XX(122)-JVS(1029)*XX(132)-JVS(1334)*XX(138)-JVS(1420)*XX(140)-JVS(1503)&
             &*XX(141)
  XX(70) = XX(70)-JVS(605)*XX(108)-JVS(766)*XX(115)-JVS(901)*XX(124)-JVS(1028)*XX(132)-JVS(1229)*XX(137)-JVS(1502)&
             &*XX(141)-JVS(1580)*XX(142)
  XX(69) = XX(69)-JVS(512)*XX(99)-JVS(685)*XX(110)-JVS(734)*XX(113)-JVS(879)*XX(123)-JVS(900)*XX(124)-JVS(1333)*XX(138)&
             &-JVS(1501)*XX(141)
  XX(68) = XX(68)-JVS(878)*XX(123)-JVS(1027)*XX(132)-JVS(1099)*XX(134)-JVS(1228)*XX(137)-JVS(1332)*XX(138)-JVS(1403)&
             &*XX(139)-JVS(1500)*XX(141)
  XX(67) = XX(67)-JVS(495)*XX(97)-JVS(604)*XX(108)-JVS(684)*XX(110)-JVS(705)*XX(111)-JVS(1026)*XX(132)-JVS(1331)*XX(138)&
             &-JVS(1499)*XX(141)
  XX(66) = XX(66)-JVS(455)*XX(93)-JVS(864)*XX(122)-JVS(1330)*XX(138)-JVS(1498)*XX(141)
  XX(65) = XX(65)-JVS(338)*XX(77)-JVS(1227)*XX(137)-JVS(1329)*XX(138)-JVS(1497)*XX(141)-JVS(1579)*XX(142)
  XX(64) = XX(64)-JVS(603)*XX(108)-JVS(1226)*XX(137)-JVS(1328)*XX(138)-JVS(1496)*XX(141)-JVS(1578)*XX(142)
  XX(63) = XX(63)-JVS(358)*XX(81)-JVS(479)*XX(95)-JVS(1168)*XX(136)-JVS(1327)*XX(138)
  XX(62) = XX(62)-JVS(537)*XX(102)-JVS(925)*XX(126)-JVS(1225)*XX(137)-JVS(1326)*XX(138)-JVS(1495)*XX(141)-JVS(1577)&
             &*XX(142)
  XX(61) = XX(61)-JVS(403)*XX(86)-JVS(536)*XX(102)-JVS(1224)*XX(137)-JVS(1325)*XX(138)-JVS(1494)*XX(141)-JVS(1576)&
             &*XX(142)
  XX(60) = XX(60)-JVS(602)*XX(108)-JVS(813)*XX(119)-JVS(827)*XX(120)-JVS(1025)*XX(132)-JVS(1324)*XX(138)-JVS(1493)&
             &*XX(141)
  XX(59) = XX(59)-JVS(366)*XX(82)-JVS(765)*XX(115)-JVS(877)*XX(123)-JVS(1024)*XX(132)-JVS(1323)*XX(138)-JVS(1492)&
             &*XX(141)
  XX(58) = XX(58)-JVS(337)*XX(77)-JVS(1137)*XX(135)-JVS(1223)*XX(137)
  XX(57) = XX(57)-JVS(1023)*XX(132)-JVS(1098)*XX(134)-JVS(1167)*XX(136)-JVS(1222)*XX(137)-JVS(1322)*XX(138)-JVS(1419)&
             &*XX(140)
  XX(56) = XX(56)-JVS(1022)*XX(132)-JVS(1166)*XX(136)-JVS(1321)*XX(138)-JVS(1418)*XX(140)
  XX(55) = XX(55)-JVS(402)*XX(86)-JVS(1491)*XX(141)
  XX(54) = XX(54)-JVS(535)*XX(102)-JVS(1221)*XX(137)-JVS(1320)*XX(138)-JVS(1490)*XX(141)-JVS(1575)*XX(142)
  XX(53) = XX(53)-JVS(374)*XX(83)-JVS(1021)*XX(132)-JVS(1220)*XX(137)-JVS(1489)*XX(141)-JVS(1574)*XX(142)
  XX(52) = XX(52)-JVS(272)*XX(65)-JVS(1219)*XX(137)-JVS(1319)*XX(138)-JVS(1488)*XX(141)-JVS(1573)*XX(142)
  XX(51) = XX(51)-JVS(253)*XX(62)-JVS(534)*XX(102)-JVS(924)*XX(126)-JVS(1318)*XX(138)-JVS(1487)*XX(141)
  XX(50) = XX(50)-JVS(733)*XX(113)-JVS(799)*XX(118)-JVS(826)*XX(120)-JVS(1218)*XX(137)-JVS(1317)*XX(138)
  XX(49) = XX(49)-JVS(511)*XX(99)-JVS(683)*XX(110)-JVS(732)*XX(113)-JVS(1316)*XX(138)-JVS(1486)*XX(141)
  XX(48) = XX(48)-JVS(208)*XX(55)-JVS(247)*XX(61)-JVS(533)*XX(102)-JVS(1315)*XX(138)-JVS(1485)*XX(141)
  XX(47) = XX(47)-JVS(419)*XX(88)-JVS(510)*XX(99)-JVS(1217)*XX(137)-JVS(1314)*XX(138)-JVS(1484)*XX(141)
  XX(46) = XX(46)-JVS(324)*XX(75)-JVS(373)*XX(83)-JVS(636)*XX(109)-JVS(1020)*XX(132)-JVS(1136)*XX(135)-JVS(1313)*XX(138)&
             &-JVS(1483)*XX(141)
  XX(45) = XX(45)-JVS(587)*XX(107)-JVS(1312)*XX(138)
  XX(44) = XX(44)-JVS(296)*XX(70)-JVS(601)*XX(108)-JVS(764)*XX(115)-JVS(1019)*XX(132)-JVS(1097)*XX(134)-JVS(1216)&
             &*XX(137)-JVS(1311)*XX(138)
  XX(43) = XX(43)-JVS(1310)*XX(138)-JVS(1417)*XX(140)
  XX(42) = XX(42)-JVS(478)*XX(95)-JVS(763)*XX(115)-JVS(1309)*XX(138)-JVS(1482)*XX(141)
  XX(41) = XX(41)-JVS(635)*XX(109)-JVS(923)*XX(126)-JVS(1165)*XX(136)-JVS(1481)*XX(141)
  XX(40) = XX(40)-JVS(191)*XX(52)-JVS(271)*XX(65)-JVS(1308)*XX(138)-JVS(1480)*XX(141)
  XX(39) = XX(39)-JVS(843)*XX(121)-JVS(1018)*XX(132)-JVS(1164)*XX(136)-JVS(1307)*XX(138)
  XX(38) = XX(38)-JVS(202)*XX(54)-JVS(532)*XX(102)-JVS(1306)*XX(138)-JVS(1479)*XX(141)
  XX(37) = XX(37)-JVS(357)*XX(81)-JVS(762)*XX(115)-JVS(1163)*XX(136)-JVS(1305)*XX(138)
  XX(36) = XX(36)-JVS(913)*XX(125)-JVS(1304)*XX(138)-JVS(1402)*XX(139)-JVS(1478)*XX(141)
  XX(35) = XX(35)-JVS(964)*XX(127)-JVS(1303)*XX(138)
  XX(34) = XX(34)-JVS(447)*XX(92)-JVS(600)*XX(108)-JVS(1302)*XX(138)-JVS(1477)*XX(141)
  XX(33) = XX(33)-JVS(1096)*XX(134)-JVS(1215)*XX(137)-JVS(1301)*XX(138)-JVS(1476)*XX(141)
  XX(32) = XX(32)-JVS(1017)*XX(132)-JVS(1300)*XX(138)-JVS(1416)*XX(140)-JVS(1475)*XX(141)
  XX(31) = XX(31)-JVS(436)*XX(91)-JVS(1474)*XX(141)
  XX(30) = XX(30)-JVS(682)*XX(110)-JVS(1016)*XX(132)-JVS(1473)*XX(141)
  XX(29) = XX(29)-JVS(761)*XX(115)-JVS(1299)*XX(138)-JVS(1472)*XX(141)
  XX(28) = XX(28)-JVS(230)*XX(58)-JVS(336)*XX(77)-JVS(1298)*XX(138)
  XX(27) = XX(27)-JVS(863)*XX(122)-JVS(1297)*XX(138)-JVS(1471)*XX(141)
  XX(26) = XX(26)-JVS(731)*XX(113)-JVS(1214)*XX(137)
  XX(25) = XX(25)-JVS(270)*XX(65)-JVS(1296)*XX(138)
  XX(24) = XX(24)-JVS(310)*XX(73)-JVS(963)*XX(127)-JVS(1295)*XX(138)-JVS(1470)*XX(141)
  XX(23) = XX(23)-JVS(730)*XX(113)-JVS(1095)*XX(134)-JVS(1294)*XX(138)-JVS(1469)*XX(141)
  XX(22) = XX(22)-JVS(229)*XX(58)-JVS(252)*XX(62)-JVS(1293)*XX(138)-JVS(1468)*XX(141)
  XX(21) = XX(21)-JVS(58)*XX(22)-JVS(190)*XX(52)-JVS(430)*XX(90)-JVS(634)*XX(109)-JVS(1162)*XX(136)-JVS(1292)*XX(138)&
             &-JVS(1467)*XX(141)
  XX(20) = XX(20)-JVS(962)*XX(127)-JVS(1291)*XX(138)
  XX(19) = XX(19)-JVS(633)*XX(109)-JVS(1161)*XX(136)-JVS(1466)*XX(141)
  XX(18) = XX(18)-JVS(201)*XX(54)-JVS(228)*XX(58)-JVS(1290)*XX(138)-JVS(1465)*XX(141)
  XX(17) = XX(17)-JVS(45)*XX(18)-JVS(48)*XX(19)-JVS(189)*XX(52)-JVS(413)*XX(87)-JVS(1289)*XX(138)-JVS(1464)*XX(141)
  XX(16) = XX(16)-JVS(729)*XX(113)-JVS(1094)*XX(134)-JVS(1213)*XX(137)
  XX(15) = XX(15)-JVS(454)*XX(93)-JVS(1288)*XX(138)-JVS(1463)*XX(141)
  XX(14) = XX(14)-JVS(200)*XX(54)-JVS(227)*XX(58)-JVS(1287)*XX(138)-JVS(1462)*XX(141)
  XX(13) = XX(13)-JVS(207)*XX(55)-JVS(632)*XX(109)-JVS(1461)*XX(141)
  XX(12) = XX(12)-JVS(31)*XX(13)-JVS(34)*XX(14)-JVS(246)*XX(61)-JVS(1286)*XX(138)-JVS(1460)*XX(141)
  XX(11) = XX(11)-JVS(372)*XX(83)-JVS(531)*XX(102)-JVS(631)*XX(109)-JVS(1285)*XX(138)-JVS(1459)*XX(141)
  XX(10) = XX(10)-JVS(269)*XX(65)-JVS(1212)*XX(137)
  XX(9) = XX(9)-JVS(961)*XX(127)-JVS(1284)*XX(138)
  XX(8) = XX(8)-JVS(89)*XX(30)-JVS(1283)*XX(138)
  XX(7) = XX(7)-JVS(1282)*XX(138)
  XX(6) = XX(6)-JVS(453)*XX(93)-JVS(1572)*XX(142)
  XX(5) = XX(5)-JVS(1281)*XX(138)-JVS(1458)*XX(141)
  XX(4) = XX(4)-JVS(446)*XX(92)-JVS(1280)*XX(138)
  XX(3) = XX(3)-JVS(477)*XX(95)-JVS(1279)*XX(138)
  XX(2) = XX(2)-JVS(1278)*XX(138)
  XX(1) = XX(1)
      
END SUBROUTINE KppSolveTR

! End of KppSolveTR function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE t1_mozcart_LinearAlgebra

