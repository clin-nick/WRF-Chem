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
! File                 : saprc99_LinearAlgebra.f90
! Time                 : Tue Apr 12 23:43:50 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/saprc99
! Equation file        : saprc99.kpp
! Output root filename : saprc99
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE saprc99_LinearAlgebra

  USE saprc99_Parameters
  USE saprc99_JacobianSP

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
  XX(3) = X(3)/JVS(8)
  XX(4) = X(4)/JVS(23)
  XX(5) = X(5)/JVS(33)
  XX(6) = X(6)/JVS(47)
  XX(7) = X(7)/JVS(50)
  XX(8) = X(8)/JVS(55)
  XX(9) = X(9)/JVS(60)
  XX(10) = (X(10)-JVS(2)*XX(1))/(JVS(69))
  XX(11) = X(11)/JVS(71)
  XX(12) = X(12)/JVS(73)
  XX(13) = X(13)/JVS(75)
  XX(14) = X(14)/JVS(77)
  XX(15) = X(15)/JVS(82)
  XX(16) = X(16)/JVS(85)
  XX(17) = X(17)/JVS(88)
  XX(18) = X(18)/JVS(91)
  XX(19) = X(19)/JVS(94)
  XX(20) = X(20)/JVS(97)
  XX(21) = X(21)/JVS(99)
  XX(22) = X(22)/JVS(101)
  XX(23) = X(23)/JVS(104)
  XX(24) = X(24)/JVS(107)
  XX(25) = X(25)/JVS(110)
  XX(26) = X(26)/JVS(113)
  XX(27) = (X(27)-JVS(78)*XX(14))/(JVS(115))
  XX(28) = X(28)/JVS(117)
  XX(29) = X(29)/JVS(122)
  XX(30) = (X(30)-JVS(9)*XX(3))/(JVS(126))
  XX(31) = (X(31)-JVS(56)*XX(8)-JVS(61)*XX(9))/(JVS(130))
  XX(32) = X(32)/JVS(135)
  XX(33) = X(33)/JVS(139)
  XX(34) = X(34)/JVS(141)
  XX(35) = X(35)/JVS(145)
  XX(36) = X(36)/JVS(149)
  XX(37) = X(37)/JVS(153)
  XX(38) = X(38)/JVS(158)
  XX(39) = (X(39)-JVS(5)*XX(2)-JVS(10)*XX(3))/(JVS(161))
  XX(40) = (X(40)-JVS(131)*XX(31))/(JVS(165))
  XX(41) = (X(41)-JVS(62)*XX(9))/(JVS(172))
  XX(42) = X(42)/JVS(177)
  XX(43) = X(43)/JVS(183)
  XX(44) = X(44)/JVS(193)
  XX(45) = X(45)/JVS(206)
  XX(46) = X(46)/JVS(229)
  XX(47) = X(47)/JVS(244)
  XX(48) = (X(48)-JVS(11)*XX(3)-JVS(207)*XX(45))/(JVS(259))
  XX(49) = (X(49)-JVS(208)*XX(45)-JVS(245)*XX(47))/(JVS(272))
  XX(50) = (X(50)-JVS(6)*XX(2)-JVS(12)*XX(3)-JVS(24)*XX(4)-JVS(57)*XX(8)-JVS(63)*XX(9)-JVS(209)*XX(45))/(JVS(281))
  XX(51) = (X(51)-JVS(166)*XX(40)-JVS(184)*XX(43)-JVS(246)*XX(47)-JVS(273)*XX(49))/(JVS(288))
  XX(52) = (X(52)-JVS(13)*XX(3)-JVS(34)*XX(5)-JVS(210)*XX(45))/(JVS(298))
  XX(53) = X(53)/JVS(310)
  XX(54) = (X(54)-JVS(14)*XX(3)-JVS(35)*XX(5)-JVS(79)*XX(14)-JVS(211)*XX(45)-JVS(230)*XX(46)-JVS(274)*XX(49)-JVS(311)&
             &*XX(53))/(JVS(332))
  XX(55) = (X(55)-JVS(15)*XX(3)-JVS(194)*XX(44)-JVS(212)*XX(45)-JVS(247)*XX(47))/(JVS(338))
  XX(56) = (X(56)-JVS(16)*XX(3)-JVS(25)*XX(4)-JVS(36)*XX(5)-JVS(213)*XX(45)-JVS(231)*XX(46)-JVS(312)*XX(53))/(JVS(344))
  XX(57) = (X(57)-JVS(17)*XX(3)-JVS(37)*XX(5)-JVS(195)*XX(44)-JVS(214)*XX(45)-JVS(248)*XX(47)-JVS(275)*XX(49))&
             &/(JVS(350))
  XX(58) = (X(58)-JVS(18)*XX(3)-JVS(26)*XX(4)-JVS(38)*XX(5)-JVS(173)*XX(41)-JVS(215)*XX(45)-JVS(232)*XX(46)-JVS(249)&
             &*XX(47)-JVS(313)*XX(53)-JVS(339)*XX(55)-JVS(351)*XX(57))/(JVS(356))
  XX(59) = (X(59)-JVS(19)*XX(3)-JVS(196)*XX(44)-JVS(216)*XX(45)-JVS(250)*XX(47)-JVS(314)*XX(53))/(JVS(363))
  XX(60) = (X(60)-JVS(217)*XX(45)-JVS(251)*XX(47))/(JVS(378))
  XX(61) = (X(61)-JVS(127)*XX(30)-JVS(218)*XX(45)-JVS(252)*XX(47))/(JVS(412))
  XX(62) = (X(62)-JVS(233)*XX(46)-JVS(315)*XX(53)-JVS(379)*XX(60)-JVS(413)*XX(61))/(JVS(435))
  XX(63) = (X(63)-JVS(64)*XX(9)-JVS(219)*XX(45)-JVS(260)*XX(48)-JVS(276)*XX(49)-JVS(282)*XX(50)-JVS(299)*XX(52)-JVS(316)&
             &*XX(53)-JVS(333)*XX(54)-JVS(340)*XX(55)-JVS(345)*XX(56)-JVS(352)*XX(57)-JVS(357)*XX(58)-JVS(364)*XX(59)&
             &-JVS(380)*XX(60)-JVS(414)*XX(61)-JVS(436)*XX(62))/(JVS(452))
  XX(64) = (X(64)-JVS(220)*XX(45)-JVS(253)*XX(47)-JVS(381)*XX(60))/(JVS(476))
  XX(65) = (X(65)-JVS(317)*XX(53)-JVS(382)*XX(60)-JVS(415)*XX(61)-JVS(477)*XX(64))/(JVS(498))
  XX(66) = (X(66)-JVS(318)*XX(53)-JVS(383)*XX(60)-JVS(416)*XX(61)-JVS(478)*XX(64)-JVS(499)*XX(65))/(JVS(519))
  XX(67) = (X(67)-JVS(7)*XX(2)-JVS(20)*XX(3)-JVS(27)*XX(4)-JVS(39)*XX(5)-JVS(65)*XX(9)-JVS(72)*XX(11)-JVS(80)*XX(14)&
             &-JVS(159)*XX(38)-JVS(162)*XX(39)-JVS(174)*XX(41)-JVS(197)*XX(44)-JVS(221)*XX(45)-JVS(234)*XX(46)-JVS(254)&
             &*XX(47)-JVS(261)*XX(48)-JVS(277)*XX(49)-JVS(283)*XX(50)-JVS(300)*XX(52)-JVS(319)*XX(53)-JVS(334)*XX(54)&
             &-JVS(341)*XX(55)-JVS(346)*XX(56)-JVS(353)*XX(57)-JVS(358)*XX(58)-JVS(365)*XX(59)-JVS(384)*XX(60)-JVS(417)&
             &*XX(61)-JVS(437)*XX(62)-JVS(453)*XX(63)-JVS(479)*XX(64)-JVS(500)*XX(65)-JVS(520)*XX(66))/(JVS(543))
  XX(68) = (X(68)-JVS(28)*XX(4)-JVS(48)*XX(6)-JVS(83)*XX(15)-JVS(289)*XX(51)-JVS(320)*XX(53)-JVS(385)*XX(60)-JVS(418)&
             &*XX(61)-JVS(521)*XX(66)-JVS(544)*XX(67))/(JVS(574))
  XX(69) = (X(69)-JVS(29)*XX(4)-JVS(40)*XX(5)-JVS(118)*XX(28)-JVS(178)*XX(42)-JVS(321)*XX(53)-JVS(419)*XX(61)-JVS(438)&
             &*XX(62)-JVS(480)*XX(64)-JVS(501)*XX(65)-JVS(522)*XX(66)-JVS(575)*XX(68))/(JVS(606))
  XX(70) = (X(70)-JVS(41)*XX(5)-JVS(51)*XX(7)-JVS(92)*XX(18)-JVS(290)*XX(51)-JVS(322)*XX(53)-JVS(386)*XX(60)-JVS(420)&
             &*XX(61)-JVS(545)*XX(67)-JVS(576)*XX(68)-JVS(607)*XX(69))/(JVS(627))
  XX(71) = (X(71)-JVS(3)*XX(1)-JVS(21)*XX(3)-JVS(66)*XX(9)-JVS(70)*XX(10)-JVS(74)*XX(12)-JVS(76)*XX(13)-JVS(81)*XX(14)&
             &-JVS(95)*XX(19)-JVS(98)*XX(20)-JVS(100)*XX(21)-JVS(105)*XX(23)-JVS(108)*XX(24)-JVS(111)*XX(25)-JVS(114)*XX(26)&
             &-JVS(116)*XX(27)-JVS(119)*XX(28)-JVS(123)*XX(29)-JVS(136)*XX(32)-JVS(140)*XX(33)-JVS(142)*XX(34)-JVS(146)&
             &*XX(35)-JVS(150)*XX(36)-JVS(154)*XX(37)-JVS(160)*XX(38)-JVS(163)*XX(39)-JVS(175)*XX(41)-JVS(179)*XX(42)&
             &-JVS(185)*XX(43)-JVS(198)*XX(44)-JVS(222)*XX(45)-JVS(235)*XX(46)-JVS(255)*XX(47)-JVS(262)*XX(48)-JVS(278)&
             &*XX(49)-JVS(284)*XX(50)-JVS(291)*XX(51)-JVS(301)*XX(52)-JVS(323)*XX(53)-JVS(335)*XX(54)-JVS(342)*XX(55)&
             &-JVS(347)*XX(56)-JVS(354)*XX(57)-JVS(359)*XX(58)-JVS(366)*XX(59)-JVS(387)*XX(60)-JVS(421)*XX(61)-JVS(439)&
             &*XX(62)-JVS(454)*XX(63)-JVS(481)*XX(64)-JVS(502)*XX(65)-JVS(523)*XX(66)-JVS(546)*XX(67)-JVS(577)*XX(68)&
             &-JVS(608)*XX(69)-JVS(628)*XX(70))/(JVS(687))
  XX(72) = (X(72)-JVS(42)*XX(5)-JVS(52)*XX(7)-JVS(86)*XX(16)-JVS(292)*XX(51)-JVS(324)*XX(53)-JVS(388)*XX(60)-JVS(422)&
             &*XX(61)-JVS(524)*XX(66)-JVS(547)*XX(67)-JVS(578)*XX(68)-JVS(609)*XX(69)-JVS(629)*XX(70)-JVS(688)*XX(71))&
             &/(JVS(713))
  XX(73) = (X(73)-JVS(22)*XX(3)-JVS(106)*XX(23)-JVS(128)*XX(30)-JVS(293)*XX(51)-JVS(325)*XX(53)-JVS(389)*XX(60)-JVS(423)&
             &*XX(61)-JVS(440)*XX(62)-JVS(455)*XX(63)-JVS(482)*XX(64)-JVS(503)*XX(65)-JVS(525)*XX(66)-JVS(548)*XX(67)&
             &-JVS(579)*XX(68)-JVS(610)*XX(69)-JVS(630)*XX(70)-JVS(689)*XX(71)-JVS(714)*XX(72))/(JVS(739))
  XX(74) = (X(74)-JVS(58)*XX(8)-JVS(67)*XX(9)-JVS(84)*XX(15)-JVS(87)*XX(16)-JVS(89)*XX(17)-JVS(93)*XX(18)-JVS(102)&
             &*XX(22)-JVS(112)*XX(25)-JVS(132)*XX(31)-JVS(137)*XX(32)-JVS(167)*XX(40)-JVS(236)*XX(46)-JVS(256)*XX(47)&
             &-JVS(294)*XX(51)-JVS(326)*XX(53)-JVS(424)*XX(61)-JVS(441)*XX(62)-JVS(456)*XX(63)-JVS(483)*XX(64)-JVS(504)&
             &*XX(65)-JVS(526)*XX(66)-JVS(549)*XX(67)-JVS(580)*XX(68)-JVS(611)*XX(69)-JVS(631)*XX(70)-JVS(690)*XX(71)&
             &-JVS(715)*XX(72)-JVS(740)*XX(73))/(JVS(782))
  XX(75) = (X(75)-JVS(59)*XX(8)-JVS(68)*XX(9)-JVS(103)*XX(22)-JVS(133)*XX(31)-JVS(155)*XX(37)-JVS(168)*XX(40)-JVS(176)&
             &*XX(41)-JVS(186)*XX(43)-JVS(199)*XX(44)-JVS(223)*XX(45)-JVS(237)*XX(46)-JVS(257)*XX(47)-JVS(263)*XX(48)&
             &-JVS(279)*XX(49)-JVS(285)*XX(50)-JVS(295)*XX(51)-JVS(302)*XX(52)-JVS(327)*XX(53)-JVS(336)*XX(54)-JVS(343)&
             &*XX(55)-JVS(348)*XX(56)-JVS(355)*XX(57)-JVS(360)*XX(58)-JVS(367)*XX(59)-JVS(390)*XX(60)-JVS(425)*XX(61)&
             &-JVS(442)*XX(62)-JVS(457)*XX(63)-JVS(484)*XX(64)-JVS(505)*XX(65)-JVS(527)*XX(66)-JVS(550)*XX(67)-JVS(581)&
             &*XX(68)-JVS(612)*XX(69)-JVS(632)*XX(70)-JVS(691)*XX(71)-JVS(716)*XX(72)-JVS(741)*XX(73)-JVS(783)*XX(74))&
             &/(JVS(823))
  XX(76) = (X(76)-JVS(30)*XX(4)-JVS(43)*XX(5)-JVS(120)*XX(28)-JVS(124)*XX(29)-JVS(328)*XX(53)-JVS(426)*XX(61)-JVS(506)&
             &*XX(65)-JVS(528)*XX(66)-JVS(582)*XX(68)-JVS(613)*XX(69)-JVS(633)*XX(70)-JVS(692)*XX(71)-JVS(717)*XX(72)&
             &-JVS(742)*XX(73)-JVS(784)*XX(74)-JVS(824)*XX(75))/(JVS(855))
  XX(77) = (X(77)-JVS(31)*XX(4)-JVS(44)*XX(5)-JVS(121)*XX(28)-JVS(180)*XX(42)-JVS(329)*XX(53)-JVS(427)*XX(61)-JVS(485)&
             &*XX(64)-JVS(507)*XX(65)-JVS(529)*XX(66)-JVS(583)*XX(68)-JVS(614)*XX(69)-JVS(634)*XX(70)-JVS(693)*XX(71)&
             &-JVS(718)*XX(72)-JVS(743)*XX(73)-JVS(785)*XX(74)-JVS(825)*XX(75)-JVS(856)*XX(76))/(JVS(899))
  XX(78) = (X(78)-JVS(32)*XX(4)-JVS(45)*XX(5)-JVS(49)*XX(6)-JVS(53)*XX(7)-JVS(96)*XX(19)-JVS(125)*XX(29)-JVS(129)*XX(30)&
             &-JVS(134)*XX(31)-JVS(138)*XX(32)-JVS(169)*XX(40)-JVS(181)*XX(42)-JVS(187)*XX(43)-JVS(258)*XX(47)-JVS(280)&
             &*XX(49)-JVS(296)*XX(51)-JVS(330)*XX(53)-JVS(428)*XX(61)-JVS(486)*XX(64)-JVS(551)*XX(67)-JVS(584)*XX(68)&
             &-JVS(615)*XX(69)-JVS(635)*XX(70)-JVS(694)*XX(71)-JVS(719)*XX(72)-JVS(744)*XX(73)-JVS(786)*XX(74)-JVS(826)&
             &*XX(75)-JVS(857)*XX(76)-JVS(900)*XX(77))/(JVS(950))
  XX(79) = (X(79)-JVS(46)*XX(5)-JVS(54)*XX(7)-JVS(90)*XX(17)-JVS(297)*XX(51)-JVS(331)*XX(53)-JVS(391)*XX(60)-JVS(429)&
             &*XX(61)-JVS(530)*XX(66)-JVS(552)*XX(67)-JVS(585)*XX(68)-JVS(616)*XX(69)-JVS(636)*XX(70)-JVS(695)*XX(71)&
             &-JVS(720)*XX(72)-JVS(745)*XX(73)-JVS(787)*XX(74)-JVS(827)*XX(75)-JVS(858)*XX(76)-JVS(901)*XX(77)-JVS(951)&
             &*XX(78))/(JVS(968))
  XX(79) = XX(79)
  XX(78) = XX(78)-JVS(967)*XX(79)
  XX(77) = XX(77)-JVS(949)*XX(78)-JVS(966)*XX(79)
  XX(76) = XX(76)-JVS(898)*XX(77)-JVS(948)*XX(78)-JVS(965)*XX(79)
  XX(75) = XX(75)-JVS(854)*XX(76)-JVS(897)*XX(77)-JVS(947)*XX(78)-JVS(964)*XX(79)
  XX(74) = XX(74)-JVS(822)*XX(75)-JVS(853)*XX(76)-JVS(896)*XX(77)-JVS(946)*XX(78)-JVS(963)*XX(79)
  XX(73) = XX(73)-JVS(781)*XX(74)-JVS(821)*XX(75)-JVS(852)*XX(76)-JVS(895)*XX(77)-JVS(945)*XX(78)-JVS(962)*XX(79)
  XX(72) = XX(72)-JVS(738)*XX(73)-JVS(780)*XX(74)-JVS(820)*XX(75)-JVS(851)*XX(76)-JVS(894)*XX(77)-JVS(944)*XX(78)&
             &-JVS(961)*XX(79)
  XX(71) = XX(71)-JVS(712)*XX(72)-JVS(737)*XX(73)-JVS(779)*XX(74)-JVS(819)*XX(75)-JVS(850)*XX(76)-JVS(893)*XX(77)&
             &-JVS(943)*XX(78)-JVS(960)*XX(79)
  XX(70) = XX(70)-JVS(686)*XX(71)-JVS(711)*XX(72)-JVS(736)*XX(73)-JVS(778)*XX(74)-JVS(818)*XX(75)-JVS(849)*XX(76)&
             &-JVS(892)*XX(77)-JVS(942)*XX(78)-JVS(959)*XX(79)
  XX(69) = XX(69)-JVS(626)*XX(70)-JVS(685)*XX(71)-JVS(710)*XX(72)-JVS(735)*XX(73)-JVS(777)*XX(74)-JVS(817)*XX(75)&
             &-JVS(848)*XX(76)-JVS(891)*XX(77)-JVS(941)*XX(78)-JVS(958)*XX(79)
  XX(68) = XX(68)-JVS(605)*XX(69)-JVS(625)*XX(70)-JVS(684)*XX(71)-JVS(709)*XX(72)-JVS(734)*XX(73)-JVS(776)*XX(74)&
             &-JVS(816)*XX(75)-JVS(847)*XX(76)-JVS(890)*XX(77)-JVS(940)*XX(78)-JVS(957)*XX(79)
  XX(67) = XX(67)-JVS(573)*XX(68)-JVS(604)*XX(69)-JVS(624)*XX(70)-JVS(683)*XX(71)-JVS(708)*XX(72)-JVS(733)*XX(73)&
             &-JVS(775)*XX(74)-JVS(815)*XX(75)-JVS(846)*XX(76)-JVS(889)*XX(77)-JVS(939)*XX(78)-JVS(956)*XX(79)
  XX(66) = XX(66)-JVS(572)*XX(68)-JVS(603)*XX(69)-JVS(682)*XX(71)-JVS(707)*XX(72)-JVS(732)*XX(73)-JVS(774)*XX(74)&
             &-JVS(814)*XX(75)-JVS(845)*XX(76)-JVS(888)*XX(77)-JVS(938)*XX(78)
  XX(65) = XX(65)-JVS(571)*XX(68)-JVS(602)*XX(69)-JVS(681)*XX(71)-JVS(706)*XX(72)-JVS(731)*XX(73)-JVS(773)*XX(74)&
             &-JVS(813)*XX(75)-JVS(844)*XX(76)-JVS(887)*XX(77)-JVS(937)*XX(78)
  XX(64) = XX(64)-JVS(570)*XX(68)-JVS(601)*XX(69)-JVS(680)*XX(71)-JVS(705)*XX(72)-JVS(772)*XX(74)-JVS(812)*XX(75)&
             &-JVS(843)*XX(76)-JVS(886)*XX(77)-JVS(936)*XX(78)
  XX(63) = XX(63)-JVS(475)*XX(64)-JVS(497)*XX(65)-JVS(518)*XX(66)-JVS(542)*XX(67)-JVS(569)*XX(68)-JVS(600)*XX(69)&
             &-JVS(623)*XX(70)-JVS(679)*XX(71)-JVS(704)*XX(72)-JVS(730)*XX(73)-JVS(771)*XX(74)-JVS(811)*XX(75)-JVS(842)&
             &*XX(76)-JVS(885)*XX(77)-JVS(935)*XX(78)-JVS(955)*XX(79)
  XX(62) = XX(62)-JVS(474)*XX(64)-JVS(496)*XX(65)-JVS(517)*XX(66)-JVS(568)*XX(68)-JVS(599)*XX(69)-JVS(678)*XX(71)&
             &-JVS(729)*XX(73)-JVS(770)*XX(74)-JVS(810)*XX(75)-JVS(841)*XX(76)-JVS(884)*XX(77)-JVS(934)*XX(78)
  XX(61) = XX(61)-JVS(677)*XX(71)-JVS(728)*XX(73)-JVS(769)*XX(74)-JVS(809)*XX(75)-JVS(933)*XX(78)
  XX(60) = XX(60)-JVS(567)*XX(68)-JVS(676)*XX(71)-JVS(768)*XX(74)-JVS(808)*XX(75)-JVS(840)*XX(76)-JVS(932)*XX(78)
  XX(59) = XX(59)-JVS(411)*XX(61)-JVS(451)*XX(63)-JVS(473)*XX(64)-JVS(495)*XX(65)-JVS(516)*XX(66)-JVS(541)*XX(67)&
             &-JVS(566)*XX(68)-JVS(598)*XX(69)-JVS(622)*XX(70)-JVS(675)*XX(71)-JVS(703)*XX(72)-JVS(727)*XX(73)-JVS(767)&
             &*XX(74)-JVS(807)*XX(75)-JVS(839)*XX(76)-JVS(883)*XX(77)-JVS(931)*XX(78)
  XX(58) = XX(58)-JVS(362)*XX(59)-JVS(377)*XX(60)-JVS(410)*XX(61)-JVS(434)*XX(62)-JVS(450)*XX(63)-JVS(472)*XX(64)&
             &-JVS(494)*XX(65)-JVS(515)*XX(66)-JVS(540)*XX(67)-JVS(565)*XX(68)-JVS(597)*XX(69)-JVS(621)*XX(70)-JVS(674)&
             &*XX(71)-JVS(702)*XX(72)-JVS(726)*XX(73)-JVS(766)*XX(74)-JVS(806)*XX(75)-JVS(838)*XX(76)-JVS(882)*XX(77)&
             &-JVS(930)*XX(78)-JVS(954)*XX(79)
  XX(57) = XX(57)-JVS(376)*XX(60)-JVS(409)*XX(61)-JVS(433)*XX(62)-JVS(471)*XX(64)-JVS(493)*XX(65)-JVS(514)*XX(66)&
             &-JVS(539)*XX(67)-JVS(564)*XX(68)-JVS(596)*XX(69)-JVS(620)*XX(70)-JVS(673)*XX(71)-JVS(701)*XX(72)-JVS(765)&
             &*XX(74)-JVS(805)*XX(75)-JVS(881)*XX(77)-JVS(929)*XX(78)
  XX(56) = XX(56)-JVS(375)*XX(60)-JVS(408)*XX(61)-JVS(432)*XX(62)-JVS(449)*XX(63)-JVS(470)*XX(64)-JVS(492)*XX(65)&
             &-JVS(513)*XX(66)-JVS(538)*XX(67)-JVS(563)*XX(68)-JVS(595)*XX(69)-JVS(672)*XX(71)-JVS(725)*XX(73)-JVS(764)&
             &*XX(74)-JVS(804)*XX(75)-JVS(837)*XX(76)-JVS(880)*XX(77)-JVS(928)*XX(78)
  XX(55) = XX(55)-JVS(407)*XX(61)-JVS(448)*XX(63)-JVS(469)*XX(64)-JVS(491)*XX(65)-JVS(537)*XX(67)-JVS(562)*XX(68)&
             &-JVS(619)*XX(70)-JVS(671)*XX(71)-JVS(700)*XX(72)-JVS(763)*XX(74)-JVS(803)*XX(75)-JVS(879)*XX(77)-JVS(927)&
             &*XX(78)
  XX(54) = XX(54)-JVS(406)*XX(61)-JVS(431)*XX(62)-JVS(447)*XX(63)-JVS(468)*XX(64)-JVS(512)*XX(66)-JVS(536)*XX(67)&
             &-JVS(561)*XX(68)-JVS(594)*XX(69)-JVS(670)*XX(71)-JVS(699)*XX(72)-JVS(724)*XX(73)-JVS(762)*XX(74)-JVS(802)&
             &*XX(75)-JVS(836)*XX(76)-JVS(878)*XX(77)-JVS(926)*XX(78)
  XX(53) = XX(53)-JVS(723)*XX(73)-JVS(761)*XX(74)-JVS(801)*XX(75)
  XX(52) = XX(52)-JVS(309)*XX(53)-JVS(337)*XX(55)-JVS(349)*XX(57)-JVS(361)*XX(59)-JVS(405)*XX(61)-JVS(446)*XX(63)&
             &-JVS(511)*XX(66)-JVS(535)*XX(67)-JVS(593)*XX(69)-JVS(618)*XX(70)-JVS(669)*XX(71)-JVS(760)*XX(74)-JVS(800)&
             &*XX(75)-JVS(835)*XX(76)-JVS(877)*XX(77)-JVS(925)*XX(78)
  XX(51) = XX(51)-JVS(404)*XX(61)-JVS(668)*XX(71)-JVS(698)*XX(72)-JVS(759)*XX(74)-JVS(799)*XX(75)-JVS(876)*XX(77)&
             &-JVS(924)*XX(78)
  XX(50) = XX(50)-JVS(374)*XX(60)-JVS(403)*XX(61)-JVS(445)*XX(63)-JVS(467)*XX(64)-JVS(490)*XX(65)-JVS(534)*XX(67)&
             &-JVS(592)*XX(69)-JVS(667)*XX(71)-JVS(798)*XX(75)-JVS(834)*XX(76)-JVS(875)*XX(77)-JVS(923)*XX(78)
  XX(49) = XX(49)-JVS(402)*XX(61)-JVS(666)*XX(71)-JVS(697)*XX(72)-JVS(758)*XX(74)-JVS(797)*XX(75)-JVS(922)*XX(78)
  XX(48) = XX(48)-JVS(271)*XX(49)-JVS(373)*XX(60)-JVS(401)*XX(61)-JVS(444)*XX(63)-JVS(466)*XX(64)-JVS(533)*XX(67)&
             &-JVS(665)*XX(71)-JVS(796)*XX(75)-JVS(833)*XX(76)-JVS(874)*XX(77)-JVS(921)*XX(78)
  XX(47) = XX(47)-JVS(664)*XX(71)-JVS(757)*XX(74)-JVS(795)*XX(75)
  XX(46) = XX(46)-JVS(308)*XX(53)-JVS(400)*XX(61)-JVS(560)*XX(68)-JVS(663)*XX(71)-JVS(832)*XX(76)
  XX(45) = XX(45)-JVS(662)*XX(71)-JVS(920)*XX(78)
  XX(44) = XX(44)-JVS(205)*XX(45)-JVS(243)*XX(47)-JVS(559)*XX(68)-JVS(661)*XX(71)-JVS(794)*XX(75)-JVS(919)*XX(78)
  XX(43) = XX(43)-JVS(242)*XX(47)-JVS(270)*XX(49)-JVS(287)*XX(51)-JVS(660)*XX(71)-JVS(793)*XX(75)-JVS(873)*XX(77)
  XX(42) = XX(42)-JVS(465)*XX(64)-JVS(659)*XX(71)-JVS(872)*XX(77)-JVS(918)*XX(78)
  XX(41) = XX(41)-JVS(241)*XX(47)-JVS(658)*XX(71)-JVS(792)*XX(75)-JVS(953)*XX(79)
  XX(40) = XX(40)-JVS(240)*XX(47)-JVS(756)*XX(74)-JVS(791)*XX(75)-JVS(917)*XX(78)
  XX(39) = XX(39)-JVS(204)*XX(45)-JVS(269)*XX(49)-JVS(399)*XX(61)-JVS(532)*XX(67)-JVS(657)*XX(71)-JVS(871)*XX(77)&
             &-JVS(916)*XX(78)
  XX(38) = XX(38)-JVS(203)*XX(45)-JVS(268)*XX(49)-JVS(464)*XX(64)-JVS(531)*XX(67)-JVS(656)*XX(71)-JVS(870)*XX(77)&
             &-JVS(915)*XX(78)
  XX(37) = XX(37)-JVS(192)*XX(44)-JVS(239)*XX(47)-JVS(286)*XX(51)-JVS(655)*XX(71)-JVS(790)*XX(75)-JVS(869)*XX(77)
  XX(36) = XX(36)-JVS(191)*XX(44)-JVS(202)*XX(45)-JVS(267)*XX(49)-JVS(307)*XX(53)-JVS(463)*XX(64)-JVS(558)*XX(68)&
             &-JVS(654)*XX(71)-JVS(868)*XX(77)-JVS(914)*XX(78)
  XX(35) = XX(35)-JVS(190)*XX(44)-JVS(201)*XX(45)-JVS(266)*XX(49)-JVS(306)*XX(53)-JVS(462)*XX(64)-JVS(557)*XX(68)&
             &-JVS(653)*XX(71)-JVS(867)*XX(77)-JVS(913)*XX(78)
  XX(34) = XX(34)-JVS(144)*XX(35)-JVS(148)*XX(36)-JVS(152)*XX(37)-JVS(157)*XX(38)-JVS(171)*XX(41)-JVS(182)*XX(43)&
             &-JVS(189)*XX(44)-JVS(265)*XX(49)-JVS(510)*XX(66)-JVS(591)*XX(69)-JVS(652)*XX(71)-JVS(866)*XX(77)-JVS(912)&
             &*XX(78)
  XX(33) = XX(33)-JVS(200)*XX(45)-JVS(228)*XX(46)-JVS(305)*XX(53)-JVS(372)*XX(60)-JVS(398)*XX(61)-JVS(461)*XX(64)&
             &-JVS(489)*XX(65)-JVS(509)*XX(66)-JVS(556)*XX(68)-JVS(590)*XX(69)-JVS(651)*XX(71)-JVS(831)*XX(76)-JVS(865)&
             &*XX(77)
  XX(32) = XX(32)-JVS(650)*XX(71)-JVS(755)*XX(74)-JVS(789)*XX(75)-JVS(911)*XX(78)
  XX(31) = XX(31)-JVS(164)*XX(40)-JVS(754)*XX(74)-JVS(910)*XX(78)
  XX(30) = XX(30)-JVS(397)*XX(61)-JVS(722)*XX(73)-JVS(753)*XX(74)-JVS(909)*XX(78)
  XX(29) = XX(29)-JVS(396)*XX(61)-JVS(649)*XX(71)-JVS(830)*XX(76)-JVS(908)*XX(78)
  XX(28) = XX(28)-JVS(395)*XX(61)-JVS(648)*XX(71)-JVS(907)*XX(78)
  XX(27) = XX(27)-JVS(143)*XX(35)-JVS(147)*XX(36)-JVS(151)*XX(37)-JVS(156)*XX(38)-JVS(170)*XX(41)-JVS(188)*XX(44)&
             &-JVS(264)*XX(49)-JVS(555)*XX(68)-JVS(589)*XX(69)-JVS(647)*XX(71)-JVS(864)*XX(77)-JVS(906)*XX(78)
  XX(26) = XX(26)-JVS(227)*XX(46)-JVS(304)*XX(53)-JVS(371)*XX(60)-JVS(394)*XX(61)-JVS(460)*XX(64)-JVS(488)*XX(65)&
             &-JVS(508)*XX(66)-JVS(588)*XX(69)-JVS(646)*XX(71)-JVS(863)*XX(77)
  XX(25) = XX(25)-JVS(226)*XX(46)-JVS(430)*XX(62)-JVS(752)*XX(74)-JVS(829)*XX(76)
  XX(24) = XX(24)-JVS(109)*XX(25)-JVS(225)*XX(46)-JVS(303)*XX(53)-JVS(370)*XX(60)-JVS(393)*XX(61)-JVS(459)*XX(64)&
             &-JVS(487)*XX(65)-JVS(587)*XX(69)-JVS(645)*XX(71)-JVS(862)*XX(77)
  XX(23) = XX(23)-JVS(644)*XX(71)-JVS(721)*XX(73)-JVS(751)*XX(74)-JVS(905)*XX(78)
  XX(22) = XX(22)-JVS(238)*XX(47)-JVS(750)*XX(74)-JVS(788)*XX(75)
  XX(21) = XX(21)-JVS(369)*XX(60)-JVS(392)*XX(61)-JVS(643)*XX(71)-JVS(861)*XX(77)-JVS(904)*XX(78)
  XX(20) = XX(20)-JVS(224)*XX(46)-JVS(458)*XX(64)-JVS(586)*XX(69)-JVS(642)*XX(71)-JVS(860)*XX(77)
  XX(19) = XX(19)-JVS(641)*XX(71)-JVS(903)*XX(78)
  XX(18) = XX(18)-JVS(617)*XX(70)-JVS(749)*XX(74)
  XX(17) = XX(17)-JVS(748)*XX(74)-JVS(952)*XX(79)
  XX(16) = XX(16)-JVS(696)*XX(72)-JVS(747)*XX(74)
  XX(15) = XX(15)-JVS(554)*XX(68)-JVS(746)*XX(74)
  XX(14) = XX(14)-JVS(553)*XX(68)
  XX(13) = XX(13)-JVS(368)*XX(60)-JVS(640)*XX(71)-JVS(859)*XX(77)
  XX(12) = XX(12)-JVS(639)*XX(71)-JVS(828)*XX(76)
  XX(11) = XX(11)-JVS(443)*XX(63)-JVS(638)*XX(71)
  XX(10) = XX(10)-JVS(637)*XX(71)-JVS(902)*XX(78)
  XX(9) = XX(9)
  XX(8) = XX(8)
  XX(7) = XX(7)
  XX(6) = XX(6)
  XX(5) = XX(5)
  XX(4) = XX(4)
  XX(3) = XX(3)
  XX(2) = XX(2)
  XX(1) = XX(1)
      
END SUBROUTINE KppSolveTR

! End of KppSolveTR function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE saprc99_LinearAlgebra

