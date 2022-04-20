! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Jacobian of Chemical Model File
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
! File                 : mozart_Jacobian.f90
! Time                 : Thu Apr 14 13:13:51 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/mozart
! Equation file        : mozart.kpp
! Output root filename : mozart
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mozart_Jacobian

  USE mozart_Parameters
  USE mozart_JacobianSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Jac_SP_Vec - function for sparse multiplication: sparse Jacobian times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JUV       - Jacobian times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Jac_SP_Vec ( JVS, UV, JUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JUV - Jacobian times user vector
  REAL(kind=dp) :: JUV(NVAR)

  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(8)+JVS(3)*UV(75)
  JUV(2) = JVS(4)*UV(2)+JVS(5)*UV(75)
  JUV(3) = JVS(6)*UV(3)+JVS(7)*UV(47)
  JUV(4) = JVS(8)*UV(4)+JVS(9)*UV(75)
  JUV(5) = JVS(10)*UV(5)+JVS(11)*UV(75)
  JUV(6) = JVS(12)*UV(6)+JVS(13)*UV(75)
  JUV(7) = JVS(14)*UV(7)+JVS(15)*UV(75)
  JUV(8) = JVS(16)*UV(8)+JVS(17)*UV(16)+JVS(18)*UV(75)+JVS(19)*UV(80)
  JUV(9) = JVS(20)*UV(9)+JVS(21)*UV(75)+JVS(22)*UV(79)
  JUV(10) = JVS(23)*UV(10)+JVS(24)*UV(75)
  JUV(11) = JVS(25)*UV(10)+JVS(26)*UV(11)+JVS(27)*UV(75)
  JUV(12) = JVS(28)*UV(12)+JVS(29)*UV(30)+JVS(30)*UV(76)
  JUV(13) = JVS(31)*UV(13)+JVS(32)*UV(73)+JVS(33)*UV(80)
  JUV(14) = JVS(34)*UV(11)+JVS(35)*UV(14)+JVS(36)*UV(73)+JVS(37)*UV(75)
  JUV(15) = JVS(38)*UV(15)+JVS(39)*UV(68)+JVS(40)*UV(75)+JVS(41)*UV(79)
  JUV(16) = JVS(42)*UV(16)+JVS(43)*UV(75)+JVS(44)*UV(80)
  JUV(17) = JVS(45)*UV(17)+JVS(46)*UV(55)+JVS(47)*UV(75)+JVS(48)*UV(77)
  JUV(18) = JVS(49)*UV(18)+JVS(50)*UV(75)+JVS(51)*UV(77)+JVS(52)*UV(79)
  JUV(19) = JVS(53)*UV(19)+JVS(54)*UV(40)+JVS(55)*UV(47)+JVS(56)*UV(72)+JVS(57)*UV(75)
  JUV(20) = JVS(58)*UV(20)+JVS(59)*UV(70)+JVS(60)*UV(75)+JVS(61)*UV(79)
  JUV(21) = JVS(62)*UV(21)+JVS(63)*UV(55)+JVS(64)*UV(75)+JVS(65)*UV(79)
  JUV(22) = JVS(66)*UV(22)+JVS(67)*UV(57)+JVS(68)*UV(75)+JVS(69)*UV(79)
  JUV(23) = JVS(70)*UV(23)+JVS(71)*UV(63)+JVS(72)*UV(75)+JVS(73)*UV(79)
  JUV(24) = JVS(74)*UV(6)+JVS(75)*UV(24)+JVS(76)*UV(75)+JVS(77)*UV(76)
  JUV(25) = JVS(78)*UV(25)+JVS(79)*UV(49)+JVS(80)*UV(75)+JVS(81)*UV(79)
  JUV(26) = JVS(82)*UV(26)+JVS(83)*UV(73)+JVS(84)*UV(75)+JVS(85)*UV(79)
  JUV(27) = JVS(86)*UV(27)+JVS(87)*UV(75)+JVS(88)*UV(81)
  JUV(28) = JVS(89)*UV(28)+JVS(90)*UV(60)+JVS(91)*UV(71)+JVS(92)*UV(75)+JVS(93)*UV(76)+JVS(94)*UV(77)+JVS(95)*UV(78)&
              &+JVS(96)*UV(80)
  JUV(29) = JVS(97)*UV(29)+JVS(98)*UV(54)+JVS(99)*UV(61)+JVS(100)*UV(75)+JVS(101)*UV(76)+JVS(102)*UV(80)
  JUV(30) = JVS(103)*UV(27)+JVS(104)*UV(30)+JVS(105)*UV(75)+JVS(106)*UV(76)
  JUV(31) = JVS(108)*UV(31)+JVS(109)*UV(44)+JVS(110)*UV(75)+JVS(111)*UV(79)
  JUV(32) = JVS(112)*UV(32)+JVS(113)*UV(74)+JVS(114)*UV(75)+JVS(115)*UV(78)+JVS(116)*UV(79)
  JUV(33) = JVS(117)*UV(33)+JVS(118)*UV(48)+JVS(119)*UV(54)+JVS(120)*UV(75)+JVS(121)*UV(76)
  JUV(34) = JVS(122)*UV(27)+JVS(123)*UV(34)+JVS(124)*UV(59)+JVS(125)*UV(61)+JVS(126)*UV(74)+JVS(127)*UV(75)+JVS(128)&
              &*UV(77)+JVS(129)*UV(78)+JVS(130)*UV(79)+JVS(131)*UV(81)
  JUV(35) = JVS(132)*UV(35)+JVS(133)*UV(50)+JVS(134)*UV(75)+JVS(135)*UV(79)
  JUV(36) = JVS(136)*UV(36)+JVS(137)*UV(73)+JVS(138)*UV(75)+JVS(139)*UV(78)
  JUV(37) = JVS(140)*UV(31)+JVS(141)*UV(37)+JVS(142)*UV(44)+JVS(143)*UV(45)+JVS(144)*UV(53)+JVS(145)*UV(75)+JVS(146)&
              &*UV(76)
  JUV(38) = JVS(148)*UV(38)+JVS(149)*UV(64)+JVS(150)*UV(75)+JVS(151)*UV(79)
  JUV(39) = JVS(152)*UV(13)+JVS(153)*UV(16)+JVS(154)*UV(39)+JVS(155)*UV(62)+JVS(156)*UV(67)+JVS(157)*UV(72)+JVS(158)&
              &*UV(73)+JVS(159)*UV(75)+JVS(160)*UV(80)
  JUV(40) = JVS(161)*UV(40)+JVS(162)*UV(47)+JVS(163)*UV(61)+JVS(164)*UV(75)+JVS(165)*UV(81)
  JUV(41) = JVS(166)*UV(41)+JVS(167)*UV(55)+JVS(168)*UV(63)+JVS(169)*UV(68)+JVS(170)*UV(70)+JVS(171)*UV(71)+JVS(172)&
              &*UV(75)+JVS(173)*UV(77)
  JUV(42) = JVS(174)*UV(42)+JVS(175)*UV(71)+JVS(176)*UV(75)+JVS(177)*UV(79)
  JUV(43) = JVS(178)*UV(43)+JVS(179)*UV(73)+JVS(180)*UV(74)+JVS(181)*UV(75)
  JUV(44) = JVS(182)*UV(10)+JVS(183)*UV(31)+JVS(184)*UV(44)+JVS(185)*UV(75)+JVS(186)*UV(76)+JVS(187)*UV(79)
  JUV(45) = JVS(188)*UV(14)+JVS(189)*UV(31)+JVS(190)*UV(44)+JVS(191)*UV(45)+JVS(192)*UV(73)+JVS(194)*UV(76)
  JUV(46) = JVS(196)*UV(46)+JVS(197)*UV(75)+JVS(198)*UV(80)+JVS(199)*UV(81)
  JUV(47) = JVS(200)*UV(3)+JVS(201)*UV(19)+JVS(202)*UV(40)+JVS(203)*UV(47)+JVS(207)*UV(81)
  JUV(48) = JVS(208)*UV(48)+JVS(209)*UV(54)+JVS(210)*UV(75)+JVS(211)*UV(79)
  JUV(49) = JVS(212)*UV(25)+JVS(213)*UV(33)+JVS(215)*UV(49)+JVS(217)*UV(75)+JVS(218)*UV(76)+JVS(219)*UV(79)
  JUV(50) = JVS(220)*UV(35)+JVS(221)*UV(50)+JVS(222)*UV(61)+JVS(223)*UV(75)+JVS(224)*UV(76)+JVS(225)*UV(79)
  JUV(51) = JVS(226)*UV(22)+JVS(227)*UV(24)+JVS(228)*UV(38)+JVS(229)*UV(48)+JVS(230)*UV(51)+JVS(231)*UV(54)+JVS(232)&
              &*UV(57)+JVS(233)*UV(64)+JVS(234)*UV(75)+JVS(235)*UV(76)+JVS(236)*UV(77)
  JUV(52) = JVS(238)*UV(47)+JVS(239)*UV(52)+JVS(242)*UV(73)+JVS(243)*UV(75)+JVS(244)*UV(79)+JVS(245)*UV(81)
  JUV(53) = JVS(246)*UV(12)+JVS(248)*UV(53)+JVS(249)*UV(68)+JVS(250)*UV(70)+JVS(251)*UV(75)+JVS(252)*UV(76)+JVS(253)&
              &*UV(77)+JVS(254)*UV(78)+JVS(255)*UV(80)
  JUV(54) = JVS(257)*UV(7)+JVS(258)*UV(48)+JVS(259)*UV(54)+JVS(260)*UV(75)+JVS(261)*UV(76)+JVS(262)*UV(79)
  JUV(55) = JVS(263)*UV(4)+JVS(264)*UV(21)+JVS(265)*UV(33)+JVS(268)*UV(55)+JVS(269)*UV(75)+JVS(270)*UV(76)+JVS(271)&
              &*UV(77)+JVS(272)*UV(79)
  JUV(56) = JVS(273)*UV(56)+JVS(274)*UV(59)+JVS(275)*UV(76)+JVS(276)*UV(79)+JVS(277)*UV(80)
  JUV(57) = JVS(278)*UV(5)+JVS(279)*UV(22)+JVS(280)*UV(57)+JVS(281)*UV(75)+JVS(282)*UV(76)+JVS(283)*UV(77)+JVS(284)&
              &*UV(79)
  JUV(58) = JVS(285)*UV(27)+JVS(286)*UV(37)+JVS(288)*UV(45)+JVS(289)*UV(53)+JVS(290)*UV(58)+JVS(291)*UV(59)+JVS(292)&
              &*UV(60)+JVS(293)*UV(61)+JVS(294)*UV(62)+JVS(295)*UV(66)+JVS(296)*UV(67)+JVS(297)*UV(68)+JVS(298)*UV(69)&
              &+JVS(299)*UV(70)+JVS(300)*UV(72)+JVS(302)*UV(75)+JVS(303)*UV(76)+JVS(304)*UV(77)+JVS(305)*UV(78)+JVS(307)&
              &*UV(80)+JVS(308)*UV(81)
  JUV(59) = JVS(309)*UV(59)+JVS(310)*UV(75)+JVS(311)*UV(80)+JVS(312)*UV(81)
  JUV(60) = JVS(313)*UV(56)+JVS(315)*UV(60)+JVS(316)*UV(70)+JVS(317)*UV(71)+JVS(318)*UV(75)+JVS(319)*UV(76)+JVS(320)&
              &*UV(79)+JVS(321)*UV(80)
  JUV(61) = JVS(323)*UV(59)+JVS(324)*UV(61)+JVS(325)*UV(69)+JVS(326)*UV(75)+JVS(327)*UV(80)+JVS(328)*UV(81)
  JUV(62) = JVS(329)*UV(17)+JVS(330)*UV(21)+JVS(331)*UV(24)+JVS(332)*UV(25)+JVS(333)*UV(35)+JVS(334)*UV(48)+JVS(335)&
              &*UV(49)+JVS(336)*UV(50)+JVS(337)*UV(54)+JVS(338)*UV(55)+JVS(339)*UV(57)+JVS(340)*UV(61)+JVS(341)*UV(62)&
              &+JVS(342)*UV(69)+JVS(343)*UV(75)+JVS(344)*UV(76)+JVS(345)*UV(77)+JVS(347)*UV(80)+JVS(348)*UV(81)
  JUV(63) = JVS(349)*UV(23)+JVS(350)*UV(51)+JVS(353)*UV(63)+JVS(355)*UV(75)+JVS(356)*UV(76)+JVS(357)*UV(77)+JVS(358)&
              &*UV(79)
  JUV(64) = JVS(359)*UV(38)+JVS(360)*UV(46)+JVS(361)*UV(64)+JVS(362)*UV(75)+JVS(363)*UV(76)+JVS(364)*UV(79)+JVS(365)&
              &*UV(80)
  JUV(65) = JVS(367)*UV(35)+JVS(368)*UV(43)+JVS(371)*UV(63)+JVS(373)*UV(65)+JVS(374)*UV(68)+JVS(376)*UV(70)+JVS(379)&
              &*UV(75)+JVS(380)*UV(76)+JVS(381)*UV(77)+JVS(382)*UV(78)+JVS(384)*UV(80)
  JUV(66) = JVS(386)*UV(38)+JVS(387)*UV(42)+JVS(388)*UV(46)+JVS(389)*UV(56)+JVS(390)*UV(59)+JVS(391)*UV(64)+JVS(392)&
              &*UV(66)+JVS(393)*UV(71)+JVS(394)*UV(75)+JVS(395)*UV(76)+JVS(396)*UV(77)+JVS(397)*UV(78)+JVS(398)*UV(79)&
              &+JVS(399)*UV(80)+JVS(400)*UV(81)
  JUV(67) = JVS(401)*UV(29)+JVS(402)*UV(31)+JVS(403)*UV(44)+JVS(404)*UV(45)+JVS(407)*UV(63)+JVS(409)*UV(65)+JVS(410)&
              &*UV(66)+JVS(411)*UV(67)+JVS(412)*UV(68)+JVS(413)*UV(69)+JVS(414)*UV(70)+JVS(418)*UV(75)+JVS(419)*UV(76)&
              &+JVS(420)*UV(77)+JVS(421)*UV(78)+JVS(423)*UV(80)+JVS(424)*UV(81)
  JUV(68) = JVS(425)*UV(15)+JVS(426)*UV(28)+JVS(427)*UV(42)+JVS(429)*UV(68)+JVS(432)*UV(75)+JVS(433)*UV(76)+JVS(434)&
              &*UV(77)+JVS(435)*UV(78)+JVS(436)*UV(79)+JVS(437)*UV(80)
  JUV(69) = JVS(439)*UV(38)+JVS(440)*UV(42)+JVS(441)*UV(46)+JVS(442)*UV(56)+JVS(443)*UV(59)+JVS(444)*UV(64)+JVS(445)&
              &*UV(69)+JVS(446)*UV(71)+JVS(447)*UV(75)+JVS(448)*UV(76)+JVS(449)*UV(77)+JVS(450)*UV(78)+JVS(451)*UV(79)&
              &+JVS(452)*UV(80)+JVS(453)*UV(81)
  JUV(70) = JVS(454)*UV(20)+JVS(455)*UV(66)+JVS(456)*UV(69)+JVS(457)*UV(70)+JVS(459)*UV(75)+JVS(460)*UV(76)+JVS(461)&
              &*UV(77)+JVS(462)*UV(78)+JVS(463)*UV(79)+JVS(464)*UV(80)
  JUV(71) = JVS(466)*UV(42)+JVS(467)*UV(59)+JVS(468)*UV(71)+JVS(469)*UV(75)+JVS(470)*UV(76)+JVS(471)*UV(77)+JVS(472)&
              &*UV(78)+JVS(473)*UV(79)+JVS(474)*UV(80)
  JUV(72) = JVS(476)*UV(12)+JVS(477)*UV(18)+JVS(478)*UV(23)+JVS(479)*UV(24)+JVS(480)*UV(27)+JVS(482)*UV(32)+JVS(483)&
              &*UV(35)+JVS(484)*UV(36)+JVS(485)*UV(40)+JVS(486)*UV(41)+JVS(487)*UV(42)+JVS(488)*UV(43)+JVS(489)*UV(47)&
              &+JVS(490)*UV(48)+JVS(491)*UV(50)+JVS(492)*UV(53)+JVS(493)*UV(54)+JVS(494)*UV(55)+JVS(495)*UV(56)+JVS(496)&
              &*UV(57)+JVS(497)*UV(59)+JVS(498)*UV(60)+JVS(499)*UV(61)+JVS(500)*UV(63)+JVS(502)*UV(65)+JVS(503)*UV(66)&
              &+JVS(504)*UV(68)+JVS(505)*UV(69)+JVS(506)*UV(70)+JVS(507)*UV(71)+JVS(508)*UV(72)+JVS(510)*UV(74)+JVS(511)&
              &*UV(75)+JVS(512)*UV(76)+JVS(513)*UV(77)+JVS(514)*UV(78)+JVS(515)*UV(79)+JVS(516)*UV(80)+JVS(517)*UV(81)
  JUV(73) = JVS(518)*UV(13)+JVS(519)*UV(14)+JVS(520)*UV(24)+JVS(521)*UV(26)+JVS(522)*UV(29)+JVS(523)*UV(30)+JVS(524)&
              &*UV(36)+JVS(525)*UV(39)+JVS(526)*UV(43)+JVS(527)*UV(44)+JVS(528)*UV(46)+JVS(529)*UV(49)+JVS(530)*UV(50)&
              &+JVS(531)*UV(52)+JVS(532)*UV(54)+JVS(533)*UV(55)+JVS(534)*UV(56)+JVS(535)*UV(57)+JVS(537)*UV(60)+JVS(540)&
              &*UV(63)+JVS(541)*UV(64)+JVS(543)*UV(68)+JVS(545)*UV(70)+JVS(546)*UV(71)+JVS(548)*UV(73)+JVS(549)*UV(74)&
              &+JVS(550)*UV(75)+JVS(551)*UV(76)+JVS(552)*UV(77)+JVS(553)*UV(78)+JVS(554)*UV(79)+JVS(555)*UV(80)+JVS(556)&
              &*UV(81)
  JUV(74) = JVS(557)*UV(20)+JVS(558)*UV(43)+JVS(559)*UV(59)+JVS(560)*UV(66)+JVS(563)*UV(73)+JVS(564)*UV(74)+JVS(565)&
              &*UV(75)+JVS(566)*UV(76)+JVS(567)*UV(77)+JVS(568)*UV(78)+JVS(569)*UV(79)+JVS(570)*UV(80)+JVS(571)*UV(81)
  JUV(75) = JVS(572)*UV(2)+JVS(573)*UV(4)+JVS(574)*UV(5)+JVS(575)*UV(6)+JVS(576)*UV(7)+JVS(577)*UV(8)+JVS(578)*UV(9)&
              &+JVS(579)*UV(10)+JVS(580)*UV(11)+JVS(581)*UV(15)+JVS(582)*UV(16)+JVS(583)*UV(17)+JVS(584)*UV(18)+JVS(585)&
              &*UV(19)+JVS(586)*UV(20)+JVS(587)*UV(21)+JVS(588)*UV(22)+JVS(589)*UV(23)+JVS(590)*UV(25)+JVS(591)*UV(26)&
              &+JVS(592)*UV(27)+JVS(593)*UV(28)+JVS(594)*UV(29)+JVS(595)*UV(31)+JVS(596)*UV(32)+JVS(597)*UV(33)+JVS(598)&
              &*UV(34)+JVS(599)*UV(35)+JVS(600)*UV(36)+JVS(601)*UV(37)+JVS(602)*UV(38)+JVS(603)*UV(39)+JVS(604)*UV(40)&
              &+JVS(605)*UV(41)+JVS(606)*UV(42)+JVS(607)*UV(43)+JVS(610)*UV(46)+JVS(611)*UV(47)+JVS(612)*UV(48)+JVS(615)&
              &*UV(51)+JVS(616)*UV(52)+JVS(617)*UV(53)+JVS(620)*UV(56)+JVS(622)*UV(58)+JVS(623)*UV(59)+JVS(624)*UV(60)&
              &+JVS(625)*UV(61)+JVS(626)*UV(62)+JVS(629)*UV(65)+JVS(630)*UV(66)+JVS(631)*UV(67)+JVS(633)*UV(69)+JVS(636)&
              &*UV(72)+JVS(637)*UV(73)+JVS(639)*UV(75)+JVS(640)*UV(76)+JVS(643)*UV(79)+JVS(644)*UV(80)+JVS(645)*UV(81)
  JUV(76) = JVS(646)*UV(3)+JVS(647)*UV(24)+JVS(648)*UV(30)+JVS(649)*UV(44)+JVS(650)*UV(47)+JVS(651)*UV(49)+JVS(652)&
              &*UV(50)+JVS(653)*UV(52)+JVS(654)*UV(54)+JVS(655)*UV(55)+JVS(656)*UV(56)+JVS(657)*UV(57)+JVS(660)*UV(63)&
              &+JVS(661)*UV(64)+JVS(662)*UV(68)+JVS(664)*UV(70)+JVS(665)*UV(71)+JVS(667)*UV(73)+JVS(668)*UV(74)+JVS(670)&
              &*UV(76)+JVS(671)*UV(77)+JVS(672)*UV(78)+JVS(673)*UV(79)+JVS(674)*UV(80)+JVS(675)*UV(81)
  JUV(77) = JVS(676)*UV(18)+JVS(677)*UV(32)+JVS(678)*UV(34)+JVS(679)*UV(36)+JVS(680)*UV(40)+JVS(681)*UV(47)+JVS(682)&
              &*UV(51)+JVS(684)*UV(55)+JVS(685)*UV(57)+JVS(687)*UV(61)+JVS(688)*UV(62)+JVS(689)*UV(63)+JVS(691)*UV(68)&
              &+JVS(692)*UV(69)+JVS(693)*UV(70)+JVS(694)*UV(71)+JVS(697)*UV(74)+JVS(698)*UV(75)+JVS(699)*UV(76)+JVS(700)&
              &*UV(77)+JVS(701)*UV(78)+JVS(702)*UV(79)+JVS(704)*UV(81)
  JUV(78) = JVS(705)*UV(23)+JVS(706)*UV(25)+JVS(707)*UV(32)+JVS(708)*UV(33)+JVS(709)*UV(36)+JVS(710)*UV(45)+JVS(712)&
              &*UV(49)+JVS(713)*UV(51)+JVS(716)*UV(62)+JVS(717)*UV(63)+JVS(719)*UV(65)+JVS(720)*UV(66)+JVS(721)*UV(67)&
              &+JVS(722)*UV(68)+JVS(723)*UV(69)+JVS(724)*UV(70)+JVS(725)*UV(71)+JVS(726)*UV(73)+JVS(727)*UV(74)+JVS(728)&
              &*UV(75)+JVS(729)*UV(76)+JVS(730)*UV(77)+JVS(731)*UV(78)+JVS(732)*UV(79)+JVS(733)*UV(80)
  JUV(79) = JVS(735)*UV(9)+JVS(736)*UV(10)+JVS(737)*UV(12)+JVS(738)*UV(14)+JVS(739)*UV(16)+JVS(740)*UV(17)+JVS(741)&
              &*UV(18)+JVS(742)*UV(19)+JVS(743)*UV(20)+JVS(744)*UV(21)+JVS(745)*UV(22)+JVS(746)*UV(24)+JVS(747)*UV(26)&
              &+JVS(748)*UV(27)+JVS(750)*UV(35)+JVS(751)*UV(37)+JVS(752)*UV(38)+JVS(753)*UV(40)+JVS(754)*UV(41)+JVS(755)&
              &*UV(42)+JVS(756)*UV(43)+JVS(757)*UV(44)+JVS(758)*UV(45)+JVS(759)*UV(46)+JVS(760)*UV(47)+JVS(761)*UV(48)&
              &+JVS(762)*UV(49)+JVS(763)*UV(50)+JVS(764)*UV(52)+JVS(765)*UV(53)+JVS(766)*UV(54)+JVS(767)*UV(55)+JVS(768)&
              &*UV(56)+JVS(769)*UV(57)+JVS(770)*UV(58)+JVS(771)*UV(59)+JVS(772)*UV(60)+JVS(773)*UV(61)+JVS(774)*UV(62)&
              &+JVS(775)*UV(63)+JVS(776)*UV(64)+JVS(777)*UV(65)+JVS(778)*UV(66)+JVS(779)*UV(67)+JVS(780)*UV(68)+JVS(781)&
              &*UV(69)+JVS(782)*UV(70)+JVS(783)*UV(71)+JVS(784)*UV(72)+JVS(785)*UV(73)+JVS(786)*UV(74)+JVS(787)*UV(75)&
              &+JVS(788)*UV(76)+JVS(789)*UV(77)+JVS(790)*UV(78)+JVS(791)*UV(79)+JVS(792)*UV(80)+JVS(793)*UV(81)
  JUV(80) = JVS(794)*UV(13)+JVS(795)*UV(16)+JVS(796)*UV(26)+JVS(797)*UV(36)+JVS(798)*UV(39)+JVS(799)*UV(43)+JVS(800)&
              &*UV(46)+JVS(801)*UV(56)+JVS(802)*UV(59)+JVS(803)*UV(60)+JVS(804)*UV(61)+JVS(805)*UV(62)+JVS(806)*UV(67)&
              &+JVS(807)*UV(68)+JVS(809)*UV(70)+JVS(810)*UV(71)+JVS(811)*UV(72)+JVS(812)*UV(73)+JVS(813)*UV(74)+JVS(814)&
              &*UV(75)+JVS(815)*UV(76)+JVS(818)*UV(79)+JVS(819)*UV(80)+JVS(820)*UV(81)
  JUV(81) = JVS(821)*UV(27)+JVS(822)*UV(46)+JVS(823)*UV(52)+JVS(824)*UV(59)+JVS(825)*UV(61)+JVS(826)*UV(66)+JVS(827)&
              &*UV(69)+JVS(830)*UV(73)+JVS(831)*UV(74)+JVS(832)*UV(75)+JVS(833)*UV(76)+JVS(835)*UV(78)+JVS(836)*UV(79)&
              &+JVS(837)*UV(80)+JVS(838)*UV(81)
      
END SUBROUTINE Jac_SP_Vec

! End of Jac_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! JacTR_SP_Vec - sparse multiplication: sparse Jacobian transposed times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JTUV      - Jacobian transposed times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE JacTR_SP_Vec ( JVS, UV, JTUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JTUV - Jacobian transposed times user vector
  REAL(kind=dp) :: JTUV(NVAR)

  JTUV(1) = JVS(1)*UV(1)
  JTUV(2) = JVS(4)*UV(2)+JVS(572)*UV(75)
  JTUV(3) = JVS(6)*UV(3)+JVS(200)*UV(47)+JVS(646)*UV(76)
  JTUV(4) = JVS(8)*UV(4)+JVS(263)*UV(55)+JVS(573)*UV(75)
  JTUV(5) = JVS(10)*UV(5)+JVS(278)*UV(57)+JVS(574)*UV(75)
  JTUV(6) = JVS(12)*UV(6)+JVS(74)*UV(24)+JVS(575)*UV(75)
  JTUV(7) = JVS(14)*UV(7)+JVS(257)*UV(54)+JVS(576)*UV(75)
  JTUV(8) = JVS(2)*UV(1)+JVS(16)*UV(8)+JVS(577)*UV(75)
  JTUV(9) = JVS(20)*UV(9)+JVS(578)*UV(75)+JVS(735)*UV(79)
  JTUV(10) = JVS(23)*UV(10)+JVS(25)*UV(11)+JVS(182)*UV(44)+JVS(579)*UV(75)+JVS(736)*UV(79)
  JTUV(11) = JVS(26)*UV(11)+JVS(34)*UV(14)+JVS(580)*UV(75)
  JTUV(12) = JVS(28)*UV(12)+JVS(246)*UV(53)+JVS(476)*UV(72)+JVS(737)*UV(79)
  JTUV(13) = JVS(31)*UV(13)+JVS(152)*UV(39)+JVS(518)*UV(73)+JVS(794)*UV(80)
  JTUV(14) = JVS(35)*UV(14)+JVS(188)*UV(45)+JVS(519)*UV(73)+JVS(738)*UV(79)
  JTUV(15) = JVS(38)*UV(15)+JVS(425)*UV(68)+JVS(581)*UV(75)
  JTUV(16) = JVS(17)*UV(8)+JVS(42)*UV(16)+JVS(153)*UV(39)+JVS(582)*UV(75)+JVS(739)*UV(79)+JVS(795)*UV(80)
  JTUV(17) = JVS(45)*UV(17)+JVS(329)*UV(62)+JVS(583)*UV(75)+JVS(740)*UV(79)
  JTUV(18) = JVS(49)*UV(18)+JVS(477)*UV(72)+JVS(584)*UV(75)+JVS(676)*UV(77)+JVS(741)*UV(79)
  JTUV(19) = JVS(53)*UV(19)+JVS(201)*UV(47)+JVS(585)*UV(75)+JVS(742)*UV(79)
  JTUV(20) = JVS(58)*UV(20)+JVS(454)*UV(70)+JVS(557)*UV(74)+JVS(586)*UV(75)+JVS(743)*UV(79)
  JTUV(21) = JVS(62)*UV(21)+JVS(264)*UV(55)+JVS(330)*UV(62)+JVS(587)*UV(75)+JVS(744)*UV(79)
  JTUV(22) = JVS(66)*UV(22)+JVS(226)*UV(51)+JVS(279)*UV(57)+JVS(588)*UV(75)+JVS(745)*UV(79)
  JTUV(23) = JVS(70)*UV(23)+JVS(349)*UV(63)+JVS(478)*UV(72)+JVS(589)*UV(75)+JVS(705)*UV(78)
  JTUV(24) = JVS(75)*UV(24)+JVS(227)*UV(51)+JVS(331)*UV(62)+JVS(479)*UV(72)+JVS(520)*UV(73)+JVS(647)*UV(76)+JVS(746)&
               &*UV(79)
  JTUV(25) = JVS(78)*UV(25)+JVS(212)*UV(49)+JVS(332)*UV(62)+JVS(590)*UV(75)+JVS(706)*UV(78)
  JTUV(26) = JVS(82)*UV(26)+JVS(521)*UV(73)+JVS(591)*UV(75)+JVS(747)*UV(79)+JVS(796)*UV(80)
  JTUV(27) = JVS(86)*UV(27)+JVS(103)*UV(30)+JVS(122)*UV(34)+JVS(285)*UV(58)+JVS(480)*UV(72)+JVS(592)*UV(75)+JVS(748)&
               &*UV(79)+JVS(821)*UV(81)
  JTUV(28) = JVS(89)*UV(28)+JVS(426)*UV(68)+JVS(593)*UV(75)
  JTUV(29) = JVS(97)*UV(29)+JVS(401)*UV(67)+JVS(522)*UV(73)+JVS(594)*UV(75)
  JTUV(30) = JVS(29)*UV(12)+JVS(104)*UV(30)+JVS(523)*UV(73)+JVS(648)*UV(76)
  JTUV(31) = JVS(108)*UV(31)+JVS(140)*UV(37)+JVS(183)*UV(44)+JVS(189)*UV(45)+JVS(402)*UV(67)+JVS(595)*UV(75)
  JTUV(32) = JVS(112)*UV(32)+JVS(482)*UV(72)+JVS(596)*UV(75)+JVS(677)*UV(77)+JVS(707)*UV(78)
  JTUV(33) = JVS(117)*UV(33)+JVS(213)*UV(49)+JVS(265)*UV(55)+JVS(597)*UV(75)+JVS(708)*UV(78)
  JTUV(34) = JVS(123)*UV(34)+JVS(598)*UV(75)+JVS(678)*UV(77)
  JTUV(35) = JVS(132)*UV(35)+JVS(220)*UV(50)+JVS(333)*UV(62)+JVS(367)*UV(65)+JVS(483)*UV(72)+JVS(599)*UV(75)+JVS(750)&
               &*UV(79)
  JTUV(36) = JVS(136)*UV(36)+JVS(484)*UV(72)+JVS(524)*UV(73)+JVS(600)*UV(75)+JVS(679)*UV(77)+JVS(709)*UV(78)+JVS(797)&
               &*UV(80)
  JTUV(37) = JVS(141)*UV(37)+JVS(286)*UV(58)+JVS(601)*UV(75)+JVS(751)*UV(79)
  JTUV(38) = JVS(148)*UV(38)+JVS(228)*UV(51)+JVS(359)*UV(64)+JVS(386)*UV(66)+JVS(439)*UV(69)+JVS(602)*UV(75)+JVS(752)&
               &*UV(79)
  JTUV(39) = JVS(154)*UV(39)+JVS(525)*UV(73)+JVS(603)*UV(75)+JVS(798)*UV(80)
  JTUV(40) = JVS(54)*UV(19)+JVS(161)*UV(40)+JVS(202)*UV(47)+JVS(485)*UV(72)+JVS(604)*UV(75)+JVS(680)*UV(77)+JVS(753)&
               &*UV(79)
  JTUV(41) = JVS(166)*UV(41)+JVS(486)*UV(72)+JVS(605)*UV(75)+JVS(754)*UV(79)
  JTUV(42) = JVS(174)*UV(42)+JVS(387)*UV(66)+JVS(427)*UV(68)+JVS(440)*UV(69)+JVS(466)*UV(71)+JVS(487)*UV(72)+JVS(606)&
               &*UV(75)+JVS(755)*UV(79)
  JTUV(43) = JVS(178)*UV(43)+JVS(368)*UV(65)+JVS(488)*UV(72)+JVS(526)*UV(73)+JVS(558)*UV(74)+JVS(607)*UV(75)+JVS(756)&
               &*UV(79)+JVS(799)*UV(80)
  JTUV(44) = JVS(109)*UV(31)+JVS(142)*UV(37)+JVS(184)*UV(44)+JVS(190)*UV(45)+JVS(403)*UV(67)+JVS(527)*UV(73)+JVS(649)&
               &*UV(76)+JVS(757)*UV(79)
  JTUV(45) = JVS(143)*UV(37)+JVS(191)*UV(45)+JVS(288)*UV(58)+JVS(404)*UV(67)+JVS(710)*UV(78)+JVS(758)*UV(79)
  JTUV(46) = JVS(196)*UV(46)+JVS(360)*UV(64)+JVS(388)*UV(66)+JVS(441)*UV(69)+JVS(528)*UV(73)+JVS(610)*UV(75)+JVS(759)&
               &*UV(79)+JVS(800)*UV(80)+JVS(822)*UV(81)
  JTUV(47) = JVS(7)*UV(3)+JVS(55)*UV(19)+JVS(162)*UV(40)+JVS(203)*UV(47)+JVS(238)*UV(52)+JVS(489)*UV(72)+JVS(611)*UV(75)&
               &+JVS(650)*UV(76)+JVS(681)*UV(77)+JVS(760)*UV(79)
  JTUV(48) = JVS(118)*UV(33)+JVS(208)*UV(48)+JVS(229)*UV(51)+JVS(258)*UV(54)+JVS(334)*UV(62)+JVS(490)*UV(72)+JVS(612)&
               &*UV(75)+JVS(761)*UV(79)
  JTUV(49) = JVS(79)*UV(25)+JVS(215)*UV(49)+JVS(335)*UV(62)+JVS(529)*UV(73)+JVS(651)*UV(76)+JVS(712)*UV(78)+JVS(762)&
               &*UV(79)
  JTUV(50) = JVS(133)*UV(35)+JVS(221)*UV(50)+JVS(336)*UV(62)+JVS(491)*UV(72)+JVS(530)*UV(73)+JVS(652)*UV(76)+JVS(763)&
               &*UV(79)
  JTUV(51) = JVS(230)*UV(51)+JVS(350)*UV(63)+JVS(615)*UV(75)+JVS(682)*UV(77)+JVS(713)*UV(78)
  JTUV(52) = JVS(239)*UV(52)+JVS(531)*UV(73)+JVS(616)*UV(75)+JVS(653)*UV(76)+JVS(764)*UV(79)+JVS(823)*UV(81)
  JTUV(53) = JVS(144)*UV(37)+JVS(248)*UV(53)+JVS(289)*UV(58)+JVS(492)*UV(72)+JVS(617)*UV(75)+JVS(765)*UV(79)
  JTUV(54) = JVS(98)*UV(29)+JVS(119)*UV(33)+JVS(209)*UV(48)+JVS(231)*UV(51)+JVS(259)*UV(54)+JVS(337)*UV(62)+JVS(493)&
               &*UV(72)+JVS(532)*UV(73)+JVS(654)*UV(76)+JVS(766)*UV(79)
  JTUV(55) = JVS(46)*UV(17)+JVS(63)*UV(21)+JVS(167)*UV(41)+JVS(268)*UV(55)+JVS(338)*UV(62)+JVS(494)*UV(72)+JVS(533)&
               &*UV(73)+JVS(655)*UV(76)+JVS(684)*UV(77)+JVS(767)*UV(79)
  JTUV(56) = JVS(273)*UV(56)+JVS(313)*UV(60)+JVS(389)*UV(66)+JVS(442)*UV(69)+JVS(495)*UV(72)+JVS(534)*UV(73)+JVS(620)&
               &*UV(75)+JVS(656)*UV(76)+JVS(768)*UV(79)+JVS(801)*UV(80)
  JTUV(57) = JVS(67)*UV(22)+JVS(232)*UV(51)+JVS(280)*UV(57)+JVS(339)*UV(62)+JVS(496)*UV(72)+JVS(535)*UV(73)+JVS(657)&
               &*UV(76)+JVS(685)*UV(77)+JVS(769)*UV(79)
  JTUV(58) = JVS(290)*UV(58)+JVS(622)*UV(75)+JVS(770)*UV(79)
  JTUV(59) = JVS(124)*UV(34)+JVS(274)*UV(56)+JVS(291)*UV(58)+JVS(309)*UV(59)+JVS(323)*UV(61)+JVS(390)*UV(66)+JVS(443)&
               &*UV(69)+JVS(467)*UV(71)+JVS(497)*UV(72)+JVS(559)*UV(74)+JVS(623)*UV(75)+JVS(771)*UV(79)+JVS(802)*UV(80)&
               &+JVS(824)*UV(81)
  JTUV(60) = JVS(90)*UV(28)+JVS(292)*UV(58)+JVS(315)*UV(60)+JVS(498)*UV(72)+JVS(537)*UV(73)+JVS(624)*UV(75)+JVS(772)&
               &*UV(79)+JVS(803)*UV(80)
  JTUV(61) = JVS(99)*UV(29)+JVS(125)*UV(34)+JVS(163)*UV(40)+JVS(222)*UV(50)+JVS(293)*UV(58)+JVS(324)*UV(61)+JVS(340)&
               &*UV(62)+JVS(499)*UV(72)+JVS(625)*UV(75)+JVS(687)*UV(77)+JVS(773)*UV(79)+JVS(804)*UV(80)+JVS(825)*UV(81)
  JTUV(62) = JVS(155)*UV(39)+JVS(294)*UV(58)+JVS(341)*UV(62)+JVS(626)*UV(75)+JVS(688)*UV(77)+JVS(716)*UV(78)+JVS(774)&
               &*UV(79)+JVS(805)*UV(80)
  JTUV(63) = JVS(71)*UV(23)+JVS(168)*UV(41)+JVS(353)*UV(63)+JVS(371)*UV(65)+JVS(407)*UV(67)+JVS(500)*UV(72)+JVS(540)&
               &*UV(73)+JVS(660)*UV(76)+JVS(689)*UV(77)+JVS(717)*UV(78)+JVS(775)*UV(79)
  JTUV(64) = JVS(149)*UV(38)+JVS(233)*UV(51)+JVS(361)*UV(64)+JVS(391)*UV(66)+JVS(444)*UV(69)+JVS(541)*UV(73)+JVS(661)&
               &*UV(76)+JVS(776)*UV(79)
  JTUV(65) = JVS(373)*UV(65)+JVS(409)*UV(67)+JVS(502)*UV(72)+JVS(629)*UV(75)+JVS(719)*UV(78)+JVS(777)*UV(79)
  JTUV(66) = JVS(295)*UV(58)+JVS(392)*UV(66)+JVS(410)*UV(67)+JVS(455)*UV(70)+JVS(503)*UV(72)+JVS(560)*UV(74)+JVS(630)&
               &*UV(75)+JVS(720)*UV(78)+JVS(778)*UV(79)+JVS(826)*UV(81)
  JTUV(67) = JVS(156)*UV(39)+JVS(296)*UV(58)+JVS(411)*UV(67)+JVS(631)*UV(75)+JVS(721)*UV(78)+JVS(779)*UV(79)+JVS(806)&
               &*UV(80)
  JTUV(68) = JVS(39)*UV(15)+JVS(169)*UV(41)+JVS(249)*UV(53)+JVS(297)*UV(58)+JVS(374)*UV(65)+JVS(412)*UV(67)+JVS(429)&
               &*UV(68)+JVS(504)*UV(72)+JVS(543)*UV(73)+JVS(662)*UV(76)+JVS(691)*UV(77)+JVS(722)*UV(78)+JVS(780)*UV(79)&
               &+JVS(807)*UV(80)
  JTUV(69) = JVS(298)*UV(58)+JVS(325)*UV(61)+JVS(342)*UV(62)+JVS(413)*UV(67)+JVS(445)*UV(69)+JVS(456)*UV(70)+JVS(505)&
               &*UV(72)+JVS(633)*UV(75)+JVS(692)*UV(77)+JVS(723)*UV(78)+JVS(781)*UV(79)+JVS(827)*UV(81)
  JTUV(70) = JVS(59)*UV(20)+JVS(170)*UV(41)+JVS(250)*UV(53)+JVS(299)*UV(58)+JVS(316)*UV(60)+JVS(376)*UV(65)+JVS(414)&
               &*UV(67)+JVS(457)*UV(70)+JVS(506)*UV(72)+JVS(545)*UV(73)+JVS(664)*UV(76)+JVS(693)*UV(77)+JVS(724)*UV(78)&
               &+JVS(782)*UV(79)+JVS(809)*UV(80)
  JTUV(71) = JVS(91)*UV(28)+JVS(171)*UV(41)+JVS(175)*UV(42)+JVS(317)*UV(60)+JVS(393)*UV(66)+JVS(446)*UV(69)+JVS(468)&
               &*UV(71)+JVS(507)*UV(72)+JVS(546)*UV(73)+JVS(665)*UV(76)+JVS(694)*UV(77)+JVS(725)*UV(78)+JVS(783)*UV(79)&
               &+JVS(810)*UV(80)
  JTUV(72) = JVS(56)*UV(19)+JVS(157)*UV(39)+JVS(300)*UV(58)+JVS(508)*UV(72)+JVS(636)*UV(75)+JVS(784)*UV(79)+JVS(811)&
               &*UV(80)
  JTUV(73) = JVS(32)*UV(13)+JVS(36)*UV(14)+JVS(83)*UV(26)+JVS(137)*UV(36)+JVS(158)*UV(39)+JVS(179)*UV(43)+JVS(192)&
               &*UV(45)+JVS(242)*UV(52)+JVS(548)*UV(73)+JVS(563)*UV(74)+JVS(637)*UV(75)+JVS(667)*UV(76)+JVS(726)*UV(78)&
               &+JVS(785)*UV(79)+JVS(812)*UV(80)+JVS(830)*UV(81)
  JTUV(74) = JVS(113)*UV(32)+JVS(126)*UV(34)+JVS(180)*UV(43)+JVS(510)*UV(72)+JVS(549)*UV(73)+JVS(564)*UV(74)+JVS(668)&
               &*UV(76)+JVS(697)*UV(77)+JVS(727)*UV(78)+JVS(786)*UV(79)+JVS(813)*UV(80)+JVS(831)*UV(81)
  JTUV(75) = JVS(3)*UV(1)+JVS(5)*UV(2)+JVS(9)*UV(4)+JVS(11)*UV(5)+JVS(13)*UV(6)+JVS(15)*UV(7)+JVS(18)*UV(8)+JVS(21)&
               &*UV(9)+JVS(24)*UV(10)+JVS(27)*UV(11)+JVS(37)*UV(14)+JVS(40)*UV(15)+JVS(43)*UV(16)+JVS(47)*UV(17)+JVS(50)&
               &*UV(18)+JVS(57)*UV(19)+JVS(60)*UV(20)+JVS(64)*UV(21)+JVS(68)*UV(22)+JVS(72)*UV(23)+JVS(76)*UV(24)+JVS(80)&
               &*UV(25)+JVS(84)*UV(26)+JVS(87)*UV(27)+JVS(92)*UV(28)+JVS(100)*UV(29)+JVS(105)*UV(30)+JVS(110)*UV(31)&
               &+JVS(114)*UV(32)+JVS(120)*UV(33)+JVS(127)*UV(34)+JVS(134)*UV(35)+JVS(138)*UV(36)+JVS(145)*UV(37)+JVS(150)&
               &*UV(38)+JVS(159)*UV(39)+JVS(164)*UV(40)+JVS(172)*UV(41)+JVS(176)*UV(42)+JVS(181)*UV(43)+JVS(185)*UV(44)&
               &+JVS(197)*UV(46)+JVS(210)*UV(48)+JVS(217)*UV(49)+JVS(223)*UV(50)+JVS(234)*UV(51)+JVS(243)*UV(52)+JVS(251)&
               &*UV(53)+JVS(260)*UV(54)+JVS(269)*UV(55)+JVS(281)*UV(57)+JVS(302)*UV(58)+JVS(310)*UV(59)+JVS(318)*UV(60)&
               &+JVS(326)*UV(61)+JVS(343)*UV(62)+JVS(355)*UV(63)+JVS(362)*UV(64)+JVS(379)*UV(65)+JVS(394)*UV(66)+JVS(418)&
               &*UV(67)+JVS(432)*UV(68)+JVS(447)*UV(69)+JVS(459)*UV(70)+JVS(469)*UV(71)+JVS(511)*UV(72)+JVS(550)*UV(73)&
               &+JVS(565)*UV(74)+JVS(639)*UV(75)+JVS(698)*UV(77)+JVS(728)*UV(78)+JVS(787)*UV(79)+JVS(814)*UV(80)+JVS(832)&
               &*UV(81)
  JTUV(76) = JVS(30)*UV(12)+JVS(77)*UV(24)+JVS(93)*UV(28)+JVS(101)*UV(29)+JVS(106)*UV(30)+JVS(121)*UV(33)+JVS(146)&
               &*UV(37)+JVS(186)*UV(44)+JVS(194)*UV(45)+JVS(218)*UV(49)+JVS(224)*UV(50)+JVS(235)*UV(51)+JVS(252)*UV(53)&
               &+JVS(261)*UV(54)+JVS(270)*UV(55)+JVS(275)*UV(56)+JVS(282)*UV(57)+JVS(303)*UV(58)+JVS(319)*UV(60)+JVS(344)&
               &*UV(62)+JVS(356)*UV(63)+JVS(363)*UV(64)+JVS(380)*UV(65)+JVS(395)*UV(66)+JVS(419)*UV(67)+JVS(433)*UV(68)&
               &+JVS(448)*UV(69)+JVS(460)*UV(70)+JVS(470)*UV(71)+JVS(512)*UV(72)+JVS(551)*UV(73)+JVS(566)*UV(74)+JVS(640)&
               &*UV(75)+JVS(670)*UV(76)+JVS(699)*UV(77)+JVS(729)*UV(78)+JVS(788)*UV(79)+JVS(815)*UV(80)+JVS(833)*UV(81)
  JTUV(77) = JVS(48)*UV(17)+JVS(51)*UV(18)+JVS(94)*UV(28)+JVS(128)*UV(34)+JVS(173)*UV(41)+JVS(236)*UV(51)+JVS(253)&
               &*UV(53)+JVS(271)*UV(55)+JVS(283)*UV(57)+JVS(304)*UV(58)+JVS(345)*UV(62)+JVS(357)*UV(63)+JVS(381)*UV(65)&
               &+JVS(396)*UV(66)+JVS(420)*UV(67)+JVS(434)*UV(68)+JVS(449)*UV(69)+JVS(461)*UV(70)+JVS(471)*UV(71)+JVS(513)&
               &*UV(72)+JVS(552)*UV(73)+JVS(567)*UV(74)+JVS(671)*UV(76)+JVS(700)*UV(77)+JVS(730)*UV(78)+JVS(789)*UV(79)
  JTUV(78) = JVS(95)*UV(28)+JVS(115)*UV(32)+JVS(129)*UV(34)+JVS(139)*UV(36)+JVS(254)*UV(53)+JVS(305)*UV(58)+JVS(382)&
               &*UV(65)+JVS(397)*UV(66)+JVS(421)*UV(67)+JVS(435)*UV(68)+JVS(450)*UV(69)+JVS(462)*UV(70)+JVS(472)*UV(71)&
               &+JVS(514)*UV(72)+JVS(553)*UV(73)+JVS(568)*UV(74)+JVS(672)*UV(76)+JVS(701)*UV(77)+JVS(731)*UV(78)+JVS(790)&
               &*UV(79)+JVS(835)*UV(81)
  JTUV(79) = JVS(22)*UV(9)+JVS(41)*UV(15)+JVS(52)*UV(18)+JVS(61)*UV(20)+JVS(65)*UV(21)+JVS(69)*UV(22)+JVS(73)*UV(23)&
               &+JVS(81)*UV(25)+JVS(85)*UV(26)+JVS(111)*UV(31)+JVS(116)*UV(32)+JVS(130)*UV(34)+JVS(135)*UV(35)+JVS(151)&
               &*UV(38)+JVS(177)*UV(42)+JVS(187)*UV(44)+JVS(211)*UV(48)+JVS(219)*UV(49)+JVS(225)*UV(50)+JVS(244)*UV(52)&
               &+JVS(262)*UV(54)+JVS(272)*UV(55)+JVS(276)*UV(56)+JVS(284)*UV(57)+JVS(320)*UV(60)+JVS(358)*UV(63)+JVS(364)&
               &*UV(64)+JVS(398)*UV(66)+JVS(436)*UV(68)+JVS(451)*UV(69)+JVS(463)*UV(70)+JVS(473)*UV(71)+JVS(515)*UV(72)&
               &+JVS(554)*UV(73)+JVS(569)*UV(74)+JVS(643)*UV(75)+JVS(673)*UV(76)+JVS(702)*UV(77)+JVS(732)*UV(78)+JVS(791)&
               &*UV(79)+JVS(818)*UV(80)+JVS(836)*UV(81)
  JTUV(80) = JVS(19)*UV(8)+JVS(33)*UV(13)+JVS(44)*UV(16)+JVS(96)*UV(28)+JVS(102)*UV(29)+JVS(160)*UV(39)+JVS(198)*UV(46)&
               &+JVS(255)*UV(53)+JVS(277)*UV(56)+JVS(307)*UV(58)+JVS(311)*UV(59)+JVS(321)*UV(60)+JVS(327)*UV(61)+JVS(347)&
               &*UV(62)+JVS(365)*UV(64)+JVS(384)*UV(65)+JVS(399)*UV(66)+JVS(423)*UV(67)+JVS(437)*UV(68)+JVS(452)*UV(69)&
               &+JVS(464)*UV(70)+JVS(474)*UV(71)+JVS(516)*UV(72)+JVS(555)*UV(73)+JVS(570)*UV(74)+JVS(644)*UV(75)+JVS(674)&
               &*UV(76)+JVS(733)*UV(78)+JVS(792)*UV(79)+JVS(819)*UV(80)+JVS(837)*UV(81)
  JTUV(81) = JVS(88)*UV(27)+JVS(131)*UV(34)+JVS(165)*UV(40)+JVS(199)*UV(46)+JVS(207)*UV(47)+JVS(245)*UV(52)+JVS(308)&
               &*UV(58)+JVS(312)*UV(59)+JVS(328)*UV(61)+JVS(348)*UV(62)+JVS(400)*UV(66)+JVS(424)*UV(67)+JVS(453)*UV(69)&
               &+JVS(517)*UV(72)+JVS(556)*UV(73)+JVS(571)*UV(74)+JVS(645)*UV(75)+JVS(675)*UV(76)+JVS(704)*UV(77)+JVS(793)&
               &*UV(79)+JVS(820)*UV(80)+JVS(838)*UV(81)
      
END SUBROUTINE JacTR_SP_Vec

! End of JacTR_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE mozart_Jacobian

