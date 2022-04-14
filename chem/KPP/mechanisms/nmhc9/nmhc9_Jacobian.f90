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
! File                 : nmhc9_Jacobian.f90
! Time                 : Tue Apr 12 23:43:14 2022
! Working directory    : /network/rit/lab/lulab/chinan/WRF/WRFV4.0/WRF/chem/KPP/mechanisms/nmhc9
! Equation file        : nmhc9.kpp
! Output root filename : nmhc9
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE nmhc9_Jacobian

  USE nmhc9_Parameters
  USE nmhc9_JacobianSP

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

  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(55)
  JUV(2) = JVS(3)*UV(2)+JVS(4)*UV(55)
  JUV(3) = JVS(5)*UV(3)+JVS(6)*UV(20)+JVS(7)*UV(55)
  JUV(4) = JVS(8)*UV(4)+JVS(9)*UV(55)
  JUV(5) = JVS(10)*UV(5)+JVS(11)*UV(52)+JVS(12)*UV(54)
  JUV(6) = JVS(13)*UV(6)+JVS(14)*UV(50)+JVS(15)*UV(54)
  JUV(7) = JVS(16)*UV(7)+JVS(17)*UV(33)+JVS(18)*UV(55)
  JUV(8) = JVS(19)*UV(8)+JVS(20)*UV(38)+JVS(21)*UV(53)+JVS(22)*UV(55)+JVS(23)*UV(57)
  JUV(9) = JVS(24)*UV(9)+JVS(25)*UV(51)+JVS(26)*UV(55)+JVS(27)*UV(57)
  JUV(10) = JVS(28)*UV(10)+JVS(29)*UV(47)+JVS(30)*UV(50)+JVS(31)*UV(51)+JVS(32)*UV(55)+JVS(33)*UV(57)
  JUV(11) = JVS(34)*UV(11)+JVS(35)*UV(48)+JVS(36)*UV(55)+JVS(37)*UV(57)
  JUV(12) = JVS(38)*UV(12)+JVS(39)*UV(50)+JVS(40)*UV(55)+JVS(41)*UV(57)
  JUV(13) = JVS(42)*UV(13)+JVS(43)*UV(51)+JVS(44)*UV(54)+JVS(45)*UV(55)
  JUV(14) = JVS(46)*UV(14)+JVS(47)*UV(49)+JVS(48)*UV(54)+JVS(49)*UV(55)
  JUV(15) = JVS(50)*UV(15)+JVS(51)*UV(19)+JVS(52)*UV(37)+JVS(53)*UV(38)+JVS(54)*UV(41)+JVS(55)*UV(53)+JVS(56)*UV(55)
  JUV(16) = JVS(57)*UV(16)+JVS(58)*UV(47)+JVS(59)*UV(55)+JVS(60)*UV(57)
  JUV(17) = JVS(61)*UV(17)+JVS(62)*UV(34)+JVS(63)*UV(55)+JVS(64)*UV(57)
  JUV(18) = JVS(65)*UV(18)+JVS(66)*UV(34)+JVS(67)*UV(55)+JVS(68)*UV(56)
  JUV(19) = JVS(69)*UV(19)+JVS(70)*UV(53)+JVS(71)*UV(55)
  JUV(20) = JVS(72)*UV(20)+JVS(73)*UV(42)+JVS(74)*UV(55)+JVS(75)*UV(57)
  JUV(21) = JVS(76)*UV(21)+JVS(77)*UV(24)+JVS(78)*UV(53)
  JUV(22) = JVS(79)*UV(22)+JVS(80)*UV(54)+JVS(81)*UV(55)+JVS(82)*UV(57)
  JUV(23) = JVS(83)*UV(23)+JVS(84)*UV(36)+JVS(85)*UV(55)+JVS(86)*UV(57)
  JUV(24) = JVS(87)*UV(21)+JVS(88)*UV(24)+JVS(89)*UV(37)+JVS(90)*UV(53)+JVS(91)*UV(55)
  JUV(25) = JVS(92)*UV(5)+JVS(93)*UV(25)+JVS(94)*UV(45)+JVS(95)*UV(46)+JVS(96)*UV(52)+JVS(97)*UV(54)+JVS(98)*UV(55)
  JUV(26) = JVS(99)*UV(26)+JVS(100)*UV(40)+JVS(101)*UV(55)+JVS(102)*UV(57)
  JUV(27) = JVS(103)*UV(27)+JVS(104)*UV(43)+JVS(105)*UV(55)+JVS(106)*UV(57)
  JUV(28) = JVS(107)*UV(17)+JVS(108)*UV(18)+JVS(109)*UV(28)+JVS(110)*UV(34)+JVS(111)*UV(50)+JVS(112)*UV(55)+JVS(113)&
              &*UV(56)
  JUV(29) = JVS(115)*UV(29)+JVS(116)*UV(49)+JVS(117)*UV(55)+JVS(118)*UV(57)
  JUV(30) = JVS(119)*UV(27)+JVS(120)*UV(30)+JVS(121)*UV(43)+JVS(122)*UV(44)+JVS(123)*UV(50)+JVS(124)*UV(55)+JVS(125)&
              &*UV(56)
  JUV(31) = JVS(127)*UV(7)+JVS(128)*UV(19)+JVS(129)*UV(29)+JVS(130)*UV(31)+JVS(132)*UV(35)+JVS(133)*UV(37)+JVS(134)&
              &*UV(38)+JVS(135)*UV(41)+JVS(136)*UV(45)+JVS(137)*UV(46)+JVS(138)*UV(49)+JVS(139)*UV(50)+JVS(140)*UV(52)&
              &+JVS(141)*UV(53)+JVS(142)*UV(55)+JVS(143)*UV(56)
  JUV(32) = JVS(145)*UV(32)+JVS(146)*UV(34)+JVS(147)*UV(36)+JVS(148)*UV(37)+JVS(149)*UV(43)+JVS(150)*UV(47)+JVS(151)&
              &*UV(48)+JVS(152)*UV(49)+JVS(153)*UV(50)+JVS(154)*UV(53)+JVS(155)*UV(55)
  JUV(33) = JVS(156)*UV(33)+JVS(157)*UV(38)+JVS(158)*UV(48)+JVS(159)*UV(52)+JVS(160)*UV(55)+JVS(161)*UV(56)
  JUV(34) = JVS(162)*UV(4)+JVS(163)*UV(17)+JVS(164)*UV(34)+JVS(165)*UV(50)+JVS(166)*UV(55)+JVS(167)*UV(56)+JVS(168)&
              &*UV(57)
  JUV(35) = JVS(169)*UV(23)+JVS(170)*UV(29)+JVS(171)*UV(35)+JVS(172)*UV(36)+JVS(173)*UV(37)+JVS(174)*UV(39)+JVS(175)&
              &*UV(41)+JVS(176)*UV(48)+JVS(177)*UV(49)+JVS(178)*UV(50)+JVS(179)*UV(53)+JVS(180)*UV(55)+JVS(181)*UV(56)
  JUV(36) = JVS(183)*UV(23)+JVS(184)*UV(28)+JVS(186)*UV(36)+JVS(187)*UV(50)+JVS(188)*UV(55)+JVS(189)*UV(56)+JVS(190)&
              &*UV(57)
  JUV(37) = JVS(191)*UV(37)+JVS(192)*UV(52)+JVS(193)*UV(53)+JVS(194)*UV(55)
  JUV(38) = JVS(195)*UV(38)+JVS(196)*UV(52)+JVS(197)*UV(53)+JVS(198)*UV(55)
  JUV(39) = JVS(199)*UV(14)+JVS(200)*UV(26)+JVS(201)*UV(29)+JVS(202)*UV(33)+JVS(203)*UV(36)+JVS(205)*UV(39)+JVS(207)&
              &*UV(48)+JVS(208)*UV(49)+JVS(209)*UV(50)+JVS(213)*UV(55)+JVS(214)*UV(56)
  JUV(40) = JVS(216)*UV(19)+JVS(217)*UV(26)+JVS(218)*UV(37)+JVS(219)*UV(40)+JVS(222)*UV(55)+JVS(223)*UV(56)+JVS(224)&
              &*UV(57)
  JUV(41) = JVS(225)*UV(11)+JVS(226)*UV(33)+JVS(227)*UV(38)+JVS(228)*UV(41)+JVS(229)*UV(48)+JVS(230)*UV(50)+JVS(232)&
              &*UV(53)+JVS(233)*UV(55)+JVS(234)*UV(56)
  JUV(42) = JVS(236)*UV(20)+JVS(237)*UV(30)+JVS(238)*UV(42)+JVS(242)*UV(55)+JVS(243)*UV(56)+JVS(244)*UV(57)
  JUV(43) = JVS(245)*UV(2)+JVS(246)*UV(27)+JVS(247)*UV(43)+JVS(248)*UV(50)+JVS(249)*UV(55)+JVS(250)*UV(56)+JVS(251)&
              &*UV(57)
  JUV(44) = JVS(252)*UV(37)+JVS(253)*UV(40)+JVS(254)*UV(42)+JVS(255)*UV(43)+JVS(256)*UV(44)+JVS(258)*UV(52)+JVS(260)&
              &*UV(55)+JVS(261)*UV(56)
  JUV(45) = JVS(263)*UV(6)+JVS(264)*UV(7)+JVS(265)*UV(11)+JVS(266)*UV(12)+JVS(267)*UV(13)+JVS(268)*UV(19)+JVS(269)&
              &*UV(21)+JVS(270)*UV(24)+JVS(271)*UV(26)+JVS(272)*UV(29)+JVS(273)*UV(32)+JVS(274)*UV(33)+JVS(275)*UV(34)&
              &+JVS(276)*UV(36)+JVS(277)*UV(37)+JVS(278)*UV(38)+JVS(279)*UV(39)+JVS(280)*UV(40)+JVS(281)*UV(41)+JVS(282)&
              &*UV(43)+JVS(283)*UV(45)+JVS(284)*UV(47)+JVS(285)*UV(48)+JVS(286)*UV(49)+JVS(287)*UV(50)+JVS(288)*UV(51)&
              &+JVS(289)*UV(52)+JVS(290)*UV(53)+JVS(292)*UV(55)+JVS(293)*UV(56)+JVS(294)*UV(57)
  JUV(46) = JVS(295)*UV(16)+JVS(296)*UV(20)+JVS(297)*UV(26)+JVS(298)*UV(27)+JVS(299)*UV(37)+JVS(300)*UV(40)+JVS(301)&
              &*UV(42)+JVS(302)*UV(43)+JVS(303)*UV(44)+JVS(304)*UV(46)+JVS(305)*UV(47)+JVS(306)*UV(50)+JVS(307)*UV(51)&
              &+JVS(308)*UV(52)+JVS(309)*UV(53)+JVS(310)*UV(55)+JVS(311)*UV(56)
  JUV(47) = JVS(313)*UV(1)+JVS(314)*UV(4)+JVS(315)*UV(16)+JVS(316)*UV(27)+JVS(317)*UV(30)+JVS(318)*UV(37)+JVS(319)&
              &*UV(43)+JVS(320)*UV(44)+JVS(321)*UV(47)+JVS(322)*UV(50)+JVS(323)*UV(51)+JVS(324)*UV(52)+JVS(325)*UV(53)&
              &+JVS(326)*UV(55)+JVS(327)*UV(56)+JVS(328)*UV(57)
  JUV(48) = JVS(329)*UV(38)+JVS(330)*UV(48)+JVS(331)*UV(50)+JVS(334)*UV(55)+JVS(335)*UV(56)+JVS(336)*UV(57)
  JUV(49) = JVS(337)*UV(14)+JVS(338)*UV(29)+JVS(339)*UV(38)+JVS(340)*UV(41)+JVS(342)*UV(49)+JVS(343)*UV(50)+JVS(345)&
              &*UV(53)+JVS(346)*UV(54)+JVS(347)*UV(55)+JVS(348)*UV(56)+JVS(349)*UV(57)
  JUV(50) = JVS(350)*UV(6)+JVS(351)*UV(9)+JVS(352)*UV(10)+JVS(353)*UV(12)+JVS(354)*UV(21)+JVS(355)*UV(24)+JVS(356)&
              &*UV(28)+JVS(357)*UV(34)+JVS(358)*UV(36)+JVS(359)*UV(37)+JVS(360)*UV(38)+JVS(361)*UV(43)+JVS(362)*UV(46)&
              &+JVS(363)*UV(47)+JVS(364)*UV(48)+JVS(365)*UV(49)+JVS(366)*UV(50)+JVS(367)*UV(51)+JVS(368)*UV(52)+JVS(369)&
              &*UV(53)+JVS(370)*UV(54)+JVS(371)*UV(55)+JVS(372)*UV(56)+JVS(373)*UV(57)
  JUV(51) = JVS(374)*UV(3)+JVS(375)*UV(9)+JVS(376)*UV(13)+JVS(377)*UV(20)+JVS(378)*UV(23)+JVS(379)*UV(28)+JVS(380)&
              &*UV(29)+JVS(381)*UV(30)+JVS(383)*UV(35)+JVS(384)*UV(36)+JVS(385)*UV(37)+JVS(386)*UV(38)+JVS(387)*UV(39)&
              &+JVS(389)*UV(41)+JVS(390)*UV(42)+JVS(393)*UV(46)+JVS(394)*UV(47)+JVS(396)*UV(49)+JVS(397)*UV(50)+JVS(398)&
              &*UV(51)+JVS(399)*UV(52)+JVS(400)*UV(53)+JVS(401)*UV(54)+JVS(402)*UV(55)+JVS(403)*UV(56)+JVS(404)*UV(57)
  JUV(52) = JVS(405)*UV(5)+JVS(406)*UV(6)+JVS(407)*UV(22)+JVS(408)*UV(25)+JVS(409)*UV(37)+JVS(410)*UV(38)+JVS(411)&
              &*UV(45)+JVS(412)*UV(46)+JVS(413)*UV(47)+JVS(416)*UV(50)+JVS(417)*UV(51)+JVS(418)*UV(52)+JVS(419)*UV(53)&
              &+JVS(420)*UV(54)+JVS(421)*UV(55)+JVS(422)*UV(56)+JVS(423)*UV(57)
  JUV(53) = JVS(424)*UV(19)+JVS(425)*UV(21)+JVS(427)*UV(37)+JVS(428)*UV(38)+JVS(429)*UV(41)+JVS(432)*UV(51)+JVS(433)&
              &*UV(52)+JVS(434)*UV(53)+JVS(435)*UV(54)+JVS(436)*UV(55)+JVS(437)*UV(56)+JVS(438)*UV(57)
  JUV(54) = JVS(439)*UV(5)+JVS(440)*UV(6)+JVS(441)*UV(7)+JVS(442)*UV(13)+JVS(443)*UV(14)+JVS(444)*UV(18)+JVS(445)*UV(22)&
              &+JVS(446)*UV(25)+JVS(447)*UV(33)+JVS(448)*UV(34)+JVS(449)*UV(36)+JVS(451)*UV(40)+JVS(452)*UV(42)+JVS(453)&
              &*UV(43)+JVS(454)*UV(44)+JVS(457)*UV(47)+JVS(458)*UV(48)+JVS(459)*UV(49)+JVS(460)*UV(50)+JVS(461)*UV(51)&
              &+JVS(462)*UV(52)+JVS(463)*UV(53)+JVS(464)*UV(54)+JVS(465)*UV(55)+JVS(466)*UV(56)+JVS(467)*UV(57)
  JUV(55) = JVS(468)*UV(1)+JVS(469)*UV(2)+JVS(470)*UV(4)+JVS(471)*UV(7)+JVS(472)*UV(8)+JVS(473)*UV(9)+JVS(474)*UV(10)&
              &+JVS(475)*UV(11)+JVS(476)*UV(12)+JVS(477)*UV(13)+JVS(478)*UV(14)+JVS(479)*UV(15)+JVS(480)*UV(16)+JVS(481)&
              &*UV(17)+JVS(482)*UV(18)+JVS(483)*UV(19)+JVS(484)*UV(20)+JVS(485)*UV(21)+JVS(486)*UV(22)+JVS(487)*UV(23)&
              &+JVS(488)*UV(24)+JVS(489)*UV(25)+JVS(490)*UV(26)+JVS(491)*UV(27)+JVS(492)*UV(28)+JVS(493)*UV(29)+JVS(494)&
              &*UV(30)+JVS(495)*UV(31)+JVS(496)*UV(32)+JVS(497)*UV(33)+JVS(499)*UV(35)+JVS(501)*UV(37)+JVS(502)*UV(38)&
              &+JVS(503)*UV(39)+JVS(505)*UV(41)+JVS(508)*UV(44)+JVS(509)*UV(45)+JVS(510)*UV(46)+JVS(516)*UV(52)+JVS(517)&
              &*UV(53)+JVS(518)*UV(54)+JVS(519)*UV(55)+JVS(520)*UV(56)+JVS(521)*UV(57)
  JUV(56) = JVS(522)*UV(34)+JVS(523)*UV(36)+JVS(524)*UV(40)+JVS(525)*UV(42)+JVS(526)*UV(43)+JVS(528)*UV(47)+JVS(529)&
              &*UV(48)+JVS(530)*UV(49)+JVS(531)*UV(50)+JVS(532)*UV(51)+JVS(533)*UV(52)+JVS(534)*UV(53)+JVS(535)*UV(54)&
              &+JVS(537)*UV(56)+JVS(538)*UV(57)
  JUV(57) = JVS(539)*UV(8)+JVS(540)*UV(11)+JVS(541)*UV(12)+JVS(542)*UV(15)+JVS(543)*UV(16)+JVS(544)*UV(17)+JVS(545)&
              &*UV(18)+JVS(546)*UV(19)+JVS(547)*UV(21)+JVS(548)*UV(22)+JVS(549)*UV(23)+JVS(550)*UV(24)+JVS(551)*UV(26)&
              &+JVS(552)*UV(27)+JVS(553)*UV(29)+JVS(554)*UV(31)+JVS(555)*UV(32)+JVS(556)*UV(33)+JVS(557)*UV(34)+JVS(558)&
              &*UV(35)+JVS(559)*UV(36)+JVS(560)*UV(37)+JVS(561)*UV(38)+JVS(562)*UV(39)+JVS(563)*UV(40)+JVS(564)*UV(41)&
              &+JVS(565)*UV(42)+JVS(566)*UV(43)+JVS(567)*UV(44)+JVS(568)*UV(45)+JVS(569)*UV(46)+JVS(570)*UV(47)+JVS(571)&
              &*UV(48)+JVS(572)*UV(49)+JVS(573)*UV(50)+JVS(574)*UV(51)+JVS(575)*UV(52)+JVS(576)*UV(53)+JVS(577)*UV(54)&
              &+JVS(578)*UV(55)+JVS(579)*UV(56)+JVS(580)*UV(57)
      
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

  JTUV(1) = JVS(1)*UV(1)+JVS(313)*UV(47)+JVS(468)*UV(55)
  JTUV(2) = JVS(3)*UV(2)+JVS(245)*UV(43)+JVS(469)*UV(55)
  JTUV(3) = JVS(5)*UV(3)+JVS(374)*UV(51)
  JTUV(4) = JVS(8)*UV(4)+JVS(162)*UV(34)+JVS(314)*UV(47)+JVS(470)*UV(55)
  JTUV(5) = JVS(10)*UV(5)+JVS(92)*UV(25)+JVS(405)*UV(52)+JVS(439)*UV(54)
  JTUV(6) = JVS(13)*UV(6)+JVS(263)*UV(45)+JVS(350)*UV(50)+JVS(406)*UV(52)+JVS(440)*UV(54)
  JTUV(7) = JVS(16)*UV(7)+JVS(127)*UV(31)+JVS(264)*UV(45)+JVS(441)*UV(54)+JVS(471)*UV(55)
  JTUV(8) = JVS(19)*UV(8)+JVS(472)*UV(55)+JVS(539)*UV(57)
  JTUV(9) = JVS(24)*UV(9)+JVS(351)*UV(50)+JVS(375)*UV(51)+JVS(473)*UV(55)
  JTUV(10) = JVS(28)*UV(10)+JVS(352)*UV(50)+JVS(474)*UV(55)
  JTUV(11) = JVS(34)*UV(11)+JVS(225)*UV(41)+JVS(265)*UV(45)+JVS(475)*UV(55)+JVS(540)*UV(57)
  JTUV(12) = JVS(38)*UV(12)+JVS(266)*UV(45)+JVS(353)*UV(50)+JVS(476)*UV(55)+JVS(541)*UV(57)
  JTUV(13) = JVS(42)*UV(13)+JVS(267)*UV(45)+JVS(376)*UV(51)+JVS(442)*UV(54)+JVS(477)*UV(55)
  JTUV(14) = JVS(46)*UV(14)+JVS(199)*UV(39)+JVS(337)*UV(49)+JVS(443)*UV(54)+JVS(478)*UV(55)
  JTUV(15) = JVS(50)*UV(15)+JVS(479)*UV(55)+JVS(542)*UV(57)
  JTUV(16) = JVS(57)*UV(16)+JVS(295)*UV(46)+JVS(315)*UV(47)+JVS(480)*UV(55)+JVS(543)*UV(57)
  JTUV(17) = JVS(61)*UV(17)+JVS(107)*UV(28)+JVS(163)*UV(34)+JVS(481)*UV(55)+JVS(544)*UV(57)
  JTUV(18) = JVS(65)*UV(18)+JVS(108)*UV(28)+JVS(444)*UV(54)+JVS(482)*UV(55)+JVS(545)*UV(57)
  JTUV(19) = JVS(51)*UV(15)+JVS(69)*UV(19)+JVS(128)*UV(31)+JVS(216)*UV(40)+JVS(268)*UV(45)+JVS(424)*UV(53)+JVS(483)&
               &*UV(55)+JVS(546)*UV(57)
  JTUV(20) = JVS(6)*UV(3)+JVS(72)*UV(20)+JVS(236)*UV(42)+JVS(296)*UV(46)+JVS(377)*UV(51)+JVS(484)*UV(55)
  JTUV(21) = JVS(76)*UV(21)+JVS(87)*UV(24)+JVS(269)*UV(45)+JVS(354)*UV(50)+JVS(425)*UV(53)+JVS(485)*UV(55)+JVS(547)&
               &*UV(57)
  JTUV(22) = JVS(79)*UV(22)+JVS(407)*UV(52)+JVS(445)*UV(54)+JVS(486)*UV(55)+JVS(548)*UV(57)
  JTUV(23) = JVS(83)*UV(23)+JVS(169)*UV(35)+JVS(183)*UV(36)+JVS(378)*UV(51)+JVS(487)*UV(55)+JVS(549)*UV(57)
  JTUV(24) = JVS(77)*UV(21)+JVS(88)*UV(24)+JVS(270)*UV(45)+JVS(355)*UV(50)+JVS(488)*UV(55)+JVS(550)*UV(57)
  JTUV(25) = JVS(93)*UV(25)+JVS(408)*UV(52)+JVS(446)*UV(54)+JVS(489)*UV(55)
  JTUV(26) = JVS(99)*UV(26)+JVS(200)*UV(39)+JVS(217)*UV(40)+JVS(271)*UV(45)+JVS(297)*UV(46)+JVS(490)*UV(55)+JVS(551)&
               &*UV(57)
  JTUV(27) = JVS(103)*UV(27)+JVS(119)*UV(30)+JVS(246)*UV(43)+JVS(298)*UV(46)+JVS(316)*UV(47)+JVS(491)*UV(55)+JVS(552)&
               &*UV(57)
  JTUV(28) = JVS(109)*UV(28)+JVS(184)*UV(36)+JVS(356)*UV(50)+JVS(379)*UV(51)+JVS(492)*UV(55)
  JTUV(29) = JVS(115)*UV(29)+JVS(129)*UV(31)+JVS(170)*UV(35)+JVS(201)*UV(39)+JVS(272)*UV(45)+JVS(338)*UV(49)+JVS(380)&
               &*UV(51)+JVS(493)*UV(55)+JVS(553)*UV(57)
  JTUV(30) = JVS(120)*UV(30)+JVS(237)*UV(42)+JVS(317)*UV(47)+JVS(381)*UV(51)+JVS(494)*UV(55)
  JTUV(31) = JVS(130)*UV(31)+JVS(495)*UV(55)+JVS(554)*UV(57)
  JTUV(32) = JVS(145)*UV(32)+JVS(273)*UV(45)+JVS(496)*UV(55)+JVS(555)*UV(57)
  JTUV(33) = JVS(17)*UV(7)+JVS(156)*UV(33)+JVS(202)*UV(39)+JVS(226)*UV(41)+JVS(274)*UV(45)+JVS(447)*UV(54)+JVS(497)&
               &*UV(55)+JVS(556)*UV(57)
  JTUV(34) = JVS(62)*UV(17)+JVS(66)*UV(18)+JVS(110)*UV(28)+JVS(146)*UV(32)+JVS(164)*UV(34)+JVS(275)*UV(45)+JVS(357)&
               &*UV(50)+JVS(448)*UV(54)+JVS(522)*UV(56)+JVS(557)*UV(57)
  JTUV(35) = JVS(132)*UV(31)+JVS(171)*UV(35)+JVS(383)*UV(51)+JVS(499)*UV(55)+JVS(558)*UV(57)
  JTUV(36) = JVS(84)*UV(23)+JVS(147)*UV(32)+JVS(172)*UV(35)+JVS(186)*UV(36)+JVS(203)*UV(39)+JVS(276)*UV(45)+JVS(358)&
               &*UV(50)+JVS(384)*UV(51)+JVS(449)*UV(54)+JVS(523)*UV(56)+JVS(559)*UV(57)
  JTUV(37) = JVS(52)*UV(15)+JVS(89)*UV(24)+JVS(133)*UV(31)+JVS(148)*UV(32)+JVS(173)*UV(35)+JVS(191)*UV(37)+JVS(218)&
               &*UV(40)+JVS(252)*UV(44)+JVS(277)*UV(45)+JVS(299)*UV(46)+JVS(318)*UV(47)+JVS(359)*UV(50)+JVS(385)*UV(51)&
               &+JVS(409)*UV(52)+JVS(427)*UV(53)+JVS(501)*UV(55)+JVS(560)*UV(57)
  JTUV(38) = JVS(20)*UV(8)+JVS(53)*UV(15)+JVS(134)*UV(31)+JVS(157)*UV(33)+JVS(195)*UV(38)+JVS(227)*UV(41)+JVS(278)&
               &*UV(45)+JVS(329)*UV(48)+JVS(339)*UV(49)+JVS(360)*UV(50)+JVS(386)*UV(51)+JVS(410)*UV(52)+JVS(428)*UV(53)&
               &+JVS(502)*UV(55)+JVS(561)*UV(57)
  JTUV(39) = JVS(174)*UV(35)+JVS(205)*UV(39)+JVS(279)*UV(45)+JVS(387)*UV(51)+JVS(503)*UV(55)+JVS(562)*UV(57)
  JTUV(40) = JVS(100)*UV(26)+JVS(219)*UV(40)+JVS(253)*UV(44)+JVS(280)*UV(45)+JVS(300)*UV(46)+JVS(451)*UV(54)+JVS(524)&
               &*UV(56)+JVS(563)*UV(57)
  JTUV(41) = JVS(54)*UV(15)+JVS(135)*UV(31)+JVS(175)*UV(35)+JVS(228)*UV(41)+JVS(281)*UV(45)+JVS(340)*UV(49)+JVS(389)&
               &*UV(51)+JVS(429)*UV(53)+JVS(505)*UV(55)+JVS(564)*UV(57)
  JTUV(42) = JVS(73)*UV(20)+JVS(238)*UV(42)+JVS(254)*UV(44)+JVS(301)*UV(46)+JVS(390)*UV(51)+JVS(452)*UV(54)+JVS(525)&
               &*UV(56)+JVS(565)*UV(57)
  JTUV(43) = JVS(104)*UV(27)+JVS(121)*UV(30)+JVS(149)*UV(32)+JVS(247)*UV(43)+JVS(255)*UV(44)+JVS(282)*UV(45)+JVS(302)&
               &*UV(46)+JVS(319)*UV(47)+JVS(361)*UV(50)+JVS(453)*UV(54)+JVS(526)*UV(56)+JVS(566)*UV(57)
  JTUV(44) = JVS(122)*UV(30)+JVS(256)*UV(44)+JVS(303)*UV(46)+JVS(320)*UV(47)+JVS(454)*UV(54)+JVS(508)*UV(55)+JVS(567)&
               &*UV(57)
  JTUV(45) = JVS(94)*UV(25)+JVS(136)*UV(31)+JVS(283)*UV(45)+JVS(411)*UV(52)+JVS(509)*UV(55)+JVS(568)*UV(57)
  JTUV(46) = JVS(95)*UV(25)+JVS(137)*UV(31)+JVS(304)*UV(46)+JVS(362)*UV(50)+JVS(393)*UV(51)+JVS(412)*UV(52)+JVS(510)&
               &*UV(55)+JVS(569)*UV(57)
  JTUV(47) = JVS(29)*UV(10)+JVS(58)*UV(16)+JVS(150)*UV(32)+JVS(284)*UV(45)+JVS(305)*UV(46)+JVS(321)*UV(47)+JVS(363)&
               &*UV(50)+JVS(394)*UV(51)+JVS(413)*UV(52)+JVS(457)*UV(54)+JVS(528)*UV(56)+JVS(570)*UV(57)
  JTUV(48) = JVS(35)*UV(11)+JVS(151)*UV(32)+JVS(158)*UV(33)+JVS(176)*UV(35)+JVS(207)*UV(39)+JVS(229)*UV(41)+JVS(285)&
               &*UV(45)+JVS(330)*UV(48)+JVS(364)*UV(50)+JVS(458)*UV(54)+JVS(529)*UV(56)+JVS(571)*UV(57)
  JTUV(49) = JVS(47)*UV(14)+JVS(116)*UV(29)+JVS(138)*UV(31)+JVS(152)*UV(32)+JVS(177)*UV(35)+JVS(208)*UV(39)+JVS(286)&
               &*UV(45)+JVS(342)*UV(49)+JVS(365)*UV(50)+JVS(396)*UV(51)+JVS(459)*UV(54)+JVS(530)*UV(56)+JVS(572)*UV(57)
  JTUV(50) = JVS(14)*UV(6)+JVS(30)*UV(10)+JVS(39)*UV(12)+JVS(111)*UV(28)+JVS(123)*UV(30)+JVS(139)*UV(31)+JVS(153)*UV(32)&
               &+JVS(165)*UV(34)+JVS(178)*UV(35)+JVS(187)*UV(36)+JVS(209)*UV(39)+JVS(230)*UV(41)+JVS(248)*UV(43)+JVS(287)&
               &*UV(45)+JVS(306)*UV(46)+JVS(322)*UV(47)+JVS(331)*UV(48)+JVS(343)*UV(49)+JVS(366)*UV(50)+JVS(397)*UV(51)&
               &+JVS(416)*UV(52)+JVS(460)*UV(54)+JVS(531)*UV(56)+JVS(573)*UV(57)
  JTUV(51) = JVS(25)*UV(9)+JVS(31)*UV(10)+JVS(43)*UV(13)+JVS(288)*UV(45)+JVS(307)*UV(46)+JVS(323)*UV(47)+JVS(367)*UV(50)&
               &+JVS(398)*UV(51)+JVS(417)*UV(52)+JVS(432)*UV(53)+JVS(461)*UV(54)+JVS(532)*UV(56)+JVS(574)*UV(57)
  JTUV(52) = JVS(11)*UV(5)+JVS(96)*UV(25)+JVS(140)*UV(31)+JVS(159)*UV(33)+JVS(192)*UV(37)+JVS(196)*UV(38)+JVS(258)&
               &*UV(44)+JVS(289)*UV(45)+JVS(308)*UV(46)+JVS(324)*UV(47)+JVS(368)*UV(50)+JVS(399)*UV(51)+JVS(418)*UV(52)&
               &+JVS(433)*UV(53)+JVS(462)*UV(54)+JVS(516)*UV(55)+JVS(533)*UV(56)+JVS(575)*UV(57)
  JTUV(53) = JVS(21)*UV(8)+JVS(55)*UV(15)+JVS(70)*UV(19)+JVS(78)*UV(21)+JVS(90)*UV(24)+JVS(141)*UV(31)+JVS(154)*UV(32)&
               &+JVS(179)*UV(35)+JVS(193)*UV(37)+JVS(197)*UV(38)+JVS(232)*UV(41)+JVS(290)*UV(45)+JVS(309)*UV(46)+JVS(325)&
               &*UV(47)+JVS(345)*UV(49)+JVS(369)*UV(50)+JVS(400)*UV(51)+JVS(419)*UV(52)+JVS(434)*UV(53)+JVS(463)*UV(54)&
               &+JVS(517)*UV(55)+JVS(534)*UV(56)+JVS(576)*UV(57)
  JTUV(54) = JVS(12)*UV(5)+JVS(15)*UV(6)+JVS(44)*UV(13)+JVS(48)*UV(14)+JVS(80)*UV(22)+JVS(97)*UV(25)+JVS(346)*UV(49)&
               &+JVS(370)*UV(50)+JVS(401)*UV(51)+JVS(420)*UV(52)+JVS(435)*UV(53)+JVS(464)*UV(54)+JVS(518)*UV(55)+JVS(535)&
               &*UV(56)+JVS(577)*UV(57)
  JTUV(55) = JVS(2)*UV(1)+JVS(4)*UV(2)+JVS(7)*UV(3)+JVS(9)*UV(4)+JVS(18)*UV(7)+JVS(22)*UV(8)+JVS(26)*UV(9)+JVS(32)&
               &*UV(10)+JVS(36)*UV(11)+JVS(40)*UV(12)+JVS(45)*UV(13)+JVS(49)*UV(14)+JVS(56)*UV(15)+JVS(59)*UV(16)+JVS(63)&
               &*UV(17)+JVS(67)*UV(18)+JVS(71)*UV(19)+JVS(74)*UV(20)+JVS(81)*UV(22)+JVS(85)*UV(23)+JVS(91)*UV(24)+JVS(98)&
               &*UV(25)+JVS(101)*UV(26)+JVS(105)*UV(27)+JVS(112)*UV(28)+JVS(117)*UV(29)+JVS(124)*UV(30)+JVS(142)*UV(31)&
               &+JVS(155)*UV(32)+JVS(160)*UV(33)+JVS(166)*UV(34)+JVS(180)*UV(35)+JVS(188)*UV(36)+JVS(194)*UV(37)+JVS(198)&
               &*UV(38)+JVS(213)*UV(39)+JVS(222)*UV(40)+JVS(233)*UV(41)+JVS(242)*UV(42)+JVS(249)*UV(43)+JVS(260)*UV(44)&
               &+JVS(292)*UV(45)+JVS(310)*UV(46)+JVS(326)*UV(47)+JVS(334)*UV(48)+JVS(347)*UV(49)+JVS(371)*UV(50)+JVS(402)&
               &*UV(51)+JVS(421)*UV(52)+JVS(436)*UV(53)+JVS(465)*UV(54)+JVS(519)*UV(55)+JVS(578)*UV(57)
  JTUV(56) = JVS(68)*UV(18)+JVS(113)*UV(28)+JVS(125)*UV(30)+JVS(143)*UV(31)+JVS(161)*UV(33)+JVS(167)*UV(34)+JVS(181)&
               &*UV(35)+JVS(189)*UV(36)+JVS(214)*UV(39)+JVS(223)*UV(40)+JVS(234)*UV(41)+JVS(243)*UV(42)+JVS(250)*UV(43)&
               &+JVS(261)*UV(44)+JVS(293)*UV(45)+JVS(311)*UV(46)+JVS(327)*UV(47)+JVS(335)*UV(48)+JVS(348)*UV(49)+JVS(372)&
               &*UV(50)+JVS(403)*UV(51)+JVS(422)*UV(52)+JVS(437)*UV(53)+JVS(466)*UV(54)+JVS(520)*UV(55)+JVS(537)*UV(56)&
               &+JVS(579)*UV(57)
  JTUV(57) = JVS(23)*UV(8)+JVS(27)*UV(9)+JVS(33)*UV(10)+JVS(37)*UV(11)+JVS(41)*UV(12)+JVS(60)*UV(16)+JVS(64)*UV(17)&
               &+JVS(75)*UV(20)+JVS(82)*UV(22)+JVS(86)*UV(23)+JVS(102)*UV(26)+JVS(106)*UV(27)+JVS(118)*UV(29)+JVS(168)&
               &*UV(34)+JVS(190)*UV(36)+JVS(224)*UV(40)+JVS(244)*UV(42)+JVS(251)*UV(43)+JVS(294)*UV(45)+JVS(328)*UV(47)&
               &+JVS(336)*UV(48)+JVS(349)*UV(49)+JVS(373)*UV(50)+JVS(404)*UV(51)+JVS(423)*UV(52)+JVS(438)*UV(53)+JVS(467)&
               &*UV(54)+JVS(521)*UV(55)+JVS(538)*UV(56)+JVS(580)*UV(57)
      
END SUBROUTINE JacTR_SP_Vec

! End of JacTR_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE nmhc9_Jacobian

