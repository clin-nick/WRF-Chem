























MODULE saprc99_Jacobian

  USE saprc99_Parameters
  USE saprc99_JacobianSP

  IMPLICIT NONE

CONTAINS












SUBROUTINE Jac_SP_Vec ( JVS, UV, JUV )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: UV(NVAR)

  REAL(kind=dp) :: JUV(NVAR)

  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(10)+JVS(3)*UV(71)
  JUV(2) = JVS(4)*UV(2)+JVS(5)*UV(39)+JVS(6)*UV(50)+JVS(7)*UV(67)
  JUV(3) = JVS(8)*UV(3)+JVS(9)*UV(30)+JVS(10)*UV(39)+JVS(11)*UV(48)+JVS(12)*UV(50)+JVS(13)*UV(52)+JVS(14)*UV(54)+JVS(15)&
             &*UV(55)+JVS(16)*UV(56)+JVS(17)*UV(57)+JVS(18)*UV(58)+JVS(19)*UV(59)+JVS(20)*UV(67)+JVS(21)*UV(71)+JVS(22)&
             &*UV(73)
  JUV(4) = JVS(23)*UV(4)+JVS(24)*UV(50)+JVS(25)*UV(56)+JVS(26)*UV(58)+JVS(27)*UV(67)+JVS(28)*UV(68)+JVS(29)*UV(69)&
             &+JVS(30)*UV(76)+JVS(31)*UV(77)+JVS(32)*UV(78)
  JUV(5) = JVS(33)*UV(5)+JVS(34)*UV(52)+JVS(35)*UV(54)+JVS(36)*UV(56)+JVS(37)*UV(57)+JVS(38)*UV(58)+JVS(39)*UV(67)&
             &+JVS(40)*UV(69)+JVS(41)*UV(70)+JVS(42)*UV(72)+JVS(43)*UV(76)+JVS(44)*UV(77)+JVS(45)*UV(78)+JVS(46)*UV(79)
  JUV(6) = JVS(47)*UV(6)+JVS(48)*UV(68)+JVS(49)*UV(78)
  JUV(7) = JVS(50)*UV(7)+JVS(51)*UV(70)+JVS(52)*UV(72)+JVS(53)*UV(78)+JVS(54)*UV(79)
  JUV(8) = JVS(55)*UV(8)+JVS(56)*UV(31)+JVS(57)*UV(50)+JVS(58)*UV(74)+JVS(59)*UV(75)
  JUV(9) = JVS(60)*UV(9)+JVS(61)*UV(31)+JVS(62)*UV(41)+JVS(63)*UV(50)+JVS(64)*UV(63)+JVS(65)*UV(67)+JVS(66)*UV(71)&
             &+JVS(67)*UV(74)+JVS(68)*UV(75)
  JUV(10) = JVS(69)*UV(10)+JVS(70)*UV(71)
  JUV(11) = JVS(71)*UV(11)+JVS(72)*UV(67)
  JUV(12) = JVS(73)*UV(12)+JVS(74)*UV(71)
  JUV(13) = JVS(75)*UV(13)+JVS(76)*UV(71)
  JUV(14) = JVS(77)*UV(14)+JVS(78)*UV(27)+JVS(79)*UV(54)+JVS(80)*UV(67)+JVS(81)*UV(71)
  JUV(15) = JVS(82)*UV(15)+JVS(83)*UV(68)+JVS(84)*UV(74)
  JUV(16) = JVS(85)*UV(16)+JVS(86)*UV(72)+JVS(87)*UV(74)
  JUV(17) = JVS(88)*UV(17)+JVS(89)*UV(74)+JVS(90)*UV(79)
  JUV(18) = JVS(91)*UV(18)+JVS(92)*UV(70)+JVS(93)*UV(74)
  JUV(19) = JVS(94)*UV(19)+JVS(95)*UV(71)+JVS(96)*UV(78)
  JUV(20) = JVS(97)*UV(20)+JVS(98)*UV(71)
  JUV(21) = JVS(99)*UV(21)+JVS(100)*UV(71)
  JUV(22) = JVS(101)*UV(22)+JVS(102)*UV(74)+JVS(103)*UV(75)
  JUV(23) = JVS(104)*UV(23)+JVS(105)*UV(71)+JVS(106)*UV(73)
  JUV(24) = JVS(107)*UV(24)+JVS(108)*UV(71)
  JUV(25) = JVS(109)*UV(24)+JVS(110)*UV(25)+JVS(111)*UV(71)+JVS(112)*UV(74)
  JUV(26) = JVS(113)*UV(26)+JVS(114)*UV(71)
  JUV(27) = JVS(115)*UV(27)+JVS(116)*UV(71)
  JUV(28) = JVS(117)*UV(28)+JVS(118)*UV(69)+JVS(119)*UV(71)+JVS(120)*UV(76)+JVS(121)*UV(77)
  JUV(29) = JVS(122)*UV(29)+JVS(123)*UV(71)+JVS(124)*UV(76)+JVS(125)*UV(78)
  JUV(30) = JVS(126)*UV(30)+JVS(127)*UV(61)+JVS(128)*UV(73)+JVS(129)*UV(78)
  JUV(31) = JVS(130)*UV(31)+JVS(131)*UV(40)+JVS(132)*UV(74)+JVS(133)*UV(75)+JVS(134)*UV(78)
  JUV(32) = JVS(135)*UV(32)+JVS(136)*UV(71)+JVS(137)*UV(74)+JVS(138)*UV(78)
  JUV(33) = JVS(139)*UV(33)+JVS(140)*UV(71)
  JUV(34) = JVS(141)*UV(34)+JVS(142)*UV(71)
  JUV(35) = JVS(143)*UV(27)+JVS(144)*UV(34)+JVS(145)*UV(35)+JVS(146)*UV(71)
  JUV(36) = JVS(147)*UV(27)+JVS(148)*UV(34)+JVS(149)*UV(36)+JVS(150)*UV(71)
  JUV(37) = JVS(151)*UV(27)+JVS(152)*UV(34)+JVS(153)*UV(37)+JVS(154)*UV(71)+JVS(155)*UV(75)
  JUV(38) = JVS(156)*UV(27)+JVS(157)*UV(34)+JVS(158)*UV(38)+JVS(159)*UV(67)+JVS(160)*UV(71)
  JUV(39) = JVS(161)*UV(39)+JVS(162)*UV(67)+JVS(163)*UV(71)
  JUV(40) = JVS(164)*UV(31)+JVS(165)*UV(40)+JVS(166)*UV(51)+JVS(167)*UV(74)+JVS(168)*UV(75)+JVS(169)*UV(78)
  JUV(41) = JVS(170)*UV(27)+JVS(171)*UV(34)+JVS(172)*UV(41)+JVS(173)*UV(58)+JVS(174)*UV(67)+JVS(175)*UV(71)+JVS(176)&
              &*UV(75)
  JUV(42) = JVS(177)*UV(42)+JVS(178)*UV(69)+JVS(179)*UV(71)+JVS(180)*UV(77)+JVS(181)*UV(78)
  JUV(43) = JVS(182)*UV(34)+JVS(183)*UV(43)+JVS(184)*UV(51)+JVS(185)*UV(71)+JVS(186)*UV(75)+JVS(187)*UV(78)
  JUV(44) = JVS(188)*UV(27)+JVS(189)*UV(34)+JVS(190)*UV(35)+JVS(191)*UV(36)+JVS(192)*UV(37)+JVS(193)*UV(44)+JVS(194)&
              &*UV(55)+JVS(195)*UV(57)+JVS(196)*UV(59)+JVS(197)*UV(67)+JVS(198)*UV(71)+JVS(199)*UV(75)
  JUV(45) = JVS(200)*UV(33)+JVS(201)*UV(35)+JVS(202)*UV(36)+JVS(203)*UV(38)+JVS(204)*UV(39)+JVS(205)*UV(44)+JVS(206)&
              &*UV(45)+JVS(207)*UV(48)+JVS(208)*UV(49)+JVS(209)*UV(50)+JVS(210)*UV(52)+JVS(211)*UV(54)+JVS(212)*UV(55)&
              &+JVS(213)*UV(56)+JVS(214)*UV(57)+JVS(215)*UV(58)+JVS(216)*UV(59)+JVS(217)*UV(60)+JVS(218)*UV(61)+JVS(219)&
              &*UV(63)+JVS(220)*UV(64)+JVS(221)*UV(67)+JVS(222)*UV(71)+JVS(223)*UV(75)
  JUV(46) = JVS(224)*UV(20)+JVS(225)*UV(24)+JVS(226)*UV(25)+JVS(227)*UV(26)+JVS(228)*UV(33)+JVS(229)*UV(46)+JVS(230)&
              &*UV(54)+JVS(231)*UV(56)+JVS(232)*UV(58)+JVS(233)*UV(62)+JVS(234)*UV(67)+JVS(235)*UV(71)+JVS(237)*UV(75)
  JUV(47) = JVS(238)*UV(22)+JVS(239)*UV(37)+JVS(240)*UV(40)+JVS(241)*UV(41)+JVS(242)*UV(43)+JVS(243)*UV(44)+JVS(244)&
              &*UV(47)+JVS(245)*UV(49)+JVS(247)*UV(55)+JVS(248)*UV(57)+JVS(251)*UV(60)+JVS(252)*UV(61)+JVS(253)*UV(64)&
              &+JVS(255)*UV(71)+JVS(256)*UV(74)+JVS(257)*UV(75)+JVS(258)*UV(78)
  JUV(48) = JVS(259)*UV(48)+JVS(260)*UV(63)+JVS(261)*UV(67)+JVS(262)*UV(71)+JVS(263)*UV(75)
  JUV(49) = JVS(264)*UV(27)+JVS(265)*UV(34)+JVS(266)*UV(35)+JVS(267)*UV(36)+JVS(268)*UV(38)+JVS(269)*UV(39)+JVS(270)&
              &*UV(43)+JVS(271)*UV(48)+JVS(272)*UV(49)+JVS(274)*UV(54)+JVS(275)*UV(57)+JVS(276)*UV(63)+JVS(277)*UV(67)&
              &+JVS(278)*UV(71)+JVS(279)*UV(75)
  JUV(50) = JVS(281)*UV(50)+JVS(282)*UV(63)+JVS(283)*UV(67)+JVS(284)*UV(71)+JVS(285)*UV(75)
  JUV(51) = JVS(286)*UV(37)+JVS(287)*UV(43)+JVS(288)*UV(51)+JVS(289)*UV(68)+JVS(290)*UV(70)+JVS(291)*UV(71)+JVS(292)&
              &*UV(72)+JVS(293)*UV(73)+JVS(294)*UV(74)+JVS(295)*UV(75)+JVS(296)*UV(78)+JVS(297)*UV(79)
  JUV(52) = JVS(298)*UV(52)+JVS(299)*UV(63)+JVS(300)*UV(67)+JVS(301)*UV(71)+JVS(302)*UV(75)
  JUV(53) = JVS(303)*UV(24)+JVS(304)*UV(26)+JVS(305)*UV(33)+JVS(306)*UV(35)+JVS(307)*UV(36)+JVS(308)*UV(46)+JVS(309)&
              &*UV(52)+JVS(310)*UV(53)+JVS(311)*UV(54)+JVS(312)*UV(56)+JVS(313)*UV(58)+JVS(314)*UV(59)+JVS(315)*UV(62)&
              &+JVS(316)*UV(63)+JVS(317)*UV(65)+JVS(318)*UV(66)+JVS(319)*UV(67)+JVS(320)*UV(68)+JVS(321)*UV(69)+JVS(322)&
              &*UV(70)+JVS(323)*UV(71)+JVS(324)*UV(72)+JVS(325)*UV(73)+JVS(327)*UV(75)+JVS(328)*UV(76)+JVS(329)*UV(77)&
              &+JVS(330)*UV(78)+JVS(331)*UV(79)
  JUV(54) = JVS(332)*UV(54)+JVS(333)*UV(63)+JVS(334)*UV(67)+JVS(335)*UV(71)+JVS(336)*UV(75)
  JUV(55) = JVS(337)*UV(52)+JVS(338)*UV(55)+JVS(339)*UV(58)+JVS(340)*UV(63)+JVS(341)*UV(67)+JVS(342)*UV(71)+JVS(343)&
              &*UV(75)
  JUV(56) = JVS(344)*UV(56)+JVS(345)*UV(63)+JVS(346)*UV(67)+JVS(347)*UV(71)+JVS(348)*UV(75)
  JUV(57) = JVS(349)*UV(52)+JVS(350)*UV(57)+JVS(351)*UV(58)+JVS(353)*UV(67)+JVS(354)*UV(71)+JVS(355)*UV(75)
  JUV(58) = JVS(356)*UV(58)+JVS(357)*UV(63)+JVS(358)*UV(67)+JVS(359)*UV(71)+JVS(360)*UV(75)
  JUV(59) = JVS(361)*UV(52)+JVS(362)*UV(58)+JVS(363)*UV(59)+JVS(364)*UV(63)+JVS(365)*UV(67)+JVS(366)*UV(71)+JVS(367)&
              &*UV(75)
  JUV(60) = JVS(368)*UV(13)+JVS(369)*UV(21)+JVS(370)*UV(24)+JVS(371)*UV(26)+JVS(372)*UV(33)+JVS(373)*UV(48)+JVS(374)&
              &*UV(50)+JVS(375)*UV(56)+JVS(376)*UV(57)+JVS(377)*UV(58)+JVS(378)*UV(60)+JVS(379)*UV(62)+JVS(380)*UV(63)&
              &+JVS(381)*UV(64)+JVS(382)*UV(65)+JVS(383)*UV(66)+JVS(384)*UV(67)+JVS(385)*UV(68)+JVS(386)*UV(70)+JVS(387)&
              &*UV(71)+JVS(388)*UV(72)+JVS(389)*UV(73)+JVS(390)*UV(75)+JVS(391)*UV(79)
  JUV(61) = JVS(392)*UV(21)+JVS(393)*UV(24)+JVS(394)*UV(26)+JVS(395)*UV(28)+JVS(396)*UV(29)+JVS(397)*UV(30)+JVS(398)&
              &*UV(33)+JVS(399)*UV(39)+JVS(400)*UV(46)+JVS(401)*UV(48)+JVS(402)*UV(49)+JVS(403)*UV(50)+JVS(405)*UV(52)&
              &+JVS(406)*UV(54)+JVS(407)*UV(55)+JVS(408)*UV(56)+JVS(409)*UV(57)+JVS(410)*UV(58)+JVS(411)*UV(59)+JVS(412)&
              &*UV(61)+JVS(413)*UV(62)+JVS(414)*UV(63)+JVS(415)*UV(65)+JVS(416)*UV(66)+JVS(417)*UV(67)+JVS(418)*UV(68)&
              &+JVS(419)*UV(69)+JVS(420)*UV(70)+JVS(421)*UV(71)+JVS(422)*UV(72)+JVS(423)*UV(73)+JVS(425)*UV(75)+JVS(426)&
              &*UV(76)+JVS(427)*UV(77)+JVS(428)*UV(78)+JVS(429)*UV(79)
  JUV(62) = JVS(430)*UV(25)+JVS(431)*UV(54)+JVS(432)*UV(56)+JVS(433)*UV(57)+JVS(434)*UV(58)+JVS(435)*UV(62)+JVS(438)&
              &*UV(69)+JVS(439)*UV(71)+JVS(440)*UV(73)+JVS(441)*UV(74)+JVS(442)*UV(75)
  JUV(63) = JVS(443)*UV(11)+JVS(444)*UV(48)+JVS(445)*UV(50)+JVS(446)*UV(52)+JVS(447)*UV(54)+JVS(448)*UV(55)+JVS(449)&
              &*UV(56)+JVS(450)*UV(58)+JVS(451)*UV(59)+JVS(452)*UV(63)+JVS(453)*UV(67)+JVS(455)*UV(73)+JVS(456)*UV(74)&
              &+JVS(457)*UV(75)
  JUV(64) = JVS(458)*UV(20)+JVS(459)*UV(24)+JVS(460)*UV(26)+JVS(461)*UV(33)+JVS(462)*UV(35)+JVS(463)*UV(36)+JVS(464)&
              &*UV(38)+JVS(465)*UV(42)+JVS(466)*UV(48)+JVS(467)*UV(50)+JVS(468)*UV(54)+JVS(469)*UV(55)+JVS(470)*UV(56)&
              &+JVS(471)*UV(57)+JVS(472)*UV(58)+JVS(473)*UV(59)+JVS(474)*UV(62)+JVS(475)*UV(63)+JVS(476)*UV(64)+JVS(477)&
              &*UV(65)+JVS(478)*UV(66)+JVS(479)*UV(67)+JVS(481)*UV(71)+JVS(484)*UV(75)
  JUV(65) = JVS(487)*UV(24)+JVS(488)*UV(26)+JVS(489)*UV(33)+JVS(490)*UV(50)+JVS(491)*UV(55)+JVS(492)*UV(56)+JVS(493)&
              &*UV(57)+JVS(494)*UV(58)+JVS(495)*UV(59)+JVS(496)*UV(62)+JVS(497)*UV(63)+JVS(498)*UV(65)+JVS(499)*UV(66)&
              &+JVS(500)*UV(67)+JVS(501)*UV(69)+JVS(502)*UV(71)+JVS(505)*UV(75)+JVS(506)*UV(76)+JVS(507)*UV(77)
  JUV(66) = JVS(508)*UV(26)+JVS(509)*UV(33)+JVS(510)*UV(34)+JVS(511)*UV(52)+JVS(512)*UV(54)+JVS(513)*UV(56)+JVS(514)&
              &*UV(57)+JVS(515)*UV(58)+JVS(516)*UV(59)+JVS(517)*UV(62)+JVS(518)*UV(63)+JVS(519)*UV(66)+JVS(520)*UV(67)&
              &+JVS(521)*UV(68)+JVS(522)*UV(69)+JVS(523)*UV(71)+JVS(524)*UV(72)+JVS(528)*UV(76)+JVS(529)*UV(77)+JVS(530)&
              &*UV(79)
  JUV(67) = JVS(531)*UV(38)+JVS(532)*UV(39)+JVS(533)*UV(48)+JVS(534)*UV(50)+JVS(535)*UV(52)+JVS(536)*UV(54)+JVS(537)&
              &*UV(55)+JVS(538)*UV(56)+JVS(539)*UV(57)+JVS(540)*UV(58)+JVS(541)*UV(59)+JVS(542)*UV(63)+JVS(543)*UV(67)&
              &+JVS(544)*UV(68)+JVS(545)*UV(70)+JVS(546)*UV(71)+JVS(547)*UV(72)+JVS(548)*UV(73)+JVS(549)*UV(74)+JVS(551)&
              &*UV(78)+JVS(552)*UV(79)
  JUV(68) = JVS(553)*UV(14)+JVS(554)*UV(15)+JVS(556)*UV(33)+JVS(557)*UV(35)+JVS(558)*UV(36)+JVS(559)*UV(44)+JVS(560)&
              &*UV(46)+JVS(561)*UV(54)+JVS(562)*UV(55)+JVS(564)*UV(57)+JVS(565)*UV(58)+JVS(566)*UV(59)+JVS(567)*UV(60)&
              &+JVS(571)*UV(65)+JVS(572)*UV(66)+JVS(573)*UV(67)+JVS(574)*UV(68)+JVS(575)*UV(69)+JVS(576)*UV(70)+JVS(577)&
              &*UV(71)+JVS(578)*UV(72)+JVS(579)*UV(73)+JVS(580)*UV(74)+JVS(581)*UV(75)+JVS(582)*UV(76)+JVS(583)*UV(77)&
              &+JVS(584)*UV(78)+JVS(585)*UV(79)
  JUV(69) = JVS(586)*UV(20)+JVS(587)*UV(24)+JVS(588)*UV(26)+JVS(589)*UV(27)+JVS(590)*UV(33)+JVS(591)*UV(34)+JVS(592)&
              &*UV(50)+JVS(593)*UV(52)+JVS(594)*UV(54)+JVS(595)*UV(56)+JVS(596)*UV(57)+JVS(597)*UV(58)+JVS(598)*UV(59)&
              &+JVS(599)*UV(62)+JVS(600)*UV(63)+JVS(601)*UV(64)+JVS(602)*UV(65)+JVS(603)*UV(66)+JVS(604)*UV(67)+JVS(605)&
              &*UV(68)+JVS(606)*UV(69)+JVS(607)*UV(70)+JVS(608)*UV(71)+JVS(609)*UV(72)+JVS(610)*UV(73)+JVS(612)*UV(75)&
              &+JVS(613)*UV(76)+JVS(614)*UV(77)+JVS(615)*UV(78)+JVS(616)*UV(79)
  JUV(70) = JVS(617)*UV(18)+JVS(618)*UV(52)+JVS(619)*UV(55)+JVS(620)*UV(57)+JVS(622)*UV(59)+JVS(623)*UV(63)+JVS(624)&
              &*UV(67)+JVS(625)*UV(68)+JVS(626)*UV(69)+JVS(627)*UV(70)+JVS(628)*UV(71)+JVS(629)*UV(72)+JVS(630)*UV(73)&
              &+JVS(631)*UV(74)+JVS(632)*UV(75)+JVS(633)*UV(76)+JVS(634)*UV(77)+JVS(635)*UV(78)+JVS(636)*UV(79)
  JUV(71) = JVS(637)*UV(10)+JVS(638)*UV(11)+JVS(639)*UV(12)+JVS(640)*UV(13)+JVS(641)*UV(19)+JVS(642)*UV(20)+JVS(643)&
              &*UV(21)+JVS(644)*UV(23)+JVS(645)*UV(24)+JVS(646)*UV(26)+JVS(647)*UV(27)+JVS(648)*UV(28)+JVS(649)*UV(29)&
              &+JVS(650)*UV(32)+JVS(651)*UV(33)+JVS(652)*UV(34)+JVS(653)*UV(35)+JVS(654)*UV(36)+JVS(655)*UV(37)+JVS(656)&
              &*UV(38)+JVS(657)*UV(39)+JVS(658)*UV(41)+JVS(659)*UV(42)+JVS(660)*UV(43)+JVS(661)*UV(44)+JVS(662)*UV(45)&
              &+JVS(663)*UV(46)+JVS(664)*UV(47)+JVS(665)*UV(48)+JVS(666)*UV(49)+JVS(667)*UV(50)+JVS(669)*UV(52)+JVS(670)&
              &*UV(54)+JVS(671)*UV(55)+JVS(672)*UV(56)+JVS(673)*UV(57)+JVS(674)*UV(58)+JVS(675)*UV(59)+JVS(676)*UV(60)&
              &+JVS(677)*UV(61)+JVS(678)*UV(62)+JVS(680)*UV(64)+JVS(681)*UV(65)+JVS(682)*UV(66)+JVS(683)*UV(67)+JVS(687)&
              &*UV(71)+JVS(689)*UV(73)+JVS(690)*UV(74)+JVS(691)*UV(75)+JVS(694)*UV(78)
  JUV(72) = JVS(696)*UV(16)+JVS(697)*UV(49)+JVS(699)*UV(54)+JVS(700)*UV(55)+JVS(701)*UV(57)+JVS(702)*UV(58)+JVS(703)&
              &*UV(59)+JVS(705)*UV(64)+JVS(706)*UV(65)+JVS(707)*UV(66)+JVS(708)*UV(67)+JVS(709)*UV(68)+JVS(710)*UV(69)&
              &+JVS(711)*UV(70)+JVS(712)*UV(71)+JVS(713)*UV(72)+JVS(714)*UV(73)+JVS(715)*UV(74)+JVS(716)*UV(75)+JVS(717)&
              &*UV(76)+JVS(718)*UV(77)+JVS(719)*UV(78)+JVS(720)*UV(79)
  JUV(73) = JVS(721)*UV(23)+JVS(722)*UV(30)+JVS(723)*UV(53)+JVS(730)*UV(63)+JVS(733)*UV(67)+JVS(734)*UV(68)+JVS(735)&
              &*UV(69)+JVS(736)*UV(70)+JVS(737)*UV(71)+JVS(738)*UV(72)+JVS(739)*UV(73)+JVS(740)*UV(74)+JVS(741)*UV(75)&
              &+JVS(742)*UV(76)+JVS(743)*UV(77)+JVS(744)*UV(78)+JVS(745)*UV(79)
  JUV(74) = JVS(746)*UV(15)+JVS(747)*UV(16)+JVS(748)*UV(17)+JVS(749)*UV(18)+JVS(750)*UV(22)+JVS(751)*UV(23)+JVS(752)&
              &*UV(25)+JVS(753)*UV(30)+JVS(754)*UV(31)+JVS(755)*UV(32)+JVS(757)*UV(47)+JVS(759)*UV(51)+JVS(760)*UV(52)&
              &+JVS(761)*UV(53)+JVS(762)*UV(54)+JVS(766)*UV(58)+JVS(770)*UV(62)+JVS(771)*UV(63)+JVS(775)*UV(67)+JVS(776)&
              &*UV(68)+JVS(777)*UV(69)+JVS(778)*UV(70)+JVS(779)*UV(71)+JVS(780)*UV(72)+JVS(781)*UV(73)+JVS(782)*UV(74)&
              &+JVS(783)*UV(75)+JVS(784)*UV(76)+JVS(785)*UV(77)+JVS(786)*UV(78)+JVS(787)*UV(79)
  JUV(75) = JVS(788)*UV(22)+JVS(789)*UV(32)+JVS(790)*UV(37)+JVS(791)*UV(40)+JVS(792)*UV(41)+JVS(793)*UV(43)+JVS(794)&
              &*UV(44)+JVS(795)*UV(47)+JVS(796)*UV(48)+JVS(797)*UV(49)+JVS(798)*UV(50)+JVS(800)*UV(52)+JVS(801)*UV(53)&
              &+JVS(802)*UV(54)+JVS(803)*UV(55)+JVS(804)*UV(56)+JVS(805)*UV(57)+JVS(806)*UV(58)+JVS(808)*UV(60)+JVS(809)&
              &*UV(61)+JVS(811)*UV(63)+JVS(812)*UV(64)+JVS(815)*UV(67)+JVS(816)*UV(68)+JVS(817)*UV(69)+JVS(818)*UV(70)&
              &+JVS(819)*UV(71)+JVS(820)*UV(72)+JVS(821)*UV(73)+JVS(822)*UV(74)+JVS(823)*UV(75)+JVS(824)*UV(76)+JVS(825)&
              &*UV(77)+JVS(826)*UV(78)+JVS(827)*UV(79)
  JUV(76) = JVS(828)*UV(12)+JVS(829)*UV(25)+JVS(830)*UV(29)+JVS(831)*UV(33)+JVS(832)*UV(46)+JVS(833)*UV(48)+JVS(834)&
              &*UV(50)+JVS(835)*UV(52)+JVS(837)*UV(56)+JVS(838)*UV(58)+JVS(839)*UV(59)+JVS(840)*UV(60)+JVS(842)*UV(63)&
              &+JVS(846)*UV(67)+JVS(847)*UV(68)+JVS(848)*UV(69)+JVS(849)*UV(70)+JVS(850)*UV(71)+JVS(851)*UV(72)+JVS(852)&
              &*UV(73)+JVS(854)*UV(75)+JVS(855)*UV(76)+JVS(856)*UV(77)+JVS(857)*UV(78)+JVS(858)*UV(79)
  JUV(77) = JVS(859)*UV(13)+JVS(860)*UV(20)+JVS(861)*UV(21)+JVS(862)*UV(24)+JVS(863)*UV(26)+JVS(864)*UV(27)+JVS(865)&
              &*UV(33)+JVS(866)*UV(34)+JVS(867)*UV(35)+JVS(868)*UV(36)+JVS(869)*UV(37)+JVS(870)*UV(38)+JVS(871)*UV(39)&
              &+JVS(872)*UV(42)+JVS(873)*UV(43)+JVS(874)*UV(48)+JVS(875)*UV(50)+JVS(877)*UV(52)+JVS(878)*UV(54)+JVS(879)&
              &*UV(55)+JVS(880)*UV(56)+JVS(881)*UV(57)+JVS(882)*UV(58)+JVS(883)*UV(59)+JVS(884)*UV(62)+JVS(885)*UV(63)&
              &+JVS(886)*UV(64)+JVS(887)*UV(65)+JVS(888)*UV(66)+JVS(889)*UV(67)+JVS(890)*UV(68)+JVS(891)*UV(69)+JVS(892)&
              &*UV(70)+JVS(893)*UV(71)+JVS(894)*UV(72)+JVS(895)*UV(73)+JVS(897)*UV(75)+JVS(898)*UV(76)+JVS(899)*UV(77)&
              &+JVS(900)*UV(78)+JVS(901)*UV(79)
  JUV(78) = JVS(902)*UV(10)+JVS(903)*UV(19)+JVS(904)*UV(21)+JVS(905)*UV(23)+JVS(906)*UV(27)+JVS(907)*UV(28)+JVS(908)&
              &*UV(29)+JVS(909)*UV(30)+JVS(910)*UV(31)+JVS(911)*UV(32)+JVS(912)*UV(34)+JVS(913)*UV(35)+JVS(914)*UV(36)&
              &+JVS(915)*UV(38)+JVS(916)*UV(39)+JVS(918)*UV(42)+JVS(919)*UV(44)+JVS(920)*UV(45)+JVS(921)*UV(48)+JVS(922)&
              &*UV(49)+JVS(923)*UV(50)+JVS(924)*UV(51)+JVS(926)*UV(54)+JVS(927)*UV(55)+JVS(928)*UV(56)+JVS(929)*UV(57)&
              &+JVS(930)*UV(58)+JVS(931)*UV(59)+JVS(932)*UV(60)+JVS(933)*UV(61)+JVS(934)*UV(62)+JVS(935)*UV(63)+JVS(936)&
              &*UV(64)+JVS(938)*UV(66)+JVS(939)*UV(67)+JVS(940)*UV(68)+JVS(941)*UV(69)+JVS(942)*UV(70)+JVS(943)*UV(71)&
              &+JVS(944)*UV(72)+JVS(945)*UV(73)+JVS(946)*UV(74)+JVS(947)*UV(75)+JVS(948)*UV(76)+JVS(949)*UV(77)+JVS(950)&
              &*UV(78)+JVS(951)*UV(79)
  JUV(79) = JVS(952)*UV(17)+JVS(953)*UV(41)+JVS(957)*UV(68)+JVS(958)*UV(69)+JVS(959)*UV(70)+JVS(960)*UV(71)+JVS(961)&
              &*UV(72)+JVS(962)*UV(73)+JVS(963)*UV(74)+JVS(964)*UV(75)+JVS(965)*UV(76)+JVS(966)*UV(77)+JVS(967)*UV(78)&
              &+JVS(968)*UV(79)
      
END SUBROUTINE Jac_SP_Vec















SUBROUTINE JacTR_SP_Vec ( JVS, UV, JTUV )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: UV(NVAR)

  REAL(kind=dp) :: JTUV(NVAR)

  JTUV(1) = JVS(1)*UV(1)
  JTUV(2) = JVS(4)*UV(2)
  JTUV(3) = JVS(8)*UV(3)
  JTUV(4) = JVS(23)*UV(4)
  JTUV(5) = JVS(33)*UV(5)
  JTUV(6) = JVS(47)*UV(6)
  JTUV(7) = JVS(50)*UV(7)
  JTUV(8) = JVS(55)*UV(8)
  JTUV(9) = JVS(60)*UV(9)
  JTUV(10) = JVS(2)*UV(1)+JVS(69)*UV(10)+JVS(637)*UV(71)+JVS(902)*UV(78)
  JTUV(11) = JVS(71)*UV(11)+JVS(443)*UV(63)+JVS(638)*UV(71)
  JTUV(12) = JVS(73)*UV(12)+JVS(639)*UV(71)+JVS(828)*UV(76)
  JTUV(13) = JVS(75)*UV(13)+JVS(368)*UV(60)+JVS(640)*UV(71)+JVS(859)*UV(77)
  JTUV(14) = JVS(77)*UV(14)+JVS(553)*UV(68)
  JTUV(15) = JVS(82)*UV(15)+JVS(554)*UV(68)+JVS(746)*UV(74)
  JTUV(16) = JVS(85)*UV(16)+JVS(696)*UV(72)+JVS(747)*UV(74)
  JTUV(17) = JVS(88)*UV(17)+JVS(748)*UV(74)+JVS(952)*UV(79)
  JTUV(18) = JVS(91)*UV(18)+JVS(617)*UV(70)+JVS(749)*UV(74)
  JTUV(19) = JVS(94)*UV(19)+JVS(641)*UV(71)+JVS(903)*UV(78)
  JTUV(20) = JVS(97)*UV(20)+JVS(224)*UV(46)+JVS(458)*UV(64)+JVS(586)*UV(69)+JVS(642)*UV(71)+JVS(860)*UV(77)
  JTUV(21) = JVS(99)*UV(21)+JVS(369)*UV(60)+JVS(392)*UV(61)+JVS(643)*UV(71)+JVS(861)*UV(77)+JVS(904)*UV(78)
  JTUV(22) = JVS(101)*UV(22)+JVS(238)*UV(47)+JVS(750)*UV(74)+JVS(788)*UV(75)
  JTUV(23) = JVS(104)*UV(23)+JVS(644)*UV(71)+JVS(721)*UV(73)+JVS(751)*UV(74)+JVS(905)*UV(78)
  JTUV(24) = JVS(107)*UV(24)+JVS(109)*UV(25)+JVS(225)*UV(46)+JVS(303)*UV(53)+JVS(370)*UV(60)+JVS(393)*UV(61)+JVS(459)&
               &*UV(64)+JVS(487)*UV(65)+JVS(587)*UV(69)+JVS(645)*UV(71)+JVS(862)*UV(77)
  JTUV(25) = JVS(110)*UV(25)+JVS(226)*UV(46)+JVS(430)*UV(62)+JVS(752)*UV(74)+JVS(829)*UV(76)
  JTUV(26) = JVS(113)*UV(26)+JVS(227)*UV(46)+JVS(304)*UV(53)+JVS(371)*UV(60)+JVS(394)*UV(61)+JVS(460)*UV(64)+JVS(488)&
               &*UV(65)+JVS(508)*UV(66)+JVS(588)*UV(69)+JVS(646)*UV(71)+JVS(863)*UV(77)
  JTUV(27) = JVS(78)*UV(14)+JVS(115)*UV(27)+JVS(143)*UV(35)+JVS(147)*UV(36)+JVS(151)*UV(37)+JVS(156)*UV(38)+JVS(170)&
               &*UV(41)+JVS(188)*UV(44)+JVS(264)*UV(49)+JVS(589)*UV(69)+JVS(647)*UV(71)+JVS(864)*UV(77)+JVS(906)*UV(78)
  JTUV(28) = JVS(117)*UV(28)+JVS(395)*UV(61)+JVS(648)*UV(71)+JVS(907)*UV(78)
  JTUV(29) = JVS(122)*UV(29)+JVS(396)*UV(61)+JVS(649)*UV(71)+JVS(830)*UV(76)+JVS(908)*UV(78)
  JTUV(30) = JVS(9)*UV(3)+JVS(126)*UV(30)+JVS(397)*UV(61)+JVS(722)*UV(73)+JVS(753)*UV(74)+JVS(909)*UV(78)
  JTUV(31) = JVS(56)*UV(8)+JVS(61)*UV(9)+JVS(130)*UV(31)+JVS(164)*UV(40)+JVS(754)*UV(74)+JVS(910)*UV(78)
  JTUV(32) = JVS(135)*UV(32)+JVS(650)*UV(71)+JVS(755)*UV(74)+JVS(789)*UV(75)+JVS(911)*UV(78)
  JTUV(33) = JVS(139)*UV(33)+JVS(200)*UV(45)+JVS(228)*UV(46)+JVS(305)*UV(53)+JVS(372)*UV(60)+JVS(398)*UV(61)+JVS(461)&
               &*UV(64)+JVS(489)*UV(65)+JVS(509)*UV(66)+JVS(556)*UV(68)+JVS(590)*UV(69)+JVS(651)*UV(71)+JVS(831)*UV(76)&
               &+JVS(865)*UV(77)
  JTUV(34) = JVS(141)*UV(34)+JVS(144)*UV(35)+JVS(148)*UV(36)+JVS(152)*UV(37)+JVS(157)*UV(38)+JVS(171)*UV(41)+JVS(182)&
               &*UV(43)+JVS(189)*UV(44)+JVS(265)*UV(49)+JVS(510)*UV(66)+JVS(591)*UV(69)+JVS(652)*UV(71)+JVS(866)*UV(77)&
               &+JVS(912)*UV(78)
  JTUV(35) = JVS(145)*UV(35)+JVS(190)*UV(44)+JVS(201)*UV(45)+JVS(266)*UV(49)+JVS(306)*UV(53)+JVS(462)*UV(64)+JVS(557)&
               &*UV(68)+JVS(653)*UV(71)+JVS(867)*UV(77)+JVS(913)*UV(78)
  JTUV(36) = JVS(149)*UV(36)+JVS(191)*UV(44)+JVS(202)*UV(45)+JVS(267)*UV(49)+JVS(307)*UV(53)+JVS(463)*UV(64)+JVS(558)&
               &*UV(68)+JVS(654)*UV(71)+JVS(868)*UV(77)+JVS(914)*UV(78)
  JTUV(37) = JVS(153)*UV(37)+JVS(192)*UV(44)+JVS(239)*UV(47)+JVS(286)*UV(51)+JVS(655)*UV(71)+JVS(790)*UV(75)+JVS(869)&
               &*UV(77)
  JTUV(38) = JVS(158)*UV(38)+JVS(203)*UV(45)+JVS(268)*UV(49)+JVS(464)*UV(64)+JVS(531)*UV(67)+JVS(656)*UV(71)+JVS(870)&
               &*UV(77)+JVS(915)*UV(78)
  JTUV(39) = JVS(5)*UV(2)+JVS(10)*UV(3)+JVS(161)*UV(39)+JVS(204)*UV(45)+JVS(269)*UV(49)+JVS(399)*UV(61)+JVS(532)*UV(67)&
               &+JVS(657)*UV(71)+JVS(871)*UV(77)+JVS(916)*UV(78)
  JTUV(40) = JVS(131)*UV(31)+JVS(165)*UV(40)+JVS(240)*UV(47)+JVS(791)*UV(75)
  JTUV(41) = JVS(62)*UV(9)+JVS(172)*UV(41)+JVS(241)*UV(47)+JVS(658)*UV(71)+JVS(792)*UV(75)+JVS(953)*UV(79)
  JTUV(42) = JVS(177)*UV(42)+JVS(465)*UV(64)+JVS(659)*UV(71)+JVS(872)*UV(77)+JVS(918)*UV(78)
  JTUV(43) = JVS(183)*UV(43)+JVS(242)*UV(47)+JVS(270)*UV(49)+JVS(287)*UV(51)+JVS(660)*UV(71)+JVS(793)*UV(75)+JVS(873)&
               &*UV(77)
  JTUV(44) = JVS(193)*UV(44)+JVS(205)*UV(45)+JVS(243)*UV(47)+JVS(559)*UV(68)+JVS(661)*UV(71)+JVS(794)*UV(75)+JVS(919)&
               &*UV(78)
  JTUV(45) = JVS(206)*UV(45)+JVS(662)*UV(71)+JVS(920)*UV(78)
  JTUV(46) = JVS(229)*UV(46)+JVS(308)*UV(53)+JVS(400)*UV(61)+JVS(560)*UV(68)+JVS(663)*UV(71)+JVS(832)*UV(76)
  JTUV(47) = JVS(244)*UV(47)+JVS(664)*UV(71)+JVS(757)*UV(74)+JVS(795)*UV(75)
  JTUV(48) = JVS(11)*UV(3)+JVS(207)*UV(45)+JVS(259)*UV(48)+JVS(271)*UV(49)+JVS(373)*UV(60)+JVS(401)*UV(61)+JVS(444)&
               &*UV(63)+JVS(466)*UV(64)+JVS(533)*UV(67)+JVS(665)*UV(71)+JVS(796)*UV(75)+JVS(833)*UV(76)+JVS(874)*UV(77)&
               &+JVS(921)*UV(78)
  JTUV(49) = JVS(208)*UV(45)+JVS(245)*UV(47)+JVS(272)*UV(49)+JVS(402)*UV(61)+JVS(666)*UV(71)+JVS(697)*UV(72)+JVS(797)&
               &*UV(75)+JVS(922)*UV(78)
  JTUV(50) = JVS(6)*UV(2)+JVS(12)*UV(3)+JVS(24)*UV(4)+JVS(57)*UV(8)+JVS(63)*UV(9)+JVS(209)*UV(45)+JVS(281)*UV(50)&
               &+JVS(374)*UV(60)+JVS(403)*UV(61)+JVS(445)*UV(63)+JVS(467)*UV(64)+JVS(490)*UV(65)+JVS(534)*UV(67)+JVS(592)&
               &*UV(69)+JVS(667)*UV(71)+JVS(798)*UV(75)+JVS(834)*UV(76)+JVS(875)*UV(77)+JVS(923)*UV(78)
  JTUV(51) = JVS(166)*UV(40)+JVS(184)*UV(43)+JVS(288)*UV(51)+JVS(759)*UV(74)+JVS(924)*UV(78)
  JTUV(52) = JVS(13)*UV(3)+JVS(34)*UV(5)+JVS(210)*UV(45)+JVS(298)*UV(52)+JVS(309)*UV(53)+JVS(337)*UV(55)+JVS(349)*UV(57)&
               &+JVS(361)*UV(59)+JVS(405)*UV(61)+JVS(446)*UV(63)+JVS(511)*UV(66)+JVS(535)*UV(67)+JVS(593)*UV(69)+JVS(618)&
               &*UV(70)+JVS(669)*UV(71)+JVS(760)*UV(74)+JVS(800)*UV(75)+JVS(835)*UV(76)+JVS(877)*UV(77)
  JTUV(53) = JVS(310)*UV(53)+JVS(723)*UV(73)+JVS(761)*UV(74)+JVS(801)*UV(75)
  JTUV(54) = JVS(14)*UV(3)+JVS(35)*UV(5)+JVS(79)*UV(14)+JVS(211)*UV(45)+JVS(230)*UV(46)+JVS(274)*UV(49)+JVS(311)*UV(53)&
               &+JVS(332)*UV(54)+JVS(406)*UV(61)+JVS(431)*UV(62)+JVS(447)*UV(63)+JVS(468)*UV(64)+JVS(512)*UV(66)+JVS(536)&
               &*UV(67)+JVS(561)*UV(68)+JVS(594)*UV(69)+JVS(670)*UV(71)+JVS(699)*UV(72)+JVS(762)*UV(74)+JVS(802)*UV(75)&
               &+JVS(878)*UV(77)+JVS(926)*UV(78)
  JTUV(55) = JVS(15)*UV(3)+JVS(194)*UV(44)+JVS(212)*UV(45)+JVS(247)*UV(47)+JVS(338)*UV(55)+JVS(407)*UV(61)+JVS(448)&
               &*UV(63)+JVS(469)*UV(64)+JVS(491)*UV(65)+JVS(537)*UV(67)+JVS(562)*UV(68)+JVS(619)*UV(70)+JVS(671)*UV(71)&
               &+JVS(700)*UV(72)+JVS(803)*UV(75)+JVS(879)*UV(77)+JVS(927)*UV(78)
  JTUV(56) = JVS(16)*UV(3)+JVS(25)*UV(4)+JVS(36)*UV(5)+JVS(213)*UV(45)+JVS(231)*UV(46)+JVS(312)*UV(53)+JVS(344)*UV(56)&
               &+JVS(375)*UV(60)+JVS(408)*UV(61)+JVS(432)*UV(62)+JVS(449)*UV(63)+JVS(470)*UV(64)+JVS(492)*UV(65)+JVS(513)&
               &*UV(66)+JVS(538)*UV(67)+JVS(595)*UV(69)+JVS(672)*UV(71)+JVS(804)*UV(75)+JVS(837)*UV(76)+JVS(880)*UV(77)&
               &+JVS(928)*UV(78)
  JTUV(57) = JVS(17)*UV(3)+JVS(37)*UV(5)+JVS(195)*UV(44)+JVS(214)*UV(45)+JVS(248)*UV(47)+JVS(275)*UV(49)+JVS(350)*UV(57)&
               &+JVS(376)*UV(60)+JVS(409)*UV(61)+JVS(433)*UV(62)+JVS(471)*UV(64)+JVS(493)*UV(65)+JVS(514)*UV(66)+JVS(539)&
               &*UV(67)+JVS(564)*UV(68)+JVS(596)*UV(69)+JVS(620)*UV(70)+JVS(673)*UV(71)+JVS(701)*UV(72)+JVS(805)*UV(75)&
               &+JVS(881)*UV(77)+JVS(929)*UV(78)
  JTUV(58) = JVS(18)*UV(3)+JVS(26)*UV(4)+JVS(38)*UV(5)+JVS(173)*UV(41)+JVS(215)*UV(45)+JVS(232)*UV(46)+JVS(313)*UV(53)&
               &+JVS(339)*UV(55)+JVS(351)*UV(57)+JVS(356)*UV(58)+JVS(362)*UV(59)+JVS(377)*UV(60)+JVS(410)*UV(61)+JVS(434)&
               &*UV(62)+JVS(450)*UV(63)+JVS(472)*UV(64)+JVS(494)*UV(65)+JVS(515)*UV(66)+JVS(540)*UV(67)+JVS(565)*UV(68)&
               &+JVS(597)*UV(69)+JVS(674)*UV(71)+JVS(702)*UV(72)+JVS(766)*UV(74)+JVS(806)*UV(75)+JVS(838)*UV(76)+JVS(882)&
               &*UV(77)+JVS(930)*UV(78)
  JTUV(59) = JVS(19)*UV(3)+JVS(196)*UV(44)+JVS(216)*UV(45)+JVS(314)*UV(53)+JVS(363)*UV(59)+JVS(411)*UV(61)+JVS(451)&
               &*UV(63)+JVS(473)*UV(64)+JVS(495)*UV(65)+JVS(516)*UV(66)+JVS(541)*UV(67)+JVS(566)*UV(68)+JVS(598)*UV(69)&
               &+JVS(622)*UV(70)+JVS(675)*UV(71)+JVS(703)*UV(72)+JVS(839)*UV(76)+JVS(883)*UV(77)+JVS(931)*UV(78)
  JTUV(60) = JVS(217)*UV(45)+JVS(251)*UV(47)+JVS(378)*UV(60)+JVS(567)*UV(68)+JVS(676)*UV(71)+JVS(808)*UV(75)+JVS(840)&
               &*UV(76)+JVS(932)*UV(78)
  JTUV(61) = JVS(127)*UV(30)+JVS(218)*UV(45)+JVS(252)*UV(47)+JVS(412)*UV(61)+JVS(677)*UV(71)+JVS(809)*UV(75)+JVS(933)&
               &*UV(78)
  JTUV(62) = JVS(233)*UV(46)+JVS(315)*UV(53)+JVS(379)*UV(60)+JVS(413)*UV(61)+JVS(435)*UV(62)+JVS(474)*UV(64)+JVS(496)&
               &*UV(65)+JVS(517)*UV(66)+JVS(599)*UV(69)+JVS(678)*UV(71)+JVS(770)*UV(74)+JVS(884)*UV(77)+JVS(934)*UV(78)
  JTUV(63) = JVS(64)*UV(9)+JVS(219)*UV(45)+JVS(260)*UV(48)+JVS(276)*UV(49)+JVS(282)*UV(50)+JVS(299)*UV(52)+JVS(316)&
               &*UV(53)+JVS(333)*UV(54)+JVS(340)*UV(55)+JVS(345)*UV(56)+JVS(357)*UV(58)+JVS(364)*UV(59)+JVS(380)*UV(60)&
               &+JVS(414)*UV(61)+JVS(452)*UV(63)+JVS(475)*UV(64)+JVS(497)*UV(65)+JVS(518)*UV(66)+JVS(542)*UV(67)+JVS(600)&
               &*UV(69)+JVS(623)*UV(70)+JVS(730)*UV(73)+JVS(771)*UV(74)+JVS(811)*UV(75)+JVS(842)*UV(76)+JVS(885)*UV(77)&
               &+JVS(935)*UV(78)
  JTUV(64) = JVS(220)*UV(45)+JVS(253)*UV(47)+JVS(381)*UV(60)+JVS(476)*UV(64)+JVS(601)*UV(69)+JVS(680)*UV(71)+JVS(705)&
               &*UV(72)+JVS(812)*UV(75)+JVS(886)*UV(77)+JVS(936)*UV(78)
  JTUV(65) = JVS(317)*UV(53)+JVS(382)*UV(60)+JVS(415)*UV(61)+JVS(477)*UV(64)+JVS(498)*UV(65)+JVS(571)*UV(68)+JVS(602)&
               &*UV(69)+JVS(681)*UV(71)+JVS(706)*UV(72)+JVS(887)*UV(77)
  JTUV(66) = JVS(318)*UV(53)+JVS(383)*UV(60)+JVS(416)*UV(61)+JVS(478)*UV(64)+JVS(499)*UV(65)+JVS(519)*UV(66)+JVS(572)&
               &*UV(68)+JVS(603)*UV(69)+JVS(682)*UV(71)+JVS(707)*UV(72)+JVS(888)*UV(77)+JVS(938)*UV(78)
  JTUV(67) = JVS(7)*UV(2)+JVS(20)*UV(3)+JVS(27)*UV(4)+JVS(39)*UV(5)+JVS(65)*UV(9)+JVS(72)*UV(11)+JVS(80)*UV(14)+JVS(159)&
               &*UV(38)+JVS(162)*UV(39)+JVS(174)*UV(41)+JVS(197)*UV(44)+JVS(221)*UV(45)+JVS(234)*UV(46)+JVS(261)*UV(48)&
               &+JVS(277)*UV(49)+JVS(283)*UV(50)+JVS(300)*UV(52)+JVS(319)*UV(53)+JVS(334)*UV(54)+JVS(341)*UV(55)+JVS(346)&
               &*UV(56)+JVS(353)*UV(57)+JVS(358)*UV(58)+JVS(365)*UV(59)+JVS(384)*UV(60)+JVS(417)*UV(61)+JVS(453)*UV(63)&
               &+JVS(479)*UV(64)+JVS(500)*UV(65)+JVS(520)*UV(66)+JVS(543)*UV(67)+JVS(573)*UV(68)+JVS(604)*UV(69)+JVS(624)&
               &*UV(70)+JVS(683)*UV(71)+JVS(708)*UV(72)+JVS(733)*UV(73)+JVS(775)*UV(74)+JVS(815)*UV(75)+JVS(846)*UV(76)&
               &+JVS(889)*UV(77)+JVS(939)*UV(78)
  JTUV(68) = JVS(28)*UV(4)+JVS(48)*UV(6)+JVS(83)*UV(15)+JVS(289)*UV(51)+JVS(320)*UV(53)+JVS(385)*UV(60)+JVS(418)*UV(61)&
               &+JVS(521)*UV(66)+JVS(544)*UV(67)+JVS(574)*UV(68)+JVS(605)*UV(69)+JVS(625)*UV(70)+JVS(709)*UV(72)+JVS(734)&
               &*UV(73)+JVS(776)*UV(74)+JVS(816)*UV(75)+JVS(847)*UV(76)+JVS(890)*UV(77)+JVS(940)*UV(78)+JVS(957)*UV(79)
  JTUV(69) = JVS(29)*UV(4)+JVS(40)*UV(5)+JVS(118)*UV(28)+JVS(178)*UV(42)+JVS(321)*UV(53)+JVS(419)*UV(61)+JVS(438)*UV(62)&
               &+JVS(501)*UV(65)+JVS(522)*UV(66)+JVS(575)*UV(68)+JVS(606)*UV(69)+JVS(626)*UV(70)+JVS(710)*UV(72)+JVS(735)&
               &*UV(73)+JVS(777)*UV(74)+JVS(817)*UV(75)+JVS(848)*UV(76)+JVS(891)*UV(77)+JVS(941)*UV(78)+JVS(958)*UV(79)
  JTUV(70) = JVS(41)*UV(5)+JVS(51)*UV(7)+JVS(92)*UV(18)+JVS(290)*UV(51)+JVS(322)*UV(53)+JVS(386)*UV(60)+JVS(420)*UV(61)&
               &+JVS(545)*UV(67)+JVS(576)*UV(68)+JVS(607)*UV(69)+JVS(627)*UV(70)+JVS(711)*UV(72)+JVS(736)*UV(73)+JVS(778)&
               &*UV(74)+JVS(818)*UV(75)+JVS(849)*UV(76)+JVS(892)*UV(77)+JVS(942)*UV(78)+JVS(959)*UV(79)
  JTUV(71) = JVS(3)*UV(1)+JVS(21)*UV(3)+JVS(66)*UV(9)+JVS(70)*UV(10)+JVS(74)*UV(12)+JVS(76)*UV(13)+JVS(81)*UV(14)&
               &+JVS(95)*UV(19)+JVS(98)*UV(20)+JVS(100)*UV(21)+JVS(105)*UV(23)+JVS(108)*UV(24)+JVS(111)*UV(25)+JVS(114)&
               &*UV(26)+JVS(116)*UV(27)+JVS(119)*UV(28)+JVS(123)*UV(29)+JVS(136)*UV(32)+JVS(140)*UV(33)+JVS(142)*UV(34)&
               &+JVS(146)*UV(35)+JVS(150)*UV(36)+JVS(154)*UV(37)+JVS(160)*UV(38)+JVS(163)*UV(39)+JVS(175)*UV(41)+JVS(179)&
               &*UV(42)+JVS(185)*UV(43)+JVS(198)*UV(44)+JVS(222)*UV(45)+JVS(235)*UV(46)+JVS(255)*UV(47)+JVS(262)*UV(48)&
               &+JVS(278)*UV(49)+JVS(284)*UV(50)+JVS(291)*UV(51)+JVS(301)*UV(52)+JVS(323)*UV(53)+JVS(335)*UV(54)+JVS(342)&
               &*UV(55)+JVS(347)*UV(56)+JVS(354)*UV(57)+JVS(359)*UV(58)+JVS(366)*UV(59)+JVS(387)*UV(60)+JVS(421)*UV(61)&
               &+JVS(439)*UV(62)+JVS(481)*UV(64)+JVS(502)*UV(65)+JVS(523)*UV(66)+JVS(546)*UV(67)+JVS(577)*UV(68)+JVS(608)&
               &*UV(69)+JVS(628)*UV(70)+JVS(687)*UV(71)+JVS(712)*UV(72)+JVS(737)*UV(73)+JVS(779)*UV(74)+JVS(819)*UV(75)&
               &+JVS(850)*UV(76)+JVS(893)*UV(77)+JVS(943)*UV(78)+JVS(960)*UV(79)
  JTUV(72) = JVS(42)*UV(5)+JVS(52)*UV(7)+JVS(86)*UV(16)+JVS(292)*UV(51)+JVS(324)*UV(53)+JVS(388)*UV(60)+JVS(422)*UV(61)&
               &+JVS(524)*UV(66)+JVS(547)*UV(67)+JVS(578)*UV(68)+JVS(609)*UV(69)+JVS(629)*UV(70)+JVS(713)*UV(72)+JVS(738)&
               &*UV(73)+JVS(780)*UV(74)+JVS(820)*UV(75)+JVS(851)*UV(76)+JVS(894)*UV(77)+JVS(944)*UV(78)+JVS(961)*UV(79)
  JTUV(73) = JVS(22)*UV(3)+JVS(106)*UV(23)+JVS(128)*UV(30)+JVS(293)*UV(51)+JVS(325)*UV(53)+JVS(389)*UV(60)+JVS(423)&
               &*UV(61)+JVS(440)*UV(62)+JVS(455)*UV(63)+JVS(548)*UV(67)+JVS(579)*UV(68)+JVS(610)*UV(69)+JVS(630)*UV(70)&
               &+JVS(689)*UV(71)+JVS(714)*UV(72)+JVS(739)*UV(73)+JVS(781)*UV(74)+JVS(821)*UV(75)+JVS(852)*UV(76)+JVS(895)&
               &*UV(77)+JVS(945)*UV(78)+JVS(962)*UV(79)
  JTUV(74) = JVS(58)*UV(8)+JVS(67)*UV(9)+JVS(84)*UV(15)+JVS(87)*UV(16)+JVS(89)*UV(17)+JVS(93)*UV(18)+JVS(102)*UV(22)&
               &+JVS(112)*UV(25)+JVS(132)*UV(31)+JVS(137)*UV(32)+JVS(167)*UV(40)+JVS(256)*UV(47)+JVS(294)*UV(51)+JVS(441)&
               &*UV(62)+JVS(456)*UV(63)+JVS(549)*UV(67)+JVS(580)*UV(68)+JVS(631)*UV(70)+JVS(690)*UV(71)+JVS(715)*UV(72)&
               &+JVS(740)*UV(73)+JVS(782)*UV(74)+JVS(822)*UV(75)+JVS(946)*UV(78)+JVS(963)*UV(79)
  JTUV(75) = JVS(59)*UV(8)+JVS(68)*UV(9)+JVS(103)*UV(22)+JVS(133)*UV(31)+JVS(155)*UV(37)+JVS(168)*UV(40)+JVS(176)*UV(41)&
               &+JVS(186)*UV(43)+JVS(199)*UV(44)+JVS(223)*UV(45)+JVS(237)*UV(46)+JVS(257)*UV(47)+JVS(263)*UV(48)+JVS(279)&
               &*UV(49)+JVS(285)*UV(50)+JVS(295)*UV(51)+JVS(302)*UV(52)+JVS(327)*UV(53)+JVS(336)*UV(54)+JVS(343)*UV(55)&
               &+JVS(348)*UV(56)+JVS(355)*UV(57)+JVS(360)*UV(58)+JVS(367)*UV(59)+JVS(390)*UV(60)+JVS(425)*UV(61)+JVS(442)&
               &*UV(62)+JVS(457)*UV(63)+JVS(484)*UV(64)+JVS(505)*UV(65)+JVS(581)*UV(68)+JVS(612)*UV(69)+JVS(632)*UV(70)&
               &+JVS(691)*UV(71)+JVS(716)*UV(72)+JVS(741)*UV(73)+JVS(783)*UV(74)+JVS(823)*UV(75)+JVS(854)*UV(76)+JVS(897)&
               &*UV(77)+JVS(947)*UV(78)+JVS(964)*UV(79)
  JTUV(76) = JVS(30)*UV(4)+JVS(43)*UV(5)+JVS(120)*UV(28)+JVS(124)*UV(29)+JVS(328)*UV(53)+JVS(426)*UV(61)+JVS(506)*UV(65)&
               &+JVS(528)*UV(66)+JVS(582)*UV(68)+JVS(613)*UV(69)+JVS(633)*UV(70)+JVS(717)*UV(72)+JVS(742)*UV(73)+JVS(784)&
               &*UV(74)+JVS(824)*UV(75)+JVS(855)*UV(76)+JVS(898)*UV(77)+JVS(948)*UV(78)+JVS(965)*UV(79)
  JTUV(77) = JVS(31)*UV(4)+JVS(44)*UV(5)+JVS(121)*UV(28)+JVS(180)*UV(42)+JVS(329)*UV(53)+JVS(427)*UV(61)+JVS(507)*UV(65)&
               &+JVS(529)*UV(66)+JVS(583)*UV(68)+JVS(614)*UV(69)+JVS(634)*UV(70)+JVS(718)*UV(72)+JVS(743)*UV(73)+JVS(785)&
               &*UV(74)+JVS(825)*UV(75)+JVS(856)*UV(76)+JVS(899)*UV(77)+JVS(949)*UV(78)+JVS(966)*UV(79)
  JTUV(78) = JVS(32)*UV(4)+JVS(45)*UV(5)+JVS(49)*UV(6)+JVS(53)*UV(7)+JVS(96)*UV(19)+JVS(125)*UV(29)+JVS(129)*UV(30)&
               &+JVS(134)*UV(31)+JVS(138)*UV(32)+JVS(169)*UV(40)+JVS(181)*UV(42)+JVS(187)*UV(43)+JVS(258)*UV(47)+JVS(296)&
               &*UV(51)+JVS(330)*UV(53)+JVS(428)*UV(61)+JVS(551)*UV(67)+JVS(584)*UV(68)+JVS(615)*UV(69)+JVS(635)*UV(70)&
               &+JVS(694)*UV(71)+JVS(719)*UV(72)+JVS(744)*UV(73)+JVS(786)*UV(74)+JVS(826)*UV(75)+JVS(857)*UV(76)+JVS(900)&
               &*UV(77)+JVS(950)*UV(78)+JVS(967)*UV(79)
  JTUV(79) = JVS(46)*UV(5)+JVS(54)*UV(7)+JVS(90)*UV(17)+JVS(297)*UV(51)+JVS(331)*UV(53)+JVS(391)*UV(60)+JVS(429)*UV(61)&
               &+JVS(530)*UV(66)+JVS(552)*UV(67)+JVS(585)*UV(68)+JVS(616)*UV(69)+JVS(636)*UV(70)+JVS(720)*UV(72)+JVS(745)&
               &*UV(73)+JVS(787)*UV(74)+JVS(827)*UV(75)+JVS(858)*UV(76)+JVS(901)*UV(77)+JVS(951)*UV(78)+JVS(968)*UV(79)
      
END SUBROUTINE JacTR_SP_Vec






END MODULE saprc99_Jacobian

