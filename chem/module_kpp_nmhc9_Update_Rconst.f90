























MODULE nmhc9_UpdateRconstWRF

  USE nmhc9_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE nmhc9_Update_RCONST(  &


PRESS, &


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


REAL(KIND=dp), INTENT(IN)  :: PRESS










  RCONST(1) = (j(Pj_o31d))
  RCONST(2) = (.7902_dp*2.1D-11*exp(115._dp/TEMP)+.20946_dp*3.2D-11*exp(70._dp/TEMP))
  RCONST(3) = (2.2D-10)
  RCONST(4) = (min(1.D-11,.20946_dp*j(Pj_O2)))
  RCONST(5) = (1.7D-12*exp(-940._dp/TEMP))
  RCONST(6) = (1.0D-14*exp(-490._dp/TEMP))
  RCONST(7) = (4.8D-11*exp(250._dp/TEMP))
  RCONST(8) = (RHO2HO2(C_M,C_H2O,TEMP))
  RCONST(9) = (j(Pj_h2o2))
  RCONST(10) = (2.9D-12*exp(-160._dp/TEMP))
  RCONST(11) = (1.57D-13+3.54D-33*C_M)
  RCONST(12) = (1.85D-20*exp(2.82_dp*log(TEMP)-987._dp/TEMP))
  RCONST(13) = (1.5D-10)
  RCONST(14) = (4.1D-13*exp(750._dp/TEMP)/(1._dp+1._dp/497.7_dp*EXP(1160._dp/TEMP)))
  RCONST(15) = (4.1D-13*exp(750._dp/TEMP)/(1._dp+497.7_dp*EXP(-1160._dp/TEMP)))
  RCONST(16) = (2.8D-12*exp(300._dp/TEMP))
  RCONST(17) = (9.5D-14*exp(390._dp/TEMP)/(1._dp+1._dp/26.2_dp*EXP(1130._dp/TEMP)))
  RCONST(18) = (9.5D-14*exp(390._dp/TEMP)/(1._dp+26.2_dp*EXP(-1130._dp/TEMP)))
  RCONST(19) = (1.3D-12)
  RCONST(20) = (j(Pj_ch3o2h))
  RCONST(21) = (RCH3OOHOH(TEMP))
  RCONST(22) = (j(Pj_ch2or))
  RCONST(23) = (j(Pj_ch2om))
  RCONST(24) = (9.52D-18*exp(2.03_dp*log(TEMP)+636._dp/TEMP))
  RCONST(25) = (3.4D-13*exp(-1900._dp/TEMP))
  RCONST(26) = (3.0D-12*exp(-1500._dp/TEMP))
  RCONST(27) = (3.5D-12*exp(250._dp/TEMP))
  RCONST(28) = (j(Pj_no2))
  RCONST(29) = (1.2D-13*exp(-2450._dp/TEMP))
  RCONST(30) = (TROE2(C_M,TEMP,.933_dp,2.85D-30,-2.67_dp,3.13D-11,363._dp))
  RCONST(31) = (RJPL(1.8D-31,3.2_dp,4.7D-12,1.4_dp,C_M,TEMP))
  RCONST(32) = (j(Pj_hno3))
  RCONST(33) = (RHNO3(C_M,TEMP))
  RCONST(34) = (j(Pj_no3o))
  RCONST(35) = (j(Pj_no3o2))
  RCONST(36) = (1.5D-11*exp(170._dp/TEMP))
  RCONST(37) = (RJPL(2.D-30,4.4_dp,1.4D-12,.7_dp,C_M,TEMP))
  RCONST(38) = (3.5D-12)
  RCONST(39) = (j(Pj_N2O5))
  RCONST(40) = (RJPL(2.D-30,4.4_dp,1.4D-12,.7_dp,C_M,TEMP)/(3.D-27*exp(10990._dp/TEMP)))
  RCONST(41) = (3.D-7*TEMP)
  RCONST(42) = (2.5D-22+C_H2O*1.8D-39)
  RCONST(43) = (j(Pj_HNO4_2))
  RCONST(44) = (RJPL(1.8D-31,3.2_dp,4.7D-12,1.4_dp,C_M,TEMP)/(2.1D-27*exp(10900._dp/TEMP)))
  RCONST(45) = (1.3D-12*exp(380._dp/TEMP))
  RCONST(46) = (5.31D-7*5.5D-12*exp(-2000._dp/TEMP))
  RCONST(47) = (7.3D-12*exp(-620._dp/TEMP))
  RCONST(48) = (RJPL(1.3D-30,4.0_dp,7.5D-12,2.0_dp,C_M,TEMP))
  RCONST(49) = (RJPL(1.3D-30,4.0_dp,7.5D-12,2.0_dp,C_M,TEMP)/(1.3D-28*exp(11200._dp/TEMP)))
  RCONST(50) = (j(Pj_HNO4_2))
  RCONST(51) = (2.54D-11*exp(410._dp/TEMP))
  RCONST(52) = (3.03D-12*exp(-446._dp/TEMP))
  RCONST(53) = (7.86D-15*exp(-1913._dp/TEMP))
  RCONST(54) = (2.22D-13*exp(1300._dp/TEMP))
  RCONST(55) = (2.54D-12*exp(360._dp/TEMP))
  RCONST(56) = (j(Pj_ch3o2h))
  RCONST(57) = (1.D-10)
  RCONST(58) = (.5_dp*(4.1D-12*exp(452._dp/TEMP)+1.9D-11*exp(175._dp/TEMP)))
  RCONST(59) = (.019_dp*j(Pj_ch2om)+.015_dp*j(Pj_MGLO))
  RCONST(60) = (.5_dp*(1.36D-15*exp(-2112._dp/TEMP)+7.51D-16*exp(-1521._dp/TEMP)))
  RCONST(61) = (2.54D-12*exp(360._dp/TEMP))
  RCONST(62) = (1.82D-13*exp(1300._dp/TEMP))
  RCONST(63) = (j(Pj_ch3o2h))
  RCONST(64) = (3.D-11)
  RCONST(65) = (1.3D-11)
  RCONST(66) = (3.7_dp*j(Pj_PAN))
  RCONST(67) = (.25_dp*RJPL(9.7D-29,5.6_dp,9.3D-12,1.5_dp,C_M,TEMP))
  RCONST(68) = (2.D-12)
  RCONST(69) = (2.D-12)
  RCONST(70) = (2.D-12)
  RCONST(71) = (2.D-12)
  RCONST(72) = (8.4D-13*exp(830._dp/TEMP))
  RCONST(73) = (j(Pj_MGLO))
  RCONST(74) = (3.D-12)
  RCONST(75) = (.074_dp*j(Pj_ch2or))
  RCONST(76) = (4.3D-13*exp(1040._dp/TEMP)/(1._dp+1._dp/37._dp*exp(660._dp/TEMP)))
  RCONST(77) = (4.3D-13*exp(1040._dp/TEMP)/(1._dp+37._dp*exp(-660._dp/TEMP)))
  RCONST(78) = (8.1D-12*exp(270._dp/TEMP))
  RCONST(79) = (RJPL(9.7D-29,5.6_dp,9.3D-12,1.5_dp,C_M,TEMP))
  RCONST(80) = (2.0D-12*exp(500._dp/TEMP)/(1._dp+1._dp/2.2D6*exp(3820._dp/TEMP)))
  RCONST(81) = (2.0D-12*exp(500._dp/TEMP)/(1._dp+2.2D6*exp(-3820._dp/TEMP)))
  RCONST(82) = (2.5D-12*exp(500._dp/TEMP))
  RCONST(83) = (4.D-12)
  RCONST(84) = (.025_dp*j(Pj_ch2or))
  RCONST(85) = (3.8D-12*exp(200._dp/TEMP))
  RCONST(86) = (2.D-14)
  RCONST(87) = (j(Pj_PAN))
  RCONST(88) = (RJPL(9.7D-29,5.6_dp,9.3D-12,1.5_dp,C_M,TEMP)/(9.D-29*exp(14000._dp/TEMP)))
  RCONST(89) = (3.2D-11)
  RCONST(90) = (RJPL(9.7D-29,5.6_dp,9.3D-12,1.5_dp,C_M,TEMP)/(9.D-29*exp(14000._dp/TEMP)))
  RCONST(91) = (j(Pj_PAN))
  RCONST(92) = (4.D-13*exp(200._dp/TEMP))
  RCONST(93) = (4.D-13)
  RCONST(94) = (5.6D-12*exp(270._dp/TEMP))
  RCONST(95) = (.19_dp*j(Pj_ch2or))
  RCONST(96) = (1.49D-17*TEMP*TEMP*exp(-499._dp/TEMP))
  RCONST(97) = (7.5D-13*exp(700._dp/TEMP))
  RCONST(98) = (1.6D-13*exp(195._dp/TEMP))
  RCONST(99) = (4.9D-12*exp(211._dp/TEMP))
  RCONST(100) = (2.7D-12*exp(350._dp/TEMP))
  RCONST(101) = (2.3D-12)
  RCONST(102) = (5.6D-12*exp(270._dp/TEMP))
  RCONST(103) = (.19_dp*j(Pj_ch2or))
  RCONST(104) = (1.4D-12*exp(-1900._dp/TEMP))
  RCONST(105) = (j(Pj_ch3o2h))
  RCONST(106) = (RCH3OOHOH(TEMP))
  RCONST(107) = (1.65D-17*TEMP*TEMP*exp(-87._dp/TEMP))
  RCONST(108) = (2.7D-12*exp(360._dp/TEMP))
  RCONST(109) = (1.9D-13*exp(1300._dp/TEMP))
  RCONST(110) = (2.0D-14*exp(-886._dp/TEMP))
  RCONST(111) = (j(Pj_ch3o2h))
  RCONST(112) = (RCH3OOHOH(TEMP))
  RCONST(113) = (3.7_dp*j(Pj_PAN))
  RCONST(114) = (6.2D-13*exp(-230._dp/TEMP))
  RCONST(115) = (1.33D-13+3.82D-11*exp(-2000._dp/TEMP))
  RCONST(116) = (j(Pj_ACET))
  RCONST(117) = (2.9D-12*exp(300._dp/TEMP))
  RCONST(118) = (8.6D-13*exp(700._dp/TEMP))
  RCONST(119) = (7.5D-13*exp(500._dp/TEMP))
  RCONST(120) = (RCH3OOHOH(TEMP))
  RCONST(121) = (j(Pj_ch3o2h))
  RCONST(122) = (RALKE(8.D-27,3.5_dp,3.D-11,0.5_dp,C_M,TEMP))
  RCONST(123) = (6.5D-15*exp(-1900._dp/TEMP))
  RCONST(124) = (4.6D-13*exp(-1155._dp/TEMP))
  RCONST(125) = (4.2D-12*exp(180._dp/TEMP))
  RCONST(126) = (6.5D-13*exp(650._dp/TEMP))
  RCONST(127) = (3.8D-12*exp(200._dp/TEMP))
  RCONST(128) = (j(Pj_ch3o2h))
  RCONST(129) = (RJPL(1.D-28,0.8_dp,8.8D-12,0._dp,C_H2O,TEMP))
  RCONST(130) = (1.2D-14*exp(-2630._dp/TEMP))
  RCONST(131) = (1.81D-17*TEMP*TEMP*exp(114._dp/TEMP))
  RCONST(132) = (2.7D-12*exp(360._dp/TEMP))
  RCONST(133) = (1.9D-13*exp(1300._dp/TEMP))
  RCONST(134) = (9.46D-14*exp(431._dp/TEMP))
  RCONST(135) = (RCH3OOHOH(TEMP))
  RCONST(136) = (j(Pj_ch3o2h))
  RCONST(137) = (.42_dp*j(Pj_ch2or))
  RCONST(138) = (1.3D-12*exp(-25._dp/TEMP))
  RCONST(139) = (2.7D-12*exp(360._dp/TEMP))
  RCONST(140) = (1.9D-13*exp(1300._dp/TEMP))
  RCONST(141) = (RCH3OOHOH(TEMP))
  RCONST(142) = (j(Pj_ch3o2h))
  RCONST(143) = (2.15_dp*j(Pj_MGLO))
  RCONST(144) = (1.7D-12)
  RCONST(145) = (3.7_dp*j(Pj_PAN))
END SUBROUTINE nmhc9_Update_RCONST













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














      REAL(KIND=dp) FUNCTION RJPL( K0300, Q, KU300, R, M, T )
      REAL(KIND=dp) :: k0300,q,ku300,r,m,t
      REAL(KIND=dp) :: tt,k0,ku,k0m,kk,lgkk,e,f

      TT= T / 3.D2
      K0= K0300 * exp(-1._dp*Q*log(TT))
      KU= KU300 * exp(-1._dp*R*log(TT))
      K0M= K0 * M
      KK= K0M / KU
      LGKK=0.43429448190324926_dp * LOG(KK) 
      E=1.D0 / ( 1.D0 + LGKK*LGKK )
      F=exp(-0.5108256237659887_dp*E)       
      RJPL = F * K0M / ( 1.D0 + KK )
      END FUNCTION

      REAL(KIND=dp) FUNCTION RALKE( K0300, Q, KU, Fc, M, T )
      REAL(KIND=dp) :: k0300,q,m,t,Fc
      real(KIND=dp) :: tt,k0,ku,k0m,kk,lgkk,e,f

      TT= T / 3.D2
      K0= K0300 * exp(-1._dp*Q*log(TT))
      K0M= K0 * M
      KK= K0M / KU
      LGKK=0.43429448190324926_dp * LOG(KK) 
      E=1.D0 / ( 1.D0 + LGKK*LGKK )
      F=exp(log(Fc)*E)
      RALKE = F * K0M / ( 1.D0 + KK )
      END FUNCTION




      REAL(KIND=dp) FUNCTION TROE2(M,T,beta,k0,k0e,kinf,Tc)
      REAL(KIND=dp) :: M,T,beta,k0,k0e,kinf,Tc
      REAL(KIND=dp) :: k0t,bcrit,Trat,dN,N,Bx,F

      k0t = k0 * exp(k0e*log(T/3.D2))
      bcrit = beta*M*k0t/kinf
      Trat = T/Tc
      dN=sign(0.1_dp-0.2605766891419492_dp*Trat,1.-bcrit) 
      N = 0.75_dp + 0.5515539920171264_dp*Trat           
      Bx = (0.43429448190324926_dp*log(bcrit)-0.12_dp) / (N+dN)
      F = exp(-1._dp *Trat/(1._dp + Bx*Bx))
      TROE2 = k0t * (beta*M/(1.+bcrit)) * F
      END FUNCTION

      REAL(KIND=dp) FUNCTION RHNO3(M,T)
      REAL(KIND=dp) :: M,T
      REAL(KIND=dp) :: K0,K2,K3


      K0=2.4D-14*EXP(460._dp/T)
      K2=2.7D-17*EXP(2199._dp/T)
      K3=M*6.5D-34*EXP(1335._dp/T)
      RHNO3 = K0 + K2 / ( 1._dp + K2/K3 )
      END FUNCTION

      REAL(KIND=dp) FUNCTION RHO2HO2(M,H2O,T)
      REAL(KIND=dp) :: M,H2O,T
      REAL(KIND=dp) :: RX1,RX2,RX3


      RX1= 1.5D-12 *EXP(19._dp/T)         
      RX2= 1.7D-33 *EXP(1000._dp/T) * M
      RX3= 1.4D-21 *EXP(2200._dp/T) * H2O
      RHO2HO2 = (RX1 + RX2)*(1._dp + RX3)
      END FUNCTION

      REAL(KIND=dp)  FUNCTION RCH3OOHOH(T) 
      REAL(KIND=dp) :: T
      RCH3OOHOH = 3.8D-12*exp(200._dp/T)
      END FUNCTION





END MODULE nmhc9_UpdateRconstWRF

