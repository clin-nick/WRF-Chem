























MODULE cbm4_UpdateRconstWRF

  USE cbm4_Parameters
  IMPLICIT NONE

CONTAINS


SUBROUTINE cbm4_Update_RCONST(  &



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










  RCONST(1) = (j(pj_no2))
  RCONST(2) = (j(pj_o33p))
  RCONST(3) = (j(pj_o31d))
  RCONST(4) = (j(pj_no3o)+j(pj_no3o2))
  RCONST(5) = (j(pj_hno2))
  RCONST(6) = (j(pj_h2o2))
  RCONST(7) = (j(pj_ch2or))
  RCONST(8) = (j(pj_ch2om))
  RCONST(9) = (4.6D-4*j(pj_no2))
  RCONST(10) = (9.04_dp*j(pj_ch2or))
  RCONST(11) = (9.64_dp*j(pj_ch2or))
  RCONST(12) = (ARR2(1.4D+3,1175.0_dp,TEMP))
  RCONST(13) = (ARR2(1.8D-12,-1370.0_dp,TEMP))
  RCONST(14) = (9.3D-12)
  RCONST(15) = (ARR2(1.6D-13,687.0_dp,TEMP))
  RCONST(16) = (ARR2(2.2D-13,602.0_dp,TEMP))
  RCONST(17) = (ARR2(1.2D-13,-2450.0_dp,TEMP))
  RCONST(18) = (ARR2(1.9D+8,390.0_dp,TEMP))
  RCONST(19) = (2.2D-10)
  RCONST(20) = (ARR2(1.6D-12,-940.0_dp,TEMP))
  RCONST(21) = (ARR2(1.4D-14,-580.0_dp,TEMP))
  RCONST(22) = (ARR2(1.3D-11,250.0_dp,TEMP))
  RCONST(23) = (ARR2(2.5D-14,-1230.0_dp,TEMP))
  RCONST(24) = (ARR2(5.3D-13,256.0_dp,TEMP))
  RCONST(25) = (1.3D-21)
  RCONST(26) = (ARR2(3.5D+14,-10897.0_dp,TEMP))
  RCONST(27) = (ARR2(1.8D-20,530.0_dp,TEMP))
  RCONST(28) = (4.4D-40)
  RCONST(29) = (ARR2(4.5D-13,806.0_dp,TEMP))
  RCONST(30) = (6.6D-12)
  RCONST(31) = (1.0D-20)
  RCONST(32) = (ARR2(1.0D-12,713.0_dp,TEMP))
  RCONST(33) = (ARR2(5.1D-15,1000.0_dp,TEMP))
  RCONST(34) = (ARR2(3.7D-12,240.0_dp,TEMP))
  RCONST(35) = (ARR2(1.2D-13,749.0_dp,TEMP))
  RCONST(36) = (ARR2(4.8D+13,-10121.0_dp,TEMP))
  RCONST(37) = (ARR2(1.3D-12,380.0_dp,TEMP))
  RCONST(38) = (ARR2(5.9D-14,1150.0_dp,TEMP))
  RCONST(39) = (ARR2(2.2D-38,5800.0_dp,TEMP))
  RCONST(40) = (ARR2(3.1D-12,-187.0_dp,TEMP))
  RCONST(41) = (2.2D-13)
  RCONST(42) = (1.0D-11)
  RCONST(43) = (ARR2(3.0D-11,-1550.0_dp,TEMP))
  RCONST(44) = (6.3D-16)
  RCONST(45) = (ARR2(1.2D-11,-986.0_dp,TEMP))
  RCONST(46) = (ARR2(7.0D-12,250.0_dp,TEMP))
  RCONST(47) = (2.5D-15)
  RCONST(48) = (ARR2(5.4D-12,250.0_dp,TEMP))
  RCONST(49) = (ARR2(8.0D-20,5500.0_dp,TEMP))
  RCONST(50) = (ARR2(9.4D+16,-14000.0_dp,TEMP))
  RCONST(51) = (2.0D-12)
  RCONST(52) = (6.5D-12)
  RCONST(53) = (ARR2(1.1D+2,-1710.0_dp,TEMP))
  RCONST(54) = (8.1D-13)
  RCONST(55) = (ARR2(1.0D+15,-8000.0_dp,TEMP))
  RCONST(56) = (1.6D+03)
  RCONST(57) = (1.5D-11)
  RCONST(58) = (ARR2(1.2D-11,-324.0_dp,TEMP))
  RCONST(59) = (ARR2(5.2D-12,504.0_dp,TEMP))
  RCONST(60) = (ARR2(1.4D-14,-2105.0_dp,TEMP))
  RCONST(61) = (7.7D-15)
  RCONST(62) = (ARR2(1.0D-11,-792.0_dp,TEMP))
  RCONST(63) = (ARR2(2.0D-12,411.0_dp,TEMP))
  RCONST(64) = (ARR2(1.3D-14,-2633.0_dp,TEMP))
  RCONST(65) = (ARR2(2.1D-12,322.0_dp,TEMP))
  RCONST(66) = (8.1D-12)
  RCONST(67) = 4.2
  RCONST(68) = (4.1D-11)
  RCONST(69) = (2.2D-11)
  RCONST(70) = (1.4D-11)
  RCONST(71) = (ARR2(1.7D-11,116.0_dp,TEMP))
  RCONST(72) = (3.0D-11)
  RCONST(73) = (ARR2(5.4D-17,-500.0_dp,TEMP))
  RCONST(74) = (1.70D-11)
  RCONST(75) = (1.80D-11)
  RCONST(76) = (9.6D-11)
  RCONST(77) = (1.2D-17)
  RCONST(78) = (3.2D-13)
  RCONST(79) = (8.1D-12)
  RCONST(80) = (ARR2(1.7D-14,1300.0_dp,TEMP))
  RCONST(81) = (6.8D-13)
END SUBROUTINE cbm4_Update_RCONST













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














END MODULE cbm4_UpdateRconstWRF

