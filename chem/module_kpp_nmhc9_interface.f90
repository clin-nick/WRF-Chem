







MODULE module_kpp_nmhc9_interf 


  USE module_state_description
  USE module_configure

  USE nmhc9_Parameters
  USE nmhc9_Precision
  USE nmhc9_UpdateRconstWRF
  USE nmhc9_Integrator

  USE module_wkppc_constants







     INTEGER, PARAMETER, PRIVATE :: Pj_o31d = 1 
     INTEGER, PARAMETER, PRIVATE :: Pj_o33p = 2 
     INTEGER, PARAMETER, PRIVATE :: Pj_no2 = 3 
     INTEGER, PARAMETER, PRIVATE :: Pj_no3o2 = 4 
     INTEGER, PARAMETER, PRIVATE :: Pj_no3o = 5 
     INTEGER, PARAMETER, PRIVATE :: Pj_hno2 = 6 
     INTEGER, PARAMETER, PRIVATE :: Pj_hno3 = 7 
     INTEGER, PARAMETER, PRIVATE :: Pj_hno4 = 8 
     INTEGER, PARAMETER, PRIVATE :: Pj_h2o2 = 9 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch2or = 10 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch2om = 11 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch3cho = 12 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch3coch3 = 13 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch3coc2h5 = 14 
     INTEGER, PARAMETER, PRIVATE :: Pj_hcocho = 15 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch3cocho = 16 
     INTEGER, PARAMETER, PRIVATE :: Pj_hcochest = 17 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch3o2h = 18 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch3coo2h = 19 
     INTEGER, PARAMETER, PRIVATE :: Pj_ch3ono2 = 20 
     INTEGER, PARAMETER, PRIVATE :: Pj_hcochob = 21 
     INTEGER, PARAMETER, PRIVATE :: Pj_macr = 22 
     INTEGER, PARAMETER, PRIVATE :: Pj_n2o5 = 23 
     INTEGER, PARAMETER, PRIVATE :: Pj_o2 = 24 
     INTEGER, PARAMETER, PRIVATE :: Pj_pan = 25 
     INTEGER, PARAMETER, PRIVATE :: Pj_acet = 26 
     INTEGER, PARAMETER, PRIVATE :: Pj_mglo = 27 
     INTEGER, PARAMETER, PRIVATE :: Pj_hno4_2 = 28 
     INTEGER, PARAMETER, PRIVATE :: Pj_clno2 = 29 
     INTEGER, PARAMETER, PRIVATE :: Pj_n2o = 30 
     INTEGER, PARAMETER, PRIVATE :: Pj_pooh = 31 
     INTEGER, PARAMETER, PRIVATE :: Pj_mpan = 32 
     INTEGER, PARAMETER, PRIVATE :: Pj_mvk = 33 
     INTEGER, PARAMETER, PRIVATE :: Pj_etooh = 34 
     INTEGER, PARAMETER, PRIVATE :: Pj_prooh = 35 
     INTEGER, PARAMETER, PRIVATE :: Pj_onitr = 36 
     INTEGER, PARAMETER, PRIVATE :: Pj_acetol = 37 
     INTEGER, PARAMETER, PRIVATE :: Pj_glyald = 38 
     INTEGER, PARAMETER, PRIVATE :: Pj_hyac = 39 
     INTEGER, PARAMETER, PRIVATE :: Pj_mek = 40 
     INTEGER, PARAMETER, PRIVATE :: Pj_open = 41 
     INTEGER, PARAMETER, PRIVATE :: Pj_gly = 42 
     INTEGER, PARAMETER, PRIVATE :: Pj_acetp = 43 
     INTEGER, PARAMETER, PRIVATE :: Pj_xooh = 44 
     INTEGER, PARAMETER, PRIVATE :: Pj_isooh = 45 
     INTEGER, PARAMETER, PRIVATE :: Pj_alkooh = 46 
     INTEGER, PARAMETER, PRIVATE :: Pj_mekooh = 47 
     INTEGER, PARAMETER, PRIVATE :: Pj_tolooh = 48 
     INTEGER, PARAMETER, PRIVATE :: Pj_terpooh = 49 
     INTEGER, PARAMETER, PRIVATE :: Pj_cl2 = 50 
     INTEGER, PARAMETER, PRIVATE :: Pj_hocl = 51 
     INTEGER, PARAMETER, PRIVATE :: Pj_fmcl = 52 
 


CONTAINS 

SUBROUTINE  nmhc9_interface( &

    chem, id, dtstepc,config_flags, & 
    p_phy,t_phy,rho_phy,moist, aero_srf_area, &
    vdrog3, ldrog, vdrog3_vbs, ldrog_vbs, &

             addt, addx, addc, etep, oltp, & 
             olip, cslp, limp, hc5p, hc8p, & 
             tolp, xylp, apip, isop, hc3p, & 
             ethp, o3p, tco3, mo2, o1d, & 
             olnn, olnd, rpho, xo2, ketp, & 
             xno2, ol2p, oln, macp, hocoo, & 
             bzno2_o, bz_o, tbu_o,  & 
             ph_o31d, ph_o33p, ph_no2, ph_no3o2, ph_no3o, & 
             ph_hno2, ph_hno3, ph_hno4, ph_h2o2, ph_ch2or, & 
             ph_ch2om, ph_ch3cho, ph_ch3coch3, ph_ch3coc2h5, ph_hcocho, & 
             ph_ch3cocho, ph_hcochest, ph_ch3o2h, ph_ch3coo2h, ph_ch3ono2, & 
             ph_hcochob, ph_macr, ph_n2o5, ph_o2, ph_pan, & 
             ph_acet, ph_mglo, ph_hno4_2, ph_clno2, ph_n2o, & 
             ph_pooh, ph_mpan, ph_mvk, ph_etooh, ph_prooh, & 
             ph_onitr, ph_acetol, ph_glyald, ph_hyac, ph_mek, & 
             ph_open, ph_gly, ph_acetp, ph_xooh, ph_isooh, & 
             ph_alkooh, ph_mekooh, ph_tolooh, ph_terpooh, ph_cl2, & 
             ph_hocl, ph_fmcl,  & 
              ids,ide, jds,jde, kds,kde,         &
              ims,ime, jms,jme, kms,kme,         &
              its,ite, jts,jte, kts,kte)


    IMPLICIT NONE


    INTEGER,      INTENT(IN   ) ::    &
                      ids,ide, jds,jde, kds,kde,      & 
                      ims,ime, jms,jme, kms,kme,      & 
                      its,ite, jts,jte, kts,kte 



    INTEGER,      INTENT(IN   ) :: id

    REAL,      INTENT(IN   ) :: dtstepc 

    TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags


    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),            &
	        INTENT(INOUT ) ::                               chem 


    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),            &
	        INTENT(IN ) ::                               moist 

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_aero_srf_area ), &
	        INTENT(INOUT) ::                             aero_srf_area 

     REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),              &
             INTENT(IN   ) ::                                      &          
                                                         p_phy,    &
                                                         t_phy,    &  
                                                         rho_phy        




    INTEGER, INTENT ( IN ) :: ldrog
  
    REAL,      INTENT(INOUT) ::                                     &
                      vdrog3(ims:ime,kms:kme-0,jms:jme,ldrog)


    INTEGER, INTENT ( IN ) :: ldrog_vbs
  
    REAL,      INTENT(INOUT) ::                                     &
                      vdrog3_vbs(ims:ime,kms:kme-0,jms:jme,ldrog_vbs)







      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),   & 
         INTENT(INOUT ) ::  & 
               addt, addx, addc, etep, oltp, & 
               olip, cslp, limp, hc5p, hc8p, & 
               tolp, xylp, apip, isop, hc3p, & 
               ethp, o3p, tco3, mo2, o1d, & 
               olnn, olnd, rpho, xo2, ketp, & 
               xno2, ol2p, oln, macp, hocoo, & 
               bzno2_o, bz_o, tbu_o




      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),   & 
         INTENT(INOUT ) ::  & 
               ph_o31d, ph_o33p, ph_no2, ph_no3o2, ph_no3o, & 
               ph_hno2, ph_hno3, ph_hno4, ph_h2o2, ph_ch2or, & 
               ph_ch2om, ph_ch3cho, ph_ch3coch3, ph_ch3coc2h5, ph_hcocho, & 
               ph_ch3cocho, ph_hcochest, ph_ch3o2h, ph_ch3coo2h, ph_ch3ono2, & 
               ph_hcochob, ph_macr, ph_n2o5, ph_o2, ph_pan, & 
               ph_acet, ph_mglo, ph_hno4_2, ph_clno2, ph_n2o, & 
               ph_pooh, ph_mpan, ph_mvk, ph_etooh, ph_prooh, & 
               ph_onitr, ph_acetol, ph_glyald, ph_hyac, ph_mek, & 
               ph_open, ph_gly, ph_acetp, ph_xooh, ph_isooh, & 
               ph_alkooh, ph_mekooh, ph_tolooh, ph_terpooh, ph_cl2, & 
               ph_hocl, ph_fmcl 




    INTEGER, PARAMETER :: njv=52
    REAL(KIND=dp), DIMENSION(njv) :: jv


    REAL(KIND=dp):: TIME_START
    REAL(KIND=dp):: TIME_END

    INTEGER, DIMENSION(20) :: ICNTRL 
    REAL(KIND=dp), DIMENSION(20) :: RCNTRL
    INTEGER, DIMENSION(20) :: ISTATUS 
    REAL(KIND=dp), DIMENSION(20) :: RSTATUS
    INTEGER :: IERR_U

    REAL(KIND=dp), DIMENSION(NREACT):: RCONST 

    REAL(KIND=dp), DIMENSION(NVAR) :: var
    REAL(KIND=dp), DIMENSION(NFIX) :: fix

     
    REAL(KIND=dp)   :: TEMP 

    REAL(KIND=dp), DIMENSION(NSPEC)  :: ATOL, RTOL
    REAL(KIND=dp), DIMENSION(NREACT) :: IRR_WRK

    REAL(KIND=dp) :: conv, oconv 
    REAL(KIND=dp) :: C_M 

    INTEGER :: i,j,k,n 
 



 

 REAL(KIND=dp)   :: PRESS



 


      DO n=1, 20
         ICNTRL(n) = 0
         RCNTRL(n) = 0._dp
      END DO


      DO n=1, NSPEC
         ATOL(n) = REAL(atols, KIND=dp)
         RTOL(n) = REAL(rtols, KIND=dp)
      END DO


      TIME_START = 0.0_dp 
      TIME_END =  REAL(dtstepc, KIND=dp) 















   ICNTRL(1) = 1


   ICNTRL(3) = 2


   RCNTRL(3) = 0.01_dp * TIME_END















    DO j=jts, jte
    DO k=kts, kte
    DO i=its, ite


      
    FIX(indf_M)  = REAL(dens2con_a * rho_phy(i,k,j), KIND=dp)
    C_M = FIX(indf_M)

      
    FIX(indf_H2O) = REAL(dens2con_w * moist(i,k,j,P_QV) * rho_phy(i,k,j), KIND=dp)


      
    TEMP = REAL(t_phy(i,k,j), KIND=dp)

      
     conv=1.E-6_dp*dens2con_a*rho_phy(i,k,j)
     oconv = 1.E0_dp/conv


    jv(Pj_o31d) = REAL(ph_o31d(i,k,j)/60., KIND=dp) 
    jv(Pj_o33p) = REAL(ph_o33p(i,k,j)/60., KIND=dp) 
    jv(Pj_no2) = REAL(ph_no2(i,k,j)/60., KIND=dp) 
    jv(Pj_no3o2) = REAL(ph_no3o2(i,k,j)/60., KIND=dp) 
    jv(Pj_no3o) = REAL(ph_no3o(i,k,j)/60., KIND=dp) 
    jv(Pj_hno2) = REAL(ph_hno2(i,k,j)/60., KIND=dp) 
    jv(Pj_hno3) = REAL(ph_hno3(i,k,j)/60., KIND=dp) 
    jv(Pj_hno4) = REAL(ph_hno4(i,k,j)/60., KIND=dp) 
    jv(Pj_h2o2) = REAL(ph_h2o2(i,k,j)/60., KIND=dp) 
    jv(Pj_ch2or) = REAL(ph_ch2or(i,k,j)/60., KIND=dp) 
    jv(Pj_ch2om) = REAL(ph_ch2om(i,k,j)/60., KIND=dp) 
    jv(Pj_ch3cho) = REAL(ph_ch3cho(i,k,j)/60., KIND=dp) 
    jv(Pj_ch3coch3) = REAL(ph_ch3coch3(i,k,j)/60., KIND=dp) 
    jv(Pj_ch3coc2h5) = REAL(ph_ch3coc2h5(i,k,j)/60., KIND=dp) 
    jv(Pj_hcocho) = REAL(ph_hcocho(i,k,j)/60., KIND=dp) 
    jv(Pj_ch3cocho) = REAL(ph_ch3cocho(i,k,j)/60., KIND=dp) 
    jv(Pj_hcochest) = REAL(ph_hcochest(i,k,j)/60., KIND=dp) 
    jv(Pj_ch3o2h) = REAL(ph_ch3o2h(i,k,j)/60., KIND=dp) 
    jv(Pj_ch3coo2h) = REAL(ph_ch3coo2h(i,k,j)/60., KIND=dp) 
    jv(Pj_ch3ono2) = REAL(ph_ch3ono2(i,k,j)/60., KIND=dp) 
    jv(Pj_hcochob) = REAL(ph_hcochob(i,k,j)/60., KIND=dp) 
    jv(Pj_macr) = REAL(ph_macr(i,k,j)/60., KIND=dp) 
    jv(Pj_n2o5) = REAL(ph_n2o5(i,k,j)/60., KIND=dp) 
    jv(Pj_o2) = REAL(ph_o2(i,k,j)/60., KIND=dp) 
    jv(Pj_pan) = REAL(ph_pan(i,k,j)/60., KIND=dp) 
    jv(Pj_acet) = REAL(ph_acet(i,k,j)/60., KIND=dp) 
    jv(Pj_mglo) = REAL(ph_mglo(i,k,j)/60., KIND=dp) 
    jv(Pj_hno4_2) = REAL(ph_hno4_2(i,k,j)/60., KIND=dp) 
    jv(Pj_clno2) = REAL(ph_clno2(i,k,j)/60., KIND=dp) 
    jv(Pj_n2o) = REAL(ph_n2o(i,k,j)/60., KIND=dp) 
    jv(Pj_pooh) = REAL(ph_pooh(i,k,j)/60., KIND=dp) 
    jv(Pj_mpan) = REAL(ph_mpan(i,k,j)/60., KIND=dp) 
    jv(Pj_mvk) = REAL(ph_mvk(i,k,j)/60., KIND=dp) 
    jv(Pj_etooh) = REAL(ph_etooh(i,k,j)/60., KIND=dp) 
    jv(Pj_prooh) = REAL(ph_prooh(i,k,j)/60., KIND=dp) 
    jv(Pj_onitr) = REAL(ph_onitr(i,k,j)/60., KIND=dp) 
    jv(Pj_acetol) = REAL(ph_acetol(i,k,j)/60., KIND=dp) 
    jv(Pj_glyald) = REAL(ph_glyald(i,k,j)/60., KIND=dp) 
    jv(Pj_hyac) = REAL(ph_hyac(i,k,j)/60., KIND=dp) 
    jv(Pj_mek) = REAL(ph_mek(i,k,j)/60., KIND=dp) 
    jv(Pj_open) = REAL(ph_open(i,k,j)/60., KIND=dp) 
    jv(Pj_gly) = REAL(ph_gly(i,k,j)/60., KIND=dp) 
    jv(Pj_acetp) = REAL(ph_acetp(i,k,j)/60., KIND=dp) 
    jv(Pj_xooh) = REAL(ph_xooh(i,k,j)/60., KIND=dp) 
    jv(Pj_isooh) = REAL(ph_isooh(i,k,j)/60., KIND=dp) 
    jv(Pj_alkooh) = REAL(ph_alkooh(i,k,j)/60., KIND=dp) 
    jv(Pj_mekooh) = REAL(ph_mekooh(i,k,j)/60., KIND=dp) 
    jv(Pj_tolooh) = REAL(ph_tolooh(i,k,j)/60., KIND=dp) 
    jv(Pj_terpooh) = REAL(ph_terpooh(i,k,j)/60., KIND=dp) 
    jv(Pj_cl2) = REAL(ph_cl2(i,k,j)/60., KIND=dp) 
    jv(Pj_hocl) = REAL(ph_hocl(i,k,j)/60., KIND=dp) 
    jv(Pj_fmcl) = REAL(ph_fmcl(i,k,j)/60., KIND=dp) 
 


    var(ind_O3) = conv  * REAL( MAX(chem(i,k,j,P_o3),0.), KIND=dp)  
    var(ind_CH4) = conv  * REAL( MAX(chem(i,k,j,P_ch4),0.), KIND=dp)  
    var(ind_CO) = conv  * REAL( MAX(chem(i,k,j,P_co),0.), KIND=dp)  
    var(ind_H2O2) = conv  * REAL( MAX(chem(i,k,j,P_h2o2),0.), KIND=dp)  
    var(ind_HNO3) = conv  * REAL( MAX(chem(i,k,j,P_hno3),0.), KIND=dp)  
    var(ind_MeOOH) = conv  * REAL( MAX(chem(i,k,j,P_OP1),0.), KIND=dp)  
    var(ind_HCHO) = conv  * REAL( MAX(chem(i,k,j,P_hcho),0.), KIND=dp)  
    var(ind_MeOH) = conv  * REAL( MAX(chem(i,k,j,P_meoh),0.), KIND=dp)  
    var(ind_MeO2NO2) = conv  * REAL( MAX(chem(i,k,j,P_meo2no2),0.), KIND=dp)  
    var(ind_ISOP) = conv  * REAL( MAX(chem(i,k,j,P_ISOPR),0.), KIND=dp)  
    var(ind_ISOOH) = conv  * REAL( MAX(chem(i,k,j,P_isooh),0.), KIND=dp)  
    var(ind_MVK) = conv  * REAL( MAX(chem(i,k,j,P_mvk),0.), KIND=dp)  
    var(ind_MVKOOH) = conv  * REAL( MAX(chem(i,k,j,P_mvkooh),0.), KIND=dp)  
    var(ind_MGLO) = conv  * REAL( MAX(chem(i,k,j,P_mglo),0.), KIND=dp)  
    var(ind_PAA) = conv  * REAL( MAX(chem(i,k,j,P_paa),0.), KIND=dp)  
    var(ind_PAN) = conv  * REAL( MAX(chem(i,k,j,P_pan),0.), KIND=dp)  
    var(ind_MPAN) = conv  * REAL( MAX(chem(i,k,j,P_mpan),0.), KIND=dp)  
    var(ind_ISON) = conv  * REAL( MAX(chem(i,k,j,P_ison),0.), KIND=dp)  
    var(ind_HCOOH) = conv  * REAL( MAX(chem(i,k,j,P_hcooh),0.), KIND=dp)  
    var(ind_CH3COOH) = conv  * REAL( MAX(chem(i,k,j,P_ch3cooh),0.), KIND=dp)  
    var(ind_ACETOL) = conv  * REAL( MAX(chem(i,k,j,P_acetol),0.), KIND=dp)  
    var(ind_NACA) = conv  * REAL( MAX(chem(i,k,j,P_naca),0.), KIND=dp)  
    var(ind_C2H6) = conv  * REAL( MAX(chem(i,k,j,P_c2h6),0.), KIND=dp)  
    var(ind_EtOOH) = conv  * REAL( MAX(chem(i,k,j,P_etooh),0.), KIND=dp)  
    var(ind_ALD) = conv  * REAL( MAX(chem(i,k,j,P_ald),0.), KIND=dp)  
    var(ind_C3H8) = conv  * REAL( MAX(chem(i,k,j,P_c3h8),0.), KIND=dp)  
    var(ind_PrOOH) = conv  * REAL( MAX(chem(i,k,j,P_prooh),0.), KIND=dp)  
    var(ind_ACET) = conv  * REAL( MAX(chem(i,k,j,P_acet),0.), KIND=dp)  
    var(ind_ACETP) = conv  * REAL( MAX(chem(i,k,j,P_acetp),0.), KIND=dp)  
    var(ind_PrONO2) = conv  * REAL( MAX(chem(i,k,j,P_prono2),0.), KIND=dp)  
    var(ind_C3H6) = conv  * REAL( MAX(chem(i,k,j,P_c3h6),0.), KIND=dp)  
    var(ind_C3H6OOH) = conv  * REAL( MAX(chem(i,k,j,P_c3h6ooh),0.), KIND=dp)  
    var(ind_C2H4) = conv  * REAL( MAX(chem(i,k,j,P_c2h4),0.), KIND=dp)  
    var(ind_C4H10) = conv  * REAL( MAX(chem(i,k,j,P_c4h10),0.), KIND=dp)  
    var(ind_C4H9OOH) = conv  * REAL( MAX(chem(i,k,j,P_c4h9ooh),0.), KIND=dp)  
    var(ind_ONIT) = conv  * REAL( MAX(chem(i,k,j,P_onit),0.), KIND=dp)  
    var(ind_MEK) = conv  * REAL( MAX(chem(i,k,j,P_mek),0.), KIND=dp)  
    var(ind_MEKOOH) = conv  * REAL( MAX(chem(i,k,j,P_mekooh),0.), KIND=dp)  
    var(ind_MeCOCO) = conv  * REAL( MAX(chem(i,k,j,P_mecoco),0.), KIND=dp)  
    var(ind_O1D) =  conv * REAL( MAX(o1d(i,k,j),0.), KIND=dp)  
    var(ind_OH) = conv  * REAL( MAX(chem(i,k,j,P_HO),0.), KIND=dp)  
    var(ind_HO2) = conv  * REAL( MAX(chem(i,k,j,P_ho2),0.), KIND=dp)  
    var(ind_NO) = conv  * REAL( MAX(chem(i,k,j,P_no),0.), KIND=dp)  
    var(ind_NO2) = conv  * REAL( MAX(chem(i,k,j,P_no2),0.), KIND=dp)  
    var(ind_NO3) = conv  * REAL( MAX(chem(i,k,j,P_no3),0.), KIND=dp)  
    var(ind_N2O5) = conv  * REAL( MAX(chem(i,k,j,P_n2o5),0.), KIND=dp)  
    var(ind_HNO4) = conv  * REAL( MAX(chem(i,k,j,P_hno4),0.), KIND=dp)  
    var(ind_MeO2) = conv  * REAL( MAX(chem(i,k,j,P_meo2),0.), KIND=dp)  
    var(ind_PA) = conv  * REAL( MAX(chem(i,k,j,P_pa),0.), KIND=dp)  
    var(ind_ISO2) = conv  * REAL( MAX(chem(i,k,j,P_iso2),0.), KIND=dp)  
    var(ind_MVKO2) = conv  * REAL( MAX(chem(i,k,j,P_mvko2),0.), KIND=dp)  
    var(ind_EtO2) = conv  * REAL( MAX(chem(i,k,j,P_eto2),0.), KIND=dp)  
    var(ind_PrO2) = conv  * REAL( MAX(chem(i,k,j,P_pro2),0.), KIND=dp)  
    var(ind_ACETO2) = conv  * REAL( MAX(chem(i,k,j,P_aceto2),0.), KIND=dp)  
    var(ind_C3H6O2) = conv  * REAL( MAX(chem(i,k,j,P_c3h6o2),0.), KIND=dp)  
    var(ind_C4H9O2) = conv  * REAL( MAX(chem(i,k,j,P_c4h9o2),0.), KIND=dp)  
    var(ind_MEKO2) = conv  * REAL( MAX(chem(i,k,j,P_meko2),0.), KIND=dp)  


     
    PRESS = REAL(p_phy(i,k,j), KIND=dp)








   CALL nmhc9_Update_Rconst(  &


PRESS, &


             jv, njv, &
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
             C_M, FIX(indf_H2O), TEMP & 

)








  CALL nmhc9_INTEGRATE(TIME_START, TIME_END, &  
          FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK, & 
          ICNTRL_U=icntrl, RCNTRL_U=rcntrl  )







    chem(i,k,j,P_o3) = MAX ( REAL (oconv * var(ind_O3), KIND=sp), 0.)  
    chem(i,k,j,P_ch4) = MAX ( REAL (oconv * var(ind_CH4), KIND=sp), 0.)  
    chem(i,k,j,P_co) = MAX ( REAL (oconv * var(ind_CO), KIND=sp), 0.)  
    chem(i,k,j,P_h2o2) = MAX ( REAL (oconv * var(ind_H2O2), KIND=sp), 0.)  
    chem(i,k,j,P_hno3) = MAX ( REAL (oconv * var(ind_HNO3), KIND=sp), 0.)  
    chem(i,k,j,P_OP1) = MAX ( REAL (oconv * var(ind_MeOOH), KIND=sp), 0.)  
    chem(i,k,j,P_hcho) = MAX ( REAL (oconv * var(ind_HCHO), KIND=sp), 0.)  
    chem(i,k,j,P_meoh) = MAX ( REAL (oconv * var(ind_MeOH), KIND=sp), 0.)  
    chem(i,k,j,P_meo2no2) = MAX ( REAL (oconv * var(ind_MeO2NO2), KIND=sp), 0.)  
    chem(i,k,j,P_ISOPR) = MAX ( REAL (oconv * var(ind_ISOP), KIND=sp), 0.)  
    chem(i,k,j,P_isooh) = MAX ( REAL (oconv * var(ind_ISOOH), KIND=sp), 0.)  
    chem(i,k,j,P_mvk) = MAX ( REAL (oconv * var(ind_MVK), KIND=sp), 0.)  
    chem(i,k,j,P_mvkooh) = MAX ( REAL (oconv * var(ind_MVKOOH), KIND=sp), 0.)  
    chem(i,k,j,P_mglo) = MAX ( REAL (oconv * var(ind_MGLO), KIND=sp), 0.)  
    chem(i,k,j,P_paa) = MAX ( REAL (oconv * var(ind_PAA), KIND=sp), 0.)  
    chem(i,k,j,P_pan) = MAX ( REAL (oconv * var(ind_PAN), KIND=sp), 0.)  
    chem(i,k,j,P_mpan) = MAX ( REAL (oconv * var(ind_MPAN), KIND=sp), 0.)  
    chem(i,k,j,P_ison) = MAX ( REAL (oconv * var(ind_ISON), KIND=sp), 0.)  
    chem(i,k,j,P_hcooh) = MAX ( REAL (oconv * var(ind_HCOOH), KIND=sp), 0.)  
    chem(i,k,j,P_ch3cooh) = MAX ( REAL (oconv * var(ind_CH3COOH), KIND=sp), 0.)  
    chem(i,k,j,P_acetol) = MAX ( REAL (oconv * var(ind_ACETOL), KIND=sp), 0.)  
    chem(i,k,j,P_naca) = MAX ( REAL (oconv * var(ind_NACA), KIND=sp), 0.)  
    chem(i,k,j,P_c2h6) = MAX ( REAL (oconv * var(ind_C2H6), KIND=sp), 0.)  
    chem(i,k,j,P_etooh) = MAX ( REAL (oconv * var(ind_EtOOH), KIND=sp), 0.)  
    chem(i,k,j,P_ald) = MAX ( REAL (oconv * var(ind_ALD), KIND=sp), 0.)  
    chem(i,k,j,P_c3h8) = MAX ( REAL (oconv * var(ind_C3H8), KIND=sp), 0.)  
    chem(i,k,j,P_prooh) = MAX ( REAL (oconv * var(ind_PrOOH), KIND=sp), 0.)  
    chem(i,k,j,P_acet) = MAX ( REAL (oconv * var(ind_ACET), KIND=sp), 0.)  
    chem(i,k,j,P_acetp) = MAX ( REAL (oconv * var(ind_ACETP), KIND=sp), 0.)  
    chem(i,k,j,P_prono2) = MAX ( REAL (oconv * var(ind_PrONO2), KIND=sp), 0.)  
    chem(i,k,j,P_c3h6) = MAX ( REAL (oconv * var(ind_C3H6), KIND=sp), 0.)  
    chem(i,k,j,P_c3h6ooh) = MAX ( REAL (oconv * var(ind_C3H6OOH), KIND=sp), 0.)  
    chem(i,k,j,P_c2h4) = MAX ( REAL (oconv * var(ind_C2H4), KIND=sp), 0.)  
    chem(i,k,j,P_c4h10) = MAX ( REAL (oconv * var(ind_C4H10), KIND=sp), 0.)  
    chem(i,k,j,P_c4h9ooh) = MAX ( REAL (oconv * var(ind_C4H9OOH), KIND=sp), 0.)  
    chem(i,k,j,P_onit) = MAX ( REAL (oconv * var(ind_ONIT), KIND=sp), 0.)  
    chem(i,k,j,P_mek) = MAX ( REAL (oconv * var(ind_MEK), KIND=sp), 0.)  
    chem(i,k,j,P_mekooh) = MAX ( REAL (oconv * var(ind_MEKOOH), KIND=sp), 0.)  
    chem(i,k,j,P_mecoco) = MAX ( REAL (oconv * var(ind_MeCOCO), KIND=sp), 0.)  
   o1d(i,k,j) = MAX (REAL (oconv * var(ind_O1D) , KIND=sp),0.) 
    chem(i,k,j,P_HO) = MAX ( REAL (oconv * var(ind_OH), KIND=sp), 0.)  
    chem(i,k,j,P_ho2) = MAX ( REAL (oconv * var(ind_HO2), KIND=sp), 0.)  
    chem(i,k,j,P_no) = MAX ( REAL (oconv * var(ind_NO), KIND=sp), 0.)  
    chem(i,k,j,P_no2) = MAX ( REAL (oconv * var(ind_NO2), KIND=sp), 0.)  
    chem(i,k,j,P_no3) = MAX ( REAL (oconv * var(ind_NO3), KIND=sp), 0.)  
    chem(i,k,j,P_n2o5) = MAX ( REAL (oconv * var(ind_N2O5), KIND=sp), 0.)  
    chem(i,k,j,P_hno4) = MAX ( REAL (oconv * var(ind_HNO4), KIND=sp), 0.)  
    chem(i,k,j,P_meo2) = MAX ( REAL (oconv * var(ind_MeO2), KIND=sp), 0.)  
    chem(i,k,j,P_pa) = MAX ( REAL (oconv * var(ind_PA), KIND=sp), 0.)  
    chem(i,k,j,P_iso2) = MAX ( REAL (oconv * var(ind_ISO2), KIND=sp), 0.)  
    chem(i,k,j,P_mvko2) = MAX ( REAL (oconv * var(ind_MVKO2), KIND=sp), 0.)  
    chem(i,k,j,P_eto2) = MAX ( REAL (oconv * var(ind_EtO2), KIND=sp), 0.)  
    chem(i,k,j,P_pro2) = MAX ( REAL (oconv * var(ind_PrO2), KIND=sp), 0.)  
    chem(i,k,j,P_aceto2) = MAX ( REAL (oconv * var(ind_ACETO2), KIND=sp), 0.)  
    chem(i,k,j,P_c3h6o2) = MAX ( REAL (oconv * var(ind_C3H6O2), KIND=sp), 0.)  
    chem(i,k,j,P_c4h9o2) = MAX ( REAL (oconv * var(ind_C4H9O2), KIND=sp), 0.)  
    chem(i,k,j,P_meko2) = MAX ( REAL (oconv * var(ind_MEKO2), KIND=sp), 0.)  



    END DO
    END DO
    END DO









END SUBROUTINE  nmhc9_interface


END MODULE module_kpp_nmhc9_interf 




