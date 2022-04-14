







MODULE module_kpp_mozcart_interf 


  USE module_state_description
  USE module_configure

  USE mozcart_Parameters
  USE mozcart_Precision
  USE mozcart_UpdateRconstWRF
  USE mozcart_Integrator

  USE module_wkppc_constants

  USE module_irr_diag






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

SUBROUTINE  mozcart_interface( &

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
              its,ite, jts,jte, kts,kte, dm,num_irr_diag,irr_rates)


    IMPLICIT NONE


    INTEGER,      INTENT(IN   ) ::    &
                      ids,ide, jds,jde, kds,kde,      & 
                      ims,ime, jms,jme, kms,kme,      & 
                      its,ite, jts,jte, kts,kte 



          INTEGER,      INTENT(IN   ) ::    dm
          INTEGER,      INTENT(IN   ) ::    num_irr_diag



    REAL,         INTENT(INOUT) :: irr_rates(ims:ime,kms:kme,jms:jme,num_irr_diag)

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
 



 

  REAL(kind=dp) :: es, qvs, rh


 


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
    var(ind_O) = conv  * REAL( MAX(chem(i,k,j,P_o),0.), KIND=dp)  
    var(ind_O1D_CB4) = conv  * REAL( MAX(chem(i,k,j,P_o1d_cb4),0.), KIND=dp)  
    var(ind_N2O) = conv  * REAL( MAX(chem(i,k,j,P_n2o),0.), KIND=dp)  
    var(ind_NO) = conv  * REAL( MAX(chem(i,k,j,P_no),0.), KIND=dp)  
    var(ind_NO2) = conv  * REAL( MAX(chem(i,k,j,P_no2),0.), KIND=dp)  
    var(ind_NO3) = conv  * REAL( MAX(chem(i,k,j,P_no3),0.), KIND=dp)  
    var(ind_NH3) = conv  * REAL( MAX(chem(i,k,j,P_nh3),0.), KIND=dp)  
    var(ind_HNO3) = conv  * REAL( MAX(chem(i,k,j,P_hno3),0.), KIND=dp)  
    var(ind_HO2NO2) = conv  * REAL( MAX(chem(i,k,j,P_HNO4),0.), KIND=dp)  
    var(ind_N2O5) = conv  * REAL( MAX(chem(i,k,j,P_n2o5),0.), KIND=dp)  
    var(ind_H2) = conv  * REAL( MAX(chem(i,k,j,P_h2),0.), KIND=dp)  
    var(ind_OH) = conv  * REAL( MAX(chem(i,k,j,P_HO),0.), KIND=dp)  
    var(ind_HO2) = conv  * REAL( MAX(chem(i,k,j,P_ho2),0.), KIND=dp)  
    var(ind_H2O2) = conv  * REAL( MAX(chem(i,k,j,P_h2o2),0.), KIND=dp)  
    var(ind_CH4) = conv  * REAL( MAX(chem(i,k,j,P_ch4),0.), KIND=dp)  
    var(ind_CO) = conv  * REAL( MAX(chem(i,k,j,P_co),0.), KIND=dp)  
    var(ind_CH3O2) = conv  * REAL( MAX(chem(i,k,j,P_ch3o2),0.), KIND=dp)  
    var(ind_CH3OOH) = conv  * REAL( MAX(chem(i,k,j,P_ch3ooh),0.), KIND=dp)  
    var(ind_CH2O) = conv  * REAL( MAX(chem(i,k,j,P_HCHO),0.), KIND=dp)  
    var(ind_CH3OH) = conv  * REAL( MAX(chem(i,k,j,P_ch3oh),0.), KIND=dp)  
    var(ind_C2H4) = conv  * REAL( MAX(chem(i,k,j,P_c2h4),0.), KIND=dp)  
    var(ind_EO) = conv  * REAL( MAX(chem(i,k,j,P_eo),0.), KIND=dp)  
    var(ind_EO2) = conv  * REAL( MAX(chem(i,k,j,P_eo2),0.), KIND=dp)  
    var(ind_CH3CHO) = conv  * REAL( MAX(chem(i,k,j,P_ALD),0.), KIND=dp)  
    var(ind_CH3COOH) = conv  * REAL( MAX(chem(i,k,j,P_ch3cooh),0.), KIND=dp)  
    var(ind_CH3COCH3) = conv  * REAL( MAX(chem(i,k,j,P_ACET),0.), KIND=dp)  
    var(ind_CH3COCHO) = conv  * REAL( MAX(chem(i,k,j,P_MGLY),0.), KIND=dp)  
    var(ind_CH3CO3) = conv  * REAL( MAX(chem(i,k,j,P_ACO3),0.), KIND=dp)  
    var(ind_CH3COOOH) = conv  * REAL( MAX(chem(i,k,j,P_PAA),0.), KIND=dp)  
    var(ind_GLYOXAL) = conv  * REAL( MAX(chem(i,k,j,P_GLY),0.), KIND=dp)  
    var(ind_PO2) = conv  * REAL( MAX(chem(i,k,j,P_po2),0.), KIND=dp)  
    var(ind_POOH) = conv  * REAL( MAX(chem(i,k,j,P_C3H6OOH),0.), KIND=dp)  
    var(ind_PAN) = conv  * REAL( MAX(chem(i,k,j,P_pan),0.), KIND=dp)  
    var(ind_MPAN) = conv  * REAL( MAX(chem(i,k,j,P_mpan),0.), KIND=dp)  
    var(ind_MCO3) = conv  * REAL( MAX(chem(i,k,j,P_mco3),0.), KIND=dp)  
    var(ind_MACR) = conv  * REAL( MAX(chem(i,k,j,P_macr),0.), KIND=dp)  
    var(ind_MACRO2) = conv  * REAL( MAX(chem(i,k,j,P_MVKO2),0.), KIND=dp)  
    var(ind_MACROOH) = conv  * REAL( MAX(chem(i,k,j,P_MVKOOH),0.), KIND=dp)  
    var(ind_MVK) = conv  * REAL( MAX(chem(i,k,j,P_mvk),0.), KIND=dp)  
    var(ind_C2H6) = conv  * REAL( MAX(chem(i,k,j,P_c2h6),0.), KIND=dp)  
    var(ind_C3H6) = conv  * REAL( MAX(chem(i,k,j,P_c3h6),0.), KIND=dp)  
    var(ind_C3H8) = conv  * REAL( MAX(chem(i,k,j,P_c3h8),0.), KIND=dp)  
    var(ind_C2H5OH) = conv  * REAL( MAX(chem(i,k,j,P_c2h5oh),0.), KIND=dp)  
    var(ind_C2H5OOH) = conv  * REAL( MAX(chem(i,k,j,P_ETOOH),0.), KIND=dp)  
    var(ind_C3H7O2) = conv  * REAL( MAX(chem(i,k,j,P_PRO2),0.), KIND=dp)  
    var(ind_C3H7OOH) = conv  * REAL( MAX(chem(i,k,j,P_PROOH),0.), KIND=dp)  
    var(ind_C10H16) = conv  * REAL( MAX(chem(i,k,j,P_c10h16),0.), KIND=dp)  
    var(ind_RO2) = conv  * REAL( MAX(chem(i,k,j,P_ACETO2),0.), KIND=dp)  
    var(ind_ROOH) = conv  * REAL( MAX(chem(i,k,j,P_ACETP),0.), KIND=dp)  
    var(ind_ONIT) = conv  * REAL( MAX(chem(i,k,j,P_onit),0.), KIND=dp)  
    var(ind_ONITR) = conv  * REAL( MAX(chem(i,k,j,P_onitr),0.), KIND=dp)  
    var(ind_ISOP) = conv  * REAL( MAX(chem(i,k,j,P_ISOPR),0.), KIND=dp)  
    var(ind_ISOPO2) = conv  * REAL( MAX(chem(i,k,j,P_ISO2),0.), KIND=dp)  
    var(ind_ISOPOOH) = conv  * REAL( MAX(chem(i,k,j,P_ISOOH),0.), KIND=dp)  
    var(ind_ISOPNO3) = conv  * REAL( MAX(chem(i,k,j,P_ISOPN),0.), KIND=dp)  
    var(ind_HYAC) = conv  * REAL( MAX(chem(i,k,j,P_ACETOL),0.), KIND=dp)  
    var(ind_GLYALD) = conv  * REAL( MAX(chem(i,k,j,P_glyald),0.), KIND=dp)  
    var(ind_HYDRALD) = conv  * REAL( MAX(chem(i,k,j,P_hydrald),0.), KIND=dp)  
    var(ind_ENEO2) = conv  * REAL( MAX(chem(i,k,j,P_eneo2),0.), KIND=dp)  
    var(ind_MEK) = conv  * REAL( MAX(chem(i,k,j,P_mek),0.), KIND=dp)  
    var(ind_MEKO2) = conv  * REAL( MAX(chem(i,k,j,P_meko2),0.), KIND=dp)  
    var(ind_C2H5O2) = conv  * REAL( MAX(chem(i,k,j,P_ETO2),0.), KIND=dp)  
    var(ind_BIGENE) = conv  * REAL( MAX(chem(i,k,j,P_bigene),0.), KIND=dp)  
    var(ind_BIGALD) = conv  * REAL( MAX(chem(i,k,j,P_OPEN),0.), KIND=dp)  
    var(ind_BIGALK) = conv  * REAL( MAX(chem(i,k,j,P_bigalk),0.), KIND=dp)  
    var(ind_ALKO2) = conv  * REAL( MAX(chem(i,k,j,P_alko2),0.), KIND=dp)  
    var(ind_ALKOOH) = conv  * REAL( MAX(chem(i,k,j,P_alkooh),0.), KIND=dp)  
    var(ind_MEKOOH) = conv  * REAL( MAX(chem(i,k,j,P_mekooh),0.), KIND=dp)  
    var(ind_TOLUENE) = conv  * REAL( MAX(chem(i,k,j,P_TOL),0.), KIND=dp)  
    var(ind_TOLO2) = conv  * REAL( MAX(chem(i,k,j,P_TO2),0.), KIND=dp)  
    var(ind_TOLOOH) = conv  * REAL( MAX(chem(i,k,j,P_tolooh),0.), KIND=dp)  
    var(ind_TERPO2) = conv  * REAL( MAX(chem(i,k,j,P_terpo2),0.), KIND=dp)  
    var(ind_TERPOOH) = conv  * REAL( MAX(chem(i,k,j,P_terpooh),0.), KIND=dp)  
    var(ind_CRESOL) = conv  * REAL( MAX(chem(i,k,j,P_CRES),0.), KIND=dp)  
    var(ind_DMS) = conv  * REAL( MAX(chem(i,k,j,P_dms),0.), KIND=dp)  
    var(ind_SO2) = conv  * REAL( MAX(chem(i,k,j,P_so2),0.), KIND=dp)  
    var(ind_SO4) = conv  * REAL( MAX(chem(i,k,j,P_SULF),0.), KIND=dp)  
    var(ind_XO2) = conv  * REAL( MAX(chem(i,k,j,P_xo2),0.), KIND=dp)  
    var(ind_XOH) = conv  * REAL( MAX(chem(i,k,j,P_xoh),0.), KIND=dp)  
    var(ind_XOOH) = conv  * REAL( MAX(chem(i,k,j,P_xooh),0.), KIND=dp)  





  es  = 1000._dp*0.6112_dp*exp(17.67_dp*(t_phy(i,k,j)-273.15_dp)/(t_phy(i,k,j)- 29.65_dp))
  qvs = es / ( p_phy(i,k,j) - es )


  rh =  moist(i,k,j,P_QV) / qvs
  rh = MIN ( MAX ( rh, 0._dp), 1._dp)






   CALL mozcart_Update_Rconst(  &


   rh, &



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








  CALL mozcart_INTEGRATE(TIME_START, TIME_END, &  
          FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK, & 
          ICNTRL_U=icntrl, RCNTRL_U=rcntrl  )






   IF( irr_option(dm) ) THEN
     DO n = param_first_scalar,num_irr_diag
       irr_rates(i,k,j,n) = real( IRR_WRK(irr_diag_ndx(n-1,dm)),kind=4 )
     ENDDO
   ENDIF

    chem(i,k,j,P_o3) = MAX ( REAL (oconv * var(ind_O3), KIND=sp), 0.)  
    chem(i,k,j,P_o) = MAX ( REAL (oconv * var(ind_O), KIND=sp), 0.)  
    chem(i,k,j,P_o1d_cb4) = MAX ( REAL (oconv * var(ind_O1D_CB4), KIND=sp), 0.)  
    chem(i,k,j,P_n2o) = MAX ( REAL (oconv * var(ind_N2O), KIND=sp), 0.)  
    chem(i,k,j,P_no) = MAX ( REAL (oconv * var(ind_NO), KIND=sp), 0.)  
    chem(i,k,j,P_no2) = MAX ( REAL (oconv * var(ind_NO2), KIND=sp), 0.)  
    chem(i,k,j,P_no3) = MAX ( REAL (oconv * var(ind_NO3), KIND=sp), 0.)  
    chem(i,k,j,P_nh3) = MAX ( REAL (oconv * var(ind_NH3), KIND=sp), 0.)  
    chem(i,k,j,P_hno3) = MAX ( REAL (oconv * var(ind_HNO3), KIND=sp), 0.)  
    chem(i,k,j,P_HNO4) = MAX ( REAL (oconv * var(ind_HO2NO2), KIND=sp), 0.)  
    chem(i,k,j,P_n2o5) = MAX ( REAL (oconv * var(ind_N2O5), KIND=sp), 0.)  
    chem(i,k,j,P_h2) = MAX ( REAL (oconv * var(ind_H2), KIND=sp), 0.)  
    chem(i,k,j,P_HO) = MAX ( REAL (oconv * var(ind_OH), KIND=sp), 0.)  
    chem(i,k,j,P_ho2) = MAX ( REAL (oconv * var(ind_HO2), KIND=sp), 0.)  
    chem(i,k,j,P_h2o2) = MAX ( REAL (oconv * var(ind_H2O2), KIND=sp), 0.)  
    chem(i,k,j,P_ch4) = MAX ( REAL (oconv * var(ind_CH4), KIND=sp), 0.)  
    chem(i,k,j,P_co) = MAX ( REAL (oconv * var(ind_CO), KIND=sp), 0.)  
    chem(i,k,j,P_ch3o2) = MAX ( REAL (oconv * var(ind_CH3O2), KIND=sp), 0.)  
    chem(i,k,j,P_ch3ooh) = MAX ( REAL (oconv * var(ind_CH3OOH), KIND=sp), 0.)  
    chem(i,k,j,P_HCHO) = MAX ( REAL (oconv * var(ind_CH2O), KIND=sp), 0.)  
    chem(i,k,j,P_ch3oh) = MAX ( REAL (oconv * var(ind_CH3OH), KIND=sp), 0.)  
    chem(i,k,j,P_c2h4) = MAX ( REAL (oconv * var(ind_C2H4), KIND=sp), 0.)  
    chem(i,k,j,P_eo) = MAX ( REAL (oconv * var(ind_EO), KIND=sp), 0.)  
    chem(i,k,j,P_eo2) = MAX ( REAL (oconv * var(ind_EO2), KIND=sp), 0.)  
    chem(i,k,j,P_ALD) = MAX ( REAL (oconv * var(ind_CH3CHO), KIND=sp), 0.)  
    chem(i,k,j,P_ch3cooh) = MAX ( REAL (oconv * var(ind_CH3COOH), KIND=sp), 0.)  
    chem(i,k,j,P_ACET) = MAX ( REAL (oconv * var(ind_CH3COCH3), KIND=sp), 0.)  
    chem(i,k,j,P_MGLY) = MAX ( REAL (oconv * var(ind_CH3COCHO), KIND=sp), 0.)  
    chem(i,k,j,P_ACO3) = MAX ( REAL (oconv * var(ind_CH3CO3), KIND=sp), 0.)  
    chem(i,k,j,P_PAA) = MAX ( REAL (oconv * var(ind_CH3COOOH), KIND=sp), 0.)  
    chem(i,k,j,P_GLY) = MAX ( REAL (oconv * var(ind_GLYOXAL), KIND=sp), 0.)  
    chem(i,k,j,P_po2) = MAX ( REAL (oconv * var(ind_PO2), KIND=sp), 0.)  
    chem(i,k,j,P_C3H6OOH) = MAX ( REAL (oconv * var(ind_POOH), KIND=sp), 0.)  
    chem(i,k,j,P_pan) = MAX ( REAL (oconv * var(ind_PAN), KIND=sp), 0.)  
    chem(i,k,j,P_mpan) = MAX ( REAL (oconv * var(ind_MPAN), KIND=sp), 0.)  
    chem(i,k,j,P_mco3) = MAX ( REAL (oconv * var(ind_MCO3), KIND=sp), 0.)  
    chem(i,k,j,P_macr) = MAX ( REAL (oconv * var(ind_MACR), KIND=sp), 0.)  
    chem(i,k,j,P_MVKO2) = MAX ( REAL (oconv * var(ind_MACRO2), KIND=sp), 0.)  
    chem(i,k,j,P_MVKOOH) = MAX ( REAL (oconv * var(ind_MACROOH), KIND=sp), 0.)  
    chem(i,k,j,P_mvk) = MAX ( REAL (oconv * var(ind_MVK), KIND=sp), 0.)  
    chem(i,k,j,P_c2h6) = MAX ( REAL (oconv * var(ind_C2H6), KIND=sp), 0.)  
    chem(i,k,j,P_c3h6) = MAX ( REAL (oconv * var(ind_C3H6), KIND=sp), 0.)  
    chem(i,k,j,P_c3h8) = MAX ( REAL (oconv * var(ind_C3H8), KIND=sp), 0.)  
    chem(i,k,j,P_c2h5oh) = MAX ( REAL (oconv * var(ind_C2H5OH), KIND=sp), 0.)  
    chem(i,k,j,P_ETOOH) = MAX ( REAL (oconv * var(ind_C2H5OOH), KIND=sp), 0.)  
    chem(i,k,j,P_PRO2) = MAX ( REAL (oconv * var(ind_C3H7O2), KIND=sp), 0.)  
    chem(i,k,j,P_PROOH) = MAX ( REAL (oconv * var(ind_C3H7OOH), KIND=sp), 0.)  
    chem(i,k,j,P_c10h16) = MAX ( REAL (oconv * var(ind_C10H16), KIND=sp), 0.)  
    chem(i,k,j,P_ACETO2) = MAX ( REAL (oconv * var(ind_RO2), KIND=sp), 0.)  
    chem(i,k,j,P_ACETP) = MAX ( REAL (oconv * var(ind_ROOH), KIND=sp), 0.)  
    chem(i,k,j,P_onit) = MAX ( REAL (oconv * var(ind_ONIT), KIND=sp), 0.)  
    chem(i,k,j,P_onitr) = MAX ( REAL (oconv * var(ind_ONITR), KIND=sp), 0.)  
    chem(i,k,j,P_ISOPR) = MAX ( REAL (oconv * var(ind_ISOP), KIND=sp), 0.)  
    chem(i,k,j,P_ISO2) = MAX ( REAL (oconv * var(ind_ISOPO2), KIND=sp), 0.)  
    chem(i,k,j,P_ISOOH) = MAX ( REAL (oconv * var(ind_ISOPOOH), KIND=sp), 0.)  
    chem(i,k,j,P_ISOPN) = MAX ( REAL (oconv * var(ind_ISOPNO3), KIND=sp), 0.)  
    chem(i,k,j,P_ACETOL) = MAX ( REAL (oconv * var(ind_HYAC), KIND=sp), 0.)  
    chem(i,k,j,P_glyald) = MAX ( REAL (oconv * var(ind_GLYALD), KIND=sp), 0.)  
    chem(i,k,j,P_hydrald) = MAX ( REAL (oconv * var(ind_HYDRALD), KIND=sp), 0.)  
    chem(i,k,j,P_eneo2) = MAX ( REAL (oconv * var(ind_ENEO2), KIND=sp), 0.)  
    chem(i,k,j,P_mek) = MAX ( REAL (oconv * var(ind_MEK), KIND=sp), 0.)  
    chem(i,k,j,P_meko2) = MAX ( REAL (oconv * var(ind_MEKO2), KIND=sp), 0.)  
    chem(i,k,j,P_ETO2) = MAX ( REAL (oconv * var(ind_C2H5O2), KIND=sp), 0.)  
    chem(i,k,j,P_bigene) = MAX ( REAL (oconv * var(ind_BIGENE), KIND=sp), 0.)  
    chem(i,k,j,P_OPEN) = MAX ( REAL (oconv * var(ind_BIGALD), KIND=sp), 0.)  
    chem(i,k,j,P_bigalk) = MAX ( REAL (oconv * var(ind_BIGALK), KIND=sp), 0.)  
    chem(i,k,j,P_alko2) = MAX ( REAL (oconv * var(ind_ALKO2), KIND=sp), 0.)  
    chem(i,k,j,P_alkooh) = MAX ( REAL (oconv * var(ind_ALKOOH), KIND=sp), 0.)  
    chem(i,k,j,P_mekooh) = MAX ( REAL (oconv * var(ind_MEKOOH), KIND=sp), 0.)  
    chem(i,k,j,P_TOL) = MAX ( REAL (oconv * var(ind_TOLUENE), KIND=sp), 0.)  
    chem(i,k,j,P_TO2) = MAX ( REAL (oconv * var(ind_TOLO2), KIND=sp), 0.)  
    chem(i,k,j,P_tolooh) = MAX ( REAL (oconv * var(ind_TOLOOH), KIND=sp), 0.)  
    chem(i,k,j,P_terpo2) = MAX ( REAL (oconv * var(ind_TERPO2), KIND=sp), 0.)  
    chem(i,k,j,P_terpooh) = MAX ( REAL (oconv * var(ind_TERPOOH), KIND=sp), 0.)  
    chem(i,k,j,P_CRES) = MAX ( REAL (oconv * var(ind_CRESOL), KIND=sp), 0.)  
    chem(i,k,j,P_dms) = MAX ( REAL (oconv * var(ind_DMS), KIND=sp), 0.)  
    chem(i,k,j,P_so2) = MAX ( REAL (oconv * var(ind_SO2), KIND=sp), 0.)  
    chem(i,k,j,P_SULF) = MAX ( REAL (oconv * var(ind_SO4), KIND=sp), 0.)  
    chem(i,k,j,P_xo2) = MAX ( REAL (oconv * var(ind_XO2), KIND=sp), 0.)  
    chem(i,k,j,P_xoh) = MAX ( REAL (oconv * var(ind_XOH), KIND=sp), 0.)  
    chem(i,k,j,P_xooh) = MAX ( REAL (oconv * var(ind_XOOH), KIND=sp), 0.)  



    END DO
    END DO
    END DO









END SUBROUTINE  mozcart_interface


END MODULE module_kpp_mozcart_interf 




