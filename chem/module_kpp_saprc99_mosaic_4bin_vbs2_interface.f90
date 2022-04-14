







MODULE module_kpp_saprc99_mosaic_4bin_vbs2_interf 


  USE module_state_description
  USE module_configure

  USE saprc99_mosaic_4bin_vbs2_Parameters
  USE saprc99_mosaic_4bin_vbs2_Precision
  USE saprc99_mosaic_4bin_vbs2_UpdateRconstWRF
  USE saprc99_mosaic_4bin_vbs2_Integrator

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

SUBROUTINE  saprc99_mosaic_4bin_vbs2_interface( &

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
    var(ind_H2O2) = conv  * REAL( MAX(chem(i,k,j,P_h2o2),0.), KIND=dp)  
    var(ind_NO) = conv  * REAL( MAX(chem(i,k,j,P_no),0.), KIND=dp)  
    var(ind_NO2) = conv  * REAL( MAX(chem(i,k,j,P_no2),0.), KIND=dp)  
    var(ind_NO3) = conv  * REAL( MAX(chem(i,k,j,P_no3),0.), KIND=dp)  
    var(ind_N2O5) = conv  * REAL( MAX(chem(i,k,j,P_n2o5),0.), KIND=dp)  
    var(ind_HONO) = conv  * REAL( MAX(chem(i,k,j,P_hono),0.), KIND=dp)  
    var(ind_HNO3) = conv  * REAL( MAX(chem(i,k,j,P_hno3),0.), KIND=dp)  
    var(ind_HNO4) = conv  * REAL( MAX(chem(i,k,j,P_hno4),0.), KIND=dp)  
    var(ind_SO2) = conv  * REAL( MAX(chem(i,k,j,P_so2),0.), KIND=dp)  
    var(ind_H2SO4) = conv  * REAL( MAX(chem(i,k,j,P_h2so4),0.), KIND=dp)  
    var(ind_CO) = conv  * REAL( MAX(chem(i,k,j,P_co),0.), KIND=dp)  
    var(ind_HCHO) = conv  * REAL( MAX(chem(i,k,j,P_hcho),0.), KIND=dp)  
    var(ind_CCHO) = conv  * REAL( MAX(chem(i,k,j,P_ccho),0.), KIND=dp)  
    var(ind_RCHO) = conv  * REAL( MAX(chem(i,k,j,P_rcho),0.), KIND=dp)  
    var(ind_ACET) = conv  * REAL( MAX(chem(i,k,j,P_acet),0.), KIND=dp)  
    var(ind_MEK) = conv  * REAL( MAX(chem(i,k,j,P_mek),0.), KIND=dp)  
    var(ind_HCOOH) = conv  * REAL( MAX(chem(i,k,j,P_hcooh),0.), KIND=dp)  
    var(ind_MEOH) = conv  * REAL( MAX(chem(i,k,j,P_meoh),0.), KIND=dp)  
    var(ind_ETOH) = conv  * REAL( MAX(chem(i,k,j,P_etoh),0.), KIND=dp)  
    var(ind_CCO_OH) = conv  * REAL( MAX(chem(i,k,j,P_cco_oh),0.), KIND=dp)  
    var(ind_RCO_OH) = conv  * REAL( MAX(chem(i,k,j,P_rco_oh),0.), KIND=dp)  
    var(ind_GLY) = conv  * REAL( MAX(chem(i,k,j,P_gly),0.), KIND=dp)  
    var(ind_MGLY) = conv  * REAL( MAX(chem(i,k,j,P_mgly),0.), KIND=dp)  
    var(ind_BACL) = conv  * REAL( MAX(chem(i,k,j,P_bacl),0.), KIND=dp)  
    var(ind_CRES) = conv  * REAL( MAX(chem(i,k,j,P_cres),0.), KIND=dp)  
    var(ind_BALD) = conv  * REAL( MAX(chem(i,k,j,P_bald),0.), KIND=dp)  
    var(ind_ISOPROD) = conv  * REAL( MAX(chem(i,k,j,P_isoprod),0.), KIND=dp)  
    var(ind_METHACRO) = conv  * REAL( MAX(chem(i,k,j,P_methacro),0.), KIND=dp)  
    var(ind_MVK) = conv  * REAL( MAX(chem(i,k,j,P_mvk),0.), KIND=dp)  
    var(ind_PROD2) = conv  * REAL( MAX(chem(i,k,j,P_prod2),0.), KIND=dp)  
    var(ind_DCB1) = conv  * REAL( MAX(chem(i,k,j,P_dcb1),0.), KIND=dp)  
    var(ind_DCB2) = conv  * REAL( MAX(chem(i,k,j,P_dcb2),0.), KIND=dp)  
    var(ind_DCB3) = conv  * REAL( MAX(chem(i,k,j,P_dcb3),0.), KIND=dp)  
    var(ind_ETHENE) = conv  * REAL( MAX(chem(i,k,j,P_ethene),0.), KIND=dp)  
    var(ind_ISOPRENE) = conv  * REAL( MAX(chem(i,k,j,P_isoprene),0.), KIND=dp)  
    var(ind_C2H6) = conv  * REAL( MAX(chem(i,k,j,P_c2h6),0.), KIND=dp)  
    var(ind_C3H8) = conv  * REAL( MAX(chem(i,k,j,P_c3h8),0.), KIND=dp)  
    var(ind_C2H2) = conv  * REAL( MAX(chem(i,k,j,P_c2h2),0.), KIND=dp)  
    var(ind_C3H6) = conv  * REAL( MAX(chem(i,k,j,P_c3h6),0.), KIND=dp)  
    var(ind_ALK3) = conv  * REAL( MAX(chem(i,k,j,P_alk3),0.), KIND=dp)  
    var(ind_ALK4) = conv  * REAL( MAX(chem(i,k,j,P_alk4),0.), KIND=dp)  
    var(ind_ALK5) = conv  * REAL( MAX(chem(i,k,j,P_alk5),0.), KIND=dp)  
    var(ind_ARO1) = conv  * REAL( MAX(chem(i,k,j,P_aro1),0.), KIND=dp)  
    var(ind_ARO2) = conv  * REAL( MAX(chem(i,k,j,P_aro2),0.), KIND=dp)  
    var(ind_OLE1) = conv  * REAL( MAX(chem(i,k,j,P_ole1),0.), KIND=dp)  
    var(ind_OLE2) = conv  * REAL( MAX(chem(i,k,j,P_ole2),0.), KIND=dp)  
    var(ind_TERP) = conv  * REAL( MAX(chem(i,k,j,P_terp),0.), KIND=dp)  
    var(ind_SESQ) = conv  * REAL( MAX(chem(i,k,j,P_sesq),0.), KIND=dp)  
    var(ind_RNO3) = conv  * REAL( MAX(chem(i,k,j,P_rno3),0.), KIND=dp)  
    var(ind_NPHE) = conv  * REAL( MAX(chem(i,k,j,P_nphe),0.), KIND=dp)  
    var(ind_PHEN) = conv  * REAL( MAX(chem(i,k,j,P_phen),0.), KIND=dp)  
    var(ind_PAN) = conv  * REAL( MAX(chem(i,k,j,P_pan),0.), KIND=dp)  
    var(ind_PAN2) = conv  * REAL( MAX(chem(i,k,j,P_pan2),0.), KIND=dp)  
    var(ind_PBZN) = conv  * REAL( MAX(chem(i,k,j,P_pbzn),0.), KIND=dp)  
    var(ind_MA_PAN) = conv  * REAL( MAX(chem(i,k,j,P_ma_pan),0.), KIND=dp)  
    var(ind_CCO_OOH) = conv  * REAL( MAX(chem(i,k,j,P_cco_ooh),0.), KIND=dp)  
    var(ind_RCO_O2) = conv  * REAL( MAX(chem(i,k,j,P_rco_o2),0.), KIND=dp)  
    var(ind_RCO_OOH) = conv  * REAL( MAX(chem(i,k,j,P_rco_ooh),0.), KIND=dp)  
    var(ind_XN) = conv  * REAL( MAX(chem(i,k,j,P_xn),0.), KIND=dp)  
    var(ind_XC) = conv  * REAL( MAX(chem(i,k,j,P_xc),0.), KIND=dp)  
    var(ind_O3P) =  conv * REAL( MAX(o3p(i,k,j),0.), KIND=dp)  
    var(ind_O1D) =  conv * REAL( MAX(o1d(i,k,j),0.), KIND=dp)  
    var(ind_OH) = conv  * REAL( MAX(chem(i,k,j,P_HO),0.), KIND=dp)  
    var(ind_HO2) = conv  * REAL( MAX(chem(i,k,j,P_ho2),0.), KIND=dp)  
    var(ind_C_O2) = conv  * REAL( MAX(chem(i,k,j,P_c_o2),0.), KIND=dp)  
    var(ind_COOH) = conv  * REAL( MAX(chem(i,k,j,P_cooh),0.), KIND=dp)  
    var(ind_ROOH) = conv  * REAL( MAX(chem(i,k,j,P_rooh),0.), KIND=dp)  
    var(ind_RO2_R) = conv  * REAL( MAX(chem(i,k,j,P_ro2_r),0.), KIND=dp)  
    var(ind_R2O2) = conv  * REAL( MAX(chem(i,k,j,P_r2o2),0.), KIND=dp)  
    var(ind_RO2_N) = conv  * REAL( MAX(chem(i,k,j,P_ro2_n),0.), KIND=dp)  
    var(ind_HOCOO) =  conv * REAL( MAX(hocoo(i,k,j),0.), KIND=dp)  
    var(ind_CCO_O2) = conv  * REAL( MAX(chem(i,k,j,P_cco_o2),0.), KIND=dp)  
    var(ind_BZCO_O2) = conv  * REAL( MAX(chem(i,k,j,P_bzco_o2),0.), KIND=dp)  
    var(ind_BZNO2_O) =  conv * REAL( MAX(bzno2_o(i,k,j),0.), KIND=dp)  
    var(ind_BZ_O) =  conv * REAL( MAX(bz_o(i,k,j),0.), KIND=dp)  
    var(ind_MA_RCO3) = conv  * REAL( MAX(chem(i,k,j,P_ma_rco3),0.), KIND=dp)  
    var(ind_TBU_O) =  conv * REAL( MAX(tbu_o(i,k,j),0.), KIND=dp)  
    var(ind_NUME) = conv  * REAL( MAX(chem(i,k,j,P_nume),0.), KIND=dp)  
    var(ind_DEN) = conv  * REAL( MAX(chem(i,k,j,P_den),0.), KIND=dp)  
    var(ind_ANT1_c) = conv  * REAL( MAX(chem(i,k,j,P_ant1_c),0.), KIND=dp)  
    var(ind_ANT1_o) = conv  * REAL( MAX(chem(i,k,j,P_ant1_o),0.), KIND=dp)  
    var(ind_BIOG1_c) = conv  * REAL( MAX(chem(i,k,j,P_biog1_c),0.), KIND=dp)  
    var(ind_BIOG1_o) = conv  * REAL( MAX(chem(i,k,j,P_biog1_o),0.), KIND=dp)  
    var(ind_PSD1) = conv  * REAL( MAX(chem(i,k,j,P_psd1),0.), KIND=dp)  
    var(ind_PSD2) = conv  * REAL( MAX(chem(i,k,j,P_psd2),0.), KIND=dp)  
    var(ind_PCG1_B_C) = conv  * REAL( MAX(chem(i,k,j,P_pcg1_b_c),0.), KIND=dp)  
    var(ind_PCG2_B_C) = conv  * REAL( MAX(chem(i,k,j,P_pcg2_b_c),0.), KIND=dp)  
    var(ind_PCG1_B_O) = conv  * REAL( MAX(chem(i,k,j,P_pcg1_b_o),0.), KIND=dp)  
    var(ind_PCG2_B_O) = conv  * REAL( MAX(chem(i,k,j,P_pcg2_b_o),0.), KIND=dp)  
    var(ind_OPCG1_B_C) = conv  * REAL( MAX(chem(i,k,j,P_opcg1_b_c),0.), KIND=dp)  
    var(ind_OPCG1_B_O) = conv  * REAL( MAX(chem(i,k,j,P_opcg1_b_o),0.), KIND=dp)  
    var(ind_PCG1_F_C) = conv  * REAL( MAX(chem(i,k,j,P_pcg1_f_c),0.), KIND=dp)  
    var(ind_PCG2_F_C) = conv  * REAL( MAX(chem(i,k,j,P_pcg2_f_c),0.), KIND=dp)  
    var(ind_PCG1_F_O) = conv  * REAL( MAX(chem(i,k,j,P_pcg1_f_o),0.), KIND=dp)  
    var(ind_PCG2_F_O) = conv  * REAL( MAX(chem(i,k,j,P_pcg2_f_o),0.), KIND=dp)  
    var(ind_OPCG1_F_C) = conv  * REAL( MAX(chem(i,k,j,P_opcg1_f_c),0.), KIND=dp)  
    var(ind_OPCG1_F_O) = conv  * REAL( MAX(chem(i,k,j,P_opcg1_f_o),0.), KIND=dp)  
    var(ind_CH4) = conv  * REAL( MAX(chem(i,k,j,P_ch4),0.), KIND=dp)  







   CALL saprc99_mosaic_4bin_vbs2_Update_Rconst(  &

   var(ind_NUME),var(ind_DEN), &



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








  CALL saprc99_mosaic_4bin_vbs2_INTEGRATE(TIME_START, TIME_END, &  
          FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK, & 
          ICNTRL_U=icntrl, RCNTRL_U=rcntrl  )







    chem(i,k,j,P_o3) = MAX ( REAL (oconv * var(ind_O3), KIND=sp), 0.)  
    chem(i,k,j,P_h2o2) = MAX ( REAL (oconv * var(ind_H2O2), KIND=sp), 0.)  
    chem(i,k,j,P_no) = MAX ( REAL (oconv * var(ind_NO), KIND=sp), 0.)  
    chem(i,k,j,P_no2) = MAX ( REAL (oconv * var(ind_NO2), KIND=sp), 0.)  
    chem(i,k,j,P_no3) = MAX ( REAL (oconv * var(ind_NO3), KIND=sp), 0.)  
    chem(i,k,j,P_n2o5) = MAX ( REAL (oconv * var(ind_N2O5), KIND=sp), 0.)  
    chem(i,k,j,P_hono) = MAX ( REAL (oconv * var(ind_HONO), KIND=sp), 0.)  
    chem(i,k,j,P_hno3) = MAX ( REAL (oconv * var(ind_HNO3), KIND=sp), 0.)  
    chem(i,k,j,P_hno4) = MAX ( REAL (oconv * var(ind_HNO4), KIND=sp), 0.)  
    chem(i,k,j,P_so2) = MAX ( REAL (oconv * var(ind_SO2), KIND=sp), 0.)  
    chem(i,k,j,P_h2so4) = MAX ( REAL (oconv * var(ind_H2SO4), KIND=sp), 0.)  
    chem(i,k,j,P_co) = MAX ( REAL (oconv * var(ind_CO), KIND=sp), 0.)  
    chem(i,k,j,P_hcho) = MAX ( REAL (oconv * var(ind_HCHO), KIND=sp), 0.)  
    chem(i,k,j,P_ccho) = MAX ( REAL (oconv * var(ind_CCHO), KIND=sp), 0.)  
    chem(i,k,j,P_rcho) = MAX ( REAL (oconv * var(ind_RCHO), KIND=sp), 0.)  
    chem(i,k,j,P_acet) = MAX ( REAL (oconv * var(ind_ACET), KIND=sp), 0.)  
    chem(i,k,j,P_mek) = MAX ( REAL (oconv * var(ind_MEK), KIND=sp), 0.)  
    chem(i,k,j,P_hcooh) = MAX ( REAL (oconv * var(ind_HCOOH), KIND=sp), 0.)  
    chem(i,k,j,P_meoh) = MAX ( REAL (oconv * var(ind_MEOH), KIND=sp), 0.)  
    chem(i,k,j,P_etoh) = MAX ( REAL (oconv * var(ind_ETOH), KIND=sp), 0.)  
    chem(i,k,j,P_cco_oh) = MAX ( REAL (oconv * var(ind_CCO_OH), KIND=sp), 0.)  
    chem(i,k,j,P_rco_oh) = MAX ( REAL (oconv * var(ind_RCO_OH), KIND=sp), 0.)  
    chem(i,k,j,P_gly) = MAX ( REAL (oconv * var(ind_GLY), KIND=sp), 0.)  
    chem(i,k,j,P_mgly) = MAX ( REAL (oconv * var(ind_MGLY), KIND=sp), 0.)  
    chem(i,k,j,P_bacl) = MAX ( REAL (oconv * var(ind_BACL), KIND=sp), 0.)  
    chem(i,k,j,P_cres) = MAX ( REAL (oconv * var(ind_CRES), KIND=sp), 0.)  
    chem(i,k,j,P_bald) = MAX ( REAL (oconv * var(ind_BALD), KIND=sp), 0.)  
    chem(i,k,j,P_isoprod) = MAX ( REAL (oconv * var(ind_ISOPROD), KIND=sp), 0.)  
    chem(i,k,j,P_methacro) = MAX ( REAL (oconv * var(ind_METHACRO), KIND=sp), 0.)  
    chem(i,k,j,P_mvk) = MAX ( REAL (oconv * var(ind_MVK), KIND=sp), 0.)  
    chem(i,k,j,P_prod2) = MAX ( REAL (oconv * var(ind_PROD2), KIND=sp), 0.)  
    chem(i,k,j,P_dcb1) = MAX ( REAL (oconv * var(ind_DCB1), KIND=sp), 0.)  
    chem(i,k,j,P_dcb2) = MAX ( REAL (oconv * var(ind_DCB2), KIND=sp), 0.)  
    chem(i,k,j,P_dcb3) = MAX ( REAL (oconv * var(ind_DCB3), KIND=sp), 0.)  
    chem(i,k,j,P_ethene) = MAX ( REAL (oconv * var(ind_ETHENE), KIND=sp), 0.)  
    chem(i,k,j,P_isoprene) = MAX ( REAL (oconv * var(ind_ISOPRENE), KIND=sp), 0.)  
    chem(i,k,j,P_c2h6) = MAX ( REAL (oconv * var(ind_C2H6), KIND=sp), 0.)  
    chem(i,k,j,P_c3h8) = MAX ( REAL (oconv * var(ind_C3H8), KIND=sp), 0.)  
    chem(i,k,j,P_c2h2) = MAX ( REAL (oconv * var(ind_C2H2), KIND=sp), 0.)  
    chem(i,k,j,P_c3h6) = MAX ( REAL (oconv * var(ind_C3H6), KIND=sp), 0.)  
    chem(i,k,j,P_alk3) = MAX ( REAL (oconv * var(ind_ALK3), KIND=sp), 0.)  
    chem(i,k,j,P_alk4) = MAX ( REAL (oconv * var(ind_ALK4), KIND=sp), 0.)  
    chem(i,k,j,P_alk5) = MAX ( REAL (oconv * var(ind_ALK5), KIND=sp), 0.)  
    chem(i,k,j,P_aro1) = MAX ( REAL (oconv * var(ind_ARO1), KIND=sp), 0.)  
    chem(i,k,j,P_aro2) = MAX ( REAL (oconv * var(ind_ARO2), KIND=sp), 0.)  
    chem(i,k,j,P_ole1) = MAX ( REAL (oconv * var(ind_OLE1), KIND=sp), 0.)  
    chem(i,k,j,P_ole2) = MAX ( REAL (oconv * var(ind_OLE2), KIND=sp), 0.)  
    chem(i,k,j,P_terp) = MAX ( REAL (oconv * var(ind_TERP), KIND=sp), 0.)  
    chem(i,k,j,P_sesq) = MAX ( REAL (oconv * var(ind_SESQ), KIND=sp), 0.)  
    chem(i,k,j,P_rno3) = MAX ( REAL (oconv * var(ind_RNO3), KIND=sp), 0.)  
    chem(i,k,j,P_nphe) = MAX ( REAL (oconv * var(ind_NPHE), KIND=sp), 0.)  
    chem(i,k,j,P_phen) = MAX ( REAL (oconv * var(ind_PHEN), KIND=sp), 0.)  
    chem(i,k,j,P_pan) = MAX ( REAL (oconv * var(ind_PAN), KIND=sp), 0.)  
    chem(i,k,j,P_pan2) = MAX ( REAL (oconv * var(ind_PAN2), KIND=sp), 0.)  
    chem(i,k,j,P_pbzn) = MAX ( REAL (oconv * var(ind_PBZN), KIND=sp), 0.)  
    chem(i,k,j,P_ma_pan) = MAX ( REAL (oconv * var(ind_MA_PAN), KIND=sp), 0.)  
    chem(i,k,j,P_cco_ooh) = MAX ( REAL (oconv * var(ind_CCO_OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rco_o2) = MAX ( REAL (oconv * var(ind_RCO_O2), KIND=sp), 0.)  
    chem(i,k,j,P_rco_ooh) = MAX ( REAL (oconv * var(ind_RCO_OOH), KIND=sp), 0.)  
    chem(i,k,j,P_xn) = MAX ( REAL (oconv * var(ind_XN), KIND=sp), 0.)  
    chem(i,k,j,P_xc) = MAX ( REAL (oconv * var(ind_XC), KIND=sp), 0.)  
   o3p(i,k,j) = MAX (REAL (oconv * var(ind_O3P) , KIND=sp),0.) 
   o1d(i,k,j) = MAX (REAL (oconv * var(ind_O1D) , KIND=sp),0.) 
    chem(i,k,j,P_HO) = MAX ( REAL (oconv * var(ind_OH), KIND=sp), 0.)  
    chem(i,k,j,P_ho2) = MAX ( REAL (oconv * var(ind_HO2), KIND=sp), 0.)  
    chem(i,k,j,P_c_o2) = MAX ( REAL (oconv * var(ind_C_O2), KIND=sp), 0.)  
    chem(i,k,j,P_cooh) = MAX ( REAL (oconv * var(ind_COOH), KIND=sp), 0.)  
    chem(i,k,j,P_rooh) = MAX ( REAL (oconv * var(ind_ROOH), KIND=sp), 0.)  
    chem(i,k,j,P_ro2_r) = MAX ( REAL (oconv * var(ind_RO2_R), KIND=sp), 0.)  
    chem(i,k,j,P_r2o2) = MAX ( REAL (oconv * var(ind_R2O2), KIND=sp), 0.)  
    chem(i,k,j,P_ro2_n) = MAX ( REAL (oconv * var(ind_RO2_N), KIND=sp), 0.)  
   hocoo(i,k,j) = MAX (REAL (oconv * var(ind_HOCOO) , KIND=sp),0.) 
    chem(i,k,j,P_cco_o2) = MAX ( REAL (oconv * var(ind_CCO_O2), KIND=sp), 0.)  
    chem(i,k,j,P_bzco_o2) = MAX ( REAL (oconv * var(ind_BZCO_O2), KIND=sp), 0.)  
   bzno2_o(i,k,j) = MAX (REAL (oconv * var(ind_BZNO2_O) , KIND=sp),0.) 
   bz_o(i,k,j) = MAX (REAL (oconv * var(ind_BZ_O) , KIND=sp),0.) 
    chem(i,k,j,P_ma_rco3) = MAX ( REAL (oconv * var(ind_MA_RCO3), KIND=sp), 0.)  
   tbu_o(i,k,j) = MAX (REAL (oconv * var(ind_TBU_O) , KIND=sp),0.) 
    chem(i,k,j,P_nume) = MAX ( REAL (oconv * var(ind_NUME), KIND=sp), 0.)  
    chem(i,k,j,P_den) = MAX ( REAL (oconv * var(ind_DEN), KIND=sp), 0.)  
    chem(i,k,j,P_ant1_c) = MAX ( REAL (oconv * var(ind_ANT1_c), KIND=sp), 0.)  
    chem(i,k,j,P_ant1_o) = MAX ( REAL (oconv * var(ind_ANT1_o), KIND=sp), 0.)  
    chem(i,k,j,P_biog1_c) = MAX ( REAL (oconv * var(ind_BIOG1_c), KIND=sp), 0.)  
    chem(i,k,j,P_biog1_o) = MAX ( REAL (oconv * var(ind_BIOG1_o), KIND=sp), 0.)  
    chem(i,k,j,P_psd1) = MAX ( REAL (oconv * var(ind_PSD1), KIND=sp), 0.)  
    chem(i,k,j,P_psd2) = MAX ( REAL (oconv * var(ind_PSD2), KIND=sp), 0.)  
    chem(i,k,j,P_pcg1_b_c) = MAX ( REAL (oconv * var(ind_PCG1_B_C), KIND=sp), 0.)  
    chem(i,k,j,P_pcg2_b_c) = MAX ( REAL (oconv * var(ind_PCG2_B_C), KIND=sp), 0.)  
    chem(i,k,j,P_pcg1_b_o) = MAX ( REAL (oconv * var(ind_PCG1_B_O), KIND=sp), 0.)  
    chem(i,k,j,P_pcg2_b_o) = MAX ( REAL (oconv * var(ind_PCG2_B_O), KIND=sp), 0.)  
    chem(i,k,j,P_opcg1_b_c) = MAX ( REAL (oconv * var(ind_OPCG1_B_C), KIND=sp), 0.)  
    chem(i,k,j,P_opcg1_b_o) = MAX ( REAL (oconv * var(ind_OPCG1_B_O), KIND=sp), 0.)  
    chem(i,k,j,P_pcg1_f_c) = MAX ( REAL (oconv * var(ind_PCG1_F_C), KIND=sp), 0.)  
    chem(i,k,j,P_pcg2_f_c) = MAX ( REAL (oconv * var(ind_PCG2_F_C), KIND=sp), 0.)  
    chem(i,k,j,P_pcg1_f_o) = MAX ( REAL (oconv * var(ind_PCG1_F_O), KIND=sp), 0.)  
    chem(i,k,j,P_pcg2_f_o) = MAX ( REAL (oconv * var(ind_PCG2_F_O), KIND=sp), 0.)  
    chem(i,k,j,P_opcg1_f_c) = MAX ( REAL (oconv * var(ind_OPCG1_F_C), KIND=sp), 0.)  
    chem(i,k,j,P_opcg1_f_o) = MAX ( REAL (oconv * var(ind_OPCG1_F_O), KIND=sp), 0.)  
    chem(i,k,j,P_ch4) = MAX ( REAL (oconv * var(ind_CH4), KIND=sp), 0.)  



    END DO
    END DO
    END DO









END SUBROUTINE  saprc99_mosaic_4bin_vbs2_interface


END MODULE module_kpp_saprc99_mosaic_4bin_vbs2_interf 




