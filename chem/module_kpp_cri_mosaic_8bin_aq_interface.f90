







MODULE module_kpp_cri_mosaic_8bin_aq_interf 


  USE module_state_description
  USE module_configure

  USE cri_mosaic_8bin_aq_Parameters
  USE cri_mosaic_8bin_aq_Precision
  USE cri_mosaic_8bin_aq_UpdateRconstWRF
  USE cri_mosaic_8bin_aq_Integrator

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

SUBROUTINE  cri_mosaic_8bin_aq_interface( &

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
 



 

  REAL( KIND = dp ) :: RO2




 


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
 


    var(ind_HONO) = conv  * REAL( MAX(chem(i,k,j,P_hono),0.), KIND=dp)  
    var(ind_O3) = conv  * REAL( MAX(chem(i,k,j,P_o3),0.), KIND=dp)  
    var(ind_HCHO) = conv  * REAL( MAX(chem(i,k,j,P_hcho),0.), KIND=dp)  
    var(ind_PAN) = conv  * REAL( MAX(chem(i,k,j,P_pan),0.), KIND=dp)  
    var(ind_C2H4) = conv  * REAL( MAX(chem(i,k,j,P_c2h4),0.), KIND=dp)  
    var(ind_CO) = conv  * REAL( MAX(chem(i,k,j,P_co),0.), KIND=dp)  
    var(ind_HNO3) = conv  * REAL( MAX(chem(i,k,j,P_hno3),0.), KIND=dp)  
    var(ind_N2O5) = conv  * REAL( MAX(chem(i,k,j,P_n2o5),0.), KIND=dp)  
    var(ind_HNO4) = conv  * REAL( MAX(chem(i,k,j,P_hno4),0.), KIND=dp)  
    var(ind_NO3) = conv  * REAL( MAX(chem(i,k,j,P_no3),0.), KIND=dp)  
    var(ind_O1D) =  conv * REAL( MAX(o1d(i,k,j),0.), KIND=dp)  
    var(ind_O3P) =  conv * REAL( MAX(o3p(i,k,j),0.), KIND=dp)  
    var(ind_OH) = conv  * REAL( MAX(chem(i,k,j,P_HO),0.), KIND=dp)  
    var(ind_HO2) = conv  * REAL( MAX(chem(i,k,j,P_ho2),0.), KIND=dp)  
    var(ind_H2O2) = conv  * REAL( MAX(chem(i,k,j,P_h2o2),0.), KIND=dp)  
    var(ind_C2H6) = conv  * REAL( MAX(chem(i,k,j,P_c2h6),0.), KIND=dp)  
    var(ind_HCOOH) = conv  * REAL( MAX(chem(i,k,j,P_hcooh),0.), KIND=dp)  
    var(ind_CH3CO3) = conv  * REAL( MAX(chem(i,k,j,P_ACO3),0.), KIND=dp)  
    var(ind_CH3OO) = conv  * REAL( MAX(chem(i,k,j,P_ch3oo),0.), KIND=dp)  
    var(ind_C2H5O2) = conv  * REAL( MAX(chem(i,k,j,P_c2h5o2),0.), KIND=dp)  
    var(ind_HSO3) = conv  * REAL( MAX(chem(i,k,j,P_hso3),0.), KIND=dp)  
    var(ind_SO3) = conv  * REAL( MAX(chem(i,k,j,P_so3),0.), KIND=dp)  
    var(ind_SO2) = conv  * REAL( MAX(chem(i,k,j,P_so2),0.), KIND=dp)  
    var(ind_NO2) = conv  * REAL( MAX(chem(i,k,j,P_no2),0.), KIND=dp)  
    var(ind_NO) = conv  * REAL( MAX(chem(i,k,j,P_no),0.), KIND=dp)  
    var(ind_C3H8) = conv  * REAL( MAX(chem(i,k,j,P_c3h8),0.), KIND=dp)  
    var(ind_NC4H10) = conv  * REAL( MAX(chem(i,k,j,P_nc4h10),0.), KIND=dp)  
    var(ind_HOCH2CH2O2) = conv  * REAL( MAX(chem(i,k,j,P_hoch2ch2o2),0.), KIND=dp)  
    var(ind_IC3H7O2) = conv  * REAL( MAX(chem(i,k,j,P_ic3h7o2),0.), KIND=dp)  
    var(ind_C5H8) = conv  * REAL( MAX(chem(i,k,j,P_c5h8),0.), KIND=dp)  
    var(ind_BENZENE) = conv  * REAL( MAX(chem(i,k,j,P_benzene),0.), KIND=dp)  
    var(ind_TOLUENE) = conv  * REAL( MAX(chem(i,k,j,P_toluene),0.), KIND=dp)  
    var(ind_OXYL) = conv  * REAL( MAX(chem(i,k,j,P_oxyl),0.), KIND=dp)  
    var(ind_NPROPOL) = conv  * REAL( MAX(chem(i,k,j,P_npropol),0.), KIND=dp)  
    var(ind_C2H2) = conv  * REAL( MAX(chem(i,k,j,P_c2h2),0.), KIND=dp)  
    var(ind_C3H6) = conv  * REAL( MAX(chem(i,k,j,P_c3h6),0.), KIND=dp)  
    var(ind_TBUT2ENE) = conv  * REAL( MAX(chem(i,k,j,P_tbut2ene),0.), KIND=dp)  
    var(ind_CH3CHO) = conv  * REAL( MAX(chem(i,k,j,P_ch3cho),0.), KIND=dp)  
    var(ind_C2H5CHO) = conv  * REAL( MAX(chem(i,k,j,P_c2h5cho),0.), KIND=dp)  
    var(ind_CH3CO2H) = conv  * REAL( MAX(chem(i,k,j,P_ch3co2h),0.), KIND=dp)  
    var(ind_CH3COCH3) = conv  * REAL( MAX(chem(i,k,j,P_KET),0.), KIND=dp)  
    var(ind_MEK) = conv  * REAL( MAX(chem(i,k,j,P_mek),0.), KIND=dp)  
    var(ind_CH3OH) = conv  * REAL( MAX(chem(i,k,j,P_ch3oh),0.), KIND=dp)  
    var(ind_C2H5OH) = conv  * REAL( MAX(chem(i,k,j,P_c2h5oh),0.), KIND=dp)  
    var(ind_IC3H7NO3) = conv  * REAL( MAX(chem(i,k,j,P_ic3h7no3),0.), KIND=dp)  
    var(ind_IPROPOL) = conv  * REAL( MAX(chem(i,k,j,P_ipropol),0.), KIND=dp)  
    var(ind_CH3NO3) = conv  * REAL( MAX(chem(i,k,j,P_ch3no3),0.), KIND=dp)  
    var(ind_C2H5NO3) = conv  * REAL( MAX(chem(i,k,j,P_c2h5no3),0.), KIND=dp)  
    var(ind_HOC2H4NO3) = conv  * REAL( MAX(chem(i,k,j,P_hoc2h4no3),0.), KIND=dp)  
    var(ind_CH3OOH) = conv  * REAL( MAX(chem(i,k,j,P_ch3ooh),0.), KIND=dp)  
    var(ind_C2H5OOH) = conv  * REAL( MAX(chem(i,k,j,P_c2h5ooh),0.), KIND=dp)  
    var(ind_IC3H7OOH) = conv  * REAL( MAX(chem(i,k,j,P_PROOH),0.), KIND=dp)  
    var(ind_CH3CO3H) = conv  * REAL( MAX(chem(i,k,j,P_PAA),0.), KIND=dp)  
    var(ind_HOC2H4OOH) = conv  * REAL( MAX(chem(i,k,j,P_hoc2h4ooh),0.), KIND=dp)  
    var(ind_RN10O2) = conv  * REAL( MAX(chem(i,k,j,P_rn10o2),0.), KIND=dp)  
    var(ind_RN13O2) = conv  * REAL( MAX(chem(i,k,j,P_rn13o2),0.), KIND=dp)  
    var(ind_RN16O2) = conv  * REAL( MAX(chem(i,k,j,P_rn16o2),0.), KIND=dp)  
    var(ind_RN19O2) = conv  * REAL( MAX(chem(i,k,j,P_rn19o2),0.), KIND=dp)  
    var(ind_RN9O2) = conv  * REAL( MAX(chem(i,k,j,P_rn9o2),0.), KIND=dp)  
    var(ind_RN12O2) = conv  * REAL( MAX(chem(i,k,j,P_rn12o2),0.), KIND=dp)  
    var(ind_RN15O2) = conv  * REAL( MAX(chem(i,k,j,P_rn15o2),0.), KIND=dp)  
    var(ind_RN18O2) = conv  * REAL( MAX(chem(i,k,j,P_rn18o2),0.), KIND=dp)  
    var(ind_NRN6O2) = conv  * REAL( MAX(chem(i,k,j,P_nrn6o2),0.), KIND=dp)  
    var(ind_NRN9O2) = conv  * REAL( MAX(chem(i,k,j,P_nrn9o2),0.), KIND=dp)  
    var(ind_NRN12O2) = conv  * REAL( MAX(chem(i,k,j,P_nrn12o2),0.), KIND=dp)  
    var(ind_CARB14) = conv  * REAL( MAX(chem(i,k,j,P_carb14),0.), KIND=dp)  
    var(ind_RN11O2) = conv  * REAL( MAX(chem(i,k,j,P_rn11o2),0.), KIND=dp)  
    var(ind_RN14O2) = conv  * REAL( MAX(chem(i,k,j,P_rn14o2),0.), KIND=dp)  
    var(ind_CARB17) = conv  * REAL( MAX(chem(i,k,j,P_carb17),0.), KIND=dp)  
    var(ind_RN8O2) = conv  * REAL( MAX(chem(i,k,j,P_rn8o2),0.), KIND=dp)  
    var(ind_RN17O2) = conv  * REAL( MAX(chem(i,k,j,P_rn17o2),0.), KIND=dp)  
    var(ind_RN10NO3) = conv  * REAL( MAX(chem(i,k,j,P_rn10no3),0.), KIND=dp)  
    var(ind_RN13NO3) = conv  * REAL( MAX(chem(i,k,j,P_rn13no3),0.), KIND=dp)  
    var(ind_RN19NO3) = conv  * REAL( MAX(chem(i,k,j,P_rn19no3),0.), KIND=dp)  
    var(ind_RN9NO3) = conv  * REAL( MAX(chem(i,k,j,P_rn9no3),0.), KIND=dp)  
    var(ind_RN12NO3) = conv  * REAL( MAX(chem(i,k,j,P_rn12no3),0.), KIND=dp)  
    var(ind_RN15NO3) = conv  * REAL( MAX(chem(i,k,j,P_rn15no3),0.), KIND=dp)  
    var(ind_RN18NO3) = conv  * REAL( MAX(chem(i,k,j,P_rn18no3),0.), KIND=dp)  
    var(ind_RN16NO3) = conv  * REAL( MAX(chem(i,k,j,P_rn16no3),0.), KIND=dp)  
    var(ind_RN10OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn10ooh),0.), KIND=dp)  
    var(ind_RN13OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn13ooh),0.), KIND=dp)  
    var(ind_RN16OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn16ooh),0.), KIND=dp)  
    var(ind_RN19OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn19ooh),0.), KIND=dp)  
    var(ind_RN8OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn8ooh),0.), KIND=dp)  
    var(ind_RN11OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn11ooh),0.), KIND=dp)  
    var(ind_RN14OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn14ooh),0.), KIND=dp)  
    var(ind_RN17OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn17ooh),0.), KIND=dp)  
    var(ind_RN9OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn9ooh),0.), KIND=dp)  
    var(ind_RN12OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn12ooh),0.), KIND=dp)  
    var(ind_RN15OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn15ooh),0.), KIND=dp)  
    var(ind_RN18OOH) = conv  * REAL( MAX(chem(i,k,j,P_rn18ooh),0.), KIND=dp)  
    var(ind_NRN6OOH) = conv  * REAL( MAX(chem(i,k,j,P_nrn6ooh),0.), KIND=dp)  
    var(ind_NRN9OOH) = conv  * REAL( MAX(chem(i,k,j,P_nrn9ooh),0.), KIND=dp)  
    var(ind_NRN12OOH) = conv  * REAL( MAX(chem(i,k,j,P_nrn12ooh),0.), KIND=dp)  
    var(ind_APINENE) = conv  * REAL( MAX(chem(i,k,j,P_apinene),0.), KIND=dp)  
    var(ind_BPINENE) = conv  * REAL( MAX(chem(i,k,j,P_bpinene),0.), KIND=dp)  
    var(ind_RN13AO2) = conv  * REAL( MAX(chem(i,k,j,P_rn13ao2),0.), KIND=dp)  
    var(ind_RN16AO2) = conv  * REAL( MAX(chem(i,k,j,P_rn16ao2),0.), KIND=dp)  
    var(ind_RN15AO2) = conv  * REAL( MAX(chem(i,k,j,P_rn15ao2),0.), KIND=dp)  
    var(ind_RN18AO2) = conv  * REAL( MAX(chem(i,k,j,P_rn18ao2),0.), KIND=dp)  
    var(ind_CARB7) = conv  * REAL( MAX(chem(i,k,j,P_carb7),0.), KIND=dp)  
    var(ind_CARB10) = conv  * REAL( MAX(chem(i,k,j,P_carb10),0.), KIND=dp)  
    var(ind_CARB13) = conv  * REAL( MAX(chem(i,k,j,P_carb13),0.), KIND=dp)  
    var(ind_CARB16) = conv  * REAL( MAX(chem(i,k,j,P_carb16),0.), KIND=dp)  
    var(ind_CARB3) = conv  * REAL( MAX(chem(i,k,j,P_carb3),0.), KIND=dp)  
    var(ind_CARB6) = conv  * REAL( MAX(chem(i,k,j,P_carb6),0.), KIND=dp)  
    var(ind_CARB9) = conv  * REAL( MAX(chem(i,k,j,P_carb9),0.), KIND=dp)  
    var(ind_CARB12) = conv  * REAL( MAX(chem(i,k,j,P_carb12),0.), KIND=dp)  
    var(ind_CARB15) = conv  * REAL( MAX(chem(i,k,j,P_carb15),0.), KIND=dp)  
    var(ind_C2H5CO3H) = conv  * REAL( MAX(chem(i,k,j,P_c2h5co3h),0.), KIND=dp)  
    var(ind_C2H5CO3) = conv  * REAL( MAX(chem(i,k,j,P_c2h5co3),0.), KIND=dp)  
    var(ind_PPN) = conv  * REAL( MAX(chem(i,k,j,P_ppn),0.), KIND=dp)  
    var(ind_HOCH2CHO) = conv  * REAL( MAX(chem(i,k,j,P_hoch2cho),0.), KIND=dp)  
    var(ind_HOCH2CO3) = conv  * REAL( MAX(chem(i,k,j,P_hoch2co3),0.), KIND=dp)  
    var(ind_HOCH2CO3H) = conv  * REAL( MAX(chem(i,k,j,P_hoch2co3h),0.), KIND=dp)  
    var(ind_PHAN) = conv  * REAL( MAX(chem(i,k,j,P_phan),0.), KIND=dp)  
    var(ind_CCARB12) = conv  * REAL( MAX(chem(i,k,j,P_ccarb12),0.), KIND=dp)  
    var(ind_RU14O2) = conv  * REAL( MAX(chem(i,k,j,P_ru14o2),0.), KIND=dp)  
    var(ind_RU12O2) = conv  * REAL( MAX(chem(i,k,j,P_ru12o2),0.), KIND=dp)  
    var(ind_CH3CL) = conv  * REAL( MAX(chem(i,k,j,P_ch3cl),0.), KIND=dp)  
    var(ind_CH2CL2) = conv  * REAL( MAX(chem(i,k,j,P_ch2cl2),0.), KIND=dp)  
    var(ind_CHCL3) = conv  * REAL( MAX(chem(i,k,j,P_chcl3),0.), KIND=dp)  
    var(ind_CH3CCL3) = conv  * REAL( MAX(chem(i,k,j,P_ch3ccl3),0.), KIND=dp)  
    var(ind_CDICLETH) = conv  * REAL( MAX(chem(i,k,j,P_cdicleth),0.), KIND=dp)  
    var(ind_TDICLETH) = conv  * REAL( MAX(chem(i,k,j,P_tdicleth),0.), KIND=dp)  
    var(ind_TRICLETH) = conv  * REAL( MAX(chem(i,k,j,P_tricleth),0.), KIND=dp)  
    var(ind_TCE) = conv  * REAL( MAX(chem(i,k,j,P_tce),0.), KIND=dp)  
    var(ind_RU10O2) = conv  * REAL( MAX(chem(i,k,j,P_ru10o2),0.), KIND=dp)  
    var(ind_UCARB12) = conv  * REAL( MAX(chem(i,k,j,P_ucarb12),0.), KIND=dp)  
    var(ind_UCARB10) = conv  * REAL( MAX(chem(i,k,j,P_ucarb10),0.), KIND=dp)  
    var(ind_RU14NO3) = conv  * REAL( MAX(chem(i,k,j,P_ru14no3),0.), KIND=dp)  
    var(ind_RU14OOH) = conv  * REAL( MAX(chem(i,k,j,P_ru14ooh),0.), KIND=dp)  
    var(ind_RU12OOH) = conv  * REAL( MAX(chem(i,k,j,P_ru12ooh),0.), KIND=dp)  
    var(ind_RU10OOH) = conv  * REAL( MAX(chem(i,k,j,P_ru10ooh),0.), KIND=dp)  
    var(ind_MPAN) = conv  * REAL( MAX(chem(i,k,j,P_mpan),0.), KIND=dp)  
    var(ind_RU12PAN) = conv  * REAL( MAX(chem(i,k,j,P_ru12pan),0.), KIND=dp)  
    var(ind_NRU14O2) = conv  * REAL( MAX(chem(i,k,j,P_nru14o2),0.), KIND=dp)  
    var(ind_NUCARB12) = conv  * REAL( MAX(chem(i,k,j,P_nucarb12),0.), KIND=dp)  
    var(ind_NRU14OOH) = conv  * REAL( MAX(chem(i,k,j,P_nru14ooh),0.), KIND=dp)  
    var(ind_NRU12O2) = conv  * REAL( MAX(chem(i,k,j,P_nru12o2),0.), KIND=dp)  
    var(ind_NRU12OOH) = conv  * REAL( MAX(chem(i,k,j,P_nru12ooh),0.), KIND=dp)  
    var(ind_NOA) = conv  * REAL( MAX(chem(i,k,j,P_noa),0.), KIND=dp)  
    var(ind_RA13O2) = conv  * REAL( MAX(chem(i,k,j,P_ra13o2),0.), KIND=dp)  
    var(ind_RA13NO3) = conv  * REAL( MAX(chem(i,k,j,P_ra13no3),0.), KIND=dp)  
    var(ind_RA13OOH) = conv  * REAL( MAX(chem(i,k,j,P_ra13ooh),0.), KIND=dp)  
    var(ind_UDCARB8) = conv  * REAL( MAX(chem(i,k,j,P_udcarb8),0.), KIND=dp)  
    var(ind_AROH14) = conv  * REAL( MAX(chem(i,k,j,P_aroh14),0.), KIND=dp)  
    var(ind_RAROH14) = conv  * REAL( MAX(chem(i,k,j,P_raroh14),0.), KIND=dp)  
    var(ind_ARNOH14) = conv  * REAL( MAX(chem(i,k,j,P_arnoh14),0.), KIND=dp)  
    var(ind_RA16O2) = conv  * REAL( MAX(chem(i,k,j,P_ra16o2),0.), KIND=dp)  
    var(ind_RA16NO3) = conv  * REAL( MAX(chem(i,k,j,P_ra16no3),0.), KIND=dp)  
    var(ind_RA16OOH) = conv  * REAL( MAX(chem(i,k,j,P_ra16ooh),0.), KIND=dp)  
    var(ind_UDCARB11) = conv  * REAL( MAX(chem(i,k,j,P_udcarb11),0.), KIND=dp)  
    var(ind_AROH17) = conv  * REAL( MAX(chem(i,k,j,P_aroh17),0.), KIND=dp)  
    var(ind_RAROH17) = conv  * REAL( MAX(chem(i,k,j,P_raroh17),0.), KIND=dp)  
    var(ind_ARNOH17) = conv  * REAL( MAX(chem(i,k,j,P_arnoh17),0.), KIND=dp)  
    var(ind_UDCARB14) = conv  * REAL( MAX(chem(i,k,j,P_udcarb14),0.), KIND=dp)  
    var(ind_RA19AO2) = conv  * REAL( MAX(chem(i,k,j,P_ra19ao2),0.), KIND=dp)  
    var(ind_RA19CO2) = conv  * REAL( MAX(chem(i,k,j,P_ra19co2),0.), KIND=dp)  
    var(ind_RA19NO3) = conv  * REAL( MAX(chem(i,k,j,P_ra19no3),0.), KIND=dp)  
    var(ind_RA19OOH) = conv  * REAL( MAX(chem(i,k,j,P_ra19ooh),0.), KIND=dp)  
    var(ind_RTN28O2) = conv  * REAL( MAX(chem(i,k,j,P_rtn28o2),0.), KIND=dp)  
    var(ind_RTN28NO3) = conv  * REAL( MAX(chem(i,k,j,P_rtn28no3),0.), KIND=dp)  
    var(ind_RTN28OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtn28ooh),0.), KIND=dp)  
    var(ind_TNCARB26) = conv  * REAL( MAX(chem(i,k,j,P_tncarb26),0.), KIND=dp)  
    var(ind_RTN26O2) = conv  * REAL( MAX(chem(i,k,j,P_rtn26o2),0.), KIND=dp)  
    var(ind_RTN26OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtn26ooh),0.), KIND=dp)  
    var(ind_NRTN28O2) = conv  * REAL( MAX(chem(i,k,j,P_nrtn28o2),0.), KIND=dp)  
    var(ind_NRTN28OOH) = conv  * REAL( MAX(chem(i,k,j,P_nrtn28ooh),0.), KIND=dp)  
    var(ind_RTN26PAN) = conv  * REAL( MAX(chem(i,k,j,P_rtn26pan),0.), KIND=dp)  
    var(ind_RTN25O2) = conv  * REAL( MAX(chem(i,k,j,P_rtn25o2),0.), KIND=dp)  
    var(ind_RTN24O2) = conv  * REAL( MAX(chem(i,k,j,P_rtn24o2),0.), KIND=dp)  
    var(ind_RTN23O2) = conv  * REAL( MAX(chem(i,k,j,P_rtn23o2),0.), KIND=dp)  
    var(ind_RTN14O2) = conv  * REAL( MAX(chem(i,k,j,P_rtn14o2),0.), KIND=dp)  
    var(ind_RTN10O2) = conv  * REAL( MAX(chem(i,k,j,P_rtn10o2),0.), KIND=dp)  
    var(ind_RTN25OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtn25ooh),0.), KIND=dp)  
    var(ind_RTN24OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtn24ooh),0.), KIND=dp)  
    var(ind_RTN23OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtn23ooh),0.), KIND=dp)  
    var(ind_RTN14OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtn14ooh),0.), KIND=dp)  
    var(ind_RTN10OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtn10ooh),0.), KIND=dp)  
    var(ind_TNCARB10) = conv  * REAL( MAX(chem(i,k,j,P_tncarb10),0.), KIND=dp)  
    var(ind_RTN25NO3) = conv  * REAL( MAX(chem(i,k,j,P_rtn25no3),0.), KIND=dp)  
    var(ind_TNCARB15) = conv  * REAL( MAX(chem(i,k,j,P_tncarb15),0.), KIND=dp)  
    var(ind_RCOOH25) = conv  * REAL( MAX(chem(i,k,j,P_rcooh25),0.), KIND=dp)  
    var(ind_RTX28O2) = conv  * REAL( MAX(chem(i,k,j,P_rtx28o2),0.), KIND=dp)  
    var(ind_RTX28NO3) = conv  * REAL( MAX(chem(i,k,j,P_rtx28no3),0.), KIND=dp)  
    var(ind_RTX28OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtx28ooh),0.), KIND=dp)  
    var(ind_TXCARB24) = conv  * REAL( MAX(chem(i,k,j,P_txcarb24),0.), KIND=dp)  
    var(ind_RTX24O2) = conv  * REAL( MAX(chem(i,k,j,P_rtx24o2),0.), KIND=dp)  
    var(ind_RTX24NO3) = conv  * REAL( MAX(chem(i,k,j,P_rtx24no3),0.), KIND=dp)  
    var(ind_RTX24OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtx24ooh),0.), KIND=dp)  
    var(ind_TXCARB22) = conv  * REAL( MAX(chem(i,k,j,P_txcarb22),0.), KIND=dp)  
    var(ind_RTX22O2) = conv  * REAL( MAX(chem(i,k,j,P_rtx22o2),0.), KIND=dp)  
    var(ind_RTX22NO3) = conv  * REAL( MAX(chem(i,k,j,P_rtx22no3),0.), KIND=dp)  
    var(ind_RTX22OOH) = conv  * REAL( MAX(chem(i,k,j,P_rtx22ooh),0.), KIND=dp)  
    var(ind_NRTX28O2) = conv  * REAL( MAX(chem(i,k,j,P_nrtx28o2),0.), KIND=dp)  
    var(ind_NRTX28OOH) = conv  * REAL( MAX(chem(i,k,j,P_nrtx28ooh),0.), KIND=dp)  
    var(ind_CARB11A) = conv  * REAL( MAX(chem(i,k,j,P_carb11a),0.), KIND=dp)  
    var(ind_ANHY) = conv  * REAL( MAX(chem(i,k,j,P_anhy),0.), KIND=dp)  
    var(ind_CH3O2NO2) = conv  * REAL( MAX(chem(i,k,j,P_ch3o2no2),0.), KIND=dp)  
    var(ind_CH4) = conv  * REAL( MAX(chem(i,k,j,P_ch4),0.), KIND=dp)  
    var(ind_H2SO4) = conv  * REAL( MAX(chem(i,k,j,P_SULF),0.), KIND=dp)  
    var(ind_HCl) = conv  * REAL( MAX(chem(i,k,j,P_hcl),0.), KIND=dp)  
    var(ind_NH3) = conv  * REAL( MAX(chem(i,k,j,P_nh3),0.), KIND=dp)  
    var(ind_RTN23NO3) = conv  * REAL( MAX(chem(i,k,j,P_rtn23no3),0.), KIND=dp)  
    var(ind_TNCARB12) = conv  * REAL( MAX(chem(i,k,j,P_tncarb12),0.), KIND=dp)  
    var(ind_TNCARB11) = conv  * REAL( MAX(chem(i,k,j,P_tncarb11),0.), KIND=dp)  
    var(ind_TM123B) = conv  * REAL( MAX(chem(i,k,j,P_tm123b),0.), KIND=dp)  
    var(ind_TM124B) = conv  * REAL( MAX(chem(i,k,j,P_tm124b),0.), KIND=dp)  
    var(ind_TM135B) = conv  * REAL( MAX(chem(i,k,j,P_tm135b),0.), KIND=dp)  
    var(ind_OETHTOL) = conv  * REAL( MAX(chem(i,k,j,P_oethtol),0.), KIND=dp)  
    var(ind_METHTOL) = conv  * REAL( MAX(chem(i,k,j,P_methtol),0.), KIND=dp)  
    var(ind_PETHTOL) = conv  * REAL( MAX(chem(i,k,j,P_pethtol),0.), KIND=dp)  
    var(ind_RA22AO2) = conv  * REAL( MAX(chem(i,k,j,P_ra22ao2),0.), KIND=dp)  
    var(ind_RA22BO2) = conv  * REAL( MAX(chem(i,k,j,P_ra22bo2),0.), KIND=dp)  
    var(ind_RA22NO3) = conv  * REAL( MAX(chem(i,k,j,P_ra22no3),0.), KIND=dp)  
    var(ind_RA22OOH) = conv  * REAL( MAX(chem(i,k,j,P_ra22ooh),0.), KIND=dp)  
    var(ind_DIME35EB) = conv  * REAL( MAX(chem(i,k,j,P_dime35eb),0.), KIND=dp)  
    var(ind_RA25O2) = conv  * REAL( MAX(chem(i,k,j,P_ra25o2),0.), KIND=dp)  
    var(ind_RA25NO3) = conv  * REAL( MAX(chem(i,k,j,P_ra25no3),0.), KIND=dp)  
    var(ind_UDCARB17) = conv  * REAL( MAX(chem(i,k,j,P_udcarb17),0.), KIND=dp)  
    var(ind_RA25OOH) = conv  * REAL( MAX(chem(i,k,j,P_ra25ooh),0.), KIND=dp)  
    var(ind_DMS) = conv  * REAL( MAX(chem(i,k,j,P_dms),0.), KIND=dp)  
    var(ind_CH3SCH2OO) = conv  * REAL( MAX(chem(i,k,j,P_ch3sch2oo),0.), KIND=dp)  
    var(ind_DMSO) = conv  * REAL( MAX(chem(i,k,j,P_dmso),0.), KIND=dp)  
    var(ind_CH3S) = conv  * REAL( MAX(chem(i,k,j,P_ch3s),0.), KIND=dp)  
    var(ind_CH3SO) = conv  * REAL( MAX(chem(i,k,j,P_ch3so),0.), KIND=dp)  
    var(ind_CH3SO2) = conv  * REAL( MAX(chem(i,k,j,P_ch3so2),0.), KIND=dp)  
    var(ind_CH3SO3) = conv  * REAL( MAX(chem(i,k,j,P_ch3so3),0.), KIND=dp)  
    var(ind_MSA) = conv  * REAL( MAX(chem(i,k,j,P_msa),0.), KIND=dp)  
    var(ind_MSIA) = conv  * REAL( MAX(chem(i,k,j,P_msia),0.), KIND=dp)  
    var(ind_DMSO2) = conv  * REAL( MAX(chem(i,k,j,P_dmso2),0.), KIND=dp)  
    var(ind_ClNO2) = conv  * REAL( MAX(chem(i,k,j,P_clno2),0.), KIND=dp)  




    RO2 = & 
        REAL(var(ind_CH3OO) + var(ind_C2H5O2) + var(ind_RN10O2) + var(ind_IC3H7O2) &
        + var(ind_RN13O2) + var(ind_RN13AO2) + var(ind_RN16AO2)  + var(ind_RN16O2) &
        + var(ind_RN19O2) + var(ind_HOCH2CH2O2) + var(ind_RN9O2) + var(ind_RN12O2) &
        + var(ind_RN15O2)+ var(ind_RN18O2)+ var(ind_RN15AO2) + var(ind_RN18AO2) &
        + var(ind_CH3CO3) + var(ind_C2H5CO3)+ var(ind_RN11O2)+ var(ind_RN14O2) &
        + var(ind_RN17O2) + var(ind_HOCH2CO3)+ var(ind_RU14O2)+ var(ind_RU12O2) &
        + var(ind_RU10O2)+ var(ind_NRN6O2)+ var(ind_NRN9O2)+ var(ind_NRN12O2) &
        + var(ind_RTN28O2) + var(ind_NRU14O2)+ var(ind_NRU12O2) + var(ind_RA13O2) &
        + var(ind_RA16O2) + var(ind_RA19AO2) + var(ind_RA19CO2)+ var(ind_RN8O2) &
        + var(ind_RTN26O2)  + var(ind_NRTN28O2) + var(ind_RTN25O2) + var(ind_RTN24O2) &
        + var(ind_RTN23O2)+ var(ind_RTN14O2)+ var(ind_RTN10O2)+ var(ind_RTX28O2) &
        + var(ind_RTX24O2)+ var(ind_RTX22O2)+ var(ind_NRTX28O2), KIND=dp)






   CALL cri_mosaic_8bin_aq_Update_Rconst(  &


RO2, &

var,nvar,ind_ru14o2,ind_no,ind_oh,ind_rtn23o2,ind_rtn26o2,ind_ho2, &
ind_rtx28o2,ind_rn19o2,ind_nrn12ooh,ind_hoch2co3,ind_no2, &



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








  CALL cri_mosaic_8bin_aq_INTEGRATE(TIME_START, TIME_END, &  
          FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK, & 
          ICNTRL_U=icntrl, RCNTRL_U=rcntrl  )







    chem(i,k,j,P_hono) = MAX ( REAL (oconv * var(ind_HONO), KIND=sp), 0.)  
    chem(i,k,j,P_o3) = MAX ( REAL (oconv * var(ind_O3), KIND=sp), 0.)  
    chem(i,k,j,P_hcho) = MAX ( REAL (oconv * var(ind_HCHO), KIND=sp), 0.)  
    chem(i,k,j,P_pan) = MAX ( REAL (oconv * var(ind_PAN), KIND=sp), 0.)  
    chem(i,k,j,P_c2h4) = MAX ( REAL (oconv * var(ind_C2H4), KIND=sp), 0.)  
    chem(i,k,j,P_co) = MAX ( REAL (oconv * var(ind_CO), KIND=sp), 0.)  
    chem(i,k,j,P_hno3) = MAX ( REAL (oconv * var(ind_HNO3), KIND=sp), 0.)  
    chem(i,k,j,P_n2o5) = MAX ( REAL (oconv * var(ind_N2O5), KIND=sp), 0.)  
    chem(i,k,j,P_hno4) = MAX ( REAL (oconv * var(ind_HNO4), KIND=sp), 0.)  
    chem(i,k,j,P_no3) = MAX ( REAL (oconv * var(ind_NO3), KIND=sp), 0.)  
   o1d(i,k,j) = MAX (REAL (oconv * var(ind_O1D) , KIND=sp),0.) 
   o3p(i,k,j) = MAX (REAL (oconv * var(ind_O3P) , KIND=sp),0.) 
    chem(i,k,j,P_HO) = MAX ( REAL (oconv * var(ind_OH), KIND=sp), 0.)  
    chem(i,k,j,P_ho2) = MAX ( REAL (oconv * var(ind_HO2), KIND=sp), 0.)  
    chem(i,k,j,P_h2o2) = MAX ( REAL (oconv * var(ind_H2O2), KIND=sp), 0.)  
    chem(i,k,j,P_c2h6) = MAX ( REAL (oconv * var(ind_C2H6), KIND=sp), 0.)  
    chem(i,k,j,P_hcooh) = MAX ( REAL (oconv * var(ind_HCOOH), KIND=sp), 0.)  
    chem(i,k,j,P_ACO3) = MAX ( REAL (oconv * var(ind_CH3CO3), KIND=sp), 0.)  
    chem(i,k,j,P_ch3oo) = MAX ( REAL (oconv * var(ind_CH3OO), KIND=sp), 0.)  
    chem(i,k,j,P_c2h5o2) = MAX ( REAL (oconv * var(ind_C2H5O2), KIND=sp), 0.)  
    chem(i,k,j,P_hso3) = MAX ( REAL (oconv * var(ind_HSO3), KIND=sp), 0.)  
    chem(i,k,j,P_so3) = MAX ( REAL (oconv * var(ind_SO3), KIND=sp), 0.)  
    chem(i,k,j,P_so2) = MAX ( REAL (oconv * var(ind_SO2), KIND=sp), 0.)  
    chem(i,k,j,P_no2) = MAX ( REAL (oconv * var(ind_NO2), KIND=sp), 0.)  
    chem(i,k,j,P_no) = MAX ( REAL (oconv * var(ind_NO), KIND=sp), 0.)  
    chem(i,k,j,P_c3h8) = MAX ( REAL (oconv * var(ind_C3H8), KIND=sp), 0.)  
    chem(i,k,j,P_nc4h10) = MAX ( REAL (oconv * var(ind_NC4H10), KIND=sp), 0.)  
    chem(i,k,j,P_hoch2ch2o2) = MAX ( REAL (oconv * var(ind_HOCH2CH2O2), KIND=sp), 0.)  
    chem(i,k,j,P_ic3h7o2) = MAX ( REAL (oconv * var(ind_IC3H7O2), KIND=sp), 0.)  
    chem(i,k,j,P_c5h8) = MAX ( REAL (oconv * var(ind_C5H8), KIND=sp), 0.)  
    chem(i,k,j,P_benzene) = MAX ( REAL (oconv * var(ind_BENZENE), KIND=sp), 0.)  
    chem(i,k,j,P_toluene) = MAX ( REAL (oconv * var(ind_TOLUENE), KIND=sp), 0.)  
    chem(i,k,j,P_oxyl) = MAX ( REAL (oconv * var(ind_OXYL), KIND=sp), 0.)  
    chem(i,k,j,P_npropol) = MAX ( REAL (oconv * var(ind_NPROPOL), KIND=sp), 0.)  
    chem(i,k,j,P_c2h2) = MAX ( REAL (oconv * var(ind_C2H2), KIND=sp), 0.)  
    chem(i,k,j,P_c3h6) = MAX ( REAL (oconv * var(ind_C3H6), KIND=sp), 0.)  
    chem(i,k,j,P_tbut2ene) = MAX ( REAL (oconv * var(ind_TBUT2ENE), KIND=sp), 0.)  
    chem(i,k,j,P_ch3cho) = MAX ( REAL (oconv * var(ind_CH3CHO), KIND=sp), 0.)  
    chem(i,k,j,P_c2h5cho) = MAX ( REAL (oconv * var(ind_C2H5CHO), KIND=sp), 0.)  
    chem(i,k,j,P_ch3co2h) = MAX ( REAL (oconv * var(ind_CH3CO2H), KIND=sp), 0.)  
    chem(i,k,j,P_KET) = MAX ( REAL (oconv * var(ind_CH3COCH3), KIND=sp), 0.)  
    chem(i,k,j,P_mek) = MAX ( REAL (oconv * var(ind_MEK), KIND=sp), 0.)  
    chem(i,k,j,P_ch3oh) = MAX ( REAL (oconv * var(ind_CH3OH), KIND=sp), 0.)  
    chem(i,k,j,P_c2h5oh) = MAX ( REAL (oconv * var(ind_C2H5OH), KIND=sp), 0.)  
    chem(i,k,j,P_ic3h7no3) = MAX ( REAL (oconv * var(ind_IC3H7NO3), KIND=sp), 0.)  
    chem(i,k,j,P_ipropol) = MAX ( REAL (oconv * var(ind_IPROPOL), KIND=sp), 0.)  
    chem(i,k,j,P_ch3no3) = MAX ( REAL (oconv * var(ind_CH3NO3), KIND=sp), 0.)  
    chem(i,k,j,P_c2h5no3) = MAX ( REAL (oconv * var(ind_C2H5NO3), KIND=sp), 0.)  
    chem(i,k,j,P_hoc2h4no3) = MAX ( REAL (oconv * var(ind_HOC2H4NO3), KIND=sp), 0.)  
    chem(i,k,j,P_ch3ooh) = MAX ( REAL (oconv * var(ind_CH3OOH), KIND=sp), 0.)  
    chem(i,k,j,P_c2h5ooh) = MAX ( REAL (oconv * var(ind_C2H5OOH), KIND=sp), 0.)  
    chem(i,k,j,P_PROOH) = MAX ( REAL (oconv * var(ind_IC3H7OOH), KIND=sp), 0.)  
    chem(i,k,j,P_PAA) = MAX ( REAL (oconv * var(ind_CH3CO3H), KIND=sp), 0.)  
    chem(i,k,j,P_hoc2h4ooh) = MAX ( REAL (oconv * var(ind_HOC2H4OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn10o2) = MAX ( REAL (oconv * var(ind_RN10O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn13o2) = MAX ( REAL (oconv * var(ind_RN13O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn16o2) = MAX ( REAL (oconv * var(ind_RN16O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn19o2) = MAX ( REAL (oconv * var(ind_RN19O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn9o2) = MAX ( REAL (oconv * var(ind_RN9O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn12o2) = MAX ( REAL (oconv * var(ind_RN12O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn15o2) = MAX ( REAL (oconv * var(ind_RN15O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn18o2) = MAX ( REAL (oconv * var(ind_RN18O2), KIND=sp), 0.)  
    chem(i,k,j,P_nrn6o2) = MAX ( REAL (oconv * var(ind_NRN6O2), KIND=sp), 0.)  
    chem(i,k,j,P_nrn9o2) = MAX ( REAL (oconv * var(ind_NRN9O2), KIND=sp), 0.)  
    chem(i,k,j,P_nrn12o2) = MAX ( REAL (oconv * var(ind_NRN12O2), KIND=sp), 0.)  
    chem(i,k,j,P_carb14) = MAX ( REAL (oconv * var(ind_CARB14), KIND=sp), 0.)  
    chem(i,k,j,P_rn11o2) = MAX ( REAL (oconv * var(ind_RN11O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn14o2) = MAX ( REAL (oconv * var(ind_RN14O2), KIND=sp), 0.)  
    chem(i,k,j,P_carb17) = MAX ( REAL (oconv * var(ind_CARB17), KIND=sp), 0.)  
    chem(i,k,j,P_rn8o2) = MAX ( REAL (oconv * var(ind_RN8O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn17o2) = MAX ( REAL (oconv * var(ind_RN17O2), KIND=sp), 0.)  
    chem(i,k,j,P_rn10no3) = MAX ( REAL (oconv * var(ind_RN10NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rn13no3) = MAX ( REAL (oconv * var(ind_RN13NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rn19no3) = MAX ( REAL (oconv * var(ind_RN19NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rn9no3) = MAX ( REAL (oconv * var(ind_RN9NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rn12no3) = MAX ( REAL (oconv * var(ind_RN12NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rn15no3) = MAX ( REAL (oconv * var(ind_RN15NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rn18no3) = MAX ( REAL (oconv * var(ind_RN18NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rn16no3) = MAX ( REAL (oconv * var(ind_RN16NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rn10ooh) = MAX ( REAL (oconv * var(ind_RN10OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn13ooh) = MAX ( REAL (oconv * var(ind_RN13OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn16ooh) = MAX ( REAL (oconv * var(ind_RN16OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn19ooh) = MAX ( REAL (oconv * var(ind_RN19OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn8ooh) = MAX ( REAL (oconv * var(ind_RN8OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn11ooh) = MAX ( REAL (oconv * var(ind_RN11OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn14ooh) = MAX ( REAL (oconv * var(ind_RN14OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn17ooh) = MAX ( REAL (oconv * var(ind_RN17OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn9ooh) = MAX ( REAL (oconv * var(ind_RN9OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn12ooh) = MAX ( REAL (oconv * var(ind_RN12OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn15ooh) = MAX ( REAL (oconv * var(ind_RN15OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rn18ooh) = MAX ( REAL (oconv * var(ind_RN18OOH), KIND=sp), 0.)  
    chem(i,k,j,P_nrn6ooh) = MAX ( REAL (oconv * var(ind_NRN6OOH), KIND=sp), 0.)  
    chem(i,k,j,P_nrn9ooh) = MAX ( REAL (oconv * var(ind_NRN9OOH), KIND=sp), 0.)  
    chem(i,k,j,P_nrn12ooh) = MAX ( REAL (oconv * var(ind_NRN12OOH), KIND=sp), 0.)  
    chem(i,k,j,P_apinene) = MAX ( REAL (oconv * var(ind_APINENE), KIND=sp), 0.)  
    chem(i,k,j,P_bpinene) = MAX ( REAL (oconv * var(ind_BPINENE), KIND=sp), 0.)  
    chem(i,k,j,P_rn13ao2) = MAX ( REAL (oconv * var(ind_RN13AO2), KIND=sp), 0.)  
    chem(i,k,j,P_rn16ao2) = MAX ( REAL (oconv * var(ind_RN16AO2), KIND=sp), 0.)  
    chem(i,k,j,P_rn15ao2) = MAX ( REAL (oconv * var(ind_RN15AO2), KIND=sp), 0.)  
    chem(i,k,j,P_rn18ao2) = MAX ( REAL (oconv * var(ind_RN18AO2), KIND=sp), 0.)  
    chem(i,k,j,P_carb7) = MAX ( REAL (oconv * var(ind_CARB7), KIND=sp), 0.)  
    chem(i,k,j,P_carb10) = MAX ( REAL (oconv * var(ind_CARB10), KIND=sp), 0.)  
    chem(i,k,j,P_carb13) = MAX ( REAL (oconv * var(ind_CARB13), KIND=sp), 0.)  
    chem(i,k,j,P_carb16) = MAX ( REAL (oconv * var(ind_CARB16), KIND=sp), 0.)  
    chem(i,k,j,P_carb3) = MAX ( REAL (oconv * var(ind_CARB3), KIND=sp), 0.)  
    chem(i,k,j,P_carb6) = MAX ( REAL (oconv * var(ind_CARB6), KIND=sp), 0.)  
    chem(i,k,j,P_carb9) = MAX ( REAL (oconv * var(ind_CARB9), KIND=sp), 0.)  
    chem(i,k,j,P_carb12) = MAX ( REAL (oconv * var(ind_CARB12), KIND=sp), 0.)  
    chem(i,k,j,P_carb15) = MAX ( REAL (oconv * var(ind_CARB15), KIND=sp), 0.)  
    chem(i,k,j,P_c2h5co3h) = MAX ( REAL (oconv * var(ind_C2H5CO3H), KIND=sp), 0.)  
    chem(i,k,j,P_c2h5co3) = MAX ( REAL (oconv * var(ind_C2H5CO3), KIND=sp), 0.)  
    chem(i,k,j,P_ppn) = MAX ( REAL (oconv * var(ind_PPN), KIND=sp), 0.)  
    chem(i,k,j,P_hoch2cho) = MAX ( REAL (oconv * var(ind_HOCH2CHO), KIND=sp), 0.)  
    chem(i,k,j,P_hoch2co3) = MAX ( REAL (oconv * var(ind_HOCH2CO3), KIND=sp), 0.)  
    chem(i,k,j,P_hoch2co3h) = MAX ( REAL (oconv * var(ind_HOCH2CO3H), KIND=sp), 0.)  
    chem(i,k,j,P_phan) = MAX ( REAL (oconv * var(ind_PHAN), KIND=sp), 0.)  
    chem(i,k,j,P_ccarb12) = MAX ( REAL (oconv * var(ind_CCARB12), KIND=sp), 0.)  
    chem(i,k,j,P_ru14o2) = MAX ( REAL (oconv * var(ind_RU14O2), KIND=sp), 0.)  
    chem(i,k,j,P_ru12o2) = MAX ( REAL (oconv * var(ind_RU12O2), KIND=sp), 0.)  
    chem(i,k,j,P_ch3cl) = MAX ( REAL (oconv * var(ind_CH3CL), KIND=sp), 0.)  
    chem(i,k,j,P_ch2cl2) = MAX ( REAL (oconv * var(ind_CH2CL2), KIND=sp), 0.)  
    chem(i,k,j,P_chcl3) = MAX ( REAL (oconv * var(ind_CHCL3), KIND=sp), 0.)  
    chem(i,k,j,P_ch3ccl3) = MAX ( REAL (oconv * var(ind_CH3CCL3), KIND=sp), 0.)  
    chem(i,k,j,P_cdicleth) = MAX ( REAL (oconv * var(ind_CDICLETH), KIND=sp), 0.)  
    chem(i,k,j,P_tdicleth) = MAX ( REAL (oconv * var(ind_TDICLETH), KIND=sp), 0.)  
    chem(i,k,j,P_tricleth) = MAX ( REAL (oconv * var(ind_TRICLETH), KIND=sp), 0.)  
    chem(i,k,j,P_tce) = MAX ( REAL (oconv * var(ind_TCE), KIND=sp), 0.)  
    chem(i,k,j,P_ru10o2) = MAX ( REAL (oconv * var(ind_RU10O2), KIND=sp), 0.)  
    chem(i,k,j,P_ucarb12) = MAX ( REAL (oconv * var(ind_UCARB12), KIND=sp), 0.)  
    chem(i,k,j,P_ucarb10) = MAX ( REAL (oconv * var(ind_UCARB10), KIND=sp), 0.)  
    chem(i,k,j,P_ru14no3) = MAX ( REAL (oconv * var(ind_RU14NO3), KIND=sp), 0.)  
    chem(i,k,j,P_ru14ooh) = MAX ( REAL (oconv * var(ind_RU14OOH), KIND=sp), 0.)  
    chem(i,k,j,P_ru12ooh) = MAX ( REAL (oconv * var(ind_RU12OOH), KIND=sp), 0.)  
    chem(i,k,j,P_ru10ooh) = MAX ( REAL (oconv * var(ind_RU10OOH), KIND=sp), 0.)  
    chem(i,k,j,P_mpan) = MAX ( REAL (oconv * var(ind_MPAN), KIND=sp), 0.)  
    chem(i,k,j,P_ru12pan) = MAX ( REAL (oconv * var(ind_RU12PAN), KIND=sp), 0.)  
    chem(i,k,j,P_nru14o2) = MAX ( REAL (oconv * var(ind_NRU14O2), KIND=sp), 0.)  
    chem(i,k,j,P_nucarb12) = MAX ( REAL (oconv * var(ind_NUCARB12), KIND=sp), 0.)  
    chem(i,k,j,P_nru14ooh) = MAX ( REAL (oconv * var(ind_NRU14OOH), KIND=sp), 0.)  
    chem(i,k,j,P_nru12o2) = MAX ( REAL (oconv * var(ind_NRU12O2), KIND=sp), 0.)  
    chem(i,k,j,P_nru12ooh) = MAX ( REAL (oconv * var(ind_NRU12OOH), KIND=sp), 0.)  
    chem(i,k,j,P_noa) = MAX ( REAL (oconv * var(ind_NOA), KIND=sp), 0.)  
    chem(i,k,j,P_ra13o2) = MAX ( REAL (oconv * var(ind_RA13O2), KIND=sp), 0.)  
    chem(i,k,j,P_ra13no3) = MAX ( REAL (oconv * var(ind_RA13NO3), KIND=sp), 0.)  
    chem(i,k,j,P_ra13ooh) = MAX ( REAL (oconv * var(ind_RA13OOH), KIND=sp), 0.)  
    chem(i,k,j,P_udcarb8) = MAX ( REAL (oconv * var(ind_UDCARB8), KIND=sp), 0.)  
    chem(i,k,j,P_aroh14) = MAX ( REAL (oconv * var(ind_AROH14), KIND=sp), 0.)  
    chem(i,k,j,P_raroh14) = MAX ( REAL (oconv * var(ind_RAROH14), KIND=sp), 0.)  
    chem(i,k,j,P_arnoh14) = MAX ( REAL (oconv * var(ind_ARNOH14), KIND=sp), 0.)  
    chem(i,k,j,P_ra16o2) = MAX ( REAL (oconv * var(ind_RA16O2), KIND=sp), 0.)  
    chem(i,k,j,P_ra16no3) = MAX ( REAL (oconv * var(ind_RA16NO3), KIND=sp), 0.)  
    chem(i,k,j,P_ra16ooh) = MAX ( REAL (oconv * var(ind_RA16OOH), KIND=sp), 0.)  
    chem(i,k,j,P_udcarb11) = MAX ( REAL (oconv * var(ind_UDCARB11), KIND=sp), 0.)  
    chem(i,k,j,P_aroh17) = MAX ( REAL (oconv * var(ind_AROH17), KIND=sp), 0.)  
    chem(i,k,j,P_raroh17) = MAX ( REAL (oconv * var(ind_RAROH17), KIND=sp), 0.)  
    chem(i,k,j,P_arnoh17) = MAX ( REAL (oconv * var(ind_ARNOH17), KIND=sp), 0.)  
    chem(i,k,j,P_udcarb14) = MAX ( REAL (oconv * var(ind_UDCARB14), KIND=sp), 0.)  
    chem(i,k,j,P_ra19ao2) = MAX ( REAL (oconv * var(ind_RA19AO2), KIND=sp), 0.)  
    chem(i,k,j,P_ra19co2) = MAX ( REAL (oconv * var(ind_RA19CO2), KIND=sp), 0.)  
    chem(i,k,j,P_ra19no3) = MAX ( REAL (oconv * var(ind_RA19NO3), KIND=sp), 0.)  
    chem(i,k,j,P_ra19ooh) = MAX ( REAL (oconv * var(ind_RA19OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rtn28o2) = MAX ( REAL (oconv * var(ind_RTN28O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtn28no3) = MAX ( REAL (oconv * var(ind_RTN28NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rtn28ooh) = MAX ( REAL (oconv * var(ind_RTN28OOH), KIND=sp), 0.)  
    chem(i,k,j,P_tncarb26) = MAX ( REAL (oconv * var(ind_TNCARB26), KIND=sp), 0.)  
    chem(i,k,j,P_rtn26o2) = MAX ( REAL (oconv * var(ind_RTN26O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtn26ooh) = MAX ( REAL (oconv * var(ind_RTN26OOH), KIND=sp), 0.)  
    chem(i,k,j,P_nrtn28o2) = MAX ( REAL (oconv * var(ind_NRTN28O2), KIND=sp), 0.)  
    chem(i,k,j,P_nrtn28ooh) = MAX ( REAL (oconv * var(ind_NRTN28OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rtn26pan) = MAX ( REAL (oconv * var(ind_RTN26PAN), KIND=sp), 0.)  
    chem(i,k,j,P_rtn25o2) = MAX ( REAL (oconv * var(ind_RTN25O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtn24o2) = MAX ( REAL (oconv * var(ind_RTN24O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtn23o2) = MAX ( REAL (oconv * var(ind_RTN23O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtn14o2) = MAX ( REAL (oconv * var(ind_RTN14O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtn10o2) = MAX ( REAL (oconv * var(ind_RTN10O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtn25ooh) = MAX ( REAL (oconv * var(ind_RTN25OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rtn24ooh) = MAX ( REAL (oconv * var(ind_RTN24OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rtn23ooh) = MAX ( REAL (oconv * var(ind_RTN23OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rtn14ooh) = MAX ( REAL (oconv * var(ind_RTN14OOH), KIND=sp), 0.)  
    chem(i,k,j,P_rtn10ooh) = MAX ( REAL (oconv * var(ind_RTN10OOH), KIND=sp), 0.)  
    chem(i,k,j,P_tncarb10) = MAX ( REAL (oconv * var(ind_TNCARB10), KIND=sp), 0.)  
    chem(i,k,j,P_rtn25no3) = MAX ( REAL (oconv * var(ind_RTN25NO3), KIND=sp), 0.)  
    chem(i,k,j,P_tncarb15) = MAX ( REAL (oconv * var(ind_TNCARB15), KIND=sp), 0.)  
    chem(i,k,j,P_rcooh25) = MAX ( REAL (oconv * var(ind_RCOOH25), KIND=sp), 0.)  
    chem(i,k,j,P_rtx28o2) = MAX ( REAL (oconv * var(ind_RTX28O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtx28no3) = MAX ( REAL (oconv * var(ind_RTX28NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rtx28ooh) = MAX ( REAL (oconv * var(ind_RTX28OOH), KIND=sp), 0.)  
    chem(i,k,j,P_txcarb24) = MAX ( REAL (oconv * var(ind_TXCARB24), KIND=sp), 0.)  
    chem(i,k,j,P_rtx24o2) = MAX ( REAL (oconv * var(ind_RTX24O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtx24no3) = MAX ( REAL (oconv * var(ind_RTX24NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rtx24ooh) = MAX ( REAL (oconv * var(ind_RTX24OOH), KIND=sp), 0.)  
    chem(i,k,j,P_txcarb22) = MAX ( REAL (oconv * var(ind_TXCARB22), KIND=sp), 0.)  
    chem(i,k,j,P_rtx22o2) = MAX ( REAL (oconv * var(ind_RTX22O2), KIND=sp), 0.)  
    chem(i,k,j,P_rtx22no3) = MAX ( REAL (oconv * var(ind_RTX22NO3), KIND=sp), 0.)  
    chem(i,k,j,P_rtx22ooh) = MAX ( REAL (oconv * var(ind_RTX22OOH), KIND=sp), 0.)  
    chem(i,k,j,P_nrtx28o2) = MAX ( REAL (oconv * var(ind_NRTX28O2), KIND=sp), 0.)  
    chem(i,k,j,P_nrtx28ooh) = MAX ( REAL (oconv * var(ind_NRTX28OOH), KIND=sp), 0.)  
    chem(i,k,j,P_carb11a) = MAX ( REAL (oconv * var(ind_CARB11A), KIND=sp), 0.)  
    chem(i,k,j,P_anhy) = MAX ( REAL (oconv * var(ind_ANHY), KIND=sp), 0.)  
    chem(i,k,j,P_ch3o2no2) = MAX ( REAL (oconv * var(ind_CH3O2NO2), KIND=sp), 0.)  
    chem(i,k,j,P_ch4) = MAX ( REAL (oconv * var(ind_CH4), KIND=sp), 0.)  
    chem(i,k,j,P_SULF) = MAX ( REAL (oconv * var(ind_H2SO4), KIND=sp), 0.)  
    chem(i,k,j,P_hcl) = MAX ( REAL (oconv * var(ind_HCl), KIND=sp), 0.)  
    chem(i,k,j,P_nh3) = MAX ( REAL (oconv * var(ind_NH3), KIND=sp), 0.)  
    chem(i,k,j,P_rtn23no3) = MAX ( REAL (oconv * var(ind_RTN23NO3), KIND=sp), 0.)  
    chem(i,k,j,P_tncarb12) = MAX ( REAL (oconv * var(ind_TNCARB12), KIND=sp), 0.)  
    chem(i,k,j,P_tncarb11) = MAX ( REAL (oconv * var(ind_TNCARB11), KIND=sp), 0.)  
    chem(i,k,j,P_tm123b) = MAX ( REAL (oconv * var(ind_TM123B), KIND=sp), 0.)  
    chem(i,k,j,P_tm124b) = MAX ( REAL (oconv * var(ind_TM124B), KIND=sp), 0.)  
    chem(i,k,j,P_tm135b) = MAX ( REAL (oconv * var(ind_TM135B), KIND=sp), 0.)  
    chem(i,k,j,P_oethtol) = MAX ( REAL (oconv * var(ind_OETHTOL), KIND=sp), 0.)  
    chem(i,k,j,P_methtol) = MAX ( REAL (oconv * var(ind_METHTOL), KIND=sp), 0.)  
    chem(i,k,j,P_pethtol) = MAX ( REAL (oconv * var(ind_PETHTOL), KIND=sp), 0.)  
    chem(i,k,j,P_ra22ao2) = MAX ( REAL (oconv * var(ind_RA22AO2), KIND=sp), 0.)  
    chem(i,k,j,P_ra22bo2) = MAX ( REAL (oconv * var(ind_RA22BO2), KIND=sp), 0.)  
    chem(i,k,j,P_ra22no3) = MAX ( REAL (oconv * var(ind_RA22NO3), KIND=sp), 0.)  
    chem(i,k,j,P_ra22ooh) = MAX ( REAL (oconv * var(ind_RA22OOH), KIND=sp), 0.)  
    chem(i,k,j,P_dime35eb) = MAX ( REAL (oconv * var(ind_DIME35EB), KIND=sp), 0.)  
    chem(i,k,j,P_ra25o2) = MAX ( REAL (oconv * var(ind_RA25O2), KIND=sp), 0.)  
    chem(i,k,j,P_ra25no3) = MAX ( REAL (oconv * var(ind_RA25NO3), KIND=sp), 0.)  
    chem(i,k,j,P_udcarb17) = MAX ( REAL (oconv * var(ind_UDCARB17), KIND=sp), 0.)  
    chem(i,k,j,P_ra25ooh) = MAX ( REAL (oconv * var(ind_RA25OOH), KIND=sp), 0.)  
    chem(i,k,j,P_dms) = MAX ( REAL (oconv * var(ind_DMS), KIND=sp), 0.)  
    chem(i,k,j,P_ch3sch2oo) = MAX ( REAL (oconv * var(ind_CH3SCH2OO), KIND=sp), 0.)  
    chem(i,k,j,P_dmso) = MAX ( REAL (oconv * var(ind_DMSO), KIND=sp), 0.)  
    chem(i,k,j,P_ch3s) = MAX ( REAL (oconv * var(ind_CH3S), KIND=sp), 0.)  
    chem(i,k,j,P_ch3so) = MAX ( REAL (oconv * var(ind_CH3SO), KIND=sp), 0.)  
    chem(i,k,j,P_ch3so2) = MAX ( REAL (oconv * var(ind_CH3SO2), KIND=sp), 0.)  
    chem(i,k,j,P_ch3so3) = MAX ( REAL (oconv * var(ind_CH3SO3), KIND=sp), 0.)  
    chem(i,k,j,P_msa) = MAX ( REAL (oconv * var(ind_MSA), KIND=sp), 0.)  
    chem(i,k,j,P_msia) = MAX ( REAL (oconv * var(ind_MSIA), KIND=sp), 0.)  
    chem(i,k,j,P_dmso2) = MAX ( REAL (oconv * var(ind_DMSO2), KIND=sp), 0.)  
    chem(i,k,j,P_clno2) = MAX ( REAL (oconv * var(ind_ClNO2), KIND=sp), 0.)  



    END DO
    END DO
    END DO









END SUBROUTINE  cri_mosaic_8bin_aq_interface


END MODULE module_kpp_cri_mosaic_8bin_aq_interf 




