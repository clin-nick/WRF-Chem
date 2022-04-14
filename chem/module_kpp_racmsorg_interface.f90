







MODULE module_kpp_racmsorg_interf 


  USE module_state_description
  USE module_configure

  USE racmsorg_Parameters
  USE racmsorg_Precision
  USE racmsorg_UpdateRconstWRF
  USE racmsorg_Integrator

  USE module_wkppc_constants



  USE module_data_sorgam, ONLY : PXYL, PTOL, PCSL1, PCSL2, PHC8, POLI1, POLI2, POLI3, POLT1, POLT2, POLT3, PAPI1, PAPI2, PAPI3, PLIM1, PLIM2, PLIM3




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

SUBROUTINE  racmsorg_interface( &

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
 



 

      REAL(KIND=dp)  :: rxylho,rtolho,rcslho,rcslno3,rhc8ho,roliho,rolino3,      &
                        rolio3,roltho,roltno3,rolto3,rapiho,rapino3,rapio3,      &
                        rlimho,rlimno3,rlimo3

      REAL(KIND=dp) , DIMENSION(ldrog)  :: PRDROG 



 


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
    var(ind_SULF) = conv  * REAL( MAX(chem(i,k,j,P_sulf),0.), KIND=dp)  
    var(ind_CO) = conv  * REAL( MAX(chem(i,k,j,P_co),0.), KIND=dp)  
    var(ind_CH4) = conv  * REAL( MAX(chem(i,k,j,P_ch4),0.), KIND=dp)  
    var(ind_ETH) = conv  * REAL( MAX(chem(i,k,j,P_eth),0.), KIND=dp)  
    var(ind_HC3) = conv  * REAL( MAX(chem(i,k,j,P_hc3),0.), KIND=dp)  
    var(ind_HC5) = conv  * REAL( MAX(chem(i,k,j,P_hc5),0.), KIND=dp)  
    var(ind_HC8) = conv  * REAL( MAX(chem(i,k,j,P_hc8),0.), KIND=dp)  
    var(ind_ETE) = conv  * REAL( MAX(chem(i,k,j,P_ete),0.), KIND=dp)  
    var(ind_OLT) = conv  * REAL( MAX(chem(i,k,j,P_olt),0.), KIND=dp)  
    var(ind_OLI) = conv  * REAL( MAX(chem(i,k,j,P_oli),0.), KIND=dp)  
    var(ind_DIEN) = conv  * REAL( MAX(chem(i,k,j,P_dien),0.), KIND=dp)  
    var(ind_ISO) = conv  * REAL( MAX(chem(i,k,j,P_iso),0.), KIND=dp)  
    var(ind_API) = conv  * REAL( MAX(chem(i,k,j,P_api),0.), KIND=dp)  
    var(ind_LIM) = conv  * REAL( MAX(chem(i,k,j,P_lim),0.), KIND=dp)  
    var(ind_TOL) = conv  * REAL( MAX(chem(i,k,j,P_tol),0.), KIND=dp)  
    var(ind_XYL) = conv  * REAL( MAX(chem(i,k,j,P_xyl),0.), KIND=dp)  
    var(ind_CSL) = conv  * REAL( MAX(chem(i,k,j,P_csl),0.), KIND=dp)  
    var(ind_HCHO) = conv  * REAL( MAX(chem(i,k,j,P_hcho),0.), KIND=dp)  
    var(ind_ALD) = conv  * REAL( MAX(chem(i,k,j,P_ald),0.), KIND=dp)  
    var(ind_KET) = conv  * REAL( MAX(chem(i,k,j,P_ket),0.), KIND=dp)  
    var(ind_GLY) = conv  * REAL( MAX(chem(i,k,j,P_gly),0.), KIND=dp)  
    var(ind_MGLY) = conv  * REAL( MAX(chem(i,k,j,P_mgly),0.), KIND=dp)  
    var(ind_DCB) = conv  * REAL( MAX(chem(i,k,j,P_dcb),0.), KIND=dp)  
    var(ind_MACR) = conv  * REAL( MAX(chem(i,k,j,P_macr),0.), KIND=dp)  
    var(ind_UDD) = conv  * REAL( MAX(chem(i,k,j,P_udd),0.), KIND=dp)  
    var(ind_HKET) = conv  * REAL( MAX(chem(i,k,j,P_hket),0.), KIND=dp)  
    var(ind_ONIT) = conv  * REAL( MAX(chem(i,k,j,P_onit),0.), KIND=dp)  
    var(ind_PAN) = conv  * REAL( MAX(chem(i,k,j,P_pan),0.), KIND=dp)  
    var(ind_TPAN) = conv  * REAL( MAX(chem(i,k,j,P_tpan),0.), KIND=dp)  
    var(ind_OP1) = conv  * REAL( MAX(chem(i,k,j,P_op1),0.), KIND=dp)  
    var(ind_OP2) = conv  * REAL( MAX(chem(i,k,j,P_op2),0.), KIND=dp)  
    var(ind_PAA) = conv  * REAL( MAX(chem(i,k,j,P_paa),0.), KIND=dp)  
    var(ind_ORA1) = conv  * REAL( MAX(chem(i,k,j,P_ora1),0.), KIND=dp)  
    var(ind_ORA2) = conv  * REAL( MAX(chem(i,k,j,P_ora2),0.), KIND=dp)  
    var(ind_HO) = conv  * REAL( MAX(chem(i,k,j,P_ho),0.), KIND=dp)  
    var(ind_HO2) = conv  * REAL( MAX(chem(i,k,j,P_ho2),0.), KIND=dp)  
    var(ind_O3P) =  conv * REAL( MAX(o3p(i,k,j),0.), KIND=dp)  
    var(ind_O1D) =  conv * REAL( MAX(o1d(i,k,j),0.), KIND=dp)  
    var(ind_MO2) =  conv * REAL( MAX(mo2(i,k,j),0.), KIND=dp)  
    var(ind_ETHP) =  conv * REAL( MAX(ethp(i,k,j),0.), KIND=dp)  
    var(ind_HC3P) =  conv * REAL( MAX(hc3p(i,k,j),0.), KIND=dp)  
    var(ind_HC5P) =  conv * REAL( MAX(hc5p(i,k,j),0.), KIND=dp)  
    var(ind_HC8P) =  conv * REAL( MAX(hc8p(i,k,j),0.), KIND=dp)  
    var(ind_ETEP) =  conv * REAL( MAX(etep(i,k,j),0.), KIND=dp)  
    var(ind_OLTP) =  conv * REAL( MAX(oltp(i,k,j),0.), KIND=dp)  
    var(ind_OLIP) =  conv * REAL( MAX(olip(i,k,j),0.), KIND=dp)  
    var(ind_ISOP) =  conv * REAL( MAX(isop(i,k,j),0.), KIND=dp)  
    var(ind_APIP) =  conv * REAL( MAX(apip(i,k,j),0.), KIND=dp)  
    var(ind_LIMP) =  conv * REAL( MAX(limp(i,k,j),0.), KIND=dp)  
    var(ind_PHO) =  conv * REAL( MAX(RPHO(i,k,j),0.), KIND=dp)  
    var(ind_ADDT) =  conv * REAL( MAX(addt(i,k,j),0.), KIND=dp)  
    var(ind_ADDX) =  conv * REAL( MAX(addx(i,k,j),0.), KIND=dp)  
    var(ind_ADDC) =  conv * REAL( MAX(addc(i,k,j),0.), KIND=dp)  
    var(ind_TOLP) =  conv * REAL( MAX(tolp(i,k,j),0.), KIND=dp)  
    var(ind_XYLP) =  conv * REAL( MAX(xylp(i,k,j),0.), KIND=dp)  
    var(ind_CSLP) =  conv * REAL( MAX(cslp(i,k,j),0.), KIND=dp)  
    var(ind_ACO3) = conv  * REAL( MAX(chem(i,k,j,P_aco3),0.), KIND=dp)  
    var(ind_TCO3) =  conv * REAL( MAX(tco3(i,k,j),0.), KIND=dp)  
    var(ind_KETP) =  conv * REAL( MAX(ketp(i,k,j),0.), KIND=dp)  
    var(ind_XO2) =  conv * REAL( MAX(xo2(i,k,j),0.), KIND=dp)  
    var(ind_OLNN) =  conv * REAL( MAX(olnn(i,k,j),0.), KIND=dp)  
    var(ind_OLND) =  conv * REAL( MAX(olnd(i,k,j),0.), KIND=dp)  







   CALL racmsorg_Update_Rconst(  &



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








  CALL racmsorg_INTEGRATE(TIME_START, TIME_END, &  
          FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK, & 
          ICNTRL_U=icntrl, RCNTRL_U=rcntrl  )










        if(p_nu0.gt.1)then
            rxylho =   RCONST(74)     
            rtolho =   RCONST(73)     
            rcslho =   RCONST(75)     
            rcslno3 =  RCONST(97)     
            rhc8ho =   RCONST(65)     
            roliho =   RCONST(68)     
            rolino3 =  RCONST(100)    
            rolio3 =   RCONST(108)    
            roltho =   RCONST(67)     
            roltno3 =  RCONST(99)     
            rolto3 =   RCONST(107)    
            rapiho =   RCONST(71)     
            rapino3 =  RCONST(103)    
            rapio3 =   RCONST(111)    
            rlimho =   RCONST(72)     
            rlimno3 =  RCONST(104)    
            rlimo3 =   RCONST(112)    
            PRDROG(PXYL)  = rxylho * var(ind_xyl)*var(ind_ho)
            PRDROG(PTOL)  = rtolho * var(ind_tol)*var(ind_ho)
            PRDROG(PCSL1) = rcslho * var(ind_csl)*var(ind_ho)
            PRDROG(PCSL2) = 0.50_dp * rcslno3* var(ind_csl)*var(ind_no3)
            PRDROG(PHC8)  = rhc8ho * var(ind_hc8)*var(ind_ho)
            PRDROG(POLI1) = roliho * var(ind_oli)*var(ind_ho)
            PRDROG(POLI2) = rolino3 * var(ind_oli)*var(ind_no3)
            PRDROG(POLI3) = rolio3 * var(ind_oli)*var(ind_o3)
            PRDROG(POLT1) = roltho * var(ind_olt)*var(ind_ho)
            PRDROG(POLT2) = roltno3 * var(ind_olt)*var(ind_no3)
            PRDROG(POLT3) = rolto3 * var(ind_olt)*var(ind_o3)
            PRDROG(PAPI1) = rapiho * var(ind_api)*var(ind_ho)
            PRDROG(PAPI2) = rapino3 * var(ind_api)*var(ind_no3)
            PRDROG(PAPI3) = rapio3 * var(ind_api)*var(ind_o3)
            PRDROG(PLIM1) = rlimho * var(ind_lim)*var(ind_ho)
            PRDROG(PLIM2) = rlimno3 * var(ind_lim)*var(ind_no3)
            PRDROG(PLIM3) = rlimo3 * var(ind_lim)*var(ind_o3)
            DO n = 1, LDROG
               VDROG3( i,k,j, n ) =  oconv * PRDROG( n ) * DTSTEPC
               VDROG3( i,k,j,n  ) = MAX( 0., VDROG3( i,k,j, n ) )
            ENDDO
        endif


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
    chem(i,k,j,P_sulf) = MAX ( REAL (oconv * var(ind_SULF), KIND=sp), 0.)  
    chem(i,k,j,P_co) = MAX ( REAL (oconv * var(ind_CO), KIND=sp), 0.)  
    chem(i,k,j,P_ch4) = MAX ( REAL (oconv * var(ind_CH4), KIND=sp), 0.)  
    chem(i,k,j,P_eth) = MAX ( REAL (oconv * var(ind_ETH), KIND=sp), 0.)  
    chem(i,k,j,P_hc3) = MAX ( REAL (oconv * var(ind_HC3), KIND=sp), 0.)  
    chem(i,k,j,P_hc5) = MAX ( REAL (oconv * var(ind_HC5), KIND=sp), 0.)  
    chem(i,k,j,P_hc8) = MAX ( REAL (oconv * var(ind_HC8), KIND=sp), 0.)  
    chem(i,k,j,P_ete) = MAX ( REAL (oconv * var(ind_ETE), KIND=sp), 0.)  
    chem(i,k,j,P_olt) = MAX ( REAL (oconv * var(ind_OLT), KIND=sp), 0.)  
    chem(i,k,j,P_oli) = MAX ( REAL (oconv * var(ind_OLI), KIND=sp), 0.)  
    chem(i,k,j,P_dien) = MAX ( REAL (oconv * var(ind_DIEN), KIND=sp), 0.)  
    chem(i,k,j,P_iso) = MAX ( REAL (oconv * var(ind_ISO), KIND=sp), 0.)  
    chem(i,k,j,P_api) = MAX ( REAL (oconv * var(ind_API), KIND=sp), 0.)  
    chem(i,k,j,P_lim) = MAX ( REAL (oconv * var(ind_LIM), KIND=sp), 0.)  
    chem(i,k,j,P_tol) = MAX ( REAL (oconv * var(ind_TOL), KIND=sp), 0.)  
    chem(i,k,j,P_xyl) = MAX ( REAL (oconv * var(ind_XYL), KIND=sp), 0.)  
    chem(i,k,j,P_csl) = MAX ( REAL (oconv * var(ind_CSL), KIND=sp), 0.)  
    chem(i,k,j,P_hcho) = MAX ( REAL (oconv * var(ind_HCHO), KIND=sp), 0.)  
    chem(i,k,j,P_ald) = MAX ( REAL (oconv * var(ind_ALD), KIND=sp), 0.)  
    chem(i,k,j,P_ket) = MAX ( REAL (oconv * var(ind_KET), KIND=sp), 0.)  
    chem(i,k,j,P_gly) = MAX ( REAL (oconv * var(ind_GLY), KIND=sp), 0.)  
    chem(i,k,j,P_mgly) = MAX ( REAL (oconv * var(ind_MGLY), KIND=sp), 0.)  
    chem(i,k,j,P_dcb) = MAX ( REAL (oconv * var(ind_DCB), KIND=sp), 0.)  
    chem(i,k,j,P_macr) = MAX ( REAL (oconv * var(ind_MACR), KIND=sp), 0.)  
    chem(i,k,j,P_udd) = MAX ( REAL (oconv * var(ind_UDD), KIND=sp), 0.)  
    chem(i,k,j,P_hket) = MAX ( REAL (oconv * var(ind_HKET), KIND=sp), 0.)  
    chem(i,k,j,P_onit) = MAX ( REAL (oconv * var(ind_ONIT), KIND=sp), 0.)  
    chem(i,k,j,P_pan) = MAX ( REAL (oconv * var(ind_PAN), KIND=sp), 0.)  
    chem(i,k,j,P_tpan) = MAX ( REAL (oconv * var(ind_TPAN), KIND=sp), 0.)  
    chem(i,k,j,P_op1) = MAX ( REAL (oconv * var(ind_OP1), KIND=sp), 0.)  
    chem(i,k,j,P_op2) = MAX ( REAL (oconv * var(ind_OP2), KIND=sp), 0.)  
    chem(i,k,j,P_paa) = MAX ( REAL (oconv * var(ind_PAA), KIND=sp), 0.)  
    chem(i,k,j,P_ora1) = MAX ( REAL (oconv * var(ind_ORA1), KIND=sp), 0.)  
    chem(i,k,j,P_ora2) = MAX ( REAL (oconv * var(ind_ORA2), KIND=sp), 0.)  
    chem(i,k,j,P_ho) = MAX ( REAL (oconv * var(ind_HO), KIND=sp), 0.)  
    chem(i,k,j,P_ho2) = MAX ( REAL (oconv * var(ind_HO2), KIND=sp), 0.)  
   o3p(i,k,j) = MAX (REAL (oconv * var(ind_O3P) , KIND=sp),0.) 
   o1d(i,k,j) = MAX (REAL (oconv * var(ind_O1D) , KIND=sp),0.) 
   mo2(i,k,j) = MAX (REAL (oconv * var(ind_MO2) , KIND=sp),0.) 
   ethp(i,k,j) = MAX (REAL (oconv * var(ind_ETHP) , KIND=sp),0.) 
   hc3p(i,k,j) = MAX (REAL (oconv * var(ind_HC3P) , KIND=sp),0.) 
   hc5p(i,k,j) = MAX (REAL (oconv * var(ind_HC5P) , KIND=sp),0.) 
   hc8p(i,k,j) = MAX (REAL (oconv * var(ind_HC8P) , KIND=sp),0.) 
   etep(i,k,j) = MAX (REAL (oconv * var(ind_ETEP) , KIND=sp),0.) 
   oltp(i,k,j) = MAX (REAL (oconv * var(ind_OLTP) , KIND=sp),0.) 
   olip(i,k,j) = MAX (REAL (oconv * var(ind_OLIP) , KIND=sp),0.) 
   isop(i,k,j) = MAX (REAL (oconv * var(ind_ISOP) , KIND=sp),0.) 
   apip(i,k,j) = MAX (REAL (oconv * var(ind_APIP) , KIND=sp),0.) 
   limp(i,k,j) = MAX (REAL (oconv * var(ind_LIMP) , KIND=sp),0.) 
   RPHO(i,k,j) = MAX (REAL (oconv * var(ind_PHO) , KIND=sp),0.) 
   addt(i,k,j) = MAX (REAL (oconv * var(ind_ADDT) , KIND=sp),0.) 
   addx(i,k,j) = MAX (REAL (oconv * var(ind_ADDX) , KIND=sp),0.) 
   addc(i,k,j) = MAX (REAL (oconv * var(ind_ADDC) , KIND=sp),0.) 
   tolp(i,k,j) = MAX (REAL (oconv * var(ind_TOLP) , KIND=sp),0.) 
   xylp(i,k,j) = MAX (REAL (oconv * var(ind_XYLP) , KIND=sp),0.) 
   cslp(i,k,j) = MAX (REAL (oconv * var(ind_CSLP) , KIND=sp),0.) 
    chem(i,k,j,P_aco3) = MAX ( REAL (oconv * var(ind_ACO3), KIND=sp), 0.)  
   tco3(i,k,j) = MAX (REAL (oconv * var(ind_TCO3) , KIND=sp),0.) 
   ketp(i,k,j) = MAX (REAL (oconv * var(ind_KETP) , KIND=sp),0.) 
   xo2(i,k,j) = MAX (REAL (oconv * var(ind_XO2) , KIND=sp),0.) 
   olnn(i,k,j) = MAX (REAL (oconv * var(ind_OLNN) , KIND=sp),0.) 
   olnd(i,k,j) = MAX (REAL (oconv * var(ind_OLND) , KIND=sp),0.) 



    END DO
    END DO
    END DO









END SUBROUTINE  racmsorg_interface


END MODULE module_kpp_racmsorg_interf 




