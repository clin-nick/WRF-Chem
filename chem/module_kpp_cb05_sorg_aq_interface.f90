







MODULE module_kpp_cb05_sorg_aq_interf 


  USE module_state_description
  USE module_configure

  USE cb05_sorg_aq_Parameters
  USE cb05_sorg_aq_Precision
  USE cb05_sorg_aq_UpdateRconstWRF
  USE cb05_sorg_aq_Integrator

  USE module_wkppc_constants






  USE module_data_sorgam, ONLY : PXYL, PTOL, PCSL1, PCSL2, PHC8, &
        POLI1, POLI2, POLI3, POLT1, POLT2, POLT3, PAPI1, PAPI2,  &
        PAPI3, PLIM1, PLIM2, PLIM3





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

SUBROUTINE  cb05_sorg_aq_interface( &

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
 


    var(ind_no2) = conv  * REAL( MAX(chem(i,k,j,P_no2),0.), KIND=dp)  
    var(ind_no) = conv  * REAL( MAX(chem(i,k,j,P_no),0.), KIND=dp)  
    var(ind_o) = conv  * REAL( MAX(chem(i,k,j,P_o),0.), KIND=dp)  
    var(ind_o3) = conv  * REAL( MAX(chem(i,k,j,P_o3),0.), KIND=dp)  
    var(ind_no3) = conv  * REAL( MAX(chem(i,k,j,P_no3),0.), KIND=dp)  
    var(ind_o1d) = conv  * REAL( MAX(chem(i,k,j,P_o1d),0.), KIND=dp)  
    var(ind_oh) = conv  * REAL( MAX(chem(i,k,j,P_oh),0.), KIND=dp)  
    var(ind_ho2) = conv  * REAL( MAX(chem(i,k,j,P_ho2),0.), KIND=dp)  
    var(ind_n2o5) = conv  * REAL( MAX(chem(i,k,j,P_n2o5),0.), KIND=dp)  
    var(ind_hno3) = conv  * REAL( MAX(chem(i,k,j,P_hno3),0.), KIND=dp)  
    var(ind_hono) = conv  * REAL( MAX(chem(i,k,j,P_hono),0.), KIND=dp)  
    var(ind_pna) = conv  * REAL( MAX(chem(i,k,j,P_pna),0.), KIND=dp)  
    var(ind_h2o2) = conv  * REAL( MAX(chem(i,k,j,P_h2o2),0.), KIND=dp)  
    var(ind_xo2) = conv  * REAL( MAX(chem(i,k,j,P_xo2),0.), KIND=dp)  
    var(ind_xo2n) = conv  * REAL( MAX(chem(i,k,j,P_xo2n),0.), KIND=dp)  
    var(ind_ntr) = conv  * REAL( MAX(chem(i,k,j,P_ntr),0.), KIND=dp)  
    var(ind_rooh) = conv  * REAL( MAX(chem(i,k,j,P_rooh),0.), KIND=dp)  
    var(ind_form) = conv  * REAL( MAX(chem(i,k,j,P_form),0.), KIND=dp)  
    var(ind_ald2) = conv  * REAL( MAX(chem(i,k,j,P_ald2),0.), KIND=dp)  
    var(ind_aldx) = conv  * REAL( MAX(chem(i,k,j,P_aldx),0.), KIND=dp)  
    var(ind_par) = conv  * REAL( MAX(chem(i,k,j,P_par),0.), KIND=dp)  
    var(ind_co) = conv  * REAL( MAX(chem(i,k,j,P_co),0.), KIND=dp)  
    var(ind_meo2) = conv  * REAL( MAX(chem(i,k,j,P_meo2),0.), KIND=dp)  
    var(ind_mepx) = conv  * REAL( MAX(chem(i,k,j,P_mepx),0.), KIND=dp)  
    var(ind_meoh) = conv  * REAL( MAX(chem(i,k,j,P_meoh),0.), KIND=dp)  
    var(ind_hco3) = conv  * REAL( MAX(chem(i,k,j,P_hco3),0.), KIND=dp)  
    var(ind_facd) = conv  * REAL( MAX(chem(i,k,j,P_facd),0.), KIND=dp)  
    var(ind_c2o3) = conv  * REAL( MAX(chem(i,k,j,P_c2o3),0.), KIND=dp)  
    var(ind_pan) = conv  * REAL( MAX(chem(i,k,j,P_pan),0.), KIND=dp)  
    var(ind_pacd) = conv  * REAL( MAX(chem(i,k,j,P_pacd),0.), KIND=dp)  
    var(ind_aacd) = conv  * REAL( MAX(chem(i,k,j,P_aacd),0.), KIND=dp)  
    var(ind_cxo3) = conv  * REAL( MAX(chem(i,k,j,P_cxo3),0.), KIND=dp)  
    var(ind_panx) = conv  * REAL( MAX(chem(i,k,j,P_panx),0.), KIND=dp)  
    var(ind_ror) = conv  * REAL( MAX(chem(i,k,j,P_ror),0.), KIND=dp)  
    var(ind_ole) = conv  * REAL( MAX(chem(i,k,j,P_ole),0.), KIND=dp)  
    var(ind_eth) = conv  * REAL( MAX(chem(i,k,j,P_eth),0.), KIND=dp)  
    var(ind_iole) = conv  * REAL( MAX(chem(i,k,j,P_iole),0.), KIND=dp)  
    var(ind_tol) = conv  * REAL( MAX(chem(i,k,j,P_tol),0.), KIND=dp)  
    var(ind_cres) = conv  * REAL( MAX(chem(i,k,j,P_cres),0.), KIND=dp)  
    var(ind_to2) = conv  * REAL( MAX(chem(i,k,j,P_to2),0.), KIND=dp)  
    var(ind_tolaer1) = conv  * REAL( MAX(chem(i,k,j,P_tolaer1),0.), KIND=dp)  
    var(ind_tolaer2) = conv  * REAL( MAX(chem(i,k,j,P_tolaer2),0.), KIND=dp)  
    var(ind_open) = conv  * REAL( MAX(chem(i,k,j,P_open),0.), KIND=dp)  
    var(ind_cro) = conv  * REAL( MAX(chem(i,k,j,P_cro),0.), KIND=dp)  
    var(ind_cslaer) = conv  * REAL( MAX(chem(i,k,j,P_cslaer),0.), KIND=dp)  
    var(ind_mgly) = conv  * REAL( MAX(chem(i,k,j,P_mgly),0.), KIND=dp)  
    var(ind_xyl) = conv  * REAL( MAX(chem(i,k,j,P_xyl),0.), KIND=dp)  
    var(ind_xylaer1) = conv  * REAL( MAX(chem(i,k,j,P_xylaer1),0.), KIND=dp)  
    var(ind_xylaer2) = conv  * REAL( MAX(chem(i,k,j,P_xylaer2),0.), KIND=dp)  
    var(ind_isop) = conv  * REAL( MAX(chem(i,k,j,P_isop),0.), KIND=dp)  
    var(ind_ispd) = conv  * REAL( MAX(chem(i,k,j,P_ispd),0.), KIND=dp)  
    var(ind_isoaer1) = conv  * REAL( MAX(chem(i,k,j,P_isoaer1),0.), KIND=dp)  
    var(ind_isoaer2) = conv  * REAL( MAX(chem(i,k,j,P_isoaer2),0.), KIND=dp)  
    var(ind_so2) = conv  * REAL( MAX(chem(i,k,j,P_so2),0.), KIND=dp)  
    var(ind_sulf) = conv  * REAL( MAX(chem(i,k,j,P_sulf),0.), KIND=dp)  
    var(ind_sulaer) = conv  * REAL( MAX(chem(i,k,j,P_sulaer),0.), KIND=dp)  
    var(ind_etoh) = conv  * REAL( MAX(chem(i,k,j,P_etoh),0.), KIND=dp)  
    var(ind_etha) = conv  * REAL( MAX(chem(i,k,j,P_etha),0.), KIND=dp)  
    var(ind_terp) = conv  * REAL( MAX(chem(i,k,j,P_terp),0.), KIND=dp)  
    var(ind_terpaer) = conv  * REAL( MAX(chem(i,k,j,P_terpaer),0.), KIND=dp)  
    var(ind_hum) = conv  * REAL( MAX(chem(i,k,j,P_hum),0.), KIND=dp)  
    var(ind_humaer) = conv  * REAL( MAX(chem(i,k,j,P_humaer),0.), KIND=dp)  
    var(ind_lim) = conv  * REAL( MAX(chem(i,k,j,P_lim),0.), KIND=dp)  
    var(ind_limaer1) = conv  * REAL( MAX(chem(i,k,j,P_limaer1),0.), KIND=dp)  
    var(ind_limaer2) = conv  * REAL( MAX(chem(i,k,j,P_limaer2),0.), KIND=dp)  
    var(ind_oci) = conv  * REAL( MAX(chem(i,k,j,P_oci),0.), KIND=dp)  
    var(ind_ociaer1) = conv  * REAL( MAX(chem(i,k,j,P_ociaer1),0.), KIND=dp)  
    var(ind_ociaer2) = conv  * REAL( MAX(chem(i,k,j,P_ociaer2),0.), KIND=dp)  
    var(ind_apin) = conv  * REAL( MAX(chem(i,k,j,P_apin),0.), KIND=dp)  
    var(ind_apinaer1) = conv  * REAL( MAX(chem(i,k,j,P_apinaer1),0.), KIND=dp)  
    var(ind_apinaer2) = conv  * REAL( MAX(chem(i,k,j,P_apinaer2),0.), KIND=dp)  
    var(ind_apinaer3) = conv  * REAL( MAX(chem(i,k,j,P_apinaer3),0.), KIND=dp)  
    var(ind_apinaer4) = conv  * REAL( MAX(chem(i,k,j,P_apinaer4),0.), KIND=dp)  
    var(ind_bpin) = conv  * REAL( MAX(chem(i,k,j,P_bpin),0.), KIND=dp)  
    var(ind_bpinaer1) = conv  * REAL( MAX(chem(i,k,j,P_bpinaer1),0.), KIND=dp)  
    var(ind_bpinaer2) = conv  * REAL( MAX(chem(i,k,j,P_bpinaer2),0.), KIND=dp)  
    var(ind_bpinaer3) = conv  * REAL( MAX(chem(i,k,j,P_bpinaer3),0.), KIND=dp)  
    var(ind_bpinaer4) = conv  * REAL( MAX(chem(i,k,j,P_bpinaer4),0.), KIND=dp)  
    var(ind_bpinaer5) = conv  * REAL( MAX(chem(i,k,j,P_bpinaer5),0.), KIND=dp)  
    var(ind_ter) = conv  * REAL( MAX(chem(i,k,j,P_ter),0.), KIND=dp)  
    var(ind_teraer1) = conv  * REAL( MAX(chem(i,k,j,P_teraer1),0.), KIND=dp)  
    var(ind_teraer2) = conv  * REAL( MAX(chem(i,k,j,P_teraer2),0.), KIND=dp)  
    var(ind_alkh) = conv  * REAL( MAX(chem(i,k,j,P_alkh),0.), KIND=dp)  
    var(ind_alkhaer1) = conv  * REAL( MAX(chem(i,k,j,P_alkhaer1),0.), KIND=dp)  
    var(ind_pah) = conv  * REAL( MAX(chem(i,k,j,P_pah),0.), KIND=dp)  
    var(ind_pahaer1) = conv  * REAL( MAX(chem(i,k,j,P_pahaer1),0.), KIND=dp)  
    var(ind_pahaer2) = conv  * REAL( MAX(chem(i,k,j,P_pahaer2),0.), KIND=dp)  
    var(ind_h2) = conv  * REAL( MAX(chem(i,k,j,P_h2),0.), KIND=dp)  
    var(ind_ch4) = conv  * REAL( MAX(chem(i,k,j,P_ch4),0.), KIND=dp)  
    var(ind_hg0) = conv  * REAL( MAX(chem(i,k,j,P_hg0),0.), KIND=dp)  
    var(ind_hg2) = conv  * REAL( MAX(chem(i,k,j,P_hg2),0.), KIND=dp)  
    var(ind_fmcl) = conv  * REAL( MAX(chem(i,k,j,P_fmcl),0.), KIND=dp)  
    var(ind_cl) = conv  * REAL( MAX(chem(i,k,j,P_cl),0.), KIND=dp)  
    var(ind_hcl) = conv  * REAL( MAX(chem(i,k,j,P_hcl),0.), KIND=dp)  
    var(ind_hocl) = conv  * REAL( MAX(chem(i,k,j,P_hocl),0.), KIND=dp)  
    var(ind_clo) = conv  * REAL( MAX(chem(i,k,j,P_clo),0.), KIND=dp)  
    var(ind_cl2) = conv  * REAL( MAX(chem(i,k,j,P_cl2),0.), KIND=dp)  







   CALL cb05_sorg_aq_Update_Rconst(  &



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








  CALL cb05_sorg_aq_INTEGRATE(TIME_START, TIME_END, &  
          FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK, & 
          ICNTRL_U=icntrl, RCNTRL_U=rcntrl  )









        if(p_nu0.gt.1)then

            rxylho =   ARR2(1.7D-11, -116.0_dp, TEMP); 
            rtolho =   ARR2(1.8D-12, -355.0_dp, TEMP); 
            rcslho =   4.1D-11;                        
            rcslno3 =  2.2D-11;
            rhc8ho =   1.97D-11;                       
            roliho =   ARR2(1.0D-11, -550.0_dp, TEMP); 
            rolino3 =  ARR2(9.6D-13, 270.0_dp, TEMP);  
            rolio3 =   ARR2(8.4D-15, 1100.0_dp, TEMP); 
            roltho =   3.2D-11;                        
            roltno3 =  ARR2(7.0D-13, 2160.0_dp, TEMP); 
            rolto3 =   ARR2(6.5D-15, 1900.0_dp, TEMP); 
            rapiho =   5.37D-11;                       
            rapino3 =  2.31D-12;                       
            rapio3 =   8.66D-17;                       
            rlimho =   1.71D-10;
            rlimno3 =  1.22D-11;
            rlimo3 =   2.00D-16;

            PRDROG(PXYL)  = rxylho * var(ind_xyl)*var(ind_oh)
            PRDROG(PTOL)  = rtolho * var(ind_tol)*var(ind_oh)
            PRDROG(PCSL1) = rcslho * var(ind_cres)*var(ind_oh)
            PRDROG(PCSL2) = rcslno3* var(ind_cres)*var(ind_no3)
            PRDROG(PHC8)  = rhc8ho * var(ind_alkh)*var(ind_oh)
            PRDROG(POLI1) = roliho * var(ind_iole)*var(ind_oh)
            PRDROG(POLI2) = rolino3 * var(ind_iole)*var(ind_no3)
            PRDROG(POLI3) = rolio3 * var(ind_iole)*var(ind_o3)
            PRDROG(POLT1) = roltho * var(ind_ole)*var(ind_oh)
            PRDROG(POLT2) = roltno3 * var(ind_ole)*var(ind_no3)
            PRDROG(POLT3) = rolto3 * var(ind_ole)*var(ind_o3)
            PRDROG(PAPI1) = rapiho * var(ind_apin)*var(ind_oh)
            PRDROG(PAPI2) = rapino3 * var(ind_apin)*var(ind_no3)
            PRDROG(PAPI3) = rapio3 * var(ind_apin)*var(ind_o3)
            PRDROG(PLIM1) = rlimho * var(ind_lim)*var(ind_oh)
            PRDROG(PLIM2) = rlimno3 * var(ind_lim)*var(ind_no3)
            PRDROG(PLIM3) = rlimo3 * var(ind_lim)*var(ind_o3)

            DO n = 1, LDROG
               VDROG3( i,k,j, n ) =  oconv * PRDROG( n ) * DTSTEPC
               VDROG3( i,k,j,n  ) = MAX( 0., VDROG3( i,k,j, n ) )
            ENDDO

        endif


    chem(i,k,j,P_no2) = MAX ( REAL (oconv * var(ind_no2), KIND=sp), 0.)  
    chem(i,k,j,P_no) = MAX ( REAL (oconv * var(ind_no), KIND=sp), 0.)  
    chem(i,k,j,P_o) = MAX ( REAL (oconv * var(ind_o), KIND=sp), 0.)  
    chem(i,k,j,P_o3) = MAX ( REAL (oconv * var(ind_o3), KIND=sp), 0.)  
    chem(i,k,j,P_no3) = MAX ( REAL (oconv * var(ind_no3), KIND=sp), 0.)  
    chem(i,k,j,P_o1d) = MAX ( REAL (oconv * var(ind_o1d), KIND=sp), 0.)  
    chem(i,k,j,P_oh) = MAX ( REAL (oconv * var(ind_oh), KIND=sp), 0.)  
    chem(i,k,j,P_ho2) = MAX ( REAL (oconv * var(ind_ho2), KIND=sp), 0.)  
    chem(i,k,j,P_n2o5) = MAX ( REAL (oconv * var(ind_n2o5), KIND=sp), 0.)  
    chem(i,k,j,P_hno3) = MAX ( REAL (oconv * var(ind_hno3), KIND=sp), 0.)  
    chem(i,k,j,P_hono) = MAX ( REAL (oconv * var(ind_hono), KIND=sp), 0.)  
    chem(i,k,j,P_pna) = MAX ( REAL (oconv * var(ind_pna), KIND=sp), 0.)  
    chem(i,k,j,P_h2o2) = MAX ( REAL (oconv * var(ind_h2o2), KIND=sp), 0.)  
    chem(i,k,j,P_xo2) = MAX ( REAL (oconv * var(ind_xo2), KIND=sp), 0.)  
    chem(i,k,j,P_xo2n) = MAX ( REAL (oconv * var(ind_xo2n), KIND=sp), 0.)  
    chem(i,k,j,P_ntr) = MAX ( REAL (oconv * var(ind_ntr), KIND=sp), 0.)  
    chem(i,k,j,P_rooh) = MAX ( REAL (oconv * var(ind_rooh), KIND=sp), 0.)  
    chem(i,k,j,P_form) = MAX ( REAL (oconv * var(ind_form), KIND=sp), 0.)  
    chem(i,k,j,P_ald2) = MAX ( REAL (oconv * var(ind_ald2), KIND=sp), 0.)  
    chem(i,k,j,P_aldx) = MAX ( REAL (oconv * var(ind_aldx), KIND=sp), 0.)  
    chem(i,k,j,P_par) = MAX ( REAL (oconv * var(ind_par), KIND=sp), 0.)  
    chem(i,k,j,P_co) = MAX ( REAL (oconv * var(ind_co), KIND=sp), 0.)  
    chem(i,k,j,P_meo2) = MAX ( REAL (oconv * var(ind_meo2), KIND=sp), 0.)  
    chem(i,k,j,P_mepx) = MAX ( REAL (oconv * var(ind_mepx), KIND=sp), 0.)  
    chem(i,k,j,P_meoh) = MAX ( REAL (oconv * var(ind_meoh), KIND=sp), 0.)  
    chem(i,k,j,P_hco3) = MAX ( REAL (oconv * var(ind_hco3), KIND=sp), 0.)  
    chem(i,k,j,P_facd) = MAX ( REAL (oconv * var(ind_facd), KIND=sp), 0.)  
    chem(i,k,j,P_c2o3) = MAX ( REAL (oconv * var(ind_c2o3), KIND=sp), 0.)  
    chem(i,k,j,P_pan) = MAX ( REAL (oconv * var(ind_pan), KIND=sp), 0.)  
    chem(i,k,j,P_pacd) = MAX ( REAL (oconv * var(ind_pacd), KIND=sp), 0.)  
    chem(i,k,j,P_aacd) = MAX ( REAL (oconv * var(ind_aacd), KIND=sp), 0.)  
    chem(i,k,j,P_cxo3) = MAX ( REAL (oconv * var(ind_cxo3), KIND=sp), 0.)  
    chem(i,k,j,P_panx) = MAX ( REAL (oconv * var(ind_panx), KIND=sp), 0.)  
    chem(i,k,j,P_ror) = MAX ( REAL (oconv * var(ind_ror), KIND=sp), 0.)  
    chem(i,k,j,P_ole) = MAX ( REAL (oconv * var(ind_ole), KIND=sp), 0.)  
    chem(i,k,j,P_eth) = MAX ( REAL (oconv * var(ind_eth), KIND=sp), 0.)  
    chem(i,k,j,P_iole) = MAX ( REAL (oconv * var(ind_iole), KIND=sp), 0.)  
    chem(i,k,j,P_tol) = MAX ( REAL (oconv * var(ind_tol), KIND=sp), 0.)  
    chem(i,k,j,P_cres) = MAX ( REAL (oconv * var(ind_cres), KIND=sp), 0.)  
    chem(i,k,j,P_to2) = MAX ( REAL (oconv * var(ind_to2), KIND=sp), 0.)  
    chem(i,k,j,P_tolaer1) = MAX ( REAL (oconv * var(ind_tolaer1), KIND=sp), 0.)  
    chem(i,k,j,P_tolaer2) = MAX ( REAL (oconv * var(ind_tolaer2), KIND=sp), 0.)  
    chem(i,k,j,P_open) = MAX ( REAL (oconv * var(ind_open), KIND=sp), 0.)  
    chem(i,k,j,P_cro) = MAX ( REAL (oconv * var(ind_cro), KIND=sp), 0.)  
    chem(i,k,j,P_cslaer) = MAX ( REAL (oconv * var(ind_cslaer), KIND=sp), 0.)  
    chem(i,k,j,P_mgly) = MAX ( REAL (oconv * var(ind_mgly), KIND=sp), 0.)  
    chem(i,k,j,P_xyl) = MAX ( REAL (oconv * var(ind_xyl), KIND=sp), 0.)  
    chem(i,k,j,P_xylaer1) = MAX ( REAL (oconv * var(ind_xylaer1), KIND=sp), 0.)  
    chem(i,k,j,P_xylaer2) = MAX ( REAL (oconv * var(ind_xylaer2), KIND=sp), 0.)  
    chem(i,k,j,P_isop) = MAX ( REAL (oconv * var(ind_isop), KIND=sp), 0.)  
    chem(i,k,j,P_ispd) = MAX ( REAL (oconv * var(ind_ispd), KIND=sp), 0.)  
    chem(i,k,j,P_isoaer1) = MAX ( REAL (oconv * var(ind_isoaer1), KIND=sp), 0.)  
    chem(i,k,j,P_isoaer2) = MAX ( REAL (oconv * var(ind_isoaer2), KIND=sp), 0.)  
    chem(i,k,j,P_so2) = MAX ( REAL (oconv * var(ind_so2), KIND=sp), 0.)  
    chem(i,k,j,P_sulf) = MAX ( REAL (oconv * var(ind_sulf), KIND=sp), 0.)  
    chem(i,k,j,P_sulaer) = MAX ( REAL (oconv * var(ind_sulaer), KIND=sp), 0.)  
    chem(i,k,j,P_etoh) = MAX ( REAL (oconv * var(ind_etoh), KIND=sp), 0.)  
    chem(i,k,j,P_etha) = MAX ( REAL (oconv * var(ind_etha), KIND=sp), 0.)  
    chem(i,k,j,P_terp) = MAX ( REAL (oconv * var(ind_terp), KIND=sp), 0.)  
    chem(i,k,j,P_terpaer) = MAX ( REAL (oconv * var(ind_terpaer), KIND=sp), 0.)  
    chem(i,k,j,P_hum) = MAX ( REAL (oconv * var(ind_hum), KIND=sp), 0.)  
    chem(i,k,j,P_humaer) = MAX ( REAL (oconv * var(ind_humaer), KIND=sp), 0.)  
    chem(i,k,j,P_lim) = MAX ( REAL (oconv * var(ind_lim), KIND=sp), 0.)  
    chem(i,k,j,P_limaer1) = MAX ( REAL (oconv * var(ind_limaer1), KIND=sp), 0.)  
    chem(i,k,j,P_limaer2) = MAX ( REAL (oconv * var(ind_limaer2), KIND=sp), 0.)  
    chem(i,k,j,P_oci) = MAX ( REAL (oconv * var(ind_oci), KIND=sp), 0.)  
    chem(i,k,j,P_ociaer1) = MAX ( REAL (oconv * var(ind_ociaer1), KIND=sp), 0.)  
    chem(i,k,j,P_ociaer2) = MAX ( REAL (oconv * var(ind_ociaer2), KIND=sp), 0.)  
    chem(i,k,j,P_apin) = MAX ( REAL (oconv * var(ind_apin), KIND=sp), 0.)  
    chem(i,k,j,P_apinaer1) = MAX ( REAL (oconv * var(ind_apinaer1), KIND=sp), 0.)  
    chem(i,k,j,P_apinaer2) = MAX ( REAL (oconv * var(ind_apinaer2), KIND=sp), 0.)  
    chem(i,k,j,P_apinaer3) = MAX ( REAL (oconv * var(ind_apinaer3), KIND=sp), 0.)  
    chem(i,k,j,P_apinaer4) = MAX ( REAL (oconv * var(ind_apinaer4), KIND=sp), 0.)  
    chem(i,k,j,P_bpin) = MAX ( REAL (oconv * var(ind_bpin), KIND=sp), 0.)  
    chem(i,k,j,P_bpinaer1) = MAX ( REAL (oconv * var(ind_bpinaer1), KIND=sp), 0.)  
    chem(i,k,j,P_bpinaer2) = MAX ( REAL (oconv * var(ind_bpinaer2), KIND=sp), 0.)  
    chem(i,k,j,P_bpinaer3) = MAX ( REAL (oconv * var(ind_bpinaer3), KIND=sp), 0.)  
    chem(i,k,j,P_bpinaer4) = MAX ( REAL (oconv * var(ind_bpinaer4), KIND=sp), 0.)  
    chem(i,k,j,P_bpinaer5) = MAX ( REAL (oconv * var(ind_bpinaer5), KIND=sp), 0.)  
    chem(i,k,j,P_ter) = MAX ( REAL (oconv * var(ind_ter), KIND=sp), 0.)  
    chem(i,k,j,P_teraer1) = MAX ( REAL (oconv * var(ind_teraer1), KIND=sp), 0.)  
    chem(i,k,j,P_teraer2) = MAX ( REAL (oconv * var(ind_teraer2), KIND=sp), 0.)  
    chem(i,k,j,P_alkh) = MAX ( REAL (oconv * var(ind_alkh), KIND=sp), 0.)  
    chem(i,k,j,P_alkhaer1) = MAX ( REAL (oconv * var(ind_alkhaer1), KIND=sp), 0.)  
    chem(i,k,j,P_pah) = MAX ( REAL (oconv * var(ind_pah), KIND=sp), 0.)  
    chem(i,k,j,P_pahaer1) = MAX ( REAL (oconv * var(ind_pahaer1), KIND=sp), 0.)  
    chem(i,k,j,P_pahaer2) = MAX ( REAL (oconv * var(ind_pahaer2), KIND=sp), 0.)  
    chem(i,k,j,P_h2) = MAX ( REAL (oconv * var(ind_h2), KIND=sp), 0.)  
    chem(i,k,j,P_ch4) = MAX ( REAL (oconv * var(ind_ch4), KIND=sp), 0.)  
    chem(i,k,j,P_hg0) = MAX ( REAL (oconv * var(ind_hg0), KIND=sp), 0.)  
    chem(i,k,j,P_hg2) = MAX ( REAL (oconv * var(ind_hg2), KIND=sp), 0.)  
    chem(i,k,j,P_fmcl) = MAX ( REAL (oconv * var(ind_fmcl), KIND=sp), 0.)  
    chem(i,k,j,P_cl) = MAX ( REAL (oconv * var(ind_cl), KIND=sp), 0.)  
    chem(i,k,j,P_hcl) = MAX ( REAL (oconv * var(ind_hcl), KIND=sp), 0.)  
    chem(i,k,j,P_hocl) = MAX ( REAL (oconv * var(ind_hocl), KIND=sp), 0.)  
    chem(i,k,j,P_clo) = MAX ( REAL (oconv * var(ind_clo), KIND=sp), 0.)  
    chem(i,k,j,P_cl2) = MAX ( REAL (oconv * var(ind_cl2), KIND=sp), 0.)  



    END DO
    END DO
    END DO









END SUBROUTINE  cb05_sorg_aq_interface


END MODULE module_kpp_cb05_sorg_aq_interf 




