MODULE module_bioemi_megan2

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

CONTAINS

  SUBROUTINE bio_emissions_megan2(id,config_flags,ktau,dtstep,         &
       curr_secs,julday,gmt,xlat,xlong,p_phy,rho_phy,dz8w,             &
       chem, ne_area,                                                  &
       current_month,                                                  &
       T2,swdown,                                                      &
       nmegan, EFmegan, msebio_isop,                                   &
       mlai,                                                           &
       pftp_bt, pftp_nt, pftp_sb, pftp_hb,                             &
       mtsa,                                                           &
       mswdown,                                                        &
       mebio_isop, mebio_apin, mebio_bpin, mebio_bcar,                 &
       mebio_acet, mebio_mbo, mebio_no,                                &
       ebio_iso,ebio_oli,ebio_api,ebio_lim,                            &
       ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,                   &
       ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_no,                   &
       ebio_c10h16,ebio_tol,ebio_bigalk, ebio_ch3oh,ebio_acet,         &
       ebio_nh3,ebio_no2,ebio_c2h5oh,ebio_ch3cooh,ebio_mek,            &
       ebio_bigene,ebio_c2h6,ebio_c2h4,ebio_c3h6,ebio_c3h8,ebio_so2,   &
       ebio_dms,ebio_hcn,                                              &
       ebio_c5h8,ebio_apinene,ebio_bpinene,ebio_toluene,               &
       ebio_ch3cho,ebio_ch3co2h,ebio_tbut2ene,ebio_c2h5cho, &
       ebio_nc4h10, &
       ebio_sesq, ebio_mbo, ebio_bpi, ebio_myrc,                       &
       ebio_alk3, ebio_alk4, ebio_alk5, ebio_ole1, ebio_ole2,          &    
       ebio_aro1, ebio_aro2, ebio_ccho, ebio_meoh,                     &    
       ebio_ethene, ebio_hcooh, ebio_terp, ebio_bald,                  &    
       ebio_cco_oh, ebio_rco_oh,                                       &    
       e_bio,                                                          &
       ids,ide, jds,jde, kds,kde,                                      &
       ims,ime, jms,jme, kms,kme,                                      &
       its,ite, jts,jte, kts,kte                                       )

    USE module_configure
    USE module_state_description
    USE module_data_megan2
    USE module_data_mgn2mech


    IMPLICIT NONE

    

    
    TYPE(grid_config_rec_type),  INTENT(IN)    :: config_flags

    
    INTEGER,   INTENT(IN   ) :: id,ktau,                               &
         ids,ide, jds,jde, kds,kde,                                    &
         ims,ime, jms,jme, kms,kme,                                    &
         its,ite, jts,jte, kts,kte

    
    INTEGER, INTENT (IN) :: julday   
    
    REAL, INTENT(IN) :: gmt,dtstep

    
    REAL(KIND=8), INTENT(IN) :: curr_secs

    
    INTEGER, INTENT (IN) :: ne_area

    
    REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                   &
         INTENT(IN) :: p_phy

    
    REAL,  DIMENSION( ims:ime , jms:jme ),                             &
         INTENT(IN) :: xlat, xlong

    
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                      &
         INTENT(IN) :: rho_phy

    
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,           &
          INTENT(IN) :: dz8w

    
    REAL,  DIMENSION( ims:ime , jms:jme ),                             &
         INTENT(IN) :: T2

    
    REAL,  DIMENSION( ims:ime , jms:jme ),                             &
         INTENT(IN) :: swdown                                    

    
    
    
    INTEGER, INTENT(IN) :: nmegan
    
    
    REAL, DIMENSION (ims:ime, jms:jme , nmegan) ,                      &
         INTENT(INOUT) :: EFmegan

    
    
    
    REAL,  DIMENSION( ims:ime , jms:jme ),                             &
         INTENT(IN ) :: msebio_isop

    
    
    REAL, DIMENSION ( ims:ime , jms:jme ),                             &
         INTENT(IN) ::                                                 &
         pftp_bt, pftp_nt, pftp_sb, pftp_hb

    
    
    REAL,  DIMENSION( ims:ime , jms:jme , 12 ),                        &
         INTENT(IN) :: mlai

    
    
    REAL,  DIMENSION( ims:ime , jms:jme , 12 ),                        &
         INTENT(IN) :: mtsa

    
    
    REAL, DIMENSION ( ims:ime , jms:jme , 12 ),                        &
         INTENT(IN) :: mswdown

    
    
    
    REAL,  DIMENSION( ims:ime , jms:jme ),                             &
         INTENT(INOUT) ::                                              &
         mebio_isop, mebio_apin, mebio_bpin, mebio_bcar,               &
         mebio_acet, mebio_mbo, mebio_no

    
    
   REAL, DIMENSION( ims:ime, jms:jme, ne_area ),                       &
         INTENT(INOUT ) :: e_bio

    
    
    
    
    REAL,  DIMENSION( ims:ime , jms:jme ),                             &
         INTENT(INOUT  ) ::                                            &
         ebio_iso,ebio_oli,ebio_api,ebio_lim,                          &
         ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,                 &
         ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_no,                 &
         ebio_c10h16,ebio_tol,ebio_bigalk, ebio_ch3oh,ebio_acet,       &
         ebio_nh3,ebio_no2,ebio_c2h5oh,ebio_ch3cooh,ebio_mek,          &
         ebio_bigene,ebio_c2h6,ebio_c2h4,ebio_c3h6,ebio_c3h8,ebio_so2, &
         ebio_dms,ebio_hcn,                                            &
         ebio_c5h8,ebio_apinene,ebio_bpinene,ebio_toluene,             &
         ebio_ch3cho,ebio_ch3co2h,ebio_tbut2ene,ebio_c2h5cho,          &
         ebio_nc4h10,                                                  &
         ebio_sesq,ebio_mbo,ebio_bpi,ebio_myrc,                        &
         ebio_alk3, ebio_alk4, ebio_alk5, ebio_ole1, ebio_ole2,        &    
         ebio_aro1, ebio_aro2, ebio_ccho, ebio_meoh,                   &    
         ebio_ethene, ebio_hcooh, ebio_terp, ebio_bald,                &    
         ebio_cco_oh, ebio_rco_oh

    
    
    
    
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),            &
         INTENT(INOUT) :: chem

    
    INTEGER, INTENT(IN) :: current_month

    

    
    REAL, PARAMETER :: min_emis = 0.001
    

    
    INTEGER, PARAMETER :: DaysInMonth(12) = (/   &
         31,28,31,30,31,30,31,31,30,31,30,31 /)
    
    REAL, PARAMETER :: PI = 3.14159 
    REAL, PARAMETER :: D2RAD = PI/180.0 


    

    CHARACTER(len=256)   ::   mesg
    INTEGER :: i,j,k,i_class, i_spc, icount, p_in_chem
    INTEGER :: previous_month

    
    
    REAL(KIND=8) :: xtime
 
    
    
    INTEGER :: ixhour
    REAL(KIND=8) :: xhour

    
    
    REAL :: xmin

    
    
    REAL :: gmtp

    
    
    REAL :: tmidh

    
    REAL :: LAIc, LAIp

    
    REAL :: tsa, tsa24, pres

    
    REAL :: lat, lon

    
    REAL :: swd, swd24
    
    
    REAL :: par, par24, pardb, pardif

    
    REAL :: zen , coszen

    
    REAL :: tstlen

    
    REAL :: epsilon

    
    
    
    REAL :: gam_LHT, gam_TMP, gam_PHO,gam_AGE, gam_SMT

    
    
    REAL :: rho

    
    REAL :: ldf

    
    REAL :: convert2

    
    REAL :: gas_emis

    

    
    
    
    
    REAL, DIMENSION(n_mgn_spc) :: adjust_factor

    
    REAL :: pft_frac(n_pft)

    
    
    REAL, DIMENSION ( n_spca_spc ) :: E_megan2

    


    
    
    
    
    
    IF ( ktau .EQ. 1 ) THEN

       
       
       
       IF ( nmegan .NE. n_spca_spc ) THEN
          WRITE(mesg,*)'namelist variable nmegan does not match n_spca_spc'
          CALL wrf_error_fatal3("<stdin>",307,&
mesg)          
       END IF

       
       
       
       IF ( (imgn_isop .NE. 1) .OR. (is_isoprene .NE. 1) ) THEN
          WRITE(mesg,*)'imgn_isop and is_isoprene in bio_emissions_megan should be 1'
          CALL wrf_error_fatal3("<stdin>",316,&
mesg)          
       END IF

    END IF


    
    ebio_iso  ( its:ite , jts:jte ) = 0.0
    ebio_oli  ( its:ite , jts:jte ) = 0.0
    ebio_api  ( its:ite , jts:jte ) = 0.0
    ebio_lim  ( its:ite , jts:jte ) = 0.0
    ebio_hc3  ( its:ite , jts:jte ) = 0.0
    ebio_ete  ( its:ite , jts:jte ) = 0.0
    ebio_olt  ( its:ite , jts:jte ) = 0.0
    ebio_ket  ( its:ite , jts:jte ) = 0.0
    ebio_ald  ( its:ite , jts:jte ) = 0.0
    ebio_hcho ( its:ite , jts:jte ) = 0.0
    ebio_eth  ( its:ite , jts:jte ) = 0.0
    ebio_ora2 ( its:ite , jts:jte ) = 0.0
    ebio_co   ( its:ite , jts:jte ) = 0.0
    ebio_no   ( its:ite , jts:jte ) = 0.0
    ebio_c10h16( its:ite , jts:jte ) = 0.0
    ebio_tol  ( its:ite , jts:jte ) = 0.0
    ebio_bigalk( its:ite , jts:jte ) = 0.0
    ebio_ch3oh ( its:ite , jts:jte ) = 0.0
    ebio_acet  ( its:ite , jts:jte ) = 0.0
    ebio_nh3   ( its:ite , jts:jte ) = 0.0
    ebio_no2   ( its:ite , jts:jte ) = 0.0
    ebio_c2h5oh( its:ite , jts:jte ) = 0.0
    ebio_ch3cooh( its:ite , jts:jte ) = 0.0
    ebio_mek   ( its:ite , jts:jte ) = 0.0
    ebio_bigene( its:ite , jts:jte ) = 0.0
    ebio_c2h4  ( its:ite , jts:jte ) = 0.0
    ebio_c2h6  ( its:ite , jts:jte ) = 0.0
    ebio_c3h6  ( its:ite , jts:jte ) = 0.0
    ebio_c3h8  ( its:ite , jts:jte ) = 0.0
    ebio_so2   ( its:ite , jts:jte ) = 0.0
    ebio_dms   ( its:ite , jts:jte ) = 0.0
    ebio_terp  ( its:ite , jts:jte ) = 0.0
    ebio_c5h8   ( its:ite , jts:jte ) = 0.0
    ebio_apinene   ( its:ite , jts:jte ) = 0.0
    ebio_bpinene   ( its:ite , jts:jte ) = 0.0
    ebio_toluene   ( its:ite , jts:jte ) = 0.0
    ebio_hcooh   ( its:ite , jts:jte ) = 0.0
    ebio_ch3cho   ( its:ite , jts:jte ) = 0.0
    ebio_c2h5oh   ( its:ite , jts:jte ) = 0.0
    ebio_ch3co2h   ( its:ite , jts:jte ) = 0.0
    ebio_tbut2ene   ( its:ite , jts:jte ) = 0.0
    ebio_c2h5cho   ( its:ite , jts:jte ) = 0.0
    ebio_nc4h10   ( its:ite , jts:jte ) = 0.0
    ebio_sesq  ( its:ite , jts:jte ) = 0.0
    ebio_mbo   ( its:ite , jts:jte ) = 0.0
    ebio_bpi   ( its:ite , jts:jte ) = 0.0
    ebio_myrc  ( its:ite , jts:jte ) = 0.0
    e_bio     ( its:ite , jts:jte , 1:ne_area) = 0.0
    
    
    
    mebio_isop ( its:ite , jts:jte ) = 0.0
    mebio_apin ( its:ite , jts:jte ) = 0.0
    mebio_bpin ( its:ite , jts:jte ) = 0.0
    mebio_bcar ( its:ite , jts:jte ) = 0.0
    mebio_acet ( its:ite , jts:jte ) = 0.0 
    mebio_mbo  ( its:ite , jts:jte ) = 0.0
    mebio_no   ( its:ite , jts:jte ) = 0.0


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    IF (current_month > 1) THEN
       previous_month = current_month -1
    ELSE
       previous_month = 12
    END IF


    
    
    
    
    
    

    
    
    

    xtime = curr_secs/60._8 + real(dtstep/120.,8)
    
    
    ixhour = int(gmt + 0.01) + int(xtime/60._8)
    xhour=real(ixhour,8)
    
    
    xmin = 60.*gmt + real(xtime-xhour*60._8,8)
    
    
    gmtp=MOD(xhour,24._8)
    
    
    tmidh= gmtp + xmin/60.

    WRITE(mesg,*) 'calculate MEGAN emissions at ktau, gmtp, tmidh = ',ktau, gmtp, tmidh
    CALL wrf_message(mesg)


    
    
    
    
    
    
    GAS_MECH_SELECT1: SELECT CASE (config_flags%chem_opt)
             
    CASE (RADM2, RADM2_KPP, RADM2SORG, RADM2SORG_AQ, RADM2SORG_AQCHEM, RADM2SORG_KPP,GOCARTRADM2)
       
       CALL get_megan2radm2_table

    CASE (RACMSORG_AQ, RACMSORG_AQCHEM_KPP, RACM_ESRLSORG_AQCHEM_KPP, RACM_ESRLSORG_KPP, RACM_KPP, GOCARTRACM_KPP, RACMSORG_KPP, &
          RACM_MIM_KPP, RACMPM_KPP)
             
       
       CALL get_megan2racm_table
    CASE (RACM_SOA_VBS_KPP,RACM_SOA_VBS_AQCHEM_KPP,RACM_SOA_VBS_HET_KPP)

        
        CALL get_megan2racmSOA_table

    CASE (CBMZ, CBMZ_BB, CBMZ_BB_KPP, CBMZ_MOSAIC_KPP, &
          CBMZ_MOSAIC_4BIN, & 
          CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
          CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, &
          CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, CBMZSORG, CBMZSORG_AQ, &
          CBMZ_CAM_MAM3_NOAQ, CBMZ_CAM_MAM3_AQ, CBMZ_CAM_MAM7_NOAQ, CBMZ_CAM_MAM7_AQ)
        
       
       CALL get_megan2cbmz_table

    CASE (CB05_SORG_AQ_KPP)
       CALL get_megan2cb05_table

    CASE ( CB05_SORG_VBS_AQ_KPP)
       CALL get_megan2cb05vbs_table

    CASE ( MOZART_KPP, MOZCART_KPP )
       
       CALL get_megan2mozcart_table
    CASE (  T1_MOZCART_KPP )
       CALL get_megan2t1_mozc_table
    CASE (  MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP )
       CALL get_megan2mozm_table

    CASE (SAPRC99_KPP,SAPRC99_MOSAIC_4BIN_VBS2_KPP, &
         SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_8BIN_VBS2_KPP)
       CALL get_megan2saprcnov_table

    CASE ( CRIMECH_KPP, CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP )
       
       CALL get_megan2crimech_table

    CASE DEFAULT
       
       CALL wrf_error_fatal3("<stdin>",493,&
'Species conversion table for MEGAN v2.04 not available. ')

    END SELECT GAS_MECH_SELECT1

    

    j_loop: DO j = jts, jte
       i_loop: DO i = its, ite


          

          tsa   = T2(i,j)                     
          pres  = 0.01*p_phy(i,kts,j)         
          lat   = xlat(i,j)                   
          lon   = xlong(i,j)                  
          swd   = swdown(i,j)                 
          LAIc  = mlai(i,j,current_month)     
          LAIp  = mlai(i,j,previous_month)    

          
          
          
          tsa24 = mtsa    (i,j,previous_month) 
          swd24 = mswdown (i,j,previous_month) 

          
          IF (tsa .LT. 200.0) THEN
             WRITE (mesg,'("temperature too low at i=",i3," ,j=",i3 )')i,j
             CALL wrf_message(mesg)
          END IF
          IF (tsa .GT. 315.0 ) THEN
             WRITE (mesg,'("temperature too high at i=",i3," ,j=",i3," ;resetting to 315K" )')i,j
             CALL wrf_message(mesg)
             tsa = 315.0
          END IF







          
          





          
          
          par = 4.766 * 0.5 * swd

          
          IF ( par .LT. 0.00 .OR. par .GT. 2600.0 ) THEN
             WRITE (mesg,'("par out of range at i=",i3," ,j=",i3," par =",f8.3 )')i,j,par
             CALL wrf_message(mesg)
          END IF

          
          
          
          par24 = swd24 * 4.766 * 0.5

          
          
          
          
          
          
          
          
          
          
          
          
          
          
          

          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          


          
          
          
          
          
          
          
          
          
          
          

          
          
          
          
          
          
          
          
          
          

          

          
          CALL GAMMA_LAI( LAIc, gam_LHT)

          
          CALL GAMMA_P( julday, tmidh, lat, lon, par, par24, gam_PHO )

          
          
          gam_SMT = 1.0

          
          

          DO i_class = 1, n_mgn_spc
             
             
             

             IF ( i_class == imgn_isop ) THEN
                CALL GAMMA_TISOP( tsa, tsa24, gam_TMP )
             ELSE
                CALL GAMMA_TNISP( i_class , tsa, gam_TMP  )
             END IF

             

             
             
             
             tstlen = REAL(DaysInMonth(previous_month),KIND(1.0))

             CALL GAMMA_A( i_class , LAIp, LAIc, TSTLEN, tsa24, gam_AGE )

             
             
             
             
             rho = rho_fct(i_class)

             
             
             
             ldf = ldf_fct(i_class)

             
             

             adjust_factor(i_class) = gam_TMP * gam_AGE * gam_LHT * gam_SMT * rho * &
                  ( (1.0-LDF) + gam_PHO*LDF )

          END DO 


          
          
          
          
          E_megan2(is_isoprene) = adjust_factor(imgn_isop)*msebio_isop(i,j)
          IF ( E_megan2(is_isoprene) .LT. min_emis ) E_megan2(is_isoprene)=0.


          
          
          
          

          
          
          
          
          DO i_spc = 2, n_spca_spc 

             
             i_class = mg20_map (i_spc)

             
             
             
             
             
             


                
                
                
                pft_frac(k_bt) = 0.01*pftp_bt(i,j)
                pft_frac(k_nt) = 0.01*pftp_nt(i,j)
                pft_frac(k_sb) = 0.01*pftp_sb(i,j)
                pft_frac(k_hb) = 0.01*pftp_hb(i,j)

                
                epsilon = 0.0
                DO k = 1, n_pft 
                   epsilon = epsilon +                             &
                        pft_frac(k)*EF(i_class,k)*EF_frac(i_spc,k)
                END DO

                
                
                
                EFmegan(i,j,i_spc) = epsilon



             
             
             
             E_megan2(i_spc) = EFmegan(i,j,i_spc)*        &
                  adjust_factor(i_class)/spca_mwt(i_spc)
             IF ( E_megan2(i_spc) .LT. min_emis ) E_megan2(i_spc)=0.

          END DO 


          
          







          mebio_isop  (i,j) = E_megan2 ( is_isoprene        )
          mebio_apin  (i,j) = E_megan2 ( is_pinene_a        )
          mebio_bpin  (i,j) = E_megan2 ( is_pinene_b        )
          mebio_bcar  (i,j) = E_megan2 ( is_caryophyllene_b )
          mebio_acet  (i,j) = E_megan2 ( is_acetone         )
          mebio_mbo   (i,j) = E_megan2 ( is_MBO_2m3e2ol     )
          mebio_no    (i,j) = E_megan2 ( is_nitric_OXD      )


          
          


          
          
          convert2 = 0.02897/(rho_phy(i,kts,j)*60.)


          
          GAS_MECH_SELECT: SELECT CASE (config_flags%chem_opt)

          CASE ( MOZART_KPP, MOZCART_KPP )

             DO icount = 1, n_megan2mozcart



                p_in_chem = p_of_mozcart(icount)
use_megan_emission : &
                IF ( p_in_chem /= non_react ) THEN



is_mozcart_species : &
                   IF ( p_in_chem >= param_first_scalar ) THEN



                      gas_emis = mozcart_per_megan(icount) * E_megan2(p_of_megan2mozcart(icount))






                      IF ( p_in_chem == p_isopr ) THEN
                         ebio_iso(i,j) = ebio_iso(i,j) + gas_emis
                         e_bio(i,j,p_isopr-1)   = e_bio(i,j,p_isopr-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_no ) THEN
                         ebio_no(i,j)  = ebio_no(i,j) + gas_emis
                         e_bio(i,j,p_no-1)   = e_bio(i,j,p_no-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_no2 ) THEN
                         ebio_no2(i,j)  = ebio_no2(i,j) + gas_emis
                         e_bio(i,j,p_no2-1)   = e_bio(i,j,p_no2-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_co ) THEN
                         ebio_co(i,j)  = ebio_co(i,j) + gas_emis
                         e_bio(i,j,p_co-1)   = e_bio(i,j,p_co-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_hcho ) THEN
                         ebio_hcho(i,j) = ebio_hcho(i,j) + gas_emis
                         e_bio(i,j,p_hcho-1)   = e_bio(i,j,p_hcho-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ald ) THEN
                         ebio_ald(i,j) = ebio_ald(i,j) + gas_emis
                         e_bio(i,j,p_ald-1)   = e_bio(i,j,p_ald-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_acet ) THEN
                         ebio_acet(i,j) = ebio_acet(i,j) + gas_emis
                         e_bio(i,j,p_acet-1)   = e_bio(i,j,p_acet-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_tol ) THEN
                         ebio_tol(i,j) = ebio_tol(i,j) + gas_emis
                         e_bio(i,j,p_tol-1)   = e_bio(i,j,p_tol-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c10h16 ) THEN
                         ebio_c10h16(i,j) = ebio_c10h16(i,j) + gas_emis
                         e_bio(i,j,p_c10h16-1)   = e_bio(i,j,p_c10h16-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_so2 ) THEN
                         ebio_so2(i,j) = ebio_so2(i,j) + gas_emis
                         e_bio(i,j,p_so2-1)   = e_bio(i,j,p_so2-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_dms ) THEN
                         ebio_dms(i,j) = ebio_dms(i,j) + gas_emis
                         e_bio(i,j,p_dms-1)   = e_bio(i,j,p_dms-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bigalk ) THEN
                         ebio_bigalk(i,j) = ebio_bigalk(i,j) + gas_emis
                         e_bio(i,j,p_bigalk-1)   = e_bio(i,j,p_bigalk-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bigene ) THEN
                         ebio_bigene(i,j) = ebio_bigene(i,j) + gas_emis
                         e_bio(i,j,p_bigene-1)   = e_bio(i,j,p_bigene-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_nh3 ) THEN
                         ebio_nh3(i,j) = ebio_nh3(i,j) + gas_emis
                         e_bio(i,j,p_nh3-1)   = e_bio(i,j,p_nh3-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ch3oh ) THEN
                         ebio_ch3oh(i,j) = ebio_ch3oh(i,j) + gas_emis
                         e_bio(i,j,p_ch3oh-1)   = e_bio(i,j,p_ch3oh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h5oh ) THEN
                         ebio_c2h5oh(i,j) = ebio_c2h5oh(i,j) + gas_emis
                         e_bio(i,j,p_c2h5oh-1)   = e_bio(i,j,p_c2h5oh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ch3cooh ) THEN
                         ebio_ch3cooh(i,j) = ebio_ch3cooh(i,j) + gas_emis
                         e_bio(i,j,p_ch3cooh-1)   = e_bio(i,j,p_ch3cooh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_mek ) THEN
                         ebio_mek(i,j) = ebio_mek(i,j) + gas_emis
                         e_bio(i,j,p_mek-1)   = e_bio(i,j,p_mek-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h4 ) THEN
                         ebio_c2h4(i,j) = ebio_c2h4(i,j) + gas_emis
                         e_bio(i,j,p_c2h4-1)   = e_bio(i,j,p_c2h4-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h6 ) THEN
                         ebio_c2h6(i,j) = ebio_c2h6(i,j) + gas_emis
                         e_bio(i,j,p_c2h6-1)   = e_bio(i,j,p_c2h6-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c3h6 ) THEN
                         ebio_c3h6(i,j) = ebio_c3h6(i,j) + gas_emis
                         e_bio(i,j,p_c3h6-1)   = e_bio(i,j,p_c3h6-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c3h8 ) THEN
                         ebio_c3h8(i,j) = ebio_c3h8(i,j) + gas_emis
                         e_bio(i,j,p_c3h8-1)   = e_bio(i,j,p_c3h8-1)  + gas_emis*convert2
                      END IF
                   END IF is_mozcart_species
                END IF use_megan_emission
             END DO

          CASE ( T1_MOZCART_KPP )

             DO icount = 1, n_megan2t1_mozc



                p_in_chem = p_of_t1_mozc(icount)
use_megan_emis_a : &
                IF ( p_in_chem /= non_react ) THEN



is_t1_mozc_species : &
                   IF ( p_in_chem >= param_first_scalar ) THEN



                      gas_emis = t1_mozc_per_megan(icount) * E_megan2(p_of_megan2t1_mozc(icount))






                      IF ( p_in_chem == p_isopr ) THEN
                         ebio_iso(i,j) = ebio_iso(i,j) + gas_emis
                         e_bio(i,j,p_isopr-1)   = e_bio(i,j,p_isopr-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_apin ) THEN
                         ebio_api(i,j) = ebio_api(i,j) + gas_emis
                         e_bio(i,j,p_apin-1)   = e_bio(i,j,p_apin-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bpin ) THEN
                         ebio_bpi(i,j) = ebio_bpi(i,j) + gas_emis
                         e_bio(i,j,p_bpin-1)   = e_bio(i,j,p_bpin-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_limon ) THEN
                         ebio_lim(i,j) = ebio_lim(i,j) + gas_emis
                         e_bio(i,j,p_limon-1)   = e_bio(i,j,p_limon-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_myrc ) THEN
                         ebio_myrc(i,j) = ebio_myrc(i,j) + gas_emis
                         e_bio(i,j,p_myrc-1)   = e_bio(i,j,p_myrc-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bcary ) THEN
                         ebio_sesq(i,j) = ebio_sesq(i,j) + gas_emis
                         e_bio(i,j,p_bcary-1)   = e_bio(i,j,p_bcary-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_mbo ) THEN
                         ebio_mbo(i,j) = ebio_mbo(i,j) + gas_emis
                         e_bio(i,j,p_mbo-1)   = e_bio(i,j,p_mbo-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ch3oh ) THEN
                         ebio_ch3oh(i,j) = ebio_ch3oh(i,j) + gas_emis
                         e_bio(i,j,p_ch3oh-1)   = e_bio(i,j,p_ch3oh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h5oh ) THEN
                         ebio_c2h5oh(i,j) = ebio_c2h5oh(i,j) + gas_emis
                         e_bio(i,j,p_c2h5oh-1)   = e_bio(i,j,p_c2h5oh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_hcho ) THEN
                         ebio_hcho(i,j) = ebio_hcho(i,j) + gas_emis
                         e_bio(i,j,p_hcho-1)   = e_bio(i,j,p_hcho-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ald ) THEN
                         ebio_ald(i,j) = ebio_ald(i,j) + gas_emis
                         e_bio(i,j,p_ald-1)   = e_bio(i,j,p_ald-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ch3cooh ) THEN
                         ebio_ch3cooh(i,j) = ebio_ch3cooh(i,j) + gas_emis
                         e_bio(i,j,p_ch3cooh-1)   = e_bio(i,j,p_ch3cooh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_hcooh ) THEN
                         ebio_hcooh(i,j) = ebio_hcooh(i,j) + gas_emis
                         e_bio(i,j,p_hcooh-1)   = e_bio(i,j,p_hcooh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_hcn ) THEN
                         ebio_hcn(i,j) = ebio_hcn(i,j) + gas_emis
                         e_bio(i,j,p_hcn-1)   = e_bio(i,j,p_hcn-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_nh3 ) THEN
                         ebio_nh3(i,j) = ebio_nh3(i,j) + gas_emis
                         e_bio(i,j,p_nh3-1)   = e_bio(i,j,p_nh3-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_co ) THEN
                         ebio_co(i,j)  = ebio_co(i,j) + gas_emis
                         e_bio(i,j,p_co-1)   = e_bio(i,j,p_co-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h4 ) THEN
                         ebio_c2h4(i,j) = ebio_c2h4(i,j) + gas_emis
                         e_bio(i,j,p_c2h4-1)   = e_bio(i,j,p_c2h4-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h6 ) THEN
                         ebio_c2h6(i,j) = ebio_c2h6(i,j) + gas_emis
                         e_bio(i,j,p_c2h6-1)   = e_bio(i,j,p_c2h6-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c3h6 ) THEN
                         ebio_c3h6(i,j) = ebio_c3h6(i,j) + gas_emis
                         e_bio(i,j,p_c3h6-1)   = e_bio(i,j,p_c3h6-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c3h8 ) THEN
                         ebio_c3h8(i,j) = ebio_c3h8(i,j) + gas_emis
                         e_bio(i,j,p_c3h8-1)   = e_bio(i,j,p_c3h8-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bigalk ) THEN
                         ebio_bigalk(i,j) = ebio_bigalk(i,j) + gas_emis
                         e_bio(i,j,p_bigalk-1)   = e_bio(i,j,p_bigalk-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bigene ) THEN
                         ebio_bigene(i,j) = ebio_bigene(i,j) + gas_emis
                         e_bio(i,j,p_bigene-1)   = e_bio(i,j,p_bigene-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_tol ) THEN
                         ebio_tol(i,j) = ebio_tol(i,j) + gas_emis
                         e_bio(i,j,p_tol-1)   = e_bio(i,j,p_tol-1)  + gas_emis*convert2
                      END IF
                   END IF is_t1_mozc_species
                END IF use_megan_emis_a
             END DO

          CASE ( MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP )

             DO icount = 1, n_megan2mozm



                p_in_chem = p_of_mozm(icount)
use_megan_emis : &
                IF ( p_in_chem /= non_react ) THEN



is_mozm_species : &
                   IF ( p_in_chem >= param_first_scalar ) THEN



                      gas_emis = mozm_per_megan(icount) * E_megan2(p_of_megan2mozm(icount))






                      IF ( p_in_chem == p_isopr ) THEN
                         ebio_iso(i,j) = ebio_iso(i,j) + gas_emis
                         e_bio(i,j,p_isopr-1)   = e_bio(i,j,p_isopr-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_no ) THEN
                         ebio_no(i,j)  = ebio_no(i,j) + gas_emis
                         e_bio(i,j,p_no-1)   = e_bio(i,j,p_no-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_no2 ) THEN
                         ebio_no2(i,j)  = ebio_no2(i,j) + gas_emis
                         e_bio(i,j,p_no2-1)   = e_bio(i,j,p_no2-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_co ) THEN
                         ebio_co(i,j)  = ebio_co(i,j) + gas_emis
                         e_bio(i,j,p_co-1)   = e_bio(i,j,p_co-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_hcho ) THEN
                         ebio_hcho(i,j) = ebio_hcho(i,j) + gas_emis
                         e_bio(i,j,p_hcho-1)   = e_bio(i,j,p_hcho-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ald ) THEN
                         ebio_ald(i,j) = ebio_ald(i,j) + gas_emis
                         e_bio(i,j,p_ald-1)   = e_bio(i,j,p_ald-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_acet ) THEN
                         ebio_acet(i,j) = ebio_acet(i,j) + gas_emis
                         e_bio(i,j,p_acet-1)   = e_bio(i,j,p_acet-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_tol ) THEN
                         ebio_tol(i,j) = ebio_tol(i,j) + gas_emis
                         e_bio(i,j,p_tol-1)   = e_bio(i,j,p_tol-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_apin ) THEN
                         ebio_api(i,j) = ebio_api(i,j) + gas_emis
                         e_bio(i,j,p_apin-1)   = e_bio(i,j,p_apin-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bpin ) THEN
                         ebio_bpi(i,j) = ebio_bpi(i,j) + gas_emis
                         e_bio(i,j,p_bpin-1)   = e_bio(i,j,p_bpin-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_limon ) THEN
                         ebio_lim(i,j) = ebio_lim(i,j) + gas_emis
                         e_bio(i,j,p_limon-1)   = e_bio(i,j,p_limon-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_mbo ) THEN
                         ebio_mbo(i,j) = ebio_mbo(i,j) + gas_emis
                         e_bio(i,j,p_mbo-1)   = e_bio(i,j,p_mbo-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_myrc ) THEN
                         ebio_myrc(i,j) = ebio_myrc(i,j) + gas_emis
                         e_bio(i,j,p_myrc-1)   = e_bio(i,j,p_myrc-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bcary ) THEN
                         ebio_sesq(i,j) = ebio_sesq(i,j) + gas_emis
                         e_bio(i,j,p_bcary-1)   = e_bio(i,j,p_bcary-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_so2 ) THEN
                         ebio_so2(i,j) = ebio_so2(i,j) + gas_emis
                         e_bio(i,j,p_so2-1)   = e_bio(i,j,p_so2-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_dms ) THEN
                         ebio_dms(i,j) = ebio_dms(i,j) + gas_emis
                         e_bio(i,j,p_dms-1)   = e_bio(i,j,p_dms-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bigalk ) THEN
                         ebio_bigalk(i,j) = ebio_bigalk(i,j) + gas_emis
                         e_bio(i,j,p_bigalk-1)   = e_bio(i,j,p_bigalk-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bigene ) THEN
                         ebio_bigene(i,j) = ebio_bigene(i,j) + gas_emis
                         e_bio(i,j,p_bigene-1)   = e_bio(i,j,p_bigene-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_nh3 ) THEN
                         ebio_nh3(i,j) = ebio_nh3(i,j) + gas_emis
                         e_bio(i,j,p_nh3-1)   = e_bio(i,j,p_nh3-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ch3oh ) THEN
                         ebio_ch3oh(i,j) = ebio_ch3oh(i,j) + gas_emis
                         e_bio(i,j,p_ch3oh-1)   = e_bio(i,j,p_ch3oh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h5oh ) THEN
                         ebio_c2h5oh(i,j) = ebio_c2h5oh(i,j) + gas_emis
                         e_bio(i,j,p_c2h5oh-1)   = e_bio(i,j,p_c2h5oh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ch3cooh ) THEN
                         ebio_ch3cooh(i,j) = ebio_ch3cooh(i,j) + gas_emis
                         e_bio(i,j,p_ch3cooh-1)   = e_bio(i,j,p_ch3cooh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_mek ) THEN
                         ebio_mek(i,j) = ebio_mek(i,j) + gas_emis
                         e_bio(i,j,p_mek-1)   = e_bio(i,j,p_mek-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h4 ) THEN
                         ebio_c2h4(i,j) = ebio_c2h4(i,j) + gas_emis
                         e_bio(i,j,p_c2h4-1)   = e_bio(i,j,p_c2h4-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h6 ) THEN
                         ebio_c2h6(i,j) = ebio_c2h6(i,j) + gas_emis
                         e_bio(i,j,p_c2h6-1)   = e_bio(i,j,p_c2h6-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c3h6 ) THEN
                         ebio_c3h6(i,j) = ebio_c3h6(i,j) + gas_emis
                         e_bio(i,j,p_c3h6-1)   = e_bio(i,j,p_c3h6-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c3h8 ) THEN
                         ebio_c3h8(i,j) = ebio_c3h8(i,j) + gas_emis
                         e_bio(i,j,p_c3h8-1)   = e_bio(i,j,p_c3h8-1)  + gas_emis*convert2
                      END IF
                   END IF is_mozm_species
                END IF use_megan_emis
             END DO

          CASE (RADM2, RADM2_KPP, RADM2SORG, RADM2SORG_AQ, RADM2SORG_AQCHEM, RADM2SORG_KPP,GOCARTRADM2)

             DO icount = 1, n_megan2radm2

                IF ( p_of_radm2(icount) .NE. non_react ) THEN
                
                   
                   
                   p_in_chem = p_of_radm2(icount)

                   
                   IF ( p_in_chem >= param_first_scalar ) THEN
                      
                      
                      gas_emis = radm2_per_megan(icount) * E_megan2(p_of_megan2radm2(icount))

                      
                      
                      
                      
                      IF ( p_in_chem .EQ. p_iso ) THEN
                         ebio_iso(i,j)        = ebio_iso(i,j)       + gas_emis
                         e_bio(i,j,p_iso-1)   = e_bio(i,j,p_iso-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_oli) THEN
                         ebio_oli(i,j)        = ebio_oli(i,j)       + gas_emis
                         e_bio(i,j,p_oli-1)   = e_bio(i,j,p_oli-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc3) THEN
                         ebio_hc3(i,j)        = ebio_hc3(i,j)       + gas_emis
                         e_bio(i,j,p_hc3-1)   = e_bio(i,j,p_hc3-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_olt) THEN
                         ebio_olt(i,j)        = ebio_olt(i,j)       + gas_emis
                         e_bio(i,j,p_olt-1)   = e_bio(i,j,p_olt-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ket) THEN
                         ebio_ket(i,j)        = ebio_ket(i,j)       + gas_emis
                         e_bio(i,j,p_ket-1)   = e_bio(i,j,p_ket-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ald) THEN
                         ebio_ald(i,j)        = ebio_ald(i,j)       + gas_emis
                         e_bio(i,j,p_ald-1)   = e_bio(i,j,p_ald-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hcho) THEN
                         ebio_hcho(i,j)       = ebio_hcho(i,j)      + gas_emis
                         e_bio(i,j,p_hcho-1)  = e_bio(i,j,p_hcho-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_eth) THEN
                         ebio_eth(i,j)        = ebio_eth(i,j)       + gas_emis
                         e_bio(i,j,p_eth-1)   = e_bio(i,j,p_eth-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ora2) THEN
                         ebio_ora2(i,j)       = ebio_ora2(i,j)      + gas_emis
                         e_bio(i,j,p_ora2-1)  = e_bio(i,j,p_ora2-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_co) THEN
                         ebio_co(i,j)         = ebio_co(i,j)        + gas_emis
                         e_bio(i,j,p_co-1)    = e_bio(i,j,p_co-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_no) THEN
                         ebio_no(i,j)         = ebio_no(i,j)        + gas_emis   
                         e_bio(i,j,p_no-1)    = e_bio(i,j,p_no-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ol2) THEN
                          e_bio(i,j,p_ol2-1)  = e_bio(i,j,p_ol2-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc5) THEN
                          e_bio(i,j,p_hc5-1)  = e_bio(i,j,p_hc5-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc8) THEN
                          e_bio(i,j,p_hc8-1)  = e_bio(i,j,p_hc8-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ora1) THEN
                          e_bio(i,j,p_ora1-1) = e_bio(i,j,p_ora1-1) + gas_emis*convert2
                      END IF

                   END IF 

                END IF 
                
             END DO

           CASE (RACMSORG_AQ, RACMSORG_AQCHEM_KPP, RACM_ESRLSORG_AQCHEM_KPP, RACM_ESRLSORG_KPP, RACM_KPP, GOCARTRACM_KPP, &
                 RACMSORG_KPP, RACM_MIM_KPP, RACMPM_KPP)

             DO icount = 1, n_megan2racm

                IF ( p_of_racm(icount) .NE. non_react ) THEN

                   
                   
                   p_in_chem = p_of_racm(icount)
                   
                   
                   IF( p_in_chem >= param_first_scalar ) THEN

                      
                      gas_emis =  racm_per_megan(icount) * E_megan2(p_of_megan2racm(icount))

                      
                      
                      
                      
                      IF ( p_in_chem .EQ. p_iso ) THEN
                         ebio_iso(i,j)        = ebio_iso(i,j)       + gas_emis
                         e_bio(i,j,p_iso-1)   = e_bio(i,j,p_iso-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_oli) THEN
                         ebio_oli(i,j)        = ebio_oli(i,j)       + gas_emis
                         e_bio(i,j,p_oli-1)   = e_bio(i,j,p_oli-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_api) THEN
                         ebio_api(i,j)        = ebio_api(i,j)       + gas_emis
                         e_bio(i,j,p_api-1)   = e_bio(i,j,p_api-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_lim) THEN
                         ebio_lim(i,j)        = ebio_lim(i,j)       + gas_emis
                         e_bio(i,j,p_lim-1)   = e_bio(i,j,p_lim-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc3) THEN
                         ebio_hc3(i,j)        = ebio_hc3(i,j)       + gas_emis
                         e_bio(i,j,p_hc3-1)   = e_bio(i,j,p_hc3-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ete) THEN
                         ebio_ete(i,j)        = ebio_ete(i,j)       + gas_emis
                         e_bio(i,j,p_ete-1)   = e_bio(i,j,p_ete-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_olt) THEN
                         ebio_olt(i,j)        = ebio_olt(i,j)       + gas_emis
                         e_bio(i,j,p_olt-1)   = e_bio(i,j,p_olt-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ket) THEN
                         ebio_ket(i,j)        = ebio_ket(i,j)       + gas_emis
                         e_bio(i,j,p_ket-1)   = e_bio(i,j,p_ket-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ald) THEN
                         ebio_ald(i,j)        = ebio_ald(i,j)       + gas_emis
                         e_bio(i,j,p_ald-1)   = e_bio(i,j,p_ald-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hcho) THEN
                         ebio_hcho(i,j)       = ebio_hcho(i,j)      + gas_emis
                         e_bio(i,j,p_hcho-1)  = e_bio(i,j,p_hcho-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_eth) THEN
                         ebio_eth(i,j)        = ebio_eth(i,j)       + gas_emis
                         e_bio(i,j,p_eth-1)   = e_bio(i,j,p_eth-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ora2) THEN
                         ebio_ora2(i,j)       = ebio_ora2(i,j)      + gas_emis
                         e_bio(i,j,p_ora2-1)  = e_bio(i,j,p_ora2-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_co) THEN
                         ebio_co(i,j)         = ebio_co(i,j)        + gas_emis
                         e_bio(i,j,p_co-1)    = e_bio(i,j,p_co-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_no) THEN
                         ebio_no(i,j)         = ebio_no(i,j)        + gas_emis   
                         e_bio(i,j,p_no-1)    = e_bio(i,j,p_no-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc5) THEN
                          e_bio(i,j,p_hc5-1)  = e_bio(i,j,p_hc5-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc8) THEN
                          e_bio(i,j,p_hc8-1)  = e_bio(i,j,p_hc8-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ora1) THEN
                          e_bio(i,j,p_ora1-1) = e_bio(i,j,p_ora1-1) + gas_emis*convert2
                      END IF

                   END IF 
                   

                END IF 

             END DO

          CASE (CB05_SORG_AQ_KPP)

             DO icount = 1, n_megan2cb05
                IF ( p_of_cb05 (icount) .NE. non_react ) THEN
                   
                   
                   p_in_chem = p_of_cb05(icount)

                   
                   
                   
                   
                   IF ( p_in_chem >= param_first_scalar ) THEN

                      
                      gas_emis = cb05_per_megan(icount) * E_megan2(p_of_megan2cb05(icount))

                      
                      
                      
                      
                      IF ( p_in_chem .EQ. p_isop ) THEN
                         ebio_iso(i,j)        = ebio_iso(i,j)       + gas_emis
                         e_bio(i,j,p_isop-1)   = e_bio(i,j,p_isop-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_aacd ) THEN
                         e_bio(i,j,p_aacd-1)  = e_bio(i,j,p_aacd-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_ald2 ) THEN
                         ebio_ald(i,j)        = ebio_ald(i,j)       + gas_emis
                         e_bio(i,j,p_ald2-1)  = e_bio(i,j,p_ald2-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_aldx ) THEN
                         ebio_ald(i,j)        = ebio_ald(i,j)       + gas_emis
                         e_bio(i,j,p_aldx-1)  = e_bio(i,j,p_aldx-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_apin ) THEN
                         ebio_api(i,j)        = ebio_api(i,j)       + gas_emis
                         e_bio(i,j,p_terp-1)  = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_bpin ) THEN
                         e_bio(i,j,p_terp-1)  = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_ch4 ) THEN
                         e_bio(i,j,p_ch4-1)   = e_bio(i,j,p_ch4-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_co ) THEN
                         ebio_co(i,j)        = ebio_co(i,j)       + gas_emis
                         e_bio(i,j,p_co-1)    = e_bio(i,j,p_co-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_eth ) THEN
                         e_bio(i,j,p_eth-1)   = e_bio(i,j,p_eth-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_etha ) THEN
                         e_bio(i,j,p_etha-1)  = e_bio(i,j,p_etha-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_etoh ) THEN
                         e_bio(i,j,p_etoh-1)  = e_bio(i,j,p_etoh-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_facd ) THEN
                         e_bio(i,j,p_facd-1)  = e_bio(i,j,p_facd-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_form ) THEN
                         ebio_hcho(i,j)        = ebio_hcho(i,j)       + gas_emis
                         e_bio(i,j,p_form-1)  = e_bio(i,j,p_form-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_hum ) THEN
                         e_bio(i,j,p_terp-1)  = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_iole ) THEN
                         e_bio(i,j,p_iole-1)  = e_bio(i,j,p_iole-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_lim ) THEN
                         ebio_lim(i,j)        = ebio_lim(i,j)       + gas_emis
                         e_bio(i,j,p_terp-1)  = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_meoh ) THEN
                         e_bio(i,j,p_meoh-1)  = e_bio(i,j,p_meoh-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_nh3 ) THEN
                         e_bio(i,j,p_nh3-1)   = e_bio(i,j,p_nh3-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_no ) THEN
                         ebio_no(i,j)        = ebio_no(i,j)       + gas_emis
                         e_bio(i,j,p_no-1)    = e_bio(i,j,p_no-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_oci ) THEN
                         e_bio(i,j,p_terp-1)  = e_bio(i,j,p_terp-1) + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_ole ) THEN
                         e_bio(i,j,p_ole-1)   = e_bio(i,j,p_ole-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_par ) THEN
                         e_bio(i,j,p_par-1)   = e_bio(i,j,p_par-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_terp ) THEN
                         ebio_terp(i,j)        = ebio_terp(i,j)       + gas_emis
                         e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_ter ) THEN
                         ebio_terp(i,j)        = ebio_terp(i,j)       + gas_emis
                         e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_tol ) THEN
                         e_bio(i,j,p_tol-1)   = e_bio(i,j,p_tol-1)  + gas_emis*convert2
                      END IF

                   END IF 

                END IF
             END DO

          CASE (CB05_SORG_VBS_AQ_KPP)

             DO icount = 1, n_megan2cb05vbs
                IF ( p_of_cb05vbs (icount) .NE. non_react ) THEN
                   
                   
                   p_in_chem = p_of_cb05vbs(icount)

                   
                   
                   
                   
                   IF ( p_in_chem >= param_first_scalar ) THEN

                      
                      gas_emis = cb05vbs_per_megan(icount) * E_megan2(p_of_megan2cb05vbs(icount))

                      
                      
                      
                      
                      IF ( p_in_chem .EQ. p_isop ) THEN
                         ebio_iso(i,j)        = ebio_iso(i,j)       + gas_emis
                         e_bio(i,j,p_isop-1)   = e_bio(i,j,p_isop-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_aacd ) THEN
                         e_bio(i,j,p_aacd-1)  = e_bio(i,j,p_aacd-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_ald2 ) THEN
                         ebio_ald(i,j)        = ebio_ald(i,j)       + gas_emis
                         e_bio(i,j,p_ald2-1)  = e_bio(i,j,p_ald2-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_aldx ) THEN
                         ebio_ald(i,j)        = ebio_ald(i,j)       + gas_emis
                         e_bio(i,j,p_aldx-1)  = e_bio(i,j,p_aldx-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_apin ) THEN
                         ebio_api(i,j)        = ebio_api(i,j)       + gas_emis
                         e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_bpin ) THEN
                         e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_ch4 ) THEN
                         e_bio(i,j,p_ch4-1)   = e_bio(i,j,p_ch4-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_co ) THEN
                         ebio_co(i,j)        = ebio_co(i,j)       + gas_emis
                         e_bio(i,j,p_co-1)    = e_bio(i,j,p_co-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_eth ) THEN
                         e_bio(i,j,p_eth-1)   = e_bio(i,j,p_eth-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_etha ) THEN
                         e_bio(i,j,p_etha-1)  = e_bio(i,j,p_etha-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_etoh ) THEN
                         e_bio(i,j,p_etoh-1)  = e_bio(i,j,p_etoh-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_facd ) THEN
                         e_bio(i,j,p_facd-1)  = e_bio(i,j,p_facd-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_form ) THEN
                         ebio_hcho(i,j)        = ebio_hcho(i,j)       + gas_emis
                         e_bio(i,j,p_form-1)  = e_bio(i,j,p_form-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_hum ) THEN
                          e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_iole ) THEN
                         e_bio(i,j,p_iole-1)  = e_bio(i,j,p_iole-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_lim ) THEN
                         ebio_lim(i,j)        = ebio_lim(i,j)       + gas_emis
                         e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_meoh ) THEN
                         e_bio(i,j,p_meoh-1)  = e_bio(i,j,p_meoh-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_nh3 ) THEN
                         e_bio(i,j,p_nh3-1)   = e_bio(i,j,p_nh3-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_no ) THEN
                         ebio_no(i,j)        = ebio_no(i,j)       + gas_emis
                         e_bio(i,j,p_no-1)    = e_bio(i,j,p_no-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_oci ) THEN
                         e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_ole ) THEN
                         e_bio(i,j,p_ole-1)   = e_bio(i,j,p_ole-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_par ) THEN
                         e_bio(i,j,p_par-1)   = e_bio(i,j,p_par-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_terp ) THEN
                         ebio_terp(i,j)        = ebio_terp(i,j)       + gas_emis
                         e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_ter ) THEN
                         ebio_terp(i,j)        = ebio_terp(i,j)       + gas_emis
                         e_bio(i,j,p_terp-1)   = e_bio(i,j,p_terp-1)  + gas_emis*convert2
                      END IF
                      IF ( p_in_chem .EQ. p_tol ) THEN
                         e_bio(i,j,p_tol-1)   = e_bio(i,j,p_tol-1)  + gas_emis*convert2
                      END IF

                   END IF 

                END IF
             END DO

          CASE (RACM_SOA_VBS_KPP,RACM_SOA_VBS_AQCHEM_KPP,RACM_SOA_VBS_HET_KPP)

          DO icount = 1, n_megan2racmSOA

                IF ( p_of_racmSOA(icount) .NE. non_react ) THEN

                   
                   
                   p_in_chem = p_of_racmSOA(icount)

                   
                   IF( p_in_chem >= param_first_scalar ) THEN

                      
                      gas_emis =  racmSOA_per_megan(icount) * E_megan2(p_of_megan2racmSOA(icount))

                      
                      
                      
                      
                      IF ( p_in_chem .EQ. p_iso ) THEN
                         ebio_iso(i,j)        = ebio_iso(i,j)       + gas_emis
                         e_bio(i,j,p_iso-1)   = e_bio(i,j,p_iso-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_oli) THEN
                         ebio_oli(i,j)        = ebio_oli(i,j)       + gas_emis
                         e_bio(i,j,p_oli-1)   = e_bio(i,j,p_oli-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_api) THEN
                         ebio_api(i,j)        = ebio_api(i,j)       + gas_emis
                         e_bio(i,j,p_api-1)   = e_bio(i,j,p_api-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_lim) THEN
                         ebio_lim(i,j)        = ebio_lim(i,j)       + gas_emis
                         e_bio(i,j,p_lim-1)   = e_bio(i,j,p_lim-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc3) THEN
                         ebio_hc3(i,j)        = ebio_hc3(i,j)       + gas_emis
                         e_bio(i,j,p_hc3-1)   = e_bio(i,j,p_hc3-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ete) THEN
                         ebio_ete(i,j)        = ebio_ete(i,j)       + gas_emis
                         e_bio(i,j,p_ete-1)   = e_bio(i,j,p_ete-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_olt) THEN
                         ebio_olt(i,j)        = ebio_olt(i,j)       + gas_emis
                         e_bio(i,j,p_olt-1)   = e_bio(i,j,p_olt-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ket) THEN
                         ebio_ket(i,j)        = ebio_ket(i,j)       + gas_emis
                         e_bio(i,j,p_ket-1)   = e_bio(i,j,p_ket-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ald) THEN
                         ebio_ald(i,j)        = ebio_ald(i,j)       + gas_emis
                         e_bio(i,j,p_ald-1)   = e_bio(i,j,p_ald-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hcho) THEN
                         ebio_hcho(i,j)       = ebio_hcho(i,j)      + gas_emis
                         e_bio(i,j,p_hcho-1)  = e_bio(i,j,p_hcho-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_eth) THEN
                         ebio_eth(i,j)        = ebio_eth(i,j)       + gas_emis
                         e_bio(i,j,p_eth-1)   = e_bio(i,j,p_eth-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ora2) THEN
                         ebio_ora2(i,j)       = ebio_ora2(i,j)      + gas_emis
                         e_bio(i,j,p_ora2-1)  = e_bio(i,j,p_ora2-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_co) THEN
                         ebio_co(i,j)         = ebio_co(i,j)        + gas_emis
                         e_bio(i,j,p_co-1)    = e_bio(i,j,p_co-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_no) THEN
                         ebio_no(i,j)         = ebio_no(i,j)        + gas_emis
                         e_bio(i,j,p_no-1)    = e_bio(i,j,p_no-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc5) THEN
                          e_bio(i,j,p_hc5-1)  = e_bio(i,j,p_hc5-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hc8) THEN
                          e_bio(i,j,p_hc8-1)  = e_bio(i,j,p_hc8-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ora1) THEN
                          e_bio(i,j,p_ora1-1) = e_bio(i,j,p_ora1-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_sesq) THEN
                          ebio_sesq(i,j)      = ebio_sesq(i,j)      + gas_emis
                          e_bio(i,j,p_sesq-1) = e_bio(i,j,p_sesq-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_mbo) THEN
                          ebio_mbo(i,j)       = ebio_mbo(i,j)        + gas_emis
                          e_bio(i,j,p_mbo-1)  = e_bio(i,j,p_mbo-1)   + gas_emis*convert2
                      END IF

                   END IF 


                END IF 

             END DO
          CASE (CBMZ, CBMZ_BB, CBMZ_BB_KPP, CBMZ_MOSAIC_KPP, &
                CBMZ_MOSAIC_4BIN, &
                CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
                CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, &
                CBMZ_MOSAIC_DMS_4BIN_AQ,CBMZ_MOSAIC_DMS_8BIN_AQ,CBMZSORG, CBMZSORG_AQ, &
                CBMZ_CAM_MAM3_NOAQ, CBMZ_CAM_MAM3_AQ, CBMZ_CAM_MAM7_NOAQ, CBMZ_CAM_MAM7_AQ)

             DO icount = 1, n_megan2cbmz

                IF ( p_of_cbmz (icount) .NE. non_react ) THEN

                   
                   
                   p_in_chem = p_of_cbmz(icount)

                   
                   
                   
                   IF( p_in_chem >= param_first_scalar ) THEN

                      
                      gas_emis = cbmz_per_megan(icount) * E_megan2(p_of_megan2cbmz(icount))


                      
                      
                      
                      
                      IF ( p_in_chem .EQ. p_iso ) THEN
                         ebio_iso(i,j)        = ebio_iso(i,j)       + gas_emis
                         e_bio(i,j,p_iso-1)   = e_bio(i,j,p_iso-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_oli) THEN
                         ebio_oli(i,j)        = ebio_oli(i,j)       + gas_emis
                         e_bio(i,j,p_oli-1)   = e_bio(i,j,p_oli-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_olt) THEN
                         ebio_olt(i,j)        = ebio_olt(i,j)       + gas_emis
                         e_bio(i,j,p_olt-1)   = e_bio(i,j,p_olt-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ket) THEN
                         ebio_ket(i,j)        = ebio_ket(i,j)       + gas_emis
                         e_bio(i,j,p_ket-1)   = e_bio(i,j,p_ket-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ald) THEN
                         ebio_ald(i,j)        = ebio_ald(i,j)       + gas_emis
                         e_bio(i,j,p_ald-1)   = e_bio(i,j,p_ald-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hcho) THEN
                         ebio_hcho(i,j)       = ebio_hcho(i,j)      + gas_emis
                         e_bio(i,j,p_hcho-1)  = e_bio(i,j,p_hcho-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_eth) THEN
                         ebio_eth(i,j)        = ebio_eth(i,j)       + gas_emis
                         e_bio(i,j,p_eth-1)   = e_bio(i,j,p_eth-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ora2) THEN
                         ebio_ora2(i,j)       = ebio_ora2(i,j)      + gas_emis
                         e_bio(i,j,p_ora2-1)  = e_bio(i,j,p_ora2-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_co) THEN
                         ebio_co(i,j)         = ebio_co(i,j)        + gas_emis
                         e_bio(i,j,p_co-1)    = e_bio(i,j,p_co-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_no) THEN
                         ebio_no(i,j)         = ebio_no(i,j)        + gas_emis   
                         e_bio(i,j,p_no-1)    = e_bio(i,j,p_no-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ol2) THEN
                          e_bio(i,j,p_ol2-1)  = e_bio(i,j,p_ol2-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ora1) THEN
                          e_bio(i,j,p_ora1-1) = e_bio(i,j,p_ora1-1) + gas_emis*convert2
                      
                      
                      
                      ELSE IF ( p_in_chem .EQ. p_par) THEN 
                         
                         e_bio(i,j,p_par-1)   = e_bio(i,j,p_par-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ch3oh) THEN	
                         ebio_ch3oh(i,j)      = ebio_ch3oh(i,j)     + gas_emis
                         e_bio(i,j,p_ch3oh-1) = e_bio(i,j,p_ch3oh-1)+ gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_c2h5oh) THEN	
                         ebio_c2h5oh(i,j)     = ebio_c2h5oh(i,j)      + gas_emis
                         e_bio(i,j,p_c2h5oh-1)= e_bio(i,j,p_c2h5oh-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_nh3) THEN	
                         ebio_nh3(i,j)        = ebio_nh3(i,j)       + gas_emis
                         e_bio(i,j,p_nh3-1)   = e_bio(i,j,p_nh3-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_tol) THEN	
                         ebio_tol(i,j)       = ebio_tol(i,j)        + gas_emis
                         e_bio(i,j,p_tol-1)  = e_bio(i,j,p_tol-1)   + gas_emis*convert2

                      END IF


                   END IF 
                   
                   
                END IF 

             END DO
            
          CASE (SAPRC99_KPP,SAPRC99_MOSAIC_4BIN_VBS2_KPP, &
               SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_8BIN_VBS2_KPP)

             DO icount = 1, n_megan2saprcnov

                IF ( p_of_saprcnov(icount) .NE. non_react ) THEN

                   
                   
                   p_in_chem = p_of_saprcnov(icount)

                   
                   IF ( p_in_chem >= param_first_scalar ) THEN

                      
                      gas_emis = saprcnov_per_megan(icount) * E_megan2(p_of_megan2saprcnov(icount))

                      
                      
                      
                      
                      IF ( p_in_chem .EQ. p_isoprene ) THEN
                         ebio_iso(i,j)        = ebio_iso(i,j)       + gas_emis
                         e_bio(i,j,p_isoprene-1)   = e_bio(i,j,p_isoprene-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_terp) THEN
                         ebio_api(i,j)       = ebio_api(i,j)      + gas_emis
                         e_bio(i,j,p_terp-1)  = e_bio(i,j,p_terp-1) + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_sesq) THEN
                         ebio_lim(i,j)         = ebio_lim(i,j)        + gas_emis
                         e_bio(i,j,p_sesq-1)    = e_bio(i,j,p_sesq-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_no) THEN
                         ebio_no(i,j)         = ebio_no(i,j)        + gas_emis
                         e_bio(i,j,p_no-1)    = e_bio(i,j,p_no-1)   + gas_emis*convert2

                      ELSE IF ( p_in_chem .EQ. p_alk3) THEN
                         ebio_alk3(i,j)         = ebio_alk3(i,j)        + gas_emis
                         e_bio(i,j,p_alk3-1)    = e_bio(i,j,p_alk3-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_alk4) THEN
                         ebio_alk4(i,j)         = ebio_alk4(i,j)        + gas_emis
                         e_bio(i,j,p_alk4-1)    = e_bio(i,j,p_alk4-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_alk5) THEN
                         ebio_alk5(i,j)         = ebio_alk5(i,j)        + gas_emis
                         e_bio(i,j,p_alk5-1)    = e_bio(i,j,p_alk5-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ole1) THEN
                         ebio_ole1(i,j)         = ebio_ole1(i,j)        + gas_emis
                         e_bio(i,j,p_ole1-1)    = e_bio(i,j,p_ole1-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ole2) THEN
                         ebio_ole2(i,j)         = ebio_ole2(i,j)        + gas_emis
                         e_bio(i,j,p_ole2-1)    = e_bio(i,j,p_ole2-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_aro1) THEN
                         ebio_aro1(i,j)         = ebio_aro1(i,j)        + gas_emis
                         e_bio(i,j,p_aro1-1)    = e_bio(i,j,p_aro1-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_aro2) THEN
                         ebio_aro2(i,j)         = ebio_aro2(i,j)        + gas_emis
                         e_bio(i,j,p_aro2-1)    = e_bio(i,j,p_aro2-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_acet) THEN
                         ebio_acet(i,j)         = ebio_acet(i,j)        + gas_emis
                         e_bio(i,j,p_acet-1)    = e_bio(i,j,p_acet-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hcho) THEN
                         ebio_hcho(i,j)         = ebio_hcho(i,j)        + gas_emis
                         e_bio(i,j,p_hcho-1)    = e_bio(i,j,p_hcho-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ccho) THEN
                         ebio_ccho(i,j)         = ebio_ccho(i,j)        + gas_emis
                         e_bio(i,j,p_ccho-1)    = e_bio(i,j,p_ccho-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_mek) THEN
                         ebio_mek(i,j)         = ebio_mek(i,j)        + gas_emis
                         e_bio(i,j,p_mek-1)    = e_bio(i,j,p_mek-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_c2h6) THEN
                         ebio_c2h6(i,j)         = ebio_c2h6(i,j)        + gas_emis
                         e_bio(i,j,p_c2h6-1)    = e_bio(i,j,p_c2h6-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_c3h6) THEN
                         ebio_c3h6(i,j)         = ebio_c3h6(i,j)        + gas_emis
                         e_bio(i,j,p_c3h6-1)    = e_bio(i,j,p_c3h6-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_c3h8) THEN
                         ebio_c3h8(i,j)         = ebio_c3h8(i,j)        + gas_emis
                         e_bio(i,j,p_c3h8-1)    = e_bio(i,j,p_c3h8-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_ethene) THEN
                         ebio_ethene(i,j)         = ebio_ethene(i,j)        + gas_emis
                         e_bio(i,j,p_ethene-1)    = e_bio(i,j,p_ethene-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_bald) THEN
                         ebio_bald(i,j)         = ebio_bald(i,j)        + gas_emis
                         e_bio(i,j,p_bald-1)    = e_bio(i,j,p_bald-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_meoh) THEN
                         ebio_meoh(i,j)         = ebio_meoh(i,j)        + gas_emis
                         e_bio(i,j,p_meoh-1)    = e_bio(i,j,p_meoh-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_hcooh) THEN
                         ebio_hcooh(i,j)         = ebio_hcooh(i,j)        + gas_emis
                         e_bio(i,j,p_hcooh-1)    = e_bio(i,j,p_hcooh-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_rco_oh) THEN
                         ebio_rco_oh(i,j)         = ebio_rco_oh(i,j)        + gas_emis
                         e_bio(i,j,p_rco_oh-1)    = e_bio(i,j,p_rco_oh-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_terp) THEN
                         ebio_terp(i,j)         = ebio_terp(i,j)        + gas_emis
                         ebio_api(i,j)         = ebio_api(i,j)        + gas_emis
                         e_bio(i,j,p_terp-1)    = e_bio(i,j,p_terp-1)   + gas_emis*convert2
                      ELSE IF ( p_in_chem .EQ. p_sesq) THEN
                         ebio_sesq(i,j)         = ebio_sesq(i,j)        + gas_emis
                         ebio_lim(i,j)         = ebio_lim(i,j)        + gas_emis
                         e_bio(i,j,p_sesq-1)    = e_bio(i,j,p_sesq-1)   + gas_emis*convert2

                      END IF

                   END IF 

                END IF 

             END DO

          CASE ( CRIMECH_KPP, CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP )

             DO icount = 1, n_megan2crimech
                IF ( p_of_crimech(icount) .NE. non_react ) THEN

                   
                   
                   p_in_chem = p_of_crimech(icount)
                   
                   
                   IF( p_in_chem >= param_first_scalar ) THEN

                      
                      gas_emis =  crimech_per_megan(icount) * E_megan2(p_of_megan2crimech(icount))

                      
                      
                      
                      
                      
                      IF ( p_in_chem == p_c5h8 ) THEN
                         ebio_c5h8(i,j) = ebio_c5h8(i,j) + gas_emis
                         e_bio(i,j,p_c5h8-1)   = e_bio(i,j,p_c5h8-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_no ) THEN
                         ebio_no(i,j)  = ebio_no(i,j) + gas_emis
                         e_bio(i,j,p_no-1)   = e_bio(i,j,p_no-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_no2 ) THEN
                         ebio_no2(i,j)  = ebio_no2(i,j) + gas_emis
                         e_bio(i,j,p_no2-1)   = e_bio(i,j,p_no2-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_co ) THEN
                         ebio_co(i,j)  = ebio_co(i,j) + gas_emis
                         e_bio(i,j,p_co-1)   = e_bio(i,j,p_co-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_hcho ) THEN
                         ebio_hcho(i,j) = ebio_hcho(i,j) + gas_emis
                         e_bio(i,j,p_hcho-1)   = e_bio(i,j,p_hcho-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ket ) THEN
                         ebio_ket(i,j) = ebio_ket(i,j) + gas_emis
                         e_bio(i,j,p_ket-1)   = e_bio(i,j,p_ket-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_toluene ) THEN
                         ebio_toluene(i,j) = ebio_toluene(i,j) + gas_emis
                         e_bio(i,j,p_toluene-1)   = e_bio(i,j,p_toluene-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_apinene ) THEN
                         ebio_apinene(i,j) = ebio_apinene(i,j) + gas_emis
                         e_bio(i,j,p_apinene-1)   = e_bio(i,j,p_apinene-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_bpinene ) THEN
                         ebio_bpinene(i,j) = ebio_bpinene(i,j) + gas_emis
                         e_bio(i,j,p_bpinene-1)   = e_bio(i,j,p_bpinene-1)  + gas_emis*convert2                         
                      ELSE IF ( p_in_chem == p_so2 ) THEN
                         ebio_so2(i,j) = ebio_so2(i,j) + gas_emis
                         e_bio(i,j,p_so2-1)   = e_bio(i,j,p_so2-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_dms ) THEN
                         ebio_dms(i,j) = ebio_dms(i,j) + gas_emis
                         e_bio(i,j,p_dms-1)   = e_bio(i,j,p_dms-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_nc4h10 ) THEN
                         ebio_nc4h10(i,j) = ebio_nc4h10(i,j) + gas_emis
                         e_bio(i,j,p_nc4h10-1)   = e_bio(i,j,p_nc4h10-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_tbut2ene ) THEN
                         ebio_tbut2ene(i,j) = ebio_tbut2ene(i,j) + gas_emis
                         e_bio(i,j,p_tbut2ene-1)   = e_bio(i,j,p_tbut2ene-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_nh3 ) THEN
                         ebio_nh3(i,j) = ebio_nh3(i,j) + gas_emis
                         e_bio(i,j,p_nh3-1)   = e_bio(i,j,p_nh3-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ch3oh ) THEN
                         ebio_ch3oh(i,j) = ebio_ch3oh(i,j) + gas_emis
                         e_bio(i,j,p_ch3oh-1)   = e_bio(i,j,p_ch3oh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h5oh ) THEN
                         ebio_c2h5oh(i,j) = ebio_c2h5oh(i,j) + gas_emis
                         e_bio(i,j,p_c2h5oh-1)   = e_bio(i,j,p_c2h5oh-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_ch3co2h ) THEN
                         ebio_ch3co2h(i,j) = ebio_ch3co2h(i,j) + gas_emis
                         e_bio(i,j,p_ch3co2h-1)   = e_bio(i,j,p_ch3co2h-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_mek ) THEN
                         ebio_mek(i,j) = ebio_mek(i,j) + gas_emis
                         e_bio(i,j,p_mek-1)   = e_bio(i,j,p_mek-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h4 ) THEN
                         ebio_c2h4(i,j) = ebio_c2h4(i,j) + gas_emis
                         e_bio(i,j,p_c2h4-1)   = e_bio(i,j,p_c2h4-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c2h6 ) THEN
                         ebio_c2h6(i,j) = ebio_c2h6(i,j) + gas_emis
                         e_bio(i,j,p_c2h6-1)   = e_bio(i,j,p_c2h6-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c3h6 ) THEN
                         ebio_c3h6(i,j) = ebio_c3h6(i,j) + gas_emis
                         e_bio(i,j,p_c3h6-1)   = e_bio(i,j,p_c3h6-1)  + gas_emis*convert2
                      ELSE IF ( p_in_chem == p_c3h8 ) THEN
                         ebio_c3h8(i,j) = ebio_c3h8(i,j) + gas_emis
                         e_bio(i,j,p_c3h8-1)   = e_bio(i,j,p_c3h8-1)  + gas_emis*convert2                         
                      ELSE IF ( p_in_chem == p_ch3cho ) THEN
                         ebio_ch3cho(i,j) = ebio_ch3cho(i,j) + gas_emis
                         e_bio(i,j,p_ch3cho-1)   = e_bio(i,j,p_ch3cho-1)  + gas_emis*convert2                        
                      ELSE IF ( p_in_chem == p_hcooh ) THEN
                         ebio_hcooh(i,j) = ebio_hcooh(i,j) + gas_emis
                         e_bio(i,j,p_hcooh-1)   = e_bio(i,j,p_hcooh-1)  + gas_emis*convert2                         
                      END IF

                   END IF 
                   

                END IF 

             END DO


 
             CASE DEFAULT

                CALL wrf_error_fatal3("<stdin>",1840,&
'Species conversion table for MEGAN v2.04 not available. ')

             END SELECT GAS_MECH_SELECT



       END DO i_loop 
    END DO j_loop    


  CONTAINS

    
    

    

    SUBROUTINE GAMMA_TISOP( TEMP, D_TEMP, gam_T )
      
      
      
      
      
      
      
      
      
      
      
      
      
      

      IMPLICIT NONE

      
      
      REAL, INTENT(IN)  :: TEMP
      

      REAL, INTENT(IN)  :: D_TEMP
      
      REAL, INTENT(OUT) :: gam_T

      
      REAL :: Eopt, Topt, X
      REAL :: AAA, BBB
      REAL, PARAMETER :: CT1 = 80.0
      REAL, PARAMETER :: CT2 = 200.0
      
      

      
      
      Eopt = 1.75 * EXP(0.08*(D_TEMP-297.0))

      
      
      Topt = 313.0 + ( 0.6*(D_TEMP-297.0) )

      
      X = ( (1.0/Topt)-(1.0/TEMP) ) / 0.00831
      AAA = Eopt*CT2*EXP(CT1*X)
      BBB = (  CT2-CT1*( 1.0-EXP(CT2*X) )  )
      gam_T = AAA/BBB

    END SUBROUTINE GAMMA_TISOP

    
    
    

    

    SUBROUTINE GAMMA_TNISP( SPCNUM, TEMP, gam_T )
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: SPCNUM               
      REAL, INTENT(IN)    :: TEMP
      REAL, INTENT(OUT)   :: gam_T
      REAL, PARAMETER     :: Ts = 303.0
      
      

      
      gam_T = EXP( TDF_PRM(SPCNUM)*(TEMP-Ts) )

    END SUBROUTINE GAMMA_TNISP


    
    

    

    SUBROUTINE GAMMA_LAI(LAI, gam_L )
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

      IMPLICIT NONE
      REAL, INTENT(IN)  ::  LAI 
      REAL, INTENT(OUT) :: gam_L

      

      
      
      gam_L = (0.49*LAI) / ( SQRT(1.0+0.2*(LAI**2)) )

      RETURN
    END SUBROUTINE GAMMA_LAI

    
    

    
    SUBROUTINE GAMMA_P(             &
         DOY_in, tmidh, LAT, LONG,  &                    
         PPFD, D_PPFD, gam_P        &
         )
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: DOY_in 

      
      
      REAL, INTENT(IN)  :: tmidh
      REAL, INTENT(IN)  ::  LAT    
      REAL, INTENT(IN)  ::  LONG   
      REAL, INTENT(IN)  ::  PPFD   
      REAL, INTENT(IN)  ::  D_PPFD 
      REAL, INTENT(OUT) ::  gam_P  


      
      INTEGER :: DOY                 
      REAL :: HOUR                   
      REAL :: AAA, BBB
      REAL :: SIN_solarangle         
      REAL :: Ptoa, Pac, PHI

      

      
      
      DOY = DOY_in
      HOUR = tmidh + long/15.
      IF ( HOUR .LT. 0.0 ) THEN
         HOUR = HOUR + 24.0
         DOY  = DOY - 1
      ENDIF

      
      
      Pac = PPFD
 
      
      CALL SOLARANGLE( DOY, HOUR, LAT, SIN_solarangle )

      
      IF ( SIN_solarangle .LE. 0.0 ) THEN
         
         gam_P = 0.0
      ELSE
         
         
         
         Ptoa = 3000.0 + 99.0 * COS( 2.*3.14*(DOY-10.)/365. )
         
         
         
         PHI = Pac/(SIN_solarangle * Ptoa)
         
         
         BBB = 1. + 0.0005*( D_PPFD-400. )
         AAA = 2.46 * BBB * PHI - 0.9 * (PHI**2)
         gam_P = SIN_solarangle * AAA

      ENDIF
      
      
      
      IF (SIN_solarangle .LE. 0.0175 .AND. gam_P .GT. 0.1) THEN
         gam_P = 0.1
      ENDIF


    END SUBROUTINE GAMMA_P

    
    

    
    SUBROUTINE GAMMA_A( i_spc, LAIp, LAIc, TSTLEN, D_TEMP, gam_A )
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      


      IMPLICIT NONE

      

      
      INTEGER, INTENT(IN) :: i_spc
      
      REAL, INTENT(IN) :: D_TEMP
      
      
      REAL, INTENT(IN) :: LAIp, LAIc
      
      REAL, INTENT(IN) ::     TSTLEN
      
      REAL, INTENT(OUT) :: gam_A

      

      
      REAL :: Fnew, Fgro, Fmat, Fold
      
      INTEGER ::  AINDX 
      
      INTEGER :: t 
      
      REAL     ti
      
      
      REAL     tm
      

      REAL     Tt                   
      

      

      
      
      

      IF (      (i_spc==imgn_acto) .OR. (i_spc==imgn_acta) .OR. (i_spc==imgn_form)   &
           .OR. (i_spc==imgn_ch4)  .OR. (i_spc==imgn_no)   .OR. (i_spc==imgn_co)     &
           ) THEN
         AINDX = 1

      ELSE IF ( (i_spc==imgn_myrc) .OR. (i_spc==imgn_sabi) .OR. (i_spc==imgn_limo)   &
           .OR. (i_spc==imgn_3car) .OR. (i_spc==imgn_ocim) .OR. (i_spc==imgn_bpin)   &
           .OR. (i_spc==imgn_apin) .OR. ( i_spc==imgn_omtp)                          &
           ) THEN
         AINDX = 2

      ELSE IF ( (i_spc==imgn_afarn) .OR. (i_spc==imgn_bcar) .OR. (i_spc==imgn_osqt)  &
           ) THEN
         AINDX = 3

      ELSE IF (i_spc==imgn_meoh) THEN
         aindx = 4

      ELSE IF ( (i_spc==imgn_isop) .OR. (i_spc==imgn_mbo) ) THEN
         aindx = 5
      ELSE
         WRITE(mesg,fmt = '("Invalid i_spc in SUBROUTINE GAMMA_A; i_spc = ", I3)') i_spc
         CALL wrf_error_fatal3("<stdin>",2254,&
mesg)
      END IF



      
      t = TSTLEN
      
      
      Tt   = D_TEMP

      
      
      IF (LAIp .EQ. LAIc) THEN
         Fnew = 0.0
         Fgro = 0.1
         Fmat = 0.8
         Fold = 0.1
      ELSEIF (LAIp .GT. LAIc) THEN
         Fnew = 0.0
         Fgro = 0.0
         Fold = ( LAIp-LAIc ) / LAIp
         Fmat = 1.0-Fold
      ELSE 
         
         
         IF (Tt .LE. 303.0) THEN
            
            ti = 5.0 + 0.7*(300.0-Tt)
         ELSE
            
            ti = 2.9
         ENDIF
         
         
         
         tm = 2.3*ti

         
         
         IF (t .LE. ti) THEN
            
            Fnew = 1.0 - (LAIp/LAIc)
         ELSE
            
            Fnew = (ti/t) * ( 1-(LAIp/LAIc) )
         ENDIF

         
         IF (t .LE. tm) THEN
            
            Fmat = LAIp/LAIc
         ELSE
            
            Fmat = (LAIp/LAIc) + ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )
         ENDIF

         Fgro = 1.0 - Fnew - Fmat
         Fold = 0.0

      ENDIF

      
      
      gam_A = Fnew*Anew(AINDX) + Fgro*Agro(AINDX)    &
           + Fmat*Amat(AINDX) + Fold*Aold(AINDX)


    END SUBROUTINE GAMMA_A

    
    

    
    SUBROUTINE SOLARANGLE( DAY, SHOUR, LAT, SIN_solarangle )
      
      
      
      
      
      
      
      
      

      IMPLICIT NONE

      
      INTEGER, INTENT(IN) :: DAY                  
      REAL, INTENT(IN)    :: SHOUR                
      REAL, INTENT(IN)    :: LAT                  
      REAL, INTENT(OUT)   :: SIN_solarangle

      
      REAL    :: sindelta, cosdelta, A, B

      

      sindelta = -SIN(0.40907) * COS( 6.28*(REAL(DAY,KIND(0.))+10.)/365. )
      cosdelta = SQRT(1.-sindelta**2.)

      A = SIN( LAT*D2RAD ) * sindelta
      B = COS( LAT*D2RAD ) * cosdelta

      SIN_solarangle = A + B * COS(2.*PI*(SHOUR-12.)/24.)


    END SUBROUTINE SOLARANGLE




  END SUBROUTINE bio_emissions_megan2

END MODULE module_bioemi_megan2
