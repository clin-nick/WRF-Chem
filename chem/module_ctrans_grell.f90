


MODULE module_ctrans_grell
USE module_cu_gd
USE module_dep_simple
USE module_state_description, only:p_co,p_qv,p_so2,p_hno3,p_hno4,p_n2o5,p_nh3,p_h2o2, &
                              p_o3,p_ora1,p_op1,p_paa,p_sulf,p_so4aj,p_nh4aj,p_no3aj, &
                              p_bc1,p_bc2,p_oc1,p_oc2,p_seas_1,p_seas_2,     &
                              p_seas_3,p_seas_4,p_dms,                       &
                              p_facd,p_mepx,p_pacd
USE module_state_description, only:p_cvasoaX,p_cvasoa1,p_cvasoa2,p_cvasoa3,p_cvasoa4,&
                              p_cvbsoaX,p_cvbsoa1,p_cvbsoa2,p_cvbsoa3,p_cvbsoa4

USE module_state_description, ONLY: mozart_mosaic_4bin_kpp, &
                              mozart_mosaic_4bin_aq_kpp, &
                              mozcart_kpp, t1_mozcart_kpp, &
                              p_hcho, p_c3h6ooh, p_onit, p_mvk, p_macr, &
                              p_etooh, p_prooh, p_acetp, p_mgly, p_mvkooh, &
                              p_onitr, p_isooh, p_ch3oh, p_c2h5oh, &
                              p_glyald, p_hydrald, p_ald, p_isopn, &
                              p_alkooh, p_mekooh, p_tolooh, p_terpooh, &
                              p_xooh, p_ch3cooh, p_hcooh, p_ch3ooh, &
                              p_so4_a01,p_no3_a01,p_smpa_a01,p_smpbb_a01,&
                              p_glysoa_r1_a01,p_glysoa_r2_a01,&
                              p_glysoa_sfc_a01,p_glysoa_nh4_a01,&
                              p_glysoa_oh_a01,&
                              p_asoaX_a01,p_asoa1_a01,p_asoa2_a01,p_asoa3_a01,p_asoa4_a01,&
                              p_bsoaX_a01,p_bsoa1_a01,p_bsoa2_a01,p_bsoa3_a01,p_bsoa4_a01,&
                              p_biog1_c_a01,p_biog1_o_a01,&
                              p_cl_a01,p_co3_a01,p_nh4_a01,p_na_a01,&
                              p_ca_a01,p_oin_a01,p_oc_a01,p_bc_a01,&
                              p_so4_a02,p_no3_a02,p_smpa_a02,p_smpbb_a02,&
                              p_glysoa_r1_a02,p_glysoa_r2_a02,&
                              p_glysoa_sfc_a02,p_glysoa_nh4_a02,&
                              p_glysoa_oh_a02,&
                              p_asoaX_a02,p_asoa1_a02,p_asoa2_a02,p_asoa3_a02,p_asoa4_a02,&
                              p_bsoaX_a02,p_bsoa1_a02,p_bsoa2_a02,p_bsoa3_a02,p_bsoa4_a02,&
                              p_biog1_c_a02,p_biog1_o_a02,&
                              p_cl_a02,p_co3_a02,p_nh4_a02,p_na_a02,&
                              p_ca_a02,p_oin_a02,p_oc_a02,p_bc_a02,&
                              p_so4_a03,p_no3_a03,p_smpa_a03,p_smpbb_a03,&
                              p_glysoa_r1_a03,p_glysoa_r2_a03,&
                              p_glysoa_sfc_a03,p_glysoa_nh4_a03,&
                              p_glysoa_oh_a03,&
                              p_asoaX_a03,p_asoa1_a03,p_asoa2_a03,p_asoa3_a03,p_asoa4_a03,&
                              p_bsoaX_a03,p_bsoa1_a03,p_bsoa2_a03,p_bsoa3_a03,p_bsoa4_a03,&
                              p_biog1_c_a03,p_biog1_o_a03,&
                              p_cl_a03,p_co3_a03,p_nh4_a03,p_na_a03,&
                              p_ca_a03,p_oin_a03,p_oc_a03,p_bc_a03,&
                              p_so4_a04,p_no3_a04,p_smpa_a04,p_smpbb_a04,&
                              p_glysoa_r1_a04,p_glysoa_r2_a04,&
                              p_glysoa_sfc_a04,p_glysoa_nh4_a04,&
                              p_glysoa_oh_a04,&
                              p_asoaX_a04,p_asoa1_a04,p_asoa2_a04,p_asoa3_a04,p_asoa4_a04,&
                              p_bsoaX_a04,p_bsoa1_a04,p_bsoa2_a04,p_bsoa3_a04,p_bsoa4_a04,&
                              p_biog1_c_a04,p_biog1_o_a04,&
                              p_cl_a04,p_co3_a04,p_nh4_a04,p_na_a04,&
                              p_ca_a04,p_oin_a04,p_oc_a04,p_bc_a04

  IMPLICIT NONE


  REAL, PARAMETER :: qcldwtr_cutoff = 1.0e-6 

  
  REAL, PARAMETER :: mwdry = 28.966  
  REAL, PARAMETER :: mwso4 = 96.00   
  REAL, PARAMETER :: mwno3 = 62.0    
  REAL, PARAMETER :: mwnh4 = 18.0985 

  REAL, PARAMETER :: mwoa  = 250.0   

  INTEGER, allocatable :: HLC_ndx(:)
  
CONTAINS


   SUBROUTINE GRELLDRVCT(DT,itimestep,DX,                           &
              rho_phy,RAINCV,chem,                                  &
              U,V,t_phy,moist,dz8w,p_phy,                           &
              XLV,CP,G,r_v,z,cu_co_ten,                             &
              wd_no3,wd_so4,                                        &
              wd_nh4,wd_oa,                                         &
              wd_so2, wd_sulf, wd_hno3, wd_nh3,                     &
              wd_cvasoa, wd_cvbsoa, wd_asoa, wd_bsoa,               &
              k22_shallow,kbcon_shallow,ktop_shallow,xmb_shallow,   &
              ishallow,num_moist,numgas,num_chem,chemopt,scalaropt, &
              conv_tr_wetscav,conv_tr_aqchem,                       &
              ids,ide, jds,jde, kds,kde,                            &
              ims,ime, jms,jme, kms,kme,                            &
              its,ite, jts,jte, kts,kte                             )



   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) ::                               &
                                  numgas,chemopt,scalaropt,     &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte,    &
                                  ishallow,num_chem,num_moist,  &
                                  conv_tr_wetscav,conv_tr_aqchem

   INTEGER,      INTENT(IN   ) :: ITIMESTEP

   REAL,         INTENT(IN   ) :: XLV, R_v
   REAL,         INTENT(IN   ) :: CP,G
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme,num_moist )         ,    &
          INTENT(IN   ) ::                              moist
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                          U,    &
                                                          V,    &
                                                      t_phy,    &
                                                      z,        &
                                                      p_phy,    &
                                                       dz8w,    &
                                                    rho_phy



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(INOUT   ) ::                                   &
                                                    cu_co_ten


   REAL, INTENT(IN   ) :: DT, DX

   REAL, DIMENSION( ims:ime , kms:kme , jms:jme, num_chem ),    &
         INTENT(INOUT) ::                                       &
                                   chem

   REAL, DIMENSION( ims:ime , jms:jme ),                        &
         INTENT(IN) ::        RAINCV, xmb_shallow
   INTEGER, DIMENSION( ims:ime , jms:jme ),                        &
         INTENT(IN) ::   k22_shallow,kbcon_shallow,ktop_shallow



   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: wd_no3,wd_so4
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: wd_nh4,wd_oa, &
                                     wd_so2, wd_sulf, wd_hno3, wd_nh3
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: &
                                    wd_cvasoa,wd_cvbsoa,wd_asoa,wd_bsoa

     real,    dimension (its:ite,kts:kte) ::                    &
        OUTT,OUTQ,OUTQC
     real,    dimension (its:ite)         ::                    &
        pret, ter11

   REAL, DIMENSION (its:ite,num_chem)         :: wetdep_1d
   REAL, DIMENSION (ims:ime,jms:jme,num_chem) :: wetdep_2d

   REAL, DIMENSION( ims:ime , jms:jme ) :: wdi_no3,wdi_so4
   REAL, DIMENSION( ims:ime , jms:jme ) :: wdi_nh4,wdi_oa, &
   wdi_so2, wdi_sulf, wdi_hno3, wdi_nh3
   REAL, DIMENSION( ims:ime , jms:jme ) :: wdi_cvasoa, wdi_cvbsoa, &
                                           wdi_asoa, wdi_bsoa





     real,    dimension (its:ite,kts:kte) ::                    &
        T,TN,q,qo,PO,P,US,VS,hstary
     real,    dimension (its:ite,kts:kte,num_chem) ::                    &
           tracer,tracert,tracert3
     real, dimension (its:ite)            ::                    &
        Z1,PSUR,AAEQ,xmb3
     integer, dimension (its:ite)            ::                    &
        ktop,k23,kbcon3,ktop3

   INTEGER :: ishallow_gf = 0
   INTEGER :: nv,i,j,k,ICLDCK,ipr,jpr,npr
   REAL    :: tcrit,dp,dq,epsilc
   INTEGER :: itf,jtf,ktf,iopt
   epsilc=1.e-30
   wetdep_1d(:,:)   = 0.0
   wetdep_2d(:,:,:) = 0.0







   ipr=0
   jpr=0
   npr=1
   if(p_co.gt.1)npr=p_co
   tcrit=258.
   iopt=0
   itf=MIN(ite,ide-1)
   ktf=MIN(kte,kde-1)
   jtf=MIN(jte,jde-1)


123  continue
     DO 100 J = jts,jtf
     if(j.eq.jpr)print *,'dt = ',dt
     DO I=ITS,ITF
         ktop(i)=0
         PSUR(I)=p_phy(I,kts,J)*.01
         TER11(I)=z(i,kts,j)
         aaeq(i)=0.




         pret(i)=raincv(i,j)/dt
         if(pret(i).le.0.)aaeq(i)=20.
     ENDDO
     if(Ishallow.eq.1)then
     DO I=ITS,ITF
         xmb3(i)=xmb_shallow(i,j)
         ktop3(i)=ktop_shallow(i,j)
         k23(i)=k22_shallow(i,j)
         kbcon3(i)=kbcon_shallow(i,j)
     ENDDO
     else
     DO I=ITS,ITF
         xmb3(i)=0.
         ktop3(i)=0
         k23(i)=0
         kbcon3(i)=0
     ENDDO
     endif
     DO K=kts,ktf
     DO I=ITS,ITF
         po(i,k)=p_phy(i,k,j)*.01
         P(I,K)=PO(i,k)
         US(I,K) =u(i,k,j)
         VS(I,K) =v(i,k,j)
         T(I,K)=t_phy(i,k,j)
         q(I,K)=moist(i,k,j,p_qv)
         IF(Q(I,K).LT.1.E-08)Q(I,K)=1.E-08
     ENDDO
     ENDDO
     do nv=2,num_chem
     DO K=kts,ktf
     DO I=ITS,ITF
         tracer(i,k,nv)=max(epsilc,chem(i,k,j,nv))
         tracert(i,k,nv)=0.
         tracert3(i,k,nv)=0.
     ENDDO
     ENDDO
     ENDDO
     DO K=kts,ktf
     DO I=ITS,ITF
         cu_co_ten(i,k,j)=0.

         if(i.eq.ipr.and.j.eq.jpr)then
          print *,k,pret(i),tracer(i,k,npr),p(i,k),z(i,k,j)
         endif
     ENDDO
     ENDDO




      CALL CUP_ct(ktop,k23,kbcon3,ktop3,xmb3,                  &
           tracer,j,AAEQ,T,Q,TER11,PRET,P,tracert,tracert3,    &
           hstary,DT,PSUR,US,VS,tcrit,wetdep_1d,               &
           xlv,r_v,cp,g,ipr,jpr,npr,num_chem,chemopt,scalaropt,&
           conv_tr_wetscav,conv_tr_aqchem,                     &
           ishallow,numgas,ids,ide, jds,jde, kds,kde,          &
           ims,ime, jms,jme, kms,kme,                          &
           its,ite, jts,jte, kts,kte                           )

            do nv=2,num_chem
            DO I=its,itf
              if(pret(i).le.0.)then
                 DO K=kts,ktf
                   tracert(i,k,nv)=0.
                 ENDDO
              endif
             enddo
             enddo
      if(ishallow.eq.1)then
        CALL neg_check_ct('shallow',pret,ktop3,epsilc,dt,tracer,tracert3,iopt,num_chem,   &
                        its,ite,kts,kte,itf,ktf,ipr,jpr,npr,j)
        do nv=2,num_chem
            DO I=its,itf
               DO K=kts,ktf
                  chem(i,k,j,nv)=max(epsilc,chem(i,k,j,nv)+tracert3(i,k,nv)*dt)
               enddo
            ENDDO
       ENDDO
     endif



      CALL neg_check_ct('deep',pret,ktop,epsilc,dt,tracer,tracert,iopt,num_chem,  &
                        its,ite,kts,kte,itf,ktf,ipr,jpr,npr,j)
        do nv=2,num_chem
            DO I=its,itf
              if(pret(i).gt.0.)then
                 DO K=kts,ktf
                   chem(i,k,j,nv)=max(epsilc,chem(i,k,j,nv)+tracert(i,k,nv)*dt)
                   if(nv.eq.npr)then
                        cu_co_ten(i,k,j)=tracert(i,k,npr)*dt

                   endif
                 ENDDO
              else
                 DO K=kts,ktf
                   tracert(i,k,nv)=0.
                   if(nv.eq.npr)cu_co_ten(i,k,j)=0.
                 enddo
              endif
            ENDDO
       ENDDO

    wetdep_2d(its:ite,j,:) = wetdep_1d(its:ite,:)

 100    continue

   if(chemopt > 0) then  
   
   
   wdi_no3(:,:) = 0.0
   wdi_so4(:,:) = 0.0
   wdi_nh4(:,:) = 0.0
   wdi_oa(:,:)  = 0.0
   wdi_so2(:,:)  = 0.0
   wdi_sulf(:,:) = 0.0
   wdi_hno3(:,:) = 0.0
   wdi_nh3(:,:)  = 0.0

   wdi_cvasoa(:,:) = 0.0
   wdi_cvbsoa(:,:) = 0.0
   wdi_asoa(:,:)   = 0.0
   wdi_bsoa(:,:)   = 0.0
   
   

   
   if (chemopt == mozart_mosaic_4bin_kpp .OR. &
       chemopt == mozart_mosaic_4bin_aq_kpp) then

     wdi_no3(its:ite,jts:jte) = wdi_no3(its:ite,jts:jte) + &
       (wetdep_2d(its:ite,jts:jte,p_no3_a01) + &
        wetdep_2d(its:ite,jts:jte,p_no3_a02) + &
        wetdep_2d(its:ite,jts:jte,p_no3_a03) + &
        wetdep_2d(its:ite,jts:jte,p_no3_a04))*dt*0.001/mwno3 

     wdi_so4(its:ite,jts:jte) = wdi_so4(its:ite,jts:jte) + &
       (wetdep_2d(its:ite,jts:jte,p_so4_a01) + &
        wetdep_2d(its:ite,jts:jte,p_so4_a02) + &
        wetdep_2d(its:ite,jts:jte,p_so4_a03) + &
        wetdep_2d(its:ite,jts:jte,p_so4_a04))*dt*0.001/mwso4 

     wdi_nh4(its:ite,jts:jte) = wdi_nh4(its:ite,jts:jte) + &
       (wetdep_2d(its:ite,jts:jte,p_nh4_a01) + &
        wetdep_2d(its:ite,jts:jte,p_nh4_a02) + &
        wetdep_2d(its:ite,jts:jte,p_nh4_a03) + &
        wetdep_2d(its:ite,jts:jte,p_nh4_a04))*dt*0.001/mwnh4 

     if (chemopt == mozart_mosaic_4bin_kpp) then

       wdi_oa(its:ite,jts:jte) = wdi_oa(its:ite,jts:jte) + &
         (wetdep_2d(its:ite,jts:jte,p_oc_a01) + &
          wetdep_2d(its:ite,jts:jte,p_oc_a02) + &
          wetdep_2d(its:ite,jts:jte,p_oc_a03) + &
          wetdep_2d(its:ite,jts:jte,p_oc_a04) + &
          wetdep_2d(its:ite,jts:jte,p_smpa_a01) + &
          wetdep_2d(its:ite,jts:jte,p_smpa_a02) + &
          wetdep_2d(its:ite,jts:jte,p_smpa_a03) + &
          wetdep_2d(its:ite,jts:jte,p_smpa_a04) + &
          wetdep_2d(its:ite,jts:jte,p_smpbb_a01) + &
          wetdep_2d(its:ite,jts:jte,p_smpbb_a02) + &
          wetdep_2d(its:ite,jts:jte,p_smpbb_a03) + &
          wetdep_2d(its:ite,jts:jte,p_smpbb_a04) + &
          wetdep_2d(its:ite,jts:jte,p_biog1_c_a01) + &
          wetdep_2d(its:ite,jts:jte,p_biog1_c_a02) + &
          wetdep_2d(its:ite,jts:jte,p_biog1_c_a03) + &
          wetdep_2d(its:ite,jts:jte,p_biog1_c_a04) + &
          wetdep_2d(its:ite,jts:jte,p_biog1_o_a01) + &
          wetdep_2d(its:ite,jts:jte,p_biog1_o_a02) + &
          wetdep_2d(its:ite,jts:jte,p_biog1_o_a03) + &
          wetdep_2d(its:ite,jts:jte,p_biog1_o_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r1_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r1_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r1_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r1_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r2_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r2_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r2_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r2_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_oh_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_oh_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_oh_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_oh_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_nh4_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_nh4_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_nh4_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_nh4_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_sfc_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_sfc_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_sfc_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_sfc_a04))*dt*0.001/mwoa 
     endif

     if (chemopt == mozart_mosaic_4bin_aq_kpp) then

       wdi_asoa(its:ite,jts:jte) = wdi_asoa(its:ite,jts:jte) + &
         (wetdep_2d(its:ite,jts:jte,p_asoaX_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoaX_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoaX_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoaX_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoa1_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoa1_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoa1_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoa1_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoa2_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoa2_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoa2_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoa2_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoa3_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoa3_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoa3_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoa3_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoa4_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoa4_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoa4_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoa4_a04))*dt*0.001/150.0 

       wdi_bsoa(its:ite,jts:jte) = wdi_bsoa(its:ite,jts:jte) + &
         (wetdep_2d(its:ite,jts:jte,p_bsoaX_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoaX_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoaX_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoaX_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa1_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa1_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa1_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa1_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa2_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa2_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa2_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa2_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa3_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa3_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa3_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa3_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa4_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa4_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa4_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa4_a04))*dt*0.001/180.0 

       wdi_cvasoa(its:ite,jts:jte) = wdi_cvasoa(its:ite,jts:jte) + &
         (wetdep_2d(its:ite,jts:jte,p_cvasoaX) + &
          wetdep_2d(its:ite,jts:jte,p_cvasoa1) + &
          wetdep_2d(its:ite,jts:jte,p_cvasoa2) + &
          wetdep_2d(its:ite,jts:jte,p_cvasoa3) + &
          wetdep_2d(its:ite,jts:jte,p_cvasoa4))*dt 

       wdi_cvbsoa(its:ite,jts:jte) = wdi_cvbsoa(its:ite,jts:jte) + &
         (wetdep_2d(its:ite,jts:jte,p_cvbsoaX) + &
          wetdep_2d(its:ite,jts:jte,p_cvbsoa1) + &
          wetdep_2d(its:ite,jts:jte,p_cvbsoa2) + &
          wetdep_2d(its:ite,jts:jte,p_cvbsoa3) + &
          wetdep_2d(its:ite,jts:jte,p_cvbsoa4))*dt 

       wdi_oa(its:ite,jts:jte) = wdi_oa(its:ite,jts:jte) + &
         (wetdep_2d(its:ite,jts:jte,p_oc_a01) + &
          wetdep_2d(its:ite,jts:jte,p_oc_a02) + &
          wetdep_2d(its:ite,jts:jte,p_oc_a03) + &
          wetdep_2d(its:ite,jts:jte,p_oc_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoaX_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoaX_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoaX_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoaX_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoa1_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoa1_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoa1_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoa1_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoa2_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoa2_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoa2_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoa2_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoa3_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoa3_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoa3_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoa3_a04) + &
          wetdep_2d(its:ite,jts:jte,p_asoa4_a01) + &
          wetdep_2d(its:ite,jts:jte,p_asoa4_a02) + &
          wetdep_2d(its:ite,jts:jte,p_asoa4_a03) + &
          wetdep_2d(its:ite,jts:jte,p_asoa4_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoaX_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoaX_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoaX_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoaX_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa1_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa1_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa1_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa1_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa2_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa2_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa2_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa2_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa3_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa3_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa3_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa3_a04) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa4_a01) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa4_a02) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa4_a03) + &
          wetdep_2d(its:ite,jts:jte,p_bsoa4_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r1_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r1_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r1_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r1_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r2_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r2_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r2_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_r2_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_oh_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_oh_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_oh_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_oh_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_nh4_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_nh4_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_nh4_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_nh4_a04) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_sfc_a01) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_sfc_a02) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_sfc_a03) + &
          wetdep_2d(its:ite,jts:jte,p_glysoa_sfc_a04))*dt*0.001/mwoa 
     endif

   else
     if (p_no3aj .gt.1) wdi_no3(its:ite,jts:jte) = wdi_no3(its:ite,jts:jte) + wetdep_2d(its:ite,jts:jte,p_no3aj)*dt*0.001/mwno3 
     if (p_so4aj .gt.1) wdi_so4(its:ite,jts:jte) = wdi_so4(its:ite,jts:jte) + wetdep_2d(its:ite,jts:jte,p_so4aj)*dt*0.001/mwso4 
   endif
   
   if (p_hno3  .gt.1) wdi_hno3(its:ite,jts:jte) = wdi_hno3(its:ite,jts:jte) + wetdep_2d(its:ite,jts:jte,p_hno3)*dt              
   if (p_hno4  .gt.1) wdi_hno3(its:ite,jts:jte) = wdi_hno3(its:ite,jts:jte) + wetdep_2d(its:ite,jts:jte,p_hno4)*dt              
   
   if (p_sulf  .gt.1) wdi_sulf(its:ite,jts:jte) = wdi_sulf(its:ite,jts:jte) + wetdep_2d(its:ite,jts:jte,p_sulf)*dt              
   if (p_so2   .gt.1) wdi_so2(its:ite,jts:jte) = wdi_so2(its:ite,jts:jte) + wetdep_2d(its:ite,jts:jte,p_so2)*dt               

   if (p_nh3   .gt.1) wdi_nh3(its:ite,jts:jte) = wdi_nh3(its:ite,jts:jte) + wetdep_2d(its:ite,jts:jte,p_nh3)*dt               
   
   
   
   wd_no3(its:ite,jts:jte) = wd_no3(its:ite,jts:jte) + wdi_no3(its:ite,jts:jte) 
   wd_so4(its:ite,jts:jte) = wd_so4(its:ite,jts:jte) + wdi_so4(its:ite,jts:jte) 
   wd_nh4(its:ite,jts:jte) = wd_nh4(its:ite,jts:jte) + wdi_nh4(its:ite,jts:jte) 
   wd_oa (its:ite,jts:jte) = wd_oa (its:ite,jts:jte) + wdi_oa (its:ite,jts:jte) 
   wd_so2 (its:ite,jts:jte) = wd_so2 (its:ite,jts:jte) + wdi_so2 (its:ite,jts:jte) 
   wd_sulf (its:ite,jts:jte) = wd_sulf (its:ite,jts:jte) + wdi_sulf (its:ite,jts:jte) 
   wd_hno3 (its:ite,jts:jte) = wd_hno3 (its:ite,jts:jte) + wdi_hno3 (its:ite,jts:jte) 
   wd_nh3 (its:ite,jts:jte) = wd_nh3 (its:ite,jts:jte) + wdi_nh3 (its:ite,jts:jte) 

   wd_asoa(its:ite,jts:jte)   = wd_asoa(its:ite,jts:jte)   + wdi_asoa(its:ite,jts:jte)   
   wd_bsoa(its:ite,jts:jte)   = wd_bsoa(its:ite,jts:jte)   + wdi_bsoa(its:ite,jts:jte)   
   wd_cvasoa(its:ite,jts:jte) = wd_cvasoa(its:ite,jts:jte) + wdi_cvasoa(its:ite,jts:jte) 
   wd_cvbsoa(its:ite,jts:jte) = wd_cvbsoa(its:ite,jts:jte) + wdi_cvbsoa(its:ite,jts:jte) 
   
   endif
   
   END SUBROUTINE GRELLDRVCT


   SUBROUTINE CUP_ct(ktop,k23,kbcon3,ktop3,xmb3,tracer,J,AAEQ,T,Q,Z1,  &
              PRE,P,tracert,tracert3,hstary,DTIME,PSUR,US,VS,TCRIT,    &
              wetdep,xl,rv,cp,g,ipr,jpr,npr,num_chem,chemopt,scalaropt,&
              conv_tr_wetscav,conv_tr_aqchem,                          &
              ishallow,numgas,ids,ide, jds,jde, kds,kde,               &
              ims,ime, jms,jme, kms,kme,                               &
              its,ite, jts,jte, kts,kte                                )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        num_chem,ids,ide, jds,jde, kds,kde,                            &
        ims,ime, jms,jme, kms,kme,scalaropt,conv_tr_wetscav,           &
        conv_tr_aqchem,ishallow,                                       &
        its,ite, jts,jte, kts,kte,ipr,jpr,npr,chemopt,numgas
     integer, intent (in   )              ::                           &
        j
  
  
  
  
  
     real,    dimension (its:ite,kts:kte,num_chem)                              &
        ,intent (inout  )                   ::                           &
        tracert,tracer,tracert3
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        pre
     integer,    dimension (its:ite)                                   &
         ,intent (inout  )                   ::                        &
          ktop,k23,kbcon3,ktop3
     integer,    dimension (its:ite)     ::                            &
        kbcon
  
  
  
  
  
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        T,P,US,VS,HSTARY
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        Z1,PSUR,AAEQ,xmb3


       real                                                            &
        ,intent (in   )                   ::                           &
        dtime,tcrit,xl,cp,rv,g


     real,    dimension (its:ite,1:3) ::                         &
        edtc










  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

     real,    dimension (its:ite,kts:kte) ::                           &
        he,hes,qes,z,pwdper,                                           &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &

        dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,zd3,      &
        clw_all3,cd3,pw3,zu3,                                          &

  
  

        cd,cdd,cdd3,scr1,DELLAH,DELLAQ,DELLAT,DELLAQC

  
  
     real,    dimension (its:ite) ::                                   &
       edt,HKB,QKB,edt3,          &
       XMB,PWAV,PWEV,BU,cap_max,cap_max_increment
     real,    dimension (its:ite,kts:kte,num_chem)       ::             &
        tr3_c,tr3_up,tr3_pw
     real,    dimension (its:ite,num_chem)         ::                   &
        trkb3
     real,    dimension (its:ite,kts:kte,num_chem)       ::             &
        tr_c,tr_up,tr_dd,tr_dd3,tre_cup,tr_pw,tr_pwd
     real,    dimension (its:ite,num_chem)         ::                   &
        trkb,wetdep
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,K22,KB,JMIN,kstabi,kstabm,                     &   
       ierr,KBMAX,ierr5

     integer                              ::                           &
       ki,I,K,KK
     real                                 ::                           &
      day,dz,mbdt,entr_rate,radius,entrd_rate,mentr_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,            &
      dh,cap_maxs,entr_rate3,mentr_rate3

     integer :: itf,jtf,ktf

     itf=MIN(ite,ide-1)
     ktf=MIN(kte,kde-1)
     jtf=MIN(jte,jde-1)


      day=86400.



      radius=12000.




      entr_rate=.2/radius
      entr_rate3=.2/200.



      mentrd_rate=0.
      mentr_rate=entr_rate
      mentr_rate3=entr_rate3



      do k=kts,ktf
      do i=its,itf
        cd(i,k)=0.1*entr_rate
        cd3(i,k)=entr_rate3
        cdd(i,k)=0.
        cdd3(i,k)=0.
        clw_all(i,k)=0.
      enddo
      enddo
      do ki=1,num_chem
      do k=kts,ktf
      do i=its,itf
        tr_dd3(i,k,ki)=0.
      enddo
      enddo
      enddo




      edtmax=.8
      edtmin=.2



      depth_min=500.




      cap_maxs=175.

      DO 7 i=its,itf
        kbmax(i)=1
        cap_max_increment(i)=0.
        edt(i)=0.
        edt3(i)=0.
        kstabm(i)=ktf-1
        IERR(i)=0
        IERR5(i)=0
        if(ktop3(i).lt.2)ierr5(i)=20
        if(xmb3(i).eq.0.)ierr5(i)=21
 7    CONTINUE
      do i=its,itf
          cap_max(i)=cap_maxs
          do k=kts,kte
           zd3(i,k)=0.
          enddo
      enddo



      zkbmax=4000.



      zcutdown=3000.



      z_detr=1250.

      mbdt=dtime*4.E-03



      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,0,xl,cp,   &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)



      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev_tr(tracer,tre_cup,num_chem,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
      if(ierr(i).eq.0)then

      do k=kts,ktf-2
        if(z_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue




      do k=kts,ktf
        if(z_cup(i,k).gt.z_detr+z1(i))then
          kdet(i)=k
          go to 26
        endif
      enddo
 26   continue

      endif
      enddo





      CALL cup_MAXIMI(HE_CUP,3,KBMAX,K22,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
       DO 36 i=its,itf
         IF(ierr(I).eq.0.)THEN
         IF(K22(I).GE.KBMAX(i))ierr(i)=2
         endif
 36   CONTINUE




      DO i=its,itf
        if(aaeq(i).ne.0.)then
           ierr(i)=20
        endif
      enddo



      call cup_kbcon(cap_max_increment,1,k22,kbcon,he_cup,hes_cup, &
           ierr,kbmax,p_cup,cap_max, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)



      CALL cup_minimi(HEs_cup,Kbcon,kstabm,kstabi,ierr,  &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
      IF(ierr(I).eq.0.)THEN
        if(kstabm(i)-1.gt.kstabi(i))then
           do k=kstabi(i),kstabm(i)-1
             cd(i,k)=cd(i,k-1)+1.5*entr_rate
             if(cd(i,k).gt.10.0*entr_rate)cd(i,k)=10.0*entr_rate
           enddo
        ENDIF
      ENDIF
      ENDDO



      call cup_up_he(k22,hkb,z_cup,cd,mentr_rate,he_cup,hc, &
           kbcon,ierr,dby,he,hes_cup, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)



      call cup_ktop(1,dby,kbcon,ktop,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      DO 37 i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(z_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,ktf
              if(z_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 go to 37
              endif
              enddo
         endif
 37   CONTINUE



      call cup_minimi(HEs_cup,K22,kzdown,JMIN,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      DO 100 i=its,ite
      IF(ierr(I).eq.0.)THEN
           if(jmin(i).le.3)then
             ierr(i)=9
             go to 100
           endif




101   continue
      if(jmin(i)-1.lt.KDET(I))kdet(i)=jmin(i)-1
      if(jmin(i).ge.Ktop(I)-1)jmin(i)=ktop(i)-2
      ki=jmin(i)
      hcd(i,ki)=hes_cup(i,ki)
      DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
      dh=dz*(HCD(i,Ki)-hes_cup(i,ki))
      dh=0.

      do k=ki-1,1,-1
         hcd(i,k)=hes_cup(i,jmin(i))
         DZ=Z_cup(i,K+1)-Z_cup(i,K)
         dh=dh+dz*(HCD(i,K)-hes_cup(i,k))
         if(dh.gt.0.)then
           jmin(i)=jmin(i)-1
           if(jmin(i).gt.3)then
             go to 101
           else if(jmin(i).le.3)then
             ierr(i)=9
             go to 100
           endif
         endif
       enddo

         IF(JMIN(I).LE.3)then
            ierr(i)=4
         endif

      ENDIF
100   continue




      do i=its,itf
      IF(ierr(I).eq.0.)THEN
           if(jmin(i).le.3)then
             ierr(i)=9
           endif
      IF(-z_cup(I,KBCON(I))+z_cup(I,KTOP(I)).LT.depth_min)then
            ierr(i)=6
      endif
      endif
      enddo




      call cup_up_nms(zu,z_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)




      call cup_dd_nms(zd,z_cup,cdd,mentrd_rate,jmin,ierr, &
           0,kdet,z1,                 &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)



      call cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,mentrd_rate, &
           jmin,ierr,he,dbyd,he_cup,  &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)



      call cup_dd_moisture(zd,hcd,hes_cup,qcd,qes_cup, &
           pwd,q_cup,z_cup,cdd,mentrd_rate,jmin,ierr,gamma_cup, &
           pwev,bu,qrcd,q,he,t_cup,1,xl, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)



      call cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav, &
           kbcon,ktop,cd,dby,mentr_rate,clw_all,      &
           q,GAMMA_cup,zu,qes_cup,k22,q_cup,xl, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)



      call cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
           pwev,edtmax,edtmin,3,edtc, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
        do i=its,itf
         if(ierr(i).eq.0)then
         edt(i)=edtc(i,2)
         endif
        enddo



        pwdper=0.
        do i=its,itf

          if(ierr(i).gt.0)pre(i)=0.
          if(ierr(i).eq.0)then
          xmb(i)=pre(i)/(pwav(i)+edt(i)*pwev(i))



          if(i.eq.ipr.and.j.eq.jpr)then
            print *,'xmb,edt,pwav = ',xmb(i),edt(i),pwav(i)
            print *,'k,pwdper(i,k),pw,pwd(i,k)',z1(i)
          endif
          do k=1,ktop(i)
          pwdper(i,k)=-edt(i)*pwd(i,k)/pwav(i)
          if(i.eq.ipr.and.j.eq.jpr)then
            print *,k,pwdper(i,k),pw(i,k),pwd(i,k)
          endif
          enddo
          endif
        enddo









       call cup_up_tracer(ierr,tcrit,t,pre,z_cup,p,tracer,tre_cup,tr_up,tr_pw, &
                   tr_c,hstary,pw,clw_all,kbcon,ktop,cd,mentr_rate,zu,k22,&
                   numgas,chemopt,scalaropt,conv_tr_wetscav,conv_tr_aqchem, &
                   ids,ide, jds,jde, kds,kde, &
                   ims,ime, jms,jme, kms,kme, &
                   its,ite, jts,jte, kts,kte, &
                   ipr,jpr,j,npr,num_chem,'deep')




       if(ishallow.eq.1)then
      call cup_up_nms(zu3,z_cup,mentr_rate3,cd3,kbcon3,ktop3,ierr5,k23, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
       call cup_up_tracer(ierr5,tcrit,t,pre,z_cup,p,tracer,tre_cup,tr3_up,tr3_pw, &
                   tr3_c,hstary,pw3,clw_all3,kbcon3,ktop3,cd3,mentr_rate3,zu3,k23,&
                   numgas,chemopt,scalaropt,conv_tr_wetscav,conv_tr_aqchem, &
                   ids,ide, jds,jde, kds,kde, &
                   ims,ime, jms,jme, kms,kme, &
                   its,ite, jts,jte, kts,kte, &
                   ipr,jpr,j,npr,num_chem,'shallow')
      call cup_dellas_tr(ierr5,z_cup,p_cup,tr_dd3,edt3,zd3,cdd3,      &
           tracer,tracert3,j,mentrd_rate,zu3,g,xmb3,                &
           cd3,tr3_up,ktop3,k23,kbcon3,mentr_rate3,jmin,tre_cup,kdet, &
           k23,ipr,jpr,npr,'shallow',num_chem,                      &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte                             )
       endif

       call cup_dd_tracer(ierr,z_cup,qrcd,tracer,tre_cup,tr_up,tr_dd, &
                    tr_pw,tr_pwd,jmin,cdd,mentrd_rate,zd,pwdper,wetdep,xmb,k22, &
                    numgas,num_chem,ids,ide, jds,jde, kds,kde, &
                    ims,ime, jms,jme, kms,kme, &
                    its,ite, jts,jte, kts,kte)
       if(j.eq.jpr)print *,'called dd_tracer'



      if(j.eq.jpr)then
        i=ipr
        print *,'in 250 loop ',edt(ipr),ierr(ipr)

        print *,k22(I),kbcon(i),ktop(i),jmin(i)
        print *,edt(i)
        do k=kts,ktf
          print *,k,z(i,k),he(i,k),hes(i,k)
        enddo
        do k=1,ktop(i)+1
          print *,zu(i,k),zd(i,k),pw(i,k),pwd(i,k)
        enddo
        print *,'tr_up(i,k,6),tr_dd(i,k,6),tr_pw(i,k,6),tr_pwd(i,k,6)'
        do k=1,ktop(i)+1
          print *,tr_up(i,k,npr),tr_dd(i,k,npr),tr_pw(i,k,npr),tr_pwd(i,k,npr)
        enddo
        endif






      call cup_dellabot_tr(ipr,jpr,tre_cup,ierr,z_cup,p,tr_dd,edt, &
           zd,cdd,tracer,tracert,j,mentrd_rate,z,g,xmb, &
           num_chem,ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)




      call cup_dellas_tr(ierr,z_cup,p_cup,tr_dd,edt,zd,cdd,      &
           tracer,tracert,j,mentrd_rate,zu,g,xmb,                &
           cd,tr_up,ktop,k22,kbcon,mentr_rate,jmin,tre_cup,kdet, &
           k22,ipr,jpr,npr,'deep',num_chem,                      &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte                             )
       if(j.eq.jpr)then
        i=ipr
        do k=kts,ktf
          print *,k,tracer(i,k,npr),tracert(i,k,npr)
        enddo
       endif









   END SUBROUTINE CUP_CT

   SUBROUTINE cup_dellabot_tr(ipr,jpr,tre_cup,ierr,z_cup,p_cup,  &
              tr_dd,edt,zd,cdd,tracer,tracert,j,mentrd_rate,z,g,xmb,     &
              num_chem,ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        num_chem,ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ipr,jpr
  
  
  
     real,    dimension (its:ite,kts:kte,1:num_chem)                   &
        ,intent (out  )                   ::                           &
        tracert
     real,    dimension (its:ite,kts:kte,1:num_chem)                   &
        ,intent (in   )                   ::                           &
        tre_cup,tracer,tr_dd
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,p_cup,zd,cdd,z
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt,xmb
     real                                                              &
        ,intent (in   )                   ::                           &
        g,mentrd_rate
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr




      integer i
      real detdo,detdo1,detdo2,entdo,dp,dz,subin,                      &
      totmas

     integer :: itf, ktf, nv, npr
     npr=24
     itf=MIN(ite,ide-1)
     ktf=MIN(kte,kde-1)


      if(j.eq.jpr)print *,'in cup dellabot '
      tracert=0.
      do 100 i=its,itf
      if(ierr(i).ne.0)go to 100
      dz=z_cup(i,2)-z_cup(i,1)
      DP=100.*(p_cup(i,1)-P_cup(i,2))
      detdo1=edt(i)*zd(i,2)*CDD(i,1)*DZ
      detdo2=edt(i)*zd(i,1)
      entdo=edt(i)*zd(i,2)*mentrd_rate*dz
      subin=-EDT(I)*zd(i,2)
      detdo=detdo1+detdo2-entdo+subin
      do nv=2,num_chem
      tracert(I,1,nv)=(detdo1*.5*(tr_dd(i,1,nv)+tr_dd(i,2,nv)) &
                 +detdo2*tr_dd(i,1,nv) &
                 +subin*tre_cup(i,2,nv) &
                 -entdo*tracer(i,1,nv))*g/dp*xmb(i)
      enddo
      if(j.eq.jpr.and.i.eq.ipr)print *,'in cup dellabot ',tracert(I,1,npr), &
        detdo1,detdo2,subin,entdo,tr_dd(i,1,npr),tr_dd(i,2,npr),tracer(i,1,npr)
 100  CONTINUE

   END SUBROUTINE cup_dellabot_tr


   SUBROUTINE cup_dellas_tr(ierr,z_cup,p_cup,tr_dd,edt,zd,cdd,             &
              tracer,tracert,j,mentrd_rate,zu,g,xmb,                       &
              cd,tr_up,ktop,k22,kbcon,mentr_rate,jmin,tre_cup,kdet,kpbl,   &
              ipr,jpr,npr,name,num_chem,                                   &
              ids,ide, jds,jde, kds,kde,                                   &
              ims,ime, jms,jme, kms,kme,                                   &
              its,ite, jts,jte, kts,kte                                    )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        num_chem,ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ipr,jpr,npr
  
  
  
     real,    dimension (its:ite,kts:kte,1:num_chem)                  &
        ,intent (inout  )                   ::                           &
        tracert
     real,    dimension (its:ite,kts:kte,1:num_chem)                  &
        ,intent (in  )                   ::                           &
        tr_up,tr_dd,tre_cup,tracer
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,p_cup,zd,cdd,cd,zu
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt,xmb
     real                                                              &
        ,intent (in   )                   ::                           &
        g,mentrd_rate,mentr_rate
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22,jmin,kdet,kpbl
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name




      integer i,k,nv
      real detdo1,detdo2,entdo,dp,dz,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas

     integer :: itf, ktf, kstart

     itf=MIN(ite,ide-1)
     ktf=MIN(kte,kde-1)
      kstart=kts+1
      if(name.eq.'shallow')kstart=kts


      i=ipr
      if(j.eq.jpr)then
        print *,'in dellas kpbl(i),k22(i),kbcon(i),ktop(i),jmin(i)'
        print *,kpbl(i),k22(i),kbcon(i),ktop(i),jmin(i)
      endif
       do nv=2,num_chem
       DO K=kstart,kte
       do i=its,itf
          tracert(i,k,nv)=0.
       enddo
       enddo
       enddo

       DO 100 k=kts+1,ktf-1
       DO 100 i=its,ite
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100




         DZ=Z_cup(I,K+1)-Z_cup(I,K)
         detdo=edt(i)*CDD(i,K)*DZ*ZD(i,k+1)
         entdo=edt(i)*mentrd_rate*dz*zd(i,k+1)
         subin=zu(i,k+1)-zd(i,k+1)*edt(i)
         entup=0.
         detup=0.
         if(k.ge.kbcon(i).and.k.lt.ktop(i))then
            entup=mentr_rate*dz*zu(i,k)
            detup=CD(i,K+1)*DZ*ZU(i,k)
         endif
         subdown=(zu(i,k)-zd(i,k)*edt(i))
         entdoj=0.
         entupk=0.
         detupk=0.

         if(k.eq.jmin(i))then
         entdoj=edt(i)*zd(i,k)
         endif

         if(k.eq.k22(i)-1)then
         entupk=zu(i,kpbl(i))
         endif

         if(k.gt.kdet(i))then
            detdo=0.
         endif

         if(k.eq.ktop(i)-0)then
         detupk=zu(i,ktop(i))
         subin=0.
         endif
         if(k.lt.kbcon(i))then
            detup=0.
         endif



         totmas=subin-subdown+detup-entup-entdo+ &
                 detdo-entupk-entdoj+detupk
          if(j.eq.jpr.and.i.eq.ipr)print *,'k,totmas,sui,sud = ',k, &
          totmas,subin,subdown




         if(abs(totmas).gt.1.e-6)then
           print *,'*********************',i,j,k,totmas,name
           print *,kpbl(i),k22(i),kbcon(i),ktop(i)




        CALL wrf_error_fatal3("<stdin>",1387,&
'cup_dellas_tr: TOTMAS > CRITICAL VALUE')
         endif
         dp=100.*(p_cup(i,k-1)-p_cup(i,k))
         do nv=2,num_chem


         tracert(i,k,nv)=(subin*tracer(i,k+1,nv) &
                    -subdown*tracer(i,k,nv) &
                    +detup*.5*(tr_up(i,K+1,nv)+tr_up(i,K,nv)) &
                    +detdo*.5*(tr_dd(i,K+1,nv)+tr_dd(i,K,nv)) &
                    -entup*tracer(i,k,nv) &
                    -entdo*tracer(i,k,nv) &
                    -entupk*tre_cup(i,k22(i),nv) &
                    -entdoj*tre_cup(i,jmin(i),nv) &
                    +detupk*tr_up(i,ktop(i),nv) &
                     )*g/dp*xmb(i)
         enddo
       if(i.eq.ipr.and.j.eq.jpr)then
         print *,k,tracert(i,k,npr),subin*tre_cup(i,k+1,npr),subdown*tre_cup(i,k,npr), &
                   detdo*.5*(tr_dd(i,K+1,npr)+tr_dd(i,K,npr))
         print *,k,detup*.5*(tr_up(i,K+1,npr)+tr_up(i,K,npr)),detupk*tr_up(i,ktop(i),npr), &
                entup*tracer(i,k,npr),entdo*tracer(i,k,npr)
         print *,k,entupk*tre_cup(i,k,npr),detupk,tr_up(i,ktop(i),npr)
       endif

 100  CONTINUE

   END SUBROUTINE cup_dellas_tr
   SUBROUTINE cup_env_clev_tr(tracer,tre_cup,num_chem,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      implicit none
     integer                                                           &
        ,intent (in   )                   ::                           &
        num_chem,ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, dimension (its:ite)                                      &
        ,intent (in)                      ::                           &
        ierr

     real,    dimension (its:ite,kts:kte,1:num_chem)                   &
        ,intent (in   )                   ::                           &
        tracer
     real,    dimension (its:ite,kts:kte,1:num_chem)                   &
        ,intent (out  )                   ::                           &
        tre_cup




     integer                              ::                           &
       i,k,nv,itf,ktf
     itf=MIN(ite,ide-1)
     ktf=MIN(kte,kde-1)
      do nv=2,num_chem
      do k=kts+1,ktf
      do i=its,ite
        if(ierr(i).eq.0)then
        tre_cup(i,k,nv)=.5*(tracer(i,k-1,nv)+tracer(i,k,nv))
        endif
      enddo
      enddo
      enddo
      do nv=2,num_chem
      do i=its,ite
        if(ierr(i).eq.0)then
        tre_cup(i,kts,nv)=tracer(i,kts,nv)
        endif
      enddo
      enddo


END subroutine cup_env_clev_tr


   SUBROUTINE  cup_up_tracer(ierr,tcrit,t,pre,z_cup,p,tracer,tre_cup,tr_up, &
                tr_pw,tr_c,hstary,cupclw,clw_all,kbcon,ktop,cd,mentr_rate,zu,k22,  &
                          numgas,chemopt,scalaropt,conv_tr_wetscav,conv_tr_aqchem, &
                          ids,ide, jds,jde, kds,kde, &
                          ims,ime, jms,jme, kms,kme, &
                          its,ite, jts,jte, kts,kte,ipr,jpr,j,npr,num_chem,name)


  USE module_state_description, only: RADM2SORG,RADM2SORG_AQ,RACMSORG_AQ,RACMSORG_KPP,   &
                                      RADM2SORG_KPP,RACM_ESRLSORG_KPP,RACM_SOA_VBS_KPP,  &
                                      RADM2SORG_AQCHEM,RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP, &
                                      RACM_SOA_VBS_AQCHEM_KPP,                           &
                                      RACM_SOA_VBS_HET_KPP,                              &
                                      CB05_SORG_VBS_AQ_KPP
  USE module_ctrans_aqchem
  USE module_input_chem_data, only: get_last_gas
  USE module_HLawConst, only: HLC
        implicit none




      INTEGER, PARAMETER :: NGAS = 12  
      INTEGER, PARAMETER :: NAER = 36  
      INTEGER, PARAMETER :: NLIQS = 41 



      INTEGER, PARAMETER :: LSO2    =  1  
      INTEGER, PARAMETER :: LHNO3   =  2  
      INTEGER, PARAMETER :: LN2O5   =  3  
      INTEGER, PARAMETER :: LCO2    =  4  
      INTEGER, PARAMETER :: LNH3    =  5  
      INTEGER, PARAMETER :: LH2O2   =  6  
      INTEGER, PARAMETER :: LO3     =  7  
      INTEGER, PARAMETER :: LFOA    =  8  
      INTEGER, PARAMETER :: LMHP    =  9  
      INTEGER, PARAMETER :: LPAA    = 10  
      INTEGER, PARAMETER :: LH2SO4  = 11  
      INTEGER, PARAMETER :: LHCL    = 12  



      INTEGER, PARAMETER :: LSO4AKN  =  1  
      INTEGER, PARAMETER :: LSO4ACC  =  2  
      INTEGER, PARAMETER :: LSO4COR  =  3  
      INTEGER, PARAMETER :: LNH4AKN  =  4  
      INTEGER, PARAMETER :: LNH4ACC  =  5  
      INTEGER, PARAMETER :: LNO3AKN  =  6  
      INTEGER, PARAMETER :: LNO3ACC  =  7  
      INTEGER, PARAMETER :: LNO3COR  =  8  
      INTEGER, PARAMETER :: LORGAAKN =  9  
      INTEGER, PARAMETER :: LORGAACC = 10  
      INTEGER, PARAMETER :: LORGPAKN = 11  
      INTEGER, PARAMETER :: LORGPACC = 12  
      INTEGER, PARAMETER :: LORGBAKN = 13  
      INTEGER, PARAMETER :: LORGBACC = 14  
      INTEGER, PARAMETER :: LECAKN   = 15  
      INTEGER, PARAMETER :: LECACC   = 16  
      INTEGER, PARAMETER :: LPRIAKN  = 17  
      INTEGER, PARAMETER :: LPRIACC  = 18  
      INTEGER, PARAMETER :: LPRICOR  = 19  
      INTEGER, PARAMETER :: LNAAKN   = 20  
      INTEGER, PARAMETER :: LNAACC   = 21  
      INTEGER, PARAMETER :: LNACOR   = 22  
      INTEGER, PARAMETER :: LCLAKN   = 23  
      INTEGER, PARAMETER :: LCLACC   = 24  
      INTEGER, PARAMETER :: LCLCOR   = 25  
      INTEGER, PARAMETER :: LNUMAKN  = 26  
      INTEGER, PARAMETER :: LNUMACC  = 27  
      INTEGER, PARAMETER :: LNUMCOR  = 28  
      INTEGER, PARAMETER :: LSRFAKN  = 29  
      INTEGER, PARAMETER :: LSRFACC  = 30  
      INTEGER, PARAMETER :: LNACL    = 31  
      INTEGER, PARAMETER :: LCACO3   = 32  
      INTEGER, PARAMETER :: LMGCO3   = 33  
      INTEGER, PARAMETER :: LA3FE    = 34  
      INTEGER, PARAMETER :: LB2MN    = 35  
      INTEGER, PARAMETER :: LK       = 36  





   

     integer                                                           &
        ,intent (in   )                   ::                           &
                         ids,ide, jds,jde, kds,kde,scalaropt,   &
                         numgas, conv_tr_wetscav,conv_tr_aqchem,        &
                         num_chem,ims,ime, jms,jme, kms,kme,chemopt,   &
                         its,ite, jts,jte, kts,kte,ipr,jpr,j,npr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,cd,zu,p,hstary,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        cupclw,clw_all
     real,    dimension (its:ite,kts:kte,1:num_chem)                  &
        ,intent (inout  )                   ::                           &
        tr_up,tr_c,tr_pw
     real,    dimension (its:ite,kts:kte,1:num_chem)                  &
        ,intent (in  )                   ::                           &
        tre_cup,tracer
     real,    dimension (its:ite)                              &
        ,intent (in  )                   ::                           &
        pre

  
     real                                                              &
        ,intent (in   )                   ::                           &
        mentr_rate,tcrit
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
   

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name


      real :: conc_equi,conc_mxr,partialp,taucld
      real :: HLCnst1, HLCnst2
      real :: HLCnst3, HLCnst4
      real :: HLCnst5, HLCnst6
      real :: kh, dk1s, dk2s, heff, tfac

     integer                              ::                           &
        iall,i,k,iwd,nv, HLndx
     real                                 ::                           &
        trcc,trch,dh,c0,dz,radius,airm,dens
     integer                              ::                           &
       itf,ktf,iaer,igas
     
     real                                 ::                           &
       frac_so4(4), frac_no3(4), frac_nh4(4), tot_so4, tot_nh4, tot_no3
     
     
     
     real aq_gas_ratio
     
     
     
     real precip 
     real, dimension (ngas) :: gas,gaswdep
     real, dimension (naer) :: aerosol,aerwdep
     real, dimension (nliqs) :: liquid
     real hpwdep
     real alfa0,alfa2,alfa3 

     real, parameter :: hion = 1.e-5
     real, parameter :: hion_inv = 1./hion
     real, parameter :: HL_t0 = 298.
     
     itf=MIN(ite,ide-1)
     ktf=MIN(kte,kde-1)


        iall=0
        c0=.002
        iwd=0






      if(name.eq.'shallow')c0=0.


      do nv=2,num_chem
      do k=kts,ktf
      do i=its,itf
        
        tr_pw(i,k,nv)=0.
        
        if(ierr(i).eq.0)tr_up(i,k,nv)=tre_cup(i,k,nv)
        
        tr_c(i,k,nv)=0.
      enddo
      enddo
      enddo
      
      do nv=2,num_chem
      do i=its,itf
      if(ierr(i).eq.0.)then
      do k=k22(i),kbcon(i)-1
        tr_up(i,k,nv)=tre_cup(i,k22(i),nv)
      enddo
      endif
      enddo
      enddo
      
      DO 100 k=kts+1,ktf-1
      DO 100 i=its,itf
      
      IF(ierr(i).ne.0)GO TO 100
      IF(K.Lt.KBCON(I))GO TO 100
      IF(K.Gt.KTOP(I)+1)GO TO 100
      
      DZ=Z_cup(i,K)-Z_cup(i,K-1)
      
      if(cupclw(i,k).le.0.)cupclw(i,k)=0.
      if(clw_all(i,k).le.0.)clw_all(i,k)=0.
      
      
      
      
      
      
      do nv=2,num_chem
        if(i.eq.ipr.and.j.eq.jpr.and.nv.eq.npr)print *,k,tr_up(i,K-1,nv),tr_up(i,K,nv),tr_pw(i,k-1,nv),clw_all(i,k),cupclw(i,k)
        tr_up(i,K,nv)=(tr_up(i,K-1,nv)*(1.-.5*CD(i,K)*DZ)+mentr_rate* &
          DZ*tracer(i,K-1,nv))/(1.+mentr_rate*DZ-.5*cd(i,k)*dz)
        tr_up(i,k,nv)=max(1.e-16,tr_up(i,K,nv))
      enddo
      
      
      
      

      if ((chemopt .EQ. RADM2SORG .OR. chemopt .EQ. RADM2SORG_AQ .OR. chemopt .EQ. RACMSORG_AQ .OR. & 
           chemopt .EQ. RACMSORG_KPP .OR. chemopt .EQ. RADM2SORG_KPP .OR. chemopt .EQ. RACM_ESRLSORG_KPP .OR. & 
           chemopt .EQ. RACM_SOA_VBS_KPP .OR. chemopt .EQ. RADM2SORG_AQCHEM .OR. chemopt .EQ. RACMSORG_AQCHEM_KPP .OR. &
           chemopt .EQ. CB05_SORG_VBS_AQ_KPP .OR.                                           &
           chemopt .EQ. RACM_SOA_VBS_HET_KPP .OR.   &
           chemopt .EQ. RACM_ESRLSORG_AQCHEM_KPP .OR. chemopt .EQ. RACM_SOA_VBS_AQCHEM_KPP) &
           .and. (conv_tr_aqchem == 1)) then
        
        
        
        
        
        
        dens = 0.1*p(i,k)/t(i,k)*mwdry/8.314472 
        
        
        airm = 1000.0*dens*dz/mwdry 
        
        
        
        GASWDEP = 0.0
        AERWDEP = 0.0
        HPWDEP  = 0.0
        
        
        
        
        precip = 0.0 
        
        alfa0 = 0.0
        alfa2 = 0.0
        alfa3 = 0.0
        
        
        
        
        gas(:) = 0.0
        
        gas(lco2)   = 380.0e-6
        
        gas(lso2)   = tr_up(i,k,p_so2)*1.0e-6
        gas(lhno3)  = tr_up(i,k,p_hno3)*1.0e-6
        gas(ln2o5)  = tr_up(i,k,p_n2o5)*1.0e-6
        gas(lnh3)   = tr_up(i,k,p_nh3)*1.0e-6
        gas(lh2o2)  = tr_up(i,k,p_h2o2)*1.0e-6
        gas(lo3)    = tr_up(i,k,p_o3)*1.0e-6
        gas(lh2so4) = tr_up(i,k,p_sulf)*1.0e-6
        if (chemopt==CB05_SORG_VBS_AQ_KPP) then
           gas(lfoa)   = tr_up(i,k,p_facd)*1.0e-6
           gas(lmhp)   = tr_up(i,k,p_mepx)*1.0e-6
           gas(lpaa)   = tr_up(i,k,p_pacd)*1.0e-6
        else
           gas(lfoa)   = tr_up(i,k,p_ora1)*1.0e-6
           gas(lmhp)   = tr_up(i,k,p_op1)*1.0e-6
           gas(lpaa)   = tr_up(i,k,p_paa)*1.0e-6
        end if

        
        
        
        
        
        
        aerosol(:) = 0.0
        
        
        
        aerosol(lso4acc) = tr_up(i,k,p_so4aj)*1.0e-9*mwdry/mwso4
        aerosol(lnh4acc) = tr_up(i,k,p_nh4aj)*1.0e-9*mwdry/mwnh4
        aerosol(lno3acc) = tr_up(i,k,p_no3aj)*1.0e-9*mwdry/mwno3
        
        
        taucld = 1800.0
      
        if (clw_all(i,k)*dens .gt. qcldwtr_cutoff) then 
          CALL AQCHEM( &
           t(i,k), &
           p(i,k)*100., &
           taucld, &
           precip, &
           clw_all(i,k)*dens, &
           clw_all(i,k)*dens, &
           airm, &
           ALFA0, &
           ALFA2, &
           ALFA3, &
           GAS, &
           AEROSOL, &
           LIQUID, &
           GASWDEP, &
           AERWDEP, &
           HPWDEP)
        endif
        
        
        
        
        tr_up(i,k,p_so2)  =  gas(lso2)*1.0e6
        tr_up(i,k,p_hno3) =  gas(lhno3)*1.0e6
        tr_up(i,k,p_n2o5) =  gas(ln2o5)*1.0e6
        tr_up(i,k,p_nh3)  =  gas(lnh3)*1.0e6
        tr_up(i,k,p_h2o2) =  gas(lh2o2)*1.0e6
        tr_up(i,k,p_o3)   =  gas(lo3)*1.0e6
        tr_up(i,k,p_sulf) =  gas(lh2so4)*1.0e6
        if (chemopt==CB05_SORG_VBS_AQ_KPP) then
           tr_up(i,k,p_facd) =  gas(lfoa)*1.0e6
           tr_up(i,k,p_mepx)  =  gas(lmhp)*1.0e6
           tr_up(i,k,p_pacd)  =  gas(lpaa)*1.0e6
        else
           tr_up(i,k,p_ora1) =  gas(lfoa)*1.0e6
           tr_up(i,k,p_op1)  =  gas(lmhp)*1.0e6
           tr_up(i,k,p_paa)  =  gas(lpaa)*1.0e6
        end if
        
        
        
        
        tr_up(i,k,p_so4aj) = aerosol(lso4acc)*1.0e9*mwso4/mwdry
        tr_up(i,k,p_nh4aj) = aerosol(lnh4acc)*1.0e9*mwnh4/mwdry
        tr_up(i,k,p_no3aj) = aerosol(lno3acc)*1.0e9*mwno3/mwdry

      else if ((chemopt .EQ. mozart_mosaic_4bin_kpp .OR. &
                chemopt .EQ. mozart_mosaic_4bin_aq_kpp)  &
               .AND. (conv_tr_aqchem == 1)) then

        
        
        

        
        dens = 0.1*p(i,k)/t(i,k)*mwdry/8.314472 

        
        airm = 1000.0*dens*dz/mwdry 

        

        GASWDEP = 0.0
        AERWDEP = 0.0
        HPWDEP  = 0.0

        
        

        precip = 0.0 

        alfa0 = 0.0
        alfa2 = 0.0
        alfa3 = 0.0

        
        

        gas(:) = 0.0

        gas(lco2)   = 380.0e-6

        gas(lso2)   = tr_up(i,k,p_so2)*1.0e-6
        gas(lhno3)  = tr_up(i,k,p_hno3)*1.0e-6
        gas(ln2o5)  = tr_up(i,k,p_n2o5)*1.0e-6
        gas(lnh3)   = tr_up(i,k,p_nh3)*1.0e-6
        gas(lh2o2)  = tr_up(i,k,p_h2o2)*1.0e-6
        gas(lo3)    = tr_up(i,k,p_o3)*1.0e-6
        gas(lfoa)   = tr_up(i,k,p_hcooh)*1.0e-6
        gas(lmhp)   = tr_up(i,k,p_ch3ooh)*1.0e-6
        gas(lpaa)   = tr_up(i,k,p_paa)*1.0e-6
        gas(lh2so4) = tr_up(i,k,p_sulf)*1.0e-6

        
        
        
        
        

        aerosol(:) = 0.0

        

        
        
        frac_so4(:) = 0.25
        frac_nh4(:) = 0.25
        frac_no3(:) = 0.25

        tot_so4     = tr_up(i,k,p_so4_a01)+tr_up(i,k,p_so4_a02)+&
                      tr_up(i,k,p_so4_a03)+tr_up(i,k,p_so4_a04)
        tot_nh4     = tr_up(i,k,p_nh4_a01)+tr_up(i,k,p_nh4_a02)+&
                      tr_up(i,k,p_nh4_a03)+tr_up(i,k,p_nh4_a04)
        tot_no3     = tr_up(i,k,p_no3_a01)+tr_up(i,k,p_no3_a02)+&
                      tr_up(i,k,p_no3_a03)+tr_up(i,k,p_no3_a04)

        if (tot_so4 > 0.0) then
          frac_so4(1) = tr_up(i,k,p_so4_a01) / tot_so4
          frac_so4(2) = tr_up(i,k,p_so4_a02) / tot_so4
          frac_so4(3) = tr_up(i,k,p_so4_a03) / tot_so4
          frac_so4(4) = tr_up(i,k,p_so4_a04) / tot_so4
          aerosol(lso4acc) = tot_so4 *1.0e-9*mwdry/mwso4
        end if

        if (tot_nh4 > 0.0) then
          frac_nh4(1) = tr_up(i,k,p_nh4_a01) / tot_nh4
          frac_nh4(2) = tr_up(i,k,p_nh4_a02) / tot_nh4
          frac_nh4(3) = tr_up(i,k,p_nh4_a03) / tot_nh4
          frac_nh4(4) = tr_up(i,k,p_nh4_a04) / tot_nh4
          aerosol(lnh4acc) = tot_nh4 *1.0e-9*mwdry/mwnh4
        end if

        if (tot_no3 > 0.0) then
          frac_no3(1) = tr_up(i,k,p_no3_a01) / tot_no3
          frac_no3(2) = tr_up(i,k,p_no3_a02) / tot_no3
          frac_no3(3) = tr_up(i,k,p_no3_a03) / tot_no3
          frac_no3(4) = tr_up(i,k,p_no3_a04) / tot_no3
          aerosol(lno3acc) = tot_no3 *1.0e-9*mwdry/mwno3
        end if

        
        taucld = 1800.0

        if (clw_all(i,k)*dens .gt. qcldwtr_cutoff) then 
          CALL AQCHEM( &
           t(i,k), &
           p(i,k)*100., &
           taucld, &
           precip, &
           clw_all(i,k)*dens, &
           clw_all(i,k)*dens, &
           airm, &
           ALFA0, &
           ALFA2, &
           ALFA3, &
           GAS, &
           AEROSOL, &
           LIQUID, &
           GASWDEP, &
           AERWDEP, &
           HPWDEP)
      endif
      
        
        

        tr_up(i,k,p_so2)    =  gas(lso2)*1.0e6
        tr_up(i,k,p_hno3)   =  gas(lhno3)*1.0e6
        tr_up(i,k,p_n2o5)   =  gas(ln2o5)*1.0e6
        tr_up(i,k,p_nh3)    =  gas(lnh3)*1.0e6
        tr_up(i,k,p_h2o2)   =  gas(lh2o2)*1.0e6
        tr_up(i,k,p_o3)     =  gas(lo3)*1.0e6
        tr_up(i,k,p_hcooh)  =  gas(lfoa)*1.0e6
        tr_up(i,k,p_ch3ooh) =  gas(lmhp)*1.0e6
        tr_up(i,k,p_paa)    =  gas(lpaa)*1.0e6
        tr_up(i,k,p_sulf)   =  gas(lh2so4)*1.0e6

        
        

        tr_up(i,k,p_so4_a01) = aerosol(lso4acc) * frac_so4(1) * 1.0e9*mwso4/mwdry
        tr_up(i,k,p_so4_a02) = aerosol(lso4acc) * frac_so4(2) * 1.0e9*mwso4/mwdry
        tr_up(i,k,p_so4_a03) = aerosol(lso4acc) * frac_so4(3) * 1.0e9*mwso4/mwdry
        tr_up(i,k,p_so4_a04) = aerosol(lso4acc) * frac_so4(4) * 1.0e9*mwso4/mwdry

        tr_up(i,k,p_nh4_a01) = aerosol(lnh4acc) * frac_nh4(1) * 1.0e9*mwnh4/mwdry
        tr_up(i,k,p_nh4_a02) = aerosol(lnh4acc) * frac_nh4(2) * 1.0e9*mwnh4/mwdry
        tr_up(i,k,p_nh4_a03) = aerosol(lnh4acc) * frac_nh4(3) * 1.0e9*mwnh4/mwdry
        tr_up(i,k,p_nh4_a04) = aerosol(lnh4acc) * frac_nh4(4) * 1.0e9*mwnh4/mwdry

        tr_up(i,k,p_no3_a01) = aerosol(lno3acc) * frac_no3(1) * 1.0e9*mwno3/mwdry
        tr_up(i,k,p_no3_a02) = aerosol(lno3acc) * frac_no3(2) * 1.0e9*mwno3/mwdry
        tr_up(i,k,p_no3_a03) = aerosol(lno3acc) * frac_no3(3) * 1.0e9*mwno3/mwdry
        tr_up(i,k,p_no3_a04) = aerosol(lno3acc) * frac_no3(4) * 1.0e9*mwno3/mwdry

      endif
      


      if (conv_tr_wetscav == 1) then
        
        
          tfac = (HL_t0 - t(i,k))/(t(i,k)*HL_t0)
          do nv=2,num_chem
            
            aq_gas_ratio = 0.0
is_moz_chm: if( chemopt == MOZCART_KPP .or. chemopt == T1_MOZCART_KPP .or. &
                chemopt == MOZART_MOSAIC_4BIN_KPP .or. chemopt == MOZART_MOSAIC_4BIN_AQ_KPP ) then
              HLndx = HLC_ndx(nv)
              if( HLndx /= 0 ) then
                HLCnst1 = HLC(HLndx)%hcnst(1) ; HLCnst2 = HLC(HLndx)%hcnst(2)
                HLCnst3 = HLC(HLndx)%hcnst(3) ; HLCnst4 = HLC(HLndx)%hcnst(4)
                HLCnst5 = HLC(HLndx)%hcnst(5) ; HLCnst6 = HLC(HLndx)%hcnst(6)
                kh = HLCnst1 * exp( HLCnst2* tfac )
                if( HLCnst3 /= 0. ) then
                  dk1s = HLCnst3 * exp( HLCnst4* tfac )
                else
                  dk1s = 0.
                endif
                if( HLCnst5 /= 0. ) then
                  dk2s = HLCnst5 * exp( HLCnst6* tfac )
                else
                  dk2s = 0.
                endif
                if( nv /= p_nh3 ) then
                  heff = kh*(1. + dk1s*hion_inv*(1. + dk2s*hion_inv))
                else
                  heff = kh*(1. + dk1s*hion/dk2s)
                endif
                aq_gas_ratio = moz_aq_frac(t(i,k), clw_all(i,k)*dens, heff )
              endif
            else is_moz_chm
            
            

            
            if (nv .eq. p_h2o2)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 8.33e+04, 7379.)
            if (nv .eq. p_hno3)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 2.6e+06, 8700.)
            if (nv .eq. p_hcho)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 6.30e+03, 6425.)
            if (nv .eq. p_ch3ooh)  aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.11e+02, 5241.)
            if (nv .eq. p_c3h6ooh) aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 2.20e+02, 5653.)
            if (nv .eq. p_paa)     aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 8.37e+02, 5308.)
            if (nv .eq. p_hno4)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 1.2e+04, 6900.) 
            if (nv .eq. p_onit)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 1.00e+03, 6000.)
            if (nv .eq. p_mvk)     aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 1.7e-03, 0.)
            if (nv .eq. p_macr)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 1.70e-03, 0.)
            if (nv .eq. p_etooh)   aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.36e+02, 5995.)
            if (nv .eq. p_prooh)   aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.36e+02, 5995.)
            if (nv .eq. p_acetp)   aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.36e+02, 5995.)
            if (nv .eq. p_mgly)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.71e+03, 7541.)
            if (nv .eq. p_mvkooh)  aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 2.6e+06, 8700.)
            if (nv .eq. p_onitr)   aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 7.51e+03, 6485.)
            if (nv .eq. p_isooh)   aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 2.6e+06, 8700.)
            if (nv .eq. p_ch3oh)   aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 2.20e+02, 4934.)
            if (nv .eq. p_c2h5oh)  aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 2.00e+02, 6500.)
            if (nv .eq. p_glyald)  aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 4.14e+04, 4630.)
            if (nv .eq. p_hydrald) aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 7.00e+01, 6000.)
            if (nv .eq. p_ald)     aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 1.14e+01, 6267.)
            if (nv .eq. p_isopn)   aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 1.00e+01, 0.)
            if (nv .eq. p_alkooh)  aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.11e+02, 5241.)
            if (nv .eq. p_mekooh)  aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.11e+02, 5241.)
            if (nv .eq. p_tolooh)  aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.11e+02, 5241.)
            if (nv .eq. p_terpooh) aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 3.11e+02, 5241.)
            if (nv .eq. p_nh3)     aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 7.40e+01, 3400.)
            if (nv .eq. p_xooh)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 90.5, 5607.)
            if (nv .eq. p_ch3cooh) aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 4.1e3, 6300.)
            if (nv .eq. p_so2)     aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 1.2, 3100.)
            if (nv .eq. p_sulf)    aq_gas_ratio = aq_frac(p(i,k)*100., t(i,k), clw_all(i,k)*dens, 1e+11, 0.) 
            ENDIF is_moz_chm
            






            if (nv.gt.numgas) aq_gas_ratio = 0.5

            if (nv.eq.p_so4aj) aq_gas_ratio = 1.0
            if (nv.eq.p_nh4aj) aq_gas_ratio = 1.0
            if (nv.eq.p_no3aj) aq_gas_ratio = 1.0

            if (nv.eq.p_bc1 .or. nv.eq.p_oc1 .or. nv.eq.p_dms) aq_gas_ratio=0.
            if (nv.eq.p_sulf .or. nv.eq.p_seas_1 .or. nv.eq.p_seas_2) aq_gas_ratio=1.
            if (nv.eq.p_seas_3 .or. nv.eq.p_seas_4) aq_gas_ratio=1.
            if (nv.eq.p_bc2 .or. nv.eq.p_oc2) aq_gas_ratio=0.8
            
            if (aq_gas_ratio > 0.0) then
              tr_c(i,k,nv)  = aq_gas_ratio*tr_up(i,k,nv) 
              trch          = tr_up(i,k,nv)-tr_c(i,k,nv) 
              trcc          = (tr_up(i,k,nv)-trch)/(1.+c0*dz*zu(i,k)) 
              tr_pw(i,k,nv) = c0*dz*trcc*zu(i,k) 
              tr_up(i,k,nv) = trcc+trch 
            endif
            
          enddo
        
     endif 
        
      
100   CONTINUE

END subroutine cup_up_tracer




REAL FUNCTION aq_frac(p, T, q, Kh_298, dHoR)

  REAL, INTENT(IN)  :: p,           & 
                       T,           & 
                       q,           & 
                       Kh_298,      & 
                       dHoR           

  REAL, PARAMETER   :: Rgas = 8.31446 

  
  REAL              :: Kh, tr_air, tr_aq


  
  Kh      = Kh_298 * exp ( dHoR * ( 1.0/T - 1.0/298 ) ) * 101.325

  
  tr_air  = 1 / (Rgas * T)

  
  tr_aq   = Kh * (q / 1000.0)

  aq_frac = min( 1.0, max( 0.0, tr_aq / (tr_aq + tr_air) ) )

END FUNCTION aq_frac

REAL FUNCTION moz_aq_frac(T, q, heff )

  REAL, INTENT(IN)  :: T,           & 
                       q,           & 
                       heff           

  REAL, PARAMETER   :: Rgas = 8.31446 

  
  REAL              :: tr_air, tr_aq

  
  tr_air  = 1. / (Rgas * T)

  
  tr_aq   = heff * (q / 1000.0)

  moz_aq_frac = min( 1.0, max( 0.0, tr_aq / (tr_aq + tr_air) ) )

END FUNCTION moz_aq_frac



  SUBROUTINE  cup_dd_tracer(ierr,z_cup,qrcd,tracer,tre_cup,tr_up,tr_dd, &
                          tr_pw,tr_pwd,jmin,cdd,entr,zd,pwdper,wetdep,xmb,k22, &
                          numgas,num_chem,ids,ide, jds,jde, kds,kde, &
                          ims,ime, jms,jme, kms,kme, &
                          its,ite, jts,jte, kts,kte)

  USE module_model_constants, only: mwdry

        implicit none



  
   

     integer                                                           &
        ,intent (in   )                   ::                           &
                               numgas,num_chem,                        &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                            &
       pwdper,zd,cdd,qrcd,z_cup 
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                            &
       xmb
     real,    dimension (its:ite,kts:kte,1:num_chem)                   &
        ,intent (inout  )                   ::                         &
        tr_dd,tr_pwd,tr_up
     real,    dimension (its:ite,kts:kte,1:num_chem)                   &
        ,intent (in  )                   ::                            &
        tre_cup,tracer,tr_pw
     real,    dimension (its:ite,1:num_chem)                           &
        ,intent (out  )                   ::                           &
         wetdep

  
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
   
  
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,k22



     integer                              ::                           &
        iall,i,k,nv,ki,kj
     real                                 ::                           &
        dh,c0,dz,radius
     integer                              ::                           &
       itf,ktf

     real                                 ::                           &
       evaporate,condensate

      itf=MIN(ite,ide-1)
      ktf=MIN(kte,kde-1)
      
      
      tr_pwd(its:ite,kts:kte,1:num_chem) = 0.0
      
      
      wetdep(its:ite,1:num_chem) = 0.0
      
      
      
      
      do i=its,itf
        
        IF(ierr(I).eq.0)then
            


          do nv=1,num_chem
            
            do k=ktf,kts,-1 
              
              
              condensate = tr_pw(i,k,nv)
              
              do kj=k,kts,-1 
                
                evaporate = condensate*pwdper(i,kj)
                
                tr_pwd(i,kj,nv) = tr_pwd(i,kj,nv) + evaporate
                
                condensate = max(0.0,condensate - evaporate)
              enddo
              
            enddo
            
            
            
            
            do k=kts,ktf
              wetdep(i,nv) = wetdep(i,nv) + tr_pw(i,k,nv) - tr_pwd(i,k,nv)
            enddo
            
            wetdep(i,nv) = max(0.0,wetdep(i,nv))
            
          enddo
          
        endif
        
      enddo
      
      
      
      do i=its,itf
        
        IF(ierr(I).eq.0)then
          
          do nv=1,num_chem
            tr_dd(i,jmin(i):ktf,nv)=tre_cup(i,jmin(i):ktf,nv) 
          enddo
          
          do ki=jmin(i)-1,1,-1
            DZ=Z_cup(i,ki+1)-Z_cup(i,ki)
            do nv=1,num_chem
              tr_dd(i,ki,nv)=(tr_dd(i,ki+1,nv)*(1.-.5*CDD(i,ki)*DZ) &
               +entr*DZ*tracer(i,ki,nv) &
                )/(1.+entr*DZ-.5*CDD(i,ki)*DZ)
            enddo
            
          enddo
          
        endif
        
      enddo
      
      
      
      do i=its,itf
        IF(ierr(I).eq.0)then
          do nv=1,num_chem
            do k=kts,ktf
              tr_dd(i,k,nv) = tr_dd(i,k,nv) + tr_pwd(i,k,nv)
            enddo
          enddo
        endif
      enddo
      
      
      
      
      do nv=1,num_chem
        if (nv <= numgas) then
          do i=its,itf
           if(ierr(I).eq.0)then
            wetdep(i,nv)=wetdep(i,nv)*xmb(i)/mwdry 
            wetdep(i,nv) = max(0.0,wetdep(i,nv))
           endif
          enddo
        else
          do i=its,itf
           if(ierr(I).eq.0)then
            wetdep(i,nv)=wetdep(i,nv)*xmb(i) 
            wetdep(i,nv) = max(0.0,wetdep(i,nv))
           endif
          enddo
        endif
      enddo

END subroutine cup_dd_tracer


   SUBROUTINE neg_check_ct(name,pret,ktop,epsilc,dt,q,outq,iopt,num_chem,    &
                           its,ite,kts,kte,itf,ktf,ipr,jpr,npr,j)

   INTEGER,      INTENT(IN   ) ::   iopt,num_chem,its,ite,kts,kte,itf,ktf,ipr,jpr,npr,j

     real, dimension (its:ite,kts:kte,num_chem  )          ,                  &
      intent(inout   ) ::                                                     &
       q,outq
     real, dimension (its:ite  )          ,                                   &
      intent(in      ) ::                                                     &
       pret
     integer, dimension (its:ite  )          ,                                &
      intent(in   ) ::                                                        &
      ktop
     real                                                                     &
        ,intent (in  )                   ::                                   &
        dt,epsilc
     real :: tracermin,tracermax,thresh,qmem,qmemf,qmem2,qtest,qmem1
     character *(*) name
     integer :: nv, i, k







      thresh=epsilc

      if(iopt.eq.0)then
      do nv=2,num_chem
      do 100 i=its,itf
         if(pret(i).le.0.and.name.eq.'deep')go to 100
         tracermin=q(i,kts,nv)
         tracermax=q(i,kts,nv)
         do k=kts+1,kte-1
           tracermin=min(tracermin,q(i,k,nv))
           tracermax=max(tracermax,q(i,k,nv))
         enddo
         tracermin=max(tracermin,thresh)
         qmemf=1.



         do k=kts,ktop(i)



            qmem=outq(i,k,nv)



            if(qmem.lt.0.)then
               qtest=q(i,k,nv)+outq(i,k,nv)*dt
               if(qtest.lt.tracermin)then



                    qmem1=outq(i,k,nv)
                    qmem2=(tracermin-q(i,k,nv))/dt
                    qmemf=min(qmemf,qmem2/qmem1)
                    if(qmemf.gt.1.)print *,'something wrong in negct_1',qmem2,qmem1
                    if(i.eq.ipr.and.j.eq.jpr.and.nv.eq.npr)then
                      print *,k,qtest,qmem2,qmem1,qmemf
                    endif
                    qmemf=max(qmemf,0.)
               endif
            endif
         enddo
         do k=kts,ktop(i)
            outq(i,k,nv)=outq(i,k,nv)*qmemf
         enddo



         qmemf=1.
         do k=kts,ktop(i)



            qmem=outq(i,k,nv)



            if(qmem.gt.0.)then
               qtest=q(i,k,nv)+outq(i,k,nv)*dt
               if(qtest.gt.tracermax)then



                    qmem1=outq(i,k,nv)
                    qmem2=(tracermax-q(i,k,nv))/dt
                    qmemf=min(qmemf,qmem2/qmem1)
                    if(qmemf.gt.1.)print *,'something wrong in negct_2',qmem2,qmem1
                    if(i.eq.ipr.and.j.eq.jpr.and.nv.eq.npr)then
                      print *,'2',k,qtest,qmem2,qmem1,qmemf
                    endif
                    qmemf=max(qmemf,0.)
               endif
            endif
         enddo
         do k=kts,ktop(i)
            outq(i,k,nv)=outq(i,k,nv)*qmemf
         enddo
 100  continue
      enddo



      elseif(iopt.eq.1)then
      do i=its,itf
      qmemf=1.
      do k=kts,ktop(i)
      do nv=2,num_chem



         qmem=outq(i,k,nv)



         if(qmem.lt.0.)then
         qtest=q(i,k,nv)+outq(i,k,nv)*dt
         if(qtest.lt.thresh)then



           qmem1=outq(i,k,nv)
           qmem2=(thresh-q(i,k,nv))/dt
           qmemf=min(qmemf,qmem2/qmem1)
           qmemf=max(0.,qmemf)
         endif
         endif
      enddo
      enddo
      do nv=2,num_chem
      do k=kts,ktop(i)
         outq(i,k,nv)=outq(i,k,nv)*qmemf
      enddo
      enddo
      enddo
      endif

   END SUBROUTINE neg_check_ct

   SUBROUTINE conv_tr_wetscav_init( numgas, num_chem )

   use module_state_description, only : param_first_scalar
   use module_scalar_tables, only : chem_dname_table
   use module_chem_utilities, only : UPCASE
   use module_HLawConst

   integer, intent(in) :: numgas, num_chem




   integer :: m, m1
   integer :: astat
   character(len=64) :: HL_tbl_name
   character(len=64) :: wrf_spc_name

is_allocated : &
   if( .not. allocated(HLC_ndx) ) then



     allocate( HLC_ndx(num_chem),stat=astat )
     if( astat /= 0 ) then
       call wrf_error_fatal3("<stdin>",2463,&
"conv_tr_wetscav_init: failed to allocate HLC_ndx")
     endif
     HLC_ndx(:) = 0
     do m = param_first_scalar,numgas
       wrf_spc_name = chem_dname_table(1,m)
       call upcase( wrf_spc_name )
       do m1 = 1,nHLC
         HL_tbl_name = HLC(m1)%name
         call upcase( HL_tbl_name )
         if( trim(HL_tbl_name) == trim(wrf_spc_name) ) then
           HLC_ndx(m) = m1
           exit
         endif
       end do
     end do
   endif is_allocated

   END SUBROUTINE conv_tr_wetscav_init


END MODULE module_ctrans_grell
