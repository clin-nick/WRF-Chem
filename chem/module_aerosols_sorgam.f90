MODULE module_aerosols_sorgam

  USE module_state_description
  USE module_data_radm2
  USE module_data_sorgam
  USE module_radm


      IMPLICIT NONE

CONTAINS
    SUBROUTINE sorgam_driver (id,ktau,dtstep,t_phy,moist,aerwrf,p8w,    &
               t8w,alt,p_phy,chem,rho_phy,dz8w,z,z_at_w,                &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,vcsulf_old,    &
               vdrog3,                                                  &
               kemit,                                                   &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )

   INTEGER,      INTENT(IN   )    ::                             &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      its,ite, jts,jte, kts,kte, &
                                      kemit,                     &
                                      id,ktau

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
         INTENT(IN ) ::                                   moist

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),         &
         INTENT(INOUT ) ::                                   chem



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(INOUT ) ::                                             &
           h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2,    &
           cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2


   REAL,  DIMENSION(ims:ime,kms:kme-0,jms:jme,ldrog),                  &
           INTENT(IN   ) ::                                            &
                                                  VDROG3               
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,           &
          INTENT(IN   ) ::                                             &
                                                      t_phy,           &
                                                        alt,           &
                                                      p_phy,           &
                                                      dz8w,            &
                                                      z    ,           &
                                              t8w,p8w,z_at_w ,         &
                                                      aerwrf ,         &
                                                    rho_phy
   REAL,  DIMENSION( ims:ime , kms:kme-0 , jms:jme )         ,         &
          INTENT(IN   ) ::                                             &
             vcsulf_old
      REAL,      INTENT(IN   ) ::                                      &
                             dtstep

      REAL drog_in(ldrog)                                    
                                                             
                                                             

      REAL condvap_in(lspcv) 
                             
                             
      REAL rgas
      DATA rgas/8.314510/
      REAL convfac,convfac2


      INTEGER blksize
      PARAMETER (blksize=1)



      INTEGER nspcsda
      PARAMETER (nspcsda=l1ae) 


      INTEGER nacv
      PARAMETER (nacv=lcva) 

      INTEGER ncv
      PARAMETER (ncv=lspcv) 

      REAL cblk(blksize,nspcsda) 
                                   
      REAL soilrat_in
                    
                    
      REAL nitrate_in
                    
      REAL nh3_in
                    
      REAL hcl_in

      REAL vsulf_in

      REAL so4rat_in
                    
      REAL epm25i(blksize),epm25j(blksize),epmcoarse(blksize)
                    
      REAL eeci_in
                    
      REAL eecj_in
                    
      REAL eorgi_in

      REAL eorgj_in
                    
                    
      REAL pres
                    
      REAL temp
                    
      REAL relhum
                    
      REAL ::p(kts:kte),t(kts:kte),rh(kts:kte)




      REAL mwso4
      PARAMETER (mwso4=96.0576)


      REAL mwhno3
      PARAMETER (mwhno3=63.01287)


      REAL mwnh3
      PARAMETER (mwnh3=17.03061)


      REAL mwhcl
      PARAMETER (mwhcl=36.46100)






      REAL mwec
      PARAMETER (mwec=12.0)


      REAL mwaro1
      PARAMETER (mwaro1=150.0)


      REAL mwaro2
      PARAMETER (mwaro2=150.0)


      REAL mwalk1
      PARAMETER (mwalk1=140.0)


      REAL mwalk2
      PARAMETER (mwalk2=140.0)



      REAL mwole1
      PARAMETER (mwole1=140.0)


      REAL mwapi1
      PARAMETER (mwapi1=200.0)


      REAL mwapi2
      PARAMETER (mwapi2=200.0)


      REAL mwlim1
      PARAMETER (mwlim1=200.0)


      REAL mwlim2
      PARAMETER (mwlim2=200.0)




   INTEGER :: i,j,k,l,debug_level




   do l=p_so4aj,num_chem
      do j=jts,jte
         do k=kts,kte
            do i=its,ite
               chem(i,k,j,l)=max(epsilc,chem(i,k,j,l)/alt(i,k,j))
            enddo
         enddo
      enddo
   enddo
      do 100 j=jts,jte
         do 100 i=its,ite
           debug_level=0
            do k=kts,kte
               t(k) = t_phy(i,k,j)
               p(k) = .001*p_phy(i,k,j)
               rh(k) = MIN( 95.,100. * moist(i,k,j,p_qv) /        &
                        (3.80*exp(17.27*(t_phy(i,k,j)-273.)/      &
                        (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j)))   )
               rh(k)=max(.1,0.01*rh(k))
            enddo
            do k=kts,kte






               cblk=0.
               do l=1,ldrog
                  drog_in(l)=0.
               enddo
               do l=1,lspcv
                  condvap_in(l)=0.
               enddo
               convfac = p(k)/rgas/t(k)*1000.
               so4rat_in=(chem(i,k,j,p_sulf)-vcsulf_old(i,k,j))/dtstep*CONVFAC*MWSO4
               soilrat_in = 0.
               nitrate_in =max(epsilc,chem(i,k,j,p_hno3)*convfac*mwhno3)
               nh3_in = max(epsilc,chem(i,k,j,p_nh3)*convfac*mwnh3)
               hcl_in = max(epsilc,chem(i,k,j,p_hcl)*convfac*mwhcl)

               vsulf_in = max(epsilc,chem(i,k,j,p_sulf)*convfac*mwso4)







               



        drog_in(PXYL ) = VDROG3(i,k,j,PXYL )
        drog_in(PTOL ) = VDROG3(i,k,j,PTOL )
        drog_in(PCSL1) = VDROG3(i,k,j,PCSL1)
        drog_in(PCSL2) = VDROG3(i,k,j,PCSL2)
        drog_in(PHC8 ) = VDROG3(i,k,j,PHC8 )
        drog_in(POLI1) = VDROG3(i,k,j,POLI1)
        drog_in(POLI2) = VDROG3(i,k,j,POLI2)
        drog_in(POLI3) = VDROG3(i,k,j,POLI3)
        drog_in(POLT1) = VDROG3(i,k,j,POLT1)
        drog_in(POLT2) = VDROG3(i,k,j,POLT2)
        drog_in(POLT3) = VDROG3(i,k,j,POLT3)

        if(p_lim.eq.1)then

            drog_in(PAPI1) = 0.
            drog_in(PAPI2) = 0.
            drog_in(PAPI3) = 0.
            drog_in(PLIM1) = 0.
            drog_in(PLIM2) = 0.
            drog_in(PLIM3) = 0.
            condvap_in(PSOAAPI1) = 0.
            condvap_in(PSOAAPI2) = 0.
            condvap_in(PSOALIM1) = 0.
            condvap_in(PSOALIM2) = 0.
        elseif(p_lim.gt.1)then

            drog_in(PAPI1) = VDROG3(i,k,j,PAPI1)
            drog_in(PAPI2) = VDROG3(i,k,j,PAPI2)
            drog_in(PAPI3) = VDROG3(i,k,j,PAPI3)
            drog_in(PLIM1) = VDROG3(i,k,j,PLIM1)
            drog_in(PLIM2) = VDROG3(i,k,j,PLIM2)
            drog_in(PLIM3) = VDROG3(i,k,j,PLIM3)
            condvap_in(PSOAAPI1) = max(epsilc,cvapi1(i,k,j))
            condvap_in(PSOAAPI2) = max(epsilc,cvapi2(i,k,j))
            condvap_in(PSOALIM1) = max(epsilc,cvlim1(i,k,j))
            condvap_in(PSOALIM2) = max(epsilc,cvlim2(i,k,j))
        endif
        condvap_in(PSOAARO1) = max(epsilc,cvaro1(i,k,j))
        condvap_in(PSOAARO2) = max(epsilc,cvaro2(i,k,j))
        condvap_in(PSOAALK1) = max(epsilc,cvalk1(i,k,j))
        condvap_in(PSOAOLE1) = max(epsilc,cvole1(i,k,j))
      cblk(1,VORGARO1J) =   chem(i,k,j,p_orgaro1j)
      cblk(1,VORGARO1I) =   chem(i,k,j,p_orgaro1i)
      cblk(1,VORGARO2J) =   chem(i,k,j,p_orgaro2j)
      cblk(1,VORGARO2I) =   chem(i,k,j,p_orgaro2i)
      cblk(1,VORGALK1J) =   chem(i,k,j,p_orgalk1j)
      cblk(1,VORGALK1I) =   chem(i,k,j,p_orgalk1i)
      cblk(1,VORGOLE1J) =   chem(i,k,j,p_orgole1j)
      cblk(1,VORGOLE1I) =   chem(i,k,j,p_orgole1i)
      cblk(1,VORGBA1J ) =   chem(i,k,j,p_orgba1j)
      cblk(1,VORGBA1I ) =   chem(i,k,j,p_orgba1i)
      cblk(1,VORGBA2J ) =   chem(i,k,j,p_orgba2j)
      cblk(1,VORGBA2I ) =   chem(i,k,j,p_orgba2i)
      cblk(1,VORGBA3J ) =   chem(i,k,j,p_orgba3j)
      cblk(1,VORGBA3I ) =   chem(i,k,j,p_orgba3i)
      cblk(1,VORGBA4J ) =   chem(i,k,j,p_orgba4j)
      cblk(1,VORGBA4I ) =   chem(i,k,j,p_orgba4i)
      cblk(1,VORGPAJ  ) =   chem(i,k,j,p_orgpaj)
      cblk(1,VORGPAI  ) =   chem(i,k,j,p_orgpai)
      cblk(1,VECJ     ) =   chem(i,k,j,p_ecj)
      cblk(1,VECI     ) =   chem(i,k,j,p_eci)
      cblk(1,VP25AJ   ) =   chem(i,k,j,p_p25j)
      cblk(1,VP25AI   ) =   chem(i,k,j,p_p25i)
      cblk(1,VANTHA   ) =   chem(i,k,j,p_antha)
      cblk(1,VSEAS    ) =   chem(i,k,j,p_seas)
      cblk(1,VSOILA   ) =   chem(i,k,j,p_soila)
      cblk(1,VH2OAJ   ) =   max(epsilc,h2oaj(i,k,j))
      cblk(1,VH2OAI   ) =   max(epsilc,h2oai(i,k,j))
      cblk(1,VNU3     ) =   max(epsilc,nu3(i,k,j))
      cblk(1,VAC3     ) =   max(epsilc,ac3(i,k,j))
      cblk(1,VCOR3    ) =   max(epsilc,cor3(i,k,j))
      cblk(1,VCVARO1  ) =   max(epsilc,cvaro1(i,k,j))
      cblk(1,VCVARO2  ) =   max(epsilc,cvaro2(i,k,j))
      cblk(1,VCVALK1  ) =   max(epsilc,cvalk1(i,k,j))
      cblk(1,VCVOLE1  ) =   max(epsilc,cvole1(i,k,j))




      cblk(1,VCVAPI1  ) =   max(epsilc,cvapi1(i,k,j))
      cblk(1,VCVAPI2  ) =   max(epsilc,cvapi2(i,k,j))
      cblk(1,VCVLIM1  ) =   max(epsilc,cvlim1(i,k,j))
      cblk(1,VCVLIM2  ) =   max(epsilc,cvlim2(i,k,j))




         epmcoarse(1) = 0.
         epm25i(1)    = 0.
         epm25j (1)   = 0.
         eeci_in      = 0.
         eecj_in      = 0.
         eorgi_in     = 0.
         eorgj_in     = 0.
         cblk(1,VSO4AJ   ) = chem(i,k,j,p_so4aj)
         cblk(1,VSO4AI   ) = chem(i,k,j,p_so4ai)
         cblk(1,VNO3AJ   ) = chem(i,k,j,p_no3aj)
         cblk(1,VNO3AI   ) = chem(i,k,j,p_no3ai)
         cblk(1,VNAAJ   )  = chem(i,k,j,p_naaj)
         cblk(1,VNAAI   )  = chem(i,k,j,p_naai)
         cblk(1,VCLAJ   )  = chem(i,k,j,p_claj)
         cblk(1,VCLAI   )  = chem(i,k,j,p_clai)































      cblk(1,vsulf) = vsulf_in
      cblk(1,vhno3) = nitrate_in
      cblk(1,vnh3) = nh3_in
      cblk(1,vhcl) = hcl_in
      cblk(1,VNH4AJ   ) =   chem(i,k,j,p_nh4aj)
      cblk(1,VNH4AI   ) =   chem(i,k,j,p_nh4ai)
      cblk(1,VNU0     ) =   max(1.e7,chem(i,k,j,p_nu0))
      cblk(1,VAC0     ) =   max(1.e7,chem(i,k,j,p_ac0))
      cblk(1,VCORN    ) =   chem(i,k,j,p_corn)


      if(debug_level.ge.1)then
     if(i.eq.its.and.j.eq.jts.and.k.eq.kts)then
        print*,'in a_mechanisms',i,j,k
        print*,'NSPCSDA, BLKSIZE',NSPCSDA, BLKSIZE
        print*,'k,DTA,PRES,TEMP,RELHUM',k,DTstep,10.*P(k),T(k),RH(k)
        print*,'nitrate_in, nh3_in, vsulf_in, so4rat_in', &
                nitrate_in, nh3_in, vsulf_in, so4rat_in
        print*,'drog_in,ldrog',drog_in,ldrog
        print*,'condvap_in,NCV,NACV',condvap_in,NCV,NACV
        print*,'eeci_in, eecj_in, eorgi_in, eorgj_in,convfac' &
            ,eeci_in, eecj_in, eorgi_in, eorgj_in,convfac
        print*,'CBLK',CBLK
       endif
    end if
      CALL rpmmod3(nspcsda,blksize,k,dtstep,10.*p(k),t(k),rh(k),nitrate_in,nh3_in, &
        vsulf_in,so4rat_in,drog_in,ldrog,condvap_in,ncv,nacv,eeci_in,eecj_in, &
        eorgi_in,eorgj_in,epm25i,epm25j,epmcoarse,soilrat_in,cblk,i,j,k)
      chem(i,k,j,p_so4aj) = cblk(1,VSO4AJ   )
      chem(i,k,j,p_so4ai) = cblk(1,VSO4AI   )
      chem(i,k,j,p_nh4aj) = cblk(1,VNH4AJ   )
      chem(i,k,j,p_nh4ai) = cblk(1,VNH4AI   )
      chem(i,k,j,p_no3aj) = cblk(1,VNO3AJ   )
      chem(i,k,j,p_no3ai) = cblk(1,VNO3AI   )
      chem(i,k,j,p_naaj) = cblk(1,VNAAJ   )
      chem(i,k,j,p_naai) = cblk(1,VNAAI   )

      chem(i,k,j,p_claj) = cblk(1,VCLAJ   )
      chem(i,k,j,p_clai) = cblk(1,VCLAI   )

      chem(i,k,j,p_orgaro1j) = cblk(1,VORGARO1J)
      chem(i,k,j,p_orgaro1i) = cblk(1,VORGARO1I)
      chem(i,k,j,p_orgaro2j) = cblk(1,VORGARO2J)
      chem(i,k,j,p_orgaro2i) = cblk(1,VORGARO2I)
      chem(i,k,j,p_orgalk1j) = cblk(1,VORGALK1J)
      chem(i,k,j,p_orgalk1i) = cblk(1,VORGALK1I)
      chem(i,k,j,p_orgole1j) = cblk(1,VORGOLE1J)
      chem(i,k,j,p_orgole1i) = cblk(1,VORGOLE1I)
      chem(i,k,j,p_orgba1j) = cblk(1,VORGBA1J )
      chem(i,k,j,p_orgba1i) = cblk(1,VORGBA1I )
      chem(i,k,j,p_orgba2j) = cblk(1,VORGBA2J )
      chem(i,k,j,p_orgba2i) = cblk(1,VORGBA2I )
      chem(i,k,j,p_orgba3j) = cblk(1,VORGBA3J )
      chem(i,k,j,p_orgba3i) = cblk(1,VORGBA3I )
      chem(i,k,j,p_orgba4j) = cblk(1,VORGBA4J )
      chem(i,k,j,p_orgba4i) = cblk(1,VORGBA4I )
      chem(i,k,j,p_orgpaj) = cblk(1,VORGPAJ  )
      chem(i,k,j,p_orgpai) = cblk(1,VORGPAI  )
      chem(i,k,j,p_ecj) = cblk(1,VECJ     )
      chem(i,k,j,p_eci) = cblk(1,VECI     )
      chem(i,k,j,p_p25j) = cblk(1,VP25AJ   )
      chem(i,k,j,p_p25i) = cblk(1,VP25AI   )
      chem(i,k,j,p_antha) =cblk(1,VANTHA   )
      chem(i,k,j,p_seas) = cblk(1,VSEAS    )
      chem(i,k,j,p_soila) = cblk(1,VSOILA   )
      chem(i,k,j,p_nu0) = max(1.e7,cblk(1,VNU0     ))
      chem(i,k,j,p_ac0) = max(1.e7,cblk(1,VAC0     ))

      chem(i,k,j,p_corn) = cblk(1,VCORN    )
      h2oaj(i,k,j) = cblk(1,VH2OAJ   )
      h2oai(i,k,j) = cblk(1,VH2OAI   )
      nu3(i,k,j) = cblk(1,VNU3     )
      ac3(i,k,j) = cblk(1,VAC3     )
      cor3(i,k,j) = cblk(1,VCOR3    )
      cvaro1(i,k,j) = cblk(1,VCVARO1  )
      cvaro2(i,k,j) = cblk(1,VCVARO2  )
      cvalk1(i,k,j) = cblk(1,VCVALK1  )
      cvole1(i,k,j) = cblk(1,VCVOLE1  )




      cvapi1(i,k,j) = cblk(1,VCVAPI1  )
      cvapi2(i,k,j) = cblk(1,VCVAPI2  )
      cvlim1(i,k,j) = cblk(1,VCVLIM1  )
      cvlim2(i,k,j) = cblk(1,VCVLIM2  )

      chem(i,k,j,p_sulf)=max(epsilc,cblk(1,vsulf)/CONVFAC/MWSO4)
      chem(i,k,j,p_hno3)=max(epsilc,cblk(1,vhno3)/CONVFAC/MWHNO3)
      chem(i,k,j,p_nh3)=max(epsilc,cblk(1,vnh3)/CONVFAC/MWNH3)
      chem(i,k,j,p_hcl)=max(epsilc,cblk(1,vhcl)/CONVFAC/MWHCL)








      enddo         
 100  continue      




      do l=p_so4aj,num_chem
         do j=jts,jte
            do k=kts,kte
               do i=its,ite
                  chem(i,k,j,l)=max(epsilc,chem(i,k,j,l)*alt(i,k,j))
               enddo
            enddo
         enddo
      enddo

    END SUBROUTINE sorgam_driver

    SUBROUTINE sum_pm_sorgam (                                         &
         alt, chem, h2oaj, h2oai,                                      &
         pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10,                   &
         dust_opt,ids,ide, jds,jde, kds,kde,                                    &
         ims,ime, jms,jme, kms,kme,                                    &
         its,ite, jts,jte, kts,kte                                     )

   INTEGER,      INTENT(IN   )    ::  dust_opt,                        &
                                      ids,ide, jds,jde, kds,kde,       &
                                      ims,ime, jms,jme, kms,kme,       &
                                      its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(IN ) :: chem

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(IN ) :: alt,h2oaj,h2oai

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(OUT) :: pm2_5_dry,pm2_5_water,pm2_5_dry_ec,pm10

   INTEGER :: i,ii,j,jj,k,n



      pm2_5_dry(its:ite, kts:kte, jts:jte)    = 0.
      pm2_5_water(its:ite, kts:kte, jts:jte)  = 0.
      pm2_5_dry_ec(its:ite, kts:kte, jts:jte) = 0.
      do j=jts,jte
         jj=min(jde-1,j)
      do k=kts,kte
      do i=its,ite
         ii=min(ide-1,i)








         do n=p_so4aj,p_p25i
            pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j)+chem(ii,k,jj,n)
         enddo
         if( p_p25cwi .gt. p_p25i) then
         do n=p_so4cwj,p_p25cwi
            pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j)+chem(ii,k,jj,n)
         enddo
         endif
         pm2_5_dry_ec(i,k,j) = pm2_5_dry_ec(i,k,j)+chem(ii,k,jj,p_ecj) &
                               + chem(ii,k,jj,p_eci)
         pm2_5_water(i,k,j) =  pm2_5_water(i,k,j)+h2oaj(i,k,j)       &
                               + h2oai(i,k,j)

         
         pm2_5_dry(i,k,j)    = pm2_5_dry(i,k,j) / alt(ii,k,jj)
         pm2_5_dry_ec(i,k,j) = pm2_5_dry_ec(i,k,j) / alt(ii,k,jj)
         pm2_5_water(i,k,j)  = pm2_5_water(i,k,j) / alt(ii,k,jj)
      enddo
      enddo
      enddo
      do j=jts,jte
         jj=min(jde-1,j)
         do k=kts,kte
            do i=its,ite
               ii=min(ide-1,i)








               pm10(i,k,j) = pm2_5_dry(i,k,j)                       &
                           + ( chem(ii,k,jj,p_antha)               &
                           + chem(ii,k,jj,p_soila)                 &
                           + chem(ii,k,jj,p_seas) ) / alt(ii,k,jj)
               if( p_p25cwi .gt. p_p25i) then
               pm10(i,k,j) = pm10(i,k,j)                       &
                           + ( chem(ii,k,jj,p_anthcw)               &
                           + chem(ii,k,jj,p_soilcw)                 &
                           + chem(ii,k,jj,p_seascw) ) / alt(ii,k,jj)
               endif

            enddo
         enddo
      enddo
    END SUBROUTINE sum_pm_sorgam

    SUBROUTINE sorgam_depdriver (id,config_flags,ktau,dtstep,           &
               ust,t_phy,moist,p8w,t8w,rmol,znt,pbl,                    &
               alt,p_phy,chem,rho_phy,dz8w,z,z_at_w,                    &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,               &
               aer_res,vgsa,                                            &
               numgas,ddflx, &
               numaer,                                                  &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )

   USE module_configure,only:  grid_config_rec_type
   TYPE (grid_config_rec_type) , INTENT (IN) :: config_flags
   INTEGER,      INTENT(IN   )    ::                             &
                                      numgas, &
                                      numaer,                    &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      its,ite, jts,jte, kts,kte, &
                                      id,ktau

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),         &
         INTENT(INOUT ) ::                                   chem



   REAL, DIMENSION( its:ite, jts:jte, numaer ),                       &
         INTENT(INOUT ) ::                                             &
         vgsa
        real, intent(inout),   &
                dimension( ims:ime, jms:jme, numgas+1:num_chem ) :: &
                ddflx

   REAL, DIMENSION( its:ite, jts:jte ),                       &
         INTENT(INOUT ) ::                                             &
         aer_res

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(INOUT ) ::                                             &
           h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2,    &
           cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                      t_phy,    &
                                                      alt,      &
                                                      p_phy,    &
                                                      dz8w,     &
                                                      z    ,    &
                                       t8w,p8w,z_at_w ,  &
                                                    rho_phy
   REAL,  DIMENSION( ims:ime ,  jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
               ust,rmol, pbl, znt
      REAL,      INTENT(IN   ) ::                               &
                             dtstep
                                                                                               
      REAL rgas
      DATA rgas/8.314510/
      REAL convfac,convfac2


      INTEGER blksize
      PARAMETER (blksize=1)



      INTEGER nspcsda
      PARAMETER (nspcsda=l1ae) 


      INTEGER nacv
      PARAMETER (nacv=lcva) 

      INTEGER ncv
      PARAMETER (ncv=lspcv) 

      REAL cblk(blksize,nspcsda) 
                                   
      REAL soilrat_in
                    
                    
      REAL nitrate_in
                    
      REAL nh3_in

      REAL hcl_in
                    
      REAL vsulf_in

      REAL so4rat_in
                    
                    
      REAL pres
                    
      REAL temp
                    
      REAL relhum
                    
      REAL ::p(kts:kte),t(kts:kte),rh(kts:kte)




      REAL mwso4
      PARAMETER (mwso4=96.0576)


      REAL mwhno3
      PARAMETER (mwhno3=63.01287)


      REAL mwnh3
      PARAMETER (mwnh3=17.03061)


      REAL mwhcl
      PARAMETER (mwhcl=36.46100)






      REAL mwec
      PARAMETER (mwec=12.0)


      REAL mwaro1
      PARAMETER (mwaro1=150.0)


      REAL mwaro2
      PARAMETER (mwaro2=150.0)


      REAL mwalk1
      PARAMETER (mwalk1=140.0)


      REAL mwalk2
      PARAMETER (mwalk2=140.0)



      REAL mwole1
      PARAMETER (mwole1=140.0)


      REAL mwapi1
      PARAMETER (mwapi1=200.0)


      REAL mwapi2
      PARAMETER (mwapi2=200.0)


      REAL mwlim1
      PARAMETER (mwlim1=200.0)


      REAL mwlim2
      PARAMETER (mwlim2=200.0)
      INTEGER NUMCELLS  


       PARAMETER( NUMCELLS = 1)


      REAL RA(BLKSIZE )             
      REAL USTAR( BLKSIZE )         
      REAL WSTAR( BLKSIZE )         
      REAL PBLH( BLKSIZE )          
      REAL ZNTT( BLKSIZE )          
      REAL RMOLM( BLKSIZE )         

      REAL BLKPRS(BLKSIZE)         
      REAL BLKTA(BLKSIZE)          
      REAL BLKDENS(BLKSIZE)        





      
      REAL XLM( BLKSIZE )           
      REAL AMU( BLKSIZE )           
      

      REAL VSED( BLKSIZE, NASPCSSED) 
      REAL VDEP( BLKSIZE, NASPCSDEP) 


      
      REAL DGNUC( BLKSIZE )         
      REAL DGACC( BLKSIZE )         
      REAL DGCOR( BLKSIZE )         
      




      
      REAL PMASSN( BLKSIZE )        
      REAL PMASSA( BLKSIZE )        
      REAL PMASSC( BLKSIZE )        



      REAL PDENSN( BLKSIZE )        
      REAL PDENSA( BLKSIZE )        
      REAL PDENSC( BLKSIZE )        



      REAL KNNUC ( BLKSIZE )        
      REAL KNACC ( BLKSIZE )        
      REAL KNCOR ( BLKSIZE )        




   INTEGER :: i,j,k,l


      do l=1,numaer
      do i=its,ite
      do j=jts,jte
      vgsa(i,j,l)=0.
      enddo
      enddo
      enddo
      vdep=0.
      do 100 j=jts,jte
         do 100 i=its,ite
            cblk=epsilc
            do k=kts,kte
               t(k) = t_phy(i,k,j)
               p(k) = .001*p_phy(i,k,j)
               rh(k) = MIN( 100.,100. * moist(i,k,j,p_qv) / &
               (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
               (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
               rh(k)=max(.05,0.01*rh(k))
            enddo

            k=kts
               convfac = p(k)/rgas/t(k)*1000.
               nitrate_in =chem(i,k,j,p_hno3)*convfac*mwhno3
               nh3_in = chem(i,k,j,p_nh3)*convfac*mwnh3
               vsulf_in = chem(i,k,j,p_sulf)*convfac*mwso4   
               hcl_in = chem(i,k,j,p_hcl)*convfac*mwhcl              
 

      BLKPRS(BLKSIZE)   = 1.e3*P(K)                
      BLKTA(BLKSIZE)   = T(K)         
      USTAR(BLKSIZE) = max(1.e-1,UST(i,j))
      WSTAR(BLKSIZE) = 0.
      pblh(blksize) = pbl(i,j)
      zntt(blksize) = znt(i,j)
      rmolm(blksize)= rmol(i,j)
      convfac2=1./alt(i,k,j)
      BLKDENS(BLKSIZE)=convfac2
      cblk(1,vsulf) = max(epsilc,vsulf_in)
      cblk(1,vhno3) = max(epsilc,nitrate_in)
      cblk(1,vnh3) = max(epsilc,nh3_in)
      cblk(1,vhcl) = max(epsilc,hcl_in)
      cblk(1,VSO4AJ   ) =   max(epsilc,chem(i,k,j,p_so4aj)*convfac2)
      cblk(1,VSO4AI   ) =   max(epsilc,chem(i,k,j,p_so4ai)*convfac2)
      cblk(1,VNH4AJ   ) =   max(epsilc,chem(i,k,j,p_nh4aj)*convfac2)
      cblk(1,VNH4AI   ) =   max(epsilc,chem(i,k,j,p_nh4ai)*convfac2)
      cblk(1,VNO3AJ   ) =   max(epsilc,chem(i,k,j,p_no3aj)*convfac2)
      cblk(1,VNO3AI   ) =   max(epsilc,chem(i,k,j,p_no3ai)*convfac2)
      if (p_naai >= param_first_scalar) &
         cblk(1,VNAAI ) =   max(epsilc,chem(i,k,j,p_naai)*convfac2)
      if (p_naaj >= param_first_scalar) &
         cblk(1,VNAAJ ) =   max(epsilc,chem(i,k,j,p_naaj)*convfac2)
      if (p_clai >= param_first_scalar) &
         cblk(1,VCLAI ) =   max(epsilc,chem(i,k,j,p_clai)*convfac2)
      if (p_claj >= param_first_scalar) &
         cblk(1,VCLAJ ) =   max(epsilc,chem(i,k,j,p_claj)*convfac2)
      cblk(1,VORGARO1J) =   max(epsilc,chem(i,k,j,p_orgaro1j)*convfac2)
      cblk(1,VORGARO1I) =   max(epsilc,chem(i,k,j,p_orgaro1i)*convfac2)
      cblk(1,VORGARO2J) =   max(epsilc,chem(i,k,j,p_orgaro2j)*convfac2)
      cblk(1,VORGARO2I) =   max(epsilc,chem(i,k,j,p_orgaro2i)*convfac2)
      cblk(1,VORGALK1J) =   max(epsilc,chem(i,k,j,p_orgalk1j)*convfac2)
      cblk(1,VORGALK1I) =   max(epsilc,chem(i,k,j,p_orgalk1i)*convfac2)
      cblk(1,VORGOLE1J) =   max(epsilc,chem(i,k,j,p_orgole1j)*convfac2)
      cblk(1,VORGOLE1I) =   max(epsilc,chem(i,k,j,p_orgole1i)*convfac2)
      cblk(1,VORGBA1J ) =   max(epsilc,chem(i,k,j,p_orgba1j)*convfac2)
      cblk(1,VORGBA1I ) =   max(epsilc,chem(i,k,j,p_orgba1i)*convfac2)
      cblk(1,VORGBA2J ) =   max(epsilc,chem(i,k,j,p_orgba2j)*convfac2)
      cblk(1,VORGBA2I ) =   max(epsilc,chem(i,k,j,p_orgba2i)*convfac2)
      cblk(1,VORGBA3J ) =   max(epsilc,chem(i,k,j,p_orgba3j)*convfac2)
      cblk(1,VORGBA3I ) =   max(epsilc,chem(i,k,j,p_orgba3i)*convfac2)
      cblk(1,VORGBA4J ) =   max(epsilc,chem(i,k,j,p_orgba4j)*convfac2)
      cblk(1,VORGBA4I ) =   max(epsilc,chem(i,k,j,p_orgba4i)*convfac2)
      cblk(1,VORGPAJ  ) =   max(epsilc,chem(i,k,j,p_orgpaj)*convfac2)
      cblk(1,VORGPAI  ) =   max(epsilc,chem(i,k,j,p_orgpai)*convfac2)
      cblk(1,VECJ     ) =   max(epsilc,chem(i,k,j,p_ecj)*convfac2)
      cblk(1,VECI     ) =   max(epsilc,chem(i,k,j,p_eci)*convfac2)
      cblk(1,VP25AJ   ) =   max(epsilc,chem(i,k,j,p_p25j)*convfac2)
      cblk(1,VP25AI   ) =   max(epsilc,chem(i,k,j,p_p25i)*convfac2)
      cblk(1,VANTHA   ) =   max(epsilc,chem(i,k,j,p_antha)*convfac2)
      cblk(1,VSEAS    ) =   max(epsilc,chem(i,k,j,p_seas)*convfac2)
      cblk(1,VSOILA   ) =   max(epsilc,chem(i,k,j,p_soila)*convfac2)
      cblk(1,VNU0     ) =   max(epsilc,chem(i,k,j,p_nu0)*convfac2)
      cblk(1,VAC0     ) =   max(epsilc,chem(i,k,j,p_ac0)*convfac2)
      cblk(1,VCORN    ) =   max(epsilc,chem(i,k,j,p_corn)*convfac2)
      cblk(1,VH2OAJ   ) =   h2oaj(i,k,j)
      cblk(1,VH2OAI   ) =   h2oai(i,k,j)
      cblk(1,VNU3     ) =   nu3(i,k,j)
      cblk(1,VAC3     ) =   ac3(i,k,j)
      cblk(1,VCOR3    ) =   cor3(i,k,j)
      cblk(1,VCVARO1  ) =   cvaro1(i,k,j)
      cblk(1,VCVARO2  ) =   cvaro2(i,k,j)
      cblk(1,VCVALK1  ) =   cvalk1(i,k,j)
      cblk(1,VCVOLE1  ) =   cvole1(i,k,j)




      cblk(1,VCVAPI1  ) =   cvapi1(i,k,j)
      cblk(1,VCVAPI2  ) =   cvapi2(i,k,j)
      cblk(1,VCVLIM1  ) =   cvlim1(i,k,j)
      cblk(1,VCVLIM2  ) =   cvlim2(i,k,j)














        CALL MODPAR(  BLKSIZE, NSPCSDA, NUMCELLS,             &
             CBLK,                                            &
             BLKTA, BLKPRS,                                   &
             PMASSN, PMASSA, PMASSC,                          &
             PDENSN, PDENSA, PDENSC,                          &
             XLM, AMU,                                        &
             DGNUC, DGACC, DGCOR,                             &
             KNNUC, KNACC,KNCOR    )                                   

        if (config_flags%aer_drydep_opt == 11) then
        CALL VDVG(  BLKSIZE, NSPCSDA, NUMCELLS,k,CBLK, &
             BLKTA, BLKDENS, AER_RES(I,j), USTAR, WSTAR,  AMU,   &   
             DGNUC, DGACC, DGCOR,                      &
             KNNUC, KNACC,KNCOR,                       &
             PDENSN, PDENSA, PDENSC,                   &
             VSED, VDEP )                                             
        else
        CALL VDVG_2(  BLKSIZE, NSPCSDA, NUMCELLS,k,CBLK, &
             BLKTA, BLKDENS, AER_RES(I,j), USTAR, PBLH,&
             ZNTT, RMOLM, AMU, DGNUC, DGACC, DGCOR,XLM,&
             KNNUC, KNACC,KNCOR,                       &
             PDENSN, PDENSA, PDENSC,                   &
             VSED, VDEP )
        endif

        VGSA(i, j, VSO4AJ   )  =  VDEP(1, VDMACC )
        VGSA(i, j, VSO4AI   )  =  VDEP(1, VDMNUC )
        VGSA(i, j, VNH4AJ   )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VNH4AI   )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VNO3AJ   )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VNO3AI   )  =  VGSA(i, j, VSO4AI )
        if (p_naai >= param_first_scalar) VGSA(i, j, VNAAI )  =  VGSA(i, j, VSO4AI )
        if (p_naaj >= param_first_scalar) VGSA(i, j, VNAAJ )  =  VGSA(i, j, VSO4AJ )
        if (p_clai >= param_first_scalar) VGSA(i, j, VCLAI )  =  VGSA(i, j, VSO4AI )
        if (p_claj >= param_first_scalar) VGSA(i, j, VCLAJ )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGARO1J)  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGARO1I)  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VORGARO2J)  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGARO2I)  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VORGALK1J)  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGALK1I)  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VORGOLE1J)  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGOLE1I)  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VORGBA1J )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGBA1I )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VORGBA2J )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGBA2I )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VORGBA3J )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGBA3I )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VORGBA4J )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGBA4I )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VORGPAJ  )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VORGPAI  )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VECJ     )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VECI     )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VP25AJ   )  =  VGSA(i, j, VSO4AJ )
        VGSA(i, j, VP25AI   )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, VANTHA   )  =  VDEP(1, VDMCOR )
        VGSA(i, j, VSEAS    )  =  VGSA(i, j, VANTHA )
        VGSA(i, j, VSOILA   )  =  VGSA(i, j, VANTHA )
        VGSA(i, j, VNU0     )  =  VDEP(1, VDNNUC )
        VGSA(i, j, VAC0     )  =  VDEP(1, VDNACC )
        VGSA(i, j, VCORN    )  =  VDEP(1, VDNCOR )
        if( config_flags%diagnostic_dep == 1) then
        ddflx(i,j,p_so4aj)=ddflx(i,j,p_so4aj)+chem(i,k,j,p_so4aj)/alt(i,k,j)*VGSA(i,j,VSO4AJ)*dtstep
        ddflx(i,j,p_so4ai)=ddflx(i,j,p_so4ai)+chem(i,k,j,p_so4ai)/alt(i,k,j)*VGSA(i,j,VSO4AI)*dtstep
        ddflx(i,j,p_nh4aj)=ddflx(i,j,p_nh4aj)+chem(i,k,j,p_nh4aj)/alt(i,k,j)*VGSA(i,j,VNH4AJ)*dtstep
        ddflx(i,j,p_nh4ai)=ddflx(i,j,p_nh4ai)+chem(i,k,j,p_nh4ai)/alt(i,k,j)*VGSA(i,j,VNH4Ai)*dtstep
        ddflx(i,j,p_no3aj)=ddflx(i,j,p_no3aj)+chem(i,k,j,p_no3aj)/alt(i,k,j)*VGSA(i,j,VNO3AJ)*dtstep
        ddflx(i,j,p_no3ai)=ddflx(i,j,p_no3ai)+chem(i,k,j,p_no3ai)/alt(i,k,j)*VGSA(i,j,VNO3AI)*dtstep
        ddflx(i,j,p_orgaro1j)=ddflx(i,j,p_orgaro1j)+chem(i,k,j,p_orgaro1j)/alt(i,k,j)*VGSA(i,j,VORGARO1J)*dtstep
        ddflx(i,j,p_orgaro1i)=ddflx(i,j,p_orgaro1i)+chem(i,k,j,p_orgaro1i)/alt(i,k,j)*VGSA(i,j,VORGARO1I)*dtstep
        ddflx(i,j,p_orgaro2j)=ddflx(i,j,p_orgaro2j)+chem(i,k,j,p_orgaro2j)/alt(i,k,j)*VGSA(i,j,VORGARO2J)*dtstep
        ddflx(i,j,p_orgaro2i)=ddflx(i,j,p_orgaro2i)+chem(i,k,j,p_orgaro2i)/alt(i,k,j)*VGSA(i,j,VORGARO2I)*dtstep
        ddflx(i,j,p_orgalk1j)=ddflx(i,j,p_orgalk1j)+chem(i,k,j,p_orgalk1j)/alt(i,k,j)*VGSA(i,j,VORGALK1J)*dtstep
        ddflx(i,j,p_orgalk1i)=ddflx(i,j,p_orgalk1i)+chem(i,k,j,p_orgalk1i)/alt(i,k,j)*VGSA(i,j,VORGALK1I)*dtstep
        ddflx(i,j,p_orgole1j)=ddflx(i,j,p_orgole1j)+chem(i,k,j,p_orgole1j)/alt(i,k,j)*VGSA(i,j,VORGOLE1J)*dtstep
        ddflx(i,j,p_orgole1i)=ddflx(i,j,p_orgole1i)+chem(i,k,j,p_orgole1i)/alt(i,k,j)*VGSA(i,j,VORGOLE1I)*dtstep
        ddflx(i,j,p_orgba1j)=ddflx(i,j,p_orgba1j)+chem(i,k,j,p_orgba1j)/alt(i,k,j)*VGSA(i,j,VORGBA1J)*dtstep
        ddflx(i,j,p_orgba1i)=ddflx(i,j,p_orgba1i)+chem(i,k,j,p_orgba1i)/alt(i,k,j)*VGSA(i,j,VORGBA1I)*dtstep
        ddflx(i,j,p_orgba2j)=ddflx(i,j,p_orgba2j)+chem(i,k,j,p_orgba2j)/alt(i,k,j)*VGSA(i,j,VORGBA2J)*dtstep
        ddflx(i,j,p_orgba2i)=ddflx(i,j,p_orgba2i)+chem(i,k,j,p_orgba2i)/alt(i,k,j)*VGSA(i,j,VORGBA2I)*dtstep
        ddflx(i,j,p_orgba3j)=ddflx(i,j,p_orgba3j)+chem(i,k,j,p_orgba3j)/alt(i,k,j)*VGSA(i,j,VORGBA3J)*dtstep
        ddflx(i,j,p_orgba3i)=ddflx(i,j,p_orgba3i)+chem(i,k,j,p_orgba3i)/alt(i,k,j)*VGSA(i,j,VORGBA3I)*dtstep
        ddflx(i,j,p_orgba4j)=ddflx(i,j,p_orgba4j)+chem(i,k,j,p_orgba4j)/alt(i,k,j)*VGSA(i,j,VORGBA4J)*dtstep
        ddflx(i,j,p_orgba4i)=ddflx(i,j,p_orgba4i)+chem(i,k,j,p_orgba4i)/alt(i,k,j)*VGSA(i,j,VORGBA4I)*dtstep
        ddflx(i,j,p_orgpaj)=ddflx(i,j,p_orgpaj)+chem(i,k,j,p_orgpaj)/alt(i,k,j)*VGSA(i,j,VORGPAJ)*dtstep
        ddflx(i,j,p_orgpai)=ddflx(i,j,p_orgpai)+chem(i,k,j,p_orgpai)/alt(i,k,j)*VGSA(i,j,VORGPAI)*dtstep
        ddflx(i,j,p_ecj)=ddflx(i,j,p_ecj)+chem(i,k,j,p_ecj)/alt(i,k,j)*VGSA(i,j,VECJ)*dtstep
        ddflx(i,j,p_eci)=ddflx(i,j,p_eci)+chem(i,k,j,p_eci)/alt(i,k,j)*VGSA(i,j,VECI)*dtstep
        ddflx(i,j,p_p25j)=ddflx(i,j,p_p25j)+chem(i,k,j,p_p25j)/alt(i,k,j)*VGSA(i,j,VP25AJ)*dtstep
        ddflx(i,j,p_p25i)=ddflx(i,j,p_p25i)+chem(i,k,j,p_p25i)/alt(i,k,j)*VGSA(i,j,VP25AI)*dtstep
        ddflx(i,j,p_naaj)=ddflx(i,j,p_naaj)+chem(i,k,j,p_naaj)/alt(i,k,j)*VGSA(i,j,VNAAJ)*dtstep
        ddflx(i,j,p_naai)=ddflx(i,j,p_naai)+chem(i,k,j,p_naai)/alt(i,k,j)*VGSA(i,j,VNAAI)*dtstep
        ddflx(i,j,p_claj)=ddflx(i,j,p_claj)+chem(i,k,j,p_claj)/alt(i,k,j)*VGSA(i,j,VCLAJ)*dtstep
        ddflx(i,j,p_clai)=ddflx(i,j,p_clai)+chem(i,k,j,p_clai)/alt(i,k,j)*VGSA(i,j,VCLAI)*dtstep
        ddflx(i,j,p_antha)=ddflx(i,j,p_antha)+chem(i,k,j,p_antha)/alt(i,k,j)*VGSA(i,j,VANTHA)*dtstep
        ddflx(i,j,p_seas)=ddflx(i,j,p_seas)+chem(i,k,j,p_seas)/alt(i,k,j)*VGSA(i,j,VSEAS)*dtstep
        ddflx(i,j,p_soila)=ddflx(i,j,p_soila)+chem(i,k,j,p_soila)/alt(i,k,j)*VGSA(i,j,VSOILA)*dtstep
        ddflx(i,j,p_nu0)=ddflx(i,j,p_nu0)+chem(i,k,j,p_nu0)/alt(i,k,j)*VGSA(i,j,VNU0)*dtstep
        ddflx(i,j,p_ac0)=ddflx(i,j,p_ac0)+chem(i,k,j,p_ac0)/alt(i,k,j)*VGSA(i,j,VAC0)*dtstep
        ddflx(i,j,p_corn)=ddflx(i,j,p_corn)+chem(i,k,j,p_corn)/alt(i,k,j)*VGSA(i,j,VCORN)*dtstep
        end if


 100  continue      
                                                                     
END SUBROUTINE sorgam_depdriver

    SUBROUTINE actcof(cat,an,gama,molnu,phimult)



































































      INTEGER xstat0
      PARAMETER (xstat0=0)

      INTEGER xstat1
      PARAMETER (xstat1=1)

      INTEGER xstat2
      PARAMETER (xstat2=2)

      INTEGER xstat3
      PARAMETER (xstat3=3)

      CHARACTER*120 xmsg




      INTEGER ncat
      PARAMETER (ncat=2)


      INTEGER nan
      PARAMETER (nan=3)




      REAL molnu

      REAL phimult
      REAL cat(ncat) 
      REAL an(nan) 
      REAL gama(ncat,nan) 



      CHARACTER*16 & 
        pname
      SAVE pname


      INTEGER ian

      INTEGER icat


      REAL fgama

      REAL i

      REAL r

      REAL s

      REAL ta

      REAL tb

      REAL tc

      REAL texpv

      REAL trm

      REAL twoi

      REAL twosri

      REAL zbar

      REAL zbar2

      REAL zot1

      REAL sri
      REAL f2(ncat)
      REAL f1(nan)
      REAL zp(ncat) 
      REAL zm(nan) 
      REAL bgama(ncat,nan)
      REAL x(ncat,nan)
      REAL m(ncat,nan) 
      REAL lgama0(ncat,nan) 
      REAL y(nan,ncat)
      REAL beta0(ncat,nan) 
      REAL beta1(ncat,nan) 
      REAL cgama(ncat,nan) 
      REAL v1(ncat,nan) 
      REAL v2(ncat,nan) 

      DATA zp/1.0, 1.0/
      DATA zm/2.0, 1.0, 1.0/
      DATA xmsg/' '/
      DATA pname/'ACTCOF'/










      DATA beta0(1,1)/2.98E-2/, beta1(1,1)/0.0/, cgama(1,1)/4.38E-2 & 
        /

      DATA beta0(1,2)/1.2556E-1/, beta1(1,2)/2.8778E-1/, cgama(1,2)/ -5.59E-3 & 
        /

      DATA beta0(1,3)/2.0651E-1/, beta1(1,3)/5.556E-1/, cgama(1,3)/0.0 & 
        /

      DATA beta0(2,1)/4.6465E-2/, beta1(2,1)/ -0.54196/, &
        cgama(2,1)/ -1.2683E-3 & 
        /

      DATA beta0(2,2)/ -7.26224E-3/, beta1(2,2)/ -1.168858/, &
        cgama(2,2)/3.51217E-5 & 
        /

      DATA beta0(2,3)/4.494E-2/, beta1(2,3)/2.3594E-1/, cgama(2,3)/ -2.962E-3 & 
        /

      DATA v1(1,1), v2(1,1)/2.0, 1.0 & 
        /
      DATA v1(2,1), v2(2,1)/2.0, 1.0 & 
        /
      DATA v1(1,2), v2(1,2)/1.0, 1.0 & 
        /
      DATA v1(2,2), v2(2,2)/1.0, 1.0 & 
        /
      DATA v1(1,3), v2(1,3)/1.0, 1.0 & 
        /
      DATA v1(2,3), v2(2,3)/1.0, 1.0 & 
        /






      i = 0.0

      DO icat = 1, ncat
        i = i + cat(icat)*zp(icat)*zp(icat)
      END DO

      DO ian = 1, nan
        i = i + an(ian)*zm(ian)*zm(ian)
      END DO

      i = 0.5*i



      IF (i==0.0) THEN

        DO ian = 1, nan
          DO icat = 1, ncat
            gama(icat,ian) = 0.0
          END DO
        END DO



        RETURN

      ELSE IF (i<0.0) THEN
        xmsg = 'WARNING: Ionic strength below zero (= negative ion concentrations) - setting ion concentrations to zero.'
        call wrf_message(xmsg)
        DO ian = 1, nan
          DO icat = 1, ncat
            gama(icat,ian) = 0.0
          END DO
        END DO
        RETURN
		
      END IF



      sri = sqrt(i)
      twosri = 2.0*sri
      twoi = 2.0*i
      texpv = 1.0 - exp(-twosri)*(1.0+twosri-twoi)
      r = 1.0 + 0.75*i
      s = 1.0 + 1.5*i
      zot1 = 0.511*sri/(1.0+sri)



      fgama = -0.392*((sri/(1.0+1.2*sri)+(2.0/1.2)*alog(1.0+1.2*sri)))

      DO icat = 1, ncat
        DO ian = 1, nan

          bgama(icat,ian) = 2.0*beta0(icat,ian) + (2.0*beta1(icat,ian)/(4.0*i) &
            )*texpv



          m(icat,ian) = (cat(icat)**v1(icat,ian)*an(ian)**v2(icat,ian))** &
            (1.0/(v1(icat,ian)+v2(icat,ian)))



          lgama0(icat,ian) = (zp(icat)*zm(ian)*fgama+m(icat,ian)*(2.0*v1(icat, &
            ian)*v2(icat,ian)/(v1(icat,ian)+v2(icat,ian))*bgama(icat, &
            ian))+m(icat,ian)*m(icat,ian)*(2.0*(v1(icat,ian)* &
            v2(icat,ian))**1.5/(v1(icat,ian)+v2(icat,ian))*cgama(icat, &
            ian)))/2.302585093

        END DO
      END DO



      DO ian = 1, nan
        DO icat = 1, ncat
          zbar = (zp(icat)+zm(ian))*0.5
          zbar2 = zbar*zbar
          y(ian,icat) = zbar2*an(ian)/i
          x(icat,ian) = zbar2*cat(icat)/i
        END DO
      END DO

      DO ian = 1, nan
        f1(ian) = 0.0
        DO icat = 1, ncat
          f1(ian) = f1(ian) + x(icat,ian)*lgama0(icat,ian) + &
            zot1*zp(icat)*zm(ian)*x(icat,ian)
        END DO
      END DO

      DO icat = 1, ncat
        f2(icat) = 0.0
        DO ian = 1, nan
          f2(icat) = f2(icat) + y(ian,icat)*lgama0(icat,ian) + &
            zot1*zp(icat)*zm(ian)*y(ian,icat)
        END DO
      END DO



      DO ian = 1, nan
        DO icat = 1, ncat

          ta = -zot1*zp(icat)*zm(ian)
          tb = zp(icat)*zm(ian)/(zp(icat)+zm(ian))
          tc = (f2(icat)/zp(icat)+f1(ian)/zm(ian))
          trm = ta + tb*tc

          IF (trm>30.0) THEN
            gama(icat,ian) = 1.0E+30


          ELSE
            gama(icat,ian) = 10.0**trm
          END IF

        END DO
      END DO

      RETURN

    END SUBROUTINE actcof




















    SUBROUTINE aeroproc(blksize,nspcsda,numcells,layer,cblk,dt,blkta,blkprs, &
        blkdens,blkrh,so4rat,orgaro1rat,orgaro2rat,orgalk1rat,orgole1rat, &
        orgbio1rat,orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog,ncv,nacv,epm25i, &
        epm25j,eorgi,eorgj,eeci,eecj,epmcoarse,esoil,eseas,xlm,amu,dgnuc, &
        dgacc,dgcor,pmassn,pmassa,pmassc,pdensn,pdensa,pdensc,knnuc,knacc, &
        kncor,fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,cgrn3,cgra3,urn00, &
        ura00,brna01,c30,deltaso4a,igrid,jgrid,kgrid)






      INTEGER blksize

      INTEGER nspcsda

      INTEGER numcells

      INTEGER layer

      INTEGER ldrog
      REAL cblk(blksize,nspcsda) 

      REAL dt



      REAL blkta(blksize) 
      REAL blkprs(blksize) 
      REAL blkdens(blksize) 
      REAL blkrh(blksize) 



      REAL so4rat(blksize) 



      INTEGER ncv

      INTEGER nacv


      REAL drog(blksize,ldrog) 


      REAL orgaro1rat(blksize)


      REAL orgaro2rat(blksize)


      REAL orgalk1rat(blksize)


      REAL orgole1rat(blksize)


      REAL orgbio1rat(blksize)


      REAL orgbio2rat(blksize)


      REAL orgbio3rat(blksize)


      REAL orgbio4rat(blksize)




      REAL epm25i(blksize) 
      REAL epm25j(blksize) 


      REAL eorgi(blksize) 
      REAL eorgj(blksize) 


      REAL eeci(blksize) 
      REAL eecj(blksize) 


      REAL esoil(blksize) 
      REAL eseas(blksize) 
      REAL epmcoarse(blksize) 






      REAL xlm(blksize) 
      REAL amu(blksize) 



      REAL dgnuc(blksize) 
      REAL dgacc(blksize) 
      REAL dgcor(blksize) 






      REAL pmassn(blksize) 
      REAL pmassa(blksize) 
      REAL pmassc(blksize) 



      REAL pdensn(blksize) 
      REAL pdensa(blksize) 
      REAL pdensc(blksize) 



      REAL knnuc(blksize) 
      REAL knacc(blksize) 
      REAL kncor(blksize) 



      REAL fconcn(blksize)
      REAL fconca(blksize)

      REAL fconcn_org(blksize)
      REAL fconca_org(blksize)





      REAL dmdt(blksize) 




      REAL dndt(blksize) 





      REAL cgrn3(blksize) 
      REAL cgra3(blksize) 





      REAL urn00(blksize) 
      REAL ura00(blksize) 




      REAL brna01(blksize) 


      REAL c30(blksize)                                                        




      REAL deltaso4a(blksize) 











      CHARACTER*16 pname
      PARAMETER (pname=' AEROPROC       ')

      INTEGER unit
      PARAMETER (unit=20)
      integer igrid,jgrid,kgrid,isorop

      isorop=1







        if(blkta(1).ge.233.15.and.blkrh(1).ge.0.1 .and. isorop.eq.1)then
           CALL eql3(blksize,nspcsda,numcells,cblk,blkta,blkrh,igrid,jgrid,kgrid)
        else if (blkta(1).ge.233.15.and.blkrh(1).ge.0.1 .and. isorop.eq.0)then
           CALL eql4(blksize,nspcsda,numcells,cblk,blkta,blkrh)
        endif



      CALL modpar(blksize,nspcsda,numcells,cblk,blkta,blkprs,pmassn,pmassa, &
        pmassc,pdensn,pdensa,pdensc,xlm,amu,dgnuc,dgacc,dgcor,knnuc,knacc, &
        kncor)



      CALL coagrate(blksize,nspcsda,numcells,cblk,blkta,pdensn,pdensa,amu, &
        dgnuc,dgacc,knnuc,knacc,urn00,ura00,brna01,c30)




      CALL nuclcond(blksize,nspcsda,numcells,cblk,dt,layer,blkta,blkprs,blkrh, &
        so4rat,orgaro1rat,orgaro2rat,orgalk1rat,orgole1rat,orgbio1rat, &
        orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog,ncv,nacv,dgnuc,dgacc, &
        fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,deltaso4a,cgrn3,cgra3)


        

      CALL aerostep(layer,blksize,nspcsda,numcells,cblk,dt,so4rat,orgaro1rat, &
        orgaro2rat,orgalk1rat,orgole1rat,orgbio1rat,orgbio2rat,orgbio3rat, &
        orgbio4rat,epm25i,epm25j,eorgi,eorgj,eeci,eecj,esoil,eseas,epmcoarse, &
        dgnuc,dgacc,fconcn,fconca,fconcn_org,fconca_org,pmassn,pmassa,pmassc, &
        dmdt,dndt,deltaso4a,urn00,ura00,brna01,c30,cgrn3,cgra3,igrid,jgrid,kgrid)




      CALL modpar(blksize,nspcsda,numcells,cblk,blkta,blkprs,pmassn,pmassa, &
        pmassc,pdensn,pdensa,pdensc,xlm,amu,dgnuc,dgacc,dgcor,knnuc,knacc, &
        kncor)

      RETURN
    END SUBROUTINE aeroproc




    SUBROUTINE aerostep(layer,blksize,nspcsda,numcells,cblk,dt,so4rat         &
       ,orgaro1rat,orgaro2rat,orgalk1rat,orgole1rat,orgbio1rat,orgbio2rat     &
       ,orgbio3rat,orgbio4rat,epm25i,epm25j,eorgi,eorgj,eeci,eecj,esoil,eseas &
       ,epmcoarse,dgnuc,dgacc,fconcn,fconca,fconcn_org,fconca_org,pmassn      &
       ,pmassa,pmassc,dmdt,dndt,deltaso4a,urn00,ura00,brna01,c30,cgrn3,cgra3, &
        igrid,jgrid,kgrid                                                     &
                                                                             )






















































      INTEGER blksize

      INTEGER numcells

      INTEGER nspcsda

      INTEGER layer
      REAL cblk(blksize,nspcsda) 
      INTEGER igrid,jgrid,kgrid
      REAL dt



      REAL so4rat(blksize) 


      REAL orgaro1rat(blksize)
      REAL orgaro2rat(blksize)


      REAL orgalk1rat(blksize)
      REAL orgole1rat(blksize)


      REAL orgbio1rat(blksize)
      REAL orgbio2rat(blksize)
      REAL orgbio3rat(blksize)
      REAL orgbio4rat(blksize)




      REAL epm25i(blksize) 
      REAL epm25j(blksize) 


      REAL eorgi(blksize) 
      REAL eorgj(blksize) 


      REAL eeci(blksize) 
      REAL eecj(blksize) 


      REAL esoil(blksize) 
      REAL eseas(blksize) 
      REAL epmcoarse(blksize) 

      REAL dgnuc(blksize) 
      REAL dgacc(blksize) 

      REAL fconcn(blksize)                                 

      REAL fconca(blksize)                                 

      REAL fconcn_org(blksize)                                 

      REAL fconca_org(blksize)                                 

      REAL dmdt(blksize)                                 

      REAL dndt(blksize)                                 

      REAL deltaso4a(blksize)                                 

      REAL urn00(blksize) 
      REAL ura00(blksize) 
      REAL brna01(blksize) 
      REAL c30(blksize)       							

      REAL cgrn3(blksize) 
      REAL cgra3(blksize) 



      REAL pmassn(blksize) 
      REAL pmassa(blksize) 
      REAL pmassc(blksize) 




      INTEGER l, lcell, & 
        spc







      REAL*8 a, b, c
      REAL*8 m1, m2, y0, y
      REAL*8 dhat, p, pexpdt, expdt
      REAL*8 loss, prod, pol, lossinv

      REAL mstrnsfr

      REAL factrans


      REAL getaf2
      REAL aaa, xnum, xm3, fnum, fm3, phnum, & 
        phm3
      REAL erf, & 
        erfc

      REAL xx














      erf(xx) = sqrt(1.0-exp(-4.0*xx*xx/pirs))
      erfc(xx) = 1.0 - erf(xx)










      DO l = 1, numcells












        a = urn00(l)
        b = brna01(l)*cblk(l,vac0)
        c = dndt(l) + factnumn*(anthfac*(epm25i(l)+eeci(l))+orgfac*eorgi(l)) 


        y0 = cblk(l,vnu0) 



        IF (c>0.0D0) THEN

          dhat = sqrt(b*b+4.0D0*a*c)

          m1 = 2.0D0*a*c/(b+dhat)

          m2 = -0.5D0*(b+dhat)

          p = -(m1-a*y0)/(m2-a*y0)

          pexpdt = p*exp(-dhat*dt)

          y = (m1+m2*pexpdt)/(a*(1.0D0+pexpdt)) 

        ELSE





          expdt = exp(-b*dt)
          IF (expdt<1.0D0) THEN
            y = b*y0*expdt/(b+a*y0*(1.0D0-expdt))
          ELSE
            y = y0
          END IF

        END IF






        cblk(l,vnu0) = max(nummin_i,y) 






        a = ura00(l)
        b = & 
          0.0D0
        c = factnuma*(anthfac*(epm25j(l)+eecj(l))+orgfac*eorgj(l)) 

        y0 = cblk(l,vac0) 





        IF (c>0.0D0) THEN

          dhat = sqrt(4.0D0*a*c)

          m1 = 2.0D0*a*c/dhat

          m2 = -0.5D0*dhat

          p = -(m1-a*y0)/(m2-a*y0)



          pexpdt = p*exp(-dhat*dt)

          y = (m1+m2*pexpdt)/(a*(1.0D0+pexpdt)) 

        ELSE

          y = y0/(1.0D0+dt*a*y0) 

          y = y0/(1.+dt*a*y0) 


        END IF

        cblk(l,vac0) = max(nummin_j,y) 



        prod = soilfac*esoil(l) + seasfac*eseas(l) + anthfac*epmcoarse(l)


        cblk(l,vcorn) = cblk(l,vcorn) + factnumc*prod*dt






        cgrn3(l) = cgrn3(l) + anthfac*(epm25i(l)+eeci(l)) + orgfac*eorgi(l) 


        cgra3(l) = cgra3(l) + c30(l) + anthfac*(epm25j(l)+eecj(l)) + &
          orgfac*eorgj(l)                                              
                                             













        loss = c30(l)/cblk(l,vnu3) 


        factrans = loss* &                             
          dt
                            



        expdt = exp(-factrans)                               


        lossinv = 1.0/ & 
          loss









        cblk(l,vsulf) = max(conmin,cblk(l,vsulf)-(deltaso4a(l)+dmdt(l)*dt))







        mstrnsfr = cblk(l,vso4ai)*factrans
        prod = deltaso4a(l)*fconcn(l)/dt + dmdt(l) 
        pol = prod*lossinv


        cblk(l,vso4ai) = pol + (cblk(l,vso4ai)-pol)*expdt

        cblk(l,vso4ai) = max(aeroconcmin,cblk(l,vso4ai))

        cblk(l,vso4aj) = cblk(l,vso4aj) + deltaso4a(l)*fconca(l) + mstrnsfr




        mstrnsfr = cblk(l,vorgaro1i)*factrans
        prod = orgaro1rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vorgaro1i) = pol + (cblk(l,vorgaro1i)-pol)*expdt

        cblk(l,vorgaro1i) = max(conmin,cblk(l,vorgaro1i))

        cblk(l,vorgaro1j) = cblk(l,vorgaro1j) + orgaro1rat(l)*fconca_org(l)*dt &
          + mstrnsfr

        mstrnsfr = cblk(l,vorgaro2i)*factrans
        prod = orgaro2rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vorgaro2i) = pol + (cblk(l,vorgaro2i)-pol)*expdt

        cblk(l,vorgaro2i) = max(conmin,cblk(l,vorgaro2i))

        cblk(l,vorgaro2j) = cblk(l,vorgaro2j) + orgaro2rat(l)*fconca_org(l)*dt &
          + mstrnsfr



        mstrnsfr = cblk(l,vorgalk1i)*factrans
        prod = orgalk1rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vorgalk1i) = pol + (cblk(l,vorgalk1i)-pol)*expdt

        cblk(l,vorgalk1i) = max(conmin,cblk(l,vorgalk1i))

        cblk(l,vorgalk1j) = cblk(l,vorgalk1j) + orgalk1rat(l)*fconca_org(l)*dt &
          + mstrnsfr

        mstrnsfr = cblk(l,vorgole1i)*factrans
        prod = orgole1rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vorgole1i) = pol + (cblk(l,vorgole1i)-pol)*expdt

        cblk(l,vorgole1i) = max(conmin,cblk(l,vorgole1i))

        cblk(l,vorgole1j) = cblk(l,vorgole1j) + orgole1rat(l)*fconca_org(l)*dt &
          + mstrnsfr



        mstrnsfr = cblk(l,vorgba1i)*factrans
        prod = orgbio1rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vorgba1i) = pol + (cblk(l,vorgba1i)-pol)*expdt

        cblk(l,vorgba1i) = max(conmin,cblk(l,vorgba1i))

        cblk(l,vorgba1j) = cblk(l,vorgba1j) + orgbio1rat(l)*fconca_org(l)*dt + &
          mstrnsfr

        mstrnsfr = cblk(l,vorgba2i)*factrans
        prod = orgbio2rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vorgba2i) = pol + (cblk(l,vorgba2i)-pol)*expdt

        cblk(l,vorgba2i) = max(conmin,cblk(l,vorgba2i))

        cblk(l,vorgba2j) = cblk(l,vorgba2j) + orgbio2rat(l)*fconca_org(l)*dt + &
          mstrnsfr


        mstrnsfr = cblk(l,vorgba3i)*factrans
        prod = orgbio3rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vorgba3i) = pol + (cblk(l,vorgba3i)-pol)*expdt

        cblk(l,vorgba3i) = max(conmin,cblk(l,vorgba3i))

        cblk(l,vorgba3j) = cblk(l,vorgba3j) + orgbio3rat(l)*fconca_org(l)*dt + &
          mstrnsfr


        mstrnsfr = cblk(l,vorgba4i)*factrans
        prod = orgbio4rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vorgba4i) = pol + (cblk(l,vorgba4i)-pol)*expdt

        cblk(l,vorgba4i) = max(conmin,cblk(l,vorgba4i))

        cblk(l,vorgba4j) = cblk(l,vorgba4j) + orgbio4rat(l)*fconca_org(l)*dt + &
          mstrnsfr



        mstrnsfr = cblk(l,vorgpai)*factrans
        prod = eorgi(l)
        pol = prod*lossinv

        cblk(l,vorgpai) = pol + (cblk(l,vorgpai)-pol)*expdt

        cblk(l,vorgpai) = max(conmin,cblk(l,vorgpai))

        cblk(l,vorgpaj) = cblk(l,vorgpaj) + eorgj(l)*dt + mstrnsfr



        mstrnsfr = cblk(l,vp25ai)*factrans
        prod = epm25i(l)
        pol = prod*lossinv

        cblk(l,vp25ai) = pol + (cblk(l,vp25ai)-pol)*expdt

        cblk(l,vp25ai) = max(conmin,cblk(l,vp25ai))

        cblk(l,vp25aj) = cblk(l,vp25aj) + epm25j(l)*dt + mstrnsfr



        mstrnsfr = cblk(l,veci)*factrans
        prod = eeci(l)
        pol = prod*lossinv

        cblk(l,veci) = pol + (cblk(l,veci)-pol)*expdt

        cblk(l,veci) = max(conmin,cblk(l,veci))

        cblk(l,vecj) = cblk(l,vecj) + eecj(l)*dt + mstrnsfr






        cblk(l,vsoila) = cblk(l,vsoila) + esoil(l)*dt
        cblk(l,vsoila) = max(conmin,cblk(l,vsoila))



        cblk(l,vseas) = cblk(l,vseas) + eseas(l)*dt
        cblk(l,vseas) = max(conmin,cblk(l,vseas))



        cblk(l,vantha) = cblk(l,vantha) + epmcoarse(l)*dt
        cblk(l,vantha) = max(conmin,cblk(l,vantha))



      END DO









      DO lcell = 1, numcells



        IF (cgrn3(lcell)>cgra3(lcell) .OR. dgnuc(lcell)>.03E-6 .AND. cblk( &
            lcell,vnu0)>cblk(lcell,vac0)) & 
            THEN


          aaa = getaf(cblk(lcell,vnu0),cblk(lcell,vac0),dgnuc(lcell), &
            dgacc(lcell),xxlsgn,xxlsga,sqrt2)






          xnum = max(aaa,xxm3)                                    
                                   
                                   


          xm3 = xnum - & 
            xxm3

          IF (xm3>0.0) & 
              THEN

            phnum = 0.5*(1.0+erf(xnum))
            phm3 = 0.5*(1.0+erf(xm3))
            fnum = 0.5*erfc(xnum)
            fm3 = 0.5*erfc(xm3)















            cblk(lcell,vac0) = cblk(lcell,vac0) + fnum*cblk(lcell,vnu0)




            cblk(lcell,vnu0) = phnum*cblk(lcell,vnu0)





            cblk(lcell,vso4aj) = cblk(lcell,vso4aj) + cblk(lcell,vso4ai)*fm3

            cblk(lcell,vnh4aj) = cblk(lcell,vnh4aj) + cblk(lcell,vnh4ai)*fm3

            cblk(lcell,vno3aj) = cblk(lcell,vno3aj) + cblk(lcell,vno3ai)*fm3

            cblk(lcell,vorgaro1j) = cblk(lcell,vorgaro1j) + &
              cblk(lcell,vorgaro1i)*fm3

            cblk(lcell,vorgaro2j) = cblk(lcell,vorgaro2j) + &
              cblk(lcell,vorgaro2i)*fm3

            cblk(lcell,vorgalk1j) = cblk(lcell,vorgalk1j) + &
              cblk(lcell,vorgalk1i)*fm3

            cblk(lcell,vorgole1j) = cblk(lcell,vorgole1j) + &
              cblk(lcell,vorgole1i)*fm3

            cblk(lcell,vorgba1j) = cblk(lcell,vorgba1j) + &
              cblk(lcell,vorgba1i)*fm3

            cblk(lcell,vorgba2j) = cblk(lcell,vorgba2j) + &
              cblk(lcell,vorgba2i)*fm3

            cblk(lcell,vorgba3j) = cblk(lcell,vorgba3j) + &
              cblk(lcell,vorgba3i)*fm3

            cblk(lcell,vorgba4j) = cblk(lcell,vorgba4j) + &
              cblk(lcell,vorgba4i)*fm3

            cblk(lcell,vorgpaj) = cblk(lcell,vorgpaj) + &
              cblk(lcell,vorgpai)*fm3

            cblk(lcell,vp25aj) = cblk(lcell,vp25aj) + cblk(lcell,vp25ai)*fm3

            cblk(lcell,vecj) = cblk(lcell,vecj) + cblk(lcell,veci)*fm3



            cblk(lcell,vso4ai) = cblk(lcell,vso4ai)*phm3


            cblk(lcell,vnh4ai) = cblk(lcell,vnh4ai)*phm3

            cblk(lcell,vno3ai) = cblk(lcell,vno3ai)*phm3

            cblk(lcell,vorgaro1i) = cblk(lcell,vorgaro1i)*phm3

            cblk(lcell,vorgaro2i) = cblk(lcell,vorgaro2i)*phm3

            cblk(lcell,vorgalk1i) = cblk(lcell,vorgalk1i)*phm3

            cblk(lcell,vorgole1i) = cblk(lcell,vorgole1i)*phm3

            cblk(lcell,vorgba1i) = cblk(lcell,vorgba1i)*phm3

            cblk(lcell,vorgba2i) = cblk(lcell,vorgba2i)*phm3

            cblk(lcell,vorgba3i) = cblk(lcell,vorgba3i)*phm3

            cblk(lcell,vorgba4i) = cblk(lcell,vorgba4i)*phm3

            cblk(lcell,vorgpai) = cblk(lcell,vorgpai)*phm3

            cblk(lcell,vp25ai) = cblk(lcell,vp25ai)*phm3

            cblk(lcell,veci) = cblk(lcell,veci)*phm3


          END IF


        END IF


      END DO



      DO spc = 1, nspcsda
        DO lcell = 1, numcells
          cblk(lcell,spc) = max(cblk(lcell,spc),conmin)
        END DO
      END DO


      RETURN


    END SUBROUTINE aerostep

    SUBROUTINE awater(irhx,mso4,mnh4,mno3,wh2o)





















































      INTEGER irhx, irh
      REAL mso4, mnh4, mno3
      REAL tso4, tnh4, tno3, wh2o, x
      REAL aw, awc

      REAL mfs0, mfs1, mfs15, mfs2
      REAL c0(4), c1(4), c15(4), c2(4)
      REAL y, y0, y1, y15, y2, y3, y40, y140, y1540, yc
      REAL kso4(6), kno3(6), mfsso4, mfsno3



      REAL mwso4, mwnh4, mwno3, mw2, mwano3


      PARAMETER (mwso4=96.0636,mwnh4=18.0985,mwno3=62.0049, &
        mw2=mwso4+2.0*mwnh4,mwano3=mwno3+mwnh4)









      DATA c1/0.9995178, -0.7952896, 0.99683673, -1.143874/
      DATA c15/1.697092, -4.045936, 5.833688, -3.463783/
      DATA c2/2.085067, -6.024139, 8.967967, -5.002934/








      DATA c0/0.798079, -1.574367, 2.536686, -1.735297/





      DATA kno3/0.2906, 6.83665, -26.9093, 46.6983, -38.803, 11.8837/
      DATA kso4/2.27515, -11.147, 36.3369, -64.2134, 56.8341, -20.0953/



      irh = irhx
      irh = max(1,irh)
      irh = min(irh,100)
      aw = float(irh)/ & 
        100.0
      tso4 = max(mso4,0.0)
      tnh4 = max(mnh4,0.0)
      tno3 = max(mno3,0.0)
      x = 0.0

      IF (tso4>0.0) THEN
        x = tnh4/tso4
      ELSE

        IF (tno3>0.0 .AND. tnh4>0.0) x = 10.0
      END IF




      IF (x<1.0) THEN

        mfs0 = poly4(c0,aw)
        mfs1 = poly4(c1,aw)
        y0 = (1.0-mfs0)/mfs0
        y1 = (1.0-mfs1)/mfs1
        y = (1.0-x)*y0 + x*y1


      ELSE IF (x<1.5) THEN

        IF (irh>=40) THEN
          mfs1 = poly4(c1,aw)
          mfs15 = poly4(c15,aw)
          y1 = (1.0-mfs1)/mfs1
          y15 = (1.0-mfs15)/mfs15
          y = 2.0*(y1*(1.5-x)+y15*(x-1.0))
        ELSE

























          awc = 0.80*(x-1.0) 
          y = 0.0
          IF (aw>=awc) & 
              THEN
            mfs1 = poly4(c1,0.40)
            mfs15 = poly4(c15,0.40)
            y140 = (1.0-mfs1)/mfs1
            y1540 = (1.0-mfs15)/mfs15
            y40 = 2.0*(y140*(1.5-x)+y1540*(x-1.0))
            yc = 2.0*y1540*(x-1.0) 
            y = y40 - (y40-yc)*(0.40-aw)/(0.40-awc)

          END IF

        END IF

      ELSE IF (x<1.9999) THEN

        y = 0.0
        IF (irh>=40) THEN
          mfs15 = poly4(c15,aw)
          mfs2 = poly4(c2,aw)
          y15 = (1.0-mfs15)/mfs15
          y2 = (1.0-mfs2)/mfs2
          y = 2.0*(y15*(2.0-x)+y2*(x-1.5))

        END IF





      ELSE






        y2 = 0.0
        y3 = 0.0
        IF (irh>=40) THEN
          mfsso4 = poly6(kso4,aw)
          mfsno3 = poly6(kno3,aw)
          y2 = (1.0-mfsso4)/mfsso4
          y3 = (1.0-mfsno3)/mfsno3

        END IF


      END IF





      IF (x<1.9999) THEN

        wh2o = y*(tso4*mwso4+mwnh4*tnh4)

      ELSE




        wh2o = y2*tso4*mw2 + y3*tno3*mwano3

      END IF

      RETURN
    END SUBROUTINE awater


    SUBROUTINE coagrate(blksize,nspcsda,numcells,cblk,blkta,pdensn,pdensa,amu, &
        dgnuc,dgacc,knnuc,knacc,urn00,ura00,brna01,c30)












































      INTEGER blksize

      INTEGER numcells

      INTEGER nspcsda


      REAL cblk(blksize,nspcsda) 
      REAL blkta(blksize) 
      REAL pdensn(blksize) 
      REAL pdensa(blksize) 
      REAL amu(blksize) 
      REAL dgnuc(blksize) 
      REAL dgacc(blksize) 
      REAL knnuc(blksize) 
      REAL knacc(blksize) 



      REAL urn00(blksize) 
      REAL ura00(blksize) 

      REAL brna01(blksize) 
      REAL c30(blksize)                                                               



      REAL*8 kncnuc, & 
        kncacc
      REAL*8 kfmnuc, & 
        kfmacc
      REAL*8 knc, & 
        kfm
      REAL*8 bencnn, & 
        bencna
      REAL*8 & 
        bencm3n
      REAL*8 befmnn, & 
        befmna
      REAL*8 & 
        befm3n
      REAL*8 betann, & 
        betana
      REAL*8 & 
        brna31
      REAL*8 & 
        s1
      REAL*8 t1, & 
        t2
      REAL*8 t16, & 
        t26
      REAL*8 rat, & 
        rin
      REAL*8 rsqt, & 
        rsq4
      REAL*8 rsqti, & 
        rsqi3
      REAL*8 & 
        dgn3
      REAL*8 & 
        dga3



      INTEGER lcell



      REAL*8 bm0
      PARAMETER (bm0=0.8D0)
      REAL*8 bm0i
      PARAMETER (bm0i=0.9D0)
      REAL*8 bm3i
      PARAMETER (bm3i=0.9D0)
      REAL*8 & 
        a
      PARAMETER (a=1.246D0)








      DO lcell = 1, & 
          numcells



        s1 = two3*boltz*blkta(lcell)/amu(lcell)



        kncnuc = s1
        kncacc = s1

        kfmnuc = sqrt(3.0*boltz*blkta(lcell)/pdensn(lcell))
        kfmacc = sqrt(3.0*boltz*blkta(lcell)/pdensa(lcell))



        knc = s1
        kfm = sqrt(6.0*boltz*blkta(lcell)/(pdensn(lcell)+pdensa(lcell)))







        dgn3 = dgnuc(lcell)**3
        dga3 = dgacc(lcell)**3

        t1 = sqrt(dgnuc(lcell))
        t2 = sqrt(dgacc(lcell))
        t16 = & 
          dgn3
        t26 = & 
          dga3




        bencnn = kncnuc*(1.0+esn08+a*knnuc(lcell)*(esn04+esn20))

        bencna = kncacc*(1.0+esa08+a*knacc(lcell)*(esa04+esa20))






        befmnn = kfmnuc*t1*(en1+esn25+2.0*esn05)*bm0

        befmna = kfmacc*t2*(ea1+esa25+2.0*esa05)*bm0







        betann = bencnn*befmnn/(bencnn+befmnn)
        betana = bencna*befmna/(bencna+befmna)



        urn00(lcell) = betann
        ura00(lcell) = betana






        rat = dgacc(lcell)/dgnuc(lcell)
        rin = 1.0D0/rat
        rsqt = sqrt(rat)
        rsq4 = rat**2

        rsqti = 1.0D0/rsqt
        rsqi3 = rin*rsqti




        bencnn = knc*(2.0+a*knnuc(lcell)*(esn04+rat*esn16*esa04)+a*knacc(lcell &
          )*(esa04+rin*esa16*esn04)+(rat+rin)*esn04*esa04)



        bencm3n = knc*dgn3*(2.0*esn36+a*knnuc(lcell)*(esn16+rat*esn04*esa04)+a &
          *knacc(lcell)*(esn36*esa04+rin*esn64*esa16)+rat*esn16*esa04+ &
          rin*esn64*esa04)












        befmnn = kfm*bm0i*t1*(en1+rsqt*ea1+2.0*rat*en1*esa04+rsq4*esn09*esa16+ &
          rsqi3*esn16*esa09+2.0*rsqti*esn04*ea1)



        befm3n = kfm*bm3i*t1*t16*(esn49+rsqt*esn36*ea1+2.0*rat*esn25*esa04+ &
          rsq4*esn09*esa16+rsqi3*esn100*esa09+2.0*rsqti*esn64*ea1)








        brna01(lcell) = bencnn*befmnn/(bencnn+befmnn)

        brna31 = bencm3n* & 
          befm3n/(bencm3n+befm3n)
        c30(lcell) = brna31*cblk(lcell,vac0)*cblk(lcell,vnu0)

                              





      END DO

      RETURN

    END SUBROUTINE coagrate







    SUBROUTINE cubic(a2,a1,a0,nr,crutes)

      INTEGER nr
      REAL*8 a2, a1, a0
      REAL crutes(3)
      REAL*8 qq, rr, a2sq, theta, sqrt3, one3rd
      REAL*8 dum1, dum2, part1, part2, part3, rrsq, phi, yy1, yy2, yy3
      REAL*8 costh, sinth
      DATA sqrt3/1.732050808/, one3rd/0.333333333/

      REAL*8 onebs
      PARAMETER (onebs=1.0)

      a2sq = a2*a2
      qq = (a2sq-3.*a1)/9.
      rr = (a2*(2.*a2sq-9.*a1)+27.*a0)/54.

      dum1 = qq*qq*qq
      rrsq = rr*rr
      dum2 = dum1 - rrsq
      IF (dum2>=0.) THEN

        phi = sqrt(dum1)
        IF (abs(phi)<1.E-20) THEN
          print *, ' cubic phi small, phi = ',phi
          crutes(1) = 0.0
          crutes(2) = 0.0
          crutes(3) = 0.0
          nr = 0
          CALL wrf_error_fatal3("<stdin>",2854,&
'PHI < CRITICAL VALUE')
        END IF
        theta = acos(rr/phi)/3.0
        costh = cos(theta)
        sinth = sin(theta)


        part1 = sqrt(qq)
        yy1 = part1*costh
        yy2 = yy1 - a2/3.0
        yy3 = sqrt3*part1*sinth
        crutes(3) = -2.0*yy1 - a2/3.0
        crutes(2) = yy2 + yy3
        crutes(1) = yy2 - yy3

        IF (crutes(1)<0.0) crutes(1) = 1.0E9
        IF (crutes(2)<0.0) crutes(2) = 1.0E9
        IF (crutes(3)<0.0) crutes(3) = 1.0E9

        crutes(1) = min(crutes(1),crutes(2),crutes(3))
        nr = 3

      ELSE

        part1 = sqrt(rrsq-dum1)
        part2 = abs(rr)
        part3 = (part1+part2)**one3rd
        crutes(1) = -sign(onebs,rr)*(part3+(qq/part3)) - a2/3.

        crutes(2) = 0.
        crutes(3) = 0.



        nr = 1
      END IF
      RETURN

    END SUBROUTINE cubic




    SUBROUTINE eql3(blksize,nspcsda,numcells,cblk,blkta,blkrh,igrid,jgrid,kgrid)



















      INTEGER blksize

      INTEGER numcells

      INTEGER nspcsda,igrid,jgrid,kgrid
      REAL cblk(blksize,nspcsda) 



      REAL blkta(blksize) 
      REAL blkrh(blksize) 



      INTEGER lcell


      REAL temp

      REAL rh

      REAL so4, no3, nh3, nh4, hno3
      REAL aso4, ano3, ah2o, anh4, gnh3, gno3

      REAL fraci

      REAL fracj



      real(kind=8) wi(5),wt(5),wt_save(5)
      real(kind=8) rhi,tempi,cntrl(2)
      real(kind=8) gas(3),aerliq(12),aersld(9),other(6)
      character*15 scasi






      DO lcell = 1, &
          numcells




        temp = blkta(lcell)
        rh = blkrh(lcell)

        rhi = amin1( rh,0.995 )
        tempi = temp
        cntrl(1) = 0.d0         
        cntrl(2) = 0.d0         

        wi(1) = (cblk(lcell,vnaaj)  + cblk(lcell,vnaai))/mw_na_aer*1.e-6      

        wi(2) = (cblk(lcell,vsulf)/(mw_so4_aer+2.) +                         &
                (cblk(lcell,vso4aj) +  cblk(lcell,vso4ai))/mw_so4_aer)*1.e-6        

        wi(3) = (cblk(lcell,vnh3)/(mw_nh4_aer-1.) +                               &
                (cblk(lcell,vnh4aj) +  cblk(lcell,vnh4ai))/mw_nh4_aer)*1.e-6        

        wi(4) = (cblk(lcell,vhno3)/(mw_no3_aer+1.) +                             &
                (cblk(lcell,vno3aj) +  cblk(lcell,vno3ai))/mw_no3_aer)*1.e-6   

       

        wi(5) = (cblk(lcell,vhcl)/(mw_cl_aer+1.) +                               &
                (cblk(lcell,vclaj)  + cblk(lcell,vclai))/mw_cl_aer)*1.e-6      


        wi(1) = max(wi(1),0.)
        wi(2) = max(wi(2),0.)
        wi(3) = max(wi(3),0.)
        wi(4) = max(wi(4),0.)
        wi(5) = max(wi(5),0.)

        wt_save(1) = wi(1) 
        wt_save(2) = wi(2) 
        wt_save(3) = wi(3) 
        wt_save(4) = wi(4) 
        wt_save(5) = wi(5) 
        if(igrid.eq.28.and.jgrid.eq.24.and.kgrid.eq.1)then
         print *,vhcl,vclai
         print *,wi(1),wi(2),wi(3),wi(4),wi(5)
        endif

       call isoropia(wi,rhi,tempi,cntrl,wt,gas,aerliq,aersld,scasi,other)








        gas(1) = min(gas(1),wt_save(3))
        gas(2) = min(gas(2),wt_save(4))
        gas(3) = min(gas(3),wt_save(5))

        gas(1) = max(gas(1),0.)
        gas(2) = max(gas(2),0.)
        gas(3) = max(gas(3),0.)




        cblk(lcell,vnh3)  = gas(1)*1.e6*(mw_nh4_aer-1.)
        cblk(lcell,vhno3) = gas(2)*1.e6*(mw_no3_aer+1.)
        cblk(lcell,vhcl) = gas(3)*1.e6*(mw_cl_aer+1.)
        if(igrid.eq.28.and.jgrid.eq.24.and.kgrid.eq.1)then
         print *,vhcl,vnh3,vhno3
         print *,cblk(lcell,vnh3),cblk(lcell,vhno3),cblk(lcell,vhcl)
        endif


        fraci = cblk(lcell,vso4ai)/(cblk(lcell,vso4aj)+cblk(lcell,vso4ai))

        fraci = min(fraci,1.0)
        fraci = max(fraci,0.0)

        fracj = 1.0 - fraci




        aerliq(8) = max(aerliq(8),0.)

        cblk(lcell,vh2oai) = fraci*aerliq(8)*18.*1.e6
        cblk(lcell,vnh4ai) = fraci*(wt_save(3) - gas(1))*mw_nh4_aer*1.e6
        cblk(lcell,vno3ai) = fraci*(wt_save(4) - gas(2))*mw_no3_aer*1.e6
        cblk(lcell,vclai)  = fraci*(wt_save(5) - gas(3))*mw_cl_aer*1.e6
        cblk(lcell,vnaai)  = fraci*wi(1)*mw_na_aer*1.e6









        cblk(lcell,vh2oaj) = fracj*aerliq(8)*18.*1.e6
        cblk(lcell,vnh4aj) = fracj*(wt_save(3) - gas(1))*mw_nh4_aer*1.e6
        cblk(lcell,vno3aj) = fracj*(wt_save(4) - gas(2))*mw_no3_aer*1.e6
        cblk(lcell,vclaj)  = fracj*(wt_save(5) - gas(3))*mw_cl_aer*1.e6
        cblk(lcell,vnaaj)  = fracj*wi(1)*mw_na_aer*1.e6






        if(igrid.eq.28.and.jgrid.eq.24.and.kgrid.eq.1)then
         print *,vh2oaj,vnh4aj,vno3aj,vclaj,vnaaj
         print *,cblk(lcell,vnh4aj),cblk(lcell,vno3aj),cblk(lcell,vclaj),aerliq(8)
        endif



      END DO


      RETURN


    END SUBROUTINE eql3




    SUBROUTINE eql4(blksize,nspcsda,numcells,cblk,blkta,blkrh)



















      INTEGER blksize

      INTEGER numcells

      INTEGER nspcsda
      REAL cblk(blksize,nspcsda) 



      REAL blkta(blksize) 
      REAL blkrh(blksize) 



      INTEGER lcell


      REAL temp

      REAL rh

      REAL so4, no3, nh3, nh4, hno3
      REAL aso4, ano3, ah2o, anh4, gnh3, gno3

      REAL fraci

      REAL fracj

      DO lcell = 1, &
          numcells




        temp = blkta(lcell)
        rh = blkrh(lcell)







        so4 = cblk(lcell,vso4aj) + cblk(lcell,vso4ai)


        no3 = cblk(lcell,vno3aj) + cblk(lcell,vno3ai)

      
        hno3 = cblk(lcell,vhno3)



        nh3 = cblk(lcell,vnh3)
        
        nh4 = cblk(lcell,vnh4aj) + cblk(lcell,vnh4ai)







        CALL rpmares_old(so4,hno3,no3,nh3,nh4,rh,temp,aso4,ano3,ah2o,anh4, &
          gnh3,gno3)



        fraci = cblk(lcell,vso4ai)/(cblk(lcell,vso4aj)+cblk(lcell,vso4ai))
        fracj = 1.0 - fraci



        cblk(lcell,vh2oai) = fraci*ah2o
        cblk(lcell,vnh4ai) = fraci*anh4
        cblk(lcell,vno3ai) = fraci*ano3



        cblk(lcell,vh2oaj) = fracj*ah2o
        cblk(lcell,vnh4aj) = fracj*anh4
        cblk(lcell,vno3aj) = fracj*ano3



        cblk(lcell,vnh3) = gnh3
        cblk(lcell,vhno3) = gno3

      END DO

      RETURN


    END SUBROUTINE eql4


    SUBROUTINE fdjac(n,x,fjac,ct,cs,imw)






























      INTEGER n
      REAL x(n) 



      REAL ct(np)
      REAL cs(np)
      REAL imw(np)

      REAL fjac(n,n)

      INTEGER i, & 
        j
      REAL a(np)
      REAL b(np)
      REAL b1(np)
      REAL b2(np)
      REAL sum_jnei

      DO i = 1, n
        a(i) = imw(i)
        sum_jnei = 0.
        DO j = 1, n
          sum_jnei = sum_jnei + x(j)*imw(j)
        END DO
        b1(i) = sum_jnei - (x(i)*imw(i))
        b2(i) = cs(i)*imw(i) - ct(i)*imw(i)
        b(i) = b1(i) + b2(i)
      END DO
      DO j = 1, n
        DO i = 1, n
          IF (i==j) THEN
            fjac(i,j) = 2.*a(i)*x(i) + b(i)
          ELSE
            fjac(i,j) = x(i)*imw(j) - ct(i)*imw(j)
          END IF
        END DO
      END DO

      RETURN
    END SUBROUTINE fdjac

    FUNCTION fmin(x,fvec,n,ct,cs,imw,m)





















      INTEGER n


      REAL ct(np)
      REAL cs(np)
      REAL imw(np)
      REAL m,fmin
      REAL x(*), fvec(np)


      INTEGER i
      REAL sum

      CALL funcv(n,x,fvec,ct,cs,imw,m)
      sum = 0.
      DO i = 1, n
        sum = sum + fvec(i)**2
      END DO
      fmin = 0.5*sum
      RETURN
    END FUNCTION fmin

    SUBROUTINE funcv(n,x,fvec,ct,cs,imw,m)













      INTEGER n
      REAL x(*)
      REAL fvec(n)



      REAL ct(np)
      REAL cs(np)
      REAL imw(np)
      REAL m

      INTEGER i, j
      REAL sum_jnei
      REAL a(np)
      REAL b(np)
      REAL c(np)

      DO i = 1, n
        a(i) = imw(i)
        sum_jnei = 0.
        DO j = 1, n
          sum_jnei = sum_jnei + x(j)*imw(j)
        END DO
        sum_jnei = sum_jnei - (x(i)*imw(i))
        b(i) = sum_jnei + cs(i)*imw(i) - ct(i)*imw(i)
        c(i) = -ct(i)*(sum_jnei+m)
        fvec(i) = a(i)*x(i)**2 + b(i)*x(i) + c(i)
      END DO

      RETURN
    END SUBROUTINE funcv
    REAL FUNCTION getaf(ni,nj,dgni,dgnj,xlsgi,xlsgj,sqrt2)


      REAL aa, bb, cc, disc, qq, alfa, l, yji
      REAL ni, nj, dgni, dgnj, xlsgi, xlsgj, sqrt2

      alfa = xlsgi/xlsgj
      yji = log(dgnj/dgni)/(sqrt2*xlsgi)
      aa = 1.0 - alfa*alfa
      l = log(alfa*nj/ni)
      bb = 2.0*yji*alfa*alfa
      cc = l - yji*yji*alfa*alfa
      disc = bb*bb - 4.0*aa*cc
      IF (disc<0.0) THEN
        getaf = - & 
          5.0
        RETURN
      END IF
      qq = -0.5*(bb+sign(1.0,bb)*sqrt(disc))
      getaf = cc/qq
      RETURN

    END FUNCTION getaf










    SUBROUTINE klpnuc(temp,rh,h2so4,ndot1,mdot1,so4rat)






      REAL temp

      REAL rh

      REAL h2so4

      REAL so4rat




      REAL ndot1

      REAL mdot1
                 
      REAL m2dot






      REAL ra


      REAL nav

      REAL nwv

      REAL nav0
                
      REAL nac


      REAL xal
      REAL nsulf, & 
        delta
      REAL*8 & 
        chi
      REAL*8 & 
        jnuc

      REAL tt, & 
        rr
      REAL pi
      PARAMETER (pi=3.14159265)

      REAL pid6
      PARAMETER (pid6=pi/6.0)


      REAL avo
      PARAMETER (avo=6.0221367E23)


      REAL rgasuniv
      PARAMETER (rgasuniv=8.314510)


      REAL atm
      PARAMETER (atm=1013.25E+02)


      REAL mwh2so4
      PARAMETER (mwh2so4=98.07948)


      REAL d35
      PARAMETER (d35=3.5E-07)
      REAL d35sq
      PARAMETER (d35sq=d35*d35)

      REAL v35
      PARAMETER (v35=pid6*d35*d35sq)


      REAL mp


                     
      REAL ugm3_ncm3

      PARAMETER (ugm3_ncm3=(avo/mwh2so4)*1.0E-12)


      REAL nc_ug
      PARAMETER (nc_ug=(1.0E6)*mwh2so4/avo)





      REAL pdens, & 
        rho_p

      REAL ad0, ad1, ad2, & 
        ad3

      PARAMETER (ad0=1.738984,ad1=-1.882301,ad2=2.951849,ad3=-1.810427) 







                
      REAL mp35

      REAL a0, a1, a2, & 
        a3
      PARAMETER (a0=1.961385E2,a1=-5.564447E2,a2=8.828801E2,a3=-5.231409E2)

      REAL ph2so4, &                         
        ph2o



      pdens(rr) = ad0 + rr*(ad1+rr*(ad2+rr*ad3))

      ph2o(tt) = exp(77.34491296-7235.4246512/tt-8.2*log(tt)+tt*5.7113E-03)

      ph2so4(tt) = exp(27.78492066-10156.0/tt)








      mp35(rr) = nc_ug*(a0+rr*(a1+rr*(a2+rr*a3)))









      nwv = rh*ph2o(temp)/(rgasuniv*temp)*avo*1.0E-6








      nav0 = ph2so4(temp)/(rgasuniv*temp)*avo*1.0E-6




      nav = ugm3_ncm3*h2so4




      nac = exp(-14.5125+0.1335*temp-10.5462*rh+1958.4*rh/temp)



      ra = nav/nav0



      delta = 1.0 + (temp-273.15)/273.14



      xal = 1.2233 - 0.0154*ra/(ra+rh) + 0.0102*log(nav) - 0.0415*log(nwv) + &
        0.0016*temp


      nsulf = log(nav/nac)



      chi = 25.1289*nsulf - 4890.8*nsulf/temp - 1743.3/temp - &
        2.2479*delta*nsulf*rh + 7643.4*xal/temp - 1.9712*xal*delta/rh

      jnuc = exp(chi) 

      ndot1 = (1.0E06)*jnuc







      rho_p = pdens(rh)





      mp = mp35(rh)                      










      mdot1 = mp*ndot1



      IF (mdot1>so4rat) THEN

        mdot1 = & 
          so4rat

        ndot1 = mdot1/ & 
          mp

      END IF


      IF (mdot1==0.) ndot1 = 0.



      m2dot = 1.0E-04*d35sq*ndot1

      RETURN

    END SUBROUTINE klpnuc
    SUBROUTINE lnsrch(ctot,n,xold,fold,g,p,x,f,stpmax,check,func, &
     fvec,ct,cs,imw,m)




























      INTEGER n
      LOGICAL check
      REAL f, fold, stpmax
      REAL g(n), p(n), x(n), xold(n)
      REAL func
      REAL ctot(n)
      REAL alf
      REAL ct(np)
      REAL cs(np)
      REAL imw(np)
      REAL fvec(n)
      REAL m

      PARAMETER (alf=1.E-04)

      EXTERNAL func

      INTEGER i
      REAL a, alam, alam2, alamin, b, disc
      REAL f2, fold2, rhs1, rhs2, slope
      REAL sum, temp, test, tmplam

      check = .FALSE.
      sum = 0.
      DO i = 1, n
        sum = sum + p(i)*p(i)
      END DO
      sum = sqrt(sum)
      IF (sum>stpmax) THEN
        DO i = 1, n
          p(i) = p(i)*stpmax/sum
        END DO
      END IF
      slope = 0.
      DO i = 1, n
        slope = slope + g(i)*p(i)
      END DO
      test = 0.
      DO i = 1, n
        temp = abs(p(i))/max(abs(xold(i)),1.)
        IF (temp>test) test = temp
      END DO
      alamin = tolx/test
      alam = 1.

10    CONTINUE




      DO i = 1, n
        x(i) = xold(i) + alam*p(i)
        IF (x(i)<=0.) x(i) = conmin
        IF (x(i)>ctot(i)) x(i) = ctot(i)
      END DO
      f = func(x,fvec,n,ct,cs,imw,m)
      IF (alam<alamin) THEN
        DO i = 1, n
          x(i) = xold(i)
        END DO
        check = .TRUE.
        RETURN
      ELSE IF (f<=fold+alf*alam*slope) THEN
        RETURN
      ELSE
        IF (alam==1.) THEN
          tmplam = -slope/(2.*(f-fold-slope))
        ELSE
          rhs1 = f - fold - alam*slope
          rhs2 = f2 - fold2 - alam2*slope
          a = (rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
          b = (-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
          IF (a==0.) THEN
            tmplam = -slope/(2.*b)
          ELSE
            disc = b*b - 3.*a*slope
            tmplam = (-b+sqrt(disc))/(3.*a)
          END IF
          IF (tmplam>0.5*alam) tmplam = 0.5*alam
        END IF
      END IF
      alam2 = alam
      f2 = f
      fold2 = fold
      alam = max(tmplam,0.1*alam)
      GO TO 10

    END SUBROUTINE lnsrch

    SUBROUTINE lubksb(a,n,np,indx,b)






















      INTEGER n, np, indx(n)
      REAL a(np,np), b(n)

      INTEGER i, ii, j, ll
      REAL sum

      ii = 0
      DO i = 1, n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        IF (ii/=0) THEN
          DO j = ii, i - 1
            sum = sum - a(i,j)*b(j)
          END DO
        ELSE IF (sum/=0) THEN
          ii = i
        END IF
        b(i) = sum
      END DO
      DO i = n, 1, -1
        sum = b(i)
        DO j = i + 1, n
          sum = sum - a(i,j)*b(j)
        END DO
        b(i) = sum/a(i,i)
      END DO

      RETURN
    END SUBROUTINE lubksb

    SUBROUTINE ludcmp(a,n,np,indx,d,klev)




























      INTEGER n, np, indx(n)
      INTEGER nmax
      PARAMETER (nmax=10) 
      REAL d, a(np,np)
      REAL tiny
      PARAMETER (tiny=1.0E-20)

      INTEGER i, imax, j, k
      REAL aamax, dum, sum, vv(nmax)
      integer klev

      d = 1
      DO i = 1, n
        aamax = 0.
        DO j = 1, n
          IF (abs(a(i,j))>aamax) aamax = abs(a(i,j))
        END DO
        IF (aamax==0) THEN
          print *, 'Singular matrix in ludcmp, klev = ',klev
          a(1,1)=epsilc

        END IF
        vv(i) = 1./aamax
      END DO
      DO j = 1, n
        DO i = 1, j - 1
          sum = a(i,j)
          DO k = 1, i - 1
            sum = sum - a(i,k)*a(k,j)
          END DO
          a(i,j) = sum
        END DO
        aamax = 0.
        DO i = j, n
          sum = a(i,j)
          DO k = 1, j - 1
            sum = sum - a(i,k)*a(k,j)
          END DO
          a(i,j) = sum
          dum = vv(i)*abs(sum)
          IF (dum>=aamax) THEN
            imax = i
            aamax = dum
          END IF
        END DO
        IF (j/=imax) THEN
          DO k = 1, n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          END DO
          d = -d
          vv(imax) = vv(j)
        END IF
        indx(j) = imax
        IF (a(j,j)==0.) a(j,j) = tiny
        IF (j/=n) THEN
          dum = 1./a(j,j)
          DO i = j + 1, n
            a(i,j) = a(i,j)*dum
          END DO
        END IF
      END DO

      RETURN
    END SUBROUTINE ludcmp



    SUBROUTINE modpar(blksize,nspcsda,numcells,cblk,blkta,blkprs,pmassn, &
        pmassa,pmassc,pdensn,pdensa,pdensc,xlm,amu,dgnuc,dgacc,dgcor,knnuc, &
        knacc,kncor)


























      INTEGER blksize

      INTEGER numcells

      INTEGER nspcsda


      REAL cblk(blksize,nspcsda) 
      REAL blkta(blksize) 
      REAL blkprs(blksize) 





      REAL dgmin
      PARAMETER (dgmin=1.0E-09)


      REAL densmin
      PARAMETER (densmin=1.0E03)

      REAL pmassn(blksize) 
      REAL pmassa(blksize) 
      REAL pmassc(blksize) 
      REAL pdensn(blksize) 
      REAL pdensa(blksize) 
      REAL pdensc(blksize) 
      REAL xlm(blksize) 
      REAL amu(blksize) 
      REAL dgnuc(blksize) 
      REAL dgacc(blksize) 
      REAL dgcor(blksize) 
      REAL knnuc(blksize) 
      REAL knacc(blksize) 
      REAL kncor(blksize) 



      INTEGER lcell





      DO lcell = 1, numcells



        cblk(lcell,vnu3) = so4fac*cblk(lcell, &
          vso4ai)+nh4fac*cblk(lcell,vnh4ai)+h2ofac*cblk(lcell, &
          vh2oai)+no3fac*cblk(lcell,vno3ai)+                   &
          nafac*cblk(lcell,vnaai)+  clfac*cblk(lcell,vclai)+   &
          orgfac*cblk(lcell, &
          vorgaro1i)+orgfac*cblk(lcell,vorgaro2i)+orgfac*cblk(lcell, &
          vorgalk1i)+orgfac*cblk(lcell,vorgole1i)+orgfac*cblk(lcell, &
          vorgba1i)+orgfac*cblk(lcell,vorgba2i)+orgfac*cblk(lcell, &
          vorgba3i)+orgfac*cblk(lcell,vorgba4i)+orgfac*cblk(lcell, &
          vorgpai)+anthfac*cblk(lcell,vp25ai)+anthfac*cblk(lcell,veci)







        cblk(lcell,vac3) = so4fac*cblk(lcell, &
          vso4aj)+nh4fac*cblk(lcell,vnh4aj)+h2ofac*cblk(lcell, &
          vh2oaj)+no3fac*cblk(lcell,vno3aj) +                  &
          nafac*cblk(lcell,vnaaj)+  clfac*cblk(lcell,vclaj)+   &
          orgfac*cblk(lcell, &
          vorgaro1j)+orgfac*cblk(lcell,vorgaro2j)+orgfac*cblk(lcell, &
          vorgalk1j)+orgfac*cblk(lcell,vorgole1j)+orgfac*cblk(lcell, &
          vorgba1j)+orgfac*cblk(lcell,vorgba2j)+orgfac*cblk(lcell, &
          vorgba3j)+orgfac*cblk(lcell,vorgba4j)+orgfac*cblk(lcell, &
          vorgpaj)+anthfac*cblk(lcell,vp25aj)+anthfac*cblk(lcell,vecj)






        cblk(lcell,vcor3) = soilfac*cblk(lcell, &
          vsoila)+seasfac*cblk(lcell,vseas)+anthfac*cblk(lcell,vantha)







        pmassn(lcell) = max(conmin,(cblk(lcell,vso4ai)+cblk(lcell, &
          vnh4ai)+cblk(lcell,vh2oai)+cblk(lcell,vno3ai)+ &
          cblk(lcell,vnaai)+cblk(lcell,vclai)+cblk(lcell, &
          vorgaro1i)+cblk(lcell,vorgaro2i)+cblk(lcell,vorgalk1i)+cblk(lcell, &
          vorgole1i)+cblk(lcell,vorgba1i)+cblk(lcell,vorgba2i)+cblk(lcell, &
          vorgba3i)+cblk(lcell,vorgba4i)+cblk(lcell,vorgpai)+cblk(lcell, &
          vp25ai)+cblk(lcell,veci)))













        pmassa(lcell) = max(conmin,(cblk(lcell,vso4aj)+cblk(lcell, &
          vnh4aj)+cblk(lcell,vh2oaj)+cblk(lcell,vno3aj)+ &
          cblk(lcell,vnaaj)+cblk(lcell,vclaj)+cblk(lcell, &
          vorgaro1j)+cblk(lcell,vorgaro2j)+cblk(lcell,vorgalk1j)+cblk(lcell, &
          vorgole1j)+cblk(lcell,vorgba1j)+cblk(lcell,vorgba2j)+cblk(lcell, &
          vorgba3j)+cblk(lcell,vorgba4j)+cblk(lcell,vorgpaj)+cblk(lcell, &
          vp25aj)+cblk(lcell,vecj)))











        pmassc(lcell) = max(conmin,cblk(lcell,vsoila)+cblk(lcell,vseas)+cblk( &
          lcell,vantha))



      END DO



      DO lcell = 1, & 
          numcells



        pdensn(lcell) = max(densmin,(f6dpim9*pmassn(lcell)/cblk(lcell,vnu3)))
        pdensa(lcell) = max(densmin,(f6dpim9*pmassa(lcell)/cblk(lcell,vac3)))
        pdensc(lcell) = max(densmin,(f6dpim9*pmassc(lcell)/cblk(lcell,vcor3)))



        xlm(lcell) = 6.6328E-8*pss0*blkta(lcell)/(tss0*blkprs(lcell))











        amu(lcell) = 1.458E-6*blkta(lcell)*sqrt(blkta(lcell))/ &
          (blkta(lcell)+110.4)


      END DO






      DO lcell = 1, & 
          numcells


        dgnuc(lcell) = max(dgmin,(cblk(lcell,vnu3)/(cblk(lcell,vnu0)*esn36))** &
          one3)


        dgacc(lcell) = max(dgmin,(cblk(lcell,vac3)/(cblk(lcell,vac0)*esa36))** &
          one3)


        dgcor(lcell) = max(dgmin,(cblk(lcell,vcor3)/(cblk(lcell,vcorn)*esc36)) &
          **one3)



      if (cw_phase > 0) then
        dgnuc(lcell) = max( dgnuc(lcell), dginin*0.2  )  
        dgnuc(lcell) = min( dgnuc(lcell), dginin*10.0 )  
        dgacc(lcell) = max( dgacc(lcell), dginia*0.2  )  
        dgacc(lcell) = min( dgacc(lcell), dginia*10.0 )  
        dgcor(lcell) = max( dgcor(lcell), dginic*0.2  )  
        dgcor(lcell) = min( dgcor(lcell), dginic*10.0 )  
      end if

      END DO

      DO lcell = 1, & 
          numcells

        knnuc(lcell) = 2.0*xlm(lcell)/dgnuc(lcell)

        knacc(lcell) = 2.0*xlm(lcell)/dgacc(lcell)

        kncor(lcell) = 2.0*xlm(lcell)/dgcor(lcell)


      END DO


      RETURN


    END SUBROUTINE modpar

    SUBROUTINE newt(layer,x,n,check,ctot,csat,imwcv,minitw,its)











































      INTEGER layer

      INTEGER n
      REAL x(n) 
      LOGICAL check
      REAL ctot(n) 
      REAL csat(n) 
      REAL imwcv(n) 

      REAL minitw



      INTEGER nn


      REAL fvec(np) 


      REAL ct(np)
      REAL cs(np)
      REAL imw(np)
      REAL m

      INTEGER i, its, j, indx(np)
      REAL d, den, f, fold, stpmax, sum, temp, test
      REAL fjac(np,np)
      REAL g(np), p(np), xold(np)





      m = minitw
      DO i = 1, n
        ct(i) = ctot(i)
        cs(i) = csat(i)
        imw(i) = imwcv(i)
      END DO

      nn = n
      f = fmin(x,fvec,nn,ct,cs,imw,m) 
      test = & 
        0.
      DO i = 1, & 
          n
        IF (abs(fvec(i))>test) test = abs(fvec(i))
      END DO
      IF (test<0.01*tolf) RETURN
      sum = & 
        0.
      DO i = 1, n
        sum = sum + x(i)**2
      END DO
      stpmax = stpmx*max(sqrt(sum),float(n))
      DO its = 1, & 
          maxits
        CALL fdjac(n,x,fjac,ct,cs,imw) 
        DO i = 1, & 
            n
          sum = 0.
          DO j = 1, n
            sum = sum + fjac(j,i)*fvec(j)
          END DO
          g(i) = sum
        END DO
        DO i = 1, & 
            n
          xold(i) = x(i)
        END DO
        fold = & 
          f
        DO i = 1, & 
            n
          p(i) = -fvec(i)
        END DO
        CALL ludcmp(fjac,n,np,indx,d,layer) 
        CALL lubksb(fjac,n,np,indx,p)
        CALL lnsrch(ctot,n,xold,fold,g, & 
          p,x,f,stpmax, & 
          check,fmin,fvec,ct,cs,imw,m) 
        test = 0.
        DO i = 1, n
          IF (abs(fvec(i))>test) test = abs(fvec(i))
        END DO
        IF (test<tolf) THEN
          check = .FALSE.
          RETURN
        END IF
        IF (check) & 
            THEN
          test = & 
            0.
          den = max(f,0.5*n)
          DO i = 1, n
            temp = abs(g(i))*max(abs(x(i)),1.)/den
            IF (temp>test) test = temp
          END DO
          IF (test<tolmin) THEN
            check = .TRUE.
          ELSE
            check = .FALSE.
          END IF
          RETURN
        END IF
        test = & 
          0.
        DO i = 1, n
          temp = (abs(x(i)-xold(i)))/max(abs(x(i)),1.)
          IF (temp>test) test = temp
        END DO
        IF (test<tolx) RETURN
      END DO


    END SUBROUTINE newt


    SUBROUTINE nuclcond(blksize,nspcsda,numcells,cblk,dt,layer,blkta,blkprs, &
        blkrh,so4rat,orgaro1rat,orgaro2rat,orgalk1rat,orgole1rat,orgbio1rat, &
        orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog,ncv,nacv,dgnuc,dgacc, &
        fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,deltaso4a,cgrn3,cgra3)



































      INTEGER blksize
      INTEGER layer

      INTEGER nspcsda

      INTEGER numcells

      INTEGER ldrog

      REAL cblk(blksize,nspcsda) 

      REAL dt
      REAL blkta(blksize) 
      REAL blkprs(blksize) 
      REAL blkrh(blksize) 
      REAL so4rat(blksize)                                       



      INTEGER ncv

      INTEGER nacv


      REAL drog(blksize,ldrog) 

      REAL orgaro1rat(blksize)                                 

      REAL orgaro2rat(blksize)                                 

      REAL orgalk1rat(blksize)                                 

      REAL orgole1rat(blksize)                                 


      REAL orgbio1rat(blksize)                                 

      REAL orgbio2rat(blksize)                                 

      REAL orgbio3rat(blksize)                                 

      REAL orgbio4rat(blksize)                                 


      REAL dgnuc(blksize) 
      REAL dgacc(blksize) 



      REAL fconcn(blksize)                                 

      REAL fconca(blksize)                                 

      REAL fconcn_org(blksize)                                 

      REAL fconca_org(blksize)                                 

      REAL dmdt(blksize)                                 

      REAL dndt(blksize)                                 

      REAL deltaso4a(blksize)                                 

      REAL cgrn3(blksize)                                 

      REAL cgra3(blksize)                                 





      INTEGER lcell



      REAL chemrat

      REAL chemrat_org
      REAL am1n, & 
        am1a
      REAL am2n, & 
        am2a
      REAL gnc3n, & 
        gnc3a
      REAL gfm3n, & 
        gfm3a

      REAL fconc

      REAL td

      REAL*8 & 
        one88
      PARAMETER (one88=1.0D0)



      REAL vapor1

      REAL vapor2


      REAL deltavap



      REAL diffcorr


      REAL csqt_org



      REAL csqt







      DO lcell = 1, & 
          numcells



        am1n = cblk(lcell,vnu0)*dgnuc(lcell)*esn04
        am1a = cblk(lcell,vac0)*dgacc(lcell)*esa04






        diffcorr = (pss0/blkprs(lcell))*(blkta(lcell)/273.16)**1.

        gnc3n = cconc*am1n*diffcorr
        gnc3a = cconc*am1a*diffcorr




        am2n = cblk(lcell,vnu0)*dgnuc(lcell)*dgnuc(lcell)*esn16
        am2a = cblk(lcell,vac0)*dgacc(lcell)*dgacc(lcell)*esa16

        csqt = ccofm*sqrt(blkta(lcell)) 



        gfm3n = csqt*am2n
        gfm3a = csqt*am2a









        fconcn(lcell) = one88*gnc3n*gfm3n/(gnc3n+gfm3n)
        fconca(lcell) = one88*gnc3a*gfm3a/(gnc3a+gfm3a)
        fconc = fconcn(lcell) + fconca(lcell)





        gnc3n = cconc_org*am1n*diffcorr
        gnc3a = cconc_org*am1a*diffcorr

        csqt_org = ccofm_org*sqrt(blkta(lcell))
        gfm3n = csqt_org*am2n
        gfm3a = csqt_org*am2a

        fconcn_org(lcell) = one88*gnc3n*gfm3n/(gnc3n+gfm3n)
        fconca_org(lcell) = one88*gnc3a*gfm3a/(gnc3a+gfm3a)






        vapor1 = cblk(lcell,vsulf) 
        vapor2 = cblk(lcell,vsulf) - so4rat(lcell)* & 
          dt

        vapor2 = max(0.0,vapor2)

        deltavap = max(0.0,(so4rat(lcell)/fconc-vapor2)*(1.0-exp(-fconc*dt)))











        deltaso4a(lcell) = max(-1.*cblk(lcell,vsulf), &
          so4rat(lcell)*dt-deltavap)




        cgrn3(lcell) = 0.0
        cgra3(lcell) = 0.0


      END DO




      IF (inucl==1) THEN













      ELSE IF (inucl==0) THEN













      ELSE IF (inucl==2) THEN




        if(blkta(1).ge.233.15.and.blkrh(1).ge.0.1)then
           CALL klpnuc(blkta(1),blkrh(1),vapor1,dndt(1),dmdt(1),so4rat(1))
        else
           dndt(1)=0.
           dmdt(1)=0.
        endif




        IF (dndt(1)==0.) dmdt(1) = 0.
        IF (dmdt(1)==0.) dndt(1) = 0.













      END IF




      CALL sorgam(layer,blkta,blkprs,orgaro1rat,orgaro2rat,orgalk1rat, &
        orgole1rat,orgbio1rat,orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog,ncv, &
        nacv,cblk,blksize,nspcsda,numcells,dt)




      DO lcell = 1, numcells




        td = 1.0/(fconcn(lcell)+fconca(lcell))
        fconcn(lcell) = td*fconcn(lcell)
        fconca(lcell) = td*fconca(lcell)

        td = 1.0/(fconcn_org(lcell)+fconca_org(lcell))
        fconcn_org(lcell) = td*fconcn_org(lcell)
        fconca_org(lcell) = td*fconca_org(lcell)

      END DO



      DO lcell = 1, & 
          numcells



        chemrat = so4fac*so4rat(lcell) 
        chemrat_org = orgfac*(orgaro1rat(lcell)+orgaro2rat(lcell)+orgalk1rat( &
          lcell)+orgole1rat(lcell)+orgbio1rat(lcell)+orgbio2rat(lcell)+ &
          orgbio3rat(lcell)+orgbio4rat(lcell)) 



        cgrn3(lcell) = so4fac*dmdt(lcell) 

        chemrat = chemrat - cgrn3(lcell)                                            


        chemrat = max(chemrat,0.0) 



        cgrn3(lcell) = cgrn3(lcell) + chemrat*fconcn(lcell) + &
          chemrat_org*fconcn_org(lcell)

        cgra3(lcell) = chemrat*fconca(lcell) + chemrat_org*fconca_org(lcell)



      END DO

      RETURN

    END SUBROUTINE nuclcond



    REAL FUNCTION poly4(a,x)
      REAL a(4), x

      poly4 = a(1) + x*(a(2)+x*(a(3)+x*(a(4))))
      RETURN
    END FUNCTION poly4
    REAL FUNCTION poly6(a,x)
      REAL a(6), x

      poly6 = a(1) + x*(a(2)+x*(a(3)+x*(a(4)+x*(a(5)+x*(a(6))))))
      RETURN
    END FUNCTION poly6






    SUBROUTINE rpmares_old(so4,hno3,no3,nh3,nh4,rh,temp,aso4,ano3,ah2o,anh4, &
        gnh3,gno3)
































































































































      REAL mwnacl
      PARAMETER (mwnacl=58.44277)


      REAL mwno3
      PARAMETER (mwno3=62.0049)


      REAL mwhno3
      PARAMETER (mwhno3=63.01287)


      REAL mwso4
      PARAMETER (mwso4=96.0576)


      REAL mwhso4
      PARAMETER (mwhso4=mwso4+1.0080)


      REAL mh2so4
      PARAMETER (mh2so4=98.07354)


      REAL mwnh3
      PARAMETER (mwnh3=17.03061)


      REAL mwnh4
      PARAMETER (mwnh4=18.03858)


      REAL mworg

      PARAMETER (mworg=175.0)



      REAL mwcl
      PARAMETER (mwcl=35.453)


      REAL mwair
      PARAMETER (mwair=28.964)


      REAL mwlct
      PARAMETER (mwlct=3.0*mwnh4+2.0*mwso4+1.0080)


      REAL mwas
      PARAMETER (mwas=2.0*mwnh4+mwso4)


      REAL mwabs
      PARAMETER (mwabs=mwnh4+mwso4+1.0080)




      REAL so4


      REAL hno3

      REAL no3

      REAL nh3

      REAL nh4

      REAL rh

      REAL temp

      REAL aso4

      REAL ano3

      REAL ah2o

      REAL anh4

      REAL gno3

      REAL gnh3




      INTEGER irh

      INTEGER nitr

      INTEGER nnn

      INTEGER nr

      REAL*8 & 
        a0
      REAL*8 & 
        a1
      REAL*8 & 
        a2

      REAL aa

      REAL bal

      REAL bb

      REAL bhat

      REAL cc

      REAL convt

      REAL dd

      REAL disc

      REAL eror

      REAL fnh3

      REAL gamaab

      REAL gamaan

      REAL gamahat

      REAL gamana

      REAL gamas1

      REAL gamas2

      REAL gamold

      REAL gasqd

      REAL hplus

      REAL k1a

      REAL k2sa

      REAL k3

      REAL kan

      REAL khat

      REAL kna

      REAL kph

      REAL kw

      REAL kw2

      REAL man

      REAL mas

      REAL mhso4

      REAL mna

      REAL mnh4

      REAL molnu

      REAL mso4

      REAL phibar

      REAL phiold

      REAL ratio

      REAL rk2sa

      REAL rkna

      REAL rknwet
      REAL rr1
      REAL rr2

      REAL stion

      REAL t1

      REAL t2

      REAL t21

      REAL t221

      REAL t3

      REAL t4

      REAL t6

      REAL tnh4

      REAL tno3

      REAL toler1

      REAL toler2

      REAL tso4

      REAL twoso4

      REAL wfrac
                                   
      REAL wh2o
                                   
                                   
                                   


      REAL wsqd

      REAL xno3

      REAL xxq

      REAL ynh4

      REAL zh2o

      REAL zso4

      REAL cat(2) 
      REAL an(3) 
      REAL crutes(3) 
      REAL gams(2,3) 

      REAL minso4
      PARAMETER (minso4=1.0E-6/mwso4)
      REAL floor
      PARAMETER (floor=1.0E-30) 







      tso4 = max(0.0,so4/mwso4)
      tno3 = max(0.0,(no3/mwno3+hno3/mwhno3))
      tnh4 = max(0.0,(nh3/mwnh3+nh4/mwnh4))




      irh = nint(100.0*rh)



      irh = max(1,irh)
      irh = min(99,irh)







      convt = 1.0/(0.082*temp)
      t6 = 0.082E-9*temp
      t1 = 298.0/temp
      t2 = alog(t1)
      t3 = t1 - 1.0
      t4 = 1.0 + t2 - t1
      kna = 2.511E+06*exp(29.17*t3+16.83*t4)*t6
      k1a = 1.805E-05*exp(-1.50*t3+26.92*t4)
      k2sa = 1.015E-02*exp(8.85*t3+25.14*t4)
      kw = 1.010E-14*exp(-22.52*t3+26.92*t4)
      kph = 57.639*exp(13.79*t3-5.39*t4)*t6

      khat = kph*k1a/kw
      kan = kna*khat




      k3 = exp(118.87-24084.0/temp-6.025*alog(temp))



      k3 = k3*convt*convt

      wh2o = 0.0
      stion = 0.0
      ah2o = 0.0
      mas = 0.0
      man = 0.0
      hplus = 0.0
      toler1 = 0.00001
      toler2 = 0.001
      nitr = 0
      nr = 0
      ratio = 0.0
      gamaan = 1.0
      gamold = 1.0


      IF (tso4>minso4) THEN
        ratio = tnh4/tso4




      ELSE

        IF (tno3==0.0) THEN



          aso4 = max(floor,aso4)
          ano3 = max(floor,ano3)
          wh2o = 0.0
          ah2o = 0.0
          gnh3 = max(floor,gnh3)
          gno3 = max(floor,gno3)
          RETURN
        END IF




        ratio = 5.0
      END IF





      IF (ratio>2.0) THEN

        gamaan = 0.1



        twoso4 = 2.0*tso4
        xno3 = 0.0
        ynh4 = twoso4







        CALL awater(irh,tso4,ynh4,tno3,ah2o) 
        wh2o = 1.0E-3*ah2o
        aso4 = tso4*mwso4
        ano3 = 0.0
        anh4 = ynh4*mwnh4
        wfrac = ah2o/(aso4+anh4+ah2o)

        IF (wfrac<0.2) THEN




          fnh3 = tnh4 - twoso4
          cc = tno3*fnh3 - k3



          IF (cc<=0.0) THEN
            xno3 = 0.0
          ELSE
            aa = 1.0
            bb = -(tno3+fnh3)
            disc = bb*bb - 4.0*cc




            IF (disc<0.0) THEN
              xno3 = 0.0
              ah2o = 1000.0*wh2o
              ynh4 = twoso4
              gno3 = tno3*mwhno3
              gnh3 = (tnh4-ynh4)*mwnh3
              aso4 = tso4*mwso4
              ano3 = 0.0
              anh4 = ynh4*mwnh4
              RETURN
            END IF



            dd = sqrt(disc)
            xxq = -0.5*(bb+sign(1.0,bb)*dd)



            xno3 = min(xxq/aa,cc/xxq)

          END IF
          ah2o = 1000.0*wh2o
          ynh4 = 2.0*tso4 + xno3
          gno3 = (tno3-xno3)*mwhno3
          gnh3 = (tnh4-ynh4)*mwnh3
          aso4 = tso4*mwso4
          ano3 = xno3*mwno3
          anh4 = ynh4*mwnh4
          RETURN

        END IF




        mas = tso4/wh2o
        man = 0.0
        xno3 = 0.0
        ynh4 = twoso4
        phiold = 1.0






        DO nnn = 1, 150
          nitr = nnn
          gasqd = gamaan*gamaan
          wsqd = wh2o*wh2o
          kw2 = kan*wsqd/gasqd
          aa = 1.0 - kw2
          bb = twoso4 + kw2*(tno3+tnh4-twoso4)
          cc = -kw2*tno3*(tnh4-twoso4)



          disc = bb*bb - 4.0*aa*cc



          IF((disc<0.0) .OR. &
             (bb>0.0 .AND. aa>0.0 .AND.cc >0.0) .OR. &
             (bb<0.0 .AND. aa<0.0 .AND.cc <0.0)) THEN
            xno3 = 0.0
            ah2o = 1000.0*wh2o
            ynh4 = twoso4
            gno3 = tno3*mwhno3
            gnh3 = (tnh4-ynh4)*mwnh3
            aso4 = tso4*mwso4
            ano3 = 0.0
            anh4 = ynh4*mwnh4

            RETURN
          END IF

          dd = sqrt(disc)
          xxq = -0.5*(bb+sign(1.0,bb)*dd)

          aa=max(aa,1.e-20)
          xxq=max(xxq,1.e-20)

          rr1 = xxq/aa
          rr2 = cc/xxq

          IF (rr1 <= 0.0 .AND. rr2 <= 0.0) THEN
            xno3 = 0.0
            ah2o = 1000.0*wh2o
            ynh4 = twoso4
            gno3 = tno3*mwhno3
            gnh3 = (tnh4-ynh4)*mwnh3
            aso4 = tso4*mwso4
            ano3 = 0.0
            anh4 = ynh4*mwnh4

            RETURN
          END IF



          IF ((rr1*rr2)<0.0) THEN
            xno3 = max(rr1,rr2)
          ELSE
            xno3 = min(rr1,rr2)
          END IF

          xno3 = min(xno3,tno3)




          CALL awater(irh,tso4,ynh4,xno3,ah2o)





          wh2o = 1.0E-3*ah2o



          man = xno3/wh2o
          mas = tso4/wh2o
          mnh4 = 2.0*mas + man
          ynh4 = mnh4*wh2o




          stion = 3.0*mas + man
          cat(1) = 0.0
          cat(2) = max(mnh4,0.0)
          an(1) = max(mas,0.0)
          an(2) = max(man,0.0)
          an(3) = 0.0
          CALL actcof(cat,an,gams,molnu,phibar)
          gamaan = gams(2,2)



          eror = abs(gamold-gamaan)/gamold
          gamold = gamaan



          IF (eror<=toler1) THEN



            aso4 = tso4*mwso4
            ano3 = xno3*mwno3
            anh4 = ynh4*mwnh4
            gno3 = (tno3-xno3)*mwhno3
            gnh3 = (tnh4-ynh4)*mwnh3
            ah2o = 1000.0*wh2o
            RETURN
          END IF

        END DO



        aso4 = tso4*mwso4
        ano3 = 0.0
        ynh4 = twoso4
        anh4 = ynh4*mwnh4
        CALL awater(irh,tso4,ynh4,xno3,ah2o)
        gno3 = tno3*mwhno3
        gnh3 = (tnh4-ynh4)*mwnh3
        RETURN

      ELSE








        wh2o = 0.0
        CALL awater(irh,tso4,tnh4,tno3,ah2o)
        wh2o = 1.0E-3*ah2o
        zh2o = ah2o



        aso4 = tso4*mwso4
        anh4 = tnh4*mwnh4
        ano3 = 0.0
        gno3 = tno3*mwhno3
        gnh3 = 0.0



        IF (wh2o==0.0) RETURN
        zso4 = tso4/wh2o









        IF (zso4>9.0) & 
            THEN
          RETURN
        END IF



        phiold = 1.0
        gamana = 1.0
        gamas1 = 1.0
        gamas2 = 1.0
        gamaab = 1.0
        gamold = 1.0



        mnh4 = tnh4/wh2o



        ynh4 = tnh4


        DO nnn = 1, 150
          nitr = nnn




          rk2sa = k2sa*gamas2*gamas2/(gamas1*gamas1*gamas1)
          rkna = kna/(gamana*gamana)
          rknwet = rkna*wh2o
          t21 = zso4 - mnh4
          t221 = zso4 + t21



          a2 = rk2sa + rknwet - t21
          a1 = rk2sa*rknwet - t21*(rk2sa+rknwet) - rk2sa*zso4 - rkna*tno3
          a0 = -(t21*rk2sa*rknwet+rk2sa*rknwet*zso4+rk2sa*rkna*tno3)


          CALL cubic(a2,a1,a0,nr,crutes)



          hplus = crutes(1)
          bal = hplus**3 + a2*hplus**2 + a1*hplus + a0
          mso4 = rk2sa*zso4/(hplus+rk2sa) 
          mhso4 = zso4 - & 
            mso4
          mna = rkna*tno3/(hplus+rknwet) 
          mna = max(0.0,mna)
          mna = min(mna,tno3/wh2o)
          xno3 = mna*wh2o
          ano3 = mna*wh2o*mwno3
          gno3 = (tno3-xno3)*mwhno3


          stion = 0.5*(hplus+mna+mnh4+mhso4+4.0*mso4)



          CALL awater(irh,tso4,ynh4,xno3,ah2o)




          wh2o = 1.0E-3*ah2o
          cat(1) = max(hplus,0.0)
          cat(2) = max(mnh4,0.0)
          an(1) = max(mso4,0.0)
          an(2) = max(mna,0.0)
          an(3) = max(mhso4,0.0)

          CALL actcof(cat,an,gams,molnu,phibar)

          gamana = gams(1,2)
          gamas1 = gams(1,1)
          gamas2 = gams(1,3)
          gamaan = gams(2,2)

          gamahat = (gamas2*gamas2/(gamaab*gamaab))
          bhat = khat*gamahat


          eror = abs(gamold-gamahat)/gamold
          gamold = gamahat




          IF (eror<=toler2) THEN



            RETURN
          END IF

        END DO



        gno3 = tno3*mwhno3
        ano3 = 0.0
        CALL awater(irh,tso4,tnh4,tno3,ah2o)
        RETURN


      END IF


    END SUBROUTINE rpmares_old











































    SUBROUTINE rpmmod3(nspcsda,blksize,layer,dtsec,pres,temp,relhum, &
        nitrate_in,nh3_in,vsulf_in,so4rat_in,drog_in,ldrog,condvap_in,ncv, &
        nacv,eeci_in,eecj_in,eorgi_in,eorgj_in,epm25i,epm25j,epmcoarse,    &
        soilrat_in,cblk,igrid,jgrid,kgrid)









      INTEGER blksize

      INTEGER numcells


      PARAMETER (numcells=1)


      INTEGER layer


      INTEGER ncell

      PARAMETER (ncell=1)


      REAL temp

      REAL relhum

      REAL pres

      REAL numnuc_in

      REAL numacc_in

      REAL numcor_in
                         
      REAL vsulf_in

                         
      REAL asulf_in


      REAL asulfi_in

      REAL nh3_in

      REAL nitrate_in

      REAL so4rat_in
                         
      REAL soilrat_in


      REAL eeci_in

      REAL eecj_in

      REAL eorgi_in

      REAL eorgj_in



      INTEGER ncv

      INTEGER nacv

      INTEGER ldrog
      REAL drog_in(ldrog)                                 

      REAL condvap_in(ncv) 
      REAL drog(blksize,ldrog)                                 





      REAL epm25i(blksize) 
      REAL epm25j(blksize) 


      REAL eorgi(blksize) 
      REAL eorgj(blksize) 


      REAL eeci(blksize) 
      REAL eecj(blksize) 



      REAL epm25(blksize) 
      REAL esoil(blksize) 
      REAL eseas(blksize) 
      REAL epmcoarse(blksize) 


      REAL dtsec



      REAL newm3

      REAL totaersulf


      INTEGER numsteps

      REAL step





      INTEGER nspcsda


      REAL cblk(blksize,nspcsda) 






      REAL blkta(blksize) 
      REAL blkprs(blksize) 
      REAL blkdens(blksize) 
      REAL blkrh(blksize) 





      REAL so4rat(blksize)                                 

      REAL orgaro1rat(blksize)                                 

      REAL orgaro2rat(blksize)                                 

      REAL orgalk1rat(blksize)                                 

      REAL orgole1rat(blksize)                                 

      REAL orgbio1rat(blksize)                                 

      REAL orgbio2rat(blksize)                                 

      REAL orgbio3rat(blksize)                                 

      REAL orgbio4rat(blksize)                                 




      REAL xlm(blksize) 
      REAL amu(blksize) 






      REAL dgnuc(blksize) 
      REAL dgacc(blksize) 
      REAL dgcor(blksize) 




      REAL pmassn(blksize) 
      REAL pmassa(blksize) 
      REAL pmassc(blksize) 



      REAL pdensn(blksize) 
      REAL pdensa(blksize) 
      REAL pdensc(blksize) 



      REAL knnuc(blksize) 
      REAL knacc(blksize) 
      REAL kncor(blksize) 



      REAL fconcn(blksize) 

      REAL fconca(blksize) 

      REAL fconcn_org(blksize)
      REAL fconca_org(blksize)





      REAL dmdt(blksize) 




      REAL dndt(blksize) 




      REAL cgrn3(blksize) 
      REAL cgra3(blksize) 





      REAL urn00(blksize) 
      REAL ura00(blksize) 




      REAL brna01(blksize) 
      REAL brna31(blksize) 



      REAL deltaso4a(blksize) 





      INTEGER unit
      PARAMETER (unit=30)

      CHARACTER*16 pname
      PARAMETER (pname=' BOX            ')




      INTEGER isp,igrid,jgrid,kgrid


      INTEGER ii, iimap(8)
      DATA iimap/1, 2, 18, 19, 21, 22, 23, 24/

















      step = & 
        dtsec
      blkta(blksize) = & 
        temp
      blkprs(blksize) = pres* & 
        100.
      blkrh(blksize) = & 
        relhum
      blkdens(blksize) = blkprs(blksize)/(rdgas*blkta(blksize)) 












      DO isp = 1, ldrog
        drog(blksize,isp) = drog_in(isp)
      END DO













      so4rat(blksize) = so4rat_in








      eorgi(blksize) = & 
        eorgi_in
      eorgj(blksize) = & 
        eorgj_in

      eeci(blksize) = & 
        eeci_in
      eecj(blksize) = & 
        eecj_in

      epm25(blksize) = & 
        0.0
      esoil(blksize) = & 
        soilrat_in
      eseas(blksize) = & 
        0.0


      dgnuc(blksize) = dginin
      dgacc(blksize) = dginia
      dgcor(blksize) = dginic
      newm3 = 0.0





      totaersulf = 0.0
      newm3 = 0.0








      CALL aeroproc(blksize,nspcsda,numcells,layer,cblk,step,blkta,blkprs, &
        blkdens,blkrh,so4rat,orgaro1rat,orgaro2rat,orgalk1rat, &
        orgole1rat,orgbio1rat,orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog,ncv, &
        nacv,epm25i,epm25j,eorgi,eorgj,eeci,eecj,epmcoarse,esoil,eseas,xlm, &
        amu,dgnuc,dgacc,dgcor,pmassn,pmassa,pmassc,pdensn,pdensa,pdensc,knnuc, &
        knacc,kncor,fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,cgrn3,cgra3, &
        urn00,ura00,brna01,brna31,deltaso4a,igrid,jgrid,kgrid)

















      RETURN


    END SUBROUTINE rpmmod3

    SUBROUTINE soa_part(layer,blkta,blkprs,orgaro1rat,orgaro2rat,orgalk1rat, &
        orgole1rat,orgbio1rat,orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog,ncv, &
        nacv,cblk,blksize,nspcsda,numcells,dt)






















































































      INTEGER layer

      INTEGER blksize

      INTEGER nspcsda

      INTEGER numcells

      INTEGER ldrog

      INTEGER ncv

      INTEGER nacv
      REAL cblk(blksize,nspcsda) 

      REAL dt
      REAL blkta(blksize) 
      REAL blkprs(blksize) 
      REAL orgaro1rat(blksize)                                       

      REAL orgaro2rat(blksize)                                       

      REAL orgalk1rat(blksize)                                       

      REAL orgole1rat(blksize)                                       

      REAL orgbio1rat(blksize) 
      REAL orgbio2rat(blksize) 
      REAL orgbio3rat(blksize) 
      REAL orgbio4rat(blksize) 
      REAL drog(blksize,ldrog) 




      REAL thrsmin
      PARAMETER (thrsmin=1.E-19)



      REAL rgas
      PARAMETER (rgas=8.314510)

      REAL tnull
      PARAMETER (tnull=298.)

      REAL mwc
      PARAMETER (mwc=12.0)

      REAL mworg
      PARAMETER (mworg=175.0)

      REAL mwso4
      PARAMETER (mwso4=96.0576)

      REAL mwnh4
      PARAMETER (mwnh4=18.03858)

      REAL mwno3
      PARAMETER (mwno3=62.01287)

      REAL rtol
      PARAMETER (rtol=1.E-04)




      INTEGER lcell
      INTEGER l, & 
        n

      REAL convfac

      REAL ttinv

      REAL minitw

      REAL mtotw

      REAL mnonow

      REAL imtotw

      REAL minit

      REAL mnono

      REAL mtot

      REAL thres

      REAL mcheck
      REAL msum(ncv) 
      REAL mwcv(ncv) 
      REAL imwcv(ncv) 
      REAL pnull(ncv) 
      REAL dhvap(ncv) 
      REAL pvap(ncv) 
      REAL ctot(ncv) 
      REAL cgas(ncv) 
      REAL caer(ncv) 
      REAL asav(ncv) 
      REAL aold(ncv) 
      REAL csat(ncv) 
      REAL alpha(ncv) 
      REAL prod(ncv) 
      REAL p(ncv) 
      REAL f(ldrog) 

      LOGICAL check

      INTEGER its








      dhvap(psoaaro1) = 156.0E03
      dhvap(psoaaro2) = 156.0E03
      dhvap(psoaalk1) = 156.0E03
      dhvap(psoaole1) = 156.0E03
      dhvap(psoaapi1) = 156.0E03
      dhvap(psoaapi2) = 156.0E03
      dhvap(psoalim1) = 156.0E03
      dhvap(psoalim2) = 156.0E03















      mwcv(psoaaro1) = 150.
      mwcv(psoaaro2) = 150.
      mwcv(psoaalk1) = 140.
      mwcv(psoaole1) = 140.
      mwcv(psoaapi1) = 184.
      mwcv(psoaapi2) = 184.
      mwcv(psoalim1) = 200.
      mwcv(psoalim2) = 200.











































      alpha(psoaaro1) = 0.039
      alpha(psoaaro2) = 0.108
      alpha(psoaalk1) = 0.048
      alpha(psoaole1) = 0.008


      alpha(psoaapi1) = & 
        0.092
      alpha(psoaapi2) = & 
        0.075
      alpha(psoalim1) = 0.163
      alpha(psoalim2) = 0.247

















      pnull(psoaaro1) = 5.7E-05
      pnull(psoaaro2) = 1.6E-03
      pnull(psoaalk1) = 5.0E-06
      pnull(psoaole1) = 5.0E-06


      pnull(psoaapi1) = & 
        2.488E-05
      pnull(psoaapi2) = & 
        2.778E-05
      pnull(psoalim1) = 2.5E-05
      pnull(psoalim2) = 1.2E-04



      f(pxyl) = & 
        1.
      f(ptol) = & 
        1.
      f(pcsl1) = & 
        1.
      f(pcsl2) = & 
        1.
      f(phc8) = & 
        1.
      f(poli1) = & 
        1.
      f(poli2) = & 
        1.
      f(poli3) = & 
        1.
      f(polt1) = & 
        1.
      f(polt2) = & 
        1.
      f(polt3) = & 
        1.



      f(papi1) = & 
        0.
      f(papi2) = & 
        0.
      f(papi3) = & 
        1.
      f(plim1) = & 
        0.228
      f(plim2) = & 
        0.
      f(plim3) = & 
        0.771



      DO lcell = 1, numcells
        DO l = 1, ldrog
          drog(lcell,l) = f(l)*drog(lcell,l)
        END DO
        ttinv = 1./tnull - 1./blkta(lcell)
        convfac = blkprs(lcell)/(rgas*blkta(lcell))
        cgas(psoaaro1) = cblk(lcell,vcvaro1)
        cgas(psoaaro2) = cblk(lcell,vcvaro2)
        cgas(psoaalk1) = cblk(lcell,vcvalk1)
        cgas(psoaole1) = cblk(lcell,vcvole1)
        cgas(psoaapi1) = cblk(lcell,vcvapi1)
        cgas(psoaapi2) = cblk(lcell,vcvapi2)
        cgas(psoalim1) = cblk(lcell,vcvlim1)
        cgas(psoalim2) = cblk(lcell,vcvlim2)
        caer(psoaaro1) = cblk(lcell,vorgaro1j) + cblk(lcell,vorgaro1i)
        caer(psoaaro2) = cblk(lcell,vorgaro2j) + cblk(lcell,vorgaro2i)
        caer(psoaalk1) = cblk(lcell,vorgalk1j) + cblk(lcell,vorgalk1i)
        caer(psoaole1) = cblk(lcell,vorgole1j) + cblk(lcell,vorgole1i)
        caer(psoaapi1) = cblk(lcell,vorgba1j) + cblk(lcell,vorgba1i)
        caer(psoaapi2) = cblk(lcell,vorgba2j) + cblk(lcell,vorgba2i)
        caer(psoalim1) = cblk(lcell,vorgba3j) + cblk(lcell,vorgba3i)
        caer(psoalim2) = cblk(lcell,vorgba4j) + cblk(lcell,vorgba4i)

        prod(psoaaro1) = drog(lcell,pxyl) + drog(lcell,ptol) + &
          drog(lcell,pcsl1) + drog(lcell,pcsl2)
        prod(psoaaro2) = drog(lcell,pxyl) + drog(lcell,ptol) + &
          drog(lcell,pcsl1) + drog(lcell,pcsl2)
        prod(psoaalk1) = drog(lcell,phc8)
        prod(psoaole1) = drog(lcell,poli1) + drog(lcell,poli2) + &

          drog(lcell,poli3) + drog(lcell,polt1) + drog(lcell,polt2) + &
          drog(lcell,polt3)
        prod(psoaapi1) = drog(lcell,papi1) + drog(lcell,papi2) + &
          drog(lcell,papi3)
        prod(psoaapi2) = drog(lcell,papi1) + drog(lcell,papi2) + &
          drog(lcell,papi3)
        prod(psoalim1) = drog(lcell,plim1) + drog(lcell,plim2) + &
          drog(lcell,plim3)
        prod(psoalim2) = drog(lcell,plim1) + drog(lcell,plim2) + &
          drog(lcell,plim3)








        thres = 0.
        mtot = 0.
        mtotw = 0.
        DO l = 1, ncv
          prod(l) = convfac*mwcv(l)*alpha(l)*prod(l)
          ctot(l) = prod(l) + cgas(l) + caer(l) 
          p(l) = prod(l)
          msum(l) = cgas(l) + caer(l) + prod(l)
          aold(l) = caer(l)
          imwcv(l) = 1./mwcv(l)
          pvap(l) = pnull(l)*exp(dhvap(l)/rgas*ttinv)
          csat(l) = pvap(l)*mwcv(l)*1.0E06/(rgas*blkta(lcell))
          thres = thres + ((cgas(l)+prod(l))/csat(l))
          mtot = mtot + caer(l)
          mtotw = mtotw + caer(l)*imwcv(l)
        END DO





        mnono = 0.0001*(cblk(lcell,vso4aj)+cblk(lcell,vnh4aj)+cblk(lcell, &
          vno3aj))
        mnono = mnono + 0.0001*(cblk(lcell,vso4ai)+cblk(lcell,vnh4ai)+cblk( &
          lcell,vno3ai))
        mnonow = 0.0001*(cblk(lcell,vso4aj)/mwso4+cblk(lcell,vnh4aj)/mwnh4+ &
          cblk(lcell,vno3aj)/mwno3)
        mnonow = mnonow + 0.0001*(cblk(lcell,vso4ai)/mwso4+cblk(lcell,vnh4ai)/ &
          mwnh4+cblk(lcell,vno3ai)/mwno3)
        mnono = max(mnono,conmin)
        mnonow = max(mnonow,conmin)




        minit = cblk(lcell,vecj) + cblk(lcell,veci) + cblk(lcell,vorgpaj) + &
          cblk(lcell,vorgpai) + mnono
        minitw = (cblk(lcell,vecj)+cblk(lcell,veci))/mwc + &
          (cblk(lcell,vorgpaj)+cblk(lcell,vorgpai))/mworg + mnonow






        minit = 0.
        minitw = 0.

        mtot = mtot + minit
        mtotw = mtotw + minitw
        imtotw = 1./mtotw



        IF ((thres>1 .AND. minitw<thrsmin) .OR. (minitw>thrsmin) .OR. &
            (mtot>thrsmin)) THEN

          DO l = 1, ncv
            ctot(l) = p(l) + cgas(l) + caer(l)
            caer(l) = ctot(l) 
          END DO




          CALL newt(layer,caer,ncv,check,ctot,csat,imwcv,minitw,its)

          IF (check) THEN

          END IF



          DO l = 1, ncv
            IF (caer(l)<=tolmin) THEN

              caer(l) = conmin
            END IF
            IF (caer(l)>ctot(l)) THEN
              IF (caer(l)-ctot(l)>tolmin) THEN

              END IF
              caer(l) = ctot(l)
            END IF
            cgas(l) = ctot(l) - caer(l)
          END DO






          cblk(lcell,vcvaro1) = max(cgas(psoaaro1),conmin)
          cblk(lcell,vcvaro2) = max(cgas(psoaaro2),conmin)
          cblk(lcell,vcvalk1) = max(cgas(psoaalk1),conmin)
          cblk(lcell,vcvole1) = max(cgas(psoaole1),conmin)
          cblk(lcell,vcvapi1) = max(cgas(psoaapi1),conmin)
          cblk(lcell,vcvapi2) = max(cgas(psoaapi2),conmin)
          cblk(lcell,vcvlim1) = max(cgas(psoalim1),conmin)
          cblk(lcell,vcvlim2) = max(cgas(psoalim2),conmin)
          orgaro1rat(lcell) = (caer(psoaaro1)-aold(psoaaro1))/dt
          orgaro2rat(lcell) = (caer(psoaaro2)-aold(psoaaro2))/dt
          orgalk1rat(lcell) = (caer(psoaalk1)-aold(psoaalk1))/dt
          orgole1rat(lcell) = (caer(psoaole1)-aold(psoaole1))/dt
          orgbio1rat(lcell) = (caer(psoaapi1)-aold(psoaapi1))/dt
          orgbio2rat(lcell) = (caer(psoaapi2)-aold(psoaapi2))/dt
          orgbio3rat(lcell) = (caer(psoalim1)-aold(psoalim1))/dt
          orgbio4rat(lcell) = (caer(psoalim2)-aold(psoalim2))/dt


        ELSE





          DO l = 1, ncv
            caer(l) = ctot(l) - csat(l)
            caer(l) = max(caer(l),0.)
            cgas(l) = ctot(l) - caer(l)
          END DO

          cblk(lcell,vcvaro1) = cgas(psoaaro1)
          cblk(lcell,vcvaro2) = cgas(psoaaro2)
          cblk(lcell,vcvalk1) = cgas(psoaalk1)
          cblk(lcell,vcvole1) = cgas(psoaole1)
          cblk(lcell,vcvapi1) = cgas(psoaapi1)
          cblk(lcell,vcvapi2) = cgas(psoaapi2)
          cblk(lcell,vcvlim1) = cgas(psoalim1)
          cblk(lcell,vcvlim2) = cgas(psoalim2)
          orgaro1rat(lcell) = (caer(psoaaro1)-aold(psoaaro1))/dt
          orgaro2rat(lcell) = (caer(psoaaro2)-aold(psoaaro2))/dt
          orgalk1rat(lcell) = (caer(psoaalk1)-aold(psoaalk1))/dt
          orgole1rat(lcell) = (caer(psoaole1)-aold(psoaole1))/dt
          orgbio1rat(lcell) = (caer(psoaapi1)-aold(psoaapi1))/dt
          orgbio2rat(lcell) = (caer(psoaapi2)-aold(psoaapi2))/dt
          orgbio3rat(lcell) = (caer(psoalim1)-aold(psoalim1))/dt
          orgbio4rat(lcell) = (caer(psoalim2)-aold(psoalim2))/dt

        END IF



        DO l = 1, ncv

          IF (cgas(l)==0. .AND. caer(l)==0. .AND. msum(l)==0) THEN
            mcheck = 1.
          ELSE
            mcheck = (cgas(l)+caer(l))/msum(l)
          END IF
          IF ((mcheck<1.-rtol) .OR. (mcheck>1.+rtol)) THEN



90020       FORMAT ('LAYER = ',I2,', L = ',I2,', MCHECK = ',E12.6,', MASS = ', &
              E12.6)
          END IF
        END DO


      END DO



      RETURN
    END SUBROUTINE soa_part
    SUBROUTINE sorgam(layer,blkta,blkprs,orgaro1rat,orgaro2rat,orgalk1rat, &
        orgole1rat,orgbio1rat,orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog,ncv, &
        nacv,cblk,blksize,nspcsda,numcells,dt)



































































      INTEGER blksize

      INTEGER nspcsda

      INTEGER numcells

      INTEGER layer

      INTEGER ldrog

      INTEGER ncv

      INTEGER nacv

      REAL dt
      REAL cblk(blksize,nspcsda) 
      REAL blkta(blksize) 
      REAL blkprs(blksize) 
      REAL orgaro1rat(blksize)                                       

      REAL orgaro2rat(blksize)                                       

      REAL orgalk1rat(blksize)                                       

      REAL orgole1rat(blksize)                                       

      REAL orgbio1rat(blksize) 
      REAL orgbio2rat(blksize) 
      REAL orgbio3rat(blksize) 
      REAL orgbio4rat(blksize) 
      REAL drog(blksize,ldrog) 



















      IF (orgaer==1) THEN



















      ELSE IF (orgaer==2) THEN








        CALL soa_part(layer,blkta,blkprs,orgaro1rat,orgaro2rat,orgalk1rat, &
          orgole1rat,orgbio1rat,orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog, &
          ncv,nacv,cblk,blksize,nspcsda,numcells,dt)
      ELSE






      END IF













90000 FORMAT ('ORGAER = ',I2)



      RETURN
    END SUBROUTINE sorgam













 
       SUBROUTINE VDVG(  BLKSIZE, NSPCSDA, NUMCELLS,           &
                     LAYER,                                    &
                     CBLK,                                     &  
                     BLKTA, BLKDENS, RA, USTAR, WSTAR,  AMU,   &
                     DGNUC, DGACC, DGCOR,                      &
                     KNNUC, KNACC,KNCOR,                       &    
                     PDENSN, PDENSA, PDENSC,                   &                 
                     VSED, VDEP )







      INTEGER BLKSIZE                  
      INTEGER NSPCSDA                  
      INTEGER NUMCELLS                
      INTEGER LAYER                   

      REAL CBLK( BLKSIZE, NSPCSDA ) 
      REAL BLKTA( BLKSIZE )         
      REAL BLKDENS(BLKSIZE) 
      REAL RA(BLKSIZE )             
      REAL USTAR( BLKSIZE )         
      REAL WSTAR( BLKSIZE )         
      REAL AMU( BLKSIZE )           
      REAL DGNUC( BLKSIZE )         
      REAL DGACC( BLKSIZE )         
      REAL DGCOR( BLKSIZE )         
      REAL KNNUC( BLKSIZE )         
      REAL KNACC( BLKSIZE )         
      REAL KNCOR( BLKSIZE )         
      REAL PDENSN( BLKSIZE )        
      REAL PDENSA( BLKSIZE )        
      REAL PDENSC( BLKSIZE )        
       



      REAL DCHAT0N( BLKSIZE), DCHAT0A(BLKSIZE), DCHAT0C(BLKSIZE)
      REAL DCHAT3N( BLKSIZE), DCHAT3A(BLKSIZE), DCHAT3C(BLKSIZE)


      
      REAL VGHAT0N( BLKSIZE), VGHAT0A(BLKSIZE), VGHAT0C(BLKSIZE)
      REAL VGHAT3N( BLKSIZE), VGHAT3A(BLKSIZE), VGHAT3C(BLKSIZE)



      REAL VDEP( BLKSIZE, NASPCSDEP) 
      REAL VSED( BLKSIZE, NASPCSSED)  
      
      
      INTEGER LCELL
      REAL DCONST1, DCONST1N, DCONST1A, DCONST1C
      REAL DCONST2, DCONST3N, DCONST3A,DCONST3C 
      REAL SC0N, SC0A, SC0C 
      REAL SC3N, SC3A, SC3C 
      REAL ST0N, ST0A, ST0C 
      REAL ST3N, ST3A, ST3C 
      REAL RD0N, RD0A, RD0C    
      REAL RD3N, RD3A, RD3C    
      REAL UTSCALE   
      REAL NU        
      REAL USTFAC      
      REAL BHAT
      PARAMETER( BHAT =  1.246 ) 




         IF ( LAYER .EQ. 1 ) THEN 

	        
         DO LCELL = 1, NUMCELLS
         
            DCONST1 = BOLTZ * BLKTA(LCELL) /                                         &
                    ( THREEPI * AMU(LCELL) )
            DCONST1N = DCONST1 / DGNUC( LCELL ) 
            DCONST1A = DCONST1 / DGACC( LCELL )
            DCONST1C = DCONST1 / DGCOR( LCELL )   
            DCONST2 = GRAV / ( 18.0 * AMU(LCELL) )
            DCONST3N = DCONST2 * PDENSN(LCELL) * DGNUC( LCELL )**2
            DCONST3A = DCONST2 * PDENSA(LCELL) * DGACC( LCELL )**2
            DCONST3C = DCONST2 * PDENSC(LCELL) * DGCOR( LCELL )**2


 
            DCHAT0N(LCELL) =  DCONST1N                             &
               * ( ESN04 + BHAT * KNNUC( LCELL ) * ESN16 )
                
            DCHAT3N(LCELL) =  DCONST1N                             &
               * ( ESNM20 + BHAT * KNNUC( LCELL ) * ESNM32 )
            
            VGHAT0N(LCELL) = DCONST3N                             &
               * ( ESN16 + BHAT * KNNUC( LCELL ) * ESN04 )
                
            VGHAT3N(LCELL) = DCONST3N                             &
               * (ESN64 + BHAT * KNNUC( LCELL ) * ESN28 )



            DCHAT0A(LCELL) =  DCONST1A                             &
              * ( ESA04 + BHAT * KNACC( LCELL ) * ESA16 )
                
            DCHAT3A(LCELL) =  DCONST1A                             &
               * ( ESAM20 + BHAT * KNACC( LCELL ) * ESAM32 )           
            
            VGHAT0A(LCELL) = DCONST3A                             &
              * ( ESA16 + BHAT * KNACC( LCELL ) * ESA04 )
                
            VGHAT3A(LCELL) = DCONST3A                             &
              * ( ESA64 + BHAT * KNACC( LCELL ) * ESA28 )




            DCHAT0C(LCELL)=  DCONST1C                             &
              * ( ESC04 + BHAT * KNCOR( LCELL ) * ESC16 )
                
            DCHAT3C(LCELL) = DCONST1C                             &
              * ( ESCM20 + BHAT * KNCOR( LCELL ) * ESCM32 )
            
            VGHAT0C(LCELL) = DCONST3C                             &
              * ( ESC16 + BHAT * KNCOR( LCELL ) * ESC04 )
                
            VGHAT3C(LCELL) = DCONST3C                             &
              * ( ESC64 + BHAT * KNCOR( LCELL ) * ESC28 )
        
        END DO
 









        DO LCELL = 1, NUMCELLS
        
         NU = AMU(LCELL) / BLKDENS(LCELL) 
         USTFAC = USTAR(LCELL) * USTAR(LCELL) / ( GRAV * NU)
         UTSCALE = USTAR(LCELL) +                             &
                 0.24 * WSTAR(LCELL) * WSTAR(LCELL) / USTAR(LCELL)


           


        SC0N = NU / DCHAT0N(LCELL)      
        ST0N = MAX( VGHAT0N(LCELL) * USTFAC , 0.01)
        RD0N = 1.0 / ( UTSCALE *                             &
                  ( SC0N**(-TWO3) + 10.0**(-3.0 / ST0N) ) ) 
      
        VDEP(LCELL, VDNNUC) = VGHAT0N(LCELL) +                             &
               1.0 / (                             &
           RA(LCELL) + RD0N + RD0N * RA(LCELL) * VGHAT0N(LCELL) )

        VSED( LCELL, VSNNUC) = VGHAT0N(LCELL) 
     


        SC0A = NU / DCHAT0A(LCELL)      
        ST0A = MAX ( VGHAT0A(LCELL) * USTFAC, 0.01)
        RD0A = 1.0 / ( UTSCALE *                             &
                  ( SC0A**(-TWO3) + 10.0**(-3.0 / ST0A) ) ) 
      
        VDEP(LCELL, VDNACC) = VGHAT0A(LCELL) +                             &
               1.0 / (                             &
           RA(LCELL) + RD0A + RD0A * RA(LCELL) * VGHAT0A(LCELL) ) 

        VSED( LCELL, VSNACC) = VGHAT0A(LCELL) 



        SC0C = NU / DCHAT0C(LCELL)      



 
         RD0C = 1.0 / ( UTSCALE *                            &
                      ( SC0C ** ( -TWO3 )  ) ) 
      
        VDEP(LCELL, VDNCOR) = VGHAT0C(LCELL) +                             &
               1.0 / (                             &
           RA(LCELL) + RD0C + RD0C * RA(LCELL) * VGHAT0C(LCELL) ) 

        VSED( LCELL, VSNCOR) = VGHAT0C(LCELL) 





        SC3N = NU / DCHAT3N(LCELL)      
        ST3N = MAX( VGHAT3N(LCELL) * USTFAC, 0.01) 
        RD3N = 1.0 / ( UTSCALE *                             &
                  ( SC3N**(-TWO3) + 10.0**(-3.0 / ST3N) ) ) 
      
        VDEP(LCELL, VDMNUC) = VGHAT3N(LCELL) +                             &
               1.0 / (                             &
           RA(LCELL) + RD3N + RD3N * RA(LCELL) * VGHAT3N(LCELL) ) 

        VSED(LCELL, VSMNUC) = VGHAT3N(LCELL)
     


        SC3A = NU / DCHAT3A(LCELL)      
        ST3A = MAX( VGHAT3A(LCELL) * USTFAC , 0.01 )
        RD3A = 1.0 / ( UTSCALE *                             &
                  ( SC3A**(-TWO3) + 10.0**(-3.0 / ST3A) ) ) 

       VDEP(LCELL, VDMACC) = VGHAT3A(LCELL) +                            &
               1.0 / (                            &
               RA(LCELL) + RD3A + RD3A * RA(LCELL) * VGHAT3A(LCELL) )
                
     









     
 







        VSED( LCELL, VSMACC ) = VGHAT3A(LCELL)



        SC3C = NU / DCHAT3C(LCELL)



   
        RD3C = 1.0 / ( UTSCALE *                            &
                     ( SC3C ** ( -TWO3 ) ) ) 
        VDEP(LCELL, VDMCOR) = VGHAT3C(LCELL) +                             &
               1.0 / (                             &
           RA(LCELL) + RD3C + RD3C * RA(LCELL) * VGHAT3C(LCELL)) 



        VSED( LCELL, VSMCOR) = VGHAT3C(LCELL) 


                                 
        END DO  
             
        ELSE   
        


         DO LCELL = 1, NUMCELLS
         
            DCONST2 = GRAV / ( 18.0 * AMU(LCELL) )
            
            DCONST3N = DCONST2 * PDENSN(LCELL) * DGNUC( LCELL )**2
            DCONST3A = DCONST2 * PDENSA(LCELL) * DGACC( LCELL )**2
            DCONST3C = DCONST2 * PDENSC(LCELL) * DGCOR( LCELL )**2

            VGHAT0N(LCELL) = DCONST3N                             &
               * ( ESN16 + BHAT * KNNUC( LCELL ) * ESN04 )
               


            VSED( LCELL, VSNNUC) = VGHAT0N(LCELL)
 
            VGHAT3N(LCELL) = DCONST3N                             &
               * (ESN64 + BHAT * KNNUC( LCELL ) * ESN28 )



	    VSED( LCELL, VSMNUC) = VGHAT3N(LCELL)

            VGHAT0A(LCELL) = DCONST3A                             &
              * ( ESA16 + BHAT * KNACC( LCELL ) * ESA04 )


     
            VSED( LCELL, VSNACC) = VGHAT0A(LCELL)      
                
            VGHAT3A(LCELL) = DCONST3A                            & 
              * ( ESA64 + BHAT * KNACC( LCELL ) * ESA28 )
     







            VSED( LCELL, VSMACC) = VGHAT3A(LCELL)     
         
            VGHAT0C(LCELL) = DCONST3C                            & 
              * ( ESC16 + BHAT * KNCOR( LCELL ) * ESC04 )


     
            VSED( LCELL, VSNCOR) = VGHAT0C(LCELL) 
       
                
            VGHAT3C(LCELL) = DCONST3C                             &
              * ( ESC64 + BHAT * KNCOR( LCELL ) * ESC28 )



            VSED( LCELL, VSMCOR) = VGHAT3C(LCELL) 
        
         END DO 
         
         END IF 
         
END SUBROUTINE VDVG


























       SUBROUTINE VDVG_2(  BLKSIZE, NSPCSDA, NUMCELLS,           &
                     LAYER,                                    &
                     CBLK, BLKTA, BLKDENS,                     &
                     RA, USTAR, PBLH, ZNTT, RMOLM,  AMU,       &
                     DGNUC, DGACC, DGCOR, XLM,                 &
                     KNNUC, KNACC,KNCOR,                       &
                     PDENSN, PDENSA, PDENSC,                   &
                     VSED, VDEP)







      INTEGER BLKSIZE                  
      INTEGER NSPCSDA                  
      INTEGER NUMCELLS                
      INTEGER LAYER                   
      INTEGER, PARAMETER :: iprnt = 0

      REAL CBLK( BLKSIZE, NSPCSDA ) 
      REAL BLKTA( BLKSIZE )         
      REAL BLKDENS(BLKSIZE) 
      REAL RA(BLKSIZE )             
      REAL USTAR( BLKSIZE )         
      REAL PBLH( BLKSIZE )          
      REAL ZNTT( BLKSIZE )          
      REAL RMOLM( BLKSIZE )         
      REAL AMU( BLKSIZE )           
      REAL XLM( BLKSIZE )           
      REAL DGNUC( BLKSIZE )         
      REAL DGACC( BLKSIZE )         
      REAL DGCOR( BLKSIZE )         
      REAL KNNUC( BLKSIZE )         
      REAL KNACC( BLKSIZE )         
      REAL KNCOR( BLKSIZE )         
      REAL PDENSN( BLKSIZE )        
      REAL PDENSA( BLKSIZE )        
      REAL PDENSC( BLKSIZE )        




      REAL VDEP( BLKSIZE, NASPCSDEP) 
      REAL VSED( BLKSIZE, NASPCSSED)  


      INTEGER LCELL,N
      REAL DCONST1, DCONST2, DCONST3, DCONST3N, DCONST3A,DCONST3C
      REAL UTSCALE,CZH   
      REAL NU        
      REAL BHAT
      PARAMETER( BHAT =  1.246 ) 
      REAL COLCTR_BIGD,COLCTR_SMALD
      PARAMETER( COLCTR_BIGD=2.E-3,COLCTR_SMALD=20.E-6)  
      REAL SUM0, SUM3, DQ, KNQ, CUNQ, VSEDQ, SCQ, STQ, RSURFQ, vdplim
      REAL Eff_dif, Eff_imp, Eff_int, RBcor
      INTEGER ISTOPvd0,IdoWesCor
      PARAMETER (ISTOPvd0 = 0)  
      PARAMETER (IdoWesCor = 0) 
      IF (ISTOPvd0.EQ.1)THEN
      RETURN
      ENDIF


      IF(iprnt.eq.1) print *,'In VDVG, Layer=',LAYER
         IF ( LAYER .EQ. 1 ) THEN 

                 
         DO LCELL = 1, NUMCELLS     
         
            DCONST1 = BOLTZ * BLKTA(LCELL) /                                         &
                    ( THREEPI * AMU(LCELL) )
            DCONST2 = GRAV / ( 18.0 * AMU(LCELL) )
            DCONST3 =  USTAR(LCELL)/(9.*AMU(LCELL)*COLCTR_BIGD)
 


         NU = AMU(LCELL) / BLKDENS(LCELL) 

         UTSCALE =  1.
        IF (IdoWesCor.EQ.1)THEN

         IF(RMOLM(LCELL).LT.0.)THEN
         CZH = -1.*PBLH(LCELL)*RMOLM(LCELL)
         IF(CZH.GT.30.0)THEN
         UTSCALE=0.45*CZH**0.6667
         ELSE
         UTSCALE=1.+(-300.*RMOLM(LCELL))**0.6667
         ENDIF
         ENDIF
        ENDIF   
         UTSCALE = USTAR(LCELL)*UTSCALE
      IF(iprnt.eq.1)THEN
      print *,'NGAUSdv,xxlsga,USTAR,UTSCALE'
      print *,NGAUSdv,xxlsga,USTAR(LCELL),UTSCALE
      print *,'DCONST2,PDENSA,DGACC,GRAV,AMU'
      print *,DCONST2,PDENSA(LCELL),DGACC(LCELL),GRAV,AMU(LCELL)
      endif
      

      
        SUM0=0.
        SUM3=0.
        DO N=1,NGAUSdv
        DQ=DGNUC(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsgn)  
        KNQ=2.*XLM(LCELL)/DQ  
        CUNQ=1.+KNQ*(1.257+.4*exp(-1.1/KNQ))  
        VSEDQ=PDENSN(LCELL)*DCONST2*CUNQ*DQ*DQ  
        SCQ=NU*DQ/DCONST1/CUNQ  
        Eff_dif=SCQ**(-TWO3)    
        STQ=DCONST3*PDENSN(LCELL)*DQ**2  
        Eff_imp=(STQ/(0.8+STQ))**2   

        Eff_int=(0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 
        RBcor=1. 
        vdplim=UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor

        vdplim=min(vdplim,.02)
        RSURFQ=RA(LCELL)+1./vdplim





        SUM0=SUM0+WGAUS(N)*(VSEDQ + 1./RSURFQ)  
        SUM3=SUM3+WGAUS(N)*(VSEDQ + 1./RSURFQ)*DQ**3  
        ENDDO
        VDEP(LCELL, VDNNUC) = SUM0/sqrtpi  
        VDEP(LCELL, VDMNUC) = SUM3/(sqrtpi*EXP((1.5*sqrt2*xxlsgn)**2)*DGNUC(LCELL)**3) 



        SUM0=0.
        SUM3=0.
        DO N=1,NGAUSdv
        DQ=DGACC(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsga)  
        KNQ=2.*XLM(LCELL)/DQ  
        CUNQ=1.+KNQ*(1.257+.4*exp(-1.1/KNQ))  
        VSEDQ=PDENSA(LCELL)*DCONST2*CUNQ*DQ*DQ  
        SCQ=NU*DQ/DCONST1/CUNQ  
        Eff_dif=SCQ**(-TWO3)    
        STQ=DCONST3*PDENSA(LCELL)*DQ**2  
        Eff_imp=(STQ/(0.8+STQ))**2   

        Eff_int=(0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 
        RBcor=1. 
        vdplim=UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor
        vdplim=min(vdplim,.02)
        RSURFQ=RA(LCELL)+1./vdplim





        SUM0=SUM0+WGAUS(N)*(VSEDQ + 1./RSURFQ)  
        SUM3=SUM3+WGAUS(N)*(VSEDQ + 1./RSURFQ)*DQ**3  
      IF(iprnt.eq.1)THEN
      print *,'N,Y_GQ,WGAUS,DQ,KNQ,CUNQ,VSEDQ,SCQ,STQ,RSURFQ'
      print *,N,Y_GQ(N),WGAUS(N),DQ,KNQ,CUNQ,VSEDQ,SCQ,STQ,RSURFQ
      print *,'N,Eff_dif,imp,int,SUM0,SUM3'
      print *,N,Eff_dif,Eff_imp,Eff_int,SUM0,SUM3
      endif
        ENDDO
        VDEP(LCELL, VDNACC) = SUM0/sqrtpi  
        VDEP(LCELL, VDMACC) = SUM3/(sqrtpi*EXP((1.5*sqrt2*xxlsga)**2)*DGACC(LCELL)**3) 
        

        
        SUM0=0.
        SUM3=0.
        DO N=1,NGAUSdv
        DQ=DGCOR(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsgc)  
        KNQ=2.*XLM(LCELL)/DQ  
        CUNQ=1.+KNQ*(1.257+.4*exp(-1.1/KNQ))  
        VSEDQ=PDENSC(LCELL)*DCONST2*CUNQ*DQ*DQ  
        SCQ=NU*DQ/DCONST1/CUNQ  
        Eff_dif=SCQ**(-TWO3)    
        STQ=DCONST3*PDENSC(LCELL)*DQ**2  
        Eff_imp=(STQ/(0.8+STQ))**2   

        Eff_int=(0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 
        EFF_int=min(1.,EFF_int)
        RBcor=exp(-2.0*(STQ**0.5)) 
        vdplim=UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor
        vdplim=min(vdplim,.02)
        vdplim=max(vdplim,1e-35) 
        RSURFQ=RA(LCELL)+1./vdplim





        SUM0=SUM0+WGAUS(N)*(VSEDQ + 1./RSURFQ)  
        SUM3=SUM3+WGAUS(N)*(VSEDQ + 1./RSURFQ)*DQ**3  
        ENDDO
        VDEP(LCELL, VDNCOR) = SUM0/sqrtpi  
        VDEP(LCELL, VDMCOR) = SUM3/(sqrtpi*EXP((1.5*sqrt2*xxlsgc)**2)*DGCOR(LCELL)**3) 
        
        END DO  
             
        ENDIF  
        


         DO LCELL = 1, NUMCELLS
         
            DCONST2 = GRAV / ( 18.0 * AMU(LCELL) )
            DCONST3N = DCONST2 * PDENSN(LCELL) * DGNUC( LCELL )**2
            DCONST3A = DCONST2 * PDENSA(LCELL) * DGACC( LCELL )**2
            DCONST3C = DCONST2 * PDENSC(LCELL) * DGCOR( LCELL )**2
               

            VSED( LCELL, VSNNUC) = DCONST3N                         &
               * ( ESN16 + BHAT * KNNUC( LCELL ) * ESN04 )
            VSED( LCELL, VSMNUC) = DCONST3N                         &
               * (ESN64 + BHAT * KNNUC( LCELL ) * ESN28 )
        

            VSED( LCELL, VSNACC) = DCONST3A                          &
              * ( ESA16 + BHAT * KNACC( LCELL ) * ESA04 )
            VSED( LCELL, VSMACC) = DCONST3A                          &
              * ( ESA64 + BHAT * KNACC( LCELL ) * ESA28 )


            VSED( LCELL, VSNCOR) = DCONST3C                          &
              * ( ESC16 + BHAT * KNCOR( LCELL ) * ESC04 )
            VSED( LCELL, VSMCOR) = DCONST3C                          &
              * ( ESC64 + BHAT * KNCOR( LCELL ) * ESC28 )

         END DO


END SUBROUTINE VDVG_2




    SUBROUTINE aerosols_sorgam_init(chem,convfac,z_at_w,             &
         pm2_5_dry,pm2_5_water,pm2_5_dry_ec,                         &
         chem_in_opt,aer_ic_opt, is_aerosol,                         &
         ids,ide, jds,jde, kds,kde,                                  &
         ims,ime, jms,jme, kms,kme,                                  &
         its,ite, jts,jte, kts,kte, config_flags                     )

    USE module_configure,only:  grid_config_rec_type

    
 
   implicit none
   INTEGER,      INTENT(IN   ) :: chem_in_opt,aer_ic_opt
   INTEGER,      INTENT(IN   ) ::                               &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte
   LOGICAL, INTENT(OUT) :: is_aerosol(num_chem)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme, num_chem ) ,     &
          INTENT(INOUT   ) ::                                      &
                              chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) ,               &
          INTENT(INOUT      ) ::                                   &
                     pm2_5_dry,pm2_5_water,pm2_5_dry_ec
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) ,               &
          INTENT(IN      ) ::                                      &
                   convfac
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) ,               &
          INTENT(IN         ) ::                                   &
                     z_at_w
   TYPE (grid_config_rec_type) , INTENT (in) :: config_flags


     integer i,j,k,l,ii,jj,kk     
     real tempfac,mwso4,zz

      REAL splitfac
                        
      REAL so4vaptoaer


      REAL m3nuc

      REAL m3acc

      REAL m3cor
      DATA splitfac/.98/
      DATA so4vaptoaer/.999/



        xxlsgn = log(sginin)
        xxlsga = log(sginia)
        xxlsgc = log(sginic)

        l2sginin = xxlsgn**2
        l2sginia = xxlsga**2
        l2sginic = xxlsgc**2

        en1 = exp(0.125*l2sginin)
        ea1 = exp(0.125*l2sginia)
        ec1 = exp(0.125*l2sginic)

        esn04 = en1**4
        esa04 = ea1**4
        esc04 = ec1**4

        esn05 = esn04*en1
        esa05 = esa04*ea1

        esn08 = esn04*esn04
        esa08 = esa04*esa04
        esc08 = esc04*esc04

        esn09 = esn04*esn05
        esa09 = esa04*esa05

        esn12 = esn04*esn04*esn04
        esa12 = esa04*esa04*esa04
        esc12 = esc04*esc04*esc04

        esn16 = esn08*esn08
        esa16 = esa08*esa08
        esc16 = esc08*esc08

        esn20 = esn16*esn04
        esa20 = esa16*esa04
        esc20 = esc16*esc04

        esn24 = esn12*esn12
        esa24 = esa12*esa12
        esc24 = esc12*esc12

        esn25 = esn16*esn09
        esa25 = esa16*esa09

        esn28 = esn20*esn08
        esa28 = esa20*esa08
        esc28 = esc20*esc08


        esn32 = esn16*esn16
        esa32 = esa16*esa16
        esc32 = esc16*esc16

        esn36 = esn16*esn20
        esa36 = esa16*esa20
        esc36 = esc16*esc20

        esn49 = esn25*esn20*esn04
        esa49 = esa25*esa20*esa04

        esn52 = esn16*esn36
        esa52 = esa16*esa36

        esn64 = esn32*esn32
        esa64 = esa32*esa32
        esc64 = esc32*esc32

        esn100 = esn36*esn64

        esnm20 = 1.0/esn20
        esam20 = 1.0/esa20
        escm20 = 1.0/esc20

        esnm32 = 1.0/esn32
        esam32 = 1.0/esa32
        escm32 = 1.0/esc32


        xxm3 = 3.0*xxlsgn/ sqrt2

        nummin_i = facatkn_min*so4fac*aeroconcmin/(dginin**3*esn36)

        nummin_j = facacc_min*so4fac*aeroconcmin/(dginia**3*esa36)

        nummin_c = anthfac*aeroconcmin/(dginic**3*esc36)








        factnumn = exp(4.5*log(sgem_i)**2)/dgvem_i**3
        factnuma = exp(4.5*log(sgem_j)**2)/dgvem_j**3
        factnumc = exp(4.5*log(sgem_c)**2)/dgvem_c**3
        ccofm = alphsulf*sqrt(pirs*rgasuniv/(2.0*mwh2so4))
        ccofm_org = alphaorg*sqrt(pirs*rgasuniv/(2.0*mworg))
        mwso4=96.03






        
        
        

        pm2_5_dry(its:ite, kts:kte-1, jts:jte)    = 0.
        pm2_5_water(its:ite, kts:kte-1, jts:jte)  = 0.
        pm2_5_dry_ec(its:ite, kts:kte-1, jts:jte) = 0.



        Y_GQ(1)=-2.651961356835233
        WGAUS(1)=0.0009717812450995
        Y_GQ(2)=-1.673551628767471
        WGAUS(2)=0.05451558281913
        Y_GQ(3)=-0.816287882858965
        WGAUS(3)=0.4256072526101
        Y_GQ(4)=-0.0
        WGAUS(4)=0.8102646175568
        Y_GQ(5)=0.816287882858965
        WGAUS(5)=WGAUS(3)
        Y_GQ(6)=1.673551628767471
        WGAUS(6)=WGAUS(2)
        Y_GQ(7)=2.651961356835233
        WGAUS(7)=WGAUS(1)





        if(chem_in_opt == 1 ) return
        do l=p_so4aj,num_chem
         chem(ims:ime,kms:kme,jms:jme,l)=epsilc
        enddo
        chem(ims:ime,kms:kme,jms:jme,p_nu0)=1.e8
        chem(ims:ime,kms:kme,jms:jme,p_ac0)=1.e8
        do j=jts,jte
        jj=min(jde-1,j)
        do k=kts,kte-1
        kk=min(kde-1,k)
        do i=its,ite
        ii=min(ide-1,i)


        if( aer_ic_opt == AER_IC_DEFAULT ) then
          chem(i,k,j,p_so4aj)=chem(ii,kk,jj,p_sulf)*CONVFAC(ii,kk,jj)*MWSO4*splitfac*so4vaptoaer
          chem(i,k,j,p_so4ai)=chem(ii,kk,jj,p_sulf)*CONVFAC(ii,kk,jj)*MWSO4* &
        (1.-splitfac)*so4vaptoaer
          chem(i,k,j,p_sulf)=chem(ii,kk,jj,p_sulf)*(1.-so4vaptoaer)
          chem(i,k,j,p_nh4aj) = 10.E-05
          chem(i,k,j,p_nh4ai) = 10.E-05
          chem(i,k,j,p_no3aj) = 10.E-05
          chem(i,k,j,p_no3ai) = 10.E-05
          chem(i,k,j,p_naaj)  = 10.E-05
          chem(i,k,j,p_naai)  = 10.E-05
          chem(i,k,j,p_claj)  = 10.E-05
          chem(i,k,j,p_clai)  = 10.E-05
        elseif( aer_ic_opt == AER_IC_PNNL ) then
           zz = (z_at_w(ii,k,jj)+z_at_w(ii,k+1,jj))*0.5
           call sorgam_init_aer_ic_pnnl(   &
                chem, zz, i,k,j, ims,ime,jms,jme,kms,kme )
        else
           call wrf_error_fatal3("<stdin>",7557,&
                "aerosols_sorgam_init: unable to parse aer_ic_opt" )
        end if


      m3nuc = so4fac*chem(i,k,j,p_so4ai) + nh4fac*chem(i,k,j,p_nh4ai) + &
        no3fac*chem(i,k,j,p_no3ai) +                                    &
        nafac*chem(i,k,j,p_naai) + clfac*chem(i,k,j,p_clai) +           &
        orgfac*chem(i,k,j,p_orgaro1i) + &
        orgfac*chem(i,k,j,p_orgaro2i) + orgfac*chem(i,k,j,p_orgalk1i) + &
        orgfac*chem(i,k,j,p_orgole1i) + orgfac*chem(i,k,j,p_orgba1i) + &
        orgfac*chem(i,k,j,p_orgba2i) + orgfac*chem(i,k,j,p_orgba3i) + &
        orgfac*chem(i,k,j,p_orgba4i) + orgfac*chem(i,k,j,p_orgpai) + &
        anthfac*chem(i,k,j,p_p25i) + anthfac*chem(i,k,j,p_eci)


      m3acc = so4fac*chem(i,k,j,p_so4aj) + nh4fac*chem(i,k,j,p_nh4aj) + &
        no3fac*chem(i,k,j,p_no3aj) +                                    &
        nafac*chem(i,k,j,p_naaj) + clfac*chem(i,k,j,p_claj) +           & 
        orgfac*chem(i,k,j,p_orgaro1j) + &
        orgfac*chem(i,k,j,p_orgaro2j) + orgfac*chem(i,k,j,p_orgalk1j) + &
        orgfac*chem(i,k,j,p_orgole1j) + orgfac*chem(i,k,j,p_orgba1j) + &
        orgfac*chem(i,k,j,p_orgba2j) + orgfac*chem(i,k,j,p_orgba3j) + &
        orgfac*chem(i,k,j,p_orgba4j) + orgfac*chem(i,k,j,p_orgpaj) + &
        anthfac*chem(i,k,j,p_p25j) + anthfac*chem(i,k,j,p_ecj)


      m3cor = soilfac*chem(i,k,j,p_soila) + seasfac*chem(i,k,j,p_seas) + &
        anthfac*chem(i,k,j,p_antha)



      chem(i,k,j,p_nu0) = m3nuc/((dginin**3)*esn36)

      chem(i,k,j,p_ac0) = m3acc/((dginia**3)*esa36)
        
      chem(i,k,j,p_corn) = m3cor/((dginic**3)*esc36)

      enddo
      enddo
      enddo


    return
    END SUBROUTINE aerosols_sorgam_init
















    SUBROUTINE sorgam_init_aer_ic_pnnl(                  &
         chem, z, i,k,j, ims,ime, jms,jme, kms,kme )

      USE module_configure,only:  num_chem, grid_config_rec_type
      implicit none

      INTEGER,INTENT(IN   ) :: i,k,j,                           &
                               ims,ime, jms,jme, kms,kme
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme, num_chem ),&
           INTENT(INOUT   ) :: chem

      real, intent(in     ) :: z
      real :: mult



















































        if( z <= 500. ) then
           mult = 1.0
        elseif( z > 500. &
             .and. z <= 1000. ) then
           mult = 1.0 - 0.001074*(z-500.)
        elseif( z > 1000. &
             .and. z <= 5000. ) then
           mult = 0.463 - 0.000111*(z-1000.)
        else
           mult = 0.019
        end if



































      chem(i,k,j,p_sulf)     = mult*conmin
      chem(i,k,j,p_so4aj)    = mult*0.300*0.97
      chem(i,k,j,p_so4ai)    = mult*0.300*0.03
      chem(i,k,j,p_nh4aj)    = mult*0.094*0.97
      chem(i,k,j,p_nh4ai)    = mult*0.094*0.03
      chem(i,k,j,p_no3aj)    = mult*0.001*0.97
      chem(i,k,j,p_no3ai)    = mult*0.001*0.03
      chem(i,k,j,p_naaj)     = 10.E-05
      chem(i,k,j,p_naai)     = 10.E-05
      chem(i,k,j,p_claj)     = 10.E-05
      chem(i,k,j,p_clai)     = 10.E-05
      chem(i,k,j,p_ecj)      = mult*0.013*0.97
      chem(i,k,j,p_eci)      = mult*0.013*0.03
      chem(i,k,j,p_p25j)     = mult*4.500*0.97
      chem(i,k,j,p_p25i)     = mult*4.500*0.03
      chem(i,k,j,p_antha)    = mult*4.500/2.0
      chem(i,k,j,p_orgpaj)   = mult*0.088*0.97
      chem(i,k,j,p_orgpai)   = mult*0.088*0.03
      chem(i,k,j,p_orgaro1j) = conmin
      chem(i,k,j,p_orgaro1i) = conmin
      chem(i,k,j,p_orgaro2j) = conmin
      chem(i,k,j,p_orgaro2i) = conmin
      chem(i,k,j,p_orgalk1j) = conmin
      chem(i,k,j,p_orgalk1i) = conmin
      chem(i,k,j,p_orgole1j) = conmin
      chem(i,k,j,p_orgole1i) = conmin
      chem(i,k,j,p_orgba1j)  = conmin
      chem(i,k,j,p_orgba1i)  = conmin
      chem(i,k,j,p_orgba2j)  = conmin
      chem(i,k,j,p_orgba2i)  = conmin
      chem(i,k,j,p_orgba3j)  = conmin
      chem(i,k,j,p_orgba3i)  = conmin
      chem(i,k,j,p_orgba4j)  = conmin
      chem(i,k,j,p_orgba4i)  = conmin
      chem(i,k,j,p_seas)     = mult*1.75


    END SUBROUTINE sorgam_init_aer_ic_pnnl



SUBROUTINE sorgam_addemiss(                                             &
     id, dtstep, u10, v10, alt, dz8w, xland, chem,                      &
     ebu,                                   &
     slai,ust,smois,ivgtyp,isltyp,                                      &
     emis_ant,dust_emiss_active,                                        &
     seasalt_emiss_active,kemit,biom,num_soil_layers,emissopt,          &
     dust_opt,ktau,p8w,u_phy,v_phy,rho_phy,g,dx,erod,                   &
     ids,ide, jds,jde, kds,kde,                                         &
     ims,ime, jms,jme, kms,kme,                                         &
     its,ite, jts,jte, kts,kte                                          )







  USE module_state_description, only:  num_chem

  INTEGER,      INTENT(IN   ) :: seasalt_emiss_active, kemit,emissopt,  &
                                 dust_emiss_active,num_soil_layers,id,  &
                                 ktau,dust_opt,                         &
                                 biom,ids,ide, jds,jde, kds,kde,        &
                                 ims,ime, jms,jme, kms,kme,             &
                                 its,ite, jts,jte, kts,kte

  REAL, INTENT(IN   ) ::    dtstep


  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),               &
       INTENT(INOUT ) ::   chem



   REAL, DIMENSION( ims:ime, kms:kemit, jms:jme,num_emis_ant ),        &
         INTENT(IN    ) ::                                             &
         emis_ant



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,num_ebu ),             &
         INTENT(IN    ) ::                                             &
         ebu


  REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                      &
       INTENT(IN   ) ::                                                 &
       alt, dz8w


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                        &
         INTENT(IN    ) :: p8w,u_phy,v_phy,rho_phy 
   REAL, INTENT(IN    ) :: dx, g
   REAL, DIMENSION( ims:ime, jms:jme, 3 ),                              &
         INTENT(IN    ) :: erod


  REAL,  DIMENSION( ims:ime , jms:jme ),                                &
       INTENT(IN   ) ::                                                 &
       u10, v10, xland, slai, ust
  INTEGER,  DIMENSION( ims:ime , jms:jme ),                             &
       INTENT(IN   ) ::   ivgtyp, isltyp
  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ),    &
       INTENT(INOUT) ::   smois


  real, dimension(its:ite,kts:kte,jts:jte) :: factor




  factor(its:ite,kts:kte,jts:jte) = alt(its:ite,kts:kte,jts:jte)*dtstep/ &
                  dz8w(its:ite,kts:kte,jts:jte)


 if (emissopt .ne. 5) then



  chem(its:ite,kts:kemit,jts:jte,p_nu0) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_nu0) +                   &
       factor(its:ite,kts:kemit,jts:jte)*factnumn*(              &
       anthfac*( emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm25i) + &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_eci)  +            &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_naai) ) +          &
       orgfac*emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi) +     &
       so4fac*emis_ant(its:ite,kts:kemit,jts:jte,p_e_so4i)   +   &
       no3fac*emis_ant(its:ite,kts:kemit,jts:jte,p_e_no3i) )


  
  chem(its:ite,kts:kemit,jts:jte,p_ac0) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_ac0) +                   &
       factor(its:ite,kts:kemit,jts:jte)*factnuma*(              &
       anthfac*( emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm25j) + &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_ecj)  +            &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_naaj) ) +          &
       orgfac*emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj) +     &
       so4fac*emis_ant(its:ite,kts:kemit,jts:jte,p_e_so4j)   +   &
       no3fac*emis_ant(its:ite,kts:kemit,jts:jte,p_e_no3j))



  chem(its:ite,kts:kemit,jts:jte,p_corn) =                       &
       chem(its:ite,kts:kemit,jts:jte,p_corn) +                  &
       factor(its:ite,kts:kemit,jts:jte)*factnumc*anthfac*                           &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm_10)



  chem(its:ite,kts:kemit,jts:jte,p_antha) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_antha) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm_10)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_p25j) =                       &
       chem(its:ite,kts:kemit,jts:jte,p_p25j) +                  &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm25j)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_p25i) =                       &
       chem(its:ite,kts:kemit,jts:jte,p_p25i) +                  &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm25i)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_ecj) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_ecj) +                   &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_ecj)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_eci) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_eci) +                   &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_eci)*factor(its:ite,kts:kemit,jts:jte)
  chem(its:ite,kts:kemit,jts:jte,p_naaj) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_naaj) +                   &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_naaj)*factor(its:ite,kts:kemit,jts:jte)
  chem(its:ite,kts:kemit,jts:jte,p_naai) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_naai) +                   &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_naai)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_orgpaj) =                     &
       chem(its:ite,kts:kemit,jts:jte,p_orgpaj) +                &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_orgpai) =                     &
       chem(its:ite,kts:kemit,jts:jte,p_orgpai) +                &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_so4aj) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_so4aj) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_so4j)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_so4ai) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_so4ai) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_so4i)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_no3aj) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_no3aj) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_no3j)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_no3ai) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_no3ai) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_no3i)*factor(its:ite,kts:kemit,jts:jte)
  elseif(emissopt == 5)then



  chem(its:ite,kts:kemit,jts:jte,p_nu0) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_nu0) +                   &
       factor(its:ite,kts:kemit,jts:jte)*factnumn*(                                  &
       anthfac*( .25*emis_ant(its:ite,kts:kemit,jts:jte,p_e_bc) ) +                      &
       orgfac*.25*emis_ant(its:ite,kts:kemit,jts:jte,p_e_oc) )


  
  chem(its:ite,kts:kemit,jts:jte,p_ac0) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_ac0) +                   &
       factor(its:ite,kts:kemit,jts:jte)*factnuma*(                                  &
       anthfac*( .75*emis_ant(its:ite,kts:kemit,jts:jte,p_e_bc) ) +                      &
       orgfac*.75*emis_ant(its:ite,kts:kemit,jts:jte,p_e_oc) )





  chem(its:ite,kts:kemit,jts:jte,p_ecj) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_ecj) +                   &
       .75*emis_ant(its:ite,kts:kemit,jts:jte,p_e_bc)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_eci) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_eci) +                   &
       .25*emis_ant(its:ite,kts:kemit,jts:jte,p_e_bc)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_orgpaj) =                     &
       chem(its:ite,kts:kemit,jts:jte,p_orgpaj) +                &
       .75*emis_ant(its:ite,kts:kemit,jts:jte,p_e_oc)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_orgpai) =                     &
       chem(its:ite,kts:kemit,jts:jte,p_orgpai) +                &
       .25*emis_ant(its:ite,kts:kemit,jts:jte,p_e_oc)*factor(its:ite,kts:kemit,jts:jte)

  endif


  if(biom == 1 )then



  chem(its:ite,kts:kte,jts:jte,p_nu0) =                        &
       chem(its:ite,kts:kte,jts:jte,p_nu0) +                   &
       factor(its:ite,kts:kte,jts:jte)*factnumn*(              &
       anthfac*( .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_pm25) +       &
              .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_bc) ) +          &
       orgfac*.25*ebu(its:ite,kts:kte,jts:jte,p_ebu_oc) )


  
  chem(its:ite,kts:kte,jts:jte,p_ac0) =                        &
       chem(its:ite,kts:kte,jts:jte,p_ac0) +                   &
       factor(its:ite,kts:kte,jts:jte)*factnuma*(              &
       anthfac*(.75*ebu(its:ite,kts:kte,jts:jte,p_ebu_pm25) +        &
      .75*ebu(its:ite,kts:kte,jts:jte,p_ebu_bc) ) +                  &
       orgfac*.75*ebu(its:ite,kts:kte,jts:jte,p_ebu_oc) )

  chem(its:ite,kts:kte,jts:jte,p_corn) =                     &
       chem(its:ite,kts:kte,jts:jte,p_corn) +                  &
       factor(its:ite,kts:kte,jts:jte)*factnumc*anthfac*       &
       ebu(its:ite,kts:kte,jts:jte,p_ebu_pm10)





  chem(its:ite,kts:kte,jts:jte,p_ecj) =                        &
       chem(its:ite,kts:kte,jts:jte,p_ecj) +                   &
       .75*ebu(its:ite,kts:kte,jts:jte,p_ebu_bc)*factor(its:ite,kts:kte,jts:jte)

  chem(its:ite,kts:kte,jts:jte,p_eci) =                        &
       chem(its:ite,kts:kte,jts:jte,p_eci) +                   &
       .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_bc)*factor(its:ite,kts:kte,jts:jte)

  chem(its:ite,kts:kte,jts:jte,p_orgpaj) =                     &
       chem(its:ite,kts:kte,jts:jte,p_orgpaj) +                &
       .75*ebu(its:ite,kts:kte,jts:jte,p_ebu_oc)*factor(its:ite,kts:kte,jts:jte)

  chem(its:ite,kts:kte,jts:jte,p_orgpai) =                     &
       chem(its:ite,kts:kte,jts:jte,p_orgpai) +                &
       .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_oc)*factor(its:ite,kts:kte,jts:jte)

  chem(its:ite,kts:kte,jts:jte,p_antha) =                      &
       chem(its:ite,kts:kte,jts:jte,p_antha) +                 &
       ebu(its:ite,kts:kte,jts:jte,p_ebu_pm10)*factor(its:ite,kts:kte,jts:jte)

  chem(its:ite,kts:kte,jts:jte,p_p25j) =                       &
       chem(its:ite,kts:kte,jts:jte,p_p25j) +                  &
       .75*ebu(its:ite,kts:kte,jts:jte,p_ebu_pm25)*factor(its:ite,kts:kte,jts:jte)

  chem(its:ite,kts:kte,jts:jte,p_p25i) =                       &
       chem(its:ite,kts:kte,jts:jte,p_p25i) +                  &
       .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_pm25)*factor(its:ite,kts:kte,jts:jte)

   endif 



  if( seasalt_emiss_active == 1 ) then
     call sorgam_seasalt_emiss(                                  &
          dtstep, u10, v10, alt, dz8w, xland, chem,              &
          ids,ide, jds,jde, kds,kde,                             &
          ims,ime, jms,jme, kms,kme,                             &
          its,ite, jts,jte, kts,kte                              )
  end if
  if( seasalt_emiss_active == 2 ) then





  end if
  if( dust_opt == 2 ) then
   
   call wrf_message("WARNING: You are calling DUSTRAN dust emission scheme with MOSAIC, which is highly experimental and not recommended for use. Please use dust_opt==13")
   

      call sorgam_dust_emiss(                                     &
           slai, ust, smois, ivgtyp, isltyp,                      &
           id, dtstep, u10, v10, alt, dz8w,                       &
           xland, num_soil_layers, chem,                          &
           ids,ide, jds,jde, kds,kde,                             &
           ims,ime, jms,jme, kms,kme,                             &
           its,ite, jts,jte, kts,kte                              )
  end if
  

  if( dust_opt == 13 ) then
  
      call sorgam_dust_gocartemis(                                &
           ktau,dtstep,num_soil_layers,alt,u_phy,                 &
           v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,        &
           ivgtyp,isltyp,xland,dx,g,                              &
           ids,ide, jds,jde, kds,kde,                             &
           ims,ime, jms,jme, kms,kme,                             &
           its,ite, jts,jte, kts,kte                              )
  end if

END SUBROUTINE sorgam_addemiss


SUBROUTINE sorgam_seasalt_emiss(                                        &
     dtstep, u10, v10, alt, dz8w, xland, chem,                          &
     ids,ide, jds,jde, kds,kde,                                         &
     ims,ime, jms,jme, kms,kme,                                         &
     its,ite, jts,jte, kts,kte                                          )






   USE module_mosaic_addemiss, only:    seasalt_emitfactors_1bin

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,            &
                                  ims,ime, jms,jme, kms,kme,            &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep


   REAL,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   u10, v10, xland


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),              &
         INTENT(INOUT ) ::   chem



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                     &
         INTENT(IN   ) ::   alt, dz8w


   integer :: i, j, k, l, l_na, l_cl, n
    integer :: p1st

    real :: dum, dumdlo, dumdhi, dumoceanfrac, dumspd10
    real :: factaa, factbb, fraccl, fracna

    real :: ssemfact_numb_i, ssemfact_numb_j, ssemfact_numb_c
    real :: ssemfact_mass_i, ssemfact_mass_j, ssemfact_mass_c








    ssemfact_numb_i = 0.
    ssemfact_mass_i = 0.







    dumdlo = 0.1e-4
    dumdhi = 1.250e-4
    call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
         ssemfact_numb_j, dum, ssemfact_mass_j )


    dumdlo = 1.25e-4
    dumdhi = 10.0e-4
    call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
         ssemfact_numb_c, dum, ssemfact_mass_c )


    ssemfact_mass_i = ssemfact_mass_i*1.0e6
    ssemfact_mass_j = ssemfact_mass_j*1.0e6
    ssemfact_mass_c = ssemfact_mass_c*1.0e6


    k = kts
    do j = jts, jte
    do i = its, ite

    
    
    
       if( xland(i,j) < 1.5 ) cycle

    
    
       dumoceanfrac = 1. 
       dumspd10 = dumoceanfrac* &
            ( (u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5*3.41) )




       factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
       factbb = factaa * dumspd10


       fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
       fraccl = 1.0 - fracna


       chem(i,k,j,p_naai) = chem(i,k,j,p_naai) +   &
                            factbb * ssemfact_mass_i * fracna
       chem(i,k,j,p_clai) = chem(i,k,j,p_clai) +   &
                            factbb * ssemfact_mass_i * fraccl
       chem(i,k,j,p_nu0)  = chem(i,k,j,p_nu0) +   &
                            factbb * ssemfact_numb_i

       chem(i,k,j,p_naaj) = chem(i,k,j,p_naaj) +   &
                            factbb * ssemfact_mass_j * fracna
       chem(i,k,j,p_claj) = chem(i,k,j,p_claj) +   &
                            factbb * ssemfact_mass_j * fraccl
       chem(i,k,j,p_ac0)  = chem(i,k,j,p_ac0) +   &
                            factbb * ssemfact_numb_j

       chem(i,k,j,p_seas) = chem(i,k,j,p_seas) +   &
                            factbb * ssemfact_mass_c
       chem(i,k,j,p_corn) = chem(i,k,j,p_corn) +   &
                            factbb * ssemfact_numb_c
    end do 
    end do 
END SUBROUTINE sorgam_seasalt_emiss


   subroutine sorgam_dust_emiss(  slai,ust, smois, ivgtyp, isltyp,         &
               id, dtstep, u10, v10, alt, dz8w, xland, num_soil_layers,    &
               chem,                                                       &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )


























   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem, param_first_scalar
   USE module_data_mosaic_asect

   IMPLICIT NONE



   INTEGER,      INTENT(IN   ) :: id,num_soil_layers,                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep


   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland, slai, ust
   INTEGER,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   ivgtyp, isltyp


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem



   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w

   REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,     &
          INTENT(INOUT) ::   smois


        integer i, j, k, l, l_oin, l_ca, l_co3, n, ii
        integer iphase, itype, izob
        integer p1st

        real dum, dumdlo, dumdhi, dumlandfrac, dumspd10
        real factaa, factbb, fracoin, fracca, fracco3, fractot
        real ustart, ustar1, ustart0
        real alphamask, f8, f50, f51, f52, wetfactor, sumdelta, ftot
        real smois_grav, wp, pclay
        real :: beta(4,7)
        real :: gamma(4), delta(4)
        real :: sz(8)
        real :: dustflux, densdust, mass1part
        real :: dp_meanvol_tmp








        beta(1,1)=0.12
        beta(2,1)=0.04
        beta(3,1)=0.04
        beta(4,1)=0.80
        beta(1,2)=0.34
        beta(2,2)=0.28
        beta(3,2)=0.28
        beta(4,2)=0.10
        beta(1,3)=0.45
        beta(2,3)=0.15
        beta(3,3)=0.15
        beta(4,3)=0.25
        beta(1,4)=0.12
        beta(2,4)=0.09
        beta(3,4)=0.09
        beta(4,4)=0.70
        beta(1,5)=0.40
        beta(2,5)=0.05
        beta(3,5)=0.05
        beta(4,5)=0.50
        beta(1,6)=0.34
        beta(2,6)=0.18
        beta(3,6)=0.18
        beta(4,6)=0.30
        beta(1,7)=0.22
        beta(2,7)=0.09
        beta(3,7)=0.09
        beta(4,7)=0.60
        gamma(1)=0.08
        gamma(2)=1.00
        gamma(3)=1.00
        gamma(4)=0.12













        sz(1)=0.0
        sz(2)=0.0
        sz(3)=0.0005
        sz(4)=0.0095
        sz(5)=0.03
        sz(6)=0.10
        sz(7)=0.18
        sz(8)=0.68


        itype = 1
        iphase = ai_phase


        k = kts
        do 1830 j = jts, jte
        do 1820 i = its, ite

    if( xland(i,j) > 1.5 ) cycle



        dumlandfrac = 1.
        dumspd10=(u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5)
        if(dumspd10 >= 5.0) then
           dumspd10 = dumlandfrac* &
         ( dumspd10*dumspd10*(dumspd10-5.0))
         else
            dumspd10=0.
         endif












         alphamask=0.001
         if (ivgtyp(i,j) .eq. 7) then
           f8=0.005
           f50=0.00
           f51=0.10
           f51=0.066
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i,j) .eq. 8) then
           f8=0.010
           f50=0.00
           f51=0.00
           f52=0.15
           f52=0.10
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i,j) .eq. 10) then
           f8=0.00
           f50=0.00
           f51=0.01
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif


















         izob=0
         if(isltyp(i,j).eq.1) izob=1
         if(isltyp(i,j).eq.2) izob=1
         if(isltyp(i,j).eq.3) izob=4
         if(isltyp(i,j).eq.4) izob=2
         if(isltyp(i,j).eq.5) izob=2
         if(isltyp(i,j).eq.6) izob=2
         if(isltyp(i,j).eq.7) izob=7
         if(isltyp(i,j).eq.8) izob=2
         if(isltyp(i,j).eq.9) izob=6
         if(isltyp(i,j).eq.10) izob=5
         if(isltyp(i,j).eq.11) izob=2
         if(isltyp(i,j).eq.12) izob=3
         if(isltyp(i,j).ge.13) izob=0
         if(izob.eq.0) goto 1840



         do ii=1,4
           delta(ii)=0.0
         enddo
         sumdelta=0.0
         do ii=1,4
           delta(ii)=beta(ii,izob)*gamma(ii)
           if(ii.lt.4) then
             sumdelta=sumdelta+delta(ii)
           endif
         enddo
         do ii=1,4
           delta(ii)=delta(ii)/sumdelta
         enddo









         pclay=beta(1,izob)*100.
         wp=0.0014*pclay*pclay+0.17*pclay
         smois_grav=(smois(i,1,j)/2.6)*100.
         if(smois_grav.gt.wp) then
           wetfactor=sqrt(1.0+1.21*(smois_grav-wp)**0.68)
         else
           wetfactor=1.0
         endif






         ustar1=ust(i,j)*100.0
         if(ustar1.gt.100.0) ustar1=100.0
         ustart0=20.0
         ustart=ustart0*wetfactor
         if(ustar1.le.ustart) then
           dustflux=0.0
         else
           dustflux=1.0e-14*(ustar1**4)*(1.0-(ustart/ustar1))
         endif
         dustflux=dustflux*10.0

         ftot=0.0
         do ii=1,2
           ftot=ftot+dustflux*alphamask*delta(ii)
         enddo

         ftot=ftot*1.0e+09


         factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
         factbb = factaa * ftot
         fracoin = 1.00


         fracca = 0.0
         fracco3 = 0.0
         fractot = fracoin + fracca + fracco3

         chem(i,k,j,p_p25j)=chem(i,k,j,p_p25j) +   &
            factbb * (sz(3)+sz(4)+sz(5)+sz(6)) * fractot

         chem(i,k,j,p_soila)=chem(i,k,j,p_soila) +   &
            factbb * (sz(7)+sz(8)) * fractot


         densdust=2.5
         dp_meanvol_tmp = 1.0e2*dginia*exp(1.5*l2sginia) 
         mass1part=0.523598*(dp_meanvol_tmp**3)*densdust*1.0e06
         chem(i,k,j,p_ac0)=chem(i,k,j,p_ac0) +   &
            factbb * (sz(3)+sz(4)+sz(5)+sz(6)) * fractot / mass1part

         dp_meanvol_tmp = 1.0e2*dginic*exp(1.5*l2sginic) 
         mass1part=0.523598*(dp_meanvol_tmp**3)*densdust*1.0e06
         chem(i,k,j,p_corn)=chem(i,k,j,p_corn) +   &
            factbb * (sz(7)+sz(8)) * fractot / mass1part


1840    continue

1820    continue
1830    continue

        return

   END subroutine sorgam_dust_emiss




  subroutine sorgam_dust_gocartemis (ktau,dt,num_soil_layers,alt,u_phy,    &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,                   &
         ivgtyp,isltyp,xland,dx,g,                                         &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_data_gocart_dust
  USE module_configure
  USE module_state_description
  USE module_model_constants, ONLY: mwdry
  USE module_data_mosaic_asect
  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ktau, num_soil_layers,           &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(INOUT) ::                               smois
   REAL,  DIMENSION( ims:ime , jms:jme, 3 )                   ,               &
          INTENT(IN   ) ::    erod
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     xland
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                     dz8w,p8w,             &
                                              u_phy,v_phy,rho_phy

  REAL, INTENT(IN   ) :: dt,dx,g



  integer :: nmx,i,j,k,ndt,imx,jmx,lmx
  integer ilwi, start_month
  real*8, DIMENSION (3) :: erodin
  real*8, DIMENSION (5) :: bems
  real*8  w10m,gwet,airden,airmas
  real*8  cdustemis,jdustemis,cdustcon,jdustcon
  real*8  cdustdens,jdustdens,mass1part,jdustdiam,cdustdiam,dp_meanvol_tmp
  real*8  dxy
  real*8  conver,converi
  real dttt
  real soilfacj,rhosoilj,rhosoilc
  real totalemis,accfrac,corfrac,rscale1,rscale2
  
  accfrac=0.07              
  corfrac=0.93              
  rscale1=1.00  
  rscale2=1.02  
  accfrac=accfrac*rscale1
  corfrac=corfrac*rscale2

  rhosoilj=2.5e3
  rhosoilc=2.6e3
  soilfacj=soilfac*rhosoilj/rhosoilc

  conver=1.e-9
  converi=1.e9


  nmx=5
  k=kts
  do j=jts,jte
  do i=its,ite


     if(xland(i,j).lt.1.5)then

     ilwi=1
     start_month = 3   
     w10m=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
     airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g   


     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))

    
    
     erodin(1)=erod(i,j,1)
     erodin(2)=erod(i,j,2)
     erodin(3)=erod(i,j,3)


     gwet=smois(i,1,j)/porosity(isltyp(i,j))
     ndt=ifix(dt)
     airden=rho_phy(i,kts,j)
     dxy=dx*dx

    call sorgam_source_du( nmx, dt,i,j, &
                     erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,start_month,g)


    
    
    totalemis=(sum(bems(1:5))/dt)*converi/dxy 
     
     
    jdustemis = totalemis*accfrac   
    cdustemis = totalemis*corfrac   

         cdustcon = sum(bems(1:5))*corfrac/airmas  
         cdustcon = cdustcon * converi   
         jdustcon = sum(bems(1:5))*accfrac/airmas  
         jdustcon = jdustcon * converi   

         chem(i,k,j,p_p25j)=chem(i,k,j,p_p25j) + jdustcon 
         chem(i,k,j,p_soila)=chem(i,k,j,p_soila) + cdustcon




       chem(i,k,j,p_ac0) =  chem(i,k,j,p_ac0) + jdustcon * factnuma*soilfacj
       chem(i,k,j,p_corn) =  chem(i,k,j,p_corn) + cdustcon * factnumc*soilfac

     endif
  enddo
  enddo

end subroutine sorgam_dust_gocartemis

  SUBROUTINE sorgam_source_du( nmx, dt1,i,j, &
                     erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,month,g0)























 USE module_data_gocart_dust

  INTEGER, INTENT(IN)    :: nmx
  REAL*8,    INTENT(IN)  :: erod(ndcls)
  INTEGER, INTENT(IN)    :: ilwi,month

  REAL*8,    INTENT(IN)    :: w10m, gwet
  REAL*8,    INTENT(IN)    :: dxy
  REAL*8,    INTENT(IN)    :: airden, airmas
  REAL*8,    INTENT(OUT)   :: bems(nmx)

  REAL*8    :: den(nmx), diam(nmx)
  REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
  REAL, intent(in)    :: g0
  REAL    :: rhoa, g,dt1
  INTEGER :: i, j, n, m, k

  
  
   ch_dust(:,:)=1.0D-9  
  
  

  
  DO n = 1, nmx
     
     den(n) = den_dust(n)*1.0D-3
     diam(n) = 2.0*reff_dust(n)*1.0D2
     g = g0*1.0E2
     
     m = ipoint(n)
     tsrc = 0.0
              rhoa = airden*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)

              
             IF (gwet < 0.5) THEN  

                 u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet))))
              ELSE
                 
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(m)*dxy  
              IF (ilwi == 1 ) THEN
                 dsrc = ch_dust(n,month)*srce*w10m**2 &
                      * (w10m - u_ts)*dt1  
              ELSE
                 dsrc = 0.0
              END IF
              IF (dsrc < 0.0) dsrc = 0.0

              
              
              bems(n) = dsrc     

  ENDDO

END SUBROUTINE sorgam_source_du














































































END Module module_aerosols_sorgam



