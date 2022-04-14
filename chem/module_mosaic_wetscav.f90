MODULE module_mosaic_wetscav

CONTAINS



   subroutine wetscav_cbmz_mosaic (id,ktau,dtstep,ktauc,config_flags,      &
               dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,        &
	       qlsink,precr,preci,precs,precg,qsrflx,                      &
	       gas_aqfrac, numgas_aqfrac,                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )







   USE module_configure, only: grid_config_rec_type
   USE module_state_description
   USE module_data_mosaic_asect


   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc, numgas_aqfrac
   REAL, INTENT(IN   ) :: dtstep,dtstepc



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(IN ) ::                                   gas_aqfrac




   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
         INTENT(IN   ) ::                                           &
                                                        alt,        &
                                                      t_phy,        &
                                                      p_phy,        &
                                                   t8w,p8w,         &
	                            qlsink,precr,preci,precs,precg, &
                                                    rho_phy,cldfra
   REAL, DIMENSION( ims:ime, jms:jme, num_chem ),          &
         INTENT(OUT ) ::                                qsrflx 

   call wetscav (id,ktau,dtstep,ktauc,config_flags,                      &
        dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,             &
        qlsink,precr,preci,precs,precg,qsrflx,                           &
        gas_aqfrac, numgas_aqfrac,                                       &
        ntype_aer, nsize_aer, ncomp_aer,                                 &
        massptr_aer, dens_aer, numptr_aer,                               &
        maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase, cw_phase, &
        volumcen_sect, volumlo_sect, volumhi_sect,                       &
        waterptr_aer, dens_water_aer,                                    &
        scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow, &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte                                        )

   end subroutine wetscav_cbmz_mosaic




   subroutine wetscav_mozart_mosaic (id,ktau,dtstep,ktauc,config_flags,      &
               dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,        &
         qlsink,precr,preci,precs,precg,qsrflx,                      &
         gas_aqfrac, numgas_aqfrac,                                  &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )







   USE module_configure, only: grid_config_rec_type
   USE module_state_description
   USE module_data_mosaic_asect


   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc, numgas_aqfrac
   REAL, INTENT(IN   ) :: dtstep,dtstepc



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(IN ) ::                                   gas_aqfrac




   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
         INTENT(IN   ) ::                                           &
                                                        alt,        &
                                                      t_phy,        &
                                                      p_phy,        &
                                                   t8w,p8w,         &
                              qlsink,precr,preci,precs,precg, &
                                                    rho_phy,cldfra
   REAL, DIMENSION( ims:ime, jms:jme, num_chem ),          &
         INTENT(OUT ) ::                                qsrflx 

   call wetscav (id,ktau,dtstep,ktauc,config_flags,                      &
        dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,             &
        qlsink,precr,preci,precs,precg,qsrflx,                           &
        gas_aqfrac, numgas_aqfrac,                                       &
        ntype_aer, nsize_aer, ncomp_aer,                                 &
        massptr_aer, dens_aer, numptr_aer,                               &
        maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase, cw_phase, &
        volumcen_sect, volumlo_sect, volumhi_sect,                       &
        waterptr_aer, dens_water_aer,                                    &
        scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow, &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte                                        )

   end subroutine wetscav_mozart_mosaic



   subroutine wetscav (id,ktau,dtstep,ktauc,config_flags,                &
        dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,             &
        qlsink,precr,preci,precs,precg,qsrflx,                           &
        gas_aqfrac, numgas_aqfrac,                                       &
        ntype_aer, nsize_aer, ncomp_aer,                                 &
        massptr_aer, dens_aer, numptr_aer,                               &
        maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase, cw_phase, &
        volumcen_sect, volumlo_sect, volumhi_sect,                       &
        waterptr_aer, dens_water_aer,                                    &
        scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow, &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte                                        )







   USE module_configure, only: grid_config_rec_type
   USE module_state_description
   USE module_model_constants, only: g, rhowater, mwdry
   USE module_dep_simple, only: is_aerosol

   IMPLICIT NONE



   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc, numgas_aqfrac
      REAL,      INTENT(IN   ) :: dtstep,dtstepc



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(IN ) ::                                   gas_aqfrac




   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
          INTENT(IN   ) ::                                          &
                                                        alt,        &
                                                      t_phy,        &
                                                      p_phy,        &
                                                    t8w,p8w,        &
                             qlsink,precr,preci,precs,precg,        &
                                             rho_phy,cldfra
   integer, intent(in) :: maxd_atype, maxd_asize, maxd_acomp, maxd_aphase
   integer, intent(in) :: ai_phase, cw_phase
   integer, intent(in) :: ntype_aer
   integer, intent(in) :: nsize_aer( maxd_atype ),   & 
      	  ncomp_aer( maxd_atype ),   & 
      	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase), & 
      	  waterptr_aer( maxd_asize, maxd_atype ), & 
      	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) 
   real, intent(in) :: dens_aer( maxd_acomp, maxd_atype ),   & 
	  dens_water_aer 
   real, intent(in) :: volumcen_sect( maxd_asize, maxd_atype ),   & 
                        volumlo_sect( maxd_asize, maxd_atype ),   & 
                        volumhi_sect( maxd_asize, maxd_atype )      

   real, intent(in) :: dlndg_nimptblgrow
   integer, intent(in) :: nimptblgrow_mind, nimptblgrow_maxd
   real, intent(in) :: scavimptblnum(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype), &
     	     scavimptblvol(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)
   REAL, DIMENSION( ims:ime, jms:jme, num_chem ),          &
         INTENT(OUT ) ::                                qsrflx 




     integer :: i,j,k,l,m,n
     integer :: lmasscw,lnumcw
     real scale
   logical :: isprx(ims:ime,kms:kme,jms:jme)

   real :: pdel(ims:ime,kms:kme,jms:jme)
   real :: delpf(kms:kme)
   real :: delpfhalf
   real :: dump
   real :: fac_pwght_ls(kms:kme)   
   real :: fapincld, fapoutcld, faptot
   real :: fapmin_ls     
   real :: fapx(ims:ime,kms:kme,jms:jme)
   real :: fracscav
   real :: pf_above, pf_below, pdel_fac
   real :: pf_ls(kms:kme)
   real :: pfoutcld
   real :: pfsmall_ls    
                             
   real :: pfsmall_min
   real :: pfx(ims:ime,kms:kme,jms:jme)
   real :: pfx_inrain(ims:ime,kms:kme,jms:jme)
   real :: sumfa, sumpf, sumpffa
   REAL :: dqdt( ims:ime, kms:kme, jms:jme, num_chem )
   logical dotend(num_chem)

   dotend(:) = .false.
   dqdt(:,:,:,:) = 0.0
   qsrflx(:,:,:) = 0.0



      do 100 j=jts,jte
      do k=kts,kte
         do i=its,ite
            pdel(i,k,j)=p8w(i,k,j)-p8w(i,k+1,j)
         end do
      end do

      do 100 k=kts,kte
      do 100 i=its,ite
         scale=1.0-dtstepc*qlsink(i,k,j)
         scale=max(0.0,min(1.0,scale)) 
         if (qlsink(i,k,j) >  0.0) then
            pdel_fac = pdel(i,k,j)/g
            do n=1,ntype_aer
               do m=1,nsize_aer(n)
                  do l=1,ncomp_aer(n)
                     lmasscw=massptr_aer(l,m,n,cw_phase)
                     if (lmasscw < param_first_scalar) cycle
                     qsrflx(i,j,lmasscw)=qsrflx(i,j,lmasscw)+chem(i,k,j,lmasscw)*(scale-1.)*pdel_fac 
                     chem(i,k,j,lmasscw)=chem(i,k,j,lmasscw)*scale
                  end do 
                  lnumcw=numptr_aer(m,n,cw_phase)
                  if (lnumcw < param_first_scalar) cycle
                  qsrflx(i,j,lnumcw)=qsrflx(i,j,lnumcw)+chem(i,k,j,lnumcw)*(scale-1.)*pdel_fac 
                  chem(i,k,j,lnumcw)=chem(i,k,j,lnumcw)*scale
               end do 
            end do 
         end if    
  100 continue




  if ( .NOT. (config_flags%chem_opt == mozart_mosaic_4bin_aq_kpp) ) then
      do 290 l = param_first_scalar, min( num_chem, numgas_aqfrac )
      if ( is_aerosol(l) ) goto 290
      do 270 j = jts,jte
      do 270 k = kts,kte
      do 270 i = its,ite
         fracscav = dtstepc*qlsink(i,k,j)*gas_aqfrac(i,k,j,l)
         if (fracscav > 0.0) then 
            fracscav = min(1.0,fracscav) 
            scale = 1.0 - fracscav
            pdel_fac = pdel(i,k,j)/(g*mwdry)
            qsrflx(i,j,l) = qsrflx(i,j,l)+chem(i,k,j,l)*(scale-1.)*pdel_fac 
            chem(i,k,j,l) = chem(i,k,j,l)*scale
         end if
270   continue
290   continue
  end if









   fapmin_ls = 1.0e-3
   pfsmall_ls = 1.0e-7

   isprx(:,:,:) = .false.
   pfx(:,:,:) = 0.0
   pfx_inrain(:,:,:) = 0.0
   fapx(:,:,:) = 0.0


      do 5900 j=jts,jte
      do 5900 i=its,ite





      pf_below = 0.0
      do k = kte,kts,-1
         pf_above = pf_below
         pf_below = precr(i,k,j) + preci(i,k,j) + precs(i,k,j) + precg(i,k,j) 
         if (pf_below < pfsmall_ls) pf_below = 0.0
         delpf(k) = pf_below - pf_above
         pf_ls(k) = 0.5*(pf_below + pf_above)
         if (pf_ls(k) < pfsmall_ls) pf_ls(k) = 0.0
      end do




      do k = kte, kts,-1
         if (k == kte) then

            delpfhalf = 0.5*delpf(k)
            if (delpfhalf > 0.0) then
               fac_pwght_ls(k) = max( cldfra(i,k,j), fapmin_ls )
               sumpffa = delpfhalf * fac_pwght_ls(k)
               sumpf = delpfhalf
            else
               fac_pwght_ls(k) = fapmin_ls
               sumpffa = 0.0
               sumpf = 0.0
	    end if
         else

            delpfhalf = 0.5*delpf(k+1)
            if (delpfhalf > 0.0) then
               sumpffa = sumpffa + delpfhalf*max( cldfra(i,k+1,j), fapmin_ls )
               sumpf = sumpf + delpfhalf
               fac_pwght_ls(k) = max( (sumpffa/sumpf), fapmin_ls )
            else
               fac_pwght_ls(k) = fac_pwght_ls(k+1)
               sumpffa = max( (sumpffa + delpfhalf*fac_pwght_ls(k)), 0.0 )
               sumpf = max( (sumpf + delpfhalf), 0.0 )
	    end if

            delpfhalf = 0.5*delpf(k)
            if (delpfhalf > 0.0) then
               sumpffa = sumpffa + delpfhalf*max( cldfra(i,k,j), fapmin_ls )
               sumpf = sumpf + delpfhalf
               fac_pwght_ls(k) = max( (sumpffa/sumpf), fapmin_ls )
            else

               sumpffa = max( (sumpffa + delpfhalf*fac_pwght_ls(k)), 0.0 )
               sumpf = max( (sumpf + delpfhalf), 0.0 )
	    end if
	 end if
      end do









      do 4900 k = kte, kts,-1








      sumpf = 0.0
      sumfa = 0.0










      if (pf_ls(k) > 0.0) then
         faptot= 0.5*fac_pwght_ls(k)
         fapincld = min( faptot, cldfra(i,k,j) )
         fapoutcld = max( 0.0, faptot-fapincld )
         pfoutcld = pf_ls(k)*(fapoutcld/faptot)
         if (pfoutcld >= pfsmall_ls) then
            isprx(i,k,j) = .true.
            pfx(i,k,j) = pfoutcld
            fapx(i,k,j) = fapoutcld
            pfx_inrain(i,k,j) = pf_ls(k)/faptot
            sumpf = sumpf + pfx(i,k,j)
            sumfa = sumfa + fapx(i,k,j)
         end if
      end if

4900 continue  

5900 continue  


      call aerimpactscav (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte, num_chem,       &
              ntype_aer, nsize_aer, ncomp_aer, massptr_aer, dens_aer, numptr_aer,  &
		      maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase,          &
		      volumcen_sect, volumlo_sect, volumhi_sect, &
		      waterptr_aer, dens_water_aer,                &
		      scavimptblvol, scavimptblnum,      &
		      nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow,     &
              dtstepc, t_phy, p_phy, pdel, chem,          &
              isprx, fapx, pfx, pfx_inrain,           &
              dqdt, dotend, qsrflx      )


  if ( .NOT. (config_flags%chem_opt == mozart_mosaic_4bin_aq_kpp) ) then
      call gasrainscav (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte, num_chem,     &
                      config_flags,      &
                      dtstepc,    t_phy,      p_phy,      pdel,  chem,      &
                      isprx,      fapx,       pfx,        pfx_inrain, &
                      dqdt,       dotend,     qsrflx      )
  end if



   do n=1,num_chem
      if(dotend(n))then
         do 6000 j=jts,jte
         do 6000 k=kts,kte
         do 6000 i=its,ite
	    chem(i,k,j,n) = chem(i,k,j,n) + dqdt(i,k,j,n)*dtstepc
 6000    continue
      end if
   end do

   end subroutine wetscav





subroutine aerimpactscav (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,num_chem,  &
     ntype_aer, nsize_aer, ncomp_aer, massptr_aer, dens_aer, numptr_aer, &
     maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase,           &
     volumcen_sect, volumlo_sect, volumhi_sect,                          &
     waterptr_aer, dens_water_aer,                                       &
     scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow, &
     deltat,     t,          pmid,       pdel, chem,                     &
     isprx,      fapx,       pfx,        pfx_inrain,                     &
     dqdt,       dotend,     qsrflx                                      )















  USE module_model_constants, only: g,rhowater, mwdry
  USE module_state_description, only: param_first_scalar

   implicit none









   integer, intent(in)  :: num_chem           
   integer, intent(in)  :: ims,ime            
   integer, intent(in)  :: kms,kme            
   integer, intent(in)  :: jms,jme            
   integer, intent(in)  :: its,ite            
   integer, intent(in)  :: kts,kte            
   integer, intent(in)  :: jts,jte            
   real, intent(in) :: deltat           

   real, intent(in) :: t(ims:ime,kms:kme,jms:jme)    
   real, intent(in) :: pmid(ims:ime,kms:kme,jms:jme) 
   real, intent(in) :: pdel(ims:ime,kms:kme,jms:jme) 
   real, intent(in) :: chem(ims:ime,kms:kme,jms:jme,num_chem) 

   logical, intent(in)  :: isprx(ims:ime,kms:kme,jms:jme) 
   real, intent(in) :: fapx(ims:ime,kms:kme,jms:jme)    
   real, intent(in) :: pfx(ims:ime,kms:kme,jms:jme)     
                                                 
   real, intent(in) :: pfx_inrain(ims:ime,kms:kme,jms:jme)  

   real, intent(inout) :: dqdt(ims:ime,kms:kme,jms:jme,num_chem) 
   logical,  intent(inout) :: dotend(num_chem)     
   real, intent(inout) :: qsrflx(ims:ime,jms:jme,num_chem) 
                        
                        
                        
   integer, intent(in)  :: maxd_atype, maxd_asize, maxd_acomp, maxd_aphase
   integer, intent(in) :: ai_phase
   integer, intent(in) :: ntype_aer
   integer, intent(in) ::  nsize_aer( maxd_atype ),   & 
      	  ncomp_aer( maxd_atype ),   & 
      	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase), & 
      	  waterptr_aer( maxd_asize, maxd_atype ), & 
      	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) 
   real, intent(in) :: volumcen_sect( maxd_asize, maxd_atype ),   & 
                        volumlo_sect( maxd_asize, maxd_atype ),   & 
                        volumhi_sect( maxd_asize, maxd_atype )      
   real, intent(in) :: dens_aer( maxd_acomp, maxd_atype ),   & 
                       dens_water_aer 

   real, intent(in) :: dlndg_nimptblgrow
   integer, intent(in) :: nimptblgrow_mind, nimptblgrow_maxd
   real, intent(in) :: scavimptblnum(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype), &
     	     scavimptblvol(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)



   integer :: i,j               
   integer :: jgrow, jp       
   integer :: k               
   integer :: l, ll, m, n        

   logical :: ispr_anywhere

   real :: dry_mass, dry_volu, dry_volu_1p
   real :: duma
   real :: dumfhi, dumflo
   real :: dumimpactamt0, dumimpactamt3, dumimpactamtsum
   real :: dumimpactratea0, dumimpactratea3
   real :: dumimpactrateb0, dumimpactrateb3
   real :: dumdgratio, dumrate, dumwetdens
   real :: dumlogdens, dumlogptot, dumlogtemp
   real :: dumnumb
   real :: dumscavratenum, dumscavratevol
   real :: pdel_dt_fac
   real :: pf_to_prmmh
   real :: scavimptbl1, scavimptbl2, scavimptbl3, scavimptbl4
   real :: xgrow
   real :: wet_mass, wet_volu
   real, save :: third = 0.33333333









   ispr_anywhere = .false.

   do 5900 i = its,ite
   do 5900 j = jts,jte

      do 4900 k = kte,kts,-1


      if ( isprx(i,k,j) ) then
        ispr_anywhere = .true.
      else
        goto 4900
      end if


	dumlogtemp = log( t(i,k,j) )

	dumlogptot = log( 10.0*pmid(i,k,j) )


       do 3900 n=1,ntype_aer
       do 3900 m=1,nsize_aer(n)
          dry_volu = 0.
          dry_mass = 0
          do ll = 1, ncomp_aer(n)
             l = massptr_aer(ll,m,n,ai_phase)
             if (l < param_first_scalar) cycle
             duma = max( chem(i,k,j,l)*1.0e-6, 0.0 )
             dry_mass = dry_mass + duma                 
             dry_volu = dry_volu + duma/dens_aer(ll,n)  
          end do


          if( dry_volu < 1.0e-26 ) goto 3900


          l = waterptr_aer(m,n)
          if (l >= param_first_scalar) then
             duma = max( chem(i,k,j,l)*1.0e-6, 0.0 )
          else
             duma = 0.0
          end if
          wet_mass = dry_mass + duma
          wet_volu = dry_volu + duma/dens_water_aer


          dumnumb = chem(i,k,j,numptr_aer(m,n,ai_phase))
          if (dry_volu > volumhi_sect(m,n)*dumnumb) then
             dry_volu_1p = volumhi_sect(m,n)
          else if (dry_volu < volumlo_sect(m,n)*dumnumb) then
             dry_volu_1p = volumlo_sect(m,n)
          else
             dry_volu_1p = dry_volu/dumnumb
          end if



          if (wet_volu <= 1.0e9*dry_volu) then
             duma = wet_volu/dry_volu
          else
             duma = 1.0e9
          end if
          dumdgratio = ((dry_volu_1p/volumcen_sect(m,n))*duma)**third


	  dumwetdens = wet_mass/wet_volu
	  dumlogdens = log( dumwetdens )
	  dumimpactamt3 = 0
	  dumimpactamt0 = 0




	  if ((dumdgratio .ge. 0.99) .and. (dumdgratio .le. 1.01)) then
	    scavimptbl1 = scavimptblvol(1,0,m,n)
	    scavimptbl2 = scavimptblvol(2,0,m,n)
	    scavimptbl3 = scavimptblvol(3,0,m,n)
	    scavimptbl4 = scavimptblvol(4,0,m,n)

	  else
	    xgrow = log( dumdgratio ) / dlndg_nimptblgrow
	    jgrow = int( xgrow )
	    if (xgrow .lt. 0.) jgrow = jgrow - 1
	    if (jgrow .lt. nimptblgrow_mind) then
		jgrow = nimptblgrow_mind
		xgrow = jgrow
	    else 
		jgrow = min( jgrow, nimptblgrow_maxd-1 )
	    end if

	    dumfhi = xgrow - jgrow
	    dumflo = 1. - dumfhi

	    scavimptbl1 = dumflo*scavimptblvol(1,jgrow,m,n) +   &
			  dumfhi*scavimptblvol(1,jgrow+1,m,n)

	    scavimptbl2 = dumflo*scavimptblvol(2,jgrow,m,n) +   &
			  dumfhi*scavimptblvol(2,jgrow+1,m,n)

	    scavimptbl3 = dumflo*scavimptblvol(3,jgrow,m,n) +   &
			  dumfhi*scavimptblvol(3,jgrow+1,m,n)

	    scavimptbl4 = dumflo*scavimptblvol(4,jgrow,m,n) +   &
			  dumfhi*scavimptblvol(4,jgrow+1,m,n)
	  end if


	  dumimpactratea3 = exp( scavimptbl1 + scavimptbl2*dumlogtemp +    &
      		scavimptbl3*dumlogptot + scavimptbl4*dumlogdens )








	dumimpactamtsum = 0.0
	    dumimpactrateb3 = dumimpactratea3 * pfx_inrain(i,k,j)
	    dumimpactamt3 =  (1. - exp(-deltat*dumimpactrateb3)) * fapx(i,k,j)
	    dumimpactamtsum = dumimpactamtsum + dumimpactamt3
	dumimpactamt3 = min( dumimpactamtsum, 1.0 )





























	if (numptr_aer(m,n,ai_phase) < param_first_scalar) then
	    dumimpactamt0 = dumimpactamt3
	    goto 3700
	end if


	if ((dumdgratio .ge. 0.99) .and. (dumdgratio .le. 1.01)) then
	    scavimptbl1 = scavimptblnum(1,0,m,n)
	    scavimptbl2 = scavimptblnum(2,0,m,n)
	    scavimptbl3 = scavimptblnum(3,0,m,n)
	    scavimptbl4 = scavimptblnum(4,0,m,n)

	else
	    scavimptbl1 = dumflo*scavimptblnum(1,jgrow,m,n) +   &
			  dumfhi*scavimptblnum(1,jgrow+1,m,n)

	    scavimptbl2 = dumflo*scavimptblnum(2,jgrow,m,n) +   &
			  dumfhi*scavimptblnum(2,jgrow+1,m,n)

	    scavimptbl3 = dumflo*scavimptblnum(3,jgrow,m,n) +   &
			  dumfhi*scavimptblnum(3,jgrow+1,m,n)

	    scavimptbl4 = dumflo*scavimptblnum(4,jgrow,m,n) +   &
			  dumfhi*scavimptblnum(4,jgrow+1,m,n)
	end if


	dumimpactratea0 = exp( scavimptbl1 + scavimptbl2*dumlogtemp +    &
      		scavimptbl3*dumlogptot + scavimptbl4*dumlogdens )

	dumimpactamt0 = 0.0
	    dumimpactrateb0 = dumimpactratea0 * pfx_inrain(i,k,j)
	    dumimpactamt0 = dumimpactamt0 +   &
		(1. - exp( -deltat*dumimpactrateb0 )) * fapx(i,k,j)
	dumimpactamt0 = min( dumimpactamt0, 1.0 )










3700	continue




	pdel_dt_fac = deltat*pdel(i,k,j)/g
	dumrate = -dumimpactamt3/(deltat*(1.0 + 1.0e-8))
	dumrate = min(0.0,max(-1.0/deltat,dumrate)) 
	do ll = 1, ncomp_aer(n)
	    l = massptr_aer(ll,m,n,ai_phase)
	    if (l < param_first_scalar) cycle
	    dqdt(i,k,j,l) = chem(i,k,j,l)*dumrate
	    qsrflx(i,j,l) = qsrflx(i,j,l) + dqdt(i,k,j,l)*pdel_dt_fac 
	end do
	l = waterptr_aer(m,n)
	if (l >= param_first_scalar) then 
	    dqdt(i,k,j,l) = chem(i,k,j,l)*dumrate
	    qsrflx(i,j,l) = qsrflx(i,j,l) + dqdt(i,k,j,l)*pdel_dt_fac 
	end if
	l = numptr_aer(m,n,ai_phase)
	if (l >= param_first_scalar) then 
	    dumrate = -dumimpactamt0/(deltat*(1.0 + 1.0e-8))
	    dqdt(i,k,j,l) = chem(i,k,j,l)*dumrate
	    qsrflx(i,j,l) = qsrflx(i,j,l) + dqdt(i,k,j,l)*pdel_dt_fac 
	end if



3900 continue  

4900 continue  

5900 continue  



   if ( ispr_anywhere ) then
      do n=1,ntype_aer
      do m=1,nsize_aer(n)
         do ll = 1, ncomp_aer(n)
            if (massptr_aer(ll,m,n,ai_phase) >= param_first_scalar) &
               dotend(massptr_aer(ll,m,n,ai_phase)) = .true.
         end do
         if (waterptr_aer(m,n) >= param_first_scalar) &
            dotend(waterptr_aer(m,n)) = .true.
         if (numptr_aer(m,n,ai_phase) >= param_first_scalar) &
            dotend(numptr_aer(m,n,ai_phase)) = .true.
      end do
      end do
   end if


   return
end subroutine aerimpactscav





  subroutine initwet( ntype_aer, nsize_aer, ncomp_aer, massptr_aer,      &
       dens_aer, numptr_aer, maxd_acomp, maxd_asize,maxd_atype,          &
       maxd_aphase, dcen_sect, sigmag_aer, waterptr_aer, dens_water_aer, &
       scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, &
       dlndg_nimptblgrow)








  implicit none

   integer, intent(in) :: maxd_atype, maxd_asize, maxd_acomp, maxd_aphase
   integer, intent(in) :: ntype_aer
   integer, intent(in) ::  nsize_aer( maxd_atype ) 
   integer, intent(in) :: ncomp_aer( maxd_atype ) 
   integer, intent(in) :: massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase) 
   integer, intent(in) :: waterptr_aer( maxd_asize, maxd_atype ) 
   integer, intent(in) :: numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) 
   real, intent(in) :: dens_aer( maxd_acomp, maxd_atype ) 
   real, intent(in) :: dens_water_aer 
   real, intent(in) :: sigmag_aer(maxd_asize, maxd_atype)
   real, intent(in) :: dcen_sect( maxd_asize, maxd_atype ) 
       
       
       
       
       

   real, intent(out) :: dlndg_nimptblgrow
   integer, intent(in) :: nimptblgrow_mind, nimptblgrow_maxd
   real, intent(out) :: scavimptblnum(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)
   real, intent(out) :: scavimptblvol(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)


	integer nnfit_maxd
	parameter (nnfit_maxd=27)

	integer i, j, m, n, jgrow, jdens, jpress, jtemp, ll, mode, nnfit
        integer lunerr

	real dg0, dg0_base, press, rhodryaero, rhowetaero, rmserr, &
     	scavratenum, scavratevol, sigmag,                &
     	temp, wetdiaratio, wetvolratio
	real aafitnum(4), xxfitnum(4,nnfit_maxd), yyfitnum(nnfit_maxd)
	real aafitvol(4), xxfitvol(4,nnfit_maxd), yyfitvol(nnfit_maxd)


	lunerr = 6
	dlndg_nimptblgrow = log( 1.25d00 )

        do 7900 n=1,ntype_aer
        do 7900 m=1,nsize_aer(n)

	sigmag = sigmag_aer(m,n)
	dg0_base = dcen_sect(m,n)*exp( -1.5*((log(sigmag))**2) )



	rhodryaero = dens_aer(1,n)

	do 7800 jgrow = nimptblgrow_mind, nimptblgrow_maxd

	wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
	dg0 = dg0_base*wetdiaratio

	wetvolratio = exp( jgrow*dlndg_nimptblgrow*3. )
	rhowetaero = 1.0 + (rhodryaero-1.0)/wetvolratio
	rhowetaero = min( rhowetaero, rhodryaero )




	nnfit = 0

	do 5900 jtemp = -1, 1
	temp = 273.16 + 25.*jtemp

	do 5900 jpress = -1, 1
	press = 0.75e6 + 0.25e6*jpress

	do 5900 jdens = 0, 2
	rhowetaero = rhodryaero**(0.5*jdens)

	call calc_1_impact_rate( &
     		dg0, sigmag, rhowetaero, temp, press, &
     		scavratenum, scavratevol, lunerr ) 

	nnfit = nnfit + 1
	if (nnfit .gt. nnfit_maxd) then
	    call wrf_error_fatal3("<stdin>",942,&
'*** subr. aerosol_impact_setup -- nnfit too big')
	end if

	xxfitnum(1,nnfit) = 1.
	xxfitnum(2,nnfit) = log( temp )
	xxfitnum(3,nnfit) = log( press )
	xxfitnum(4,nnfit) = log( rhowetaero )
	yyfitnum(nnfit) = log( scavratenum )

	xxfitvol(1,nnfit) = 1.
	xxfitvol(2,nnfit) = log( temp )
	xxfitvol(3,nnfit) = log( press )
	xxfitvol(4,nnfit) = log( rhowetaero )
	yyfitvol(nnfit) = log( scavratevol )

5900	continue





	call mlinft( xxfitnum, yyfitnum, aafitnum, nnfit, 4, 4, rmserr )
	call mlinft( xxfitvol, yyfitvol, aafitvol, nnfit, 4, 4, rmserr )

	do i = 1, 4
	    scavimptblnum(i,jgrow,m,n) = aafitnum(i)
	    scavimptblvol(i,jgrow,m,n) = aafitvol(i)
	end do


7800	continue
7900	continue

        return
      end subroutine initwet





	subroutine calc_1_impact_rate(             &
     		dg0, sigmag, rhoaero, temp, press, &
     		scavratenum, scavratevol, lunerr )













	implicit none


	integer lunerr
	real dg0, sigmag, rhoaero, temp, press, scavratenum, scavratevol


	integer nrainsvmax
	parameter (nrainsvmax=50)
	real rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
     		vfallrainsv(nrainsvmax)

	integer naerosvmax
	parameter (naerosvmax=51)
	real aaerosv(naerosvmax), &
     	ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

	integer i, ja, jr, na, nr
	real a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
     	real anumsum, avolsum, boltzmann, cair, chi
     	real d, dr, dum, dumfuchs, dx
     	real ebrown, eimpact, eintercept, etotal, freepath, gravity
     	real pi, precip, precipmmhr, precipsum
     	real r, rainsweepout, reynolds, rhi, rhoair, rhowater, rlo, rnumsum
     	real scavsumnum, scavsumnumbb
     	real scavsumvol, scavsumvolbb
     	real schmidt, sqrtreynolds, sstar, stokes, sx
     	real taurelax, vfall, vfallstp
     	real x, xg0, xg3, xhi, xlo, xmuwaterair


	rlo = .005
	rhi = .250
	dr = 0.005
	nr = 1 + nint( (rhi-rlo)/dr )
	if (nr .gt. nrainsvmax) then
	    call wrf_error_fatal3("<stdin>",1035,&
'*** subr. calc_1_impact_rate -- nr > nrainsvmax')
	end if

	precipmmhr = 1.0
	precip = precipmmhr/36000. 

	ag0 = dg0/2.
	sx = log( sigmag )
	xg0 = log( ag0 )
	xg3 = xg0 + 3.*sx*sx

	xlo = xg3 - 4.*sx
	xhi = xg3 + 4.*sx
	dx = 0.2*sx

	dx = max( 0.2*sx, 0.01 )
	xlo = xg3 - max( 4.*sx, 2.*dx )
	xhi = xg3 + max( 4.*sx, 2.*dx )

	na = 1 + nint( (xhi-xlo)/dx )
	if (na .gt. naerosvmax) then
	    call wrf_error_fatal3("<stdin>",1057,&
'*** subr. calc_1_impact_rate -- na > naerosvmax')
	end if


	cair = press/(8.31436e7*temp)

	rhoair = 28.966*cair

	freepath = 2.8052e-10/cair

	boltzmann = 1.3807e-16

	rhowater = 1.0

	gravity = 980.616

	airdynvisc = 1.8325e-4 * (416.16/(temp+120.)) *    &
                                        ((temp/296.16)**1.5)

	airkinvisc = airdynvisc/rhoair

	xmuwaterair = 60.0

	pi = 3.1415926536








	precipsum = 0.
	do i = 1, nr
	    r = rlo + (i-1)*dr
	    rrainsv(i) = r
	    xnumrainsv(i) = exp( -r/2.7e-2 )

	    d = 2.*r
	    if (d .le. 0.007) then
		vfallstp = 2.88e5 * d**2.
	    else if (d .le. 0.025) then
		vfallstp = 2.8008e4 * d**1.528
	    else if (d .le. 0.1) then
		vfallstp = 4104.9 * d**1.008
	    else if (d .le. 0.25) then
		vfallstp = 1812.1 * d**0.638
	    else
		vfallstp = 1069.8 * d**0.235
	    end if

	    vfall = vfallstp * sqrt(1.204e-3/rhoair)
	    vfallrainsv(i) = vfall
	    precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
	end do
	precipsum = precipsum*pi*1.333333

	rnumsum = 0.
	do i = 1, nr
	    xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
	    rnumsum = rnumsum + xnumrainsv(i)
	end do







	anumsum = 0.
	avolsum = 0.
	do i = 1, na
	    x = xlo + (i-1)*dx
	    a = exp( x )
	    aaerosv(i) = a
	    dum = (x - xg0)/sx
	    ynumaerosv(i) = exp( -0.5*dum*dum )
	    yvolaerosv(i) = ynumaerosv(i)*1.3333*pi*a*a*a
	    anumsum = anumsum + ynumaerosv(i)
	    avolsum = avolsum + yvolaerosv(i)
	end do

	do i = 1, na
	    ynumaerosv(i) = ynumaerosv(i)/anumsum
	    yvolaerosv(i) = yvolaerosv(i)/avolsum
	end do





	scavsumnum = 0.
	scavsumvol = 0.



	do 5900 jr = 1, nr

	r = rrainsv(jr)
	vfall = vfallrainsv(jr)

	reynolds = r * vfall / airkinvisc
	sqrtreynolds = sqrt( reynolds )




	scavsumnumbb = 0.
	scavsumvolbb = 0.

	do 5500 ja = 1, na

	a = aaerosv(ja)

	chi = a/r

	dum = freepath/a
	dumfuchs = 1. + 1.246*dum + 0.42*dum*exp(-0.87/dum)
	taurelax = 2.*rhoaero*a*a*dumfuchs/(9.*rhoair*airkinvisc)

	aeromass = 4.*pi*a*a*a*rhoaero/3.
	aerodiffus = boltzmann*temp*taurelax/aeromass

	schmidt = airkinvisc/aerodiffus
	stokes = vfall*taurelax/r

	ebrown = 4.*(1. + 0.4*sqrtreynolds*(schmidt**0.3333333)) /  &
     			(reynolds*schmidt)

	dum = (1. + 2.*xmuwaterair*chi) /         &
     			(1. + xmuwaterair/sqrtreynolds)
	eintercept = 4.*chi*(chi + dum)

	dum = log( 1. + reynolds )
	sstar = (1.2 + dum/12.) / (1. + dum)
	eimpact = 0.
	if (stokes .gt. sstar) then
	    dum = stokes - sstar
	    eimpact = (dum/(dum+0.6666667)) ** 1.5
	end if

	etotal = ebrown + eintercept + eimpact
	etotal = min( etotal, 1.0 )

	rainsweepout = xnumrainsv(jr)*4.*pi*r*r*vfall

	scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
	scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

5500	continue

	scavsumnum = scavsumnum + scavsumnumbb
	scavsumvol = scavsumvol + scavsumvolbb
5900	continue

	scavratenum = scavsumnum*3600.
	scavratevol = scavsumvol*3600.



  return
  end subroutine calc_1_impact_rate





subroutine gasrainscav (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,num_chem,  &
                      config_flags,   &
                      deltat,     t,          pmid,       pdel,       &
                      chem,                                  &
                      isprx,      fapx,       pfx,        pfx_inrain, &
                      dqdt,       dotend,     qsrflx      )















  USE module_model_constants, only: g,rhowater, mwdry
  use module_configure, only:  grid_config_rec_type,   &
		param_first_scalar,   &
		p_so2, p_h2o2, p_sulf,p_h2so4, p_msa,   &
		p_hno3, p_hcl, p_nh3

   implicit none











   integer, intent(in)  :: num_chem           
   integer, intent(in)  :: ims,ime            
   integer, intent(in)  :: kms,kme            
   integer, intent(in)  :: jms,jme            
   integer, intent(in)  :: its,ite            
   integer, intent(in)  :: kts,kte            
   integer, intent(in)  :: jts,jte            
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   real, intent(in) :: deltat           

   real, intent(in) :: t(ims:ime,kms:kme,jms:jme)    
   real, intent(in) :: pmid(ims:ime,kms:kme,jms:jme) 
   real, intent(in) :: pdel(ims:ime,kms:kme,jms:jme) 
   real, intent(in) :: chem(ims:ime,kms:kme,jms:jme,num_chem) 

   logical, intent(in)  :: isprx(ims:ime,kms:kme,jms:jme) 
   real, intent(in) :: fapx(ims:ime,kms:kme,jms:jme)    
   real, intent(in) :: pfx(ims:ime,kms:kme,jms:jme)     
                                                 
   real, intent(in) :: pfx_inrain(ims:ime,kms:kme,jms:jme)  

   real, intent(out) :: dqdt(ims:ime,kms:kme,jms:jme,num_chem) 
   logical,  intent(inout) :: dotend(num_chem)     
   real, intent(inout) :: qsrflx(ims:ime,jms:jme,num_chem)
                        
                        
                        



   integer :: i, j, k         
   integer :: jp              
   integer :: p1st
   logical :: ispr_anywhere

   integer, parameter :: ng = 7
   integer, parameter :: ig_so2   = 1
   integer, parameter :: ig_h2o2  = 2
   integer, parameter :: ig_h2so4 = 3
   integer, parameter :: ig_msa   = 4
   integer, parameter :: ig_hno3  = 5
   integer, parameter :: ig_hcl   = 6
   integer, parameter :: ig_nh3   = 7
   integer :: ig, lg, lg_ptr(ng)

   real :: amtscav(ng), amtscav_sub(ng)
   real :: fracgas(ng)
   real :: fracscav(ng), fracscav_sub(ng)
   real :: deltatinv
   real :: dum, dumamt, dumprecipmmh, dumpress, dumtemp
   real :: pdel_fac
   real :: r_gc(ng)
   real :: scavrate_hno3
   real :: scavrate(ng), scavrate_factor(ng)









   ispr_anywhere = .false.
   deltatinv = 1.0/(deltat*(1.0d0 + 1.0d-15))

   p1st = param_first_scalar

   lg_ptr(ig_so2  ) = p_so2  
   lg_ptr(ig_h2o2 ) = p_h2o2 
   lg_ptr(ig_h2so4) = p_sulf
   lg_ptr(ig_h2so4) = p_h2so4
   lg_ptr(ig_msa  ) = p_msa  
   lg_ptr(ig_hno3 ) = p_hno3 
   lg_ptr(ig_hcl  ) = p_hcl  
   lg_ptr(ig_nh3  ) = p_nh3  

   scavrate_factor(ig_so2  ) = 1.08
   scavrate_factor(ig_h2o2 ) = 1.38
   scavrate_factor(ig_h2so4) = 0.80
   scavrate_factor(ig_msa  ) = 0.80
   scavrate_factor(ig_hno3 ) = 1.00
   scavrate_factor(ig_hcl  ) = 1.15
   scavrate_factor(ig_nh3  ) = 1.59


   do 5900 j = jts,jte
   do 5900 i = its,ite

      do 4900 k = kts,kte


      if ( isprx(i,k,j)  ) then
        ispr_anywhere = .true.
      else
        goto 4900
      end if



	dumtemp = t(i,k,j)
	if (dumtemp .le. 273.16) goto 4900
	dumpress = 10.0*pmid(i,k,j)

	do ig = 1, ng
	    fracscav(ig) = 0.0
	    fracgas(ig) = 1.0
	    lg = lg_ptr(ig)
	    if (lg .ge. p1st) then
		r_gc(ig) = max( chem(i,k,j,lg), 0.0 )



	    else
		r_gc(ig) = 0.0
	    end if
	end do
		
	if ( .not. isprx(i,k,j) ) goto 3600


	dumprecipmmh = pfx_inrain(i,k,j)*3600.0



	scavrate_hno3 = 6.262e-5*(dumprecipmmh**0.7366)   &
      		* ((dumtemp/298.0)**1.12)   &
      		* ((1.013e6/dumpress)**.75)

	do ig = 1, ng
	    scavrate(ig) = scavrate_hno3*scavrate_factor(ig)
	    fracscav_sub(ig) = (1. - exp(-scavrate(ig)*deltat))   &
			*fracgas(ig)*fapx(i,k,j)
	    amtscav_sub(ig) = r_gc(ig)*min( fracscav_sub(ig), 1.0 )
	end do



	dumamt = min( amtscav_sub(ig_so2), amtscav_sub(ig_h2o2) )
	fracscav_sub(ig_so2 ) = dumamt/max( r_gc(ig_so2 ), 1.0e-30 )
	fracscav_sub(ig_h2o2) = dumamt/max( r_gc(ig_h2o2), 1.0e-30 )
	amtscav_sub(ig_so2 ) = r_gc(ig_so2 )*min( fracscav_sub(ig_so2 ), 1.0 )
	amtscav_sub(ig_h2o2) = r_gc(ig_h2o2)*min( fracscav_sub(ig_h2o2), 1.0 )


	dumamt = 2.0*amtscav_sub(ig_so2)   &
		+ 2.0*amtscav_sub(ig_h2so4) + amtscav_sub(ig_msa)   &
		+ amtscav_sub(ig_hno3) + amtscav_sub(ig_hcl)
	dumamt = min( dumamt, amtscav_sub(ig_nh3) )
	fracscav_sub(ig_nh3) = dumamt/max( r_gc(ig_nh3), 1.0e-30 )
	amtscav_sub(ig_nh3 ) = r_gc(ig_nh3 )*min( fracscav_sub(ig_nh3 ), 1.0 )

	do ig = 1, ng
	    fracscav(ig) = fracscav(ig) + fracscav_sub(ig)
	end do






















3600	continue




	pdel_fac = (pdel(i,k,j)/(g*mwdry))

	do ig = 1, ng
	    fracscav(ig) = max(0.0,min(1.0,fracscav(ig))) 
	    amtscav(ig)  = fracscav(ig)*r_gc(ig)
	    lg = lg_ptr(ig)
	    if (lg .ge. p1st) then
		dqdt(i,k,j,lg) = -deltatinv*amtscav(ig)  
		qsrflx(i,j,lg) = qsrflx(i,j,lg) - pdel_fac*amtscav(ig) 
	    end if
	end do

4900 continue  

5900 continue  



   if ( ispr_anywhere ) then
       do ig = 1, ng
           if (lg_ptr(ig) .ge. p1st) dotend(lg_ptr(ig)) = .true.
       end do
   end if


   return
end subroutine gasrainscav





	subroutine mlinft( x, y, a, n, m, mmaxd, rmserr )
















	implicit none


	integer n, m, mmaxd
	real x(mmaxd,n), y(n), a(mmaxd), rmserr


	integer i, j, jflag, k
	real aa(10,10), bb(10), errsq, resid, ydum

	if (n .le. 1) then
	    a(1) = 1.e30
	    rmserr = 0.
	    return
	end if

	do 2900 i = 1, m
	    do 2100 j = 1, m
		aa(i,j) = 0.0
2100	    continue
	    bb(i) = 0.0

	    do 2500 k = 1, n
		do 2300 j = 1, m
		    aa(i,j) = aa(i,j) + x(i,k)*x(j,k)
2300		continue
		bb(i) = bb(i) + x(i,k)*y(k)
2500	    continue

2900	continue







	call linsolv( aa, a, bb, m, 10, 10, jflag )


	errsq = 0.
	do 3300 k = 1, n
	    ydum = 0.0
	    do 3100 i = 1, m
		ydum = ydum + a(i)*x(i,k)
3100	    continue
	    resid = ydum - y(k)
	    errsq = errsq + resid*resid
3300	continue
	rmserr = sqrt( errsq/n )

	return
	end subroutine mlinft





	subroutine linsolv( a, x, b, n, m1, m2, jflag )


















	implicit none


	integer n, m1, m2, jflag
	real a(m1,m2), x(n), b(n)


	integer i, imax, iup, j, k
	real amax, asmall, dmy, rsmall
	parameter (rsmall = 1.0e-16)

	jflag = 0




	do 1900 k = 1, n




	    imax = k
	    amax = abs( a(imax,k) )
	    do 1200 i = k+1, n
		if (abs(a(i,k)) .gt. amax) then
		    imax = i
		    amax = abs(a(i,k))
		end if
1200	    continue
	    if (amax .eq. 0.) return

	    if (imax .ne. k) then
		do 1400 j = k, n
		    dmy = a(imax,j)
		    a(imax,j) = a(k,j)
		    a(k,j) = dmy
1400		continue
		dmy = b(imax)
		b(imax) = b(k)
		b(k) = dmy
	    end if




	    asmall = abs(a(k,k))
	    do 1700 i = k+1, n
		if (a(i,k) .ne. 0.0) then
		    if (asmall .le. abs(rsmall*a(i,k))) return
		    dmy = a(i,k)/a(k,k)
		    a(i,k) = 0.0
		    do 1600 j = k+1, n
			a(i,j) = a(i,j) - dmy*a(k,j)
1600		    continue
		    b(i) = b(i) - dmy*b(k)
		end if
1700	    continue

1900	continue




	do 2900 iup = 1, n
	    i = n + 1 - iup
	    dmy = b(i)
	    do 2500 j = i+1, n
		dmy = dmy - a(i,j)*x(j)
2500	    continue
	    asmall = abs(a(i,i))
	    if (abs(a(i,i)) .le. abs(rsmall*dmy)) return
	    x(i) = dmy/a(i,i)
2900	continue

	jflag = 1

	return
	end subroutine linsolv

END MODULE module_mosaic_wetscav
