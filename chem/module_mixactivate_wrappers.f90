














MODULE module_mixactivate_wrappers

CONTAINS



      subroutine mosaic_mixactivate (                        &
           id, ktau, dtstep, config_flags, idrydep_onoff,    &
           rho_phy, t_phy, w, cldfra, cldfra_old,            &
           ddvel, z, dz8w, p_at_w, t_at_w, exch_h,           &
           qv, qc, qi, qndrop3d, f_qc, f_qi, chem,           &
	       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
	       qsrflx, &
           ids,ide, jds,jde, kds,kde,                        &
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte                         )

    USE module_configure, only: grid_config_rec_type
	use module_state_description, only:  num_chem
	use module_data_mosaic_asect
	use module_mixactivate, only: mixactivate



	implicit none


	integer, intent(in) ::               &
         id, ktau,                       &
         ids, ide, jds, jde, kds, kde,   &
         ims, ime, jms, jme, kms, kme,   &
         its, ite, jts, jte, kts, kte,   &
         idrydep_onoff

	real, intent(in) :: dtstep

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		rho_phy, t_phy, w,   &
		z, dz8w, p_at_w, t_at_w, exch_h

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: cldfra, cldfra_old

	real, intent(in),   &
		dimension( its:ite, jts:jte, num_chem ) :: ddvel

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		qv, qc, qi

    LOGICAL, intent(in) :: f_qc, f_qi

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		qndrop3d

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem
      real, intent(out), dimension(ims:ime,kms:kme,jms:jme) :: nsource,&
	     ccn1,ccn2,ccn3,ccn4,ccn5,ccn6  

	type(grid_config_rec_type), intent(in) :: config_flags
 real, intent(out) :: qsrflx(ims:ime, jms:jme, num_chem) 


	real sumhygro,sumvol
	integer i,j,k,l,m,n
	real hygro( its:ite, kts:kte, jts:jte, maxd_asize, maxd_atype ) 

      qsrflx(:,:,:) = 0.0


      do 100 j=jts,jte
      do 100 k=kts,kte
      do 100 i=its,ite
       do n=1,ntype_aer
       do m=1,nsize_aer(n)
	       sumhygro=0.
	       sumvol=0.
	       do l=1,ncomp_aer(n)
	          sumhygro = sumhygro+hygro_aer(l,n)*   &
                   chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	          sumvol = sumvol+chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	       end do 
           hygro(i,k,j,m,n)=sumhygro/sumvol
	end do 
	end do 
  100 continue



      call mixactivate(  msectional, &
           chem, num_chem, qv, qc, qi, qndrop3d,   &
           t_phy, w, ddvel, idrydep_onoff,  &
           maxd_acomp, maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, dlo_sect, dhi_sect, sigmag_aer, dcen_sect,  &
           dens_aer, mw_aer,           &
           waterptr_aer, hygro,  ai_phase, cw_phase,                &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old, qsrflx, &
	       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
	       id, ktau, dtstep, &
           f_qc, f_qi               )

      end subroutine mosaic_mixactivate





      subroutine mosaic_mixactivate_init(                    &
           config_flags, chem, scalar,                       &
           chem_in_opt,                                      & 
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte                         )

      USE module_configure, only: grid_config_rec_type
      use module_state_description, only:  num_chem, num_scalar, p_qndrop
      use module_data_mosaic_asect

	implicit none


      type(grid_config_rec_type), intent(in) :: config_flags

      integer, intent(in) ::               &
           ims, ime, jms, jme, kms, kme,   &
           its, ite, jts, jte, kts, kte
      INTEGER,      INTENT(IN   ) :: chem_in_opt 
      real, intent(inout),   &
           dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
           chem

      real, intent(inout),   &
           dimension( ims:ime, kms:kme, jms:jme, 1:num_scalar ) :: &
           scalar

      integer :: i, j, k, m, n, l

      do j=jts,jte
         do k=kts,kte
            do i=its,ite
               scalar(i,k,j,p_qndrop) = 0.               
            end do
         end do
      end do

      if( cw_phase > 0 ) then   
                                
         if (config_flags%chem_in_opt == 1) then 
            do n=1,ntype_aer
               do m=1,nsize_aer(n)
                  chem(its:ite,kts:kte,jts:jte,numptr_aer(m,n,cw_phase)) = 0.
                  do l=1,ncomp_aer(n)
                     if( ai_phase > 0 ) then
                        
                        chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,ai_phase))= &
                             chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,ai_phase)) + &
                             chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,cw_phase))
                        
                     endif 
                     chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,cw_phase)) = 0.
                  end do              
               end do                 
            end do         
         else
            do n=1,ntype_aer
               do m=1,nsize_aer(n)
                  chem(its:ite,kts:kte,jts:jte,numptr_aer(m,n,cw_phase)) = 0.
                  do l=1,ncomp_aer(n)
                     chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,cw_phase)) = 0.
                  end do              
               end do                 
            end do                 
         endif 
      end if

      end subroutine mosaic_mixactivate_init






      subroutine sorgam_mixactivate (                        &
           id, ktau, dtstep, config_flags, idrydep_onoff,    &
           rho_phy, t_phy, w, cldfra, cldfra_old,            &
           ddvel, z, dz8w, p_at_w, t_at_w, exch_h,           &
           qv, qc, qi, qndrop3d, f_qc, f_qi, chem,           &
	       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
           ids,ide, jds,jde, kds,kde,                        &
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte                         )

    USE module_configure, only: grid_config_rec_type
	use module_state_description, only:  num_chem
	use module_data_sorgam
	use module_mixactivate, only: mixactivate



	implicit none


	integer, intent(in) ::                  &
		id, ktau,                       &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte,   &
                idrydep_onoff

	real, intent(in) :: dtstep

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		rho_phy, t_phy, w,   &
		z, dz8w, p_at_w, t_at_w, exch_h

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: cldfra, cldfra_old

	real, intent(in),   &
		dimension( its:ite, jts:jte, num_chem ) :: ddvel

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		qv, qc, qi

    LOGICAL, intent(in) :: f_qc, f_qi

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme  ) :: &
		qndrop3d

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem
      real, intent(out), dimension(ims:ime,kms:kme,jms:jme) :: nsource, &
	     ccn1,ccn2,ccn3,ccn4,ccn5,ccn6  

	type(grid_config_rec_type), intent(in) :: config_flags


	real qsrflx(ims:ime, jms:jme, num_chem) 
	real sumhygro,sumvol
	integer i,j,k,l,m,n
	real hygro( its:ite, kts:kte, jts:jte,maxd_asize, maxd_atype )



      do 100 j=jts,jte
      do 100 k=kts,kte
      do 100 i=its,ite
       do n=1,ntype_aer
       do m=1,nsize_aer(n)
	       sumhygro=0
	       sumvol=0
	       do l=1,ncomp_aer(n)
	          sumhygro = sumhygro+hygro_aer(l,n)*   &
                   chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	          sumvol = sumvol+chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	       end do 
               hygro(i,k,j,m,n)=sumhygro/sumvol
	end do 
	end do 
  100 continue




      call mixactivate(  msectional, &
           chem, num_chem, qv, qc, qi, qndrop3d,   &
           t_phy, w, ddvel, idrydep_onoff,  &
           maxd_acomp, maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, dlo_sect, dhi_sect, sigmag_aer, dcen_sect,  &
           dens_aer, mw_aer,           &
           waterptr_aer, hygro,  ai_phase, cw_phase,                 &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old, qsrflx,                      &
	       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
	       id, ktau, dtstep, &
           f_qc, f_qi               )

      end subroutine sorgam_mixactivate


      subroutine soa_vbs_mixactivate (                       &
           id, ktau, dtstep, config_flags, idrydep_onoff,    &
           rho_phy, t_phy, w, cldfra, cldfra_old,            &
           ddvel, z, dz8w, p_at_w, t_at_w, exch_h,           &
           qv, qc, qi, qndrop3d, f_qc, f_qi, chem,           &
           ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
           ids,ide, jds,jde, kds,kde,                        &
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte                         )

    USE module_configure, only: grid_config_rec_type
        use module_state_description, only:  num_chem
        use module_data_soa_vbs
        use module_mixactivate, only: mixactivate



        implicit none


        integer, intent(in) ::                  &
                id, ktau,                       &
                ids, ide, jds, jde, kds, kde,   &
                ims, ime, jms, jme, kms, kme,   &
                its, ite, jts, jte, kts, kte,   &
                idrydep_onoff

        real, intent(in) :: dtstep

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                rho_phy, t_phy, w,   &
                z, dz8w, p_at_w, t_at_w, exch_h

        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: cldfra, cldfra_old

        real, intent(in),   &
                dimension( its:ite, jts:jte, num_chem ) :: ddvel

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                qv, qc, qi

    LOGICAL, intent(in) :: f_qc, f_qi

        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme  ) :: &
                qndrop3d

        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
                chem
      real, intent(out), dimension(ims:ime,kms:kme,jms:jme) :: nsource, &
             ccn1,ccn2,ccn3,ccn4,ccn5,ccn6  

        type(grid_config_rec_type), intent(in) :: config_flags 
        real qsrflx(ims:ime, jms:jme, num_chem) 
        real sumhygro,sumvol
        integer i,j,k,l,m,n
        real hygro( its:ite, kts:kte, jts:jte,maxd_asize, maxd_atype )



      do 100 j=jts,jte
      do 100 k=kts,kte
      do 100 i=its,ite
       do n=1,ntype_aer
       do m=1,nsize_aer(n)
               sumhygro=0
               sumvol=0
               do l=1,ncomp_aer(n)
                  sumhygro = sumhygro+hygro_aer(l,n)*   &
                             chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
                  sumvol = sumvol+chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
               end do 
               hygro(i,k,j,m,n)=sumhygro/sumvol
        end do 
        end do 
  100 continue




      call mixactivate(  msectional, &
           chem, num_chem, qv, qc, qi, qndrop3d,   &
           t_phy, w, ddvel, idrydep_onoff,  &
           maxd_acomp, maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, dlo_sect, dhi_sect, sigmag_aer, dcen_sect, &
           dens_aer, mw_aer,           &
           waterptr_aer, hygro,  ai_phase, cw_phase,                 &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old, qsrflx,                      &
               ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
               id, ktau, dtstep, &
           f_qc, f_qi               )

      end subroutine soa_vbs_mixactivate

      subroutine sorgam_vbs_mixactivate (                        &
           id, ktau, dtstep, config_flags, idrydep_onoff,    &
           rho_phy, t_phy, w, cldfra, cldfra_old,            &
           ddvel, z, dz8w, p_at_w, t_at_w, exch_h,           &
           qv, qc, qi, qndrop3d, f_qc, f_qi, chem,           &
           ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
           ids,ide, jds,jde, kds,kde,                        &
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte                         )

    USE module_configure, only: grid_config_rec_type
        use module_state_description, only:  num_chem
        use module_data_sorgam_vbs
        use module_mixactivate, only: mixactivate



        implicit none


        integer, intent(in) ::                  &
                id, ktau,                       &
                ids, ide, jds, jde, kds, kde,   &
                ims, ime, jms, jme, kms, kme,   &
                its, ite, jts, jte, kts, kte,   &
                idrydep_onoff

        real, intent(in) :: dtstep

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                rho_phy, t_phy, w,   &
                z, dz8w, p_at_w, t_at_w, exch_h

        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: cldfra, cldfra_old

        real, intent(in),   &
                dimension( its:ite, jts:jte, num_chem ) :: ddvel

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                qv, qc, qi

    LOGICAL, intent(in) :: f_qc, f_qi

        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme  ) :: &
                qndrop3d

        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
                chem
      real, intent(out), dimension(ims:ime,kms:kme,jms:jme) :: nsource, &
             ccn1,ccn2,ccn3,ccn4,ccn5,ccn6  

        type(grid_config_rec_type), intent(in) :: config_flags


        real qsrflx(ims:ime, jms:jme, num_chem) 
        real sumhygro,sumvol
        integer i,j,k,l,m,n
        real hygro( its:ite, kts:kte, jts:jte,maxd_asize, maxd_atype )



      do 100 j=jts,jte
      do 100 k=kts,kte
      do 100 i=its,ite
       do n=1,ntype_aer
       do m=1,nsize_aer(n)
               sumhygro=0
               sumvol=0
               do l=1,ncomp_aer(n)
                  sumhygro = sumhygro+hygro_aer(l,n)*   &
                   chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
                  sumvol = sumvol+chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
               end do 
               hygro(i,k,j,m,n)=sumhygro/sumvol
        end do 
        end do 
  100 continue




      call mixactivate(  msectional, &
           chem, num_chem, qv, qc, qi, qndrop3d,   &
           t_phy, w, ddvel, idrydep_onoff,  &
           maxd_acomp, maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, dlo_sect, dhi_sect, sigmag_aer, dcen_sect,  &
           dens_aer, mw_aer,           &
           waterptr_aer, hygro,  ai_phase, cw_phase,                 &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old, qsrflx,                      &
               ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
               id, ktau, dtstep, &
           f_qc, f_qi               )

      end subroutine sorgam_vbs_mixactivate

END MODULE module_mixactivate_wrappers
