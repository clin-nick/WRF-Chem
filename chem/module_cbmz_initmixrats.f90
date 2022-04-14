








      module module_cbmz_initmixrats


      use module_peg_util


      integer, parameter :: cbmz_init_wrf_mixrats_flagaa = 1
                            


      contains



      subroutine cbmz_init_wrf_mixrats(   &
               config_flags,              &
               z_at_w, g,                 &
               chem, numgas,              &
               ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte  )

















   USE module_configure, only:  grid_config_rec_type, num_chem, &
	p_o3, p_ald, p_hc3, p_hc5, p_hc8, p_ket, p_oli, p_olt, p_ora2, &
    p_hcl, p_par
   USE module_state_description, only:  param_first_scalar,   &
                        gas_ic_pnnl
   USE module_input_chem_data, only:  bdy_chem_value

   IMPLICIT NONE





   INTEGER, INTENT(IN) :: numgas, &
                          ids,ide, jds,jde, kds,kde, &
                          ims,ime, jms,jme, kms,kme, &
                          its,ite, jts,jte, kts,kte

   REAL, INTENT(IN) :: g


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), &
         INTENT(IN) :: z_at_w


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), &
         INTENT(INOUT) :: chem

   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags



	integer i, j, k, kp1
    real, dimension( its:ite, kts:kte, jts:jte ) :: z

	if (cbmz_init_wrf_mixrats_flagaa <= 0) return




    do j = jts, min(jte,jde-1)
       do k = kts, kte
          kp1 = min(k+1, kte)
          do i = its, min(ite,ide-1)
             z(i,k,j) = (z_at_w(i,k,j)+z_at_w(i,kp1,j))*0.5
          end do
       end do
    end do















	do j = jts, min(jte,jde-1)
	do k = kts, kte
	do i = its, min(ite,ide-1)
       if( z(i,k,j) <= 1000. ) then
          chem(i,k,j,p_hcl) = 0.4*1e-3
       elseif( z(i,k,j) > 1000. &
            .and. z(i,k,j) <= 2500. ) then
          chem(i,k,j,p_hcl) = (0.4*1e-3) + (z(i,k,j)-1000.)* &
               ((0.1*1e-3)-(0.4*1e-3)) / (2500.-1000.)
       else
          chem(i,k,j,p_hcl) = 0.1*1e-3
       end if
    end do
    end do
    end do



















	return
	end subroutine cbmz_init_wrf_mixrats  




	end module module_cbmz_initmixrats




	subroutine bdy_chem_value_cbmz ( chem_bv, z, nch, numgas )












	use module_configure, only:  grid_config_rec_type,   &
	    p_o3, p_ald, p_hc3, p_hc5, p_hc8, p_ket, p_oli,  &
	    p_olt, p_ora2, p_hcl, p_par
	use module_input_chem_data, only:  bdy_chem_value

	implicit none


	REAL,    INTENT(OUT)  :: chem_bv    
	REAL,    INTENT(IN)   :: z          
	INTEGER, INTENT(IN)   :: nch        
    INTEGER, INTENT(IN)   :: numgas     

	real chem_bv_ald, chem_bv_hc3, chem_bv_hc5,   &
	     chem_bv_hc8, chem_bv_ket, chem_bv_oli,   &
	     chem_bv_olt, chem_bv_ora2
	real, parameter :: chem_bv_def = 1.0e-20


    if( nch == p_hcl ) then
       
       
       if( z <= 1000. ) then
          chem_bv = 0.4*1e-3
       elseif( z > 1000. &
            .and. z <= 2500. ) then
          chem_bv = (0.4*1e-3) + (z-1000.)* &
               ((0.1*1e-3)-(0.4*1e-3)) / (2500.-1000.)
       else
          chem_bv = 0.1*1e-3
       end if

    else if( nch == p_par ) then
       call bdy_chem_value( chem_bv_ald,  z, p_ald, numgas )

       
       call bdy_chem_value( chem_bv_hc3,  z, numgas+1+1, numgas )
       call bdy_chem_value( chem_bv_hc5,  z, numgas+1+2, numgas )
       call bdy_chem_value( chem_bv_hc8,  z, numgas+1+3, numgas )

       call bdy_chem_value( chem_bv_ket,  z, p_ket, numgas )
       call bdy_chem_value( chem_bv_oli,  z, p_oli, numgas )
       call bdy_chem_value( chem_bv_olt,  z, p_olt, numgas )
       call bdy_chem_value( chem_bv_ora2, z, p_ora2, numgas )

       chem_bv = 0.4*chem_bv_ald + 2.9*chem_bv_hc3    &
            + 4.8*chem_bv_hc5 + 7.9*chem_bv_hc8       &
            + 0.9*chem_bv_ket + 2.8*chem_bv_oli       &
            + 1.8*chem_bv_olt + 1.0*chem_bv_ora2

    else
       
       chem_bv = chem_bv_def

    end if

	return
	end subroutine bdy_chem_value_cbmz
