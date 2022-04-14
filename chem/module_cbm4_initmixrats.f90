module module_cbm4_initmixrats
  
  use module_peg_util
  
  
  integer, parameter :: cbm4_init_wrf_mixrats_flagaa = 1
  
  
  
contains
  
  
  
  subroutine cbm4_init_wrf_mixrats(   &
       config_flags,              &
       z_at_w, g,                 &
       chem, numgas,              &
       ids,ide, jds,jde, kds,kde, &
       ims,ime, jms,jme, kms,kme, &
       its,ite, jts,jte, kts,kte  )

















    USE module_configure, only:  grid_config_rec_type, num_chem, &
         p_o3, p_ald2, p_par
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

	if (cbm4_init_wrf_mixrats_flagaa <= 0) return




    do j = jts, min(jte,jde-1)
       do k = kts, kte
          kp1 = min(k+1, kte)
          do i = its, min(ite,ide-1)
             z(i,k,j) = (z_at_w(i,k,j)+z_at_w(i,kp1,j))*0.5
          end do
       end do
    end do





	if ( (config_flags%chem_in_opt == 0) .and.   &
	     (config_flags%gas_ic_opt == gas_ic_pnnl) ) then
	    do j = jts, min(jte,jde-1)
	    do k = kts, kte
	    do i = its, min(ite,ide-1)
           call bdy_chem_value( chem(i,k,j,p_o3),z(i,k,j), p_o3, numgas )
	    end do
	    end do
	    end do
	end if





















	return
      end subroutine cbm4_init_wrf_mixrats







	subroutine bdy_chem_value_cbm4 ( chem_bv, z, nch, numgas )












	use module_configure, only:  grid_config_rec_type,   &
	    p_par, p_ald2
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















    if( nch == p_par ) then
       call bdy_chem_value( chem_bv_ald,  z, p_ald2, numgas )

       
       call bdy_chem_value( chem_bv_hc3,  z, numgas+1+1, numgas )
       call bdy_chem_value( chem_bv_hc5,  z, numgas+1+2, numgas )
       call bdy_chem_value( chem_bv_hc8,  z, numgas+1+3, numgas )

       call bdy_chem_value( chem_bv_ket,  z, numgas+1+4, numgas )
       call bdy_chem_value( chem_bv_oli,  z, numgas+1+5, numgas )
       call bdy_chem_value( chem_bv_olt,  z, numgas+1+6, numgas )
       call bdy_chem_value( chem_bv_ora2, z, numgas+1+7, numgas )

       chem_bv = 0.4*chem_bv_ald + 2.9*chem_bv_hc3    &
            + 4.8*chem_bv_hc5 + 7.9*chem_bv_hc8       &
            + 0.9*chem_bv_ket + 2.8*chem_bv_oli       &
            + 1.8*chem_bv_olt + 1.0*chem_bv_ora2

    else
       
       chem_bv = chem_bv_def

    end if

	return
      end subroutine bdy_chem_value_cbm4
   


    end module module_cbm4_initmixrats
