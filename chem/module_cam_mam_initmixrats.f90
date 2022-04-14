




      module module_cam_mam_initmixrats

      private
      public :: bdy_chem_value_cam_mam
      public :: cam_mam_initmixrats

      contains



	subroutine cam_mam_initmixrats(         &
		iflagaa, numgas, config_flags,  &
		chem, convfac, alt, z_at_w, g,  &
		ids,ide, jds,jde, kds,kde,      &
		ims,ime, jms,jme, kms,kme,      &
		its,ite, jts,jte, kts,kte       )








	use module_configure, only:  grid_config_rec_type, num_chem, &
		p_ncl_a1, p_ncl_a2, p_nh4_a1, p_nh4_a2, &
		p_so4_a1, p_so4_a2, p_sulf
	use module_state_description, only:  num_chem, param_first_scalar,   &
		aer_ic_default, aer_ic_pnnl
	use module_data_cam_mam_asect, only:  ai_phase, dens_aer, &
		massptr_aer, ncomp_aer, nsize_aer, ntype_aer, numptr_aer, &
		volumcen_sect
	use module_cam_mam_init, only:  pom_icbc_1p4_factor, so4_icbc_1p2_factor

	use module_data_sorgam, only:  conmin

	implicit none



	type(grid_config_rec_type), intent(in) :: config_flags

	integer, intent(in) ::   &
		iflagaa, numgas,   &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem

	real, intent(in),      &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		alt, convfac, z_at_w

        real, intent(in) :: g


        integer :: i, ic, iphase, isize, itype
        integer :: j
        integer :: k
        integer :: l, ll

        real, parameter :: splitfac = 0.98   
        real, parameter :: so4vaptoaer = 0.999   
        real :: tmpa
        real :: zz


        do ic = numgas+1, num_chem
           chem(its:ite,kts:kte,jts:jte,ic) = conmin
        end do

        do j = jts, jte
        do k = kts, kte-1
        do i = its, ite

        if (config_flags%aer_ic_opt == aer_ic_default) then

           chem(i,k,j,p_so4_a1)=chem(i,k,j,p_sulf)*convfac(i,k,j)* &
              96.0*splitfac*so4vaptoaer * so4_icbc_1p2_factor
           chem(i,k,j,p_so4_a2)=chem(i,k,j,p_sulf)*convfac(i,k,j)* &
              96.0*(1.-splitfac)*so4vaptoaer * so4_icbc_1p2_factor
           chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)*(1.-so4vaptoaer)
           chem(i,k,j,p_ncl_a1) = 20.e-05
           chem(i,k,j,p_ncl_a2) = 20.e-05

        else if (config_flags%aer_ic_opt == aer_ic_pnnl) then
           zz = (z_at_w(i,k,j)+z_at_w(i,k+1,j))*0.5
           call cam_mam_init_aer_ic_pnnl(   &
                chem, zz, i,k,j, ims,ime,jms,jme,kms,kme )

        else

           call wrf_error_fatal3("<stdin>",101,&
                "cam_mam_initmixrats: bad value for aer_ic_opt" )
        end if


        iphase = ai_phase
        do itype = 1, ntype_aer
        do isize = 1, nsize_aer(itype)
           tmpa = 0.0
           do ll = 1, ncomp_aer(itype)
              l = massptr_aer(ll,isize,itype,iphase)
              if ((l >= param_first_scalar) .and. (l <= num_chem)) &
                 tmpa = tmpa + chem(i,k,j,l)/dens_aer(ll,itype)
           end do 
           tmpa = tmpa*1.0e-6   
           l = numptr_aer(isize,itype,iphase)
           chem(i,k,j,l) = tmpa/volumcen_sect(isize,itype)
        end do 
        end do 

        end do 
        end do 
        end do 




        do ic = numgas+1, num_chem
        do j = jts,jte
        do k = kts,kte-1
        do i = its,ite
            chem(i,k,j,ic) = chem(i,k,j,ic)*alt(i,k,j)
        end do 
        end do 
        end do 
        end do 


        do ic = numgas+1, num_chem
        do j = jts,jte
        do i = its,ite
            chem(i,kte,j,ic) = chem(i,kte-1,j,ic)
        end do 
        end do 
        end do 


	return
	end subroutine cam_mam_initmixrats











      subroutine cam_mam_init_aer_ic_pnnl(                  &
          chem, z, i,k,j, ims,ime, jms,jme, kms,kme )

      use module_configure, only:  grid_config_rec_type, num_chem, &
          p_bc_a1, p_bc_a3, p_dst_a1, p_dst_a3, p_dst_a5, p_dst_a7, &
          p_ncl_a1, p_ncl_a2, p_ncl_a3, p_ncl_a4, p_ncl_a6, &
          p_nh4_a1, p_nh4_a2, p_pom_a1, p_pom_a3, &
          p_soa_a1, p_soa_a2, p_so4_a1, p_so4_a2, p_sulf
      use module_data_sorgam, only:  conmin
      use module_cam_mam_init, only:  pom_icbc_1p4_factor, so4_icbc_1p2_factor

      implicit none

      integer,intent(in   ) :: i,k,j,                           &
                               ims,ime, jms,jme, kms,kme
      real,  dimension( ims:ime , kms:kme , jms:jme, num_chem ),&
           intent(inout   ) :: chem

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



      chem(i,k,j,p_so4_a1)   = mult*0.300*0.97 * so4_icbc_1p2_factor
      chem(i,k,j,p_so4_a2)   = mult*0.300*0.03 * so4_icbc_1p2_factor











      chem(i,k,j,p_ncl_a2)   = 20.e-05
      chem(i,k,j,p_ncl_a1)   = 20.e-05



      chem(i,k,j,p_bc_a1)    = mult*0.013




      chem(i,k,j,p_dst_a1)   = mult*4.500
      chem(i,k,j,p_dst_a3)   = mult*4.500/2.0



      chem(i,k,j,p_pom_a1)   = mult*0.088 * pom_icbc_1p4_factor

















      chem(i,k,j,p_soa_a1)   = conmin
      chem(i,k,j,p_soa_a2)   = conmin


      chem(i,k,j,p_ncl_a3)   = mult*1.75


      end subroutine cam_mam_init_aer_ic_pnnl










      subroutine bdy_chem_value_cam_mam( chem, z, nch, config_flags, &
                                         alt, convfac, g )

      use module_configure, only:  grid_config_rec_type
      use module_state_description, only:  aer_bc_default, aer_bc_pnnl
      use module_data_sorgam, only: conmin

      implicit none

      integer, intent(in)   :: nch        
      real, intent(out) :: chem
      real, intent(in)  :: z          
      real, intent(in)  :: alt, convfac, g
      type (grid_config_rec_type), intent(in) :: config_flags




      if (config_flags%aer_bc_opt == aer_bc_pnnl) then
         call cam_mam_set_aer_bc_pnnl( chem, z, nch, config_flags )
         return
      else if (config_flags%aer_bc_opt == aer_bc_default) then
         chem = conmin
         return
      else
         call wrf_error_fatal3("<stdin>",296,&
            "bdy_chem_value_cam_mam -- unable to parse aer_bc_opt" )
      end if

      end subroutine bdy_chem_value_cam_mam











      subroutine cam_mam_set_aer_bc_pnnl( chem, z, nch, config_flags )

      use module_configure, only:  grid_config_rec_type, num_chem, &
          p_bc_a1, p_dst_a1, p_dst_a3, p_dst_a5, p_dst_a7, &
          p_ncl_a3, p_ncl_a6, p_nh4_a1, p_nh4_a2, &
          p_num_a3, p_num_a7, p_pom_a1, p_so4_a1, p_so4_a2

      use module_state_description, only:  param_first_scalar

      use module_data_sorgam, only:  conmin

      use module_data_cam_mam_asect, only:  ai_phase, dens_aer, massptr_aer, &
          ncomp_aer, nsize_aer, ntype_aer, numptr_aer, volumcen_sect

      use module_cam_mam_init, only:  pom_icbc_1p4_factor, so4_icbc_1p2_factor

      implicit none

      integer,intent(in) :: nch
      real,intent(in   ) :: z
      real,intent(inout) :: chem
      type (grid_config_rec_type), intent (in) :: config_flags


      integer :: inumber, isize, itype
      integer :: l, ll
      real :: bv_so4ai, bv_so4aj,         &
              bv_nh4ai, bv_nh4aj,         &
              bv_no3ai, bv_no3aj,         &
              bv_eci,   bv_ecj,           &
              bv_p25i,  bv_p25j,          &
              bv_orgpai,bv_orgpaj,        &
              bv_antha, bv_seas, bv_soila
      real :: mult
      real :: tmpchem(num_chem)
      real :: tmpvol

      character(len=160) :: msg


      chem = conmin

      if ((nch < p_so4_a1) .or. (nch > p_num_a3)) return


















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


      bv_so4aj = mult*0.300*0.97 * so4_icbc_1p2_factor
      bv_so4ai = mult*0.300*0.03 * so4_icbc_1p2_factor
      bv_nh4aj = mult*0.094*0.97
      bv_nh4ai = mult*0.094*0.03
      bv_no3aj = mult*0.001*0.97
      bv_no3ai = mult*0.001*0.03
      bv_ecj   = mult*0.013*0.97
      bv_eci   = mult*0.013*0.03
      bv_p25j  = mult*4.500*0.97
      bv_p25i  = mult*4.500*0.03
      bv_antha = mult*4.500/2.0
      bv_orgpaj = mult*0.088*0.97
      bv_orgpai = mult*0.088*0.03
      bv_seas   = mult*1.75
      bv_soila  = conmin




      tmpchem(:) = conmin

      tmpchem(p_so4_a1) = bv_so4aj
      tmpchem(p_so4_a2) = bv_so4ai


      tmpchem(p_pom_a1) = (bv_orgpai + bv_orgpaj) * pom_icbc_1p4_factor

      tmpchem(p_bc_a1 ) = bv_eci + bv_ecj

      tmpchem(p_dst_a1) = bv_p25i + bv_p25j
      tmpchem(p_dst_a3) = bv_antha + bv_soila

      tmpchem(p_ncl_a3) = bv_seas





      inumber = 0
itype_loop01: &
      do itype = 1, ntype_aer
      do isize = 1, nsize_aer(itype)
         if (nch == numptr_aer(isize,itype,ai_phase)) then

            tmpvol = 0.0
            do ll = 1, ncomp_aer(itype)
               l = massptr_aer(ll,isize,itype,ai_phase)
               tmpvol = tmpchem(l)/dens_aer(ll,itype)
            end do




            chem = tmpvol*1.0e-6/volumcen_sect(isize,itype)
            inumber = 1
            exit itype_loop01
         end if
      end do
      end do itype_loop01

      if (inumber <= 0) then

         chem = tmpchem(nch)
      end if


      return
      end subroutine cam_mam_set_aer_bc_pnnl




      end module module_cam_mam_initmixrats

