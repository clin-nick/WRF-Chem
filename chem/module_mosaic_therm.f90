









      module module_mosaic_therm



      use module_data_mosaic_therm
      use module_peg_util



      implicit none

      intrinsic max, min

      contains




































































      subroutine aerchemistry( iclm, jclm, kclm_calcbgn, kclm_calcend,   &
                               dtchem_sngl, idiagaa,vbs_nbin )

      use module_data_mosaic_asect
      use module_data_mosaic_other
      use module_mosaic_movesect, only:  move_sections







      integer iclm, jclm, kclm_calcbgn, kclm_calcend, idiagaa,vbs_nbin(1)
      real dtchem_sngl

      real(kind=8) :: dtchem
      integer k, m



      dtchem = dtchem_sngl

      lunerr_aer = lunerr
      ncorecnt_aer = ncorecnt


      call aerchem_boxtest_output( 1, iclm, jclm, 0, 0, dtchem )

      iclm_aer = iclm
      jclm_aer = jclm
      kclm_aer_calcbgn = kclm_calcbgn
      kclm_aer_calcend = kclm_calcend


      do 200 m = 1, nsubareas
        mclm_aer = m

        do 100 k = kclm_aer_calcbgn, kclm_aer_calcend
          kclm_aer = k
          if (afracsubarea(k,m) .lt. 1.e-4) goto 100

          istat_mosaic_fe1 = 1

          call mosaic( k, m, dtchem,vbs_nbin)

          if (istat_mosaic_fe1 .lt. 0) then
             nfe1_mosaic_cur = nfe1_mosaic_cur + 1
             nfe1_mosaic_tot = nfe1_mosaic_tot + 1
             if (iprint_mosaic_fe1 .gt. 0) then
                write(6,*) 'mosaic aerchemistry fatal error - i/j/k/m =',   &
                   iclm_aer, jclm_aer, kclm_aer, mclm_aer
                call print_input
                if (iprint_mosaic_fe1 .ge. 10)   &
                   call mosaic_aerchem_error_dump( 0, 0, lunerr_aer,   &
                      'aerchemistry fatal error' )
             end if
             goto 100
          end if

          call specialoutaa( iclm, jclm, k, m, 'befor_movesect' )
          call move_sections( 1, iclm, jclm, k, m)
          call specialoutaa( iclm, jclm, k, m, 'after_movesect' )

100     continue	

200   continue		



      call aerchem_boxtest_output( 3, iclm, jclm, 0, 0, dtchem )

      return
      end subroutine aerchemistry
















      subroutine mosaic(k, m, dtchem,vbs_nbin)

      use module_data_mosaic_asect
      use module_data_mosaic_other






      integer k, m,vbs_nbin(1)
      real(kind=8) dtchem

      real(kind=8) yh2o, dumdum
      integer iclm_debug, jclm_debug, kclm_debug, ncnt_debug




      iclm_debug=-28; jclm_debug=1; kclm_debug=9; ncnt_debug=6



      if(iclm_aer .eq. iclm_debug .and.   &
         jclm_aer .eq. jclm_debug .and.   &
         kclm_aer .eq. kclm_debug  .and.   &
         ncorecnt_aer .eq. ncnt_debug)then
        dumdum = 0.0
      endif



         if(1.eq.0)then
           call hijack_input(k,m)
         endif


          t_k = rsub(ktemp,k,m)			
          p_atm = ptotclm(k) /1.032d6		
          yh2o = rsub(kh2o,k,m)			
          rh_pc = 100.*relhumclm(k)		
          ah2o = relhumclm(k)			


          call load_mosaic_parameters		

          call initialize_mosaic_variables

          call update_thermodynamic_constants(vbs_nbin)	

          call map_mosaic_species(k, m, 0)


          call overall_massbal_in 
          iprint_input = myes     


          call mosaic_dynamic_solver( dtchem,vbs_nbin )
          if (istat_mosaic_fe1 .lt. 0) return


          call overall_massbal_out(0) 

          call map_mosaic_species(k, m, 1)



      return
      end subroutine mosaic


















      subroutine mosaic_dynamic_solver( dtchem,vbs_nbin )




      real(kind=8) dtchem

      integer ibin, iv, k, m,vbs_nbin(1)
      real(kind=8) xt, dumdum











      do 500 ibin = 1, nbin_a

        call check_aerosol_mass(ibin)
        if(jaerosolstate(ibin) .eq. no_aerosol)goto 500

        call conform_electrolytes(jtotal,ibin,xt) 	

        call check_aerosol_mass(ibin) 			
        if(jaerosolstate(ibin) .eq. no_aerosol)goto 500	

        call conform_aerosol_number(ibin)   		

500   continue






      call save_pregrow_props

      call specialoutaa( iclm_aer, jclm_aer, kclm_aer, 77,   &
      		'after_conform' )




      if(mgas_aer_xfer .eq. mon)then

        call astem(dtchem,vbs_nbin)

      endif













      do 600 ibin = 1, nbin_a
        if(jaerosolstate(ibin).eq.no_aerosol) goto 600

        if(jhyst_leg(ibin) .eq. jhyst_lo)then
          water_a_hyst(ibin) = 0.0
        elseif(jhyst_leg(ibin) .eq. jhyst_up)then
          water_a_up(ibin)   = aerosol_water_up(ibin)	
          water_a_hyst(ibin) = water_a_up(ibin)
        endif

        call calc_dry_n_wet_aerosol_props(ibin)		
600   continue

      return
      end subroutine mosaic_dynamic_solver














      subroutine hijack_input(k, m)

      use module_data_mosaic_asect
      use module_data_mosaic_other








      integer k, m

      integer ibin, igas, iphase, isize, itype
      real(kind=8) t_kdum, p_atmdum, rhdum, cairclmdum
      real(kind=8) gasdum(4), aerdum(14,8)





      open(92, file = 'box.txt')

      read(92,*)t_kdum, p_atmdum, rhdum, cairclmdum

        read(92,*)gasdum(1),gasdum(2),gasdum(3),gasdum(4)


      do ibin = 1, nbin_a
        read(92,*)aerdum(1,ibin),aerdum(2,ibin),aerdum(3,ibin),   &
                  aerdum(4,ibin),aerdum(5,ibin),aerdum(6,ibin),   &
                  aerdum(7,ibin),aerdum(8,ibin),aerdum(9,ibin),   &
                  aerdum(10,ibin),aerdum(11,ibin),aerdum(12,ibin),   &
                  aerdum(13,ibin),aerdum(14,ibin)
      enddo

      close(92)




      rsub(ktemp,k,m) = t_kdum			
      ptotclm(k)      = p_atmdum*1.032d6
      relhumclm(k)    = rhdum/100.0		
      cairclm(k)      = cairclmdum		




      cair_mol_m3 = cairclm(k)*1.e6	
      cair_mol_cc = cairclm(k)



      conv1a = cair_mol_m3*1.e9		
      conv1b = 1./conv1a		
      conv2a = cair_mol_m3*18.*1.e-3	
      conv2b = 1./conv2a		




        rsub(kh2so4,k,m) = gasdum(1)
        rsub(khno3,k,m)  = gasdum(2)
        rsub(khcl,k,m)   = gasdum(3)
        rsub(knh3,k,m)   = gasdum(4)



        iphase = ai_phase
        ibin = 0
        do 10 itype = 1, ntype_aer
        do 10 isize = 1, nsize_aer(itype)
        ibin = ibin + 1

        rsub(lptr_so4_aer(isize,itype,iphase),k,m) = aerdum(1,ibin)
        rsub(lptr_no3_aer(isize,itype,iphase),k,m) = aerdum(2,ibin)
        rsub(lptr_cl_aer(isize,itype,iphase),k,m)  = aerdum(3,ibin)
        rsub(lptr_nh4_aer(isize,itype,iphase),k,m) = aerdum(4,ibin)
        rsub(lptr_oc_aer(isize,itype,iphase),k,m)  = aerdum(5,ibin)
        rsub(lptr_co3_aer(isize,itype,iphase),k,m) = aerdum(6,ibin)
        rsub(lptr_msa_aer(isize,itype,iphase),k,m) = aerdum(7,ibin)
        rsub(lptr_bc_aer(isize,itype,iphase),k,m)  = aerdum(8,ibin)
        rsub(lptr_na_aer(isize,itype,iphase),k,m)  = aerdum(9,ibin)
        rsub(lptr_ca_aer(isize,itype,iphase),k,m)  = aerdum(10,ibin)
        rsub(lptr_oin_aer(isize,itype,iphase),k,m) = aerdum(11,ibin)

        rsub(hyswptr_aer(isize,itype),k,m) = aerdum(12,ibin) 
        rsub(waterptr_aer(isize,itype),k,m)       = aerdum(13,ibin)	
        rsub(numptr_aer(isize,itype,iphase),k,m)          = aerdum(14,ibin)	
10    continue

      return
      end subroutine hijack_input











      subroutine initialize_mosaic_variables



      integer iaer, ibin, iv, ja, jc, je



      do iv = 1, ngas_ioa
          gas(iv)           = 0.0
      enddo


      do ibin = 1, nbin_a

        num_a(ibin)          = 0.0
        mass_dry_a(ibin)     = 0.0
        mass_soluble_a(ibin) = 0.0

        do iaer = 1, naer
          aer(iaer,jtotal,ibin)  = 0.0
          aer(iaer,jsolid,ibin)  = 0.0
          aer(iaer,jliquid,ibin) = 0.0
        enddo

        do je = 1, nelectrolyte
          electrolyte(je,jtotal,ibin)  = 0.0
          electrolyte(je,jsolid,ibin)  = 0.0
          electrolyte(je,jliquid,ibin) = 0.0
          epercent(je,jtotal,ibin)     = 0.0	
          epercent(je,jsolid,ibin)     = 0.0	
          epercent(je,jliquid,ibin)    = 0.0	
          activity(je,ibin)            = 0.0
          gam(je,ibin)                 = 0.0
        enddo

          gam_ratio(ibin)   = 0.0

        do iv = 1, ngas_ioa
          flux_s(iv,ibin)   = 0.0
          flux_l(iv,ibin)   = 0.0
          kg(iv,ibin)       = 0.0

          phi_volatile_s(iv,ibin) = 0.0
          phi_volatile_l(iv,ibin) = 0.0
          df_gas_s(iv,ibin)   = 0.0
          df_gas_l(iv,ibin)   = 0.0
          volatile_s(iv,ibin) = 0.0
        enddo


        jaerosolstate(ibin) = -1	
        jphase(ibin) = 0

        do jc = 1, ncation
          mc(jc,ibin) = 0.0
        enddo

        do ja = 1, nanion
          ma(ja,ibin) = 0.0
        enddo

      enddo	


      return
      end subroutine initialize_mosaic_variables












      subroutine map_mosaic_species(k, m, imap)

      use module_data_mosaic_asect
      use module_data_mosaic_other
      use module_state_description, only:  param_first_scalar









      integer k, m, imap

      integer ibin, iphase, isize, itsi, itype, l, p1st



      p1st = param_first_scalar



      cair_mol_m3 = cairclm(k)*1.e6	
      cair_mol_cc = cairclm(k)



      conv1a = cair_mol_m3*1.d9		
      conv1b = 1.d0/conv1a		
      conv2a = cair_mol_m3*18.*1.d-3	
      conv2b = 1.d0/conv2a		








      if(imap.eq.0)then    

	if (kh2so4 .ge. p1st) then
	    gas(ih2so4_g) = rsub(kh2so4,k,m)*conv1a	
	else
	    gas(ih2so4_g) = 0.0
	end if
	if (khno3 .ge. p1st) then
	    gas(ihno3_g)  = rsub(khno3,k,m)*conv1a
	else
	    gas(ihno3_g) = 0.0
	end if
	if (khcl .ge. p1st) then
	    gas(ihcl_g)   = rsub(khcl,k,m)*conv1a
	else
	    gas(ihcl_g) = 0.0
	end if
	if (knh3 .ge. p1st) then
	    gas(inh3_g)   = rsub(knh3,k,m)*conv1a
	else
	    gas(inh3_g) = 0.0
	end if
	if (kn2o5 .ge. p1st) then
	    gas(in2o5_g)   = rsub(kn2o5,k,m)*conv1a
	else
	    gas(in2o5_g) = 0.0
	end if
	if (kclno2 .ge. p1st) then
	    gas(iclno2_g)   = rsub(kclno2,k,m)*conv1a
	else
	    gas(iclno2_g) = 0.0
	end if


        if (kpcg1_b_c .ge. p1st) then
            gas(ipcg1_b_c_g)   = rsub(kpcg1_b_c,k,m)*conv1a
        else
        gas(ipcg1_b_c_g) = 0.0
        end if
        if (kpcg2_b_c .ge. p1st) then
            gas(ipcg2_b_c_g)   = rsub(kpcg2_b_c,k,m)*conv1a
        else
        gas(ipcg2_b_c_g) = 0.0
        end if
        if (kpcg3_b_c .ge. p1st) then
            gas(ipcg3_b_c_g)   = rsub(kpcg3_b_c,k,m)*conv1a
        else
        gas(ipcg3_b_c_g) = 0.0
        end if
        if (kpcg4_b_c .ge. p1st) then
            gas(ipcg4_b_c_g)   = rsub(kpcg4_b_c,k,m)*conv1a
        else
        gas(ipcg4_b_c_g) = 0.0
        end if
        if (kpcg5_b_c .ge. p1st) then
            gas(ipcg5_b_c_g)   = rsub(kpcg5_b_c,k,m)*conv1a
        else
        gas(ipcg5_b_c_g) = 0.0
        end if
        if (kpcg6_b_c .ge. p1st) then
            gas(ipcg6_b_c_g)   = rsub(kpcg6_b_c,k,m)*conv1a
        else
        gas(ipcg6_b_c_g) = 0.0
        end if
        if (kpcg7_b_c .ge. p1st) then
            gas(ipcg7_b_c_g)   = rsub(kpcg7_b_c,k,m)*conv1a
        else
        gas(ipcg7_b_c_g) = 0.0
        end if
        if (kpcg8_b_c .ge. p1st) then
            gas(ipcg8_b_c_g)   = rsub(kpcg8_b_c,k,m)*conv1a
        else
        gas(ipcg8_b_c_g) = 0.0
        end if
        if (kpcg9_b_c .ge. p1st) then
            gas(ipcg9_b_c_g)   = rsub(kpcg9_b_c,k,m)*conv1a
        else
        gas(ipcg9_b_c_g) = 0.0
        end if
        if (kpcg1_b_o .ge. p1st) then
            gas(ipcg1_b_o_g)   = rsub(kpcg1_b_o,k,m)*conv1a
        else
        gas(ipcg1_b_o_g) = 0.0
        end if
        if (kpcg2_b_o .ge. p1st) then
            gas(ipcg2_b_o_g)   = rsub(kpcg2_b_o,k,m)*conv1a
        else
        gas(ipcg2_b_o_g) = 0.0
        end if
        if (kpcg3_b_o .ge. p1st) then
            gas(ipcg3_b_o_g)   = rsub(kpcg3_b_o,k,m)*conv1a
        else
        gas(ipcg3_b_o_g) = 0.0
        end if
        if (kpcg4_b_o .ge. p1st) then
            gas(ipcg4_b_o_g)   = rsub(kpcg4_b_o,k,m)*conv1a
        else
        gas(ipcg4_b_o_g) = 0.0
        end if
        if (kpcg5_b_o .ge. p1st) then
            gas(ipcg5_b_o_g)   = rsub(kpcg5_b_o,k,m)*conv1a
        else
        gas(ipcg5_b_o_g) = 0.0
        end if
        if (kpcg6_b_o .ge. p1st) then
            gas(ipcg6_b_o_g)   = rsub(kpcg6_b_o,k,m)*conv1a
        else
        gas(ipcg6_b_o_g) = 0.0
        end if
        if (kpcg7_b_o .ge. p1st) then
            gas(ipcg7_b_o_g)   = rsub(kpcg7_b_o,k,m)*conv1a
        else
        gas(ipcg7_b_o_g) = 0.0
        end if
        if (kpcg8_b_o .ge. p1st) then
            gas(ipcg8_b_o_g)   = rsub(kpcg8_b_o,k,m)*conv1a
        else
        gas(ipcg8_b_o_g) = 0.0
        end if
        if (kpcg9_b_o .ge. p1st) then
            gas(ipcg9_b_o_g)   = rsub(kpcg9_b_o,k,m)*conv1a
        else
        gas(ipcg9_b_o_g) = 0.0
        end if
        if (kopcg1_b_c .ge. p1st) then
            gas(iopcg1_b_c_g)   = rsub(kopcg1_b_c,k,m)*conv1a
        else
        gas(iopcg1_b_c_g) = 0.0
        end if
        if (kopcg2_b_c .ge. p1st) then
            gas(iopcg2_b_c_g)   = rsub(kopcg2_b_c,k,m)*conv1a
        else
        gas(iopcg2_b_c_g) = 0.0
        end if
        if (kopcg3_b_c .ge. p1st) then
            gas(iopcg3_b_c_g)   = rsub(kopcg3_b_c,k,m)*conv1a
        else
        gas(iopcg3_b_c_g) = 0.0
        end if
        if (kopcg4_b_c .ge. p1st) then
            gas(iopcg4_b_c_g)   = rsub(kopcg4_b_c,k,m)*conv1a
        else
        gas(iopcg4_b_c_g) = 0.0
        end if
        if (kopcg5_b_c .ge. p1st) then
            gas(iopcg5_b_c_g)   = rsub(kopcg5_b_c,k,m)*conv1a
        else
        gas(iopcg5_b_c_g) = 0.0
        end if
        if (kopcg6_b_c .ge. p1st) then
            gas(iopcg6_b_c_g)   = rsub(kopcg6_b_c,k,m)*conv1a
        else
        gas(iopcg6_b_c_g) = 0.0
        end if
        if (kopcg7_b_c .ge. p1st) then
            gas(iopcg7_b_c_g)   = rsub(kopcg7_b_c,k,m)*conv1a
        else
        gas(iopcg7_b_c_g) = 0.0
        end if
        if (kopcg8_b_c .ge. p1st) then
            gas(iopcg8_b_c_g)   = rsub(kopcg8_b_c,k,m)*conv1a
        else
        gas(iopcg8_b_c_g) = 0.0
        end if
        if (kopcg1_b_o .ge. p1st) then
            gas(iopcg1_b_o_g)   = rsub(kopcg1_b_o,k,m)*conv1a
        else
        gas(iopcg1_b_o_g) = 0.0
        end if
        if (kopcg2_b_o .ge. p1st) then
            gas(iopcg2_b_o_g)   = rsub(kopcg2_b_o,k,m)*conv1a
        else
        gas(iopcg2_b_o_g) = 0.0
        end if
        if (kopcg3_b_o .ge. p1st) then
            gas(iopcg3_b_o_g)   = rsub(kopcg3_b_o,k,m)*conv1a
        else
        gas(iopcg3_b_o_g) = 0.0
        end if
        if (kopcg4_b_o .ge. p1st) then
            gas(iopcg4_b_o_g)   = rsub(kopcg4_b_o,k,m)*conv1a
        else
        gas(iopcg4_b_o_g) = 0.0
        end if
        if (kopcg5_b_o .ge. p1st) then
            gas(iopcg5_b_o_g)   = rsub(kopcg5_b_o,k,m)*conv1a
        else
        gas(iopcg5_b_o_g) = 0.0
        end if
        if (kopcg6_b_o .ge. p1st) then
            gas(iopcg6_b_o_g)   = rsub(kopcg6_b_o,k,m)*conv1a
        else
        gas(iopcg6_b_o_g) = 0.0
        end if
        if (kopcg7_b_o .ge. p1st) then
            gas(iopcg7_b_o_g)   = rsub(kopcg7_b_o,k,m)*conv1a
        else
        gas(iopcg7_b_o_g) = 0.0
        end if
        if (kopcg8_b_o .ge. p1st) then
            gas(iopcg8_b_o_g)   = rsub(kopcg8_b_o,k,m)*conv1a
        else
        gas(iopcg8_b_o_g) = 0.0
        end if
        if (kpcg1_f_c .ge. p1st) then
            gas(ipcg1_f_c_g)   = rsub(kpcg1_f_c,k,m)*conv1a
        else
        gas(ipcg1_f_c_g) = 0.0
        end if
        if (kpcg2_f_c .ge. p1st) then
            gas(ipcg2_f_c_g)   = rsub(kpcg2_f_c,k,m)*conv1a
        else
        gas(ipcg2_f_c_g) = 0.0
        end if
        if (kpcg3_f_c .ge. p1st) then
            gas(ipcg3_f_c_g)   = rsub(kpcg3_f_c,k,m)*conv1a
        else
        gas(ipcg3_f_c_g) = 0.0
        end if
        if (kpcg4_f_c .ge. p1st) then
            gas(ipcg4_f_c_g)   = rsub(kpcg4_f_c,k,m)*conv1a
        else
        gas(ipcg4_f_c_g) = 0.0
        end if
        if (kpcg5_f_c .ge. p1st) then
            gas(ipcg5_f_c_g)   = rsub(kpcg5_f_c,k,m)*conv1a
        else
        gas(ipcg5_f_c_g) = 0.0
        end if
        if (kpcg6_f_c .ge. p1st) then
            gas(ipcg6_f_c_g)   = rsub(kpcg6_f_c,k,m)*conv1a
        else
        gas(ipcg6_f_c_g) = 0.0
        end if
        if (kpcg7_f_c .ge. p1st) then
            gas(ipcg7_f_c_g)   = rsub(kpcg7_f_c,k,m)*conv1a
        else
        gas(ipcg7_f_c_g) = 0.0
        end if
        if (kpcg8_f_c .ge. p1st) then
            gas(ipcg8_f_c_g)   = rsub(kpcg8_f_c,k,m)*conv1a
        else
        gas(ipcg8_f_c_g) = 0.0
        end if
        if (kpcg9_f_c .ge. p1st) then
            gas(ipcg9_f_c_g)   = rsub(kpcg9_f_c,k,m)*conv1a
        else
        gas(ipcg9_f_c_g) = 0.0
        end if
        if (kpcg1_f_o .ge. p1st) then
            gas(ipcg1_f_o_g)   = rsub(kpcg1_f_o,k,m)*conv1a
        else
        gas(ipcg1_f_o_g) = 0.0
        end if
        if (kpcg2_f_o .ge. p1st) then
            gas(ipcg2_f_o_g)   = rsub(kpcg2_f_o,k,m)*conv1a
        else
        gas(ipcg2_f_o_g) = 0.0
        end if
        if (kpcg3_f_o .ge. p1st) then
            gas(ipcg3_f_o_g)   = rsub(kpcg3_f_o,k,m)*conv1a
        else
        gas(ipcg3_f_o_g) = 0.0
        end if
        if (kpcg4_f_o .ge. p1st) then
            gas(ipcg4_f_o_g)   = rsub(kpcg4_f_o,k,m)*conv1a
        else
        gas(ipcg4_f_o_g) = 0.0
        end if
        if (kpcg5_f_o .ge. p1st) then
            gas(ipcg5_f_o_g)   = rsub(kpcg5_f_o,k,m)*conv1a
        else
        gas(ipcg5_f_o_g) = 0.0
        end if
        if (kpcg6_f_o .ge. p1st) then
            gas(ipcg6_f_o_g)   = rsub(kpcg6_f_o,k,m)*conv1a
        else
        gas(ipcg6_f_o_g) = 0.0
        end if
        if (kpcg7_f_o .ge. p1st) then
            gas(ipcg7_f_o_g)   = rsub(kpcg7_f_o,k,m)*conv1a
        else
        gas(ipcg7_f_o_g) = 0.0
        end if
        if (kpcg8_f_o .ge. p1st) then
            gas(ipcg8_f_o_g)   = rsub(kpcg8_f_o,k,m)*conv1a
        else
        gas(ipcg8_f_o_g) = 0.0
        end if
        if (kpcg9_f_o .ge. p1st) then
            gas(ipcg9_f_o_g)   = rsub(kpcg9_f_o,k,m)*conv1a
        else
        gas(ipcg9_f_o_g) = 0.0
        end if
        if (kopcg1_f_c .ge. p1st) then
            gas(iopcg1_f_c_g)   = rsub(kopcg1_f_c,k,m)*conv1a
        else
        gas(iopcg1_f_c_g) = 0.0
        end if
        if (kopcg2_f_c .ge. p1st) then
            gas(iopcg2_f_c_g)   = rsub(kopcg2_f_c,k,m)*conv1a
        else
        gas(iopcg2_f_c_g) = 0.0
        end if
        if (kopcg3_f_c .ge. p1st) then
            gas(iopcg3_f_c_g)   = rsub(kopcg3_f_c,k,m)*conv1a
        else
        gas(iopcg3_f_c_g) = 0.0
        end if
        if (kopcg4_f_c .ge. p1st) then
            gas(iopcg4_f_c_g)   = rsub(kopcg4_f_c,k,m)*conv1a
        else
        gas(iopcg4_f_c_g) = 0.0
        end if
        if (kopcg5_f_c .ge. p1st) then
            gas(iopcg5_f_c_g)   = rsub(kopcg5_f_c,k,m)*conv1a
        else
        gas(iopcg5_f_c_g) = 0.0
        end if
        if (kopcg6_f_c .ge. p1st) then
            gas(iopcg6_f_c_g)   = rsub(kopcg6_f_c,k,m)*conv1a
        else
        gas(iopcg6_f_c_g) = 0.0
        end if
        if (kopcg7_f_c .ge. p1st) then
            gas(iopcg7_f_c_g)   = rsub(kopcg7_f_c,k,m)*conv1a
        else
        gas(iopcg7_f_c_g) = 0.0
        end if
        if (kopcg8_f_c .ge. p1st) then
            gas(iopcg8_f_c_g)   = rsub(kopcg8_f_c,k,m)*conv1a
        else
        gas(iopcg8_f_c_g) = 0.0
        end if
        if (kopcg1_f_o .ge. p1st) then
            gas(iopcg1_f_o_g)   = rsub(kopcg1_f_o,k,m)*conv1a
        else
        gas(iopcg1_f_o_g) = 0.0
        end if
        if (kopcg2_f_o .ge. p1st) then
            gas(iopcg2_f_o_g)   = rsub(kopcg2_f_o,k,m)*conv1a
        else
        gas(iopcg2_f_o_g) = 0.0
        end if
        if (kopcg3_f_o .ge. p1st) then
            gas(iopcg3_f_o_g)   = rsub(kopcg3_f_o,k,m)*conv1a
        else
        gas(iopcg3_f_o_g) = 0.0
        end if
        if (kopcg4_f_o .ge. p1st) then
            gas(iopcg4_f_o_g)   = rsub(kopcg4_f_o,k,m)*conv1a
        else
        gas(iopcg4_f_o_g) = 0.0
        end if
        if (kopcg5_f_o .ge. p1st) then
            gas(iopcg5_f_o_g)   = rsub(kopcg5_f_o,k,m)*conv1a
        else
        gas(iopcg5_f_o_g) = 0.0
        end if
        if (kopcg6_f_o .ge. p1st) then
            gas(iopcg6_f_o_g)   = rsub(kopcg6_f_o,k,m)*conv1a
        else
        gas(iopcg6_f_o_g) = 0.0
        end if
        if (kopcg7_f_o .ge. p1st) then
            gas(iopcg7_f_o_g)   = rsub(kopcg7_f_o,k,m)*conv1a
        else
        gas(iopcg7_f_o_g) = 0.0
        end if
        if (kopcg8_f_o .ge. p1st) then
            gas(iopcg8_f_o_g)   = rsub(kopcg8_f_o,k,m)*conv1a
        else
        gas(iopcg8_f_o_g) = 0.0
        end if

       if (ksmpa .ge. p1st) then
            gas(ismpa_g)   = rsub(ksmpa,k,m)*conv1a
        else
        gas(ismpa_g) = 0.0
        end if
        if (ksmpbb .ge. p1st) then
            gas(ismpbb_g)   = rsub(ksmpbb,k,m)*conv1a
        else
        gas(ismpbb_g) = 0.0
        end if
       if (kgly .ge. p1st) then
            gas(igly)   = rsub(kgly,k,m)*conv1a
        else
        gas(igly) = 0.0
        end if
        if (koh .ge. p1st) then
            gas(iho)   = rsub(koh,k,m)*conv1a
        else
        gas(koh) = 0.0
        end if


        if (kant1_c .ge. p1st) then
            gas(iant1_c_g)   = rsub(kant1_c,k,m)*conv1a
        else
        gas(iant1_c_g) = 0.0
        end if
        if (kant2_c .ge. p1st) then
            gas(iant2_c_g)   = rsub(kant2_c,k,m)*conv1a
        else
        gas(iant2_c_g) = 0.0
        end if
        if (kant3_c .ge. p1st) then
            gas(iant3_c_g)   = rsub(kant3_c,k,m)*conv1a
        else
        gas(iant3_c_g) = 0.0
        end if
        if (kant4_c .ge. p1st) then
            gas(iant4_c_g)   = rsub(kant4_c,k,m)*conv1a
        else
        gas(iant4_c_g) = 0.0
        end if

        if (kant1_o .ge. p1st) then
            gas(iant1_o_g)   = rsub(kant1_o,k,m)*conv1a
        else
        gas(iant1_o_g) = 0.0
        end if
        if (kant2_o .ge. p1st) then
            gas(iant2_o_g)   = rsub(kant2_o,k,m)*conv1a
        else
        gas(iant2_o_g) = 0.0
        end if
        if (kant3_o .ge. p1st) then
            gas(iant3_o_g)   = rsub(kant3_o,k,m)*conv1a
        else
        gas(iant3_o_g) = 0.0
        end if
        if (kant4_o .ge. p1st) then
            gas(iant4_o_g)   = rsub(kant4_o,k,m)*conv1a
        else
        gas(iant4_o_g) = 0.0
        end if

        if (kbiog1_c .ge. p1st) then
            gas(ibiog1_c_g)   = rsub(kbiog1_c,k,m)*conv1a
        else
        gas(ibiog1_c_g) = 0.0
        end if
        if (kbiog2_c .ge. p1st) then
            gas(ibiog2_c_g)   = rsub(kbiog2_c,k,m)*conv1a
        else
        gas(ibiog2_c_g) = 0.0
        end if
        if (kbiog3_c .ge. p1st) then
            gas(ibiog3_c_g)   = rsub(kbiog3_c,k,m)*conv1a
        else
        gas(ibiog3_c_g) = 0.0
        end if
        if (kbiog4_c .ge. p1st) then
            gas(ibiog4_c_g)   = rsub(kbiog4_c,k,m)*conv1a
        else
        gas(ibiog4_c_g) = 0.0
        end if

        if (kbiog1_o .ge. p1st) then
            gas(ibiog1_o_g)   = rsub(kbiog1_o,k,m)*conv1a
        else
        gas(ibiog1_o_g) = 0.0
        end if
        if (kbiog2_o .ge. p1st) then
            gas(ibiog2_o_g)   = rsub(kbiog2_o,k,m)*conv1a
        else
        gas(ibiog2_o_g) = 0.0
        end if
        if (kbiog3_o .ge. p1st) then
            gas(ibiog3_o_g)   = rsub(kbiog3_o,k,m)*conv1a
        else
        gas(ibiog3_o_g) = 0.0
        end if
        if (kbiog4_o .ge. p1st) then
            gas(ibiog4_o_g)   = rsub(kbiog4_o,k,m)*conv1a
        else
        gas(ibiog4_o_g) = 0.0
        end if

        if (kasoaX .ge. p1st) then
            gas(iasoaX_g)   = rsub(kasoaX,k,m)*conv1a
        else
        gas(iasoaX_g) = 0.0
        end if

        if (kasoa1 .ge. p1st) then
            gas(iasoa1_g)   = rsub(kasoa1,k,m)*conv1a
        else
        gas(iasoa1_g) = 0.0
        end if

        if (kasoa2 .ge. p1st) then
            gas(iasoa2_g)   = rsub(kasoa2,k,m)*conv1a
        else
        gas(iasoa2_g) = 0.0
        end if

        if (kasoa3 .ge. p1st) then
            gas(iasoa3_g)   = rsub(kasoa3,k,m)*conv1a
        else
        gas(iasoa3_g) = 0.0
        end if

        if (kasoa4 .ge. p1st) then
            gas(iasoa4_g)   = rsub(kasoa4,k,m)*conv1a
        else
        gas(iasoa4_g) = 0.0
        end if

        if (kbsoaX .ge. p1st) then
            gas(ibsoaX_g)   = rsub(kbsoaX,k,m)*conv1a
        else
        gas(ibsoaX_g) = 0.0
        end if

        if (kbsoa1 .ge. p1st) then
            gas(ibsoa1_g)   = rsub(kbsoa1,k,m)*conv1a
        else
        gas(ibsoa1_g) = 0.0
        end if

        if (kbsoa2 .ge. p1st) then
            gas(ibsoa2_g)   = rsub(kbsoa2,k,m)*conv1a
        else
        gas(ibsoa2_g) = 0.0
        end if

        if (kbsoa3 .ge. p1st) then
            gas(ibsoa3_g)   = rsub(kbsoa3,k,m)*conv1a
        else
        gas(ibsoa3_g) = 0.0
        end if

        if (kbsoa4 .ge. p1st) then
            gas(ibsoa4_g)   = rsub(kbsoa4,k,m)*conv1a
        else
        gas(ibsoa4_g) = 0.0
        end if






        iphase = ai_phase
        ibin = 0
        do 10 itype = 1, ntype_aer
        do 10 isize = 1, nsize_aer(itype)
        ibin = ibin + 1






        l = lptr_so4_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(iso4_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(iso4_a,jtotal,ibin)=0.0
        end if

        l = lptr_no3_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(ino3_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(ino3_a,jtotal,ibin)=0.0
        end if

        l = lptr_cl_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(icl_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(icl_a,jtotal,ibin)=0.0
        end if

        l = lptr_nh4_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(inh4_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(inh4_a,jtotal,ibin)=0.0
        end if

        l = lptr_oc_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(ioc_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(ioc_a,jtotal,ibin)=0.0
        end if

        l = lptr_bc_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(ibc_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(ibc_a,jtotal,ibin)=0.0
        end if

        l = lptr_na_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(ina_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(ina_a,jtotal,ibin)=0.0
        end if

        l = lptr_oin_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(ioin_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(ioin_a,jtotal,ibin)=0.0
        end if

        l = lptr_msa_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(imsa_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(imsa_a,jtotal,ibin)=0.0
        end if

        l = lptr_co3_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(ico3_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(ico3_a,jtotal,ibin)=0.0
        end if

        l = lptr_ca_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            aer(ica_a,jtotal,ibin)=rsub(l,k,m)*conv1a
        else
            aer(ica_a,jtotal,ibin)=0.0
        end if



       l = lptr_pcg1_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg1_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg1_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg2_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg2_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg2_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg3_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg3_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg3_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg4_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg4_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg4_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg5_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg5_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg5_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg6_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg6_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg6_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg7_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg7_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg7_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg8_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg8_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg8_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg9_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg9_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg9_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg1_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg1_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg1_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg2_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg2_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg2_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg3_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg3_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg3_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg4_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg4_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg4_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg5_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg5_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg5_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg6_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg6_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg6_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg7_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg7_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg7_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg8_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg8_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg8_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg9_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg9_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg9_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg1_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg1_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg1_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg2_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg2_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg2_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg3_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg3_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg3_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg4_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg4_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg4_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg5_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg5_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg5_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg6_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg6_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg6_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg7_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg7_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg7_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg8_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg8_b_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg8_b_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg1_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg1_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg1_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg2_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg2_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg2_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg3_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg3_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg3_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg4_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg4_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg4_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg5_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg5_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg5_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg6_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg6_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg6_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg7_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg7_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg7_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg8_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg8_b_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg8_b_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg1_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg1_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg1_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg2_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg2_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg2_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg3_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg3_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg3_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg4_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg4_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg4_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg5_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg5_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg5_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg6_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg6_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg6_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg7_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg7_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg7_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg8_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg8_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg8_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg9_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg9_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg9_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg1_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg1_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg1_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg2_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg2_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg2_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg3_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg3_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg3_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg4_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg4_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg4_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg5_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg5_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg5_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg6_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg6_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg6_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg7_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg7_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg7_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg8_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg8_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg8_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_pcg9_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ipcg9_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ipcg9_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg1_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg1_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg1_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg2_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg2_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg2_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg3_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg3_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg3_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg4_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg4_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg4_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg5_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg5_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg5_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg6_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg6_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg6_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg7_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg7_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg7_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg8_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg8_f_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg8_f_c_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg1_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg1_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg1_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg2_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg2_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg2_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg3_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg3_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg3_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg4_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg4_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg4_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg5_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg5_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg5_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg6_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg6_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg6_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg7_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg7_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg7_f_o_a,jtotal,ibin)=0.0
       end if
       l = lptr_opcg8_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iopcg8_f_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iopcg8_f_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_smpa_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ismpa_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ismpa_a,jtotal,ibin)=0.0
       end if
       l = lptr_smpbb_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ismpbb_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ismpbb_a,jtotal,ibin)=0.0
       end if

       l = lptr_glysoa_r1_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iglysoa_r1_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iglysoa_r1_a,jtotal,ibin)=0.0
       end if

       l = lptr_glysoa_r2_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iglysoa_r2_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iglysoa_r2_a,jtotal,ibin)=0.0
       end if

       l = lptr_glysoa_sfc_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iglysoa_sfc_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iglysoa_sfc_a,jtotal,ibin)=0.0
       end if

       l = lptr_glysoa_nh4_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iglysoa_nh4_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iglysoa_nh4_a,jtotal,ibin)=0.0
       end if

       l = lptr_glysoa_oh_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iglysoa_oh_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iglysoa_oh_a,jtotal,ibin)=0.0
       end if

       l = lptr_ant1_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iant1_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iant1_c_a,jtotal,ibin)=0.0
       end if

       l = lptr_ant2_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iant2_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iant2_c_a,jtotal,ibin)=0.0
       end if

       l = lptr_ant3_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iant3_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iant3_c_a,jtotal,ibin)=0.0
       end if

       l = lptr_ant4_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iant4_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iant4_c_a,jtotal,ibin)=0.0
       end if

       l = lptr_ant1_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iant1_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iant1_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_ant2_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iant2_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iant2_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_ant3_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iant3_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iant3_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_ant4_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iant4_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iant4_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_biog1_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibiog1_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibiog1_c_a,jtotal,ibin)=0.0
       end if

       l = lptr_biog2_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibiog2_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibiog2_c_a,jtotal,ibin)=0.0
       end if

       l = lptr_biog3_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibiog3_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibiog3_c_a,jtotal,ibin)=0.0
       end if

       l = lptr_biog4_c_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibiog4_c_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibiog4_c_a,jtotal,ibin)=0.0
       end if

       l = lptr_biog1_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibiog1_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibiog1_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_biog2_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibiog2_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibiog2_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_biog3_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibiog3_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibiog3_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_biog4_o_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibiog4_o_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibiog4_o_a,jtotal,ibin)=0.0
       end if

       l = lptr_asoaX_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iasoaX_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iasoaX_a,jtotal,ibin)=0.0
       end if

       l = lptr_asoa1_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iasoa1_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iasoa1_a,jtotal,ibin)=0.0
       end if

       l = lptr_asoa2_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iasoa2_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iasoa2_a,jtotal,ibin)=0.0
       end if

       l = lptr_asoa3_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iasoa3_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iasoa3_a,jtotal,ibin)=0.0
       end if

       l = lptr_asoa4_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(iasoa4_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(iasoa4_a,jtotal,ibin)=0.0
       end if

       l = lptr_bsoaX_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibsoaX_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibsoaX_a,jtotal,ibin)=0.0
       end if

       l = lptr_bsoa1_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibsoa1_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibsoa1_a,jtotal,ibin)=0.0
       end if

       l = lptr_bsoa2_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibsoa2_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibsoa2_a,jtotal,ibin)=0.0
       end if

       l = lptr_bsoa3_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibsoa3_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibsoa3_a,jtotal,ibin)=0.0
       end if

       l = lptr_bsoa4_aer(isize,itype,iphase)
       if (l .ge. p1st) then
           aer(ibsoa4_a,jtotal,ibin)=rsub(l,k,m)*conv1a
       else
            aer(ibsoa4_a,jtotal,ibin)=0.0
       end if


        l = hyswptr_aer(isize,itype)
        if (l .ge. p1st) then
            water_a_hyst(ibin)=rsub(l,k,m)*conv2a
        else
            water_a_hyst(ibin)=0.0
        end if


        l = waterptr_aer(isize,itype)
        if (l .ge. p1st) then
            water_a(ibin)=rsub(l,k,m)*conv2a
        else
            water_a(ibin)=0.0
        end if


        l = numptr_aer(isize,itype,iphase)
        num_a(ibin) = rsub(l,k,m)*cair_mol_cc


        sigmag_a(ibin)	= 1.02

10      continue







      else                 




	if (kh2so4 .ge. p1st)   &
	    rsub(kh2so4,k,m) = gas(ih2so4_g)*conv1b
	if (khno3 .ge. p1st)   &
	    rsub(khno3,k,m)  = gas(ihno3_g)*conv1b
	if (khcl .ge. p1st)   &
	    rsub(khcl,k,m)   = gas(ihcl_g)*conv1b
	if (knh3 .ge. p1st)   &
	    rsub(knh3,k,m)   = gas(inh3_g)*conv1b
	if (kn2o5 .ge. p1st)   &
	    rsub(kn2o5,k,m)   = gas(in2o5_g)*conv1b
	if (kclno2 .ge. p1st)   &
	    rsub(kclno2,k,m)   = gas(iclno2_g)*conv1b


        if (kpcg1_b_c .ge. p1st)   &
            rsub(kpcg1_b_c,k,m)   = gas(ipcg1_b_c_g)*conv1b
        if (kpcg2_b_c .ge. p1st)   &
            rsub(kpcg2_b_c,k,m)   = gas(ipcg2_b_c_g)*conv1b
        if (kpcg3_b_c .ge. p1st)   &
            rsub(kpcg3_b_c,k,m)   = gas(ipcg3_b_c_g)*conv1b
        if (kpcg4_b_c .ge. p1st)   &
            rsub(kpcg4_b_c,k,m)   = gas(ipcg4_b_c_g)*conv1b
        if (kpcg5_b_c .ge. p1st)   &
            rsub(kpcg5_b_c,k,m)   = gas(ipcg5_b_c_g)*conv1b
        if (kpcg6_b_c .ge. p1st)   &
            rsub(kpcg6_b_c,k,m)   = gas(ipcg6_b_c_g)*conv1b
        if (kpcg7_b_c .ge. p1st)   &
            rsub(kpcg7_b_c,k,m)   = gas(ipcg7_b_c_g)*conv1b
        if (kpcg8_b_c .ge. p1st)   &
            rsub(kpcg8_b_c,k,m)   = gas(ipcg8_b_c_g)*conv1b
        if (kpcg9_b_c .ge. p1st)   &
            rsub(kpcg9_b_c,k,m)   = gas(ipcg9_b_c_g)*conv1b
        if (kpcg1_b_o .ge. p1st)   &
            rsub(kpcg1_b_o,k,m)   = gas(ipcg1_b_o_g)*conv1b
        if (kpcg2_b_o .ge. p1st)   &
            rsub(kpcg2_b_o,k,m)   = gas(ipcg2_b_o_g)*conv1b
        if (kpcg3_b_o .ge. p1st)   &
            rsub(kpcg3_b_o,k,m)   = gas(ipcg3_b_o_g)*conv1b
        if (kpcg4_b_o .ge. p1st)   &
            rsub(kpcg4_b_o,k,m)   = gas(ipcg4_b_o_g)*conv1b
        if (kpcg5_b_o .ge. p1st)   &
            rsub(kpcg5_b_o,k,m)   = gas(ipcg5_b_o_g)*conv1b
        if (kpcg6_b_o .ge. p1st)   &
            rsub(kpcg6_b_o,k,m)   = gas(ipcg6_b_o_g)*conv1b
        if (kpcg7_b_o .ge. p1st)   &
            rsub(kpcg7_b_o,k,m)   = gas(ipcg7_b_o_g)*conv1b
        if (kpcg8_b_o .ge. p1st)   &
            rsub(kpcg8_b_o,k,m)   = gas(ipcg8_b_o_g)*conv1b
        if (kpcg9_b_o .ge. p1st)   &
            rsub(kpcg9_b_o,k,m)   = gas(ipcg9_b_o_g)*conv1b
        if (kopcg1_b_c .ge. p1st)   &
            rsub(kopcg1_b_c,k,m)   = gas(iopcg1_b_c_g)*conv1b
        if (kopcg2_b_c .ge. p1st)   &
            rsub(kopcg2_b_c,k,m)   = gas(iopcg2_b_c_g)*conv1b
        if (kopcg3_b_c .ge. p1st)   &
            rsub(kopcg3_b_c,k,m)   = gas(iopcg3_b_c_g)*conv1b
        if (kopcg4_b_c .ge. p1st)   &
            rsub(kopcg4_b_c,k,m)   = gas(iopcg4_b_c_g)*conv1b
        if (kopcg5_b_c .ge. p1st)   &
            rsub(kopcg5_b_c,k,m)   = gas(iopcg5_b_c_g)*conv1b
        if (kopcg6_b_c .ge. p1st)   &
            rsub(kopcg6_b_c,k,m)   = gas(iopcg6_b_c_g)*conv1b
        if (kopcg7_b_c .ge. p1st)   &
            rsub(kopcg7_b_c,k,m)   = gas(iopcg7_b_c_g)*conv1b
        if (kopcg8_b_c .ge. p1st)   &
            rsub(kopcg8_b_c,k,m)   = gas(iopcg8_b_c_g)*conv1b
        if (kopcg1_b_o .ge. p1st)   &
            rsub(kopcg1_b_o,k,m)   = gas(iopcg1_b_o_g)*conv1b
        if (kopcg2_b_o .ge. p1st)   &
            rsub(kopcg2_b_o,k,m)   = gas(iopcg2_b_o_g)*conv1b
        if (kopcg3_b_o .ge. p1st)   &
            rsub(kopcg3_b_o,k,m)   = gas(iopcg3_b_o_g)*conv1b
        if (kopcg4_b_o .ge. p1st)   &
            rsub(kopcg4_b_o,k,m)   = gas(iopcg4_b_o_g)*conv1b
        if (kopcg5_b_o .ge. p1st)   &
            rsub(kopcg5_b_o,k,m)   = gas(iopcg5_b_o_g)*conv1b
        if (kopcg6_b_o .ge. p1st)   &
            rsub(kopcg6_b_o,k,m)   = gas(iopcg6_b_o_g)*conv1b
        if (kopcg7_b_o .ge. p1st)   &
            rsub(kopcg7_b_o,k,m)   = gas(iopcg7_b_o_g)*conv1b
        if (kopcg8_b_o .ge. p1st)   &
            rsub(kopcg8_b_o,k,m)   = gas(iopcg8_b_o_g)*conv1b
        if (kpcg1_f_c .ge. p1st)   &
            rsub(kpcg1_f_c,k,m)   = gas(ipcg1_f_c_g)*conv1b
        if (kpcg2_f_c .ge. p1st)   &
            rsub(kpcg2_f_c,k,m)   = gas(ipcg2_f_c_g)*conv1b
        if (kpcg3_f_c .ge. p1st)   &
            rsub(kpcg3_f_c,k,m)   = gas(ipcg3_f_c_g)*conv1b
        if (kpcg4_f_c .ge. p1st)   &
            rsub(kpcg4_f_c,k,m)   = gas(ipcg4_f_c_g)*conv1b
        if (kpcg5_f_c .ge. p1st)   &
            rsub(kpcg5_f_c,k,m)   = gas(ipcg5_f_c_g)*conv1b
        if (kpcg6_f_c .ge. p1st)   &
            rsub(kpcg6_f_c,k,m)   = gas(ipcg6_f_c_g)*conv1b
        if (kpcg7_f_c .ge. p1st)   &
            rsub(kpcg7_f_c,k,m)   = gas(ipcg7_f_c_g)*conv1b
        if (kpcg8_f_c .ge. p1st)   &
            rsub(kpcg8_f_c,k,m)   = gas(ipcg8_f_c_g)*conv1b
        if (kpcg9_f_c .ge. p1st)   &
            rsub(kpcg9_f_c,k,m)   = gas(ipcg9_f_c_g)*conv1b
        if (kpcg1_f_o .ge. p1st)   &
            rsub(kpcg1_f_o,k,m)   = gas(ipcg1_f_o_g)*conv1b
        if (kpcg2_f_o .ge. p1st)   &
            rsub(kpcg2_f_o,k,m)   = gas(ipcg2_f_o_g)*conv1b
        if (kpcg3_f_o .ge. p1st)   &
            rsub(kpcg3_f_o,k,m)   = gas(ipcg3_f_o_g)*conv1b
        if (kpcg4_f_o .ge. p1st)   &
            rsub(kpcg4_f_o,k,m)   = gas(ipcg4_f_o_g)*conv1b
        if (kpcg5_f_o .ge. p1st)   &
            rsub(kpcg5_f_o,k,m)   = gas(ipcg5_f_o_g)*conv1b
        if (kpcg6_f_o .ge. p1st)   &
            rsub(kpcg6_f_o,k,m)   = gas(ipcg6_f_o_g)*conv1b
        if (kpcg7_f_o .ge. p1st)   &
            rsub(kpcg7_f_o,k,m)   = gas(ipcg7_f_o_g)*conv1b
        if (kpcg8_f_o .ge. p1st)   &
            rsub(kpcg8_f_o,k,m)   = gas(ipcg8_f_o_g)*conv1b
        if (kpcg9_f_o .ge. p1st)   &
            rsub(kpcg9_f_o,k,m)   = gas(ipcg9_f_o_g)*conv1b
        if (kopcg1_f_c .ge. p1st)   &
            rsub(kopcg1_f_c,k,m)   = gas(iopcg1_f_c_g)*conv1b
        if (kopcg2_f_c .ge. p1st)   &
            rsub(kopcg2_f_c,k,m)   = gas(iopcg2_f_c_g)*conv1b
        if (kopcg3_f_c .ge. p1st)   &
            rsub(kopcg3_f_c,k,m)   = gas(iopcg3_f_c_g)*conv1b
        if (kopcg4_f_c .ge. p1st)   &
            rsub(kopcg4_f_c,k,m)   = gas(iopcg4_f_c_g)*conv1b
        if (kopcg5_f_c .ge. p1st)   &
            rsub(kopcg5_f_c,k,m)   = gas(iopcg5_f_c_g)*conv1b
        if (kopcg6_f_c .ge. p1st)   &
            rsub(kopcg6_f_c,k,m)   = gas(iopcg6_f_c_g)*conv1b
        if (kopcg7_f_c .ge. p1st)   &
            rsub(kopcg7_f_c,k,m)   = gas(iopcg7_f_c_g)*conv1b
        if (kopcg8_f_c .ge. p1st)   &
            rsub(kopcg8_f_c,k,m)   = gas(iopcg8_f_c_g)*conv1b
        if (kopcg1_f_o .ge. p1st)   &
            rsub(kopcg1_f_o,k,m)   = gas(iopcg1_f_o_g)*conv1b
        if (kopcg2_f_o .ge. p1st)   &
            rsub(kopcg2_f_o,k,m)   = gas(iopcg2_f_o_g)*conv1b
        if (kopcg3_f_o .ge. p1st)   &
            rsub(kopcg3_f_o,k,m)   = gas(iopcg3_f_o_g)*conv1b
        if (kopcg4_f_o .ge. p1st)   &
            rsub(kopcg4_f_o,k,m)   = gas(iopcg4_f_o_g)*conv1b
        if (kopcg5_f_o .ge. p1st)   &
            rsub(kopcg5_f_o,k,m)   = gas(iopcg5_f_o_g)*conv1b
        if (kopcg6_f_o .ge. p1st)   &
            rsub(kopcg6_f_o,k,m)   = gas(iopcg6_f_o_g)*conv1b
        if (kopcg7_f_o .ge. p1st)   &
            rsub(kopcg7_f_o,k,m)   = gas(iopcg7_f_o_g)*conv1b
        if (kopcg8_f_o .ge. p1st)   &
            rsub(kopcg8_f_o,k,m)   = gas(iopcg8_f_o_g)*conv1b
        if (ksmpa .ge. p1st)   &
            rsub(ksmpa,k,m)   = gas(ismpa_g)*conv1b
        if (kgly .ge. p1st)   &
            rsub(kgly,k,m)   = gas(igly)*conv1b
        


        if (ksmpbb .ge. p1st)   &
            rsub(ksmpbb,k,m)   = gas(ismpbb_g)*conv1b
        if (kant1_c .ge. p1st)   &
            rsub(kant1_c,k,m)   = gas(iant1_c_g)*conv1b
        if (kant2_c .ge. p1st)   &
            rsub(kant2_c,k,m)   = gas(iant2_c_g)*conv1b
        if (kant3_c .ge. p1st)   &
            rsub(kant3_c,k,m)   = gas(iant3_c_g)*conv1b
        if (kant4_c .ge. p1st)   &
            rsub(kant4_c,k,m)   = gas(iant4_c_g)*conv1b
        if (kant1_o .ge. p1st)   &
            rsub(kant1_o,k,m)   = gas(iant1_o_g)*conv1b
        if (kant2_o .ge. p1st)   &
            rsub(kant2_o,k,m)   = gas(iant2_o_g)*conv1b
        if (kant3_o .ge. p1st)   &
            rsub(kant3_o,k,m)   = gas(iant3_o_g)*conv1b
        if (kant4_o .ge. p1st)   &
            rsub(kant4_o,k,m)   = gas(iant4_o_g)*conv1b
        if (kbiog1_c .ge. p1st)   &
            rsub(kbiog1_c,k,m)   = gas(ibiog1_c_g)*conv1b
        if (kbiog2_c .ge. p1st)   &
            rsub(kbiog2_c,k,m)   = gas(ibiog2_c_g)*conv1b
        if (kbiog3_c .ge. p1st)   &
            rsub(kbiog3_c,k,m)   = gas(ibiog3_c_g)*conv1b
        if (kbiog4_c .ge. p1st)   &
            rsub(kbiog4_c,k,m)   = gas(ibiog4_c_g)*conv1b
        if (kbiog1_o .ge. p1st)   &
            rsub(kbiog1_o,k,m)   = gas(ibiog1_o_g)*conv1b
        if (kbiog2_o .ge. p1st)   &
            rsub(kbiog2_o,k,m)   = gas(ibiog2_o_g)*conv1b
        if (kbiog3_o .ge. p1st)   &
            rsub(kbiog3_o,k,m)   = gas(ibiog3_o_g)*conv1b
        if (kbiog4_o .ge. p1st)   &
            rsub(kbiog4_o,k,m)   = gas(ibiog4_o_g)*conv1b
        if (kasoaX .ge. p1st)   &
            rsub(kasoaX,k,m)   = gas(iasoaX_g)*conv1b
        if (kasoa1 .ge. p1st)   &
            rsub(kasoa1,k,m)   = gas(iasoa1_g)*conv1b
        if (kasoa2 .ge. p1st)   &
            rsub(kasoa2,k,m)   = gas(iasoa2_g)*conv1b
        if (kasoa3 .ge. p1st)   &
            rsub(kasoa3,k,m)   = gas(iasoa3_g)*conv1b
        if (kasoa4 .ge. p1st)   &
            rsub(kasoa4,k,m)   = gas(iasoa4_g)*conv1b
        if (kbsoaX .ge. p1st)   &
            rsub(kbsoaX,k,m)   = gas(ibsoaX_g)*conv1b
        if (kbsoa1 .ge. p1st)   &
            rsub(kbsoa1,k,m)   = gas(ibsoa1_g)*conv1b
        if (kbsoa2 .ge. p1st)   &
            rsub(kbsoa2,k,m)   = gas(ibsoa2_g)*conv1b
        if (kbsoa3 .ge. p1st)   &
            rsub(kbsoa3,k,m)   = gas(ibsoa3_g)*conv1b
        if (kbsoa4 .ge. p1st)   &
            rsub(kbsoa4,k,m)   = gas(ibsoa4_g)*conv1b


        iphase = ai_phase
        ibin = 0
        do 20 itype = 1, ntype_aer
        do 20 isize = 1, nsize_aer(itype)
        ibin = ibin + 1




        l = lptr_so4_aer(isize,itype,iphase)
        rsub(l,k,m) = aer(iso4_a,jtotal,ibin)*conv1b

        l = lptr_no3_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(ino3_a,jtotal,ibin)*conv1b

        l = lptr_cl_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(icl_a,jtotal,ibin)*conv1b

        l = lptr_nh4_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(inh4_a,jtotal,ibin)*conv1b

        l = lptr_oc_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(ioc_a,jtotal,ibin)*conv1b

        l = lptr_bc_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(ibc_a,jtotal,ibin)*conv1b

        l = lptr_na_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(ina_a,jtotal,ibin)*conv1b

        l = lptr_oin_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(ioin_a,jtotal,ibin)*conv1b

        l = lptr_msa_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(imsa_a,jtotal,ibin)*conv1b

        l = lptr_co3_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(ico3_a,jtotal,ibin)*conv1b

        l = lptr_ca_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) = aer(ica_a,jtotal,ibin)*conv1b



       l = lptr_pcg1_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg1_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg2_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg2_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg3_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg3_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg4_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg4_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg5_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg5_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg6_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg6_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg7_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg7_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg8_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg8_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg9_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg9_b_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg1_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg1_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg2_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg2_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg3_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg3_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg4_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg4_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg5_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg5_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg6_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg6_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg7_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg7_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg8_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg8_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg9_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg9_b_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg1_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg1_b_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg2_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg2_b_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg3_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg3_b_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg4_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg4_b_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg5_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg5_b_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg6_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg6_b_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg7_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg7_b_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg8_b_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg8_b_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg1_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg1_b_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg2_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg2_b_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg3_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg3_b_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg4_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg4_b_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg5_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg5_b_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg6_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg6_b_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg7_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg7_b_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg8_b_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg8_b_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg1_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg1_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg2_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg2_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg3_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg3_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg4_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg4_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg5_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg5_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg6_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg6_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg7_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg7_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg8_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg8_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg9_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg9_f_c_a,jtotal,ibin)*conv1b
       l = lptr_pcg1_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg1_f_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg2_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg2_f_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg3_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg3_f_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg4_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg4_f_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg5_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg5_f_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg6_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg6_f_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg7_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg7_f_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg8_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg8_f_o_a,jtotal,ibin)*conv1b
       l = lptr_pcg9_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ipcg9_f_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg1_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg1_f_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg2_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg2_f_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg3_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg3_f_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg4_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg4_f_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg5_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg5_f_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg6_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg6_f_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg7_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg7_f_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg8_f_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg8_f_c_a,jtotal,ibin)*conv1b
       l = lptr_opcg1_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg1_f_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg2_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg2_f_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg3_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg3_f_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg4_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg4_f_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg5_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg5_f_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg6_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg6_f_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg7_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg7_f_o_a,jtotal,ibin)*conv1b
       l = lptr_opcg8_f_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iopcg8_f_o_a,jtotal,ibin)*conv1b

       l = lptr_smpa_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ismpa_a,jtotal,ibin)*conv1b
       l = lptr_smpbb_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ismpbb_a,jtotal,ibin)*conv1b
       l = lptr_glysoa_r1_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iglysoa_r1_a,jtotal,ibin)*conv1b
       l = lptr_glysoa_r2_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iglysoa_r2_a,jtotal,ibin)*conv1b
       l = lptr_glysoa_sfc_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iglysoa_sfc_a,jtotal,ibin)*conv1b
       l = lptr_glysoa_nh4_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iglysoa_nh4_a,jtotal,ibin)*conv1b
       l = lptr_glysoa_oh_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iglysoa_oh_a,jtotal,ibin)*conv1b

       l = lptr_ant1_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iant1_c_a,jtotal,ibin)*conv1b
       l = lptr_ant2_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iant2_c_a,jtotal,ibin)*conv1b
       l = lptr_ant3_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iant3_c_a,jtotal,ibin)*conv1b
       l = lptr_ant4_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iant4_c_a,jtotal,ibin)*conv1b
       l = lptr_ant1_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iant1_o_a,jtotal,ibin)*conv1b
       l = lptr_ant2_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iant2_o_a,jtotal,ibin)*conv1b
       l = lptr_ant3_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iant3_o_a,jtotal,ibin)*conv1b
       l = lptr_ant4_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iant4_o_a,jtotal,ibin)*conv1b
       l = lptr_biog1_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibiog1_c_a,jtotal,ibin)*conv1b
       l = lptr_biog2_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibiog2_c_a,jtotal,ibin)*conv1b
       l = lptr_biog3_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibiog3_c_a,jtotal,ibin)*conv1b
       l = lptr_biog4_c_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibiog4_c_a,jtotal,ibin)*conv1b
       l = lptr_biog1_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibiog1_o_a,jtotal,ibin)*conv1b
       l = lptr_biog2_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibiog2_o_a,jtotal,ibin)*conv1b
       l = lptr_biog3_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibiog3_o_a,jtotal,ibin)*conv1b
       l = lptr_biog4_o_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibiog4_o_a,jtotal,ibin)*conv1b

       l = lptr_asoaX_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iasoaX_a,jtotal,ibin)*conv1b
       l = lptr_asoa1_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iasoa1_a,jtotal,ibin)*conv1b
       l = lptr_asoa2_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iasoa2_a,jtotal,ibin)*conv1b
       l = lptr_asoa3_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iasoa3_a,jtotal,ibin)*conv1b
       l = lptr_asoa4_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(iasoa4_a,jtotal,ibin)*conv1b
       l = lptr_bsoaX_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibsoaX_a,jtotal,ibin)*conv1b
       l = lptr_bsoa1_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibsoa1_a,jtotal,ibin)*conv1b
       l = lptr_bsoa2_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibsoa2_a,jtotal,ibin)*conv1b
       l = lptr_bsoa3_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibsoa3_a,jtotal,ibin)*conv1b
       l = lptr_bsoa4_aer(isize,itype,iphase)
       if (l .ge. p1st) rsub(l,k,m) = aer(ibsoa4_a,jtotal,ibin)*conv1b




        l = hyswptr_aer(isize,itype)
        if (l .ge. p1st) rsub(l,k,m) = water_a_hyst(ibin)*conv2b

        l = waterptr_aer(isize,itype)
        if (l .ge. p1st) rsub(l,k,m) = water_a(ibin)*conv2b

        l = numptr_aer(isize,itype,iphase)
        if (l .ge. p1st) rsub(l,k,m) =  num_a(ibin)/cair_mol_cc


        drymass_aftgrow(isize,itype) = mass_dry_a(ibin)/cair_mol_cc 
        if(jaerosolstate(ibin) .eq. no_aerosol) then
	    drydens_aftgrow(isize,itype) = -1.
	else
            drydens_aftgrow(isize,itype) = dens_dry_a(ibin)         
	end if

20      continue

      endif

      return
      end subroutine map_mosaic_species





      subroutine isize_itype_from_ibin( ibin, isize, itype )





      use module_data_mosaic_asect
      use module_data_mosaic_other, only:  lunerr



      integer ibin, isize, itype

      integer jdum_bin, jdum_size, jdum_type
      character*80 msg

      isize = -999888777
      itype = -999888777

      jdum_bin = 0
      do jdum_type = 1, ntype_aer
      do jdum_size = 1, nsize_aer(jdum_type)
          jdum_bin = jdum_bin + 1
          if (ibin .eq. jdum_bin) then
              isize = jdum_size
              itype = jdum_type
          end if
      end do
      end do

      if (isize .le. 0) then
          write(msg,'(a,1x,i5)')   &
              '*** subr isize_itype_from_ibin - bad ibin =', ibin
          call peg_error_fatal( lunerr, msg )
      end if

      return
      end subroutine isize_itype_from_ibin




      subroutine overall_massbal_in

      use module_data_mosaic_asect
      use module_data_mosaic_other



      integer ibin

      tot_so4_in = gas(ih2so4_g)
      tot_no3_in = gas(ihno3_g)
      tot_cl_in  = gas(ihcl_g)
      tot_nh4_in = gas(inh3_g)
      tot_na_in  = 0.0
      tot_ca_in  = 0.0


      do ibin = 1, nbin_a
        tot_so4_in = tot_so4_in + aer(iso4_a,jtotal,ibin)
	tot_no3_in = tot_no3_in + aer(ino3_a,jtotal,ibin)
        tot_cl_in  = tot_cl_in  + aer(icl_a, jtotal,ibin)
        tot_nh4_in = tot_nh4_in + aer(inh4_a,jtotal,ibin)
        tot_na_in  = tot_na_in  + aer(ina_a,jtotal,ibin)
        tot_ca_in  = tot_ca_in  + aer(ica_a,jtotal,ibin)
      enddo


        total_species(inh3_g) = tot_nh4_in
        total_species(ihno3_g)= tot_no3_in
        total_species(ihcl_g) = tot_cl_in


      return
      end subroutine overall_massbal_in



      subroutine overall_massbal_out(mbin)








      integer mbin

      integer ibin



        tot_so4_out = gas(ih2so4_g)
	tot_no3_out = gas(ihno3_g)
        tot_cl_out  = gas(ihcl_g)
        tot_nh4_out = gas(inh3_g)
        tot_na_out  = 0.0
        tot_ca_out  = 0.0

	do ibin = 1, nbin_a
          tot_so4_out = tot_so4_out + aer(iso4_a,jtotal,ibin)
	  tot_no3_out = tot_no3_out + aer(ino3_a,jtotal,ibin)
          tot_cl_out  = tot_cl_out  + aer(icl_a,jtotal,ibin)
          tot_nh4_out = tot_nh4_out + aer(inh4_a,jtotal,ibin)
          tot_na_out  = tot_na_out  + aer(ina_a,jtotal,ibin)
          tot_ca_out  = tot_ca_out  + aer(ica_a,jtotal,ibin)
	enddo

        diff_so4 = tot_so4_out - tot_so4_in
	diff_no3 = tot_no3_out - tot_no3_in
        diff_cl  = tot_cl_out  - tot_cl_in
        diff_nh4 = tot_nh4_out - tot_nh4_in
        diff_na  = tot_na_out  - tot_na_in
        diff_ca  = tot_ca_out  - tot_ca_in


        reldiff_so4 = 0.0
	if(tot_so4_in .gt. 1.e-25 .or. tot_so4_out .gt. 1.e-25)then
	  reldiff_so4 = diff_so4/max(tot_so4_in, tot_so4_out)
	endif

        reldiff_no3 = 0.0
	if(tot_no3_in .gt. 1.e-25 .or. tot_no3_out .gt. 1.e-25)then
	  reldiff_no3 = diff_no3/max(tot_no3_in, tot_no3_out)
	endif

        reldiff_cl = 0.0
	if(tot_cl_in .gt. 1.e-25 .or. tot_cl_out .gt. 1.e-25)then
	  reldiff_cl = diff_cl/max(tot_cl_in, tot_cl_out)
	endif

        reldiff_nh4 = 0.0
	if(tot_nh4_in .gt. 1.e-25 .or. tot_nh4_out .gt. 1.e-25)then
	  reldiff_nh4 = diff_nh4/max(tot_nh4_in, tot_nh4_out)
	endif

        reldiff_na = 0.0
	if(tot_na_in .gt. 1.e-25 .or. tot_na_out .gt. 1.e-25)then
	  reldiff_na = diff_na/max(tot_na_in, tot_na_out)
	endif

        reldiff_ca = 0.0
	if(tot_ca_in .gt. 1.e-25 .or. tot_ca_out .gt. 1.e-25)then
	  reldiff_ca = diff_ca/max(tot_ca_in, tot_ca_out)
	endif



      if(  abs(reldiff_so4) .gt. 1.e-4 .or.   &
           abs(reldiff_no3) .gt. 1.e-4 .or.   &
           abs(reldiff_cl)  .gt. 1.e-4 .or.   &
           abs(reldiff_nh4) .gt. 1.e-4 .or.   &
           abs(reldiff_na)  .gt. 1.e-4 .or.   &
           abs(reldiff_ca)  .gt. 1.e-4)then


        if (iprint_mosaic_diag1 .gt. 0) then
          if (iprint_input .eq. myes) then
            write(6,*)'*** mbin = ', mbin, '  isteps = ', isteps_ASTEM
            write(6,*)'reldiff_so4 = ', reldiff_so4
            write(6,*)'reldiff_no3 = ', reldiff_no3
            write(6,*)'reldiff_cl  = ', reldiff_cl
            write(6,*)'reldiff_nh4 = ', reldiff_nh4
            write(6,*)'reldiff_na  = ', reldiff_na
            write(6,*)'reldiff_ca  = ', reldiff_ca
            call print_input
            iprint_input = mno
          endif
        endif

      endif


      return
      end subroutine overall_massbal_out







      subroutine print_input

      use module_data_mosaic_asect
      use module_data_mosaic_other








      integer k, m

      integer ibin, iphase, isize, itype
      integer ipasstmp, luntmp



        if (iprint_mosaic_input_ok .le. 0) return
        if (iprint_input .ne. myes) return
        iprint_input = mno

        k = kclm_aer
        m = mclm_aer


        tot_so4_out = gas(ih2so4_g)
        tot_no3_out = gas(ihno3_g)
        tot_cl_out  = gas(ihcl_g)
        tot_nh4_out = gas(inh3_g)
        tot_na_out  = 0.0
        tot_ca_out  = 0.0

	do ibin = 1, nbin_a
          tot_so4_out = tot_so4_out + aer(iso4_a,jtotal,ibin)
          tot_no3_out = tot_no3_out + aer(ino3_a,jtotal,ibin)
          tot_cl_out  = tot_cl_out  + aer(icl_a,jtotal,ibin)
          tot_nh4_out = tot_nh4_out + aer(inh4_a,jtotal,ibin)
          tot_na_out  = tot_na_out  + aer(ina_a,jtotal,ibin)
          tot_ca_out  = tot_ca_out  + aer(ica_a,jtotal,ibin)
	enddo

        diff_so4 = tot_so4_out - tot_so4_in
	diff_no3 = tot_no3_out - tot_no3_in
        diff_cl  = tot_cl_out  - tot_cl_in
        diff_nh4 = tot_nh4_out - tot_nh4_in
        diff_na  = tot_na_out  - tot_na_in
        diff_ca  = tot_ca_out  - tot_ca_in


        reldiff_so4 = 0.0
	if(tot_so4_in .gt. 1.e-25 .or. tot_so4_out .gt. 1.e-25)then
	  reldiff_so4 = diff_so4/max(tot_so4_in, tot_so4_out)
	endif

        reldiff_no3 = 0.0
	if(tot_no3_in .gt. 1.e-25 .or. tot_no3_out .gt. 1.e-25)then
	  reldiff_no3 = diff_no3/max(tot_no3_in, tot_no3_out)
	endif

        reldiff_cl = 0.0
	if(tot_cl_in .gt. 1.e-25 .or. tot_cl_out .gt. 1.e-25)then
	  reldiff_cl = diff_cl/max(tot_cl_in, tot_cl_out)
	endif

        reldiff_nh4 = 0.0
	if(tot_nh4_in .gt. 1.e-25 .or. tot_nh4_out .gt. 1.e-25)then
	  reldiff_nh4 = diff_nh4/max(tot_nh4_in, tot_nh4_out)
	endif

        reldiff_na = 0.0
	if(tot_na_in .gt. 1.e-25 .or. tot_na_out .gt. 1.e-25)then
	  reldiff_na = diff_na/max(tot_na_in, tot_na_out)
	endif

        reldiff_ca = 0.0
	if(tot_ca_in .gt. 1.e-25 .or. tot_ca_out .gt. 1.e-25)then
	  reldiff_ca = diff_ca/max(tot_ca_in, tot_ca_out)
	endif


        do 2900 ipasstmp = 1, 2

        if (ipasstmp .eq. 1) then
           luntmp = 6     
        else
           luntmp = 67    

        endif


          write(luntmp,*)'+++++++++++++++++++++++++++++++++++++++++'
          write(luntmp,*)'i j k n = ', iclm_aer, jclm_aer, kclm_aer,   &
                                  ncorecnt_aer
          write(luntmp,*)'relative so4 mass bal = ', reldiff_so4
	  write(luntmp,*)'relative no3 mass bal = ', reldiff_no3
          write(luntmp,*)'relative cl  mass bal = ', reldiff_cl
          write(luntmp,*)'relative nh4 mass bal = ', reldiff_nh4
          write(luntmp,*)'relative na  mass bal = ', reldiff_na
          write(luntmp,*)'relative ca  mass bal = ', reldiff_ca
          write(luntmp,*)'inputs:'
          write(luntmp,*)'t (k), p (atm), rh (%), cair (mol/cc) = '
          write(luntmp,44) t_k, p_atm, rh_pc, cairclm(k)
	  write(luntmp,*)'gas h2so4, hno3, hcl, nh3 (mol/mol)'
	  write(luntmp,44)rsub(kh2so4,k,m), rsub(khno3,k,m),   &
                          rsub(khcl,k,m), rsub(knh3,k,m)


	  iphase = ai_phase
          ibin = 0
          do itype = 1, ntype_aer
          do isize = 1, nsize_aer(itype)
          ibin = ibin + 1

	  write(luntmp,44) rsub(lptr_so4_aer(ibin,itype,iphase),k,m),   &
                      rsub(lptr_no3_aer(ibin,itype,iphase),k,m),   &
                      rsub(lptr_cl_aer(ibin,itype,iphase),k,m),   &
                      rsub(lptr_nh4_aer(ibin,itype,iphase),k,m),   &
                      rsub(lptr_oc_aer(ibin,itype,iphase),k,m),	   &  
                      rsub(lptr_co3_aer(ibin,itype,iphase),k,m),   &
                      rsub(lptr_msa_aer(ibin,itype,iphase),k,m),   &
                      rsub(lptr_bc_aer(ibin,itype,iphase),k,m),	   &  
                      rsub(lptr_na_aer(ibin,itype,iphase),k,m),   &
                      rsub(lptr_ca_aer(ibin,itype,iphase),k,m),   &
                      rsub(lptr_oin_aer(ibin,itype,iphase),k,m),	   &
                      rsub(hyswptr_aer(ibin,itype),k,m),   &
                      rsub(waterptr_aer(ibin,itype),k,m),   &
                      rsub(numptr_aer(ibin,itype,iphase),k,m)
          enddo
          enddo

          write(luntmp,*)'+++++++++++++++++++++++++++++++++++++++++'

2900    continue


44      format(14e20.10)



      return
      end subroutine print_input

























      subroutine check_aerosol_mass(ibin)



      integer ibin

      integer iaer
      real(kind=8) drymass, aer_H

      mass_dry_a(ibin) = 0.0

      aer_H = (2.*aer(iso4_a,jtotal,ibin) +  &
                  aer(ino3_a,jtotal,ibin) +  &
                  aer(icl_a,jtotal,ibin)  +  &
                  aer(imsa_a,jtotal,ibin) +  &
               2.*aer(ico3_a,jtotal,ibin))-  &
              (2.*aer(ica_a,jtotal,ibin)  +  &
                  aer(ina_a,jtotal,ibin)  +  &
                  aer(inh4_a,jtotal,ibin))


      do iaer = 1, naer
        mass_dry_a(ibin) = mass_dry_a(ibin) +   &
                           aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)	
      enddo
      mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H

      drymass = mass_dry_a(ibin)			
      mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15	

      if(drymass .lt. mass_cutoff)then			
        jaerosolstate(ibin) = no_aerosol
        jphase(ibin) = 0
        if(drymass .eq. 0.)num_a(ibin) = 0.0
      endif

      return
      end subroutine check_aerosol_mass

















      subroutine conform_aerosol_number(ibin)

      use module_data_mosaic_asect







      integer ibin

      integer je, l, iaer, isize, itype
      real(kind=8) num_at_dlo, num_at_dhi, numold
      real(kind=8) aer_H

      vol_dry_a(ibin)  = 0.0		

      if(jaerosolstate(ibin) .eq. no_aerosol) return

      aer_H = (2.*aer(iso4_a,jtotal,ibin) +  &
                  aer(ino3_a,jtotal,ibin) +  &
                  aer(icl_a,jtotal,ibin)  +  &
                  aer(imsa_a,jtotal,ibin) +  &
               2.*aer(ico3_a,jtotal,ibin))-  &
              (2.*aer(ica_a,jtotal,ibin)  +  &
                  aer(ina_a,jtotal,ibin)  +  &
                  aer(inh4_a,jtotal,ibin))

      do iaer = 1, naer
        vol_dry_a(ibin) = vol_dry_a(ibin) +   &
        aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  
      enddo
      vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

      vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15	


      call isize_itype_from_ibin( ibin, isize, itype )
      num_at_dlo = vol_dry_a(ibin)/volumlo_sect(isize,itype)
      num_at_dhi = vol_dry_a(ibin)/volumhi_sect(isize,itype)

      numold = num_a(ibin)
      num_a(ibin) = min(num_a(ibin), num_at_dlo) 
      num_a(ibin) = max(num_a(ibin), num_at_dhi) 













      return
      end subroutine conform_aerosol_number











      subroutine aerosol_phase_state(ibin)



      integer ibin

      integer js, je, iaer, iv, iter_kelvin
      real(kind=8) ah2o_a_new, rel_err

      real(kind=8) kelvin_toler, term
      real(kind=8) aer_H


      ah2o = rh_pc*0.01
      ah2o_a(ibin) = ah2o
      kelvin(ibin) = 1.0
      do iv = 1, ngas_volatile+ngas_het
        kel(iv,ibin) = 1.0
      enddo

      if(rh_pc .le. 99)then
        kelvin_toler = 1.e-2
      else
        kelvin_toler = 1.e-6
      endif


      mass_dry_a(ibin) = 0.0		
      vol_dry_a(ibin)  = 0.0		

      aer_H = (2.*aer(iso4_a,jtotal,ibin) +  &
                  aer(ino3_a,jtotal,ibin) +  &
                  aer(icl_a,jtotal,ibin)  +  &
                  aer(imsa_a,jtotal,ibin) +  &
               2.*aer(ico3_a,jtotal,ibin))-  &
              (2.*aer(ica_a,jtotal,ibin)  +  &
                  aer(ina_a,jtotal,ibin)  +  &
                  aer(inh4_a,jtotal,ibin))

      do iaer = 1, naer
        mass_dry_a(ibin) = mass_dry_a(ibin) +   &
                           aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)	
        vol_dry_a(ibin)  = vol_dry_a(ibin) +   &
        aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  	
      enddo
      mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
      vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

      mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15			
      vol_dry_a(ibin)  = vol_dry_a(ibin)*1.e-15				


      mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3		
      vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3		


      water_a_up(ibin) = aerosol_water_up(ibin)	

      iter_kelvin = 0

10    iter_kelvin = iter_kelvin + 1
      do je = 1, nelectrolyte
        molality0(je) = bin_molality(je,ibin)	
      enddo

      call mesa(ibin)
      if(jaerosolstate(ibin) .eq. all_solid)then
        return
      endif
      if (istat_mosaic_fe1 .lt. 0) return


      mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3		
      vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3		

      call calculate_kelvin(ibin)

      ah2o_a_new = rh_pc*0.01/kelvin(ibin)

      rel_err = abs( (ah2o_a_new - ah2o_a(ibin))/ah2o_a(ibin))

      if(rel_err .gt. kelvin_toler .and. iter_kelvin.le.10)then
        ah2o_a(ibin) = ah2o_a_new
        goto 10
      endif

      if(jaerosolstate(ibin) .eq. all_liquid)jhyst_leg(ibin) = jhyst_up


      do iv = 1,  ngas_volatile+ngas_het
        term = 4.*sigma_soln(ibin)*partial_molar_vol(iv)/  &
                       (8.3144e7*T_K*DpmV(ibin))
        kel(iv,ibin) = 1. + term*(1. + 0.5*term*(1. + term/3.))
      enddo


      return
      end subroutine aerosol_phase_state












      subroutine calculate_kelvin(ibin)



      integer ibin

      real(kind=8) term



      volume_a(ibin) = vol_wet_a(ibin) 					
      dpmv(ibin)=(6.*volume_a(ibin)/(num_a(ibin)*3.1415926))**(1./3.)	
      sigma_soln(ibin) = sigma_water + 49.0*(1. - ah2o_a(ibin)) 	
      term = 72.*sigma_soln(ibin)/(8.3144e7*t_k*dpmv(ibin))		

      kelvin(ibin) = 1. + term*(1. + 0.5*term*(1. + term/3.))


      return
      end subroutine calculate_kelvin

























      subroutine mesa(ibin)	



      integer ibin


      integer idissolved, j_index, jdum, js, je		
      real(kind=8) crh, solids, sum_soluble, sum_insoluble, xt


      real(kind=8) h_ion, sum_dum				



      sum_dum = 0.0
      do je = 1, nelectrolyte
        sum_dum = sum_dum + electrolyte(je,jtotal,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      do je = 1, nelectrolyte
        epercent(je,jtotal,ibin) = 100.*electrolyte(je,jtotal,ibin)/sum_dum
      enddo


      call calculate_xt(ibin,jtotal,xt)

      crh = 0.35  


      if( (ah2o_a(ibin) .lt. crh)   .and. &
          (xt.gt.1.0 .or. xt.lt.0.) .and. &
          (epercent(jcano3,jtotal,ibin) .le. ptol_mol_astem) .and. &
          (epercent(jcacl2,jtotal,ibin) .le. ptol_mol_astem) )then     
        jaerosolstate(ibin) = all_solid
        jphase(ibin)    = jsolid
        jhyst_leg(ibin) = jhyst_lo
        call adjust_solid_aerosol(ibin)
        return
      endif







      if (mhyst_method .eq. mhyst_uporlo_waterhyst) then
         jdum = 0
         if (water_a_hyst(ibin) .gt. 0.5*water_a_up(ibin)) jdum = 1
      else if (mhyst_method .eq. mhyst_force_up) then
         jdum = 1
      else 
         jdum = 0
      end if

      if (jdum .eq. 1) then

        call do_full_deliquescence(ibin)


























          jaerosolstate(ibin) = all_liquid
          jhyst_leg(ibin) = jhyst_up
          jphase(ibin) = jliquid
          water_a(ibin) = aerosol_water(jtotal,ibin)

          if(water_a(ibin) .lt. 0.0)then    
            jaerosolstate(ibin) = all_solid 
            jphase(ibin)    = jsolid
            jhyst_leg(ibin) = jhyst_lo
            call adjust_solid_aerosol(ibin)
          else
            call adjust_liquid_aerosol(ibin)
            call compute_activities(ibin)
          endif


          mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3	
          vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3	
          growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)	

          return



      endif





      if(xt .lt. 1. .and. xt .gt. 0. )goto 10	

      jdum = 0
      do js = 1, nsalt
        jsalt_present(js) = 0			

        if(epercent(js,jtotal,ibin) .gt. ptol_mol_astem)then
          jsalt_present(js) = 1			
          jdum = jdum + jsalt_index(js)
        endif
      enddo

      if(jdum .eq. 0)then
        jaerosolstate(ibin) = all_solid 
        jphase(ibin) = jsolid
        call adjust_solid_aerosol(ibin)
        return
      endif

      if(xt .ge. 2.0 .or. xt .lt. 0.0)then
        j_index = jsulf_poor(jdum)
      else
        j_index = jsulf_rich(jdum)
      endif

      mdrh(ibin) = mdrh_t(j_index)

      if(ah2o_a(ibin)*100. .lt. mdrh(ibin)) then
        jaerosolstate(ibin) = all_solid
        jphase(ibin) = jsolid
        jhyst_leg(ibin) = jhyst_lo
        call adjust_solid_aerosol(ibin)
        return
      endif



10    call do_full_deliquescence(ibin)
      call mesa_ptc(ibin)	
      if (istat_mosaic_fe1 .lt. 0) return



      return
      end subroutine mesa


















      subroutine do_full_deliquescence(ibin)	



      integer ibin

      integer js





      do js = 1, nelectrolyte
       electrolyte(js,jsolid,ibin)  = 0.0
       electrolyte(js,jliquid,ibin) = electrolyte(js,jtotal,ibin)
      enddo


      electrolyte(jcaco3,jsolid,ibin) = electrolyte(jcaco3,jtotal,ibin)
      electrolyte(jcaso4,jsolid,ibin) = electrolyte(jcaso4,jtotal,ibin)
      electrolyte(jcaco3,jliquid,ibin)= 0.0
      electrolyte(jcaso4,jliquid,ibin)= 0.0




      aer(iso4_a,jsolid,ibin) = electrolyte(jcaso4,jsolid,ibin)
      aer(ino3_a,jsolid,ibin) = 0.0
      aer(icl_a, jsolid,ibin) = 0.0
      aer(inh4_a,jsolid,ibin) = 0.0
      aer(ioc_a, jsolid,ibin) = aer(ioc_a,jtotal,ibin)
      aer(imsa_a,jsolid,ibin) = 0.0
      aer(ico3_a,jsolid,ibin) = aer(ico3_a,jtotal,ibin)
      aer(ina_a, jsolid,ibin) = 0.0
      aer(ica_a, jsolid,ibin) = electrolyte(jcaco3,jsolid,ibin) +   &
                                electrolyte(jcaso4,jsolid,ibin)
      aer(ibc_a, jsolid,ibin) = aer(ibc_a,jtotal,ibin)
      aer(ioin_a,jsolid,ibin) = aer(ioin_a,jtotal,ibin)
      aer(ipcg1_b_c_a,jsolid,ibin)= aer(ipcg1_b_c_a,jtotal,ibin)
      aer(ipcg2_b_c_a,jsolid,ibin)= aer(ipcg2_b_c_a,jtotal,ibin)
      aer(ipcg3_b_c_a,jsolid,ibin)= aer(ipcg3_b_c_a,jtotal,ibin)
      aer(ipcg4_b_c_a,jsolid,ibin)= aer(ipcg4_b_c_a,jtotal,ibin)
      aer(ipcg5_b_c_a,jsolid,ibin)= aer(ipcg5_b_c_a,jtotal,ibin)
      aer(ipcg6_b_c_a,jsolid,ibin)= aer(ipcg6_b_c_a,jtotal,ibin)
      aer(ipcg7_b_c_a,jsolid,ibin)= aer(ipcg7_b_c_a,jtotal,ibin)
      aer(ipcg8_b_c_a,jsolid,ibin)= aer(ipcg8_b_c_a,jtotal,ibin)
      aer(ipcg9_b_c_a,jsolid,ibin)= aer(ipcg9_b_c_a,jtotal,ibin)
      aer(ipcg1_b_o_a,jsolid,ibin)= aer(ipcg1_b_o_a,jtotal,ibin)
      aer(ipcg2_b_o_a,jsolid,ibin)= aer(ipcg2_b_o_a,jtotal,ibin)
      aer(ipcg3_b_o_a,jsolid,ibin)= aer(ipcg3_b_o_a,jtotal,ibin)
      aer(ipcg4_b_o_a,jsolid,ibin)= aer(ipcg4_b_o_a,jtotal,ibin)
      aer(ipcg5_b_o_a,jsolid,ibin)= aer(ipcg5_b_o_a,jtotal,ibin)
      aer(ipcg6_b_o_a,jsolid,ibin)= aer(ipcg6_b_o_a,jtotal,ibin)
      aer(ipcg7_b_o_a,jsolid,ibin)= aer(ipcg7_b_o_a,jtotal,ibin)
      aer(ipcg8_b_o_a,jsolid,ibin)= aer(ipcg8_b_o_a,jtotal,ibin)
      aer(ipcg9_b_o_a,jsolid,ibin)= aer(ipcg9_b_o_a,jtotal,ibin)
      aer(iopcg1_b_c_a,jsolid,ibin)= aer(iopcg1_b_c_a,jtotal,ibin)
      aer(iopcg2_b_c_a,jsolid,ibin)= aer(iopcg2_b_c_a,jtotal,ibin)
      aer(iopcg3_b_c_a,jsolid,ibin)= aer(iopcg3_b_c_a,jtotal,ibin)
      aer(iopcg4_b_c_a,jsolid,ibin)= aer(iopcg4_b_c_a,jtotal,ibin)
      aer(iopcg5_b_c_a,jsolid,ibin)= aer(iopcg5_b_c_a,jtotal,ibin)
      aer(iopcg6_b_c_a,jsolid,ibin)= aer(iopcg6_b_c_a,jtotal,ibin)
      aer(iopcg7_b_c_a,jsolid,ibin)= aer(iopcg7_b_c_a,jtotal,ibin)
      aer(iopcg8_b_c_a,jsolid,ibin)= aer(iopcg8_b_c_a,jtotal,ibin)
      aer(iopcg1_b_o_a,jsolid,ibin)= aer(iopcg1_b_o_a,jtotal,ibin)
      aer(iopcg2_b_o_a,jsolid,ibin)= aer(iopcg2_b_o_a,jtotal,ibin)
      aer(iopcg3_b_o_a,jsolid,ibin)= aer(iopcg3_b_o_a,jtotal,ibin)
      aer(iopcg4_b_o_a,jsolid,ibin)= aer(iopcg4_b_o_a,jtotal,ibin)
      aer(iopcg5_b_o_a,jsolid,ibin)= aer(iopcg5_b_o_a,jtotal,ibin)
      aer(iopcg6_b_o_a,jsolid,ibin)= aer(iopcg6_b_o_a,jtotal,ibin)
      aer(iopcg7_b_o_a,jsolid,ibin)= aer(iopcg7_b_o_a,jtotal,ibin)
      aer(iopcg8_b_o_a,jsolid,ibin)= aer(iopcg8_b_o_a,jtotal,ibin)
      aer(ipcg1_f_c_a,jsolid,ibin)= aer(ipcg1_f_c_a,jtotal,ibin)
      aer(ipcg2_f_c_a,jsolid,ibin)= aer(ipcg2_f_c_a,jtotal,ibin)
      aer(ipcg3_f_c_a,jsolid,ibin)= aer(ipcg3_f_c_a,jtotal,ibin)
      aer(ipcg4_f_c_a,jsolid,ibin)= aer(ipcg4_f_c_a,jtotal,ibin)
      aer(ipcg5_f_c_a,jsolid,ibin)= aer(ipcg5_f_c_a,jtotal,ibin)
      aer(ipcg6_f_c_a,jsolid,ibin)= aer(ipcg6_f_c_a,jtotal,ibin)
      aer(ipcg7_f_c_a,jsolid,ibin)= aer(ipcg7_f_c_a,jtotal,ibin)
      aer(ipcg8_f_c_a,jsolid,ibin)= aer(ipcg8_f_c_a,jtotal,ibin)
      aer(ipcg9_f_c_a,jsolid,ibin)= aer(ipcg9_f_c_a,jtotal,ibin)
      aer(ipcg1_f_o_a,jsolid,ibin)= aer(ipcg1_f_o_a,jtotal,ibin)
      aer(ipcg2_f_o_a,jsolid,ibin)= aer(ipcg2_f_o_a,jtotal,ibin)
      aer(ipcg3_f_o_a,jsolid,ibin)= aer(ipcg3_f_o_a,jtotal,ibin)
      aer(ipcg4_f_o_a,jsolid,ibin)= aer(ipcg4_f_o_a,jtotal,ibin)
      aer(ipcg5_f_o_a,jsolid,ibin)= aer(ipcg5_f_o_a,jtotal,ibin)
      aer(ipcg6_f_o_a,jsolid,ibin)= aer(ipcg6_f_o_a,jtotal,ibin)
      aer(ipcg7_f_o_a,jsolid,ibin)= aer(ipcg7_f_o_a,jtotal,ibin)
      aer(ipcg8_f_o_a,jsolid,ibin)= aer(ipcg8_f_o_a,jtotal,ibin)
      aer(ipcg9_f_o_a,jsolid,ibin)= aer(ipcg9_f_o_a,jtotal,ibin)
      aer(iopcg1_f_c_a,jsolid,ibin)= aer(iopcg1_f_c_a,jtotal,ibin)
      aer(iopcg2_f_c_a,jsolid,ibin)= aer(iopcg2_f_c_a,jtotal,ibin)
      aer(iopcg3_f_c_a,jsolid,ibin)= aer(iopcg3_f_c_a,jtotal,ibin)
      aer(iopcg4_f_c_a,jsolid,ibin)= aer(iopcg4_f_c_a,jtotal,ibin)
      aer(iopcg5_f_c_a,jsolid,ibin)= aer(iopcg5_f_c_a,jtotal,ibin)
      aer(iopcg6_f_c_a,jsolid,ibin)= aer(iopcg6_f_c_a,jtotal,ibin)
      aer(iopcg7_f_c_a,jsolid,ibin)= aer(iopcg7_f_c_a,jtotal,ibin)
      aer(iopcg8_f_c_a,jsolid,ibin)= aer(iopcg8_f_c_a,jtotal,ibin)
      aer(iopcg1_f_o_a,jsolid,ibin)= aer(iopcg1_f_o_a,jtotal,ibin)
      aer(iopcg2_f_o_a,jsolid,ibin)= aer(iopcg2_f_o_a,jtotal,ibin)
      aer(iopcg3_f_o_a,jsolid,ibin)= aer(iopcg3_f_o_a,jtotal,ibin)
      aer(iopcg4_f_o_a,jsolid,ibin)= aer(iopcg4_f_o_a,jtotal,ibin)
      aer(iopcg5_f_o_a,jsolid,ibin)= aer(iopcg5_f_o_a,jtotal,ibin)
      aer(iopcg6_f_o_a,jsolid,ibin)= aer(iopcg6_f_o_a,jtotal,ibin)
      aer(iopcg7_f_o_a,jsolid,ibin)= aer(iopcg7_f_o_a,jtotal,ibin)
      aer(iopcg8_f_o_a,jsolid,ibin)= aer(iopcg8_f_o_a,jtotal,ibin)
      aer(ismpa_a,jsolid,ibin)= aer(ismpa_a,jtotal,ibin)
      aer(ismpbb_a,jsolid,ibin)= aer(ismpbb_a,jtotal,ibin)
      aer(iglysoa_r1_a,jsolid,ibin)= aer(iglysoa_r1_a,jtotal,ibin)
      aer(iglysoa_r2_a,jsolid,ibin)= aer(iglysoa_r2_a,jtotal,ibin)
      aer(iglysoa_sfc_a,jsolid,ibin)= aer(iglysoa_sfc_a,jtotal,ibin)
      aer(iglysoa_nh4_a,jsolid,ibin)= aer(iglysoa_nh4_a,jtotal,ibin)
      aer(iglysoa_oh_a,jsolid,ibin)= aer(iglysoa_oh_a,jtotal,ibin)
      aer(iant1_c_a,jsolid,ibin)= aer(iant1_c_a,jtotal,ibin)
      aer(iant2_c_a,jsolid,ibin)= aer(iant2_c_a,jtotal,ibin)
      aer(iant3_c_a,jsolid,ibin)= aer(iant3_c_a,jtotal,ibin)
      aer(iant4_c_a,jsolid,ibin)= aer(iant4_c_a,jtotal,ibin)
      aer(iant1_o_a,jsolid,ibin)= aer(iant1_o_a,jtotal,ibin)
      aer(iant2_o_a,jsolid,ibin)= aer(iant2_o_a,jtotal,ibin)
      aer(iant3_o_a,jsolid,ibin)= aer(iant3_o_a,jtotal,ibin)
      aer(iant4_o_a,jsolid,ibin)= aer(iant4_o_a,jtotal,ibin)
      aer(ibiog1_c_a,jsolid,ibin)= aer(ibiog1_c_a,jtotal,ibin)
      aer(ibiog2_c_a,jsolid,ibin)= aer(ibiog2_c_a,jtotal,ibin)
      aer(ibiog3_c_a,jsolid,ibin)= aer(ibiog3_c_a,jtotal,ibin)
      aer(ibiog4_c_a,jsolid,ibin)= aer(ibiog4_c_a,jtotal,ibin)
      aer(ibiog1_o_a,jsolid,ibin)= aer(ibiog1_o_a,jtotal,ibin)
      aer(ibiog2_o_a,jsolid,ibin)= aer(ibiog2_o_a,jtotal,ibin)
      aer(ibiog3_o_a,jsolid,ibin)= aer(ibiog3_o_a,jtotal,ibin)
      aer(ibiog4_o_a,jsolid,ibin)= aer(ibiog4_o_a,jtotal,ibin)
      aer(iasoaX_a,jsolid,ibin)= aer(iasoaX_a,jtotal,ibin)
      aer(iasoa1_a,jsolid,ibin)= aer(iasoa1_a,jtotal,ibin)
      aer(iasoa2_a,jsolid,ibin)= aer(iasoa2_a,jtotal,ibin)
      aer(iasoa3_a,jsolid,ibin)= aer(iasoa3_a,jtotal,ibin)
      aer(iasoa4_a,jsolid,ibin)= aer(iasoa4_a,jtotal,ibin)
      aer(ibsoaX_a,jsolid,ibin)= aer(ibsoaX_a,jtotal,ibin)
      aer(ibsoa1_a,jsolid,ibin)= aer(ibsoa1_a,jtotal,ibin)
      aer(ibsoa2_a,jsolid,ibin)= aer(ibsoa2_a,jtotal,ibin)
      aer(ibsoa3_a,jsolid,ibin)= aer(ibsoa3_a,jtotal,ibin)
      aer(ibsoa4_a,jsolid,ibin)= aer(ibsoa4_a,jtotal,ibin)


      aer(iso4_a,jliquid,ibin) = aer(iso4_a,jtotal,ibin) -   &
                                 electrolyte(jcaso4,jsolid,ibin)
      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jtotal,ibin)
      aer(icl_a, jliquid,ibin) = aer(icl_a,jtotal,ibin)
      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jtotal,ibin)
      aer(ioc_a, jliquid,ibin) = 0.0
      aer(imsa_a,jliquid,ibin) = aer(imsa_a,jtotal,ibin)
      aer(ico3_a,jliquid,ibin) = 0.0
      aer(ina_a, jliquid,ibin) = aer(ina_a,jtotal,ibin)
      aer(ica_a, jliquid,ibin) = electrolyte(jcano3,jtotal,ibin) +   &
                                 electrolyte(jcacl2,jtotal,ibin)
      aer(ibc_a, jliquid,ibin) = 0.0
      aer(ioin_a,jliquid,ibin) = 0.0
      aer(ipcg1_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg2_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg3_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg4_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg5_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg6_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg7_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg8_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg9_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg1_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg2_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg3_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg4_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg5_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg6_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg7_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg8_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg9_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg1_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg2_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg3_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg4_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg5_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg6_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg7_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg8_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg1_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg2_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg3_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg4_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg5_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg6_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg7_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg8_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg1_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg2_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg3_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg4_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg5_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg6_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg7_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg8_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg9_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg1_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg2_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg3_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg4_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg5_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg6_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg7_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg8_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg9_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg1_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg2_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg3_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg4_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg5_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg6_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg7_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg8_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg1_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg2_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg3_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg4_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg5_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg6_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg7_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg8_f_o_a,jliquid,ibin)= 0.0
      aer(ismpa_a,jliquid,ibin)= 0.0
      aer(ismpbb_a,jliquid,ibin)= 0.0
      aer(iglysoa_r1_a,jliquid,ibin)= 0.0
      aer(iglysoa_r2_a,jliquid,ibin)= 0.0
      aer(iglysoa_sfc_a,jliquid,ibin)= 0.0
      aer(iglysoa_nh4_a,jliquid,ibin)= 0.0
      aer(iglysoa_oh_a,jliquid,ibin)= 0.0
      aer(iant1_c_a,jliquid,ibin)= 0.0
      aer(iant2_c_a,jliquid,ibin)= 0.0
      aer(iant3_c_a,jliquid,ibin)= 0.0
      aer(iant4_c_a,jliquid,ibin)= 0.0
      aer(iant1_o_a,jliquid,ibin)= 0.0
      aer(iant2_o_a,jliquid,ibin)= 0.0
      aer(iant3_o_a,jliquid,ibin)= 0.0
      aer(iant4_o_a,jliquid,ibin)= 0.0
      aer(ibiog1_c_a,jliquid,ibin)= 0.0
      aer(ibiog2_c_a,jliquid,ibin)= 0.0
      aer(ibiog3_c_a,jliquid,ibin)= 0.0
      aer(ibiog4_c_a,jliquid,ibin)= 0.0
      aer(ibiog1_o_a,jliquid,ibin)= 0.0
      aer(ibiog2_o_a,jliquid,ibin)= 0.0
      aer(ibiog3_o_a,jliquid,ibin)= 0.0
      aer(ibiog4_o_a,jliquid,ibin)= 0.0
      aer(iasoaX_a,jliquid,ibin)= 0.0
      aer(iasoa1_a,jliquid,ibin)= 0.0
      aer(iasoa2_a,jliquid,ibin)= 0.0
      aer(iasoa3_a,jliquid,ibin)= 0.0
      aer(iasoa4_a,jliquid,ibin)= 0.0
      aer(ibsoaX_a,jliquid,ibin)= 0.0
      aer(ibsoa1_a,jliquid,ibin)= 0.0
      aer(ibsoa2_a,jliquid,ibin)= 0.0
      aer(ibsoa3_a,jliquid,ibin)= 0.0
      aer(ibsoa4_a,jliquid,ibin)= 0.0





      return
      end subroutine do_full_deliquescence































      subroutine mesa_ptc(ibin)		



      integer ibin

      integer iaer, iconverge, iconverge_flux, iconverge_mass,   &
           idissolved, itdum, js, je, jp			
      real(kind=8) tau_p(nsalt), tau_d(nsalt)
      real(kind=8) hsalt_min
      real(kind=8) phi_prod, alpha_fac, sum_dum		
      real(kind=8) aer_H






      itdum = 0		
      hsalt_max = 1.e25



      do js = 1, nsalt
        hsalt(js)     = 0.0
        sat_ratio(js) = 0.0
        phi_salt(js)  = 0.0
        flux_sl(js)   = 0.0
      enddo



      sum_dum = 0.0
      do je = 1, nelectrolyte
        sum_dum = sum_dum + electrolyte(je,jtotal,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      do je = 1, nelectrolyte
        epercent(je,jtotal,ibin) = 100.*electrolyte(je,jtotal,ibin)/sum_dum
      enddo



      do js = 1, nsalt
        jsalt_present(js) = 0			
        if(epercent(js,jtotal,ibin) .gt. 1.0)then
          jsalt_present(js) = 1			
        endif
      enddo


      mass_dry_a(ibin) = 0.0

      aer_H = (2.*aer(iso4_a,jtotal,ibin) +  &
                  aer(ino3_a,jtotal,ibin) +  &
                  aer(icl_a,jtotal,ibin)  +  &
                  aer(imsa_a,jtotal,ibin) +  &
               2.*aer(ico3_a,jtotal,ibin))-  &
              (2.*aer(ica_a,jtotal,ibin)  +  &
                  aer(ina_a,jtotal,ibin)  +  &
                  aer(inh4_a,jtotal,ibin))
      aer_H = max(aer_H, 0.0d0)		

      do iaer = 1, naer
       mass_dry_a(ibin) = mass_dry_a(ibin) +  &
          aer(iaer,jtotal,ibin)*mw_aer_mac(iaer) 	
        vol_dry_a(ibin)  = vol_dry_a(ibin) +  &
          aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  	
      enddo
      mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
      vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

      mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15			
      vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15				

      mass_dry_salt(ibin) = 0.0		
      do je = 1, nsalt
        mass_dry_salt(ibin) = mass_dry_salt(ibin) +  &
              electrolyte(je,jtotal,ibin)*mw_electrolyte(je)*1.e-15	
      enddo

      nmesa_call = nmesa_call + 1



      do 500 itdum = 1, nmax_mesa



      call mesa_flux_salt(ibin)
      if (istat_mosaic_fe1 .lt. 0) return



      call mesa_convergence_criterion(ibin,      &
                                      iconverge_mass,   &
                                      iconverge_flux,   &
                                      idissolved)

      if(iconverge_mass .eq. myes)then
        iter_mesa(ibin) = iter_mesa(ibin) + itdum
        niter_mesa = niter_mesa + itdum
        niter_mesa_max = max(niter_mesa_max, itdum)
        jaerosolstate(ibin) = all_solid
        call adjust_solid_aerosol(ibin)
        jhyst_leg(ibin) = jhyst_lo
        growth_factor(ibin) = 1.0
        return
      elseif(iconverge_flux .eq. myes)then
        iter_mesa(ibin) = iter_mesa(ibin)+ itdum
        niter_mesa = niter_mesa + itdum
        niter_mesa_max = max(niter_mesa_max, itdum)
        mass_wet_a(ibin)    = mass_dry_a(ibin) + water_a(ibin)*1.e-3	
        vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3		
        growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)		

        if(idissolved .eq. myes)then
          jaerosolstate(ibin) = all_liquid

        else
          jaerosolstate(ibin) = mixed
          jhyst_leg(ibin) = jhyst_lo
        endif


        sum_dum = 0.0
        jp = jsolid
        do je = 1, nelectrolyte
          electrolyte(je,jp,ibin) = max(0.D0,electrolyte(je,jp,ibin)) 
          sum_dum = sum_dum + electrolyte(je,jp,ibin)
        enddo
        electrolyte_sum(jp,ibin) = sum_dum
        if(sum_dum .eq. 0.)sum_dum = 1.0
        do je = 1, nelectrolyte
          epercent(je,jp,ibin) = 100.*electrolyte(je,jp,ibin)/sum_dum
        enddo

        return
      endif



      hsalt_min = 1.e25
      do js = 1, nsalt

        phi_prod = phi_salt(js) * phi_salt_old(js)

        if(itdum .gt. 1 .and. phi_prod .gt. 0.0)then
          phi_bar(js) = (abs(phi_salt(js))-abs(phi_salt_old(js)))/   &
                                    alpha_salt(js)
        else
          phi_bar(js) = 0.0			
        endif

        if(phi_bar(js) .lt. 0.0)then		
          phi_bar(js) = max(phi_bar(js), -10.0D0)
          alpha_fac = 3.0*exp(phi_bar(js))
          alpha_salt(js) = min(alpha_fac*abs(phi_salt(js)), 0.9D0)
        elseif(phi_bar(js) .gt. 0.0)then	
           alpha_salt(js) = min(abs(phi_salt(js)), 0.5D0)
        else					
           alpha_salt(js) = min(abs(phi_salt(js))/3.0, 0.5D0)
        endif



        phi_salt_old(js) = phi_salt(js)		


        if(flux_sl(js) .gt. 0.)then

          tau_p(js) = eleliquid(js)/flux_sl(js)	
          if(tau_p(js) .eq. 0.0)then
            hsalt(js) = 1.e25
            flux_sl(js) = 0.0
            phi_salt(js)= 0.0
          else
            hsalt(js) = alpha_salt(js)*tau_p(js)
          endif

        elseif(flux_sl(js) .lt. 0.)then

          tau_p(js) = -eleliquid(js)/flux_sl(js)	
          tau_d(js) = -electrolyte(js,jsolid,ibin)/flux_sl(js) 
          if(tau_p(js) .eq. 0.0)then
            hsalt(js) = alpha_salt(js)*tau_d(js)
          else
            hsalt(js) = alpha_salt(js)*min(tau_p(js),tau_d(js))
          endif

        else

          hsalt(js) = 1.e25

        endif

          hsalt_min = min(hsalt(js), hsalt_min)

      enddo




      do js = 1, nsalt
        electrolyte(js,jsolid,ibin) =    &
                         electrolyte(js,jsolid,ibin)  +   &
                         hsalt(js) * flux_sl(js)
      enddo



      call electrolytes_to_ions(jsolid,ibin)



      do iaer = 1, naer
        aer(iaer,jliquid,ibin) = aer(iaer,jtotal,ibin) -   &
                                       aer(iaer,jsolid,ibin)
      enddo





500   continue	

      nmesa_fail = nmesa_fail + 1
      iter_mesa(ibin) = iter_mesa(ibin) + itdum
      niter_mesa = niter_mesa + itdum
      jaerosolstate(ibin) = mixed
      jhyst_leg(ibin) = jhyst_lo
      mass_wet_a(ibin)    = mass_dry_a(ibin) + water_a(ibin)*1.e-3	
      vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3		
      growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)		

      return
      end subroutine mesa_ptc
















      subroutine mesa_flux_salt(ibin)	



      integer ibin

      integer js, je						
      real(kind=8) xt, calcium, sum_salt, sum_dum	



      call ions_to_electrolytes(jliquid,ibin,xt)
      if (istat_mosaic_fe1 .lt. 0) return
      call compute_activities(ibin)
      activity(jna3hso4,ibin)   = 0.0

      if(water_a(ibin) .le. 0.0)then
        do js = 1, nsalt
         flux_sl(js) = 0.0
        enddo
        return
      endif


      call mesa_estimate_eleliquid(ibin,xt)

      calcium = aer(ica_a,jliquid,ibin)




      sum_dum = 0.0
      do je = 1, nelectrolyte
        sum_dum = sum_dum + electrolyte(je,jliquid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      do je = 1, nelectrolyte
        epercent(je,jliquid,ibin) = 100.*electrolyte(je,jliquid,ibin)/sum_dum
      enddo




      sum_salt = 0.0
      do js = 1, nsalt
        sum_salt = sum_salt + electrolyte(js,jsolid,ibin)
      enddo
      electrolyte_sum(jsolid,ibin) = sum_salt
      if(sum_salt .eq. 0.0)sum_salt = 1.0
      do js = 1, nsalt
        frac_salt_solid(js) = electrolyte(js,jsolid,ibin)/sum_salt
        frac_salt_liq(js)   = epercent(js,jliquid,ibin)/100.
      enddo




      do js = 1, nsalt		


        sat_ratio(js) = activity(js,ibin)/keq_sl(js)

        phi_salt(js)  = (sat_ratio(js) - 1.0)/max(sat_ratio(js),1.0D0)


        if(sat_ratio(js)       .lt. 1.00 .and.   &
           frac_salt_solid(js) .lt. 0.01 .and.   &
           frac_salt_solid(js) .gt. 0.0)then
          call mesa_dissolve_small_salt(ibin,js)
          call mesa_estimate_eleliquid(ibin,xt)
          sat_ratio(js) = activity(js,ibin)/keq_sl(js)
        endif


        flux_sl(js) = sat_ratio(js) - 1.0


        if( (sat_ratio(js)               .lt. 1.0 .and.   &
             electrolyte(js,jsolid,ibin) .eq. 0.0) .or.   &
            (calcium .gt. 0.0 .and. frac_salt_liq(js).lt.0.01).or.   &
            (calcium .gt. 0.0 .and. jsalt_present(js).eq.0) )then
          flux_sl(js) = 0.0
          phi_salt(js)= 0.0
        endif

      enddo



      sat_ratio(jcano3) = 1.0
      phi_salt(jcano3)  = 0.0
      flux_sl(jcano3)   = 0.0

      sat_ratio(jcacl2) = 1.0
      phi_salt(jcacl2)  = 0.0
      flux_sl(jcacl2)   = 0.0


      return
      end subroutine mesa_flux_salt






















      subroutine mesa_estimate_eleliquid(ibin,xt)	



      integer ibin, jp
      real(kind=8) xt

      integer iaer, je, jc, ja, icase
      real(kind=8) store(naer), sum_dum, sum_naza, sum_nczc, sum_na_nh4,   &
           f_nh4, f_na, xh, xb, xl, xs, xt_d, xna_d, xnh4_d,   &
           xdum, dum, cat_net
      real(kind=8) nc(ncation), na(nanion)
      real(kind=8) dum_ca, dum_no3, dum_cl, cano3, cacl2




      do iaer =  1, naer
      aer(iaer,jliquid,ibin) = max(0.0D0, aer(iaer,jliquid,ibin))
      enddo



      call calculate_xt(ibin,jliquid,xt)

      if(xt .ge. 2.0 .or. xt.lt.0.)then
       icase = 1	
      else
       icase = 2	
      endif



      do je = 1, nelectrolyte
        eleliquid(je) = 0.0
      enddo




      jp = jliquid

      if(icase.eq.1)then 

        dum_ca  = aer(ica_a,jp,ibin)
        dum_no3 = aer(ino3_a,jp,ibin)
        dum_cl  = aer(icl_a,jp,ibin)

        cano3   = min(dum_ca, 0.5*dum_no3)
        dum_ca  = max(0.D0, dum_ca - cano3)
        dum_no3 = max(0.D0, dum_no3 - 2.*cano3)

        cacl2   = min(dum_ca, 0.5*dum_cl)
        dum_ca  = max(0.D0, dum_ca - cacl2)
        dum_cl  = max(0.D0, dum_cl - 2.*cacl2)

        na(ja_hso4)= 0.0
        na(ja_so4) = aer(iso4_a,jp,ibin)
        na(ja_no3) = aer(ino3_a,jp,ibin)
        na(ja_cl)  = aer(icl_a, jp,ibin)
        na(ja_msa) = aer(imsa_a,jp,ibin)

        nc(jc_ca)  = aer(ica_a, jp,ibin)
        nc(jc_na)  = aer(ina_a, jp,ibin)
        nc(jc_nh4) = aer(inh4_a,jp,ibin)

        cat_net =     &
            ( 2.d0*na(ja_so4)+na(ja_no3)+na(ja_cl)+na(ja_msa) ) -  &
            ( nc(jc_h)+2.d0*nc(jc_ca) +nc(jc_nh4)+nc(jc_na) )

        if(cat_net .lt. 0.0)then

          nc(jc_h) = 0.0

        else  

          nc(jc_h) = cat_net

        endif



      sum_naza = 0.0
      do ja = 1, nanion
        sum_naza = sum_naza + na(ja)*za(ja)
      enddo

      sum_nczc = 0.0
      do jc = 1, ncation
        sum_nczc = sum_nczc + nc(jc)*zc(jc)
      enddo

      if(sum_naza .eq. 0. .or. sum_nczc .eq. 0.)then
        if (iprint_mosaic_diag1 .gt. 0) then
          write(6,*)'subroutine mesa_estimate_eleliquid'
          write(6,*)'ionic concentrations are zero'
          write(6,*)'sum_naza = ', sum_naza
          write(6,*)'sum_nczc = ', sum_nczc
        endif
        return
      endif

      do ja = 1, nanion
        xeq_a(ja) = na(ja)*za(ja)/sum_naza
      enddo

      do jc = 1, ncation
        xeq_c(jc) = nc(jc)*zc(jc)/sum_nczc
      enddo

      na_ma(ja_so4) = na(ja_so4) *mw_a(ja_so4)
      na_ma(ja_no3) = na(ja_no3) *mw_a(ja_no3)
      na_ma(ja_cl)  = na(ja_cl)  *mw_a(ja_cl)
      na_ma(ja_hso4)= na(ja_hso4)*mw_a(ja_hso4)
      na_Ma(ja_msa) = na(ja_msa) *MW_a(ja_msa)

      nc_mc(jc_ca)  = nc(jc_ca) *mw_c(jc_ca)
      nc_mc(jc_na)  = nc(jc_na) *mw_c(jc_na)
      nc_mc(jc_nh4) = nc(jc_nh4)*mw_c(jc_nh4)
      nc_mc(jc_h)   = nc(jc_h)  *mw_c(jc_h)



      eleliquid(jna2so4) = (xeq_c(jc_na) *na_ma(ja_so4) +  &
                            xeq_a(ja_so4)*nc_mc(jc_na))/   &
                             mw_electrolyte(jna2so4)

      eleliquid(jnahso4) = (xeq_c(jc_na) *na_ma(ja_hso4) +  &
                            xeq_a(ja_hso4)*nc_mc(jc_na))/   &
                             mw_electrolyte(jnahso4)

      eleliquid(jnamsa)  = (xeq_c(jc_na) *na_ma(ja_msa) + &
                            xeq_a(ja_msa)*nc_mc(jc_na))/  &
                             mw_electrolyte(jnamsa)

      eleliquid(jnano3)  = (xeq_c(jc_na) *na_ma(ja_no3) +  &
                            xeq_a(ja_no3)*nc_mc(jc_na))/   &
                             mw_electrolyte(jnano3)

      eleliquid(jnacl)   = (xeq_c(jc_na) *na_ma(ja_cl) +   &
                            xeq_a(ja_cl) *nc_mc(jc_na))/   &
                             mw_electrolyte(jnacl)

      eleliquid(jnh4so4) = (xeq_c(jc_nh4)*na_ma(ja_so4) +   &
                            xeq_a(ja_so4)*nc_mc(jc_nh4))/   &
                             mw_electrolyte(jnh4so4)

      eleliquid(jnh4hso4)= (xeq_c(jc_nh4)*na_ma(ja_hso4) +   &
                            xeq_a(ja_hso4)*nc_mc(jc_nh4))/   &
                             mw_electrolyte(jnh4hso4)

      eleliquid(jnh4msa) = (xeq_c(jc_nh4) *na_ma(ja_msa) +  &
                            xeq_a(ja_msa)*nc_mc(jc_nh4))/   &
                             mw_electrolyte(jnh4msa)

      eleliquid(jnh4no3) = (xeq_c(jc_nh4)*na_ma(ja_no3) +   &
                            xeq_a(ja_no3)*nc_mc(jc_nh4))/   &
                             mw_electrolyte(jnh4no3)

      eleliquid(jnh4cl)  = (xeq_c(jc_nh4)*na_ma(ja_cl) +   &
                            xeq_a(ja_cl) *nc_mc(jc_nh4))/  &
                             mw_electrolyte(jnh4cl)

      eleliquid(jcano3)  = (xeq_c(jc_ca) *na_ma(ja_no3) +  &
                            xeq_a(ja_no3)*nc_mc(jc_ca))/   &
                             mw_electrolyte(jcano3)

      eleliquid(jcamsa2) = (xeq_c(jc_ca) *na_ma(ja_msa) +  &
                            xeq_a(ja_msa)*nc_mc(jc_ca))/   &
                             mw_electrolyte(jcamsa2)

      eleliquid(jcacl2)  = (xeq_c(jc_ca) *na_ma(ja_cl) +   &
                            xeq_a(ja_cl) *nc_mc(jc_ca))/   &
                             mw_electrolyte(jcacl2)

      eleliquid(jh2so4)  = (xeq_c(jc_h)  *na_ma(ja_hso4) + &
                            xeq_a(ja_hso4)*nc_mc(jc_h))/   &
                             mw_electrolyte(jh2so4)

      eleliquid(jhno3)   = (xeq_c(jc_h)  *na_ma(ja_no3) +  &
                            xeq_a(ja_no3)*nc_mc(jc_h))/    &
                             mw_electrolyte(jhno3)

      eleliquid(jhcl)    = (xeq_c(jc_h) *na_ma(ja_cl) +   &
                            xeq_a(ja_cl)*nc_mc(jc_h))/    &
                             mw_electrolyte(jhcl)

      eleliquid(jmsa)    = (xeq_c(jc_h)  *na_ma(ja_msa) + &
                            xeq_a(ja_msa)*nc_mc(jc_h))/   &
                             mw_electrolyte(jmsa)



      elseif(icase.eq.2)then 

        jp = jliquid

        store(iso4_a) = aer(iso4_a,jp,ibin)
        store(imsa_a) = aer(imsa_a,jp,ibin)
        store(inh4_a) = aer(inh4_a,jp,ibin)
        store(ina_a)  = aer(ina_a, jp,ibin)
        store(ica_a)  = aer(ica_a, jp,ibin)

        call form_camsa2(store,jp,ibin)

        sum_na_nh4 = store(ina_a) + store(inh4_a)
        if(sum_na_nh4 .gt. 0.0)then
          f_nh4 = store(inh4_a)/sum_na_nh4
          f_na  = store(ina_a)/sum_na_nh4
        else
          f_nh4 = 0.0
          f_na  = 0.0
        endif


        if(sum_na_nh4 .gt. store(imsa_a))then
          eleliquid(jnh4msa) = f_nh4*store(imsa_a)
          eleliquid(jnamsa)  = f_na *store(imsa_a)
          store(inh4_a)= store(inh4_a)-eleliquid(jnh4msa) 
          store(ina_a) = store(ina_a) -eleliquid(jnamsa)  
        else
          eleliquid(jnh4msa) = store(inh4_a)
          eleliquid(jnamsa)  = store(ina_a)
          eleliquid(jmsa)    = store(imsa_a) - sum_na_nh4
          store(inh4_a)= 0.0  
          store(ina_a) = 0.0  
        endif

        if(store(iso4_a).eq.0.0)goto 10

        xt_d  = xt
        xna_d = 1. + 0.5*aer(ina_a,jp,ibin)/aer(iso4_a,jp,ibin)
        xdum = aer(iso4_a,jp,ibin) - aer(inh4_a,jp,ibin)

        dum = 2.d0*aer(iso4_a,jp,ibin) - aer(ina_a,jp,ibin)
        if(aer(inh4_a,jp,ibin) .gt. 0.0 .and. dum .gt. 0.0)then
          xnh4_d = 2.*aer(inh4_a,jp,ibin)/   &
                  (2.*aer(iso4_a,jp,ibin) - aer(ina_a,jp,ibin))
        else
          xnh4_d = 0.0
        endif


        if(aer(inh4_a,jp,ibin) .gt. 0.0)then


        if(xt_d .ge. xna_d)then
          eleliquid(jna2so4) = 0.5*aer(ina_a,jp,ibin)

          if(xnh4_d .ge. 5./3.)then
            eleliquid(jnh4so4) = 1.5*aer(ina_a,jp,ibin)   &
                               - 3.*xdum - aer(inh4_a,jp,ibin)
            eleliquid(jlvcite) = 2.*xdum + aer(inh4_a,jp,ibin)   &
                               - aer(ina_a,jp,ibin)
          elseif(xnh4_d .ge. 1.5)then
            eleliquid(jnh4so4) = aer(inh4_a,jp,ibin)/5.
            eleliquid(jlvcite) = aer(inh4_a,jp,ibin)/5.
          elseif(xnh4_d .ge. 1.0)then
            eleliquid(jnh4so4) = aer(inh4_a,jp,ibin)/6.
            eleliquid(jlvcite) = aer(inh4_a,jp,ibin)/6.
            eleliquid(jnh4hso4)= aer(inh4_a,jp,ibin)/6.
          endif

        elseif(xt_d .gt. 1.0)then
          eleliquid(jnh4so4)  = aer(inh4_a,jp,ibin)/6.
          eleliquid(jlvcite)  = aer(inh4_a,jp,ibin)/6.
          eleliquid(jnh4hso4) = aer(inh4_a,jp,ibin)/6.
          eleliquid(jna2so4)  = aer(ina_a,jp,ibin)/3.
          eleliquid(jnahso4)  = aer(ina_a,jp,ibin)/3.
        elseif(xt_d .le. 1.0)then
          eleliquid(jna2so4)  = aer(ina_a,jp,ibin)/4.
          eleliquid(jnahso4)  = aer(ina_a,jp,ibin)/2.
          eleliquid(jlvcite)  = aer(inh4_a,jp,ibin)/6.
          eleliquid(jnh4hso4) = aer(inh4_a,jp,ibin)/2.
        endif

        else

        if(xt_d .gt. 1.0)then
          eleliquid(jna2so4) = aer(ina_a,jp,ibin) - aer(iso4_a,jp,ibin)
          eleliquid(jnahso4) = 2.*aer(iso4_a,jp,ibin) -   &
                                  aer(ina_a,jp,ibin)
        else
          eleliquid(jna2so4) = aer(ina_a,jp,ibin)/4.
          eleliquid(jnahso4) = aer(ina_a,jp,ibin)/2.
        endif


        endif



      endif


10    return
      end subroutine mesa_estimate_eleliquid
















      subroutine mesa_dissolve_small_salt(ibin,js)



      integer ibin, js, jp

      jp = jsolid


      if(js .eq. jnh4so4)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                           2.*electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                         2.*electrolyte(jnh4so4,jp,ibin) +   &
                         3.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jnh4msa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jna2so4,jp,ibin) +   &
                         2.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnh4so4,jp,ibin) +   &
                         2.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jlvcite)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                           3.*electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                           2.*electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                         2.*electrolyte(jnh4so4,jp,ibin) +   &
                         3.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jnh4msa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jna2so4,jp,ibin) +   &
                         2.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnh4so4,jp,ibin) +   &
                         2.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jnh4hso4)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                             electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                         2.*electrolyte(jnh4so4,jp,ibin) +   &
                         3.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jnh4msa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jna2so4,jp,ibin) +   &
                         2.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnh4so4,jp,ibin) +   &
                         2.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jna2so4)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                           2.*electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                            electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jna2so4,jp,ibin) +   &
                         3.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jna2so4,jp,ibin) +   &
                         2.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnh4so4,jp,ibin) +   &
                         2.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jna3hso4)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                           3.*electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                           2.*electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                            electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jna2so4,jp,ibin) +   &
                         3.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jna2so4,jp,ibin) +   &
                         2.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnh4so4,jp,ibin) +   &
                         2.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jnahso4)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                            electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jna2so4,jp,ibin) +   &
                         3.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jna2so4,jp,ibin) +   &
                         2.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnh4so4,jp,ibin) +   &
                         2.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jnh4no3)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                         2.*electrolyte(jnh4so4,jp,ibin) +   &
                         3.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jnh4msa,jp,ibin)

        aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
                         2.*electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jhno3,jp,ibin)
        return
      endif


      if(js .eq. jnh4cl)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                         2.*electrolyte(jnh4so4,jp,ibin) +   &
                         3.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jnh4msa,jp,ibin)

        aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                            electrolyte(jhcl,jp,ibin)
        return
      endif


      if(js .eq. jnano3)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                            electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jna2so4,jp,ibin) +   &
                         3.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)

        aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
                         2.*electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jhno3,jp,ibin)
        return
      endif


      if(js .eq. jnacl)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                            electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jna2so4,jp,ibin) +   &
                         3.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)

        aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                            electrolyte(jhcl,jp,ibin)
        return
      endif


      if(js .eq. jcano3)then
        aer(ica_a,jliquid,ibin)  = aer(ica_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
                            2.*electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jcaco3,jp,ibin)  +   &
                            electrolyte(jcamsa2,jp,ibin)

        aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
                         2.*electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jhno3,jp,ibin)
        return
      endif


      if(js .eq. jcacl2)then
        aer(ica_a,jliquid,ibin) = aer(ica_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(icl_a,jliquid,ibin) = aer(icl_a,jliquid,ibin) +   &
                            2.*electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jcaco3,jp,ibin)  +   &
                            electrolyte(jcamsa2,jp,ibin)

        aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                            electrolyte(jhcl,jp,ibin)
        return
      endif



      return
      end subroutine mesa_dissolve_small_salt













      subroutine mesa_convergence_criterion(ibin,  &  
                                       iconverge_mass,    &
                                       iconverge_flux,    &
                                       idissolved)



      integer ibin, iconverge_mass, iconverge_flux, idissolved

      integer je, js, iaer
      real(kind=8) mass_solid, mass_solid_salt, frac_solid, xt, h_ion, &
           crustal_solids, sumflux


      idissolved = mno		


      iconverge_mass = mno	








      mass_solid_salt = 0.0
      do je = 1, nsalt
        mass_solid_salt = mass_solid_salt + &
             electrolyte(je,jsolid,ibin)*mw_electrolyte(je)*1.e-15	
      enddo














      if( mass_dry_salt(ibin) < 1e-30 ) then
         iconverge_mass = myes
         return
      else
         frac_solid = mass_solid_salt/mass_dry_salt(ibin)         
         if(frac_solid .ge. 0.98)then
            iconverge_mass = myes
            return
         endif
      end if



      iconverge_flux = myes
      do js = 1, nsalt
        if(abs(phi_salt(js)).gt. rtol_mesa)then
          iconverge_flux = mno
          return
        endif
      enddo





      sumflux = 0.0
      do js = 1, nsalt
        sumflux = sumflux + abs(flux_sl(js))
      enddo





      crustal_solids = electrolyte(jcaco3,jsolid,ibin)*mw_electrolyte(jcaco3) +  &
                       electrolyte(jcaso4,jsolid,ibin)*mw_electrolyte(jcaso4) +  &
                       aer(ioin_a,jsolid,ibin)*mw_aer_mac(ioin_a)


      if ( sumflux .eq. 0.0 .and. &
           crustal_solids .le. xhyst_up_crustal_thresh*(mass_dry_a(ibin)*1.0e15) ) then
         
        idissolved = myes
      endif



      return
      end subroutine mesa_convergence_criterion














      subroutine adjust_solid_aerosol(ibin)



      integer ibin

      integer iaer, je


      jphase(ibin)    = jsolid
      jhyst_leg(ibin) = jhyst_lo	
      water_a(ibin)   = 0.0


      do iaer = 1, naer
        aer(iaer, jsolid, ibin) = aer(iaer,jtotal,ibin)
        aer(iaer, jliquid,ibin) = 0.0
      enddo


      do je = 1, nelectrolyte
        electrolyte(je,jliquid,ibin) = 0.0
        epercent(je,jliquid,ibin)    = 0.0
        electrolyte(je,jsolid,ibin)  = electrolyte(je,jtotal,ibin)
        epercent(je,jsolid,ibin)     = epercent(je,jtotal,ibin)
      enddo


      aer(inh4_a,jtotal,ibin) = aer(inh4_a,jsolid,ibin)
      aer(ino3_a,jtotal,ibin) = aer(ino3_a,jsolid,ibin)
      aer(icl_a,jtotal,ibin)  = aer(icl_a,jsolid,ibin)


      do je = 1, nelectrolyte
        electrolyte(je,jtotal,ibin) = electrolyte(je,jsolid,ibin)
        epercent(je,jtotal,ibin)    = epercent(je,jsolid,ibin)
      enddo

      return
      end subroutine adjust_solid_aerosol















      subroutine adjust_liquid_aerosol(ibin)



      integer ibin

      integer je




      jphase(ibin)    = jliquid
      jhyst_leg(ibin) = jhyst_up	


      do je = 1, nelectrolyte
        electrolyte(je,jsolid,ibin)  = 0.0
        epercent(je,jsolid,ibin)     = 0.0
        electrolyte(je,jliquid,ibin) = electrolyte(je,jtotal,ibin)
        epercent(je,jliquid,ibin)    = epercent(je,jtotal,ibin)
      enddo

      electrolyte(jcaco3,jsolid,ibin) = electrolyte(jcaco3,jtotal,ibin)
      electrolyte(jcaso4,jsolid,ibin) = electrolyte(jcaso4,jtotal,ibin)
      epercent(jcaco3,jsolid,ibin)    = epercent(jcaco3,jtotal,ibin)
      epercent(jcaso4,jsolid,ibin)    = epercent(jcaso4,jtotal,ibin)
      electrolyte(jcaco3,jliquid,ibin)= 0.0
      electrolyte(jcaso4,jliquid,ibin)= 0.0
      epercent(jcaco3,jliquid,ibin)   = 0.0
      epercent(jcaso4,jliquid,ibin)   = 0.0




      aer(iso4_a,jsolid,ibin) = electrolyte(jcaso4,jsolid,ibin)
      aer(ino3_a,jsolid,ibin) = 0.0
      aer(icl_a,jsolid,ibin)  = 0.0
      aer(inh4_a,jsolid,ibin) = 0.0
      aer(ioc_a,jsolid,ibin)  = aer(ioc_a,jtotal,ibin)
      aer(imsa_a,jsolid,ibin) = 0.0
      aer(ico3_a,jsolid,ibin) = aer(ico3_a,jtotal,ibin)
      aer(ina_a,jsolid,ibin)  = 0.0
      aer(ica_a,jsolid,ibin)  = electrolyte(jcaco3,jsolid,ibin) + &
                                electrolyte(jcaso4,jsolid,ibin)
      aer(ibc_a,jsolid,ibin)  = aer(ibc_a,jtotal,ibin)
      aer(ioin_a,jsolid,ibin) = aer(ioin_a,jtotal,ibin)
      aer(ipcg1_b_c_a,jsolid,ibin)= aer(ipcg1_b_c_a,jtotal,ibin)
      aer(ipcg2_b_c_a,jsolid,ibin)= aer(ipcg2_b_c_a,jtotal,ibin)
      aer(ipcg3_b_c_a,jsolid,ibin)= aer(ipcg3_b_c_a,jtotal,ibin)
      aer(ipcg4_b_c_a,jsolid,ibin)= aer(ipcg4_b_c_a,jtotal,ibin)
      aer(ipcg5_b_c_a,jsolid,ibin)= aer(ipcg5_b_c_a,jtotal,ibin)
      aer(ipcg6_b_c_a,jsolid,ibin)= aer(ipcg6_b_c_a,jtotal,ibin)
      aer(ipcg7_b_c_a,jsolid,ibin)= aer(ipcg7_b_c_a,jtotal,ibin)
      aer(ipcg8_b_c_a,jsolid,ibin)= aer(ipcg8_b_c_a,jtotal,ibin)
      aer(ipcg9_b_c_a,jsolid,ibin)= aer(ipcg9_b_c_a,jtotal,ibin)
      aer(ipcg1_b_o_a,jsolid,ibin)= aer(ipcg1_b_o_a,jtotal,ibin)
      aer(ipcg2_b_o_a,jsolid,ibin)= aer(ipcg2_b_o_a,jtotal,ibin)
      aer(ipcg3_b_o_a,jsolid,ibin)= aer(ipcg3_b_o_a,jtotal,ibin)
      aer(ipcg4_b_o_a,jsolid,ibin)= aer(ipcg4_b_o_a,jtotal,ibin)
      aer(ipcg5_b_o_a,jsolid,ibin)= aer(ipcg5_b_o_a,jtotal,ibin)
      aer(ipcg6_b_o_a,jsolid,ibin)= aer(ipcg6_b_o_a,jtotal,ibin)
      aer(ipcg7_b_o_a,jsolid,ibin)= aer(ipcg7_b_o_a,jtotal,ibin)
      aer(ipcg8_b_o_a,jsolid,ibin)= aer(ipcg8_b_o_a,jtotal,ibin)
      aer(ipcg9_b_o_a,jsolid,ibin)= aer(ipcg9_b_o_a,jtotal,ibin)
      aer(iopcg1_b_c_a,jsolid,ibin)= aer(iopcg1_b_c_a,jtotal,ibin)
      aer(iopcg2_b_c_a,jsolid,ibin)= aer(iopcg2_b_c_a,jtotal,ibin)
      aer(iopcg3_b_c_a,jsolid,ibin)= aer(iopcg3_b_c_a,jtotal,ibin)
      aer(iopcg4_b_c_a,jsolid,ibin)= aer(iopcg4_b_c_a,jtotal,ibin)
      aer(iopcg5_b_c_a,jsolid,ibin)= aer(iopcg5_b_c_a,jtotal,ibin)
      aer(iopcg6_b_c_a,jsolid,ibin)= aer(iopcg6_b_c_a,jtotal,ibin)
      aer(iopcg7_b_c_a,jsolid,ibin)= aer(iopcg7_b_c_a,jtotal,ibin)
      aer(iopcg8_b_c_a,jsolid,ibin)= aer(iopcg8_b_c_a,jtotal,ibin)
      aer(iopcg1_b_o_a,jsolid,ibin)= aer(iopcg1_b_o_a,jtotal,ibin)
      aer(iopcg2_b_o_a,jsolid,ibin)= aer(iopcg2_b_o_a,jtotal,ibin)
      aer(iopcg3_b_o_a,jsolid,ibin)= aer(iopcg3_b_o_a,jtotal,ibin)
      aer(iopcg4_b_o_a,jsolid,ibin)= aer(iopcg4_b_o_a,jtotal,ibin)
      aer(iopcg5_b_o_a,jsolid,ibin)= aer(iopcg5_b_o_a,jtotal,ibin)
      aer(iopcg6_b_o_a,jsolid,ibin)= aer(iopcg6_b_o_a,jtotal,ibin)
      aer(iopcg7_b_o_a,jsolid,ibin)= aer(iopcg7_b_o_a,jtotal,ibin)
      aer(iopcg8_b_o_a,jsolid,ibin)= aer(iopcg8_b_o_a,jtotal,ibin)
      aer(ipcg1_f_c_a,jsolid,ibin)= aer(ipcg1_f_c_a,jtotal,ibin)
      aer(ipcg2_f_c_a,jsolid,ibin)= aer(ipcg2_f_c_a,jtotal,ibin)
      aer(ipcg3_f_c_a,jsolid,ibin)= aer(ipcg3_f_c_a,jtotal,ibin)
      aer(ipcg4_f_c_a,jsolid,ibin)= aer(ipcg4_f_c_a,jtotal,ibin)
      aer(ipcg5_f_c_a,jsolid,ibin)= aer(ipcg5_f_c_a,jtotal,ibin)
      aer(ipcg6_f_c_a,jsolid,ibin)= aer(ipcg6_f_c_a,jtotal,ibin)
      aer(ipcg7_f_c_a,jsolid,ibin)= aer(ipcg7_f_c_a,jtotal,ibin)
      aer(ipcg8_f_c_a,jsolid,ibin)= aer(ipcg8_f_c_a,jtotal,ibin)
      aer(ipcg9_f_c_a,jsolid,ibin)= aer(ipcg9_f_c_a,jtotal,ibin)
      aer(ipcg1_f_o_a,jsolid,ibin)= aer(ipcg1_f_o_a,jtotal,ibin)
      aer(ipcg2_f_o_a,jsolid,ibin)= aer(ipcg2_f_o_a,jtotal,ibin)
      aer(ipcg3_f_o_a,jsolid,ibin)= aer(ipcg3_f_o_a,jtotal,ibin)
      aer(ipcg4_f_o_a,jsolid,ibin)= aer(ipcg4_f_o_a,jtotal,ibin)
      aer(ipcg5_f_o_a,jsolid,ibin)= aer(ipcg5_f_o_a,jtotal,ibin)
      aer(ipcg6_f_o_a,jsolid,ibin)= aer(ipcg6_f_o_a,jtotal,ibin)
      aer(ipcg7_f_o_a,jsolid,ibin)= aer(ipcg7_f_o_a,jtotal,ibin)
      aer(ipcg8_f_o_a,jsolid,ibin)= aer(ipcg8_f_o_a,jtotal,ibin)
      aer(ipcg9_f_o_a,jsolid,ibin)= aer(ipcg9_f_o_a,jtotal,ibin)
      aer(iopcg1_f_c_a,jsolid,ibin)= aer(iopcg1_f_c_a,jtotal,ibin)
      aer(iopcg2_f_c_a,jsolid,ibin)= aer(iopcg2_f_c_a,jtotal,ibin)
      aer(iopcg3_f_c_a,jsolid,ibin)= aer(iopcg3_f_c_a,jtotal,ibin)
      aer(iopcg4_f_c_a,jsolid,ibin)= aer(iopcg4_f_c_a,jtotal,ibin)
      aer(iopcg5_f_c_a,jsolid,ibin)= aer(iopcg5_f_c_a,jtotal,ibin)
      aer(iopcg6_f_c_a,jsolid,ibin)= aer(iopcg6_f_c_a,jtotal,ibin)
      aer(iopcg7_f_c_a,jsolid,ibin)= aer(iopcg7_f_c_a,jtotal,ibin)
      aer(iopcg8_f_c_a,jsolid,ibin)= aer(iopcg8_f_c_a,jtotal,ibin)
      aer(iopcg1_f_o_a,jsolid,ibin)= aer(iopcg1_f_o_a,jtotal,ibin)
      aer(iopcg2_f_o_a,jsolid,ibin)= aer(iopcg2_f_o_a,jtotal,ibin)
      aer(iopcg3_f_o_a,jsolid,ibin)= aer(iopcg3_f_o_a,jtotal,ibin)
      aer(iopcg4_f_o_a,jsolid,ibin)= aer(iopcg4_f_o_a,jtotal,ibin)
      aer(iopcg5_f_o_a,jsolid,ibin)= aer(iopcg5_f_o_a,jtotal,ibin)
      aer(iopcg6_f_o_a,jsolid,ibin)= aer(iopcg6_f_o_a,jtotal,ibin)
      aer(iopcg7_f_o_a,jsolid,ibin)= aer(iopcg7_f_o_a,jtotal,ibin)
      aer(iopcg8_f_o_a,jsolid,ibin)= aer(iopcg8_f_o_a,jtotal,ibin)
      aer(ismpa_a,jsolid,ibin)= aer(ismpa_a,jtotal,ibin)
      aer(ismpbb_a,jsolid,ibin)= aer(ismpbb_a,jtotal,ibin)
      aer(iglysoa_r1_a,jsolid,ibin)= aer(iglysoa_r1_a,jtotal,ibin)
      aer(iglysoa_r2_a,jsolid,ibin)= aer(iglysoa_r2_a,jtotal,ibin)
      aer(iglysoa_sfc_a,jsolid,ibin)= aer(iglysoa_sfc_a,jtotal,ibin)
      aer(iglysoa_nh4_a,jsolid,ibin)= aer(iglysoa_nh4_a,jtotal,ibin)
      aer(iglysoa_oh_a,jsolid,ibin)= aer(iglysoa_oh_a,jtotal,ibin)
      aer(iant1_c_a,jsolid,ibin)= aer(iant1_c_a,jtotal,ibin)
      aer(iant2_c_a,jsolid,ibin)= aer(iant2_c_a,jtotal,ibin)
      aer(iant3_c_a,jsolid,ibin)= aer(iant3_c_a,jtotal,ibin)
      aer(iant4_c_a,jsolid,ibin)= aer(iant4_c_a,jtotal,ibin)
      aer(iant1_o_a,jsolid,ibin)= aer(iant1_o_a,jtotal,ibin)
      aer(iant2_o_a,jsolid,ibin)= aer(iant2_o_a,jtotal,ibin)
      aer(iant3_o_a,jsolid,ibin)= aer(iant3_o_a,jtotal,ibin)
      aer(iant4_o_a,jsolid,ibin)= aer(iant4_o_a,jtotal,ibin)
      aer(ibiog1_c_a,jsolid,ibin)= aer(ibiog1_c_a,jtotal,ibin)
      aer(ibiog2_c_a,jsolid,ibin)= aer(ibiog2_c_a,jtotal,ibin)
      aer(ibiog3_c_a,jsolid,ibin)= aer(ibiog3_c_a,jtotal,ibin)
      aer(ibiog4_c_a,jsolid,ibin)= aer(ibiog4_c_a,jtotal,ibin)
      aer(ibiog1_o_a,jsolid,ibin)= aer(ibiog1_o_a,jtotal,ibin)
      aer(ibiog2_o_a,jsolid,ibin)= aer(ibiog2_o_a,jtotal,ibin)
      aer(ibiog3_o_a,jsolid,ibin)= aer(ibiog3_o_a,jtotal,ibin)
      aer(ibiog4_o_a,jsolid,ibin)= aer(ibiog4_o_a,jtotal,ibin)
      aer(iasoaX_a,jsolid,ibin)= aer(iasoaX_a,jtotal,ibin)
      aer(iasoa1_a,jsolid,ibin)= aer(iasoa1_a,jtotal,ibin)
      aer(iasoa2_a,jsolid,ibin)= aer(iasoa2_a,jtotal,ibin)
      aer(iasoa3_a,jsolid,ibin)= aer(iasoa3_a,jtotal,ibin)
      aer(iasoa4_a,jsolid,ibin)= aer(iasoa4_a,jtotal,ibin)
      aer(ibsoaX_a,jsolid,ibin)= aer(ibsoaX_a,jtotal,ibin)
      aer(ibsoa1_a,jsolid,ibin)= aer(ibsoa1_a,jtotal,ibin)
      aer(ibsoa2_a,jsolid,ibin)= aer(ibsoa2_a,jtotal,ibin)
      aer(ibsoa3_a,jsolid,ibin)= aer(ibsoa3_a,jtotal,ibin)
      aer(ibsoa4_a,jsolid,ibin)= aer(ibsoa4_a,jtotal,ibin)





      aer(iso4_a,jliquid,ibin) = aer(iso4_a,jtotal,ibin) - &
                                 aer(iso4_a,jsolid,ibin)
      aer(iso4_a,jliquid,ibin) = max(0.D0, aer(iso4_a,jliquid,ibin))
      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jtotal,ibin)
      aer(icl_a,jliquid,ibin)  = aer(icl_a,jtotal,ibin)
      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jtotal,ibin)
      aer(ioc_a,jliquid,ibin)  = 0.0
      aer(imsa_a,jliquid,ibin) = aer(imsa_a,jtotal,ibin)
      aer(ico3_a,jliquid,ibin) = 0.0
      aer(ina_a,jliquid,ibin)  = aer(ina_a,jtotal,ibin)
      aer(ica_a,jliquid,ibin)  = aer(ica_a,jtotal,ibin) - &
                                 aer(ica_a,jsolid,ibin)
      aer(ica_a,jliquid,ibin)  = max(0.D0, aer(ica_a,jliquid,ibin))
      aer(ibc_a,jliquid,ibin)  = 0.0
      aer(ioin_a,jliquid,ibin) = 0.0
      aer(ipcg1_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg2_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg3_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg4_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg5_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg6_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg7_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg8_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg9_b_c_a,jliquid,ibin)= 0.0
      aer(ipcg1_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg2_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg3_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg4_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg5_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg6_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg7_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg8_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg9_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg1_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg2_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg3_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg4_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg5_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg6_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg7_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg8_b_c_a,jliquid,ibin)= 0.0
      aer(iopcg1_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg2_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg3_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg4_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg5_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg6_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg7_b_o_a,jliquid,ibin)= 0.0
      aer(iopcg8_b_o_a,jliquid,ibin)= 0.0
      aer(ipcg1_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg2_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg3_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg4_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg5_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg6_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg7_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg8_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg9_f_c_a,jliquid,ibin)= 0.0
      aer(ipcg1_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg2_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg3_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg4_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg5_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg6_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg7_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg8_f_o_a,jliquid,ibin)= 0.0
      aer(ipcg9_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg1_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg2_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg3_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg4_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg5_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg6_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg7_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg8_f_c_a,jliquid,ibin)= 0.0
      aer(iopcg1_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg2_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg3_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg4_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg5_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg6_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg7_f_o_a,jliquid,ibin)= 0.0
      aer(iopcg8_f_o_a,jliquid,ibin)= 0.0
      aer(ismpa_a,jliquid,ibin)= 0.0
      aer(ismpbb_a,jliquid,ibin)= 0.0
      aer(iglysoa_r1_a,jliquid,ibin)= 0.0
      aer(iglysoa_r2_a,jliquid,ibin)= 0.0
      aer(iglysoa_sfc_a,jliquid,ibin)= 0.0
      aer(iglysoa_nh4_a,jliquid,ibin)= 0.0
      aer(iglysoa_oh_a,jliquid,ibin)= 0.0
      aer(iant1_c_a,jliquid,ibin)= 0.0
      aer(iant2_c_a,jliquid,ibin)= 0.0
      aer(iant3_c_a,jliquid,ibin)= 0.0
      aer(iant4_c_a,jliquid,ibin)= 0.0
      aer(iant1_o_a,jliquid,ibin)= 0.0
      aer(iant2_o_a,jliquid,ibin)= 0.0
      aer(iant3_o_a,jliquid,ibin)= 0.0
      aer(iant4_o_a,jliquid,ibin)= 0.0
      aer(ibiog1_c_a,jliquid,ibin)= 0.0
      aer(ibiog2_c_a,jliquid,ibin)= 0.0
      aer(ibiog3_c_a,jliquid,ibin)= 0.0
      aer(ibiog4_c_a,jliquid,ibin)= 0.0
      aer(ibiog1_o_a,jliquid,ibin)= 0.0
      aer(ibiog2_o_a,jliquid,ibin)= 0.0
      aer(ibiog3_o_a,jliquid,ibin)= 0.0
      aer(ibiog4_o_a,jliquid,ibin)= 0.0
      aer(iasoaX_a,jliquid,ibin)= 0.0
      aer(iasoa1_a,jliquid,ibin)= 0.0
      aer(iasoa2_a,jliquid,ibin)= 0.0
      aer(iasoa3_a,jliquid,ibin)= 0.0
      aer(iasoa4_a,jliquid,ibin)= 0.0
      aer(ibsoaX_a,jliquid,ibin)= 0.0
      aer(ibsoa1_a,jliquid,ibin)= 0.0
      aer(ibsoa2_a,jliquid,ibin)= 0.0
      aer(ibsoa3_a,jliquid,ibin)= 0.0
      aer(ibsoa4_a,jliquid,ibin)= 0.0





      return
      end subroutine adjust_liquid_aerosol























      subroutine ASTEM(dtchem,vbs_nbin)

      USE module_mosaic_gly, only : glysoa_complex, glysoa_simple





      real(kind=8) dtchem

      integer ibin
      real(kind=8) dumdum
      integer vbs_nbin(1)
      integer start_svoc, Nsoa



      
      integer, save :: iclm_debug, jclm_debug, kclm_debug, ncnt_debug
      data iclm_debug /25/
      data jclm_debug /1/
      data kclm_debug /9/
      data ncnt_debug /2/



      if(iclm_aer .eq. iclm_debug .and.   &
         jclm_aer .eq. jclm_debug .and.   &
         kclm_aer .eq. kclm_debug  .and.   &
         ncorecnt_aer .eq. ncnt_debug)then
        dumdum = 0.0
      endif




      nASTEM_call  = nASTEM_call + 1


      iprint_input = mYES





      do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .ne. no_aerosol)then
          call aerosol_phase_state(ibin)
          if (istat_mosaic_fe1 .lt. 0) return
          call calc_dry_n_wet_aerosol_props(ibin)
        endif

      enddo









      call aerosolmtc(vbs_nbin)
      if (istat_mosaic_fe1 .lt. 0) return


      call ASTEM_non_volatiles(dtchem)	
      if (istat_mosaic_fe1 .lt. 0) return



	  call overall_massbal_in 



      call ASTEM_semi_volatiles(dtchem)	
      if (istat_mosaic_fe1 .lt. 0) return

      if (glysoa_param == glysoa_param_simple)  call glysoa_simple(dtchem)
      if (glysoa_param == glysoa_param_complex) call glysoa_complex(dtchem)


      if (istat_mosaic_fe1 .lt. 0) return

      start_svoc = 1
      Nsoa       = 0
      
      if (vbs_nbin(1).eq.0) then
        start_svoc = ismpa_g
      
      else if (vbs_nbin(1).eq.4) then
        start_svoc = iasoaX_g
      
      else
        start_svoc = ipcg1_b_c_g

      end if
      Nsoa       = ngas_ioa + ngas_soa - start_svoc + 1

      call equilibrium(start_svoc,Nsoa)



















      return
      end subroutine astem









      subroutine print_mosaic_stats( iflag1 )



      integer iflag1

      integer ibin
      real(kind=8) p_mesa_fails, p_astem_fails, dumcnt


      if (iflag1 .le. 0) goto 2000



      dumcnt = float(max(nmesa_call,1))
      p_mesa_fails  = 100.*float(nmesa_fail)/dumcnt
      niter_mesa_avg = float(niter_mesa)/dumcnt

      dumcnt = float(max(nastem_call,1))
      p_astem_fails = 100.*float(nastem_fail)/dumcnt
      nsteps_astem_avg = float(nsteps_astem)/dumcnt


      if (iprint_mosaic_perform_stats .gt. 0) then
        write(6,*)'------------------------------------------------'
        write(6,*)'     astem performance statistics'
        write(6,*)'number of astem calls=', nastem_call
        write(6,*)'percent astem fails  =', nastem_fail
        write(6,*)'avg steps per dtchem =', nsteps_astem_avg
        write(6,*)'max steps per dtchem =', nsteps_astem_max
        write(6,*)'  '
        write(6,*)'     mesa performance statistics'
        write(6,*)'number of mesa calls =', nmesa_call
        write(6,*)'total mesa fails     =', nmesa_fail
        write(6,*)'percent mesa fails   =', p_mesa_fails
        write(6,*)'avg iterations/call  =', niter_mesa_avg
        write(6,*)'max iterations/call  =', niter_mesa_max
        write(6,*)'  '
      endif

      if (iprint_mosaic_fe1 .gt. 0) then
         if ((nfe1_mosaic_cur .gt. 0) .or.   &
             (iprint_mosaic_fe1 .ge. 100)) then
            write(6,*)'-----------------------------------------'
            write(6,*)'mosaic failure count (current step) =',   &
               nfe1_mosaic_cur
            write(6,*)'mosaic failure count (all step tot) =',   &
               nfe1_mosaic_tot
            write(6,*)'  '
         endif
      endif

      if (nfe1_mosaic_tot .gt. 9999) then
         write(6,'(a)') "MOSAIC FAILURE COUNT > 9999 -- SOMETHING IS SERIOUSLY WRONG !!!"
         call peg_error_fatal( lunerr_aer, &
              "---> MOSAIC FAILURE COUNT > 9999 -- SOMETHING IS SERIOUSLY WRONG !!!" )
      endif

2000  continue


      nfe1_mosaic_cur = 0

      nmesa_call   = 0
      nmesa_fail   = 0
      niter_mesa   = 0.0
      niter_mesa_max = 0

      nastem_call = 0
      nastem_fail = 0

      nsteps_astem = 0.0
      nsteps_astem_max = 0.0


      return
      end subroutine print_mosaic_stats






        subroutine  equilibrium(start_ind,N)








        implicit none
        real(kind=8), parameter :: tinys=1.0d-15
        integer, intent(in) :: start_ind, N

        integer, parameter :: itermax=2000
        integer idxfresh(N),idxaged(N)   
        real(kind=8) :: dq,frqfresh(nbin_a),frqaged(nbin_a)
        real(kind=8) :: frqtotfresh,frqtotaged,frt
        real(kind=8) :: xsumfresh(nbin_a),xsumaged(nbin_a)
        real(kind=8) :: mnkfresh,mxkfresh,mnkaged,mxkaged
        real betak

        real(kind=8) ::  Csatfresh(N), Ctotfresh(N)
        real(kind=8) ::  Cgasfresh(N),Caerfresh(N) 

        real(kind=8) ::    Csataged(N), Ctotaged(N)
        real(kind=8) ::  Cgasaged(N),Caeraged(N)
        integer nsolfresh,nsolaged,ntrack,icontfresh,icontaged 
        real(kind=8) :: cpxfresh,cpxaged 
        integer ibin,iter 

        integer iv, jp
         real(kind=8) :: dum, sum_dum, sum_soa, small_oc




        real(kind=8) :: cpx 
        real(kind=8) :: Ctot(N),Caer(N),Cgas(N),Csat(N)
        real(kind=8) :: Paer(ngas_volatile)
        integer :: i

        jp=jtotal
        iter=0
         cpxaged=0.0
        cpxfresh=0.0 
        nsolfresh=0
         nsolaged=0
         icontfresh=0
         icontaged=0
         dq=0.0

          do iv=1,ngas_volatile
           Paer(iv)=0.0
          enddo

          do i=1,N
             flagsoap(i)=1
             Ctot(i) = 0.0
             Ctotaged(i) = 0.0
             Ctotfresh(i) = 0.0
             Caer(i) = 0.0
             Caeraged(i) = 0.0
             Caerfresh(i) = 0.0
             Cgas(i) = 0.0
             Cgasaged(i) = 0.0
             Cgasfresh(i) = 0.0
             Csat(i) = 0.0
             Csataged(i) = 0.0
             Csatfresh(i) = 0.0
          enddo



              do iv = start_ind, (start_ind + N - 1)
        total_species(iv) = gas(iv)
        do ibin = 1, nbin_a
          total_species(iv) = total_species(iv) + aer(iv,jtotal,ibin)
           Paer(iv)=Paer(iv)+aer(iv,jtotal,ibin)
        enddo
      enddo

        do ibin=1,nbin_a
        cpxaged= cpxaged+aer(ioc_a,jp,ibin)
         enddo


        do i=1,N
           Ctot(i)=total_species(start_ind+i-1)
           Caer(i)=Paer(start_ind+i-1)
           Csat(i)=sat_soa(start_ind+i-1)
           Cgas(i)=gas(start_ind+i-1)
         enddo


          do i=1,N
            idxfresh(i)=0
            idxaged(i)=0
          enddo



         do i=1,N
            flagsoap(i)=1 
         enddo






























      do i=1,N
         if (flagsoap(i).eq.2) then 
           icontfresh=icontfresh+1  
            idxfresh(icontfresh) = i  
            Csatfresh(icontfresh)=Csat(i)
            Ctotfresh(icontfresh)=Ctot(i)
            Caerfresh(icontfresh)=Caer(i)
            Cgasfresh(icontfresh)=Cgas(i)
            nsolfresh=nsolfresh+1
         elseif (flagsoap(i).eq.1) then                       
            icontaged=icontaged+1
            idxaged(icontaged) = i
            Csataged(icontaged)=Csat(i)
            Ctotaged(icontaged)=Ctot(i)
            Caeraged(icontaged)=Caer(i)
            Cgasaged(icontaged)=Cgas(i)
            nsolaged=nsolaged+1
         endif
      enddo






      if (nsolfresh.gt.0)  call soap(nsolfresh,Ctotfresh, &
                    Csatfresh,Caerfresh,Cgasfresh,cpxfresh)
      if (nsolaged.gt.0)  call soap(nsolaged,Ctotaged, &
                  Csataged,Caeraged,Cgasaged,cpxaged)



        ntrack=0
       do i=1,N 
         if (idxfresh(i).gt.0) then
         Caer(idxfresh(i))= Caerfresh(i)
         Cgas(idxfresh(i))= Cgasfresh(i)
         Ctot(idxfresh(i))= Ctotfresh(i)
         ntrack=ntrack+1
         endif
         if (idxaged(i).gt.0) then
         Caer(idxaged(i))= Caeraged(i)
         Cgas(idxaged(i))= Cgasaged(i)
         Ctot(idxaged(i))= Ctotaged(i)
         ntrack=ntrack+1
         endif
       enddo

        if (ntrack.ne.N) then
        call wrf_error_fatal3("<stdin>",5373,&
'Error in mapping fresh and primary species arrays')
        endif






         do ibin=1,nbin_a
           xsumfresh(ibin)=0.0
           xsumaged(ibin)=0.0
              xsumaged(ibin)= xsumaged(ibin)+aer(ioc_a,jp,ibin)

         do iv = start_ind, (start_ind + N - 1)
           if (flagsoap(iv-start_ind+1).eq.2) then
               xsumfresh(ibin)= xsumfresh(ibin)+aer(iv,jtotal,ibin)
           elseif (flagsoap(iv-start_ind+1).eq.1) then
              xsumaged(ibin)= xsumaged(ibin)+aer(iv,jtotal,ibin)
                elseif (flagsoap(iv-start_ind+1).eq.0) then
                 print *, 'Error in mapping flagsoap to start_ind'
           endif
         enddo











          if (xsumfresh(ibin).eq.0.0) xsumfresh(ibin)=tinys
          if (xsumaged(ibin).eq.0.0) xsumaged(ibin)=tinys
        enddo





          do iv = start_ind, (start_ind + N - 1)
           if (Ctot(iv-start_ind+1).lt.1d-10) goto 120 
            dq=gas(iv)-Cgas(iv-start_ind+1) 



           frqtotfresh=0.0d0
           frqtotaged=0.0d0
           mnkfresh=0.0d0
           mnkaged=0.0d0
           mxkfresh=0.0d0
           mxkaged=0.0d0
             do ibin=1,nbin_a



           if (flagsoap(iv-start_ind+1).eq.2) then
              frqfresh(ibin)= kg(iv,ibin)*(gas(iv) & 
              -(aer(iv,jtotal,ibin))/xsumfresh(ibin) &
              *Csat(iv-start_ind+1))
          endif

           if (flagsoap(iv-start_ind+1).eq.1) then
              frqaged(ibin)= kg(iv,ibin)*(gas(iv) & 
             -(aer(iv,jtotal,ibin))/xsumaged(ibin) &
              *Csat(iv-start_ind+1))
          endif












            mnkfresh=min(mnkfresh,frqfresh(ibin))
            mnkaged=min(mnkaged,frqaged(ibin))

            mxkfresh=max(mxkfresh,frqfresh(ibin))
            mxkaged=max(mxkaged,frqaged(ibin))
          enddo 

            if (flagsoap(iv-start_ind+1).eq.2) then


          if(dq.gt.0.and.mnkfresh.lt.0.and.mxkfresh.gt.0) then
             do ibin=1,nbin_a
               frqfresh(ibin)=max(frqfresh(ibin)-mnkfresh,0.0d0)
              enddo

          elseif(dq.lt.0.and.mxkfresh.gt.0.and.mnkfresh.lt.0) then
              do ibin=1,nbin_a
              frqfresh(ibin)=min(frqfresh(ibin)-mxkfresh,0.0d0)
              enddo
           endif
           do ibin=1,nbin_a
            frqtotfresh=frqtotfresh+frqfresh(ibin)
           enddo




          do ibin=1,nbin_a
           frqfresh(ibin)=frqfresh(ibin)/frqtotfresh
           enddo
 
            elseif(flagsoap(iv-start_ind+1).eq.1) then

          if(dq.gt.0.and.mnkaged.lt.0.and.mxkaged.gt.0) then
             do ibin=1,nbin_a
               frqaged(ibin)=max(frqaged(ibin)-mnkaged,0.0d0)
              enddo
          elseif(dq.lt.0.and.mxkaged.gt.0.and.mnkaged.lt.0) then
              do ibin=1,nbin_a
              frqaged(ibin)=min(frqaged(ibin)-mxkaged,0.0d0)
              enddo
           endif

           do ibin=1,nbin_a
            frqtotaged=frqtotaged+frqaged(ibin)
           enddo

           do ibin=1,nbin_a
           frqaged(ibin)=frqaged(ibin)/frqtotaged
           enddo

           endif 

           if(dq.gt.0.0d0) then

            
             do ibin=1,nbin_a
                 if (flagsoap(iv-start_ind+1).eq.2) then
                 aer(iv,jtotal,ibin)= aer(iv,jtotal,ibin)+dq*frqfresh(ibin)
                 endif
                if (flagsoap(iv-start_ind+1).eq.1) then
                aer(iv,jtotal,ibin)= aer(iv,jtotal,ibin)+dq*frqaged(ibin)
                endif
             enddo

                gas(iv)=Cgas(iv-start_ind+1)













         elseif(dq.lt.0.0d0) then
            iter=0
100         frt=1.0d0
               do ibin=1,nbin_a
                   if (flagsoap(iv-start_ind+1).eq.2) then


                 if(frqfresh(ibin).gt.0.0d0) &
         frt=MAX(MIN(aer(iv,jtotal,ibin)/abs(-dq*frqfresh(ibin)),frt),0.0d0)

                  elseif(flagsoap(iv-start_ind+1).eq.1) then
               if(frqaged(ibin).gt.0.0d0) &
         frt=MAX(MIN(aer(iv,jtotal,ibin)/abs(-dq*frqaged(ibin)),frt),0.0d0)
                  endif 
               enddo 



         frqtotfresh=0.0d0
         frqtotaged=0.0d0

             do ibin=1,nbin_a
        if (flagsoap(iv-start_ind+1).eq.2) then

               aer(iv,jtotal,ibin)= &

               MAX(aer(iv,jtotal,ibin)+frt*dq*frqfresh(ibin),0.0d0)
         if(aer(iv,jtotal,ibin).lt.tinys) frqfresh(ibin)=0.0d0
              frqtotfresh=frqtotfresh+frqfresh(ibin)

        elseif (flagsoap(iv-start_ind+1).eq.1) then
               aer(iv,jtotal,ibin)= &
               MAX(aer(iv,jtotal,ibin)+frt*dq*frqaged(ibin),0.0d0)
         if(aer(iv,jtotal,ibin).lt.tinys) frqaged(ibin)=0.0d0
              frqtotaged=frqtotaged+frqaged(ibin)
         endif 
             enddo 


          dq=(1.0d0-frt)*dq

         if (flagsoap(iv-start_ind+1).eq.2) then
           if(dq.lt.-1.d-8) then 
             if(frqtotfresh.gt.tinys) then 
              if(iter.le.itermax) then 
                iter = iter + 1
                do ibin = 1,nbin_a
                  frqfresh(ibin) = frqfresh(ibin) / frqtotfresh
                enddo 
             goto 100
            endif 
          endif 
           endif 

          elseif (flagsoap(iv-start_ind+1).eq.1) then
           if(dq.lt.-1.d-8) then
             if(frqtotaged.gt.tinys) then 
              if(iter.le.itermax) then 
                iter = iter + 1
                do ibin = 1,nbin_a
                  frqaged(ibin) = frqaged(ibin) / frqtotaged
                enddo
               goto 100
          endif
            endif
            endif

            
            
            
           endif 



           gas(iv)=Ctot(iv-start_ind+1)
             do ibin=1,nbin_a
               gas(iv)=gas(iv)-aer(iv,jtotal,ibin)
             enddo
        endif 

120       continue
           enddo 

       end subroutine equilibrium














        subroutine spfcn(N,Ctot,Csat,Ca,cpx,tom,fval)


      implicit none
       real(kind=8):: Ctot(N),Csat(N),Ca(N),tom,fval,cpx

         integer i,N
        fval=0.0
         do i=1, N
         Ca(i)=Ctot(i)*tom/(tom+Csat(i)/1)
        fval=fval+Ca(i)/1 
        enddo
          fval=fval+cpx-tom
        return

       end subroutine spfcn


        subroutine soap(N,Ctot,Csat,Ca,Cgas,cpx)





        real(kind=8),  parameter :: xtol = 5.0e-5
          real(kind=8):: Ctot(N),Csat(N),cpx,Ca(N),Cgas(N)
          real(kind=8):: xend,dx,xmid,fend,fmid,sun
         integer i,N,znum
        
         sun=0.0
          do i=1,N
            if (Csat(i).gt.0) then
            sun=sun+Ctot(i)/Csat(i) 
           else
           endif
          enddo
         if(cpx.lt.1e-9.and.sun.le.1.0) then 
           do i=1,N
             Cgas(i)=Ctot(i)
             Ca(i)=0.0
           enddo
         goto 900
        endif

       xend=0.0
       do i=1, N
         xend=xend+Ctot(i)/1 
         enddo
        xend=xend+cpx 
       if (xend.gt.1e-10) then 
           call spfcn(N,Ctot,Csat,Ca,cpx,xend,fend) 
        else

              goto 100
      endif
          if(abs(fend).le.xtol*xend) goto 99 
          if (fend.gt.0.0) then 
         write (2,*) "Error in SOAP"
         goto 50
        endif
           dx=xend-cpx
        do znum=1,200
        dx=0.5*dx
         xmid=xend-dx 
           call spfcn (N,Ctot,Csat,Ca,cpx,xmid,fmid) 
          if(abs(fmid).le.xtol*xmid.or.dx.le.xtol*xmid) goto 100 
           if (fmid.lt.0.0) xend=xmid
         enddo
50       call wrf_message("Error in SOAP")
         call wrf_error_fatal3("<stdin>",5700,&
"Error: max number of iterations reached")


99     xmid=xend
100    continue
        do i=1, N
        Ca(i)=min(Ctot(i), Ca(i))
        Cgas(i)=Ctot(i)-Ca(i)
       enddo
900   continue
        


     return

       end subroutine soap









      subroutine ASTEM_semi_volatiles(dtchem)




      real(kind=8) dtchem

      integer ibin, iv, jp
      real(kind=8) dtmax, t_new, t_old, t_out, xt
      real(kind=8) sum1, sum2, sum3, sum4, sum4a, sum4b, h_flux_s



      t_old = 0.0
      t_out = dtchem


      isteps_ASTEM = 0
      do ibin = 1, nbin_a
        iter_MESA(ibin) = 0
      enddo




10    isteps_ASTEM = isteps_ASTEM + 1


      phi_nh4no3_s = 0.0
      phi_nh4cl_s  = 0.0
      ieqblm_ASTEM = mYES			

      do 501 ibin = 1, nbin_a

        idry_case3a(ibin) = mNO			

        do iv = 1, ngas_ioa
          sfc_a(iv)                  = gas(iv)
          df_gas_s(iv,ibin)          = 0.0
          df_gas_l(iv,ibin)          = 0.0
          flux_s(iv,ibin)            = 0.0
          flux_l(iv,ibin)            = 0.0
          Heff(iv,ibin)              = 0.0
          volatile_s(iv,ibin)        = 0.0
          phi_volatile_s(iv,ibin)    = 0.0
          phi_volatile_l(iv,ibin)    = 0.0
          integrate(iv,jsolid,ibin)  = mNO	
          integrate(iv,jliquid,ibin) = mNO	
        enddo


        if(jaerosolstate(ibin) .eq. all_solid)then
          jphase(ibin) = jsolid
          call ASTEM_flux_dry(ibin)
        elseif(jaerosolstate(ibin) .eq. all_liquid)then
          jphase(ibin) = jliquid
          call ASTEM_flux_wet(ibin)
        elseif(jaerosolstate(ibin) .eq. mixed)then

          if( electrolyte(jnh4no3,jsolid,ibin).gt. 0.0 .or. &
              electrolyte(jnh4cl, jsolid,ibin).gt. 0.0 )then
            call ASTEM_flux_mix(ibin)	
          else
            jphase(ibin) = jliquid
            call ASTEM_flux_wet(ibin)
          endif

        endif

501   continue

      if(ieqblm_ASTEM .eq. mYES)goto 30	





11    call ASTEM_calculate_dtmax(dtchem, dtmax)     
      t_new = t_old + dtmax	
      if(t_new .gt. t_out)then	
        dtmax = t_out - t_old
        t_new = t_out*1.01
      endif





      do 20 iv = 2, 4

        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        sum4 = 0.0
        sum4a= 0.0
        sum4b= 0.0

        do 21 ibin = 1, nbin_a
          if(jaerosolstate(ibin) .eq. no_aerosol)goto 21

          jp = jliquid
          sum1 = sum1 + aer(iv,jp,ibin)/ &
          (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin)*integrate(iv,jp,ibin))

          sum2 = sum2 + kg(iv,ibin)*integrate(iv,jp,ibin)/ &
          (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin)*integrate(iv,jp,ibin))

          jp = jsolid
          sum3 = sum3 + aer(iv,jp,ibin)

          if(flux_s(iv,ibin) .gt. 0.)then
            h_flux_s = dtmax*flux_s(iv,ibin)
            sum4a = sum4a + h_flux_s
            aer(iv,jp,ibin) = aer(iv,jp,ibin) + h_flux_s
          elseif(flux_s(iv,ibin) .lt. 0.)then
            h_flux_s = min(h_s_i_m(iv,ibin),dtmax)*flux_s(iv,ibin)
            sum4b = sum4b + h_flux_s
            aer(iv,jp,ibin) = aer(iv,jp,ibin) + h_flux_s
            aer(iv,jp,ibin) = max(aer(iv,jp,ibin), 0.0D0)
          endif
          
21      continue

        sum4 = sum4a + sum4b



        gas(iv) = (total_species(iv) - (sum1 + sum3 + sum4) )/ &
                              (1. + dtmax*sum2)
        gas(iv) = max(gas(iv), 0.0D0)


        

        do 22 ibin = 1, nbin_a

          if(integrate(iv,jliquid,ibin) .eq. mYES)then
            aer(iv,jliquid,ibin) =  &
             (aer(iv,jliquid,ibin) + dtmax*kg(iv,ibin)*gas(iv))/ &
                  (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin))

          endif

22      continue


20    continue








      do 40 ibin = 1, nbin_a
        if(jaerosolstate(ibin) .eq. no_aerosol)goto 40

        if(jphase(ibin) .eq. jsolid)then
          call form_electrolytes(jsolid,ibin,XT)  
        elseif(jphase(ibin) .eq. jliquid)then
          call form_electrolytes(jliquid,ibin,XT) 
        elseif(jphase(ibin) .eq. jtotal)then
          call form_electrolytes(jsolid,ibin,XT)  
          call form_electrolytes(jliquid,ibin,XT) 
        endif



        do iv = 2, ngas_ioa
          aer(iv,jtotal,ibin)=aer(iv,jsolid,ibin)+aer(iv,jliquid,ibin)
        enddo



        call form_electrolytes(jtotal,ibin,XT)	




        if(jhyst_leg(ibin) .eq. jhyst_lo)then
          call ASTEM_update_phase_eqblm(ibin)
        else
          call do_full_deliquescence(ibin)		
        endif
      

40    continue



      t_old = t_new
    

      if(isteps_astem .ge. nmax_astem)then
        nastem_fail = nastem_fail + 1
        write(6,*)'ASTEM internal steps exceeded', nmax_astem
        if(iprint_input .eq. mYES)then
          write(67,*)'ASTEM internal steps exceeded', nmax_astem
          call print_input
          iprint_input = mNO
        endif
        goto 30
      elseif(t_new .lt. t_out)then
        goto 10
      endif



      if(t_new .lt. 0.9999*t_out) goto 10

30    nsteps_astem = nsteps_astem + isteps_astem		
      nsteps_astem_max = max(nsteps_astem_max, isteps_astem)	
























      return
      end subroutine ASTEM_semi_volatiles
     

















      subroutine ASTEM_calculate_dtmax(dtchem, dtmax)
       use module_data_mosaic_other, only:  lunerr



      real(kind=8) dtchem, dtmax

      integer ibin, iv   
      real(kind=8) alpha, h_gas, h_sub_max,  &
           h_gas_i(ngas_ioa), h_gas_l, h_gas_s,  &
           sum_kg_phi, sumflux_s


      h_sub_max = 100.0	







      h_gas_s = 2.e16

      do 5 iv = 2, ngas_ioa  
        h_gas_i(iv) = 1.e16
        sumflux_s = 0.0
        do ibin = 1, nbin_a
          if(flux_s(iv,ibin) .gt. 0.0)then
            sumflux_s = sumflux_s + flux_s(iv,ibin)
          endif        
        enddo
        
        if(sumflux_s .gt. 0.0)then
          h_gas_i(iv) = 0.1*gas(iv)/sumflux_s     
          h_gas_s     = min(h_gas_s, h_gas_i(iv))
        endif

5     continue
      




      h_gas_l = 2.e16

      do 6 iv = 2, ngas_ioa  
        h_gas_i(iv) = 1.e16
        sum_kg_phi = 0.0
        do ibin = 1, nbin_a
          if(integrate(iv,jliquid,ibin) .eq. mYES)then
          sum_kg_phi = sum_kg_phi +  &
                       abs(phi_volatile_l(iv,ibin))*kg(iv,ibin)
          endif        
        enddo
        
        if(sum_kg_phi .gt. 0.0)then
          h_gas_i(iv) = alpha_astem/sum_kg_phi
          h_gas_l     = min(h_gas_l, h_gas_i(iv))
        endif

6     continue

      h_gas = min(h_gas_s, h_gas_l)
      h_gas = min(h_gas, h_sub_max)







      do ibin = 1, nbin_a

        volatile_s(ino3_a,ibin) = electrolyte(jnh4no3,jsolid,ibin)
        volatile_s(inh4_a,ibin) = electrolyte(jnh4cl,jsolid,ibin) +  &
                                  electrolyte(jnh4no3,jsolid,ibin)

        if(idry_case3a(ibin) .eq. mYES)then
          volatile_s(icl_a,ibin)  = aer(icl_a,jsolid,ibin)
        else
          volatile_s(icl_a,ibin)  = electrolyte(jnh4cl,jsolid,ibin)
        endif

      enddo



      do iv = 2, ngas_ioa

        sum_bin_s(iv) = 0.0
        sum_vdf_s(iv) = 0.0
        sum_vol_s(iv) = 0.0

        do ibin = 1, nbin_a
          if(flux_s(iv,ibin) .lt. 0.)then	
            sum_bin_s(iv) = sum_bin_s(iv) + 1.0
            sum_vdf_s(iv) = sum_vdf_s(iv) +  &
                            volatile_s(iv,ibin)*df_gas_s(iv,ibin)
            sum_vol_s(iv) = sum_vol_s(iv) + volatile_s(iv,ibin)
          endif
        enddo

        if(sum_vol_s(iv) .gt. 0.0)then
          avg_df_gas_s(iv) = sum_vdf_s(iv)/sum_vol_s(iv)
        else
          avg_df_gas_s(iv) = 1.0 
        endif

      enddo





      do 20 ibin = 1, nbin_a
        
        if(jaerosolstate(ibin) .eq. no_aerosol) goto 20        
        
        do 10 iv = 2, ngas_ioa

          if(flux_s(iv,ibin) .lt. 0.)then				

            alpha = abs(avg_df_gas_s(iv))/  &
                   (volatile_s(iv,ibin)*sum_bin_s(iv))
            alpha = min(alpha, 1.0D0)

            if(idry_case3a(ibin) .eq. mYES)alpha = 1.0D0

            h_s_i_m(iv,ibin) =  &
                 -alpha*volatile_s(iv,ibin)/flux_s(iv,ibin)

          endif

10      continue
        

20    continue
      

      dtmax = min(dtchem, h_gas)


      if(dtmax .eq. 0.0)then
        write(6,*)' dtmax = ', dtmax
        write(67,*)' dtmax = ', dtmax
        call print_input
        iprint_input = mNO
        call peg_error_fatal( lunerr, " " )
      endif

      return
      end subroutine astem_calculate_dtmax






















      subroutine ASTEM_update_phase_eqblm(ibin)	



      integer ibin

      integer jdum, js, j_index, je	
      real(kind=8) XT, sum_dum	
      



      sum_dum = 0.0
      do je = 1, nelectrolyte
        sum_dum = sum_dum + electrolyte(je,jtotal,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      do je = 1, nelectrolyte
        epercent(je,jtotal,ibin) = 100.*electrolyte(je,jtotal,ibin)/sum_dum
      enddo



      call calculate_XT(ibin,jtotal,XT)		
      

      if(XT .lt. 1. .and. XT .gt. 0. )goto 10	
      
      jdum = 0
      do js = 1, nsalt
        jsalt_present(js) = 0			
        
        if(epercent(js,jtotal,ibin) .gt. ptol_mol_astem)then
          jsalt_present(js) = 1			
          jdum = jdum + jsalt_index(js)
        endif
      enddo
      
      if(jdum .eq. 0)then
        jaerosolstate(ibin) = all_solid 
        jphase(ibin) = jsolid
        call adjust_solid_aerosol(ibin)      
        return
      endif
      
      if(XT .ge. 2.0 .or. XT .lt. 0.0)then
        j_index = jsulf_poor(jdum)
      else
        j_index = jsulf_rich(jdum)
      endif
      
      MDRH(ibin) = MDRH_T(j_index)
      
      if(aH2O*100. .lt. MDRH(ibin)) then
        jaerosolstate(ibin) = all_solid
        jphase(ibin) = jsolid
        call adjust_solid_aerosol(ibin)
        return
      endif



10    if(jphase(ibin) .eq. jsolid)then
        call do_full_deliquescence(ibin)
        call MESA_PTC(ibin)
      else
        call MESA_PTC(ibin)
      endif



      return
      end subroutine ASTEM_update_phase_eqblm






















      subroutine ASTEM_flux_wet(ibin)
      use module_data_mosaic_other, only:  lunerr



      integer ibin

      integer iv, iadjust, iadjust_intermed
      real(kind=8) xt, g_nh3_hno3, g_nh3_hcl, a_nh4_no3, a_nh4_cl



      call ions_to_electrolytes(jliquid,ibin,XT)  	
      call compute_activities(ibin)

      if(water_a(ibin) .eq. 0.0)then
	write(6,*)'Water is zero in liquid phase'
        call peg_error_fatal( lunerr, "Stopping in ASTEM_flux_wet" )
      endif




      if(electrolyte(jcaco3,jsolid,ibin) .gt. 0.0)then
        call ASTEM_flux_wet_case1(ibin)
        return
      endif




      if(XT.lt.1.9999 .and. XT.ge.0.)then
        call ASTEM_flux_wet_case2(ibin)
        return
      endif



      if( (gas(inh3_g)+aer(inh4_a,jliquid,ibin)) .lt. 1.e-25)goto 10  





      iadjust = mNO		
      iadjust_intermed = mNO	


      g_nh3_hno3 = gas(inh3_g)*gas(ihno3_g)
      a_nh4_no3  = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)

      if(g_nh3_hno3 .gt. 0. .and. a_nh4_no3 .eq. 0.)then
        call absorb_tiny_nh4no3(ibin)
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	
        iadjust_intermed = mNO	
      endif


      g_nh3_hcl = gas(inh3_g)*gas(ihcl_g)
      a_nh4_cl  = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)

      if(g_nh3_hcl .gt. 0. .and. a_nh4_cl .eq. 0.)then
        call absorb_tiny_nh4cl(ibin)
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	
      endif
    
      if(iadjust .eq. mYES)then
        call compute_activities(ibin)			
      endif





      kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
      Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	

      kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
      Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl	

      call ASTEM_flux_wet_case3(ibin)

      return






10    iadjust = mNO		
      iadjust_intermed = mNO	


      if(gas(ihno3_g).gt.0. .and. aer(ino3_a,jliquid,ibin).eq.0. .and. &
         aer(icl_a,jliquid,ibin) .gt. 0.0)then
        call absorb_tiny_hno3(ibin)	
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	
        iadjust_intermed = mNO	
      endif


      if(gas(ihcl_g).gt.0. .and. aer(icl_a,jliquid,ibin).eq.0. .and. &
         aer(ino3_a,jliquid,ibin) .gt. 0.0)then
        call absorb_tiny_hcl(ibin)	
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	
      endif

      if(iadjust .eq. mYES)then
        call compute_activities(ibin)			
      endif
      


      call ASTEM_flux_wet_case4(ibin)


      return
      end subroutine ASTEM_flux_wet





















      subroutine ASTEM_flux_wet_case1(ibin)



      integer ibin

      integer iv
      
      mc(jc_h,ibin) = sqrt(Keq_ll(3))


      if(gas(ihno3_g) .gt. 1.e-5)then
        sfc_a(ihno3_g) = 0.0
        df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
        phi_volatile_s(ihno3_g,ibin) = 1.0
        flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
        integrate(ihno3_g,jsolid,ibin) = mYES
        jphase(ibin) = jsolid
        ieqblm_ASTEM = mNO
      endif

      if(gas(ihcl_g) .gt. 1.e-5)then
        sfc_a(ihcl_g)  = 0.0
        df_gas_s(ihcl_g,ibin) = gas(ihcl_g)
        phi_volatile_s(ihcl_g,ibin) = 1.0
        flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
        integrate(ihcl_g,jsolid,ibin)  = mYES
        jphase(ibin) = jsolid
        ieqblm_ASTEM = mNO
      endif

      return
      end subroutine ASTEM_flux_wet_case1






      subroutine ASTEM_flux_wet_case2(ibin)



      integer ibin

      real(kind=8) dum_hno3, dum_hcl, dum_nh3


      sfc_a(inh3_g)  = kel(inh3_g,ibin)* &
                       gam_ratio(ibin)*mc(jc_nh4,ibin)*Keq_ll(3)/ &
                        (mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))

      sfc_a(ihno3_g) = kel(ihno3_g,ibin)* &
                   mc(jc_h,ibin)*ma(ja_no3,ibin)*gam(jhno3,ibin)**2/ &
                   Keq_gl(3)

      sfc_a(ihcl_g)  = kel(ihcl_g,ibin)* &
                   mc(jc_h,ibin)*ma(ja_cl,ibin)*gam(jhcl,ibin)**2/ &
                   Keq_gl(4)

      dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
      dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
      dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))



      if(dum_hno3 .gt. 0.0)then
        df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
        phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
      else
        phi_volatile_l(ihno3_g,ibin)= 0.0
      endif

      if(dum_hcl .gt. 0.0)then
        df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
        phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
      else
        phi_volatile_l(ihcl_g,ibin) = 0.0
      endif

      if(dum_nh3 .gt. 0.0)then
        df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
        phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
      else
        phi_volatile_l(inh3_g,ibin) = 0.0
      endif


      if(phi_volatile_l(ihno3_g,ibin) .le. rtol_eqb_astem .and. &
         phi_volatile_l(ihcl_g,ibin)  .le. rtol_eqb_astem .and. &
         phi_volatile_l(inh3_g,ibin)  .le. rtol_eqb_astem)then

        return

      endif



      if(dum_hno3 .gt. 0.0)then
        Heff(ihno3_g,ibin)=  &
          kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/ &
                       (water_a(ibin)*Keq_gl(3))
        integrate(ihno3_g,jliquid,ibin)= mYES
        ieqblm_ASTEM = mNO
      endif

      if(dum_hcl .gt. 0.0)then
        Heff(ihcl_g,ibin)=  &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/ &
                       (water_a(ibin)*Keq_gl(4))
        integrate(ihcl_g,jliquid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif

      if(dum_nh3 .gt. 0.0)then
        Heff(inh3_g,ibin) =  &
             kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/ &
             (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
        integrate(inh3_g,jliquid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif


      return
      end subroutine ASTEM_flux_wet_case2











      subroutine ASTEM_flux_wet_case3(ibin)



      integer ibin

      real(kind=8) a, b, c, dum_hno3, dum_hcl, dum_nh3



      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)  &
          + kg(ihno3_g,ibin)*gas(ihno3_g)  &
          + kg(ihcl_g,ibin)*gas(ihcl_g)
      c = -(kg(ihno3_g,ibin)*Keq_nh4no3 + kg(ihcl_g,ibin)*Keq_nh4cl)

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihno3_g) = Keq_nh4no3/max(sfc_a(inh3_g),1.D-20)
      sfc_a(ihcl_g)  = Keq_nh4cl/max(sfc_a(inh3_g),1.D-20)



      if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/ &
        (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
      elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/ &
        (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
      else
        call equilibrate_acids(ibin)	
        mc(jc_h,ibin)  = max(mc(jc_h,ibin), sqrt(Keq_ll(3)))

        sfc_a(inh3_g)  = kel(inh3_g,ibin)* &
                         gam_ratio(ibin)*mc(jc_nh4,ibin)*Keq_ll(3)/ &
                        (mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))

        sfc_a(ihno3_g) = kel(ihno3_g,ibin)* &
                   mc(jc_h,ibin)*ma(ja_no3,ibin)*gam(jhno3,ibin)**2/ &
                   Keq_gl(3)
        sfc_a(ihcl_g)  = kel(ihcl_g,ibin)* &
                   mc(jc_h,ibin)*ma(ja_cl,ibin)*gam(jhcl,ibin)**2/ &
                   Keq_gl(4)
      endif



      dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
      dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
      dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))


      if(dum_hno3 .gt. 0.0)then
        df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
        phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
      else
        phi_volatile_l(ihno3_g,ibin)= 0.0
      endif

      if(dum_hcl .gt. 0.0)then
        df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
        phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
      else
        phi_volatile_l(ihcl_g,ibin) = 0.0
      endif

      if(dum_nh3 .gt. 0.0)then
        df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
        phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
      else
        phi_volatile_l(inh3_g,ibin) = 0.0
      endif



      if(phi_volatile_l(ihno3_g,ibin) .le. rtol_eqb_astem .and. &
         phi_volatile_l(ihcl_g,ibin)  .le. rtol_eqb_astem .and. &
         phi_volatile_l(inh3_g,ibin)  .le. rtol_eqb_astem)then

        return

      endif



      if(dum_hno3 .gt. 0.0)then
        Heff(ihno3_g,ibin)=  &
          kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/ &
                       (water_a(ibin)*Keq_gl(3))
        integrate(ihno3_g,jliquid,ibin)= mYES
        ieqblm_ASTEM = mNO
      endif

      if(dum_hcl .gt. 0.0)then
        Heff(ihcl_g,ibin)=  &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/ &
                       (water_a(ibin)*Keq_gl(4))
        integrate(ihcl_g,jliquid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif

      if(dum_nh3 .gt. 0.0)then
        Heff(inh3_g,ibin) =  &
             kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/ &
             (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
        integrate(inh3_g,jliquid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif



      return
      end subroutine ASTEM_flux_wet_case3












      subroutine ASTEM_flux_wet_case3a(ibin)	



      integer ibin

      real(kind=8) a, b, c, dum_hno3, dum_nh3




      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g) &
          + kg(ihno3_g,ibin)*gas(ihno3_g) 
      c = -(kg(ihno3_g,ibin)*Keq_nh4no3)

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihno3_g) = Keq_nh4no3/sfc_a(inh3_g)



      if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/ &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
      else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
      endif



      dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
      dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))


      if(dum_hno3 .gt. 0.0)then
        df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
        phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
      else
        phi_volatile_l(ihno3_g,ibin)= 0.0
      endif

      if(dum_nh3 .gt. 0.0)then
        df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
        phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
      else
        phi_volatile_l(inh3_g,ibin) = 0.0
      endif


      if(phi_volatile_l(ihno3_g,ibin) .le. rtol_eqb_astem .and. &
         phi_volatile_l(inh3_g,ibin)  .le. rtol_eqb_astem)then

        return

      endif



      Heff(ihno3_g,ibin)=  &
        kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/ &
                     (water_a(ibin)*Keq_gl(3))
      integrate(ihno3_g,jliquid,ibin)= mYES


      Heff(inh3_g,ibin) =  &
           kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/ &
           (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
      integrate(inh3_g,jliquid,ibin) = mYES


      ieqblm_ASTEM = mNO


      return
      end subroutine ASTEM_flux_wet_case3a












      subroutine ASTEM_flux_wet_case3b(ibin)	



      integer ibin

      real(kind=8) a, b, c, dum_hcl, dum_nh3



      
      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g) &
          + kg(ihcl_g,ibin)*gas(ihcl_g)  
      c = -(kg(ihcl_g,ibin)*Keq_nh4cl)
        
      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihcl_g)  = Keq_nh4cl /sfc_a(inh3_g)



      if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/ &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
      else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
      endif



      dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
      dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))



      if(dum_hcl .gt. 0.0)then
        df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
        phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
      else
        phi_volatile_l(ihcl_g,ibin) = 0.0
      endif

      if(dum_nh3 .gt. 0.0)then
        df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
        phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
      else
        phi_volatile_l(inh3_g,ibin) = 0.0
      endif



      if(phi_volatile_l(ihcl_g,ibin)  .le. rtol_eqb_astem .and. &
         phi_volatile_l(inh3_g,ibin)  .le. rtol_eqb_astem)then

        return

      endif




      Heff(ihcl_g,ibin)=  &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/ &
                       (water_a(ibin)*Keq_gl(4))
      integrate(ihcl_g,jliquid,ibin) = mYES


      Heff(inh3_g,ibin) =  &
             kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/ &
             (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
      integrate(inh3_g,jliquid,ibin) = mYES


      ieqblm_ASTEM = mNO



      return
      end subroutine ASTEM_flux_wet_case3b












      subroutine ASTEM_flux_wet_case4(ibin)



      integer ibin

      real(kind=8) dum_numer, dum_denom, gas_eqb_ratio, dum_hno3, dum_hcl
      

      dum_numer = kel(ihno3_g,ibin)*Keq_gl(4)*ma(ja_no3,ibin)* &
                  gam(jhno3,ibin)**2
      dum_denom = kel(ihcl_g,ibin)*Keq_gl(3)*ma(ja_cl ,ibin)* &
                  gam(jhcl,ibin)**2


      if(dum_denom .eq. 0.0 .or. dum_numer .eq. 0.0)then
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
        return
      endif

      gas_eqb_ratio = dum_numer/dum_denom	
     


      sfc_a(ihcl_g) =  &
       ( kg(ihno3_g,ibin)*gas(ihno3_g)+kg(ihcl_g,ibin)*gas(ihcl_g) )/ &
           ( kg(ihcl_g,ibin) + gas_eqb_ratio*kg(ihno3_g,ibin) )
      sfc_a(ihno3_g)= gas_eqb_ratio*sfc_a(ihcl_g)



      if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/ &
        (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
      elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/ &
        (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
      else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
      endif



      dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g)) 
      dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))   


      if(dum_hno3 .gt. 0.0)then
        df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
        phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
      else
        phi_volatile_l(ihno3_g,ibin)= 0.0
      endif

      if(dum_hcl .gt. 0.0)then
        df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
        phi_volatile_l(ihcl_g,ibin)= df_gas_l(ihcl_g,ibin)/dum_hcl
      else
        phi_volatile_l(ihcl_g,ibin)= 0.0
      endif


      if(phi_volatile_l(ihno3_g,ibin) .le. rtol_eqb_astem .and. &
         phi_volatile_l(ihcl_g,ibin)  .le. rtol_eqb_astem)then

        return

      endif




      Heff(ihno3_g,ibin)=  &
          kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/ &
                       (water_a(ibin)*Keq_gl(3))
      integrate(ihno3_g,jliquid,ibin)= mYES


      Heff(ihcl_g,ibin)=  &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/ &
                       (water_a(ibin)*Keq_gl(4))
      integrate(ihcl_g,jliquid,ibin) = mYES


      ieqblm_ASTEM = mNO



      return
      end subroutine ASTEM_flux_wet_case4

























      subroutine ASTEM_flux_dry(ibin)



      integer ibin

      integer iv
      real(kind=8) XT, prod_nh4no3, prod_nh4cl, volatile_cl
     
     
     
      
      call calculate_XT(ibin,jsolid,XT)
      



      if(electrolyte(jcaco3,jsolid,ibin) .gt. 0.0)then
        
        call ASTEM_flux_dry_case1(ibin)
      
        return
      endif




      if(XT.lt.1.9999 .and. XT.ge.0.)then	

	call ASTEM_flux_dry_case2(ibin)
     
        return
      endif




      volatile_cl  = electrolyte(jnacl,jsolid,ibin) + &
                     electrolyte(jcacl2,jsolid,ibin)
      

      if(volatile_cl .gt. 0.0 .and. gas(ihno3_g).gt. 0.0 )then
     
        call ASTEM_flux_dry_case3a(ibin)

        Keq_nh4cl_0  = min(Kp_nh4cl_0,  Keq_sg(2))	

        prod_nh4cl = max( (gas(inh3_g)*gas(ihcl_g)-Keq_nh4cl_0), 0.0d0) +   &
                     electrolyte(jnh4cl, jsolid,ibin)	

        if(prod_nh4cl .gt. 0.0)then
          call ASTEM_flux_dry_case3b(ibin)
        endif

        return
      endif




      Keq_nh4no3_0 = min(Kp_nh4no3_0, Keq_sg(1))	
      Keq_nh4cl_0  = min(Kp_nh4cl_0,  Keq_sg(2))	

      prod_nh4no3 = max( (gas(inh3_g)*gas(ihno3_g)-Keq_nh4no3_0), 0.0d0) +   &
                    electrolyte(jnh4no3,jsolid,ibin)	
      prod_nh4cl  = max( (gas(inh3_g)*gas(ihcl_g) -Keq_nh4cl_0), 0.0d0) +   &
                    electrolyte(jnh4cl, jsolid,ibin)	

      if(prod_nh4no3 .gt. 0.0 .or. prod_nh4cl .gt. 0.0)then
        call ASTEM_flux_dry_case4(ibin)
        return
      endif
      


      return                                  
      end subroutine ASTEM_flux_dry
      























      subroutine ASTEM_flux_dry_case1(ibin)



      integer ibin


      if(gas(ihno3_g) .gt. 1.e-5)then
        sfc_a(ihno3_g) = 0.0
        df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
        phi_volatile_s(ihno3_g,ibin) = 1.0
        flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
        integrate(ihno3_g,jsolid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif

      if(gas(ihcl_g) .gt. 1.e-5)then
        sfc_a(ihcl_g)  = 0.0
        df_gas_s(ihcl_g,ibin) = gas(ihcl_g)
        phi_volatile_s(ihcl_g,ibin) = 1.0
        flux_s(ihcl_g,ibin)  = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
        integrate(ihcl_g,jsolid,ibin)  = mYES
        ieqblm_ASTEM = mNO
      endif


      return
      end subroutine ASTEM_flux_dry_case1






      subroutine ASTEM_flux_dry_case2(ibin) 



      integer ibin
      

      if(gas(inh3_g).gt.1.e-5)then
        sfc_a(inh3_g) = 0.0
        df_gas_s(inh3_g,ibin) = gas(inh3_g)
        phi_volatile_s(inh3_g,ibin)  = 1.0
        flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*gas(inh3_g)
        integrate(inh3_g,jsolid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif
      

      return
      end subroutine ASTEM_flux_dry_case2







      subroutine ASTEM_flux_dry_case3a(ibin)



      integer ibin
      

      if(gas(ihno3_g) .gt. 1.e-5)then
        sfc_a(ihno3_g) = 0.0
        sfc_a(ihcl_g)  = gas(ihcl_g) + aer(icl_a,jsolid,ibin)

        df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
        df_gas_s(ihcl_g,ibin)  = -aer(icl_a,jsolid,ibin)
    
        flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*gas(ihno3_g)
        flux_s(ihcl_g,ibin)  = -flux_s(ihno3_g,ibin)

        phi_volatile_s(ihno3_g,ibin) = 1.0
        phi_volatile_s(ihcl_g,ibin)=df_gas_s(ihcl_g,ibin)/sfc_a(ihcl_g)

        integrate(ihno3_g,jsolid,ibin) = mYES
        integrate(ihcl_g,jsolid,ibin)  = mYES

        idry_case3a(ibin) = mYES
        ieqblm_ASTEM = mNO
      endif

      return
      end subroutine ASTEM_flux_dry_case3a







      subroutine ASTEM_flux_dry_case3b(ibin)	



      integer ibin

      integer iactive_nh4cl, js	
      real(kind=8) a, b, c, sum_dum	





      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jnh4cl,jsolid,ibin) = 100.*electrolyte(jnh4cl,jsolid,ibin)/sum_dum






      iactive_nh4cl  = 1



      phi_nh4cl_s = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/ &
                    max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))





      if( abs(phi_nh4cl_s) .lt. rtol_eqb_ASTEM )then
        iactive_nh4cl = 0
      elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and. &
             epercent(jnh4cl, jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4cl = 0
        if(epercent(jnh4cl, jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4cl(ibin)
        endif
      endif



      if(iactive_nh4cl .eq. 0)return

            



      
      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g) &
          + kg(ihcl_g,ibin)*gas(ihcl_g)  
      c = -(kg(ihcl_g,ibin)*Keq_sg(2))
        
      sfc_a(inh3_g) = quadratic(a,b,c)
      sfc_a(ihcl_g) = Keq_sg(2)/sfc_a(inh3_g)

      df_gas_s(ihcl_g,ibin) = gas(ihcl_g) - sfc_a(ihcl_g)
      df_gas_s(inh3_g,ibin) = gas(inh3_g) - sfc_a(inh3_g)
      
      flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*df_gas_s(inh3_g,ibin)
      flux_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin) + flux_s(inh3_g,ibin)

      phi_volatile_s(inh3_g,ibin) = phi_nh4cl_s

      if(flux_s(ihcl_g,ibin) .gt. 0.0)then
        df_gas_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin)/kg(ihcl_g,ibin)	
        phi_volatile_s(ihcl_g,ibin) = phi_nh4cl_s
      else
        sfc_a(ihcl_g)  = gas(ihcl_g) + aer(icl_a,jsolid,ibin)
        df_gas_s(ihcl_g,ibin) = -aer(icl_a,jsolid,ibin)
        phi_volatile_s(ihcl_g,ibin)=df_gas_s(ihcl_g,ibin)/sfc_a(ihcl_g)  
      endif

      integrate(inh3_g,jsolid,ibin) = mYES
      integrate(ihcl_g,jsolid,ibin) = mYES	
            
      ieqblm_ASTEM = mNO

      return
      end subroutine ASTEM_flux_dry_case3b







      subroutine ASTEM_flux_dry_case4(ibin)	



      integer ibin

      integer iactive_nh4no3, iactive_nh4cl, iactive, js	
      real(kind=8) a, b, c, sum_dum					






      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jnh4no3,jsolid,ibin) = 100.*electrolyte(jnh4no3,jsolid,ibin)/sum_dum
      epercent(jnh4cl, jsolid,ibin) = 100.*electrolyte(jnh4cl, jsolid,ibin)/sum_dum





      iactive_nh4no3 = 1
      iactive_nh4cl  = 2



      phi_nh4no3_s = (gas(inh3_g)*gas(ihno3_g) - Keq_sg(1))/ &
                     max(gas(inh3_g)*gas(ihno3_g),Keq_sg(1))
      phi_nh4cl_s  = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/ &
                     max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))






      if( abs(phi_nh4no3_s) .lt. rtol_eqb_ASTEM )then
        iactive_nh4no3 = 0
      elseif(gas(inh3_g)*gas(ihno3_g) .lt. Keq_sg(1) .and. &
             epercent(jnh4no3,jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4no3 = 0
        if(epercent(jnh4no3,jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4no3(ibin)
        endif
      endif


      if( abs(phi_nh4cl_s) .lt. rtol_eqb_ASTEM )then
        iactive_nh4cl = 0
      elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and. &
             epercent(jnh4cl, jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4cl = 0
        if(epercent(jnh4cl, jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4cl(ibin)
        endif
      endif

              
      iactive = iactive_nh4no3 + iactive_nh4cl


      if(iactive .eq. 0)return


      goto (1,2,3),iactive



1     call ASTEM_flux_dry_case4a(ibin)

      return
      
            


2     call ASTEM_flux_dry_case4b(ibin)
            
      return

      


3     call ASTEM_flux_dry_case4ab(ibin)




      return
      end subroutine ASTEM_flux_dry_case4










      subroutine ASTEM_flux_dry_case4a(ibin) 



      integer ibin

      real(kind=8) a, b, c





      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)  &
          + kg(ihno3_g,ibin)*gas(ihno3_g) 
      c = -(kg(ihno3_g,ibin)*Keq_nh4no3_0)	

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihno3_g) = Keq_nh4no3_0/sfc_a(inh3_g) 

      integrate(ihno3_g,jsolid,ibin) = mYES
      integrate(inh3_g,jsolid,ibin)  = mYES

      df_gas_s(ihno3_g,ibin)=gas(ihno3_g)-sfc_a(ihno3_g)
      df_gas_s(inh3_g,ibin) =gas(inh3_g) -sfc_a(inh3_g)
      
      phi_volatile_s(ihno3_g,ibin)= phi_nh4no3_s
      phi_volatile_s(inh3_g,ibin) = phi_nh4no3_s

      flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
      flux_s(inh3_g,ibin)  = flux_s(ihno3_g,ibin)

      ieqblm_ASTEM = mNO

      return
      end subroutine ASTEM_flux_dry_case4a







      subroutine ASTEM_flux_dry_case4b(ibin) 



      integer ibin

      real(kind=8) a, b, c




      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g) &
          + kg(ihcl_g,ibin)*gas(ihcl_g)  
      c = -(kg(ihcl_g,ibin)*Keq_nh4cl_0)	
        
      sfc_a(inh3_g) = quadratic(a,b,c)
      sfc_a(ihcl_g) = Keq_nh4cl_0 /sfc_a(inh3_g)	

      integrate(ihcl_g,jsolid,ibin) = mYES
      integrate(inh3_g,jsolid,ibin) = mYES

      df_gas_s(ihcl_g,ibin) = gas(ihcl_g)-sfc_a(ihcl_g)
      df_gas_s(inh3_g,ibin) = gas(inh3_g)-sfc_a(inh3_g)

      phi_volatile_s(ihcl_g,ibin) = phi_nh4cl_s
      phi_volatile_s(inh3_g,ibin) = phi_nh4cl_s

      flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
      flux_s(inh3_g,ibin) = flux_s(ihcl_g,ibin)

      ieqblm_ASTEM = mNO

      return
      end subroutine ASTEM_flux_dry_case4b







      subroutine ASTEM_flux_dry_case4ab(ibin)	



      integer ibin

      real(kind=8) a, b, c, &
           flux_nh3_est, flux_nh3_max, ratio_flux



      call ASTEM_flux_dry_case4a(ibin)
      call ASTEM_flux_dry_case4b(ibin)




      flux_nh3_est = flux_s(ihno3_g,ibin)+flux_s(ihcl_g,ibin)
      flux_nh3_max = kg(inh3_g,ibin)*gas(inh3_g)


      if(flux_nh3_est .le. flux_nh3_max)then

        flux_s(inh3_g,ibin) = flux_nh3_est			
        sfc_a(inh3_g)       = gas(inh3_g) -  &			
                              flux_s(inh3_g,ibin)/kg(inh3_g,ibin)
        phi_volatile_s(inh3_g,ibin) = max(abs(phi_nh4no3_s), &
                                          abs(phi_nh4cl_s))

      else			
     
        ratio_flux          = flux_nh3_max/flux_nh3_est
        flux_s(inh3_g,ibin) = flux_nh3_max
        flux_s(ihno3_g,ibin)= flux_s(ihno3_g,ibin)*ratio_flux
        flux_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin) *ratio_flux

        sfc_a(inh3_g) = 0.0
        sfc_a(ihno3_g)= gas(ihno3_g) -  &	
                        flux_s(ihno3_g,ibin)/kg(ihno3_g,ibin)
        sfc_a(ihcl_g) = gas(ihcl_g) -   &	
                        flux_s(ihcl_g,ibin)/kg(ihcl_g,ibin)

        df_gas_s(inh3_g,ibin) =gas(inh3_g) -sfc_a(inh3_g)
        df_gas_s(ihno3_g,ibin)=gas(ihno3_g)-sfc_a(ihno3_g)
        df_gas_s(ihcl_g,ibin) =gas(ihcl_g) -sfc_a(ihcl_g)

        phi_volatile_s(inh3_g,ibin) = max(abs(phi_nh4no3_s), &
                                          abs(phi_nh4cl_s))

      endif

      ieqblm_ASTEM = mNO

      return
      end subroutine ASTEM_flux_dry_case4ab






















      subroutine ASTEM_flux_mix(ibin)
      use module_data_mosaic_other, only:  lunerr



      integer ibin

      integer iv, iadjust, iadjust_intermed, js		
      real(kind=8) XT, g_nh3_hno3, g_nh3_hcl, &
           a_nh4_no3, a_nh4_cl, a_no3, a_cl, &
           prod_nh4no3, prod_nh4cl
      real(kind=8) volatile_cl, sum_dum			
     

      call ions_to_electrolytes(jliquid,ibin,XT)  	
      call compute_activities(ibin)

      if(water_a(ibin) .eq. 0.0)then
	write(6,*)'Water is zero in liquid phase'
        call peg_error_fatal( lunerr, "Stopping in ASTEM_flux_wet" )
      endif
      



      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jcaco3,jsolid,ibin) = 100.*electrolyte(jcaco3,jsolid,ibin)/sum_dum



	Keq_nh4no3_0 = Keq_sg(1)	
	Keq_nh4cl_0  = Keq_sg(2)	




      if(epercent(jcaco3,jsolid,ibin) .gt. 0.0)then
        jphase(ibin) = jliquid
        call ASTEM_flux_wet_case1(ibin)
        return
      endif




      if(XT.lt.1.9999 .and. XT.ge.0.)then	
        jphase(ibin) = jliquid
	call ASTEM_flux_wet_case2(ibin)
        return
      endif




      volatile_cl  = electrolyte(jnacl,jsolid,ibin) +   &
                     electrolyte(jcacl2,jsolid,ibin)


      if(volatile_cl .gt. 0.0 .and. gas(ihno3_g).gt. 0.0 )then

        call ASTEM_flux_dry_case3a(ibin)

        prod_nh4cl = max( (gas(inh3_g)*gas(ihcl_g)-Keq_sg(2)), 0.0d0) +   &
                     electrolyte(jnh4cl, jsolid,ibin)

        if(prod_nh4cl .gt. 0.0)then
          call ASTEM_flux_dry_case3b(ibin)
        endif

        jphase(ibin) = jsolid

        return
      endif




      if( electrolyte(jnh4no3,jsolid,ibin).gt.0. .and. &
          electrolyte(jnh4cl,jsolid,ibin) .gt.0. )then
        jphase(ibin) = jsolid
        call ASTEM_flux_dry_case4(ibin)

        if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/ &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
        elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/ &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

        return

      elseif( electrolyte(jnh4no3,jsolid,ibin).gt.0. )then

        g_nh3_hcl= gas(inh3_g)*gas(ihcl_g)
        a_nh4_cl = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)

        iadjust = mNO		
        if(g_nh3_hcl .gt. 0.0 .and. a_nh4_cl .eq. 0.0)then
          call absorb_tiny_nh4cl(ibin)
          iadjust = mYES
        elseif(g_nh3_hcl .eq. 0.0 .and. a_nh4_cl .gt. 0.0)then
          call degas_tiny_nh4cl(ibin)
          iadjust = mYES
        endif
    
        if(iadjust .eq. mYES)then
          call ions_to_electrolytes(jliquid,ibin,XT)  
          call compute_activities(ibin)			
        endif

        call ASTEM_flux_mix_case4a(ibin)	
        jphase(ibin) = jtotal
        return

      elseif( electrolyte(jnh4cl,jsolid,ibin).gt.0.)then

        g_nh3_hno3= gas(inh3_g)*gas(ihno3_g)
        a_nh4_no3 = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)

        iadjust = mNO		
        if(g_nh3_hno3 .gt. 0.0 .and. a_nh4_no3 .eq. 0.0)then
          call absorb_tiny_nh4no3(ibin)
          iadjust = mYES
        elseif(g_nh3_hno3 .eq. 0.0 .and. a_nh4_no3 .gt. 0.0)then
          call degas_tiny_nh4no3(ibin)
          iadjust = mYES
        endif

        if(iadjust .eq. mYES)then
          call ions_to_electrolytes(jliquid,ibin,XT)  	
          call compute_activities(ibin)			
        endif

        kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
        Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	

        call ASTEM_flux_mix_case4b(ibin)	
        jphase(ibin) = jtotal
        return
      endif




      if( (gas(inh3_g)+aer(inh4_a,jliquid,ibin)) .lt. 1.e-25)goto 10  





      iadjust = mNO		
      iadjust_intermed = mNO	


      g_nh3_hno3 = gas(inh3_g)*gas(ihno3_g)
      a_nh4_no3  = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)

      if(g_nh3_hno3 .gt. 0. .and. a_nh4_no3 .eq. 0.)then
        call absorb_tiny_nh4no3(ibin)
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	
        iadjust_intermed = mNO	
      endif


      g_nh3_hcl = gas(inh3_g)*gas(ihcl_g)
      a_nh4_cl  = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)

      if(g_nh3_hcl .gt. 0. .and. a_nh4_cl .eq. 0.)then
        call absorb_tiny_nh4cl(ibin)
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	
      endif

      if(iadjust .eq. mYES)then
        call compute_activities(ibin)			
      endif





      kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
      Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	

      kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
      Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl	

      call ASTEM_flux_wet_case3(ibin)
      jphase(ibin) = jliquid

      return






10    iadjust = mNO		
      iadjust_intermed = mNO	


      if(gas(ihno3_g).gt.0. .and. aer(ino3_a,jliquid,ibin).eq.0. .and.   &
         aer(icl_a,jliquid,ibin) .gt. 0.0)then
        call absorb_tiny_hno3(ibin)	
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	
        iadjust_intermed = mNO	
      endif


      if(gas(ihcl_g).gt.0. .and. aer(icl_a,jliquid,ibin) .eq. 0. .and.   &
         aer(ino3_a,jliquid,ibin) .gt. 0.0)then
        call absorb_tiny_hcl(ibin)			
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	
      endif

      if(iadjust .eq. mYES)then
        call compute_activities(ibin)			
      endif



      call ASTEM_flux_wet_case4(ibin)
      jphase(ibin) = jliquid

     

      return
      end subroutine ASTEM_flux_mix
      












      subroutine ASTEM_flux_mix_case4a(ibin)	



      integer ibin

      integer iactive_nh4no3, iactive_nh4cl, js	
      real(kind=8) sum_dum				



      iactive_nh4no3 = mYES
      iactive_nh4cl  = mYES



      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jnh4no3,jsolid,ibin) = 100.*electrolyte(jnh4no3,jsolid,ibin)/sum_dum




      phi_nh4no3_s = (gas(inh3_g)*gas(ihno3_g) - Keq_sg(1))/ &
                     max(gas(inh3_g)*gas(ihno3_g),Keq_sg(1))


      kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
      Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl	





      if( abs(phi_nh4no3_s) .le. rtol_eqb_ASTEM )then
        iactive_nh4no3 = mNO
      elseif(gas(inh3_g)*gas(ihno3_g) .lt. Keq_sg(1) .and. &
             epercent(jnh4no3,jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4no3 = mNO
        if(epercent(jnh4no3,jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4no3(ibin)
        endif
      endif


      if( gas(inh3_g)*gas(ihcl_g).eq.0. .or. Keq_nh4cl.eq.0. )then
        iactive_nh4cl = mNO
      endif
              


      if(iactive_nh4no3 .eq. mYES)then

        jphase(ibin) = jsolid
        call ASTEM_flux_dry_case4a(ibin)	

        if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/ &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
        elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/ &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

      endif 


      if(iactive_nh4cl .eq. mYES)then

        jphase(ibin) = jliquid
        call ASTEM_flux_wet_case3b(ibin)	

        if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/ &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

      endif


      if(iactive_nh4cl .eq. mYES .and. iactive_nh4no3 .eq. mYES)then
        jphase(ibin) = jtotal
      endif


            
      return
      end subroutine ASTEM_flux_mix_case4a











      subroutine ASTEM_flux_mix_case4b(ibin)	



      integer ibin

      integer iactive_nh4no3, iactive_nh4cl, js	
	real(kind=8) sum_dum				



      iactive_nh4cl  = mYES
      iactive_nh4no3 = mYES



      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jnh4cl,jsolid,ibin) = 100.*electrolyte(jnh4cl,jsolid,ibin)/sum_dum




      phi_nh4cl_s  = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/ &
                     max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))


      kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
      Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	





      if( abs(phi_nh4cl_s) .le. rtol_eqb_ASTEM )then
        iactive_nh4cl = mNO
      elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and. &
             epercent(jnh4cl,jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4cl = mNO
        if(epercent(jnh4cl,jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4cl(ibin)
        endif
      endif


      if( gas(inh3_g)*gas(ihno3_g).eq.0. .or. Keq_nh4no3.eq.0. )then
        iactive_nh4no3 = mNO
      endif



      if(iactive_nh4cl .eq. mYES)then
      
        jphase(ibin) = jsolid
        call ASTEM_flux_dry_case4b(ibin)	

        if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/ &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
        elseif(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/ &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

      endif


      if(iactive_nh4no3 .eq. mYES)then

        jphase(ibin) = jliquid
        call ASTEM_flux_wet_case3a(ibin)	

        if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/ &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

      endif


      if(iactive_nh4cl .eq. mYES .and. iactive_nh4no3 .eq. mYES)then
        jphase(ibin) = jtotal
      endif

                 

      return
      end subroutine ASTEM_flux_mix_case4b


















      subroutine ASTEM_non_volatiles(dtchem) 



      real(kind=8) dtchem

      integer ibin, iupdate_phase_state
      real(kind=8) decay_h2so4, decay_msa,   &
           delta_h2so4, delta_tmsa, delta_nh3, delta_hno3, delta_hcl, &
           delta_so4(nbin_a), delta_msa(nbin_a), &
           delta_nh4(nbin_a)
		
		
      real(kind=8) :: decay_n2o5,   &
           delta_n2o5, delta_clno2, &
           delta_no3_rct1(nbin_a), delta_no3_rct2(nbin_a)
      real(kind=8) XT
    



      sumkg_h2so4 = 0.0
      sumkg_msa   = 0.0
      sumkg_nh3   = 0.0
      sumkg_hno3  = 0.0
      sumkg_hcl   = 0.0
      do ibin = 1, nbin_a
        sumkg_h2so4 = sumkg_h2so4 + kg(ih2so4_g,ibin)
        sumkg_msa   = sumkg_msa   + kg(imsa_g,ibin)
        sumkg_nh3   = sumkg_nh3   + kg(inh3_g,ibin)
        sumkg_hno3  = sumkg_hno3  + kg(ihno3_g,ibin)
        sumkg_hcl   = sumkg_hcl   + kg(ihcl_g,ibin)
      enddo
		
      sumkg_n2o5  = 0.0
      do ibin = 1, nbin_a
        sumkg_n2o5 = sumkg_n2o5 + kg(in2o5_g,ibin)
      enddo





      if(gas(ih2so4_g) .gt. 1.e-14)then


        decay_h2so4   = exp(-sumkg_h2so4*dtchem)
        delta_h2so4   = gas(ih2so4_g)*(1.0 - decay_h2so4)
        gas(ih2so4_g) = gas(ih2so4_g)*decay_h2so4



        do ibin = 1, nbin_a
          if(jaerosolstate(ibin) .ne. no_aerosol)then
            delta_so4(ibin) = delta_h2so4*kg(ih2so4_g,ibin)/sumkg_h2so4
            aer(iso4_a,jtotal,ibin) = aer(iso4_a,jtotal,ibin) + &
                                      delta_so4(ibin)
          endif
        enddo

      else

        delta_h2so4 = 0.0
        do ibin = 1, nbin_a
            delta_so4(ibin) = 0.0
        enddo

      endif






      if(gas(imsa_g) .gt. 1.e-14)then


        decay_msa   = exp(-sumkg_msa*dtchem)
        delta_tmsa  = gas(imsa_g)*(1.0 - decay_msa)
        gas(imsa_g) = gas(imsa_g)*decay_msa


        do ibin = 1, nbin_a
          if(jaerosolstate(ibin) .ne. no_aerosol)then
            delta_msa(ibin) = delta_tmsa*kg(imsa_g,ibin)/sumkg_msa
            aer(imsa_a,jtotal,ibin) = aer(imsa_a,jtotal,ibin) + &
                                      delta_msa(ibin)
          endif
        enddo

      else

        delta_tmsa = 0.0
        do ibin = 1, nbin_a
            delta_msa(ibin) = 0.0
        enddo

      endif





	if(n2o5_flag .gt. 0) then
		
		
		
		
		
		if(gas(in2o5_g) .gt. 1.e-14 .and. sumkg_n2o5 .gt. 0.0)then

			
			decay_n2o5   = exp(-sumkg_n2o5*dtchem)
			delta_n2o5   = gas(in2o5_g)*(1.0 - decay_n2o5)
			gas(in2o5_g) = gas(in2o5_g)*decay_n2o5


			
			do ibin = 1, nbin_a
				if(jaerosolstate(ibin) .ne. no_aerosol)then
					delta_no3_rct1(ibin) = delta_n2o5*frac_n2o5_h2o(ibin)*kg(in2o5_g,ibin)/sumkg_n2o5
					delta_no3_rct2(ibin) = delta_n2o5*(1.0-frac_n2o5_h2o(ibin))*kg(in2o5_g,ibin)/sumkg_n2o5

					aer(ino3_a,jtotal,ibin) = aer(ino3_a,jtotal,ibin) + &
										  (2.0*delta_no3_rct1(ibin)+delta_no3_rct2(ibin))
					
					
					if(aer(icl_a,jtotal,ibin).ge.delta_no3_rct2(ibin))then
						aer(icl_a,jtotal,ibin)  = aer(icl_a,jtotal,ibin) - &
												  delta_no3_rct2(ibin)
						gas(iclno2_g)           = gas(iclno2_g) + &
												  delta_no3_rct2(ibin)
					else
						aer(ino3_a,jtotal,ibin) = aer(ino3_a,jtotal,ibin) + &
												  (delta_no3_rct2(ibin)-aer(icl_a,jtotal,ibin))
						gas(iclno2_g)           = gas(iclno2_g) + &
												  aer(icl_a,jtotal,ibin)

						
						
						delta_no3_rct1(ibin) = delta_no3_rct1(ibin) + (delta_no3_rct2(ibin)-aer(icl_a,jtotal,ibin))
						delta_no3_rct2(ibin) = aer(icl_a,jtotal,ibin)

						aer(icl_a,jtotal,ibin)  = 0.0
					endif
				endif
			enddo

		else

			delta_n2o5 = 0.0
			do ibin = 1, nbin_a
				delta_no3_rct1(ibin) = 0.0
				delta_no3_rct2(ibin) = 0.0
			enddo

		endif
	else
		delta_n2o5 = 0.0	
		do ibin = 1, nbin_a
			delta_no3_rct1(ibin) = 0.0
			delta_no3_rct2(ibin) = 0.0
		enddo
	endif





      delta_nh3 = gas(inh3_g) *(1.0 - exp(-sumkg_nh3*dtchem))
      delta_hno3= gas(ihno3_g)*(1.0 - exp(-sumkg_hno3*dtchem))
      delta_hcl = gas(ihcl_g) *(1.0 - exp(-sumkg_hcl*dtchem))
      

      do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .ne. no_aerosol)then
          delta_nh3_max(ibin) = delta_nh3*kg(inh3_g,ibin)/sumkg_nh3
          delta_hno3_max(ibin)= delta_hno3*kg(ihno3_g,ibin)/sumkg_hno3
          delta_hcl_max(ibin) = delta_hcl*kg(ihcl_g,ibin)/sumkg_hcl
        endif
      enddo


      if(delta_h2so4 .eq. 0.0 .and. delta_tmsa .eq. 0.0 .and. delta_n2o5 .eq. 0.0)then
        iupdate_phase_state = mNO
        goto 100
      endif



      do ibin = 1, nbin_a

        if(epercent(jnacl,jtotal,ibin)  .eq. 0.0 .and. &
           epercent(jcacl2,jtotal,ibin) .eq. 0.0 .and. &
           epercent(jnano3,jtotal,ibin) .eq. 0.0 .and. &
           epercent(jcano3,jtotal,ibin) .eq. 0.0 .and. &
           epercent(jcaco3,jtotal,ibin) .eq. 0.0 .and. &
           jaerosolstate(ibin) .ne. no_aerosol)then
        
          delta_nh4(ibin)=min( (2.*delta_so4(ibin)+delta_msa(ibin)+2.*delta_no3_rct1(ibin)+delta_no3_rct2(ibin)), &
                                delta_nh3_max(ibin) )
     
          aer(inh4_a,jtotal,ibin) = aer(inh4_a,jtotal,ibin) +        &	
                                    delta_nh4(ibin)

          gas(inh3_g) = gas(inh3_g) - delta_nh4(ibin)		

        else

          delta_nh4(ibin) = 0.0

        endif

      enddo

      iupdate_phase_state = mYES



100   if(iupdate_phase_state .eq. mYES)then
        do ibin = 1, nbin_a
          if(jaerosolstate(ibin) .ne. no_aerosol)then
            call conform_electrolytes(jtotal,ibin,XT)
            call aerosol_phase_state(ibin)
          endif
        enddo
      endif

      return
      end subroutine ASTEM_non_volatiles














      subroutine aerosolmtc(vbs_nbin)

      use module_data_mosaic_asect





      integer nghq,vbs_nbin(1)
      integer start_ind
      parameter (nghq = 2)		
      integer ibin, iq, iv
      real(kind=8) tworootpi, root2, beta
      parameter (tworootpi = 3.5449077, root2 = 1.4142135, beta = 2.0)
      real(kind=8) cdum, dp, dp_avg, fkn, kn, lnsg, lndpgn, lndp, speed,   &
           sumghq
      real(kind=8) xghq(nghq), wghq(nghq)			
      real(kind=8) mw_vol(ngas_volatile+ngas_het), v_molar(ngas_volatile+ngas_het), 		     &  
           freepath(ngas_volatile+ngas_het), accom(ngas_volatile+ngas_het),   &
           dg(ngas_volatile+ngas_het) 				









      mw_vol(ih2so4_g) = 98.0
      mw_vol(ihno3_g)  = 63.0
      mw_vol(ihcl_g)   = 36.5
      mw_vol(inh3_g)   = 17.0
      mw_vol(in2o5_g)  = 108.0
      mw_vol(iclno2_g) = 81.5
      mw_vol(imsa_g)   = 96.0
      mw_vol(ipcg1_b_c_g) =250.0
      mw_vol(ipcg2_b_c_g) =250.0
      mw_vol(ipcg3_b_c_g)=250.0
      mw_vol(ipcg4_b_c_g)=250.0
      mw_vol(ipcg5_b_c_g)=250.0
      mw_vol(ipcg6_b_c_g)=250.0
      mw_vol(ipcg7_b_c_g)=250.0
      mw_vol(ipcg8_b_c_g)=250.0
      mw_vol(ipcg9_b_c_g)=250.0
      mw_vol(iopcg1_b_c_g)=250.0
      mw_vol(iopcg2_b_c_g)=250.0
      mw_vol(iopcg3_b_c_g)=250.0
      mw_vol(iopcg4_b_c_g)=250.0
      mw_vol(iopcg5_b_c_g)=250.0
      mw_vol(iopcg6_b_c_g)=250.0
      mw_vol(iopcg7_b_c_g)=250.0
      mw_vol(iopcg8_b_c_g)=250.0
      mw_vol(ipcg1_b_o_g)=250.0
      mw_vol(ipcg2_b_o_g)=250.0
      mw_vol(ipcg3_b_o_g)=250.0
      mw_vol(ipcg4_b_o_g)=250.0
      mw_vol(ipcg5_b_o_g)=250.0
      mw_vol(ipcg6_b_o_g)=250.0
      mw_vol(ipcg7_b_o_g)=250.0
      mw_vol(ipcg8_b_o_g)=250.0
      mw_vol(ipcg9_b_o_g)=250.0
      mw_vol(iopcg1_b_o_g)=250.0
      mw_vol(iopcg2_b_o_g)=250.0
      mw_vol(iopcg3_b_o_g)=250.0
      mw_vol(iopcg4_b_o_g)=250.0
      mw_vol(iopcg5_b_o_g)=250.0
      mw_vol(iopcg6_b_o_g)=250.0
      mw_vol(iopcg7_b_o_g)=250.0
      mw_vol(iopcg8_b_o_g)=250.0
      mw_vol(ipcg1_f_c_g) =250.0
      mw_vol(ipcg2_f_c_g) =250.0
      mw_vol(ipcg3_f_c_g)=250.0
      mw_vol(ipcg4_f_c_g)=250.0
      mw_vol(ipcg5_f_c_g)=250.0
      mw_vol(ipcg6_f_c_g)=250.0
      mw_vol(ipcg7_f_c_g)=250.0
      mw_vol(ipcg8_f_c_g)=250.0
      mw_vol(ipcg9_f_c_g)=250.0
      mw_vol(iopcg1_f_c_g)=250.0
      mw_vol(iopcg2_f_c_g)=250.0
      mw_vol(iopcg3_f_c_g)=250.0
      mw_vol(iopcg4_f_c_g)=250.0
      mw_vol(iopcg5_f_c_g)=250.0
      mw_vol(iopcg6_f_c_g)=250.0
      mw_vol(iopcg7_f_c_g)=250.0
      mw_vol(iopcg8_f_c_g)=250.0
      mw_vol(ipcg1_f_o_g)=250.0
      mw_vol(ipcg2_f_o_g)=250.0
      mw_vol(ipcg3_f_o_g)=250.0
      mw_vol(ipcg4_f_o_g)=250.0
      mw_vol(ipcg5_f_o_g)=250.0
      mw_vol(ipcg6_f_o_g)=250.0
      mw_vol(ipcg7_f_o_g)=250.0
      mw_vol(ipcg8_f_o_g)=250.0
      mw_vol(ipcg9_f_o_g)=250.0
      mw_vol(iopcg1_f_o_g)=250.0
      mw_vol(iopcg2_f_o_g)=250.0
      mw_vol(iopcg3_f_o_g)=250.0
      mw_vol(iopcg4_f_o_g)=250.0
      mw_vol(iopcg5_f_o_g)=250.0
      mw_vol(iopcg6_f_o_g)=250.0
      mw_vol(iopcg7_f_o_g)=250.0
      mw_vol(iopcg8_f_o_g)=250.0
      mw_vol(ismpa_g)=250.0
      mw_vol(ismpbb_g)=250.0
      mw_vol(igly)=58.0
      mw_vol(iho)=17.0
      mw_vol(iant1_c_g)=250.0
      mw_vol(iant2_c_g)=250.0
      mw_vol(iant3_c_g)=250.0
      mw_vol(iant4_c_g)=250.0
      mw_vol(iant1_o_g)=250.0
      mw_vol(iant2_o_g)=250.0
      mw_vol(iant3_o_g)=250.0
      mw_vol(iant4_o_g)=250.0
      mw_vol(ibiog1_c_g)=250.0
      mw_vol(ibiog2_c_g)=250.0
      mw_vol(ibiog3_c_g)=250.0
      mw_vol(ibiog4_c_g)=250.0
      mw_vol(ibiog1_o_g)=250.0
      mw_vol(ibiog2_o_g)=250.0
      mw_vol(ibiog3_o_g)=250.0
      mw_vol(ibiog4_o_g)=250.0
      mw_vol(iasoaX_g)=250.0
      mw_vol(iasoa1_g)=250.0
      mw_vol(iasoa2_g)=250.0
      mw_vol(iasoa3_g)=250.0
      mw_vol(iasoa4_g)=250.0
      mw_vol(ibsoaX_g)=250.0
      mw_vol(ibsoa1_g)=250.0
      mw_vol(ibsoa2_g)=250.0
      mw_vol(ibsoa3_g)=250.0
      mw_vol(ibsoa4_g)=250.0





      v_molar(ih2so4_g)= 42.88
      v_molar(ihno3_g) = 24.11
      v_molar(ihcl_g)  = 21.48
      v_molar(inh3_g)  = 14.90
      v_molar(imsa_g)  = 58.00
      v_molar(in2o5_g) = 60.40
      v_molar(iclno2_g)= 52.70


      accom(ih2so4_g)  = 0.1
      accom(ihno3_g)   = 0.1
      accom(ihcl_g)    = 0.1
      accom(inh3_g)    = 0.1
      accom(in2o5_g)   = 0.1  
      accom(iclno2_g)  = 0.1  
      accom(imsa_g)    = 0.1
      accom(ipcg1_b_c_g) =0.1
      accom(ipcg2_b_c_g) =0.1
      accom(ipcg3_b_c_g)=0.1
      accom(ipcg4_b_c_g)=0.1
      accom(ipcg5_b_c_g)=0.1
      accom(ipcg6_b_c_g)=0.1
      accom(ipcg7_b_c_g)=0.1
      accom(ipcg8_b_c_g)=0.1
      accom(ipcg9_b_c_g)=0.1
      accom(iopcg1_b_c_g)=0.1
      accom(iopcg2_b_c_g)=0.1
      accom(iopcg3_b_c_g)=0.1
      accom(iopcg4_b_c_g)=0.1
      accom(iopcg5_b_c_g)=0.1
      accom(iopcg6_b_c_g)=0.1
      accom(iopcg7_b_c_g)=0.1
      accom(iopcg8_b_c_g)=0.1
      accom(ipcg1_b_o_g)=0.1
      accom(ipcg2_b_o_g)=0.1
      accom(ipcg3_b_o_g)=0.1
      accom(ipcg4_b_o_g)=0.1
      accom(ipcg5_b_o_g)=0.1
      accom(ipcg6_b_o_g)=0.1
      accom(ipcg7_b_o_g)=0.1
      accom(ipcg8_b_o_g)=0.1
      accom(ipcg9_b_o_g)=0.1
      accom(iopcg1_b_o_g)=0.1
      accom(iopcg2_b_o_g)=0.1
      accom(iopcg3_b_o_g)=0.1
      accom(iopcg4_b_o_g)=0.1
      accom(iopcg5_b_o_g)=0.1
      accom(iopcg6_b_o_g)=0.1
      accom(iopcg7_b_o_g)=0.1
      accom(iopcg8_b_o_g)=0.1
      accom(ipcg1_f_c_g) =0.1
      accom(ipcg2_f_c_g) =0.1
      accom(ipcg3_f_c_g)=0.1
      accom(ipcg4_f_c_g)=0.1
      accom(ipcg5_f_c_g)=0.1
      accom(ipcg6_f_c_g)=0.1
      accom(ipcg7_f_c_g)=0.1
      accom(ipcg8_f_c_g)=0.1
      accom(ipcg9_f_c_g)=0.1
      accom(iopcg1_f_c_g)=0.1
      accom(iopcg2_f_c_g)=0.1
      accom(iopcg3_f_c_g)=0.1
      accom(iopcg4_f_c_g)=0.1
      accom(iopcg5_f_c_g)=0.1
      accom(iopcg6_f_c_g)=0.1
      accom(iopcg7_f_c_g)=0.1
      accom(iopcg8_f_c_g)=0.1
      accom(ipcg1_f_o_g)=0.1
      accom(ipcg2_f_o_g)=0.1
      accom(ipcg3_f_o_g)=0.1
      accom(ipcg4_f_o_g)=0.1
      accom(ipcg5_f_o_g)=0.1
      accom(ipcg6_f_o_g)=0.1
      accom(ipcg7_f_o_g)=0.1
      accom(ipcg8_f_o_g)=0.1
      accom(ipcg9_f_o_g)=0.1
      accom(iopcg1_f_o_g)=0.1
      accom(iopcg2_f_o_g)=0.1
      accom(iopcg3_f_o_g)=0.1
      accom(iopcg4_f_o_g)=0.1
      accom(iopcg5_f_o_g)=0.1
      accom(iopcg6_f_o_g)=0.1
      accom(iopcg7_f_o_g)=0.1
      accom(iopcg8_f_o_g)=0.1
      accom(ismpa_g)=0.1
      accom(ismpbb_g)=0.1
      
      accom(igly)=0.1
      accom(iho)=0.1
      accom(iant1_c_g)=0.1
      accom(iant2_c_g)=0.1
      accom(iant3_c_g)=0.1
      accom(iant4_c_g)=0.1
      accom(iant1_o_g)=0.1
      accom(iant2_o_g)=0.1
      accom(iant3_o_g)=0.1
      accom(iant4_o_g)=0.1
      accom(ibiog1_c_g)=0.1
      accom(ibiog2_c_g)=0.1
      accom(ibiog3_c_g)=0.1
      accom(ibiog4_c_g)=0.1
      accom(ibiog1_o_g)=0.1
      accom(ibiog2_o_g)=0.1
      accom(ibiog3_o_g)=0.1
      accom(ibiog4_o_g)=0.1
      accom(iasoaX_g)=0.1
      accom(iasoa1_g)=0.1
      accom(iasoa2_g)=0.1
      accom(iasoa3_g)=0.1
      accom(iasoa4_g)=0.1
      accom(ibsoaX_g)=0.1
      accom(ibsoa1_g)=0.1
      accom(ibsoa2_g)=0.1
      accom(ibsoa3_g)=0.1
      accom(ibsoa4_g)=0.1





      xghq(1) =  0.70710678
      xghq(2) = -0.70710678
      wghq(1) =  0.88622693
      wghq(2) =  0.88622693





      do iv = 1, ngas_ioa
        speed  = mean_molecular_speed(t_k,mw_vol(iv))	
        dg(iv) = gas_diffusivity(t_k,p_atm,mw_vol(iv),v_molar(iv)) 
        freepath(iv) = 3.*dg(iv)/speed			
      enddo


      start_ind = 1
      if(vbs_nbin(1) .eq. 0) then
        start_ind = ismpa_g
      else if (vbs_nbin(1) .eq. 4) then
        start_ind = iasoaX_g
      else
        start_ind = ipcg1_b_c_g
      end if
      
      
      do iv = start_ind, ngas_ioa + ngas_soa+2
        speed = mean_molecular_speed(t_k,mw_vol(iv))    
        dg(iv) = 0.1                                    
        freepath(iv) = 3.*dg(iv)/speed
      enddo


      do iv = (ngas_volatile+1), (ngas_volatile+ngas_het)
        speed = mean_molecular_speed(t_k,mw_vol(iv))	
		dg(iv) = gas_diffusivity(t_k,p_atm,mw_vol(iv),v_molar(iv)) 
		freepath(iv) = 3.*dg(iv)/speed			
      enddo




      if (msize_framework .eq. mmodal) then


      do 10 ibin = 1, nbin_a

        if(jaerosolstate(ibin) .eq. no_aerosol)goto 10
        call calc_dry_n_wet_aerosol_props(ibin)

        dpgn_a(ibin) = dp_wet_a(ibin)	

        lnsg   = log(sigmag_a(ibin))
        lndpgn = log(dpgn_a(ibin))
        cdum   = tworootpi*num_a(ibin)*   &
                 exp(beta*lndpgn + 0.5*(beta*lnsg)**2)

        do 20 iv = 1, ngas_volatile + ngas_het

		  if(iv.eq.in2o5_g)then	
		  						
		  						
		  	if(n2o5_flag.gt.0)then
		  		accom(iv) = acc_n2o5_bert_thorn(water_a(ibin),&
		  							aer(ino3_a,jtotal,ibin),&
		  							aer(icl_a,jtotal,ibin),&
		  							vol_wet_a(ibin))
		  	else
		  		accom(iv) = 0.0
		  	endif
		  end if

          sumghq = 0.0
          do 30 iq = 1, nghq	
            lndp = lndpgn + beta*lnsg**2 + root2*lnsg*xghq(iq)
            dp = exp(lndp)
            kn = 2.*freepath(iv)/dp
            fkn = fuchs_sutugin(kn,accom(iv))
            sumghq = sumghq + wghq(iq)*dp*fkn/(dp**beta)
30        continue

        kg(iv,ibin) = cdum*dg(iv)*sumghq		
20      continue
		
		if(n2o5_flag.gt.0)then
			
	  		
	  		frac_n2o5_h2o(ibin) = split_n2o5_bert_thorn(water_a(ibin),&
		  							aer(icl_a,jtotal,ibin),&
		  							vol_wet_a(ibin))
		else
			frac_n2o5_h2o(ibin) = 0.0
		endif
			
10    continue

      elseif(msize_framework .eq. msection)then


      do 11 ibin = 1, nbin_a

        if(jaerosolstate(ibin) .eq. no_aerosol)goto 11

        call calc_dry_n_wet_aerosol_props(ibin)

        dp_avg = dp_wet_a(ibin)
        cdum  = 6.283185*dp_avg*num_a(ibin)

        do 21 iv = 1, ngas_volatile+ngas_het
		  if(iv.eq.in2o5_g)then	
		  						
		  						
		  	if(n2o5_flag.gt.0)then
		  		accom(iv) = acc_n2o5_bert_thorn(water_a(ibin),&
		  							aer(ino3_a,jtotal,ibin),&
		  							aer(icl_a,jtotal,ibin),&
		  							vol_wet_a(ibin))
		  	else
		  		accom(iv) = 0.0
		  	end if
		  end if
          kn = 2.*freepath(iv)/dp_avg
          fkn = fuchs_sutugin(kn,accom(iv))
          kg(iv,ibin) = cdum*dg(iv)*fkn		

21      continue
		if(n2o5_flag.gt.0)then
			
	  		
	  		frac_n2o5_h2o(ibin) = split_n2o5_bert_thorn(water_a(ibin),&
		  							aer(icl_a,jtotal,ibin),&
		  							vol_wet_a(ibin))	
		else
			frac_n2o5_h2o(ibin) = 0.0
		end if
		
11    continue

      else

        if (iprint_mosaic_fe1 .gt. 0) then
          write(6,*)'error in the choice of msize_framework'
          write(6,*)'mosaic fatal error in subr. aerosolmtc'
        endif

        istat_mosaic_fe1 = -1900
        return

      endif


      return
      end subroutine aerosolmtc


















      subroutine calc_dry_n_wet_aerosol_props(ibin)

      use module_data_mosaic_asect





      integer ibin

      integer jc, je, iaer, isize, itype
      real(kind=8) aer_H
      complex(kind=8) ri_dum



      mass_dry_a(ibin) = 0.0		
      vol_dry_a(ibin)  = 0.0		
      area_dry_a(ibin) = 0.0		

      if(jaerosolstate(ibin) .ne. no_aerosol)then

        aer_H = (2.*aer(iso4_a,jtotal,ibin) +  &
                    aer(ino3_a,jtotal,ibin) +  &
                    aer(icl_a,jtotal,ibin)  +  &
                    aer(imsa_a,jtotal,ibin) +  &
                 2.*aer(ico3_a,jtotal,ibin))-  &
                (2.*aer(ica_a,jtotal,ibin)  +  &
                    aer(ina_a,jtotal,ibin)  +  &
                    aer(inh4_a,jtotal,ibin))

      do iaer = 1, naer
        mass_dry_a(ibin) = mass_dry_a(ibin) +   &
                           aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)	
        vol_dry_a(ibin) = vol_dry_a(ibin) +   &
        aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  	
      enddo
        mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
        vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

      mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15			
      vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15				


        mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3	
        vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3	


        dens_dry_a(ibin) = mass_dry_a(ibin)/vol_dry_a(ibin) 
        dens_wet_a(ibin) = mass_wet_a(ibin)/vol_wet_a(ibin) 


        dp_dry_a(ibin)=(1.90985*vol_dry_a(ibin)/num_a(ibin))**0.3333333	
        dp_wet_a(ibin)=(1.90985*vol_wet_a(ibin)/num_a(ibin))**0.3333333 


        area_dry_a(ibin)= 3.14159*num_a(ibin)*dp_dry_a(ibin)**2	
        area_wet_a(ibin)= 3.14159*num_a(ibin)*dp_wet_a(ibin)**2	



        do je = 1, nelectrolyte
          comp_a(je)=electrolyte(je,jtotal,ibin)*mw_comp_a(je)*1.e-15	
        enddo
        comp_a(joc)  = aer(ioc_a,jtotal,ibin)*mw_comp_a(je)*1.e-15	
        comp_a(jbc)  = aer(ibc_a,jtotal,ibin)*mw_comp_a(je)*1.e-15	
        comp_a(join) = aer(ioin_a,jtotal,ibin)*mw_comp_a(je)*1.e-15	
         comp_a(jpcg1_b_c)= aer(ipcg1_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg2_b_c)= aer(ipcg2_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg3_b_c)= aer(ipcg3_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg4_b_c)= aer(ipcg4_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg5_b_c)= aer(ipcg5_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg6_b_c)= aer(ipcg6_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg7_b_c)= aer(ipcg7_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg8_b_c)= aer(ipcg8_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg9_b_c)= aer(ipcg9_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg1_b_c)= aer(iopcg1_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg2_b_c)= aer(iopcg2_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg3_b_c)= aer(iopcg3_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg4_b_c)= aer(iopcg4_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg5_b_c)= aer(iopcg5_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg6_b_c)= aer(iopcg6_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg7_b_c)= aer(iopcg7_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg8_b_c)= aer(iopcg8_b_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg1_b_o)= aer(ipcg1_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg2_b_o)= aer(ipcg2_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg3_b_o)= aer(ipcg3_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg4_b_o)= aer(ipcg4_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg5_b_o)= aer(ipcg5_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg6_b_o)= aer(ipcg6_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg7_b_o)= aer(ipcg7_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg8_b_o)= aer(ipcg8_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg9_b_o)= aer(ipcg9_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg1_b_o)= aer(iopcg1_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg2_b_o)= aer(iopcg2_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg3_b_o)= aer(iopcg3_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg4_b_o)= aer(iopcg4_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg5_b_o)= aer(iopcg5_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg6_b_o)= aer(iopcg6_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg7_b_o)= aer(iopcg7_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg8_b_o)= aer(iopcg8_b_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg1_f_c)= aer(ipcg1_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg2_f_c)= aer(ipcg2_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg3_f_c)= aer(ipcg3_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg4_f_c)= aer(ipcg4_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg5_f_c)= aer(ipcg5_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg6_f_c)= aer(ipcg6_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg7_f_c)= aer(ipcg7_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg8_f_c)= aer(ipcg8_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg9_f_c)= aer(ipcg9_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg1_f_c)= aer(iopcg1_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg2_f_c)= aer(iopcg2_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg3_f_c)= aer(iopcg3_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg4_f_c)= aer(iopcg4_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg5_f_c)= aer(iopcg5_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg6_f_c)= aer(iopcg6_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg7_f_c)= aer(iopcg7_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg8_f_c)= aer(iopcg8_f_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg1_f_o)= aer(ipcg1_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg2_f_o)= aer(ipcg2_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg3_f_o)= aer(ipcg3_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg4_f_o)= aer(ipcg4_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg5_f_o)= aer(ipcg5_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg6_f_o)= aer(ipcg6_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg7_f_o)= aer(ipcg7_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg8_f_o)= aer(ipcg8_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jpcg9_f_o)= aer(ipcg9_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg1_f_o)= aer(iopcg1_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg2_f_o)= aer(iopcg2_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg3_f_o)= aer(iopcg3_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg4_f_o)= aer(iopcg4_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg5_f_o)= aer(iopcg5_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg6_f_o)= aer(iopcg6_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg7_f_o)= aer(iopcg7_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jopcg8_f_o)= aer(iopcg8_f_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jsmpa)= aer(ismpa_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jsmpbb)= aer(ismpbb_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jglysoa_r1)= aer(iglysoa_r1_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jglysoa_r2)= aer(iglysoa_r2_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jglysoa_sfc)= aer(iglysoa_sfc_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jglysoa_nh4)= aer(iglysoa_nh4_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jglysoa_oh)= aer(iglysoa_oh_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jant1_c)= aer(iant1_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jant2_c)= aer(iant2_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jant3_c)= aer(iant3_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jant4_c)= aer(iant4_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jant1_o)= aer(iant1_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jant2_o)= aer(iant2_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jant3_o)= aer(iant3_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jant4_o)= aer(iant4_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbiog1_c)= aer(ibiog1_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbiog2_c)= aer(ibiog2_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbiog3_c)= aer(ibiog3_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbiog4_c)= aer(ibiog4_c_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbiog1_o)= aer(ibiog1_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbiog2_o)= aer(ibiog2_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbiog3_o)= aer(ibiog3_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbiog4_o)= aer(ibiog4_o_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jasoaX)= aer(iasoaX_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jasoa1)= aer(iasoa1_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jasoa2)= aer(iasoa2_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jasoa3)= aer(iasoa3_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jasoa4)= aer(iasoa4_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbsoaX)= aer(ibsoaX_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbsoa1)= aer(ibsoa1_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbsoa2)= aer(ibsoa2_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbsoa3)= aer(ibsoa3_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 
         comp_a(jbsoa4)= aer(ibsoa4_a,jtotal,ibin)*mw_comp_a(je)*1.e-15 



        comp_a(jh2o) = water_a(ibin)*1.e-3				

        ri_dum = (0.0,0.0)
        do jc = 1, naercomp
          if (dens_comp_a(jc).gt.0) then
          ri_dum = ri_dum + ref_index_a(jc)*comp_a(jc)/dens_comp_a(jc)
          endif
        enddo

        ri_avg_a(ibin) = ri_dum/vol_wet_a(ibin)

      else	

        dens_dry_a(ibin) = 1.0	 
        dens_wet_a(ibin) = 1.0	 

        call isize_itype_from_ibin( ibin, isize, itype )
        dp_dry_a(ibin) = dcen_sect(isize,itype)	
        dp_wet_a(ibin) = dcen_sect(isize,itype)	

        ri_avg_a(ibin) = (1.5,0.0)
      endif


      return
      end subroutine calc_dry_n_wet_aerosol_props


























      subroutine compute_activities(ibin)



      integer ibin

      integer jp, ja
      real(kind=8) xt, xmol(nelectrolyte), sum_elec, dumK, c_bal, a_c
      real(kind=8) quad, aq, bq, cq, xq, dum




      water_a(ibin) = aerosol_water(jliquid,ibin)	
      if(water_a(ibin) .eq. 0.0)return


      call calculate_xt(ibin,jliquid,xt)

      if(xt.gt.2.0 .or. xt.lt.0.)then




      ma(ja_so4,ibin)  = 1.e-9*aer(iso4_a,jliquid,ibin)/water_a(ibin)
      ma(ja_hso4,ibin) = 0.0
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)
      ma(ja_msa,ibin)  = 1.e-9*aer(imsa_a,jliquid,ibin)/water_a(ibin)


      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)
      a_c              = ( 2.d0*ma(ja_so4,ibin)+  &
                                ma(ja_no3,ibin)+  &
                                ma(ja_cl,ibin) +  &
                                ma(ja_msa,ibin) ) - &
                         ( 2.d0*mc(jc_ca,ibin) +  &
                                mc(jc_nh4,ibin)+  &
                                mc(jc_na,ibin) )
      mc(jc_h,ibin) = 0.5*a_c + sqrt(a_c**2 + 4.*Keq_ll(3))

      if(mc(jc_h,ibin) .eq. 0.0)then
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
      endif


      jp = jliquid
      
      
      sum_elec = 2.*electrolyte(jnh4no3,jp,ibin) +  &
                 2.*electrolyte(jnh4cl,jp,ibin)  +  &
                 3.*electrolyte(jnh4so4,jp,ibin) +  &
                 3.*electrolyte(jna2so4,jp,ibin) +  &
                 2.*electrolyte(jnano3,jp,ibin)  +  &
                 2.*electrolyte(jnacl,jp,ibin)   +  &
                 3.*electrolyte(jcano3,jp,ibin)  +  &
                 3.*electrolyte(jcacl2,jp,ibin)  +  &
                 2.*electrolyte(jhno3,jp,ibin)   +  &
                 2.*electrolyte(jhcl,jp,ibin)

      if(sum_elec .eq. 0.0)then
        do ja = 1, nelectrolyte
          gam(ja,ibin) = 1.0
        enddo
        goto 10
      endif
     
     

      xmol(jnh4no3) = 2.*electrolyte(jnh4no3,jp,ibin)/sum_elec
      xmol(jnh4cl)  = 2.*electrolyte(jnh4cl,jp,ibin) /sum_elec
      xmol(jnh4so4) = 3.*electrolyte(jnh4so4,jp,ibin)/sum_elec
      xmol(jna2so4) = 3.*electrolyte(jna2so4,jp,ibin)/sum_elec
      xmol(jnano3)  = 2.*electrolyte(jnano3,jp,ibin) /sum_elec
      xmol(jnacl)   = 2.*electrolyte(jnacl,jp,ibin)  /sum_elec
      xmol(jcano3)  = 3.*electrolyte(jcano3,jp,ibin) /sum_elec
      xmol(jcacl2)  = 3.*electrolyte(jcacl2,jp,ibin) /sum_elec
      xmol(jhno3)   = 2.*electrolyte(jhno3,jp,ibin)  /sum_elec
      xmol(jhcl)    = 2.*electrolyte(jhcl,jp,ibin)   /sum_elec


      ja = jnh4so4
      if(xmol(ja).gt.0.0)then
      log_gam(ja) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnh4so4,ibin) = mc(jc_nh4,ibin)**2*ma(ja_so4,ibin)* &
                               gam(jnh4so4,ibin)**3
      endif



      jA = jnh4no3
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnh4no3,ibin) = mc(jc_nh4,ibin)*ma(ja_no3,ibin)* &
                               gam(jnh4no3,ibin)**2
      endif


      jA = jnh4cl
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnh4cl,ibin)  = mc(jc_nh4,ibin)*ma(ja_cl,ibin)* &
                               gam(jnh4cl,ibin)**2
      endif
      
     
      jA = jna2so4
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jna2so4,ibin) = mc(jc_na,ibin)**2*ma(ja_so4,ibin)* &
                               gam(jna2so4,ibin)**3
      endif


      jA = jnano3
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnano3,ibin)  = mc(jc_na,ibin)*ma(ja_no3,ibin)* &
                               gam(jnano3,ibin)**2
      endif



      jA = jnacl
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnacl,ibin)   = mc(jc_na,ibin)*ma(ja_cl,ibin)* &
                               gam(jnacl,ibin)**2
      endif










     






      jA = jcano3
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jcano3,ibin)  = mc(jc_ca,ibin)*ma(ja_no3,ibin)**2* &
                               gam(jcano3,ibin)**3
      endif


     
      jA = jcacl2
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jcacl2,ibin)  = mc(jc_ca,ibin)*ma(ja_cl,ibin)**2* &
                               gam(jcacl2,ibin)**3
      endif

     
      jA = jhno3
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jhno3,ibin)   = mc(jc_h,ibin)*ma(ja_no3,ibin)* &
                               gam(jhno3,ibin)**2


      jA = jhcl
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +  &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +  &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +  &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +  &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +  &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +  &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +  &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +  &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jhcl,ibin)    = mc(jc_h,ibin)*ma(ja_cl,ibin)* &
                               gam(jhcl,ibin)**2


10    gam(jlvcite,ibin) = 1.0
     
      gam(jnh4hso4,ibin)= 1.0

      gam(jnh4msa,ibin) = 1.0

      gam(jna3hso4,ibin) = 1.0
     
      gam(jnahso4,ibin) = 1.0

      gam(jnamsa,ibin)  = 1.0

      gam(jcamsa2,ibin) = 1.0  

      activity(jlvcite,ibin) = 0.0

      activity(jnh4hso4,ibin)= 0.0

      activity(jnh4msa,ibin) = mc(jc_nh4,ibin)*ma(ja_msa,ibin)* &
                               gam(jnh4msa,ibin)**2
     
      activity(jna3hso4,ibin)= 0.0

      activity(jnahso4,ibin) = 0.0

      activity(jnamsa,ibin) = mc(jc_na,ibin)*ma(ja_msa,ibin)* &  
                               gam(jnamsa,ibin)**2
      
      activity(jcamsa2,ibin) = mc(jc_ca,ibin) * ma(ja_msa,ibin)**2 * &  
                               gam(jcamsa2,ibin)**3

      gam_ratio(ibin) = gam(jnh4no3,ibin)**2/gam(jhno3,ibin)**2


      else


      jp = jliquid
            
      sum_elec = 3.*electrolyte(jh2so4,jp,ibin)    +  &
                 2.*electrolyte(jnh4hso4,jp,ibin)  +  &
                 5.*electrolyte(jlvcite,jp,ibin)   +  &
                 3.*electrolyte(jnh4so4,jp,ibin)   +  &
                 2.*electrolyte(jnahso4,jp,ibin)   +  &
                 5.*electrolyte(jna3hso4,jp,ibin)  +  &
                 3.*electrolyte(jna2so4,jp,ibin)   +  &
                 2.*electrolyte(jhno3,jp,ibin)     +  &
                 2.*electrolyte(jhcl,jp,ibin)
     

      if(sum_elec .eq. 0.0)then
        do jA = 1, nelectrolyte
          gam(jA,ibin) = 1.0
        enddo
        goto 20
      endif
      

      xmol(jh2so4)  = 3.*electrolyte(jh2so4,jp,ibin)/sum_elec
      xmol(jnh4hso4)= 2.*electrolyte(jnh4hso4,jp,ibin)/sum_elec
      xmol(jlvcite) = 5.*electrolyte(jlvcite,jp,ibin)/sum_elec
      xmol(jnh4so4) = 3.*electrolyte(jnh4so4,jp,ibin)/sum_elec
      xmol(jnahso4) = 2.*electrolyte(jnahso4,jp,ibin)/sum_elec
      xmol(jna3hso4)= 5.*electrolyte(jna3hso4,jp,ibin)/sum_elec
      xmol(jna2so4) = 3.*electrolyte(jna2so4,jp,ibin)/sum_elec
      xmol(jhno3)   = 2.*electrolyte(jhno3,jp,ibin)/sum_elec
      xmol(jhcl)    = 2.*electrolyte(jhcl,jp,ibin)/sum_elec
            
      

      jA = jh2so4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)

      

      jA = jhhso4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      
      

      jA = jnh4hso4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      
      

      jA = jlvcite
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      
      

      jA = jnh4so4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      
      

      jA = jnahso4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      


      jA = jna3hso4










      gam(jA,ibin) = 1.0



      jA = jna2so4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)



      jA = jhno3
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      
      

      jA = jhcl
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +  &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+  &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +  &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +  &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +  &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+  &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +  &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +  &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


20    gam(jnh4no3,ibin) = 1.0
      gam(jnh4cl,ibin)  = 1.0
      gam(jnano3,ibin)  = 1.0
      gam(jnacl,ibin)   = 1.0
      gam(jcano3,ibin)  = 1.0
      gam(jcacl2,ibin)  = 1.0

      gam(jnh4msa,ibin) = 1.0
      gam(jnamsa,ibin)  = 1.0
      gam(jcamsa2,ibin) = 1.0  




      mc(jc_ca,ibin)   = 0.0	
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)


      mSULF            = 1.e-9*aer(iso4_a,jliquid,ibin)/water_a(ibin)
      ma(ja_hso4,ibin) = 0.0
      ma(ja_so4,ibin)  = 0.0
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)
      ma(ja_msa,ibin)  = 1.e-9*aer(imsa_a,jliquid,ibin)/water_a(ibin)

      gam_ratio(ibin)  = gam(jnh4hso4,ibin)**2/gam(jhhso4,ibin)**2
      dumK = Keq_ll(1)*gam(jhhso4,ibin)**2/gam(jh2so4,ibin)**3
      
      c_bal =  mc(jc_nh4,ibin) + mc(jc_na,ibin) + 2.*mc(jc_ca,ibin) & 
         - ma(ja_no3,ibin) - ma(ja_cl,ibin) - mSULF - ma(ja_msa,ibin)
      
      aq = 1.0
      bq = dumK + c_bal
      cq = dumK*(c_bal - mSULF)



        if(bq .ne. 0.0)then
        xq = 4.*(1./bq)*(cq/bq)
        else
        xq = 1.e+6
        endif
                
        if(abs(xq) .lt. 1.e-6)then
          dum = xq*(0.5 + xq*(0.125 + xq*0.0625))
          quad = (-0.5*bq/aq)*dum
          if(quad .lt. 0.)then
            quad = -bq/aq - quad
          endif
        else
          quad = 0.5*(-bq+sqrt(bq*bq - 4.*cq))
        endif      


      mc(jc_h,ibin) = max(quad, 1.D-7)
      ma(ja_so4,ibin) = mSULF*dumK/(mc(jc_h,ibin) + dumK)
      ma(ja_hso4,ibin)= mSULF - ma(ja_so4,ibin)


      activity(jcamsa2,ibin) = mc(jc_ca,ibin) * ma(ja_msa,ibin)**2 * & 
                               gam(jcamsa2,ibin)**3

      activity(jnh4so4,ibin) = mc(jc_nh4,ibin)**2*ma(ja_so4,ibin)* &
                               gam(jnh4so4,ibin)**3
     
      activity(jlvcite,ibin) = mc(jc_nh4,ibin)**3*ma(ja_hso4,ibin)* &
                               ma(ja_so4,ibin) * gam(jlvcite,ibin)**5

      activity(jnh4hso4,ibin)= mc(jc_nh4,ibin)*ma(ja_hso4,ibin)* & 
                               gam(jnh4hso4,ibin)**2

      activity(jnh4msa,ibin) = mc(jc_nh4,ibin)*ma(ja_msa,ibin)* &
                               gam(jnh4msa,ibin)**2
     
      activity(jna2so4,ibin) = mc(jc_na,ibin)**2*ma(ja_so4,ibin)* &
                               gam(jna2so4,ibin)**3

      activity(jnahso4,ibin) = mc(jc_na,ibin)*ma(ja_hso4,ibin)* & 
                               gam(jnahso4,ibin)**2

      activity(jnamsa,ibin)  = mc(jc_na,ibin)*ma(ja_msa,ibin)* &
                               gam(jnamsa,ibin)**2
     



      activity(jna3hso4,ibin)= 0.0
     
      activity(jhno3,ibin)   = mc(jc_h,ibin)*ma(ja_no3,ibin)* &
                               gam(jhno3,ibin)**2
      
      activity(jhcl,ibin)    = mc(jc_h,ibin)*ma(ja_cl,ibin)* &
                               gam(jhcl,ibin)**2

      activity(jmsa,ibin)    = mc(jc_h,ibin)*ma(ja_msa,ibin)* &
                               gam(jmsa,ibin)**2
      


      activity(jnh4no3,ibin) = 0.0
     
      activity(jnh4cl,ibin)  = 0.0

      activity(jnano3,ibin)  = 0.0
      
      activity(jnacl,ibin)   = 0.0
     
      activity(jcano3,ibin)  = 0.0
      
      activity(jcacl2,ibin)  = 0.0


      endif




      return
      end subroutine compute_activities






















      subroutine mtem_compute_log_gamz



      integer ja





      ja = jhno3
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)


      ja = jhcl
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)


      ja = jnh4so4
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)


      ja = jnh4no3
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jnh4cl
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jna2so4
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)


      ja = jnano3
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jnacl
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jcano3
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jcacl2
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnh4no3) = fnlog_gamz(ja,jnh4no3)
      log_gamz(ja,jnh4cl)  = fnlog_gamz(ja,jnh4cl)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jnano3)  = fnlog_gamz(ja,jnano3)
      log_gamz(ja,jnacl)   = fnlog_gamz(ja,jnacl)
      log_gamz(ja,jcano3)  = fnlog_gamz(ja,jcano3)
      log_gamz(ja,jcacl2)  = fnlog_gamz(ja,jcacl2)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)



      ja = jh2so4
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jhhso4
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jnh4hso4
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jlvcite
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jnahso4
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)


      ja = jna3hso4
      log_gamz(ja,jh2so4)  = fnlog_gamz(ja,jh2so4)
      log_gamz(ja,jnh4hso4)= fnlog_gamz(ja,jnh4hso4)
      log_gamz(ja,jlvcite) = fnlog_gamz(ja,jlvcite)
      log_gamz(ja,jnh4so4) = fnlog_gamz(ja,jnh4so4)
      log_gamz(ja,jnahso4) = fnlog_gamz(ja,jnahso4)
      log_gamz(ja,jna3hso4)= fnlog_gamz(ja,jna3hso4)
      log_gamz(ja,jna2so4) = fnlog_gamz(ja,jna2so4)
      log_gamz(ja,jhno3)   = fnlog_gamz(ja,jhno3)
      log_gamz(ja,jhcl)    = fnlog_gamz(ja,jhcl)

      return
      end subroutine mtem_compute_log_gamz


































      subroutine calculate_xt(ibin,jp,xt)



      integer ibin, jp
      real(kind=8) xt


      if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
        xt   = ( aer(inh4_a,jp,ibin) +   &
     &           aer(ina_a,jp,ibin)  +   &
     &        2.*aer(ica_a,jp,ibin) )/   &
     &         (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
      else
        xt   = -1.0
      endif


      return
      end subroutine calculate_xt











      subroutine electrolytes_to_ions(jp,ibin)



      integer jp, ibin


      aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jna2so4,jp,ibin) +   &
                         2.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnh4so4,jp,ibin) +   &
                         2.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jh2so4,jp,ibin)

      aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
                         2.*electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jhno3,jp,ibin)

      aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                            electrolyte(jhcl,jp,ibin)

      aer(imsa_a,jp,ibin) = electrolyte(jnh4msa,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)  +   &
                         2.*electrolyte(jcamsa2,jp,ibin) +   &
                            electrolyte(jmsa,jp,ibin)

      aer(ico3_a,jp,ibin) = electrolyte(jcaco3,jp,ibin)

      aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jcaco3,jp,ibin)  +   &
                            electrolyte(jcamsa2,jp,ibin)

      aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                            electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jna2so4,jp,ibin) +   &
                         3.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)

      aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                         2.*electrolyte(jnh4so4,jp,ibin) +   &
                         3.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jnh4msa,jp,ibin)


      return
      end subroutine electrolytes_to_ions




















      subroutine ions_to_electrolytes(jp,ibin,xt)



      integer ibin, jp
      real(kind=8) xt

      integer iaer, je, jc, ja, icase
      real(kind=8) store(naer), sum_dum, sum_naza, sum_nczc, sum_na_nh4,   &
           f_nh4, f_na, xh, xb, xl, xs, cat_net, rem_nh4, rem_na
      real(kind=8) nc(ncation), na(nanion)




      if(jp .ne. jliquid)then
        if (iprint_mosaic_fe1 .gt. 0) then
          write(6,*)' jp must be jliquid'
          write(6,*)' in ions_to_electrolytes sub'
          write(6,*)' wrong jp = ', jp
          write(6,*)' mosaic fatal error in ions_to_electrolytes'
        endif

        istat_mosaic_fe1 = -2000
        return
      endif


      do iaer = 1, naer
      aer(iaer,jp,ibin) = max(0.0D0, aer(iaer,jp,ibin))
      enddo



      store(ica_a)  = aer(ica_a, jp,ibin)
      store(iso4_a) = aer(iso4_a,jp,ibin)

      call form_caso4(store,jp,ibin)

      if(jp .eq. jliquid)then 
        aer(ica_a,jliquid,ibin) = aer(ica_a,jliquid,ibin) -   &
                                  electrolyte(jcaso4,jliquid,ibin)

        aer(iso4_a,jliquid,ibin)= aer(iso4_a,jliquid,ibin)-   &
                                  electrolyte(jcaso4,jliquid,ibin)

        aer(ica_a,jsolid,ibin)  = aer(ica_a,jsolid,ibin) +   &
                                  electrolyte(jcaso4,jliquid,ibin)

        aer(iso4_a,jsolid,ibin) = aer(iso4_a,jsolid,ibin) +   &
                                  electrolyte(jcaso4,jliquid,ibin)

        electrolyte(jcaso4,jsolid,ibin)=electrolyte(jcaso4,jsolid,ibin) &
                                       +electrolyte(jcaso4,jliquid,ibin)
        electrolyte(jcaso4,jliquid,ibin)= 0.0
      endif



      call calculate_xt(ibin,jp,xt)

      if(xt .ge. 1.9999 .or. xt.lt.0.)then
       icase = 1	
      else
       icase = 2	
      endif



      do je = 1, nelectrolyte
        electrolyte(je,jp,ibin) = 0.0
      enddo




      if(icase.eq.1)then 

        na(ja_hso4)= 0.0
        na(ja_so4) = aer(iso4_a,jp,ibin)
        na(ja_no3) = aer(ino3_a,jp,ibin)
        na(ja_cl)  = aer(icl_a, jp,ibin)
        na(ja_msa) = aer(imsa_a,jp,ibin)

        nc(jc_ca)  = aer(ica_a, jp,ibin)
        nc(jc_na)  = aer(ina_a, jp,ibin)
        nc(jc_nh4) = aer(inh4_a,jp,ibin)

        cat_net =&
                 ( 2.*na(ja_so4)+na(ja_no3)+na(ja_cl)+na(ja_msa) )- &
                 ( 2.*nc(jc_ca) +nc(jc_nh4)+nc(jc_na) )

        if(cat_net .lt. 0.0)then

          nc(jc_h) = 0.0

        else  

          nc(jc_h) = cat_net

        endif



      sum_naza = 0.0
      do ja = 1, nanion
        sum_naza = sum_naza + na(ja)*za(ja)
      enddo

      sum_nczc = 0.0
      do jc = 1, ncation
        sum_nczc = sum_nczc + nc(jc)*zc(jc)
      enddo

      if(sum_naza .eq. 0. .or. sum_nczc .eq. 0.)then
        if (iprint_mosaic_diag1 .gt. 0) then
          write(6,*)'mosaic ions_to_electrolytes'
          write(6,*)'ionic concentrations are zero'
          write(6,*)'sum_naza = ', sum_naza
          write(6,*)'sum_nczc = ', sum_nczc
        endif
        return
      endif

      do ja = 1, nanion
        xeq_a(ja) = na(ja)*za(ja)/sum_naza
      enddo

      do jc = 1, ncation
        xeq_c(jc) = nc(jc)*zc(jc)/sum_nczc
      enddo

      na_ma(ja_so4) = na(ja_so4) *mw_a(ja_so4)
      na_ma(ja_no3) = na(ja_no3) *mw_a(ja_no3)
      na_ma(ja_cl)  = na(ja_cl)  *mw_a(ja_cl)
      na_ma(ja_msa) = na(ja_msa) *mw_a(ja_msa)
      na_ma(ja_hso4)= na(ja_hso4)*mw_a(ja_hso4)

      nc_mc(jc_ca)  = nc(jc_ca) *mw_c(jc_ca)
      nc_mc(jc_na)  = nc(jc_na) *mw_c(jc_na)
      nc_mc(jc_nh4) = nc(jc_nh4)*mw_c(jc_nh4)
      nc_mc(jc_h)   = nc(jc_h)  *mw_c(jc_h)



      if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_so4) .gt. 0.)then
        electrolyte(jna2so4,jp,ibin) = (xeq_c(jc_na) *na_ma(ja_so4) + &
                                        xeq_a(ja_so4)*nc_mc(jc_na))/  &
                                         mw_electrolyte(jna2so4)
      endif

      electrolyte(jnahso4,jp,ibin) = 0.0

      if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
        electrolyte(jnamsa,jp,ibin)  = (xeq_c(jc_na) *na_Ma(ja_msa) + &
                                        xeq_a(ja_msa)*nc_Mc(jc_na))/  &
                                         mw_electrolyte(jnamsa)
      endif

      if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
        electrolyte(jnano3, jp,ibin) = (xeq_c(jc_na) *na_ma(ja_no3) + &
                                        xeq_a(ja_no3)*nc_mc(jc_na))/  &
                                         mw_electrolyte(jnano3)
      endif

      if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
        electrolyte(jnacl,  jp,ibin) = (xeq_c(jc_na) *na_ma(ja_cl) +  &
                                        xeq_a(ja_cl) *nc_mc(jc_na))/  &
                                         mw_electrolyte(jnacl)
      endif

      if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_so4) .gt. 0.)then
        electrolyte(jnh4so4,jp,ibin) = (xeq_c(jc_nh4)*na_ma(ja_so4) + &
                                        xeq_a(ja_so4)*nc_mc(jc_nh4))/ &
                                         mw_electrolyte(jnh4so4)
      endif

      electrolyte(jnh4hso4,jp,ibin)= 0.0

      if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
        electrolyte(jnh4msa,jp,ibin) = (xeq_c(jc_nh4)*na_Ma(ja_msa) + &
                                        xeq_a(ja_msa)*nc_Mc(jc_nh4))/ &
                                         mw_electrolyte(jnh4msa)
      endif

      if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
        electrolyte(jnh4no3,jp,ibin) = (xeq_c(jc_nh4)*na_ma(ja_no3) + &
                                        xeq_a(ja_no3)*nc_mc(jc_nh4))/ &
                                         mw_electrolyte(jnh4no3)
      endif

      if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
        electrolyte(jnh4cl, jp,ibin) = (xeq_c(jc_nh4)*na_ma(ja_cl) +  &
                                        xeq_a(ja_cl) *nc_mc(jc_nh4))/ &
                                         mw_electrolyte(jnh4cl)
      endif

      if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.0)then
        electrolyte(jcano3, jp,ibin) = (xeq_c(jc_ca) *na_ma(ja_no3) + &
                                        xeq_a(ja_no3)*nc_mc(jc_ca))/  &
                                         mw_electrolyte(jcano3)
      endif

      if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
        electrolyte(jcacl2, jp,ibin) = (xeq_c(jc_ca) *na_ma(ja_cl) +  &
                                        xeq_a(ja_cl) *nc_mc(jc_ca))/  &
                                         mw_electrolyte(jcacl2)
      endif

      if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
        electrolyte(jcamsa2,jp,ibin) = (xeq_c(jc_ca) *na_Ma(ja_msa) + &
                                        xeq_a(ja_msa) *nc_Mc(jc_ca))/ &
                                         mw_electrolyte(jcamsa2)
      endif

      electrolyte(jh2so4, jp,ibin) = 0.0

      if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
      electrolyte(jhno3,  jp,ibin) = (xeq_c(jc_h)  *na_ma(ja_no3) +   &
                                      xeq_a(ja_no3)*nc_mc(jc_h))/     &
                                       mw_electrolyte(jhno3)
      endif

      if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
        electrolyte(jhcl,   jp,ibin) = (xeq_c(jc_h) *na_ma(ja_cl) +   &
                                        xeq_a(ja_cl)*nc_mc(jc_h))/    &
                                         mw_electrolyte(jhcl)
      endif

      if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
        electrolyte(jmsa,jp,ibin)    = (xeq_c(jc_h) *na_ma(ja_msa) +  &
                                        xeq_a(ja_msa)*nc_mc(jc_h))/   &
                                         mw_electrolyte(jmsa)
      endif



      elseif(icase.eq.2)then 

        store(imsa_a) = aer(imsa_a,jp,ibin)
        store(ica_a)  = aer(ica_a, jp,ibin)
        
        call form_camsa2(store,jp,ibin)

        sum_na_nh4 = aer(ina_a,jp,ibin) + aer(inh4_a,jp,ibin)

        if(sum_na_nh4 .gt. 0.0)then
          f_nh4 = aer(inh4_a,jp,ibin)/sum_na_nh4
          f_na  = aer(ina_a,jp,ibin)/sum_na_nh4
        else
          f_nh4 = 0.0
          f_na  = 0.0
        endif


        if(sum_na_nh4 .gt. store(imsa_a))then
          electrolyte(jnamsa,jp,ibin)  = f_na *store(imsa_a)
          electrolyte(jnh4msa,jp,ibin) = f_nh4*store(imsa_a)
          rem_na = aer(ina_a,jp,ibin) - electrolyte(jnamsa,jp,ibin)  
          rem_nh4= aer(inh4_a,jp,ibin)- electrolyte(jnh4msa,jp,ibin) 
        else
          electrolyte(jnamsa,jp,ibin)  = aer(ina_a,jp,ibin)
          electrolyte(jnh4msa,jp,ibin) = aer(inh4_a,jp,ibin)
          electrolyte(jmsa,jp,ibin)    = store(imsa_a) - sum_na_nh4
          rem_nh4 = 0.0  
          rem_na  = 0.0  
        endif



        if(aer(iso4_a,jp,ibin).gt.0.0)then
          xt = (rem_nh4 + rem_na)/aer(iso4_a,jp,ibin)
        else
          goto 10
        endif

        if(xt .le. 1.0)then	
          xh = (1.0 - xt)
          xb = xt
          electrolyte(jh2so4,jp,ibin)   = xh*aer(iso4_a,jp,ibin)
          electrolyte(jnh4hso4,jp,ibin) = xb*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jnahso4,jp,ibin)  = xb*f_na *aer(iso4_a,jp,ibin)
        elseif(xt .le. 1.5)then	
          xb = 3.0 - 2.0*xt
          xl = xt - 1.0
          electrolyte(jnh4hso4,jp,ibin) = xb*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jnahso4,jp,ibin)  = xb*f_na *aer(iso4_a,jp,ibin)
          electrolyte(jlvcite,jp,ibin)  = xl*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna3hso4,jp,ibin) = xl*f_na *aer(iso4_a,jp,ibin)
        else			
          xl = 2.0 - xt
          xs = 2.0*xt - 3.0
          electrolyte(jlvcite,jp,ibin)  = xl*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna3hso4,jp,ibin) = xl*f_na *aer(iso4_a,jp,ibin)
          electrolyte(jnh4so4,jp,ibin)  = xs*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna2so4,jp,ibin)  = xs*f_na *aer(iso4_a,jp,ibin)
        endif

        electrolyte(jhno3,jp,ibin) = aer(ino3_a,jp,ibin)
        electrolyte(jhcl,jp,ibin)  = aer(icl_a,jp,ibin)

      endif















10    return
      end subroutine ions_to_electrolytes

































      subroutine conform_electrolytes(jp,ibin,xt)



      integer ibin, jp
      real(kind=8) xt

      integer i, ixt_case, je
      real(kind=8) sum_dum, xna_prime, xnh4_prime, xt_prime
      real(kind=8) store(naer)


      do i=1,naer
      aer(i,jp,ibin) = max(0.0D0, aer(i,jp,ibin))
      enddo


      call calculate_xt(ibin,jp,xt)

      if(xt .ge. 1.9999 .or. xt.lt.0.)then
       ixt_case = 1	
      else
       ixt_case = 2	
      endif




      store(iso4_a) = aer(iso4_a,jp,ibin)
      store(ino3_a) = aer(ino3_a,jp,ibin)
      store(icl_a)  = aer(icl_a, jp,ibin)
      store(imsa_a) = aer(imsa_a,jp,ibin)
      store(ico3_a) = aer(ico3_a,jp,ibin)
      store(inh4_a) = aer(inh4_a,jp,ibin)
      store(ina_a)  = aer(ina_a, jp,ibin)
      store(ica_a)  = aer(ica_a, jp,ibin)

      do je=1,nelectrolyte
      electrolyte(je,jp,ibin) = 0.0
      enddo



      if(ixt_case.eq.1)then



        call form_caso4(store,jp,ibin)
        call form_camsa2(store,jp,ibin)
        call form_na2so4(store,jp,ibin)
        call form_namsa(store,jp,ibin)
        call form_cano3(store,jp,ibin)
        call form_nano3(store,jp,ibin)
        call form_nacl(store,jp,ibin)
        call form_cacl2(store,jp,ibin)
        call form_caco3(store,jp,ibin)
        call form_nh4so4(store,jp,ibin)
        call form_nh4msa(store,jp,ibin)
        call form_nh4no3(store,jp,ibin)
        call form_nh4cl(store,jp,ibin)
        call form_msa(store,jp,ibin)
        call degas_hno3(store,jp,ibin)
        call degas_hcl(store,jp,ibin)
        call degas_nh3(store,jp,ibin)

      elseif(ixt_case.eq.2)then



        call form_caso4(store,jp,ibin)
        call form_camsa2(store,jp,ibin)
        call form_namsa(store,jp,ibin)
        call form_nh4msa(store,jp,ibin)
        call form_msa(store,jp,ibin)

        if(store(iso4_a).eq.0.0)goto 10


        xt_prime =(store(ina_a)+store(inh4_a))/   &
                        store(iso4_a)
        xna_prime=0.5*store(ina_a)/store(iso4_a) + 1.

        if(xt_prime.ge.xna_prime)then
          call form_na2so4(store,jp,ibin)
          xnh4_prime = 0.0
          if(store(iso4_a).gt.1.e-15)then
            xnh4_prime = store(inh4_a)/store(iso4_a)
          endif

          if(xnh4_prime .ge. 1.5)then
            call form_nh4so4_lvcite(store,jp,ibin)
          else
            call form_lvcite_nh4hso4(store,jp,ibin)
          endif

        elseif(xt_prime.ge.1.)then
          call form_nh4hso4(store,jp,ibin)
          call form_na2so4_nahso4(store,jp,ibin)
        elseif(xt_prime.lt.1.)then
          call form_nahso4(store,jp,ibin)
          call form_nh4hso4(store,jp,ibin)
          call form_h2so4(store,jp,ibin)
        endif

10    call degas_hno3(store,jp,ibin)
      call degas_hcl(store,jp,ibin)
      call degas_nh3(store,jp,ibin)

      endif 



      call electrolytes_to_ions(jp, ibin)
















      return
      end subroutine conform_electrolytes

















      subroutine form_electrolytes(jp,ibin,xt)



      integer ibin, jp
      real(kind=8) xt

      integer i, ixt_case, j, je
      real(kind=8) sum_dum, xna_prime, xnh4_prime, xt_prime
      real(kind=8) store(naer)


      do i=1,naer
      aer(i,jp,ibin) = max(0.0D0, aer(i,jp,ibin))
      enddo


      call calculate_xt(ibin,jp,xt)

      if(xt .ge. 1.9999 .or. xt.lt.0.)then
       ixt_case = 1	
      else
       ixt_case = 2	
      endif




      store(iso4_a) = aer(iso4_a,jp,ibin)
      store(ino3_a) = aer(ino3_a,jp,ibin)
      store(icl_a)  = aer(icl_a, jp,ibin)
      store(imsa_a) = aer(imsa_a,jp,ibin)
      store(ico3_a) = aer(ico3_a,jp,ibin)
      store(inh4_a) = aer(inh4_a,jp,ibin)
      store(ina_a)  = aer(ina_a, jp,ibin)
      store(ica_a)  = aer(ica_a, jp,ibin)

      do j=1,nelectrolyte
      electrolyte(j,jp,ibin) = 0.0
      enddo



      if(ixt_case.eq.1)then


        call form_caso4(store,jp,ibin)
        call form_camsa2(store,jp,ibin)
        call form_na2so4(store,jp,ibin)
        call form_namsa(store,jp,ibin)
        call form_cano3(store,jp,ibin)
        call form_nano3(store,jp,ibin)
        call form_nacl(store,jp,ibin)
        call form_cacl2(store,jp,ibin)
        call form_caco3(store,jp,ibin)
        call form_nh4so4(store,jp,ibin)
        call form_nh4msa(store,jp,ibin)
        call form_nh4no3(store,jp,ibin)
        call form_nh4cl(store,jp,ibin)
        call form_msa(store,jp,ibin)

        if(jp .eq. jsolid)then
          call degas_hno3(store,jp,ibin)
          call degas_hcl(store,jp,ibin)
          call degas_nh3(store,jp,ibin)
        else
          call form_hno3(store,jp,ibin)
          call form_hcl(store,jp,ibin)
          call degas_nh3(store,jp,ibin)
        endif



      elseif(ixt_case.eq.2)then



        call form_caso4(store,jp,ibin)
        call form_camsa2(store,jp,ibin)
        call form_namsa(store,jp,ibin)
        call form_nh4msa(store,jp,ibin)
        call form_msa(store,jp,ibin)

        if(store(iso4_a).eq.0.0)goto 10


        xt_prime =(store(ina_a)+store(inh4_a))/   &
                        store(iso4_a)
        xna_prime=0.5*store(ina_a)/store(iso4_a) + 1.

        if(xt_prime.ge.xna_prime)then
          call form_na2so4(store,jp,ibin)
          xnh4_prime = 0.0
          if(store(iso4_a).gt.1.e-15)then
            xnh4_prime = store(inh4_a)/store(iso4_a)
          endif

          if(xnh4_prime .ge. 1.5)then
            call form_nh4so4_lvcite(store,jp,ibin)
          else
            call form_lvcite_nh4hso4(store,jp,ibin)
          endif

        elseif(xt_prime.ge.1.)then
          call form_nh4hso4(store,jp,ibin)
          call form_na2so4_nahso4(store,jp,ibin)
        elseif(xt_prime.lt.1.)then
          call form_nahso4(store,jp,ibin)
          call form_nh4hso4(store,jp,ibin)
          call form_h2so4(store,jp,ibin)
        endif

10      if(jp .eq. jsolid)then
          call degas_hno3(store,jp,ibin)
          call degas_hcl(store,jp,ibin)
          call degas_nh3(store,jp,ibin)
        else
          call form_hno3(store,jp,ibin)
          call form_hcl(store,jp,ibin)
          call degas_nh3(store,jp,ibin)
        endif

      endif 



      call electrolytes_to_ions(jp, ibin)
















      return
      end subroutine form_electrolytes




















      subroutine form_caso4(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jcaso4,jp,ibin) = min(store(ica_a),store(iso4_a))
      store(ica_a)  = store(ica_a) - electrolyte(jcaso4,jp,ibin)
      store(iso4_a) = store(iso4_a) - electrolyte(jcaso4,jp,ibin)
      store(ica_a)  = max(0.D0, store(ica_a))
      store(iso4_a) = max(0.D0, store(iso4_a))

      return
      end subroutine form_caso4



      subroutine form_camsa2(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)
      
      electrolyte(jcamsa2,jp,ibin) = min(store(ica_a),0.5*store(imsa_a))
      store(ica_a)  = store(ica_a) - electrolyte(jcamsa2,jp,ibin)
      store(imsa_a) = store(imsa_a) - 2.d0*electrolyte(jcamsa2,jp,ibin)
      store(ica_a)  = max(0.D0, store(ica_a))
      store(imsa_a) = max(0.D0, store(imsa_a))

      return
      end subroutine form_camsa2



      subroutine form_cano3(store,jp,ibin)	



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jcano3,jp,ibin) = min(store(ica_a),0.5*store(ino3_a))

      store(ica_a)  = store(ica_a) - electrolyte(jcano3,jp,ibin)
      store(ino3_a) = store(ino3_a) - 2.*electrolyte(jcano3,jp,ibin)
      store(ica_a)  = max(0.D0, store(ica_a))
      store(ino3_a) = max(0.D0, store(ino3_a))

      return
      end subroutine form_cano3


      subroutine form_cacl2(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jcacl2,jp,ibin) = min(store(ica_a),0.5*store(icl_a))

      store(ica_a)  = store(ica_a) - electrolyte(jcacl2,jp,ibin)
      store(icl_a)  = store(icl_a) - 2.*electrolyte(jcacl2,jp,ibin)
      store(ica_a)  = max(0.D0, store(ica_a))
      store(icl_a)  = max(0.D0, store(icl_a))

      return
      end subroutine form_cacl2


      subroutine form_caco3(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      if(jp.eq.jtotal .or. jp.eq.jsolid)then
      electrolyte(jcaco3,jp,ibin) = store(ica_a)

      aer(ico3_a,jp,ibin)= electrolyte(jcaco3,jp,ibin)	

      store(ica_a) = 0.0
      store(ico3_a)= 0.0
      endif

      return
      end subroutine form_caco3


      subroutine form_na2so4(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jna2so4,jp,ibin) = min(.5*store(ina_a),   &
                                            store(iso4_a))
      store(ina_a) = store(ina_a) - 2.*electrolyte(jna2so4,jp,ibin)
      store(iso4_a)= store(iso4_a) - electrolyte(jna2so4,jp,ibin)
      store(ina_a) = max(0.D0, store(ina_a))
      store(iso4_a)= max(0.D0, store(iso4_a))

      return
      end subroutine form_na2so4



      subroutine form_nahso4(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnahso4,jp,ibin) = min(store(ina_a),   &
                                         store(iso4_a))
      store(ina_a)  = store(ina_a) - electrolyte(jnahso4,jp,ibin)
      store(iso4_a) = store(iso4_a) - electrolyte(jnahso4,jp,ibin)
      store(ina_a)  = max(0.D0, store(ina_a))
      store(iso4_a) = max(0.D0, store(iso4_a))

      return
      end subroutine form_nahso4



      subroutine form_namsa(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnamsa,jp,ibin) = min(store(ina_a), &
                                        store(imsa_a))
      store(ina_a)  = store(ina_a) - electrolyte(jnamsa,jp,ibin)
      store(imsa_a) = store(imsa_a) - electrolyte(jnamsa,jp,ibin)
      store(ina_a)  = max(0.D0, store(ina_a))
      store(imsa_a) = max(0.D0, store(imsa_a))

      return
      end subroutine form_namsa



      subroutine form_nano3(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnano3,jp,ibin)=min(store(ina_a),store(ino3_a))
      store(ina_a)  = store(ina_a) - electrolyte(jnano3,jp,ibin)
      store(ino3_a) = store(ino3_a) - electrolyte(jnano3,jp,ibin)
      store(ina_a)  = max(0.D0, store(ina_a))
      store(ino3_a) = max(0.D0, store(ino3_a))

      return
      end subroutine form_nano3



      subroutine form_nacl(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnacl,jp,ibin) = store(ina_a)

      store(ina_a) = 0.0
      store(icl_a) = store(icl_a) - electrolyte(jnacl,jp,ibin)
     
      if(store(icl_a) .lt. 0.)then 				
        aer(icl_a,jp,ibin)= aer(icl_a,jp,ibin)- store(icl_a)	

        if(jp .ne. jtotal)then
          aer(icl_a,jtotal,ibin)= aer(icl_a,jliquid,ibin)+ &		
                                  aer(icl_a,jsolid,ibin) 
        endif

        gas(ihcl_g) = gas(ihcl_g) + store(icl_a)			

        if(gas(ihcl_g) .lt. 0.0)then
          total_species(ihcl_g) = total_species(ihcl_g) - gas(ihcl_g)	
          tot_cl_in = tot_cl_in - gas(ihcl_g)				
        endif

        gas(ihcl_g) = max(0.D0, gas(ihcl_g))				
        store(icl_a) = 0.        				

      endif
     
      store(icl_a) = max(0.D0, store(icl_a))

      return
      end subroutine form_nacl



      subroutine form_nh4so4(store,jp,ibin)	



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnh4so4,jp,ibin)= min(.5*store(inh4_a),   &
                                           store(iso4_a))
      store(inh4_a)= store(inh4_a) - 2.*electrolyte(jnh4so4,jp,ibin)
      store(iso4_a)= store(iso4_a) - electrolyte(jnh4so4,jp,ibin)
      store(inh4_a) = max(0.D0, store(inh4_a))
      store(iso4_a) = max(0.D0, store(iso4_a))

      return
      end subroutine form_nh4so4



      subroutine form_nh4hso4(store,jp,ibin)	



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnh4hso4,jp,ibin) = min(store(inh4_a),   &
                                          store(iso4_a))
      store(inh4_a)= store(inh4_a) - electrolyte(jnh4hso4,jp,ibin)
      store(iso4_a)= store(iso4_a) - electrolyte(jnh4hso4,jp,ibin)
      store(inh4_a) = max(0.D0, store(inh4_a))
      store(iso4_a) = max(0.D0, store(iso4_a))

      return
      end subroutine form_nh4hso4



      subroutine form_nh4msa(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnh4msa,jp,ibin) = min(store(inh4_a), &
                                         store(imsa_a))
      store(inh4_a) = store(inh4_a) - electrolyte(jnh4msa,jp,ibin)
      store(imsa_a) = store(imsa_a) - electrolyte(jnh4msa,jp,ibin)
      store(inh4_a) = max(0.D0, store(inh4_a))
      store(imsa_a) = max(0.D0, store(imsa_a))

      return
      end subroutine form_nh4msa



      subroutine form_nh4cl(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnh4cl,jp,ibin) = min(store(inh4_a),   &
                                        store(icl_a))
      store(inh4_a) = store(inh4_a) - electrolyte(jnh4cl,jp,ibin)
      store(icl_a)  = store(icl_a) - electrolyte(jnh4cl,jp,ibin)
      store(inh4_a) = max(0.D0, store(inh4_a))
      store(icl_a)  = max(0.D0, store(icl_a))

      return
      end subroutine form_nh4cl



      subroutine form_nh4no3(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnh4no3,jp,ibin) = min(store(inh4_a),   &
                                         store(ino3_a))
      store(inh4_a) = store(inh4_a) - electrolyte(jnh4no3,jp,ibin)
      store(ino3_a) = store(ino3_a) - electrolyte(jnh4no3,jp,ibin)
      store(inh4_a) = max(0.D0, store(inh4_a))
      store(ino3_a) = max(0.D0, store(ino3_a))

      return
      end subroutine form_nh4no3



      subroutine form_nh4so4_lvcite(store,jp,ibin) 



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jnh4so4,jp,ibin)= 2.*store(inh4_a) - 3.*store(iso4_a)
      electrolyte(jlvcite,jp,ibin)= 2.*store(iso4_a) - store(inh4_a)
      electrolyte(jnh4so4,jp,ibin)= max(0.D0,   &
                                    electrolyte(jnh4so4,jp,ibin))
      electrolyte(jlvcite,jp,ibin)= max(0.D0,   &
                                    electrolyte(jlvcite,jp,ibin))
      store(inh4_a) = 0.
      store(iso4_a) = 0.

      return
      end subroutine form_nh4so4_lvcite



      subroutine form_lvcite_nh4hso4(store,jp,ibin) 



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jlvcite,jp,ibin) = store(inh4_a) - store(iso4_a)
      electrolyte(jnh4hso4,jp,ibin)= 3.*store(iso4_a) - 2.*store(inh4_a)
      electrolyte(jlvcite,jp,ibin) = max(0.D0,   &
                                      electrolyte(jlvcite,jp,ibin))
      electrolyte(jnh4hso4,jp,ibin)= max(0.D0,   &
                                      electrolyte(jnh4hso4,jp,ibin))
      store(inh4_a) = 0.
      store(iso4_a) = 0.

      return
      end subroutine form_lvcite_nh4hso4



      subroutine form_na2so4_nahso4(store,jp,ibin) 



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jna2so4,jp,ibin)= store(ina_a) - store(iso4_a)
      electrolyte(jnahso4,jp,ibin)= 2.*store(iso4_a) - store(ina_a)
      electrolyte(jna2so4,jp,ibin)= max(0.D0,   &
                                    electrolyte(jna2so4,jp,ibin))
      electrolyte(jnahso4,jp,ibin)= max(0.D0,   &
                                    electrolyte(jnahso4,jp,ibin))
      store(ina_a)  = 0.
      store(iso4_a) = 0.



      return
      end subroutine form_na2so4_nahso4




      subroutine form_h2so4(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jh2so4,jp,ibin) = max(0.0D0, store(iso4_a))
      store(iso4_a) = 0.0

      return
      end subroutine form_h2so4




      subroutine form_msa(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jmsa,jp,ibin) = max(0.0D0, store(imsa_a))
      store(imsa_a) = 0.0

      return
      end subroutine form_msa



      subroutine form_hno3(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jhno3,jp,ibin) = max(0.0D0, store(ino3_a))
      store(ino3_a) = 0.0

      return
      end subroutine form_hno3




      subroutine form_hcl(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      electrolyte(jhcl,jp,ibin) = max(0.0D0, store(icl_a))
      store(icl_a) = 0.0

      return
      end subroutine form_hcl




      subroutine degas_hno3(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      store(ino3_a) = max(0.0D0, store(ino3_a))
      gas(ihno3_g) = gas(ihno3_g) + store(ino3_a)
      aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) - store(ino3_a)
      aer(ino3_a,jp,ibin) = max(0.0D0,aer(ino3_a,jp,ibin))


      if(jp .ne. jtotal)then
        aer(ino3_a,jtotal,ibin) = aer(ino3_a,jsolid, ibin) +   &
                                  aer(ino3_a,jliquid,ibin)
      endif

      electrolyte(jhno3,jp,ibin) = 0.0
      store(ino3_a) = 0.0

      return
      end subroutine degas_hno3



      subroutine degas_hcl(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      store(icl_a) = max(0.0D0, store(icl_a))
      gas(ihcl_g) = gas(ihcl_g) + store(icl_a)
      aer(icl_a,jp,ibin) = aer(icl_a,jp,ibin) - store(icl_a)
      aer(icl_a,jp,ibin) = max(0.0D0,aer(icl_a,jp,ibin))


      if(jp .ne. jtotal)then
        aer(icl_a,jtotal,ibin) = aer(icl_a,jsolid, ibin) +   &
                                 aer(icl_a,jliquid,ibin)
      endif

      electrolyte(jhcl,jp,ibin) = 0.0
      store(icl_a) = 0.0

      return
      end subroutine degas_hcl



      subroutine degas_nh3(store,jp,ibin)



      integer jp, ibin
      real(kind=8) store(naer)

      store(inh4_a) = max(0.0D0, store(inh4_a))
      gas(inh3_g) = gas(inh3_g) + store(inh4_a)
      aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) - store(inh4_a)
      aer(inh4_a,jp,ibin) = max(0.0D0,aer(inh4_a,jp,ibin))


      if(jp .ne. jtotal)then
        aer(inh4_a,jtotal,ibin)= aer(inh4_a,jsolid, ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      endif

      store(inh4_a) = 0.0

      return
      end subroutine degas_nh3









      subroutine degas_acids(jp,ibin,xt)



      integer jp, ibin
      real(kind=8) xt

      real(kind=8) ehno3, ehcl



      if(jp .ne. jliquid)then
        if (iprint_mosaic_diag1 .gt. 0) then
          write(6,*)'mosaic - error in degas_acids'
          write(6,*)'wrong jp'
        endif
      endif

      ehno3 = electrolyte(jhno3,jp,ibin)
      ehcl  = electrolyte(jhcl,jp,ibin)


      gas(ihno3_g) = gas(ihno3_g) + ehno3
      gas(ihcl_g)  = gas(ihcl_g)  + ehcl


      aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) - ehno3
      aer(icl_a, jp,ibin) = aer(icl_a, jp,ibin) - ehcl


      aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
                                aer(ino3_a,jsolid, ibin)

      aer(icl_a,jtotal,ibin)  = aer(icl_a,jliquid,ibin) +   &
                                aer(icl_a,jsolid, ibin)

      electrolyte(jhno3,jp,ibin) = 0.0
      electrolyte(jhcl,jp,ibin)  = 0.0

      return
      end subroutine degas_acids






















      subroutine degas_solid_nh4no3(ibin)



      integer ibin

      integer jp
      real(kind=8) a, b, c, xgas, xt



      jp = jsolid

      a = 1.0
      b = gas(inh3_g) + gas(ihno3_g)
      c = gas(inh3_g)*gas(ihno3_g) - keq_sg(1)
      xgas = quadratic(a,b,c)

      if(xgas .ge. electrolyte(jnh4no3,jp,ibin))then 

          gas(inh3_g) = gas(inh3_g)  + electrolyte(jnh4no3,jp,ibin)
          gas(ihno3_g)= gas(ihno3_g) + electrolyte(jnh4no3,jp,ibin)
          aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) -   &
                                electrolyte(jnh4no3,jp,ibin)
          aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) -   &
                                electrolyte(jnh4no3,jp,ibin)

      else	

          gas(inh3_g) = gas(inh3_g)  + xgas
          gas(ihno3_g)= gas(ihno3_g) + xgas
          aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) - xgas
          aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) - xgas
      endif



      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)

      return
      end subroutine degas_solid_nh4no3










      subroutine degas_solid_nh4cl(ibin)



      integer ibin

      integer jp
      real(kind=8) a, b, c, xgas, xt



      jp = jsolid

      a = 1.0
      b = gas(inh3_g) + gas(ihcl_g)
      c = gas(inh3_g)*gas(ihcl_g) - keq_sg(2)
      xgas = quadratic(a,b,c)

      if(xgas .ge. electrolyte(jnh4cl,jp,ibin))then 

          gas(inh3_g) = gas(inh3_g) + electrolyte(jnh4cl,jp,ibin)
          gas(ihcl_g) = gas(ihcl_g) + electrolyte(jnh4cl,jp,ibin)
          aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) -   &
                                electrolyte(jnh4cl,jp,ibin)
          aer(icl_a,jp,ibin)  = aer(icl_a,jp,ibin) -   &
                                electrolyte(jnh4cl,jp,ibin)

      else	

          gas(inh3_g) = gas(inh3_g) + xgas
          gas(ihcl_g) = gas(ihcl_g) + xgas
          aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) - xgas
          aer(icl_a,jp,ibin)  = aer(icl_a,jp,ibin)  - xgas

      endif



      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
                                 aer(icl_a,jliquid,ibin)

      return
      end subroutine degas_solid_nh4cl



















      subroutine absorb_tiny_nh4no3(ibin)



      integer ibin

      real(kind=8) small_aer, small_gas, small_amt
      integer je					



      electrolyte_sum(jtotal,ibin) = 0.0	
      do je = 1, nelectrolyte
        electrolyte_sum(jtotal,ibin) = electrolyte_sum(jtotal,ibin) + &
                                       electrolyte(je,jtotal,ibin)
      enddo


      small_gas = 0.01 * min(gas(inh3_g), gas(ihno3_g))
      small_aer = 0.01 * electrolyte_sum(jtotal,ibin)
      if(small_aer .eq. 0.0)small_aer = small_gas

      small_amt = min(small_gas, small_aer)

      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) + small_amt
      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) + small_amt


      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)


      gas(inh3_g)    = gas(inh3_g) - small_amt
      gas(ihno3_g)   = gas(ihno3_g) - small_amt

      return
      end subroutine absorb_tiny_nh4no3








      subroutine absorb_tiny_nh4cl(ibin)



      integer ibin

      real(kind=8) small_aer, small_gas, small_amt
	integer je					



      electrolyte_sum(jtotal,ibin) = 0.0	
      do je = 1, nelectrolyte
        electrolyte_sum(jtotal,ibin) = electrolyte_sum(jtotal,ibin) + &
                                       electrolyte(je,jtotal,ibin)
      enddo


      small_gas = 0.01 * min(gas(inh3_g), gas(ihcl_g))
      small_aer = 0.01 * electrolyte_sum(jtotal,ibin)
      if(small_aer .eq. 0.0)small_aer = small_gas

      small_amt = min(small_gas, small_aer)

      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) + small_amt
      aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin)  + small_amt


      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
                                 aer(icl_a,jliquid,ibin)


      gas(inh3_g)   = gas(inh3_g) - small_amt
      gas(ihcl_g)   = gas(ihcl_g) - small_amt

      return
      end subroutine absorb_tiny_nh4cl















      subroutine degas_tiny_nh4no3(ibin)



      integer ibin

      real(kind=8) small_amt

      small_amt = 0.01 * electrolyte(jnh4no3,jliquid,ibin)

      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) - small_amt
      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) - small_amt


      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)


      gas(inh3_g)  = gas(inh3_g)  + small_amt
      gas(ihno3_g) = gas(ihno3_g) + small_amt

      return
      end subroutine degas_tiny_nh4no3






      subroutine degas_tiny_nh4cl(ibin)



      integer ibin

      real(kind=8) small_amt


      small_amt = 0.01 * electrolyte(jnh4cl,jliquid,ibin)

      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) - small_amt
      aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin) - small_amt


      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
                                 aer(icl_a,jliquid,ibin)


      gas(inh3_g) = gas(inh3_g) + small_amt
      gas(ihcl_g) = gas(ihcl_g) + small_amt

      return
      end subroutine degas_tiny_nh4cl









      subroutine absorb_tiny_hcl(ibin)	



      integer ibin

      real(kind=8) small_aer, small_amt, small_gas

      small_gas = 0.01 * gas(ihcl_g)
      small_aer = 0.01 * aer(ino3_a,jliquid,ibin)

      small_amt = min(small_gas, small_aer)


      aer(icl_a,jliquid,ibin)= aer(icl_a,jliquid,ibin) + small_amt
      aer(icl_a,jtotal,ibin) = aer(icl_a,jsolid,ibin) +   &
                               aer(icl_a,jliquid,ibin)
      gas(ihcl_g) = gas(ihcl_g) - small_amt


      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) - small_amt
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)


      gas(ihno3_g) = gas(ihno3_g) + small_amt

      return
      end subroutine absorb_tiny_hcl





      subroutine absorb_tiny_hno3(ibin)	



      integer ibin

      real(kind=8) small_aer, small_amt, small_gas

      small_gas = 0.01 * gas(ihno3_g)
      small_aer = 0.01 * aer(icl_a,jliquid,ibin)

      small_amt = min(small_gas, small_aer)


      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) + small_amt
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)
      gas(ihno3_g) = gas(ihno3_g) - small_amt


      aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin) - small_amt
      aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin) +   &
                                 aer(icl_a,jliquid,ibin)


      gas(ihcl_g) = gas(ihcl_g) + small_amt

      return
      end subroutine absorb_tiny_hno3















      subroutine equilibrate_acids(ibin)



      integer ibin



      if(gas(ihcl_g)*gas(ihno3_g) .gt. 0.)then
        call equilibrate_hcl_and_hno3(ibin)
      elseif(gas(ihcl_g) .gt. 0.)then
        call equilibrate_hcl(ibin)
      elseif(gas(ihno3_g) .gt. 0.)then
        call equilibrate_hno3(ibin)
      endif


      return
      end subroutine equilibrate_acids









      subroutine equilibrate_hcl(ibin)



      integer ibin

      real(kind=8) a, aerh, aerhso4, aerso4, b, c, dum, kdash_hcl, mh, tcl,   &
        w, xt, z


      aerso4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
      aerhso4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

      tcl = aer(icl_a,jliquid,ibin) + gas(ihcl_g)		
      kdash_hcl = keq_gl(4)*1.e+18/gam(jhcl,ibin)**2	
      z = (   aer(ina_a, jliquid,ibin) + 		   &  
              aer(inh4_a,jliquid,ibin) +   &
           2.*aer(ica_a, jliquid,ibin) ) -   &
          (2.*aerso4  +   &
              aerhso4 +   &
              aer(ino3_a,jliquid,ibin) )


      w     = water_a(ibin)				

      kdash_hcl = keq_gl(4)*1.e+18/gam(jhcl,ibin)**2	
      a = 1.0
      b = (kdash_hcl*w + z/w)*1.e-9
      c = kdash_hcl*(z - tcl)*1.e-18


      dum = b*b - 4.*a*c
      if (dum .lt. 0.) return		


      if(c .lt. 0.)then
        mh = quadratic(a,b,c)	
        aerh = mh*w*1.e+9
        aer(icl_a,jliquid,ibin) = aerh + z
      else
        mh = sqrt(keq_ll(3))
      endif

      call form_electrolytes(jliquid,ibin,xt)


      gas(ihcl_g) = tcl - aer(icl_a,jliquid,ibin)



      ma(ja_so4,ibin)  = 1.e-9*aerso4/water_a(ibin)
      ma(ja_hso4,ibin) = 1.e-9*aerhso4/water_a(ibin)
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

      mc(jc_h,ibin)    = mh
      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)



      activity(jhcl,ibin)    = mc(jc_h,ibin)  *ma(ja_cl,ibin)  *   &
                               gam(jhcl,ibin)**2

      activity(jhno3,ibin)   = mc(jc_h,ibin)  *ma(ja_no3,ibin) *   &
                               gam(jhno3,ibin)**2

      activity(jnh4cl,ibin)  = mc(jc_nh4,ibin)*ma(ja_cl,ibin) *   &
                               gam(jnh4cl,ibin)**2



      aer(icl_a,jtotal,ibin) = aer(icl_a,jliquid,ibin) +   &
                               aer(icl_a,jsolid,ibin)

      electrolyte(jhcl,jtotal,ibin) = electrolyte(jhcl,jliquid,ibin)

      return
      end subroutine equilibrate_hcl





      subroutine equilibrate_hno3(ibin)



      integer ibin

      real(kind=8) a, aerh, aerhso4, aerso4, b, c, dum, kdash_hno3, mh,   &
        tno3, w, xt, z


      aerso4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
      aerhso4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

      tno3 = aer(ino3_a,jliquid,ibin) + gas(ihno3_g)	
      kdash_hno3 = keq_gl(3)*1.e+18/gam(jhno3,ibin)**2	
      z = (   aer(ina_a, jliquid,ibin) + 		   &  
              aer(inh4_a,jliquid,ibin) +   &
           2.*aer(ica_a, jliquid,ibin) ) -   &
          (2.*aerso4  +   &
              aerhso4 +   &
              aer(icl_a,jliquid,ibin) )


      w     = water_a(ibin)				

      kdash_hno3 = keq_gl(3)*1.e+18/gam(jhno3,ibin)**2	
      a = 1.0
      b = (kdash_hno3*w + z/w)*1.e-9
      c = kdash_hno3*(z - tno3)*1.e-18

      dum = b*b - 4.*a*c
      if (dum .lt. 0.) return		



      if(c .lt. 0.)then
        mh = quadratic(a,b,c)	
        aerh = mh*w*1.e+9
        aer(ino3_a,jliquid,ibin) = aerh + z
      else
        mh = sqrt(keq_ll(3))
      endif

      call form_electrolytes(jliquid,ibin,xt)


      gas(ihno3_g)= tno3 - aer(ino3_a,jliquid,ibin)



      ma(ja_so4,ibin)  = 1.e-9*aerso4/water_a(ibin)
      ma(ja_hso4,ibin) = 1.e-9*aerhso4/water_a(ibin)
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

      mc(jc_h,ibin)    = mh
      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)



      activity(jhcl,ibin)    = mc(jc_h,ibin)  *ma(ja_cl,ibin)  *   &
                               gam(jhcl,ibin)**2

      activity(jhno3,ibin)   = mc(jc_h,ibin)  *ma(ja_no3,ibin) *   &
                               gam(jhno3,ibin)**2

      activity(jnh4no3,ibin) = mc(jc_nh4,ibin)*ma(ja_no3,ibin) *   &
                               gam(jnh4no3,ibin)**2



      aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
                                aer(ino3_a,jsolid,ibin)

      electrolyte(jhno3,jtotal,ibin) = electrolyte(jhno3,jliquid,ibin)

      return
      end subroutine equilibrate_hno3











      subroutine equilibrate_hcl_and_hno3(ibin)



      integer ibin

      real(kind=8) aerh, aerhso4, aerso4, kdash_hcl, kdash_hno3,   &
        mh, p, q, r, tcl, tno3, w, xt, z



      aerso4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
      aerhso4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

      tcl  = aer(icl_a,jliquid,ibin)  + gas(ihcl_g)	
      tno3 = aer(ino3_a,jliquid,ibin) + gas(ihno3_g)	

      kdash_hcl  = keq_gl(4)*1.e+18/gam(jhcl,ibin)**2	
      kdash_hno3 = keq_gl(3)*1.e+18/gam(jhno3,ibin)**2	

      z = (   aer(ina_a, jliquid,ibin) + 		   &  
              aer(inh4_a,jliquid,ibin) +   &
           2.*aer(ica_a, jliquid,ibin) ) -   &
          (2.*aerso4 + aerhso4 )


      w = water_a(ibin)

      kdash_hcl  = keq_gl(4)*1.e+18/gam(jhcl,ibin)**2	
      kdash_hno3 = keq_gl(3)*1.e+18/gam(jhno3,ibin)**2	

      p = (z/w + w*(kdash_hcl + kdash_hno3))*1.e-9

      q = 1.e-18*kdash_hcl*kdash_hno3*w**2  +   &
          1.e-18*z*(kdash_hcl + kdash_hno3) -   &
          1.e-18*kdash_hcl*tcl -   &
          1.e-18*kdash_hno3*tno3

      r = 1.e-18*kdash_hcl*kdash_hno3*w*(z - tcl - tno3)*1.e-9

      mh = cubic(p,q,r)

      if(mh .gt. 0.0)then
        aerh = mh*w*1.e+9
        aer(ino3_a,jliquid,ibin) = kdash_hno3*w*w*tno3/   &
                                  (aerh + kdash_hno3*w*w)
        aer(icl_a, jliquid,ibin) = kdash_hcl*w*w*tcl/   &
                                  (aerh + kdash_hcl*w*w)
      else
        mh = sqrt(keq_ll(3))
      endif

      call form_electrolytes(jliquid,ibin,xt)


      gas(ihno3_g)= tno3 - aer(ino3_a,jliquid,ibin)
      gas(ihcl_g) = tcl  - aer(icl_a,jliquid,ibin)



      ma(ja_so4,ibin)  = 1.e-9*aerso4/water_a(ibin)
      ma(ja_hso4,ibin) = 1.e-9*aerhso4/water_a(ibin)
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

      mc(jc_h,ibin)    = mh
      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)



      activity(jhcl,ibin)    = mc(jc_h,ibin)*ma(ja_cl,ibin)   *   &
                               gam(jhcl,ibin)**2

      activity(jhno3,ibin)   = mc(jc_h,ibin)*ma(ja_no3,ibin)  *   &
                               gam(jhno3,ibin)**2

      activity(jnh4no3,ibin) = mc(jc_nh4,ibin)*ma(ja_no3,ibin)*   &
                               gam(jnh4no3,ibin)**2

      activity(jnh4cl,ibin)  = mc(jc_nh4,ibin)*ma(ja_cl,ibin) *   &
                               gam(jnh4cl,ibin)**2



      aer(icl_a,jtotal,ibin)  = aer(icl_a,jliquid,ibin) +   &
                                aer(icl_a,jsolid,ibin)

      aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
                                aer(ino3_a,jsolid,ibin)

      electrolyte(jhno3,jtotal,ibin) = electrolyte(jhno3,jliquid,ibin)
      electrolyte(jhcl, jtotal,ibin) = electrolyte(jhcl, jliquid,ibin)

      return
      end subroutine equilibrate_hcl_and_hno3




















      subroutine load_mosaic_parameters




      integer iaer, je, ja, j_index, ibin



      logical, save :: first = .true.



      if(first)then
        first=.false.



      msize_framework = msection	
      mgas_aer_xfer   = myes		


      nmax_astem      = 200		
      alpha_astem     = 0.05		

      rtol_eqb_astem  = 0.01		
      ptol_mol_astem  = 0.01		


      nmax_mesa       = 80		
      rtol_mesa       = 0.01		







      ih2so4_g	= 1	
      ihno3_g	= 2	
      ihcl_g	= 3	
      inh3_g	= 4	
      imsa_g	= 5	
      ipcg1_b_c_g =6
      ipcg2_b_c_g =7
      ipcg3_b_c_g =8
      ipcg4_b_c_g =9
      ipcg5_b_c_g =10
      ipcg6_b_c_g =11
      ipcg7_b_c_g =12
      ipcg8_b_c_g =13
      ipcg9_b_c_g =14
      ipcg1_b_o_g =15
      ipcg2_b_o_g =16
      ipcg3_b_o_g =17
      ipcg4_b_o_g =18
      ipcg5_b_o_g =19
      ipcg6_b_o_g =20
      ipcg7_b_o_g =21
      ipcg8_b_o_g =22
      ipcg9_b_o_g =23
      iopcg1_b_c_g =24
      iopcg2_b_c_g = 25
      iopcg3_b_c_g =26
      iopcg4_b_c_g =27
      iopcg5_b_c_g =28
      iopcg6_b_c_g =29
      iopcg7_b_c_g =30
      iopcg8_b_c_g =31
      iopcg1_b_o_g =32
      iopcg2_b_o_g =33
      iopcg3_b_o_g =34
      iopcg4_b_o_g =35
      iopcg5_b_o_g =36
      iopcg6_b_o_g =37
      iopcg7_b_o_g =38
      iopcg8_b_o_g =39
      ipcg1_f_c_g =40
      ipcg2_f_c_g =41
      ipcg3_f_c_g =42
      ipcg4_f_c_g =43
      ipcg5_f_c_g =44
      ipcg6_f_c_g =45
      ipcg7_f_c_g =46
      ipcg8_f_c_g =47
      ipcg9_f_c_g =48
      ipcg1_f_o_g =49
      ipcg2_f_o_g =50
      ipcg3_f_o_g =51
      ipcg4_f_o_g =52
      ipcg5_f_o_g =53
      ipcg6_f_o_g =54
      ipcg7_f_o_g =55
      ipcg8_f_o_g =56
      ipcg9_f_o_g =57
      iopcg1_f_c_g =58
      iopcg2_f_c_g =59
      iopcg3_f_c_g =60
      iopcg4_f_c_g =61
      iopcg5_f_c_g =62
      iopcg6_f_c_g =63
      iopcg7_f_c_g =64
      iopcg8_f_c_g =65
      iopcg1_f_o_g =66
      iopcg2_f_o_g =67
      iopcg3_f_o_g =68
      iopcg4_f_o_g =69
      iopcg5_f_o_g =70
      iopcg6_f_o_g =71
      iopcg7_f_o_g =72
      iopcg8_f_o_g =73
      ismpa_g =74
      ismpbb_g =75
      iant1_c_g =76
      iant2_c_g =77
      iant3_c_g =78
      iant4_c_g =79
      iant1_o_g =80
      iant2_o_g =81
      iant3_o_g =82
      iant4_o_g =83
      ibiog1_c_g =84
      ibiog2_c_g =85
      ibiog3_c_g =86
      ibiog4_c_g =87
      ibiog1_o_g =88
      ibiog2_o_g =89
      ibiog3_o_g =90
      ibiog4_o_g =91





      iasoaX_g=92
      iasoa1_g=93
      iasoa2_g=94
      iasoa3_g=95
      iasoa4_g=96
      ibsoaX_g=97
      ibsoa1_g=98
      ibsoa2_g=99
      ibsoa3_g=100
      ibsoa4_g=101
      in2o5_g    =102  
      iclno2_g   =103  

      igly       =104
      iho        =105





      iso4_a	=  1	
      ino3_a	=  2	
      icl_a	=  3	
      inh4_a	=  4	
      imsa_a	=  5	
      ipcg1_b_c_a =6
      ipcg2_b_c_a =7
      ipcg3_b_c_a =8
      ipcg4_b_c_a =9
      ipcg5_b_c_a =10
      ipcg6_b_c_a =11
      ipcg7_b_c_a =12
      ipcg8_b_c_a =13
      ipcg9_b_c_a =14
      ipcg1_b_o_a =15
      ipcg2_b_o_a =16
      ipcg3_b_o_a =17
      ipcg4_b_o_a =18
      ipcg5_b_o_a =19
      ipcg6_b_o_a =20
      ipcg7_b_o_a =21
      ipcg8_b_o_a =22
      ipcg9_b_o_a =23
      iopcg1_b_c_a =24
      iopcg2_b_c_a = 25
      iopcg3_b_c_a =26
      iopcg4_b_c_a =27
      iopcg5_b_c_a =28
      iopcg6_b_c_a =29
      iopcg7_b_c_a =30
      iopcg8_b_c_a =31
      iopcg1_b_o_a =32
      iopcg2_b_o_a = 33
      iopcg3_b_o_a =34
      iopcg4_b_o_a =35
      iopcg5_b_o_a =36
      iopcg6_b_o_a =37
      iopcg7_b_o_a =38
      iopcg8_b_o_a =39
      ipcg1_f_c_a =40
      ipcg2_f_c_a =41
      ipcg3_f_c_a =42
      ipcg4_f_c_a =43
      ipcg5_f_c_a =44
      ipcg6_f_c_a =45
      ipcg7_f_c_a =46
      ipcg8_f_c_a =47
      ipcg9_f_c_a =48
      ipcg1_f_o_a =49
      ipcg2_f_o_a =50
      ipcg3_f_o_a =51
      ipcg4_f_o_a =52
      ipcg5_f_o_a =53
      ipcg6_f_o_a =54
      ipcg7_f_o_a =55
      ipcg8_f_o_a =56
      ipcg9_f_o_a =57
      iopcg1_f_c_a =58
      iopcg2_f_c_a =59
      iopcg3_f_c_a =60
      iopcg4_f_c_a =61
      iopcg5_f_c_a =62
      iopcg6_f_c_a =63
      iopcg7_f_c_a =64
      iopcg8_f_c_a =65
      iopcg1_f_o_a =66
      iopcg2_f_o_a =67
      iopcg3_f_o_a =68
      iopcg4_f_o_a =69
      iopcg5_f_o_a =70
      iopcg6_f_o_a =71
      iopcg7_f_o_a =72
      iopcg8_f_o_a =73
      ismpa_a =74
      ismpbb_a =75
      iant1_c_a =76
      iant2_c_a =77
      iant3_c_a =78
      iant4_c_a =79
      iant1_o_a =80
      iant2_o_a =81
      iant3_o_a =82
      iant4_o_a =83
      ibiog1_c_a =84
      ibiog2_c_a =85
      ibiog3_c_a =86
      ibiog4_c_a =87
      ibiog1_o_a =88
      ibiog2_o_a =89
      ibiog3_o_a =90
      ibiog4_o_a =91








      iasoaX_a=92
      iasoa1_a=93
      iasoa2_a=94
      iasoa3_a=95
      iasoa4_a=96
      ibsoaX_a=97
      ibsoa1_a=98
      ibsoa2_a=99
      ibsoa3_a=100
      ibsoa4_a=101
      iglysoa_r1_a = 102
      iglysoa_r2_a = 103
      iglysoa_sfc_a = 104
      iglysoa_nh4_a = 105
      iglysoa_oh_a = 106

      ico3_a    = 107    
      ina_a     = 108
      ica_a     = 109
      ioin_a    = 110
      ioc_a     = 111
      ibc_a     = 112



      
      jnh4so4	=  1	
      jlvcite	=  2	
      jnh4hso4	=  3	
      jnh4msa	=  4	
      jnh4no3	=  5	
      jnh4cl	=  6	
      jna2so4	=  7	
      jna3hso4	=  8	
      jnahso4	=  9	
      jnamsa	= 10	
      jnano3	= 11	
      jnacl	= 12	
      jcano3	= 13	
      jcacl2	= 14	
      jcamsa2	= 15	
      jh2so4	= 16	
      jmsa	= 17	
      jhno3	= 18	
      jhcl	= 19	
      jhhso4	= 20	
      jcaso4	= 21	
      jcaco3	= 22	
      joc	= 23	
      jbc	= 24	
      join	= 25	
      jpcg1_b_c =26
      jpcg2_b_c =27
      jpcg3_b_c =28
      jpcg4_b_c =29
      jpcg5_b_c =30
      jpcg6_b_c =31
      jpcg7_b_c =32
      jpcg8_b_c =33
      jpcg9_b_c =34
      jpcg1_b_o =35
      jpcg2_b_o =36
      jpcg3_b_o =37
      jpcg4_b_o =38
      jpcg5_b_o =39
      jpcg6_b_o =40
      jpcg7_b_o =41
      jpcg8_b_o =42
      jpcg9_b_o =43
      jopcg1_b_c =44
      jopcg2_b_c =45
      jopcg3_b_c =46
      jopcg4_b_c =47
      jopcg5_b_c =48
      jopcg6_b_c =49
      jopcg7_b_c =50
      jopcg8_b_c =51
      jopcg1_b_o =52
      jopcg2_b_o =53
      jopcg3_b_o =54
      jopcg4_b_o =55
      jopcg5_b_o =56
      jopcg6_b_o =57
      jopcg7_b_o =58
      jopcg8_b_o =59
      jpcg1_f_c =60
      jpcg2_f_c =61
      jpcg3_f_c =62
      jpcg4_f_c =63
      jpcg5_f_c =64
      jpcg6_f_c =65
      jpcg7_f_c =66
      jpcg8_f_c =67
      jpcg9_f_c =68
      jpcg1_f_o =69
      jpcg2_f_o =70
      jpcg3_f_o =71
      jpcg4_f_o =72
      jpcg5_f_o =73
      jpcg6_f_o =74
      jpcg7_f_o =75
      jpcg8_f_o =76
      jpcg9_f_o =77
      jopcg1_f_c =78
      jopcg2_f_c =79
      jopcg3_f_c =80
      jopcg4_f_c =81
      jopcg5_f_c =82
      jopcg6_f_c =83
      jopcg7_f_c =84
      jopcg8_f_c =85
      jopcg1_f_o =86
      jopcg2_f_o =87
      jopcg3_f_o =88
      jopcg4_f_o =89
      jopcg5_f_o =90
      jopcg6_f_o =91
      jopcg7_f_o =92
      jopcg8_f_o =93
      jsmpa =94
      jsmpbb =95
      jant1_c =96
      jant2_c =97
      jant3_c =98
      jant4_c =99
      jant1_o =100
      jant2_o =101
      jant3_o =102
      jant4_o =103
      jbiog1_c =104
      jbiog2_c =105
      jbiog3_c =106
      jbiog4_c =107
      jbiog1_o =108
      jbiog2_o =109
      jbiog3_o =110
      jbiog4_o =111

      jasoaX=112
      jasoa1=113
      jasoa2=114
      jasoa3=115
      jasoa4=116
      jbsoaX=117
      jbsoa1=118
      jbsoa2=119
      jbsoa3=120
      jbsoa4=121
      jglysoa_r1  = 122
      jglysoa_r2  = 123
      jglysoa_sfc = 124
      jglysoa_nh4 = 125
      jglysoa_oh  = 126
      jh2o  = 127 



      jc_h	=  1
      jc_nh4	=  2
      jc_na	=  3
      jc_ca	=  4


      ja_hso4	=  1
      ja_so4  	=  2
      ja_no3  	=  3
      ja_cl   	=  4
      ja_msa	=  5










      aer_name(iso4_a) = 'so4'
      aer_name(ino3_a) = 'no3'
      aer_name(icl_a)  = 'cl '
      aer_name(inh4_a) = 'nh4'
      aer_name(ioc_a)  = 'oc '
      aer_name(imsa_a) = 'msa'
      aer_name(ico3_a) = 'co3'
      aer_name(ina_a)  = 'na '
      aer_name(ica_a)  = 'ca '
      aer_name(ibc_a)  = 'bc '
      aer_name(ioin_a) = 'oin'
      aer_name(ipcg1_b_c_a)="pcg1_b_c"
      aer_name(ipcg2_b_c_a)="pcg2_b_c"
      aer_name(ipcg3_b_c_a)="pcg3_b_c"
      aer_name(ipcg4_b_c_a)="pcg4_b_c"
      aer_name(ipcg5_b_c_a)="pcg5_b_c"
      aer_name(ipcg6_b_c_a)="pcg6_b_c"
      aer_name(ipcg7_b_c_a)="pcg7_b_c"
      aer_name(ipcg8_b_c_a)="pcg8_b_c"
      aer_name(ipcg9_b_c_a)="pcg9_b_c"
      aer_name(iopcg1_b_c_a)="opcg1_b_c"
      aer_name(iopcg2_b_c_a)="opcg2_b_c"
      aer_name(iopcg3_b_c_a)="opcg3_b_c"
      aer_name(iopcg4_b_c_a)="opcg4_b_c"
      aer_name(iopcg5_b_c_a)="opcg5_b_c"
      aer_name(iopcg6_b_c_a)="opcg6_b_c"
      aer_name(iopcg7_b_c_a)="opcg7_b_c"
      aer_name(iopcg8_b_c_a)="opcg8_b_c"
      aer_name(ipcg1_b_o_a)="pcg1_b_o"
      aer_name(ipcg2_b_o_a)="pcg2_b_o"
      aer_name(ipcg3_b_o_a)="pcg3_b_o"
      aer_name(ipcg4_b_o_a)="pcg4_b_o"
      aer_name(ipcg5_b_o_a)="pcg5_b_o"
      aer_name(ipcg6_b_o_a)="pcg6_b_o"
      aer_name(ipcg7_b_o_a)="pcg7_b_o"
      aer_name(ipcg8_b_o_a)="pcg8_b_o"
      aer_name(ipcg9_b_o_a)="pcg9_b_o"
      aer_name(iopcg1_b_o_a)="opcg1_b_o"
      aer_name(iopcg2_b_o_a)="opcg2_b_o"
      aer_name(iopcg3_b_o_a)="opcg3_b_o"
      aer_name(iopcg4_b_o_a)="opcg4_b_o"
      aer_name(iopcg5_b_o_a)="opcg5_b_o"
      aer_name(iopcg6_b_o_a)="opcg6_b_o"
      aer_name(iopcg7_b_o_a)="opcg7_b_o"
      aer_name(iopcg8_b_o_a)="opcg8_b_o"
      aer_name(ipcg1_f_c_a)="pcg1_f_c"
      aer_name(ipcg2_f_c_a)="pcg2_f_c"
      aer_name(ipcg3_f_c_a)="pcg3_f_c"
      aer_name(ipcg4_f_c_a)="pcg4_f_c"
      aer_name(ipcg5_f_c_a)="pcg5_f_c"
      aer_name(ipcg6_f_c_a)="pcg6_f_c"
      aer_name(ipcg7_f_c_a)="pcg7_f_c"
      aer_name(ipcg8_f_c_a)="pcg8_f_c"
      aer_name(ipcg9_f_c_a)="pcg9_f_c"
      aer_name(iopcg1_f_c_a)="opcg1_f_c"
      aer_name(iopcg2_f_c_a)="opcg2_f_c"
      aer_name(iopcg3_f_c_a)="opcg3_f_c"
      aer_name(iopcg4_f_c_a)="opcg4_f_c"
      aer_name(iopcg5_f_c_a)="opcg5_f_c"
      aer_name(iopcg6_f_c_a)="opcg6_f_c"
      aer_name(iopcg7_f_c_a)="opcg7_f_c"
      aer_name(iopcg8_f_c_a)="opcg8_f_c"
      aer_name(ipcg1_f_o_a)="pcg1_f_o"
      aer_name(ipcg2_f_o_a)="pcg2_f_o"
      aer_name(ipcg3_f_o_a)="pcg3_f_o"
      aer_name(ipcg4_f_o_a)="pcg4_f_o"
      aer_name(ipcg5_f_o_a)="pcg5_f_o"
      aer_name(ipcg6_f_o_a)="pcg6_f_o"
      aer_name(ipcg7_f_o_a)="pcg7_f_o"
      aer_name(ipcg8_f_o_a)="pcg8_f_o"
      aer_name(ipcg9_f_o_a)="pcg9_f_o"
      aer_name(iopcg1_f_o_a)="opcg1_f_o"
      aer_name(iopcg2_f_o_a)="opcg2_f_o"
      aer_name(iopcg3_f_o_a)="opcg3_f_o"
      aer_name(iopcg4_f_o_a)="opcg4_f_o"
      aer_name(iopcg5_f_o_a)="opcg5_f_o"
      aer_name(iopcg6_f_o_a)="opcg6_f_o"
      aer_name(iopcg7_f_o_a)="opcg7_f_o"
      aer_name(iopcg8_f_o_a)="opcg8_f_o"
      aer_name(ismpa_a)="smpa"
      aer_name(ismpbb_a)="smpbb"
      aer_name(iglysoa_r1_a)="glysoa_r1"
      aer_name(iglysoa_r2_a)="glysoa_r2"
      aer_name(iglysoa_sfc_a)="glysoa_sfc"
      aer_name(iglysoa_nh4_a)="glysoa_nh4"
      aer_name(iglysoa_oh_a)="glysoa_oh"
      aer_name(iant1_c_a)="ant1_c"
      aer_name(iant2_c_a)="ant2_c"
      aer_name(iant3_c_a)="ant3_c"
      aer_name(iant4_c_a)="ant4_c"
      aer_name(iant1_o_a)="ant1_o"
      aer_name(iant2_o_a)="ant2_o"
      aer_name(iant3_o_a)="ant3_o"
      aer_name(iant4_o_a)="ant4_o"
      aer_name(ibiog1_c_a)="biog1_c"
      aer_name(ibiog2_c_a)="biog2_c"
      aer_name(ibiog3_c_a)="biog3_c"
      aer_name(ibiog4_c_a)="biog4_c"
      aer_name(ibiog1_o_a)="biog1_o"
      aer_name(ibiog2_o_a)="biog2_o"
      aer_name(ibiog3_o_a)="biog3_o"
      aer_name(ibiog4_o_a)="biog4_o"
      aer_name(iasoaX_a)="asoaX"
      aer_name(iasoa1_a)="asoa1"
      aer_name(iasoa2_a)="asoa2"
      aer_name(iasoa3_a)="asoa3"
      aer_name(iasoa4_a)="asoa4"
      aer_name(ibsoaX_a)="bsoaX"
      aer_name(ibsoa1_a)="bsoa1"
      aer_name(ibsoa2_a)="bsoa2"
      aer_name(ibsoa3_a)="bsoa3"
      aer_name(ibsoa4_a)="bsoa4"


      gas_name(ih2so4_g) = 'h2so4'
      gas_name(ihno3_g)  = 'hno3 '
      gas_name(ihcl_g)   = 'hcl  '
      gas_name(inh3_g)   = 'nh3  '
      gas_name(imsa_g)   = "msa  "
      gas_name(ipcg1_b_c_g)="pcg1_b_c"
      gas_name(ipcg2_b_c_g)="pcg2_b_c"
      gas_name(ipcg3_b_c_g)="pcg3_b_c"
      gas_name(ipcg4_b_c_g)="pcg4_b_c"
      gas_name(ipcg5_b_c_g)="pcg5_b_c"
      gas_name(ipcg6_b_c_g)="pcg6_b_c"
      gas_name(ipcg7_b_c_g)="pcg7_b_c"
      gas_name(ipcg8_b_c_g)="pcg8_b_c"
      gas_name(ipcg9_b_c_g)="pcg9_b_c"
      gas_name(iopcg1_b_c_g)="opcg1_b_c"
      gas_name(iopcg2_b_c_g)="opcg2_b_c"
      gas_name(iopcg3_b_c_g)="opcg3_b_c"
      gas_name(iopcg4_b_c_g)="opcg4_b_c"
      gas_name(iopcg5_b_c_g)="opcg5_b_c"
      gas_name(iopcg6_b_c_g)="opcg6_b_c"
      gas_name(iopcg7_b_c_g)="opcg7_b_c"
      gas_name(iopcg8_b_c_g)="opcg8_b_c"
      gas_name(ipcg1_b_o_g)="pcg1_b_o"
      gas_name(ipcg2_b_o_g)="pcg2_b_o"
      gas_name(ipcg3_b_o_g)="pcg3_b_o"
      gas_name(ipcg4_b_o_g)="pcg4_b_o"
      gas_name(ipcg5_b_o_g)="pcg5_b_o"
      gas_name(ipcg6_b_o_g)="pcg6_b_o"
      gas_name(ipcg7_b_o_g)="pcg7_b_o"
      gas_name(ipcg8_b_o_g)="pcg8_b_o"
      gas_name(ipcg9_b_o_g)="pcg9_b_o"
      gas_name(iopcg1_b_o_g)="opcg1_b_o"
      gas_name(iopcg2_b_o_g)="opcg2_b_o"
      gas_name(iopcg3_b_o_g)="opcg3_b_o"
      gas_name(iopcg4_b_o_g)="opcg4_b_o"
      gas_name(iopcg5_b_o_g)="opcg5_b_o"
      gas_name(iopcg6_b_o_g)="opcg6_b_o"
      gas_name(iopcg7_b_o_g)="opcg7_b_o"
      gas_name(iopcg8_b_o_g)="opcg8_b_o"
      gas_name(ipcg1_f_c_g)="pcg1_f_c"
      gas_name(ipcg2_f_c_g)="pcg2_f_c"
      gas_name(ipcg3_f_c_g)="pcg3_f_c"
      gas_name(ipcg4_f_c_g)="pcg4_f_c"
      gas_name(ipcg5_f_c_g)="pcg5_f_c"
      gas_name(ipcg6_f_c_g)="pcg6_f_c"
      gas_name(ipcg7_f_c_g)="pcg7_f_c"
      gas_name(ipcg8_f_c_g)="pcg8_f_c"
      gas_name(ipcg9_f_c_g)="pcg9_f_c"
      gas_name(iopcg1_f_c_g)="opcg1_f_c"
      gas_name(iopcg2_f_c_g)="opcg2_f_c"
      gas_name(iopcg3_f_c_g)="opcg3_f_c"
      gas_name(iopcg4_f_c_g)="opcg4_f_c"
      gas_name(iopcg5_f_c_g)="opcg5_f_c"
      gas_name(iopcg6_f_c_g)="opcg6_f_c"
      gas_name(iopcg7_f_c_g)="opcg7_f_c"
      gas_name(iopcg8_f_c_g)="opcg8_f_c"
      gas_name(ipcg1_f_o_g)="pcg1_f_o"
      gas_name(ipcg2_f_o_g)="pcg2_f_o"
      gas_name(ipcg3_f_o_g)="pcg3_f_o"
      gas_name(ipcg4_f_o_g)="pcg4_f_o"
      gas_name(ipcg5_f_o_g)="pcg5_f_o"
      gas_name(ipcg6_f_o_g)="pcg6_f_o"
      gas_name(ipcg7_f_o_g)="pcg7_f_o"
      gas_name(ipcg8_f_o_g)="pcg8_f_o"
      gas_name(ipcg9_f_o_g)="pcg9_f_o"
      gas_name(iopcg1_f_o_g)="opcg1_f_o"
      gas_name(iopcg2_f_o_g)="opcg2_f_o"
      gas_name(iopcg3_f_o_g)="opcg3_f_o"
      gas_name(iopcg4_f_o_g)="opcg4_f_o"
      gas_name(iopcg5_f_o_g)="opcg5_f_o"
      gas_name(iopcg6_f_o_g)="opcg6_f_o"
      gas_name(iopcg7_f_o_g)="opcg7_f_o"
      gas_name(iopcg8_f_o_g)="opcg8_f_o"
      gas_name(ismpa_g)="smpa"
      gas_name(ismpbb_g)="smpbb"
      gas_name(iant1_c_g)="ant1_c"
      gas_name(iant2_c_g)="ant2_c"
      gas_name(iant3_c_g)="ant3_c"
      gas_name(iant4_c_g)="ant4_c"
      gas_name(iant1_o_g)="ant1_o"
      gas_name(iant2_o_g)="ant2_o"
      gas_name(iant3_o_g)="ant3_o"
      gas_name(iant4_o_g)="ant4_o"
      gas_name(ibiog1_c_g)="biog1_c"
      gas_name(ibiog2_c_g)="biog2_c"
      gas_name(ibiog3_c_g)="biog3_c"
      gas_name(ibiog4_c_g)="biog4_c"
      gas_name(ibiog1_o_g)="biog1_o"
      gas_name(ibiog2_o_g)="biog2_o"
      gas_name(ibiog3_o_g)="biog3_o"
      gas_name(ibiog4_o_g)="biog4_o"
      gas_name(in2o5_g) = "n2o5 "
      gas_name(iclno2_g)= "clno2"
      gas_name(iasoaX_g)="asoaX"
      gas_name(iasoa1_g)="asoa1"
      gas_name(iasoa2_g)="asoa2"
      gas_name(iasoa3_g)="asoa3"
      gas_name(iasoa4_g)="asoa4"
      gas_name(ibsoaX_g)="bsoaX"
      gas_name(ibsoa1_g)="bsoa1"
      gas_name(ibsoa2_g)="bsoa2"
      gas_name(ibsoa3_g)="bsoa3"
      gas_name(ibsoa4_g)="bsoa4"
      gas_name(igly)="gly"
      gas_name(iho)="ho" 
      

      ename(jnh4so4) = 'amso4'
      ename(jlvcite) = '(nh4)3h(so4)2'
      ename(jnh4hso4)= 'nh4hso4'
      ename(jnh4msa) = "ch3so3nh4"
      ename(jnh4no3) = 'nh4no3'
      ename(jnh4cl)  = 'nh4cl'
      ename(jnacl)   = 'nacl'
      ename(jnano3)  = 'nano3'
      ename(jna2so4) = 'na2so4'
      ename(jna3hso4)= 'na3h(so4)2'
      ename(jnamsa)  = "ch3so3na"
      ename(jnahso4) = 'nahso4'
      ename(jcaso4)  = 'caso4'
      ename(jcamsa2) = "(ch3so3)2ca"
      ename(jcano3)  = 'ca(no3)2'
      ename(jcacl2)  = 'cacl2'
      ename(jcaco3)  = 'caco3'
      ename(jh2so4)  = 'h2so4'
      ename(jhhso4)  = 'hhso4'
      ename(jhno3)   = 'hno3'
      ename(jhcl)    = 'hcl'
      ename(jmsa)    = "ch3so3h"


      mw_electrolyte(jnh4so4) = 132.0
      mw_electrolyte(jlvcite) = 247.0
      mw_electrolyte(jnh4hso4)= 115.0
      mw_electrolyte(jnh4msa) = 113.0
      mw_electrolyte(jnh4no3) = 80.0
      mw_electrolyte(jnh4cl)  = 53.5
      mw_electrolyte(jnacl)   = 58.5
      mw_electrolyte(jnano3)  = 85.0
      mw_electrolyte(jna2so4) = 142.0
      mw_electrolyte(jna3hso4)= 262.0
      mw_electrolyte(jnahso4) = 120.0
      mw_electrolyte(jnamsa)  = 118.0
      mw_electrolyte(jcaso4)  = 136.0
      mw_electrolyte(jcamsa2) = 230.0
      mw_electrolyte(jcano3)  = 164.0
      mw_electrolyte(jcacl2)  = 111.0
      mw_electrolyte(jcaco3)  = 100.0
      mw_electrolyte(jh2so4)  = 98.0
      mw_electrolyte(jhno3)   = 63.0
      mw_electrolyte(jhcl)    = 36.5
      mw_electrolyte(jmsa)    = 96.0



      mw_c(jc_h)  =  1.0
      mw_c(jc_nh4)= 18.0
      mw_c(jc_na) = 23.0
      mw_c(jc_ca) = 40.0

      mw_a(ja_so4) = 96.0
      mw_a(ja_hso4)= 97.0
      mw_a(ja_no3) = 62.0
      mw_a(ja_cl)  = 35.5
      MW_a(ja_msa) = 95.0



      zc(jc_h)   = 1
      zc(jc_nh4) = 1
      zc(jc_na)  = 1
      zc(jc_ca)  = 2

      za(ja_hso4)= 1
      za(ja_so4) = 2
      za(ja_no3) = 1
      za(ja_cl)  = 1
      za(ja_msa) = 1



      dens_electrolyte(jnh4so4)  = 1.8
      dens_electrolyte(jlvcite)  = 1.8
      dens_electrolyte(jnh4hso4) = 1.8
      dens_electrolyte(jnh4msa)  = 1.8 
      dens_electrolyte(jnh4no3)  = 1.8
      dens_electrolyte(jnh4cl)   = 1.8
      dens_electrolyte(jnacl)    = 2.2
      dens_electrolyte(jnano3)   = 2.2
      dens_electrolyte(jna2so4)  = 2.2
      dens_electrolyte(jna3hso4) = 2.2
      dens_electrolyte(jnahso4)  = 2.2
      dens_electrolyte(jnamsa)   = 2.2 
      dens_electrolyte(jcaso4)   = 2.6
      dens_electrolyte(jcamsa2)  = 2.6	
      dens_electrolyte(jcano3)   = 2.6
      dens_electrolyte(jcacl2)   = 2.6
      dens_electrolyte(jcaco3)   = 2.6
      dens_electrolyte(jh2so4)   = 1.8
      dens_electrolyte(jhhso4)   = 1.8
      dens_electrolyte(jhno3)    = 1.8
      dens_electrolyte(jhcl)     = 1.8
      dens_electrolyte(jmsa)     = 1.8 



      dens_comp_a(jnh4so4)  = 1.8
      dens_comp_a(jlvcite)  = 1.8
      dens_comp_a(jnh4hso4) = 1.8
      dens_comp_a(jnh4msa)  = 1.8	
      dens_comp_a(jnh4no3)  = 1.7
      dens_comp_a(jnh4cl)   = 1.5
      dens_comp_a(jnacl)    = 2.2
      dens_comp_a(jnano3)   = 2.2
      dens_comp_a(jna2so4)  = 2.2
      dens_comp_a(jna3hso4) = 2.2
      dens_comp_a(jnahso4)  = 2.2
      dens_comp_a(jnamsa)   = 2.2	
      dens_comp_a(jcaso4)   = 2.6
      dens_comp_a(jcamsa2)  = 2.6	
      dens_comp_a(jcano3)   = 2.6
      dens_comp_a(jcacl2)   = 2.6
      dens_comp_a(jcaco3)   = 2.6
      dens_comp_a(jh2so4)   = 1.8
      dens_comp_a(jhhso4)   = 1.8
      dens_comp_a(jhno3)    = 1.8
      dens_comp_a(jhcl)     = 1.8
      dens_comp_a(jmsa)     = 1.8	
      dens_comp_a(joc)      = 1.0
      dens_comp_a(jbc)      = 1.8
      dens_comp_a(join)     = 2.6
      dens_comp_a(jh2o)     = 1.0
      dens_comp_a(ipcg1_b_c_a) =1.0
      dens_comp_a(ipcg2_b_c_a) =1.0
      dens_comp_a(ipcg3_b_c_a)=1.0
      dens_comp_a(ipcg4_b_c_a)=1.0
      dens_comp_a(ipcg5_b_c_a)=1.0
      dens_comp_a(ipcg6_b_c_a)=1.0
      dens_comp_a(ipcg7_b_c_a)=1.0
      dens_comp_a(ipcg8_b_c_a)=1.0
      dens_comp_a(ipcg9_b_c_a)=1.0
      dens_comp_a(iopcg1_b_c_a)=1.0
      dens_comp_a(iopcg2_b_c_a)=1.0
      dens_comp_a(iopcg3_b_c_a)=1.0
      dens_comp_a(iopcg4_b_c_a)=1.0
      dens_comp_a(iopcg5_b_c_a)=1.0
      dens_comp_a(iopcg6_b_c_a)=1.0
      dens_comp_a(iopcg7_b_c_a)=1.0
      dens_comp_a(iopcg8_b_c_a)=1.0
      dens_comp_a(ipcg1_b_o_a)=1.0
      dens_comp_a(ipcg2_b_o_a)=1.0
      dens_comp_a(ipcg3_b_o_a)=1.0
      dens_comp_a(ipcg4_b_o_a)=1.0
      dens_comp_a(ipcg5_b_o_a)=1.0
      dens_comp_a(ipcg6_b_o_a)=1.0
      dens_comp_a(ipcg7_b_o_a)=1.0
      dens_comp_a(ipcg8_b_o_a)=1.0
      dens_comp_a(ipcg9_b_o_a)=1.0
      dens_comp_a(iopcg1_b_o_a)=1.0
      dens_comp_a(iopcg2_b_o_a)=1.0
      dens_comp_a(iopcg3_b_o_a)=1.0
      dens_comp_a(iopcg4_b_o_a)=1.0
      dens_comp_a(iopcg5_b_o_a)=1.0
      dens_comp_a(iopcg6_b_o_a)=1.0
      dens_comp_a(iopcg7_b_o_a)=1.0
      dens_comp_a(iopcg8_b_o_a)=1.0
      dens_comp_a(ipcg1_f_c_a) =1.0
      dens_comp_a(ipcg2_f_c_a) =1.0
      dens_comp_a(ipcg3_f_c_a)=1.0
      dens_comp_a(ipcg4_f_c_a)=1.0
      dens_comp_a(ipcg5_f_c_a)=1.0
      dens_comp_a(ipcg6_f_c_a)=1.0
      dens_comp_a(ipcg7_f_c_a)=1.0
      dens_comp_a(ipcg8_f_c_a)=1.0
      dens_comp_a(ipcg9_f_c_a)=1.0
      dens_comp_a(iopcg1_f_c_a)=1.0
      dens_comp_a(iopcg2_f_c_a)=1.0
      dens_comp_a(iopcg3_f_c_a)=1.0
      dens_comp_a(iopcg4_f_c_a)=1.0
      dens_comp_a(iopcg5_f_c_a)=1.0
      dens_comp_a(iopcg6_f_c_a)=1.0
      dens_comp_a(iopcg7_f_c_a)=1.0
      dens_comp_a(iopcg8_f_c_a)=1.0
      dens_comp_a(ipcg1_f_o_a)=1.0
      dens_comp_a(ipcg2_f_o_a)=1.0
      dens_comp_a(ipcg3_f_o_a)=1.0
      dens_comp_a(ipcg4_f_o_a)=1.0
      dens_comp_a(ipcg5_f_o_a)=1.0
      dens_comp_a(ipcg6_f_o_a)=1.0
      dens_comp_a(ipcg7_f_o_a)=1.0
      dens_comp_a(ipcg8_f_o_a)=1.0
      dens_comp_a(ipcg9_f_o_a)=1.0
      dens_comp_a(iopcg1_f_o_a)=1.0
      dens_comp_a(iopcg2_f_o_a)=1.0
      dens_comp_a(iopcg3_f_o_a)=1.0
      dens_comp_a(iopcg4_f_o_a)=1.0
      dens_comp_a(iopcg5_f_o_a)=1.0
      dens_comp_a(iopcg6_f_o_a)=1.0
      dens_comp_a(iopcg7_f_o_a)=1.0
      dens_comp_a(iopcg8_f_o_a)=1.0
      dens_comp_a(ismpa_a)=1.0
      dens_comp_a(ismpbb_a)=1.0
      dens_comp_a(iglysoa_r1_a)=1.0
      dens_comp_a(iglysoa_r2_a)=1.0
      dens_comp_a(iglysoa_sfc_a)=1.0
      dens_comp_a(iglysoa_nh4_a)=1.0
      dens_comp_a(iglysoa_oh_a)=1.0
      dens_comp_a(iant1_c_a)=1.0
      dens_comp_a(iant2_c_a)=1.0
      dens_comp_a(iant3_c_a)=1.0
      dens_comp_a(iant4_c_a)=1.0
      dens_comp_a(iant1_o_a)=1.0
      dens_comp_a(iant2_o_a)=1.0
      dens_comp_a(iant3_o_a)=1.0
      dens_comp_a(iant4_o_a)=1.0
      dens_comp_a(ibiog1_c_a)=1.0
      dens_comp_a(ibiog2_c_a)=1.0
      dens_comp_a(ibiog3_c_a)=1.0
      dens_comp_a(ibiog4_c_a)=1.0
      dens_comp_a(ibiog1_o_a)=1.0
      dens_comp_a(ibiog2_o_a)=1.0
      dens_comp_a(ibiog3_o_a)=1.0
      dens_comp_a(ibiog4_o_a)=1.0
      dens_comp_a(iasoaX_a)=1.5
      dens_comp_a(iasoa1_a)=1.5
      dens_comp_a(iasoa2_a)=1.5
      dens_comp_a(iasoa3_a)=1.5
      dens_comp_a(iasoa4_a)=1.5
      dens_comp_a(ibsoaX_a)=1.5
      dens_comp_a(ibsoa1_a)=1.5
      dens_comp_a(ibsoa2_a)=1.5
      dens_comp_a(ibsoa3_a)=1.5
      dens_comp_a(ibsoa4_a)=1.5


      mw_aer_mac(iso4_a) = 96.0
      mw_aer_mac(ino3_a) = 62.0
      mw_aer_mac(icl_a)  = 35.5
      mw_aer_mac(imsa_a) = 95.0 
      mw_aer_mac(ico3_a) = 60.0
      mw_aer_mac(inh4_a) = 18.0
      mw_aer_mac(ina_a)  = 23.0
      mw_aer_mac(ica_a)  = 40.0
      mw_aer_mac(ioin_a) = 1.0          
      mw_aer_mac(ibc_a)  = 1.0          
      mw_aer_mac(ioc_a)  = 250.0  
      mw_aer_mac(ipcg1_b_c_a) =250.0
      mw_aer_mac(ipcg2_b_c_a) =250.0
      mw_aer_mac(ipcg3_b_c_a)=250.0
      mw_aer_mac(ipcg4_b_c_a)=250.0
      mw_aer_mac(ipcg5_b_c_a)=250.0
      mw_aer_mac(ipcg6_b_c_a)=250.0
      mw_aer_mac(ipcg7_b_c_a)=250.0
      mw_aer_mac(ipcg8_b_c_a)=250.0
      mw_aer_mac(ipcg9_b_c_a)=250.0
      mw_aer_mac(iopcg1_b_c_a)=250.0
      mw_aer_mac(iopcg2_b_c_a)=250.0
      mw_aer_mac(iopcg3_b_c_a)=250.0
      mw_aer_mac(iopcg4_b_c_a)=250.0
      mw_aer_mac(iopcg5_b_c_a)=250.0
      mw_aer_mac(iopcg6_b_c_a)=250.0
      mw_aer_mac(iopcg7_b_c_a)=250.0
      mw_aer_mac(iopcg8_b_c_a)=250.0
      mw_aer_mac(ipcg1_b_o_a)=250.0
      mw_aer_mac(ipcg2_b_o_a)=250.0
      mw_aer_mac(ipcg3_b_o_a)=250.0
      mw_aer_mac(ipcg4_b_o_a)=250.0
      mw_aer_mac(ipcg5_b_o_a)=250.0
      mw_aer_mac(ipcg6_b_o_a)=250.0
      mw_aer_mac(ipcg7_b_o_a)=250.0
      mw_aer_mac(ipcg8_b_o_a)=250.0
      mw_aer_mac(ipcg9_b_o_a)=250.0
      mw_aer_mac(iopcg1_b_o_a)=250.0
      mw_aer_mac(iopcg2_b_o_a)=250.0
      mw_aer_mac(iopcg3_b_o_a)=250.0
      mw_aer_mac(iopcg4_b_o_a)=250.0
      mw_aer_mac(iopcg5_b_o_a)=250.0
      mw_aer_mac(iopcg6_b_o_a)=250.0
      mw_aer_mac(iopcg7_b_o_a)=250.0
      mw_aer_mac(iopcg8_b_o_a)=250.0
      mw_aer_mac(ipcg1_f_c_a) =250.0
      mw_aer_mac(ipcg2_f_c_a) =250.0
      mw_aer_mac(ipcg3_f_c_a)=250.0
      mw_aer_mac(ipcg4_f_c_a)=250.0
      mw_aer_mac(ipcg5_f_c_a)=250.0
      mw_aer_mac(ipcg6_f_c_a)=250.0
      mw_aer_mac(ipcg7_f_c_a)=250.0
      mw_aer_mac(ipcg8_f_c_a)=250.0
      mw_aer_mac(ipcg9_f_c_a)=250.0
      mw_aer_mac(iopcg1_f_c_a)=250.0
      mw_aer_mac(iopcg2_f_c_a)=250.0
      mw_aer_mac(iopcg3_f_c_a)=250.0
      mw_aer_mac(iopcg4_f_c_a)=250.0
      mw_aer_mac(iopcg5_f_c_a)=250.0
      mw_aer_mac(iopcg6_f_c_a)=250.0
      mw_aer_mac(iopcg7_f_c_a)=250.0
      mw_aer_mac(iopcg8_f_c_a)=250.0
      mw_aer_mac(ipcg1_f_o_a)=250.0
      mw_aer_mac(ipcg2_f_o_a)=250.0
      mw_aer_mac(ipcg3_f_o_a)=250.0
      mw_aer_mac(ipcg4_f_o_a)=250.0
      mw_aer_mac(ipcg5_f_o_a)=250.0
      mw_aer_mac(ipcg6_f_o_a)=250.0
      mw_aer_mac(ipcg7_f_o_a)=250.0
      mw_aer_mac(ipcg8_f_o_a)=250.0
      mw_aer_mac(ipcg9_f_o_a)=250.0
      mw_aer_mac(iopcg1_f_o_a)=250.0
      mw_aer_mac(iopcg2_f_o_a)=250.0
      mw_aer_mac(iopcg3_f_o_a)=250.0
      mw_aer_mac(iopcg4_f_o_a)=250.0
      mw_aer_mac(iopcg5_f_o_a)=250.0
      mw_aer_mac(iopcg6_f_o_a)=250.0
      mw_aer_mac(iopcg7_f_o_a)=250.0
      mw_aer_mac(iopcg8_f_o_a)=250.0
      mw_aer_mac(ismpa_a) = 250.0
      mw_aer_mac(ismpbb_a) = 250.0
      mw_aer_mac(iglysoa_r1_a) = 250.0
      mw_aer_mac(iglysoa_r2_a) = 250.0
      mw_aer_mac(iglysoa_sfc_a) = 250.0
      mw_aer_mac(iglysoa_nh4_a) = 250.0
      mw_aer_mac(iglysoa_oh_a) = 250.0
      mw_aer_mac(iant1_c_a) = 250.0
      mw_aer_mac(iant2_c_a) = 250.0
      mw_aer_mac(iant3_c_a) = 250.0
      mw_aer_mac(iant4_c_a) = 250.0
      mw_aer_mac(iant1_o_a) = 250.0
      mw_aer_mac(iant2_o_a) = 250.0
      mw_aer_mac(iant3_o_a) = 250.0
      mw_aer_mac(iant4_o_a) = 250.0
      mw_aer_mac(ibiog1_c_a) = 250.0
      mw_aer_mac(ibiog2_c_a) = 250.0
      mw_aer_mac(ibiog3_c_a) = 250.0
      mw_aer_mac(ibiog4_c_a) = 250.0
      mw_aer_mac(ibiog1_o_a) = 250.0
      mw_aer_mac(ibiog2_o_a) = 250.0
      mw_aer_mac(ibiog3_o_a) = 250.0
      mw_aer_mac(ibiog4_o_a) = 250.0
      mw_aer_mac(iasoaX_a) = 250.0
      mw_aer_mac(iasoa1_a) = 250.0
      mw_aer_mac(iasoa2_a) = 250.0
      mw_aer_mac(iasoa3_a) = 250.0
      mw_aer_mac(iasoa4_a) = 250.0
      mw_aer_mac(ibsoaX_a) = 250.0
      mw_aer_mac(ibsoa1_a) = 250.0
      mw_aer_mac(ibsoa2_a) = 250.0
      mw_aer_mac(ibsoa3_a) = 250.0
      mw_aer_mac(ibsoa4_a) = 250.0



      mw_comp_a(jnh4so4) = 132.0
      mw_comp_a(jlvcite) = 247.0
      mw_comp_a(jnh4hso4)= 115.0
      mw_comp_a(jnh4msa) = 113.0
      mw_comp_a(jnh4no3) = 80.0
      mw_comp_a(jnh4cl)  = 53.5
      mw_comp_a(jnacl)   = 58.5
      mw_comp_a(jnano3)  = 85.0
      mw_comp_a(jna2so4) = 142.0
      mw_comp_a(jna3hso4)= 262.0
      mw_comp_a(jnahso4) = 120.0
      mw_comp_a(jnamsa)  = 118.0
      mw_comp_a(jcaso4)  = 136.0
      mw_comp_a(jcamsa2) = 230.0
      mw_comp_a(jcano3)  = 164.0
      mw_comp_a(jcacl2)  = 111.0
      mw_comp_a(jcaco3)  = 100.0
      mw_comp_a(jh2so4)  = 98.0
      mw_comp_a(jhhso4)  = 98.0
      mw_comp_a(jhno3)   = 63.0
      mw_comp_a(jhcl)    = 36.5
      mw_comp_a(jmsa)    = 96.0
      mw_comp_a(joc)	 = 250.0
      mw_comp_a(jbc)	 = 1.0
      mw_comp_a(join)    = 1.0
      mw_comp_a(jh2o)    = 18.0
      mw_comp_a(jpcg1_b_c) =250.0
      mw_comp_a(jpcg2_b_c) =250.0
      mw_comp_a(jpcg3_b_c)=250.0
      mw_comp_a(jpcg4_b_c)=250.0
      mw_comp_a(jpcg5_b_c)=250.0
      mw_comp_a(jpcg6_b_c)=250.0
      mw_comp_a(jpcg7_b_c)=250.0
      mw_comp_a(jpcg8_b_c)=250.0
      mw_comp_a(jpcg9_b_c)=250.0
      mw_comp_a(jopcg1_b_c)=250.0
      mw_comp_a(jopcg2_b_c)=250.0
      mw_comp_a(jopcg3_b_c)=250.0
      mw_comp_a(jopcg4_b_c)=250.0
      mw_comp_a(jopcg5_b_c)=250.0
      mw_comp_a(jopcg6_b_c)=250.0
      mw_comp_a(jopcg7_b_c)=250.0
      mw_comp_a(jopcg8_b_c)=250.0
      mw_comp_a(jpcg1_b_o)=250.0
      mw_comp_a(jpcg2_b_o)=250.0
      mw_comp_a(jpcg3_b_o)=250.0
      mw_comp_a(jpcg4_b_o)=250.0
      mw_comp_a(jpcg5_b_o)=250.0
      mw_comp_a(jpcg6_b_o)=250.0
      mw_comp_a(jpcg7_b_o)=250.0
      mw_comp_a(jpcg8_b_o)=250.0
      mw_comp_a(jpcg9_b_o)=250.0
      mw_comp_a(jopcg1_b_o)=250.0
      mw_comp_a(jopcg2_b_o)=250.0
      mw_comp_a(jopcg3_b_o)=250.0
      mw_comp_a(jopcg4_b_o)=250.0
      mw_comp_a(jopcg5_b_o)=250.0
      mw_comp_a(jopcg6_b_o)=250.0
      mw_comp_a(jopcg7_b_o)=250.0
      mw_comp_a(jopcg8_b_o)=250.0
      mw_comp_a(jpcg1_f_c) =250.0
      mw_comp_a(jpcg2_f_c) =250.0
      mw_comp_a(jpcg3_f_c)=250.0
      mw_comp_a(jpcg4_f_c)=250.0
      mw_comp_a(jpcg5_f_c)=250.0
      mw_comp_a(jpcg6_f_c)=250.0
      mw_comp_a(jpcg7_f_c)=250.0
      mw_comp_a(jpcg8_f_c)=250.0
      mw_comp_a(jpcg9_f_c)=250.0
      mw_comp_a(jopcg1_f_c)=250.0
      mw_comp_a(jopcg2_f_c)=250.0
      mw_comp_a(jopcg3_f_c)=250.0
      mw_comp_a(jopcg4_f_c)=250.0
      mw_comp_a(jopcg5_f_c)=250.0
      mw_comp_a(jopcg6_f_c)=250.0
      mw_comp_a(jopcg7_f_c)=250.0
      mw_comp_a(jopcg8_f_c)=250.0
      mw_comp_a(jpcg1_f_o)=250.0
      mw_comp_a(jpcg2_f_o)=250.0
      mw_comp_a(jpcg3_f_o)=250.0
      mw_comp_a(jpcg4_f_o)=250.0
      mw_comp_a(jpcg5_f_o)=250.0
      mw_comp_a(jpcg6_f_o)=250.0
      mw_comp_a(jpcg7_f_o)=250.0
      mw_comp_a(jpcg8_f_o)=250.0
      mw_comp_a(jpcg9_f_o)=250.0
      mw_comp_a(jopcg1_f_o)=250.0
      mw_comp_a(jopcg2_f_o)=250.0
      mw_comp_a(jopcg3_f_o)=250.0
      mw_comp_a(jopcg4_f_o)=250.0
      mw_comp_a(jopcg5_f_o)=250.0
      mw_comp_a(jopcg6_f_o)=250.0
      mw_comp_a(jopcg7_f_o)=250.0
      mw_comp_a(jopcg8_f_o)=250.0
      mw_comp_a(jsmpa)=250.0
      mw_comp_a(jsmpbb)=250.0
      mw_comp_a(jglysoa_r1)=250.0
      mw_comp_a(jglysoa_r2)=250.0
      mw_comp_a(jglysoa_sfc)=250.0
      mw_comp_a(jglysoa_nh4)=250.0
      mw_comp_a(jglysoa_oh)=250.0
      mw_comp_a(jant1_c)=250.0
      mw_comp_a(jant2_c)=250.0
      mw_comp_a(jant3_c)=250.0
      mw_comp_a(jant4_c)=250.0
      mw_comp_a(jant1_o)=250.0
      mw_comp_a(jant2_o)=250.0
      mw_comp_a(jant3_o)=250.0
      mw_comp_a(jant4_o)=250.0
      mw_comp_a(jbiog1_c)=250.0
      mw_comp_a(jbiog2_c)=250.0
      mw_comp_a(jbiog3_c)=250.0
      mw_comp_a(jbiog4_c)=250.0
      mw_comp_a(jbiog1_o)=250.0
      mw_comp_a(jbiog2_o)=250.0
      mw_comp_a(jbiog3_o)=250.0
      mw_comp_a(jbiog4_o)=250.0
      mw_comp_a(jasoaX)=250.0
      mw_comp_a(jasoa1)=250.0
      mw_comp_a(jasoa2)=250.0
      mw_comp_a(jasoa3)=250.0
      mw_comp_a(jasoa4)=250.0
      mw_comp_a(jbsoaX)=250.0
      mw_comp_a(jbsoa1)=250.0
      mw_comp_a(jbsoa2)=250.0
      mw_comp_a(jbsoa3)=250.0
      mw_comp_a(jbsoa4)=250.0


      dens_aer_mac(iso4_a) = 1.8	
      dens_aer_mac(ino3_a) = 1.8	
      dens_aer_mac(icl_a)  = 2.2	
      dens_aer_mac(imsa_a) = 1.8	
      dens_aer_mac(ico3_a) = 2.6	
      dens_aer_mac(inh4_a) = 1.8	
      dens_aer_mac(ina_a)  = 2.2	
      dens_aer_mac(ica_a)  = 2.6	
      dens_aer_mac(ioin_a) = 2.6	
      dens_aer_mac(ioc_a)  = 1.0	
      dens_aer_mac(ibc_a)  = 1.7	
      dens_aer_mac(ipcg1_b_c_a) =1.0
      dens_aer_mac(ipcg2_b_c_a) =1.0
      dens_aer_mac(ipcg3_b_c_a)=1.0
      dens_aer_mac(ipcg4_b_c_a)=1.0
      dens_aer_mac(ipcg5_b_c_a)=1.0
      dens_aer_mac(ipcg6_b_c_a)=1.0
      dens_aer_mac(ipcg7_b_c_a)=1.0
      dens_aer_mac(ipcg8_b_c_a)=1.0
      dens_aer_mac(ipcg9_b_c_a)=1.0
      dens_aer_mac(iopcg1_b_c_a)=1.0
      dens_aer_mac(iopcg2_b_c_a)=1.0
      dens_aer_mac(iopcg3_b_c_a)=1.0
      dens_aer_mac(iopcg4_b_c_a)=1.0
      dens_aer_mac(iopcg5_b_c_a)=1.0
      dens_aer_mac(iopcg6_b_c_a)=1.0
      dens_aer_mac(iopcg7_b_c_a)=1.0
      dens_aer_mac(iopcg8_b_c_a)=1.0
      dens_aer_mac(ipcg1_b_o_a)=1.0
      dens_aer_mac(ipcg2_b_o_a)=1.0
      dens_aer_mac(ipcg3_b_o_a)=1.0
      dens_aer_mac(ipcg4_b_o_a)=1.0
      dens_aer_mac(ipcg5_b_o_a)=1.0
      dens_aer_mac(ipcg6_b_o_a)=1.0
      dens_aer_mac(ipcg7_b_o_a)=1.0
      dens_aer_mac(ipcg8_b_o_a)=1.0
      dens_aer_mac(ipcg9_b_o_a)=1.0
      dens_aer_mac(iopcg1_b_o_a)=1.0
      dens_aer_mac(iopcg2_b_o_a)=1.0
      dens_aer_mac(iopcg3_b_o_a)=1.0
      dens_aer_mac(iopcg4_b_o_a)=1.0
      dens_aer_mac(iopcg5_b_o_a)=1.0
      dens_aer_mac(iopcg6_b_o_a)=1.0
      dens_aer_mac(iopcg7_b_o_a)=1.0
      dens_aer_mac(iopcg8_b_o_a)=1.0
      dens_aer_mac(ipcg1_f_c_a) =1.0
      dens_aer_mac(ipcg2_f_c_a) =1.0
      dens_aer_mac(ipcg3_f_c_a)=1.0
      dens_aer_mac(ipcg4_f_c_a)=1.0
      dens_aer_mac(ipcg5_f_c_a)=1.0
      dens_aer_mac(ipcg6_f_c_a)=1.0
      dens_aer_mac(ipcg7_f_c_a)=1.0
      dens_aer_mac(ipcg8_f_c_a)=1.0
      dens_aer_mac(ipcg9_f_c_a)=1.0
      dens_aer_mac(iopcg1_f_c_a)=1.0
      dens_aer_mac(iopcg2_f_c_a)=1.0
      dens_aer_mac(iopcg3_f_c_a)=1.0
      dens_aer_mac(iopcg4_f_c_a)=1.0
      dens_aer_mac(iopcg5_f_c_a)=1.0
      dens_aer_mac(iopcg6_f_c_a)=1.0
      dens_aer_mac(iopcg7_f_c_a)=1.0
      dens_aer_mac(iopcg8_f_c_a)=1.0
      dens_aer_mac(ipcg1_f_o_a)=1.0
      dens_aer_mac(ipcg2_f_o_a)=1.0
      dens_aer_mac(ipcg3_f_o_a)=1.0
      dens_aer_mac(ipcg4_f_o_a)=1.0
      dens_aer_mac(ipcg5_f_o_a)=1.0
      dens_aer_mac(ipcg6_f_o_a)=1.0
      dens_aer_mac(ipcg7_f_o_a)=1.0
      dens_aer_mac(ipcg8_f_o_a)=1.0
      dens_aer_mac(ipcg9_f_o_a)=1.0
      dens_aer_mac(iopcg1_f_o_a)=1.0
      dens_aer_mac(iopcg2_f_o_a)=1.0
      dens_aer_mac(iopcg3_f_o_a)=1.0
      dens_aer_mac(iopcg4_f_o_a)=1.0
      dens_aer_mac(iopcg5_f_o_a)=1.0
      dens_aer_mac(iopcg6_f_o_a)=1.0
      dens_aer_mac(iopcg7_f_o_a)=1.0
      dens_aer_mac(iopcg8_f_o_a)=1.0
      dens_aer_mac(ismpa_a)=1.0
      dens_aer_mac(ismpbb_a)=1.0
      dens_aer_mac(iglysoa_r1_a)=1.0
      dens_aer_mac(iglysoa_r2_a)=1.0
      dens_aer_mac(iglysoa_sfc_a)=1.0
      dens_aer_mac(iglysoa_nh4_a)=1.0
      dens_aer_mac(iglysoa_oh_a)=1.0
      dens_aer_mac(iant1_c_a)=1.0
      dens_aer_mac(iant2_c_a)=1.0
      dens_aer_mac(iant3_c_a)=1.0
      dens_aer_mac(iant4_c_a)=1.0
      dens_aer_mac(iant1_o_a)=1.0
      dens_aer_mac(iant2_o_a)=1.0
      dens_aer_mac(iant3_o_a)=1.0
      dens_aer_mac(iant4_o_a)=1.0
      dens_aer_mac(ibiog1_c_a)=1.0
      dens_aer_mac(ibiog2_c_a)=1.0
      dens_aer_mac(ibiog3_c_a)=1.0
      dens_aer_mac(ibiog4_c_a)=1.0
      dens_aer_mac(ibiog1_o_a)=1.0
      dens_aer_mac(ibiog2_o_a)=1.0
      dens_aer_mac(ibiog3_o_a)=1.0
      dens_aer_mac(ibiog4_o_a)=1.0
      dens_aer_mac(iasoaX_a)=1.5
      dens_aer_mac(iasoa1_a)=1.5
      dens_aer_mac(iasoa2_a)=1.5
      dens_aer_mac(iasoa3_a)=1.5
      dens_aer_mac(iasoa4_a)=1.5
      dens_aer_mac(ibsoaX_a)=1.5
      dens_aer_mac(ibsoa1_a)=1.5
      dens_aer_mac(ibsoa2_a)=1.5
      dens_aer_mac(ibsoa3_a)=1.5
      dens_aer_mac(ibsoa4_a)=1.5


      partial_molar_vol(ih2so4_g) = 51.83
      partial_molar_vol(ihno3_g)  = 31.45
      partial_molar_vol(ihcl_g)   = 20.96
      partial_molar_vol(inh3_g)   = 24.03
      partial_molar_vol(imsa_g)   = 53.33
      partial_molar_vol(ipcg1_b_c_g) =250.0
      partial_molar_vol(ipcg2_b_c_g) =250.0
      partial_molar_vol(ipcg3_b_c_g)=250.0
      partial_molar_vol(ipcg4_b_c_g)=250.0
      partial_molar_vol(ipcg5_b_c_g)=250.0
      partial_molar_vol(ipcg6_b_c_g)=250.0
      partial_molar_vol(ipcg7_b_c_g)=250.0
      partial_molar_vol(ipcg8_b_c_g)=250.0
      partial_molar_vol(ipcg9_b_c_g)=250.0
      partial_molar_vol(iopcg1_b_c_g)=250.0
      partial_molar_vol(iopcg2_b_c_g)=250.0
      partial_molar_vol(iopcg3_b_c_g)=250.0
      partial_molar_vol(iopcg4_b_c_g)=250.0
      partial_molar_vol(iopcg5_b_c_g)=250.0
      partial_molar_vol(iopcg6_b_c_g)=250.0
      partial_molar_vol(iopcg7_b_c_g)=250.0
      partial_molar_vol(iopcg8_b_c_g)=250.0
      partial_molar_vol(ipcg1_b_o_g)=250.0
      partial_molar_vol(ipcg2_b_o_g)=250.0
      partial_molar_vol(ipcg3_b_o_g)=250.0
      partial_molar_vol(ipcg4_b_o_g)=250.0
      partial_molar_vol(ipcg5_b_o_g)=250.0
      partial_molar_vol(ipcg6_b_o_g)=250.0
      partial_molar_vol(ipcg7_b_o_g)=250.0
      partial_molar_vol(ipcg8_b_o_g)=250.0
      partial_molar_vol(ipcg9_b_o_g)=250.0
      partial_molar_vol(iopcg1_b_o_g)=250.0
      partial_molar_vol(iopcg2_b_o_g)=250.0
      partial_molar_vol(iopcg3_b_o_g)=250.0
      partial_molar_vol(iopcg4_b_o_g)=250.0
      partial_molar_vol(iopcg5_b_o_g)=250.0
      partial_molar_vol(iopcg6_b_o_g)=250.0
      partial_molar_vol(iopcg7_b_o_g)=250.0
      partial_molar_vol(iopcg8_b_o_g)=250.0
      partial_molar_vol(ipcg1_f_c_g) =250.0
      partial_molar_vol(ipcg2_f_c_g) =250.0
      partial_molar_vol(ipcg3_f_c_g)=250.0
      partial_molar_vol(ipcg4_f_c_g)=250.0
      partial_molar_vol(ipcg5_f_c_g)=250.0
      partial_molar_vol(ipcg6_f_c_g)=250.0
      partial_molar_vol(ipcg7_f_c_g)=250.0
      partial_molar_vol(ipcg8_f_c_g)=250.0
      partial_molar_vol(ipcg9_f_c_g)=250.0
      partial_molar_vol(iopcg1_f_c_g)=250.0
      partial_molar_vol(iopcg2_f_c_g)=250.0
      partial_molar_vol(iopcg3_f_c_g)=250.0
      partial_molar_vol(iopcg4_f_c_g)=250.0
      partial_molar_vol(iopcg5_f_c_g)=250.0
      partial_molar_vol(iopcg6_f_c_g)=250.0
      partial_molar_vol(iopcg7_f_c_g)=250.0
      partial_molar_vol(iopcg8_f_c_g)=250.0
      partial_molar_vol(ipcg1_f_o_g)=250.0
      partial_molar_vol(ipcg2_f_o_g)=250.0
      partial_molar_vol(ipcg3_f_o_g)=250.0
      partial_molar_vol(ipcg4_f_o_g)=250.0
      partial_molar_vol(ipcg5_f_o_g)=250.0
      partial_molar_vol(ipcg6_f_o_g)=250.0
      partial_molar_vol(ipcg7_f_o_g)=250.0
      partial_molar_vol(ipcg8_f_o_g)=250.0
      partial_molar_vol(ipcg9_f_o_g)=250.0
      partial_molar_vol(iopcg1_f_o_g)=250.0
      partial_molar_vol(iopcg2_f_o_g)=250.0
      partial_molar_vol(iopcg3_f_o_g)=250.0
      partial_molar_vol(iopcg4_f_o_g)=250.0
      partial_molar_vol(iopcg5_f_o_g)=250.0
      partial_molar_vol(iopcg6_f_o_g)=250.0
      partial_molar_vol(iopcg7_f_o_g)=250.0
      partial_molar_vol(iopcg8_f_o_g)=250.0
      partial_molar_vol(ismpa_g)=250.0
      partial_molar_vol(ismpbb_g)=250.0
      partial_molar_vol(iant1_c_g)=250.0
      partial_molar_vol(iant2_c_g)=250.0
      partial_molar_vol(iant3_c_g)=250.0
      partial_molar_vol(iant4_c_g)=250.0
      partial_molar_vol(iant1_o_g)=250.0
      partial_molar_vol(iant2_o_g)=250.0
      partial_molar_vol(iant3_o_g)=250.0
      partial_molar_vol(iant4_o_g)=250.0
      partial_molar_vol(ibiog1_c_g)=250.0
      partial_molar_vol(ibiog2_c_g)=250.0
      partial_molar_vol(ibiog3_c_g)=250.0
      partial_molar_vol(ibiog4_c_g)=250.0
      partial_molar_vol(ibiog1_o_g)=250.0
      partial_molar_vol(ibiog2_o_g)=250.0
      partial_molar_vol(ibiog3_o_g)=250.0
      partial_molar_vol(ibiog4_o_g)=250.0
      partial_molar_vol(in2o5_g)  = 200.0	
      partial_molar_vol(iclno2_g) = 200.0	
      partial_molar_vol(iasoaX_g)=250.0
      partial_molar_vol(iasoa1_g)=250.0
      partial_molar_vol(iasoa2_g)=250.0
      partial_molar_vol(iasoa3_g)=250.0
      partial_molar_vol(iasoa4_g)=250.0
      partial_molar_vol(ibsoaX_g)=250.0
      partial_molar_vol(ibsoa1_g)=250.0
      partial_molar_vol(ibsoa2_g)=250.0
      partial_molar_vol(ibsoa3_g)=250.0
      partial_molar_vol(ibsoa4_g)=250.0
      partial_molar_vol(igly)=58.0
      partial_molar_vol(iho)=17.0


      ref_index_a(jnh4so4) = cmplx(1.52,0.)
      ref_index_a(jlvcite) = cmplx(1.50,0.)
      ref_index_a(jnh4hso4)= cmplx(1.47,0.)
      ref_index_a(jnh4msa) = cmplx(1.50,0.)	
      ref_index_a(jnh4no3) = cmplx(1.50,0.)
      ref_index_a(jnh4cl)  = cmplx(1.50,0.)
      ref_index_a(jnacl)   = cmplx(1.45,0.)
      ref_index_a(jnano3)  = cmplx(1.50,0.)
      ref_index_a(jna2so4) = cmplx(1.50,0.)
      ref_index_a(jna3hso4)= cmplx(1.50,0.)
      ref_index_a(jnahso4) = cmplx(1.50,0.)
      ref_index_a(jnamsa)  = cmplx(1.50,0.)	
      ref_index_a(jcaso4)  = cmplx(1.56,0.006)
      ref_index_a(jcamsa2) = cmplx(1.56,0.006)	
      ref_index_a(jcano3)  = cmplx(1.56,0.006)
      ref_index_a(jcacl2)  = cmplx(1.52,0.006)
      ref_index_a(jcaco3)  = cmplx(1.68,0.006)
      ref_index_a(jh2so4)  = cmplx(1.43,0.)
      ref_index_a(jhhso4)  = cmplx(1.43,0.)
      ref_index_a(jhno3)   = cmplx(1.50,0.)
      ref_index_a(jhcl)    = cmplx(1.50,0.)
      ref_index_a(jmsa)    = cmplx(1.43,0.)	
      ref_index_a(joc)	   = cmplx(1.45,0.)
      ref_index_a(jbc)	   = cmplx(1.82,0.74)
      ref_index_a(join)    = cmplx(1.55,0.006)
      ref_index_a(jh2o)    = cmplx(1.33,0.)


      jsalt_index(jnh4so4) = 5		
      jsalt_index(jlvcite) = 2		
      jsalt_index(jnh4hso4)= 1		
      jsalt_index(jnh4no3) = 2		
      jsalt_index(jnh4cl)  = 1		
      jsalt_index(jna2so4) = 60		
      jsalt_index(jnahso4) = 10		
      jsalt_index(jnano3)  = 40		
      jsalt_index(jnacl)   = 10		
      jsalt_index(jcano3)  = 120	
      jsalt_index(jcacl2)  = 80		
      jsalt_index(jnh4msa) = 0		
      jsalt_index(jnamsa)  = 0		
      jsalt_index(jcamsa2) = 0		







      jsulf_poor(1)   = 	1	
      jsulf_poor(2)   = 	2	
      jsulf_poor(5)   = 	3	
      jsulf_poor(10)  = 	4	
      jsulf_poor(40)  = 	5	
      jsulf_poor(60)  = 	6	
      jsulf_poor(80)  = 	7	
      jsulf_poor(120) = 	8	
      jsulf_poor(3)   = 	9	
      jsulf_poor(6)   = 	10	
      jsulf_poor(7)   = 	11	
      jsulf_poor(8)   =  	12	
      jsulf_poor(11)  = 	13	
      jsulf_poor(41)  = 	14	
      jsulf_poor(42)  = 	15	
      jsulf_poor(43)  = 	16	
      jsulf_poor(50)  = 	17	
      jsulf_poor(51)  = 	18	
      jsulf_poor(61)  = 	19	
      jsulf_poor(62)  = 	20	
      jsulf_poor(63)  = 	21	
      jsulf_poor(65)  = 	22	
      jsulf_poor(66)  = 	23	
      jsulf_poor(67)  = 	24	
      jsulf_poor(68)  = 	25	
      jsulf_poor(70)  = 	26	
      jsulf_poor(71)  = 	27	
      jsulf_poor(100) = 	28	
      jsulf_poor(101) = 	29	
      jsulf_poor(102) = 	30	
      jsulf_poor(103) = 	31	
      jsulf_poor(110) = 	32	
      jsulf_poor(111) = 	33	
      jsulf_poor(81)  = 	34	
      jsulf_poor(90)  = 	35	
      jsulf_poor(91)  = 	36	
      jsulf_poor(121) = 	37	
      jsulf_poor(122) = 	38	
      jsulf_poor(123) = 	39	
      jsulf_poor(130) = 	40	
      jsulf_poor(131) = 	41	
      jsulf_poor(160) = 	42	
      jsulf_poor(161) = 	43	
      jsulf_poor(162) = 	44	
      jsulf_poor(163) = 	45	
      jsulf_poor(170) = 	46	
      jsulf_poor(171) = 	47	
      jsulf_poor(200) = 	48	
      jsulf_poor(201) = 	49	
      jsulf_poor(210) = 	50	
      jsulf_poor(211) = 	51	


      jsulf_rich(1)   = 	52	
      jsulf_rich(2)   = 	53	
      jsulf_rich(10)  = 	54	
      jsulf_rich(3)   = 	55	
      jsulf_rich(7)   = 	56	
      jsulf_rich(70)  = 	57	
      jsulf_rich(62)  = 	58	
      jsulf_rich(67)  = 	59	
      jsulf_rich(61)  = 	60	
      jsulf_rich(63)  = 	61	
      jsulf_rich(11)  = 	62	
      jsulf_rich(71)  = 	63	
      jsulf_rich(5)   = 	3	
      jsulf_rich(60)  = 	6	
      jsulf_rich(65)  = 	22	










      je = jnh4so4
      a_zsr(1,je)  =  1.30894
      a_zsr(2,je)  = -7.09922
      a_zsr(3,je)  =  20.62831
      a_zsr(4,je)  = -32.19965
      a_zsr(5,je)  =  25.17026
      a_zsr(6,je)  = -7.81632
      aw_min(je)   = 0.1


      je = jlvcite
      a_zsr(1,je)  =  1.10725
      a_zsr(2,je)  = -5.17978
      a_zsr(3,je)  =  12.29534
      a_zsr(4,je)  = -16.32545
      a_zsr(5,je)  =  11.29274
      a_zsr(6,je)  = -3.19164
      aw_min(je)   = 0.1


      je = jnh4hso4
      a_zsr(1,je)  =  1.15510
      a_zsr(2,je)  = -3.20815
      a_zsr(3,je)  =  2.71141
      a_zsr(4,je)  =  2.01155
      a_zsr(5,je)  = -4.71014
      a_zsr(6,je)  =  2.04616
      aw_min(je)   = 0.1


      je = jnh4msa
      a_zsr(1,je)  =  1.15510
      a_zsr(2,je)  = -3.20815
      a_zsr(3,je)  =  2.71141
      a_zsr(4,je)  =  2.01155
      a_zsr(5,je)  = -4.71014
      a_zsr(6,je)  =  2.04616
      aw_min(je)   = 0.1


      je = jnh4no3
      a_zsr(1,je)  =  0.43507
      a_zsr(2,je)  =  6.38220
      a_zsr(3,je)  = -30.19797
      a_zsr(4,je)  =  53.36470
      a_zsr(5,je)  = -43.44203
      a_zsr(6,je)  =  13.46158
      aw_min(je)   = 0.1


      je = jnh4cl
      a_zsr(1,je)  =  0.45309
      a_zsr(2,je)  =  2.65606
      a_zsr(3,je)  = -14.7730
      a_zsr(4,je)  =  26.2936
      a_zsr(5,je)  = -20.5735
      a_zsr(6,je)  =  5.94255
      aw_min(je)   = 0.1


      je = jnacl
      a_zsr(1,je)  =  0.42922
      a_zsr(2,je)  = -1.17718
      a_zsr(3,je)  =  2.80208
      a_zsr(4,je)  = -4.51097
      a_zsr(5,je)  =  3.76963
      a_zsr(6,je)  = -1.31359
      aw_min(je)   = 0.1


      je = jnano3
      a_zsr(1,je)  =  1.34966
      a_zsr(2,je)  = -5.20116
      a_zsr(3,je)  =  11.49011
      a_zsr(4,je)  = -14.41380
      a_zsr(5,je)  =  9.07037
      a_zsr(6,je)  = -2.29769
      aw_min(je)   = 0.1


      je = jna2so4
      a_zsr(1,je)  =  0.39888
      a_zsr(2,je)  = -1.27150
      a_zsr(3,je)  =  3.42792
      a_zsr(4,je)  = -5.92632
      a_zsr(5,je)  =  5.33351
      a_zsr(6,je)  = -1.96541
      aw_min(je)   = 0.1


      je = jna3hso4
      a_zsr(1,je)  =  0.31480
      a_zsr(2,je)  = -1.01087
      a_zsr(3,je)  =  2.44029
      a_zsr(4,je)  = -3.66095
      a_zsr(5,je)  =  2.77632
      a_zsr(6,je)  = -0.86058
      aw_min(je)   = 0.1


      je = jnahso4
      a_zsr(1,je)  =  0.62764
      a_zsr(2,je)  = -1.63520
      a_zsr(3,je)  =  4.62531
      a_zsr(4,je)  = -10.06925
      a_zsr(5,je)  =  10.33547
      a_zsr(6,je)  = -3.88729
      aw_min(je)   = 0.1


      je = jnamsa
      a_zsr(1,je)  =  0.62764
      a_zsr(2,je)  = -1.63520
      a_zsr(3,je)  =  4.62531
      a_zsr(4,je)  = -10.06925
      a_zsr(5,je)  =  10.33547
      a_zsr(6,je)  = -3.88729
      aw_min(je)   = 0.1


      je = jcano3
      a_zsr(1,je)  =  0.38895
      a_zsr(2,je)  = -1.16013
      a_zsr(3,je)  =  2.16819
      a_zsr(4,je)  = -2.23079
      a_zsr(5,je)  =  1.00268
      a_zsr(6,je)  = -0.16923
      aw_min(je)   = 0.1


      je = jcacl2
      a_zsr(1,je)  =  0.29891
      a_zsr(2,je)  = -1.31104
      a_zsr(3,je)  =  3.68759
      a_zsr(4,je)  = -5.81708
      a_zsr(5,je)  =  4.67520
      a_zsr(6,je)  = -1.53223
      aw_min(je)   = 0.1


      je = jh2so4
      a_zsr(1,je) =  0.32751
      a_zsr(2,je) = -1.00692
      a_zsr(3,je) =  2.59750
      a_zsr(4,je) = -4.40014
      a_zsr(5,je) =  3.88212
      a_zsr(6,je) = -1.39916
      aw_min(je)  = 0.1


      je = jmsa
      a_zsr(1,je) =  0.32751
      a_zsr(2,je) = -1.00692
      a_zsr(3,je) =  2.59750
      a_zsr(4,je) = -4.40014
      a_zsr(5,je) =  3.88212
      a_zsr(6,je) = -1.39916
      aw_min(je)  = 0.1


      je = jhhso4
      a_zsr(1,je) =  0.32751
      a_zsr(2,je) = -1.00692
      a_zsr(3,je) =  2.59750
      a_zsr(4,je) = -4.40014
      a_zsr(5,je) =  3.88212
      a_zsr(6,je) = -1.39916
      aw_min(je)  = 1.0


      je = jhno3
      a_zsr(1,je) =  0.75876
      a_zsr(2,je) = -3.31529
      a_zsr(3,je) =  9.26392
      a_zsr(4,je) = -14.89799
      a_zsr(5,je) =  12.08781
      a_zsr(6,je) = -3.89958
      aw_min(je)  = 0.1


      je = jhcl
      a_zsr(1,je) =  0.31133
      a_zsr(2,je) = -0.79688
      a_zsr(3,je) =  1.93995
      a_zsr(4,je) = -3.31582
      a_zsr(5,je) =  2.93513
      a_zsr(6,je) = -1.07268
      aw_min(je)  = 0.1


      je = jcaso4
      a_zsr(1,je)  =  0.0
      a_zsr(2,je)  =  0.0
      a_zsr(3,je)  =  0.0
      a_zsr(4,je)  =  0.0
      a_zsr(5,je)  =  0.0
      a_zsr(6,je)  =  0.0
      aw_min(je)   = 1.0


      je = jcamsa2
      a_zsr(1,je)  =  0.38895
      a_zsr(2,je)  = -1.16013
      a_zsr(3,je)  =  2.16819
      a_zsr(4,je)  = -2.23079
      a_zsr(5,je)  =  1.00268
      a_zsr(6,je)  = -0.16923
      aw_min(je)   = 0.1


      je = jcaco3
      a_zsr(1,je)  =  0.0
      a_zsr(2,je)  =  0.0
      a_zsr(3,je)  =  0.0
      a_zsr(4,je)  =  0.0
      a_zsr(5,je)  =  0.0
      a_zsr(6,je)  =  0.0
      aw_min(je)   = 1.0







      b_zsr(jnh4so4)  = 28.0811


      b_zsr(jlvcite)  = 14.7178


      b_zsr(jnh4hso4) = 29.4779


      b_zsr(jnh4msa)  = 29.4779 


      b_zsr(jnh4no3)  = 33.4049


      b_zsr(jnh4cl)   = 30.8888


      b_zsr(jnacl)    = 29.8375


      b_zsr(jnano3)   = 32.2756


      b_zsr(jna2so4)  = 27.6889


      b_zsr(jna3hso4) = 14.2184


      b_zsr(jnahso4)  = 28.3367


      b_zsr(jnamsa)   = 28.3367 


      b_zsr(jcano3)   = 18.3661


      b_zsr(jcacl2)   = 20.8792


      b_zsr(jh2so4)   = 26.7347


      b_zsr(jhhso4)   = 26.7347


      b_zsr(jhno3)    = 28.8257


      b_zsr(jhcl)     = 27.7108


      b_zsr(jmsa)     = 26.7347 


      b_zsr(jcaso4)   = 0.0


      b_zsr(jcamsa2)  = 18.3661 


      b_zsr(jcaco3)   = 0.0













      ja = jnh4so4


      je = jnh4so4
      b_mtem(1,ja,je) = -2.94685
      b_mtem(2,ja,je) = 17.3328
      b_mtem(3,ja,je) = -64.8441
      b_mtem(4,ja,je) = 122.7070
      b_mtem(5,ja,je) = -114.4373
      b_mtem(6,ja,je) = 41.6811


      je = jnh4no3
      b_mtem(1,ja,je) = -2.7503
      b_mtem(2,ja,je) = 4.3806
      b_mtem(3,ja,je) = -1.1110
      b_mtem(4,ja,je) = -1.7005
      b_mtem(5,ja,je) = -4.4207
      b_mtem(6,ja,je) = 5.1990


      je = jnh4cl
      b_mtem(1,ja,je) = -2.06952
      b_mtem(2,ja,je) = 7.1240
      b_mtem(3,ja,je) = -24.4274
      b_mtem(4,ja,je) = 51.1458
      b_mtem(5,ja,je) = -54.2056
      b_mtem(6,ja,je) = 22.0606


      je = jna2so4
      b_mtem(1,ja,je) = -2.17361
      b_mtem(2,ja,je) = 15.9919
      b_mtem(3,ja,je) = -69.0952
      b_mtem(4,ja,je) = 139.8860
      b_mtem(5,ja,je) = -134.9890
      b_mtem(6,ja,je) = 49.8877


      je = jnano3
      b_mtem(1,ja,je) = -4.4370
      b_mtem(2,ja,je) = 24.0243
      b_mtem(3,ja,je) = -76.2437
      b_mtem(4,ja,je) = 128.6660
      b_mtem(5,ja,je) = -110.0900
      b_mtem(6,ja,je) = 37.7414


      je = jnacl
      b_mtem(1,ja,je) = -1.5394
      b_mtem(2,ja,je) = 5.8671
      b_mtem(3,ja,je) = -22.7726
      b_mtem(4,ja,je) = 47.0547
      b_mtem(5,ja,je) = -47.8266
      b_mtem(6,ja,je) = 18.8489


      je = jhno3
      b_mtem(1,ja,je) = -0.35750
      b_mtem(2,ja,je) = -3.82466
      b_mtem(3,ja,je) = 4.55462
      b_mtem(4,ja,je) = 5.05402
      b_mtem(5,ja,je) = -14.7476
      b_mtem(6,ja,je) = 8.8009


      je = jhcl
      b_mtem(1,ja,je) = -2.15146
      b_mtem(2,ja,je) = 5.50205
      b_mtem(3,ja,je) = -19.1476
      b_mtem(4,ja,je) = 39.1880
      b_mtem(5,ja,je) = -39.9460
      b_mtem(6,ja,je) = 16.0700


      je = jh2so4
      b_mtem(1,ja,je) = -2.52604
      b_mtem(2,ja,je) = 9.76022
      b_mtem(3,ja,je) = -35.2540
      b_mtem(4,ja,je) = 71.2981
      b_mtem(5,ja,je) = -71.8207
      b_mtem(6,ja,je) = 28.0758



      je = jnh4hso4
      b_mtem(1,ja,je) = -4.13219
      b_mtem(2,ja,je) = 13.8863
      b_mtem(3,ja,je) = -34.5387
      b_mtem(4,ja,je) = 56.5012
      b_mtem(5,ja,je) = -51.8702
      b_mtem(6,ja,je) = 19.6232



      je = jlvcite
      b_mtem(1,ja,je) = -2.53482
      b_mtem(2,ja,je) = 12.3333
      b_mtem(3,ja,je) = -46.1020
      b_mtem(4,ja,je) = 90.4775
      b_mtem(5,ja,je) = -88.1254
      b_mtem(6,ja,je) = 33.4715



      je = jnahso4
      b_mtem(1,ja,je) = -3.23425
      b_mtem(2,ja,je) = 18.7842
      b_mtem(3,ja,je) = -78.7807
      b_mtem(4,ja,je) = 161.517
      b_mtem(5,ja,je) = -154.940
      b_mtem(6,ja,je) = 56.2252



      je = jna3hso4
      b_mtem(1,ja,je) = -1.25316
      b_mtem(2,ja,je) = 7.40960
      b_mtem(3,ja,je) = -34.8929
      b_mtem(4,ja,je) = 72.8853
      b_mtem(5,ja,je) = -72.4503
      b_mtem(6,ja,je) = 27.7706




      ja = jnh4no3


      je = jnh4so4
      b_mtem(1,ja,je) = -3.5201
      b_mtem(2,ja,je) = 21.6584
      b_mtem(3,ja,je) = -72.1499
      b_mtem(4,ja,je) = 126.7000
      b_mtem(5,ja,je) = -111.4550
      b_mtem(6,ja,je) = 38.5677


      je = jnh4no3
      b_mtem(1,ja,je) = -2.2630
      b_mtem(2,ja,je) = -0.1518
      b_mtem(3,ja,je) = 17.0898
      b_mtem(4,ja,je) = -36.7832
      b_mtem(5,ja,je) = 29.8407
      b_mtem(6,ja,je) = -7.9314


      je = jnh4cl
      b_mtem(1,ja,je) = -1.3851
      b_mtem(2,ja,je) = -0.4462
      b_mtem(3,ja,je) = 8.4567
      b_mtem(4,ja,je) = -11.5988
      b_mtem(5,ja,je) = 2.9802
      b_mtem(6,ja,je) = 1.8132


      je = jna2so4
      b_mtem(1,ja,je) = -1.7602
      b_mtem(2,ja,je) = 10.4044
      b_mtem(3,ja,je) = -35.5894
      b_mtem(4,ja,je) = 64.3584
      b_mtem(5,ja,je) = -57.8931
      b_mtem(6,ja,je) = 20.2141


      je = jnano3
      b_mtem(1,ja,je) = -3.24346
      b_mtem(2,ja,je) = 16.2794
      b_mtem(3,ja,je) = -48.7601
      b_mtem(4,ja,je) = 79.2246
      b_mtem(5,ja,je) = -65.8169
      b_mtem(6,ja,je) = 22.1500


      je = jnacl
      b_mtem(1,ja,je) = -1.75658
      b_mtem(2,ja,je) = 7.71384
      b_mtem(3,ja,je) = -22.7984
      b_mtem(4,ja,je) = 39.1532
      b_mtem(5,ja,je) = -34.6165
      b_mtem(6,ja,je) = 12.1283


      je = jcano3
      b_mtem(1,ja,je) = -0.97178
      b_mtem(2,ja,je) = 6.61964
      b_mtem(3,ja,je) = -26.2353
      b_mtem(4,ja,je) = 50.5259
      b_mtem(5,ja,je) = -47.6586
      b_mtem(6,ja,je) = 17.5074


      je = jcacl2
      b_mtem(1,ja,je) = -0.41515
      b_mtem(2,ja,je) = 6.44101
      b_mtem(3,ja,je) = -26.4473
      b_mtem(4,ja,je) = 49.0718
      b_mtem(5,ja,je) = -44.2631
      b_mtem(6,ja,je) = 15.3771


      je = jhno3
      b_mtem(1,ja,je) = -1.20644
      b_mtem(2,ja,je) = 5.70117
      b_mtem(3,ja,je) = -18.2783
      b_mtem(4,ja,je) = 31.7199
      b_mtem(5,ja,je) = -27.8703
      b_mtem(6,ja,je) = 9.7299


      je = jhcl
      b_mtem(1,ja,je) = -0.680862
      b_mtem(2,ja,je) = 3.59456
      b_mtem(3,ja,je) = -10.7969
      b_mtem(4,ja,je) = 17.8434
      b_mtem(5,ja,je) = -15.3165
      b_mtem(6,ja,je) = 5.17123




      ja = jnh4cl


      je = jnh4so4
      b_mtem(1,ja,je) = -2.8850
      b_mtem(2,ja,je) = 20.6970
      b_mtem(3,ja,je) = -70.6810
      b_mtem(4,ja,je) = 124.3690
      b_mtem(5,ja,je) = -109.2880
      b_mtem(6,ja,je) = 37.5831


      je = jnh4no3
      b_mtem(1,ja,je) = -1.9386
      b_mtem(2,ja,je) = 1.3238
      b_mtem(3,ja,je) = 11.8500
      b_mtem(4,ja,je) = -28.1168
      b_mtem(5,ja,je) = 21.8543
      b_mtem(6,ja,je) = -5.1671


      je = jnh4cl
      b_mtem(1,ja,je) = -0.9559
      b_mtem(2,ja,je) = 0.8121
      b_mtem(3,ja,je) = 4.3644
      b_mtem(4,ja,je) = -8.9258
      b_mtem(5,ja,je) = 4.2362
      b_mtem(6,ja,je) = 0.2891


      je = jna2so4
      b_mtem(1,ja,je) = 0.0377
      b_mtem(2,ja,je) = 6.0752
      b_mtem(3,ja,je) = -30.8641
      b_mtem(4,ja,je) = 63.3095
      b_mtem(5,ja,je) = -61.0070
      b_mtem(6,ja,je) = 22.1734


      je = jnano3
      b_mtem(1,ja,je) = -1.8336
      b_mtem(2,ja,je) = 12.8160
      b_mtem(3,ja,je) = -42.3388
      b_mtem(4,ja,je) = 71.1816
      b_mtem(5,ja,je) = -60.5708
      b_mtem(6,ja,je) = 20.5853


      je = jnacl
      b_mtem(1,ja,je) = -0.1429
      b_mtem(2,ja,je) = 2.3561
      b_mtem(3,ja,je) = -10.4425
      b_mtem(4,ja,je) = 20.8951
      b_mtem(5,ja,je) = -20.7739
      b_mtem(6,ja,je) = 7.9355


      je = jcano3
      b_mtem(1,ja,je) = 0.76235
      b_mtem(2,ja,je) = 3.08323
      b_mtem(3,ja,je) = -23.6772
      b_mtem(4,ja,je) = 53.7415
      b_mtem(5,ja,je) = -55.4043
      b_mtem(6,ja,je) = 21.2944


      je = jcacl2
      b_mtem(1,ja,je) = 1.13864
      b_mtem(2,ja,je) = -0.340539
      b_mtem(3,ja,je) = -8.67025
      b_mtem(4,ja,je) = 22.8008
      b_mtem(5,ja,je) = -24.5181
      b_mtem(6,ja,je) = 9.3663


      je = jhno3
      b_mtem(1,ja,je) = 2.42532
      b_mtem(2,ja,je) = -14.1755
      b_mtem(3,ja,je) = 38.804
      b_mtem(4,ja,je) = -58.2437
      b_mtem(5,ja,je) = 43.5431
      b_mtem(6,ja,je) = -12.5824


      je = jhcl
      b_mtem(1,ja,je) = 0.330337
      b_mtem(2,ja,je) = 0.0778934
      b_mtem(3,ja,je) = -2.30492
      b_mtem(4,ja,je) = 4.73003
      b_mtem(5,ja,je) = -4.80849
      b_mtem(6,ja,je) = 1.78866




      ja = jna2so4


      je = jnh4so4
      b_mtem(1,ja,je) = -2.6982
      b_mtem(2,ja,je) = 22.9875
      b_mtem(3,ja,je) = -98.9840
      b_mtem(4,ja,je) = 198.0180
      b_mtem(5,ja,je) = -188.7270
      b_mtem(6,ja,je) = 69.0548


      je = jnh4no3
      b_mtem(1,ja,je) = -2.4844
      b_mtem(2,ja,je) = 6.5420
      b_mtem(3,ja,je) = -9.8998
      b_mtem(4,ja,je) = 11.3884
      b_mtem(5,ja,je) = -13.6842
      b_mtem(6,ja,je) = 7.7411


      je = jnh4cl
      b_mtem(1,ja,je) = -1.3325
      b_mtem(2,ja,je) = 13.0406
      b_mtem(3,ja,je) = -56.1935
      b_mtem(4,ja,je) = 107.1170
      b_mtem(5,ja,je) = -97.3721
      b_mtem(6,ja,je) = 34.3763


      je = jna2so4
      b_mtem(1,ja,je) = -1.2832
      b_mtem(2,ja,je) = 12.8526
      b_mtem(3,ja,je) = -62.2087
      b_mtem(4,ja,je) = 130.3876
      b_mtem(5,ja,je) = -128.2627
      b_mtem(6,ja,je) = 48.0340


      je = jnano3
      b_mtem(1,ja,je) = -3.5384
      b_mtem(2,ja,je) = 21.3758
      b_mtem(3,ja,je) = -70.7638
      b_mtem(4,ja,je) = 121.1580
      b_mtem(5,ja,je) = -104.6230
      b_mtem(6,ja,je) = 36.0557


      je = jnacl
      b_mtem(1,ja,je) = 0.2175
      b_mtem(2,ja,je) = -0.5648
      b_mtem(3,ja,je) = -8.0288
      b_mtem(4,ja,je) = 25.9734
      b_mtem(5,ja,je) = -32.3577
      b_mtem(6,ja,je) = 14.3924


      je = jhno3
      b_mtem(1,ja,je) = -0.309617
      b_mtem(2,ja,je) = -1.82899
      b_mtem(3,ja,je) = -1.5505
      b_mtem(4,ja,je) = 13.3847
      b_mtem(5,ja,je) = -20.1284
      b_mtem(6,ja,je) = 9.93163


      je = jhcl
      b_mtem(1,ja,je) = -0.259455
      b_mtem(2,ja,je) = -0.819366
      b_mtem(3,ja,je) = -4.28964
      b_mtem(4,ja,je) = 16.4305
      b_mtem(5,ja,je) = -21.8546
      b_mtem(6,ja,je) = 10.3044


      je = jh2so4
      b_mtem(1,ja,je) = -1.84257
      b_mtem(2,ja,je) = 7.85788
      b_mtem(3,ja,je) = -29.9275
      b_mtem(4,ja,je) = 61.7515
      b_mtem(5,ja,je) = -63.2308
      b_mtem(6,ja,je) = 24.9542


      je = jnh4hso4
      b_mtem(1,ja,je) = -1.05891
      b_mtem(2,ja,je) = 2.84831
      b_mtem(3,ja,je) = -21.1827
      b_mtem(4,ja,je) = 57.5175
      b_mtem(5,ja,je) = -64.8120
      b_mtem(6,ja,je) = 26.1986


      je = jlvcite
      b_mtem(1,ja,je) = -1.16584
      b_mtem(2,ja,je) = 8.50075
      b_mtem(3,ja,je) = -44.3420
      b_mtem(4,ja,je) = 97.3974
      b_mtem(5,ja,je) = -98.4549
      b_mtem(6,ja,je) = 37.6104


      je = jnahso4
      b_mtem(1,ja,je) = -1.95805
      b_mtem(2,ja,je) = 6.62417
      b_mtem(3,ja,je) = -31.8072
      b_mtem(4,ja,je) = 77.8603
      b_mtem(5,ja,je) = -84.6458
      b_mtem(6,ja,je) = 33.4963


      je = jna3hso4
      b_mtem(1,ja,je) = -0.36045
      b_mtem(2,ja,je) = 3.55223
      b_mtem(3,ja,je) = -24.0327
      b_mtem(4,ja,je) = 54.4879
      b_mtem(5,ja,je) = -56.6531
      b_mtem(6,ja,je) = 22.4956




      ja = jnano3


      je = jnh4so4
      b_mtem(1,ja,je) = -2.5888
      b_mtem(2,ja,je) = 17.6192
      b_mtem(3,ja,je) = -63.2183
      b_mtem(4,ja,je) = 115.3520
      b_mtem(5,ja,je) = -104.0860
      b_mtem(6,ja,je) = 36.7390


      je = jnh4no3
      b_mtem(1,ja,je) = -2.0669
      b_mtem(2,ja,je) = 1.4792
      b_mtem(3,ja,je) = 10.5261
      b_mtem(4,ja,je) = -27.0987
      b_mtem(5,ja,je) = 23.0591
      b_mtem(6,ja,je) = -6.0938


      je = jnh4cl
      b_mtem(1,ja,je) = -0.8325
      b_mtem(2,ja,je) = 3.9933
      b_mtem(3,ja,je) = -15.3789
      b_mtem(4,ja,je) = 30.4050
      b_mtem(5,ja,je) = -29.4204
      b_mtem(6,ja,je) = 11.0597


      je = jna2so4
      b_mtem(1,ja,je) = -1.1233
      b_mtem(2,ja,je) = 8.3998
      b_mtem(3,ja,je) = -31.9002
      b_mtem(4,ja,je) = 60.1450
      b_mtem(5,ja,je) = -55.5503
      b_mtem(6,ja,je) = 19.7757


      je = jnano3
      b_mtem(1,ja,je) = -2.5386
      b_mtem(2,ja,je) = 13.9039
      b_mtem(3,ja,je) = -42.8467
      b_mtem(4,ja,je) = 69.7442
      b_mtem(5,ja,je) = -57.8988
      b_mtem(6,ja,je) = 19.4635


      je = jnacl
      b_mtem(1,ja,je) = -0.4351
      b_mtem(2,ja,je) = 2.8311
      b_mtem(3,ja,je) = -11.4485
      b_mtem(4,ja,je) = 22.7201
      b_mtem(5,ja,je) = -22.4228
      b_mtem(6,ja,je) = 8.5792


      je = jcano3
      b_mtem(1,ja,je) = -0.72060
      b_mtem(2,ja,je) = 5.64915
      b_mtem(3,ja,je) = -23.5020
      b_mtem(4,ja,je) = 46.0078
      b_mtem(5,ja,je) = -43.8075
      b_mtem(6,ja,je) = 16.1652


      je = jcacl2
      b_mtem(1,ja,je) = 0.003928
      b_mtem(2,ja,je) = 3.54724
      b_mtem(3,ja,je) = -18.6057
      b_mtem(4,ja,je) = 38.1445
      b_mtem(5,ja,je) = -36.7745
      b_mtem(6,ja,je) = 13.4529


      je = jhno3
      b_mtem(1,ja,je) = -1.1712
      b_mtem(2,ja,je) = 7.20907
      b_mtem(3,ja,je) = -22.9215
      b_mtem(4,ja,je) = 38.1257
      b_mtem(5,ja,je) = -32.0759
      b_mtem(6,ja,je) = 10.6443


      je = jhcl
      b_mtem(1,ja,je) = 0.738022
      b_mtem(2,ja,je) = -1.14313
      b_mtem(3,ja,je) = 0.32251
      b_mtem(4,ja,je) = 0.838679
      b_mtem(5,ja,je) = -1.81747
      b_mtem(6,ja,je) = 0.873986




      ja = jnacl


      je = jnh4so4
      b_mtem(1,ja,je) = -1.9525
      b_mtem(2,ja,je) = 16.6433
      b_mtem(3,ja,je) = -61.7090
      b_mtem(4,ja,je) = 112.9910
      b_mtem(5,ja,je) = -101.9370
      b_mtem(6,ja,je) = 35.7760


      je = jnh4no3
      b_mtem(1,ja,je) = -1.7525
      b_mtem(2,ja,je) = 3.0713
      b_mtem(3,ja,je) = 4.8063
      b_mtem(4,ja,je) = -17.5334
      b_mtem(5,ja,je) = 14.2872
      b_mtem(6,ja,je) = -3.0690


      je = jnh4cl
      b_mtem(1,ja,je) = -0.4021
      b_mtem(2,ja,je) = 5.2399
      b_mtem(3,ja,je) = -19.4278
      b_mtem(4,ja,je) = 33.0027
      b_mtem(5,ja,je) = -28.1020
      b_mtem(6,ja,je) = 9.5159


      je = jna2so4
      b_mtem(1,ja,je) = 0.6692
      b_mtem(2,ja,je) = 4.1207
      b_mtem(3,ja,je) = -27.3314
      b_mtem(4,ja,je) = 59.3112
      b_mtem(5,ja,je) = -58.7998
      b_mtem(6,ja,je) = 21.7674


      je = jnano3
      b_mtem(1,ja,je) = -1.17444
      b_mtem(2,ja,je) = 10.9927
      b_mtem(3,ja,je) = -38.9013
      b_mtem(4,ja,je) = 66.8521
      b_mtem(5,ja,je) = -57.6564
      b_mtem(6,ja,je) = 19.7296


      je = jnacl
      b_mtem(1,ja,je) = 1.17679
      b_mtem(2,ja,je) = -2.5061
      b_mtem(3,ja,je) = 0.8508
      b_mtem(4,ja,je) = 4.4802
      b_mtem(5,ja,je) = -8.4945
      b_mtem(6,ja,je) = 4.3182


      je = jcano3
      b_mtem(1,ja,je) = 1.01450
      b_mtem(2,ja,je) = 2.10260
      b_mtem(3,ja,je) = -20.9036
      b_mtem(4,ja,je) = 49.1481
      b_mtem(5,ja,je) = -51.4867
      b_mtem(6,ja,je) = 19.9301


      je = jcacl2
      b_mtem(1,ja,je) = 1.55463
      b_mtem(2,ja,je) = -3.20122
      b_mtem(3,ja,je) = -0.957075
      b_mtem(4,ja,je) = 12.103
      b_mtem(5,ja,je) = -17.221
      b_mtem(6,ja,je) = 7.50264


      je = jhno3
      b_mtem(1,ja,je) = 2.46187
      b_mtem(2,ja,je) = -12.6845
      b_mtem(3,ja,je) = 34.2383
      b_mtem(4,ja,je) = -51.9992
      b_mtem(5,ja,je) = 39.4934
      b_mtem(6,ja,je) = -11.7247


      je = jhcl
      b_mtem(1,ja,je) = 1.74915
      b_mtem(2,ja,je) = -4.65768
      b_mtem(3,ja,je) = 8.80287
      b_mtem(4,ja,je) = -12.2503
      b_mtem(5,ja,je) = 8.668751
      b_mtem(6,ja,je) = -2.50158




      ja = jcano3


      je = jnh4no3
      b_mtem(1,ja,je) = -1.86260
      b_mtem(2,ja,je) = 11.6178
      b_mtem(3,ja,je) = -30.9069
      b_mtem(4,ja,je) = 41.7578
      b_mtem(5,ja,je) = -33.7338
      b_mtem(6,ja,je) = 12.7541


      je = jnh4cl
      b_mtem(1,ja,je) = -1.1798
      b_mtem(2,ja,je) = 25.9608
      b_mtem(3,ja,je) = -98.9373
      b_mtem(4,ja,je) = 160.2300
      b_mtem(5,ja,je) = -125.9540
      b_mtem(6,ja,je) = 39.5130


      je = jnano3
      b_mtem(1,ja,je) = -1.44384
      b_mtem(2,ja,je) = 13.6044
      b_mtem(3,ja,je) = -54.4300
      b_mtem(4,ja,je) = 100.582
      b_mtem(5,ja,je) = -91.2364
      b_mtem(6,ja,je) = 32.5970


      je = jnacl
      b_mtem(1,ja,je) = -0.099114
      b_mtem(2,ja,je) = 2.84091
      b_mtem(3,ja,je) = -16.9229
      b_mtem(4,ja,je) = 37.4839
      b_mtem(5,ja,je) = -39.5132
      b_mtem(6,ja,je) = 15.8564


      je = jcano3
      b_mtem(1,ja,je) = 0.055116
      b_mtem(2,ja,je) = 4.58610
      b_mtem(3,ja,je) = -27.6629
      b_mtem(4,ja,je) = 60.8288
      b_mtem(5,ja,je) = -61.4988
      b_mtem(6,ja,je) = 23.3136


      je = jcacl2
      b_mtem(1,ja,je) = 1.57155
      b_mtem(2,ja,je) = -3.18486
      b_mtem(3,ja,je) = -3.35758
      b_mtem(4,ja,je) = 18.7501
      b_mtem(5,ja,je) = -24.5604
      b_mtem(6,ja,je) = 10.3798


      je = jhno3
      b_mtem(1,ja,je) = 1.04446
      b_mtem(2,ja,je) = -3.19066
      b_mtem(3,ja,je) = 2.44714
      b_mtem(4,ja,je) = 2.07218
      b_mtem(5,ja,je) = -6.43949
      b_mtem(6,ja,je) = 3.66471


      je = jhcl
      b_mtem(1,ja,je) = 1.05723
      b_mtem(2,ja,je) = -1.46826
      b_mtem(3,ja,je) = -1.0713
      b_mtem(4,ja,je) = 4.64439
      b_mtem(5,ja,je) = -6.32402
      b_mtem(6,ja,je) = 2.78202




      ja = jcacl2


      je = jnh4no3
      b_mtem(1,ja,je) = -1.43626
      b_mtem(2,ja,je) = 13.6598
      b_mtem(3,ja,je) = -38.2068
      b_mtem(4,ja,je) = 53.9057
      b_mtem(5,ja,je) = -44.9018
      b_mtem(6,ja,je) = 16.6120


      je = jnh4cl
      b_mtem(1,ja,je) = -0.603965
      b_mtem(2,ja,je) = 27.6027
      b_mtem(3,ja,je) = -104.258
      b_mtem(4,ja,je) = 163.553
      b_mtem(5,ja,je) = -124.076
      b_mtem(6,ja,je) = 37.4153


      je = jnano3
      b_mtem(1,ja,je) = 0.44648
      b_mtem(2,ja,je) = 8.8850
      b_mtem(3,ja,je) = -45.5232
      b_mtem(4,ja,je) = 89.3263
      b_mtem(5,ja,je) = -83.8604
      b_mtem(6,ja,je) = 30.4069


      je = jnacl
      b_mtem(1,ja,je) = 1.61927
      b_mtem(2,ja,je) = 0.247547
      b_mtem(3,ja,je) = -18.1252
      b_mtem(4,ja,je) = 45.2479
      b_mtem(5,ja,je) = -48.6072
      b_mtem(6,ja,je) = 19.2784


      je = jcano3
      b_mtem(1,ja,je) = 2.36667
      b_mtem(2,ja,je) = -0.123309
      b_mtem(3,ja,je) = -24.2723
      b_mtem(4,ja,je) = 65.1486
      b_mtem(5,ja,je) = -71.8504
      b_mtem(6,ja,je) = 28.3696


      je = jcacl2
      b_mtem(1,ja,je) = 3.64023
      b_mtem(2,ja,je) = -12.1926
      b_mtem(3,ja,je) = 20.2028
      b_mtem(4,ja,je) = -16.0056
      b_mtem(5,ja,je) = 1.52355
      b_mtem(6,ja,je) = 2.44709


      je = jhno3
      b_mtem(1,ja,je) = 5.88794
      b_mtem(2,ja,je) = -29.7083
      b_mtem(3,ja,je) = 78.6309
      b_mtem(4,ja,je) = -118.037
      b_mtem(5,ja,je) = 88.932
      b_mtem(6,ja,je) = -26.1407


      je = jhcl
      b_mtem(1,ja,je) = 2.40628
      b_mtem(2,ja,je) = -6.16566
      b_mtem(3,ja,je) = 10.2851
      b_mtem(4,ja,je) = -12.9035
      b_mtem(5,ja,je) = 7.7441
      b_mtem(6,ja,je) = -1.74821




      ja = jhno3


      je = jnh4so4
      b_mtem(1,ja,je) = -3.57598
      b_mtem(2,ja,je) = 21.5469
      b_mtem(3,ja,je) = -77.4111
      b_mtem(4,ja,je) = 144.136
      b_mtem(5,ja,je) = -132.849
      b_mtem(6,ja,je) = 47.9412


      je = jnh4no3
      b_mtem(1,ja,je) = -2.00209
      b_mtem(2,ja,je) = -3.48399
      b_mtem(3,ja,je) = 34.9906
      b_mtem(4,ja,je) = -68.6653
      b_mtem(5,ja,je) = 54.0992
      b_mtem(6,ja,je) = -15.1343


      je = jnh4cl
      b_mtem(1,ja,je) = -0.63790
      b_mtem(2,ja,je) = -1.67730
      b_mtem(3,ja,je) = 10.1727
      b_mtem(4,ja,je) = -14.9097
      b_mtem(5,ja,je) = 7.67410
      b_mtem(6,ja,je) = -0.79586


      je = jnacl
      b_mtem(1,ja,je) = 1.3446
      b_mtem(2,ja,je) = -2.5578
      b_mtem(3,ja,je) = 1.3464
      b_mtem(4,ja,je) = 2.90537
      b_mtem(5,ja,je) = -6.53014
      b_mtem(6,ja,je) = 3.31339


      je = jnano3
      b_mtem(1,ja,je) = -0.546636
      b_mtem(2,ja,je) = 10.3127
      b_mtem(3,ja,je) = -39.9603
      b_mtem(4,ja,je) = 71.4609
      b_mtem(5,ja,je) = -63.4958
      b_mtem(6,ja,je) = 22.0679


      je = jna2so4
      b_mtem(1,ja,je) = 1.35059
      b_mtem(2,ja,je) = 4.34557
      b_mtem(3,ja,je) = -35.8425
      b_mtem(4,ja,je) = 80.9868
      b_mtem(5,ja,je) = -81.6544
      b_mtem(6,ja,je) = 30.4841


      je = jcano3
      b_mtem(1,ja,je) = 0.869414
      b_mtem(2,ja,je) = 2.98486
      b_mtem(3,ja,je) = -22.255
      b_mtem(4,ja,je) = 50.1863
      b_mtem(5,ja,je) = -51.214
      b_mtem(6,ja,je) = 19.2235


      je = jcacl2
      b_mtem(1,ja,je) = 1.42800
      b_mtem(2,ja,je) = -1.78959
      b_mtem(3,ja,je) = -2.49075
      b_mtem(4,ja,je) = 10.1877
      b_mtem(5,ja,je) = -12.1948
      b_mtem(6,ja,je) = 4.64475


      je = jhno3
      b_mtem(1,ja,je) = 0.22035
      b_mtem(2,ja,je) = 2.94973
      b_mtem(3,ja,je) = -12.1469
      b_mtem(4,ja,je) = 20.4905
      b_mtem(5,ja,je) = -17.3966
      b_mtem(6,ja,je) = 5.70779


      je = jhcl
      b_mtem(1,ja,je) = 1.55503
      b_mtem(2,ja,je) = -3.61226
      b_mtem(3,ja,je) = 6.28265
      b_mtem(4,ja,je) = -8.69575
      b_mtem(5,ja,je) = 6.09372
      b_mtem(6,ja,je) = -1.80898


      je = jh2so4
      b_mtem(1,ja,je) = 1.10783
      b_mtem(2,ja,je) = -1.3363
      b_mtem(3,ja,je) = -1.83525
      b_mtem(4,ja,je) = 7.47373
      b_mtem(5,ja,je) = -9.72954
      b_mtem(6,ja,je) = 4.12248


      je = jnh4hso4
      b_mtem(1,ja,je) = -0.851026
      b_mtem(2,ja,je) = 12.2515
      b_mtem(3,ja,je) = -49.788
      b_mtem(4,ja,je) = 91.6215
      b_mtem(5,ja,je) = -81.4877
      b_mtem(6,ja,je) = 28.0002


      je = jlvcite
      b_mtem(1,ja,je) = -3.09464
      b_mtem(2,ja,je) = 14.9303
      b_mtem(3,ja,je) = -43.0454
      b_mtem(4,ja,je) = 72.6695
      b_mtem(5,ja,je) = -65.2140
      b_mtem(6,ja,je) = 23.4814


      je = jnahso4
      b_mtem(1,ja,je) = 1.22973
      b_mtem(2,ja,je) = 2.82702
      b_mtem(3,ja,je) = -17.5869
      b_mtem(4,ja,je) = 28.9564
      b_mtem(5,ja,je) = -23.5814
      b_mtem(6,ja,je) = 7.91153


      je = jna3hso4
      b_mtem(1,ja,je) = 1.64773
      b_mtem(2,ja,je) = 0.94188
      b_mtem(3,ja,je) = -19.1242
      b_mtem(4,ja,je) = 46.9887
      b_mtem(5,ja,je) = -50.9494
      b_mtem(6,ja,je) = 20.2169




      ja = jhcl


      je = jnh4so4
      b_mtem(1,ja,je) = -2.93783
      b_mtem(2,ja,je) = 20.5546
      b_mtem(3,ja,je) = -75.8548
      b_mtem(4,ja,je) = 141.729
      b_mtem(5,ja,je) = -130.697
      b_mtem(6,ja,je) = 46.9905


      je = jnh4no3
      b_mtem(1,ja,je) = -1.69063
      b_mtem(2,ja,je) = -1.85303
      b_mtem(3,ja,je) = 29.0927
      b_mtem(4,ja,je) = -58.7401
      b_mtem(5,ja,je) = 44.999
      b_mtem(6,ja,je) = -11.9988


      je = jnh4cl
      b_mtem(1,ja,je) = -0.2073
      b_mtem(2,ja,je) = -0.4322
      b_mtem(3,ja,je) = 6.1271
      b_mtem(4,ja,je) = -12.3146
      b_mtem(5,ja,je) = 8.9919
      b_mtem(6,ja,je) = -2.3388


      je = jnacl
      b_mtem(1,ja,je) = 2.95913
      b_mtem(2,ja,je) = -7.92254
      b_mtem(3,ja,je) = 13.736
      b_mtem(4,ja,je) = -15.433
      b_mtem(5,ja,je) = 7.40386
      b_mtem(6,ja,je) = -0.918641


      je = jnano3
      b_mtem(1,ja,je) = 0.893272
      b_mtem(2,ja,je) = 6.53768
      b_mtem(3,ja,je) = -32.3458
      b_mtem(4,ja,je) = 61.2834
      b_mtem(5,ja,je) = -56.4446
      b_mtem(6,ja,je) = 19.9202


      je = jna2so4
      b_mtem(1,ja,je) = 3.14484
      b_mtem(2,ja,je) = 0.077019
      b_mtem(3,ja,je) = -31.4199
      b_mtem(4,ja,je) = 80.5865
      b_mtem(5,ja,je) = -85.392
      b_mtem(6,ja,je) = 32.6644


      je = jcano3
      b_mtem(1,ja,je) = 2.60432
      b_mtem(2,ja,je) = -0.55909
      b_mtem(3,ja,je) = -19.6671
      b_mtem(4,ja,je) = 53.3446
      b_mtem(5,ja,je) = -58.9076
      b_mtem(6,ja,je) = 22.9927


      je = jcacl2
      b_mtem(1,ja,je) = 2.98036
      b_mtem(2,ja,je) = -8.55365
      b_mtem(3,ja,je) = 15.2108
      b_mtem(4,ja,je) = -15.9359
      b_mtem(5,ja,je) = 7.41772
      b_mtem(6,ja,je) = -1.32143


      je = jhno3
      b_mtem(1,ja,je) = 3.8533
      b_mtem(2,ja,je) = -16.9427
      b_mtem(3,ja,je) = 45.0056
      b_mtem(4,ja,je) = -69.6145
      b_mtem(5,ja,je) = 54.1491
      b_mtem(6,ja,je) = -16.6513


      je = jhcl
      b_mtem(1,ja,je) = 2.56665
      b_mtem(2,ja,je) = -7.13585
      b_mtem(3,ja,je) = 14.8103
      b_mtem(4,ja,je) = -21.8881
      b_mtem(5,ja,je) = 16.6808
      b_mtem(6,ja,je) = -5.22091


      je = jh2so4
      b_mtem(1,ja,je) = 2.50179
      b_mtem(2,ja,je) = -6.69364
      b_mtem(3,ja,je) = 11.6551
      b_mtem(4,ja,je) = -13.6897
      b_mtem(5,ja,je) = 7.36796
      b_mtem(6,ja,je) = -1.33245


      je = jnh4hso4
      b_mtem(1,ja,je) = 0.149955
      b_mtem(2,ja,je) = 11.8213
      b_mtem(3,ja,je) = -53.9164
      b_mtem(4,ja,je) = 101.574
      b_mtem(5,ja,je) = -91.4123
      b_mtem(6,ja,je) = 31.5487


      je = jlvcite
      b_mtem(1,ja,je) = -2.36927
      b_mtem(2,ja,je) = 14.8359
      b_mtem(3,ja,je) = -44.3443
      b_mtem(4,ja,je) = 73.6229
      b_mtem(5,ja,je) = -65.3366
      b_mtem(6,ja,je) = 23.3250


      je = jnahso4
      b_mtem(1,ja,je) = 2.72993
      b_mtem(2,ja,je) = -0.23406
      b_mtem(3,ja,je) = -10.4103
      b_mtem(4,ja,je) = 13.1586
      b_mtem(5,ja,je) = -7.79925
      b_mtem(6,ja,je) = 2.30843


      je = jna3hso4
      b_mtem(1,ja,je) = 3.51258
      b_mtem(2,ja,je) = -3.95107
      b_mtem(3,ja,je) = -11.0175
      b_mtem(4,ja,je) = 38.8617
      b_mtem(5,ja,je) = -48.1575
      b_mtem(6,ja,je) = 20.4717




      ja = jh2so4


      je = jh2so4
      b_mtem(1,ja,je) = 0.76734
      b_mtem(2,ja,je) = -1.12263
      b_mtem(3,ja,je) = -9.08728
      b_mtem(4,ja,je) = 30.3836
      b_mtem(5,ja,je) = -38.4133
      b_mtem(6,ja,je) = 17.0106


      je = jnh4hso4
      b_mtem(1,ja,je) = -2.03879
      b_mtem(2,ja,je) = 15.7033
      b_mtem(3,ja,je) = -58.7363
      b_mtem(4,ja,je) = 109.242
      b_mtem(5,ja,je) = -102.237
      b_mtem(6,ja,je) = 37.5350


      je = jlvcite
      b_mtem(1,ja,je) = -3.10228
      b_mtem(2,ja,je) = 16.6920
      b_mtem(3,ja,je) = -59.1522
      b_mtem(4,ja,je) = 113.487
      b_mtem(5,ja,je) = -110.890
      b_mtem(6,ja,je) = 42.4578


      je = jnh4so4
      b_mtem(1,ja,je) = -3.43885
      b_mtem(2,ja,je) = 21.0372
      b_mtem(3,ja,je) = -84.7026
      b_mtem(4,ja,je) = 165.324
      b_mtem(5,ja,je) = -156.101
      b_mtem(6,ja,je) = 57.3101


      je = jnahso4
      b_mtem(1,ja,je) = 0.33164
      b_mtem(2,ja,je) = 6.55864
      b_mtem(3,ja,je) = -33.5876
      b_mtem(4,ja,je) = 65.1798
      b_mtem(5,ja,je) = -63.2046
      b_mtem(6,ja,je) = 24.1783


      je = jna3hso4
      b_mtem(1,ja,je) = 3.06830
      b_mtem(2,ja,je) = -3.18408
      b_mtem(3,ja,je) = -19.6332
      b_mtem(4,ja,je) = 61.3657
      b_mtem(5,ja,je) = -73.4438
      b_mtem(6,ja,je) = 31.2334


      je = jna2so4
      b_mtem(1,ja,je) = 2.58649
      b_mtem(2,ja,je) = 0.87921
      b_mtem(3,ja,je) = -39.3023
      b_mtem(4,ja,je) = 101.603
      b_mtem(5,ja,je) = -109.469
      b_mtem(6,ja,je) = 43.0188


      je = jhno3
      b_mtem(1,ja,je) = 1.54587
      b_mtem(2,ja,je) = -7.50976
      b_mtem(3,ja,je) = 12.8237
      b_mtem(4,ja,je) = -10.1452
      b_mtem(5,ja,je) = -0.541956
      b_mtem(6,ja,je) = 3.34536


      je = jhcl
      b_mtem(1,ja,je) = 0.829757
      b_mtem(2,ja,je) = -4.11316
      b_mtem(3,ja,je) = 3.67111
      b_mtem(4,ja,je) = 3.6833
      b_mtem(5,ja,je) = -11.2711
      b_mtem(6,ja,je) = 6.71421




      ja = jhhso4


      je = jh2so4
      b_mtem(1,ja,je) = 2.63953
      b_mtem(2,ja,je) = -6.01532
      b_mtem(3,ja,je) = 10.0204
      b_mtem(4,ja,je) = -12.4840
      b_mtem(5,ja,je) = 7.78853
      b_mtem(6,ja,je) = -2.12638


      je = jnh4hso4
      b_mtem(1,ja,je) = -0.77412
      b_mtem(2,ja,je) = 14.1656
      b_mtem(3,ja,je) = -53.4087
      b_mtem(4,ja,je) = 93.2013
      b_mtem(5,ja,je) = -80.5723
      b_mtem(6,ja,je) = 27.1577


      je = jlvcite
      b_mtem(1,ja,je) = -2.98882
      b_mtem(2,ja,je) = 14.4436
      b_mtem(3,ja,je) = -40.1774
      b_mtem(4,ja,je) = 67.5937
      b_mtem(5,ja,je) = -61.5040
      b_mtem(6,ja,je) = 22.3695


      je = jnh4so4
      b_mtem(1,ja,je) = -1.15502
      b_mtem(2,ja,je) = 8.12309
      b_mtem(3,ja,je) = -38.4726
      b_mtem(4,ja,je) = 80.8861
      b_mtem(5,ja,je) = -80.1644
      b_mtem(6,ja,je) = 30.4717


      je = jnahso4
      b_mtem(1,ja,je) = 1.99641
      b_mtem(2,ja,je) = -2.96061
      b_mtem(3,ja,je) = 5.54778
      b_mtem(4,ja,je) = -14.5488
      b_mtem(5,ja,je) = 14.8492
      b_mtem(6,ja,je) = -5.1389


      je = jna3hso4
      b_mtem(1,ja,je) = 2.23816
      b_mtem(2,ja,je) = -3.20847
      b_mtem(3,ja,je) = -4.82853
      b_mtem(4,ja,je) = 20.9192
      b_mtem(5,ja,je) = -27.2819
      b_mtem(6,ja,je) = 11.8655


      je = jna2so4
      b_mtem(1,ja,je) = 2.56907
      b_mtem(2,ja,je) = 1.13444
      b_mtem(3,ja,je) = -34.6853
      b_mtem(4,ja,je) = 87.9775
      b_mtem(5,ja,je) = -93.2330
      b_mtem(6,ja,je) = 35.9260


      je = jhno3
      b_mtem(1,ja,je) = 2.00024
      b_mtem(2,ja,je) = -4.80868
      b_mtem(3,ja,je) = 8.29222
      b_mtem(4,ja,je) = -11.0849
      b_mtem(5,ja,je) = 7.51262
      b_mtem(6,ja,je) = -2.07654


      je = jhcl
      b_mtem(1,ja,je) = 2.8009
      b_mtem(2,ja,je) = -6.98416
      b_mtem(3,ja,je) = 14.3146
      b_mtem(4,ja,je) = -22.0068
      b_mtem(5,ja,je) = 17.5557
      b_mtem(6,ja,je) = -5.84917




      ja = jnh4hso4


      je = jh2so4
      b_mtem(1,ja,je) = 0.169160
      b_mtem(2,ja,je) = 2.15094
      b_mtem(3,ja,je) = -9.62904
      b_mtem(4,ja,je) = 18.2631
      b_mtem(5,ja,je) = -17.3333
      b_mtem(6,ja,je) = 6.19835


      je = jnh4hso4
      b_mtem(1,ja,je) = -2.34457
      b_mtem(2,ja,je) = 12.8035
      b_mtem(3,ja,je) = -35.2513
      b_mtem(4,ja,je) = 53.6153
      b_mtem(5,ja,je) = -42.7655
      b_mtem(6,ja,je) = 13.7129


      je = jlvcite
      b_mtem(1,ja,je) = -2.56109
      b_mtem(2,ja,je) = 11.1414
      b_mtem(3,ja,je) = -30.2361
      b_mtem(4,ja,je) = 50.0320
      b_mtem(5,ja,je) = -44.1586
      b_mtem(6,ja,je) = 15.5393


      je = jnh4so4
      b_mtem(1,ja,je) = -0.97315
      b_mtem(2,ja,je) = 7.06295
      b_mtem(3,ja,je) = -29.3032
      b_mtem(4,ja,je) = 57.6101
      b_mtem(5,ja,je) = -54.9020
      b_mtem(6,ja,je) = 20.2222


      je = jnahso4
      b_mtem(1,ja,je) = -0.44450
      b_mtem(2,ja,je) = 3.33451
      b_mtem(3,ja,je) = -15.2791
      b_mtem(4,ja,je) = 30.1413
      b_mtem(5,ja,je) = -26.7710
      b_mtem(6,ja,je) = 8.78462


      je = jna3hso4
      b_mtem(1,ja,je) = -0.99780
      b_mtem(2,ja,je) = 4.69200
      b_mtem(3,ja,je) = -16.1219
      b_mtem(4,ja,je) = 29.3100
      b_mtem(5,ja,je) = -26.3383
      b_mtem(6,ja,je) = 9.20695


      je = jna2so4
      b_mtem(1,ja,je) = -0.52694
      b_mtem(2,ja,je) = 7.02684
      b_mtem(3,ja,je) = -33.7508
      b_mtem(4,ja,je) = 70.0565
      b_mtem(5,ja,je) = -68.3226
      b_mtem(6,ja,je) = 25.2692


      je = jhno3
      b_mtem(1,ja,je) = 0.572926
      b_mtem(2,ja,je) = -2.04791
      b_mtem(3,ja,je) = 2.1134
      b_mtem(4,ja,je) = 0.246654
      b_mtem(5,ja,je) = -3.06019
      b_mtem(6,ja,je) = 1.98126


      je = jhcl
      b_mtem(1,ja,je) = 0.56514
      b_mtem(2,ja,je) = 0.22287
      b_mtem(3,ja,je) = -2.76973
      b_mtem(4,ja,je) = 4.54444
      b_mtem(5,ja,je) = -3.86549
      b_mtem(6,ja,je) = 1.13441




      ja = jlvcite


      je = jh2so4
      b_mtem(1,ja,je) = -1.44811
      b_mtem(2,ja,je) = 6.71815
      b_mtem(3,ja,je) = -25.0141
      b_mtem(4,ja,je) = 50.1109
      b_mtem(5,ja,je) = -50.0561
      b_mtem(6,ja,je) = 19.3370


      je = jnh4hso4
      b_mtem(1,ja,je) = -3.41707
      b_mtem(2,ja,je) = 13.4496
      b_mtem(3,ja,je) = -34.8018
      b_mtem(4,ja,je) = 55.2987
      b_mtem(5,ja,je) = -48.1839
      b_mtem(6,ja,je) = 17.2444


      je = jlvcite
      b_mtem(1,ja,je) = -2.54479
      b_mtem(2,ja,je) = 11.8501
      b_mtem(3,ja,je) = -39.7286
      b_mtem(4,ja,je) = 74.2479
      b_mtem(5,ja,je) = -70.4934
      b_mtem(6,ja,je) = 26.2836


      je = jnh4so4
      b_mtem(1,ja,je) = -2.30561
      b_mtem(2,ja,je) = 14.5806
      b_mtem(3,ja,je) = -55.1238
      b_mtem(4,ja,je) = 103.451
      b_mtem(5,ja,je) = -95.2571
      b_mtem(6,ja,je) = 34.2218


      je = jnahso4
      b_mtem(1,ja,je) = -2.20809
      b_mtem(2,ja,je) = 13.6391
      b_mtem(3,ja,je) = -57.8246
      b_mtem(4,ja,je) = 117.907
      b_mtem(5,ja,je) = -112.154
      b_mtem(6,ja,je) = 40.3058


      je = jna3hso4
      b_mtem(1,ja,je) = -1.15099
      b_mtem(2,ja,je) = 6.32269
      b_mtem(3,ja,je) = -27.3860
      b_mtem(4,ja,je) = 55.4592
      b_mtem(5,ja,je) = -54.0100
      b_mtem(6,ja,je) = 20.3469


      je = jna2so4
      b_mtem(1,ja,je) = -1.15678
      b_mtem(2,ja,je) = 8.28718
      b_mtem(3,ja,je) = -37.3231
      b_mtem(4,ja,je) = 76.6124
      b_mtem(5,ja,je) = -74.9307
      b_mtem(6,ja,je) = 28.0559


      je = jhno3
      b_mtem(1,ja,je) = 0.01502
      b_mtem(2,ja,je) = -3.1197
      b_mtem(3,ja,je) = 3.61104
      b_mtem(4,ja,je) = 3.05196
      b_mtem(5,ja,je) = -9.98957
      b_mtem(6,ja,je) = 6.04155


      je = jhcl
      b_mtem(1,ja,je) = -1.06477
      b_mtem(2,ja,je) = 3.38801
      b_mtem(3,ja,je) = -12.5784
      b_mtem(4,ja,je) = 25.2823
      b_mtem(5,ja,je) = -25.4611
      b_mtem(6,ja,je) = 10.0754




      ja = jnahso4


      je = jh2so4
      b_mtem(1,ja,je) = 0.68259
      b_mtem(2,ja,je) = 0.71468
      b_mtem(3,ja,je) = -5.59003
      b_mtem(4,ja,je) = 11.0089
      b_mtem(5,ja,je) = -10.7983
      b_mtem(6,ja,je) = 3.82335


      je = jnh4hso4
      b_mtem(1,ja,je) = -0.03956
      b_mtem(2,ja,je) = 4.52828
      b_mtem(3,ja,je) = -25.2557
      b_mtem(4,ja,je) = 54.4225
      b_mtem(5,ja,je) = -52.5105
      b_mtem(6,ja,je) = 18.6562


      je = jlvcite
      b_mtem(1,ja,je) = -1.53503
      b_mtem(2,ja,je) = 8.27608
      b_mtem(3,ja,je) = -28.9539
      b_mtem(4,ja,je) = 55.2876
      b_mtem(5,ja,je) = -51.9563
      b_mtem(6,ja,je) = 18.6576


      je = jnh4so4
      b_mtem(1,ja,je) = -0.38793
      b_mtem(2,ja,je) = 7.14680
      b_mtem(3,ja,je) = -38.7201
      b_mtem(4,ja,je) = 84.3965
      b_mtem(5,ja,je) = -84.7453
      b_mtem(6,ja,je) = 32.1283


      je = jnahso4
      b_mtem(1,ja,je) = -0.41982
      b_mtem(2,ja,je) = 4.26491
      b_mtem(3,ja,je) = -20.2351
      b_mtem(4,ja,je) = 42.6764
      b_mtem(5,ja,je) = -40.7503
      b_mtem(6,ja,je) = 14.2868


      je = jna3hso4
      b_mtem(1,ja,je) = -0.32912
      b_mtem(2,ja,je) = 1.80808
      b_mtem(3,ja,je) = -8.01286
      b_mtem(4,ja,je) = 15.5791
      b_mtem(5,ja,je) = -14.5494
      b_mtem(6,ja,je) = 5.27052


      je = jna2so4
      b_mtem(1,ja,je) = 0.10271
      b_mtem(2,ja,je) = 5.09559
      b_mtem(3,ja,je) = -30.3295
      b_mtem(4,ja,je) = 66.2975
      b_mtem(5,ja,je) = -66.3458
      b_mtem(6,ja,je) = 24.9443


      je = jhno3
      b_mtem(1,ja,je) = 0.608309
      b_mtem(2,ja,je) = -0.541905
      b_mtem(3,ja,je) = -2.52084
      b_mtem(4,ja,je) = 6.63297
      b_mtem(5,ja,je) = -7.24599
      b_mtem(6,ja,je) = 2.88811


      je = jhcl
      b_mtem(1,ja,je) = 1.98399
      b_mtem(2,ja,je) = -4.51562
      b_mtem(3,ja,je) = 8.36059
      b_mtem(4,ja,je) = -12.4948
      b_mtem(5,ja,je) = 9.67514
      b_mtem(6,ja,je) = -3.18004




      ja = jna3hso4


      je = jh2so4
      b_mtem(1,ja,je) = -0.83214
      b_mtem(2,ja,je) = 4.99572
      b_mtem(3,ja,je) = -20.1697
      b_mtem(4,ja,je) = 41.4066
      b_mtem(5,ja,je) = -42.2119
      b_mtem(6,ja,je) = 16.4855


      je = jnh4hso4
      b_mtem(1,ja,je) = -0.65139
      b_mtem(2,ja,je) = 3.52300
      b_mtem(3,ja,je) = -22.8220
      b_mtem(4,ja,je) = 56.2956
      b_mtem(5,ja,je) = -59.9028
      b_mtem(6,ja,je) = 23.1844


      je = jlvcite
      b_mtem(1,ja,je) = -1.31331
      b_mtem(2,ja,je) = 8.40835
      b_mtem(3,ja,je) = -38.1757
      b_mtem(4,ja,je) = 80.5312
      b_mtem(5,ja,je) = -79.8346
      b_mtem(6,ja,je) = 30.0219


      je = jnh4so4
      b_mtem(1,ja,je) = -1.03054
      b_mtem(2,ja,je) = 8.08155
      b_mtem(3,ja,je) = -38.1046
      b_mtem(4,ja,je) = 78.7168
      b_mtem(5,ja,je) = -77.2263
      b_mtem(6,ja,je) = 29.1521


      je = jnahso4
      b_mtem(1,ja,je) = -1.90695
      b_mtem(2,ja,je) = 11.6241
      b_mtem(3,ja,je) = -50.3175
      b_mtem(4,ja,je) = 105.884
      b_mtem(5,ja,je) = -103.258
      b_mtem(6,ja,je) = 37.6588


      je = jna3hso4
      b_mtem(1,ja,je) = -0.34780
      b_mtem(2,ja,je) = 2.85363
      b_mtem(3,ja,je) = -17.6224
      b_mtem(4,ja,je) = 38.9220
      b_mtem(5,ja,je) = -39.8106
      b_mtem(6,ja,je) = 15.6055


      je = jna2so4
      b_mtem(1,ja,je) = -0.75230
      b_mtem(2,ja,je) = 10.0140
      b_mtem(3,ja,je) = -50.5677
      b_mtem(4,ja,je) = 106.941
      b_mtem(5,ja,je) = -105.534
      b_mtem(6,ja,je) = 39.5196


      je = jhno3
      b_mtem(1,ja,je) = 0.057456
      b_mtem(2,ja,je) = -1.31264
      b_mtem(3,ja,je) = -1.94662
      b_mtem(4,ja,je) = 10.7024
      b_mtem(5,ja,je) = -14.9946
      b_mtem(6,ja,je) = 7.12161


      je = jhcl
      b_mtem(1,ja,je) = 0.637894
      b_mtem(2,ja,je) = -2.29719
      b_mtem(3,ja,je) = 0.765361
      b_mtem(4,ja,je) = 4.8748
      b_mtem(5,ja,je) = -9.25978
      b_mtem(6,ja,je) = 4.91773










      j_index = 1
      d_mdrh(j_index,1) = -58.00268351
      d_mdrh(j_index,2) = 2.031077573
      d_mdrh(j_index,3) = -0.008281218
      d_mdrh(j_index,4) = 1.00447e-05


      j_index = 2
      d_mdrh(j_index,1) = 1039.137773
      d_mdrh(j_index,2) = -11.47847095
      d_mdrh(j_index,3) = 0.047702786
      d_mdrh(j_index,4) = -6.77675e-05


      j_index = 3
      d_mdrh(j_index,1) = 115.8366357
      d_mdrh(j_index,2) = 0.491881663
      d_mdrh(j_index,3) = -0.00422807
      d_mdrh(j_index,4) = 7.29274e-06


      j_index = 4
      d_mdrh(j_index,1) = 253.2424151
      d_mdrh(j_index,2) = -1.429957864
      d_mdrh(j_index,3) = 0.003727554
      d_mdrh(j_index,4) = -3.13037e-06


      j_index = 5
      d_mdrh(j_index,1) = -372.4306506
      d_mdrh(j_index,2) = 5.3955633
      d_mdrh(j_index,3) = -0.019804438
      d_mdrh(j_index,4) = 2.25662e-05


      j_index = 6
      d_mdrh(j_index,1) = 286.1271416
      d_mdrh(j_index,2) = -1.670787758
      d_mdrh(j_index,3) = 0.004431373
      d_mdrh(j_index,4) = -3.57757e-06


      j_index = 7
      d_mdrh(j_index,1) = -1124.07059
      d_mdrh(j_index,2) = 14.26364209
      d_mdrh(j_index,3) = -0.054816822
      d_mdrh(j_index,4) = 6.70107e-05


      j_index = 8
      d_mdrh(j_index,1) = 1855.413934
      d_mdrh(j_index,2) = -20.29219473
      d_mdrh(j_index,3) = 0.07807482
      d_mdrh(j_index,4) = -1.017887858e-4


      j_index = 9
      d_mdrh(j_index,1) = 1761.176886
      d_mdrh(j_index,2) = -19.29811062
      d_mdrh(j_index,3) = 0.075676987
      d_mdrh(j_index,4) = -1.0116959e-4


      j_index = 10
      d_mdrh(j_index,1) = 122.1074303
      d_mdrh(j_index,2) = 0.429692122
      d_mdrh(j_index,3) = -0.003928277
      d_mdrh(j_index,4) = 6.43275e-06


      j_index = 11
      d_mdrh(j_index,1) = 2424.634678
      d_mdrh(j_index,2) = -26.54031307
      d_mdrh(j_index,3) = 0.101625387
      d_mdrh(j_index,4) = -1.31544547798e-4


      j_index = 12
      d_mdrh(j_index,1) = 2912.082599
      d_mdrh(j_index,2) = -31.8894185
      d_mdrh(j_index,3) = 0.121185849
      d_mdrh(j_index,4) = -1.556534623e-4


      j_index = 13
      d_mdrh(j_index,1) = 172.2596493
      d_mdrh(j_index,2) = -0.511006195
      d_mdrh(j_index,3) = 4.27244597e-4
      d_mdrh(j_index,4) = 4.12797e-07


      j_index = 14
      d_mdrh(j_index,1) = 1596.184935
      d_mdrh(j_index,2) = -16.37945565
      d_mdrh(j_index,3) = 0.060281218
      d_mdrh(j_index,4) = -7.6161e-05


      j_index = 15
      d_mdrh(j_index,1) = 1916.072988
      d_mdrh(j_index,2) = -20.85594868
      d_mdrh(j_index,3) = 0.081140141
      d_mdrh(j_index,4) = -1.07954274796e-4


      j_index = 16
      d_mdrh(j_index,1) = 1467.165935
      d_mdrh(j_index,2) = -16.01166196
      d_mdrh(j_index,3) = 0.063505582
      d_mdrh(j_index,4) = -8.66722e-05


      j_index = 17
      d_mdrh(j_index,1) = 158.447059
      d_mdrh(j_index,2) = -0.628167358
      d_mdrh(j_index,3) = 0.002014448
      d_mdrh(j_index,4) = -3.13037e-06


      j_index = 18
      d_mdrh(j_index,1) = 1115.892468
      d_mdrh(j_index,2) = -11.76936534
      d_mdrh(j_index,3) = 0.045577399
      d_mdrh(j_index,4) = -6.05779e-05


      j_index = 19
      d_mdrh(j_index,1) = 269.5432407
      d_mdrh(j_index,2) = -1.319963885
      d_mdrh(j_index,3) = 0.002592363
      d_mdrh(j_index,4) = -1.44479e-06


      j_index = 20
      d_mdrh(j_index,1) = 2841.334784
      d_mdrh(j_index,2) = -31.1889487
      d_mdrh(j_index,3) = 0.118809274
      d_mdrh(j_index,4) = -1.53007e-4


      j_index = 21
      d_mdrh(j_index,1) = 2199.36914
      d_mdrh(j_index,2) = -24.11926569
      d_mdrh(j_index,3) = 0.092932361
      d_mdrh(j_index,4) = -1.21774e-4


      j_index = 22
      d_mdrh(j_index,1) = 395.0051604
      d_mdrh(j_index,2) = -2.521101657
      d_mdrh(j_index,3) = 0.006139319
      d_mdrh(j_index,4) = -4.43756e-06


      j_index = 23
      d_mdrh(j_index,1) = 386.5150675
      d_mdrh(j_index,2) = -2.4632138
      d_mdrh(j_index,3) = 0.006139319
      d_mdrh(j_index,4) = -4.98796e-06


      j_index = 24
      d_mdrh(j_index,1) = 3101.538491
      d_mdrh(j_index,2) = -34.19978105
      d_mdrh(j_index,3) = 0.130118605
      d_mdrh(j_index,4) = -1.66873e-4


      j_index = 25
      d_mdrh(j_index,1) = 2307.579403
      d_mdrh(j_index,2) = -25.43136774
      d_mdrh(j_index,3) = 0.098064728
      d_mdrh(j_index,4) = -1.28301e-4


      j_index = 26
      d_mdrh(j_index,1) = 291.8309602
      d_mdrh(j_index,2) = -1.828912974
      d_mdrh(j_index,3) = 0.005053148
      d_mdrh(j_index,4) = -4.57516e-06


      j_index = 27
      d_mdrh(j_index,1) = 188.3914345
      d_mdrh(j_index,2) = -0.631345031
      d_mdrh(j_index,3) = 0.000622807
      d_mdrh(j_index,4) = 4.47196e-07


      j_index = 28
      d_mdrh(j_index,1) = -167.1252839
      d_mdrh(j_index,2) = 2.969828002
      d_mdrh(j_index,3) = -0.010637255
      d_mdrh(j_index,4) = 1.13175e-05


      j_index = 29
      d_mdrh(j_index,1) = 1516.782768
      d_mdrh(j_index,2) = -15.7922661
      d_mdrh(j_index,3) = 0.058942209
      d_mdrh(j_index,4) = -7.5301e-05


      j_index = 30
      d_mdrh(j_index,1) = 1739.963163
      d_mdrh(j_index,2) = -19.06576022
      d_mdrh(j_index,3) = 0.07454963
      d_mdrh(j_index,4) = -9.94302e-05


      j_index = 31
      d_mdrh(j_index,1) = 2152.104877
      d_mdrh(j_index,2) = -23.74998008
      d_mdrh(j_index,3) = 0.092256654
      d_mdrh(j_index,4) = -1.21953e-4


      j_index = 32
      d_mdrh(j_index,1) = 221.9976265
      d_mdrh(j_index,2) = -1.311331272
      d_mdrh(j_index,3) = 0.004406089
      d_mdrh(j_index,4) = -5.88235e-06


      j_index = 33
      d_mdrh(j_index,1) = 1205.645615
      d_mdrh(j_index,2) = -12.71353459
      d_mdrh(j_index,3) = 0.048803922
      d_mdrh(j_index,4) = -6.41899e-05


      j_index = 34
      d_mdrh(j_index,1) = 506.6737879
      d_mdrh(j_index,2) = -3.723520818
      d_mdrh(j_index,3) = 0.010814242
      d_mdrh(j_index,4) = -1.21087e-05


      j_index = 35
      d_mdrh(j_index,1) = -1123.523841
      d_mdrh(j_index,2) = 14.08345977
      d_mdrh(j_index,3) = -0.053687823
      d_mdrh(j_index,4) = 6.52219e-05


      j_index = 36
      d_mdrh(j_index,1) = -1159.98607
      d_mdrh(j_index,2) = 14.44309169
      d_mdrh(j_index,3) = -0.054841073
      d_mdrh(j_index,4) = 6.64259e-05


      j_index = 37
      d_mdrh(j_index,1) = 756.0747916
      d_mdrh(j_index,2) = -8.546826257
      d_mdrh(j_index,3) = 0.035798677
      d_mdrh(j_index,4) = -5.06629e-05


      j_index = 38
      d_mdrh(j_index,1) = 338.668191
      d_mdrh(j_index,2) = -2.971223403
      d_mdrh(j_index,3) = 0.012294866
      d_mdrh(j_index,4) = -1.87558e-05


      j_index = 39
      d_mdrh(j_index,1) = -53.18033508
      d_mdrh(j_index,2) = 0.663911748
      d_mdrh(j_index,3) = 9.16326e-4
      d_mdrh(j_index,4) = -6.70354e-06


      j_index = 40
      d_mdrh(j_index,1) = 3623.831129
      d_mdrh(j_index,2) = -39.27226457
      d_mdrh(j_index,3) = 0.144559515
      d_mdrh(j_index,4) = -1.78159e-4


      j_index = 41
      d_mdrh(j_index,1) = 3436.656743
      d_mdrh(j_index,2) = -37.16192684
      d_mdrh(j_index,3) = 0.136641377
      d_mdrh(j_index,4) = -1.68262e-4


      j_index = 42
      d_mdrh(j_index,1) = 768.608476
      d_mdrh(j_index,2) = -8.051517149
      d_mdrh(j_index,3) = 0.032342332
      d_mdrh(j_index,4) = -4.52224e-05


      j_index = 43
      d_mdrh(j_index,1) = 33.58027951
      d_mdrh(j_index,2) = -0.308772182
      d_mdrh(j_index,3) = 0.004713639
      d_mdrh(j_index,4) = -1.19658e-05


      j_index = 44
      d_mdrh(j_index,1) = 57.80183041
      d_mdrh(j_index,2) = 0.215264604
      d_mdrh(j_index,3) = 4.11406e-4
      d_mdrh(j_index,4) = -4.30702e-06


      j_index = 45
      d_mdrh(j_index,1) = -234.368984
      d_mdrh(j_index,2) = 2.721045204
      d_mdrh(j_index,3) = -0.006688341
      d_mdrh(j_index,4) = 2.31729e-06


      j_index = 46
      d_mdrh(j_index,1) = 3879.080557
      d_mdrh(j_index,2) = -42.13562874
      d_mdrh(j_index,3) = 0.155235005
      d_mdrh(j_index,4) = -1.91387e-4


      j_index = 47
      d_mdrh(j_index,1) = 3600.576985
      d_mdrh(j_index,2) = -39.0283489
      d_mdrh(j_index,3) = 0.143710316
      d_mdrh(j_index,4) = -1.77167e-4


      j_index = 48
      d_mdrh(j_index,1) = -1009.729826
      d_mdrh(j_index,2) = 12.9145339
      d_mdrh(j_index,3) = -0.049811146
      d_mdrh(j_index,4) = 6.09563e-05


      j_index = 49
      d_mdrh(j_index,1) = -577.0919514
      d_mdrh(j_index,2) = 8.020324227
      d_mdrh(j_index,3) = -0.031469556
      d_mdrh(j_index,4) = 3.82181e-05


      j_index = 50
      d_mdrh(j_index,1) = -728.9983499
      d_mdrh(j_index,2) = 9.849458215
      d_mdrh(j_index,3) = -0.03879257
      d_mdrh(j_index,4) = 4.78844e-05


      j_index = 51
      d_mdrh(j_index,1) = -803.7026845
      d_mdrh(j_index,2) = 10.61881494
      d_mdrh(j_index,3) = -0.041402993
      d_mdrh(j_index,4) = 5.08084e-05




      j_index = 52
      d_mdrh(j_index,1) = -493.6190458
      d_mdrh(j_index,2) = 6.747053851
      d_mdrh(j_index,3) = -0.026955267
      d_mdrh(j_index,4) = 3.45118e-05


      j_index = 53
      d_mdrh(j_index,1) = 53.37874093
      d_mdrh(j_index,2) = 1.01368249
      d_mdrh(j_index,3) = -0.005887513
      d_mdrh(j_index,4) = 8.94393e-06


      j_index = 54
      d_mdrh(j_index,1) = 206.619047
      d_mdrh(j_index,2) = -1.342735684
      d_mdrh(j_index,3) = 0.003197691
      d_mdrh(j_index,4) = -1.93603e-06


      j_index = 55
      d_mdrh(j_index,1) = -493.6190458
      d_mdrh(j_index,2) = 6.747053851
      d_mdrh(j_index,3) = -0.026955267
      d_mdrh(j_index,4) = 3.45118e-05


      j_index = 56
      d_mdrh(j_index,1) = 53.37874093
      d_mdrh(j_index,2) = 1.01368249
      d_mdrh(j_index,3) = -0.005887513
      d_mdrh(j_index,4) = 8.94393e-06


      j_index = 57
      d_mdrh(j_index,1) = 206.619047
      d_mdrh(j_index,2) = -1.342735684
      d_mdrh(j_index,3) = 0.003197691
      d_mdrh(j_index,4) = -1.93603e-06


      j_index = 58
      d_mdrh(j_index,1) = 41.7619047
      d_mdrh(j_index,2) = 1.303872053
      d_mdrh(j_index,3) = -0.007647908
      d_mdrh(j_index,4) = 1.17845e-05


      j_index = 59
      d_mdrh(j_index,1) = 41.7619047
      d_mdrh(j_index,2) = 1.303872053
      d_mdrh(j_index,3) = -0.007647908
      d_mdrh(j_index,4) = 1.17845e-05


      j_index = 60
      d_mdrh(j_index,1) = -369.7142842
      d_mdrh(j_index,2) = 5.512878771
      d_mdrh(j_index,3) = -0.02301948
      d_mdrh(j_index,4) = 3.0303e-05


      j_index = 61
      d_mdrh(j_index,1) = -369.7142842
      d_mdrh(j_index,2) = 5.512878771
      d_mdrh(j_index,3) = -0.02301948
      d_mdrh(j_index,4) = 3.0303e-05


      j_index = 62
      d_mdrh(j_index,1) = -162.8095232
      d_mdrh(j_index,2) = 2.399326592
      d_mdrh(j_index,3) = -0.009336219
      d_mdrh(j_index,4) = 1.17845e-05


      j_index = 63
      d_mdrh(j_index,1) = -735.4285689
      d_mdrh(j_index,2) = 8.885521857
      d_mdrh(j_index,3) = -0.033488456
      d_mdrh(j_index,4) = 4.12458e-05

      call load_kappa_nonelectro

      endif 

      return
      end subroutine load_mosaic_parameters




      subroutine load_kappa_nonelectro

      use module_data_mosaic_asect, only: &
         hygro_oin_aer, hygro_oc_aer, hygro_bc_aer,  &
         hygro_pcg1_b_c_aer,  hygro_pcg2_b_c_aer,  hygro_pcg3_b_c_aer,  &
         hygro_pcg4_b_c_aer,  hygro_pcg5_b_c_aer,  hygro_pcg6_b_c_aer,  &
         hygro_pcg7_b_c_aer,  hygro_pcg8_b_c_aer,  hygro_pcg9_b_c_aer,  &
         hygro_pcg1_b_o_aer,  hygro_pcg2_b_o_aer,  hygro_pcg3_b_o_aer,  &
         hygro_pcg4_b_o_aer,  hygro_pcg5_b_o_aer,  hygro_pcg6_b_o_aer,  &
         hygro_pcg7_b_o_aer,  hygro_pcg8_b_o_aer,  hygro_pcg9_b_o_aer,  &
         hygro_opcg1_b_c_aer, hygro_opcg2_b_c_aer, hygro_opcg3_b_c_aer,  &
         hygro_opcg4_b_c_aer, hygro_opcg5_b_c_aer, hygro_opcg6_b_c_aer,  &
         hygro_opcg7_b_c_aer, hygro_opcg8_b_c_aer,  &
         hygro_opcg1_b_o_aer, hygro_opcg2_b_o_aer, hygro_opcg3_b_o_aer,  &
         hygro_opcg4_b_o_aer, hygro_opcg5_b_o_aer, hygro_opcg6_b_o_aer,  &
         hygro_opcg7_b_o_aer, hygro_opcg8_b_o_aer,  &
         hygro_pcg1_f_c_aer,  hygro_pcg2_f_c_aer,  hygro_pcg3_f_c_aer,  &
         hygro_pcg4_f_c_aer,  hygro_pcg5_f_c_aer,  hygro_pcg6_f_c_aer,  &
         hygro_pcg7_f_c_aer,  hygro_pcg8_f_c_aer,  hygro_pcg9_f_c_aer,  &
         hygro_pcg1_f_o_aer,  hygro_pcg2_f_o_aer,  hygro_pcg3_f_o_aer,  &
         hygro_pcg4_f_o_aer,  hygro_pcg5_f_o_aer,  hygro_pcg6_f_o_aer,  &
         hygro_pcg7_f_o_aer,  hygro_pcg8_f_o_aer,  hygro_pcg9_f_o_aer,  &
         hygro_opcg1_f_c_aer, hygro_opcg2_f_c_aer, hygro_opcg3_f_c_aer,  &
         hygro_opcg4_f_c_aer, hygro_opcg5_f_c_aer, hygro_opcg6_f_c_aer,  &
         hygro_opcg7_f_c_aer, hygro_opcg8_f_c_aer,  &
         hygro_opcg1_f_o_aer, hygro_opcg2_f_o_aer, hygro_opcg3_f_o_aer,  &
         hygro_opcg4_f_o_aer, hygro_opcg5_f_o_aer, hygro_opcg6_f_o_aer,  &
         hygro_opcg7_f_o_aer, hygro_opcg8_f_o_aer,  &
         hygro_ant1_c_aer,  hygro_ant2_c_aer,  hygro_ant3_c_aer,  hygro_ant4_c_aer,  &
         hygro_ant1_o_aer,  hygro_ant2_o_aer,  hygro_ant3_o_aer,  hygro_ant4_o_aer,  &
         hygro_biog1_c_aer, hygro_biog2_c_aer, hygro_biog3_c_aer, hygro_biog4_c_aer,  &
         hygro_biog1_o_aer, hygro_biog2_o_aer, hygro_biog3_o_aer, hygro_biog4_o_aer,  &
         hygro_smpa_aer, hygro_smpbb_aer,  &
         hygro_glysoa_r1_aer,  hygro_glysoa_r2_aer,  hygro_glysoa_oh_aer,  &
         hygro_glysoa_nh4_aer, hygro_glysoa_sfc_aer,  &
         hygro_asoaX_aer, hygro_asoa1_aer, hygro_asoa2_aer,  &
         hygro_asoa3_aer, hygro_asoa4_aer,  &
         hygro_bsoaX_aer, hygro_bsoa1_aer, hygro_bsoa2_aer,  &
         hygro_bsoa3_aer, hygro_bsoa4_aer

      if (ioin_a        > 0) kappa_nonelectro(ioin_a       ) = hygro_oin_aer
      if (ioc_a         > 0) kappa_nonelectro(ioc_a        ) = hygro_oc_aer
      if (ibc_a         > 0) kappa_nonelectro(ibc_a        ) = hygro_bc_aer

      if (ipcg1_b_c_a   > 0) kappa_nonelectro(ipcg1_b_c_a  ) = hygro_pcg1_b_c_aer
      if (ipcg2_b_c_a   > 0) kappa_nonelectro(ipcg2_b_c_a  ) = hygro_pcg2_b_c_aer
      if (ipcg3_b_c_a   > 0) kappa_nonelectro(ipcg3_b_c_a  ) = hygro_pcg3_b_c_aer
      if (ipcg4_b_c_a   > 0) kappa_nonelectro(ipcg4_b_c_a  ) = hygro_pcg4_b_c_aer
      if (ipcg5_b_c_a   > 0) kappa_nonelectro(ipcg5_b_c_a  ) = hygro_pcg5_b_c_aer
      if (ipcg6_b_c_a   > 0) kappa_nonelectro(ipcg6_b_c_a  ) = hygro_pcg6_b_c_aer
      if (ipcg7_b_c_a   > 0) kappa_nonelectro(ipcg7_b_c_a  ) = hygro_pcg7_b_c_aer
      if (ipcg8_b_c_a   > 0) kappa_nonelectro(ipcg8_b_c_a  ) = hygro_pcg8_b_c_aer
      if (ipcg9_b_c_a   > 0) kappa_nonelectro(ipcg9_b_c_a  ) = hygro_pcg9_b_c_aer
      if (ipcg1_b_o_a   > 0) kappa_nonelectro(ipcg1_b_o_a  ) = hygro_pcg1_b_o_aer
      if (ipcg2_b_o_a   > 0) kappa_nonelectro(ipcg2_b_o_a  ) = hygro_pcg2_b_o_aer
      if (ipcg3_b_o_a   > 0) kappa_nonelectro(ipcg3_b_o_a  ) = hygro_pcg3_b_o_aer
      if (ipcg4_b_o_a   > 0) kappa_nonelectro(ipcg4_b_o_a  ) = hygro_pcg4_b_o_aer
      if (ipcg5_b_o_a   > 0) kappa_nonelectro(ipcg5_b_o_a  ) = hygro_pcg5_b_o_aer
      if (ipcg6_b_o_a   > 0) kappa_nonelectro(ipcg6_b_o_a  ) = hygro_pcg6_b_o_aer
      if (ipcg7_b_o_a   > 0) kappa_nonelectro(ipcg7_b_o_a  ) = hygro_pcg7_b_o_aer
      if (ipcg8_b_o_a   > 0) kappa_nonelectro(ipcg8_b_o_a  ) = hygro_pcg8_b_o_aer
      if (ipcg9_b_o_a   > 0) kappa_nonelectro(ipcg9_b_o_a  ) = hygro_pcg9_b_o_aer
      if (iopcg1_b_c_a  > 0) kappa_nonelectro(iopcg1_b_c_a ) = hygro_opcg1_b_c_aer
      if (iopcg2_b_c_a  > 0) kappa_nonelectro(iopcg2_b_c_a ) = hygro_opcg2_b_c_aer
      if (iopcg3_b_c_a  > 0) kappa_nonelectro(iopcg3_b_c_a ) = hygro_opcg3_b_c_aer
      if (iopcg4_b_c_a  > 0) kappa_nonelectro(iopcg4_b_c_a ) = hygro_opcg4_b_c_aer
      if (iopcg5_b_c_a  > 0) kappa_nonelectro(iopcg5_b_c_a ) = hygro_opcg5_b_c_aer
      if (iopcg6_b_c_a  > 0) kappa_nonelectro(iopcg6_b_c_a ) = hygro_opcg6_b_c_aer
      if (iopcg7_b_c_a  > 0) kappa_nonelectro(iopcg7_b_c_a ) = hygro_opcg7_b_c_aer
      if (iopcg8_b_c_a  > 0) kappa_nonelectro(iopcg8_b_c_a ) = hygro_opcg8_b_c_aer
      if (iopcg1_b_o_a  > 0) kappa_nonelectro(iopcg1_b_o_a ) = hygro_opcg1_b_o_aer
      if (iopcg2_b_o_a  > 0) kappa_nonelectro(iopcg2_b_o_a ) = hygro_opcg2_b_o_aer
      if (iopcg3_b_o_a  > 0) kappa_nonelectro(iopcg3_b_o_a ) = hygro_opcg3_b_o_aer
      if (iopcg4_b_o_a  > 0) kappa_nonelectro(iopcg4_b_o_a ) = hygro_opcg4_b_o_aer
      if (iopcg5_b_o_a  > 0) kappa_nonelectro(iopcg5_b_o_a ) = hygro_opcg5_b_o_aer
      if (iopcg6_b_o_a  > 0) kappa_nonelectro(iopcg6_b_o_a ) = hygro_opcg6_b_o_aer
      if (iopcg7_b_o_a  > 0) kappa_nonelectro(iopcg7_b_o_a ) = hygro_opcg7_b_o_aer
      if (iopcg8_b_o_a  > 0) kappa_nonelectro(iopcg8_b_o_a ) = hygro_opcg8_b_o_aer
      if (ipcg1_f_c_a   > 0) kappa_nonelectro(ipcg1_f_c_a  ) = hygro_pcg1_f_c_aer
      if (ipcg2_f_c_a   > 0) kappa_nonelectro(ipcg2_f_c_a  ) = hygro_pcg2_f_c_aer
      if (ipcg3_f_c_a   > 0) kappa_nonelectro(ipcg3_f_c_a  ) = hygro_pcg3_f_c_aer
      if (ipcg4_f_c_a   > 0) kappa_nonelectro(ipcg4_f_c_a  ) = hygro_pcg4_f_c_aer
      if (ipcg5_f_c_a   > 0) kappa_nonelectro(ipcg5_f_c_a  ) = hygro_pcg5_f_c_aer
      if (ipcg6_f_c_a   > 0) kappa_nonelectro(ipcg6_f_c_a  ) = hygro_pcg6_f_c_aer
      if (ipcg7_f_c_a   > 0) kappa_nonelectro(ipcg7_f_c_a  ) = hygro_pcg7_f_c_aer
      if (ipcg8_f_c_a   > 0) kappa_nonelectro(ipcg8_f_c_a  ) = hygro_pcg8_f_c_aer
      if (ipcg9_f_c_a   > 0) kappa_nonelectro(ipcg9_f_c_a  ) = hygro_pcg9_f_c_aer
      if (ipcg1_f_o_a   > 0) kappa_nonelectro(ipcg1_f_o_a  ) = hygro_pcg1_f_o_aer
      if (ipcg2_f_o_a   > 0) kappa_nonelectro(ipcg2_f_o_a  ) = hygro_pcg2_f_o_aer
      if (ipcg3_f_o_a   > 0) kappa_nonelectro(ipcg3_f_o_a  ) = hygro_pcg3_f_o_aer
      if (ipcg4_f_o_a   > 0) kappa_nonelectro(ipcg4_f_o_a  ) = hygro_pcg4_f_o_aer
      if (ipcg5_f_o_a   > 0) kappa_nonelectro(ipcg5_f_o_a  ) = hygro_pcg5_f_o_aer
      if (ipcg6_f_o_a   > 0) kappa_nonelectro(ipcg6_f_o_a  ) = hygro_pcg6_f_o_aer
      if (ipcg7_f_o_a   > 0) kappa_nonelectro(ipcg7_f_o_a  ) = hygro_pcg7_f_o_aer
      if (ipcg8_f_o_a   > 0) kappa_nonelectro(ipcg8_f_o_a  ) = hygro_pcg8_f_o_aer
      if (ipcg9_f_o_a   > 0) kappa_nonelectro(ipcg9_f_o_a  ) = hygro_pcg9_f_o_aer
      if (iopcg1_f_c_a  > 0) kappa_nonelectro(iopcg1_f_c_a ) = hygro_opcg1_f_c_aer
      if (iopcg2_f_c_a  > 0) kappa_nonelectro(iopcg2_f_c_a ) = hygro_opcg2_f_c_aer
      if (iopcg3_f_c_a  > 0) kappa_nonelectro(iopcg3_f_c_a ) = hygro_opcg3_f_c_aer
      if (iopcg4_f_c_a  > 0) kappa_nonelectro(iopcg4_f_c_a ) = hygro_opcg4_f_c_aer
      if (iopcg5_f_c_a  > 0) kappa_nonelectro(iopcg5_f_c_a ) = hygro_opcg5_f_c_aer
      if (iopcg6_f_c_a  > 0) kappa_nonelectro(iopcg6_f_c_a ) = hygro_opcg6_f_c_aer
      if (iopcg7_f_c_a  > 0) kappa_nonelectro(iopcg7_f_c_a ) = hygro_opcg7_f_c_aer
      if (iopcg8_f_c_a  > 0) kappa_nonelectro(iopcg8_f_c_a ) = hygro_opcg8_f_c_aer
      if (iopcg1_f_o_a  > 0) kappa_nonelectro(iopcg1_f_o_a ) = hygro_opcg1_f_o_aer
      if (iopcg2_f_o_a  > 0) kappa_nonelectro(iopcg2_f_o_a ) = hygro_opcg2_f_o_aer
      if (iopcg3_f_o_a  > 0) kappa_nonelectro(iopcg3_f_o_a ) = hygro_opcg3_f_o_aer
      if (iopcg4_f_o_a  > 0) kappa_nonelectro(iopcg4_f_o_a ) = hygro_opcg4_f_o_aer
      if (iopcg5_f_o_a  > 0) kappa_nonelectro(iopcg5_f_o_a ) = hygro_opcg5_f_o_aer
      if (iopcg6_f_o_a  > 0) kappa_nonelectro(iopcg6_f_o_a ) = hygro_opcg6_f_o_aer
      if (iopcg7_f_o_a  > 0) kappa_nonelectro(iopcg7_f_o_a ) = hygro_opcg7_f_o_aer
      if (iopcg8_f_o_a  > 0) kappa_nonelectro(iopcg8_f_o_a ) = hygro_opcg8_f_o_aer

      if (iant1_c_a     > 0) kappa_nonelectro(iant1_c_a    ) = hygro_ant1_c_aer
      if (iant2_c_a     > 0) kappa_nonelectro(iant2_c_a    ) = hygro_ant2_c_aer
      if (iant3_c_a     > 0) kappa_nonelectro(iant3_c_a    ) = hygro_ant3_c_aer
      if (iant4_c_a     > 0) kappa_nonelectro(iant4_c_a    ) = hygro_ant4_c_aer
      if (iant1_o_a     > 0) kappa_nonelectro(iant1_o_a    ) = hygro_ant1_o_aer
      if (iant2_o_a     > 0) kappa_nonelectro(iant2_o_a    ) = hygro_ant2_o_aer
      if (iant3_o_a     > 0) kappa_nonelectro(iant3_o_a    ) = hygro_ant3_o_aer
      if (iant4_o_a     > 0) kappa_nonelectro(iant4_o_a    ) = hygro_ant4_o_aer
      if (ibiog1_c_a    > 0) kappa_nonelectro(ibiog1_c_a   ) = hygro_biog1_c_aer
      if (ibiog2_c_a    > 0) kappa_nonelectro(ibiog2_c_a   ) = hygro_biog2_c_aer
      if (ibiog3_c_a    > 0) kappa_nonelectro(ibiog3_c_a   ) = hygro_biog3_c_aer
      if (ibiog4_c_a    > 0) kappa_nonelectro(ibiog4_c_a   ) = hygro_biog4_c_aer
      if (ibiog1_o_a    > 0) kappa_nonelectro(ibiog1_o_a   ) = hygro_biog1_o_aer
      if (ibiog2_o_a    > 0) kappa_nonelectro(ibiog2_o_a   ) = hygro_biog2_o_aer
      if (ibiog3_o_a    > 0) kappa_nonelectro(ibiog3_o_a   ) = hygro_biog3_o_aer
      if (ibiog4_o_a    > 0) kappa_nonelectro(ibiog4_o_a   ) = hygro_biog4_o_aer

      if (ismpa_a       > 0) kappa_nonelectro(ismpa_a      ) = hygro_smpa_aer
      if (ismpbb_a      > 0) kappa_nonelectro(ismpbb_a     ) = hygro_smpbb_aer
      if (iglysoa_r1_a  > 0) kappa_nonelectro(iglysoa_r1_a ) = hygro_glysoa_r1_aer
      if (iglysoa_r2_a  > 0) kappa_nonelectro(iglysoa_r2_a ) = hygro_glysoa_r2_aer
      if (iglysoa_oh_a  > 0) kappa_nonelectro(iglysoa_oh_a ) = hygro_glysoa_oh_aer
      if (iglysoa_nh4_a > 0) kappa_nonelectro(iglysoa_nh4_a) = hygro_glysoa_nh4_aer
      if (iglysoa_sfc_a > 0) kappa_nonelectro(iglysoa_sfc_a) = hygro_glysoa_sfc_aer
      if (iasoaX_a      > 0) kappa_nonelectro(iasoaX_a     ) = hygro_asoaX_aer
      if (iasoa1_a      > 0) kappa_nonelectro(iasoa1_a     ) = hygro_asoa1_aer
      if (iasoa2_a      > 0) kappa_nonelectro(iasoa2_a     ) = hygro_asoa2_aer
      if (iasoa3_a      > 0) kappa_nonelectro(iasoa3_a     ) = hygro_asoa3_aer
      if (iasoa4_a      > 0) kappa_nonelectro(iasoa4_a     ) = hygro_asoa4_aer
      if (ibsoaX_a      > 0) kappa_nonelectro(ibsoaX_a     ) = hygro_bsoaX_aer
      if (ibsoa1_a      > 0) kappa_nonelectro(ibsoa1_a     ) = hygro_bsoa1_aer
      if (ibsoa2_a      > 0) kappa_nonelectro(ibsoa2_a     ) = hygro_bsoa2_aer
      if (ibsoa3_a      > 0) kappa_nonelectro(ibsoa3_a     ) = hygro_bsoa3_aer
      if (ibsoa4_a      > 0) kappa_nonelectro(ibsoa4_a     ) = hygro_bsoa4_aer

      return
      end subroutine load_kappa_nonelectro










      subroutine update_thermodynamic_constants(vbs_nbin)



      integer iv, j_index, ibin, je,vbs_nbin(1)
      integer start_ind
      real(kind=8) :: tr, rt, term
      real(kind=8) :: gam_nh4no3_0, gam_nh4cl_0, m_nh4no3_0, m_nh4cl_0  




      tr = 298.15			
      rt = 82.056*t_k/(1.e9*1.e6)	


      keq_gl(1)= 1.0				         
      keq_gl(2)= fn_keq(57.64d0 , 13.79d0, -5.39d0,t_k)*rt     
      keq_gl(3)= fn_keq(2.63d6, 29.17d0, 16.83d0,t_k)*rt     
      keq_gl(4)= fn_keq(2.00d6, 30.20d0, 19.91d0,t_k)*rt     


      keq_ll(1)= fn_keq(1.0502d-2, 8.85d0, 25.14d0,t_k)      
      keq_ll(2)= fn_keq(1.805d-5, -1.50d0, 26.92d0,t_k)      
      keq_ll(3)= fn_keq(1.01d-14,-22.52d0, 26.92d0,t_k)      


      kp_nh3   = keq_ll(3)/(keq_ll(2)*keq_gl(2))
      kp_nh4no3= kp_nh3/keq_gl(3)
      kp_nh4cl = kp_nh3/keq_gl(4)



      keq_sg(1)= fn_keq(4.72d-17,-74.38d0,6.12d0,t_k)/rt**2  
      keq_sg(2)= fn_keq(8.43d-17,-71.00d0,2.40d0,t_k)/rt**2  



      keq_sl(jnh4so4) = fn_keq(1.040d0,-2.65d0, 38.57d0, t_k)  
      keq_sl(jlvcite) = fn_keq(11.8d0, -5.19d0, 54.40d0, t_k)  
      keq_sl(jnh4hso4)= fn_keq(117.0d0,-2.87d0, 15.83d0, t_k)  
      keq_sl(jnh4msa) = 1.e15				 
      keq_sl(jnh4no3) = fn_keq(12.21d0,-10.4d0, 17.56d0, t_k)  
      keq_sl(jnh4cl)  = fn_keq(17.37d0,-6.03d0, 16.92d0, t_k)  
      keq_sl(jna2so4) = fn_keq(0.491d0, 0.98d0, 39.75d0, t_k)  
      keq_sl(jnahso4) = fn_keq(313.0d0, 0.8d0,  14.79d0, t_k)  
      keq_sl(jna3hso4)= 1.e15		 	         
      keq_sl(jnamsa)  = 1.e15				 
      keq_sl(jnano3)  = fn_keq(11.95d0,-8.22d0, 16.01d0, t_k)  
      keq_sl(jnacl)   = fn_keq(38.28d0,-1.52d0, 16.89d0, t_k)  
      keq_sl(jcacl2)  = fn_keq(8.0d11,32.84d0,44.79d0, t_k)*1.e5  
      keq_sl(jcano3)  = fn_keq(4.31d5, 7.83d0,42.01d0, t_k)*1.e5  
      keq_sl(jcamsa2) = 1.e15				 

      start_ind = 1
      if (vbs_nbin(1).eq.0) then
        start_ind = ismpa_g
      else if (vbs_nbin(1) .eq. 4) then
        start_ind = iasoaX_g
      else
        start_ind = ipcg1_b_c_g
      endif
      
      do iv = start_ind, ngas_ioa + ngas_soa
        sat_soa(iv) = 0.0       
      enddo

       if (vbs_nbin(1).eq.9) then

      po_soa(ipcg1_b_c_g) = fn_po(9.91d-8, 112.0d0, T_K) 
      po_soa(ipcg2_b_c_g) = fn_po(9.91d-7, 106.0d0, T_K) 
      po_soa(ipcg3_b_c_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(ipcg4_b_c_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(ipcg5_b_c_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(ipcg6_b_c_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(ipcg7_b_c_g) = fn_po(9.91d-2, 76.0d0, T_K) 
      po_soa(ipcg8_b_c_g) = fn_po(9.91d-1, 70.0d0, T_K) 
      po_soa(ipcg9_b_c_g) = fn_po(9.91d0, 64.0d0, T_K) 
      po_soa(iopcg1_b_c_g) = fn_po(9.91d-8, 112.0d0, T_K) 
      po_soa(iopcg2_b_c_g) = fn_po(9.91d-7, 106.0d0, T_K) 
      po_soa(iopcg3_b_c_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(iopcg4_b_c_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(iopcg5_b_c_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(iopcg6_b_c_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(iopcg7_b_c_g) = fn_po(9.91d-2, 76.0d0, T_K) 
      po_soa(iopcg8_b_c_g) = fn_po(9.91d-1, 70.0d0, T_K) 
      po_soa(ipcg1_b_o_g) = fn_po(9.91d-8, 112.0d0, T_K) 
      po_soa(ipcg2_b_o_g) = fn_po(9.91d-7, 106.0d0, T_K) 
      po_soa(ipcg3_b_o_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(ipcg4_b_o_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(ipcg5_b_o_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(ipcg6_b_o_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(ipcg7_b_o_g) = fn_po(9.91d-2, 76.0d0, T_K) 
      po_soa(ipcg8_b_o_g) = fn_po(9.91d-1, 70.0d0, T_K) 
      po_soa(ipcg9_b_o_g) = fn_po(9.91d0, 64.0d0, T_K) 
      po_soa(iopcg1_b_o_g) = fn_po(9.91d-8, 112.0d0, T_K) 
      po_soa(iopcg2_b_o_g) = fn_po(9.91d-7, 106.0d0, T_K) 
      po_soa(iopcg3_b_o_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(iopcg4_b_o_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(iopcg5_b_o_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(iopcg6_b_o_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(iopcg7_b_o_g) = fn_po(9.91d-2, 76.0d0, T_K) 
      po_soa(iopcg8_b_o_g) = fn_po(9.91d-1, 70.0d0, T_K) 
      po_soa(ipcg1_f_c_g) = fn_po(9.91d-8, 112.0d0, T_K) 
      po_soa(ipcg2_f_c_g) = fn_po(9.91d-7, 106.0d0, T_K) 
      po_soa(ipcg3_f_c_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(ipcg4_f_c_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(ipcg5_f_c_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(ipcg6_f_c_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(ipcg7_f_c_g) = fn_po(9.91d-2, 76.0d0, T_K) 
      po_soa(ipcg8_f_c_g) = fn_po(9.91d-1, 70.0d0, T_K) 
      po_soa(ipcg9_f_c_g) = fn_po(9.91d0, 64.0d0, T_K) 
      po_soa(iopcg1_f_c_g) = fn_po(9.91d-8, 112.0d0, T_K) 
      po_soa(iopcg2_f_c_g) = fn_po(9.91d-7, 106.0d0, T_K) 
      po_soa(iopcg3_f_c_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(iopcg4_f_c_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(iopcg5_f_c_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(iopcg6_f_c_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(iopcg7_f_c_g) = fn_po(9.91d-2, 76.0d0, T_K) 
      po_soa(iopcg8_f_c_g) = fn_po(9.91d-1, 70.0d0, T_K) 
      po_soa(ipcg1_f_o_g) = fn_po(9.91d-8, 112.0d0, T_K) 
      po_soa(ipcg2_f_o_g) = fn_po(9.91d-7, 106.0d0, T_K) 
      po_soa(ipcg3_f_o_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(ipcg4_f_o_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(ipcg5_f_o_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(ipcg6_f_o_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(ipcg7_f_o_g) = fn_po(9.91d-2, 76.0d0, T_K) 
      po_soa(ipcg8_f_o_g) = fn_po(9.91d-1, 70.0d0, T_K) 
      po_soa(ipcg9_f_o_g) = fn_po(9.91d0, 64.0d0, T_K) 
      po_soa(iopcg1_f_o_g) = fn_po(9.91d-8, 112.0d0, T_K) 
      po_soa(iopcg2_f_o_g) = fn_po(9.91d-7, 106.0d0, T_K) 
      po_soa(iopcg3_f_o_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(iopcg4_f_o_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(iopcg5_f_o_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(iopcg6_f_o_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(iopcg7_f_o_g) = fn_po(9.91d-2, 76.0d0, T_K) 
      po_soa(iopcg8_f_o_g) = fn_po(9.91d-1, 70.0d0, T_K) 

      po_soa(iant1_c_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(iant2_c_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(iant3_c_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(iant4_c_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(iant1_o_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(iant2_o_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(iant3_o_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(iant4_o_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(ibiog1_c_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(ibiog2_c_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(ibiog3_c_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(ibiog4_c_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      po_soa(ibiog1_o_g) = fn_po(9.91d-6, 100.0d0, T_K) 
      po_soa(ibiog2_o_g) = fn_po(9.91d-5, 94.0d0, T_K) 
      po_soa(ibiog3_o_g) = fn_po(9.91d-4, 88.0d0, T_K) 
      po_soa(ibiog4_o_g) = fn_po(9.91d-3, 82.0d0, T_K) 
      endif

      if (vbs_nbin(1).eq.4) then
        po_soa(iasoaX_g) = fn_po(9.91d-10, 40.0d0, T_K) 
        po_soa(iasoa1_g) = fn_po(9.91d-6, dhr_approx(0.0d0), T_K) 
        po_soa(iasoa2_g) = fn_po(9.91d-5, dhr_approx(1.0d0), T_K) 
        po_soa(iasoa3_g) = fn_po(9.91d-4, dhr_approx(2.0d0), T_K) 
        po_soa(iasoa4_g) = fn_po(9.91d-3, dhr_approx(3.0d0), T_K) 
        po_soa(ibsoaX_g) = fn_po(9.91d-10, 40.0d0, T_K) 
        po_soa(ibsoa1_g) = fn_po(9.91d-6, dhr_approx(0.0d0), T_K) 
        po_soa(ibsoa2_g) = fn_po(9.91d-5, dhr_approx(1.0d0), T_K) 
        po_soa(ibsoa3_g) = fn_po(9.91d-4, dhr_approx(2.0d0), T_K) 
        po_soa(ibsoa4_g) = fn_po(9.91d-3, dhr_approx(3.0d0), T_K) 
      endif

      if (vbs_nbin(1).eq.3) then


        po_soa(ipcg1_b_c_g)  = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ipcg2_b_c_g)  = fn_po(9.91d-1, 83.0d0, T_K) 
        po_soa(iopcg1_b_c_g) = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ipcg1_b_o_g)  = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ipcg2_b_o_g)  = fn_po(9.91d-1, 83.0d0, T_K) 
        po_soa(iopcg1_b_o_g) = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ipcg1_f_c_g)  = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ipcg2_f_c_g)  = fn_po(9.91d-1, 83.0d0, T_K) 
        po_soa(iopcg1_f_c_g) = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ipcg1_f_o_g)  = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ipcg2_f_o_g)  = fn_po(9.91d-1, 83.0d0, T_K) 
        po_soa(iopcg1_f_o_g) = fn_po(9.91d-8, 83.0d0, T_K) 

        po_soa(iant1_c_g)    = fn_po(9.91d-7, 106.0d0, T_K) 
        po_soa(iant2_c_g)    = fn_po(9.91d-6, 100.0d0, T_K) 
        po_soa(iant3_c_g)    = fn_po(9.91d-5,  94.0d0, T_K) 
        po_soa(iant4_c_g)    = fn_po(9.91d-4,  88.0d0, T_K) 
        po_soa(ibiog1_c_g)   = fn_po(9.91d-7, 106.0d0, T_K) 
        po_soa(ibiog2_c_g)   = fn_po(9.91d-6, 100.0d0, T_K) 
        po_soa(ibiog3_c_g)   = fn_po(9.91d-5,  94.0d0, T_K) 
        po_soa(ibiog1_o_g)   = fn_po(9.91d-7, 106.0d0, T_K) 
        po_soa(ibiog2_o_g)   = fn_po(9.91d-6, 100.0d0, T_K) 
      endif

      if (vbs_nbin(1).eq.2) then
      po_soa(ipcg1_b_c_g) = fn_po(9.91d-8, 83.0d0, T_K) 
      po_soa(ipcg2_b_c_g) = fn_po(9.91d-1, 83.0d0, T_K) 
      po_soa(iopcg1_b_c_g) = fn_po(9.91d-8, 83.0d0, T_K) 
      po_soa(ipcg1_b_o_g) = fn_po(9.91d-8, 83.0d0, T_K) 
      po_soa(ipcg2_b_o_g) = fn_po(9.91d-1, 83.0d0, T_K) 
      po_soa(iopcg1_b_o_g) = fn_po(9.91d-8, 83.0d0, T_K) 
      po_soa(ipcg1_f_c_g) = fn_po(9.91d-8, 83.0d0, T_K) 
      po_soa(ipcg2_f_c_g) = fn_po(9.91d-1, 83.0d0, T_K) 
      po_soa(iopcg1_f_c_g) = fn_po(9.91d-8, 83.0d0, T_K) 
      po_soa(ipcg1_f_o_g) = fn_po(9.91d-8, 83.0d0, T_K) 
      po_soa(ipcg2_f_o_g) = fn_po(9.91d-1, 83.0d0, T_K) 
      po_soa(iopcg1_f_o_g) = fn_po(9.91d-8, 83.0d0, T_K) 
      po_soa(iant1_c_g) = fn_po(9.91d-6, 83.0d0, T_K) 
      po_soa(iant1_o_g) = fn_po(9.91d-6, 83.0d0, T_K) 
      po_soa(ibiog1_c_g) = fn_po(9.91d-6, 83.0d0, T_K) 
      po_soa(ibiog1_o_g) = fn_po(9.91d-6, 83.0d0, T_K) 
      endif
      if (vbs_nbin(1).eq.0) then
        po_soa(ismpa_g) = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ismpbb_g) = fn_po(9.91d-8, 83.0d0, T_K) 
        po_soa(ibiog1_c_g) = fn_po(9.91d-6, 83.0d0, T_K) 
        po_soa(ibiog1_o_g) = fn_po(9.91d-6, 83.0d0, T_K) 
      endif

      start_ind = 1
      if (vbs_nbin(1).eq.0) then
        start_ind = ismpa_g
      else if (vbs_nbin(1).eq.4) then
        start_ind = iasoaX_g
      else
        start_ind = ipcg1_b_c_g
      end if

      do iv = start_ind, ngas_ioa + ngas_soa
        sat_soa(iv) = 1.e9*po_soa(iv)/(8.314*t_k)	
      enddo


      term = (647.15 - t_k)/647.15
      sigma_water = 0.2358*term**1.256 * (1. - 0.625*term) 


      do j_index = 1, 63
        mdrh_t(j_index) = drh_mutual(j_index)
      enddo




      do ibin = 1, nbin_a
        ah2o_a(ibin) = ah2o			
      enddo

      call mtem_compute_log_gamz		


      gam_nh4no3_0 = 10.**log_gamZ(jnh4no3,jnh4no3)
      gam_nh4cl_0  = 10.**log_gamZ(jnh4cl,jnh4cl)

      m_nh4no3_0   = molality_0(jnh4no3)
      m_nh4cl_0    = molality_0(jnh4cl)

      Kp_nh4no3_0  = Kp_nh4no3*(m_nh4no3_0*gam_nh4no3_0)**2
      Kp_nh4cl_0   = Kp_nh4cl *(m_nh4cl_0 *gam_nh4cl_0 )**2




      return
      end subroutine update_thermodynamic_constants

      
      
      
      
      real(kind=8) function dhr_approx(log10_Csat_298)

        real(kind=8), intent(in) :: log10_Csat_298

        dhr_approx = -11.0 * log10_Csat_298 + 131.0 

      end function dhr_approx













      real(kind=8) function fn_keq(keq_298, a, b, t)


      real(kind=8) keq_298, a, b, t

      real(kind=8) tt


        tt = 298.15/t
        fn_keq = keq_298*exp(a*(tt-1.)+b*(1.+log(tt)-tt))

      return
      end function fn_keq







      real(kind=8) function fn_po(po_298, dh, t)	


      real(kind=8) po_298, dh, t


        fn_po = po_298*exp(-(dh/8.314e-3)*(1./t - 3.354016435e-3))

      return
      end function fn_po







      real(kind=8) function drh_mutual(j_index)



      integer j_index

      integer j


      j = j_index

      if(j_index .eq. 7 .or. j_index .eq. 8 .or.   &
        (j_index.ge. 34 .and. j_index .le. 51))then

        drh_mutual = 10.0  

      else

        drh_mutual =  d_mdrh(j,1) + t_k*   &
                     (d_mdrh(j,2) + t_k*   &
                     (d_mdrh(j,3) + t_k*   &
                      d_mdrh(j,4) )) + 1.0

      endif


      return
      end function drh_mutual










      real(kind=8) function aerosol_water_up(ibin) 



      integer ibin

      integer jp, je
      real(kind=8) dum




      jp = jtotal
      dum = 0.0

      do je = 1, (nsalt+4)	
        dum = dum + 1.e-9*electrolyte(je,jp,ibin)/bin_molality_60(je)
      enddo

      aerosol_water_up = dum

      return
      end function aerosol_water_up












      real(kind=8) function aerosol_water(jp,ibin) 



      integer jp, ibin

      integer ja, je
      real(kind=8) dum, tmpa




      dum = 0.0
      do je = 1, (nsalt+4)	
        dum = dum + electrolyte(je,jp,ibin)/bin_molality(je,ibin)
      enddo

      if (mwater_kappa_nonelectro > 0) then
         tmpa = 0.0
         do ja = 1, naer
            if (kappa_nonelectro(ja) > 0.0) then
               tmpa = tmpa + (aer(ja,jtotal,ibin)*mw_aer_mac(ja)/dens_aer_mac(ja))*kappa_nonelectro(ja)
            end if
         end do
         dum = dum + 1.0e-3*tmpa*aH2O_a(ibin)/(1.0-aH2O_a(ibin))
      end if

      aerosol_water = dum*1.e-9  

      if(aerosol_water .le. 0.0)then
        if (iprint_mosaic_diag1 .gt. 0) then
          write(6,*)'mosaic aerosol_water - water .le. 0'
          write(6,*)'iclm  jclm  ibin  jp = ',   &
                     iclm_aer, jclm_aer, ibin, jp
          write(6,*)'ah2o, water = ', ah2o, aerosol_water
          write(6,*)'dry mass = ', mass_dry_a(ibin)
          write(6,*)'soluble mass = ', mass_soluble_a(ibin)
          write(6,*)'number = ', num_a(ibin)
          do je = 1, nsoluble
            write(6,44)ename(je), electrolyte(je,jp,ibin)
          enddo
          write(6,*)'error in water calculation'
          write(6,*)'ibin = ', ibin
          write(6,*)'water content cannot be negative or zero'
          write(6,*)'setting jaerosolstate to all_solid'
        endif

        call print_input

        jaerosolstate(ibin) = all_solid
        jphase(ibin)    = jsolid
        jhyst_leg(ibin) = jhyst_lo



      endif

44    format(a7, 2x, e11.3)


      return
      end function aerosol_water







      real(kind=8) function bin_molality(je,ibin)



      integer je, ibin

      real(kind=8) aw, xm


      aw = max(ah2o_a(ibin), aw_min(je))
      aw = min(aw, 0.999999D0)


      if(aw .lt. 0.97)then

        xm =     a_zsr(1,je) +   &
             aw*(a_zsr(2,je) +   &
             aw*(a_zsr(3,je) +   &
             aw*(a_zsr(4,je) +   &
             aw*(a_zsr(5,je) +   &
             aw* a_zsr(6,je) ))))

        bin_molality = 55.509*xm/(1. - xm)

      else

        bin_molality = -b_zsr(je)*log(aw)

      endif


      return
      end function bin_molality







      real(kind=8) function bin_molality_60(je)



      integer je

      real(kind=8) aw, xm


      aw = 0.6

        xm =  a_zsr(1,je) + aw*   &
             (a_zsr(2,je) + aw*   &
             (a_zsr(3,je) + aw*   &
             (a_zsr(4,je) + aw*   &
             (a_zsr(5,je) + aw*   &
              a_zsr(6,je) ))))

      bin_molality_60 = 55.509*xm/(1. - xm)

      return
      end function bin_molality_60





      real(kind=8) function molality_0(je)


      integer je

      real(kind=8) :: aw, xm


      aw = max(ah2o, aw_min(je))
      aw = min(aw, 0.999999d0)


      if(aw .lt. 0.97)then

        xm =     a_zsr(1,je) +   &
             aw*(a_zsr(2,je) +   &
             aw*(a_zsr(3,je) +   &
             aw*(a_zsr(4,je) +   &
             aw*(a_zsr(5,je) +   &
             aw* a_zsr(6,je) ))))

        molality_0 = 55.509*xm/(1. - xm)

      else

        molality_0 = -b_zsr(je)*log(aw)

      endif


      return
      end function molality_0





      real(kind=8) function fnlog_gamz(ja,je)	



      integer ja, je

      real(kind=8) aw


      aw = max(ah2o, aw_min(je))

      fnlog_gamz = b_mtem(1,ja,je) + aw*   &
                  (b_mtem(2,ja,je) + aw*   &
                  (b_mtem(3,ja,je) + aw*   &
                  (b_mtem(4,ja,je) + aw*   &
                  (b_mtem(5,ja,je) + aw*   &
                   b_mtem(6,ja,je) ))))

      return
      end function fnlog_gamz






      real(kind=8) function mean_molecular_speed(t, mw)	


      real(kind=8) t, mw	

        mean_molecular_speed = 1.455e4 * sqrt(t/mw)

      return
      end function mean_molecular_speed






      real(kind=8) function gas_diffusivity(t, p, mw, vm)	


      real(kind=8) mw, vm, t, p	


      gas_diffusivity = (1.0e-3 * t**1.75 * sqrt(1./mw + 0.035))/   &
                             (p * (vm**0.333333 + 2.7189)**2)


      return
      end function gas_diffusivity






      real(kind=8) function fuchs_sutugin(rkn,a)


      real(kind=8) rkn, a

      real(kind=8) rnum, denom


      rnum  = 0.75*a*(1. + rkn)
      denom = rkn**2 + rkn + 0.283*rkn*a + 0.75*a
      fuchs_sutugin = rnum/denom

      return
      end function fuchs_sutugin





    real(kind=8) function acc_n2o5_bert_thorn(mass_h2o,mol_no3,mol_cl,vol)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	real(kind=8), parameter :: A_bt  = 3.2e-8
	real(kind=8), parameter :: b_bt  = 1.15e6
	real(kind=8), parameter :: d_bt  = 1.3e-1
	real(kind=8), parameter :: k3_bt = 6.0e-2
	real(kind=8), parameter :: k4_bt = 29e0

	
	real(kind=8), parameter :: nmol_mol = 1e-9	
	real(kind=8), parameter :: m3_litre = 1e3	
	real(kind=8), parameter :: mm_h2o   = 18e-3	

	
	real(kind=8) :: mass_h2o	
	real(kind=8) :: mol_no3		
	real(kind=8) :: mol_cl		
	real(kind=8) :: vol			

	
	real(kind=8) :: part_step
	real(kind=8) :: aer_h2o, aer_no3, aer_cl	


	
	aer_h2o = mass_h2o / (mm_h2o*vol*m3_litre) 
	aer_no3 = mol_no3*nmol_mol / (vol*m3_litre)
	aer_cl  = mol_cl*nmol_mol / (vol*m3_litre)
	
	if(n2o5_flag.eq.1)then 
		aer_cl = 0.0
	end if
	
	if(aer_h2o .ne. 0.0)then
		part_step =  b_bt - b_bt * exp(-d_bt*aer_h2o)
		if(aer_no3 .ne. 0.0)then
			acc_n2o5_bert_thorn = A_bt * part_step *  &		
					(1.0 - 1.0 / (                    &
						1.0 +                         &
						(k3_bt*aer_h2o/aer_no3) +     &
						(k4_bt*aer_cl/aer_no3)        &
					))
		else
			acc_n2o5_bert_thorn = A_bt * part_step
		endif
	else 
		acc_n2o5_bert_thorn = 0.0
	endif

	return
	end function acc_n2o5_bert_thorn




	real(kind=8) function split_n2o5_bert_thorn(mass_h2o,mol_cl,vol)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	real(kind=8), parameter :: k3_bt = 6.0e-2
	real(kind=8), parameter :: k4_bt = 29e0

	
	real(kind=8), parameter :: nmol_mol = 1e-9	
	real(kind=8), parameter :: m3_litre = 1e3	
	real(kind=8), parameter :: mm_h2o   = 18e-3	

	
	real(kind=8) :: mass_h2o	
	real(kind=8) :: mol_cl		
	real(kind=8) :: vol			

	
	real(kind=8) :: part_step
	real(kind=8) :: aer_h2o, aer_cl	


	
	aer_h2o = mass_h2o / (mm_h2o*vol*m3_litre) 
	aer_cl  = mol_cl*nmol_mol / (vol*m3_litre)

	if(n2o5_flag.eq.1)then 
		aer_cl = 0.0
	end if

	if(aer_h2o .ne. 0.0)then
		split_n2o5_bert_thorn = 1e0 / &
					( 1e0 + (k4_bt*aer_cl)/(k3_bt*aer_h2o) )	
	else
		split_n2o5_bert_thorn = 0.0
	endif



	return
	end function split_n2o5_bert_thorn










      real(kind=8) function cubic( p, q, r )


      real(kind=8), intent(in) :: p, q, r

      real(kind=8) a, b, d, m, n, third, y
      real(kind=8) k, phi, thesign, x(3), duma
      integer icase, kk

      third = 1.d0/3.d0

      a = (1.d0/3.d0)*((3.d0*q) - (p*p))
      b = (1.d0/27.d0)*((2.d0*p*p*p) - (9.d0*p*q) + (27.d0*r))

      d = ( ((a*a*a)/27.d0) + ((b*b)/4.d0) )

      if(d .gt. 0.)then	
        icase = 1
      elseif(d .eq. 0.)then 
        icase = 2
      else	
        icase = 3
      endif


      goto (1,2,3), icase


1     thesign = 1.
      if(b .gt. 0.)then
        b = -b
        thesign = -1.
      endif

      m = thesign*((-b/2.d0) + (sqrt(d)))**(third)
      n = thesign*((-b/2.d0) - (sqrt(d)))**(third)

      cubic = real( (m) + (n) - (p/3.d0) )
      return


2     thesign = 1.
      if(b .gt. 0.)then
        b = -b
        thesign = -1.
      endif

      m = thesign*(-b/2.d0)**third
      n = m

      x(1) = real( (m) + (n) - (p/3.d0) )
      x(2) = real( (-m/2.d0) + (-n/2.d0) - (p/3.d0) )
      x(2) = real( (-m/2.d0) + (-n/2.d0) - (p/3.d0) )

      cubic = 0.
      do kk = 1, 3
        if(x(kk).gt.cubic) cubic = x(kk)
      enddo
      return


3     if(b.gt.0.)then
        thesign = -1.
      elseif(b.lt.0.)then
        thesign = 1.
      endif



      duma = thesign*sqrt( (b*b/4.d0)/(-a*a*a/27.d0) )
      duma = min( duma, +1.0D0 )
      duma = max( duma, -1.0D0 )
      phi  = acos( duma )	


      cubic = 0.
      do kk = 1, 3
        k = kk-1
        y = 2.*sqrt(-a/3.)*cos(phi + 120.*k*0.017453293)
        x(kk) = real((y) - (p/3.d0))
        if(x(kk).gt.cubic) cubic = x(kk)
      enddo
      return

      end function cubic






      real(kind=8) function quadratic(a,b,c)


      real(kind=8) a, b, c

      real(kind=8) x, dum, quad1, quad2


        if(b .ne. 0.0)then
        x = 4.*(a/b)*(c/b)
        else
        x = 1.e+6
        endif

        if(abs(x) .lt. 1.e-6)then
          dum = (0.5*x) +   &
                (0.125*x**2) +   &
                (0.0625*x**3)

          quadratic = (-0.5*b/a)*dum

          if(quadratic .lt. 0.)then
            quadratic = -b/a - quadratic
          endif

        else
          quad1 = (-b+sqrt(b*b-4.*a*c))/(2.*a)
          quad2 = (-b-sqrt(b*b-4.*a*c))/(2.*a)

          quadratic = max(quad1, quad2)
        endif

      return
      end function quadratic








 
      subroutine quadratix(a,b,c, qx1,qx2)


      real(kind=8) a, b, c, qx1, qx2

      real(kind=8) x, dum


      if(b .ne. 0.0)then
        x = 4.*(a/b)*(c/b)
        else
        x = 1.e+6
      endif

      if(abs(x) .lt. 1.e-6)then
        dum = (0.5*x) +   &
              (0.125*x**2) +   &
              (0.0625*x**3)

        qx1 = (-0.5*b/a)*dum
        qx2 = -b/a - qx1

      else

        qx1 = (-b+sqrt(b*b - 4.*a*c))/(2.*a)
        qx2 = (-b-sqrt(b*b - 4.*a*c))/(2.*a)

      endif

      return
      end subroutine quadratix





























      subroutine save_pregrow_props

      use module_data_mosaic_asect
      use module_data_mosaic_other










      integer ibin, isize, itype



      cair_mol_cc = cairclm(kclm_aer)


      do ibin = 1, nbin_a

      call calc_dry_n_wet_aerosol_props( ibin )

      call isize_itype_from_ibin( ibin, isize, itype )
      drymass_pregrow(isize,itype) = mass_dry_a(ibin)/cair_mol_cc	
      if(jaerosolstate(ibin) .eq. no_aerosol) then
          drydens_pregrow(isize,itype) = -1.
      else
          drydens_pregrow(isize,itype) = dens_dry_a(ibin)		
      end if

      end do

      return
      end subroutine save_pregrow_props












	subroutine specialoutaa( iclm, jclm, kclm, msub, fromwhere )



	integer iclm, jclm, kclm, msub
	character*(*) fromwhere

	return
	end subroutine specialoutaa









	subroutine aerchem_boxtest_output(   &
      		iflag, iclm, jclm, kclm, msub, dtchem )

	use module_data_mosaic_asect
	use module_data_mosaic_other






	integer iflag, iclm, jclm, kclm, msub
	real(kind=8) dtchem


	integer lun
	parameter (lun=83)
	integer, save :: ientryno = -13579
	integer icomp, iphase, isize, itype, k, l, m, n

	real(kind=8) dtchem_sv1
	save dtchem_sv1
	real(kind=8) rsub_sv1(l2maxd,kmaxd,nsubareamaxd)



	if (maerchem_boxtest_output .le. 0) return






	itype = 1
	iphase = ai_phase


	if (ientryno .ne. -13579) goto 1000

	ientryno = +1
	call peg_message( lunerr, '***' )
	call peg_message( lunerr, '*** doing initial aerchem_boxtest_output' )
	call peg_message( lunerr, '***' )

	write(lun) ltot, ltot2, itot, jtot, ktot
	write(lun) (name(l), l=1,ltot2)

	write(lun) maerocoag, maerchem, maeroptical
	write(lun) msectional, maerosolincw

	write(lun) nsize_aer(itype), ntot_mastercomp_aer

	do icomp = 1, ntot_mastercomp_aer
	    write(lun)   &
      		name_mastercomp_aer(icomp)
	    write(lun)   &
      		dens_mastercomp_aer(icomp),     mw_mastercomp_aer(icomp)
	end do

	do isize = 1, nsize_aer(itype)
	    write(lun)   &
      		ncomp_plustracer_aer(itype),   &
		ncomp_aer(itype),   &
      		waterptr_aer(isize,itype),   &
		numptr_aer(isize,itype,iphase),   &
      		mprognum_aer(isize,itype,iphase)
	    write(lun)   &
      	      ( mastercompptr_aer(l,itype),   &
		massptr_aer(l,isize,itype,iphase),   &
      		l=1,ncomp_plustracer_aer(itype) )
	    write(lun)   &
      		volumcen_sect(isize,itype),   &
		volumlo_sect(isize,itype),   &
      		volumhi_sect(isize,itype),   &
		dcen_sect(isize,itype),   &
      		dlo_sect(isize,itype),   &
		dhi_sect(isize,itype)
	    write(lun)   &
      		lptr_so4_aer(isize,itype,iphase),   &
      		lptr_msa_aer(isize,itype,iphase),   &
      		lptr_no3_aer(isize,itype,iphase),   &
      		lptr_cl_aer(isize,itype,iphase),   &
      		lptr_co3_aer(isize,itype,iphase),   &
      		lptr_nh4_aer(isize,itype,iphase),   &
      		lptr_na_aer(isize,itype,iphase),   &
      		lptr_ca_aer(isize,itype,iphase),   &
      		lptr_oin_aer(isize,itype,iphase),   &
      		lptr_oc_aer(isize,itype,iphase),   &
      		lptr_bc_aer(isize,itype,iphase),   &
      		hyswptr_aer(isize,itype)
	end do




1000	continue
	if (iflag .eq. 1) goto 1010
	if (iflag .eq. 2) goto 2000
	if (iflag .eq. 3) goto 3000
	return




1010	continue
	dtchem_sv1 = dtchem
	do m = 1, nsubareas
	do k = 1, ktot
	do l = 1, ltot2
	    rsub_sv1(l,k,m) = rsub(l,k,m)
	end do
	end do
	end do

	return





2000	continue
	return





3000	continue
	do m = 1, nsubareas
	do k = 1, ktot

	write(lun) iymdcur, ihmscur, iclm, jclm, k, m, nsubareas
	write(lun) t, dtchem_sv1, cairclm(k), relhumclm(k),   &
      		ptotclm(k), afracsubarea(k,m)

	write(lun) (rsub_sv1(l,k,m), rsub(l,k,m), l=1,ltot2)

	end do
	end do


	return
	end subroutine aerchem_boxtest_output








	subroutine mosaic_aerchem_error_dump( istop, ibin, luna, msga )




	use module_data_mosaic_asect
	use module_data_mosaic_other



	integer istop, ibin, luna
	character*(*) msga


	integer icomp, iphase, isize, itype, k, l, lunb, m, n
	real(kind=8) dtchem_sv1





	itype = 1


	lunb = luna
	if (lunb .le. 0) lunb = 6

9000	format( a )
9010	format( 7i10 )
9020	format( 3(1pe19.11) )

	write(lunb,9000)
	write(lunb,9000) 'begin mosaic_aerchem_error_dump - msga ='
	write(lunb,9000) msga
	write(lunb,9000) 'i, j, k, msub,ibin ='
	write(lunb,9010) iclm_aer, jclm_aer, kclm_aer, mclm_aer, ibin

	write(lunb,9010) ltot, ltot2, itot, jtot, ktot
	write(lunb,9000) (name(l), l=1,ltot2)

	write(lunb,9010) maerocoag, maerchem, maeroptical
	write(lunb,9010) msectional, maerosolincw

	write(lunb,9010) nsize_aer(itype), ntot_mastercomp_aer

	do icomp = 1, ntot_mastercomp_aer
	    write(lunb,9000)   &
      		name_mastercomp_aer(icomp)
	    write(lunb,9020)   &
      		dens_mastercomp_aer(icomp),     mw_mastercomp_aer(icomp)
	end do

	do isize = 1, nsize_aer(itype)
	    write(lunb,9010)   &
      		ncomp_plustracer_aer(itype),   &
		ncomp_aer(itype),   &
      		waterptr_aer(isize,itype),   &
		numptr_aer(isize,itype,iphase),   &
      		mprognum_aer(isize,itype,iphase)
	    write(lunb,9010)   &
      	      ( mastercompptr_aer(l,itype),   &
		massptr_aer(l,isize,itype,iphase),   &
      		l=1,ncomp_plustracer_aer(itype) )
	    write(lunb,9020)   &
      		volumcen_sect(isize,itype),   &
		volumlo_sect(isize,itype),   &
      		volumhi_sect(isize,itype),   &
		dcen_sect(isize,itype),   &
      		dlo_sect(isize,itype),   &
		dhi_sect(isize,itype)
	    write(lunb,9010)   &
      		lptr_so4_aer(isize,itype,iphase),   &
      		lptr_msa_aer(isize,itype,iphase),   &
      		lptr_no3_aer(isize,itype,iphase),   &
      		lptr_cl_aer(isize,itype,iphase),   &
      		lptr_co3_aer(isize,itype,iphase),   &
      		lptr_nh4_aer(isize,itype,iphase),   &
      		lptr_na_aer(isize,itype,iphase),   &
      		lptr_ca_aer(isize,itype,iphase),   &
      		lptr_oin_aer(isize,itype,iphase),   &
      		lptr_oc_aer(isize,itype,iphase),   &
      		lptr_bc_aer(isize,itype,iphase),   &
      		hyswptr_aer(isize,itype)
	end do


	dtchem_sv1 = -1.0
	do m = 1, nsubareas
	do k = 1, ktot

	write(lunb,9010) iymdcur, ihmscur, iclm_aer, jclm_aer, k, m, nsubareas
	write(lunb,9020) t, dtchem_sv1, cairclm(k), relhumclm(k),   &
      		ptotclm(k), afracsubarea(k,m)

	write(lunb,9020) (rsub(l,k,m), l=1,ltot2)

	end do
	end do

	write(lunb,9000) 'end mosaic_aerchem_error_dump'


	if (istop .gt. 0) call peg_error_fatal( luna, msga )

	return
	end subroutine mosaic_aerchem_error_dump


      end module module_mosaic_therm

