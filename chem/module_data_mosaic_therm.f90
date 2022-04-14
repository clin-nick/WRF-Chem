

















      module module_data_mosaic_therm


      use module_data_mosaic_aero, only:  nbin_a_max, nbin_a


      implicit none



















      integer, parameter :: nbin_a_maxd = nbin_a_max







      integer ngas_ioa, ngas_soa, ngas_volatile, ngas_het,	&
              naer, naercomp, nelectrolyte, nsalt,	&
              nsoluble, ncation, nanion
      parameter(ngas_ioa = 5)	






      parameter(ngas_soa = 68+2+16+10)	

      parameter(ngas_volatile = ngas_ioa + ngas_soa + 1 + 1)
      parameter(ngas_het = 2)   

      parameter(naer = 11+68+2+16+5+10)	
      parameter(naercomp = 26+68+2+16+5+10)	
      parameter(nelectrolyte = 22) 
      parameter(nsalt    = 15)	
      parameter(nsoluble = 20)	
      parameter(ncation = 4)	
      parameter(nanion  = 5)	

      integer nrxn_aer_gl, nrxn_aer_ll, nrxn_aer_sg, nrxn_aer_sl
      parameter(nrxn_aer_gl = 4) 
      parameter(nrxn_aer_ll = 3) 
      parameter(nrxn_aer_sg = 2) 
      parameter(nrxn_aer_sl = nsalt)

      integer mmodal, msection,   &
              mon, moff, myes, mno
      parameter(mmodal  = 1)	
      parameter(msection= 2)	
      parameter(mon     = 1)	
      parameter(moff    = 0)    
      parameter(myes	= mon)	
      parameter(mno	= moff)	


      integer jtotal, jsolid, jliquid
      parameter(jsolid = 1)
      parameter(jliquid= 2)
      parameter(jtotal = 3)

      integer jhyst_lo, jhyst_up
      parameter(jhyst_lo = 0)	
      parameter(jhyst_up = 1) 	

  
      
      integer, parameter :: mhyst_uporlo_waterhyst = 2	
      integer, parameter :: mhyst_force_up = 3	        
      integer, parameter :: mhyst_force_lo = 4	        
      integer, parameter :: mhyst_method = mhyst_uporlo_waterhyst
  
      real(kind=8), parameter :: xhyst_up_crustal_thresh = 0.30
      
      
      
      
      
  
      integer, parameter :: mwater_kappa_nonelectro = 1  
  
      integer no_aerosol, all_solid, all_liquid, mixed
      parameter(no_aerosol = 0)	
      parameter(all_solid  = 1) 
      parameter(all_liquid = 2) 
      parameter(mixed      = 3)	

      integer soluble, insoluble
      parameter(soluble   = 1)  
      parameter(insoluble = 2)  

      real(kind=8) mass_cutoff
      parameter(mass_cutoff = 1.d-15)	






      integer, save ::   &
       ih2so4_g,     ihno3_g,      ihcl_g,      inh3_g,   &
       imsa_g, in2o5_g, iclno2_g

      integer, save ::   &
      ipcg1_b_c_g,ipcg2_b_c_g,ipcg3_b_c_g,ipcg4_b_c_g, &
      ipcg5_b_c_g,ipcg6_b_c_g,ipcg7_b_c_g,ipcg8_b_c_g, &
      ipcg9_b_c_g,ipcg1_b_o_g,ipcg2_b_o_g,ipcg3_b_o_g, &
      ipcg4_b_o_g,ipcg5_b_o_g,ipcg6_b_o_g,ipcg7_b_o_g, &
      ipcg8_b_o_g,ipcg9_b_o_g,iopcg1_b_c_g,iopcg2_b_c_g,&
      iopcg3_b_c_g, iopcg4_b_c_g,iopcg5_b_c_g,iopcg6_b_c_g,&
      iopcg7_b_c_g,iopcg8_b_c_g,iopcg1_b_o_g,iopcg2_b_o_g,&
      iopcg3_b_o_g,iopcg4_b_o_g,iopcg5_b_o_g,iopcg6_b_o_g,&
      iopcg7_b_o_g,iopcg8_b_o_g,&
      ipcg1_f_c_g,ipcg2_f_c_g,ipcg3_f_c_g,ipcg4_f_c_g, &
      ipcg5_f_c_g,ipcg6_f_c_g,ipcg7_f_c_g,ipcg8_f_c_g, &
      ipcg9_f_c_g,ipcg1_f_o_g,ipcg2_f_o_g,ipcg3_f_o_g, &
      ipcg4_f_o_g,ipcg5_f_o_g,ipcg6_f_o_g,ipcg7_f_o_g, &
      ipcg8_f_o_g,ipcg9_f_o_g,iopcg1_f_c_g,iopcg2_f_c_g,&
      iopcg3_f_c_g, iopcg4_f_c_g,iopcg5_f_c_g,iopcg6_f_c_g,&
      iopcg7_f_c_g,iopcg8_f_c_g,iopcg1_f_o_g,iopcg2_f_o_g,&
      iopcg3_f_o_g,iopcg4_f_o_g,iopcg5_f_o_g,iopcg6_f_o_g,&
      iopcg7_f_o_g,iopcg8_f_o_g,iant1_c_g,iant2_c_g,iant3_c_g, &
      iant4_c_g,ibiog1_c_g,ibiog2_c_g,ibiog3_c_g,ibiog4_c_g, &
      iant1_o_g,iant2_o_g,iant3_o_g, &
      iant4_o_g,ibiog1_o_g,ibiog2_o_g,ibiog3_o_g,ibiog4_o_g, &
      ismpa_g,ismpbb_g, &
      iasoaX_g, iasoa1_g, iasoa2_g, iasoa3_g, iasoa4_g, &
      ibsoaX_g, ibsoa1_g, ibsoa2_g, ibsoa3_g, ibsoa4_g, &
      igly, iho



      integer, save ::   &
       iso4_a,     ino3_a,     icl_a,     inh4_a,     ico3_a,  &
       imsa_a,     ina_a,      ica_a,     ioc_a,      ibc_a,   &
       ioin_a

      integer, save ::   &
      ipcg1_b_c_a,ipcg2_b_c_a,ipcg3_b_c_a,ipcg4_b_c_a, &
      ipcg5_b_c_a,ipcg6_b_c_a,ipcg7_b_c_a,ipcg8_b_c_a, &
      ipcg9_b_c_a,ipcg1_b_o_a,ipcg2_b_o_a,ipcg3_b_o_a, &
      ipcg4_b_o_a,ipcg5_b_o_a,ipcg6_b_o_a,ipcg7_b_o_a, &
      ipcg8_b_o_a,ipcg9_b_o_a,iopcg1_b_c_a,iopcg2_b_c_a,&
      iopcg3_b_c_a, iopcg4_b_c_a,iopcg5_b_c_a,iopcg6_b_c_a,&
      iopcg7_b_c_a,iopcg8_b_c_a,iopcg1_b_o_a,iopcg2_b_o_a,&
      iopcg3_b_o_a,iopcg4_b_o_a,iopcg5_b_o_a,iopcg6_b_o_a,&
      iopcg7_b_o_a,iopcg8_b_o_a,&
      ipcg1_f_c_a,ipcg2_f_c_a,ipcg3_f_c_a,ipcg4_f_c_a, &
      ipcg5_f_c_a,ipcg6_f_c_a,ipcg7_f_c_a,ipcg8_f_c_a, &
      ipcg9_f_c_a,ipcg1_f_o_a,ipcg2_f_o_a,ipcg3_f_o_a, &
      ipcg4_f_o_a,ipcg5_f_o_a,ipcg6_f_o_a,ipcg7_f_o_a, &
      ipcg8_f_o_a,ipcg9_f_o_a,iopcg1_f_c_a,iopcg2_f_c_a,&
      iopcg3_f_c_a, iopcg4_f_c_a,iopcg5_f_c_a,iopcg6_f_c_a,&
      iopcg7_f_c_a,iopcg8_f_c_a,iopcg1_f_o_a,iopcg2_f_o_a,&
      iopcg3_f_o_a,iopcg4_f_o_a,iopcg5_f_o_a,iopcg6_f_o_a,&
      iopcg7_f_o_a,iopcg8_f_o_a, &
      ismpa_a,ismpbb_a, &
      iglysoa_r1_a, iglysoa_r2_a, iglysoa_oh_a, iglysoa_sfc_a, iglysoa_nh4_a, &
      iant1_c_a,iant2_c_a,iant3_c_a, &
      iant4_c_a,ibiog1_c_a,ibiog2_c_a,ibiog3_c_a,ibiog4_c_a, &
      iant1_o_a,iant2_o_a,iant3_o_a, &
      iant4_o_a,ibiog1_o_a,ibiog2_o_a,ibiog3_o_a,ibiog4_o_a, &
      iasoaX_a, iasoa1_a,iasoa2_a,iasoa3_a,iasoa4_a,&
      ibsoaX_a, ibsoa1_a,ibsoa2_a,ibsoa3_a,ibsoa4_a



      integer, save ::   &
       jnh4so4,    jlvcite,    jnh4hso4,   jnh4no3,    jnh4cl,  &
       jna2so4,    jna3hso4,   jnahso4,    jnano3,     jnacl,   &
       jcaso4,     jcano3,     jcacl2,     jcaco3,     jh2so4,  &
       jhno3,      jhcl,       jhhso4,                          &
       jnh4msa,    jnamsa,     jcamsa2,    jmsa,                &
       joc,        jbc,        join,       jh2o

      integer, save ::   &
      jpcg1_b_c,jpcg2_b_c,jpcg3_b_c,jpcg4_b_c, &
      jpcg5_b_c,jpcg6_b_c,jpcg7_b_c,jpcg8_b_c, &
      jpcg9_b_c,jpcg1_b_o,jpcg2_b_o,jpcg3_b_o, &
      jpcg4_b_o,jpcg5_b_o,jpcg6_b_o,jpcg7_b_o, &
      jpcg8_b_o,jpcg9_b_o,jopcg1_b_c,jopcg2_b_c,&
      jopcg3_b_c, jopcg4_b_c,jopcg5_b_c,jopcg6_b_c,&
      jopcg7_b_c,jopcg8_b_c,jopcg1_b_o,jopcg2_b_o,&
      jopcg3_b_o,jopcg4_b_o,jopcg5_b_o,jopcg6_b_o,&
      jopcg7_b_o,jopcg8_b_o,&
      jpcg1_f_c,jpcg2_f_c,jpcg3_f_c,jpcg4_f_c, &
      jpcg5_f_c,jpcg6_f_c,jpcg7_f_c,jpcg8_f_c, &
      jpcg9_f_c,jpcg1_f_o,jpcg2_f_o,jpcg3_f_o, &
      jpcg4_f_o,jpcg5_f_o,jpcg6_f_o,jpcg7_f_o, &
      jpcg8_f_o,jpcg9_f_o,jopcg1_f_c,jopcg2_f_c,&
      jopcg3_f_c, jopcg4_f_c,jopcg5_f_c,jopcg6_f_c,&
      jopcg7_f_c,jopcg8_f_c,jopcg1_f_o,jopcg2_f_o,&
      jopcg3_f_o,jopcg4_f_o,jopcg5_f_o,jopcg6_f_o,&
      jopcg7_f_o,jopcg8_f_o, &
      jsmpa,jsmpbb, &
      jglysoa_r1, jglysoa_r2, jglysoa_oh, jglysoa_sfc, jglysoa_nh4, &
      jant1_c,jant2_c,jant3_c, &
      jant4_c,jbiog1_c,jbiog2_c,jbiog3_c,jbiog4_c, &
      jant1_o,jant2_o,jant3_o, &
      jant4_o,jbiog1_o,jbiog2_o,jbiog3_o,jbiog4_o, &
      jasoaX,jasoa1,jasoa2,jasoa3,jasoa4,&
      jbsoaX,jbsoa1,jbsoa2,jbsoa3,jbsoa4



      integer, save ::   			&
       jc_h,    jc_nh4, jc_na,  jc_ca,		&
       ja_hso4, ja_so4, ja_no3, ja_cl, ja_msa     




      integer, save ::			&
	iclm_aer,			&  
	jclm_aer,			&  
	kclm_aer,			&  
	kclm_aer_calcbgn,		&  
	kclm_aer_calcend,		&  
	mclm_aer,			&  
	mgas_aer_xfer,			&  
	mdynamic_solver,		&  
	msize_framework,		&  
	jaerosolstate(nbin_a_maxd),	&  
	jphase(nbin_a_maxd),		&  
	jhyst_leg(nbin_a_maxd),		&  
	iprint_input,			&  
	lunerr_aer,			&  
	ncorecnt_aer,       &  
	n2o5_flag				



      integer, save :: istat_mosaic_fe1       
                       
                       
      integer, save :: nfe1_mosaic_cur = 0
                       
      integer, save :: nfe1_mosaic_tot = 0
                       
      integer, save :: iprint_mosaic_fe1 = 1
                       
                       
                       
      integer, save :: iprint_mosaic_perform_stats = 1 
                       
      integer, save :: iprint_mosaic_diag1 = 1 
                       
      integer, save :: iprint_mosaic_input_ok = 1 
                       
                       


      real(kind=8), save ::		&
      	num_a(nbin_a_maxd),		&  
      	dpgn_a(nbin_a_maxd),		&  
      	dp_dry_a(nbin_a_maxd),		&  
      	dp_wet_a(nbin_a_maxd),		&  
      	area_dry_a(nbin_a_maxd),	&  
      	area_wet_a(nbin_a_maxd),	&  
	mass_dry_salt(nbin_a_maxd),	&  
      	mass_dry_a(nbin_a_maxd),	&  
      	mass_wet_a(nbin_a_maxd),	&  
      	mass_soluble_a(nbin_a_maxd),	&  
      	vol_dry_a(nbin_a_maxd),		&  
      	vol_wet_a(nbin_a_maxd),		&  
      	dens_dry_a(nbin_a_maxd),	&  
      	dens_wet_a(nbin_a_maxd),	&  
      	sigmag_a(nbin_a_maxd),		&  
      	water_a(nbin_a_maxd), 		&  
      	water_a_hyst(nbin_a_maxd),	&  
      	water_a_up(nbin_a_maxd),	&  
      	ph(nbin_a_maxd),		&  
        c_as(nbin_a_maxd),          & 
        c_an(nbin_a_maxd),          & 
        a_nh4(nbin_a_maxd),         & 
      	aer(naer,3,nbin_a_maxd),	&  
	aer_sum(3,nbin_a_maxd),		&  
      	aer_percent(naer,3,nbin_a_maxd), &  
      	comp_a(naercomp),		&  
      	electrolyte(nelectrolyte,3,nbin_a_maxd),   &  
      	electrolyte_sum(nelectrolyte,nbin_a_maxd), &  
      	epercent(nelectrolyte,3,nbin_a_maxd),	   &  
      	gas(ngas_volatile+ngas_het),		&  
      	ah2o,				&  
      	ah2o_a(nbin_a_maxd),		&  
      	dpmv(nbin_a_maxd),		&  
      	volume_a(nbin_a_maxd),		&  
	volume_bin(nbin_a_maxd),	&  
      	kelvin(nbin_a_maxd),		&  
	kel(ngas_volatile+ngas_het,nbin_a_maxd),	&  
	kelvin_nh4no3,			&  
	kelvin_nh4cl,			&  
	total_species(ngas_volatile)	   




      integer, save ::			&
	idry_case3a(nbin_a_maxd),	&  
	ieqblm_bin(nbin_a_maxd),	&  
	ieqblm_astem,			&  
        ieqblm_soa,                     &  
	nastem_call,			&  
	nastem_fail,			&  
	isteps_astem,			&  
	nsteps_astem,			&  
        isteps_SOA,                     &
	nsteps_astem_max,		&  
	nmax_ASTEM,			&  
        flagsoap(ngas_soa),             &       
	integrate(ngas_volatile,3,nbin_a_maxd)  


      real(kind=8), save ::			&
	po_soa(ngas_volatile),			&  
	sat_soa(ngas_volatile),			&  
	x_soa(naer),				&  
	sfc_a(ngas_volatile),			&  
	Heff(ngas_volatile,nbin_a_maxd),	&  
	kg(ngas_volatile+ngas_het,nbin_a_maxd),		&  
        fraceq(ngas_volatile,nbin_a_maxd),      &  
	df_gas_s(ngas_volatile,nbin_a_maxd),	&  
	df_gas_l(ngas_volatile,nbin_a_maxd),	&  
        df_gas_o(ngas_volatile,nbin_a_maxd),     &  
	flux_s(ngas_volatile,nbin_a_maxd),	&  
	flux_l(ngas_volatile,nbin_a_maxd),	&  
        flux_o(ngas_volatile,nbin_a_maxd),      &  
	sumkg_h2so4,				&  
	sumkg_msa,				&  
	sumkg_nh3,				&  
	sumkg_hno3,				&  
	sumkg_hcl,				&  
	sumkg_n2o5,				&  
	delta_nh3_max(nbin_a_maxd),		&  
	delta_hno3_max(nbin_a_maxd),		&  
	delta_hcl_max(nbin_a_maxd),		&  
	keq_nh4no3,				&  
	keq_nh4cl,				&  
	Keq_nh4no3_0,				&  
	Keq_nh4cl_0,				&  
	volatile_s(ngas_volatile,nbin_a_maxd), 	&  
	phi_volatile_s(ngas_volatile,nbin_a_maxd),	&  
	phi_volatile_l(ngas_volatile,nbin_a_maxd),	&  
        phi_volatile_o(ngas_volatile,nbin_a_maxd),      &  
	phi_nh4no3_s,				&  
	phi_nh4cl_s,				&  
	sum_vdf_s(ngas_volatile),		&  
	sum_vol_s(ngas_volatile),		&  
	sum_bin_s(ngas_volatile),		&  
	avg_df_gas_s(ngas_volatile),		&  
	h_s_i_m(ngas_volatile,nbin_a_maxd),	&  
	alpha_astem,				&  
	rtol_eqb_astem,				&  
	ptol_mol_astem,				&  
	nsteps_astem_avg			   

      integer, parameter :: glysoa_param_off     = 0, &
                            glysoa_param_simple  = 1, &
                            glysoa_param_complex = 2
      integer, save      :: glysoa_param



      integer, save ::   		&
      	jsalt_index(nsalt),		&
      	jsulf_poor(211),		&
      	jsulf_rich(71),			&
      	jsalt_present(nsalt),		&
      	nmax_mesa,			&
      	nmesa_call,			&
      	nmesa_fail,			&
	iter_mesa(nbin_a_maxd),		&
	niter_mesa,			&
      	niter_mesa_max


      real(kind=8), save ::   		&
	eleliquid(nelectrolyte),	&
	flux_sl(nsalt),			&
	phi_salt(nsalt),		&
	phi_salt_old(nsalt),		&
	phi_bar(nsalt),			&
	alpha_salt(nsalt),		&
	sat_ratio(nsalt),		&
	hsalt(nsalt),			&
	hsalt_max,			&
	frac_salt_liq(nsalt),		&
	frac_salt_solid(nsalt),		&
	growth_factor(nbin_a_maxd),	&
	d_mdrh(63,4),			&  
	mdrh(nbin_a_maxd),		&
	mdrh_t(63),			&
	molality0(nelectrolyte),	&
	rtol_mesa,			&
	niter_mesa_avg




      character(len=8), save ::		&
	ename(nelectrolyte),		&  
	aer_name(naer),			&  
	gas_name(ngas_volatile+ngas_het)		   

      character(len=6), save ::		&
	phasestate(4)


      real(kind=8), save ::			&
	t_k,				&  
	p_atm,				&  
	rh_pc,				&  
	cair_mol_cc,			&  
	cair_mol_m3,			&  
	conv1a,				&
	conv1b,				&
	conv2a,				&
	conv2b,				&
	mw_electrolyte(nelectrolyte),	&  
	mw_aer_mac(naer),		&  
	mw_comp_a(naercomp),		&  
	mw_c(ncation),			&  
	mw_a(nanion),			&  
	dens_electrolyte(nelectrolyte),	&  
	dens_aer_mac(naer),		&  
	dens_comp_a(naercomp),		&  
	kappa_nonelectro(naer),		&  
	partial_molar_vol(ngas_volatile+ngas_het), & 
	sigma_water,			&  
	sigma_soln(nbin_a_maxd),	&  
	keq_gl(nrxn_aer_gl),		&  
	keq_ll(nrxn_aer_ll),		&  
	keq_sg(nrxn_aer_sg),		&  
	keq_sl(nrxn_aer_sl), 		&  
	kp_nh3, 			&  
	kp_nh4no3, 			&  
	kp_nh4no3_0, 			&  
	kp_nh4cl,	                &  
	kp_nh4cl_0,			&   
	frac_n2o5_h2o(nbin_a_maxd)	

      complex, save ::			&
      		ref_index_a(naercomp),	&  
      		ri_avg_a(nbin_a_maxd)	   





      real(kind=8), save ::			&
	mc(ncation,nbin_a_maxd),		&  
	ma(nanion,nbin_a_maxd),			&  
	msulf,					&  
	zc(ncation),				&  
	za(nanion),				&  
	gam(nelectrolyte,nbin_a_maxd),		&
	gam_ratio(nbin_a_maxd),			&
	log_gamz(nelectrolyte,nelectrolyte),	&
	log_gam(nelectrolyte),			&
	activity(nelectrolyte,nbin_a_maxd),	&
	xeq_a(nanion),				&
	xeq_c(ncation),				&
	na_ma(nanion),				&
	nc_mc(ncation),				&
	a_zsr(6,nelectrolyte),			&  
	b_zsr(nelectrolyte),			&  
	aw_min(nelectrolyte),			&  
	b_mtem(6,nelectrolyte,nelectrolyte)	   




      real(kind=8), save ::	&
	tot_so4_in,	&
	tot_no3_in,	&
	tot_cl_in,	&
	tot_nh4_in,	&
	tot_na_in,	&
	tot_ca_in,	&
	tot_so4_out,	&
	tot_no3_out,	&
	tot_cl_out,	&
	tot_nh4_out,	&
	tot_na_out,	&
	tot_ca_out,	&
	diff_so4,	&
	diff_no3,	&
	diff_cl,	&
	diff_nh4,	&
	diff_na,	&
	diff_ca,	&
	reldiff_so4,	&
	reldiff_no3,	&
	reldiff_cl,	&
	reldiff_nh4,	&
	reldiff_na,	&
	reldiff_ca




      end module module_data_mosaic_therm
