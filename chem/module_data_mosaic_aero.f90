module module_data_mosaic_aero








  
  use module_data_mosaic_kind, only:  r8
  
  implicit none
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  

  integer, parameter :: nbin_a_max = 8  

  integer, save :: nbin_a     = -999888777  
  
  
  integer, parameter :: ngas_ioa = 4+1	




  integer, parameter :: ngas_soa = 68+2+16+10	

  integer, parameter :: ngas_het = 2+2
  integer, parameter :: ngas_volatile = ngas_ioa + ngas_soa + 2 
  integer, parameter :: ngas_aerchtot = ngas_volatile + ngas_het 

  integer, parameter :: naer = 11+68+2+16+5+10	

  integer, parameter :: naercomp = 26+68+2+16+5+10	

  integer, parameter :: nelectrolyte = 22+1 
  integer, parameter :: nsalt   = 12+3	
  integer, parameter :: nsoluble= 16+4  
  integer, parameter :: ncation = 4	
  integer, parameter :: nanion  = 4+1	
  
  integer, parameter :: nrxn_aer_gl = 4 
  integer, parameter :: nrxn_aer_ll = 3 
  integer, parameter :: nrxn_aer_sg = 2 
  integer, parameter :: nrxn_aer_sl = nsalt
  
  integer, parameter :: mASTEM = 1	
  integer, parameter :: mLSODE = 2	
  integer, parameter :: mMODAL  = 1	
  integer, parameter :: mUNSTRUCTURED = 2 
  
  integer, parameter :: mSECTIONAL = 3	
  integer, parameter :: mON     = 1	
  integer, parameter :: mOFF    = 0    
  integer, parameter :: mYES	= mON	
  integer, parameter :: mNO	= mOFF	
  
  integer, parameter :: jsolid = 1
  integer, parameter :: jliquid= 2
  integer, parameter :: jtotal = 3
  
  integer, parameter :: jhyst_lo = 0	
  integer, parameter :: jhyst_up = 1 	
  integer, parameter :: jhyst_undefined = -1	
  
  
  integer, parameter :: mhyst_uporlo_jhyst = 1	
  
  
  
  integer, parameter :: mhyst_uporlo_waterhyst = 2	
  
  integer, parameter :: mhyst_force_up = 3	
  integer, parameter :: mhyst_force_lo = 4	
  
  real(r8), parameter :: xhyst_up_crustal_thresh = 0.30_r8
  
  
  
  
  

  integer, parameter :: no_aerosol = 0	
  integer, parameter :: all_solid  = 1 
  integer, parameter :: all_liquid = 2 
  integer, parameter :: mixed      = 3	
  
  integer, parameter :: soluble   = 1  
  integer, parameter :: insoluble = 2  

  integer, parameter :: MDRH_T_NUM     = 63     
  integer, parameter :: jsulf_poor_NUM = 211    
  integer, parameter :: jsulf_rich_NUM = 71     
  integer, parameter :: d_mdrh_DIM2    = 4      
  
  real(r8), parameter :: mass_cutoff = 1.e-6	
  
  real(r8), parameter :: density_min_allow = 1.0	
  real(r8), parameter :: density_max_allow = 3.0	
  real(r8), parameter :: ah2o_max = 0.99                

  
  
  type :: mosaic_vars_aa_type
     integer :: it_host
     integer :: it_mosaic
     integer, dimension(6) :: hostgridinfo(6)
     integer :: idiagbb_host
     integer :: f_mos_fail
     integer :: isteps_astem
     integer :: isteps_astem_max
     integer :: jastem_call
     integer :: jastem_fail
     integer :: jmesa_call
     integer :: jmesa_fail
     integer :: niter_mesa_max
     integer :: nmax_astem
     integer :: nmax_mesa
     integer :: fix_astem_negative
     logical :: flag_itr_kel
     logical :: zero_water_flag
     real(r8) :: cumul_steps_astem
     real(r8) :: niter_mesa
     real(r8) :: swdown
     real(r8), dimension(5,4) :: xnerr_astem_negative
     integer, dimension(:), allocatable :: iter_mesa
  end type mosaic_vars_aa_type

  
  
  
  
  
  
  integer, save :: isoa_first
  
  integer, save :: jsoa_first

  
  integer, save ::   &
       ih2so4_g,     ihno3_g,      ihcl_g,      inh3_g,   &
       imsa_g,       in2o5_g,      iclno2_g,   &
       iaro1_g,      iaro2_g,      ialk1_g,     iole1_g,   &   
       iapi1_g,      iapi2_g,      ilim1_g,     ilim2_g        
  
  
  integer, save ::   &
       iso4_a,     ino3_a,     icl_a,     inh4_a,     ico3_a,   &
       imsa_a,     ina_a,      ica_a,     ioc_a,      ibc_a,   &
       ioin_a,   &
       iaro1_a,    iaro2_a,   ialk1_a,    iole1_a,   &   
       iapi1_a,    iapi2_a,    ilim1_a,   ilim2_a        
  
  
  integer, save ::   &
       jnh4so4,    jlvcite,    jnh4hso4,   jnh4no3,    jnh4cl,   &
       jna2so4,    jna3hso4,   jnahso4,    jnano3,     jnacl,   &
       jcaso4,     jcano3,     jcacl2,     jcaco3,     jh2so4,   &
       jhno3,      jhcl,       jhhso4,   &
       jnh4msa,    jnamsa,     jcamsa2,    jmsa,   &
       joc,        jbc,        join,       jh2o,   &
       jaro1,      jaro2,      jalk1,      jole1,   &   
       japi1,      japi2,      jlim1,      jlim2        
  
  
  integer, save ::   &
       jc_h,    jc_nh4, jc_na,  jc_ca,   &
       ja_hso4, ja_so4, ja_no3, ja_cl, ja_msa     
  
  
      integer, save ::                                        &
      ipcg1_b_c_g, ipcg2_b_c_g, ipcg3_b_c_g, ipcg4_b_c_g, &
      ipcg5_b_c_g, ipcg6_b_c_g, ipcg7_b_c_g, ipcg8_b_c_g, ipcg9_b_c_g, &
      ipcg1_b_o_g, ipcg2_b_o_g, ipcg3_b_o_g, ipcg4_b_o_g, &
      ipcg5_b_o_g, ipcg6_b_o_g, ipcg7_b_o_g, ipcg8_b_o_g, ipcg9_b_o_g, &
      iopcg1_b_c_g, iopcg2_b_c_g, iopcg3_b_c_g, iopcg4_b_c_g, &
      iopcg5_b_c_g, iopcg6_b_c_g, iopcg7_b_c_g, iopcg8_b_c_g, &
      iopcg1_b_o_g, iopcg2_b_o_g, iopcg3_b_o_g, iopcg4_b_o_g, &
      iopcg5_b_o_g, iopcg6_b_o_g, iopcg7_b_o_g, iopcg8_b_o_g, &
      ipcg1_f_c_g, ipcg2_f_c_g, ipcg3_f_c_g, ipcg4_f_c_g, &
      ipcg5_f_c_g, ipcg6_f_c_g, ipcg7_f_c_g, ipcg8_f_c_g, ipcg9_f_c_g, &
      ipcg1_f_o_g, ipcg2_f_o_g, ipcg3_f_o_g, ipcg4_f_o_g, &
      ipcg5_f_o_g, ipcg6_f_o_g, ipcg7_f_o_g, ipcg8_f_o_g, ipcg9_f_o_g, &
      iopcg1_f_c_g, iopcg2_f_c_g, iopcg3_f_c_g, iopcg4_f_c_g, &
      iopcg5_f_c_g, iopcg6_f_c_g, iopcg7_f_c_g, iopcg8_f_c_g, &
      iopcg1_f_o_g, iopcg2_f_o_g, iopcg3_f_o_g, iopcg4_f_o_g, &
      iopcg5_f_o_g, iopcg6_f_o_g, iopcg7_f_o_g, iopcg8_f_o_g, &
      iant1_c_g, iant2_c_g, iant3_c_g, iant4_c_g, &
      iant1_o_g, iant2_o_g, iant3_o_g, iant4_o_g, &
      ibiog1_c_g, ibiog2_c_g, ibiog3_c_g, ibiog4_c_g, &
      ibiog1_o_g, ibiog2_o_g, ibiog3_o_g, ibiog4_o_g, &
      ismpa_g, ismpbb_g, &
      iasoaX_g, iasoa1_g, iasoa2_g, iasoa3_g, iasoa4_g, &
      ibsoaX_g, ibsoa1_g, ibsoa2_g, ibsoa3_g, ibsoa4_g, &
      igly, iho

  
      integer, save ::   &
      ipcg1_b_c_a, ipcg2_b_c_a, ipcg3_b_c_a, ipcg4_b_c_a, &
      ipcg5_b_c_a, ipcg6_b_c_a, ipcg7_b_c_a, ipcg8_b_c_a, ipcg9_b_c_a, &
      ipcg1_b_o_a, ipcg2_b_o_a, ipcg3_b_o_a, ipcg4_b_o_a, &
      ipcg5_b_o_a, ipcg6_b_o_a, ipcg7_b_o_a, ipcg8_b_o_a, ipcg9_b_o_a, &
      iopcg1_b_c_a, iopcg2_b_c_a, iopcg3_b_c_a, iopcg4_b_c_a, &
      iopcg5_b_c_a, iopcg6_b_c_a, iopcg7_b_c_a, iopcg8_b_c_a, &
      iopcg1_b_o_a, iopcg2_b_o_a, iopcg3_b_o_a, iopcg4_b_o_a, &
      iopcg5_b_o_a, iopcg6_b_o_a, iopcg7_b_o_a, iopcg8_b_o_a,&
      ipcg1_f_c_a, ipcg2_f_c_a, ipcg3_f_c_a, ipcg4_f_c_a, &
      ipcg5_f_c_a, ipcg6_f_c_a, ipcg7_f_c_a, ipcg8_f_c_a, ipcg9_f_c_a, &
      ipcg1_f_o_a, ipcg2_f_o_a, ipcg3_f_o_a, ipcg4_f_o_a, &
      ipcg5_f_o_a, ipcg6_f_o_a, ipcg7_f_o_a, ipcg8_f_o_a, ipcg9_f_o_a, &
      iopcg1_f_c_a, iopcg2_f_c_a, iopcg3_f_c_a, iopcg4_f_c_a, &
      iopcg5_f_c_a, iopcg6_f_c_a, iopcg7_f_c_a, iopcg8_f_c_a, &
      iopcg1_f_o_a, iopcg2_f_o_a, iopcg3_f_o_a, iopcg4_f_o_a, &
      iopcg5_f_o_a, iopcg6_f_o_a, iopcg7_f_o_a, iopcg8_f_o_a, &
      ismpa_a, ismpbb_a, &
      iglysoa_r1_a, iglysoa_r2_a, iglysoa_oh_a, iglysoa_sfc_a, iglysoa_nh4_a, &
      iant1_c_a, iant2_c_a, iant3_c_a, iant4_c_a,&
      iant1_o_a, iant2_o_a, iant3_o_a, iant4_o_a,&
      ibiog1_c_a, ibiog2_c_a, ibiog3_c_a, ibiog4_c_a, &
      ibiog1_o_a, ibiog2_o_a, ibiog3_o_a, ibiog4_o_a, &
      iasoaX_a, iasoa1_a, iasoa2_a, iasoa3_a, iasoa4_a, &
      ibsoaX_a, ibsoa1_a, ibsoa2_a, ibsoa3_a, ibsoa4_a

  
      integer, save ::   &
      jpcg1_b_c, jpcg2_b_c, jpcg3_b_c, jpcg4_b_c, &
      jpcg5_b_c, jpcg6_b_c, jpcg7_b_c, jpcg8_b_c, jpcg9_b_c, &
      jpcg1_b_o, jpcg2_b_o, jpcg3_b_o, jpcg4_b_o, &
      jpcg5_b_o, jpcg6_b_o, jpcg7_b_o, jpcg8_b_o, jpcg9_b_o, &
      jopcg1_b_c, jopcg2_b_c, jopcg3_b_c, jopcg4_b_c, &
      jopcg5_b_c, jopcg6_b_c, jopcg7_b_c, jopcg8_b_c, &
      jopcg1_b_o, jopcg2_b_o, jopcg3_b_o, jopcg4_b_o, &
      jopcg5_b_o, jopcg6_b_o, jopcg7_b_o, jopcg8_b_o, &
      jpcg1_f_c, jpcg2_f_c, jpcg3_f_c, jpcg4_f_c, &
      jpcg5_f_c, jpcg6_f_c, jpcg7_f_c, jpcg8_f_c, jpcg9_f_c, &
      jpcg1_f_o, jpcg2_f_o, jpcg3_f_o, jpcg4_f_o, &
      jpcg5_f_o, jpcg6_f_o, jpcg7_f_o, jpcg8_f_o, jpcg9_f_o, &
      jopcg1_f_c, jopcg2_f_c, jopcg3_f_c, jopcg4_f_c, &
      jopcg5_f_c, jopcg6_f_c, jopcg7_f_c, jopcg8_f_c, &
      jopcg1_f_o, jopcg2_f_o, jopcg3_f_o, jopcg4_f_o, &
      jopcg5_f_o, jopcg6_f_o, jopcg7_f_o, jopcg8_f_o, &
      jsmpa, jsmpbb, &
      jglysoa_r1, jglysoa_r2, jglysoa_oh, jglysoa_sfc, jglysoa_nh4, &
      jant1_c, jant2_c, jant3_c, jant4_c,  &
      jant1_o, jant2_o, jant3_o, jant4_o,  &
      jbiog1_c, jbiog2_c, jbiog3_c, jbiog4_c, &
      jbiog1_o, jbiog2_o, jbiog3_o, jbiog4_o, &
      jasoaX, jasoa1, jasoa2, jasoa3, jasoa4, &
      jbsoaX, jbsoa1, jbsoa2, jbsoa3, jbsoa4

  
  

  
  
  
  
  
  
  integer, save :: use_cam5mam_soa_params  = 0   
  integer, save :: use_cam5mam_accom_coefs = 0   

  integer, save :: 	&
       
       
       
       
       
       mclm_aer,			   &  
       mGAS_AER_XFER,			   &  
       mDYNAMIC_SOLVER,		   &  
       mSIZE_FRAMEWORK,		   &  
       mhyst_method,			   &  
       maersize_init_flag1,		   &  
       mcoag_flag1,			   &  
       mmovesect_flag1,		   &  
       mnewnuc_flag1,			   &  
       msectional_flag1,		   &  
       msectional_flag2,		   &  
       method_bcfrac,      		   &  
       method_kappa,       		   &  
       ifreq_coag,            		   &  
       ipmcmos_aero,                       &
       maeroptic_aero,      		   &
       msoa_flag1,      		   &  
       msoa_vbs_info(9)    		      

       
       
  
  
  
  
  
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
  
  real(r8), save :: 	&
       dlo_aersize_init, 		   &  
       dhi_aersize_init, 		   &  
       xcutlo_atype_md1_init, 		   &  
       xcuthi_atype_md1_init, 		   &  
       xcutlo_atype_md2_init, 		   &  
       xcuthi_atype_md2_init  		      
  
  integer, save :: 	&
       method_atype_md1_init, 		   &  
       method_atype_md2_init  		      
  
  
  
  
  integer, save ::	   &
  
  
  
  
  
  
  
  
  
  nmax_ASTEM
  
  
  integer, save :: m_gas2bin_uptk_flag = 1


  integer, save :: i_gas2bin_uptk_flag(ngas_aerchtot,nbin_a_max)




  real(r8), save ::	&
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  alpha_ASTEM,			     &  
  rtol_eqb_ASTEM,			     &  
  ptol_mol_ASTEM
  
  
  
  
  
  integer, save ::	&
       jsalt_index(nsalt),   &
       jsulf_poor(jsulf_poor_NUM),   &
       jsulf_rich(jsulf_rich_NUM),   &
       
       Nmax_mesa
       
       
       
       

  real(r8), save ::	&
       
       
       
       
       
       
       
       
       
       
       
       
       d_mdrh(MDRH_T_NUM,d_mdrh_DIM2),		   &  
       
       
       
       rtol_mesa
       
       
       
       
  
  
  
  character(len= 6), save :: phasestate(0:4)
  character(len=16), save :: ename(nelectrolyte)	
  character(len=16), save :: aer_name(naer)		
  character(len=16), save :: gas_name(ngas_aerchtot)	
  
  real(r8), save ::      &
       
       
       
       
       
       
       
       
       
       mw_electrolyte(nelectrolyte),		   &  
       mw_aer_mac(naer),			   &  
       mw_comp_a(naercomp),			   &  
       mw_c(ncation),				   &  
       mw_a(nanion),				   &  
       mw_gas(ngas_aerchtot),			   &  
       dens_electrolyte(nelectrolyte),		   &  
       dens_aer_mac(naer),			   &  
       dens_comp_a(naercomp),			   &  
       kappa_aer_mac(naer),			   &  
       partial_molar_vol(ngas_aerchtot),	   &  
       v_molar_gas(ngas_aerchtot)		      
       
       
       
       
       
       
       
       
       
  
  complex   &
       ref_index_a(naercomp)
  
  
  
  
  
  
  
  real(r8), save ::	&
       
       
       
       zc(Ncation),			     &  
       za(Nanion),			     &  
       
       
       
       
       
       
       
       
       
       a_zsr(6,nelectrolyte),		     &  
       b_zsr(nelectrolyte),		     &  
       aw_min(nelectrolyte),		     &  
       b_mtem(6,nelectrolyte,nelectrolyte)     
  
  
  

























  
  

end module module_data_mosaic_aero
