module module_sorgam_vbs_aqchem
  
  
  
  
  
  
  REAL, PARAMETER :: epsilc = 1.0E-16
  
  REAL, PARAMETER :: qcldwtr_cutoff = 1.0e-6 
  
  REAL, PARAMETER :: mwdry = 28.966  
  REAL, PARAMETER :: mwso4 = 96.00   
  REAL, PARAMETER :: mwno3 = 62.0    
  REAL, PARAMETER :: mwnh4 = 18.0985 
  REAL, PARAMETER :: mwna  = 22.990  
  REAL, PARAMETER :: mwcl  = 35.453  
  
  
  
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
	
  

  INTEGER, PARAMETER :: LACL        =  1  
  INTEGER, PARAMETER :: LNH4L       =  2  
  INTEGER, PARAMETER :: LCAL        =  3  
  INTEGER, PARAMETER :: LNAACCL     =  4  
  INTEGER, PARAMETER :: LOHL        =  5  
  INTEGER, PARAMETER :: LSO4ACCL    =  6  
  INTEGER, PARAMETER :: LHSO4ACCL   =  7  
  INTEGER, PARAMETER :: LSO3L       =  8  
  INTEGER, PARAMETER :: LHSO3L      =  9  
  INTEGER, PARAMETER :: LSO2L       = 10  
  INTEGER, PARAMETER :: LCO3L       = 11  
  INTEGER, PARAMETER :: LHCO3L      = 12  
  INTEGER, PARAMETER :: LCO2L       = 13  
  INTEGER, PARAMETER :: LNO3ACCL    = 14  
  INTEGER, PARAMETER :: LNH3L       = 15  
  INTEGER, PARAMETER :: LCLACCL     = 16  
  INTEGER, PARAMETER :: LH2O2L      = 17  
  INTEGER, PARAMETER :: LO3L        = 18  
  INTEGER, PARAMETER :: LFEL        = 19  
  INTEGER, PARAMETER :: LMNL        = 20  
  INTEGER, PARAMETER :: LAL         = 21  
  INTEGER, PARAMETER :: LFOAL       = 22  
  INTEGER, PARAMETER :: LHCO2L      = 23  
  INTEGER, PARAMETER :: LMHPL       = 24  
  INTEGER, PARAMETER :: LPAAL       = 25  
  INTEGER, PARAMETER :: LHCLL       = 26  
  INTEGER, PARAMETER :: LPRIML      = 27  
  INTEGER, PARAMETER :: LMGL        = 28  
  INTEGER, PARAMETER :: LKL         = 29  
  INTEGER, PARAMETER :: LBL         = 30  
  INTEGER, PARAMETER :: LHNO3L      = 31  
  INTEGER, PARAMETER :: LPRIMCORL   = 32  
  INTEGER, PARAMETER :: LNUMCORL    = 33  
  INTEGER, PARAMETER :: LTS6CORL    = 34  
  INTEGER, PARAMETER :: LNACORL     = 35  
  INTEGER, PARAMETER :: LCLCORL     = 36  
  INTEGER, PARAMETER :: LNO3CORL    = 37  
  INTEGER, PARAMETER :: LORGAL      = 38  
  INTEGER, PARAMETER :: LORGPL      = 39  
  INTEGER, PARAMETER :: LORGBL      = 40  
  INTEGER, PARAMETER :: LECL        = 41  
  
  contains
  

  
	subroutine sorgam_vbs_aqchem_driver(   &
	    id, ktau, ktauc, dtstepc, config_flags,   &
	    p_phy, t_phy, rho_phy, alt, dz8w,  &
	    moist, chem,   &
	    gas_aqfrac, numgas_aqfrac,   &
	    ids,ide, jds,jde, kds,kde,   &
	    ims,ime, jms,jme, kms,kme,   &
	    its,ite, jts,jte, kts,kte )
    
    use module_ctrans_aqchem, only: aqchem
  	
  	use module_configure, only: grid_config_rec_type
    
    use module_state_description, only: &
      num_chem, &
  		num_moist, &

      p_so2, &
      p_sulf, &
      p_nh3, &
      p_h2o2, &
      p_o3, &
      p_op1, &
      p_ora1, &
      p_paa, &
      p_hno3, &
      p_n2o5, &
      p_so4cwi, &
      p_nh4cwi, &
      p_no3cwi, &
      p_so4cwj, &
      p_nh4cwj, &
      p_no3cwj, &
      p_nacwi, &
      p_nacwj, &
      p_clcwi, &
      p_clcwj, &


      p_qv, &
      p_qc, &

      p_facd, &
      p_mepx, &
      p_pacd, &
      CB05_SORG_VBS_AQ_KPP
    
    use module_data_sorgam_vbs, only: cw_phase, nphase_aer
    
    implicit none
    
    
    
    
    
  	
    
    
    
    


    
    

    
    
    
  	integer, intent(in) ::   &
  		id, ktau, ktauc,   &
  		numgas_aqfrac,   &
  		ids, ide, jds, jde, kds, kde,   &
  		ims, ime, jms, jme, kms, kme,   &
  		its, ite, jts, jte, kts, kte
    
    
    type(grid_config_rec_type), intent(in) :: config_flags
    
    
    real, intent(in) :: dtstepc
    
    
    
    
    
    
    
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: &
      p_phy, t_phy, rho_phy, alt, dz8w
    
    
    
    
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: moist
    
    
    
    
    real, intent(inout), dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: chem
    
    
    
    real, intent(inout), dimension( ims:ime, kms:kme, jms:jme, numgas_aqfrac ) :: gas_aqfrac
    
    
    
    
    
    real, dimension (ngas)  :: gas     
    real, dimension (naer)  :: aerosol 
    real, dimension (nliqs) :: liquid  
    
    real, dimension (ngas) :: gaswdep 
    real, dimension (naer) :: aerwdep 
    real                   :: hpwdep  
    
    real :: precip    
    real :: airm      
    real :: rho_dry   
    real :: h2o_aq    
    real :: h2o_total 
    
    real :: alfa0 
    real :: alfa2 
    real :: alfa3 

    
    
    
    
    integer :: it, jt, kt
    
    real :: conv_factor
    
    

    if ((cw_phase .le. 0) .or. (cw_phase .gt. nphase_aer)) then
      write(*,*) '*** module_sorgam_aqchem - cw_phase not active'
      return
    endif
    
  	write(*,'(a,8(1x,i6))') 'entering module_sorgam_aqchem - ktau =', ktau
    
    
    
    
    precip = 0.0 
    
    alfa0 = 0.0
    alfa2 = 0.0
    alfa3 = 0.0
    
    
    
    gaswdep(:) = 0.0
    aerwdep(:) = 0.0
    hpwdep  = 0.0
    
    
    
    do jt = jts, jte
  	do it = its, ite
  	do kt = kts, kte
    
      if (moist(it,kt,jt,p_qc).gt.qcldwtr_cutoff) then
        
        
        airm = 1000.0*rho_phy(it,kt,jt)*dz8w(it,kt,jt)/mwdry 
        
        
        rho_dry = 1.0/alt(it,kt,jt) 
        
        
        h2o_aq = moist(it,kt,jt,p_qc)*rho_dry 
        
        
        h2o_total = (moist(it,kt,jt,p_qc)+moist(it,kt,jt,p_qv))*rho_dry 
        
        
        
        
        gas(:) = 0.0
        



          gas(lco2) = 380.0e-6

        
        if (p_so2 .gt. 1)  gas(lso2)   = chem(it,kt,jt,p_so2)*1.0e-6
        if (p_hno3 .gt. 1) gas(lhno3)  = chem(it,kt,jt,p_hno3)*1.0e-6
        if (p_n2o5 .gt. 1) gas(ln2o5)  = chem(it,kt,jt,p_n2o5)*1.0e-6
        if (p_nh3 .gt. 1)  gas(lnh3)   = chem(it,kt,jt,p_nh3)*1.0e-6
        if (p_h2o2 .gt. 1) gas(lh2o2)  = chem(it,kt,jt,p_h2o2)*1.0e-6
        if (p_o3 .gt. 1)   gas(lo3)    = chem(it,kt,jt,p_o3)*1.0e-6



        if (p_sulf .gt. 1) gas(lh2so4) = chem(it,kt,jt,p_sulf)*1.0e-6


        if (config_flags%chem_opt==CB05_SORG_VBS_AQ_KPP) then
        if (p_facd .gt. 1) gas(lfoa)   = chem(it,kt,jt,p_facd)*1.0e-6
        if (p_mepx .gt. 1)  gas(lmhp)   = chem(it,kt,jt,p_mepx)*1.0e-6
        if (p_pacd .gt. 1)  gas(lpaa)   = chem(it,kt,jt,p_pacd)*1.0e-6

        else
        if (p_ora1 .gt. 1) gas(lfoa)   = chem(it,kt,jt,p_ora1)*1.0e-6
        if (p_op1 .gt. 1)  gas(lmhp)   = chem(it,kt,jt,p_op1)*1.0e-6
        if (p_paa .gt. 1)  gas(lpaa)   = chem(it,kt,jt,p_paa)*1.0e-6        
        end if
        
        
        
        
        
        
        
        aerosol(:) = 0.0
        
        aerosol(lso4akn)  = chem(it,kt,jt,p_so4cwi)*1.0e-9*mwdry/mwso4 
        aerosol(lnh4akn)  = chem(it,kt,jt,p_nh4cwi)*1.0e-9*mwdry/mwnh4 
        aerosol(lno3akn)  = chem(it,kt,jt,p_no3cwi)*1.0e-9*mwdry/mwno3 


        
        aerosol(lorgaakn) = 0.0                                        
        aerosol(lorgpakn) = 0.0                                        
        aerosol(lorgbakn) = 0.0                                        
        aerosol(lecakn)   = 0.0                                        
        aerosol(lpriakn)  = 0.0                                        
        
        aerosol(lso4acc)  = chem(it,kt,jt,p_so4cwj)*1.0e-9*mwdry/mwso4 
        aerosol(lnh4acc)  = chem(it,kt,jt,p_nh4cwj)*1.0e-9*mwdry/mwnh4 
        aerosol(lno3acc)  = chem(it,kt,jt,p_no3cwj)*1.0e-9*mwdry/mwno3 
        aerosol(lnaacc)   = chem(it,kt,jt,p_nacwj)*1.0e-9*mwdry/mwna   
        aerosol(lclacc)   = chem(it,kt,jt,p_clcwj)*1.0e-9*mwdry/mwcl   
        
        aerosol(lorgaacc) = 0.0                                        
        aerosol(lorgpacc) = 0.0                                        
        aerosol(lorgbacc) = 0.0                                        
        aerosol(lecacc)   = 0.0                                        
        aerosol(lpriacc)  = 0.0                                        
        


        aerosol(lnacor)   = 0.0                                        
        aerosol(lclcor)   = 0.0                                        
        aerosol(lpricor)  = 0.0                                        


        aerosol(LA3FE) = 0.01*alt(it,kt,jt)*1.0e-9*mwdry/55.8
        aerosol(LB2MN) = 0.005*alt(it,kt,jt)*1.0e-9*mwdry/54.9       
 
        
        
        liquid(:) = 0.0
        
        call aqchem( &
          t_phy(it,kt,jt), &
          p_phy(it,kt,jt), &
          dtstepc, &
          precip, &
          h2o_aq, &
          h2o_total, &
          airm, &
          alfa0, &
          alfa2, &
          alfa3, &
          gas, &
          aerosol, &
          liquid, &
          gaswdep, &
          aerwdep, &
          hpwdep)
        
        
        
        

        if (p_so2 .gt. 1)  chem(it,kt,jt,p_so2)  = gas(lso2)*1.0e6
        if (p_hno3 .gt. 1) chem(it,kt,jt,p_hno3) = gas(lhno3)*1.0e6
        if (p_n2o5 .gt. 1) chem(it,kt,jt,p_n2o5) = gas(ln2o5)*1.0e6
        if (p_nh3 .gt. 1)  chem(it,kt,jt,p_nh3)  = gas(lnh3)*1.0e6
        if (p_h2o2 .gt. 1) chem(it,kt,jt,p_h2o2) = gas(lh2o2)*1.0e6
        if (p_o3 .gt. 1)   chem(it,kt,jt,p_o3)   = gas(lo3)*1.0e6



        if (p_sulf .gt. 1) chem(it,kt,jt,p_sulf) = gas(lh2so4)*1.0e6


        if (config_flags%chem_opt==CB05_SORG_VBS_AQ_KPP) then
        if (p_facd .gt. 1) chem(it,kt,jt,p_facd) = gas(lfoa)*1.0e6
        if (p_mepx .gt. 1)  chem(it,kt,jt,p_mepx)  = gas(lmhp)*1.0e6
        if (p_pacd .gt. 1)  chem(it,kt,jt,p_pacd)  = gas(lpaa)*1.0e6

        else
        if (p_ora1 .gt. 1) chem(it,kt,jt,p_ora1) = gas(lfoa)*1.0e6
        if (p_op1 .gt. 1)  chem(it,kt,jt,p_op1)  = gas(lmhp)*1.0e6
        if (p_paa .gt. 1)  chem(it,kt,jt,p_paa)  = gas(lpaa)*1.0e6
        end if
        
        
        
        
        chem(it,kt,jt,p_so4cwi) = aerosol(lso4akn) *1.0e9/mwdry*mwso4 
        chem(it,kt,jt,p_nh4cwi) = aerosol(lnh4akn) *1.0e9/mwdry*mwnh4 
        chem(it,kt,jt,p_no3cwi) = aerosol(lno3akn) *1.0e9/mwdry*mwno3 
        chem(it,kt,jt,p_nacwi)  = aerosol(lnaakn)  *1.0e9/mwdry*mwna  
        chem(it,kt,jt,p_clcwi)  = aerosol(lclakn)  *1.0e9/mwdry*mwcl  
        





        
        chem(it,kt,jt,p_so4cwj) = aerosol(lso4acc) *1.0e9/mwdry*mwso4 
        chem(it,kt,jt,p_nh4cwj) = aerosol(lnh4acc) *1.0e9/mwdry*mwnh4 
        chem(it,kt,jt,p_no3cwj) = aerosol(lno3acc) *1.0e9/mwdry*mwno3 
        chem(it,kt,jt,p_nacwj)  = aerosol(lnaacc)  *1.0e9/mwdry*mwna  
        chem(it,kt,jt,p_clcwj)  = aerosol(lclacc)  *1.0e9/mwdry*mwcl  
        





                                  





        
        
        
        gas_aqfrac(it,kt,jt,:) = 0.0
        
        conv_factor = 1.0E-3*moist(it,kt,jt,p_qc)*mwdry 
        

        if (p_so2  .gt. 1 .and. gas(lso2)  .gt. epsilc) gas_aqfrac(it,kt,jt,p_so2)  = conv_factor*liquid(lso2l)/gas(lso2)
        if (p_nh3  .gt. 1 .and. gas(lnh3)  .gt. epsilc) gas_aqfrac(it,kt,jt,p_nh3)  = conv_factor*liquid(lnh3l)/gas(lnh3)
        if (p_hno3 .gt. 1 .and. gas(lhno3) .gt. epsilc) gas_aqfrac(it,kt,jt,p_hno3) = conv_factor*liquid(lhno3l)/gas(lhno3)
        if (p_h2o2 .gt. 1 .and. gas(lh2o2) .gt. epsilc) gas_aqfrac(it,kt,jt,p_h2o2) = conv_factor*liquid(lh2o2l)/gas(lh2o2)
        if (p_o3   .gt. 1 .and. gas(lo3)   .gt. epsilc) gas_aqfrac(it,kt,jt,p_o3)   = conv_factor*liquid(lo3l)/gas(lo3)


        if (config_flags%chem_opt==CB05_SORG_VBS_AQ_KPP) then
        if (p_facd .gt. 1 .and. gas(lfoa)  .gt. epsilc) gas_aqfrac(it,kt,jt,p_facd) = conv_factor*liquid(lfoal)/gas(lfoa)
        if (p_mepx  .gt. 1 .and. gas(lmhp)  .gt. epsilc) gas_aqfrac(it,kt,jt,p_mepx)  = conv_factor*liquid(lmhpl)/gas(lmhp)
        if (p_pacd  .gt. 1 .and. gas(lpaa)  .gt. epsilc) gas_aqfrac(it,kt,jt,p_pacd)  = conv_factor*liquid(lpaal)/gas(lpaa)

        else
        if (p_ora1 .gt. 1 .and. gas(lfoa)  .gt. epsilc) gas_aqfrac(it,kt,jt,p_ora1) = conv_factor*liquid(lfoal)/gas(lfoa)
        if (p_op1  .gt. 1 .and. gas(lmhp)  .gt. epsilc) gas_aqfrac(it,kt,jt,p_op1)  = conv_factor*liquid(lmhpl)/gas(lmhp)
        if (p_paa  .gt. 1 .and. gas(lpaa)  .gt. epsilc) gas_aqfrac(it,kt,jt,p_paa)  = conv_factor*liquid(lpaal)/gas(lpaa)
        end if




        
      endif
    
    enddo
    enddo
    enddo
    
	end subroutine sorgam_vbs_aqchem_driver
  
end module module_sorgam_vbs_aqchem
