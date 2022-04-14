

module module_cam_mam_gas_wetdep_driver

  use shr_kind_mod,       only: r8 => shr_kind_r8
  use module_cam_support, only: pcnst =>pcnst_runtime, pcols, pver

  implicit none
  save

  private
  public :: cam_mam_gas_wetdep_driver
  public :: cam_mam_gas_wetdep_inti

contains

  subroutine cam_mam_gas_wetdep_inti()
    use mo_sethet, only : sethet_inti
    implicit none

    
    call  sethet_inti
    
  end subroutine cam_mam_gas_wetdep_inti
  
  subroutine cam_mam_gas_wetdep_driver(          &
       
       chem,                                     &
       
       dtstepc, config_flags, ht, XLAT, nevapr3d,&
       rprdsh3d, rprddp3d, prain3d, z_sea_level, & 
       p_phy, t_phy, alt, moist, scalar,         &
       ids,ide, jds,jde, kds,kde,                &
       ims,ime, jms,jme, kms,kme,                &
       its,ite, jts,jte, kts,kte                 )


    use module_state_description,  only: num_moist, num_chem, p_qv, p_qc,    &
         p_qi, p_qnc, p_qni, param_first_scalar, num_scalar, CAMUWSHCUSCHEME,&
         CAMZMSCHEME

    use module_cam_support,        only: pcnst =>pcnst_runtime, pcols, pver, &
         pcnst_non_chem => pcnst_non_chem_modal_aero, nfs,                   &
         gas_pcnst => gas_pcnst_modal_aero, hetcnt => gas_pcnst_modal_aero 

    use module_data_cam_mam_asect, only: lptr_chem_to_q, lptr_chem_to_qqcw,  &
         factconv_chem_to_q, mw_q_array
    use physconst,                 only: mwdry, avogad, gravit, rga

    use mo_sethet,                 only: sethet
    use infnan,                    only: nan
    use module_configure,          only: grid_config_rec_type

    
    implicit none

    TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

    
    integer, intent(in) ::          &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         its,ite, jts,jte, kts,kte

    real, intent(in) :: dtstepc       

    
    real, dimension( ims:ime, jms:jme ), intent(in) :: ht     
    real, dimension( ims:ime, jms:jme ), intent(in) :: XLAT   

    
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: prain3d     
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: p_phy       
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: t_phy       
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: alt         
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: z_sea_level 
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: nevapr3d    
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: rprdsh3d    
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: rprddp3d    

    
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme, 1:num_moist )  :: moist  
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme, 1:num_scalar ) :: scalar 

    
    real, intent(inout), dimension( ims:ime, kms:kme, jms:jme, 1:num_chem )   :: chem 


    
    
    

    integer            :: lchnk                         
    integer            :: ncol                          
    real(r8)           :: delt                          
    real(r8)           :: tfld(pcols,kte)               
    real(r8)           :: pmid(pcols,kte)               

    
    
    
    real(r8), parameter :: m2km  = 1.e-3_r8

    integer  :: i, k, m, n
    real(r8) :: zsurf(pcols)                          
    real(r8) :: phis(pcols)    
    real(r8) :: cmfdqr(pcols, kte)
    real(r8) :: nevapr(pcols, kte)
    real(r8) :: rprdsh(pcols, kte)
    real(r8) :: rprddp(pcols, kte)
    real(r8) :: state_zm(pcols,kte)
    real(r8) :: zmid(pcols,kte)
    real(r8) :: invariants(pcols,kte,nfs)
    real(r8) :: vmr(pcols,kte,gas_pcnst)              
    real(r8) :: het_rates(pcols,kte,hetcnt)           
    real(r8) :: prain(pcols,kte)

    
    integer  ::  icol, itsm1, itile_len
    integer  ::  iw, jw, kw, ktep1, kflip, l, l2, l3
    integer  ::  p1st, ichem
    real(r8) ::  dp, multFrc, fconv
    real(r8) ::  rlat(pcols) 
    real(r8) ::  xhnm(pcols,kte)
    real(r8) ::  state_q(pcols,kte,pcnst)


    
    delt = dtstepc 
    pver = kte

    
    p1st = param_first_scalar 

    ncol = pcols
    icol = ncol 

    
    if(ncol .NE. 1) then
       call wrf_error_fatal3("<stdin>",140,&
'Number of CAM Columns (NCOL) in CAM_MAM_CLOUDCHEM scheme must be 1')
    endif

    
    
    
    itsm1     = its - 1
    itile_len = ite - itsm1
    do jw     = jts , jte
       do iw  = its , ite

          lchnk   = (jw - jts) * itile_len + (iw - itsm1)             
          ktep1   = kte + 1

          phis(1) = ht(iw,jw)  * gravit 

          
          do kw  = kts, kte
             kflip             = ktep1 - kw
             pmid(1,kflip)     = p_phy(iw,kw,jw)                   
             tfld(1,kflip)     = t_phy(iw,kw,jw)                   
             state_zm(1,kflip) = z_sea_level(iw,kw,jw) - ht(iw,jw) 
             rlat(1)           = xlat(iw,jw)


             
             multFrc              = 1._r8/(1._r8 + moist(iw,kw,jw,P_QV))
             state_q(1,kflip,1)   = max( moist(iw,kw,jw,P_QV)*multFrc, 1.0e-30_r8 ) 
             state_q(1,kflip,2)   = moist(iw,kw,jw,P_QC)*multFrc                    
             state_q(1,kflip,3)   = moist(iw,kw,jw,P_QI)*multFrc                    
             state_q(1,kflip,4)   = scalar(iw,kw,jw,P_QNC)*multFrc                  
             state_q(1,kflip,5)   = scalar(iw,kw,jw,P_QNI)*multFrc                  

             
             
             do l = p1st, num_chem
                l2 = lptr_chem_to_q(l)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   state_q(1,kflip,l2) = chem(iw,kw,jw,l)*factconv_chem_to_q(l)
                end if
             end do 
             prain(1,kflip)        = prain3d(iw,kw,jw) 
             nevapr(1,kflip)       = nevapr3d(iw,kw,jw)                           

             rprdsh(1,kflip)            = 0.0_r8
             
             if(config_flags%shcu_physics==CAMUWSHCUSCHEME) then
                
                rprdsh(1,kflip)         = rprdsh3d(iw,kw,jw)                     
             endif
             
             rprddp(1,kflip)            = 0.0_r8
             if(config_flags%cu_physics==CAMZMSCHEME)then
                
                rprddp(1,kflip)         = rprddp3d(iw,kw,jw)                     
             endif

             
             xhnm(1,kflip) = (avogad*1.0e-6_r8)/(mwdry*alt(iw,kw,jw)) 

             
             invariants(1,kflip,:) = nan

             do ichem = p1st , num_chem
                l2 = lptr_chem_to_q(ichem)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   l3                   = l2 - pcnst_non_chem
                   fconv                = mwdry/mw_q_array(l2)
                   vmr(icol,kflip,l3)   = state_q(icol,kflip,l2)*fconv
                endif
             end do

          enddo 

          zsurf(:ncol) = rga * phis(:ncol)
          do k = 1,pver
             zmid(:ncol,k) = m2km * (state_zm(:ncol,k) + zsurf(:ncol))
          end do


          cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)
          
          call sethet( het_rates, pmid, zmid, phis, tfld, &
                   cmfdqr, prain, nevapr, delt, xhnm, & 
                   vmr, ncol, lchnk,rlat)
          
          do l2 = pcnst_non_chem+1, pcnst
             l3                 = l2 - pcnst_non_chem    
             het_rates(icol,:,l3) = min(max(0._r8, het_rates(icol,:,l3)),(1._r8-1.e-5_r8)/delt) 
             state_q(icol,:,l2) = state_q(icol,:,l2) - het_rates(icol,:,l3)*state_q(icol,:,l2)*delt
          end do


          do kw = kts , kte
             kflip = kte-kw+1

             do l = p1st, num_chem
                l2 = lptr_chem_to_q(l)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   chem(iw,kw,jw,l) = state_q(1,kflip,l2)/factconv_chem_to_q(l)             
                end if
             end do 
          end do 


       enddo 
    enddo    



  end subroutine cam_mam_gas_wetdep_driver

end module module_cam_mam_gas_wetdep_driver
