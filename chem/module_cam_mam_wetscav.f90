
module module_cam_mam_wetscav
  
  
  
  
  

  
  
  

  
  
  

  
  
  
  

  
  
  
  
  
  
  
  
  
  
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  
  implicit none
  private
  save
  
  
  
  
  public :: wetscav_cam_mam_driver_init 
  public :: wetscav_cam_mam_driver      
  
  integer      :: itf, jtf, ktf, pverp
  real(r8)     :: dp1, dp2
  
  
contains
  
  
  subroutine wetscav_cam_mam_driver(itimestep,p_hyd,p8w,t_phy,             &
       dgnum4d,dgnumwet4d,dlf3d,dlf2_3d,dtstep,qme3d,                      &
       prain3d,nevapr3d,rate1ord_cw2pr_st3d,shfrc3d,cmfmc3d,cmfmc2_3d,     &
       evapcsh3d,icwmrsh3d,rprdsh3d,evapcdp3d,icwmrdp3d,rprddp3d,          &
       qs_curr, f_ice_phy,f_rain_phy,config_flags,cldfra_mp_all,cldfrai,   &
       cldfral,cldfra,is_CAMMGMP_used,                                     &
       ids,ide, jds,jde, kds,kde,                                          &
       ims,ime, jms,jme, kms,kme,                                          &
       its,ite, jts,jte, kts,kte,                                          &
       
       qv_curr,qc_curr,qi_curr,ni3d,nc3d,chem,                             &
       
       fracis3D                                                            )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    use module_mp_cammgmp_driver,  only: physics_update, physics_ptend_init
    use module_cam_support,        only: pcnst =>pcnst_runtime, pcols, pver
    use module_state_description,  only: num_chem, param_first_scalar,F_QC, F_QI, F_QV, F_QS, &
         CAMZMSCHEME, CAMUWSHCUSCHEME
    use module_data_cam_mam_asect, only: lptr_chem_to_q, lptr_chem_to_qqcw, factconv_chem_to_q, waterptr_aer
    use wetdep,                    only: clddiag
    use mz_aerosols_intr,          only: mz_aero_wet_intr
    use modal_aero_data,           only: ntot_amode
    use module_configure,          only: grid_config_rec_type
    use infnan,                    only: nan
    
    implicit none
    
    
    
    
    
    TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

    logical, intent(in) :: is_CAMMGMP_used

    integer, intent(in) :: itimestep                           
    integer, intent(in) :: ids,ide, jds,jde, kds,kde
    integer, intent(in) :: ims,ime, jms,jme, kms,kme
    integer, intent(in) :: its,ite, jts,jte, kts,kte

    real, intent(in)    :: dtstep                              
    
    
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: dlf3d          
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: dlf2_3d        
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: p8w            
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: p_hyd          
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: t_phy          
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: qme3d          
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: prain3d        
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: nevapr3d       
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: rate1ord_cw2pr_st3d 
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: shfrc3d        
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: cmfmc3d        
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: cmfmc2_3d      
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: evapcsh3d      
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: icwmrsh3d      
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: rprdsh3d       
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: evapcdp3d      
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: icwmrdp3d      
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: rprddp3d       
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: qs_curr        
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: F_ICE_PHY      
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: F_RAIN_PHY     
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: cldfra_mp_all
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: cldfrai
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: cldfral
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) :: cldfra
    
    
    real, dimension( ims:ime, kms:kme, jms:jme, ntot_amode ), intent(in) :: dgnum4d    
    real, dimension( ims:ime, kms:kme, jms:jme, ntot_amode ), intent(in) :: dgnumwet4d 

    
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(inout) :: qv_curr     
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(inout) :: qc_curr     
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(inout) :: qi_curr     
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(inout) :: ni3d        
    real, dimension( ims:ime, kms:kme, jms:jme ), intent(inout) :: nc3d        


    
    real, dimension( ims:ime, kms:kme, jms:jme, num_chem ),   intent(inout) :: chem     

    
    real, dimension( ims:ime, kms:kme, jms:jme, pcnst ), intent(out) :: fracis3d      

    
    real(r8) :: dt                      
    integer  :: nstep                   
    real(r8) :: clds(pcols,kte)         
    
    integer  i, k
    
    character*24 :: ptend_name            
    
    logical      :: ptend_ls              
    logical      :: ptend_lq(pcnst)       
    
    integer      :: iw, kw, jw, itsm1
    integer      :: itile_len, ktep1
    integer      :: kflip, l, l2, p1st
    integer      :: imode, kcam
    integer      :: ncol                  
    integer      :: lchnk                 
    
    real(r8)     :: dp, multFrc
    
    real(r8)     :: cldc(pcols,kte)           
    real(r8)     :: cldv(pcols,kte)           
    real(r8)     :: cldvcu(pcols,kte)         
    real(r8)     :: cldvst(pcols,kte)         
    real(r8)     :: conicw(pcols, kte)
    real(r8)     :: cmfdqr(pcols, kte)
    real(r8)     :: evapc(pcols, kte)         
    real(r8)     :: rainmr(pcols, kte)        
    real(r8)     :: dlf(pcols,kte)            
    real(r8)     :: dlf2(pcols,kte)           
    real(r8)     :: cmfmc(pcols,kte+1)
    real(r8)     :: cmfmc2(pcols,kte+1)
    real(r8)     :: calday                    
    
    


    real(r8)     :: state_t(pcols,kte)
    real(r8)     :: state_q(pcols,kte,pcnst)
    real(r8)     :: state_pmid(pcols,kte)
    real(r8)     :: state_pdel(pcols,kte)
    
    real(r8)     :: ptend_s(pcols,kte)         
    real(r8)     :: state_s(pcols,kte)         
    
    real(r8)     :: ptend_q(pcols,kte,pcnst)
    
    
    real(r8)     :: cam_out

    
    
    
    real(r8), dimension(pcols,kte)       :: cldn       
    real(r8), dimension(pcols,kte)       :: cme        
    real(r8), dimension(pcols,kte)       :: prain      
    real(r8), dimension(pcols,kte)       :: evapr      
    real(r8), dimension(pcols,kte)       :: icwmrdp    
    real(r8), dimension(pcols,kte)       :: rprddp     
    real(r8), dimension(pcols,kte)       :: icwmrsh    
    real(r8), dimension(pcols,kte)       :: rprdsh     
    real(r8), dimension(pcols,kte,pcnst) :: fracis     
    
    
    real(r8), dimension(pcols,kte)       ::  sh_frac   
    real(r8), dimension(pcols,kte)       ::  dp_frac   
    real(r8), dimension(pcols,kte)       ::  evapcsh   
    real(r8), dimension(pcols,kte)       ::  evapcdp   
    
    
    real(r8), dimension(pcols,kte,ntot_amode) :: dgnum_pbuf, dgnumwet_pbuf 
    
    
    real(r8), dimension(pcols,kte,pcnst)      :: qqcw     
    real(r8), dimension(pcols,kte,ntot_amode) :: qaerwat  
    real(r8), dimension(pcols,kte)            :: rate1ord_cw2pr_st   
    
    nstep = itimestep
    if(itimestep == 1) then
       if(config_flags%shcu_physics .NE. CAMUWSHCUSCHEME) call wrf_message('WARNING: sh_frac,evapcsh,icwmrsh,rprdsh,cmfmc,cmfmc2  are set to zero in CAM_MAM_WETSCAV')
       if(config_flags%cu_physics   .NE. CAMZMSCHEME)     call wrf_message('WARNING: evapcdp,icwmrdp,rprddp,dlf,dlf2 are set to zero in CAM_MAM_WETSCAV')
    endif
    
    call physics_ptend_init(ptend_name,ptend_q,ptend_s,ptend_lq,ptend_ls,pcnst)
    
    
    cme(:,:) = nan
    
    
    
    cam_out = nan
    
    
    

    state_s(:,:) = nan
    ptend_ls     = .FALSE.
    ptend_s(:,:) = nan
    
    
    p1st  = param_first_scalar 
    dt    = real(dtstep,r8)    
    
    ncol  = pcols
    
    if(ncol .NE. 1) then
       call wrf_error_fatal3("<stdin>",259,&
'Number of CAM Columns (NCOL) in CAM_MAM_WETSCAV scheme must be 1')
    endif
    
    
    
    
    calday = nan
    
    
    

    
    itsm1     = its - 1 
    itile_len = ite - itsm1
    do jw     = jts , jte 
       do iw  = its , ite 
          
          lchnk   = (jw - jts) * itile_len + (iw - itsm1)             
          ktep1   = kte + 1
          
          
          do kw  = kts, kte
             kflip                = ktep1 - kw
             
             state_pmid(1,kflip)  = p_hyd(iw,kw,jw)                   
             dp                   = p8w(iw,kw,jw) - p8w(iw,kw+1,jw)   
             state_pdel(1,kflip)  = dp
             state_t(1,kflip)     = t_phy(iw,kw,jw)                   
             

             
             multFrc              = 1._r8/(1._r8 + qv_curr(iw,kw,jw))
             state_q(1,kflip,1)   = max( qv_curr(iw,kw,jw)*multFrc, 1.0e-30_r8 ) 
             state_q(1,kflip,2)   = qc_curr(iw,kw,jw)*multFrc                    
             state_q(1,kflip,3)   = qi_curr(iw,kw,jw)*multFrc                    
             state_q(1,kflip,4)   = nc3d(iw,kw,jw)*multFrc                       
             state_q(1,kflip,5)   = ni3d(iw,kw,jw)*multFrc                       
             
             
             
             do l = p1st, num_chem
                l2 = lptr_chem_to_q(l)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   state_q(1,kflip,l2) = chem(iw,kw,jw,l)*factconv_chem_to_q(l)
                end if
                l2 = lptr_chem_to_qqcw(l)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   qqcw(1,kflip,l2) = chem(iw,kw,jw,l)*factconv_chem_to_q(l)     
                end if
             end do 
             
             
             
             do imode = 1 , ntot_amode
                dgnum_pbuf(1,kflip,imode)    = dgnum4D(iw,kw,jw,imode)           
                dgnumwet_pbuf(1,kflip,imode) = dgnumwet4D(iw,kw,jw,imode) 
                
                l = waterptr_aer(1,imode)
                if ((l >= p1st) .and. (l <= num_chem)) then
                   qaerwat(1,kflip,imode) = chem(iw,kw,jw,l)*factconv_chem_to_q(l)
                endif
             enddo
             

             
             
             cme(1,kflip)         = nan                                          
             prain(1,kflip)       = prain3d(iw,kw,jw)                            
             evapr(1,kflip)       = nevapr3d(iw,kw,jw)                           
             
             rate1ord_cw2pr_st(1,kflip) = rate1ord_cw2pr_st3d(iw,kw,jw)          
             if(is_CAMMGMP_used) then
                cldn(1,kflip)              = cldfra_mp_all(iw,kw,jw)                
             else
                cldn(1,kflip)              = cldfra(iw,kw,jw)
             endif
             cldn(1,kflip)              = min(max(cldn(1,kflip),0._r8),1._r8)
             
             sh_frac(1,kflip)           = 0.0_r8
             evapcsh(1,kflip)           = 0.0_r8
             icwmrsh(1,kflip)           = 0.0_r8
             rprdsh(1,kflip)            = 0.0_r8
             
             if(config_flags%shcu_physics==CAMUWSHCUSCHEME) then
                
                sh_frac(1,kflip)        = shfrc3d(iw,kw,jw)                      
                evapcsh(1,kflip)        = evapcsh3d(iw,kw,jw)                    
                icwmrsh(1,kflip)        = icwmrsh3d(iw,kw,jw)                    
                rprdsh(1,kflip)         = rprdsh3d(iw,kw,jw)                     
             endif
             
             evapcdp(1,kflip)           = 0.0_r8
             icwmrdp(1,kflip)           = 0.0_r8
             rprddp(1,kflip)            = 0.0_r8
             dlf(1,kflip)               = 0.0_r8
             dlf2(1,kflip)              = 0.0_r8
             
             if(config_flags%cu_physics==CAMZMSCHEME)then
                
                evapcdp(1,kflip)        = evapcdp3d(iw,kw,jw)                    
                icwmrdp(1,kflip)        = icwmrdp3d(iw,kw,jw)                    
                rprddp(1,kflip)         = rprddp3d(iw,kw,jw)                     
                dlf(1,kflip)            = dlf3d(iw,kw,jw)                        
                dlf2(1,kflip)           = dlf2_3d(iw,kw,jw)                      
             endif
          enddo
          
          do kw = kts, kte+1
             kflip = kte - kw + 2
             
             cmfmc(1,kflip)      = 0.0_r8
             cmfmc2(1,kflip)     = 0.0_r8
             if(config_flags%shcu_physics==CAMUWSHCUSCHEME) then
                cmfmc(1,kflip)   = cmfmc3d(iw,kw,jw)    
                cmfmc2(1,kflip)  = cmfmc2_3d(iw,kw,jw)  
             endif
          end do
          
          do kcam = 1, kte
             
             dp_frac(1,kcam)         = max(0.0_r8,min(dp1*log(1.0_r8+dp2*(cmfmc(1,kcam+1)-cmfmc2(1,kcam+1))),0.60_r8))
          end do
          
          
          cldc(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:) 
          evapc(:ncol,:) = evapcsh(:ncol,:) + evapcdp(:ncol,:) 
          clds(:ncol,:)  = cldn(:ncol,:) - cldc(:ncol,:)       
          
          
          
          conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
               max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
          
          cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)
          
          
          
          
          
          call clddiag( state_t, state_pmid, state_pdel, cmfdqr, evapc, cldn, cldc, clds, cme, evapr, prain, &
               cldv, cldvcu, cldvst, rainmr, ncol )
          
          ptend_name = 'wetdep'
          
          
          
          
          
          

          
          fracis(:,:,:) = 1.0_r8
          call mz_aero_wet_intr (lchnk, ncol, state_q,                 &
               state_pdel, state_pmid, state_t, ptend_name,            &
               ptend_lq, ptend_q, nstep, dt, cme, prain, evapr, cldv,  &
               cldvcu, cldvst, cldc, cldn, fracis, calday, cmfdqr,     &
               evapc, conicw, rainmr,                                  &
               rate1ord_cw2pr_st,                                      &   
               dgnumwet_pbuf, qqcw, qaerwat, cam_out, dlf              )

          call physics_update(lchnk,dt,state_q,ptend_q,state_s,ptend_s,ptend_name,ptend_lq,ptend_ls, pcnst)
          

          do kw=kts,kte
             kflip = kte-kw+1
             do imode = 1,  ntot_amode
                l = waterptr_aer(1,imode)
                if ((l >= p1st) .and. (l <= num_chem)) then
                   chem(iw,kw,jw,l) = qaerwat(1,kflip,imode)/factconv_chem_to_q(l)
                endif
             end do 
             
             
             qv_curr(iw,kw,jw)       = state_q(1,kflip,1) / (1.0_r8 - state_q(1,kflip,1)) 
             multFrc                 = 1._r8 + qv_curr(iw,kw,jw)
             
             qc_curr(iw,kw,jw)       = state_q(1,kflip,2) * multFrc
             qi_curr(iw,kw,jw)       = state_q(1,kflip,3) * multFrc 
             nc3d(iw,kw,jw)          = state_q(1,kflip,4) * multFrc  
             ni3d(iw,kw,jw)          = state_q(1,kflip,5) * multFrc
             do l = 1 ,5
                fracis3d(iw,kw,jw,l)     = fracis(1,kflip,l)          
             enddo
             do l = p1st, num_chem
                l2 = lptr_chem_to_q(l)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   chem(iw,kw,jw,l) = state_q(1,kflip,l2)/factconv_chem_to_q(l)
                   fracis3d(iw,kw,jw,l2)     = fracis(1,kflip,l2)          
                end if
                l2 = lptr_chem_to_qqcw(l)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   chem(iw,kw,jw,l) = qqcw(1,kflip,l2)/factconv_chem_to_q(l)
                end if
             end do 
          end do
          
       enddo 
    enddo 
    return
    
  end subroutine wetscav_cam_mam_driver
  
  
  subroutine wetscav_cam_mam_driver_init(ids,ide, jds,jde, kds,kde, &
       ims,ime, jms,jme, kms,kme,                           &
       its,ite, jts,jte, kts,kte                            )
    
    
    
    
    
    
    use module_cam_support,        only: pver
    use mz_aerosols_intr,          only: modal_aero_bcscavcoef_init, mz_aero_initialize 
    implicit none
    integer, intent(in) :: ids,ide, jds,jde, kds,kde
    integer, intent(in) :: ims,ime, jms,jme, kms,kme
    integer, intent(in) :: its,ite, jts,jte, kts,kte
    
    jtf   = min(jte,jde-1)
    ktf   = min(kte,kde-1)
    itf   = min(ite,ide-1)

    
    pver  = ktf - kts + 1 
    pverp = pver + 1

    
    dp1   = 0.10_r8 
    dp2   = 500.0_r8
    
    
    call  mz_aero_initialize 
    
  end subroutine wetscav_cam_mam_driver_init
  
end module module_cam_mam_wetscav
