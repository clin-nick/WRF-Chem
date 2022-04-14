





module module_cam_mam_cloudchem

  use shr_kind_mod,       only: r8 => shr_kind_r8
  use module_cam_support, only: pcnst =>pcnst_runtime, pcols, pver, fieldname_len, &
       gas_pcnst => gas_pcnst_modal_aero,iam

  implicit none
  save

  private
  public :: cam_mam_cloudchem_driver
  public :: cam_mam_cloudchem_inti

  integer :: synoz_ndx, so4_ndx, h2o_ndx, o2_ndx, o_ndx, hno3_ndx, dst_ndx, cldice_ndx
  integer :: o3_ndx
  integer :: het1_ndx
  logical :: inv_o3, inv_oh, inv_no3, inv_ho2
  integer :: id_o3, id_oh, id_no3, id_ho2
  integer :: dgnum_idx       = 0
  integer :: dgnumwet_idx    = 0
  integer :: wetdens_ap_idx  = 0

contains

  subroutine cam_mam_cloudchem_inti()
    use mo_setsox, only : sox_inti
    implicit none

    
    call  sox_inti
    
  end subroutine cam_mam_cloudchem_inti
  
  subroutine cam_mam_cloudchem_driver(           &
       
       dvmrdt_sv13d,dvmrcwdt_sv13d,              &
       
       chem,                                     &
       
       moist, scalar, p8w, prain3d, p_phy,       &
       t_phy, dtstepc, ktau,alt, f_ice_phy,      &
       f_rain_phy,cldfra, cldfra_mp_all,         &
       cldfrai, cldfral, is_CAMMGMP_used,        & 
       ids,ide, jds,jde, kds,kde,                &
       ims,ime, jms,jme, kms,kme,                &
       its,ite, jts,jte, kts,kte                 )
    
    
    
    
    
    use module_configure,          only: grid_config_rec_type
    use module_state_description,  only: num_moist, num_chem, p_qv, p_qc,    &
         p_qi, p_qs, p_qnc, p_qni, param_first_scalar, num_scalar, f_qc,     &
         f_qi,  f_qv, f_qs 
    use module_cam_support,        only: pcnst =>pcnst_runtime, pcols, pver, &
         pcnst_non_chem => pcnst_non_chem_modal_aero, nfs,                   &
         gas_pcnst => gas_pcnst_modal_aero,                                  &
         gas_pcnst_pos => gas_pcnst_modal_aero_pos
    use constituents,              only: cnst_get_ind
    use module_data_cam_mam_asect, only: lptr_chem_to_q, lptr_chem_to_qqcw,  &
         factconv_chem_to_q, mw_q_array
    use physconst,                 only: mwdry, avogad

    use mo_setsox,                 only: setsox, has_sox
    use modal_aero_data,           only: ntot_amode
    use infnan,                    only: nan

    
    implicit none

    logical, intent(in) :: is_CAMMGMP_used
    
    integer, intent(in) :: ktau       
    integer, intent(in) ::          &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         its,ite, jts,jte, kts,kte

    real, intent(in) :: dtstepc       

    
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: p8w     
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: prain3d 
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: p_phy   
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: t_phy   
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: alt
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: F_ICE_PHY   
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: F_RAIN_PHY  
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: cldfra   
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: cldfra_mp_all   
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: cldfrai
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: cldfral
    
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme, 1:num_moist )  :: moist  
    real, intent(in), dimension( ims:ime, kms:kme, jms:jme, 1:num_scalar ) :: scalar 
    
    real, intent(inout), dimension( ims:ime, kms:kme, jms:jme, 1:num_chem )   :: chem 
    
    real, intent(out), dimension( ims:ime, kms:kme, jms:jme, 1:gas_pcnst_pos )   :: dvmrdt_sv13d,dvmrcwdt_sv13d 


    
    
    

    integer            :: lchnk                         
    integer            :: ncol                          
    integer            :: imozart                       
    real(r8)           :: delt                          
    real(r8)           :: tfld(pcols,kte)               
    real(r8)           :: pmid(pcols,kte)               
    real(r8)           :: pdel(pcols,kte)               
    real(r8)           :: cldw(pcols,kte)               
    real(r8)           :: ncldwtr(pcols,kte)            

    
    
    
    integer      ::  i, k, m, n
    real(r8)     ::  invariants(pcols,kte,nfs)
    real(r8)     ::  vmr(pcols,kte,gas_pcnst)              
    real(r8)     ::  vmrcw(pcols,kte,gas_pcnst)            
    real(r8), dimension(pcols,kte) :: &
         mbar                                           
    real(r8), dimension(pcols,kte) :: &
         cwat, &                                           
         cldnum, &                                         
         cldfr, &                                          
         prain

    real(r8) :: qh2o(pcols,kte)               
    real(r8) :: dvmrdt_sv1(pcols,pver,gas_pcnst_pos)
    real(r8) :: dvmrcwdt_sv1(pcols,pver,gas_pcnst_pos)

    

    logical, parameter :: do_cam_sulfchem = .FALSE. 

    integer      ::  imozart_m1, icol, itsm1, itile_len
    integer      ::  iw, jw, kw, ktep1, kflip, l, l2, l3
    integer      ::  l_mo_h2so4, l_mo_soag, p1st, ichem
    real(r8)     ::  dp, multFrc, fconv
    real(r8)     ::  xhnm(pcols,kte)
    real(r8)     ::  state_q(pcols,kte,pcnst)
    real(r8)     ::  qqcw(pcols,kte,pcnst)      


    
    delt = dtstepc
    pver = kte


    
    imozart_m1 = pcnst_non_chem
    imozart = imozart_m1 + 1

    
    call cnst_get_ind( 'h2so4', l_mo_h2so4, .false. )
    l_mo_h2so4 = l_mo_h2so4 - imozart_m1
    if ((l_mo_h2so4 < 1) .or. (l_mo_h2so4 > gas_pcnst)) &
         call wrf_error_fatal3("<stdin>",168,&
'cam_mam_cloudchem error -- no h2so4 species' )
    write(*,*) 'l_mo_h2so4 = ', l_mo_h2so4

    call cnst_get_ind( 'soag', l_mo_soag, .false. )
    l_mo_soag = l_mo_soag - imozart_m1
    if ((l_mo_soag < 1) .or. (l_mo_soag > gas_pcnst)) &
         call  wrf_error_fatal3("<stdin>",175,&
'cam_mam_cloudchem error -- no soag species' )
    write(*,*) 'l_mo_soag = ', l_mo_soag


    
    p1st = param_first_scalar 


    ncol = pcols
    icol = ncol 

    
    if(ncol .NE. 1) then
       call wrf_error_fatal3("<stdin>",189,&
'Number of CAM Columns (NCOL) in CAM_MAM_CLOUDCHEM scheme must be 1')
    endif

    
    
    
    itsm1     = its - 1
    itile_len = ite - itsm1
    do jw     = jts , jte
       do iw  = its , ite

          lchnk   = (jw - jts) * itile_len + (iw - itsm1)             
          ktep1   = kte + 1

          
          do kw  = kts, kte
             kflip          = ktep1 - kw
             pmid(1,kflip)  = p_phy(iw,kw,jw)                   
             dp             = p8w(iw,kw,jw) - p8w(iw,kw+1,jw)   
             pdel(1,kflip)  = dp
             tfld(1,kflip)  = t_phy(iw,kw,jw)                   


             
             multFrc              = 1._r8/(1._r8 + moist(iw,kw,jw,P_QV))
             state_q(1,kflip,1)   = max( moist(iw,kw,jw,P_QV)*multFrc, 1.0e-30_r8 ) 
             state_q(1,kflip,2)   = moist(iw,kw,jw,P_QC)*multFrc                    
             state_q(1,kflip,3)   = moist(iw,kw,jw,P_QI)*multFrc                    
             state_q(1,kflip,4)   = scalar(iw,kw,jw,P_QNC)*multFrc                  
             state_q(1,kflip,5)   = scalar(iw,kw,jw,P_QNI)*multFrc                  

             
             qh2o(1,kflip)        = state_q(1,kflip,1)
             cwat(1,kflip)        = state_q(1,kflip,2) + state_q(1,kflip,3)
             cldnum(1,kflip)      = state_q(1,kflip,4)

             
             
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

             prain(1,kflip)        = prain3d(iw,kw,jw)                      
             if(is_CAMMGMP_used) then
                cldfr(1,kflip)        = cldfral(iw,kw,jw)                       
             else
                cldfr(1,kflip)        = cldfra(iw,kw,jw)                        
             endif
             cldfr(1,kflip)        = min(max(cldfr(1,kflip),0._r8),1._r8)
             
             invariants(1,kflip,:) = nan
             mbar(1,kflip)         = mwdry

             
             xhnm(1,kflip) = (avogad*1.0e-6_r8)/(mwdry*alt(iw,kw,jw)) 
             
             dvmrdt_sv13d(iw,kw,jw,:)   = 0.0
             dvmrcwdt_sv13d(iw,kw,jw,:) = 0.0

             
             vmr(icol,kflip,:)          = 0.0_r8
             vmrcw(icol,kflip,:)        = 0.0_r8

             
             
             
             

             
             do ichem = p1st , num_chem
                l2 = lptr_chem_to_q(ichem)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   l3                   = l2 - pcnst_non_chem
                   fconv                = mwdry/mw_q_array(l2)
                   vmr(icol,kflip,l3)   = state_q(icol,kflip,l2)*fconv
                endif

                if (iw*jw == 1 .and. kw == kts .and. l3 == l_mo_h2so4) then
                   write(*,'(a,1p,2e11.3)') '1,1,1 h2so4 q8 & vmr8', state_q(icol,pver,l2), vmr(icol,pver,l3)
                endif
                if (iw*jw == 1 .and. kw == kts .and. l3 == l_mo_soag) then
                   write(*,'(a,1p,2e11.3)') '1,1,1 soag  q8 & vmr8', state_q(icol,pver,l2), vmr(icol,pver,l3)
                endif                

                l2 = lptr_chem_to_qqcw(ichem)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   l3                   = l2 - pcnst_non_chem
                   fconv                = mwdry/mw_q_array(l2)
                   vmrcw(icol,kflip,l3) = qqcw(icol,kflip,l2)*fconv
                endif
             end do

          enddo 

          dvmrdt_sv1 = vmr
          dvmrcwdt_sv1 = vmrcw
          
          
          if( has_sox .and. (.not. do_cam_sulfchem) ) then
             call setsox( ncol,   &
                  pmid,   &
                  delt,   &
                  tfld,   &
                  qh2o,   &
                  cwat,   &
                  lchnk,  &
                  pdel,   &
                  mbar,   &
                  prain,  &
                  cldfr,  &
                  cldnum, &
                  vmrcw,    &
                  imozart-1,&
                  xhnm, & 
                  vmr, &
                  invariants )

          endif
          dvmrdt_sv1 = (vmr - dvmrdt_sv1)/delt
          dvmrcwdt_sv1 = (vmrcw - dvmrcwdt_sv1)/delt


          do l2 = pcnst_non_chem+1, pcnst
             l3                 = l2 - pcnst_non_chem
             fconv              = mw_q_array(l2)/mwdry
             state_q(icol,:,l2) = vmr(icol,:,l3)*fconv
             qqcw(icol,:,l2)    = vmrcw(icol,:,l3)*fconv
             if (iw*jw == 1 .and. kw == kts .and. l3 == l_mo_h2so4) &
                  write(*,'(a,1p,2e11.3)') '1,1,1 h2so4 q8 & vmr8', state_q(icol,pver,l2), vmr(icol,pver,l3)
             if (iw*jw == 1 .and. kw == kts .and. l3 == l_mo_soag) &
                  write(*,'(a,1p,2e11.3)') '1,1,1 soag  q8 & vmr8', state_q(icol,pver,l2), vmr(icol,pver,l3)
          end do


          do kw = kts , kte
             kflip = kte-kw+1

             do l = p1st, num_chem
                l2 = lptr_chem_to_q(l)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   chem(iw,kw,jw,l)         = state_q(1,kflip,l2)/factconv_chem_to_q(l)             
                end if
                l2 = lptr_chem_to_qqcw(l)
                if ((l2 >= 1) .and. (l2 <= pcnst)) then
                   chem(iw,kw,jw,l) = qqcw(1,kflip,l2)/factconv_chem_to_q(l)
                end if
             end do 
             do l = 1 , gas_pcnst
                dvmrdt_sv13d(iw,kw,jw,l)   = dvmrdt_sv1(1,kflip,l) 
                dvmrcwdt_sv13d(iw,kw,jw,l) = dvmrcwdt_sv1(1,kflip,l) 
             enddo

          end do 


       enddo 
    enddo    



  end subroutine cam_mam_cloudchem_driver

end module module_cam_mam_cloudchem
