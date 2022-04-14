
module mz_aerosols_intr

  use shr_kind_mod,only: r8 => shr_kind_r8
  use module_cam_support, only: masterproc, pcnst =>pcnst_runtime, pcols, pver, endrun
  use constituents,only: cnst_name, cnst_get_ind 
    use modal_aero_data, only: ntot_amode
  use module_cam_support, only: iulog

  implicit none

  private          

  save

  
  
  

  public :: mz_aero_initialize                           
  public :: mz_aero_wet_intr                             
  public :: aer_wetdep_list
  public :: do_cam_sulfchem
  public :: modal_aero_bcscavcoef_init  
  public :: modal_aero_bcscavcoef_get   

  integer :: num_mz_aerosols
  integer, allocatable :: mz_aerosol_ids(:) 
  character(len=16), allocatable, dimension(:) :: aer_wetdep_list
  logical :: use_cam_sulfchem = .false.
  logical :: do_cam_sulfchem
  logical :: inv_o3, inv_oh, inv_no3, inv_ho2
  integer, pointer :: id_so2, id_so4, id_dms, id_o3, id_h2o2, id_oh, id_no3, id_ho2
  integer, target  :: spc_ids(8)

  character(len=16), allocatable,dimension(:) :: aer_drydep_lis

 integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12
 real(r8) dlndg_nimptblgrow
 real(r8) scavimptblnum(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)
 real(r8) scavimptblvol(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)

contains
  
  subroutine mz_aero_initialize( )
    use module_cam_support, only: addfld, add_default, phys_decomp

    implicit none

    integer :: m                                
    integer :: mm                               
    integer :: astat, id

    logical :: history_aerosol      

    
    history_aerosol = .FALSE.

    num_mz_aerosols = 0
    

    
    num_mz_aerosols = pcnst
  end subroutine mz_aero_initialize

  
  subroutine mz_aero_wet_intr ( lchnk_in, ncol_in, state_q, state_pdel,          &
       state_pmid, state_t, ptend_name, ptend_lq, ptend_q, nstep, dt, cme, prain,    &
       evapr, cldv, cldvcu, cldvst, cldc, cldn, fracis, calday, cmfdqr, evapc, conicw, rainmr  &
     , rate1ord_cw2pr_st                                                & 
     , dgncur_awet,  qqcw, qaerwat                                       &
     , cam_out , dlf                                                         )
    
    
    use module_cam_support, only: outfld
    use wetdep,        only: wetdepa_v1, wetdepa_v2
    use physconst,     only: gravit
    use constituents,  only: cnst_mw
    use physconst,     only: amass => mwdry         
    use physconst,     only: boltz                  
    use modal_aero_data

    
    implicit none
    
    
    
    
    real(r8),            intent(in)  :: dt             
    integer,  intent(in) :: lchnk_in,ncol_in
    real(r8), intent(in) :: state_q(pcols,pver,pcnst) 
    real(r8), intent(in) :: state_pdel(pcols,pver)
    real(r8), intent(in) :: state_t(pcols,pver)
    real(r8), intent(in) :: state_pmid(pcols,pver) 
    integer,  intent(in) :: nstep
    real(r8), intent(in) :: cme(pcols,pver)            
    real(r8), intent(in) :: prain(pcols,pver)          
    real(r8), intent(in) :: evapr(pcols,pver)          
    real(r8), intent(in) :: cldn(pcols,pver)           
    real(r8), intent(in) :: cldc(pcols,pver)           
    real(r8), intent(in) :: cldv(pcols,pver)           
    real(r8), intent(in) :: cldvcu(pcols,pver)         
    real(r8), intent(in) :: cldvst(pcols,pver)         
    real(r8), intent(in) :: conicw(pcols, pver)
    real(r8), intent(in) :: cmfdqr(pcols, pver)
    real(r8), intent(in) :: evapc(pcols, pver)         
    real(r8), intent(in) :: rainmr(pcols, pver)         
    real(r8), intent(in) :: dlf(pcols,pver)            
    character(*), intent(inout) :: ptend_name 
    logical,      intent(inout) :: ptend_lq(pcnst)
    real(r8),     intent(inout) :: ptend_q(pcols,pver,pcnst) 
    real(r8), target, intent(inout) :: qqcw(pcols,pver,pcnst) 
    real(r8),            intent(inout) :: fracis(pcols,pver,pcnst) 
    real(r8), intent(in) :: calday        
    real(r8), intent(in) :: rate1ord_cw2pr_st(pcols,pver) 
    real(r8), intent(in) :: dgncur_awet(pcols,pver,ntot_amode) 
                                                               
    real(r8), intent(inout) :: qaerwat(pcols,pver,ntot_amode)

    real(r8), intent(in) :: cam_out

    
    
    
    integer  :: m                                  
    integer  :: ixcldice, ixcldliq
    integer  :: lchnk                              
    integer  :: ncol                               
    real(r8) :: obuf(1)
    real(r8) :: iscavt(pcols, pver)
    real(r8) :: totcond(pcols, pver) 
    integer  :: yr, mon, day, ncsec
    integer  :: ncdate
    integer  :: mm
    integer  :: nphob
    real(r8), dimension(pcols,pver) :: h2o23d,o3,oh,no3,ho2,so2,so4,dms,h2o2,ekso2,ekh2o2

    real(r8), dimension(pcols,pver) ::&
         ph,          &
         hion          

    integer ::&
         indcp(pcols,pver),    &
         ncldypts(pver)         

    real(r8) :: icwmr2(pcols,pver) 
    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factic(pcols,pver)
    real(r8) :: sol_factbi, sol_factii, sol_factiic
    real(r8) :: xhnm(pcols, pver)
    real(r8) :: sflx(pcols)            
    real(r8) :: dicorfac(pcols)        
    integer :: i,k
    real(r8) :: scavcoef(pcols,pver) 
    integer :: jnv                     
    integer :: lphase                  
    integer :: lspec                   
    integer :: lcoardust, lcoarnacl    
    real(r8) :: dqdt_tmp(pcols,pver)   
    real(r8) :: f_act_conv(pcols,pver) 
    real(r8) :: f_act_conv_coarse(pcols,pver) 
    real(r8) :: f_act_conv_coarse_dust, f_act_conv_coarse_nacl                                         
    real(r8) :: fracis_cw(pcols,pver)
    real(r8) :: hygro_sum_old(pcols,pver)  
    real(r8) :: hygro_sum_del(pcols,pver)  
    real(r8) :: hygro_sum_old_ik, hygro_sum_new_ik
    real(r8) :: prec(pcols)                
    real(r8) :: q_tmp(pcols,pver)          
    real(r8) :: qqcw_tmp(pcols,pver)       
    real(r8) :: scavcoefnv(pcols,pver,0:2) 
                                           
                                           
    real(r8) :: tmpa, tmpb
    real(r8) :: tmpdust, tmpnacl
    real(r8) :: water_old, water_new   
    logical :: isprx(pcols,pver) 
    real(r8) :: aerdepwetis(pcols,pcnst)  
    real(r8) :: aerdepwetcw(pcols,pcnst)  
    real(r8), pointer :: fldcw(:,:)

    
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    lchnk = lchnk_in
    ncol  = ncol_in

    totcond(:ncol, :) = state_q(:ncol,:,ixcldliq) + state_q(:ncol,:,ixcldice)

    if ( num_mz_aerosols > 0 ) then

       
       ptend_name  = ptend_name//'+mz_aero_wetdep'

     prec(:ncol)=0._r8
     do k=1,pver
        where (prec(:ncol) >= 1.e-7)
            isprx(:ncol,k) = .true.
        elsewhere
            isprx(:ncol,k) = .false.
        endwhere
        prec(:ncol) = prec(:ncol) + (prain(:ncol,k) + cmfdqr(:ncol,k) - evapr(:ncol,k)) &
                    *state_pdel(:ncol,k)/gravit
     end do




     f_act_conv_coarse(:,:) = 0.60_r8   
     f_act_conv_coarse_dust = 0.40_r8   
     f_act_conv_coarse_nacl = 0.80_r8   
     if (modeptr_coarse > 0) then
        lcoardust = lptr_dust_a_amode(modeptr_coarse)
        lcoarnacl = lptr_nacl_a_amode(modeptr_coarse)
        if ((lcoardust > 0) .and. (lcoarnacl > 0)) then
           do k = 1, pver
           do i = 1, ncol
              tmpdust = max( 0.0_r8, state_q(i,k,lcoardust) + ptend_q(i,k,lcoardust)*dt )
              tmpnacl = max( 0.0_r8, state_q(i,k,lcoarnacl) + ptend_q(i,k,lcoarnacl)*dt )
              if ((tmpdust+tmpnacl) > 1.0e-30_r8) then

                 f_act_conv_coarse(i,k) = (f_act_conv_coarse_dust*tmpdust &
                                         + f_act_conv_coarse_nacl*tmpnacl)/(tmpdust+tmpnacl)  
              end if
           end do
           end do
        end if
     end if


    scavcoefnv(:,:,0) = 0.0_r8   

    do m = 1, ntot_amode   

       do lphase = 1, 2   






















          if (lphase == 1) then   
             hygro_sum_old(:,:) = 0.0_r8
             hygro_sum_del(:,:) = 0.0_r8
             call modal_aero_bcscavcoef_get( m, ncol, isprx, dgncur_awet,   &
                                             scavcoefnv(:,:,1), scavcoefnv(:,:,2) )

             sol_factb  = 0.1_r8   


             sol_facti  = 0.0_r8   

             sol_factic = 0.4_r8      

             if (m == modeptr_pcarbon) then

                f_act_conv = 0.0_r8   
             else if ((m == modeptr_finedust) .or. (m == modeptr_coardust)) then

                f_act_conv = 0.4_r8   
             else

                f_act_conv = 0.8_r8   
             end if

          else   
             sol_factb  = 0.0_r8   

             sol_facti  = 1.0_r8   

             sol_factic = 0.0_r8   
                                   
             f_act_conv = 0.0_r8   

          end if














          sol_factbi  = sol_factb
          sol_factii  = sol_facti
          sol_factiic = sol_factic(1,1)


          do lspec = 0, nspec_amode(m)+1   

             if (lspec == 0) then   
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jnv = 1
                else
                   mm = numptrcw_amode(m)
                   jnv = 0
                endif
             else if (lspec <= nspec_amode(m)) then   
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jnv = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   jnv = 0
                endif
             else   

                cycle
                if (lphase == 1) then
                   mm = 0

                   jnv = 2
                else
                   mm = 0
                   jnv = 0
                endif
             endif

          if (mm <= 0) cycle








          if ((lphase == 1) .and. (m == modeptr_coarse)) then

             f_act_conv = f_act_conv_coarse   
             if (lspec > 0) then
                if (lmassptr_amode(lspec,m) == lptr_dust_a_amode(m)) then

                   f_act_conv = f_act_conv_coarse_dust   
                else if (lmassptr_amode(lspec,m) == lptr_nacl_a_amode(m)) then

                   f_act_conv = f_act_conv_coarse_nacl   
                end if
             end if
          end if
          if ((lphase == 1) .and. (lspec <= nspec_amode(m))) then
             ptend_lq(mm) = .TRUE.
             dqdt_tmp(:,:) = 0.0_r8

             q_tmp(1:ncol,:) = state_q(1:ncol,:,mm) + ptend_q(1:ncol,:,mm)*dt
             fldcw => qqcw(:,:,mm)
             call wetdepa_v2( state_t, state_pmid, state_q, state_pdel,  &
                  cldn, cldc, cmfdqr, evapc, conicw, prain, cme,                     &
                  evapr, totcond, q_tmp(:,:), dt,            &
                  dqdt_tmp(:,:), iscavt, cldv, cldvcu, cldvst, dlf, fracis(:,:,mm), sol_factb, ncol, &
                  scavcoefnv(:,:,jnv), &
                  .false., rate1ord_cw2pr_st, fldcw, f_act_conv, &     
                  icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                  sol_facti_in=sol_facti, sol_factbi_in=sol_factbi, sol_factii_in=sol_factii, &   
                  sol_factic_in=sol_factic, sol_factiic_in=sol_factiic )                          

             ptend_q(1:ncol,:,mm) = ptend_q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

             call outfld( trim(cnst_name(mm))//'WET', dqdt_tmp(:,:), pcols, lchnk)
             call outfld( trim(cnst_name(mm))//'SIC', icscavt, pcols, lchnk)
             call outfld( trim(cnst_name(mm))//'SIS', isscavt, pcols, lchnk)
             call outfld( trim(cnst_name(mm))//'SBC', bcscavt, pcols, lchnk)
             call outfld( trim(cnst_name(mm))//'SBS', bsscavt, pcols, lchnk)

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+dqdt_tmp(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name(mm))//'SFWET', sflx, pcols, lchnk)
             aerdepwetis(:ncol,mm) = sflx(:ncol)

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+icscavt(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name(mm))//'SFSIC', sflx, pcols, lchnk)
             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+isscavt(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name(mm))//'SFSIS', sflx, pcols, lchnk)
             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+bcscavt(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name(mm))//'SFSBC', sflx, pcols, lchnk)
             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+bsscavt(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name(mm))//'SFSBS', sflx, pcols, lchnk)

             if (lspec > 0) then
                tmpa = spechygro(lspectype_amode(lspec,m))/ &
                       specdens_amode(lspectype_amode(lspec,m))
                tmpb = tmpa*dt
                hygro_sum_old(1:ncol,:) = hygro_sum_old(1:ncol,:) &
                                         + tmpa*q_tmp(1:ncol,:)
                hygro_sum_del(1:ncol,:) = hygro_sum_del(1:ncol,:) &
                                         + tmpb*dqdt_tmp(1:ncol,:)
             end if

          else if ((lphase == 1) .and. (lspec == nspec_amode(m)+1)) then









             do k = 1, pver
             do i = 1, ncol

                water_old = max( 0.0_r8, qaerwat(i,k,mm) )
                hygro_sum_old_ik = max( 0.0_r8, hygro_sum_old(i,k) )
                hygro_sum_new_ik = max( 0.0_r8, hygro_sum_old_ik+hygro_sum_del(i,k) )
                if (hygro_sum_new_ik >= 10.0_r8*hygro_sum_old_ik) then
                   water_new = 10.0_r8*water_old
                else
                   water_new = water_old*(hygro_sum_new_ik/hygro_sum_old_ik)
                end if

                qaerwat(i,k,mm) = water_new
             end do
             end do













          else   
             dqdt_tmp(:,:) = 0.0_r8
             qqcw_tmp(:,:) = 0.0_r8    
             fldcw => qqcw(:,:,mm)
             call wetdepa_v2( state_t, state_pmid, state_q, state_pdel,  &
                  cldn, cldc, cmfdqr, evapc, conicw, prain, cme,                     &
                  evapr, totcond, fldcw, dt,            &
                  dqdt_tmp(:,:), iscavt, cldv, cldvcu, cldvst, dlf, fracis_cw(:,:), sol_factb, ncol, &
                  scavcoefnv(:,:,jnv), &
                  .true., rate1ord_cw2pr_st, qqcw_tmp, f_act_conv, &     
                  icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                  sol_facti_in=sol_facti, sol_factbi_in=sol_factbi, sol_factii_in=sol_factii, &   
                  sol_factic_in=sol_factic, sol_factiic_in=sol_factiic )                          

             fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+dqdt_tmp(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name_cw(mm))//'SFWET', sflx, pcols, lchnk)
             aerdepwetcw(:ncol,mm) = sflx(:ncol)

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+icscavt(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name_cw(mm))//'SFSIC', sflx, pcols, lchnk)
             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+isscavt(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name_cw(mm))//'SFSIS', sflx, pcols, lchnk)
             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+bcscavt(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name_cw(mm))//'SFSBC', sflx, pcols, lchnk)
             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+bsscavt(i,k)*state_pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(cnst_name_cw(mm))//'SFSBS', sflx, pcols, lchnk)

          endif

          enddo   
       enddo   
    enddo   

    endif
    return

  end subroutine mz_aero_wet_intr



  
  subroutine modal_aero_bcscavcoef_get( m, ncol, isprx, dgn_awet, scavcoefnum, scavcoefvol )

    use modal_aero_data
    
    implicit none

	  integer,intent(in) :: m, ncol
	  logical,intent(in):: isprx(pcols,pver)
	  real(r8), intent(in) :: dgn_awet(pcols,pver,ntot_amode)
	  real(r8), intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)

	  integer i, k, jgrow
	  real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum


      do k = 1, pver
      do i = 1, ncol


         if ( isprx(i,k)  ) then



	    dumdgratio = dgn_awet(i,k,m)/dgnum_amode(m)

	    if ((dumdgratio .ge. 0.99) .and. (dumdgratio .le. 1.01)) then
	        scavimpvol = scavimptblvol(0,m)
	        scavimpnum = scavimptblnum(0,m)
	    else
	        xgrow = log( dumdgratio ) / dlndg_nimptblgrow
	        jgrow = int( xgrow )
	        if (xgrow .lt. 0.) jgrow = jgrow - 1
	        if (jgrow .lt. nimptblgrow_mind) then
		    jgrow = nimptblgrow_mind
		    xgrow = jgrow
	        else
	            jgrow = min( jgrow, nimptblgrow_maxd-1 )
	        end if

	        dumfhi = xgrow - jgrow
	        dumflo = 1. - dumfhi

	        scavimpvol = dumflo*scavimptblvol(jgrow,m) +   &
			  dumfhi*scavimptblvol(jgrow+1,m)
	        scavimpnum = dumflo*scavimptblnum(jgrow,m) +   &
			  dumfhi*scavimptblnum(jgrow+1,m)

	    end if


	    scavcoefvol(i,k) = exp( scavimpvol  )

	    scavcoefnum(i,k) = exp( scavimpnum  )






          else
            scavcoefvol(i,k) = 0._r8
	    scavcoefnum(i,k) = 0._r8
          end if

      end do
      end do

        return
	end subroutine modal_aero_bcscavcoef_get


  
  subroutine modal_aero_bcscavcoef_init









  use shr_kind_mod,only: r8 => shr_kind_r8
  use modal_aero_data
  use module_cam_support, only: endrun

  implicit none



	integer nnfit_maxd
	parameter (nnfit_maxd=27)

	integer i, jgrow, jdens, jpress, jtemp, ll, mode, nnfit
        integer lunerr

	real(r8) dg0, dg0_cgs, press, &
	rhodryaero, rhowetaero, rhowetaero_cgs, rmserr, &
	scavratenum, scavratevol, sigmag,                &
	temp, wetdiaratio, wetvolratio
	real(r8) aafitnum(1), xxfitnum(1,nnfit_maxd), yyfitnum(nnfit_maxd)
	real(r8) aafitvol(1), xxfitvol(1,nnfit_maxd), yyfitvol(nnfit_maxd)
        
	lunerr = 6
	dlndg_nimptblgrow = log( 1.25d00 )

	do 7900 mode = 1, ntot_amode

	sigmag = sigmag_amode(mode)

	ll = lspectype_amode(1,mode)
	rhodryaero = specdens_amode(ll)

	do 7800 jgrow = nimptblgrow_mind, nimptblgrow_maxd

	wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
	dg0 = dgnum_amode(mode)*wetdiaratio

	wetvolratio = exp( jgrow*dlndg_nimptblgrow*3. )
	rhowetaero = 1.0 + (rhodryaero-1.0)/wetvolratio
	rhowetaero = min( rhowetaero, rhodryaero )




	nnfit = 0

	temp = 273.16
	press = 0.75e6   
	rhowetaero = rhodryaero

	dg0_cgs = dg0*1.0e2   
	rhowetaero_cgs = rhowetaero*1.0e-3   
	call calc_1_impact_rate( &
     		dg0_cgs, sigmag, rhowetaero_cgs, temp, press, &
     		scavratenum, scavratevol, lunerr )

	nnfit = nnfit + 1
	if (nnfit .gt. nnfit_maxd) then
	    write(lunerr,9110)
	    call endrun()
	end if
9110	format( '*** subr. modal_aero_bcscavcoef_init -- nnfit too big' )

	xxfitnum(1,nnfit) = 1.
	yyfitnum(nnfit) = log( scavratenum )

	xxfitvol(1,nnfit) = 1.
	yyfitvol(nnfit) = log( scavratevol )

5900	continue
















	scavimptblnum(jgrow,mode) = yyfitnum(1)
	scavimptblvol(jgrow,mode) = yyfitvol(1)



7800	continue
7900	continue

        return
        end subroutine modal_aero_bcscavcoef_init


  
	subroutine calc_1_impact_rate(             &
     		dg0, sigmag, rhoaero, temp, press, &
     		scavratenum, scavratevol, lunerr )













     use shr_kind_mod, only: r8 => shr_kind_r8

	implicit none


	integer lunerr
	real(r8) dg0, sigmag, rhoaero, temp, press, scavratenum, scavratevol


	integer nrainsvmax
	parameter (nrainsvmax=50)
	real(r8) rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
     		vfallrainsv(nrainsvmax)

	integer naerosvmax
	parameter (naerosvmax=51)
	real(r8) aaerosv(naerosvmax), &
     	ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

	integer i, ja, jr, na, nr
	real(r8) a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
     	real(r8) anumsum, avolsum, boltzmann, cair, chi
     	real(r8) d, dr, dum, dumfuchs, dx
     	real(r8) ebrown, eimpact, eintercept, etotal, freepath, gravity 
     	real(r8) pi, precip, precipmmhr, precipsum
     	real(r8) r, rainsweepout, reynolds, rhi, rhoair, rhowater, rlo, rnumsum
     	real(r8) scavsumnum, scavsumnumbb
     	real(r8) scavsumvol, scavsumvolbb
     	real(r8) schmidt, sqrtreynolds, sstar, stokes, sx              
     	real(r8) taurelax, vfall, vfallstp
     	real(r8) x, xg0, xg3, xhi, xlo, xmuwaterair                     


	rlo = .005
	rhi = .250
	dr = 0.005
	nr = 1 + nint( (rhi-rlo)/dr )
	if (nr .gt. nrainsvmax) then
	    write(lunerr,9110)
	    call endrun()
	end if
9110	format( '*** subr. calc_1_impact_rate -- nr > nrainsvmax' )

	precipmmhr = 1.0
	precip = precipmmhr/36000.

	ag0 = dg0/2.
	sx = log( sigmag )
	xg0 = log( ag0 )
	xg3 = xg0 + 3.*sx*sx

	xlo = xg3 - 4.*sx
	xhi = xg3 + 4.*sx
	dx = 0.2*sx

	dx = max( 0.2_r8*sx, 0.01_r8 )
	xlo = xg3 - max( 4._r8*sx, 2._r8*dx )
	xhi = xg3 + max( 4._r8*sx, 2._r8*dx )

	na = 1 + nint( (xhi-xlo)/dx )
	if (na .gt. naerosvmax) then
	    write(lunerr,9120)
	    call endrun()
	end if
9120	format( '*** subr. calc_1_impact_rate -- na > naerosvmax' )


	cair = press/(8.31436e7*temp)

	rhoair = 28.966*cair

	freepath = 2.8052e-10/cair

	boltzmann = 1.3807e-16

	rhowater = 1.0

	gravity = 980.616

	airdynvisc = 1.8325e-4 * (416.16/(temp+120.)) *    &
                                        ((temp/296.16)**1.5)

	airkinvisc = airdynvisc/rhoair

	xmuwaterair = 60.0

	pi = 3.1415926536








	precipsum = 0.
	do i = 1, nr
	    r = rlo + (i-1)*dr
	    rrainsv(i) = r
	    xnumrainsv(i) = exp( -r/2.7e-2 )

	    d = 2.*r
	    if (d .le. 0.007) then
		vfallstp = 2.88e5 * d**2.
	    else if (d .le. 0.025) then
		vfallstp = 2.8008e4 * d**1.528
	    else if (d .le. 0.1) then
		vfallstp = 4104.9 * d**1.008
	    else if (d .le. 0.25) then
		vfallstp = 1812.1 * d**0.638
	    else
		vfallstp = 1069.8 * d**0.235
	    end if

	    vfall = vfallstp * sqrt(1.204e-3/rhoair)
	    vfallrainsv(i) = vfall
	    precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
	end do
	precipsum = precipsum*pi*1.333333

	rnumsum = 0.
	do i = 1, nr
	    xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
	    rnumsum = rnumsum + xnumrainsv(i)
	end do







	anumsum = 0.
	avolsum = 0.
	do i = 1, na
	    x = xlo + (i-1)*dx
	    a = exp( x )
	    aaerosv(i) = a
	    dum = (x - xg0)/sx
	    ynumaerosv(i) = exp( -0.5*dum*dum )
	    yvolaerosv(i) = ynumaerosv(i)*1.3333*pi*a*a*a
	    anumsum = anumsum + ynumaerosv(i)
	    avolsum = avolsum + yvolaerosv(i)
	end do

	do i = 1, na
	    ynumaerosv(i) = ynumaerosv(i)/anumsum
	    yvolaerosv(i) = yvolaerosv(i)/avolsum
	end do





	scavsumnum = 0.
	scavsumvol = 0.



	do 5900 jr = 1, nr

	r = rrainsv(jr)
	vfall = vfallrainsv(jr)

	reynolds = r * vfall / airkinvisc
	sqrtreynolds = sqrt( reynolds )




	scavsumnumbb = 0.
	scavsumvolbb = 0.

	do 5500 ja = 1, na

	a = aaerosv(ja)

	chi = a/r

	dum = freepath/a
	dumfuchs = 1. + 1.246*dum + 0.42*dum*exp(-0.87/dum)
	taurelax = 2.*rhoaero*a*a*dumfuchs/(9.*rhoair*airkinvisc)

	aeromass = 4.*pi*a*a*a*rhoaero/3.
	aerodiffus = boltzmann*temp*taurelax/aeromass

	schmidt = airkinvisc/aerodiffus
	stokes = vfall*taurelax/r

	ebrown = 4.*(1. + 0.4*sqrtreynolds*(schmidt**0.3333333)) /  &
     			(reynolds*schmidt)

	dum = (1. + 2.*xmuwaterair*chi) /         &
     			(1. + xmuwaterair/sqrtreynolds)
	eintercept = 4.*chi*(chi + dum)

	dum = log( 1. + reynolds )
	sstar = (1.2 + dum/12.) / (1. + dum)
	eimpact = 0.
	if (stokes .gt. sstar) then
	    dum = stokes - sstar
	    eimpact = (dum/(dum+0.6666667)) ** 1.5
	end if

	etotal = ebrown + eintercept + eimpact
	etotal = min( etotal, 1.0_r8 )

	rainsweepout = xnumrainsv(jr)*4.*pi*r*r*vfall

	scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
	scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

5500	continue

	scavsumnum = scavsumnum + scavsumnumbb
	scavsumvol = scavsumvol + scavsumvolbb
5900	continue

	scavratenum = scavsumnum*3600.
	scavratevol = scavsumvol*3600.


   
  return
  end subroutine calc_1_impact_rate




end module mz_aerosols_intr
