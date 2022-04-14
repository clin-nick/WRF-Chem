















   module modal_aero_calcsize




   implicit none
   private
   save
                                                                                                                             

   public modal_aero_calcsize_sub, modal_aero_calcsize_init
                                                                                                                             

   logical  :: do_adjust
   logical  :: do_aitacc_transfer

                                                                                                                             












                                                                                                                             
                                                                                                                             
  contains
                                                                                                                             
                                                                                                                             



    subroutine modal_aero_calcsize_sub(            &
         lchnk,   ncol,    nstep,   icalcaer_flag, &
         loffset,                                  &
         aero_mmr_flag,    deltat,                 &
         t,       pmid,    pdel,    q,             &
         dqdt,    dotend,                          &
         qqcw,                                     &
         dgncur_a       )

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

      use modal_aero_data,  only: numptrcw_amode, mprognum_amode, numptr_amode, &
      lmassptrcw_amode, &
           alnsg_amode, voltonumblo_amode, lmassptr_amode, modeptr_accum, ntot_amode, modeptr_aitken, ntot_aspectype, &
           dgnum_amode, lspectype_amode, specmw_amode, specdens_amode, voltonumb_amode, nspec_amode, dgnumlo_amode, &
           cnst_name_cw, voltonumbhi_amode, dgnumhi_amode
      use modal_aero_rename, only: lspectooa_renamexf, lspecfrma_renamexf, lspectooc_renamexf, lspecfrmc_renamexf, &
           modetoo_renamexf, nspecfrm_renamexf, npair_renamexf, modefrm_renamexf

      use shr_kind_mod,  only:  r8 => shr_kind_r8
      use physconst,     only:  avogad, gravit, mwdry, rair
      use physconst,     only: pi
      use module_cam_support, only: pcols, pver, pcnst =>pcnst_runtime, &
           endrun,fieldname_len, outfld, masterproc
      use constituents,  only: cnst_name

      implicit none

      
      
      integer, intent(in)  :: lchnk                
      integer, intent(in)  :: ncol                 

      integer, intent(in)  :: nstep                
      integer, intent(in)  :: icalcaer_flag
      
      integer, intent(in)  :: loffset
      

      logical, intent(in)  :: aero_mmr_flag        
      

      real(r8), intent(in) :: deltat               
      real(r8), intent(in) :: t(pcols,pver)        
      real(r8), intent(in) :: pmid(pcols,pver)     
      real(r8), intent(in) :: pdel(pcols,pver)     
      real(r8), intent(in) :: q(pcols,pver,pcnst)  

      logical,  intent(inout) :: dotend(pcnst)          
      real(r8), intent(inout) :: dqdt(pcols,pver,pcnst) 
      real(r8), target, intent(inout) :: qqcw(pcols,pver,pcnst) 
      
      
      
      
      
      real(r8), intent(inout) :: dgncur_a(pcols,pver,ntot_amode)

      
      
      integer  :: i, icol_diag, iduma, ipair, iq
      integer  :: ixfer_acc2ait, ixfer_ait2acc
      integer  :: ixfer_acc2ait_sv(pcols,pver), ixfer_ait2acc_sv(pcols,pver)
      integer  :: j, jac, jsrflx, k 
      integer  :: l, l1, la, lc, lna, lnc, lsfrm, lstoo
      integer  :: n, nacc, nait
      integer  :: lat(pcols), lon(pcols)

      integer, save  :: idiagaa = 1

      logical  :: dotendqqcw(pcnst)
      logical  :: noxf_acc2ait(ntot_aspectype)

      character(len=fieldname_len)   :: tmpnamea, tmpnameb
      character(len=fieldname_len+3) :: fieldname

      real(r8), parameter :: third = 1.0_r8/3.0_r8
      real(r8), pointer :: fldcw(:,:)
      real(r8) :: delnum_a2, delnum_c2            
      real(r8) :: delnum_a3, delnum_c3, delnum_t3 
      real(r8) :: deltatinv                     
      real(r8) :: dgncur_c(pcols,pver,ntot_amode)
      real(r8) :: dgnyy, dgnxx                  
      real(r8) :: dqqcwdt(pcols,pver,pcnst)     
      real(r8) :: drv_a, drv_c, drv_t           
      real(r8) :: drv_t0
      real(r8) :: drv_a_noxf, drv_c_noxf, drv_t_noxf 
      real(r8) :: drv_a_acc, drv_c_acc
      real(r8) :: drv_a_accsv(pcols,pver), drv_c_accsv(pcols,pver)
      real(r8) :: drv_a_aitsv(pcols,pver), drv_c_aitsv(pcols,pver)
      real(r8) :: drv_a_sv(pcols,pver,ntot_amode), drv_c_sv(pcols,pver,ntot_amode)
      real(r8) :: dryvol_a(pcols,pver)          
      
      real(r8) :: dryvol_c(pcols,pver)          
      real(r8) :: duma, dumb, dumc, dumd        
      real(r8) :: dumfac, dummwdens             
      real(r8) :: frelaxadj                     
      
      real(r8) :: fracadj                       
      real(r8) :: num_a0, num_c0, num_t0        
      real(r8) :: num_a1, num_c1                
      real(r8) :: num_a2, num_c2, num_t2        
      real(r8) :: num_a, num_c, num_t           
      real(r8) :: num_t_noxf
      real(r8) :: numbnd                        
      real(r8) :: num_a_acc, num_c_acc
      real(r8) :: num_a_accsv(pcols,pver), num_c_accsv(pcols,pver)
      real(r8) :: num_a_aitsv(pcols,pver), num_c_aitsv(pcols,pver)
      real(r8) :: num_a_sv(pcols,pver,ntot_amode), num_c_sv(pcols,pver,ntot_amode)
      real(r8) :: pdel_fac                      
      real(r8) :: tadj                          
      real(r8) :: tadjinv                       
      real(r8) :: v2ncur_a(pcols,pver,ntot_amode)
      real(r8) :: v2ncur_c(pcols,pver,ntot_amode)
      real(r8) :: v2nyy, v2nxx, v2nzz           
      real(r8) :: v2nyyrl, v2nxxrl              
      real(r8) :: xfercoef
      real(r8) :: xfercoef_num_acc2ait, xfercoef_vol_acc2ait
      real(r8) :: xfercoef_num_ait2acc, xfercoef_vol_ait2acc
      real(r8) :: xferfrac_num_acc2ait, xferfrac_vol_acc2ait
      real(r8) :: xferfrac_num_ait2acc, xferfrac_vol_ait2acc
      real(r8) :: xfertend, xfertend_num(2,2)

      integer, parameter :: nsrflx = 4    
      real(r8) :: qsrflx(pcols,pcnst,nsrflx,2)
      
      
      
      
      
      
      
      

      dotendqqcw(:) = .false.
      dqqcwdt(:,:,:) = 0.0_r8
      qsrflx(:,:,:,:) = 0.0_r8

      nait = modeptr_aitken
      nacc = modeptr_accum

      deltatinv = 1.0/(deltat*(1.0d0 + 1.0d-15))
      
      
      tadj = deltat
      tadj = 86400
      tadj = max( tadj, deltat )
      tadjinv = 1.0/(tadj*(1.0d0 + 1.0d-15))
      fracadj = deltat*tadjinv
      fracadj = max( 0.0_r8, min( 1.0_r8, fracadj ) )


      
      
      
      
      
      
      
      
      do n = 1, ntot_amode


         
         do k=1,pver
            do i=1,ncol
               
               
               dgncur_a(i,k,n) = dgnum_amode(n)
               dgncur_c(i,k,n) = dgnum_amode(n)
               v2ncur_a(i,k,n) = voltonumb_amode(n)
               v2ncur_c(i,k,n) = voltonumb_amode(n)
               dryvol_a(i,k) = 0.0_r8
               dryvol_c(i,k) = 0.0_r8
            end do
         end do

         
         
         do l1 = 1, nspec_amode(n)
            if ( aero_mmr_flag ) then
               
               dummwdens = 1.0_r8 / specdens_amode(lspectype_amode(l1,n))
            else
               
               dummwdens = specmw_amode(lspectype_amode(l1,n))   &
                    / specdens_amode(lspectype_amode(l1,n))
            end if
            la = lmassptr_amode(l1,n) - loffset
            do k=1,pver
               do i=1,ncol
                  dryvol_a(i,k) = dryvol_a(i,k)    &
                       + max(0.0_r8,q(i,k,la))*dummwdens
               end do
            end do
            
            fldcw => qqcw(:,:,lmassptrcw_amode(l1,n))
            do k=1,pver
               do i=1,ncol
                  dryvol_c(i,k) = dryvol_c(i,k)    &
                       + max(0.0_r8,fldcw(i,k))*dummwdens
               end do
            end do
         end do

         
         lna = numptr_amode(n) - loffset
         lnc = numptrcw_amode(n) - loffset
         fldcw => qqcw(:,:,numptrcw_amode(n))


         
         if (mprognum_amode(n) <= 0) then

            
            
            
            
            if (lna > 0) then
               dotend(lna) = .true.
               do k=1,pver
                  do i=1,ncol
                     dqdt(i,k,lna) = (dryvol_a(i,k)*voltonumb_amode(n)   &
                          - q(i,k,lna)) * deltatinv
                  end do
               end do
            end if
            if (lnc > 0) then
               dotendqqcw(lnc) = .true.
               do k=1,pver
                  do i=1,ncol
                     dqqcwdt(i,k,lnc) = (dryvol_c(i,k)*voltonumb_amode(n)   &
                          - fldcw(i,k)) * deltatinv
                  end do
               end do
            end if
         else


            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
         end if
         frelaxadj = 27.0
         dumfac = exp(4.5*alnsg_amode(n)**2)*pi/6.0
         v2nxx = voltonumbhi_amode(n)
         v2nyy = voltonumblo_amode(n)
         v2nxxrl = v2nxx/frelaxadj
         v2nyyrl = v2nyy*frelaxadj
         dgnxx = dgnumhi_amode(n)
         dgnyy = dgnumlo_amode(n)
         if ( do_aitacc_transfer ) then
            if (n == nait) v2nxx = v2nxx/1.0e6
            if (n == nacc) v2nyy = v2nyy*1.0e6
            v2nxxrl = v2nxx/frelaxadj   
            v2nyyrl = v2nyy*frelaxadj   
         end if

         if (do_adjust) then
            dotend(lna) = .true.
            dotendqqcw(lnc) = .true.
         end if

         do  k = 1, pver
            do  i = 1, ncol

               drv_a = dryvol_a(i,k)
               num_a0 = q(i,k,lna)
               num_a = max( 0.0_r8, num_a0 )
               drv_c = dryvol_c(i,k)
               num_c0 = fldcw(i,k)
               num_c = max( 0.0_r8, num_c0 )

               if ( do_adjust) then

                  
                  
                  
                  
                  
                  
                  
                  if ((drv_a <= 0.0_r8) .and. (drv_c <= 0.0_r8)) then
                     
                     
                     num_a = 0.0_r8
                     dqdt(i,k,lna) = -num_a0*deltatinv
                     num_c = 0.0_r8
                     dqqcwdt(i,k,lnc) = -num_c0*deltatinv
                  else if (drv_c <= 0.0_r8) then
                     
                     
                     num_c = 0.0_r8
                     dqqcwdt(i,k,lnc) = -num_c0*deltatinv
                     num_a1 = num_a
                     numbnd = max( drv_a*v2nxx, min( drv_a*v2nyy, num_a1 ) )
                     num_a  = num_a1 + (numbnd - num_a1)*fracadj
                     dqdt(i,k,lna) = (num_a - num_a0)*deltatinv

                  else if (drv_a <= 0.0_r8) then
                     
                     num_a = 0.0_r8
                     dqdt(i,k,lna) = -num_a0*deltatinv
                     num_c1 = num_c
                     numbnd = max( drv_c*v2nxx, min( drv_c*v2nyy, num_c1 ) )
                     num_c  = num_c1 + (numbnd - num_c1)*fracadj
                     dqqcwdt(i,k,lnc) = (num_c - num_c0)*deltatinv
                  else
                     
                     
                     
                     num_a1 = num_a
                     num_c1 = num_c
                     
                     
                     
                     
                     numbnd = max( drv_a*v2nxxrl, min( drv_a*v2nyyrl, num_a1 ) )
                     delnum_a2 = (numbnd - num_a1)*fracadj
                     num_a2 = num_a1 + delnum_a2
                     numbnd = max( drv_c*v2nxxrl, min( drv_c*v2nyyrl, num_c1 ) )
                     delnum_c2 = (numbnd - num_c1)*fracadj
                     num_c2 = num_c1 + delnum_c2
                     if ((delnum_a2 == 0.0) .and. (delnum_c2 /= 0.0)) then
                        num_a2 = max( drv_a*v2nxxrl, min( drv_a*v2nyyrl,   &
                             num_a1-delnum_c2 ) )
                     else if ((delnum_a2 /= 0.0) .and. (delnum_c2 == 0.0)) then
                        num_c2 = max( drv_c*v2nxxrl, min( drv_c*v2nyyrl,   &
                             num_c1-delnum_a2 ) )
                     end if
                     
                     
                     drv_t = drv_a + drv_c
                     num_t2 = num_a2 + num_c2
                     delnum_a3 = 0.0_r8
                     delnum_c3 = 0.0_r8
                     if (num_t2 < drv_t*v2nxx) then
                        delnum_t3 = (drv_t*v2nxx - num_t2)*fracadj
                        
                        
                        if ((num_a2 < drv_a*v2nxx) .and. (num_c2 < drv_c*v2nxx)) then
                           delnum_a3 = delnum_t3*(num_a2/num_t2)
                           delnum_c3 = delnum_t3*(num_c2/num_t2)
                        else if (num_c2 < drv_c*v2nxx) then
                           delnum_c3 = delnum_t3
                        else if (num_a2 < drv_a*v2nxx) then
                           delnum_a3 = delnum_t3
                        end if
                     else if (num_t2 > drv_t*v2nyy) then
                        delnum_t3 = (drv_t*v2nyy - num_t2)*fracadj
                        
                        
                        if ((num_a2 > drv_a*v2nyy) .and. (num_c2 > drv_c*v2nyy)) then
                           delnum_a3 = delnum_t3*(num_a2/num_t2)
                           delnum_c3 = delnum_t3*(num_c2/num_t2)
                        else if (num_c2 > drv_c*v2nyy) then
                           delnum_c3 = delnum_t3
                        else if (num_a2 > drv_a*v2nyy) then
                           delnum_a3 = delnum_t3
                        end if
                     end if
                     num_a = num_a2 + delnum_a3
                     dqdt(i,k,lna) = (num_a - num_a0)*deltatinv
                     num_c = num_c2 + delnum_c3
                     dqqcwdt(i,k,lnc) = (num_c - num_c0)*deltatinv
                  end if

               end if 

               
               
               
               if (drv_a > 0.0_r8) then
                  if (num_a <= drv_a*v2nxx) then
                     dgncur_a(i,k,n) = dgnxx
                     v2ncur_a(i,k,n) = v2nxx
                  else if (num_a >= drv_a*v2nyy) then
                     dgncur_a(i,k,n) = dgnyy
                     v2ncur_a(i,k,n) = v2nyy
                  else
                     dgncur_a(i,k,n) = (drv_a/(dumfac*num_a))**third
                     v2ncur_a(i,k,n) = num_a/drv_a
                  end if
               end if
               pdel_fac = pdel(i,k)/gravit   
               jac = 1
               qsrflx(i,lna,1,jac) = qsrflx(i,lna,1,jac) + max(0.0_r8,dqdt(i,k,lna))*pdel_fac
               qsrflx(i,lna,2,jac) = qsrflx(i,lna,2,jac) + min(0.0_r8,dqdt(i,k,lna))*pdel_fac

               if (drv_c > 0.0_r8) then
                  if (num_c <= drv_c*v2nxx) then
                     dgncur_c(i,k,n) = dgnumhi_amode(n)
                     v2ncur_c(i,k,n) = v2nxx
                  else if (num_c >= drv_c*v2nyy) then
                     dgncur_c(i,k,n) = dgnumlo_amode(n)
                     v2ncur_c(i,k,n) = v2nyy
                  else
                     dgncur_c(i,k,n) = (drv_c/(dumfac*num_c))**third
                     v2ncur_c(i,k,n) = num_c/drv_c
                  end if
               end if
               jac = 2
               qsrflx(i,lnc,1,jac) = qsrflx(i,lnc,1,jac) + max(0.0_r8,dqqcwdt(i,k,lnc))*pdel_fac
               qsrflx(i,lnc,2,jac) = qsrflx(i,lnc,2,jac) + min(0.0_r8,dqqcwdt(i,k,lnc))*pdel_fac


               
               if ( do_aitacc_transfer ) then
                  if (n == nait) then
                     drv_a_aitsv(i,k) = drv_a
                     num_a_aitsv(i,k) = num_a
                     drv_c_aitsv(i,k) = drv_c
                     num_c_aitsv(i,k) = num_c
                  else if (n == nacc) then
                     drv_a_accsv(i,k) = drv_a
                     num_a_accsv(i,k) = num_a
                     drv_c_accsv(i,k) = drv_c
                     num_c_accsv(i,k) = num_c
                  end if
               end if
               drv_a_sv(i,k,n) = drv_a
               num_a_sv(i,k,n) = num_a
               drv_c_sv(i,k,n) = drv_c
               num_c_sv(i,k,n) = num_c

            end do
         end do


         
         
         
         
      end do


      


      
      
      
      
      
      
      
      
      
      
      
      
      
      ixfer_ait2acc_sv(:,:) = 0
      ixfer_acc2ait_sv(:,:) = 0
      if ( do_aitacc_transfer ) then

         
         
         
         if (npair_renamexf .le. 0) then
            npair_renamexf = 0
            
            if (npair_renamexf .le. 0) then
               write( 6, '(//a//)' )   &
                    '*** modal_aero_calcaersize_sub error -- npair_renamexf <= 0'
               call endrun( 'modal_aero_calcaersize_sub error' )
            end if
         end if

         
         ipair = 1
         if ((modefrm_renamexf(ipair) .ne. nait) .or.   &
              (modetoo_renamexf(ipair) .ne. nacc)) then
            write( 6, '(//2a//)' )   &
                 '*** modal_aero_calcaersize_sub error -- ',   &
                 'modefrm/too_renamexf(1) are wrong'
            call endrun( 'modal_aero_calcaersize_sub error' )
         end if

         
         do iq = 1, nspecfrm_renamexf(ipair)
            lsfrm = lspecfrma_renamexf(iq,ipair) - loffset
            lstoo = lspectooa_renamexf(iq,ipair) - loffset
            if ((lsfrm > 0) .and. (lstoo > 0)) then
               dotend(lsfrm) = .true.
               dotend(lstoo) = .true.
            end if
            lsfrm = lspecfrmc_renamexf(iq,ipair) - loffset
            lstoo = lspectooc_renamexf(iq,ipair) - loffset
            if ((lsfrm > 0) .and. (lstoo > 0)) then
               dotendqqcw(lsfrm) = .true.
               dotendqqcw(lstoo) = .true.
            end if
         end do

         
         noxf_acc2ait(:) = .true.
         do l1 = 1, nspec_amode(nacc)
            la = lmassptr_amode(l1,nacc)
            do iq = 1, nspecfrm_renamexf(ipair)
               if (lspectooa_renamexf(iq,ipair) == la) then
                  noxf_acc2ait(l1) = .false.
               end if
            end do
         end do

         
         
         v2nzz = sqrt(voltonumb_amode(nait)*voltonumb_amode(nacc))

         
         do  k = 1, pver
            do  i = 1, ncol

               pdel_fac = pdel(i,k)/gravit   
               xfertend_num(:,:) = 0.0_r8

               
               ixfer_ait2acc = 0
               xfercoef_num_ait2acc = 0.0_r8
               xfercoef_vol_ait2acc = 0.0_r8

               drv_t = drv_a_aitsv(i,k) + drv_c_aitsv(i,k)
               num_t = num_a_aitsv(i,k) + num_c_aitsv(i,k)
               if (drv_t > 0.0_r8) then
                  if (num_t < drv_t*v2nzz) then
                     ixfer_ait2acc = 1
                     if (num_t < drv_t*voltonumb_amode(nacc)) then
                        xferfrac_num_ait2acc = 1.0_r8
                        xferfrac_vol_ait2acc = 1.0_r8
                     else
                        xferfrac_vol_ait2acc = ((num_t/drv_t) - v2nzz)/   &
                             (voltonumb_amode(nacc) - v2nzz)
                        xferfrac_num_ait2acc = xferfrac_vol_ait2acc*   &
                             (drv_t*voltonumb_amode(nacc)/num_t)
                        if ((xferfrac_num_ait2acc <= 0.0_r8) .or.   &
                             (xferfrac_vol_ait2acc <= 0.0_r8)) then
                           xferfrac_num_ait2acc = 0.0_r8
                           xferfrac_vol_ait2acc = 0.0_r8
                        else if ((xferfrac_num_ait2acc >= 1.0_r8) .or.   &
                             (xferfrac_vol_ait2acc >= 1.0_r8)) then
                           xferfrac_num_ait2acc = 1.0_r8
                           xferfrac_vol_ait2acc = 1.0_r8
                        end if
                     end if
                     xfercoef_num_ait2acc = xferfrac_num_ait2acc*tadjinv
                     xfercoef_vol_ait2acc = xferfrac_vol_ait2acc*tadjinv
                     xfertend_num(1,1) = num_a_aitsv(i,k)*xfercoef_num_ait2acc
                     xfertend_num(1,2) = num_c_aitsv(i,k)*xfercoef_num_ait2acc
                  end if
               end if

               
               
               
               
               
               
               ixfer_acc2ait = 0
               xfercoef_num_acc2ait = 0.0_r8
               xfercoef_vol_acc2ait = 0.0_r8

               drv_t = drv_a_accsv(i,k) + drv_c_accsv(i,k)
               num_t = num_a_accsv(i,k) + num_c_accsv(i,k)
               drv_a_noxf = 0.0_r8
               drv_c_noxf = 0.0_r8
               if (drv_t > 0.0_r8) then
                  if (num_t > drv_t*v2nzz) then
                     do l1 = 1, nspec_amode(nacc)

                        if ( noxf_acc2ait(l1) ) then
                           if ( aero_mmr_flag ) then
                              
                              dummwdens = 1.0_r8   &
                                   / specdens_amode(lspectype_amode(l1,nacc))
                           else
                              
                              dummwdens = specmw_amode(lspectype_amode(l1,nacc))   &
                                   / specdens_amode(lspectype_amode(l1,nacc))
                           end if
                           la = lmassptr_amode(l1,nacc) - loffset
                           drv_a_noxf = drv_a_noxf    &
                                + max(0.0_r8,q(i,k,la))*dummwdens
                           lc = lmassptrcw_amode(l1,nacc) - loffset
                            fldcw => qqcw(:,:,lmassptrcw_amode(l1,nacc))
                           drv_c_noxf = drv_c_noxf    &
                                + max(0.0_r8,fldcw(i,k))*dummwdens
                        end if
                     end do
                     drv_t_noxf = drv_a_noxf + drv_c_noxf
                     num_t_noxf = drv_t_noxf*voltonumblo_amode(nacc)
                     num_t0 = num_t
                     drv_t0 = drv_t
                     num_t = max( 0.0_r8, num_t - num_t_noxf )
                     drv_t = max( 0.0_r8, drv_t - drv_t_noxf )
                  end if
               end if

               if (drv_t > 0.0_r8) then
                  if (num_t > drv_t*v2nzz) then
                     ixfer_acc2ait = 1
                     if (num_t > drv_t*voltonumb_amode(nait)) then
                        xferfrac_num_acc2ait = 1.0
                        xferfrac_vol_acc2ait = 1.0
                     else
                        xferfrac_vol_acc2ait = ((num_t/drv_t) - v2nzz)/   &
                             (voltonumb_amode(nait) - v2nzz)
                        xferfrac_num_acc2ait = xferfrac_vol_acc2ait*   &
                             (drv_t*voltonumb_amode(nait)/num_t)
                        if ((xferfrac_num_acc2ait <= 0.0_r8) .or.   &
                             (xferfrac_vol_acc2ait <= 0.0_r8)) then
                           xferfrac_num_acc2ait = 0.0_r8
                           xferfrac_vol_acc2ait = 0.0_r8
                        else if ((xferfrac_num_acc2ait >= 1.0_r8) .or.   &
                             (xferfrac_vol_acc2ait >= 1.0_r8)) then
                           xferfrac_num_acc2ait = 1.0_r8
                           xferfrac_vol_acc2ait = 1.0_r8
                        end if
                     end if
                     duma = 1.0e-37
                     xferfrac_num_acc2ait = xferfrac_num_acc2ait*   &
                          num_t/max( duma, num_t0 )
                     xfercoef_num_acc2ait = xferfrac_num_acc2ait*tadjinv
                     xfercoef_vol_acc2ait = xferfrac_vol_acc2ait*tadjinv
                     xfertend_num(2,1) = num_a_accsv(i,k)*xfercoef_num_acc2ait
                     xfertend_num(2,2) = num_c_accsv(i,k)*xfercoef_num_acc2ait
                  end if
               end if

               
               if (ixfer_ait2acc+ixfer_acc2ait > 0) then
                  ixfer_ait2acc_sv(i,k) = ixfer_ait2acc
                  ixfer_acc2ait_sv(i,k) = ixfer_acc2ait

                  
                  
                  
                  
                  do n = nait, nacc, (nacc-nait)
                     if (n .eq. nait) then
                        duma = (xfertend_num(1,1) - xfertend_num(2,1))*deltat
                        num_a     = max( 0.0_r8, num_a_aitsv(i,k) - duma )
                        num_a_acc = max( 0.0_r8, num_a_accsv(i,k) + duma )
                        duma = (drv_a_aitsv(i,k)*xfercoef_vol_ait2acc -   &
                             (drv_a_accsv(i,k)-drv_a_noxf)*xfercoef_vol_acc2ait)*deltat
                        drv_a     = max( 0.0_r8, drv_a_aitsv(i,k) - duma )
                        drv_a_acc = max( 0.0_r8, drv_a_accsv(i,k) + duma )
                        duma = (xfertend_num(1,2) - xfertend_num(2,2))*deltat
                        num_c     = max( 0.0_r8, num_c_aitsv(i,k) - duma )
                        num_c_acc = max( 0.0_r8, num_c_accsv(i,k) + duma )
                        duma = (drv_c_aitsv(i,k)*xfercoef_vol_ait2acc -   &
                             (drv_c_accsv(i,k)-drv_c_noxf)*xfercoef_vol_acc2ait)*deltat
                        drv_c     = max( 0.0_r8, drv_c_aitsv(i,k) - duma )
                        drv_c_acc = max( 0.0_r8, drv_c_accsv(i,k) + duma )
                     else
                        num_a = num_a_acc
                        drv_a = drv_a_acc
                        num_c = num_c_acc
                        drv_c = drv_c_acc
                     end if

                     if (drv_a > 0.0_r8) then
                        if (num_a <= drv_a*voltonumbhi_amode(n)) then
                           dgncur_a(i,k,n) = dgnumhi_amode(n)
                           v2ncur_a(i,k,n) = voltonumbhi_amode(n)
                        else if (num_a >= drv_a*voltonumblo_amode(n)) then
                           dgncur_a(i,k,n) = dgnumlo_amode(n)
                           v2ncur_a(i,k,n) = voltonumblo_amode(n)
                        else
                           dgncur_a(i,k,n) = (drv_a/(dumfac*num_a))**third
                           v2ncur_a(i,k,n) = num_a/drv_a
                        end if
                     else
                        dgncur_a(i,k,n) = dgnum_amode(n)
                        v2ncur_a(i,k,n) = voltonumb_amode(n)
                     end if

                     if (drv_c > 0.0_r8) then
                        if (num_c <= drv_c*voltonumbhi_amode(n)) then
                           dgncur_c(i,k,n) = dgnumhi_amode(n)
                           v2ncur_c(i,k,n) = voltonumbhi_amode(n)
                        else if (num_c >= drv_c*voltonumblo_amode(n)) then
                           dgncur_c(i,k,n) = dgnumlo_amode(n)
                           v2ncur_c(i,k,n) = voltonumblo_amode(n)
                        else
                           dgncur_c(i,k,n) = (drv_c/(dumfac*num_c))**third
                           v2ncur_c(i,k,n) = num_c/drv_c
                        end if
                     else
                        dgncur_c(i,k,n) = dgnum_amode(n)
                        v2ncur_c(i,k,n) = voltonumb_amode(n)
                     end if

                  end do


                  
                  
                  

                  if ( masterproc ) then
                     if (idiagaa > 0) then
                        do j = 1, 2
                           do iq = 1, nspecfrm_renamexf(ipair)
                              do jac = 1, 2
                                 if (j .eq. 1) then
                                    if (jac .eq. 1) then
                                       lsfrm = lspecfrma_renamexf(iq,ipair)
                                       lstoo = lspectooa_renamexf(iq,ipair)
                                    else
                                       lsfrm = lspecfrmc_renamexf(iq,ipair)
                                       lstoo = lspectooc_renamexf(iq,ipair)
                                    end if
                                 else
                                    if (jac .eq. 1) then
                                       lsfrm = lspectooa_renamexf(iq,ipair)
                                       lstoo = lspecfrma_renamexf(iq,ipair)
                                    else
                                       lsfrm = lspectooc_renamexf(iq,ipair)
                                       lstoo = lspecfrmc_renamexf(iq,ipair)
                                    end if
                                 end if
                                 write( 6, '(a,3i3,2i4)' ) 'calcsize j,iq,jac, lsfrm,lstoo',   &
                                      j,iq,jac, lsfrm,lstoo
                              end do
                           end do
                        end do
                     end if
                  end if
                  idiagaa = -1


                  
                  do  j = 1, 2

                     if ((j .eq. 1 .and. ixfer_ait2acc > 0) .or. &
                          (j .eq. 2 .and. ixfer_acc2ait > 0)) then

                        jsrflx = j+2
                        if (j .eq. 1) then
                           xfercoef = xfercoef_vol_ait2acc
                        else
                           xfercoef = xfercoef_vol_acc2ait
                        end if

                        do  iq = 1, nspecfrm_renamexf(ipair)

                           
                           do  jac = 1, 2

                              
                              
                              
                              
                              if (j .eq. 1) then
                                 if (jac .eq. 1) then
                                    lsfrm = lspecfrma_renamexf(iq,ipair)
                                    lstoo = lspectooa_renamexf(iq,ipair)
                                 else
                                    lsfrm = lspecfrmc_renamexf(iq,ipair)
                                    lstoo = lspectooc_renamexf(iq,ipair)
                                 end if
                              else
                                 if (jac .eq. 1) then
                                    lsfrm = lspectooa_renamexf(iq,ipair)
                                    lstoo = lspecfrma_renamexf(iq,ipair)
                                 else
                                    lsfrm = lspectooc_renamexf(iq,ipair)
                                    lstoo = lspecfrmc_renamexf(iq,ipair)
                                 end if
                              end if

                              lsfrm = lsfrm - loffset

                              lstoo = lstoo - loffset
                              if ((lsfrm > 0) .and. (lstoo > 0)) then
                                 if (jac .eq. 1) then
                                    if (iq .eq. 1) then
                                       xfertend = xfertend_num(j,jac)
                                    else
                                       xfertend = max(0.0_r8,q(i,k,lsfrm))*xfercoef
                                    end if
                                    dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xfertend
                                    dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xfertend
                                 else
                                    if (iq .eq. 1) then
                                       xfertend = xfertend_num(j,jac)
                                    else
                                       fldcw => qqcw(:,:,lsfrm+loffset)
                                       xfertend = max(0.0_r8,fldcw(i,k))*xfercoef
                                    end if
                                    dqqcwdt(i,k,lsfrm) = dqqcwdt(i,k,lsfrm) - xfertend
                                    dqqcwdt(i,k,lstoo) = dqqcwdt(i,k,lstoo) + xfertend
                                 end if
                                 qsrflx(i,lsfrm,jsrflx,jac) = qsrflx(i,lsfrm,jsrflx,jac) - xfertend*pdel_fac
                                 qsrflx(i,lstoo,jsrflx,jac) = qsrflx(i,lstoo,jsrflx,jac) + xfertend*pdel_fac
                              end if

                           end do
                        end do
                     end if
                  end do

               end if
            end do
         end do


      end if  
      lsfrm = -123456789   

      

      
      
      
      do l = 1, pcnst
         lc = l - loffset
         if ( lc>0 .and. dotendqqcw(lc) ) then
            fldcw=> qqcw(:,:,l)
            do k = 1, pver
               do i = 1, ncol
                  fldcw(i,k) = max( 0.0_r8,   &
                       (fldcw(i,k) + dqqcwdt(i,k,lc)*deltat) )
               end do
            end do
         end if
      end do

      
      
      

      
      if ( .not. do_adjust ) return

      do n = 1, ntot_amode 
         if (mprognum_amode(n) <= 0) cycle

         do jac = 1, 2
            if (jac == 1) then
               l = numptr_amode(n)
               tmpnamea = cnst_name(l)
            else
               l = numptrcw_amode(n)
               tmpnamea = cnst_name_cw(l)
            end if
            l = l - loffset
            fieldname = trim(tmpnamea) // '_sfcsiz1'
            call outfld( fieldname, qsrflx(:,l,1,jac), pcols, lchnk)

            fieldname = trim(tmpnamea) // '_sfcsiz2'
            call outfld( fieldname, qsrflx(:,l,2,jac), pcols, lchnk)
         end do   

      end do   


      
      if ( .not. do_aitacc_transfer ) return

      do iq = 1, nspecfrm_renamexf(ipair)

         
         do jac = 1, 2

            
            
            if (jac .eq. 1) then
               lsfrm = lspecfrma_renamexf(iq,ipair)
               lstoo = lspectooa_renamexf(iq,ipair)
            else
               lsfrm = lspecfrmc_renamexf(iq,ipair)
               lstoo = lspectooc_renamexf(iq,ipair)
            end if
            if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

            if (jac .eq. 1) then
               tmpnamea = cnst_name(lsfrm)
               tmpnameb = cnst_name(lstoo)
            else
               tmpnamea = cnst_name_cw(lsfrm)
               tmpnameb = cnst_name_cw(lstoo)
            end if
            lsfrm = lsfrm - loffset
            lstoo = lstoo - loffset
            if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

            fieldname = trim(tmpnamea) // '_sfcsiz3'
            call outfld( fieldname, qsrflx(:,lsfrm,3,jac), pcols, lchnk)

            fieldname = trim(tmpnameb) // '_sfcsiz3'
            call outfld( fieldname, qsrflx(:,lstoo,3,jac), pcols, lchnk)

            fieldname = trim(tmpnamea) // '_sfcsiz4'
            call outfld( fieldname, qsrflx(:,lsfrm,4,jac), pcols, lchnk)

            fieldname = trim(tmpnameb) // '_sfcsiz4'
            call outfld( fieldname, qsrflx(:,lstoo,4,jac), pcols, lchnk)

         end do   
      end do   
    end subroutine modal_aero_calcsize_sub
 




subroutine modal_aero_calcsize_init












use modal_aero_data
use modal_aero_rename
use module_cam_support, only: endrun, addfld, add_default, fieldname_len, phys_decomp, &
     pcnst =>pcnst_runtime, masterproc
use constituents,  only :  cnst_name


implicit none






   integer  :: ipair, iq
   integer  :: jac
   integer  :: lsfrm, lstoo
   integer  :: n, nacc, nait

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit

   logical                        :: history_aerosol      
 
   
      history_aerosol = .FALSE.



      do_adjust = .true.




      nait = modeptr_aitken
      nacc = modeptr_accum
      do_aitacc_transfer = .false.
      if ((modeptr_aitken > 0) .and.   &
          (modeptr_accum  > 0) .and.   &
          (modeptr_aitken /= modeptr_accum)) then
         do_aitacc_transfer = .true.
         if (mprognum_amode(nait) <= 0) do_aitacc_transfer = .false.
         if (mprognum_amode(nacc) <= 0) do_aitacc_transfer = .false.
      end if




      if ( .not. do_adjust ) return
      do n = 1, ntot_amode 
         if (mprognum_amode(n) <= 0) cycle

         do jac = 1, 2
            if (jac == 1) then
               tmpnamea = cnst_name(numptr_amode(n))
            else
               tmpnamea = cnst_name_cw(numptrcw_amode(n))
            end if
            unit = '#/m2/s'
            fieldname = trim(tmpnamea) // '_sfcsiz1'
            long_name = trim(tmpnamea) // ' calcsize number-adjust column source'
            call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
            endif

            if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

            fieldname = trim(tmpnamea) // '_sfcsiz2'
            long_name = trim(tmpnamea) // ' calcsize number-adjust column sink'
            call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
            endif
            if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname
         end do   

      end do   



      if ( .not. do_aitacc_transfer ) return


      ipair = 1
      if ((modefrm_renamexf(ipair) .ne. nait) .or.   &
          (modetoo_renamexf(ipair) .ne. nacc)) then
         write( 6, '(//2a//)' )   &
            '*** modal_aero_calcaersize_init error -- ',   &
            'modefrm/too_renamexf(1) are wrong'
         call endrun( 'modal_aero_calcaersize_init error' )
      end if

      do iq = 1, nspecfrm_renamexf(ipair)


         do jac = 1, 2



            if (jac .eq. 1) then
               lsfrm = lspecfrma_renamexf(iq,ipair)
               lstoo = lspectooa_renamexf(iq,ipair)
            else
               lsfrm = lspecfrmc_renamexf(iq,ipair)
               lstoo = lspectooc_renamexf(iq,ipair)
            end if
            if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

            if (jac .eq. 1) then
               tmpnamea = cnst_name(lsfrm)
               tmpnameb = cnst_name(lstoo)
            else
               tmpnamea = cnst_name_cw(lsfrm)
               tmpnameb = cnst_name_cw(lstoo)
            end if

            unit = 'kg/m2/s'
            if ((tmpnamea(1:3) == 'num') .or. &
                (tmpnamea(1:3) == 'NUM')) unit = '#/m2/s'
            fieldname = trim(tmpnamea) // '_sfcsiz3'
            long_name = trim(tmpnamea) // ' calcsize aitken-to-accum adjust column tendency'
            call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
            endif
            if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

            fieldname = trim(tmpnameb) // '_sfcsiz3'
            long_name = trim(tmpnameb) // ' calcsize aitken-to-accum adjust column tendency'
            call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
            endif
            if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

            fieldname = trim(tmpnamea) // '_sfcsiz4'
            long_name = trim(tmpnamea) // ' calcsize accum-to-aitken adjust column tendency'
            call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
            endif
            if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

            fieldname = trim(tmpnameb) // '_sfcsiz4'
            long_name = trim(tmpnameb) // ' calcsize accum-to-aitken adjust column tendency'
            call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
            endif
            if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         end do   
      end do   


      return
      end subroutine modal_aero_calcsize_init




   end module modal_aero_calcsize



