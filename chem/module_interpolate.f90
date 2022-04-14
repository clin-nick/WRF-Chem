module module_interpolate












  implicit none

  private

  integer,parameter :: r8 = selected_real_kind(12) 




  public :: interp_type, lininterp, vertinterp, bilin, get_timeinterp_factors
  public :: lininterp_init, lininterp_finish
  type interp_type
     real(r8), pointer :: wgts(:)
     real(r8), pointer :: wgtn(:)
     integer, pointer  :: jjm(:)
     integer, pointer  :: jjp(:)
  end type interp_type
  interface lininterp
     module procedure lininterp_original
     module procedure lininterp_full1d
     module procedure lininterp1d
     module procedure lininterp2d2d
     module procedure lininterp2d1d
     module procedure lininterp3d2d
  end interface
  

contains
  subroutine lininterp_full1d (arrin, yin, nin, arrout, yout, nout)
    integer, intent(in) :: nin, nout
    real(r8), intent(in) :: arrin(nin), yin(nin), yout(nout)
    real(r8), intent(out) :: arrout(nout)
    type (interp_type) :: interp_wgts

    call lininterp_init(yin, nin, yout, nout, 1, interp_wgts)
    call lininterp1d(arrin, nin, arrout, nout, interp_wgts)
    call lininterp_finish(interp_wgts)

  end subroutine lininterp_full1d

  subroutine lininterp_init(yin, nin, yout, nout, extrap_method, interp_wgts, &
       cyclicmin, cyclicmax)
















    integer, intent(in) :: nin
    integer, intent(in) :: nout
    real(r8), intent(in) :: yin(:)           
    real(r8), intent(in) :: yout(:)         
    integer, intent(in) :: extrap_method       
                                               
                                               
    real(r8), intent(in), optional :: cyclicmin, cyclicmax

    type (interp_type), intent(out) :: interp_wgts

    real(r8) :: cmin, cmax
    real(r8) :: extrap
    real(r8) :: dyinwrap
    real(r8) :: ratio
    real(r8) :: avgdyin
    integer :: i, j, icount
    integer :: jj
    real(r8), pointer :: wgts(:)
    real(r8), pointer :: wgtn(:)
    integer, pointer :: jjm(:)
    integer, pointer :: jjp(:)
    logical :: increasing
    character(len=132) :: message
    
    
    
    
    if (nin.lt.2) then
       call wrf_abort  
    end if
    if(present(cyclicmin)) then
       cmin=cyclicmin
    else
       cmin=0_r8
    end if
    if(present(cyclicmax)) then
       cmax=cyclicmax
    else
       cmax=360_r8
    end if
    if(cmax<=cmin) then
       call wrf_abort  
    end if
    increasing=.true.
    icount = 0
    do j=1,nin-1
       if (yin(j).gt.yin(j+1)) icount = icount + 1
    end do
    if(icount.eq.nin-1) then
       increasing = .false.
       icount=0
    endif
    if (icount.gt.0) then
       call wrf_abort  
    end if
    allocate(interp_wgts%jjm(nout), &
         interp_wgts%jjp(nout), &
         interp_wgts%wgts(nout), &
         interp_wgts%wgtn(nout))

    jjm => interp_wgts%jjm
    jjp => interp_wgts%jjp
    wgts =>  interp_wgts%wgts
    wgtn =>  interp_wgts%wgtn

    
    
    
    jjm = 0
    jjp = 0

    extrap = 0.
    if(extrap_method.eq.0) then
       
       
       
       
       do j=1,nout
          if(increasing) then
             if (yout(j).lt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 0.
                wgtn(j) = 0.
                extrap = extrap + 1.
             else if (yout(j).gt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 0.
                wgtn(j) = 0.
                extrap = extrap + 1.
             end if
          else
             if (yout(j).gt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 0.
                wgtn(j) = 0.
                extrap = extrap + 1.
             else if (yout(j).lt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 0.
                wgtn(j) = 0.
                extrap = extrap + 1.
             end if
          end if
       end do
    else if(extrap_method.eq.1) then
       
       
       
       
       do j=1,nout
          if(increasing) then
             if (yout(j).le.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 1.
                wgtn(j) = 0.
                extrap = extrap + 1.
             else if (yout(j).gt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 1.
                wgtn(j) = 0.
                extrap = extrap + 1.
             end if
          else
             if (yout(j).gt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 1.
                wgtn(j) = 0.
                extrap = extrap + 1.
             else if (yout(j).le.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 1.
                wgtn(j) = 0.
                extrap = extrap + 1.
             end if
          end if
       end do
    else if(extrap_method.eq.2) then
       
       
       
       
       dyinwrap = yin(1) + (cmax-cmin) - yin(nin)
       avgdyin = abs(yin(nin)-yin(1))/(nin-1.)
       ratio = dyinwrap/avgdyin
       if (ratio < 0.9 .or. ratio > 1.1) then
          write(message,*) 'Lininterp: Bad dyinwrap value =',dyinwrap,&
               ' avg=', avgdyin, yin(1),yin(nin)
          call wrf_message( trim(message) )
          call wrf_abort  
       end if

       do j=1,nout
          if(increasing) then
             if (yout(j) <= yin(1)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
             else if (yout(j) > yin(nin)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)-yin(nin))/dyinwrap
             end if
          else
             if (yout(j) > yin(1)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
             else if (yout(j) <= yin(nin)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin)-yin(nin))/dyinwrap
             end if

          endif
       end do
    end if

    
    
    
    if(increasing) then
       do j=1,nout
          do jj=1,nin-1
             if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    else
       do j=1,nout
          do jj=1,nin-1
             if (yout(j).le.yin(jj) .and. yout(j).gt.yin(jj+1)) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    end if

    
    
    
    extrap = 100.*extrap/real(nout,r8)
    if (extrap.gt.50. ) then
       write(message,*) 'interpolate_data:','yout=',minval(yout),maxval(yout),increasing,nout
       call wrf_message( trim(message) )
       write(message,*) 'interpolate_data:','yin=',yin(1),yin(nin)
       call wrf_message( trim(message) )
       write(message,*) 'interpolate_data:',extrap,' % of output grid will have to be extrapolated'
       call wrf_message( trim(message) )
       call wrf_abort  
    end if

    
    
    
    icount = 0
    do j=1,nout
       if (jjm(j).eq.0 .or. jjp(j).eq.0) icount = icount + 1
       ratio=wgts(j)+wgtn(j)
       if((ratio<0.9.or.ratio>1.1).and.extrap_method.ne.0) then
          write(message,*) j, wgts(j),wgtn(j),jjm(j),jjp(j), increasing,extrap_method
          call wrf_message( trim(message) )
          call wrf_abort  
       end if
    end do
    if (icount.gt.0) then
       call wrf_abort  
    end if

  end subroutine lininterp_init

  subroutine lininterp1d (arrin, n1, arrout, m1, interp_wgts)
    
    
    
    
    
    
    
    
    
    
    implicit none
    
    
    
    
    integer, intent(in) :: n1                 
    integer, intent(in) :: m1                
                                                                                                  
    real(r8), intent(in) :: arrin(n1)    
    type(interp_type), intent(in) :: interp_wgts
    real(r8), intent(out) :: arrout(m1) 
                                                                                                  
    
    
    
    integer j                
    integer, pointer :: jjm(:)
    integer, pointer :: jjp(:)
                                                                                                  
    real(r8), pointer :: wgts(:)
    real(r8), pointer :: wgtn(:)
                                                                                                  
                                                                                                  
    jjm => interp_wgts%jjm
    jjp => interp_wgts%jjp
    wgts =>  interp_wgts%wgts
    wgtn =>  interp_wgts%wgtn
                                                                                                  
    
    
    
    do j=1,m1
      arrout(j) = arrin(jjm(j))*wgts(j) + arrin(jjp(j))*wgtn(j)
    end do
                                                                                                  
    return
  end subroutine lininterp1d
                                                                                                  
  subroutine lininterp2d2d(arrin, n1, n2, arrout, m1, m2, wgt1, wgt2)
    implicit none
    
    
    
    
    integer, intent(in) :: n1, n2, m1, m2
    real(r8), intent(in) :: arrin(n1,n2)    
    type(interp_type), intent(in) :: wgt1, wgt2
    real(r8), intent(out) :: arrout(m1,m2) 
    
    
    
    integer i,j                
    integer, pointer :: iim(:), jjm(:)
    integer, pointer :: iip(:), jjp(:)
                                                                                                  
    real(r8), pointer :: wgts1(:), wgts2(:)
    real(r8), pointer :: wgtn1(:), wgtn2(:)
                                                                                                  
    real(r8) :: arrtmp(n1,m2)
                                                                                                  
                                                                                                  
    jjm => wgt2%jjm
    jjp => wgt2%jjp
    wgts2 => wgt2%wgts
    wgtn2 => wgt2%wgtn
                                                                                                  
    iim => wgt1%jjm
    iip => wgt1%jjp
    wgts1 => wgt1%wgts
    wgtn1 => wgt1%wgtn

    do j=1,m2
      do i=1,n1
        arrtmp(i,j) = arrin(i,jjm(j))*wgts2(j) + arrin(i,jjp(j))*wgtn2(j)
      end do
    end do

    do j=1,m2
      do i=1,m1
        arrout(i,j) = arrtmp(iim(i),j)*wgts1(i) + arrtmp(iip(i),j)*wgtn1(i)
      end do
    end do
                                                                                                  
  end subroutine lininterp2d2d
  subroutine lininterp2d1d(arrin, n1, n2, arrout, m1, wgt1, wgt2, fldname)
    implicit none
    
    
    
    
    integer, intent(in) :: n1, n2, m1
    real(r8), intent(in) :: arrin(n1,n2)    
    type(interp_type), intent(in) :: wgt1, wgt2
    real(r8), intent(out) :: arrout(m1) 
    character(len=*), intent(in), optional :: fldname(:)
    
    
    
    integer i                
    integer, pointer :: iim(:), jjm(:)
    integer, pointer :: iip(:), jjp(:)
                                                                                                  
    real(r8), pointer :: wgts(:), wgte(:)
    real(r8), pointer :: wgtn(:), wgtw(:)

    jjm => wgt2%jjm
    jjp => wgt2%jjp
    wgts => wgt2%wgts
    wgtn => wgt2%wgtn
                                                                                                  
    iim => wgt1%jjm
    iip => wgt1%jjp
    wgtw => wgt1%wgts
    wgte => wgt1%wgtn

    do i=1,m1
       arrout(i) = arrin(iim(i),jjm(i))*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i))*wgte(i)*wgts(i) + &
                   arrin(iim(i),jjp(i))*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i))*wgte(i)*wgtn(i)
    end do

                                                                                                  
  end subroutine lininterp2d1d
  subroutine lininterp3d2d(arrin, n1, n2, n3, arrout, m1, len1, wgt1, wgt2)
    implicit none
    
    
    
    
    integer, intent(in) :: n1, n2, n3, m1, len1   
    real(r8), intent(in) :: arrin(n1,n2,n3)    
    type(interp_type), intent(in) :: wgt1, wgt2
    real(r8), intent(out) :: arrout(len1, n3) 

    
    
    
    integer i, k               
    integer, pointer :: iim(:), jjm(:)
    integer, pointer :: iip(:), jjp(:)
                                                                                                  
    real(r8), pointer :: wgts(:), wgte(:)
    real(r8), pointer :: wgtn(:), wgtw(:)

    jjm => wgt2%jjm
    jjp => wgt2%jjp
    wgts => wgt2%wgts
    wgtn => wgt2%wgtn
                                                                                                  
    iim => wgt1%jjm
    iip => wgt1%jjp
    wgtw => wgt1%wgts
    wgte => wgt1%wgtn

    do k=1,n3
       do i=1,m1
          arrout(i,k) = arrin(iim(i),jjm(i),k)*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i),k)*wgte(i)*wgts(i) + &
               arrin(iim(i),jjp(i),k)*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i),k)*wgte(i)*wgtn(i)
       end do
    end do
                                                                                                  
  end subroutine lininterp3d2d
                                                                                                  



  subroutine lininterp_finish(interp_wgts)
    type(interp_type) :: interp_wgts

    deallocate(interp_wgts%jjm, &
         interp_wgts%jjp, &
         interp_wgts%wgts, &
         interp_wgts%wgtn)

    nullify(interp_wgts%jjm, &
         interp_wgts%jjp, &
         interp_wgts%wgts, &
         interp_wgts%wgtn)
  end subroutine lininterp_finish

  subroutine lininterp_original (arrin, yin, nlev, nlatin, arrout, &
       yout, nlatout)
    
    
    
    
    
    
    
    
    
    
    
    
    
    implicit none
    
    
    
    
    integer, intent(in) :: nlev                   
    integer, intent(in) :: nlatin                 
    integer, intent(in) :: nlatout                

    real(r8), intent(in) :: arrin(nlev,nlatin)    
    real(r8), intent(in) :: yin(nlatin)           
    real(r8), intent(in) :: yout(nlatout)         

    real(r8), intent(out) :: arrout(nlev,nlatout) 
    
    
    
    integer j, jj              
    integer jjprev             
    integer k                  
    integer icount             

    real(r8) extrap            
    
    
    
    integer :: jjm(nlatout)
    integer :: jjp(nlatout)

    real(r8) :: wgts(nlatout)
    real(r8) :: wgtn(nlatout)
    
    
    
    
    if (nlatin.lt.2) then
       call wrf_abort  
    end if

    icount = 0
    do j=1,nlatin-1
       if (yin(j).gt.yin(j+1)) icount = icount + 1
    end do


    if (icount.gt.0) then
       call wrf_abort  
    end if
    
    
    
    do j=1,nlatout
       jjm(j) = 0
       jjp(j) = 0
    end do
    
    
    
    
    extrap = 0.

    do j=1,nlatout
       if (yout(j).le.yin(1)) then
          jjm(j) = 1
          jjp(j) = 1
          wgts(j) = 1.
          wgtn(j) = 0.
          extrap=extrap+1.
       else if (yout(j).gt.yin(nlatin)) then
          jjm(j) = nlatin
          jjp(j) = nlatin
          wgts(j) = 1.
          wgtn(j) = 0.
          extrap=extrap+1.
       endif
    end do

    
    
    
    do j=1,nlatout
       do jj=1,nlatin-1
          if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
             jjm(j) = jj
             jjp(j) = jj + 1
             wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
             wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
             exit
          end if
       end do
    end do
    
    
    
    icount = 0
    do j=1,nlatout
       if (jjm(j).eq.0 .or. jjp(j).eq.0) then
          icount = icount + 1
       end if
    end do
    if (icount.gt.0) then
       call wrf_abort  
    end if
    
    
    
    do j=1,nlatout
       do k=1,nlev
          arrout(k,j) = arrin(k,jjm(j))*wgts(j) + arrin(k,jjp(j))*wgtn(j)
       end do
    end do

    return
  end subroutine lininterp_original


  subroutine bilin (arrin, xin, yin, nlondin, nlonin, &
       nlevdin, nlev, nlatin, arrout, xout, &
       yout, nlondout, nlonout, nlevdout, nlatout)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    
    implicit none
    
    integer,parameter :: r8 = selected_real_kind(12) 
    
    
    
    integer, intent(in) :: nlondin                        
    integer, intent(in) :: nlonin                         
    integer, intent(in) :: nlevdin                        
    integer, intent(in) :: nlev                           
    integer, intent(in) :: nlatin                         
    integer, intent(in) :: nlatout                        
    integer, intent(in) :: nlondout                       
    integer, intent(in) :: nlonout(nlatout)               
    integer, intent(in) :: nlevdout                       

    real(r8), intent(in) :: arrin(nlondin,nlevdin,nlatin) 
    real(r8), intent(in) :: xin(nlondin)                  
    real(r8), intent(in) :: yin(nlatin)                   
    real(r8), intent(in) :: xout(nlondout,nlatout)        
    real(r8), intent(in) :: yout(nlatout)                 
    
    
    
    real(r8), intent(out) :: arrout(nlondout,nlevdout,nlatout) 
    
    
    
    integer :: i, ii, iw, ie, iiprev 
    integer :: j, jj, js, jn, jjprev 
    integer :: k                     
    integer :: icount                

    real(r8) :: extrap               
    real(r8) :: dxinwrap             
    real(r8) :: avgdxin              
    real(r8) :: ratio                
    real(r8) :: sum                  
    
    
    
    integer :: iim(nlondout)         
    integer :: iip(nlondout)         
    integer :: jjm(nlatout)          
    integer :: jjp(nlatout)          

    real(r8) :: wgts(nlatout)        
    real(r8) :: wgtn(nlatout)        
    real(r8) :: wgte(nlondout)       
    real(r8) :: wgtw(nlondout)       
    real(r8) :: igrid(nlonin)        

    character(len=132) :: message
    
    
    
    
    if (nlonin < 2 .or. nlatin < 2) then
       call wrf_abort  
    end if

    if (xin(1) < 0._r8 .or. xin(nlonin) > 360._r8) then
       call wrf_abort  
    end if

    icount = 0
    do i=1,nlonin-1
       if (xin(i) >= xin(i+1)) icount = icount + 1
    end do

    do j=1,nlatin-1
       if (yin(j) >= yin(j+1)) icount = icount + 1
    end do

    do j=1,nlatout-1
       if (yout(j) >= yout(j+1)) icount = icount + 1
    end do

    do j=1,nlatout
       do i=1,nlonout(j)-1
          if (xout(i,j) >= xout(i+1,j)) icount = icount + 1
       end do
    end do

    if (icount > 0) then
       call wrf_abort  
    end if

    if (yout(nlatout) <= yin(1) .or. yout(1) >= yin(nlatin)) then
       call wrf_abort  
    end if

    do j=1,nlatout
       if (xout(1,j) < 0._r8 .or. xout(nlonout(j),j) > 360._r8) then
          call wrf_abort  
       end if

       if (xout(nlonout(j),j) <= xin(1) .or.  &
            xout(1,j)          >= xin(nlonin)) then
          call wrf_abort  
       end if
    end do
    
    
    
    do j=1,nlatout
       jjm(j) = 0
       jjp(j) = 0
    end do
    
    
    
    
    do js=1,nlatout
       if (yout(js) > yin(1)) exit
       jjm(js) = 1
       jjp(js) = 1
       wgts(js) = 1._r8
       wgtn(js) = 0._r8
    end do

    do jn=nlatout,1,-1
       if (yout(jn) <= yin(nlatin)) exit
       jjm(jn) = nlatin
       jjp(jn) = nlatin
       wgts(jn) = 1._r8
       wgtn(jn) = 0._r8
    end do
    
    
    
    jjprev = 1
    do j=js,jn
       do jj=jjprev,nlatin-1
          if (yout(j) > yin(jj) .and. yout(j) <= yin(jj+1)) then
             jjm(j) = jj
             jjp(j) = jj + 1
             wgts(j) = (yin(jj+1) - yout(j)) / (yin(jj+1) - yin(jj))
             wgtn(j) = (yout(j)   - yin(jj)) / (yin(jj+1) - yin(jj))
             goto 30
          end if
       end do
       call wrf_abort  
30     jjprev = jj
    end do

    dxinwrap = xin(1) + 360._r8 - xin(nlonin)
    
    
    
    avgdxin = (xin(nlonin)-xin(1))/(nlonin-1._r8)
    ratio = dxinwrap/avgdxin
    if (ratio < 0.9_r8 .or. ratio > 1.1_r8) then
       write(message,*)'BILIN: Insane dxinwrap value =',dxinwrap,' avg=', avgdxin
       call wrf_message( trim(message) )
       call wrf_abort  
    end if
    
    
    
    extrap = 100._r8*((js - 1._r8) + real(nlatout - jn,r8))/nlatout
    if (extrap > 20._r8) then
       write(message,*)'BILIN:',extrap,' % of N/S output grid will have to be extrapolated'
       call wrf_message( trim(message) )
    end if
    
    
    
    
    icount = 0
    do j=1,nlatout
       if (jjm(j) == 0 .or. jjp(j) == 0) icount = icount + 1
       sum = wgts(j) + wgtn(j)
       if (sum < 0.99999_r8 .or. sum > 1.00001_r8) icount = icount + 1
       if (wgts(j) < 0._r8 .or. wgts(j) > 1._r8) icount = icount + 1
       if (wgtn(j) < 0._r8 .or. wgtn(j) > 1._r8) icount = icount + 1
    end do

    if (icount > 0) then
       call wrf_abort  
    end if
    
    
    
    do j=1,nlatout
       
       
       
       do i=1,nlondout
          iim(i) = 0
          iip(i) = 0
       end do
       
       
       
       
       do iw=1,nlonout(j)
          if (xout(iw,j) > xin(1)) exit
          iim(iw) = nlonin
          iip(iw) = 1
          wgtw(iw) = (xin(1)        - xout(iw,j))   /dxinwrap
          wgte(iw) = (xout(iw,j)+360._r8 - xin(nlonin))/dxinwrap
       end do

       do ie=nlonout(j),1,-1
          if (xout(ie,j) <= xin(nlonin)) exit
          iim(ie) = nlonin
          iip(ie) = 1
          wgtw(ie) = (xin(1)+360._r8 - xout(ie,j))   /dxinwrap
          wgte(ie) = (xout(ie,j)    - xin(nlonin))/dxinwrap
       end do
       
       
       
       iiprev = 1
       do i=iw,ie
          do ii=iiprev,nlonin-1
             if (xout(i,j) > xin(ii) .and. xout(i,j) <= xin(ii+1)) then
                iim(i) = ii
                iip(i) = ii + 1
                wgtw(i) = (xin(ii+1) - xout(i,j)) / (xin(ii+1) - xin(ii))
                wgte(i) = (xout(i,j) - xin(ii))   / (xin(ii+1) - xin(ii))
                goto 60
             end if
          end do
          call wrf_abort  
60        iiprev = ii
       end do

       icount = 0
       do i=1,nlonout(j)
          if (iim(i) == 0 .or. iip(i) == 0) icount = icount + 1
          sum = wgtw(i) + wgte(i)
          if (sum < 0.99999_r8 .or. sum > 1.00001_r8) icount = icount + 1
          if (wgtw(i) < 0._r8 .or. wgtw(i) > 1._r8) icount = icount + 1
          if (wgte(i) < 0._r8 .or. wgte(i) > 1._r8) icount = icount + 1
       end do

       if (icount > 0) then
          write(message,*)'BILIN: j=',j,' Something bad in longitude indices or weights'
          call wrf_message( trim(message) )
          call wrf_abort  
       end if
       
       
       
       do k=1,nlev
          do i=1,nlonin
             igrid(i) = arrin(i,k,jjm(j))*wgts(j) + arrin(i,k,jjp(j))*wgtn(j)
          end do

          do i=1,nlonout(j)
             arrout(i,k,j) = igrid(iim(i))*wgtw(i) + igrid(iip(i))*wgte(i)
          end do
       end do
    end do


    return
  end subroutine bilin

  subroutine vertinterp(ncol, ncold, nlev, pmid, pout, arrin, arrout)

    
    
    
    
    
    
    
    
    
    
    

    implicit none

    
    integer , intent(in)  :: ncol              
    integer , intent(in)  :: ncold             
    integer , intent(in)  :: nlev              
    real(r8), intent(in)  :: pmid(ncold,nlev)  
    real(r8), intent(in)  :: pout              
    real(r8), intent(in)  :: arrin(ncold,nlev) 
    real(r8), intent(out) :: arrout(ncold)     
    

    
    integer i,k               
    integer kupper(ncold)     
    real(r8) dpu              
    real(r8) dpl              
    logical found(ncold)      
    logical error             
    
    
    
    
    do i=1,ncol
       found(i)  = .false.
       kupper(i) = 1
    end do
    error = .false.
    
    
    
    
    
    do k=1,nlev-1
       do i=1,ncol
          if ((.not. found(i)) .and. pmid(i,k)<pout .and. pout<=pmid(i,k+1)) then
             found(i) = .true.
             kupper(i) = k
          end if
       end do
    end do
    

    
    
    
    do i=1,ncol
       if (pout <= pmid(i,1)) then
          arrout(i) = arrin(i,1)
       else if (pout >= pmid(i,nlev)) then
          arrout(i) = arrin(i,nlev)
       else if (found(i)) then
          dpu = pout - pmid(i,kupper(i))
          dpl = pmid(i,kupper(i)+1) - pout
          arrout(i) = (arrin(i,kupper(i)  )*dpl + arrin(i,kupper(i)+1)*dpu)/(dpl + dpu)
       else
          error = .true.
       end if
    end do
    
    
    
    if (error) then
       call wrf_abort  
    end if

    return
  end subroutine vertinterp

  subroutine get_timeinterp_factors (cycflag, np1, cdayminus, cdayplus, cday, &
       fact1, fact2, str)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    implicit none
    
    
    
    logical, intent(in) :: cycflag             

    integer, intent(in) :: np1                 

    real(r8), intent(in) :: cdayminus          
    real(r8), intent(in) :: cdayplus           
    real(r8), intent(in) :: cday               
    real(r8), intent(out) :: fact1             
    real(r8), intent(out) :: fact2             

    character(len=*), intent(in) :: str        
    
    
    
    real(r8) :: deltat                         
    real(r8), parameter :: daysperyear = 365.  

    character(len=132) :: message
    
    
    
    if (np1 == 1 .and. .not. cycflag) then
       call wrf_abort  
    end if

    if (np1 < 1) then
       call wrf_abort  
    end if

    if (cycflag) then
       if ((cday < 1.) .or. (cday > (daysperyear+1.))) then
          write(message,*) 'GETFACTORS:', str, ' bad cday=',cday
          call wrf_message( trim(message) )
          call wrf_abort  
       end if
    else
       if (cday < 1.) then
          write(message,*) 'GETFACTORS:', str, ' bad cday=',cday
          call wrf_message( trim(message) )
          call wrf_abort  
       end if
    end if
    
    
    
    
    if (cycflag .and. np1 == 1) then                     
       deltat = cdayplus + daysperyear - cdayminus
       if (cday > cdayplus) then                         
          fact1 = (cdayplus + daysperyear - cday)/deltat
          fact2 = (cday - cdayminus)/deltat
       else                                              
          fact1 = (cdayplus - cday)/deltat
          fact2 = (cday + daysperyear - cdayminus)/deltat
       end if
    else
       deltat = cdayplus - cdayminus
       fact1 = (cdayplus - cday)/deltat
       fact2 = (cday - cdayminus)/deltat
    end if

    if (.not. valid_timeinterp_factors (fact1, fact2)) then
       write(message,*) 'GETFACTORS: ', str, ' bad fact1 and/or fact2=', fact1, fact2
       call wrf_message( trim(message) )
       call wrf_abort  
    end if

    return
  end subroutine get_timeinterp_factors

  logical function valid_timeinterp_factors (fact1, fact2)
    
    
    
    
    
    implicit none

    real(r8), intent(in) :: fact1, fact2           

    valid_timeinterp_factors = .true.

    
    if (abs(fact1+fact2-1.) > 1.e-6 .or. &
         fact1 > 1.000001 .or. fact1 < -1.e-6 .or. &
         fact2 > 1.000001 .or. fact2 < -1.e-6 .or. &
         fact1 .ne. fact1 .or. fact2 .ne. fact2) then

       valid_timeinterp_factors = .false.
    end if

    return
  end function valid_timeinterp_factors

end module module_interpolate
