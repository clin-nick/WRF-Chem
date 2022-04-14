













   module modal_aero_wateruptake


   use shr_kind_mod, only : r8 => shr_kind_r8
   use module_cam_support, only: iulog

   implicit none
   private
   save
                                                                                                                             

   public modal_aero_wateruptake_sub
                                                                                                                             



                                                                                                                             












                                                                                                                             
                                                                                                                             
   contains
                                                                                                                             
                                                                                                                             

      subroutine modal_aero_wateruptake_sub(                &
                 lchnk, ncol, nstep,                        &
                 iwaterup_flag, loffset,                    &
                 aero_mmr_flag, h2o_mmr_flag,               &
                 deltat, h2ommr, t, pmid, pdel, cldn,       &
                 raer, raertend, dotend, qaerwat,           &
                 dgncur_a, dgncur_awet, wetdens             &
                                                            )












      use modal_aero_data
      use physconst,     only: cpair, epsilo, gravit, mwdry, mwh2o,   &
                               rair, rga, rhoh2o, rh2o, latvap
      use wv_saturation, only: aqsat, qsat_water
      use module_cam_support, only: pcnst => pcnst_runtime, &
                                    pcols, pver, &
                                    endrun, masterproc, outfld
      use physconst,     only: pi

      implicit none

      integer,  intent(in)  :: lchnk              
      integer,  intent(in)  :: ncol               
      integer,  intent(in)  :: nstep              
      integer,  intent(in)  :: iwaterup_flag       
                               
      integer,  intent(in)  :: loffset            

      logical,  intent(in)  :: aero_mmr_flag      
                                                  
      logical,  intent(in)  :: h2o_mmr_flag       
      logical,  intent(inout)::dotend(pcnst)
                               

      real(r8), intent(in)  :: deltat             
      real(r8), intent(in)  :: h2ommr(pcols,pver) 
      real(r8), intent(in)  :: t(pcols,pver)      
      real(r8), intent(in)  :: pmid(pcols,pver)   
      real(r8), intent(in)  :: pdel(pcols,pver)   
      real(r8), intent(in)  :: cldn(pcols,pver)   
      real(r8), intent(in)  :: raer(pcols,pver,pcnst)
                               
      real(r8), intent(inout)::raertend(pcols,pver,pcnst)
                               
                               
      real(r8), intent(out)   :: qaerwat(pcols,pver,ntot_amode)
      real(r8), intent(in)    :: dgncur_a(pcols,pver,ntot_amode)
      real(r8), intent(out)   :: dgncur_awet(pcols,pver,ntot_amode)
      real(r8), intent(out)   :: wetdens(pcols,pver,ntot_amode)













      integer i,k,m
      integer icol_diag
      integer l 
      integer lmass 
      integer ltype 

      integer  lat(pcols), lon(pcols)      

      real(r8) density_water                   
      real(r8) drydens(ntot_amode)   
      real(r8) drymass(ntot_amode)   
      real(r8) dryrad(pcols,pver,ntot_amode) 
      real(r8) dryvol(ntot_amode)    
      real(r8) dryvolmr(ntot_amode)  
      real(r8) duma, dumb
      real(r8) es(pcols,pver)        
      real(r8) hygro(ntot_amode)     
      real(r8) hystfac(ntot_amode)   
      real(r8) pi43
      real(r8) qs(pcols,pver)        
      real(r8) qwater                
      real(r8) rh(pcols,pver)        
      real(r8) third
      real(r8) v2ncur_a(pcols,pver,ntot_amode)
      real(r8) wtrvol(ntot_amode)    
      real(r8) wetvol(ntot_amode)    

      real(r8) :: maer(pcols,pver,ntot_amode)
                              
      real(r8) :: naer(pcols,pver,ntot_amode)
                              
      real(r8) :: wetrad(pcols,pver,ntot_amode)  
                              

      character(len=3) :: trnum       









      do m=1,ntot_amode



      end do


      third=1./3.
      pi43 = pi*4.0/3.0
      density_water = rhoh2o   






      do m=1,ntot_amode
           hystfac(m) = 1.0 / max( 1.0e-5_r8,   &
                              (rhdeliques_amode(m)-rhcrystal_amode(m)) )
      enddo




      do k=1,pver
      do i=1,ncol

         qs(i,k)=qsat_water(t(i,k),pmid(i,k))

         if ( h2o_mmr_flag ) then
            rh(i,k) = h2ommr(i,k)/qs(i,k)
         else
            rh(i,k) = h2ommr(i,k)*mwh2o/(mwdry*qs(i,k))
         end if
         rh(i,k) = max(rh(i,k),0.0_r8)
         rh(i,k) = min(rh(i,k),0.98_r8)
         if (cldn(i,k) .lt. 1.0_r8) then
           rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_r8 - cldn(i,k))  
         end if
         rh(i,k) = max(rh(i,k),0.0_r8)
           


         do m=1,ntot_amode

            maer(i,k,m)=0.
            dryvolmr(m)=0.
            hygro(m)=0.
            do l = 1, nspec_amode(m)
               lmass = lmassptr_amode(l,m) - loffset
               ltype = lspectype_amode(l,m)
               if ( aero_mmr_flag ) then
                  duma = raer(i,k,lmass)
               else
                  duma = raer(i,k,lmass)*(specmw_amode(ltype)/mwdry)
               end if
               maer(i,k,m) = maer(i,k,m) + duma
               dumb = duma/specdens_amode(ltype)
               dryvolmr(m) = dryvolmr(m) + dumb
               hygro(m) = hygro(m) + dumb*spechygro(ltype)
            enddo
            if (dryvolmr(m) > 1.0e-30_r8) then
               hygro(m) = hygro(m)/dryvolmr(m)
            else
               hygro(m) = spechygro( lspectype_amode(1,m) )
            end if



            v2ncur_a(i,k,m) = 1. / ( (pi/6.)*                            &
                (dgncur_a(i,k,m)**3.)*exp(4.5*alnsg_amode(m)**2.) )
            naer(i,k,m) = dryvolmr(m)*v2ncur_a(i,k,m)
         enddo   






         do m=1,ntot_amode
            if (maer(i,k,m) .gt. 1.0e-31) then
               drydens(m) = maer(i,k,m)/dryvolmr(m)
            else
               drydens(m) = 1.0
            end if
            dryvol(m) = 1.0/v2ncur_a(i,k,m)
            drymass(m) = drydens(m)*dryvol(m)
            dryrad(i,k,m) = (dryvol(m)/pi43)**third
         enddo


         do m=1,ntot_amode
            call modal_aero_kohler(                          &
                    dryrad(i,k,m), hygro(m), rh(i,k),        &
                    wetrad(i,k,m), 1, 1                      )

            wetrad(i,k,m)=max(wetrad(i,k,m),dryrad(i,k,m))
            dgncur_awet(i,k,m) = dgncur_a(i,k,m)*   &
                                     (wetrad(i,k,m)/dryrad(i,k,m))
            wetvol(m) = pi43*wetrad(i,k,m)*wetrad(i,k,m)*wetrad(i,k,m)
            wetvol(m) = max(wetvol(m),dryvol(m))
            wtrvol(m) = wetvol(m) - dryvol(m)
            wtrvol(m) = max( wtrvol(m), 0.0_r8 )




            if (rh(i,k) < rhcrystal_amode(m)) then
               wetrad(i,k,m) = dryrad(i,k,m)
               wetvol(m) = dryvol(m)
               wtrvol(m) = 0.0_r8
            else if (rh(i,k) < rhdeliques_amode(m)) then
               wtrvol(m) = wtrvol(m)*hystfac(m)   &
                                    *(rh(i,k) - rhcrystal_amode(m))
               wtrvol(m) = max( wtrvol(m), 0.0_r8 )
               wetvol(m) = dryvol(m) + wtrvol(m)
               wetrad(i,k,m) = (wetvol(m)/pi43)**third
            end if




            if ( aero_mmr_flag ) then
               duma = 1.0_r8
            else
               duma = mwdry/mwh2o
            end if
            qwater = density_water*naer(i,k,m)*wtrvol(m)*duma








            qaerwat(i,k,m) = qwater


            if (wetvol(m) > 1.0e-30_r8) then
               wetdens(i,k,m) = (drymass(m) + density_water*wtrvol(m))/wetvol(m)
            else
               wetdens(i,k,m) = specdens_amode( lspectype_amode(1,m) )
            end if

         enddo


      end do   
      end do   



      do m = 1, ntot_amode
         write( trnum, '(i3.3)' ) m
         call outfld( 'wat_a'//trnum(3:3),  qaerwat(:,:,m),     pcols, lchnk)

         call outfld( 'dgnd_a'//trnum(2:3), dgncur_a(:,:,m),    pcols, lchnk)
         call outfld( 'dgnw_a'//trnum(2:3), dgncur_awet(:,:,m), pcols, lchnk)
      end do


      return
      end subroutine modal_aero_wateruptake_sub



      subroutine modal_aero_kohler(   &
          rdry_in, hygro, s, rwet_out, im, imx )







      implicit none


      integer :: im         
      integer :: imx        
      real(r8) :: rdry_in(imx)    
      real(r8) :: hygro(imx)      
      real(r8) :: s(imx)          
      real(r8) :: rwet_out(imx)   


      integer, parameter :: imax=200
      integer :: i, n, nsol

      real(r8) :: a, b
      real(r8) :: p40(imax),p41(imax),p42(imax),p43(imax) 
      real(r8) :: p30(imax),p31(imax),p32(imax) 
      real(r8) :: p
      real(r8) :: r3, r4
      real(r8) :: r(imx)        
      real(r8) :: rdry(imax)    
      real(r8) :: ss            
      real(r8) :: slog(imax)    
      real(r8) :: vol(imax)     
      real(r8) :: xi, xr

      complex(r8) :: cx4(4,imax),cx3(3,imax)

      real(r8), parameter :: eps = 1.e-4
      real(r8), parameter :: mw = 18.
      real(r8), parameter :: pi = 3.14159
      real(r8), parameter :: rhow = 1.
      real(r8), parameter :: surften = 76.
      real(r8), parameter :: tair = 273.
      real(r8), parameter :: third = 1./3.
      real(r8), parameter :: ugascon = 8.3e7



      a=2.e4*mw*surften/(ugascon*tair*rhow)

      do i=1,im
           rdry(i) = rdry_in(i)*1.0e6   
           vol(i) = rdry(i)**3          
           b = vol(i)*hygro(i)


           ss=min(s(i),1.-eps)
           ss=max(ss,1.e-10_r8)
           slog(i)=log(ss)
           p43(i)=-a/slog(i)
           p42(i)=0.
           p41(i)=b/slog(i)-vol(i)
           p40(i)=a*vol(i)/slog(i)

           p32(i)=0.
           p31(i)=-b/a
           p30(i)=-vol(i)
      end do


       do 100 i=1,im


        if(vol(i).le.1.e-12)then
           r(i)=rdry(i)
           go to 100
        endif

        p=abs(p31(i))/(rdry(i)*rdry(i))
        if(p.lt.eps)then

           r(i)=rdry(i)*(1.+p*third/(1.-slog(i)*rdry(i)/a))
        else
           call makoh_quartic(cx4(1,i),p43(i),p42(i),p41(i),p40(i),1)

           r(i)=1000.*rdry(i)
           nsol=0
           do n=1,4
              xr=real(cx4(n,i))
              xi=imag(cx4(n,i))
              if(abs(xi).gt.abs(xr)*eps) cycle  
              if(xr.gt.r(i)) cycle  
              if(xr.lt.rdry(i)*(1.-eps)) cycle  
              if(xr.ne.xr) cycle  
              r(i)=xr
              nsol=n
           end do  
           if(nsol.eq.0)then
              write(iulog,*)   &
               'ccm kohlerc - no real(r8) solution found (quartic)'
              call wrf_message(iulog)
              write(iulog,*)'roots =', (cx4(n,i),n=1,4)
              call wrf_message(iulog)
              write(iulog,*)'p0-p3 =', p40(i), p41(i), p42(i), p43(i)
              call wrf_message(iulog)
              write(iulog,*)'rh=',s(i)
              call wrf_message(iulog)
              write(iulog,*)'setting radius to dry radius=',rdry(i)
              call wrf_message(iulog)
              r(i)=rdry(i)

           endif
        endif

        if(s(i).gt.1.-eps)then

           r4=r(i)

           p=abs(p31(i))/(rdry(i)*rdry(i))
           if(p.lt.eps)then
              r(i)=rdry(i)*(1.+p*third)
           else
              call makoh_cubic(cx3,p32,p31,p30,im)

              r(i)=1000.*rdry(i)
              nsol=0
              do n=1,3
                 xr=real(cx3(n,i))
                 xi=imag(cx3(n,i))
                 if(abs(xi).gt.abs(xr)*eps) cycle  
                 if(xr.gt.r(i)) cycle  
                 if(xr.lt.rdry(i)*(1.-eps)) cycle  
                 if(xr.ne.xr) cycle  
                 r(i)=xr
                 nsol=n
              end do  
              if(nsol.eq.0)then
                 write(iulog,*)   &
                  'ccm kohlerc - no real(r8) solution found (cubic)'
                 call wrf_message(iulog)
                 write(iulog,*)'roots =', (cx3(n,i),n=1,3)
                 call wrf_message(iulog)
                 write(iulog,*)'p0-p2 =', p30(i), p31(i), p32(i)
                 call wrf_message(iulog)
                 write(iulog,*)'rh=',s(i)
                 call wrf_message(iulog)
                 write(iulog,*)'setting radius to dry radius=',rdry(i)
                 call wrf_message(iulog)
                 r(i)=rdry(i)

              endif
           endif
           r3=r(i)

           r(i)=(r4*(1.-s(i))+r3*(s(i)-1.+eps))/eps
        endif

  100 continue


      do i=1,im
         r(i) = min(r(i),30._r8) 
         rwet_out(i) = r(i)*1.e-6
      end do

      return
      end subroutine modal_aero_kohler



      subroutine makoh_cubic( cx, p2, p1, p0, im )




      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx)
      complex(r8) :: cx(3,imx)

      integer :: i
      real(r8) :: eps, q(imx), r(imx), sqrt3, third
      complex(r8) :: ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

      save eps
      data eps/1.e-20/

      third=1./3.
      ci=dcmplx(0.,1.)
      sqrt3=sqrt(3.)
      cw=0.5*(-1+ci*sqrt3)
      cwsq=0.5*(-1-ci*sqrt3)

      do i=1,im
      if(p1(i).eq.0.)then

         cx(1,i)=(-p0(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
      else
         q(i)=p1(i)/3.
         r(i)=p0(i)/2.
         crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
         crad(i)=sqrt(crad(i))

         cy(i)=r(i)-crad(i)
         if (abs(cy(i)).gt.eps) cy(i)=cy(i)**third
         cq=q(i)
         cz(i)=-cq/cy(i)

         cx(1,i)=-cy(i)-cz(i)
         cx(2,i)=-cw*cy(i)-cwsq*cz(i)
         cx(3,i)=-cwsq*cy(i)-cw*cz(i)
      endif
      enddo

      return
      end subroutine makoh_cubic



      subroutine makoh_quartic( cx, p3, p2, p1, p0, im )




      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx), p3(imx)
      complex(r8) :: cx(4,imx)

      integer :: i
      real(r8) :: third, q(imx), r(imx)
      complex(r8) :: cb(imx), cb0(imx), cb1(imx),   &
                     crad(imx), cy(imx), czero


      czero=cmplx(0.0,0.0)
      third=1./3.

      do 10 i=1,im

      q(i)=-p2(i)*p2(i)/36.+(p3(i)*p1(i)-4*p0(i))/12.
      r(i)=-(p2(i)/6)**3+p2(i)*(p3(i)*p1(i)-4*p0(i))/48.   &
       +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16

      crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
      crad(i)=sqrt(crad(i))

      cb(i)=r(i)-crad(i)
      if(cb(i).eq.czero)then

         cx(1,i)=(-p1(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
         cx(4,i)=cx(1,i)
      else
         cb(i)=cb(i)**third

         cy(i)=-cb(i)+q(i)/cb(i)+p2(i)/6

         cb0(i)=sqrt(cy(i)*cy(i)-p0(i))
         cb1(i)=(p3(i)*cy(i)-p1(i))/(2*cb0(i))

         cb(i)=p3(i)/2+cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)+cb0(i))
         crad(i)=sqrt(crad(i))
         cx(1,i)=(-cb(i)+crad(i))/2.
         cx(2,i)=(-cb(i)-crad(i))/2.

         cb(i)=p3(i)/2-cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)-cb0(i))
         crad(i)=sqrt(crad(i))
         cx(3,i)=(-cb(i)+crad(i))/2.
         cx(4,i)=(-cb(i)-crad(i))/2.
      endif
   10 continue

      return
      end subroutine makoh_quartic



   end module modal_aero_wateruptake


