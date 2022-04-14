
module wetdep







use shr_kind_mod, only: r8 => shr_kind_r8
use physconst,    only: gravit, rair, tmelt
use module_cam_support, only: pcols, pver, iulog

implicit none
save
private

public :: wetdepa_v1  
public :: wetdepa_v2  
public :: clddiag     

real(r8), parameter :: cmftau = 3600._r8
real(r8), parameter :: rhoh2o = 1000._r8            
real(r8), parameter :: molwta = 28.97_r8            


contains


subroutine clddiag(t, pmid, pdel, cmfdqr, evapc, &
                   cldt, cldcu, cldst, cme, evapr, &
                   prain, cldv, cldvcu, cldvst, rain, &
                   ncol)

   
   
   
   
   
   
   
   
   

   
   real(r8), intent(in) :: t(pcols,pver)        
   real(r8), intent(in) :: pmid(pcols,pver)     
   real(r8), intent(in) :: pdel(pcols,pver)     
   real(r8), intent(in) :: cmfdqr(pcols,pver)   
   real(r8), intent(in) :: evapc(pcols,pver)    
   real(r8), intent(in) :: cldt(pcols,pver)    
   real(r8), intent(in) :: cldcu(pcols,pver)    
   real(r8), intent(in) :: cldst(pcols,pver)    
   real(r8), intent(in) :: cme(pcols,pver)      
   real(r8), intent(in) :: evapr(pcols,pver)    
   real(r8), intent(in) :: prain(pcols,pver)    
   integer, intent(in) :: ncol

   
   real(r8), intent(out) :: cldv(pcols,pver)     
   real(r8), intent(out) :: cldvcu(pcols,pver)   
   real(r8), intent(out) :: cldvst(pcols,pver)   
   real(r8), intent(out) :: rain(pcols,pver)     

   
   integer  i, k
   real(r8) convfw         
   real(r8) sumppr(pcols)        
   real(r8) sumpppr(pcols)       
   real(r8) cldv1(pcols)         
   real(r8) lprec                
   real(r8) lprecp               
   real(r8) rho                  
   real(r8) vfall
   real(r8) sumppr_cu(pcols)     
   real(r8) sumpppr_cu(pcols)    
   real(r8) cldv1_cu(pcols)      
   real(r8) lprec_cu             
   real(r8) lprecp_cu            
   real(r8) sumppr_st(pcols)     
   real(r8) sumpppr_st(pcols)    
   real(r8) cldv1_st(pcols)      
   real(r8) lprec_st             
   real(r8) lprecp_st            
   

   convfw = 1.94_r8*2.13_r8*sqrt(rhoh2o*gravit*2.7e-4_r8)
   do i=1,ncol
      sumppr(i) = 0._r8
      cldv1(i) = 0._r8
      sumpppr(i) = 1.e-36_r8
      sumppr_cu(i)  = 0._r8
      cldv1_cu(i)   = 0._r8
      sumpppr_cu(i) = 1.e-36_r8
      sumppr_st(i)  = 0._r8
      cldv1_st(i)   = 0._r8
      sumpppr_st(i) = 1.e-36_r8
   end do

   do k = 1,pver
      do i = 1,ncol
         cldv(i,k) = &
            max(min(1._r8, &
            cldv1(i)/sumpppr(i) &
            )*sumppr(i)/sumpppr(i), &
            cldt(i,k) &
            )
         lprec = pdel(i,k)/gravit &
            *(prain(i,k)+cmfdqr(i,k)-evapr(i,k))
         lprecp = max(lprec,1.e-30_r8)
         cldv1(i) = cldv1(i)  + cldt(i,k)*lprecp
         sumppr(i) = sumppr(i) + lprec
         sumpppr(i) = sumpppr(i) + lprecp

         
         cldvcu(i,k)   = max(min(1._r8,cldv1_cu(i)/sumpppr_cu(i))*(sumppr_cu(i)/sumpppr_cu(i)),0._r8)
         lprec_cu      = (pdel(i,k)/gravit)*(cmfdqr(i,k)-evapc(i,k))
         lprecp_cu     = max(lprec_cu,1.e-30_r8)
         cldv1_cu(i)   = cldv1_cu(i) + cldcu(i,k)*lprecp_cu
         sumppr_cu(i)  = sumppr_cu(i) + lprec_cu
         sumpppr_cu(i) = sumpppr_cu(i) + lprecp_cu

         
         cldvst(i,k)   = max(min(1._r8,cldv1_st(i)/sumpppr_st(i))*(sumppr_st(i)/sumpppr_st(i)),0._r8)
         lprec_st      = (pdel(i,k)/gravit)*(prain(i,k)-evapr(i,k))
         lprecp_st     = max(lprec_st,1.e-30_r8)
         cldv1_st(i)   = cldv1_st(i) + cldst(i,k)*lprecp_st
         sumppr_st(i)  = sumppr_st(i) + lprec_st
         sumpppr_st(i) = sumpppr_st(i) + lprecp_st

         rain(i,k) = 0._r8
         if(t(i,k) .gt. tmelt) then
            rho = pmid(i,k)/(rair*t(i,k))
            vfall = convfw/sqrt(rho)
            rain(i,k) = sumppr(i)/(rho*vfall)
            if (rain(i,k).lt.1.e-14_r8) rain(i,k) = 0._r8
         endif
      end do
   end do

end subroutine clddiag





subroutine wetdepa_v2(t, p, q, pdel, &
                   cldt, cldc, cmfdqr, evapc, conicw, precs, conds, &
                       evaps, cwat, tracer, deltat, &
                       scavt, iscavt, cldv, cldvcu, cldvst, dlf, fracis, sol_fact, ncol, &
                       scavcoef,  &
                       is_strat_cloudborne, rate1ord_cw2pr_st, qqcw, f_act_conv, &     
                       icscavt, isscavt, bcscavt, bsscavt, &
                       sol_facti_in, sol_factbi_in, sol_factii_in, &
                       sol_factic_in, sol_factiic_in )

      
      
      
      
      
      
      
      



      implicit none

      real(r8), intent(in) ::&
         t(pcols,pver),        &
         p(pcols,pver),        &
         q(pcols,pver),        &
         pdel(pcols,pver),     &
         cldt(pcols,pver),    &
         cldc(pcols,pver),     &
         cmfdqr(pcols,pver),   &

         evapc(pcols,pver),    &

         conicw(pcols,pver),   &
         cwat(pcols,pver),     &
         precs(pcols,pver),    &
         conds(pcols,pver),    &
         evaps(pcols,pver),    &
         cldv(pcols,pver),     &

         cldvcu(pcols,pver),   &
         cldvst(pcols,pver),   &
         dlf(pcols,pver),      &

         deltat,               &
         tracer(pcols,pver)     
      
            
      
            
            
         real(r8), intent(in) :: sol_fact 
         integer, intent(in) :: ncol
         real(r8), intent(in) :: scavcoef(pcols,pver) 
         
         
         logical, intent(in), optional :: is_strat_cloudborne   
         
         real(r8), intent(in), optional  :: rate1ord_cw2pr_st(pcols,pver)
         
         real(r8), intent(in), optional  :: qqcw(pcols,pver)
         
         real(r8), intent(in), optional  :: f_act_conv(pcols,pver)
         

         real(r8), intent(in), optional :: sol_facti_in   
         real(r8), intent(in), optional :: sol_factbi_in  
         real(r8), intent(in), optional :: sol_factii_in  
         real(r8), intent(in), optional :: sol_factic_in(pcols,pver)  
         real(r8), intent(in), optional :: sol_factiic_in 
         
      real(r8), intent(out) ::&
         scavt(pcols,pver),    &
         iscavt(pcols,pver),   &
         fracis(pcols,pver)     

      real(r8), intent(out), optional ::    icscavt(pcols,pver)     
      real(r8), intent(out), optional ::    isscavt(pcols,pver)     
      real(r8), intent(out), optional ::    bcscavt(pcols,pver)     
      real(r8), intent(out), optional ::    bsscavt(pcols,pver)     

      

      integer i                 
      integer k                 

      real(r8) adjfac               
      real(r8) aqfrac               
      real(r8) cwatc                
      real(r8) cwats                
      real(r8) cwatp                
      real(r8) fracev(pcols)        

      real(r8) fracev_cu(pcols)     

      real(r8) fracp                
      real(r8) gafrac               
      real(r8) hconst               
                                
      real(r8) mpla                 
      real(r8) mplb                 
      real(r8) omsm                 
      real(r8) part                 
      real(r8) patm                 
      real(r8) pdog                 
      real(r8) precabc(pcols)       
      real(r8) precabs(pcols)       
      real(r8) precbl               
      real(r8) precmin              
      real(r8) rat(pcols)           
      real(r8) scavab(pcols)        
      real(r8) scavabc(pcols)       
      real(r8) srcc                 
      real(r8) srcs                 
      real(r8) srct(pcols)          
      real(r8) tracab(pcols)        

      real(r8) fins                 
      real(r8) finc                 
      real(r8) srcs1                
      real(r8) srcs2                
      real(r8) tc                   
      real(r8) weight               
      real(r8) cldmabs(pcols)       
      real(r8) cldmabc(pcols)       
      real(r8) odds                 
      real(r8) dblchek(pcols)
      logical :: found

    
    
    

      real(r8) tracer_incu
      real(r8) tracer_mean

    

      real(r8) sol_facti,  sol_factb  
      real(r8) sol_factii, sol_factbi 
      real(r8) sol_factic(pcols,pver)             
      real(r8) sol_factiic            
      
      
      
      
      

      

      omsm = 1._r8-2*epsilon(1._r8) 
      precmin =  0.1_r8/8.64e4_r8      

      adjfac = deltat/(max(deltat,cmftau)) 

      

	

      sol_facti = sol_fact
      sol_factb = sol_fact
      sol_factii = sol_fact
      sol_factbi = sol_fact

      if ( present(sol_facti_in) )  sol_facti = sol_facti_in
      if ( present(sol_factii_in) )  sol_factii = sol_factii_in
      if ( present(sol_factbi_in) )  sol_factbi = sol_factbi_in

      sol_factic  = sol_facti
      sol_factiic = sol_factii
      if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in
      if ( present(sol_factiic_in) )  sol_factiic = sol_factiic_in

      
      
      
      
      
      
      
      

      do i = 1,pcols
         precabs(i) = 0
         precabc(i) = 0
         scavab(i) = 0
         scavabc(i) = 0
         tracab(i) = 0
         cldmabs(i) = 0
         cldmabc(i) = 0

       
       
       
       
       
       
       

      end do

      do k = 1,pver
         do i = 1,ncol
            tc     = t(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) 
            weight = 0._r8                                 

            pdog = pdel(i,k)/gravit

            
            
            
            fracev(i) = evaps(i,k)*pdel(i,k)/gravit &
                     /max(1.e-12_r8,precabs(i))

            
            fracev(i) = max(0._r8,min(1._r8,fracev(i)))


            fracev_cu(i) = evapc(i,k)*pdel(i,k)/gravit/max(1.e-12_r8,precabc(i))
            fracev_cu(i) = max(0._r8,min(1._r8,fracev_cu(i)))

            
            

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

            
            
            
            


            
              fracp = cmfdqr(i,k)*deltat/max(1.e-12_r8,cldc(i,k)*conicw(i,k)+(cmfdqr(i,k)+dlf(i,k))*deltat) 
            
            
            
            fracp = max(min(1._r8,fracp),0._r8)
            




            
            if ( is_strat_cloudborne ) then
               
               srcs1 = 0._r8
            else
               tracer_incu = f_act_conv(i,k)*(tracer(i,k)+& 
                             min(qqcw(i,k),tracer(i,k)*((cldt(i,k)-cldc(i,k))/max(0.01_r8,(1._r8-(cldt(i,k)-cldc(i,k)))))))              
               srcs1 = sol_factic(i,k)*cldc(i,k)*fracp*tracer_incu*(1._r8-weight)/deltat &  
                     + sol_factiic    *cldc(i,k)*fracp*tracer_incu*(weight)/deltat          
            end if


            

            

            
            
            
            
            cldmabc(i) = cldvcu(i,k)

            
            
            

            
            if ( is_strat_cloudborne ) then
               
               srcs2 = 0._r8
            else
               tracer_mean = tracer(i,k)*(1._r8-cldc(i,k)*f_act_conv(i,k))-cldc(i,k)*f_act_conv(i,k)*&
                             min(qqcw(i,k),tracer(i,k)*((cldt(i,k)-cldc(i,k))/max(0.01_r8,(1._r8-(cldt(i,k)-cldc(i,k))))))
               tracer_mean = max(0._r8,tracer_mean) 
               odds  = max(min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8)*scavcoef(i,k)*deltat),0._r8) 
               srcs2 = sol_factb *cldmabc(i)*odds*tracer_mean*(1._r8-weight)/deltat & 
                     + sol_factbi*cldmabc(i)*odds*tracer_mean*(weight)/deltat         
            end if





            srcc = srcs1 + srcs2  
            finc = srcs1/(srcc + 1.e-36_r8) 

            
            

            

            
            if ( is_strat_cloudborne ) then
               
               
               
             
             
               fracp = precs(i,k)*deltat/max(cwat(i,k)+precs(i,k)*deltat,1.e-12_r8) 
               fracp = max(0._r8,min(1._r8,fracp))
               srcs1 = sol_facti *fracp*tracer(i,k)/deltat*(1._r8-weight) &  
                     + sol_factii*fracp*tracer(i,k)/deltat*(weight)          
            else
               
               srcs1 = 0._r8
            end if
            


            




            cldmabs(i) = cldvst(i,k) 

            
            
            

            
            if ( is_strat_cloudborne ) then
               
               srcs2 = 0._r8
            else
               odds = precabs(i)/max(cldmabs(i),1.e-5_r8)*scavcoef(i,k)*deltat
               odds = max(min(1._r8,odds),0._r8)
               srcs2 = sol_factb *cldmabs(i)*odds*tracer_mean*(1._r8-weight)/deltat & 
                     + sol_factbi*cldmabs(i)*odds*tracer_mean*(weight)/deltat         
            end if


            srcs = srcs1 + srcs2             
            fins=srcs1/(srcs + 1.e-36_r8)    

            
            
            rat(i) = tracer(i,k)/max(deltat*(srcc+srcs),1.e-36_r8)
            if (rat(i).lt.1._r8) then
               srcs = srcs*rat(i)
               srcc = srcc*rat(i)
            endif
            srct(i) = (srcc+srcs)*omsm

            
            
            
            fracp = deltat*srct(i)/max(cldmabs(i)*tracer(i,k),1.e-36_r8)  
            fracp = max(0._r8,min(1._r8,fracp))
            fracis(i,k) = 1._r8 - fracp

            
            
         
            scavt(i,k) = -srct(i) + (fracev(i)*scavab(i)+fracev_cu(i)*scavabc(i))*gravit/pdel(i,k)
            iscavt(i,k) = -(srcc*finc + srcs*fins)*omsm

            if ( present(icscavt) ) icscavt(i,k) = -(srcc*finc) * omsm
            if ( present(isscavt) ) isscavt(i,k) = -(srcs*fins) * omsm
            if ( present(bcscavt) ) bcscavt(i,k) = -(srcc * (1-finc)) * omsm +  &
                 fracev_cu(i)*scavabc(i)*gravit/pdel(i,k)
            if ( present(bsscavt) ) bsscavt(i,k) = -(srcs * (1-fins)) * omsm +  &
                 fracev(i)*scavab(i)*gravit/pdel(i,k)

            dblchek(i) = tracer(i,k) + deltat*scavt(i,k)

            
            scavab(i) = scavab(i)*(1-fracev(i)) + srcs*pdel(i,k)/gravit
            precabs(i) = precabs(i) + (precs(i,k) - evaps(i,k))*pdel(i,k)/gravit
            scavabc(i) = scavabc(i)*(1-fracev_cu(i)) + srcc*pdel(i,k)/gravit
            precabc(i) = precabc(i) + (cmfdqr(i,k) - evapc(i,k))*pdel(i,k)/gravit
            tracab(i) = tracab(i) + tracer(i,k)*pdel(i,k)/gravit

       
       

       

       
       

       
       

       

         end do 

         found = .false.
         do i = 1,ncol
            if ( dblchek(i) < 0._r8 ) then
               found = .true.
               exit
            end if
         end do

         if ( found ) then
            do i = 1,ncol
               if (dblchek(i) .lt. 0._r8) then
                  write(iulog,*) ' wetdapa: negative value ', i, k, tracer(i,k), &
                       dblchek(i), scavt(i,k), srct(i), rat(i), fracev(i)
                  call wrf_message(iulog)
               endif
            end do
         endif

      end do 

   end subroutine wetdepa_v2







   subroutine wetdepa_v1( t, p, q, pdel, &
                       cldt, cldc, cmfdqr, conicw, precs, conds, &
                       evaps, cwat, tracer, deltat, &
                       scavt, iscavt, cldv, fracis, sol_fact, ncol, &
                       scavcoef,icscavt, isscavt, bcscavt, bsscavt, &
                       sol_facti_in, sol_factbi_in, sol_factii_in, &
                       sol_factic_in, sol_factiic_in )

      
      
      
      
      
      
      

      implicit none

      real(r8), intent(in) ::&
         t(pcols,pver),        &
         p(pcols,pver),        &
         q(pcols,pver),        &
         pdel(pcols,pver),     &
         cldt(pcols,pver),    &
         cldc(pcols,pver),     &
         cmfdqr(pcols,pver),   &
         conicw(pcols,pver),   &
         cwat(pcols,pver),     &
         precs(pcols,pver),    &
         conds(pcols,pver),    &
         evaps(pcols,pver),    &
         cldv(pcols,pver),     &
         deltat,               &
         tracer(pcols,pver)     
      
            
      
            
            
         real(r8), intent(in) :: sol_fact 
         real(r8), intent(in), optional :: sol_facti_in   
         real(r8), intent(in), optional :: sol_factbi_in  
         real(r8), intent(in), optional :: sol_factii_in  
         real(r8), intent(in), optional :: sol_factic_in(pcols,pver)  
         real(r8), intent(in), optional :: sol_factiic_in 
         real(r8), intent(in) :: scavcoef(pcols,pver) 
         
      integer, intent(in) :: ncol

      real(r8), intent(out) ::&
         scavt(pcols,pver),    &
         iscavt(pcols,pver),   &
         fracis(pcols,pver)     

      real(r8), intent(out), optional ::    icscavt(pcols,pver)     
      real(r8), intent(out), optional ::    isscavt(pcols,pver)     
      real(r8), intent(out), optional ::    bcscavt(pcols,pver)     
      real(r8), intent(out), optional ::    bsscavt(pcols,pver)     

      

      integer i                 
      integer k                 

      real(r8) adjfac               
      real(r8) aqfrac               
      real(r8) cwatc                
      real(r8) cwats                
      real(r8) cwatp                
      real(r8) fracev(pcols)        
      real(r8) fracp                
      real(r8) gafrac               
      real(r8) hconst               
                                
      real(r8) mpla                 
      real(r8) mplb                 
      real(r8) omsm                 
      real(r8) part                 
      real(r8) patm                 
      real(r8) pdog                 
      real(r8) precabc(pcols)       
      real(r8) precabs(pcols)       
      real(r8) precbl               
      real(r8) precmin              
      real(r8) rat(pcols)           
      real(r8) scavab(pcols)        
      real(r8) scavabc(pcols)       
      real(r8) srcc                 
      real(r8) srcs                 
      real(r8) srct(pcols)          
      real(r8) tracab(pcols)        

      real(r8) fins                 
      real(r8) finc                 
      real(r8) srcs1                
      real(r8) srcs2                
      real(r8) tc                   
      real(r8) weight               
      real(r8) cldmabs(pcols)       
      real(r8) cldmabc(pcols)       
      real(r8) odds                 
      real(r8) dblchek(pcols)
      logical :: found

      real(r8) sol_facti,  sol_factb  
      real(r8) sol_factii, sol_factbi 
      real(r8) sol_factic(pcols,pver)             
      real(r8) sol_factiic            
      
      
      
      
      

      

      omsm = 1._r8-2*epsilon(1._r8) 
      precmin =  0.1_r8/8.64e4_r8      

      adjfac = deltat/(max(deltat,cmftau)) 

      

	

      sol_facti = sol_fact
      sol_factb = sol_fact
      sol_factii = sol_fact
      sol_factbi = sol_fact

      if ( present(sol_facti_in) )  sol_facti = sol_facti_in
      if ( present(sol_factii_in) )  sol_factii = sol_factii_in
      if ( present(sol_factbi_in) )  sol_factbi = sol_factbi_in

      sol_factic  = sol_facti
      sol_factiic = sol_factii
      if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in
      if ( present(sol_factiic_in) )  sol_factiic = sol_factiic_in

      
      
      
      
      
      
      
      

      do i = 1,pcols
         precabs(i) = 0
         precabc(i) = 0
         scavab(i) = 0
         scavabc(i) = 0
         tracab(i) = 0
         cldmabs(i) = 0
         cldmabc(i) = 0
      end do

      do k = 1,pver
         do i = 1,ncol
            tc     = t(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) 
            weight = 0._r8                                 

            pdog = pdel(i,k)/gravit

            
            
            
            fracev(i) = evaps(i,k)*pdel(i,k)/gravit &
                     /max(1.e-12_r8,precabs(i))

            
            fracev(i) = max(0._r8,min(1._r8,fracev(i)))

            
            

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

            
            fracp = cmfdqr(i,k)*deltat/max(1.e-8_r8,conicw(i,k))
            
            fracp = max(min(1._r8,fracp),0._r8)
            



            srcs1 = sol_factic(i,k)*cldt(i,k)*fracp*tracer(i,k)*(1._r8-weight)/deltat &  
                 +  sol_factiic*cldt(i,k)*fracp*tracer(i,k)*(weight)/deltat      


            

            

            
            
            cldmabc(i) = max(cldv(i,k),cldmabc(i))
            cldmabc(i) = cldv(i,k)

            odds=max( &
                 min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8) &
                 *scavcoef(i,k)*deltat),0._r8) 
            srcs2 = sol_factb*cldmabc(i)*odds*tracer(i,k)*(1._r8-weight)/deltat & 
                 +  sol_factbi*cldmabc(i)*odds*tracer(i,k)*(weight)/deltat    



            srcc = srcs1 + srcs2  
            finc = srcs1/(srcc + 1.e-36_r8) 

            
            

            

            
            fracp =  precs(i,k)*deltat/max(cwat(i,k),1.e-12_r8)
            fracp = max(0._r8,min(1._r8,fracp))


            
            
            
            

            srcs1 = sol_facti*cldt(i,k)*fracp*tracer(i,k)/deltat*(1._r8-weight) &  
                 + sol_factii*cldt(i,k)*fracp*tracer(i,k)/deltat*(weight)       


            


            cldmabs(i) = cldv(i,k)   


            odds = precabs(i)/max(cldmabs(i),1.e-5_r8)*scavcoef(i,k)*deltat
            odds = max(min(1._r8,odds),0._r8)
            srcs2 =sol_factb*(cldmabs(i)*odds) *tracer(i,k)*(1._r8-weight)/deltat & 
                 + sol_factbi*(cldmabs(i)*odds) *tracer(i,k)*(weight)/deltat       



            srcs = srcs1 + srcs2             
            fins=srcs1/(srcs + 1.e-36_r8)    

            
            
            rat(i) = tracer(i,k)/max(deltat*(srcc+srcs),1.e-36_r8)
            if (rat(i).lt.1._r8) then
               srcs = srcs*rat(i)
               srcc = srcc*rat(i)
            endif
            srct(i) = (srcc+srcs)*omsm

            
            
            
            fracp = deltat*srct(i)/max(cldmabs(i)*tracer(i,k),1.e-36_r8)  
            fracp = max(0._r8,min(1._r8,fracp))
            fracis(i,k) = 1._r8 - fracp

            
            scavt(i,k) = -srct(i) + fracev(i)*scavab(i)*gravit/pdel(i,k)
            iscavt(i,k) = -(srcc*finc + srcs*fins)*omsm

            if ( present(icscavt) ) icscavt(i,k) = -(srcc*finc) * omsm
            if ( present(isscavt) ) isscavt(i,k) = -(srcs*fins) * omsm
            if ( present(bcscavt) ) bcscavt(i,k) = -(srcc * (1-finc)) * omsm
            if ( present(bsscavt) ) bsscavt(i,k) = -(srcs * (1-fins)) * omsm +  &
                 fracev(i)*scavab(i)*gravit/pdel(i,k)

            dblchek(i) = tracer(i,k) + deltat*scavt(i,k)

            
            scavab(i) = scavab(i)*(1-fracev(i)) + srcs*pdel(i,k)/gravit
            precabs(i) = precabs(i) + (precs(i,k) - evaps(i,k))*pdel(i,k)/gravit
            scavabc(i) = scavabc(i) + srcc*pdel(i,k)/gravit
            precabc(i) = precabc(i) + (cmfdqr(i,k))*pdel(i,k)/gravit
            tracab(i) = tracab(i) + tracer(i,k)*pdel(i,k)/gravit

         end do

         found = .false.
         do i = 1,ncol
            if ( dblchek(i) < 0._r8 ) then
               found = .true.
               exit
            end if
         end do

         if ( found ) then
            do i = 1,ncol
               if (dblchek(i) .lt. 0._r8) then
                  write(iulog,*) ' wetdapa: negative value ', i, k, tracer(i,k), &
                       dblchek(i), scavt(i,k), srct(i), rat(i), fracev(i)
                  call wrf_message(iulog)
               endif
            end do
         endif

      end do

   end subroutine wetdepa_v1


end module wetdep
