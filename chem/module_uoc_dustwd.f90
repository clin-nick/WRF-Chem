MODULE module_uoc_dustwd









USE module_state_description   
USE module_model_constants     
USE physconst                  
USE module_data_gocart_dust         

 CONTAINS

subroutine uoc_dustwd_driver(precr,chem,p_phy,t_phy,                       &
                             ids,ide, jds,jde, kds,kde,                    &
                             ims,ime, jms,jme, kms,kme,                    &
                             its,ite, jts,jte, kts,kte,                    &
                             dtstepc,                                      &
                             dustwd_1, dustwd_2,                           &
                             dustwd_3, dustwd_4,                           &
                             dustwd_5,                                     &
                             wetdep_1, wetdep_2,                           &
                             wetdep_3, wetdep_4,                           &
                             wetdep_5,                                     &
                             dustwdload_1, dustwdload_2,                   &
                             dustwdload_3, dustwdload_4,                   &
                             dustwdload_5,                                 &
                             alt, dz8w, epsilc                             ) 
                             
IMPLICIT NONE

 INTEGER :: debug_level 
 CHARACTER*(100) :: text
 
real :: dustold                           

integer :: d, i, j, k, l                  

integer, parameter :: nbins=5                 
integer, parameter :: nbinsa=5                
integer :: bins(nbinsa)                       


real, dimension(nbins)    :: dbinm             
real, dimension(nbins)    :: dbinmm            


real :: conver9, converi9                 
real :: conver6, converi6                 

real :: wt                                

real :: rmin, rmax                        

integer, parameter        :: nrbins=30    
real, dimension(nrbins)   :: raind        
real, dimension(nrbins+1) :: rend         
real, dimension(nrbins)   :: dsd_rn       
real, dimension(nrbins)   :: vt           
real                      :: z            

real :: visca                             
real :: tw, viscw                         

real :: colece                            

real :: scrate, scavn                     
real :: delR, delv, rnflx, carea          
real :: tair, pair                        
real :: rhoair                            

INTEGER,  INTENT(IN)    ::  ids,ide, jds,jde, kds,kde,    & 
                            ims,ime, jms,jme, kms,kme,    &
                            its,ite, jts,jte, kts,kte
                            
REAL, INTENT(IN)        :: dtstepc, & 
                           epsilc
                            
REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN)    ::                                   precr


REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                &
          INTENT(INOUT) ::                                   chem

 
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                           &
          INTENT(IN)    ::                                   p_phy, t_phy

          
real, dimension( ims:ime, kms:kme, jms:jme ) ::              rnrate


REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) ::                     &
dustwd_1, dustwd_2, dustwd_3, dustwd_4, dustwd_5


REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ) ::           wdbins


REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::                     &
dustwdload_1, dustwdload_2, dustwdload_3, dustwdload_4, dustwdload_5


REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::                     &
wetdep_1, wetdep_2, wetdep_3, wetdep_4, wetdep_5


REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) :: alt, dz8w


den_dust(:) = 2650.

dbinm(:) = 2.0d0*reff_dust(:)

CALL get_wrf_debug_level(debug_level)

rmin =  100.*1.e-6
rmax = 5800.*1.e-6              

 conver6  = 1.e-6               
 converi6 = 1.e6                

 conver9  = 1.e-9               
 converi9 = 1.e9                

bins(1)=p_dust_1                
bins(2)=p_dust_2                
bins(3)=p_dust_3                
bins(4)=p_dust_4                
bins(5)=p_dust_5                

dustwd_1(its:ite,kts:kte,jts:jte) = 0.
dustwd_2(its:ite,kts:kte,jts:jte) = 0.
dustwd_3(its:ite,kts:kte,jts:jte) = 0.
dustwd_4(its:ite,kts:kte,jts:jte) = 0.
dustwd_5(its:ite,kts:kte,jts:jte) = 0.
wetdep_1(its:ite,jts:jte) = 0.
wetdep_2(its:ite,jts:jte) = 0.
wetdep_3(its:ite,jts:jte) = 0.
wetdep_4(its:ite,jts:jte) = 0.
wetdep_5(its:ite,jts:jte) = 0.
dustwdload_1(its:ite,jts:jte) = 0.
dustwdload_2(its:ite,jts:jte) = 0.
dustwdload_3(its:ite,jts:jte) = 0.
dustwdload_4(its:ite,jts:jte) = 0.
dustwdload_5(its:ite,jts:jte) = 0.

wdbins(its:ite,kts:kte,jts:jte,1)=dustwd_1(its:ite,kts:kte,jts:jte)             
wdbins(its:ite,kts:kte,jts:jte,2)=dustwd_2(its:ite,kts:kte,jts:jte)             
wdbins(its:ite,kts:kte,jts:jte,3)=dustwd_3(its:ite,kts:kte,jts:jte)             
wdbins(its:ite,kts:kte,jts:jte,4)=dustwd_4(its:ite,kts:kte,jts:jte)             
wdbins(its:ite,kts:kte,jts:jte,5)=dustwd_5(its:ite,kts:kte,jts:jte)             

rnrate = precr*3600.            

dbinmm = dbinm*conver6          

tw=0                            
call dviscw(tw,viscw)           


do l = 1,nrbins+1
 z = alog10(rmin*converi6)+(alog10(rmax*converi6)-alog10(rmin*converi6))*real(l-1)/real(nrbins)
 rend(l) = real(int(10.**z))*conver6
enddo
do l = 1,nrbins
 z = (alog10(rend(l)*converi6)+alog10(rend(l+1)*converi6))*0.5
 raind(l) = real(int(10.**z))*conver6
enddo






do j = jts, jte
 do k = kts, kte
  do i = its, ite
   if ( precr(i,k,j).gt.0.) then 

    call dsd_rain(rnrate(i,k,j),raind,nrbins,dsd_rn,3) 

    tair=t_phy(i,k,j)
    pair=p_phy(i,k,j)
    rhoair = pair/(rair*tair)     

    call dvisca(tair,visca)

    do l = 1,nrbins
     call fallv(raind(l),tair,pair,rhoair,visca,vt(l))
    enddo

    do d=1,nbins 

     if ( chem(i,k,j,bins(d)).gt.0.) then 

      call w_t(wt,dbinmm(d),den_dust(d),rhoair)

      scrate=0.
      scavn=0.
      do l = 1,nrbins
       
       call coleff(dbinmm(d),raind(l),den_dust(d),pair,tair,wt,vt(l),rhoair,visca,viscw,colece)
       
       delR = (rend(l+1)-rend(l))*1.e2 
       delv = (vt(l)-wt)*1.e2 
       delv = amax1(delv,0.)
       rnflx = dsd_rn(l)*delv
       carea = pi*0.25*(raind(l)*1.e2+dbinmm(d)*1.e2)**2.
       scavn = scavn+carea*colece*rnflx*delR
      enddo
      scrate = scavn

      dustold = chem(i,k,j,bins(d))
      chem(i,k,j,bins(d))=chem(i,k,j,bins(d))-max(0.,chem(i,k,j,bins(d))*scrate*dtstepc)
      chem(i,k,j,bins(d))=max(chem(i,k,j,bins(d)),epsilc)
      wdbins(i,k,j,d)=max(epsilc,dustold-chem(i,k,j,bins(d)))
      
     endif
    enddo






    
   endif

  enddo
 enddo
enddo

do j=jts,jte
 do i=its,ite
  do k=kts,kte
   dustwdload_1(i,j)= max(epsilc,dustwdload_1(i,j) + wdbins(i,k,j,1)/alt(i,k,j) * dz8w(i,k,j))
   dustwdload_2(i,j)= max(epsilc,dustwdload_2(i,j) + wdbins(i,k,j,2)/alt(i,k,j) * dz8w(i,k,j))
   dustwdload_3(i,j)= max(epsilc,dustwdload_3(i,j) + wdbins(i,k,j,3)/alt(i,k,j) * dz8w(i,k,j))
   dustwdload_4(i,j)= max(epsilc,dustwdload_4(i,j) + wdbins(i,k,j,4)/alt(i,k,j) * dz8w(i,k,j))
   dustwdload_5(i,j)= max(epsilc,dustwdload_5(i,j) + wdbins(i,k,j,5)/alt(i,k,j) * dz8w(i,k,j))
   dustwd_1(i,k,j)=max(epsilc,wdbins(i,k,j,1))
   dustwd_2(i,k,j)=max(epsilc,wdbins(i,k,j,2))
   dustwd_3(i,k,j)=max(epsilc,wdbins(i,k,j,3))
   dustwd_4(i,k,j)=max(epsilc,wdbins(i,k,j,4))
   dustwd_5(i,k,j)=max(epsilc,wdbins(i,k,j,5))
  enddo
   wetdep_1(i,j)= max(epsilc,wdbins(i,kts,j,1)/alt(i,kts,j) * dz8w(i,kts,j) / dtstepc)
   wetdep_2(i,j)= max(epsilc,wdbins(i,kts,j,2)/alt(i,kts,j) * dz8w(i,kts,j) / dtstepc)
   wetdep_3(i,j)= max(epsilc,wdbins(i,kts,j,3)/alt(i,kts,j) * dz8w(i,kts,j) / dtstepc)
   wetdep_4(i,j)= max(epsilc,wdbins(i,kts,j,4)/alt(i,kts,j) * dz8w(i,kts,j) / dtstepc)
   wetdep_5(i,j)= max(epsilc,wdbins(i,kts,j,5)/alt(i,kts,j) * dz8w(i,kts,j) / dtstepc)  
 enddo
enddo
                  

end subroutine uoc_dustwd_driver






      subroutine coleff(d,rain,prho,p,t,pv,rv,arho,visca,viscw,ecolec)








      real, intent (in)    :: d                                             
      real, intent (in)    :: rain                                          
      real(8), intent (in) :: prho                                          
      real, intent (in)    :: p                                             
      real, intent (in)    :: t                                             
      real, intent (in)    :: pv                                            
      real, intent (in)    :: rv                                            
      real, intent (in)    :: arho                                          
      real, intent (in)    :: visca, viscw



      real, intent (out) ::  ecolec



      real    :: beta
      real    :: rmin
      real(8) :: E_im,E_in,E_br
      real    :: dd, raind, rhop, pr, ta, vp, vs, rhoa

      CHARACTER*(100) :: text
      integer :: debug_level 

      data rhow/1/                                                       

      CALL get_wrf_debug_level(debug_level)

      dd = d*1.e6       
      raind = rain*1.e6 
      rhop = prho*1.e-3 
      pr = p*1.e-2      
      ta = t-273.15     
      vp = pv*1.e2      
      vs = rv*1.e2      
      rhoa = arho*1.e-3 
      
      
      ecol=0

      if(dd.le.2) then
        call eslinn(raind,dd,pr,ta,rhop,vp,vs,rhoa,visca,viscw,ecol)      
      elseif(dd.gt.2 .and. dd.le.40) then
        call eimpact(rhop,dd,raind,ecol)
      elseif(dd.gt.40) then
        ecol=0
      endif

      ecolec=ecol


      end subroutine coleff









      subroutine eslinn(dr,dp,pair,tair,rhop,vp,vs,rhoa,visca,viscw,eslin)



      real, intent (in) :: dp                                             
      real, intent (in) :: dr                                             
      real, intent (in) :: rhop                                           
      real, intent (in) :: pair                                           
      real, intent (in) :: tair                                           
      real, intent (in) :: vp                                             
      real, intent (in) :: vs                                             
      real, intent (in) :: rhoa                                           
      real, intent (in) :: visca, viscw



      real, intent (out) ::  eslin
      
      real :: pi

      CHARACTER*(100) :: text
      integer :: debug_level 

      real :: gnarf, ren

      E_br=0
      E_in=0
      E_im=0

      drc=dr*1.e-4                                        
      dpc=dp*1.e-4                                        

      pr=pair
      ta=tair

      pi=acos(-1.0)
      xk=1.381e-16                                        




      omg=viscw/visca
      ren=0.5*drc*vs*rhoa/visca                           
      sstar=(1.2+alog(1+ren)/12)/(1+alog(1+ren))          

      xlambda=0.0651*1.e-4 
      cc=1+(2*xlambda/dpc)*(1.257+0.4*exp(-1.1*dpc/(2*xlambda))) 

      tau=dpc**2*rhop*cc/(18*visca)
      stn=2*tau*(vs-vp)/drc                               
      diff=(xk*(ta+273.15)*cc)/(3*pi*visca*dpc)
      scn=visca/(diff*rhoa)                               
      phi=dpc/drc
      pen=ren*scn                                         



      E_br=4*(1+0.4*sqrt(ren)*scn**(1.0/3.0)+0.16*sqrt(ren)*sqrt(scn))/pen



      E_in=4*phi*(1/omg+(1+2*sqrt(ren))*phi)



      if(stn.ge.sstar) then
          E_im=((stn-sstar)/(stn-sstar+2.0/3.0))**1.5
      endif

      eslin=E_br+E_in


      end subroutine eslinn






subroutine eimpact(rhop,dp,raind,e_im)

    USE module_data_uoc_wd   

      integer, parameter :: lp=17, lr=17





      real, intent (in) :: dp                                             
      real, intent (in) :: raind                                          
      real, intent (in) :: rhop                                           



      real, intent (out) ::  e_im

      real  ::  ee(lr)
      real  ::  u(lr+1)
      real  ::  enew(lr+1)
      real  ::  e_r(lp)
      real  ::  e_c(lp)
      real  ::  f(lp)
      real  ::  v(lp+1)
      real  ::  dnew(lp+1)

      CHARACTER*(100) :: text
      integer :: debug_level 
 
      CALL get_wrf_debug_level(debug_level)





      if(raind.lt.100 .or. raind.gt.6000) then
       text = 'raindrop diameter out of range'
       text=trim(trim(text)//" - ERROR - UoC dust wet deposition")       
       call wrf_debug (debug_level,text)
      endif       
      if(dp.lt.2 .or. dp.gt.40) then
       text = 'particle diameter out of range'
       text=trim(trim(text)//" - ERROR - UoC dust wet deposition")       
       call wrf_debug (debug_level,text)
      endif       
      


      rr=raind*0.5
      rp=dp*0.5

      do i=1,lp


        jj=0
        do j=1,lr
          ee(j)=edatawd(i,j)
          if(rr.ge.rainwd(j)) jj=j
        enddo
        if(jj.eq.0) then
         text = 'error in raindrop radius'
         text=trim(trim(text)//" - ERROR - UoC dust wet deposition")       
         call wrf_debug (debug_level,text)
        endif       
        if(rr.gt.rainwd(jj)) then
           lrp1=lr+1
           do l=1,jj
             u(l)=rainwd(l)
           enddo
             u(jj+1)=rr
           do l=jj+1,lr
             u(l+1)=rainwd(l)
           enddo
           call intrpl(lr,log(rainwd),ee,lrp1,log(u),enew)
           e_r(i)=enew(jj+1)
        elseif(rr.eq.rainwd(jj)) then
           e_r(i)=ee(jj)
        endif
      enddo



      rmax=sqrt(1/rhop)

      do i=1,lp
        if(dustwd(i).le.1) then
           f(i)=1
        elseif(dustwd(i).gt.1 .and. dustwd(i).lt.3.98) then
           f(i)=(1-rmax)*(dustwd(i)-3.98)**2/(1-3.98)**2+rmax
        elseif(dustwd(i).ge.3.98) then
           f(i)=rmax
        endif
      enddo

      do i=1,lp
        e_c(i)=e_r(i)*f(i)
      enddo

      ii=0
      e_d=0
      do i=1,lp
        if(rp.ge.dustwd(i)) ii=i
      enddo
      if(ii.eq.0) then
       text = 'wrong particle radius'
       text=trim(trim(text)//" - ERROR - UoC dust wet deposition")       
       call wrf_debug (debug_level,text)
      endif       
      if(rp.gt.dustwd(ii)) then
         do l=1,ii
           v(l)=dustwd(l)
         enddo
           v(ii+1)=rp
         do l=ii+1,lp
           v(l+1)=dustwd(l)
         enddo
           lpp1=lp+1
           call intrpl(lp,dustwd,e_c,lpp1,v,dnew)
           e_d=dnew(ii+1)
      elseif(rp.eq.dustwd(ii)) then
           e_d=e_c(ii)
      endif

      e_im=e_d


      end subroutine eimpact






      subroutine dvisca(tair,dvisc)
      



      integer, parameter :: l=9, n=1

      real, dimension(l) :: ta, visca
      real, dimension(l) :: x, y
      real, dimension(n) :: u, v
     
      data ta/-173,-73,0,20,25,27,127,227,327/
      data visca/0.71e-4,1.33e-4,1.72e-4,1.797e-4,1.818e-4, &
     &           1.86e-4,2.31e-4,2.71e-4,3.08e-4/



      do i = 1,l
         x(i)=ta(i)+273
         y(i)=visca(i)
      enddo

      u(1)=tair

      call intrpl(l,x,y,n,u,v)

      dvisc=v(1)


      end subroutine dvisca







      subroutine dviscw(twt,dvisc)

      integer, parameter :: l=11, n=1

      real, dimension(l) :: tw, viscw
      real, dimension(l) :: x, y
      real, dimension(l) :: u, v

      data tw/0,10,20,30,40,50,60,70,80,90,100/
      data viscw/1.7930e-2,1.3070e-2,1.002e-2,0.7977e-2,0.6532e-2, &
     &           0.5470e-2,0.4665e-2,0.404e-2,0.3544e-2,0.3145e-2, &
     &           0.2818e-2/
     






      do i = 1,l
         x(i)=tw(i)+273.15
         y(i)=viscw(i)
      enddo

      u(1)=twt+273.15

      call intrpl(l,x,y,n,u,v)

      dvisc=v(1)


      end subroutine dviscw







subroutine  intrpl(l,x,y,n,u,v)                                


























      integer, intent(in) :: l, n
      real, dimension(l)  :: x, y
      real, dimension(n)  :: u, v
      equivalence  (p0,x3),(q0,y3),(q1,t3)
      real                :: m1,m2,m3,m4,m5
      equivalence  (uk,dx),(imn,x2,a1,m1),(imx,x5,a5,m5), &
     &             (j,sw,sa),(y2,w2,w4,q2),(y5,w3,q3)
      real    :: a1, a2, a3, a4, a5
      integer :: debug_level 
      CALL get_wrf_debug_level(debug_level)


   10 l0=l
      lm1=l0-1
      lm2=lm1-1
      lp1=l0+1
      n0=n
      if(lm2.lt.0)        go to 90
      if(n0.le.0)         go to 91
      do 11  i=2,l0
        if(x(i-1)-x(i))   11,95,96
   11   continue
      ipv=0

      do 80  k=1,n0
        uk=u(k)

   20   if(lm2.eq.0)      go to 27
        if(uk.ge.x(l0))   go to 26
        if(uk.lt.x(1))    go to 25
        imn=2
        imx=l0
   21   i=(imn+imx)/2
        if(uk.ge.x(i))    go to 23
   22   imx=i
        go to 24
   23   imn=i+1
   24   if(imx.gt.imn)    go to 21
        i=imx
        go to 30
   25   i=1
        go to 30
   26   i=lp1
        go to 30
   27   i=2

   30   if(i.eq.ipv)      go to 70
        ipv=i


   40   j=i
        if(j.eq.1)        j=2
        if(j.eq.lp1)      j=l0
        x3=x(j-1)
        y3=y(j-1)
        x4=x(j)
        y4=y(j)
        a3=x4-x3
        m3=(y4-y3)/a3
        if(lm2.eq.0)      go to 43
        if(j.eq.2)        go to 41
        x2=x(j-2)
        y2=y(j-2)
        a2=x3-x2
        m2=(y3-y2)/a2
        if(j.eq.l0)       go to 42
   41   x5=x(j+1)
        y5=y(j+1)
        a4=x5-x4
        m4=(y5-y4)/a4
        if(j.eq.2)        m2=m3+m3-m4
        go to 45
   42   m4=m3+m3-m2
        go to 45
   43   m2=m3
        m4=m3
   45   if(j.le.3)        go to 46
        a1=x2-x(j-3)
        m1=(y2-y(j-3))/a1
        go to 47
   46   m1=m2+m2-m3
   47   if(j.ge.lm1)      go to 48
        a5=x(j+2)-x5
        m5=(y(j+2)-y5)/a5
        go to 50
   48   m5=m4+m4-m3

   50   if(i.eq.lp1)      go to 52
        w2=abs(m4-m3)
        w3=abs(m2-m1)
        sw=w2+w3
        if(sw.ne.0.0)     go to 51
        w2=0.5
        w3=0.5
        sw=1.0
   51   t3=(w2*m2+w3*m3)/sw
        if(i.eq.1)        go to 54
   52   w3=abs(m5-m4)
        w4=abs(m3-m2)
        sw=w3+w4
        if(sw.ne.0.0)     go to 53
        w3=0.5
        w4=0.5
        sw=1.0
   53   t4=(w3*m3+w4*m4)/sw
        if(i.ne.lp1)      go to 60
        t3=t4
        sa=a2+a3
        t4=0.5*(m4+m5-a2*(a2-a3)*(m2-m3)/(sa*sa))
        x3=x4
        y3=y4
        a3=a2
        m3=m4
        go to 60
   54   t4=t3
        sa=a3+a4
        t3=0.5*(m1+m2-a4*(a3-a4)*(m3-m4)/(sa*sa))
        x3=x3-a4
        y3=y3-m2*a4
        a3=a4
        m3=m2

   60   q2=(2.0*(m3-t3)+m3-t4)/a3
        q3=(-m3-m3+t3+t4)/(a3*a3)

   70   dx=uk-p0
   80   v(k)=q0+dx*(q1+dx*(q2+dx*q3))
      return

   90 call wrf_debug (debug_level,'error 90 in subroutine intrpl - ERROR - UoC dust wet deposition') 
      go to 99
   91 call wrf_debug (debug_level,'error 91 in subroutine intrpl - ERROR - UoC dust wet deposition') 
      go to 99
   95 call wrf_debug (debug_level,'error 95 in subroutine intrpl - ERROR - UoC dust wet deposition') 
      go to 97
   96 call wrf_debug (debug_level,'error 96 in subroutine intrpl - ERROR - UoC dust wet deposition') 
   97 call wrf_debug (debug_level,'error 97 in subroutine intrpl - ERROR - UoC dust wet depositionk') 
   99 call wrf_debug (debug_level,'error 98 in subroutine intrpl - ERROR - UoC dust wet deposition') 
      return








      end subroutine intrpl






      subroutine fallv(d,t,p,rhoair,visc,vt)












real :: gg, pr0, t0, visc0, rhow, rhoa, dp, dl, cl, csc, vt, ta, dropd, pr
real :: a0, a1, a2, a3, a4, a5, a6
real :: b0, b1, b2, b3, b4, b5
real :: c1, c2, dan, x, vy, rey, c3, bon, phn, rhoair, p, t, d

 CHARACTER*(100) :: text
 integer :: debug_level 
 CALL get_wrf_debug_level(debug_level)




        dropd = d*1.e6 
        ta = t-273.15  
        pr = p*1.e-2   

        call sfctens(ta,sig)              

        gg=980                                     
        pr0=1013.25                               
        t0=293.15                                 
        visc0=1.818e-4                            
        rhow=1                                    

        rhoa = rhoair*1.e-3 



        vt=0

        dp=dropd
        dl=dropd*1.e-4                            

        c1=(rhow-rhoa)*gg/(18*visc)
        cl=6.62e-6*(visc/visc0)*(pr0/pr)*sqrt((ta+273.15)/t0)  
        csc=1+2.51*cl/dl

        if(dp.lt.19.) then

         vt=c1*csc*dl**2
         vt=vt*1.e-2 

        elseif(dp.ge.19. .and. dp.lt.1070.) then

         a0=-0.318657e+1
         a1=+0.992696
         a2=-0.153193e-2
         a3=-0.987059e-3
         a4=-0.578878e-3
         a5=+0.855176e-4
         a6=-0.327815e-5
     
         c2=4*rhoa*(rhow-rhoa)*gg/(3*visc*visc)
         dan=c2*dl**3
         x=alog(dan)
         vy=a0+a1*x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6
         rey=csc*exp(vy)
         vt=visc*rey/(rhoa*dl)
         vt=vt*1.e-2 

        elseif(dp.ge.1070. .and. dp.le.7000.) then

         b0=-0.500015e+1
         b1=+0.523778e+1
         b2=-0.204914e+1
         b3=+0.475294
         b4=-0.542819e-1
         b5=+0.238449e-2

         c3=4*(rhow-rhoa)*gg/(3*sig)
         bon=c3*dl*dl                               
         phn=sig**3*rhoa**2/(visc**4*gg*(rhow-rhoa)) 

         y=alog(bon*phn**(1.0/6.0))
         vy=b0+b1*y+b2*y**2+b3*y**3+b4*y**4+b5*y**5
         rey=phn**(1.0/6.0)*exp(vy)
         vt=visc*rey/(rhoa*dl)
         vt=vt*1.e-2 

        endif


        end subroutine fallv






        subroutine sfctens(temp,sig)



        integer, parameter :: n=14
        real, dimension(n) :: t, sfctn
        data t/0,5,10,15,18,20,25,30,40,50,60,70,80,100/
        data sfctn/75.6,74.9,74.22,73.49,73.05,72.75,71.97,71.18,69.56, &
     &             67.91,66.18,64.4,62.6,58.9/
     
        real :: ta

        kk = 1
        xsfc = 0
        ta=max(temp-273.15,0.0)

        do i = 1, 14
           if(i.lt.14) then
              if(ta.ge.t(i) .and. ta.lt.t(i+1)) then
                 kk=i 
                 goto 1
              endif
           else
              if(ta.ge.t(i)) then
                 kk=14
                 goto 1
              endif
           endif
        enddo 

1       continue
       
        if(kk.lt.14) then
             sig = sfctn(kk)+ &
     &       (sfctn(kk+1)-sfctn(kk))*(ta-t(kk))/(t(kk+1)-t(kk))    
        else
             sig = sfctn(kk)+ &
     &       (sfctn(kk)-sfctn(kk-1))*(ta-t(kk))/(t(kk)-t(kk-1))
        endif


        end subroutine sfctens






        subroutine dsd_rain(rnrate,raind,nrbin,dsd_rn,ldsd)


        real, intent(in) :: raind(nrbin), rnrate 
        integer, intent(in) :: ldsd, nrbin

        real :: dsd_mp(nrbin)
        real :: dsd_ss(nrbin)
        real :: dsd_wt(nrbin)
        real :: dsd_fl(nrbin)
        real :: ng
        real :: nt
        real :: pi, dr, wc, d0, dcm, const, dumm

        real, intent(out) :: dsd_rn(nrbin)


        d0=0
        do nk=1,nrbin

          dr=raind(nk)*1.e3                            
       
          if(ldsd.eq.1) then


          dsd_mp(nk)=8.e3*exp(-4.1*dr*rnrate**(-0.21))  
          dsd_mp(nk)=dsd_mp(nk)*1.e-5                   
          dsd_rn(nk)=dsd_mp(nk)

          elseif(ldsd.eq.2) then


          dsd_ss(nk)=0.07*rnrate**0.37*exp(-3.8*dr/rnrate**0.14) 
          dsd_rn(nk)=dsd_ss(nk)

          elseif(ldsd.eq.3) then


          wc=0.062*rnrate**0.913                        
          d0=0.1571*wc**0.1681                          
          ng=512.85*wc*1.e-6/d0**4*d0**(-2.16)
          dcm=dr*0.1                                    
          dsd_wt(nk)=ng*dcm**2.16*exp(-5.5880*dcm/d0)   
          dsd_rn(nk)=dsd_wt(nk)

          elseif(ldsd.eq.4) then
 

          pi=acos(-1.0)
          const=sqrt(2*pi)*alog(1.43)
          nt=172*rnrate**0.22                           
          dg=0.72*rnrate**0.21
          dumm=0.5*alog(dr/dg)**2/alog(1.43)**2
          dsd_fl(nk)=nt/const/dr*exp(-dumm)             
          dsd_fl(nk)=dsd_fl(nk)*1.e-5                   
          dsd_rn(nk)=dsd_fl(nk)

          endif

        enddo

        return
        end subroutine dsd_rain








subroutine w_p(d,rho,wt) 










IMPLICIT NONE 

      REAL, INTENT(OUT) :: wt
      REAL              :: lambda
      REAL              :: rden, err, test, x, xl, xh
      REAL              :: f, fm, fl, fh
      INTEGER           :: i, isign
      real, intent(in)  :: rho, d
      real              :: cc

      REAL, PARAMETER :: visc=1.5E-05
      INTEGER :: debug_level 

      CHARACTER*(100) :: text




       CALL get_wrf_debug_level(debug_level)

       rden=rho 
       lambda=0.0651*1.e-6 
       cc=1+(2*lambda/d)*(1.257+0.4*exp(-1.1*d/(2*lambda))) 

       err=1E-6
       xl=1E-12
       xh=1E3

       wt = xl
       fl = CDSPH(wt*d/visc)*wt**2/cc - 4.0*d*g*rden/3.0

       wt = xh
       fh = CDSPH(wt*d/visc)*wt**2/cc - 4.0*d*g*rden/3.0

       IF (fh*fl.GT.0) call wrf_debug (debug_level,'w_t: F(XL), F(XH) HAVE SAME SIGN - ERROR - UoC dust wet deposition')
       IF (fh-fl.EQ.0) call wrf_debug (debug_level,'w_t: F(XL) = F(XH) - ERROR - UoC dust wet deposition')

       isign=1
       IF (fh.LT.fl) isign=-1

       fh=isign*fh
       fl=isign*fl

       DO 1 i=1,100
         x=(xl+xh)/2
         wt = x
         f = CDSPH(wt*d/visc)*wt**2/cc - 4.0*d*g*rden/3.0
         fm=isign*f
         IF (fm.GT.fh.OR.fm.LT.fl) call wrf_debug (debug_level,'w_t: F(X) NON-MONOTONIC - ERROR - UoC dust wet deposition')
         IF (fm.GT.0) xh=x
         IF (fm.LT.0) xl=x
         IF (fm.EQ.0) then
          return 
         ENDIF
         test=ABS(xh-xl)/ABS(x)
         IF (test.LT.ABS(err)) then
          return
         ENDIF
1      CONTINUE
       call wrf_debug (debug_level,'w_t: NO SOLUTION IN 100 LOOPS - ERROR - UoC dust wet deposition')

END subroutine w_p







SUBROUTINE w_t (wt, dmm, rhop, rhoa)













      REAL(8), INTENT(IN) :: rhop
      REAL,    INTENT(IN) :: rhoa      

      REAL :: wt, dmm
      REAL :: RDEN, D, ERR, TEST, X, XL, XH
      REAL :: F, FM, FL, FH
      REAL, PARAMETER :: VISC=1.5E-05
      INTEGER :: I, ISIGN

      CHARACTER*(100) :: text
      REAL            :: lambda
      real            :: cc

      


      D = dmm
      RDEN=rhop/rhoa 
      lambda=0.0651*1.e-6 
      cc=1+(2*lambda/D)*(1.257+0.4*exp(-1.1*D/(2*lambda))) 

      ERR=1E-6
      XL=1E-12
      XH=1E3
      WT = XL
      FL = CDSPH(WT*D/VISC)*WT**2/cc - 4.0*D*G*RDEN/3.0
      WT = XH
      FH = CDSPH(WT*D/VISC)*WT**2/cc - 4.0*D*G*RDEN/3.0
      IF (FH*FL.GT.0) call wrf_debug (debug_level,'w_t: F(XL), F(XH) HAVE SAME SIGN - ERROR - UoC dust wet deposition')
      IF (FH-FL.EQ.0) call wrf_debug (debug_level,'w_t: F(XL) = F(XH) - ERROR - UoC dust wet deposition')
      ISIGN=1
      IF (FH.LT.FL) ISIGN=-1
      FH=ISIGN*FH
      FL=ISIGN*FL
      DO 1 I=1,100
        X=(XL+XH)/2
        WT = X
        F = CDSPH(WT*D/VISC)*WT**2/cc - 4.0*D*G*RDEN/3.0
        FM=ISIGN*F
        IF (FM.GT.FH.OR.FM.LT.FL) call wrf_debug (debug_level,'w_t: F(X) NON-MONOTONIC - ERROR - UoC dust wet deposition')
        IF (FM.GT.0) XH=X
        IF (FM.LT.0) XL=X
        IF (FM.EQ.0) then
         return 
        ENDIF
        TEST=ABS(XH-XL)/ABS(X)
        IF (TEST.LT.ABS(ERR)) then
         return
        ENDIF
1     CONTINUE
      call wrf_debug (debug_level,'w_t: NO SOLUTION IN 100 LOOPS - ERROR - UoC dust wet deposition')

END subroutine w_t





real FUNCTION CDSPH(RE)














IMPLICIT NONE 

real :: RE

IF (RE.LE.0) THEN
 WRITE(6,*) RE
 STOP 'CDSPH: RE.LE.0'
ELSE IF (RE.LE.0.1) THEN
 CDSPH = 24/RE
ELSE IF (RE.LE.1) THEN
 CDSPH = 22.73/RE + 0.0903/RE**2 + 3.69
ELSE IF (RE.LE.10) THEN
 CDSPH = 29.1667/RE - 3.8889/RE**2 + 1.222
ELSE IF (RE.LE.100) THEN
 CDSPH = 46.5/RE - 116.67/RE**2 + 0.6167
ELSE IF (RE.LE.1000) THEN
 CDSPH = 98.33/RE - 2778/RE**2 + 0.3644
ELSE IF (RE.LE.5000) THEN
 CDSPH = 148.62/RE - 47500/RE**2 + 0.357
ELSE IF (RE.LE.10000) THEN
 CDSPH = -490.546/RE + 578700/RE**2 + 0.46
ELSE IF (RE.LE.50000) THEN
 CDSPH = -1662.5/RE + 5416700/RE**2 + 0.5191
ELSE
 CDSPH = 0.48802
END IF

END FUNCTION CDSPH



END MODULE module_uoc_dustwd
