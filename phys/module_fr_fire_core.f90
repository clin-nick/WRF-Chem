



module module_fr_fire_core

use module_fr_fire_phys, only: fire_params , fire_ros
use module_fr_fire_util








type ignition_line_type
  REAL  ros, &          
        stop_time, &    
        wind_red,  &    
        wrdist,   &     
        wrupwind, &     
        start_x, &      
        start_y, &      
        end_x, &        
        end_y, &        
        start_time, &   
        end_time, &     
        radius          
end type ignition_line_type

contains




    
subroutine init_no_fire(&
    ifds,ifde,jfds,jfde, &
    ifms,ifme,jfms,jfme, &
    ifts,ifte,jfts,jfte, &
    fdx,fdy,time_now,    & 
    fuel_frac,fire_area,lfn,tign)    
implicit none
             



integer, intent(in):: ifds,ifde,jfds,jfde   
integer, intent(in):: ifts,ifte,jfts,jfte   
integer, intent(in):: ifms,ifme,jfms,jfme   
real, intent(in) :: fdx,fdy,time_now        
real, intent(out), dimension (ifms:ifme,jfms:jfme) :: & 
                   fuel_frac,fire_area,lfn,tign       


intrinsic epsilon
                                                

integer:: i,j
real lfn_init,time_init

lfn_init = 2*max((ifde-ifds+1)*fdx,(jfde-jfds+1)*fdy)      
time_init=time_now + max(time_now,1.0)*epsilon(time_now) 
 
do j=jfts,jfte
    do i=ifts,ifte
        fuel_frac(i,j)=1.          
        fire_area(i,j)=0.          
        tign(i,j) = time_init      
        lfn(i,j) = lfn_init        
    enddo
enddo
call message('init_no_fire: state set to no fire')

end subroutine init_no_fire




 

subroutine ignite_fire( ifds,ifde,jfds,jfde,                    & 
                        ifms,ifme,jfms,jfme,                      &
                        ifts,ifte,jfts,jfte,                      &
                        ignition_line,                            &
                        start_ts,end_ts,                    &
                        coord_xf,coord_yf,                &     
                        unit_xf,unit_yf,                  &
                        lfn,tign,ignited)
implicit none














integer, intent(in):: ifds,ifde,jfds,jfde   
integer, intent(in):: ifts,ifte,jfts,jfte   
integer, intent(in):: ifms,ifme,jfms,jfme   
type(ignition_line_type), intent(in):: ignition_line    
real, intent(in):: start_ts,end_ts          
real, dimension(ifms:ifme, jfms:jfme), intent(in):: & 
    coord_xf,coord_yf                       
real, intent(in):: unit_xf,unit_yf          
real, intent(inout), dimension (ifms:ifme,jfms:jfme) :: & 
                   lfn, tign                
integer, intent(out):: ignited              
                        

integer:: i,j
real::lfn_new,time_ign,ax,ay,rels,rele,d
real:: sx,sy                    
real:: ex,ey                    
real:: st,et                    
character(len=128):: msg
real::cx2,cy2,dmax,axmin,axmax,aymin,aymax,dmin
real:: start_x,start_y          
real:: end_x,end_y              
real:: radius                   
real:: start_time,end_time      
real:: ros,tos                  




start_x    = ignition_line%start_x 
start_y    = ignition_line%start_y 
end_x      = ignition_line%end_x   
end_y      = ignition_line%end_y   
start_time = ignition_line%start_time 
end_time   = ignition_line%end_time
radius     = ignition_line%radius  
ros        = ignition_line%ros     
tos        = radius/ros            
st         = start_time            
et         = min(end_ts,end_time)  


if(start_ts>et+tos .or. end_ts<st)return   

if(start_time < end_time)then  
        
        
	
	

        sx = start_x
        sy = start_y
        rele =  (et - start_time) / (end_time - start_time)    
	ex = start_x + rele * (end_x - start_x)
	ey = start_y + rele * (end_y - start_y)
else
        


	sx = start_x
	sy = start_y
	ex = end_x
	ey = end_y
endif


cx2=unit_xf*unit_xf
cy2=unit_yf*unit_yf

axmin=coord_xf(ifts,jfts)
aymin=coord_yf(ifts,jfts)
axmax=coord_xf(ifte,jfte)
aymax=coord_yf(ifte,jfte)
!$OMP CRITICAL(FIRE_CORE_CRIT)
write(msg,'(a,2f11.6,a,2f11.6)')'IGN from ',sx,sy,' to ',ex,ey
call message(msg)
write(msg,'(a,2f10.2,a,2f10.2,a)')'IGN timestep [',start_ts,end_ts,'] in [',start_time,end_time,']'
call message(msg)
write(msg,'(a,2g13.6,a,2g13.6)')'IGN tile coord from  ',axmin,aymin,' to ',axmax,aymax
call message(msg)
!$OMP END CRITICAL(FIRE_CORE_CRIT)
ignited=0
dmax=0
dmin=huge(dmax)
11      format('IGN ',6(a,g17.7,1x)) 
12      format('IGN ',4(a,2g13.7,1x))
do j=jfts,jfte   
    do i=ifts,ifte
        ax=coord_xf(i,j)
        ay=coord_yf(i,j)

        
        
        call nearest(d,time_ign,ax,ay,sx,sy,st,ex,ey,et,cx2,cy2)
        dmax=max(d,dmax)
        dmin=min(d,dmin)

        lfn_new=d - min( radius, ros*(end_ts - time_ign) )  
        if(fire_print_msg.ge.3)then
!$OMP CRITICAL(FIRE_CORE_CRIT)
            write(msg,*)'IGN1 i,j=',i,j,' lfn(i,j)=',lfn(i,j),' tign(i,j)=',tign(i,j)
            call message(msg)
            write(msg,*)'IGN2 i,j=',i,j,' lfn_new= ',lfn_new, ' time_ign= ',time_ign,' d=',d
            call message(msg)
!$OMP END CRITICAL(FIRE_CORE_CRIT)
        endif
        if(.not.lfn_new>0.) then
            ignited=ignited+1   
        endif
        if(lfn(i,j)>0. .and. .not. lfn_new > 0.) then 
            tign(i,j)=time_ign + d/ros  
            if(fire_print_msg.ge.3)then
!$OMP CRITICAL(FIRE_CORE_CRIT)
                write(msg,'(a,2i6,a,2g13.6,a,f10.2,a,2f10.2,a)')'IGN ignited cell ',i,j,' at',ax,ay, &
                    ' time',tign(i,j),' in [',start_ts,end_ts,']'
                call message(msg)
                write(msg,'(a,g10.3,a,f10.2,a,2f10.2,a)')'IGN distance',d,' from ignition line at',time_ign,' in [',st,et,']'
                call message(msg)
!$OMP END CRITICAL(FIRE_CORE_CRIT)
            endif
            if(tign(i,j) < start_ts .or. tign(i,j) > end_ts )then
!$OMP CRITICAL(FIRE_CORE_CRIT)
                write(msg,'(a,2i6,a,f11.6,a,2f11.6,a)')'WARNING ',i,j, &
                ' fixing ignition time ',tign(i,j),' outside of the time step [',start_ts,end_ts,']'
                call message (msg)
!$OMP END CRITICAL(FIRE_CORE_CRIT)
                tign(i,j) = min(max(tign(i,j),start_ts),end_ts)
            endif
        endif
        lfn(i,j)=min(lfn(i,j),lfn_new)  
        if(fire_print_msg.ge.3)then
!$OMP CRITICAL(FIRE_CORE_CRIT)
            write(msg,*)'IGN3 i,j=',i,j,' lfn(i,j)=',lfn(i,j),' tign(i,j)=',tign(i,j)
            call message(msg)
!$OMP END CRITICAL(FIRE_CORE_CRIT)
        endif
    enddo
enddo
!$OMP CRITICAL(FIRE_CORE_CRIT)
write(msg,'(a,2g13.2,a,g10.2,a,g10.2)')'IGN units ',unit_xf,unit_yf,' m max dist ',dmax,' min',dmin
call message(msg)
write(msg,'(a,f6.1,a,f8.1,a,i10)')'IGN radius ',radius,' time of spread',tos,' ignited nodes',ignited
call message(msg)
!$OMP END CRITICAL(FIRE_CORE_CRIT)

return
99 continue

end subroutine ignite_fire


!DEC$ ATTRIBUTES FORCEINLINE
SUBROUTINE nearest(d,t,ax,ay,sx,sy,st,ex,ey,et,cx2,cy2)
        implicit none

        real, intent(out):: d,t
        real, intent(in):: ax,ay,sx,sy,st,ex,ey,et,cx2,cy2
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        real:: mx,my,dam2,dames,am_es,cos2,dmc2,mcrel,mid_t,dif_t,des2,cx,cy
        character(len=128):: msg


11      format('IGN ',6(a,g17.7,1x))
12      format('IGN ',4(a,2g13.7,1x))

        
        mx = (sx + ex)*0.5
        my = (sy + ey)*0.5
        dam2=(ax-mx)*(ax-mx)*cx2+(ay-my)*(ay-my)*cy2      
        des2 = (ex-sx)*(ex-sx)*cx2+(ey-sy)*(ey-sy)*cy2          
        dames = dam2*des2
        am_es=(ax-mx)*(ex-sx)*cx2+(ay-my)*(ey-sy)*cy2       
        if(dames>0)then
            cos2 = (am_es*am_es)/dames                  
        else 
            cos2 = 0.
        endif
        dmc2 = dam2*cos2                                
        if(4.*dmc2 < des2)then                          
            
            mcrel = sign(sqrt(4.*dmc2/des2),am_es)      
        elseif(am_es>0)then                             
            mcrel = 1.0 
        else                                            
            mcrel = -1.0
        endif
	cx = (ex + sx)*0.5 + mcrel*(ex - sx)*0.5     
	cy = (ey + sy)*0.5 + mcrel*(ey - sy)*0.5     
        d=sqrt((ax-cx)*(ax-cx)*cx2+(ay-cy)*(ay-cy)*cy2) 
	t = (et + st)*0.5 + mcrel*(et - st)*0.5     
        if(fire_print_msg.ge.3)then
!$OMP CRITICAL(FIRE_CORE_CRIT)
            write(msg,12)'find nearest to [',ax,ay,'] from [',sx,sy,'] [',ex,ey,']' 
            call message(msg)
            write(msg,12)'end times',st,et,' scale squared',cx2,cy2 
            call message(msg)
            write(msg,11)'nearest at mcrel=',mcrel,'from the midpoint, t=',t 
            call message(msg)
            write(msg,12)'nearest is [',cx,cy,'] d=',d 
            call message(msg)
            write(msg,11)'dam2=',dam2,'des2=',des2,'dames=',dames
            call message(msg)
            write(msg,11)'am_es=',am_es,'cos2=',cos2,'dmc2=',dmc2 
            call message(msg)
!$OMP END CRITICAL(FIRE_CORE_CRIT)
        endif
END SUBROUTINE nearest






subroutine fuel_left(&
    ims,ime,jms,jme, &
    its,ite,jts,jte, &
    ifs,ife,jfs,jfe, &
    lfn, tign, fuel_time, time_now, fuel_frac, fire_area)
implicit none






integer, intent(in) :: its,ite,jts,jte,ims,ime,jms,jme,ifs,ife,jfs,jfe
real, intent(in), dimension(ims:ime,jms:jme)::lfn,tign,fuel_time
real, intent(in):: time_now
real, intent(out), dimension(ifs:ife,jfs:jfe)::fuel_frac
real, intent(out), dimension(ims:ime,jms:jme):: fire_area













integer::i,j,ir,jr,icl,jcl,isubcl,jsubcl,i2,j2,ii,jj
real::fmax,frat,helpsum1,helpsum2,fuel_left_ff,fire_area_ff,rx,ry,tignf(2,2)

real::lffij,lffi1j,lffij1,lffi1j1,tifij,tifi1j,tifij1,tifi1j1,tx,ty,txx,tyy

character(len=128)::msg
integer::m,omp_get_thread_num
     


ir=fuel_left_irl
jr=fuel_left_jrl

if ((ir.ne.2).or.(jr.ne.2)) then 
   call crash('fuel_left: ir.ne.2 or jr.ne.2 ')
endif

rx=1./ir 
ry=1./jr
















do icl=its,ite
  do jcl=jts,jte
    helpsum1=0
    helpsum2=0

    do isubcl=1,ir
      do jsubcl=1,jr 
        i=(icl-its)*ir+isubcl
        j=(jcl-jts)*jr+jsubcl


        if ((isubcl.eq.1).and.(jsubcl.eq.1)) then
           i2=icl-1
           j2=jcl-1
           ty=0.5
           tx=0.5
           tyy=1.0
           txx=1.0
        else if ((isubcl.eq.2).and.(jsubcl.eq.1)) then
           i2=icl
           j2=jcl-1
           ty=0.5
           tx=0
           tyy=1.0
           txx=0.5
        else if ((isubcl.eq.1).and.(jsubcl.eq.2)) then
           i2=icl-1
           j2=jcl
           tx=0.5
           ty=0
           txx=1.0
           tyy=0.5
        else if ((isubcl.eq.2).and.(jsubcl.eq.2)) then
           i2=icl
           j2=jcl
           tx=0
           ty=0
           txx=0.5
           tyy=0.5
        else
           call crash('fuel_left: isubcl,jsubcl should be only 1 or 2')
        endif 


        lffij=                             &    
                  (1-tx)*(1-ty)*lfn(i2,j2)      &
             +    (1-tx)*ty  *lfn(i2,j2+1)      &
             +     tx*(1-ty)*lfn(i2+1,j2)       &
             +       tx*ty  *lfn(i2+1,j2+1)
        lffi1j=                            &
                    (1-txx)*(1-ty)*lfn(i2,j2)   &
             +      (1-txx)*ty  *lfn(i2,j2+1)   &
             +      (txx)*(1-ty)*lfn(i2+1,j2)   &
             +      (txx)*ty  *lfn(i2+1,j2+1)
        lffij1=                            &
                    (1-tx)*(1-tyy)*lfn(i2,j2)   &
             +      (1-tx)*(tyy)  *lfn(i2,j2+1) &
             +      tx*(1-tyy)*lfn(i2+1,j2)     &
             +      tx*(tyy)  *lfn(i2+1,j2+1)
        lffi1j1 =                               &
                      (1-txx)*(1-tyy)*lfn(i2,j2)     &
             +      (1-txx)*(tyy)  *lfn(i2,j2+1)   &        
             +      (txx)*(1-tyy)*lfn(i2+1,j2)     &
             +      (txx)*(tyy)  *lfn(i2+1,j2+1)

        
        do ii=1,2
          do jj=1,2
            tignf(ii,jj)=tign(i2+ii-1,j2+jj-1)
          enddo
        enddo
        tifij=                                 &
                   (1-tx)*(1-ty)*tignf(1,1)        &
             +     (1-tx)*ty*tignf(1,1+1)          &
             +     tx*(1-ty)*tignf(1+1,1)          &
             +     tx*ty*tignf(1+1,1+1)
        tifi1j=                               &
                   (1-txx)*(1-ty)*tignf(1,1)      &
             +     (1-txx)*ty*tignf(1,1+1)        &
             +     (txx)*(1-ty)*tignf(1+1,1)      &
             +     (txx)*(ty)*tignf(1+1,1+1)            
        tifij1=                               &
                   (1-tx)*(1-tyy)*tignf(1,1)      &
             +     (1-tx)*(tyy)*tignf(1,1+1)      &
             +      tx*(1-tyy)*tignf(1+1,1)       &
             +      tx*(tyy)*tignf(1+1,1+1)
        tifi1j1=                               &
                   (1-txx)*(1-tyy)*tignf(1,1)     &
             +     (1-txx)*(tyy)*tignf(1,1+1)     &
             +     (txx)*(1-tyy)*tignf(1+1,1)     &
             +     (txx)*(tyy)*tignf(1+1,1+1) 

         
        if(fuel_left_method.eq.1)then
          call fuel_left_cell_1( fuel_left_ff, fire_area_ff, &
             lffij,lffij1,lffi1j,lffi1j1,&
             tifij,tifij1,tifi1j,tifi1j1,&
             time_now, fuel_time(icl,jcl))
        elseif(fuel_left_method.eq.2)then
          fire_area_ff=0  
          fuel_left_ff=fuel_left_cell_2( &
             lffij,lffij1,lffi1j,lffi1j1,&
             tifij,tifij1,tifi1j,tifi1j1,&
             time_now, fuel_time(icl,jcl)) 
        else
          call crash('fuel_left: unknown fuel_left_method')
        endif

        
        if(fire_area_ff.lt.-1e-6 .or.  &
          (fire_area_ff.eq.0. .and. fuel_left_ff.lt.1.-1e-6))then
!$OMP CRITICAL(FIRE_CORE_CRIT)
           write(msg,'(a,2i6,2(a,f11.8))')'fuel_left: at node',i,j, &
              ' of refined mesh fuel burnt',1-fuel_left_ff,' fire area',fire_area_ff
!$OMP END CRITICAL(FIRE_CORE_CRIT)
           call crash(msg)
        endif

        helpsum1=helpsum1+fuel_left_ff
        helpsum2=helpsum2+fire_area_ff
      enddo
    enddo
    fuel_frac(icl,jcl)=helpsum1 
    fire_area(icl,jcl)=helpsum2
  enddo 
enddo
  





do j=jts,jte
    do i=its,ite        
        fuel_frac(i,j) = fuel_frac(i,j) /(ir*jr) 
        fire_area(i,j) = fire_area(i,j) /(ir*jr) 
    enddo
enddo


fmax=0
do j=jts,jte
    do i=its,ite        
       if(fire_area(i,j).eq.0.)then
           if(fuel_frac(i,j).lt.1.-1e-6)then
!$OMP CRITICAL(FIRE_CORE_CRIT)
               write(msg,'(a,2i6,2(a,f11.8))')'fuel_left: at node',i,j, &
                   ' fuel burnt',1-fuel_frac(i,j),' but fire area',fire_area(i,j)
!$OMP END CRITICAL(FIRE_CORE_CRIT)
               call crash(msg)
           endif
       else
           frat=(1-fuel_frac(i,j))/fire_area(i,j)
           fmax=max(fmax,frat)
       endif
    enddo
enddo
!$OMP CRITICAL(FIRE_CORE_CRIT)
write(msg,'(a,4i6,a,f10.7)')'fuel_left: tile',its,ite,jts,jte,' max fuel burnt/area',fmax 
!$OMP END CRITICAL(FIRE_CORE_CRIT)
call message(msg)
return


end subroutine fuel_left





subroutine fuel_left_cell_1( fuel_frac_left, fire_frac_area, &
    lfn00,lfn01,lfn10,lfn11, &
    tign00,tign01,tign10,tign11,&
    time_now, fuel_time_cell)

implicit none

real, intent(out):: fuel_frac_left, fire_frac_area 
real, intent(in)::lfn00,lfn01,lfn10,lfn11    
real, intent(in)::tign00,tign01,tign10,tign11
real, intent(in)::time_now                   
real, intent(in)::fuel_time_cell            









































































intrinsic tiny


real::ps,aps,area,ta,out
real::t00,t01,t10,t11
real,parameter::safe=tiny(aps)
character(len=128)::msg







t00=tign00-time_now
if(lfn00>0. .or. t00>0.)t00=0.
t01=tign01-time_now
if(lfn01>0. .or. t01>0.)t01=0.
t10=tign10-time_now
if(lfn10>0. .or. t10>0.)t10=0.
t11=tign11-time_now
if(lfn11>0. .or. t11>0.)t11=0.


ps = lfn00+lfn01+lfn10+lfn11   
aps = abs(lfn00)+abs(lfn01)+abs(lfn10)+abs(lfn11)
aps=max(aps,safe)
area =(-ps/aps+1.)/2.
area = max(area,0.) 
area = min(area,1.)
    

ta=0.25*(t00+t01+t10+t11)


out=1.

if(area>0)out=area*exp(ta/fuel_time_cell) + (1. - area)

if(out>1.)then
!$OMP CRITICAL(FIRE_CORE_CRIT)
    write(msg,*)'out=',out,'>1 area=',area,' ta=',ta
    call message(msg)
    write(msg,*)'tign=',tign00,tign01,tign10,tign11,' time_now=',time_now
!$OMP END CRITICAL(FIRE_CORE_CRIT)
    call message(msg)

    call crash('fuel_left_cell_1: fuel fraction > 1')
endif




fuel_frac_left = out
fire_frac_area = area

end subroutine fuel_left_cell_1





real function fuel_left_cell_2(  &
    lfn00,lfn01,lfn10,lfn11, &
    tign00,tign01,tign10,tign11,&
    time_now, fuel_time_cell)

implicit none

real, intent(in)::lfn00,lfn01,lfn10,lfn11    
real, intent(in)::tign00,tign01,tign10,tign11
real, intent(in)::time_now                   
real, intent(in)::fuel_time_cell            








































































call crash('fuel_left_cell_2: not implemented, please use fire_fuel_left_method=1')
fuel_left_cell_2=0.  
end function fuel_left_cell_2

subroutine prop_ls_rk3(id,               &             
                ifds,ifde,jfds,jfde,     &                 
                ifms,ifme,jfms,jfme,     &                    
                ifps,ifpe,jfps,jfpe,     &            
                ifts,ifte,jfts,jfte,     &                    
                ts,dt,dx,dy,             &                     
                tbound,                  &                  
                lfn_in,                  &
                lfn_0,lfn_1,lfn_2,       & 
                lfn_out,tign,ros,        &           
                fp,                      &
                grid,                    &       
                ids,ide,jds,jde,kds,kde, & 
                ims,ime,jms,jme,kms,kme, & 
                ips,ipe,jps,jpe,kps,kpe  &
                   )


    USE module_dm        , ONLY : ntasks_x,ntasks_y,local_communicator,mytask,ntasks
    USE module_comm_dm , ONLY : halo_fire_lfn_1_sub, halo_fire_lfn_2_sub, halo_fire_lfn_0_sub
USE module_domain , only: domain

implicit none









































  
















type(domain) , target :: grid                     
integer, intent(in):: ids,ide,jds,jde,kds,kde, &
    ims,ime,jms,jme,kms,kme, &
    ips,ipe,jps,jpe,kps,kpe


integer,intent(in)::id,ifms,ifme,jfms,jfme,ifds,ifde,jfds,jfde,ifts,ifte,jfts,jfte,ifps,ifpe,jfps,jfpe 
real,dimension(ifms:ifme,jfms:jfme),intent(inout)::lfn_in,tign
real,dimension(ifms:ifme,jfms:jfme),intent(inout)::lfn_1,lfn_2,lfn_0 
real,dimension(ifms:ifme,jfms:jfme),intent(out)::lfn_out,ros
real,intent(in)::dx,dy,ts,dt
real,intent(out)::tbound
type(fire_params),intent(in)::fp

real,dimension(ifms:ifme,jfms:jfme):: tend

real::grad2,rr,tbound2,tbound3

real::gradx,grady,aspeed,err,aerr,time_now
integer::i,j,its1,ite1,jts1,jte1,k,kk,id1
character(len=128)::msg
integer::nfirenodes,nfireline
real::sum_err,min_err,max_err,sum_aerr,min_aerr,max_aerr   


integer,parameter :: mstep=1000, printl=1
real, parameter:: zero=0.,one=1.,eps=epsilon(zero),tol=100*eps, &
    safe=2.,rmin=safe*tiny(zero),rmax=huge(zero)/safe



intrinsic max,min,sqrt,nint,epsilon,tiny,huge



!$OMP CRITICAL(FIRE_CORE_CRIT)
write(msg,'(a8,i5,a6,i5,3(a1,i5))')'prop_ls:',id,' tile ',ifts,':',ifte,',',jfts,':',jfte
!$OMP END CRITICAL(FIRE_CORE_CRIT)
call message(msg)





    do j=jfts,jfte
        do i=ifts,ifte
            lfn_0(i,j) = lfn_in(i,j)
        enddo
    enddo







CALL HALO_FIRE_LFN_0_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


    id1 = id  
    if(id1.ne.0)id1=id1+1000
    call tend_ls(id1,    &
    ifds,ifde,jfds,jfde, &                       
    ifps,ifpe,jfps,jfpe, &                       
    ifts,ifte,jfts,jfte, &                       
    ifms,ifme,jfms,jfme, &                       
    ts,dt,dx,dy,         &                       
    lfn_0,               &                       
    tbound,              &                       
    tend, ros,           &                       
    fp                   &                       
)

    do j=jfts,jfte 
        do i=ifts,ifte 
            lfn_1(i,j) = lfn_0(i,j) + (dt/3.0)*tend(i,j) 
        enddo
    enddo







CALL HALO_FIRE_LFN_1_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )






    if(id1.ne.0)id1=id1+1000
    call tend_ls(id1,    &
    ifds,ifde,jfds,jfde, &                     
    ifps,ifpe,jfps,jfpe, &                    
    ifts,ifte,jfts,jfte, &                    
    ifms,ifme,jfms,jfme, &               
    ts+dt,dt,dx,dy,      &                  
    lfn_1,               &                             
    tbound2,             &                              
    tend,ros,            &                           
    fp                   &
)

    do j=jfts,jfte
        do i=ifts,ifte
            lfn_2(i,j) = lfn_0(i,j) + (dt/2.0)*tend(i,j)
        enddo
    enddo     







CALL HALO_FIRE_LFN_2_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )






    if(id1.ne.0)id1=id1+1000
    call tend_ls(id1,    &
    ifds,ifde,jfds,jfde, &                      
    ifps,ifpe,jfps,jfpe, &                     
    ifts,ifte,jfts,jfte, &                    
    ifms,ifme,jfms,jfme, &                 
    ts+dt,dt,dx,dy,      &                  
    lfn_2,               &                            
    tbound3,             &                            
    tend,ros,            &                           
    fp                   &
)

    do j=jfts,jfte
        do i=ifts,ifte
            lfn_out(i,j) = lfn_0(i,j) + dt*tend(i,j)
            lfn_2(i,j) = lfn_out(i,j) 
        enddo
    enddo     







CALL HALO_FIRE_LFN_2_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )




    tbound=min(tbound,tbound2,tbound3)

!$OMP CRITICAL(FIRE_CORE_CRIT)
    write(msg,'(a,f10.2,4(a,f7.2))')'prop_ls: time',ts,' dt=',dt,' bound',min(tbound,999.99), &
        ' dx=',dx,' dy=',dy
!$OMP END CRITICAL(FIRE_CORE_CRIT)
    call message(msg)
    if(dt>tbound)then
!$OMP CRITICAL(FIRE_CORE_CRIT)
        write(msg,'(2(a,f10.2))')'prop_ls: WARNING: time step ',dt, &
        ' > bound =',tbound
!$OMP END CRITICAL(FIRE_CORE_CRIT)
        call message(msg)
    endif
    
end subroutine prop_ls_rk3





subroutine reinit_ls_rk3(id,                       &
                ifts,ifte,jfts,jfte,               &                     
                ifms,ifme,jfms,jfme,               &                     
                ifds,ifde,jfds,jfde,               &                     
                ifps,ifpe,jfps,jfpe,               &                      
                ts,dt,dx,dy,                       &                    
                lfn_in,                            & 
                lfn_2,lfn_s0,lfn_s1,lfn_s2,lfn_s3, & 
                lfn_out,tign,                      &             
                grid,                              &          
                ids,ide,jds,jde,kds,kde,           &
                ims,ime,jms,jme,kms,kme,           & 
                ips,ipe,jps,jpe,kps,kpe            &
                ) 


    USE module_dm        , ONLY : ntasks_x,ntasks_y,local_communicator,mytask,ntasks
    USE module_comm_dm , ONLY : halo_fire_lfn_s1_sub,halo_fire_lfn_s2_sub,halo_fire_lfn_s3_sub
USE module_domain , only: domain
USE module_configure, only: grid_config_rec_type

implicit none













type(domain) , target :: grid
integer, intent(in):: ids,ide,jds,jde,kds,kde, & 
                      ims,ime,jms,jme,kms,kme, & 
                      ips,ipe,jps,jpe,kps,kpe    

integer,intent(in)::ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,id 
integer,intent(in)::ifds,ifde,jfds,jfde,ifps,ifpe,jfps,jfpe  
real,dimension(ifms:ifme,jfms:jfme),intent(inout)::lfn_in,tign
real,dimension(ifms:ifme,jfms:jfme),intent(inout)::lfn_2,lfn_s0,lfn_s1,lfn_s2,lfn_s3
real,dimension(ifms:ifme,jfms:jfme),intent(inout)::lfn_out
real,intent(in)::dx,dy,ts,dt
real::dt_s
real,dimension(ifts:ifte,jfts:jfte):: tend_1,tend_2,tend_3
real::diffLx,diffLy,diffRx,diffRy,diff2x,diff2y,grad,time_now
integer::nts,i,j,k,kk
intrinsic epsilon
character(len=128)::msg
integer::itso,iteo,jtso,jteo
real::threshold_HLl,threshold_HLu
integer,parameter::bdy_eno1=10


threshold_HLl=-fire_lsm_band_ngp*dx 
threshold_HLu=fire_lsm_band_ngp*dx



    do j=jfts,jfte 
        do i=ifts,ifte 
            lfn_s0(i,j) = lfn_2(i,j)/sqrt(lfn_2(i,j)**2.0+dx**2.0) 
            lfn_s3(i,j) = lfn_2(i,j)
        enddo
    enddo






CALL HALO_FIRE_LFN_S3_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


    call continue_at_boundary(1,1,fire_lfn_ext_up, &   
    ifms,ifme,jfms,jfme, &                             
    ifds,ifde,jfds,jfde, &                             
    ifps,ifpe,jfps,jfpe, &                             
    ifts,ifte,jfts,jfte, &                             
    itso,iteo,jtso,jteo, &                             
    lfn_s3) 

    
dt_s=0.01*dx
do nts=1,fire_lsm_reinit_iter 
    
    
    call advance_ls_reinit( &
    ifms,ifme,jfms,jfme,    &           
    ifds,ifde,jfds,jfde,    &           
    ifts,ifte,jfts,jfte,    &               
    dx,dy,dt_s,bdy_eno1,threshold_HLu, &               
    lfn_s0,lfn_s3,lfn_s3,lfn_s1,1.0/3.0)   







CALL HALO_FIRE_LFN_S1_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


    call continue_at_boundary(1,1,fire_lfn_ext_up, & 
    ifms,ifme,jfms,jfme, &                          
    ifds,ifde,jfds,jfde, &                           
    ifps,ifpe,jfps,jfpe, &                           
    ifts,ifte,jfts,jfte, &                           
    itso,iteo,jtso,jteo, &               
    lfn_s1) 

    

    call advance_ls_reinit( &
    ifms,ifme,jfms,jfme,    &
    ifds,ifde,jfds,jfde,    &
    ifts,ifte,jfts,jfte,    &
    dx,dy,dt_s,bdy_eno1,threshold_HLu, &
    lfn_s0,lfn_s3,lfn_s1,lfn_s2,1.0/2.0)







CALL HALO_FIRE_LFN_S2_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


    call continue_at_boundary(1,1,fire_lfn_ext_up, &  
    ifms,ifme,jfms,jfme, &              
    ifds,ifde,jfds,jfde, &                 
    ifps,ifpe,jfps,jfpe, &                 
    ifts,ifte,jfts,jfte, &               
    itso,iteo,jtso,jteo, & 
    lfn_s2) 

    

    call advance_ls_reinit( &
    ifms,ifme,jfms,jfme,    &
    ifds,ifde,jfds,jfde,    &
    ifts,ifte,jfts,jfte,    &
    dx,dy,dt_s,bdy_eno1,threshold_HLu, &
    lfn_s0,lfn_s3,lfn_s2,lfn_s3,1.0)  







CALL HALO_FIRE_LFN_S3_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


    call continue_at_boundary(1,1,fire_lfn_ext_up, &  
    ifms,ifme,jfms,jfme, &              
    ifds,ifde,jfds,jfde, &                  
    ifps,ifpe,jfps,jfpe, &                
    ifts,ifte,jfts,jfte, &              
    itso,iteo,jtso,jteo, &           
    lfn_s3) 

enddo 

    do j=jfts,jfte 
        do i=ifts,ifte 
          lfn_out(i,j)=lfn_s3(i,j) 
          lfn_out(i,j)=min(lfn_out(i,j),lfn_in(i,j)) 
        enddo
    enddo

end subroutine reinit_ls_rk3





subroutine advance_ls_reinit(ifms,ifme,jfms,jfme, &
                          ifds,ifde,jfds,jfde,    &
                          ifts,ifte,jfts,jfte,    &
                          dx,dy,dt_s,bdy_eno1,threshold_HLu, &
                          lfn_s0,lfn_ini,lfn_curr,lfn_fin,rk_coeff)









implicit none

integer,intent(in)::ifms,ifme,jfms,jfme,ifts,ifte,jfts,jfte,ifds,ifde,jfds,jfde,bdy_eno1
real,dimension(ifms:ifme,jfms:jfme),intent(in)::lfn_s0,lfn_ini,lfn_curr
real,dimension(ifms:ifme,jfms:jfme),intent(inout)::lfn_fin
real,intent(in)::dx,dy,dt_s,threshold_HLu,rk_coeff

integer::i,j
real::diffLx,diffLy,diffRx,diffRy,diff2x,diff2y,grad,tend_r


    do j=jfts,jfte 
        do i=ifts,ifte 
          if (i.lt.ifds+bdy_eno1 .OR. i.gt.ifde-bdy_eno1 .OR. j.lt.jfds+bdy_eno1 .OR. j.gt.jfde-bdy_eno1) then 
            diffLx=(lfn_curr(i,j)-lfn_curr(i-1,j))/dx
            diffLy=(lfn_curr(i,j)-lfn_curr(i,j-1))/dy
            diffRx=(lfn_curr(i+1,j)-lfn_curr(i,j))/dx
            diffRy=(lfn_curr(i,j+1)-lfn_curr(i,j))/dy
            diff2x=select_eno(diffLx,diffRx)
            diff2y=select_eno(diffLy,diffRy)
          else
            select case(fire_upwinding_reinit)
            case(1)
             diff2x=select_4th(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i+1,j),lfn_curr(i+2,j))
             diff2y=select_4th(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j+1),lfn_curr(i,j+2))
             diff2x=select_weno3(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i+1,j),lfn_curr(i+2,j),lfn_s0(i,j)*diff2x)
             diff2y=select_weno3(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j+1),lfn_curr(i,j+2),lfn_s0(i,j)*diff2y)
            case(2)
             diff2x=select_4th(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i+1,j),lfn_curr(i+2,j))
             diff2y=select_4th(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j+1),lfn_curr(i,j+2))
             diff2x=select_weno5(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i-3,j),lfn_curr(i+1,j),lfn_curr(i+2,j),lfn_curr(i+3,j),lfn_s0(i,j)*diff2x)
             diff2y=select_weno5(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j-3),lfn_curr(i,j+1),lfn_curr(i,j+2),lfn_curr(i,j+3),lfn_s0(i,j)*diff2y)
            case(3)
             if (lfn_curr(i,j).lt.threshold_HLu) then
               diff2x=select_4th(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i+1,j),lfn_curr(i+2,j))
               diff2y=select_4th(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j+1),lfn_curr(i,j+2))
               diff2x=select_weno3(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i+1,j),lfn_curr(i+2,j),lfn_s0(i,j)*diff2x)
               diff2y=select_weno3(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j+1),lfn_curr(i,j+2),lfn_s0(i,j)*diff2y)
             else
               diffLx=(lfn_curr(i,j)-lfn_curr(i-1,j))/dx
               diffLy=(lfn_curr(i,j)-lfn_curr(i,j-1))/dy
               diffRx=(lfn_curr(i+1,j)-lfn_curr(i,j))/dx
               diffRy=(lfn_curr(i,j+1)-lfn_curr(i,j))/dy
               diff2x=select_eno(diffLx,diffRx)
               diff2y=select_eno(diffLy,diffRy)
             endif
            case(4)
             if (lfn_curr(i,j).lt.threshold_HLu) then
               diff2x=select_4th(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i+1,j),lfn_curr(i+2,j))
               diff2y=select_4th(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j+1),lfn_curr(i,j+2))
               diff2x=select_weno5(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i-3,j),lfn_curr(i+1,j),lfn_curr(i+2,j),lfn_curr(i+3,j),lfn_s0(i,j)*diff2x)
               diff2y=select_weno5(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j-3),lfn_curr(i,j+1),lfn_curr(i,j+2),lfn_curr(i,j+3),lfn_s0(i,j)*diff2y)
             else
               diffLx=(lfn_curr(i,j)-lfn_curr(i-1,j))/dx
               diffLy=(lfn_curr(i,j)-lfn_curr(i,j-1))/dy
               diffRx=(lfn_curr(i+1,j)-lfn_curr(i,j))/dx
               diffRy=(lfn_curr(i,j+1)-lfn_curr(i,j))/dy
               diff2x=select_eno(diffLx,diffRx)
               diff2y=select_eno(diffLy,diffRy)
             endif
            case default
             if (lfn_curr(i,j).lt.threshold_HLu) then
               diff2x=select_4th(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i+1,j),lfn_curr(i+2,j))
               diff2y=select_4th(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j+1),lfn_curr(i,j+2))
               diff2x=select_weno5(dx,lfn_curr(i,j),lfn_curr(i-1,j),lfn_curr(i-2,j),lfn_curr(i-3,j),lfn_curr(i+1,j),lfn_curr(i+2,j),lfn_curr(i+3,j),lfn_s0(i,j)*diff2x)
               diff2y=select_weno5(dy,lfn_curr(i,j),lfn_curr(i,j-1),lfn_curr(i,j-2),lfn_curr(i,j-3),lfn_curr(i,j+1),lfn_curr(i,j+2),lfn_curr(i,j+3),lfn_s0(i,j)*diff2y)
             else
               diffLx=(lfn_curr(i,j)-lfn_curr(i-1,j))/dx
               diffLy=(lfn_curr(i,j)-lfn_curr(i,j-1))/dy
               diffRx=(lfn_curr(i+1,j)-lfn_curr(i,j))/dx
               diffRy=(lfn_curr(i,j+1)-lfn_curr(i,j))/dy
               diff2x=select_eno(diffLx,diffRx)
               diff2y=select_eno(diffLy,diffRy)
             endif
            end select
          endif
            grad=sqrt(diff2x*diff2x+diff2y*diff2y)
            tend_r=lfn_s0(i,j)*(1.0-grad)
            lfn_fin(i,j)=lfn_ini(i,j)+(dt_s*rk_coeff)*tend_r
        enddo
    enddo

end subroutine advance_ls_reinit





subroutine tign_update(ifts,ifte,jfts,jfte, &              
                       ifms,ifme,jfms,jfme, &              
                       ifds,jfds,ifde,jfde, &              
                       ts,dt,               &              
                       lfn_in,lfn_out,tign  &              
                      )

implicit none

integer,intent(in)::ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,ifds,jfds,ifde,jfde
real,dimension(ifms:ifme,jfms:jfme),intent(inout)::tign
real,dimension(ifms:ifme,jfms:jfme),intent(in)::lfn_in,lfn_out
real,intent(in)::ts,dt
real::time_now
integer::i,j,k,kk
intrinsic epsilon
character(len=128)::msg











    
time_now=ts+dt
time_now = time_now + abs(time_now)*epsilon(time_now)*2.
do j=jfts,jfte
    do i=ifts,ifte
        
        if (.not. lfn_out(i,j)>0 .and. lfn_in(i,j)>0)then
            tign(i,j) = ts + dt * lfn_in(i,j) / (lfn_in(i,j) - lfn_out(i,j))
    endif
        
        if(lfn_out(i,j)>0.)tign(i,j)=time_now
    enddo
enddo

do j=jfts,jfte
  if (j.le.boundary_guard .or. j.gt.(jfde-boundary_guard)) then
    do i=ifts,ifte
      if (lfn_out(i,j).lt.0.) then 
        WRITE(msg,*)'j-boundary reached'
        CALL wrf_debug(0,msg)
        WRITE(msg,*)'i,j,lfn_out=',i,j,lfn_out(i,j)
        CALL wrf_debug(0,msg)
        call crash('wrf: SUCCESS COMPLETE WRF. Fire has reached domain boundary.')
      endif
    enddo
  endif
enddo
do i=ifts,ifte
  if (i.le.boundary_guard .or. i.gt.(ifde-boundary_guard)) then
    do j=jfts,jfte
      if (lfn_out(i,j).lt.0.) then 
        WRITE(msg,*)'i-boundary reached'
        CALL wrf_debug(0,msg)
        WRITE(msg,*)'i,j,lfn_out=',i,j,lfn_out(i,j)
        CALL wrf_debug(0,msg)
        call crash('wrf: SUCCESS COMPLETE WRF. Fire has reached domain boundary.')
      endif
    enddo
  endif
enddo
    
call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, &
               lfn_out,'prop_ls: lfn out')

end subroutine tign_update





subroutine calc_flame_length(ifts,ifte,jfts,jfte, &
                             ifms,ifme,jfms,jfme, &
                             ros,iboros,flame_length,ros_fl,fire_area)




implicit none

integer,intent(in)::ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme
real,dimension(ifms:ifme,jfms:jfme),intent(in)::ros,iboros,fire_area
real,dimension(ifms:ifme,jfms:jfme),intent(out)::flame_length,ros_fl

integer::i,j

do j=jfts,jfte
  do i=ifts,ifte
    if (fire_area(i,j).gt.0.0 .AND. fire_area(i,j).lt.1.0) then
      flame_length(i,j)=0.0775*(iboros(i,j)*ros(i,j))**0.46 
      ros_fl(i,j)=ros(i,j)
    else
      flame_length(i,j)=0.0
      ros_fl(i,j)=0.0
    endif
  enddo
enddo

end subroutine calc_flame_length





subroutine tend_ls( id, &
    ids,ide,jds,jde,    &              
    ips,ipe,jps,jpe,    &              
    its,ite,jts,jte,    &              
    ims,ime,jms,jme,    &              
    t,dt,dx,dy,         &              
    lfn,                &              
    tbound,             &              
    tend, ros,          &              
    fp                  &
)

implicit none




integer,intent(in)::id 
integer,intent(in)::ims,ime,jms,jme,its,ite,jts,jte
integer, intent(in)::ids,ide,jds,jde,ips,ipe,jps,jpe 
real,intent(in)::t                                 
real,intent(in)::dt,dx,dy                          
real,dimension(ims:ime,jms:jme),intent(inout)::lfn 
real,dimension(ims:ime,jms:jme),intent(out)::tend  
real,dimension(ims:ime,jms:jme),intent(out)::ros   
real,intent(out)::tbound                           
type(fire_params),intent(in)::fp


real:: te,diffLx,diffLy,diffRx,diffRy, & 
   diffCx,diffCy,diff2x,diff2y,grad,rr, &
   ros_base,ros_wind,ros_slope,ros_back,advx,advy,scale,nvx,nvy, &
   speed,tanphi
integer::i,j,itso,iteo,jtso,jteo
character(len=128)msg


real, parameter:: eps=epsilon(0.0)

real, parameter:: zero=0.,one=1.,tol=100*eps, &
    safe=2.,rmin=safe*tiny(zero),rmax=huge(zero)/safe

real :: diff2xn, diff2yn
real :: a_valor, signo_x, signo_y
real :: threshold_HLl,threshold_HLu
real :: threshold_AV,fire_viscosity_var
integer,parameter :: bdy_eno1=10



intrinsic max,min,sqrt,nint,tiny,huge




threshold_HLl=-fire_lsm_band_ngp*dx 
threshold_HLu=fire_lsm_band_ngp*dx
threshold_AV=fire_viscosity_ngp*dx
    
    call continue_at_boundary(1,1,fire_lfn_ext_up, &   
    ims,ime,jms,jme, &                                 
    ids,ide,jds,jde, &                                 
    ips,ipe,jps,jpe, &                                 
    its,ite,jts,jte, &                                 
    itso,iteo,jtso,jteo, &                             
    lfn)                                               

    tbound=0    
    do j=jts,jte
        do i=its,ite
            
            diffRx = (lfn(i+1,j)-lfn(i,j))/dx
            diffLx = (lfn(i,j)-lfn(i-1,j))/dx
            diffRy = (lfn(i,j+1)-lfn(i,j))/dy
            diffLy = (lfn(i,j)-lfn(i,j-1))/dy
            diffCx = diffLx+diffRx   
            diffCy = diffLy+diffRy

            if (i.lt.ids+bdy_eno1 .OR. i.gt.ide-bdy_eno1 .OR. j.lt.jds+bdy_eno1 .OR. j.gt.jde-bdy_eno1) then 
              diff2x=select_eno(diffLx,diffRx)
              diff2y=select_eno(diffLy,diffRy)
              grad=sqrt(diff2x*diff2x + diff2y*diff2y)
            else
              select case(fire_upwinding)
              case(0)  
                  grad=sqrt(diffCx**2 + diffCy**2)
              case(1) 
                  diff2x=select_upwind(diffLx,diffRx)
                  diff2y=select_upwind(diffLy,diffRy)
                  grad=sqrt(diff2x*diff2x + diff2y*diff2y)
              case(2) 
                  diff2x=select_godunov(diffLx,diffRx)
                  diff2y=select_godunov(diffLy,diffRy)
                  grad=sqrt(diff2x*diff2x + diff2y*diff2y)
              case(3) 
                  diff2x=select_eno(diffLx,diffRx)
                  diff2y=select_eno(diffLy,diffRy)
                  grad=sqrt(diff2x*diff2x + diff2y*diff2y)
              case(4) 
                  grad=sqrt(max(diffLx,0.)**2+min(diffRx,0.)**2   &
                          + max(diffLy,0.)**2+min(diffRy,0.)**2)
              case(5) 
                  diff2x=select_2nd(rr,dx,lfn(i,j),lfn(i-1,j),lfn(i+1,j))
                  diff2y=select_2nd(rr,dy,lfn(i,j),lfn(i,j-1),lfn(i,j+1))
                  grad=sqrt(diff2x*diff2x + diff2y*diff2y)
              case(6) 
                  a_valor=select_4th(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j))*fp%vx(i,j)+ &
                          select_4th(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))*fp%vy(i,j)
                  signo_x=a_valor*select_4th(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j))
                  signo_y=a_valor*select_4th(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))
                  diff2x=select_weno3(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j),signo_x)
                  diff2y=select_weno3(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2),signo_y)
                  grad=sqrt(diff2x*diff2x + diff2y*diff2y)
              case(7) 
                  a_valor=select_4th(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j))*fp%vx(i,j)+ &
                          select_4th(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))*fp%vy(i,j)
                  signo_x=a_valor*select_4th(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j))
                  signo_y=a_valor*select_4th(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))
                  diff2x=select_weno5(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i-3,j),lfn(i+1,j),lfn(i+2,j),lfn(i+3,j),signo_x)
                  diff2y=select_weno5(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j-3),lfn(i,j+1),lfn(i,j+2),lfn(i,j+3),signo_y)
                  grad=sqrt(diff2x*diff2x + diff2y*diff2y)
              case(8) 
                  if (abs(lfn(i,j)).lt.threshold_HLu) then
                    a_valor=select_4th(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j))*fp%vx(i,j)+ &
                            select_4th(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))*fp%vy(i,j)
                    signo_x=a_valor*select_4th(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j))
                    signo_y=a_valor*select_4th(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))
                    diff2x=select_weno3(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j),signo_x)
                    diff2y=select_weno3(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2),signo_y)
                    grad=sqrt(diff2x*diff2x + diff2y*diff2y)
                  else
                    diff2x=select_eno(diffLx,diffRx)
                    diff2y=select_eno(diffLy,diffRy)
                    grad=sqrt(diff2x*diff2x + diff2y*diff2y)
                  endif
              case(9) 
                  if (abs(lfn(i,j)).lt.threshold_HLu) then
                    a_valor=select_4th(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j))*fp%vx(i,j)+ &
                            select_4th(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))*fp%vy(i,j)
                    signo_x=a_valor*select_4th(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i+1,j),lfn(i+2,j))
                    signo_y=a_valor*select_4th(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))
                    diff2x=select_weno5(dx,lfn(i,j),lfn(i-1,j),lfn(i-2,j),lfn(i-3,j),lfn(i+1,j),lfn(i+2,j),lfn(i+3,j),signo_x)
                    diff2y=select_weno5(dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j-3),lfn(i,j+1),lfn(i,j+2),lfn(i,j+3),signo_y)
                    grad=sqrt(diff2x*diff2x + diff2y*diff2y)
                  else
                    diff2x=select_eno(diffLx,diffRx)
                    diff2y=select_eno(diffLy,diffRy)
                    grad=sqrt(diff2x*diff2x + diff2y*diff2y)
                  endif
              case default
                   grad=0.
              end select
            endif

            scale=sqrt(grad**2.0+eps) 
            nvx=diff2x/scale 
            nvy=diff2y/scale 

            
            speed =  fp%vx(i,j)*nvx + fp%vy(i,j)*nvy
            
            call fire_ros(ros_base,ros_wind,ros_slope, &
            nvx,nvy,i,j,fp)
            rr=ros_base + ros_wind + fire_slope_factor*ros_slope

            if(fire_grows_only.gt.0)rr=max(rr,0.)
            
            if(i.ge.its.and.i.le.ite.and.j.ge.jts.and.j.le.jte)ros(i,j)=rr

            if(fire_upwind_split.eq.0)then

                
                te = -rr*grad   

            else

                
                te = - ros_base*grad

		
                if (abs(speed)> eps) then
                    advx=fp%vx(i,j)*ros_wind/speed
                    advy=fp%vy(i,j)*ros_wind/speed
                else 
                    advx=0
                    advy=0
                endif

                tanphi =  fp%dzdxf(i,j)*nvx + fp%dzdyf(i,j)*nvy
		
                if(abs(tanphi)>eps) then
                    advx=advx+fp%dzdxf(i,j)*ros_slope/tanphi
                    advy=advy+fp%dzdyf(i,j)*ros_slope/tanphi
                endif

                if(fire_upwind_split.eq.1)then   

                    
                    te = te - max(advx,0.)*diffLx - min(advx,0.)*diffRy &
                            - max(advy,0.)*diffLy - min(advy,0.)*diffRy


                elseif(fire_upwind_split.eq.2)then   
 
                    
                    call crash('prop_ls: bad fire_upwind_split, Lax-Friedrichs not done yet')

                else

                    call crash('prop_ls: bad fire_upwind_split')

                endif
            endif

            
            if (grad > 0.) then
                 tbound = max(tbound,rr*(abs(diff2x)/dx+abs(diff2y)/dy)/grad)
            endif

            
            if (abs(lfn(i,j)).lt.threshold_AV .AND. (i.gt.ids+bdy_eno1 .and. i.lt.ide-bdy_eno1) .AND. (j.gt.jds+bdy_eno1 .and. j.lt.jde-bdy_eno1)) then 
              fire_viscosity_var=fire_viscosity_bg
            elseif (abs(lfn(i,j)).ge.threshold_AV .AND. abs(lfn(i,j)).lt.threshold_AV*(1.0+fire_viscosity_band) .AND. (i.gt.ids+10 .and. i.lt.ide-10) .AND. (j.gt.jds+10 .and. j.lt.jde-10)) then
              fire_viscosity_var=min(fire_viscosity_bg+(fire_viscosity-fire_viscosity_bg)*(abs(lfn(i,j))-threshold_AV)/(fire_viscosity_band*threshold_AV),fire_viscosity)
            else
              fire_viscosity_var=fire_viscosity
            endif
            te=te + fire_viscosity_var*abs(rr)*((diffRx-diffLx)+(diffRy-diffLy))
            tend(i,j)=te
        enddo
    enddo        

    call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme, & 
                   tend,'tend_ls: tend out')

    
    tbound = 1/(tbound+tol)

end subroutine tend_ls






real function select_upwind(diffLx,diffRx)
implicit none
real, intent(in):: diffLx, diffRx
real diff2x



diff2x=0
if (diffLx>0.and.diffRx>0.)diff2x=diffLx
if (diffLx<0.and.diffRx<0.)diff2x=diffRx

select_upwind=diff2x
end function select_upwind






real function select_godunov(diffLx,diffRx)
implicit none
real, intent(in):: diffLx, diffRx
real diff2x,diffCx






diff2x=0
diffCx=diffRx+diffLx
if (diffLx>0.and..not.diffCx<0)diff2x=diffLx
if (diffRx<0.and.     diffCx<0)diff2x=diffRx

select_godunov=diff2x
end function select_godunov





real function select_eno(diffLx,diffRx)
implicit none
real, intent(in):: diffLx, diffRx
real diff2x



if    (.not.diffLx>0 .and. .not.diffRx>0)then
    diff2x=diffRx
elseif(.not.diffLx<0 .and. .not.diffRx<0)then
    diff2x=diffLx
elseif(.not.diffLx<0 .and. .not.diffRx>0)then
    if(.not. abs(diffRx) < abs(diffLx))then
        diff2x=diffRx
    else
        diff2x=diffLx
    endif
else
    diff2x=0.
endif

select_eno=diff2x
end function select_eno
      





real function select_2nd(ros,dx,lfn_i,lfn_im1,lfn_ip1)
implicit none
real, intent(in):: lfn_i, lfn_im1, lfn_ip1
real, intent(in):: ros, dx
real:: diff2x_p, diff2x_m 



diff2x_p=0.
diff2x_m=0.
diff2x_p=(lfn_ip1+lfn_i)/(2.*dx)
diff2x_m=(lfn_i+lfn_im1)/(2.*dx)
select_2nd=diff2x_p-diff2x_m
end function select_2nd






real function select_4th(dx,lfn_i,lfn_im1,lfn_im2,lfn_ip1,lfn_ip2)
implicit none
real, intent(in):: lfn_i,lfn_im1,lfn_im2,lfn_ip1,lfn_ip2
real, intent(in):: dx
real:: diff2x_p,diff2x_m 



diff2x_p=0.
diff2x_m=0.
diff2x_p=(7.0*lfn_ip1+7.0*lfn_i-lfn_ip2-lfn_im1)/(12.*dx)
diff2x_m=(7.0*lfn_i+7.0*lfn_im1-lfn_ip1-lfn_im2)/(12.*dx)
select_4th=diff2x_p-diff2x_m
end function select_4th






real function select_weno3(dx,lfn_it,lfn_im1t,lfn_im2t,lfn_ip1t,lfn_ip2t,uf)
implicit none
real, intent(in)::lfn_it,lfn_im1t,lfn_im2t,lfn_ip1t,lfn_ip2t,uf
real:: lfn_i,lfn_im1,lfn_im2,lfn_ip1,lfn_ip2
real, intent(in):: dx
real:: flux_p,flux_m 
real, parameter:: gamma1=1./3.,gamma2=2./3.,tol=1e-6
real:: w1,w2,w1t,w2t,beta1,beta2 
real:: fh_1,fh_2



flux_p=0.
flux_m=0.

if (uf .ge. 0.0 ) then
  lfn_i=lfn_it
  lfn_im1=lfn_im1t
  lfn_im2=lfn_im2t
  lfn_ip1=lfn_ip1t
else
  lfn_i=lfn_it
  lfn_im1=lfn_ip1t
  lfn_im2=lfn_ip2t
  lfn_ip1=lfn_im1t
endif


fh_1=-0.5*lfn_im1+1.5*lfn_i
fh_2=0.5*lfn_i+0.5*lfn_ip1
beta1=(lfn_i-lfn_im1)**2
beta2=(lfn_ip1-lfn_i)**2
w1t=gamma1/(beta1+tol)**2
w2t=gamma2/(beta2+tol)**2
w1=w1t/(w1t+w2t)
w2=w2t/(w1t+w2t)
flux_p=w1*fh_1+w2*fh_2

fh_1=-0.5*lfn_im2+1.5*lfn_im1
fh_2=0.5*lfn_im1+0.5*lfn_i
beta1=(lfn_im1-lfn_im2)**2
beta2=(lfn_i-lfn_im1)**2
w1t=gamma1/(beta1+tol)**2
w2t=gamma2/(beta2+tol)**2
w1=w1t/(w1t+w2t)
w2=w2t/(w1t+w2t)
flux_m=w1*fh_1+w2*fh_2

if (uf .ge. 0.0 ) then
  select_weno3=(flux_p-flux_m)/dx
else
  select_weno3=(flux_m-flux_p)/dx
endif
end function select_weno3






real function select_weno5(dx,lfn_it,lfn_im1t,lfn_im2t,lfn_im3t,lfn_ip1t,lfn_ip2t,lfn_ip3t,uf)
implicit none
real, intent(in):: lfn_it,lfn_im1t,lfn_im2t,lfn_im3t,lfn_ip1t,lfn_ip2t,lfn_ip3t,uf
real:: lfn_i,lfn_im1,lfn_im2,lfn_im3,lfn_ip1,lfn_ip2,lfn_ip3
real, intent(in):: dx
real:: flux_p,flux_m 
real, parameter:: gamma1=1./10.,gamma2=3./5.,gamma3=3./10.,tol=1e-6
real:: w1,w2,w3,w1t,w2t,w3t,beta1,beta2,beta3 
real:: fh_1,fh_2,fh_3



flux_p=0.
flux_m=0.

if (uf .ge. 0.0 ) then
  lfn_i=lfn_it
  lfn_im1=lfn_im1t
  lfn_im2=lfn_im2t
  lfn_im3=lfn_im3t
  lfn_ip1=lfn_ip1t
  lfn_ip2=lfn_ip2t
else
  lfn_i=lfn_it
  lfn_im1=lfn_ip1t
  lfn_im2=lfn_ip2t
  lfn_im3=lfn_ip3t
  lfn_ip1=lfn_im1t
  lfn_ip2=lfn_im2t
endif


fh_1=(2.0*lfn_im2-7.0*lfn_im1+11.0*lfn_i)/6.0
fh_2=(-1.0*lfn_im1+5.0*lfn_i+2.0*lfn_ip1)/6.0
fh_3=(2.0*lfn_i+5.0*lfn_ip1-1.0*lfn_ip2)/6.0
beta1=(13.0/12.0)*(lfn_im2-2.0*lfn_im1+lfn_i)**2+0.25*(lfn_im2-4.0*lfn_im1+3.0*lfn_i)**2
beta2=(13.0/12.0)*(lfn_im1-2.0*lfn_i+lfn_ip1)**2+0.25*(lfn_im1-lfn_ip1)**2
beta3=(13.0/12.0)*(lfn_i-2.0*lfn_ip1+lfn_ip2)**2+0.25*(3.0*lfn_i-4.0*lfn_ip1+lfn_ip2)**2
w1t=gamma1/(beta1+tol)**2
w2t=gamma2/(beta2+tol)**2
w3t=gamma3/(beta3+tol)**2
w1=w1t/(w1t+w2t+w3t)
w2=w2t/(w1t+w2t+w3t)
w3=w3t/(w1t+w2t+w3t)
flux_p=w1*fh_1+w2*fh_2+w3*fh_3

fh_1=(2.0*lfn_im3-7.0*lfn_im2+11.0*lfn_im1)/6.0
fh_2=(-1.0*lfn_im2+5.0*lfn_im1+2.0*lfn_i)/6.0
fh_3=(2.0*lfn_im1+5.0*lfn_i-1.0*lfn_ip1)/6.0
beta1=(13.0/12.0)*(lfn_im3-2.0*lfn_im2+lfn_im1)**2+0.25*(lfn_im3-4.0*lfn_im2+3.0*lfn_im1)**2
beta2=(13.0/12.0)*(lfn_im2-2.0*lfn_im1+lfn_i)**2+0.25*(lfn_im2-lfn_i)**2
beta3=(13.0/12.0)*(lfn_im1-2.0*lfn_i+lfn_ip1)**2+0.25*(3.0*lfn_im1-4.0*lfn_i+lfn_ip1)**2
w1t=gamma1/(beta1+tol)**2
w2t=gamma2/(beta2+tol)**2
w3t=gamma3/(beta3+tol)**2
w1=w1t/(w1t+w2t+w3t)
w2=w2t/(w1t+w2t+w3t)
w3=w3t/(w1t+w2t+w3t)
flux_m=w1*fh_1+w2*fh_2+w3*fh_3

if (uf .ge. 0.0 ) then
  select_weno5=(flux_p-flux_m)/dx
else
  select_weno5=(flux_m-flux_p)/dx
endif
end function select_weno5





real function speed_func(diffCx,diffCy,dx,dy,i,j,fp)


implicit none

real, intent(in)::diffCx,diffCy  
real, intent(in)::dx,dy  
integer, intent(in)::i,j         
type(fire_params),intent(in)::fp

real::scale,nvx,nvy,r
real::ros_base , ros_wind , ros_slope
real, parameter:: eps=epsilon(0.0)

            
            scale=sqrt(diffCx*diffCx+diffCy*diffCy+eps) 
            nvx=diffCx/scale
            nvy=diffCy/scale
                      
            

            call fire_ros(ros_base,ros_wind,ros_slope, &
            nvx,nvy,i,j,fp)

            r=ros_base + ros_wind + ros_slope
            if(fire_grows_only.gt.0)r=max(r,0.)
            speed_func=r

end function speed_func

end module module_fr_fire_core
