



module module_fr_fire_atm

use module_model_constants, only: cp,xlv
use module_fr_fire_util
use module_state_description, only: num_tracer 
use module_state_description, only: p_fire_smoke

contains


subroutine add_fire_tracer_emissions(    &
       tracer_opt,dt,dx,dy,              &
       ifms,ifme,jfms,jfme,              &
       ifts,ifte,jtfs,jfte,              &
       ids,ide,kds,kde,jds,jde,          &
       ims,ime,kms,kme,jms,jme,          &
       its,ite,kts,kte,jts,jte,          &
       rho,dz8w,                         &
       burnt_area_dt,fgip,               &
       tracer,fire_tracer_smoke          &
)

implicit none

integer,intent(in)::tracer_opt
real,intent(in)::fire_tracer_smoke
real,intent(in)::dt,dx,dy
integer,intent(in)::ifms,ifme,jfms,jfme,ifts,ifte,jtfs,jfte,ids,ide,kds,kde,jds,jde,ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte
real,intent(in)::rho(ims:ime,kms:kme,jms:jme),dz8w(ims:ime,kms:kme,jms:jme)
real,intent(in),dimension(ifms:ifme,jfms:jfme)::burnt_area_dt,fgip
real,intent(inout)::tracer(ims:ime,kms:kme,jms:jme,num_tracer)

integer::isz1,jsz1,isz2,jsz2,ir,jr
integer::i,j,ibase,jbase,i_f,ioff,j_f,joff
real::avgw,emis,conv

isz1 = ite-its+1
jsz1 = jte-jts+1
isz2 = ifte-ifts+1
jsz2 = jfte-jtfs+1
ir=isz2/isz1
jr=jsz2/jsz1
avgw = 1.0/(ir*jr)

do j=max(jds+1,jts),min(jte,jde-2)
    jbase=jtfs+jr*(j-jts)
    do i=max(ids+1,its),min(ite,ide-2)
       ibase=ifts+ir*(i-its)
       do joff=0,jr-1
           j_f=joff+jbase
           do ioff=0,ir-1
               i_f=ioff+ibase
               if (num_tracer >0)then
                     emis=avgw*fire_tracer_smoke*burnt_area_dt(i_f,j_f)*fgip(i_f,j_f)*1000/(rho(i,kts,j)*dz8w(i,kts,j)) 
                     tracer(i,kts,j,p_fire_smoke)=tracer(i,kts,j,p_fire_smoke)+emis
                endif
           enddo
       enddo
    enddo
enddo

end subroutine add_fire_tracer_emissions





SUBROUTINE fire_tendency( &
    ids,ide, kds,kde, jds,jde,   & 
    ims,ime, kms,kme, jms,jme,   &
    its,ite, kts,kte, jts,jte,   &
    grnhfx,grnqfx,canhfx,canqfx, & 
    alfg,alfc,z1can,             & 
    zs,z_at_w,dz8w,mu,c1h,c2h,rho, &
    rthfrten,rqvfrten)             








   IMPLICIT NONE



   INTEGER , INTENT(in) :: ids,ide, kds,kde, jds,jde, &
                           ims,ime, kms,kme, jms,jme, &
                           its,ite, kts,kte, jts,jte

   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: grnhfx,grnqfx  
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: canhfx,canqfx  
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: zs  
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: mu  
   REAL, INTENT(in), DIMENSION( kms:kme         ) :: c1h, c2h 

   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: z_at_w 
   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: dz8w   
   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: rho    

   REAL, INTENT(in) :: alfg 
   REAL, INTENT(in) :: alfc 
   REAL, INTENT(in) :: z1can    



   REAL, INTENT(out), DIMENSION( ims:ime,kms:kme,jms:jme ) ::   &
       rthfrten, & 
       rqvfrten    


   INTEGER :: i,j,k
   INTEGER :: i_st,i_en, j_st,j_en, k_st,k_en

   REAL :: cp_i
   REAL :: rho_i
   REAL :: xlv_i
   REAL :: z_w
   REAL :: fact_g, fact_c
   REAL :: alfg_i, alfc_i

   REAL, DIMENSION( its:ite,kts:kte,jts:jte ) :: hfx,qfx
   


        do j=jts,jte
            do k=kts,min(kte+1,kde)
               do i=its,ite
                   rthfrten(i,k,j)=0.
                   rqvfrten(i,k,j)=0.
               enddo
            enddo
        enddo



   

   cp_i = 1./cp     
   xlv_i = 1./xlv   
   alfg_i = 1./alfg
   alfc_i = 1./alfc




   call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,grnhfx,'fire_tendency:grnhfx')
   call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,grnqfx,'fire_tendency:grnqfx')



   i_st = MAX(its,ids+1)
   i_en = MIN(ite,ide-1)
   k_st = kts
   k_en = MIN(kte,kde-1)
   j_st = MAX(jts,jds+1)
   j_en = MIN(jte,jde-1)



   DO j = j_st,j_en
      DO k = k_st,k_en
         DO i = i_st,i_en

            

            z_w = z_at_w(i,k,j) - zs(i,j) 

            

            fact_g = cp_i * EXP( - alfg_i * z_w )
            IF ( z_w < z1can ) THEN
               fact_c = cp_i
            ELSE
               fact_c = cp_i * EXP( - alfc_i * (z_w - z1can) )
            END IF
            hfx(i,k,j) = fact_g * grnhfx(i,j) + fact_c * canhfx(i,j)





            

            fact_g = xlv_i * EXP( - alfg_i * z_w )
            IF (z_w < z1can) THEN
               fact_c = xlv_i
            ELSE
               fact_c = xlv_i * EXP( - alfc_i * (z_w - z1can) )
            END IF
            qfx(i,k,j) = fact_g * grnqfx(i,j) + fact_c * canqfx(i,j)
            






         END DO
      END DO
   END DO






   DO j = j_st,j_en
      DO k = k_st,k_en-1
         DO i = i_st,i_en

            rho_i = 1./rho(i,k,j)

            rthfrten(i,k,j) = - (c1h(k)*mu(i,j)+c2h(k)) * rho_i * (hfx(i,k+1,j)-hfx(i,k,j)) / dz8w(i,k,j)
            rqvfrten(i,k,j) = - (c1h(k)*mu(i,j)+c2h(k)) * rho_i * (qfx(i,k+1,j)-qfx(i,k,j)) / dz8w(i,k,j)

         END DO
      END DO
   END DO

   call print_3d_stats(its,ite,kts,kte,jts,jte,ims,ime,kms,kme,jms,jme,rthfrten,'fire_tendency:rthfrten')
   call print_3d_stats(its,ite,kts,kte,jts,jte,ims,ime,kms,kme,jms,jme,rqvfrten,'fire_tendency:rqvfrten')

   RETURN

END SUBROUTINE fire_tendency





end module module_fr_fire_atm
