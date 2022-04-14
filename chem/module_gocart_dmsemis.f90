module module_dms_emis
  USE module_data_gocartchem
contains

        subroutine gocart_dmsemis(dt,config_flags,alt,t_phy,u_phy,  &
         v_phy,chem,rho_phy,dz8w,u10,v10,p8w,dms_0,tsk,                  &
         ivgtyp,isltyp,xland,dx,g, &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_configure
  USE module_state_description
  USE module_model_constants, only : mwdry
  IMPLICIT NONE
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, jms:jme),                                     &
         INTENT(IN ) ::    dms_0,tsk
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::                              alt,               &
                                                      t_phy,               &
                                                      dz8w,                &
                                              p8w,u_phy,v_phy,rho_phy
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     xland
   real, intent(in) :: dx,g,dt
  



  integer :: i,j,k,ndt,imx,jmx,lmx,nmx
  integer,dimension (1,1) :: ilwi
  real*8, DIMENSION (1,1,1,1) :: tc
  real*8, DIMENSION (1,1,1) :: bems,airmas
  real*8, DIMENSION (1,1) :: emsdms
  real*8, dimension (1,1) :: w10m,gwet,airden,tskin,dmso
  real*8, dimension (1) :: dxy



  imx=1
  jmx=1
  lmx=1
  nmx=1
  k=kts
  ndt=ifix(dt)
  do j=jts,jte
  do i=its,ite



     if(xland(i,j).gt.1.5)then
     ilwi(1,1)=0

    tc(1,1,1,1)=chem(i,kts,j,p_dms)*1.d-6
    dmso(1,1)=dms_0(i,j)
     w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
  
     airmas(1,1,1)=-1.0*(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g
     dxy(1)=dx*dx
     tskin(1,1)=tsk(i,j)



     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))

    call  srcdms(imx, jmx, lmx, nmx, ndt, tc, mwdry,&
                  tskin, ilwi, dmso, w10m, airmas, dxy, emsdms, bems)
    chem(i,kts,j,p_dms)=tc(1,1,1,1)*1.e6
  endif
  enddo
  enddo

end subroutine gocart_dmsemis

SUBROUTINE srcdms(imx, jmx, lmx, nmx, ndt1, tc,airmw, &
                  tskin, ilwi, dmso, w10m, airmas, dxy, emsdms, bems)
 
  
  
  
  
  


  IMPLICIT NONE
  REAL,    PARAMETER :: dms_mw = 62.00
  REAL,    PARAMETER :: tcmw(1)=dms_mw
  INTEGER, INTENT(IN)    :: imx, jmx, lmx, nmx, ndt1
  REAL*8,    INTENT(IN)    :: tskin(imx,jmx), dmso(imx,jmx)
  INTEGER, INTENT(IN)    :: ilwi(imx,jmx)
  REAL*8,    INTENT(IN)    :: dxy(jmx), w10m(imx,jmx)
  REAL,      INTENT(IN)    :: airmw
  REAL*8,    INTENT(IN)    :: airmas(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(INOUT) :: emsdms(imx,jmx)
  REAL*8,    INTENT(OUT)   :: bems(imx,jmx,nmx)

  INTEGER :: i,j
  REAL*8    :: sst, sc, conc, w10, scco2, akw, erate, dmssrc, c

  
  
  
  
  

  

  lat_loop: DO j = 1,jmx
     lon_loop: DO i = 1,imx
        
        sst = tskin(i,j) - 273.15
        if_water: IF (ilwi(i,j) == 0) THEN
           
           
           sc = 2674.0 - 147.12*sst + 3.726*(sst**2) - 0.038*(sst**3)














           conc = dmso(i,j)

           w10  = w10m(i,j)















           
           scco2 = 600.0
           IF (w10 <= 3.6) THEN
              akw = 0.17 * w10
           ELSE IF (w10 <= 13.0) THEN
              akw = 2.85 * w10 - 9.65
           ELSE
              akw = 5.90 * w10 - 49.3
           END IF
           
           
           IF (w10 <= 3.6) THEN
              akw = akw * ((scco2/sc) ** 0.667)
           ELSE
              akw = akw * SQRT(scco2/sc)
           END IF











           erate  = akw/100.0/3600.0*conc*1.0E-12*1000.0*REAL(ndt1)*tcmw(NDMS)
           dmssrc = erate * dxy(j)

        ELSE   

           dmssrc = 0.0
           
        END IF if_water 





        
        c = dmssrc / airmas(i,j,1) * airmw / tcmw(NDMS)
        tc(i,j,1,NDMS) = tc(i,j,1,NDMS) + c

        
        
        
        emsdms(i,j) = emsdms(i,j) + dmssrc * smw / tcmw(NDMS) 

        bems(i,j,NDMS) = dmssrc  

     END DO lon_loop
  END DO lat_loop

END SUBROUTINE srcdms

END module module_dms_emis
