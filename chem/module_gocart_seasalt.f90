MODULE GOCART_SEASALT

CONTAINS
  subroutine gocart_seasalt_driver(ktau,dt,config_flags,julday,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,u10,v10,p8w,z_at_w,                  &
         xland,xlat,xlong,dx,g,emis_seas,seasin, &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_configure
  USE module_state_description
  USE module_model_constants, ONLY: mwdry
  IMPLICIT NONE
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: julday, ktau,                     &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_seas),OPTIONAL,&
         INTENT(INOUT ) ::                                                 &
         emis_seas
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     xland,                &
                                                     xlat,                 &
                                                     xlong
   REAL,  DIMENSION( ims:ime , jms:jme, 5 ),                        &
          INTENT(INOUT   ) ::              seasin
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                      t_phy,               &
                                               dz8w,p8w,z_at_w,            &
                                              u_phy,v_phy,rho_phy

  REAL, INTENT(IN   ) :: dt,dx,g



  integer :: ipr,nmx,i,j,k,ndt,imx,jmx,lmx
  integer,dimension (1,1) :: ilwi
  real*8, DIMENSION (4) :: tc,bems
  real*8, dimension (1,1) :: w10m,gwet,airden,airmas
  real*8, dimension (1) :: dxy
  real*8 conver,converi
  conver=1.d-9
  converi=1.d9



  imx=1
  jmx=1
  lmx=1
  nmx=4
  k=kts
  do j=jts,jte
  do i=its,ite



     if(xland(i,j).gt.1.5.and.z_at_w(i,kts,j).lt.1.e-3)then
     ilwi(1,1)=0
     if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
        tc(:)=1.e-16*conver
     else
     tc(1)=chem(i,kts,j,p_seas_1)*conver
     tc(2)=chem(i,kts,j,p_seas_2)*conver
     tc(3)=chem(i,kts,j,p_seas_3)*conver
     tc(4)=chem(i,kts,j,p_seas_4)*conver
     endif
     w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
     airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g



     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))

     dxy(1)=dx*dx
       ipr=0

    call source_ss( imx,jmx,lmx,nmx, dt, tc,ilwi, dxy, w10m, airmas, bems,ipr)

    if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
      seasin(i,j,1:4)=tc(1:4)*converi
     else
     chem(i,kts,j,p_seas_1)=tc(1)*converi
     chem(i,kts,j,p_seas_2)=tc(2)*converi
     chem(i,kts,j,p_seas_3)=tc(3)*converi
     chem(i,kts,j,p_seas_4)=tc(4)*converi
     endif

     emis_seas(i,1,j,p_eseas1)=bems(1)
     emis_seas(i,1,j,p_eseas2)=bems(2)
     emis_seas(i,1,j,p_eseas3)=bems(3)
     emis_seas(i,1,j,p_eseas4)=bems(4)
     endif
  enddo
  enddo


end subroutine gocart_seasalt_driver

SUBROUTINE source_ss(imx,jmx,lmx,nmx, dt1, tc, &
                     ilwi, dxy, w10m, airmas, &
                     bems,ipr)
































 

  USE module_data_gocart_seas

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nmx,imx,jmx,lmx,ipr
  INTEGER, INTENT(IN)    :: ilwi(imx,jmx)
  REAL*8,    INTENT(IN)    :: dxy(jmx), w10m(imx,jmx)
  REAL*8,    INTENT(IN)    :: airmas(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(OUT)   :: bems(imx,jmx,nmx)

  REAL*8 :: c0(5), b0(2)


  
  REAL*8, PARAMETER :: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/) 
  REAL*8, PARAMETER :: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
  REAL*8, PARAMETER :: b_old(2)=(/0.380, 0.650/)
  REAL*8, PARAMETER :: b_new(2)=(/0.433, 0.433/)
  REAL*8, PARAMETER :: dr=5.0D-2 
  REAL*8, PARAMETER :: theta=30.0
  

  REAL*8,    PARAMETER :: frh = 2.d0
  LOGICAL, PARAMETER :: old=.TRUE., new=.FALSE.
  REAL*8 :: rho_d, r0, r1, r, r_w, a, b, dfn, r_d, dfm, src
  INTEGER :: i, j, n, nr, ir
  REAL :: dt1


  REAL*8                  :: tcmw(nmx), ar(nmx), tcvv(nmx)
  REAL*8                  :: ar_wetdep(nmx), kc(nmx)
  CHARACTER(LEN=20)     :: tcname(nmx), tcunits(nmx)
  LOGICAL               :: aerosol(nmx)


  REAL*8 :: tc1(imx,jmx,lmx,nmx)
  REAL*8, TARGET :: tcms(imx,jmx,lmx,nmx) 
  REAL*8, TARGET :: tcgm(imx,jmx,lmx,nmx) 

  
  
  



  
  
  
  REAL*8 :: e_an(imx,jmx,2,nmx), e_bb(imx,jmx,nmx), &
          e_ac(imx,jmx,lmx,nmx)

  
  
  








  
  REAL*8 :: tmas0(nmx), tmas1(nmx)
  REAL*8 :: dtems(nmx), dttrp(nmx), dtdif(nmx), dtcnv(nmx)
  REAL*8 :: dtwet(nmx), dtdry(nmx), dtstl(nmx)
  REAL*8 :: dtems2(nmx), dttrp2(nmx), dtdif2(nmx), dtcnv2(nmx)
  REAL*8 :: dtwet2(nmx), dtdry2(nmx), dtstl2(nmx)

  
  REAL*8, TARGET :: ems_an(imx,jmx,nmx),    ems_bb(imx,jmx,nmx), ems_tp(imx,jmx)
  REAL*8, TARGET :: ems_ac(imx,jmx,lmx,nmx)
  REAL*8, TARGET :: ems_co(imx,jmx,nmx)


  

  DO n = 1,nmx

     bems(:,:,n) = 0.0
     rho_d = den_seas(n)
     r0 = ra(n)*frh
     r1 = rb(n)*frh
     r = r0
     nr = INT((r1-r0)/dr+.001)

     DO ir = 1,nr
        r_w = r + dr*0.5
        r = r + dr
        IF (new) THEN
           a = 4.7*(1.0 + theta*r_w)**(-0.017*r_w**(-1.44))
           c0 = c_new
           b0 = b_new
        ELSE
           a = 3.0
           c0 = c_old
           b0 = b_old
        END IF
        
        b = (b0(1) - LOG10(r_w))/b0(2)
        dfn = (c0(1)/r_w**a)*(1.0 + c0(3)*r_w**c0(4))* &
             10**(c0(5)*EXP(-(b**2)))
        
        r_d = r_w/frh*1.0D-6  
        dfm = 4.0/3.0*pi*r_d**3*rho_d*frh*dfn*dr*dt1
        DO i = 1,imx
           DO j = 1,jmx

              IF (ilwi(i,j) == 0) THEN

                 src = dfm*dxy(j)*w10m(i,j)**c0(2)

                 if(src < 0.0 ) src=0.
                 tc(i,j,1,n) = tc(i,j,1,n) + src/airmas(i,j,1)

              ELSE
                 src = 0.0
              END IF
              bems(i,j,n) = bems(i,j,n) + src
           END DO  
        END DO 
     END DO 
  END DO 

END SUBROUTINE source_ss
END MODULE GOCART_SEASALT
