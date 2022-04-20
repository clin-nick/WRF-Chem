


MODULE module_dust_emis
CONTAINS






  subroutine bulk_dust_emis (ktau,dt,num_soil_layers,u_phy,v_phy,          &
         rho_phy,alt,u10,v10,p8w,dz8w,smois,erod,                          &
         ivgtyp,isltyp,vegfra,albbck,xland,dx,g,                           &
         nifa2d,                                                           &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )





































  USE module_data_gocart_dust

  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ktau, num_soil_layers,           &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) ::                        &
                                                     vegfra,               &
                                                     albbck
   REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(INOUT) ::                               smois
   REAL,  DIMENSION( ims:ime , jms:jme, 3 )                   ,               &
          INTENT(IN   ) ::    erod
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     xland
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                     dz8w,p8w,             &
                                              u_phy,v_phy,rho_phy

  REAL, INTENT(IN   ) :: dt,dx,g
  REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: NIFA2D



  integer :: i,j,k, nmx, month, n, m, maxflux
  real*8  w10m,gwet,den,diam,airden,airmas,erodin
  real*8  dxy
  real*8  converi
  REAL*8  u_ts0, u_ts, dsrc, srce
  real rhoa, g0, totalemis, pi, dustmas

  converi=1.e9
  g0 = g*1.0E2
  pi = 3.14159
  dxy=dx*dx
  maxflux = 0

  do j=jts,jte
  do i=its,ite
     nifa2d(i,j) = 0.0
  enddo
  enddo

  nmx = 5

  do n = 1, nmx

  den = den_dust(n)*1.0D-3            
  diam = 2.0*reff_dust(n)*1.0D2       
  dustmas = (pi/6.) * den * (diam)**3 
  ch_dust(n,:)=1.0D-9  
  month = 3   
  m = ipoint(n)  

  k=kts
  do j=jts,jte
  do i=its,ite

    dsrc = 0.0




    if(xland(i,j).lt.1.5) then

      w10m=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
      airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dxy/g   


      if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))


      w10m = w10m*1.33

      erodin=erod(i,j,m)


      if (erodin .gt. 1.E-8 .AND. albbck(i,j).gt.0.175 .and. vegfra(i,j).lt.12.5) then
         erodin = MIN(0.5, erodin + 0.1*albbck(i,j))
      endif


      gwet=smois(i,1,j)/porosity(isltyp(i,j))
      airden=rho_phy(i,kts,j)
      rhoa = airden*1.0D-3


      u_ts0 = 0.13*1.0D-2*SQRT(den*g0*diam/rhoa)* &
              SQRT(1.0+0.006/den/g0/(diam)**2.5)/ &
              SQRT(1.928*(1331.0*(diam)**1.56+0.38)**0.092-1.0)


      IF (gwet < 0.5) THEN  

         u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet))))
      ELSE

         u_ts = 100.0
      END IF
      srce = frac_s(n)*erodin*dxy  

      dsrc = MAX(0.0, ch_dust(n,month)*srce*w10m*w10m *(w10m - u_ts)*dt) 








      nifa2d(i,j) = nifa2d(i,j) + dsrc/dt * 1000.0/dustmas/airmas

    endif 
    maxflux = MAX(maxflux,int(nifa2d(i,j)))
  enddo 
  enddo 

  enddo 


end subroutine bulk_dust_emis

END MODULE module_dust_emis
