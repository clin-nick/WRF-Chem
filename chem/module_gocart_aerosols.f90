MODULE module_gocart_aerosols
  USE module_model_constants, only: mwdry
    INTEGER, PARAMETER ::NBC1=1, NOC1=2, NBC2=3, NOC2=4

CONTAINS
  subroutine gocart_aerosols_driver(ktau,dt,config_flags,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,dx,g,         &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_configure
  USE module_state_description
  IMPLICIT NONE
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: ktau,                     &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                     t_phy,               &
                                            dz8w,p8w,      &
                                              rho_phy

  REAL, INTENT(IN   ) :: dt,dx,g
  integer :: ndt1,nmx,i,j,k,imx,jmx,lmx
  real*8, DIMENSION (1,1,1) :: tmp,airden,airmas
  REAL*8    :: chmlos(1,1,1,4)
  REAL*8    :: bchmlos(1,1,4)
  REAL*8 :: pc2(1,1,1,2)
  REAL*8 :: tc(4),tt1,tt2
  real, parameter :: mw_c = 12.

        INTRINSIC abs

       imx=1
       jmx=1
       lmx=1
       nmx=4
       ndt1=ifix(dt)


       do j=jts,jte
       do k=kts,kte-1
       do i=its,ite
          airmas(1,1,1)=-(p8w(i,k+1,j)-p8w(i,k,j))*dx*dx/g
          pc2(1,1,1,1)=0.
          pc2(1,1,1,2)=0.
          tc(1)=chem(i,k,j,p_bc1)/mw_c*mwdry*1.d-9
          tc(2)=chem(i,k,j,p_oc1)/mw_c*mwdry*1.d-9
          tc(3)=chem(i,k,j,p_bc2)/mw_c*mwdry*1.d-9
          tc(4)=chem(i,k,j,p_oc2)/mw_c*mwdry*1.d-9
          tt1=tc(3)

         CALL chem_1(imx,jmx,lmx, nmx, ndt1, airmas, tc, &
              chmlos, bchmlos, pc2)
         CALL chem_2(imx,jmx,lmx, nmx, ndt1, airmas, tc, pc2)




          tt2=tc(3)-tt1
          chem(i,k,j,p_bc1)=tc(1)/mwdry*mw_c*1.e9
          chem(i,k,j,p_oc1)=tc(2)/mwdry*mw_c*1.e9
          chem(i,k,j,p_bc2)=tc(3)/mwdry*mw_c*1.e9
          chem(i,k,j,p_oc2)=(tc(4)+8.*tt2)/mwdry*mw_c*1.e9

       enddo
       enddo
       enddo
end subroutine gocart_aerosols_driver
       subroutine sum_pm_gocart (                                      &
            alt, chem,pm2_5_dry, pm2_5_dry_ec, pm10,                   &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )
   USE module_configure
   USE module_state_description
   USE module_data_gocartchem, only: nh4_mfac,oc_mfac
   IMPLICIT NONE
   REAL, PARAMETER :: mwso4 = 96.0576
   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                          &
         INTENT(INOUT ) :: pm2_5_dry, pm2_5_dry_ec, pm10     
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                          &
         INTENT(IN ) :: alt
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(IN ) ::                                   chem
   real d_2_5,s_2_5,d_10,sulfate
   integer i,j,k,ii,jj,n
    d_2_5=0.380       
    s_2_5=0.834       
    d_10=0.737        




      pm2_5_dry(its:ite, kts:kte, jts:jte)    = 0.
      pm10(its:ite, kts:kte, jts:jte)    = 0.
      pm2_5_dry_ec(its:ite, kts:kte, jts:jte) = 0.
      do j=jts,jte                    
         jj=min(jde-1,j)              
         do k=kts,kte
         do i=its,ite
            ii=min(ide-1,i)
            sulfate=chem(ii,k,jj,p_sulf)*mwso4/mwdry*1.e3
            do n=p_p25,p_dust_1
               pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j)+chem(ii,k,jj,n)        
            enddo
            pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j)+chem(ii,k,jj,p_dust_2)*d_2_5     &
                                            +chem(ii,k,jj,p_seas_1)           &        
                                            +chem(ii,k,jj,p_seas_2)*s_2_5     &        
                                            +sulfate*nh4_mfac                 &
                                            +(chem(ii,k,jj,p_oc1)+chem(ii,k,jj,p_oc2))*(oc_mfac-1.)
   
            
            pm2_5_dry(i,k,j)    = pm2_5_dry(i,k,j) / alt(ii,k,jj)
      enddo
      enddo
      enddo
      do j=jts,jte
         jj=min(jde-1,j)
         do k=kts,kte
            do i=its,ite
               ii=min(ide-1,i)
               sulfate=chem(ii,k,jj,p_sulf)*mwso4/mwdry*1.e3
               do n=p_p25,p_dust_3
                  pm10(i,k,j) = pm10(i,k,j)+chem(ii,k,jj,n)        
               enddo
               do n=p_seas_1,p_seas_3
                  pm10(i,k,j) = pm10(i,k,j)+chem(ii,k,jj,n)        
               enddo
               pm10(i,k,j) = pm10(i,k,j) + sulfate*nh4_mfac        &
                            +chem(ii,k,jj,p_dust_4)*d_10           &
                            +chem(ii,k,jj,p_p10)                &
                            +(chem(ii,k,jj,p_oc1)+chem(ii,k,jj,p_oc2))*(oc_mfac-1.)
               pm10(i,k,j) = pm10(i,k,j)/ alt(ii,k,jj)
            enddo
         enddo
      enddo

end subroutine sum_pm_gocart

SUBROUTINE chem_1(imx,jmx,lmx, nmx, &
     ndt1, airm, tc, chmlos, bchmlos, pc2)













  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: lmx, nmx,imx,jmx, ndt1
  REAL*8,    INTENT(IN)    :: airm(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(INOUT) :: chmlos(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(INOUT) :: bchmlos(imx,jmx,nmx)
  REAL*8,    INTENT(OUT)   :: pc2(imx,jmx,lmx,2)

  REAL*8 :: r1, c0, r2, rkt, c1
  INTEGER :: np, n, i, j, l

  

  r1 = 4.63D-6
  
  DO n = 1,nmx
     IF (n == NBC1 .OR. n == NOC1) THEN
        IF (n == NBC1) np = 1
        IF (n == NOC1) np = 2
        DO l = 1,lmx
           DO j = 1,jmx
              DO i = 1,imx
                 
                 c0 = tc(i,j,l,n)
                 r2 = 0.0 
                 rkt = (r1 + r2) * REAL(ndt1)
                 
                 c1 = c0 *  EXP(-rkt)
                 c1 = MAX(c1, 1.0D-32)
                 tc(i,j,l,n) = c1
                 
                 pc2(i,j,l,np) = (c0 - c1) * r1/(r1 + r2)
                 
                 

                 chmlos(i,j,l,n) = chmlos(i,j,l,n) + pc2(i,j,l,np)*airm(i,j,l)
                 
              END DO
           END DO
        END DO

        DO j = 1,jmx
           DO i = 1,imx
              bchmlos(i,j,n) = bchmlos(i,j,n) + SUM(chmlos(i,j,:,n))
           END DO
        END DO

     END IF
  END DO
  
END SUBROUTINE chem_1

SUBROUTINE chem_2(imx,jmx,lmx, nmx, &
                  ndt1, airm, tc,  pc2)








  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: lmx,imx,jmx, nmx, ndt1
  REAL*8,    INTENT(IN)    :: airm(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(IN)    :: pc2(imx,jmx,lmx,2)

  INTEGER :: np, n, i, j, l
  REAL*8  :: c0, pp, rkt, c1

  
  
  DO n = 1,nmx
     IF (n == NBC2 .OR. n == NOC2) THEN
        IF (n == NBC2) np = 1
        IF (n == NOC2) np = 2

        DO l = 1,lmx
           DO j = 1,jmx
              DO i = 1,imx
                 
                 c0 = tc(i,j,l,n)

                 pp = pc2(i,j,l,np)
                 c1 = c0 + pp
                 
                 c1 = MAX(c1, 1.0D-32)
                 tc(i,j,l,n) = c1
                 
                 
              END DO
           END DO
        END DO
     END IF
  END DO
  
END SUBROUTINE chem_2

END MODULE module_gocart_aerosols
