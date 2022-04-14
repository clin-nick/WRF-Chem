















 MODULE module_lightning_nox_decaria

 IMPLICIT NONE

 CONTAINS












 SUBROUTINE lightning_nox_decaria ( &
                          
                            dx, dy, xland, ht, t, rho, z, p,      &
                            ic_flashrate, cg_flashrate,           & 
                          
                            refl,                                 &
                          
                            N_IC, N_CG,                           &
                            ltng_temp_upper,ltng_temp_lower,      &
                            cellcount_method,                     &
                          
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          
                            lnox_ic_tend, lnox_cg_tend            & 
                          )


 USE module_state_description


 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_max_real, wrf_dm_min_real, wrf_dm_sum_real


 USE module_lightning_driver, only: countCells

 IMPLICIT NONE



 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: t, rho, z, p
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: ic_flashrate  , cg_flashrate 



 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl


 REAL,    INTENT(IN   )    ::       N_IC, N_CG
 REAL,    INTENT(IN   )    ::       ltng_temp_lower, ltng_temp_upper
 INTEGER, INTENT(IN   )    ::       cellcount_method


 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe


 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) :: lnox_ic_tend,lnox_cg_tend


 INTEGER :: i,k,j
 INTEGER :: ktop,kbtm,kupper,klower
 REAL :: ic_fr, cg_fr, delta 
 REAL :: reflmax, cellmax
 REAL :: term2, B, B_denom
 CHARACTER (LEN=250) :: message
 REAL, DIMENSION( kps:kpe ) :: cellcount
 REAL, DIMENSION( kps:kpe ) :: z_average, t_average, p_average, rho_average, conv
 REAL, DIMENSION( kps:kpe ) :: fd, fd2, dz 

 REAL, PARAMETER :: refl_threshold = 20.



 lnox_ic_tend (ips:ipe,kps:kpe,jps:jpe ) = 0.
 lnox_cg_tend (ips:ipe,kps:kpe,jps:jpe ) = 0.


 CALL countCells( &
      
        refl, refl_threshold, cellcount_method,     &
      
        ids, ide, jds, jde, kds, kde,              &
        ims, ime, jms, jme, kms, kme,              &
        ips, ipe, jps, jpe, kps, kpe,              &
      
        cellcount )


 ic_fr = sum(ic_flashrate(ips:ipe,jps:jpe))
 cg_fr = sum(cg_flashrate(ips:ipe,jps:jpe))
 if ( cellcount_method .eq. 2 ) then
    ic_fr = wrf_dm_sum_real(ic_fr)
    cg_fr = wrf_dm_sum_real(cg_fr)
 ENDIF
 reflmax = maxval(refl(ips:ipe,kps:kpe,jps:jpe))
 cellmax = maxval(cellcount(kps:kpe))
 WRITE(message, * ) ' LNOx tracer: max_refl, max_cellcount, ic_fr = ',  reflmax, cellmax, ic_fr
 CALL wrf_debug ( 100, message )




 CALL horizontalAverage( z( ips:ipe,kps:kpe,jps:jpe ), ips, ipe, kps, kpe, jps, jpe, z_average )
 CALL horizontalAverage( t( ips:ipe,kps:kpe,jps:jpe ), ips, ipe, kps, kpe, jps, jpe, t_average )
 CALL horizontalAverage( p( ips:ipe,kps:kpe,jps:jpe ), ips, ipe, kps, kpe, jps, jpe, p_average )
 CALL horizontalAverage( rho( ips:ipe,kps:kpe,jps:jpe ), ips, ipe, kps, kpe, jps, jpe, rho_average )




 conv(kps:kpe) = 8.314 *t_average(kps:kpe) / (dx * dy)                


 CALL  kfind ( cellcount, t_average,            &
               ltng_temp_upper,ltng_temp_lower, cellcount_method, &
                ips, ipe, jps, jpe, kps, kpe,              &
              
                ktop,kbtm,kupper,klower )



 IF (( ic_fr > 0 ) .and. (( ktop > klower ) .and. (kbtm < klower) ) )THEN
   call bellcurve(kbtm,ktop,klower,z_average, kps,kpe, fd, dz)
   if (ktop .gt. kupper) then
     call bellcurve(kbtm,ktop,kupper,z_average, kps,kpe, fd2, dz)
     fd(kbtm:ktop) = 0.5*( fd(kbtm:ktop) + fd2(kbtm:ktop) )         
   endif

   B_denom = DOT_PRODUCT( fd(kbtm:ktop),p_average(kbtm:ktop) )      


   DO k=kbtm,ktop
     if ( cellcount(k) .gt. 0. ) THEN

        
        
        
        

        delta = (ic_fr * N_IC / cellcount(k)) * fd(k) / B_denom * conv(k)/dz(k) * 1E6




        where(refl(ips:ipe,k,jps:jpe) .gt. refl_threshold )
          lnox_ic_tend(ips:ipe,k,jps:jpe) = delta
        endwhere
     ENDIF
   ENDDO

 ENDIF 




 IF ((cg_fr > 0 ) .and. (( ktop > klower ) .and. (kbtm < klower) ) ) THEN
   call bellcurve(kps,ktop,klower,z_average, kps,kpe, fd, dz)


   B_denom = DOT_PRODUCT( fd(kbtm:ktop),p_average(kbtm:ktop) )      

   k = ktop

   DO WHILE (k .ge. kps)
     IF (cellcount(k) .gt. 0) THEN


      delta = (cg_fr * N_CG / cellcount(k)) * fd(k) / B_denom * conv(k)/dz(k) * 1E6




       where( refl(ips:ipe,k,jps:jpe) .gt. refl_threshold )
          lnox_cg_tend(ips:ipe,k,jps:jpe) = delta
       endwhere

     ENDIF

     k = k - 1      

   ENDDO

 ENDIF

 END SUBROUTINE lightning_nox_decaria















 SUBROUTINE bellcurve ( k_min, k_max, k_mu, z, kps,kpe, f, dz )


 IMPLICIT NONE
 INTEGER,                      INTENT(IN   ) :: k_min, k_max, k_mu
 REAL,   DIMENSION( kps:kpe ), INTENT(IN   ) :: z       
 INTEGER,                      INTENT(IN   ) :: kps,kpe

 REAL,   DIMENSION( kps:kpe ), INTENT(  OUT) :: f, dz

 INTEGER :: i,j,k
 REAL, DIMENSION( kps:kpe ) :: ex
 REAL :: sigma, z_mu, cuml_f_dist
 REAL, PARAMETER :: two_pi = 6.2831854



 f(kps:kpe) = 0.
 z_mu = z(k_mu)
 sigma = AMIN1(z(k_max)-z_mu,z_mu-z(k_min))/3.0

 
 ex(k_min:k_max) = (z(k_min:k_max)-z_mu)/sigma
 
 
 f(k_min:k_max) = (1.0/(sqrt(two_pi)*sigma))*exp(-ex(k_min:k_max)*ex(k_min:k_max)/2.0)




 dz(kps) = (z(kps+1) - z(kps))*.5             
 dz(kpe) = (z(kpe) - z(kpe-1))*.5             
 DO k=kps+1,kpe-1

   dz(k) = (z(k+1) - z(k-1))*.5
 ENDDO

 
 cuml_f_dist = DOT_PRODUCT(dz(k_min:k_max),f(k_min:k_max))
 f(k_min:k_max) = f(k_min:k_max)*dz(k_min:k_max)/cuml_f_dist

 END SUBROUTINE bellcurve









 SUBROUTINE kfind ( &
              
                cellcount, t,                         &
              
                ltng_temp_upper,ltng_temp_lower,      &
                cellcount_method,                     &
              
                ips, ipe, jps, jpe, kps, kpe,          &
              
                ktop,kbtm,kupper,klower               &
            )


 USE module_state_description


 USE module_model_constants

 USE module_dm, only: wrf_dm_max_real, wrf_dm_min_real, wrf_dm_sum_real

 IMPLICIT NONE



 REAL, DIMENSION( kps:kpe ), INTENT(IN   ) :: cellcount
 REAL, DIMENSION( kps:kpe ), INTENT(IN   ) :: t


 REAL,    INTENT(IN   )    ::       ltng_temp_lower, ltng_temp_upper
 INTEGER, INTENT(IN   )    ::       cellcount_method


 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe


 INTEGER, INTENT(  OUT)    ::       ktop,kbtm,kupper,klower


 CHARACTER (LEN=250) :: message
 REAL    :: ktop_r, kbtm_r, kupper_r, klower_r
 INTEGER :: k


 ktop = kps
 kbtm = kps
 kupper = kps
 klower = kps

 
 k = kpe
 DO WHILE ( cellcount(k) .eq. 0 .and. k .gt. kps)
  k = k-1
 ENDDO
 ktop = k

 
 k = kps
 DO WHILE( cellcount(k) .eq. 0 .and. k .le. ktop )
  k = k+1
 ENDDO
 kbtm = k


 
 k = kps
 DO WHILE ( t(k) .gt. ltng_temp_lower + 273.15 .and. k .lt. kpe )
   k = k + 1
 ENDDO
 klower = k

 DO WHILE ( t(k) .gt. ltng_temp_upper + 273.15 .and. k .lt. kpe )
   k = k + 1
 ENDDO
 kupper = k

 WRITE(message, * ) ' LNOx_driver: kbtm, ktop, klower, kupper = ', kbtm, ktop, klower, kupper
 CALL wrf_debug ( 100, message )
 
 IF ( cellcount_method .eq. 2 ) THEN
   kbtm_r = real(kbtm)
   ktop_r = real(ktop)
   klower_r = real(klower)
   kupper_r = real(kupper)
   kbtm = nint(wrf_dm_min_real(kbtm_r))
   ktop = nint(wrf_dm_max_real(ktop_r))
   klower = nint(wrf_dm_max_real(klower_r))
   kupper = nint(wrf_dm_max_real(kupper_r))

   WRITE(message, * ) ' lightning_driver: kbtm, ktop, klower, kupper = ', kbtm, ktop, klower, kupper
   CALL wrf_debug ( 100, message )
 endif

 END SUBROUTINE kfind




 SUBROUTINE horizontalAverage( array3D, ips, ipe, kps, kpe, jps, jpe, array1D )

 IMPLICIT NONE
 REAL, DIMENSION(ips:ipe,kps:kpe,jps:jpe), INTENT(IN) :: array3D
 INTEGER, INTENT(IN) :: ips,ipe,kps,kpe,jps,jpe

 INTEGER :: k
 REAL, DIMENSION(kps:kpe), INTENT(OUT) :: array1D

 DO k=kps,kpe
   array1D(k) = sum(array3D(ips:ipe,k,jps:jpe))/((ipe-ips+1)*(jpe-jps+1))
 ENDDO
    
 END SUBROUTINE horizontalAverage

 END MODULE module_lightning_nox_decaria
