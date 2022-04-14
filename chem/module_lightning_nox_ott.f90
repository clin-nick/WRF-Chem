
















 MODULE module_lightning_nox_ott

 IMPLICIT NONE

 CONTAINS















 SUBROUTINE lightning_nox_ott ( &
                          
                            dx, dy, xlat, xland, ht, rho, z,      &
                            ic_flashrate, cg_flashrate,           & 
                          
                            N_IC, N_CG,                           &
                          
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte,         &
                          
                            lnox_total_tend                       & 
                          )


 USE module_state_description


 USE module_model_constants
 USE module_wrf_error

 IMPLICIT NONE



 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xlat, xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: rho, z
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: ic_flashrate  , cg_flashrate 


 REAL,    INTENT(IN   )    ::       N_IC, N_CG


 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       its,ite, jts,jte, kts,kte


 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) :: lnox_total_tend


 REAL, PARAMETER :: subtrop_midlat = 35.
 REAL, PARAMETER :: trop_subtrop = 20.


  INTEGER,                  PARAMETER :: vds = 0
  INTEGER,                  PARAMETER :: vde = 16
                             
  REAL, DIMENSION(vds:vde), PARAMETER :: &
       ott_subtrop(vde+1) = (/ .010,.020,.039,.058,.077,.093,.105,.110,.110,.104,.092,.075,.055,.034,.015,.002,.000 /)
  REAL, DIMENSION(vds:vde), PARAMETER :: &
       ott_midlat(vde+1)  = (/ .024,.050,.074,.093,.106,.114,.115,.110,.099,.083,.063,.042,.022,.005,.000,.000,.000 /)
  REAL, DIMENSION(vds:vde), PARAMETER :: & 
       ott_trpcon(vde+1)  = (/ .002,.005,.006,.014,.027,.040,.050,.062,.086,.103,.116,.124,.127,.124,.076,.030,.008 /)
  REAL, DIMENSION(vds:vde), PARAMETER :: & 
       ott_trpmar(vde+1)  = (/ .006,.015,.029,.043,.054,.067,.066,.085,.096,.102,.105,.102,.082,.065,.045,.022,.005 /)


 INTEGER :: i,k,j
 INTEGER :: v       
 REAL    :: total_NO, mass_of_air, dA
 REAL, DIMENSION( kts:kte ):: zkm     
 REAL, DIMENSION( vds:vde ):: NOperkm 



 lnox_total_tend(its:ite,kts:kte,jts:jte) = 0.
 dA = dx * dy

 DO J=jts,jte
   DO I=its,ite

     
     total_NO = ic_flashrate(I,J)*N_IC + cg_flashrate(I,J)*N_CG
     IF ( total_NO .eq. 0. ) CONTINUE

     
     IF ( xlat(I,J) .gt. subtrop_midlat ) THEN
       NOperkm(:) = ott_midlat(:) * total_NO
     ELSE IF ( xlat(I,J) .gt. trop_subtrop ) THEN
       NOperkm(:) = ott_subtrop(:) * total_NO
     ELSE IF ( xland(I,J) .gt. 1.5 ) THEN
       NOperkm(:) = ott_trpcon(:) * total_NO
     ELSE
       NOperkm(:) = ott_trpmar(:) * total_NO
     ENDIF

     
     
     
     
     k = kts
     zkm(kts:kte) = ( z(i,kts:kte,j) - ht(i,j) ) / 1000.
     v = MAX( vds, int(zkm(k)) )
     DO WHILE ( (v .le. vde) .and. (k .le. kte) )
       mass_of_air = rho(i,k,j) * 1E3 * dA / .02897       
       lnox_total_tend(i,k,j) = NOperkm(v)/mass_of_air * 1E6  

       k = k + 1
       IF ( int(zkm(k)) .gt. v ) v = int( zkm(k) )
     ENDDO
   ENDDO
 ENDDO

 END SUBROUTINE lightning_nox_ott


 END MODULE module_lightning_nox_ott
