MODULE module_bioemi_simple





  USE module_data_radm2
      INTEGER, PARAMETER ::  nlu = 25,  &
        iswater_temp = 16,isice_temp = 24
      REAL :: aefiso(nlu), aefmter(nlu), aefovoc(nlu), aef_n(nlu)
      CHARACTER (4),PARAMETER :: mminlu_loc = 'USGS'
      INTEGER :: ixxxlu(nlu)


    CONTAINS
      SUBROUTINE bio_emissions(id,ktau,dtstep,DX,                         &
               config_flags,                                              &
               gmt,julday,t_phy,moist,p8w,t8w,                            &
               e_bio,p_phy,chem,rho_phy,dz8w,ne_area,                     &
               ivgtyp,gsw,vegfra,rmol,ust,znt,xlat,xlong,z_at_w,          &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
  USE module_configure
  USE module_state_description
  IMPLICIT NONE
   INTEGER,      INTENT(IN   ) :: id,julday, ne_area,                     &
                                  ids,ide, jds,jde, kds,kde,              &
                                  ims,ime, jms,jme, kms,kme,              &
                                  its,ite, jts,jte, kts,kte
   INTEGER,      INTENT(IN   ) ::                                         &
                                  ktau
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),               &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, jms:jme, ne_area ),                          &
         INTENT(INOUT ) ::                               e_bio
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,              &
          INTENT(IN   ) ::                                                &
                                                      t_phy,              &
                                                      p_phy,              &
                                                      dz8w,               &
                                              t8w,p8w,z_at_w ,            &
                                                    rho_phy
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,              &
          INTENT(IN   ) ::                                                &
                                                     ivgtyp
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,              &
          INTENT(IN   ) ::                                                &
                                                     gsw,                 &
                                                  vegfra,                 &
                                                     rmol,                &
                                                     ust,                 &
                                                     xlat,                &
                                                     xlong,               &
                                                     znt
      REAL,      INTENT(IN   ) ::                                         &
                             dtstep,dx,gmt


   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags







      REAL :: emiss_bio(ne_area) 
      LOGICAL :: highnh3, rainflag, vegflag, wetflag
      CHARACTER (4) :: luse_typ


      REAL ::  clwchem,eiso,eisoc,emter,emterc,eovoc,eovocc,e_n,e_nn,  &
        pa,rad, rhchem, ta, ustar, vegfrac, vocsc, xtimin, z1,zntt
      INTEGER :: i,j,iland, iprt, iseason, n, nr, ipr,jpr,nvr



      INTRINSIC max, min

      luse_typ=mminlu_loc

      iseason=1
      if(julday.lt.90.or.julday.gt.270)then
        iseason=2
        CALL wrf_debug(100,'setting iseason in bio_emi to 2')
      endif
        
                         

                         

                         
      CALL bioemiin(iseason,luse_typ,vegflag) 
      do 100 j=jts,jte  
      do 100 i=its,ite  
      iland = ivgtyp(i,j)
      ta = t_phy(i,kts,j)      
      rad = gsw(i,j)
      vegfrac = vegfra(i,j)
      pa = .01*p_phy(i,kts,j)
      clwchem = moist(i,kts,j,p_qc)
      ustar = ust(i,j) 
      zntt = znt(i,j)                                                 
      z1 = z_at_w(i,kts+1,j)-z_at_w(i,kts,j)                          
                                                                      

      rainflag = .FALSE.                                              
      wetflag = .FALSE.                                               
      highnh3 = .FALSE.                                               
                                                                      
      if(moist(i,kts,j,p_qr).gt.0.)rainflag = .true.                  

                                                                      

      rhchem = MIN( 100.,100. * moist(i,kts,j,p_qv) / &               
               (3.80*exp(17.27*(t_phy(i,kts,j)-273.)/(t_phy(i,kts,j)-36.))/pa))     
      rhchem = max(rhchem,5.)
      if (rhchem >= 95.) wetflag = .true.                             

      if(chem(i,kts,j,p_nh3).gt.2.*chem(i,kts,j,p_so2))highnh3 = .true.
      iseason = 1                                                     

      emiss_bio=0.
      CALL biogen(iland,ta,rad,eiso,emter,eovoc,e_n,vocsc,eisoc,emterc,eovocc, &
        e_nn,pa,luse_typ,iseason,vegflag)








      CALL biosplit(iland,eiso,emter,eovoc,e_n,emiss_bio,ne_area,vegfrac, &
        config_flags, luse_typ,vegflag)




      DO n = 1, ne_area-2
        e_bio(i,j,n) = emiss_bio(n)

      END DO
 100  continue
END SUBROUTINE bio_emissions



      SUBROUTINE bioemiin(isn,mminlu,vegflag)





















































































        LOGICAL :: vegflag
        CHARACTER (4) :: mminlu
        INTEGER :: isn







        INTEGER :: sum




        IF (mminlu=='OLD ') THEN

          aefiso(1) = 0.

          aefiso(2) = 8.

          aefiso(3) = 0.

          aefiso(4) = 4400.

          aefiso(5) = 780.

          aefiso(6) = 5775.

          aefiso(7) = 0.

          aefiso(8) = 0.

          aefiso(9) = 0.

          aefiso(10) = 70.

          aefiso(11) = 0.

          aefiso(12) = 3100.

          aefiso(13) = 0
        END IF
        IF (mminlu=='USGS') THEN

          aefiso(1) = 0.

          aefiso(2) = 8.

          aefiso(3) = 8.

          aefiso(4) = 8.

          aefiso(5) = 4.

          aefiso(6) = 2204.

          aefiso(7) = 0.

          aefiso(8) = 0.

          aefiso(9) = 0.

          aefiso(10) = 0.

          aefiso(11) = 4400.

          aefiso(12) = 780.

          aefiso(13) = 4400.

          aefiso(14) = 780.

          aefiso(15) = 5775.

          aefiso(16) = 0.

          aefiso(17) = 0.

          aefiso(18) = 5775.

          aefiso(19) = 0.

          aefiso(20) = 70.

          aefiso(21) = 70.

          aefiso(22) = 70.

          aefiso(23) = 0.

          aefiso(24) = 0.

          aefiso(25) = 0.
        END IF
        IF (mminlu=='SiB ') THEN

          aefiso(1) = 4400.

          aefiso(2) = 4400.

          aefiso(3) = 4400.

          aefiso(4) = 780.

          aefiso(5) = 780.

          aefiso(6) = 0.

          aefiso(7) = 0.

          aefiso(8) = 0.

          aefiso(9) = 0.

          aefiso(10) = 0.

          aefiso(11) = 0.

          aefiso(12) = 8.

          aefiso(13) = 0.

          aefiso(14) = 0.

          aefiso(15) = 0.

          aefiso(16) = 0.

          aefiso(17) = 0.
        END IF



        IF (mminlu=='OLD ') THEN

          aefmter(1) = 0.

          aefmter(2) = 20.

          aefmter(3) = 20.

          aefmter(4) = 385.

          aefmter(5) = 1380.

          aefmter(6) = 1001.

          aefmter(7) = 0.

          aefmter(8) = 0.

          aefmter(9) = 0.

          aefmter(10) = 0.

          aefmter(11) = 0.

          aefmter(12) = 270.

          aefmter(13) = 0
        END IF
        IF (mminlu=='USGS') THEN

          aefmter(1) = 0.

          aefmter(2) = 20.

          aefmter(3) = 20.

          aefmter(4) = 20.

          aefmter(5) = 20.

          aefmter(6) = 202.5

          aefmter(7) = 20.

          aefmter(8) = 20.

          aefmter(9) = 20.

          aefmter(10) = 0

          aefmter(11) = 385.

          aefmter(12) = 1380.

          aefmter(13) = 385.

          aefmter(14) = 1380.

          aefmter(15) = 1001.

          aefmter(16) = 0.

          aefmter(17) = 0.

          aefmter(18) = 1001.

          aefmter(19) = 0.

          aefmter(20) = 0.

          aefmter(21) = 0.

          aefmter(22) = 0.

          aefmter(23) = 0.

          aefmter(24) = 0.

          aefmter(25) = 0.
        END IF
        IF (mminlu=='SiB ') THEN

          aefmter(1) = 385.

          aefmter(2) = 385.

          aefmter(3) = 385.

          aefmter(4) = 1380.

          aefmter(5) = 1380.

          aefmter(6) = 20.

          aefmter(7) = 20.

          aefmter(8) = 20.

          aefmter(9) = 20.

          aefmter(10) = 20.

          aefmter(11) = 0.

          aefmter(12) = 20.

          aefmter(13) = 0.

          aefmter(14) = 0.

          aefmter(15) = 0.

          aefmter(16) = 0.

          aefmter(17) = 0.
        END IF



        IF (mminlu=='OLD ') THEN

          aefovoc(1) = 0.

          aefovoc(2) = 12.

          aefovoc(3) = 80.

          aefovoc(4) = 715.

          aefovoc(5) = 840.

          aefovoc(6) = 924.

          aefovoc(7) = 0.

          aefovoc(8) = 0.

          aefovoc(9) = 0.

          aefovoc(10) = 0.

          aefovoc(11) = 0.

          aefovoc(12) = 0.

          aefovoc(13) = 0
        END IF
        IF (mminlu=='USGS') THEN

          aefovoc(1) = 0.

          aefovoc(2) = 12.

          aefovoc(3) = 12.

          aefovoc(4) = 12.

          aefovoc(5) = 46.

          aefovoc(6) = 363.5

          aefovoc(7) = 80.

          aefovoc(8) = 80.

          aefovoc(9) = 80.

          aefovoc(10) = 0

          aefovoc(11) = 715.

          aefovoc(12) = 840.

          aefovoc(13) = 715.

          aefovoc(14) = 840.

          aefovoc(15) = 924.

          aefovoc(16) = 0.

          aefovoc(17) = 0.

          aefovoc(18) = 924.

          aefovoc(19) = 0.

          aefovoc(20) = 0.

          aefovoc(21) = 0.

          aefovoc(22) = 0.

          aefovoc(23) = 0.

          aefovoc(24) = 0.

          aefovoc(25) = 0.
        END IF
        IF (mminlu=='SiB ') THEN

          aefovoc(1) = 715.

          aefovoc(2) = 715.

          aefovoc(3) = 715.

          aefovoc(4) = 840.

          aefovoc(5) = 840.

          aefovoc(6) = 80.

          aefovoc(7) = 80.

          aefovoc(8) = 80.

          aefovoc(9) = 80.

          aefovoc(10) = 80.

          aefovoc(11) = 0.

          aefovoc(12) = 12.

          aefovoc(13) = 0.

          aefovoc(14) = 0.

          aefovoc(15) = 0.

          aefovoc(16) = 0.

          aefovoc(17) = 0.
        END IF



        IF (mminlu=='OLD ') THEN

          aef_n(1) = 0.

          aef_n(2) = 9.

          aef_n(3) = 0.9

          aef_n(4) = 0.07

          aef_n(5) = 0.07

          aef_n(6) = 0.07

          aef_n(7) = 0.

          aef_n(8) = 0.

          aef_n(9) = 0.

          aef_n(10) = 0.

          aef_n(11) = 0.

          aef_n(12) = 1.78

          aef_n(13) = 0
        END IF
        IF (mminlu=='USGS') THEN

          aef_n(1) = 0.

          aef_n(2) = 9.

          aef_n(3) = 9.

          aef_n(4) = 9.

          aef_n(5) = 4.95

          aef_n(6) = 4.535

          aef_n(7) = 0.9

          aef_n(8) = 0.07

          aef_n(9) = 0.07

          aef_n(10) = 0.

          aef_n(11) = 0.07

          aef_n(12) = 0.07

          aef_n(13) = 0.07

          aef_n(14) = 0.07

          aef_n(15) = 0.07

          aef_n(16) = 0.

          aef_n(17) = 0.

          aef_n(18) = 0.07

          aef_n(19) = 0.

          aef_n(20) = 0.

          aef_n(21) = 0.

          aef_n(22) = 0.

          aef_n(23) = 0.

          aef_n(24) = 0.

          aef_n(25) = 0.
        END IF
        IF (mminlu=='SiB ') THEN

          aef_n(1) = 0.07

          aef_n(2) = 0.07

          aef_n(3) = 0.07

          aef_n(4) = 0.07

          aef_n(5) = 0.07

          aef_n(6) = 0.07

          aef_n(7) = 0.9

          aef_n(8) = 0.07

          aef_n(9) = 0.07

          aef_n(10) = 0.07

          aef_n(11) = 0.

          aef_n(12) = 9.

          aef_n(13) = 0.

          aef_n(14) = 0.

          aef_n(15) = 0.

          aef_n(16) = 0.

          aef_n(17) = 0.
        END IF















        IF (mminlu=='OLD ') THEN
          ixxxlu(1) = 1
          ixxxlu(2) = 2
          ixxxlu(3) = 3
          ixxxlu(4) = 4
          ixxxlu(5) = 5
          ixxxlu(6) = 5
          ixxxlu(7) = 0
          ixxxlu(8) = 6
          ixxxlu(9) = 1
          ixxxlu(10) = 6
          ixxxlu(11) = 0
          ixxxlu(12) = 4
          ixxxlu(13) = 6
        END IF
        IF (mminlu=='USGS') THEN
          ixxxlu(1) = 1
          ixxxlu(2) = 2
          ixxxlu(3) = 2
          ixxxlu(4) = 2
          ixxxlu(5) = 2
          ixxxlu(6) = 4
          ixxxlu(7) = 3
          ixxxlu(8) = 6
          ixxxlu(9) = 3
          ixxxlu(10) = 6
          ixxxlu(11) = 4
          ixxxlu(12) = 5
          ixxxlu(13) = 4
          ixxxlu(14) = 5
          ixxxlu(15) = 5
          ixxxlu(16) = 0
          ixxxlu(17) = 6
          ixxxlu(18) = 4
          ixxxlu(19) = 1
          ixxxlu(20) = 6
          ixxxlu(21) = 4
          ixxxlu(22) = 6
          ixxxlu(23) = 1
          ixxxlu(24) = 0
          ixxxlu(25) = 1
        END IF
        IF (mminlu=='SiB ') THEN
          ixxxlu(1) = 4
          ixxxlu(2) = 4
          ixxxlu(3) = 4
          ixxxlu(4) = 5
          ixxxlu(5) = 5
          ixxxlu(6) = 6
          ixxxlu(7) = 3
          ixxxlu(8) = 6
          ixxxlu(9) = 6
          ixxxlu(10) = 6
          ixxxlu(11) = 1
          ixxxlu(12) = 2
          ixxxlu(13) = 6
          ixxxlu(14) = 1
          ixxxlu(15) = 0
          ixxxlu(16) = 0
          ixxxlu(17) = 1
        END IF










        IF (mminlu=='OLD ') THEN

          IF (isn==2) THEN

            aefiso(2) = 0.

            aefiso(4) = 0.

            aefiso(6) = 5775./2.

            aefiso(10) = 0.

            aefmter(2) = 0.

            aefmter(4) = 0.

            aefmter(6) = 1001./2.

            aefovoc(2) = 0.

            aefovoc(4) = 0.

            aefovoc(6) = 924./2.
          END IF
        END IF

        IF (mminlu=='USGS') THEN

          sum = 0.





          IF (sum>1) THEN
            vegflag = .TRUE.
          ELSE
            vegflag = .FALSE.
          END IF

          IF (( .NOT. vegflag) .AND. (isn==2)) THEN




            aefiso(2) = 0.

            aefiso(3) = 0.

            aefiso(4) = 0.

            aefiso(5) = 0.

            aefiso(6) = 0.

            aefiso(11) = 0.

            aefiso(12) = 0.

            aefiso(15) = 5775./2.

            aefiso(18) = 5775./2.

            aefiso(20) = 0.

            aefiso(21) = 0.

            aefiso(22) = 0.

            aefmter(2) = 0.

            aefmter(3) = 0.

            aefmter(4) = 0.

            aefmter(5) = 10.

            aefmter(6) = 0.

            aefmter(11) = 0.

            aefmter(12) = 0.

            aefmter(15) = 1001./2.

            aefmter(18) = 1001./2.

            aefovoc(2) = 0.

            aefovoc(3) = 0.

            aefovoc(4) = 0.

            aefovoc(5) = 40.

            aefovoc(6) = 0.

            aefovoc(11) = 0.

            aefovoc(12) = 0.

            aefovoc(15) = 924./2.

            aefovoc(18) = 924./2.
          END IF
        END IF

        IF (mminlu=='SiB ') THEN

          IF (isn==2) THEN

            aefiso(1) = 0.

            aefiso(2) = 0.

            aefiso(3) = 0.

            aefiso(12) = 0.

            aefmter(1) = 0.

            aefmter(2) = 0.

            aefmter(3) = 0.

            aefmter(12) = 0.

            aefovoc(1) = 0.

            aefovoc(2) = 0.

            aefovoc(3) = 0.

            aefovoc(12) = 0.
          END IF
        END IF

      END SUBROUTINE bioemiin



      SUBROUTINE biogen(iland,ta,rad,eiso,emter,eovoc,e_n,vocsc,eisoc,emterc, &
          eovocc,e_nn,pa,mminlu,isn,vegflag)










































        REAL :: eiso, eisoc, emter, emterc, eovoc, eovocc, e_n, e_nn, pa, rad, &
          ta, vocsc
        INTEGER :: iland, isn
        LOGICAL :: vegflag
        CHARACTER (4) :: mminlu






        REAL :: alpha, beta, cl, cl1, coniso, conn, conovoc, conter, ct, ct1, &
          ct2, ecf_iso, ecf_mter, ecf_n, ecf_ovoc, par, r, rat, tm, ts, tsoil


        INTRINSIC exp, sqrt


        alpha = 0.0027

        cl1 = 1.066

        r = 8.314

        ct1 = 95000

        ct2 = 230000

        tm = 314.

        ts = 303.

        beta = 0.09












        IF ((ixxxlu(iland)==4) .OR. (ixxxlu(iland)==5)) THEN

          par = 2.0*rad






          cl = alpha*cl1*par/sqrt(1+alpha*alpha*par*par)
          ct = exp(ct1*(ta-ts)/(r*ts*ta))/(1+exp(ct2*(ta-tm)/(r*ts*ta)))

          ecf_iso = cl*ct

          ecf_mter = exp(beta*(ta-ts)) 
          ecf_ovoc = ecf_mter

          tsoil = 0.84*(ta-273.15) + 3.6
          ecf_n = exp(0.071*tsoil)

        END IF




        IF (ixxxlu(iland)==2) THEN
          ecf_iso = exp(0.1*(ta-30.-273.15)) 
          ecf_mter = ecf_iso
          ecf_ovoc = ecf_iso

          tsoil = 0.72*(ta-273.15) + 5.8
          ecf_n = exp(0.071*tsoil)
        END IF




        IF ((ixxxlu(iland)==3) .OR. (ixxxlu(iland)==6)) THEN
          ecf_iso = exp(0.1*(ta-30.-273.15)) 
          ecf_mter = ecf_iso
          ecf_ovoc = ecf_iso

          tsoil = 0.66*(ta-273.15) + 8.8
          ecf_n = exp(0.071*tsoil)
        END IF




        IF ((ixxxlu(iland)==1) .OR. (iland==iswater_temp) .OR. (iland==isice_temp)) THEN
          ecf_iso = 0.
          ecf_mter = 0.
          ecf_ovoc = 0.
          ecf_n = 0.
        END IF











        rat = ta/pa



        coniso = rat*2.3095E-5
        eisoc = aefiso(iland)*ecf_iso
        eiso = coniso*eisoc




        conter = rat*1.1548E-5
        emterc = aefmter(iland)*ecf_mter
        emter = conter*emterc






        conovoc = rat*1.4435E-5
        eovocc = aefovoc(iland)*ecf_ovoc
        eovoc = conovoc*eovocc



        vocsc = eisoc + emterc + eovocc











        conn = rat*3.5633E-4
        e_nn = aef_n(iland)*ecf_n
        e_n = conn*e_nn


      END SUBROUTINE biogen



      SUBROUTINE biosplit(iland,eiso,emter,eovoc,e_n,emiss_bio,ne_area, &
          vegfrc, &
          config_flags, mminlu,vegflag)




















































  USE module_configure
  USE module_state_description

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags



        REAL :: eiso, emter, eovoc, e_n, vegfrc
        INTEGER :: iland, ne_area




        REAL :: emiss_bio(ne_area)



        LOGICAL :: vegflag
        CHARACTER (4) :: mminlu



        IF ((mminlu=='USGS') .AND. (vegflag)) THEN
          eiso = eiso*vegfrc/100.
          emter = emter*vegfrc/100.
          eovoc = eovoc*vegfrc/100.
        END IF




        emiss_bio(liso) = eiso
        emiss_bio(lno) = emiss_bio(lno) + e_n

      if (config_flags%chem_opt == CB05_SORG_AQ_KPP .OR. &
          config_flags%chem_opt == CB05_SORG_VBS_AQ_KPP ) then

      emiss_bio(ltpan) = emter 




        IF (ixxxlu(iland)==2) THEN
          emiss_bio(lhc5) = emiss_bio(lhc5) + 0.16*eovoc
          emiss_bio(lhc8) = emiss_bio(lhc8) + 0.27*eovoc
          emiss_bio(lolt) = emiss_bio(lolt) + 0.05*eovoc
          emiss_bio(loli) = emiss_bio(loli) + 0.37*eovoc
          emiss_bio(lket) = emiss_bio(lket) + 0.03*eovoc
          emiss_bio(lald) = emiss_bio(lald) + 0.12*eovoc
        END IF




        IF (ixxxlu(iland)==3) THEN
          emiss_bio(lhc5) = emiss_bio(lhc5) + 0.09*eovoc
          emiss_bio(lolt) = emiss_bio(lolt) + 0.07*eovoc
          emiss_bio(loli) = emiss_bio(loli) + 0.51*eovoc
          emiss_bio(lket) = emiss_bio(lket) + 0.15*eovoc
          emiss_bio(lald) = emiss_bio(lald) + 0.18*eovoc
        END IF




        IF (ixxxlu(iland)==4) THEN
          emiss_bio(lhcho) = emiss_bio(lhcho) + 0.19*eovoc
          emiss_bio(lald) = emiss_bio(lald) + 0.13*eovoc
          emiss_bio(lxyl) = emiss_bio(lxyl) + 0.04*emter
          emiss_bio(lhc5) = emiss_bio(lhc5) + 0.03*eovoc
          emiss_bio(loli) = emiss_bio(loli) + 0.07*eovoc
          emiss_bio(lora1) = emiss_bio(lora1) + 0.23*eovoc
          emiss_bio(lora2) = emiss_bio(lora2) + 0.35*eovoc
        END IF





        IF (ixxxlu(iland)==5) THEN
          emiss_bio(lhcho) = emiss_bio(lhcho) + 0.04*eovoc
          emiss_bio(lhcho) = emiss_bio(lhcho) + 0.04*eovoc
          emiss_bio(lald) = emiss_bio(lald) + 0.14*eovoc
          emiss_bio(lhc3) = emiss_bio(lhc3) + 0.07*eovoc
          emiss_bio(lhc5) = emiss_bio(lhc5) + 0.07*eovoc
          emiss_bio(lolt) = emiss_bio(lolt) + 0.07*eovoc
          emiss_bio(loli) = emiss_bio(loli) + 0.50*eovoc
          emiss_bio(lket) = emiss_bio(lket) + 0.03*eovoc
          emiss_bio(lora1) = emiss_bio(lora1) + 0.03*eovoc
          emiss_bio(lora2) = emiss_bio(lora2) + 0.05*eovoc
        END IF

      else
      




        IF (ixxxlu(iland)==2) THEN
          emiss_bio(loli) = emiss_bio(loli) + 0.80*emter
          emiss_bio(liso) = emiss_bio(liso) + 0.20*emter
          emiss_bio(lhc5) = emiss_bio(lhc5) + 0.16*eovoc
          emiss_bio(lhc8) = emiss_bio(lhc8) + 0.27*eovoc
          emiss_bio(lolt) = emiss_bio(lolt) + 0.05*eovoc
          emiss_bio(loli) = emiss_bio(loli) + 0.37*eovoc
          emiss_bio(lket) = emiss_bio(lket) + 0.03*eovoc
          emiss_bio(lald) = emiss_bio(lald) + 0.12*eovoc
        END IF




        IF (ixxxlu(iland)==3) THEN
          emiss_bio(loli) = emiss_bio(loli) + 0.98*emter
          emiss_bio(liso) = emiss_bio(liso) + 0.02*emter
          emiss_bio(lhc5) = emiss_bio(lhc5) + 0.09*eovoc
          emiss_bio(lolt) = emiss_bio(lolt) + 0.07*eovoc
          emiss_bio(loli) = emiss_bio(loli) + 0.51*eovoc
          emiss_bio(lket) = emiss_bio(lket) + 0.15*eovoc
          emiss_bio(lald) = emiss_bio(lald) + 0.18*eovoc
        END IF




        IF (ixxxlu(iland)==4) THEN
          emiss_bio(loli) = emiss_bio(loli) + 0.94*emter
          emiss_bio(liso) = emiss_bio(liso) + 0.02*emter
          emiss_bio(lhcho) = emiss_bio(lhcho) + 0.19*eovoc
          emiss_bio(lald) = emiss_bio(lald) + 0.13*eovoc
          emiss_bio(lxyl) = emiss_bio(lxyl) + 0.04*emter
          emiss_bio(lhc5) = emiss_bio(lhc5) + 0.03*eovoc
          emiss_bio(loli) = emiss_bio(loli) + 0.07*eovoc
          emiss_bio(lora1) = emiss_bio(lora1) + 0.23*eovoc
          emiss_bio(lora2) = emiss_bio(lora2) + 0.35*eovoc
        END IF





        IF (ixxxlu(iland)==5) THEN
          emiss_bio(loli) = emiss_bio(loli) + 0.85*emter
          emiss_bio(liso) = emiss_bio(liso) + 0.15*emter
          emiss_bio(lhcho) = emiss_bio(lhcho) + 0.04*eovoc
          emiss_bio(lald) = emiss_bio(lald) + 0.14*eovoc
          emiss_bio(lhc3) = emiss_bio(lhc3) + 0.07*eovoc
          emiss_bio(lhc5) = emiss_bio(lhc5) + 0.07*eovoc
          emiss_bio(lolt) = emiss_bio(lolt) + 0.07*eovoc
          emiss_bio(loli) = emiss_bio(loli) + 0.50*eovoc
          emiss_bio(lket) = emiss_bio(lket) + 0.03*eovoc
          emiss_bio(lora1) = emiss_bio(lora1) + 0.03*eovoc
          emiss_bio(lora2) = emiss_bio(lora2) + 0.05*eovoc
        END IF




        IF ((mminlu=='OLD ') .AND. (iland==12)) THEN
          emiss_bio(loli) = emiss_bio(loli) + emter
        END IF

      end if

      END SUBROUTINE biosplit

    END MODULE module_bioemi_simple
