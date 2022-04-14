





































MODULE module_input_chem_bioemiss

    USE module_io_domain
    USE module_domain
    USE module_driver_constants
    USE module_state_description
    USE module_configure
    USE module_date_time
    USE module_wrf_error
    USE module_timing
    USE module_data_radm2
    USE module_aerosols_sorgam
    USE module_get_file_names


CONTAINS


SUBROUTINE input_ext_chem_beis3_file (grid)


   IMPLICIT NONE

   TYPE(domain)           ::  grid

   INTEGER ::  i,j,n,numfil,status,system

   INTEGER :: ids, ide, jds, jde, kds, kde,    &
              ims, ime, jms, jme, kms, kme,    &
              ips, ipe, jps, jpe, kps, kpe

   REAL, ALLOCATABLE, DIMENSION(:,:) :: emiss




      PARAMETER(numfil=19)

   CHARACTER (LEN=80) :: message

   TYPE (grid_config_rec_type)              :: config_flags












      CHARACTER*100 onefil
      CHARACTER*12 emfil(numfil)
      DATA emfil/'ISO','OLI','API','LIM','XYL','HC3','ETE','OLT',  &
        'KET','ALD','HCHO','ETH','ORA2','CO','NR',                 &
        'NOAG_GROW','NOAG_NONGROW','NONONAG','ISOP'/



       
       CALL get_ijk_from_grid (  grid ,                        &
                                 ids, ide, jds, jde, kds, kde,    &
                                 ims, ime, jms, jme, kms, kme,    &
                                 ips, ipe, jps, jpe, kps, kpe    )

     WRITE( message , FMT='(A,4I5)' ) ' DIMS: ',ids,ide-1,jds,jde-1
     CALL  wrf_message ( message )

     ALLOCATE( emiss(ids:ide-1,jds:jde-1) )



      DO n=1,numfil


       status=system('rm -f scratem*')


       IF(n.LE.15)THEN 
        onefil='../../run/BIOREF_'//             &
         TRIM(ADJUSTL(emfil(n)))//'.gz'

       ELSE IF(n.GE.16.AND.n.LE.18)THEN 
        onefil='../../run/AVG_'//                &
         TRIM(ADJUSTL(emfil(n)))//'.gz'

       ELSE
        onefil='../../run/LAI_'//                &
         TRIM(ADJUSTL(emfil(n)))//'S.gz'
       ENDIF


       status=system('cp '//TRIM(ADJUSTL(onefil))//' scratem.gz')


       status=system('gunzip scratem')


       OPEN(26,FILE='scratem',FORM='FORMATTED')
       IF(n.EQ. 1) then
             READ(26,'(12E9.2)') emiss
             grid%sebio_iso(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ. 2)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_oli(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ. 3)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_api(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ. 4)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_lim(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ. 5)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_xyl(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ. 6)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_hc3(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ. 7)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_ete(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ. 8)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_olt(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ. 9)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_ket(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.10)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_ald(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.11)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_hcho(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.12)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_eth(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.13)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_ora2(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.14)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_co(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.15)then
              READ(26,'(12E9.2)') emiss
              grid%sebio_nr(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.16)then
              READ(26,'(12E9.2)') emiss
              grid%noag_grow(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.17)then
              READ(26,'(12E9.2)') emiss
              grid%noag_nongrow(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.18)then
              READ(26,'(12E9.2)') emiss
              grid%nononag(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       IF(n.EQ.19)then
              READ(26,'(12E9.2)') emiss
              grid%slai(ids:ide-1,jds:jde-1) = emiss
       ENDIF
       CLOSE(26)

      ENDDO


    DEALLOCATE( emiss )

END SUBROUTINE input_ext_chem_beis3_file 



SUBROUTINE input_ext_chem_megan2_file (grid)

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



   IMPLICIT NONE          

   TYPE(domain)           ::  grid

   INTEGER ::  i,j,v,status,system, itmp, jtmp

   INTEGER :: ids, ide, jds, jde, kds, kde,    &
              ims, ime, jms, jme, kms, kme,    &
              ips, ipe, jps, jpe, kps, kpe

   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: emiss

   CHARACTER (LEN=80) :: message

   TYPE (grid_config_rec_type)              :: config_flags


   
   
   integer, parameter :: n_mgnin = 41
   integer, parameter ::        & 
        &  mgnin_isop     =  1  & 
        & ,mgnin_lai01    =  2  & 
        & ,mgnin_lai02    =  3  & 
        & ,mgnin_lai03    =  4  & 
        & ,mgnin_lai04    =  5  & 
        & ,mgnin_lai05    =  6  & 
        & ,mgnin_lai06    =  7  & 
        & ,mgnin_lai07    =  8  & 
        & ,mgnin_lai08    =  9  & 
        & ,mgnin_lai09    = 10  & 
        & ,mgnin_lai10    = 11  & 
        & ,mgnin_lai11    = 12  & 
        & ,mgnin_lai12    = 13  & 
        & ,mgnin_pftp_bt  = 14  & 
        & ,mgnin_pftp_nt  = 15  & 
        & ,mgnin_pftp_sb  = 16  & 
        & ,mgnin_pftp_hb  = 17  & 
        & ,mgnin_tsa01    = 18  & 
        & ,mgnin_tsa02    = 19  & 
        & ,mgnin_tsa03    = 20  & 
        & ,mgnin_tsa04    = 21  & 
        & ,mgnin_tsa05    = 22  & 
        & ,mgnin_tsa06    = 23  & 
        & ,mgnin_tsa07    = 24  & 
        & ,mgnin_tsa08    = 25  & 
        & ,mgnin_tsa09    = 26  & 
        & ,mgnin_tsa10    = 27  & 
        & ,mgnin_tsa11    = 28  & 
        & ,mgnin_tsa12    = 29  & 
        & ,mgnin_swdown01 = 30  & 
        & ,mgnin_swdown02 = 31  & 
        & ,mgnin_swdown03 = 32  & 
        & ,mgnin_swdown04 = 33  & 
        & ,mgnin_swdown05 = 34  & 
        & ,mgnin_swdown06 = 35  & 
        & ,mgnin_swdown07 = 36  & 
        & ,mgnin_swdown08 = 37  & 
        & ,mgnin_swdown09 = 38  & 
        & ,mgnin_swdown10 = 39  & 
        & ,mgnin_swdown11 = 40  & 
        & ,mgnin_swdown12 = 41    

      CHARACTER*100 onefil



       
       CALL get_ijk_from_grid (  grid ,                           &
                                 ids, ide, jds, jde, kds, kde,    &
                                 ims, ime, jms, jme, kms, kme,    &
                                 ips, ipe, jps, jpe, kps, kpe    )

     WRITE( message , FMT='(A,4I5)' ) ' in input_ext_chem_megan2_file, DIMS: ',ids,ide-1,jds,jde-1
     CALL  wrf_message ( message )

     ALLOCATE( emiss(ids:ide-1,jds:jde-1,n_mgnin) )

     



     
     onefil='MEGAN_input_WRFchem.txt'


     

     OPEN(26,FILE=trim(onefil),FORM='FORMATTED', status='old')

     

     do i = ids, ide-1
        do j = jds, jde-1
           read (26, FMT='(2(I5,1x),41(ES11.2,1x))') itmp, jtmp, (emiss(i,j,v),v=1,n_mgnin)
           
           if ( (i /= itmp) .or. j /= jtmp ) then
              WRITE( message , FMT='(A,I3,I3,A,I3,I3)' ) 'Something is wrong (i,j) = ',i,j,"itmp, jtmp = ",itmp,jtmp
              call wrf_error_fatal3("<stdin>",349,&
message)
           end if
        end do
     end do


     
     grid%msebio_isop(ids:ide-1,jds:jde-1) = emiss(ids:ide-1,jds:jde-1,mgnin_isop)
     
     grid%mlai    (ids:ide-1,jds:jde-1,01) = emiss(ids:ide-1,jds:jde-1,mgnin_lai01)
     grid%mlai    (ids:ide-1,jds:jde-1,02) = emiss(ids:ide-1,jds:jde-1,mgnin_lai02)
     grid%mlai    (ids:ide-1,jds:jde-1,03) = emiss(ids:ide-1,jds:jde-1,mgnin_lai03)
     grid%mlai    (ids:ide-1,jds:jde-1,04) = emiss(ids:ide-1,jds:jde-1,mgnin_lai04)
     grid%mlai    (ids:ide-1,jds:jde-1,05) = emiss(ids:ide-1,jds:jde-1,mgnin_lai05)
     grid%mlai    (ids:ide-1,jds:jde-1,06) = emiss(ids:ide-1,jds:jde-1,mgnin_lai06)
     grid%mlai    (ids:ide-1,jds:jde-1,07) = emiss(ids:ide-1,jds:jde-1,mgnin_lai07)
     grid%mlai    (ids:ide-1,jds:jde-1,08) = emiss(ids:ide-1,jds:jde-1,mgnin_lai08)
     grid%mlai    (ids:ide-1,jds:jde-1,09) = emiss(ids:ide-1,jds:jde-1,mgnin_lai09)
     grid%mlai    (ids:ide-1,jds:jde-1,10) = emiss(ids:ide-1,jds:jde-1,mgnin_lai10)
     grid%mlai    (ids:ide-1,jds:jde-1,11) = emiss(ids:ide-1,jds:jde-1,mgnin_lai11)
     grid%mlai    (ids:ide-1,jds:jde-1,12) = emiss(ids:ide-1,jds:jde-1,mgnin_lai12)
     
     grid%pftp_bt  (ids:ide-1,jds:jde-1) = emiss(ids:ide-1,jds:jde-1,mgnin_pftp_bt)
     grid%pftp_nt  (ids:ide-1,jds:jde-1) = emiss(ids:ide-1,jds:jde-1,mgnin_pftp_nt)
     grid%pftp_sb  (ids:ide-1,jds:jde-1) = emiss(ids:ide-1,jds:jde-1,mgnin_pftp_sb)
     grid%pftp_hb  (ids:ide-1,jds:jde-1) = emiss(ids:ide-1,jds:jde-1,mgnin_pftp_hb)
     
     
     
     grid%mtsa    (ids:ide-1,jds:jde-1,01) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa01)
     grid%mtsa    (ids:ide-1,jds:jde-1,02) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa02)
     grid%mtsa    (ids:ide-1,jds:jde-1,03) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa03)
     grid%mtsa    (ids:ide-1,jds:jde-1,04) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa04)
     grid%mtsa    (ids:ide-1,jds:jde-1,05) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa05)
     grid%mtsa    (ids:ide-1,jds:jde-1,06) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa06)
     grid%mtsa    (ids:ide-1,jds:jde-1,07) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa07)
     grid%mtsa    (ids:ide-1,jds:jde-1,08) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa08)
     grid%mtsa    (ids:ide-1,jds:jde-1,09) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa09)
     grid%mtsa    (ids:ide-1,jds:jde-1,10) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa10)
     grid%mtsa    (ids:ide-1,jds:jde-1,11) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa11)
     grid%mtsa    (ids:ide-1,jds:jde-1,12) = emiss(ids:ide-1,jds:jde-1,mgnin_tsa12)
     
     
     grid%mswdown (ids:ide-1,jds:jde-1,01) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown01)
     grid%mswdown (ids:ide-1,jds:jde-1,02) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown02)
     grid%mswdown (ids:ide-1,jds:jde-1,03) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown03)
     grid%mswdown (ids:ide-1,jds:jde-1,04) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown04)
     grid%mswdown (ids:ide-1,jds:jde-1,05) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown05)
     grid%mswdown (ids:ide-1,jds:jde-1,06) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown06)
     grid%mswdown (ids:ide-1,jds:jde-1,07) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown07)
     grid%mswdown (ids:ide-1,jds:jde-1,08) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown08)
     grid%mswdown (ids:ide-1,jds:jde-1,09) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown09)
     grid%mswdown (ids:ide-1,jds:jde-1,10) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown10)
     grid%mswdown (ids:ide-1,jds:jde-1,11) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown11)
     grid%mswdown (ids:ide-1,jds:jde-1,12) = emiss(ids:ide-1,jds:jde-1,mgnin_swdown12)




    DEALLOCATE( emiss )

  end SUBROUTINE input_ext_chem_megan2_file




END MODULE module_input_chem_bioemiss

