





































MODULE module_input_gocart_dms

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


SUBROUTINE input_ext_chem_gocart_dms (grid)


   IMPLICIT NONE

   TYPE(domain)           ::  grid

   INTEGER ::  i,j,n,numfil,status,system

   INTEGER :: ids, ide, jds, jde, kds, kde,    &
              ims, ime, jms, jme, kms, kme,    &
              ips, ipe, jps, jpe, kps, kpe

   REAL, ALLOCATABLE, DIMENSION(:,:) :: tmp



      PARAMETER(numfil=19)

   CHARACTER (LEN=80) :: message

   TYPE (grid_config_rec_type)              :: config_flags


      CHARACTER*100 onefil



       
       CALL get_ijk_from_grid (  grid ,                        &
                                 ids, ide, jds, jde, kds, kde,    &
                                 ims, ime, jms, jme, kms, kme,    &
                                 ips, ipe, jps, jpe, kps, kpe    )

     WRITE( message , FMT='(A,4I5)' ) ' DIMS: ',ids,ide-1,jds,jde-1
     CALL  wrf_message ( message )

     ALLOCATE(   tmp(ids:ide-1,jds:jde-1) )



     OPEN(19,FILE='wrf_gocart_dms.txt',FORM='FORMATTED')
     READ(19,'(12E15.5)') tmp
     grid%dms_0(ids:ide-1,jds:jde-1) = tmp
     CLOSE(19)



      DEALLOCATE( tmp )


END SUBROUTINE input_ext_chem_gocart_dms 


END MODULE module_input_gocart_dms

