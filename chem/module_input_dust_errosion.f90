





































MODULE module_input_dust_errosion

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


SUBROUTINE input_ext_chem_erod_file (grid)


   IMPLICIT NONE

   TYPE(domain)           ::  grid

   INTEGER ::  i,j,n,numfil,status,system

   INTEGER :: ids, ide, jds, jde, kds, kde,    &
              ims, ime, jms, jme, kms, kme,    &
              ips, ipe, jps, jpe, kps, kpe

   REAL, ALLOCATABLE, DIMENSION(:,:) :: tmp
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: erod


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

     ALLOCATE( tmp(ids:ide-1,jds:jde-1) )
     ALLOCATE( erod(ids:ide-1,jds:jde-1,3) )


     OPEN(19,FILE='wrf_dust_erod.txt',FORM='FORMATTED')
     READ(19,'(12E15.5)') tmp
     grid%erod(ids:ide-1,jds:jde-1,1) = tmp
     READ(19,'(12E15.5)') tmp
     grid%erod(ids:ide-1,jds:jde-1,2) = tmp
     READ(19,'(12E15.5)') tmp
     grid%erod(ids:ide-1,jds:jde-1,3) = tmp
     CLOSE(19)



      DEALLOCATE( erod )
      DEALLOCATE( tmp )

END SUBROUTINE input_ext_chem_erod_file 


END MODULE module_input_dust_errosion

