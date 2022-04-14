module module_mosaic_support
  
  
  

  use module_data_mosaic_asecthp, only: lunerr
  use module_peg_util,            only: peg_error_fatal, peg_message

  implicit none
  private

  public:: mosaic_warn_mess
  public:: mosaic_err_mess
  
  
contains

  subroutine mosaic_warn_mess(message)
    

    character(len=*), intent(in) :: message

    
    character(len=16), parameter :: warn_str = 'MOSAIC WARNING: ' 
    character(len=500) :: str_to_prnt 



    str_to_prnt = warn_str // message
    call peg_message( lunerr, trim(str_to_prnt) )


    
  end subroutine mosaic_warn_mess


  subroutine mosaic_err_mess(message)
    
    character(len=*), intent(in) :: message

    
    character(len=14), parameter :: err_str = 'MOSAIC ERROR: ' 
    character(len=500) :: str_to_prnt 

    write(str_to_prnt,*)err_str,message
    



    str_to_prnt = err_str // message
    call peg_error_fatal( lunerr, trim(str_to_prnt) )
  end subroutine mosaic_err_mess



end module module_mosaic_support
