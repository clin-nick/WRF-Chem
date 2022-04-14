







	module module_peg_util


	contains



	subroutine peg_debugmsg( lun, level, str )




	implicit none

	integer, intent(in) :: lun, level
	character(len=*), intent(in) :: str

	integer n

	n = max( 1, len_trim(str) )
	if (lun .gt. 0) then
	    write(lun,'(a)') str(1:n)
	else
	    call wrf_debug( level, str(1:n) )
	end if
	return
	end subroutine peg_debugmsg



	subroutine peg_message( lun, str )




	implicit none

	integer, intent(in) :: lun
	character(len=*), intent(in) :: str

	integer n

	n = max( 1, len_trim(str) )
	if (lun .gt. 0) then
	    write(lun,'(a)') str(1:n)
	else
	    call wrf_message( str(1:n) )
	end if
	return
	end subroutine peg_message



	subroutine peg_error_fatal( lun, str )




	implicit none

	integer, intent(in) :: lun
	character(len=*), intent(in) :: str

	integer n

	n = max( 1, len_trim(str) )
	if (lun .gt. 0) write(lun,'(a)') str(1:n)
	call wrf_error_fatal3("<stdin>",76,&
str(1:n) )
	return
	end subroutine peg_error_fatal



	end module module_peg_util
