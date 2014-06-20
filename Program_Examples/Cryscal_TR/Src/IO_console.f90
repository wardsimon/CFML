!     Last change:  TR   11 Dec 2007    5:36 pm
module IO_module

 implicit none

 PRIVATE
 PUBLIC    :: read_input_line, write_info

 contains


 subroutine read_input_line(string_info)
  character (len=256) , intent(out) :: string_info

  do
   read(*, '(a)')   string_info
   if(len_trim(string_info)/=0) exit
  end do

  return
 end subroutine read_input_line

 !-----------------------------------------------------------
 subroutine write_info(string_info)
   character (len=*) , intent(in) :: string_info
   integer                        :: i_error

  if(len_trim(string_info) == 0) then
   write(*,'(2x,a)') ' '
   write(2,'(2x,a)', iostat=i_error) ' '
   if(i_error /=0) then
    write(*,*) ' Unable to write in cryscal.log file. Please check permissions in the current folder !'
	write(*,*) ' Program will be stopped !'
	stop
   endif
  else
   write(*,'(2x,a)') string_info
   write(2,'(2x,a)', iostat=i_error) string_info
   if(i_error /=0) then
    write(*,*) ' Unable to write in cryscal.log file. Please check permissions in the current folder !'
	write(*,*) ' Program will be stopped !'
	stop
   endif
  endif
  return

 end subroutine write_info

end module IO_module
