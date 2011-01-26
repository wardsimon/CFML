!     Last change:  TR   12 Dec 2007    9:38 am
module IO_module

 use Realwin, only: message_box, write_scroll_line, SCROLL_TEXT, create_window,select_font, &
                    read_scroll_line

 implicit none
 PRIVATE
 PUBLIC    :: read_input_line, write_info

 integer, private             :: cryscal_window = -1
 logical, private             :: wscroll = .false.
 character (len=256), private :: scroll_line
 integer, private             :: windows_status

 contains


 subroutine read_input_line(string_info)
  character (len=256) , intent(out) :: string_info

  !do
  ! read(*, '(a)')   string_info
  ! if(len_trim(string_info)/=0) exit
  !end do

  do
   call read_scroll_line(window= cryscal_window, text=scroll_line, stat = windows_status)
   IF(windows_status /=0) stop
   if(len_trim(scroll_line) /=0) exit
  end do
  string_info = scroll_line

  return
 end subroutine read_input_line

 !-----------------------------------------------------------
 subroutine write_info(string_info)
  !use cryscal_RW_module, only : main_window
  ! ecrit un message dans la fenetre main_window et dans le fichier CRYSCAL.OUT
   character (len=*) , intent(in) :: string_info

!  if(len_trim(string_info) == 0) then
!   call write_scroll_line(window = main_window, text = ' ')
!  else
!   call write_scroll_line(window = main_window, text = string_info)
!  endif
!  write(2, '(a)') string_info
!
!  return


  !---- Creation de la fenetre cryscal_window si necessaire ----!
  if (.not. wscroll) then
     cryscal_window = create_window(window_name="CRYSCAL (TR/CDIFX/2009)",x=0.1,y=0.01,    &
                                    width= 0.8, height=0.45,                               &
                                    paint_code=SCROLL_TEXT,                                &
                                    text_font=select_font(typeface='courier',point=8))
      wscroll=.true.
  end if
  call write_scroll_line(cryscal_window,text=string_info)

  write(2, '(2x,a)') string_info
 end subroutine write_info

end module IO_module
