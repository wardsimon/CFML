!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: CFML_IO_MESSAGES
!!----   INFO: Input / Output General Messages. It is convenient to use these intermediate procedures instead of
!!----         Fortran Write(*,*) or Print*, because it is much more simple to modify a program for making a GUI.
!!----         Usually GUI tools and libraries need special calls to dialog boxes for screen messages. These
!!----         calls may be implemented within this module using the same name procedures. The subroutines
!!----         ERROR_MESSAGE and INFO_MESSAGE are just wrappers for the actual calls.
!!--..
!!--.. WINTERACTER ZONE
!!--..
!!---- HISTORY
!!----
!!----    Update: February - 2005
!!----            June    - 1999   Updated by JGP
!!----
!!---- DEPENDENCIES
!!--++    Winteracter or X/Winteracter Library
!!----
!!---- VARIABLES
!!----    WIN_CONSOLE
!!--++    WSCROLL                      [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       CLOSE_SCROLL_WINDOW
!!----       ERROR_MESSAGE
!!----       INFO_MESSAGE
!!----       QUESTION_MESSAGE
!!----       STOP_MESSAGE
!!----       WARNING_MESSAGE
!!----       WRITE_SCROLL_TEXT
!!----
!!
 Module CFML_IO_Messages
    !---- Use Modules ----!
    use Winteracter, only: YesNo,OKOnly,CommonYes, CommonOK, Modeless, ViewOnly,         &
                           WordWrap, NoMenu, NoToolbar, SystemFixed, EditTextLength,     &
                           ScreenHeight, StopIcon,InformationIcon, ExclamationIcon,      &
                           QuestionIcon, WMessageBox, WindowCloseChild, WindowOpenChild, &
                           WEditFile, WEditPutTextPart, WindowSelect, WInfoEditor,       &
                           CommandLine,WInfoScreen,CourierNew,win_message

    !---- Definitions ----!
    implicit none

    private

    !---- List of public subroutines ----!
    public :: Close_scroll_window, error_message, info_message, question_message, warning_message, &
              stop_message, write_scroll_text

    !---- Definitions ----!
    !!----
    !!---- WIN_CONSOLE
    !!----    integer, private :: icwindow
    !!----
    !!----    Code number for Scroll Window
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: win_console= -1

    !!--++
    !!--++ WSCROLL
    !!--++    logical, private :: wscroll
    !!--++
    !!--++    Logical variable to indicate if the Scroll Window is
    !!--++    active or not.
    !!--++
    !!--++ Update: March - 2005
    !!
    logical, private :: wscroll = .false.

 Contains
    !!----
    !!---- SUBROUTINE CLOSE_SCROLL_WINDOW()
    !!----
    !!----    Close the Scroll Window
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Close_Scroll_Window()

       if (wscroll) call WindowCloseChild(win_console)
       win_console= -1
       wscroll=.false.

       return
    End Subroutine Close_Scroll_Window

    !!----
    !!---- Subroutine Error_Message(Line, Iunit)
    !!----    character(len=*), intent(in)           :: Line    !  In -> Error information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an error message on the screen and in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Error_Message(line,iunit)
       !---- Arguments ----!
       character(len=*), intent(in) :: line
       integer, optional,intent(in) :: iunit

       call WMessageBox(OKOnly, ExclamationIcon, CommonOK, line,"Error Message")

       if (present(iunit)) then
          write(unit=iunit,fmt="(1x,a)") "****"
          write(unit=iunit,fmt="(1x,a)") "**** ERROR: "//line
          write(unit=iunit,fmt="(1x,a)") "****"
          write(unit=iunit,fmt="(1x,a)") " "
       end if

       return
    End Subroutine Error_Message

    !!----
    !!---- Subroutine Info_Message(Line, Iunit)
    !!----    character(len=*), intent(in)           :: Line    !  In -> Info information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an message on the screen or in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Info_Message(line, iunit)
       !---- Arguments ----!
       character(len=*), intent(in)           :: line
       integer,          intent(in), optional :: iunit

       call WMessageBox(OKOnly, InformationIcon, CommonOK, line,"Information Message")

       if (present(iunit) ) then
          write(unit=iunit,fmt="(1x,a)") " "
          write(unit=iunit,fmt="(1x,a)") " "//line
          write(unit=iunit,fmt="(1x,a)") " "
       end if

       return
    End Subroutine Info_Message

    !!----
    !!---- SUBROUTINE Question_Message(line,titl)
    !!----    character(len=*)  :: line
    !!----
    !!----    Print an question on the screen
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Question_Message(line,titl)
       !---- Argument ----!
       character (len=*),           intent(in) :: line
       character (len=*), optional, intent(in) :: Titl

       !---- Variable ----!
       character(len=80) :: ch_title

       ch_title=' '
       if (present(titl)) ch_title=titl
       call WMessageBox(YesNo,QuestionIcon,CommonYes,line,ch_title)

       return
    End Subroutine Question_Message

    !!----
    !!---- SUBROUTINE Stop_Message(line,Titl)
    !!----    character(len=*)  :: line
    !!----
    !!----    Print an Stop on the screen
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Stop_Message(line,Titl)
       !---- Argument ----!
       character (len=*),           intent(in) :: line
       character (len=*), optional, intent(in) :: Titl

       !---- Variable ----!
       character(len=80) :: ch_title

       ch_title=' '
       if (present(titl)) ch_title=titl
       call WMessageBox(YesNo,StopIcon,CommonYes,line,ch_title)

       return
    End Subroutine Stop_Message

    !!----
    !!---- SUBROUTINE WARNING_MESSAGE(Line, Iunit)
    !!----    character(len=*), intent(in)           :: Line    !  In -> Info information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an message on the screen or in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Warning_Message(line, iunit)
       !---- Arguments ----!
       character(len=*), intent(in) :: line
       integer, optional,intent(in) :: iunit

       call WMessageBox(OKOnly,ExclamationIcon,CommonOK, line,"Warning Message")

       if (present(iunit) ) then
          write(unit=iunit,fmt="(1x,a)") "****"
          write(unit=iunit,fmt="(1x,a)") "**** WARNING: "//line
          write(unit=iunit,fmt="(1x,a)") "****"
       end if

       return
    End Subroutine Warning_Message

    !!----
    !!---- Subroutine Write_Scroll_Text(Line)
    !!----    character(len=*), intent(in)           :: Line
    !!----
    !!----    Print the string in a the scroll window
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Write_Scroll_Text(Line,icode)
       !---- Argument ----!
       character(len=*), intent(in) :: line
       integer, optional,intent(in) :: icode

       !---- Local variables ----!
       character(len=2), parameter :: newline = char(13)//char(10)
       integer                     :: iendbuf, ic

       ic=0
       if (present(icode)) ic=icode
       !---- Open the Scroll Window if necessary ----!
       if (.not. wscroll) then
          call WindowOpenChild(win_console,height=nint(WInfoScreen(ScreenHeight)*0.5), &
                               title='Info Window')
          select case (ic)
             case (0)
                call WEditFile(" ",IFlags=ViewOnly+NoMenu+NoToolbar+CommandLine,             &
                             IFont=CourierNew)
             case (1)
                call WEditFile(" ",IFlags=ViewOnly+NoMenu+NoToolbar,                         &
                             IFont=CourierNew)
          end select
          wscroll=.true.
       end if

       call WindowSelect(win_console)
       iendbuf=WInfoEditor(win_console,EditTextLength)+1
       call WEditPutTextPart(trim(line)//newline,iendbuf)
       call WindowSelect(0)

       return
    End Subroutine Write_Scroll_Text

 End Module CFML_IO_Messages
