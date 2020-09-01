!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----          Nebil Ayape Katcho      (ILL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
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
!!----    Update: 02/03/2011
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

    use CFML_GlobalDeps, only: OPS

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
    !!---- Subroutine Error_Message(Mess, Iunit, Routine, Fatal)
    !!----    character(len=*), intent(in)           :: Mess    !  In -> Error information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----    character(len=*), intent(in), optional :: routine !  Added for consistency with the CFML_IO_Mess.f90 version.
    !!----    logical,          intent(in), optional :: fatal   !  Added for consistency with the CFML_IO_Mess.f90 version.
    !!----
    !!----    Print an error message on the screen and in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Error_Message(Mess, Iunit, Routine, Fatal)
       !---- Arguments ----!
       character(len=*),            intent(in) :: Mess
       integer, optional,           intent(in) :: iunit
       Character(Len =*), Optional, Intent(In) :: Routine
       Logical, Optional,           Intent(In) :: Fatal

       call WMessageBox(OKOnly, ExclamationIcon, CommonOK, Mess,"Error Message")

       if (present(iunit)) then
          write(unit=iunit,fmt="(tr1,a)") "****"
          write(unit=iunit,fmt="(tr1,a)") "**** ERROR: "//trim(Mess)
          write(unit=iunit,fmt="(tr1,a)") "****"
          write(unit=iunit,fmt="(tr1,a)") " "
          If (Present(Routine)) Then
            Write(Unit = iunit, Fmt = "(tr1,a)") "**** Subroutine: "//trim(Routine)
          End If
          If (Present(Fatal)) Then
           If (Fatal) Then
               Write(Unit = iunit, Fmt = "(/tr1,a)") "**** The Program Will Stop Here."
               Stop
           End If
          End If
       end if

       return
    End Subroutine Error_Message

    !!----
    !!---- Subroutine Info_Message(Mess, Iunit)
    !!----    character(len=*), intent(in)           :: Mess    !  In -> Info information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an message on the screen or in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Info_Message(Mess, iunit)
       !---- Arguments ----!
       character(len=*), intent(in)           :: Mess
       integer,          intent(in), optional :: iunit

       call WMessageBox(OKOnly, InformationIcon, CommonOK, Mess,"Information Message")

       if (present(iunit) ) then
          write(unit=iunit,fmt="(tr1,a)") " "
          write(unit=iunit,fmt="(tr1,a)") " "//trim(Mess)
          write(unit=iunit,fmt="(tr1,a)") " "
       end if

       return
    End Subroutine Info_Message

    !!----
    !!---- Subroutine Question_Message(Mess,Title
    !!----    character(len=*)  :: Mess
    !!----
    !!----    Print an question on the screen
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Question_Message(Mess,Title)
       !---- Argument ----!
       character (len=*),           intent(in) :: Mess
       character (len=*), optional, intent(in) :: Title

       !---- Variable ----!
       character(len=80) :: ch_title

       ch_title=' '
       if (present(title)) ch_title=title
       call WMessageBox(YesNo,QuestionIcon,CommonYes,Mess,ch_title)

       return
    End Subroutine Question_Message

    !!----
    !!---- Subroutine Stop_Message(Mess,Title)
    !!----    character(len=*)  :: Mess
    !!----
    !!----    Print an Stop on the screen
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Stop_Message(Mess,Title)
       !---- Argument ----!
       character (len=*),           intent(in) :: Mess
       character (len=*), optional, intent(in) :: Title

       !---- Variable ----!
       character(len=80) :: ch_title

       ch_title=' '
       if (present(title)) ch_title=title
       call WMessageBox(YesNo,StopIcon,CommonYes,Mess,ch_title)

       return
    End Subroutine Stop_Message

    !!----
    !!---- SUBROUTINE WARNING_MESSAGE(Mess, Iunit)
    !!----    character(len=*), intent(in)           :: Mess    !  In -> Info information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an message on the screen or in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Warning_Message(Mess, Iunit)
       !---- Arguments ----!
       character(len=*), intent(in) :: Mess
       integer, optional,intent(in) :: iunit

       call WMessageBox(OKOnly,ExclamationIcon,CommonOK, Mess,"Warning Message")

       if (present(iunit) ) then
          write(unit=iunit,fmt="(tr1,a)") "****"
          write(unit=iunit,fmt="(tr1,a)") "**** WARNING: "//trim(Mess)
          write(unit=iunit,fmt="(tr1,a)") "****"
       end if

       return
    End Subroutine Warning_Message

    !!----
    !!---- Subroutine Write_Scroll_Text(Mess, ICmd)
    !!----    character(len=*), intent(in)  :: Mess    ! Message to write
    !!----    integer, optional,intent(in)  :: ICmd    ! Define the type of the Editor opened
    !!----
    !!----    Print the string in a the window
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Write_Scroll_Text(Mess,ICmd)
       !---- Argument ----!
       character(len=*), intent(in) :: Mess
       integer, optional,intent(in) :: ICmd

       !---- Local variables ----!
       character(len=2), save :: newline = char(13)//char(10)
       integer                :: iendbuf, ic

       ic=0
       if (present(ICmd)) ic=ICmd
       !---- Open the Scroll Window if necessary ----!
       if (.not. wscroll) then
          if(OPS /= 1) newline=char(10)
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
       call WEditPutTextPart(trim(Mess)//newline,iendbuf)
       call WindowSelect(0)

       return
    End Subroutine Write_Scroll_Text

 End Module CFML_IO_Messages
