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
!!----
!!--..
!!--..         REALWIN ZONE
!!--..
!!---- HISTORY
!!----
!!----    Update: 02/03/2011
!!----
!!---- DEPENDENCIES
!!--++    RealWin Library
!!----
!!---- VARIABLES
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ERROR_MESSAGE
!!----       INFO_MESSAGE
!!----       WRITE_SCROLL_TEXT
!!----
!!
 Module CFML_IO_Messages
    !---- Use Modules ----!
    use Realwin, only: message_box, Write_Scroll_Line,Scroll_Text, Create_Window,Select_Font

    !---- Definitions ----!
    implicit none

    !---- List private variables ----!
    integer, private :: messval

    !---- List of public subroutines ----!
    public :: error_message, info_message, write_scroll_text

    !---- Definitions ----!
    !!--++
    !!--++ ICWINDOW
    !!--++    integer, private :: icwindow
    !!--++
    !!--++    Code number for Scroll Window
    !!--++
    !!--++ Update: March - 2005
    !!
    integer, private :: icwindow= -1

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
    !!---- Subroutine Error_Message(Mess, Iunit, Routine, Faltal)
    !!----    character(len=*), intent(in)           :: Mess    !  In -> Error information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an error message on the screen and in 'Iunit' if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Error_Message(Mess, Iunit, Routine, Fatal)
       !---- Arguments ----!
       character(len=*),            intent(in) :: Mess
       integer, optional,           intent(in) :: iunit
       Character(Len =*), Optional, Intent(In) :: Routine
       Logical, Optional,           Intent(In) :: Fatal

       messval=message_box(text="Error: "//trim(mess),title="Stop/Warning Box")  !rw

       if (present(iunit)) then
          write(unit=iunit,fmt="(1x,a)") "***"
          write(unit=iunit,fmt="(1x,a)") "*** ERROR: "//mess
          write(unit=iunit,fmt="(1x,a)") "***"
          write(unit=iunit,fmt="(1x,a)") " "
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
    !!---- Subroutine Info_Message(Mess, Iunit, Scroll_Window)
    !!----    character(len=*),  intent(in) :: Mess           !  In -> Info information
    !!----    integer, optional, intent(in) :: Iunit          !  In -> Write information on Iunit unit
    !!----    integer, optional, intent(in) :: Scroll_Window  !  In -> Write information on scroll windows
    !!----
    !!----    Print an message on the screen and in 'Iunit' if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Info_Message(Mess, Iunit, Scroll_Window)
       character(len=*), intent(in)           :: Mess
       integer,          intent(in), optional :: iunit
       integer,          intent(in), optional :: scroll_window

       if(present(scroll_window)) then
         call write_scroll_line(scroll_window,text=trim(Mess))
       else
         messval=message_box(text=trim(Mess),title="Info Message")  !rw
       end if
       if (present(iunit) ) then
          write(unit=iunit,fmt="(1x,a)") " "
          write(unit=iunit,fmt="(1x,a)") " "//trim(Mess)
          write(unit=iunit,fmt="(1x,a)") " "
       end if
       return
    End Subroutine Info_Message

    !!----
    !!---- Subroutine Write_Scroll_Text(Mess)
    !!----    character(len=*), intent(in)           :: Mess
    !!----
    !!----    Print the string in a the scroll window
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Write_Scroll_Text(Mess)
       !---- Argument ----!
       character(len=*), intent(in) :: Mess


       !---- Open the Scroll Window if necessary ----!
       if (.not. wscroll) then
          icwindow = create_window(window_name="Scroll Text Window",x=0.15,y=0.07, &
                            width= 0.7, height=0.42, &
                            paint_code=SCROLL_TEXT,&
                            text_font=select_font(typeface='courier',point=8))
         wscroll=.true.
       end if
       call write_scroll_line(icwindow,text=trim(mess))

       return
    End Subroutine Write_Scroll_Text

 End Module CFML_IO_Messages
