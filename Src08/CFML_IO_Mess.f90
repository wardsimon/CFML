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
!!--..
!!--..         NON-GRAPHICS ZONE
!!--..
!!---- HISTORY
!!----
!!----    Update: 11/07/2015
!!----
!!----
!!
 Module CFML_IO_Messages
    !---- Use Modules ----!

    !---- Definitions ----!
    implicit none

    public :: error_message, info_message, print_message, wait_message, write_scroll_text

 Contains

    !!----
    !!---- SUBROUTINE ERROR_MESSAGE
    !!----
    !!----    Print an error message on the screen or in "Iunit" if present
    !!----    If "routine" is present the subroutine where the occured will be also displayed.
    !!----    If "fatal" is present and .True. the program will stop after the printing.
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Error_Message(Mess, Iunit, Routine, Fatal)
       !---- Arguments ----!
       Character(Len=*),           Intent(In) :: Mess       ! Error information
       Integer,          optional, Intent(In) :: Iunit      ! Write information on Iunit unit
       Character(Len=*), optional, Intent(In) :: Routine    ! The subroutine where the error occured
       Logical,          optional, Intent(In) :: Fatal      ! Should the program stop here ?

       !---- Local Variables ----!
       Integer          :: Lun, Lenm, Lenr
       character(len=1) :: ent

       Lun = 6
       If (Present(Iunit)) Lun = Iunit

       If (Present(Routine)) Then
           Lenr = Len_Trim(Routine)
           Write(Unit = Lun, Fmt = "(A)") " => Error in Subroutine: "//Routine(1:Lenr)
       End If

       Lenm = Len_Trim(Mess)
       Write(Unit = Lun, Fmt = "(A)") " => Message Error: "//Mess(1:Lenm)

       If (Present(Fatal)) Then
           If (Fatal) Then
               Write(Unit = Lun, Fmt = "(A)") " => Fatal error: the program stops here."
               Write(Unit= *, Fmt="(/,a)") " => Press <enter> to finish "
               Read (Unit= *, Fmt="(a)") ent
               Stop
           End If
       End If

       Return
    End Subroutine Error_Message

    !!----
    !!---- SUBROUTINE INFO_MESSAGE
    !!----
    !!----    Print an message on the screen or in "Iunit" if present
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Info_Message(Mess, Iunit)
       !---- Arguments ----!
       character(len=*), intent(in)           :: Mess       ! Info information
       integer,          intent(in), optional :: iunit      ! Write information on Iunit unit

       !---- Local Variables ----!
       integer :: lun

       lun=6
       if (present(iunit)) lun=iunit
       write(unit=lun,fmt="(a)") " => "//trim(Mess)

       return
    End Subroutine Info_Message

    !!----
    !!---- SUBROUTINE PRINT_MESSAGE
    !!----
    !!----    Print an message on the screen
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Print_Message(Mess)
       !---- Arguments ----!
       character(len=*),intent(in) ::  Mess ! Print information

       !---- Local Variables ----!
       integer :: lon

       lon=len_trim(mess)
       if (lon == 0) then
          write(unit=*,fmt="(a)") "  "
       else
          if (mess(1:1) == "=" .or. mess(2:2) == "=") then
             write(unit=*,fmt="(a)") mess(1:lon)
          else
             write(unit=*,fmt="(a,a)")" =>", mess(1:lon)
          end if
       end if

       return
    End Subroutine Print_Message

    !!----
    !!---- SUBROUTINE WAIT_MESSAGE
    !!----
    !!----    Similar to Pause for Console version
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Wait_Message(Mess)
       !---- Argument ----!
       character(len=*), optional, intent(in) :: Mess    ! Message

       !---- Local variable ----!
       character(len=1) :: car

       write(unit=*,fmt="(a)") " "
       if (present(mess)) write(unit=*,fmt="(a)", advance="no") mess
       read(unit=*,fmt="(a)") car
       if( len_trim(car) == 0) return

       return
    End Subroutine Wait_Message

    !!----
    !!---- SUBROUTINE WRITE_SCROLL_TEXT
    !!----
    !!----    Print the string in a default output unit
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Write_Scroll_Text(Mess)
       !---- Argument ----!
       character(len=*), intent(in) :: Mess   ! Message

       write(unit=*, fmt="(a)") trim(mess)

       return
    End Subroutine Write_Scroll_Text

 End Module CFML_IO_Messages
