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
!!----    Update: 02/03/2011
!!----
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ERROR_MESSAGE
!!----       INFO_MESSAGE
!!----       PRINT_MESSAGE
!!----       WAIT_MESSAGE
!!----       WRITE_SCROLL_TEXT
!!----
!!
 Module CFML_IO_Messages
    !---- Use Modules ----!

    !---- Definitions ----!
    implicit none

    !---- List of public subroutines ----!
    public :: Info_Message, Error_Message, Print_Message, Wait_Message, Write_Scroll_Text


 Contains

    !!----
    !!---- Subroutine Error_Message(Mess, Iunit, Routine, Fatal)
    !!----    character(len=*), intent(in)           :: Mess          !  In -> Error information
    !!----    integer,          intent(in), optional :: Iunit         !  In -> Write information on Iunit unit
    !!----    character(len=*), intent(in), optional :: Routine       !  In -> The subroutine where the error occured
    !!----    logical,          intent(in), optional :: Fatal         !  In -> Should the program stop here ?
    !!----
    !!----    Print an error message on the screen or in "Iunit" if present
    !!----    If "routine" is present the subroutine where the occured will be also displayed.
    !!----    If "fatal" is present and .True. the program will stop after the printing.
    !!----
    !!---- Update: January - 2010
    !!
    Subroutine Error_Message(Mess, Iunit, Routine, Fatal)
       !---- Arguments ----!
       Character ( Len = * ), Intent(In)           :: Mess
       Integer,               Intent(In), Optional :: Iunit
       Character ( Len = * ), Intent(In), Optional :: Routine
       Logical,               Intent(In), Optional :: Fatal

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
    !!---- Subroutine Info_Message(Mess, Iunit)
    !!----    character(len=*), intent(in)           :: Mess          !  In -> Info information
    !!----    integer,          intent(in), optional :: Iunit         !  In -> Write information on Iunit unit
    !!----
    !!----    Print an message on the screen or in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Info_Message(Mess, iunit)
       !---- Arguments ----!
       character(len=*), intent(in)           :: Mess
       integer,          intent(in), optional :: iunit

       !---- Local Variables ----!
       integer :: lun

       lun=6
       if (present(iunit)) lun=iunit
       write(unit=lun,fmt="(a)") " => "//Mess

       return
    End Subroutine Info_Message

    !!----
    !!---- Subroutine Print_Message(Mess)
    !!----    character(len=*), intent(in)  :: Mess    !  In -> Print information
    !!----
    !!----    Print an message on the screen
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Print_Message(Mess)
       !---- Arguments ----!
       character(len=*),intent(in) ::  Mess

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
    !!---- Subroutine Wait_Message(Mess)
    !!----    character(len=*), optional, intent(in) :: Mess
    !!----
    !!----    Similar to Pause for Console version
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Wait_Message(Mess)
       !---- Argument ----!
       character(len=*), optional, intent(in) :: Mess

       !---- Local variable ----!
       character(len=1) :: car

       write(unit=*,fmt="(a)") " "
       if (present(mess)) write(unit=*,fmt="(a)", advance="no") mess
       read(unit=*,fmt="(a)") car
       if( len_trim(car) == 0) return

       return
    End Subroutine Wait_Message

    !!----
    !!---- Subroutine Write_Scroll_Text(Mess)
    !!----    character(len=*), intent(in) :: Mess
    !!----
    !!----    Print the string in a default output unit
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Write_Scroll_Text(Mess)
       !---- Argument ----!
       character(len=*), intent(in) :: Mess

       write(unit=*, fmt="(a)") trim(mess)

       return
    End Subroutine Write_Scroll_Text

 End Module CFML_IO_Messages
