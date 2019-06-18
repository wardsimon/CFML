!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2019  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!----               Ross Angel         (University of Pavia) 
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
!!----
!!
 Module CFML_Messages
    !---- Use Modules ----!
    use Winteracter, only: YesNo,OKOnly,CommonYes, CommonOK, Modeless, ViewOnly,         &
                           WordWrap, NoMenu, NoToolbar, SystemFixed, EditTextLength,     &
                           ScreenHeight, StopIcon,InformationIcon, ExclamationIcon,      &
                           QuestionIcon, WMessageBox, WindowCloseChild, WindowOpenChild, &
                           WEditFile, WEditPutTextPart, WindowSelect, WInfoEditor,       &
                           CommandLine,WInfoScreen,CourierNew,win_message

    use CFML_GlobalDeps, only: OPS, NEWLINE

    !---- Definitions ----!
    implicit none

    private

    !---- List of public subroutines ----!
    public :: Close_scroll_window, error_message, info_message, question_message, warning_message, &
              stop_message, write_scroll_text

    !---- Definitions ----!
    integer, public :: win_console= -1       ! Number of the Scroll window
    logical         :: wscroll    = .false.  ! Logical variable to indicate if the Scroll Window is active

    !---- Interfaces ----!
    Interface
       Module Subroutine Close_Scroll_Window()
          !---- Arguments ----!
       End Subroutine Close_Scroll_Window
       
       Module Subroutine Error_Message(Mess, Iunit, Routine, Fatal)
          !---- Arguments ----!
          character(len=*),            intent(in) :: Mess         
          integer, optional,           intent(in) :: iunit        
          Character(Len =*), Optional, Intent(In) :: Routine      
          Logical,           Optional, Intent(In) :: Fatal        
       End Subroutine Error_Message
       
       Module Subroutine Info_Message(Mess, iunit)
          !---- Arguments ----!
          character(len=*), intent(in)           :: Mess        
          integer,          intent(in), optional :: iunit       
       End Subroutine Info_Message
       
       Module Subroutine Question_Message(Mess, Title)
          !---- Argument ----!
          character (len=*),           intent(in) :: Mess        
          character (len=*), optional, intent(in) :: Title    
       End Subroutine Question_Message
       
       Module Subroutine Stop_Message(Mess, Title)
          !---- Argument ----!
          character (len=*),           intent(in) :: Mess      
          character (len=*), optional, intent(in) :: Title    
       End Subroutine Stop_Message 
       
       Module Subroutine Warning_Message(Mess, Iunit)
          !---- Arguments ----!
          character(len=*), intent(in) :: Mess             
          integer, optional,intent(in) :: iunit           
       End Subroutine Warning_Message
       
       Module Subroutine Write_Scroll_Text(Mess,ICmd)
          !---- Argument ----!
          character(len=*), intent(in) :: Mess      
          integer, optional,intent(in) :: ICmd      
       End Subroutine Write_Scroll_Text
    
    End Interface

 End Module CFML_Messages
