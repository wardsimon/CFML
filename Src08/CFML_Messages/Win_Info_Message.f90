!!----
SubModule (CFML_Messages) Wmess_Infos
  Implicit none
   Contains

   !!----
   !!---- INFO_MESSAGE
   !!----
   !!----    Print an message on the screen or in "Iunit" if present
   !!----
   !!---- 03/05/2019
   !!
   Module Subroutine Info_Message(Mess, iunit)
      !---- Arguments ----!
      character(len=*), intent(in)           :: Mess        ! Info information
      integer,          intent(in), optional :: iunit       ! Write information on Iunit unit

      call WMessageBox(OKOnly, InformationIcon, CommonOK, trim(Mess),"CrysFML: Information Message")

      if (present(iunit) ) then
         write(unit=iunit,fmt="(tr1,a)") " "
         write(unit=iunit,fmt="(tr1,a)") " "//trim(Mess)
         write(unit=iunit,fmt="(tr1,a)") " "
      end if

      return
   End Subroutine Info_Message

End SubModule Wmess_Infos