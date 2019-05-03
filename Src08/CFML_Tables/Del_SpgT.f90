!!----
!!----
!!----
SubModule (CFML_Symmetry_Tables) Del_Routines
  
   Contains
   !!----
   !!---- Remove_Spgr_Info()
   !!----
   !!----    Deallocating SPGR_INFO Data
   !!----
   !!---- 23/04/2019 
   !!
   Module Subroutine Remove_Spgr_Info()

      if (allocated(spgr_info)) deallocate(spgr_info)

      return
   End Subroutine Remove_Spgr_Info

   !!----
   !!---- Remove_System_Equiv()
   !!----
   !!----    Deallocating SPGR_INFO Data
   !!----
   !!---- 23/04/2019 
   !!
   Module Subroutine Remove_System_Equiv()

      if (allocated(System_Equiv)) deallocate(System_Equiv)

      return
   End Subroutine Remove_System_Equiv

   !!----
   !!---- Remove_Wyckoff_Info()
   !!----
   !!----    Deallocating WYCKOFF_INFO Data
   !!----
   !!---- 23/04/2019 
   !!
   Module Subroutine Remove_Wyckoff_Info()

      if (allocated(wyckoff_info)) deallocate(wyckoff_info)

      return
   End Subroutine Remove_Wyckoff_Info 
   
End SubModule Del_Routines 