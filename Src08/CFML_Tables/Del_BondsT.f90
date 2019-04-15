!!----
!!----
!!----
SubModule (CFML_Bonds_Tables) Del_Routines
  
  Contains
   
   !!----
   !!---- REMOVE_BONDS_TABLE
   !!----    Deallocating Bond_Length_Table
   !!----
   !!---- 15/04/2019 
   !!
   Module Subroutine Remove_Bonds_Table()

      if (allocated(Bond_Length_Table)) deallocate(Bond_Length_Table)
      Set_BT_Variable =.false.

      return
   End Subroutine Remove_Bonds_Table
   
End SubModule Del_Routines 