!!----
!!----
!!----
!!
SubModule (CFML_Groups) CFML_Grp_014
   Contains
   
   !!----
   !!---- SYMBOL_OPERATOR
   !!----
   !!---- 19/04/2019 
   !!
   Module Function Get_Symb_from_Oper(Op, Strcode) Result(symb)
      !---- Arguments ----!
      type(Symm_Oper_Type),       intent(in) :: Op
      character(len=*), optional, intent(in) :: Strcode
      character(len=80)                      :: symb
      
      !> Init
      symb=" "
      if (present(strcode)) then
         Symb=Get_Symb_from_Mat(Op%Mat,strcode,Op%time_inv)
      else
         Symb=Get_Symb_from_Mat(Op%Mat,"xyz",Op%time_inv)
      end if
      
      return
   End Function Get_Symb_from_Oper

End SubModule CFML_Grp_014   
   
