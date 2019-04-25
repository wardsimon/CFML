!!----
!!----
!!----
!!
SubModule (CFML_SpaceG) SPG_011
   Contains
   
   !!----
   !!---- GET_SYMB_FROM_OP
   !!----
   !!---- 19/04/2019 
   !!
   Module Function Get_Symb_from_Op(Op, Strcode) Result(symb)
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
   End Function Get_Symb_from_Op

End SubModule SPG_011
   
