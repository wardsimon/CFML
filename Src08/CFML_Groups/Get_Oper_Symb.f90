!!----
!!----
!!----
!!
SubModule (CFML_Groups) CFML_Grp_015
   Contains
   
   !!----
   !!---- OPERATOR_FROM_SYMBOL
   !!----
   !!---- 19/04/19
   !!
   Module Function Get_Oper_from_Symb(Symb) Result(Op)
      !---- Arguments ----!
      character(len=*),     intent(in) :: symb
      type(Symm_Oper_Type)             :: Op
      
      !---- Local variables ----!
      integer :: d
      
      d=Get_Dimension_Gener(Symb)
      
      allocate(Op%Mat(d,d))
      call Get_Mat_From_Symb(Symb,Op%Mat,Op%time_inv)
      Op%dt=rdet(Op%Mat(1:3,1:3))
      
      return
   End Function Get_Oper_from_Symb

End SubModule CFML_Grp_015   
   
