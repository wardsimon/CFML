!!----
SubModule (CFML_gSpaceGroups) SPG_Apply_OP
   implicit none
   Contains
   !!----
   !!---- APPLY_OP
   !!----    Apply a symmetry operator to a vector:  Vp = Apply_OP(Op,v)
   !!----
   !!---- 30/05/2019
   !!
   Module Function Apply_OP(Op, V) Result(S)
      !---- Arguments ----!
      Type(Symm_Oper_Type),         intent(in) :: Op    ! Symmetry Operator
      real(kind=cp), dimension(3),  intent(in) :: v     ! Vector
      real(kind=cp), dimension(3)              :: S     ! Output vector

      !---- Local Variables ----!
      type(rational), dimension(4) :: tr

      tr(1:3)=v
      tr(4)=1_LI//1_LI

      tr=matmul(Op%Mat,Tr)
      S=tr(1:3)
   End Function Apply_OP

End SubModule SPG_Apply_OP