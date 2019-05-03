SubModule (CFML_SpaceG) Oper_Mult
   Contains
   
   !!----
   !!---- MULTIPLY_SYMM_OPER
   !!----
   !!---- 19/04/2019
   !!
   Module Pure Function Multiply_Symm_Oper(Op1,Op2) Result (Op3)
      !---- Arguments ----!
      type(Symm_Oper_Type), intent(in) :: Op1,Op2
      type(Symm_Oper_Type)             :: Op3
      
      !---- Local Variables ----!
      integer :: n,d,i
      
      n=size(Op1%Mat,dim=1)
      allocate(Op3%Mat(n,n))
      
      d=n-1
      Op3%Mat=matmul(Op1%Mat,Op2%Mat) 
      Op3%Mat(1:d,n)=mod(Op3%Mat(1:d,n),1_LI)
      do i=1,d
         do
            if (Op3%Mat(i,n) < 0_LI//1_LI) then
               Op3%Mat(i,n) = Op3%Mat(i,n) + 1_LI
            else
               exit
            end if
         end do
      end do
      Op3%time_inv=Op1%time_inv*Op2%time_inv
      Op3%dt=Op1%dt*Op2%dt
    End Function Multiply_Symm_Oper
    
End SubModule Oper_Mult