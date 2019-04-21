!!----
!!----
!!----
!!
SubModule (CFML_Groups) CFML_Grp_018
   Contains
   
   !!----
   !!---- GET_MULTIPLICATION_TABLE
   !!----
   !!----   This subroutine construct the Cayley table of a Group
   !!----   defined by the rational operators Op. It is assumed that
   !!----   the first element is the identity.
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Get_Multiplication_Table(Op,Table)
      !---- Arguments ----!
      type(Symm_Oper_Type), dimension(:), intent(in) :: Op
      integer, dimension(:,:),allocatable,intent(out):: Table
      
      !---- Local Variables ----!
      integer:: i,j,m,Multip
      type(Symm_Oper_Type) :: Opm
      
      multip=size(Op)
      allocate(Table(multip,multip))
      Table=0
      Table(1,:) = [(i,i=1,multip)] !It is supposed that the first operator is the identity
      Table(:,1) = [(i,i=1,multip)]
      do i=2,multip
         do j=2,multip
            Opm=Op(i)*Op(j)
            do m=1,multip
               if (Opm == Op(m)) then
                  Table(i,j)=m
                  exit
               end if
            end do
         end do
      end do
      
      return
   End Subroutine Get_Multiplication_Table

End SubModule CFML_Grp_018   
   
