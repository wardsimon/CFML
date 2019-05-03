!!----
!!----
!!----
!!
SubModule (CFML_SpaceG) SPG_015
   Contains
   
   !!----
   !!---- GET_MULTIP_OP_TABLE
   !!----
   !!----   This subroutine construct the Cayley table of a Group
   !!----   defined by the rational operators Op. It is assumed that
   !!----   the first element is the identity.
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Get_Multip_OP_Table(Op,Table)
      !---- Arguments ----!
      type(Symm_Oper_Type), dimension(:),   intent(in) :: Op
      integer, dimension(:,:), allocatable, intent(out):: Table
      
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
   End Subroutine Get_Multip_OP_Table

End SubModule SPG_015   
   
