!!----
!!----
!!----
!!
SubModule (CFML_Groups) CFML_GRP_024
   Contains
   
   !!----
   !!---- GET_GROUP_FROM_TABLE
   !!----
   !!---- 20/04/19  TO BE TESTED
   !!
   Module Subroutine Get_Group_From_Table(N, G, Table, Ord)  
      !---- Arguments ----! 
      integer,                 intent(in)     :: n   !Number of initial elements in G
      integer, dimension(:),   intent(in out) :: G   !Vector containing the different elements of the groups
      integer, dimension(:,:), intent(in)     :: Table
      integer,                 intent(out)    :: ord !order or the final group
      
      !---- Local Variables ----!
      integer :: i,j,k,m,nt,mult,ij
      logical, dimension(size(Table,1),size(Table,1)) :: done
      
      mult=size(Table,1)
      done=.false.

      nt=n
      do_ext:do
         m=nt
         do i=1,m
            do_j:do j=1,m
               if (done(i,j)) cycle do_j
               ij=table(i,j)
               do k=1,nt
                  if (ij == G(k)) then
                     done(i,j)=.true.
                     cycle do_j
                  end if
               end do
               
               done(i,j)=.true.
               nt=nt+1
               G(nt)=ij
               if (nt > mult) then
                  nt=nt-1
                  exit do_ext
               end if
            end do do_j
         end do  !i
         if ( m == nt) exit do_ext
      end do do_ext
      ord=nt
      
      return
   End Subroutine Get_Group_From_Table

End SubModule CFML_GRP_024   
   
