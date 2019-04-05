!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Maths) Sort
 Contains
 
    !!----
    !!---- SORT_PR_PARTITION
    !!----
    !!----    (Private)
    !!----    Utilised by Sort_Strings.
    !!----
    !!---- 03/04/2019 
    !!
    Module Pure Subroutine Sort_PR_Partition(A, Marker)    
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in out) :: A
       integer,                        intent(   out) :: marker

       !---- Local variables ----!
       integer                  :: i, j
       character(len=len(A(1))) :: temp
       character(len=len(A(1))) :: x      ! pivot point

       x = A(1)
       i= 0
       j= size(A) + 1

       do
          j = j-1
          do
             if (A(j) <= x) exit
             j = j-1
          end do
          i = i+1
          do
             if (A(i) >= x) exit
             i = i+1
          end do
          if (i < j) then
             !---- exchange A(i) and A(j)
             temp = A(i)
             A(i) = A(j)
             A(j) = temp
          else if (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          end if
       end do

       return
    End Subroutine Sort_PR_Partition
 
    !!----
    !!---- SORT_I
    !!----    Sort an array such the arr(indx(j)) is in ascending
    !!----    order for j=1,2,...,N.
    !!----
    !!---- 03/04/2019 
    !!
    Module Function Sort_I(arr,n) Result(indx)    
       !---- Arguments ----!
       integer, dimension(:), intent(in ) :: arr       ! Vector
       integer              , intent(in ) :: n         ! Dimension
       integer, dimension(n)              :: indx      ! Index

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer                      :: i,indxt,ir,itemp,j,jstack,k,l
       integer                      :: a


       !> Init
       do j=1,n
          indx(j)=j
       end do

       istack=0
       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if (arr(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l+1)) > arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
             do
                i=i+1
                if (arr(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (arr(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                err_cfml%ierr=1
                err_cfml%msg=" NSTACK too small in SORT_I procedure"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Function Sort_I
 
    !!---- 
    !!---- SORT_R
    !!----
    !!----    (OVERLOADED)
    !!----    Sort an array such the arr(indx(j)) is in ascending
    !!----    order for j=1,2,...,N.
    !!----
    !!---- 03/04/2019 
    !!
    Module Function Sort_R(arr,n) Result(indx)    
       !---- Arguments ----!
       real(kind=cp),dimension(:), intent(in) :: arr       ! Input Vector
       integer,                    intent(in) :: n         ! Dimension
       integer,      dimension(n)             :: indx

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer                      :: i,indxt,ir,itemp,j,jstack,k,l
       real(kind=cp)                :: a

       !> Init
       do j=1,n
          indx(j)=j
       end do

       istack=0
       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if (arr(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l+1)) > arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
             do
                i=i+1
                if (arr(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (arr(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                err_cfml%ierr=1
                err_cfml%msg=" NSTACK too small in SORT_R procedure"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Function Sort_R
 
    !!---- 
    !!---- SORT_STRINGS
    !!----    Sort an array of string
    !!----
    !!---- 03/04/2019 
    !!
    Module Recursive Subroutine Sort_Strings(Str)    
       !---- Argument ----!
       character(len=*), dimension(:), intent(in out) :: Str

       !---- Local variables ----!
       integer :: iq

       if (size(Str) > 1) then
          call Sort_PR_Partition(Str, iq)
          call Sort_Strings(Str(:iq-1))
          call Sort_Strings(Str(iq:))
       end if

       return
    End Subroutine Sort_Strings
 
End Submodule Sort
