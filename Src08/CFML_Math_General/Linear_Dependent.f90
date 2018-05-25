!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Linear_Dependent
 Contains
 
    !!--++ SUBROUTINE LINEAR_DEPENDENTC
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++
    !!--++    For the case of complex vectors in Cn the problem can be reduced to real vectors
    !!--++    of dimension R2n. Each complex vector contributes as two real vectors of dimension
    !!--++    2n: (R,I) and (-I,R). A complex vector V is linearly dependent on n complex vectors
    !!--++    if V can be written as: V = Sigma{j=1,n}(Cj.Vj), with Cj complex numbers and Vj
    !!--++    having n complex components. One may write:
    !!--++
    !!--++     V = Sigma{j=1,n}(Cj.Vj)
    !!--++     (R,I) = Sigma{j=1,n} (Cjr Vj + i Cji Vj) = Sigma{j=1,n} (Cjr (Rj,Ij) +  Cji (-Ij,Rj) )
    !!--++     (R,I) = Sigma{j=1,n} (aj (Rj,Ij) + bj (-Ij,Rj) )  = Sigma{j=1,2n} (Aj.Uj)
    !!--++     Were Uj=(Rj,Ij) and U(j+1)= (-Ij,Rj)
    !!--++
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Linear_DependentC(A,na,B,nb,mb,info)    
       !---- Arguments ----!
       complex, dimension(:),   intent(in)  :: a
       complex, dimension(:,:), intent(in)  :: b
       integer,                 intent(in)  :: na,nb,mb
       logical,                 intent(out) :: info

       !---- Local variables ----!
       integer                                                     :: r,n1
       real(kind=dp), parameter                                    :: tol= 100.0_dp*deps
       real(kind=dp), dimension(2*max(nb+1,mb+1),2*max(nb+1,mb+1)) :: c

       c=0.0_dp
       call clear_error()

       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_CFML%state=.true.
          err_cfml%flag=2
          err_cfml%msg=" Linear_DependentC: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then
          n1=2*nb+1
          if(n1+1 > 2*mb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(b(1:nb,1:mb))
          c(1:nb,     mb+1:mb+na) = aimag(b(1:nb,1:mb))
          c(nb+1:2*nb,      1:mb) =-aimag(b(1:nb,1:mb))
          c(nb+1:2*nb,mb+1:mb+na) =  real(b(1:nb,1:mb))
          c(n1,             1:mb) =  real(a(1:na))
          c(n1,      mb+1:mb+na ) = aimag(a(1:na))
          c(n1+1,           1:mb) =-aimag(a(1:na))
          c(n1+1,    mb+1:mb+na ) =  real(a(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*mb)) info=.false.

       else if( na == nb) then
          n1=2*mb+1
          if(n1+1 > 2*nb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(b(1:nb,1:mb))
          c(nb+1:nb+na,     1:mb) = aimag(b(1:nb,1:mb))
          c(1:nb,      mb+1:2*mb) =-aimag(b(1:nb,1:mb))
          c(nb+1:nb+na,mb+1:2*mb) =  real(b(1:nb,1:mb))
          c(1:na,             n1) =  real(a(1:na))
          c(nb+1:nb+na,       n1) = aimag(a(1:na))
          c(1:na,           1+n1) =-aimag(a(1:na))
          c(nb+1:nb+na,     1+n1) =  real(a(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*nb)) info=.false.

       else
          ERR_CFML%state=.true.
          err_cfml%flag=2
          err_cfml%msg=" Linear_DependentC: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentC
 
    !!--++ SUBROUTINE Linear_DependentI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Linear_DependentI(A,na,B,nb,mb,info)    
       !---- Arguments ----!
       integer, dimension(:),   intent(in)  :: a
       integer, dimension(:,:), intent(in)  :: b
       integer,                 intent(in)  :: na,nb,mb
       logical,                 intent(out) :: info

       !---- Local variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       c=0.0_dp
       call clear_error()
       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_CFML%state=.true.
          err_cfml%flag=2
          err_cfml%msg=" Linear_DependentI: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(b(1:nb,1:mb))
          c(n1,  1:mb)=real(a(1:na))      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) info=.false.

       else if( na == nb) then
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(b(1:nb,1:mb))
          c(1:nb,  n1)=real(a(1:na))     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) info=.false.

       else
          ERR_CFML%state=.true.
          err_cfml%flag=2
          err_cfml%msg=" Linear_DependentI: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentI
 
    !!--++ SUBROUTINE LINEAR_DEPENDENTR
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Linear_DependentR(A,na,B,nb,mb,info)    
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in)  :: a
       real(kind=cp), dimension(:,:), intent(in)  :: b
       integer,                       intent(in)  :: na,nb,mb
       logical,                       intent(out) :: info

       !---- Local Variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       c=0.0_dp
       call clear_error()
       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_CFML%state=.true.
          err_cfml%flag=2
          err_cfml%msg=" Linear_DependentR: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then    !Vector added as an additional row
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=b(1:nb,1:mb)
          c(n1,  1:mb)=a(1:na)      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) info=.false.

       else if( na == nb) then   !Vector added as an additional column
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=b(1:nb,1:mb)
          c(1:nb,  n1)=a(1:na)     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) info=.false.
       else
          ERR_CFML%state=.true.
          err_cfml%flag=2
          err_cfml%msg=" Linear_DependentR: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentR
 
   
End Submodule Linear_Dependent
