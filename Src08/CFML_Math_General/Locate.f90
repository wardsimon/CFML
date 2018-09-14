!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Locate
 Contains

    !!--++ FUNCTION LOCATE_I
    !!--++
    !!--++    Subroutine for locating the index J of an array V(:)
    !!--++    satisfying:
    !!--++
    !!--++               V(J) <= X < V(J+1)
    !!--++
    !!--++ Update: June - 2011
    !!
    Module Function Locate_I(V,x) Result(j)
       !---- Argument ----!
       integer, dimension(:), intent(in):: v  ! Input vector
       integer,               intent(in):: x  ! Value
       integer                          :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm, i1,i2

       i1=lbound(v,dim=1)
       i2=ubound(v,dim=1)

       if (x <= v(i1)) then
          j=i1
          return
       end if
       if(x >= v(i2)) then
         j=i2
         return
       end if

       jl=i1-1
       ju=i2+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((v(i2) > v(i1)) .eqv. (x > v(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_I

    !!--++ FUNCTION LOCATE_R
    !!--++
    !!--++    Function for locating the index J of an array XX(:)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++ Update: June - 2011
    !!
    Module Function Locate_R(V,x) Result(j)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in):: v
       real(kind=cp),               intent(in):: x
       integer                                :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm, i1,i2

       i1=lbound(v,dim=1)
       i2=ubound(v,dim=1)

       if (x <= v(i1)) then
         j=i1
         return
       end if
       if(x >= v(i2)) then
         j=i2
         return
       end if

       jl=i1-1
       ju=i2+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((v(i2) > v(i1)) .eqv. (x > v(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_R


End Submodule Locate
