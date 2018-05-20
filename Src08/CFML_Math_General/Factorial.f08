!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Factorial
 Contains
 
    Module Elemental Function Factorial(N) Result(Fact)    
       !---- Argument ----!
       integer, intent(in) :: N        ! Factorial of N
       integer             :: Fact

       !---- Local variables ----!
       integer, parameter :: N_LIM = 12
       integer            :: i

       !> Check the current limits
       if ( n > N_LIM) then
          Fact=0
          return
       end if

       if (n ==0) then
          Fact=1
       else
          Fact=1
          do i=1,n
             Fact=Fact*i
          end do
       end if

       return
    End Function Factorial
    
    !!---- FUNCTION FACTORIAL_DP
    !!----
    !!----    Factorial of N but the value returned is a double
    !!----
    !!---- Update: January - 2017
    !!
    Module Elemental Function Factorial_DP(N) Result(Fact)    
       !---- Arguments ----!
       integer, intent(in) :: N      ! Factorial of N
       real(kind=dp)       :: Fact

       !---- Local Variables ----!
       integer       :: i

       if (n == 0) then
          Fact = 1.0_dp
       else
          Fact=1.0_dp
          do i=1,n
             Fact = Fact * real(i,kind=dp)
          end do
       end if

       return
    End Function Factorial_DP
    
    !!---- FUNCTION FACTORIAL_SP
    !!----
    !!----    Factorial of N but the value returned is a real number
    !!----
    !!---- Update: January - 2017
    !!
    Module Elemental Function Factorial_SP(N) Result(Fact)    
       !---- Arguments ----!
       integer, intent(in) :: N    ! Factorial of N
       real(kind=sp)       :: Fact

       !---- Local Variables ----!
       integer       :: i

       if (n == 0) then
          Fact = 1.0_sp
       else
          Fact=1.0_sp
          do i=1,n
             Fact = Fact * real(i)
          end do
       end if

       return
    End Function Factorial_SP
 
End Submodule Factorial
