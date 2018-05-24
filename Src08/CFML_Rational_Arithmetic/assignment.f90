Submodule (CFML_Rational_Arithmetic) Assignment

 Contains
    !!----
    !!---- FUNCTION ASSIGN_RATIONAL_INT_IL 
    !!----
    !!----
    !! 
    Module Pure Subroutine Assign_Rational_Int_IL(Res, I)
       !---- Arguments ----!
       type(rational),   intent (out) :: res  ! volatile
       integer(kind=il), intent (in)  :: i
      
       res = i // 1_il
       
       return
    End Subroutine Assign_Rational_Int_IL
    
    !!----
    !!---- FUNCTION ASSIGN_RATIONAL_INT_IT 
    !!----
    !!----
    !! 
    Module Pure Subroutine Assign_Rational_Int(Res, I)
       !---- Arguments ----!
       type(rational),  intent (out) :: res  ! volatile
       integer,         intent (in)  :: i
      
       res = i // 1
       
       return 
    End Subroutine Assign_Rational_Int

    !!----
    !!---- SUBROUTINE ASSIGN_RATIONAL_REAL_CP 
    !!----
    !!----
    !! 
    Module Elemental Subroutine Assign_Rational_Real_Cp(Res, Xr)
       !---- Arguments ----!
       type(rational), intent(out) :: res  ! volatile
       real(kind=cp),  intent (in) :: xr
      
       !---- Local variables ----!      
       integer(kind=il)                     :: maxden,ai,t,si
       integer(kind=il), dimension(0:1,0:1) :: m
       real(kind=cp)                        :: x,rai,eps !,startx,er1,er2

       maxden=maximum_denominator; eps=1.0e-6_cp
       m = 0; m(0,0)=1_il; m(1,1)=1_il
       si=sign(1.0,xr)
       x=abs(xr)
      
       !startx=x
       !First (more precise) option
       do
          ai=int(x)
          if ( m(1,0)*ai+m(1,1) > maxden) exit
          t = m(0,0) * ai + m(0,1)
          m(0,1) = m(0,0); m(0,0) = t
          t = m(1,0) * ai + m(1,1)
          m(1,1) = m(1,0)
          m(1,0) =  t
          rai=real(ai,kind=cp)
          if ( abs(x - rai) < eps) exit !division by zero
          x = 1.0_cp/(x - rai)
       end do
      
       res= si*m(0,0)// m(1,0)

       !er1=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)
       !Second option
       ! ai = (maxden - m(2,2)) / m(1,2)
       ! m(1,1) = m(1,1) * ai + mm(2,1)
       ! m(1,2) = mm(1,2) * ai + mm(2,2)
       ! er2=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)

       return
    End Subroutine Assign_Rational_Real_Cp
    
    !!----
    !!---- SUBROUTINE ASSIGN_RATIONAL_REAL_DP 
    !!----
    !!----
    !! 
    Module Elemental Subroutine Assign_Rational_Real_Dp(Res,Xr)
       !---- Arguments ----!
       type(rational), intent(out) :: res
       real(kind=dp),  intent (in) :: xr
      
       !---- Local variables ----!
       integer(kind=il)                      :: maxden,ai,t,si
       integer(kind=il), dimension(0:1,0:1)  :: m
       real(kind=dp)                         :: x,rai,eps

       maxden=maximum_denominator; eps=1.0e-8_dp
       m = 0; m(0,0)=1_il; m(1,1)=1_il
       si=sign(1.0_dp,xr)
       x=abs(xr)

       do
          ai=int(x)
          if ( m(1,0)*ai+m(1,1) > maxden) exit
          t = m(0,0) * ai + m(0,1)
          m(0,1) = m(0,0); m(0,0) = t
          t = m(1,0) * ai + m(1,1)
          m(1,1) = m(1,0)
          m(1,0) =  t
          rai=real(ai,kind=dp)
          if ( abs(x - rai ) < eps) exit !division by zero
          x = 1.0_dp/(x - rai)
       end do
       res= si*m(0,0)// m(1,0)
       
       return
    End Subroutine Assign_Rational_Real_Dp
    
    !!----
    !!---- SUBROUTINE ASSIGN_INT_RATIONAL 
    !!----
    !!----
    !! 
    Module Elemental Subroutine Assign_Int_Rational(I, Res)
       !---- Arguments ----!
       type(rational), intent (in)   :: res  !, volatile
       integer,        intent (out)  :: i
      
       i= nint(real(res%numerator,kind=dp)/real(res%denominator,kind=dp))
       
       return
    End Subroutine Assign_Int_Rational
    
    !!----
    !!---- SUBROUTINE ASSIGN_INTIL_RATIONAL 
    !!----
    !!----
    !! 
    Module Elemental Subroutine Assign_Intil_Rational(I, Res)
       !---- Arguments ----!
       type(rational),  intent (in)   :: res  !, volatile
       integer(kind=il),intent (out)  :: i
      
       i= nint(real(res%numerator,kind=dp)/real(res%denominator,kind=dp))
       
       return
    End Subroutine Assign_Intil_Rational

    !!----
    !!---- SUBROUTINE ASSIGN_REAL_RATIONAL_CP 
    !!----
    !!----
    !! 
    Module Elemental Subroutine Assign_Real_Rational_Cp(X, Res)
       !---- Arguments ----!
       type(rational), intent(in)   :: res
       real(kind=cp),  intent (out) :: x
      
       x=real(res%numerator,kind=cp)/real(res%denominator,kind=cp)
       
       return
    End Subroutine Assign_Real_Rational_Cp
    
    !!----
    !!---- SUBROUTINE ASSIGN_REAL_RATIONAL_DP 
    !!----
    !!----
    !!
    Module Elemental Subroutine Assign_Real_Rational_Dp(X, Res)
       !---- Arguments ----!
       type(rational), intent(in)   :: res
       real(kind=dp),  intent (out) :: x
      
       x=real(res%numerator,kind=dp)/real(res%denominator,kind=dp)
    
       return   
    End Subroutine Assign_Real_Rational_Dp

End Submodule Assignment