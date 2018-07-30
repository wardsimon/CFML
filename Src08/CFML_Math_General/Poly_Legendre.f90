!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) CFML_MG_08
 Contains
 
    !!---- FUNCTION POLY_LEGENDRE
    !!----
    !!----    Compute the Associated Legendre Polynomial Pml(x).
    !!----
    !!---- Here m (order) and l (degree) are integers satisfying
    !!----    0 <= m <= l
    !!----   -1 <= x <= 1.
    !!----
    !!---- Update: 11/07/2015
    !!
    Module Elemental Function Poly_Legendre(l,m,x) result(plmx)    
       !---- Arguments ----!
       integer,      intent (in) :: l
       integer,      intent (in) :: m
       real(kind=cp),intent (in) :: x
       real(kind=cp)             :: plmx

       !---- Local variables ----!
       integer       :: i, ll
       real(kind=cp) :: fact, pll, pmm, pmmp1, somx2


       !> Check limits
       if (m < 0 .or. m > l .or. abs(x) > 1.0) then
          plmx=0.0
          return
       end if

       !> Calculation
       pmm=1.0
       if (m > 0) then                  !Compute P(m,m)
          somx2=sqrt((1.0_cp-x)*(1.0_cp+x))
          fact=1.0_cp
          do i=1,m
             pmm=-pmm*fact*somx2
             fact=fact+2.0_cp
          end do
       end if

       if (l == m) then
          plmx=pmm
       else
          pmmp1=x*real(2*m+1,kind=cp)*pmm           !Compute P(m,m+1)
          if (l == m+1) then
             plmx=pmmp1                 !Compute P(m,L), L > m+1
          else
             do ll=m+2,l
                pll=(x*real(2*ll-1,kind=cp)*pmmp1-real(ll+m-1,kind=cp)*pmm)/real(ll-m,kind=cp)
                pmm=pmmp1
                pmmp1=pll
             end do
             plmx=pll
          end if
       end if

       return
    End Function Poly_Legendre
    
End Submodule CFML_MG_08