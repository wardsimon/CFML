!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Sph_Jn
 
 Contains
 
    !!---- SUBROUTINE SPH_JN
    !!----
    !!----    Compute spherical Bessel functions jn(x) and their derivatives
    !!----
    !!---- Update: 11/07/2015
    !!
    Module Subroutine Sph_Jn(n,x,nm,jn,djn)    
       !---- Arguments ----!
       integer,                       intent(in)  :: n   !Order of jn(x) (n=0,1,2,3,...)
       real(kind=dp),                 intent(in)  :: x   !Argument of jn(x)
       integer,                       intent(out) :: nm  !Highest order computed
       real(kind=dp), dimension(0:n), intent(out) :: jn  !array with spherical Bessel functions jn(x)
       real(kind=dp), dimension(0:n), intent(out) :: djn !array with derivatives jn'(x)

       !---- Local variables ----!
       integer       :: k,m
       real(kind=dp) :: sa,sb, f,f1,f0, cs

       nm=n
       if (abs(x) <= 1.0e-30_dp) then
          do k=0,n
             jn(k) = 0.0_dp
             djn(k)= 0.0_dp
          end do
          jn(0)=1.0_dp
          djn(1)=1.0_dp/3.0_dp
          return
       end if

       jn(0)=sin(x)/x
       jn(1)=(jn(0)-cos(x))/x
       if (n >= 2) then
          sa=jn(0)
          sb=jn(1)
          m=Sphjn_PR_Start1(x,200)
          if (m < n) then
             nm=m
          else
             m=Sphjn_PR_Start2(x,n,15)
          end if
          f0=0.0_dp
          f1=1.0e-100_dp
          do k=m,0,-1
             f=(2.0_dp*k+3.0_dp)*f1/x-f0
             if (k <= nm) jn(k)=f
             f0=f1
             f1=f
          end do
          if (abs(sa) > abs(sb)) cs=sa/f
          if (abs(sa) <= abs(sb)) cs=sb/f0
          do k=0,nm
             jn(k)=cs*jn(k)
          end do
       end if

       djn(0)=(cos(x)-sin(x)/x)/x
       do k=1,nm
          djn(k)=jn(k-1)-(k+1.0_dp)*jn(k)/x
       end do

       return
    End Subroutine Sph_Jn
    
    !!--++ FUNCTION SPHJN_PR_ENVJ
    !!--++
    !!--++    (PRIVATE)
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Module Pure Function Sphjn_PR_Envj(N,X) Result(Y)    
       !---- Arguments ----!
       integer,       intent(in) :: n
       real(Kind=dp), intent(in) :: x
       real(Kind=dp)             :: y

       y=0.5_dp*log10(6.28_dp*real(n,kind=dp))-n*log10(1.36_dp*x/real(n,kind=dp))

       return
    End Function Sphjn_PR_Envj
 
    !!--++ FUNCTION SPHN_PR_START1
    !!--++
    !!--++    (PRIVATE)
    !!--++    Determine the starting point for backward
    !!--++    recurrence such that the magnitude of Jn(x) at that point is
    !!--++    about 10^(-MP).
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Module Pure Function Sphjn_PR_Start1(X,Mp) Result (Start)    
       !---- Arguments ----!
       real(kind=dp), intent(in) :: x         ! Argument of Jn(x)
       integer, intent(in)       :: mp        ! Value of magnitude
       integer                   :: start

       !---- Local variables ----!
       integer      :: n1,n0,nn, it
       real(kind=dp):: f,f0,f1,a0

       a0=abs(x)
       n0=int(1.1_dp*a0)+1
       f0=Sphjn_Pr_Envj(n0,a0)-mp
       n1=n0+5
       f1=Sphjn_Pr_Envj(n1,a0)-mp
       do it=1,20
          nn=n1-(n1-n0)/(1.0_dp-f0/f1)
          f=Sphjn_Pr_Envj(nn,a0)-mp
          if (abs(nn-n1) < 1) exit
          n0=n1
          f0=f1
          n1=nn
          f1=f
       end do
       start=nn

       return
    End Function Sphjn_PR_Start1
 
    !!--++ FUNCTION SPHJN_PR_START2
    !!--++
    !!--++    (PRIVATE)
    !!--++    Determine the starting point for backward
    !!--++    recurrence such that all Jn(x) has MP significants digits
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Module Pure Function Sphjn_PR_Start2(X,N,Mp) Result(Start)    
       !---- Arguments ----!
       real(kind=dp), intent(in) :: x     ! Argument of Jn(x)
       integer,       intent(in) :: n     ! Order of Jn(x)
       integer,       intent(in) :: mp    ! Significant digit
       integer                   :: start

       !---- Local variables ----!
       real(kind=dp) :: a0, hmp, ejn, obj,f,f0,f1
       integer       :: n0,n1,nn, it

       a0=abs(x)
       hmp=0.5_dp*mp
       ejn=Sphjn_Pr_Envj(n,a0)
       if (ejn <= hmp) then
          obj=mp
          n0=int(1.1_dp*a0)+1  ! +1 was absent in the original version ... this solves the problem of
       else                    ! Intel, gfortran and g95 compilers ... Lahey was calculating well event if n0=0!
          obj=hmp+ejn
          n0=n
       end if
       f0=Sphjn_Pr_Envj(n0,a0)-obj
       n1=n0+5
       f1=Sphjn_Pr_Envj(n1,a0)-obj
       do it=1,20
          nn=n1-(n1-n0)/(1.0_dp-f0/f1)
          f=Sphjn_Pr_Envj(nn,a0)-obj
          if (abs(nn-n1) < 1) exit
          n0=n1
          f0=f1
          n1=nn
          f1=f
       end do
       start=nn+10

       return
    End Function Sphjn_PR_Start2
 
End Submodule Sph_Jn
