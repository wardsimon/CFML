!!----
!!----
!!----
!!
SubModule (CFML_Profiles) PRF_003
   Contains
   !!----
   !!---- FUNCTION GAUSSIAN
   !!----
   !!----
   !!---- Update: October - 2005
   !!
   Pure Module Function Gaussian(X,Par) Result (Gauss_Val)
      !---- Arguments ----!
      real(kind=cp),              intent(in) :: x
      real(kind=cp), dimension(:),intent(in) :: par
      real(kind=cp)                          :: gauss_val

      !--- Local variables ---!
      real(kind=cp) :: H,ag,bg

      H=par(1)
      ag= 0.93943727869965133377234032841018/H
      bg= 2.7725887222397812376689284858327/(H*H)
      gauss_val = ag* exp(-bg*x*x)

      return
   End Function Gaussian

   !!----
   !!---- Subroutine Gaussian_Der
   !!----
   !!---- Update: October - 2005
   !!
   Pure Module Subroutine Gaussian_Der(X,Par,Gauss_Val,Dpar)
      !---- Arguments ----!
      real(kind=cp),                       intent(in) :: x
      real(kind=cp),          dimension(:),intent(in) :: par
      real(kind=cp),                       intent(out):: gauss_val
      real(kind=cp), optional,dimension(:),intent(out):: dpar

      !--- Local variables ---!
      real(kind=cp) :: H,ag,bg, gaussp,dgaussH,x2

      H=par(1)
      ag= 0.93943727869965133377234032841018/H
      bg= 2.7725887222397812376689284858327/(H*H)
      x2=x*x
      gauss_val = ag* exp(-bg*x2)
      if (present(dpar)) then
         gaussp  = -2.0 * gauss_val * bg * x
         dgaussH = (2.0*bg*x2-1.0)*gauss_val/H
         dpar(1:2)=(/gaussp,dgaussH/)
      end if

      return
   End Subroutine Gaussian_Der

   !!--++
   !!--++ GAUSSIAN
   !!--++    Return value of Gaussian at 'Pos' for peak at 'Pos0' and 'Gamma'
   !!--++
   !!--++ 21/04/2019
   !!
   Module Subroutine Prof_Gaussian(Pos , Pos0 , Gamma , Dgdt , Dgdg, Gauss )
      !---- Arguments ----!
      real(kind=cp), intent(in)  :: pos
      real(kind=cp), intent(in)  :: pos0
      real(kind=cp), intent(in)  :: gamma
      real(kind=cp), intent(out) :: dgdt       !is derivative of G wrt Pos0
      real(kind=cp), intent(out) :: dgdg       !is derivative of G wrt Gamma
      real(kind=cp), intent(out) :: gauss

      !---- Local Variables ----!
      real(kind=cp),parameter    ::   c= 1.6651092223154_cp   ! 2*sqrt(ln2)
      real(kind=cp),parameter    ::  cg= 0.93943727869965_cp  ! 2*sqrt(ln2/pi)
      real(kind=cp)              ::  delp , ex

      delp = pos - pos0
      if (abs(delp)/gamma > 6.0 ) then
         Gauss = 0.0
         dgdt = 0.0
         dgdg = 0.0
      else
         ex=(delp * c /gamma)**2
         Gauss = cg * exp(-ex)/gamma
         dgdg = Gauss * ( -1.0 + 2.0 * ex) / gamma
         dgdt = 2.0 * c*c * delp * Gauss/(gamma*gamma)
      end if

      return
   End Subroutine Prof_Gaussian

End SubModule PRF_003