!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) CFML_M3D_01
 Contains
 
    !!----  SUBROUTINE MATRIX_DIAGEIGEN
    !!----
    !!----    Diagonalize the matrix Array, put eigenvalues in EigenValues and
    !!----    eigenvectors in EigenVec
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine Matrix_DiagEigen(Array,EigenValues,EigenVec)    
       !---- Arguments ----!
       real(kind=cp), intent(in)  , dimension(3,3)    :: array
       real(kind=cp), intent(out) , dimension(3)      :: EigenValues
       real(kind=cp), intent(out) , dimension(3,3)    :: EigenVec

       !---- Local Variables ----!
       integer, parameter            :: N=3, ITMAX=50
       real(kind=cp), parameter      :: EPS1=1.e-7 , EPS2=1.e-7 , EPS3=1.e-7

       integer                       :: i, j, k, nm1, ip1, iter
       real(kind=cp), dimension(3)   :: u
       real(kind=cp), dimension(3,3) :: e
       real(kind=cp)                 :: sigma1, offdsq, p, q, spq, csa, sna
       real(kind=cp)                 :: holdik, holdki, sigma2

       !> Init
       call clear_error()

       nm1=n-1
       do i=1,n
          do j=1,n
             e(i,j)=array(i,j)
             EigenVec(i,j)=0.0_cp
             if (j < i) e(i,j)=0.0_cp
          end do
       end do
       sigma1=0.0_cp
       offdsq=0.0_cp

       do i=1,n
          sigma1=sigma1+e(i,i)**2
          EigenVec(i,i)=1.0_cp
          ip1=i+1
          if (i >= n) exit
          do j=ip1,n
             offdsq=offdsq+e(i,j)**2
          end do
       end do

       do iter=1,itmax
          do i=1,nm1
             ip1=i+1
             do j=ip1,n
                q=abs(e(i,i)-e(j,j))
                if (q <= eps1) then
                   csa=1.0/sqrt(2.0)
                   sna=csa
                else
                   if (abs(e(i,j)) <= eps2) then
                      e(i,j)=0.0_cp
                      cycle
                   end if
                   p=2.0*e(i,j)*q/(e(i,i)-e(j,j))
                   spq=sqrt(p*p+q*q)
                   csa=sqrt((1.0+q/spq)/2.0)
                   sna=p/(2.0*csa*spq)
                end if
                do k=1,n
                   holdki=EigenVec(k,i)
                   EigenVec(k,i)=holdki*csa + EigenVec(k,j)*sna
                   EigenVec(k,j)=holdki*sna - EigenVec(k,j)*csa
                end do
                do k=i,n
                   if (k > j) then
                      holdik=e(i,k)
                      e(i,k)=csa*holdik + sna*e(j,k)
                      e(j,k)=sna*holdik - csa*e(j,k)
                   else
                      u(k)=e(i,k)
                      e(i,k)=csa*u(k)+sna*e(k,j)
                      if (k /= j) cycle
                      e(j,k)=sna*u(k)-csa*e(j,k)
                   end if
                end do
                u(j)=sna*u(i)-csa*u(j)
                do k=1,j
                   if (k <= i)  then
                      holdki=e(k,i)
                      e(k,i)=csa*holdki+sna*e(k,j)
                      e(k,j)=sna*holdki-csa*e(k,j)
                   else
                      e(k,j)=sna*u(k)-csa*e(k,j)
                   end if
                end do
                e(i,j)=0.0_cp
             end do
          end do
          sigma2=0.0_cp
          do i=1,n
             EigenValues(i)=e(i,i)
             sigma2=sigma2+EigenValues(i)*EigenValues(i)
          end do
          if (1.0-sigma1/sigma2 <= eps3) return
          sigma1=sigma2
       end do

       !> Error zone
       Err_CFML%state=.true.
       err_cfml%flag=2
       err_cfml%msg =" Convergence not reached in diagonalization "

       return
    End Subroutine Matrix_DiagEigen
 
End Submodule CFML_M3D_01
