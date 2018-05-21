!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Co_Prime
 Contains
 
    !!---- FUNCTION CO_PRIME
    !!----
    !!---- Provides the value .TRUE. if the array V contains co-prime
    !!----
    !!---- integers: there is no common divisor for all the integers.
    !!---- Only the first 1000 prime numbers are stored in the module array "primes"
    !!---- imax is the maximum prime number to be tested. It is calculated if not given.
    !!----
    !!
    Module Function Co_Prime(v,imax) result(cop)    
       !---- Arguments ----!
       integer, dimension(:),           intent(in) :: v          ! Input vector of numbers
       integer,               optional, intent(in) :: imax       ! Maximun prime number to be tested
       Logical                                     :: cop

       !---- Local variables ----!
       integer :: i,j,im,k,dimv,imaxv,maxv

       !> init
       cop=.true.

       maxv=maxval(abs(v))
       if (present(imax)) then
          imaxv=imax
       else
          imaxv=maxv
       end if

       !---- If the maximum value of the indices is 1 they are coprimes
       if (maxv == 1) return
       if (maxv == 0) then
          cop=.false.
          return
       end if

       !---- Search the maximum prime number to be tested
       if (imaxv > 7919) then
          im=1000
       else
          do i=1,1000
             if (imaxv > primes(i)) cycle
             im=i
             exit
          end do
       end if

       !---- Indices greater than 1
       dimv=size(v)
       do_p: do i=1,im
          k=primes(i)
          do j=1,dimv
             if( mod(v(j),k) /= 0) cycle do_p
          end do
          cop=.false.
          exit
       end  do do_p

       return
    End Function Co_Prime
    
    !!---- SUBROUTINE CO_PRIME_VECTOR
    !!----
    !!----     Calculates the co-prime vector (cop) parallel to the input vector (v)
    !!----     It uses the list of the first thousand prime numbers.
    !!----
    !!----   Updated: January 2012 (JRC), copied from Nodal_Indices (Laue_Mod) in July 2013 (JRC)
    !!----
    Module Subroutine Co_Prime_Vector(V,Cop,f)    
       !---- Arguments ----!
       integer, dimension(:),           intent(in)  :: v              !input integer vector
       integer, dimension(:),           intent(out) :: cop            !Output co-prime vector
       integer,               optional, intent(out) :: f              !Common multiplicative factor

       !---- Local variables ----!
       integer :: i,j,max_ind,k,im,dimv,n

       cop=v
       n=1
       if (present(f)) f=1
       max_ind=maxval(abs(cop))

       !> If the maximum value of the indices is 1 they are already coprimes
       if (max_ind <= 1) return

       !>- Indices greater than 1
       dimv=size(v)
       im=0
       do i=1,size(primes)
          if (primes(i) > max_ind) then  !primes is an array within this module
             im=i
             exit
          end if
       end do
       if(im == 0) return

       do_p: do i=1,im
          k=primes(i)
          do
             do j=1,dimv
                if( mod(cop(j),k) /= 0) cycle do_p
             end do
             n=n*k
             cop=cop/k
          end do
       end do do_p

       if (present(f)) f=n

       return
    End Subroutine Co_Prime_vector
 
End Submodule Co_Prime
