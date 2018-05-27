!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) CFML_MG_01
 Contains
 
    !!---- SUBROUTINE MEDIAN_QS
    !!----
    !!---- Subroutine calculating the median of a real array
    !!----      Find the median of X(1), ... , X(N), using as much of the quicksort
    !!----      algorithm as is needed to isolate it.
    !!----
    !!----      N.B. On exit, the array x is partially ordered.
    !!----      Based in Alan Miller's median.f90 code.
    !!----
    !!
    Module Pure Subroutine Median_QS(x, n, xmed)    
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out) :: x      ! In: Vector  Out: Sorted vector
       integer,                     intent(in)     :: n      ! Number of data in X
       real(kind=cp),               intent(out)    :: xmed   ! Media of consided data

       !---- Local variables ----!
       real    :: temp, xhi, xlo, xmax, xmin
       logical :: odd
       integer :: hi, lo, nby2, nby2p1, mid, i, j, k

       nby2 = n / 2
       nby2p1 = nby2 + 1
       odd = .true.

       !> HI & LO are position limits encompassing the median.
       if (n == 2 * nby2) odd = .false.
       lo = 1
       hi = n
       if (n < 3) then
          if (n < 1) then
             xmed = 0.0
             return
          end if
          xmed = x(1)
          if (n == 1) return
          xmed = 0.5*(xmed + x(2))
          return
       end if

       !> Find median of 1st, middle & last values.
       do
          mid = (lo + hi)/2
          xmed = x(mid)
          xlo = x(lo)
          xhi = x(hi)
          if (xhi < xlo) then          ! swap xhi & xlo
             temp = xhi
             xhi = xlo
             xlo = temp
          end if
          if (xmed > xhi) then
             xmed = xhi
          else if (xmed < xlo) then
             xmed = xlo
          end if

          !> The basic quicksort algorithm to move all values <= the sort key (XMED)
          !> to the left-hand end, and all higher values to the other end.
          i = lo
          j = hi
          do
              do
                 if (x(i) >= xmed) exit
                 i = i + 1
              end do
              do
                 if (x(j) <= xmed) exit
                 j = j - 1
              end do
              if (i < j) then
                 temp = x(i)
                 x(i) = x(j)
                 x(j) = temp
                 i = i + 1
                 j = j - 1
                 !> Decide which half the median is in.
                 if (i <= j) cycle
              end if
              exit
          end do
          if (.not. odd) then
             if (j == nby2 .and. i == nby2p1) then
                !> Special case, N even, J = N/2 & I = J + 1, so the median is
                !> between the two halves of the series.   Find max. of the first
                !> half & min. of the second half, then average.
                xmax = x(1)
                do k = lo, j
                   xmax = max(xmax, x(k))
                end do
                xmin = x(n)
                do k = i, hi
                   xmin = min(xmin, x(k))
                end do
                xmed = 0.5*(xmin + xmax)
                return
             end if
             if (j < nby2) lo = i
             if (i > nby2p1) hi = j
             if (i /= j) then
                if (lo < hi - 1) cycle
                exit
             end if

             if (i == nby2) lo = nby2
             if (j == nby2p1) hi = nby2p1
          else
             if (j < nby2p1) lo = i
             if (i > nby2p1) hi = j
             if (i /= j) then
                if (lo < hi - 1) cycle
                exit
             end if
             !> test whether median has been isolated.
             if (i == nby2p1) return
          end if
          if (lo < hi - 1) cycle
          exit
       end do

       if (.not. odd) then
          xmed = 0.5*(x(nby2) + x(nby2p1))
          return
       end if
       temp = x(lo)
       if (temp > x(hi)) then
          x(lo) = x(hi)
          x(hi) = temp
       end if
       xmed = x(nby2p1)

       return
    End Subroutine Median_QS
 
End Submodule CFML_MG_01
