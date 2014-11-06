  Module diffax_calc

   use CFML_GlobalDeps,           only : sp
   use CFML_String_Utilities,     only : number_lines, reading_lines,  init_findfmt, findfmt , &
                                         iErr_fmt, getword, err_string, err_string_mess, getnum, Ucase
   use CFML_Optimization_General, only : Opt_Conditions_Type
   use CFML_LSQ_TypeDef,          only : LSQ_Conditions_type
   use read_data,                 only : crys_2d_type, read_structure_file, length, opti, cond
   use diffax_mod

   implicit none
   private
   public ::  salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, &
              optimz,point,  gospec, gostrk, gointr,gosadp, getfnm, get_sym,&
              chk_sym, overlp, nmcoor


   !type (Opt_Conditions_Type),     save  :: opti
   !type (LSQ_Conditions_type),     save  :: cond

   contains
 ! ______________________________________________________________________
! Title: AGLQ16
! Author: MWD
! Date: 18 Aug 1988
! Description:  This routine does adaptive 16-point Gauss-Legendre
! quadrature in an interval (h,k,a) to (h,k,b) of reciprocal space.
! ok is returned as .TRUE. if GLQ16 hasn't blown its stack. The
! integrated result is returned in AGLQ16.

!      ARGUMENTS:
!            h   -  reciprocal lattice vector h-component. (input).
!            k   -  reciprocal lattice vector k-component. (input).
!            a   -  l-value of the lower bound of reciprocal
!                   lattice integration region. (input).
!            b   -  l-value of the upper bound of reciprocal
!                   lattice integration region. (input).
!            ok  -  logical flag indicating all went well. (output).

!      AGLQ16 returns the adaptively integrated value.
! ______________________________________________________________________

      REAL*8 FUNCTION aglq16(h, k, a, b, ok)
!   Utiliza las variables : DIFFaX.par , h, k ,a, b ,  ok ,maxstk= 200, stp, n, n2 , sum
!                           sum1, sum2, sum3, epsilon= FIVE * eps4, epsilon2, GLQ16,  stk(maxstk),
!                           d1, d2, d3, x
!   Utiliza las funciones: GLQ16(h, k, a, b, ok) externa  ,
!   Utiliza las subrutinas:

      Integer,      Intent(In)  :: h
      Integer,      Intent(In)  :: k
      Real(Kind=8), Intent(In)  :: a
      Real(Kind=8), Intent(In)  :: b
      Logical,      Intent(Out) :: ok
      Integer  :: stp, n, n2
      Integer, Parameter :: maxstk = 200
      Real(Kind=8)            :: sum, sum1, sum2, sum3,  epsilon2,  stk(maxstk), d1, d2, d3, x
      Real(Kind=8), Parameter :: epsilon = five * eps4

! external function
      aglq16 = zero
! initalize stack; top is at highest index
      stp = maxstk
! get first integration
      sum3 = glq16(h, k, a, b, ok)


      IF(.NOT.ok) GO TO 999
      n2 = 2
      n = 0
      d1 = a
      d3 = b
      sum = zero
      20 d2 = half * (d1 + d3)
      n = n + 2
      sum1 = glq16(h, k, d1, d2, ok)
      IF(.NOT. ok) GO TO 999
      sum2 = glq16(h, k, d2, d3, ok)
      IF(.NOT. ok) GO TO 999
      x = sum1 + sum2
! determine figure of merit
      epsilon2 = MAX(epsilon, ABS(x) * epsilon)
      IF(ABS(x - sum3) > epsilon2) THEN
! the area of these two panels is not accurately known
! check for stack overflow
        IF(stp < 3) GO TO 980
! push right one onto stack
        stk(stp) = sum2
        stk(stp-1) = d2
        stk(stp-2) = d3
        stp = stp - 3
        d3 = d2
        sum3 = sum1
      ELSE
! this panel has been accurately integrated; save its area
        sum = sum + x
! get next panel
! check for stack underflow--happens when no panels left to pop off
        IF(stp == maxstk) GO TO 30
        d3 = stk(stp+1)
        d1 = stk(stp+2)
        sum3 = stk(stp+3)
        stp = stp + 3
      END IF
      IF(n == n2) THEN
        n2 = n2 * 2
      END IF
      GO TO 20
      30 aglq16 = sum
      RETURN
      980 WRITE(op,402) 'Stack overflow in AGLQ16.'
      RETURN
      999 WRITE(op,402) 'GLQ16 returned an error to AGLQ16.'
      WRITE(op,403) h, k, d1, d2
      RETURN
      402 FORMAT(1X, a)
      403 FORMAT(3X,'at: h = ',i3,' k = ',i3,' l1 = ',g12.5,' l2 = ',g12.5)
      END FUNCTION aglq16

! ______________________________________________________________________
! Title: APPR_F
! Author: MMJT
! Date: 11 Feb 1989
! Description: This subroutine returns polynomial approximations of f
!  at 16 points, for use in GLQ16. The f's are returned in a MAX_L x 16
!  array. The order of the polynomial is n-1, where n is input. The 16
!  points in reciprocal space for which f's are needed are given by
!  (h,k,ag_l(i)), where h, k and ag_l are input.
!  The ll are the n sampling points, whose f's must be calculated
!  by direct calls to GET_F, and from which the interpolations are
!  calculated. list contains the indices of those n points. There is no
!  need to interpolate for these values since we know them already!

!      ARGUMENTS:
!            f      -  Array of layer form factors. (output).
!            h      -  reciprocal lattice vector h-component. (input).
!            k      -  reciprocal lattice vector k-component. (input).
!            ll     -  an array of n l-values to be used in the
!                      interpolation. (input).
!            ag_l   -  an array of 16 l_values at which interpolated
!                      form factors are required. (input).
!            n      -  the order of the polynomial approximation.
!                                                              (input).
!            list   -  A list of the indices of the n ll points
!                      entered. Interpolation is not needed at these
!                      input values. (input).
!            ok     -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:     a0, b0, c0, d0, n_layers
! ______________________________________________________________________

      SUBROUTINE appr_f(f,h,k,ll,ag_l,n,list,ok)
!   Utiliza las variables :  DIFFaX.par,DIFFaX.inc, n, list(n) ,ll(n), ag_l(16)
!                            h, k , f(MAX_L,16) , ok , know_f,  i, j, m, p, max_poly =10 ,
!                            Q2, l  ,ff(MAX_L,max_poly), fa(MAX_L), f_ans, f_error
!
!   Utiliza las funciones :   Q2(h,k,l)
!   Utiliza las subrutinas : POLINT  ,GET_F   ,

      COMPLEX*16, INTENT(IN OUT)               :: f(max_l,16)
      INTEGER*4, INTENT(IN )                   :: h
      INTEGER*4, INTENT(IN )                   :: k
      REAL*8, INTENT(IN OUT)                   :: ll(n)
      REAL*8, INTENT(IN OUT)                   :: ag_l(16)
      INTEGER*4, INTENT(IN)                    :: n
      INTEGER*4, INTENT(IN)                    :: list(n)
      LOGICAL, INTENT(IN OUT)                  :: ok

      LOGICAL :: know_f
      INTEGER*4 i, j, m, p
      INTEGER*4, PARAMETER :: max_poly = 10
      REAL*8 q2, l
      COMPLEX*16 ff(max_l,max_poly), fa(max_l), f_ans, f_error

! external subroutines (Some compilers need them declared external)
!      external POLINT, GET_F

! statement function
! Q2 is the value of 1/d**2 at hkl
      q2(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0

! sample GET_F n times
      DO  i = 1, n
        CALL get_f(ff(:,i), q2(h,k,ll(i)), ll(i))
      END DO

! store the n sampled f values for each layer i in fa.
! call POLINT for each layer in turn
! do this 16 times.
      DO  m = 1, 16
! check to see that we haven't just calculated f at l(m)
        know_f = .false.
        DO  i = 1, n
          IF(m == list(i)) THEN
            p = i
            know_f = .true.
          END IF
        END DO
! if we have, then simply copy it
        IF(know_f) THEN
          DO  i = 1, n_layers
            f(i,m) = ff(i,p)
          END DO
        ELSE
! else, use polynomial interpolation.
          DO  i = 1, n_layers
            DO  j = 1, n
              fa(j) = ff(i,j)
            END DO
            CALL polint(ll,fa,n,ag_l(m),f_ans,f_error,ok)
            IF(.NOT. ok) GO TO 999
            f(i,m) = f_ans
          END DO
        END IF
      END DO

      RETURN
      999 WRITE(op,100) 'POLINT returned an error to APPR_F.'
      RETURN
      100 FORMAT(1X, a)
      END SUBROUTINE appr_f

! ______________________________________________________________________
! Title: ATOMS
! Authors: MMJT
! Date: 18 Feb 90
! Description: This routine lists the legal atom names accepted by
! DIFFaX. These names are taken from the file 'sfname'.

!      ARGUMENTS:
!           No arguments are used.

!      COMMON VARIABLES:
!            uses:  sfname

!        modifies:  No COMMON variables are modified.
! ______________________________________________________________________

      SUBROUTINE atoms()
!   Utiliza las variables :  DIFFaX.par,DIFFaX.inc,  i,  atom_name  , atomline
!                            line,
!   Utiliza las subrutinas :

      INTEGER*4 i
      CHARACTER (LEN=4) :: atom_name
      CHARACTER (LEN=80) :: atomline
      CHARACTER (LEN=120) :: line

! sfname has been opened once already,
! no need for elaborate error checking
      OPEN(UNIT = sf, FILE = sfname, STATUS = 'old', ERR = 999)

      WRITE(op,200) 'Legal atom types are: '

! skip first few lines which contain no data, until we come to Hydrogen
      1 READ(sf,'(a)') line
      IF(line(1:4) /= 'H   ') GO TO 1
      BACKSPACE(UNIT = sf, ERR = 500)

! initialize atomline
      atomline = ' '
! write 10 atom names to a line
      10 i = 0
      20 READ(sf, '(a)', END = 100) line
      atom_name = line(1:4)
      WRITE(atomline(i*7+1:i*7+7),'(3a)') '''', atom_name, ''' '
      i = i + 1
      IF(MOD(i,10) == 0) THEN
        WRITE(op,210) atomline
        GO TO 10
      END IF
      GO TO 20
      100 CONTINUE

      CLOSE(UNIT = sf, ERR = 600)

      999 RETURN
      500 WRITE(op,220) 'Unable to backspace scattering factor file ''',  &
          sfname(1:length(sfname)), '''.'
      RETURN
      600 WRITE(op,220) 'Unable to close scattering factor file ''',  &
          sfname(1:length(sfname)), '''.'
      RETURN
      200 FORMAT(1X, a)
      210 FORMAT(2X, a)
      220 FORMAT(1X, 'ERROR in ATOMS: ', 3A)
      END SUBROUTINE atoms

!!S Bfiles
! ______________________________________________________________________
! Title: BINPOW
! Author: MMJT
! Date: 18 Mar 1990
! Description:  This function breaks down the number 'n' into its
! binary representation. The result is stored in the global array 'pow'.
! This is used for efficiently multiplying a square matrix by itself
! n times if the RECURSIVE option was chosen for a finite number of
! layers.

! n must be such that n <= RCSV_MAX+1 <= 2**(MAX_BIN+1) - 1

!      ARGUMENTS:
!            n   -  number to binarize. (input).

!      COMMON VARIABLES:
!            uses:  max_pow

!        modifies:  max_pow, pow

! BINPOW returns logical TRUE if no problems were encountered.
! ______________________________________________________________________

      LOGICAL FUNCTION binpow(n)
!       Utiliza las variables:'DIFFaX.par , DIFFaX.inc, n , i, j, itmp,
!       Utiliza las funciones:
!       Utiliza las subrutinas:

        REAL*4               :: n
        INTEGER*4 i, j, itmp

        binpow = .false.
        itmp = n
        max_pow = 0
        i = max_bin + 1
        do
          i = i - 1
          j = itmp / 2**(i-1)
          ! j should be either 1 or zero if n was within bounds
          IF(j == 1) THEN
            pow(i) = 1
            itmp = itmp - 2**(i-1)
            IF(max_pow == 0) max_pow = i
          ELSE IF(j == 0) THEN
            pow(i) = 0
          ELSE
            WRITE(op,"(a, i8)") '  ERROR in BINPOW: Invalid exponent ', n
            WRITE(op,"(a, i8)") '  Maximum No. of layers allowed = ', rcsv_max
            WRITE(op,"(a, i8)") '  Maximum supported by DIFFaX is', 2**(max_bin+1)-1
            RETURN
          END IF
          IF(i <= 1) exit
        end do
        binpow = .true.
        RETURN
      END FUNCTION binpow

! ______________________________________________________________________
! Title: BOUNDS
! Authors: MMJT
! Date: 24 Feb 1990
! Description: This function translates the value x so that it lies
! within the range 0 to 1.
!      ARGUMENTS:
!                  x  -  real number to convert. (input).
!      COMMON VARIABLES:
!                  No COMMON variables are used
!      BOUNDS returns the translated value of x. x is not modified.
! ______________________________________________________________________

      REAL*8 FUNCTION bounds(x)
!       Utiliza las variables:'DIFFaX.par , x, y
!       Utiliza las funciones:
!       Utiliza las subrutinas:
        REAL*8, INTENT(IN)    :: x
        REAL*8 y

        y = x - INT(x) + one
        y = y - INT(y)

        ! allow for rounding error
        IF((one - y) < eps5) y = zero
        bounds = y
        RETURN
      END FUNCTION bounds
!________________________________________________________________________
!Title:CALC_BCKG
!
!Linear interpolation of the background
!
!________________________________________________________________________

 ! SUBROUTINE calc_bckg(x,y)

!  REAL*8 ,dimension(:), allocatable   ::  x,y
!!S Cfiles
! ______________________________________________________________________
! Title: CHK_SYM
! Author: MMJT
! Date: 15 Aug 1989; 21 Jan 1995
! Checks the user's assertions in the data file about the symmetry of
! the diffraction pattern. The symmetry is used to define the smallest
! possible volume in reciprocal space over which integration can
! occur and be representative of the whole pattern. CHK_SYM gives the
! user a crude 'goodness of fit', or tolerance, within which a random
! sampling of intensities fit the asserted symmetry. CHK_SYM does
! not override the user's judgment if the fit is poor, unless the cell
! dimensions or angle are inconsistent with the requested symmetry. The
! intention of CHK_SYM is to alert the user that the data does not
! conform to the expected symmetry. Useful for debugging datafiles.

!      ARGUMENTS:
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  cell_a, cell_b, cell_gamma, max_angle, SymGrpNo,
!                   pnt_grp, tolerance, PI, PI2, PS_VGT, RAD2DEG

!        modifies:  max_var, h_mirror, k_mirror, check_sym, SymGrpNo
! ______________________________________________________________________

      SUBROUTINE chk_sym(ok)


!     Utiliza las variables:'DIFFaX.par , DIFFaX.inc,ok , diad, triad, tetrad
!                            TST_ROT, TST_MIR   , cell90, cell120, eq_sides ,
!                             GET_SYM, LENGTH, idum  , tmp
!
!     Utiliza las funciones:  TST_ROT, TST_MIR, GET_SYM, LENGTH  todas externas
!     Utiliza las subrutinas:

      LOGICAL, INTENT(OUT)     :: ok
      LOGICAL :: diad, triad, tetrad
      LOGICAL :: cell90, cell120, eq_sides
      INTEGER*4  idum
      REAL*8 tmp

! external functions

! external subroutine (Some compilers need them declared external)

! reinitialize random numbers in RAN3
      idum = -1
      diad   = .false.
      triad  = .false.
      tetrad = .false.
      cell90 =  ABS(cell_gamma - half*pi) < half*pi*eps6
      cell120 = ABS(cell_gamma - pi2/three) < pi2*eps6/three

! sample reciprocal space to get an idea of the sort of intensities
! that are out there.
! 360 degrees, no symmetry (-1)
! there is nothing to check. The center of symmetry is a given.
      IF(symgrpno == 1) THEN
        max_var = zero
        GO TO 900
      END IF

! 180 degrees, rotation only (2/M, 1st setting)
      IF(symgrpno == 2) THEN
        diad = tst_rot(2, idum, ok)
        IF(.NOT.ok) GO TO 999
        GO TO 900
      END IF

! 180 degrees, vertical mirror (2/M, 2nd setting)
      IF(symgrpno == 3) THEN
        h_mirror = tst_mir(1, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = max_var
        max_var = zero
        k_mirror = tst_mir(2, idum, ok)
        IF(.NOT.ok) GO TO 999
        max_var = MAX(tmp, max_var)
        GO TO 900
      END IF

! 90 degrees, vertical mirrors (MMM)
      IF(symgrpno == 4) THEN
        IF(.NOT.cell90) GO TO 910
        diad = tst_rot(2, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = max_var
        max_var = zero
        h_mirror = tst_mir(1, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = MAX(tmp, max_var)
        max_var = zero
        k_mirror = tst_mir(2, idum, ok)
        IF(.NOT.ok) GO TO 999
        max_var = MAX(tmp, max_var)
        GO TO 900
      END IF

! the following point groups require equi-sided cells
      eq_sides = ABS(cell_a - cell_b) <= half*eps6*(cell_a + cell_b)
      IF(.NOT.eq_sides) GO TO 920


! 120 degrees, rotation only (-3)
      IF(symgrpno == 5) THEN
        IF(.NOT.cell120) GO TO 910
        triad = tst_rot(3, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = max_var
        max_var = zero
        h_mirror = tst_mir(1, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = MAX(tmp, max_var)
        max_var = zero
        hk_mirror = tst_mir(3, idum, ok)
        IF(.NOT.ok) GO TO 999
        max_var = MAX(tmp, max_var)
        GO TO 900
      END IF

! 60 degrees, vertical mirrors (-3M)
      IF(symgrpno == 6) THEN
        IF(.NOT.cell120) GO TO 910
        triad = tst_rot(3, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = max_var
        max_var = zero
        h_mirror = tst_mir(1, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = MAX(tmp, max_var)
        max_var = zero
        hk_mirror = tst_mir(3, idum, ok)
        IF(.NOT.ok) GO TO 999
        max_var = MAX(tmp, max_var)
        GO TO 900
      END IF

! 90 degrees, rotation (4/M)
      IF(symgrpno == 7) THEN
        IF(.NOT.cell90) GO TO 910
        tetrad = tst_rot(4, idum, ok)
        IF(.NOT.ok) GO TO 999
        GO TO 900
      END IF

! 45 degrees, vertical mirrors (4/MMM)
      IF(symgrpno == 8) THEN
        IF(.NOT.cell90) GO TO 910
        tetrad = tst_rot(4, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = max_var
        max_var = zero
        h_mirror = tst_mir(1, idum, ok)
        IF(.NOT.ok) GO TO 999
        IF(.NOT.h_mirror) THEN
          tmp = MAX(tmp, max_var)
          max_var = zero
          hk_mirror = tst_mir(3, idum, ok)
          IF(.NOT.ok) GO TO 999
        END IF
        max_var = MAX(tmp, max_var)
        GO TO 900
      END IF

! 60 degrees, rotation (6/M)
      IF(symgrpno == 9) THEN
        IF(.NOT.cell120) GO TO 910
        diad = tst_rot(2, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = max_var
        max_var = zero
        triad = tst_rot(3, idum, ok)
        IF(.NOT.ok) GO TO 999
        max_var = MAX(tmp, max_var)
        GO TO 900
      END IF

! 30 degrees, vertical mirrors (6/MMM)
      IF(symgrpno == 10) THEN
        IF(.NOT.cell120) GO TO 910
        diad = tst_rot(2, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = max_var
        max_var = zero
        triad = tst_rot(3, idum, ok)
        IF(.NOT.ok) GO TO 999
        tmp = MAX(tmp, max_var)
        max_var = zero
        h_mirror = tst_mir(1, idum, ok)
        IF(.NOT.ok) GO TO 999
        IF(.NOT.h_mirror) THEN
          tmp = MAX(tmp, max_var)
          max_var = zero
          hk_mirror = tst_mir(3, idum, ok)
          IF(.NOT.ok) GO TO 999
        END IF
        max_var = MAX(tmp, max_var)
      END IF

      IF(symgrpno > 10) GO TO 500

      900 WRITE(op,*) '=> The diffraction data fits the point group symmetry ', pnt_grp(1:length(pnt_grp)),''''
      IF(max_var > eps6 .AND. max_var <= eps1) THEN
        WRITE(op,203) '   with a tolerance of one part in ', nint(one / max_var)
      ELSE IF(max_var > eps1) THEN
        WRITE(op,204) '   with a tolerance of one part in ', one / max_var
      ELSE
        WRITE(op,100) '   with a tolerance better than one part in a million.'
      END IF

      500 RETURN

! The user's guess is inconsistent with cell_gamma.
! Override the user.
      910 WRITE(op,200) 'The cell angle of',cell_gamma * rad2deg, ' degrees,'
      WRITE(op,202) ' is inconsistent with point group symmetry ''',  &
          pnt_grp(1:length(pnt_grp)),''''
      WRITE(op,300)
      check_sym = .false.
      symgrpno = get_sym(ok)
      IF(.NOT.ok) GO TO 999
      WRITE(op,205) pnt_grp(1:length(pnt_grp))
      IF(tolerance > eps6 .AND. tolerance <= eps1) THEN
        WRITE(op,203) '  with a tolerance of one part in ', nint(one / tolerance)
      ELSE IF(tolerance > eps1) THEN
        WRITE(op,204) '  with a tolerance of one part in ', one / tolerance
      ELSE
        WRITE(op,100) '  with a tolerance better than one part in a million.'
      END IF
      RETURN

! The user's guess is inconsistent with cell dimensions.
! Override the user.
      920 WRITE(op,201) 'The cell a and b dimensions, ',  &
          cell_a,' Angstroms by ',cell_b,' Angstroms,'
      WRITE(op,202) '   are inconsistent with point group symmetry ''',  &
          pnt_grp(1:length(pnt_grp)),''''
      WRITE(op,300)
! reset check_sym flag, since we are now evaluating from scratch
      check_sym = .false.
      max_var = zero
      symgrpno = get_sym(ok)
      IF(.NOT.ok) GO TO 999
      WRITE(op,205) pnt_grp(1:length(pnt_grp))

      IF(tolerance > eps6 .AND. tolerance <= eps1) THEN
        WRITE(op,203) '  with a tolerance of one part in ', nint(one / tolerance)
      ELSE IF(tolerance > eps1) THEN
        WRITE(op,204) '  with a tolerance of one part in ', one / tolerance
      ELSE
        WRITE(op,100) '  with a tolerance better than one part in a million.'
      END IF

      RETURN

      999 WRITE(op,100) 'ERROR in CHK_SYM'
      RETURN
      100 FORMAT(1X, a)
      200 FORMAT(1X, a, f7.3, a)
      201 FORMAT(1X, 2(a, f7.3), a)
      202 FORMAT(1X, 3A)
      203 FORMAT(1X, a, i6)
      204 FORMAT(1X, a, f4.1)
      205 FORMAT(1X, 'Diffraction point symmetry is found to be ''',a,'''')
      300 FORMAT(1X, 'Re-evaluating diffraction symmetry')
      END SUBROUTINE chk_sym


! ______________________________________________________________________
! Title: CHWDTH
! Author: MMJT
! Date: 6 Mar 1995; 31st Oct 1996
! Description:  This routine adds shape broadening caused by the finite
! lateral width of the layers. This routine does not add the broadening
! caused by finite width in the stacking direction. That is handled in
! a different manner by the routine INTEN2 and associated routines.
! The broadening handled here is added phenomenologically by applying a
! Lorentzian profile to the computed intensity at each peak in the h-k
! plane. If we are on the 00l axis the shape broadening is not
! symmetrical. For a Lorentzian shape broadening intensity profile in
! the h-k plane, the curvature of the Ewald sphere ensures a sharp
! onset of intensity with a slowly decaying tail at higher l values. For
! a symmetrical disk, the integrated intensity decays logarithmically.
! In the event that the crystal width is different in the h, and h-
! perpendicular directions, the disk of confusion becomes elongated
! into a streak perpendicular to l, and the tail becomes more
! Lorentzian-like. This is modelled in a phenomenological fashion, by
! mixing Lorentzian and logarithmic terms in a manner that depends
! on the ratio Wa/Wb. The costly logarithm function is avoided by
! using its derivative.
! When off the 00l axis, the broadening is modelled as a symmetric
! Lorentzian whose half-width depends on the angle that the Ewald
! sphere intercepts the disk of confusion (controlled by l). If
! the lateral dimensions are not equal, then the half width is also
! dependent on h and k. The Lorentzian is pre-computed in OPTIMZ to gain
! computational speed.
! This routine is called by GETSPC.

!      ARGUMENTS:
!            h       -  Reciprocal space index. (input).
!            k       -  Reciprocal space index. (input).
!            l0      -  Lower bound of the l reciprocal space index
!                       that is being integrated over. (input).
!            l1      -  Lower bound of the l reciprocal space index
!                       that is being integrated over. (input).
!            x       -  The integrated intensity value along the
!                       line defined by h,k,l0,l1. (input).
!            m       -  The current index of the array 'spec'
!                       corresponding to h,k,l0,l1. (input).
!          max_indx  -  Maximum array value in spec that will be
!                       accessed. (input).

!      COMMON VARIABLES:
!            uses:      spec, brd_spec, FFACT_SIZE, formfactor, d_theta
!                       ffact_scale, ffhkcnst, ffwdth

!        modifies:      spec, brd_spec
! ______________________________________________________________________

      SUBROUTINE chwdth(h,k,l0,l1,x,m,max_indx)
!     Utiliza las variables: DIFFaX.par , DIFFaX.inc,  h, k, m, max_indx ,n, p, i, indx
!                            S, h_wdth, n_hw, d_hk, norm, l, scale, avg, xx, dx, tmp
!     Utiliza las funciones:   S(h,k,l)
!     Utiliza las subrutinas:

      Integer,       Intent(In) :: h
      Integer,       Intent(In) :: k
      Real(Kind=8),  Intent(In) :: l0
      Real(Kind=8),  Intent(In) :: l1
      Real(Kind=8),  Intent(In) :: x
      Integer,       Intent(In) :: m
      Integer,       Intent(In) :: max_indx

      Integer      :: n, p, i, indx
      Real(Kind=8) :: s, h_wdth, n_hw, d_hk, norm, l, scale, avg, xx, dx, tmp

! indx indexes into the arrays spec and brd_spec
! n indexes into the array formfactor
! p is the index of the centroid of the array formfactor
! h_wdth contains the effective half-width of the size broadening
! ffhkcnst depends only on h and k, and was computed in GETSPC
! scale contains the calibration term for the formfactor array
! d_hk is the average radius in reciprocal Angstroms that the
! Ewald sphere of radius theta+d_theta intercepts the 00l plane
! brd_spc is used for temporary storage.
! We are only accessing half of the symmetrical formfactor array

! statement functions
! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0

!***********************************************************************
! special case if we are on the l axis
! Don't bother if the intensity is negligible in the first place
      IF(h == 0 .AND. k == 0 .AND. x > ten*tiny_inty) THEN
        l = half*(l0+l1)
        d_hk = two*l*d_theta / (lambda*cell_c*ffwdth*ffwdth)
        norm = zero
        indx = 0
        xx = wa / wb
        IF(xx > one) xx = one / xx
        10   indx = indx + 1
        tmp = one / (one + indx*d_hk)
! balance streak, versus disk of confusion (phenomenological treatment)
! xx=1 means disk, xx=0 means streak. Intensity falls off more
! rapidly for a streak
        dx = ((one-xx)*SQRT(DBLE(indx))*tmp + xx)*tmp
        IF(m+indx-1 <= max_indx) brd_spc(m+indx-1) = dx
        norm = norm + dx
! eps5 is reasonable. However, it may be worth experimenting more.
        IF(dx < eps5) GO TO 20
        GO TO 10
        20   CONTINUE
        norm = x / norm
        DO  i = 0, indx - 1
          IF(m+i <= max_indx) spec(m+i) = spec(m+i) + norm*brd_spc(m+i)
        END DO
! We were on the 00l axis, we can exit now
        RETURN
      END IF

!***********************************************************************
! We are not on the l-axis. Broadening is handled differently.
! scale relates the formfactor array indices to the spec array indices.
! h_wdth and ffact_scale should never be zero. In case they are,
! make scale large enough so that n.ge.p-1 and the loop below is
! exited early
      h_wdth = ffhkcnst / SQRT(s(h,k,half*(l0+l1)))
      IF(h_wdth > zero .AND. ffact_scale > zero) THEN
        scale = d_theta / (ffact_scale * h_wdth)
      ELSE
        scale = ffact_size
      END IF
      p = ffact_size/2 + 1
      norm = one
      brd_spc(m) = one
      indx = 0
      40 indx = indx + 1
      n_hw = indx*scale
      n = INT(n_hw)
      IF(n >= p-1) GO TO 50
! linear interpolation of the pre-computed pseudo-Lorentzian
      xx = n_hw - n
      avg = (one-xx)*formfactor(p+n) + xx*formfactor(p+n+1)
      IF(m+indx <= max_indx) brd_spc(m+indx) = avg
      IF(m-indx > 0)        brd_spc(m-indx) = avg
! intensity x is being redistributed. We will need to normalize later
      norm = norm + two*avg
      GO TO 40
      50 CONTINUE
      norm = x / norm
      spec(m) = spec(m) + norm*brd_spc(m)
      DO  i = 1, indx - 1
        IF(m+i <= max_indx) spec(m+i) = spec(m+i) + norm*brd_spc(m+i)
        IF(m-i > 0)        spec(m-i) = spec(m-i) + norm*brd_spc(m-i)
      END DO
      RETURN
      END SUBROUTINE chwdth

!***********************************************************************
!***************************LINPACK ROUTINES****************************
!***********************************************************************

! The following are the standard Linpack routines for solving complex
! simultaneous equations. They were found to reduce DIFFaX run time by
! significant amount (30% in one case) compared with the Numerical
! Recipes routines LUDCMP and LUBKSB. The only changes are

!                         complex -> complex*16
!                         real    -> real*8
!                         real()  -> dble()
!                         aimag   -> dimag

!***********************************************************************
! ______________________________________________________________________
! Title: CGESL (LINPACK ROUTINE)
! Author: cleve moler, university of new mexico, argonne national lab.
! Date: linpack. this version dated 08/14/78
! Description:
!     CGESL solves the complex system
!     a * x = b  or  ctrans(a) * x = b
!     using the factors computed by cgeco or CGEFA.

!     on entry

!        a       complex(lda, n)
!                the output from cgeco or CGEFA.

!        lda     integer
!                the leading dimension of the array  a .

!        n       integer
!                the order of the matrix  a .

!        ipvt    integer(n)
!                the pivot vector from cgeco or CGEFA.

!        b       complex(n)
!                the right hand side vector.

!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  ctrans(a)*x = b  where
!                            ctrans(a)  is the conjugate transpose.

!     on return

!        b       the solution vector  x .

!     error condition

!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if cgeco has set rcond .gt. 0.0
!        or CGEFA has set info .eq. 0 .

!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call cgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call CGESL(a,lda,n,ipvt,c(1,j),0)
!        10 continue

!     subroutines and functions

!     blas CAXPY,CDOTC
!     fortran conjg

! ______________________________________________________________________

      SUBROUTINE cgesl(a,lda,n,ipvt,b,job)
!     Utiliza las variables:lda,n,ipvt(:),job  , a(:,:),b(:) , CDOTC,t
!                           k,kb,l,nm1
!     Utiliza las funciones:  CDOTC externa,
!     Utiliza las subrutinas: CAXPY

      IMPLICIT NONE
      COMPLEX*16, INTENT(IN OUT)               :: a(:,:) !a(lda,1)
      INTEGER                                  :: lda
      INTEGER, INTENT(IN)                      :: n
      INTEGER, INTENT(IN)                      :: ipvt(:)
      COMPLEX*16, INTENT(IN OUT)               :: b(:)
      INTEGER                                  :: job
      COMPLEX*16 t
      INTEGER :: k,kb,l,nm1
! MMJT: external subroutine
! external function
      nm1 = n - 1
      IF (job /= 0) GO TO 50

!        job = 0 , solve  a * x = b
!        first solve  l*y = b

      IF (nm1 >= 1) THEN
        DO  k = 1, nm1
          l = ipvt(k)
          t = b(l)
          IF (l /= k) then
            b(l) = b(k)
            b(k) = t
          end if
          CALL caxpy(n-k,t,a(k+1:,k),1,b(k+1:),1)
        END DO
      END IF

!        now solve  u*x = y

      DO  kb = 1, n
        k = n + 1 - kb
        b(k) = b(k)/a(k,k)
        t = -b(k)
        CALL caxpy(k-1,t,a(1:,k),1,b(1:),1)
      END DO
      GO TO 100
      50 CONTINUE

!        job = nonzero, solve  ctrans(a) * x = b
!        first solve  ctrans(u)*y = b

      DO  k = 1, n
        t = cdotc(k-1,a(1:,k),1,b(1:),1)
        b(k) = (b(k) - t)/CONJG(a(k,k))
      END DO

!        now solve ctrans(l)*x = y

      IF (nm1 >= 1) THEN
        DO  kb = 1, nm1
          k = n - kb
          b(k) = b(k) + cdotc(n-k,a(k+1:,k),1,b(k+1:),1)
          l = ipvt(k)
          IF (l /= k) then
              t  = b(l)
            b(l) = b(k)
            b(k) = t
          END IF
        END DO
      END IF

      100 CONTINUE
      RETURN
      END SUBROUTINE cgesl
! ______________________________________________________________________
! Title: CAXPY (LINPACK ROUTINE)
! Author: jack dongarra
! Date: linpack, 3/11/78
! Description: constant times a vector plus a vector.
! ______________________________________________________________________

      SUBROUTINE caxpy(n,ca,cx,incx,cy,incy)

!     Utiliza las variables: cx(:),cy(:),ca , n,incx,incy , i,ix,iy
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      IMPLICIT NONE
      INTEGER, INTENT(IN)                      :: n
      COMPLEX*16, INTENT(IN)                   :: ca
      COMPLEX*16, INTENT(IN)                   :: cx(:)
      INTEGER, INTENT(IN)                      :: incx
      COMPLEX*16, INTENT(IN OUT)               :: cy(:)
      INTEGER, INTENT(IN)                      :: incy

      INTEGER :: i,ix,iy

      IF(n <= 0)RETURN
      IF (ABS(DBLE(ca)) + ABS(DIMAG(ca)) == 0.0 ) RETURN
      IF(incx == 1 .AND. incy == 1) then
!        code for both increments equal to 1
         DO  i = 1,n
           cy(i) = cy(i) + ca*cx(i)
         END DO
      ELSE
!        code for unequal increments or equal increments
!          not equal to 1
         ix = 1
         iy = 1
         IF(incx < 0)ix = (-n+1)*incx + 1
         IF(incy < 0)iy = (-n+1)*incy + 1
         DO  i = 1,n
           cy(iy) = cy(iy) + ca*cx(ix)
           ix = ix + incx
           iy = iy + incy
         END DO
      END IF

      RETURN
      END SUBROUTINE caxpy
! ______________________________________________________________________
! Title: CDOTC (LINPACK ROUTINE)
! Author: jack dongarra
! Date: linpack,  3/11/78.
! Description:
!     forms the dot product of two vectors, conjugating the first
!     vector.
! ______________________________________________________________________

      COMPLEX*16 FUNCTION cdotc(n,cx,incx,cy,incy)
      IMPLICIT NONE

!     Utiliza las variables: cx(:),cy(:), n,incx,incy , i,ix,iy,ctemp
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      INTEGER, INTENT(IN)                      :: n
      COMPLEX*16, INTENT(IN OUT)               :: cx(:)
      INTEGER, INTENT(IN)                      :: incx
      COMPLEX*16, INTENT(IN)                   :: cy(:)
      INTEGER, INTENT(IN)                      :: incy
      COMPLEX*16 ctemp
      INTEGER :: i,ix,iy

      ctemp = (0.0,0.0)
      cdotc = (0.0,0.0)
      IF(n <= 0) RETURN
      IF(incx == 1 .AND. incy == 1)GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      ix = 1
      iy = 1
      IF(incx < 0)ix = (-n+1)*incx + 1
      IF(incy < 0)iy = (-n+1)*incy + 1
      DO  i = 1,n
        ctemp = ctemp + CONJG(cx(ix))*cy(iy)
        ix = ix + incx
        iy = iy + incy
      END DO
      cdotc = ctemp
      RETURN

!        code for both increments equal to 1

      20 DO  i = 1,n
        ctemp = ctemp + CONJG(cx(i))*cy(i)
      END DO
      cdotc = ctemp
      RETURN
      END FUNCTION cdotc
! ______________________________________________________________________
! Title: CGEFA (LINPACK ROUTINE)
! Author: cleve moler, university of new mexico, argonne national lab.
! Date: linpack. this version dated 08/14/78
! Description:

!     CGEFA factors a complex matrix by gaussian elimination.

!     CGEFA is usually called by cgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for cgeco) = (1 + 9/n)*(time for CGEFA) .

!     on entry

!        a       complex(lda, n)
!                the matrix to be factored.

!        lda     integer
!                the leading dimension of the array  a .

!        n       integer
!                the order of the matrix  a .

!     on return

!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.

!        ipvt    integer(n)
!                an integer vector of pivot indices.

!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that CGESL or cgedi will divide by zero
!                     if called.  use  rcond  in cgeco for a reliable
!                     indication of singularity.

!     subroutines and functions

!     blas CAXPY,CSCAL,ICAMAX
!     fortran abs,aimag,real
! ______________________________________________________________________

      SUBROUTINE cgefa(a,lda,n,ipvt,info)

!     Utiliza las variables:  lda,n,ipvt(:),info  , a(:,:) , t, ICAMAX,j,k,kp1,l,nm1
!                             zdum ,  cabs1,
!     Utiliza las funciones:  ICAMAX  externa,  cabs1(zdum)
!     Utiliza las subrutinas: CSCAL, CAXPY
      IMPLICIT NONE
      COMPLEX*16, INTENT(IN OUT)               :: a(lda,n)
      INTEGER, INTENT(IN)                      :: lda
      INTEGER, INTENT(IN OUT)                  :: n
      INTEGER, INTENT(IN OUT)                  :: ipvt(:) !ipvt(1)
      INTEGER, INTENT(IN OUT)                  :: info
      COMPLEX*16 t
      INTEGER :: j,k,kp1,l,nm1

      COMPLEX*16 zdum
      REAL*8 cabs1

! MMJT: external subroutine
! statement function
      cabs1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))

!     gaussian elimination with partial pivoting

      info = 0
      nm1 = n - 1
      IF (nm1 < 1) GO TO 70
      DO  k = 1, nm1
        kp1 = k + 1

!        find l = pivot index

        l = icamax(n-k+1,a(k,k),1) + k - 1
        ipvt(k) = l

!        zero pivot implies this column already triangularized

        IF (cabs1(a(l,k)) == 0.0E0) GO TO 40

!           interchange if necessary

        IF (l == k) GO TO 10
        t = a(l,k)
        a(l,k) = a(k,k)
        a(k,k) = t
        10       CONTINUE

!           compute multipliers

        t = -(1.0E0,0.0E0)/a(k,k)
        CALL cscal(n-k,t,a(k+1,k),1)

!           row elimination with column indexing

        DO  j = kp1, n
          t = a(l,j)
          IF (l == k) GO TO 20
          a(l,j) = a(k,j)
          a(k,j) = t
          20          CONTINUE
          CALL caxpy(n-k,t,a(k+1:,k),1,a(k+1:,j),1)
        END DO
        GO TO 50
        40    CONTINUE
        info = k
        50    CONTINUE
      END DO
      70 CONTINUE
      ipvt(n) = n
      IF (cabs1(a(n,n)) == 0.0E0) info = n
      RETURN
      END SUBROUTINE cgefa
! ______________________________________________________________________
! Title: CSCAL (LINPACK ROUTINE)
! Author: jack dongarra
! Date: linpack,  3/11/78.
! Description: scales a vector by a constant.
! ______________________________________________________________________

      SUBROUTINE  cscal(n,ca,cx,incx)
!     Utiliza las variables: ca,cx(1), incx,n , i,nincx
!     Utiliza las funciones:
!     Utiliza las subrutinas: ,
      IMPLICIT NONE
      INTEGER                                  :: n
      COMPLEX*16                               :: ca
      COMPLEX*16, INTENT(IN OUT)               :: cx(1)
      INTEGER                                  :: incx
      INTEGER :: i,nincx

      IF(n <= 0)RETURN
      IF(incx == 1)GO TO 20
!        code for increment not equal to 1
      nincx = n*incx
      DO  i = 1,nincx,incx
        cx(i) = ca*cx(i)
      END DO
      RETURN
!        code for increment equal to 1
      20 DO  i = 1,n
        cx(i) = ca*cx(i)
      END DO
      RETURN
      END SUBROUTINE  cscal
! ______________________________________________________________________
! Title: ICAMAX (LINPACK ROUTINE)
! Author: jack dongarra
! Date: linpack, 3/11/78
! Description:
!     finds the index of element having max. absolute value.
! ______________________________________________________________________

      INTEGER FUNCTION icamax(n,cx,incx)
!     Utiliza las variables: cx(1) , incx,n  ,smax  , i,ix ,  zdum, cabs1
!     Utiliza las funciones: cabs1(zdum)
!     Utiliza las subrutinas:
      IMPLICIT NONE
      INTEGER                                  :: n
      COMPLEX*16, INTENT(IN)                   :: cx(1)
      INTEGER                                  :: incx
      REAL*8 smax
      INTEGER :: i,ix
      COMPLEX*16 zdum

! statement function
      REAL*8 cabs1
      cabs1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))

      icamax = 0
      IF( n < 1 ) RETURN
      icamax = 1
      IF(n == 1)RETURN
      IF(incx == 1) GO TO 20

!        code for increment not equal to 1

      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      DO  i = 2,n
        IF(cabs1(cx(ix)) > smax) then
          icamax = i
          smax = cabs1(cx(ix))
        end if
        ix = ix + incx
      END DO
      RETURN
!        code for increment equal to 1
      20 smax = cabs1(cx(1))
      DO  i = 2,n
        IF(cabs1(cx(i)) <= smax) CYCLE
        icamax = i
        smax = cabs1(cx(i))
      END DO
      RETURN
      END FUNCTION icamax

!!S Dfiles
! ______________________________________________________________________
! Title: DETUN
! Author: MMJT
! Date: 27 Jan 1989
! Description: This subroutine detunes the sharp peaks so that
! they can be integrated. In effect, it modifies the stacking
! probability matrix such that the determinant will never be
! exactly zero.
!      ARGUMENTS:
!            No input arguments.
!      COMMON VARIABLES:
!            uses:    n_layers
!        modifies:    detune
! ______________________________________________________________________

      SUBROUTINE detun()
!     Utiliza las variables: DIFFaX.par , DIFFaX.inc, i, j,delta ,
!     Utiliza las funciones:
!     Utiliza las subrutinas:
      INTEGER*4 i, j
      REAL*8 delta

      delta = eps3
! A value of delta = 0.001 is found to be optimum.
! If preferred, user can specify 'delta' interactively
! using the following three lines.
!   30 write(op,400) 'Enter detune parameter'
!      read(cntrl,*,err=30,end=999) delta
!      if(CFile) write(op,401) delta
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          detune(j,i) = one - ABS(delta)
        END DO
      END DO

      RETURN
!  999 stop 'ERROR: Bad delta value. DIFFaX aborted.'
!  400 format(1x, a)
!  401 format(1x, g12.5)
      END SUBROUTINE detun

! ______________________________________________________________________
! Title: DUMP
! Author: MWD and MMJT
! Date: 18 Aug 1988; 15 Mar 1995
! Description: This subroutine prints out the data read in from
! the data file. It reassures the user that the data was read in
! correctly.

!      ARGUMENTS:
!            infile  -  The name of the input data file. (input).
!            ok      -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:      CFile, GAUSS, LORENZ, NEUTRN, PI2, PS_VGT,
!                       RAD2DEG, SymGrpNo, X_RAY, a_B, a_name,
!                       a_number, a_occup, a_pos, a_type, blurring,
!                       cell_a, cell_b, cell_c, cell_gamma, cntrl
!                       d_theta, e_sf, FWHM, inf_thick, l_actual
!                       l_alpha, l_cnt, l_g, l_n_atoms, l_r, l_symmetry
!                       lambda, l_seq, n_actual, n_atoms, n_layers
!                       pnt_grp, pv_gamma, pv_u, pv_v, pv_w, r_B11
!                       r_B12, r_B22, r_B23, r_B31, r_B33, rad_type
!                       recrsv, rndm, th2_max, th2_min
!                       tolerance, xplcit

!        modifies:      no COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE dump(infile,dmp,ok)

!     Utiliza las variables: DIFFaX.par , DIFFaX.inc, infile, ok, scale, atom_cnt(MAX_TA),
!                            cum_atom_cnt(MAX_TA), norm  , i, i2, j, n, type, num_types, tot_types
!                            LENGTH, print_width ,dmpfile , list,
!     Utiliza las funciones: LENGTH  externa,
!     Utiliza las subrutinas:GETFNM


      CHARACTER (LEN=*), INTENT(IN) :: infile
      integer, intent(in)           :: dmp
      LOGICAL, INTENT(OUT)          :: ok

      REAL(kind=8) :: scale, atom_cnt(max_ta), cum_atom_cnt(max_ta), norm
      INTEGER :: i, i2, j, n, TYPE, num_types, tot_types
      INTEGER :: print_width
      CHARACTER(LEN=80) :: list(5)


      i2 = 1
      !WRITE(op,200) 'Enter 1 for full atomic position dump: '
      !READ(cntrl,*,ERR=99,END=9999) i2
      !IF(cfile) WRITE(op,'(1x,i3)') i2
      WRITE(dmp,"(/a)") ' ---------------------------------------------------------------------'
      WRITE(dmp,"(a)")  ' Writing input data parameter read in file '//trim(infile)
      WRITE(dmp,"(a/)") ' ---------------------------------------------------------------------'

!title
      WRITE(dmp,"(a)") 'TITLE:', ttl
! sundry details about layers
      WRITE(dmp,125) 'Number of layers = ', n_layers
      WRITE(dmp,125) 'Number of unique layers = ', n_actual
      WRITE(dmp,125) 'Number of different atom types = ', n_atoms
! cell dimensions
      WRITE(dmp,"(a,3f12.6,f12.2)") ' cell_a,      cell_b,      cell_c,      cell_gamma ',  &
          cell_a, cell_b, cell_c, rad2deg * cell_gamma


! in-plane layer widths
      if(finite_width) then
        if(Wa < inf_width) then
          write(dmp,126) 'width along a', Wa
            if (Wb < inf_width) then
              write(dmp,126) 'width along b', Wb
            end if
        end if
      else
          write(*,"(a)") ' => Layers are to be treated as having infinite lateral width.'
          write(dmp,200) 'Layers are to be treated as having infinite lateral width.'
      end if


! radiation type
      list(1) = 'X-RAY '
      list(2) = 'NEUTRON '
      list(3) = 'ELECTRON '
      WRITE(dmp,100) 'radiation type = ', list(rad_type+1)
! wavelength
      WRITE(dmp,170) 'Wavelength lambda = ', lambda, lambda2, ratio
! symmetry
      WRITE(dmp,100) 'Diffraction point group symmetry specified = ',pnt_grp
      IF(symgrpno == UNKNOWN) WRITE(dmp,175)  &
          'Tolerance on intensities in symmetry evaluation = +/-',  &
          tolerance * hundred, ' %'
! instrumental broadening to simulate
      list(1) = 'NONE'
      list(2) = 'GAUSSIAN'
      list(3) = 'LORENTZIAN'
      list(4) = 'PSEUDO-VOIGT'
      i = blurring + 1
! see if it's a pure Gaussian pseudo-Voigt
      IF(blurring == pv_gss) i = 2
! see if it's a pure Lorentzian pseudo-Voigt
      IF(blurring == pv_lrn) i = 3
      WRITE(dmp,100) 'Instrumental broadening to be simulated = ', list(i)
      IF(blurring == gauss .OR. blurring == lorenz) WRITE(dmp,170)  &
          'Full width half-maximum = ', fwhm
      IF(blurring == ps_vgt) THEN
        WRITE(dmp,200) 'Pseudo-Voigt parameters:   u,       v,       w,       x,      Dg,          Dl'
        WRITE(dmp,180) pv_u, pv_v, pv_w, pv_x, pv_dg, pv_dl
      END IF
      IF(blurring == pv_gss) THEN
        WRITE(dmp,200) 'Gaussian parameters:      u,       v,       w'
        WRITE(dmp,185) pv_u, pv_v, pv_w
      END IF
      IF(blurring == pv_lrn) THEN
        WRITE(dmp,200) 'Lorentzian parameters:    u,       v,       w'
        WRITE(dmp,185) pv_u, pv_v, pv_w
      END IF
! are we to trim out the origin?
      IF(blurring /= NONE) THEN
        IF(trim_origin) THEN
          WRITE(dmp,200) '  Intensity near origin will be ignored'
        ELSE
          WRITE(dmp,200) '  Intensity near origin will be included'
        END IF
      END IF

! equivalent layers
      WRITE(dmp,200) ' '
      DO  i = 1, n_layers
        WRITE(dmp,110) 'LAYER ',i, ' is equivalent to fundamental LAYER ',l_actual(i)
        WRITE(dmp,170) '   Existence probability = ', l_g(i)
      END DO

      DO  j = 1, max_ta
        cum_atom_cnt(j) = zero
      END DO

! layer structure
      WRITE(dmp,200) ' '
      DO  i = 1, n_actual
        WRITE(dmp,130) 'fundamental LAYER ',i
        list(1) = 'NONE '
        list(2) = 'CENTROSYMMETRIC '
        WRITE(dmp,100) 'symmetry = ', list(l_symmetry(i)+1)
! write out the layer composition
        WRITE(dmp,200) 'LAYER composition is:'
        DO  j = 1, max_ta
          atom_cnt(j) = zero
        END DO
        num_types = zero
        scale = one
        IF(l_symmetry(i) == centro) scale = two
        DO  j = 1, l_n_atoms(i)
          TYPE = a_type(j,i)
          IF(TYPE > num_types) num_types = TYPE
          atom_cnt(TYPE) = atom_cnt(TYPE) + a_occup(j,i) * scale
        END DO
       ! IF(num_types > tot_types) tot_types = num_types
! accumulate weighted atom count
        DO  j = 1, max_ta
          cum_atom_cnt(j) = cum_atom_cnt(j) + atom_cnt(j) * l_g(i)
        END DO
        DO  j = 1, num_types
          WRITE(dmp,280) atom_l(j), ':', atom_cnt(j)
        END DO
        WRITE(dmp,200) ' '

! do we want all the details about each atom?
        IF(i2 /= 1) CYCLE
! yes, go for it.
        DO  j = 1, l_n_atoms(i)
          WRITE(dmp,200) ' '
          WRITE(dmp,120) 'ATOM number', a_number(j,i), ' in layer',i,  &
              ' is ''', a_name(j,i),''''
          WRITE(dmp,121)
          WRITE(dmp,122) a_pos(1,j,i), a_pos(2,j,i), a_pos(3,j,i),  &
              a_b(j,i), a_occup(j,i)
          IF(rad_type == x_ray) THEN
            WRITE(dmp,140) 'X-ray scattering factor data Ai, Bi',  &
                (x_sf(n,a_type(j,i)),n=1,9)
          ELSE IF(rad_type == neutrn) THEN
            WRITE(dmp,140) 'Neutron scattering factor data', n_sf(a_type(j,i))
          ELSE IF(rad_type == electn) THEN
            WRITE(dmp,140)'Electron scattering factor data Ai, Bi and Z',  &
                (x_sf(n,a_type(j,i)),n=1,9), e_sf(a_type(j,i))
          ELSE
            WRITE(op,200) 'ERROR: Illegal radiation type in DUMP'
          END IF
        END DO
        WRITE(dmp,200) ' '
      END DO

! Print out average composition
      IF(.NOT.xplcit) THEN
        norm = zero
        DO  j = 1, num_types
          norm = norm + cum_atom_cnt(j)
        END DO
        WRITE(dmp,200) ' '
        WRITE(dmp,200) 'Average crystal composition is:'
        DO  j = 1, num_types
          WRITE(dmp,281) atom_l(j), ':', cum_atom_cnt(j) / norm
        END DO
        WRITE(dmp,200) ' '
      END IF

! stacking details
      WRITE(dmp,200) ' '
      WRITE(dmp,200) 'STACKING'
      IF(recrsv) THEN
        WRITE(dmp,210) 'RECURSIVELY'
        IF(inf_thick) THEN
          WRITE(dmp,220) 'INFINITE'
        ELSE
          WRITE(dmp,230) nint(l_cnt)
        END IF
      ELSE IF(xplcit) THEN
        IF(randm) THEN
          WRITE(dmp,240)
        ELSEIF(semirandm) THEN
          WRITE(dmp,245)
        ELSE
          WRITE(dmp,250)
        END IF
        WRITE(dmp,260) 'Sequence for ', nint(l_cnt), ' layers is:'
        print_width = 25
        j = 1
        n = print_width
        35   IF(n > nint(l_cnt)) n = nint(l_cnt)
        WRITE(dmp,270) (l_seq(i), i = j, n)
        IF(n < nint(l_cnt)) THEN
          j = j + print_width
          n = n + print_width
          GO TO 35
        END IF
      END IF

! stacking transition data
      WRITE(dmp,200) ' '
      WRITE(dmp,200) ' '
      WRITE(dmp,200) 'TRANSITIONS'
      WRITE(dmp,200) 'Layer stacking probabilities and stacking vectors.'
      WRITE(dmp,200) '  i  j |   alpha_ij      Rx_ij        Ry_ij        Rz_ij'
      WRITE(dmp,200) '  -----|-------------------------------------------------'
      DO  i = 1, n_layers
        WRITE(dmp,200) '       |'
        DO  j = 1, n_layers
          WRITE(dmp,150) i, j, l_alpha(j,i), l_r(1,j,i), l_r(2,j,i), l_r(3,j,i)
        END DO
      END DO

      WRITE(dmp,200) ' '
      WRITE(dmp,200) 'Anisotropic Layer stacking ''UNCERTAINTY'' FACTORS.'
! MMJT: 3/15/95. Ordering of B23 and B31 swapped. B31 -> C13
      WRITE(dmp,200) '  i  j |     C11      C22      C33      C12      C13      C23'
      WRITE(dmp,200)  &
          '  -----|-------------------------------------------------------'
      DO  i = 1, n_layers
        WRITE(dmp,200) '       |'
        DO  j = 1, n_layers
! MMJT: 3/15/95. Ordering of B23 and B31 swapped
          WRITE(dmp,151) i, j, r_b11(j,i), r_b22(j,i), r_b33(j,i),  &
              r_b12(j,i), r_b31(j,i), r_b23(j,i)
        END DO
      END DO

      WRITE(dmp,200) ' '   !Per deixar espai
      WRITE(dmp,200) ' '

      write(dmp,"(a)")     " CALCULATION  "
        if (opt == 0) then
          write(dmp,"(a)")          " SIMULATION"
          write(dmp,"(3f10.4)")     th2_min, th2_max, d_theta
        elseif (opt == 3) then
          write(dmp,"(2a)")          " LOCAL_OPTIMIZER   ", opti%method
          write(dmp,"(a,i4)")          " MXFUN  ", opti%mxfun
          write(dmp,"(a,f10.4)")          " EPS  ", opti%eps
          write(dmp,"(a, i2)")          " IOUT  ", opti%iout
          write(dmp,"(a,f10.7)")          " ACC  ", opti%acc
        elseif (opt == 4) then
          write(dmp,"(a)")          " LMQ"
          if (Cond%constr) write(dmp,"(a,f5.2)")          " BOXP    " , Cond%percent
          write(dmp,"(a,i4)")    " CORRMAX    ", cond%corrmax
          write(dmp,"(a,i4)")    " MAXFUN     ", cond%icyc
          write(dmp,"(a,f10.4)")    " TOL     ", cond%tol
          write(dmp,"(a,i2)")    " Nprint     ", cond%nprint
        else
          write(*,*) "ERROR writing *.dmp file: Problem with calculation section"
        end if

       WRITE(dmp,200) ' '   !Per deixar espai
       WRITE(dmp,200) ' '

       if(opt == 3 .or. opt == 4) then
         write(dmp,"(a)")              "  "
         write(dmp,"(a)")          " EXPERIMENTAL"
         write(dmp,"(2a)")         " FILE  ", dfile
         if (nexcrg /= 0) then
           write(dmp,"(a, i2)")    " EXCLUDED_REGIONS  ",  nexcrg
           do i=1,nexcrg
             write(dmp,"(2f10.4)")  alow(i),ahigh(i)
           end do

         end if
         write(dmp,"(2a)")         " FFORMAT  ",fmode
         write(dmp,"(2a)")         " BGR  ",background_file
         write(dmp,"(2a)")         " BCALC  ",mode
       end if

      999 RETURN
      9999 ok = .false.
      WRITE(op,200) 'DUMP aborted'
      RETURN
      100 FORMAT(1X, 2A)
      110 FORMAT(1X, a, i4, a, i4)
      120 FORMAT(1X, a, i3, a, i2, 3A)
      121 FORMAT(4X,'x_rel',8X,'y_rel',8X,'z_rel',10X,'DW',10X,'Occ')
      122 FORMAT(1X, 5G13.5)
      125 FORMAT(1X, a40, i4)
      126 FORMAT(1X, 'layer characteristic ', a, ' = ', f9.2, ' Angstroms')
      127 FORMAT(1X, 'layer characteristic ', a, ' = INFINITY Angstroms')
      128 FORMAT(1X, '   which is equivalent to ', f9.2, ' unit cells')
      130 FORMAT(1X, a, i4)
      140 FORMAT(1X, a, / 5G13.5, / 4G13.5, i4)
      150 FORMAT(1X, 2I3, ' |', 4(1X, g12.5))
      151 FORMAT(1X, 2I3, ' |', 6(1X, f8.3))
      170 FORMAT(1X, a, 3(g12.5))
      175 FORMAT(1X, a, g12.5, a)
      180 FORMAT(21X, 4(2X, f7.3), 2(4x, f9.3))
      185 FORMAT(21X, 3(2X, f7.3))
      200 FORMAT(1X, a)
      210 FORMAT(1X, 'Stacking is to be treated ', a, ' by FAULTS.')
      220 FORMAT(1X, 'Number of layers along the fault direction is ', a)
      230 FORMAT(1X, 'There are ',i5 ,' layers along the fault direction')
      240 FORMAT(1X, 'Sequencing is defined RANDOMLY by FAULTS.')
      245 FORMAT(1X, 'Sequencing is defined SEMIRANDOMLY by FAULTS.')
      250 FORMAT(1X, 'Sequencing is defined EXPLICITLY by the user.')
      260 FORMAT(1X, a, i4, a)
      270 FORMAT(1X, 30I3)
      280 FORMAT(23X, 2A, f7.2)
      281 FORMAT(23X, 2A, f7.4)
      300 FORMAT(1X, 3A)
      END SUBROUTINE dump

!!S Efiles
! ______________________________________________________________________
! Title: EQUALB
! Author: MMJT
! Date: 13 Mar 1990; 21 July 1997
! Description:  This routine determines if all of the stacking
! uncertainty parameters are identical. There are six arrays to be
! tested, namely r_B11, r_B22, r_B33, r_B12, r_B23 and r_B31. These are
! passed one at a time from OPTIMZ as r_B. The average value of the
! r_B parameters is returned in a_B.

!      ARGUMENTS:
!            r_B  -  Array of stacking uncertainty parameters. (input).
!            av_B  -  Average of r_B. (output).

!      COMMON VARIABLES:
!            uses the array 'there'. n_layers
! ______________________________________________________________________

      LOGICAL FUNCTION equalb(r_b, av_b)
!     Utiliza las variables: DIFFaX.par , DIFFaX.inc,r_B(MAX_L,MAX_L), av_B, i, j, m , error
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      REAL*8, INTENT(IN OUT)                   :: r_b(max_l,max_l)
      REAL*8, INTENT(IN OUT)                   :: av_b
      INTEGER*4 i, j, m
      REAL*8 error

      av_b = zero
      m = 0
      DO  i = 1, n_layers
        DO  j = 1, n_layers
! Examine only those transitions that actually occur
          IF(there(j,i)) THEN
            m = m + 1
            av_b = av_b + ABS(r_b(j,i))
          END IF
        END DO
      END DO

! Take average
      IF(m /= 0) av_b = av_b / DBLE(m)

      error = zero
! find absolute deviation of coefficients
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          IF(there(j,i)) error = error + ABS(r_b(j,i) - av_b)
        END DO
      END DO
! get relative error
      IF(av_b /= zero) error = error / av_b

      equalb = ABS(error) <= eps3

      RETURN
      END FUNCTION equalb


!!S GETfiles1
! ______________________________________________________________________
! ______________________________________________________________________
! Title: GETFNM
! Author: JRC
! Date: 08 Sept 2010
! Description: This subroutine makes sure that we do not write over
! an existing file, and creates a suitable filename (name2). name2
! equals 'name1.append' append is supplied by the calling routine.
! If 'name1.append' is the name of an existing file, then append is
! set to 'append#' where #=1,2.
!
!      ARGUMENTS:
!            name1   -  The name of the input data file. (input).
!            name2   -  A derivative filename. (output).
!            append  -  A token to be appended to name2. (input).
!            ok      -  logical flag signalling all went well.
!                                                      (output).
! ______________________________________________________________________

      SUBROUTINE getfnm(name1,name2,APPEND,ok)
        CHARACTER (LEN=*), INTENT(IN)        :: name1
        CHARACTER (LEN=*), INTENT(OUT)       :: name2
        CHARACTER (LEN=*), INTENT(IN)        :: APPEND
        LOGICAL, INTENT(IN OUT)              :: ok

        LOGICAL :: fexist
        CHARACTER (LEN=len(name1)) :: bname
        integer :: ln, i

        ok=.false.

        i=index(name1,".",back=.true.)
        if(i == 0) then
           bname=name1
        else
           bname=name1(1:i-1)
        end if
        !ln=len_trim(bname) + len_trim(APPEND) + 3
        !if(ln > length(name2) ) return

        i=1
        WRITE(name2,"(a,i1)") trim(bname),i
        name2=trim(name2)//trim(append)


        do
          INQUIRE(FILE = name2, EXIST = fexist)
          if(fexist) then
            i=i+1
            name2=" "
            Select Case (i)
              Case(1:9)
                  WRITE(name2,"(a,i1)") trim(bname),i
              Case(10:99)
                  WRITE(name2,"(a,i2)") trim(bname),i
              Case(100:999)
                  WRITE(name2,"(a,i3)") trim(bname),i
            End Select
            name2=trim(name2)//trim(append)
          else
            exit
          end if
        end do


        ok=.true.
        return
      END SUBROUTINE getfnm
! ______________________________________________________________________
! Title: GETLAY
! Author: MMJT
! Date: 4 Oct 1989
! Description: This function generates a random sequence of layers
! weighted by the stacking probabilities. Needed only when the
! 'EXPLICIT' and 'RANDOM' options were specified in the 'STACKING'
! description.

!      ARGUMENTS:
!           No arguments are used. All data is in 'COMMON'.

!      COMMON VARIABLES:
!            uses:      l_cnt, l_g, n_layers, l_seq, l_alpha

!        modifies:      l_seq

!      GETLAY returns logical .true. if all went well.
! ______________________________________________________________________

      LOGICAL FUNCTION getlay()


!     Utiliza las variables: DIFFaX.par , 'DIFFaX.inc' , okay , messge, i, j, idum
!                             , RAN3, x, sum ,
!     Utiliza las funciones: RAN3 externa
!     Utiliza las subrutinas:

      LOGICAL :: okay
      CHARACTER (LEN=80) :: messge
      INTEGER*4 i, j, idum
      REAL*8  x, sum!, rnd
! external function

      call random_seed()

      getlay = .false.
      okay = .true.

! initialize random numbers in RAN3
     ! call random_number (rnd)
     ! idum = int (rnd*(-100))

      WRITE(*,*) 'Generating a random sequence of ', l_cnt, ' layers.'

! Get the first layer. Even though some stacking transition
! probabilities to and from the layer are non-zero, the layer itself
! may not actually exist anywhere in the crystal! Must use the
! existence probabilities, l_g.
!en definitiva, si capa 1 existeix, comenca per aquesta

    IF (l_seq(1)==0) then
      call random_number (x)
      IF(x == one) x = x - eps7 !perque no sigui exactament 1
      sum = zero
      i = 1
      10 sum = sum + l_g(i)
      IF(x > sum) THEN
        i = i + 1
        IF(i <= n_layers) GO TO 10
! if all goes well, we should not reach here.
        messge = 'GETLAY could not generate the first layer.$'
        GO TO 999
      END IF
      l_seq(1) = i


      ! Now generate the remaining l_cnt-1 layers
      j = 2
   21 if  (l_seq(j)==0)  then
         call random_number (x)
         IF(x == one) x = x - eps7
         sum = zero
          i = 1
          31 sum = sum + l_alpha(i,l_seq(j-1))

          IF(x > sum) THEN                   !perque la probabilitat no sigui zero
            i = i + 1
           IF(i <= n_layers) GO TO 31
! if all goes well, we should not reach here.
             WRITE(messge,101) 'GETLAY could not generate layer ', j, '.$'
             GO TO 999
          END IF
           l_seq(j) = i
           j = j + 1
           IF(j <= nint(l_cnt)) GO TO 21


      elseif  (l_seq(j)/=0) then

          if ( l_alpha(l_seq(j),l_seq(j-1))== 0) then    !cal canviar l_seq j-1

              do

                if ( l_seq(j)/= 0) then
                  l_seq(j-1) = l_seq(j)
                  if  ( l_seq(j+1)== 0) l_seq(j)=0
                 ! if ( l_alpha(j-1,l_seq(j-2))== 0) then
                 ! l_seq(j-2) = l_seq(j-1)
                 ! end if

                  j=j+1
                  cycle
                else
                  exit
                end if
              end do
           j= j -1
           IF(j <= nint(l_cnt)) GO TO 21
          end if

         j = j + 1
         IF(j <= nint(l_cnt)) GO TO 21



      end if
 !****************************************************************************************
    else

! Now generate the remaining l_cnt-1 layers
      j = 2
      do

           if  (l_seq(j)/=0)  then
                 j=j+1
                 cycle
           else

              call random_number (x)
              IF(x == one) x = x - eps7
                sum = zero
                i = 1
             131 sum = sum + l_alpha(i,l_seq(j-1))

                IF(x > sum) THEN                   !perque la probabilitat no sigui zero
                 i = i + 1
                 IF(i <= n_layers) GO TO 131
        ! if all goes well, we should not reach here.
                 WRITE(messge,101) 'GETLAY could not generate layer ', j, '.$'
                 GO TO 999
                END IF
                  l_seq(j) = i
                  j = j + 1

            IF(j .GT. nint(l_cnt)) exit

              41 if  (l_seq(j)==0)  then
                    call random_number (x)
                    IF(x == one) x = x - eps7
                    sum = zero
                    i = 1
                    51 sum = sum + l_alpha(i,l_seq(j-1))

                    IF(x > sum) THEN                   !perque la probabilitat no sigui zero
                       i = i + 1
                       IF(i <= n_layers) GO TO 51
                       WRITE(messge,101) 'GETLAY could not generate layer ', j, '.$'
                       GO TO 999
                    END IF
                    l_seq(j) = i
                    j = j + 1
                    IF(j <= nint(l_cnt)) GO TO 41


                 elseif  (l_seq(j)/=0) then

                    if ( l_alpha(l_seq(j),l_seq(j-1))== 0) then    !cal canviar l_seq j-1
                        do
                           if ( l_seq(j)/= 0) then
                               l_seq(j-1) = l_seq(j)
                               if  ( l_seq(j+1)== 0) l_seq(j)=0
                               j=j+1
                               cycle
                           else
                               exit
                           end if
                        end do
                           j= j -1
                        IF(j <= nint(l_cnt)) GO TO 41
                 end if

                    j = j + 1
                    IF(j .GT. nint(l_cnt)) exit

            end if
            exit
        end if

 !****************************************************************************************
      end do
    end if



      getlay = .true.
      RETURN
      999 WRITE(op,102) messge(1:INDEX(messge,'$')-1)
      RETURN
      100 FORMAT(1X, a, i4, a)
      101 FORMAT(a, i4, a)
      102 FORMAT(1X, 'ERROR: ', a)
      END FUNCTION getlay

! Title: GETSAD
! Author: MMJT
! Date: 29 Oct 1989; 16th April 1999
! Description: This subroutine generates the selected area diffraction
! pattern (SADP). GETSAD returns .true. if all went well. The pattern
! is stored in the linear array 'spec'. 'spec' is read by WRTSAD
! to write the diffraction pattern image to disk as a binary file.
!
!      ARGUMENTS:
!            FN      -  Function name passed by reference. The
!                       choice is between GLQ16 (non-adaptive
!                       Gauss-Legendre integration), and AGLQ16
!                       (adaptive Gauss-Legendre integration). (input).
!            view    -  Choice of beam direction (input).
!                              1  =  normal to the plane k = 0
!                              2  =  normal to the plane h = 0
!                              3  =  normal to the plane h = k
!                              4  =  normal to the plane h = -k
!            l_upper -  Upper limit of l. (input).
!            hk_lim  -  Upper limit of h (or k). (output).
!            infile  -  The name of the input data file. (input).
!            ok      -  logical flag indicating all went well.
!                                                      (output).
!
!      COMMON VARIABLES:
!            uses:      a0, b0, c0, d0, lambda, has_l_mirror, sadblock,
!                       spec, loglin, brightness, X_RAY, rad_type
!
!        modifies:      scaleint
! ______________________________________________________________________
!
      SUBROUTINE getsad(fn, view, l_upper, hk_lim, infile, ok)

!     Utiliza las variables: DIFFaX.par , 'DIFFaX.inc' , FN, l_upper , view, hk_lim, infile ,ok,
!                            h, k, i, j, n, info_step, info, cnt, LENGTH, origin  , x, S,
!                            S_value, ANGLE, W4, PNTINT, theta, Q2, l, l_lower, dl, high1, high2, intervals =TWENTY
!     Utiliza las funciones:  FN, LENGTH, PNTINT  externas ,  S(h,k,l), ANGLE(h,k,l) , W4(theta),
!     Utiliza las subrutinas: XYPHSE, PRE_MAT


      INTEGER*4, INTENT(IN OUT)                :: view
      REAL*8, INTENT(IN out)                   :: l_upper
      INTEGER*4, INTENT(IN OUT)                :: hk_lim
      CHARACTER (LEN=*), INTENT(IN OUT)        :: infile
      LOGICAL, INTENT(IN OUT)                  :: ok

      INTEGER*4 h, k, i, j, n, info_step, info, cnt,  origin, li, n_step
      REAL*8 x, s, s_value, angle, w4,  theta, q2, l
      REAL*8 l_lower, dl, high1, high2
      REAL*8, PARAMETER :: intervals = twenty

! external functions (FN is either GLQ16 or AGLQ16)
      REAL*8 fn
      EXTERNAL fn

! external subroutines (Some compilers need them declared external)
!      external XYPHSE, PRE_MAT

! statement functions
! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
! ANGLE is the Bragg angle (in radians) of the h,k,l plane
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))
! W4 is the X-ray polarization factor
      w4(theta) = half * (one + (COS(two*theta))**2)

      q2 = four / (lambda**2)

! check angles are meaningful
      s_value = s(0,0,l_upper)
      IF(s_value <= zero) THEN
        WRITE(op,101) 'ERROR: Illegal value in GETSAD: 1/d**2 = ',s_value
        ok = .false.
        GO TO 990
      END IF
      IF(s_value > q2) THEN
        WRITE(op,100)
        l_upper = two / (lambda*SQRT(c0)) - eps10
        WRITE(op,101) 'Upper bound reduced to ', l_upper
        s_value = s(0,0,l_upper)
      END IF

! Set h and k limit
      IF(view == 1) THEN
        hk_lim =  INT(l_upper * SQRT(c0 / a0))
      ELSE IF(view == 2) THEN
        hk_lim =  INT(l_upper * SQRT(c0 / b0))
      ELSE IF(view == 3) THEN
        hk_lim =  INT(l_upper * SQRT(c0 / (a0 + b0 + d0)))
      ELSE IF(view == 4) THEN
        hk_lim =  INT(l_upper * SQRT(c0 / (a0 + b0 - d0)))
      END IF

! Get l increment
      dl = l_upper / DBLE(sadsize/2)
! reset l_upper so that our integration window, of width dl,
! straddles l = 0.
      l_upper = l_upper - half*dl

! Get diffraction pattern.
! First check scan limits, and make sure we are not going to overflow
! the array 'spec' which will be used to hold the scan information
! prior to writing to file in WRTSAD.
      IF(has_l_mirror) THEN
        l_lower = zero - half*dl
        sadblock = sadsize/2
        IF((hk_lim + 1)*(sadsize/2) > max_sp) THEN
          hk_lim = max_sp/(sadsize/2) - 1
        END IF
      ELSE
        l_lower = -l_upper
        sadblock = sadsize
        IF((hk_lim + 1)*sadsize > max_sp) THEN
          hk_lim = max_sp/sadsize - 1
        END IF
      END IF
      info_step = nint( l_upper / (dl * intervals) )
      IF(info_step <= 0) info_step = 1

      cnt = 0
      DO  i = 0, hk_lim
! set h and k
        IF(view == 1) THEN
          h =  i
          k =  0
        ELSE IF(view == 2) THEN
          h =  0
          k =  i
        ELSE IF(view == 3) THEN
          h =  i
          k =  i
        ELSE IF(view == 4) THEN
          h =  i
          k = -i
        ELSE
          WRITE(op,105) 'ERROR: Illegal view value in GETSAD', view
          GO TO 990
        END IF

! write progress to screen
        WRITE(op,102) h, k, infile(1:length(infile))

        CALL xyphse(h, k)

        CALL pre_mat(h, k)

! The following piece of monkey business is here because if there is
! no mirror then the l-loop only calculates SADSIZE-1 values. If there
! is a mirror then the l-loop makes SADSIZE/2 calculations. We wish
! the final binary data file to be in a square array SADSIZE x SADSIZE,
! so we pad out the first value with a zero.

        IF(.NOT. has_l_mirror) THEN
          cnt = cnt + 1
          spec(cnt) = zero
        END IF

        info = 0
        n_step=nint((l_upper-eps10-l_lower)/dl+1.0)
        do  li=1, n_step
          l=l_lower+(li-1)*dl
        !DO  l = l_lower, l_upper-eps10, dl
! If we are at the origin, don't integrate, we'll probably overflow.
          IF(i == 0 .AND. ABS(l+dl) <= dl+eps10) THEN
            x = zero
            origin = cnt + 1
          ELSE IF(s(h,k,l+dl) > q2) THEN
            x = zero
          ELSE
            x = fn(h,k,l,l+dl,ok)
            IF(.NOT.ok) GO TO 999
            IF(rad_type == x_ray) x = x*w4(angle(h,k,l+half*dl))
          END IF
          cnt = cnt + 1
! make sure we do not overflow
          IF(cnt > max_sp) GO TO 998
          spec(cnt) = x
          IF(MOD(info,info_step) == 0) THEN
            IF(loglin == 0) THEN
              IF(one+x > zero) THEN
! write out log(1+x), since x may get to be quite small
                WRITE(op,103) 'l = ',l,' log(intensity) = ',LOG(one+x)
              ELSE
                WRITE(op,103) 'l = ', l, ' log(intensity) = ', zero
              END IF
            ELSE
              WRITE(op,103) 'l = ', l, ' intensity = ', x
            END IF
          END IF
          info = info + 1
        END DO
      END DO
! check cnt
      IF(cnt < 2) GO TO 980

! patch the intensity at the origin to put a well-defined peak there
      IF(has_l_mirror) THEN
! origin = 1, make origin slightly bigger than proceeding value
        spec(origin) = (one+eps4) * spec(origin+1)
      ELSE
! make it slightly larger than the biggest adjacent value
        spec(origin) = (one+eps4) * MAX(spec(origin-1),spec(origin+1))
      END IF

! we need to find the second highest peak intensity so as to scale the
! data. The highest intensity should be at the origin, which we discard.
      high1 = zero
      high2 = zero
      DO  i = 0, hk_lim
        DO  j = 1, sadblock-1
          n = i*sadblock + j
! check scaling type. If logarithmic, make sure values are always +ve
          IF(loglin == 0) THEN
            IF(one+spec(n) > zero) THEN
              spec(n) = LOG(one + spec(n))
            ELSE
              spec(n) = zero
            END IF
          END IF
! check if origin is the first value. If so, preset high value.
          IF(n == 1 .AND. origin == 1) THEN
            high1 = spec(origin)
            CYCLE
          END IF
          x = spec(n)
          IF(j == 1) THEN
            IF(x > spec(n+1)) THEN
              IF(x > high1) THEN
                high2 = high1
                high1 = x
              ELSE IF(x > high2) THEN
                high2 = x
              END IF
            END IF
          ELSE
            IF(x > spec(n-1) .AND. x > spec(n+1)) THEN
              IF(x > high1) THEN
                high2 = high1
                high1 = x
              ELSE IF(x > high2) THEN
                high2 = x
              END IF
            END IF
          END IF
        END DO
      END DO

      IF(loglin /= 0 .AND. high2 <= zero) THEN
        WRITE(op,101)  &
            'ERROR in intensity scaling in GETSAD. ''SCALE FACTOR'' = ', high2
        ok = .false.
        GO TO 990
      END IF

      IF(loglin == 0 .AND. high1 <= zero) THEN
        WRITE(op,101)  &
            'ERROR in intensity scaling in GETSAD. ''SCALE FACTOR'' = ', high1
        ok = .false.
        GO TO 990
      END IF

! If logarithmic, scale to the brightest peak
! If linear, scale to the 2nd brightest peak
! Note: Intensity scaling can be modified by the user-defined
! brightness value
      if(loglin.eq.0) then
        scaleint = brightness * (maxsad - ONE) / high1
      else
        scaleint = brightness * (maxsad - ONE) / high2
      endif

      990 RETURN
      980 WRITE(op,105) 'Error in GETSAD: loop counter is too small. cnt = ', cnt
      ok = .false.
      RETURN
      998 WRITE(op,104) 'ERROR in GETSAD: spectrum array overflow at h = ',  &
          h,', k = ',k,', l = ',l
      ok = .false.
      RETURN
      999 WRITE(op,104) 'ERROR in GETSAD at h = ',h,', k = ',k,', l = ',l
      RETURN
      100 FORMAT(1X, 'Upper bound exceeds 180 degrees!')
      101 FORMAT(1X, a, g12.5)
      102 FORMAT(1X, 'h = ', i3, ' k = ', i3, 10X, '''', a, '''')
      103 FORMAT(1X, a, f10.5, a, g12.5)
      104 FORMAT(1X, 2(a, i3), a, f10.5)
      105 FORMAT(1X, a, i3)
      END SUBROUTINE getsad

! ______________________________________________________________________
! Title: GETSPC
! Authors: MWD and MMJT
! Date: 17 Mar 1989; Last tweaked on 7 Mar 1995
! Description: This subroutine calculates the spectrum.

!      ARGUMENTS:
!            FN      -  Function name passed by reference. The
!                       choice is between GLQ16 (non-adaptive
!                       Gauss-Legendre), and AGLQ16
!                       (adaptive Gauss-Legendre). (input).
!            infile  -  The name of the input data file. (input).

!      COMMON VARIABLES:
!            uses:      CFile, ELECTN, NEUTRN, SymGrpNo, X_RAY, a0
!                       any_sharp, b0, bnds_wt, c0, cntrl, d0, d_theta
!                       full_brd, lambda, mltplcty, rad_type, rot_only,
!                       spec, th2_max, th2_min, theta1, theta2

!        modifies:      full_shrp

!      GETSPC returns logical .true. if all went well.
! ______________________________________________________________________

      LOGICAL FUNCTION getspc(fn,infile)

!     Utiliza las variables: DIFFaX.par , 'DIFFaX.inc' , FN , infile,  ok, SHARP, on_bndry, l_axis, shrp
!                            h, k, h_lower, h_upper, k_lower, k_upper ,  m, i, max_indx , LENGTH, S, Q,
!                            theta, tmp, tmp2, tmp3, fact, h_val, k_val , HKANGL, LL, ANGLE, AGLQ16 ,
!                             l, hk_th, x, GLQ16, l_max, min_th, max_th, W1, l1, l0, d_l, INTENS, L_STEP,
!                            W2, W3, l00,f(MAX_L)
!     Utiliza las funciones:  FN, GLQ16, AGLQ16, INTENS, SHARP, L_STEP, LENGTH  externas todas  , S(h,k,l)
!                              LL(theta,h,k) , ANGLE(h,k,l), HKANGL(k_val,h_val),   W1(theta) ,  W2(theta) ,
!                               W3(theta),
!     Utiliza las subrutinas:  XYPHSE, PRE_MAT, GET_F, CHWDTH


      CHARACTER (LEN=*), INTENT(IN OUT)        :: infile

      LOGICAL :: ok,  on_bndry, l_axis, shrp
      INTEGER*4 h, k, h_lower, h_upper, k_lower, k_upper, i_th, i_thm
      INTEGER*4 m, i, max_indx, lz, lzf
      REAL*8 s, q, theta, tmp, tmp2, tmp3, fact, h_val, k_val, tmpa, tmpb, tmpc, tmpd, tmpe, tmpf , tmpg, tmph
      REAL*8  :: hkangl, ll, angle , angles
      REAL*8 l, hk_th, x,  l_max, min_th, max_th
      REAL*8 w1, l1, l0, d_l,   w2, w3, l00
      COMPLEX*16 f(max_l)

! external functions
      real*8 fn
      EXTERNAL fn
! external subroutines (Some compilers need them declared external)
!      external XYPHSE, PRE_MAT, GET_F, CHWDTH

! statement functions
! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
! LL is the maximum allowable l-value for a given h,k and theta
      ll(theta,h,k) = SQRT((fact * (SIN(theta))**2 - h*h*a0 - k*k*b0 - h*k*d0)/ c0)
! ANGLE is the Bragg angle (in radians) of the h,k,l plane
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))
! HKANGL is the angle between the vector (h_val,k_val,0) and (1,0,0)
      hkangl(k_val,h_val) = ATAN2(k_val*SQRT(a0*b0 - d0*d0*quarter),  &
          (h_val*a0 + k_val*d0*half))
! These factors are for powder diffraction patterns
! W1 is the polarization factor for weighting x-ray intensities
! it also incorporates the Lorentz factor
      w1(theta) = (one + COS(two*theta) * COS(two*theta)) /  &
          (SIN(theta) * SIN(two*theta))
! W2 is the neutron weighting factor--the Lorentz factor
      w2(theta) = one / (SIN(theta) * SIN(two*theta))
! W3 is the electron weighting factor--the Lorentz factor
      w3(theta) = one / (SIN(theta) * SIN(two*theta))

      getspc = .false.
      ok = .true.

! Make sure we are within bounds. If not, adjust.
      min_th = half * th2_min
      max_th = half * th2_max
      max_indx = INT((max_th - min_th) / d_theta + 1)
      IF(max_indx > max_sp) THEN
        d_theta = (max_th - min_th) / (max_sp - 1)
        max_indx = INT((max_th - min_th) / d_theta + 1)
        WRITE(op,300) ''
        WRITE(op,250) 'd_theta is too small and has been adjusted to ',  &
            two*d_theta*rad2deg
      END IF

    !  1234 WRITE(op,300) 'Enter 1 for adaptive quadrature over all l values'
    !  WRITE(op,300) 'on rows with "sharp" spots: '
    !  READ(cntrl,*,ERR=1234) full_shrp
    !  IF(cfile) WRITE(op,400) full_shrp

       full_shrp = 1
! zero out spectra
      DO  i = 1, max_sp
        spec(i) = zero
      END DO
! See if there is a chance of any sharp peaks.
! If so, an appropriate step along l is found, and any_sharp is .true.
      d_l = l_step(ok)
      IF(.NOT. ok) GO TO 999
! If no sharp peaks were unambiguously detected, override user.
      IF(d_l == zero) full_shrp = 1
! determine extreme values of h
      q = two * SIN(max_th) / lambda
      fact = two / lambda
      fact = fact * fact
! h_upper is the largest value that the index h can have,
! consistent with the value of Q. (When cell_gamma is not 90 degrees,
! k is not necessarily zero at this extreme).
! In case the incredible happens, immune system to rescue.
      tmp3 = four * a0 * b0 - d0 * d0
      IF(tmp3 <= zero) GO TO 990
      tmp3 = two * q * SQRT(one / tmp3)
      h_upper = INT(tmp3 * SQRT(b0))
      h_lower = -h_upper
      k_upper = INT(tmp3 * SQRT(a0))
      k_lower = -k_upper
      h_min = h_lower
      h_max = h_upper
      k_min = k_lower
      k_max = k_upper


      n_hkl=(h_upper -h_lower+1)* (k_upper -k_lower+1)* 80  !
      if(allocated(dos_theta))   deallocate(dos_theta)
      allocate(dos_theta(n_hkl))
      dos_theta=0.0
      if(allocated(hkl_list))   deallocate(hkl_list)
      allocate(hkl_list(3,n_hkl))
      hkl_list=0
      n_hkl=0

! scan along h-axis from h_lower to h_upper
      DO  h = h_lower, h_upper
! determine limits along k for a given h
        DO  k = k_lower, k_upper
! if out of bounds, cycle
          IF(s(h,k,zero) > q*q) CYCLE
          l_axis = h == 0 .AND. k == 0
          hk_th = theta1
          IF(.NOT. l_axis) hk_th = hkangl(DBLE(k), DBLE(h))
! see if in wedge to be scanned
          IF((theta1-hk_th)*(theta2-hk_th) <= eps3 .OR. symgrpno == 1) THEN
! if rotational symmetry only, do not take points on upper wedge plane
            IF(rot_only .AND. (theta2-hk_th) <= eps3 .AND. symgrpno /= 1) CYCLE
            IF(symgrpno == 11 .AND. .NOT. l_axis) CYCLE
         !   WRITE(op,200) 'Integrating along l at ',h,k,  &
         !       '''',infile(1:length(infile)),''''
            on_bndry = ABS(hk_th-theta1) <= eps3 .OR. ABS(hk_th-theta2) <= eps3
! set up the phases in the structure factors and stacking vectors
! which depend on h and k only
            CALL xyphse(h, k)
            CALL pre_mat(h, k)
! assign a corrected shape-broadening half-width

            IF(finite_width) THEN
              tmp2 = (h + k*COS(pi-cell_gamma))/(wa*SIN(pi-cell_gamma))
              tmp3 = k / wb
              ffhkcnst = quarter*lambda*SQRT(a0*tmp2*tmp2 + b0*tmp3*tmp3)
            END IF
! get starting value of l
            IF(l_axis) THEN   ! h==0 and k==0
              tmp = MIN(d_theta, max_th)
              IF(tmp < min_th) tmp = min_th
              l1 = ll(tmp, h, k)
              shrp = any_sharp

            ELSE

              tmp = angle(h, k, zero)
              IF(tmp < min_th) THEN
                l1 = ll(min_th,h,k)
                tmp = angle(h, k, l1)
              ELSE
                l1 = zero
              END IF
              IF(any_sharp .AND. full_shrp /= 1) THEN
                shrp = sharp(h, k, d_l)
                IF(.NOT. ok) GO TO 999
              ELSE
                shrp = any_sharp
              END IF

            END IF
            call update_reflections(h,k,l1)

! m indexes into the array spec
            m = INT((tmp - min_th) / d_theta) + 1
            IF(.NOT. shrp .OR. full_shrp == 1) THEN
! broad streak or full adaptive integration over sharp spots
             ! IF(full_shrp == 1 .OR. full_brd == 1)  &
               !   WRITE(op,300) 'Full adaptive integration'
! integrate each d_theta's range of reciprocal space
              !DO  theta = tmp, max_th-eps10, d_theta
              i_thm=nint((max_th-eps10-tmp)/d_theta+1.0d0)
              DO  i_th = 1,i_thm
                theta=tmp+(i_th-1)*d_theta
                l0 = l1
                tmp2 = MIN(d_theta, max_th-theta)
                l1 = ll(theta+tmp2, h, k)
! sharp spots; do not use knowledge of where they are
                IF(shrp) THEN
                  x = aglq16(h, k, l0, l1, ok)
                ELSE
! broad streaks
                  x = fn(h, k, l0, l1, ok)
                END IF

                IF(.NOT. ok) GO TO 110

! include weighting factors for radiation type
                IF(rad_type == x_ray) THEN
                  x = two * x * w1(theta + half * tmp2)
                ELSE IF(rad_type == neutrn) THEN
                  x = two * x * w2(theta + half * tmp2)
                ELSE IF(rad_type == electn) THEN
                  x = two * x * w3(theta + half * tmp2)
                ELSE
                  ok = .false.
                  GO TO 130
                END IF

! see if not on l-axis
                IF(.NOT. l_axis) THEN
! apply multiplicity factor
                  x = x * mltplcty
! if on boundary, apply appropriate weighting (mirror vs rotation only)
                  IF(on_bndry) x = x * bnds_wt
                END IF

                call update_reflections(h,k,l1)

                IF(finite_width) THEN
                  CALL chwdth(h,k,l0,l1,x,m,max_indx)
                ELSE
                  spec(m) = spec(m) + x
                END IF
                m = m + 1
              END DO

            ELSE
! line of sharp spots--detuned delta functions
! use knowledge of where spots are
! make sure we do all l values a multiple of d_l
! starting with spot on hk-plane
              WRITE(op,300) 'which is a line of sharp spots'
              l00 = zero
              IF(l_axis) THEN
                l00 = d_l
                do
                 IF(l00 >= th2_min) exit
                 l00 = l00 + d_l
                end do
              END IF
              l_max = ll(max_th, h, k)
! avoid trouble by ignoring l = l_max

              lzf=nint((l_max-l00)/d_l+1.0d0)
              !DO  l = l00, l_max, d_l
              DO  lz = 1, lzf
                l=l00+(lz-1)*d_l
                IF(l == l_max) CYCLE
                theta = angle(h,k,l)
                CALL get_f(f, s(h,k,l), l)
                tmp = intens(f, h, k, l, ok) * eps8
! find width of peak
                x = eps10
                80  IF(.NOT. ok) GO TO 120
                x = two * x
                CALL get_f(f, s(h,k,l+x), l+x)
                IF(intens(f, h, k, l+x, ok) > tmp .AND. x <= eps2*d_l) GO TO 80
                IF(.NOT. ok) GO TO 120
                l0 = MAX(l - x, zero)
                l1 = MIN(l + x, l_max)
                x = aglq16(h, k, l0, l1, ok)
                IF(.NOT. ok) GO TO 110

! include weighting factors for radiation type
                IF(rad_type == x_ray) THEN
                  x = two * x * w1(theta)
                ELSE IF(rad_type == neutrn) THEN
                  x = two * x * w2(theta)
                ELSE IF(rad_type == electn) THEN
                  x = two * x * w3(theta)
                ELSE
                  ok = .false.
                  GO TO 130
                END IF

! see if not on l-axis
                IF(.NOT. l_axis) THEN
! apply multiplicity factor
                  x = x * mltplcty

! if on boundary, apply appropriate weighting (mirror vs rotation only)
                  IF(on_bndry) x = x * bnds_wt
                END IF
                m = INT(theta / d_theta) + 1

                call update_reflections(h,k,l)

                IF(finite_width) THEN
                  CALL chwdth(h,k,l0,l1,x,m,max_indx)
                ELSE
                  spec(m) = spec(m) + x
                END IF
              END DO
            END IF
          END IF
        END DO
      END DO

      getspc = .true.
      RETURN
      110 WRITE(op,300) 'GLQ16 returned error in GETSPC.'
      RETURN
      120 WRITE(op,300) 'INTENS returned error in GETSPC.'
      RETURN
      130 WRITE(op,300) 'ERROR: Radiation type is undefined in GETSPC'
      999 RETURN
      990 WRITE(op,300) 'Illegal cell parameters in GETSPC.'
      WRITE(op,250) '4*a0*b0-d0*d0 = ', four * a0 * b0 - d0 * d0
      RETURN
      200 FORMAT(1X, a, 2I4, 6X, 3A)
      250 FORMAT(1X, a, g12.5)
      300 FORMAT(1X, a)
      400 FORMAT(1X, i3)
      END FUNCTION getspc

      Subroutine update_reflections(h,k,l)
        integer,       intent(in) :: h,k
        real (kind=8), intent(in) :: l
        real (kind=8) :: delta, ela, s, angle
        integer :: j
        integer, dimension(3) :: hkl
        logical :: esta
        ! S is the value of 1/d**2 at hkl
        s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
        ! ANGLE is the Bragg angle (in radians) of the h,k,l plane
        angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))

        delta=abs(real(nint(l),kind=8) - l)
        if( delta < eps2 ) then
        !if( l < eps2 ) then
          ela=nint(l)   ! ela is now zero
          hkl=(/ h, k, nint(ela) /)
          esta=.false.
          do j=1,n_hkl
            if(hkl(1) == hkl_list(1,j)  .and. hkl(2) == hkl_list(2,j) .and. hkl(3) == hkl_list(3,j)) then
              esta=.true.
              exit
            end if
          end do
          if(.not. esta) then
            n_hkl=n_hkl+1
            dos_theta(n_hkl)= angle(h,k,ela) * rad2deg *2
            hkl_list(:,n_hkl)= hkl
          end if
        end if
        return
      End Subroutine update_reflections

!!S GETfiles2
! ______________________________________________________________________
! Title: GET_BDS
! Author: MMJT
! Date: 1 July 1989
! Description:  This routine assigns reciprocal space vectors in the
! h,k plane within which integration can be confined. Weightings are
! assigned to off-axis spot intensities to allow for their
! multiplicity relative to spots on the 00l axis. Spots that occur on
! one of the boundaries are assigned special weighting depending on
! whether or not the boundary is also a mirror plane.

!      ARGUMENTS:
!           No arguments are used. All data is in 'COMMON'.

!      COMMON VARIABLES:
!            uses:      SymGrpNo, h_mirror, k_mirror

!        modifies:      h_start, k_start, h_end, k_end, mltplcty,
!                       bnds_wt, rot_only, pnt_grp
! ______________________________________________________________________

      SUBROUTINE get_bds()
!     Utiliza las variables: DIFFaX.par , 'DIFFaX.inc' ,
!     Utiliza las funciones:
!     Utiliza las subrutinas:
! 360 degrees, no symmetry (-1)
! note, the scan vectors are not used in this instance since there is
! an ambiguity between 0 and 360 degrees. We assign values anyway, so
! that the variables are at least initialized in a controlled way.
      IF(symgrpno == 1) THEN
        h_start  =  one
        k_start  =  zero
        h_end    =  one
        k_end    =  zero
        mltplcty =  one
        bnds_wt  =  one
        rot_only = .true.
      END IF

! 180 degrees, rotation only (2/M, 1st setting)
      IF(symgrpno == 2) THEN
        h_start  =  one
        k_start  =  zero
        h_end    = -one
        k_end    =  zero
        mltplcty =  two
        bnds_wt  =  one
        rot_only = .true.
      END IF

! 180 degrees, vertical mirror (2/M, 2nd setting)
! we need to know which mirror plane.
      IF(symgrpno == 3) THEN
        IF(h_mirror .AND. .NOT.k_mirror) THEN
          h_start  =  one
          k_start  =  zero
          h_end    = -one
          k_end    =  zero
          mltplcty =  two
          bnds_wt  =  half
          rot_only = .false.
        ELSE IF(k_mirror .AND. .NOT.h_mirror) THEN
          h_start  =  zero
          k_start  =  one
          h_end    =  zero
          k_end    = -one
          mltplcty =  two
          bnds_wt  =  half
          rot_only = .false.
        ELSE
! In theory, the following should never be needed, but just in case,
! let's bolster DIFFaX's immune system.
          WRITE(op,400) 'DIFFaX is confused about vertical mirrors.'
          WRITE(op,400) 'To be safe, symmetry is being set to -1'
          symgrpno = 1
          pnt_grp = '-1'
          h_start  =  one
          k_start  =  zero
          h_end    =  one
          k_end    =  zero
          mltplcty =  one
          bnds_wt  =  one
          rot_only = .true.
        END IF
      END IF

! 90 degrees, vertical mirrors (MMM)
      IF(symgrpno == 4) THEN
        h_start  =  one
        k_start  =  zero
        h_end    =  zero
        k_end    =  one
        mltplcty =  four
        bnds_wt  =  half
        rot_only = .false.
      END IF

! 120 degrees, rotation only (-3)
      IF(symgrpno == 5) THEN
        h_start  =  one
        k_start  =  zero
        h_end    = -one
        k_end    =  one
        mltplcty =  three
        bnds_wt  =  one
        rot_only = .true.
      END IF

! 60 degrees, vertical mirrors (-3M)
      IF(symgrpno == 6) THEN
        h_start  =  one
        k_start  =  zero
        h_end    =  zero
        k_end    =  one
        mltplcty =  six
        bnds_wt  =  half
        rot_only = .false.
      END IF

! 90 degrees, rotation (4/M)
      IF(symgrpno == 7) THEN
        h_start  =  one
        k_start  =  zero
        h_end    =  zero
        k_end    =  one
        mltplcty =  four
        bnds_wt  =  one
        rot_only = .true.
      END IF

! 45 degrees, vertical mirrors (4/MMM)
      IF(symgrpno == 8) THEN
        h_start  =  one
        k_start  =  zero
        h_end    =  one
        k_end    =  one
        mltplcty =  eight
        bnds_wt  =  half
        rot_only = .false.
      END IF

! 60 degrees, rotation (6/M)
      IF(symgrpno == 9) THEN
        h_start  =  one
        k_start  =  zero
        h_end    =  zero
        k_end    =  one
        mltplcty =  six
        bnds_wt  =  one
        rot_only = .true.
      END IF

! 30 degrees, vertical mirrors (6/MMM)
      IF(symgrpno == 10) THEN
        h_start  =  one
        k_start  =  zero
        h_end    =  one
        k_end    =  one
        mltplcty =  twelve
        bnds_wt  =  half
        rot_only = .false.
      END IF

! integrate along 0 0 l only (axial)
! the following are somewhat arbitrary in this case. Assign values
! anyway just to make sure they are initialized
      IF(symgrpno == 11) THEN
        h_start  =  one
        k_start  =  zero
        h_end    =  one
        k_end    =  zero
        mltplcty =  one
        bnds_wt  =  one
        rot_only = .true.
      END IF

      RETURN
      400 FORMAT(1X, a)
      END SUBROUTINE get_bds

! ______________________________________________________________________
! Title: GET_F
! Author: MWD and MMJT
! Date: 22 Mar 1989
! Description: This routine calculates the form factors for each layer.
! Since this routine is the main holdup for complex structures, it
! attempts to make use of shortcuts detected in OPTIMZ. The Debye-
! Waller diffuse background is assumed to be small and is not
! calculated in this version.

!      ARGUMENTS:
!            f   -  Array of layer form factors. (output).
!            Q2  -  Value of 1/d**2 at h,k,l. (input).
!            l   -  reciprocal lattice vector l-component. (input).

!      COMMON VARIABLES:
!            uses:  x_sf, n_sf, e_sf, rad_type, n_atoms, l_symmetry,
!                   one_B, l_n_atoms, a_type, hx_ky, a_pos, a_occup,
!                   a_B, l_actual, CENTRO, ELECTN, NEUTRN, X_RAY
!                   n_actual, n_layers

!        modifies:  no COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE get_f(f, s2, l)
!     Utiliza las variables: DIFFaX.par , 'DIFFaX.inc' ,  S2, l ,  f(MAX_L),i, j, m, n, type
!                            fact(MAX_TA), tmp(MAX_TA), tmp_sum, dot, e_factor = 0.023934D0, Q2
!                             ctmp(MAX_TA), f_uniq(MAX_L), ctmp_sum
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      COMPLEX*16,dimension(max_l), INTENT(IN OUT)   :: f !max_l
      REAL*8,                          INTENT(IN)   :: s2
      REAL*8,                          INTENT(IN)   :: l
      INTEGER*4 i, j, m, n, TYPE
      REAL*8 fact(max_ta), tmp(max_ta), tmp_sum, dot,  q2
      REAL*8, PARAMETER :: e_factor = 0.023934D0
      COMPLEX*16 ctmp(max_ta), f_uniq(max_l), ctmp_sum

      q2 = quarter * s2
! Q2 = sin(theta)**2 / lamba**2
! determine scattering factors for each atom type
      IF(rad_type == x_ray .OR. rad_type == electn) THEN
        DO  i = 1, n_atoms
! This empirical formula comes from p. 71 of
! "International Tables for X-ray Crystallography, Vol. IV"
! (The Kynoch Press: Birmingham, England), 1974.
          fact(i) = x_sf(1,i) * EXP(-x_sf(2,i) * q2) +  &
              x_sf(3,i) * EXP(-x_sf(4,i) * q2) + x_sf(5,i) * EXP(-x_sf(6,i) * q2) +  &
              x_sf(7,i) * EXP(-x_sf(8,i) * q2) + x_sf(9,i)
        END DO
      ELSE IF(rad_type == neutrn) THEN
        fact(1:n_atoms) = n_sf(1:n_atoms)
      END IF

! get electron scattering factor from x-ray scattering factor
! s = 4 pi sin(theta) / lambda
! f_electron(s) = (8 pi**2 m e**2 / h**2) {Z - fx(s)} / s**2

!       = 0.023934 lambda**2 {Z - fx(s)} / sin(theta)**2
      IF(rad_type == electn) THEN
        DO  i = 1, n_atoms
          fact(i) = e_factor * ( DBLE(e_sf(i)) - fact(i) ) / q2
        END DO
      END IF

      DO  m = 1, n_actual
        tmp_sum  = zero
        ctmp_sum = c_zero
        DO  n = 1, n_atoms
          tmp(n)  = zero
          ctmp(n) = c_zero
        END DO

! First calculate the scattering factors of the unique layers.
! Check to see if f_uniq(m) will all be real and if Debye-Waller is
! invariant
! Note: hx_ky(j,m) contains h*a_pos(1,j,m) + k*a_pos(2,j,m)

        IF(l_symmetry(m) == centro .AND. one_b(m)) THEN
          DO  j = 1, l_n_atoms(m)
            TYPE = a_type(j,m)
            dot = hx_ky(j,m) + l*a_pos(3,j,m)
            tmp(TYPE) = tmp(TYPE) + a_occup(j,m) * COS(dot)
          END DO
          DO  j = 1, n_atoms
            tmp_sum = tmp_sum + tmp(j) * fact(j)
          END DO
! NOTE: twice since centrosymmetric
          f_uniq(m) = two * EXP(-a_b(1,m) * q2) * tmp_sum

! Debye-Waller is not invariant
        ELSE IF(l_symmetry(m) == centro) THEN
          DO  j = 1, l_n_atoms(m)
            TYPE = a_type(j,m)
            dot = hx_ky(j,m) + l*a_pos(3,j,m)
            tmp(TYPE) = tmp(TYPE) + a_occup(j,m) * EXP(-a_b(j,m) * q2) * COS(dot)

          END DO
          DO  j = 1, n_atoms
            tmp_sum = tmp_sum + tmp(j) * fact(j)
          END DO
! NOTE: twice since centrosymmetric
          f_uniq(m) = two * tmp_sum

! check if Debye-Waller is the only invariant
! f(i) will be complex
        ELSE IF(one_b(m)) THEN
          DO  j = 1, l_n_atoms(m)
            TYPE = a_type(j,m)
            dot = hx_ky(j,m) + l*a_pos(3,j,m)
            ctmp(TYPE) = ctmp(TYPE) + a_occup(j,m) * DCMPLX(COS(dot), SIN(dot))
          END DO
          DO  j = 1, n_atoms
            ctmp_sum = ctmp_sum + ctmp(j) * fact(j)
          END DO
          f_uniq(m) = EXP(-a_b(1,m) * q2) * ctmp_sum


! Nothing is invariant
        ELSE
          DO  j = 1, l_n_atoms(m)
            TYPE = a_type(j,m)
            dot = hx_ky(j,m) + l*a_pos(3,j,m)
            ctmp(TYPE) = ctmp(TYPE) + a_occup(j,m) *  &
                EXP(-a_b(j,m) * q2) * DCMPLX(COS(dot), SIN(dot))

          END DO
          DO  j = 1, n_atoms
            ctmp_sum = ctmp_sum + ctmp(j) * fact(j)
          END DO
          f_uniq(m) = ctmp_sum
        END IF
      END DO

! Now assign scattering factors to all the layers
      DO  i = 1, n_layers
        f(i) = f_uniq(l_actual(i))
      END DO

      RETURN
      END SUBROUTINE get_f

! ______________________________________________________________________
! Title: GET_G
! Author: MWD and MMJT
! Date: 18 Aug 1988; 15 Mar 1995
! Description: This function determines g_i, the a-priori probability
! that a layer of type i, will actually occur somewhere within the
! crystal.
! 'cnt' counts the l_alpha(i,i) = 1.0 terms. Then, l_g(i) = 1.0/cnt.

!      ARGUMENTS:
!            No input arguments.

!      COMMON VARIABLES:
!            uses:   n_layers, l_alpha

!        modifies:   l_g

!      GET_G returns logical .true. if all went well.
! ______________________________________________________________________

      LOGICAL FUNCTION get_g()
!     Utiliza las variables: DIFFaX.par,'DIFFaX.inc,    singular, LUDCMP ,
!                             i, j, cnt, index(MAX_L),  sum, g_mat(MAX_L,MAX_L), Det
!     Utiliza las funciones:LUDCMP    externa
!     Utiliza las subrutinas: LUBKSB

      LOGICAL :: singular
      INTEGER::  i, j, cnt, INDX(max_l)
      REAL(kind=8) :: suma, g_mat(max_l,max_l), det

! external function

! external subroutine (Some compilers need them declared external)
!      external LUBKSB

      get_g = .false.
! set up matrix that represents the probabilities
! only n-1 equations are independent

      DO  i = 1, n_layers - 1
        l_g(i) = zero
        suma = zero
        DO  j = 1, n_layers
          suma = suma + l_alpha(j,i)
        END DO
        suma = one / suma
! suma should actually be ONE
        DO  j = 1, n_layers
          g_mat(i,j) = suma * l_alpha(i,j)
        END DO
        g_mat(i,i) = g_mat(i,i) - one
      END DO

      l_g(n_layers) = one

! the sum of the g's must be 1
      DO  i = 1, n_layers
        g_mat(n_layers, i) = one
      END DO

! before we invert the matrix, let's catch the pathological values
      cnt = 0
      DO  i = 1, n_layers
        IF(l_alpha(i,i) == one) cnt = cnt + 1
      END DO

      IF(cnt /= 0) THEN
        singular = .true.
      ELSE
        singular = .false.
      END IF

      IF(singular) THEN
        DO  i = 1, n_layers

          IF(l_alpha(i,i) == one) THEN
! arbitrarily assume such layers occur with equal probability
            l_g(i) = one / DBLE(cnt)
          ELSE
            l_g(i) = zero
          END IF
        END DO
      ELSE
! solve the matrix
        IF(.NOT. ludcmp(g_mat,INDX,n_layers,max_l,det)) GO TO 100
        CALL lubksb(g_mat,l_g,INDX,n_layers,max_l)
      END IF


! If we are here then all went well
      get_g = .true.

! Are some of the layers non-existent?
      DO  i = 1, n_layers
        IF(l_g(i) < eps6) THEN
          WRITE(op,410) 'WARNING: Layer ', i,  &
              ' does not occur in any significant quantity.'
        END IF
      END DO

      RETURN
      100 WRITE(op,400)  &
          'ERROR: Stacking probabilities give a singular matrix in GET_G'
      RETURN
      400 FORMAT(1X, a)
      410 FORMAT(1X, a, i2, a)
      END FUNCTION get_g

! ______________________________________________________________________
! Title: GETALPHA
! Author: MCC
! Date: 4 Feb 2004
! Description: This function generates a random sequence of stacking probabilities
! ______________________________________________________________________

      LOGICAL FUNCTION get_alpha()

        INTEGER     ::  i, j, k
        REAL*8      ::  suma ,sumo, ran

        get_alpha = .false.

        DO i=1, n_layers
           sumo = zero
           suma = zero
           DO j=1, n_layers-1
             l_alpha(j,i) = zero
             call random_number(ran)
             l_alpha(j,i)= ran
             sumo = sumo + l_alpha(j,i)
             IF(ABS(sumo - one) < eps6)  l_alpha(j,i)=zero
             IF   (sumo > one )  l_alpha(j,i)= one - suma
             suma = suma +l_alpha(j,i)
           END DO
           l_alpha( n_layers, i)= one - suma
        END DO
        get_alpha = .true.
        RETURN
      END FUNCTION get_alpha
! ______________________________________________________________________
! Title: GET_MAT
! Author: MMJT
! Date: 21 Mar 1990; 14th July 1995; 21st July 1997
! Description:  This subroutine calculates the elements of 'mat', the
! stacking transition matrix ie. {alpha * exp(2*pi*u.R)}ij. The h and k
! components were pre-calculated in PRE_MAT. GET_MAT calculates the
! l-component to complete the matrix. However, the mat(i,i) terms
! have yet to have 1 subtracted from them. This is done in GET_S.

!      ARGUMENTS:
!            h   -  reciprocal vector h-component. (input).
!            k   -  reciprocal vector k-component. (input).
!            l   -  reciprocal lattice vector l-component. (input).

!      COMMON VARIABLES:
!            uses:  same_rz, n_layers, same_Bs, Bs_zero, mat1, c0, bc0,
!                   ca0, r_B33, r_B23, r_B31, l_r, PI2, there,
!                   fatsWalla_hk

!        modifies:  mat
! ______________________________________________________________________

      SUBROUTINE get_mat(h, k, l)


!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',  h, k , l , dot, twopi_l, fatsWaller, i, j
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      Integer,       Intent(In) :: h
      Integer,       Intent(In) :: k
      Real(Kind=8),  Intent(In) :: l

      Real(Kind=8) :: dot, twopi_l, fatswaller
      Integer      :: i, j

! set up matrix that represents the sequences
! Note: mat is in 'i,j' order.
      twopi_l = pi2 * l
      IF(same_bs) THEN
        IF(all_bs_zero) THEN
! all the fatsWalla_hk terms equal 1
          DO  i = 1, n_layers
            DO  j = 1, n_layers
              IF(there(j,i)) THEN
                dot = twopi_l * l_r(3,j,i)
                mat(i,j) = mat1(i,j) * DCMPLX( COS(dot), SIN(dot) )
              ELSE
                mat(i,j) = c_zero
              END IF
            END DO
          END DO
        ELSE
! all the fatsWalla terms are identical, but are less than 1
! fatsWalla_hk already contains the h-k component computed in PRE_MAT
          fatswaller = fatswalla_hk * EXP(-(l*(quarter*a_b33*c0*l  &
              + half*(a_b23*bc0*k + a_b31*ca0*h))))
          DO  i = 1, n_layers
            DO  j = 1, n_layers
              IF(there(j,i)) THEN
                dot = twopi_l * l_r(3,j,i)
                mat(i,j) = mat1(i,j) * fatswaller * DCMPLX( COS(dot), SIN(dot) )
              ELSE
                mat(i,j) = c_zero
              END IF
            END DO
          END DO
        END IF
      ELSE
! the fatsWalla terms differ. mat1 already contains the h-k component
        DO  i = 1, n_layers
          DO  j = 1, n_layers
            IF(there(j,i)) THEN
              dot = twopi_l * l_r(3,j,i)
              IF(bs_zero(j,i)) THEN
                mat(i,j) = mat1(i,j) * DCMPLX( COS(dot), SIN(dot) )
              ELSE
                mat(i,j) = mat1(i,j) * DCMPLX( COS(dot), SIN(dot) )  &
                    * EXP( -(l*(quarter*r_b33(j,i)*c0*l  &
                    + half*(r_b23(j,i)*bc0*k + r_b31(j,i)*ca0*h))) )
              END IF
            ELSE
              mat(i,j) = c_zero
            END IF
          END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE get_mat

! ______________________________________________________________________
! Title: GET_S
! Author: MWD and MMJT
! Date: 5th Aug 1991
! Description:  This function determines the S's--the average scattered
! wave functions of each layer at h, k, l for crystals with an
! infinite number of layers. GET_MAT should be called prior to GET_S.
! Note, since GET_S is called billions of times from deep within the
! inner loops of DIFFaX's bowels, part of the matrix mat1 has been
! precalculated in PRE_MAT in order to speed things up a bit.

!      ARGUMENTS:
!            f   -  Array of layer form factors. (input).
!            s   -  Array of average layer wavefunctions. (output).
!            h   -  reciprocal vector h-component. (input).
!            k   -  reciprocal vector k-component. (input).
!            l   -  reciprocal lattice vector l-component. (input).

!      COMMON VARIABLES:
!            uses:  n_layers

!        modifies:  mat

!      GET_S returns logical .true. if all went well.
! ______________________________________________________________________

      LOGICAL FUNCTION get_s(f, s, h, k, l)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',  h, k , l ,  f(MAX_L), s(MAX_L)
!                             i_ok, index(MAX_L), i, Det, s_tmp(2)
!     Utiliza las funciones:
!     Utiliza las subrutinas: CGEFA, CGESL

      COMPLEX*16, INTENT(IN)                   :: f(max_l)
      COMPLEX*16, INTENT(IN OUT)               :: s(max_l)
      INTEGER*4, INTENT(IN )                :: h
      INTEGER*4, INTENT(IN )                :: k
      REAL*8, INTENT(IN )                   :: l

! i_ok is used by Linpack routines
      INTEGER :: i_ok, INDEX(max_l)
      INTEGER*4 i
      COMPLEX*16 det, s_tmp(2)
! external subroutines (Some compilers need them declared external)
! CGEFA and CGESL are Linpack routines


      get_s = .false.

! subtract identity matrix (need do this for diagonal terms only).
      DO  i = 1, n_layers
        mat(i,i) = mat(i,i) - c_one
      END DO

! Now solve the system of equations.
      IF(n_layers > 2) THEN
! now call LINPACK routines
        CALL cgefa(mat, max_l, n_layers, INDEX, i_ok)
        IF(i_ok /= 0) GO TO 999
        CALL cgesl(mat, max_l, n_layers, INDEX, s, 0)
      ELSE IF(n_layers == 2) THEN
! its a simple 2 x 2, so solve it directly
        det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
        IF(det == c_zero) GO TO 999
! copy s (remember, if infinitely thick, s = -f)
        s_tmp(1) = -s(1)
        s_tmp(2) = -s(2)
        s(1) = ( mat(1,2) * s_tmp(2) - mat(2,2) * s_tmp(1) ) / det
        s(2) = ( mat(2,1) * s_tmp(1) - mat(1,1) * s_tmp(2) ) / det
      ELSE IF(n_layers == 1) THEN
! only one layer, so solve it immediately
        s(1) = -f(1) / mat(1,1)
      END IF

      get_s = .true.
      RETURN

  999 CONTINUE
      WRITE(op,"(1X,a)") 'GET_S: Solving for sequence produces a singular matrix.'
      WRITE(unit=op,fmt="(1X,a,i3, a, i3, a, g12.5)") 'h = ',h, ' k = ',k, ' l = ',l
      DO  i = 1, n_layers
        WRITE(op,"(5X,a,i2,a,g12.5,a, g12.5,a)") 'f(', i, ') = (', real(f(i)), ',' ,aimag(f(i)),')'
      END DO
      RETURN
      END FUNCTION get_s

! ______________________________________________________________________
! Title: GET_S2
! Author: MMJT
! Date: 5 Feb 1990
! Description:  This function determines the S's--the average scattered
! wave functions of each layer at h, k, l for crystals with only a
! finite number of layers, l_cnt. The equation being solved is

!   inv(Ident-T) * ((N+1)*Ident - inv(Ident-T)*(Ident-T**(N+1)) * F / N

!  where N = l_cnt, and T is the stacking probability matrix

!      ARGUMENTS:
!            f   -  Array of layer form factors. (input).
!            s   -  Array of average layer wavefunctions. (output).
!            h   -  reciprocal vector h-component. (input).
!            k   -  reciprocal vector k-component. (input).
!            l   -  reciprocal lattice vector l-component. (input).

!      COMMON VARIABLES:
!            uses:  n_layers, l_cnt

!        modifies:  mat

!      GET_S2 returns logical .true. if all went well.
! ______________________________________________________________________

      LOGICAL FUNCTION get_s2(f, s, h, k, l)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',  h, k , l , f(MAX_L), s(MAX_L)
!                              ok, GET_S, MAT2N ,  i, j  , ctmp, mat_n(MAX_L,MAX_L), tmp_mat(MAX_L,MAX_L)
!     Utiliza las funciones:   GET_S, MAT2N externas
!     Utiliza las subrutinas:

      Complex(Kind=8), Intent(In)      :: f(max_l)
      Complex(Kind=8), Intent(In Out)  :: s(max_l)
      Integer,         Intent(In)      :: h
      Integer,         Intent(In)      :: k
      Real(Kind=8),    Intent(In)      :: l

      Logical         :: ok
      Integer         ::  i, j
      Complex(Kind=8) :: ctmp, mat_n(max_l,max_l), tmp_mat(max_l,max_l)
! external functions
      get_s2 = .false.

! get matrix mat multiplied by itself l_cnt+1 times
      ok = mat2n(mat_n)
      IF(.NOT.ok) GO TO 990

! subtract identity matrix, and make a copy of mat.
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          tmp_mat(j,i) = mat(j,i)
        END DO
        mat_n(i,i) = mat_n(i,i) - c_one
      END DO

! Multiply out -(Ident - T**(l_cnt+1))F.
      DO  i = 1, n_layers
        ctmp = c_zero
        DO  j = 1, n_layers
          ctmp = ctmp + mat_n(i,j) * f(j)
        END DO
        s(i) = ctmp
      END DO
! Next, solve. ie. inv(Ident - T) * (Ident - T**(l_cnt+1))*F
! where now s = (Ident - T**(l_cnt+1))*F
      ok = get_s(f, s, h, k, l)
      IF(.NOT.ok) GO TO 999

! Use result to build a new vector, and solve again.
! First, reconstruct mat.
      DO  i = 1, n_layers
        s(i) = (s(i) - f(i) * DBLE(l_cnt + 1)) / DBLE(l_cnt)
        DO  j = 1, n_layers
          mat(j,i) = tmp_mat(j,i)
        END DO
      END DO
! Solve with new RHS vector
      ok = get_s(f, s, h, k, l)
      IF(.NOT.ok) GO TO 999

      get_s2 = .true.
      RETURN
      990 WRITE(op,400) 'MAT2N returned an error in GET_S2.'
      WRITE(op,401) h, k, l
      RETURN
      999 WRITE(op,400) 'Solving for sequence produces a singular matrix.'
      WRITE(op,401) h, k, l
      DO  i = 1, n_layers
        WRITE(op,402) i, f(i)
      END DO
      RETURN
      400 FORMAT(1X, 'GET_S2:', a)
      401 FORMAT(1X, 'h = ', i3, ' k = ', i3, ' l = ', g12.5)
      402 FORMAT(5X, 'f(', i2, ') = (', g12.5, ',', g12.5, ')')
      END FUNCTION get_s2

! ______________________________________________________________________
! Title: GET_SYM
! Author: MMJT
! Date: 19 June 1990
! Determines the symmetry of the diffraction pattern, thereby
! defining the smallest volume of reciprocal space which needs to
! be integrated over. There are only 10 kinematical diffraction point
! groups in the presence of streaking. Friedel's law ensures there
! is always a center of symmetry, and the possibility of streaking
! prohibits the cubic point groups. The point group symmetry number
! is returned as GET_SYM.
! The 10 groups are:

!              GET_SYM          point group
!           -------------------------------
!                 1:       |        -1
!                 2:       |        2/M(1) (diad along streaks)
!                 3:       |        2/M(2) (mirror contains streaks)
!                 4:       |        MMM
!                 5:       |        -3
!                 6:       |        -3M
!                 7:       |        4/M
!                 8:       |        4/MMM
!                 9:       |        6/M
!                10:       |        6/MMM

! The point group symbol is returned in the global character string
! 'pnt_grp'.

!      ARGUMENTS:
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  cell_gamma, DoSymDump, cell_a, cell_b, PI, PI2
!                   RAD2DEG

!        modifies:  max_var, pnt_grp, h_mirror, k_mirror

!      GET_SYM returns one of the ten symmetry flags listed above.
! ______________________________________________________________________

      INTEGER*4 FUNCTION get_sym(ok)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',ok , diad, triad, tetrad, hexad
!                              TST_ROT, TST_MIR , cell90, cell120, eq_sides, idum, rot_sym, tmp_var
!     Utiliza las funciones:   TST_ROT, TST_MIR    externas
!     Utiliza las subrutinas:

      LOGICAL, INTENT(IN OUT)                  :: ok
      LOGICAL :: diad, triad, tetrad, hexad
      LOGICAL :: cell90, cell120, eq_sides
      INTEGER*4 idum, rot_sym
      REAL*8 tmp_var

! external functions

! external subroutine (Some compilers need them declared external)

! initialize random numbers in RAN3
      idum = -1

! initialize function
      get_sym = 0

      max_var = zero
      diad   = .false.
      triad  = .false.
      tetrad = .false.
      hexad  = .false.
      cell90 =  ABS(cell_gamma - half*pi) < half*pi*eps6
      cell120 = ABS(cell_gamma - pi2/three) < pi2*eps6/three

! sample reciprocal space to get an idea of the sort of intensities
! that are out there.
! test for vertical mirror (equivalent to a 2-fold, Friedel's law)
      tmp_var = max_var
      rot_sym = 2
      diad = tst_rot(rot_sym, idum, ok)
      IF(.NOT.ok) GO TO 997
      IF(.NOT.diad) max_var = tmp_var

! if the cell angle is neither 90 nor 120 degrees, then the symmetry
! has to be either -1 or 2/M.
      IF( .NOT.cell90 .AND. .NOT.cell120 ) THEN
        IF(dosymdump) THEN
          WRITE(sy,'(a)') ' '
          WRITE(sy,220) cell_gamma * rad2deg
          WRITE(sy,221)
        END IF

        IF(diad) THEN
          get_sym = 2
          pnt_grp = '2/M(1)'
        ELSE
          get_sym = 1
          pnt_grp = '-1'
        END IF
        GO TO 900
      END IF
      eq_sides = ABS(cell_a - cell_b) <= half*eps6*(cell_a + cell_b)
      IF(.NOT.eq_sides .AND. dosymdump) THEN
        WRITE(sy,'(a)') ' '
        WRITE(sy,225) cell_a, cell_b
        WRITE(sy,226)
        WRITE(sy,221)
      END IF

! cell_gamma is either 90 or 120 degrees.
! if cell_a = cell_b, higher rotational symmetry is possible
      IF(eq_sides) THEN
! Note, it is quite possible for an oblique cell (whose cell_gamma is
! not equal to 120 degrees) to have 3-fold symmetry. We do not test
! for this.
        tmp_var = max_var
        IF(cell120) THEN
          rot_sym = 3
          triad = tst_rot(rot_sym, idum, ok)
          IF(.NOT.ok) GO TO 997
        ELSE
          triad = .false.
        END IF
        IF(.NOT.triad) max_var = tmp_var
        hexad = diad.AND.triad
        IF(hexad.AND.dosymdump) WRITE(sy,200)
        IF(diad.AND.cell90 .AND. (.NOT.triad)) THEN
          tmp_var = max_var
          rot_sym = 4
          tetrad = tst_rot(rot_sym, idum, ok)
          IF(.NOT.ok) GO TO 997
          IF(.NOT.tetrad) max_var = tmp_var
        ELSE
          tetrad = .false.
        END IF
      END IF

! Now test for mirrors.
      tmp_var = max_var
      h_mirror = tst_mir(1, idum, ok)
      IF(.NOT.ok) GO TO 998
      IF(.NOT.h_mirror) max_var = tmp_var
      tmp_var = max_var
      k_mirror = tst_mir(2, idum, ok)
      IF(.NOT.ok) GO TO 999
      IF(.NOT.k_mirror) max_var = tmp_var
      tmp_var = max_var
      hk_mirror = tst_mir(3, idum, ok)
      IF(.NOT.ok) GO TO 999
      IF(.NOT.hk_mirror) max_var = tmp_var

! If h_mirror does not equal k_mirror, then there cannot be a higher
! rotation symmetry than a diad. If, by some bizarre freak, this is
! inconsistent with the result of TST_ROT, choose the lower symmetry.
      IF(h_mirror.NEQV.k_mirror) THEN
        IF(diad) THEN
          get_sym = 2
          pnt_grp = '2/M(1)'
        ELSE
          get_sym = 3
          pnt_grp = '2/M(2)'
        END IF
        GO TO 900
      END IF
! Now check for combinations of mirrors and rotation axes.

! 6-fold
      IF(hexad) THEN
        IF(h_mirror.OR.hk_mirror) THEN
          get_sym = 10
          pnt_grp = '6/MMM'
        ELSE
          get_sym = 9
          pnt_grp = '6/M'
        END IF
        GO TO 900
      END IF
! 4-fold
      IF(tetrad) THEN
        IF(h_mirror.OR.hk_mirror) THEN
          get_sym = 8
          pnt_grp = '4/MMM'
        ELSE
          get_sym = 7
          pnt_grp = '4/M'
        END IF
        GO TO 900
      END IF
! 3-fold
      IF(triad.AND.(.NOT.diad)) THEN
        IF(h_mirror.OR.hk_mirror) THEN
          get_sym = 6
          pnt_grp = '-3M'
        ELSE
          get_sym = 5
          pnt_grp = '-3'
        END IF
        GO TO 900
      END IF
! 2-fold
      IF(diad) THEN
! handle special case of non-orthogonal mesh which has a diad,
! no triad, and one mirror. Diad prevails, vertical mirrors ignored.
        IF((h_mirror.OR.hk_mirror).AND.cell90) THEN
          get_sym = 4
          pnt_grp = 'MMM'
        ELSE
          get_sym = 2
          pnt_grp = '2/M(1)'
        END IF
        GO TO 900
      END IF
! if no symmetry has been detected opt for lowest symmetry
      get_sym = 1
      pnt_grp = '-1'

      900 RETURN
      997 WRITE(op,228) 'ERROR in GET_SYM: error returned by TST_ROT'
      WRITE(op,229) '      while testing for ', rot_sym, '-fold axis.'
      RETURN
      998 WRITE(op,228) 'ERROR in GET_SYM: error returned by TST_MIR'
      WRITE(op,228) '   while testing for mirror about the a - c plane'
      RETURN
      999 WRITE(op,228) 'ERROR in GET_SYM: error returned by TST_MIR'
      WRITE(op,228) '   while testing for mirror about the b - c plane'
      RETURN
      200 FORMAT(1X,'THE 2-FOLD AND 3-FOLD IMPLY 6-FOLD ROTATION SYMMETRY')
      220 FORMAT(1X, 'Cell_gamma = ', g12.5, ' degrees')
      221 FORMAT(1X, 'Rotational symmetry higher than 2-fold is unlikely.')
      225 FORMAT(1X, 'cell-a = ', g12.5, ' cell-b = ', g12.5)
      226 FORMAT(1X, 'Cell sides are not equal.')
      228 FORMAT(1X, a)
      229 FORMAT(1X, a, i2, a)
      END FUNCTION get_sym

!!S Gfiles1
! ______________________________________________________________________
! Title: GAUSSN
! Author: MMJT
! Date: 17 Feb 1990; 7 Mar 1995
! Description: This subroutine simulates Gaussian instrumental
! broadening. std_dev is in degrees. The algorithm used does not
! conserve intensity when std_dev is comparable to d_theta. Intensities
! at the extreme ends of the spectrum are corrupted slightly.

!      ARGUMENTS:
!            th2_low  -  lowest 2theta angle to consider. (input).

!      COMMON VARIABLES:
!            uses:  NONE, PI2, RAD2DEG, blurring, d_theta, FWHM
!                   th2_max

!        modifies:  brd_spc, spec

! ______________________________________________________________________

      SUBROUTINE gaussn(th2_low)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',th2_low, i, j, n_low, n_high, m ,
!                            k1, k2, k3, const, gss, std_dev, tmp, tmp1, tmp2
!
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      REAL*8, INTENT(IN OUT)                       :: th2_low

      INTEGER*4 i, j, n_low,  m
      REAL*8 k1, k2, k3, const, gss, std_dev, tmp, tmp1, tmp2

      IF(fwhm <= zero) GO TO 999
      std_dev = fwhm / SQRT(eight * LOG(two))

! check that cut-off is reasonable
      IF(th2_low < zero .OR. th2_low >= th2_max) THEN
        WRITE(op,101) 'GAUSSN: Cut-off angle ', th2_low,  &
            ' is out of bounds. Angle reset to zero.'
        th2_low = zero
      END IF

! th2_low is the angle relative to th2_min
! 2*d_theta is the angular step size
      n_low  = INT(half*th2_low/d_theta) + 1

      const = two * rad2deg * d_theta
      k1 = const / (SQRT(pi2) * std_dev)
      k2 = half * (const / std_dev)**2

      DO  i = 1, n_high
        brd_spc(i) = zero
      END DO

! go out to 40 standard deviations, or to the end of the spectrum
      m = nint(two * twenty * std_dev / const)
      IF(m > n_high) m = n_high
      DO  i = 0, m
        k3 = k2*DBLE(i*i)
        gss = k1*EXP(-k3)
        DO  j = n_low+1, n_high
          tmp1 = zero
          tmp2 = zero
          IF((j-i) > n_low)  tmp1 = spec(j-i)
          IF((j+i) <= n_high) tmp2 = spec(j+i)
          tmp = tmp1 + tmp2
          IF(i == 0) tmp = half * tmp
          brd_spc(j) = brd_spc(j) + gss * tmp
        END DO
      END DO
      RETURN
      999 WRITE(op,101) 'Illegal FWHM ', fwhm, ' in GAUSSN.'
      WRITE(op,100)'Gaussian instrumental broadening not added'
! kill blurring option
      blurring = NONE
      RETURN
      100 FORMAT(1X, a)
      101 FORMAT(1X, a, g12.5, a)
      END SUBROUTINE gaussn

!!S Gfiles2
! ______________________________________________________________________
! Title: GLQ16
! Authors: MWD and MMJT
! Date: 6 April 1989; 13th July 95
!  This routine performs 16-point Gauss-Legendre quadrature on an
!  interval in reciprocal space. The interval is (h,k,a) to (h,k,b).
!  The routine calls INTENS at each of the 16 points, and if all goes
!  well, returns .true..
!  In the interests of speed, this routine has the option of calling
!  APPR_F, which returns interpolated f values for each of the 16
!  points. This modification slows down the procedure for structures
!  with few atoms per layer, but speeds it up if there are many atoms
!  per layer.

!      ARGUMENTS:
!            h   -  reciprocal lattice vector h-component. (input).
!            k   -  reciprocal lattice vector k-component. (input).
!            a   -  l-value of the lower bound of reciprocal
!                   lattice integration region. (input).
!            b   -  l-value of the upper bound of reciprocal
!                   lattice integration region. (input).
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, recrsv

!        modifies:  intp_F

!      GLQ16 returns the integrated value.
! ______________________________________________________________________

      REAL*8 FUNCTION glq16(h, k, a, b, ok)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc', ok, h, k , a, b  , o, too_close ,intens, inten2
!                            o, too_close, c1, c2, x1 = 0.095012509837637440185D0,
!                            x2 = 0.281603550779258913230D0 , x3 = 0.458016777657227386342D0 ...
!                            x8 = 0.989400934991649932596D0, w1 = 0.189450610455068496285D0...
!                            w8 = 0.027152459411754094852D0 , i, j , n=3, list(n),  Q2, l , ag_l(16), samp_l(n)
!                            f(MAX_L,16)
!     Utiliza las funciones:   INTENS, INTEN2   externas  , Q2(h,k,l)
!     Utiliza las subrutinas: APPR_F, GET_F  ,
      INTEGER*4               :: h
      INTEGER*4               :: k
      REAL*8                  :: a
      REAL*8                  :: b
      LOGICAL, INTENT(IN OUT) :: ok

      LOGICAL :: o, too_close

      REAL*8 c1, c2

      REAL*8, PARAMETER :: x1 = 0.095012509837637440185D0
      REAL*8, PARAMETER :: x2 = 0.281603550779258913230D0
      REAL*8, PARAMETER :: x3 = 0.458016777657227386342D0
      REAL*8, PARAMETER :: x4 = 0.617876244402643748447D0
      REAL*8, PARAMETER :: x5 = 0.755404408355003033895D0
      REAL*8, PARAMETER :: x6 = 0.865631202387831743880D0
      REAL*8, PARAMETER :: x7 = 0.944575023073232576078D0
      REAL*8, PARAMETER :: x8 = 0.989400934991649932596D0

      REAL*8, PARAMETER :: w1 = 0.189450610455068496285D0
      REAL*8, PARAMETER :: w2 = 0.182603415044923588867D0
      REAL*8, PARAMETER :: w3 = 0.169156519395002538189D0
      REAL*8, PARAMETER :: w4 = 0.149595988816576732081D0
      REAL*8, PARAMETER :: w5 = 0.124628971255533872052D0
      REAL*8, PARAMETER :: w6 = 0.095158511682492784810D0
      REAL*8, PARAMETER :: w7 = 0.062253523938647892863D0
      REAL*8, PARAMETER :: w8 = 0.027152459411754094852D0

      INTEGER*4 i, j
! f is approximated by a polynomial of order (n-1)

      INTEGER*4, PARAMETER :: n = 3
      INTEGER*4 list(n)

      REAL*8 q2, l
      REAL*8 ag_l(16), samp_l(n)
      COMPLEX*16 f(max_l,16)

! external functions

! external subroutines (Some compilers need them declared external)
! statement function
! Q2 is the value of 1/d**2 at hkl
      q2(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0

! initialize, to suppress compiler warnings
      glq16 = zero

! check that the integration range is legitimate
      IF(b < a) GO TO 999
! catch pathological values of a and b - i.e. zero integration range!
      IF(b == a) THEN
        ok = .true.
        GO TO 900
      END IF

! let's interpolate ('hard-wired' in this version)
      intp_f = .true.

      c1 = half*(b - a)
      c2 = c1 + a

! set up the 16 l-values to sample
      ag_l(1) = -c1*x8 + c2
      ag_l(2) = -c1*x7 + c2
      ag_l(3) = -c1*x6 + c2
      ag_l(4) = -c1*x5 + c2
      ag_l(5) = -c1*x4 + c2
      ag_l(6) = -c1*x3 + c2
      ag_l(7) = -c1*x2 + c2
      ag_l(8) = -c1*x1 + c2

      ag_l( 9) = c1*x1 + c2
      ag_l(10) = c1*x2 + c2
      ag_l(11) = c1*x3 + c2
      ag_l(12) = c1*x4 + c2
      ag_l(13) = c1*x5 + c2
      ag_l(14) = c1*x6 + c2
      ag_l(15) = c1*x7 + c2
      ag_l(16) = c1*x8 + c2

      IF(intp_f) THEN

! choose special values to sample (3 point interpolation)
        list(1) = 1
        list(2) = 8
        list(3) = 16
        samp_l(1) = ag_l(list(1))
        samp_l(2) = ag_l(list(2))
        samp_l(3) = ag_l(list(3))

! Deal with very rare cases when the spread in l values is too small
        too_close = (samp_l(1) == samp_l(3)) .OR. (samp_l(1) == samp_l(2)) .OR.  &
            (samp_l(2) == samp_l(3))

        IF(.NOT.too_close) THEN
          CALL appr_f(f, h, k, samp_l, ag_l, n, list, ok)
          IF(.NOT.ok) GO TO 990
        ELSE
! sample f values once, and set all 16 values over l-range to be equal
          CALL get_f(f(:,1), q2(h,k,ag_l(1)), ag_l(1))
          DO  j = 1, n_layers
            DO  i = 2, 16
              f(j,i) = f(j,1)
            END DO
          END DO
        END IF

      ELSE
! do not interpolate
        DO  i = 1, 16
          CALL get_f(f(:,i), q2(h,k,ag_l(i)), ag_l(i))
        END DO

      END IF

! sum intensities

      o = .true.
      IF(recrsv) THEN
        glq16 = c1 * (  &
            w8*(intens(f(1,1),h,k,ag_l(1),o)+intens(f(1,16),h,k,ag_l(16),o))+  &
            w7*(intens(f(1,2),h,k,ag_l(2),o)+intens(f(1,15),h,k,ag_l(15),o))+  &
            w6*(intens(f(1,3),h,k,ag_l(3),o)+intens(f(1,14),h,k,ag_l(14),o))+  &
            w5*(intens(f(1,4),h,k,ag_l(4),o)+intens(f(1,13),h,k,ag_l(13),o))+  &
            w4*(intens(f(1,5),h,k,ag_l(5),o)+intens(f(1,12),h,k,ag_l(12),o))+  &
            w3*(intens(f(1,6),h,k,ag_l(6),o)+intens(f(1,11),h,k,ag_l(11),o))+  &
            w2*(intens(f(1,7),h,k,ag_l(7),o)+intens(f(1,10),h,k,ag_l(10),o))+  &
            w1*(intens(f(1,8),h,k,ag_l(8),o)+intens(f(1, 9),h,k,ag_l( 9),o)))

      ELSE

        glq16 = c1 * (  &
            w8*(inten2(f(1,1),h,k,ag_l(1),o)+inten2(f(1,16),h,k,ag_l(16),o))+  &
            w7*(inten2(f(1,2),h,k,ag_l(2),o)+inten2(f(1,15),h,k,ag_l(15),o))+  &
            w6*(inten2(f(1,3),h,k,ag_l(3),o)+inten2(f(1,14),h,k,ag_l(14),o))+  &
            w5*(inten2(f(1,4),h,k,ag_l(4),o)+inten2(f(1,13),h,k,ag_l(13),o))+  &
            w4*(inten2(f(1,5),h,k,ag_l(5),o)+inten2(f(1,12),h,k,ag_l(12),o))+  &
            w3*(inten2(f(1,6),h,k,ag_l(6),o)+inten2(f(1,11),h,k,ag_l(11),o))+  &
            w2*(inten2(f(1,7),h,k,ag_l(7),o)+inten2(f(1,10),h,k,ag_l(10),o))+  &
            w1*(inten2(f(1,8),h,k,ag_l(8),o)+inten2(f(1, 9),h,k,ag_l( 9),o)))

      END IF
      ok = o
      900 RETURN
      999 ok = .false.
      WRITE(op,100) 'GLQ16: Illegal integration interval!'
      WRITE(op,101) h, k, a, b
      RETURN
      990 ok = .false.
      WRITE(op,100) 'GLQ16: ERROR returned from APPR_F'
      WRITE(op,101) h, k, a, b
      RETURN
      100 FORMAT(1X, a)
      101 FORMAT(1X, 'h = ',i4,', k = ',i4,', l0 = ',g12.5,', l1 = ',g12.5)
      END FUNCTION glq16

! ______________________________________________________________________
! Title: GOINTR
! Author: MMJT
! Date: 23 Oct 1989
! Description: This subroutine sets up the parameters necessary to
! integrate over an interval of reciprocal space.

!      ARGUMENTS:
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  cntrl, CFile

!        modifies:  no COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE gointr(ok)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc', ok, i, AGLQ16, GLQ16 ,
!
!     Utiliza las funciones: AGLQ16, GLQ16   externas
!     Utiliza las subrutinas:INTEGR

      LOGICAL, INTENT(IN OUT)                  :: ok
      INTEGER*4 i
! external functions

! external subroutine (Some compilers need them declared external)
!      external INTEGR

      1 WRITE(op,100) ' '
      WRITE(op,100) 'INTEGRATING INTENSITY IN AN INTERVAL ALONG l. . .'
      WRITE(op,100) 'Enter 1 for adaptive quadrature.'

      READ(cntrl,*,ERR=1) i
      IF(cfile) WRITE(op,101) i
      IF(i == 1) THEN
        CALL integr(aglq16,ok)
      ELSE
        CALL integr(glq16,ok)
      END IF
      IF(.NOT.ok) GO TO 999
      i = 1
      2 WRITE(op,100) 'Enter 1 to integrate another interval.'
      READ(cntrl,*,ERR=2) i
      IF(cfile) WRITE(op,101) i
      IF(i == 1) GO TO 1

      999 RETURN
      100 FORMAT(1X, a)
      101 FORMAT(1X, i3)
      END SUBROUTINE gointr

! ______________________________________________________________________
! Title: GOSADP
! Author: MMJT
! Date: 23 Oct 1989; 1st Mar 2000
! Description: This subroutine sets up a file whose name is given
! by 'outfile', which is derived from the input filename 'infile'.
! The selected area diffraction pattern (SADP) is calculated and is
! then written in binary format to 'outfile'.

!      ARGUMENTS:
!            infile   -  name of input file (input).
!            outfile  -  name of output file (output).
!            ok       -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  cntrl, CFile

!        modifies:  loglin, brightness
! ______________________________________________________________________

      SUBROUTINE gosadp(infile,outfile,ok)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',  infile, outfile , ok ,i_adapt,
!                              i_plane, io_err, hk_lim ,AGLQ16, GLQ16, l_upper
!     Utiliza las funciones: AGLQ16, GLQ16  externas
!     Utiliza las subrutinas: STREAK, GETFNM, GETSAD, WRTSAD

      CHARACTER (LEN=*), INTENT(IN OUT)        :: infile
      CHARACTER (LEN=*), INTENT(OUT)           :: outfile
      LOGICAL, INTENT(IN OUT)                  :: ok

      INTEGER*4 i_adapt, i_plane, io_err, hk_lim
      REAL*8   l_upper

! external functions

! external subroutines (Some compilers need them declared external)
!      external STREAK, GETFNM, GETSAD, WRTSAD

! Get a suitable file name, and open file.
      CALL getfnm(infile,outfile,'sadp',ok)
      IF(.NOT.ok) GO TO 999
! Open unformatted for binary write.
      IF(sa /= op) OPEN(UNIT = sa, FILE = outfile, STATUS = 'new',  &
          FORM = 'unformatted', ERR = 990, IOSTAT = io_err)
! some compilers need the following added for true binary output
!     |       ,recordtype = 'stream'

      WRITE(op,100) ' '
    !  WRITE(op,100)'CALCULATING SELECTED AREA DIFFRACTION PATTERN. . .'

      1234 WRITE(op,100) 'Enter 1 for adaptive quadrature.'
      READ(cntrl,*,ERR=1234) i_adapt
      IF(cfile) WRITE(op,101) i_adapt

      1 WRITE(op,100) 'Choose a plane in reciprocal space to view.'
      WRITE(op,100) '       the l-axis is included by default.'
      WRITE(op,100) '1: k = 0.   2: h = 0.   3: h = k.   4: h = -k.'
      READ(cntrl,*,ERR=1) i_plane
      IF(cfile) WRITE(op,101) i_plane
      IF(i_plane < 1 .OR. i_plane > 4) THEN
        IF(cfile) THEN
          WRITE(op,100) 'Illegal choice. Default is k = 0.'
          i_plane = 1
        ELSE
          WRITE(op,100) 'Illegal choice. Choose again.'
          GO TO 1
        END IF
      END IF

! Get upper bounds. The SADP will be square, so we only need one bound.
! Choose l-axis, since this is common to all 4 views and makes scaling
! simpler if more than one SAD pattern is requested.
      2345 WRITE(op,100) 'Enter maximum value of l.'
      READ(cntrl,*,ERR=2345) l_upper
      IF(cfile) WRITE(op,104) l_upper

! 8-bit images or 16-bit images?
 3456 write(op,100) 'Choose the bit-depth of the image.'
      write(op,100) '   8: - 8-bits  16: - 16-bits    :  '
      read(cntrl,*,err=3456) bitdepth
      if(CFile) write(op,101) bitdepth
      if(bitdepth /= 8 .and. bitdepth /= 16) then
        write(op,100) 'Illegal bit-depth.'
        if(CFile) then
          write(op,100) 'Using 8-bit as the default'
          bitdepth = 8
        else
          write(op,100) 'Re-enter. . .  '
          goto 3456
        endif
      endif
      if(bitdepth == 16) then
! Bypass issue of signed or unsigned format. Use range 0 - 32767 only.
        write(op,100) 'File will be saved in unsigned 16-bit format.'
        maxsad = FIFTEENBITS - ONE
      else
        write(op,100) 'File will be saved in unsigned 8-bit format.'
        maxsad = EIGHTBITS - ONE
      endif

! Logarithmic or linear intensity scaling?
      2 WRITE(op,100) 'Choose intensity scaling'
      WRITE(op,100) '   0: - Logarithmic  1: - Linear'
      READ(cntrl,*,ERR=2) loglin
      IF(cfile) WRITE(op,101) loglin
      IF(loglin < 0 .OR. loglin > 1) THEN
        WRITE(op,100) 'Illegal intensity scaling type.'
        IF(cfile) THEN
          WRITE(op,100) 'Using linear as the default'
          loglin = 1
        ELSE
          WRITE(op,100) 'Re-enter. . . '
          GO TO 2
        END IF
      END IF

! By how much do we wish to saturate?
    3 write(op,100) 'Enter a brightness (must be +ve)'
      if(bitdepth == 16) then
        write(op,100)'  1 - scales intensities to the range 0 - 32767'
        write(op,100)' 10 - scales intensities to the range 0 - 327670,'
        write(op,100)'      but values above 65535 will be saturated'
      else
        write(op,100)'  1 - scales intensities to the range 0 - 255'
        write(op,100)' 10 - scales intensities to the range 0 - 2550,'
        write(op,100)'      but values above 255 will be saturated'
      endif
      READ(cntrl,*,ERR=3) brightness
      IF(cfile) WRITE(op,104) brightness
      IF(brightness <= zero) THEN
        WRITE(op,100) 'Illegal value for brightness. Must be positive'
        IF(cfile) THEN
          IF(loglin == 0) THEN
            WRITE(op,100) 'Using default of 1 for logarithmic scaling'
            brightness = one
          ELSE
            WRITE(op,100) 'Using default of 10 for linear scaling'
            brightness = ten
          END IF
        ELSE
          IF(loglin == 0) THEN
            WRITE(op,100) '1 is a good default for logarithmic scaling'
          ELSE
            WRITE(op,100) '10 is a good default for linear scaling'
          END IF
          WRITE(op,100) 'Re-enter. . . '
          GO TO 3
        END IF
      END IF

! Generate intensity data for the SAD pattern.
      IF(i_adapt == 1) THEN
        CALL getsad(aglq16, i_plane, l_upper, hk_lim, infile, ok)
      ELSE
        CALL getsad( glq16, i_plane, l_upper, hk_lim, infile, ok)
      END IF

      IF(ok) CALL wrtsad(outfile, i_plane, l_upper, hk_lim, ok)

      999 RETURN
      990 WRITE(op,102) 'Problems opening output file ', outfile
      WRITE(op,103) 'IOSTAT = ', io_err
      ok = .false.
      RETURN
      100 FORMAT(1X, a)
      101 FORMAT(1X, i3)
      102 FORMAT(1X, 2A)
      103 FORMAT(1X, a, i5)
      104 FORMAT(1X, g12.5)
      END SUBROUTINE gosadp

! ______________________________________________________________________
! Title: GOSPEC
! Author: MMJT
! Date: 23 Oct 1989
! Description: This subroutine sets up a file whose name is given
! by 'outfile', which is derived from the input filename 'infile'.
! The powder pattern data is then written to 'outfile', which is closed
! by the subroutine 'WRTSPC'.

!      ARGUMENTS:
!            infile   -  name of input file (input)
!            outfile  -  name of output file (output)
!            ok       -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  blurring, CFile, GAUSS, LORENZ, NONE, PS_VGT
!                   PV_GSS, PV_LRN, cntrl

!        modifies:  full_brd
! ______________________________________________________________________

      SUBROUTINE gospec(infile,outfile,ok)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',  infile, outfile ,ok , GETSPC, TRMSPC, RDRNGE ,
!                             io_err , AGLQ16, GLQ16, cut_off
!     Utiliza las funciones: RDRNGE, GETSPC, AGLQ16, GLQ16, TRMSPC      externas
!     Utiliza las subrutinas: GETFNM, PV, GAUSSN, LORNZN, WRTSPC

      CHARACTER (LEN=*), INTENT(IN OUT)        :: infile
      CHARACTER (LEN=*)                        :: outfile
      LOGICAL, INTENT(IN OUT)                  :: ok

      INTEGER*4 io_err, i
      REAL*8   cut_off

! get angular range and step
   !   ok = rdrnge()
   !   IF(.NOT.ok) GO TO 999

! create a new filename to save spectrum data to
   !   CALL getfnm(infile, outfile, 'spc', ok)

    !  IF(.NOT.ok) GO TO 999
    !  IF(unit_sp /= op) OPEN(UNIT = unit_sp, FILE = outfile, STATUS = 'new',  &
    !      ERR = 990, IOSTAT = io_err)
      full_brd = 1
     ! 3456 WRITE(op,100) 'Enter 1 for adaptive quadrature on broad peaks: '
     ! READ(cntrl,*,ERR=3456) full_brd
     ! IF(cfile) WRITE(op,101) full_brd

      if (conv_g == 1 .or. numcal == 0 ) then
          IF(full_brd == 1) THEN
           ok = getspc(aglq16, infile)
          ELSE
            ok = getspc(glq16, infile)
          END IF

      end if


! suppress the huge peak near the origin if required
          cut_off = zero
      IF(ok .AND. (th2_min == zero) .AND. trim_origin) ok = trmspc(cut_off)


      IF(ok) THEN
        IF(blurring == gauss) THEN
          !WRITE(op,104) 'Gaussian'
          CALL gaussn(cut_off)
        ELSE IF(blurring == lorenz) THEN
          !WRITE(op,104) 'Lorentzian'
          CALL lornzn(cut_off)
        ELSE IF(blurring == ps_vgt) THEN
          !WRITE(op,104) 'pseudo-Voigt'
          CALL pv(cut_off)
        ELSE IF(blurring == pv_gss) THEN
          !WRITE(op,104) 'Gaussian'
          CALL pv(cut_off)
        ELSE IF(blurring == pv_lrn) THEN
          !WRITE(op,104) 'Lorentzian'
          CALL pv(cut_off)
        ELSE IF(blurring /= NONE) THEN
          WRITE(op,100) 'Instrumental broadening type is undefined in GOSPEC.'
        END IF
      END IF
  !    IF(ok) CALL wrtspc(outfile, ok)


      999 RETURN
     ! 990 WRITE(op,102) 'Problems opening output file ', outfile
     ! WRITE(op,103) 'IOSTAT = ', io_err
     ! ok = .false.
     ! RETURN
      100 FORMAT(1X, a)
      101 FORMAT(1X, i3)
      102 FORMAT(1X, 2A)
      103 FORMAT(1X, a, i5)
      104 FORMAT(1X, '=> Adding ', a, ' instrumental broadening . . .')
      END SUBROUTINE gospec

! ______________________________________________________________________
! Title: GOSTRK
! Author: MMJT
! Date: 23 Oct 1989
! Description: This subroutine sets up a file whose name is given
! by 'outfile', which is derived from the input filename 'infile'.
! The streak intensity is then written to 'outfile', which is closed
! by the subroutine 'STREAK'.

!      ARGUMENTS:
!            infile   -  name of input file. (input).
!            outfile  -  name of output file. (output).
!            ok       -  logical flag indicating all went well.
!                                                       (output).

!      COMMON VARIABLES:
!            uses:  cntrl, CFile

!        modifies:  no COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE gostrk(infile,outfile,ok)
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',  infile, outfile,ok, i, io_err, AGLQ16, GLQ16
!     Utiliza las funciones:  AGLQ16, GLQ16     externas
!     Utiliza las subrutinas:  STREAK, GETFNM

      CHARACTER (LEN=*), INTENT(IN OUT)        :: infile
      CHARACTER (LEN=*), INTENT(OUT)           :: outfile
      LOGICAL, INTENT(IN OUT)                  :: ok

      INTEGER*4 i, io_err
! external functions

! external subroutines (Some compilers need them declared external)
!      external STREAK, GETFNM

      CALL getfnm(infile, outfile, 'str', ok)
      IF(.NOT.ok) GO TO 999
      IF(sk /= op) OPEN(UNIT = sk, FILE = outfile, STATUS = 'new',  &
          ERR = 990, IOSTAT = io_err)
      4567 WRITE(op,100) ' '
    !  WRITE(op,100) 'CALCULATING INTENSITY ALONG A STREAK. . .'
      WRITE(op,100) 'Enter 1 for adaptive quadrature: '
      READ(cntrl,*,ERR=4567) i
      IF(cfile) WRITE(op,101) i
! select integration function
      IF(i == 1) THEN
        CALL streak(aglq16, outfile, ok)
      ELSE
        CALL streak( glq16, outfile, ok)
      END IF

      999 RETURN
      990 WRITE(op,102) 'Problems opening output file ', outfile
      WRITE(op,103) 'IOSTAT = ', io_err
      ok = .false.
      RETURN
      100 FORMAT(1X, a)
      101 FORMAT(1X, i3)
      102 FORMAT(1X, 2A)
      103 FORMAT(1X, a, i5)
      END SUBROUTINE gostrk

!!S Hfiles
! ______________________________________________________________________
! Title: HKL_LIM
! Author: MMJT
! Date: 2 August 1989
! Obtains upper limits of h, k, and l given the global variable
! 'max_angle' (in radians).
! The limits are returned in the global variables h_bnd, k_bnd and
! l_bnd. HKL_LIM may need to decrease the value of lambda if h_bnd
! and k_bnd are too small to allow adequate off-axis symmetry testing.
! lambda is restored in OPTIMZ after symmetry testing.

!      ARGUMENTS:
!           No arguments are used. All data is in 'COMMON'.

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, lambda

!        modifies:  h_bnd, k_bnd, l_bnd, lambda
! ______________________________________________________________________

      SUBROUTINE hkl_lim()
!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',  HKLUP, h, k ,l, y ,
!     Utiliza las funciones:  HKLUP(y)
!     Utiliza las subrutinas:

      INTEGER*4 hklup
      REAL*8  y

! HKLUP returns the maximum value of h, k or l given 'max_angle'
      hklup(y) =  INT(two * SIN(half*max_angle) / (lambda*SQRT(y)))
! define upper h, k, l values consistent with max_angle
      1 h_bnd = hklup(a0)
      k_bnd = hklup(b0)
     l_bnd = hklup(c0)
!write(*,*) "hkl_lim2", lambda, h_bnd, k_bnd
! make sure bounds are not too small. This could occur for a small
! unit cell, a small value of max_angle, or a long wavelength.
      IF(h_bnd < 2 .OR. k_bnd < 2) THEN
        lambda = half * lambda
        GO TO 1
      END IF
! Make sure bounds are not too large either
      IF(h_bnd > 10) h_bnd = 10
      IF(k_bnd > 10) k_bnd = 10
      IF(l_bnd > 10) l_bnd = 10
      RETURN
      END SUBROUTINE hkl_lim

!!S Ifiles
! ______________________________________________________________________
! Title: INTEGR
! Author: MMJT and MWD
! Date: 15 Feb 1990; 7 Mar 1995; 28th May 1996
! Description:  This routine integrates intensity from
!               h,k,l0 to h,k,l1.

!      ARGUMENTS:
!            FN   -  Function name passed by reference. The
!                    choice is between GLQ16 (non-adaptive
!                    Gauss-Legendre integration), and AGLQ16
!                    (adaptive Gauss-Legendre integration). (input).
!            ok   -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:   a0, b0, c0, d0, lambda, cntrl, CFile, xplcit,
!                    X_RAY, rad_type, th2_max

!        modifies:   no COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE integr(fn, ok)


!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc', FN, ok, divided  ,h, k
!                            l0, l1, max_th, x, W4, ANGLE, theta, l, S, Q2, t1,
!                            sum, tmp, LL, fact, d_th, l_tmp
!     Utiliza las funciones:  FN   externa  S(h,k,l) , ANGLE(h,k,l)
!                             LL(theta,h,k), W4(theta)
!     Utiliza las subrutinas: XYPHSE, PRE_MAT, GET_F



      LOGICAL, INTENT(IN OUT)                  :: ok



      LOGICAL :: divided
      INTEGER*4 h, k, i_th, i_thm
      REAL*8 l0, l1, max_th, x, w4, angle, theta, l, s, q2
      REAL*8 t1, sum, tmp, ll, fact, d_th, l_tmp

! external function, passed by reference
      REAL*8 fn
      EXTERNAL fn
! external subroutines (Some compilers need them declared external)
!      external XYPHSE, PRE_MAT, GET_F

! statement functions
! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
! ANGLE is the Bragg angle (in radians) of the h,k,l plane
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))
! LL is the maximum l value for a given h,k
      ll(theta,h,k) = SQRT((q2 * SIN(theta) * SIN(theta)  &
          - h*h*a0 - k*k*b0 - h*k*d0)/ c0)
! W4 is the X-ray polarization factor
      w4(theta) = half * (one + (COS(two*theta))**2)

      max_th = half * th2_max
      fact = two / lambda
      q2 = fact * fact
! set increment to a safe value, in case l1-l0 is too large
      d_th = eps3 * deg2rad

      10 WRITE(op,400) 'Enter h, k, l0, l1: '
      READ(cntrl,*,ERR=10)  h, k, l0, l1
      IF(cfile) WRITE(op,401) h, k, l0, l1
! check values
      IF(l1 == l0) THEN
        WRITE(op,400) 'Illegal input: l1 equals l0'
        IF(cfile) THEN
          l1 = l0 + d_th
          WRITE(op,403) 'l1 is set to ', l1
        ELSE
          GO TO 10
        END IF
      END IF
! make sure we are not going to blow up at the origin
      IF(h == 0 .AND. k == 0 .AND. rad_type == electn) THEN
        IF(l0*l1 <= zero) THEN
          WRITE(op,400) 'Cannot integrate across the origin for electron radiation'
          WRITE(op,400) 'Re-enter. . . '
          GO TO 10
        END IF
      END IF
! Finally, check angles are meaningful
      IF(s(h,k,l0) > q2 .OR. s(h,k,l1) > q2) THEN
        IF(s(h,k,l0) > q2) WRITE(op,402) h, k, l0,  &
            ' exceeds 180 degree scattering angle!'
        IF(s(h,k,l1) > q2) WRITE(op,402) h, k, l1,  &
            ' exceeds 180 degree scattering angle!'
        GO TO 10
      END IF
! get angles corresponding to start and stop
! check if we need to break the integration region into two parts
! because h,k,l, and h,k,-l may subtend the same angle
      divided  = .false.
      IF(l0 <= zero .AND. l1 <= zero) THEN
! use Friedel's law, and keep l +ve. Swap l0 and l1
        h = -h
        k = -k
        tmp = -l0
        l0 = -l1
        l1 = tmp
      ELSE IF(l0 < zero .AND. l1 > zero) THEN
        h = -h
        k = -k
        l_tmp = l1
        l1 = -l0
        l0 = zero
        divided = .true.
      END IF
! swap if in reverse order
      IF(l0 > l1) THEN
        tmp = l0
        l0 = l1
        l1 = tmp
      END IF

      sum = zero
      30 max_th = angle(h,k,l1)
      t1     = angle(h,k,l0)
      l1 = l0

      CALL xyphse(h, k)

      CALL pre_mat(h, k)
! integrate each d_th's range of reciprocal space

      !DO  theta = t1, max_th-eps14, d_th
      i_thm=nint((max_th-eps14-t1)/d_th+1.0d0)
      DO  i_th = 1, i_thm
        theta=t1+(i_th-1)*d_th
        l0 = l1
        tmp = MIN(d_th, max_th-theta)
        l1 = ll(theta+tmp,h,k)
        x = fn(h,k,l0,l1,ok)
        IF(.NOT. ok) GO TO 999
        IF(rad_type == x_ray) x = x * w4(theta + half*tmp)
        sum = sum + x
      END DO
! do we need to integrate the other part?
      IF(divided) THEN
! goto -h,-k,-l and continue
        h = -h
        k = -k
        l0 = zero
        l1 = l_tmp
        divided = .false.
        GO TO 30
      END IF
      WRITE(op,403) 'Integrated intensity = ', sum
      999 RETURN
      400 FORMAT(1X, a)
      401 FORMAT(1X, 2I3, 2G12.5)
      402 FORMAT(1X, 2I3, g12.5, a)
      403 FORMAT(1X, a, g15.8)
      END SUBROUTINE integr

! ______________________________________________________________________
! Title: INTEN2
! Author: MMJT
! Date: 8 Apr 1992
! Description: This routine determines the intensity at (h,k,l)
! in reciprocal space, for an explicitly defined stacking sequence.

!      ARGUMENTS:
!            f   -  Layer form factors. (input).
!            h   -  reciprocal lattice vector h-component. (input).
!            k   -  reciprocal lattice vector k-component. (input).
!            l   -  reciprocal lattice vector l-component. (input).
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  same_layer, l_seq, l_r, n_layers, l_phi, l_cnt
!                   PI2

!        modifies:  wavefn

!      INTEN2 returns the intensity at h, k, l
! ______________________________________________________________________

      REAL*8 FUNCTION inten2(f, h, k, l, ok)



!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc',ok=.true.  ,  h, k , l, f(MAX_L)
!                             i, j, m , twopi_l, dot, tmp, phi(MAX_L, MAX_L), z, z_to_n
!     Utiliza las funciones:
!     Utiliza las subrutinas:


      COMPLEX*16, INTENT(IN OUT)               :: f(max_l)
      INTEGER*4, INTENT(IN OUT)                :: h
      INTEGER*4, INTENT(IN)                    :: k
      REAL*8, INTENT(IN OUT)                   :: l
      LOGICAL, INTENT(IN OUT)                  :: ok





      INTEGER*4 i, j, m
      REAL*8 twopi_l, dot, tmp
      COMPLEX*16 phi(max_l, max_l), z, z_to_n

      twopi_l = pi2 * l

! 'ok' is not used, but is included for compatibility with INTENS
      ok = .true.

! Is there only one layer? If so, let's get it over with.
      IF(nint(l_cnt) == 1) THEN
        wavefn = f(l_seq(1))
        GO TO 900
      END IF

! Check for an obvious optimization when all layers are identical.
      IF(same_layer) THEN
        i = l_seq(1)
        dot = pi2*(h*l_r(1,i,i) + k*l_r(2,i,i) + l*l_r(3,i,i))
        z = DCMPLX( COS(dot), SIN(dot) )
        tmp = dot / pi2
! check we are not about to execute 0.0 / 0.0
        IF(ABS(tmp - nint(tmp)) <= eps5) THEN
          wavefn = f(i) * l_cnt
        ELSE
! sum the series
          dot = dot * l_cnt
          z_to_n = DCMPLX( COS(dot), SIN(dot) )
          wavefn = f(i) * (c_one - z_to_n) / (c_one - z)
        END IF
        GO TO 900
      END IF

! Else do it the long way
! Get phases
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          dot = twopi_l * l_r(3,j,i)
          phi(j,i) = l_phi(j,i) * DCMPLX( COS(dot), SIN(dot) )
        END DO
      END DO
! Count down to the first layer (we know l_cnt is greater than 1)
! Initialize wavefunction to the scattering factor of the last layer
      wavefn = f(l_seq(nint(l_cnt)))

      DO  m = nint(l_cnt) - 1, 1, -1
        i = l_seq(m)
        j = l_seq(m+1)
        wavefn = f(i) + wavefn * phi(j,i)
      END DO

! Normalize to the number of layers
      900 inten2 = wavefn * CONJG(wavefn) / l_cnt

      RETURN
      END FUNCTION inten2

! ______________________________________________________________________
! Title: INTENS
! Author: MWD and MMJT
! Date: 10 Feb 1989
! Description: This routine determines the intensity at (h,k,l)
! in reciprocal space, for recursive stacking. For this function
! to be called, 'rcrsv' must be TRUE.
! Note: The diffuse background is handled automatically by the
! recursion algorithm.

!      ARGUMENTS:
!            f   -  Layer form factors. (input).
!            h   -  reciprocal lattice vector h-component. (input).
!            k   -  reciprocal lattice vector k-component. (input).
!            l   -  reciprocal lattice vector l-component. (input).
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  only_real, l_g, n_layers, inf_thick

!        modifies:  no COMMON variables are modified

!      INTENS returns the intensity at h, k, l
! ______________________________________________________________________

      REAL*8 FUNCTION intens(f, h, k, l, ok)


!     Utiliza las variables:  DIFFaX.par', DIFFaX.inc' ,  ok ,  h, k , l, f(MAX_L) ,  GET_S, GET_S2
!                             i , sum, x , s(MAX_L)
!     Utiliza las funciones: GET_S, GET_S2         externas
!     Utiliza las subrutinas: GET_MAT

      Complex(Kind=8), Intent(In) :: f(max_l)
      Integer,         Intent(In) :: h
      Integer,         Intent(In) :: k
      Real(Kind=8)                :: l
      Logical,         Intent(Out):: ok



      INTEGER*4 i
      REAL*8 sum, x
      COMPLEX*16 s(max_l)

! external function

! external subroutines (Some compilers need them declared external)
!      external GET_MAT
      sum = zero
      CALL get_mat(h, k, l)
      IF(inf_thick) THEN

! initialize s to -f, since mat is -(1 - T)
        DO  i = 1, n_layers
          s(i) = - f(i)
        END DO


        ok = get_s(f, s, h, k, l)
      ELSE
! s is initialized inside GET_S2, from where GET_S is called
        ok = get_s2(f, s, h, k, l)
      END IF
      IF(ok) THEN
! only use real part of f(i) if that's all that's there
        IF(only_real) THEN
          DO  i = 1, n_layers
            sum = sum + l_g(i) * DBLE(f(i)) * DBLE(s(i))
          END DO
          sum = two * sum
          DO  i = 1, n_layers
            x = DBLE(f(i))
            sum = sum - l_g(i) * x*x
          END DO
        ELSE
! must use complex part of f(i)
          DO  i = 1, n_layers
            sum = sum + l_g(i) * DBLE(CONJG(f(i)) * s(i))
          END DO
          sum = two * sum
          DO  i = 1, n_layers
            x = ABS(f(i))
            sum = sum - l_g(i) * x*x
          END DO
        END IF
      END IF

      intens = sum

      RETURN
      END FUNCTION intens

!!S Lfiles
! ______________________________________________________________________
! Title: LORNZN
! Author: MMJT
! Date: 17 Feb 1990; 7 Mar 1995
! Description: This subroutine performs the Lorentzian
! instrumental broadening. FWHM is in degrees. Does not conserve
! intensity well when FWHM is comparable to d_theta. Data at the
! extreme ends of the spectrum are corrupted slightly.

!      ARGUMENTS:
!            th2_low  -  lowest 2theta angle to consider. (input).

!      COMMON VARIABLES:
!            uses:  th2_max, d_theta, FWHM, spec, NONE, PI, RAD2DEG
!                   blurring

!        modifies:  brd_spc
! ______________________________________________________________________

      SUBROUTINE lornzn(th2_low)


!     Utiliza las variables: DIFFaX.par' , 'DIFFaX.inc' , th2_low , i, j, n_low, n_high
!                             k1, k2, k3, const, lrnz, tmp, tmp1, tmp2
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      REAL*8, INTENT(IN OUT)                       :: th2_low


      INTEGER*4 i, j, n_low
      REAL*8 k1, k2, k3, const, lrnz, tmp, tmp1, tmp2

      IF(fwhm <= zero) GO TO 999
! check that cut-off is reasonable
      IF(th2_low < zero .OR. th2_low >= th2_max) THEN
        WRITE(op,101) 'LORNZN: Cut-off angle ', th2_low,  &
            ' is out of bounds. Angle reset to zero.'
        th2_low = zero
      END IF

! th2_low is the angle relative to th2_min
! 2*d_theta is the angular step size
      n_low  = INT(half*th2_low/d_theta) + 1

      const = two * rad2deg * d_theta
      k1 = const * two / (pi * fwhm)
      k2 = (const * two / fwhm)**2

      DO  i = 1, n_high
        brd_spc(i) = zero
      END DO

      DO  i = 0, n_high
        k3 = one + k2*i*i
        lrnz = k1 / k3
        DO  j = n_low+1, n_high
          tmp1 = zero
          tmp2 = zero
          IF((j-i) > n_low)  tmp1 = spec(j-i)
          IF((j+i) <= n_high) tmp2 = spec(j+i)
          tmp = tmp1 + tmp2
          IF(i == 0) tmp = half * tmp
          brd_spc(j) = brd_spc(j) + lrnz * tmp
        END DO
      END DO
      RETURN
      999 WRITE(op,101) 'Illegal FWHM ', fwhm, ' in LORNZN()'
      WRITE(op,100) 'Lorentzian instrumental broadening not added'
! kill blurring option
      blurring = NONE
      RETURN
      100 FORMAT(1X, a)
      101 FORMAT(1X, a, g12.5, a)
      END SUBROUTINE lornzn

! ______________________________________________________________________
! Title: LUBKSB
! Author: MWD, adapted from
! "Numerical Recipes: The Art of Scientific Computing."
! Date: 18 Aug 1988
!  Description: Solves the set of linear equations ax = b,
!  where a, x and b contain real variables.
!  Here, a is input as the LU-decomposition of a, determined by the
!  routine LUDCMP. index is input as the permutation vector returned
!  by LUDCMP. b is input as the right hand side vector, and returns
!  with the solution x. a, n, MAX_N and index are not modified by this
!  routine. In DIFFaX, LUDCMP and LUBKSB are used to solve for l_g(i),
!  the a-priori probability that layer i exists within the crystal.

!      ARGUMENTS:
!            a      -  LU-decomposed square matrix of real numbers.
!                                                           (input).
!            b      -  vector of real numbers, the right hand side of
!                      ax = b, is input. The solution x is output.
!            index  -  vector containing the record of the row
!                      permutation effected by the partial pivoting
!                      carried out in CLUDCM (input).
!            n      -  size of the square matrix. (input).
!            MAX_N  -  physical dimension of a (MAX_N x MAX_N).(input).
! ______________________________________________________________________

      SUBROUTINE lubksb(a,b,INDEX,n,max_n)
!     save


!     Utiliza las variables: DIFFaX.par' , n, MAX_N, index(MAX_N), a(MAX_N,MAX_N), b(MAX_N) , i, i2, j, row ,sum
!     Utiliza las funciones:
!     Utiliza las subrutinas:


      REAL*8, INTENT(IN)                       :: a(max_n,max_n)
      REAL*8, INTENT(OUT)                      :: b(max_n)
      INTEGER*4, INTENT(IN)                    :: INDEX(max_n)
      INTEGER*4, INTENT(IN)                    :: n
      INTEGER*4, INTENT(IN)                    :: max_n



      INTEGER*4 i, i2, j, row
      REAL*8 sum


      i2 = 0
      DO  i = 1, n
        row = INDEX(i)
        sum = b(row)
        b(row) = b(i)
        IF(i2 /= 0) THEN
          DO  j = i2, i-1
            sum = sum - a(i,j) * b(j)
          END DO
        ELSE IF(ABS(sum) /= zero) THEN
          i2 = i
        END IF
        b(i) = sum
      END DO
      DO  i = n, 1, -1
        sum = b(i)
        DO  j = i+1, n
          sum = sum - a(i,j) * b(j)
        END DO
        b(i) = sum / a(i,i)
      END DO
      RETURN
      END SUBROUTINE lubksb

! ______________________________________________________________________
! Title: LUDCMP
! Author: MWD, adapted from
! "Numerical Recipes: The Art of Scientific Computing."
! Date: 18 Aug 1988
!  Description: This is an LU decomposition routine, and accepts
!  real*8 variables.
!  Given an n x n matrix a, with physical dimension MAX_N, this
!  routine replaces it by the LU decomposition of a rowwise permutation
!  of itself. a and n are input. a is the LU decomposed output; index
!  is an output vector which records the row permutation affected by
!  the partial pivoting; Det is the determinant of a. This routine is
!  used in combination with LUBKSB to solve linear equations. In DIFFaX,
!  these routines are used to solve for l_g(i), the a-priori
!  probability that layer i exists within the crystal.
!  LUDCMP returns .false. if the matrix turns out to be singular.

!      ARGUMENTS:
!            a      -  Square matrix of real numbers to LU-decompose is
!                      input. a is then replaced by the
!                      LU-decomposed result.
!            index  -  output vector which records the row permutation
!                      affected by the partial pivoting. (output).
!            n      -  size of the square matrix. (input).
!            MAX_N  -  physical dimension of a (MAX_N x MAX_N).(input).
!            Det    -  determinant of the matrix. (output).

!      LUDCMP returns logical .true. if all proceeded happily.
! ______________________________________________________________________

      LOGICAL FUNCTION ludcmp(a,INDEX,n,max_n,det)
!     save

!     Utiliza las variables: DIFFaX.par' , n, MAX_N, index(MAX_N), a(MAX_N,MAX_N),Det
!                            L_MAX=100, i, j, m, row , tiny = 1.0D-20, tmp(L_MAX), sum, max, tmp2
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      REAL*8, INTENT(IN OUT)                   :: a(max_n,max_n)
      INTEGER*4, INTENT(IN OUT)                :: INDEX(max_n)
      INTEGER*4, INTENT(IN OUT)                :: n
      INTEGER*4                                :: max_n
      REAL*8, INTENT(IN OUT)                   :: det



      INTEGER*4  i, j, m, row
      INTEGER*4, PARAMETER :: l_max = 100
      REAL*8  tmp(l_max), sum, MAX, tmp2
      REAL*8, PARAMETER :: tiny = 1.0D-20

      ludcmp = .false.
      det = one
      IF(n > l_max) THEN
        WRITE(op,400) 'Matrix too large for LUDCMP'
        RETURN
      END IF
      DO  i = 1, n
        MAX = zero
        DO  j = 1, n
          IF(ABS(a(i,j)) > MAX) MAX = ABS(a(i,j))
        END DO
        IF(MAX == zero) CYCLE
        tmp(i) = one / MAX
      END DO
      DO  j = 1, n
        DO  i = 1, j-1
          sum = a(i,j)
          DO  m = 1, i-1
            sum = sum - a(i,m) * a(m,j)
          END DO
          a(i,j) = sum
        END DO
        MAX = zero
        DO  i = j, n
          sum = a(i,j)
          DO  m = 1, j-1
            sum = sum - a(i,m) * a(m,j)
          END DO
          a(i,j) = sum
          tmp2 = tmp(i)*ABS(sum)
          IF(ABS(tmp2) >= MAX) THEN
            row = i
            MAX = tmp2
          END IF
        END DO
        IF(j /= row) THEN
          DO  m = 1, n
            tmp2 = a(row,m)
            a(row,m) = a(j,m)
            a(j,m) = tmp2
          END DO
          det = -det
          tmp(row) = tmp(j)
        END IF
        INDEX(j) = row
        IF(ABS(a(j,j)) == zero) a(j,j) = tiny
        tmp2 = one / a(j,j)
        DO  i = j+1, n
          a(i,j) = a(i,j) * tmp2
        END DO
        det = det * a(j,j)
      END DO
      ludcmp = .true.
      RETURN
      200 CONTINUE
      RETURN
      400 FORMAT(1X, a)
      END FUNCTION ludcmp

! ______________________________________________________________________
! Title: L_STEP
! Author: MMJT
! Date: 13 Aug 1989
! Description:  L_STEP attempts to determine whether or not there
! are any sharp peaks, and if so, their spacing. The algorithm inspects
! the sum of certain permutations of the stacking vector z_components.
! The components are chosen so that their values should be independent
! of the layer origins chosen by the user. ie. Rz(i,i) and
! Rz(i,j) + Rz(j,i) and Rz(i,j) + Rz(j,k) + Rz(k,i), etc...
! In the interests of clarity, the permutations of the indices are
! written out explicity, instead of being generated by a subroutine
! (which would considerably shorten the code). We stop searching after
! combinations of 4 stacking vectors, since then the structure is too
! complicated to be attempting to trick DIFFaX into believing that we
! have found sharp peaks. Under these circumstances, work out the
! diffraction pattern the long, but reliable, way.

!      ARGUMENTS:
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  n_layers, there, l_r

!        modifies:  any_sharp

!      L_STEP returns the l-increment that sharp peaks are likely to
!      be found at.
! ______________________________________________________________________

      REAL*8 FUNCTION l_step(ok)



!     Utiliza las variables: DIFFaX.par' , 'DIFFaX.inc', ok   , tmp, z_step , YRDSTK, resonant, decided
!                              i1, i2, i3, i4
!     Utiliza las funciones: YRDSTK externa
!     Utiliza las subrutinas:


      LOGICAL, INTENT(IN OUT)                  :: ok


      REAL*8 tmp, z_step
      LOGICAL ::  resonant, decided
      INTEGER*4 i1, i2, i3, i4

! external function


! initialize return value
      l_step = zero

      resonant = .true.
      decided  = .false.
      z_step   = zero

! Check z-components of Rii stacking vectors
! if any of the transitions do not exist, set resonant to false.
      DO  i1 = 1, n_layers
        IF(resonant) THEN
          IF(there(i1,i1)) THEN
            decided = .true.
            tmp = l_r(3,i1,i1)
            resonant = resonant .AND. yrdstk(z_step, tmp, ok)
            IF(.NOT. ok) GO TO 990
          END IF
        END IF
      END DO

! Rii terms do not occur (ie. no layer i will stack to another layer i)
! We must therefore check z-components of Rij + Rji sequences (i.ne.j).
      IF((n_layers > 1) .AND. .NOT.decided) THEN
        DO  i1 = 1, n_layers
          DO  i2 = i1 + 1, n_layers
            IF(resonant) THEN
              IF(there(i2,i1) .AND. there(i1,i2)) THEN
                decided = .true.
                tmp = l_r(3,i2,i1) + l_r(3,i1,i2)
                resonant = resonant .AND. yrdstk(z_step, tmp, ok)
                IF(.NOT. ok) GO TO 991
              END IF
            END IF
          END DO
        END DO
      END IF

! No Rij + Rji sequences occur.
! Check z-components of Rij + Rjk + Rki sequences (where i.ne.j.ne.k).
      IF((n_layers > 2) .AND. .NOT. decided) THEN
        DO  i1 = 1, n_layers
          DO  i2 = i1 + 1, n_layers
            DO  i3 = i2 + 1, n_layers
              IF(resonant) THEN
! There are 2 permutations
                IF(there(i2,i1).AND. there(i3,i2).AND.  &
                      there(i1,i3)) THEN
                  decided = .true.
                  tmp = l_r(3,i2,i1) + l_r(3,i3,i2) + l_r(3,i1,i3)
                  resonant = resonant .AND. yrdstk(z_step, tmp, ok)
                  IF(.NOT. ok) GO TO 992
                END IF
                IF(there(i3,i1).AND. there(i2,i3).AND.  &
                      there(i1,i2) .AND. resonant) THEN
                  decided = .true.
                  tmp = l_r(3,i3,i1) + l_r(3,i2,i3) + l_r(3,i1,i2)
                  resonant = resonant .AND. yrdstk(z_step, tmp, ok)
                  IF(.NOT. ok) GO TO 993
                END IF
              END IF
            END DO
          END DO
        END DO
      END IF

! No Rij + Rjk + Rki sequences occur.
! Check z-components of Rij + Rjk + Rkl + Rli sequences
! (where i.ne.j.ne.k.ne.l).
      IF((n_layers > 3) .AND. .NOT.decided) THEN
        DO  i1 = 1, n_layers
          DO  i2 = i1 + 1, n_layers
            DO  i3 = i2 + 1, n_layers
              DO  i4 = i3 + 1, n_layers
                IF(resonant) THEN
! There are 6 permutations
                  IF(there(i2,i1).AND. there(i3,i2).AND.  &
                        there(i4,i3).AND. there(i1,i4)) THEN
                    decided = .true.
                    tmp = l_r(3,i2,i1) + l_r(3,i3,i2) + l_r(3,i4,i3) + l_r(3,i1,i4)
                    resonant = resonant.AND.yrdstk(z_step,tmp,ok)
                    IF(.NOT. ok) GO TO 994
                  END IF
                  IF(there(i2,i1).AND. there(i4,i2).AND.  &
                        there(i3,i4).AND. there(i1,i3) .AND. resonant) THEN
                    decided = .true.
                    tmp = l_r(3,i2,i1) + l_r(3,i4,i2) + l_r(3,i3,i4) + l_r(3,i1,i3)
                    resonant = resonant .AND. yrdstk(z_step,tmp,ok)
                    IF(.NOT. ok) GO TO 995
                  END IF
                  IF(there(i3,i1).AND. there(i2,i3).AND.  &
                        there(i4,i2).AND. there(i1,i4) .AND. resonant) THEN
                    decided = .true.
                    tmp = l_r(3,i3,i1) + l_r(3,i2,i3) + l_r(3,i4,i2) + l_r(3,i1,i4)
                    resonant = resonant .AND. yrdstk(z_step,tmp,ok)
                    IF(.NOT. ok) GO TO 996
                  END IF
                  IF(there(i3,i1).AND. there(i4,i3).AND.  &
                        there(i2,i4).AND. there(i1,i2) .AND. resonant) THEN
                    decided = .true.
                    tmp = l_r(3,i3,i1) + l_r(3,i4,i3) + l_r(3,i2,i4) + l_r(3,i1,i2)
                    resonant = resonant .AND. yrdstk(z_step,tmp,ok)
                    IF(.NOT. ok) GO TO 997
                  END IF
                  IF(there(i4,i1).AND. there(i2,i4).AND.  &
                        there(i3,i2).AND. there(i1,i3) .AND. resonant) THEN
                    decided = .true.
                    tmp = l_r(3,i4,i1) + l_r(3,i2,i4) + l_r(3,i3,i2) + l_r(3,i1,i3)
                    resonant = resonant .AND. yrdstk(z_step,tmp,ok)
                    IF(.NOT. ok) GO TO 998
                  END IF
                  IF(there(i4,i1).AND. there(i3,i4).AND.  &
                        there(i2,i3).AND. there(i1,i2).AND.resonant) THEN
                    decided = .true.
                    tmp = l_r(3,i4,i1) + l_r(3,i3,i4) + l_r(3,i2,i3) + l_r(3,i1,i2)
                    resonant = resonant .AND. yrdstk(z_step,tmp,ok)
                    IF(.NOT.ok) GO TO 999
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
      END IF

! If there is no stacking sequence that can bring us back to a layer
! similar to that at the origin after 4 layers, then this
! structure is sufficiently complicated that we may be better
! off doing adaptive integration anyway. (d_l still equals 0.0.)

      IF(decided .AND. resonant .AND. (tmp /= zero)) THEN
        l_step = one / tmp
        any_sharp = .true.
      ELSE
        l_step = zero
        any_sharp = .false.
      END IF

      RETURN
      990 WRITE(op,200)
      WRITE(op,201) i1, i1
      RETURN
      991 WRITE(op,202)
      WRITE(op,203) i1, i2, i2, i1
      RETURN
      992 WRITE(op,202)
      WRITE(op,204) i1, i2, i2, i3, i3, i1
      RETURN
      993 WRITE(op,202)
      WRITE(op,204) i1, i3, i3, i2, i2, i1
      RETURN
      994 WRITE(op,202)
      WRITE(op,205) i1, i2, i2, i3, i3, i4, i4, i1
      RETURN
      995 WRITE(op,202)
      WRITE(op,205) i1, i2, i2, i4, i4, i3, i3, i1
      RETURN
      996 WRITE(op,202)
      WRITE(op,205) i1, i3, i3, i2, i2, i4, i4, i1
      RETURN
      997 WRITE(op,202)
      WRITE(op,205) i1, i3, i3, i4, i4, i2, i2, i1
      RETURN
      998 WRITE(op,202)
      WRITE(op,205) i1, i4, i4, i2, i2, i3, i3, i1
      RETURN
      999 WRITE(op,202)
      WRITE(op,205) i1, i4, i4, i3, i3, i2, i2, i1
      RETURN
      200 FORMAT(1X,'L_STEP: Non-physical z-component of stacking vector.')
      201 FORMAT(1X,'Rz(',i2,',',i2,') = 0.0')
      202 FORMAT(1X,'L_STEP:Non-physical z-components of stacking vectors')
      203 FORMAT(1X,'Rz(',i2,',',i2,') + Rz(',i2,',',i2,') = 0.0')
      204 FORMAT(1X,'Rz(',i2,',',i2,')',2(' + Rz(',i2,',',i2,')'),' = 0.0')
      205 FORMAT(1X,'Rz(',i2,',',i2,')',3(' + Rz(',i2,',',i2,')'),' = 0.0')
      END FUNCTION l_step

!!S Mfiles
! ______________________________________________________________________
! Title: MATMUL
! Author: MMJT
! Date: 5 Feb 1990
! Description:  This subroutine multiplies the complex matrices
! 'a' and 'b', of logical size n x n. Result is returned in 'a'.

!      ARGUMENTS:
!            a   -  Complex*16 array to store result. (input and output).
!            b   -  Complex*16 array. (input).
!            n   -  Logical size of matrices. (input).

! ______________________________________________________________________

      SUBROUTINE matmul(a, b, n)
!     save

!     Utiliza las variables: DIFFaX.par' , n  , a(MAX_L,MAX_L), b(MAX_L,MAX_L),  i, j, m  , c(MAX_L,MAX_L), ctmp
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      COMPLEX*16, INTENT(IN OUT)               :: a(max_l,max_l)
      COMPLEX*16, INTENT(IN)                   :: b(max_l,max_l)
      INTEGER*4, INTENT(IN OUT)                :: n



      INTEGER*4 i, j, m
      COMPLEX*16 c(max_l,max_l), ctmp

! first copy a into c
      DO  j = 1, n
        DO  i = 1, n
          c(i,j) = a(i,j)
        END DO
      END DO

      DO  j = 1, n
        DO  i = 1, n
          ctmp = c_zero
          DO  m = 1, n
            ctmp = ctmp + c(i,m) * b(m,j)
          END DO
          a(i,j) = ctmp
        END DO
      END DO

      RETURN
      END SUBROUTINE matmul

! ______________________________________________________________________
! Title: MAT2N
! Author: MMJT
! Date: 5 Feb 1990
! Description:  This function multiplies the matrix 'mat' by itself
! n = l_cnt+1 times. In order to speed up the process, n has been broken
! down into its binary representation (by the function BINPOW, which is
! called once in OPTIMZ), and the result is given by

!         mat**n = mat**n0 * mat**n1 * mat**n2 * mat**n3 * etc

! where ni = 2**i
! and n = n0 + n1 + n2 + n3 + . . .

! mat**ni is given by (mat**n(i-1)) * (mat**n(i-1))
! ie. mat**8 = (mat**4) * (mat**4)
! and similarly mat**4 = (mat**2) * (mat**2)

! n must be such that n <= RCSV_MAX+1 <= 2**(MAX_BIN+1) - 1

!      ARGUMENTS:
!            a   -  Complex*16 array to store result. (output).

!      COMMON VARIABLES:
!            uses:  n_layers, pow, max_pow

!        modifies:  No COMMON variable are modified.

!      MAT2N returns TRUE if all went well.
! ______________________________________________________________________

      LOGICAL FUNCTION mat2n(a)


!     Utiliza las variables: DIFFaX.par' ,'DIFFaX.inc', a(MAX_L,MAX_L),i, j ,tmp_mat(MAX_L,MAX_L,MAX_BIN)
!     Utiliza las funciones:
!     Utiliza las subrutinas:  MATSQR, MATMUL

      COMPLEX*16, INTENT(IN OUT)               :: a(max_l,max_l)


      INTEGER*4 i, j
      COMPLEX*16 tmp_mat(max_l,max_l,max_bin)

! external subroutine (Some compilers need them declared external)
!      external MATSQR, MATMUL

      mat2n = .false.

! copy mat into the first 2-dimensional tmp_mat array. Initialize a
! to be the identity matrix.
      DO  j = 1, n_layers
        DO  i = 1, n_layers
          a(i,j) = c_zero
          tmp_mat(i,j,1) = mat(i,j)
        END DO
        a(j,j) = c_one
      END DO

      DO  i = 1, max_pow-1
        IF(pow(i) == 1) THEN
          CALL matmul(a, tmp_mat(1,1,i), n_layers)
        END IF
        CALL matsqr(tmp_mat(1,1,i+1), tmp_mat(1,1,i), n_layers)
      END DO
      IF(pow(max_pow) == 1) CALL matmul(a, tmp_mat(1,1,max_pow), n_layers)

      mat2n = .true.

      RETURN
      END FUNCTION mat2n

! ______________________________________________________________________
! Title: MATSQR
! Author: MMJT
! Date: 5 Feb 1990
! Description:  This subroutine multiplies the complex matrix
! 'b', of logical size n x n, by itself. Result is returned in 'a'.

!      ARGUMENTS:
!            a   -  Complex*16 array to store result. (output).
!            b   -  Complex*16 array to be 'squared'. (input).
!            n   -  Logical size of matrices. (input).

! ______________________________________________________________________

      SUBROUTINE matsqr(a, b, n)
!     save


!     Utiliza las variables: DIFFaX.par' ,  n , a(MAX_L,MAX_L), b(MAX_L,MAX_L) , i, j, m ,ctmp
!     Utiliza las funciones:
!     Utiliza las subrutinas:


      COMPLEX*16, INTENT(IN OUT)               :: a(max_l,max_l)
      COMPLEX*16, INTENT(IN)                   :: b(max_l,max_l)
      INTEGER*4, INTENT(IN OUT)                :: n



      INTEGER*4 i, j, m
      COMPLEX*16 ctmp

      DO  j = 1, n
        DO  i = 1, n
          ctmp = c_zero
          DO  m = 1, n
            ctmp = ctmp + b(i,m) * b(m,j)
          END DO
          a(i,j) = ctmp
        END DO
      END DO

      RETURN
      END SUBROUTINE matsqr

!!S Nfiles
! ______________________________________________________________________
! Title: NAMER
! Author: MMJT
! Date: 13 Oct 1988
! Description: This subroutine creates the appropriate filenames for
! the various files needed.
! It generates the filename by reading the structure data file name.
! It scans the name to see if a period ('.') is already within
! the name (ie. if the data file is called 'Structure.dat'). It deletes
! the characters after the period and appends whatever is in append.
! 'name1' and 'name2' are assumed to be the same physical length.
! This subroutine is called by GETFNM

!      ARGUMENTS:
!            name1   -  The name of the input data file. (input).
!            name2   -  A derivative filename. (output).
!            append  -  A token to be appended to name2. (input).
! ______________________________________________________________________

      SUBROUTINE namer(name1,name2,APPEND)
!     save

!     Utiliza las variables: DIFFaX.par' ,  name1, name2, append ,applen, namlen, i, idot ,LENGTH
!     Utiliza las funciones: LENGTH     externa
!     Utiliza las subrutinas:

      CHARACTER (LEN=*), INTENT(IN)        :: name1
      CHARACTER (LEN=*), INTENT(IN OUT)    :: name2
      CHARACTER (LEN=*), INTENT(IN OUT)    :: APPEND


      INTEGER*4 applen, namlen, i, idot


! get length of the string holding the filename
      namlen = LEN(name1)
      name2  = ' '
      applen = length(APPEND)

      idot = INDEX(name1,'.')

      IF(idot == 0) idot = INDEX(name1,' ')

      IF(idot == 0) THEN
        IF(namlen < max_nam) THEN
          idot = namlen + 1
        ELSE
          idot = max_nam
          namlen = max_nam
        END IF
      END IF
! truncate root filename so that the appendage will appear correctly
     ! IF(idot+applen > namlen) idot = namlen - applen

      DO  i = 1, idot-1
        WRITE(name2(i:i),'(a)') name1(i:i)
      END DO

      WRITE(name2(idot:idot),'(a)') '.'

      DO  i = 1, applen
        WRITE(name2((idot+i):(idot+i)),'(a)') APPEND(i:i)
      END DO


      DO  i = idot+applen+1, namlen
        WRITE(name2(i:i),'(a)') ' '
      END DO

      RETURN
      END SUBROUTINE namer

! ______________________________________________________________________
! Title: NMCOOR
! Author: MWD
! Date: 18 Aug 1988
! Description:  This subroutine multiplies the relative coordinates
! of each atom by 2*pi, which is the useful form for calculating phases.

!      ARGUMENTS:
!            No input arguments.

!      COMMON VARIABLES:
!            uses:  n_actual, l_n_atoms, a_pos, PI2

!        modifies:  a_pos
! ______________________________________________________________________

      SUBROUTINE nmcoor()


!     Utiliza las variables: DIFFaX.par' , 'DIFFaX.inc'  , i, j
!     Utiliza las funciones:
!     Utiliza las subrutinas:


      INTEGER*4 i, j

      DO  i = 1, n_actual
        DO  j = 1, l_n_atoms(i)
          a_pos(1,j,i) = a_pos(1,j,i) * pi2
          a_pos(2,j,i) = a_pos(2,j,i) * pi2
          a_pos(3,j,i) = a_pos(3,j,i) * pi2
        END DO
      END DO

      RETURN
      END SUBROUTINE nmcoor

! ______________________________________________________________________
!!S Ofiles
! ______________________________________________________________________
! Title: OPTIMZ
! Author: MWD and MMJT
! Date: 8 Apr 1992; 15 Mar 1995; 24 July 1997
! Description:  This routine determines if any shortcuts can be taken
! in the calculations.

!      ARGUMENTS:
!            rootnam  -  The name of the input data file. (input).
!            ok       -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, d0, n_layers, l_alpha, l_actual, l_n_atoms,
!                   a_B, r_B11, r_B22, r_B33, r_B12, r_B23, r_B31,
!                   l_symmetry, l_seq, DoSymDump, SymGrpNo, lambda,
!                   l_cnt, CENTRO, PI, l_r, n_actual, finite_width

!        modifies:  there, one_B, Bs_zero, same_Bs, only_real, l_rz,
!                   same_layer, max_var, no_trials, tolerance, SymGrpNo,
!                   lambda, check_sym, has_l_mirror, theta1, theta2
!                   h_end, h_start, k_end, k_start, pnt_grp, same_rz
!                   xplcit, formfactor, all_Bs_zero,
!                   a_B11,a_B22,a_B33,a_B12,a_B23,a_B31
! ______________________________________________________________________

      SUBROUTINE optimz(rootnam, ok)


!     Utiliza las variables: DIFFaX.par ,DIFFaX.inc',rootnam ,ok ,sym_fnam  ,EQUALB, BINPOW
!                           GET_SYM, i, j, j2, m, n, LENGTH ,HKANGL, h_val, k_val ,x, error,
!                           tmp, incr, z, old_lambda , did_it(MAX_L,MAX_L)
!     Utiliza las funciones: GET_SYM, LENGTH, EQUALB, BINPOW externas    ,HKANGL(k_val,h_val)
!     Utiliza las subrutinas:GETFNM, GET_BDS, CHK_SYM, NMCOOR, OVERLP , thresh(ok)

      CHARACTER (LEN=*), INTENT(IN OUT)        :: rootnam
      LOGICAL, INTENT(IN OUT)                  :: ok



      CHARACTER (LEN=31) :: sym_fnam

      INTEGER*4  i, j, j2, m, n
      REAL*8 hkangl, h_val, k_val
      REAL*8 x, error, tmp, incr, z, old_lambda
      LOGICAL :: did_it(max_l,max_l)

! external functions

! external subroutines (Some compilers need them declared external)
!      external GETFNM, GET_BDS, CHK_SYM, NMCOOR, OVERLP

! statement function
! HKANGL is the angle between the vector (h_val,k_val,0) and (1,0,0)
      hkangl(k_val,h_val) = ATAN2(k_val*SQRT(a0*b0 - d0*d0*quarter),  &
          (h_val*a0 + k_val*d0*half))

! set up logic table for stacking transitions
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          there(j,i) = l_alpha(j,i) >= eps7
        END DO
      END DO

! see if there are any overlapping atoms
    !  CALL overlp()

! multiply all atom coordinates by 2*PI
   !   CALL nmcoor()

! If calculation is to be recursive for a finite number of layers,
! then store the binary form of l_cnt+1 in an array for efficient
! matrix multiplication.
      IF(recrsv .AND. .NOT. inf_thick) THEN
        ok = binpow(l_cnt+1)
        IF(.NOT. ok) THEN
          WRITE(op,202) 'ERROR returned by BINPOW to OPTIMZ'
          WRITE(op,201) 'The argument passed was l_cnt+1 = ', nint(l_cnt)+1
          GO TO 999
        END IF
      END IF

! see if Debye-Waller coefficients are same for all atoms in a layer
      DO  i = 1, n_layers
        x = zero
        j2 = l_actual(i)
        m = l_n_atoms(j2)
        j = 0
        DO  j = 1, m
          x = x + a_b(j,j2)
        END DO
        x = x / DBLE(m)
        error = zero
! find absolute deviation of coefficients
        DO  j = 1, m
          error = error + ABS(a_b(j,j2) - x)
        END DO
! get relative error
        IF(x /= zero) error = error / (x * m)
        one_b(j2) = ABS(error) <= eps3
      END DO

     ! Check that the layer uncertainty factors are physically reasonable.
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          IF(there(j,i)) THEN
           ! check on r_B12
            x = r_b11(j,i)*r_b22(j,i)*a0*b0 - r_b12(j,i)*r_b12(j,i)*ab0*ab0
            IF(x < zero) THEN
              WRITE(op,500) 'C12'
              WRITE(op,501) i, j
              WRITE(op,502) 'C11 and C22'
              ok = .false.
              GO TO 999
            END IF
           ! check on r_B23
            x = r_b22(j,i)*r_b33(j,i)*b0*c0 - r_b23(j,i)*r_b23(j,i)*bc0*bc0
            IF(x < zero) THEN
              WRITE(op,500) 'C23'
              WRITE(op,501) i, j
              WRITE(op,502) 'C22 and C33'
              ok = .false.
              GO TO 999
            END IF
            ! check on r_B31
            x = r_b11(j,i)*r_b33(j,i)*a0*c0 - r_b31(j,i)*r_b31(j,i)*ca0*ca0
            IF(x < zero) THEN
              WRITE(op,500) 'C13'
              WRITE(op,501) i, j
              WRITE(op,502) 'C11 and C33'
              ok = .false.
              GO TO 999
            END IF
          END IF
        END DO
      END DO

      ! see if stacking 'uncertainty' coefficients are same for all layers.
      ! Special flag if they are all zero.
      all_bs_zero = .true.
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          bs_zero(j,i) = r_b11(j,i) == zero .AND. r_b22(j,i) == zero .AND.  &
              r_b33(j,i) == zero .AND. r_b12(j,i) == zero .AND.  &
              r_b23(j,i) == zero .AND. r_b31(j,i) == zero
          all_bs_zero = all_bs_zero .AND. bs_zero(j,i)
        END DO
      END DO

      ! run through all 6 coefficients
      ! until it is clear that some are different
      same_bs = equalb(r_b11, a_b11)
      IF(same_bs) same_bs = equalb(r_b22, a_b22)
      IF(same_bs) same_bs = equalb(r_b33, a_b33)
      IF(same_bs) same_bs = equalb(r_b12, a_b12)
      IF(same_bs) same_bs = equalb(r_b23, a_b23)
      IF(same_bs) same_bs = equalb(r_b31, a_b31)

      ! see if all layers are centrosymmetric
      only_real = .true.
      DO  i = 1, n_actual
        only_real = only_real .AND. (l_symmetry(i) == centro)
      END DO

      ! see if all z-components of the stacking vectors are the same
      l_rz = zero
      n = 0
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          IF(there(j,i)) THEN
            l_rz = l_rz + l_r(3,j,i)
            n = n + 1
          END IF
        END DO
      END DO
      l_rz = l_rz / DBLE(n)
      error = zero
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          IF(there(j,i)) error = error + ABS(l_r(3,j,i) - l_rz)
        END DO
      END DO
      same_rz = ABS(error) <= eps4

      ! If the stacking is explicit, check to see if all the layers are
      ! the same
      same_layer = .false.
      IF(xplcit) THEN
        IF(nint(l_cnt) == 1) GO TO 140
        same_layer = .true.
        j = l_seq(1)
        i = 2
        150   IF(l_seq(i) == j) THEN
          i = i + 1
          IF(i <= nint(l_cnt)) GO TO 150
        ELSE
          same_layer = .false.
        END IF

        ! Check if any of the layer transitions have non-zero probability
        ! initialize flags so that we do not produce duplicate error messages
        DO  i = 1, n_layers
          DO  j = 1, n_layers
            did_it(j,i) = .false.
          END DO
        END DO

        ! now check for legal transitions
        WRITE(op,400)
        DO  n = 1, nint(l_cnt)-1
          i = l_seq(n)
          j = l_seq(n+1)
          IF(.NOT.there(j,i)) THEN
            ok = .false.
            IF(.NOT.did_it(j,i)) THEN
              did_it(j,i) = .true.
              WRITE(op,401) j, i
              STOP
            END IF
          END IF
        END DO
        140   CONTINUE
      END IF
      IF(.NOT.ok) GO TO 999

      ! Pre-compute a pseudo-Lorentzian form factor for lateral
      ! planar widths.
      IF(finite_width) THEN
        ! ideally, FFACT_SIZE is odd
        m = ffact_size/2
        ! incr contains the correct increment size such that the first and
        ! last array elements contain zero.
        incr = (three*n_sigmas*n_sigmas + one) / (two*n_sigmas*m)
        ffact_scale = incr
        z = one + n_sigmas*n_sigmas
        tmp = one / (z*z)
        ! put the peak in the middle (if FFACT_SIZE is odd)
        formfactor(m+1) = one
        DO   n = 1, m
          z = n*incr
          IF(z <= DBLE(n_sigmas)) THEN
            ! Lorentzian part
            x = one / (one + z*z)
          ELSE
            ! Linear part
            x = tmp*(three*n_sigmas*n_sigmas + one - two*n_sigmas*z)
          END IF
          IF(m+n+1 <= ffact_size) formfactor(m+n+1) = x
          IF(m-n+1 > 0)          formfactor(m-n+1) = x
        END DO
        ! compute half width in reciprocal Angstroms for use in computing the
        ! subtle impact of a-b shape broadening on 00l reflections
        tmp = wa*SIN(pi-cell_gamma)
        ffwdth = SQRT(one/(tmp*tmp) + one/(wb*wb))
      END IF
      ! get diffraction symmetry

      ! establish some bounds
      no_trials = 25
      max_var = zero
      ! save lambda, HKL_LIM (called by THRESH), may change it.
      old_lambda = lambda
      CALL thresh(ok)
      IF(.NOT. ok) GO TO 999

      IF(symgrpno == 11)   WRITE(op,202) 'Axial integration only selected.'

      IF(dosymdump) THEN
        CALL getfnm(rootnam, sym_fnam, 'sym', ok)
        IF(.NOT. ok) THEN
          WRITE(op,202) 'OPTIMZ: ERROR in creating symmetry dumpfile.'
          GO TO 999
        END IF
        IF(sy /= op) OPEN(UNIT = sy, FILE = sym_fnam, STATUS = 'new')
        WRITE(op,204) 'Writing symmetry data to file ''',  &
            sym_fnam(1:length(sym_fnam)),'''. . .'
        WRITE(sy,303) '''', rootnam(1:length(rootnam)),''''
        WRITE(sy,203) no_trials
        WRITE(sy,206) tiny_inty
      END IF

      ! restore lambda.
      lambda = old_lambda

      has_l_mirror = symgrpno /= 1 .AND. symgrpno /= 3 .AND.  &
          symgrpno /= 5 .AND. symgrpno /= 6 .AND. symgrpno /= 11

      IF(dosymdump) THEN
        WRITE(sy,202) ' '
        WRITE(sy,204) 'The diffraction data fits the point group symmetry ''',  &
            pnt_grp(1:length(pnt_grp)),''''
        IF(symgrpno /= 1 .AND. symgrpno /= 11) THEN
          IF(max_var > eps6 .AND. max_var <= eps1) THEN
            WRITE(sy,201) '  with a tolerance of one part in ', nint(one / max_var)
          ELSE IF(max_var > eps1) THEN
            WRITE(sy,205) '  with a tolerance of one part in ', one / max_var
          ELSE
            WRITE(sy,202) '  with a tolerance better than one part in a million.'
          END IF
        ELSE
          WRITE(sy,202)  &
              'By definition, all diffraction data has a center of symmetry'
          WRITE(sy,202) 'thus, there is no need to test for inversion.'
        END IF
        ! close file, unless the output was to the default device
        IF(sy /= op) CLOSE (UNIT = sy)
      END IF
      ! establish integration limits and weighting factors
      CALL get_bds()
      ! compute angles of scanning vectors relative to 1 0 0
      theta1 = hkangl(k_start, h_start)
      theta2 = hkangl(k_end, h_end)
      ! resolve ambiguity in the case of -1 point symmetry
      IF(symgrpno == 1 .OR. symgrpno == 11) THEN
        theta1 = -pi
        theta2 =  pi
      END IF
      999 RETURN
      200 FORMAT(1X, 2A)
      201 FORMAT(1X, a, i6)
      202 FORMAT(1X, a)
      203 FORMAT(1X, 'Number of trials per symmetry element = ', i4)
      204 FORMAT(1X, 3A)
      205 FORMAT(1X, a, f4.1)
      206 FORMAT(1X, "Threshold intensity = ", g22.6)
      302 FORMAT(1X, 'SYMMETRY EVALUATIONS FOR DATA IN FILE ''', a, '''')
      303 FORMAT(1X, 'SYMMETRY EVALUATIONS FOR DATA IN FILE ', 3A)
      400 FORMAT(1X,'Checking for conflicts in layer stackings . . .')
      401 FORMAT(1X,'ERROR: Layer ',i2,' cannot stack after layer ',i2)
      500 FORMAT(1X,'ERROR: Non-physical interlayer uncertainty factor ',a)
      501 FORMAT(1X,'       for stacking from layer ', i2, ' to ', i2)
      502 FORMAT(1X,'       It is too large relative to ', a)
      END SUBROUTINE optimz

! ______________________________________________________________________
! Title: OVERLP
! Authors: MMJT
! Date: 24 Feb 1990
! Description: This function compares the coordinates of atoms in
! adjacent layers, and searches for any overlap. If it finds any
! overlap, it checks the summed occupancy to check that it does
! not exceed 1. If OVERLP finds atoms too close it provides a warning
! message only.

!      ARGUMENTS:
!            No arguments are passed.

!      COMMON VARIABLES:
!            uses:  n_layers, there, l_n_atoms, l_r, a_pos

!        modifies:  No COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE overlp()

!     Utiliza las variables: DIFFaX.par ,DIFFaX.inc',  invert    ,txt, i, j, m, n, nn,
!                            j2, err_no, fact, at_num  , prune ,max_err = 100
!                             lay(3,2*max_a),x1, y1, z1, x2, y2, z2, sum_occ,  tmp, bounds , tol = eps1
!     Utiliza las funciones: bounds, prune   externas
!     Utiliza las subrutinas:


      LOGICAL :: invert
      CHARACTER (LEN=33) :: txt
      INTEGER*4 i, j, m, n, nn, j2, err_no, max_err, fact, at_num
      PARAMETER(max_err = 100)
      REAL*8 lay(3,2*max_a)
      REAL*8 x1, y1, z1, x2, y2, z2, sum_occ, tol, tmp
      PARAMETER(tol = eps1)

      ! external functions
     !EXTERNAL bounds

      WRITE(op,fmt="(a)") ' => Checking for conflicts in atomic positions . . .'

      err_no = 0
      DO  i = 1, n_layers
        fact = 1
        invert = l_symmetry(i) == centro
        IF(invert) fact = 2
        at_num = l_n_atoms(i)
        DO  j2 = 1, at_num
          lay(1,j2) = bounds(a_pos(1,j2,i))
          lay(2,j2) = bounds(a_pos(2,j2,i))
         ! Remember, only the a and b directions are truly periodic.
         ! The scaling along c is arbitrary. Furthermore, c is perpendicular
         ! to a and b, and is not necessarily parallel to Rii, which (if it
         ! exists) would define the third cell-repeat direction. In other words,
         ! for the general case, we cannot use the BOUNDS function along c.
          lay(3,j2) = a_pos(3,j2,i)
        END DO
        IF(invert) THEN
          DO  j2 = 1, at_num
            lay(1,at_num+j2) = bounds(-a_pos(1,j2,i))
            lay(2,at_num+j2) = bounds(-a_pos(2,j2,i))
            lay(3,at_num+j2) = -a_pos(3,j2,i)
          END DO
        END IF

        DO  m = 1, at_num
          x1 = lay(1,m)
          y1 = lay(2,m)
          z1 = lay(3,m)
          DO  n = m + 1, fact * at_num
            IF(n > at_num) THEN
              nn = n - at_num
            ELSE
              nn = n
            END IF
            x2 = lay(1,n)
            y2 = lay(2,n)
            z2 = lay(3,n)
            IF(ABS(x1-x2)*cell_a > tol) CYCLE
            IF(ABS(y1-y2)*cell_b > tol) CYCLE
            IF(ABS(z1-z2)*cell_c <= tol) THEN
              sum_occ = a_occup(nn,i) + a_occup(m,i)
              IF((sum_occ - one) > eps4) THEN
                IF(n <= at_num) THEN
                  txt = 'are too close in layer'
                ELSE
                  txt = '(inverted) are too close in layer'
                END IF
                WRITE(op,410) 'Atom ', a_name(nn,i), a_number(nn,i),  &
                    ' and atom ', a_name(m,i), a_number(m,i)
                WRITE(op,412) trim(txt), i
                WRITE(op,420) 'Their combined occupancy is ', sum_occ
                err_no = err_no + 1
                IF(err_no > max_err) GO TO 999
               END IF
            END IF
          END DO
        END DO
      END DO

      ! now let's look at i-j layer transitions and generate a simple warning
      ! message if it seems that the layers are intertwined.
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          IF(there(j,i) .AND. i /= j) THEN
            tmp = l_r(3,j,i) + low_atom(l_actual(j)) - high_atom(l_actual(i))
            IF(tmp*cell_c <= -tol) THEN
              WRITE(op,430) 'Atoms from layer', j, ' extend into layer', i
            END IF
          END IF
        END DO
      END DO

      if(err_no == 0) then
      WRITE(op,fmt="(a)") ' No overlap of atoms has been detected'
      end if

      RETURN
      999 WRITE(op,401) 'WARNING: Number of errors exceeds ', max_err
      RETURN
      400 FORMAT(1X, a)
      401 FORMAT(1X, a, i5)
      410 FORMAT(1X, 'WARNING: ', 2(2A, ' (number ', i3, ')' ) )
      412 FORMAT(10X, a, 1X, i2)
      420 FORMAT(1X, a, g12.5)
      430 FORMAT(1X, 'WARNING: ', 2(a, i3))
      END SUBROUTINE overlp


!!S Pfiles
! ______________________________________________________________________
! ______________________________________________________________________
! Title: PNTINT
! Author: MMJT
! Date: 30 July 1989
! Gets the intensity at the point h, k, l. This differs from the
! subroutine POINT only in that PNTINT does not interact with the user.

!      ARGUMENTS:
!            h   -  reciprocal vector h-component. (input).
!            k   -  reciprocal vector k-component. (input).
!            l   -  reciprocal lattice vector l-component. (input).
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, lambda, recrsv, rad_type, X_RAY

!        modifies:  no COMMON variables are modified

!      PNTINT returns the intensity at h, k, l
! ______________________________________________________________________

      REAL*8 FUNCTION pntint(h, k, l, ok)


!     Utiliza las variables: par.f90' ,'inc.f90' ,h,k,l,ok  ,s, angle, w4, intens, inten2, theta, x
!                            f(max_l)
!     Utiliza las funciones: intens, inten2         externas ,  s(h,k,l) , angle(h,k,l), w4(theta)
!     Utiliza las subrutinas: GET_F, XYPHSE, PRE_MAT

      INTEGER*4                    :: h
      INTEGER*4                   :: k
      REAL*8                   :: l
      LOGICAL, INTENT(IN OUT)                  :: ok




      REAL*8 s, angle, w4,   theta, x
      COMPLEX*16 f(max_l)

      ! external functions

      ! external subroutines (Some compilers need them declared external)
      !      external GET_F, XYPHSE, PRE_MAT

      ! statement functions
      ! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
      ! ANGLE is the Bragg angle (in radians) of the h,k,l plane
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))
      ! W4 is the X-ray polarization factor
      w4(theta) = half * (one + (COS(two*theta))**2)


      CALL xyphse(h, k)

      CALL pre_mat(h, k)
      CALL get_f(f, s(h,k,l), l)
      IF(recrsv) THEN
        x = intens(f, h, k, l, ok)
      ELSE
        x = inten2(f, h, k, l, ok)
      END IF
      IF(.NOT.ok) THEN
        IF(recrsv) THEN
          WRITE(op,200) 'INTENS'
        ELSE
          WRITE(op,200) 'INTEN2'
        END IF
        WRITE(op,201) 'h,k,l = ', h,',', k,',', l
      END IF
      IF(rad_type == x_ray)  x = x * w4(angle(h,k,l))

      pntint = x

      RETURN
      200 FORMAT(1X, 'ERROR returned from ', a, ' in subroutine PNTINT')
      201 FORMAT(1X, 2(a, i3), a, g12.5)
      END FUNCTION pntint

! ______________________________________________________________________
! Title: POINT
! Author: MWD and MMJT
! Date: 18 Aug 1988
! Description:  This subroutine prompts the user for h, k, l values
! and displays the intensity at that point.

!      ARGUMENTS:
!            ok  -  logical flag indicating all went well. (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, lambda, cntrl, CFile, xplcit,
!                   recrsv, rad_type, n_layers, inf_thick, wavefn

!        modifies:  no COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE point(ok)


!     Utiliza las variables: par.f90' ,'inc.f90' ,   get_s, get_s2 ,h, k, i ,chr ,l, intens, inten2,
!                            shkl, x, angle, ss, w4, theta, q2 , f(max_l), s(max_l)
!     Utiliza las funciones: intens, inten2, get_s, get_s2      externas ,  ss(h,k,l) , angle(h,k,l), w4(theta)
!     Utiliza las subrutinas:XYPHSE, PRE_MAT, GET_MAT, GET_F

      LOGICAL, INTENT(IN OUT)                  :: ok



      INTEGER*4 h, k, i
      CHARACTER (LEN=1) :: chr
      REAL*8 l,  shkl, x, angle, ss, w4, theta, q2
      COMPLEX*16 f(max_l), s(max_l)

     ! external functions

     ! external subroutines (Some compilers need them declared external)
     !      external XYPHSE, PRE_MAT, GET_MAT, GET_F

      ! statement functions
      ! SS is the value of 1/d**2 at hkl
      ss(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
      ! ANGLE is the Bragg angle (in radians) of the h,k,l plane
      angle(h,k,l) = ASIN(half * lambda * SQRT(ss(h,k,l)))
      ! W4 is the X-ray polarization factor
      w4(theta) = half * (one + (COS(two*theta))**2)

      1 WRITE(op,400) ' '
      !  WRITE(op,400) 'CALCULATING INTENSITY AT A POINT. . .'
      10 WRITE(op,400) 'Enter h, k, l : '
      READ(cntrl,*,ERR=10)  h, k, l
      IF(cfile) WRITE(op,401) h, k, l
      ! check angles are meaningful
      q2 = four / (lambda**2)
      shkl = ss(h,k,l)
      ! Check that SS(h,k,l) is legal
      IF(shkl < zero) THEN
        WRITE(op,403) 'ERROR: In POINT(), 1/d**2 = ', shkl
        ok = .false.
        GO TO 999
      END IF
      IF(shkl > q2) THEN
        WRITE(op,402) h, k, l,' exceeds 180 degree scattering angle!'
        WRITE(op,"(a)") ' Re-enter. . . '
        GO TO 10
      END IF
      ! make sure we are not going to blow up at the origin
      IF(h == 0 .AND. k == 0 .AND. rad_type == electn) THEN
        IF(shkl <= eps4) THEN
          WRITE(op,"(a)") ' Cannot integrate across the origin for electron radiation'
          WRITE(op,"(a)") ' Re-enter. . .'
          GO TO 10
        END IF
      END IF
      ! get intensity

      CALL xyphse(h, k)

      CALL pre_mat(h, k)
      CALL get_f(f, shkl, l)
      IF(recrsv) THEN
        x = intens(f, h, k, l, ok)
        IF(.NOT.ok) GO TO 999
        ! set up mat again to re-call GET_S (or GET_S2)
        CALL pre_mat(h, k)
        CALL get_mat(h, k, l)
        IF(inf_thick) THEN
          ! initialize s to -f, since mat is -(ident - T)
          DO  i = 1, n_layers
            s(i) = - f(i)
          END DO
          ok = get_s(f, s, h, k, l)
        ELSE
          ! s is initialized within GET_S2, which then calls GET_S
          ok = get_s2(f, s, h, k, l)
        END IF
      ELSE
        x = inten2(f, h, k, l, ok)
      END IF
      IF(.NOT.ok) GO TO 999
      ! for diagnostics write out the s values as well
      IF(rad_type == x_ray)  x = x * w4(angle(h,k,l))
      WRITE(op,404) '2theta = ', rad2deg * two * angle(h,k,l), ' degrees'
      IF(shkl > zero) WRITE(op,403) 'd = ', one / SQRT(shkl)
      IF(shkl >= zero) WRITE(op,403) '1/d = ', SQRT(shkl)
      WRITE(op,"(a)") ' '
      WRITE(op,"(a)") ' Layer scattering factors'
      DO  i = 1, n_layers
        WRITE(op,405) i, f(i)
      END DO
      WRITE(op,"(a)") ' '
      IF(recrsv) THEN
        WRITE(op,"(a)") ' Average scattered wavefunctions'
        DO  i = 1, n_layers
          WRITE(op,407) i, s(i)
        END DO
      ELSE
        WRITE(op,408) 'Crystal wavefunction', wavefn
      END IF
      WRITE(op,"(3a)") ' '
      WRITE(op,406) 'Intensity at ', h, k, l, ' = ', x
      WRITE(op,"(3a)") 'Another point? (y/n)'
      READ(cntrl,'(a)') chr
      IF(chr(1:1) == 'y' .OR. chr(1:1) == 'Y') GO TO 1

      999 RETURN
      400 FORMAT(1X, a)
      401 FORMAT(1X, 2I3, g12.5)
      402 FORMAT(1X, 2I3, g12.5, a)
      403 FORMAT(1X, a, g12.5)
      404 FORMAT(1X, a, g12.5, a)
      405 FORMAT(1X, 'f(', i2, ') = (', g14.7, ',', g14.7, ')')
      406 FORMAT(1X, a, 2I3, g12.5, a, g14.7)
      407 FORMAT(1X, 'psi(', i2,') = (', g14.7, ',', g14.7, ')')
      408 FORMAT(1X, a, ' psi = (', g14.7, ',', g14.7, ')')
      END SUBROUTINE point

! ______________________________________________________________________
! Title: POLINT
! Author: MMJT, adapted from Numerical Recipes Software
! Date: Copyright (C) 1985
! Description: Given arrays xa and ya, each of length n, and given
! a value x, this routine returns an interpolated estimate for y(x),
! and an error estimate dy. If P(x) is the polynomial of degree n-1,
! such that P(xa_i) = ya_i, i=1,....N, then the returned value
! y = P(x). This routine is called by APPR_F.

!      ARGUMENTS:
!            xa  -  Array of length n containing the x values
!                   corresponding to the n ya values. (input).
!            ya  -  Array of length n containing input y values. (input).
!            n   -  Length of the arrays entered. (input).
!            x   -  x position at which an interpolated value is
!                   required. (input).
!            y   -  interpolated y value. (output).
!            dy  -  estimate of the error in y. (output).
!            ok  -  logical flag indicating all went well. (output).
! ______________________________________________________________________

      SUBROUTINE polint(xa,ya,n,x,y,dy,ok)
!     save

!     Utiliza las variables: par.f90' ,xa(n),ya(n) , n,x,y,dy,ok ,i, m, ns ,nmax = 10,dif, dift, ho, hp
!                            c(nmax), d(nmax), w, den
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      REAL*8, INTENT(IN)                       :: xa(n)
      COMPLEX*16, INTENT(IN)                   :: ya(n)
      INTEGER*4, INTENT(IN)                    :: n
      REAL*8, INTENT(IN OUT)                   :: x
      COMPLEX*16, INTENT(IN OUT)               :: y
      COMPLEX*16, INTENT(IN OUT)               :: dy
      LOGICAL, INTENT(IN OUT)                  :: ok





      INTEGER*4 :: i, m, ns
      INTEGER*4, PARAMETER :: nmax = 10
      REAL*8 dif, dift, ho, hp
      COMPLEX*16 c(nmax), d(nmax), w, den

      ns = 1
      dif = ABS(x - xa(1))
      DO  i = 1, n
        dift = ABS(x - xa(i))
        IF(dift < dif) THEN
          ns = i
          dif = dift
        END IF
        c(i) = ya(i)
        d(i) = ya(i)
      END DO
      y = ya(ns)
      ns = ns - 1
      DO  m = 1, n - 1
        DO  i = 1, n - m
          ho = xa(i) - x
          hp = xa(i + m) - x
          w = c(i+1) - d(i)
          den = DCMPLX(ho - hp, zero)
          IF(den == c_zero) GO TO 999
          den = w / den
          d(i) = hp * den
          c(i) = ho * den
        END DO
        IF(2*ns < n-m) THEN
          dy = c(ns + 1)
        ELSE
          dy = d(ns)
          ns = ns - 1
        END IF
        y = y + dy
      END DO

      RETURN
      999 ok = .false.
      WRITE(op,100) 'ERROR: Zero denominator in POLINT.'
      RETURN
      100 FORMAT(1X, a)
      END SUBROUTINE polint

! ______________________________________________________________________
! Title: PRE_MAT
! Author: MMJT
! Date: 21 Mar 1990; 21 July 1997
! Description:  This subroutine premultiplies some of the factors
! needed to calculate the matrix mat ( = Identity - alphaij * Rij)
! at (h,k,0) which remain constant when integrating along a streak. The
! pre-multiplied factors are stored in mat1, and fatsWalla_hk.

!      ARGUMENTS:
!            h   -  reciprocal lattice vector h-component. (input).
!            k   -  reciprocal lattice vector k-component. (input).

!      COMMON VARIABLES:
!            uses:   n_layers, l_r, there, detune, a0, b0, ab0, r_B11,
!                    r_B22, r_B12, same_Bs, Bs_zero, PI2, l_alpha

!        modifies:   l_phi, mat1, fatsWalla_hk
! ______________________________________________________________________

      SUBROUTINE pre_mat(h, k)


!     Utiliza las variables: par.f90' , 'inc.f90', h, k, dot, i, j,
!     Utiliza las funciones:
!     Utiliza las subrutinas:


      INTEGER*4, INTENT(IN)                :: h
      INTEGER*4, INTENT(IN)                    :: k


      REAL*8 dot
      INTEGER*4 i, j

      ! Set up matrix that represents the sequences
      ! For the matrix inversion routines, 'mat' and 'mat1' have to be
      ! in 'i,j' format rather than the quicker 'j,i' format
      DO  i = 1, n_layers
        DO  j = 1, n_layers
          IF(there(j,i)) THEN
            dot = pi2*(h*l_r(1,j,i) + k*l_r(2,j,i))
            l_phi(j,i) = DCMPLX(COS(dot),SIN(dot))
            IF(same_bs.OR.bs_zero(j,i)) THEN
              mat1(i,j) =  detune(j,i) * l_alpha(j,i) * l_phi(j,i)
            ELSE
              ! h-k components only. l-components are handled later by GET_MAT.
              mat1(i,j) = detune(j,i) * l_alpha(j,i) * l_phi(j,i)  &
                  * EXP(-quarter*(r_b11(j,i)*a0*h*h + r_b22(j,i)*b0*k*k)  &
                  + half*r_b12(j,i)*ab0*h*k )
            END IF
          ELSE
            mat1(i,j) = c_zero
          END IF
        END DO
      END DO

      ! Are all the uncertainty factors identical?
      ! Here, we compute only the h-k components.
      ! The l-components are handled later by GET_MAT.
      IF(same_bs) THEN
        IF(all_bs_zero) THEN
          ! This initialization is not actually necessary, since if we are here,
          ! fatsWalla_hk will not be needed by GET_MAT. However, let's be safe.
          fatswalla_hk = one
        ELSE
          fatswalla_hk = EXP(-(quarter*(a_b11*a0*h*h + a_b22*b0*k*k) +  &
              half*a_b12*ab0*h*k) )
        END IF
      END IF

      RETURN
      END SUBROUTINE pre_mat

! ______________________________________________________________________
! ______________________________________________________________________
! Title: PV
! Author: MMJT
! Date: 21st Dec 1990; 7 Mar 1995; 9 June 1998
! Description: This subroutine performs the pseudo-Voigt
! instrumental broadening. Expects the standard u, v and w and gamma
! parameters to be available in 'common'. PV does not conserve
! intensity well when the broadening width is comparable to d_theta.
! Data at the extreme ends of the spectrum are corrupted slightly.

!      ARGUMENTS:
!            th2_low  -  lowest 2theta angle to consider. (input).

!      COMMON VARIABLES:
!            uses:  pv_u, pv_v, pv_w, pv_gamma, d_theta, spec
!                   NONE, PI, RAD2DEG, blurring, th2_max

!        modifies:  brd_spc
! ______________________________________________________________________


      SUBROUTINE pv(th2_low)


!     Utiliza las variables: par.f90' , 'inc.f90' ,th2_low  ,i, j, n_low, n_high, indx
!                            th_rng, tn_th, c00, hk_inv, th0,k1, k2, k3, k4, k5, pvoigt, const, tmp, speci
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      REAL*8, INTENT(IN OUT)    :: th2_low


      INTEGER*4 i, j, n_low,  indx
      REAL      :: th_rng,th_rng2, tn_th, c00, hk_inv, th0 , temp1, temp2, temp3, temp4, delta_lambda
      real(kind = sp ) ::  pv_hg, pv_hl, tmp
      REAL      :: k1, k2, k3, k4, k5, pvoigt, const, speci,pv1,pv2

      ! first check the numbers
      IF(pv_u == zero .AND. pv_v == zero .AND. pv_w == zero) THEN
        WRITE(op,"(a)") '  pseudo-Voigt parameters are zero in PV()'
        WRITE(op,"(a)") '  pseudo-Voigt instrumental broadening not added'
        blurring = NONE
        RETURN
      END IF
      IF(pv_gamma < zero .OR. pv_gamma > one) THEN
        WRITE(op,"(a)") '  Illegal pseudo-Voigt gamma value in PV()'
        WRITE(op,"(a)") '  pseudo-Voigt instrumental broadening not added'
        blurring = NONE
        RETURN
      END IF

      IF(th2_low < zero .OR. th2_low >= th2_max) THEN
        WRITE(op,"(a, g12.5, a)") '  PV: Cut-off angle ', th2_low,  &
            ' is out of bounds. Angle reset to zero.'
        th2_low = zero
      END IF


      ! th2_low is the angle relative to th2_min
      ! 2*d_theta is the angular step size
      n_low  = INT(half*th2_low/d_theta) + 1
      delta_lambda = lambda2 - lambda
      if ( delta_lambda <= 0 ) delta_lambda = delta_lambda * (-1)
      c00 = four * LOG(two)
      const = two * rad2deg * d_theta
      DO  i = 1, n_high
        brd_spc(i) = zero
      END DO
      k3 = -c00
      th0 = half*th2_min
      temp2 = (4.0 * log(2.0) *lambda * lambda ) * rad2deg*rad2deg
      temp1 = pi * pv_dg * pv_dg
      temp4 = (2*lambda/pi) * rad2deg

      DO  i = n_low, n_high
        ! get tan((2theta)/2)
        tn_th = TAN(i*d_theta + th0)
        temp3 = (COS(i*d_theta + th0))
        pv_hg = (pv_u * tn_th * tn_th) + pv_v * tn_th  +  pv_w + (temp2 / (temp1 * temp3* temp3))
        pv_hg = sqrt(pv_hg)
        pv_hl = pv_x * tn_th + temp4/(pv_dl * temp3)
        call tch(pv_hg,pv_hl,tmp,pv_gamma)
        k1 = pv_gamma*two/pi
        k2 = (one - pv_gamma)*SQRT(c00/pi)
        tmp=tmp*tmp
        IF(tmp <= zero) GO TO 995
        hk_inv = one / SQRT(tmp)
        tmp = (const * hk_inv)**2

        k4 = k1 * hk_inv * const
        k5 = k2 * hk_inv * const
        speci = spec(i)

       if (lambda2 /= 0 ) then
        DO  j = n_low - i, n_high - i
          th_rng = tmp * j*j
          th_rng2 =  (((-2.0*rad2deg*delta_lambda * tn_th / lambda)+ j*const) * hk_inv) **2
          pv1=(k4/(one+four*th_rng) + k5*EXP(k3*th_rng))* speci
          pv2=(k4/(one+four*th_rng2) + k5*EXP(k3*th_rng2))*speci
          pvoigt = pv1 + ratio * pv2
          indx = i + j
          brd_spc(indx) = brd_spc(indx) + pvoigt
        END DO

       else

        DO  j = n_low - i, n_high - i
          th_rng = tmp * j*j
          pvoigt = (k4/(one+four*th_rng) + k5*EXP(k3*th_rng)) * speci
          indx = i + j
          brd_spc(indx) = brd_spc(indx) + pvoigt
        END DO

       end if
      END DO


      RETURN

      995 WRITE(op,"(a, g12.5)")  &
          ' ERROR: pseudo-Voigt spread function is complex at theta = ', i*d_theta
      WRITE(op,"(a)") '   u, v, w parameters are illegal.'
      WRITE(op,"(a)") '   pseudo-Voigt instrumental broadening not added'
      blurring = NONE
      RETURN
      END SUBROUTINE pv

!------------------------------------------------------------------------


      SUBROUTINE TCH(hg,hl,fwhm,eta)
!-------------------------------------------------------------------------
! Calculation of eta and FWHM of the pV-function for the
! T-C-H representation.
!-------------------------------------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(kind=sp), INTENT(IN)     :: hg
      real(kind=sp), INTENT(IN)     :: hl
      real(kind=sp), INTENT(OUT)    :: fwhm
      real(kind=sp), INTENT(OUT)    :: eta
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(kind=sp), PARAMETER :: o1=2.69269, o2=2.42843, o3=4.47163, o4= 0.07842
      real(kind=sp), PARAMETER :: e1= 1.36603, e2=0.47719, e3=0.11116
      real(kind=sp) :: ctl, tlr
!-----------------------------------------------
! There is no exception handling because it is supposed to be perfomed before calling TCH
         ctl=hg**5.0+o1*hg**4.0*hl+o2*hg**3.0*hl**2.0+o3*hg**2.0*hl**3.0+  &
             o4*hg*hl**4.0+hl**5.0
         fwhm=ctl**0.2
         tlr = hl/fwhm
         eta = max(1.0e-06, e1*tlr-e2*tlr*tlr+e3*tlr**3.0)  !eta
        RETURN
    END SUBROUTINE TCH
!---------------------------------------------------------------------------
!!S RDfiles1
!!S Rfiles
! ______________________________________________________________________
! Title: RAN3
! Authors: Press, Flannery, Teukolsky and Vetterling
! Date: Copyright (C) 1985
! Returns a uniform random deviate between 0.0 and 1.0. Set 'idum'
! to any negative value to initialize or reinitialize the sequence.
! This version is modified to return real*8 values, and enforces static
! storage of all local variables by use of the 'save' statement
! (In fact 'seed' is the important variable to save, but we save all
! anyway).

!      ARGUMENTS:
!            idum       -  Set -ve to initialize. (input).

!      RAN3 returns a real random number between 0 and 1
! ______________________________________________________________________

      REAL*8 FUNCTION ran3(idum)

      IMPLICIT NONE
      INTEGER*4, intent (in out)             :: idum
      SAVE


!     Utiliza las variables: idum ,big=4000000.0D0 ,seed=1618033.0D0  ,mz=0.0D0 ,fac=2.5D-7,ma(55)
!                            mj, mk ,iff, ii, i, j, inext, inextp
!     Utiliza las funciones:
!     Utiliza las subrutinas:


      REAL*8, PARAMETER :: big=4000000.0D0
      REAL*8, PARAMETER :: seed=1618033.0D0
      REAL*8, PARAMETER :: mz=0.0D0
      REAL*8, PARAMETER :: fac=2.5D-7
      REAL*8 ma(55)
      REAL*8 mj, mk
      INTEGER*4 iff, ii, i, j, inext, inextp

      DATA iff /0/

      IF(idum < 0 .OR. iff == 0) THEN
        iff = 1
        mj = seed - IABS(idum)
        mj = MOD(mj,big)
        ma(55) = mj
        mk = 1.0D0
        DO  i = 1, 54
          ii = MOD(21*i,55)
          ma(ii) = mk
          mk = mj - mk
          IF(mk < mz) mk = mk + big
          mj = ma(ii)
        END DO
        DO  j = 1, 4
          DO  i = 1, 55
            ma(i) = ma(i) - ma(1 + MOD(i + 30,55))
            IF(ma(i) < mz) ma(i) = ma(i) + big
          END DO
        END DO
        inext = 0
        inextp = 31
        idum = 1
      END IF

      inext = inext + 1
      IF(inext == 56) inext = 1
      inextp = inextp + 1
      IF(inextp == 56) inextp = 1
      mj = ma(inext) - ma(inextp)
      IF(mj < mz) mj = mj + big
      ma(inext) = mj
      ran3 = mj * fac

      RETURN
      END FUNCTION ran3

!!S Sfiles
! ______________________________________________________________________
! Title: SALUTE
! Author: MMJT
! Date: 16th Apr, 2002
! Draws the DIFFaX flag. Saluting is optional.

!      ARGUMENTS:
!           No arguments are needed.
! ______________________________________________________________________

      SUBROUTINE salute(i_write)
      integer, optional, intent(in) :: i_write
      integer :: iw
      if(present(i_write)) then
        iw=i_write
      else
        iw=op
      end if
      WRITE(iw,"(a)") ' _________________________________________________'
      WRITE(iw,"(a)") ' _________________________________________________'
      WRITE(iw,"(a)") '            _______ FAULTS 2014 _______           '
      WRITE(iw,"(a)") ' _________________________________________________'
      WRITE(iw,"(a)") ' _________________________________________________'
      WRITE(iw,"(a)") '                                                  '
      WRITE(iw,"(a)") '     A computer program based in DIFFax for       '
      WRITE(iw,"(a)") '     refining faulted layered structures          '
      WRITE(iw,"(a)") '     Authors: M.Casas-Cabanas  (CIC energiGUNE)   '
      WRITE(iw,"(a)") '              J.Rodriguez-Carvajal (ILL)          '
      WRITE(iw,"(a)") '                                                  '
      WRITE(iw,"(a)") '                [version: Nov. 2014]              '
      WRITE(iw,"(a)") ' _________________________________________________'
      WRITE(iw,"(a)")
      WRITE(iw,"(a)")

      RETURN
      END SUBROUTINE salute

! ______________________________________________________________________
! Title: SFC
! Authors: MWD and MMJT
! Date: 11 Feb 90; 7 Mar 1995
! Description: This routine reads the atomic scattering factor constants
! from the file 'sfname' (usually given as "data.sfc"). The scattering
! factors are used to generate the layer form factors. SFC returns
! .false. if there are too many distinct atom types (ie. more than
! MAX_TA). The value of MAX_TA can be reset in 'DIFFaX.par'.

!      ARGUMENTS:
!           No arguments are used. All data is in 'COMMON'.

!      COMMON VARIABLES:
!            uses:  n_actual, l_n_atoms, a_name, n_atoms, rad_type
!                   NEUTRN, X_RAY, rad_type, sfname

!        modifies:  n_atoms, x_sf, n_sf, e_sf, a_type

!      SFC returns logical .true. if it read the scattering factor data
!      safely.
! ______________________________________________________________________

      LOGICAL FUNCTION sfc()

!     Utiliza las variables: par.f90' ,'inc.f90', ok, new_atom, our_atom, done ,list(max_ta+1)
!                            i, j, n, m, length,NAME ,line ,tmp
!     Utiliza las funciones: length externa
!     Utiliza las subrutinas: ATOMS, touppr(NAME)

      LOGICAL :: ok, new_atom, our_atom, done
      LOGICAL :: list(max_ta+1)
      INTEGER*4 i, j, n, m
      CHARACTER (LEN=4)   :: NAMEa
      CHARACTER (LEN=120) :: line
      REAL*4 tmp

      ! external subroutine (Some compilers need them declared external)
      !      external ATOMS

      WRITE(op,"(3a)") ' => Reading scattering factor datafile''',  &
          sfname(1:length(sfname)),'''. . .'

      sfc = .false.
      DO  i = 1, max_ta
        list(i) = .false.
        atom_l(i) = ' '
      END DO

      m = 0

      ! determine all distinct atom types
      DO  i = 1, n_actual
        DO  j = 1, l_n_atoms(i)
          NAMEa = a_name(j,i)
          CALL Ucase(NAMEa)
          ! see if this name has been seen already
          n = 0
          do
            n = n + 1
            new_atom = NAMEa /= atom_l(n)
            IF(new_atom .AND. n < m) cycle
            exit
          end do
          IF(new_atom) THEN
            m = m + 1
            IF(m > max_ta) GO TO 220
            atom_l(m) = NAMEa
            a_type(j,i) = m

          ELSE
            a_type(j,i) = n

          END IF
        END DO
      END DO
      n_atoms = m
      ! now find data for each atom type in file
      ! pass through file only once
   60 READ(sf, 300, END=90, ERR=200) line
      NAMEa = line(1:4)
      CALL ucase(NAMEa)
      ! see if this is one of the distinct atoms
      i = 0
      70   i = i + 1
      our_atom = NAMEa == atom_l(i)
      IF(.NOT. our_atom .AND. i < n_atoms) GO TO 70
      ! have we read all that we need to know?
      done = .true.
      DO  n = 1, n_atoms
        done = done .AND. list(n)
      END DO
      ! If we're done, close file.
      IF(done) GO TO 90
      IF(our_atom .AND. .NOT. list(i)) THEN
        ! mark this atom's data as read in
        list(i) = .true.
        ! and read it in
        IF(rad_type == x_ray) THEN
          READ(line,310,ERR=210) (x_sf(j,i), j=1, 9)
        ELSE IF(rad_type == neutrn) THEN
          READ(line,320,ERR=210) n_sf(i)
        ELSE
          READ(line,310,ERR=210) (x_sf(j,i), j=1, 9), tmp, e_sf(i)
        END IF
      END IF
      GO TO 60

      90 CLOSE(sf,ERR=240)

      ! see if all the data for each atom has been read in
      ok = .true.
      DO  i = 1, n_atoms
        IF(.NOT.list(i)) THEN
          ok = .false.
          WRITE(op,330) 'ERROR: Data for atom ''', atom_l(i),  &
              ''' NOT FOUND IN FILE ''', sfname(1:length(sfname)), ''''
          CALL atoms()
        END IF
      END DO
      IF(ok) THEN
        sfc = .true.
        WRITE(op,"(3a)") ' => Scattering factor data read in.'
      END IF
      RETURN
      200 WRITE(op,"(3a)") ' => Scattering factor file ''',  &
          sfname(1:length(sfname)), ''' DEFECTIVE.'
      CLOSE(sf,ERR=240)
      RETURN
      210 WRITE(op,"(3a)") ' => ERROR reading scattering factor data.'
      CLOSE(sf,ERR=240)
      RETURN
      220 WRITE(op,402) ' => There are too many types of atoms in layer ', i
      WRITE(op,"(3a)") ' => Atoms recorded so far are:'
      DO  j = 1, max_ta
        WRITE(op,403) '       Atom type ', j, '   ', atom_l(j)
      END DO
      WRITE(op,402) ' => Maximum number of types allowed is ', max_ta
      RETURN
      240 WRITE(op,"(3a)") ' => Unable to close scattering factor file ''',  &
          sfname(1:length(sfname)), '''.'
      RETURN
      300 FORMAT(a)
      310 FORMAT(t5, 10F11.6, i3)
      320 FORMAT(t104, f11.6)
      330 FORMAT(1X, 5A)
      400 FORMAT(1X, a)
      401 FORMAT(1X, 2A)
      402 FORMAT(1X, a, i2)
      403 FORMAT(1X, a, i2, 2A)
      END FUNCTION sfc

! ______________________________________________________________________
! Title: SHARP
! Author: MMJT
! Date: 29 Aug 1991
! Description: This subroutine determines whether or not
! spots on a given h, k row are sharp. It does this by examining
! the intensity at l = 0 (or, if absent, at l = d_l) and comparing
! with the intensity at a nearby l-value. If the peak is sharp, there
! will be a large change in intensity.

!      ARGUMENTS:
!            h    -  reciprocal vector h-component. (input).
!            k    -  reciprocal vector k-component. (input).
!            d_l  -  steps in l-component of the reciprocal lattice
!                    vector that sharp spots are likely to be found at.
!                                                            (input).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, lambda, PI

!        modifies:  no COMMON variables are modified

!      SHARP returns logical .true. if it thinks h, k contains sharp
!      spots.
! ______________________________________________________________________

      LOGICAL FUNCTION sharp(h, k, d_l)


!     Utiliza las variables: par.f90' ,'inc.f90', h,k,d_l, ok,s, pntint, angle, ll,
!                            l, i1, i2, x, theta, l_next
!     Utiliza las funciones:  pntint externas   ,  s(h,k,l)  ,angle(h,k,l) ,ll(theta,h,k)
!     Utiliza las subrutinas: GET_F

      Integer,      Intent(In)  :: h
      Integer,      Intent(In)  :: k
      Real(kind=8), Intent(In)  :: d_l



      Logical       :: ok
      Real(kind=8)  :: s,  angle, ll, l, i1, i2, x, theta, l_next

     ! external subroutine (Some compilers need them declared external)
     !      external GET_F
     ! external functions


      ! statement functions
      ! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
      ! ANGLE is the Bragg angle (in radians) of the h,k,l plane
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))
      ! LL is the maximum allowable l value for a given h,k and theta
      ll(theta,h,k) =  SQRT( ( (two*SIN(theta)/lambda)**2  &
          - h*h*a0 - k*k*b0 - h*k*d0) / c0 )

      sharp = .false.

      ! get the intensity at hkl, with l = 0 initially
      l = zero
      10 i1 = pntint(h, k, l, ok)
      IF(.NOT. ok) GO TO 100
      ! If there is an extinction at hkl, try again at l = l + d_l
      IF(i1 < eps4) THEN
        l = l + d_l
        IF(angle(h,k,l) < half*pi) THEN
          GO TO 10
        ELSE
          GO TO 100
        END IF
      END IF

      ! Define a spot to be sharp if intensity is halved at l = l + d_l/100
      theta = angle(h, k, l)
      x = MIN(d_theta, half*th2_max-theta)
      l_next = ll(theta+x, h, k)
      i2 = pntint(h, k, l+l_next*eps2, ok)
      IF(.NOT.ok) GO TO 100

      sharp = i1 > two*i2

      100 RETURN
      END FUNCTION sharp

! ______________________________________________________________________
! Title: SMUDGE
! Author: MMJT
! Date: 30 Oct 1989; 16 April 1999
! Description: This subroutine convolutes 'array', of length
! 'arrsize' with a Gaussian of standard deviation 'sigma'.
! NOTE: This routine does not conserve area under the convoluted curve.
! It is designed so that the peak heights are unchanged.

!      ARGUMENTS:
!            array    -  The name of the input array. (input).
!            arrsize  -  Size of the input array. (input).
!            sigma    -  Standard deviation of the Gaussian as a
!                        multiple of the array sampling step. (input).
!            ok       -  logical flag indicating all went well.
!                                                      (output).
! ______________________________________________________________________

      SUBROUTINE smudge(array, arrsize, sigma, ok)
!     save

!     Utiliza las variables: par.f90' ,array(arrsize) ,arrsize ,sigma ,ok ,m, i, j ,tmparr(sadsize)
!                             k1, k2, tmp, tmp1, tmp2, gss, normalize
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      INTEGER,                         INTENT(IN )     :: arrsize
      REAL(kind=8),dimension(arrsize), INTENT(IN OUT)  :: array
      REAL(kind=8),                    INTENT(IN OUT)  :: sigma
      LOGICAL,                         INTENT(IN OUT)  :: ok




      INTEGER ::  m, i, j
      REAL(kind=8) :: tmparr(sadsize)
      REAL(kind=8) :: k1, k2, tmp, tmp1, tmp2, gss, normalize

      IF(sigma == zero) RETURN

      DO  i = 1, arrsize
        if(array(i) > maxsad) array(i) = maxsad
        tmparr(i) = array(i)
        array(i) = zero
      END DO

! go out to 5 standard deviations, or to the end of the spectrum
      m = nint(five * sigma)
      k1 = half / ( sigma * sigma )
      IF(m > arrsize) m = arrsize
! get normalization constant. We wish peak heights to be unchanged.
      normalize = one
      DO  i = 1, m
        normalize = normalize + two * EXP( -k1 * DBLE(i*i) )
      END DO

      IF(normalize == zero) THEN
        WRITE(op,100) 'ERROR in SMUDGE: Zero normalization constant.'
        ok = .false.
        GO TO 999
      END IF
      normalize = one / normalize

      DO  i = 0, m
        k2 = k1*DBLE(i*i)
        gss = EXP(-k2)
        DO  j = 1, arrsize
          tmp1 = zero
          tmp2 = zero
          IF((j-i) > 0) tmp1 = tmparr(j-i)
          IF((j+i) <= arrsize) tmp2 = tmparr(j+i)
          tmp = tmp1 + tmp2
          IF(i == 0) tmp = half * tmp
          array(j) = array(j) + gss * tmp * normalize
        END DO
      END DO

      999 RETURN
      100 FORMAT(1X, a)
      END SUBROUTINE smudge

! ______________________________________________________________________
! Title: SPHCST
! Author: MWD
! Date: 18 Aug 1988
! Description: This subroutine determines the constants used in
! determining the magnitude of reciprocal lattice vectors.

!      ARGUMENTS:
!            No input arguments.

!      COMMON VARIABLES:
!            uses:  cell_a, cell_b, cell_c, cell_gamma

!        modifies:  a0, b0, c0, d0, ab0, bc0, ca0
! ______________________________________________________________________

      SUBROUTINE sphcst()


!     Utiliza las variables: par.f90', inc.f90'
!     Utiliza las funciones:
!     Utiliza las subrutinas:


      a0 = one / (cell_a * SIN(cell_gamma))**2
      b0 = one / (cell_b * SIN(cell_gamma))**2
      c0 = one / (cell_c)**2
      d0 = -two * COS(cell_gamma) / (cell_a * cell_b * SIN(cell_gamma)**2)

      ab0 = SQRT(a0 * b0)
      bc0 = SQRT(b0 * c0)
      ca0 = SQRT(c0 * a0)



      RETURN
      END SUBROUTINE sphcst

! ______________________________________________________________________
! Title: STREAK
! Author: MWD & MMJT
! Date: 18 Aug 1988, revised 22 June 1989.
! Description:  This routine outputs the integrated (ie. averaged)
! intensitites from h,k,l0 to h,k,l1, in steps of dl, to the file
! strkfile.

!      ARGUMENTS:
!            FN        -  Function name passed by reference. The
!                         choice is between GLQ16 (non-adaptive
!                         Gauss-Legendre integration), and AGLQ16
!                         (adaptive Gauss-Legendre integration).
!                                                            (input).
!            strkfile  -  The name of the file that streak data is
!                         to be output to. (input).
!            ok        -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, lambda, cntrl, CFile, xplcit,
!                   rad_type, X_RAY

!        modifies:  no COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE streak(fn, strkfile, ok)

!     Utiliza las variables: par.f90', inc.f90'  ,strkfile ,ok ,its_hot, h, k, length, i, i_step
!                            l, theta, x, angle, s, q2, l0, l1 ,dl, w4 ,intervals = twenty ,fn,
!     Utiliza las funciones:  fn, length externas    s(h,k,l) ,angle(h,k,l), w4(theta)
!     Utiliza las subrutinas: XYPHSE, PRE_MAT, GET_F

      CHARACTER (LEN=*), INTENT(IN OUT)        :: strkfile
      LOGICAL, INTENT(IN OUT)                  :: ok

      LOGICAL :: its_hot
      INTEGER*4 h, k, i, i_step, lz, lzf
      REAL*8 l, theta, x, angle, s, q2, l0, l1
      REAL*8 dl, w4
      REAL*8, PARAMETER :: intervals = twenty

! external functions
      REAL*8 fn
      EXTERNAL fn
! external subroutines (Some compilers need them declared external)
!      external XYPHSE, PRE_MAT, GET_F

! statement functions
! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
! ANGLE is the Bragg angle (in radians) of the h,k,l plane
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))
! W4 is the X-ray polarization factor
      w4(theta) = half * (one + (COS(two*theta))**2)

      q2 = four / (lambda**2)
      10 WRITE(op,"(a)") ' Enter h, k, l0, l1, delta l : '
      READ(cntrl,*,ERR=10) h, k, l0, l1, dl
      IF(cfile) WRITE(op,401) h, k, l0, l1, dl
! check input
      IF(l1 == l0) THEN
        WRITE(op,400) 'Illegal input: l0 equals l1'
        GO TO 999
      ELSE IF(dl == zero) THEN
        WRITE(op,400) 'Illegal zero value of dl entered'
        WRITE(op,402)'A value of ',(l1-l0)/(five*hundred),' is assumed'
      ELSE IF(l1 > l0 .AND. dl < zero) THEN
        WRITE(op,400) 'l1 is greater than l0. +ve dl assumed'
        dl = -dl
      ELSE IF(l1 < l0 .AND. dl > zero) THEN
        WRITE(op,400) 'l0 is greater than l1. -ve dl assumed'
        dl = -dl
      END IF
! The origin may be hotter than hell! Let's check first.
      its_hot = (h == 0 .AND. k == 0) .AND. l0*l1 <= zero .AND. rad_type == electn
      IF(its_hot) THEN
        WRITE(op,400) 'Cannot scan the origin with electron radiation'
        WRITE(op,400) 'Origin will be skipped.'
      END IF
! check angles are meaningful
      IF(s(h,k,l0) > q2 .OR. s(h,k,l1) > q2) THEN
        IF(s(h,k,l0) > q2) WRITE(op,403) h, k, l0,  &
            ' exceeds 180 degree scattering angle!'
        IF(s(h,k,l1) > q2) WRITE(op,403) h, k, l1,  &
            ' exceeds 180 degree scattering angle!'
        GO TO 10
      END IF

      WRITE(op,404) 'Writing streak data to file ''',  &
          strkfile(1:length(strkfile)),'''. . .'

      CALL xyphse(h, k)

      CALL pre_mat(h, k)

      i_step = nint( (l1 - l0) / (dl * intervals) )
! Immune system to the rescue!
      IF(i_step <= 0) i_step = 10

      i = 0
!      DO  l = l0, l1, dl
      lzf=nint((l1-l0)/dl+1.0D0)

      DO  lz = 1, lzf
        l=l0+(lz-1)*dl
! If we are dealing with electrons, make sure we avoid the origin
        IF(its_hot .AND. l*(l+dl) <= zero) THEN
          x = zero
          GO TO 30
        END IF
        i = i + 1
        IF(MOD(i,i_step) == 0) WRITE(op,405) 'Reached l = ',l
        x = fn(h,k,l,l+dl,ok)
        IF(.NOT.ok) GO TO 999
! note: since this is streak data, only the X_RAY input needs
! correcting for polarization.
        IF(rad_type == x_ray)  x = x * w4(angle(h,k,l+half*dl))
        30   WRITE(sk,406,ERR=100) l, CHAR(9), x
      END DO
      IF(sk /= op) CLOSE(sk,ERR=110)
      WRITE(op,404) 'Streak data file, ''',  &
          strkfile(1:length(strkfile)),''' WRITTEN TO DISK.'
      RETURN
      100 WRITE(op,404) 'ERROR writing to streak data file ''',  &
          strkfile(1:length(strkfile)),''''
      IF(sk /= op) CLOSE(sk,ERR=110)
      RETURN
      110 WRITE(op,404) 'Unable to close streak data file ''',  &
          strkfile(1:length(strkfile)),''''
      RETURN
      999 WRITE(op,405) 'ERROR encountered in streak integration at l = ',l
      RETURN
      400 FORMAT(1X, a)
      401 FORMAT(1X, 2I3, 3F10.5)
      402 FORMAT(1X, a, f10.5, a)
      403 FORMAT(1X, 2I3, f10.5, a)
      404 FORMAT(1X, 3A)
      405 FORMAT(1X, a, f10.5)
      406 FORMAT(1X, e12.5, a, e14.6)
      END SUBROUTINE streak

!!S Tfiles
! ______________________________________________________________________
! Title: THRESH
! Author: MMJT
! Date: 21 Jan 1995
! Samples intensities in reciprocal space to get a "feeling" for what
! kind of values are out there prior to testing diffraction symmetry.
! CHK_SYM and GET_SYM measure the fractional deviations of intensity
! from (potentially) symmetry-related points in reciprocal space.
! This method runs into problems when the intensity being measured
! is close to zero. Miniscule intensity variations can appear to be
! huge relative to zero!
! This function is needed in order to obtain a (crude) estimate of
! which intensity values are too small to worry about, even if the
! relative intensity variations seem to be large.
! This function will be of no use if there are no streaks, that is,
! if the crystal is perfect.

!      ARGUMENTS:
!            ok       -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, lambda, no_trials,
!                   max_angle

!        modifies:  ok, max_angle, h_bnd, k_bnd, l_bnd, tiny_inty

! ______________________________________________________________________

      SUBROUTINE thresh(ok)


!     Utiliza las variables: par.f90', inc.f90'  ok, i, h, k, idum, ran3, pntint, s, angle , l, tot_int
!     Utiliza las funciones: ran3, pntint  externas    s(h,k,l) ,angle(h,k,l)
!     Utiliza las subrutinas:  hkl_lim()

      LOGICAL, INTENT(IN OUT)                  :: ok


      INTEGER*4 i, h, k, idum
      REAL*8   s, angle
      REAL*8 l, tot_int

! external functions


! statement functions
! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))
! initialize random numbers in RAN3
      idum = -1

! First define angular range to sample. h_bnd, k_bnd, l_bnd
! (defined in HKL_LIM) and max_angle, are used later on
! in GET_SYM and CHK_SYM.
      max_angle = quarter * pi
      CALL hkl_lim()
! Sample the typical intensities that are out there
      tot_int = zero
      DO  i = 1, no_trials
        20   h = INT( DBLE(h_bnd + 1) * ran3(idum) )
        k = INT( DBLE(k_bnd + 1) * ran3(idum) )
        l = l_bnd * ran3(idum)
! make sure we are not sampling at too high an angle
        IF(two * angle(h, k, l) > max_angle) GO TO 20
! get I(h,k,l)
        tot_int = tot_int + pntint(h, k, l, ok)
        IF(.NOT.ok) GO TO 999
      END DO

! Estimate some suitable fraction of the average intensity. This
! fraction defines a baseline intensity equivalent to zero for
! the purposes of the symmetry testing later on.

      tiny_inty = tot_int * eps5 / no_trials

      900 RETURN
      999 WRITE(op,100) 'ERROR in intensity calculation in THRESH'
      WRITE(op,200) '   at h,k,l = ', h,',', k,',', l
      RETURN
      100 FORMAT(1X, a)
      200 FORMAT(1X, a, i3, a, i3, a, f7.2)
      END SUBROUTINE thresh


! ______________________________________________________________________
! Title: TRMSPC
! Author: MMJT
! Date: 20 April 1989; 7 Mar 1995
! Description:  This function locates a suitable cut-off angle,
! below which we can ignore the huge intensities which may occur
! close to the origin when the full adaptive integration option is used.
! TRMSPC is called only when theta = 0 is in the integrated range.

!      ARGUMENTS:
!            th2_low  -  a 2theta cut-off angle to be found (output)

!      COMMON VARIABLES:
!            uses:  th2_min, th2_max, d_theta, spec, RAD2DEG

!        modifies:  spec

!            TRMSPC returns logical .true. if all went well.
! ______________________________________________________________________

      LOGICAL FUNCTION trmspc(th2_low)



!     Utiliza las variables: par.f90, inc.f90,  th2_low, i, i_min, i_max
!     Utiliza las funciones:
!     Utiliza las subrutinas:


      REAL*8, INTENT(IN OUT)                   :: th2_low


      INTEGER*4 i, i_min, i_max

      i_max = INT(half*(th2_max - th2_min) / d_theta) + 1
! spec(1) corresponds to the intensity at the origin and is always zero.
      i = 2
      20 i = i + 1
      IF(i >= i_max+1) THEN
        WRITE(op,100) 'No peaks were found in spectrum.'
        write(*,*) th2_max, th2_min, d_theta
        GO TO 900
      END IF
! locate the first minimum after the huge peak at the origin
      IF(spec(i) <= spec(i-1)) GO TO 20
      i_min = i - 1

! NOTE: The absolute angle is th2_low + th2_min
      th2_low = i_min * d_theta

      900 trmspc = .true.

      RETURN
      100 FORMAT(1X, a)
      END FUNCTION trmspc

! ______________________________________________________________________
! Title: TST_MIR
! Author: MMJT
! Date: 30 July 1989; 22 Feb 1995
! Identifies the presence of mirror planes which contain the streaky
! axis. The integer argument 'mir_sym' can have one of two values,
! 1 or 2.
! mir_sym = 1 is the plane containing the cell sides h and l.
! mir_sym = 2 is the plane containing the cell sides k and l.
! For rotational symmetry greater than 2, one mirror implies the
! presence of the other. For diffraction symmetry group No 3 (2/M),
! only one or the other can occur, but we must test for both.
! TST_MIR returns '.true.' if a mirror was found, '.false.' otherwise.

!      ARGUMENTS:
!            mir_sym  -  Plane across which we wish to test for mirror
!                        symmetry. Takes the values 1 or 2. (input).
!            idum     -  parameter used by RAN3. Is -ve if RAN3
!                        is to be reset. (input).
!            ok       -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, lambda, DoSymDump, no_trials,
!                   max_angle, PI, PI2, RAD2DEG, cell_gamma, check_sym
!                   h_bnd, k_bnd, l_bnd, tolerance, tiny_inty

!        modifies:  max_var

!      TST_MIR returns logical .true. if the diffraction intensities
!      have the mirror symmetry about the plane requested.
! ______________________________________________________________________

      LOGICAL FUNCTION tst_mir(mir_sym, idum, ok)


!     Utiliza las variables: par.f90, inc.f90,mir_sym , idum , ok  ,is_good, match, eq_sides
!                             i, h, k, h_tmp, k_tmp  ,cell90, cell120 ,ran3, pntint, s, angle,l
!                             i_avg, tol, i1, i2, variance, rel_var ,tiny = five * eps4
!     Utiliza las funciones:   ran3, pntint externas    s(h,k,l) ,angle(h,k,l)
!     Utiliza las subrutinas:

      INTEGER*4, INTENT(IN)                :: mir_sym
      INTEGER*4, INTENT(IN out)                :: idum
      LOGICAL,   INTENT(OUT)               :: ok



      LOGICAL :: is_good, match, eq_sides
      INTEGER*4 i, h, k, h_tmp, k_tmp
      LOGICAL :: cell90, cell120
      REAL*8   s, angle
      REAL*8   l
      REAL*8 i_avg, tol, i1, i2, variance, rel_var
      REAL*8, PARAMETER :: tiny = five * eps4

! external functions


! statement functions
! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
! ANGLE is the Bragg angle (in radians) of the h,k,l plane
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))

      cell90  = ABS(cell_gamma - half*pi) < half*pi*tiny
      cell120 = ABS(cell_gamma - pi2/three) < pi2*tiny/three
      eq_sides = ABS(cell_a - cell_b) <= half*eps6*(cell_a + cell_b)
      tst_mir = .false.
      is_good = .false.
      IF(mir_sym < 1 .OR. mir_sym > 3) GO TO 900
      IF(dosymdump) THEN
        WRITE(sy,200)
        IF(.NOT.cell90 .AND. .NOT.cell120) THEN
          WRITE(sy,230) 'cell angle = ', cell_gamma * rad2deg,  &
              ' degrees. NO HORIZONTAL MIRRORS ARE LIKELY.'
          GO TO 900
        END IF
      END IF

      IF(dosymdump) THEN
        IF(mir_sym == 1) THEN
          WRITE(sy,199) 'Testing for mirror about the h-l plane'
        ELSE IF(mir_sym == 2) THEN
          WRITE(sy,199) 'Testing for mirror about the k-l plane'
        ELSE IF(mir_sym == 3) THEN
          WRITE(sy,199) 'Testing for mirror about the h=k,l plane'
        END IF
        WRITE(sy,200)
        WRITE(sy,210)
      END IF
      DO  i = 1, no_trials
! get usable h,k >= 0
        20   IF(mir_sym == 1) THEN
          h_tmp = INT( DBLE(h_bnd + 1) * ran3(idum) )
          30     k_tmp = INT( DBLE(k_bnd + 1) * ran3(idum) )
          IF(k_tmp == 0) GO TO 30
        ELSE IF(mir_sym == 2) THEN
          k_tmp = INT( DBLE(k_bnd + 1) * ran3(idum) )
          40     h_tmp = INT( DBLE(h_bnd + 1) * ran3(idum) )
          IF(h_tmp == 0) GO TO 40
        ELSE IF(mir_sym == 3) THEN
          k_tmp = INT( DBLE(k_bnd + 1) * ran3(idum) )
          45     h_tmp = INT( DBLE(h_bnd + 1) * ran3(idum) )
          IF(h_tmp == k_tmp) GO TO 45
        END IF
! get usable l > 0
        50   l = l_bnd * ran3(idum)
        IF(ABS(l) <= eps2) GO TO 50
! make sure we are not sampling at too high an angle
        IF(two * angle(h_tmp, k_tmp, l) > max_angle) GO TO 20
! I(h,k,l)
        h = h_tmp
        k = k_tmp
        i1 = pntint(h, k, l, ok)
        IF(.NOT.ok) GO TO 999
        IF(dosymdump) WRITE(sy,220) h, k, l, i1
! mirror on h-l plane
        IF(mir_sym == 1) THEN
! is the cell angle equal to 90 degrees
          IF(cell90) THEN
! I(h,-k,l), rectangular cell
            h =  h_tmp
            k = -k_tmp
            i2 = pntint(h, k, l, ok)
            IF(.NOT.ok) GO TO 999
            IF(dosymdump) WRITE(sy,220) h, k, l, i2
          ELSE IF(cell120) THEN
! I(h+k,-k,l), hexagonal cell
            h =  h_tmp + k_tmp
            k = -k_tmp
            i2 = pntint(h, k, l, ok)
            IF(.NOT.ok) GO TO 999
            IF(dosymdump) WRITE(sy,220) h, k, l, i2
          END IF
! else mirror on k-l plane, mir = 2
        ELSE IF(mir_sym == 2) THEN
! is the cell angle equal to 90 degrees
          IF(cell90) THEN
! I(-h,k,l), rectangular cell
            h = -h_tmp
            k =  k_tmp
            i2 = pntint(h, k, l, ok)
            IF(.NOT.ok) GO TO 999
            IF(dosymdump) WRITE(sy,220) h, k, l, i2
          ELSE IF(cell120) THEN
! I(-h,h+k,l), hexagonal cell
            h = -h_tmp
            k =  h_tmp + k_tmp
            i2 = pntint(h, k, l, ok)
            IF(.NOT.ok) GO TO 999
            IF(dosymdump) WRITE(sy,220) h, k, l, i2
          END IF
! else mirror on hk-l plane, mir = 3
        ELSE IF(mir_sym == 3) THEN
! is the cell square
          IF(cell90 .AND. eq_sides) THEN
! I(-h,k,l), square cell
            h = k_tmp
            k = h_tmp
            i2 = pntint(h, k, l, ok)
            IF(.NOT.ok) GO TO 999
            IF(dosymdump) WRITE(sy,220) h, k, l, i2
          ELSE IF(cell120) THEN
! I(-h,h+k,l), hexagonal cell
            h = k_tmp
            k = h_tmp
            i2 = pntint(h, k, l, ok)
            IF(.NOT.ok) GO TO 999
            IF(dosymdump) WRITE(sy,220) h, k, l, i2
          END IF
        END IF
! compare mirrored intensities
        i_avg = half * (i1 + i2)
        variance = half * (ABS(i_avg-i1) + ABS(i_avg-i2))
! Be careful intensities are not actually zero
        IF(i_avg < tiny_inty) THEN
          tol = tiny_inty
        ELSE
          tol = i_avg * tolerance
          rel_var = variance / i_avg
          IF(rel_var > max_var) max_var = rel_var
        END IF
        match = ABS(i_avg-i1) < tol .AND. ABS(i_avg-i2) < tol
        is_good = (i == 1 .OR. is_good) .AND. match
        IF(dosymdump) THEN
          WRITE(sy,270)
          WRITE(sy,240) i_avg
          WRITE(sy,360) variance, hundred * variance / i_avg
!          if(.not.check_sym) then
          WRITE(sy,260) tol
          IF(match) THEN
            IF(mir_sym == 1) THEN
              WRITE(sy,250)
            ELSE IF(mir_sym == 2) THEN
              WRITE(sy,280)
            ELSE IF(mir_sym == 3) THEN
              WRITE(sy,285)
            END IF
          ELSE
            IF(mir_sym == 1) THEN
              WRITE(sy,290)
            ELSE IF(mir_sym == 2) THEN
              WRITE(sy,300)
            ELSE IF(mir_sym == 3) THEN
              WRITE(sy,305)
            END IF
          END IF
!          endif
          WRITE(sy,200)
        END IF
      END DO
      tst_mir = is_good

      IF(dosymdump) THEN
!        if(.not.check_sym) then
        IF(mir_sym == 1) THEN
          IF(is_good) THEN
            WRITE(sy,310)
          ELSE
            WRITE(sy,320)
          END IF
        ELSE IF(mir_sym == 2) THEN
          IF(is_good) THEN
            WRITE(sy,330)
          ELSE
            WRITE(sy,340)
          END IF
        ELSE IF(mir_sym == 3) THEN
          IF(is_good) THEN
            WRITE(sy,345)
          ELSE
            WRITE(sy,346)
          END IF
        END IF
!        endif
        WRITE(sy,200)
      END IF

      900 RETURN
      999 WRITE(op,199) 'ERROR in intensity calculation in TST_MIR'
      WRITE(op,350) '   at h,k,l = ', h,',', k,',', l
      RETURN
      199 FORMAT(1X, a)
      200 FORMAT(' ')
      210 FORMAT(1X, '  h', 5X, 'k', 7X, 'l', 20X, 'Intensity')
      220 FORMAT(1X, i3, 3X, i3, 2X, f9.4, 5X, f22.6)
      230 FORMAT(1X, a, f6.2, a)
      240 FORMAT(6X, 'Average Intensity = ', f22.6)
      250 FORMAT(1X, 'Intensities are consistent with an h-l mirror plane')
      260 FORMAT(1X, 'Intensity tolerance = +/-', f22.6)
      270 FORMAT(26X, '----------------------')
      280 FORMAT(1X, 'Intensities are consistent with a k-l mirror plane')
      285 FORMAT(1X, 'Intensities are consistent with a h=k,l mirror plane')
      290 FORMAT(1X, 'Intensities not consistent with an h-l mirror plane')
      300 FORMAT(1X, 'Intensities not consistent with a k-l mirror plane')
      305 FORMAT(1X, 'Intensities not consistent with a h=k,l mirror plane')
      310 FORMAT(1X, 'THERE IS A MIRROR ABOUT THE H-L PLANE')
      320 FORMAT(1X, 'THERE IS NO MIRROR ABOUT THE H-L PLANE')
      330 FORMAT(1X, 'THERE IS A MIRROR ABOUT THE K-L PLANE')
      340 FORMAT(1X, 'THERE IS NO MIRROR ABOUT THE K-L PLANE')
      345 FORMAT(1X, 'THERE IS A MIRROR ABOUT THE H=K,L PLANE')
      346 FORMAT(1X, 'THERE IS NO MIRROR ABOUT THE H=K,L PLANE')
      350 FORMAT(1X, a, i3, a, i3, a, f7.2)
      360 FORMAT(1X,'  Average variation = +/-', f22.6,'  (+/-',g9.2,'%)')
      END FUNCTION tst_mir

! ______________________________________________________________________
! Title: TST_ROT
! Author: MMJT
! Date: 30 July 1989; 22 Feb 1995
! Identifies the rotational symmetry of the diffraction pattern.
! The rotational symmetry to test for is passed as 'rot_sym'.
! The values accepted for 'rot_sym' are 2, 3 and 4. If the diffraction
! intensity has the symmetry requested, TST_ROT returns '.true.'. If
! the pattern does not have the requested symmetry, or if an illegal
! value was passed (i.e. rot_sym = 5), TST_ROT returns '.false.'.
! NOTE. A 6-fold axis is found by calling TST_ROT twice: once for
! rot_sym = 2 and again for rot_sym = 3.

!      ARGUMENTS:
!            rot_sym  -  Rotational symmetry to test for.
!                        Accepts the values 2, 3 or 4. (input).
!            idum     -  parameter used by RAN3. Is -ve if RAN3
!                        is to be reset. (input).
!            ok       -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, lambda, DoSymDump, no_trials,
!                   max_angle, check_sym, h_bnd, k_bnd, l_bnd
!                   tolerance

!        modifies:  max_var

!      TST_ROT returns logical .true. if the diffraction intensities
!      have the requested rotational symmetry.
! ______________________________________________________________________

      LOGICAL FUNCTION tst_rot(rot_sym, idum, ok)


!     Utiliza las variables: par.f90, inc.f90, rot_sym ,idum , ok ,is_good, match ,i, h, k, h_tmp, k_tmp
!                            ran3, pntint, s, angle  ,l, i_avg, tol  ,i1, i2, i3, i4, variance, rel_var
!     Utiliza las funciones:   ran3, pntint externas    s(h,k,l) ,angle(h,k,l)
!     Utiliza las subrutinas:


      INTEGER*4, INTENT(IN)                    :: rot_sym
      INTEGER*4, INTENT(IN OUT)                :: idum
      LOGICAL, INTENT(IN OUT)                  :: ok



      LOGICAL :: is_good, match
      INTEGER*4 i, h, k, h_tmp, k_tmp
      REAL*8   s, angle
      REAL*8 l, i_avg, tol
      REAL*8 i1, i2, i3, i4, variance, rel_var

! external functions


! statement functions
! S is the value of 1/d**2 at hkl
      s(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
      angle(h,k,l) = ASIN(half * lambda * SQRT(s(h,k,l)))

      tst_rot = .false.
      is_good = .false.
! Is rot valid?
      IF(rot_sym < 2 .OR. rot_sym > 4) GO TO 900
! Now test for rotational symmetry.
! 2-fold and 4-fold
! Avoid both h, k = 0. Also, avoid l = 0, since Friedel's law will
! create a pseudo 2-fold.
      IF(rot_sym == 2 .OR. rot_sym == 4) THEN
        IF(dosymdump) THEN
          WRITE(sy,210)
          WRITE(sy,330) 'Testing for ', rot_sym, '-fold axis'
          WRITE(sy,210)
        END IF
        DO  i = 1, no_trials
          IF(dosymdump) WRITE(sy,220)
! get usable h,k >= 0
          20     h_tmp = INT( DBLE(h_bnd + 1) * ran3(idum) )
          k_tmp = INT( DBLE(k_bnd + 1) * ran3(idum) )
          IF(h_tmp == 0 .AND. k_tmp == 0) GO TO 20
! get usable l > 0
          30     l = l_bnd * ran3(idum)
! keep l off the l = 0 plane, else we might confuse the inversion
! with a 2-fold
          IF(ABS(l) <= eps2) GO TO 30
! make sure we are not sampling at too high an angle
          IF(two * angle(h_tmp, k_tmp, l) > max_angle) GO TO 20
! I(h,k,l)
          h = h_tmp
          k = k_tmp
          i1 = pntint(h, k, l, ok)
          IF(.NOT.ok) GO TO 999
          IF(dosymdump) WRITE(sy,230) h, k, l, i1
! I(-h,-k,l)
          h = -h_tmp
          k = -k_tmp
          i2 = pntint(h, k, l, ok)
          IF(.NOT.ok) GO TO 999
          IF(dosymdump) WRITE(sy,230) h, k, l, i2
! compare 2-fold intensities
          IF(rot_sym == 2) THEN
            i_avg = half * (i1 + i2)
            variance = half * (ABS(i_avg-i1) + ABS(i_avg-i2))
! Be careful intensities are not actually zero
            IF(i_avg < tiny_inty) THEN
              tol = tiny_inty
            ELSE
              tol = i_avg * tolerance
              rel_var = variance / i_avg
              IF(rel_var > max_var) max_var = rel_var
            END IF
            match = ABS(i_avg-i1) < tol .AND. ABS(i_avg-i2) < tol
            is_good = (i == 1 .OR. is_good) .AND. match
          ELSE
! I(-k,h,l)
            h = -k_tmp
            k =  h_tmp
            i3 = pntint(h, k, l, ok)
            IF(.NOT.ok) GO TO 999
            IF(dosymdump) WRITE(sy,230) h, k, l, i3
! I(k,-h,l)
            h =  k_tmp
            k = -h_tmp
            i4 = pntint(h, k, l, ok)
            IF(.NOT.ok) GO TO 999
            IF(dosymdump) WRITE(sy,230) h, k, l, i4
! compare 4-fold intensities
            i_avg = quarter * (i1 + i2 + i3 + i4)
            variance = quarter * (ABS(i_avg-i1) + ABS(i_avg-i2) +  &
                ABS(i_avg-i3) + ABS(i_avg-i4))
! Be careful intensities are not actually zero
            IF(i_avg < tiny_inty) THEN
              tol = tiny_inty
            ELSE
              tol = i_avg * tolerance
              rel_var = variance / i_avg
              IF(rel_var > max_var) max_var = rel_var
            END IF
            match = ABS(i_avg-i1) < tol .AND. ABS(i_avg-i2) < tol  &
                .AND. ABS(i_avg-i3) < tol .AND. ABS(i_avg-i4) < tol
            is_good = (i == 1.OR.is_good) .AND. match
          END IF

          IF(dosymdump) THEN
            WRITE(sy,240)
            WRITE(sy,250) i_avg
            WRITE(sy,260) variance, hundred * variance / i_avg
!            if(.not.check_sym) then
            WRITE(sy,270) tol
            IF(match) THEN
              WRITE(sy,280) rot_sym
            ELSE
              WRITE(sy,290) rot_sym
            END IF
!            endif
            WRITE(sy,210)
          END IF
        END DO
        tst_rot = is_good
        GO TO 900
      END IF
! 3-fold
! Avoid both h, k = 0.
      IF(rot_sym == 3) THEN
        IF(dosymdump) THEN
          WRITE(sy,200) rot_sym
          WRITE(sy,210)
          WRITE(sy,220)
        END IF
        DO  i = 1, no_trials
! get usable h,k >= 0
          50     h_tmp = INT( DBLE(h_bnd + 1) * ran3(idum) )
          k_tmp = INT( DBLE(k_bnd + 1) * ran3(idum) )
          IF(h_tmp == 0 .AND. k_tmp == 0) GO TO 50
! get l (l=0 is allowed)
          l = l_bnd * ran3(idum)
! make sure we are not sampling at too high an angle
          IF(two * angle(h_tmp, k_tmp, l) > max_angle) GO TO 50
! I(h,k,l)
          h = h_tmp
          k = k_tmp
          i1 = pntint(h, k, l, ok)
          IF(.NOT.ok) GO TO 999
          IF(dosymdump) WRITE(sy,230) h, k, l, i1
! I(-h-k,h,l)
          h = -(h_tmp + k_tmp)
          k = h_tmp
          i2 = pntint(h, k, l, ok)
          IF(.NOT.ok) GO TO 999
          IF(dosymdump) WRITE(sy,230) h, k, l, i2
! I(k,-h-k,l)
          h = k_tmp
          k = -(h_tmp + k_tmp)
          i3 = pntint(h, k, l, ok)
          IF(dosymdump) WRITE(sy,230) h, k, l, i3
! compare intensities
          i_avg = (i1 + i2 + i3) / three
          variance = (ABS(i_avg-i1) + ABS(i_avg-i2) + ABS(i_avg-i3)) / three
! Be careful intensities are not actually zero
          IF(i_avg < tiny_inty) THEN
            tol = tiny_inty
          ELSE
            tol = i_avg * tolerance
            rel_var = variance / i_avg
            IF(rel_var > max_var) max_var = rel_var
          END IF
          match = ABS(i_avg-i1) < tol .AND. ABS(i_avg-i2) < tol  &
              .AND. ABS(i_avg-i3) < tol
          is_good = (i == 1 .OR. is_good) .AND. match
          IF(dosymdump) THEN
            WRITE(sy,240)
            WRITE(sy,250) i_avg
            WRITE(sy,260) variance, hundred * variance / i_avg
!            if(.not.check_sym) then
            WRITE(sy,270) tol
            IF(match) THEN
              WRITE(sy,280) rot_sym
            ELSE
              WRITE(sy,290) rot_sym
            END IF
!            endif
            WRITE(sy,210)
          END IF
        END DO
        tst_rot = is_good
      END IF

      900 IF(dosymdump) THEN
!        if(.not.check_sym) then
        IF(is_good) THEN
          WRITE(sy,300) rot_sym
        ELSE
          WRITE(sy,310) rot_sym
        END IF
!        endif
        WRITE(sy,210)
      END IF
      RETURN

      999 WRITE(op,400) 'ERROR in intensity calculation in TST_ROT'
      WRITE(op,320) '   at h,k,l = ', h,',', k,',', l
      RETURN
      200 FORMAT(1X, 'Testing for a ', i1, '-fold axis')
      210 FORMAT(' ')
      220 FORMAT(1X, '  h', 5X, 'k', 7X, 'l', 20X, 'Intensity')
      230 FORMAT(1X, i3, 3X, i3, 2X, f9.4, 5X, f22.6)
      240 FORMAT(26X, '----------------------')
      250 FORMAT(6X, 'Average Intensity = ', f22.6)
      260 FORMAT(1X, '  Average variation = +/-',f22.6,'  (+/-',g9.2,'%)')
      270 FORMAT(1X, 'Intensity tolerance = +/-', f22.6)
      280 FORMAT(1X, 'Intensities are consistent with a ', i1, '-fold')
      290 FORMAT(1X, 'Intensities are not consistent with a ', i1, '-fold')
      300 FORMAT(1X, 'INTENSITY DISTRIBUTION HAS A ', i1, '-FOLD AXIS')
      310 FORMAT(1X, 'INTENSITY DISTRIBUTION HAS NO ', i1, '-FOLD AXIS')
      320 FORMAT(1X, a, i3, a, i3, a, f7.2)
      330 FORMAT(1X, a, i1, a)
      400 FORMAT(1X, a)
      END FUNCTION tst_rot


! ______________________________________________________________________
! Title: WRTSAD
! Author: MMJT
! Date: 29 Oct 1989; 1st Mar 2000
! Description: This subroutine writes the selected area diffraction
! pattern (SADP) data to a binary file.

!      ARGUMENTS:
!            outfile  -  The name of the file that the binary
!                        selected area diffraction data is to be
!                        written to. (input).
!            view     -  Choice of beam direction. (input).
!                              1  =  normal to the plane k = 0
!                              2  =  normal to the plane h = 0
!                              3  =  normal to the plane h = k
!                              4  =  normal to the plane h = -k
!            l_upper  -  Upper limit of l. (input).
!            hk_lim   -  Upper limit of h (or k). (input).
!            ok       -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  a0, b0, c0, d0, sadblock, scaleint, has_l_mirror

!        modifies:  spec
! ______________________________________________________________________

      SUBROUTINE wrtsad(outfile, view, l_upper, hk_lim, ok)


!     Utiliza las variables: par.f90, inc.f90' , outfile ,view , l_upper, hk_lim ,ok
!                            i, j, p1, p2, length ,rowint(sadsize), incr, sigma, x
!     Utiliza las funciones:length      externa
!     Utiliza las subrutinas:  SMUDGE


      CHARACTER (LEN=*), INTENT(IN OUT)        :: outfile
      INTEGER*4, INTENT(IN OUT)                :: view
      REAL*8, INTENT(IN OUT)                   :: l_upper
      INTEGER*4, INTENT(IN OUT)                :: hk_lim
      LOGICAL, INTENT(IN OUT)                  :: ok

      INTEGER*4 i, j, p1, p2, irow(sadsize)
      REAL*8 rowint(sadsize), incr, sigma, x
      integer(kind=2), dimension(sadsize) :: short_data
      integer(kind=2),          parameter :: zero_2bytes=0

! external subroutine (Some compilers need them declared external)
!      external SMUDGE

      WRITE(op,102) 'Writing SADP data to file ''',  &
          outfile(1:length(outfile)), '''. . .'

! set standard deviation
      sigma = zero

! Establish scaling in reciprocal space, the number of pixels per unit
! first get number of pixels per l
      incr = DBLE(sadsize/2) / l_upper
      IF(view == 1) THEN
! number of pixels per unit h
        incr = incr * SQRT(a0 / c0)
      ELSE IF(view == 2) THEN
! number of pixels per unit k
        incr = incr * SQRT(b0 / c0)
      ELSE IF(view == 3) THEN
! number of pixels per unit along (h = k)
        incr = incr * SQRT((a0 + b0 + d0) / c0)
      ELSE IF(view == 4) THEN
! number of pixels per unit along (h = -k)
        incr = incr * SQRT((a0 + b0 - d0) / c0)
      END IF

      DO  j = sadblock - 1, 0, -1
        DO  i = 1, sadsize
          rowint(i) = zero
        END DO
        DO  i = 0, hk_lim
          p1 = sadsize/2 + nint(i*incr) + 1
          p2 = sadsize/2 - nint(i*incr) + 1

! cycle if we have gone out of bounds
          IF(p1 > sadsize .OR. p1 < 0) THEN
            CYCLE
          END IF
          IF(p2 > sadsize .OR. p2 < 0) THEN
            CYCLE
          END IF

          x = spec(i*sadblock + j + 1) * scaleint
! handle saturated pixels
          IF(x > maxsad) x = maxsad
          rowint(p1) = x
          IF(has_l_mirror) THEN
            rowint(p2) = x
          ELSE
            x = spec(i*sadblock + sadblock - j) * scaleint
! handle saturated pixels
            IF(x > maxsad) x = maxsad
            rowint(p2) = x
          END IF
        END DO
        CALL smudge(rowint, sadsize, sigma, ok)
        IF(.NOT. ok) GO TO 999
        IF(bitdepth == 8) THEN
          WRITE(sa,ERR=900) (CHAR(nint(rowint(i))), i = 1, sadsize)
        ELSE
          do  i = 1, SADSIZE
            irow(i)=nint(rowint(i))
            if( irow(i) > 65535) then
              short_data(i) = -1
            else if( irow(i) > 32767) then
              short_data(i) = irow(i) - 65536
            else
              short_data(i) = irow(i)
            end if
          end do
          write(sa,err=900) (short_data(i), i = 1, SADSIZE)
        END IF
      END DO

! Repeat, in reverse for the bottom half if data had a mirror.
      IF(has_l_mirror) THEN
        DO  j = 1, sadblock - 1
          DO  i = 1, sadsize
            rowint(i) = zero
          END DO
          DO  i = 0, hk_lim
            p1 = sadsize/2 + nint(i*incr) + 1
            p2 = sadsize/2 - nint(i*incr) + 1

! cycle if we have gone out of bounds
            IF(p1 > sadsize .OR. p1 < 0) THEN
              CYCLE
            END IF
            IF(p2 > sadsize .OR. p2 < 0) THEN
              CYCLE
            END IF

            x = spec(i*sadblock + j + 1) * scaleint
! handle saturated pixels
            IF(x > maxsad) x = maxsad
            rowint(p1) = x
            rowint(p2) = x
          END DO
          CALL smudge(rowint, sadsize, sigma, ok)
          IF(.NOT. ok) GO TO 999
          IF(bitdepth == 8) THEN
            WRITE(sa,ERR=900) (CHAR(nint(rowint(i))), i = 1, sadsize)
          ELSE
            do  i = 1, SADSIZE
              irow(i)=nint(rowint(i))
              if( irow(i) > 65535) then
                short_data(i) = -1
              else if( irow(i) > 32767) then
                short_data(i) = irow(i) - 65536
              else
                short_data(i) = irow(i)
              end if
            end do
            write(sa,err=900) (short_data(i), i = 1, SADSIZE)
          END IF
        END DO
! write a blank last line to make the array SADSIZE x SADSIZE square
        IF(bitdepth == 8) THEN
          WRITE(sa,ERR=900) (CHAR(0), i = 1, sadsize)
        ELSE
          write(sa,err=900) (zero_2bytes, i = 1, 2*SADSIZE)
        END IF
      END IF

      IF(sa /= op) CLOSE(UNIT = sa, ERR = 990)

      WRITE(op,103) sadsize, ' x ', sadsize, ' pixels: ', bitdepth, ' bits deep.'

      RETURN
      900 WRITE(op,100) 'ERROR: problems writing to binary SADP file.'
      ok = .false.
      RETURN
      990 WRITE(op,100) 'ERROR: problems closing binary SADP file.'
      ok = .false.
      RETURN
      999 WRITE(op,100) 'ERROR: SMUDGE returned error to WRTSAD.'
      WRITE(op,101) i, j
      RETURN
      100 FORMAT(1X, a)
      101 FORMAT(5X, 'with local variables i = ', i5, ' j = ', i5)
      102 FORMAT(1X, 3A)
      103 FORMAT(1X, i4, a, i4, a, i2, a)
      END SUBROUTINE wrtsad

! ______________________________________________________________________
! Title: WRTSPC
! Author: MWD and MMJT
! Date: 18 Aug 1988; 7 Mar 1995
! Description:  This routine writes the spectrum arrays to file.

!      ARGUMENTS:
!            spcfile  -  The name of the output data file. (input).
!            ok       -  logical flag indicating all went well.
!                                                      (output).

!      COMMON VARIABLES:
!            uses:  th2_min, th2_max, d_theta, spec, brd_spc, blurring
!                   NONE, RAD2DEG

!        modifies:  no COMMON variables are modified
! ______________________________________________________________________

      SUBROUTINE wrtspc(spcfile, ok)

!     Utiliza las variables: par.f90, inc.f90' , spcfile ,ok ,tab ,i, n_low, n_high, length, theta
!     Utiliza las funciones:length      externa
!     Utiliza las subrutinas:

      CHARACTER (LEN=*), INTENT(IN OUT)        :: spcfile
      LOGICAL, INTENT(IN OUT)                  :: ok
      CHARACTER (LEN=1) :: tab
      INTEGER*4 i, n_low
      REAL*8 theta

      n_low  = 1
      tab = CHAR(9)
      theta = th2_min * rad2deg
      WRITE(op,200) 'Writing spectrum data to file ''',      spcfile(1:length(spcfile)), '''. . .'
!      do 10 i = int(HALF*th2_min / d_theta) + 1,
!     |                  int(HALF*th2_max / d_theta) + 1
      IF(blurring == NONE) THEN
        DO  i = n_low, n_high
         ! WRITE(sp,101,ERR=50) theta, tab, spec(i)
      !    WRITE(unit_sp,"(1X, e12.5, a, g13.6,a)",ERR=50) theta, "    ", spec(i),"    1"
          theta = theta + two * rad2deg * d_theta
        END DO
      ELSE
        DO  i = n_low, n_high
         ! WRITE(sp,100,ERR=50) theta, tab, spec(i), tab, brd_spc(i)
     !     WRITE(unit_sp,"(1X, e12.5, 2(a, g13.6),a)",ERR=50) &
     !              theta, "    ", spec(i), "   ", brd_spc(i),"    1"
          theta = theta + two * rad2deg * d_theta
        END DO
      END IF

    !  IF(unit_sp /= op) CLOSE(sp,ERR=60)
      WRITE(op,202) 'Spectrum written.'
      WRITE(op,202)
      RETURN
      50 WRITE(op,202) 'Unable to write to spectrum file.'
    !  IF(unit_sp /= op) CLOSE(sp,ERR=60)
      ok = .false.
      RETURN
      60 WRITE(op,202) 'Unable to close spectrum file.'
      ok = .false.
      RETURN
      100 FORMAT(1X, e12.5, 2(a, g13.6))
      101 FORMAT(1X, e12.5, a, g13.6)
      200 FORMAT(1X, 3A)
      201 FORMAT(1X, a, i2, 2A)
      202 FORMAT(1X, a)
      END SUBROUTINE wrtspc

!!S XYfiles
! ______________________________________________________________________
! Title: XYPHSE
! Author: MMJT
! Date: 8 June 1990
! Description:  This routine pre-calculates the h and k components of
! phases for each atom. This is called when l = 0, since h and k
! are held constant while l is varied along each streak.

!      ARGUMENTS:
!            h  -  reciprocal lattice vector h-component. (input).
!            k  -  reciprocal lattice vector k-component. (input).

!      COMMON VARIABLES:
!            uses:  n_actual, l_n_atoms, a_pos, l_actual, n_layers

!        modifies:  hx_ky
! ______________________________________________________________________

      Subroutine Xyphse(H, K)
!     Utiliza las variables: par.f90, inc.f90' ,h,k, i, m
!     Utiliza las funciones:
!     Utiliza las subrutinas:

      Integer, Intent(In)  :: h
      Integer, Intent(In)  :: k

      Integer :: i, m

      Do  m = 1, n_actual
        Do  i = 1, l_n_atoms(m)
          hx_ky(i,m) = h*a_pos(1,i,m) + k*a_pos(2,i,m)
        End Do
      End Do

      Return
      End Subroutine xyphse

! ______________________________________________________________________
! Title: YRDSTK
! Author: MMJT
! Date: 12 Aug 1989
! Description: YRDSTK checks that y is a multiple of x.

!      ARGUMENTS:
!            x   -  reference value. (input).
!            y   -  value to be tested. (input)
!            ok  -  logical flag set to .false. if y is zero. (output).

!      YRDSTK returns logical .true. if y is a multiple of x.
! ______________________________________________________________________

      LOGICAL FUNCTION yrdstk(x, y, ok)
!       Utiliza las variables: par.f90, x,y, ok, tmp
!       Utiliza las funciones:
!       Utiliza las subrutinas:

        REAL*8, INTENT(IN OUT)      :: x
        REAL*8, INTENT(IN)          :: y
        LOGICAL, INTENT(IN OUT)     :: ok
        REAL*8 :: tmp

        yrdstk = .false.
        IF(y == zero) then
         ok = .false.
         RETURN
        END IF
        IF(x == zero) THEN
        ! This is the first visit to YRDSTK. Set things up.
          x = y
          yrdstk = .true.
        ELSE
          tmp = y / x
          IF(ABS(nint(tmp)-tmp) <= eps3*tmp) yrdstk = .true.
        END IF
        RETURN
      END FUNCTION yrdstk

  End Module diffax_calc
