!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!----
!!----
!!---- MODULE: CFML_Optimization_LSQ
!!----   INFO: Module implementing several algorithms for non-linear least-squares.
!!----         At present only two versions of the Levenberg-Marquardt method are implemented.
!!--..
!!--..         General Information Concerning The Least-Squares Procedures
!!--..
!!--..         There are two high level procedures contained in CFML_Optimization_LSQ based
!!--..         in the Levenberg-Marquardt method. The first procedure, called "Marquardt_Fit",
!!--..         is a simple implementation of the method and the second one is a Fortran 90
!!--..         version of the MINPACK Fortran 77 LMxxx subroutines, accessible through the
!!--..         general name "Levenberg_Marquardt_Fit". The second one is, in principle, more
!!--..         robust for general LSQ problems.
!!--..
!!--..
!!--..         Marquardt_Fit (overloaded subroutine)
!!--..         =====================================
!!--..         This procedure needs analytical derivatives to be provided by the user. This
!!--..         procedure minimizes the following cost function:
!!--..
!!--..         C =  Sum(i) { w(i) ( y(x(i)) - yc(x(i),a))^2 }  =
!!--..           =  Sum(i) {  [( y(x(i)) - yc(x(i),a))/sig(i)]^2 }
!!--..
!!--..         Where x(i),y(i),w(i) are observational values (normally w(i)=1/variance(y(i)))
!!--..         and yc is a model to the data depending on the vector a containing the
!!--..         parameters of the model. The user should provide a subroutine for calculating
!!--..         yc(x(i),a) and invoke the procedure passing as the first argument the name
!!--..         of the subroutine.
!!--..
!!--..         Two interfaces are available for calling this procedure:
!!--..
!!--..         First interface:
!!--..         ----------------
!!--..         Subroutine Marquardt_Fit(Model_Functn, X, Y, W, Yc, Nobs, c, vs, Ipr, Chi2, scroll_lines)
!!--..            real(kind=cp), dimension(:),   intent(in)             :: x      !Vector of x-values
!!--..            real(kind=cp), dimension(:),   intent(in)             :: y      !Vector of observed y-values
!!--..            real(kind=cp), dimension(:),   intent(in out)         :: w      !Vector of weights-values (1/variance)
!!--..            real(kind=cp), dimension(:),   intent(   out)         :: yc     !Vector of calculated y-values
!!--..            integer                    ,   intent(in)             :: nobs   !Number of effective components of x,y,w,yc
!!--..            type(LSQ_conditions_type),     intent(in out)         :: c      !Conditions for the algorithm
!!--..            type(LSQ_State_Vector_type),   intent(in out)         :: vs     !State vector for the model calculation
!!--..            integer                    ,   intent(in)             :: Ipr    !Logical unit for writing
!!--..            real(kind=cp),                 intent(out)            :: chi2   !Reduced Chi-2
!!--..            character(len=*),dimension(:), intent(out), optional  :: scroll_lines  !If present, part of the output is stored
!!--..                                                                            !in this text for treatment in the calling program
!!--..
!!--..            Model_functn                                                    !Name of the subroutine calculating yc(i) for point x(i)
!!--..            Interface                                      !Interface for the Model_Functn subroutine
!!--..              Subroutine Model_Functn(iv,Xv,ycalc,a,der)
!!--..                   use CFML_GlobalDeps, only: cp
!!--..                   integer,                             intent(in) :: iv     !Number of the component: "i" in x(i)
!!--..                   real(kind=cp),                       intent(in) :: xv     !Value of x(i)
!!--..                   real(kind=cp),                       intent(out):: ycalc  !Value of yc at point x(i) => ycalc=yc(i)
!!--..                   real(kind=cp),dimension(:),          intent(in) :: a      !Vector of free parameters
!!--..                   real(kind=cp),dimension(:),optional, intent(out):: der    !Derivatives of the function w.r.t. free parameters
!!--..              End Subroutine Model_Functn                                    !at the given point
!!--..           End Interface
!!--..
!!--..
!!--..           The model function subroutine calculates the value of the function "ycalc"
!!--..           at the provided single point "Xv" (component i of vector X).It contains
!!--..           as an argument the vector "a" of "free parameters", whereas the total number
!!--..           of parameters in the model is stored in the structure "Vs" of type
!!--..           LSQ_State_Vector_type (see the definition of the type for details). It is
!!--..           responsibility of the user to write the Model_Functn subroutine with access
!!--..           to Vs (via a module variable) and do the necessary connection between the
!!--..           "a" and "Vs%pv" vectors, as well as the calculation of the derivatives of
!!--..           "ycalc" w.r.t "a" stored in the optional array "der"
!!--..
!!--..         Second interface:
!!--..         -----------------
!!--..         Subroutine Marquardt_Fit(Model_Functn, d, c, vs, Ipr, Chi2, scroll_lines)
!!--..            Type(LSQ_Data_Type),        intent(in out) :: d      !Data type for LSQ
!!--..            type(LSQ_conditions_type),  intent(in out) :: c      !Conditions for the algorithm
!!--..            type(LSQ_State_Vector_type),intent(in out) :: vs     !State vector for the model calculation
!!--..            integer                    ,intent(in)     :: Ipr    !Logical unit for writing
!!--..            real(kind=cp),              intent(out)    :: chi2   !Reduced Chi-2
!!--..            character(len=*),dimension(:), intent(out), optional  :: scroll_lines  !If present, part of the output is stored
!!--..                                                                                   !in this text for treatment in the calling program
!!--..
!!--..           Model_functn                                            !Name of the subroutine calculating yc(i) for point x(i)
!!--..           Interface                                               !Interface for the Model_Functn subroutine
!!--..               Subroutine Model_Functn(iv,Xv,ycalc,Vsa,calder)
!!--..                    use CFML_GlobalDeps, only: cp
!!--..                    integer,                    intent(in)     :: iv     !Number of the component: "i" in x(i)
!!--..                    real(kind=cp),              intent(in)     :: xv     !Value of x(i)
!!--..                    real(kind=cp),              intent(out)    :: ycalc  !Value of yc at point x(i) => ycalc=yc(i)
!!--..                    Type(LSQ_State_Vector_type),intent(in out) :: Vsa    !Parameters, codes, and derivatives
!!--..                    logical,optional,           intent(in)     :: calder !If present, derivatives, stored in Vsa%dpv, are calculated
!!--..               End Subroutine Model_Functn
!!--..           End Interface
!!--..
!!--..           This second interface is simpler than the first one, it encapsulates the
!!--..           handling of the connection between free and total number of parameters.
!!--..           The structure "d", of type LSQ_Data_Type, contains the information that
!!--..           was explicitly used in the first interface (X, Y, W, Yc, Nobs) in five
!!--..           arguments.
!!--..           The model function uses directly as an argument the object "Vsa", of type
!!--..           LSQ_State_Vector_type in terms of which can be explicitly written the
!!--..           calculation of the model function "ycalc" at the given point "iv,Xv". The
!!--..           calculation of the derivatives (stored in Vsa%dpv vector) are needed only
!!--..           if the argument "calder" is present (irrespective of its value .true. or
!!--..           .false.).
!!--..
!!--..
!!--..         Levenberg_Marquardt_Fit (overloaded subroutine), Based on MINPACK
!!--..         =================================================================
!!--..         This procedure needs either analytical (provided by the user) or numerical
!!--..         derivatives (automatically calculated by the procedure from the provided
!!--..         model function subroutine).The procedure minimizes a cost function of the
!!--..         form:
!!--..
!!--..             C =  Sum(i) { fvec(i)^2 }  =   Sum(i) {  [( y(i) - yc(i,a))/sig(i)]^2 }
!!--..
!!--..         The user-provided subroutine should have as arguments the vector of free
!!--..         parameters x(1:n) and the vector or residuals fvec(m).
!!--..
!!--..         Interface No_Fderivatives
!!--..           Subroutine Model_Functn(m, n, x, fvec, iflag)             !Model Function subroutine
!!--..             Use CFML_GlobalDeps, Only: cp
!!--..             Integer,                       Intent(In)    :: m, n    !Number of observations and free parameters
!!--..             Real (Kind=cp),Dimension(:),   Intent(In)    :: x       !Array with the values of free parameters: x(1:n)
!!--..             Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec    !Array of residuals fvec=(y-yc)/sig : fvec(1:m)
!!--..             Integer,                       Intent(In Out):: iflag   !If iflag=1 calculate only fvec without changing fjac
!!--..           End Subroutine Model_Functn                               !If iflag=2 calculate only fjac keeping fvec fixed
!!--..         End Interface No_Fderivatives
!!--..
!!--..
!!--..         If analytical derivatives can be calculated the subroutine should have also
!!--..         as an additional argument the Jacobian rectangular matrix fjac(m,n)
!!--..
!!--..         Interface  FDerivatives
!!--..           Subroutine Model_Functn(m, n, x, fvec, fjac, iflag)       !Model Function subroutine
!!--..             Use CFML_GlobalDeps, Only: cp
!!--..             Integer,                       Intent(In)    :: m, n    !Number of observations and free parameters
!!--..             Real (Kind=cp),Dimension(:),   Intent(In)    :: x       !Array with the values of free parameters: x(1:n)
!!--..             Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec    !Array of residuals fvec=(y-yc)/sig : fvec(1:m)
!!--..             Real (Kind=cp),Dimension(:,:), Intent(Out)   :: fjac    !Jacobian Dfvec/Dx(i,j) = [ dfvec(i)/dx(j) ] : fjac(1:m,1:n)
!!--..             Integer,                       Intent(In Out):: iflag   !If iflag=1 calculate only fvec without changing fjac
!!--..           End Subroutine Model_Functn                               !If iflag=2 calculate only fjac keeping fvec fixed
!!--..         End Interface FDerivatives
!!--..
!!--..
!!--..         The five available interfaces are related to these two options:
!!--..
!!--..         First/second Interface(s) (no derivatives required)
!!--..         ---------------------------------------------------
!!--..
!!--..         Subroutine Levenberg_Marquard_Fit(Model_Functn, m, c, Vs, chi2, infout,residuals)
!!--..            Integer,                     Intent(In)      :: m        !Number of observations
!!--..            type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
!!--..            type(LSQ_State_Vector_type), Intent(In Out)  :: Vs       !State vector
!!--..            Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
!!--..            character(len=*),            Intent(out)     :: infout   !Information about the refinement (min length 256)
!!--..            Real (Kind=cp), dimension(:),optional, intent(out) :: residuals
!!--..         End Subroutine
!!--..
!!--..
!!--..         The second interface corresponds to the original version in MINPACK
!!--..         (called lmdif1)
!!--..
!!--..         Subroutine Levenberg_Marquard_Fit(Model_Functn, m, n, x, fvec, tol, info, iwa)
!!--..            Integer,                     Intent(In)      :: m
!!--..            Integer,                     Intent(In)      :: n
!!--..            Real (Kind=cp),Dimension(:), Intent(In Out)  :: x
!!--..            Real (Kind=cp),Dimension(:), Intent(Out)     :: fvec
!!--..            Real (Kind=cp),              Intent(In)      :: tol
!!--..            Integer,                     Intent(Out)     :: info
!!--..            Integer,Dimension(:),        Intent(Out)     :: iwa
!!--..         End Subroutine
!!--..
!!--..         The interface for Model_Functn is identical to the interface No_Fderivatives
!!--..
!!--..
!!--..         Third/fourth/fifth Interface(s) (analytical derivatives required)
!!--..         -----------------------------------------------------------
!!--..
!!--..         Subroutine Levenberg_Marquard_Fit(Model_Functn, m, c, Vs, chi2, calder, infout,residuals)
!!--..            Integer,                     Intent(In)      :: m        !Number of observations
!!--..            type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
!!--..            type(LSQ_State_Vector_type), Intent(In Out)  :: Vs       !State vector
!!--..            Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
!!--..            logical,                     Intent(in)      :: calder   !logical (should be .true.) used only for purposes
!!--..                                                                     !of making unambiguous the generic procedure
!!--..            character(len=*),            Intent(out)     :: infout   !Information about the refinement (min length 256)
!!--..
!!--..         This third interface is effectively called LM_DER
!!--..
!!--..         The fourth interface corresponds to the original version in MINPACK
!!--..         (called lmder1)
!!--..
!!--..        Subroutine lmder1(Model_Functn, m, n, x, fvec, fjac, tol, info, ipvt)
!!--..          Integer,                        Intent(In)      :: m
!!--..          Integer,                        Intent(In)      :: n
!!--..          Real (Kind=cp),Dimension(:),    Intent(In Out)  :: x
!!--..          Real (Kind=cp),Dimension(:),    Intent(Out)     :: fvec
!!--..          Real (Kind=cp),Dimension(:,:),  Intent(In Out)  :: fjac
!!--..          Real (Kind=cp),                 Intent(In)      :: tol
!!--..          Integer,                        Intent(Out)     :: info
!!--..          Integer,Dimension(:),           Intent(In Out)  :: ipvt
!!--..
!!--..          The interface for Model_Functn is identical to the previous interface
!!--..          Fderivatives
!!--..
!!--..         The fifth interface corresponds to a mixed version that uses the
!!--..         conditions of refinement as the third interface, but the state vector
!!--..         is directly an 1D array. This cal use directly ve state vector as
!!--..         defined in CFML_Refcodes
!!--..         This fifth interface is effectively called LM_DERV
!!--..
!!--..         Subroutine Levenberg_Marquard_Fit(Model_Functn, m, c, Vect, chi2, calder, infout,residuals)
!!--..            Integer,                     Intent(In)      :: m        !Number of observations
!!--..            type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
!!--..            real, dimension(:),          Intent(In Out)  :: Vect     !State vector like in CFML_Refcodes V_Vec
!!--..            Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
!!--..            logical,                     Intent(in)      :: calder   !logical (should be .true.) used only for purposes
!!--..                                                                     !of making unambiguous the generic procedure
!!--..            character(len=*),            Intent(out)     :: infout   !Information about the refinement (min length 256)
!!--..
!!--..          The interface for Model_Functn is identical to the previous interfaces
!!--..          with derivatives
!!--..
!!--..
!!---- HISTORY
!!----    Update: 04/03/2011
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,    only : Cp, Dp
!!--++    Use CFML_Math_General,  only : Invert_Matrix, enorm => Euclidean_Norm
!!--++    Use CFML_LSQ_TypeDef
!!----
!!---- VARIABLES
!!--++    CH                      [Private]
!!--++    CORREL                  [Private]
!!--++    CURV_MAT                [Private]
!!----    ERR_LSQ
!!----    ERR_LSQ_MESS
!!----    INFO_LSQ_MESS
!!--++    NAMFREE                 [Private]
!!--++    PN                      [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       FCHISQ
!!----
!!----    Subroutines:
!!--++       BOX_CONSTRAINTS           [Private]
!!--++       CURFIT_v1                 [Private]
!!--++       CURFIT_v2                 [Private]
!!--++       FDJAC2                    [Private]
!!----       LEVENBERG_MARQUARDT_FIT   [Overloaded]   {LM_DER, LM_DERV, LM_DIF, LMDER1, LMDIF1}
!!--++       LM_DER                    [Private]
!!--++       LM_DERV                   [Private]
!!--++       LMDER                     [Private]
!!--++       LM_DIF                    [Private]
!!--++       LMDIF                     [Private]
!!--++       LMPAR                     [Private]
!!----       MARQUARDT_FIT             [Overloaded]  {MARQUARDT_FIT_v1, MARQUARDT_FIT_v2}
!!--++       MARQUARDT_FIT_v1          [Private]
!!--++       MARQUARDT_FIT_v2          [Private]
!!--++       OUTPUT_CYC                [Private]
!!----       INFO_LSQ_LM               [Overloaded] {INFO_LSQ_LM_V, INFO_LSQ_LM_VS}
!!--++       INFO_LSQ_LM_V             [Private]
!!--++       INFO_LSQ_LM_VS            [Private]
!!----       INFO_LSQ_OUTPUT
!!--++       QRFAC                     [Private]
!!--++       QRSOLV                    [Private]
!!----
!!----
!!
 Module CFML_Optimization_LSQ
    !---- Use Files ----!
    Use CFML_GlobalDeps,   only: Cp, Dp
    Use CFML_LSQ_TypeDef
    Use CFML_Math_General, only: Invert_Matrix, enorm => Euclidean_Norm, &
                                 Upper_Triangular,SVDcmp

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public  :: fchisq
    !---- List of public overloaded procedures: functions ----!

    !---- List of public subroutines ----!
    public  :: Marquardt_Fit, Levenberg_Marquardt_Fit, Info_LSQ_Output, Info_LSQ_LM

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!

    !---- List of private subroutines ----!
    private :: marquardt_fit_v1,marquardt_fit_v2,curfit_v1,curfit_v2,box_constraints,output_cyc,&
               LM_Dif, lmdif, LM_Der, LM_DerV, lmder, fdjac2, lmpar, qrfac, qrsolv, lmder1, lmdif1

    !---- Definitions ----!

    !!--++
    !!--++ CH
    !!--++    real(kind=cp),     dimension(Max_Free_Par), private   :: ch
    !!--++
    !!--++    (PRIVATE)
    !!--++    Vector holding the change in the values of parameters (ch = pn - pv)
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), dimension(Max_Free_Par), private   :: ch         ! ch = pn - pv

    !!--++
    !!--++ CORREL
    !!--++    real(kind=cp), dimension(Max_Free_Par,Max_Free_Par), public  :: correl
    !!--++
    !!--++    Variance/covariance/correlation matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), dimension(Max_Free_Par,Max_Free_Par), private  :: correl     !Variance/covariance/correlation matrix


    !!--++
    !!--++ CURV_MAT
    !!--++    real(kind=cp), dimension(Max_Free_Par,Max_Free_Par), public  :: curv_mat
    !!--++
    !!--++    Curvature matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), dimension(Max_Free_Par,Max_Free_Par), private  :: curv_mat   !Curvature matrix

    !!----
    !!---- ERR_LSQ
    !!----    logical, public  :: ERR_LSQ
    !!----
    !!----    Logical variable. The vaule .true. indicates that an error condition occurs
    !!----
    !!---- Update: February - 2005
    !!
    logical, public   :: ERR_Lsq =.false.

    !!----
    !!---- ERR_LSQ_MESS
    !!----    Character(len=150), public  :: Err_Lsq_Mess
    !!----
    !!----    Character variable containing the error message associated to the
    !!----    ocurrence of an error condition
    !!----
    !!---- Update: February - 2005
    !!
    Character(len=150), public  :: ERR_Lsq_Mess

    !!----
    !!---- INFO_LSQ_MESS
    !!----    Character(len=150), public  :: Info_Lsq_Mess
    !!----
    !!----    Character variable containing the information message associated to the
    !!----    exit parameter "info" of the Levenberg_Marquardt_Fit procedure.
    !!----
    !!---- Update: August - 2005
    !!
    Character(len=150), public  :: Info_Lsq_Mess


    !!--++
    !!--++ NAMFREE
    !!--++    character(len=40), dimension(Max_Free_Par), private   :: namfree
    !!--++
    !!--++  (PRIVATE)
    !!--++  Names of refined parameters
    !!--++
    !!--++ Update: January - 2014
    !!
    character(len=40), dimension(Max_Free_Par), private,save   :: namfree    !Names of refined parameters

    !!--++
    !!--++ PN
    !!--++    real(kind=cp), dimension(Max_Free_Par), public   :: pn
    !!--++
    !!--++    Vector with new values of parameters
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), dimension(Max_Free_Par), private   :: pn   !Vector with new values of parameters

    !---- Interfaces ----!

    Interface Marquardt_fit
      Module Procedure Marquardt_Fit_v1
      Module Procedure Marquardt_Fit_v2
    End Interface

    Interface Levenberg_Marquardt_Fit
      Module Procedure LM_Dif
      Module Procedure lmdif1
      Module Procedure LM_Der
      Module Procedure LM_DerV
      Module Procedure lmder1
    End Interface Levenberg_Marquardt_Fit

    Interface Info_LSQ_LM
      Module Procedure Info_LSQ_LM_VS
      Module Procedure Info_LSQ_LM_V
    End Interface Info_LSQ_LM

 Contains

    !---- Functions ----!

    !!----
    !!---- Function Fchisq(Nfr,Nobs,Y,W,Yc) Result(Chisq)
    !!----    integer, intent(in)                    :: Nfr
    !!----    integer, intent(in)                    :: Nobs
    !!----    real(kind=cp), intent(in),dimension(:) :: y
    !!----    real(kind=cp), intent(in),dimension(:) :: w
    !!----    real(kind=cp), intent(in),dimension(:) :: yc
    !!----
    !!----    Evaluate reduced chi2
    !!----
    !!---- Update: February - 2005
    !!
    Function Fchisq(Nfr,Nobs,Y,W,Yc) Result(Chisq)
       !---- Arguments ----!
       integer,                    intent(in) :: nfr
       integer,                    intent(in) :: nobs
       real(kind=cp),dimension(:), intent(in) :: y
       real(kind=cp),dimension(:), intent(in) :: w
       real(kind=cp),dimension(:), intent(in) :: yc
       real(kind=cp)                          :: chisq

       !---- Local variables ----!
       integer :: i

       chisq=0.0
       do i=1,nobs
          chisq=chisq+w(i)*(y(i)-yc(i))**2
       end do
       chisq=chisq/real(nfr)

       return
    End Function Fchisq

    !---- Subroutines ----!

    !!--++
    !!--++ Subroutine Box_Constraints(A,Sa,Fixed,c,vs)
    !!--++    real(kind=cp), dimension (Max_Free_Par), intent(in out) :: a
    !!--++    real(kind=cp), dimension (Max_Free_Par), intent(in out) :: sa
    !!--++    logical                        ,         intent(   out) :: fixed
    !!--++    Type(LSQ_Conditions_type),               intent(in)     :: c     !conditions for refinement
    !!--++    Type(LSQ_State_Vector_type),             intent(in)     :: vs    !State vector
    !!--++
    !!--++    (PRIVATE)
    !!--++    This subroutine avoid a peak-position parameter undergoing a change
    !!--++    greater than percent% w.r.t. the initial value. This is probably enough
    !!--++    to avoid lack of convergence when a peak makes an excursion outside
    !!--++    the current angular range. The parameter is fixed for the next cycle
    !!--++    The intensity parameters are also constrained to be strictly positive
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Box_Constraints(A,Sa,Fixed,c,vs)
       !---- Arguments ----!
       real(kind=cp), dimension (Max_Free_Par), intent(in out) :: a
       real(kind=cp), dimension (Max_Free_Par), intent(in out) :: sa
       logical                        ,         intent(   out) :: fixed
       Type(LSQ_Conditions_type),               intent(in)     :: c     !conditions for refinement
       Type(LSQ_State_Vector_type),             intent(in)     :: vs    !State vector

       !---- Local variables ----!
       integer       :: i,ncount
       real(kind=cp) :: per
       logical       :: ifi

       fixed=.false.
       ncount=0
       per=c%percent*0.01
       do i=1,vs%np
          ifi=.false.
          if (vs%code(i) /= 0) then
             ncount=ncount+1
             if (abs(vs%pv(i)) > 0.0000001) then
                if (abs(a(ncount)-vs%pv(i))/vs%pv(i) > per )  ifi=.true.
                if (ifi) then
                   fixed=.true.
                   a(ncount) = vs%pv(i)
                   namfree(ncount) = trim(namfree(ncount))//"-*"
                   sa(ncount) = 0.0
                   pn(i) = vs%pv(i)
                   ch(i) = 0.0
                end if
             end if
          end if
       end do

       return
    End Subroutine Box_Constraints

    !!--++
    !!--++  Subroutine Curfit_v1(Model_Functn, X, Y, W, Nobs, c, A, Sa, Fl, Yc, Chir, Ifail)
    !!--++     real(kind=cp),    dimension(:),      intent(in)      :: x     !vector with abcisae
    !!--++     real(kind=cp),    dimension(:),      intent(in)      :: y     !Observed values
    !!--++     real(kind=cp),    dimension(:),      intent(in out)  :: w     !weight of observations
    !!--++     integer,                             intent(in)      :: nobs  !number of observations
    !!--++     Type(LSQ_Conditions_type),           intent(in)      :: c     !conditions for refinement
    !!--++     real(kind=cp),dimension(:),          intent(in out)  :: a     !vector of parameter
    !!--++     real(kind=cp),dimension(:),          intent(in out)  :: sa    !estimated standard deviations
    !!--++     real(kind=cp),                       intent(in out)  :: fl    !Marquardt LAMBDA value
    !!--++     real(kind=cp),dimension(:),          intent(out)     :: yc    !Calculated
    !!--++     real(kind=cp),                       intent(out)     :: chir
    !!--++     integer,                    intent(out)     :: ifail
    !!--++
    !!--++     Interface
    !!--++      Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
    !!--++         use CFML_GlobalDeps, only: cp
    !!--++         integer,                             intent(in) :: iv
    !!--++         real(kind=cp),                       intent(in) :: xv
    !!--++         real(kind=cp),dimension(:),          intent(in) :: aa
    !!--++         real(kind=cp),                       intent(out):: ycalc
    !!--++         real(kind=cp),dimension(:),optional, intent(out):: der
    !!--++      End Subroutine Model_Functn
    !!--++     End Interface
    !!--++
    !!--++ Update: February - 2003
    !!
    Subroutine Curfit_v1(Model_Functn, X, Y, W, Nobs, c, A, Sa, Fl, Yc, Chir, Ifail,nt)
       !---- Arguments ----!
       real(kind=cp),    dimension(:),      intent(in)      :: x     !vector with abcisae
       real(kind=cp),    dimension(:),      intent(in)      :: y     !Observed values
       real(kind=cp),    dimension(:),      intent(in out)  :: w     !weight of observations
       integer,                             intent(in)      :: nobs  !number of observations
       Type(LSQ_Conditions_type),           intent(in)      :: c     !conditions for refinement
       real(kind=cp),dimension(:),          intent(in out)  :: a     !vector of parameter
       real(kind=cp),dimension(:),          intent(in out)  :: sa    !estimated standard deviations
       real(kind=cp),                       intent(in out)  :: fl    !Marquardt LAMBDA value
       real(kind=cp),dimension(:),          intent(out)     :: yc    !Calculated
       real(kind=cp),                       intent(out)     :: chir
       integer,                             intent(out)     :: ifail,nt

       Interface
        Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
           use CFML_GlobalDeps, only: cp
           integer,                             intent(in) :: iv
           real(kind=cp),                       intent(in) :: xv
           real(kind=cp),dimension(:),          intent(in) :: aa
           real(kind=cp),                       intent(out):: ycalc
           real(kind=cp),dimension(:),optional, intent(out):: der
        End Subroutine Model_Functn
       End Interface

       !---- Local variables ----!
       logical                          :: change_par
       logical                          :: singular
       integer                          :: ntrials
       integer                          :: nfr,j,i,k,ntr
       real(kind=cp)                    :: chisq1
       real(kind=cp), dimension(c%npvar):: beta, b, der


       ntrials=40

       !---- Optimization with MARQUARDT ALGORITHM ----!
       nfr=nobs-c%npvar
       nt=0
       ifail=0
       change_par=.false.
       if (nobs > 2000) ntrials= 20

       !---- Evaluate ALFA and BETA matrices       ----!
       !---- (BETA=DP*ALFA)--- DP=BETA*ALFA(-1)    ----!
       beta(:) =0.0
       curv_mat(:,:)=0.0


       do i=1,nobs
          call model_functn(i,x(i),yc(i),a,der)
          if (c%iw == 1 .and. yc(i) /= 0.0) w(i)=1.0/yc(i)
          do j=1,c%npvar
             beta(j)=beta(j)+w(i)*(y(i)-yc(i))*der(j)
             do k=1,j
                curv_mat(j,k)=curv_mat(j,k)+w(i)*der(j)*der(k)
             end do
          end do
       end do

       do i=1,c%npvar
          if (curv_mat(i,i) < 1.0e-30) then
             Err_lsq =.true.
             write(unit=Err_Lsq_Mess,fmt="(a,i5,a,a)")  &
                  " => Singular matrix!!, problem with parameter no.:",i," -> ",trim(namfree(i))
             ifail=2
             return
          end if
       end do

       do j=1,c%npvar
          do k=1,j
             curv_mat(k,j)=curv_mat(j,k)
          end do
       end do

       !---- Evaluate CHI2 at starting point ----!
       chisq1=fchisq(nfr,nobs,y,w,yc)  !chi2 before inverting and applying shifts

       !---- Invert modified curvature matrix to find new parameters
       do ntr=1,ntrials
          do j=1,c%npvar
             do k=1,c%npvar
                correl(j,k)=curv_mat(j,k)/sqrt(curv_mat(j,j)*curv_mat(k,k))
             end do
             correl(j,j)=1.0+fl
          end do

          call Invert_Matrix(correl(1:c%npvar,1:c%npvar),correl(1:c%npvar,1:c%npvar),singular)

          if (singular) then
             Err_lsq =.true.
             do i=1,c%npvar
                if (correl(i,i) < 1.0e-20) then
                   j=i !iperm(i)
                   exit
                end if
             end do
             write(unit=Err_Lsq_Mess,fmt="(a,i5,a,a)")  &
                  " => Singular matrix!!, problem with parameter no.:",j," -> ",trim(namfree(j))
             ifail=2
             return
          end if

          if (fl < 1.0e-20)  then
             chir=chisq1
             exit
          end if

          do j=1,c%npvar
             b(j)=a(j)
             do k=1,c%npvar
                b(j)=b(j)+beta(k)*correl(j,k)/sqrt(curv_mat(j,j)*curv_mat(k,k))
             end do
          end do


          do i=1,nobs
             call model_functn(i,x(i),yc(i),b)
          end do

          chir=fchisq(nfr,nobs,y,w,yc)

          if (ntr == ntrials .or. fl > 1.0e+20 ) then    !no improvement
             ifail=1
             change_par=.false.
             exit
          end if

          !---- If Chi2 increase, increase fl and try again ----!
          if (chisq1-chir < 0.0) then  !Chi2 increases
             fl=10.0*fl                !Increase fl
             cycle
          end if
          change_par=.true.
          exit
       end do  ! ntrials
       nt=ntr
       !---- Evaluate parameters and uncertainties ----!
       do j=1,c%npvar
          sa(j)=sqrt(abs(correl(j,j)/curv_mat(j,j)))
       end do
       if (change_par) then
          do j=1,c%npvar
             a(j)=b(j)
          end do
       end if
       fl=fl/10.0

       return
    End Subroutine Curfit_v1

    !!--++
    !!--++  Subroutine Curfit_v2(Model_Functn, d, c, vs, Fl, Chir, Ifail)
    !!--++     Type(LSQ_Data_type)                  intent(in out)  :: d     !Data
    !!--++     Type(LSQ_Conditions_type),           intent(in)      :: c     !conditions for refinement
    !!--++     Type(LSQ_State_Vector_type),         intent(in out)  :: vs    !State Vector with model parameters
    !!--++     real(kind=cp),                       intent(in out)  :: fl    !Marquardt LAMBDA value
    !!--++     real(kind=cp),dimension(:),          intent(out)     :: yc    !Calculated
    !!--++     real(kind=cp),                       intent(out)     :: chir
    !!--++     integer,                    intent(out)     :: ifail
    !!--++
    !!--++     Interface
    !!--++      Subroutine Model_Functn(iv,xv,ycalc,Vsa,der)
    !!--++         use CFML_GlobalDeps, only: cp
    !!--++         use CFML_LSQ_TypeDef,  only: LSQ_State_Vector_type
    !!--++         integer,                             intent(in) :: iv
    !!--++         real(kind=cp),                       intent(in) :: xv
    !!--++         real(kind=cp),                       intent(out):: ycalc
    !!--++         Type(LSQ_State_Vector_type),         intent(in) :: Vsa
    !!--++         real(kind=cp),dimension(:),optional, intent(out):: der
    !!--++      End Subroutine Model_Functn
    !!--++     End Interface
    !!--++
    !!--++ Update: February - 2003
    !!
    Subroutine Curfit_v2(Model_Functn, d, c, vs, Fl, Chir, Ifail,nt)
       !---- Arguments ----!
       Type(LSQ_Data_Type),                 intent(in out)  :: d
       Type(LSQ_Conditions_type),           intent(in)      :: c     !conditions for refinement
       Type(LSQ_State_Vector_type),         intent(in out)  :: vs
       real(kind=cp),                       intent(in out)  :: fl    !Marquardt LAMBDA value
       real(kind=cp),                       intent(out)     :: chir
       integer,                             intent(out)     :: ifail,nt

       Interface
        Subroutine Model_Functn(iv,xv,ycalc,Vsa,calder)
           use CFML_GlobalDeps, only: cp
           use CFML_LSQ_TypeDef,  only: LSQ_State_Vector_type
           integer,                    intent(in)     :: iv
           real(kind=cp),              intent(in)     :: xv
           real(kind=cp),              intent(out)    :: ycalc
           Type(LSQ_State_Vector_type),intent(in out) :: Vsa
           Logical, optional,          intent(in)     :: calder
        End Subroutine Model_Functn
       End Interface

       !---- Local variables ----!
       logical                          :: change_par
       logical                          :: singular
       integer                          :: ntrials,ncount
       integer                          :: nfr,j,i,k,ntr
       real(kind=cp)                    :: chisq1
       real(kind=cp), dimension(c%npvar):: beta, der
       real(kind=cp), dimension(c%npvar):: b     !vector of parameter
       real(kind=cp), dimension(c%npvar):: sb    !estimated standard deviations
       Type(LSQ_State_Vector_type)      :: lvs   !local state vector


       ntrials=40
       nt=0
       lvs=vs         !Save current vector state for calling the function

       !---- Optimization with MARQUARDT ALGORITHM ----!
       nfr=d%nobs-c%npvar
       ifail=0
       change_par=.false.
       if (d%nobs > 2000) ntrials= 20

       !---- Evaluate ALFA and BETA matrices       ----!
       !---- (BETA=DP*ALFA)--- DP=BETA*ALFA(-1)    ----!
       beta(:) =0.0
       curv_mat(:,:)=0.0


       do i=1,d%nobs
          call model_functn(i,d%x(i),d%yc(i),lvs,.true.)
          ncount=0
          do j=1,lvs%np
             if (lvs%code(j) == 0) cycle
             ncount=ncount+1
             der(ncount)=lvs%dpv(j)  !Extracting derivatives of free parameters
          end do
          if (c%iw == 1 .and. d%yc(i) /= 0.0) d%sw(i)=1.0/d%yc(i)
          do j=1,c%npvar
             beta(j)=beta(j)+d%sw(i)*(d%y(i)-d%yc(i))*der(j)
             do k=1,j
                curv_mat(j,k)=curv_mat(j,k)+d%sw(i)*der(j)*der(k)
             end do
          end do
       end do

       do i=1,c%npvar
          if (curv_mat(i,i) < 1.0e-30) then
             Err_lsq =.true.
             write(unit=Err_Lsq_Mess,fmt="(a,i5,a,a)")  &
                  " => Singular matrix!!, problem with parameter no.:",i," -> ",trim(namfree(i))
             ifail=2
             return
          end if
       end do

       do j=1,c%npvar
          do k=1,j
             curv_mat(k,j)=curv_mat(j,k)
          end do
       end do

       !---- Evaluate CHI2 at starting point ----!
       chisq1=fchisq(nfr,d%nobs,d%y,d%sw,d%yc)  !chi2 before inverting and applying shifts

       !---- Invert modified curvature matrix to find new parameters
       do ntr=1,ntrials   !Apply the Marquardt procedure without re-calculating derivatives in the call to the function
          nt=nt+1
          ncount=0
          do j=1,lvs%np
             if (lvs%code(j) == 0) cycle
             ncount=ncount+1
             b(ncount)=lvs%pv(j)     !Select and store in b the varying parameters
          end do
          do j=1,c%npvar
             do k=1,c%npvar
                correl(j,k)=curv_mat(j,k)/sqrt(curv_mat(j,j)*curv_mat(k,k))
             end do
             correl(j,j)=1.0+fl
          end do

          call Invert_Matrix(correl(1:c%npvar,1:c%npvar),correl(1:c%npvar,1:c%npvar),singular)

          if (singular) then
             Err_lsq =.true.
             do i=1,c%npvar
                if (correl(i,i) < 1.0e-20) then
                   j=i !iperm(i)
                   exit
                end if
             end do
             write(unit=Err_Lsq_Mess,fmt="(a,i5,a,a)")  &
                  " => Singular matrix!!, problem with parameter no.:",j," -> ",trim(namfree(j))
             ifail=2
             return
          end if

          if (fl < 1.0e-20)  then
             chir=chisq1
             exit
          end if
          !Update parameters
          do j=1,c%npvar
             do k=1,c%npvar
                b(j)=b(j)+beta(k)*correl(j,k)/sqrt(curv_mat(j,j)*curv_mat(k,k))
             end do
          end do
          ncount=0
          do i=1,lvs%np
             if (lvs%code(i) == 0) cycle
             ncount=ncount+1
             lvs%pv(i)=b(ncount)   !Provisional change in local vector state
          end do
          do i=1,d%nobs    !Calculation of the function for all points with the new parameters
             call model_functn(i,d%x(i),d%yc(i),lvs)
          end do

          chir=fchisq(nfr,d%nobs,d%y,d%sw,d%yc)

          if (ntr == ntrials .or. fl > 1.0e+20 ) then    !no improvement
             ifail=1
             change_par=.false.
             exit
          end if

          !---- If Chi2 increase, increase fl and try again ----!
          if (chisq1-chir < 0.0) then  !Chi2 increases
             fl=10.0*fl                !Increase fl
             lvs=vs                    !Restore the value of the local state vector to the input value
             cycle
          end if
          change_par=.true.
          exit
       end do  ! ntrials

       !---- Evaluate parameters and uncertainties ----!
       do j=1,c%npvar
          sb(j)=sqrt(abs(correl(j,j)/curv_mat(j,j)))
       end do
       fl=fl/10.0
       ncount=0
       do i=1,vs%np            !Complete update with standard deviations and evaluate the change
          if (vs%code(i) == 0) then
             vs%spv(i)=0.0
             pn(i)=lvs%pv(i)   !Change in Vs done in the calling routine
             ch(i)=0.0
          else
             ncount=ncount+1
             pn(i)=b(ncount)      !Change in Vs done in the calling routine
             ch(i)=pn(i)-vs%pv(i) !change w.r.t input parameters
             if(change_par) vs%spv(i)=sb(ncount)
          end if
       end do
       return
    End Subroutine Curfit_v2

    !!--++
    !!--++ Subroutine Fdjac2(Fcn, M, N, X, Fvec, Fjac, Iflag, Epsfcn)
    !!--++   Integer,                       Intent(In)      :: m
    !!--++   Integer,                       Intent(In)      :: n
    !!--++   Real (Kind=cp),dimension(n),   Intent(In Out)  :: x
    !!--++   Real (Kind=cp),dimension(m),   Intent(In)      :: fvec
    !!--++   Real (Kind=cp),dimension(:,:), Intent(Out)     :: fjac     ! fjac(ldfjac,n)
    !!--++   Integer,                       Intent(In Out)  :: iflag
    !!--++   Real (Kind=cp),                Intent(In)      :: epsfcn
    !!--++
    !!--++   Interface
    !!--++     Subroutine fcn(m, n, x, fvec, iflag)
    !!--++       Use CFML_GlobalDeps, Only: cp
    !!--++       Integer,                     Intent(In)      :: m, n
    !!--++       Real (Kind=cp),Dimension(:), Intent(In)      :: x
    !!--++       Real (Kind=cp),Dimension(:), Intent(In Out)  :: fvec
    !!--++       Integer,                     Intent(In Out)  :: iflag
    !!--++     End Subroutine fcn
    !!--++   End Interface
    !!--++
    !!--..   Original documentation:
    !!--..
    !!--..  This subroutine computes a forward-difference approximation to the m by n
    !!--..  jacobian matrix associated with a specified problem of m functions in
    !!--..  n variables.
    !!--..  The subroutine statement is
    !!--..
    !!--..    subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
    !!--..
    !!--..  where
    !!--..
    !!--..    fcn   is the name of the user-supplied subroutine which calculates the
    !!--..          functions.  fcn must be declared in an external statement in the
    !!!-..+         user calling program, and should be written as follows.
    !!--..
    !!--..          subroutine fcn(m,n,x,fvec,iflag)
    !!--..             integer m,n,iflag
    !!--..             real (kind=cp) x(n),fvec(m)
    !!--..             ----------
    !!--..             calculate the functions at x and
    !!--..             return this vector in fvec.
    !!--..             ----------
    !!--..             return
    !!--..          end
    !!--..
    !!--..          the value of iflag should not be changed by fcn unless the user
    !!--..          wants to terminate execution of fdjac2. In this case set iflag
    !!--..          to a negative integer.
    !!--..
    !!--..    m     is a positive integer input variable set to the number of functions.
    !!--..
    !!--..    n     is a positive integer input variable set to the number of variables.
    !!--..          n must not exceed m.
    !!--..
    !!--..    x     is an input array of length n.
    !!--..
    !!--..    fvec  is an input array of length m which must contain the functions
    !!--..          evaluated at x.
    !!--..
    !!--..    fjac  is an output m by n array which contains the approximation to the
    !!--..          jacobian matrix evaluated at x.
    !!--..
    !!--..   ldfjac is a positive integer input variable not less than m which specifies
    !!--..          the leading dimension of the array fjac.
    !!--..
    !!--..    iflag is an integer variable which can be used to terminate the execution
    !!--..          of fdjac2.  see description of fcn.
    !!--..
    !!--..   epsfcn is an input variable used in determining a suitable step length
    !!--..          for the forward-difference approximation.  This approximation assumes
    !!--..          that the relative errors in the functions are of the order of epsfcn.
    !!--..          If epsfcn is less than the machine precision, it is assumed that the
    !!--..          relative errors in the functions are of the order of the machine
    !!--..          precision.
    !!--..
    !!--..    wa    is a work array of length m.
    !!--..
    !!--..
    !!--..  Argonne National Laboratory. minpack project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!--++
    !!--++
    !!--++ Update: February - 2009
    !!
    Subroutine Fdjac2(Fcn, M, N, X, Fvec, Fjac, Iflag, Epsfcn)
       !---- Arguments ----!
       Integer,                       Intent(In)      :: m
       Integer,                       Intent(In)      :: n
       Real (Kind=cp),dimension(n),   Intent(In Out)  :: x
       Real (Kind=cp),dimension(m),   Intent(In)      :: fvec
       Real (Kind=cp),dimension(:,:), Intent(Out)     :: fjac     ! fjac(ldfjac,n)
       Integer,                       Intent(In Out)  :: iflag
       Real (Kind=cp),                Intent(In)      :: epsfcn

       Interface
          Subroutine fcn(m, n, x, fvec, iflag)
             Use CFML_GlobalDeps, Only: cp
             Integer,                     Intent(In)      :: m, n
             Real (Kind=cp),Dimension(:), Intent(In)      :: x
             Real (Kind=cp),Dimension(:), Intent(In Out)  :: fvec
             Integer,                     Intent(In Out)  :: iflag
          End Subroutine fcn
       End Interface

       !---- Local Variables ----!
       Integer                   :: j
       Real (Kind=cp)            :: eps, epsmch, h, temp, wa(m)
       Real (Kind=cp), Parameter :: zero = 0.0_CP

       ! epsmch is the machine precision.
       epsmch = Epsilon(zero)

       eps = Sqrt(Max(epsfcn, epsmch))
       Do j = 1, n
          temp = x(j)
          h = eps*Abs(temp)
          If (h == zero) h = eps
          x(j) = temp + h
          Call fcn(m, n, x, wa, iflag)
          If (iflag < 0) Exit
          x(j) = temp
          fjac(1:m,j) = (wa(1:m) - fvec(1:m))/h
       End Do

       Return
    End Subroutine Fdjac2

    !!--..
    !!--..  Subroutine Info_LSQ_LM_V(Chi2,Lun,c,v,vstd,vnam)
    !!--..   real(kind=cp),                 intent(in) :: chi2         !Final Chi2
    !!--..   integer,                       intent(in) :: lun          !Logical unit for output
    !!--..   type(LSQ_conditions_type),     intent(in) :: c            !Conditions of the refinement
    !!--..   real(kind=cp),   dimension(:), intent(in) :: v,vstd       !State vector and standad deviations (parameters of the model)
    !!--..   character(len=*),dimension(:), intent(in) :: vnam         !Names of the refined parameters
    !!--..
    !!--..  Subroutine for output information at the end of refinement of a Levenberg-Marquard fit
    !!--..
    !!---- Update: November 1 - 2013
    !!
    Subroutine Info_LSQ_LM_V(Chi2,Lun,c,v,vstd,vnam)
       !---- Arguments ----!
       real(kind=cp),                 intent(in) :: chi2
       integer,                       intent(in) :: lun
       type(LSQ_conditions_type),     intent(in) :: c
       real(kind=cp),   dimension(:), intent(in) :: v,vstd
       character(len=*),dimension(:), intent(in) :: vnam
       !---- Local variables ----!
       integer       :: i,j,inum

       !---- Correlation matrix ----!
       !Here the correlation matrix is already well calculated

       write(unit=lun,fmt="(/,a,/)")   " => Correlation Matrix: "

       do i=1,c%npvar
          do j=i,c%npvar
             correl(i,j)=correl(i,j)*100.0
             correl(j,i)=correl(i,j)
          end do
       end do

       inum=0
       do i=1,c%npvar-1
          do j=i+1,c%npvar
             if (correl(i,j) > real(c%corrmax) ) then
                write(unit=lun,fmt="(a,i4,a,i2,4a)") "    Correlation:",nint(min(correl(i,j),100.0)),  &
                     " > ",c%corrmax,"% for parameters:   ", adjustr(vnam(i))," & ", vnam(j)
                inum=inum+1
             end if
          end do
       end do
       if (inum == 0) then
          write(unit=lun,fmt="(/,a,i2,a)") " => There is no correlations greater than ",c%corrmax,"% "
       else
          write(unit=lun,fmt="(/,a,i3,a,i2,a)") " => There are ",inum," values of Correlation > ",c%corrmax,"%"
       end if

       write(unit=lun,fmt="(/,/,a,/,a,/)") "      FINAL LIST OF REFINED PARAMETERS AND STANDARD DEVIATIONS",&
                                           "      --------------------------------------------------------"
       write(unit=lun,fmt="(/,a,/)") &
       "    Parameter name     No.(LSQ)         Final-Value   Standard Deviation"
       do i=1,c%npvar
        write(unit=lun,fmt="(a,i6,2f20.5)") "    "//vnam(i),i,v(i),vstd(i)
       end do
       write(unit=lun,fmt="(/,a,g13.5)")  " => Final value of Chi2: ",chi2
       write(unit=lun,fmt="(a)")  " => Well-behaved LSQ-problems give Tikhonov regularization equal to zero (convergence status above)"
       write(unit=lun,fmt="(a)")  &
       " => In ill-behaved LSQ-problems Tikhonov regularization is selected equal to 10^-6*Maximum(Singular Value) in SVDecomposition"
       if(chi2 > 1.0) then
          write(unit=lun,fmt="(a)") " => The LSQ-standard deviations have been mutiplied by SQRT(Chi2)"
       end if

       return
    End Subroutine Info_LSQ_LM_V

    !!--..
    !!--..  Subroutine Info_LSQ_LM_VS(Chi2,Lun,c,vs)
    !!--..   real(kind=cp),              intent(in)     :: chi2       !Final Chi2
    !!--..   integer,                    intent(in)     :: lun        !Logical unit for output
    !!--..   type(LSQ_conditions_type),  intent(in)     :: c          !Conditions of the refinement
    !!--..   type(LSQ_State_Vector_type),intent(in)     :: vs         !State vector (parameters of the model)
    !!--..
    !!--..  Subroutine for output information at the end of refinement of a Levenberg-Marquard fit
    !!--..
    !!--.. Update: November 1 - 2013
    !!
    Subroutine Info_LSQ_LM_VS(Chi2,Lun,c,vs)
       !---- Arguments ----!
       real(kind=cp),              intent(in)     :: chi2
       integer,                    intent(in)     :: lun
       type(LSQ_conditions_type),  intent(in)     :: c
       type(LSQ_State_Vector_type),intent(in)     :: vs

       !---- Local variables ----!
       integer       :: i,j,inum
       real(kind=cp) :: g2


       !---- Correlation matrix ----!
       write(unit=lun,fmt="(/,a,/)")   " => Correlation Matrix: "

       do i=1,c%npvar
          do j=i,c%npvar
             correl(i,j)=correl(i,j)/sqrt(curv_mat(i,i)*curv_mat(j,j))
             correl(j,i)=correl(i,j)
          end do
       end do
       do i=1,c%npvar
          g2=sqrt(correl(i,i))
          do j=i,c%npvar
             correl(i,j)=correl(i,j)/g2/sqrt(correl(j,j))*100.0
             correl(j,i)=correl(i,j)
          end do
       end do

       inum=0
       do i=1,c%npvar-1
          do j=i+1,c%npvar
             if (correl(i,j) > real(c%corrmax) ) then
                write(unit=lun,fmt="(a,i4,a,i2,4a)") "    Correlation:",nint(min(correl(i,j),100.0)),  &
                     " > ",c%corrmax,"% for parameters:   ", adjustr(namfree(i))," & ", namfree(j)
                inum=inum+1
             end if
          end do
       end do
       if (inum == 0) then
          write(unit=lun,fmt="(/,a,i2,a)") " => There is no correlations greater than ",c%corrmax,"% "
       else
          write(unit=lun,fmt="(/,a,i3,a,i2,a)") " => There are ",inum," values of Correlation > ",c%corrmax,"%"
       end if

       write(unit=lun,fmt="(/,/,a,/,a,/)") "      FINAL LIST OF REFINED PARAMETERS AND STANDARD DEVIATIONS",&
                                           "      --------------------------------------------------------"
       write(unit=lun,fmt="(/,a,/)") &
       "    #   Parameter name                       No.(Model)         Final-Value   Standard Deviation"

       inum=0
       do i=1,vs%np
          if (vs%code(i)/=0) then
            inum=inum+1
            write(unit=lun,fmt="(i5,a,i6,2f20.5)") inum,"    "//vs%nampar(i),i,vs%pv(i),vs%spv(i)
          end if
       end do

       write(unit=lun,fmt="(/,a,g13.5)") " => Final value of Chi2: ",chi2


       return
    End Subroutine Info_LSQ_LM_VS

    !!----
    !!----  Subroutine Info_LSQ_Output(Chi2,FL,Nobs,X,Y,Yc,W,Lun,c,vs,out_obscal)
    !!----   real(kind=cp),              intent(in)     :: chi2       !Final Chi2
    !!----   real(kind=cp),              intent(in)     :: FL         !Final Marquardt lambda
    !!----   integer,                    intent(in)     :: nobs       !Number of data points
    !!----   real(kind=cp),dimension(:), intent(in)     :: x          !Array with abcisae of Data points
    !!----   real(kind=cp),dimension(:), intent(in)     :: y          !Array with data point values
    !!----   real(kind=cp),dimension(:), intent(in)     :: yc         !Array with calculated values
    !!----   real(kind=cp),dimension(:), intent(in)     :: w          !Array with weight factors
    !!----   integer,                    intent(in)     :: lun        !Logical unit for output
    !!----   type(LSQ_conditions_type),  intent(in)     :: c          !Conditions of the refinement
    !!----   type(LSQ_State_Vector_type),intent(in)     :: vs         !State vector (parameters of the model)
    !!----   character(len=*), optional, intent(in)     :: out_obscal !If present the vectors X,Y,Yc,Sig(=sqrt(1/w))
    !!----                                                            !Are output in a file called LM_fit.xy
    !!----
    !!----  Subroutine for output information at the end of refinement
    !!----
    !!---- Update: August - 2009
    !!
    Subroutine Info_LSQ_Output(Chi2,FL,Nobs,X,Y,Yc,W,Lun,c,vs,out_obscal)
       !---- Arguments ----!
       real(kind=cp),              intent(in)     :: chi2
       real(kind=cp),              intent(in)     :: FL
       integer,                    intent(in)     :: nobs
       real(kind=cp),dimension(:), intent(in)     :: x
       real(kind=cp),dimension(:), intent(in)     :: y
       real(kind=cp),dimension(:), intent(in)     :: yc
       real(kind=cp),dimension(:), intent(in)     :: w
       integer,                    intent(in)     :: lun
       type(LSQ_conditions_type),  intent(in)     :: c
       type(LSQ_State_Vector_type),intent(in)     :: vs
       character(len=*), optional, intent(in)     :: out_obscal

       !---- Local variables ----!
       integer       :: i,j,inum, lob=22
       real(kind=dp) :: rfact,rwfact,riobs,rex
       real(kind=cp) :: del,g2

       !---- Final calculation and writings R-Factors calculations ----!
       rfact=0.0
       rwfact=0.0
       riobs=0.0
       do i=1,nobs
          riobs=riobs+y(i)
          del=y(i)-yc(i)
          rfact=rfact+abs(del)
          rwfact=rwfact+w(i)*del*del
       end do
       rfact=rfact/riobs*100.0
       rwfact=sqrt(rwfact/riobs)*100.0
       rex=sqrt(real(nobs-c%npvar)/riobs)*100.0
       write(unit=lun,fmt="(/,(3(a,f8.3)))") "  Rfact= ",rfact,"   Rwfact= ",rwfact,"   Rex= ",rex
       write(unit=lun,fmt="(/,a,F16.3)") "  Final value of Marquardt F-Lambda = ",FL

       !---- Correlation matrix ----!
       write(unit=lun,fmt="(/,a,/)")   " => Correlation Matrix: "

       do i=1,c%npvar
          do j=i,c%npvar
             correl(i,j)=correl(i,j)/sqrt(curv_mat(i,i)*curv_mat(j,j))
             correl(j,i)=correl(i,j)
          end do
       end do
       do i=1,c%npvar
          g2=sqrt(correl(i,i))
          do j=i,c%npvar
             correl(i,j)=correl(i,j)/g2/sqrt(correl(j,j))*100.0
             correl(j,i)=correl(i,j)
          end do
       end do

       inum=0
       do i=1,c%npvar-1
          do j=i+1,c%npvar
             if (correl(i,j) > real(c%corrmax) ) then
                write(unit=lun,fmt="(a,i4,a,i2,4a)") "    Correlation:",nint(min(correl(i,j),100.0)),  &
                     " > ",c%corrmax,"% for parameters:   ", adjustr(namfree(i))," & ", namfree(j)
                inum=inum+1
             end if
          end do
       end do
       if (inum == 0) then
          write(unit=lun,fmt="(/,a,i2,a)") " => There is no correlations greater than ",c%corrmax,"% "
       else
          write(unit=lun,fmt="(/,a,i3,a,i2,a)") " => There are ",inum," values of Correlation > ",c%corrmax,"%"
       end if

       write(unit=lun,fmt="(/,/,a,/,a,/)") "      FINAL LIST OF REFINED PARAMETERS AND STANDARD DEVIATIONS",&
                                           "      --------------------------------------------------------"
       write(unit=lun,fmt="(/,a,/)") &
       "    #   Parameter name                       No.(Model)         Final-Value   Standard Deviation"
       inum=0
       do i=1,vs%np
          if (vs%code(i)/=0) then
            inum=inum+1
            write(unit=lun,fmt="(i5,a,i6,2f20.5)") inum,"    "//vs%nampar(i),i,vs%pv(i),vs%spv(i)
          end if
       end do
       write(unit=lun,fmt="(/,a,g13.5)") " => Final value of Chi2: ",chi2

       if (present(out_obscal)) then
          !---- Output of a file with the observed and calculated curves ----!
          open(unit=lob,file="LM_fit.xy",status="replace", action="write")
          write(unit=lob,fmt="(a)") "!        X             Y-obs          Y-calc           Sigma"
          do i=1,nobs
             write(unit=lob,fmt="(4f16.4)") x(i),y(i),yc(i),sqrt(1.0/w(i))
          end do
          close(unit=lob)
       end if

       return
    End Subroutine Info_LSQ_Output


    !!--++
    !!--++   Subroutine LM_Der(Model_Functn, m, c, Vs, chi2, calder, infout,residuals)
    !!--++     Integer,                     Intent(In)      :: m        !Number of observations
    !!--++     type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
    !!--++     type(LSQ_State_Vector_type), Intent(In Out)  :: Vs       !State vector
    !!--++     Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
    !!--++     logical,                     Intent(in)      :: calder   !logical (should be .true.) used only for purposes
    !!--++                                                              !of making unambiguous the generic procedure
    !!--++     character(len=*),            Intent(out)     :: infout   !Information about the refinement (min length 256)
    !!--++     Real (Kind=cp), dimension(:),optional, intent(out) :: residuals
    !!--++
    !!--++     !--- Local Variables ---!
    !!--++     Interface
    !!--++       Subroutine Model_Functn(m, n, x, fvec, fjac, iflag)       !Model Function subroutine
    !!--++         Use CFML_GlobalDeps, Only: cp
    !!--++         Integer,                       Intent(In)    :: m, n    !Number of obserations and free parameters
    !!--++         Real (Kind=cp),Dimension(:),   Intent(In)    :: x       !Array with the values of free parameters: x(1:n)
    !!--++         Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec    !Array of residuals fvec=(y-yc)/sig : fvec(1:m)
    !!--++         Real (Kind=cp),Dimension(:,:), Intent(Out)   :: fjac    !Jacobian Dfvec/Dx(i,j) = [ dfvec(i)/dx(j) ] : fjac(1:m,1:n)
    !!--++         Integer,                       Intent(In Out):: iflag   !If iflag=1 calculate only fvec without changing fjac
    !!--++       End Subroutine Model_Functn                               !If iflag=2 calculate only fjac keeping fvec fixed
    !!--++     End Interface
    !!--++
    !!--++
    !!--++     This interface to lmder has been modified with respect to the original
    !!--++     lmder1 (also available in the overloaded procedure) in order to use the
    !!--++     LSQ types of this module.
    !!--++     The residuals vector and the Jacobian matrix have been eliminated from the
    !!--++     arguments in order to expose to the user only the important information
    !!--++     that are stored in Vs of type LSQ_Vector_State_Type.
    !!--++     The "ipvt" argument has also been suppressed.
    !!--++     If one uses Vs%code_comp=.true.  be careful with the calculations inside
    !!--++     Model_Functn, for making properly the constraints one has to work with
    !!--++     shift of parameters instead of the parameters themselves.
    !!--++
    !!--++
    !!--++  Updated: January - 2014
    !!
    Subroutine LM_Der(Model_Functn, m, c, Vs, chi2, calder, infout,residuals)
       !---- Arguments ----!
       Integer,                               Intent(In)      :: m
       type(LSQ_conditions_type),             Intent(In Out)  :: c
       type(LSQ_State_Vector_type),           Intent(In Out)  :: Vs
       Real (Kind=cp),                        Intent(out)     :: chi2
       logical,                               Intent(in)      :: calder
       character(len=*),                      Intent(out)     :: infout
       real (Kind=cp), dimension(:),optional, Intent(out)     :: residuals

       Interface
         Subroutine Model_Functn(m, n, x, fvec, fjac, iflag)
            Use CFML_GlobalDeps, Only: cp
            Integer,                       Intent(In)    :: m, n
            Real (Kind=cp),Dimension(:),   Intent(In)    :: x
            Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec
            Real (Kind=cp),Dimension(:,:), Intent(Out)   :: fjac
            Integer,                       Intent(In Out):: iflag
         End Subroutine Model_Functn
       End Interface

       !--- Local Variables ---!
       Integer                              :: i,j, maxfev, mode, nfev, njev, nprint, n, &
                                               iflag,info
       Integer,        dimension(c%npvar)   :: ipvt
       Integer, dimension(c%npvar,c%npvar)  :: p,id
       Real (Kind=cp), dimension(c%npvar)   :: x, sx,xo
       Real (Kind=cp), dimension(vs%np)     :: vp
       Real (Kind=cp), dimension(m)         :: fvec  !Residuals
       Real (Kind=cp), dimension(m,c%npvar) :: fjac
       Real (Kind=cp)                       :: ftol, gtol, xtol,iChi2,deni,denj
       Real (Kind=cp), Parameter            :: factor = 100.0_CP, zero = 0.0_CP
       Logical                              :: singular

       info = 0
       n=c%npvar
       c%reached=.false.
       id=0
       do i=1,n
          id(i,i) = 1
       end do
       Chi2=1.0e35
       infout=" "

       !     check the input parameters for errors.
       If ( n <= 0 .OR. m < n .OR. c%tol < zero ) then
         write(unit=infout,fmt="(2a,i5,a,i5,a,f10.7)") "Improper input parameters in LM-optimization (n,m,tol):", &
                        " nb of free (refined) parameters: ", n, &
                        ", number of observations: ", m, &
                        ", tolerance: ", c%tol
          return
       end if

       j=0
       if(vs%code_comp) then
          !save the initial values of all parameters
          vp=vs%pv(1:vs%np)
          do j=1,n
            do i=1,vs%np
               if(vs%code(i) == j) then
                  x(j)=Vs%pv(i)
                  namfree(j)=vs%nampar(i)
                  exit
               end if
            end do
          end do
          xo=x !saving the initial free parameter
       else
          do i=1,vs%np
             if (vs%code(i) == 0) cycle
             j=j+1
             x(j)=Vs%pv(i)
             namfree(j)=vs%nampar(i)
          end do
       end if

       ! Initial calculation of Chi2
       iflag=1
       if(calder) iflag=1
       Call Model_Functn(m, n, x, fvec, fjac, iflag)
       ichi2=enorm(m,fvec)
       if (ichi2 < 1.0e15) then
          ichi2=ichi2*ichi2/max(1.0_cp,real(m-n,kind=cp))
       end if
       !If ( c%icyc < 5*n) then
       !   maxfev = 100*(n + 1)
       !   c%icyc= maxfev
       !Else
          maxfev = c%icyc
       !End if
       ftol = c%tol
       xtol = c%tol
       gtol = zero
       mode = 1
       nprint = c%nprint

       ! Call to the core procedure of the MINPACK implementation of the Levenberg-Marquardt algorithm
       Call lmder(Model_Functn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev,  &
                  mode, factor, nprint, info, nfev, njev, ipvt)

       !Calculate final Chi2 or norm
       chi2=enorm(m,fvec)
       if (chi2 < 1.0e15) then
          chi2=chi2*chi2/max(1.0_cp,real(m-n,kind=cp))
       end if
       sx=zero

       !Extract the curvature matrix side of equation below from fjac
       !        t     t           t
       !       p *(jac *jac)*p = r *r
       !           t            t     t
       ! cvm = (jac *jac) = p*(r *r)*p     -> Permutation matrices are orthogonal
       ! See the documentation of lmder1 and lmdif1 for why we use only the (n,n) submatrix of fjac
       !curv_mat(1:n,1:n)=matmul( transpose( fjac(1:n,1:n) ) , fjac(1:n,1:n) )
       curv_mat(1:n,1:n)=matmul( transpose( Upper_Triangular(fjac,n)) , Upper_Triangular(fjac,n) )
       do j=1,n
          p(1:n,j) = id(1:n,ipvt(j))
       end do
       curv_mat(1:n,1:n) = matmul(  p, matmul( curv_mat(1:n,1:n),transpose(p) )  )

       !Make an additional direct final call to the model function in order to calculate the Jacobian J,
       !curvature matrix A=Jt.J, invert it and take the diagonal elements for calculating the standard
       !deviations at the final point.
       !iflag=2
       !Call Model_Functn(m, n, x, fvec, fjac, iflag)
       !curv_mat(1:n,1:n) = matmul (transpose(fjac),fjac)
       do j=1,n
          denj=curv_mat(j,j)
          if( denj <= zero) denj=1.0
          do i=1,n
             deni=curv_mat(i,i)
             if( deni <= zero) deni=1.0
             correl(j,i)=curv_mat(j,i)/sqrt(denj*deni)
          end do
          correl(j,j)=1.00001
       end do
       Call Invert_Matrix(correl(1:n,1:n),correl(1:n,1:n),singular)
       If (.not. singular) then
          Do i=1,n
             deni = curv_mat(i,i)
             if( deni <= zero) deni=1.0
             sx(i) = sqrt(abs(correl(i,i)/deni))         !sqrt(abs(correl(i,i)*Chi2))
          End Do
       Else
          Do i=1,n
             if (correl(i,i) <= zero) then
                info=0
                write(unit=infout,fmt="(a,i3)") "Final Singular Matrix!, problem with parameter #",i
                exit
             end if
          End Do
          Err_lsq =.true.
          Err_Lsq_Mess=" Singular Matrix at the end of LM_Der for calculating std deviations ..."
       End if

       !Update the State vector Vs
       n=0
       if(vs%code_comp) then
          Do i=1,Vs%np
             if (Vs%code(i) == 0) then
                 Vs%spv(i)=0.0
             Else
                n=vs%code(i)
                Vs%pv(i) = vp(i)+(x(n)-xo(n))*Vs%mul(i)
                Vs%spv(i)=sx(n)*Vs%mul(i)
             End if
          End do
       else
          Do i=1,Vs%np
             if (Vs%code(i) == 0) then
                 Vs%spv(i)=0.0
             Else
                n=n+1
                Vs%pv(i) = x(n)
                Vs%spv(i)=sx(n)
             End if
          End do
       end if

       Select Case (info)
          Case(0)
             c%reached=.false.
          Case(1)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Initial Chi2:", ichi2, &
               "     Convergence reached: The relative error in the sum of squares is at most ",c%tol

          Case(2)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Initial Chi2:", ichi2,&
                "     Convergence reached: The relative error between x and the solution is at most ",c%tol

          Case(3)
            c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Initial Chi2:", ichi2,&
             "     Convergence reached: The relative error "// &
             "in the sum of squares and the difference between x and the solution are both at most ",c%tol

          Case(4,8)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a)") "Initial Chi2:", ichi2, &
              "     Convergence reached: Residuals vector is orthogonal to the columns of the Jacobian to machine precision"

          Case(5)
             c%reached=.false.
             write(unit=infout,fmt="(a,i6)") &
             "Convergence NOT reached: Number of calls to model function with iflag=1 has reached ",maxfev

          Case(6)
             c%reached=.false.
             write(unit=infout,fmt="(a,e12.5,a)") "Convergence NOT reached: Provided tolerance ",c%tol,&
                                            " is too small! No further reduction in the sum of squares is possible"
          Case(7)
             c%reached=.false.
             write(unit=infout,fmt="(a,e12.5,a)") "Convergence NOT reached: Provided tolerance ",c%tol,&
                                            " is too small! No further improvement in the approximate solution x is possible"
       End Select
       if (present(residuals)) residuals(1:m)=fvec(1:m)

       Return
    End Subroutine LM_Der

    !!--++
    !!--++   Subroutine LM_DerV(Model_Functn, m, c, V, Vstd, chi2, calder, infout,residuals)
    !!--++     Integer,                     Intent(In)      :: m        !Number of observations
    !!--++     type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
    !!--++     Real (Kind=cp),dimension(:), Intent(In Out)  :: V        !State vector
    !!--++     Real (Kind=cp),dimension(:), Intent(In Out)  :: Vstd     !Standard deviations vector
    !!--++     Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
    !!--++     logical,                     Intent(in)      :: calder   !logical (should be .true.) used only for purposes
    !!--++                                                              !of making unambiguous the generic procedure
    !!--++     character(len=*),            Intent(out)     :: infout   !Information about the refinement (min length 256)
    !!--++     Real (Kind=cp), dimension(:),optional, intent(out) :: residuals
    !!--++
    !!--++     !--- Local Variables ---!
    !!--++     Interface
    !!--++       Subroutine Model_Functn(m, n, x, fvec, fjac, iflag)       !Model Function subroutine
    !!--++         Use CFML_GlobalDeps, Only: cp
    !!--++         Integer,                       Intent(In)    :: m, n    !Number of obserations and free parameters
    !!--++         Real (Kind=cp),Dimension(:),   Intent(In)    :: x       !Array with the values of free parameters: x(1:n)
    !!--++         Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec    !Array of residuals fvec=(y-yc)/sig : fvec(1:m)
    !!--++         Real (Kind=cp),Dimension(:,:), Intent(Out)   :: fjac    !Jacobian Dfvec/Dx(i,j) = [ dfvec(i)/dx(j) ] : fjac(1:m,1:n)
    !!--++         Integer,                       Intent(In Out):: iflag   !If iflag=1 calculate only fvec without changing fjac
    !!--++       End Subroutine Model_Functn                               !If iflag=2 calculate only fjac keeping fvec fixed
    !!--++     End Interface
    !!--++
    !!--++
    !!--++     This interface to lmder has been modified with respect to the original
    !!--++     lmder1 (also available in the overloaded procedure) in order to use the
    !!--++     LSQ types of this module. This is similar to LM_DER
    !!--++     In this case we use only the LSQ_conditions_type 'c' derived type.
    !!--++     The state vector is similar to V_vec in CFML_Refcodes.
    !!--++
    !!--++
    !!--++ Update: November 1 - 2013
    !!
    Subroutine LM_DerV(Model_Functn, m, c, V, Vstd, chi2, calder, infout,residuals)
       !---- Arguments ----!
       Integer,                               Intent(In)      :: m
       type(LSQ_conditions_type),             Intent(In Out)  :: c
       Real (Kind=cp), dimension(:),          Intent(In Out)  :: V,Vstd
       Real (Kind=cp),                        Intent(out)     :: chi2
       logical,                               Intent(in)      :: calder
       character(len=*),                      Intent(out)     :: infout
       real (Kind=cp), dimension(:),optional, Intent(out)     :: residuals

       Interface
         Subroutine Model_Functn(m, n, x, fvec, fjac, iflag)
            Use CFML_GlobalDeps, Only: cp
            Integer,                       Intent(In)    :: m, n
            Real (Kind=cp),Dimension(:),   Intent(In)    :: x
            Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec
            Real (Kind=cp),Dimension(:,:), Intent(Out)   :: fjac
            Integer,                       Intent(In Out):: iflag
         End Subroutine Model_Functn
       End Interface

       !--- Local Variables ---!
       Integer                                    :: i,j, maxfev, mode, nfev, njev, nprint, n, &
                                                     iflag,info
       Integer,        dimension(c%npvar)         :: ipvt
       Integer,        dimension(c%npvar,c%npvar) :: p,id
       Real (Kind=cp), dimension(c%npvar,c%npvar) :: U,VS,Sigm
       Real (Kind=cp), dimension(c%npvar)         :: Sig
       Real (Kind=cp), dimension(m)               :: fvec  !Residuals
       Real (Kind=cp), dimension(m,c%npvar)       :: fjac
       Real (Kind=cp)                             :: ftol, gtol, xtol,iChi2,deni,denj, &
                                                     Tikhonov
       Real (Kind=cp), Parameter                  :: factor = 100.0_CP, zero = 0.0_CP
       Logical                                    :: singular

       info = 0
       n=c%npvar
       c%reached=.false.
       id=0
       do i=1,n
          id(i,i) = 1
       end do
       Chi2=1.0e35
       infout=" "

       !     check the input parameters for errors.
       If ( n <= 0 .OR. m < n .OR. c%tol < zero ) then
         write(unit=infout,fmt="(2a,i5,a,i5,a,f10.7)") "Improper input parameters in LM-optimization (n,m,tol):", &
                        " nb of free (refined) parameters: ", n, &
                        ", number of observations: ", m, &
                        ", tolerance: ", c%tol
          return
       end if

       ! Initial calculation of Chi2
       iflag=1
       if(calder) iflag=1  !Just for for purposes of making unambiguous the generic procedure
       Call Model_Functn(m, n, v, fvec, fjac, iflag)
       ichi2=enorm(m,fvec)
       if (ichi2 < 1.0e15) then
          ichi2=ichi2*ichi2/max(1.0_cp,real(m-n,kind=cp))
       end if
       If ( c%icyc < 5*n) then
          maxfev = 100*(n + 1)
          c%icyc= maxfev
       Else
          maxfev = c%icyc
       End if
       ftol = c%tol
       xtol = c%tol
       gtol = zero
       mode = 1
       nprint = c%nprint

       ! Call to the core procedure of the MINPACK implementation of the Levenberg-Marquardt algorithm
       Call lmder(Model_Functn, m, n, v, fvec, fjac, ftol, xtol, gtol, maxfev,  &
                  mode, factor, nprint, info, nfev, njev, ipvt)

       !Calculate final Chi2 or norm
       chi2=enorm(m,fvec)
       if (chi2 < 1.0e15) then
          chi2=chi2*chi2/max(1.0_cp,real(m-n,kind=cp))
       end if
       Vstd=zero

       !Extract the curvature matrix side of equation below from fjac
       !        t     t           t
       !       p *(jac *jac)*p = r *r
       !           t            t     t
       ! cvm = (jac *jac) = p*(r *r)*p     -> Permutation matrices are orthogonal
       ! See the documentation of lmder1 and lmdif1 for why we use only the (n,n) submatrix of fjac
       !curv_mat(1:n,1:n)=matmul( transpose( fjac(1:n,1:n) ) , fjac(1:n,1:n) )  !this is Rt.R

        curv_mat(1:n,1:n)=matmul( transpose( Upper_Triangular(fjac,n)) , Upper_Triangular(fjac,n) )
        do j=1,n
           p(1:n,j) = id(1:n,ipvt(j))
        end do
        curv_mat(1:n,1:n) = matmul(  p, matmul( curv_mat(1:n,1:n),transpose(p) )  )
        Call Invert_Matrix(curv_mat(1:n,1:n),correl(1:n,1:n),singular)
        !If the final curvature matrix is singular perform a Tikhonov
        !regularization (This increases the error bars!)
        Tikhonov=0.0
        if(singular) then
          Err_lsq =.true.
          Err_Lsq_Mess="Regularization (SVD,Tikhonov) of the final Curvature Matrix unsuccessfull (no standard deviations) ..."
          j=0
          do
            !first do a SVD decomposition
            U=curv_mat(1:n,1:n)
            call SVDcmp(U,Sig,VS)
            Tikhonov=maxval(Sig)*1.0e-6
            Sig=Sig+Tikhonov
            Sigm=0.0
            do i=1,n
               Sigm(i,i)=Sig(i)
            end do
            curv_mat(1:n,1:n)=matmul(U,matmul(Sigm,transpose(VS)))
            Call Invert_Matrix(curv_mat(1:n,1:n),correl(1:n,1:n),singular)
            if(.not. singular) exit
            j=j+1
            if(j > 3) exit
          end do
          if(j <= 3) then
            Err_Lsq_Mess="Regularization (SVD,Tikhonov) of the final Curvature Matrix OK! ... bigger Standard Deviations"
          end if
        end if
       !
       !Alternatively, make an additional direct final call to the model function in order to calculate the Jacobian J,
       !curvature matrix A=Jt.J, invert it and take the diagonal elements for calculating the standard
       !deviations at the final point.
       !iflag=2
       !Call Model_Functn(m, n, v, fvec, fjac, iflag)
       !curv_mat(1:n,1:n) = matmul (transpose(fjac),fjac) !Least squares matrix

       If (.not. singular) then
          !Here correl is the inverse of the curvature matrix (Variance-co-variance matrix)
          deni=max(chi2,1.0)
          Do i=1,n
             Vstd(i) = sqrt(deni*abs(correl(i,i)))
          End Do
          !Now correl is the true correlation matrix
          forall (i=1:n,j=1:n)
              correl(i,j)=correl(i,j)/sqrt(correl(i,i)*correl(j,j))
          end forall
       Else
          !Here correl is the "pseudo"-inverse of the curvature-matrix
          Do i=1,n
             if( deni <= zero) then
               Vstd(i) = 99999.9
             else
               Vstd(i) = sqrt(chi2*abs(correl(i,i)/deni))         !sqrt(abs(correl(i,i)*Chi2))
             end if
          End Do
          !Pseudo-correlation matrix
          do j=1,n
             denj=correl(j,j)
             if( denj <= zero) denj=1.0
             do i=1,n
                deni=correl(i,i)
                if( deni <= zero) deni=1.0
                correl(j,i)=correl(j,i)/sqrt(denj*deni)
             end do
          end do
       End if


       Select Case (info)
          Case(0)
             c%reached=.false.
          Case(1)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Tikhonov regularization:", Tikhonov, &
               "     Convergence reached: The relative error in the sum of squares is at most ",c%tol

          Case(2)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Tikhonov regularization:", Tikhonov,&
                "     Convergence reached: The relative error between x and the solution is at most ",c%tol

          Case(3)
            c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Tikhonov regularization:", Tikhonov,&
             "     Convergence reached: The relative error "// &
             "in the sum of squares and the difference between x and the solution are both at most ",c%tol

          Case(4,8)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a)") "Tikhonov regularization:", Tikhonov, &
              "     Convergence reached: Residuals vector is orthogonal to the columns of the Jacobian to machine precision"

          Case(5)
             c%reached=.false.
             write(unit=infout,fmt="(a,i6)") &
             "Convergence NOT reached: Number of calls to model function with iflag=1 has reached ",maxfev

          Case(6)
             c%reached=.false.
             write(unit=infout,fmt="(a,e12.5,a)") "Convergence NOT reached: Provided tolerance ",c%tol,&
                                            " is too small! No further reduction in the sum of squares is possible"
          Case(7)
             c%reached=.false.
             write(unit=infout,fmt="(a,e12.5,a)") "Convergence NOT reached: Provided tolerance ",c%tol,&
                                            " is too small! No further improvement in the approximate solution x seems to be possible"
       End Select
       if (present(residuals)) residuals(1:m)=fvec(1:m)

       Return
    End Subroutine LM_DerV

    !!--++
    !!--++   Subroutine LM_Dif(Model_Functn, m, c, Vs, chi2, infout,residuals,idebug)
    !!--++     Integer,                     Intent(In)      :: m        !Number of observations
    !!--++     type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
    !!--++     type(LSQ_State_Vector_type), Intent(In Out)  :: Vs       !State vector
    !!--++     Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
    !!--++     character(len=*),            Intent(out)     :: infout   !Information about the refinement (min length 256)
    !!--++     Real (Kind=cp), dimension(:),optional, intent(out) :: residuals
    !!--++     integer,                     optional, intent(in)  :: idebug !logical unit for writing results
    !!--++     !--- Local Variables ---!
    !!--++     Interface
    !!--++       Subroutine Model_Functn(m, n, x, fvec, iflag)             !Model Function subroutine
    !!--++         Use CFML_GlobalDeps, Only: cp
    !!--++         Integer,                       Intent(In)    :: m, n    !Number of observations and free parameters
    !!--++         Real (Kind=cp),Dimension(:),   Intent(In)    :: x       !Array with the values of free parameters: x(1:n)
    !!--++         Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec    !Array of residuals fvec=(y-yc)/sig : fvec(1:m)
    !!--++         Integer,                       Intent(In Out):: iflag   !If iflag=1 calculate only fvec without changing fjac
    !!--++       End Subroutine Model_Functn                               !If iflag=2 calculate only fjac keeping fvec fixed
    !!--++     End Interface
    !!--++
    !!--++
    !!--++
    !!--++
    !!--++  This interface to lmdif has been modified with respect to the original lmdif1 (also available in
    !!--++  the overloaded procedure) in order to use the LSQ types of this module.
    !!--++  The residuals vector has been eliminated from the arguments in order to expose to the
    !!--++  user only the important information that are stored in Vs of type LSQ_Vector_State_Type.
    !!--++  The "iwa" argument has also been suppressed.
    !!--++  If one uses Vs%code_comp=.true.  be careful with the calculations inside
    !!--++  Model_Functn, for making properly the constraints one has to work with
    !!--++  shift of parameters instead of the parameters themselves.
    !!--++
    !!--++ Updated: January - 2014
    !!
    Subroutine LM_Dif(Model_Functn, m, c, Vs, chi2, infout,residuals,idebug)
       !---- Arguments ----!
       Integer,                     Intent(In)      :: m        !Number of observations
       type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
       type(LSQ_State_Vector_type), Intent(In Out)  :: Vs       !State vector
       Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
       character(len=*),            Intent(out)     :: infout   !Information about the refinement  (min length 256)
       Real (Kind=cp), dimension(:),optional, intent(out) :: residuals
       integer,                     optional, intent(in)  :: idebug !logical unit for writing results

       Interface
         Subroutine Model_Functn(m, n, x, fvec, iflag)             !Model Function subroutine
            Use CFML_GlobalDeps, Only: cp
            Integer,                       Intent(In)    :: m, n    !Number of observations and free parameters
            Real (Kind=cp),Dimension(:),   Intent(In)    :: x       !Array with the values of free parameters: x(1:n)
            Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec    !Array of residuals fvec=(y-yc)/sig : fvec(1:m)
            Integer,                       Intent(In Out):: iflag   !If iflag=1 calculate only fvec without changing fjac
         End Subroutine Model_Functn                               !If iflag=2 calculate only fjac keeping fvec fixed
       End Interface

       !--- Local Variables ---!
       Integer                              :: i,j, maxfev, mode, nfev, nprint, n, &
                                               iflag,info
       Integer,        dimension(c%npvar)   :: ipvt
       Integer, dimension(c%npvar,c%npvar)  :: p,id
       Real (Kind=cp), dimension(c%npvar)   :: x, sx,xo
       Real (Kind=cp), dimension(vs%np)     :: vp
       Real (Kind=cp), dimension(m)         :: fvec  !Residuals
       Real (Kind=cp), dimension(m,c%npvar) :: fjac
       Real (Kind=cp)                       :: epsfcn,ftol, gtol, xtol,iChi2,deni,denj
       Real (Kind=cp), Parameter            :: factor = 100.0_CP, zero = 0.0_CP
       Logical                              :: singular
       character(len=20)                    :: formt

       info = 0
       n=c%npvar
       c%reached=.false.
       id=0
       do i=1,n
          id(i,i) = 1
       end do
       Chi2=1.0e35
       infout=" "
       !     check the input parameters for errors.
       If ( n <= 0 .OR. m < n .OR. c%tol < zero ) then
         write(unit=infout,fmt="(2a,i5,a,i5,a,f10.7)") "Improper input parameters in LM-optimization (n,m,tol):", &
                        " nb of free (refined) parameters: ", n, &
                        ", number of observations: ", m, &
                        ", tolerance: ", c%tol
          Return
       End if


       j=0
       if(vs%code_comp) then
          !save the initial values of all parameters
          vp=vs%pv(1:vs%np)
          do j=1,n
            do i=1,vs%np
               if(vs%code(i) == j) then
                  x(j)=Vs%pv(i)
                  namfree(j)=vs%nampar(i)
                  exit
               end if
            end do
          end do
          xo=x !saving the initial free parameter
       else
          do i=1,vs%np
             if (vs%code(i) == 0) cycle
             j=j+1
             x(j)=Vs%pv(i)
             namfree(j)=vs%nampar(i)
          end do
       end if

       !Initial calculation of Chi2
       iflag=1
       Call Model_Functn(m, n, x, fvec, iflag)
       ichi2=enorm(m,fvec)
       if (ichi2 < 1.0e15) then
          ichi2=ichi2*ichi2/max(1.0_cp,real(m-n,kind=cp))
       end if
       If ( c%icyc < 5*n) then
          maxfev = 200*(n + 1)
          c%icyc= maxfev
       Else
          maxfev = c%icyc
       End if
       ftol = c%tol
       xtol = c%tol
       gtol = zero
       epsfcn = zero
       mode = 1
       nprint = c%nprint

       ! Call to the core procedure of the MINPACK implementation of the Levenberg-Marquardt algorithm
       Call lmdif(Model_Functn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,   &
                  mode, factor, nprint, info, nfev, fjac, ipvt)

       !Calculate final Chi2 or norm
       chi2=enorm(m,fvec)
       if (chi2 < 1.0e15) then
          chi2=chi2*chi2/max(1.0_cp,real(m-n,kind=cp))
       end if
       sx=zero
       !Extract the curvature matrix side of equation below from fjac
       !        t     t           t
       !       p *(jac *jac)*p = r *r
       !           t            t     t
       ! cvm = (jac *jac) = p*(r *r)*p     -> Permutation matrices are orthogonal
       ! See the documentation of lmder1 and lmdif1 for why we use only the (n,n) submatrix of fjac
       !curv_mat(1:n,1:n)=matmul( transpose( fjac(1:n,1:n) ) , fjac(1:n,1:n) )
       curv_mat(1:n,1:n)=matmul( transpose( Upper_Triangular(fjac,n)) , Upper_Triangular(fjac,n) )
       do j=1,n
          p(1:n,j) = id(1:n,ipvt(j))
       end do
       curv_mat(1:n,1:n) = matmul(  p, matmul( curv_mat(1:n,1:n),transpose(p) )  )

       do j=1,n
          denj=curv_mat(j,j)
          if( denj <= zero) denj=1.0
          do i=1,n
             deni=curv_mat(i,i)
             if( deni <= zero) deni=1.0
             correl(j,i)=curv_mat(j,i)/sqrt(denj*deni)
          end do
          correl(j,j)=1.0001
       end do

       if(present(idebug)) then
         formt="(100g10.3)"
         !write(unit=formt(2:5),fmt="(i4)") n/12
         write(unit=idebug,fmt="(a,i5,a,/)") "  Output pseudo-Jacobian for ",n," free parameters"
         do j=1,n
           write(unit=idebug,fmt=formt) fjac(j,1:n)
         end do
         write(unit=idebug,fmt="(/,a,i5,a,/)") "  Diagonal of the calculated Curvature matrix for ",n," free parameters"
         write(unit=idebug,fmt=formt) (curv_mat(j,j),j=1,n)
         write(unit=idebug,fmt="(a,g14.4)") "  Maxvalue of Correl before inversion: ",maxval(correl)
         write(unit=idebug,fmt="(a,g14.4)") "  Minvalue of Correl before inversion: ",minval(correl)
       end if

       Call Invert_Matrix(correl(1:n,1:n),correl(1:n,1:n),singular)

       if(present(idebug)) then
         write(unit=idebug,fmt="(/,a,i5,a,/)") "  Diagonal of the Inverse Correlation matrix for ",n," free parameters"
         write(unit=idebug,fmt=formt) (correl(j,j),j=1,n)
         write(unit=idebug,fmt="(/,a,i5,a,/)") "  Variances for ",n," free parameters"
         write(unit=idebug,fmt=formt) (correl(j,j)/curv_mat(j,j),j=1,n)

       end if

       If (.not. singular) then
          Do i=1,n
             deni = curv_mat(i,i)
             if( deni <= zero) deni=1.0
             sx(i) = sqrt(abs(correl(i,i)/deni))         !sqrt(abs(correl(i,i)*Chi2))
          End Do
       Else
          Do i=1,n
             deni = curv_mat(i,i)
             if( deni <= zero) then
               sx(i) = 99999.9
             else
               sx(i) = sqrt(abs(correl(i,i)/deni))         !sqrt(abs(correl(i,i)*Chi2))
             end if
          End Do
          Err_lsq =.true.
          Err_Lsq_Mess=" => Singular Matrix at the end of LM_Der for calculating std deviations ... provided sigmas are dubious!"
       End if

       !Update the State vector Vs
       n=0
       if(vs%code_comp) then
          Do i=1,Vs%np
             if (Vs%code(i) == 0) then
                 Vs%spv(i)=0.0
             Else
                n=vs%code(i)
                Vs%pv(i) = vp(i)+(x(n)-xo(n))*Vs%mul(i)
                Vs%spv(i)=sx(n)*Vs%mul(i)
             End if
          End do
       else
          Do i=1,Vs%np
             if (Vs%code(i) == 0) then
                 Vs%spv(i)=0.0
             Else
                n=n+1
                Vs%pv(i) = x(n)
                Vs%spv(i)=sx(n)
             End if
          End do
       end if

       Select Case (info)
          Case(0)
             c%reached=.false.
          Case(1)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Initial Chi2:", ichi2, &
                   "     Convergence reached: The relative error in the sum of squares is at most ",c%tol

          Case(2)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Initial Chi2:", ichi2,&
                    "     Convergence reached: The relative error between x and the solution is at most ",c%tol

          Case(3)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a,e12.5)") "Initial Chi2:", ichi2, &
             "    Convergence reached: The relative error"// &
             " in the sum of squares and the difference between x and the solution are both at most ",c%tol

          Case(4,8)
             c%reached=.true.
             write(unit=infout,fmt="(a,g12.5,a)") "Initial Chi2:", ichi2, &
              "    Convergence reached: Residuals vector is orthogonal to the columns of the Jacobian to machine precision"

          Case(5)
             c%reached=.false.
             write(unit=infout,fmt="(a,i6)") &
             "Convergence NOT reached: Number of calls to model function with iflag=1 has reached ",maxfev

          Case(6)
             c%reached=.false.
             write(unit=infout,fmt="(a,e12.5,a)") "Convergence NOT reached: Provided tolerance ",c%tol,&
                                            " is too small! No further reduction in the sum of squares is possible"
          Case(7)
             c%reached=.false.
             write(unit=infout,fmt="(a,e12.5,a)") "Convergence NOT reached: Provided tolerance ",c%tol,&
                                            " is too small! No further improvement in the approximate solution x is possible"
       End Select
       if (present(residuals)) residuals(1:m)=fvec(1:m)

       Return
    End Subroutine LM_Dif

    !!--++
    !!--++  Subroutine lmder(fcn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev, &
    !!--++                   mode, factor, nprint, info, nfev, njev, ipvt)
    !!--++    Integer,                        Intent(In)      :: m
    !!--++    Integer,                        Intent(In)      :: n
    !!--++    Real (Kind=cp), Dimension(:),   Intent(In Out)  :: x
    !!--++    Real (Kind=cp), Dimension(m),   Intent(Out)     :: fvec
    !!--++    Real (Kind=cp), Dimension(:,:), Intent(Out)     :: fjac    ! fjac(ldfjac,n)
    !!--++    Real (Kind=cp),                 Intent(In)      :: ftol
    !!--++    Real (Kind=cp),                 Intent(In)      :: xtol
    !!--++    Real (Kind=cp),                 Intent(In Out)  :: gtol
    !!--++    Integer,                        Intent(In Out)  :: maxfev
    !!--++    Integer,                        Intent(In)      :: mode
    !!--++    Real (Kind=cp),                 Intent(In)      :: factor
    !!--++    Integer,                        Intent(In)      :: nprint
    !!--++    Integer,                        Intent(Out)     :: info
    !!--++    Integer,                        Intent(Out)     :: nfev
    !!--++    Integer,                        Intent(Out)     :: njev
    !!--++    Integer,        Dimension(:),   Intent(Out)     :: ipvt
    !!--++
    !!--++    Interface
    !!--++      Subroutine fcn(m, n, x, fvec, fjac, iflag)
    !!--++        Use CFML_GlobalDeps, Only: cp
    !!--++        Integer,                       Intent(In)      :: m, n
    !!--++        Real (Kind=cp),Dimension(:),   Intent(In)      :: x
    !!--++        Real (Kind=cp),Dimension(:),   Intent(In Out)  :: fvec
    !!--++        Real (Kind=cp),Dimension(:,:), Intent(Out)     :: fjac
    !!--++        Integer,                       Intent(In Out)  :: iflag
    !!--++      End Subroutine fcn
    !!--++    End Interface
    !!--++
    !!--++
    !!--..  Original Documentation
    !!--..
    !!--..  The purpose of lmder is to minimize the sum of the squares of m nonlinear
    !!--..  functions in n variables by a modification of the levenberg-marquardt
    !!--..  algorithm. the user must provide a subroutine which calculates the functions
    !!--..  and the jacobian.
    !!--..
    !!--..  The subroutine statement is
    !!--..
    !!--..      subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, &
    !!--..                       maxfev,diag,mode,factor,nprint,info,nfev,  &
    !!--..                       njev,ipvt,qtf,wa1,wa2,wa3,wa4)
    !!--..
    !!--..  where
    !!--..
    !!--..    fcn      is the name of the user-supplied subroutine which calculates the
    !!--..             functions and the jacobian. fcn must be declared in an external
    !!--..             statement in the user calling program, and should be written as
    !!--..             follows.
    !!--..
    !!--..             subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    !!--..                integer m,n,ldfjac,iflag
    !!--..                real (kind=cp) x(:),fvec(m),fjac(ldfjac,n)
    !!--..                ----------
    !!--..                if iflag = 1 calculate the functions at x and
    !!--..                return this vector in fvec. do not alter fjac.
    !!--..                if iflag = 2 calculate the jacobian at x and
    !!--..                return this matrix in fjac.  Do not alter fvec.
    !!--..                ----------
    !!--..                return
    !!--..             end subroutine
    !!--..
    !!--..             The value of iflag should not be changed by fcn unless the user
    !!--..             wants to terminate execution of lmder. In this case set iflag to
    !!--..             a negative integer.
    !!--..
    !!--..    m        is a positive integer input variable set to the number of functions.
    !!--..
    !!--..    n        is a positive integer input variable set to the number of variables.
    !!--..             n must not exceed m.
    !!--..
    !!--..    x        is an array of length n. on input x must contain an initial estimate
    !!--..             of the solution vector. on output x contains the final estimate of
    !!--..             the solution vector.
    !!--..
    !!--..    fvec     is an output array of length m which contains the functions evaluated
    !!--..             at the output x.
    !!--..
    !!--..    fjac     is an output m by n array. the upper n by n submatrix of fjac
    !!--..             contains an upper triangular matrix r with diagonal elements of
    !!--..             nonincreasing magnitude such that
    !!--..
    !!--..                  t     t           t
    !!--..                 p *(jac *jac)*p = r *r
    !!--..
    !!--..             where p is a permutation matrix and jac is the final calculated
    !!--..             jacobian.  Column j of p is column ipvt(j) (see below) of the
    !!--..             identity matrix.  The lower trapezoidal part of fjac contains
    !!--..             information generated during the computation of r.
    !!--..
    !!--..    ldfjac   is a positive integer input variable not less than m which specifies
    !!--..             the leading dimension of the array fjac.
    !!--..
    !!--..    ftol     is a nonnegative input variable.  Termination occurs when both
    !!--..             the actual and predicted relative reductions in the sum of squares
    !!--..             are at most ftol.   Therefore, ftol measures the relative error
    !!--..             desired in the sum of squares.
    !!--..
    !!--..    xtol     is a nonnegative input variable. termination occurs when the relative
    !!--..             error between two consecutive iterates is at most xtol. therefore,
    !!--..             xtol measures the relative error desired in the approximate solution.
    !!--..
    !!--..    gtol     is a nonnegative input variable.  Termination occurs when the cosine
    !!--..             of the angle between fvec and any column of the jacobian is at most
    !!--..             gtol in absolute value.  Therefore, gtol measures the orthogonality
    !!--..             desired between the function vector and the columns of the jacobian.
    !!--..
    !!--..    maxfev   is a positive integer input variable.  Termination occurs when the
    !!--..             number of calls to fcn with iflag = 1 has reached maxfev.
    !!--..
    !!--..    diag     is an array of length n.  If mode = 1 (see below), diag is internally
    !!--..             set.  If mode = 2, diag must contain positive entries that serve as
    !!--..             multiplicative scale factors for the variables.
    !!--..
    !!--..    mode     is an integer input variable.  if mode = 1, the variables will be
    !!--..             scaled internally.  if mode = 2, the scaling is specified by the
    !!--..             input diag.  other values of mode are equivalent to mode = 1.
    !!--..
    !!--..    factor   is a positive input variable used in determining the initial step
    !!--..             bound. this bound is set to the product of factor and the euclidean
    !!--..             norm of diag*x if nonzero, or else to factor itself. in most cases
    !!--..             factor should lie in the interval (.1,100.).100. is a generally
    !!--..             recommended value.
    !!--..
    !!--..    nprint   is an integer input variable that enables controlled printing of
    !!--..             iterates if it is positive.  In this case, fcn is called with
    !!--..             iflag = 0 at the beginning of the first iteration and every nprint
    !!--..             iterations thereafter and immediately prior to return, with x, fvec,
    !!--..             and fjac available for printing.  fvec and fjac should not be
    !!--..             altered.  If nprint is not positive, no special calls of fcn with
    !!--..             iflag = 0 are made.
    !!--..
    !!--..    info     is an integer output variable.  If the user has terminated execution,
    !!--..             info is set to the (negative) value of iflag. See description of fcn.  Otherwise, info is set as follows.
    !!--..
    !!--..             info = 0  improper input parameters.
    !!--..
    !!--..             info = 1  both actual and predicted relative reductions
    !!--..                       in the sum of squares are at most ftol.
    !!--..
    !!--..             info = 2  relative error between two consecutive iterates
    !!--..                       is at most xtol.
    !!--..
    !!--..             info = 3  conditions for info = 1 and info = 2 both hold.
    !!--..
    !!--..             info = 4  the cosine of the angle between fvec and any column of
    !!--..                       the jacobian is at most gtol in absolute value.
    !!--..
    !!--..             info = 5  number of calls to fcn with iflag = 1 has reached maxfev.
    !!--..
    !!--..             info = 6  ftol is too small.  No further reduction in
    !!--..                       the sum of squares is possible.
    !!--..
    !!--..             info = 7  xtol is too small.  No further improvement in
    !!--..                       the approximate solution x is possible.
    !!--..
    !!--..             info = 8  gtol is too small.  fvec is orthogonal to the
    !!--..                       columns of the jacobian to machine precision.
    !!--..
    !!--..    nfev     is an integer output variable set to the number of calls to fcn with
    !!--..             iflag = 1.
    !!--..
    !!--..    njev     is an integer output variable set to the number of calls to fcn with
    !!--..             iflag = 2.
    !!--..
    !!--..    ipvt     is an integer output array of length n.  ipvt defines a permutation
    !!--..             matrix p such that jac*p = q*r, where jac is the final calculated
    !!--..             jacobian, q is orthogonal (not stored), and r is upper triangular
    !!--..             with diagonal elements of nonincreasing magnitude. column j of p is
    !!--..             column ipvt(j) of the identity matrix.
    !!--..
    !!--..    qtf      is an output array of length n which contains the first n elements
    !!--..             of the vector (q transpose)*fvec.
    !!--..
    !!--..    wa1..3   are work arrays of length n.
    !!--..
    !!--..    wa4      is a work array of length m.
    !!--..
    !!--..  Argonne National Laboratory. minpack project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!--++
    !!--++
    !!--++ Update: February - 2009
    !!
    Subroutine lmder(fcn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev, &
                     mode, factor, nprint, info, nfev, njev, ipvt)
       !---- Arguments ----!
       Integer,                        Intent(In)      :: m
       Integer,                        Intent(In)      :: n
       Real (Kind=cp), Dimension(:),   Intent(In Out)  :: x
       Real (Kind=cp), Dimension(m),   Intent(Out)     :: fvec
       Real (Kind=cp), Dimension(:,:), Intent(Out)     :: fjac    ! fjac(ldfjac,n)
       Real (Kind=cp),                 Intent(In)      :: ftol
       Real (Kind=cp),                 Intent(In)      :: xtol
       Real (Kind=cp),                 Intent(In Out)  :: gtol
       Integer,                        Intent(In Out)  :: maxfev
       Integer,                        Intent(In)      :: mode
       Real (Kind=cp),                 Intent(In)      :: factor
       Integer,                        Intent(In)      :: nprint
       Integer,                        Intent(Out)     :: info
       Integer,                        Intent(Out)     :: nfev
       Integer,                        Intent(Out)     :: njev
       Integer,        Dimension(:),   Intent(Out)     :: ipvt

       Interface
         Subroutine fcn(m, n, x, fvec, fjac, iflag)
            Use CFML_GlobalDeps, Only: cp
            Integer,                       Intent(In)      :: m, n
            Real (Kind=cp),Dimension(:),   Intent(In)      :: x
            Real (Kind=cp),Dimension(:),   Intent(In Out)  :: fvec
            Real (Kind=cp),Dimension(:,:), Intent(Out)     :: fjac
            Integer,                       Intent(In Out)  :: iflag
         End Subroutine fcn
       End Interface

       !--- Local Variables ---!
       Integer                      :: i, iflag, iter, j, l
       Real (Kind=cp)               :: actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm,  &
                                       par, pnorm, prered, ratio, suma, temp, temp1, temp2, xnorm
       Real (Kind=cp), dimension(n) :: diag, qtf, wa1, wa2, wa3
       Real (Kind=cp), dimension(m) :: wa4
       Real (Kind=cp), Parameter    :: one = 1.0_CP, p1 = 0.1_CP, p5 = 0.5_CP,  &
                                       p25 = 0.25_CP, p75 = 0.75_CP, p0001 = 0.0001_CP, &
                                       zero = 0.0_CP
       Logical                      :: lmerror

       !     epsmch is the machine precision.
       epsmch = Epsilon(zero)
       lmerror=.false.
       info = 0
       iflag = 0
       nfev = 0
       njev = 0

       !     check the input parameters for errors.
       If (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
           .OR. maxfev <= 0 .OR. factor <= zero) lmerror=.true.

       If (mode == 2 .and. .not. lmerror) Then
         Do j = 1, n
            If (diag(j) <= zero) Then
               lmerror=.true.
               exit
            End If
         End Do
       End If

       !     Evaluate the function at the starting point and calculate its norm.
       If (.not. lmerror) then
         iflag = 1
         Call fcn(m, n, x, fvec, fjac, iflag)
         nfev = 1
         If (iflag < 0) lmerror=.true.
       end if

       If (.not. lmerror) then
         fnorm = enorm(m, fvec)

         !     Initialize Levenberg-Marquardt parameter and iteration counter.
         par = zero
         iter = 1
         Do_Out: Do !     Beginning of the outer loop.

            !        Calculate the Jacobian matrix.
            iflag = 2
            Call fcn(m, n, x, fvec, fjac, iflag)
            njev = njev + 1
            If (iflag < 0) Exit Do_Out

            !        If Requested, call fcn to enable printing of iterates.
            If (nprint > 0) Then
               iflag = 0
               If (Mod(iter-1,nprint) == 0) Call fcn(m, n, x, fvec, fjac, iflag)
               If (iflag < 0) Exit Do_Out
            End If

            !        Compute the qr factorization of the jacobian.
            Call qrfac(m, n, fjac, .TRUE., ipvt, wa1, wa2)

            !        On the first iteration and if mode is 1, scale according
            !        to the norms of the columns of the initial jacobian.
            If (iter == 1) Then
               If (mode /= 2) Then
                  Do j = 1, n
                     diag(j) = wa2(j)
                     If (wa2(j) == zero) diag(j) = one
                  End Do
               End If

               !        On the first iteration, calculate the norm of the scaled x
               !        and initialize the step bound delta.
               wa3(1:n) = diag(1:n)*x(1:n)
               xnorm = enorm(n,wa3)
               delta = factor*xnorm
               If (delta == zero) delta = factor
            End If

            !        Form (q transpose)*fvec and store the first n components in qtf.
            wa4(1:m) = fvec(1:m)
            Do j = 1, n
               If (fjac(j,j) /= zero) Then
                  suma = Dot_Product( fjac(j:m,j), wa4(j:m) )
                  temp = -suma/fjac(j,j)
                  Do i = j, m
                     wa4(i) = wa4(i) + fjac(i,j)*temp
                  End Do
               End If
               fjac(j,j) = wa1(j)
               qtf(j) = wa4(j)
            End Do

            !        Compute the norm of the scaled gradient.
            gnorm = zero
            If (fnorm /= zero) Then
               Do j = 1, n
                  l = ipvt(j)
                  If (wa2(l) == zero) Cycle
                  suma = zero
                  Do i = 1, j
                     suma = suma + fjac(i,j)*(qtf(i)/fnorm)
                  End Do
                  gnorm = Max(gnorm,Abs(suma/wa2(l)))
               End Do
            End If

            !        Test for convergence of the gradient norm.
            If (gnorm <= gtol) info = 4
            If (info /= 0) Exit Do_Out

            !        Rescale if necessary.
            If (mode /= 2) Then
               Do j = 1, n
                  diag(j) = Max(diag(j), wa2(j))
               End Do
            End If

            !        Beginning of the inner loop.
            Do_In: Do

               !           Determine the Levenberg-Marquardt parameter.
               Call lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)

               !           Store the direction p and x + p. calculate the norm of p.
               Do j = 1, n
                  wa1(j) = -wa1(j)
                  wa2(j) = x(j) + wa1(j)
                  wa3(j) = diag(j)*wa1(j)
               End Do
               pnorm = enorm(n, wa3)

               !           On the first iteration, adjust the initial step bound.
               If (iter == 1) delta = Min(delta,pnorm)

               !           Evaluate the function at x + p and calculate its norm.
               iflag = 1
               Call fcn(m, n, wa2, wa4, fjac, iflag)
               nfev = nfev + 1
               If (iflag < 0) Exit Do_Out
               fnorm1 = enorm(m, wa4)

               !           Compute the scaled actual reduction.
               actred = -one
               If (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

               !           Compute the scaled predicted reduction and
               !           the scaled directional derivative.
               Do j = 1, n
                  wa3(j) = zero
                  l = ipvt(j)
                  temp = wa1(l)
                  wa3(1:j) = wa3(1:j) + fjac(1:j,j)*temp
               End Do
               temp1 = enorm(n,wa3)/fnorm
               temp2 = (Sqrt(par)*pnorm)/fnorm
               prered = temp1**2 + temp2**2/p5
               dirder = -(temp1**2 + temp2**2)

               !           Compute the ratio of the actual to the predicted reduction.
               ratio = zero
               If (prered /= zero) ratio = actred/prered

               !           Update the step bound.
               If (ratio <= p25) Then
                  If (actred >= zero) temp = p5
                  If (actred < zero) temp = p5*dirder/(dirder + p5*actred)
                  If (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
                  delta = temp*Min(delta, pnorm/p1)
                  par = par/temp
               Else
                  If (.not. (par /= zero .AND. ratio < p75) ) Then
                     delta = pnorm/p5
                     par = p5*par
                  End If
               End If

               !           Test for successful iteration.
               If (ratio >= p0001) Then
                  !           Successful iteration. Update x, fvec, and their norms.
                  Do j = 1, n
                     x(j) = wa2(j)
                     wa2(j) = diag(j)*x(j)
                  End Do
                  fvec(1:m) = wa4(1:m)
                  xnorm = enorm(n,wa2)
                  fnorm = fnorm1
                  iter = iter + 1
               End If

               !           Tests for convergence.

               If (Abs(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
               If (delta <= xtol*xnorm) info = 2
               If (Abs(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one .AND. info == 2) info = 3
               If (info /= 0) Exit Do_Out

               !           Tests for termination and stringent tolerances.
               If (nfev >= maxfev) info = 5
               If (Abs(actred) <= epsmch .AND. prered <= epsmch .AND. p5*ratio <= one) info = 6
               If (delta <= epsmch*xnorm) info = 7
               If (gnorm <= epsmch) info = 8
               If (info /= 0) Exit Do_Out

               !           End of the inner loop. repeat if iteration unsuccessful.
               If (ratio >= p0001) Exit Do_In
            End Do Do_In
            !        End of the outer loop.
         End Do Do_Out
       End If

       !     Termination, either normal or user imposed.
       If (iflag < 0) info = iflag
       iflag = 0
       If (nprint > 0) Call fcn(m, n, x, fvec, fjac, iflag)

       Return
    End Subroutine lmder

    !!--++
    !!--++   Subroutine lmder1(Model_Functn, m, n, x, fvec, fjac, tol, nprint, info, ipvt)
    !!--++     Integer,                        Intent(In)      :: m
    !!--++     Integer,                        Intent(In)      :: n
    !!--++     Real (Kind=cp),Dimension(:),    Intent(In Out)  :: x
    !!--++     Real (Kind=cp),Dimension(:),    Intent(Out)     :: fvec
    !!--++     Real (Kind=cp),Dimension(:,:),  Intent(In Out)  :: fjac
    !!--++     Real (Kind=cp),                 Intent(In)      :: tol
    !!--++     Integer,                        Intent(In)      :: nprint
    !!--++     Integer,                        Intent(Out)     :: info
    !!--++     Integer,Dimension(:),           Intent(In Out)  :: ipvt
    !!--++     !--- Local Variables ---!
    !!--++     Interface
    !!--++       Subroutine Model_Functn(m, n, x, fvec, fjac, iflag)
    !!--++         Use CFML_GlobalDeps, Only: cp
    !!--++         Integer,                       Intent(In)    :: m, n
    !!--++         Real (Kind=cp),Dimension(:),   Intent(In)    :: x
    !!--++         Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec
    !!--++         Real (Kind=cp),Dimension(:,:), Intent(Out)   :: fjac
    !!--++         Integer,                       Intent(In Out):: iflag
    !!--++       End Subroutine Model_Functn
    !!--++     End Interface
    !!--++
    !!--++     The purpose of lmder1 is to minimize the sum of the squares of m nonlinear
    !!--++     functions in n variables by a modification of the levenberg-marquardt
    !!--++     algorithm.  This is done by using the more general least-squares solver
    !!--++     lmder.  The user must provide a subroutine which calculates the functions
    !!--++     and the jacobian.
    !!--++
    !!--..     Original documentation
    !!--..
    !!--..     The subroutine statement is
    !!--..
    !!--..         subroutine lmder1(fcn, m, n, x, fvec, fjac, tol, nprint, info, ipvt)
    !!--..
    !!--..     where
    !!--..
    !!--..     fcn     is the name of the user-supplied subroutine which calculates the
    !!--..             functions and the jacobian.  fcn must be declared in an interface
    !!--..             statement in the user calling program, and should be written as
    !!--..             follows.
    !!--..
    !!--..             subroutine fcn(m, n, x, fvec, fjac, iflag)
    !!--..                integer   :: m, n, ldfjac, iflag
    !!--..                real (kind=cp) :: x(:), fvec(:), fjac(:,:)
    !!--..                ----------
    !!--..                if iflag = 1 calculate the functions at x and
    !!--..                return this vector in fvec. do not alter fjac.
    !!--..                if iflag = 2 calculate the jacobian at x and
    !!--..                return this matrix in fjac. do not alter fvec.
    !!--..                ----------
    !!--..                return
    !!--..             end subroutine
    !!--..
    !!--..             The value of iflag should not be changed by fcn unless the user
    !!--..             wants to terminate execution of LM_Der. In this case set iflag
    !!--..             to a negative integer.
    !!--..
    !!--..     m       is a positive integer input variable set to the number of functions.
    !!--..
    !!--..     n       is a positive integer input variable set to the number of variables.
    !!--..             n must not exceed m.
    !!--..
    !!--..     x       is an array of length n. on input x must contain an initial estimate
    !!--..             of the solution vector. on output x contains the final estimate of
    !!--..             the solution vector.
    !!--..
    !!--..     fvec    is an output array of length m which contains the functions evaluated
    !!--..             at the output x.
    !!--..
    !!--..     fjac    is an output m by n array. the upper n by n submatrix of fjac
    !!--..             contains an upper triangular matrix r with diagonal elements of
    !!--..             non increasing magnitude such that
    !!--..
    !!--..                  t     t           t
    !!--..                 p *(jac *jac)*p = r *r,
    !!--..
    !!--..             where p is a permutation matrix and jac is the final calculated
    !!--..             Jacobian.  Column j of p is column ipvt(j) (see below) of the
    !!--..             identity matrix.  The lower trapezoidal part of fjac contains
    !!--..             information generated during the computation of r.
    !!--..
    !!--..    ldfjac   is a positive integer input variable not less than m which specifies
    !!--..             the leading dimension of the array fjac.
    !!--..
    !!--..    tol      is a nonnegative input variable. termination occurs when the
    !!--..             algorithm estimates either that the relative error in the sum of
    !!--..             squares is at most tol or that the relative error between x and the
    !!--..             solution is at most tol.
    !!--..
    !!--..    info     is an integer output variable.  If the user has terminated execution,
    !!--..             info is set to the (negative) value of iflag. See description of
    !!--..             fcn.  Otherwise, info is set as follows.
    !!--..
    !!--..             info = 0  improper input parameters.
    !!--..
    !!--..             info = 1  algorithm estimates that the relative error
    !!--..                       in the sum of squares is at most tol.
    !!--..
    !!--..             info = 2  algorithm estimates that the relative error
    !!--..                       between x and the solution is at most tol.
    !!--..
    !!--..             info = 3  conditions for info = 1 and info = 2 both hold.
    !!--..
    !!--..             info = 4  fvec is orthogonal to the columns of the
    !!--..                       jacobian to machine precision.
    !!--..
    !!--..             info = 5  number of calls to fcn with iflag = 1 has reached 100*(n+1).
    !!--..
    !!--..             info = 6  tol is too small.  No further reduction in
    !!--..                       the sum of squares is possible.
    !!--..
    !!--..             info = 7  tol is too small.  No further improvement in
    !!--..                       the approximate solution x is possible.
    !!--..
    !!--..    ipvt     is an integer output array of length n. ipvt defines a permutation
    !!--..             matrix p such that jac*p = q*r, where jac is the final calculated
    !!--..             jacobian, q is orthogonal (not stored), and r is upper triangular
    !!--..             with diagonal elements of nonincreasing magnitude. column j of p
    !!--..             is column ipvt(j) of the identity matrix.
    !!--..
    !!--..    wa       is a work array of length lwa.
    !!--..
    !!--..    lwa      is a positive integer input variable not less than 5*n+m.
    !!--..
    !!--..  Argonne National Laboratory. minpack project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!--++
    !!--++ Update: February - 2009
    !!
    Subroutine lmder1(Model_Functn, m, n, x, fvec, fjac, tol, nprint, info, ipvt)
       !---- Arguments ----!
       Integer,                        Intent(In)      :: m
       Integer,                        Intent(In)      :: n
       Real (Kind=cp),Dimension(:),    Intent(In Out)  :: x
       Real (Kind=cp),Dimension(:),    Intent(Out)     :: fvec
       Real (Kind=cp),Dimension(:,:),  Intent(In Out)  :: fjac
       Real (Kind=cp),                 Intent(In)      :: tol
       Integer,                        Intent(in)      :: nprint
       Integer,                        Intent(Out)     :: info
       Integer,Dimension(:),           Intent(In Out)  :: ipvt

       Interface
         Subroutine Model_Functn(m, n, x, fvec, fjac, iflag)
            Use CFML_GlobalDeps, Only: cp
            Integer,                       Intent(In)    :: m, n
            Real (Kind=cp),Dimension(:),   Intent(In)    :: x
            Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec
            Real (Kind=cp),Dimension(:,:), Intent(Out)   :: fjac
            Integer,                       Intent(In Out):: iflag
         End Subroutine Model_Functn
       End Interface

       !--- Local Variables ---!
       Integer                   :: maxfev, mode, nfev, njev
       Real (Kind=cp)            :: ftol, gtol, xtol
       Real (Kind=cp), Parameter :: factor = 100.0_CP, zero = 0.0_CP

       info = 0
       !     check the input parameters for errors.
       If ( n <= 0 .OR. m < n .OR. tol < zero ) Return
       maxfev = 100*(n + 1)
       ftol = tol
       xtol = tol
       gtol = zero
       mode = 1
       Call lmder(Model_Functn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev,  &
                  mode, factor, nprint, info, nfev, njev, ipvt)
       If (info == 8) info = 4

       Return
    End Subroutine lmder1

    !!--++
    !!--++   Subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,  &
    !!--++                    mode, factor, nprint, info, nfev, fjac, ipvt)
    !!--++     INTEGER,                       INTENT(IN)    :: m
    !!--++     INTEGER,                       INTENT(IN)    :: n
    !!--++     REAL (kind=cp),dimension(:),   INTENT(IN OUT):: x
    !!--++     REAL (kind=cp),dimension(:),   INTENT(OUT)   :: fvec
    !!--++     REAL (kind=cp),                INTENT(IN)    :: ftol
    !!--++     REAL (kind=cp),                INTENT(IN)    :: xtol
    !!--++     REAL (kind=cp),                INTENT(IN OUT):: gtol
    !!--++     INTEGER,                       INTENT(IN OUT):: maxfev
    !!--++     REAL (kind=cp),                INTENT(IN OUT):: epsfcn
    !!--++     INTEGER,                       INTENT(IN)    :: mode
    !!--++     REAL (kind=cp),                INTENT(IN)    :: factor
    !!--++     INTEGER,                       INTENT(IN)    :: nprint
    !!--++     INTEGER,                       INTENT(OUT)   :: info
    !!--++     INTEGER,                       INTENT(OUT)   :: nfev
    !!--++     REAL (kind=cp),dimension(:,:), INTENT(OUT)   :: fjac     ! fjac(ldfjac,n)
    !!--++     INTEGER,dimension(:),          INTENT(OUT)   :: ipvt(:)
    !!--++
    !!--++     INTERFACE
    !!--++       SUBROUTINE fcn(m, n, x, fvec, iflag)
    !!--++         Use CFML_GlobalDeps, only: cp
    !!--++         INTEGER,                     INTENT(IN)      :: m, n
    !!--++         REAL (kind=cp),dimension(:), INTENT(IN)      :: x
    !!--++         REAL (kind=cp),dimension(:), INTENT(IN OUT)  :: fvec
    !!--++         INTEGER,                     INTENT(IN OUT)  :: iflag
    !!--++       END SUBROUTINE fcn
    !!--++     END INTERFACE
    !!--++
    !!--++
    !!--..  Original documentation
    !!--..
    !!--..  The purpose of lmdif is to minimize the sum of the squares of m nonlinear
    !!--..  functions in n variables by a modification of the Levenberg-Marquardt
    !!--..  algorithm.  The user must provide a subroutine which calculates the
    !!--..  functions.  The jacobian is then calculated by a forward-difference
    !!--..  approximation.
    !!--..
    !!--..  The subroutine statement is
    !!--..
    !!--..    subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,
    !!--..                     diag, mode, factor, nprint, info, nfev, fjac,
    !!--..                     ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)
    !!--..
    !!--..  where
    !!--..
    !!--..    fcn       is the name of the user-supplied subroutine which calculates the
    !!--..              functions.  fcn must be declared in an external statement in the
    !!--..              user calling program, and should be written as follows.
    !!--..
    !!--..              subroutine fcn(m, n, x, fvec, iflag)
    !!--..                 integer m, n, iflag
    !!--..                 real (kind=cp) x(:), fvec(m)
    !!--..                 ----------
    !!--..                 calculate the functions at x and return this vector in fvec.
    !!--..                 ----------
    !!--..                 return
    !!--..              end
    !!--..
    !!--..              The value of iflag should not be changed by fcn unless the user
    !!--..              wants to terminate execution of lmdif. In this case set iflag to
    !!--..              a negative integer.
    !!--..
    !!--..    m         is a positive integer input variable set to the number of functions.
    !!--..
    !!--..    n         is a positive integer input variable set to the number of variables.
    !!--..              n must not exceed m.
    !!--..
    !!--..    x         is an array of length n.  On input x must contain an initial
    !!--..              estimate of the solution vector. On output x contains the final
    !!--..              estimate of the solution vector.
    !!--..
    !!--..    fvec      is an output array of length m which contains the functions evaluated
    !!--..              at the output x.
    !!--..
    !!--..    ftol      is a nonnegative input variable. Termination occurs when both the
    !!--..              actual and predicted relative reductions in the sum of squares are
    !!--..              at most ftol.  Therefore, ftol measures the relative error desired
    !!--..              in the sum of squares.
    !!--..
    !!--..    xtol      is a nonnegative input variable.  Termination occurs when both the
    !!--..              actual and predicted relative reductions in the sum of squares are
    !!--..              at most ftol.  Therefore, ftol measures the relative error desired
    !!--..              in the sum of squares.
    !!--..
    !!--..    gtol      is a nonnegative input variable.  Termination occurs when both the
    !!--..              actual and predicted relative reductions in the sum of squares are
    !!--..              at most ftol.  Therefore, ftol measures the relative error desired
    !!--..              in the sum of squares.
    !!--..
    !!--..    maxfev    is a positive integer input variable.  Termination occurs when the
    !!--..              number of calls to fcn is at least maxfev by the end of an iteration.
    !!--..
    !!--..    epsfcn    is an input variable used in determining a suitable step length
    !!--..              for the forward-difference approximation.  This approximation
    !!--..              assumes that the relative errors in the functions are of the order
    !!--..              of epsfcn.
    !!--..              If epsfcn is less than the machine precision, it is assumed that
    !!--..              the relative errors in the functions are of the order of the
    !!--..              machine precision.
    !!--..
    !!--..    diag      is an array of length n.  If mode = 1 (see below), diag is
    !!--..              internally set.  If mode = 2, diag must contain positive entries
    !!--..              that serve as multiplicative scale factors for the variables.
    !!--..
    !!--..    mode      is an integer input variable.  If mode = 1, the variables will be
    !!--..              scaled internally.  If mode = 2, the scaling is specified by the
    !!--..              input diag. other values of mode are equivalent to mode = 1.
    !!--..
    !!--..    factor    is a positive input variable used in determining the initial step
    !!--..              bound.  This bound is set to the product of factor and the euclidean
    !!--..              norm of diag*x if nonzero, or else to factor itself.  In most cases
    !!--..              factor should lie in the interval (.1,100.). 100. is a generally
    !!--..              recommended value.
    !!--..
    !!--..    nprint    is an integer input variable that enables controlled printing of
    !!--..              iterates if it is positive.  In this case, fcn is called with
    !!--..              iflag = 0 at the beginning of the first iteration and every nprint
    !!--..              iterations thereafter and immediately prior to return, with x and
    !!--..              fvec available for printing.  If nprint is not positive, no special
    !!--..              calls of fcn with iflag = 0 are made.
    !!--..
    !!--..    info      is an integer output variable.  If the user has terminated execution,
    !!--..              info is set to the (negative) value of iflag. See description of fcn.
    !!--..              Otherwise, info is set as follows.
    !!--..
    !!--..              info = 0  improper input parameters.
    !!--..
    !!--..              info = 1  both actual and predicted relative reductions
    !!--..                        in the sum of squares are at most ftol.
    !!--..
    !!--..              info = 2  relative error between two consecutive iterates <= xtol.
    !!--..
    !!--..              info = 3  conditions for info = 1 and info = 2 both hold.
    !!--..
    !!--..              info = 4  the cosine of the angle between fvec and any column of
    !!--..                        the Jacobian is at most gtol in absolute value.
    !!--..
    !!--..              info = 5  number of calls to fcn has reached or exceeded maxfev.
    !!--..
    !!--..              info = 6  ftol is too small. no further reduction in
    !!--..                        the sum of squares is possible.
    !!--..
    !!--..              info = 7  xtol is too small. no further improvement in
    !!--..                        the approximate solution x is possible.
    !!--..
    !!--..              info = 8  gtol is too small. fvec is orthogonal to the
    !!--..                        columns of the jacobian to machine precision.
    !!--..
    !!--..    nfev      is an integer output variable set to the number of calls to fcn.
    !!--..
    !!--..    fjac      is an output m by n array. the upper n by n submatrix of fjac
    !!--..              contains an upper triangular matrix r with diagonal elements of
    !!--..              nonincreasing magnitude such that
    !!--..
    !!--..                      t     t           t
    !!--..                     p *(jac *jac)*p = r *r,
    !!--..
    !!--..              where p is a permutation matrix and jac is the final calculated
    !!--..              Jacobian.  Column j of p is column ipvt(j) (see below) of the
    !!--..              identity matrix. the lower trapezoidal part of fjac contains
    !!--..              information generated during the computation of r.
    !!--..
    !!--..    ldfjac    is a positive integer input variable not less than m which i
    !!--..              specifies the leading dimension of the array fjac.
    !!--..
    !!--..    ipvt      is an integer output array of length n.  ipvt defines a permutation
    !!--..              matrix p such that jac*p = q*r, where jac is the final calculated
    !!--..              jacobian, q is orthogonal (not stored), and r is upper triangular
    !!--..              with diagonal elements of nonincreasing magnitude.
    !!--..              Column j of p is column ipvt(j) of the identity matrix.
    !!--..
    !!--..    qtf       is an output array of length n which contains the first n elements i
    !!--..              of the vector (q transpose)*fvec.
    !!--..
    !!--..    wa1...3   are work arrays of length n.
    !!--..
    !!--..    wa4       is a work array of length m.
    !!--..
    !!--..  Argonne National Laboratory. minpack project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!--++
    !!--++ Updated: February - 2009
    !!
    Subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,  &
                     mode, factor, nprint, info, nfev, fjac, ipvt)
       !---- Arguments ----!
       Integer,                        Intent(In)    :: m
       Integer,                        Intent(In)    :: n
       Real (Kind=cp), Dimension(:),   Intent(In Out):: x
       Real (Kind=cp), Dimension(:),   Intent(Out)   :: fvec
       Real (Kind=cp),                 Intent(In)    :: ftol
       Real (Kind=cp),                 Intent(In)    :: xtol
       Real (Kind=cp),                 Intent(In Out):: gtol
       Integer,                        Intent(In Out):: maxfev
       Real (Kind=cp),                 Intent(In Out):: epsfcn
       Integer,                        Intent(In)    :: mode
       Real (Kind=cp),                 Intent(In)    :: factor
       Integer,                        Intent(In)    :: nprint
       Integer,                        Intent(Out)   :: info
       Integer,                        Intent(Out)   :: nfev
       Real (Kind=cp), Dimension(:,:), Intent(Out)   :: fjac     ! fjac(ldfjac,n)
       Integer,        Dimension(:),   Intent(Out)   :: ipvt(:)

       Interface
         Subroutine fcn(m, n, x, fvec, iflag)
            Use CFML_GlobalDeps, Only: cp
            Integer,                     Intent(In)      :: m, n
            Real (Kind=cp),Dimension(:), Intent(In)      :: x
            Real (Kind=cp),Dimension(:), Intent(In Out)  :: fvec
            Integer,                     Intent(In Out)  :: iflag
         End Subroutine fcn
       End Interface

       !--- Local Variables ---!
       Integer                     :: i, iflag, iter, j, l
       Real (Kind=cp)              :: actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm,  &
                                      par, pnorm, prered, ratio, suma, temp, temp1, temp2, xnorm
       Real (Kind=cp), Dimension(n):: diag, qtf, wa1, wa2, wa3
       Real (Kind=cp), Dimension(m):: wa4
       Real (Kind=cp), Parameter   :: one = 1.0_CP, p1 = 0.1_CP, p5 = 0.5_CP,  &
                                      p25 = 0.25_CP, p75 = 0.75_CP, p0001 = 0.0001_CP, &
                                      zero = 0.0_CP
       Logical                     :: lmerror

       !     epsmch is the machine precision.
       epsmch = Epsilon(zero)
       lmerror=.false.
       info = 0
       iflag = 0
       nfev = 0

       !     check the input parameters for errors.
       If (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
           .OR. maxfev <= 0 .OR. factor <= zero) lmerror=.true.

       If (mode == 2 .and. .not. lmerror) Then
          Do j = 1, n
             If (diag(j) <= zero) Then
                lmerror=.true.
                exit
             End If
          End Do
       End If

       !     evaluate the function at the starting point and calculate its norm.
       if (.not. lmerror) then
          iflag = 1
          Call fcn(m, n, x, fvec, iflag)
          nfev = 1
          If (iflag < 0) lmerror=.true.
       end if

       if (.not. lmerror) then
          fnorm = enorm(m, fvec)

          !     initialize levenberg-marquardt parameter and iteration counter.
          par = zero
          iter = 1

          !     beginning of the outer loop.
          Do_out: Do  !

             !        calculate the jacobian matrix.
             iflag = 2
             Call fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)
             nfev = nfev + n
             If (iflag < 0) Exit Do_out

             !        If requested, call fcn to enable printing of iterates.
             If (nprint > 0) Then
                iflag = 0
                If (Mod(iter-1,nprint) == 0) Call fcn(m, n, x, fvec, iflag)
                If (iflag < 0) Exit Do_out
             End If

             !        Compute the qr factorization of the jacobian.
             Call qrfac(m, n, fjac, .TRUE., ipvt, wa1, wa2)

             !        On the first iteration and if mode is 1, scale according
             !        to the norms of the columns of the initial jacobian.
             If (iter == 1) Then
                If (mode /= 2) Then
                   Do j = 1, n
                      diag(j) = wa2(j)
                      If (wa2(j) == zero) diag(j) = one
                   End Do
                End If

                !        On the first iteration, calculate the norm of the scaled x
                !        and initialize the step bound delta.
                wa3(1:n) = diag(1:n)*x(1:n)
                xnorm = enorm(n, wa3)
                delta = factor*xnorm
                If (delta == zero) delta = factor
             End If

             !        Form (q transpose)*fvec and store the first n components in qtf.
             wa4(1:m) = fvec(1:m)
             Do j = 1, n
                If (fjac(j,j) /= zero) Then
                   suma = Dot_Product( fjac(j:m,j), wa4(j:m) )
                   temp = -suma/fjac(j,j)
                   Do i = j, m
                      wa4(i) = wa4(i) + fjac(i,j)*temp
                   End Do
                End If
                fjac(j,j) = wa1(j)
                qtf(j) = wa4(j)
             End Do

             !        compute the norm of the scaled gradient.
             gnorm = zero
             If (fnorm /= zero) Then
                Do j = 1, n
                   l = ipvt(j)
                   If (wa2(l) == zero) Cycle
                   suma = zero
                   Do i = 1, j
                      suma = suma + fjac(i,j)*(qtf(i)/fnorm)
                   End Do
                   gnorm = Max(gnorm, Abs(suma/wa2(l)))
                End Do
             End If

             !        test for convergence of the gradient norm.
             If (gnorm <= gtol) info = 4
             If (info /= 0) Exit Do_out

             !        rescale if necessary.
             If (mode /= 2) Then
                Do j = 1, n
                   diag(j) = Max(diag(j), wa2(j))
                End Do
             End If

             !        beginning of the inner loop.
             Do_In: Do

                !           determine the Levenberg-Marquardt parameter.
                Call lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)

                !           store the direction p and x + p. calculate the norm of p.
                Do j = 1, n
                   wa1(j) = -wa1(j)
                   wa2(j) = x(j) + wa1(j)
                   wa3(j) = diag(j)*wa1(j)
                End Do
                pnorm = enorm(n, wa3)

                !           on the first iteration, adjust the initial step bound.
                If (iter == 1) delta = Min(delta, pnorm)

                !           evaluate the function at x + p and calculate its norm.
                iflag = 1
                Call fcn(m, n, wa2, wa4, iflag)
                nfev = nfev + 1
                If (iflag < 0) Exit Do_Out
                fnorm1 = enorm(m, wa4)

                !           compute the scaled actual reduction.
                actred = -one
                If (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

                !           Compute the scaled predicted reduction and
                !           the scaled directional derivative.
                Do j = 1, n
                   wa3(j) = zero
                   l = ipvt(j)
                   temp = wa1(l)
                   Do i = 1, j
                      wa3(i) = wa3(i) + fjac(i,j)*temp
                   End Do
                End Do
                temp1 = enorm(n,wa3)/fnorm
                temp2 = (Sqrt(par)*pnorm)/fnorm
                prered = temp1**2 + temp2**2/p5
                dirder = -(temp1**2 + temp2**2)

                !           compute the ratio of the actual to the predicted reduction.
                ratio = zero
                If (prered /= zero) ratio = actred/prered

                !           update the step bound.
                If (ratio <= p25) Then
                   If (actred >= zero) temp = p5
                   If (actred < zero) temp = p5*dirder/(dirder + p5*actred)
                   If (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
                   delta = temp*Min(delta,pnorm/p1)
                   par = par/temp
                Else
                   If (.not. (par /= zero .AND. ratio < p75) ) Then
                      delta = pnorm/p5
                      par = p5*par
                   End If
                End If

                !           test for successful iteration.
                If (ratio >= p0001) Then

                   !           successful iteration. update x, fvec, and their norms.
                   Do j = 1, n
                      x(j) = wa2(j)
                      wa2(j) = diag(j)*x(j)
                   End Do
                   fvec(1:m) = wa4(1:m)
                   xnorm = enorm(n, wa2)
                   fnorm = fnorm1
                   iter = iter + 1
                End If

                !           tests for convergence.
                If (Abs(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
                If (delta <= xtol*xnorm) info = 2
                If (Abs(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one .AND. info == 2) info = 3
                If (info /= 0) Exit Do_out

                !           tests for termination and stringent tolerances.
                If (nfev >= maxfev) info = 5
                If (Abs(actred) <= epsmch .AND. prered <= epsmch .AND. p5*ratio <= one) info = 6
                If (delta <= epsmch*xnorm) info = 7
                If (gnorm <= epsmch) info = 8
                If (info /= 0) Exit Do_out

                !           end of the inner loop. repeat if iteration unsuccessful.
                If (ratio >= p0001) Exit Do_In
             End Do Do_in
             !        end of the outer loop.
          End Do  Do_out
       End If !.not. lmerror

       !     termination, either normal or user imposed.
       If (iflag < 0) info = iflag
       iflag = 0
       If (nprint > 0) Call fcn(m, n, x, fvec, iflag)

       Return
    End Subroutine lmdif

    !!--++
    !!--++  Subroutine lmdif1(fcn, m, n, x, fvec, tol, nprint, info, iwa)
    !!--++     Integer,                     Intent(In)      :: m
    !!--++     Integer,                     Intent(In)      :: n
    !!--++     Real (Kind=cp),Dimension(:), Intent(In Out)  :: x
    !!--++     Real (Kind=cp),Dimension(:), Intent(Out)     :: fvec
    !!--++     Real (Kind=cp),              Intent(In)      :: tol
    !!--++     Integer,                     Intent(In)      :: nprint
    !!--++     Integer,                     Intent(Out)     :: info
    !!--++     Integer,Dimension(:),        Intent(Out)     :: iwa
    !!--++
    !!--++     Interface
    !!--++       Subroutine fcn(m, n, x, fvec, iflag)
    !!--++         Use CFML_GlobalDeps, Only: cp
    !!--++         Integer,                     Intent(In)      :: m, n
    !!--++         Real (Kind=cp),Dimension(:), Intent(In)      :: x
    !!--++         Real (Kind=cp),Dimension(:), Intent(In Out)  :: fvec
    !!--++         Integer,                     Intent(In Out)  :: iflag
    !!--++       End Subroutine fcn
    !!--++     End Interface
    !!--++
    !!--..  Original documentation
    !!--..
    !!--..  The purpose of lmdif1 is to minimize the sum of the squares of m nonlinear
    !!--..  functions in n variables by a modification of the Levenberg-Marquardt
    !!--..  algorithm.  This is done by using the more general least-squares solver lmdif.
    !!--..  The user must provide a subroutine which calculates the functions.  The
    !!--..  jacobian is then calculated by a forward-difference approximation.
    !!--..
    !!--..  The subroutine statement is
    !!--..
    !!--..      subroutine lmdif1(fcn, m, n, x, fvec, tol, nprint, info, iwa)
    !!--..
    !!--..  where
    !!--..
    !!--..    fcn      is the name of the user-supplied subroutine which calculates
    !!--..             the functions.  fcn must be declared in an external statement in
    !!--..             the user calling program, and should be written as follows.
    !!--..
    !!--..             subroutine fcn(m, n, x, fvec, iflag)
    !!--..                integer m, n, iflag
    !!--..                real (kind=cp) x(n), fvec(m)
    !!--..                ----------
    !!--..                calculate the functions at x and return this vector in fvec.
    !!--..                ----------
    !!--..                return
    !!--..             end
    !!--..
    !!--..             The value of iflag should not be changed by fcn unless the user
    !!--..             wants to terminate execution of LM_Dif. In this case set iflag
    !!--..             to a negative integer.
    !!--..
    !!--..    m        is a positive integer input variable set to the number of functions.
    !!--..
    !!--..    n        is a positive integer input variable set to the number of variables.
    !!--..             n must not exceed m.
    !!--..
    !!--..    x        is an array of length n.  On input x must contain an initial
    !!--..             estimate of the solution vector.  On output x contains the final
    !!--..             estimate of the solution vector.
    !!--..
    !!--..    fvec     is an output array of length m which contains the functions evaluated
    !!--..             at the output x.
    !!--..
    !!--..    tol      is a nonnegative input variable.  Termination occurs when the
    !!--..             algorithm estimates either that the relative error in the sum of
    !!--..             squares is at most tol or that the relative error between x and the
    !!--..             solution is at most tol.
    !!--..
    !!--..    info     is an integer output variable.  If the user has terminated execution,
    !!--..             info is set to the (negative) value of iflag.  See description of
    !!--..             fcn. Otherwise, info is set as follows.
    !!--..
    !!--..             info = 0  improper input parameters.
    !!--..
    !!--..             info = 1  algorithm estimates that the relative error
    !!--..                       in the sum of squares is at most tol.
    !!--..
    !!--..             info = 2  algorithm estimates that the relative error
    !!--..                       between x and the solution is at most tol.
    !!--..
    !!--..             info = 3  conditions for info = 1 and info = 2 both hold.
    !!--..
    !!--..             info = 4  fvec is orthogonal to the columns of the
    !!--..                       jacobian to machine precision.
    !!--..
    !!--..             info = 5  number of calls to fcn has reached or exceeded 200*(n+1).
    !!--..
    !!--..             info = 6  tol is too small. no further reduction in
    !!--..                       the sum of squares is possible.
    !!--..
    !!--..             info = 7  tol is too small.  No further improvement in
    !!--..                       the approximate solution x is possible.
    !!--..
    !!--..    iwa      is an integer work array of length n.
    !!--..
    !!--..    wa       is a work array of length lwa.
    !!--..
    !!--..    lwa      is a positive integer input variable not less than m*n+5*n+m.
    !!--..
    !!--..  Argonne National Laboratory. minpack project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!--++
    !!--++ Update: February - 2009
    !!
    Subroutine lmdif1(Model_Functn, m, n, x, fvec, tol, nprint, info, iwa)
       !---- Argument ----!
       Integer,                     Intent(In)      :: m
       Integer,                     Intent(In)      :: n
       Real (Kind=cp),Dimension(:), Intent(In Out)  :: x
       Real (Kind=cp),Dimension(:), Intent(Out)     :: fvec
       Real (Kind=cp),              Intent(In)      :: tol
       Integer,                     Intent(In)      :: nprint
       Integer,                     Intent(Out)     :: info
       Integer,Dimension(:),        Intent(Out)     :: iwa

       Interface
         Subroutine Model_Functn(m, n, x, fvec, iflag)
            Use CFML_GlobalDeps, Only: cp
            Integer,                     Intent(In)      :: m, n
            Real (Kind=cp),Dimension(:), Intent(In)      :: x
            Real (Kind=cp),Dimension(:), Intent(In Out)  :: fvec
            Integer,                     Intent(In Out)  :: iflag
         End Subroutine Model_Functn
       End Interface

       !--- Local Variables ---!
       Integer                       :: maxfev, mode, nfev
       Real (Kind=cp)                :: epsfcn, ftol, gtol, xtol
       Real (Kind=cp),Dimension(m,n) :: fjac
       Real (Kind=cp), Parameter     :: factor = 100.0_CP, zero = 0.0_CP

       info = 0
       !     check the input parameters for errors.
       If (n <= 0 .OR. m < n .OR. tol < zero) Return
       !     call lmdif.
       maxfev = 200*(n + 1)
       ftol = tol
       xtol = tol
       gtol = zero
       epsfcn = zero
       mode = 1

       Call lmdif(Model_Functn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,   &
                  mode, factor, nprint, info, nfev, fjac, iwa)
       If (info == 8) info = 4

       return
    End Subroutine lmdif1

    !!--++
    !!--++   Subroutine lmpar(n, r, ipvt, diag, qtb, delta, par, x, sdiag)
    !!--++      Integer,                        Intent(In)      :: n
    !!--++      Real (Kind=cp), dimension(:,:), Intent(In Out)  :: r
    !!--++      Integer,        dimension(:),   Intent(In)      :: ipvt
    !!--++      Real (Kind=cp), dimension(:),   Intent(In)      :: diag
    !!--++      Real (Kind=cp), dimension(:),   Intent(In)      :: qtb
    !!--++      Real (Kind=cp),                 Intent(In)      :: delta
    !!--++      Real (Kind=cp),                 Intent(Out)     :: par
    !!--++      Real (Kind=cp), dimension(:),   Intent(Out)     :: x
    !!--++      Real (Kind=cp), dimension(:),   Intent(Out)     :: sdiag
    !!--++
    !!--..  Original documentation:
    !!--..
    !!--..  Given an m by n matrix a, an n by n nonsingular diagonal matrix d,
    !!--..  an m-vector b, and a positive number delta, the problem is to determine a
    !!--..  value for the parameter par such that if x solves the system
    !!--..
    !!--..        a*x = b ,     sqrt(par)*d*x = 0 ,
    !!--..
    !!--..  in the least squares sense, and dxnorm is the euclidean
    !!--..  norm of d*x, then either par is zero and
    !!--..
    !!--..        (dxnorm-delta) <= 0.1*delta ,
    !!--..
    !!--..  or par is positive and
    !!--..
    !!--..        abs(dxnorm-delta) <= 0.1*delta .
    !!--..
    !!--..  This subroutine completes the solution of the problem if it is provided
    !!--..  with the necessary information from the r factorization, with column
    !!--..  qpivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
    !!--..  q has orthogonal columns, and r is an upper triangular matrix with diagonal
    !!--..  elements of nonincreasing magnitude, then lmpar expects the full upper
    !!--..  triangle of r, the permutation matrix p, and the first n components of
    !!--..  (q transpose)*b.
    !!--..  On output lmpar also provides an upper triangular matrix s such that
    !!--..
    !!--..         t   t                   t
    !!--..        p *(a *a + par*d*d)*p = s *s .
    !!--..
    !!--..  s is employed within lmpar and may be of separate interest.
    !!--..
    !!--..  Only a few iterations are generally needed for convergence of the algorithm.
    !!--..  If, however, the limit of 10 iterations is reached, then the output par
    !!--..  will contain the best value obtained so far.
    !!--..
    !!--..  the subroutine statement is
    !!--..
    !!--..    subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag, wa1,wa2)
    !!--..
    !!--..  where
    !!--..
    !!--..    n is a positive integer input variable set to the order of r.
    !!--..
    !!--..    r is an n by n array. on input the full upper triangle
    !!--..      must contain the full upper triangle of the matrix r.
    !!--..      On output the full upper triangle is unaltered, and the
    !!--..      strict lower triangle contains the strict upper triangle
    !!--..      (transposed) of the upper triangular matrix s.
    !!--..
    !!--..    ldr is a positive integer input variable not less than n
    !!--..      which specifies the leading dimension of the array r.
    !!--..
    !!--..    ipvt is an integer input array of length n which defines the
    !!--..      permutation matrix p such that a*p = q*r. column j of p
    !!--..      is column ipvt(j) of the identity matrix.
    !!--..
    !!--..    diag is an input array of length n which must contain the
    !!--..      diagonal elements of the matrix d.
    !!--..
    !!--..    qtb is an input array of length n which must contain the first
    !!--..      n elements of the vector (q transpose)*b.
    !!--..
    !!--..    delta is a positive input variable which specifies an upper
    !!--..      bound on the euclidean norm of d*x.
    !!--..
    !!--..    par is a nonnegative variable. on input par contains an
    !!--..      initial estimate of the levenberg-marquardt parameter.
    !!--..      on output par contains the final estimate.
    !!--..
    !!--..    x is an output array of length n which contains the least
    !!--..      squares solution of the system a*x = b, sqrt(par)*d*x = 0,
    !!--..      for the output par.
    !!--..
    !!--..    sdiag is an output array of length n which contains the
    !!--..      diagonal elements of the upper triangular matrix s.
    !!--..
    !!--..    wa1...2 are work arrays of length n.
    !!--..
    !!--..  Argonne National Laboratory. minpack project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!--++
    !!--++
    !!--++ Update: February - 2009
    !!
    Subroutine lmpar(n, r, ipvt, diag, qtb, delta, par, x, sdiag)
       !---- Argument ----!
       Integer,                        Intent(In)      :: n
       Real (Kind=cp), dimension(:,:), Intent(In Out)  :: r
       Integer,        dimension(:),   Intent(In)      :: ipvt
       Real (Kind=cp), dimension(:),   Intent(In)      :: diag
       Real (Kind=cp), dimension(:),   Intent(In)      :: qtb
       Real (Kind=cp),                 Intent(In)      :: delta
       Real (Kind=cp),                 Intent(Out)     :: par
       Real (Kind=cp), dimension(:),   Intent(Out)     :: x
       Real (Kind=cp), dimension(:),   Intent(Out)     :: sdiag

       !--- Local Variables ---!
       Integer                      :: iter, j, jm1, jp1, k, l, nsing
       Real (Kind=cp)               :: dxnorm, dwarf, fp, gnorm, parc, parl, paru, suma, temp
       Real (Kind=cp), dimension(n) :: wa1, wa2
       Real (Kind=cp), Parameter    :: p1 = 0.1_CP, p001 = 0.001_CP, zero = 0.0_CP
       logical                      :: lmerror

       !     dwarf is the smallest positive magnitude.
       dwarf = Tiny(zero)
       lmerror=.false.

       !     compute and store in x the gauss-newton direction. if the
       !     jacobian is rank-deficient, obtain a least squares solution.
       nsing = n
       Do j = 1, n
          wa1(j) = qtb(j)
          If (r(j,j) == zero .AND. nsing == n) nsing = j - 1
          If (nsing < n) wa1(j) = zero
       End Do

       Do k = 1, nsing
          j = nsing - k + 1
          wa1(j) = wa1(j)/r(j,j)
          temp = wa1(j)
          jm1 = j - 1
          wa1(1:jm1) = wa1(1:jm1) - r(1:jm1,j)*temp
       End Do

       Do j = 1, n
          l = ipvt(j)
          x(l) = wa1(j)
       End Do

       !     initialize the iteration counter.
       !     evaluate the function at the origin, and test
       !     for acceptance of the gauss-newton direction.
       iter = 0
       wa2(1:n) = diag(1:n)*x(1:n)
       dxnorm = enorm(n, wa2)
       fp = dxnorm - delta
       If (fp <= p1*delta) lmerror=.true.
       If (.not. lmerror) then

          !     if the jacobian is not rank deficient, the newton
          !     step provides a lower bound, parl, for the zero of
          !     the function.  Otherwise set this bound to zero.
          parl = zero
          If (nsing >= n) Then
             Do j = 1, n
                l = ipvt(j)
                wa1(j) = diag(l)*(wa2(l)/dxnorm)
             End Do
             Do j = 1, n
                suma = Dot_Product( r(1:j-1,j), wa1(1:j-1) )
                wa1(j) = (wa1(j) - suma)/r(j,j)
             End Do
             temp = enorm(n,wa1)
             parl = ((fp/delta)/temp)/temp
          End If

          !     calculate an upper bound, paru, for the zero of the function.
          Do j = 1, n
             suma = Dot_Product( r(1:j,j), qtb(1:j) )
             l = ipvt(j)
             wa1(j) = suma/diag(l)
          End Do
          gnorm = enorm(n,wa1)
          paru = gnorm/delta
          If (paru == zero) paru = dwarf/Min(delta,p1)

          !     if the input par lies outside of the interval (parl,paru),
          !     set par to the closer endpoint.
          par = Max(par,parl)
          par = Min(par,paru)
          If (par == zero) par = gnorm/dxnorm

          !     beginning of an iteration.
          Do_Iter: Do
             iter = iter + 1

             !        evaluate the function at the current value of par.
             If (par == zero) par = Max(dwarf, p001*paru)
             temp = Sqrt(par)
             wa1(1:n) = temp*diag(1:n)
             Call qrsolv(n, r, ipvt, wa1, qtb, x, sdiag)
             wa2(1:n) = diag(1:n)*x(1:n)
             dxnorm = enorm(n, wa2)
             temp = fp
             fp = dxnorm - delta

             !        if the function is small enough, accept the current value
             !        of par. also test for the exceptional cases where parl
             !        is zero or the number of iterations has reached 10.
             If (Abs(fp) <= p1*delta .OR. parl == zero .AND. fp <= temp  &
                 .AND. temp < zero .OR. iter == 10) Exit Do_Iter

             !        compute the newton correction.
             Do j = 1, n
                l = ipvt(j)
                wa1(j) = diag(l)*(wa2(l)/dxnorm)
             End Do
             Do j = 1, n
                wa1(j) = wa1(j)/sdiag(j)
                temp = wa1(j)
                jp1 = j + 1
                wa1(jp1:n) = wa1(jp1:n) - r(jp1:n,j)*temp
             End Do
             temp = enorm(n,wa1)
             parc = ((fp/delta)/temp)/temp

             !        depending on the sign of the function, update parl or paru.
             If (fp > zero) parl = Max(parl,par)
             If (fp < zero) paru = Min(paru,par)

             !        compute an improved estimate for par.
             par = Max(parl, par+parc)

             !        end of an iteration.
          End Do Do_iter
          !     termination.
       End If ! .not. lmerror
       If (iter == 0) par = zero

       Return
    End Subroutine lmpar

    !!----
    !!---- Subroutine Marquardt_Fit(Model_Functn, X, Y, W, Yc, Nobs, c, vs, Ipr, Chi2, scroll_lines)
    !!----    real(kind=cp), dimension(:),             intent(in)     :: x             !Vector of x-values
    !!----    real(kind=cp), dimension(:),             intent(in)     :: y             !Vector of observed y-values
    !!----    real(kind=cp), dimension(:),             intent(in out) :: w             !Vector of weights-values (1/variance)
    !!----    real(kind=cp), dimension(:),             intent(   out) :: yc            !Vector of calculated y-values
    !!----    integer                    ,             intent(in)     :: nobs          !Number of effective components of x,y,w,yc
    !!----    type(LSQ_conditions_type),               intent(in out) :: c             !Conditions for the algorithm
    !!----    type(LSQ_State_Vector_type),             intent(in out) :: vs            !State vector for the model calculation
    !!----    integer                    ,             intent(in)     :: Ipr           !Logical unit for writing
    !!----    real(kind=cp),                           intent(out)    :: chi2          !Reduced Chi-2
    !!----    character(len=*),dimension(:), optional, intent(out)    :: scroll_lines  !If present, part of the output is stored
    !!----                                                                             !in this text for treatment in the calling program
    !!----
    !!----    Model_functn                                                             !Name of the subroutine calculating yc(i) for point x(i)
    !!----    Interface                                                                !Interface for the Model_Functn subroutine
    !!----       Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
    !!----          use CFML_GlobalDeps, only: cp
    !!----          integer,                             intent(in) :: iv     !Number of the component: "i" in x(i)
    !!----          real(kind=cp),                       intent(in) :: xv     !Value of x(i)
    !!----          real(kind=cp),                       intent(out):: ycalc  !Value of yc at point x(i) => ycalc=yc(i)
    !!----          real(kind=cp),dimension(:),          intent(in) :: aa     !Vector of parameters
    !!----          real(kind=cp),dimension(:),optional, intent(out):: der    !Derivatives of the function w.r.t. free parameters
    !!----       End Subroutine Model_Functn
    !!----    End Interface
    !!----
    !!----     or
    !!----
    !!---- Subroutine Marquardt_Fit(Model_Functn, d, c, vs, Ipr, Chi2, scroll_lines)
    !!----    Type(LSQ_Data_Type),                     intent(in out) :: d             !Vector of x-values
    !!----    type(LSQ_conditions_type),               intent(in out) :: c             !Conditions for the algorithm
    !!----    type(LSQ_State_Vector_type),             intent(in out) :: vs            !State vector for the model calculation
    !!----    integer                    ,             intent(in)     :: Ipr           !Logical unit for writing
    !!----    real(kind=cp),                           intent(out)    :: chi2          !Reduced Chi-2
    !!----    character(len=*),dimension(:), optional, intent(out)    :: scroll_lines  !If present, part of the output is stored
    !!----                                                                             !in this text for treatment in the calling program
    !!----
    !!----    Model_functn                                                             !Name of the subroutine calculating yc(i) for point x(i)
    !!----    Interface                                                                !Interface for the Model_Functn subroutine
    !!----       Subroutine Model_Functn(iv,Xv,ycalc,Vsa,calder)
    !!----          use CFML_GlobalDeps, only: cp
    !!----          integer,                    intent(in)     :: iv     !Number of the component: "i" in x(i)
    !!----          real(kind=cp),              intent(in)     :: xv     !Value of x(i)
    !!----          real(kind=cp),              intent(out)    :: ycalc  !Value of yc at point x(i) => ycalc=yc(i)
    !!----          Type(LSQ_State_Vector_type),intent(in out) :: Vsa    !Parameters, codes, and derivatives
    !!----          logical,optional,           intent(in)     :: calder !If present, derivatives, stored in Vsa%dpv, are calculated
    !!----       End Subroutine Model_Functn
    !!----    End Interface
    !!----
    !!----
    !!---- Update: January 2010
    !!

    !!--++
    !!--++ Subroutine Marquardt_Fit_v1(Model_Functn, X, Y, W, Yc, Nobs, c, vs, Ipr, Chi2, scroll_lines)
    !!--++    real(kind=cp), dimension(:),intent(in)     :: x      !Vector of x-values
    !!--++    real(kind=cp), dimension(:),intent(in)     :: y      !Vector of observed y-values
    !!--++    real(kind=cp), dimension(:),intent(in out) :: w      !Vector of weights-values (1/variance)
    !!--++    real(kind=cp), dimension(:),intent(   out) :: yc     !Vector of calculated y-values
    !!--++    integer                    ,intent(in)     :: nobs   !Number of effective components of x,y,w,yc
    !!--++    type(LSQ_conditions_type),  intent(in out) :: c      !Conditions for the algorithm
    !!--++    type(LSQ_State_Vector_type),intent(in out) :: vs     !State vector for the model calculation
    !!--++    integer                    ,intent(in)     :: Ipr    !Logical unit for writing
    !!--++    real(kind=cp),              intent(out)    :: chi2   !Reduced Chi-2
    !!--++    character(len=*),dimension(:), intent(out), optional  :: scroll_lines  !If present, part of the output is stored
    !!--++                                                                           !in this text for treatment in the calling program
    !!--++
    !!--++    Model_functn                                            !Name of the subroutine calculating yc(i) for point x(i)
    !!--++    Interface                                               !Interface for the Model_Functn subroutine
    !!--++       Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
    !!--++            use CFML_GlobalDeps, only: cp
    !!--++            integer,                             intent(in) :: iv     !Number of the component: "i" in x(i)
    !!--++            real(kind=cp),                       intent(in) :: xv     !Value of x(i)
    !!--++            real(kind=cp),                       intent(out):: ycalc  !Value of yc at point x(i) => ycalc=yc(i)
    !!--++            real(kind=cp),dimension(:),          intent(in) :: aa     !Vector of parameters
    !!--++            real(kind=cp),dimension(:),optional, intent(out):: der    !Derivatives of the function w.r.t. free parameters
    !!--++       End Subroutine Model_Functn
    !!--++    End Interface
    !!--++
    !!--++
    !!--++    Subroutine for applying the Levenberg-Marquardt method for Least-Squares.
    !!--++    The user must provide a model function according to the interface above.
    !!--++    The model function should use at least some of the public variables of the
    !!--++    present module in order to set the derivatives with respect to the model
    !!--++    parameters. Examples of using this module are given in the program
    !!--++    templates CW_fit and TOF_fit. This subroutine ignores vs%code_comp!
    !!--++
    !!--..    INFORMATION
    !!--..        Author: Juan Rodriguez-Carvajal (based in text descriptions of the literature)
    !!--..                Translated, extracted and modified from the Fortran 77 code: XRFIT (JRC 1986)
    !!--..        Module implementing the Marquardt method for non-linear least-squares.
    !!--..        For using this module, the user must provide:
    !!--..            1: Number of cycles (c%icyc), type of weighting scheme (c%iw), number of
    !!--..               model parameters (vs%np), constraint conditions (c%constr, c%percent)
    !!--..            2: A character(len=15) name for all possible parameters of the
    !!--..               model stored in the array vs%nampar(:).
    !!--..            3: A set of initial value for all parameters stored in vs%pv(:)
    !!--..            4: A set of flags values to refine or fix the parameters vs%code(:)
    !!--..            5: The model function (subroutine model_functn) must be provided by the user
    !!--..               according to the interface described below. The actual name of the model function
    !!--..               is arbitrary and it is passed to the only public procedure "marquardt_fit" as a
    !!--..               dummy argument.
    !!--..
    !!--..       The values of all possible refinable parameters are stored in the array vs%pv(:).
    !!--..       The derivatives must be calculated within model_functn, by using the array vs%der(:)
    !!--..       The actual refined parameters a(:) are selected from the larger vs%pv(:) array by
    !!--..       the integer array of flags: vs%code(:). A value vs%code(j) /= 0 means that the j-th
    !!--..       parameter is to be varied. A value vs%code(k) = 0 means that the k-th parameter is
    !!--..       kept fixed through the refinement cycles.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Marquardt_Fit_v1(Model_Functn,X,Y,W,Yc,Nobs,c,vs,Ipr,Chi2,scroll_lines)
       !---- Arguments ----!
       real(kind=cp),   dimension(:), intent(in)             :: x,y
       real(kind=cp),   dimension(:), intent(in out)         :: w
       real(kind=cp),   dimension(:), intent(   out)         :: yc
       integer,                       intent(in)             :: nobs,Ipr
       type(LSQ_conditions_type),     intent(in out)         :: c
       type(LSQ_State_Vector_type),   intent(in out)         :: vs
       real(kind=cp),                 intent(   out)         :: chi2
       character(len=*),dimension(:), intent(out), optional  :: scroll_lines

       Interface
        Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
           use CFML_GlobalDeps, only: cp
           integer,                             intent(in) :: iv
           real(kind=cp),                       intent(in) :: xv
           real(kind=cp),dimension(:),          intent(in) :: aa
           real(kind=cp),                       intent(out):: ycalc
           real(kind=cp),dimension(:),optional, intent(out):: der
        End Subroutine Model_Functn
       End Interface

       !---- local variables ----!
       logical                                :: fixed, write_cyc
       character(len=180)                     :: line
       integer                                :: ifail,nt
       integer                                :: i,j,ncount,li,ntex
       real(kind=cp)                          :: fl,chi1
       real(kind=cp), dimension(Max_Free_Par) :: a,sa

       fixed = .false.
       write_cyc=.true.
       if (c%icyc > 50) write_cyc=.false.

       !---- Beginning the iterations ----!
       fl=0.001
       chi1=1.0e+30
       c%reached =.false.
       write(unit=ipr,fmt="(/,a)")   " -------------------------------------------------"
       write(unit=ipr,fmt="(a,i3,a)")" => Marquardt Least Squares Fitting for ",c%icyc," cycles"
       write(unit=ipr,fmt="(a,/)")   " -------------------------------------------------"

       write(unit=ipr,fmt="(a,i6)")        " => Number of observations   : ",nobs
       ncount=0
       do i=1,vs%np
          if (vs%code(i) == 0) cycle
          ncount=ncount+1
          j=vs%code(i)
          namfree(j)=vs%nampar(i)
       end do
       c%npvar=ncount
       write(unit=ipr,fmt="(a,i6)")        " => Number of free parameters: ",ncount
       write(unit=ipr,fmt="(2(a,f12.4),a)")" => Range of variable X      : (",x(1),",",x(nobs),")"

       ntex=0
       if(present(scroll_lines)) scroll_lines=" "  !Clear the stored text

       do li=1,c%icyc                      !Loop over icyc cycles
          ncount=0
          do i=1,vs%np
             if (vs%code(i) == 0) cycle
             ncount=ncount+1
             a(ncount)=vs%pv(i)
          end do

          call curfit_v1(Model_Functn,x,y,w,nobs,c,a,sa,fl,yc,chi2,ifail,nt)

          if (ifail /= 0) then
             line= " "
             write(unit=line,fmt="(a,i3,a)") " => IFAIL /= 0 : ",ifail, " "//trim(ERR_Lsq_Mess)
             if(present(scroll_lines)) then
                ntex=ntex+1
                scroll_lines(ntex)=line
             else
                 if(ifail /= 1) write(unit=*,fmt="(a)") trim(line)
             end if
             if (ABS(chi2-chi1)*100/Chi2 <= 0.001 .and. ifail==1) then
                c%reached=.true.
                line= " "
                write(unit=line,fmt="(a)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.00001 * Chi2(old)"
                if(present(scroll_lines)) then
                    ntex=ntex+1
                    scroll_lines(ntex)=line
                else
                    write(unit=*,fmt="(a)") trim(line)
                end if
                write(unit=Ipr,fmt="(/,a,/)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.00001 * Chi2(old)"
                exit
             end if

             if(present(scroll_lines)) then
                 ntex=ntex+1
                 scroll_lines(ntex)=" => Marquardt convergence not reached (F-lambda increased too much, see output file)."
             else
                 write(unit=*,fmt="(a)") " => Marquardt convergence not reached (F-lambda increased too much, see output file)."
             end if

             write(unit=ipr,fmt="(/,a)")" => Marquardt convergence not reached (F-lambda increased too much)."
             write(unit=ipr,fmt="(a)" ) "    If Chi2 is reasonably low, you may have reached the minimum"
             write(unit=ipr,fmt="(a)" ) "    within the precision of the machine.                       "
             write(unit=ipr,fmt="(a)" ) "    Otherwise, look at the output file to search for the reason "
             write(unit=ipr,fmt="(a,/)")"    of anomalous behaviour and try to start from new parameters"

             write(unit=ipr,fmt="(/,a,/)")" =>  Diagonal elements of LSQ-matrix before inversion"
             do i=1,c%npvar
                write(unit=ipr,fmt="(a30,i3,a,i3,a,f16.6)") "    "//trim(namfree(i))//"  (",i,",",i,")", curv_mat(i,i)
             end do

             exit
          end if  !ifail/=0

          if (abs(chi2-chi1)*100.0/chi2 >= 0.0001 ) then
             c%reached=.false.
             chi1=chi2
             ncount=0
             do i=1,vs%np
                if (vs%code(i) == 0) then
                   vs%spv(i)=0.0
                   pn(i)=vs%pv(i)
                   ch(i)=0.0
                else
                   ncount=ncount+1
                   pn(i)=a(ncount)
                   ch(i)=pn(i)-vs%pv(i)
                   vs%spv(i)=sa(ncount)
                end if
             end do

             if (c%constr)  then
                call box_constraints(a,sa,fixed,c,vs)
                if (fixed) then
                   write(unit=ipr,fmt="(/,a,/,a,i2,a)") " => Some parameters have been restored to their initial value ", &
                                                        "    because of local divergence (change > ",nint(c%percent),"%)"
                   write(unit=ipr,fmt="(a,/)")          "    Look at the above list for parameters having a '*'"
                end if  !fixed
             end if

             line= " "
             write(unit=line,fmt="(a,i5,2(a,f14.6),a,i3)")" => Cycle No.:",li,"  Chi2 =",chi2, "  Marquardt-Lambda: ",fl, &
                                                          "  Ntrials: ",nt
             if(present(scroll_lines)) then
                 ntex=ntex+1
                 scroll_lines(ntex)=line
             else
                 write(unit=*,fmt="(a)") trim(line)
             end if

             if (write_cyc) call output_cyc(li,Ipr,chi2,vs)
             vs%pv(1:vs%np)=pn(1:vs%np)
             cycle
          else
             line= " "
             write(unit=line,fmt="(a)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.000001 * Chi2(old)"
             if(present(scroll_lines)) then
                 ntex=ntex+1
                 scroll_lines(ntex)=line
             else
                 write(unit=*,fmt="(a)") trim(line)
             end if
             write(unit=Ipr,fmt="(/,a,/)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.000001 * Chi2(old)"
             c%reached=.true.
             exit
          end if
       end do  !li=1,c%icyc

       if (c%reached) then
          if (fl <= 10.0) then
             fl=0.0
          else
             fl=0.5e-20
          end if
          call curfit_v1(Model_Functn,x,y,w,nobs,c,a,sa,fl,yc,chi2,ifail,nt)
          ncount=0
          do i=1,vs%np
             if (vs%code(i) == 0) then
                vs%spv(i)=0.0
                pn(i)=vs%pv(i)
                ch(i)=0.0
             else
                ncount=ncount+1
                namfree(ncount)=vs%nampar(i)
                pn(i)=a(ncount)
                ch(i)=pn(i)-vs%pv(i)
                vs%spv(i)=sa(ncount)
                vs%pv(i)=a(ncount)
             end if
          end do
       end if

       call Info_LSQ_Output(chi2,FL,nobs,x,y,yc,w,Ipr,c,vs)

       return
    End Subroutine Marquardt_Fit_v1

    !!--++
    !!--++ Subroutine Marquardt_Fit_v2(Model_Functn, d, c, vs, Ipr, Chi2, scroll_lines)
    !!--++    Type(LSQ_Data_Type),                     intent(in out) :: d             !Vector of x-values
    !!--++    type(LSQ_conditions_type),               intent(in out) :: c             !Conditions for the algorithm
    !!--++    type(LSQ_State_Vector_type),             intent(in out) :: vs            !State vector for the model calculation
    !!--++    integer                    ,             intent(in)     :: Ipr           !Logical unit for writing
    !!--++    real(kind=cp),                           intent(out)    :: chi2          !Reduced Chi-2
    !!--++    character(len=*),dimension(:), optional, intent(out)    :: scroll_lines  !If present, part of the output is stored
    !!--++                                                                             !in this text for treatment in the calling program
    !!--++
    !!--++    Model_functn                                                             !Name of the subroutine calculating yc(i) for point x(i)
    !!--++    Interface                                                                !Interface for the Model_Functn subroutine
    !!--++       Subroutine Model_Functn(iv,Xv,ycalc,Vsa,calder)
    !!--++          use CFML_GlobalDeps, only: cp
    !!--++          integer,                    intent(in)     :: iv     !Number of the component: "i" in x(i)
    !!--++          real(kind=cp),              intent(in)     :: xv     !Value of x(i)
    !!--++          real(kind=cp),              intent(out)    :: ycalc  !Value of yc at point x(i) => ycalc=yc(i)
    !!--++          Type(LSQ_State_Vector_type),intent(in out) :: Vsa    !Parameters, codes, and derivatives
    !!--++          logical,optional,           intent(in)     :: calder !If present, derivatives, stored in Vsa%dpv, are calculated
    !!--++       End Subroutine Model_Functn
    !!--++    End Interface
    !!--++
    !!--++
    !!--++    Subroutine for applying the Levenberg-Marquardt method for Least-Squares.
    !!--++    The user must provide a model function according to the interface above.
    !!--++    The model function should use at least some of the public variables of the
    !!--++    present module in order to set the derivatives with respect to the model
    !!--++    parameters. Examples of using this module are given in the program
    !!--++    templates CW_fit and TOF_fit.
    !!--++
    !!--..    INFORMATION
    !!--..        Author: Juan Rodriguez-Carvajal (based in text descriptions of the literature)
    !!--..                Translated, extracted and modified from the Fortran 77 code: XRFIT (JRC 1986)
    !!--..        Module implementing the Marquardt method for non-linear least-squares.
    !!--..        For using this module, the user must provide:
    !!--..            1: Number of cycles (c%icyc), type of weighting scheme (c%iw), number of
    !!--..               model parameters (vs%np), constraint conditions (c%constr, c%percent)
    !!--..            2: A character(len=15) name for all possible parameters of the
    !!--..               model stored in the array vs%nampar(:).
    !!--..            3: A set of initial value for all parameters stored in vs%pv(:)
    !!--..            4: A set of flags values to refine or fix the parameters vs%code(:)
    !!--..            5: The model function (subroutine model_functn) must be provided by the user
    !!--..               according to the interface described below. The actual name of the model function
    !!--..               is arbitrary and it is passed to the only public procedure "marquardt_fit" as a
    !!--..               dummy argument.
    !!--..
    !!--..       The values of all possible refinable parameters are stored in the array vs%pv(:).
    !!--..       The derivatives must be calculated within model_functn, by using the array vs%dpv(:)
    !!--..       The actual refined parameters are selected from the larger vs%pv(:) array by
    !!--..       the integer array of flags: vs%code(:). A value vs%code(j)=1 means that the j-th
    !!--..       parameter is to be varied. A value vs%code(k)=0 means that the k-th parameter is
    !!--..
    !!--..       Subroutine Model_Functn(iv,Xv,ycalc,Vsa,calder)
    !!--..       It is recommended that Model_Functn be stored in a Module. The integer iv (counter
    !!--..       of the loop for all observations for which the subroutine is invoked) is passed
    !!--..       because in many cases part of the calculations and derivatives may be calculated only
    !!--..       for iv=1 and stored in local variables that should have the SAVE attribute or be private
    !!--..       module variables accessible by host association. The only output values of the subroutine
    !!--..       are ycalc and Vsa%dpv(:) that vary for each "iv" point.
    !!--..       In this version of the algorithm the derivatives with respect to the free parameters
    !!--..       in Model_Functn should be calculated before exiting by using the following loop:
    !!--..
    !!--..       vs%dpv=0.0          !Should be always be initialized for each point
    !!--..       if(present(calder)) then
    !!--..         Do i=1,vs%np
    !!--..           if(vs%code(i) == 0) cycle
    !!--..           vs%dpv(i) = derivative_wrt(i) !symbolic form for calculating the derivative w.r.t. parameter "i"
    !!--..         end do                          !at the current point xv (observation "iv")
    !!--..       end if
    !!--++
    !!--++ Update: August - 2009
    !!
    Subroutine Marquardt_Fit_v2(Model_Functn,d,c,vs,Ipr,Chi2,scroll_lines)
       !---- Arguments ----!
       Type(LSQ_Data_Type),           intent(in out)         :: d
       type(LSQ_conditions_type),     intent(in out)         :: c
       type(LSQ_State_Vector_type),   intent(in out)         :: vs
       integer,                       intent(in)             :: Ipr
       real(kind=cp),                 intent(   out)         :: chi2
       character(len=*),dimension(:), intent(out), optional  :: scroll_lines

       Interface
        Subroutine Model_Functn(iv,xv,ycalc,Vsa,calder)
           use CFML_GlobalDeps, only: cp
           use CFML_LSQ_TypeDef,  only: LSQ_State_Vector_type
           integer,                     intent(in)     :: iv
           real(kind=cp),               intent(in)     :: xv
           real(kind=cp),               intent(out)    :: ycalc
           Type(LSQ_State_Vector_type), intent(in out) :: Vsa
           logical,optional,            intent(in)     :: calder
        End Subroutine Model_Functn
       End Interface

       !---- local variables ----!
       logical                                :: fixed, write_cyc
       character(len=180)                     :: line
       integer                                :: ifail
       integer                                :: i,ncount,li,ntex,nt
       real(kind=cp)                          :: fl,chi1
       real(kind=cp), dimension(Max_Free_Par) :: a,sa

       fixed = .false.
       write_cyc=.true.
       if (c%icyc > 50) write_cyc=.false.

       !---- Beginning the iterations ----!
       fl=0.001
       chi1=1.0e+30
       c%reached =.false.
       write(unit=ipr,fmt="(/,a)")   " -------------------------------------------------"
       write(unit=ipr,fmt="(a,i3,a)")" => Marquardt Least Squares Fitting for ",c%icyc," cycles"
       write(unit=ipr,fmt="(a,/)")   " -------------------------------------------------"

       write(unit=ipr,fmt="(a,i6)")        " => Number of observations   : ",d%nobs
       ncount=0
       do i=1,vs%np
          if (vs%code(i) == 0) cycle
          ncount=ncount+1
          namfree(ncount)=vs%nampar(i)
       end do
       c%npvar=ncount
       write(unit=ipr,fmt="(a,i6)")        " => Number of free parameters: ",ncount
       write(unit=ipr,fmt="(2(a,f12.4),a)")" => Range of variable X      : (",d%x(1),",",d%x(d%nobs),")"

       ntex=0
       if(present(scroll_lines)) scroll_lines=" "  !Clear the stored text
       do li=1,c%icyc                      !Loop over icyc cycles

          call curfit_v2(Model_Functn,d,c,vs,fl,chi2,ifail,nt)

          if (ifail /= 0) then
             line= " "
             write(unit=line,fmt="(a,i3,a)") " => IFAIL /= 0 : ",ifail, " "//trim(ERR_Lsq_Mess)
             if(present(scroll_lines)) then
                ntex=ntex+1
                scroll_lines(ntex)=line
              else
                 if(ifail /= 1) write(unit=*,fmt="(a)") trim(line)
            end if
             if (ABS(chi2-chi1)*100/Chi2 <= 0.001 .and. ifail==1) then
                c%reached=.true.
                line= " "
                write(unit=line,fmt="(a)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.00001 * Chi2(old)"
                if(present(scroll_lines)) then
                    ntex=ntex+1
                    scroll_lines(ntex)=line
                else
                    write(unit=*,fmt="(a)") trim(line)
                end if
                write(unit=Ipr,fmt="(/,a,/)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.00001 * Chi2(old)"
                exit
             end if

             do i=1,d%nobs            !Calculation of the function for all points with the previous parameters
                call model_functn(i,d%x(i),d%yc(i),vs)
             end do
             chi2=fchisq(d%nobs-c%npvar,d%nobs,d%y,d%sw,d%yc)
             line= " "
             write(unit=line,fmt="(a,i5,2(a,f14.6),a,i3)") " => Marquardt convergence not reached! Cycle No.:",&
                                            li,"  Last Chi2 =",chi2, "  Marquardt-Lambda: ",fl, "  Ntrials: ",nt
             if(present(scroll_lines)) then
                ntex=ntex+1
                scroll_lines(ntex)=line
              else
                write(unit=*,fmt="(a)") trim(line)
             end if

             write(unit=ipr,fmt="(/,a)") trim(line)
             write(unit=ipr,fmt="(a)" ) "    If Chi2 is reasonably low, you may have reached the minimum"
             write(unit=ipr,fmt="(a)" ) "    within the precision of the machine.                       "
             write(unit=ipr,fmt="(a)" ) "    Otherwise, look at the output file to search for the reason "
             write(unit=ipr,fmt="(a,/)")"    of anomalous behaviour and try to start from new parameters"

             write(unit=ipr,fmt="(/,a,/)")" =>  Diagonal elements of LSQ-matrix before inversion"
             do i=1,c%npvar
                write(unit=ipr,fmt="(a30,i3,a,i3,a,f16.6)") "    "//trim(namfree(i))//"  (",i,",",i,")", curv_mat(i,i)
             end do

             exit
          end if  !ifail/=0

          if (abs(chi2-chi1)*100.0/chi2 >= 0.0001 ) then
             c%reached=.false.
             chi1=chi2
              if (c%constr)  then
                !Extract the vectors a and sa from pn
                ncount=0
                do i=1,vs%np
                   if (vs%code(i) == 0) cycle
                   ncount=ncount+1
                   a(ncount)=pn(i)
                   sa(ncount)=vs%spv(i)
                end do
                call box_constraints(a,sa,fixed,c,vs)
                if (fixed) then
                   write(unit=ipr,fmt="(/,a,/,a,i2,a)") " => Some parameters have been restored to their initial value ", &
                                                        "    because of local divergence (change > ",nint(c%percent),"%)"
                   write(unit=ipr,fmt="(a,/)")          "    Look at the above list for parameters having a '*'"
                end if  !fixed
             end if

             line= " "
             write(unit=line,fmt="(a,i5,2(a,f14.6),a,i3)") &
                " => Cycle No.:",li,"  Chi2 =",chi2, "  Marquardt-Lambda: ",fl, "  Ntrials: ",nt
             if(present(scroll_lines)) then
                 ntex=ntex+1
                 scroll_lines(ntex)=line
             else
                 write(unit=*,fmt="(a)") trim(line)
             end if

             if (write_cyc) call output_cyc(li,Ipr,chi2,vs)
             vs%pv(1:vs%np)=pn(1:vs%np)  !Change is done here
             cycle
          else
             line= " "
             write(unit=line,fmt="(a)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.000001 * Chi2(old)"
             if(present(scroll_lines)) then
                 ntex=ntex+1
                 scroll_lines(ntex)=line
             else
                 write(unit=*,fmt="(a)") trim(line)
             end if
             write(unit=Ipr,fmt="(/,a,/)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.000001 * Chi2(old)"
             c%reached=.true.
             exit
          end if
       end do  !li=1,c%icyc

       if (c%reached) then
          if (fl <= 10.0) then
             fl=0.0
          else
             fl=0.5e-20
          end if
          call curfit_v2(Model_Functn,d,c,vs,fl,chi2,ifail,nt)
       end if


       call Info_LSQ_Output(chi2,FL,d%nobs,d%x,d%y,d%yc,d%sw,Ipr,c,vs)

       return
    End Subroutine Marquardt_Fit_v2

    !!--++
    !!--++  Subroutine Output_Cyc(Ic,Lun,Chi2,vs)
    !!--++   integer,                    intent(in) :: ic  !cycle number
    !!--++   integer,                    intent(in) :: lun !logical number of the output file
    !!--++   real(kind=cp),              intent(in) :: chi2
    !!--++   type(LSQ_State_Vector_type),intent(in) :: vs
    !!--++
    !!--++  (PRIVATE)
    !!--++  Subroutine for output information on each cycle
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Output_Cyc(Ic,Lun,Chi2,vs)
       !---- Arguments ----!
       integer,                    intent(in) :: ic  !cycle number
       integer,                    intent(in) :: lun !logical number of the output file
       real(kind=cp),              intent(in) :: chi2
       type(LSQ_State_Vector_type),intent(in) :: vs

       !---- Local variables ----!
       integer       :: i,j
       real(kind=cp) :: rat

       !---- Writing during cycles ----!
       write(unit=lun,fmt="(/,/,a,i5,a,f14.6)")" => Cycle No.:",ic,"  Chi2 =",chi2
       write(unit=lun,fmt="(/,/,a,/)") &
       "              Name-Par       No.      Old-Value          Change        New-Value         Sigma        Change/Sigma"
       j=0
       do i=1,vs%np
          if (vs%code(i)/=0) then
             j=j+1
             if (vs%spv(i) > 1.0e-36) then
                rat=ch(i)/vs%spv(i)
             else
                rat=0.0
             end if
             write(unit=lun,fmt="(a25,i6,a,5f16.6)") " "//trim(namfree(j)),i," ",vs%pv(i),ch(i),pn(i),vs%spv(i),rat
          end if
       end do

       return
    End Subroutine Output_Cyc

    !!--++
    !!--++  Subroutine qrfac(m, n, a, pivot, ipvt, rdiag, acnorm)
    !!--++     Integer,                       Intent(In)      :: m
    !!--++     Integer,                       Intent(In)      :: n
    !!--++     Real (Kind=cp),Dimension(:,:), Intent(In Out)  :: a
    !!--++     Logical,                       Intent(In)      :: pivot
    !!--++     Integer,Dimension(:),          Intent(Out)     :: ipvt
    !!--++     Real (Kind=cp),Dimension(:),   Intent(Out)     :: rdiag
    !!--++     Real (Kind=cp),Dimension(:),   Intent(Out)     :: acnorm
    !!--++
    !!--..  Original documentation:
    !!--..
    !!--..  This subroutine uses Householder transformations with column pivoting
    !!--..  (optional) to compute a qr factorization of the m by n matrix a. That is,
    !!--..  qrfac determines an orthogonal matrix q, a permutation matrix p, and an
    !!--..  upper trapezoidal matrix r with diagonal elements of nonincreasing magnitude,
    !!--..  such that a*p = q*r.  The householder transformation for column k,
    !!--..  k = 1,2,...,min(m,n), is of the form
    !!--..
    !!--..                        t
    !!--..        i - (1/u(k))*u*u
    !!--..
    !!--..  where u has zeros in the first k-1 positions.  The form of this transformation
    !!--..  and the method of pivoting first appeared in the corresponding linpack i
    !!--..  subroutine.
    !!--..
    !!--..  The subroutine statement is
    !!--..
    !!--..      subroutine qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa)
    !!--..
    !!--..  where
    !!--..
    !!--..    m     is a positive integer input variable set to the number of rows of a.
    !!--..
    !!--..    n     is a positive integer input variable set to the number of columns of a.
    !!--..
    !!--..    a     is an m by n array.  On input a contains the matrix for which the qr
    !!--..          factorization is to be computed.  On output the strict upper trapezoidal
    !!--..          part of a contains the strict upper trapezoidal part of r, and the lower
    !!--..          trapezoidal part of a contains a factored form of q (the non-trivial
    !!--..          elements of the u vectors described above).
    !!--..
    !!--..    lda   is a positive integer input variable not less than m which specifies
    !!--..          the leading dimension of the array a.
    !!--..
    !!--..    pivot is a logical input variable.  If pivot is set true, then column pivoting
    !!--..          is enforced.  If pivot is set false, then no column pivoting is done.
    !!--..
    !!--..    ipvt  is an integer output array of length lipvt.  ipvt defines the permutationi
    !!--..          matrix p such that a*p = q*r. Column j of p is column ipvt(j) of the i
    !!--..          identity matrix. If pivot is false, ipvt is not referenced.
    !!--..
    !!--..    lipvt is a positive integer input variable.  If pivot is false, then lipvt may
    !!--..          be as small as 1.  If pivot is true, then lipvt must be at least n.
    !!--..
    !!--..    rdiag is an output array of length n which contains the diagonal elements of
    !!--..          r.
    !!--..
    !!--..   acnorm is an output array of length n which contains the norms of the
    !!--..          corresponding columns of the input matrix a. If this information is not
    !!--..          needed, then acnorm can coincide with rdiag.
    !!--..
    !!--..    wa    is a work array of length n.  If pivot is false, then wa can coincide i
    !!--..          with rdiag.
    !!--..
    !!--..  Argonne National Laboratory. minpack project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!--++
    !!--++
    !!--++ Update: August - 2009
    !!
    Subroutine qrfac(m, n, a, pivot, ipvt, rdiag, acnorm)
       !---- Arguments ----!
       Integer,                       Intent(In)      :: m
       Integer,                       Intent(In)      :: n
       Real (Kind=cp),Dimension(:,:), Intent(In Out)  :: a
       Logical,                       Intent(In)      :: pivot
       Integer,Dimension(:),          Intent(Out)     :: ipvt
       Real (Kind=cp),Dimension(:),   Intent(Out)     :: rdiag
       Real (Kind=cp),Dimension(:),   Intent(Out)     :: acnorm

       !--- Local Variables ---!
       Integer                      :: i, j, jp1, k, kmax, minmn
       Real (Kind=cp)               :: ajnorm, epsmch, suma, temp
       Real (Kind=cp), Dimension(n) :: wa
       Real (Kind=cp), Parameter    :: one = 1.0_CP, p05 = 0.05_CP, zero = 0.0_CP

       !     epsmch is the machine precision.
       epsmch = Epsilon(zero)

       !     compute the initial column norms and initialize several arrays.
       Do j = 1, n
          acnorm(j) = enorm(m,a(1:,j))
          rdiag(j) = acnorm(j)
          wa(j) = rdiag(j)
          If (pivot) ipvt(j) = j
       End Do

       !     Reduce a to r with Householder transformations.
       minmn = Min(m,n)
       Do j = 1, minmn
          If (pivot) Then

             !Bring the column of largest norm into the pivot position.
             kmax = j
             Do k = j, n
                If (rdiag(k) > rdiag(kmax)) kmax = k
             End Do

             If (kmax /= j) Then
                Do i = 1, m
                   temp = a(i,j)
                   a(i,j) = a(i,kmax)
                   a(i,kmax) = temp
                End Do
                rdiag(kmax) = rdiag(j)
                wa(kmax) = wa(j)
                k = ipvt(j)
                ipvt(j) = ipvt(kmax)
                ipvt(kmax) = k
             End If
          End If

          !     Compute the Householder transformation to reduce the
          !     j-th column of a to a multiple of the j-th unit vector.
          ajnorm = enorm(m-j+1, a(j:,j))
          If (ajnorm == zero) Cycle
          If (a(j,j) < zero) ajnorm = -ajnorm
          a(j:m,j) = a(j:m,j)/ajnorm
          a(j,j) = a(j,j) + one

          !     Apply the transformation to the remaining columns and update the norms.
          jp1 = j + 1
          Do k = jp1, n
             suma = Dot_Product( a(j:m,j), a(j:m,k) )
             temp = suma/a(j,j)
             a(j:m,k) = a(j:m,k) - temp*a(j:m,j)
             If (.NOT. pivot .OR. rdiag(k) == zero) Cycle
             temp = a(j,k)/rdiag(k)
             rdiag(k) = rdiag(k)*Sqrt(Max(zero, one-temp**2))
             If (p05*(rdiag(k)/wa(k))**2 > epsmch) Cycle
             rdiag(k) = enorm(m-j, a(jp1:,k))
             wa(k) = rdiag(k)
          End Do
          rdiag(j) = -ajnorm
       End Do

       Return
    End Subroutine qrfac

    !!--++
    !!--++ Subroutine qrsolv(n, r, ipvt, diag, qtb, x, sdiag)
    !!--++    Integer,                        Intent(In)      :: n
    !!--++    Real (Kind=cp), dimension(:,:), Intent(In Out)  :: r
    !!--++    Integer, dimension(:),          Intent(In)      :: ipvt
    !!--++    Real (Kind=cp), dimension(:),   Intent(In)      :: diag
    !!--++    Real (Kind=cp), dimension(:),   Intent(In)      :: qtb
    !!--++    Real (Kind=cp), dimension(:),   Intent(Out)     :: x
    !!--++    Real (Kind=cp), dimension(:),   Intent(Out)     :: sdiag
    !!--++
    !!--..    Original documentation:
    !!--..
    !!--..    Given an m by n matrix a, an n by n diagonal matrix d, and an m-vector b,
    !!--..    the problem is to determine an x which solves the system
    !!--..
    !!--..        a*x = b ,     d*x = 0 ,
    !!--..
    !!--..    in the least squares sense.
    !!--..
    !!--..    This subroutine completes the solution of the problem if it is provided
    !!--..    with the necessary information from the qr factorization, with column
    !!--..    pivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
    !!--..    q has orthogonal columns, and r is an upper triangular matrix with diagonal
    !!--..    elements of nonincreasing magnitude, then qrsolv expects the full upper
    !!--..    triangle of r, the permutation matrix p, and the first n components of
    !!--..    (q transpose)*b.  The system a*x = b, d*x = 0, is then equivalent to
    !!--..
    !!--..               t       t
    !!--..        r*z = q *b ,  p *d*p*z = 0 ,
    !!--..
    !!--..    where x = p*z. if this system does not have full rank,
    !!--..    then a least squares solution is obtained.  On output qrsolv
    !!--..    also provides an upper triangular matrix s such that
    !!--..
    !!--..         t   t               t
    !!--..        p *(a *a + d*d)*p = s *s .
    !!--..
    !!--..    s is computed within qrsolv and may be of separate interest.
    !!--..
    !!--..    The subroutine statement is
    !!--..
    !!--..        subroutine qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)
    !!--..
    !!--..    where
    !!--..
    !!--..    n     is a positive integer input variable set to the order of r.
    !!--..
    !!--..    r     is an n by n array.  On input the full upper triangle must contain
    !!--..          the full upper triangle of the matrix r.
    !!--..          On output the full upper triangle is unaltered, and the strict lower
    !!--..          triangle contains the strict upper triangle (transposed) of the
    !!--..          upper triangular matrix s.
    !!--..
    !!--..    ldr   is a positive integer input variable not less than n which specifies
    !!--..          the leading dimension of the array r.
    !!--..
    !!--..    ipvt   is an integer input array of length n which defines the permutation
    !!--..           matrix p such that a*p = q*r.  Column j of p is column ipvt(j) of
    !!--..           the identity matrix.
    !!--..
    !!--..    diag   is an input array of length n which must contain the diagonal elements
    !!--..           of the matrix d.
    !!--..
    !!--..    qtb    is an input array of length n which must contain the first n elements
    !!--..           of the vector (q transpose)*b.
    !!--..
    !!--..    x      is an output array of length n which contains the least squares
    !!--..           solution of the system a*x = b, d*x = 0.
    !!--..
    !!--..    sdiag  is an output array of length n which contains the diagonal elements
    !!--..           of the upper triangular matrix s.
    !!--..
    !!--..    wa     is a work array of length n.
    !!--..
    !!--..  Argonne National Laboratory. minpack project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!--++
    !!--++ Update: August - 2009
    !!
    Subroutine qrsolv(n, r, ipvt, diag, qtb, x, sdiag)
       !---- Arguments ----!
       Integer,                        Intent(In)      :: n
       Real (Kind=cp), dimension(:,:), Intent(In Out)  :: r
       Integer, dimension(:),          Intent(In)      :: ipvt
       Real (Kind=cp), dimension(:),   Intent(In)      :: diag
       Real (Kind=cp), dimension(:),   Intent(In)      :: qtb
       Real (Kind=cp), dimension(:),   Intent(Out)     :: x
       Real (Kind=cp), dimension(:),   Intent(Out)     :: sdiag

       !--- Local Variables ---!
       Integer                      :: i, j, k, kp1, l, nsing
       Real (Kind=cp)               :: cosi, cotan, qtbpj, sini, suma, tani, temp
       Real (Kind=cp), dimension(n) ::  wa
       Real (Kind=cp), Parameter :: p5 = 0.5_CP, p25 = 0.25_CP, zero = 0.0_CP

       !     Copy r and (q transpose)*b to preserve input and initialize s.
       !     In particular, save the diagonal elements of r in x.

       Do j = 1, n
          r(j:n,j) = r(j,j:n)
          x(j) = r(j,j)
          wa(j) = qtb(j)
       End Do

       !     Eliminate the diagonal matrix d using a givens rotation.
       Do j = 1, n

          !        Prepare the row of d to be eliminated, locating the
          !        diagonal element using p from the qr factorization.
          l = ipvt(j)
          If (diag(l) == zero) Cycle
          sdiag(j:n) = zero
          sdiag(j) = diag(l)

          !     The transformations to eliminate the row of d modify only a single
          !     element of (q transpose)*b beyond the first n, which is initially zero.
          qtbpj = zero
          Do k = j, n

             !        Determine a givens rotation which eliminates the
             !        appropriate element in the current row of d.
             If (sdiag(k) == zero) Cycle
             If (Abs(r(k,k)) < Abs(sdiag(k))) Then
                cotan = r(k,k)/sdiag(k)
                sini = p5/Sqrt(p25 + p25*cotan**2)
                cosi = sini*cotan
             Else
                tani = sdiag(k)/r(k,k)
                cosi = p5/Sqrt(p25 + p25*tani**2)
                sini = cosi*tani
             End If

             !        Compute the modified diagonal element of r and
             !        the modified element of ((q transpose)*b,0).
             r(k,k) = cosi*r(k,k) + sini*sdiag(k)
             temp = cosi*wa(k) + sini*qtbpj
             qtbpj = -sini*wa(k) + cosi*qtbpj
             wa(k) = temp

             !        Accumulate the tranformation in the row of s.
             kp1 = k + 1
             Do i = kp1, n
                temp = cosi*r(i,k) + sini*sdiag(i)
                sdiag(i) = -sini*r(i,k) + cosi*sdiag(i)
                r(i,k) = temp
             End Do
          End Do

          !     Store the diagonal element of s and restore
          !     the corresponding diagonal element of r.
          sdiag(j) = r(j,j)
          r(j,j) = x(j)
       End Do

       !     Solve the triangular system for z.  If the system is singular,
       !     then obtain a least squares solution.
       nsing = n
       Do j = 1, n
          If (sdiag(j) == zero .AND. nsing == n) nsing = j - 1
          If (nsing < n) wa(j) = zero
       End Do

       Do k = 1, nsing
          j = nsing - k + 1
          suma = Dot_Product( r(j+1:nsing,j), wa(j+1:nsing) )
          wa(j) = (wa(j) - suma)/sdiag(j)
       End Do

       !     Permute the components of z back to components of x.
       Do j = 1, n
          l = ipvt(j)
          x(l) = wa(j)
       End Do

       Return
    End Subroutine qrsolv

 End Module CFML_Optimization_LSQ
