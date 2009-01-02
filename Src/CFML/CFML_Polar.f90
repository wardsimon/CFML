!!----
!!---- Copyleft(C) 1999 - 2009,              Version: 4.0
!!---- Juan Rodriguez-Carvajal & Marc Janoschek
!!----
!!---- MODULE: CFML_Polarimetry
!!----   INFO: Subroutines and Functions to calculate the polarisation tensor
!!----   as it will be measured. It uses matrices defined in CFML_Crystal_Metrics in
!!----   order to calculate the polar tensor with respect to the coordinate
!!----   frame defined in the Blume equations (Phys. Rev. Vol. 130 p.1670-1676,
!!----   1963, see also the definitions below in magn_Inter_Vec_PF). As input
!!----   the nuclear structure factor, the magnetic interaction vector with
!!----   respect to the crystal frame and the matrices defined in CFML_Crystal_Metrics
!!----   for the crystal frame are needed.
!!----
!!---- HISTORY
!!----    Update: April - 2008
!!----            April - 2005: Created by MJ and revised by JRC
!!----            December - 2006: Added function Write_Polar_line for more convenient
!!----                             output of matrices of many reflections in one file
!!----
!!---- DEPENDENCIES
!!--++    use CFML_Crystal_Metrics, only: Set_Crystal_Cell, Crystal_Cell_type, Cart_Vector
!!--++    Use CFML_Constants,       only: cp,tpi
!!--++    USE CFML_Math_3D,         only: Cross_Product
!!----
!!---- VARIABLES
!!----    POLAR_INFO_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!--++       IM_NM_Y                   [Private]
!!--++       IM_NM_Z                   [Private]
!!--++       MAG_Y                     [Private]
!!--++       MAG_Z                     [Private]
!!--++       MAGN_INTER_VEC_PF         [Private]
!!--++       MM                        [Private]
!!--++       NUC_CONTR                 [Private]
!!--++       REAL_NM_Y                 [Private]
!!--++       REAL_NM_Z                 [Private]
!!--++       TCHIRAL                   [Private]
!!----
!!----    Subroutines:
!!----       SET_POLAR_INFO
!!----       WRITE_POLAR_INFO
!!----       WRITE_POLAR_LINE
!!----
!!
 Module CFML_Polarimetry
    !---- Used External Modules ----!
    Use CFML_Constants,       only: cp, tpi
    Use CFML_Crystal_Metrics, only: Set_Crystal_Cell, Crystal_Cell_type, Cart_Vector
    Use CFML_Math_3D,         only: Cross_Product

    !---- Variables ----!
    implicit none

    private

    !---- List of public overloaded operators ----!

    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public  :: Set_Polar_Info, Write_Polar_Info, Write_Polar_line

    !---- List of private Operators ----!

    !---- List of private functions ----!
    private :: Magn_Inter_Vec_PF, Nuc_Contr, Mag_Y, Mag_Z, Real_Nm_Y, Real_Nm_Z, &
               Im_Nm_Y, Im_Nm_Z, Tchiral, Mm

    !---- List of private subroutines ----!


    !---- Definitions ----!

    !!----
    !!---- TYPE :: POLAR_INFO_TPYE
    !!--..
    !!---- Type, public :: Polar_Info_type
    !!----     real(kind=cp), dimension (3)    :: H     ! Scattering vector in hkl
    !!----     real(kind=cp), dimension (3)    :: SPV   ! Second vector in Scattering plane apart of scattering vector to define plane
    !!----     type(crystal_cell_type)         :: Cell  ! Unit Cell of Crystal
    !!----     real(kind=cp)                   :: P     ! magnitude of initial polarisation vector
    !!----     complex, dimension (3)          :: MIV   ! magnetic interaction vector
    !!----     complex                         :: NSF   ! nuclear structure factor
    !!----     real(kind=cp)                   :: NC    ! nuclear scattering contribution
    !!----     real(kind=cp)                   :: MY    ! magnetic contribution along y
    !!----     real(kind=cp)                   :: MZ    ! magnetic contribution along z
    !!----     real(kind=cp)                   :: RY    ! real part of nuclear magnetic interference term along y
    !!----     real(kind=cp)                   :: RZ    ! real part of nuclear magnetic interference term along z
    !!----     real(kind=cp)                   :: IY    ! imaginary part of nuclear magnetic interference term along y
    !!----     real(kind=cp)                   :: IZ    ! imaginary part of nuclear magnetic interference term along y
    !!----     real(kind=cp)                   :: TC    ! chiral contribution
    !!----     real(kind=cp)                   :: MM    ! magnetic-magnetic interference term
    !!----     real(kind=cp), dimension (3)    :: CS    ! the three different elastic cross-sections depending on the direction of the initial polar vector
    !!----     real(kind=cp), dimension (3,3)  :: Pij   ! the polarisation tensor
    !!---- End Type Polar_Info_type
    !!----
    !!---- Update: April 2008
    !!
    Type, public :: Polar_Info_type
       real(kind=cp), DIMENSION (3)     :: H
       real(kind=cp), DIMENSION (3)     :: SPV
       TYPE(Crystal_Cell_Type)          :: Cell
       real(kind=cp)                    :: P
       COMPLEX, DIMENSION (3)           :: MIV
       COMPLEX                          :: NSF
       real(kind=cp)                    :: NC
       real(kind=cp)                    :: MY
       real(kind=cp)                    :: MZ
       real(kind=cp)                    :: RY
       real(kind=cp)                    :: RZ
       real(kind=cp)                    :: IY
       real(kind=cp)                    :: IZ
       real(kind=cp)                    :: TC
       real(kind=cp)                    :: MM
       real(kind=cp), DIMENSION (3)     :: CS
       real(kind=cp), DIMENSION (3,3)   :: Pij
    End Type Polar_Info_type

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!
    
    !!--++
    !!--++ Function Im_Nm_Y(Nsf, Miv_Pf) Result(I_Nm_Y)
    !!--++    Complex,               intent( in)  :: NSF      !  In  -> Nuclear Structure Factor
    !!--++    Complex, dimension(3), intent( in)  :: MIV_PF   !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                       :: I_NM_Y   !  Out -> Imaginary part of nuclear magnetic interference contribution along Y
    !!--++
    !!--++    (Private)
    !!--++    Calculates the imaginary part of the nuclear magnetic interference contribution along Y
    !!--++    to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Im_Nm_Y(Nsf, Miv_Pf) Result(I_Nm_Y)
       !---- Argument ----!
       Complex,               intent( in)  :: NSF
       Complex, dimension(3), intent( in)  :: MIV_PF
       real(kind=cp)                       :: I_NM_Y
       
       I_NM_Y = AIMAG(NSF * CONJG(MIV_PF(2)) - CONJG(NSF) * MIV_PF(2))
       
       return
    End Function  Im_Nm_Y
    
    !!--++
    !!--++ Real Function Im_Nm_Z(Nsf, Miv_Pf) Result(I_Nm_Z)
    !!--++    Complex,               intent(in) :: NSF     !  In -> Nuclear Structure Factor
    !!--++    Complex, dimension(3), intent(in) :: MIV_PF  !  In -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                     :: I_NM_Z  !  Out-> Imaginary part of nuclear magnetic interference contribution along Z
    !!--++
    !!--++    (Private)
    !!--++    Calculates the imaginary part of the nuclear magnetic interference contribution along Z
    !!--++    to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Im_Nm_Z(Nsf, Miv_Pf) Result(I_Nm_Z)
       !---- Argument ----!
       Complex,               intent(in) :: NSF
       Complex, dimension(3), intent(in) :: MIV_PF
       real(kind=cp)                     :: I_NM_Z
       
       !---- Local variables ----!

       I_NM_Z = AIMAG(NSF * CONJG(MIV_PF(3)) - CONJG(NSF) * MIV_PF(3))
       
       return
    End Function  Im_Nm_Z
    
    !!--++
    !!--++ Function Mag_Y(Miv_Pf) Result(My)
    !!--++    Complex, dimension(3), intent( in) :: MIV_PF !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                      :: MY     !  Out -> Magnetic contribution along Y
    !!--++
    !!--++    (Private)
    !!--++    Calculates the magnetic contribution along Y to scattering in the polarisation
    !!--++    coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Mag_Y(Miv_Pf) Result(My)
       !---- Argument ----!
       Complex, dimension(3), intent( in)  :: MIV_PF
       real(kind=cp)                       :: MY
       
       !---- Local variables ----!

       MY = MIV_PF(2) * CONJG(MIV_PF(2))
       
       return
    End Function  Mag_Y
    
    !!--++
    !!--++ Function Mag_Z(Miv_Pf) Result(Mz)
    !!--++    Complex, dimension(3), intent( in):: MIV_PF  !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                     :: MZ      !  Out -> Magnetic contribution along Z
    !!--++
    !!--++    (Private)
    !!--++    Calculates the magnetic contribution along Z to scattering in the polarisation
    !!--++    coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Mag_Z(Miv_Pf) Result(Mz)
       !---- Argument ----!
       Complex, dimension(3), intent( in)       :: MIV_PF
       real(kind=cp)                            :: MZ
       
       !---- Local variables ----!

       MZ = MIV_PF(3) * CONJG(MIV_PF(3))
       
       return
    End Function  Mag_Z

    !!--++
    !!--++ Function Magn_Inter_Vec_Pf(Miv,H,Spv, Cell) Result(Miv_Pf)
    !!--++    Complex, dimension(3), intent( in)       :: MIV            !  In -> Magnetic Interaction Vector
    !!--++    real(kind=cp), dimension(3), intent( in) :: H              !  In -> Scattering Vector in hkl
    !!--++    real(kind=cp), dimension(3), intent( in) :: SPV            !  In -> Second Scattering plane vector in hkl
    !!--++    Type (Crystal_Cell_Type),  intent(in)    :: Cell           !  In -> Cell variable which holds transformation matrices
    !!--++    Complex, dimension(3)                    :: MIV_PF         !  Out -> Magnetic Interaction Vector in polarisation coordinate frame
    !!--++
    !!--++    (Private)
    !!--++    Calculates the magnetic interaction vector in the polarisation coordinate frame according to the Blume equations
    !!--++    and therefore depending on the scattering vector.
    !!--++
    !!--++    Polarisation coordinate frame according to Blume
    !!--++    X  || scattering vector Q     (where Q is the scattering Vector in cartesian real space coordinates, it will be calculated from H and matrices in Cell)
    !!--++    Y _|_ scattering vector Q in scattering plane
    !!--++    Z _|_ scattering vector Q out of scattering plane (ATTENTION: This is choice not non-ambiguous, there are always two possible choices
    !!--++                                                       for a right handed coordinate frame which will fullfils this condition!!!)
    !!--++
    !!--++                           Y
    !!--++                          /|\
    !!--++                           |
    !!--++                   Q       | Z
    !!--++            ____________ _\o_____\ X
    !!--++            \             /      /
    !!--++             \           /
    !!--++              \         /
    !!--++               \       /
    !!--++                \     /  K_f
    !!--++             K_i \   /
    !!--++                 _\|/_
    !!--++
    !!--++    Therefore the right handed coordinate frame will be explicitly chosen like this:
    !!--++    X := Q/|Q|               where Q is the scattering Vector in cartesian real space coordinates
    !!--++    Z := (Q x SV)/|(Q x SV)| where SV is a second vector in the scattering plane in cartesian real space coordinates
    !!--++    Y := (Z x X)
    !!--++
    !!--++    ATTENTION: Be aware that the choice of SV with respect to Q will decide which of the two possible right handed coordinates fullfilling
    !!--++               the conditions above wille be used!!!
    !!--++
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Magn_Inter_Vec_Pf(Miv,H,Spv, Cell) Result(Miv_Pf)
       !---- Argument ----!
       Complex, dimension(3),       intent(in) :: MIV
       real(kind=cp), dimension(3), intent(in) :: H
       real(kind=cp), dimension(3), intent(in) :: SPV
       Type (Crystal_Cell_Type),    intent(in) :: Cell
       Complex, dimension(3)                   :: MIV_PF

       !---- Local variables ----!
       real(kind=cp), DIMENSION (3)            :: QSV,Q,SV,X,Y,Z
       real(kind=cp), DIMENSION (3,3)          :: M
       INTEGER                                 :: i

       Q = Cart_Vector("R",H,Cell)
       SV = Cart_Vector("R",SPV,Cell)

       X = Q / sqrt(dot_product(Q,Q))
       QSV= Cross_Product(Q,SV)
       Z = QSV / sqrt(dot_product(QSV,QSV))
       Y = Cross_Product(Z,X)

       DO i = 1,3
          M(1,i) = 0.0    ! X-Component of Magnetic Interaction Vector is always equal to ZERO because
                          ! of MRI = Q x (M(Q) x Q); where M(Q) is the Fourier Transform of Magnetic Density of the sample
          M(2,i) = Y(i)
          M(3,i) = Z(i)
       END DO

       MIV_PF = MATMUL(M, MIV)
       
       return
    End Function  Magn_Inter_Vec_Pf
    
    !!--++
    !!--++ Function Mm(Nsf, Miv_Pf) Result(Mmc)
    !!--++    Complex, dimension(3), intent( in) :: MIV_PF !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                      :: MMC    !  Out -> magnetic magnetic interference terme
    !!--++
    !!--++    (Private)
    !!--++    Calculates the magnetic magnetic interfernce contribution to scattering in the
    !!--++    polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update April - 2005
    !!
    Function Mm(Miv_Pf) Result(Mmc)
       !---- Argument ----!
       Complex, dimension(3), intent( in) :: MIV_PF
       real(kind=cp)                      :: MMC
       
       !---- Local variables ----!

       MMC = REAL(MIV_PF(2) * CONJG(MIV_PF(3)) + CONJG(MIV_PF(2)) * MIV_PF(3))
       
       return
    End Function  Mm

    !!--++
    !!--++ Function Nuc_Contr(Nsf, Miv_Pf) Result(Nsc)
    !!--++    Complex,               intent( in):: NSF     !  In  -> Nuclear Structure Factor
    !!--++    Complex, dimension(3), intent( in):: MIV_PF  !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                     :: NSC     !  Out -> nuclear scattering contribution
    !!--++
    !!--++    (Private)
    !!--++    Calculates the nuclear contribution to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Nuc_Contr(Nsf) Result(Nsc)
       !---- Argument ----!
       Complex, intent( in)                     :: NSF
       real(kind=cp)                            :: NSC
       
       !---- Local variables ----!

       NSC = NSF * CONJG(NSF)
       
       return
    End Function  Nuc_Contr

    !!--++
    !!--++ Function Real_Nm_Y(Nsf, Miv_Pf) Result(R_Nm_Y)
    !!--++    Complex,               intent(in)  :: NSF     !  In  -> Nuclear Structure Factor
    !!--++    Complex, dimension(3), intent(in)  :: MIV_PF  !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                      :: R_NM_Y  !  Out -> real part of nuclear magnetic interference contribution along Y
    !!--++
    !!--++    (Private)
    !!--++    Calculates the real part of the nuclear magnetic interference contribution along Y to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Real_Nm_Y(Nsf, Miv_Pf) Result(R_Nm_Y)
       !---- Argument ----!
       Complex,               intent(in)  :: NSF
       Complex, dimension(3), intent(in)  :: MIV_PF
       real(kind=cp)                      :: R_NM_Y
       
       !---- Local variables ----!

       R_NM_Y = REAL(NSF * CONJG(MIV_PF(2)) + CONJG(NSF) * MIV_PF(2))
       
       return
    End Function  Real_Nm_Y

    !!--++
    !!--++ Function Real_Nm_Z(Nsf, Miv_Pf) Result(R_Nm_Z)
    !!--++    Complex, intent( in)               :: NSF            !  In  -> Nuclear Structure Factor
    !!--++    Complex, dimension(3), intent( in) :: MIV_PF         !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                      :: R_NM_Z         !  Out -> nuclear real part of magnetic interference contribution along Z
    !!--++
    !!--++    (Private)
    !!--++    Calculates the real part of the nuclear magnetic interference contribution along Z
    !!--++    to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Real_Nm_Z(Nsf, Miv_Pf) Result(R_Nm_Z)
       !---- Argument ----!
       Complex, intent( in)                     :: NSF
       Complex, dimension(3), intent( in)       :: MIV_PF
       real(kind=cp)                            :: R_NM_Z
       
       !---- Local variables ----!

       R_NM_Z = REAL(NSF * CONJG(MIV_PF(3)) + CONJG(NSF) * MIV_PF(3))
       
       return
    End Function  Real_Nm_Z

    !!--++
    !!--++ Function Tchiral(Nsf, Miv_Pf) Result(Tc)
    !!--++    Complex, dimension(3), intent( in)  :: MIV_PF !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                       :: TC     !  Out -> chiral contribution
    !!--++
    !!--++    (Private)
    !!--++    Calculates the chiral contribution to scattering in the polarisation coordinate frame
    !!--++    according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Tchiral(Miv_Pf) Result(Tc)
       !---- Argument ----!
       Complex, dimension(3), intent( in) :: MIV_PF
       real(kind=cp)                      :: TC
       
       !---- Local variables ----!

       TC = - AIMAG(MIV_PF(2) * CONJG(MIV_PF(3)) - CONJG(MIV_PF(2)) * MIV_PF(3))
       
       return
    End Function  Tchiral

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Set_Polar_Info(Cell, H, Spv, Pin, Nsf, Miv, Polari)
    !!----    Type (Crystal_Cell_Type),  intent(in)    :: Cell   !  In -> Cell variable
    !!----    real(kind=cp), DIMENSION (3),intent(in)  :: H      !  In -> Scattering vector in hkl
    !!----    real(kind=cp), dimension(3), intent( in) :: SPV    !  In -> Second Scattering plane vector in hkl
    !!----    real(kind=cp), intent( in)               :: Pin    !  In -> magnitude of initial polarisation
    !!----    COMPLEX, intent( in)                     :: NSF    !  In -> Nuclear Scattering Factor
    !!----    COMPLEX, dimension(3), intent( in)       :: MIV    !  In -> Magnetic interaction vector
    !!----    Type (Polar_Info_type), intent( out)     :: Polari !  Out ->type with all information about polarisation in
    !!----                                                                 one point hkl
    !!----
    !!----    Initializes the polarisation info type
    !!----
    !!---- Update: April - 2008
    !!
    Subroutine Set_Polar_Info(Cell, H, Spv, Pin, Nsf, Miv, Polari)
       !---- Arguments ----!
       Type (Crystal_Cell_Type),  intent(in)         :: Cell
       real(kind=cp), DIMENSION (3),intent(in)       :: H
       real(kind=cp), dimension(3), intent( in)      :: SPV
       real(kind=cp), intent( in)                    :: Pin
       COMPLEX, intent( in)                          :: NSF
       COMPLEX, dimension(3), intent( in)            :: MIV
       Type (Polar_Info_type), intent( out)          :: Polari

       !---- Local variables ----!
       real(kind=cp), DIMENSION (3)     :: sigma        ! elastic cross for different inicdent polarisation directions
       real(kind=cp)                    :: nc, my, mz, rnmy, rnmz, inmy, inmz, tc, mmc, A !the different contribution to cross-section
       COMPLEX, DIMENSION (3)           :: MIV_PF       !MIV in polarisation frame


       A = tpi**3/Cell%CellVol
       
       !First store given info in Polari
       Polari%H = H
       Polari%SPV = SPV
       Polari%Cell = Cell
       Polari%P = Pin
       Polari%NSF = NSF

       !Calculate the rest and also store it in Polari

       !magnetic interaction in polarisation frame
       MIV_PF = Magn_Inter_Vec_PF(MIV,H,SPV, Cell)
       Polari%MIV = MIV_PF

       !the different contributions to the scattering cross-section
       nc = nuc_contr(NSF)
       Polari%NC = A * nc

       my = mag_y(MIV_PF)
       Polari%MY = A * my

       mz = mag_z(MIV_PF)
       Polari%MZ = A * mz

       rnmy = real_nm_y(NSF, MIV_PF)
       Polari%RY = A * rnmy

       rnmz = real_nm_z(NSF, MIV_PF)
       Polari%RZ = A * rnmz

       inmy = im_nm_y(NSF, MIV_PF)
       Polari%IY = A * inmy

       inmz = im_nm_z(NSF, MIV_PF)
       Polari%IZ = A * inmz

       tc = tchiral(MIV_PF)
       Polari%TC = tc

       mmc = mm(MIV_PF)
       Polari%MM = A * mmc

       !scattering cross-section for the different initial polarisation vectors
       sigma = (/ nc + my + mz - Pin * tc, nc + my + mz + Pin * rnmy, nc + my + mz + Pin * rnmz /)
       Polari%CS = sigma

       !the polar matrix
       Polari%Pij(1,1) = ((nc - my - mz)* Pin + tc)/sigma(1)
       Polari%Pij(1,2) = (inmz * Pin + tc)/sigma(2)
       Polari%Pij(1,3) = (-inmy * Pin + tc)/sigma(3)
       Polari%Pij(2,1) = (-inmz * Pin + rnmy)/sigma(1)
       Polari%Pij(2,2) = ((nc + my - mz) * Pin + rnmy)/sigma(2)
       Polari%Pij(2,3) = (mmc * Pin + rnmy)/sigma(3)
       Polari%Pij(3,1) = (inmy * Pin + rnmz)/sigma(1)
       Polari%Pij(3,2) = (mmc * Pin + rnmz)/sigma(2)
       Polari%Pij(3,3) = ((nc - my + mz) * Pin + rnmz)/sigma(3)

       return
    End Subroutine Set_Polar_Info

    !!----
    !!---- Subroutine Write_Polar_Info(Polari, Lun, Info)
    !!----    Type (Polar_Info_type),     intent( in) :: Polrari !  in ->type with all information about polarisation in one point hkl
    !!----    integer,          optional, intent(in)  :: lun     !  In -> Unit to write
    !!----    character(len=*), optional, intent( in) :: info    !  in -> if info "P" also print information about coordinate frame
    !!----                                                       !        if info "C" also print information about crystal
    !!----                                                       !        if info "B" also print information about both
    !!----
    !!----    Outputs the polarisation info type in nice form
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Write_Polar_Info(Polari, Lun, info)
       !---- Arguments ----!
       Type (Polar_Info_type),    intent( in)  :: Polari !
       integer, optional,         intent(in)   :: Lun
       character(len=*),OPTIONAL, intent( in)  :: info
       !---- Local variables ----!
       integer            :: iunit

       iunit=6
       if (present(lun)) iunit=lun

       Write(unit=iunit,fmt="(/,a)")            "        Polar information:"
       Write(unit=iunit,fmt="(a,/)")            "        -------------------"
       Write(unit=iunit,fmt="(a,/)")            " => Initial parameters:"
       Write(unit=iunit,fmt="(3(a,f12.3), a)")  "    Scattering vector  (Qh,Qh,Ql) =  (", Polari%h(1) ,", ", Polari%h(2) , &
                                                ", ", Polari%h(3), ")"
       Write(unit=iunit,fmt="(3(a,f12.3), a)")  "    Add. in-plane vector      SPV =  (", Polari%SPV(1) ,", ", Polari%SPV(2) ,&
                                                ", ", Polari%SPV(3), ")"
       Write(unit=iunit,fmt="(a,f12.3)")        "              Polarisation degree =   ", Polari%P
       Write(unit=iunit,fmt="(a,/)")            " "
       IF (present(info)) THEN
         IF (info == "C" .OR. info == "c" .OR. info == "B" .OR. info == "b") THEN
           Write(unit=iunit,fmt="(a,/)")        " => Crystal information:"
           Write(unit=iunit,fmt="(3(a,f12.4))") "      a = ", Polari%Cell%cell(1),"      b = ", Polari%Cell%cell(2), "      c = ",&
                                                Polari%Cell%cell(3)
           Write(unit=iunit,fmt="(3(a,f12.3))") "  alpha = ", Polari%Cell%ang(1) ,"   beta = ", Polari%Cell%ang(2) , "  gamma = ",&
                                                Polari%Cell%ang(3)
           Write(unit=iunit,fmt="(a,f12.4)")    "                     Direct Cell Volume = ",Polari%Cell%CellVol
           Write(unit=iunit,fmt="(a,/)")        ""
         END IF
         IF (info == "P" .OR. info == "p" .OR. info == "B" .OR. info == "b") THEN
           Write(unit=iunit,fmt="(a,/)")  " => Polarisation coordinate frame:"
           Write(unit=iunit,fmt="(a,/)")  " "
           Write(unit=iunit,fmt="(a,/)")  "    Polarisation coordinate frame according to Blume"
           Write(unit=iunit,fmt="(a,a,/)")"    X  || scattering vector Q    ",&
                                          " (where Q is the scattering Vector in cartesian real space coordinates)"
           Write(unit=iunit,fmt="(a,/)")  "    Y _|_ scattering vector Q in scattering plane"
           Write(unit=iunit,fmt="(a,/)")  "    Z _|_ scattering vector Q out of scattering plane"
           Write(unit=iunit,fmt="(a,/)")  "    (ATTENTION: This choice is not non-ambiguous, there are always two possible choices"
           Write(unit=iunit,fmt="(a,/)")  "    for a right handed coordinate frame which will fullfils this condition!!!)"
           Write(unit=iunit,fmt="(a,/)")  " "
           Write(unit=iunit,fmt="(a,/)")  "    Therefore the right handed coordinate frame is explicitly chosen like this:"
           Write(unit=iunit,fmt="(a,a,/)")"    X := Q/|Q|                 where Q is the scattering Vector in ", &
                                          "cartesian real space coordinates"
           Write(unit=iunit,fmt="(a,a,/)")"    Z := (Q x SVP)/|(Q x SVP)| where SVP is a second vector in the scattering ",&
                                          "plane in cartesian real space coordinates"
           Write(unit=iunit,fmt="(a,/)")  "    Y := (Z x X)"
           Write(unit=iunit,fmt="(a,/)")  " "
           Write(unit=iunit,fmt="(a)")    "                            Y"
           Write(unit=iunit,fmt="(a)")    "                          /|\"
           Write(unit=iunit,fmt="(a)")    "                           |"
           Write(unit=iunit,fmt="(a)")    "                   Q       | Z "
           Write(unit=iunit,fmt="(a)")    "            ____________ _\o_____\ X"
           Write(unit=iunit,fmt="(a)")    "            \             /      /"
           Write(unit=iunit,fmt="(a)")    "             \           / "
           Write(unit=iunit,fmt="(a)")    "              \         /"
           Write(unit=iunit,fmt="(a)")    "               \       /"
           Write(unit=iunit,fmt="(a)")    "                \     /  K_f"
           Write(unit=iunit,fmt="(a)")    "             K_i \   /"
           Write(unit=iunit,fmt="(a)")    "                 _\|/_"
           Write(unit=iunit,fmt="(a,/)")  ""
         END IF
       END IF
       
       Write(unit=iunit,fmt="(a,/)")         " => Interaction potentials:"
       Write(unit=iunit,fmt="(2(a,f7.3))")   "       NSF = ", real(Polari%NSF)," + i ", AIMAG(Polari%NSF)
       Write(unit=iunit,fmt="(2(a,3f7.3),a)")"       MiV = (", REAL(Polari%MIV), ") + i (", AIMAG(Polari%MIV),")"
       Write(unit=iunit,fmt="(a,/)")         " "
       Write(unit=iunit,fmt="(a,/)")         " => Different contributions to the cross-section:"
       Write(unit=iunit,fmt="(a,f12.3)")     "         Nuclear Contribution = ", Polari%NC
       Write(unit=iunit,fmt="(a,f12.3)")     "             Magnetic along y = ", Polari%MY
       Write(unit=iunit,fmt="(a,f12.3)")     "             Magnetic along z = ", Polari%MZ
       Write(unit=iunit,fmt="(a,f12.3)")     "  Real nuclear magnetic al. y = ", Polari%RY
       Write(unit=iunit,fmt="(a,f12.3)")     "  Real nuclear magnetic al. z = ", Polari%RZ
       Write(unit=iunit,fmt="(a,f12.3)")     "  Imag nuclear magnetic al. y = ", Polari%IY
       Write(unit=iunit,fmt="(a,f12.3)")     "  Imag nuclear magnetic al. z = ", Polari%IZ
       Write(unit=iunit,fmt="(a,f12.3)")     "          Chiral contribution = ", Polari%TC
       Write(unit=iunit,fmt="(a,f12.3)")     "            Magnetic Magnetic = ", Polari%MM
       Write(unit=iunit,fmt="(a,/)")         " "
       Write(unit=iunit,fmt="(a,/)")           " => Cross-section for initial polar vector:"
       Write(unit=iunit,fmt="(a,f12.3)")       "                      along x = ", Polari%CS(1)
       Write(unit=iunit,fmt="(a,f12.3)")       "                      along y = ", Polari%CS(2)
       Write(unit=iunit,fmt="(a,f12.3)")       "                      along z = ", Polari%CS(3)
       Write(unit=iunit,fmt="(a,/)")           " "
       Write(unit=iunit,fmt="(a,/)")           " => Polarisation tensor as it will be measured:"
       Write(unit=iunit,fmt="(a,/)")           " "
       Write(unit=iunit,fmt="(3(a,f12.4), a)") "        /", Polari%Pij(1,1),"  ", Polari%Pij(1,2) , "  ", Polari%Pij(1,3), "  \"
       Write(unit=iunit,fmt="(3(a,f12.4), a)") "  PT  = |", Polari%Pij(2,1),"  ", Polari%Pij(2,2) , "  ", Polari%Pij(2,3), "  |"
       Write(unit=iunit,fmt="(3(a,f12.4), a)") "        \", Polari%Pij(3,1),"  ", Polari%Pij(3,2) , "  ", Polari%Pij(3,3), "  /"

       return
    End Subroutine Write_Polar_Info

    !!----
    !!---- Subroutine Write_Polar_Line(Polari, Lun)
    !!----    Type (Polar_Info_type), intent( in)     :: Polrari !  in ->type with all information about polarization in one point hkl
    !!----    integer, optional,      intent(in)      :: lun     !  In -> Unit to write
    !!----
    !!----    Outputs the polarization info type in line form, so you can write it to a file
    !!----
    !!---- Update: May - 2005
    !!
    Subroutine Write_Polar_Line(Polari, Lun)
       !---- Arguments ----!
       Type (Polar_Info_type), intent( in)     :: Polari !
       integer, optional,      intent(in)      :: Lun
       
       !---- Local variables ----!
       integer            :: iunit

       iunit=6
       if (present(lun)) iunit=lun

       Write(unit=iunit,fmt="(/,a)")        "     H         K         L         NSF^2     NSF_r     NSF_i"
       Write(unit=iunit,fmt="(6(f10.6))")   Polari%H(1), Polari%H(2), Polari%H(3), Polari%NC, real(Polari%NSF), AIMAG(Polari%NSF)
       Write(unit=iunit,fmt="(a)")          "    MIV^2     MIV_xr    MIV_xi    MIV_yr    MIV_yi    MIV_zr    MIV_zi"
       Write(unit=iunit,fmt="(7(f10.6))")   (Polari%MY+Polari%MZ), Polari%MIV(:)
       Write(unit=iunit,fmt="(/,a)")        "      CSx       CSy       CSz       NC        My        Mz"
       Write(unit=iunit,fmt="(6(f10.6))")   Polari%CS(1:3), Polari%NC, Polari%MY, Polari%MZ
       Write(unit=iunit,fmt="(a)")          "      RNMy      RNMz      INMy      INMz       CC       MM"
       Write(unit=iunit,fmt="(6(f10.6))")   Polari%RY, Polari%RZ, Polari%IY, Polari%IZ, Polari%TC, Polari%MM
       Write(unit=iunit,fmt="(/,a)")        "       Pix       Piy       Piz       Pfx       Pfy       Pfz"
       Write(unit=iunit,fmt="(6(f10.3))")   Polari%P, 0.0, 0.0, Polari%Pij(1,1), Polari%Pij(2,1), Polari%Pij(3,1)
       Write(unit=iunit,fmt="(6(f10.3))")   0.0,Polari%P, 0.0, Polari%Pij(1,2), Polari%Pij(2,2), Polari%Pij(3,2)
       Write(unit=iunit,fmt="(6(f10.3))")   0.0, 0.0, Polari%P, Polari%Pij(1,3), Polari%Pij(2,3), Polari%Pij(3,3)
       Write(unit=iunit,fmt="(a,/)")        "-------------------------------------------------------------------------"

       return
    End Subroutine Write_Polar_line

 End Module CFML_Polarimetry