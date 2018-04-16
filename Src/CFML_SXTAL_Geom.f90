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
!!---- MODULE: CFML_Geometry_SXTAL
!!----   INFO: Module for making geometrical calculations in single crystal instruments.
!!----         All subroutines gathered in this module have been taken from packages
!!----         at ILL maintained informally by different people. When the original authors
!!----         have made some comments in the header of the routines these have been kept.
!!----
!!---- HISTORY
!!----    Update: 06/03/2011
!!----
!!--..
!!--..    The original routines have been transformed to Fortran 95 and the names
!!--..    have been changed to indicate more clearly the purpose.
!!--..
!!--..    Code converted using MODIF90 by Juan Rodriguez-Carvajal
!!--..    from manually modified files produced by a pre-treatment with the
!!--..    program To_F90, by Alan Miller.
!!--..    Date: 2005-03-06  Time: 22:07:48
!!--..
!!--..    The persons that, to my knowledge have contributed to maintain these procedures
!!--..    in some way are: R.F.D Stansfield, A. Barthelemy, Garry McIntyre, John Allibon,
!!--..    S.A. Mason, M.Thomas and Alain Filhol
!!--..
!!--..    Juan Rodriguez-Carvajal
!!--..
!!--..    The module is under test, I'm entirely responsible for all the errors and bugs.
!!--..
!!--..
!!--..    Some words about notation and conventions
!!--..
!!--..    A reflection in r.l.u. units is represented by a real(kind=cp) vector h
!!--..    The same vector referred to the Cartesian laboratory system when all
!!--..    motor angles are set to zero is noted as z1, so that [z1]= [UB] [h]
!!--..    where [] denotes a column matrices for 3-D vectors.
!!--..
!!--..    The rotational matrices [PHI], [CHI], [OMG], [PSI] are calculated
!!--..    according to conventions in Busing & Levy
!!--..    After applying a [PHI] matrix to vector [z1] we obtain [z2] ... and so on
!!--..     [z2] = [PHI] [z1]
!!--..     [z3] = [CHI] [z2]
!!--..     [z4] = [OMG] [z3]
!!--..
!!--..
!!--..
!!--..    The default instrument cartesian frame is that defined by
!!--..    W.R.Busing & H.A.Levy (Acta Cryst. 22,457-464 (1967)
!!--..    The reciprocal Cartesian Frame used by Busing-Levy has "x" along a*
!!--..    "y" in the (a*,b*) plane and "z" perpendicular to that plane. The
!!--..    matrix passing from (h,k,l) in r.l.u. to the cartesian frame is
!!--..           | a*   b*.cos(gamma*)   c*.cos(beta*)            |
!!--..    [B] =  | 0    b*.sin(gamma*)  -c*.sin(beta*).cos(alpha) |
!!--..           ! 0          0          1/c                      |
!!--..
!!--..        [hc] = [B] [h]
!!--..
!!--..    Let U be the orthogonal matrix relating the phi-axis system
!!--..    attached to the phi-shaft of the instrument with the Cartesian crystal system
!!--..    then     [h-phi]= [U] [hc]      [h-phi] = [U] [B] [h] = [UB] [h]  --> z1
!!--..
!!--..    Let us define three more cartesian systems attached to the Chi,Omega and
!!--..    Theta axes, coincident with phi-axis when all instrument angles are set
!!--..    to zero :
!!--..       [h-chi] = [PHI] [h-phi]     --> z2
!!--..       [h-omg] = [CHI] [h-chi]     --> z3
!!--..       [h-tet] = [OMG] [h-omg]     --> z4
!!--..
!!--..             |  cos(phi)      sin(phi)      0  |             |zL-axis
!!--..     [PHI] = | -sin(phi)      cos(phi)      0  |             |
!!--..             |     0            0           1  |             |
!!--..                                           incoming beam ->__|_________yL-axis
!!--..             |  cos(chi)       0     sin(chi)  |            /
!!--..     [CHI] = |     0           1       0       |           /
!!--..             | -sin(chi)       0     cos(chi)  |          /xL-axis
!!--..
!!--..             |  cos(omg)      sin(omg)      0  |
!!--..     [OMG] = | -sin(omg)      cos(omg)      0  |
!!--..             |     0            0           1  |
!!--..
!!--..    Let us define 2 more Cartesian frames: The laboratory fixed system
!!--..    and the 2Theta-axis system attached to the 2theta shaft. All coincide
!!--..    for all angles equal to zero.
!!--..     [h-Lab] = [THE]  [h-tet] = [N] [h-omg]
!!--..     [h-tte] = [THE]t [h-tet] = [M] [h-omg]
!!--..    Where:
!!--..             |  cos(the)      sin(the)      0  |
!!--..     [THE] = | -sin(the)      cos(the)      0  |
!!--..             |     0            0           1  |
!!--..
!!--..                         |  cos(nu)      sin(nu)      0  |
!!--..     [N] = [THE] [OMG] = | -sin(nu)      cos(nu)      0  |  with nu=omg+tet
!!--..                         |     0            0         1  |
!!--..
!!--..                          |  cos(mu)      sin(mu)      0  |
!!--..     [M] = [THE]t [OMG] = | -sin(mu)      cos(mu)      0  |  with mu=omg-tet
!!--..                          |     0            0         1  |
!!--..    All angles, except Chi, are left-handed rotations about their respective axes.
!!--..
!!--..    In this module the origin of the Omega angle is different from that of BL.
!!--..    The Omega angle is always evaluated from the YL laboratory axis instead of
!!--..    taking the origin from the YTheta axis. This means that the bisecting geometry
!!--..    condition is Omega=Theta instead of Omega=0. This makes a change in some BL
!!--..    equations to which we have to add the value of Theta.
!!--..
!!--..    Basic Diffractometer equations:
!!--..          q = |h| = sqrt(dot_product(hc,hc)), with [hc]=[B][h]
!!--..    Bragg equation: sin(theta) = Lambda.q/2
!!--..    The reflection [h] is in the diffraction position if
!!--..         [h-tet] = [OMG] [CHI] [PHI] [U] [B] [h]
!!--..    has the form [h-tet]t = (q, 0, 0)
!!--..
!!--..    The orientation matrix U can be obtained from two non-colinear reflections h1 & h2,
!!--..    provided the cell parameters are known
!!--..
!!--..    Two kind of Busing-Levy (BL) frames are possible:
!!--..    For both frames the origin is at sample position and the y-axis
!!--..    positive sense is along the secondary beam (from monochromator to sample)
!!--..    a)  z upward
!!--..    b)  z downward
!!--..    The x-axis makes a right-handed frame with the other axes. Positive 2theta
!!--..    angles are from y-axis towards x-axis
!!--..    A change of geometry is always possible by providing the components of the
!!--..    axes {e1,e2,e3} with respect to the standard {i,j,k} BL-frame
!!--..
!!--..    igeom=1: Bisecting (PSI=0)
!!--..    igeom=2: Bisecting - HiCHI
!!--..    igeom=3: Normal beam (Chi=-90)
!!--..    igeom=4: Parallel    (Chi=90, PSI=90)
!!--..
!!----
!!---- DEPENDENCIES
!!----
!!----
!!---- VARIABLES
!!--..    Types
!!----    PSD_VAL_TYPE
!!----    SXD_VAL_TYPE
!!--++    PSD_SET            [Private]
!!----    SXD
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ANGS_4C_BISECTING
!!----       CALANG
!!----       CALC_OM_CHI_PHI
!!----       CALC_PSI
!!----       CELL_FR_UB
!!----       CHI_MAT
!!----       PSD_CONVERT
!!----       D19PSD
!!----       DSPACE
!!----       EQUATORIAL_CHI_PHI
!!----       FIXDNU
!!----       FLAT_CONE_VERTDET
!!----       GENB
!!----       GENUB
!!----       GET_ANGS_NB
!!----       GET_DSPACING_THETA
!!----       GET_GAOMNU_FRCHIPHI
!!----       GET_UB_FROM_HKL_HKL_OMEGA
!!----       GET_UB_FROM_UVW_HKL_OMEGA
!!----       GET_WAVEGANU_FRZ4
!!----       GET_Z1_D9ANGLS
!!----       GET_Z1_FROM_PIXEL
!!----       NORMAL
!!----       NORMAL_BEAM_ANGLES
!!----       PHI_MAT
!!----       PSI_MAT
!!----       REFVEC
!!----       S4CNB
!!----       SET_PSD
!!----       SNB4C
!!----       SXDPSD
!!----       TRIPLE
!!----       Z1FRFC
!!----       Z1FRMD
!!----       Z1FRNB
!!----       Z1FRZ2
!!----       Z1FRZ3
!!----       Z1FRZ4
!!----       Z2FRZ1
!!----       Z3FRZ1
!!----       Z4FRGN
!!----       Z4FRZ1
!!----
!!
 Module CFML_Geometry_SXTAL
    !---- Use Modules ----!
    Use CFML_GlobalDeps,        Only: Cp,To_Deg,To_Rad
    Use CFML_Math_General,      Only: cosd,sind,atan2d,acosd,asind,tand,co_linear,in_limits
    Use CFML_Math_3D,           Only: Cross_Product, invert => invert_A
    Use CFML_Crystal_Metrics,   Only: Crystal_Cell_Type, Zone_Axis_type, Get_basis_from_uvw, &
                                      Rot_Matrix
    Use CFML_Geometry_Calc,     Only: Get_OmegaChiPhi, Get_Matrix_moving_v_to_u,ERR_Geom_Mess,&
                                      ERR_Geom
    Use CFML_ILL_Instrm_data,   Only: Current_Orient, Current_Instrm, SXTAL_Numor_type

    !---- Variables ----!
    Implicit None

    Private

    !---- List of public functions ----!
    Public :: Chkin180

    !---- List of public subroutines ----!
    Public :: Angs_4C_bisecting, Equatorial_Chi_Phi, Get_dspacing_theta,                        &
              Get_GaOmNu_frChiPhi, Chi_mat, Phi_mat, Psi_mat, Get_Angs_NB,                      &
              Calc_Om_Chi_Phi, Calc_Psi, d19psd, dspace, fixdnu, Normal_Beam_Angles,            &
              s4cnb, snb4c, Flat_Cone_vertDet, Get_WaveGaNu_frZ4, normal, refvec, sxdpsd,       &
              triple, z3frz1, z2frz1, z1frfc, z1frnb, z1frmd, z1frz4, z1frz3, z1frz2, z4frgn,   &
              z4frz1, calang, genb, genub, cell_fr_UB, set_psd, get_z1_from_pixel,              &
              Get_z1_D9angls, psd_convert, Get_UB_from_uvw_hkl_omega, Get_UB_from_hkl_hkl_omega,&
              Get_FlatCone_Angles_D10

    !---- Definitions ----!

    !!----
    !!---- ERR_SXTGEOM
    !!----    logical, public :: err_sxtgeom
    !!----
    !!----    Logical Variable indicating an error in CFML_Geometry_SXTAL module
    !!----
    !!---- Update: February - 2010
    !!
    logical, public :: ERR_SXTGeom

    !!----
    !!---- ERR_SXTGEOM_MESS
    !!----    character(len=150), public :: ERR_SXTGeom_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2010
    !!
    character(len=150), public :: ERR_SXTGeom_Mess

    !!----
    !!---- TYPE :: PSD_VAL_TYPE
    !!--..
    !!---- Type, public :: Psd_Val_Type       !Set of caracteristic parameters of the PSD
    !!----    Character(len=12) :: name_inst  ! Name of the instrument
    !!----    real(kind=cp)     :: xoff       !
    !!----    real(kind=cp)     :: zoff       !
    !!----    real(kind=cp)     :: radius     !
    !!----    real(kind=cp)     :: yoff       !
    !!----    real(kind=cp)     :: cgap       !
    !!----    real(kind=cp)     :: agap       !
    !!----    Integer           :: ncat       !
    !!----    Integer           :: nano       !
    !!----    Integer           :: ipsd       !
    !!---- End Type Psd_Val_Type
    !!----
    !!---- Update: July - 2010
    !!
    Type, Public :: Psd_Val_Type
       Character(len=12) :: name_inst ! Name of the instrument
       real(kind=cp)     :: xoff
       real(kind=cp)     :: zoff
       real(kind=cp)     :: radius
       real(kind=cp)     :: yoff
       real(kind=cp)     :: cgap
       real(kind=cp)     :: agap
       Integer           :: ncat
       Integer           :: nano
       Integer           :: ipsd
    End Type Psd_Val_Type

    !!----
    !!---- PSD
    !!----    Type(Psd_Val_Type), Public :: psd
    !!----
    !!----
    !!---- Update: January - 2009
    !!
    Type(Psd_Val_Type), Public :: psd

    !!----
    !!---- TYPE :: Sxd_Val_Type
    !!--..
    !!---- Type, public :: Sxd_Val_Type        !Set of caracteristic parameters of the diffractometr SXD (old?) at ISIS
    !!----    real(kind=cp)      :: distms         !Not tested! (used only in the subroutine sxdpsd)
    !!----    real(kind=cp)      :: distsd         !
    !!----    real(kind=cp)      :: dimx           !
    !!----    real(kind=cp)      :: dimz           !
    !!----    real(kind=cp)      :: xoff           !
    !!----    real(kind=cp)      :: yoff           !
    !!----    real(kind=cp)      :: zoff           !
    !!----    real(kind=cp)      :: toff           !
    !!----    real(kind=cp)      :: velcon         !
    !!----    Integer            :: nxcel          !
    !!----    Integer            :: nzcel          !
    !!---- End Type Sxd_Val_Type
    !!----
    !!---- Update: April - 2008
    !!
    Type, Public  :: Sxd_Val_Type
       real(kind=cp)      :: distms
       real(kind=cp)      :: distsd
       real(kind=cp)      :: dimx
       real(kind=cp)      :: dimz
       real(kind=cp)      :: xoff
       real(kind=cp)      :: yoff
       real(kind=cp)      :: zoff
       real(kind=cp)      :: toff
       real(kind=cp)      :: velcon
       Integer            :: nxcel
       Integer            :: nzcel
    End Type Sxd_Val_Type

    !!--++
    !!--++ PSD_SET
    !!--++    logical, private :: psd_set
    !!--++
    !!--++    ...No comments!!!
    !!--++
    !!--++ Update: April - 2008
    !!
    logical, private :: psd_set=.false.

    !!----
    !!---- SXD
    !!----    Type(Sxd_Val_Type), Public :: sxd
    !!----
    !!----    ...No comments!!!
    !!----
    !!---- Update: April - 2008
    !!
    Type(Sxd_Val_Type), Public :: sxd

 Contains

    !!---- Elemental Function chkin180(angle) result(angle_in)
    !!----   real(kind=cp), intent(in) :: angle
    !!----   real(kind=cp)             :: angle_in
    !!----
    !!----  Elemental function (can be used with arrays) putting the
    !!----  input angle (in degrees) into the interval (-180,180)
    !!----
    !!----  Created: March 2013 (JRC)
    !!----
    Elemental Function chkin180(angle) result(angle_in)
      real(kind=cp), intent(in) :: angle
      real(kind=cp)             :: angle_in

      angle_in=angle
      do
        if(angle_in > 180.0_cp) then
          angle_in=angle_in-360.0_cp
        else if(angle_in < -180.0_cp) then
          angle_in=angle_in+360.0_cp
        else
          exit
        end if
      end do
      return
    End Function chkin180

    !!----
    !!---- Subroutine Angs_4c_bisecting(wave,z1,tth,om,ch,ph,ierr)
    !!----    real(kind=cp), Intent(In)                :: wave
    !!----    real(kind=cp), Intent(In), Dimension(3)  :: z1
    !!----    real(kind=cp), Intent(Out)               :: tth,om,ch,ph
    !!----    Integer, Intent(Out)                     :: ierr
    !!----
    !!----    Calculate 2-THETA, OMEGA (=THETA), CHI, PHI to put the
    !!----    vector Z1 in the bisecting diffraction condition. The
    !!----    reciprocal vector Z1 is given in cartesian components
    !!----    with respect to the laboratory system. This geometry
    !!----    corresponds to the bisecting PSI=0 (igeom=1).
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Angs_4c_bisecting(wave,z1,tth,om,ch,ph,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                :: wave
       real(kind=cp), Intent(In), Dimension(3)  :: z1
       real(kind=cp), Intent(Out)               :: tth,om,ch,ph
       Integer, Intent(Out)                     :: ierr

       !---- Local Variables ----!
       real(kind=cp) :: th, ds

       ierr=0
       Call Equatorial_Chi_Phi(z1,ch,ph) !Eq. 38 of Busing & Levy
       Call Get_dspacing_theta(wave,z1,ds,th,ierr)
       If(ierr == 0) Then
          om=th             !Bisecting position
          tth=th*2.0
       Else
          tth=0.0  !Case of error (ierr=1, null z1; ierr=2 outside resolution sphere)
          om=0.0
          ch=0.0
          ph=0.0
       End If

       Return
    End Subroutine Angs_4c_bisecting

    !!----
    !!---- Subroutine calang(h,tteta,om,ch,ph,ierr,wav,ubm,geom)
    !!----    real(kind=cp),Dimension(3),             Intent(In) :: h
    !!----    real(kind=cp),                          Intent(Out):: tteta,om,ch,ph
    !!----    Integer,                                Intent(Out):: ierr
    !!----    real(kind=cp),                optional, intent(in) :: wav
    !!----    real(kind=cp), dimension(3,3),optional, intent(in) :: ubm
    !!----    integer,                      optional, intent(in) :: geom
    !!----
    !!--<<    Original comments:
    !!----    SUBROUTINE *** CALANG ***
    !!----    VERSION DU 6/10/1976
    !!----    CALCULE LES ANGLES 2THETA OMEGA CHI PHI  :
    !!----           IPARA=1  GEOMETRIE BISSECTRICE (PSI=0)
    !!----           IPARA=2  GEOMETRIE BISECTING - HICHI  (D9, displex,cryostat)
    !!----           IPARA=3  GEOMETRIE NORMAL BEAM
    !!----           IPARA=4  GEOMETRIE PARALLELE (PSI=90)  (D15,D16)
    !!----      POUR REFLECTION H,K,L ET MATRICE UB DONNEES
    !!----
    !!----    This subroutine is a more general variant of Angs_4c_bisecting.
    !!----    If the optional arguments are given, the corresponding values
    !!-->>    are adopted instead of those of the current instrument.
    !!----
    !!---- Update: July 2008
    !!
    Subroutine calang(h,tteta,om,ch,ph,ierr,wav,ubm,geom)
       !---- Arguments ----!
       real(kind=cp),Dimension(3),             Intent(In) :: h
       real(kind=cp),                          Intent(Out):: tteta,om,ch,ph
       Integer,                                Intent(Out):: ierr
       real(kind=cp),                optional, intent(in) :: wav
       real(kind=cp), dimension(3,3),optional, intent(in) :: ubm
       integer,                      optional, intent(in) :: geom

       !--- Local Variables ---!
       real(kind=cp)                 :: ttmax,ttmin,sint,wave,chmax,diff,ds,theta
       real(kind=cp), Dimension(3)   :: z1
       real(kind=cp), Dimension(3,3) :: ub
       Integer                       :: igeom

       !---- Set local variables with the current instrument
       if (present(wav) .and. present(ubm) .and. present(geom)) then
          wave = wav
          ub   = ubm
          ttmin= 0.0
          ttmax= 180.0
          igeom= geom
          if(igeom == 2) then  !high-chi
             chmax= 200.0
          else
             chmax= 170.0
          end if
       else
          wave = Current_Orient%wave
          ub   = Current_Orient%ub
          ttmin= Current_Instrm%Ang_Limits(1,1)  !    1       2      3     4   ...
          ttmax= Current_Instrm%Ang_Limits(1,2)  ! 2theta   Omega   Chi   Phi  ...
          igeom= Current_Instrm%igeom
          chmax= Current_Instrm%Ang_Limits(3,2)
       end if
       if(present(geom)) igeom=geom
       ierr=0

       z1=Matmul(ub,h)
       sint=0.5*wave*Sqrt(Dot_Product(z1,z1))
       If (abs(sint) > 1.0) Then
          ierr=1
          Return
       End If
       theta =asind(sint)
       om=theta       !Theta = omega in bisecting geometry
       tteta=2.0*theta

       If (tteta < ttmin .or. tteta > ttmax) Then
          ierr=1
          Return
       End If

       ds = Sqrt(z1(1)*z1(1)+z1(2)*z1(2))
       Select Case(igeom)
          Case(1,2)    !1:Bisecting Geometry  (PSI=0), 2:Bisecting Geometry - High CHI
             ph = atan2d(z1(2),z1(1))  !Eqs 38 BL
             ch = atan2d (z1(3),ds)
             If(igeom == 2) Then
               ch=180.0-ch  !Eqs 40 BL
               ph=180.0+ph
             End If

             If(ph >  180.0) ph=ph-360.0
             If(ph < -180.0) ph=ph+360.0
             If(ch >  180.0) ch=ch-360.0
             If(ch < -180.0) ch=ch+360.0
             If (chmax > 180.0) Then
                diff=chmax-180.0
                If (ch < -180.0 + diff .AND. ch >= -180.0)  ch=ch+360.0
             End If

          Case(3)    !3:  Geometry NORMAL BEAM
             ch=-90.0
             ph=-atan2d(z1(1),z1(2))
             om=theta+90.0-atan2d(z1(3),ds)

          Case(4)    !4:  Geometry parallel (PSI=90) ..D15-D16  (PSI=90)
             ch=+90.0
             ph=atan2d(z1(1),-z1(2))     ! Eqs 41 BL (with different definition of omega
             om=theta+atan2d(-ds,z1(3))  ! Omega=Theta+atan2d(-ds,z1(3))

       End Select

       Return
    End Subroutine calang

    !!----
    !!---- Subroutine Calc_Om_Chi_Phi(vhkl,vlab1,psi,ub,om,ch,ph,ierr)
    !!----    real(kind=cp), Intent(In),Dimension(3)     :: vhkl,vlab1
    !!----    real(kind=cp), Intent(In)                  :: psi
    !!----    real(kind=cp), Intent(In),Dimension(3,3)   :: ub
    !!----    real(kind=cp), Intent(In Out)              :: om,ch,ph
    !!----    Integer, Intent(Out)                       :: ierr
    !!----
    !!----    Calculate OM, CH, PH for diffraction vector at azimuthal angle PSI
    !!----    from the diffraction condition expressed as :
    !!----    [TZ]=[OM].[CH].[PH].[TS].[PSI]-1 if [R]=[OM].[CH].[PH]
    !!----    then [R]'=[TZ].[PSI].[TS]-1
    !!----    The om,ch,ph angles are provided, on input, to calculate the components
    !!----    of the vector vlab1 in the Theta-system for Psi=0
    !!----    Used only in the Flat_Cone_vertDet subroutine.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Calc_Om_Chi_Phi(vhkl,vlab1,psi,ub,om,ch,ph,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In),Dimension(3)     :: vhkl,vlab1
       real(kind=cp), Intent(In)                  :: psi
       real(kind=cp), Intent(In),Dimension(3,3)   :: ub
       real(kind=cp), Intent(In Out)              :: om,ch,ph
       Integer, Intent(Out)                       :: ierr

       !--- Local Variables ---!
       real(kind=cp),Dimension(3)   :: z1,z4,vs,vz
       real(kind=cp),Dimension(3,3) :: r,dum1,dum2,ts,tz

       z1=vlab1
       Call z4frz1(z1,om,ch,ph,z4)
       Call refvec(vhkl,ub,vs,vz,ierr)
       Call triple(z1,vs,ts,ierr)
       Call triple(z4,vz,tz,ierr)
       If(ierr /= 0) Return
       r=Transpose(ts)
       Call Psi_mat(psi,dum2)
       dum1=Matmul(dum2,r)
       r=Matmul(tz,dum1)

       !---- Solve R for new Eulerian components
       If (Abs(r(3,3)) < 1.0) Then

          !---- We choose 0.0<= CH <=180.0 but remember that an equivalent
          !---- solution is CH'=-CH, PH'=180.0+PH, OM'=180.0+OM
          om=atan2d(-r(2,3),r(1,3))
          ch=atan2d( Sqrt(r(1,3)*r(1,3)+r(2,3)*r(2,3)), r(3,3) )
          ph=atan2d( -r(3,2), -r(3,1) )
       Else
          !---- The singular case of CH=0; OM and PH no longer independent
          !---- Choose OM=90.0+THE angle that bisects projected INC. and DIF. beams
          Call z4frz1(vlab1,om,ch,ph,z4)
          om=90.0+atan2d(-z4(2),z4(1))
          If(r(3,3) > 0.0) ch=0.0
          If(r(3,3) < 0.0) ch=180.0
          ph=Sign(1.0_cp,r(3,3))*(atan2d(r(1,2),r(2,2))-om)

          !---- Convert -180.0 <= PH <= +180.0
          ph=ph-360.0_cp*Int((Sign(180.0_cp,ph)+ph)/360.0_cp)
       End If

       Return
    End Subroutine Calc_Om_Chi_Phi

    !!----
    !!---- Subroutine Calc_Psi(vhkl,vlab1,om,ch,ph,ub,psi,ierr)
    !!----    real(kind=cp), Intent(In),Dimension(3)        :: vhkl,vlab1
    !!----    real(kind=cp), Intent(In)                     :: om,ch,ph
    !!----    real(kind=cp), Intent(In),Dimension(3,3)      :: ub
    !!----    real(kind=cp), Intent(Out)                    :: psi
    !!----    Integer, Intent(Out)                          :: ierr
    !!----
    !!--<<    Calculate PSI for the reflection VHKL positioned at OM, CH and PH.
    !!----    The value of PSI is taken to be zero when (the values of OM, CH, PH
    !!----    are such that) the reflection H0 lies in the plane of VHKL and VZ,
    !!----    on the same side of VHKL as VZ.
    !!----    The reference vectors H0 and VZ are defined in subroutine REFVEC.
    !!----    there, the vector VZ is the z-axis of the fixed laboratory system
    !!----    (Busing and Levy Convention, Y along beam, X in positive 2-THETA
    !!----    direction). H0 is (0,0,1) for all VHKL except when VHKL is parallel
    !!----    to (0,0,1), in which case (0,1,0) is chosen.
    !!----    Construct two orthonormal-vector-triplets, TS in the PHI-axis
    !!----    system with normalised VLAB1, (VLAB1 x TS2) x VLAB1, VLAB1 x TS2 as
    !!----    columns; and TZ in the fixed-lab system with normalised [R].VLAB1,
    !!----    ([R].VLAB1 x TZ2) x [R].VLAB1, [R].VLAB1 x TZ2 as columns.
    !!----    Then the diffraction condition can be expressed as:
    !!----    [TZ]=[OM].[CH].[PH].[TS].[PSI]-1  (all are orthonormal matrices)
    !!----    hence [PSI]=[TZ]-1.[R].[TS] where [R]=[OM].[CH].[PH]
    !!----                                             R.F.D.Stansfield May-83
    !!-->>    PSI matrix and equations corrected 1/9/83 - results not changed
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Calc_Psi(vhkl,vlab1,om,ch,ph,ub,psi,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In),Dimension(3)        :: vhkl,vlab1
       real(kind=cp), Intent(In)                     :: om,ch,ph
       real(kind=cp), Intent(In),Dimension(3,3)      :: ub
       real(kind=cp), Intent(Out)                    :: psi
       Integer, Intent(Out)                          :: ierr

       !--- Local Variables ---!
       real(kind=cp),Dimension(3)   ::  z1,z4,vs,vz
       real(kind=cp),Dimension(3,3) :: dum1,dum2,dum3,ts,tz

       z1=vlab1
       Call z4frz1(z1,om,ch,ph,z4)
       Call refvec(vhkl,ub,vs,vz,ierr)
       Call triple(z1,vs,ts,ierr)
       Call triple(z4,vz,tz,ierr)
       If(ierr /= 0) Return
       Call Phi_mat(ph,dum2)
       dum1=Matmul(dum2,ts)
       Call Chi_mat(ch,dum3)
       dum2=Matmul(dum3,dum1)
       Call Phi_mat(om,dum1)
       dum3=Matmul(dum1,dum2)
       dum2=Transpose(tz)
       dum1=Matmul(dum2,dum3)
       psi=atan2d(-dum1(2,3),dum1(2,2))

       Return
    End Subroutine Calc_Psi

    !!----
    !!---- Subroutine cell_fr_UB(ub,ipr,dcel,rcel)
    !!---- real(kind=cp),Dimension(3,3),         Intent(In) :: ub
    !!---- Integer, optional,                    Intent(In)  :: ipr
    !!---- real(kind=cp),Dimension(6), optional, Intent(out) :: dcel,rcel
    !!----
    !!----    Calculate and print cell parameters from UB-matrix
    !!----
    !!---- Update: May 2011
    !!
    Subroutine cell_fr_UB(ub,ipr,dcel,rcel)
       !---- Arguments ----!
       real(kind=cp),Dimension(3,3),         Intent(In)  :: ub
       Integer, optional,                    Intent(In)  :: ipr
       real(kind=cp),Dimension(6), optional, Intent(out) :: dcel,rcel

       !--- Local Variables ---!
       real(kind=cp), Dimension(3,3) :: g,sg,ubinv,angle
       real(kind=cp), Dimension(3)   :: acal,angcal,cala,calang
       integer                       :: i,j,k,jn,kn
       real(kind=cp)                 :: x

       g=Matmul(Transpose(ub),ub)
       sg=g
       !..inverse matrix g=b'b  to get direct cell parameters
       !..On first run through loop calculation of reciprocal cell, second do real cell
       Do k=1,2
           If(k==2) sg=invert(sg)
           g=sg
           Do  i=1,3
             acal(i)=Sqrt(g(i,i))
           End Do
           Do  i=1,3
             j=i
             jn=Mod(j,3)+1
             kn=Mod(jn,3)+1
             angcal(i)=acosd(g(jn,kn)/(acal(jn)*acal(kn)))
           End Do
           If(k==2) Exit
           cala   = acal   !store reciprocal latice in the first pass
           calang = angcal
       End Do
       if(present(dcel)) then
         dcel(1:3)=acal
         dcel(4:6)=angcal
       end if
       if(present(rcel)) then
         rcel(1:3)=cala
         rcel(4:6)=calang
       end if
       if(present(ipr)) then
         !.....Now invert UB to obtain the hkl's along the orthogonal diffractometer axes
         ubinv=invert(ub)
         Write(Unit=ipr,Fmt="(/,a)")                " => Parameters deduced from the UB matrix "
         Write(Unit=ipr,Fmt="(a,3f12.5,tr4,3f9.4)") " => Direct     cell dimensions: ",acal,angcal
         Write(Unit=ipr,Fmt="(a,3f12.8,tr4,3f9.4)") " => Reciprocal cell dimensions: ",cala,calang
         Write(Unit=ipr,Fmt="(/,a,/)") " =>  Inverse of UB-Matrix "
         Write(Unit=ipr,Fmt="(a)")            "             X(PH=0,CH=0)  Y(PH=90,CH=0)    Z(CHI=90)"
         Write(Unit=ipr,Fmt="(a,3f14.8)")     "          H",ubinv(1,:)
         Write(Unit=ipr,Fmt="(a,3f14.8)")     "          K",ubinv(2,:)
         Write(Unit=ipr,Fmt="(a,3f14.8)")     "          L",ubinv(3,:)
         !.....Now calculate angles between recip axes and orthogonal diffract. axes
         Do  i=1,3
           Do  j=1,3
             x = ub(i,j)/cala(j)
             If (x > 1.0) Then
                Write (ipr,"(a,3e12.4)") " Error x >1.0! Values of x,ub,acal: ",x,ub(i,j),acal(i)
             Else
               angle(i,j) = acosd(x)
             End If
           End Do
         End Do
         Write(Unit=ipr,Fmt="(/,a)")"    With all diffractometer angles set to 0, the angles between the "
         Write(Unit=ipr,Fmt="(a)")  "    reciprocal (A*,B*,C*) and diffractometer (X,Y,Z) axes are ..."
         Write(Unit=ipr,Fmt="(a,3f12.4,a)") "    (A*-X  B*-X  C*-X)  (",angle(1,:),")"
         Write(Unit=ipr,Fmt="(a,3f12.4,a)") "    (A*-Y  B*-Y  C*-Y)  (",angle(2,:),")"
         Write(Unit=ipr,Fmt="(a,3f12.4,a)") "    (A*-Z  B*-Z  C*-Z)  (",angle(3,:),")"
       end if
       Return
    End Subroutine cell_fr_UB

    !!----
    !!---- Subroutine Chi_mat(chi,dum)
    !!----    real(kind=cp), Intent(In)                   :: chi
    !!----    real(kind=cp), Intent(Out), Dimension(3,3)  :: dum
    !!----
    !!----    Calculate the Busing and Levy conventional rotation matrix for CHI (eq. 9)
    !!----    The CHI angle must be provided in degrees.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Chi_mat(chi,dum)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                   :: chi
       real(kind=cp), Intent(Out), Dimension(3,3)  :: dum

       dum=0.0
       dum(1,1)=cosd(chi)
       dum(1,3)=sind(chi)
       dum(2,2)=1.0
       dum(3,1)=-dum(1,3)
       dum(3,3)= dum(1,1)

       Return
    End Subroutine Chi_mat

    !!----
    !!---- Subroutine psd_convert(mpsd,gamm,gamp,nup,xobs,zobs,cath,anod,ierr)
    !!----    Integer, Intent(In)            :: mpsd
    !!----    real(kind=cp), Intent(In)      :: gamm
    !!----    real(kind=cp), Intent(In Out)  :: gamp
    !!----    real(kind=cp), Intent(In Out)  :: nup
    !!----    real(kind=cp), Intent(Out)     :: xobs
    !!----    real(kind=cp), Intent(Out)     :: zobs
    !!----    real(kind=cp), Intent(in Out)  :: cath
    !!----    real(kind=cp), Intent(in Out)  :: anod
    !!----    Integer, Intent(Out)           :: ierr
    !!----
    !!----    Original name: d19amd, now generalized to a series of PSDs
    !!----    Subroutine for getting Gamma and Nu of a reflections spot (GamP,NuP), given the
    !!----    gamma angle of the detector (GamM) and the pixel values (cath,anod). This is
    !!----    calculated when mpsd > 0, otherwise the inverse calculation is done. In both
    !!----    cases the detector coordinates (xobs,zobs) in mm are also calculated.
    !!----    The caracteristics of the detector are accessed via de global variable PSD of
    !!----    Type(Psd_Val_Type), that should be set by the calling program.
    !!----
    !!--<<    Original Comment:
    !!----    Specifically for D19A bannana detector, 4 x 64 degrees - 16 x 512
    !!----    cells and vertically curved.             R.F.D. STANSFIELD SEP-83
    !!----
    !!----    Modified for general case GamM .NE. GamP                     Feb-84
    !!----
    !!----    MPSD +VE - Find GamP, NuP given GamM, Cath and Anod
    !!----    MPSD -VE - Find Cath, Anod given GamM, GamP and NuP
    !!-->>
    !!----
    !!----    Extended for D19 cylindrical banana (from MJ Turner, peak find ...)
    !!----
    !!---- Update: July 2010
    !!
    Subroutine psd_convert(mpsd,gamm,gamp,nup,xobs,zobs,cath,anod,ierr)
       !---- Arguments ----!
       Integer, Intent(In)               :: mpsd
       real(kind=cp),    Intent(In)      :: gamm
       real(kind=cp),    Intent(In Out)  :: gamp
       real(kind=cp),    Intent(In Out)  :: nup
       real(kind=cp),    Intent(Out)     :: xobs
       real(kind=cp),    Intent(Out)     :: zobs
       real(kind=cp),    Intent(in Out)  :: cath
       real(kind=cp),    Intent(in Out)  :: anod
       Integer, Intent(Out)              :: ierr

       !--- Local Variables ---!
       real(kind=cp) :: cmid, amid, delga, td, tn, e, f, g, cd, a, b, z, d

       cmid=real(psd%ncat-1)/2.0
       amid=real(psd%nano-1)/2.0
       ierr = 0

       If (mpsd < 0) Then   !Find Cath, Anod given GamM, GamP and NuP
          delga=gamp - gamm
          td=tand(delga)
          tn=tand(nup)

          Select Case(psd%ipsd)

            Case(1)                               ! Vertically curved banana detector
             e=Sqrt(1.0 + td*td)
             f=psd%yoff*tn*e - psd%zoff
             g=psd%radius*Sqrt(1.0 + tn*tn + tn*tn*td*td)
             zobs=psd%radius*( Atan(tn*e) + Asin(f/g) )
             xobs=td*( psd%radius*Cos(zobs/psd%radius) + psd%yoff ) - psd%xoff
             anod=amid - zobs/psd%agap                          ! D19A
             cath=cmid - xobs/psd%cgap                          ! D19A

            Case(2)                                             ! Flat detector
             cd=cosd(delga)                                     ! D9, D10, Db21
             xobs=(psd%radius+psd%yoff)*td    - psd%xoff        ! D16, etc.
             zobs=(psd%radius+psd%yoff)*tn/cd - psd%zoff
             anod=amid + zobs/psd%agap
             cath=cmid - xobs/psd%cgap

            Case(3)                         ! Now checked!
             !xobs=(psd%radius + psd%yoff)*delga*to_rad - psd%xoff    ! D19-like Horizontal banana
             !zobs=tn*(psd%radius + psd%yoff) - psd%zoff
             e=sqrt(1 + tn*tn)
             f=psd%yoff*td*e - psd%xoff
             g=psd%radius*sqrt(1.0 + td*td + td*td*tn*tn)
             xobs=psd%radius*( atan(td*e) + asin(f/g) )
             zobs=tn*( psd%radius*cos(xobs/psd%radius) + psd%yoff ) - psd%zoff
             anod=amid - zobs/psd%agap
             cath=cmid - xobs/psd%cgap

          End Select

          If (anod > 0.0 .and. anod < real(psd%nano-1) .AND.  &
              cath > 0.0 .and. cath < real(psd%ncat-1)) Return
             ierr=-1

       Else   ! Find GamP, NuP given GamM, Cath and Anod

          Select Case(psd%ipsd)

            Case(1)                                  ! Vertically curved detector
             xobs=(cmid-cath)*psd%cgap               ! D19A
             zobs=(amid-anod)*psd%agap               ! D19A
             a=xobs                            + psd%xoff
             b=psd%radius*Cos(zobs/psd%radius) + psd%yoff
             z=psd%radius*Sin(zobs/psd%radius) + psd%zoff
             d=Sqrt(a*a + b*b)
             gamp=gamm + atan2d(a,b)
             nup=atan2d(z,d)

            Case(2)                            ! Flat detector
             xobs=  (cmid-cath)*psd%cgap       ! D9, D10, Db21
             zobs= -(amid-anod)*psd%agap       ! D16, etc.
             a=xobs       + psd%xoff
             b=psd%radius + psd%yoff
             z=zobs       + psd%zoff
             d=Sqrt(a*a + b*b)
             gamp=gamm + atan2d(a,b)
             nup=atan2d(z,d)

            Case(3)

              xobs=(cmid-cath)*psd%cgap     ! D19-like Horizontal banana (like in retreat)
              zobs=(amid-anod)*psd%agap     ! Need to check yoff at large xobs
              a= psd%radius*sin(xobs/psd%radius) + psd%xoff
              b= psd%radius*cos(xobs/psd%radius) + psd%yoff
              z=    zobs   + psd%zoff
              gamp=gamm + Atan2d(a,b)
              d=sqrt(a*a + b*b)
              nup=Atan2d(z,d)

          End Select
       End If

       Return
    End Subroutine psd_convert

    !!----
    !!---- Subroutine d19psd(mpsd,ga,nu,cath,anod,ierr)
    !!----    Integer, Intent(In Out)            :: mpsd
    !!----    real(kind=cp), Intent(In Out)      :: ga
    !!----    real(kind=cp), Intent(In Out)      :: nu
    !!----    real(kind=cp), Intent(In Out)      :: cath
    !!----    real(kind=cp), Intent(In Out)      :: anod
    !!----    Integer, Intent(In Out)            :: ierr
    !!----
    !!--<<    Original Comment:
    !!----    Specifically for D19A bannana detector, 4 X 64 degrees - 16 X 512
    !!----    cells and vertically curved.                RFD STANSFIELD SEP-83
    !!----
    !!----    MPSD +VE - Calculate delta GAMMA and NU from the cathode, anode
    !!----               co-ordinates
    !!----    MPSD -VE - Calculate the anode co-ordinate from NU
    !!----
    !!----    Some of the variables making reference to the characteristics of the
    !!-->>    detector are provisionally stored in a public type(Psd_Val_Type):: PSD
    !!----
    !!---- Update: April 2008
    !!
    Subroutine d19psd(mpsd,ga,nu,cath,anod,ierr)
       !---- Arguments ----!
       Integer, Intent(In)                :: mpsd
       real(kind=cp), Intent(In Out)      :: ga
       real(kind=cp), Intent(In Out)      :: nu
       real(kind=cp), Intent(In Out)      :: cath
       real(kind=cp), Intent(In Out)      :: anod
       Integer, Intent(In Out)            :: ierr

       !--- Local Variables ---!
       real(kind=cp) :: nurad,cmid,amid
       real(kind=cp) :: xo,yo,a,b,z,d

       cmid=real(psd%ncat-1)/2.0
       amid=real(psd%nano-1)/2.0

       If (mpsd >= 0) Then
          xo=(cmid-cath)*psd%cgap
          yo=(amid-anod)*psd%agap
          a=xo                            + psd%xoff
          b=psd%radius*Cos(yo/psd%radius) + psd%yoff
          z=psd%radius*Sin(yo/psd%radius) + psd%zoff
          d=Sqrt(a*a + b*b)
          ga=ga + atan2d(a,b)
          nu=atan2d(z,d)
       Else
          xo=-psd%xoff
          cath=cmid-xo/psd%cgap

          If (cath > 0.0 .AND. cath < real(psd%ncat-1)) Then
             nurad=nu*To_Rad
             yo=psd%radius*(nurad - Asin((psd%zoff*cosd(nu)-psd%yoff*sind(nu))/psd%radius))
             anod=amid-yo/psd%agap
             If ( .Not. (anod > 0.0 .AND. anod < real(psd%nano-1)) ) Then
                ierr=-1
                xo=0.0
                yo=0.0
             End If
          Else
             ierr=-1
             xo=0.0
             yo=0.0
          End If
       End If

       Return
    End Subroutine d19psd
    !!----
    !!---- Subroutine dspace(wave,vhkl,cell,ds,th,ierr)
    !!----    real(kind=cp), Intent(In)               :: wave
    !!----    real(kind=cp), Intent(In),Dimension(3)  :: vhkl
    !!----    real(kind=cp), Intent(In),Dimension(6)  :: cell
    !!----    real(kind=cp), Intent(Out)              :: ds,th
    !!----    Integer, Intent(Out)                    :: ierr
    !!----
    !!----    Calculate d-spacing and theta from cell parameters
    !!----    and wavelength assume triclinic symmetry. The reflection
    !!----    vector vhkl is provided in reciprocal lattice components
    !!----
    !!---- Update: April 2008
    !!
    Subroutine dspace(wave,vhkl,cell,ds,th,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In)               :: wave
       real(kind=cp), Intent(In),Dimension(3)  :: vhkl
       real(kind=cp), Intent(In),Dimension(6)  :: cell
       real(kind=cp), Intent(Out)              :: ds,th
       Integer, Intent(Out)                    :: ierr

       !--- Local Variables ---!
       real(kind=cp) :: a,b,c,al,be,ga,d,e,f,g, sint

       a=vhkl(1)/cell(1)
       b=vhkl(2)/cell(2)
       c=vhkl(3)/cell(3)
       al=cosd(cell(4))
       be=cosd(cell(5))
       ga=cosd(cell(6))
       d=a*(1.0-al*al) + b*(al*be-ga) + c*(ga*al-be)
       e=b*(1.0-be*be) + c*(be*ga-al) + a*(al*be-ga)
       f=c*(1.0-ga*ga) + a*(ga*al-be) + b*(be*ga-al)
       g=1.0 - al*al - be*be - ga*ga + 2.0*al*be*ga
       If (g /= 0.0) Then
          ds=(a*d + b*e + c*f)/g
          If (ds > 0.0) Then
             ds=1.0/Sqrt(ds)
             sint=wave/(2.0*ds)
             If (Abs(sint) <= 1.0)Then
                th=asind(sint)
             Else
                ierr=2
                th=0.0
             End If
          Else
             ierr=1
             th=0.0
          End If
       Else
          ierr=-1
          ds=0.0
          th=0.0
       End If

       Return
    End Subroutine dspace

    !!----
    !!---- Subroutine Equatorial_Chi_Phi(z1,ch,ph)
    !!----    real(kind=cp), Intent(In), Dimension(3)  :: z1
    !!----    real(kind=cp), Intent(Out)               :: ch,ph
    !!----
    !!----    Calculate CHI, PHI to put the vector Z1 in the equatorial plane
    !!----    The reciprocal vector Z1 is given in cartesian components with
    !!----    respect to the laboratory system (Eq. 38 of Busing & Levy)
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Equatorial_Chi_Phi(z1,ch,ph)
       !---- Arguments ----!
       real(kind=cp), Intent(In), Dimension(3)  :: z1
       real(kind=cp), Intent(Out)               :: ch,ph

       !---- Local Variables ----!
       real(kind=cp) :: d

       If (z1(1) /= 0.0 .or. z1(2) /= 0.0)  Then
          ph=atan2d(z1(2),z1(1))
          d=Sqrt(z1(1)*z1(1)+z1(2)*z1(2))
          ch=atan2d(z1(3),d)
       Else
          ph=0.0
          ch=90.0
          If(z1(3) < 0.0) ch=-ch
       End If

       Return
    End Subroutine Equatorial_Chi_Phi

    !!----
    !!---- Subroutine fixdnu(wave,z1,ch,ph,ga,om,nu,ierr)
    !!----    real(kind=cp), Intent(In)               :: wave
    !!----    real(kind=cp), Intent(In), Dimension(3) :: z1
    !!----    real(kind=cp), Intent(In)               :: nu
    !!----    real(kind=cp), Intent(Out)              :: ch,ph,ga,om
    !!----    Integer, Intent(Out)                    :: ierr
    !!----
    !!----    Calculate a setting CH,PH,GA,OM to put the diffracted beam at NU.
    !!----    PH puts the diffraction vector Z1 into the CHI circle (as for
    !!----    bisecting geometry), CH brings the vector to the appropriate NU
    !!----    and OM then positions the beam at GA.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine fixdnu(wave,z1,nu,ch,ph,ga,om,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In)               :: wave
       real(kind=cp), Intent(In), Dimension(3) :: z1
       real(kind=cp), Intent(In)               :: nu
       real(kind=cp), Intent(Out)              :: ch,ph,ga,om
       Integer, Intent(Out)                    :: ierr

       !--- Local Variables ---!
       real(kind=cp), Dimension(3) ::  z4
       real(kind=cp) :: ch1,ch2,twoth,theta, cosga

       ierr=0
       !---- Get CH1 and PH as for bisecting geometry
       Call Angs_4c_bisecting(wave,z1,twoth,theta,ch1,ph,ierr)
       If (ierr == 0) Then
          If (Abs(cosd(nu)) <= 0.0001) Then
             !---- One unique vector can diffract at ABS(NU)=90. (for which THETA=45)
             If (.Not. (theta > 44.99 .and. theta < 45.01) )  Then
                ierr=-1
             Else
                ga=90.0
                om=90.0
                ch2=Sign(1.0_cp,nu)*45.0_cp
                ch=ch1-ch2
                !---- Convert to -180.0 <= CH <= +180.0
                ch=ch-360.0_cp*Int((Sign(180.0_cp,ch)+ch)/360.0_cp)
             End If
          Else
             cosga=cosd(twoth)/cosd(nu)  !General expression coming from the definition of the angles
             If (Abs(cosga) > 1.0_cp) Then
                ierr=-2
             Else
                ga=acosd(cosga)
                !- --- Diffraction vector in lab system
                Call z4frgn(wave,ga,nu,z4)
                om=atan2d(-z4(2),z4(1))
                !---- CH2 is the angle between the scattering vector and horizontal plane
                ch2=asind(z4(3)*wave/(2.0*sind(theta)))
                ch=ch1-ch2
                !---- Convert to -180.0<= CH <=+180.0
                ch=ch-360.0_cp*Int((Sign(180.0_cp,ch)+ch)/360.0_cp)
             End if
          End If
       End If

       If (ierr /= 0) Then
          ch=0.0_cp
          ph=0.0_cp
          ga=0.0_cp
          om=0.0_cp
       End If

       Return
    End Subroutine fixdnu
    !!----
    !!---- Subroutine Flat_Cone_vertDet(wave,z1,ub,vrho,rho,ch,ph,ga,om,nu,ierr)
    !!----    real(kind=cp), Intent(In)                   :: wave
    !!----    real(kind=cp), Intent(In),  Dimension(3)    :: z1
    !!----    real(kind=cp), Intent(In),Dimension(3,3)    :: ub
    !!----    real(kind=cp), Intent(In Out), Dimension(3) :: vrho
    !!----    real(kind=cp), Intent(Out)                  :: rho,ch,ph,ga,om,nu
    !!----    Integer, Intent(Out)                        :: ierr
    !!----
    !!----    Calculate RHO (=PSI) about given rotation vector VRHO (and the
    !!----    corresponding angles OM, CH, PH) to put the vector Z1 into the
    !!----    flat-cone diffracting position. Calculate resulting GAMMA and NU.
    !!----
    !!--<<    V = [UB].H = [UB].[G].D where D is a real-space vector.
    !!-->>        [G]-1 = [B]T.[B] = [UB]T.[UB] hence V = [[UB]T]-1.D
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Flat_Cone_vertDet(wave,z1,ub,vrho,rho,ch,ph,ga,om,nu,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                      :: wave
       real(kind=cp), Intent(In),     Dimension(3)    :: z1
       real(kind=cp), Intent(In),   Dimension(3,3)    :: ub
       real(kind=cp), Intent(In Out), Dimension(3)    :: vrho
       real(kind=cp), Intent(Out)                     :: rho,ch,ph,ga,om,nu
       Integer, Intent(Out)                           :: ierr

       !---- Local Variables ----!
       real(kind=cp), Dimension(3)   ::  z3, z4, v1
       real(kind=cp), Dimension(3,3) ::  dum1
       real(kind=cp)                 ::  pstar, theta, duf, sinmu, cosnu, psi1,csi1,csi2

       dum1=Transpose(ub)
       dum1=invert(dum1)
       v1=Matmul(dum1,vrho)  ! Components of vrho in reciprocal lattice
       Call normal(v1,ierr)
       If (ierr == 0) Then
          !---- Pstar is the reciprocal spacing of this layer
          pstar= Dot_Product(z1,v1)
          If (pstar < 0.0_cp) Then
             v1=-v1
             vrho=-vrho
             pstar=-pstar
          End If

          !---- Get GA (=MU) for this flat-cone, then NU.
          !---- Choose 0.0_cp <= GA <=180.0_cp and 0.0_cp <= NU <=90.0_cp  (GA,NU = 180.0_cp+GA,180.0_cp-NU)
          Call Get_dspacing_theta(wave,z1,duf,theta,ierr)
          If (ierr == 0) Then
             sinmu=wave*pstar
             If (sinmu > 1.0_cp) Then
                ierr=-1
             Else
                ga=asind(sinmu)
                !---- TH > 45.0_cp is equivalent to PSTAR > (SQRT(2.0_cp))/WAVE
                If (Abs(theta) > 45.0_cp) ga=180.0_cp-ga
                !---- TH = 45.0 is equivalent to MU=90.0, and NU is indeterminate.
                If (sinmu >= 1.0_cp) Then
                   nu=0.0_cp
                Else
                   cosnu=cosd(2.0_cp*theta)/cosd(ga)
                   If (Abs(cosnu) > 1.0_cp) Then
                      ierr=-2
                   Else
                      nu=acosd(cosnu)
                   End If
                End If

                If (ierr == 0) Then
                   !---- Get CH,PH to put rec.-lat.-plane-normal V1 in equatorial plane
                   !---- OM=GA puts V1 in flat-cone setting
                   Call Equatorial_Chi_Phi(v1,ch,ph)
                   om=ga
                   !---- Get RHO, the PSI angle about V1 to make Z1 diffract!
                   !---- The present OM,CH,PH define PSI1 about V1.
                   Call Calc_Psi(vrho,v1,om,ch,ph,ub,psi1,ierr)
                   If (ierr == 0) Then
                      !---- The angle (-CSI1+CSI2), about V1 (coincident with the x-omega axis),
                      !---- puts the present Z3 into the diffracting position.
                      Call z3frz1(z1,ch,ph,z3)
                      csi1=atan2d(z3(3),z3(2))
                      Call z4frgn(wave,ga,nu,z4)
                      Call z1frz2(z4,om,z3)
                      csi2=atan2d(z3(3),z3(2))
                      rho=psi1-csi1+csi2
                      rho=rho-360.0_cp*Int((Sign(180.0_cp,rho)+rho)/360.0_cp)
                      Call Calc_Om_Chi_Phi(vrho,v1,rho,ub,om,ch,ph,ierr)
                      If (ierr == 0) Return
                   End If
                End If !ierr==0
             End If
          End If
       End If
       If (ierr /= 0) Then
          rho=0.0_cp
          ch=0.0_cp
          ph=0.0_cp
          ga=0.0_cp
          om=0.0_cp
          nu=0.0_cp
       End If

       Return
    End Subroutine Flat_Cone_vertDet

    !!----
    !!---- Subroutine genb(c,b)
    !!----    Type(Crystal_Cell_Type), Intent(In)   :: c
    !!----    real(kind=cp), Intent(Out), Dimension(3,3)     :: b
    !!----
    !!--<<   Calculation of [B] matrix
    !!----   BUSING&LEVY ACTA CRYST.(1967)22,457-464  EQUATION 3
    !!----   WOOSTER R. SCI. INSTRUM. (1962)39,103
    !!----      C%cell  : Direct cell parameters
    !!----      C%rcell : Reciprocal cell parameters
    !!----      C%ang   : Direct cell angles
    !!-->>      C%rcell : Reciprocal cell angles
    !!----
    !!---- Update: April 2008
    !!
    Subroutine genb(c,b)
       !---- Arguments ----!
       Type(Crystal_Cell_Type),      Intent(In)   :: c
       real(kind=cp), Dimension(3,3),Intent(Out)  :: b

       b(1,1)=c%rcell(1)
       b(1,2)=c%rcell(2)*cosd(c%rang(3))
       b(1,3)=c%rcell(3)*cosd(c%rang(2))
       b(2,2)=c%rcell(2)*sind(c%rang(3))
       b(2,3)=-(c%rcell(3)*sind(c%rang(2))*cosd(c%ang(1)))
       b(3,3)=1.0_cp/c%cell(3)
       b(2,1)=0.0_cp
       b(3,1)=0.0_cp
       b(3,2)=0.0_cp

       Return
    End Subroutine genb

    !!----
    !!---- Subroutine GenUB(b,h1,h2,h1o,h2o,ub, ierr)
    !!----    real(kind=cp), Dimension(3,3), Intent(In)  :: B        !Busing-Levy B-matrix
    !!----    real(kind=cp), Dimension(3),   Intent(In)  :: h1,h2    !Miller indices
    !!----    real(kind=cp), Dimension(3),   Intent(In)  :: h1o,h2o  !Components in Lab system
    !!----    real(kind=cp), Dimension(3,3), Intent(Out) :: ub
    !!----
    !!----    Original from   A.Filhol  25/05/84
    !!----    Given the [B] matrix, the Miller indices of two reflections, h1 & h2,
    !!----    and the components of these two reflections, h1o & h2o, in the laboratory
    !!----    system, this subroutine provides the matrix UB. Only the direction in the
    !!----    laboratory system of reflections are needed, e.g. h1o and h2o may be unitary
    !!----    vectors or whatever other vector along these directions.
    !!----
    !!----    [hc] : Reflection H,K,L in the reciprocal lattice orthonormal system
    !!----    [hc] = [B] [h]
    !!----    [ho] : Reflection H,K,L in the Cartesian laboratory system
    !!----    [ho]=[UB]*[hc]
    !!----
    !!---- Update: April 2008
    !!
    Subroutine GenUB(B,h1,h2,h1o,h2o,UB, ierr)
       !---- Arguments ----!
       real(kind=cp), Dimension(3,3), Intent(In)  :: B        !Busing-Levy B-matrix
       real(kind=cp), Dimension(3),   Intent(In)  :: h1,h2    !Miller indices
       real(kind=cp), Dimension(3),   Intent(In)  :: h1o,h2o  !Components in Lab system
       real(kind=cp), Dimension(3,3), Intent(Out) :: UB       !Busing-Levy UB-matrix
       Integer,                       Intent(Out) :: ierr

       !--- Local Variables ---!
       real(kind=cp), Dimension(3)  :: h1c,h2c,v1,v2
       real(kind=cp), Dimension(3,3):: trpc,trpo,u

       ierr=0
       !Calculation of the reciprocal components in the cartesian reciprocal axes
       h1c=Matmul(B,h1)
       h2c=Matmul(B,h2)
       v1=h1o
       v2=h2o
       Call triple(h1c,h2c,trpc,ierr) !Orthonormal frame attached to h1c,h2c
       If (ierr /= 0) Then
          ub=0.0_cp
          Return
       End If
       Call triple(v1,v2,trpo,ierr) !Orthonormal frame attached to h1o,h2o
       If (ierr /= 0) Then
          ub=0.0_cp
          Return
       End If
       !..... MATRIX [U]  *** Equation 27 (B-L)
       U=Matmul(trpo,Transpose(trpc))

       !..... MATRIX [UB]
       UB=Matmul(U,B)

       Return
    End Subroutine GenUB

    !!----
    !!---- Subroutine Get_Angs_NB(wave,z1,ga,om,nu,ierr)
    !!----    real(kind=cp), Intent(In)                      :: wave
    !!----    real(kind=cp), Intent(In),Dimension(3)         :: z1
    !!----    real(kind=cp), Intent(Out)                     :: ga,om,nu
    !!----    Integer, Intent(Out)                           :: ierr
    !!----
    !!----    Calculate normal-beam angles GAMMA, OMEGA, NU to put the
    !!----    vector Z1 into the diffracting condition.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Get_Angs_NB(wave,z1,ga,om,nu,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                      :: wave
       real(kind=cp), Intent(In),Dimension(3)         :: z1
       real(kind=cp), Intent(Out)                     :: ga,om,nu
       Integer, Intent(Out)                           :: ierr

       !--- Local Variables ---!
       real(kind=cp), Dimension(3)   ::  znew
       real(kind=cp), Dimension(3,3) ::  dum
       real(kind=cp)                 :: theta,d,a,sint,b

       Call Get_dspacing_theta(wave,z1,d,theta,ierr)
       ga=0.0_cp
       om=0.0_cp
       nu=0.0_cp
       If (ierr == 0) Then
          a=Sqrt(z1(1)*z1(1)+z1(2)*z1(2))
          If (a <= 1.0e-10_cp) Then
             !---- Anything on the omega axis is blind
             ierr=-1
          Else
             sint=sind(theta)
             b=2.0_cp*sint*sint/(wave*a)
             If (b > 1.0_cp) Then
                ierr=-2
             Else
                a=-atan2d(z1(2),-z1(1))
                b=-asind(b)
                om=a+b
                Call Phi_mat(om,dum)
                znew=Matmul(dum,z1)
                If (znew(1) <= 0.0_cp) Then
                   om=om-2.0_cp*atan2d(-znew(1),-znew(2))
                End If
                om=om-360.0_cp*Int((Sign(180.0_cp,om)+om)/360.0_cp)
                nu=asind(wave*z1(3))
                ga=acosd(cosd(2.0_cp*theta)/cosd(nu))
             End If
          End If
       End If

       Return
    End Subroutine Get_Angs_NB

    !!----
    !!---- Subroutine Get_dspacing_theta(wave,z1,ds,th,ierr)
    !!----    real(kind=cp), Intent(In)                 :: wave
    !!----    real(kind=cp), Intent(In), Dimension(3)   :: z1
    !!----    real(kind=cp), Intent(Out)                :: ds,th
    !!----    Integer, Intent(Out)                      :: ierr
    !!----
    !!----    Calculate D-spacing (real space) and THETA from the length of Z1
    !!----    The reciprocal vector Z1 is given in cartesian components with
    !!----    respect to the laboratory system. If ierr=1 the calculated d-spacing
    !!----    is is fixed to 0.0 as well as theta. This error condition appears when
    !!----    the length of the reciprocal vector z1 is lower or equal to 0.0001
    !!----    If ierr=2 the reflection is outside the resolution sphere.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Get_dspacing_theta(wave,z1,ds,th,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                 :: wave
       real(kind=cp), Intent(In), Dimension(3)   :: z1
       real(kind=cp), Intent(Out)                :: ds,th
       Integer, Intent(Out)                      :: ierr

       !---- Local Variables ---!
       real(kind=cp) :: dstar, sint

       ierr=0
       dstar=Sqrt(Dot_Product(z1,z1))
       If (dstar > 0.0001_cp) Then
          ds=1.0/dstar
          sint=wave*dstar/2.0_cp
          If (Abs(sint) <= 1.0_cp) Then
             th=asind(sint)
          Else
             ierr=2
             th=0.0_cp
          End If
       Else
          ierr=1
          ds=0.0_cp
          th=0.0_cp
       End If

       Return
    End Subroutine Get_dspacing_theta


    !!---- Subroutine Get_FlatCone_Angles_D10(z,mu_fc,psi1,psi2,n_ang,limits,psi,omega,chi,phi,inlim,Rot_base)
    !!----   real(kind=cp), dimension(3),             intent(in)  :: z
    !!----   real(kind=cp),                           intent(in)  :: mu_fc,psi1,psi2
    !!----   integer,                                 intent(in)  :: n_ang
    !!----   real(kind=cp), dimension(:,:),           intent(in)  :: limits
    !!----   real(kind=cp), dimension(:),             intent(out) :: psi,omega,chi,phi
    !!----   logical,       dimension(:),             intent(out) :: inlim
    !!----   real(kind=cp), dimension(3,3), optional, intent(out) :: Rot_base
    !!----
    !!----  Subroutine for getting all accessible angles for making a flat-cone data collection with
    !!----  a geometry like D10
    !!----
    !!----
    !!----
    !!----  Created: March 2013 (JRC)
    !!----
    Subroutine Get_FlatCone_Angles_D10(z,mu_fc,psi1,psi2,n_ang,limits,psi,omega,chi,phi,inlim,Rot_base)
      real(kind=cp), dimension(3),             intent(in)  :: z
      real(kind=cp),                           intent(in)  :: mu_fc,psi1,psi2
      integer,                                 intent(in)  :: n_ang
      real(kind=cp), dimension(:,:),           intent(in)  :: limits
      real(kind=cp), dimension(:),             intent(out) :: psi,omega,chi,phi
      logical,       dimension(:),             intent(out) :: inlim
      real(kind=cp), dimension(3,3), optional, intent(out) :: Rot_base

      !--- Local variables ---!
      real(kind=cp),dimension(3)    :: dL,angl
      real(kind=cp),dimension(3,3)  :: Rot,R_psi,Rt
      real(kind=cp)                 :: ps,om,ch,ph,oms,chs,phs,step
      integer                       :: i
      logical                       :: ok

      !Initializing all output arrays
      inlim=.false.; psi=0.0_cp; omega=0.0_cp; chi=0.0_cp; phi=0.0_cp

      step=(psi2-psi1)/real(n_ang-1,kind=cp)
      dL=(/0.0_cp, sind(mu_fc), cosd(mu_fc)/)    !unitary vector normal to detector plane in L-system for D10 z-down
      call Get_Matrix_moving_v_to_u(z,dL,Rt)     !One of the infinite matrices

      !Saving the base rotation matrix corresponding to psi=0.0 if needed
      if(present(Rot_base))   Rot_base=Rt

      ok=.false.
      do i=1,n_ang
        ps=psi1+real(i-1,kind=cp)*step
        R_psi= Rot_Matrix(dL,ps)
        Rot=Matmul(R_psi,Rt)   !Full rotation matrix

        !Get the Euler angles corresponding to matrix Rot
        Call Get_OmegaChiPhi(Rot,om,ch,ph,"D")
        if(ERR_Geom) then
            write(*,"(a)") trim(ERR_Geom_Mess)
        else
            angl=(/om,ch,ph/)
            if(in_limits(3,limits,angl)) then
              omega(i)=om; chi(i)=ch; phi(i)=ph; psi(i)=ps
              inlim(i)=.true.
              cycle
            end if
            oms=om+180.0_cp
            chs=-ch
            phs=ph+180_cp
            angl=Chkin180((/oms,chs,phs/))
            if(in_limits(3,limits,angl)) then
              omega(i)=oms; chi(i)=chs; phi(i)=phs; psi(i)=ps
              inlim(i)=.true.
              cycle
            end if
            omega(i)=om; chi(i)=ch; phi(i)=ph
        end if
      end do
      return
    End Subroutine Get_FlatCone_Angles_D10

    !!----
    !!---- Subroutine Get_GaOmNu_frChiPhi(wave,z1,ch,ph,ga,om,nu,ierr)
    !!----    real(kind=cp), Intent(In)                :: wave
    !!----    real(kind=cp), Intent(In), Dimension(3)  :: z1
    !!----    real(kind=cp), Intent(In)                :: ch,ph
    !!----    real(kind=cp), Intent(Out)               :: ga,om,nu
    !!----    Integer, Intent(Out)                     :: ierr
    !!----
    !!----    Given CHI & PHI, calculate normal-beam angles GAMMA, OMEGA, NU
    !!----    to put the vector Z1 into the diffraction condition.
    !!----    converts VLAB2=[CHI0].[PHI0].VLAB1 and relies on subroutine
    !!----    Get_Angs_NB to do the rest!
    !!----    The reciprocal vector Z1 is given in Cartesian components with
    !!----    respect to the laboratory system.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Get_GaOmNu_frChiPhi(wave,z1,ch,ph,ga,om,nu,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                :: wave
       real(kind=cp), Intent(In), Dimension(3)  :: z1
       real(kind=cp), Intent(In)                :: ch,ph
       real(kind=cp), Intent(Out)               :: ga,om,nu
       Integer, Intent(Out)                     :: ierr

       !--- Local Variables ---!
       real(kind=cp),Dimension(3) :: z3

       Call z3frz1(z1,ch,ph,z3)
       Call Get_Angs_NB(wave,z3,ga,om,nu,ierr)

       Return
    End Subroutine Get_GaOmNu_frChiPhi

    !!---- Subroutine Get_UB_from_hkl_hkl_omega(wave,Cell,h1,h2,omega,UB,ok,mess)
    !!----   real(kind=cp),                 intent(in)    :: wave  !Vavelength
    !!----   type (Crystal_Cell_Type),      intent(in)    :: Cell  !Unit cell object
    !!----   real(kind=cp), dimension(3),   intent(in)    :: h1    !Indices of the known first reflection in plane
    !!----   real(kind=cp), dimension(3),   intent(in)    :: h2    !Indices of the known second reflection in plane
    !!----   real(kind=cp),                 intent(in)    :: omega !Value of the omega motor for the second reflection (vertical spindle)
    !!----   real(kind=cp), dimension(3,3), intent(out)   :: UB    !Generated Busing-Levy UB-matrix
    !!----   logical,                       intent(out)   :: ok    !If .true. everything has gone well
    !!----   character(len=*),              intent(out)   :: mess  !Error message in case ok=.false.
    !!----
    !!----   This subroutine generates a UB matrix when two reflections in the horizontal plane
    !!----   are known (indices hkl and h'h'l') and the second reflection has been measured and
    !!----   its omega angle is known.
    !!----
    !!----   Updated: June-2012 (JRC)
    !!----
    Subroutine Get_UB_from_hkl_hkl_omega(wave,Cell,h1,h2,omega,UB,ok,mess)
      real(kind=cp),                 intent(in)    :: wave
      type (Crystal_Cell_Type),      intent(in)    :: Cell
      real(kind=cp), dimension(3),   intent(in)    :: h1
      real(kind=cp), dimension(3),   intent(in)    :: h2
      real(kind=cp),                 intent(in)    :: omega
      real(kind=cp), dimension(3,3), intent(out)   :: UB
      logical,                       intent(out)   :: ok
      character(len=*),              intent(out)   :: mess
      ! Local variables
      integer                       :: ierr
      real(kind=cp)                 :: theta1,theta2,alpha,del_omega,d1s,d2s
      real(kind=cp), dimension(3)   :: ho1,ho2,s1,s2
      real(kind=cp), dimension(3,3) :: Rot

      ok=.true.
      mess= " "
      if(Co_linear(h1,h2,3)) then
        ok=.false.
        mess="The two provided reflections are co-linear, no UB-matrix can be calculated"
        return
      end if
      !
      !Determination of the Bragg angles of two reflections in the scattering plane
      !
      d1s=sqrt(dot_product(h1,matmul(Cell%GR,h1))) !d1*
      d2s=sqrt(dot_product(h2,matmul(Cell%GR,h2))) !d2*
      theta1=asind(0.5*wave*d1s)
      theta2=asind(0.5*wave*d2s)
      alpha=acosd(dot_product(h1,matmul(Cell%GR,h2))/d1s/d2s) !Angle between the two reciprocal vectors
      del_omega=theta1-theta2+alpha  !Variation in omega to put the first reflection in diffraction position
      !
      !Calculation of the Cartesian components of the two reflections in the scattering plane
      !
      s2=d2s*(/cosd(Theta2),-sind(Theta2),0.0_cp/)   !z4   diffraction position
      call Phi_Mat(omega,Rot)
      ho2=matmul(transpose(rot),s2)               !z1   zero motor angles
      s1=d1s*(/cosd(Theta1),-sind(Theta1),0.0_cp/)   !z4   diffraction position
      call Phi_mat(omega+del_omega,Rot)
      ho1=matmul(transpose(rot),s1)               !z1   zero motor angles
      !
      ! Generate UB-matrix
      !
      call GenUB(Cell%BL_M,h1,h2,ho1,ho2,UB, ierr)
      if(ierr /= 0) then
        ok=.false.
        mess = "Error in the calculation of UB-matrix "
      end if
      return
    End Subroutine Get_UB_from_hkl_hkl_omega

    !!---- Subroutine Get_UB_from_uvw_hkl_omega(wave,Cell,Zone_Axis,h1,omega,UB,ok,mess)
    !!----   real(kind=cp),                 intent(in)    :: wave      !Vavelength
    !!----   type (Crystal_Cell_Type),      intent(in)    :: Cell      !Unit cell object
    !!----   Type (Zone_Axis_type),         intent(in out):: Zone_Axis !Zone axis (See CFML_Crystal_Metrics module)
    !!----   real(kind=cp), dimension(3),   intent(in)    :: h1        !Indices of the known reflection in plane
    !!----   real(kind=cp),                 intent(in)    :: omega     !Value of the omega motor (vertical spindle)
    !!----   real(kind=cp), dimension(3,3), intent(out)   :: UB        !Generated Busing-Levy UB-matrix
    !!----   logical,                       intent(out)   :: ok        !If .true. everything has gone well
    !!----   character(len=*),              intent(out)   :: mess      !Error message in case ok=.false.
    !!----
    !!----   This subroutine generates a UB matrix when the vertical zone axis of the crystal
    !!----   is known and a reflection in the horizonal plane has been measured with known
    !!----   indices an value of the omega angle.
    !!----
    !!----   Updated: June-2012 (JRC)
    !!----
    Subroutine Get_UB_from_uvw_hkl_omega(wave,Cell,Zone_Axis,h1,omega,UB,ok,mess)
      real(kind=cp),                 intent(in)    :: wave
      type (Crystal_Cell_Type),      intent(in)    :: Cell
      Type (Zone_Axis_type),         intent(in out):: Zone_Axis
      real(kind=cp), dimension(3),   intent(in)    :: h1
      real(kind=cp),                 intent(in)    :: omega
      real(kind=cp), dimension(3,3), intent(out)   :: UB
      logical,                       intent(out)   :: ok
      character(len=*),              intent(out)   :: mess
      ! Local variables
      integer                        :: ierr
      real(kind=cp)                  :: theta1,theta2,alpha,del_omega,d1s,d2s
      real(kind=cp), dimension(3)    :: h2,ho1,ho2,s1,s2
      real(kind=cp), dimension(3,3)  :: Rot

      ok=.true.
      mess= " "
      !First check that the provided reflection is perpendicular to uvw
      if(dot_product(real(Zone_Axis%uvw),h1) > 0.01_cp) then
        ok=.false.
        write(unit=mess,fmt="(2(a,f8.3))") "The given reflection: ",h1," is not perpendicular to: ",real(Zone_Axis%uvw)
        return
      end if
      call Get_basis_from_uvw(1.0_cp,Zone_Axis%uvw,Cell,zone_axis,ok) !Here we use a dmin=1.0 angstrom
      if(.not. ok) then
        mess = "Error in the calculation of reflection plane "
        return
      end if
      h2=real(zone_axis%rx)       !Second reflection in the plane
      if(Co_linear(h1,h2,3)) h2=real(zone_axis%ry)
      !
      !Determination of the Bragg angles of two reflections in the scattering plane
      !
      d1s=sqrt(dot_product(h1,matmul(Cell%GR,h1))) !d1*
      d2s=sqrt(dot_product(h2,matmul(Cell%GR,h2))) !d2*
      theta1=asind(0.5*wave*d1s)
      theta2=asind(0.5*wave*d2s)
      alpha=acosd(dot_product(h1,matmul(Cell%GR,h2))/d1s/d2s) !Angle between reciprocal vectors
      del_omega=theta2-theta1+alpha  !Variation in omega to put the second reflection in diffraction position
      !
      !Calculation of the Cartesian components of the two reflections in the scattering plane
      !
      s1=d1s*(/cosd(Theta1),-sind(Theta1),0.0_cp/)
      call Phi_Mat(omega,Rot)
      ho1=matmul(transpose(rot),s1)
      s2=d2s*(/cosd(Theta2),-sind(Theta2),0.0_cp/)
      call Phi_mat(omega+del_omega,Rot)
      ho2=matmul(transpose(rot),s2)
      !
      ! Generate UB-matrix
      !
      call GenUB(Cell%BL_M,h1,h2,ho1,ho2,UB, ierr)
      if(ierr /= 0) then
        ok=.false.
        mess = "Error in the calculation of UB-matrix "
      end if
      return
    End Subroutine Get_UB_from_uvw_hkl_omega

    !!----
    !!---- Subroutine Get_WaveGaNu_frZ4(z4,wave,ga,nu,ierr)
    !!----    real(kind=cp),    Intent(In), Dimension(3)  :: z4
    !!----    real(kind=cp),    Intent(Out)               :: wave,ga,nu
    !!----    Integer, Intent(Out)               :: ierr
    !!----
    !!----    Calculate GA, NU and wavelength for diffraction
    !!----    vector in laboratory system
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Get_WaveGaNu_frZ4(z4,wave,ga,nu,ierr)
       !---- Arguments ----!
       real(kind=cp),    Intent(In), Dimension(3)  :: z4
       real(kind=cp),    Intent(Out)               :: wave,ga,nu
       Integer, Intent(Out)                        :: ierr

       !--- Local Variables ---!
       real(kind=cp) :: a,b,sinnu, cosnu,cosga, singa

       a=-2.0_cp*z4(2)
       b=Dot_Product(z4,z4)
       If (.Not. (a <= 0.0_cp .or. b <= 0.0_cp) ) Then
          wave=a/b
          sinnu=wave*z4(3)
          If (Abs(sinnu) > 1.0_cp) Then
             ierr=1
             nu=0.0_cp
             ga=0.0_cp
          Else If(Abs(sinnu) < 1.0_cp) Then
             nu=asind(sinnu)
             cosnu=Sqrt(1.0_cp - sinnu*sinnu)
             singa=wave*z4(1)/cosnu
             If (Abs(singa) > 1.0_cp) Then
                ierr=1
                nu=0.0_cp
                ga=0.0_cp
                Return
             End If
             cosga=(wave*z4(2)+1.0_cp)/cosnu
             ga=atan2d(singa,cosga)
          Else
             nu=90.0_cp
             ga=0.0_cp
          End If
       Else
          wave=0.0_cp
          ierr=1
          nu=0.0_cp
          ga=0.0_cp
       End If

       Return
    End Subroutine Get_WaveGaNu_frZ4

    !!----
    !!---- Subroutine Get_z1_D9angls(wave,ttheta,om,ch,ph,z1)
    !!----    real(kind=cp),               intent (in)  :: wave,ttheta,om,ch,ph
    !!----    real(kind=cp), dimension(3), intent (out) :: z1
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Get_z1_D9angls(wave,ttheta,om,ch,ph,z1)
       !---- Arguments ----!
       real(kind=cp),               intent (in)  :: wave,ttheta,om,ch,ph
       real(kind=cp), dimension(3), intent (out) :: z1

       !--- Local Variables ---!
       real(kind=cp) :: th,dsp,ome,cosom,coski,cosfi,sinom,sinfi

       th=ttheta*0.5_cp
       dsp= 2.0_cp*sind(th)/wave
       ome=om-th
       cosom=cosd(ome)
       sinom=sind(ome)
       coski=cosd(ch)
       sinfi=sind(ph)
       cosfi=cosd(ph)
       z1(1)=cosom*coski*cosfi-sinom*sinfi
       z1(2)=cosom*coski*sinfi+sinom*cosfi
       z1(3)=cosom*sind(ch)
       z1=z1*dsp

       return
    End Subroutine Get_z1_D9angls

    !!----
    !!---- Subroutine Get_z1_from_pixel(npx,npz,ifr,snum,z1)
    !!----    integer,                intent(in)  :: npx,npz,ifr (pixel and frame)
    !!----    type(SXTAL_Numor_type), intent(in)  :: snum
    !!----    real(kind=cp), dimension(3),     intent(out) :: z1
    !!----
    !!---- Update: May 2011
    !!
    Subroutine Get_z1_from_pixel(npx,npz,ifr,snum,z1)
       !---- Arguments ----!
       integer,                      intent(in)  :: npx,npz,ifr
       type(SXTAL_Numor_type),       intent(in)  :: snum
       real(kind=cp), dimension(3),  intent(out) :: z1

       !---- Local Variables ----!
       integer        :: ier, mpsd
       real(kind=cp)  :: gamm,gamp,nup,xobs,zobs,cath,anod, wave,chim,phim,omem

       mpsd  = 1  !Find Gamma_Pixel and Nu_Pixel given GamM, Cath and Anod in PSD_Convert
      phim  = snum%angles(1)  !Angles corresponding to hmin,kmin,lmin
      chim  = snum%angles(2)
      omem  = snum%angles(3)
      gamm  = snum%angles(4)
      if(snum%scantype == 'omega') then
         omem  = snum%tmc_ang(4,ifr)
      else if(snum%scantype == 'phi') then
         phim  = snum%tmc_ang(4,ifr)
      else if(snum%scantype == 'q-scan') then  ! Angles: gamma, omega, Chi,phi, psi?
         gamm  = snum%tmc_ang(4,ifr)
         omem  = snum%tmc_ang(5,ifr)
         chim  = snum%tmc_ang(6,ifr)
         phim  = snum%tmc_ang(7,ifr)
      end if

       anod  = npx
       cath  = npz
       wave  = Current_Orient%wave

       ! Find GAMP and NUP for this pixel
       Call Psd_Convert(mpsd,gamm,gamp,nup,xobs,zobs,cath,anod,ier)

       ! Find the scattering vector in cartesian coordinates for this pixel
       ! from GAMP, NUP, CHIM, OMEGM and PHIM
       Call z1frmd(wave,chim,phim,gamp,omem,nup,z1)

       return
    End Subroutine Get_z1_from_pixel

    !!----
    !!---- Subroutine normal(v,ierr)
    !!----    real(kind=cp), Intent(In Out), Dimension(3)   :: v
    !!----    Integer, Intent(Out)                          :: ierr
    !!----
    !!----    Normalise vector V (in Cartesian components)
    !!----
    !!---- Update: April 2008
    !!
    Subroutine normal(v,ierr)
       !---- Argument ----!
       real(kind=cp), Intent(In Out), Dimension(3)   :: v
       Integer, Intent(Out)                          :: ierr

       !--- Local Variables ---!
       real(kind=cp) :: d

       d=Dot_Product(v,v)
       If (d <= 0.0_cp) Then
          ierr=-1
       Else
          ierr=0
          d=Sqrt(d)
          v=v/d
       End If

       Return
    End Subroutine normal

    !!----
    !!---- Subroutine Normal_Beam_Angles(wav,ub,h,sig,anbcal,ier,zer,nusign)
    !!----    real(kind=cp),                 Intent(In)         :: wav    ! Wavelength
    !!----    real(kind=cp), dimension(3,3), Intent(In)         :: ub     ! [UB] matrix
    !!----    real(kind=cp), dimension(3),   Intent(In)         :: h      ! Miller indices of reflection h
    !!----    Integer,                       Intent(In Out)     :: sig    ! -1 for negative gammas
    !!----    real(kind=cp), dimension(:),   Intent(Out)        :: anbcal ! Normal beam angles
    !!----    Integer,                       Intent(In Out)     :: ier    ! Zero correction on input, Error flag on output
    !!----    real(kind=cp), dimension(3),   Intent(In),optional:: zer    ! zero corrections of NB angles in degrees
    !!----    integer,                       Intent(In),optional:: nusign ! Sign of the nu-angle (=-1 if nu is positive when elevation towards z < 0)
    !!----
    !!--<<    Based on the subroutine CANNB from  A.FILHOL 01-Avr-1981 & 08-Aou-1984
    !!----    Calculation of the normal-beam diffraction angles for reflection h()
    !!----    from the [UB} matrix of the crystal.
    !!----    On input:
    !!----    sig      : 1/-1 defines the sign for "GAMMA"
    !!----    If zer is provided, the calculated angles are corrected for zero shifts
    !!----
    !!----    On output:
    !!----      Angles(calculated) = Angles(theoretical) + zero shifts
    !!----      ANBCAL(1:4) -> GAMMA, OMEGA(NB), NU, THETA   (in degrees)
    !!----     -> Error Flag
    !!----        ier  : 0   all angles are calculable
    !!----             : 1   reflection hors sphere d'Ewald
    !!----             : 2   reflection dans zone aveugle  ==> angle "NU"
    !!----             : 3   reflection dans zone aveugle  ==> angle "GAMMA"
    !!-->>             : 4   H,K,L tous les trois nuls
    !!----
    !!---- Created: April 2008, Updated: April 2013 (JRC)
    !!
    Subroutine Normal_Beam_Angles(wav,ub,h,sig,anbcal,ier,zer,nusign)
       !---- Arguments ----!
       real(kind=cp),                 Intent(In)          :: wav
       real(kind=cp), dimension(3,3), Intent(In)          :: ub
       real(kind=cp), dimension(3),   Intent(In)          :: h
       Integer,                       Intent(In Out)      :: sig
       real(kind=cp), dimension(:),   Intent(Out)         :: anbcal
       Integer,                       Intent(Out)         :: ier
       real(kind=cp), dimension(3),   Intent(In),optional :: zer
       integer,                       Intent(In),optional :: nusign

       !---- Local Variables ----!
       real(kind=cp), Dimension(3) ::  z1,zo1
       real(kind=cp)               ::  ds,dpl, thc, sthc,snu,cnu,cga,sga,oma,omb,tom1,tom2,ome,wav2,sn

       ier=0
       anbcal(:)=0.0_cp
       sn=1.0_cp
       if(present(nusign)) sn=real(nusign)

       wav2=wav*wav
       If(sig == 0) sig=1

       !---- Calculation of angle THETA and tests, Ewald sphere limits , etc.

       z1=matmul(ub,h)                ! h in laboratory system
       ds=sqrt(dot_product(z1,z1))    ! |h|=dstar
       zo1=z1/ds                      ! unit vector along z1
       dpl= z1(1)*z1(1) + z1(2)*z1(2) ! Projection of z1 on equatorial plane

       If (ds < 0.000001_cp) Then
          ier=4
          Return
       End If

       !---- test limit Ewald sphere
       sthc = ds*wav*0.5_cp
       If(Abs(sthc) > 1.0_cp) Then
         ier=1
         Return
       End If

       !-- Bragg angle THETA
       thc = Asind(sthc)
       thc = real(sig)*thc

       !-- NU   (Elevation angle of detector)
       snu=sn*wav*z1(3)         !Sin(nu) sn=1 for nu > 0 when elevation w.r.t. positive z-axis
       If (Abs(snu) > 1.0_cp  .or. dpl < 0.00001_cp) Then
          ier=2
          Return
       End If
       cnu=Sqrt(1.0_cp-snu*snu)   !Cos(nu)

       !-- GAM    (Projection of 2THETA on the equatorial plane)
       cga=((1.0_cp+cnu*cnu)-wav2*dpl)
       cga=cga/(2.0_cp*cnu)
       If (Abs(cga) > 1.0_cp) Then
          ier=3
          Return
       End If
       sga=sig*Sqrt(1.0_cp-cga*cga)

       !-- OME(NB)   (Rotation angle of the sample)
       oma=sig*cnu*sga
       omb=1.0_cp-cnu*cga
       tom1=omb*z1(1)+oma*z1(2)
       tom2=oma*z1(1)-omb*z1(2)
       ome=Atan2d(sig*tom1,tom2)

       !- Final angle values
       anbcal(1)=sig*Acosd(cga)
       anbcal(2)=ome
       anbcal(3)=Asind(snu)
       anbcal(4)=thc
       if (present(zer)) then
          anbcal(1:3) = anbcal(1:3) + zer(1:3)
       end if

       Return
    End Subroutine Normal_Beam_Angles

    !!----
    !!---- Subroutine Phi_mat(phi,dum)
    !!----    real(kind=cp), Intent(In)                   :: phi
    !!----    real(kind=cp), Intent(Out), Dimension(3,3)  :: dum
    !!----
    !!----    Calculate the Busing and Levy conventional rotation matrix for PHI
    !!----    or OMEGA [eq. 8 & 10]. The PHI/OMEGA angle must be provided in degrees.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Phi_mat(phi,dum)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                 :: phi
       real(kind=cp), Intent(Out), Dimension(3,3):: dum

       dum=0.0_cp
       dum(1,1)= cosd(phi)
       dum(1,2)= sind(phi)
       dum(2,1)=-dum(1,2)
       dum(2,2)= dum(1,1)
       dum(3,3)= 1.0_cp

       Return
    End Subroutine Phi_mat

    !!----
    !!---- Subroutine Psi_mat(psi,dum)
    !!----    real(kind=cp), Intent(In)                   :: psi
    !!----    real(kind=cp), Intent(Out), Dimension(3,3)  :: dum
    !!----
    !!----    Calculate the Busing and Levy  [eq. 46] conventional rotation matrix for PSI
    !!----    The PSI angle must be provided in degrees.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Psi_mat(psi,dum)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                 :: psi
       real(kind=cp), Intent(Out), Dimension(3,3):: dum

       dum=0.0_cp
       dum(1,1)=1.0_cp
       dum(2,2)= cosd(psi)
       dum(2,3)=-sind(psi)
       dum(3,2)=-dum(2,3)
       dum(3,3)= dum(2,2)

       Return
    End Subroutine Psi_mat

    !!----
    !!---- Subroutine refvec(vhkl,ub,vs,vz,ierr)
    !!----    real(kind=cp), Intent(In),  Dimension(3)    :: vhkl
    !!----    real(kind=cp), Intent(In),  Dimension(3,3)  :: ub
    !!----    real(kind=cp), Intent(Out), Dimension(3)    :: vs,vz
    !!----    Integer,       Intent(Out)                  :: ierr
    !!----
    !!----    Calculate vs,vz as reference vectors for defining Psi=0
    !!----    The B-L convention is that Psi=0 when the reflection hkl is
    !!----    in diffraction position and the c* is in the plane defined
    !!----    by vhkl and vz (z-axis of the laboratory system) for all
    !!----    reflections except when vhkl is parallel to c* in which case
    !!----    the vector b* plays the role of c* in the above prescription.
    !!----    The vector vhkl is provided with components in the reciprocal
    !!----    lattice.
    !!----
    !!---- Update: July 2008
    !!
    Subroutine refvec(vhkl,ub,vs,vz,ierr)
       !---- Arguments ----!
       real(kind=cp), Intent(In),  Dimension(3)    :: vhkl
       real(kind=cp), Intent(In),  Dimension(3,3)  :: ub
       real(kind=cp), Intent(Out), Dimension(3)    :: vs,vz
       Integer,       Intent(Out)                  :: ierr

       !--- Local Variables ---!
       real(kind=cp),Dimension(3) :: hn,h0
       real(kind=cp),Dimension(3) :: h1=(/0.0_cp,0.0_cp,1.0_cp/),h2=(/0.0_cp,1.0_cp,0.0_cp/),v0=(/0.0_cp,0.0_cp,1.0_cp/)

       !---- Test if VHKL is parallel to H1
       hn=vhkl
       !Call normal(hn,ierr)  Not needed at all => misleading because vhkl is in reciprocal lattice
       !Check that the vector is non-null
       if(sum(abs(hn)) < 0.0001_cp) then
         ierr=-1
         return
       else
         ierr=0
       end if
       h0=Cross_Product(hn,h1)
       If (Sum(Abs(h0)) > 0.0001_cp) Then
          h0=h1         !vhkl IS NOT parallel to c*(=h1), so h1 can be used as reference
       Else
          h0=h2         !vhkl IS parallel to c*(=h1), so h2=b* is used as reference
       End If
       vs=Matmul(ub,h0) !Put the reciprocal vector c* (or b*) in the laboratory system
       vz=v0            !Z-xis of the laboratory system

       Return
    End Subroutine refvec

    !!----
    !!---- Subroutine s4cnb(angl_4C,angl_NB,ierr)
    !!----    real(kind=cp), dimension(4), Intent(In)     :: angl_4C
    !!----    real(kind=cp), dimension(3), Intent(Out)    :: angl_NB
    !!----    Integer,                     Intent(Out)    :: ierr
    !!----
    !!--<<    Subroutine *** S4CNB ***   A.Filhol & M.Thomas
    !!----    Conversion of diffraction angles from the geometry
    !!----    4-CERCLES === TO ==> NORMAL-BEAM
    !!----        "ANGL_4C()"                       "ANGL_NB()"
    !!----   (/ 2Theta, Omega, Chi, Phi /)   (/ Gamma, Omega_nb, nu /)
    !!----    Error messages:
    !!----      IERR=0  O.K.
    !!----      IERR=1  calculation of nu impossible
    !!----      IERR=2  calculation of gamma impossible
    !!-->>      IERR=3  calculation of phi impossible
    !!----
    !!---- Update: April 2008
    !!
    Subroutine s4cnb(angl_4C,angl_NB,ierr)
       !---- Arguments ----!
       real(kind=cp), dimension(4), Intent(In)     :: angl_4C
       real(kind=cp), dimension(3), Intent(Out)    :: angl_NB
       Integer,                     Intent(Out)    :: ierr

       !--- Local variables ---!
       real(kind=cp) :: si, sint, sint2, sinnu, phi

       si=1.0_cp
       ierr=0

       !angl_4C=(/ 2Theta, Omega, Chi, Phi /)
       !angl_NB=(/ Gamma, Omega_nb, nu /)

       !.....NU NB
       sint=Sind(angl_4C(1)/2.0_cp)
       sint2=sint*sint
       angl_NB(3)=2.0_cp*Sind(angl_4C(3)) * sint
       If (Abs(angl_NB(3)) > 1.0_cp ) Then
          ierr=1
          return
       End If
       angl_NB(3)= Asind(angl_NB(3))
       angl_NB(3)= Abs(angl_NB(3)) *Sign(si,angl_4C(3))
       sinnu = Sind(angl_NB(3))

       !.....GAMMA NB
       angl_NB(1)=Cosd(angl_4C(1))/Cosd(angl_NB(3))
       If (Abs(angl_NB(1)) > 1.0_cp ) Then
          ierr=2
          Return
       End If
       angl_NB(1)=Acosd(angl_NB(1))
       angl_NB(1)=angl_NB(1)*Sign(si,angl_4C(1))

       !.....OMEGA NB
       phi=2.0_cp*sint2/Sqrt(4.0_cp*sint2-sinnu**2)
       If (Abs(phi) > 1.0_cp ) Then
          ierr=3
          Return
       End If
       phi=Acosd(phi)-90.0_cp
       phi=phi*Sign(si,angl_4C(1))
       angl_NB(2)=angl_4C(4)-phi

       return
    End Subroutine s4cnb

    !!---
    !!--- Subroutine Set_PSD()
    !!---
    !!---- Update: July 2008
    !!
    Subroutine Set_PSD(dist,cg,ag,nh,nv,ip)
       real(kind=cp), optional, intent(in) :: dist,cg,ag
       integer,       optional, intent(in) :: nh,nv,ip

       if(present(dist) .and. present(cg) .and. present(ag) .and. present(nh) &
                        .and. present(nv) .and. present(ip)) then
          psd%xoff   = 0.0_cp; psd%yoff=0.0_cp; psd%zoff=0.0_cp
          psd%radius = dist
          psd%cgap   = cg
          psd%agap   = ag
          psd%ncat   = nh
          psd%nano   = nv
          psd%ipsd   = ip
       else
          psd%xoff   = Current_Instrm%det_offsets(1)
          psd%yoff   = Current_Instrm%det_offsets(2)
          psd%zoff   = Current_Instrm%det_offsets(3)
          psd%radius = Current_Instrm%dist_samp_detector
          psd%cgap   = Current_Instrm%cgap
          psd%agap   = Current_Instrm%agap
          psd%ncat   = Current_Instrm%np_horiz
          psd%nano   = Current_Instrm%np_vert
          psd%ipsd   = Current_Instrm%ipsd
       end if
       psd_set    = .true.

       return
    End Subroutine Set_PSD

    !!----
    !!---- Subroutine snb4c(angl_NB,angl_4C)
    !!----    real(kind=cp), dimension(3), Intent(In )  :: angl_NB(3)
    !!----    real(kind=cp), dimension(4), Intent(Out)  :: angl_4C(4)
    !!--<<
    !!----    Subroutine ***SNB4C ***      A.Filhol & M.Thomas (I.L.L.)
    !!----   Conversion of diffraction angles from the geometry
    !!----   NORMAL-BEAM  === TO ==>  4-CERCLES
    !!-->>      "ANGL_NB()"                         "ANGL_4C()"
    !!----    (/ Gamma, Omega_nb, nu /)  (/ 2Theta, Omega, Chi, Phi /)
    !!----
    !!---- Update: April 2008
    !!
    Subroutine snb4c(angl_NB,angl_4C)
       !---- Arguments ----!
       real(kind=cp), dimension(3), Intent(In )  :: angl_NB(3)
       real(kind=cp), dimension(4), Intent(Out)  :: angl_4C(4)

       !--- Local variables ---!
       real(kind=cp) :: si, sint, sint2, sinnu, phi

       si=1.0_cp
       !.....2THETA 4C
       angl_4C(1)=Acosd(Cosd(angl_NB(1))*Cosd(angl_NB(3)))
       angl_4C(1)=angl_4C(1)*Sign(si,angl_NB(1))

       !.....OMEGA 4C
       angl_4C(2)=angl_4C(1)/2.0_cp

       !.....CHI 4C
       sint=Sind(angl_4C(2))
       sint2=sint*sint
       sinnu=Sind(angl_NB(3))
       angl_4C(3)=Asind(sinnu/(2.0_cp* sint))

       ! ***** On D15, both NU & 2-THETA might be positive or negative
       angl_4C(3)= Abs(angl_4C(3)) *Sign(si,angl_NB(3))

       !.....PHI 4C
       phi=2.0_cp*sint2/Sqrt(4.0_cp*sint2-sinnu**2)
       phi=Acosd(phi)-90.0_cp
       phi=phi*Sign(si,angl_NB(1))
       angl_4C(4)=angl_NB(2)+phi

       Return
    End Subroutine snb4c

    !!----
    !!---- Subroutine sxdpsd(mpsd,gamm,wave,nup,gamp,xobs,zobs, xcel,time,zcel,ierr)
    !!----    Integer, Intent(In)          :: mpsd
    !!----    real(kind=cp),    Intent(In)          :: gamm
    !!----    real(kind=cp),    Intent(In Out)      :: wave,nup,gamp
    !!----    real(kind=cp),    Intent(Out)         :: xobs,zobs, xcel,time,zcel
    !!----    Integer, Intent(Out)         :: ierr
    !!----
    !!--<<    Original comments:
    !!----
    !!----
    !!----    SPECIFICALLY FOR SXD ON SNS AT RAL.       R.F.D STANSFIELD DEC-84
    !!----
    !!----    The coordinate system adopted, whether the origin is at the
    !!----    moderator, at the sample (called the fixed laboratory system),
    !!----    or at the surface of the PSD when positioned at 0 degrees;
    !!----    is Y parallel to the beam, X in the horizontal plane on the
    !!----    diffraction side, and Z vertical. Hence if neutrons are diffracted
    !!----    to the left, Z is vertically down.
    !!----    The PSD is driven to an angle GamM. When GamM=0, the direct
    !!----    beam strikes a perfectly aligned PSD at its centre C.
    !!----    Call this point in space A. A is where we define NuP=0, GamP=0.
    !!----    The coordinates of A with respect to the sample are (0, Distsd, 0);
    !!----    and with respect to the moderator are (0, Distms+Distsd, 0).
    !!----    For a mis-aligned detector, the coordinates of A with respect to
    !!----    C are the translational offsets (Xoff, Yoff, Zoff) in mm.
    !!----    With the PSD at a general position GamM, the point where the
    !!----    direct beam struck it is rotated to O, where now NuP=0, GamP=GamM.
    !!----    For convenience, define a new cartesian system by rotating the
    !!----    axes of the fixed laboratory system about the vertical, such that
    !!----    X is now along the line joining the sample and O on the PSD.
    !!----    In this system, the coordinates of a Bragg peak P with respect to
    !!----    C are (Xobs, 0, Zobs) in mm, or (Xcel, 0, Zcel) in pixels.
    !!----    Hence the coordinates of P with respect to O, in this system, are:
    !!----     (x)     (Xobs   + Xoff)
    !!----     (y)  =  (Distsd + Yoff)  and: Tan(GamP - GamM) = x/y
    !!----     (z)     (Zobs   + Zoff)       Tan(NuP)         = z/SQRT(x*x + y*y)
    !!----    The PSD front surface measures Dimx by Dimy mm, and is divided
    !!----    into Nxcel by Nzcel pixels.
    !!----
    !!----                 -------------------------       At GamM=0 the direct
    !!----                 I          PSD          I       beam hits the detector
    !!----                 I     front surface     I       at A (NuP=0, GamP=0)
    !!----    Gamma________I_______________.O______I_________________________.A
    !!----           Zoff  I               !       I                         !
    !!----        X________I___________.C  !       I                         !
    !!----                 I           !   !       I                         !
    !!----           Zobs  I           !   !       I                         !
    !!----         ________I____.P     !   !       I                         !
    !!----                 I    !      !   !       I                         !
    !!----                 -----!------!---!--------                         !
    !!----                      ! Xobs !Xoff                                 !
    !!----                      !      !   !                                 !
    !!----                             Z  Nu                              GamM=0
    !!----
    !!----    Time is the time coordinate (bin) relative to an elapsed time
    !!----    Toff after the emission of a pulse at the moderator.
    !!----    The effect of moderator thickness on Time is NOT included yet.
    !!----    Distot is the total distance travelled from the moderator to a
    !!----    particular pixel on the PSD surface, in a total time Timtot.
    !!----    Wave = Velcon * Timtot / Distot, where Velcon is the velocity of a
    !!----    1 Angstrom neutron in Km/s or mm/mms.
    !!----
    !!----    MPSD +VE -> Find Wave, NuP, GamP; given GamM, Xcel, Zcel and Time
    !!----    MPSD -VE -> Find Xcel, Zcel, Time; given Wave, NuP, GamP and GamM
    !!----    (Routine not tested, probably obsolete for present SXD!!!)
    !!-->>
    !!----
    !!---- Update: April 2008
    !!
    Subroutine sxdpsd(mpsd,gamm,wave,nup,gamp,xobs,zobs, xcel,time,zcel,ierr)
       !---- Arguments ----!
       Integer, Intent(In)                   :: mpsd
       real(kind=cp),    Intent(In)          :: gamm
       real(kind=cp),    Intent(In Out)      :: wave,nup,gamp
       real(kind=cp),    Intent(Out)         :: xobs,zobs, xcel,time,zcel
       Integer, Intent(Out)                  :: ierr

       !--- Local variables ---!
       real(kind=cp) :: xmax,zmax, a, b, c, d, e, timtot, distot, delga, td, f, tn

       xmax=real(sxd%nxcel/2)
       zmax=real(sxd%nzcel/2)
       If (mpsd < 0) Then
          !---- Xcel and Zcel are the coords. in pixels of P from C, on the PSD
          !---- NuP, GamP are the angular coordinates of P from the sample
          delga=gamp - gamm
          If (Abs(delga) >= 90.0) then
             ierr=-1
             return
          End if
          b=sxd%distsd+sxd%yoff
          td=tand(delga)
          a=b*td
          xobs=a - sxd%xoff
          xcel=xobs*real(sxd%nxcel)/sxd%dimx
          f=Sqrt(1.0 + td*td)
          tn=tand(nup)
          c=b*f*tn
          zobs=c -sxd%zoff
          zcel=zobs*real(sxd%nzcel)/sxd%dimz
          If (Abs(xcel) > xmax .or. Abs(zcel) > zmax) then
             ierr=-1
             return
          End If
          e=Sqrt(a*a + b*b + c*c)
          distot=sxd%distms+e
          timtot=wave*distot/sxd%velcon
          time=timtot-sxd%toff
       Else
          If (Abs(xcel) > xmax.or.abs(zcel) > zmax) Then
             ierr=-1
             return
          End If

          !---- Xobs and Zobs are the coords. in mm of P from C on the PSD surface
          !---- (A,B,C) are the coordinates of P from the sample in the lab system

          xobs=xcel*sxd%dimx/real(sxd%nxcel)
          zobs=zcel*sxd%dimz/real(sxd%nzcel)
          a=xobs  +sxd%xoff
          b=sxd%distsd+sxd%yoff
          c=zobs  +sxd%zoff
          d=Sqrt(a*a + b*b)
          e=Sqrt(c*c + d*d)
          gamp=gamm + atan2d(a,b)
          nup=atan2d(c,d)
          timtot=time+sxd%toff
          distot=sxd%distms+e
          wave=sxd%velcon*timtot/distot
       End If

       Return
    End Subroutine sxdpsd

    !!----
    !!---- Subroutine triple(v1,v2,tv,ierr)
    !!----    real(kind=cp),    Intent(In Out), Dimension(3)  :: v1,v2
    !!----    real(kind=cp),    Intent(Out),    Dimension(3,3):: tv
    !!----    Integer, Intent(Out)                            :: ierr
    !!----
    !!----    Construct orthonormal triplet matrix TV, with column vectors :
    !!----    V1, (V1 x V2) x V1, V1 x V2.
    !!----
    !!---- Update: July 2008
    !!
    Subroutine triple(v1,v2,tv,ierr)
       !---- Arguments ----!
       real(kind=cp),    Intent(In Out), Dimension(3)  :: v1,v2  !they come back normalized and V2 perp. to V1
       real(kind=cp),    Intent(Out),    Dimension(3,3):: tv
       Integer, Intent(Out)                            :: ierr

       !--- Local Variables ---!
       real(kind=cp),Dimension(3) :: v3

       Call normal(v1,ierr)
       v3=Cross_Product(v1,v2)
       Call normal(v3,ierr)
       v2=Cross_Product(v3,v1)
       Call normal(v2,ierr)
       If(ierr /= 0) Return
       tv(:,1)=v1(:)
       tv(:,2)=v2(:)
       tv(:,3)=v3(:)

       Return
    End Subroutine triple

    !!----
    !!---- Subroutine z1frfc(wave,tth,om,ch,ph,z1)
    !!----    real(kind=cp), Intent(In)               :: wave,tth,om,ch,ph
    !!----    real(kind=cp), Intent(Out),Dimension(3) :: z1
    !!----
    !!----    Z1 from Four Circle diffractometer angles
    !!----    Calculate diffraction vector Z1 from TTH, OM, CH, PH (Need not
    !!----    be bisecting, but Z4 is assumed to be in the equatorial plane)
    !!----
    !!---- Update: July 2008
    !!
    Subroutine z1frfc(wave,tth,om,ch,ph,z1)
       !---- Argument ----!
       real(kind=cp), Intent(In)               :: wave,tth,om,ch,ph
       real(kind=cp), Intent(Out),Dimension(3) :: z1

       !--- Local Variables ---!
       real(kind=cp),Dimension(3) ::  z4
       !real(kind=cp) :: th  !Not needed

                                    ! Original: th=(tth/2.0)
       z4(1)=sind(tth)/wave         ! ( 2.0*sind(th)*cosd(th))/wave
       z4(2)=(cosd(tth)-1.0)/wave   ! (-2.0*sind(th)*sind(th))/wave
       z4(3)=0.0                    ! 0.0
       Call z1frz4(z4,om,ch,ph,z1)

       Return
    End Subroutine z1frfc

    !!----
    !!---- Subroutine z1frmd(wave,ch,ph,ga,om,nu,z1)
    !!----    real(kind=cp), Intent(In)                 :: wave,ch,ph,ga,om,nu
    !!----    real(kind=cp), Intent(Out), Dimension(3)  :: z1
    !!----
    !!----    Z1 from Four Circle + PSD multi-detector diffractometer angles
    !!----    Calculate diffraction vector Z1 from CH, PH, GA, OM, NU
    !!----    for a multi-detector. The angles Chi,Phi,Gamma,Omega and Nu for
    !!----    the equatorial plane are Chi,Phi,2Theta and Omega (Nu=0)
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z1frmd(wave,ch,ph,ga,om,nu,z1)
       !---- Argument ----!
       real(kind=cp), Intent(In)                 :: wave,ch,ph,ga,om,nu
       real(kind=cp), Intent(Out), Dimension(3)  :: z1

       !--- Local Variables ---!
       real(kind=cp), Dimension(3)   ::  z3
       Call z1frnb(wave,ga,om,nu,z3)
       Call z1frz3(z3,ch,ph,z1)

       Return
    End Subroutine z1frmd

    !!----
    !!---- Subroutine z1frnb(wave,ga,om,nu,z1)
    !!----    real(kind=cp), Intent(In)                 :: wave,ga,om,nu
    !!----    real(kind=cp), Intent(Out), Dimension(3)  :: z1
    !!----
    !!----    Z1 from Normal Beam diffractometer angles
    !!----    Calculate diffraction vector Z1 from GA, OM, NU, assuming CH=PH=0
    !!----    This is is the normal beam geometry for a Lifting arm detector or
    !!----    for a PSD with a single Omega axis for the sample.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z1frnb(wave,ga,om,nu,z1)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                 :: wave,ga,om,nu
       real(kind=cp), Intent(Out), Dimension(3)  :: z1

       !--- Local Variables ---!
       real(kind=cp), Dimension(3,3) ::  omg
       real(kind=cp), Dimension(3)   ::  z4

       Call z4frgn(wave,ga,nu,z4)
       Call Phi_mat(om,omg)
       z1=Matmul(Transpose(omg),z4)

       Return
    End Subroutine z1frnb

    !!----
    !!---- Subroutine z1frz2(z2,ph,z1)
    !!----    real(kind=cp), Intent(In ), Dimension(3)  :: z2
    !!----    real(kind=cp), Intent(In )                :: ph
    !!----    real(kind=cp), Intent(Out), Dimension(3)  :: z1
    !!----
    !!----    Calculate Z1 = [PHI]T.Z2
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z1frz2(z2,ph,z1)
       !---- Arguments ----!
       real(kind=cp), Intent(In ), Dimension(3)  :: z2
       real(kind=cp), Intent(In )                :: ph
       real(kind=cp), Intent(Out), Dimension(3)  :: z1

       !--- Local Variables ---!
       real(kind=cp),Dimension(3,3)  :: phim

       Call Phi_mat(ph,phim)
       z1=Matmul(Transpose(phim),z2)

       Return
    End Subroutine z1frz2

    !!----
    !!----   Subroutine z1frz3(z3,ch,ph,z1)
    !!----      real(kind=cp), Intent(In), Dimension(3)   :: z3
    !!----      real(kind=cp), Intent(In)                 :: ch,ph
    !!----      real(kind=cp), Intent(Out), Dimension(3)  :: z1
    !!----
    !!----    CALCULATE Z1 = [PHI]T.[CHI]T.Z3
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z1frz3(z3,ch,ph,z1)
       !---- Arguments ----!
       real(kind=cp), Intent(In), Dimension(3)   :: z3
       real(kind=cp), Intent(In)                 :: ch,ph
       real(kind=cp), Intent(Out), Dimension(3)  :: z1

       !--- Local Variables ---!
       real(kind=cp),Dimension(3)    :: z2
       real(kind=cp),Dimension(3,3)  :: chim

       !----
       Call Chi_mat(ch,chim)
       z2=Matmul(Transpose(chim),z3)
       Call z1frz2(z2,ph,z1)

       Return
    End Subroutine z1frz3

    !!----
    !!---- Subroutine z1frz4(z4,om,ch,ph,z1)
    !!----    real(kind=cp), Intent(In), Dimension(3)   :: z4
    !!----    real(kind=cp), Intent(In)                 :: om,ch,ph
    !!----    real(kind=cp), Intent(Out), Dimension(3)  :: z1
    !!----
    !!----    Calculate Z1 = [PHI]T.[CHI]T.[OM]T.Z4
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z1frz4(z4,om,ch,ph,z1)
       !---- Arguments ----!
       real(kind=cp), Intent(In), Dimension(3)   :: z4
       real(kind=cp), Intent(In)                 :: om,ch,ph
       real(kind=cp), Intent(Out), Dimension(3)  :: z1

       !--- Local Variables ---!
       real(kind=cp),Dimension(3)    :: z3
       real(kind=cp),Dimension(3,3)  :: phim

       Call Phi_mat(om,phim)
       z3=Matmul(Transpose(phim),z4)
       Call z1frz3(z3,ch,ph,z1)

       Return
    End Subroutine z1frz4

    !!----
    !!---- Subroutine z2frz1(z1,ph,z2)
    !!----    real(kind=cp), Intent(In),Dimension(3)    :: z1
    !!----    real(kind=cp), Intent(In)                 :: ph
    !!----    real(kind=cp), Intent(Out),Dimension(3)   :: z2
    !!----
    !!----    Calculate Z2 = [PHI].Z1
    !!----    The reciprocal vector Z1 is given in cartesian components with
    !!----    respect to the laboratory system.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z2frz1(z1,ph,z2)
       !---- Arguments ----!
       real(kind=cp), Intent(In),Dimension(3)    :: z1
       real(kind=cp), Intent(In)                 :: ph
       real(kind=cp), Intent(Out),Dimension(3)   :: z2

       !--- Local Variables ---!
       real(kind=cp), Dimension(3,3) ::  dum

       Call Phi_mat(ph,dum)
       z2=Matmul(dum,z1)

       Return
    End Subroutine z2frz1

    !!----
    !!---- Subroutine z3frz1(z1,ch,ph,z3)
    !!----    real(kind=cp), Intent(In), Dimension(3)    :: z1
    !!----    real(kind=cp), Intent(In)                  :: ch,ph
    !!----    real(kind=cp), Intent(Out), Dimension(3)   :: z3
    !!----
    !!----    Calculate Z3 = [CHI].Z2 = [CHI].[PHI].Z1
    !!----    The reciprocal vector Z1 is given in cartesian components with
    !!----    respect to the laboratory system.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z3frz1(z1,ch,ph,z3)
       !---- Arguments ----!
       real(kind=cp), Intent(In), Dimension(3)    :: z1
       real(kind=cp), Intent(In)                  :: ch,ph
       real(kind=cp), Intent(Out), Dimension(3)   :: z3

       !--- Local Variables ---!
       real(kind=cp), Dimension(3)   ::  z2
       real(kind=cp), Dimension(3,3) ::  dum

       Call z2frz1(z1,ph,z2)
       Call Chi_mat(ch,dum)
       z3=Matmul(dum,z2)

       Return
    End Subroutine z3frz1

    !!----
    !!---- Subroutine z4frgn(wave,ga,nu,z4)
    !!----    real(kind=cp), Intent(In)                 :: wave
    !!----    real(kind=cp), Intent(In)                 :: ga,nu
    !!----    real(kind=cp), Intent(Out), Dimension(3)  :: z4
    !!----
    !!----    Calculates diffraction vector of a reflection in the
    !!----    laboratory system from the angles GA and NU.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z4frgn(wave,ga,nu,z4)
       !---- Arguments ----!
       real(kind=cp), Intent(In)                 :: wave
       real(kind=cp), Intent(In)                 :: ga,nu
       real(kind=cp), Intent(Out), Dimension(3)  :: z4

       z4(1)=( sind(ga)*cosd(nu)     )/wave
       z4(2)=( cosd(ga)*cosd(nu)-1.0 )/wave  ! without 1 => Cartesian system on the
       z4(3)=( sind(nu)              )/wave  ! centre of Ewald sphere

       Return
    End Subroutine z4frgn

    !!----
    !!---- Subroutine z4frz1(z1,om,ch,ph,z4)
    !!----    real(kind=cp), Intent(In), Dimension(3)     :: z1
    !!----    real(kind=cp), Intent(In)                   :: om,ch,ph
    !!----    real(kind=cp), Intent(Out), Dimension(3)    :: z4
    !!----
    !!----    Calculate Z4 = [OM].[CHI].[PHI].Z1
    !!----
    !!---- Update: April 2008
    !!
    Subroutine z4frz1(z1,om,ch,ph,z4)
       !---- Arguments ----!
       real(kind=cp), Intent(In), Dimension(3)     :: z1
       real(kind=cp), Intent(In)                   :: om,ch,ph
       real(kind=cp), Intent(Out), Dimension(3)    :: z4

       !--- Local Variables ---!
       real(kind=cp), Dimension(3)  ::  z3
       real(kind=cp), Dimension(3,3)::  omg

       Call z3frz1(z1,ch,ph,z3)
       Call Phi_mat(om,omg)
       z4=Matmul(omg,z3)

       Return
    End Subroutine z4frz1

 End Module CFML_Geometry_SXTAL
