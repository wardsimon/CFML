!!----  Program MCMAG, written by Philippe Lacorre in 1988
!!----  Adapted to Fortran 90 and some more modifications by J.Rodriguez-Carvajal
!!----  A new input file containing the whole set of instructions for running
!!----  the program has been implemented. Output for FullProf Studio and for
!!----  calculating the single crystal and powder neutron diffraction patterns
!!----  is being implemented (May 2015, JRC)
!!----
    Module mcm_inc
    !Global variables of the program MCMAG (no more COMMONs!)
    implicit none
    integer, parameter, public :: sp = selected_real_kind(6,30)
    integer, parameter, public :: dp = selected_real_kind(14,150)
    real(kind=dp), parameter, public ::  pi = 3.141592653589793238463_dp
    integer, parameter, public ::    &
              numax=30,    &    ! maximal number of cells along the U axis
              nvmax=30,    &    ! maximal number of cells along the V axis
              nwmax=30,    &    ! maximal number of cells along the W axis
              NAMAX=100,   &    ! maximal number of sites per cell
              NNMAX=50,    &    ! maximal number of neighbours per site
              NMOMAX=200,  &    ! maximal number of isotropic coupling
                                !   constants which can be modified
                                !   interactively
              nal=8,       &    ! maximal number of site moments to be
                                !   written on a single line of the output
                                !   file MYFILE.LIS
              I01=31,      &    !  input unit number for file MYFILE.DAT
              I02=32,      &    !  input unit number for file MYFILE.INI
              I03=33,      &    ! output unit number for file MYFILE.RES
              I04=34,      &    ! output unit number for file MYFILE.SPI
              IBA=35,      &    ! output unit number for file MYFILE.REC
              LPT=7,       &    ! output unit number for file MYFILE.LIS
              ITTI=5,      &    ! input unit number for terminal
                                !   (interactive mode)
              ITTO=6            ! output unit number for terminal
                                !   (interactive mode)

          real(kind=sp), dimension(namax)                         :: d,ddix,ddiy,ddiz, dx,dy,dz, d0, amps, amps2, &
                                                                     fjsx,fjsy,fjsz,fjsm
          character(len=6), dimension(namax) :: namt=" ",scattf=" "
          character(len=2), dimension(namax) :: elem=" "
          real(kind=sp), dimension(3,namax)  :: xyz=0.0
          logical                            :: coordinates_read=.false.,rotated=.false.
          logical                            :: batch=.false.
          real(kind=sp), dimension(namax,numax+1,nvmax+1,nwmax+1) :: sx, sy, sz
          integer,       dimension(namax,nnmax,3)                 :: M !Connectivity matrix with relative translations
          integer,       dimension(namax,nnmax)                   :: nav !
          integer,       dimension(namax)                         :: isigx,isigy,isigz,isigm,nu,nv,nw

          integer,       dimension(namax)                         :: Nvs, &!Total number of neighbours of a site
                                                                     nat !
          real(kind=sp), dimension(namax,nnmax)                   :: ajjx,ajjy,ajjz,ajxy,ajxz,ajyx,ajyz,ajzx,ajzy, &
                                                                     ajxx,ajyy,ajzz
          real(kind=sp), dimension(nnmax)                         :: pt
          real(kind=dp), dimension(3)                             :: tco

          real(kind=sp), dimension(namax)                         :: ajsx,ajsy,ajsz

          real(kind=sp),dimension(namax,numax,nvmax,nwmax)        :: somx, somx2,somy2,somz2, somrx,somry,somrz, &
                                                                     somrx2,somry2,somrz2
          real(kind=sp),dimension(namax,numax+1,nvmax+1,nwmax+1)  :: somy,somz,comx,tmom

          integer,       dimension(nmomax)  ::ium,ivm,iwm,nm,numb
          real(kind=sp), dimension(nmomax)  ::ajmod


          integer(kind=8) :: npaj ! Number of accepted jumps

          integer ::      lu,lv,lw,   &! Number of unit cells along u,v,w
                          it,         &! Number of Montecarlo cycles
                          nom,        &! Number of thermalization cycles
                          ianb,       &! Boundary condirions (1: periodic, 0:free, -1:mixed)
                          na,         &! Total number of sites in the basic unit cell (if negative no single ion anisotropy)
                          iani,       &! Index for single ion anisotropy (set to -1 if na < 0, before puttin na=|na|)
                          iran,       &! Index for the spin-model (1:Ising, 2:planar XY, 3:Heisenberg, 4: q-state Potts planar)
                          jcod,       &! Code of interactions (0: isotropic, 1: anisotropic diagonal tensor, 2: anisotropic)
                          nmot,       &! Number of formula units per cell
                          nma,        &! Total number of unit cells: nma=lu*lv*lw
                          nq,         &! Number of q-states in q-Potts models
                          infa,nfa,   &! Option for averages (1:sample, 2:unit cell, 3: site, 4:mole)
                          nta,        &! Number of units for averages(1:sample, nma: cell, nta:site, nma*nmot: mole)
                          j1,j2,j3,   &! Indices for three sites (output)
                          last,       &! =1 for the last iteration
                          iwmx,iwmxy, &! Printing controls
                          nor,        &! Ordinal number of the moment to be re-oriented
                          im,iwrt=0,  &! Index, counter for cpu_time
                          ntemp        ! Total number of temperatures

          real(kind=sp):: t,           &! Initial temperature
                        coef,          &! Kirpatrick coefficient (>1 heating, <1 cooling)
                        tfin,          &! Final temperature
                        hk,            &! Strength of magnetic field in kelvins
                        hx0,hy0,hz0,   &! Direction of magnetic field
                        hsta,hste,hfin,&! Initial, step and final magnetic field
                        hx,hy,hz,      &! Effective components of the magnetic field
                        hdix,hdiy,hdiz,&! Intermediate components of the magnetic field
                        snom,          &! snom=FLOAT(it-nom)
                        amx,amy,amz,cpu,cpu_par,remaining
        real(kind=sp) :: delta, & ! Current energy difference between the old and new configuration
                         x,y,z, & ! Randomly generated components of a spin
                         dorix,doriy,doriz,dorx,dory,dorz  !Direction to rotate a spin (X-tal and Cartesian systems)
        real(kind=sp) :: fconst,rc  !Constraint functions
        integer       :: nfip,  &! if zero no random spin S=(x,y,z) is generated before calculating delta
                         nmult   ! Multiplier (2) applied to H and D in order to calculate the energy E1=1/2 E0

          character (len=132):: titl1,titl2,titl3,itit1,itit2,itit3
          character (len=50) :: titl
          character (len=10) :: bcond
          character (len=120):: nfil
          character (len=9)  :: avchar
          character (len=12) :: spimod
          character (len=1)  :: ansi,   &! Answer for correct file
                                ansmod, &! Answer for modifying locally the exchange interactions
                                ians,   &! R of I (Random initial configuration or read from input file)
                                ansnev, &! Printout of the configuration of magnetic moments (N: Never, E: End, A: All temperatures)
                                ansrf,  &! Output of averaged characteristics in file *.res
                                ansc,   &! Continue or not
                                ansh,   &! Interation of temperature (T) or Field (F)
                                cod,    &! Code for multiplier (M) or for additive coefficient (A)
                                ansm,   &! Printout of the coomponents of magnetic moments in crystallographic system
                                codt     ! Code '+' or 'x' depending on cod.
        real(kind=dp) :: tcojx,tcojy,tcojz, tcojxy,tcojxz,tcojyz,tcojyx,tcojzx,tcojzy !Tests of internal compatibility variables
        real(kind=dp) :: en,enm,en2m,mag2
        real(kind=dp), dimension(3,3) :: aa,ac,ad
        real(kind=sp) :: aax,aay,aaz,aaa,aab,aag !Cell parameters

  contains
      subroutine conv_matrix()
        !  This subroutine calls the following function:
        !     INVERT
         real(kind=dp) :: degre,aar
         degre=pi/180.0_dp
         aaa=aaa*degre
         aab=aab*degre
         aag=aag*degre

         aa(:,:)=0.0_dp
         ad(:,:)=0.0_dp

         aar=SQRT( 1.0-COS(aaa)*COS(aaa) -COS(aab)*COS(aab)  &
             -COS(aag)*COS(aag) +2*COS(aaa)*COS(aab)*COS(aag) )
         aa(1,1)=1.0
         aa(1,2)=COS(aag)
         aa(1,3)=COS(aab)
         aa(2,2)=SIN(aag)
         aa(2,3)=(COS(aab)*COS(aag)-COS(aaa))/SIN(aag)
         aa(3,3)=aar/SIN(aag)

         ad(1,1)=aax
         ad(1,2)=aay*COS(aag)
         ad(1,3)=aaz*COS(aab)
         ad(2,2)=aay*SIN(aag)
         ad(2,3)=aaz*(COS(aab)*COS(aag)-COS(aaa))/SIN(aag)
         ad(3,3)=aaz*aar/SIN(aag)

         aaa=aaa/degre
         aab=aab/degre
         aag=aag/degre
         !Calculates AC=INV(AA)
         ac=Invert(aa)
      end subroutine conv_matrix

      Function Invert(a) Result(b)
         !---- Arguments ----!
         real(kind=dp),dimension(3,3), intent(in) :: a
         real(kind=dp),dimension(3,3)             :: b

         !---- Local variables ----!
         real(kind=dp)  :: dmat

         b(1,1) =   a(2,2)*a(3,3)-a(2,3)*a(3,2)
         b(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
         b(3,1) =   a(2,1)*a(3,2)-a(2,2)*a(3,1)
         b(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
         b(2,2) =   a(1,1)*a(3,3)-a(1,3)*a(3,1)
         b(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
         b(1,3) =   a(1,2)*a(2,3)-a(1,3)*a(2,2)
         b(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
         b(3,3) =   a(1,1)*a(2,2)-a(1,2)*a(2,1)
         dmat = a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1) !determinant of A

         if (abs(dmat) > tiny(dmat)) then
            b= b/dmat
         else
            b=0.0_dp
         end if

         return
      End Function Invert


!***********************************************************************
!  The subroutines RANnD(RADIUS,X1,Y1,,) generate a pseudo-random
!  spin configuration that is a random vector suitable for the selected
!  spin model:
!      RAN1D(RADIUS,X1)         for    Ising spins
!      RAN2D(RADIUS,X1,Y1)      for    XY spins
!      RAN3D(RADIUS,X1,Y1,Z1)   for    Heisenberg spins
!      RANqD(RADIUS,X1,Y1)      for    q-state Potts planar (XY) spins
!  The subroutines RAN2D(RADIUS,X1,Y1) and RAN3D(RADIUS,X1,Y1,Z1)
!  originate from the Program Library of CERN Computer Center
!  (CERN - Geneva), with original names RAN2VS and RAN3VS.
!  These subroutines call no subroutine.
!  They call the following function :
!  RANF() : this portable real function [Program library of CERN
!           Computer Centre (CERN - Geneva)] is used to generate random
!           numbers between 0 and 1. This function may be slower than
!           dedicated random number generators. Therefore, in order to
!           save CPU time, users should replace RANF() by a call to
!           the equivalent function on the computer they use.
!           For instance, DEC-computer users should call the random
!           generator RAN(I1,I2). For sake of convenience, the calls to
!           this function are included in the current version of the
!           program as comment cards. After implementing and testing
!           the program (see manual), the use of RAN(I1,I2) or of a
!           faster computer random generator, is recommended.

!***********************************************************************
      SUBROUTINE ran1d(radius,x1,y1,z1,n)
      real(kind=sp), intent(in)  :: radius
      real(kind=sp), intent(out) :: x1,y1,z1
      integer,       intent(in)  :: n
!***********************************************************************
!  For a general description of the subroutines RANnD(RADIUS,X1,Y1...),
!  see the comments just above.
!  RAN1D(RADIUS,X1) generates pseudo-randomly a number +RADIUS or
!  -RADIUS.
!***********************************************************************
         IF (ranf()-0.5 > 0) THEN
           x1 = radius*dx(n)
           y1 = radius*dy(n)
           z1 = radius*dz(n)
         ELSE
           x1 = -radius*dx(n)
           y1 = -radius*dy(n)
           z1 = -radius*dz(n)
         END IF
         RETURN
      END SUBROUTINE ran1d


      SUBROUTINE ran2d(radius,x1,y1)
      real(kind=sp), intent(in)  :: radius
      real(kind=sp), intent(out) :: x1,y1
!***********************************************************************
!  For a general description of the subroutines RANnD(RADIUS,X1,Y1...),
!  see the comments just before the subroutine RAN1D(RADIUS,X1).
!  RAN2D(RADIUS,X1,Y1) generates pseudo-randomly a point on a circle of
!  radius RADIUS.
!  It originates from the Program Library of CERN Computer Center
!  (CERN - Geneva) with the original name RAN2VS.
!***********************************************************************
         real(kind=sp) :: x,y,rr,scale
         DO
           x  =  2.0*ranf() - 1.0
           y  =  2.0*ranf() - 1.0
           rr =  x*x + y*y
           IF(rr <= 1.0)  EXIT
         END DO
         scale  =  radius / SQRT(rr)
         x1  =  x * scale
         y1  =  y * scale
         RETURN
      END SUBROUTINE ran2d


      SUBROUTINE ran3d(radius,x1,y1,z1)
      real(kind=sp), intent(in)  :: radius
      real(kind=sp), intent(out) :: x1,y1,z1
      real(kind=sp) :: x,y,z,rr,scale
!***********************************************************************
!  For a general description of the subroutines RANnD(RADIUS,X1,Y1...),
!  see the comments just before the subroutine RAN1D(RADIUS,X1).
!  RAN3D(RADIUS,X1,Y1,Z1) generates pseudo-randomly a point at the
!  surface of a sphere of radius RADIUS.
!  It originates from the Program Library of CERN Computer Center
!  (CERN - Geneva) with the original name RAN3VS.
!***********************************************************************
         DO
           x  =  2.0*ranf() - 1.0
           y  =  2.0*ranf() - 1.0
           z  =  2.0*ranf() - 1.0
           rr =  x*x + y*y + z*z
           IF(rr <= 1.0)  EXIT
         END DO
         scale  =  radius / SQRT(rr)
         x1  =  x * scale
         y1  =  y * scale
         z1  =  z * scale
         RETURN
      END SUBROUTINE ran3d

      SUBROUTINE ranqd(radius,x1,y1)
      real(kind=sp), intent(in)  :: radius
      real(kind=sp), intent(out) :: x1,y1
      real(kind=sp) :: aq, uran
      integer :: i
!***********************************************************************
!  For a general description of the subroutines RANnD(RADIUS,X1,Y1...),
!  see the comments just before the subroutine RAN1D(RADIUS,X1).
!  RANqD(RADIUS,X1,Y1) generates pseudo-randomly a point on a circle of
!  radius RADIUS among q points spaced out of an angle 360/q degrees.
!***********************************************************************
         uran=ranf()
         DO i=1,nq
           aq = (FLOAT(i))/(FLOAT(nq))
           IF (uran < aq) THEN
             x1 = radius*COS(2*pi*(i-1)/nq)
             y1 = radius*SIN(2*pi*(i-1)/nq)
             EXIT
           END IF
         END DO
         RETURN
      END SUBROUTINE ranqd


      FUNCTION ranf() result(r)
        real(kind=sp)   :: r
        call random_number(r)
      END FUNCTION ranf

      Subroutine Write_Date_Time(lun,dtim)
        integer,         optional,intent(in) :: lun
        character(len=*),optional,intent(out):: dtim
        !--- Local variables ----!
        character (len=10) :: dat
        character (len=10) :: tim
        call date_and_time(date=dat,time=tim)
        if(present(lun)) &
        write(unit=lun,fmt="(/,4a)") &
          " => Date: ",dat(7:8)//"/"//dat(5:6)//"/"//dat(1:4),      &
            "  Time: ",tim(1:2)//":"//tim(3:4)//":"//tim(5:10)
        if(present(dtim)) &
         dtim="#   Date: "//dat(7:8)//"/"//dat(5:6)//"/"//dat(1:4)//      &
               "  Time: "//tim(1:2)//":"//tim(3:4)//":"//tim(5:10)
        return
      End Subroutine Write_Date_Time

      Function L_Case(Text) Result (Mtext)
         !---- Argument ----!
         character (len=*), intent(in) :: text
         character (len=len(text))     :: mtext

         !---- Local variables ----!
         integer, parameter :: inc = ICHAR("A") - ICHAR("a")
         integer            :: leng, pos

         mtext=text
         leng=len_trim(mtext)
         do pos=1,leng
            if (mtext(pos:pos) >= "A" .and. mtext(pos:pos) <= "Z")           &
                mtext(pos:pos) = CHAR ( ICHAR(mtext(pos:pos)) - inc )
         end do

         return
      End Function L_Case

      Function U_Case(Text) Result (Mtext)
         !---- Argument ----!
         character (len=*), intent(in) :: text
         character (len=len(text))     :: mtext

         !---- Local variables ----!
         integer, parameter :: inc = ICHAR("A") - ICHAR("a")
         integer            :: leng, pos

         mtext=text
         leng=len_trim(mtext)
         do pos=1,leng
            if (mtext(pos:pos) >= "a" .and. mtext(pos:pos) <= "z")           &
                mtext(pos:pos) = CHAR ( ICHAR(mtext(pos:pos)) + inc )
         end do

         return
      End Function U_Case

  End Module mcm_inc

!                  M   M   CCC  M   M   AA    GGG
!                  MM MM  C     MM MM  A  A  G
!                  M M M  C     M M M  AAAA  G GGG
!                  M   M  C     M   M  A  A  G   G
!  *************** M   M   CCC  M   M  A  A   GGG ****************
! *                                                               *
! *   MCMAG - A FORTRAN PROGRAM TO SIMULATE MAGNETIC STRUCTURES   *
! *                  Version May 2015 (JRC-ILL)                   *
! *                                                               *
! *                    written by P. Lacorre                      *
! *                modified by J. Rodriguez-Carvajal              *
! *  Institut Laue-Langevin, 156X, 38042 Grenoble Cedex - FRANCE  *
! *                             and                               *
! *  Laboratoire des Fluorures URA CNRS 449, Faculte des Sciences *
! *         Universite du Maine, Avenue Olivier-Messiaen          *
! *                 72017 Le Mans Cedex - FRANCE                  *
! *                                                               *
! *****************************************************************

!      MCMAG is a Fortran program to simulate magnetic models from
!  postulated coupling constants. Ground state configurations are
!  obtained by a simulated annealing algorithm [S. Kirkpatrick,
!  C.D. Gellatt Jr. and M.P. Vecchi, Science, 220 (1983) 671] and
!  the configuration space is explored with random walks by using
!  Metropolis importance sampling [N. Metropolis, A.W. Rosenbluth,
!  M.N. Rosenbluth, A.H. Teller and E. Teller, J. Chem. Phys. 21
!  (1953) 1087]. Any kind of 1D, 2D or 3D periodic model, magnetic
!  cluster as well as real magnetic structure can be simulated
!  without limitation of interacting distance.

!  Further details on the method used for this purpose are given
!  in the paper :
!  MCMAG: a computer program to simulate magnetic structures.
!  P. Lacorre and J. Pannetier, J. Magn. Magn. Mat. 71 (1987) 63.

!  Topological and magnetic parameters must be stored in a file
!  MYFILE.DAT (where MYFILE is a 6 characters filename given by
!  user). This input file contains a list of the magnetic
!  sites of the structure, together with a list of their neighbours
!  and corresponding exchange integrals, as well as anisotropy
!  coefficients and spin amplitudes (for a complete description see
!  the manual MCMAG.MAN). No symmetry element is taken into account
!  by the program: all spins are independent.
!  The so-called basic unit cell (b.u.c.) which is described in
!  the input file represents the smallest sample that can be used
!  for a simulation. Larger size samples are built up from the
!  stacking of such unit cells along the three crystallographic
!  axes and defined by the number of cells along these axes (namely
!  LU, LV and LW for axes U, V and W respectively). Three kinds of
!  boundary conditions are allowed:
!       - free edges
!       - periodic boundary conditions
!       - mixed conditions (<=> free edges along W and periodic
!         boundary conditions along U and V)

!  The upper limit in sample size is imposed by the array dimension
!  of the program. In order to allow an easy modification of these
!  limits (which could be convenient for the study of specific cases)
!  they have been introduced as parameters in the aside file
!  MCM.INC, which is included to the source file during the compi-
!  lation.These upper limits are called :
!       - NUMAX : maximal number of cells along the U axis
!       - NVMAX : maximal number of cells along the V axis
!       - NWMAX : maximal number of cells along the W axis
!       - NAMAX : maximal number of sites per cell
!       - NNMAX : maximal number of neighbours per site

!  The file MCM.INC contains the following parameters:
!     PARAMETER(
!    1          NUMAX=20,      ! maximal number of cells along the U axis
!    1          NVMAX=20,      ! maximal number of cells along the U axis
!    1          NWMAX=20,      ! maximal number of cells along the U axis
!    1          NAMAX=100,    ! maximal number of sites per cell
!    1          NNMAX=50,     ! maximal number of neighbours per site
!    1          NMOMAX=200,   ! maximal number of isotropic coupling
!    1                        !   constants which can be modified
!    1                        !   interactively
!    1          NAL=8,        ! maximal number of site moments to be
!    1                        !   written on a single line of the output
!    1                        !   file MYFILE.LIS
!    1          I01=31,       ! input unit number for file MYFILE.MCM
!    1          I02=32,       ! input unit number for file MYFILE.INI
!    1          I03=33,       ! output unit number for file MYFILE.RES
!    1          I04=34,       ! output unit number for file MYFILE.SPI
!    1          IBA=35,       ! output unit number for file MYFILE.REC
!    1          LPT=7,        ! output unit number for file MYFILE.LIS
!    1          ITTI=5,       ! input unit number for terminal
!    1                        !   (interactive mode)
!    1          ITTO=6        ! output unit number for terminal
!    1                        !   (interactive mode)
!    1                   )

!  The initial configuration of spins is either randomly generated
!  by the program or can be imposed (file MYFILE.INI).
!  During the simulation, new orientations of spins are drawn at
!  random and adopted or rejected (according to Boltzmann statistics)
!  in order to progressively minimize the energy of the system.
!  The Hamiltonian used for the energy calculation includes three
!  terms :
!       - a general two-spin coupling term with facility for anisotropic
!         asymmetric exchange
!       - a single-ion anisotropy term
!       - an applied magnetic field term
!  The energy formula may be expressed as follows:
!  (with the convention that a vector V' represents the transposed
!  form of a vector V)

!       E = - 0.5*Si'[Jij]Sj + ((Di'Si)**2)/|Di| - H'Si

!       where :
!          Si,Sj  are vectors representing two interacting spins
!          [Jij]  is the coupling tensor:
!                                                         Jxx Jxy Jxz
!                 for every <i,j> pair, [Jij] represents  Jyx Jyy Jyz
!                                                         Jzx Jzy Jzz
!          Di     is the vector of single-ion anisotropy
!          H      is the vector of applied magnetic field

!  The reference used for spins is totally independent of the crys-
!  tallographic one: spins are always generated and expressed in an
!  orthonormal coordinate system during the simulation process. The
!  same orthonormal coordinate system is used for the anisotropy
!  coefficients and for the applied magnetic field if any. For sake
!  of convenience, a possibility is offered at the end of the
!  simulation, to project the magnetic moments on the real crystal
!  axes and to rotate the whole system to align one moment along a
!  desired direction of the crystal frame of reference. This option
!  facilitates the comparison between the result of simulation and
!  real or theoretical arrangements of magnetic moments.

!  The external conditions of the simulation (cooling rate, varying
!  field ...) are introduced in an interactive way, as sequences of
!  iteration on temperature at H=cst or on magnetic field at T=cst.
!  The arrangements of magnetic moments and last spin configuration
!  are stored during the process and at the end of simulation, respec-
!  tively. Averaged characteristics such as magnetic susceptibility
!  or specific heat and also energy and magnetization are collected
!  at every step of the iterations.

!  The following files are generated during or at the end of the
!  simulation :

!  - MYFILE.REC : contains the sequence of instructions given by
!                 user during the run of MCMAG (can be used as ins-
!                 truction file for further batch runs)

!  - MYFILE.LIS : contains all informations about the arrangement of
!                 magnetic moments and the averaged characteristics
!                 collected at every step of the simulation process

!  - MYFILE.RES : created at user's request, this file stores averaged
!                 characteristics collected at every temperature or
!                 field step, for further analysis or plotting

!  - MYFILE.SPI : created at the end of simulation, it contains the
!                 last generated spin configuration of the simulation.
!                 The file format is the same as that of file MYFILE.INI,
!                 which allows to split up a simulation in several runs

!       The main program calls the following subroutines :

!  READ_FILE   : read input data and simulation parameters
!  COMP_TEST   : check internal consistency of data
!  MOD_INT     : if activated, modify some selected isotropic coupling
!  GEN_SPI     : generate or read initial spin configuration
!  MC_CYCLE    : execute the Monte Carlo loops
!  EN_TOT      : compute the total energy of the system
!  CONST_FUNCT : compute the constraint function (isotropic coupling)
!  OUT_RES     : output averaged characteristics
!  OUT_CONF    : output the configuration of magnetic moments
!  PRO_ROTA    : if selected, project moments on crystal cell (rotate)

!  These subroutines call other subroutines or function :

!  BOUND_COND : compute the selected boundary conditions
!  ENn_CALC   : compute energy difference between new and old spin conf.
!  RANnD      : generate a pseudo-random vector
!  RANF       : generate a pseudo-random number between 0 and 1

!  Four subroutines ENn_CALC and four subroutines RANnD are included in
!  the program (n = 1, 2, 3 or q). Only one among the four is activated
!  during a simulation, depending on the selected spin model, namely :
!  - Ising spin model           (n = 1)
!  - XY spin model              (n = 2)
!  - Heisenberg spin model      (n = 3)
!  - q-states planar Potts model(n = q)

!***********************************************************************

!  The following lines give a list of the main parameters and arrays
!  used in MCMAG (main program and subroutines) :

!    NAV                   : indices of the neighbours of a site

!    PT                    : pointer of NAV

!    NU,NV,NW              : absolute indices of neighbouring cells

!    NAT                   : indices of sites in a cell

!    NVS                   : total number of neighbours of a site

!    M                     : relative indices of neighbouring cells

!    AA,AD                 : matrices for inversion (to rotate moments)

!    TCO                   : test of internal compatibility
!                            (topological self-consistency)

!    TCOJX,JY,JZ,JXY,JXZ   : test of internal compatibility (interactions)
!    ...JYX,JYZ,JZX,JZY

!    AMPS,AMPS2            : spin amplitude

!    SX,SY,SZ              : spin components

!    SOMX,SOMY,SOMZ        : components of magnetic moment
!    COMX                  : storage of SOMX when SOMX is used for
!                            another purpose

!    TMOM                  : total magnetic moment

!    SOMRX,SOMRY,SOMRZ     : components of magnetic moment after rotation

!    SOMX2,SOMY2,SOMZ2     : variance and standard deviation of magnetic
!        moment

!    SOMRX2,SOMRY2,SOMRZ2  : variance and standard deviation of magnetic
!        moment after rotation

!    ISIGX,ISIGY,ISIGZ     : 1000*(standard deviation of components)
!    ISIGM                 : 1000*(standard deviation of total moment)

!    AJSX,AJSY,AJSZ        : components of the local molecular field

!    AJXX,AJXY,AJXZ        : elements of the exchange tensor
!    AJYX,AJYY,AJYZ
!    AJZX,AJZY,AJZZ

!    AJJX,AJJY,AJJZ        : diagonal elements of the exchange tensor
!                            or storage of isotropic (AJXX,AJYY,AJZZ)
!                            before eventual modifications

!    D                     : anisotropy coefficient

!    DDIX,DDIY,DDIZ        : direction of anisotropy
!    D0                    : normalizing factor for the direction of
!                            anisotropy

!    DX,DY,DZ              : normalized direction of anisotropy


!***********************************************************************
!                  *******************************
!                  *                             *
!                  *        MAIN PROGRAM         *
!                  *                             *
!                  *******************************
      PROGRAM MCMAG

      Use mcm_inc
      implicit none
      integer :: nres,icou,i,j,n,ou3,ier
      real(kind=sp) :: hm0,hi,hki,tci,ta,t_ini,t_fin
      real(kind=dp) :: energy,e_aver,temp
      real(kind=sp),dimension(3) :: pos, spv,spr
      integer :: u,v,w,narg
      character(len=20)  :: atm_text=" ",keyword
      character(len=120) :: line=" "

      narg=COMMAND_ARGUMENT_COUNT()
      if(narg > 0) then
           nfil=" "
           call GET_COMMAND_ARGUMENT(1,nfil)
           i=index(nfil,".")
           if(i /= 0) nfil=nfil(1:i-1)
           batch=.true.
      end if

      nres=0
      Write(itto,"(a)")'  ----------------------------------------------------------------------------'
      Write(itto,"(a)")'                               Program MCMAG'
      Write(itto,"(a)")'           SIMULATION OF MAGNETIC STRUCTURES BY A MONTE-CARLO METHOD'
      Write(itto,"(a)")'           MCMAG: a computer program to simulate magnetic structures'
      Write(itto,"(a)")'          P. Lacorre and J. Pannetier, J. Magn. Magn. Mat. 71 (1987) 63.'
      Write(itto,"(a)")'                   Fortran-90 Version May 2015 (JRC-ILL)'
      Write(itto,"(a)")'            Energy Calculated According to the Following Expression :'
      Write(itto,"(a)")"               E = - 0.5*(Si'[Jij]Sj) + Di(ui'Si)**2 - H'Si "
      Write(itto,"(a)")'  ----------------------------------------------------------------------------'
      Write(itto,"(a)")' '
      Write(itto,"(a)")' '
      if(.not. batch) then
         WRITE(itto,"(a)",advance="no")' => Input file name ? : '
         READ(itti,"(a)") nfil
      end if
      Open(Unit=i01,File=trim(nfil)//'.mcm',Status='OLD',iostat=ier)
      if(ier /= 0) then
        write(*,"(a)") " => Error opening the file: "//trim(nfil)//'.mcm'
        stop
      end if
      READ(i01,"(a)")titl1    ! Read First Title of the Input File
      READ(i01,"(a)")titl2    ! Read Second Title of the Input File
      Open(Unit=lpt,STATUS='unknown',File=trim(nfil)//'.lis')

      CALL read_file()   ! Read the Data File and
                         ! Parameters of the Simulation

      CALL comp_test()   ! Test of Internal Compatibility
                         ! Concerning Topology and Interactions


      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(l_case(line))
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (SpinModel): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line,"spinmodel")
         if(j /= 0) then
           keyword=adjustl(line(j+9:))
           j=index(keyword," ")
           Select Case(keyword(1:j-1))
            Case("ising")
              iran=1
            Case("xy")
              iran=2
            Case("heisenberg")
              iran=3
            Case("q-potts")
              iran=4
              read(keyword(j:),*) nq
            Case default
              write(*,"(a)") " => The keyword SpinModel should followed by: Ising, XY, Heisenberg or q-Potts !"
              stop
           End Select
         else
           write(*,"(a)") " => The keyword: SpinModel, followed by Ising, XY, Heisen, or q-Potts, is needed!"
           stop
         end if
      else
         WRITE(itto,"(a)",advance='no') ' => Ising(1), XY(2) Heisenberg(3) or q-State Potts planar(4) spins?: '
                                        ! ISING, XY, HEISENBERG OR q-STATE
                                        ! POTTS PLANAR SPINS ?
         READ(itti,708)iran
      end if
      IF ((iran < 1).OR.(iran > 4)) then
        write(*,"(a)") " => Model ",iran," not supported!"
        stop
      end if

      IF (iran == 4 .and. .not. batch) THEN            ! q-STATE POTTS PLANAR SPINS
        do
          WRITE(itto,"(a)",advance='no')  ' => q-Value (Integer)?: '
          READ(itti,*) nq
          IF (nq < 2) cycle
          exit
        end do
      END IF

      Select Case(iran)
        Case(1)            ! Ising Spins
           spimod='ISING'
           icou = 0
           DO i=1,na
             IF (d0(i) == 0) THEN
               dx(i)=1.
               icou = icou+1
             END IF
           END DO
           IF (icou == na) iwmx=1

        Case(2)      ! XY Spins
           spimod='XY PLANAR'
           iwmxy=1

        Case(3)           ! Heisenberg Spins
        spimod='HEISENBERG'

        Case(4)       ! q-State Potts Planar Spins
           spimod='-STATE POTTS'
           iwmxy=1
      End Select

      if(.not. batch) then
         Open(Unit=Iba,Status='unknown',FILE=trim(nfil)//'.rec')
         Write(Iba,4003) Nfil,iran
         If (Iran == 4) Then
           Write(Iba,*)Nq
         End If
      End If

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Title): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         titl="Simulation using MCMAG"
         if(trim(keyword) == "title") titl=line(j+1:)
      else
         Write(itto,"(a)",advance='no')  ' => Title of simulation ?: '
         Read(itti,"(a)") titl
         Write(iba,"(a)") trim(titl)
      end if

      titl3 = spimod//' Spins - '//trim(titl)
      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Ncells): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "ncells") then
            read(line(j+1:),*,iostat=ier) lu,lv,lw
            if(ier /= 0) then
              write(*,"(a)") " => The keyword Ncells should be followed by the number of cells along a, b and c!"
              stop
             end if
         else
            write(*,"(a)") " => The keyword Ncells is neede here !"
            stop
         end if
      else
         do
           Write(itto,"(a)")           ' => Size of the sample box along the three axes? '
           Write(itto,"(a,3i4,a)",advance='no') '        ( Maximal values: ',numax,nvmax,nwmax,'):  '  ! SIZE OF THE SAMPLE ?
           Read(itti,*,iostat=ier)lu,lv,lw
           If ((lu > numax).OR.(lv > nvmax).OR.(lw > nwmax) .or. ier/=0) cycle
           If ((lu < 1).OR.(lv < 1).OR.(lw < 1)) Then
             Write(itto,888)
             cycle
           End If
           exit
         end do
         Write(iba,*)lu,lv,lw
      end if
      nma=lu*lv*lw
      nta=na*nma
      ansmod='N'
      If (jcod == 0 .and. .not. batch) Then
        do
           WritE(itto,"(a)",advance='no')' => Modify some interactions in the sample box (Y or N)?: '
           Read(itti,"(a)")ansmod
           If ((ansmod /= 'Y').AND.(ansmod /= 'y').AND.  &
               (ansmod /= 'N').AND.(ansmod /= 'n')) cycle
           exit
        end do
        Write(iba,"(a)")ansmod
        If ((ansmod == 'Y').OR.(ansmod == 'y')) Then

          Call mod_int()    ! Possibility to modify locally
                            ! some interactions (for isotropic exchange only)
        End If
      End If

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (InitConf): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "initconf") then
           ians=adjustl(line(j+1:))
         else
            write(*,"(a)") " => The keyword InitConf (followed by R or I) is neede here !"
            stop
         end if
      else
         do
           WriTE(itto,"(a)",advance='no') ' => Initial spin configuration :'            ! INITIAL SPIN CONFIGURATION ?
           WriTE(itto,"(a,a,a)",advance="no") ' => Random(R) or from the input file ',trim(nfil)//'.ini',' (I)?:  '
           ReaD(itti,"(a)")ians
           If ((ians /= 'R').AND.(ians /= 'r').AND.  &
               (ians /= 'I').AND.(ians /= 'i')) cycle
           exit
         end do
         WriTE(iba,"(a)") ians
      end if

      IF ((ians == 'I') .OR. (ians == 'i')) THEN
        OPEN(UNIT=i02,STATUS='OLD',FILE=trim(nfil)//'.ini',iostat=ier)
        if(ier /= 0) then
          write(*,"(a)") " => Error opening the file "//trim(nfil)//'.ini'
          stop
        end if
        READ(i02,2222)itit1,itit2,itit3
      END IF

      CALL gen_spi()  ! Read or Generate the Initial
                      ! Configuration of Spins

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Boundary): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "boundary") then
           bcond=adjustl(l_case(line(j+1:)))
           if(index(bcond,"free") /= 0) ianb=0
           if(index(bcond,"peri") /= 0) ianb=1
           if(index(bcond,"mixe") /= 0) ianb=-1
         else
            write(*,"(a)") " => The keyword Boundary (followed by Free, Periodic or Mixed) is neede here !"
            stop
         end if
      else
         Write(itto,"(a)",advance='no') ' => Free edges(0), Periodic(1) or mixed boundary conditions(-1)?: '
         Read(itti,889,iostat=ier) ianb
         Write(iba,889)ianb
      end if

      IF (ianb < 0) THEN            ! MIXED BOUNDARY CONDITIONS
        bcond = 'MIXED'

      ELSE IF (ianb == 0) THEN       ! FREE EDGES
        bcond = 'FREE EDGES'

      ELSE IF (ianb > 0) THEN       ! PERIODIC BOUNDARY CONDITIONS
        bcond = 'PERIODIC'

      END IF

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Scale): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "scale") then
           avchar=adjustl(l_case(line(j+1:)))
           infa=2
           if(index(avchar,"samp") /= 0) infa=1
           if(index(avchar,"cell") /= 0) infa=2
           if(index(avchar,"site") /= 0) infa=3
           if(index(avchar,"mole") /= 0) infa=4
         else
            write(*,"(a)") " => The keyword Scale (followed by Samp, Cell, Site or Mole) is neede here !"
            stop
         end if
      else
        do
          Write(itto,"(a)") ' => Scaling of averaged characteristics of sample ? :'
          Write(itto,"(a)") '          1 = per sample'
          Write(itto,"(a)") '          2 = per basic unit cell (input file)'
          Write(itto,"(a)") '          3 = per site'
          Write(itto,"(a)") '          4 = per mole'
          Write(itto,"(a)",advance='no') ' => Option: '
          Read(itti,708,iostat=ier) infa
          if(ier /= 0) cycle
          If ((infa < 1).OR.(infa > 4)) cycle
          exit
        end do
        Write(iba,708) infa
      end if


      IF (infa == 1) THEN            ! 1 = PER SAMPLE
        nfa=1
        avchar = 'SAMPLE'
      ELSE IF (infa == 2) THEN       ! 2 = PER UNIT CELL
        nfa=nma
        avchar = 'UNIT CELL'
      ELSE IF (infa == 3) THEN       ! 3 = PER SITE
        nfa=nta
        avchar = 'SITE'
      ELSE IF (infa == 4) THEN       ! 4 = PER MOLE
        nfa=nma*nmot
        avchar = 'MOLE'
      END IF

      CALL en_tot(1, energy)                 ! CALCULATE THE TOTAL ENERGY
! WITH THE CURRENT MAGNETIC FIELD


!***** SUMMARIZE SIMULATION PARAMETERS *****

      WRITE(lpt,596)
      WRITE(lpt,2003)titl
      WRITE(lpt,2015)lu,lv,lw,na,nmot,avchar,trim(bcond)
      IF (nq /= 0) THEN
        WRITE(lpt,2014)nq,spimod
      ELSE
        WRITE(lpt,2016) spimod
      END IF

!*****


      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Sites): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "sites") then
           read(line(j+1:),*,iostat=ier) j1,j2,j3
           if(ier /= 0) then
             j1=1; j2=2; j3=3
           end if
         else
            write(*,"(a)") " => The keyword Sites followed by three integers is neede here !"
            stop
         end if
      else
         WRITE(itto,1262,advance='no')        ! Indices of Three Sites ?
                                              ! (For Output of Resulting Moments)
         READ(itti,*)j1,j2,j3
         WRITE(iba,*)j1,j2,j3

         IF((j1 <= 0).OR.(j1 > na))j1=1
         IF((j2 <= 0).OR.(j2 > na))j2=2
         IF((j3 <= 0).OR.(j3 > na))j3=3
      end if

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Schedule): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "schedule") then
           read(line(j+1:),*,iostat=ier) t,coef,tfin
           if(ier /= 0) then
            write(*,"(a)") " => Error reading T_ini, Coeff and  T_final !"
            stop
           end if
         else
            write(*,"(a)") " => The keyword Schedule followed by three reals (T_ini, Coeff, T_final) is neede here !"
            stop
         end if
      else
         do
            Write(itto,"(a)") ' => Give initial temperature, HEATING (>1) or COOLING (<1) multiplier'
            Write(itto,"(a)",advance='no') '    and final temperature (in Kelvin) : '

            !Kirkpatrick's Cooling(Heating) Scheme : T(n+1) = COEF*T(n)

            Read(itti,*,iostat=ier)t,coef,tfin
            if(ier /= 0) cycle
            exit
         end do
         Write(iba,*) t,coef,tfin
      end if
      !Calculating the expected number of temperatures
      temp=t
      ntemp=0
      if( coef < 1.0) then
         do
           ntemp=ntemp+1
           temp=temp*coef
           if(temp < tfin) exit
         end do
      else
         do
           ntemp=ntemp+1
           temp=temp*coef
           if(temp > tfin) exit
         end do
      end if
      write(lpt,"(a,i5)") "  => Expected number of temperatures: ",ntemp
      iwrt=ntemp
      ansh = 'T'     ! Ansh=T => Iteration On Temperature
      cod  = 'M'     ! Cod=M => Multiplicative Coefficient
      codt = 'x'

  706 if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Hfield): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "hfield") then
           read(line(j+1:),*,iostat=ier) hk,hx0,hy0,hz0
           if(ier /= 0) then
            write(*,"(a)") " => Error reading H, hx,hy,hz !"
            stop
           end if
           If ( (hk == 0).AND.(hx0 == 0).AND. (hy0 == 0).AND.(hz0 == 0) ) hx0 = 1.0
         else
            write(*,"(a)") " => The keyword Hfield followed by four reals (H-Strength, hx,hy,hz) is neede here !"
            stop
         end if
      else
         do
            Write(itto,"(a)",advance='no') ' => Strength and direction of the magnetic field ?: '               ! MAGNETIC FIELD ?
            Read(itti,*,iostat=ier)hk,hx0,hy0,hz0
            if(ier /= 0) cycle
            If ( (hk /= 0).AND.(hx0 == 0).AND. (hy0 == 0).AND.(hz0 == 0) ) cycle
            If ( (hk == 0).AND.(hx0 == 0).AND. (hy0 == 0).AND.(hz0 == 0) ) hx0 = 1.0
            Write(iba,*)hk,hx0,hy0,hz0
            exit
         end do
      end if

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Mcyc): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "mcyc") then
           read(line(j+1:),*,iostat=ier) it,nom
           if(ier /= 0) then
            write(*,"(a)") " => Error reading Number of Montecarlo cycles !"
            stop
           end if
         else
            write(*,"(a)") " => The keyword Mcyc followed by two integers (Num_MCS, Num_Thermal) is neede here !"
            stop
         end if
      else
         do
            WRITE(itto,"(a/a)",advance='no') '   Give : - Total number of MCS/S', &
                                '          - Number of thermalisation MCS/S: '
            READ(itti,*,iostat=ier)it,nom
            if(ier /= 0) cycle
            WRITE(iba,*)it,nom
            exit
         end do
      end if

      snom=FLOAT(it-nom)

      WRITE(lpt,2112)                ! PARAMETERS FOR CURRENT SCHEME
        2112 FORMAT(8X,  &
            ' _______________________________________________________ '  &
            /8X,'|                                                       |'  &
            /8X,'|   SIMULATION PARAMETERS USED FOR THE FOLLOWING RUN    |')
      WRITE(lpt,2013)it,nom,t,tfin,codt,coef,hk,hk
        2113 FORMAT(8X,  &
            '|_______________________________________________________|'  &
            /8X,'|                                                       |'  &
            /8X,'| at any H :',i6,'MCS/S  with',i6,' omitted for average |'  &
            /8X,'|     From T =',f8.3,'  to T =',f8.3,'                  |'  &
            /8X,'|     Temperature program :        none                 |'  &
            /8X,'|     From H =',f9.2,' to H =',f9.2,'  kOe            |'  &
            /8X,'|     Field program       : H(n+1) = H(n) + ',f8.2,'    |'  &
            /8X,'|     Field in the direction  :  ',3(f6.2,1X),'  |'  &
            /8X,'|_______________________________________________________|')
      If (hk /= 0) Then
        Write(lpt,2031)hx0,hy0,hz0
      End If
      Write(lpt,2030)

  243 hm0=SQRT(hx0*hx0+hy0*hy0+hz0*hz0)
      hdix=hx0/hm0
      hdiy=hy0/hm0
      hdiz=hz0/hm0
      hi=hk/14.88732                 ! CONVERSION HK(Kelvin)->HI(kOe)
      hx=hdix*hi                     ! HX,HY,HZ: EFFECTIVE COMPONENTS
      hy=hdiy*hi                     ! OF THE MAGNETIC FIELD
      hz=hdiz*hi
      WRITE(lpt,"(/)")

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#"  .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Print): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "print") then
           read(line(j+1:),*,iostat=ier) ansnev
           if(ier /= 0) ansnev="N"
         else
            write(*,"(a)") " => The keyword Print followed by N, E or A, is neede here !"
            stop
         end if
      else
         do
            Write(itto,"(a)") ' => Printout of the configuration of magnetic moments:'
            Write(itto,"(a)") '     N = never'
            Write(itto,"(a)") '     E = at the final temperature/field only'
            Write(itto,"(a)") '     A = at every temperature/field'
            Write(itto,"(a)",advance='no') ' => Option: '
            Read(itti,"(a)")ansnev
            IF ((ansnev /= 'A').AND.(ansnev /= 'a').AND.  &
                (ansnev /= 'E').AND.(ansnev /= 'e').AND.  &
                (ansnev /= 'N').AND.(ansnev /= 'n')) cycle
            exit
         end do
         Write(iba,"(a)")ansnev
      end if

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#"  .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Averages): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "averages") then
           ansrf='Y'
         else
            ansrf='N'
            backspace(unit=i01)
         end if
      else
         do
           WRITE(itto,"(3a)",advance='no') ' => Output of averaged characteristics in file ',trim(nfil)//'.res',' (Y or N)?: '
           READ(itti,"(a)") ansrf
           IF ((ansrf /= 'Y').AND.(ansrf /= 'y').AND.  &
              (ansrf /= 'N').AND.(ansrf /= 'n')) cycle
           exit
         end do
         Write(iba,"(a)")ansrf
      end if


!***** OPENING AND FIRST LINES OF THE AVERAGED CHARARACTERISTICS FILE *****


      IF ((ansrf == 'Y') .OR. (ansrf == 'y')) THEN
         OPEN(UNIT=i03,FILE=trim(nfil)//'.res',STATUS='unknown')
         ou3=1
         WRITE(itto,"(3a)") ' => File ',trim(nfil)//'.res',' created'
         nres=nres+1
         IF (nres == 1) THEN
           WRITE(i03,"(3a)") '[ Input file : ',trim(nfil)//'.mcm ]  ',titl1
           WRITE(i03,"(a)") "   "//titl2
           WRITE(i03,"(a)") "   "//titl
           IF (nq /= 0) THEN
             WRITE(i03,"(4(a,i2),2a)") "   ",nq,spimod//' Spins - Sample size : ',lu,' X',lv,' X',lw,' - results per ',avchar
           ELSE
             WRITE(i03,"(3(a,i2),2a)") "   "//spimod,lu,' X',lv,' X',lw,' - results per ',avchar
           END IF
         END IF
         WRITE(i03,"(1x,i2,i10,a,i5,a,3(f6.2,2X))")  &
                 nres,it,' (-',nom,' )  MCS/S   with H and/or M// along : ',hx0,hy0,hz0
         WRITE(i03,"(a,3(i2,a))") '       T      H      1/SUS    SP. HT.      EN     MTOT  M(', &
                                  j1,') M(',j2,') M(',j3,')     M//   Fc'
      END IF

      !***** Averaged Characteristics File Opened *****

      call cpu_time(t_ini)
      cpu=t_ini
      do   !Loop over possible different conditions (interactive)
         If (hk /= 0) Then
           Write(lpt,484)t,hk,hi,hx0,hy0,hz0
           Write(itto,484)t,hk,hi,hx0,hy0,hz0
           484 FORMAT(1X,'***********************',/,' **  T =  ',f8.3,  &
               '    **************************************************',/,  &
               ' **  H =',f9.2,' kOe ( ',f7.2,' [K] )  along : ',3(1X,f6.3),  &
               '  **',/,' ************************************************',  &
               '***********************')
         Else
           WRITE(lpt,494)t,hk,hi
           WRITE(itto,494)t,hk,hi
           494 FORMAT(1X,'***********************',/,' **  T =  ',f8.3,  &
               '    *******************',/,' **  H =',f9.2,  &
               ' kOe ( ',f7.2,' [K] )  **',/, ' ****************************************')
         End If

         Write(itto,"(a)") ' => Starting monte-carlo loops'

         Call en_tot(1,energy)          ! Calculate the Total Energy
                                        ! with the Current Magnetic Field
         Call mc_cycle()                ! Execute the Monte-Carlo Loops


         !*****    TEST IF THE CURRENT ITERATION IS THE LAST ONE *****


         IF ((ansh == 'F').OR.(ansh == 'f')) THEN
           hki=hk+hste
           IF ((hste > 0).AND.(hki > hfin)) last=1
           IF ((hste < 0).AND.(hki < hfin)) last=1

         ELSE
           IF ((cod == 'M').OR.(cod == 'm')) THEN
             tci=t*coef
             IF ((coef > 1).AND.(tci > tfin)) last=1
             IF ((coef < 1).AND.(tci < tfin)) last=1
           ELSE
             tci=t+coef
             IF ((coef > 0).AND.(tci > tfin)) last=1
             IF ((coef < 0).AND.(tci < tfin)) last=1
           END IF

         END IF


         !*****    END OF LAST ITERATION TEST *****


         CALL out_conf()         ! Output Configuration Of Moments

         IF (jcod == 0) THEN     ! Calculate the Constraint Function
           CALL const_funct()    ! of the System (For Isotropic
         END IF                  ! Exchange Only)

         CALL out_res()          ! Output Averaged Characteristics

         CALL en_tot(2,energy)   ! Calculate the Total Energy
                                 ! With the Current Magnetic Field
         !*****    INCREMENT FIELD OR TEMPERATURE *****


         IF ((ansh == 'F').OR.(ansh == 'f')) THEN
                                        ! ANSH=F => ITERATION ON FIELD
           hk=hk+hste                   ! HK = NEW MAGNETIC FIELD
           hi=hk/14.88732
           hx=hdix*hi
           hy=hdiy*hi
           hz=hdiz*hi

           IF ((hste > 0).AND.(hk <= hfin)) cycle
           IF ((hste < 0).AND.(hk >= hfin)) cycle

         ELSE                           ! ANSH=T => ITERATION ON TEMPERATURE
           IF ((cod == 'M').OR.(cod == 'm')) THEN
                                        ! COD=M => MULTIPLICATIVE COEFFICIENT
                                        ! KIRKPATRICK'S COOLING(HEATING) SCHEME
             ta=t*coef
             IF ((coef > 1).AND.(ta <= tfin)) THEN
               t=ta
               cycle
             END IF
             IF ((coef < 1).AND.(ta >= tfin)) THEN
               t=ta
               cycle
             END IF
           ELSE                         ! COD=A => ADDITIVE COEFFICIENT
                                        ! LINEAR COOLING(HEATING) SCHEME
             ta=t+coef
             IF ((coef > 0).AND.(ta <= tfin)) THEN
               t=ta
               cycle
             END IF
             IF ((coef < 0).AND.(ta >= tfin)) THEN
               t=ta
               cycle
             END IF
           END IF
         END IF
         exit
      end do

!***** END OF INCREMENTATION *****


      last=0
      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#"  .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Continue): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "continue") then
            ansc='Y'
         else
            ansc='N'
            backspace(unit=i01)
         end if
      else
        do
           WRITE(itto,"(a)",advance='no') ' => Do you want to continue (y or n) ?: '
           READ(itti,"(a)")ansc
           WRITE(lpt,"(a)")
           IF ((ansc /= 'Y').AND.(ansc /= 'y').AND.  &
               (ansc /= 'N').AND.(ansc /= 'n')) cycle
           WRITE(iba,"(a)")ansc
           exit
        end do
      end if
        IF ((ansc == 'Y').OR.(ansc == 'y')) THEN
          do
             WRITE(itto,"(a)",advance='no') ' => Modify temperature (T) or magnetic field (F) ?: '
             READ(itti,"(a)")ansh
             IF ((ansh /= 'F').AND.(ansh /= 'f').AND.  &
                 (ansh /= 'T').AND.(ansh /= 't')) cycle
             WRITE(iba,"(a)")ansh
             exit
          end do
          IF ((ansh == 'T').OR.(ansh == 't')) THEN
            do
               WRITE(itto,"(a)")  ' => Define code for HEATING/COOLING factor :'
               WRITE(itto,"(a)")  '              M = multiply'
               WRITE(itto,"(a)")  '              A = add'
               WRITE(itto,"(a)",advance='no')  ' => Option: '
               READ(itti,"(a)")cod
               IF ((cod /= 'M').AND.(cod /= 'm').AND.  &
                   (cod /= 'A').AND.(cod /= 'a')) cycle
               WRITE(iba,"(a)")cod
               exit
            end do
            do
              WRITE(itto,"(a)") ' => Heating/cooling factor and final temperature ?: '

              !   Kirkpatrick's Cooling(Heating) Scheme : T(n+1) = COEF*T(n)
              !   Linear Cooling(Heating) Scheme : T(n+1) = COEF+T(n)

              READ(itti,*,iostat=ier)coef,tfin
              if(ier /= 0) cycle
              WRITE(iba,*)coef,tfin
              exit
            end do

            IF ((cod == 'M').OR.(cod == 'm')) THEN
              t=t*coef     ! FIRST TEMPERATURE OF THE NEW ITERATION
              codt = 'x'

            ELSE
              t=t+coef     ! FIRST TEMPERATURE OF THE NEW ITERATION
              codt = '+'

            END IF

            GO TO 706

          END IF

   1116   WRITE(itto,"(a)",advance='no') ' => Orientation of the magnetic field h ?: '
          READ(itti,*,ERR=1116)hx0,hy0,hz0
          IF ( (hx0 == 0) .AND. (hy0 == 0).AND.(hz0 == 0) ) GO TO 1116
          if(.not. batch) WRITE(iba,*) hx0,hy0,hz0
   1117   WRITE(itto,"(a)",advance='no')   ' => Starting, increment and final fields ?: '
          READ(itti,*,ERR=1117)hsta,hste,hfin
          if(.not. batch) WRITE(iba,*)hsta,hste,hfin
          hk=hsta

    519   WRITE(itto,"(a/a)",advance='no') '   GIVE : - TOTAL NUMBER OF MCS/S', &
                                           '   -     NUMBER OF THERMALISATION MCS/S: '
          READ(itti,*,ERR=519)it,nom
          if(.not. batch) WRITE(iba,*)it,nom
          snom=FLOAT(it-nom)

          WRITE(lpt,2112)              ! PARAMETERS FOR CURRENT SCHEME
          WRITE(lpt,2113)it,nom,t,t,hsta,hfin,hste,hx0,hy0,hz0

          GO TO 243

        END IF

        IF ((ansrf == 'Y').OR.(ansrf == 'y')) THEN
          WRITE(i03,509)
        END IF

      if(batch) then
         do
           read(i01,"(a)",iostat=ier) line
           if(ier == 0) then
             write(lpt,"(a)") "    mcm-line => "//trim(line)
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#"  .or. len_trim(line) == 0) cycle
             exit
           else
             write(*,"(a)") " => Premature end of file (Cryst): "//trim(nfil)//".mcm"
             stop
           end if
         end do       !123456789
         j=index(line," ")
         keyword=l_case(line(1:j-1))
         if(trim(keyword) == "cryst") then
            ansm='Y'
            read(unit=line(j+1:),fmt=*,iostat=ier)  nor,dorx,dory,dorz
            if(ier /= 0) then
              write(*,"(a)") " => Error reading site and orienting direction after keyword CRYST ! "
              stop
            end if
         else
            ansm='N'
            backspace(unit=i01)
         end if
      else
          do
             Write(itto,"(a)") ' => Do you want a printout of the coomponents of magnetic moments'
             Write(itto,"(a)",advance='no') '    in the Crystallographic System of Coordinates ?: '
             Read(itti,"(a)")ansm
             If ((ansm /= 'Y').AND.(ansm /= 'y').AND.  &
                 (ansm /= 'N').AND.(ansm /= 'n'))cycle
             exit
          end do
          Write(iba,"(a)")ansm
      end if
        IF ((ansm == 'Y').OR.(ansm == 'y')) THEN
          CALL pro_rota()              ! Projection and Rotation of
          rotated=.true.               ! Moments Inside the Crystal
                                       ! Coordinates System
        END IF
!***** OUTPUT OF THE LAST GENERATED SPIN CONFIGURATION (ALSO IN FST format) *****


        OPEN(UNIT=i04,STATUS='unknown',FILE=trim(nfil)//'.spi' )
        WRITE(i04,"(a)")titl1
        WRITE(i04,"(a)")titl2

        IF (nq /= 0) THEN
          WRITE(i04,"(i3,a)")nq," "//titl3
        ELSE
          WRITE(i04,"(a)")titl3
        END IF

        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              WRITE(i04,798)
              WRITE(i04,*)u,v,w
              WRITE(i04,799)
              DO n=1,na
                WRITE(i04,*)n,sx(n,u,v,w),sy(n,u,v,w),sz(n,u,v,w)
              END DO
            END DO
          END DO
        END DO
        close(unit=i04)
        if(coordinates_read) then
           OPEN(UNIT=i04,STATUS='replace',FILE=trim(nfil)//'.fst' )
           WRITE(i04,"(a)")"!  File generated by MCMAG"
           WRITE(i04,"(a)")"TITLE "//trim(titl1)//"  "//trim(titl2)
           WRITE(i04,"(a)")"SPACEG P 1 "
           WRITE(i04,"(a,3f12.6,3f12.4)")"CELL ",lu*aax,lv*aay,lw*aaz,aaa,aab,aag
           WRITE(i04,"(a)")"  "
           WRITE(i04,"(a)") 'BOX   -0.15  1.15   -0.15  1.15    -0.15  1.15'
           DO w=1,lw
             DO v=1,lv
               DO u=1,lu
                 DO n=1,na
                   atm_text=" "
                   write(atm_text,"(a,3i2.2)") trim(namt(n))//"_",u-1,v-1,w-1
                   pos=(real([u-1,v-1,w-1])+xyz(:,n))/real([lu,lv,lw])
                   write(i04,"(a,3f10.5,a)") "Atom "//trim(atm_text)//"  "//elem(n),pos
                 END DO
               END DO
             END DO
           END DO
           WRITE(i04,"(a)")"  "
           WRITE(i04,"(a)")'Mag_Structure'
           WRITE(i04,"(a)")'Lattice P'
           WRITE(i04,"(a)")'Kvect    0.00   0.00   0.00'
           WRITE(i04,"(a)")'SYMM   x, y, z'
           WRITE(i04,"(a)")'MSYM   u, v, w, 0.00'
           DO w=1,lw
             DO v=1,lv
               DO u=1,lu
                 DO n=1,na
                   atm_text=" "
                   spv=[sx(n,u,v,w),sy(n,u,v,w),sz(n,u,v,w)] !Spins of the last configuration in Cartesian Coordinates
                   spr=matmul(ac,spv)                        !Converted to unitary system along a,b,c
                   write(atm_text,"(a,3i2.2)") trim(namt(n))//"_",u-1,v-1,w-1
                   pos=(real([u-1,v-1,w-1])+xyz(:,n))/real([lu,lv,lw])
                   write(i04,"(a,3f10.5,a)") "Matom "//trim(atm_text)//"  "//elem(n),pos,'  SCALE 1.0 GROUP'
                   write(i04,"(a20,3f10.5,a)") "SkP    1   1 ",spr,'  0.00  0.00  0.00    0.00  '
                 END DO
               END DO
             END DO
           END DO
           WRITE(i04,"(a)")'End_Mag_Structure'
           close(unit=i04)

           if(rotated) then
             !Write another fst-file with rotated moments
             OPEN(UNIT=i04,STATUS='replace',FILE=trim(nfil)//'.cfl' )
             WRITE(i04,"(a,i3,a,3f8.3,a)")"!  File generated by MCMAG: Moments rotated by aligning moment: ",&
                                         nor, " along [",dorix,doriy,doriz," ]"
             WRITE(i04,"(a)")"Title "//trim(titl1)//"  "//trim(titl2)
             WRITE(i04,"(a)")"Spaceg P 1 "
             WRITE(i04,"(a,3f12.6,3f12.4)")"Cell ",lu*aax,lv*aay,lw*aaz,aaa,aab,aag
             WRITE(i04,"(a)")"  "
             WRITE(i04,"(a)") 'Box   -0.15  1.15   -0.15  1.15    -0.15  1.15'
             WRITE(i04,"(a)")"  "
             WRITE(i04,"(a)")'Mag_Structure'
             WRITE(i04,"(a)")'Lattice P'
             WRITE(i04,"(a)")'Kvect    0.00   0.00   0.00'
             WRITE(i04,"(a)")'Symm   x, y, z'
             WRITE(i04,"(a)")'Msym   u, v, w, 0.00'
             DO w=1,lw
               DO v=1,lv
                 DO u=1,lu
                   DO n=1,na
                     atm_text=" "  !(COMX,SOMY,SOMZ)<- output of pro_rota
                     !spr=[comx(n,u,v,w),somy(n,u,v,w),somz(n,u,v,w)] !Average Spins in Crystallographic Coordinates
                     spr=[somrx(n,u,v,w),somry(n,u,v,w),somrz(n,u,v,w)] !Average Spins in Crystallographic Coordinates
                     write(atm_text,"(a,3i2.2)") trim(namt(n))//"_",u-1,v-1,w-1
                     pos=(real([u-1,v-1,w-1])+xyz(:,n))/real([lu,lv,lw])
                     write(i04,"(a,3f10.5,a)") "Matom "//trim(atm_text)//"  "//scattf(n),pos,'  SCALE 1.0 GROUP'
                     write(i04,"(a20,3f10.5,a)") "SkP    1   1 ",spr,'  0.00  0.00  0.00    0.00  '
                     spv=matmul(aa,spr)  !Cartesian components
                     sx(n,u,v,w)=spv(1)  !Copy the average Cartesian moments into spin variables
                     sy(n,u,v,w)=spv(2)
                     sz(n,u,v,w)=spv(3)
                     !write(*,"(4i4,3f10.5,tr5,3f10.5)") n,u,v,w,spr,spv
                   END DO
                 END DO
               END DO
             END DO
             WRITE(i04,"(a)")'End_Mag_Structure'
             close(unit=i04)
           else    !!Copy the average moments into spin variables even without rotation
             DO w=1,lw
               DO v=1,lv
                 DO u=1,lu
                   DO n=1,na
                     spr=[comx(n,u,v,w),somy(n,u,v,w),somz(n,u,v,w)] !Moments of the last cycle in Crystallographic Coordinates
                     spv=matmul(aa,spr)  !Cartesian components
                     sx(n,u,v,w)=spv(1)  !Copy the average moments into spin variables
                     sy(n,u,v,w)=spv(2)  !Before calculating the total energy
                     sz(n,u,v,w)=spv(3)
                   END DO
                 END DO
               END DO
             END DO
           end if !Rotated
           call en_tot(3,e_aver)
           delta=e_aver-energy

           write(lpt,"(a/a,f12.5/)") " => Difference betwen average energy and the energy", &
                          "    of the last spin configuration:",delta
           write(*,"(a/a,f12.5/)") " => Difference betwen average energy and the energy", &
                          "    of the last spin configuration:",delta
        end if

! ***** OUTPUT COMMENTS AT THE END OF THE PROGRAM *****
      call cpu_time(t_fin)
      write(unit=*,fmt="(a,f12.4,a)")  " => CPU time: ",t_fin-t_ini, " seconds"


      write(itto,"(a)") " => The following files have been created: "
      write(itto,"(a)") "    "//trim(nfil)//'.rec    with the sequence of instruction you gave'
      write(itto,"(a)") "    "//trim(nfil)//'.lis    with all the results of simulation'
      write(itto,"(a)") "    "//trim(nfil)//'.spi    with the last generated spin configuration (Cartesian)'
      if(coordinates_read) then
         write(itto,"(a)") "    "//trim(nfil)//'.fst    with the last generated spin configuration (Crystallographic)'
        if(rotated) write(itto,"(a)") "    "//trim(nfil)//'.cfl    with the average moments after rotations (Crystallographic)'
      end if

      write(lpt,"(a)") " => The following files have been created: "
      write(lpt,"(a)") "    "//trim(nfil)//'.rec    with the sequence of instruction you gave'
      write(lpt,"(a)") "    "//trim(nfil)//'.lis    with all the results of simulation'
      write(lpt,"(a)") "    "//trim(nfil)//'.spi    with the last generated spin configuration (Cartesian)'
      if(coordinates_read) then
         write(lpt,"(a)") "    "//trim(nfil)//'.fst    with the last generated spin configuration (Crystallographic)'
        if(rotated) write(lpt,"(a)") "    "//trim(nfil)//'.cfl    with the average moments after rotations (Crystallographic)'
      end if

      IF(ou3 /= 0) THEN
        write(lpt,"(a)")  "    "//trim(nfil)//'.res    with condensed list of averaged characteristics'
        write(itto,"(a)") "    "//trim(nfil)//'.res    with condensed list of averaged characteristics'
      END IF

      WRITE(lpt,597)                 ! END COMMENTS
      STOP


      2 FORMAT(a6)
      2003 FORMAT(1X,' ------------------------------------------',  &
          '----------------------------------',// ,12X,a50,//,  &
          1X,' ------------------------------------------',  &
          '----------------------------------')
      4003 FORMAT(a6,/,'Y',/,i1)
      2222 FORMAT(a80/a80/a80)
      2022 FORMAT(2X,a80)
      888 FORMAT(' => Sample size must be greater than or equal to 1 along any axis')
        708 FORMAT(i1)
        889 FORMAT(i2)
        1262 FORMAT(' => Give the indices of 3 sites whose moments will be collected: ' )
        2015 FORMAT(10X,'  _________________________________________________'  &
            ,/,11X,'|                                                 |'  &
            ,/,11X,'|        GENERAL PARAMETERS OF SIMULATION         |'  &
            ,/,11X,'|_________________________________________________|'  &
            ,/,11X,'|                                                 |'  &
            ,/,11X,'|  SAMPLE SIZE    :    ',i2,' x',i2,' x',i2,'  CELLS'  &
            ,10X,'|',/,11X,'|  with  ',i2,' SITES per CELL and ',i2,' MOLES per CELL  |'  &
            ,/,11X,'|  AVERAGED CHARACTERISTICS   :    per ',a9,  &
            '  |',/,11X,'|  BOUNDARY CONDITIONS        :     ',a10,'    |')
        2014 FORMAT(11X,'|  SPIN MODEL     :     ',i2,a12,'            |'  &
            ,/,11X,'|_________________________________________________|'//)
        2016 FORMAT(11X,'|  SPIN MODEL     :        ',a12,'           |'  &
            ,/,11X,'|_________________________________________________|'//)
        2013 FORMAT(8X,  &
            '|_______________________________________________________|'  &
            /8X,'|                                                       |'  &
            /8X,'| at every T :',i6,'MCS/S  (',i6,' omitted for average) |'  &
            /8X,'|     From T =',f8.3,'  to T =',f8.3,'                  |'  &
            /8X,'|     Temperature program : T(n+1) = T(n) ',a1,f8.3,'     |'  &
            /8X,'|     From H =',f9.2,' to H =',f9.2,'  kOe            |'  &
            /8X,'|     Field program       :        none                 |')
        2031 FORMAT(8X, '|     Field in the direction  :  ',3(f6.2,1X),'  |')
        2030 FORMAT(8X, '|_______________________________________________________|')
        2018 FORMAT(3X,i3,a77)
        1018 FORMAT(1X,'!! Temperature cannot be less than or eq. to zero !!')
          509 FORMAT(' 99')
            798 FORMAT(1X,'CELL')
            799 FORMAT(9X,'Site',5X,'Sx',13X,'Sy',13X,'Sz')
            3000 FORMAT(1X,' ------------------------------------------'  &
                ,'----------------------------------')
            596 FORMAT(  '  ****************************** START SIMULATION *******************************',//)
          597 FORMAT(///,'  ************************************ END **************************************')
          595 FORMAT(1X)

      END PROGRAM MCMAG

      subroutine write_fst(ier)
        Use mcm_inc
        implicit none
        integer,intent(out):: ier
        integer:: n,u,v,w,i_fst=34
        character(len=20) :: atm_text
        real(kind=sp), dimension(3) :: spr,spv,pos
        if(coordinates_read) then
           OPEN(unit=i_fst,FILE='simann.fst',STATUS='replace',action="write",iostat=ier )
           if(ier /= 0) return
           WRITE(i_fst,"(a)")"!  File generated by MCMAG"
           WRITE(i_fst,"(a)")"TITLE "//trim(titl1)//"  "//trim(titl2)
           WRITE(i_fst,"(a)")"SPACEG P 1 "
           WRITE(i_fst,"(a,3f12.6,3f12.4)")"CELL ",lu*aax,lv*aay,lw*aaz,aaa,aab,aag
           WRITE(i_fst,"(a)")"  "
           WRITE(i_fst,"(a)") 'BOX   -0.15  1.15   -0.15  1.15    -0.15  1.15'
           WRITE(i_fst,"(a)")"  "
           WRITE(i_fst,"(a)")'{'
           WRITE(i_fst,"(a)")'Lattice P'
           WRITE(i_fst,"(a)")'K    0.00   0.00   0.00'
           WRITE(i_fst,"(a)")'SYMM   x, y, z'
           WRITE(i_fst,"(a)")'MSYM   u, v, w, 0.00'
           DO w=1,lw
             DO v=1,lv
               DO u=1,lu
                 DO n=1,na
                   atm_text=" "
                   spv=[sx(n,u,v,w),sy(n,u,v,w),sz(n,u,v,w)] !Spins of the last configuration in Cartesian Coordinates
                   spr=matmul(ac,spv)                        !Converted to unitary system along a,b,c
                   write(atm_text,"(a,3i2.2)") trim(namt(n))//"_",u-1,v-1,w-1
                   pos=(real([u-1,v-1,w-1])+xyz(:,n))/real([lu,lv,lw])
                   write(i_fst,"(a,3f10.5,a)") "Matom "//trim(atm_text)//"  "//elem(n),pos,'  SCALE 1.0 GROUP'
                   write(i_fst,"(a20,3f10.5,a)") "SkP    1   1 ",spr,'  0.00  0.00  0.00    0.00  '
                 END DO
               END DO
             END DO
           END DO
           WRITE(i_fst,"(a)")'}'
           call flush(i_fst)
           close(unit=i_fst)
        end if
      end subroutine write_fst


      SUBROUTINE read_file()

!***********************************************************************

!  This subroutine reads the input file MYFILE.DAT which contains the
!  topological and magnetic parameters.

!  This subroutine calls no subroutine.

!***********************************************************************


      Use mcm_inc
        implicit none
        Integer :: nin,nfi,j,i,ier,naco,k
        real(kind=sp) :: smo
        Character (Len=34)  :: Nfila
        Character (Len=132) :: Line
        Character (Len=6)   :: aux_scf

        naco=0
        Write(lpt,"(a)")'  ----------------------------------------------------------------------------'
        Write(lpt,"(a)")'                               Program MCMAG'
        Write(lpt,"(a)")'           SIMULATION OF MAGNETIC STRUCTURES BY A MONTE-CARLO METHOD'
        Write(lpt,"(a)")'           MCMAG: a computer program to simulate magnetic structures'
        Write(lpt,"(a)")'          P. Lacorre and J. Pannetier, J. Magn. Magn. Mat. 71 (1987) 63.'
        Write(lpt,"(a)")'                   Fortran-90 Version May 2015 (JRC-ILL)'
        Write(lpt,"(a)")'            Energy Calculated According to the Following Expression :'
        Write(lpt,"(a)")"               E = - 0.5*(Si'[Jij]Sj) + Di(ui'Si)**2 - H'Si "
        Write(lpt,"(a)")'  ----------------------------------------------------------------------------'
        call Write_Date_Time(lpt)
        WRITE(lpt,"(//a//)") ' Input file : '//trim(nfil)//'.mcm  contains the following data : '
        WRITE(*,*) " => Input File: ",trim(nfil)//".mcm"
        WRITE(lpt,"(1x,a)")titl1
        WRITE(lpt,"(1x,a)")titl2
        WRITE(lpt,"(a)") " "
        nfila = trim(nfil)//'.mcm'
        do
          read(i01,"(a)") line
          line=adjustl(line)
          if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
          read(line,*,iostat=ier) na,jcod,nmot
          if(ier /= 0) then
            write(*,"(a)") " => Error reading the MCM file at the items NA,JCOD,NV !"
            write(*,"(a)") " => Error in Line: "//trim(line)
            stop
          end if
          exit
        end do

        if (((jcod /= 0).and.(jcod /= 1).and.(jcod /= 2))  &
            .or.(na == 0).or.(nmot < 1))   then
            write(*,"(a/a)") " =>  NA and NMOT must be integers and at least equal to 1",&
                             ' =>  JCOD must be equal to 0, 1 or 2'
            stop
        end if

        iani=1

        IF (na < 0) THEN
          iani=-1
          na=-na
        END IF


        DO i=1,na
          do
             read(i01,"(a)",iostat=ier) line
             line=adjustl(line)
             if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
             exit
          end do
          IF (iani == -1) THEN
            READ(line,*,iostat=ier)nat(i),nvs(nat(i))
            if(ier /= 0) nat(i)=0
          ELSE
            READ(line,*,iostat=ier) nat(i),nvs(nat(i)),d(nat(i)),  &
                ddix(nat(i)),ddiy(nat(i)),ddiz(nat(i))
            if(ier /= 0) nat(i)=0
          END IF
          if ((nat(i) < 1).OR.(nvs(nat(i)) < 1).OR. (nat(i) > na))  then
               write(*,"(a/a/a)")' => NAT must be an integer between 1 and NA',  &
                                 '    NV must be an integer at least equal to 1 ',&
                                 "    or an error has been produced in reading single ion anisotropy"
               write(*,"(a)") " => Error in Line: "//trim(line)
               stop
          end if
          !Reading the name and coordinates of the atoms
          j=index(line,"::")
          if(j /= 0) then
             line=adjustl(line(j+2:))
             j=index(line," ")
             namt(i)=line(1:j)
             elem(i)=line(1:2)
             if(index("0123456789_-()",line(2:2)) /= 0) elem(i)(2:2)=" "
             read(line(j+1:),*,iostat=ier) xyz(:,i)
             if(ier /= 0) then
                coordinates_read=.false.
             else
                coordinates_read=.true.
                write(lpt,"(a,i3,a,3f10.5)") "    Atom #",i,"   "//namt(i)//"  "//elem(i),xyz(:,i)
             end if
          end if

          write(lpt,"(/,1X,'SITE ',i2,' [ D = ',f8.3,' along the direction ',3(f6.3,2X),']',/,9X,'with ',i2,' neighbours :')") &
              nat(i),d(nat(i)), ddix(nat(i)),ddiy(nat(i)),ddiz(nat(i)), nvs(nat(i))

          IF (jcod == 0) THEN          ! JCOD=0 => ISOTROPIC EXCHANGE
            DO j=1,nvs(nat(i))
              do
                 read(i01,"(a)",iostat=ier) line
                 line=adjustl(line)
                 if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
                 exit
              end do
              read(line,*,iostat=ier) nav(nat(i),j),(m(nat(i),j,k),k=1,3), ajjx(nat(i),j)
              if(ier /= 0 .or. (nav(nat(i),j) < 1).OR.(nav(nat(i),j) > na) ) then
                 write(*,"(a/a/a)")' => NAV must be an integer between 1 and NA',  &
                                   "    or an error is produced at reading a connetivity line: NAV,Av,Bv,Cv,J"
                 write(*,"(a)") " => Error in Line: "//trim(line)
                 stop
              end if

              write(lpt,"(9X,'Site ',i2,'  in cell (',i2,',',i2,',',i2,')','  with J(',i2,',',i2,') =',f10.3,' K')") &
                    nav(nat(i),j),(m(nat(i),j,k),k=1,3), nat(i),nav(nat(i),j),ajjx(nat(i),j)
              ajjy(nat(i),j)=ajjx(nat(i),j)
              ajjz(nat(i),j)=ajjx(nat(i),j)
            END DO


          ELSE IF (jcod == 1) THEN     ! JCOD=1 => DIAGONAL TENSOR
            DO j=1,nvs(nat(i))
              do
                 read(i01,"(a)",iostat=ier) line
                 line=adjustl(line)
                 if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
                 exit
              end do
              read(line,*,iostat=ier)nav(nat(i),j),(m(nat(i),j,k),k=1,3),  &
                  ajjx(nat(i),j),ajjy(nat(i),j),ajjz(nat(i),j)
              if (ier /= 0 .or. nav(nat(i),j) < 1) then
                 write(*,"(a/a/a)")' => NAV must be an integer between 1 and NA',  &
                                   "    or an error is produced at reading a connectivity line: NAV,Av,Bv,Cv,Jx,Jy,Jz"
                 write(*,"(a)") " => Error in Line: "//trim(line)
                 stop
              end if

              WRITE(lpt,92)nav(nat(i),j),(m(nat(i),j,k),k=1,3),  &
                  nat(i),nav(nat(i),j), ajjx(nat(i),j),ajjy(nat(i),j),ajjz(nat(i),j)
            END DO


          ELSE IF (jcod == 2) THEN     ! JCOD=2 => ANISOTROPIC TENSOR
            DO j=1,nvs(nat(i))
              do
                 read(i01,"(a)",iostat=ier) line
                 line=adjustl(line)
                 if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
                 exit
              end do
              read(line,*,iostat=ier)nav(nat(i),j),(m(nat(i),j,k),k=1,3),  &
                  ajjx(nat(i),j),ajjy(nat(i),j),ajjz(nat(i),j),  &
                  ajxy(nat(i),j),ajxz(nat(i),j),ajyz(nat(i),j),  &
                  ajyx(nat(i),j),ajzx(nat(i),j),ajzy(nat(i),j)
              IF (ier/= 0 .or. nav(nat(i),j) < 1) then
                 write(*,"(a/a/a)")' => NAV must be an integer between 1 and NA',  &
                                   "    or an error is produced at reading a connectivity line:", &
                                   "    NAV,Av,Bv,Cv,Jxx,Jyy,Jzz,Jxy,Jxz,Jyz,Jyx,Jzx,Jzy"
                 write(*,"(a)") " => Error in Line: "//trim(line)
                 stop
              end if

              write(lpt,192)nav(nat(i),j),(m(nat(i),j,k),k=1,3),  &
                  nat(i),nav(nat(i),j), ajjx(nat(i),j),ajxy(nat(i),j),ajxz(nat(i),j),  &
                  ajyx(nat(i),j),ajjy(nat(i),j),ajyz(nat(i),j),  &
                  ajzx(nat(i),j),ajzy(nat(i),j),ajjz(nat(i),j)
            END DO

          END IF

        END DO

        WRITE(lpt,"(a)")" "
        WRITE(lpt,"(a)")" "
        do
           do
              read(i01,"(a)",iostat=ier) line
              line=adjustl(line)
              if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
              exit
           end do
           read(line,*,iostat=ier) nin,nfi,smo
           i=index(line,"M")
           if(i == 0) then
              i=index(line,"J")
              if(i == 0) then
                aux_scf=" "
              else
                aux_scf=line(i:)
              end if
           else
             aux_scf=line(i:)
           end if
           if (ier/= 0 .or. (nin < 1).OR.(nfi < 1).OR.(nfi < nin).OR. (smo < 0.0)) then
              write(*,"(5(a,/))")" => NI and NF must be integers and at least equal to 1",  &
                                 "    NF must always be larger than or equal to NI", &
                                 "    The largest NF must be equal to NAT", &
                                 "    A spin amplitude must be attributed to any site", &
                                 "    Spin amplitudes must be positive or equal to 0"
              write(*,"(a)") " => Error in Line: "//trim(line)
              stop
           end if
           naco = naco+nfi-nin+1

           DO i=nin,nfi
             amps(i)=smo
             if(len_trim(aux_scf) /= 0 ) then
                scattf(i)= aux_scf
             else
                scattf(i)="M"//u_case(elem(i))//"3"
             end if
             write(lpt,"(22x,a,i2,a,f6.3)")' Spin amplitude ( ',i,' ) = ',amps(i)
          END DO

           IF (naco /= na) cycle
           exit
        end do

        DO i=1,na
          d0(i)=SQRT(ddix(i)*ddix(i)+ddiy(i)*ddiy(i)+ddiz(i)*ddiz(i))
          IF(d0(i) /= 0) THEN
            dx(i)=ddix(i)/d0(i)
            dy(i)=ddiy(i)/d0(i)
            dz(i)=ddiz(i)/d0(i)
          END IF
        END DO

        !Read a,b,c,alpha,beta,gamma
        do
           read(i01,"(a)",iostat=ier) line
           line=adjustl(line)
           if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
           exit
        end do
        read(line,*,iostat=ier) aax,aay,aaz,aaa,aab,aag
        if(ier /= 0) then
            coordinates_read=.false. !fst file cannot be written
            write(*,"(a)") " => Error reading the cell parameters"
            write(*,"(a)") " => Error in Line: "//trim(line)
            stop
        end if
        ! * * * COMPUTATION OF TRANSFORMATION MATRICES AA(I,J), AC(I,J) AND AD(I,J) * * *
        !     Vc= AD Vx  <= crystallographic directions transformed to Cartesian direction
        !     Mx= AC Mc  <= Cartesian moments to crystallographic unitary system (AC=AA^-1)
        call conv_matrix()   !Calculate the matrix AC for converting moment to crystallographic components in any case
        !Printing the matrices passing from orthogonal to crystallographic frames and viceversa
        write(lpt,"(/,a)")  "        Matrices for conversions Cartesian <-> Crystallographic "
        write(lpt,"(/,a)")  "    M(Cart.) = AA . M(Cryst.)               M(Cryst.)= AC .  M(Cart.)"
        do i=1,3
          write(lpt,"(3f10.5,tr10,3f10.5)")  aa(i,:),ac(i,:)
        end do
        write(lpt,"(a)")" "


        WRITE(itto,"(a)") " "
        92 FORMAT(9X,'Site ',i2,'  in cell (',i2,',',i2,',',i2,')',  &
            '  with coupling (',i2,',',i2,') :',/,15X,  &
            'Jx = ',f8.3,' K , Jy = ',f8.3,' K , Jz = ',f8.3,' K')
        192 FORMAT(9X,'Site ',i2,'  in cell (',i2,',',i2,',',i2,')',/,  &
            15X,'with J tensor (',i2,',',i2,') =   ',3(f8.3,1X),/,  &
            41X,3(f8.3,1X),/,41X,3(f8.3,1X))
        9876 FORMAT(1X,'D(',i2,') = ',f8.3,2X,'( direction ',3(f6.3,2X), ')')
        RETURN
      END SUBROUTINE read_file


      SUBROUTINE comp_test()
!***********************************************************************

!  This subroutine tests the validity of connectivity between spins.
!  Every interaction is described twice in the input file (spin i with
!  neighbour j, then spin j with neighbour i).
!  The subroutine checks that the two descriptions are compatible
!  as well from the topological viewpoint as from the coupling
!  constant viewpoint.
!  Any time an anomaly is detected, an error message is generated,
!  that gives the location of the corresponding error.

!  The test is designed as follows :
!  - TOPOLOGY
!    For any coupling Jij, the consistency of the cell relative
!    indexing (M(i,j,k), k=1,2,3) of neighbouring spins is checked
!    by calculating TCO(k) = M(i,j,k) + M(j,i,k). If any anomaly is
!    detected (non-zero TCO(k)), an error message is generated, that
!    gives the location of corresponding error (namely M(i,j,k))
!  - INTERACTIONS
!    Following the same scheme, the consistency of interactions Jij
!    is checked by calculating TCOJ* = Jij - Jji for any component
!    of the coupling tensor (* means X,Y,Z,XY,XZ,YX,YZ,ZX or ZY).
!    If any anomaly is detected (non-zero TCOJ*), an error message
!    is generated, that gives the location of corresponding error
!    (namely J(i,j))
!  Remark : For specific purpose, such as the study of non symmetrical
!  coupling (Jij different from Jji), the last part of the test (at
!  least) should be bypassed. If you bypass the test, make sure that
!  your input data are correct !

!  This subroutine calls no subroutine.

!***********************************************************************

      Use mcm_inc

      nh=1
      nhp=1

      DO i=1,na

        DO nj=1,nnmax
          pt(nj)=0
        END DO

   do_j:DO  j=1,nvs(nat(i))
          IF((jcod == 2).AND.(nat(i) == nav(nat(i),j)))GO TO 4446
          DO ni=1,nnmax
            IF ( nav(nat(i),j) == pt(ni) ) cycle do_j
          END DO
          DO k=1,3
            tco(k)=m(nat(i),j,k)
          END DO
          tcojx=ajjx(nat(i),j)
          IF (jcod > 0) THEN
            tcojy=ajjy(nat(i),j)
            tcojz=ajjz(nat(i),j)
            IF (jcod == 2) THEN
              tcojxy=ajxy(nat(i),j)
              tcojxz=ajxz(nat(i),j)
              tcojyx=ajyx(nat(i),j)
              tcojyz=ajyz(nat(i),j)
              tcojzx=ajzx(nat(i),j)
              tcojzy=ajzy(nat(i),j)
            END IF
          END IF

          nh=1

          DO n=j+1,nvs(nat(i))
            IF ( nav(nat(i),n) == nav(nat(i),j) ) THEN
              DO k=1,3
                tco(k)=tco(k)+m(nat(i),n,k)
              END DO
              nh=-nh
              IF(nav(nat(i),n) /= nat(i))nh=1
              tcojx=tcojx+nh*ajjx(nat(i),n)
              IF (jcod > 0) THEN
                tcojy=tcojy+nh*ajjy(nat(i),n)
                tcojz=tcojz+nh*ajjz(nat(i),n)
                IF (jcod == 2) THEN
                  tcojxy=tcojxy+nh*ajxy(nat(i),n)
                  tcojxz=tcojxz+nh*ajxz(nat(i),n)
                  tcojyx=tcojyx+nh*ajyx(nat(i),n)
                  tcojyz=tcojyz+nh*ajyz(nat(i),n)
                  tcojzx=tcojzx+nh*ajzx(nat(i),n)
                  tcojzy=tcojzy+nh*ajzy(nat(i),n)
                END IF
              END IF
              pt(n)=nav(nat(i),n)
            END IF
          END DO
          ii=nav(nat(i),j)
          nhp=1

          DO l=1,nvs(ii)
            IF ( nav(ii,l) == nat(i) ) THEN
              DO k=1,3
                tco(k)=tco(k)+m(nav(nat(i),j),l,k)
              END DO
              nhp=-nhp
              IF (ii /= nat(i)) nhp=1
              tcojx=tcojx-nhp*ajjx(nav(nat(i),j),l)
              IF (jcod > 0) THEN
                tcojy=tcojy-nhp*ajjy(nav(nat(i),j),l)
                tcojz=tcojz-nhp*ajjz(nav(nat(i),j),l)
                IF (jcod > 0) THEN
                  tcojxy=tcojxy-nhp*ajyx(nav(nat(i),j),l)
                  tcojxz=tcojxz-nhp*ajzx(nav(nat(i),j),l)
                  tcojyx=tcojyx-nhp*ajxy(nav(nat(i),j),l)
                  tcojyz=tcojyz-nhp*ajzy(nav(nat(i),j),l)
                  tcojzx=tcojzx-nhp*ajxz(nav(nat(i),j),l)
                  tcojzy=tcojzy-nhp*ajyz(nav(nat(i),j),l)
                END IF
              END IF
            END IF
          END DO

          DO k=1,3
            IF (jcod == 0) THEN
              IF ( (tco(k) /= 0).OR.(tcojx /= 0) ) THEN
                WRITE(lpt,131)
                WRITE(itto,131)
                WRITE(lpt,109)nat(i),nav(nat(i),j),k,  &
                    tco(k),nat(i),nav(nat(i),j), tcojx,nat(i),nav(nat(i),j)
                WRITE(itto,109)nat(i),nav(nat(i),j),k,  &
                    tco(k),nat(i),nav(nat(i),j), tcojx,nat(i),nav(nat(i),j)
                STOP
              END IF

            ELSE IF (jcod == 1) THEN
              IF ( (tco(k) /= 0).OR.(tcojx /= 0).OR.  &
                    (tcojy /= 0).OR.(tcojz /= 0) ) THEN
                WRITE(lpt,131)
                WRITE(itto,131)
                WRITE(lpt,709)nat(i),nav(nat(i),j),k,  &
                    tco(k),nat(i),nav(nat(i),j), tcojx,nat(i),nav(nat(i),j),  &
                    tcojy,nat(i),nav(nat(i),j), tcojz,nat(i),nav(nat(i),j)
                WRITE(itto,709)nat(i),nav(nat(i),j),k,  &
                    tco(k),nat(i),nav(nat(i),j), tcojx,nat(i),nav(nat(i),j),  &
                    tcojy,nat(i),nav(nat(i),j), tcojz,nat(i),nav(nat(i),j)
                STOP
              END IF

            ELSE IF (jcod == 2) THEN
              IF ( (tco(k) /= 0).OR.(tcojx /= 0).OR.  &
                    (tcojy /= 0).OR.(tcojz /= 0).OR.  &
                    (tcojxy /= 0).OR.(tcojxz /= 0).OR.  &
                    (tcojyx /= 0).OR.(tcojyz /= 0).OR.  &
                    (tcojzx /= 0).OR.(tcojzy /= 0) ) THEN
                WRITE(lpt,131)
                WRITE(itto,131)
                WRITE(lpt,809)nat(i),nav(nat(i),j),k,  &
                    tco(k),nat(i),nav(nat(i),j),  &
                    tcojx,tcojxy,tcojxz,nat(i),nav(nat(i),j),  &
                    tcojyx,tcojy,tcojyz,nat(i),nav(nat(i),j),  &
                    tcojzx,tcojzy,tcojz,nat(i),nav(nat(i),j)
                WRITE(itto,809)nat(i),nav(nat(i),j),k,  &
                    tco(k),nat(i),nav(nat(i),j),  &
                    tcojx,tcojxy,tcojxz,nat(i),nav(nat(i),j),  &
                    tcojyx,tcojy,tcojyz,nat(i),nav(nat(i),j),  &
                    tcojzx,tcojzy,tcojz,nat(i),nav(nat(i),j)
                STOP
              END IF
            END IF
          END DO

        END DO  do_j

      END DO

      WRITE(lpt,132)
      WRITE(itto,133)
      RETURN

 4446 WRITE(itto,"(a/a)") ' => Please DO NOT connect sites to themselves with asymmetric coupling',&
                          '    You should double your cell'

      131 FORMAT(/'  ----------------------------------------------------------------------------',/,  &
                  '                         COMPATIBILITY TEST : BAD !!',/,&
                  '                    ****** VERIFY YOUR DATA FILE ! ******',/,  &
                  '  ----------------------------------------------------------------------------')
      132 FORMAT(/'  ----------------------------------------------------------------------------',/,  &
                  '                     COMPATIBILITY TEST : DATA FILE LOOKS OK',/,  &
                  '  ----------------------------------------------------------------------------')
      133 FORMAT( '  ----------------------------------------------------------------------------',/,  &
                  '                     COMPATIBILITY TEST : DATA FILE LOOKS OK',/,  &
                  '  ----------------------------------------------------------------------------')
      109 FORMAT('   Test sym. on  M (',3I2,' ) = ',f5.2,/,'     Test on J (',2I2,' ) = ',f8.3,/,  &
                 ' Dubious connectivity or exchange constant between site ',i2,' and site ',i2,/)
      709 FORMAT('  Test sym. on  M (',3I2,' ) = ',f5.2,/,  &
          '    Test on Jx (',2I2,' ) = ',f8.3,/,  &
          '    Test on Jy (',2I2,' ) = ',f8.3,/,  &
          '    Test on Jz (',2I2,' ) = ',f8.3,/,  &
          ' Dubious connectivity or exchange constant between site ',i2,' and site ',i2,/)
      809 FORMAT('  Test sym. on  M (',3I2,' ) = ',f5.2,/,  &
          '    Test on Jxx,Jxy,Jxz (',2I2,' ) = ',3(f8.3,1X),/,  &
          '    Test on Jyx,Jyy,Jyz (',2I2,' ) = ',3(f8.3,1X),/,  &
          '    Test on Jzx,Jzy,Jzz (',2I2,' ) = ',3(f8.3,1X),/,  &
          ' Dubious connectivity or exchange constant between site ',i2,' and site ',i2,/)


      END SUBROUTINE comp_test


      SUBROUTINE mod_int()

!***********************************************************************

!  This subroutine allows to modify one or more coupling constant(s)
!  within the sample box, for instance to test the perturbation of
!  the magnetic configuration due to the introduction of impurities
!  in the sample. The subroutine is activated only if JCOD=0
!  (isotropic exchange).

!  This subroutine calls no subroutine.

!***********************************************************************

      Use mcm_inc
      implicit none
      integer :: i,j,k,ja,jb,jba, maj1,maj2,maj3
      im=1
 7706 WRITE(itto,7701)
      READ(itti,*)nm(im),ium(im),ivm(im),iwm(im)
      WRITE(iba,*)nm(im),ium(im),ivm(im),iwm(im)

      IF (nm(im) /= 0) THEN
        WRITE(itto,7702)nm(im),nvs(nm(im))
        DO j=1,nvs(nm(im))
          WRITE(itto,7703)j,nav(nm(im),j),(m(nm(im),j,k),k=1,3), ajjx(nm(im),j)
        END DO
        7708   WRITE(itto,7705)
        READ(itti,*)numb(im),ajmod(im)
        WRITE(iba,*)numb(im),ajmod(im)
        IF(numb(im) == 0)GO TO 7706
        ja=nav(nm(im),numb(im))

        DO jb=1,nvs(ja)
          maj1=m(ja,jb,1)+m(nm(im),numb(im),1)
          maj2=m(ja,jb,2)+m(nm(im),numb(im),2)
          maj3=m(ja,jb,3)+m(nm(im),numb(im),3)
          IF ( (nav(ja,jb) == nm(im)).AND.(maj1 == 0).AND.  &
                (maj2 == 0).AND.(maj3 == 0) ) THEN
            jba=jb
          END IF
        END DO

        numb(im+1)=jba
        ajmod(im+1)=ajmod(im)
        nm(im+1)=ja
        ium(im+1)=ium(im)+m(nm(im),numb(im),1)
        IF(ium(im+1) <= 0) ium(im+1)=ium(im+1)+lu
        IF(ium(im+1) > lu) ium(im+1)=ium(im+1)-lu
        ivm(im+1)=ivm(im)+m(nm(im),numb(im),2)
        IF(ivm(im+1) <= 0) ivm(im+1)=ivm(im+1)+lv
        IF(ivm(im+1) > lv) ivm(im+1)=ivm(im+1)-lv
        iwm(im+1)=iwm(im)+m(nm(im),numb(im),3)
        IF(iwm(im+1) <= 0) iwm(im+1)=iwm(im+1)+lw
        IF(iwm(im+1) > lw) iwm(im+1)=iwm(im+1)-lw
        im=im+2
        nm(im)=nm(im-2)
        ium(im)=ium(im-2)
        ivm(im)=ivm(im-2)
        iwm(im)=iwm(im-2)
        GO TO 7708

      END IF

      im=im-1
      WRITE(itto,7711)
      WRITE(lpt,7711)
      DO i=1,im
        WRITE(itto,7712)nm(i),ium(i),ivm(i),iwm(i),numb(i),ajmod(i)
        WRITE(lpt,7712)nm(i),ium(i),ivm(i),iwm(i),numb(i),ajmod(i)
      END DO
      WRITE(itto,"(/)")
      WRITE(lpt,"(/)")

      7701 FORMAT(1X,' Identify site (N) and cell (U,V,W) [0 0 0 0' ,' to stop] :')
      7702 FORMAT(1X,' Site number',i2,' has ',i2,' neighbours :')
      7703 FORMAT(1X,' Neigh. numb.',i2,' = site',i2,' in relative cell'  &
          , ' (',i2,',',i2,',',i2,') with J =',f8.3,' K')
      7705 FORMAT(1X,' Enter modification (Neigh. numb. and new J) [0 0'  &
          ,' to stop] :')
      7711 FORMAT(/,1X,' The following modifications have been applied'  &
          , ' for the current simulation :',/,' (Original file is preserved)')
      7712 FORMAT(1X,' From site',i2,' in cell (',i2,',',i2,',',i2,')'  &
          ,' to its neigh. numb.',i2,' : new J =',f8.3,' K')

      RETURN

      END SUBROUTINE mod_int


      SUBROUTINE gen_spi()

!***********************************************************************

!  This subroutine generates or reads the input spin configuration.
!  In the former case, spins are randomly generated.
!  In the latter case, spin values are read from the file MYFILE.INI.
!  The input format of this file is identical to the format of the
!  file MYFILE.SPI which is generated at the end of every simulation and
!  contains the last generated spin configuration of the simulation. Thus
!  a simulation can be split up in several runs.

!  This subroutine calls the following subroutines:
!     RANnD(RADIUS,X1,Y1...)

!***********************************************************************

      Use mcm_inc
      implicit none
      integer :: u,v,w,ju,jv,jw,n,npassi,npf,npi,jn

      IF(iran == 1) THEN

!*****     FOR ISING SPINS     *****

        DO w=1,lw
          DO v=1,lv
            DO u=1,lu

              IF ((ians == 'I').OR.(ians == 'i')) THEN
                READ(i02,*)
                READ(i02,*)ju,jv,jw
                READ(i02,*)
              END IF

              DO n=1,na
                IF ((ians == 'R') .OR. (ians == 'r')) THEN
                  CALL ran1d(amps(n),x,y,z,n)
                  sx(n,u,v,w)=x
                  sy(n,u,v,w)=y
                  sz(n,u,v,w)=z
                ELSE
                  IF (iani == -1) THEN
                    READ(i02,*)jn,sx(jn,ju,jv,jw)
                  ELSE
                    READ(i02,*)jn,sx(jn,ju,jv,jw),sy(jn,ju,jv,jw), sz(jn,ju,jv,jw)
                  END IF
                END IF
              END DO

            END DO
          END DO
        END DO


      ELSE IF(iran == 2) THEN

!*****     FOR XY SPINS     *****

        DO w=1,lw
          DO v=1,lv
            DO u=1,lu

              IF ((ians == 'I').OR.(ians == 'i')) THEN
                READ(i02,*)
                READ(i02,*)ju,jv,jw
                READ(i02,*)
              END IF

              DO n=1,na
                IF ((ians == 'R').OR.(ians == 'r')) THEN
                  CALL ran2d(amps(n),x,y)
                  sx(n,u,v,w)=x
                  sy(n,u,v,w)=y
                ELSE
                  READ(i02,*)jn,sx(jn,ju,jv,jw),sy(jn,ju,jv,jw)
                END IF
              END DO

            END DO
          END DO
        END DO


      ELSE IF(iran == 3) THEN

!*****     FOR HEISENBERG SPINS     *****

        DO w=1,lw
          DO v=1,lv
            DO u=1,lu

              IF ((ians == 'I').OR.(ians == 'i')) THEN
                READ(i02,*)
                READ(i02,*)ju,jv,jw
                READ(i02,*)
                WRITE(itto,*)ju,jv,jw
              END IF

              DO n=1,na
                IF ((ians == 'R').OR.(ians == 'r')) THEN
                  CALL ran3d(amps(n),x,y,z)
                  sx(n,u,v,w)=x
                  sy(n,u,v,w)=y
                  sz(n,u,v,w)=z
                ELSE
                  READ(i02,*)jn,sx(jn,ju,jv,jw),sy(jn,ju,jv,jw), sz(jn,ju,jv,jw)
                END IF
              END DO

            END DO
          END DO
        END DO


      ELSE IF (iran == 4) THEN

!*****     FOR q-STATE POTTS PLANAR SPINS     *****

        DO w=1,lw
          DO v=1,lv
            DO u=1,lu

              IF ((ians == 'I').OR.(ians == 'i')) THEN
                READ(i02,*)
                READ(i02,*)ju,jv,jw
                READ(i02,*)
              END IF

              DO n=1,na
                IF ((ians == 'R').OR.(ians == 'r')) THEN
                  CALL ranqd(amps(n),x,y)
                  sx(n,u,v,w)=x
                  sy(n,u,v,w)=y
                ELSE
                  READ(i02,*)jn,sx(jn,ju,jv,jw),sy(jn,ju,jv,jw)
                END IF
              END DO

            END DO
          END DO
        END DO


      END IF


!*****     END OF GENERATION AND READING     *****


!*****     OUTPUT     *****

      IF ((ians == 'R') .OR. (ians == 'r')) THEN
        WRITE(lpt,"(///a//)")  '  Random Initial Spin Configuration '
      ELSE
        WRITE(lpt,"(///a,a//)")'  Initial Spin Configuration Read from File: ',trim(nfil)//'.ini'
      END IF

      DO w=1,lw
        DO v=1,lv
          DO u=1,lu

            npi=1
            npf=na
            WRITE(lpt ,"(a,i2,a,i2,a,i2,a)")'  CELL (',u,',',v,',',w,')'
            WRITE(itto,"(a,i2,a,i2,a,i2,a)")'  CELL (',u,',',v,',',w,')'
            WRITE(lpt,"(a)")
            WRITE(itto,"(a)")

            npassi=0
 4412       IF((na-nal*npassi) > nal) npf=npi+nal-1
            npassi=npassi+1
            WRITE(itto,491)
            WRITE(lpt,491)
            DO n=npi,npf
              WRITE(itto,496)n
              WRITE(lpt,496)n
            END DO
            WRITE(lpt,497)(sx(n,u,v,w),n=npi,npf)
            WRITE(itto,497)(sx(n,u,v,w),n=npi,npf)
            IF (iwmx == 1) GO TO 999
            WRITE(lpt,493)(sy(n,u,v,w),n=npi,npf)
            WRITE(itto,493)(sy(n,u,v,w),n=npi,npf)
            IF (iwmxy == 1) GO TO 999
            WRITE(lpt,494)(sz(n,u,v,w),n=npi,npf)
            WRITE(itto,494)(sz(n,u,v,w),n=npi,npf)
  999       WRITE(lpt,595)
            WRITE(itto,595)
            IF (npf >= npi+nal-1) THEN
              IF ((npi /= 1).OR.(npf /= na)) THEN
                npi=npf+1
                npf=na
                IF (npi <= npf) GO TO 4412
              END IF
            END IF

          END DO
        END DO
      END DO


      496 FORMAT('Site',i2,'   ',$)
      497 FORMAT(/1X,' Sx ',50(1X,f6.3,2X))
      493 FORMAT(1X,' Sy ',50(1X,f6.3,2X))
      494 FORMAT(1X,' Sz ',50(1X,f6.3,2X))
      491 FORMAT('+     ',$)
      595 FORMAT(/)

      RETURN

      END SUBROUTINE gen_spi


      SUBROUTINE mc_cycle()
!***********************************************************************
!  This subroutine executes the Monte-Carlo loops according to the
!  following scheme:
!
!               |
!               V
!               |
!               |--------->--------.
!               |                  |
!               |---->---.         |
!               |        |         |
!               |   -----------    |
!               |  | spin loop |   |
!               |  | on sample |   |
!               |   -----------    |
!               |        |         |
!               |----<---'    -----------
!               |            | IT loops  |
!               |            | of MCS/S  |
!               |             -----------
!               |                  |
!               |---------<--------'
!               |
!               V
!               |
!  This subroutine calls the following subroutines and function:
!  ENn_CALC(N,U,V,W)
!  RANF() : this portable real function [Program library of CERN
!           Computer Centre (CERN - Geneva)] is used to generate random
!           numbers between 0 and 1. This function may be slower than
!           dedicated random number generators. Therefore, in order to
!           save CPU time, users should replace RANF() by a call to
!           the equivalent function on the computer they use.
!           For instance, DEC-computer users should call the random
!           generator RAN(I1,I2). For sake of convenience, the calls to
!           this function are included in the current version of the
!           program as comment cards. After implementing and testing
!           the program (see manual), the use of RAN(I1,I2) or of a
!           faster computer random generator, is recommended.
!***********************************************************************
      Use mcm_inc
      implicit none
      real(kind=sp) :: omx,omy,omz,asup,ex,ainf,ax,ay,az
      integer :: i,ini,itim,n,k,u,v,w
      omx=0.0
      omy=0.0
      omz=0.0
      DO w=1,lw
        DO v=1,lv
          DO u=1,lu
            DO k=1,na
              somx(k,u,v,w)=0.0
              somy(k,u,v,w)=0.0
              somz(k,u,v,w)=0.0
              somx2(k,u,v,w)=0.0
              somy2(k,u,v,w)=0.0
              somz2(k,u,v,w)=0.0
            END DO
          END DO
        END DO
      END DO
      ini=0
      npaj=0

      IF (iran == 1) THEN          ! ISING SPINS

!***** MONTE CARLO LOOPS FOR ISING SPINS *****

        DO  itim=1,it
          DO w=1,lw
            DO v=1,lv
              DO u=1,lu
                DO n=1,na
                  CALL en1_calc(n,u,v,w)
                  asup=-delta/t-88.029
                  IF(asup <= 0.0) THEN
                    ainf=-delta/t+89.415
                    IF(ainf < 0) CYCLE
                    ex=EXP(-delta/t)
                    IF(ex-ranf() < 0.0) CYCLE
                  END IF
                  en=en+delta
                  sx(n,u,v,w)=x
                  sy(n,u,v,w)=y
                  sz(n,u,v,w)=z
                  IF (itim > nom) npaj=npaj+1
                END DO
              END DO
            END DO
          END DO

! * * * MEAN VALUES AND VARIANCE OF MAGNETIC MOMENTS AND ENERGY * * *
! FIRST PART OF CALCULATION: SOMMATION OVER THE SNOM=IT-NOM LAST CYCLES

          IF(itim-nom > 0) THEN
            ini=ini+1
          ELSE
            CYCLE
          END IF
          DO w=1,lw
            DO v=1,lv
              DO u=1,lu
                DO k=1,na
                  omx=omx+sx(k,u,v,w)
                  omy=omy+sy(k,u,v,w)
                  omz=omz+sz(k,u,v,w)
                  somx(k,u,v,w)=somx(k,u,v,w)+sx(k,u,v,w)
                  somy(k,u,v,w)=somy(k,u,v,w)+sy(k,u,v,w)
                  somz(k,u,v,w)=somz(k,u,v,w)+sz(k,u,v,w)
                  somx2(k,u,v,w)=somx2(k,u,v,w)+sx(k,u,v,w)**2
                  somy2(k,u,v,w)=somy2(k,u,v,w)+sy(k,u,v,w)**2
                  somz2(k,u,v,w)=somz2(k,u,v,w)+sz(k,u,v,w)**2
                END DO
              END DO
            END DO
          END DO
          enm=enm+en
          en2m=en2m+en*en
          mag2=mag2+omx**2+omy**2+omz**2
          omx=0.0
          omy=0.0
          omz=0.0
        END DO

!***** END OF MONTE CARLO LOOPS FOR ISING SPINS *****

      ELSE IF (iran == 2) THEN     ! XY SPINS

!***** MONTE CARLO LOOPS FOR XY SPINS *****

      DO  itim=1,it
        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              DO n=1,na
                CALL en2_calc(n,u,v,w)
                asup=-delta/t-88.029
                IF(asup <= 0.0) THEN
                  ainf=-delta/t+89.415
                  IF(ainf < 0) CYCLE
                  ex=EXP(-delta/t)
                  IF(ex-ranf() < 0.0) CYCLE
                END IF
                en=en+delta
                sx(n,u,v,w)=x
                sy(n,u,v,w)=y
                IF (itim > nom) npaj=npaj+1
              END DO
            END DO
          END DO
        END DO

! * * * MEAN VALUES AND VARIANCE OF MAGNETIC MOMENTS AND ENERGY * * *
! FIRST PART OF CALCULATION: SOMMATION OVER THE SNOM=IT-NOM LAST CYCLES

        IF(itim-nom > 0) THEN
          ini=ini+1
        ELSE
          CYCLE
        END IF
        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              DO k=1,na
                omx=omx+sx(k,u,v,w)
                omy=omy+sy(k,u,v,w)
                somx(k,u,v,w)=somx(k,u,v,w)+sx(k,u,v,w)
                somy(k,u,v,w)=somy(k,u,v,w)+sy(k,u,v,w)
                somx2(k,u,v,w)=somx2(k,u,v,w)+sx(k,u,v,w)**2
                somy2(k,u,v,w)=somy2(k,u,v,w)+sy(k,u,v,w)**2
              END DO
            END DO
          END DO
        END DO
        enm=enm+en
        en2m=en2m+en*en
        mag2=mag2+omx**2+omy**2
        omx=0.0
        omy=0.0
        omz=0.0
      END DO

!***** END OF MONTE CARLO LOOPS FOR XY SPINS *****

      ELSE IF (iran == 3) THEN     ! HEISENBERG SPINS

!***** MONTE CARLO LOOPS FOR HEISENBERG SPINS *****

      DO  itim=1,it
        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              DO n=1,na
                CALL en3_calc(n,u,v,w)
                asup=-delta/t-88.029
                IF(asup <= 0) THEN
                  ainf=-delta/t+89.415
                  IF(ainf < 0)  CYCLE
                  ex=EXP(-delta/t)
                  IF(ex-ranf() < 0.0) CYCLE
                ENDIF
                en=en+delta
                sx(n,u,v,w)=x
                sy(n,u,v,w)=y
                sz(n,u,v,w)=z
                IF (itim > nom) npaj=npaj+1
              END DO
            END DO
          END DO
        END DO

! * * * MEAN VALUES AND VARIANCE OF MAGNETIC MOMENTS AND ENERGY * * *
! FIRST PART OF CALCULATION: SOMMATION OVER THE SNOM=IT-NOM LAST CYCLES

        IF(itim-nom > 0) THEN
          ini=ini+1
        ELSE
          CYCLE
        END IF
        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              DO k=1,na
                omx=omx+sx(k,u,v,w)
                omy=omy+sy(k,u,v,w)
                omz=omz+sz(k,u,v,w)
                somx(k,u,v,w)=somx(k,u,v,w)+sx(k,u,v,w)
                somy(k,u,v,w)=somy(k,u,v,w)+sy(k,u,v,w)
                somz(k,u,v,w)=somz(k,u,v,w)+sz(k,u,v,w)
                somx2(k,u,v,w)=somx2(k,u,v,w)+sx(k,u,v,w)**2
                somy2(k,u,v,w)=somy2(k,u,v,w)+sy(k,u,v,w)**2
                somz2(k,u,v,w)=somz2(k,u,v,w)+sz(k,u,v,w)**2
              END DO
            END DO
          END DO
        END DO
        enm=enm+en
        en2m=en2m+en*en
        mag2=mag2+omx**2+omy**2+omz**2
        omx=0.0
        omy=0.0
        omz=0.0
      END DO


!***** END OF MONTE CARLO LOOPS FOR HEISENBERG SPINS *****

      ELSE IF (iran == 4) THEN     ! q-STATE POTTS PLANAR SPINS

!***** MONTE CARLO LOOPS FOR q-STATE PLANAR POTTS SPINS *****

      DO  itim=1,it
        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              DO n=1,na
                CALL enq_calc(n,u,v,w)
                asup=-delta/t-88.029
                IF(asup <= 0.0) THEN
                  ainf=-delta/t+89.415
                  IF(ainf < 0) CYCLE
                  ex=EXP(-delta/t)
                  IF(ex-ranf() < 0.0) CYCLE
                END IF
                en=en+delta
                sx(n,u,v,w)=x
                sy(n,u,v,w)=y
                IF (itim > nom) npaj=npaj+1
              END DO
            END DO
          END DO
        END DO

! * * * MEAN VALUES AND VARIANCE OF MAGNETIC MOMENTS AND ENERGY * * *
! FIRST PART OF CALCULATION: SOMMATION OVER THE SNOM=IT-NOM LAST CYCLES

        IF(itim-nom > 0) THEN
          ini=ini+1
        ELSE
          CYCLE
        END IF
        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              DO k=1,na
                omx=omx+sx(k,u,v,w)
                omy=omy+sy(k,u,v,w)
                somx(k,u,v,w)=somx(k,u,v,w)+sx(k,u,v,w)
                somy(k,u,v,w)=somy(k,u,v,w)+sy(k,u,v,w)
                somx2(k,u,v,w)=somx2(k,u,v,w)+sx(k,u,v,w)**2
                somy2(k,u,v,w)=somy2(k,u,v,w)+sy(k,u,v,w)**2
              END DO
            END DO
          END DO
        END DO
        enm=enm+en
        en2m=en2m+en*en
        mag2=mag2+omx**2+omy**2
        omx=0.0
        omy=0.0
        omz=0.0
      END DO


!***** END OF MONTE CARLO LOOPS FOR q-STATE POTTS SPINS *****


      END IF

      DO w=1,lw
        DO v=1,lv
          DO u=1,lu
            DO k=1,na
              somx(k,u,v,w)=somx(k,u,v,w)/snom
              comx(k,u,v,w)=somx(k,u,v,w)
              somy(k,u,v,w)=somy(k,u,v,w)/snom
              somz(k,u,v,w)=somz(k,u,v,w)/snom
              somx2(k,u,v,w)=somx2(k,u,v,w)/snom
              somy2(k,u,v,w)=somy2(k,u,v,w)/snom
              somz2(k,u,v,w)=somz2(k,u,v,w)/snom
              ax=somx2(k,u,v,w)-somx(k,u,v,w)**2
              IF(ax < 0) ax=0
              ay=somy2(k,u,v,w)-somy(k,u,v,w)**2
              IF(ay < 0) ay=0
              az=somz2(k,u,v,w)-somz(k,u,v,w)**2
              IF(az < 0) az=0
              somx2(k,u,v,w)=SQRT(ax)
              somy2(k,u,v,w)=SQRT(ay)
              somz2(k,u,v,w)=SQRT(az)
            END DO
          END DO
        END DO
      END DO

      RETURN

      END SUBROUTINE mc_cycle


      SUBROUTINE en_tot(L,en1)
!***********************************************************************
!  This subroutine computes the total energy of the system.
!  This subroutine calls the following subroutine:
!     EN3_CALC(N,U,V,W)
!***********************************************************************
      Use mcm_inc
      implicit none
      integer, intent(in) :: L
      real(kind=dp),intent(out) :: en1
      real(kind=dp) :: en0
      integer :: i,n,k,u,v,w,ier
      real(kind=sp) :: aux
      character(len=80) :: time_line
      nmult=2
      DO i=1,na
        amps2(i)=0
      END DO
      !For L=3 the copy the average magnetic moments into spin variables
      !has been done before calling the subroutine
      x=0.0; y=0.0; z=0.0 !Initialize the spin components
      en0=0               !Needed for calculating the total energy with the
      nfip=0              !current spin configuration [Sx,Sy,Sz](n,u,v,w)
      DO w=1,lw
        DO v=1,lv
          DO u=1,lu
            DO n=1,na
              CALL en3_calc(n,u,v,w)
              en0=en0-delta
            END DO
          END DO
        END DO
      END DO
      en=en0/2
      en1=en/nfa
      IF (l == 1) THEN
        WRITE(lpt ,"(a,f12.5)")'       *  Energy of the starting spin configuration = ',en1
        WRITE(itto,"(a,f12.5)")'       *  Energy of the starting spin configuration = ',en1
      ELSE IF (l == 2) then
        WRITE(lpt ,"(a,f12.5)")'       *  Energy of the last spin configuration = ',en1
        WRITE(itto,"(a,f12.5)")'       *  Energy of the last spin configuration = ',en1
        do i=1,10
         call write_fst(ier)
         if(ier == 0) exit
        end do
        iwrt=iwrt-1
        call cpu_time(cpu_par)
        cpu= cpu_par-cpu
        remaining=iwrt*cpu/60.0/60.0  !remaining time in decimal hours
        time_line(1:80)=" "
        write(unit=time_line(1:10),fmt="(i3,a)") int(remaining)," hours,"
        remaining=(remaining-int(remaining))*60.0 !remaining time in decimal minutes
        write(unit=time_line(12:21),fmt="(i2,a)") int(remaining)," minutes"
        aux=(remaining - int(remaining))*60.0 !Number of seconds remaining
        write(unit=time_line(23:50),fmt="(f8.3,a)") aux," seconds"
        write(unit=*,fmt="(a,f8.2,a,a,/)")" => Partial time:",cpu*60.0, &
               " seconds  ->  Approx. remaining time: ",trim(time_line)
        cpu=cpu_par
      ELSE IF (l == 3) then
        WRITE(lpt ,"(a,f12.5)")'       *  Energy of the average spin configuration = ',en1
        WRITE(itto,"(a,f12.5)")' => Energy of the average spin configuration = ',en1
      END IF
      nfip=1
      DO i=1,na
        amps2(i)=amps(i)
      END DO
      nmult=1
      RETURN

      END SUBROUTINE en_tot


      SUBROUTINE const_funct()

!***********************************************************************
!  This subroutine calculates the constraint function Fc of the system
!  (for isotropic exchange only), according to the formula:
!  (with the convention that a vector V' represents the transposed
!  form of a vector V)
!
!       Fc = -E/Eb  with  E  = - JijSi'Sj             (energy)
!                    and  Eb = -|Jij||Si||Sj|      (basic energy)
!
!  Further details about the constraint function are given in the
!  following paper:
!  The constaint functions: an attempt to evaluate the constraint
!  rate inside compounds that undergo ordered magnetic frustration.
!  P. Lacorre, J. Phys. C: Solid State Phys. 20 (1987) L775-L781
!
!  This subroutine calls the following subroutine:
!     BOUND_COND
!
!***********************************************************************


      Use mcm_inc
      implicit none
      real(kind=sp) :: fcnu,fcde
      integer :: k,n,u,v,w
      fcnu = 0
      fcde = 0
      DO w=1,lw
        DO v=1,lv
          DO u=1,lu
            DO n=1,na
              fjsx(n)=0
              fjsy(n)=0
              fjsz(n)=0
              fjsm(n)=0
              DO k=1,nvs(n)

                CALL bound_cond(n,k,u,v,w)

                fjsx(n)=fjsx(n)+ajxx(n,k)*comx(nav(n,k),nu(k),nv(k),nw(k))
                fjsy(n)=fjsy(n)+ajyy(n,k)*somy(nav(n,k),nu(k),nv(k),nw(k))
                fjsz(n)=fjsz(n)+ajzz(n,k)*somz(nav(n,k),nu(k),nv(k),nw(k))
                fjsm(n)=fjsm(n)+ABS(ajxx(n,k)) *tmom(nav(n,k),nu(k),nv(k),nw(k))
              END DO

! * * *   COMPUTATION OF THE NUMERATOR (FCNU=2*(MEAN ENERGY))   * * *
! * * *   AND DENOMINATOR (FCDE=-2*(BASIC MEAN ENERGY)) OF FC   * * *

             fcnu=fcnu-comx(n,u,v,w)*fjsx(n)-somy(n,u,v,w)*fjsy(n)-somz(n,u,v,w)*fjsz(n)
             fcde=fcde+tmom(n,u,v,w)*fjsm(n)
            END DO
          END DO
        END DO
      END DO
      IF (fcde < 0.000001) fcde=0.000001
      fconst=fcnu/fcde
      rc=50*(1+fconst)

      RETURN
      END SUBROUTINE const_funct


      SUBROUTINE pro_rota()
!***********************************************************************
!  This subroutine computes, at the end of simulation, the projection
!  of magnetic moments onto the crystal axes, and rotates the whole
!  system in order to align a given magnetic moment along the desired
!  direction.
!  Calculation for rotation are performed within the orthonormal coord-
!  inate system, then magnetic moments are projected onto the crystal-
!  lographic reference frame.
!  The output file MYFILE.LIS stores :
!  - the components of magnetic moments inside the crystallographic
!    reference frame before rotation
!  - the components of magnetic moments inside the crystallographic
!    reference frame after rotation
!
!  The following parameters are used in the subroutine :
!
!                    TRANSFORMATION MATRICES
!         between orthonormal coordinate system (o.c.s.)
!         and crystallographic reference frame (c.r.f.)
!  ( ref. : V.A.LIOPO_Sov. Phys. Crystallogr. 30 (1985) 687)
!
!  AA : transformation matrix (c.r.f.) ---> (o.c.s.) for magnetic moments
!  AC=INV(AA) : transformation matrix (o.c.s.) ---> (c.r.f.)
!  AD : transformation matrix (c.r.f.) ---> (o.c.s.) for direction
!
!                           VECTORS
!  SORX,SORY,SORZ : components of the magnetic moment taken as reference
!  DORX,DORY,DORZ : new direction of the magnetic moment
!  UORX,UORY,UORZ : rotation axis ( UOR = SOR^DOR )
!  OM = angle(SOR,DOR)  [ SIOM = SIN(OM) , COOM = COS(OM) ]
!
!                       ROTATION FORMULA
!                       ----------------
!  ( ref.: J.BASS_Cours de Mathematiques T.I f.1 MASSON(1977) )
!                       ( p 130 - ex 36 )
!
!  if S = arbitrary vector
!     S'= transformation of S by the rotation
!     U = normalised UOR
!
!         S' = (1-COOM)*(U.S)*U + COOM*S + SIOM*(U^S)
!         -------------------------------------------
!
!
!  This subroutine calls the following subroutine:
!     conv_matrix
!
!***********************************************************************


      Use mcm_inc
      implicit none
      integer       :: ier,k,n,npi,npf,npassa,npassb,u,v,w
      real(kind=sp) :: sorx,sory,sorz,dor,sorm, sor,coom, uorx, &
                       uory,uorz, uor, siom, paj,ammc,ampara,ammca,pscal
      real(kind=dp) :: enma,enm2
      real(kind=dp),dimension(3) :: vecto,vectx

      !READ(i01,*)aax,aay,aaz,aaa,aab,aag  !Removed from here for reading cell event if PRO_ROTA is not invoked


! * * * NOR = INDEX OF MOMENT TAKEN AS REFERENCE (CELL (1,1,1)) * * *
      if(.not. batch) then
         DO
           WRITE(itto,"(a)") ' => Do you want to redefine the orientation of magnetic moments'
           WRITE(itto,"(a)") '    in the crystal cell ?'
           WRITE(itto,"(a)") '              0    : do not'
           WRITE(itto,"(a)") '           N =< NA : ordinal number of the magnetic moment to be'
           WRITE(itto,"(a)") '                     used to calculate the new orientation'
           WRITE(itto,"(a)",advance='no') ' => Option: '
           READ(itti,*,iostat=ier) nor
           if(ier /= 0) cycle
           IF ((nor < 0) .OR. (nor > na)) CYCLE
           EXIT
         END DO
         WRITE(iba,*)nor
      end if
      !--------------------
      IF (nor /= 0) THEN
      !--------------------
        sorx=comx(nor,1,1,1)
        sory=somy(nor,1,1,1)
        sorz=somz(nor,1,1,1)
        if(.not. batch) then
           Do
             Write(itto,"(a,i2,a)",advance='no')' => Give new orientation of magnetic moment of SITE ',nor, &
                              ' in the crystal cell: '
             Read(itti,*,iostat=ier) dorx,dory,dorz
             If(ier /= 0) cycle
             Exit
           End Do
           Write(iba,*)dorx,dory,dorz
       end if
! DORX,DORY,DORZ : NEW DIRECTION OF MOMENT NOR INSIDE THE CRYSTAL. REFERENCE

        vectx=[dorx,dory,dorz]  !Crystallographic direction
        dorix=dorx
        doriy=dory
        doriz=dorz
        vecto=matmul(ad,vectx)  !In Cartesian frame  Vc= AD Vx

        dorx=vecto(1) !dorx*ad(1,1)+dory*ad(1,2)+dorz*ad(1,3)   This is an error
        dory=vecto(2) !dorx*ad(2,1)+dory*ad(2,2)+dorz*ad(2,3)
        dorz=vecto(2) !dorx*ad(3,1)+dory*ad(3,2)+dorz*ad(3,3)

! DORX,DORY,DORZ : NEW DIRECTION OF MOMENT NOR INSIDE THE ORTHO. SYSTEM

! * * * CALCULATION INSIDE THE ORTHONORMAL COORDINATE SYSTEM * * *

        dor=SQRT(dorx*dorx+dory*dory+dorz*dorz)
        sorm=SQRT(sorx*sorx+sory*sory+sorz*sorz)
        sor=dorx*sorx+dory*sory+dorz*sorz
        coom=sor/(dor*sorm)
        uorx=sory*dorz-sorz*dory
        uory=sorz*dorx-sorx*dorz
        uorz=sorx*dory-sory*dorx
        uor=SQRT(uorx*uorx+uory*uory+uorz*uorz)
        siom=uor/(dor*sorm)
        uorx=uorx/uor
        uory=uory/uor
        uorz=uorz/uor

        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              DO k=1,na
                pscal= uorx*comx(k,u,v,w) +uory*somy(k,u,v,w)+uorz*somz(k,u,v,w)
                somrx(k,u,v,w)= (1-coom)*pscal*uorx + coom*comx(k,u,v,w)  &
                                + siom*( uory*somz(k,u,v,w) -uorz*somy(k,u,v,w) )
                somry(k,u,v,w)= (1-coom)*pscal*uory + coom*somy(k,u,v,w)  &
                                + siom*( uorz*comx(k,u,v,w) -uorx*somz(k,u,v,w) )
                somrz(k,u,v,w)= (1-coom)*pscal*uorz + coom*somz(k,u,v,w)  &
                                + siom*( uorx*somy(k,u,v,w) -uory*comx(k,u,v,w) )

!   *   SOMRX           COMPONENTS OF MAGNETIC MOMENTS   *
!   *   SOMRY   ------> INSIDE THE ORTHONORMAL SYSTEM    *
!   *   SOMRZ           AFTER ROTATION                   *


                somrx2(k,u,v,w)= (1-coom)*pscal*uorx + coom*somx2(k,u,v,w)  &
                                 + siom*( uory*somz2(k,u,v,w) -uorz*somy2(k,u,v,w) )
                somry2(k,u,v,w)= (1-coom)*pscal*uory + coom*somy2(k,u,v,w)  &
                                 + siom*( uorz*somx2(k,u,v,w) -uorx*somz2(k,u,v,w) )
                somrz2(k,u,v,w)= (1-coom)*pscal*uorz + coom*somz2(k,u,v,w)  &
                                 + siom*( uorx*somy2(k,u,v,w) -uory*somx2(k,u,v,w) )

!   *   SOMRX2          VARIANCE OF THE COMPONENTS OF       *
!   *   SOMRY2  ------> MAGNETIC MOMENTS INSIDE THE         *
!   *   SOMRZ2          ORTHONORMAL SYSTEM AFTER ROTATION   *


              END DO
            END DO
          END DO
        END DO


! * * * END OF CALCULATION INSIDE THE ORTHONORMAL SYSTEM * * *

! * * * COMPUTATION INSIDE THE CRISTALLOGRAPHIC REFERENCE FRAME * * *

!   *   PROJECTION OF MAGNETIC MOMENTS AFTER ROTATION   *
!                  (SOMRX,SOMRY,SOMRZ)

        DO w=1,lw
          DO v=1,lv
            DO u=1,lu
              DO k=1,na
                 vecto=[somrx(k,u,v,w),somry(k,u,v,w),somrz(k,u,v,w)]  !Mx= AC Mc
                 vectx=matmul(ac,vecto)  !Below there was an error
                 somrx(k,u,v,w)= vectx(1)!ac(1,1)* somrx(k,u,v,w) +ac(1,2)* somry(k,u,v,w)+ac(1,3)* somrz(k,u,v,w)
                 somry(k,u,v,w)= vectx(2)!ac(2,1)* somrx(k,u,v,w) +ac(2,2)* somry(k,u,v,w)+ac(2,3)* somrz(k,u,v,w)
                 somrz(k,u,v,w)= vectx(3)!ac(3,1)* somrx(k,u,v,w) +ac(3,2)* somry(k,u,v,w)+ac(3,3)* somrz(k,u,v,w)
                 vecto=[somrx2(k,u,v,w),somry2(k,u,v,w),somrz2(k,u,v,w)]
                 vectx=matmul(ac,vecto)
                                         !Below there was an error
                somrx2(k,u,v,w)= vectx(1)!ac(1,1)*somrx2(k,u,v,w) +ac(1,2)*somry2(k,u,v,w)+ac(1,3)*somrz2(k,u,v,w)
                somry2(k,u,v,w)= vectx(2)!ac(2,1)*somrx2(k,u,v,w) +ac(2,2)*somry2(k,u,v,w)+ac(2,3)*somrz2(k,u,v,w)
                somrz2(k,u,v,w)= vectx(3)!ac(3,1)*somrx2(k,u,v,w) +ac(3,2)*somry2(k,u,v,w)+ac(3,3)*somrz2(k,u,v,w)
              END DO
            END DO
          END DO
        END DO
      !--------------------
      END IF  !nor /= 0
      !--------------------

!   *   PROJECTION OF MAGNETIC MOMENTS BEFORE ROTATION   *
!                      (COMX,SOMY,SOMZ)

      DO w=1,lw
        DO v=1,lv
          DO u=1,lu
            DO k=1,na
               vecto=[comx(k,u,v,w),somy(k,u,v,w),somz(k,u,v,w)] !Cartesian
               vectx=matmul(ac,vecto) ! Xtallographic: Below there was an error
               comx(k,u,v,w)= vectx(1)!ac(1,1)* comx(k,u,v,w) +ac(1,2)* somy(k,u,v,w)+ac(1,3)* somz(k,u,v,w)
               somy(k,u,v,w)= vectx(2)!ac(2,1)* comx(k,u,v,w) +ac(2,2)* somy(k,u,v,w)+ac(2,3)* somz(k,u,v,w)
               somz(k,u,v,w)= vectx(3)!ac(3,1)* comx(k,u,v,w) +ac(3,2)* somy(k,u,v,w)+ac(3,3)* somz(k,u,v,w)

               vecto=[somx2(k,u,v,w),somy2(k,u,v,w),somz2(k,u,v,w)]
               vectx=matmul(ac,vecto)
                                      !Below there was an error
              somx2(k,u,v,w)= vectx(1)!ac(1,1)*somx2(k,u,v,w) +ac(1,2)*somy2(k,u,v,w)+ac(1,3)*somz2(k,u,v,w)
              somy2(k,u,v,w)= vectx(2)!ac(2,1)*somx2(k,u,v,w) +ac(2,2)*somy2(k,u,v,w)+ac(2,3)*somz2(k,u,v,w)
              somz2(k,u,v,w)= vectx(3)!ac(3,1)*somx2(k,u,v,w) +ac(3,2)*somy2(k,u,v,w)+ac(3,3)*somz2(k,u,v,w)
            END DO
          END DO
        END DO
      END DO

! ***** END OF CALCULATION FOR PROJECTION AND ROTATION *****

! ***** OUTPUT RESULTS OF PROJECTION AND ROTATION *****

      WRITE(lpt,517)aax,aay,aaz,aaa,aab,aag
      WRITE(itto,"(/)")
      WRITE(lpt,"(/)")
      WRITE(itto,618)
      WRITE(lpt,618)
      DO w=1,lw
        DO v=1,lv
          DO u=1,lu
            DO n=1,na
              isigx(n)=NINT(1000*somx2(n,u,v,w))
              isigy(n)=NINT(1000*somy2(n,u,v,w))
              isigz(n)=NINT(1000*somz2(n,u,v,w))
            END DO
            npi=1
            npf=na
            WRITE(lpt,498)u,v,w
            WRITE(itto,498)u,v,w
            WRITE(lpt,"(a)")
            WRITE(itto,"(a)")

            npassa=0
 3412       IF((na-nal*npassa) > nal) npf=npi+nal-1
            npassa=npassa+1
            WRITE(itto,491)
            WRITE(lpt,491)
            DO n=npi,npf
              WRITE(itto,496)n
              WRITE(lpt,496)n
            END DO
            WRITE(itto,"(a)")
            WRITE(lpt,"(a)")
            WRITE(lpt,497)(comx(n,u,v,w),n=npi,npf)
            WRITE(lpt,487)(isigx(n),n=npi,npf)
            WRITE(lpt,493)(somy(n,u,v,w),n=npi,npf)
            WRITE(lpt,487)(isigy(n),n=npi,npf)
            WRITE(lpt,494)(somz(n,u,v,w),n=npi,npf)
            WRITE(lpt,487)(isigz(n),n=npi,npf)
            WRITE(itto,497)(comx(n,u,v,w),n=npi,npf)
            WRITE(itto,487)(isigx(n),n=npi,npf)
            WRITE(itto,493)(somy(n,u,v,w),n=npi,npf)
            WRITE(itto,487)(isigy(n),n=npi,npf)
            WRITE(itto,494)(somz(n,u,v,w),n=npi,npf)
            WRITE(itto,487)(isigz(n),n=npi,npf)
            WRITE(itto,"(/)")
            WRITE(lpt,"(/)")
            IF (npf >= npi+nal-1) THEN
              IF ( (npi /= 1).OR.(npf /= na) ) THEN
                npi=npf+1
                npf=na
                IF (npi < npf) GO TO 3412
              END IF
            END IF
          END DO
        END DO
      END DO

      IF(nor == 0) return
      WRITE(itto,"(/)")
      WRITE(lpt,"(/)")
      WRITE(itto,638)nor,dorix,doriy,doriz
      WRITE(lpt,638)nor,dorix,doriy,doriz
      WRITE(itto,"(/)")
      WRITE(lpt,"(/)")
      WRITE(itto,628)
      WRITE(lpt,628)
      DO w=1,lw
        DO v=1,lv
          DO u=1,lu
            npi=1
            npf=na
            WRITE(lpt,498)u,v,w
            WRITE(itto,498)u,v,w
            WRITE(lpt,"(a)")
            WRITE(itto,"(a)")

            npassb=0
 2412       IF((na-nal*npassb) > nal) npf=npi+nal-1
            npassb=npassb+1
            WRITE(itto,491)
            WRITE(lpt,491)
            DO n=npi,npf
              WRITE(itto,496)n
              WRITE(lpt,496)n
            END DO
            WRITE(itto,"(a)")
            WRITE(lpt,"(a)")
            WRITE(lpt,497)(somrx(n,u,v,w),n=npi,npf)
            WRITE(lpt,493)(somry(n,u,v,w),n=npi,npf)
            WRITE(lpt,494)(somrz(n,u,v,w),n=npi,npf)
            WRITE(itto,497)(somrx(n,u,v,w),n=npi,npf)
            WRITE(itto,493)(somry(n,u,v,w),n=npi,npf)
            WRITE(itto,494)(somrz(n,u,v,w),n=npi,npf)
            WRITE(itto,"(/)")
            WRITE(lpt,"(/)")
            IF (npf >= npi+nal-1) THEN
              IF ( (npi /= 1).OR.(npf /= na) ) THEN
                npi=npf+1
                npf=na
                IF (npi > npf) GO TO 2412
              END IF
            END IF
          END DO
        END DO
      END DO


      613 FORMAT(1X,3(f8.4,1X),6X,3(f8.4,1X))
      496 FORMAT('Site',i2,'   ',$)
      497 FORMAT(1X,' Mx ',50(1X,f6.3,2X))
      493 FORMAT(1X,' My ',50(1X,f6.3,2X))
      494 FORMAT(1X,' Mz ',50(1X,f6.3,2X))
      487 FORMAT(5X,50(3X,i4,2X))
      491 FORMAT('+     ',$)
      498 FORMAT(1X,' CELL (',i2,',',i2,',',i2,')')
        517 FORMAT(1X,' Crystal CELL parameters :',/,5X,' a = ',f9.4,  &
            /,5X,' b = ',f9.4,/,5X,' c = ',f9.4,/,5X,' alpha = ',f9.4,  &
            /,5X,' beta  = ',f9.4,/,5X,' gamma = ',f9.4)
        618 FORMAT(1X,' COMPONENTS OF THE MAGNETIC MOMENTS IN THE CRYSTAL'  &
            , ' CELL ',/,'   *************************************'  &
            , '***************',/)
        638 FORMAT(1X,' ***** MOMENTS ARE ROTATED TO ALIGN MOMENT OF SITE ',  &
            i2,' ALONG THE DIRECTION *****',/,13X,'***** ',  &
            3(f7.3,1X),'OF THE CRYSTAL CELL  *****')
        628 FORMAT(1X,' COMPONENTS OF THE MAGNETIC MOMENTS IN THE CRYSTAL'  &
            , ' CELL AFTER ROTATION',/,'   ***********************'  &
            , '********************************************',/)

            RETURN

      END SUBROUTINE pro_rota


      SUBROUTINE out_res()
!***********************************************************************
!  This subroutine outputs the averaged characteristics in files
!  MYFILE.LIS and MYFILE.RES (if activated) at the end of every
!  temperature/field loop.
!
!  The calculated averaged characteristics are the following :
!  - mean total energy (ENMA)
!  - mean total magnetization (AMMCA)
!  - mean total magnetization along the applied field (AMPARA)
!  - mean magnetization on three selected sites of cell (1,1,1)
!  - specific heat (CVA)
!  - magnetic susceptibility (SUSA)
!  - inverse susceptibility (SUAM1)
!
!  The specific heat and magnetic susceptibility are computed according
!  to the formula :
!  CVA = (<E**2>-<E>**2)/T**2    (where E represents the energy)
!  SUSA = (<M**2>-<M>**2)/T      (where M represents the magnetization)
!
!  This subroutine calls no subroutine.
!***********************************************************************

      Use mcm_inc
      implicit none

      real(kind=sp) :: paj, pscal
      real(kind=sp) :: xm,ym,zm,susa,suam1,cva, &
                       ammc,ampara,ammca
      real(kind=dp) :: enma,enm2

      paj=100.0*real(npaj)/(snom*real(nta))

      enm=enm/snom
      enm2=enm*enm
      en2m=en2m/snom
      mag2=mag2/snom
      ammc=SQRT(amx*amx+amy*amy+amz*amz)

      enma=enm/nfa
      ampara=hdix*amx+hdiy*amy+hdiz*amz
      ammca=ammc/nfa
      ampara=ampara/nfa
      IF (ammc /= 0.0) THEN
        xm=amx/ammc
        ym=amy/ammc
        zm=amz/ammc
      ELSE
        xm=0.0
        ym=0.0
        zm=0.0
      END IF
      susa=(mag2-ammc**2)/t
      susa=susa/nfa
      IF (susa > 0.0) THEN
        suam1=1.0/susa
        IF(ABS(suam1) > 9999.999) suam1=9999.999
      END IF
      cva=(en2m-enm2)/(t*t)
      cva=cva/nfa
      IF(cva > 9999.999) cva=9999.999
      enm=0
      en2m=0
      mag2=0
      ammc=0
      IF (last == 1) last=0
      IF (jcod /= 0) THEN
        WRITE(lpt,993)paj,xm,ym,zm,hx0,hy0,hz0
        WRITE(itto,993)paj,xm,ym,zm,hx0,hy0,hz0
      ELSE
        WRITE(lpt,693)paj,fconst,rc,xm,ym,zm,hx0,hy0,hz0
        WRITE(itto,693)paj,fconst,rc,xm,ym,zm,hx0,hy0,hz0
      693 FORMAT(7X,'*  Rate of accepted jumps : ',f7.3,' %',/,  &
          7X,'*  Fc = ',f8.5,'   -   Rc = ',f7.3,' %',/,  &
          7X,'*  MTOT   along the direction ',3(1X,f6.3),/,  &
          7X,'*  M//    along the direction ',3(1X,f6.3))
      END IF
      WRITE(itto,263)j1,j2,j3
      WRITE(lpt,263)j1,j2,j3
      WRITE(itto,264)susa,suam1,cva,enma,ammca,  &
          tmom(j1,1,1,1),tmom(j2,1,1,1),tmom(j3,1,1,1), ampara
      WRITE(lpt,264)susa,suam1,cva,enma,ammca,  &
          tmom(j1,1,1,1),tmom(j2,1,1,1),tmom(j3,1,1,1), ampara
      IF ((ansrf == 'Y').OR.(ansrf == 'y')) THEN
        WRITE(i03,265)t,hk,suam1,cva,enma,ammca,  &
            tmom(j1,1,1,1),tmom(j2,1,1,1),tmom(j3,1,1,1), ampara,fconst
      END IF


      993 FORMAT(7X,'*  Rate of accepted jumps : ',f7.3,' %',/,  &
          7X,'*  MTOT   along the direction ',3(1X,f6.3),/,  &
          7X,'*  M//    along the direction ',3(1X,f6.3))
      264 FORMAT(7X,'*',f9.4,1X,f8.3,1X,f8.3,1X,f10.2,1X,f7.3,3(1X,f6.3), 1X,f7.3)
      265 FORMAT(3X,f7.2,1X,f7.2,1X,f8.3,1X,f8.3,1X,f10.2,1X,f7.3,  &
          3(1X,f6.3),1X,f7.3,1X,f6.3)
      263 FORMAT(7X,'*    SUSC',4X,'1/SUS',4X,'SP. HT.',6X,'EN',5X,'MTOT',  &
          1X,3(1X,'M(',i2,')'),4X,'M//')

      RETURN

      END SUBROUTINE out_res


      SUBROUTINE out_conf()
!***********************************************************************
!  This subroutine output the configuration of magnetic moments at
!  current temperature and field.
!  This subroutine calls no subroutine.
!***********************************************************************
      Use mcm_inc
      implicit none
      integer :: npi,npf,k,n,l,u,v,w
      logical :: prt,prta,prte

      amx=0.0
      amy=0.0
      amz=0.0
      prt=((ansnev == 'E').OR.(ansnev == 'e'))
      prta=((ansnev == 'A').OR.(ansnev == 'a'))
      prte= (prt .or. prta) .AND. (last == 1)

      IF ( prte) then

        WRITE(lpt,"(/a/a/)") &
          '                    Mean Configuration of Magnetic Moments',   &
          '                            ---------------------'
      END IF

      DO w=1,lw
        DO v=1,lv
          DO u=1,lu

            DO k=1,na
              tmom(k,u,v,w)=SQRT( somx(k,u,v,w)**2 +somy(k,u,v,w)**2  &
                  +somz(k,u,v,w)**2 )
            END DO
            npi=1
            npf=na

            ! Modification of the original IF block!  (caused by illegal goto 1412!)

            IF ( prta .OR. prt)  THEN

                if( prte) then !Print only at the end
                   WRITE(itto,"(a,i2,a,i2,a,i2,a,/)")'  CELL (', u, ',', v, ',', w, ')'
                   WRITE( lpt,"(a,i2,a,i2,a,i2,a,/)")'  CELL (', u, ',', v, ',', w, ')'

                   WRITE(itto,"(3x,50(a,i2))")("   Site",n,n=npi,npf)
                   WRITE( lpt,"(3x,50(a,i2))")("   Site",n,n=npi,npf)
                end if

              DO k=npi,npf
                isigx(k)=NINT(1000*somx2(k,u,v,w))
                isigy(k)=NINT(1000*somy2(k,u,v,w))
                isigz(k)=NINT(1000*somz2(k,u,v,w))
              END DO

              if(prte) then !Print only at the end
                 WRITE(lpt,"(a,50(1X,f6.3,2X))") "   Mx",(somx(n,u,v,w),n=npi,npf)
                 WRITE(lpt,"(5X,50(3X,i4,2X))")          (isigx(n),n=npi,npf)
                 WRITE(itto,"(a,50(1X,f6.3,2X))")"   Mx",(somx(n,u,v,w),n=npi,npf)
                 WRITE(itto,"(5X,50(3X,i4,2X))")         (isigx(n),n=npi,npf)
                 IF(iwmx /= 1) then
                    WRITE(lpt,"(a,50(1X,f6.3,2X))") "   My",(somy(n,u,v,w),n=npi,npf)
                    WRITE(lpt,"(5X,50(3X,i4,2X))")          (isigy(n),n=npi,npf)
                    WRITE(itto,"(a,50(1X,f6.3,2X))")"   My",(somy(n,u,v,w),n=npi,npf)
                    WRITE(itto,"(5X,50(3X,i4,2X))")         (isigy(n),n=npi,npf)
                    IF(iwmxy /= 1) then
                       WRITE(lpt,"(a,50(1X,f6.3,2X))") "   Mz",(somz(n,u,v,w),n=npi,npf)
                       WRITE(lpt,"(5X,50(3X,i4,2X))")          (isigz(n),n=npi,npf)
                       WRITE(itto,"(a,50(1X,f6.3,2X))")"   Mz",(somz(n,u,v,w),n=npi,npf)
                       WRITE(itto,"(5X,50(3X,i4,2X))")         (isigz(n),n=npi,npf)
                    end if
                 end if
              end if
              do k=npi,npf
                if (tmom(k,u,v,w) > 0.00001) then
                  somx(k,u,v,w)=( ABS(somx(k,u,v,w))*somx2(k,u,v,w)  &
                      +ABS(somy(k,u,v,w))*somy2(k,u,v,w)  &
                      +ABS(somz(k,u,v,w))*somz2(k,u,v,w) ) /tmom(k,u,v,w)
                else
                  somx(k,u,v,w)=amps(k)
                end if
                isigm(k)=NINT(1000*somx(k,u,v,w))
              end do
              if(prte) then !Print only at the end
                IF(iwmx /= 1) then
                   Write(lpt,"(a,50(1X,f6.3,2X))") " Mtot",(tmom(n,u,v,w),n=npi,npf)
                   Write(lpt,"(5X,50(3X,i4,2X))")          (isigm(n),n=npi,npf)
                   Write(itto,"(a,50(1X,f6.3,2X))")" Mtot",(tmom(n,u,v,w),n=npi,npf)
                   Write(itto,"(5X,50(3X,i4,2X))")         (isigm(n),n=npi,npf)
                end if
                write(itto,"(a)")
                write(lpt,"(a)")
              end if

              DO l=npi,npf
                amx=amx+comx(l,u,v,w)
                amy=amy+somy(l,u,v,w)
                amz=amz+somz(l,u,v,w)
              END DO

            END IF

          END DO
        END DO
      END DO

      RETURN

      END SUBROUTINE out_conf


      SUBROUTINE bound_cond(n,k,u,v,w)
      Use mcm_inc
      implicit none
      integer, intent(in) :: n,k,u,v,w


!***********************************************************************
!  This subroutine implements the boundary conditions for spins at
!  the edges of the sample. Three kinds of boundary conditions are
!  available:
!
!  - free edges:
!    ----------
!    spins at any edge of the sample are considered as free of nei-
!    ghbours outside of the sample. The sample is then representa-
!    tive of a real small crystallite delimited by surfaces
!
!  - periodic boundary conditions (p.b.c.):
!    ----------------------------
!    spins at one edge of the sample are connected to spins at the
!    opposite edge of the sample. This condition imposes a 3D perio-
!    dicity to the system. Then the system may be thought of as an
!    infinite size crystal with an elementary cell equal to the
!    sample dimension.
!
!  - mixed boundary conditions:
!    -------------------------
!    this new kind of boundary conditions has been introduced to
!    simulate helimagnets with a propagation axis of the helix.
!    It corresponds to free edges along one direction (axis W) and
!    p.b.c along the other two directions.
!
!     This subroutine calls no subroutine.
!***********************************************************************
      nu(k)=u+m(n,k,1)
      nv(k)=v+m(n,k,2)
      nw(k)=w+m(n,k,3)
      IF(ianb < 0) THEN
        GO TO    97
      ELSE IF (ianb == 0) THEN
        GO TO   199
      ELSE
        GO TO    99
      END IF

!   * MIXED BOUNDARY CONDITIONS   *

      97   IF(nu(k) > 0) THEN
        GO TO   910
      END IF
      900     nu(k)=nu(k)+lu
      GO TO 930
      910     IF(lu-nu(k) < 0) THEN
        GO TO   920
      ELSE
        GO TO   930
      END IF
      920       nu(k)=nu(k)-lu
      930       IF(nv(k) > 0) THEN
        GO TO   950
      END IF
      940         nv(k)=nv(k)+lv
      GO TO 970
      950         IF(lv-nv(k) < 0) THEN
        GO TO   960
      ELSE
        GO TO   970
      END IF
      960           nv(k)=nv(k)-lv
      970           IF(nw(k) > 0) THEN
        GO TO   990
      END IF
      980             nw(k)=nwmax+1
      GO TO 911
      990             IF(lw-nw(k) < 0) THEN
        GO TO   901
      ELSE
        GO TO   911
      END IF
      901               nw(k)=nwmax+1
      911               GO TO 98

!   * PERIODIC BOUNDARY CONDITIONS   *

      99 IF(nu(k) > 0) THEN
        GO TO   110
      END IF
      100     nu(k)=nu(k)+lu
      GO TO 130
      110     IF(lu-nu(k) < 0) THEN
        GO TO   120
      ELSE
        GO TO   130
      END IF
      120       nu(k)=nu(k)-lu
      130       IF(nv(k) > 0) THEN
        GO TO   150
      END IF
      140         nv(k)=nv(k)+lv
      GO TO 170
      150         IF(lv-nv(k) < 0) THEN
        GO TO   160
      ELSE
        GO TO   170
      END IF
      160           nv(k)=nv(k)-lv
      170           IF(nw(k) > 0) THEN
        GO TO   190
      END IF
      180             nw(k)=nw(k)+lw
      GO TO 210
      190             IF(lw-nw(k) < 0) THEN
        GO TO   200
      ELSE
        GO TO   210
      END IF
      200               nw(k)=nw(k)-lw
      210               GO TO 98

!   * FREE EDGES   *

      199   IF(nu(k) > 0) THEN
        GO TO   310
      END IF
      300     nu(k)=numax+1
      GO TO 330
      310     IF(lu-nu(k) < 0) THEN
        GO TO   320
      ELSE
        GO TO   330
      END IF
      320       nu(k)=numax+1
      330       IF(nv(k) > 0) THEN
        GO TO   350
      END IF
      340         nv(k)=nvmax+1
      GO TO 370
      350         IF(lv-nv(k) < 0) THEN
        GO TO   360
      ELSE
        GO TO   370
      END IF
      360           nv(k)=nvmax+1
      370           IF(nw(k) > 0) THEN
        GO TO   390
      END IF
      380            nw(k)=nwmax+1
      GO TO 98
      390             IF(lw-nw(k) < 0) THEN
        GO TO   400
      ELSE
        GO TO    98
      END IF
      400               nw(k)=nwmax+1
      98               CONTINUE

      RETURN

      END SUBROUTINE bound_cond

!***********************************************************************
!       General description of subroutines ENn_CALC(N,U,V,W)
!
!  The subroutines ENn_CALC(N,U,V,W) calculate the energy difference
!  between the new randomly generated configuration (NEW) and the old
!  configuration (OLD), expressed as DELTA = E(OLD) - E(NEW)
!  according to the energy formula :
!  (with the convention that a vector V' represents the transposed
!  form of a vector V)
!
!  E = - 0.5*Si'[Jij]Sj + Di(ui'Si)**2) - H'Si
!                                              AJXX AJXY AJXZ
!  [Jij] is the exchange tensor with elements  AJYX AJYY AJYZ
!                                              AJZX AJZY AJZZ
!
!  The exchange term Si'[Jij]Sj can be written Si'Jsi where
!  Jsi = [Jij]Sj is the vector (AJSX,AJSY,AJSZ).
!  The calculation of DELTA is made in two steps:
!      - first is calculated the term Jsi, which corresponds to the
!        microscopic molecular field at site i, and that is obviously
!        the same for the two configurations of spin i
!      _ then DELTA difference is calculated according to the above
!        formula.
!
!  The four subroutines ENn_CALC(N,U,V,W) concern different spin
!  models, namely:
!      the Ising model                       with  EN1_CALC(N,U,V,W)
!      the XY model                          with  EN2_CALC(N,U,V,W)
!      the Heisenberg model                  with  EN3_CALC(N,U,V,W)
!      the q-state Potts planar (XY) model   with  ENq_CALC(N,U,V,W)
!  Only one of them is activated at every run of the program.
!
!  These subroutines call the following subroutines:
!     BOUND_COND
!     RANnD(RADIUS,X1,Y1,,)
!***********************************************************************


      SUBROUTINE en1_calc(n,u,v,w)
      Use mcm_inc
      implicit none
      integer, intent(in) :: n,u,v,w

!***********************************************************************
!  See the general description of subroutines ENn_CALC(N,U,V,W)
!  just above.
!  EN1_CALC(N,U,V,W) deals with Ising spins.
!***********************************************************************
      integer :: i,k

      ajsx(n)=0
      ajsy(n)=0
      ajsz(n)=0
      DO k=1,nvs(n)
        ajxx(n,k)=ajjx(n,k)
        ajyy(n,k)=ajjy(n,k)
        ajzz(n,k)=ajjz(n,k)

        IF ((ansmod == 'Y').OR.(ansmod == 'y')) THEN


!*****    LOCAL MODIFICATION OF COUPLING (IF ANY)    ******


          DO i=1,im
            IF ( (w == iwm(i)).AND.(v == ivm(i)).AND.(u == ium(i))  &
                  .AND.(n == nm(i)).AND.(k == numb(i)) ) THEN
              ajxx(n,k)=ajmod(i)
              ajyy(n,k)=ajmod(i)
              ajzz(n,k)=ajmod(i)
            END IF
          END DO


!*****    END OF MODIFICATION    *****


        END IF

        CALL bound_cond(n,k,u,v,w)


!*****    COMPUTE DELTA    *****

        IF ( (nav(n,k) /= n).OR.(nu(k) /= u).OR.(nv(k) /= v)  &
              .OR.(nw(k) /= w) ) THEN
          ajsx(n)=ajsx(n)+ajxx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajxy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajxz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))
          ajsy(n)=ajsy(n)+ajyy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajyx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajyz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))
          ajsz(n)=ajsz(n)+ajzz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajzx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajzy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))
        END IF
      END DO

      IF (nfip /= 0) THEN
        CALL ran1d (amps2(n),x,y,z,n)
      END IF

      delta = (sx(n,u,v,w)-x)*(ajsx(n)+(nmult*hx))  &
          +(sy(n,u,v,w)-y)*(ajsy(n)+(nmult*hy)) +(sz(n,u,v,w)-z)*(ajsz(n)+(nmult*hz))  &
          -nmult*d(n) * ( ( dx(n)*sx(n,u,v,w)  &
          +dy(n)*sy(n,u,v,w) +dz(n)*sz(n,u,v,w) )**2  &
          -( dx(n)*x + dy(n)*y + dz(n)*z )**2 )

      RETURN

      END SUBROUTINE en1_calc



      SUBROUTINE en2_calc(n,u,v,w)
      Use mcm_inc
      implicit none
      integer, intent(in) :: n,u,v,w
!***********************************************************************
!  See the general description of subroutines EN*_CALC(N,U,V,W)
!  just before the subroutine EN1_CALC(N,U,V,W).
!  EN2_CALC(N,U,V,W) concerns XY spins.
!***********************************************************************
      integer :: k,i
      ajsx(n)=0.0
      ajsy(n)=0.0
      ajsz(n)=0.0

      DO k=1,nvs(n)
        ajxx(n,k)=ajjx(n,k)
        ajyy(n,k)=ajjy(n,k)
        ajzz(n,k)=ajjz(n,k)


!*****    LOCAL MODIFICATION OF COUPLING (IF ANY)    ******


        IF ((ansmod == 'Y').OR.(ansmod == 'y')) THEN
          DO i=1,im
            IF ( (w == iwm(i)).AND.(v == ivm(i)).AND.(u == ium(i))  &
                  .AND.(n == nm(i)).AND.(k == numb(i)) ) THEN
              ajxx(n,k)=ajmod(i)
              ajyy(n,k)=ajmod(i)
              ajzz(n,k)=ajmod(i)
            END IF
          END DO


!*****    END OF MODIFICATION    *****


        END IF

        CALL bound_cond(n,k,u,v,w)


!*****    COMPUTE DELTA    *****

        IF ( (nav(n,k) /= n).OR.(nu(k) /= u).OR.(nv(k) /= v)  &
              .OR.(nw(k) /= w) ) THEN
          ajsx(n)=ajsx(n)+ajxx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajxy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajxz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))
          ajsy(n)=ajsy(n)+ajyy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajyx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajyz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))
          ajsz(n)=ajsz(n)+ajzz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajzx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajzy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))
        END IF
      END DO

      IF (nfip /= 0) THEN
        CALL ran2d(amps2(n),x,y)
      END IF

      delta = (sx(n,u,v,w)-x)*(ajsx(n)+(nmult*hx))  &
          +(sy(n,u,v,w)-y)*(ajsy(n)+(nmult*hy)) +(sz(n,u,v,w)-z)*(ajsz(n)+(nmult*hz))  &
          -nmult*d(n) * ( ( dx(n)*sx(n,u,v,w)  &
          +dy(n)*sy(n,u,v,w) +dz(n)*sz(n,u,v,w) )**2  &
          -( dx(n)*x + dy(n)*y + dz(n)*z )**2 )

      RETURN
      END SUBROUTINE en2_calc



      SUBROUTINE en3_calc(n,u,v,w)
      Use mcm_inc
      implicit none
      integer, intent(in) :: n,u,v,w
!***********************************************************************
!  See the general description of subroutines EN*_CALC(N,U,V,W)
!  just before the subroutine EN1_CALC(N,U,V,W).
!  EN3_CALC(N,U,V,W) concerns Heisenberg spins.
!***********************************************************************
      integer :: k,i
      ajsx(n)=0
      ajsy(n)=0
      ajsz(n)=0

      DO k=1,nvs(n)
        ajxx(n,k)=ajjx(n,k)
        ajyy(n,k)=ajjy(n,k)
        ajzz(n,k)=ajjz(n,k)


!*****    LOCAL MODIFICATION OF COUPLING (IF ANY)    ******


        IF ((ansmod == 'Y').OR.(ansmod == 'y')) THEN
          DO i=1,im
            IF ( (w == iwm(i)).AND.(v == ivm(i)).AND.(u == ium(i))  &
                  .AND.(n == nm(i)).AND.(k == numb(i)) ) THEN
              ajxx(n,k)=ajmod(i)
              ajyy(n,k)=ajmod(i)
              ajzz(n,k)=ajmod(i)
            END IF
          END DO


!*****    END OF MODIFICATION    *****


        END IF

        CALL bound_cond(n,k,u,v,w)


!*****    COMPUTE DELTA    *****

        IF ( (nav(n,k) /= n).OR.(nu(k) /= u).OR.(nv(k) /= v)  &
              .OR.(nw(k) /= w) ) THEN
          ajsx(n)=ajsx(n)+ajxx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajxy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajxz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))
          ajsy(n)=ajsy(n)+ajyy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajyx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajyz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))
          ajsz(n)=ajsz(n)+ajzz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajzx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajzy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))
        END IF
      END DO

      IF (nfip /= 0) THEN
        CALL ran3d (amps2(n),x,y,z)
      END IF

      delta = (sx(n,u,v,w)-x)*(ajsx(n)+(nmult*hx))  &
          +(sy(n,u,v,w)-y)*(ajsy(n)+(nmult*hy)) +(sz(n,u,v,w)-z)*(ajsz(n)+(nmult*hz))  &
          -nmult*d(n) * ( ( dx(n)*sx(n,u,v,w)  &
          +dy(n)*sy(n,u,v,w) +dz(n)*sz(n,u,v,w) )**2  &
          -( dx(n)*x + dy(n)*y + dz(n)*z )**2 )

      RETURN
      END SUBROUTINE en3_calc


      SUBROUTINE enq_calc(n,u,v,w)
!***********************************************************************
!  See the general description of subroutines EN*_CALC(N,U,V,W)
!  just before the subroutine EN1_CALC(N,U,V,W).
!  ENq_CALC(N,U,V,W) concerns q-state Potts planar (XY) spins.
!***********************************************************************
      Use mcm_inc
      implicit none
      integer, intent(in) :: n,u,v,w
      integer :: k,i

      ajsx(n)=0
      ajsy(n)=0
      ajsz(n)=0

      DO k=1,nvs(n)
        ajxx(n,k)=ajjx(n,k)
        ajyy(n,k)=ajjy(n,k)
        ajzz(n,k)=ajjz(n,k)


!*****    LOCAL MODIFICATION OF COUPLING (IF ANY)    ******


        IF ((ansmod == 'Y').OR.(ansmod == 'y')) THEN
          DO i=1,im
            IF ( (w == iwm(i)).AND.(v == ivm(i)).AND.(u == ium(i))  &
                  .AND.(n == nm(i)).AND.(k == numb(i)) ) THEN
              ajxx(n,k)=ajmod(i)
              ajyy(n,k)=ajmod(i)
              ajzz(n,k)=ajmod(i)
            END IF
          END DO


!*****    END OF MODIFICATION    *****


        END IF

        CALL bound_cond(n,k,u,v,w)


!*****    COMPUTE DELTA    *****

        IF ( (nav(n,k) /= n).OR.(nu(k) /= u).OR.(nv(k) /= v)  &
              .OR.(nw(k) /= w) ) THEN
          ajsx(n)=ajsx(n)+ajxx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajxy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajxz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))
          ajsy(n)=ajsy(n)+ajyy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajyx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajyz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))
          ajsz(n)=ajsz(n)+ajzz(n,k)*sz(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajzx(n,k)*sx(nav(n,k),nu(k),nv(k),nw(k))  &
              +ajzy(n,k)*sy(nav(n,k),nu(k),nv(k),nw(k))
        END IF

      END DO

      IF (nfip /= 0) THEN
        CALL ranqd (amps2(n),x,y)
      END IF

      delta = (sx(n,u,v,w)-x)*(ajsx(n)+(nmult*hx))  &
          +(sy(n,u,v,w)-y)*(ajsy(n)+(nmult*hy)) +(sz(n,u,v,w)-z)*(ajsz(n)+(nmult*hz))  &
          -nmult*d(n) * ( ( dx(n)*sx(n,u,v,w)  &
          +dy(n)*sy(n,u,v,w) +dz(n)*sz(n,u,v,w) )**2  &
          -( dx(n)*x + dy(n)*y + dz(n)*z )**2 )

      RETURN

      END SUBROUTINE enq_calc
