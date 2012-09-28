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
!!---- MODULE: CFML_Propagation_Vectors
!!----   INFO: Series of procedures handling operation with Propagation
!!----         vectors
!!----
!!---- HISTORY
!!----    Update: 06/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,                only: Cp
!!--++    Use CFML_Math_General,              only: Zbelong
!!--++    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type
!!--++    Use CFML_Reflections_Utilities,     only: Hkl_R, Hkl_Equal
!!----
!!---- VARIABLES
!!--++    EPS                      [Private]
!!----    GROUP_K_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       HK_EQUIV
!!----       K_EQUIV
!!----       K_EQUIV_MINUS_K
!!----
!!----    Subroutines:
!!----       K_STAR
!!----       WRITE_GROUP_K
!!----
!!
 Module CFML_Propagation_Vectors

    !---- Use Modules ----!
    Use CFML_GlobalDeps,                only: Cp, Eps
    Use CFML_Math_General,              only: Zbelong
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type
    Use CFML_Reflections_Utilities,     only: Hkl_R, Hkl_Equal

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public :: Hk_Equiv, K_Equiv , K_Equiv_Minus_K

    !---- List of public subroutines ----!
    public :: K_Star, Write_Group_K

    !---- Definitions ----!

    !!----
    !!---- TYPE :: GROUP_K_TYPE
    !!--..
    !!---- Type, Public :: Group_K_Type
    !!----    type(Space_Group_Type)        :: G0             !Initial space group
    !!----    integer                       :: ngk            !Number of elements of G_k
    !!----    logical                       :: k_equiv_minusk !true if k equiv -k
    !!----    integer,      dimension(192)  :: p              !Pointer to operations of G0 that changes/fix k
    !!----                                                    !First indices: G_k, last indices: Stark
    !!----    integer                       :: nk             !Number of star arms
    !!----    real(kind=cp),dimension(3,24) :: stark          !Star of the wave vector k
    !!---- End Type Group_K_Type
    !!----
    !!--<<
    !!----    The integer pointer p is used as follows:
    !!----    If we defined the object G as ->  type(Group_K_Type) :: G
    !!----    G%p(1:ngk) gives the numeral of the symmetry operators of G%G0
    !!----    belonging to G_k.
    !!----    G%p(192:193-nk) gives the numeral of the the symmetry operators of G%G0 that
    !!----    transform the initial k-vector to the other arms of the star.
    !!----    G%co(:,kk) gives also the numerals of the the symmetry operators of G%G0 that
    !!----    transform the initial k-vector to the arm kk of the star to the representative
    !!----    of the coset decomposition of G%G0 with respect to G_k.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Type, Public :: Group_k_Type
       type(Space_Group_Type)           :: G0
       integer                          :: ngk=0
       logical                          :: k_equiv_minusk=.false.
       integer, dimension(192)          :: p=0
       integer, dimension(48,48)        :: co=0
       integer                          :: nk=0
       real(kind=cp),dimension(3,48)    :: stark=0.0
    End Type Group_k_Type

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!----
    !!---- Logical Function  Hk_Equiv(H,K, Spacegk ,Friedel)
    !!----    real(kind=cp), dimension(3), intent(in) :: h
    !!----    real(kind=cp), dimension(3), intent(in) :: k
    !!----    Type (Group_k_Type),         intent(in) :: SpaceGk
    !!----    logical, optional,           intent(in) :: Friedel
    !!----
    !!----    Calculate if two real reflections are equivalent
    !!----
    !!---- Update: February - 2005
    !!
    Function Hk_Equiv(H,K,Spacegk,Friedel) Result (Info)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent (in) :: h,k
       Type (Group_k_Type),         intent (in) :: SpaceGk
       logical, optional,           intent(in)  :: Friedel
       logical                                  :: info

       !---- Local Variables ----!
       integer                    :: i,j
       real(kind=cp), dimension(3):: hh

       info=.false.
       do i=1, SpaceGk%ngk
          j=SpaceGk%p(i)
          hh = hkl_R(h,SpaceGk%g0%SymOp(j))
          if (hkl_equal(k,hh)) then
             info=.true.
             exit
          end if
          if (present(Friedel) .and. SpaceGk%g0%centred /= 2 .and. SpaceGk%k_equiv_minusk) then
             if (Friedel) then
                if (hkl_equal(k,-hh)) then
                   info=.true.
                   exit
                end if
             end if
          end if
       end do

       return
    End Function Hk_Equiv

    !!----
    !!---- Logical Function  K_Equiv(H,K,Latyp)  Result (Info)
    !!----    real(kind=cp), dimension(3),    intent (in) :: h      !  In ->
    !!----    real(kind=cp), dimension(3),    intent (in) :: k      !  In ->
    !!----    character (len=*),              intent (in) :: latyp  !  In ->
    !!----
    !!----    Calculate if two k-vectors are equivalent in the sense
    !!----    that "h" is equivalent to "k" if "h-k" belongs to the reciprocal
    !!----    lattice. Only lattice type is needed.
    !!----
    !!---- Update: February - 2005
    !!
    Function K_Equiv(H,K,Latyp) Result (Info)
       !---- Arguments ----!
       real(kind=cp), dimension(3),intent (in) :: h,k
       character (len=*),          intent (in) :: latyp
       logical                                 :: info

       !---- Local Variables ----!
       character(len=1)             :: latypc
       real(kind=cp), dimension(3)  :: vec
       integer, dimension(3)        :: ivec

       info=.false.
       vec=h-k
       if (.not. Zbelong(vec) ) return
       ivec=nint(vec)
       latypc=adjustl(latyp)
       select case (Latypc)
          case("P","p")
             info=.true.
          case("A","a")
             if (mod(ivec(2)+ivec(3),2) == 0) info=.true.
          case("B","b")
             if (mod(ivec(1)+ivec(3),2) == 0) info=.true.
          case("C","c")
             if (mod(ivec(1)+ivec(2),2) == 0) info=.true.
          case("I","i")
             if (mod(ivec(1)+ivec(2)+ivec(3),2) == 0) info=.true.
          case("R","r")
             if (mod(-ivec(1)+ivec(2)+ivec(3),3) == 0) info=.true.
          case("F","f")
             if (mod(ivec(2)+ivec(3),2) == 0              &
                 .and. mod(ivec(1)+ivec(3),2) == 0        &
                 .and. mod(ivec(1)+ivec(2),2) == 0 ) info=.true.
       end select

       return
    End Function K_Equiv

    !!----
    !!---- Logical Function K_Equiv_Minus_K(Vec,Lat) Result(Equiv)
    !!----    real(kind=cp), dimension(3), intent(in) :: vec      !  In ->
    !!----    character (len=*),           intent(in) :: Lat      !  In ->
    !!----    logical                                 :: equiv
    !!----
    !!----    Determine whether a k-vector is equivalent to -k.
    !!----
    !!---- Update: February - 2005
    !!
    Function K_Equiv_Minus_K(Vec,Lat) result(equiv)
       !---- Argument ----!
       real(kind=cp), dimension(3), intent(in) :: vec
       character(len=*),            intent(in) :: Lat
       logical                                 :: equiv

       !---- Local variables ----!
       character(len=1)      :: latc
       logical               :: integralv
       integer, dimension(3) :: ivec

       integralv=.false.

       if ( Zbelong(2.0*vec) )  integralv=.true.
       equiv=.false.
       latc=adjustl(lat)

       if (integralv) then
          ivec = nint (2.0 * vec)
          select case (latc)
             case("P","p")
                equiv=.true.

             case("A","a")
                if (mod(ivec(2)+ivec(3),2) == 0) equiv=.true.

             case("B","b")
                if (mod(ivec(1)+ivec(3),2) == 0) equiv=.true.

             case("C","c")
                if (mod(ivec(1)+ivec(2),2) == 0) equiv=.true.

             case("I","i")
                if (mod(ivec(1)+ivec(2)+ivec(3),2) == 0) equiv=.true.

             case("R","r")
                if (mod(-ivec(1)+ivec(2)+ivec(3),3) == 0) equiv=.true.

             case("F","f")
                if (mod(ivec(2)+ivec(3),2) == 0        &
                   .and. mod(ivec(1)+ivec(3),2) == 0        &
                   .and. mod(ivec(1)+ivec(2),2) == 0 ) equiv=.true.

             case default
                equiv=.true.

          end select
       end if

       return
    End Function K_Equiv_Minus_K

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine K_Star(K,Spacegroup,Gk)
    !!----    integer, dimension(3),   intent(in)  :: k           !  In ->
    !!----    type (Space_Group_Type), intent(in)  :: SpaceGroup  !  In ->
    !!----    Type (Group_k_Type),     intent(out) :: Groupk      ! Out ->
    !!----
    !!----    Calculate the star of the propagation vector
    !!----    and the group of the vector k. Construct the object
    !!----    "Groupk" that has as a component SpaceGroup + Pointers to the
    !!----    operators belonging to Gk and the star of k.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine K_Star(K,Spacegroup,Gk)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent (in) :: k
       Type (Space_Group_Type),     intent (in) :: SpaceGroup
       Type (Group_k_Type),         intent(out) :: Gk

       !---- Local Variables ----!
       real(kind=cp), dimension(3):: h
       integer, dimension(48)     :: ind=1
       integer                    :: i, j, ng

       ng=SpaceGroup%numops
       Gk%g0=SpaceGroup  !<- Copy the whole space group to the component G0 of Gk
       Gk%p(1)= 1    !<- pointer to the identity that always belong to G_k
       Gk%p(192)= 1  !<- pointer to the identity that always is a coset representative of G0[Gk]
       Gk%nk  = 1
       Gk%ngk = 1
       Gk%stark(:,1) = k  !<- First arm of the star of k
       Gk%k_equiv_minusk = .true. !it is supposed that k equiv -k
       if (SpaceGroup%Centred /= 1) then
          j=ng+1
          h = hkl_R(k,SpaceGroup%SymOp(j))
          !---- k not equivalent to -k
          if (.not. k_EQUIV(h,k, SpaceGroup%SPG_lat))  Gk%k_equiv_minusk = .false.
          ng=ng*2
       else
          if (.not. k_EQUIV(k,-k, SpaceGroup%SPG_lat))  Gk%k_equiv_minusk = .false.
       end if

       do_ng: do i=2,ng
          h = hkl_R(k,SpaceGroup%SymOp(i))

          !---- this operation belong to Gk
          if (k_EQUIV(h,k, SpaceGroup%SPG_lat)) then
             Gk%ngk = Gk%ngk +1
             Gk%p(Gk%ngk)      = i
             Gk%co(Gk%ngk,1)   = i
             cycle do_ng
          end if
          !---- Passing here means that h is not equivalent to k, so there is eventually a
          !-----new k-vector (=h) to be added to the list of the star.
          !---- Look if the generated h-vector is already equivalent to
          !---- one of the list
          do j=1, Gk%nk
             if (k_EQUIV(h,Gk%stark(:,j), SpaceGroup%SPG_lat) ) then
                ind(j)=ind(j)+1
                Gk%co(ind(j),j)   = i
                cycle do_ng
             end if
          end do
          !---- Passing here is to add a new vector to the list
          Gk%nk = Gk%nk +1           !this construct the star of k
          Gk%stark(:,Gk%nk) = h
          Gk%p(193-Gk%nk)   = i        !This points to operators changing k
          Gk%co(1,Gk%nk)    = i
       end do do_ng

       return
    End Subroutine K_Star

    !!----
    !!---- Subroutine Write_Group_K(Gk,Lun)
    !!----    Type (Group_k_Type),   intent(in) :: Gk    !  In ->
    !!----    integer, optional,     intent(in) :: lun   !  In ->
    !!----
    !!----    Subroutine to write the operators of the propagation vector
    !!----    group and the list of all vectors {k}, belonging to the star
    !!----    of k.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Group_k(Gk,lun)
       !---- Arguments ----!
       Type (Group_k_Type),   intent(in) :: Gk
       Integer, optional,     intent(in) :: lun

       !---- Local Variables ----!
       real(kind=cp), dimension(3):: h,k
       character(len=4)           :: kvec
       character(len=22)          :: formato
       integer                    :: i, j, lu, l, m

       lu=6
       if (present(lun)) lu=lun

       write(unit=lu,fmt="(/,a)") "      PROPAGATION VECTOR GROUP INFORMATION"
       write(unit=lu,fmt="(a,/)") "      ===================================="

       h=Gk%stark(:,1)
       write(unit=lu,fmt="(a,3f8.4,a)") " => The input propagation vector is: K=(",h," )"
       if (Gk%k_equiv_minusk) then
          write(unit=lu,fmt="(a,i3,a)")    " => K .. IS .. equivalent to -K "
       else
          write(unit=lu,fmt="(a,i3,a)")    " => K .. IS NOT .. equivalent to -K "
       end if
       write(unit=lu,fmt="(a)")       " => The operators following the k-vectors constitute the co-set decomposition G[Gk]"
       write(unit=lu,fmt="(a)")       "    The list of equivalent k-vectors are also given on the right of operators.     "
       write(unit=lu,fmt="(a,i3,a,/)")" => The star of K is formed by the following ",Gk%nk," vectors:"

       m=0
       do i=1,Gk%G0%multip
          if(len_trim(Gk%G0%SymopSymb(i)) > m ) m=len_trim(Gk%G0%SymopSymb(i))
       end do
       formato="(a,i3,a,a  ,a,3f8.4,a)"
       write(unit=formato(10:11),fmt="(i2)") m
       do i=1, Gk%nk
          j=Gk%p(193-i)
          kvec="    "
          if (i < 10) then
             write(unit=kvec,fmt="(a,i1)")"k_",i
          else
             write(unit=kvec,fmt="(a,i2)")"k_",i
          end if
          k= Gk%stark(:,i)
          if (.not. k_EQUIV(-h,k, Gk%G0%SPG_lat)) then
             write(unit=lu,fmt="(/,a,3f8.4,a,i3,a,a)") "           "//kvec//" = (",k,  &
                   " )    Op: (",j,") ",trim(Gk%G0%SymopSymb(j))
          else
             write(unit=lu,fmt="(/,a,3f8.4,a,i3,a,a)") "  Eqv. -K: "//kvec//" = (",k,  &
                   " )    Op: (",j,") ",trim(Gk%G0%SymopSymb(j))
          end if
          do l=2,48
             j=Gk%co(l,i)
             if (j == 0) exit
             k = hkl_R(h,Gk%G0%SymOp(j))
             write(unit=lu,fmt=formato)"                                                 Op: (",  &
                                        j,") ",     Gk%G0%SymopSymb(j) ,"  -> ( ",k," )"
          end do
       end do

       write(unit=lu,fmt="(/,a,/)")    " => G_k has the following symmetry operators: "

       do i=1, Gk%ngk
          j=Gk%p(i)
          write(unit=lu,fmt="(2(a,i3),a,a)") "    ",i," SYMM(",j,") = ", trim(Gk%G0%SymopSymb(j))
       end do

       if (.not. Gk%k_equiv_minusk .and. Gk%G0%Centred == 1) then
          write(unit=lu,fmt="(/,a)")    " => The original space group is non-centrosymmetric and K is not equivalent to -K"
          write(unit=lu,fmt="(a  )")    "    If the vector k=-K does not appear within the arm of the K-vector, the propagation"
          write(unit=lu,fmt="(a  )")    "    vector k=-K should be included in applications where final vector components have "
          write(unit=lu,fmt="(a,/)")    "    to be real values.   "
       end if

       return
    End Subroutine Write_Group_k

 End Module CFML_Propagation_Vectors

