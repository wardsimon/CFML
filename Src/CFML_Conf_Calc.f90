!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2018  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!----               Nebil A. Katcho    (ILL)
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
!!---- MODULE: CFML_BVS_Energy_Calc
!!----   INFO: Subroutines related to calculations of energy or
!!----         configuration properties depending on the crystal structure: BVS, Energy,....
!!----
!!---- HISTORY
!!----    Updated: 16/12/2014, 20/02/2018
!!----
!!----
!!---- DEPENDENCIES
!!----
!!----
!!---- VARIABLES
!!----    ATOM_CONF_LIST
!!----    ERR_CONF
!!----    ERR_CONF_MESS
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ALLOCATE_ATOMS_CONF_LIST
!!--++       BOND_VALENCE                  [Private]
!!----       CALC_BVS
!!----       CALC_MAP_BVEL
!!----       CALC_MAP_BVS
!!--++       COMPLETE_TABLE                [Private]
!!----       COST_BVS
!!----       COST_BVS_COULOMBREP
!!----       DEALLOCATE_ATOMS_CONF_LIST
!!----       EWALD
!!----       INIT_ERR_CONF
!!----       SET_TABLE_D0_B
!!----       SPECIES_ON_LIST
!!----       SET_FORMAL_CHARGES
!!----
!!
 Module CFML_BVS_Energy_Calc
    !---- Use Files ----!
    Use CFML_GlobalDeps,                 only: Sp,Cp,dp,pi,tpi
    Use CFML_Math_3D,                    only: cross_product
    Use CFML_Math_General,               only: Sort_Strings,Sort,cosd,set_epsg
    use CFML_String_Utilities,           only: Getword, U_Case,pack_string, get_logunit
    Use CFML_Scattering_Chemical_Tables, only: Get_Ionic_Radius, Get_Chemsymb, Get_Covalent_Radius, Set_Chem_Info
    use CFML_Crystal_Metrics,            only: Crystal_Cell_Type
    use CFML_Crystallographic_Symmetry,  only: Space_Group_Type, get_orbit
    use CFML_Atom_TypeDef,               only: Atom_type, Init_Atom_type, Write_Atom_List, Atom_list_Type, Allocate_Atom_List, &
                                               Deallocate_Atom_List, AtList1_ExtenCell_AtList2, Atoms_Cell_Type, Allocate_Atoms_Cell, &
                                               Deallocate_Atoms_Cell
    use CFML_Geometry_Calc,              only: Coord_Info, Distance, calc_dist_angle_sigma, calc_dist_angle
    use CFML_Export_VTK,                 only: write_grid_VESTA
    use CFML_BVSpar

    !---- Variables ----!
    implicit none

    private

    !---- List of public subroutines ----!
    public :: Allocate_Atoms_Conf_List, Calc_BVS, Deallocate_Atoms_Conf_List,           &
              Init_Err_Conf, Set_Table_d0_b, Species_on_List, Set_Table_BVEL_Params,    &
              Calc_Map_BVS, Cost_BVS, Cost_BVS_CoulombRep, Calc_Map_BVEL, Ewald,        &
              Set_Formal_Charges

    !---- List of public private ----!
    private :: Bond_Valence, Complete_Table, Complete_Table_BVEL, get_soft_covalent_radius

    !---- Definitions ----!

    !!----
    !!---- TYPE :: ATOMS_CONF_LIST_TYPE
    !!--..
    !!---- Type, public :: Atoms_Conf_List_Type
    !!----    integer                                     :: natoms    ! Total number of atoms in the list
    !!----    integer                                     :: N_Spec    ! Number of different species in the list
    !!----    integer                                     :: N_Anions  ! Number of anions in the list
    !!----    integer                                     :: N_Cations ! Number of cations in the list
    !!----    real(kind=cp)                               :: Tol       ! Tolerance(%) for sum of radii conditions
    !!----    real(kind=cp)                               :: totatoms  ! Total number of atoms in the unit cell
    !!----    character(len=4), dimension(:), allocatable :: Species   ! Symbol + valence
    !!----    real(kind=cp),    dimension(:), allocatable :: Radius    !ionic/atomic radius of species
    !!----    type(Atom_Type),  dimension(:), allocatable :: atom
    !!---- End Type Atoms_Conf_List_Type
    !!----
    !!---- Update: March - 2005
    !!
    Type, public :: Atoms_Conf_List_Type
       integer                                     :: natoms    ! Total number of atoms in the list
       integer                                     :: N_Spec    ! Number of different species in the list
       integer                                     :: N_Anions  ! Number of anions in the list
       integer                                     :: N_Cations ! Number of cations in the list
       real(kind=cp)                               :: Tol       ! Tolerance(%) for sum of radii conditions
       real(kind=cp)                               :: totatoms  ! Total number of atoms in the unit cell
       character(len=4), dimension(:), allocatable :: Species
       real(kind=cp),    dimension(:), allocatable :: Radius    !ionic/atomic radius of species
       type(Atom_Type),  dimension(:), allocatable :: Atom
    End type Atoms_Conf_List_Type

    !!----
    !!---- ERR_CONF
    !!----    logical, public  :: err_conf
    !!----
    !!----    Logical Variable taking the value .true. if an error in the module
    !!----    CONFIGURATIONS_CALCULATIONS occurs.
    !!----
    !!---- Update: March - 2005
    !!
    logical, public  :: err_conf

    !!----
    !!---- ERR_CONF_MESS
    !!----    character(len=150), public:: ERR_Conf_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: March - 2005
    !!
    character(len=150), public :: ERR_Conf_Mess

    !!----
    !!---- ERR_CHAR
    !!----    logical, public  :: err_char
    !!----
    !!----    Logical Variable taking the value .true. if an error  in formal charges subroutine
    !!----
    !!---- Update: July - 2016
    !!
    Logical, Public :: ERR_Char

    !!----
    !!---- ERR_CHAR_MESS
    !!----    character(len=150), public:: ERR_Char_Mess
    !!----
    !!----    String containing information about the last error in formal charges subroutine
    !!----
    !!---- Update: July - 2016
    !!
    Character(len=150), Public :: ERR_Char_Mess

    !!----
    !!---- WARN_CHAR
    !!----    logical, public  :: warn_char
    !!----
    !!----    Logical Variable taking the value .true. if a warning  in formal charges subroutine
    !!----
    !!---- Update: July - 2016
    !!
    Logical, Public :: WARN_Char

    !!----
    !!---- WARN_CHAR_MESS
    !!----    character(len=150), public:: ERR_Warn_Mess
    !!----
    !!----    String containing information about the last warning
    !!----
    !!---- Update: July - 2016
    !!
    Character(len=150), Public :: WARN_Char_Mess

 Contains

    !!----
    !!---- Subroutine Allocate_Atoms_Conf_List(N,A)
    !!----    integer, intent(in)                         :: n    !  In -> Atoms in asymmetric unit
    !!----    type (Atoms_Conf_List_Type), intent(in out) :: A    !  In -> Objet to be allocated
    !!----
    !!----    Allocation of objet A of type atom_list_Conf. This subroutine
    !!----    should be called before using an object of type atom_list_Conf.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_Atoms_Conf_List(n,A)
       !---- Arguments ----!
       integer,                     intent(in)       :: N  !N. atoms in asymmetric unit
       type (Atoms_Conf_List_Type), intent(in out)   :: A  !Objet to be allocated

       !---- Local variables ----!
       integer :: i

       A%natoms   = n
       A%n_spec   = 0
       A%n_anions = 0
       A%n_cations= 0
       A%tol      = 20.0 ! 20% tolerance by default
       A%totatoms = 0.0
       if(.not. allocated(A%atom)) allocate (A%atom(n))
       do i=1,n
          call init_atom_type(A%atom(i))
       end do

       return
    End Subroutine Allocate_Atoms_Conf_List

    !!--++
    !!--++ Subroutine Bond_Valence(d0,b0,d,sd,bv,sbv)
    !!--++    real(kind=cp), intent(in)  :: d0    ! BVS parameter
    !!--++    real(kind=cp), intent(in)  :: b0    ! BVS parameter
    !!--++    real(kind=cp), intent(in)  :: d     ! Distance
    !!--++    real(kind=cp), intent(in)  :: sd    ! Sigma(d)
    !!--++    real(kind=cp), intent(out) :: bv    ! BVS
    !!--++    real(kind=cp), intent(out) :: sbv   ! Sigma(bv)
    !!--++
    !!--++    (Private)
    !!--++    Zachariasen exponential expression of Bond Valence
    !!--++
    !!--++ Update: Dic - 2007
    !!
    Subroutine Bond_Valence(D0,B0,D,Sd,Bv,Sbv)
       !---- Arguments ----!
       real(kind=cp),  intent(in)  :: d0,b0  !Bond-valence parameters
       real(kind=cp),  intent(in)  :: d,sd   !Bond distance and sigma
       real(kind=cp),  intent(out) :: bv,sbv !Bond-valence and sigma

       bv=EXP((d0-d)/b0)
       sbv=bv*sd/b0

       return
    End Subroutine Bond_Valence

    !!----
    !!---- Subroutine Calc_BVS(A, Ipr, N_BVSm, BVS_M, Filecod)
    !!----    type (Atoms_Conf_List_type),              intent(in)   :: A            !  In -> Object of Atoms_Conf_List_type
    !!----    integer,                        optional, intent(in)   :: Ipr
    !!----    integer,                        optional, intent(in)   :: n_bvsm       !  In -> Number of modifications
    !!----    character(len=*), dimension(:), optional, intent(in)   :: BVS_M        ! In -> Text with BVS parameters
    !!----    character(len=*),               optional, intent(in)   :: Filecod
    !!----
    !!----    Subroutine to calculate Bond-Valence sums.
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Calc_Dist_Angle_Sigma" in order
    !!----    to update the internal private variables related to distance/angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    Control for error is present.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Calc_BVS(A, Ipr, N_Bvsm, Bvs_M, Filecod)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),            intent(in)  :: A      !  In -> Object of Atoms_Conf_List_type
       integer,                      optional, intent(in)  :: Ipr
       integer,                      optional, intent(in)  :: N_bvsm
       character(len=*),dimension(:),optional, intent(in)  :: bvs_m
       character(len=*),             optional, intent(in)  :: Filecod

       !---- Local variables ----!
       integer           :: i,j,k,n1,n2,icm,l,icn,isoc,kk,is0, &
                            isstr,isdav,isigt,ibvs,sig1,sig2
       real(kind=cp)     :: tol,fact,del2,s0,q2,  &
                            dd,sigtot,efcn,sums,dav,sdav,q1,d2,  &
                            str2,sstr,ric,r_2,del,perc,spred,disp,  &
                            str,rg1,dist,gii_a,gii_b,gii_c
       character(len=4)  :: rnatom

       call init_err_conf()

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20


       if (present(ipr)) then
          write(unit=ipr,fmt="(/,a)")   "  ------------------------------------------------"
          write(unit=ipr,fmt="(a)")     "  {--- BOND-VALENCE AND POLYHEDRA DISTORTIONS ---}"
          write(unit=ipr,fmt="(a,/)")   "  ------------------------------------------------"
       end if


       if (present(n_bvsm).and. present(ipr)) then
          write(unit=ipr,fmt="(a,/,a,/)")  &
          " Bond-Valence parameters (d0,B0) for Zachariasen formula:  s= exp{(d0-d)/B0}", &
                "   (read from external conditions)"
       else
          write(unit=ipr,fmt="(a,/,a,/)")  &
          " Bond-Valence parameters (d0,B0) for Zachariasen formula:  s= exp{(d0-d)/B0}", &
                "   (data read from internal table)"
       end if

       !---- Each line: Cation Anion d0 b
       if (present(n_bvsm) .and. present(bvs_m)) then
          call Complete_Table(A,N_bvsm,bvs_m)
       end if

       do n1=1,A%N_Cations
          do j=1,A%N_Anions
             n2=A%N_Cations+j
             if (present(ipr)) then
                write(unit=ipr,fmt="(2(a,i3,a,a4),2(a,f6.3),a)") &
                      "   Type",n1,": ",A%Species(n1)," with type",n2,": ",A%Species(n2),&
                      " d0=",Table_d0(n1,n2),"    B0=",Table_b(n1,n2), "   => Reference: "//trim(references(Table_ref(n1,n2)))
                write(unit=ipr,fmt="(2(a,a,a,f6.3,a),/)") &
                      "   Cation (Eff. radius): ",A%Species(n1),"(",A%Radius(n1),")   ", &
                      "   Anion  (Eff. radius): ",A%Species(n2),"(",A%Radius(n2),")"
             end if
          end do
       end do

       del2=0.0
       call Get_LogUnit(ibvs)

       if (present(filecod)) then
          open(unit=ibvs,file=trim(filecod)//"_sum.bvs",status="replace",action="write", position="append")
          write(unit=ibvs,fmt="(a)") "  Subroutine Calc_BVS (JRC-LLB, version: March-2005)"
          write(unit=ibvs,fmt="(a)") "  Title: Summary of Bond-Valence calculations for file: "//trim(filecod)//".cfl"
          write(unit=ibvs,fmt="(a)") "   Atom      Coord  D_aver Sigm   Distort(x10-4)    Valence    BVSum(Sigma)"
       else
          open(unit=ibvs,file="summary.bvs",status="replace",action="write", position="append")
          write(unit=ibvs,fmt="(a)") "  Subroutine Calc_BVS (JRC-LLB, version: March-2005)"
          write(unit=ibvs,fmt="(a)") "  Title: Summary of Bond-Valence calculations "
          write(unit=ibvs,fmt="(a)") "   Atom      Coord  D_aver Sigm   Distort(x10-4)    Valence    BVSum(Sigma)"
       end if

       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii_a=0.0
       gii_b=0.0
       gii_c=0.0
       do i=1,A%natoms
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0_cp,q1)
          icn=0
          efcn=0.0
          sums=0.0
          sigtot=0.0
          dav=0.0
          sdav=0.0
          str2=0.0
          d2=0.0
          isoc=INT(A%atom(i)%VarF(2)*1000.0+0.5)
          if (present(ipr)) then
             write(unit=ipr,fmt="(/,/,a,/,a,a4,a,f6.3,a,i3,a,/,a,/)") &
                  "    -------------------------------------------------------------------",  &
                  " => Bond-valence and coordination of atom: ",A%atom(i)%lab ," occupancy: ",A%atom(i)%VarF(1),"(",isoc,")",  &
                  "    -------------------------------------------------------------------"
          end if
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(j,i))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(j,i))%charge
             sig2= SIGN(1.0_cp,q2)
             if (sig1 == sig2) cycle
             dd=coord_info%dist(j,i)
             if (dd > (A%Radius(l)+A%Radius(k))*(1.0+tol)) cycle
             icn=icn+1
             efcn=efcn+A%Atom(coord_info%n_cooatm(j,i))%VarF(1) !/A%Atom(i)%VarF(1)
             kk=k
             d2=d2+dd*dd
             s0=coord_info%s_dist(j,i)
             dav=dav+dd
             sdav=sdav+s0*s0
             rnatom=A%atom(coord_info%n_cooatm(j,i))%lab
             call Bond_valence(Table_d0(l,k),Table_b(l,k),dd,s0,str,sstr)
             str=str*A%Atom(coord_info%n_cooatm(j,i))%VarF(1)             !/A%Atom(i)%VarF(1)
             sstr=sstr*A%Atom(coord_info%n_cooatm(j,i))%VarF(1)           !/A%Atom(i)%VarF(1)
             sums=sums+str
             str2=str2+str*str
             sigtot=sigtot+sstr*sstr
             is0=nint(s0*10000)
             isstr=nint(sstr*1000)
             if (present(ipr)) then
                write(unit=ipr,fmt="(a,a4,a,a4,a,f8.4,a,i4,a,f6.3,a,i3,a)")  &
                     " (",A%atom(i)%lab ,")-(",rnatom,") :",dd,"(",is0,")  ",str,"(",isstr,")"
             end if
          end do

          ric=real(icn)
          if (icn == 0) then
            if (present(ipr) ) then
                write(unit=ipr,fmt=*) " => Warning!! atom: ",A%atom(i)%lab ," is non-coordinated"
                write(unit=ipr,fmt=*) "    Increase the tolerance factor for ionic radii"
             end if
             cycle
          end if
          d2=d2/ric
          sigtot=SQRT(sigtot)
          dav=dav/ric
          sdav=SQRT(sdav)/ric

          isdav=INT(sdav*10000+0.5)
          isigt=INT(sigtot*1000+0.5)
          dist=10000.0*(d2/(dav*dav)-1.0)
          r_2=sums/ric
          r_2=SQRT(ABS(str2/ric-r_2*r_2))
          del=sums-ABS(q1)
          del2=del2+del*del
          perc=100.0*ABS(del/q1)

          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii_a=gii_a+abs(del)*fact
          gii_b=gii_b+perc*fact
          gii_c=gii_c+del*del*fact

          !efcn=efcn/A%Atom(i)%VarF(1)                    !Division by the site occupancy
          spred=ABS(q1)/efcn                              !Predicted valence of a single bond
          disp=table_d0(l,kk)-table_b(l,kk)*LOG(spred)    !Pred. distance
          if (present(ipr)) then
             write(unit=ipr,fmt="(/,a,i5,a,f5.2,a,a4)") " Coordination number: ",icn, &
                 "      Eff.Coor. number: ",efcn,"  for atom: ",A%atom(i)%lab
             write(unit=ipr,fmt="(a,f8.4,a,i4,a,f8.3,a)")  &
                  " Average distance  :",dav,"(",isdav,")  Distortion:  ",dist," xE-04"
             write(unit=ipr,fmt="(a,f8.4,a,f6.3)") " Predicted distance:",disp, &
                  "        Single bond-valence S=",spred
             write(unit=ipr,fmt="(a,f8.3,/,a,f8.3,a,i3,a)") &
                  "                                      Valence: ",q1,  &
                  "                                         Sums: ",sums, "(",isigt,")"
             write(unit=ipr,fmt="(a,2f8.3,/,a)") &
                  " Deviation from the Valence Sum Rule (r1,%dev):",del,perc,  &
                  " {r1=Sumj(sij)-Vi, %dev=100abs(r1)/Vi} "
             write(unit=ipr,fmt="(a,f8.3,/,a)")  &
                  " Deviation from the Equal Valence Rule    (r2):",r_2,  &
                  " {r2=<sij-<sij>>rms}"
             write(unit=ibvs,fmt="(tr4,a4,tr4,f6.2,f8.4,a,i4,a,f14.3,2f12.3,a,i3,a)") &
                  A%atom(i)%lab,efcn,dav,"(",isdav,")",dist,q1,sums, "(",isigt,")"

          end if
       end do

       rg1=SQRT(del2/real(A%natoms))*100.0
       gii_a=gii_a*100.0/A%totatoms
       gii_b=gii_b/A%totatoms  !*100.0 already multiplied
       gii_c=sqrt(gii_c/A%totatoms)*100.0

       if (present(ipr)) then
          write(unit=ipr,fmt="(/,6(a,/))")  &
             " => Lines concerning predicted average distances and single",  &
             "    bond-valence values, as well as the deviations from the",  &
             "    Equal Valence Rule,  apply only to those central atoms",  &
             "    having N coordination-atoms of the same chemical species. ",  &
             "    (The term 'single bond-valence' refers to the valence value",  &
             "    of a single bond for a regular polyhedron, so S=Valence/N)"
           write(unit=ipr,fmt="(/,4(a,/))")  &
             " => The Old Global Instability Index (GII) is calculated with the atoms of the asymetric unit (Num_Atoms).",&
             "    The normalized GII(a,b,c) below are calculated using the sum over asymmetric unit but multiplying ",&
             "    differences by the multiplicity of the site. N_Atoms_UCell is the total number of atoms in the ", &
             "    conventional unit cell. In all cases the result of the different expressions is multiplied by 100.0"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Old Global Instability Index (  GII=SQRT{SUM{|BVS-abs(q)|^2}/Num_Atoms} ) =", &
               rg1," /100"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(a)=       SUM {|BVS-abs(q)|  *mult}       /N_Atoms_UCell =", &
               gii_a," /100"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(b)=       SUM {|BVS-abs(q)|  *mult/abs(q)}/N_Atoms_UCell =", &
               gii_b," %"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(c)= SQRT{ SUM {|BVS-abs(q)|^2*mult}       /N_Atoms_UCell}=", &
               gii_c," /100"
       end if

       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Old Global Instability Index (  GII=SQRT{SUM{|BVS-abs(q)|^2}/Num_Atoms} ) =", &
            rg1," /100"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(a)=       SUM {|BVS-abs(q)|  *mult}       /N_Atoms_UCell =", &
            gii_a," /100"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(b)=       SUM {|BVS-abs(q)|  *mult/abs(q)}/N_Atoms_UCell =", &
            gii_b," %"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(c)= SQRT{ SUM {|BVS-abs(q)|^2*mult}       /N_Atoms_UCell}=", &
            gii_c," /100"

       call flush(ibvs)
       close (unit=ibvs)

       return
    End Subroutine Calc_BVS

    !!----
    !!---- Subroutine Calc_Map_BVEL(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta,vol,emin,npix,outp,bvel_map)
    !!----    !---- Arguments ----!
    !!----    type (Atoms_Conf_List_type), intent(in) :: A
    !!----    type (Space_Group_Type),     intent(in) :: SpG
    !!----    Type (Crystal_Cell_Type),    intent(in) :: Cell
    !!----    character(len=*),            intent(in) :: Filecod
    !!----    integer,                     intent(in) :: ndimx
    !!----    integer,                     intent(in) :: ndimy
    !!----    integer,                     intent(in) :: ndimz
    !!----    character(len=*),            intent(in) :: atname  !Species of the probe atom e.g. LI+1
    !!----    real(kind=cp),               intent(in) :: drmax
    !!----    real(kind=cp), optional,     intent(in) :: delta
    !!----    real(kind=cp), optional,     intent(out):: vol
    !!----    real(kind=cp), optional,     intent(out):: emin
    !!----    integer,       optional,     intent(out):: npix
    !!----    real(kind=cp), optional,  dimension(:,:,:), allocatable, intent(out) :: bvel_map
    !!----
    !!----    Calculate Bond-Valence energy landscape map where the energy at each point of the grid
    !!----    is determined by a species representative defined in atname. The BV-site Energy value
    !!----    is evaluated for distances below drmax value. If delta is present only the points with
    !!----    energy below EMin+delta .
    !!----
    !!---- Created: January 2015
    !!
    Subroutine Calc_Map_BVEL(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta,vol,emin,npix,outp,bvel_map)
       !---- Arguments ----!
       type (Atoms_Conf_List_type), intent(in) :: A
       type (Space_Group_Type),     intent(in) :: SpG
       Type (Crystal_Cell_Type),    intent(in) :: Cell
       character(len=*),            intent(in) :: Filecod
       integer,                     intent(in) :: ndimx
       integer,                     intent(in) :: ndimy
       integer,                     intent(in) :: ndimz
       character(len=*),            intent(in) :: atname  !Species of the probe atom e.g. LI+1
       real(kind=cp),               intent(in) :: drmax
       real(kind=cp), optional,     intent(in) :: delta
       real(kind=cp), optional,     intent(out):: vol
       real(kind=cp), optional,     intent(out):: emin
       integer,       optional,     intent(out):: npix
       logical,       optional,     intent(in) :: outp
       real(kind=cp), optional,  dimension(:,:,:), allocatable, intent(out) :: bvel_map

       !---- Local variables ----!
       character(len=4)                             :: car,atm
       integer                                      :: i,j,k,n,n1,n2,np,jbvs
       integer                                      :: nx1,nx2,ny1,ny2,nz1,nz2
       integer                                      :: i1,j1,k1,sig1,sig2,ncont
       real(kind=cp)                                :: rx1,ry1,rz1,qval,q1,q2,rep,p,s,cose
       real(kind=cp)                                :: sbvs, dd, occ, radius, rho, dmin, d_cutoff, &
                                                       dzero, alpha,c_rep,c_atr
       real(kind=cp), dimension(3)                  :: pto,pta,step,extend
       real(kind=cp),   dimension(:,:,:), allocatable :: map_bvs
       integer(kind=2), dimension(:,:,:), allocatable :: contrib
       type (Atom_list_Type)                        :: At1, At2
       logical                                      :: anion,cation
       !Coulomb constant (1/4pie0) to get energies in eV, use electron charges and distances in angstroms
       real(kind=dp), parameter :: ke=14.399644850445155254866066
       !Principal quantum numbers of the test=ion  and the species of all the atoms of the list
       real(kind=cp) :: n_tion, ferfc
       real(kind=cp), dimension(:), allocatable :: n_j
       logical :: all_present


       !---- Initial conditions ----!
       if (A%natoms <= 0) return
       if (ndimx <= 0 .or. ndimy <= 0 .or. ndimz <= 0) return

       all_present=present(delta) .and. present(vol) .and. present(npix) .and. present(emin)
       !---- Preparing the Variables and calculations ----!
       call Allocate_Atom_List(A%natoms,At1)
       At1%atom=A%atom  !Atom list A%Atom(:)%ind(1) contains the pointer to the species
                        !coming from A (Atoms_Conf_List_type, set in the main program)
       atm=u_case(atname)
       n1=0
       allocate(n_j(A%N_Spec))
       n_j=0.0
       if (.not. allocated(Ap_Table)) call Set_Atomic_Properties()
       do j=1,A%N_Spec
         car=u_case(A%species(j))
         do i=1,Ap_species_n
           if(Ap_Table(i)%Symb == car) then
              n_j(j)=real(Ap_Table(i)%n)
              exit
           end if
         end do
         if(car == atm) n_tion=n_j(j)
       end do

       do i=1,At1%natoms  !This exclude all atoms with the same species as Atname
          n2=A%Atom(i)%ind(1)
          if (trim(u_case(A%species(n2))) == trim(atm)) then
             At1%atom(i)%active=.false.
             q1=At1%atom(i)%charge
             sig1=nint(SIGN(1.0_cp,q1))
             n1=n2
             radius=A%radius(n1)
             car=A%Species(n1)
          end if
       end do
       if (n1 ==0) then
          err_conf=.true.
          ERR_Conf_Mess="The point atom "//atname//" is not in the Species Table"
          return
       end if

       call AtList1_ExtenCell_AtList2(Spg,At1,At2,.true.)
       !call Deallocate_atom_list(At1)
       !check that all species are well set in the list
       !Write(unit=*,fmt="(a)") " => List of atoms for calculating BVEL"
       do n=1,At2%natoms
           n2=At2%Atom(n)%ind(1)
           !write(unit=*,fmt="(a,a,3f12.5)") At2%Atom(n)%Lab,At2%Atom(n)%SfacSymb,At2%Atom(n)%x
           if (n2 ==0) then
              err_conf=.true.
              ERR_Conf_Mess="The atom "//trim(At2%Atom(n)%lab)//" is not in the Species Table"
              return
           end if
       end do

       allocate(map_bvs(ndimx,ndimy,ndimz),contrib(ndimx,ndimy,ndimz))
       map_bvs=0.0
       contrib=0
       ! Determination of the valence expected for a good position
       i=index(car,"+")
       cation=.false.
       anion=.false.
       if( i /= 0) then
         read(unit=car(i+1:),fmt=*) qval
         cation=.true.
       else
         i=index(car,"-")
         read(unit=car(i+1:),fmt=*) qval
         anion=.true.
       end if
       step=(/ 1.0/real(ndimx),  1.0/real(ndimy), 1.0/real(ndimz) /)
       extend=(/ drmax/cell%cell(1),  drmax/cell%cell(2), drmax/cell%cell(3) /)

       np=0; npix=0
       !---- Map Calculation ----!
       do k=1,ndimz
          pto(3)=(k-1)*step(3)
          do j=1,ndimy
             pto(2)=(j-1)*step(2)
             do i=1,ndimx
                pto(1)=(i-1)*step(1)

                rx1=pto(1)-extend(1)
                if (rx1 <= 0.0) then
                    nx1=int(rx1)-1
                else
                    nx1=int(rx1)
                end if
                nx2=int(pto(1)+extend(1))

                ry1=pto(2)-extend(2)
                if (ry1 <= 0.0) then
                    ny1=int(ry1)-1
                else
                    ny1=int(ry1)
                end if
                ny2=int(pto(2)+extend(2))

                rz1=pto(3)-extend(3)
                if (rz1 <= 0.0) then
                    nz1=int(rz1)-1
                else
                    nz1=int(rz1)
                end if
                nz2=int(pto(3)+extend(3))

                sbvs=0.0
                rep=0.0
                ncont=0
                do n=1,At2%natoms
                   n2=At2%Atom(n)%ind(1)
                   q2=At2%Atom(n)%charge
                   sig2= nint(SIGN(1.0_cp,q2))
                   rho=(radius+A%radius(n2))*0.74
                   dzero=Table_Dzero(n1,n2)
                    dmin=Table_Rmin(n1,n2)
                   alpha=Table_Alpha(n1,n2)
                   d_cutoff=Table_Rcutoff(n1,n2)
                   occ=At2%Atom(n)%VarF(1)
                   c_rep=occ*q1*q2/sqrt(n_tion*n_j(n2))
                   c_atr=occ*dzero
                   ferfc=erfc(drmax/rho)/drmax !always below 10^(-9) when drmax/rho > 4.2
                   do k1=nz1,nz2
                      do j1=ny1,ny2
                         do i1=nx1,nx2
                            pta=At2%Atom(n)%x+real((/i1,j1,k1/))
                            dd=max(Distance(pto,pta,Cell),0.0001) !To avoid division by zero
                            if (dd > drmax) cycle
                            contrib(i,j,k)=contrib(i,j,k)+1
                            if (sig1 == sig2) then
                                rep=rep + c_rep*(erfc(dd/rho)/dd-ferfc)
                            else
                               !if(dd > d_cutoff) cycle
                               sbvs=sbvs+ c_atr*((exp(alpha*(dmin-dd))-1.0)**2-1.0)
                            end if
                         end do
                      end do
                   end do

                end do
                !Multiply the repulsion term by the Coulomb constant to convert to eV
                map_bvs(i,j,k)=sbvs+ke*rep
             end do
          end do
       end do
       !Calculation of the volume available for mobility path
       if(all_present) then
          emin=minval(map_bvs)
          npix=count(map_bvs <= emin+delta)
          p=1.0; s=1.0
          do i=1,3
             cose=cosd(Cell%ang(i))
             p=p*cose
             s=s-cose*cose
          end do
          vol=sqrt(abs(s+2.0*p))
          do i=1,3
             p=Cell%cell(i)*step(i)
             vol=vol*p
          end do
          vol=vol*real(npix)
       end if

       !---- Export a File ----!
       if(present(outp)) then
          call Get_LogUnit(jbvs)
          open(unit=jbvs,file=trim(filecod)//".map",status="replace",action="write")

          write (unit=jbvs, fmt='(a)') "BVEL Map Calculations using Bond_STR Program"
          write (unit=jbvs, fmt='(a)') "BVEL Map for species "//trim(car)
          write (unit=jbvs, fmt='(a,3f12.4,3x,3f8.3)') "CELL ",cell%cell,cell%ang
          write (unit=jbvs, fmt='(a)') "SPGR  "//trim(SpG%spg_symb)
          write (unit=jbvs, fmt='(a,f10.4,a)')" => Global distance cutoff:",drmax," angstroms"
          write (unit=jbvs,fmt="(a,/,a,/)")  &
          " Bond-Valence Energy parameters (D0,Rmin,alpha) for Morse Potential:  D0*[{exp(alpha(dmin-d))-1}^2-1]", &
                "   (data read from internal table or provided by the user)"
          do n1=1,A%N_Cations
             do j=1,A%N_Anions
                n2=A%N_Cations+j
                write(unit=jbvs,fmt="(2(a,i3,a,a4),/,3(a,f9.5),/,3(a,f9.5),a)")           &
                      "   Type",n1,": ",A%Species(n1)," with type",n2,": ",A%Species(n2), &
                      "    D0  =",Table_Dzero(n1,n2),"       Rmin =",Table_Rmin(n1,n2),   &
                      "  Alpha =",Table_Alpha(n1,n2),"  Av. Coord.=",Table_Avcoor(n1,n2), &
                      "    R0  =",Table_Rzero(n1,n2),"   R-cutoff =",Table_Rcutoff(n1,n2),&
                      "   => Reference: "//trim(references(Table_ref(n1,n2)))
                write(unit=jbvs,fmt="(2(a,a,a,f6.3,a),/)") &
                      "   Cation (Eff. radius): ",A%Species(n1),"(",A%Radius(n1),")   ",  &
                      "   Anion  (Eff. radius): ",A%Species(n2),"(",A%Radius(n2),")"
             end do
          end do
          write(unit=jbvs,fmt="(a,i3,a)") &
           "   Principal Quantum Number (test ion): ",nint(n_tion)," for species: "//trim(atname)
          do j=1,A%N_Spec
                write(unit=jbvs,fmt="(a,i3,a)") &
                      "   Principal Quantum Number: ",nint(n_j(j))," for species: "//A%Species(j)
          end do
          write (unit=jbvs, fmt='(a)')     " List of atoms  "
          write (unit=jbvs, fmt='(a)')     "      Label  Species       X           Y           Z         Biso     OccFactor    Occupancy "
          do i=1,At1%natoms
            if(At1%atom(i)%active) then
                write(unit=jbvs, fmt='(a5,2a6,t20,6f12.5)')"Atom ",trim(A%Atom(i)%lab)," "//A%Atom(i)%SfacSymb, &
                                       A%Atom(i)%x, A%Atom(i)%biso, A%Atom(i)%occ, A%Atom(i)%VarF(1)
            end if
          end do
          write (unit=jbvs, fmt='(/,a)')     " Expanded List of atoms  "
          write (unit=jbvs, fmt='(a)')     "    #      Label  Species     X           Y           Z         Biso     OccFactor    Occupancy "
          do i=1,At2%natoms
                write(unit=jbvs, fmt='(i5,3a6,t23,6f12.5)') i," Atom ",trim(At2%Atom(i)%lab)," "//At2%Atom(i)%SfacSymb, &
                                       At2%Atom(i)%x, At2%Atom(i)%biso, At2%Atom(i)%occ, At2%Atom(i)%VarF(1)
          end do
          if(present(delta) .and. present(vol)) then
             write (unit=jbvs, fmt='(/,a,f10.4,a)')   "Value of delta for volumen calculation:",delta," eV"
             write (unit=jbvs, fmt='(a,f10.4,a)')     "Available volume for ion mobility in the unit cell:",vol," angstroms^3"
             write (unit=jbvs, fmt='(a,f10.4,a)')     "Volume  fraction for ion mobility in the unit cell:",vol/Cell%CellVol*100.0, "%"
             write (unit=jbvs, fmt='(a,f10.4,a,i8)')  "Minum Energy (in eV):", emin,"  Number of pixels with Emin < Energy < Emin+Delta: ",npix
          end if
          write (unit=jbvs, fmt='(a)')     "! Grid values: ndimx,ndimy,ndimz [0,ndimx-1] [0,ndimy-1] [0,ndimz-1]  "
          write (unit=jbvs, fmt='(a,9i6)') "GRID ",ndimx,ndimy,ndimz,0,ndimx-1,0,ndimy-1,0,ndimz-1
          if(.not. outp) then
              write (unit=jbvs, fmt='(a)')     "    X         Y         Z        Energy(eV)   #Contr.Terms"
              do k=1,ndimz
                 pto(3)=(k-1)*step(3)
                 do j=1,ndimy
                    pto(2)=(j-1)*step(2)
                    do i=1,ndimx
                       pto(1)=(i-1)*step(1)
                       ncont=contrib(i,j,k)
                       write(unit=jbvs, fmt='(3f10.5,f14.6,i12)') pto,map_bvs(i,j,k),ncont
                    end do
                 end do
              end do
          else
             write (unit=jbvs, fmt='(a)')     "DENSITY_MAP"
             write (unit=jbvs,fmt='(10g14.5)') map_bvs
          end if
         close(unit=jbvs)
       end if
       if(present(bvel_map)) then
         allocate(bvel_map(ndimx,ndimy,ndimz))
         bvel_map=map_bvs
       end if
       call write_grid_VESTA(map_bvs,cell,"Bond Valence Energy Landscape",trim(filecod)//"_bvel","P")
       return
    End Subroutine Calc_Map_BVEL

    !!----
    !!---- Subroutine Calc_Map_BVS(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta,vol)
    !!----    type (Atoms_Conf_List_type), intent(in) :: A
    !!----    type (Space_Group_Type),     intent(in) :: SpG
    !!----    Type (Crystal_Cell_Type),    intent(in) :: Cell
    !!----    character(len=*),            intent(in) :: Filecod
    !!----    integer,                     intent(in) :: ndimx
    !!----    integer,                     intent(in) :: ndimy
    !!----    integer,                     intent(in) :: ndimz
    !!----    character(len=*),            intent(in) :: atname
    !!----    real(kind=cp),               intent(in) :: drmax
    !!----    real(kind=cp), optional,     intent(in) :: delta !Tolerance in v.u. for output the map
    !!----    real(kind=cp), optional,     intent(out):: vol
    !!----
    !!----    Calculate a map of BVS values at each point of the grid is determined
    !!----    by a species representative defined in atname. The BVS value is evaluated
    !!----    for distances below drmax value. If delta is present only the points with valence
    !!----    close ( q-delta < BVS < q+delta) to that of the atname are output in the map.
    !!----
    !!---- Update: December - 2007, December 2014 ( JRC,change the order of indices and VESTA output)
    !!
    Subroutine Calc_Map_BVS(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta,vol)
       !---- Arguments ----!
       type (Atoms_Conf_List_type), intent(in) :: A
       type (Space_Group_Type),     intent(in) :: SpG
       Type (Crystal_Cell_Type),    intent(in) :: Cell
       character(len=*),            intent(in) :: Filecod
       integer,                     intent(in) :: ndimx
       integer,                     intent(in) :: ndimy
       integer,                     intent(in) :: ndimz
       character(len=*),            intent(in) :: atname  !Species of the probe atom e.g. LI+1
       real(kind=cp),               intent(in) :: drmax
       real(kind=cp), optional,     intent(in) :: delta
       real(kind=cp), optional,     intent(out):: vol

       !---- Local variables ----!
       character(len=4)                             :: car,atm
       integer                                      :: i,j,k,n,n1,n2,np,L,nL
       integer                                      :: nx1,nx2,ny1,ny2,nz1,nz2
       integer                                      :: i1,j1,k1,sig1,sig2
       integer                                      :: jbvs=63,npix
       integer,      dimension(:),     allocatable  :: ind
       real(kind=cp)                                :: rx1,ry1,rz1,qval,dif,q1,q2,rep,p,s
       real(kind=cp)                                :: sbvs, dd, occ, radius,sig, cose
       real(kind=cp), dimension(3)                  :: pto,pta,step,extend
       real(kind=cp), dimension(:,:,:), allocatable :: map_bvs
       real(kind=cp), dimension(:,:),   allocatable :: peaks
       real(kind=cp), dimension(:),     allocatable :: VD_peaks
       type (Atom_list_Type)                        :: At1, At2
       logical                                      :: new_peak,anion,cation

       !---- Initial conditions ----!
       if (A%natoms <= 0) return
       if (ndimx <= 0 .or. ndimy <= 0 .or. ndimz <= 0) return

       !---- Preparing the Variables and calculations ----!
       call Allocate_Atom_List(A%natoms,At1)
       At1%atom=A%atom  !Atom list A%Atom(:)%ind(1) contains the pointer to the species
                        !coming from A (Atoms_Conf_List_type, set in the main program)
       atm=u_case(atname)
       n1=0
       do i=1,At1%natoms  !This exclude all atoms with the same species as Atname
          n2=A%Atom(i)%ind(1)
          if (trim(u_case(A%species(n2))) == trim(atm)) then
             At1%atom(i)%active=.false.
             q1=At1%atom(i)%charge
             sig1=nint(SIGN(1.0_cp,q1))
             n1=n2
             radius=A%radius(n1)
             car=A%Species(n1)
          end if
       end do
       if (n1 ==0) then
          err_conf=.true.
          ERR_Conf_Mess="The point atom "//atname//" is not in the Species Table"
          return
       end if

       call AtList1_ExtenCell_AtList2(Spg,At1,At2,.true.)
       call Deallocate_atom_list(At1)
       !check that all species are well set in the list
       !write(unit=*,fmt="(a)") " => List of atoms for calculating BVS map"
       do n=1,At2%natoms
           n2=At2%Atom(n)%ind(1)
           !write(unit=*,fmt="(a,a,3f12.5)") At2%Atom(n)%Lab,At2%Atom(n)%SfacSymb,At2%Atom(n)%x
           if (n2 ==0) then
              err_conf=.true.
              ERR_Conf_Mess="The atom "//trim(At2%Atom(n)%lab)//" is not in the Species Table"
              return
           end if
       end do

       allocate(map_bvs(ndimx,ndimy,ndimz))
       map_bvs=0.0

       if(.not. present(delta)) then
          allocate(peaks(4,ndimx)) ! A maximum of ndimx peaks will be stored
          peaks=0.0
          allocate(VD_peaks(ndimx)) ! Vector holding the differences of bond-valences
          VD_peaks=0.0
          allocate(ind(ndimx))      ! Vector pointer for ordering the peaks
          ind=0
       end if

       ! Determination of the valence expected for a good position
       i=index(car,"+")
       cation=.false.
       anion=.false.
       if( i /= 0) then
         read(unit=car(i+1:),fmt=*) qval
         cation=.true.
       else
         i=index(car,"-")
         read(unit=car(i+1:),fmt=*) qval
         anion=.true.
       end if
       step=(/ 1.0/real(ndimx),  1.0/real(ndimy), 1.0/real(ndimz) /)
       extend=(/ drmax/cell%cell(1),  drmax/cell%cell(2), drmax/cell%cell(3) /)

       np=0; npix=0
       !---- Map Calculation ----!
       do k=1,ndimz
          pto(3)=(k-1)*step(3)
          do j=1,ndimy
             pto(2)=(j-1)*step(2)
             do i=1,ndimx
                pto(1)=(i-1)*step(1)

                rx1=pto(1)-extend(1)
                if (rx1 <= 0.0) then
                    nx1=int(rx1)-1
                else
                    nx1=int(rx1)
                end if
                nx2=int(pto(1)+extend(1))

                ry1=pto(2)-extend(2)
                if (ry1 <= 0.0) then
                    ny1=int(ry1)-1
                else
                    ny1=int(ry1)
                end if
                ny2=int(pto(2)+extend(2))

                rz1=pto(3)-extend(3)
                if (rz1 <= 0.0) then
                    nz1=int(rz1)-1
                else
                    nz1=int(rz1)
                end if
                nz2=int(pto(3)+extend(3))

                sbvs=0.0
                rep=0.0
                do n=1,At2%natoms
                   n2=At2%Atom(n)%ind(1)
                   q2=At2%Atom(n)%charge
                   sig2= nint(SIGN(1.0_cp,q2))
                   sig=(radius+A%radius(n2))*0.99
                   do k1=nz1,nz2
                      do j1=ny1,ny2
                         do i1=nx1,nx2
                            pta=At2%Atom(n)%x+real((/i1,j1,k1/))
                            occ=At2%Atom(n)%VarF(1)
                            dd=max(Distance(pto,pta,Cell),0.0001) !To avoid division by zero
                            if (dd > drmax) cycle
                            if (sig1 == sig2) then
                               rep=rep + (sig/dd)**18
                               cycle
                            end if
                            sbvs=sbvs+occ*exp((table_d0(n1,n2)-dd)/table_b(n1,n2))
                         end do
                      end do
                   end do
                end do
                dif=abs(sbvs-qval)
                if(present(delta)) then
                  if(dif > delta .or. rep > 0.01) cycle
                else
                  !Algorithm for selecting the optimum positions
                  if(dif < 0.15 .and. np < ndimx) then
                    new_peak=.true.
                    do L=1,np
                       dd=Distance(peaks(1:3,L),pto,Cell)
                       if( dd < 0.8 ) then
                         new_peak=.false.
                         nL=L
                         exit
                       end if
                    end do
                    if(new_peak) then
                      np=np+1
                      peaks(1:3,np)= pto
                      peaks(4,np)  = sbvs
                      VD_peaks(np) = dif
                    else  !now compare with the peak stored at nL and interchange them if sbvs is more favourable
                      if( dif < abs(qval-peaks(4,nL)) ) then
                        peaks(1:3,nL) = pto
                        peaks(  4,nL) = sbvs
                        VD_peaks(nL)  = dif
                      end if
                    end if
                  end if
                end if
                !end of peaks construction
                npix=npix+1
                map_bvs(i,j,k)=sbvs
             end do
          end do
       end do
       !Calculation of the volume available for mobility path
       if(present(delta) .and. present(vol))then
          p=1.0; s=1.0
          do i=1,3
             cose=cosd(Cell%ang(i))
             p=p*cose
             s=s-cose*cose
          end do
          vol=sqrt(abs(s+2.0*p))
          do i=1,3
             p=Cell%cell(i)*step(i)
             vol=vol*p
          end do
          vol=vol*real(npix)
       end if
       !Sorting the favourable peak positions for inserting the species of atom Atname
       if(.not. present(delta)) call sort(VD_peaks,np,ind)
       !---- Export a File ----!
       !call Get_LogUnit(jbvs)
       open(unit=jbvs,file=trim(filecod)//".map",status="replace",action="write")
       write (unit=jbvs, fmt='(a)') "BVS Map Calculations using Bond_STR Program"
       write (unit=jbvs, fmt='(a)') "BVS Map for species "//trim(car)
       write (unit=jbvs, fmt='(a,3f12.4,3x,3f8.3)') "CELL ",cell%cell,cell%ang
       write (unit=jbvs, fmt='(a)')     "SPGR  "//trim(SpG%spg_symb)
       write (unit=jbvs, fmt='(a,f10.4,a)')" => Global distance cutoff:",drmax," angstroms"
       write (unit=jbvs,fmt="(/,a,/,a,/)")  &
       " Bond-Valence parameters (d0,B0) for Zachariasen formula:  s= exp{(d0-d)/B0}", &
             "   (data read from internal table or provided by the user)"
       do n1=1,A%N_Cations
          do j=1,A%N_Anions
             n2=A%N_Cations+j
                write(unit=jbvs,fmt="(2(a,i3,a,a4),2(a,f6.3),a)") &
                      "   Type",n1,": ",A%Species(n1)," with type",n2,": ",A%Species(n2),&
                      " d0=",Table_d0(n1,n2),"    B0=",Table_b(n1,n2), "   => Reference: "//trim(references(Table_ref(n1,n2)))
                write(unit=jbvs,fmt="(2(a,a,a,f6.3,a),/)") &
                      "   Cation (Eff. radius): ",A%Species(n1),"(",A%Radius(n1),")   ", &
                      "   Anion  (Eff. radius): ",A%Species(n2),"(",A%Radius(n2),")"
          end do
       end do
       write (unit=jbvs, fmt='(a)')     "! List ot atoms  "
       do i=1,A%natoms
         write(unit=jbvs, fmt='(a,t20,5f12.5)')"Atom "//trim(A%Atom(i)%lab)//"  "//A%Atom(i)%SfacSymb, &
                                                A%Atom(i)%x, A%Atom(i)%biso, A%Atom(i)%occ
       end do
       if(present(delta) .and. present(vol))then
          write (unit=jbvs, fmt='(/,a,f10.4,a)')  "Value of delta for volumen calculation:",delta," eV"
          write (unit=jbvs, fmt='(a,f10.4,a)')    "Available volume for ion mobility in the unit cell:",vol," angstroms^3"
          write (unit=jbvs, fmt='(a,f10.4,a)')    "Volume  fraction for ion mobility in the unit cell:",vol/Cell%CellVol*100.0, "%"
       end if
       if(.not. present(delta)) then
         write (unit=jbvs, fmt='(a)')     "! List ot favourable positions for inserting the species "//trim(car)
         do i=1,np
           j=ind(i)
           write(unit=jbvs, fmt='(a,i4,a,3f10.5,a,f8.3)')"#",i,"  Position: (",peaks(1:3,j),")  Valence: ",peaks(4,j)
         end do
       end if
       write (unit=jbvs, fmt='(a)')     "! Grid values: ndimx,ndimy,ndimz [0,ndimx-1] [0,ndimy-1] [0,ndimz-1]  "
       write (unit=jbvs, fmt='(a,9i6)') "GRID ",ndimx,ndimy,ndimz,0,ndimx-1,0,ndimy-1,0,ndimz-1
       write (unit=jbvs, fmt='(a)')     "DENSITY_MAP"
       write (unit=jbvs,fmt='(10g14.5)') map_bvs
       close(unit=jbvs)
       call write_grid_VESTA(map_bvs,cell,"Bond Valence Map",trim(filecod)//"_bvs","P")

       !---- End Procedure ----!
       if (allocated(map_bvs)) deallocate(map_bvs)
       call deallocate_atom_list(At2)

       return
    End Subroutine Calc_Map_BVS

    !!----
    !!---- Subroutine Complete_Table(A,N_bvsm,bvs_m)
    !!----    type (Atoms_Conf_List_type),   intent(in) :: A      !  In -> Object of Atoms_Conf_List_type
    !!----    Integer,                       intent(in) :: N_bvsm
    !!----    Character(len=*),dimension(:), intent(in) :: bvs_m
    !!----
    !!----    Sets up the table of BV parameters (add provided external values)
    !!----    completing the table when the user gives his/her own BV parameters
    !!----
    !!---- Update: January - 2008
    !!
    Subroutine Complete_Table(A,N_bvsm,bvs_m)
       type (Atoms_Conf_List_type),   intent(in) :: A      !  In -> Object of Atoms_Conf_List_type
       Integer,                       intent(in) :: N_bvsm
       Character(len=*),dimension(:), intent(in) :: bvs_m
       !---- Local variables ----!
       Integer       :: i,j,k,ic,icat,ian
       real(kind=cp) :: d0_n, b_n
       character(len=10), dimension(5)        :: dire

          do k=1,N_bvsm
                 dire=" "
                 call getword(bvs_m(k),dire,ic)
                 if (ic <= 0) cycle
                 icat=0
                 ian=0
                 if (ic < 3 ) then
                          err_conf=.true.
                          ERR_Conf_Mess="Cation-Anion D0,b parameters must be provided"
                          return
                 end if

                 do i=1,A%N_Spec
                          if (u_case(dire(1)(1:4)) /= A%Species(i)) cycle
                          icat=i
                 end do
                 do i=1,A%N_Spec
                          if (u_case(dire(2)(1:4)) /= A%Species(i)) cycle
                          ian=i
                 end do
                 if (icat == 0 .or. ian == 0) then
                    err_conf=.true.
                    ERR_Conf_Mess="The given Cation or Anion cannot be found in the atom list"
                    return
                 end if
                 if (icat > ian) then
                          j=icat
                          icat=ian
                          ian=j
                 end if
                 if (icat > A%N_Cations) then
                    err_conf=.true.
                    ERR_Conf_Mess="One of the given cations is not found in the atom list"
                    return
                 end if

                 read(unit=dire(3),fmt=*) d0_n
                 if (ic > 3) read(unit=dire(4),fmt=*) b_n

                 Table_d0(icat,ian)=d0_n
                 Table_b(icat,ian)=b_n
                 Table_ref(icat,ian)=0

                 Table_d0(ian,icat)=d0_n
                 Table_b(ian,icat)=b_n
                 Table_ref(ian,icat)=0

          end do

    End Subroutine Complete_Table

    !!----
    !!---- Subroutine Complete_Table_BVEL(A,N_bvel,bvel)
    !!----    type (Atoms_Conf_List_type),   intent(in) :: A      !  In -> Object of Atoms_Conf_List_type
    !!----    Integer,                       intent(in) :: N_bvel
    !!----    Character(len=*),dimension(:), intent(in) :: bvel
    !!----
    !!----    Sets up the table of BVEL parameters (add provided external values)
    !!----    Completing the table when the user gives his/her own BVEL parameters
    !!----
    !!---- Created: December - 2014
    !!
    Subroutine Complete_Table_BVEL(A,N_bvel,bvel)
       type (Atoms_Conf_List_type),   intent(in) :: A      !  In -> Object of Atoms_Conf_List_type
       Integer,                       intent(in) :: N_bvel
       Character(len=*),dimension(:), intent(in) :: bvel
       !---- Local variables ----!
       Integer       :: i,j,k,ic,icat,ian
       real(kind=cp) :: Avcoor,Rzero,Rcutoff,Dzero,Rmin,alpha
       character(len=12), dimension(9)  :: dire

          do k=1,N_bvel
                 dire=" "
                 call getword(bvel(k),dire,ic)
                 if (ic <= 0) cycle
                 icat=0
                 ian=0
                 if (ic < 8 ) then
                          err_conf=.true.
                          ERR_Conf_Mess="Cation-Anion {Nc,R0,D0,Rmin,alpha} parameters must be provided"
                          return
                 end if

                 do i=1,A%N_Spec
                          if (u_case(dire(1)(1:4)) /= A%Species(i)) cycle
                          icat=i
                 end do
                 do i=1,A%N_Spec
                          if (u_case(dire(2)(1:4)) /= A%Species(i)) cycle
                          ian=i
                 end do
                 if (icat == 0 .or. ian == 0) then
                    err_conf=.true.
                    ERR_Conf_Mess="The given Cation or Anion cannot be found in the atom list"
                    return
                 end if
                 if (icat > ian) then
                          j=icat
                          icat=ian
                          ian=j
                 end if
                 if (icat > A%N_Cations) then
                    err_conf=.true.
                    ERR_Conf_Mess="One of the given cations is not found in the atom list"
                    return
                 end if

                 read(unit=dire(3),fmt=*) Avcoor
                 read(unit=dire(4),fmt=*) Rzero
                 read(unit=dire(5),fmt=*) Rcutoff
                 read(unit=dire(6),fmt=*) Dzero
                 read(unit=dire(7),fmt=*) Rmin
                 read(unit=dire(8),fmt=*) alpha

                 Table_Avcoor  (icat,ian)=Avcoor
                 Table_Rzero   (icat,ian)=Rzero
                 Table_Rcutoff (icat,ian)=Rcutoff
                 Table_Dzero   (icat,ian)=Dzero
                 Table_Rmin    (icat,ian)=Rmin
                 Table_alpha   (icat,ian)=alpha
                 Table_ref     (icat,ian)=0

                 Table_Avcoor  (ian,icat)=Avcoor
                 Table_Rzero   (ian,icat)=Rzero
                 Table_Rcutoff (ian,icat)=Rcutoff
                 Table_Dzero   (ian,icat)=Dzero
                 Table_Rmin    (ian,icat)=Rmin
                 Table_alpha   (ian,icat)=alpha
                 Table_ref     (ian,icat)=0

          end do

    End Subroutine Complete_Table_BVEL
    !!----
    !!---- Subroutine Cost_BVS(A, GII, ERep, gic)
    !!----    type (Atoms_Conf_List_type),  intent(in out) :: A    !  In  -> Object of Atoms_Conf_List_type
    !!----    real(kind=cp),                intent(out)    :: GII  !  OUT -> Global instability index
    !!----    real(kind=cp),      optional, intent(out)    :: ERep !  OUT -> Repulstion term from soft spheres
    !!----    character(len=*),   optional, intent(in)     :: gic  ! If present GII_c is put in GII
    !!----
    !!----    Subroutine to calculate the Global Instability Index.
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Set_TDist_Coordination" in order
    !!----    to update the internal private variables related to distance/angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    All items corresponding to the bond-valence parameters contained in A have to
    !!----    be properly set before calling this procedure.
    !!----
    !!---- Update: April - 2005, July -2104 (JRC, calculation only for the subset having VarF(5) = 0)
    !!
    Subroutine Cost_BVS(A, GII, ERep,gic)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),  intent(in out)  :: A    !  In -> Object of Atoms_Conf_List_type
       real(kind=cp),                intent(out)     :: GII  !  GII_a
       real(kind=cp),      optional, intent(out)     :: ERep !  Repulsion term due to shoft spheres
       character(len=*),   optional, intent(in)      :: gic  !  If present GII_c is put in GII

       !---- Local variables ----!
       integer       :: i,j,k,icm,l,sig1,sig2
       real(kind=cp) :: tol,fact,q2,dd,sums,q1, del, bv,gii_a,gii_c,efcn,sig,rep

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20
       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii_a=0.0
       gii_c=0.0
       rep=0.0
       do i=1,A%natoms
          if(A%Atom(i)%VarF(5) > 0.01) cycle
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0_cp,q1)
          sums=0.0
          efcn=0.0
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(j,i))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(j,i))%charge
             sig2= SIGN(1.0_cp,q2)
             sig=A%radius(l)+A%radius(k)
             dd=coord_info%dist(j,i)
             if (sig1 == sig2) then
                rep=rep + (0.8*sig/dd)**18
                cycle
             end if
             if (dd > sig*(1.0+tol)) cycle
             efcn=efcn+A%Atom(coord_info%n_cooatm(j,i))%VarF(1)
             bv=EXP((Table_d0(l,k)-dd)/Table_b(l,k))
             bv=bv*A%Atom(coord_info%n_cooatm(j,i))%VarF(1) !Occupacy
             sums=sums+bv
          end do
          A%Atom(i)%varF(3)=efcn
          del=sums-ABS(q1)

          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii_a=gii_a+abs(del)*fact
          gii_c=gii_c+del*del*fact

       end do
       gii_a=gii_a*100.0/A%totatoms
       gii_c=sqrt(gii_c/A%totatoms)*100.0
       GII=gii_a
       if(present(ERep)) ERep=rep
       if(present(gic)) GII=gii_c
       return
    End Subroutine Cost_BVS

    !!----
    !!---- Subroutine Cost_BVS_CoulombRep(A, GII, ERep)
    !!----    type (Atoms_Conf_List_type),  intent(in out):: A     !  In  -> Object of Atoms_Conf_List_type
    !!----    real(kind=cp),                intent(out)   :: GII   !  OUT -> Global instability index
    !!----    real(kind=cp),                intent(out)   :: ERep   !  Pseudo Repulsion Coulomb "energy"
    !!----
    !!----
    !!----    Subroutine to calculate the Global Instability Index Gii_a and
    !!----    a pseudo Coulomb repulsion energy useful to avoid cation-cation and
    !!----    anion-anion overlap when using this cost function for predicting or
    !!----    solving a ionic crystal structure. It was used in the old program PiXSA,
    !!----    by J. Pannetier, J. Bassas-Alsina, J.Rodriguez-Carvajal and V. Caignaert,
    !!----    in "Prediction of Crystal Structures from Crystal Chemistry Rules by
    !!----    Simulated Annealing", Nature 346, 343-345 (1990).
    !!----    An additional repulsive term (soft sphere) of the form Esph=(sig/d)**18,
    !!----    with sig equal to the sum of ionic radii of the pair, has been introduced.
    !!----
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Set_TDist_Coordination" in order
    !!----    to update the internal Coord_Info variables related to distance and
    !!----    angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    All items corresponding to the bond-valence parameters contained in
    !!----    "A" have to be properly set before calling this procedure.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Cost_BVS_CoulombRep(A, GII, ERep)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),  intent(in out):: A      !  In -> Object of Atoms_Conf_List_type
       real(kind=cp),                intent(out)   :: GII    !  GII_a
       real(kind=cp),                intent(out)   :: ERep   !  Pseudo Repulsion Coulomb "energy"

       !---- Local variables ----!
       integer        :: i,j,k,icm,l,sig1,sig2
       real(kind=cp)  :: tol,fact,q2,dd,sums,q1, del, bv,efcn,sig

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20

       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii=0.0
       Erep=0.0
       do i=1,A%natoms
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0_cp,q1)
          sums=0.0
          efcn=0.0
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(j,i))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(j,i))%charge
             sig2= SIGN(1.0_cp,q2)
             dd=coord_info%dist(j,i)
             sig=A%radius(l)+A%radius(k) !sum of ionic radii of the current pair
             if (sig1 == sig2) then
                Erep=Erep+ q1*q2/dd + (sig/dd)**18  !Coulomb potential + soft sphere repulsion (avoid short distances)
                cycle
             end if
             if (dd > sig*(1.0+tol)) cycle
             efcn=efcn+A%Atom(coord_info%n_cooatm(j,i))%VarF(1)
             bv=EXP((Table_d0(l,k)-dd)/Table_b(l,k))
             bv=bv*A%Atom(coord_info%n_cooatm(j,i))%VarF(1)
             sums=sums+bv
          end do

          del=sums-ABS(q1)
          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii=gii+abs(del)*fact
          A%Atom(i)%varF(3)=efcn
       end do
       gii=gii*100.0/A%totatoms

       return
    End Subroutine Cost_BVS_CoulombRep

    !!----
    !!---- Subroutine Deallocate_Atoms_Conf_List(A)
    !!----    type (Atoms_Conf_List_Type), intent(in out)   :: A  ! In/ Out -> Objet to be deallocated
    !!----
    !!----    De-allocation of objet A of type Atoms_Conf_List. This subroutine
    !!----    should be after using an object of type Atoms_Conf_List that is no
    !!----    more needed.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Deallocate_Atoms_Conf_List(A)
       !---- Arguments ----!
       type (Atoms_Conf_List_Type), intent(in out)   :: A  !Objet to be deallocated

       if (allocated(A%atom)) deallocate (A%atom)
       A%natoms   = 0
       A%n_spec   = 0
       A%n_anions = 0
       A%n_cations= 0

       return
    End Subroutine Deallocate_Atoms_Conf_List

    !!----
    !!---- Subroutine Ewald(Lattvec,Vol,Ac,e)
    !!----
    !!----   Real(kind=cp), Dimension(3,3),  Intent(in)   :: Lattvec !Transformation matrix from Crystallographic to Cartesian components
    !!----   Real(kind=cp),                  Intent(in)   :: Vol     !Volume of the unit cell
    !!----   Type (Atoms_Cell_Type),         Intent(in)   :: Ac      !List of atoms in the unit cell
    !!----   Real(kind=dp),                  Intent(out)  :: e       !Electrostatic lattice energy
    !!----
    !!----   This subroutine computes the electrostatic energy using the Ewald method.
    !!----
    !!----  Update: February - 2018
    !!
    Subroutine Ewald(Lattvec,Vol,Ac,e)
      Real(kind=cp), Dimension(3,3),  Intent(in)   :: Lattvec
      Real(kind=cp),                  Intent(in)   :: Vol
      Type (Atoms_Cell_Type),         Intent(in)   :: Ac
      Real(kind=8),                   Intent(out)  :: e

      Integer                                         i, j, k, h, l, ik, nk
      Integer, Dimension(3)                        :: ncell, nrcell

      Real(kind=8), Parameter                      :: k_coul = 14.39963977d0   ! eV * angstrom
      Real(kind=8), Parameter                      :: k_boltz = 8.617332358e-5 ! eV / K
      Real(kind=8), Parameter                      :: beta = 0.4
      Real(kind=8), Parameter                      :: rcut = 12.               ! angstrom
      Real(kind=8), Parameter                      :: kcut =  4.               ! angstrom ^ -1
      Real(kind=8)                                    rcut2, kcut2, four_beta2, fourpi
      Real(kind=8)                                 :: r2, e_s, e_l, e_self, fase, sk2
      Real(kind=8),  Dimension(2)                  :: sk
      Real(kind=8),  Dimension(3)                  :: latt, rlatt, rdmin, kdmin, r0, r, a, b, c
      Real(kind=cp), Dimension(3,3)                :: rlattvec
      Real(kind=cp), Dimension(:), Allocatable     :: expok
      Real(kind=cp), Dimension(:,:), Allocatable   :: xyz, errorfc, g

      ! Constants

      rcut2         = rcut ** 2
      kcut2         = kcut ** 2
      fourpi      = 4.0 * pi
      four_beta2  = 4.0 * beta ** 2

      ! Compute reciprocal lattice vectors rlattvec

      rlattvec(1:3,1) = cross_product(lattvec(1:3,2),lattvec(1:3,3))
      rlattvec(1:3,2) = cross_product(lattvec(1:3,3),lattvec(1:3,1))
      rlattvec(1:3,3) = cross_product(lattvec(1:3,1),lattvec(1:3,2))

      rlattvec = rlattvec * tpi / Vol

      ! Compute module of lattice and reciprocal lattice vectors

      Do i = 1 , 3
         latt(i)  = Dot_product(lattvec(1:3,i),lattvec(1:3,i))
         rlatt(i) = Dot_product(rlattvec(1:3,i),rlattvec(1:3,i))
         latt(i)  = Sqrt(latt(i))
         rlatt(i) = Sqrt(rlatt(i))
      End Do

      Do i = 1 , 3
         rdmin(i)  = Dot_product(lattvec(1:3,i),rlattvec(1:3,i))
         kdmin(i)  = rdmin(i) / latt(i)
         rdmin(i)  = rdmin(i) / rlatt(i)
         ncell(i)  = Int(rcut / rdmin(i)) + 1
         nrcell(i) = Int(kcut / kdmin(i)) + 1
      End Do

      ! Direct to Cartesian Coordinates

      Allocate(xyz(3,Ac%Nat))
      xyz=0.0
      Do i = 1 , Ac%Nat
        xyz(:,i)= matmul(lattvec,Ac%XYZ(:,i))
      End Do

      ! - Arrays construction: ERRORFC, G, EXPOK

      ! - ERROFC

      Allocate(errorfc(Ac%Nat,Ac%Nat))

      Do i = 1 , Ac%Nat
         Do j = 1 , Ac%Nat
            errorfc(i,j) = 0.
            r0 = xyz(:,j) - xyz(:,i)
            Do h = -ncell(1) , ncell(1)
               Do k = -ncell(2) , ncell(2)
                  Do l = -ncell(3) , ncell(3)
                     If (i /= j .Or. &
                         h /= 0 .Or. &
                         k /= 0 .Or. &
                         l /= 0) Then
                         r = r0 + h * lattvec(:,1) + k * lattvec(:,2) + l * lattvec(:,3)
                        r2 = Dot_Product(r,r)
                        If (r2 <= rcut2) Then
                           r(1)         = Sqrt(r2)
                           errorfc(i,j) = errorfc(i,j) + erfc(beta * r(1)) / r(1)
                        End If

                     End If
                  End Do
               End Do
            End Do
         End Do
      End Do

      ! - G:  reciprocal vectors used in the sum

      nk = (2 * nrcell(1) + 1) * (2 * nrcell(2) + 1) * (2 * nrcell(3) + 1)

      Allocate(g(1:nk,1:4))

      nk = 0
      i  = 0
      g  = 0

      Do h = - nrcell(1) , nrcell(1)
         a = h * rlattvec(:,1)
         Do k = - nrcell(2) , nrcell(2)
            b = k * rlattvec(:,2)
            Do l = - nrcell(3) , nrcell(3)
               c = l * rlattvec(:,3)
               r = a  + b + c
               r2 = Dot_Product(r,r)
               If (r2 <= kcut2 .And. r2 > 0.0d0) Then
                  nk        = nk + 1
                  g(nk,1:3) = r
                  g(nk,4)   = r2
               End If
            End Do
         End Do
      End Do

      ! - EXPOK

      Allocate(expok(1:nk))

      Do ik = 1 , nk
         expok(ik) = Exp((-g(ik,4) / four_beta2)) / g(ik,4)
      end Do

      ! - Compute energy

      e_self        = 0.0d0
      e_s           = 0.0d0
      e_l           = 0.0d0

      ! - E_SELF

      Do i = 1 , Ac%Nat
         e_self = e_self + Ac%Charge(i) * Ac%Charge(i)
      End Do

      ! - E_S (Real space contribution)

      Do i = 1 , Ac%Nat
         Do j = 1 , Ac%Nat
            e_s = e_s + Ac%Charge(i) * Ac%Charge(j) * errorfc(i,j)
         End Do
      End Do

      ! E_L (Reciprocal space contribution)

      Do ik = 1 , nk
         sk(1:2) = 0.0d0
         Do j = 1 , Ac%Nat
            fase = Dot_product(xyz(1:3,j),g(ik,1:3))
            sk(1) = sk(1) + Ac%Charge(j) * cos(fase)
            sk(2) = sk(2) + Ac%Charge(j) * Sin(fase)
         End Do
         sk2 = sk(1) ** 2 + sk(2) ** 2
         e_l = e_l + expok(ik) * sk2
      End Do
      e_s    = e_s * k_coul / 2.
      e_l    = e_l * k_coul * fourpi / 2. / Vol
      e_self = e_self * k_coul * beta  / Sqrt(pi)
      e      = (e_s + e_l - e_self) / Ac%Nat
    End Subroutine Ewald

    !!----
    !!---- Subroutine get_soft_covalent_radius(nam,radius)
    !!----    character(len=*), intent(in) :: nam    ! In -> Name of the species containing the charge
    !!----    real(kind=cp),    intent(out):: radius ! Out-> Covalent radius according to Ap_Table
    !!----
    !!----    Getting the covalent radius according to Ap_Table
    !!----
    !!---- Created: January - 2015
    !!
    Subroutine get_soft_covalent_radius(nam,radius)
      character(len=*), intent(in) :: nam
      real(kind=cp),    intent(out):: radius
      !--- Local variables ---!
      integer :: i
      radius=0.0
      do i=1,Ap_species_n
        if(nam == Ap_Table(i)%Symb) then
          radius=Ap_table(i)%Rc
          exit
        end if
      end do
      if(radius < 0.01) radius=1.0
    End Subroutine get_soft_covalent_radius

    !!----
    !!---- Subroutine Init_Err_Conf()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Init_Err_Conf()

       err_conf=.false.
       ERR_Conf_Mess=" "

       return
    End Subroutine Init_Err_Conf

    !!----
    !!---- Subroutine Set_Charges(SpGr,Cell,A,filcod)
    !!----     Type (Space_Group_Type),     Intent(in)    :: SpGr
    !!----     Type (Crystal_Cell_Type),    Intent(in)    :: Cell
    !!----     Type (Atom_List_Type),       Intent(inout) :: A
    !!----     Real(kind=cp), Dimension(*), Intent(out)   :: charges
    !!----     Character(len=*), Optional,  Intent(in)    :: filcod
    !!----
    !!----    This subroutine sets atomic charges for object A using
    !!----    bond-valence sums. Output is stored in Formal_Charges.fch
    !!----    or filcod.fch if filcod is provided.
    !!----
    !!---- Update: July - 2016
    !!
    Subroutine Set_Formal_Charges(SpGr,Cell,A,eps_val,iwrt)
      !---- Arguments ----!
      Type (Space_Group_Type),     Intent(in)    :: SpGr
      Type (Crystal_Cell_Type),    Intent(in)    :: Cell
      Type (Atom_List_Type),       Intent(in out):: A
      Real(kind=cp),    Optional,  Intent(in)    :: eps_val
      Integer,          Optional,  Intent(in)    :: iwrt
      !---- Local variables ----!
      Type (Atoms_Cell_Type)                     :: Ac,Bc

      Character(len=2)                           :: label

      Integer                                    :: i, j, k, l, n, m, icomb, p, diff, penalty_aux, orbit_mult,lun
      Integer                                    :: nelements, nanions_candidate, ncombinations, nsolutions, nunusual_aux
      Integer                                    :: nsol_usual, nelemmix, nequiv_sol, npenalty_min, penalty_min
      Integer, Parameter                         :: massive_limit = 100
      Integer, Dimension(:), Allocatable         :: anion_candidate, first_neighbor, iontype, nvalences, contador, nunusual
      Integer, Dimension(:), Allocatable         :: atom2elem, elemmix
      Integer, Dimension(:), Allocatable         :: penalty, Ac2Bc
      Integer, Dimension(:,:), Allocatable       :: valence

      Real(kind=cp)                              :: rmin, total_charge, energy_aux, norm, epsv
      Real(kind=cp), Parameter                   :: ZEROCHARGE = 1.0e-4
      Real(kind=cp)                              :: Dmax=4.0, Dangle=0.0
      Real(kind=dp),  Dimension(:), Allocatable  :: energy
      Real(kind=cp), Dimension(:), Allocatable   :: sumocc, solution_aux
      Real(kind=cp), Dimension(:,:), Allocatable :: solution, combination, solution_new

      Logical                                    :: error_pairs=.False.,error_element=.False.
      Logical                                    :: warning_valence
      Logical                                    :: newelement, incremento, usual, ordered, xyz_stored, newsolution, identical

      Type Element
         Character(len=2)                        :: ChemSymb
         Character(len=10)                       :: Lab
         Integer                                 :: Z
         Integer                                 :: Type
      End Type Element

      Type(Element), Dimension(:), Allocatable   :: Elem

      ERR_Char  = .False.
      WARN_Char = .False.
      epsv=0.01
      lun=6
      if(present(iwrt)) lun=iwrt
      if(present(eps_val)) epsv=eps_val

      Write(unit=lun,fmt="(/,/,7(a,/))")                           &
           "             ============================"           , &
           "             ====== FORMAL CHARGES ======"           , &
           "             ============================"           , &
           "    ***********************************************" , &
           "    *  Formal charges from  *.cfl or *.cif files  *" , &
           "    ***********************************************" , &
           "      (Nebil A. Katcho - ILL, version: July 2016 )"
      Write(unit=lun,fmt="(a,/)") " "

      ! -------------------
      ! Build the unit cell
      ! -------------------
      call Deallocate_Atoms_Cell(Ac)
      Call Allocate_Atoms_Cell(A%Natoms,SpGr%Multip,0.0,Ac)
      Call Set_Epsg(epsv) !needed for well controlling the calculation of multiplicities
      Write(unit=lun,fmt="(a,i3)")  "    Space Group Number: ", SpGr%NumSpg
      Write(unit=lun,fmt="(a,a,/)") "    Space Group Symbol: ", SpGr%Spg_Symb

      Ac%Nat = 0
      j      = 1

      Do i = 1 , A%Natoms
         Ac%Nat = Ac%Nat + A%Atom(i)%Mult
         Ac%Noms(j:Ac%Nat) = A%Atom(i)%Lab
         Call Get_Orbit(A%Atom(i)%X,SpGr,orbit_mult,Ac%XYZ(:,j:Ac%Nat))
         If (orbit_mult /= A%Atom(i)%Mult) Then
            ERR_Char = .True.
            ERR_Char_Mess = "Error in multiplicities. Check space group (or increase epsg)"
            Return
         End If
         j = j + A%Atom(i)%Mult
      End Do

      ! We build Bc. It is build from Ac, we only avoid to repeat lattice sites. Bc and Ac
      ! are different when two different atomic species occupy the same lattice sites. Bc
      ! is needed to compute the electrostatic energy.
      Allocate(Bc%XYZ(3,Ac%Nat),Bc%Charge(Ac%Nat))
      Allocate(Ac2Bc(Ac%Nat)) ! maps atoms of Ac into Bc
      Allocate(sumocc(Ac%Nat))

      sumocc(:) = 0.
      k         = 0
      Do i = 1 , A%Natoms
         Do m = 1 , A%Atom(i)%Mult
            xyz_stored = .False.
            j = 1
            k = k + 1
            Do   !While (.Not. xyz_stored .And. j <= Bc%Nat)
               If (xyz_stored .Or. j > Bc%Nat) Exit
               If (Abs(Bc%XYZ(1,j)-Ac%XYZ(1,k)) < 0.01 .And. &
                   Abs(Bc%XYZ(2,j)-Ac%XYZ(2,k)) < 0.01 .And. &
                   Abs(Bc%XYZ(3,j)-Ac%XYZ(3,k)) < 0.01) Then
                  xyz_stored = .True.
               Else
                  j = j + 1
               End If
            End Do
            If (.Not. xyz_stored) Then
               Bc%Nat           = Bc%Nat + 1
               Bc%XYZ(:,Bc%Nat) = Ac%XYZ(:,k)
            End If
            Ac2Bc(k)       = Bc%Nat
            sumocc(Bc%Nat) = sumocc(Bc%Nat) + A%Atom(i)%Occ
         End Do

      End Do

      ! We compute the number of elements and the number of elements that can be considered
      ! as anions

      Call Set_Chem_Info()
      Call Set_Oxidation_States_Table()
      Call Set_Common_Oxidation_States_Table()

      nelements         = 0
      nanions_candidate = 0

      Allocate(anion_candidate(A%natoms),atom2elem(A%natoms))

      Do i = 1 , A%natoms
         newelement         = .True.
         anion_candidate(i) = 0
         Label = A%Atom(i)%ChemSymb
         Call Get_ChemSymb(label,A%Atom(i)%ChemSymb,A%Atom(i)%Z)
         Do j = 1 , i - 1
            If (A%Atom(i)%ChemSymb == A%Atom(j)%ChemSymb) Then
               newelement = .False.
               atom2elem(i) = atom2elem(j)
            End If
         End Do

         j = 1
         Do
            If (j > 11) Exit
            If (anion_candidate(i) /= 0) Exit
            If (OxStates_Table(j,A%Atom(i)%Z) < 0) Then
               anion_candidate(i) = -1
               nanions_candidate  = nanions_candidate + 1
            End If
            j = j + 1
         End Do

         If (anion_candidate(i) == -1) nanions_candidate = nanions_candidate + 1
         If (newelement) Then
            nelements    = nelements + 1
            atom2elem(i) = nelements
         End If
      End Do

      Write(unit=lun,fmt="(a,i2/)") " => Number of elements found: ", nelements
      Write(unit=lun,fmt="(a,/)") " => Searching ions that can be anions...."
      Write(unit=lun,fmt="(tr4,a)") "----------------"
      Write(unit=lun,fmt="(tr4,a)") "Anion candidates"
      Write(unit=lun,fmt="(tr4,a)") "----------------"

      Do i = 1 , A%natoms
         If (anion_candidate(i) == 0) Then
            Write(unit=lun,fmt="(tr4,a4,tr9,a)") A%Atom(i)%Lab, " NO"
         Else
            Write(unit=lun,fmt="(tr4,a4,tr9,a)") A%Atom(i)%Lab, "YES"
         End If
      End Do

      Write(unit=lun,fmt="(tr4,a,/)") "----------------"

      If (nelements == 1) Then
         ! We have an element, not a compound
         ! Formal charges = 0.
         Do i = 1 , A%nAtoms
            A%Atom(i)%Charge = 0.0
         End Do
         Write(unit=lun,fmt="(/,a)") " => This is an element. All charges are have been set to zero"

      Else If (nanions_candidate == 0) Then
         ! No anions
         ! Formal charges = 0.0
         Do i = 1 , A%nAtoms
            A%Atom(i)%Charge = 0.0_cp
         End Do
         Write(unit=lun,fmt="(/,a)") " => No anions candidates have been found"
         Write(unit=lun,fmt="(tr4,a)") "All charges have been set to zero"

      Else
         ! Here we determine which ion is a cation, which is an anion
         !     1. we determine the first neighbor of each atom
         !     2. each pair must be a cation-anion pair

         Call Calc_Dist_Angle(Dmax,Dangle,Cell,SpGr,A)
         Allocate(first_neighbor(A%natoms))
         first_neighbor = 0

         ! First neighbors

         Write(unit=lun,fmt="(a,/)") " => Computing first neighbor distance, needed to define cations and anions..."
         Write(unit=lun,fmt="(tr4,a)") "-----------------------"
         Write(unit=lun,fmt="(tr4,a)")  "First neighbor distance"
         Write(unit=lun,fmt="(tr4,a)")  "-----------------------"

         Do i = 1 , A%natoms
            rmin = 1.0e8_cp
            Do j = 1 , coord_info%Max_coor
               If (coord_info%N_CooAtm(j,i) /= 0) Then
                  If (coord_info%Dist(j,i) < rmin) Then
                     rmin = coord_info%Dist(j,i)
                     first_neighbor(i) = coord_info%N_CooAtm(j,i)
                  End If
               End If
            End Do
            If (first_neighbor(i) == 0) Then
               ERR_Char = .True.
               ERR_Char_Mess = " => Error! Coordination for atom: "//A%Atom(i)%ChemSymb//" is zero. Increase Dmax!"
               Return
            End If
            Write(unit=lun,fmt="(tr4,a4,tr2,a4,tr7,f6.3)") A%Atom(i)%Lab,A%Atom(first_neighbor(i))%Lab,rmin
         End Do
         Write(unit=lun,fmt="(tr4,a,/)") "-----------------------"

         ! Setting anion/cation nature

         Call Set_Pauling_Electronegativity()

         Allocate(iontype(A%natoms))

         iontype(:) = 0

         ! Ions with anion_candidate = 0 be cations
         ! Ions with anion_candidate greater than zero can be cations or anions
         !     - if the first neighbor anion_candidate = 0, it is an anion
         !     - if the first neighbor anion_candidate /= 0,
         !          it is a cation if it has lower electronegativity
         !          it is an anion otherwise

         Do i = 1 , A%natoms
            If (anion_candidate(i) == 0) Then
               iontype(i) = 1
            Else
               If (anion_candidate(first_neighbor(i)) == 0) Then
                  iontype(i) = -1
               Else
                  If (PaulingX(A%Atom(i)%Z) < PaulingX(A%Atom(first_neighbor(i))%Z)) Then
                     iontype(i) = 1
                  Else
                     iontype(i) = -1
                  End If
               End If
            End If
         End Do

         Write(unit=lun,fmt="(a,/)") " => Anion / Cation assignment...."
         Write(unit=lun,fmt="(tr4,a)") "----------------"
         Write(unit=lun,fmt="(tr4,a)") "Anions / Cations"
         Write(unit=lun,fmt="(tr4,a)") "----------------"

         Do i = 1 , A%natoms
            If (iontype(i) == 1) Then
               Write(unit=lun,fmt="(tr4,a4,tr6,a)") A%Atom(i)%Lab, "CATION"
            Else If (iontype(i) == -1) Then
               Write(unit=lun,fmt="(tr4,a4,tr6,a)") A%Atom(i)%Lab, " ANION"
            End If
         End Do

         Write(unit=lun,fmt="(tr4,a,/)") "----------------"

         ! We check every first-neighbor pair is a cation-anion pair

         Do i = 1 , A%natoms
            If (iontype(i) == iontype(first_neighbor(i))) Then
               error_pairs=.True.
               If (iontype(i) == 1) Then
                  Write(unit=lun,fmt="(a)") " => Error! The pair "//Trim(A%Atom(i)%Lab)//"-"&
                       //Trim(A%Atom(first_neighbor(i))%Lab)//" is a CATION-CATION pair"
               Else
                  Write(unit=lun,fmt="(a)") " => Error! The pair "//Trim(A%Atom(i)%Lab)//"-"&
                       //Trim(A%Atom(first_neighbor(i))%Lab)//" is an ANION-ANION pair"
               End If
            End If
         End Do

         If (error_pairs) Then
            ERR_Char = .True.
            ERR_Char_Mess = " => Error in first-neighbor:  Cation-Cation or Anion-Anion first-neighbor pairs have been found"
            Return
         End If

         ! We check all atoms of a given element are all cations or all anions

         Allocate(elem(nelements))

         Do i = 1 , nelements
            j = 1
            Do
               If (atom2elem(j) == i) Then
                  elem(i)%Z        = A%Atom(j)%Z
                  elem(i)%Lab      = A%Atom(j)%Lab
                  elem(i)%type     = iontype(j)
                  elem(i)%ChemSymb = A%Atom(j)%ChemSymb
                  Exit
               Else
                  j = j + 1
               End If
            End Do
            k = j + 1
            Do l = k , A%Natoms
               If (atom2elem(l) == i .And. iontype(l) /= iontype(j)) Then
                  Write(unit=lun,fmt="(2(a,i4),a)") " => Error! Atom ", j, " and Atom ", l, " must be both cations or anions"
                  error_element = .True.
               End If
            End Do
         End Do

         If (error_element) Then
            ERR_Char = .True.
            ERR_Char_Mess = " => Error in element:  Cations and anions found for an element"
            Return
         End If

         ! We extract the possible valences for each element

         Allocate(nvalences(nelements))
         Allocate(valence(1:11,nelements))
         nvalences = 0

         Do i = 1 , nelements
            j = 1
            Do
               If (j > 11) Exit
               If (OxStates_Table(j,elem(i)%Z) == 0) Exit

               If (elem(i)%type == 1 .And. OxStates_Table(j,elem(i)%Z) > 0) Then
                  nvalences(i)            = nvalences(i) + 1
                  valence(nvalences(i),i) = OxStates_Table(j,elem(i)%Z)
               Else If (elem(i)%type == -1 .And. OxStates_Table(j,elem(i)%Z) < 0) Then
                  nvalences(i)            = nvalences(i) + 1
                  valence(nvalences(i),i) = OxStates_Table(j,elem(i)%Z)
               End If
               j = j + 1
            End Do
            If (nvalences(i) == 0) Then ! He, some noble gas elements
               nvalences(i) = 1
               valence(1,i) = 0
            End If
         End Do

         Write(unit=lun,fmt="(a,/)") " => Extracting valences from Oxidation States Table...."
         Write(unit=lun,fmt="(tr4,a)") "----------------------------------------------"
         Write(unit=lun,fmt="(tr4,a)") "Ion, Usual / (Reported) valences"
         Write(unit=lun,fmt="(tr4,a)") "----------------------------------------------"

         Do i = 1 , nelements
            If (elem(i)%Z == 8) nvalences(i) = 1
            Write(unit=lun,fmt="(tr4,a4)",advance='no') elem(i)%ChemSymb
            Do j = 1 , nvalences(i)
               k = 1
               usual = .False.
               Do !While (.Not. usual .And. k <= 8)
                  If (valence(j,i) == Common_OxStates_Table(k,elem(i)%Z)) usual = .True.
                  k = k + 1
                  If (usual .Or. k > 8) Exit
               End Do
               If (usual) Then
                  Write(unit=lun,fmt="(tr3,i2,tr1)",advance='no') valence(j,i)
               Else
                  Write(unit=lun,fmt="(tr2,a1,i2,a1)",advance='no') "(", valence(j,i), ")"
               End If
            End Do
            Write(unit=lun,fmt=*)
         End Do

         Write(unit=lun,fmt="(tr4,a,/)") "----------------------------------------------"

         Write(unit=lun,fmt="(a,/)") " => Searching for combinations that fulfill the electroneutrality rule..."

         ncombinations=1

         Do i = 1 , nelements
            ncombinations = ncombinations * nvalences(i)
         End Do

         Allocate(solution(nelements,ncombinations),combination(nelements,ncombinations))
         Allocate(contador(nelements))

         nsolutions  = 0
         contador(:) = 1
         icomb       = 1

         Do icomb = 1 , ncombinations
            total_charge = 0.
            Do j = 1 , A%natoms
               combination(atom2elem(j),icomb) = valence(contador(atom2elem(j)),atom2elem(j))
               total_charge = total_charge + valence(contador(atom2elem(j)),atom2elem(j)) * A%atom(j)%occ
            End Do

            If (Abs(total_charge) < ZEROCHARGE) Then
               nsolutions = nsolutions + 1
               Do j = 1 , A%Natoms
                  solution(atom2elem(j),nsolutions) = valence(contador(atom2elem(j)),atom2elem(j))
               End Do
            End If
            incremento = .False.
            j = nelements

            Do !While (.Not. incremento .And. icomb < ncombinations)
               If (incremento .Or. icomb == ncombinations) Exit
               If (contador(j) < nvalences(j)) Then
                  contador(j) = contador(j) + 1
                  incremento  = .True.
               Else
                  contador(j) = 1
                  j = j - 1
               End If
            End Do
         End Do

         If (nsolutions == 0) Then

            Write(unit=lun,fmt="(a)") "    No solution found with integer valences "
            Write(unit=lun,fmt="(a)") "    Trying solutions with non-integer valences (mixed valence)"
            Write(unit=lun,fmt="(a)")

            ! We try with mixed valences (non-integer valences) for transition and earth-rare metals, if any

            Allocate(elemmix(nelements))

            nelemmix   = 0
            elemmix(:) = 0

            Do i = 1 , nelements
               If ( (elem(i)%Z > 20 .And. elem(i)%Z < 31) .Or. &
                    (elem(i)%Z > 38 .And. elem(i)%Z < 49) .Or. &
                    (elem(i)%Z > 56 .And. elem(i)%Z < 81) .Or. &
                    (elem(i)%Z > 88)) Then
                  nelemmix   = nelemmix + 1
                  elemmix(i) = nelemmix
               End If
            End Do

            If (nelemmix == 0) Then
               ERR_Char = .True.
               ERR_Char_Mess = " => No transition / earth-rare ion found. Unable to set charges"
               Return
            Else
               norm = 0.
               nsolutions = ncombinations * nelemmix

               If (nelemmix > 1) Then
                  Deallocate(solution)
                  Allocate(solution(nelements,nsolutions))
               End If

               nsolutions = 0

               Do n = 1 , nelemmix
                  Do i = 1 , A%Natoms
                     If (elemmix(atom2elem(i)) == n ) norm = norm + A%Atom(i)%occ
                  End Do

                  Do i = 1 , ncombinations
                     total_charge = 0.
                     Do j = 1 , A%Natoms
                        If (elemmix(atom2elem(j)) /= n) &
                             total_charge  = total_charge + combination(atom2elem(j),i) * A%Atom(j)%occ
                     End Do
                     Do j = 1 , nelements
                        If (elemmix(j) /= n) Then
                           solution(j,nsolutions+1) = combination(j,i)
                        Else
                           solution(j,nsolutions+1) = - total_charge / norm
                        End If
                     End Do

                     ! check if this solution is new

                     newsolution = .True.
                     m  = 1

                     Do !While (newsolution .And. m <= nsolutions)
                        If (newsolution .Or. m > nsolutions) Exit
                        identical = .True.

                        Do k = 1 , nelements
                           If (solution(k,nsolutions+1) /= solution(k,m)) identical = .False.
                        End Do

                        If (identical) newsolution = .False.
                        m = m + 1
                     End Do

                     If (newsolution) nsolutions = nsolutions + 1
                  End Do
               End Do
            End If
         End If

         If (nsolutions == 0) Then
            ERR_Char = .True.
            ERR_Char_Mess = " => No solution found. Unable to set charges"
            Return
          End If

         ! We compute how many unusual valences are in each solution, and set a penalty according to this

         Allocate(penalty(nsolutions))
         Allocate(nunusual(nsolutions))

         nunusual(:)    = 0
         penalty(:)     = 0.

         Do i = 1 , nsolutions
            Do j = 1 , nelements
               k = 1
               usual = .False.
               Do !While (.Not. usual .And. k <= 8)
                  If (solution(j,i) == Common_OxStates_Table(k,elem(j)%Z)) usual = .True.
                  k = k + 1
                  If (usual .Or. k > 8) Exit
               End Do
               If (.Not. usual) Then
                  nunusual(i)    = nunusual(i) + 1
                  ! Compute penalization
                  p = 100
                  k = 1
                  Do
                     If (k > 8) Exit
                     If (Common_OxStates_Table(k,elem(j)%Z) == 0) Exit
                     If (elem(j)%type == 1 .And. Common_OxStates_Table(k,elem(j)%Z) > 0) Then
                        diff = Abs(solution(j,i) - Common_OxStates_Table(k,elem(j)%Z))
                        If (p > diff) p = diff
                     Else If (elem(j)%type == -1 .And. Common_OxStates_Table(k,elem(j)%Z) < 0) Then
                        diff = Abs(solution(j,i) - Common_OxStates_Table(k,elem(j)%Z))
                        If (p > diff) p = diff
                     End If
                     k = k + 1
                  End Do
                  penalty(i) = penalty(i) + p
               End If
            End Do
         End Do

         penalty_min = Minval(penalty)
         nsol_usual = 0.

         Do i = 1 , nsolutions
            If (nunusual(i) == 0) nsol_usual = nsol_usual + 1
         End Do

         Write(unit=lun,fmt="(tr4,i10,a)") nsolutions, " solutions found"
         Write(unit=lun,fmt="(tr4,i10,a,/)") nsol_usual, " solutions with usual formal charges (penalty = 0)"

         If (nsolutions > massive_limit) Then
            Write(unit=lun,fmt="(tr4,a)") "number of solutions too large"
            Write(unit=lun,fmt="(tr4,a,/)") "work only with solutions with the minimum penalty"

            npenalty_min = 0
            Allocate(solution_aux(nelements))

            Do i = 1 , nsolutions

               If (penalty(i) == penalty_min) Then
                  npenalty_min = npenalty_min + 1
                  solution_aux(:) = solution(:,i)
                  solution(:,npenalty_min) = solution_aux(:)
               End If

            End Do

            Deallocate(solution_aux)
            Allocate(solution_new(nelements,npenalty_min))
            solution_new(:,:) = solution(:,1:npenalty_min)
            Deallocate(solution)
            Allocate(solution(nelements,npenalty_min))
            solution(:,:) = solution_new(:,:)
            Deallocate(solution_new)
            Deallocate(penalty)
            Allocate(penalty(npenalty_min))
            penalty(:) = penalty_min
            nsolutions = npenalty_min
         End If

         ! We sort solutions according to their penalty

         Allocate(solution_aux(nelements))
         ordered = .False.

         Do !While (.Not. ordered)
            ordered = .True.
            Do i = 1 , nsolutions - 1
               If (penalty(i+1) < penalty(i)) Then
                  ordered = .False.
                  penalty_aux     = penalty(i+1)
                  nunusual_aux    = nunusual(i+1)
                  solution_aux(:) = solution(:,i+1)
                  penalty(i+1)    = penalty(i)
                  nunusual(i+1)   = nunusual(i)
                  solution(:,i+1) = solution(:,i)
                  penalty(i)      = penalty_aux
                  nunusual(i)     = nunusual_aux
                  solution(:,i)   = solution_aux(:)
               End If
            End Do
            If (ordered) Exit
         End Do

         ! Compute Energies

         Write(unit=lun,fmt="(a,/)") " => Computing energies..."

         Allocate(energy(nsolutions))
         energy=0.0
         Do i = 1 , nsolutions
            ! Assign charges to every atom and compute the energy
            k = 0
            Bc%Charge(:) = 0.
            Do j = 1 , A%Natoms
               Do m = 1 , A%Atom(j)%Mult
                  k = k + 1
                  Bc%Charge(Ac2Bc(k)) =  Bc%Charge(Ac2Bc(k)) + solution(atom2elem(j),i) * A%Atom(j)%Occ / sumocc(Ac2Bc(k))
               End Do
            End Do
            Call Ewald(Cell%CR_Orth_Cel,Cell%CellVol,Bc,energy(i))
         End Do
         ! For each penalty value, we sort solutions according to their energy
         m = 1
         Do !While (m < nsolutions)
            If (m >= nsolutions) Exit
            n = m
            j = m + 1
            Do       !While (j <= nsolutions)
               If (j > nsolutions) Exit
               If (penalty(m) == penalty(j)) n = n + 1
               j = j + 1
            End Do
            If (n /= m) Then
               ordered = .False.
               Do     !While (.Not. ordered)
                  ordered = .True.
                  Do i = m , n - 1
                     If (energy(i+1) < energy(i)) Then
                        ordered         = .False.
                        energy_aux      = energy(i+1)
                        nunusual_aux    = nunusual(i+1)
                        solution_aux(:) = solution(:,i+1)
                        energy(i+1)     = energy(i)
                        nunusual(i+1)   = nunusual(i)
                        solution(:,i+1) = solution(:,i)
                        energy(i)       = energy_aux
                        nunusual(i)     = nunusual_aux
                        solution(:,i)   = solution_aux(:)
                     End If
                  End Do
                  If (ordered) Exit
               End Do
               m = n
            Else
               m = m + 1
            End If
         End Do

         Write(unit=lun,fmt="(tr5)",advance="no")
         Do i = 1 , nelements
            Write(unit=lun,fmt="(a)",advance='no') "---------"
         End Do
         Write(unit=lun,fmt="(a)") "---------------------"
         Do i = 1 , nelements
            Write(unit=lun,fmt="(tr5,a4)",advance="no") elem(i)%ChemSymb
         End Do
         Write(unit=lun,fmt="(a)") "  Energy (eV/atom) Penalty"
         Write(unit=lun,fmt="(tr5)",advance="no")
         Do i = 1 , nelements
            Write(unit=lun,fmt="(a)",advance='no') "---------"
         End Do
         Write(unit=lun,fmt="(a)") "---------------------"

         Do i = 1 , nsolutions
            Write(unit=lun,fmt="(tr4)",advance="no")
            Do j = 1 , nelements
               Write(unit=lun,fmt="(f5.2,tr4)",advance='no') solution(j,i)
            End Do
            Write(unit=lun,fmt="(tr5,f8.3,i9)") energy(i), penalty(i)
         End Do

         Write(unit=lun,fmt="(tr5)",advance="no")

         Do i = 1 , nelements
            Write(unit=lun,fmt="(a)",advance='no') "---------"
         End Do
         Write(unit=lun,fmt="(a)") "---------------------"
         Write(unit=lun,fmt="()")

         ! Determine the number of solutions with the minimum penalty and the same energy
         !         minimum = penalty(1)

         nequiv_sol = 1
         Do i = 2 , nsolutions
            If (penalty(i) == penalty(1) .And. Abs(energy(i) - energy(1)) < 0.001) nequiv_sol = nequiv_sol + 1
         End Do

         If (nequiv_sol > 1) Then
            Write(unit=lun,fmt="(tr5,i2,a)") nequiv_sol, " equivalent solutions found."
            Write(unit=lun,fmt="(tr8,a,/)") "final solution = average of equivalent solutions"
         End If

         warning_valence = .False.
         Do i = 1 , A%Natoms
            A%Atom(i)%Charge = 0.
            Do j = 1 , nequiv_sol
               A%Atom(i)%Charge = A%Atom(i)%Charge + solution(atom2elem(i),j)
            End Do
            A%Atom(i)%Charge = A%Atom(i)%Charge / nequiv_sol
            If (Abs(Abs(A%Atom(i)%Charge) - Nint(Abs(A%Atom(i)%Charge))) > 0.001) warning_valence = .True.
         End Do

      End If

      If (warning_valence) Then
         Write(unit=lun,fmt="(a,/)") " => WARNING: Non-integer valences found"
      End If

    End Subroutine Set_Formal_Charges

    !!----
    !!---- Subroutine Set_Table_BVEL_Params(A,N_bvel,bvel)
    !!---- type (Atoms_Conf_List_type),            intent(in)  :: A
    !!---- integer,                      optional, intent(in)  :: N_bvel   !Number of bvel strings with externally provided values
    !!---- character(len=*),dimension(:),optional, intent(in)  :: bvel     !bvs strings with externally provided values
    !!---- logical,                      optional, intent(in)  :: soft     !Calculate D0 and Rmin
    !!----
    !!----
    !!---- Created: December - 2014
    !!
    Subroutine Set_Table_BVEL_Params(A,N_bvel,bvel,soft)
       type (Atoms_Conf_List_type),            intent(in)  :: A
       integer,                      optional, intent(in)  :: N_bvel
       character(len=*),dimension(:),optional, intent(in)  :: bvel
       logical,                      optional, intent(in)  :: soft

       !---- Local Variables ----!
       integer :: i,j,k,L,ia,ic,ac,an
       real(kind=cp) :: rmin,d0,cn,b0,r0
       real(kind=cp),parameter :: f1=0.9185, f2=0.2285 !in eV units

       if (A%N_Spec == 0) then
          err_conf=.true.
          ERR_Conf_Mess=" The number of different species is zero, tables cannot be set"
          return
       end if

       if (allocated(Table_Avcoor))   deallocate(Table_Avcoor)
       if (allocated(Table_Rzero))    deallocate(Table_Rzero)
       if (allocated(Table_Rcutoff))  deallocate(Table_Rcutoff)
       if (allocated(Table_Dzero))    deallocate(Table_Dzero)
       if (allocated(Table_Rmin))     deallocate(Table_Rmin)
       if (allocated(Table_Alpha))    deallocate(Table_Alpha)
       if (allocated(Table_Ref))      deallocate(Table_Ref)

       allocate(Table_Avcoor(A%N_Spec,A%N_Spec))  ; Table_Avcoor  =0.0
       allocate(Table_Rzero(A%N_Spec,A%N_Spec))   ; Table_Rzero   =0.0
       allocate(Table_Rcutoff(A%N_Spec,A%N_Spec)) ; Table_Rcutoff =0.0
       allocate(Table_Dzero(A%N_Spec,A%N_Spec))   ; Table_Dzero   =0.0
       allocate(Table_Rmin(A%N_Spec,A%N_Spec))    ; Table_Rmin    =0.0
       allocate(Table_Alpha(A%N_Spec,A%N_Spec))   ; Table_Alpha   =0.0
       allocate(Table_ref(A%N_Spec,A%N_Spec))     ; Table_ref     =0

       if(present(soft)) then  !Calculate Rmin and D0 from softBVS parameters and softness
          call Set_Atomic_Properties()
          call Set_SBVS_Table()
          do i=1,A%N_Cations
             ic=0; ac=0
             do j=1,sbvs_species_n
                if (trim(A%Species(i)) == trim(sBVS_Table(j)%Symb)) then
                   ic=j
                   do k=1,Ap_species_n
                     if (trim(A%Species(i)) == trim(Ap_Table(k)%Symb)) then
                       ac=k
                       exit
                     end if
                   end do
                   exit
                end if
             end do
             if (ic == 0 .or. ac == 0) then
                if(.not. present(N_bvel)) then
                  err_conf=.true.
                  ERR_Conf_Mess=" Cation not found on the internal list: "//A%Species(i)
                  return
                else
                   call Complete_Table_BVEL(A,N_bvel,bvel)
                   if(err_conf) then
                        return
                   else
                        cycle
                   end if
                end if
             end if

             do k=1,A%N_Anions
                ia=0; an=0
                do j=1,bvs_anions_n
                   if (trim(A%Species(A%N_Cations+k)) == trim(bvs_anions(j))) then
                      ia=j
                      do L=1,Ap_species_n
                        if (A%Species(i) == Ap_Table(L)%Symb) then
                          an=L
                          exit
                        end if
                      end do
                      exit
                   end if
                end do
                if (ia == 0 .or. an == 0) then
                   err_conf=.true.
                   ERR_Conf_Mess=" Anion not found on the internal list: "//A%Species(A%N_Cations+k)
                   return
                end if

                b0=sbvs_table(ic)%b_par(ia)
                r0=sbvs_table(ic)%d0(ia)
                cn=sbvs_table(ic)%cn(ia)
                !write(*,"(a,3f12.4)") " => b0,r0,cn: ", b0,r0,cn
                if( b0+r0 < 0.00001) then
                   err_conf=.true.
                   ERR_Conf_Mess=" Bond-Valence parameters not found for pair: "//trim(A%Species(i))//"-"//bvs_anions(ia)
                   return
                end if
                if( cn < 0.1 ) cn=6.0
                Table_Avcoor (i,A%N_Cations+k)=cn
                Table_Rzero  (i,A%N_Cations+k)=r0
                Table_Rcutoff(i,A%N_Cations+k)=sbvs_table(ic)%ctoff(ia)
                Table_Alpha  (i,A%N_Cations+k)=1.0/b0
                Table_ref    (i,A%N_Cations+k)=sbvs_table(ic)%refnum(ia)

                !Formula 24 of ref. 34
                Rmin=r0*(f1+f2*abs(Ap_Table(ac)%sigma-Ap_Table(an)%sigma))-b0*log(real(Ap_Table(ac)%oxs)/cn)
                Select Case(Ap_Table(ac)%b)
                   Case(0,1)
                      D0=14.4*real(Ap_Table(ac)%oxs*Ap_Table(an)%oxs)*b0*b0/(2.0*rmin*sqrt(real(Ap_Table(an)%n*Ap_Table(ac)%n)))
                   Case(2:)
                      D0=14.4*sqrt(real(Ap_Table(ac)%oxs*Ap_Table(an)%oxs))*b0*b0/(rmin*sqrt(real(Ap_Table(an)%n*Ap_Table(ac)%n)))
                End Select
                !write(*,"(a,2f12.4)") " => D0 & Rmin: ",D0,Rmin
                Table_Dzero  (i,A%N_Cations+k)=D0
                Table_Rmin   (i,A%N_Cations+k)=Rmin

                Table_Avcoor (A%N_Cations+k,i)=cn
                Table_Rzero  (A%N_Cations+k,i)=r0
                Table_Rcutoff(A%N_Cations+k,i)=sbvs_table(ic)%ctoff(ia)
                Table_Dzero  (A%N_Cations+k,i)=D0
                Table_Rmin   (A%N_Cations+k,i)=Rmin
                Table_Alpha  (A%N_Cations+k,i)=1.0/b0
                Table_ref    (A%N_Cations+k,i)=sbvs_table(ic)%refnum(ia)

             end do  !Anions
          end do  !Cations
          call deallocate_SBVS_Table()
          call deallocate_Ap_table()

       else  !Just use the tabulated values

          call Set_BVEL_Table()

          do i=1,A%N_Cations
             ic=0
             do j=1,bvel_species_n
                if (trim(A%Species(i)) == trim(BVEL_Table(j)%Symb)) then
                   ic=j
                   exit
                end if
             end do
             if (ic == 0) then
                if(.not. present(N_bvel)) then
                  err_conf=.true.
                  ERR_Conf_Mess=" Cation not found on the internal list: "//A%Species(i)
                  return
                else
                   call Complete_Table_BVEL(A,N_bvel,bvel)
                   if(err_conf) then
                        return
                   else
                        cycle
                   end if
                end if
             end if

             do k=1,A%N_Anions
                ia=0
                do j=1,bvel_anions_n
                   if (A%Species(A%N_Cations+k) == bvel_anions(j)) then
                      ia=j
                      exit
                   end if
                end do
                if (ia == 0) then
                   err_conf=.true.
                   ERR_Conf_Mess=" Anion not found on the internal list: "//A%Species(A%N_Cations+k)
                   return
                end if
                Table_Avcoor (i,A%N_Cations+k)=bvel_table(ic)%Avcoor(ia)
                Table_Rzero  (i,A%N_Cations+k)=bvel_table(ic)%Rzero(ia)
                Table_Rcutoff(i,A%N_Cations+k)=bvel_table(ic)%Rcutoff(ia)
                Table_Dzero  (i,A%N_Cations+k)=bvel_table(ic)%Dzero(ia)
                Table_Rmin   (i,A%N_Cations+k)=bvel_table(ic)%Rmin(ia)
                Table_Alpha  (i,A%N_Cations+k)=bvel_table(ic)%Alpha(ia)
                Table_ref    (i,A%N_Cations+k)=bvel_table(ic)%refnum(ia)

                Table_Avcoor (A%N_Cations+k,i)=bvel_table(ic)%Avcoor(ia)
                Table_Rzero  (A%N_Cations+k,i)=bvel_table(ic)%Rzero(ia)
                Table_Rcutoff(A%N_Cations+k,i)=bvel_table(ic)%Rcutoff(ia)
                Table_Dzero  (A%N_Cations+k,i)=bvel_table(ic)%Dzero(ia)
                Table_Rmin   (A%N_Cations+k,i)=bvel_table(ic)%Rmin(ia)
                Table_Alpha  (A%N_Cations+k,i)=bvel_table(ic)%Alpha(ia)
                Table_ref    (A%N_Cations+k,i)=bvel_table(ic)%refnum(ia)

             end do
          end do
          call Deallocate_BVEL_Table()
       end if

       return
    End Subroutine Set_Table_BVEL_Params

    !!----
    !!---- Subroutine Set_Table_d0_b(A,N_bvsm,bvs_m)
    !!---- type (Atoms_Conf_List_type),            intent(in)  :: A
    !!---- integer,                      optional, intent(in)  :: N_bvsm   !Number of bvs strings with externally provided values
    !!---- character(len=*),dimension(:),optional, intent(in)  :: bvs_m    !bvs strings with externally provided values
    !!----
    !!----
    !!----
    !!---- Updated: March - 2005, December-2014
    !!
    Subroutine Set_Table_d0_b(A,N_bvsm,bvs_m,soft)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),            intent(in)  :: A
       integer,                      optional, intent(in)  :: N_bvsm
       character(len=*),dimension(:),optional, intent(in)  :: bvs_m
       logical,                      optional, intent(in)  :: soft
       !---- Local Variables ----!
       integer :: i,j,k,ia,ic
       logical :: found

       if (A%N_Spec == 0) then
          err_conf=.true.
          ERR_Conf_Mess=" The number of different species is zero, tables cannot be set"
          return
       end if

       if (allocated(Table_d0))   deallocate(Table_d0)
       if (allocated(Table_b))    deallocate(Table_b)
       if (allocated(Table_ref))  deallocate(Table_ref)

       allocate(Table_d0(A%N_Spec,A%N_Spec))
       allocate(Table_b(A%N_Spec,A%N_Spec))
       allocate(Table_ref(A%N_Spec,A%N_Spec))

       Table_d0=0.0
       Table_b=0.37
       Table_ref = 0

       call Set_BVS_Table()
       if(present(soft)) call set_SBVS_Table()

       do i=1,A%N_Cations
          ic=0; found=.false.
          if(present(soft)) then
             do j=1,sbvs_species_n
                if (trim(A%Species(i)) == trim(sBVS_Table(j)%Symb)) then
                   ic=j
                   found=.true.
                   exit
                end if
             end do
          end if
          if(ic == 0) then
            do j=1,bvs_species_n
               if (trim(A%Species(i)) == trim(BVS_Table(j)%Symb)) then
                  ic=j
                  exit
               end if
            end do
          end if
          if (ic == 0) then
             if(.not. present(N_bvsm)) then
               err_conf=.true.
               ERR_Conf_Mess=" Cation not found on the internal list: "//A%Species(i)
               return
             else
                call Complete_Table(A,N_bvsm,bvs_m)
                if(err_conf) then
                     return
                else
                     cycle
                end if
             end if
          end if

          do k=1,A%N_Anions
             ia=0
             do j=1,bvs_anions_n
                if (trim(A%Species(A%N_Cations+k)) == trim(bvs_anions(j))) then
                   ia=j
                   exit
                end if
             end do
             if (ia == 0) then
                err_conf=.true.
                ERR_Conf_Mess=" Anion not found on the internal list: "//A%Species(A%N_Cations+k)
                return
             end if
             if(found) then
                Table_d0 (i,A%N_Cations+k)=sbvs_table(ic)%d0(ia)
                Table_b  (i,A%N_Cations+k)=sbvs_table(ic)%b_par(ia)
                Table_ref(i,A%N_Cations+k)=sbvs_table(ic)%refnum(ia)

                Table_d0 (A%N_Cations+k,i)=sbvs_table(ic)%d0(ia)
                Table_b  (A%N_Cations+k,i)=sbvs_table(ic)%b_par(ia)
                Table_ref(A%N_Cations+k,i)=sbvs_table(ic)%refnum(ia)
             else
                Table_d0 (i,A%N_Cations+k)=bvs_table(ic)%d0(ia)
                Table_b  (i,A%N_Cations+k)=bvs_table(ic)%b_par(ia)
                Table_ref(i,A%N_Cations+k)=bvs_table(ic)%refnum(ia)

                Table_d0 (A%N_Cations+k,i)=bvs_table(ic)%d0(ia)
                Table_b  (A%N_Cations+k,i)=bvs_table(ic)%b_par(ia)
                Table_ref(A%N_Cations+k,i)=bvs_table(ic)%refnum(ia)
             end if

          end do
       end do

       call Deallocate_BVS_Table()
       if(present(soft)) call Deallocate_BVS_Table()

       return
    End Subroutine Set_Table_d0_b


    !!----
    !!---- Subroutine Species_on_List(A,MulG,tol, covalent,softbvs)
    !!----    type (Atoms_Conf_List_Type), intent(in out) :: A
    !!----    Integer, optional,           intent(in)     :: MulG
    !!----    real(kind=cp), optional,     intent(in)     :: tol
    !!----    logical,       optional,     intent(in)     :: covalent
    !!----    logical,       optional,     intent(in)     :: softbvs
    !!----
    !!----    Determines the different species in the List and,
    !!----    optionally, sets the tolerance factor for ionic radii
    !!----    conditions and provides "corrected" occupation factors
    !!----    (mult/MulG) when the user is using a multiplier. The
    !!----    general multiplicity of the space group MulG must be
    !!----    provided in such a case. This first free variable of the
    !!----    Atom-type A%Atom%VFree(1) is set to the corrected
    !!----    occupation. The first atom in the list must completely
    !!----    occupy its site. If covalent is provided the covalent
    !!----    radius is used instead of the ionic radius. If softbvs is
    !!----    present the Morse parameters are calculated from softBVS
    !!----    parameters
    !!----
    !!---- Update: March - 2005, December 2014, January 2015
    !!
    Subroutine Species_on_List(A,MulG, tol, covalent,softbvs)
       !---- Arguments ----!
       type (Atoms_Conf_List_Type), intent(in out) :: A
       Integer, optional,           intent(in)     :: MulG
       real(kind=cp), optional,     intent(in)     :: tol
       logical,       optional,     intent(in)     :: covalent
       logical,       optional,     intent(in)     :: softbvs

       !---- Local variables ----!
       character(len=4), dimension(50) :: cation,anion,spec
       character(len=2)                :: car,cv
       character(len=4)                :: canio
       integer                         :: i,im,j,v,ns,nc,na
       real(kind=cp)                   :: fac1,fact


       if (A%natoms == 0) return

       ns=0
       spec  = " "
       nc=0
       cation= " "
       na=0
       anion = " "

       if(present(tol)) A%tol=tol

       if(present(MulG)) then

         fac1=A%atom(1)%Occ*real(MulG)/real(A%atom(1)%mult)
         fac1=1.0/fac1
         A%totatoms=0.0
         do i=1,a%natoms
            fact=real(MulG)/real(a%atom(i)%mult)
            A%Atom(i)%VarF(1)=A%atom(i)%occ*fact*fac1      !Site Occupancy (=1, full occupation)
            A%Atom(i)%VarF(2)=A%atom(i)%occ_std*fact*fac1  !standard deviation of Site Occupancy
            A%totatoms=A%totatoms + A%Atom(i)%VarF(1)*real(a%atom(i)%mult) !total number of atoms/conventional unit cell
         end do

       else

         A%totatoms=0.0
         do i=1,a%natoms
            A%Atom(i)%VarF(1)=A%atom(i)%occ  !The user has given site occupancy
            A%Atom(i)%VarF(2)=A%atom(i)%occ_std
            A%totatoms=A%totatoms + A%Atom(i)%VarF(1)*real(a%atom(i)%mult) !total number of atoms/conventional unit cell
         end do

       end if

       if(present(softbvs) .and. .not. allocated(Ap_Table)) call Set_Atomic_Properties()

       loop1:do i=1, A%natoms
          car=u_case(a%atom(i)%ChemSymb)
          v=nint(a%atom(i)%charge)
          if (v == 0) then
             err_conf=.true.
             ERR_Conf_Mess=" The Atom "//a%atom(i)%lab//"has not charge"
             return
          end if
          write(unit=cv,fmt="(i2)") v
          if (v > 0) cv(1:1)="+"
          canio=car//cv
          canio=pack_string(canio)

          if (v > 0) then
             do j=1,nc
                if (canio == cation(j)) cycle loop1
             end do
             nc=nc+1
             cation(nc)=canio
          else
             do j=1,na
                if (canio == anion(j)) cycle loop1
             end do
             na=na+1
             anion(na)=canio
          end if
          if (na+nc == 50) exit
       end do loop1

       ns=nc+na
       A%N_Spec    = ns
       A%N_Anions  = na
       A%N_Cations = nc

       !---- Order the Species vector ----!
       call sort_strings(cation(1:nc))
       call sort_strings(anion(1:na))
       spec(1:nc)=cation(1:nc)
       spec(nc+1:nc+na)=anion(1:na)

       if (allocated(A%Species)) deallocate(A%Species)
       allocate(A%Species(ns))
       A%Species=spec(1:ns)

       if (allocated(A%Radius)) deallocate(A%Radius)
       allocate(A%Radius(ns))

       do i=1,nc
          im=index(A%Species(i),"+")
          car=A%Species(i)(im+1:im+2)
          read(unit=car,fmt="(i1)") j
          car=A%Species(i)(1:im-1)
          if(present(covalent) .and. present(softbvs)) then
             call get_soft_covalent_radius(A%Species(i),A%Radius(i))
          else if(present(covalent)) then
             call get_covalent_radius(car,A%Radius(i))
          else
             call get_ionic_radius(car,j,A%Radius(i))
          end if
          if (A%Radius(i) < 0.01) A%Radius(i)=0.8
       end do

       do i=1,A%N_Anions
          do j=1,bvs_anions_n
             if (A%Species(nc+i) == bvs_anions(j)) then
                if(present(covalent) .and. present(softbvs)) then
                    call get_soft_covalent_radius(A%Species(nc+i),A%Radius(nc+i))
                else if(present(covalent)) then
                    call get_covalent_radius(A%Species(nc+i),A%Radius(nc+i))
                else
                    A%Radius(nc+i) = bvs_anions_rion(j)
                end if
                exit
             end if
          end do
       end do

       !---- Fix the index on Atom_type pointing to the Species vector ----!
       do i=1, A%natoms
          do j=1,ns
             im=index(A%Species(j),"+")
             if (im == 0) im=index(A%Species(j),"-")
             car=A%Species(j)(1:im-1)
             cv=A%Species(j)(im:im+1)
             if (cv(1:1)=="+") cv(1:1)=" "
                read(unit=cv,fmt="(i2)") v
                if (u_case(A%Atom(i)%ChemSymb) == car .and. nint(A%Atom(i)%charge) == v) then
                A%atom(i)%ind(1)=j
                exit
             end if
          end do
       end do

       return
    End Subroutine Species_on_List

 End Module CFML_BVS_Energy_Calc
