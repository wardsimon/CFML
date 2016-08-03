!!----
!!---- MODULE: CFML_Formal_Charges
!!----   INFO: Contains the subroutine to set atomic charges using Bond-Valence Theory
!!----
!!---- HISTORY
!!----    Update: 28/07/2016
!!----
!!---- DEPENDENCIES
!!----   Use CFML_GlobalDeps,                  Only: Cp                   
!!----   Use CFML_Math_General,                Only: Set_Epsg
!!----   Use CFML_Crystal_Metrics,             Only: Crystal_Cell_Type
!!----   Use CFML_Crystallographic_Symmetry,   Only: Space_Group_Type, Get_Orbit
!!----   Use CFML_Atom_TypeDef,                Only: Atom_List_Type, Atoms_Cell_Type, Allocate_Atoms_Cell
!!----   Use CFML_Scattering_Chemical_Tables,  Only: Get_ChemSymb, Set_Chem_Info
!!----   Use CFML_Geometry_Calc
!!----   Use CFML_BVSpar
!!----   Use CFML_BVS_Energy_Calc
!!----   Use CFML_Atomic_Data
!!----
!!---- VARIABLES
!!----    ERR_CHAR
!!----    ERR_CHAR_MESS
!!----    WARN_CHAR
!!----    WARN_CHAR_MESS
!!----    
!!---- PROCEDURES
!!----    Subroutines:
!!----       SET_CHARGES
!!----
!!
 Module CFML_Formal_Charges

   !---- Use Modules ----!
   Use CFML_GlobalDeps,                  Only: Cp                   
   Use CFML_Math_General,                Only: Set_Epsg
   Use CFML_Crystal_Metrics,             Only: Crystal_Cell_Type
   Use CFML_Crystallographic_Symmetry,   Only: Space_Group_Type,Get_Orbit
   Use CFML_Atom_TypeDef,                Only: Atom_List_Type,Atoms_Cell_Type,Allocate_Atoms_Cell
   Use CFML_Scattering_Chemical_Tables,  Only: Get_ChemSymb,Set_Chem_Info
   Use CFML_Geometry_Calc
   Use CFML_BVSpar
   Use CFML_BVS_Energy_Calc
   Use CFML_Atomic_Data

   !---- Definitions ----!

   !!----
   !!---- ERR_CHAR
   !!----    logical, public  :: err_char
   !!----
   !!----    Logical Variable taking the value .true. if an error in the module
   !!----
   !!---- Update: July - 2016
   !!
   Logical, Public :: ERR_Char

   !!----
   !!---- ERR_CHAR_MESS
   !!----    character(len=150), public:: ERR_Char_Mess
   !!----
   !!----    String containing information about the last error
   !!----
   !!---- Update: July - 2016
   !!
   Character(len=150), Public :: ERR_Char_Mess

   !!----
   !!---- WARN_CHAR
   !!----    logical, public  :: warn_char
   !!----
   !!----    Logical Variable taking the value .true. if a warning in the module
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
   
   !---------------------!
   !---- Subroutines ----!
   !---------------------!

    !!----
    !!---- Subroutine Set_Charges(SpGr,Cell,A,charges,filcod)
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
   Subroutine Set_Charges(SpGr,Cell,A,charges,filcod)
     !---- Arguments ----!
     Type (Space_Group_Type),     Intent(in)    :: SpGr
     Type (Crystal_Cell_Type),    Intent(in)    :: Cell
     Type (Atom_List_Type),       Intent(inout) :: A
     Real(kind=cp), Dimension(*), Intent(out)   :: charges
     Character(len=*), Optional,  Intent(in)    :: filcod

     !---- Local variables ----!
     Type (Atoms_Cell_Type)                     :: Ac
     Type (Atoms_Conf_List_Type)                :: Acl
     
     Character(len=2)                           :: label

     Integer                                    :: i,j,k,orbit_mult,icomb,l,is
     Integer                                    :: nmax_val=10
     Integer                                    :: nelements,nsites,nanions,nanions_candidate,ncombinations,nsolutions
     Integer                                    :: maxv,minv
     Integer                                    :: icm,icn,sig1,sig2
     Integer                                    :: lun=1
     Integer, Dimension(:), Allocatable         :: anion_candidate,nvalences_all,nvalences,first_neighbor,iontype,contador
     Integer, Dimension(:), Allocatable         :: atom2elem,atom2site,tmre
     Integer, Dimension(:,:), Allocatable       :: valence,valence_all
     
     Real(kind=cp)                              :: tol,s0,q2,dd,q1,str,dav
     Real(kind=cp)                              :: rmin,total_charge,sumq_min
     Real(kind=cp)                              :: penalty_aux,sumq_aux
     Real(kind=cp)                              :: Dmax=4.,Dangle=0.,zerocharge=0.1,ttol=30.,DQ=0.1
     Real(kind=cp), Dimension(:), Allocatable   :: siteocc,sumq,sumq_old,penalty
     Real(kind=cp), Dimension(:), Allocatable   :: bvsum_aux,solution_aux
     Real(kind=cp), Dimension(:,:), Allocatable :: combination,solution,bvsum
     
     Logical                                    :: soft
     Logical                                    :: newelement,newsite,incremento,electroneutrality,ordered
     
     Type Element
        Character(len=2)                        :: ChemSymb
        Character(len=10)                       :: Lab
        Integer                                 :: Z
        Integer                                 :: Type
     End Type Element

     Type(Element), Dimension(:), Allocatable   :: Elem

     ERR_Char  = .False.
     WARN_Char = .False.

     If (Present(filcod)) Then
        Open(unit=lun,file=Trim(filcod)//".fch",status="replace",action="write")
     Else
        Open(unit=lun,file="Formal_Charges.fch", status="replace",action="write")
     End If

     Write(unit=lun,fmt="(/,/,7(a,/))")                              &
          "             ============================"           , &
          "             ====== FORMAL CHARGES ======"           , &
          "             ============================"           , &
          "    ***********************************************" , &
          "    *  Formal charges from  *.cfl or *.cif files  *" , &
          "    ***********************************************" , &
          "          (JRC - ILL, version: January 2016 )"
     Write(unit=lun,fmt="(a,/)") " "

     ! -------------------
     ! Build the unit cell
     ! -------------------

     Call Allocate_Atoms_Cell(A%Natoms,SpGr%Multip,0.,Ac)
     Call Set_Epsg(0.001) !needed for well controlling the calculation of multiplicities

     Write(unit=lun,fmt="(a,i3)") "    Space Group Number: ", SpGr%NumSpg
     Write(unit=lun,fmt="(a,a,/)") "    Space Group Symbol: ", SpGr%Spg_Symb

     Ac%Nat = 0
     j      = 1

     Do i = 1 , A%Natoms
        Ac%Nat = Ac%Nat + A%Atom(i)%Mult
        Ac%Noms(j:Ac%Nat) = A%Atom(i)%Lab
        Call Get_Orbit(A%Atom(i)%X,SpGr,orbit_mult,Ac%XYZ(:,j:Ac%Nat))

        If (orbit_mult .Ne. A%Atom(i)%Mult) Then
           ERR_Char = .True.
           ERR_Char_Mess = "Error in multiplicities. Check space group"
           Return
        End If
     
        j = j + A%Atom(i)%Mult
     End Do

     ! -----------------------------------
     ! Classify ions as cations and anions
     ! -----------------------------------

     Call Set_Chem_Info()
     Call Set_Atomic_Properties()

     ! - Determine elements

     nelements         = 0
     nanions_candidate = 0

     Allocate(anion_candidate(A%natoms),atom2elem(A%natoms))

     Do i = 1 , A%Natoms
        newelement         = .True.
        anion_candidate(i) = 0
        label = A%Atom(i)%ChemSymb
        Call Get_ChemSymb(label,A%Atom(i)%ChemSymb,A%Atom(i)%Z)
        Do j = 1 , i - 1
           If (A%Atom(i)%ChemSymb == A%Atom(j)%ChemSymb) Then 
              newelement = .False.
              atom2elem(i) = atom2elem(j)
           End If
        End Do
        If (newelement) Then 
           nelements    = nelements + 1
           atom2elem(i) = nelements
        End If
     End Do

     Write(unit=lun,fmt="(a,i2/)") " => Number of elements found: ", nelements

     ! - Determine number of sites

     nsites = 0
  
     Allocate(atom2site(A%natoms))
  
     Do i = 1 , A%Natoms
        newsite = .True.
        Do j = 1 , i -1
           If (Abs(A%Atom(i)%X(1)-A%Atom(j)%X(1)) < 0.01 .And. &
               Abs(A%Atom(i)%X(2)-A%Atom(j)%X(2)) < 0.01 .And. &
               Abs(A%Atom(i)%X(3)-A%Atom(j)%X(3)) < 0.01) Then
              newsite = .False.
              atom2site(i) = atom2site(j)
           End If
        End Do
        If (newsite) Then
           nsites = nsites + 1
           atom2site(i) = nsites
        End If
     End Do

     Allocate(siteocc(nsites))
  
     siteocc(:) = 0.

     Do i = 1 , A%Natoms
        siteocc(atom2site(i)) = siteocc(atom2site(i)) + A%Atom(i)%Occ 
     End Do

     Write(unit=lun,fmt="(a,i2/)") " => Number of crystallographic sites found: ", nsites

     ! - Determine which elements are transition metal or rare earths

     Allocate(tmre(A%Natoms))
  
     tmre(:) = 0

     Do i = 1 , A%Natoms
        If ( (A%Atom(i)%Z > 20 .And. A%Atom(i)%Z < 31) .Or. &
             (A%Atom(i)%Z > 38 .And. A%Atom(i)%Z < 49) .Or. &
             (A%Atom(i)%Z > 56 .And. A%Atom(i)%Z < 81) .Or. &
             (A%Atom(i)%Z > 88)) tmre(i) = 1
     End Do

     ! - Extract all possible valences 

     Write(unit=lun,fmt="(a,/)") " => Extracting valences from Atomic Properties Table...."

     Allocate(valence_all(nmax_val,A%Natoms),nvalences_all(A%Natoms))

     nvalences_all(:) = 0
  
     Do i = 1 , A%Natoms
        j = 1
        Do
           If (Ap_Table(j)%Z == A%Atom(i)%Z) Then
              Exit
           Else
              j = j + 1
           End If
        End Do
        Do
           If (Ap_Table(j)%Z == A%Atom(i)%Z) Then
              nvalences_all(i) = nvalences_all(i) + 1
              valence_all(nvalences_all(i),i) = Ap_Table(j)%oxs
              j = j + 1
           Else
              Exit
           End If
        End Do
     End Do

     ! - Determine anions candidates

     Write(unit=lun,fmt="(a,/)") " => Searching ions that can be anions...."
     
     Do i = 1 , A%Natoms
        j = 1
        Do 
           If (valence_all(j,i) < 0) Then 
              anion_candidate(i) = -1
              nanions_candidate  = nanions_candidate + 1
           End If
           j = j + 1
           If (j > nvalences_all(i) .Or. &
                anion_candidate(i) == -1) Exit
        End Do
     End Do

     If (nanions_candidate == 0) Then
        ERR_Char = .True.
        ERR_Char_Mess = "Error in anions candidates. No anion candidate has been found"
        Return
     End If

     Write(unit=lun,fmt="(4x,a)") "----------------"
     Write(unit=lun,fmt="(4x,a)") "Anion candidates"
     Write(unit=lun,fmt="(4x,a)") "----------------"

     Do i = 1 , A%natoms
        If (anion_candidate(i) == 0) Then
           Write(unit=lun,fmt="(4x,a4,9x,a)") A%Atom(i)%Lab, " NO"
        Else
           Write(unit=lun,fmt="(4x,a4,9x,a)") A%Atom(i)%Lab, "YES"
        End If
     End Do

     Write(unit=lun,fmt="(4x,a,/)") "----------------"

     ! - Determine first neighbor of anions candidates

     Allocate(first_neighbor(A%Natoms))
     first_neighbor(:) = 0

     ! First neighbors

     Write(unit=lun,fmt="(a,/)") " => Computing first neighbor of anion candidates (needed to define cations and anions)...."
     Write(unit=lun,fmt="(4x,a)")  "-----------------------"
     Write(unit=lun,fmt="(4x,a)")  "First neighbor distance"
     Write(unit=lun,fmt="(4x,a)")  "-----------------------"

     Call Calc_Dist_Angle_Sigma(Dmax,Dangle,Cell,SpGr,A)  
  
     Do i = 1 , A%natoms
        rmin = 1d8
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
           ERR_Char_Mess = "Error in coordination. Coordination for atom "//Trim(A%Atom(i)%Lab)//" is zero. Increase Dmax"   
           Return
        End If
        Write(unit=lun,fmt="(4x,a4,2x,a4,7x,f6.3)") A%Atom(i)%Lab,A%Atom(first_neighbor(i))%Lab,rmin
     End Do

     Write(unit=lun,fmt="(4x,a,/)") "-----------------------"

     ! - Classify ions as anions and cations

     Call Set_Pauling_Electronegativity()

     Allocate(iontype(A%natoms))
     
     iontype(:) = 0

     ! Ions with anion_candidate = 0 be cations
     ! Ions with anion_candidate < 0 can be cations or anions
     !     - if the first neighbor anion_candidate = 0, it is an anion
     !     - if the first neighbor anion_candidate .ne. 0, 
     !          it is a cation if it has lower electronegativity
     !          it is an anion otherwise

     Do i = 1 , A%Natoms
        If (anion_candidate(i) == 0) iontype(i) = 1
     End Do

     Do i = 1 , A%Natoms
        If (anion_candidate(i) == -1) Then
           If (anion_candidate(first_neighbor(i)) == 0) Then
              If (iontype(i) <= 0) Then
                 iontype(i) = -1
              Else
                 ERR_Char = .True.
                 ERR_Char_Mess = "Error in classification. Atom "//Trim(A%Atom(i)%Lab)//" has been classified as anion and cation"
                 Return
              End If
           Else
              If (PaulingX(A%Atom(i)%Z) < PaulingX(A%Atom(first_neighbor(i))%Z)) Then
                 If (iontype(i) >= 0) Then
                    iontype(i) = 1
                 Else
                    ERR_Char = .True.
                    ERR_Char_Mess = "Error in classification. Atom "&
                         //Trim(A%Atom(i)%Lab)//" has been classified as anion and cation"
                    Return
                 End If
                 If (iontype(first_neighbor(i)) == 1) Then
                    ERR_Char = .True.
                    ERR_Char_Mess = "Error in classification. Atom "&
                         //Trim(A%Atom(first_neighbor(i))%Lab)//" has been classified as anion and cation"
                    Return
                 Else 
                    iontype(first_neighbor(i)) = -1
                 End If
              Else
                 If (iontype(i) <= 0) Then
                    iontype(i) = -1
                 Else
                    ERR_Char = .True.
                    ERR_Char_Mess = "Error in classification. Atom "&
                         //Trim(A%Atom(i)%Lab)//" has been classified as anion and cation"
                    Return
                 End If
                 If (iontype(first_neighbor(i)) == -1) Then
                    ERR_Char = .True.
                    ERR_Char_Mess = "Error in classification. Atom "&
                         //Trim(A%Atom(first_neighbor(i))%Lab)//" has been classified as anion and cation"
                    Return
                 Else 
                    iontype(first_neighbor(i)) = 1
                 End If
              End If
           End If
        End If
     End Do

     Write(unit=lun,fmt="(a,/)") " => Anion / Cation assignment...."
     Write(unit=lun,fmt="(4x,a)") "----------------"
     Write(unit=lun,fmt="(4x,a)") "Anions / Cations"
     Write(unit=lun,fmt="(4x,a)") "----------------"

     Do i = 1 , A%natoms
        
        If (iontype(i) == 1) Then
           Write(unit=lun,fmt="(4x,a4,6x,a)") A%Atom(i)%Lab, "CATION"
        Else If (iontype(i) == -1) Then
           Write(unit=lun,fmt="(4x,a4,6x,a)") A%Atom(i)%Lab, " ANION"
        End If
      
     End Do

     Write(unit=lun,fmt="(4x,a,/)") "----------------"

     ! - Get cationic and anionic valences

     Allocate(valence(nmax_val,A%Natoms),nvalences(A%Natoms))

     nvalences(:) = 0

     Do i = 1 , A%Natoms
        Do j = 1 , nvalences_all(i)
           If ( (iontype(i) ==  1 .And. valence_all(j,i) > 0) .Or. &
                (iontype(i) == -1 .And. valence_all(j,i) < 0) ) Then
              nvalences(i) = nvalences(i) + 1
              valence(nvalences(i),i) = valence_all(j,i)
           End If
        End Do
        If (A%Atom(i)%Z == 8) nvalences(i) = 1 ! Correction for oxygen, there are two valences = -2 in Atomic Properties Table
     End Do

     Write(unit=lun,fmt="(a,/)") " => Possible atomic valences...."
     Write(unit=lun,fmt="(4x,a)") "-----------------------------------------"
     Write(unit=lun,fmt="(4x,a)") "Ion, Valences"
     Write(unit=lun,fmt="(4x,a)") "-----------------------------------------"

     Do i = 1 , A%Natoms
        Write(unit=lun,fmt="(4x,a4)",advance='no') A%Atom(i)%ChemSymb
        Do j = 1 , nvalences(i)
           Write(unit=lun,fmt="(3x,i2,1x)",advance='no') valence(j,i)
        End Do
        Write(unit=lun,fmt=*)
     End Do

     Write(unit=lun,fmt="(4x,a,/)") "-----------------------------------------"

     ! ----------------------------------------------
     ! Find solutions which fulfill electroneutrality
     ! ----------------------------------------------

     Write(unit=lun,fmt="(a,/)") " => Searching for combinations that fulfill the electroneutrality rule..."   

     ncombinations=1
      
     Do i = 1 , A%Natoms
        ncombinations = ncombinations * nvalences(i)
     End Do

     Allocate(combination(A%Natoms,ncombinations))
     Allocate(sumq(ncombinations))
     Allocate(contador(A%Natoms))

     nsolutions  = 0
     sumq_min    = 1d8
     contador(:) = 1

     Do icomb = 1 , ncombinations
        total_charge = 0.
        Do j = 1 , A%Natoms
           combination(j,icomb) = valence(contador(j),j)
           total_charge         = total_charge + valence(contador(j),j) * A%atom(j)%occ
        End Do
        sumq(icomb) = Abs(total_charge)
        If (Abs(total_charge) < sumq_min) sumq_min = Abs(total_charge)
        If (icomb < ncombinations) Then
           incremento = .False.
           j = A%Natoms
           Do
              If (contador(j) < nvalences(j)) Then
                 contador(j) = contador(j) + 1
                 incremento  = .True.
              Else 
                 contador(j) = 1
                 j = j - 1
              End If
              If (incremento) Exit
           End Do
        End If
     End Do

     If (sumq_min <= zerocharge) Then
        electroneutrality = .True.
     Else
        electroneutrality = .False.
        Write(unit=lun,fmt="(8x,a)") "No solution found fulfilling electroneutrality. Possible problems:"
        Write(unit=lun,fmt="(12x,a)") " - wrong compositions"
        Write(unit=lun,fmt="(12x,a)") " - fractional charge --> mixed valence oxidation state"
        Write(unit=lun,fmt="(12x,a)") " - unusual oxidation state (not tabulated)"
        WARN_Char = .True.
        WARN_Char_Mess = "WARNING! No solution found fulfilling electroneutrality"
     End If

     Write(unit=lun,fmt="(a,/)") " => Selecting combinations...."

     nsolutions = 0
     Do i = 1 , ncombinations
        If (sumq(i) <= sumq_min + DQ) nsolutions = nsolutions + 1
     End Do

     Allocate(sumq_old(ncombinations))
     sumq_old(:) = sumq(:)
     Deallocate(sumq)
     Allocate(solution(A%Natoms,nsolutions),sumq(nsolutions))

     nsolutions = 0
     Do i = 1 , ncombinations
        If (sumq_old(i) <= sumq_min+0.01) Then
           nsolutions = nsolutions + 1
           solution(:,nsolutions) = combination(:,i)
           sumq(nsolutions) = sumq_old(i)
        End If
     End Do

     Deallocate(combination)
     Deallocate(sumq_old)

     Write(unit=lun,fmt="(4x,i10,2x,a,/)") nsolutions, "solutions found"
     Write(unit=lun,fmt="(4x,a,/)") "Computing penalizations...."
     Write(unit=lun,fmt="(8x,a,/)") "For each ion, atomic valence and (bond valence sum)"
     
     ! ---------------------
     ! Compute penalizations
     ! ---------------------

     Call Allocate_Atoms_Conf_List(A%Natoms,Acl)
     Allocate(bvsum(A%Natoms,nsolutions),penalty(nsolutions))

     ! if for a given crystallographic site, there are 
     ! two or more transition/rare earth metals, we apply
     ! a penalization if the difference between the maximum
     ! and minimum valence >= 3
     
     Do is = 1 , nsolutions
        penalty(is) = 0.
        Do i = 1 , nsites
           maxv = 0
           minv = 0
           k    = 0
           Do j = 1 , A%Natoms
              If (tmre(j) == 1 .And. atom2site(j) == i) Then
                 k = k + 1
                 If (k == 1) Then
                    maxv = solution(j,is)
                    minv = solution(j,is)
                 Else
                    If (solution(j,is) > maxv) maxv = solution(j,is)
                    If (solution(j,is) < minv) minv = solution(j,is)
                 End If
              End If
           End Do
           If ((maxv-minv) > 2) penalty(is) = penalty(is) + maxv - minv
        End Do
     End Do

     Do is = 1 , nsolutions
        Acl%atom(1:A%natoms)=A%atom(:)
        Do j = 1 , A%Natoms
           Acl%Atom(j)%Charge = solution(j,is)
        End Do

        soft=.False.
        Call Init_Err_Conf
        Call Species_on_List(Acl,SpGr%Multip,ttol)     
        Call Set_Table_d0_b(Acl)
        If (Err_Conf) Then
           soft=.True.
           Call Init_Err_Conf
           Call Species_on_List(Acl,SpGr%Multip,ttol,soft)     
           Call Set_Table_d0_b(Acl,soft=soft)
        End If
        If (Err_Conf) Return

        ! Compute BVS sums

        tol = Acl%tol*0.01
        If (tol <= 0.001) tol=0.20

        Do i = 1 , Acl%natoms
           icm         = coord_info%coord_num(i)
           l           = Acl%Atom(i)%ind(1)
           q1          = Acl%Atom(i)%charge
           sig1        = Sign(1.0_cp,q1)
           icn         = 0
           dav         = 0.
           bvsum(i,is) = 0.0
           Do j = 1 , icm
              k           = Acl%Atom(coord_info%n_cooatm(j,i))%ind(1)
              q2          = Acl%Atom(coord_info%n_cooatm(j,i))%charge
              sig2        = Sign(1.0_cp,q2)
              If (sig1 == sig2) Cycle
              dd          = coord_info%dist(j,i)
              If (dd > (Acl%Radius(l)+Acl%Radius(k))*(1.0+tol)) Cycle
              icn         = icn + 1
              dav         = dav + dd
              s0          = coord_info%s_dist(j,i)
              str         = Exp((Table_d0(l,k)-dd)/Table_b(l,k))
              str         = str*Acl%Atom(coord_info%n_cooatm(j,i))%VarF(1)
              bvsum(i,is) = bvsum(i,is)+str
           End Do
           penalty(is) = penalty(is) + &
                (Abs(bvsum(i,is) - Abs(Acl%Atom(i)%Charge))) * Acl%Atom(i)%Occ / siteocc(atom2site(i))
        End Do
     End Do

     ! - Order solutions according to the penalty

     Allocate(solution_aux(A%Natoms),bvsum_aux(A%Natoms))
     ordered = .False.

     Do 
        ordered = .True.
        Do i = 1 , nsolutions - 1
           If (penalty(i+1) < penalty(i)) Then
              ordered = .False.
              penalty_aux     = penalty(i+1)
              sumq_aux        = sumq(i+1)
              solution_aux(:) = solution(:,i+1)
              bvsum_aux(:)    = bvsum(:,i+1)
              penalty(i+1)    = penalty(i)
              sumq(i+1)       = sumq(i)
              solution(:,i+1) = solution(:,i)
              bvsum(:,i+1)    = bvsum(:,i)
              penalty(i)      = penalty_aux
              sumq(i)         = sumq_aux
              solution(:,i)   = solution_aux(:)
              bvsum(:,i)      = bvsum_aux(:)
           End If
        End Do
        If (ordered) Exit
     End Do

     Write(unit=lun,fmt="(5x)",advance="no")

     Do i = 1 , A%Natoms
        Write(unit=lun,fmt="(a)",advance='no') "------------------"
     End Do
     Write(unit=lun,fmt="(a)") "--------------------------"
     
     Do i = 1 , A%Natoms
        Write(unit=lun,fmt="(5x,a4,9x)",advance="no") A%Atom(i)%Lab
     End Do
     
     Write(unit=lun,fmt="(10x,a)") "Penalty  Total_Charge"
     Write(unit=lun,fmt="(5x)",advance="no")
     
     Do i = 1 , A%Natoms
        Write(unit=lun,fmt="(a)",advance='no') "------------------"
     End Do
     Write(unit=lun,fmt="(a)") "--------------------------"
     
     Do i = 1 , nsolutions
        Write(unit=lun,fmt="(4x)",advance="no")
        
        Do j = 1 , A%Natoms
           Write(unit=lun,fmt="(f5.2,1x,a1,f6.4,a1,4x)",advance='no') &
                solution(j,i),"(",bvsum(j,i),")"
        End Do
        
        Write(unit=lun,fmt="(5x,f8.3,f14.5)") penalty(i),sumq(i)
     End Do
     
     Write(unit=lun,fmt="(5x)",advance="no")
     
     Do i = 1 , A%Natoms
        Write(unit=lun,fmt="(a)",advance='no') "------------------"
     End Do
     Write(unit=lun,fmt="(a,/)") "--------------------------"
     
     charges(1:A%Natoms) = solution(1:A%Natoms,1)

   End Subroutine Set_Charges

 End Module CFML_Formal_Charges

