!!----
!!---- Program Formal_Charges
!!----
!!---- Using CrysFML versions above 4.01
!!----
Program Formal_Charges

  !---- Use Modules ----!
  Use CFML_GlobalDeps
  Use CFML_String_Utilities
  Use CFML_Math_General,                Only: Set_Epsg
  Use CFML_Crystal_Metrics,             Only: Crystal_Cell_Type
  Use CFML_Crystallographic_Symmetry,   Only: Space_Group_Type, Get_Orbit
  Use CFML_Atom_TypeDef,                Only: Atom_List_Type, Atoms_Cell_Type, Allocate_Atoms_Cell, Atom_List_To_Cell
  Use CFML_IO_Formats,                  Only: Readn_set_Xtal_Structure,ERR_Form_Mess,err_form,File_List_Type,Write_CFL
  Use CFML_Scattering_Chemical_Tables,  Only: Get_ChemSymb, Set_Chem_Info
  Use CFML_Geometry_Calc
  Use CFML_BVSpar
  Use CFML_Atomic_Data

  !---- Variables ----!
  Implicit None

  Type (Space_Group_Type)                    :: SpGr
  Type (Crystal_Cell_Type)                   :: Cell
  Type (Atom_List_Type)                      :: A
  Type (Atoms_Cell_Type)                     :: Ac, Bc
  Type (File_List_Type)                      :: Fich_cfl
 
  Character(len=2)                           :: label
  Character(len=256)                         :: filcod

  Integer                                    :: i, j, k, l, n, m, icomb, p, diff, penalty_aux, orbit_mult
  Integer                                    :: narg, nelements, nanions, nanions_candidate, ncombinations, nsolutions, nunusual_aux
  Integer                                    :: nsol_usual, ntotal_atoms, nelemmix, nequiv_sol, npenalty_min, penalty_min
  Integer                                    :: lun=1,i_cfl=2,err=3
  Integer, Parameter                         :: massive_limit = 100
  Integer, Dimension(:), Allocatable         :: anion_candidate, first_neighbor, iontype, nvalences, contador, nunusual
  Integer, Dimension(:), Allocatable         :: atom2elem, elemmix
  Integer, Dimension(:), Allocatable         :: penalty, Ac2Bc
  Integer, Dimension(:,:), Allocatable       :: valence

  Real(kind=cp)                              :: rmin, total_charge, energy_aux, norm
  Real(kind=cp), Parameter                   :: ZEROCHARGE = 1d-4
  Real(kind=cp)                              :: Dmax=4.0, Dangle=0.
  Real(kind=8),  Dimension(:), Allocatable   :: energy
  Real(kind=cp), Dimension(:), Allocatable   :: sumocc, solution_aux
  Real(kind=cp), Dimension(:,:), Allocatable :: solution, combination, solution_new
  
  Logical                                    :: arggiven=.False.,esta,cif=.False.,out_cif=.False.
  Logical                                    :: error_pairs=.False.,error_element=.False.
  Logical                                    :: warning_valence
  Logical                                    :: newelement, incremento, usual, ordered, xyz_stored, newsolution, identical

  Type Element
     Character(len=2)                        :: ChemSymb
     Character(len=10)                       :: Lab
     Integer                                 :: Z
     Integer                                 :: type
  End Type Element

  Type(Element), Dimension(:), Allocatable   :: Elem

  ! Arguments on the command line

  narg=Command_Argument_Count()

  If (narg > 0) Then
     Call GET_COMMAND_ARGUMENT(1,filcod)
     i=Index(filcod,".cfl")
     If(i /= 0) filcod=filcod(1:i-1)
     i=Index(filcod,".cif")
     If(i /= 0) filcod=filcod(1:i-1)
     arggiven=.True.
  End If

  Write(unit=*,fmt="(/,/,7(a,/))")                                &
       "             ============================"           , &
       "             ====== FORMAL CHARGES ======"           , &
       "             ============================"           , &
       "    ***********************************************" , &
       "    *  Formal charges from  *.cfl or *.cif files  *" , &
       "    ***********************************************" , &
       "          (JRC - ILL, version: January 2016 )"
   
  Write(unit=*,fmt=*) " "

  If(.Not. arggiven) Then
     Write(unit=*,fmt="(a)",advance='no') " => Code of the file xx.cfl(cif) (give xx): "
     Read(unit=*,fmt="(a)") filcod
     If(Len_trim(filcod) == 0) Stop
  End If

  Inquire(file=Trim(filcod)//".cfl",exist=esta)

  If(esta) Then
     Call Readn_set_Xtal_Structure(Trim(filcod)//".cfl",Cell,SpGr,A,Mode="CFL",file_list=fich_cfl)
     cif=.False.
     out_cif=.True.
  Else
     Inquire(file=Trim(filcod)//".cif",exist=esta)
     If(.Not. esta) Then
        Write(unit=*,fmt="(a)") " File: "//Trim(filcod)//".cfl (or .cif) does'nt exist!"
        stop
     End If
     Call Readn_set_Xtal_Structure(Trim(filcod)//".cif",Cell,SpGr,A,Mode="CIF",file_list=fich_cfl)
     !Call Readn_set_Xtal_Structure(Trim(filcod)//".cif",Cell,SpGr,A,Mode="CIF")
     cif=.True.
  End If

  !
  If (err_form) Then
     Write(unit=*,fmt="(a)") Trim(ERR_Form_Mess)
     Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
     Stop
  End If
  
  Open(unit=lun,file=Trim(filcod)//".fch", status="replace",action="write")
  Write(unit=lun,fmt="(/,/,7(a,/))")                              &
       "             ============================"           , &
       "             ====== FORMAL CHARGES ======"           , &
       "             ============================"           , &
       "    ***********************************************" , &
       "    *  Formal charges from  *.cfl or *.cif files  *" , &
       "    ***********************************************" , &
       "          (JRC - ILL, version: January 2016 )"
  Write(unit=lun,fmt="(a,/)") " "

  If (cif) Then
     Write(unit=lun,fmt="(a,/)") " => Data obtained from CIF file: "//Trim(filcod)//".cif"     
  Else
     Write(unit=lun,fmt="(a,/)") " => Content of the input file: "
     Do i=1,fich_cfl%nlines
        Write(unit=lun,fmt="(tr10,a)") fich_cfl%line(i)
     End Do
  End If

  ! We build the unit cell

  Call Allocate_Atoms_Cell(A%Natoms,SpGr%Multip,0.,Ac)
  Call Set_Epsg(0.01) !needed for well controlling the calculation of multiplicities

  Write(unit=*,fmt="(a,i3)") "    Space Group Number: ", SpGr%NumSpg
  Write(unit=*,fmt="(a,a,/)") "    Space Group Symbol: ", SpGr%Spg_Symb

  Write(unit=lun,fmt="(a,i3)") "    Space Group Number: ", SpGr%NumSpg
  Write(unit=lun,fmt="(a,a,/)") "    Space Group Symbol: ", SpGr%Spg_Symb

  Ac%Nat = 0
  j      = 1

  Do i = 1 , A%Natoms
     Ac%Nat = Ac%Nat + A%Atom(i)%Mult
     Ac%Noms(j:Ac%Nat) = A%Atom(i)%Lab
     Call Get_Orbit(A%Atom(i)%X,SpGr,orbit_mult,Ac%XYZ(:,j:Ac%Nat))

     If (orbit_mult .Ne. A%Atom(i)%Mult) Then
        Open(unit=err,file=Trim(filcod)//".err", status="replace",action="write")
        Write(unit=err,fmt="(a)") " => Error in multiplicities. Increase epsg"
        Write(unit=err,fmt="(a,i4,a1)")  "    Atom ", i, ":"
        Write(unit=err,fmt="(a,i4)")     "                Multiplicity in cfl(cif) file:",   A%Atom(i)%Mult
        Write(unit=err,fmt="(a,i4)")     "                Multiplicity found by Get_Orbit:", orbit_mult
        Close(unit=err)
        Close(unit=lun)
        Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
        Write(unit=*,fmt="(a)")   " => Check "//Trim(filcod)//".err file"
        Stop
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
        j         = 1
        k         = k + 1

        Do !While (.Not. xyz_stored .And. j <= Bc%Nat)
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
  Write(unit=*,fmt="(a,/)") " => Searching ions that can be anions...."
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
  
  If (nelements == 1) Then

     ! We have an element, not a compound
     ! Formal charges = 0.

     Do i = 1 , A%nAtoms
        A%Atom(i)%Charge = 0.
     End Do

     Write(unit=lun,fmt="(/,a)") " => This is an element. All charges are have been set to zero"
     Write(unit=*,fmt="(/,a)") " => This is an element. All charges are have been set to zero"

  Else If (nanions_candidate == 0) Then

     ! No anions 
     ! Formal charges = 0.
     
     Do i = 1 , A%nAtoms
        A%Atom(i)%Charge = 0.
     End Do

     Write(unit=lun,fmt="(/,a)") " => No anions candidates have been found"
     Write(unit=lun,fmt="(4x,a)") "All charges have been set to zero"
     Write(unit=*,fmt="(/,a)") " => No anions candidates have been found"
     Write(unit=*,fmt="(4x,a)") "All charges have been set to zero"
      
  Else
      
     ! Here we determine which ion is a cation, which is an anion
     !     1. we determine the first neighbor of each atom
     !     2. each pair must be a cation-anion pair

     Call Calc_Dist_Angle(Dmax,Dangle,Cell,SpGr,A)
     Allocate(first_neighbor(A%natoms))
     first_neighbor = 0

     ! First neighbors

     Write(unit=lun,fmt="(a,/)") " => Computing first neighbor distance, needed to define cations and anions..."
     Write(unit=*,fmt="(a,/)") " => Computing first neighbor distance, needed to define cations and anions..."
     Write(unit=lun,fmt="(4x,a)") "-----------------------"
     Write(unit=lun,fmt="(4x,a)")  "First neighbor distance"
     Write(unit=lun,fmt="(4x,a)")  "-----------------------"
      
     Do i = 1 , A%natoms
        rmin = 1d8

        Do j = 1 , coord_info%Max_coor

           If (coord_info%N_CooAtm(j,i) .Ne. 0) Then 

              If (coord_info%Dist(j,i) < rmin) Then
                 rmin = coord_info%Dist(j,i)
                 first_neighbor(i) = coord_info%N_CooAtm(j,i)
              End If

           End If

        End Do
        
        If (first_neighbor(i) == 0) Then
           Open(unit=err,file=Trim(filcod)//".err", status="replace",action="write")
           Write(unit=err,fmt="(a)") " => Error in coordination "
           Write(unit=err,fmt="(a,1x,a2,1x,a)") "    Coordination for atom", A%Atom(i)%ChemSymb, "is zero. Increase Dmax" 
           Close(unit=err)
           Close(unit=lun)
           Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
           Write(unit=*,fmt="(a)")   " => Check "//Trim(filcod)//".err file"
           Stop
        End If

        Write(unit=lun,fmt="(4x,a4,2x,a4,7x,f6.3)") A%Atom(i)%Lab,A%Atom(first_neighbor(i))%Lab,rmin
     End Do

     Write(unit=lun,fmt="(4x,a,/)") "-----------------------"
     
     ! Setting anion/cation nature

     Call Set_Pauling_Electronegativity()

     Allocate(iontype(A%natoms))
   
     iontype(:) = 0

     ! Ions with anion_candidate = 0 be cations
     ! Ions with anion_candidate greater than zero can be cations or anions
     !     - if the first neighbor anion_candidate = 0, it is an anion
     !     - if the first neighbor anion_candidate .ne. 0, 
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
     Write(unit=*,fmt="(a,/)") " => Anion / Cation assignment...."
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
        Open(unit=err,file=Trim(filcod)//".err", status="replace",action="write")
        Write(unit=err,fmt="(a)") " => Error in first-neighbor "
        Write(unit=err,fmt="(a)") "    Cation-Cation or Anion-Anion first-neighbor pairs has been found"
        Close(unit=err)
        Close(unit=lun)
        Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
        Write(unit=*,fmt="(a)")   " => Check "//Trim(filcod)//".err file"
        Stop
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
              Write(unit=lun,fmt="(a,i4,1x,a,i4,1x,a)") " => Error! Atom", j, "and Atom", l, "must be both cations or anions"
              error_element = .True.
           End If

        End Do

     End Do

     If (error_element) Then 
        Open(unit=err,file=Trim(filcod)//".err", status="replace",action="write")
        Write(unit=err,fmt="(a)") " => Error in element"
        Write(unit=err,fmt="(a)") "    Cations and anions found for an element"
        Close(unit=err)
        Close(unit=lun)
        Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
        Write(unit=*,fmt="(a)")   " => Check "//Trim(filcod)//".err file"
        Stop
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
     Write(unit=*,fmt="(a,/)") " => Extracting valences from Oxidation States Table...."
     Write(unit=lun,fmt="(4x,a)") "----------------------------------------------"
     Write(unit=lun,fmt="(4x,a)") "Ion, Usual / (Reported) valences"
     Write(unit=lun,fmt="(4x,a)") "----------------------------------------------"

     Do i = 1 , nelements
        If (elem(i)%Z == 8) nvalences(i) = 1

        Write(unit=lun,fmt="(4x,a4)",advance='no') elem(i)%ChemSymb

        Do j = 1 , nvalences(i)
           k = 1 
           usual = .False.

           Do !While (.Not. usual .And. k <= 8)
              If (valence(j,i) == Common_OxStates_Table(k,elem(i)%Z)) usual = .True.
              k = k + 1
              If (usual .Or. k > 8) Exit
           End Do

           If (usual) Then 
              Write(unit=lun,fmt="(3x,i2,1x)",advance='no') valence(j,i)
           Else
              Write(unit=lun,fmt="(2x,a1,i2,a1)",advance='no') "(", valence(j,i), ")"
           End If

        End Do
         
        Write(unit=lun,fmt=*)
     End Do

     Write(unit=lun,fmt="(4x,a,/)") "----------------------------------------------"
     
     Write(unit=lun,fmt="(a,/)") " => Searching for combinations that fulfill the electroneutrality rule..."   
     Write(unit=*,fmt="(a,/)") " => Searching for combinations that fulfill the electroneutrality rule..."   

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
           total_charge         = total_charge + valence(contador(atom2elem(j)),atom2elem(j)) * A%atom(j)%occ
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
        !Open(unit=err,file=Trim(filcod)//".err", status="replace",action="write")
        !Write(unit=err,fmt="(a)") " => Error in valences "
        !Write(unit=err,fmt="(a)") "    No solution found with integer valences "
        !Close(unit=err)
        !Close(unit=lun)
        !Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
        !Write(unit=*,fmt="(a)")   " => Check "//Trim(filcod)//".err file"
        !Stop

        Write(unit=*,fmt="(a)")         "    No solution found with integer valences "
        Write(unit=*,fmt="(a,i3,2x,a)") "    Trying solutions with non-integer valences (mixed valence)"
        Write(unit=*,fmt="()")
        Write(unit=lun,fmt="(a)")         "    No solution found with integer valences "
        Write(unit=lun,fmt="(a,i3,2x,a)") "    Trying solutions with non-integer valences (mixed valence)"
        Write(unit=lun,fmt="()")

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
           Write(unit=*,fmt="(a,i3,2x,a)") "    No transition / earth-rare ion found. Unable to set charges"        
           Stop
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
                 m           = 1 
                 
                 Do !While (newsolution .And. m <= nsolutions)
                    If (newsolution .Or. m > nsolutions) Exit
                    identical = .True.

                    Do k = 1 , nelements
                       If (solution(k,nsolutions+1) .Ne. solution(k,m)) identical = .False.
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
        Write(unit=*,fmt="(a,i3,2x,a)") " => No solution found"
        Write(unit=*,fmt="(a,i3,2x,a)") "    Unable to set charges"
        Write(unit=*,fmt="()")

        Write(unit=lun,fmt="(a,i3,2x,a)") "    No solution found"
        Write(unit=lun,fmt="(a,i3,2x,a)") "    Unable to set charges"
        Write(unit=lun,fmt="()")
        Stop
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

     Write(unit=*,fmt="(4x,i10,2x,a)") , nsolutions, "sets of formal charges found"
     Write(unit=*,fmt="(4x,i10,2x,a,/)") nsol_usual, "sets with usual formal charges"
     Write(unit=*,fmt="(a,/)") " => Sets have been written in "//Trim(filcod)//".fch"         
     
     Write(unit=lun,fmt="(4x,i10,2x,a)") nsolutions, "solutions found"
     Write(unit=lun,fmt="(4x,i10,2x,a,/)") nsol_usual, "solutions with usual formal charges (penalty = 0)"

     If (nsolutions > massive_limit) Then
        Write(unit=lun,fmt="(4x,a)") "number of solutions too large"
        Write(unit=lun,fmt="(4x,a,/)") "work only with solutions with the minimum penalty"
        
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
     Write(unit=*,fmt="(a,/)") " => Computing energies..."   

     Allocate(energy(nsolutions))

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

        Do !While (j <= nsolutions)
           If (j > nsolutions) Exit
           If (penalty(m) == penalty(j)) n = n + 1
           j = j + 1
        End Do
        
        If (n .Ne. m) Then
           ordered = .False.

           Do !While (.Not. ordered)
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
     
     Write(unit=lun,fmt="(5x)",advance="no")

     Do i = 1 , nelements
        Write(unit=lun,fmt="(a)",advance='no') "---------"
     End Do
     Write(unit=lun,fmt="(a)") "---------------------"
     
     Do i = 1 , nelements
        Write(unit=lun,fmt="(5x,a4)",advance="no") elem(i)%ChemSymb
     End Do
     
     Write(unit=lun,fmt="(2x,a)") "Energy (eV/atom) Penalty"
     
     Write(unit=lun,fmt="(5x)",advance="no")
     
     Do i = 1 , nelements
        Write(unit=lun,fmt="(a)",advance='no') "---------"
     End Do
     Write(unit=lun,fmt="(a)") "---------------------"
     
     Do i = 1 , nsolutions
        Write(unit=lun,fmt="(4x)",advance="no")
        
        Do j = 1 , nelements
           Write(unit=lun,fmt="(f5.2,4x)",advance='no') solution(j,i)
        End Do
        
        Write(unit=lun,fmt="(5x,f8.3,i9)") energy(i), penalty(i)
     End Do
     
     Write(unit=lun,fmt="(5x)",advance="no")
     
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
        Write(unit=lun,fmt="(5x,i2,1x,a)") nequiv_sol, "equivalent solutions found."
        Write(unit=lun,fmt="(8x,a,/)") "final solution = average of equivalent solutions"
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
     Write(unit=*,fmt="(a,/)") " => WARNING: Non-integer valences found"
     Write(unit=lun,fmt="(a,/)") " => WARNING: Non-integer valences found"
  End If

  Write(unit=*,fmt="(a,/)") " => Formal Charges: "
  Write(unit=lun,fmt="(a,/)") " => Formal Charges: "
  
  Do i = 1 , A%Natoms
     Write(unit=*,fmt="(4x,a4,a1,1x,f5.2)") A%Atom(i)%Lab, ":", A%Atom(i)%Charge
     Write(unit=lun,fmt="(4x,a4,a1,1x,f5.2)") A%Atom(i)%Lab, ":", A%Atom(i)%Charge
  End Do
  
  Write(unit=*,fmt="()")

  Open(unit=i_cfl,file="CFL_file.cfl",status="replace",action="write")
  Write(unit=i_cfl,fmt="(a)") "Title  CFL-file generated from by Formal_Charges.f90"
  Call Write_CFL(i_cfl,Cell,SpGr,A)
  Close(unit=i_cfl)
  
  Write(unit=*,fmt="(a,/)") " => A CFL-file has been generated -> CFL_file.cfl"

End Program Formal_Charges

Subroutine Ewald(Lattvec,Vol,Ac,e)

  !---- Use Modules ----!
  Use CFML_GlobalDeps,                  Only: Cp
  Use CFML_String_Utilities
  Use CFML_Math_3D
  Use CFML_Crystal_Metrics,             Only: Crystal_Cell_Type
  Use CFML_Atom_TypeDef,                Only: Atoms_Cell_Type

  !---- Variables ----!
  Implicit None

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
  Real(kind=8)                                    rcut2, kcut2, four_beta2, pi, twopi, fourpi
  Real(kind=8)                                 :: r2, e_s, e_l, e_self, fase, sk2
  Real(kind=8),  Dimension(2)                  :: sk
  Real(kind=8),  Dimension(3)                  :: latt, rlatt, rdmin, kdmin, r0, r, a, b, c
  Real(kind=cp), Dimension(3,3)                :: rlattvec
  Real(kind=cp), Dimension(:), Allocatable     :: expok
  Real(kind=cp), Dimension(:,:), Allocatable   :: xyz, errorfc, g

  ! Constants

  rcut2         = rcut ** 2
  kcut2         = kcut ** 2
  pi            = Acos(-1.0d0)
  twopi         = 2. * pi
  fourpi      = 4. * pi
  four_beta2  = 4. * beta ** 2

  ! Compute reciprocal lattice vectors rlattvec
  
  rlattvec(1:3,1) = cross_product(lattvec(1:3,2),lattvec(1:3,3))
  rlattvec(1:3,2) = cross_product(lattvec(1:3,3),lattvec(1:3,1))
  rlattvec(1:3,3) = cross_product(lattvec(1:3,1),lattvec(1:3,2))
  
  rlattvec = rlattvec * twopi / Vol

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

  Do i = 1 , Ac%Nat

     Do j = 1 , 3
        xyz(j,i) = 0

        Do k = 1 , 3
           xyz(j,i) = xyz(j,i) + lattvec(j,k) * Ac%XYZ(k,i)
        End Do
           
     End Do

  End Do

  ! - Arrays construction: ERRORFC, G, EXPOK

  ! - ERROFC

  Allocate(errorfc(Ac%Nat,Ac%Nat))

  Do i = 1 , Ac%Nat 

     Do j = 1 , Ac%Nat
        errorfc(i,j) = 0.
        r0(1:3) = xyz(:,j) - xyz(:,i)
        
        Do h = -ncell(1) , ncell(1)
           
           Do k = -ncell(2) , ncell(2)

              Do l = -ncell(3) , ncell(3)

                 If (i .Ne. j .Or. &
                     h .Ne. 0 .Or. &
                     k .Ne. 0 .Or. &
                     l .Ne. 0) Then
           
                    r(1:3) = r0(1:3) + &
                             h * lattvec(1:3,1) + &
                             k * lattvec(1:3,2) + &
                             l * lattvec(1:3,3)
                 
                    r2 = Dot_Product(r(1:3),r(1:3))

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
     a(1:3) = h * rlattvec(1:3,1)

     Do k = - nrcell(2) , nrcell(2)
        b(1:3) = k * rlattvec(1:3,2)

        Do l = - nrcell(3) , nrcell(3)
           c(1:3) = l * rlattvec(1:3,3)
           r(1:3) = a(1:3) + b(1:3) + c(1:3)
              
           r2 = Dot_Product(r(1:3),r(1:3))

           If (r2 <= kcut2 .And. r2 > 0.0d0) Then
              nk        = nk + 1
              g(nk,1:3) = r(1:3)
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
