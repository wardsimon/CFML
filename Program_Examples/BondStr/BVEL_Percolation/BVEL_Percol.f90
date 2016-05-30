Module Common

  Use CFML_GlobalDeps,       Only:  cp

  Type :: node
     Integer                     :: nbonds, nbonds_pbc    
     Integer                     :: type                ! 1: interior, -1: face, -2: edge, -3: vertex
     Integer                     :: state               ! 0: non-visited, 1: visited 
     Integer                     :: cc, cc_pbc          ! connected component to which the node belongs, before and after pbc
     Integer, Dimension(3)       :: face
     Integer, Dimension(12)      :: bond, bond_pbc
     Real(kind=cp), Dimension(3) :: xyz, xyz_aux
  End Type node

  Type :: component
     Integer                            :: nbonds
     Integer                            :: state
     Integer                            :: cc
     Integer                            :: pbc          ! 0: apply pbc, 1: don't apply
     Integer, Dimension(:), Allocatable :: bond
  End Type component

  Type(node), Dimension(:), Allocatable      :: nodes
  Type(component), Dimension(:), Allocatable :: cts

End Module Common

Program BVEL_Percol

  !---- Use Modules ----!
  Use CFML_GlobalDeps
  Use CFML_Crystal_Metrics,                 Only: Crystal_Cell_Type
  Use CFML_Maps_Calculations,               Only: Calculate_Mesh, Max_Points
  Use Common

  !---- Variables ----!
  Implicit None

  Type (Crystal_Cell_Type)                     :: Cell

  Integer                                      :: i, j, k, l, m, n, v, ic, jc, kc, cc_aux
  Integer                                      :: narg, nlevel, np_max, ncc, nbonds_aux, nbonds_in, nbonds_out, ncc_pbc
  Integer                                      :: lun=1,err=3
  Integer(kind=4)                              :: gType, fType, nval, ndim, nASYM
  Integer(kind=4), Dimension(3)                :: ngrid
  Integer(kind=4), Dimension(4)                :: version
  Integer,         Dimension(0:4)              :: nnodes
  Integer,         Dimension(12)               :: bond_aux
  Integer,         Dimension(2,0:2)            :: tricon
  Integer,         Dimension(:), Allocatable   :: npoints, id

  Real(kind=cp)                                :: Emin, Efin, Eini, dr, dE, dE_max
  Real(kind=cp), Dimension(3)                  :: drmin, t, E_percol
  Real(kind=cp), Dimension(:),     Allocatable :: levels
  Real(kind=cp), Dimension(:,:),   Allocatable :: xyz
  Real(kind=cp), Dimension(:,:,:), Allocatable :: rho
  
  Character(len=1)                             :: end_line=Char(10)
  Character(len=2)                             :: MC_method
  Character(len=79)                            :: title
  Character(len=256)                           :: filnam, filcod, arg2, arg3, arg4

  Logical                                      :: newbond, allvisited, pbc
  Logical, Dimension(3)                        :: percolation
  Logical, Dimension(:), Allocatable           :: newnode

  Write(unit=*,fmt="(/,/,7(a,/))")                                &
       "            ============================="           , &
       "            ======== BVEL_PERCOL ========"           , &
       "            ============================="           , &
       "       ***************************************     " , &
       "       * BVEL Percolation from *.pgrid files *     " , &
       "       ***************************************     " , &
       "          (JRC - ILL, version: January 2016 )"
   
  Write(unit=*,fmt=*) " "

  ! Arguments on the command line

  narg = Command_Argument_Count()

  If (narg > 0) Then
     Call GET_COMMAND_ARGUMENT(1,filnam)
     Call GET_COMMAND_ARGUMENT(2,arg2)
     Call GET_COMMAND_ARGUMENT(3,arg3)
     Call GET_COMMAND_ARGUMENT(4,arg4)

     Read(arg2,*) Eini
     Read(arg3,*) Efin
     Read(arg4,*) dE

     i=Index(filnam,".pgrid")
     If(i /= 0) filcod=filnam(1:i-1)
  Else
     Write(unit=*,fmt="(a)") 'Error! A pgrid file must be given on the command line'
     Stop
  End If

  ! -----------------
  ! Read the BVEL map 
  ! -----------------

  Open(11, file = filnam, &
       form   = 'UNFORMATTED', &
       access = 'STREAM',      &
       action = 'READ')

  Read(11) version
  Read(11) title(1:79),end_line
  Read(11) gType
  Read(11) fType !=0                                                                                                             
  Read(11) nval  !=1                                                                                                             
  Read(11) ndim  !=3                                                                                                             
  Read(11) ngrid
  Read(11) nASYM !                                                                                                               
  Read(11) Cell%cell,Cell%ang
  
  Allocate(rho(ngrid(1),ngrid(2),ngrid(3)))

  Read(11) rho
  
  Close(11)

  Emin = Minval(rho)

  Write(unit=*,fmt="(a,2x,es20.10,/)") "    Emin:", Emin

  rho(:,:,:)  = rho(:,:,:) - Emin

  Emin = 0.

  ! -------------------------------------------------------
  ! Set parameters and arrays for the marching cube routine
  ! and the percolation analysis
  ! -------------------------------------------------------

  ! Parameters for the marching cube method
  
  nlevel    = 1
  MC_method = "TR"
  np_max    = ngrid(1) * ngrid(2) * ngrid(3) * 12 / 8

  Allocate(npoints(nlevel),levels(nlevel),xyz(4,np_max))

  levels(1)      = Eini
  percolation(:) = .False.

  Open(unit=lun,file=Trim(filcod)//".percol", status="replace",action="write")
  Write(unit=lun,fmt="(/,/,7(a,/))")                              &
       "             ==========================="           , &
       "             ======= BVEL PERCOL ======="           , &
       "             ==========================="           , &
       "         ***********************************" , &
       "         *  Percolation from*.pgrid files  *" , &
       "         ***********************************" , &
       "          (JRC - ILL, version: January 2016)"
  Write(unit=lun,fmt="(a,/)") " "

  Do 

     ! ----------------------
     ! Compute the isosurface
     ! ----------------------

     Write(unit=*,fmt="(4x,a,f6.2,a3,/)")  "Energy above the minimum = ", levels(1) - emin, " eV"
     Write(unit=lun,fmt="(4x,a,f6.2,a3,/)")  "Energy above the minimum = ", levels(1) - emin, " eV"

     Write(unit=*,fmt="(8x,a,/)") 'Computing Isosurface....'
     Write(unit=lun,fmt="(8x,a,/)") 'Computing Isosurface....'

     Call Calculate_Mesh(rho,ngrid,nlevel,levels,MC_Method,npoints,xyz)

     If (npoints(1) == Max_points) Then
        Open(unit=err,file=Trim(filcod)//".err", status="replace",action="write")
        Write(unit=err,fmt='(a)') "Error! npoints = Max_Points. Increase Max_Points in CFML_Maps.f90"
        Close(unit=err)
     End If

     ! ---------------------------------------------
     ! Extract graph nodes from xyz
     !
     !    1. Determine how many nodes have the graph
     !    2. Store the nodes
     !    3. Store conectivity (bonds)
     ! ---------------------------------------------

     Write(unit=*,fmt="(8x,a,/)") 'Extracting Graph from the Isosurface....'
     Write(unit=lun,fmt="(8x,a,/)") 'Extracting Graph from the Isosurface....'

     Allocate(id(npoints(1)),newnode(npoints(1)))

     !  1.)  

     Write(unit=*,fmt="(10x,a)") 'Computing nodes....'
     Write(unit=lun,fmt="(10x,a)") 'Computing nodes....'

     Do i = 1 , 3
        !  [1. / real(ngrid(i))] is a measure of the density of points along dimension i
        !  drmin will be used to determine when two nodes are equivalent by translation symmetry
        drmin(i) = 0.01 / Real(ngrid(i))
     End Do

     Do i = 1 , npoints(1)

        ! Put every point inside the unit cell

        Do j = 1 , 3
           If (xyz(j,i) < 0.) Then
              xyz(j,i) = xyz(j,i) + 1.
           Else If (xyz(j,i) > 1.) Then
              xyz(j,i) = xyz(j,i) - 1.
           End If
        End Do

     End Do

     nnodes(:)  = 0

     Do i = 1 , npoints(1)
        newnode(i) = .True.

        Do j = i - 1 , 1 , -1 ! This is the most time consuming part
                              ! Running the loop in this way is much faster
           If (Abs(xyz(1,i)-xyz(1,j)) < drmin(1) .And. &
               Abs(xyz(2,i)-xyz(2,j)) < drmin(2) .And. &
               Abs(xyz(3,i)-xyz(3,j)) < drmin(3)) Then
              newnode(i) = .False.
              Exit
           End If
        End Do
        
        If (newnode(i)) Then 
           nnodes(4) = nnodes(4) + 1
           id(i)  = nnodes(4)
        Else
           id(i)  = id(j)
        End If

     End Do

     ! 2.)

     Write(unit=*,fmt="(10x,a)") 'Storing nodes....'
     Write(unit=lun,fmt="(10x,a)") 'Storing nodes....'

     Allocate(nodes(nnodes(4)))

     Do i = 1 , npoints(1)

        If (newnode(i)) Then 
           nodes(id(i))%xyz(:)     = xyz(1:3,i)
           nodes(id(i))%state      = 0
           nodes(id(i))%cc         = 0
           nodes(id(i))%cc_pbc     = 0
           nodes(id(i))%nbonds     = 0
           nodes(id(i))%nbonds_pbc = 0
           nodes(id(i))%face(:)    = 0

           ! Classify the node: 1 (interior), -1 (face), -2 (edge), -3 (vertex)
        
           nodes(id(i))%type = 0

           Do j = 1 , 3

              If (Abs(nodes(id(i))%xyz(j) - 0.0d0) < 1e-6) Then 
                 nodes(id(i))%type    = nodes(id(i))%type + 1
                 nodes(id(i))%face(j) = -1
              Else If (Abs(nodes(id(i))%xyz(j) - 1.0d0) < 1e-6) Then
                 nodes(id(i))%type    = nodes(id(i))%type + 1
                 nodes(id(i))%face(j) =  1
              End If

           End Do

           nnodes(nodes(id(i))%type) = nnodes(nodes(id(i))%type) + 1
        End If

     End Do

     ! 3.)

     Write(unit=*,fmt="(10x,a,/)") 'Computing bonds....'
     Write(unit=lun,fmt="(10x,a,/)") 'Computing bonds....'

     i          = 1
     nbonds_in  = 0
     nbonds_out = 0
  
     Do j = 1 , npoints(1) / 3 ! run over each triangle stored in xyz
     
        ! tricon stores the triangle conectivity
     
        !        i ------ i+1
        !          \    /
        !           \  /
        !            \/
        !             i+2

        tricon(1,0) = i + 1
        tricon(2,0) = i + 2
        tricon(1,1) = i
        tricon(2,1) = i + 2
        tricon(1,2) = i
        tricon(2,2) = i + 1
     
        Do k = 0 , 2 ! run over the three vertex of the triangle
           v = i + k
        
           Do l = 1 , 2 ! store the two bonds if they haven't previously be stored
              newbond = .True.

              Do m = 1 , nodes(id(v))%nbonds

                 If (nodes(id(v))%bond(m) == id(tricon(l,k))) Then 
                    newbond = .False.
                    Exit
                 End If

              End Do

              If (newbond) Then

                 Do m = 1 , nodes(id(v))%nbonds_pbc

                    If (nodes(id(v))%bond_pbc(m) == id(tricon(l,k))) Then 
                       newbond = .False.
                       Exit
                    End If

                 End Do

              End If

              If (newbond) Then
                 pbc = .False.
              
                 Do m = 1 , 3

                    If (Abs(nodes(id(v))%xyz(m) - nodes(id(tricon(l,k)))%xyz(m)) > 0.5) Then
                       pbc = .True.
                       Exit
                    End If

                 End Do

                 If (.Not. pbc) Then
                    nbonds_in                              = nbonds_in + 1
                    nodes(id(v))%nbonds                    = nodes(id(v))%nbonds + 1
                    nodes(id(v))%bond(nodes(id(v))%nbonds) = id(tricon(l,k))
                 Else
                    nbonds_out                                     = nbonds_out + 1
                    nodes(id(v))%nbonds_pbc                        = nodes(id(v))%nbonds_pbc + 1
                    nodes(id(v))%bond_pbc(nodes(id(v))%nbonds_pbc) = id(tricon(l,k))
                 End If

              End If
           
           End Do
           
        End Do

        i = i + 3
     End Do

     Deallocate(id)
     Deallocate(newnode)

     Write(unit=*,fmt="(12x,a,i6)")   "Total    nodes: ", nnodes(4)
     Write(unit=*,fmt="(12x,a,i6)")   "Interior nodes: ", nnodes(0)
     Write(unit=*,fmt="(12x,a,i6)")   "Face     nodes: ", nnodes(1)
     Write(unit=*,fmt="(12x,a,i6)")   "Edge     nodes: ", nnodes(2)
     Write(unit=*,fmt="(12x,a,i6)")   "Vertex   nodes: ", nnodes(3)
     Write(unit=*,fmt="(12x,a,i6)")   "Interior bonds: ", nbonds_in
     Write(unit=*,fmt="(12x,a,i6,/)") "Exterior bonds: ", nbonds_out

     Write(unit=lun,fmt="(12x,a,i6)")   "Total    nodes: ", nnodes(4)
     Write(unit=lun,fmt="(12x,a,i6)")   "Interior nodes: ", nnodes(0)
     Write(unit=lun,fmt="(12x,a,i6)")   "Face     nodes: ", nnodes(1)
     Write(unit=lun,fmt="(12x,a,i6)")   "Edge     nodes: ", nnodes(2)
     Write(unit=lun,fmt="(12x,a,i6)")   "Vertex   nodes: ", nnodes(3)
     Write(unit=lun,fmt="(12x,a,i6)")   "Interior bonds: ", nbonds_in
     Write(unit=lun,fmt="(12x,a,i6,/)") "Exterior bonds: ", nbonds_out

     ! ----------------------------
     ! Compute connected components
     ! ----------------------------

     Write(unit=*,fmt="(8x,a,/)") 'Computing connected components (Depth First Search algorithm)....'
     Write(unit=lun,fmt="(8x,a,/)") 'Computing connected components (Depth First Search algorithm)....'
  
     m          = 0
     ncc        = 0

     Do 
        allvisited = .True.
        m          = m + 1

        Do i = m , nnodes(4)
           
           If (nodes(i)%state == 0) Then
              allvisited      = .False.
              ncc             = ncc + 1
              nodes(i)%state  = 1
              nodes(i)%cc     = ncc

              Do j = 1 , nodes(i)%nbonds
                 k = nodes(i)%bond(j)
                 If (nodes(k)%state == 0) Call DFS_N(k,ncc)
              End Do

           End If

        End Do

        If (allvisited) Exit
     End Do

     Write(unit=*,fmt="(12x,a,i6)")   "Number of connected components found before applying boundary conditions: ", ncc
     Write(unit=lun,fmt="(12x,a,i6)")   "Number of connected components found before applying boundary conditions: ", ncc

     ! Build the graph for components

     Allocate(cts(ncc))

     Do i = 1 , ncc
        Allocate(cts(i)%bond(ncc))
        cts(i)%cc     = 0
        cts(i)%pbc    = 0
        cts(i)%state  = 0
        cts(i)%nbonds = 0
     End Do

     Do i  = 1 , nnodes(4)
        !Write(1000+nodes(i)%cc,*) nodes(i)%xyz
        ic = nodes(i)%cc
        Do j = 1 , nodes(i)%nbonds_pbc
           k  = nodes(i)%bond_pbc(j)
           kc = nodes(k)%cc
           If (ic /= kc) Then
              newbond = .True.
              Do l = 1 , cts(ic)%nbonds
                 If (cts(ic)%bond(l) == kc) newbond = .False.
              End Do
              If (newbond) Then
                 cts(ic)%nbonds               = cts(ic)%nbonds + 1
                 cts(ic)%bond(cts(ic)%nbonds) = kc
              End If
           End If
        End Do
     End Do

     m       = 0
     ncc_pbc = 0

     Do 
        allvisited = .True.
        m          = m + 1

        Do i = m , ncc
           
           If (cts(i)%state == 0) Then
              allvisited      = .False.
              ncc_pbc         = ncc_pbc + 1
              cts(i)%state  = 1
              cts(i)%cc     = ncc_pbc

              Do j = 1 , cts(i)%nbonds
                 k = cts(i)%bond(j)
                 If (cts(k)%state == 0) Call DFS_C(k,ncc_pbc)
              End Do

           End If

        End Do

        If (allvisited) Exit
     End Do
     
     Write(unit=*,fmt="(12x,a,i6,/)")   "Number of connected components found after  applying boundary conditions: ", ncc_pbc
     Write(unit=lun,fmt="(12x,a,i6,/)")   "Number of connected components found after  applying boundary conditions: ", ncc_pbc

     ! For each component_pbc, we must select one component as origin, that is, the component
     ! for which we do not apply pbc

     Do ic = 1 , ncc_pbc
        jc = 1
        Do
           If (cts(jc)%cc == ic) Exit
           jc = jc + 1
        End Do
        cts(jc)%pbc = 1
     End Do
     
     ! Apply boundary conditions

     Do i = 1 , ncc_pbc
        ic = 1
        Do 
           If (cts(ic)%cc == i .And. cts(ic)%pbc == 0) Then
              Do j = 1 , nnodes(4)
                 If (nodes(j)%cc == ic .And. nodes(j)%nbonds_pbc > 0) Then
                    Do l = 1 , nodes(j)%nbonds_pbc
                       k  = nodes(j)%bond_pbc(l)
                       kc = nodes(k)%cc
                       If (cts(kc)%pbc == 1) Then
                          Do v = 1 , 3
                             dr = nodes(j)%xyz(v) - nodes(k)%xyz(v)
                             If (dr > 0.5) Then
                                t(v) = -1. 
                             Else If (dr < -0.5) Then
                                t(v) =  1. 
                             Else 
                                t(v) =  0.
                             End If
                          End Do
                          cts(ic)%pbc = 1
                          Do m = 1 , nnodes(4)
                             If (nodes(m)%cc == ic) &
                                  nodes(m)%xyz(:) = nodes(m)%xyz(:) + t(:) 
                          End Do
                       End If
                    End Do
                 End If
              End Do
           End If
           pbc = .True.
           Do jc = 1 , ncc
              If (cts(jc)%cc == i .And. cts(jc)%pbc == 0) pbc = .False.
           End Do
           If (pbc) Exit
           ic = ic + 1
           If (ic > ncc) ic = 1
        End Do
     End Do
     
     ! --------------------
     ! Percolation analysis
     ! --------------------

     Write(unit=*,fmt="(8x,a,/)") "Percolation analysis...."
     Write(unit=lun,fmt="(8x,a,/)") "Percolation analysis...."
     Write(*,*) nnodes(4)

     !Do j = 1 , nnodes(4)
     !   Write(2000+cts(nodes(j)%cc)%cc,*) nodes(j)%xyz
     !   Write(3000+nodes(j)%cc,*) nodes(j)%xyz
     !End Do

     Do i = 1 , ncc_pbc
        Write(unit=*,fmt="(12x,a,i3)") "Analysing component ", i
        Write(unit=lun,fmt="(12x,a,i3)") "Analysing component ", i
     
        Do j = 1 , nnodes(4)
           If (cts(nodes(j)%cc)%cc == i .And. nodes(j)%nbonds_pbc > 0) Then
           
              Do k = 1 , nodes(j)%nbonds_pbc
                 l = nodes(j)%bond_pbc(k)
                 
                 Do v = 1 , 3
                    dr = nodes(j)%xyz(v) - nodes(l)%xyz(v)
                    If (Abs(dr) > 0.5 .And. .Not. percolation(v)) Then 
                       percolation(v) = .True.
                       E_percol(v)    = levels(1) - emin
                    End If
                 End Do
              
              End Do

           End If

        End Do

     End Do

     Write(unit=*,fmt="()")
     Write(unit=lun,fmt="()")

     If (percolation(1)) Then
        Write(unit=*,fmt="(12x,a,f6.2,a3)") "Percolation along a: Yes, Activation energy: ", E_percol(1), " eV"
        Write(unit=lun,fmt="(12x,a,f6.2,a3)") "Percolation along a: Yes, Activation energy: ", E_percol(1), " eV"
     Else
        Write(unit=*,fmt="(12x,a)") "Percolation along a: No"
        Write(unit=lun,fmt="(12x,a)") "Percolation along a: No"
     End If
     
     If (percolation(2)) Then
        Write(unit=*,fmt="(12x,a,f6.2,a3)") "Percolation along b: Yes, Activation energy: ", E_percol(2), " eV"
        Write(unit=lun,fmt="(12x,a,f6.2,a3)") "Percolation along b: Yes, Activation energy: ", E_percol(2), " eV"
     Else
        Write(unit=*,fmt="(12x,a)") "Percolation along b: No"
        Write(unit=lun,fmt="(12x,a)") "Percolation along b: No"
     End If
     
     If (percolation(3)) Then
        Write(unit=*,fmt="(12x,a,f6.2,a3/)") "Percolation along c: Yes, Activation energy: ", E_percol(3), " eV"
        Write(unit=lun,fmt="(12x,a,f6.2,a3/)") "Percolation along c: Yes, Activation energy: ", E_percol(3), " eV"
     Else
        Write(unit=*,fmt="(12x,a,/)") "Percolation along c: No"
        Write(unit=lun,fmt="(12x,a,/)") "Percolation along c: No"
     End If

     levels(1) = levels(1) + dE

     If (dE > 0) Then
        If ((percolation(1) .And. percolation(2) .And. percolation(3)) .Or. levels(1) > efin) Exit
     Else
        If ((percolation(1) .And. percolation(2) .And. percolation(3)) .Or. levels(1) < efin) Exit
     End If

     Deallocate(nodes)
     Deallocate(cts)

  End Do

  Close(lun)

End Program BVEL_Percol

Recursive Subroutine DFS_N(id,ncc)

  Use Common

  Implicit None

  Integer, Intent(in) :: id, ncc
  Integer             :: j, k

  nodes(id)%state  = 1
  nodes(id)%cc     = ncc

  Do j = 1 , nodes(id)%nbonds
     k = nodes(id)%bond(j)
     If (nodes(k)%state == 0) Call DFS_N(k,ncc)
  End Do

End Subroutine DFS_N

Recursive Subroutine DFS_C(id,ncc)

  Use Common

  Implicit None

  Integer, Intent(in) :: id, ncc
  Integer             :: j, k

  cts(id)%state  = 1
  cts(id)%cc     = ncc

  Do j = 1 , cts(id)%nbonds
     k = cts(id)%bond(j)
     If (cts(k)%state == 0) Call DFS_C(k,ncc)
  End Do

End Subroutine DFS_C
