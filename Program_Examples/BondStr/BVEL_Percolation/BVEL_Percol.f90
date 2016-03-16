Module Common

  Use CFML_GlobalDeps,       Only:  cp

  Type :: node
     Integer                     :: nbonds, nbonds_pbc    
     Integer                     :: type                ! 1: interior, -1: face, -2: edge, -3: vertex
     Integer                     :: state               ! 0: non-visited, 1: visited 
     Integer                     :: cc                  ! connected component to which the node belongs
     Integer, Dimension(3)       :: face
     Integer, Dimension(12)      :: bond, bond_pbc
     Real(kind=cp), Dimension(3) :: xyz, xyz_aux
  End Type node

  Type(node), Dimension(:), Allocatable        :: nodes
     
End Module Common

Program BVEL_Percol

  !---- Use Modules ----!
  Use CFML_GlobalDeps
  Use CFML_Crystal_Metrics,                 Only: Crystal_Cell_Type
  Use CFML_Maps_Calculations,               Only: Calculate_Mesh
  Use Common

  !---- Variables ----!
  Implicit None

  Type (Crystal_Cell_Type)                     :: Cell

  Integer                                      :: i, j, k, l, m, n, v, cc_aux
  Integer                                      :: narg, nlevel, np_max, ncc, nbonds_aux, nbonds_in, nbonds_out, ncc_pbc
  Integer(kind=4)                              :: gType, fType, nval, ndim, nASYM
  Integer(kind=4), Dimension(3)                :: ngrid
  Integer(kind=4), Dimension(4)                :: version
  Integer,         Dimension(0:4)              :: nnodes
  Integer,         Dimension(12)               :: bond_aux
  Integer,         Dimension(2,0:2)            :: tricon
  Integer,         Dimension(:), Allocatable   :: npoints, id

  Real(kind=cp)                                :: Emin, Emax, dr
  Real(kind=cp), Parameter                     :: dE = 0.01, dE_max = 1.0
  Real(kind=cp), Dimension(3)                  :: drmin, t, E_percol
  Real(kind=cp), Dimension(:),     Allocatable :: levels
  Real(kind=cp), Dimension(:,:),   Allocatable :: xyz
  Real(kind=cp), Dimension(:,:,:), Allocatable :: rho
  
  Character(len=1)                             :: end_line=Char(10)
  Character(len=2)                             :: MC_method
  Character(len=79)                            :: title
  Character(len=256)                           :: filnam, filcod

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
     i=Index(filnam,".pgrid")
     If(i /= 0) filcod=filcod(1:i-1)
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

  emin = Minval(rho)
  emax = Maxval(rho)

  ! -------------------------------------------------------
  ! Set parameters and arrays for the marching cube routine
  ! and the percolation analysis
  ! -------------------------------------------------------

  ! Parameters for the marching cube method
  
  nlevel    = 1
  MC_method = "TR"
  np_max    = ngrid(1) * ngrid(2) * ngrid(3) * 12 / 8

  Allocate(npoints(nlevel),levels(nlevel),xyz(4,np_max))

  levels(1)      = emin
  percolation(:) = .False.
  
  Do 

     ! ----------------------
     ! Compute the isosurface
     ! ----------------------

     Write(unit=*,fmt="(4x,a,f6.2,a3,/)")  "Energy above the minimum = ", levels(1) - emin, " eV"

     Write(unit=*,fmt="(8x,a,/)") 'Computing Isosurface....'

     Call Calculate_Mesh(rho,ngrid,nlevel,levels,MC_Method,npoints,xyz)

     ! ---------------------------------------------
     ! Extract graph nodes from xyz
     !
     !    1. Determine how many nodes have the graph
     !    2. Store the nodes
     !    3. Store conectivity (bonds)
     ! ---------------------------------------------

     Write(unit=*,fmt="(8x,a,/)") 'Extracting Graph from the Isosurface....'

     Allocate(id(npoints(1)),newnode(npoints(1)))

     !  1.)  

     Write(unit=*,fmt="(10x,a)") 'Computing nodes....'

     Do i = 1 , 3
        !  [1. / real(ngrid(i))] is a measure of the density of points along dimension i
        !  drmin will be used to determine when two nodes are equivalent by translation symmetry
        drmin(i) = 0.01 / Real(ngrid(i))
     End Do

     Do i = 1 , npoints(1)

        ! Put every point inside the unit cell

        Do j = 1 , 3
           If (xyz(j,i) .Lt. 0.) Then
              xyz(j,i) = xyz(j,i) + 1.
           Else If (xyz(j,i) .Gt. 1.) Then
              xyz(j,i) = xyz(j,i) - 1.
           End If
        End Do

     End Do

     nnodes(:)  = 0

     Do i = 1 , npoints(1)
        newnode(i) = .True.

        Do j = i - 1 , 1 , -1 ! This is the most time consuming part
                              ! Running the loop in this way is much faster
           If (Abs(xyz(1,i)-xyz(1,j)) .Lt. drmin(1) .And. &
               Abs(xyz(2,i)-xyz(2,j)) .Lt. drmin(2) .And. &
               Abs(xyz(3,i)-xyz(3,j)) .Lt. drmin(3)) Then
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

     Allocate(nodes(nnodes(4)))

     Do i = 1 , npoints(1)

        If (newnode(i)) Then 
           nodes(id(i))%xyz(:)     = xyz(1:3,i)
           nodes(id(i))%state      = 0
           nodes(id(i))%cc         = 0
           nodes(id(i))%nbonds     = 0
           nodes(id(i))%nbonds_pbc = 0
           nodes(id(i))%face(:)    = 0

           ! Classify the node: 1 (interior), -1 (face), -2 (edge), -3 (vertex)
        
           nodes(id(i))%type = 0

           Do j = 1 , 3

              If (Abs(nodes(id(i))%xyz(j) - 0.0d0) .Lt. 1e-6) Then 
                 nodes(id(i))%type    = nodes(id(i))%type + 1
                 nodes(id(i))%face(j) = -1
              Else If (Abs(nodes(id(i))%xyz(j) - 1.0d0) .Lt. 1e-6) Then
                 nodes(id(i))%type    = nodes(id(i))%type + 1
                 nodes(id(i))%face(j) =  1
              End If

           End Do

           nnodes(nodes(id(i))%type) = nnodes(nodes(id(i))%type) + 1
        End If

     End Do

     ! 3.)

     Write(unit=*,fmt="(10x,a,/)") 'Computing bonds....'

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

     ! ----------------------------
     ! Compute connected components
     ! ----------------------------

     Write(unit=*,fmt="(8x,a,/)") 'Computing connected components (Depth First Search algorithm)....'
  
     m          = 0
     ncc        = 0

     Do 
        allvisited = .True.
        m          = m + 1

        Do i = m , nnodes(4)
           
           If (nodes(i)%state .Eq. 0) Then
              allvisited     = .False.
              ncc            = ncc + 1
              nodes(i)%state = 1
              nodes(i)%cc    = ncc

              Do j = 1 , nodes(i)%nbonds
                 k = nodes(i)%bond(j)
                 If (nodes(k)%state == 0) Call DFS(k,ncc)
              End Do

           End If

        End Do

        If (allvisited) Exit
     End Do

     Write(unit=*,fmt="(12x,a,i6)")   "Number of connected components found before applying boundary conditions: ", ncc

     ! Apply periodic boundary conditions
     ! Here we connect components through the boundaries of the unit cell

     ncc_pbc = 0
  
     Do i = 1 , ncc
        n = 0

        Do j = 1 , nnodes(4)
        
           If (nodes(j)%cc == i) Then
              n           = n + 1
              If (n .Eq. 1) ncc_pbc = ncc_pbc + 1
              nodes(j)%cc = ncc_pbc
           
              Do k = j + 1 , nnodes(4)

                 If (nodes(k)%cc == i) Then 
                    n           = n + 1
                    nodes(k)%cc = ncc_pbc
                 End If

              End Do

              Do k = j , nnodes(4)

                 If (nodes(k)%cc == ncc_pbc) Then 
              
                    Do l = 1 , nodes(k)%nbonds_pbc
                       m = nodes(k)%bond_pbc(l)
                       
                       If (nodes(m)%cc /= ncc_pbc) Then
                          
                          Do v = 1 , 3
                             dr = nodes(m)%xyz(v) - nodes(k)%xyz(v)
                             
                             If (dr .Gt. 0.5) Then
                                t(v) = -1.
                             Else If (dr .Lt. -0.5) Then
                                t(v) =  1.
                             Else 
                                t(v) =  0.
                             End If

                          End Do
                             
                          cc_aux = nodes(m)%cc

                          Do v = 1 , nnodes(4)

                             If (nodes(v)%cc == cc_aux) Then
                                n = n + 1
                                nodes(v)%cc     = ncc_pbc
                                nodes(v)%xyz(:) = nodes(v)%xyz(:) + t(:) 
                             End If

                          End Do

                       End If
                       
                    End Do
                    
                 End If

              End Do

           End If

        End Do

     End Do

     Write(unit=*,fmt="(12x,a,i6,/)")   "Number of connected components found after  applying boundary conditions: ", ncc_pbc
     
     ! --------------------
     ! Percolation analysis
     ! --------------------

     Write(unit=*,fmt="(8x,a,/)") "Percolation analysis...."

     Do i = 1 , ncc_pbc
        Write(unit=*,fmt="(12x,a,i3)") "Analysing component ", i
     
        Do j = 1 , nnodes(4)

           If (nodes(j)%cc == i .And. nodes(j)%nbonds_pbc > 0) Then
           
              Do k = 1 , nodes(j)%nbonds_pbc
                 l = nodes(j)%bond_pbc(k)
                 
                 Do v = 1 , 3
                    dr = nodes(j)%xyz(v) - nodes(l)%xyz(v)
                    If (Abs(dr) .Gt. 0.5 .And. .Not. percolation(v)) Then 
                       percolation(v) = .True.
                       E_percol(v)    = levels(1) - emin
                       Emax           = levels(1) + dE_max
                    End If
                 End Do
              
              End Do

           End If

        End Do

     End Do

     Write(unit=*,fmt="()")

     If (percolation(1)) Then
        Write(unit=*,fmt="(12x,a,f6.2,a3)") "Percolation through bc: Yes, Activation energy: ", E_percol(1), " eV"
     Else
        Write(unit=*,fmt="(12x,a)") "Percolation through bc: No"
     End If
     
     If (percolation(2)) Then
        Write(unit=*,fmt="(12x,a,f6.2,a3)") "Percolation through ac: Yes, Activation energy: ", E_percol(2), " eV"
     Else
        Write(unit=*,fmt="(12x,a)") "Percolation through ac: No"
     End If
     
     If (percolation(3)) Then
        Write(unit=*,fmt="(12x,a,f6.2,a3/)") "Percolation through ab: Yes, Activation energy: ", E_percol(3), " eV"
     Else
        Write(unit=*,fmt="(12x,a,/)") "Percolation through ab: No"
     End If

     levels(1) = levels(1) + dE
     If ((percolation(1) .And. percolation(2) .And. percolation(3)) .Or. levels(1) > emax) Exit

     Deallocate(nodes)

  End Do

End Program BVEL_Percol

Recursive Subroutine DFS(id,ncc)

  Use Common

  Implicit None

  Integer, Intent(in) :: id, ncc
  Integer             :: j, k

  nodes(id)%state = 1
  nodes(id)%cc    = ncc

  Do j = 1 , nodes(id)%nbonds
     k = nodes(id)%bond(j)
     If (nodes(k)%state == 0) Call DFS(k,ncc)
  End Do

End Subroutine DFS
