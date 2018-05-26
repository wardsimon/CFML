Submodule (CFML_Crystal_Metrics) GenMetrics

 Contains
    !!----
    !!---- FUNCTION SIGMAV_CELLTYPE
    !!----
    !!----    Calculates the standard deviation of the unit cell volume
    !!----    from the standard deviations of cell parameters. The input
    !!----    variable is of type Crytal_Cell_Type, if the standard deviations of
    !!----    of both cell axes and cell angles are zero the result is sigma=0.0,
    !!----    otherwise the calculation is performed
    !!----    It is assumed that there is no correlation (covariance terms) between
    !!----    the standard deviations of the different cell parameters.
    !!----
    !!---- Updated: January - 2013 (JRC)
    !!
    Module Pure Function SigmaV_CellType(Cell) Result(sigma)
       !---- Arguments ----!
       class(CrysCell_Type), intent(in) :: Cell      ! Cell Parameters
       real(kind=cp)                    :: sigma     ! Sigma

       !--- Local variables ---!
       integer                     :: i
       real(kind=cp)               :: cr1, cr2 
       real(kind=cp)               :: q,ca,cb,cc,vc,sa,sb,sc
       real(kind=cp), dimension(3) :: var_ang

       !> Init
       sigma=0.0_cp
       cr1=sum(abs(cell%scell))
       cr2=sum(abs(cell%sang))
       if (cr1 < eps .and. cr2 < eps ) return

       vc=0.0_cp
       do i=1,3
          q=cell%scell(i)/cell%cell(i)
          vc=vc+q*q
       end do
       
       if (cr2 > eps) then
          ca=cosd(cell%ang(1)) ;  sa=sind(cell%ang(1))
          cb=cosd(cell%ang(2)) ;  sb=sind(cell%ang(2))
          cc=cosd(cell%ang(3)) ;  sc=sind(cell%ang(3))
          q=1.0-ca*ca-cb*cb-cc*cc+2.0*ca*cb*cc
          
          var_ang = (cell%sang * TO_RAD)**2/q
          vc=vc+ (ca-cb*cc)*(ca-cb*cc)*sa*sa * var_ang(1)
          vc=vc+ (cb-ca*cc)*(cb-ca*cc)*sb*sb * var_ang(2)
          vc=vc+ (cc-ca*cb)*(cc-ca*cb)*sc*sc * var_ang(3)
       end if

       sigma=cell%vol*sqrt(vc)

       return
    End Function SigmaV_CellType
    
    !!----
    !!---- FUNCTION SIGMAV_CELLTYPE
    !!----
    !!----    Calculates the standard deviation of the unit cell volume
    !!----    from the standard deviations of cell parameters. The input
    !!----    variable is of type Crytal_Cell_Type, if the standard deviations of
    !!----    of both cell axes and cell angles are zero the result is sigma=0.0,
    !!----    otherwise the calculation is performed
    !!----    It is assumed that there is no correlation (covariance terms) between
    !!----    the standard deviations of the different cell parameters.
    !!----
    !!---- Updated: January - 2013 (JRC)
    !!
    Module Pure Function SigmaV_Cell(cell,ang,scell,sang) Result(sigma)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: cell      ! Cell parameters
       real(kind=cp), dimension(3), intent(in) :: ang 
       real(kind=cp), dimension(3), intent(in) :: scell     ! standard deviation for cell parameters
       real(kind=cp), dimension(3), intent(in) :: sang
       real(kind=cp)                           :: sigma     ! Sigma

       !--- Local variables ---!
       integer                     :: i
       real(kind=cp)               :: cr1, cr2 
       real(kind=cp)               :: q,ca,cb,cc,vc,sa,sb,sc
       real(kind=cp), dimension(3) :: var_ang

       !> Init
       sigma=0.0_cp
       cr1=sum(abs(scell))
       cr2=sum(abs(sang))
       if (cr1 < eps .and. cr2 < eps ) return

       vc=0.0_cp
       do i=1,3
          q=scell(i)/cell(i)
          vc=vc+q*q
       end do
       
       if (cr2 > eps) then
          ca=cosd(ang(1)) ;  sa=sind(ang(1))
          cb=cosd(ang(2)) ;  sb=sind(ang(2))
          cc=cosd(ang(3)) ;  sc=sind(ang(3))
          q=1.0-ca*ca-cb*cb-cc*cc+2.0*ca*cb*cc
          
          var_ang = (sang * TO_RAD)**2/q
          vc=vc+ (ca-cb*cc)*(ca-cb*cc)*sa*sa * var_ang(1)
          vc=vc+ (cb-ca*cc)*(cb-ca*cc)*sb*sb * var_ang(2)
          vc=vc+ (cc-ca*cb)*(cc-ca*cb)*sc*sc * var_ang(3)
       end if

       sigma=volume_cell(cell,ang)*sqrt(vc)

       return
    End Function SigmaV_Cell
    
    !!----
    !!---- FUNCTION VOLUME_CELL
    !!----
    !!---- Calculate the volume of a Cell
    !!----
    !!
    Module Pure Function Volume_Cell(cell,ang) Result(Vol)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: cell
       real(kind=cp), dimension(3), intent(in) :: ang
       real(kind=cp)                           :: vol
       
       !---- Local Variables ----!
       integer       :: i
       real(kind=cp) :: p,s,cs
       !> Init
       p=1.0_cp; s=1.0_cp
       
       do i=1,3
          cs=cosd(ang(i))
          p=p*cs
          s=s-cs*cs
       end do
       vol=sqrt(abs(s + 2.0_cp*p))

       do i=1,3
          vol=vol*cell(i)
       end do
       
       return
    End Function Volume_Cell  
    
    !!--++
    !!--++ FUNCTION METRICS
    !!--++
    !!--++    (PRIVATE)
    !!--++    Constructs the metric tensor
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Pure Function Metrics(cell,ang) Result(G)
       !---- Arguments ----!
       real(kind=cp), dimension(3)  , intent(in ) :: cell  ! Cell Parameters
       real(kind=cp), dimension(3)  , intent(in ) :: ang
       real(kind=cp), dimension(3,3)              :: G     ! Metric Tensor

       !---- Local Variables ----!
       integer :: i

       G(1,2)= cell(1)*cell(2)*cosd(ang(3))
       G(1,3)= cell(1)*cell(3)*cosd(ang(2))
       G(2,3)= cell(2)*cell(3)*cosd(ang(1))

       do i=1,3
          G(i,i)= cell(i)*cell(i)
       end do

       G(2,1)=G(1,2)
       G(3,1)=G(1,3)
       G(3,2)=G(2,3)

       return
    End Function Metrics
    
    !!--++
    !!--++ SUBROUTINE GET_CRYST_ORTHOG_MATRIX
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtains the matrix giving the crystallographic basis in
    !!--++    direct space in terms of a Cartesian basis. The output matrix
    !!--++    can be directly used for transforming crystallographic components
    !!--++    to Cartesian components of the components of a vector considered
    !!--++    as a column vector:   XC = CrystOrt X.
    !!--++
    !!--++    If CartType is not present, or if it is not equal to 'A',
    !!--++    the cartesian system is defined as:
    !!--++          z // c; y is in the bc-plane; x is y ^ z
    !!--++    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
    !!--++    b = (         0         ,     b sinalpha      , b cosalpha)
    !!--++    c = (         0         ,         0           , c         )
    !!--++
    !!--++    If CartType = 'A', the Cartesian system is defined as:
    !!--++         x // a; y is in the ab-plane; z is x ^ z
    !!--++    a = (       a   ,         0           ,       0             )
    !!--++    b = ( b cosgamma,    b singamma       ,       0             )
    !!--++    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )
    !!--++
    !!--++    The output matrix is the tranposed of the above one(s) so that the
    !!--++    matrix can directly be used for transforming "components" given
    !!--++    in a crystallographic basis to "components" in cartesian basis
    !!--++    when the components are used as "column" vectors.
    !!--++
    !!--++      [a] = C [e] , In [a],[e] basis vectors are in column form
    !!--++      (a) = (e) CT, In (a),(e) basis vectors are in row form
    !!--++      CrystOrt = CT  => (a) = (e) CystOrt, in ITC: (a) = (e) P
    !!--++
    !!--++    Remember that  C.CT = GD (direct cell metrics)
    !!--++
    !!--++
    !!--++      Xc = CrystOrt X (Xc Cartesian components, X crystallographic components)
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Pure Subroutine Get_Cryst_Orthog_Matrix(Cell,Ang, Mat,CarType)
       !---- Arguments ----!
       real(kind=cp), dimension(3  ), intent (in ) :: cell,ang   ! Cell Parameters
       real(kind=cp), dimension(3,3), intent (out) :: Mat        ! Convsersion matrix
       character(len=*), optional,    intent (in ) :: CarType    ! Type of Cartesian axes

       !---- Local Variables ----!
       character(len=1) :: car
       real(kind=cp)    :: cosgas, singas

       car='C'
       if (present(CarType)) car=adjustl(u_case(cartype))
       
       select case (car)
          case ('A')   ! x//a
             !  Transponse of the following matrix:
             !    a = (       a   ,         0           ,       0             )
             !    b = ( b cosgamma,    b singamma       ,       0             )
             !    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )
             cosgas =(cosd(ang(3))*cosd(ang(2))-cosd(ang(1)))/(sind(ang(3))*sind(ang(2)))
             singas = sqrt(1.0-cosgas**2)
             Mat(1,1) = cell(1)
             Mat(1,2) = cell(2)*cosd(ang(3))
             Mat(1,3) = cell(3)*cosd(ang(2))
             Mat(2,1) = 0.0
             Mat(2,2) = cell(2)*sind(ang(3))
             Mat(2,3) =-cell(3)*sind(ang(2))*cosgas
             Mat(3,1) = 0.0
             Mat(3,2) = 0.0
             Mat(3,3) = cell(3)*sind(ang(2))*singas
             
          case default  
             !  By default, the cartesian frame is such as z//c
             !  Transponse of the following matrix:
             !    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
             !    b = (         0         ,     b sinalpha      , b cosalpha)
             !    c = (         0         ,         0           , c         )
             cosgas =(cosd(ang(1))*cosd(ang(2))-cosd(ang(3)))/(sind(ang(1))*sind(ang(2)))
             singas = sqrt(1.0-cosgas**2)
             Mat(1,1) = cell(1)*sind(ang(2))*singas
             Mat(1,2) = 0.0
             Mat(1,3) = 0.0
             Mat(2,1) =-cell(1)*sind(ang(2))*cosgas
             Mat(2,2) = cell(2)*sind(ang(1))
             Mat(2,3) = 0.0
             Mat(3,1) = cell(1)*cosd(ang(2))
             Mat(3,2) = cell(2)*cosd(ang(1))
             Mat(3,3) = cell(3)
             
       end select 
       
       return
    End Subroutine Get_Cryst_Orthog_Matrix 
    
    !!--++
    !!--++ Subroutine ReciprocalCell
    !!--++
    !!--++    (PRIVATE)
    !!--++    Calculates the reciprocal lattice vectors 
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Pure Subroutine ReciprocalCell(cell,ang,rcell,rang,rVol)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in ) :: cell,ang
       real(kind=cp), dimension(3), intent(out) :: rcell,rang
       real(kind=cp),               intent(out) :: rvol

       !---- Local Variables ----!
       integer        :: i
       real(kind=cp)  :: vol,s,p,cose

       !> Init
       rcell=0.0_cp; rang=0.0_cp; rvol=0.0_cp     ! Initializing return values
       rvol=0.0_cp
       
       vol=Volume_Cell(cell,ang)
       if (vol <= eps .or. err_CFML%state) return
       
       rvol=1.0_cp/vol

       rcell(1)=cell(2)*Cell(3)*sind(ang(1))*rvol
       rcell(2)=cell(3)*Cell(1)*sind(ang(2))*rvol
       rcell(3)=cell(1)*Cell(2)*sind(ang(3))*rvol
       rang(1)=(cosd(ang(2))*cosd(ang(3))-cosd(ang(1)))/(sind(ang(2))*sind(ang(3)))
       rang(2)=(cosd(ang(1))*cosd(ang(3))-cosd(ang(2)))/(sind(ang(1))*sind(ang(3)))
       rang(3)=(cosd(ang(2))*cosd(ang(1))-cosd(ang(3)))/(sind(ang(2))*sind(ang(1)))
       rang=acosd(rang)

       return
    End Subroutine ReciprocalCell
    
    !!----
    !!---- SUBROUTINE SET_CRYSTAL_CELL
    !!----
    !!----    Constructs the object "Celda" of type Crystal_Cell
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine Set_Crystal_Cell(VCell,VAng,Cell,Cartype,Vscell,Vsang)
       !---- Arguments ----!
       real(kind=cp), dimension (3),        intent(in ) :: Vcell, Vang    ! Cell parameters
       class(CrysCell_Type),                intent(out) :: Cell           ! Cell Object
       character (len=*),          optional,intent(in ) :: CarType        ! Orientation in Cartesian
       real(kind=cp), dimension(3),optional,intent(in ) :: Vscell, Vsang  ! Standard deviations

       !---- Local Variables ----!
       integer :: ifail

       !> init
       call clear_error()

       !> a,b,c
       cell%cell=vcell
       if (present(Vscell)) cell%scell=Vscell 
       
       !> angles
       cell%ang=vang
       if (present(Vsang)) cell%scell=Vsang 
       where(cell%ang < eps) cell%ang =90.0_cp
       
       select type (cell)
          class is (Cryscell_M_Type)
             if (present(Cartype)) cell%cartType=adjustl(u_case(cartype)) 
       end select

       !> Volume
       cell%vol=Volume_Cell(Vcell, Vang)
       if (cell%vol <= tiny(0.0)) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          err_CFML%Msg=" Volume of the current cell parameters was zero. Please, check it!"
          return
       end if   
       cell%svol=sigmaVolume(cell)
       
       select type (cell)
          class is (crysCell_M_Type)
             !> Reciprocal 
             call reciprocalcell(Vcell,Vang,cell%rcell,cell%rang,cell%rvol)
             
             !> GD and GR
             cell%GD=metrics(Vcell,Vang)
             cell%GR=metrics(cell%rcell,cell%rang) 
             
             !> Cartesian conversion
             if (cell%cartType /= 'A') cell%cartType='C'
             call Get_Cryst_Orthog_matrix(Vcell,Vang,Cell%Cr_Orth_cel,cell%CartType)
             Cell%Orth_Cr_cel=invert_array3x3(Cell%Cr_Orth_cel)
                  
             !> Busing-Levy matrix component 
             if (cell%CartType == "C") then
                cell%bl_m=Transpose(Cell%Orth_Cr_cel)
                cell%bl_minv=Transpose(Cell%Cr_Orth_cel)
                
             else
                Cell%bl_m(1,1)=cell%rcell(1)
                Cell%bl_m(1,2)=cell%rcell(2)*cosd(cell%rang(3))
                Cell%bl_m(1,3)=cell%rcell(3)*cosd(cell%rang(2))
                Cell%bl_m(2,2)=cell%rcell(2)*sind(cell%rang(3))
                Cell%bl_m(2,3)=-(cell%rcell(3)*sind(cell%rang(2))*cosd(cell%ang(1)))
                Cell%bl_m(3,3)=1.0_cp/cell%cell(3)
                Cell%bl_m(2,1)=0.0_cp
                Cell%bl_m(3,1)=0.0_cp
                Cell%bl_m(3,2)=0.0_cp
                Cell%bl_Minv=invert_array3x3(Cell%bl_M)
             end if
                     
       end select 
       
       return
    End Subroutine Set_Crystal_Cell
    
    !!----
    !!---- SUBROUTINE GET_CRYST_FAMILY
    !!----
    !!----    Obtain the Crystal Family, Symbol and System from cell parameters
    !!----
    !!---- Update: May - 2005
    !!----
    Module Subroutine Get_Cryst_Family(Cell, Str_Family, Str_Symbol, Str_System)
       !---- Arguments ----!
       class(CrysCell_Type),   intent(in ) :: Cell
       character(len=*),       intent(out) :: Str_Family
       character(len=*),       intent(out) :: Str_Symbol
       character(len=*),       intent(out) :: Str_System

       !---- Local variables ----!
       integer, dimension(3) :: icodp, icoda
       integer               :: n1,n2

       !> Init
       Str_Family=" "
       Str_Symbol=" "
       Str_System=" "

       icodp=0; icoda=0

       !> Codification 

       !> a
       icodp(1)=1

       !> b 
       if (abs(cell%cell(2)-cell%cell(1)) <= eps) then
          icodp(2)=icodp(1)
       else
          icodp(2)=2
       end if

       !> c
       if (abs(cell%cell(3)-cell%cell(1)) <= eps) then
          icodp(3)=icodp(1)
       else
          icodp(3)=3
       end if

       !> alpha
       icoda(1)=1

       !> beta
       if (abs(cell%ang(2)-cell%ang(1)) <= eps) then
          icoda(2)=icoda(1)
       else
          icoda(2)=2
       end if

       !>gamma 
       if (abs(cell%ang(3)-cell%ang(1)) <= eps) then
          icoda(3)=icoda(1)
       else
          icoda(3)=3
       end if

       n1=count(icoda==icoda(1))  ! angles
       n2=count(icodp==icodp(1))  ! parameters
       
       select case (n1)
          case (1) ! All angles are differents
             Str_Family="Triclinic"
             Str_Symbol ="a"
             Str_System ="Triclinic"

          case (2) ! two angles are equal
             if (icoda(1) == icoda(2)) then
                if (abs(cell%ang(3)-120.0) <= eps) then
                   if (icodp(1)==icodp(2)) then
                      !> Hexagonal 
                      Str_Family="Hexagonal"
                      Str_Symbol ="h"
                      Str_System ="Hexagonal"
                   end if
                else
                   !> Monoclinic 
                   Str_Family="Monoclinic"
                   Str_Symbol ="m"
                   Str_System ="Monoclinic"
                end if

             else
                !> Monoclic b-unique setting
                if (abs(cell%ang(1)-90.0) <= eps) then
                   Str_Family="Monoclinic"
                   Str_Symbol ="m"
                   Str_System ="Monoclinic"
                end if
             end if

          case (3) ! all angles are equals
             if (abs(cell%ang(1) - 90.000) <= eps) then
                select case (n2)
                   case (1) 
                      !> Orthorhombic 
                      Str_Family="Orthorhombic"
                      Str_Symbol ="o"
                      Str_System ="Orthorhombic"

                   case (2)
                      !> Tetragonal 
                      if (icodp(1)==icodp(2)) then
                         Str_Family="Tetragonal"
                         Str_Symbol ="t"
                         Str_System ="Tetragonal"
                      end if   

                   case (3)
                      !> Cubic
                      Str_Family="Cubic"
                      Str_Symbol ="c"
                      Str_System ="Cubic"
                end select

             else
                if (n2 == 3) then
                   !> Hexagonal with rhombohedral axes 
                   Str_Family="Hexagonal"
                   Str_Symbol ="h"
                   Str_System ="Trigonal"
                end if
             end if

       end select ! n
       
       !> Error?
       if (len_trim(Str_Family) <= 0) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          err_CFML%Msg=" Error obtaining the Crystal Family. Please, check the cell parameters!"
       end if 

       return
    End Subroutine Get_Cryst_Family
    
    !!----
    !!---- SUBROUTINE GET_DERIV_ORTH_CELL
    !!----
    !!----    Subroutine to get derivative matrix of the transformation matrix
    !!----    to orthogonal frame. Useful for calculations of standard deviations
    !!----    of distances and angles. The specialised subroutine calculating
    !!----    sigmas of distances "distance_and_sigma" is in Atom_mod.
    !!----    The output matrices "de_Orthcell" are the derivatives of, with
    !!----    respect to a(1),b(2),c(3),alpha(4),beta(5) and gamma(6) of the
    !!----    matrix   "Cellp%Cr_Orth_cel".
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine Get_Deriv_Orth_Cell(Cell,De_Orthcell,Cartype)
       !---- Arguments ----!
       class(CrysCell_type),            intent(in ) :: cell
       real(kind=cp), dimension(3,3,6), intent(out) :: de_Orthcell
       character (len=1), optional,     intent(in ) :: CarType

       !---- Local Variables ----!
       real(kind=cp)    ::  ca,cb,cg,sa,sb,sg,f,g, fa,fb,fc,ga,gb,gc
       character(len=1) :: car

       !> Init
       de_Orthcell=0.0_cp
       
       car=" "
       if (present(CarType)) car=adjustl(u_case(Cartype))

       ca=cosd(cell%ang(1))
       cb=cosd(cell%ang(2))
       cg=cosd(cell%ang(3))
       sa=sind(cell%ang(1))
       sb=sind(cell%ang(2))
       sg=sind(cell%ang(3))

       select case (car)
          case ('A')
             f=(ca-cb*cg)/sg    !-cosgas*sinbeta
             g=SQRT(sb*sb-f*f)  ! singas*sinbeta
             fa=-sa/sg          ! df/dalpha
             fb=sb*cg/sg        ! df/dbeta
             fc=cb/sg**2        ! df/dgamma
             ga=-f*fa/g         ! dg/dalpha
             gb=(sb*cb-f*fb)/g  ! dg/dbeta
             gc=f/g*fc          ! dg/dgamma

             ! M: Transponse of the following matrix:
             !    a = (       a   ,         0           ,       0             )
             !    b = ( b cosgamma,    b singamma       ,       0             )
             !    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )

             !
             !        (   a         b*cg        c*cb )
             !    M = (   0         b*sg        c*f  )
             !        (   0          0          c*g  )
             !
             !           (   1      0      0 )
             !  dM_da =  (   0      0      0 )
             !           (   0      0      0 )
             de_Orthcell(1,1,1) = 1.0

             !           (   0      cg     0 )
             !  dM_db =  (   0      sg     0 )
             !           (   0      0      0 )
             de_Orthcell(1,2,2) = cg
             de_Orthcell(2,2,2) = sg

             !
             !            (   0          0          cb )
             !  dM_dc =   (   0          0          f  )
             !            (   0          0          g  )
             de_Orthcell(1,3,3) = cb
             de_Orthcell(2,3,3) = f
             de_Orthcell(3,3,3) = g

             !
             !             (   0          0           0   )
             ! dM_dalpha=  (   0          0          c*fa )
             !             (   0          0          c*ga )
             !
             de_Orthcell(2,3,4) = cell%cell(3)*fa
             de_Orthcell(3,3,4) = cell%cell(3)*ga

             !
             !             (   0          0         -c*sb )
             ! dM_dbeta =  (   0          0          c*fb )
             !             (   0          0          c*gb )
             !
             de_Orthcell(1,3,5) = -cell%cell(3)*sb
             de_Orthcell(2,3,5) =  cell%cell(3)*fb
             de_Orthcell(3,3,5) =  cell%cell(3)*gb

             !
             !              (   0        -b*sg         0   )
             ! dM_dgamma =  (   0         b*cg        c*fc )
             !              (   0          0          c*gc )
             !
             de_Orthcell(1,2,6) = -cell%cell(2)*sg
             de_Orthcell(2,2,6) =  cell%cell(2)*cg
             de_Orthcell(2,3,6) =  cell%cell(3)*fc
             de_Orthcell(3,3,6) =  cell%cell(3)*gc
             
          case default 
             !  By default, the cartesian frame is such as z//c
             !  Transponse of the following matrix:
             !    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
             !    b = (         0         ,     b sinalpha      , b cosalpha)
             !    c = (         0         ,         0           , c         )
             
             !         ( a sinbeta singamma*          0             0 )
             !    M =  (-a sinbeta cosgamma*      b sinalpha        0 )
             !         ( a cosbeta                b cosalpha        c )
             
             f=(cg-ca*cb)/sa    !-sinbeta . cosgamma*
             g=SQRT(sb*sb-f*f)  ! sinbeta . singamma*
             fa= cb/sa**2       ! df/dalpha
             fb=sb*ca/sa        ! df/dbeta
             fc=-sb/sa          ! df/dgamma
             ga=-f*fa/g         ! dg/dalpha
             gb=(sb*cb-f*fb)/g  ! dg/dbeta
             gc=f/g*fc          ! dg/dgamma
             
             !         ( a*g        0         0 )
             !    M =  ( a*f      b*sa        0 )
             !         ( a*cb     b*ca        c )
             
             !
             !           (   g       0      0 )
             !  dM_da =  (   f       0      0 )
             !           (   cb      0      0 )
             de_Orthcell(1,1,1) = g
             de_Orthcell(1,2,1) = f
             de_Orthcell(1,3,1) = cb
             
             !           (   0      0      0 )
             !  dM_db =  (   0      sa     0 )
             !           (   0      ca     0 )
             de_Orthcell(1,2,2) = sa
             de_Orthcell(3,2,2) = ca
             
             !
             !            (   0      0      0  )
             !  dM_dc =   (   0      0      0  )
             !            (   0      0      1  )
             de_Orthcell(3,3,3) = 1
             
             !
             !             ( a*ga         0          0 )
             ! dM_dalpha=  ( a*fa       -b*ca        0 )
             !             (   0         b*sa        0 )
             !
             de_Orthcell(1,1,4) = cell%cell(1)*ga
             de_Orthcell(2,1,4) = cell%cell(1)*fa
             de_Orthcell(2,2,4) =-cell%cell(2)*ca
             de_Orthcell(3,2,4) = cell%cell(2)*sa
             
             !
             !             (  a*gb        0         0 )
             ! dM_dbeta =  (  a*fb        0         0 )
             !             ( -a*sb        0         0 )
             !
             de_Orthcell(1,1,5) = cell%cell(1)*gb
             de_Orthcell(2,1,5) = cell%cell(1)*fb
             de_Orthcell(3,1,5) =-cell%cell(1)*sb
             
             !
             !              (  a*gc     0      0   )
             ! dM_dgamma =  (  a*fc     0      0   )
             !              (   0       0      0   )
             !
             de_Orthcell(1,1,6) = cell%cell(1)*gc
             de_Orthcell(2,1,6) = cell%cell(1)*fc 
       end select 
       
       return
    End Subroutine Get_Deriv_Orth_Cell
    
    !!----
    !!---- SUBROUTINE CHANGE_SETTING_CELL
    !!----
    !!---- Calculates a new cell giving the transformation matrix.
    !!---- The input matrix can be given as the S-matrix in International
    !!---- Tables or its transposed (default) that corresponds to the matrix
    !!---- relating formal column matrices containing the basis vectors.
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine Change_Setting_Cell(Cell,Mat,Celln,Matkind)
       !---- Arguments ----!
       class(CrysCell_Type),          intent( in)     :: Cell
       real(kind=cp), dimension (3,3),intent( in)     :: Mat
       class(CrysCell_Type),          intent(out)     :: Celln
       character(len=*),  optional,   intent (in)     :: Matkind

       !--- Local variables ---!
       integer                       :: i
       real(kind=cp), dimension(3)   :: cellv,angl
       real(kind=cp), dimension(3,3) :: S,GN,ST


       if (present(matkind)) then
          if (matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
             S=Mat
            ST=transpose(Mat)
          else
             S=transpose(Mat)
            ST=Mat
          end if
       else
          S=transpose(Mat)
         ST=Mat
       end if

       select type (cell)
          class is (CrysCell_M_Type)
             !> Get the new metric tensor
             !> GDN= Mat GD MatT  or GDN= ST GD S
             GN=matmul(ST,matmul(Cell%GD,S))

             !> Calculate new cell parameters from the new metric tensor
             do i=1,3
                cellv(i)=sqrt(GN(i,i))
             end do
             angl(1)=acosd(GN(2,3)/(cellv(2)*cellv(3)))
             angl(2)=acosd(GN(1,3)/(cellv(1)*cellv(3)))
             angl(3)=acosd(GN(1,2)/(cellv(1)*cellv(2)))

             call Set_Crystal_Cell(cellv,angl,Celln)
          
          class default
             err_CFML%state=.true.
             err_CFML%Flag=2
             err_CFML%Msg="The input cell is an incompatible variable for Change_Setting_Cell procedure."
       end select      

       return
    End Subroutine Change_Setting_Cell
    
    !!----
    !!---- SUBROUTINE GET_PRIMITIVE_CELL
    !!----
    !!----    Subroutine for getting the primitive cell from a centred cell
    !!----    On input Lat_type is the lattice type: P,A,B,C,I,R or F
    !!----    Centred_cell is the Crystal_Cell_Type of the input lattice
    !!----    The subroutine calculates the transformation matric "transfm"
    !!----
    !!---- Update: April - 2008
    !!
    Module Subroutine Get_Primitive_Cell(Lat_Type,Centred_Cell,Primitive_Cell,Transfm)
       !---- Arguments ----!
       character(len=*),              intent(in)  :: lat_type          ! Lattice type
       class(CrysCell_Type),          intent(in)  :: centred_cell      ! Input Cell Object
       class(CrysCell_Type),          intent(out) :: primitive_cell    ! Output Cell Object
       real(kind=cp), dimension(3,3), intent(out) :: transfm           ! Transformation Matrix between Cell objects 

       !---- Local variables ----!
       integer                       :: i
       real(kind=cp), dimension(3)   :: celp,celang
       real(kind=cp), dimension(3,3) :: cart,metric
       character(len=1)              :: lat

       lat=adjustl(lat_type)
       Select Case(lat)
          case("a","A")
             transfm= reshape((/1.0,0.0,0.0,  0.0,0.5,0.5,  0.0,-0.5,0.5/),(/3,3/))
          
          case("b","B")
             transfm= reshape((/0.5,0.0,0.5,  0.0,1.0,0.0, -0.5, 0.0,0.5/),(/3,3/))
          
          case("c","C")
             transfm= reshape((/0.5,0.5,0.0, -0.5,0.5,0.0,  0.0, 0.0,1.0/),(/3,3/))
          
          case("i","I")
             transfm= reshape((/1.0,0.0,0.0,  0.0,1.0,0.0,  0.5, 0.5,0.5/),(/3,3/))
          
          case("r","R")
             transfm= reshape((/2.0/3.0, 1.0/3.0, 1.0/3.0,  &
                               -1.0/3.0, 1.0/3.0, 1.0/3.0,  &
                               -1.0/3.0,-2.0/3.0, 1.0/3.0/),(/3,3/))
          
          case("f","F")
             transfm= reshape((/0.5,0.0,0.5,  0.5,0.5,0.0,  0.0, 0.5,0.5/),(/3,3/))
          
          case default  !assumed primitive
             transfm= reshape((/1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0/),(/3,3/))
       End Select
       
       select type (Centred_Cell)
          class is (CrysCell_M_Type)
             transfm=transpose(transfm)
             cart=matmul(transfm,transpose(Centred_Cell%Cr_Orth_cel))
             metric=matmul(cart,transpose(cart))

             !> Calculate new cell parameters from the new metric tensor
             do i=1,3
                Celp(i)=sqrt(metric(i,i))
             end do

             celang(1)=acosd(metric(2,3)/(celp(2)*celp(3)))
             celang(2)=acosd(metric(1,3)/(celp(1)*celp(3)))
             celang(3)=acosd(metric(1,2)/(celp(1)*celp(2)))
             
          class default
             celp=centred_cell%cell
             celang=centred_cell%ang
       end select         
       call Set_Crystal_Cell(celp,celang,primitive_cell)

       return
    End Subroutine Get_Primitive_Cell
    
    !!----
    !!---- SUBROUTINE GET_TRANSFM_MATRIX
    !!----
    !!----    Subroutine for getting the transformation matrix between two
    !!----    primitive unit cells (the range of indices is fixed to -2 to 2)
    !!----
    !!---- Update: January - 2011
    !!
    Module Subroutine Get_Transfm_Matrix(cella,cellb,trm,tol)
       !---- Arguments ----!
       class(CrysCell_Type),          intent(in) :: cella  ! Cell object
       class(CrysCell_Type),          intent(in) :: cellb  ! Cell object
       real(kind=cp), dimension(3,3), intent(out):: trm    ! Transformation matrix
       real(kind=cp), optional,       intent(in) :: tol    ! Tolerabce

       !---- Local variables ----!
       type(CrysCell_Type)     :: Cellt
       integer,dimension(3,3)  :: Nu
       integer                 :: j,i1,i2,i3,i4,i5,i6,i7,i8,i9
       real(kind=cp)           :: tolt
       Logical                 :: ok

       !> Init
       tolt=0.3
       if(present(tol)) tolt=tol
       Trm=0.0_cp
       
       select type(cella)
          type is (CrysCell_type)
             err_CFML%state=.false.
             err_CFML%Flag=2
             err_CFML%Msg=" The input cell type is incompatible to use on Get_Transfm_Matrix"
             return
             
          type is (CrysCell_LS_type)  
             err_CFML%state=.false.
             err_CFML%Flag=2
             err_CFML%Msg=" The input cell type is incompatible to use on Get_Transfm_Matrix"
             return 
             
       end select      
             
       ok=.false.
       dox: do i1=-2,2                     !         |i1  i4  i7|
          do i2=-2,2                       !    Nu = |i2  i5  i8|
             do i3=-2,2                    !         |i3  i6  i9|
                do i4=-2,2
                   do i5=-2,2
                      do i6=-2,2
                         do i7=-2,2
                            do i8=-2,2
                               do i9=-2,2
                                  j=i1*i5*i9+i4*i8*i3+i2*i6*i7-i3*i5*i7-i8*i6*i1-i2*i4*i9     !determinant (much faster than calling determ_A)
                                  if ( j /= 1) cycle
                                  Nu=reshape((/i1,i2,i3,i4,i5,i6,i7,i8,i9/),(/3,3/))
                                  Trm=real(Nu)
                                  call Change_Setting_Cell(Cella,Trm,Cellt)
                                  if (Sum(abs(Cellt%cell(:)-Cellb%cell(:)))+Sum(abs(Cellt%ang(:)-Cellb%ang(:))) < tolt  ) then
                                     ok=.true.
                                     exit dox
                                  end if
                               end do    !i9
                            end do     !i8
                         end do      !i7
                      end do       !i6
                   end do        !i5
                end do         !i4
             end do          !i3
          end do           !i2
       end do  dox       !i1
       
       if (.not. ok) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          err_CFML%Msg="Get the transformation Matrix between Cells failed!"
       end if 

       return
    End Subroutine Get_Transfm_Matrix
    
    !!----
    !!---- SUBROUTINE GET_BASIS_FROM_UVW
    !!----
    !!----  Subroutine to construct ZA of type Zone_Axis. This subroutine picks up two reciprocal
    !!----  lattice vectors satisfying the equation
    !!----                            hu+kv+lw=0
    !!----  The two reciprocal lattice vectors have no coprime factors and
    !!----  constitute the basis of a reciprocal lattice plane. They are
    !!----  obtained as the shortest two reciprocal lattice vectors satisfying
    !!----  the above equation. If mode is provided and mode="R", we interpret
    !!----  that the input zone axis is a reciprocal lattice vector and what we
    !!----  obtain is the basis of a direct plane in terms of lattice vectors.
    !!----  If mode="R", dmin corresponds n(uvw)max
    !!----  This subroutine has been imported from resvis_proc.f90.
    !!----
    !!
    Module Subroutine Get_basis_from_uvw(dmin,u,cell,ZoneB,mode)
       !--- Arguments ---!
       real(kind=cp),              intent(in) :: dmin      ! Minimum d-spacing (smax=1/dmin)
       integer, dimension(3),      intent(in) :: u         ! Zone axis indices
       class(CrysCell_Type),       intent(in) :: cell      ! Cell object
       type (Zone_Axis_Type),      intent(out):: ZoneB     ! !Object containing u and basis vector in the plane
       character(len=*), optional, intent(in) :: mode

       !--- Local Variables ---!
       integer                :: n,ik,il,um,iv,i1,i2,i,coun01,coun02,coun1,coun2
       integer                :: kmin,kmax,lmin,lmax
       integer,dimension(3)   :: au,h,mu
       real, dimension(2)     :: rm
       real, dimension(3,3)   :: mat
       integer,dimension(3,2) :: bas
       real                   :: rv,s2max
       logical                :: ok

       !> Init
       ZoneB%nlayer=0
       ZoneB%uvw=u
       
       select type(cell)
          class is (CrysCell_M_Type)
             if (present(mode)) then 
                Mat=cell%GD
             else
                Mat=cell%GR
             end if
                   
          class default
             err_CFML%state=.true.
             err_CFML%Flag=2
             err_CFML%Msg=" The input cell type is incompatible for this procedure Get_basis_from_uvw" 
             return 
       end select
       
       ok=.false.

       au=abs(u)
       um=3*maxval(au)
       i=maxloc(au,dim=1)
       iv=u(i)
       mu=u
       if (iv < 0) then
          mu=-u
          iv=-iv
       end if

       Select Case (i)
          Case(1)
             i1=2; i2=3
             
          Case(2)
             i1=1; i2=3
             
          Case(3)
             i1=1; i2=2
       End Select

       rm(1)=100000.0; rm(2)=rm(1)
       bas(:,1) = (/ 71,121, 113/)
       bas(:,2) = (/117, 91,-111/)

       if (present(mode)) then    ! Direct
          s2max=dmin*dmin   !here dmin is really n_max
          kmax=nint(dmin/Cell%cell(i1)+1.0)
          lmax=nint(dmin/Cell%cell(i2)+1.0)
          kmax=min(um,kmax)
          lmax=min(um,lmax)
          !mat=cell%GD
          
       else
          s2max=1.0/(dmin*dmin)
          kmax=nint(Cell%cell(i1)/dmin+1.0)
          lmax=nint(Cell%cell(i2)/dmin+1.0)
          kmax=min(um,kmax)
          lmax=min(um,lmax)
          !mat=cell%gr
       end if

       kmin=-kmax; lmin=-lmax
       coun1=0; coun2=0
       do ik=kmax,kmin,-1
          do il=lmax,lmin,-1
             if (ik == 0 .and. il == 0) cycle
             n=-ik*mu(i1)-il*mu(i2)
             if (mod(n,iv) == 0) then               !n is multiple of iv
                h(i)= n/iv ; h(i1)=ik ; h(i2) = il  !h is solution of hu+kv+lw=0
                rv=dot_product(real(h),matmul(mat,real(h)))
                if (rv > s2max  .or. rv < 1.0e-20) cycle
                if (rv < rm(1)) then
                   if (.not. co_linear(bas(:,1),h,3) ) then
                      bas(:,2)=bas(:,1)
                      rm(2) = rm(1)
                      if (coun1 >=1) coun2=coun2+1
                   end if
                   bas(:,1)=h
                   rm(1) = rv
                   coun1=coun1+1
                   
                else if (rv < rm(2) .and. .not. co_linear(bas(:,1),h,3) ) then
                   bas(:,2)=h
                   rm(2) = rv
                   coun2=coun2+1
                end if
             end if
          end do
       end do
       
       ZoneB%rx=bas(:,1)
       ZoneB%ry=bas(:,2)
       if (coun1 >= 1 .and. coun2 >=1) ok=.true.
       
       coun01=0; coun02=0; coun1=0; coun2=0
       do i=1,3
          if (ZoneB%rx(i) < 0) coun1=coun1+1
          if (ZoneB%ry(i) < 0) coun2=coun2+1
          if (ZoneB%rx(i) == 0) coun01=coun01+1
          if (ZoneB%ry(i) == 0) coun02=coun02+1
       end do
       if (coun1 >= 2 .or. (coun1 == 1 .and. coun01 == 2)) ZoneB%rx=-ZoneB%rx
       if (coun2 >= 2 .or. (coun2 == 1 .and. coun02 == 2)) ZoneB%ry=-ZoneB%ry
       
       if (.not. ok) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          Err_CFML%Msg=" Problems in the Procedure of Get_Basis_From_Uvw"
       end if 

       return
    End Subroutine Get_Basis_From_Uvw
    
    !!----
    !!---- SUBROUTINE GET_TWOFOLD_AXES
    !!----
    !!----    Subroutine for getting the possible two-fold axes (within an
    !!----    angular tolerance tol) existing in the lattice generated by the
    !!----    unit cell "Celln". Strictly independent two-fold axes are stored
    !!----    in the variable "twofold" that is of type Twofold_Axes_Type
    !!----    The output order of the two-fold axes is ascending in their
    !!----    modulus. Shorter vectors appears before longer ones.
    !!----    The conditions for a reciprocal or direct row to be a two-fold
    !!----    axis are discussed by Y. Le Page in J.Appl.Cryst. 15, 255 (1982).
    !!----
    !!----
    !!---- Update: November - 2008
    !!
    Module Subroutine Get_TwoFold_Axes(Cell,Tol,Twofold)
       !---- Arguments ----!
       class(CrysCell_M_Type),   intent (in) :: Cell     ! Cell object
       real(kind=cp),           intent (in) :: Tol      ! angular tolerance in degrees
       Type(twofold_axes_type), intent(out) :: Twofold

       !---- Local variables ----!
       integer                        :: i,j,n,m, ih,ik,il,iu,iv,iw,imax,ntwo
       real(kind=cp), dimension(3)    :: dv, rv, a, b, c, as, bs, cs, cross
       real(kind=cp), dimension(  12) :: maxes,crossa
       integer, dimension(  12)       :: dota,ind
       real(kind=cp), dimension(3,12) :: caxes
       integer, dimension(3,12)       :: dtw,rtw
       integer, dimension(3)          :: v,h
       real(kind=cp)                  :: dot,crossm

       maxes=0.0; crossa=0.0; dota=0; caxes=0.0; dtw=0; rtw=0
       a=Cell%Cr_Orth_cel(:,1)
       b=Cell%Cr_Orth_cel(:,2)
       c=Cell%Cr_Orth_cel(:,3)
       twofold%a=a
       twofold%b=b
       twofold%c=c
       as=cross_product(b,c)/Cell%Vol !Reciprocal lattice vectors in
       bs=cross_product(c,a)/Cell%Vol !Cartesian components
       cs=cross_product(a,b)/Cell%Vol
       ntwo=0
       imax=2   !Is inough if the input cell is the Buerger or Niggli cell

       do_iu: do iu=imax, 0,-1
          do iv=imax,-imax,-1
             do iw=imax,-imax,-1
                v=(/iu,iv,iw/)
                if (.not. Co_Prime(v,2)) cycle
                do ih=imax,0,-1
                   do ik=imax,-imax,-1
                      do_il:do il=imax,-imax,-1
                         h=(/ih,ik,il/)
                         if (.not. Co_Prime(h,2)) cycle
                         n=abs(ih*iu+ik*iv+il*iw)
                         if ( n == 2 .or. n == 1) then
                            dv=real(iu)*a+real(iv)*b+real(iw)*c
                            rv=real(ih)*as+real(ik)*bs+real(il)*cs
                            cross=cross_product(dv,rv)
                            dot=sqrt(dot_product(cross,cross))
                            crossm=atand(dot/real(n))
                            if (abs(crossm) <= tol) then
                               do m=1,ntwo
                                  if (determ_Vec((/17,41,71/),v,dtw(:,m) ) == 0) cycle do_il
                               end do
                               ntwo=ntwo+1
                               dtw(:,ntwo)= v
                               dv=v(1)*a+v(2)*b+v(3)*c
                               caxes(:,ntwo)=dv
                               maxes(ntwo)=sqrt(dot_product(dv,dv))
                               rtw(:,ntwo)= h
                               dota(ntwo)=n
                               crossa(ntwo)=crossm
                            end if
                            if (ntwo == 12) exit do_iu
                         end if
                      end do do_il
                   end do
                end do
             end do
          end do
       end do do_iu
       
       call sort(maxes,ntwo,ind)
       do i=1,ntwo
          j=ind(i)
          twofold%dtwofold(:,i)= dtw(:,j)
          twofold%caxes(:,i)= caxes(:,j)
          twofold%maxes(i)= maxes(j)
          twofold%rtwofold(:,i)= rtw(:,j)
          twofold%dot(i)= dota(j)
          twofold%cross(i)= crossa(j)
       End do
       twofold%ntwo=ntwo
       twofold%tol=tol

       return
    End Subroutine Get_TwoFold_Axes
    
    !!----
    !!---- SUBROUTINE GET_CONVENTIONAL_CELL
    !!----
    !!----  This subroutine provides the "conventional" (or quasi! being still tested )
    !!----  from the supplied object "twofold" that has been obtained from a previous
    !!----  call to Get_TwoFold_Axes. The conventional unit cell can be deduced from
    !!----  the distribution of two-fold axes in the lattice. The cell produced in this
    !!----  procedure applies some rules for obtaining the conventional cell, for instance
    !!----  in monoclinic lattices (a single two-fold axis) the two-fold axis is along
    !!----  b and the final cell is right handed with a <= c and beta >= 90. It may be
    !!----  A,C or I centred. The convertion to the C-centred setting in the A and I
    !!----  centring, is not attempted. The angular tolerance for accepting a two-fold
    !!----  axis, or higher order axes, as such has been previously set into twofold%tol
    !!----  component. The output Tr-matrix is the transpose of the IT convention.
    !!----  It corresponds to the transformation between formal column matrices containing
    !!----  the basis vectors.
    !!----  The tolerance for comparing distances in angstroms told is optional.
    !!----- By default the used tolerance is 0.2 angstroms.
    !!----
    !!---- Update: November - 2008
    !!----
    Module Subroutine Get_Conventional_Cell(Twofold,Cell,Tr,Message,told)
       !---- Arguments ----!
       Type(Twofold_Axes_Type), intent(in)  :: Twofold
       class(CrysCell_Type),    intent(out) :: Cell
       integer, dimension(3,3), intent(out) :: tr
       character(len=*),        intent(out) :: message
       real(kind=cp), optional, intent(in)  :: told

       !---- Local variables ----!
       integer, dimension(1)          :: ix
       integer, dimension(2)          :: ab
       integer, dimension(3)          :: rw,h1,h2
       integer, dimension(66)         :: inp
       integer, dimension(3,48)       :: row
       real(kind=cp), dimension(3)    :: u,v1,v2,v3,a,b,c,vec,vi,vj,vk
       real(kind=cp), dimension(48)   :: mv
       real(kind=cp), dimension(66)   :: ang
       integer                        :: iu,iv,iw,nax,i,j,k,m,namina,naminb,naminc,ia
       real(kind=cp)                  :: dot,ep,domina,dominb,dominc,aij,aik,ajk
       real(kind=cp)                  :: delt,tola
       logical                        :: ok,hexap, hexac

       !> Init
       a=twofold%a; b=twofold%b; c=twofold%c
       delt=twofold%tol
       ep=cosd(90.0-delt)
       domina=9.0e+30; dominc=domina
       tr=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
       ab=0; mv=0.0; ang=0.0; row=0; inp=0
       Call Set_Crystal_Cell([1.0,1.0,1.0],[90.0,90.0,90.0],Cell)
       message=" "
       call clear_error()
       
       ok=.true.
       tola=0.2
       if(present(told)) tola=told

       Select Case(twofold%ntwo)
          Case (1)    !Monoclinic n-2foldaxes=1
             v2=twofold%caxes(:,1)
             u = v2/twofold%maxes(1)
             tr(2,:)=twofold%dtwofold(:,1)
             nax=0
             do iu=-3,3
                do iv=-3,3
                   do_iw: do iw=0,3
                      rw=(/iu,iv,iw/)
                      !> if(iu == 0 .and. iv == 0 .and. iw == 0) cycle
                      if (.not. Co_prime(rw,3)) cycle
                      vec=real(iu)*a+real(iv)*b+real(iw)*c
                      dot=sqrt(dot_product(vec,vec))
                      vec=vec/dot
                      if (abs(dot_product(u,vec)) < ep) then
                         do m=1,nax
                            if(co_linear(rw,row(:,m),3)) cycle do_iw
                         end do
                         nax=nax+1
                         row(:,nax) = rw
                         mv(nax) = dot
                         if (dot < domina) then
                            domina=dot
                            namina=nax
                            tr(1,:)=rw
                            v1=real(iu)*a+real(iv)*b+real(iw)*c
                         end if
                      end if
                   end do do_iw
                end do
             end do

             do i=1,nax
                if (i == namina) cycle
                if (mv(i) < dominc) then
                   dominc=mv(i)
                   naminc=i
                end if
             end do
             tr(3,:)=row(:,naminc)
             v3=row(1,naminc)*a+row(2,naminc)*b+row(3,naminc)*c

             !> Length of the three basis vectors should be stored in mv(1),mv(2),mv(3)
             mv(1)=sqrt(dot_product(v1,v1))
             mv(2)=sqrt(dot_product(v2,v2))
             mv(3)=sqrt(dot_product(v3,v3))

             !>The two shortest vectors perpendicular to the primary twofold axis have been found
             !> and the transformation matrix has been constructed
             namina=determ_3x3(tr)
             if (namina < 0) then   !right handed system
                tr(2,:)=-tr(2,:)
                v2=-v2
                namina=-namina
             end if

             !> Test if beta is lower than 90 in such a case invert c and b
             dominb=dot_product(v1/mv(1),v3/mv(3))
             if (dominb > 0.0) then  !angle beta < 90
                tr(2,:)=-tr(2,:)
                v2=-v2
                tr(3,:)=-tr(3,:)
                v3=-v3
             end if

             Select Case (namina)
                Case(1)
                   message="Monoclinic, primitive cell"
                   
                Case(2)
                   rw=matmul((/0,1,1/),tr)
                   if (.not. co_prime(rw,3)) then
                      message="Monoclinic, A-centred cell"
                      
                   else
                      rw=matmul((/1,1,1/),tr)
                      if (.not. co_prime(rw,3)) then
                         message="Monoclinic, I-centred cell"
                      else
                         rw=matmul((/1,1,0/),tr)
                         if(.not. co_prime(rw,3)) message="Monoclinic, C-centred cell"
                      end if
                   end if

                Case(3:)
                   message="Error in monoclinic cell"
                   ok=.false.
                   
                   err_CFML%state=.true.
                   err_CFML%Flag=2
                   err_CFML%Msg="Error for monoclinic cell in Get_Conventional_Cell"
                   return
             End Select

          Case (3)    !Orthorhombic/Trigonal n-2foldaxes=3
             u(1:3)=twofold%maxes(1:3)
             ix=minloc(u)
             namina=ix(1)
             ix=maxloc(u)
             naminc=ix(1)
             if (naminc == namina) then
                namina=1; naminb=2; naminc=3
             else
                do i=1,3
                   if(i == namina) cycle
                   if(i == naminc) cycle
                   naminb=i
                   exit
                end do
             end if
             tr(1,:) = twofold%dtwofold(:,namina)
             tr(2,:) = twofold%dtwofold(:,naminb)
             tr(3,:) = twofold%dtwofold(:,naminc)
             v1 = twofold%caxes(:,namina)
             v2 = twofold%caxes(:,naminb)
             v3 = twofold%caxes(:,naminc)
             mv(1)=twofold%maxes(namina)
             mv(2)=twofold%maxes(naminb)
             mv(3)=twofold%maxes(naminc)

             !> Check the system by verifying that the two-fold axes form 90 (orthorhombic)
             !> or 120 degrees (Trigonal)
             domina=dot_product(v2/mv(2),v3/mv(3))
             dominb=dot_product(v1/mv(1),v3/mv(3))
             dominc=dot_product(v1/mv(1),v2/mv(2))

             if (abs(domina) < ep .and. abs(dominb) < ep .and. abs(dominc) < ep) then !orthorhombic
                namina=determ_3x3(tr)
                if (namina < 0) then
                   tr(3,:)=-tr(3,:)
                   v3=-v3
                   namina=-namina
                end if
                Select Case (namina)
                   Case(1)
                      message="Orthorhombic, Primitive cell"

                   Case(2)
                      rw=matmul((/0,1,1/),tr)
                      if (.not. co_prime(rw,3)) then
                         message="Orthorhombic, A-centred cell"
                      else
                         rw=matmul((/1,1,1/),tr)
                         if (.not. co_prime(rw,3)) then
                            message="Orthorhombic, I-centred cell"
                         else
                            rw=matmul((/1,1,0/),tr)
                            if (.not. co_prime(rw,3)) then
                               message="Orthorhombic, C-centred cell"
                            else
                               rw=matmul((/1,0,1/),tr)
                               if (.not. co_prime(rw,3)) message="Orthorhombic, B-centred cell"
                            end if
                         end if
                      end if

                   Case(3:)
                      message="Orthorhombic, F-centred cell"
                End Select

             else !Rhombohedral/Trigonal

                !> In the Trigonal system the two-fold axes are in the plane perpendicular to
                !> the three-fold axis, and valid a,b, vectors can be chosen among any two two-fold
                !> axes forming an angle of 120 degrees
                !> verify that 1 and 2 form 120
                ang(1)=acosd(domina)    !2-3
                ang(2)=acosd(dominb)    !1-3
                ang(3)=acosd(dominc)    !1-2
                dot=1.0
                iu=1
                j=0
                do i=1,3
                   if (abs(ang(i)-120.0) < delt) then
                      j=i
                      exit
                   end if
                end do

                if ( j == 0) then
                   do i=1,3
                      if (abs(ang(i)-60.0) < delt) then
                         j=i
                         dot=-1.0
                         iu=-1
                         exit
                      end if
                   end do
                End if

                if ( j == 0) then
                   message="Trigonal/Rhombohedral test failed! Supply only one two-fold axis"
                   ok=.false.
                   
                   err_CFML%state=.true.
                   err_CFML%Flag=2
                   err_CFML%Msg=trim(message)
                   return    
                else
                   Select Case (j)
                      case(1)
                         vi=v2
                         vj=dot*v3
                         h1=tr(2,:); h2=iu*tr(3,:)
                         tr(3,:)=tr(1,:)
                         tr(1,:)=h1
                         tr(2,:)=h2

                      case(2)
                         vi=v1
                         vj=dot*v3
                         h2=iu*tr(3,:)
                         tr(3,:)=tr(2,:)
                         tr(2,:)=h2

                      case(3)
                         vi=v1
                         vj=dot*v2
                         tr(2,:)=iu*tr(2,:)

                   End Select

                   v1 = vi
                   v2 = vj
                   mv(1)=sqrt(dot_product(v1,v1))
                   mv(2)=sqrt(dot_product(v2,v2))
                   vi=v1/mv(1)
                   vj=v2/mv(2)
                   
                   ok=.false.

                   do_iu: do iu=-3,3
                      do iv=-3,3
                         do iw=0,3
                            rw=(/iu,iv,iw/)
                            if (.not. Co_prime(rw,3)) cycle
                            vec=real(iu)*a+real(iv)*b+real(iw)*c
                            dot=sqrt(dot_product(vec,vec))
                            vec=vec/dot
                            if (abs(dot_product(vi,vec)) < ep  .and. abs(dot_product(vj,vec)) < ep) then
                               tr(3,:)=rw
                               ok=.true.
                               exit do_iu
                            end if
                         end do
                      end do
                   end do do_iu

                   If (ok) then
                      namina=determ_3x3(tr)
                      if (namina < 0) then
                         tr(3,:)=-tr(3,:)
                         namina=-namina
                      end if
                      v3 = tr(3,1)*a+tr(3,2)*b+tr(3,3)*c
                      mv(3)=sqrt(dot_product(v3,v3))
                      Select Case (namina)
                         case(1)
                            message="Primitive hexagonal cell"
                         case(3)
                            rw=matmul((/2,1,1/),tr)
                            if (.not. co_prime(rw,3)) then
                               message="Rhombohedral, obverse setting cell"
                            else
                               message="Rhombohedral, reverse setting cell"
                            end if
                      End Select

                   Else
                      message="Trigonal/Rhombohedral test failed! Supply only one two-fold axis"
                      ok=.false.
                      
                      err_CFML%state=.true.
                      err_CFML%Flag=2
                      err_CFML%Msg=trim(message)
                      return  
                   End if
                End if !j==0
             End if  !orthorhombic test

          Case (5)    !Tetragonal n-2foldaxes=5
             m=0
             inp=0
             mv(1:5)=twofold%maxes(1:5)
             do i=1,twofold%ntwo-1
                vi=twofold%caxes(:,i)/twofold%maxes(i)
                do j=i+1,twofold%ntwo
                   vj=twofold%caxes(:,j)/twofold%maxes(j)
                   m=m+1
                   ang(m)=acosd(dot_product(vi,vj))
                   if (abs(ang(m)-45.0) < delt .or. abs(ang(m)-135.0) < delt) then
                      inp(i)=1
                      inp(j)=1
                      if (mv(i) > mv(j)) then
                         ia=j
                      else
                         ia=i
                      end if
                      if (ab(1) == 0) then
                         ab(1) = ia
                      else
                         ab(2) = ia
                      end if
                   end if
                end do
             end do

             !Determination of the c-axis (that making 90 degree with all the others)
             ix=minloc(inp)
             naminc=ix(1)

             !The two axes forming a,b are those of indices ab(1) and ab(2)
             namina=ab(1)
             naminb=ab(2)
             if (namina == 0 .or. naminb == 0) then
                ok=.false.
                message="Basis vectors a-b not found!"
                
                err_CFML%state=.true.
                err_CFML%Flag=2
                err_CFML%Msg=trim(message)
                return  
             else   
                tr(1,:) = twofold%dtwofold(:,namina)
                tr(2,:) = twofold%dtwofold(:,naminb)
                tr(3,:) = twofold%dtwofold(:,naminc)
                v1 = twofold%caxes(:,namina)
                v2 = twofold%caxes(:,naminb)
                v3 = twofold%caxes(:,naminc)
                mv(1)=twofold%maxes(namina)
                mv(2)=twofold%maxes(naminb)
                mv(3)=twofold%maxes(naminc)
                namina=determ_3x3(tr)
                if (namina < 0) then
                   tr(3,:)=-tr(3,:)
                   v3=-v3
                   namina=-namina
                end if
                
                Select Case (namina)
                   Case(1)
                      message="Tetragonal, Primitive cell"
                   Case(2)
                      message="Tetragonal, I-centred cell"
                   Case(3:)
                      message="Error in tetragonal cell"
                      ok=.false.
                      
                      err_CFML%state=.true.
                      err_CFML%Flag=2
                      err_CFML%Msg=trim(message)
                      return  
                End Select
             end if

          Case (7)    !Hexagonal n-2foldaxes=7

             m=0
             inp=0
             mv(1:7)=twofold%maxes(1:7)
             hexap=.false.;  hexac=.false.

             !Search tha a-b plane
             do_ii:do i=1,twofold%ntwo-1
                vi=twofold%caxes(:,i)/twofold%maxes(i)
                do j=i+1,twofold%ntwo
                   vj=twofold%caxes(:,j)/twofold%maxes(j)
                   aij=acosd(dot_product(vi,vj))
                   if (abs(aij-120.0) < delt) then
                      if (abs(mv(i)-mv(j)) < tola .and. .not. hexap ) then
                         rw(1)=i; rw(2)=j
                         u(1)=mv(i); u(2)=mv(j)
                         hexap=.true.
                         exit do_ii
                      end if
                   end if
                end do
             end do do_ii

             if (hexap) then ! Search the c-axis, it should be also a two-fold axis!
                             ! because Op(6).Op(6).Op(6)=Op(2)
                v1 = twofold%caxes(:,rw(1))
                v2 = twofold%caxes(:,rw(2))
                vj=v1/u(1)
                vk=v2/u(2)
                do i=1,twofold%ntwo
                   vi=twofold%caxes(:,i)/twofold%maxes(i)
                   aij=acosd(dot_product(vi,vj))
                   aik=acosd(dot_product(vi,vk))
                   if (abs(aij-90.0) < delt .and. abs(aik-90.0) < delt ) then
                      rw(3)=i
                      u(3)= mv(i)
                      hexac=.true.
                      exit
                   end if
                end do
             else
                ok=.false.
                
                err_CFML%state=.true.
                err_CFML%Flag=2
                err_CFML%Msg="Error in Hexagonal n-2fold axes=7"
                return
             end if

             if (hexac) then
                do i=1,3
                   tr(i,:) = twofold%dtwofold(:,rw(i))
                   mv(i)=u(i)
                end do
                v3 = twofold%caxes(:,rw(3))
                namina=determ_3x3(tr)
                if (namina < 0) then
                   tr(3,:)=-tr(3,:)
                   v3=-v3
                   namina=-namina
                end if

                Select Case (namina)
                   Case(1)
                      message="Hexagonal, Primitive cell"
                   Case(2:)
                      message="Hexagonal, centred cell? possible mistake"
                End Select

             else
                ok=.false.
                message="The c-axis of a hexagonal cell was not found!"
                
                err_CFML%state=.true.
                err_CFML%Flag=2
                err_CFML%Msg=trim(message)
                return
             end if

          Case (9)   !Cubic n-2foldaxes=9
             m=0
             inp=0
             mv(1:9)=twofold%maxes(1:9)
             do_i:do i=1,twofold%ntwo-2
                vi=twofold%caxes(:,i)/twofold%maxes(i)
                do j=i+1,twofold%ntwo-1
                   vj=twofold%caxes(:,j)/twofold%maxes(j)
                   do k=j+1,twofold%ntwo
                      vk=twofold%caxes(:,k)/twofold%maxes(k)
                      aij=acosd(dot_product(vi,vj))
                      aik=acosd(dot_product(vi,vk))
                      ajk=acosd(dot_product(vj,vk))
                      if (abs(aij-90.0) < delt .and. abs(aik-90.0) < delt .and. abs(ajk-90.0) < delt ) then
                         if (abs(mv(i)-mv(j)) < tola .and. abs(mv(i)-mv(k)) < tola .and. abs(mv(j)-mv(k)) < tola ) then
                            rw(1)=i; rw(2)=j; rw(3)=k
                            u(1)=mv(i); u(2)=mv(j); u(3)=mv(k)
                            exit do_i
                         end if
                      end if
                   end do
                end do
             end do do_i

             do i=1,3
                tr(i,:) = twofold%dtwofold(:,rw(i))
                mv(i)=u(i)
             end do
             v1 = twofold%caxes(:,rw(1))
             v2 = twofold%caxes(:,rw(2))
             v3 = twofold%caxes(:,rw(3))
             namina=determ_3x3(tr)
             if (namina < 0) then
                tr(3,:)=-tr(3,:)
                v3=-v3
                namina=-namina
             end if

             Select Case (namina)
                Case(0)
                  message="Pseudo-cubic but tolerance too small ... "
                  ok=.false.
                  
                  err_CFML%state=.true.
                  err_CFML%Flag=2
                  err_CFML%Msg=trim(message)
                  return
                Case(1)
                   message="Cubic, Primitive cell"
                Case(2)
                   rw=matmul((/0,1,1/),tr)
                   if (.not. co_prime(rw,3)) then
                      message="Cubic, A-centred cell"
                   else
                      rw=matmul((/1,1,1/),tr)
                      if (.not. co_prime(rw,3)) then
                         message="Cubic, I-centred cell"
                      else
                         rw=matmul((/1,1,0/),tr)
                         if (.not. co_prime(rw,3)) then
                            message="Cubic, C-centred cell"
                         else
                            rw=matmul((/1,0,1/),tr)
                            if (.not. co_prime(rw,3)) message="Cubic, B-centred cell"
                         end if
                      end if
                   end if

                Case(3:)
                  message="Cubic, F-centred cell"
             End Select

          case default
             write(unit=message,fmt="(a,i3)") "Wrong number of two-fold axes! ",twofold%ntwo
             ok=.false.
             
             err_CFML%state=.true.
             err_CFML%Flag=2
             err_CFML%Msg=trim(message)
             return

       End Select

      !> Calculation of the new cell
      ang(1)=acosd(dot_product(v2/mv(2),v3/mv(3)))
      ang(2)=acosd(dot_product(v1/mv(1),v3/mv(3)))
      ang(3)=acosd(dot_product(v1/mv(1),v2/mv(2)))
      Call Set_Crystal_Cell(mv(1:3),ang(1:3),Cell)

      return
    End Subroutine Get_Conventional_Cell
        
 
End Submodule GenMetrics 