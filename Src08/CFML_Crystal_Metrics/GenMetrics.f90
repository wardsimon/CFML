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
        
 
End Submodule GenMetrics 