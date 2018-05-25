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
        
 
End Submodule GenMetrics 