 SubModule (CFML_Geometry_SXTAL) SXTAL_UB

  Contains

    !!----
    !!---- Module Subroutine cell_fr_UB(ub,ipr,dcel,rcel)
    !!---- real(kind=cp),Dimension(3,3),         Intent(In)  :: ub
    !!---- Integer, optional,                    Intent(In)  :: ipr
    !!---- real(kind=cp),Dimension(6), optional, Intent(out) :: dcel,rcel
    !!----
    !!----    Calculate and print cell parameters from UB-matrix
    !!----
    !!---- Update: May 2011, June 2020
    !!
    Module Subroutine cell_fr_UB(ub,ipr,dcel,rcel)
       !---- Arguments ----!
       real(kind=cp),Dimension(3,3),         Intent(In)  :: ub
       Integer, optional,                    Intent(In)  :: ipr
       real(kind=cp),Dimension(6), optional, Intent(out) :: dcel,rcel

       !--- Local Variables ---!
       real(kind=dp), Dimension(3,3) :: g,ubinv,angle
       real(kind=cp), Dimension(3,3) :: sg
       real(kind=dp), Dimension(3)   :: acal,angcal,cala,calang
       integer                       :: i,j,k,jn,kn
       real(kind=dp)                 :: x

       g=Matmul(Transpose(ub),ub)
       sg=g
       !..inverse matrix g=b'b  to get direct cell parameters
       !..On first run through loop calculation of reciprocal cell, second do real cell
       Do k=1,2
           If(k==2) sg=invert(sg)
           g=sg
           Do  i=1,3
             acal(i)=Sqrt(g(i,i))
           End Do
           Do  i=1,3
             j=i
             jn=Mod(j,3)+1
             kn=Mod(jn,3)+1
             angcal(i)=acosd(g(jn,kn)/(acal(jn)*acal(kn)))
           End Do
           If(k==2) Exit
           cala   = acal   !store reciprocal latice in the first pass
           calang = angcal
       End Do
       if(present(dcel)) then
         dcel(1:3)=acal
         dcel(4:6)=angcal
       end if
       if(present(rcel)) then
         rcel(1:3)=cala
         rcel(4:6)=calang
       end if
       if(present(ipr)) then
         !.....Now invert UB to obtain the hkl's along the orthogonal diffractometer axes
         ubinv=invert(ub)
         Write(Unit=ipr,Fmt="(/,a)")                " => Parameters deduced from the UB matrix "
         Write(Unit=ipr,Fmt="(a,3f12.5,tr4,3f9.4)") " => Direct     cell dimensions: ",acal,angcal
         Write(Unit=ipr,Fmt="(a,3f12.8,tr4,3f9.4)") " => Reciprocal cell dimensions: ",cala,calang
         Write(Unit=ipr,Fmt="(/,a)")" =>               UB-Matrix                                          Inverse of UB-Matrix "
         Write(Unit=ipr,Fmt="(a)")    "           A*            B*            C*                 X(PH=0,CH=0)  Y(PH=90,CH=0)    Z(CHI=90)"
         Write(Unit=ipr,Fmt="(a,3f14.8,a,3f14.8)") "  X",ub(1,:),"          H",ubinv(1,:)
         Write(Unit=ipr,Fmt="(a,3f14.8,a,3f14.8)") "  Y",ub(2,:),"          K",ubinv(2,:)
         Write(Unit=ipr,Fmt="(a,3f14.8,a,3f14.8)") "  Z",ub(3,:),"          L",ubinv(3,:)
         !.....Now calculate angles between recip axes and orthogonal diffract. axes
         Do  i=1,3
           Do  j=1,3
             x = ub(i,j)/cala(j)
             If (x > 1.0) Then
                Write (ipr,"(a,3e12.4)") " Error x >1.0! Values of x,ub,acal: ",x,ub(i,j),acal(i)
             Else
               angle(i,j) = acosd(x)
             End If
           End Do
         End Do
         Write(Unit=ipr,Fmt="(/,a)")"    With all diffractometer angles set to 0, the angles between the "
         Write(Unit=ipr,Fmt="(a)")  "    reciprocal (A*,B*,C*) and diffractometer (X,Y,Z) axes are ..."
         Write(Unit=ipr,Fmt="(a,3f12.4,a)") "    (A*-X  B*-X  C*-X)  (",angle(1,:),")"
         Write(Unit=ipr,Fmt="(a,3f12.4,a)") "    (A*-Y  B*-Y  C*-Y)  (",angle(2,:),")"
         Write(Unit=ipr,Fmt="(a,3f12.4,a)") "    (A*-Z  B*-Z  C*-Z)  (",angle(3,:),")"
       end if
    End Subroutine cell_fr_UB

    !!----
    !!---- Module Function genb(c) Result (b)
    !!----    Type(Cell_G_Type), Intent(In)   :: c
    !!----    real(kind=cp), Dimension(3,3)   :: b
    !!----
    !!--<<   Calculation of [B] matrix
    !!----   Busing&Levy Acta Cryst.(1967)22,457-464  Equation 3
    !!----   Wooster R. Sci. Instrum. (1962)39,103
    !!----      C%cell  : Direct cell parameters
    !!----      C%rcell : Reciprocal cell parameters
    !!----      C%ang   : Direct cell angles
    !!-->>      C%rcell : Reciprocal cell angles
    !!----
    !!---- Update: April 2008, June 2020
    !!
    Module Function genb(c) Result (b)
       !---- Arguments ----!
       Type(Cell_G_Type), Intent(In)  :: c
       real(kind=cp), Dimension(3,3)  :: b

       b(1,1)=c%rcell(1)
       b(1,2)=c%rcell(2)*cosd(c%rang(3))
       b(1,3)=c%rcell(3)*cosd(c%rang(2))
       b(2,2)=c%rcell(2)*sind(c%rang(3))
       b(2,3)=-(c%rcell(3)*sind(c%rang(2))*cosd(c%ang(1)))
       b(3,3)=1.0_cp/c%cell(3)
       b(2,1)=0.0_cp
       b(3,1)=0.0_cp
       b(3,2)=0.0_cp
    End Function genb

    !!----
    !!---- Module Function GenUB(b,h1,h2,h1o,h2o) Result (ub)
    !!----    real(kind=cp), Dimension(3,3), Intent(In)  :: B        !Busing-Levy B-matrix
    !!----    real(kind=cp), Dimension(3),   Intent(In)  :: h1,h2    !Miller indices
    !!----    real(kind=cp), Dimension(3),   Intent(In)  :: h1o,h2o  !Components in Lab system
    !!----
    !!----    Original from   A.Filhol  25/05/84
    !!----    Given the [B] matrix, the Miller indices of two reflections, h1 & h2,
    !!----    and the components of these two reflections, h1o & h2o, in the laboratory
    !!----    system, this subroutine provides the matrix UB. Only the direction in the
    !!----    laboratory system of reflections are needed, e.g. h1o and h2o may be unitary
    !!----    vectors or whatever other vector along these directions.
    !!----
    !!----    [hc] : Reflection H,K,L in the reciprocal lattice orthonormal system
    !!----    [hc] = [B] [h]
    !!----    [ho] : Reflection H,K,L in the Cartesian laboratory system
    !!----    [ho]=[UB]*[hc]
    !!----
    !!---- Update: April 2008, June 2020
    !!
    Module Function GenUB(b,h1,h2,h1o,h2o) Result (UB)
       !---- Arguments ----!
       real(kind=cp), Dimension(3,3), Intent(In)  :: B        !Busing-Levy B-matrix
       real(kind=cp), Dimension(3),   Intent(In)  :: h1,h2    !Miller indices
       real(kind=cp), Dimension(3),   Intent(In)  :: h1o,h2o  !Components in Lab system
       real(kind=cp), Dimension(3,3)              :: UB       !Busing-Levy UB-matrix

       !--- Local Variables ---!
       real(kind=cp), Dimension(3)  :: h1c,h2c,v1,v2
       real(kind=cp), Dimension(3,3):: trpc,trpo,u

       call clear_error()
       !Calculation of the reciprocal components in the cartesian reciprocal axes
       h1c=Matmul(B,h1)
       h2c=Matmul(B,h2)
       v1=h1o
       v2=h2o
       Call triple(h1c,h2c,trpc) !Orthonormal frame attached to h1c,h2c
       If (err_CFML%Ierr /= 0) Then
          ub=0.0_cp
          Return
       End If
       Call triple(v1,v2,trpo) !Orthonormal frame attached to h1o,h2o
       If (err_CFML%Ierr /= 0) Then
          ub=0.0_cp
          Return
       End If
       !..... MATRIX [U]  *** Equation 27 (B-L)
       U=Matmul(trpo,Transpose(trpc))

       !..... MATRIX [UB]
       UB=Matmul(U,B)
    End Function GenUB

    !!---- Module Function Get_UB_from_hkl_hkl_omega(wave,Cell,h1,h2,omega) Result(UB)
    !!----   real(kind=cp),                 intent(in)    :: wave  !Vavelength
    !!----   type (Cell_G_Type),            intent(in)    :: Cell  !Unit cell object
    !!----   real(kind=cp), dimension(3),   intent(in)    :: h1    !Indices of the known first reflection in plane
    !!----   real(kind=cp), dimension(3),   intent(in)    :: h2    !Indices of the known second reflection in plane
    !!----   real(kind=cp),                 intent(in)    :: omega !Value of the omega motor for the second reflection (vertical spindle)
    !!----   real(kind=cp), dimension(3,3)                :: UB    !Generated Busing-Levy UB-matrix
    !!----
    !!----   This subroutine generates a UB matrix when two reflections in the horizontal plane
    !!----   are known (indices hkl and h'h'l') and the second reflection has been measured and
    !!----   its omega angle is known.
    !!----
    !!----   Updated: June-2012 (JRC), June 2020
    !!----
    Module Function Get_UB_from_hkl_hkl_omega(wave,Cell,h1,h2,omega) Result(UB)
      real(kind=cp),                 intent(in)    :: wave
      type (Cell_G_Type),            intent(in)    :: Cell
      real(kind=cp), dimension(3),   intent(in)    :: h1
      real(kind=cp), dimension(3),   intent(in)    :: h2
      real(kind=cp),                 intent(in)    :: omega
      real(kind=cp), dimension(3,3)                :: UB
      ! Local variables
      integer                       :: ierr
      real(kind=cp)                 :: theta1,theta2,alpha,del_omega,d1s,d2s
      real(kind=cp), dimension(3)   :: ho1,ho2,s1,s2
      real(kind=cp), dimension(3,3) :: Rot

      call clear_error()
      if(Co_linear(h1,h2,3)) then
        Err_CFML%Ierr=1
        Err_CFML%Msg="The two provided reflections are co-linear, no UB-matrix can be calculated"
        return
      end if
      !
      !Determination of the Bragg angles of two reflections in the scattering plane
      !
      d1s=sqrt(dot_product(h1,matmul(Cell%GR,h1))) !d1*
      d2s=sqrt(dot_product(h2,matmul(Cell%GR,h2))) !d2*
      theta1=asind(0.5*wave*d1s)
      theta2=asind(0.5*wave*d2s)
      alpha=acosd(dot_product(h1,matmul(Cell%GR,h2))/d1s/d2s) !Angle between the two reciprocal vectors
      del_omega=theta1-theta2+alpha  !Variation in omega to put the first reflection in diffraction position
      !
      !Calculation of the Cartesian components of the two reflections in the scattering plane
      !
      s2=d2s*(/cosd(Theta2),-sind(Theta2),0.0_cp/)   !z4   diffraction position
      Rot= Phi_Mat(omega)
      ho2=matmul(transpose(rot),s2)                  !z1   zero motor angles
      s1=d1s*(/cosd(Theta1),-sind(Theta1),0.0_cp/)   !z4   diffraction position
      Rot= Phi_mat(omega+del_omega)
      ho1=matmul(transpose(rot),s1)                  !z1   zero motor angles
      !
      ! Generate UB-matrix
      !
      UB= GenUB(Cell%BL_M,h1,h2,ho1,ho2)
      if(Err_CFML%Ierr /= 0) then
        Err_CFML%Msg = "Error in the calculation of UB-matrix hkl_hkl_omega: "//trim(Err_CFML%Msg)
      end if

    End Function Get_UB_from_hkl_hkl_omega

    !!---- Module Function Get_UB_from_uvw_hkl_omega(wave,Cell,Zone_Axis,h1,omega) Result(UB)
    !!----   real(kind=cp),                 intent(in)    :: wave      !Vavelength
    !!----   type (Cell_G_Type),            intent(in)    :: Cell      !Unit cell object
    !!----   Type (Zone_Axis_type),         intent(in out):: Zone_Axis !Zone axis (See CFML_Crystal_Metrics module)
    !!----   real(kind=cp), dimension(3),   intent(in)    :: h1        !Indices of the known reflection in plane
    !!----   real(kind=cp),                 intent(in)    :: omega     !Value of the omega motor (vertical spindle)
    !!----   real(kind=cp), dimension(3,3)                :: UB        !Generated Busing-Levy UB-matrix
    !!----
    !!----   This subroutine generates a UB matrix when the vertical zone axis of the crystal
    !!----   is known and a reflection in the horizonal plane has been measured with known
    !!----   indices an value of the omega angle.
    !!----
    !!----   Updated: June-2012 (JRC), June 2020
    !!----
    Module Function Get_UB_from_uvw_hkl_omega(wave,Cell,Zone_Axis,h1,omega) Result(UB)
      real(kind=cp),                 intent(in)    :: wave
      type (Cell_G_Type),            intent(in)    :: Cell
      Type (Zone_Axis_type),         intent(in out):: Zone_Axis
      real(kind=cp), dimension(3),   intent(in)    :: h1
      real(kind=cp),                 intent(in)    :: omega
      real(kind=cp), dimension(3,3)                :: UB
      ! Local variables
      integer                        :: ierr
      real(kind=cp)                  :: theta1,theta2,alpha,del_omega,d1s,d2s
      real(kind=cp), dimension(3)    :: h2,ho1,ho2,s1,s2
      real(kind=cp), dimension(3,3)  :: Rot

      call clear_error()

      !First check that the provided reflection is perpendicular to uvw
      if(dot_product(real(Zone_Axis%uvw),h1) > 0.01_cp) then
        Err_CFML%Ierr=1
        write(unit=Err_CFML%Msg,fmt="(2(a,f8.3))") "The given reflection: ",h1," is not perpendicular to: ",real(Zone_Axis%uvw)
        return
      end if
      zone_axis = Get_basis_from_uvw(1.0_cp,Zone_Axis%uvw,Cell) !Here we use a dmin=1.0 angstrom
      if(Err_CFML%Ierr /= 0) then
        Err_CFML%Msg = "Error in the calculation of reflection plane: "//trim(Err_CFML%Msg )
        return
      end if
      h2=real(zone_axis%rx)       !Second reflection in the plane
      if(Co_linear(h1,h2,3)) h2=real(zone_axis%ry)
      !
      !Determination of the Bragg angles of two reflections in the scattering plane
      !
      d1s=sqrt(dot_product(h1,matmul(Cell%GR,h1))) !d1*
      d2s=sqrt(dot_product(h2,matmul(Cell%GR,h2))) !d2*
      theta1=asind(0.5*wave*d1s)
      theta2=asind(0.5*wave*d2s)
      alpha=acosd(dot_product(h1,matmul(Cell%GR,h2))/d1s/d2s) !Angle between reciprocal vectors
      del_omega=theta2-theta1+alpha  !Variation in omega to put the second reflection in diffraction position
      !
      !Calculation of the Cartesian components of the two reflections in the scattering plane
      !
      s1=d1s*(/cosd(Theta1),-sind(Theta1),0.0_cp/)
      Rot= Phi_Mat(omega)
      ho1=matmul(transpose(rot),s1)
      s2=d2s*(/cosd(Theta2),-sind(Theta2),0.0_cp/)
      Rot= Phi_mat(omega+del_omega)
      ho2=matmul(transpose(rot),s2)
      !
      ! Generate UB-matrix
      !
      UB= GenUB(Cell%BL_M,h1,h2,ho1,ho2)
      if(Err_CFML%Ierr /= 0) then
        Err_CFML%Msg = "Error in the calculation of UB-matrix: "//trim(Err_CFML%Msg)
      end if
    End Function Get_UB_from_uvw_hkl_omega


    !!----
    !!---- Module Subroutine normal(v,ierr)
    !!----    real(kind=cp), Intent(In Out), Dimension(3)   :: v
    !!----    Integer, Intent(Out)                          :: ierr
    !!----
    !!----    Normalise vector V (in Cartesian components)
    !!----
    !!---- Update: April 2008
    !!
    Module Subroutine normal(v)
       !---- Argument ----!
       real(kind=cp), Intent(In Out), Dimension(3)   :: v

       !--- Local Variables ---!
       real(kind=cp) :: d

       d=Dot_Product(v,v)
       If (d <= 0.0_cp) Then
          Err_CFML%Ierr=-1
       Else
          ierr=0
          d=Sqrt(d)
          v=v/d
       End If
    End Subroutine normal

    !!----
    !!---- Module Subroutine refvec(vhkl,ub,vs,vz)
    !!----    real(kind=cp), Intent(In),  Dimension(3)    :: vhkl
    !!----    real(kind=cp), Intent(In),  Dimension(3,3)  :: ub
    !!----    real(kind=cp), Intent(Out), Dimension(3)    :: vs,vz
    !!----
    !!----    Calculate vs,vz as reference vectors for defining Psi=0
    !!----    The B-L convention is that Psi=0 when the reflection hkl is
    !!----    in diffraction position and the c* is in the plane defined
    !!----    by vhkl and vz (z-axis of the laboratory system) for all
    !!----    reflections except when vhkl is parallel to c* in which case
    !!----    the vector b* plays the role of c* in the above prescription.
    !!----    The vector vhkl is provided with components in the reciprocal
    !!----    lattice.
    !!----
    !!---- Update: July 2008, June 2020
    !!
    Module Subroutine refvec(vhkl,ub,vs,vz)
       !---- Arguments ----!
       real(kind=cp), Intent(In),  Dimension(3)    :: vhkl
       real(kind=cp), Intent(In),  Dimension(3,3)  :: ub
       real(kind=cp), Intent(Out), Dimension(3)    :: vs,vz

       !--- Local Variables ---!
       real(kind=cp),Dimension(3) :: hn,h0
       real(kind=cp),Dimension(3) :: h1=(/0.0_cp,0.0_cp,1.0_cp/),h2=(/0.0_cp,1.0_cp,0.0_cp/),v0=(/0.0_cp,0.0_cp,1.0_cp/)

       !---- Test if VHKL is parallel to H1
       hn=vhkl
       !Check that the vector is non-null
       if(sum(abs(hn)) < 0.0001_cp) then
         Err_CFML%Ierr=-1
         Err_CFML%Msg = " Error: Null input vector @ refvec"
         return
       else
         ierr=0
       end if
       h0=Cross_Product(hn,h1)
       If (Sum(Abs(h0)) > 0.0001_cp) Then
          h0=h1         !vhkl IS NOT parallel to c*(=h1), so h1 can be used as reference
       Else
          h0=h2         !vhkl IS parallel to c*(=h1), so h2=b* is used as reference
       End If
       vs=Matmul(ub,h0) !Put the reciprocal vector c* (or b*) in the laboratory system
       vz=v0            !Z-xis of the laboratory system
    End Subroutine refvec

    !!----
    !!---- Module Subroutine triple(v1,v2,tv)
    !!----    real(kind=cp),    Intent(In Out), Dimension(3)  :: v1,v2
    !!----    real(kind=cp),    Intent(Out),    Dimension(3,3):: tv
    !!----
    !!----    Construct orthonormal triplet matrix TV, with column vectors :
    !!----    V1, (V1 x V2) x V1, V1 x V2.
    !!----
    !!---- Update: July 2008, June 2020
    !!
    Module Subroutine triple(v1,v2,tv)
       !---- Arguments ----!
       real(kind=cp),    Intent(In Out), Dimension(3)  :: v1,v2  !they come back normalized and V2 perp. to V1
       real(kind=cp),    Intent(Out),    Dimension(3,3):: tv

       !--- Local Variables ---!
       real(kind=cp),Dimension(3) :: v3
       call clear_error()
       Call normal(v1)
       v3=Cross_Product(v1,v2)
       Call normal(v3)
       v2=Cross_Product(v3,v1)
       Call normal(v2)
       If(Err_CFML%Ierr /= 0) Return
       tv(:,1)=v1(:)
       tv(:,2)=v2(:)
       tv(:,3)=v3(:)
    End Subroutine triple

 End SubModule SXTAL_UB
