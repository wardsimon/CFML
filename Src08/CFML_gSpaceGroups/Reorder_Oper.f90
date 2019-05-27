!!----
!!----
!!----
!!
SubModule (CFML_gSpaceGroups) SPG_009
   Contains
   
   !!----
   !!---- REORDER_OPERATORS
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Reorder_Operators(Multip, Op, Centred, Centre_Coord, Anticentred, &
                                       Anticentre_Coord, Numops, Num_Lat, Num_Alat,    &
                                       Lat_Tr, Alat_Tr, Mag_Type)
      !---- Arguments ----!
      integer,                            intent(in)     :: multip
      type(Symm_Oper_Type), dimension(:), intent(in out) :: Op
      integer,                            intent(out)    :: num_lat,num_alat,Numops, centred, anticentred, mag_type
      type(rational),dimension(:,:),      intent(out)    :: Lat_tr
      type(rational),dimension(:,:),      intent(out)    :: aLat_tr
      type(rational),dimension(:),        intent(out)    :: centre_coord
      type(rational),dimension(:),        intent(out)    :: anticentre_coord
      
      !--- Local variables ---!
      integer                          :: i,j,L,n,m,Ng,invt,i_centre,d
      real(kind=cp), dimension(multip) :: tr   !Sum of absolute values of Translations components associated to the array of operators
      logical,       dimension(multip) :: nul  !Logical to control the exclusion of an operator from the independet list
      type(rational), dimension(:,:),allocatable:: identity,invers,imat !(d,d)
      type(Symm_Oper_Type), dimension(multip)  :: Opr,Op_Lat,Op_aLat
      type(Symm_Oper_Type)                     :: Op_aux,Op_aux1,Op_aux2, &
                                                  Op_centre,Op_identp
      real(kind=cp), parameter :: loc_eps=0.001_cp
      real(kind=cp)            :: tmin
      type(rational)           :: ZERO, ONE, ONE_HALF
      
      type(rational), dimension(:), allocatable :: atr,nulo 

      !> Init
      call clear_error()
      ZERO=0//1;  ONE=1//1; ONE_HALF=1//2

      !> dimension of the full square matrices
      n=size(Op(1)%Mat,1)  
      
      !> dimension of the square matrix containing rotational operator
      d=n-1
      allocate(identity(d,d),invers(d,d),imat(d,d))
      call Allocate_Symm_Op(n,Op_identp)   ! {1|0}'
      
      identity=ZERO; nul=.false.; mag_type=1
      do i=1,d
         Op_identp%Mat(i,i)=ONE
         identity(i,i)=ONE
      end do
      Op_identp%time_inv=-1
      invers=-identity 
      centred=1 
      anticentred=1

      !> Insertion sort putting the negative determinants at the bottom
      call Sort_Oper(multip,Op(1:multip),"det")
      do i=1,Multip
         tr(i)=sum(abs(Op(i)%Mat(1:d,n)))
         call Allocate_Symm_Op(n,Opr(i))
      end do

      !> Check if the group is paramagnetic
      j=0
      do i=2,Multip
         if (Op(i)== Op_identp) then
            j=i
            exit
         end if
      end do
      
      if (j /= 0) Then
         if (Op(j)%time_inv < 0) then
            do i=2,Multip   !Nullify all primed operators
               if (Op(i)%time_inv < 0) nul(i) = .true.
            end do
            mag_type=2
         end if
      end if

      !> Look for centre of symmetry, and centring translations
      num_lat=0; num_alat=0; tmin=1.0e8; i_centre=0
      do j=2,Multip
         if (nul(j)) cycle
         invt= Op(j)%time_inv
         imat=Op(j)%Mat(1:d,1:d)
         if (rational_equal(identity,imat) .and. invt == 1) then
            num_lat=num_lat+1
            Lat_tr(:,num_lat)=Op(j)%Mat(1:d,n)
            Op_Lat(num_lat)=Op(j)
            nul(j)=.true.   !Nullify centring translations
            cycle
         end if
         if (rational_equal(imat,invers) .and. invt == 1 ) then
            nul(j) = .true.
            if (tr(j) < tmin) Then
               tmin= tr(j)
               i_centre=j
            end if
         end if
      end do

      if (i_centre /= 0) then
         Op_centre=Op(i_centre)
         centre_coord=ONE_HALF*Op(i_centre)%Mat(1:d,n)
         if (tr(i_centre) < loc_eps) then
            centred = 2
         else
            centred = 0
         end if
      end if

      !> Nullify operators deduced by lattice translations and centre of symmetry
      do j=2,Multip-1   
         if (nul(j)) cycle

         if (mag_type ==2) then
            Op_aux=Op(j)*Op_identp
            do i=j+1,Multip
               if (nul(i)) cycle
               if (Op_aux == Op(i)) then
                  nul(i)=.true.
               end if
            end do
         end if

         if (num_lat > 0 .and. i_centre /= 0) then
            Op_aux=Op(j)*Op_centre
            do L=1,num_lat
               Op_aux1=Op(j)*Op_Lat(L)
               Op_aux2=Op_aux1*Op_centre
               do i=j+1,Multip
                  if (nul(i)) cycle
                  if (Op_aux == Op(i) .or. Op_aux1 == Op(i) .or. Op_aux2 == Op(i)) then
                     nul(i)=.true.
                  end if
               end do
            end do
            cycle
         end if

         if (num_lat > 0) then
            do L=1,num_lat
               Op_aux=Op(j)*Op_Lat(L)
               do i=j+1,Multip
                  if (nul(i)) cycle
                  if (Op_aux == Op(i)) then
                     nul(i)=.true.
                  end if
               end do
            end do
         end if

         if (i_centre /= 0) then
            Op_aux=Op(j)*Op_centre
            do i=j+1,Multip
               if (nul(i)) cycle
               if (Op_aux == Op(i) ) then
                  nul(i)=.true.
               end if
            end do
         end if
      end do

      !> Determine the lattice anti-translations
      if (allocated(atr)) deallocate(atr, nulo)
      allocate(atr(d), nulo(d))
      nulo=zero
      
      do_ext: do j=2,Multip
         !if(nul(j)) cycle
         invt= Op(j)%time_inv
         imat=Op(j)%Mat(1:d,1:d)
         atr=Op(j)%Mat(1:d,n)
         !> Modification JGP
         !if (rational_equal(identity,imat) .and. invt == -1) then 
             !if (rational_equal(atr,nulo)) cycle
         if (rational_equal(identity,imat) .and. invt == -1 .and. mag_type /= 2) then
            num_alat=num_alat+1
            aLat_tr(:,num_alat)=atr
            Op_aLat(num_alat)=Op(j)
         end if
      end do  do_ext !j=2,Multip

      if (num_alat > 0) then
         if (mag_type /= 2) then
            mag_type=4
         end if
      else
         if (mag_type /= 2) then
            if (any(Op(:)%time_inv < 0)) mag_type=3
         end if
      end if

      !> Determine the reduced set of symmetry operators"
      j=0
      do i=1,Multip
         if (nul(i)) cycle
         j=j+1
         Opr(j) = Op(i)
      end do
      Numops=j

      !> Promote the reduced set of symmetry operators to the top of the list
      Op(1:j)=Opr(1:j)

      !> Reorder the reduced set putting primed elements at the bottom
      call Sort_oper(Numops,Op(1:Numops),"tim")

      m=Numops*(num_lat+1)*cent(centred)
      if (mag_type == 2) m=m*2
      if ( m /= Multip) then !Check that it is OK
         Err_CFML%Ierr = 1
         write(unit=Err_CFML%Msg,fmt="(2(a,i4))") "Reorder_Operators@SPACEG: Warning in Multip=",Multip, " Calculated Multip: ",m
         return
      end if

      !> Re-Construct, in an ordered way, all the symmetry operators
      !> starting with the reduced set
      m=Numops
      ng=m
      if (i_centre /= 0) then   !First apply the centre of symmetry
         do i=1,Numops
            m=m+1
            Op(m) = Op(i) * Op_centre
         end do
      end if

      ng=m  ! Number or symmetry operators including centre of symmetry
      if (Num_Lat > 0) then  !Fourth apply the lattice centring translations
         do L=1,Num_Lat
            do i=1,ng
               m=m+1
               Op(m)=Op_Lat(L)*Op(i)
            end do
         end do
      end if

      if (mag_type == 2) then
         ng=m
         do i=1,ng
            m=m+1
            Op(m)=Op_identp*Op(i)
         end do
      end if

      !> Normally here the number of operators should be equal to multiplicity
      !> Test that everything is OK
      ng=m
      if (ng /= Multip) then
         Err_CFML%Ierr = 1
         write(unit=Err_CFML%Msg,fmt="(2(a,i3))") "Reorder_Operators@SPACEG: Error in the multiplicity ",&
                                                  Multip," has not been recovered, value of ng=",ng
         return
      end if

      !> Determine if the group contain the primed inversion
      do i = 1 , ng
         if (rational_equal(Op(i)%Mat(1:d,1:d),invers) .and. Op(i)%time_inv == -1 ) then
            anticentre_coord(1:d)=ONE_HALF*Op(i)%Mat(1:d,n)
            if (tr(i) < loc_eps) then
               anticentred = 2
            else
               anticentred = 0
            end if
            exit
         end if
      end do
   End Subroutine Reorder_Operators

End SubModule SPG_009
   
