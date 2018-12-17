! Template algorithm for re-ordering an arbitrary group characterized by matrices of
! whatever kind and dimensions. The direct operations (det=1) are placed on top, followed
! by those related with the first block by the centre of symmetry if it exists and finally
! followed by the new block plus centring translations. In case of paramagnetic groups
!(the operation {1|0}' is present) the primed operation are put at the bottom.

      !Initializing
      n=size(Op(1)%Mat,1) !dimension of the full square matrices
      d=n-1               !here d is the dimension of the square matrix containing rotational operator
      allocate(identity(d,d),invers(d,d),imat(d,d))
      call Allocate_Operator(n,Op_identp)   ! {1|0}'
      identity=ZERO; nul=.false.; mag_type=1
      do i=1,d
        Op_identp%Mat(i,i)=ONE
        identity(i,i)=ONE
      end do
      Op_identp%time_inv=-1
      invers=-identity !Inversion
      centred=1 !Default value for non-centrosymmetric groups

      !Insertion sort putting the negative determinants at the bottom
      call sort_Op(multip,Op(1:multip),"det")
      do i=1,Multip
        tr(i)=sum(abs(Op(i)%Mat(1:d,n)))
        call Allocate_Operator(n,Opr(i))
      end do


      !Testing
      !write(*,*) " "
      !write(*,*) "List of operators after re-ordering by determinant: "
      !m=0
      !do i=1,Multip
      !  if(nul(i)) cycle
      !  m=m+1
      !  write(*,"(2i5,a30,f12.5,i6,tr2,L)") m,i,trim(Symbol_Operator(Op(i))),tr(i),Op(i)%dt,nul(i)
      !end do
      !end Testing


      !Check if the group is paramagnetic
      j=0
      do i=2,Multip
        if(Op(i)== Op_identp) then
          j=i
          exit
        end if
      end do
      if(j /= 0) Then
        if(Op(j)%time_inv < 0) then
          do i=2,Multip   !Nullify all primed operators
            if(Op(i)%time_inv < 0) nul(i) = .true.
          end do
          mag_type=2
        end if
      end if
      !----End intial re-ordering

      !Look for centre of symmetry, and centring translations
      num_lat=0; num_alat=0; tmin=1.0e8; i_centre=0
      do j=2,Multip
         if(nul(j)) cycle
         invt= Op(j)%time_inv
         imat=Op(j)%Mat(1:d,1:d)
         if(equal_matrix(identity,imat) .and. invt == 1) then
            num_lat=num_lat+1
            Lat_tr(:,num_lat)=Op(j)%Mat(1:d,n)
            Op_Lat(num_lat)=Op(j)
            nul(j)=.true.   !Nullify centring translations
            cycle
         end if
         if(equal_matrix(imat,invers) .and. invt == 1 ) then
             nul(j) = .true.
             if(tr(j) < tmin) Then
               tmin= tr(j)
               i_centre=j
             end if
         end if
      end do

      if(i_centre /= 0) then
         Op_centre=Op(i_centre)
         centre_coord=ONE_HALF*Op(i_centre)%Mat(1:d,n)
         if(tr(i_centre) < loc_eps) then
           centred = 2
         else
           centred = 0
         end if
      end if

      do j=2,Multip-1   !Nullify operators deduced by lattice translations and centre of symmetry
         if(nul(j)) cycle

           if(mag_type ==2) then
             Op_aux=Op(j)*Op_identp
             do i=j+1,Multip
                if(nul(i)) cycle
                if(Op_aux == Op(i)) then
                   nul(i)=.true.
                end if
             end do
           end if

           if(num_lat > 0 .and. i_centre /= 0) then
              Op_aux=Op(j)*Op_centre
              do L=1,num_lat
                 Op_aux1=Op(j)*Op_Lat(L)
                 Op_aux2=Op_aux1*Op_centre
                 do i=j+1,Multip
                    if(nul(i)) cycle
                    if(Op_aux == Op(i) .or. Op_aux1 == Op(i) .or. Op_aux2 == Op(i)) then
                       nul(i)=.true.
                    end if
                 end do
              end do
              cycle
           end if

           if(num_lat > 0) then
              do L=1,num_lat
                 Op_aux=Op(j)*Op_Lat(L)
                 do i=j+1,Multip
                    if(nul(i)) cycle
                    if(Op_aux == Op(i)) then
                       nul(i)=.true.
                    end if
                 end do
              end do
           end if

           if(i_centre /= 0) then
              Op_aux=Op(j)*Op_centre
              do i=j+1,Multip
                 if(nul(i)) cycle
                 if(Op_aux == Op(i) ) then
                    nul(i)=.true.
                 end if
              end do
           end if
      end do

      !Determine the lattice anti-translations
      do_ext: do j=2,Multip
        if(nul(j)) cycle
        invt= Op(j)%time_inv
        imat=Op(j)%Mat(1:d,1:d)
        if(equal_matrix(identity,imat) .and. invt == -1) then
          num_alat=num_alat+1
          aLat_tr(:,num_alat)=Op(j)%Mat(1:d,n)
          Op_aLat(num_alat)=Op(j)
        end if
      end do  do_ext !j=2,Multip


      if(num_alat > 0) then
        if(mag_type /= 2) then
          mag_type=4
        end if
      else
        if(mag_type /= 2) then
          if(any(Op(:)%time_inv < 0)) mag_type=3
        end if
      end if

      ! => Determine the reduced set of symmetry operators"
      j=0
      do i=1,Multip
        if(nul(i)) cycle
        j=j+1
        Opr(j) = Op(i)
      end do
      Numops=j

      !Promote the reduced set of symmetry operators to the top of the list
      Op(1:j)=Opr(1:j)

      !Reorder the reduced set putting primed elements at the bottom
      call sort_op(Numops,Op(1:Numops),"tim")

      !Testing
      !write(*,*) " "
      !write(*,*) "List of reduced set of operators after re-ordering by time inversion: "
      !m=0
      !do i=1,Numops
      !  if(nul(i)) cycle
      !  m=m+1
      !  write(*,"(2i5,a30,f12.5,i6,tr2,L)") m,i,trim(Symbol_Operator(Op(i))),tr(i),Op(i)%time_inv,nul(i)
      !end do
      !end Testing



      m=Numops*(num_lat+1)*cent(centred)
      if(mag_type == 2) m=m*2
      if( m /= Multip) then !Check that it is OK
        write(unit=Err_group_mess,fmt="(2(a,i4))") " Warning! Multip=",Multip, " Calculated Multip: ",m
        Err_group=.true.
        return
      end if

      !Re-Construct, in an ordered way, all the symmetry operators
      !starting with the reduced set
      m=Numops
      ng=m
      if(i_centre /= 0) then   !First apply the centre of symmetry
        do i=1,Numops
          m=m+1
          Op(m) = Op(i) * Op_centre
        end do
      end if

      ng=m  ! Number or symmetry operators including centre of symmetry

      if(Num_Lat > 0) then  !Fourth apply the lattice centring translations
        do L=1,Num_Lat
           do i=1,ng
             m=m+1
             Op(m)=Op_Lat(L)*Op(i)
           end do
        end do
      end if

      if(mag_type == 2) then
         ng=m
         do i=1,ng
           m=m+1
           Op(m)=Op_identp*Op(i)
         end do
      end if

      !Normally here the number of operators should be equal to multiplicity
      !Test that everything is OK
      ng=m
      if(ng /= Multip) then
        Err_group=.true.
        write(unit=Err_group_mess,fmt="(2(a,i3))") " => Problem! the multiplicity ",Multip," has not been recovered, value of ng=",ng
      end if
