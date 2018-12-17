! Template for constructing subgroups
       !---- Construct first the generators of centring translations ----!
       ng=0; nc=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings & centre of symmetry
       if (SpG%centred /= 1) then
          nop=nop*2 !!number of symmetry operators excluding lattice centrings
          nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
          gen_cent=SpG%Symb_Op(nc)
          call Allocate_Operator(SpG%d,Op_cent)
          Op_cent=SpG%Op(nc)
       end if
       if(SpG%num_lat > 0) then
         do i=1,SpG%num_lat
            ng=ng+1
            gen_lat(ng)= SpG%Symb_Op(1+nop*i)
            Op_lat(ng)= SpG%Op(1+nop*i)
         end do
       end if
       !First work with the Numops operators to determine the subgroups, the other subgroups
       !will be obtained adding progressively the rest of generators (centre of symmetry and
       !lattice centrings.
       L=0
       !---- Determine first the groups with only one rotational generator
       ng=1
       do i=2,SpG%numops
          gen(1) = SpG%Symb_Op(i)
          L=L+1
          if (L > maxg) then
             nsg=maxg
             return
          end if
          newg=.true.
          !write(*,*) (trim(gen(j))//" ; ",j=1,ng)
          call Group_Constructor(gen(1:ng),SubG(L))
          do k=1,L-1
            if (SubG(L) == SubG(k)) then
               newg=.false.
               exit
            end if
          end do
          if (.not. newg) L=L-1
       end do
       ns_1=L
       !---- Determine now the groups with two rotational generators
       if(SpG%numops > 2) then
         ng=2
         do i=2,SpG%numops-1
            gen(1) = SpG%Symb_Op(i)
            do j=i+1,SpG%numops
              gen(2)=SpG%Symb_Op(j)
              L=L+1
              if (L > maxg) then
                 nsg=maxg
                 return
              end if
              newg=.true.
              call Group_Constructor(gen(1:ng),SubG(L))
              do k=1,L-1
                if (SubG(L) == SubG(k)) then
                   newg=.false.
                   exit
                end if
              end do
              if (.not. newg) L=L-1
            end do
         end do
         ns_2=L-ns_1
       end if
       nsg=L
       n_nc_group=L
       !write(*,*) " Number of subgroups of the first Numops elements: ",n_nc_group

       !---- Determine now the new groups adding a centre of symmetry (without lattice centring) if it exists
       if (SpG%centred /= 1) then !This doubles the number of groups
         do i=1,n_nc_group
           L=L+1
           call Allocate_Group(SpG%d,2*SubG(i)%multip,SubG(L))
           if(allocated(SubG(L)%centre_coord)) deallocate(SubG(L)%centre_coord)
           allocate(SubG(L)%centre_coord(d-1))
           if(SubG(i)%num_alat /= 0) then
               if(allocated(SubG(L)%aLat_tr)) deallocate(SubG(L)%aLat_tr)
               allocate(SubG(L)%aLat_tr(d-1,SubG(i)%num_alat))
               SubG(L)%aLat_tr=SubG(i)%aLat_tr
               SubG(L)%num_alat=SubG(i)%num_alat
           end if

           k=SubG(i)%numops
           do j=1,SubG(i)%numops
             SubG(L)%Op(j)=SubG(i)%Op(j)
             SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
             k=k+1
             SubG(L)%Op(k)=SubG(i)%Op(j)*Op_cent
             SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
           end do
           SubG(L)%Numops=SubG(i)%Numops
           SubG(L)%mag_type=SubG(i)%mag_type
           SubG(L)%centred=SpG%centred
           SubG(L)%centre_coord=SpG%centre_coord
           SubG(L)%Multip= 2*SubG(i)%multip
         end do
       end if
       nsg=L
       n_nc_group=L
       !if (SpG%centred /= 1) write(*,*) " Number of subgroups of adding a centre of symmetry: ",n_nc_group

       !Determine now the rest of groups adding the lattice translations if they exist in the
       !original space group
       if(SpG%num_lat > 0) then
         Select Case (SpG%num_lat)
           case(1)
             do i=1,n_nc_group
               L=L+1
               call Allocate_Group(SpG%d,2*SubG(i)%multip,SubG(L))
               if(allocated(SubG(L)%centre_coord)) deallocate(SubG(L)%centre_coord)
               allocate(SubG(L)%centre_coord(d-1))
               if(allocated(SubG(L)%Lat_tr)) deallocate(SubG(L)%Lat_tr)
               allocate(SubG(L)%Lat_tr(d-1,1))
               SubG(L)%Lat_tr(:,1)=SpG%Lat_tr(:,1)
               if(SubG(i)%num_alat /= 0) then
                   if(allocated(SubG(L)%aLat_tr)) deallocate(SubG(L)%aLat_tr)
                   allocate(SubG(L)%aLat_tr(d-1,SubG(i)%num_alat))
                   SubG(L)%aLat_tr=SubG(i)%aLat_tr
                   SubG(L)%num_alat=SubG(i)%num_alat
               end if
               SubG(L)%num_lat=1
               k=SubG(i)%numops*cent(SubG(i)%centred)
               kb=k
               do j=1,kb
                 SubG(L)%Op(j)=SubG(i)%Op(j)
                 SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(1)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               SubG(L)%Multip=2*SubG(i)%Multip
               SubG(L)%Numops=SubG(i)%Numops
               SubG(L)%mag_type=SubG(i)%mag_type
               SubG(L)%centred=SubG(i)%centred
               if(SubG(L)%centred /= 1) SubG(L)%centre_coord=SpG%centre_coord
             end do
           case(2)
             do i=1,n_nc_group
               L=L+1
               call Allocate_Group(SpG%d,3*SubG(i)%multip,SubG(L))
               if(allocated(SubG(L)%centre_coord)) deallocate(SubG(L)%centre_coord)
               allocate(SubG(L)%centre_coord(d-1))
               if(allocated(SubG(L)%Lat_tr)) deallocate(SubG(L)%Lat_tr)
               allocate(SubG(L)%Lat_tr(d-1,2))
               SubG(L)%Lat_tr(:,:)=SpG%Lat_tr(:,:)
               SubG(L)%num_lat=2
               if(SubG(i)%num_alat /= 0) then
                   if(allocated(SubG(L)%aLat_tr)) deallocate(SubG(L)%aLat_tr)
                   allocate(SubG(L)%aLat_tr(d-1,SubG(i)%num_alat))
                   SubG(L)%aLat_tr=SubG(i)%aLat_tr
                   SubG(L)%num_alat=SubG(i)%num_alat
               end if
               k=SubG(i)%numops*cent(SubG(i)%centred)
               kb=k
               do j=1,kb
                 SubG(L)%Op(j)=SubG(i)%Op(j)
                 SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(1)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               do j=1,kb
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(2)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               SubG(L)%Multip=3*SubG(i)%Multip
               SubG(L)%Numops=SubG(i)%Numops
               SubG(L)%mag_type=SubG(i)%mag_type
               SubG(L)%centred=SubG(i)%centred
               if(SubG(L)%centred /= 1) SubG(L)%centre_coord=SpG%centre_coord
             end do
           case(3)
             do i=1,n_nc_group
               L=L+1
               call Allocate_Group(SpG%d,4*SubG(i)%multip,SubG(L))
               if(allocated(SubG(L)%centre_coord)) deallocate(SubG(L)%centre_coord)
               allocate(SubG(L)%centre_coord(d-1))
               if(allocated(SubG(L)%Lat_tr)) deallocate(SubG(L)%Lat_tr)
               allocate(SubG(L)%Lat_tr(d-1,3))
               SubG(L)%Lat_tr(:,:)=SpG%Lat_tr(:,:)
               SubG(L)%num_lat=3
               if(SubG(i)%num_alat /= 0) then
                   if(allocated(SubG(L)%aLat_tr)) deallocate(SubG(L)%aLat_tr)
                   allocate(SubG(L)%aLat_tr(d-1,SubG(i)%num_alat))
                   SubG(L)%aLat_tr=SubG(i)%aLat_tr
                   SubG(L)%num_alat=SubG(i)%num_alat
               end if
               k=SubG(i)%numops*cent(SubG(i)%centred)
               kb=k
               do j=1,kb
                 SubG(L)%Op(j)=SubG(i)%Op(j)
                 SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(1)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               do j=1,kb
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(2)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               do j=1,kb
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(3)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               SubG(L)%Numops=SubG(i)%Numops
               SubG(L)%mag_type=SubG(i)%mag_type
               SubG(L)%centred=SubG(i)%centred
               if(SubG(L)%centred /= 1) SubG(L)%centre_coord=SpG%centre_coord
               SubG(L)%Multip=4*SubG(i)%Numops
             end do
         End Select

       end if

       !Include now the case in which the original group is paramagnetic
       !if(SpG%mag_type == 2) then
       !  nsg=L
       !  do i=1,nsg
       !  end do
       !end if

       nsg=L
       if(present(point)) then
         point=.false.
         do j=1,nsg
           L=1
           do i=1,SpG%multip
              do k=L,SubG(j)%multip
               if(SubG(j)%Symb_Op(k) == SpG%Symb_Op(i)) then
                  point(i,j) = .true.
                  L=k+1
                  exit
               end if
              end do
           end do
         end do
       end if