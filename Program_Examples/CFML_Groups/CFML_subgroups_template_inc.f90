! Template for constructing subgroups
       !---- Construct first the generators of centring translations ----!
       ng=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings
       if (SpG%centred /= 1) nop=nop*2
       if(SpG%num_lat > 0) then
         do i=1,SpG%num_lat
            ng=ng+1
            gen(ng)= SpG%Symb_Op(1+nop*i)
         !   Op_lat(ng)= SpG%Op(1+nop*i))
         end do
       end if
       nla=ng
       mp=SpG%num_lat+1 !<-- Index for passing generators (fixed!)
       nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
       L=0
       !---- Determine first the triclinic subgroups
       cen_added=.false.
       ng=mp
       do nla=0,SpG%num_lat
         do
             L=L+1
             newg=.true.
             call Group_Constructor(gen(mp-nla:ng),SubG(L))
             do j=1,L-1
                if (SubG(L) == SubG(j)) then
                   newg=.false.
                   exit
                end if
             end do
             if (.not. newg) L=L-1
             if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
                ng=ng+1
                gen(ng)=SpG%Symb_Op(nc)
                cen_added=.true.
             else
                exit
             end if
         end do
       end do

       !---- Determine first the groups with only one rotational generator
       ng=mp+1
       do nla=0,SpG%num_lat
         do i=2,nop
            gen(ng) = SpG%Symb_Op(i)
            cen_added=.false.
            do
               L=L+1
               if (L > maxg) then
                  nsg=maxg
                  return
               end if
               newg=.true.
               write(*,*) " Section of gen:",mp-nla,ng
               call Group_Constructor(gen(mp-nla:ng),SubG(L))
               do j=1,L-1
                 if (SubG(L) == SubG(j)) then
                    newg=.false.
                    exit
                 end if
               end do
               if (.not. newg) L=L-1
               if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
                  ng=ng+1
                  gen(ng)=SpG%Symb_Op(nc)
                  cen_added=.true.
               else
                  exit
               end if
            end do
         end do
       end do

       !---- Determine now the groups with two rotational generator ----!

      !ng=mp+2
      !
      !do nla=0,SpG%num_lat
      !  do i1=2,nop-1
      !     gen(ng-1) = SpG%Symb_Op(i1)
      !     do i2 = i1+1,nop
      !        gen(ng) = SpG%Symb_Op(i2)
      !        !ng=nla+2
      !        cen_added=.false.
      !        do
      !           L=L+1
      !           if (L > maxg) then
      !              nsg=maxg
      !              return
      !           end if
      !           newg=.true.
      !           call Group_Constructor(gen(mp-nla:ng),SubG(L))
      !           do j=1,L-1
      !              if (SubG(L) == SubG(j)) then
      !                 newg=.false.
      !                 exit
      !              end if
      !           end do
      !           if (.not. newg) L=L-1
      !           if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
      !              ng=ng+1
      !              gen(ng)=SpG%Symb_Op(nc)
      !              cen_added=.true.
      !           else
      !              exit
      !           end if
      !        end do
      !     end do
      !  end do
      !end do

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
