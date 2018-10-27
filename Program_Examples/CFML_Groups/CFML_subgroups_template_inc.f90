! Template for constructing subgroups
       !---- Construct first the generators of centring translations ----!
       ng=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings
       if (SpG%centred /= 1) nop=nop*2
       do i=2,SpG%num_lat
          ng=ng+1
          gen(ng)= SpG%Symb_Op(1+nop*(i-1))
       end do

       nla=ng
       nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
       L=0
       !---- Determine first the triclinic subgroups
       cen_added=.false.
       do
           L=L+1
           newg=.true.
           call Group_Constructor(gen(1:ng),SubG(L))
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

       !---- Determine first the groups with only one rotational generator
       do i=2,nop
          ng=nla+1
          gen(ng) = SpG%Symb_Op(i)
          cen_added=.false.
          do
             L=L+1
             if (L > maxg) then
                nsg=maxg
                return
             end if
             newg=.true.
             write(*,*) (trim(gen(j))//" ; ",j=1,ng)
             call Group_Constructor(gen(1:ng),SubG(L))
             do j=1,L-1
               write(*,"(a,4i6)") " j, L,mulj,mulL:",j,l,SubG(j)%multip,SubG(L)%multip
               call print_group(SubG(j))
               call print_group(SubG(L))
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

       !---- Determine now the groups with two rotational generator ----!

       do i1=2,nop-1
          gen(nla+1) = SpG%Symb_Op(i1)
          do i2 = i1+1,nop
             gen(nla+2) = SpG%Symb_Op(i2)
             ng=nla+2
             cen_added=.false.
             do
                L=L+1
                if (L > maxg) then
                   nsg=maxg
                   return
                end if
                newg=.true.
                call Group_Constructor(gen(1:ng),SubG(L))
                if(mod(nop,SubG(L)%Numops) /= 0 .or. SubG(L)%multip == 0) then
                  L=L-1
                  newg=.false.
                else
                  do j=1,L-1
                     if (SubG(L) == SubG(j)) then
                        newg=.false.
                        exit
                     end if
                  end do
                  if (.not. newg) L=L-1
                end if
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
