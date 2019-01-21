! Template for constructing subgroups
       !---- Construct first the generators of centring translations ----!
       ng=0; nc=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings & centre of symmetry
       if (SpG%centred /= 1) then
          nop=nop*2 !!number of symmetry operators excluding lattice centrings
          nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
          gen_cent=SpG%Symb_Op(nc)
          call Allocate_Operator(SpG%d,Op_cent)
          Op_cent=SpG%Op(nc)  !Operator corresponding to the centre of symmetry
       end if
       nla=0
       if(SpG%num_lat > 0) then
         do i=1,SpG%num_lat
            ng=ng+1
            gen_lat(ng)= SpG%Symb_Op(1+nop*i)
            Op_lat(ng)= SpG%Op(1+nop*i)  !Operators corresponding to the lattice centring vectors
         end do
         nla=ng
       end if
       !First work with the Numops operators to determine the subgroups, the other subgroups
       !will be obtained adding progressively the rest of generators (centre of symmetry and
       !lattice centrings.
       L=0
       !---- Determine first the groups with three rotational generators
       if(nop > 3) then
         ng=3
         do i=2, SpG%numops    !SpG%numops-1
            gen(1) = SpG%Symb_Op(i)
            do j=i+1,nop-1  !SpG%numops
              gen(2)=SpG%Symb_Op(j)
              do m=j+1,nop
                gen(3)=SpG%Symb_Op(m)
                L=L+1
                if (L > maxg) then
                   nsg=maxg
                   return
                end if
                newg=.true.; Err_group=.false.
                !write(*,"(a,i5,a)") "  Group # ",L, "  "//trim(gen(1))//"  "//trim(gen(2))//"  "//trim(gen(3))
                call Group_Constructor(gen(1:ng),SubG(L))
                index_sg(L)=SpG%multip/SubG(L)%multip
                if(index_sg(L) == 1) then
                   L=L-1
                   cycle
                end if
                if(present(indexg)) then
                   if(index_sg(L) /= indexg) then
                     L=L-1
                     cycle
                   end if
                end if
                if(.not. Err_group) then
                  do k=L-1,1,-1
                    if(index_sg(L) /= index_sg(k)) cycle
                    if (SubG(L) == SubG(k)) then
                       newg=.false.
                       exit
                    end if
                  end do
                  if (.not. newg) then
                    L=L-1
                  end if
                else
                  L=L-1
                end if
              end do
            end do
         end do
       end if
       ns_3=L
       !---- Determine now the groups with two rotational generators
       !if(SpG%numops > 2) then
       if(nop > 2) then
         ng=2
         do i=2, SpG%numops    !SpG%numops-1
            gen(1) = SpG%Symb_Op(i)
            do j=i+1,nop   !SpG%numops
              gen(2)=SpG%Symb_Op(j)
              L=L+1
              if (L > maxg) then
                 nsg=maxg
                 return
              end if
              newg=.true.; Err_group=.false.
              call Group_Constructor(gen(1:ng),SubG(L))
              index_sg(L)=SpG%multip/SubG(L)%multip
              if(index_sg(L) == 1) then
                 L=L-1
                 cycle
              end if
              if(present(indexg)) then
                 if(index_sg(L) /= indexg) then
                   L=L-1
                   cycle
                 end if
              end if
              if(.not. Err_group) then
                do k=L-1,ns_3,-1
                  if(index_sg(L) /= index_sg(k)) cycle
                  if (SubG(L) == SubG(k)) then
                     newg=.false.
                     exit
                  end if
                end do
                if (.not. newg) L=L-1
              else
                L=L-1
              end if
            end do
         end do
       end if
       ns_2=L-ns_3
       nsg=L
       !---- Determine finally the groups with only one rotational generator (cyclic groups)
       ng=1
       do i=2,nop !SpG%numops
          gen(1) = SpG%Symb_Op(i)
          L=L+1
          if (L > maxg) then
             nsg=maxg
             return
          end if
          newg=.true.; Err_group=.false.
          !write(*,*) (trim(gen(j))//" ; ",j=1,ng)
          call Group_Constructor(gen(1:ng),SubG(L))
          index_sg(L)=SpG%multip/SubG(L)%multip
          if(index_sg(L) == 1) then
             L=L-1
             cycle
          end if
          if(present(indexg)) then
             if(index_sg(L) /= indexg) then
               L=L-1
               cycle
             end if
          end if
          if(.not. Err_group) then
             do k=L-1,nsg,-1
               if(index_sg(L) /= index_sg(k)) cycle
               if (SubG(L) == SubG(k)) then
                  newg=.false.
                  exit
               end if
             end do
             if (.not. newg) L=L-1
          else
             L=L-1
          end if
       end do
       nsg=L
       n_nc_group=L

       !---- Determine now the new groups adding a centre of symmetry (without lattice centring) if it exists
       if (SpG%centred /= 1) then !This doubles the number of groups
         do i=1,n_nc_group
           if(SubG(i)%centred == 0 .or. SubG(i)%centred == 2) cycle
           L=L+1
           if (L > maxg) then
              nsg=maxg
              return
           end if
           !List of generators for the new group
           aux_string=adjustl(SubG(i)%generators_list//";"//trim(gen_cent))
           newg=.true.; Err_group=.false.
           !write(*,"(a,i5,a)") "  Group # ",L,"  "//trim(aux_string)
           call Group_Constructor(aux_string,SubG(L))
           index_sg(L)=SpG%multip/SubG(L)%multip
           if(index_sg(L) == 1) then
              L=L-1
              cycle
           end if
           if(present(indexg)) then
              if(index_sg(L) /= indexg) then
                L=L-1
                cycle
              end if
           end if
           if(.not. Err_group) then
             do k=L-1,n_nc_group,-1
               if(index_sg(L) /= index_sg(k)) cycle
               if (SubG(L) == SubG(k)) then
                  newg=.false.
                  exit
               end if
             end do
             if (.not. newg) L=L-1
           else
             L=L-1
           end if
         end do
       end if
       n_nc_group=L

       !Determine now the rest of groups adding the lattice translations if they exist in the
       !original space group

       if(SpG%num_lat > 0) then
         do j=1,nla
             do i=1,n_nc_group
               L=L+1
               if (L > maxg) then
                  nsg=maxg
                  return
               end if
               aux_string=adjustl(SubG(i)%generators_list//";"//trim(gen_lat(j)))
               newg=.true.; Err_group=.false.
               call Group_Constructor(aux_string,SubG(L))
               index_sg(L)=SpG%multip/SubG(L)%multip
               if(index_sg(L) == 1) then
                  L=L-1
                  cycle
               end if
               if(present(indexg)) then
                  if(index_sg(L) /= indexg) then
                    L=L-1
                    cycle
                  end if
               end if
               if(.not. Err_group) then
                 do k=L-1,n_nc_group,-1
                   if(index_sg(L) /= index_sg(k)) cycle
                   if (SubG(L) == SubG(k)) then
                      newg=.false.
                      exit
                   end if
                 end do
                 if (.not. newg) L=L-1
               else
                 L=L-1
               end if
             end do
         end do
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