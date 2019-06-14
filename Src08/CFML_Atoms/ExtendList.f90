!!----
!!----
!!----
!!
SubModule (CFML_Atoms) Atm_005
   Contains
   
   !!----
   !!---- EXTEND_LIST
   !!----
   !!----    Subroutine to generate atoms in the primitive (conven=.false.) or the conventional
   !!----    unit cell (conven=.true.), Excluding atoms with A%atom(:)%active=.false.
   !!----
   !!---- Update: February - 2005
   !!
   Module Subroutine Extend_List(A, B, Spg, Conven)
      !---- Arguments ----!
      class(alist_type),   intent(in)     :: A         ! Atom list (asymmetric unit)
      class(alist_type),   intent(in out) :: B         ! Atom list into the unit cell
      type(SpG_Type),      intent(in)     :: SpG       ! SpaceGroup
      logical, optional,   intent(in)     :: Conven    ! If present and .true. using the whole conventional unit cell

      !---- Local Variables ----!
      character(len=4)                      :: fmm
      integer                               :: k,j,l,nt,npeq,n
      real(kind=cp),dimension(3)            :: xo,xx
      real(kind=cp),dimension(3,Spg%multip) :: u
      logical                               :: ccell
      
      type(Atm_List_Type)      :: c1
      type(Atm_std_List_Type)  :: c2
      type(MAtm_std_List_Type) :: c3
      type(Atm_ref_List_Type)  :: c4

      !> Init
      call clear_error()
      ccell=.false.
      if (present(conven)) ccell=conven
      
      npeq=SpG%numops
      if (SpG%centred == 2) npeq=npeq*2
      if (ccell) npeq=SpG%multip
      
      n=0
      select type (A)
         type is (atm_list_type)
            call allocate_atom_list(npeq*A%natoms,c1)
            do k=1,A%natoms
               if (.not. A%active(k)) cycle
               l=1
               n=n+1
               c1%Atom(n)=A%Atom(k)
               xo= modulo_lat(A%atom(k)%x)
               u(:,l)= xo
               c1%Atom(n)%x=xo

               loop1:do j=2,npeq
                  xx=Apply_OP(SpG%Op(j),xo)
                  xx=modulo_lat(xx)
                  do nt=1,l
                     if (equal_vector(u(:,nt),xx,3)) then
                        c1%atom(n-(l-nt))%occ = c1%atom(n-(l-nt))%occ + A%atom(k)%occ
                        cycle loop1
                     end if
                  end do
                  l=l+1
                  u(:,l)=xx(:)
                  n=n+1
            
                  select case (l)
                     case(:9)
                        write(unit=fmm,fmt="(i1)") l
                     case(10:99)
                        write(unit=fmm,fmt="(i2)") l
                     case(100:999)
                        write(unit=fmm,fmt="(i3)") l
                  end select
                  c1%Atom(n)     = A%Atom(k)
                  c1%Atom(n)%lab = trim(A%Atom(k)%lab)//"_"//adjustl(fmm)
                  c1%Atom(n)%x   = xx
                  c1%Active(n)   = .true.
                  c1%Atom(n)%Mult= 1.0_cp
               end do loop1
            end do
            
         type is (atm_std_list_type)
            call allocate_atom_list(npeq*A%natoms,c2)
            do k=1,A%natoms
               if (.not. A%active(k)) cycle
               l=1
               n=n+1
               c2%Atom(n)=A%Atom(k)
               xo= modulo_lat(A%atom(k)%x)
               u(:,l)= xo
               c2%Atom(n)%x=xo

               loop2:do j=2,npeq
                  xx=Apply_OP(SpG%Op(j),xo)
                  xx=modulo_lat(xx)
                  do nt=1,l
                     if (equal_vector(u(:,nt),xx,3)) then
                        c2%atom(n-(l-nt))%occ = c2%atom(n-(l-nt))%occ + A%atom(k)%occ
                        cycle loop2
                     end if
                  end do
                  l=l+1
                  u(:,l)=xx(:)
                  n=n+1
            
                  select case (l)
                     case(:9)
                        write(unit=fmm,fmt="(i1)") l
                     case(10:99)
                        write(unit=fmm,fmt="(i2)") l
                     case(100:999)
                        write(unit=fmm,fmt="(i3)") l
                  end select
                  c2%Atom(n)     = A%Atom(k)
                  c2%Atom(n)%lab = trim(A%Atom(k)%lab)//"_"//adjustl(fmm)
                  c2%Atom(n)%x   = xx
                  c2%Active(n)   = .true.
                  c2%Atom(n)%Mult= 1.0_cp
               end do loop2
            end do
            
         type is (Matm_std_list_type)
            call allocate_atom_list(npeq*A%natoms,c3)
            do k=1,A%natoms
               if (.not. A%active(k)) cycle
               l=1
               n=n+1
               c3%Atom(n)=A%Atom(k)
               xo= modulo_lat(A%atom(k)%x)
               u(:,l)= xo
               c3%Atom(n)%x=xo

               loop3:do j=2,npeq
                  xx=Apply_OP(SpG%Op(j),xo)
                  xx=modulo_lat(xx)
                  do nt=1,l
                     if (equal_vector(u(:,nt),xx,3)) then
                        c3%atom(n-(l-nt))%occ = c3%atom(n-(l-nt))%occ + A%atom(k)%occ
                        cycle loop3
                     end if
                  end do
                  l=l+1
                  u(:,l)=xx(:)
                  n=n+1
            
                  select case (l)
                     case(:9)
                        write(unit=fmm,fmt="(i1)") l
                     case(10:99)
                        write(unit=fmm,fmt="(i2)") l
                     case(100:999)
                        write(unit=fmm,fmt="(i3)") l
                  end select
                  c3%Atom(n)     = A%Atom(k)
                  c3%Atom(n)%lab = trim(A%Atom(k)%lab)//"_"//adjustl(fmm)
                  c3%Atom(n)%x   = xx
                  c3%Active(n)   = .true.
                  c3%Atom(n)%Mult= 1.0_cp
               end do loop3
            end do
            
         type is (atm_ref_list_type)
            call allocate_atom_list(npeq*A%natoms,c4) 
            do k=1,A%natoms
               if (.not. A%active(k)) cycle
               l=1
               n=n+1
               c4%Atom(n)=A%Atom(k)
               xo= modulo_lat(A%atom(k)%x)
               u(:,l)= xo
               c4%Atom(n)%x=xo

               loop4:do j=2,npeq
                  xx=Apply_OP(SpG%Op(j),xo)
                  xx=modulo_lat(xx)
                  do nt=1,l
                     if (equal_vector(u(:,nt),xx,3)) then
                        c4%atom(n-(l-nt))%occ = c4%atom(n-(l-nt))%occ + A%atom(k)%occ
                        cycle loop4
                     end if
                  end do
                  l=l+1
                  u(:,l)=xx(:)
                  n=n+1
            
                  select case (l)
                     case(:9)
                        write(unit=fmm,fmt="(i1)") l
                     case(10:99)
                        write(unit=fmm,fmt="(i2)") l
                     case(100:999)
                        write(unit=fmm,fmt="(i3)") l
                  end select
                  c4%Atom(n)     = A%Atom(k)
                  c4%Atom(n)%lab = trim(A%Atom(k)%lab)//"_"//adjustl(fmm)
                  c4%Atom(n)%x   = xx
                  c4%Active(n)   = .true.
                  c4%Atom(n)%Mult= 1.0_cp
               end do loop4
            end do        
      end select
      
      call allocate_atom_list(n,B)
      select type (B)
         type is (atm_list_type)
            B%active(1:n)=c1%active(1:n)
            B%atom(1:n)=C1%atom(1:n)
            call allocate_atom_list(0,c1)
            
         type is (atm_std_list_type)
            B%active(1:n)=c2%active(1:n)
            B%atom(1:n)=C2%atom(1:n)
            call allocate_atom_list(0,c2)
            
         type is (Matm_std_list_type)
            B%active(1:n)=c3%active(1:n)
            B%atom(1:n)=C3%atom(1:n)
            call allocate_atom_list(0,c3)
            
         type is (atm_ref_list_type)
            B%active(1:n)=c4%active(1:n)
            B%atom(1:n)=C4%atom(1:n) 
            call allocate_atom_list(0,c4)       
      end select
      
   End Subroutine Extend_List
End SubModule   