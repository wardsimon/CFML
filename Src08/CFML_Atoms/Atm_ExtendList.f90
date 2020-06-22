!!----
!!----
!!----
!!
SubModule (CFML_Atoms) Generating_Atoms_inCell
   Contains

   !!----
   !!---- EXTEND_LIST
   !!----
   !!----    Subroutine to generate atoms in the primitive (conven=.false.) or the conventional
   !!----    unit cell (conven=.true.), Excluding atoms with A%atom(:)%active=.false.
   !!----
   !!---- Update: February - 2005
   !!
   Module Subroutine Extend_List(A, B, Spg, Type_Atm,Conven)
      !---- Arguments ----!
      type(atlist_type),   intent(in)     :: A         ! Atom list (asymmetric unit)
      type(atlist_type),   intent(in out) :: B         ! Atom list into the unit cell
      type(SpG_Type),      intent(in)     :: SpG       ! SpaceGroup
      character(len=*),    intent(in)     :: Type_Atm  ! !Atomic type: Atm, Atm_Std, MAtm_Std, Atm_Ref, MAtm_Ref
      logical, optional,   intent(in)     :: Conven    ! If present and .true. using the whole conventional unit cell

      !---- Local Variables ----!
      character(len=4)                      :: fmm
      integer                               :: k,j,l,nt,npeq,n
      real(kind=cp),dimension(3)            :: xo,xx
      real(kind=cp),dimension(3,Spg%multip) :: u
      logical                               :: ccell

      type(atlist_type) :: c_atm

      !> Init
      call clear_error()
      ccell=.false.
      if (present(conven)) ccell=conven

      !> Check arguments
      if (.not. extends_type_of(B,A) .or. (.not. same_type_as(B,A)) ) then
         err_CFML%IErr=1
         err_CFML%Msg="Extend_List@CFML_Atoms: Incompatible arguments!"
         return
      end if

      npeq=SpG%numops
      if (SpG%centred == 2) npeq=npeq*2
      if (ccell) npeq=SpG%multip
      call allocate_atom_list(npeq*A%natoms,c_atm,Type_Atm,3)
      n=0
      do k=1,A%natoms
         l=1       ! Number of symetric atom
         n=n+1     ! Number of atom in the list

         xo= modulo_lat(A%atom(k)%x)
         u(:,l)= xo
         c_atm%atom(n)%x=xo

         loop:do j=2,npeq
            xx=Apply_OP(SpG%Op(j),xo)
            xx=modulo_lat(xx)
            do nt=1,l
               if (equal_vector(u(:,nt),xx,3)) then
                  c_atm%atom(n-(l-nt))%occ = c_atm%atom(n-(l-nt))%occ + A%atom(k)%occ
                  cycle loop
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
            !c_atm%Atom(n)     = A%atom(k)     ! Valid for Intel and not for GFortran
            c_atm%Atom(n)%lab     = trim(A%atom(k)%lab)//"_"//adjustl(fmm)
            c_atm%Atom(n)%Chemsymb=A%atom(k)%Chemsymb
            c_atm%Atom(n)%Z       =A%atom(k)%Z
            c_atm%Atom(n)%Mult    = real(l)/SpG%multip
            c_atm%Atom(n)%Charge  =A%atom(k)%Charge
            c_atm%Atom(n)%x       = xx
            c_atm%Atom(n)%UType   =A%atom(k)%UType
            c_atm%Atom(n)%ThType  =A%atom(k)%ThType
            c_atm%Atom(n)%U_iso   =A%atom(k)%U_iso
            c_atm%Atom(n)%Magnetic=A%atom(k)%Magnetic
            c_atm%Atom(n)%Moment  =A%atom(k)%Moment
            c_atm%Atom(n)%SfacSymb=A%atom(k)%SfacSymb
            c_atm%Atom(n)%Ind_ff  =A%atom(k)%Ind_ff

            !select type(atm => A%atom(k))
            !   class is (Atm_Std_Type)
            !      c_atm%Atom(n)%u_iso_std=atm%u_iso_std
            !      !c_atm%Atom(n)%x_std=A%atom(k)%x_std
            !      !c_atm%Atom(n)%occ_std=A%atom(k)%occ_std
            !end select

            !select type(atm => A%atom(k))
            !   type is (MAtm_Std_Type)
            !      c_atm%Atom(n)%wyck   = atm%wyck
            !      c_atm%Atom(n)%n_mc   = atm%n_mc
            !      c_atm%Atom(n)%n_dc   = atm%n_dc
            !      !c_atm%Atom(n)%mcs    = A%atom(k)%mcs
            !      !c_atm%Atom(n)%mcs_std= A%atom(k)%mcs_std
            !      !c_atm%Atom(n)%dcs    = A%atom(k)%dcs
            !      !c_atm%Atom(n)%dcs_std= A%atom(k)%dcs_std
            !end select

            c_atm%Active(n)   = A%Active(k)

         end do loop
      end do

      if (n == 0) then
         err_CFML%IErr=1
         err_CFML%Msg="Extend_List@CFML_Atoms: Number of extended atoms is zero!"
         call allocate_atom_list(0,c_atm,Type_Atm,3)
         return
      end if

      call allocate_atom_list(n,B,Type_Atm,3)
      B%Active=c_atm%active(1:n)
      B%Atom=c_atm%atom(1:n)

      !> DEallocate
      call allocate_atom_list(0,c_atm,Type_Atm,3)

   End Subroutine Extend_List

   !!----
   !!----    Module Subroutine Set_Atom_Equiv_List(SpG,cell,A,Ate,lun)
   !!----      type(SpG_Type),             intent(in) :: SpG
   !!----      type(Cell_G_Type),          intent(in) :: Cell
   !!----      type(Atlist_Type),          intent(in) :: A
   !!----      type(Atom_Equiv_List_Type), intent(out):: Ate
   !!----      integer, optional,          intent(in) :: lun
   !!----
   !!---- Subroutine constructing the list of all atoms in the unit cell.
   !!---- The atoms are in a structure of type "Atom_Equiv_List_Type" containing
   !!---- just the fractional coordinates of all the atoms in the cell.
   !!---- This a simplified version of the Extend_List Subroutine useful for geometric
   !!---- calculations, using the type Atom_Equiv_List_Type, without the burden of
   !!---- all components of Aton_Type
   !!----
   !!---- Updated: May 2020
   !!
   Module Subroutine Set_Atom_Equiv_List(SpG,cell,A,Ate,lun)
     type(SpG_Type),             intent(in) :: SpG
     type(Cell_G_Type),          intent(in) :: Cell
     type(Atlist_Type),          intent(in) :: A
     type(Atom_Equiv_List_Type), intent(out):: Ate
     integer, optional,          intent(in) :: lun

     ! local variables
     real(kind=cp),  dimension(3)            :: xx,xo,v,xc
     real(kind=cp),  dimension(3,SpG%Multip) :: u
     character(len=20),dimension(SpG%Multip) :: label
     integer                                 :: k,j,L,nt
     character (len=6)                       :: fmm
     character (len=20)                      :: nam
     real(kind=cp), parameter                :: epsi = 0.002

     if (.not. allocated (Ate%atm)) allocate(Ate%atm(A%natoms))
     ate%nauas=A%natoms
     if (present(lun))  then
        write(unit=lun,fmt="(/,a)") "     LIST OF ATOMS INSIDE THE CONVENTIONAL UNIT CELL "
        write(unit=lun,fmt="(a,/)") "     =============================================== "
     end if
     do k=1,A%natoms
        ate%atm(k)%ChemSymb = A%atom(k)%ChemSymb
        xo(:) =Modulo_Lat(A%atom(k)%x(:))
        L=1
        u(:,L)=xo(:)
        xc =matmul(cell%Cr_Orth_cel,xo)
        if (present(lun))then
         write(unit=lun,fmt="(/,a,a)") " => Equivalent positions of atom: ",A%atom(k)%lab
         write(unit=lun,fmt="(a)")  &
         "                                    x         y         z          Xc        Yc        Zc"
        end if
        fmm="(a,i1)"
        write(unit=label(L),fmt=fmm) trim(A%Atom(k)%lab)//"_",L
        nam=label(L)
        if (present(lun)) write(unit=lun,fmt="(3a,3f10.5,a,3f10.5)") "       ",nam,"  ", xo,"  ", xc

        do_eq:DO j=2,SpG%multip
           xx=Apply_OP(SpG%Op(j),xo)
           xx=modulo_lat(xx)
           do nt=1,L
              v=u(:,nt)-xx(:)
              if (sum(abs((v))) < epsi ) cycle do_eq
           end do
           L=L+1
           u(:,L)=xx(:)
           if ( L > 9 .and. L < 100)  fmm="(a,i2)"
           if ( L >= 100 )  fmm="(a,i3)"
           write(unit=label(L),fmt=fmm) trim(A%Atom(k)%lab)//"_",L
           nam=Label(L)
           xc=matmul(cell%Cr_Orth_cel,xx)
           if (present(lun)) write(unit=lun,fmt="(3a,3f10.5,a,3f10.5)") "       ",nam,"  ", xx,"  ", xc
        end do do_eq

        if(allocated(Ate%Atm(k)%Lab)) deallocate(Ate%Atm(k)%Lab)
        allocate(Ate%Atm(k)%lab(L))
        if(allocated(Ate%Atm(k)%x)) deallocate(Ate%Atm(k)%x)
        allocate(Ate%Atm(k)%x(3,L))

        Ate%Atm(k)%mult=L
        do j=1,Ate%Atm(k)%mult
          Ate%Atm(k)%lab(j)=Label(j)
          Ate%Atm(k)%x(:,j)=u(:,j)
        end do
     end do
     if (present(lun))  write(unit=lun,fmt="(/)")
   End Subroutine Set_Atom_Equiv_List

End SubModule Generating_Atoms_inCell