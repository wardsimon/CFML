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
      type(atlist_type),   intent(in)     :: A         ! Atom list (asymmetric unit)
      type(atlist_type),   intent(in out) :: B         ! Atom list into the unit cell
      type(SpG_Type),      intent(in)     :: SpG       ! SpaceGroup
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
      call allocate_atom_list(npeq*A%natoms,c_atm)
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
         err_CFML%Msg="Extend_List@CFML_Atoms: Number of extended atoms was zero!"
         call allocate_atom_list(0,c_atm)
         return
      end if   
      
      call allocate_atom_list(n,B)
      B%Active=c_atm%active(1:n)
      B%Atom=c_atm%atom(1:n)
      
      !> DEallocate
      call allocate_atom_list(0,c_atm)
      
   End Subroutine Extend_List
End SubModule   