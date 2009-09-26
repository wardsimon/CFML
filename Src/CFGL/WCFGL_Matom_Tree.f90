module WCFGL_matom_tree
!------------------------------------------------
! Written by Laurent C.Chapon
! September 2004.
! Updated : Now linked list of magnetic atoms
!------------------------------------------------
  use OPENGL
  use WCFGL_constant,                  only : PI
  use WCFGL_atomic_table
  use WCFGL_objects_definition,        only : moment_definition, init_moment
  use WCFGL_metrix
  use CFML_Propagation_Vectors,        only : k_equiv_minus_k
  use CFML_Math_General,               only : modulo_lat
  use CFML_Crystallographic_Symmetry,  only : read_xsym, read_msymm, &
                                              set_spacegroup, space_group_type
  implicit none

  character(len=1),  public, pointer                 :: current_magnetic_lattice => null()
  real             , dimension(:,:), public, pointer :: current_k_list => null()
  character(len=80), dimension(:),   public, pointer :: current_xsym_list => null()
  character(len=80), dimension(:,:), public, pointer :: current_msym_list => null()


!-------------------------------------------------------------------------------
  type :: gl_moment
    integer                             :: ktype    !
    integer                             :: msymtype !
    complex   ,dimension(3)             :: Skj      ! Fourier coeficients
    real                                :: phikj    ! Phas in units of 2pi
    complex   ,dimension(:,:), pointer  :: Skj_eq   ! Skj equivalents
    real      ,dimension(4)             :: color    ! Color
  end type gl_moment
!------------------------------------------------------------------------------
  type :: matom_list
    character(len=10)                          :: Label    ! Label
    character(len=2)                           :: Symbol   ! Chemical Symbol
    real             , dimension(3)            :: xf       ! Atomic position
    real             , dimension(:,:), pointer :: xf_eq    ! List of equivalent positions
    type(gl_moment)  , dimension(:)  , pointer :: moment   ! A gl_moment in the list
    logical                                    :: dead     ! If true then is not displayed
    logical                                    :: group    !
    logical                                    :: envelop  ! If true then the envelop of the modulation on the site is shown
    real             , dimension(4)            :: envelop_color ! Envelop Color
    logical                                    :: edges    ! If edges are visible on the envelop
    real             , dimension(4)            :: edges_color  ! Color of the envelop edges
    real                                       :: edges_radius ! Radius of the envelop edges
    real                                       :: scal     ! Scale factor
    integer                                    :: dlist    ! Store the integer refering to the display list
    type(matom_list), pointer                  :: next
  end type matom_list
!------------------------------------------------------------------------------

  type(matom_list) , pointer, private       :: matom_list_head => null(), &
                                                current_node =>   null(),   &
                                                previous_node =>  null(),  &
                                                temp => null()

  integer, pointer, save, private            :: main_matom_list => null()

  real,             parameter, private       :: eps=0.0001
  contains
!------------------------------------------------------------------------------
  subroutine define_magnetic_phase(lattice,k_list,xsym_list,msym_list,nvk,nxsym,nmsym)
    character(len=*),                 intent(in) :: lattice
    real            , dimension(:,:), intent(in) :: k_list
    character(len=*), dimension(:)  , intent(in) :: xsym_list
    character(len=*), dimension(:,:), intent(in) :: msym_list
    integer                         , intent(in) :: nvk,nxsym,nmsym

  if (associated(current_magnetic_lattice)) deallocate(current_magnetic_lattice)
  if (.not.(associated(current_magnetic_lattice))) allocate(current_magnetic_lattice)
  current_magnetic_lattice=lattice

  if (associated(current_k_list)) deallocate(current_k_list)
  if (.not.(associated(current_k_list))) allocate(current_k_list(3,nvk))
  current_k_list=k_list(:,1:nvk)

  if (associated(current_xsym_list)) deallocate(current_xsym_list)
  if (.not.(associated(current_xsym_list))) allocate(current_xsym_list(nxsym))
  current_xsym_list=xsym_list(1:nxsym)

  if (associated(current_msym_list)) deallocate(current_msym_list)
  if (.not.(associated(current_msym_list))) allocate(current_msym_list(nxsym,nmsym))
  current_msym_list=msym_list(1:nxsym,1:nmsym)

  return

  end subroutine define_magnetic_phase
!------------------------------------------------------------------------------
  subroutine erase_magnetic_phase()

  if (associated(current_magnetic_lattice)) deallocate(current_magnetic_lattice)
  if (associated(current_k_list))           deallocate(current_k_list)
  if (associated(current_xsym_list))        deallocate(current_xsym_list)
  if (associated(current_msym_list))        deallocate(current_msym_list)

  return

  end subroutine erase_magnetic_phase
!------------------------------------------------------------------------------
  subroutine append_k_vector(k)
  real, dimension(3), intent(in) :: k
  real, dimension(:,:), allocatable :: temp_k_list
  integer :: num_k

  if (associated(current_k_list)) then
    num_k=size(current_k_list,DIM=2)
    allocate(temp_k_list(3,num_k+1))
    temp_k_list(:,1:num_k)=current_k_list
    temp_k_list(:,num_k+1)=k
    deallocate(current_k_list)
    allocate(current_k_list(3,num_k+1))
    current_k_list=temp_k_list
    deallocate(temp_k_list)
  else
    allocate(current_k_list(3,1))
    current_k_list(:,1)=k
  end if

  return
  end subroutine append_k_vector
!------------------------------------------------------------------------------
  subroutine erase_k_vector(k_id)
  integer, intent(in) :: k_id
  real, dimension(:,:), allocatable :: temp_k_list
  integer :: num_k

  if (associated(current_k_list)) then
    num_k=size(current_k_list,DIM=2)
    allocate(temp_k_list(3,num_k-1))
    if (k_id==1) then
      temp_k_list=current_k_list(:,2:num_k)
    else
      temp_k_list(:,1:k_id-1)=current_k_list(:,1:k_id-1)
      temp_k_list(:,k_id:num_k-1)=current_k_list(:,k_id+1:num_k)
    end if
      deallocate(current_k_list)
      allocate(current_k_list(3,num_k-1))
      current_k_list=temp_k_list
      deallocate(temp_k_list)
  end if

  return
  end subroutine erase_k_vector
!------------------------------------------------------------------------------
  subroutine modify_k_vector(k_id,new_k)
    integer,            intent(in) :: k_id
    real, dimension(3), intent(in) :: new_k

    if(.not.(associated(current_k_list))) return
    if (k_id<1.or.k_id>size(current_k_list,DIM=2)) then
      return
    else
      current_k_list(:,k_id)=new_k
    end if

    return
  end subroutine modify_k_vector
!------------------------------------------------------------------------------
  subroutine push_matom(label, symbol, xf, point2k, point2msym, Skj, Phikj, &
                         color,scal,dead,group,num_contrib,envelop,envelop_color,edges,edges_color,edges_radius)
    character(len=*),                 intent(in) :: label
    character(len=*),                 intent(in) :: symbol
    real ,   dimension(3),            intent(in) :: xf
    integer, dimension(:),            intent(in) :: point2k, point2msym
    complex, dimension(:,:),          intent(in) :: Skj
    real   , dimension(:)   ,         intent(in) :: phikj
    real , dimension(:,:) ,           intent(in) :: color
    real,                             intent(in) :: scal
    logical,                          intent(in) :: dead
    logical    ,                      intent(in) :: group
    integer    ,                      intent(in) :: num_contrib
    logical    ,                      intent(in), optional :: envelop
    real , dimension(4) ,             intent(in), optional :: envelop_color
    logical    ,                      intent(in), optional :: edges
    real , dimension(4) ,             intent(in), optional :: edges_color
    real                  ,           intent(in), optional :: edges_radius
    ! Local variables
    type(space_group_type) :: space
    integer                :: numlat, nsym, i, j, l, eq_at
    integer :: sim(3,3)
    real    :: tt(3), p_mag, tphas(3), rtrans(3), va(3),vb(3), xm(3)
    complex :: rephas1   !, rephas2
    real    :: skjrealpart(3),skjimpart(3)
    logical :: envelop_flag, edges_flag

    envelop_flag=.false.
    edges_flag=.false.
    ! First check input--------

    ! End of check input-------

    if (.not.(associated(matom_list_head))) then
      allocate(matom_list_head)
      matom_list_head%next =>null()
      current_node => matom_list_head
      if (.not.(associated(moment_definition))) call init_moment()
    else
      allocate(current_node%next)
      current_node%next%next => null()
      current_node => current_node%next
    end if

    current_node%label   = label
    current_node%symbol  = symbol
    xm                   = modulo_lat(xf)
    current_node%xf      = xm
    current_node%group   = group
    current_node%dead    = dead

    current_node%scal    = scal


    call set_spacegroup(current_magnetic_lattice,space)
    numlat=space%numlat

    nsym=size(current_xsym_list)

    ! Allocations
    allocate(current_node%xf_eq(3,nsym*numlat))
    allocate(current_node%moment(num_contrib))

    rtrans=current_node%xf-xf ! Returning translation due to modulo_lat on position

    do l=1, num_contrib
      allocate(current_node%moment(l)%skj_eq(3,nsym*numlat))
      current_node%moment(l)%skj=skj(:,l)
      current_node%moment(l)%phikj=phikj(l)+dot_product(current_k_list(:,point2k(l)),rtrans)
      current_node%moment(l)%color=color(:,l)
      current_node%moment(l)%ktype=point2k(l)
      current_node%moment(l)%msymtype=point2msym(l)
    end do
    ! End of allocations

    do i=1, numlat     ! Sum over number of lattice translations
      do j=1, nsym     ! Sum over number of symmetry operators
        call read_xsym(current_xsym_list(j),1,sim,tt)
        va=matmul(sim,xm)+tt                          !va=matmul(sim,xf)+tt  <---- Usign xf is a bug
        vb=modulo_lat(va+space%latt_trans(:,i))
        current_node%xf_eq(:,j+(i-1)*nsym)=vb
        tphas=vb-va
        do l=1, num_contrib   ! Sum over all contribs
          skjrealpart= real(current_node%moment(l)%skj)
          skjimpart  =AIMAG(current_node%moment(l)%skj)
          rephas1=exp(-2.0*PI*CMPLX(0.0,1.0)*dot_product(current_k_list(:,point2k(l)),tphas))
          call read_msymm(current_msym_list(j,point2msym(l)),sim,p_mag)
          current_node%moment(l)%Skj_eq(:,j+(i-1)*nsym)=CMPLX(matmul(sim,Skjrealpart),&
          matmul(sim,Skjimpart)) &
          *exp(-2.0*PI*CMPLX(0.0,p_mag))*rephas1
        end do
      end do
    end do

    eq_at=size(current_node%xf_eq, DIM=2)

    if (present(envelop)) envelop_flag=envelop
    if (present(edges))  edges_flag=edges

    if (.not.(envelop_flag)) then
        current_node%envelop=.false.
        current_node%envelop_color=(/0.4,0.4,0.,0.5/)
        current_node%edges=.false.
        current_node%edges_color=(/0.,0.4,0.,1./)
        current_node%edges_radius=1.0
                current_node%dlist=glgenLists(1) !  Only one list
    else
                current_node%envelop=envelop_flag
                if (present(envelop_color)) then
                current_node%envelop_color=envelop_color
                else
                current_node%envelop_color(1:3)=color(1:3,1)
                current_node%envelop_color(4)=0.5
                end if
        current_node%edges=edges_flag
        if (present(edges_color)) then
                current_node%edges_color=edges_color
                else
                current_node%edges_color=(/0.0,0.0,0.0,1.0/)
                end if
                if (present(edges_radius)) then
        current_node%edges_radius=edges_radius
        else
        current_node%edges_radius=2.0
        end if
        if (.not.(edges_flag)) then
                current_node%dlist=glgenLists(1+eq_at) !   One list for moment and one list /envelop
        else
                current_node%dlist=glgenLists(1+2*eq_at) !   One list for moment and one list /envelop and 1 list /edges
        end if
    end if

    call construct_matom_dlist(current_node)

    call create_main_matom_list()

    return

  end subroutine push_matom
!------------------------------------------------------------------------------
  subroutine pop_matom(matom_id)
    integer, intent(in) :: matom_id
    ! Local variables---------------------
    integer :: id_count, l
    integer :: nlist_to_destroy, eq_at
    temp => matom_list_head                      ! Back to top
    id_count=0
    do
      if (.not. associated(temp)) exit          ! Tail of the list
      id_count=id_count+1
      if (id_count==matom_id) exit                ! We have a match
      previous_node => temp
      temp => temp%next
    end do

    if (associated(temp)) then
      if (associated(temp,matom_list_head)) then   ! If id is the head of the list
         matom_list_head => temp%next
       else                                        ! If it is not
        previous_node%next => temp%next
      endif
      if (associated(temp,current_node)) current_node => previous_node
      nlist_to_destroy=1
      eq_at=size(temp%xf_eq, DIM=2)
      if (temp%envelop) then
        if (.not.(temp%edges)) then
          nlist_to_destroy=1+eq_at
        else
          nlist_to_destroy=1+2*eq_at
        end if
      end if
        if (associated(temp%xf_eq))  deallocate(temp%xf_eq)
        if (associated(temp%moment)) then
          do l=1, size(temp%moment)
            if (associated(temp%moment(l)%skj_eq)) deallocate(temp%moment(l)%skj_eq)
          end do
          deallocate(temp%moment)
        end if
        call gldeletelists(temp%dlist,nlist_to_destroy)          ! Delete the display list to free memory
        deallocate(temp)                                         ! Deallocate the pointer
    end if

    return

  end subroutine pop_matom
!------------------------------------------------------------------------------
  subroutine construct_matom_dlist(matom_pointer)
    type(matom_list), intent(in) :: matom_pointer
    ! Local variables---------------------
    integer       :: i, j, k, l, m, eq_at, ei, envLoop
    real          :: kdotRl
    complex       :: expterm
    real(glfloat) :: pos(3), trans(3), Rl(3), mom(3), spher(3)
    complex, dimension(3) :: skj, skjstar

    if (.not.(associated(current_box))) call init_box()

    eq_at=size(matom_pointer%xf_eq, DIM=2)

    if (matom_pointer%envelop) then ! IF envelop keyword is true
        do envLoop=1,eq_at
                if (.not.(matom_pointer%edges)) then
                        call glnewlist(matom_pointer%dlist+envLoop,GL_COMPILE) ! Begin envelop
                else
                        call glnewlist(matom_pointer%dlist+envLoop*2-1,GL_COMPILE) ! Begin envelop
                end if
                call glmaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,matom_pointer%envelop_color)
                call glBegin(GL_TRIANGLE_FAN)
                        call glVertex3f(0.0,0.0,0.0)
                do ei=0, 200 ! 200 points for the envelop
                        mom=0.0
                        do m=1, size(matom_pointer%moment)
                                kdotRl=-2.0*PI*ei/200.0
                                expterm=exp(CMPLX(0.0,kdotRl))
                                skj=matom_pointer%moment(m)%skj_eq(:,envLoop)
                                skjstar=conjg(skj)
                                if (k_equiv_minus_k(current_k_list(:,matom_pointer%moment(m)%ktype),current_magnetic_lattice)) then
                                        mom=mom+real(skj*expterm*exp(CMPLX(0.0,-2.0*PI*matom_pointer%moment(m)%phikj)))
                                else
                                        mom=mom+0.5*real((skj*expterm*exp(CMPLX(0.0,-2.0*PI*matom_pointer%moment(m)%phikj))) &
                                        +(skjstar*conjg(expterm)*exp(CMPLX(0.0,2.0*PI*matom_pointer%moment(m)%phikj))))
                                end if
                        end do
                        mom=matom_pointer%scal*mom*0.5 ! if scale is applied
                        if (associated(current_cell)) mom=matmul(current_cell%t_mtx_unit,mom)
                        call glVertex3fv(mom)
                end do
                call glEnd()
                call glEndList()  !Finish envelop
        end do

        if (matom_pointer%edges) then
                do envLoop=1,eq_at
                                call glnewlist(matom_pointer%dlist+envLoop*2,GL_COMPILE) ! Begin envelop
                                !call glpushattrib(GL_LIGHTING)
                        !call gldisable(GL_LIGHTING)
                                call glLineWidth(matom_pointer%edges_radius)
                        call glBegin(GL_LINE_STRIP)
                        call glColor4fv(matom_pointer%edges_color)
                        do ei=0, 200 ! 200 points for the envelop
                                mom=0.0
                                do m=1, size(matom_pointer%moment)
                                        kdotRl=-2.0*PI*ei/200.0
                                        expterm=exp(CMPLX(0.0,kdotRl))
                                        skj=matom_pointer%moment(m)%skj_eq(:,envLoop)
                                        skjstar=conjg(skj)
                                        if (k_equiv_minus_k(current_k_list(:,matom_pointer%moment(m)%ktype),current_magnetic_lattice)) then
                                                mom=mom+real(skj*expterm*exp(CMPLX(0.0,-2.0*PI*matom_pointer%moment(m)%phikj)))
                                        else
                                                mom=mom+0.5*real((skj*expterm*exp(CMPLX(0.0,-2.0*PI*matom_pointer%moment(m)%phikj))) &
                                        +(skjstar*conjg(expterm)*exp(CMPLX(0.0,2.0*PI*matom_pointer%moment(m)%phikj))))
                                        end if
                                end do
                                mom=matom_pointer%scal*mom*0.5 ! if scale is applied
                                if (associated(current_cell)) mom=matmul(current_cell%t_mtx_unit,mom)
                                call glVertex3fv(mom)
                        end do
                        call glEnd()
                        call glEndList()  !Finish envelop
                        !call glpopattrib()
                end do
        end if
     end if


     call glnewlist(matom_pointer%dlist,GL_COMPILE)

     do l=1, eq_at
      do i=int(current_box(1))-1,int(current_box(2))+1
        do j=int(current_box(3))-1,int(current_box(4))+1
          do k=int(current_box(5))-1, int(current_box(6))+1
            Rl=real((/i,j,k/))
            pos=matom_pointer%xf_eq(:,l)+Rl
            if (in_limit(pos)) then
              trans=f2c(pos)
              call glpushmatrix()
              call gltranslatef(trans(1),trans(2),trans(3))
              if (.not.(matom_pointer%group)) then
                do m=1, size(matom_pointer%moment)
                  kdotRl=-2.0*PI*dot_product(current_k_list(:,matom_pointer%moment(m)%ktype),Rl)
                  expterm=exp(CMPLX(0.0,kdotRl))
                  skj=matom_pointer%moment(m)%skj_eq(:,l)
                  skjstar=conjg(skj)
                  if (k_equiv_minus_k(current_k_list(:,matom_pointer%moment(m)%ktype),current_magnetic_lattice)) then
                    mom=real(skj*expterm*exp(CMPLX(0.0,-2.0*PI*matom_pointer%moment(m)%phikj)))
                  else
                    mom=0.5*real((skj*expterm*exp(CMPLX(0.0,-2.0*PI*matom_pointer%moment(m)%phikj)))&
                    +(skjstar*conjg(expterm)*exp(CMPLX(0.0,2.0*PI*matom_pointer%moment(m)%phikj))))
                  end if
                  if (associated(current_cell)) mom=matmul(current_cell%t_mtx_unit,mom)
                  spher=cart2spher(mom)
                  call glpushmatrix()
                  call glrotatef(spher(2),0.0,0.0,1.0)
                  call glrotatef(-spher(3),0.0,1.0,0.0)
                  if (spher(1)>eps) then
                    call glscalef(spher(1)*matom_pointer%scal,1.0,1.0)
                    call glmaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,matom_pointer%moment(m)%color)
                    call glcalllist(moment_definition)
                  end if
                  call glpopmatrix()
                end do
              else
                mom=0.0
                do m=1, size(matom_pointer%moment)
                  kdotRl=-2.0*PI*dot_product(current_k_list(:,matom_pointer%moment(m)%ktype),Rl)
                  expterm=exp(CMPLX(0.0,kdotRl))
                  skj=matom_pointer%moment(m)%skj_eq(:,l)
                  skjstar=conjg(skj)
                  if (k_equiv_minus_k(current_k_list(:,matom_pointer%moment(m)%ktype),current_magnetic_lattice)) then
                    mom=mom+real(skj*expterm*exp(CMPLX(0.0,-2.0*PI*matom_pointer%moment(m)%phikj)))
                  else
                    mom=mom+0.5*real((skj*expterm*exp(CMPLX(0.0,-2.0*PI*matom_pointer%moment(m)%phikj))) &
                    +(skjstar*conjg(expterm)*exp(CMPLX(0.0,2.0*PI*matom_pointer%moment(m)%phikj))))
                  end if
                end do
                if (associated(current_cell)) mom=matmul(current_cell%t_mtx_unit,mom)
                spher=cart2spher(mom)
                call glpushmatrix()
                call glrotatef(spher(2),0.0,0.0,1.0)
                call glrotatef(-spher(3),0.0,1.0,0.0)
                if (spher(1)>eps) then
                  call glscalef(spher(1)*matom_pointer%scal,1.0,1.0)
                  call glmaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,matom_pointer%moment(1)%color)
                  call glcalllist(moment_definition)
                end if
                call glpopmatrix()
                if (matom_pointer%envelop) then
                        call glmaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,matom_pointer%envelop_color)
                        if (matom_pointer%edges) then
                                call glcalllist(matom_pointer%dlist+2*l-1)
                                call glpushattrib(GL_LIGHTING)
                                call gldisable(GL_LIGHTING)
                                call glColor4fv(matom_pointer%edges_color)
                                call glcalllist(matom_pointer%dlist+2*l)
                                call glPopAttrib()
                        else
                                call glcalllist(matom_pointer%dlist+l)
                        end if
                end if
              end if
              call glpopmatrix()
            end if
          end do
        end do
      end do
    end do

    call glendlist()

    return

  end subroutine construct_matom_dlist
!------------------------------------------------------------------------------
subroutine construct_matom_envelop(matom_pointer)
type(matom_list), intent(in) :: matom_pointer

if (.not.(associated(current_box))) call init_box()


end subroutine construct_matom_envelop
!------------------------------------------------------------------------------
  subroutine change_matoms_plotting_status(status)
    logical, intent(in) :: status

    if (.not.(associated(matom_list_head))) return

    temp => matom_list_head
    do
      temp%dead=status
      if (.not.(associated(temp%next))) exit
      temp => temp%next
    end do

    call create_main_matom_list()

    return

  end subroutine change_matoms_plotting_status
!----------------------------------------------------------------------------------
  subroutine create_main_matom_list()

    if (.not.(associated(matom_list_head))) return

   if (.not.(associated(main_matom_list))) then
     allocate(main_matom_list)
     main_matom_list=glgenlists(1)
   end if

    call glnewlist(main_matom_list,GL_COMPILE)

    temp => matom_list_head
    do
      if (.not.(temp%dead)) then
        call glcalllist(temp%dlist)
      end if
      if (.not.(associated(temp%next))) exit
      temp => temp%next
    end do
    call glendlist()

  return

  end subroutine create_main_matom_list
!----------------------------------------------------------------------------------
subroutine update_all_matoms()

   if (.not.(associated(matom_list_head))) return

     temp => matom_list_head
     do
       call construct_matom_dlist(temp)
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do

   return

end subroutine update_all_matoms
!----------------------------------------------------------------------------------

  subroutine draw_matoms()

    if (.not.(associated(main_matom_list))) return
    if (new_box .or. new_cell) call update_all_matoms()
      call glcalllist(main_matom_list)
    return

  end subroutine draw_matoms
!----------------------------------------------------------------------------------
subroutine empty_matom_list()

  do while(associated(matom_list_head))
    call pop_matom(1)
  end do

  if (associated(main_matom_list)) deallocate(main_matom_list)

  return

end subroutine empty_matom_list
!----------------------------------------------------------------------------------
  subroutine write_magnetic_structure(unit_no)
    integer, intent(in) :: unit_no
    character(len=10) :: dead_string
    character(len=6)  :: group_string
    integer :: i, j


  write(unit=unit_no,fmt='(a)') "{"

  if (associated(current_magnetic_lattice)) then
    write(unit=unit_no,fmt='(2a)') "LATTICE ", current_magnetic_lattice
  end if

  if (associated(current_k_list)) then
    do i=1, size(current_k_list,DIM=2)
    write(unit=unit_no,fmt='(a,3f8.4)') "K ", current_k_list(:,i)
    end do
  end if

  if (associated(current_xsym_list).and. associated(current_msym_list)) then
    do i=1, size(current_xsym_list)
      write(unit=unit_no,fmt='(2a)') "SYMM ", current_xsym_list(i)
      do j=1, size(current_msym_list,DIM=2)
        write(unit=unit_no,fmt='(a)') current_msym_list(i,j)
      end do
    end do
  end if

  if (.not.(associated(matom_list_head))) then
    write(unit=unit_no,fmt='(a)') "}"
    return
  end if

  temp => matom_list_head
     do

       if (temp%dead) then
         dead_string=" NODISPLAY"
       else
         dead_string=" "
       end if

       if (temp%group) then
         group_string=" GROUP"
       else
         group_string=" "
       end if
       if (.not.(temp%envelop)) then
       write(unit=unit_no,fmt='(3a,3F9.5,a,F5.2,2a)') "MATOM ", temp%label, &
       temp%symbol, temp%xf, " SCALE ", temp%scal, group_string, dead_string
       else
        if (.not.(temp%edges)) then
        write(unit=unit_no,fmt='(3a,3F9.5,a,F5.2,2a,4F5.2,a)') "MATOM ", temp%label, &
        temp%symbol, temp%xf, " SCALE ", temp%scal, group_string, " ENVELOP ENVELOPCOL ",temp%envelop_color, dead_string
        else
        write(unit=unit_no,fmt='(3a,3F9.5,a,F5.2,2a,4F5.2,a,4F5.2,a,F5.2,a)') "MATOM ", temp%label, &
        temp%symbol, temp%xf, " SCALE ", temp%scal, group_string, " ENVELOP ENVELOPCOL ",temp%envelop_color, " EDGE EDGECOL ",temp%edges_color, " RADIUS ", temp%edges_radius, dead_string
        end if
       end if
       do i=1, size(temp%moment)
         write(unit=unit_no,fmt='(a,2i3,7f9.3,a,4f6.3)') "SKP ", temp%moment(i)%ktype, &
         temp%moment(i)%msymtype,&
         real(temp%moment(i)%skj(1)),&
         real(temp%moment(i)%skj(2)),&
         real(temp%moment(i)%skj(3)),&
         AIMAG(temp%moment(i)%skj(1)),&
         AIMAG(temp%moment(i)%skj(2)),&
         AIMAG(temp%moment(i)%skj(3)),&
         temp%moment(i)%phikj, " COLOR ", temp%moment(i)%color
       end do

       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do

  write(unit=unit_no,fmt='(a)') "}"

  return

  end subroutine write_magnetic_structure
!----------------------------------------------------------------------------------
  subroutine write_magnetic_moment_list(unit_no)
    integer, intent(in)   :: unit_no
    integer               :: i, j, k, l, m, eq_at
    real                  :: kdotRl, modulus
    complex               :: expterm
    real(glfloat)         :: pos(3), Rl(3), mom(3), tot_mom(3)
    complex, dimension(3) :: skj, skjstar
    character(len=80)     :: Rlstring

    if (.not.(associated(matom_list_head))) then
      write(unit=unit_no,fmt='(a)') "No moments found"
      return
    end if

    if (associated(current_magnetic_lattice)) &
    write(unit=unit_no,fmt='(2a)') " Magnetic lattice type : ", current_magnetic_lattice

    write(unit=unit_no,fmt='(a)') " "

    if (associated(current_k_list)) then
        write(unit=unit_no,fmt='(a)') " Magnetic k-vectors :"
      do i=1, size(current_k_list,DIM=2)
        write(unit=unit_no,fmt='(3f8.3)') current_k_list(:,i)
      end do
    end if

    write(unit=unit_no,fmt='(a)') " "

    if (associated(current_xsym_list).and. associated(current_msym_list)) then
      write(unit=unit_no,fmt='(a)') " Symmetry operations :"
      do i=1, size(current_xsym_list)
        write(unit=unit_no,fmt='(2a)') "SYMM ", current_xsym_list(i)
        do j=1, size(current_msym_list,DIM=2)
          write(unit=unit_no,fmt='(a)') current_msym_list(i,j)
        end do
      end do
    end if

    write(unit=unit_no,fmt='(a)') " "

    temp => matom_list_head

    do
      write(unit=unit_no,fmt='(3a)') "Atom : ", temp%label, temp%symbol
      write(unit=unit_no,fmt='(a)') "------------------------------------------------------------------------------------"
      write(unit=unit_no,fmt='(a)') "        x        y        z        Translation    k MSYM     m(a)     m(b)     m(c)     Mtot"
      eq_at=size(temp%xf_eq, DIM=2)
      do l=1, eq_at
        write(unit=unit_no,fmt='(3f9.3)') temp%xf_eq(:,l)
        do i=int(current_box(1))-1,int(current_box(2))+1
          do j=int(current_box(3))-1,int(current_box(4))+1
            do k=int(current_box(5))-1, int(current_box(6))+1
              Rl=real((/i,j,k/))
              pos=temp%xf_eq(:,l)+Rl
              if (in_limit(pos)) then
                write(unit=Rlstring,fmt='(a,i4,a,i4,a,i4,a)') "(",i,",",j,",",k,")"
                write(unit=unit_no,fmt='(2a)') "                              ", Rlstring
                tot_mom=0.0
                do m=1, size(temp%moment)
                  kdotRl=-2.0*PI*dot_product(current_k_list(:,temp%moment(m)%ktype),Rl)
                  expterm=exp(CMPLX(0.0,kdotRl))
                  skj=temp%moment(m)%skj_eq(:,l)
                  skjstar=conjg(skj)
                  if (k_equiv_minus_k(current_k_list(:,temp%moment(m)%ktype),current_magnetic_lattice)) then
                    mom=real(skj*expterm*exp(CMPLX(0.0,-2.0*PI*temp%moment(m)%phikj)))
                    tot_mom=tot_mom+mom
                  else
                    mom=0.5*real((skj*expterm*exp(CMPLX(0.0,-2.0*PI*temp%moment(m)%phikj)))&
                    +(skjstar*conjg(expterm)*exp(CMPLX(0.0,2.0*PI*temp%moment(m)%phikj))))
                    tot_mom=tot_mom+mom
                  end if
                  write(unit=unit_no,fmt='(a,2i5,3f9.3)') "                                              ", &
                  temp%moment(m)%ktype,temp%moment(m)%msymtype, mom(1), mom(2), mom(3)
                end do
                  write(unit=unit_no,fmt='(a)')     "                                                        ------------------------------------"
                  mom=matmul(current_cell%t_mtx_unit,tot_mom) !Transform to cartesian frame
                  modulus = sqrt( dot_product( mom, mom ) )
                  write(unit=unit_no,fmt='(a,4f9.3)')"                                                        ",tot_mom, modulus
              end if
            end do
          end do
        end do
      end do

      if (.not.(associated(temp%next))) exit
      temp => temp%next
    end do

     return

  end subroutine write_magnetic_moment_list
!----------------------------------------------------------------------------------
end module WCFGL_matom_tree
