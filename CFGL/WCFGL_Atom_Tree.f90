module WCFGL_atom_tree
!------------------------------------------------
! Written by Laurent C.Chapon
! July 2004.
! Updated :
!------------------------------------------------

  use OPENGL
  use WCFGL_geometry
  use WCFGL_objects_definition, only : atom_definition, init_atom
  use WCFGL_atomic_table
  use WCFGL_glatom
  use WCFGL_metrix
  use WCFGL_trackball

  implicit none

!------------------------------------------------------------------------------
  type :: atom_list
    type(gl_atom)            :: atom  ! A gl_Atom in the list
    logical                  :: dead  ! If true then is not display
    integer                  :: dlist ! Store the integer refering to the display list
    type(atom_list), pointer :: next
  end type atom_list
!------------------------------------------------------------------------------
  type(atom_list) , pointer, private       :: atom_list_head => null(), &
                                              current_node =>   null(),   &
                                              previous_node =>  null(),  &
                                              temp => null()

  integer, pointer, save, private          :: main_atom_list => null()


  real(glfloat), parameter, private        :: eps=0.0002

  contains
!------------------------------------------------------------------------------
  subroutine draw_atoms()

    if (new_box .or. new_cell .or. new_spacegroup) call update_all_atoms()

    if (associated(main_atom_list)) then
      call glcalllist(main_atom_list)
    end if

    return

  end subroutine draw_atoms
!-----------------------------------------------------------------------------------
  subroutine push_atom(label, symbol, xf, Biso, radius,dead,color)
    character(len=10),      intent(in)                 :: label
    character(len=2),       intent(in)                 :: symbol
    real         ,          intent(in)                 :: xf(3)
    real         ,          intent(in)                 :: Biso
    real         ,          intent(in)                 :: radius
    logical      ,          intent(in)                 :: dead
    real         ,          intent(in), optional       :: color(4)
    ! Local variables
    integer :: ier
    real, dimension(4) :: my_color

    if (.not.(associated(atom_list_head))) then
      allocate(atom_list_head)
      atom_list_head%next =>null()
      current_node => atom_list_head
      if (.not.(associated(atom_definition))) call init_atom()
    else
      allocate(current_node%next)
      current_node%next%next => null()
      current_node => current_node%next
    end if

    if (.not.(associated(current_space_group))) call init_spacegroup()

    if (.not.(present(color))) then
      call new_gl_atom(label, symbol, xf, Biso, radius, get_color_from_symbol(symbol), current_space_group, current_node%atom)
      current_node%atom%colorfromtable=.true.
    else
      call new_gl_atom(label, symbol, xf, Biso, radius, color, current_space_group, current_node%atom)
      current_node%atom%colorfromtable=.false.
    end if

    current_node%dead=dead
    current_node%dlist=glgenlists(1)                                 ! Associate a display list.

    call construct_atom_dlist(current_node)

    call create_main_atom_list()

    return

  end subroutine push_atom
!----------------------------------------------------------------------------------
  subroutine construct_atom_dlist(atom_pointer)
    type(atom_list), intent(inout) :: atom_pointer
    ! Local variables
    real     :: pos(3), trans(3)
    integer  :: i, j, k, l, eq_at

    if (.not.(associated(current_box))) call init_box()

    if (new_spacegroup) call calculate_eq_positions(atom_pointer%atom,current_space_group)

    eq_at=size(atom_pointer%atom%xf_eq, DIM=2)

    call glnewlist(atom_pointer%dlist,GL_COMPILE)                    ! Create display list
    call glmaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,atom_pointer%atom%color(1:3))
    do i=int(current_box(1))-1,int(current_box(2))+1
      do j=int(current_box(3))-1,int(current_box(4))+1
        do k=int(current_box(5))-1, int(current_box(6))+1
          do l=1, eq_at
            pos=atom_pointer%atom%xf_eq(:,l)+real((/i,j,k/))
            if (in_limit(pos)) then
              trans=f2c(pos)
              call glpushmatrix()
              call gltranslatef(trans(1),trans(2),trans(3))
              call glscalef(atom_pointer%atom%radius,atom_pointer%atom%radius,atom_pointer%atom%radius)
              call glcalllist(atom_definition)
              call glpopmatrix()
            end if
          end do
        end do
      end do
    end do
  call glendlist()

  return

  end subroutine construct_atom_dlist
!----------------------------------------------------------------------------------
  subroutine pop_atom(atom_id)
    integer, intent(in) :: atom_id
    integer :: id_count

    temp => atom_list_head                      ! Back to top
    id_count=0
    do
      if (.not. associated(temp)) exit          ! Tail of the list
      id_count=id_count+1
      if (id_count==atom_id) exit                ! We have a match
      previous_node => temp
      temp => temp%next
    end do

    if (associated(temp)) then
      if (associated(temp,atom_list_head)) then   ! If id is the head of the list
         atom_list_head => temp%next
       else                                        ! If it is not
        previous_node%next => temp%next
      endif
      if (associated(temp,current_node)) current_node => previous_node
        if (associated(temp%atom%xf_eq)) deallocate(temp%atom%xf_eq)
        call gldeletelists(temp%dlist,1)          ! Delete the display list to free memory
        deallocate(temp)                          ! Deallocate the pointer
    endif

    return

  end subroutine pop_atom
!----------------------------------------------------------------------------------
subroutine create_main_atom_list()

   if (.not.(associated(atom_list_head))) return

   if (.not.(associated(main_atom_list))) then
     allocate(main_atom_list)
     main_atom_list=glgenlists(1)
   end if

   call glnewlist(main_atom_list,GL_COMPILE)

     temp => atom_list_head
     do
       if (.not.(temp%dead)) then
         call glcalllist(temp%dlist)
       end if
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do
   call glendlist()

   return

end subroutine create_main_atom_list
!----------------------------------------------------------------------------------
subroutine update_all_atoms()

   if (.not.(associated(atom_list_head))) return

     temp => atom_list_head
     do
       call construct_atom_dlist(temp)
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do

   return

end subroutine update_all_atoms
!----------------------------------------------------------------------------------
subroutine search_atom_list_by_label(label,atomlist,success)
  character(len=10)                        , intent(in)   :: label
  type(gl_atom), dimension (:), allocatable, intent(out)  :: atomlist
  logical                                  , intent(out)  :: success
  type(gl_atom), dimension(400) :: templist
  integer :: count, i


   count=0
   success=.false.

   if (.not.(associated(atom_list_head))) return

     temp => atom_list_head
     do
       if (temp%atom%label==label) then
         success=.true.
         count=count+1
         templist(count)=temp%atom
       end if
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do

   allocate(atomlist(count))
   atomlist(:)=templist(1:count)

   return

end subroutine search_atom_list_by_label
!----------------------------------------------------------------------------------
subroutine change_atoms_plotting_status(status)
  logical, intent(in) :: status

  if (.not.(associated(atom_list_head))) return

  temp => atom_list_head
  do
    temp%dead=status
    if (.not.(associated(temp%next))) exit
    temp => temp%next
  end do
  call create_main_atom_list()

  return

end subroutine change_atoms_plotting_status
!----------------------------------------------------------------------------------
subroutine change_Hatoms_plotting_status(status)
  logical, intent(in) :: status

  if (.not.(associated(atom_list_head))) return

  temp => atom_list_head
  do
    if (temp%atom%Symbol=="H " .or. temp%atom%Symbol=="D ") then
      temp%dead=status
    end if
    if (.not.(associated(temp%next))) exit
    temp => temp%next
  end do

  call create_main_atom_list()

  return

end subroutine change_Hatoms_plotting_status
!----------------------------------------------------------------------------------
subroutine search_atom_list_by_symbol(symbol,atomlist,success)
  character(len=2)                         , intent(in)   :: symbol
  type(gl_atom), dimension (:), allocatable, intent(out)  :: atomlist
  logical                                  , intent(out)  :: success
  type(gl_atom), dimension(400) :: templist
  integer :: count, i


   count=0
   success=.false.

   if (.not.(associated(atom_list_head))) return

     temp => atom_list_head
     do
       if (temp%atom%symbol==symbol) then
         success=.true.
         count=count+1
         templist(count)=temp%atom
       end if
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do

   allocate(atomlist(count))
   atomlist(:)=templist(1:count)

   return

end subroutine search_atom_list_by_symbol
!----------------------------------------------------------------------------------
subroutine empty_atom_list()

  do while(associated(atom_list_head))
    call pop_atom(1)
  end do

  if (associated(main_atom_list)) deallocate(main_atom_list)

  return

end subroutine empty_atom_list
!----------------------------------------------------------------------------------
subroutine write_atom_list(unit_no)
  integer, intent(in) :: unit_no
  ! Local variables
  character(len=10) :: dead_string=" "
  if (.not.(associated(atom_list_head))) return

  temp => atom_list_head
     do

       if (temp%dead) then
         dead_string=" NODISPLAY"
       else
         dead_string=" "
       end if

       if (temp%atom%colorfromtable) then
         write(unit=unit_no,fmt='(3a,3F9.5,a,F6.3,a)') "ATOM ", temp%atom%label, &
         temp%atom%symbol, temp%atom%xf, " RADIUS",&
         temp%atom%radius, dead_string
       else
         write(unit=unit_no,fmt='(3a,3F9.5,a,4F6.3,a,F6.3,a)') "ATOM ", temp%atom%label, &
         temp%atom%symbol, temp%atom%xf, " COLOR", temp%atom%color, " RADIUS",&
         temp%atom%radius, dead_string
       end if
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do
  return

end subroutine write_atom_list
!----------------------------------------------------------------------------------
subroutine generate_atom_label_array(atom_num,label_array)
  integer                                     , intent(out) :: atom_num
  character(len=13), dimension(:), allocatable, intent(out) :: label_array
  ! Local variables
  character(len=13), dimension(300)            :: temp_string=" "
  integer :: count


  if (allocated(label_array)) deallocate(label_array)

  if (.not.(associated(atom_list_head))) then
    allocate(label_array(1))
    atom_num=1
    label_array(1)="-New atom-"
    return
  end if

  count = 0
  temp => atom_list_head
     do
       count=count+1
       temp_string(count)=trim(temp%atom%symbol)//" "//trim(temp%atom%label)
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do


    atom_num=count+1
    allocate(label_array(atom_num))
    label_array(atom_num)="-New atom-"
    label_array(1:atom_num-1)=temp_string(1:atom_num-1)

  return

end subroutine generate_atom_label_array
!----------------------------------------------------------------------------------
subroutine get_atom_info(atom_id,label,symbol,xf,Biso,radius,dead,color,multip,colorfromtable)
    integer, intent(in) :: atom_id
    character(len=10),      intent(out)                 :: label
    character(len=2),       intent(out)                 :: symbol
    real         ,          intent(out)                 :: xf(3)
    real         ,          intent(out)                 :: Biso
    real         ,          intent(out)                 :: radius
    logical      ,          intent(out)                 :: dead
    real         ,          intent(out)                 :: color(4)
    integer      ,          intent(out)                 :: multip
    logical      ,          intent(out)                 :: colorfromtable
    integer :: id_count

    temp => atom_list_head                      ! Back to top
    id_count=0
    do
      if (.not. associated(temp)) exit          ! Tail of the list
      id_count=id_count+1
      if (id_count==atom_id) exit                ! We have a match
      previous_node => temp
      temp => temp%next
    end do

    if (associated(temp)) then
      label=temp%atom%label
      symbol=temp%atom%symbol
      xf=temp%atom%xf
      Biso=temp%atom%Biso
      radius=temp%atom%radius
      dead=temp%dead
      color=temp%atom%color
      multip=temp%atom%multip
      colorfromtable=temp%atom%colorfromtable
    else
      label="???"
      symbol="???"
      xf=0.0
      Biso=0.0
      radius=1.0
      dead=.false.
      color=0.0
      multip=0
      colorfromtable=.true.
    end if

    return

end subroutine get_atom_info
!----------------------------------------------------------------------------------
subroutine modify_atom(atom_id,label,symbol,xf,Biso,radius,dead,color,colorfromtable)
    integer, intent(in) :: atom_id
    character(len=10),      intent(in), optional :: label
    character(len=2),       intent(in), optional :: symbol
    real         ,          intent(in), optional :: xf(3)
    real         ,          intent(in), optional :: Biso
    real         ,          intent(in), optional :: radius
    logical      ,          intent(in), optional :: dead
    real         ,          intent(in), optional :: color(4)
    logical      ,          intent(in), optional :: colorfromtable
    integer :: id_count

    temp => atom_list_head                      ! Back to top
    id_count=0
    do
      if (.not. associated(temp)) exit          ! Tail of the list
      id_count=id_count+1
      if (id_count==atom_id) exit                ! We have a match
      previous_node => temp
      temp => temp%next
    end do

    if (associated(temp)) then
      if (present(label))  temp%atom%label=label
      if (present(symbol)) temp%atom%symbol=symbol

      if (present(dead)) then
        temp%dead=dead
        call create_main_atom_list()
      end if

      if (present(xf).or.present(radius).or.present(color).or.present(colorfromtable)) then
        if (present(radius)) temp%atom%radius=radius
        if (present(xf)) then
          temp%atom%xf=xf
          call calculate_eq_positions(temp%atom,current_space_group)
        end if
        if (present(color)) then
          temp%atom%color=color
          temp%atom%colorfromtable=.false.
        end if
        if(present(colorfromtable)) then
          temp%atom%color=get_color_from_symbol(temp%atom%symbol)
          temp%atom%colorfromtable=.true.
        end if
        call construct_atom_dlist(temp)
      end if
    end if

    return

end subroutine modify_atom
!----------------------------------------------------------------------------------
end module WCFGL_atom_tree

