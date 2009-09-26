module WCFGL_bond_tree

  use OPENGL
  use WCFGL_geometry
  use WCFGL_metrix
  use WCFGL_glatom
  use WCFGL_glbond
  use WCFGL_objects_definition, only : bond_definition, init_bond
  use WCFGL_atom_tree

  implicit none
!------------------------------------------------------------------------------
  type :: bond_list
    character(len=10)                    :: label1, label2
    real                                 :: minvalue, maxvalue
    logical                              :: conn
    logical                              :: colorfromatoms
    real                                 :: radius
    type(gl_bond), dimension(:), pointer :: bond
    logical                              :: dead
    integer                              :: dlist
    type(bond_list), pointer             :: next
  end type bond_list
!------------------------------------------------------------------------------
  type(bond_list) , pointer, private       :: bond_list_head => null(), &
                                              current_node =>   null(),   &
                                              previous_node =>  null(),  &
                                              temp => null()

  integer, pointer, save, private          :: main_bond_list => null()

  integer, parameter, private              :: eps=0.01

  real, public                             :: maxval_poly  ! Maximum value to consider a bond for a polyhedron
                                                           ! Initialised to 1.9 in read_fst_file
                                                           ! and updated if there are CONN with higher values
                                                           ! for maximum distance
                                                           ! This remove the interference of Bond instructions with
                                                           ! polyhedra plotting
  contains
!------------------------------------------------------------------------------
  subroutine draw_bonds()

    if (new_box .or. new_cell .or. new_spacegroup) call update_all_bonds()

    if (associated(main_bond_list)) then
      call glcalllist(main_bond_list)
    end if

    return

  end subroutine draw_bonds
!-----------------------------------------------------------------------------------
  subroutine push_bond_by_label(label1,label2,minvalue,maxvalue,radius,dead,color)
    character(len=10), intent(in)           :: label1, label2
    real             , intent(in)           :: minvalue
    real             , intent(in)           :: maxvalue
    real             , intent(in)           :: radius
    logical          , intent(in)           :: dead
    real             , intent(in), optional :: color(4)
    ! Local variables.
    !real    :: rad
    !real    :: col(4)
    logical :: success1, success2
    type(gl_atom), dimension(:), allocatable :: list1, list2
    type(gl_bond), dimension(:), allocatable :: bondbuffer
    !integer :: i,j

    success1=.false.
    success2=.false.

    call search_atom_list_by_label(label1,list1,success1)
    call search_atom_list_by_label(label2,list2,success2)

    if (success1.and.success2) then
      if (present(color)) then
        call new_gl_bonds(list1,list2,minvalue,maxvalue,bondbuffer,.false.,color)
      else
        call new_gl_bonds(list1,list2,minvalue,maxvalue,bondbuffer,.false.)
      end if

      if (.not.(associated(bond_list_head))) then
        allocate(bond_list_head)
        bond_list_head%next =>null()
        current_node => bond_list_head
        if (.not.(associated(bond_definition))) call init_bond()
      else
        allocate(current_node%next)
        current_node%next%next => null()
        current_node => current_node%next
      end if

      current_node%radius=radius
      current_node%label1=label1
      current_node%label2=label2
      current_node%minvalue=minvalue
      current_node%maxvalue=maxvalue
      current_node%conn=.false.

      if (present(color)) then
        current_node%colorfromatoms=.true.
      else
        current_node%colorfromatoms=.false.
      end if

      allocate(current_node%bond(size(bondbuffer)))
      current_node%bond(:)=bondbuffer(:)
      current_node%dead=dead
      current_node%dlist=glgenlists(1)

      call construct_bond_dlist(current_node)
      call create_main_bond_dlist()

      if (allocated(bondbuffer)) deallocate(bondbuffer)
      if (allocated(list1))      deallocate(list1)
      if (allocated(list2))      deallocate(list2)

    end if

    return

  end subroutine push_bond_by_label
!-----------------------------------------------------------------------------------
  subroutine push_bond_by_symbol(symbol1,symbol2,minvalue,maxvalue,radius,dead,color)
    character(len=2) , intent(in)           :: symbol1, symbol2
    real             , intent(in)           :: minvalue
    real             , intent(in)           :: maxvalue
    real             , intent(in)           :: radius
    logical          , intent(in)           :: dead
    real             , intent(in), optional :: color(4)
    ! Local variables.
    !real    :: rad
    !real    :: col(4)
    logical :: success1, success2
    type(gl_atom), dimension(:), allocatable :: list1, list2
    type(gl_bond), dimension(:), allocatable :: bondbuffer
    !integer :: i,j

    success1=.false.
    success2=.false.

    call search_atom_list_by_symbol(symbol1,list1,success1)
    call search_atom_list_by_symbol(symbol2,list2,success2)

    if (success1.and.success2) then
      if (present(color)) then
        call new_gl_bonds(list1,list2,minvalue,maxvalue,bondbuffer,.true.,color)
      else
        call new_gl_bonds(list1,list2,minvalue,maxvalue,bondbuffer,.true.)
      end if


      if (.not.(associated(bond_list_head))) then
        allocate(bond_list_head)
        bond_list_head%next =>null()
        current_node => bond_list_head
        if (.not.(associated(bond_definition))) call init_bond()
      else
        allocate(current_node%next)
        current_node%next%next => null()
        current_node => current_node%next
      end if

      current_node%radius=radius
      current_node%label1=symbol1
      current_node%label2=symbol2
      current_node%minvalue=minvalue
      current_node%maxvalue=maxvalue
      current_node%conn=.true.

      if (present(color)) then
        current_node%colorfromatoms=.true.
      else
        current_node%colorfromatoms=.false.
      end if

      allocate(current_node%bond(size(bondbuffer)))
      current_node%bond(:)=bondbuffer(:)
      current_node%dead=dead
      current_node%dlist=glgenlists(1)
      call construct_bond_dlist(current_node)
      call create_main_bond_dlist()

      if (allocated(bondbuffer)) deallocate(bondbuffer)
      if (allocated(list1))      deallocate(list1)
      if (allocated(list2))      deallocate(list2)

    end if

    return

  end subroutine push_bond_by_symbol
!----------------------------------------------------------------------------------
  subroutine construct_bond_dlist(bond_pointer)
    type(bond_list), intent(inout) :: bond_pointer
    ! Local variables
    real     :: pos1(3), pos2(3), middle(3), spher(3)
    integer  :: i, j, k, l, eq_at, m

    if (.not.(associated(current_box))) call init_box()

    call glnewlist(bond_pointer%dlist,GL_COMPILE)

    do m=1, size(bond_pointer%bond)
      eq_at=size(bond_pointer%bond(m)%xf, DIM=2)

      do i=int(current_box(1))-1,int(current_box(2))
        do j=int(current_box(3))-1,int(current_box(4))
          do k=int(current_box(5))-1,int(current_box(6))
            do l=1, eq_at
              pos1=bond_pointer%bond(m)%xf(1:3,l)+real((/i,j,k/))
              pos2=bond_pointer%bond(m)%xf(4:6,l)+real((/i,j,k/))
              if ((in_limit(pos1).and.in_limit(pos2)).and.(.norm.(pos1-pos2)>eps)) then
                pos1=f2c(pos1)
                pos2=f2c(pos2)
                middle=(pos1+pos2)/2.0
                spher=cart2spher(pos2-pos1)
                call glpushmatrix()
                  call gltranslatef(pos1(1),pos1(2),pos1(3))
                  call glrotatef(spher(2),0.0,0.0,1.0)
                  call glrotatef(-spher(3),0.0,1.0,0.0)
                  call glmaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,bond_pointer%bond(m)%color1)
                  call glscalef(spher(1)/2.0,bond_pointer%radius,bond_pointer%radius)
                  call glcalllist(bond_definition)
                call glpopmatrix()
                call glpushmatrix()
                  call gltranslatef(middle(1),middle(2),middle(3))
                  call glrotatef(spher(2),0.0,0.0,1.0)
                  call glrotatef(-spher(3),0.0,1.0,0.0)
                  call glmaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,bond_pointer%bond(m)%color2)
                  call glscalef(spher(1)/2.0,bond_pointer%radius,bond_pointer%radius)
                  call glcalllist(bond_definition)
                call glpopmatrix()
              end if
            end do
          end do
        end do
      end do
    end do

    call glendlist()

    return

  end subroutine construct_bond_dlist
!----------------------------------------------------------------------------------
  subroutine pop_bond(bond_id)
    integer, intent(in) :: bond_id
    integer :: id_count, i

    temp => bond_list_head                      ! Back to top
    id_count=0
    do
      if (.not. associated(temp)) exit          ! Tail of the list
      id_count=id_count+1
      if (id_count==bond_id) exit                ! We have a match
      previous_node => temp
      temp => temp%next
    end do

    if (associated(temp)) then
      if (associated(temp,bond_list_head)) then   ! If id is the head of the list
         bond_list_head => temp%next
       else                                        ! If it is not
        previous_node%next => temp%next
      endif
      if (associated(temp,current_node)) current_node => previous_node
        call gldeletelists(temp%dlist,1)          ! Delete the display list to free memory
        do i=1, size(temp%bond)
          if (associated(temp%bond(i)%xf)) deallocate(temp%bond(i)%xf)
        end do
        deallocate(temp%bond)                     ! Deallocate the pointer
        deallocate(temp)
    endif

    return

  end subroutine pop_bond
!----------------------------------------------------------------------------------
subroutine create_main_bond_dlist()

   if (.not.(associated(bond_list_head))) return

   if (.not.(associated(main_bond_list))) then
     allocate(main_bond_list)
     main_bond_list=glgenlists(1)
   end if

   call glnewlist(main_bond_list,GL_COMPILE)
     temp => bond_list_head
     do
       if (.not.(temp%dead)) then
         call glcalllist(temp%dlist)
       end if
         if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do
   call glendlist()

   return

end subroutine create_main_bond_dlist
!----------------------------------------------------------------------------------
subroutine update_all_bonds()

   if (.not.(associated(bond_list_head))) return

     temp => bond_list_head
     do
       if (new_box) call construct_bond_dlist(temp)
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do

   return

end subroutine update_all_bonds
!----------------------------------------------------------------------------------
subroutine change_bonds_plotting_status(status)
  logical, intent(in) :: status

  if (.not.(associated(bond_list_head))) return

  temp => bond_list_head
  do
    temp%dead=status
    if (.not.(associated(temp%next))) exit
    temp => temp%next
  end do

  call create_main_bond_dlist()

  return

end subroutine change_bonds_plotting_status
!----------------------------------------------------------------------------------
subroutine empty_bond_list()

  do while(associated(bond_list_head))
    call pop_bond(1)
  end do

  if (associated(main_bond_list)) deallocate(main_bond_list)

  return

end subroutine empty_bond_list
!----------------------------------------------------------------------------------
subroutine search_bond_list_by_label(label,bondlist,success)
  character(len=10)                        , intent(in)   :: label
  type(gl_bond), dimension (:), allocatable, intent(out)  :: bondlist
  logical                                  , intent(out)  :: success
  type(gl_bond), dimension(1000) :: templist
  integer :: count, taille !, i


   count=0
   success=.false.

   if (.not.(associated(bond_list_head))) return

     temp => bond_list_head
     do
!       if (temp%label1==label) then
       if (temp%label1 == label .and. temp%label2 /= label .and. temp%minvalue < maxval_poly) then
         success=.true.
         taille=size(temp%bond)
         templist(count+1:count+taille)=temp%bond(:)
         count=count+taille
       end if
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do

   allocate(bondlist(count))
   bondlist(:)=templist(1:count)

   return

end subroutine search_bond_list_by_label
!----------------------------------------------------------------------------------
subroutine write_bond_list(unit_no)
  integer, intent(in) :: unit_no
  ! Local variables
  character(len=10) :: dead_string=" "
  character(len=5)  :: bond_string=" "

  if (.not.(associated(bond_list_head))) return

  temp => bond_list_head
     do

       if (temp%dead) then
         dead_string=" NODISPLAY"
       else
         dead_string=" "
       end if

       if (temp%conn) then
         bond_string="CONN "
       else
         bond_string="BOND "
       end if
       if (temp%colorfromatoms) then
         write(unit=unit_no,fmt='(3a,2F7.3,a,f7.3,a,4f6.3,a)') bond_string, temp%label1,&
                                                       temp%label2, temp%minvalue,&
                                                       temp%maxvalue," RADIUS", temp%radius,&
                                                       " COLOR", temp%bond(1)%color1, dead_string
       else
         write(unit=unit_no,fmt='(3a,2F6.3,a,f6.3,a)') bond_string, temp%label1,&
                                               temp%label2, temp%minvalue,&
                                               temp%maxvalue," RADIUS", temp%radius, dead_string
       end if
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do
  return

end subroutine write_bond_list
!----------------------------------------------------------------------------------
subroutine generate_bond_label_array(bond_num,label_array)
  integer                                     , intent(out) :: bond_num
  character(len=22), dimension(:), allocatable, intent(out) :: label_array
  ! Local variables
  character(len=22), dimension(300)            :: temp_string=" "
  integer :: count


  if (allocated(label_array)) deallocate(label_array)

  if (.not.(associated(bond_list_head))) then
    allocate(label_array(1))
    bond_num=1
    label_array(1)="-New bond-"
    return
  end if

  count = 0
  temp => bond_list_head
     do
       count=count+1
       temp_string(count)=trim(temp%label1)//"--"//trim(temp%label2)
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do


    bond_num=count+1
    allocate(label_array(bond_num))
    label_array(bond_num)="-New bond-"
    label_array(1:bond_num-1)=temp_string(1:bond_num-1)

  return

end subroutine generate_bond_label_array
!----------------------------------------------------------------------------------
subroutine get_bond_info(bond_id,label1,label2,minvalue,maxvalue,radius,conn,colorfromatoms,dead,color)
    integer,                intent(in)                  :: bond_id
    character(len=10),      intent(out)                 :: label1
    character(len=10),      intent(out)                 :: label2
    real         ,          intent(out)                 :: minvalue
    real         ,          intent(out)                 :: maxvalue
    real         ,          intent(out)                 :: radius
    logical      ,          intent(out)                 :: conn
    logical      ,          intent(out)                 :: colorfromatoms
    logical      ,          intent(out)                 :: dead
    real         ,          intent(out)                 :: color(4)
    integer :: id_count

    temp => bond_list_head                      ! Back to top
    id_count=0
    do
      if (.not. associated(temp)) exit          ! Tail of the list
      id_count=id_count+1
      if (id_count==bond_id) exit                ! We have a match
      previous_node => temp
      temp => temp%next
    end do

    if (associated(temp)) then
      label1=temp%label1
      label2=temp%label2
      minvalue=temp%minvalue
      maxvalue=temp%maxvalue
      radius=temp%radius
      conn=temp%conn
      colorfromatoms=temp%colorfromatoms
      dead=temp%dead
      color=temp%bond(1)%color1
    else
      label1="???"
      label2="???"
      minvalue=0.0
      maxvalue=2.0
      radius=1.0
      conn=.true.
      colorfromatoms=.true.
      dead=.false.
      color=(/1.0,1.0,1.0,1.0/)
    end if

    return

end subroutine get_bond_info
!----------------------------------------------------------------------------------
end module WCFGL_bond_tree
