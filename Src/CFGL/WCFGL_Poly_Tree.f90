module WCFGL_poly_tree

  use OPENGL
  use WCFGL_geometry
  use WCFGL_metrix
  use WCFGL_glpoly
  use WCFGL_chull3D

  implicit none
!------------------------------------------------------------------------------
  type :: poly_list
    character(len=10)                    :: label
    type(gl_poly), dimension(:), pointer :: poly
    logical                              :: dead
    integer                              :: dlist
    type(poly_list), pointer             :: next
  end type poly_list
!------------------------------------------------------------------------------
  type(poly_list) , pointer, private       :: poly_list_head => null(), &
                                              current_node =>   null(),   &
                                              previous_node =>  null(),  &
                                              temp => null()

  integer, pointer, save, private          :: main_poly_list => null()

  integer, parameter, private              :: eps=0.0001

  contains
!------------------------------------------------------------------------------
  subroutine draw_polys()

    if (new_box .or. new_cell .or. new_spacegroup) call update_all_polys()
    if (associated(main_poly_list)) then
      call glcalllist(main_poly_list)
    end if

    return

  end subroutine draw_polys
!-----------------------------------------------------------------------------------
  subroutine push_poly(label,dead,color,show_edges,edge_color,edge_radius)
    character(len=10), intent(in)           :: label
    logical          , intent(in)           :: dead
    real             , intent(in), optional :: color(4)
    logical          , intent(in), optional :: show_edges
    real             , intent(in), optional :: edge_color(3)
    real             , intent(in), optional :: edge_radius
    ! Local variables.
    logical :: success
    type(gl_poly), dimension(:), allocatable :: polybuffer
    logical :: show_edges_on
    real :: ecolor(3), eradius

    eradius=1.0
    if (present(edge_radius)) eradius=edge_radius

    show_edges_on=.true.
    if (present(show_edges))  show_edges_on=show_edges

    ecolor=(/0,0,0/)
    if (present(edge_color))  ecolor=edge_color

    if (present(color)) then
      call new_gl_polys(label,polybuffer,success,show_edges_on,ecolor,eradius,color)
    else
      call new_gl_polys(label,polybuffer,success,show_edges_on,ecolor,eradius)
    end if

    if (success) then
    if (.not.(associated(poly_list_head))) then
      allocate(poly_list_head)
      poly_list_head%next =>null()
      current_node => poly_list_head
    else
      allocate(current_node%next)
      current_node%next%next => null()
      current_node => current_node%next
    end if

      current_node%dead=dead
      current_node%label=label

      if (size(polybuffer)>0) then
      allocate(current_node%poly(size(polybuffer)))
      current_node%poly(:)=polybuffer(:)
      end if
      current_node%dlist=glgenlists(1)

      call construct_poly_dlist(current_node)

      call create_main_poly_dlist()

      if (allocated(polybuffer)) deallocate(polybuffer)
    end if
    return

  end subroutine push_poly
!----------------------------------------------------------------------------------
  subroutine construct_poly_dlist(temp)
    type(poly_list), intent(inout) :: temp
    ! Local variables
    integer :: ii, jj, kk, l, m
    real, dimension(3) :: trans  !pos, 

    if (.not.(associated(current_box))) call init_box()

      call glNewList(temp%dlist,GL_COMPILE) ! New display list.
      do ii=int(current_box(1))-1,int(current_box(2)) ! Loop on cell x
        do jj=int(current_box(3))-1,int(current_box(4)) ! Loop on cell y
          do kk=int(current_box(5))-1,int(current_box(6)) ! Loop on cell z
            do l=1, size(temp%poly)                       ! Loop on poly
              do m=1, size(temp%poly(l)%x_center, DIM=2)        ! Loop on eq positions
              if (in_limit(temp%poly(l)%x_center(:,m)+real((/ii,jj,kk/)))) then
                      trans=f2c(real((/ii,jj,kk/)))
                      call glPushMatrix()
                      call gltranslatef(trans(1),trans(2),trans(3))
                      call glCallList(temp%poly(l)%dlist(m))
                      call glPopMatrix()
              end if
              end do ! End on eq positions
            end do ! End om poly
          end do ! End on loop on cell z
        end do ! End loop on cell y
      end do ! End loop on cell x
      call glEndList()
        return
  end subroutine construct_poly_dlist
!----------------------------------------------------------------------------------
  subroutine pop_poly(poly_id)
    integer, intent(in) :: poly_id
    integer :: id_count, i,j

    temp => poly_list_head                      ! Back to top
    id_count=0
    do
      if (.not. associated(temp)) exit          ! Tail of the list
      id_count=id_count+1
      if (id_count==poly_id) exit                ! We have a match
      previous_node => temp
      temp => temp%next
    end do

    if (associated(temp)) then
      if (associated(temp,poly_list_head)) then   ! If id is the head of the list
         poly_list_head => temp%next
       else                                        ! If it is not
        previous_node%next => temp%next
      endif
        call gldeletelists(temp%dlist,1)          ! Delete the display list to free memory
        do i=1, size(temp%poly)
          if (associated(temp%poly(i)%xf))       deallocate(temp%poly(i)%xf)
          if (associated(temp%poly(i)%x_center)) deallocate(temp%poly(i)%x_center)
          do j=1, size(temp%poly(i)%dlist)
          call gldeletelists(temp%poly(i)%dlist(j),1) ! Make sure grpahic memory is cleared-up for
          end do
          if (associated(temp%poly(i)%dlist))    deallocate(temp%poly(i)%dlist)

        end do
        if (associated(temp%poly)) deallocate(temp%poly)                     ! Deallocate the pointer
        deallocate(temp)
    endif

    return

  end subroutine pop_poly
!----------------------------------------------------------------------------------
subroutine create_main_poly_dlist()

   if (.not.(associated(poly_list_head))) return

   if (.not.(associated(main_poly_list))) then
     allocate(main_poly_list)
     main_poly_list=glgenlists(1)
   end if
   call glnewlist(main_poly_list,GL_COMPILE)
     temp => poly_list_head
     do
       if (.not.(temp%dead)) call glcalllist(temp%dlist)
         if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do
   call glendlist()

   return

end subroutine create_main_poly_dlist
!----------------------------------------------------------------------------------
subroutine empty_poly_list()

  do while(associated(poly_list_head))
    call pop_poly(1)
  end do

  if (associated(main_poly_list)) deallocate(main_poly_list)

  return

end subroutine empty_poly_list
!----------------------------------------------------------------------------------
subroutine update_all_polys()

   if (.not.(associated(poly_list_head))) return

     temp => poly_list_head
     do
       call construct_poly_dlist(temp)
       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do

   return

end subroutine update_all_polys
!----------------------------------------------------------------------------------
subroutine write_poly_list(unit_no)
  integer, intent(in) :: unit_no
  ! Local variables
  character(len=10) :: dead_string=" "
  integer           :: i

  if (.not.(associated(poly_list_head))) return

  temp => poly_list_head
     do

       if (temp%dead) then
         dead_string=" NODISPLAY"
       else
         dead_string=" "
       end if

       do i=1, size(temp%poly)
       write(unit=unit_no,fmt='(3a,4F6.3,a)') "POLY ", temp%label &
         ," COLOR", temp%poly(i)%poly_color, dead_string
       end do

       if (.not.(associated(temp%next))) exit
       temp => temp%next
     end do
  return

end subroutine write_poly_list
!----------------------------------------------------------------------------------
end module WCFGL_poly_tree
