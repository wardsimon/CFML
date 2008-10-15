module WCFGL_chull3D
!------------------------------------------------
! Written by Laurent C.Chapon 
! December 2004.
! Compute the 3D convex hull based on an incremental algorithm 
! described by Joseph O'Rourke in "Computational Geometry in C"
! Translation of the algorithm and data structures in C. 
! Updated : 
!------------------------------------------------ 
  use OPENGL 
  use WCFGL_geometry 
  
  implicit none 
  
  type :: vertex
    real, dimension(3)    :: v            ! Coordinates of the points
    integer               :: num          ! Number of the point
    logical               :: onhull, mark ! True if point on hull, / processed
    type(edge)  , pointer :: duplicate => null() ! Point to duplicated edge
    type(vertex), pointer :: prev => null() , next => null()
  end type vertex
!------------------------------------------------  
  type :: edge
    type(face),   pointer :: adjface1 => null(), adjface2 => null() ! Adjacent faces
    type(vertex), pointer :: endpts1 => null(), endpts2 => null()   ! Endpoints
    type(face),   pointer :: newface => null()                      ! New face
    logical               :: delete                                 ! True if should be deleted
    type(edge),   pointer :: prev => null(), next=> null() 
  end type edge
!------------------------------------------------  
  type :: face
    type(edge),   pointer  :: edge1 => null(), edge2 => null(), edge3 => null() ! 3 edges
    type(vertex), pointer  :: vertex1 => null(),vertex2 => null(),vertex3 => null() ! 3 vertex
    logical                :: visible ! True if visible from new point 
    type(face),   pointer  :: prev => null(), next => null()
  end type face
!------------------------------------------------  
! Global variables definitions   
  type(vertex), pointer, private :: vertices => null()
  type(edge),   pointer, private :: edges    => null()
  type(face),   pointer, private :: faces    => null()
  logical, private, parameter    :: onhull=.true., removed=.true.,&
                                    visible=.true., processed=.true.
  real, parameter, private       :: eps=0.0001
  logical,            public     :: err_chull3D=.false.
  character(len=120), public     :: err_chull3D_mess="No error in chull3D"
!------------------------------------------------  

  interface  ADD
    module procedure addvertex
    module procedure addedge
    module procedure addface
  end interface

  interface  DELETE
    module procedure delvertex
    module procedure deledge
    module procedure delface
  end interface

  contains 
!------------------------------------------------ 
  subroutine read_vertices(points)
    real, dimension(:,:), intent(in) :: points
    ! Local variables
    integer :: num
    type(vertex), pointer :: v
    
    num=0
    err_chull3D=.false.
    
    if (size(points, DIM=1)/=3) then
      err_chull3D     =.true.
      err_chull3D_mess=" Points not of dimension 3" 
    end if   
    if (size(points, DIM=2)<4) then 
      err_chull3D     =.true.
      err_chull3D_mess=" At least 4 points are required for 3D convex Hull" 
    end if

    do num=1, size(points, DIM=2)
      call makenullvertex(v)
      v%v=points(:,num)
      v%num=num
    end do 
    
    return 
  end subroutine read_vertices
!------------------------------------------------ 
  subroutine doubletriangle()
    type(vertex), pointer :: v1, v2, v3, v4
    type(face), pointer   :: f1, f2
    real :: vol 
    
    f2=>null()
    v1=>vertices
    do 
      if (.not.(colinear(v1,v1%next,v1%next%next))) exit
      v1=>v1%next
      if (associated(v1,vertices)) then
        err_chull3D=.true.
        err_chull3D_mess="All points are colinear."
        return 
      end if 
    end do
    
    v2=>v1%next
    v3=>v1%next%next
    
    v1%mark=processed;v2%mark=processed;v3%mark=processed
    
    call makeface(v1,v2,v3,f2,f1)
    call makeface(v3,v2,v1,f1,f2)

    f1%edge1%adjface2 => f2;f1%edge2%adjface2 => f2;f1%edge3%adjface2 => f2
    f2%edge1%adjface2 => f1;f2%edge2%adjface2 => f1;f2%edge3%adjface2 => f1
	
    v4 => v3%next
    do 
      if (abs(volumesign(f1,v4))>0.00001) exit
      v4=>v4%next
      if (associated(v4,v1)) then
        err_chull3D=.true.
        err_chull3D_mess="The Hull is coplanar."
        return 
      end if 
    end do
    vertices => v4
    return
  end subroutine doubletriangle
!------------------------------------------------
  subroutine constructhull()
    type(vertex), pointer :: v, vnext
    real :: vol

      v => vertices
      do 
        vnext => v%next
        if (.not.(v%mark)) then 
          v%mark = PROCESSED
	      call addone(v)
	      call cleanup(vnext)
        end if
        v => vnext
        if (associated(v,vertices)) exit
      end do
      return
   
  end subroutine constructhull
!------------------------------------------------
  subroutine cleanup(v)
    type(vertex), pointer :: v 
    
    call cleanedges()
    call cleanfaces()
    call cleanvertices(v)
    return 
  end subroutine cleanup
!------------------------------------------------
  subroutine addone(p)
    type(vertex), pointer :: p
    type(face), pointer :: f
    type(edge), pointer :: e, temp
    real    :: vol
    logical :: vis, visible1, visible2
  

    vis=.false.

    f => faces
    do 
      vol = volumesign(f,p)
      if (vol<0) then
	    f%visible = visible  
	    vis = .true.  
      end if
      f => f%next
      if (associated(f,faces)) exit
    end do
    
    if (.not.(vis)) then 
      p%onhull = .not.(onhull)
      return
    end if
    
    e => edges
    do 
      temp => e%next
      
      if (associated(e%adjface1)) then 
        visible1=e%adjface1%visible
      else
        visible1=.false.
      end if
      
      if (associated(e%adjface2)) then 
        visible2=e%adjface2%visible
      else
        visible2=.false.
      end if

      if (visible1 .and. visible2) then
	    e%delete = removed
      else if (visible1 .or. visible2) then 
	    call makeconeface(e,p)
      end if
      e => temp
      if (associated(e,edges)) exit
   end do
   return 
  end subroutine addone
!------------------------------------------------
  subroutine makeconeface(e,p)
   type(edge),   pointer :: e
   type(vertex), pointer :: p
   type(edge),   pointer ::  new_edge1, new_edge2
   type(face),   pointer ::  new_face
      

        if (.not.(associated(e%endpts1%duplicate))) then                     
	      call makenulledge(new_edge1)
	      new_edge1%endpts1 => e%endpts1
	      new_edge1%endpts2 => p
	      e%endpts1%duplicate => new_edge1
        else
          new_edge1 => e%endpts1%duplicate
        end if 
        
        if (.not.(associated(e%endpts2%duplicate))) then
          call makenulledge(new_edge2)
	      new_edge2%endpts1 => e%endpts2
	      new_edge2%endpts2 => p
	      e%endpts2%duplicate => new_edge2
        else 
          new_edge2 => e%endpts2%duplicate
        end if 
        
      call makenullface(new_face)   
      new_face%edge1 => e
      new_face%edge2 => new_edge1
      new_face%edge3 => new_edge2
      call makeccw(new_face,e,p) 
      
      if (.not.(associated(new_edge1%adjface1))) then
        new_edge1%adjface1 => new_face
      else if (.not.(associated(new_edge1%adjface2))) then
        new_edge1%adjface2 => new_face   
      end if 
         
      if (.not.(associated(new_edge2%adjface1))) then
        new_edge2%adjface1 => new_face
      else if (.not.(associated(new_edge2%adjface2))) then
        new_edge2%adjface2 => new_face
      end if
      
    e%newface=> new_face
      
    return
     
  end subroutine makeconeface
!------------------------------------------------
  subroutine makeccw(f,e,p)
    type(face),   pointer :: f, fv
    type(edge),   pointer :: e, s
    type(vertex), pointer :: p
    
    if (associated(e%adjface1)) then
      if (e%adjface1%visible) then      
        fv => e%adjface1
      else 
        fv => e%adjface2
      end if
    else
        fv => e%adjface2
    end if
       
    if ((associated(fv%vertex1,e%endpts1))) then
      if (.not.(associated(fv%vertex2,e%endpts2))) then
        f%vertex1 => e%endpts2  
        f%vertex2 => e%endpts1    
      else                                
        f%vertex1 => e%endpts1  
        f%vertex2 => e%endpts2     
        call swap(s,f%edge2,f%edge3)
      end if
    end if
    
    if ((associated(fv%vertex2,e%endpts1))) then
      if (.not.(associated(fv%vertex3,e%endpts2))) then
        f%vertex1 => e%endpts2  
        f%vertex2 => e%endpts1    
      else                                
        f%vertex1 => e%endpts1  
        f%vertex2 => e%endpts2     
        call swap(s,f%edge2,f%edge3)
      end if
    end if

    if ((associated(fv%vertex3,e%endpts1))) then
      if (.not.(associated(fv%vertex1,e%endpts2))) then
        f%vertex1 => e%endpts2  
        f%vertex2 => e%endpts1    
      else                                
        f%vertex1 => e%endpts1  
        f%vertex2 => e%endpts2     
        call swap(s,f%edge2,f%edge3)
      end if
    end if

    f%vertex3 => p
    
    return   
 
  end subroutine makeccw
!------------------------------------------------
  subroutine cleanedges()
    type(edge), pointer :: e, t
		
    e => edges
    do
      if (associated(e%newface)) then 
	    if (e%adjface1%visible) then
	      e%adjface1 => e%newface 
	    else
          e%adjface2 => e%newface
        end if
	    e%newface => null()
      end if
      e => e%next
      if (associated(e,edges)) exit
    end do

    do while (associated(edges) .and. edges%delete)  
       e => edges
       call delete(edges,e)
    end do
   
    e => edges%next
    do 
      if (e%delete) then 
	    t => e
	    e => e%next
	    call delete(edges,t)
      else
        e => e%next
      end if 
      if (associated(e,edges)) exit
    end do
    
    return 
    
  end subroutine cleanedges
!------------------------------------------------
  subroutine cleanfaces()
    type(face), pointer :: f, t
    
    do while (associated(faces) .and. faces%visible)  
      f => faces
      call delete(faces,f)
    end do
   
    f => faces%next
    do 
      if (f%visible ) then 
	    t => f
	    f => f%next
	    call delete(faces,t)
      else
        f => f%next
      end if 
      if (associated(f,faces)) exit
    end do
   
    return
  end subroutine cleanfaces
!------------------------------------------------
  subroutine cleanvertices(pvnext)
    type(vertex), pointer :: pvnext
    type(edge),   pointer :: e
    type(vertex), pointer :: v, t
    
  
    e => edges
    
    do 
      e%endpts1%onhull = onhull
      e%endpts2%onhull = onhull
      e%endpts1%duplicate => null()
      e%endpts2%duplicate => null()
      e => e%next
      if (associated(e,edges)) exit
    end do
	
    do while (associated(vertices) .and. (vertices%mark) .and. (.not.(vertices%onhull)))  
      if (associated(v,pvnext)) pvnext => v%next
      v => vertices
      call delete(vertices,v)
    end do
     
    v => vertices%next
    do 
      if (v%mark .and. (.not.(v%onhull))) then    
	    t => v 
	    v => v%next
	    call delete(vertices,t)
      else
        v => v%next
      end if
      if (associated(v,vertices)) exit
    end do
	
    v => vertices
    do
      v%duplicate => null() 
      v%onhull    = .not.(onhull) 
      v => v%next
      if (associated(v,vertices)) exit
    end do
   
    return
  end subroutine cleanvertices
!------------------------------------------------
  subroutine makeface(v1,v2,v3,fold,f) 
    type(vertex), pointer :: v1,v2,v3
    type(face)  , pointer :: fold
    type(edge), pointer :: e1, e2, e3
    type(face), pointer :: f
    
    if (.not.(associated(fold))) then 
      call makenulledge(e1);call makenulledge(e2);call makenulledge(e3)
    else
      e1=>fold%edge3;e2=>fold%edge2;e3=>fold%edge1
    end if 
    
    e1%endpts1 => v1;e1%endpts2 => v2
    e2%endpts1 => v2;e2%endpts2 => v3
    e3%endpts1 => v3;e3%endpts2 => v1
	
    call makenullface(f)
    f%edge1 => e1;f%edge2 => e2;f%edge3 => e3
    f%vertex1 => v1;f%vertex2 => v2;f%vertex3 => v3
    e1%adjface1 => f;e2%adjface1 => f;e3%adjface1 => f
        
    return 
   end subroutine makeface
!------------------------------------------------
  subroutine swap(t,x,y)
    type(edge), pointer :: t,x,y
      t => x
      x => y
      y => t
    return
  end subroutine swap
!------------------------------------------------ 
  function colinear(v0,v1,v2) result(col)
    type(vertex), intent(in) :: v0,v1,v2
    logical                  :: col
    real, dimension(3) :: vect1, vect2 
    real :: a,b,c
    
    vect1=v1%v-v0%v
    vect2=v2%v-v0%v
    
    a=vect1(2)*vect2(3)-vect1(3)*vect2(2)
    b=vect1(3)*vect2(1)-vect1(1)*vect2(3)
    c=vect1(1)*vect2(2)-vect1(2)*vect2(1)
    
    if (abs(a)<eps .and. abs(b)<eps .and. abs(c)<eps) then 
      col=.true.
    else
      col=.false.
    end if
    
    return
  end function colinear 
!------------------------------------------------ 
  subroutine makenullvertex(p) 
    type(vertex), pointer :: p
    
    allocate(p)
    p%onhull=.not.(onhull)
    p%mark  =.not.(processed)
    call add(vertices,p)
    return   
  end subroutine makenullvertex
!------------------------------------------------  
  subroutine makenulledge(e)
    type(edge), pointer :: e
   
    allocate(e)
    e%delete = .not.(REMOVED)
    call add(edges,e)
    return
  end subroutine makenulledge
!------------------------------------------------ 
  subroutine makenullface(f)
    type(face), pointer :: f 
    
    allocate(f)
    f%visible = .not.(visible)
    call add(faces,f)
    return
  end subroutine makenullface
!------------------------------------------------ 
  subroutine addvertex(head,p)
    type(vertex), pointer :: head, p 
      if (associated(head)) then
        p%next=>head
        p%prev=>head%prev
        head%prev=>p
        p%prev%next=>p
      else
        head=>p
        head%prev=>p
        head%next=>p
      end if 
    return
  end subroutine addvertex
!------------------------------------------------  
  subroutine delvertex(head,p)
    type(vertex), pointer :: head, p 
      if (associated(head)) then
        if (associated(head,head%next)) then
          head => null()
        else if (associated(p,head)) then
          head=>head%next
        end if
        p%next%prev=>p%prev
        p%prev%next=>p%next
        p => null()
      end if
    return
  end subroutine delvertex
!------------------------------------------------ 
  subroutine addedge(head,p)
    type(edge), pointer :: head, p 
      if (associated(head)) then
        p%next=>head
        p%prev=>head%prev
        head%prev=>p
        p%prev%next=>p
      else
        head=>p
        head%prev=>p
        head%next=>p
      end if 
    return
  end subroutine addedge
!------------------------------------------------  
  subroutine deledge(head,p)
    type(edge), pointer :: head, p 
      if (associated(head)) then
        if (associated(head,head%next)) then
          head => null()
        else if (associated(p,head)) then
          head=>head%next
        end if
        p%next%prev=>p%prev
        p%prev%next=>p%next
        if (associated(p%adjface1)) then
        if (associated(p%adjface1%edge1,p)) then 
          nullify(p%adjface1%edge1)
        end if
        if (associated(p%adjface1%edge2,p)) then 
          nullify(p%adjface1%edge2)
        end if
        if (associated(p%adjface1%edge3,p)) then
          nullify(p%adjface1%edge3)
        end if
        end if
       if (associated(p%adjface2)) then
        if (associated(p%adjface2%edge1,p)) then 
          nullify(p%adjface2%edge1)
        end if
        if (associated(p%adjface2%edge2,p)) then 
          nullify(p%adjface2%edge2)
        end if
        if (associated(p%adjface2%edge3,p)) then 
          nullify(p%adjface2%edge3)
        end if
        end if
        p => null()
      end if
    return
  end subroutine deledge
!------------------------------------------------ 
  subroutine addface(head,p)
    type(face), pointer :: head, p 
      if (associated(head)) then
        p%next=>head
        p%prev=>head%prev
        head%prev=>p
        p%prev%next=>p
      else
        head=>p
        head%prev=>p
        head%next=>p
      end if 
    return
  end subroutine addface
!------------------------------------------------  
  subroutine delface(head,p)
    type(face), pointer :: head, p 
      if (associated(head)) then
        if (associated(head,head%next)) then
          head => null()
        else if (associated(p,head)) then
          head=>head%next
        end if
        p%next%prev=>p%prev
        p%prev%next=>p%next
        
        if (associated(p%edge1)) then
          if (associated(p%edge1%adjface1,p)) then
            nullify(p%edge1%adjface1)
          end if 
        end if
        if (associated(p%edge1)) then 
          if (associated(p%edge1%adjface2,p)) then 
            nullify(p%edge1%adjface2)
          end if 
        end if 
        if (associated(p%edge2)) then 
          if (associated(p%edge2%adjface1,p)) then 
            nullify(p%edge2%adjface1)
          end if 
        end if
        if (associated(p%edge2)) then 
          if (associated(p%edge2%adjface2,p)) then 
            nullify(p%edge2%adjface2)
          end if 
        end if 
        if (associated(p%edge3)) then
          if (associated(p%edge3%adjface1,p)) then 
            nullify(p%edge3%adjface1)
          end if 
        end if
        if (associated(p%edge3)) then
          if (associated(p%edge3%adjface2,p)) then 
            nullify(p%edge3%adjface2)
          end if 
        end if 
        p => null()
      end if
    return
  end subroutine delface
!------------------------------------------------ 
  function volumesign(f,p) result(vol)
    type(face),   pointer :: f
    type(vertex), pointer :: p
    real                  :: vol
    real, dimension(3) :: a, b, c
    
   a = f%vertex1%v-p%v;b = f%vertex2%v-p%v;c = f%vertex3%v-p%v

   vol=a(1)*(b(2)*c(3)-b(3)*c(2))+a(2)*(b(3)*c(1)-b(1)*c(3))+a(3)*(b(1)*c(2)-b(2)*c(1))
   
   return 
   
   end function volumesign 
!------------------------------------------------  
  subroutine construct_poly(npoints,fcolor,with_edges,ecolor,eradius)
    real, dimension(:,:),    intent(in)  :: npoints
    logical,                 intent(in)  :: with_edges
    real,                    intent(in)  :: fcolor(4), ecolor(3), eradius 
    type(face),  pointer :: temp 
    type(edge),  pointer :: te
    real, dimension(3) :: p1, p2, normal, p11, p12, p21, p22, n1, n2, pv
    
    call clean_up_hull()        ! Make sure Hull is cleaned-up before starting
    call read_vertices(npoints) ! Read the vertices
    if (err_chull3D) return     ! Error in reading vertices
    call doubletriangle()       ! Make the first triangle, double-faced
    call constructhull()        ! Construct Hull 
        
    temp=> faces
    do
      p1=temp%vertex2%v-temp%vertex1%v
      p2=temp%vertex3%v-temp%vertex2%v 
      normal=p1.vect.p2  
      call glmaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,fcolor)
      call glBegin(GL_TRIANGLES)
      call glNormal3fv(normal)
      call glVertex3fv(temp%vertex1%v)
      call glVertex3fv(temp%vertex2%v)
      call glVertex3fv(temp%vertex3%v)
      call glEnd()
      temp => temp%next
      if (associated(faces,temp)) exit
    end do
    
    if (with_edges) then
      call glpushattrib(GL_LIGHTING)
      call gldisable(GL_LIGHTING) 
      call glLineWidth(eradius)
      call glColor3fv(ecolor)
      te=> edges
      do
        if (associated(te%adjface1) .and. associated(te%adjface2)) then
        p11=te%adjface1%vertex2%v-te%adjface1%vertex1%v
      	p12=te%adjface1%vertex3%v-te%adjface1%vertex2%v
      	n1=.unit.(p11.vect.p12)
     	p21=te%adjface2%vertex2%v-te%adjface2%vertex1%v
     	p22=te%adjface2%vertex3%v-te%adjface2%vertex2%v
     	n2=.unit.(p21.vect.p22)
     	if (abs(n1(1)-n2(1))>0.01 .or. abs(n1(2)-n2(2))>0.01 .or. abs(n1(3)-n2(3))>0.01) then
          call glBegin(GL_LINES)
      	  call glVertex3fv(te%endpts1%v)
      	  call glVertex3fv(te%endpts2%v)
      	  call glEnd()
        end if
        end if 
        te => te%next
        if (associated(edges,te)) exit
      end do
    call glpopattrib()
    end if
    
    return 
  end subroutine construct_poly
!------------------------------------------------  
subroutine clean_up_hull()
type(vertex), pointer :: v
type(edge), pointer   :: e 
type(face), pointer   :: f
integer :: c1

    do while (associated(edges))  
       e => edges
       e%delete=.true.
       call delete(edges,e)
    end do
    
    do while (associated(faces))  
       f => faces
       call delete(faces,f)
    end do
    
    do while (associated(vertices))  
       v => vertices
       call delete(vertices,v)
    end do
    
	
	return
end subroutine clean_up_hull
  
end module WCFGL_chull3D
!------------------------------------------------  
