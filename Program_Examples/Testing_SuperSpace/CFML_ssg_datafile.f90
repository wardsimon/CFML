  Module CFML_ssg_datafile
      ! Read data for (3+d)=D-dimensional superspace groups (d=1,2,3)
      !
    implicit none
    !private

    public  :: Read_SSG, Read_single_SSG
    logical, public            :: ssg_database_allocated=.false.
    !logical, public            :: Err_ssg_database=.false. 
    !character,len(80),public   :: Err_ssg_database_mess
    integer, parameter, public :: m_cen=16, m_ncl=322, m_ngs=16697, m_ops=48, &
                                  m_dim=6, m_cond=50, m_qv=3 !D=3+d (d=3)
    !
    integer                              :: nclasses=0          ! number of Bravais classes
    integer,   dimension(:), allocatable :: iclass_nmod         ! for each Bravais class: number of modulation q vectors
    integer,   dimension(:), allocatable :: iclass_number       ! class number
    integer,   dimension(:), allocatable :: iclass_spacegroup   ! basic space group of lattice
    integer,   dimension(:), allocatable :: iclass_nstars       ! number of different stars of q
    integer, dimension(:,:), allocatable :: iclass_nmodstar     ! number of modulation q vectors for each star

    character(len=5),  dimension(:)       , allocatable :: class_nlabel       ! class number label: 1.1, 1.2, 1.3, etc.
    character(len=50), dimension(:)       , allocatable :: class_label        !   class label
    integer,           dimension(:,:,:,:) , allocatable :: iclass_qvec        !   q vectors
    integer,           dimension(:)       , allocatable :: iclass_ncentering  ! number of centering translations
    integer,           dimension(:,:,:)   , allocatable :: iclass_centering   ! centering translations: D=3+d integers followed by a common denominator
    integer                                             :: ngroups            ! number of superspace groups
    integer,           dimension(:)       , allocatable :: igroup_number      ! for each superspace group,  group number
    integer,           dimension(:)       , allocatable :: igroup_class       ! Bravais class
    integer,           dimension(:)       , allocatable :: igroup_spacegroup  ! Basic space group
    character(len=13), dimension(:)       , allocatable :: group_nlabel       !   group number label: 1.1.1.1, 2,1,1,1, etc.
    character(len=60), dimension(:)       , allocatable :: group_label        !   group label
    integer,           dimension(:)       , allocatable :: igroup_nops        !   number of operators
    !   (D+1)x(D+1) augmented matrix for each operator in supercentered setting
    !   common denominator in element (D+1,D+1)
    integer, dimension(:,:,:,:), allocatable :: igroup_ops 
    integer, dimension(:)      , allocatable :: igroup_nconditions  ! number of reflection conditions
    integer, dimension(:,:,:,:), allocatable :: igroup_condition1   ! matrix representation of righthand side
    integer, dimension(:,:,:)  , allocatable :: igroup_condition2   ! vector representation of lefthand side
    integer, dimension(:)      , allocatable :: pos_group           !position in the file of the groups
    integer, dimension(:)      , allocatable :: pos_class           !position in the file of the Bravais classes
    integer, dimension(:)      , allocatable :: pos_all             !position in the file of all

    contains
      	
      Subroutine Allocate_SSG_database()
        if(ssg_database_allocated) return
        if(.not. allocated(iclass_nmod))        allocate(iclass_nmod (m_ncl))                         ; iclass_nmod=0      
        if(.not. allocated(iclass_number))      allocate(iclass_number(m_ncl))                        ; iclass_number=0        
        if(.not. allocated(iclass_spacegroup))  allocate(iclass_spacegroup(m_ncl))                    ; iclass_spacegroup=0
        if(.not. allocated(iclass_nstars))      allocate(iclass_nstars(m_ncl))                        ; iclass_nstars=0        
        if(.not. allocated(iclass_nmodstar))    allocate(iclass_nmodstar(3,m_ncl))                    ; iclass_nmodstar=0    
        if(.not. allocated(class_nlabel))       allocate(class_nlabel(m_ncl))                         ; class_nlabel=" " 
        if(.not. allocated(class_label ))       allocate(class_label (m_ncl))                         ; class_label =" "  
        if(.not. allocated(iclass_qvec ))       allocate(iclass_qvec (3,3,m_qv,m_ncl))                ; iclass_qvec =0    
        if(.not. allocated(iclass_ncentering))  allocate(iclass_ncentering(m_ncl))                    ; iclass_ncentering=0 
        if(.not. allocated(iclass_centering ))  allocate(iclass_centering(m_dim+1,m_cen,m_ncl))       ; iclass_centering =0  
        
        if(.not. allocated(igroup_number    ))  allocate(igroup_number(m_ngs))                        ; igroup_class=0                 
        if(.not. allocated(igroup_class     ))  allocate(igroup_class(m_ngs))                         ; igroup_class=0                 
        if(.not. allocated(igroup_spacegroup))  allocate(igroup_spacegroup(m_ngs))                    ; igroup_spacegroup=0            
        if(.not. allocated(group_nlabel     ))  allocate(group_nlabel(m_ngs))                         ; group_nlabel=" " 
        if(.not. allocated(group_label      ))  allocate(group_label(m_ngs))                          ; group_label=" "  
        if(.not. allocated(igroup_nops      ))  allocate(igroup_nops(m_ngs))                          ; igroup_nops=0      
        if(.not. allocated(igroup_ops        )) allocate(igroup_ops(m_dim+1,m_dim+1,m_ops,m_ngs))     ; igroup_ops        =0
        if(.not. allocated(igroup_nconditions)) allocate(igroup_nconditions(m_ngs))                   ; igroup_nconditions=0  
        if(.not. allocated(igroup_condition1 )) allocate(igroup_condition1(m_dim,m_dim,m_cond,m_ngs)) ; igroup_condition1 =0  
        if(.not. allocated(igroup_condition2 )) allocate(igroup_condition2(m_dim+1,m_cond,m_ngs))     ; igroup_condition2 =0  
        if(.not. allocated(pos_group ))         allocate(pos_group(m_ngs))                            ; pos_group=0       
        if(.not. allocated(pos_class ))         allocate(pos_class(m_ncl))                            ; pos_class=0             
        if(.not. allocated(pos_all   ))         allocate(pos_all(m_ncl+m_ngs))                        ; pos_all  =0           
        ssg_database_allocated=.true.
      End Subroutine Allocate_SSG_database
      
      Subroutine deAllocate_SSG_database()
        if( .not. ssg_database_allocated) return
        if( allocated(iclass_nmod))        deallocate(iclass_nmod)       
        if( allocated(iclass_number))      deallocate(iclass_number)     
        if( allocated(iclass_spacegroup))  deallocate(iclass_spacegroup) 
        if( allocated(iclass_nstars))      deallocate(iclass_nstars)     
        if( allocated(iclass_nmodstar))    deallocate(iclass_nmodstar)   
        if( allocated(class_nlabel))       deallocate(class_nlabel)      
        if( allocated(class_label ))       deallocate(class_label )      
        if( allocated(iclass_qvec ))       deallocate(iclass_qvec)       
        if( allocated(iclass_ncentering))  deallocate(iclass_ncentering) 
        if( allocated(iclass_centering ))  deallocate(iclass_centering)  
        if( allocated(igroup_number))      deallocate(igroup_number)      
        if( allocated(igroup_class))       deallocate(igroup_class)      
        if( allocated(igroup_spacegroup))  deallocate(igroup_spacegroup) 
        if( allocated(group_nlabel))       deallocate(group_nlabel)      
        if( allocated(group_label ))       deallocate(group_label)       
        if( allocated(igroup_nops ))       deallocate(igroup_nops)       
        if( allocated(igroup_ops))         deallocate(igroup_ops)        
        if( allocated(igroup_nconditions)) deallocate(igroup_nconditions)
        if( allocated(igroup_condition1 )) deallocate(igroup_condition1) 
        if( allocated(igroup_condition2 )) deallocate(igroup_condition2) 
        if( allocated(pos_group))          deallocate(pos_group)         
        if( allocated(pos_class))          deallocate(pos_class)         
        if( allocated(pos_all  ))          deallocate(pos_all)           
        ssg_database_allocated=.false.
      End Subroutine deAllocate_SSG_database

      Subroutine Read_SSG(ok,mess)
        logical,          intent(out) :: ok
        character(len=*), intent(out) :: mess
        !
        integer :: i,j,k,m,n,imax,nmod,iclass
        integer :: i_db, ier,L
        character(len=512) :: ssg_file
        character(len=4)   :: line
         !imax=0
         if(.not. ssg_database_allocated) then
         	 call Allocate_SSG_database()
         end if
         i_db=1; ier=0
         ok=.true.
         mess=" "
         ! open data file
         ssg_file='ssg_datafile.txt'
         open(unit=i_db,file=ssg_file,status='old',action='read',position='rewind',iostat=ier)
         if(ier /= 0) then
           ok=.false.
           mess= 'Error opening the database file: '//trim(ssg_file)
           return
         end if   
         L=0
         do i=1,2526
           read(i_db,"(a)") line
           if(line(1:1) == '"') then
             L=L+1
             pos_class(L)= i-1
           end if
         end do
         L=0

         do i=2527,300000
           read(i_db,"(a)",iostat=ier) line
           if (ier /= 0) exit
           if(line(1:1) == '"') then
             L=L+1
             pos_group(L)= i-1
           end if
         end do
         rewind(i_db)
         ! skip heading
         read(i_db,*)
         read(i_db,*)
         read(i_db,*)
         ! read number of Bravais classes
         read(i_db,*) nclasses
         ! read each Bravais class
         do m=1,nclasses
           read(i_db,*)n,iclass_nmod(m),iclass_number(m), iclass_spacegroup(m),iclass_nstars(m), &
                     (iclass_nmodstar(i,m),i=1,iclass_nstars(m))
           nmod=iclass_nmod(m)
           if(n /= m) then
             ok=.false.
             write(mess,"(a,i3)") 'Error in ssg_datafile @reading Bravais class #: ',m
             return
           end if

           read(i_db,*)class_nlabel(m),class_label(m)
           read(i_db,*)(((iclass_qvec(i,j,k,m),i=1,3),j=1,3),k=1,nmod)
           read(i_db,*)iclass_ncentering(m)
           read(i_db,*)((iclass_centering(i,j,m),i=1,nmod+4),j=1,iclass_ncentering(m))
         end do
 
         ! read number of superspace groups
         read(i_db,*)ngroups
         ! read each superspace group
         do m=1,ngroups
           !write(6,'(i5)')m
           read(i_db,*)n,igroup_number(m),igroup_class(m),igroup_spacegroup(m)
           if(n /= m)then
             ok=.false.
             write(mess,"(a,i3)") 'Error in ssg_datafile @reading group#: ',m
             return
           end if
           iclass=igroup_class(m)
           nmod=iclass_nmod(iclass)
           read(i_db,*)group_nlabel(m),group_label(m)
           read(i_db,*)igroup_nops(m)
           read(i_db,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))
           read(i_db,*)igroup_nconditions(m)
           if(igroup_nconditions(m) > 0)then
               read(i_db,*)(((igroup_condition1(i,j,k,m),i=1,nmod+3),j=1,nmod+3),(igroup_condition2(j,k,m),j=1,nmod+4), &
                          k=1,igroup_nconditions(m))
           end if
         end do
         ssg_database_allocated=.true.
      End Subroutine Read_SSG

      Subroutine Read_single_SSG(num,ok,Mess)
        integer,          intent(in)  :: num
        logical,          intent(out) :: ok
        character(len=*), intent(out) :: mess
        !
        integer :: i,j,k,n,m,i_pos,n_skip,nmod,i_db,ier,iclass
        character(len=512) ssg_file,pos_file
         ok=.true.
         mess=" "
         ! open data file
         ssg_file='ssg_datafile.txt'
         open(newunit=i_db,file=ssg_file,status='old',action='read',position='rewind',iostat=ier)
         if(ier /= 0) then
           ok=.false.
           mess= 'Error opening the database file: '//trim(ssg_file)
           return
         end if
         pos_file='class+group_pos.txt'
         open(newunit=i_pos,file=pos_file,status='old',action='read',position='rewind',iostat=ier)
         if(ier /= 0) then
           ok=.false.
           mess= 'Error opening the database file: '//trim(pos_file)
           return
         end if
         read(unit=i_pos,fmt=*) !skip class line
         read(unit=i_pos,fmt=*) pos_class
         read(unit=i_pos,fmt=*) !skip group line
         read(unit=i_pos,fmt=*) pos_group
         close(unit=i_pos)
         read(i_db,*)
         read(i_db,*)
         read(i_db,*)
         ! read number of Bravais classes
         read(i_db,*) nclasses
         ! read each Bravais class
         do m=1,nclasses
           read(i_db,*)n,iclass_nmod(m),iclass_number(m), iclass_spacegroup(m),iclass_nstars(m), &
                     (iclass_nmodstar(i,m),i=1,iclass_nstars(m))
           nmod=iclass_nmod(m)
           read(i_db,*)class_nlabel(m),class_label(m)
           read(i_db,*)(((iclass_qvec(i,j,k,m),i=1,3),j=1,3),k=1,nmod)
           read(i_db,*)iclass_ncentering(m)
           read(i_db,*)((iclass_centering(i,j,m),i=1,nmod+4),j=1,iclass_ncentering(m))
         end do
         rewind(i_db)
         !write(*,"(10i8)") pos_group
         n_skip=pos_group(num)-1
         !write(*,"(a,i12)") "Skipping ",n_skip
         do i=1,n_skip
           read(unit=i_db,fmt=*)
         end do
         m=num
         read(i_db,*)n,igroup_number(m),igroup_class(m),igroup_spacegroup(m)
         if(n /= m)then
           ok=.false.
           write(mess,"(a,2i5)") 'Error in ssg_datafile @reading group#: ',m,n
           return
         end if
         iclass=igroup_class(m)
         nmod=iclass_nmod(iclass)
         read(i_db,*)group_nlabel(m),group_label(m)
         read(i_db,*)igroup_nops(m)
         read(i_db,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))  
         call deAllocate_SSG_database()
      End Subroutine Read_single_SSG

  End Module CFML_ssg_datafile
