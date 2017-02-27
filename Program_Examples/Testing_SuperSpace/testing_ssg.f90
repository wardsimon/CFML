    Program read_ssg_datafile
      use CFML_ssg_datafile
      use CFML_SuperSpaceGroups
      use CFML_Rational_Arithmetic

      implicit none

      integer :: iclass,nmod,i,j,k,m,multip,Dp, i1,i2,i3, ng
      character(len=20)  :: forma
      character(len=130) :: message,symb
      logical :: ok
      type(SuperSpaceGroup_Type) :: SSpaceGroup
      type(SuperSpaceGroup_Type),   dimension(200) :: spg
      type(rational),   dimension(:,:),allocatable :: Mat
      integer,          dimension(:,:),allocatable :: table
      character(len=40),dimension(:,:),allocatable :: matrix
      character(len=60) :: Operator_Symbol

      call Read_SSG(ok,message)
      if(.not. ok) then
        write(*,"(a)") "   !!! "//message//" !!!"
        stop
      end if
      !!
      !open(unit=1,file="class+group_pos.txt",status="replace",action="write")
      !write(1,"(a,i6)") "Bravais Classes positions in database, Number of classes ",m_ncl
      !write(1,"(10i7)") pos_class
      !write(1,"(a,i6)") "Space Group positions in database, Number of groups ",m_ngs
      !write(1,"(10i7)") pos_group
      !close(unit=1)
      !stop
      do
        !write(*,"(a)",advance="no")  " => Enter the dimension d of the matrix: "
        !read(*,*) Dp
        !if(Dp < 4) exit
        !if(allocated(mat)) deallocate(Mat)
        !allocate(Mat(Dp,Dp))
        !if(allocated(matrix)) deallocate(Matrix)
        !allocate(Matrix(Dp,Dp))
        !forma="( a8)"
        !write(forma(2:2),"(i1)") Dp
        !do
        !   write(*,"(a)",advance="no")  " => Enter the symbol of the operator: "
        !   read(*,"(a)") symb
        !   if(len_trim(symb) == 0) exit
        !   call Get_Mat_From_SSymSymb(Symb,Mat)
        !   matrix=print_rational(Mat)
        !   write(unit=*,fmt="(a)") "  Rational Matrix corresponding to "//trim(symb)
        !   do j=1,Dp
        !      write(unit=*,fmt=forma) (trim( Matrix(j,k))//" ",k=1,Dp)
        !   end do
        !   !Retransformation to a symbol
        !   write(unit=*,fmt="(a)")" "
        !   call Get_SSymSymb_from_Mat(Mat,Symb,"xyz")
        !   write(unit=*,fmt="(a)") "     xyz_type: "//trim(Symb)
        !   call Get_SSymSymb_from_Mat(Mat,Symb,"x1x2x3")
        !   write(unit=*,fmt="(a)") "  x1x2x3_type: "//trim(Symb)
        !   call Get_SSymSymb_from_Mat(Mat,Symb,"abc")
        !   write(unit=*,fmt="(a)") "     abc_type: "//trim(Symb)
        !end do

        write(*,"(a)",advance="no") " => Enter the number of the SSG: "
        read(*,*) m
        if(m <= 0) exit
        if(m > 16697) then
          write(*,"(a)") " => There are only 16697 superspace groups in the database! "
          cycle
        end if
        !Call Read_single_SSG(m,ok,Message)
        !Call Read_SSG(ok,Message)
        !if(.not. ok) then
        !  write(*,"(a)") "   !!! "//trim(message)//" !!!"
        !  stop
        !end if
        !write(*,"(3(a,i5))") " Group number:",igroup_number(m), " Bravais class:",igroup_class(m), "  Basic Group #:",igroup_spacegroup(m)
        !iclass=igroup_class(m)
        !nmod=iclass_nmod(iclass)
        !write(*,"(a,tr4,a)")  group_nlabel(m), group_label(m)
        !write(*,"(a,i3)") " Number of operators:", igroup_nops(m)
        !do k=1,igroup_nops(m)
        !  write(*,"(a,i3)") " Operator #",k
        !  forma="(  i4)"
        !  write(forma(2:3),"(i2)") nmod+4
        !  do i=1,nmod+4
        !    write(*,forma) igroup_ops(i,1:nmod+4,k,m)
        !    !write(*,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))
        !  end do
        !end do

        call Set_SSG_Reading_Database(m,SSpaceGroup,ok,message)
        if(.not. ok) then
          write(*,"(a)") "   !!! "//message//" !!!"
          stop
        end if

        call Write_SSG(SSpaceGroup,full=.true.)
        Dp=size(SSpaceGroup%SymOp(1)%Mat,dim=1)
        !if(SSpaceGroup%Centred == 2) then
        !  SSpaceGroup%SymOp(4)%Mat = -SSpaceGroup%SymOp(1)%Mat
        !  SSpaceGroup%SymOp(4)%Mat(Dp,Dp) = 1
        !end if

        !Generate subgroups 
        i1=0; i2=0; i3=0
        if(SSpaceGroup%Num_Lat > 1) then
          i1=SSpaceGroup%NumOps+1
          if(SSpaceGroup%Num_Lat > 2) then
            i2=2*SSpaceGroup%NumOps+1
            if(SSpaceGroup%Num_Lat > 3) then
              i3=3*SSpaceGroup%NumOps+1
            end if
          end if
        end if
        do i=1,200
          call Allocate_SSG_SymmOps(Dp-4,4*SSpaceGroup%multip,spg(i)%SymOp)
        end do
        ng=4
        do i=2,SSpaceGroup%numops-1
          spg(i)%SymOp(1)=SSpaceGroup%SymOp(1)
          spg(i)%SymOp(2)=SSpaceGroup%SymOp(i)
          spg(i)%SymOp(3)=SSpaceGroup%SymOp(i+1)
          spg(i)%SymOp(3)%Mat(Dp,Dp)=-1//1
          spg(i)%SymOp(4)=spg(i)%SymOp(1)  !Operator 4 is {1'|0001/2}
          spg(i)%SymOp(4)%Mat(Dp,Dp)=-1//1
          spg(i)%SymOp(4)%Mat(Dp-1,Dp)=1//2
          if(SSpaceGroup%Num_Lat > 1) then
           spg(i)%SymOp(5)=SSpaceGroup%SymOp(i1)
           ng=4+1
           if(SSpaceGroup%Num_Lat > 2) then
             spg(i)%SymOp(6)=SSpaceGroup%SymOp(i2)
             ng=4+2
             if(SSpaceGroup%Num_Lat > 3) then
               spg(i)%SymOp(7)=SSpaceGroup%SymOp(i3)
               ng=4+3
             end if
           end if
          end if
          Call Gen_Group(ng,spg(i)%SymOp,multip,table)
        !Writing of the rational operator matrices
          Write(*,"(/,2(a,i4))") " Number of generated operators for subgroup # ",i,":",multip
        !if(allocated(matrix)) deallocate(Matrix)
        !allocate(Matrix(Dp,Dp))
        !forma="( a8)"
        !write(forma(2:2),"(i1)") Dp
          do j=1,multip
            call Get_SSymSymb_from_Mat(spg(i)%SymOp(j)%Mat,Operator_Symbol,"xyz")
            write(unit=*,fmt="(a,i3,a)") "  Operator # ",j,"  "//trim(Operator_Symbol)
            !matrix=print_rational(SSpaceGroup%SymOp(i)%Mat)
            !write(unit=*,fmt="(a,i3)") "  Rational Operator #",i
            !do j=1,Dp
            !   write(unit=*,fmt=forma) (trim( Matrix(j,k))//" ",k=1,Dp)
            !end do
          end do
          write(unit=*,fmt="(/,a)") " Writing the multiplication table "
          do j=1,multip
            write(unit=*,fmt="(192i3)") table(j,1:multip)
          end do

        end do
        !write(*,"(a,i3)") " Reflection conditions: ",igroup_nconditions(m)
        !if(igroup_nconditions(m) > 0)then
        ! do k=1,igroup_nconditions(m)
        !   write(*,*)((igroup_condition1(i,j,k,m),i=1,nmod+3),j=1,nmod+3),(igroup_condition2(j,k,m),j=1,nmod+4)
        ! end do
        !end if
      end do

    End Program read_ssg_datafile