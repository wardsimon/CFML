  Module utilities
      implicit none
      private

      public :: File_To_FileList, Get_LogUnit

    !!----
    !!---- TYPE :: FILE_LIST_TYPE
    !!--..
    !!---- Type,public :: File_List_Type
    !!----    integer                                       :: nlines ! Number of lines in the file
    !!----    character(len=256), allocatable, dimension(:) :: line   ! Content of the lines
    !!---- End Type file_list_type
    !!----
    !!---- Updated: February - 2005, November 2012
    !!
    Type,public :: File_List_Type
       integer                                       :: nlines
       character(len=256), allocatable, dimension(:) :: line
    End Type File_List_Type

    Logical,          public  :: ERR_Util=.false.
    Character(len=80),public  :: ERR_Util_Mess=" "

    contains
    
    !!----
    !!---- Character Function U_Case(Text) Result (Mtext)
    !!----    character (len=*), intent(in) :: text   !  In -> String:"Input Line"
    !!----    character (len=len(text))     :: mtext  ! Out -> String:"INPUT LINE"
    !!----
    !!----    Conversion to upper case, text is not modified
    !!----
    !!---- Update: February - 2005
    !!
    Function U_Case(Text) Result (Mtext)
       !---- Argument ----!
       character (len=*), intent(in) :: text
       character (len=len(text))     :: mtext

       !---- Local variables ----!
       integer, parameter :: inc = ICHAR("A") - ICHAR("a")
       integer            :: leng, pos

       mtext=text
       leng=len_trim(mtext)
       do pos=1,leng
          if (mtext(pos:pos) >= "a" .and. mtext(pos:pos) <= "z")           &
              mtext(pos:pos) = CHAR ( ICHAR(mtext(pos:pos)) + inc )
       end do

       return
    End Function U_Case    

    !!----
    !!---- Subroutine File_To_FileList(File_dat,File_list)
    !!----   character(len=*),     intent( in) :: file_dat  !Input data file
    !!----   type(file_list_type), intent(out) :: file_list !File list structure
    !!----
    !!----    Charge an external file to an object of File_List_Type.
    !!----
    !!---- Update: August - 2008
    !!
    Subroutine File_To_FileList(File_dat,File_list)
       !---- Arguments ----!
       character(len=*),      intent( in) :: file_dat
       type(file_list_type),  intent(out) :: file_list

       !---- Local Variables ----!
       integer                           :: nlines

       !---- Number of Lines in the input file ----!
       call Number_Lines(trim(File_dat), nlines)

       if (nlines==0) then
          ERR_Util=.true.
          ERR_Util_Mess="The file "//trim(File_dat)//" contains nothing"
          return
       else
          file_list%nlines=nlines
          if (allocated(file_list%line)) deallocate(file_list%line)
          allocate(file_list%line(nlines))
          call reading_Lines(trim(File_dat),nlines,file_list%line)
       end if

       return
    End Subroutine File_To_FileList

    !!----
    !!---- Subroutine Get_LogUnit(lun)
    !!----   integer,     intent(out) :: lun !First logical unit available
    !!----
    !!----   Provides the number of the first logical unit that is not opened.
    !!----   Useful for getting a logical unit to a file that should be opened
    !!----   of the flight.
    !!----
    !!----   Update: February - 2005
    !!
    Subroutine Get_LogUnit(lun)
       !---- Arguments ----!
       integer,  intent(out) :: lun

       !---- Local variables ----!
       logical :: op
       integer, parameter :: max_iunits=500

       lun=1
       do
          inquire(unit=lun,opened=op)
          if (.not. op) exit
          lun=lun+1
          if (lun == max_iunits) then
             lun=-1
             exit
          end if
       end do

       return
    End Subroutine Get_LogUnit

    Subroutine Init_ERR_Util()
       ERR_Util=.false.
       ERR_Util_Mess=" "
    End Subroutine Init_ERR_Util

    !!----
    !!---- Subroutine Number_Lines(Filename,n, input_string)
    !!----    character(len=*), intent(in) :: Filename     !  In -> Name of the file
    !!----    integer        , intent(out) :: N            ! Out -> Number of lines in the file
    !!----    character(len=*), optional,intent(in) :: input_string   ! In -> String to exit
    !!----
    !!----    Return the number of lines contained in a file. The file will be opened and closed before
    !!----    returning to the calling unit.
    !!----    If 'input_string' is present, return the number of lines until 'input_string' is founded
    !!----    as first string in the line
    !!----    (example : input_string =='END' : avoid Q peaks in a SHELX file)
    !!----
    !!---- Update: February - 2005, March-2014 (removing the "opened" inquire, JRC)
    !!
    Subroutine Number_Lines(filename,n, input_string)
       !---- Arguments ----!
       character(len=*), intent(in)  :: filename
       integer,          intent(out) :: n
       character(len=*), optional, intent(in) :: input_string       ! TR may 2013

       !---- Local Variables ----!
       logical            :: info
       integer            :: lun,cond
       character (len=256):: read_line                             ! TR may 2013
       integer            :: long                                  ! TR may 2013

       !---- Init ----!
       info=.false.
       call get_logunit(lun)
       n=0
       cond=0

       if(present(input_string)) long = len_trim(input_string)    ! TR may 2013

       !---- Exist filename ? ----!
       inquire (file=filename,exist=info)
       if (.not. info) return

       open(unit=lun,file=filename, status="old",action="read", position="rewind")

       !---- Counting lines ----!
       do
          read(unit=lun,fmt="(a)",iostat=cond) read_line
          if (cond /= 0) exit
          read_line=adjustl(read_line)
          if(present(input_string)) then                                         ! TR may 2013
            if(u_case(read_line(1:long)) == u_case(input_string(1:long))) exit
          end if
          n=n+1
       end do

       close(unit=lun)

       return
    End Subroutine Number_Lines

    !!----
    !!---- Subroutine Reading_Lines(Filename,Nlines,Filevar)
    !!----    character(len= *), intent(in)                :: Filename   !  In -> Filename
    !!----    integer,           intent(in)                :: Nlines     !  In -> Number of lines to read
    !!----    character(len= *), dimension(:), intent(out) :: Filevar    ! Out -> String vector
    !!----
    !!----    Read nlines of the file and put the information on Filevar. The file is opened to read the
    !!----    lines and closed before returning to the calling unit.
    !!----
    !!---- Update: February - 2005, March-2014 (eliminating the "opened" inquire,JRC)
    !!
    Subroutine Reading_Lines(filename,nlines,filevar)
       !---- Arguments ----!
       character(len=*),               intent( in) :: filename
       integer,                        intent( in) :: nlines
       character(len=*), dimension(:), intent(out) :: filevar

       !---- Local Variables ----!
       logical :: info
       integer :: lun,i

       !---- Init ----!
       call init_ERR_Util()
       info=.false.
       call get_logunit(lun)

       !---- Exist filename ? ----!
       inquire (file=filename,exist=info)
       if (.not. info) then
          ERR_Util=.true.
          ERR_Util_Mess="The file"//trim(filename)//" does not exist "
          return
       end if

       open(unit=lun,file=filename, status="old",action="read", position="rewind")

       !---- Reading... ----!
       do i=1,nlines
          read(unit=lun,fmt="(a)") filevar(i)
       end do

       close(unit=lun)

       return
    End Subroutine Reading_Lines


  End Module utilities

  Module Qsort_Mod
    implicit none
    private

    public :: QSort, quick_sort, simple_sort

    Type, public :: group
       integer   :: order    ! original order of unsorted data
       real      :: value    ! values to be sorted by
    End Type group

    contains

    Recursive Subroutine QSort(a,na)
       type (group), dimension(nA), intent(in out) :: A
       integer,                     intent(in)     :: nA
       ! Local Variables
       integer :: left, right
       real :: random
       real :: pivot
       type (group) :: temp
       integer :: marker

       if (nA > 1) then

           call random_number(random)
           pivot = A(int(random*real(nA-1))+1)%value   ! random pivor (not best performance, but avoids worst-case)
           left = 0
           right = nA + 1

           do while (left < right)
               right = right - 1
               do while (A(right)%value > pivot)
                   right = right - 1
               end do
               left = left + 1
               do while (A(left)%value < pivot)
                   left = left + 1
               end do
               if (left < right) then
                   temp = A(left)
                   A(left) = A(right)
                   A(right) = temp
               end if
           end do

           if (left == right) then
               marker = left + 1
           else
               marker = left
           end if

           call QSort(A(:marker-1),marker-1)
           call QSort(A(marker:),nA-marker+1)

       end if

    End Subroutine QSort
    
    Recursive Subroutine quick_sort(list, order)    
       ! Quick sort routine from:
       ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
       ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
       ! Modified by Alan Miller to include an associated integer array which gives
       ! the positions of the elements in the original order.
       
       Real,    Dimension (:), Intent(In out)  :: list
       Integer, Dimension (:), Intent(Out)     :: order
       
       ! Local variable
       Integer :: i
       
       Do i = 1, SIZE(list)
         order(i) = i
       End Do
       
       Call quick_sort_1(1, SIZE(list))
       
       Contains
    
         Recursive Subroutine quick_sort_1(left_end, right_end)
         
            Integer, Intent(In) :: left_end, right_end
            
            !     Local Variables
            Integer             :: i, j, itemp
            Real                :: reference, temp
            Integer, Parameter  :: max_simple_sort_size = 6
            
            If (right_end < left_end + max_simple_sort_size) Then
              ! Use interchange sort for small lists
              Call interchange_sort(left_end, right_end)
            
            Else
              ! Use partition ("quick") sort
              reference = list((left_end + right_end)/2)
              i = left_end - 1; j = right_end + 1
            
              Do
                ! Scan list from left end until element >= reference is found
                Do
                  i = i + 1
                  If (list(i) >= reference) Exit
                End Do
                ! Scan list from right end until element <= reference is found
                Do
                  j = j - 1
                  If (list(j) <= reference) Exit
                End Do
            
            
                If (i < j) Then
                  ! Swap two out-of-order elements
                  temp = list(i); list(i) = list(j); list(j) = temp
                  itemp = order(i); order(i) = order(j); order(j) = itemp
                Else If (i == j) Then
                  i = i + 1
                  Exit
                Else
                  Exit
                End If
              End Do
            
              If (left_end < j) Call quick_sort_1(left_end, j)
              If (i < right_end) CalL quick_sort_1(i, right_end)
            End If
         
         End Subroutine quick_sort_1
         
         
         Subroutine interchange_sort(left_end, right_end)
         
            Integer, Intent(In) :: left_end, right_end
            
            !     Local variables
            Integer             :: i, j, itemp
            Real                :: temp
            
            Do i = left_end, right_end - 1
              Do j = i+1, right_end
                If (list(i) > list(j)) Then
                  temp = list(i); list(i) = list(j); list(j) = temp
                  itemp = order(i); order(i) = order(j); order(j) = itemp
                End If
              End Do
            End Do
         
         End Subroutine interchange_sort
    
    End Subroutine quick_sort
    
    Subroutine simple_sort(arr, order)
       real,    dimension(:), intent(in)  :: arr
       integer, dimension(:), intent(out) :: order            
       !     Local variables
       real, dimension(size(arr)) :: list
       Integer                    :: i, j, itemp
       Real                       :: temp
       
       list=arr
       Do i=1,size(list)
       	 order(i)=i
       End Do
       Do i = 1, size(list) - 1
         Do j = i+1, size(list)
           If (list(i) > list(j)) Then
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
           End If
         End Do
       End Do
    
    End Subroutine simple_sort    
  
  End Module Qsort_Mod
  
  Program reorder
    use utilities
    use qsort_mod
    implicit none
    logical :: arggiven
    integer :: narg, i, j, k, nsol, n, lun
    real    :: ini,fin
    type(file_list_type) :: file_lines
    character(len=120)   :: filcod
    integer,     dimension(:), allocatable :: init_line, end_line, ind
   !Type(group), dimension(:), allocatable :: merit
    real,        dimension(:), allocatable :: merit


    narg=COMMAND_ARGUMENT_COUNT()


    if(narg > 0) then
        call GET_COMMAND_ARGUMENT(1,filcod)
        arggiven=.true.
        i=index(filcod,'.ind',back=.true.)
        if( i /= 0) filcod=filcod(1:i-1)
    end if
    call cpu_time(ini)
    call File_To_FileList(trim(filcod)//".ind",file_lines)
    if(ERR_util) then
      write(unit=*,fmt="(a)") trim(ERR_util_Mess)
      stop
    end if
    ! First pass
    ! Look of keywords in the input file to determine the number of solutions
    nsol=0
    do i=1, file_lines%nlines
      if(index(file_lines%line(i),"DIRECT PARAMETERS :") /= 0) nsol=nsol+1
    end do
    if(nsol > 0) then
      allocate(init_line(nsol),end_line(nsol),ind(nsol),merit(nsol))
      !do i=1,nsol
      !	merit(i)%order=i
      !end do
    else
      write(unit=*,fmt="(a)") " => There's no solution!"
      stop
    end if
    ! Second pass
    ! Determine the initial and end line for each solution
    n=0
    do i=1, file_lines%nlines
      if(index(file_lines%line(i),"DIRECT PARAMETERS :") /= 0) then
       n=n+1
       init_line(n)= i-1  !precedent line containing the crystal system
      end if
      if(index(file_lines%line(i),"1.- M(") /= 0) then !Figure of merit
       j=index( file_lines%line(i),"=")
       read(unit=file_lines%line(i)(j+1:),fmt=*) merit(n) !%value
      end if
      if(index(file_lines%line(i),"             ---------") /= 0) then !end_line
        end_line(n)= i
      end if
    end do
    if( n /= nsol) then !Test
      write(unit=*,fmt="(a)") " => Warning! n /= nsol!"
    end if
    call cpu_time(fin)
    write(unit=*,fmt="(a,f10.6)")  " => Partial CPU time before ordering: ",fin-ini
    !Order by figure of merit in ascending order
    !call QSort(merit,nsol)
    !call quick_sort(merit,ind)
    call simple_sort(merit,ind)    
    !Write the *.ord file
    Call Get_LogUnit(lun)
    open(unit=lun,file=trim(filcod)//".ord", status="replace",action="write")
    !Copy the first lines up to "SELECTED OPTION"
    do i=1, file_lines%nlines
    	write(unit=lun,fmt="(a)") trim(file_lines%line(i))
    	if(index(file_lines%line(i),"SELECTED OPTION") /= 0) exit    	
    end do
    do i=1,nsol
    	!j=merit(nsol-i+1)%order  !this is the index of the hihger figure of merit
    	j=ind(nsol-i+1)   !this is the index of the hihger figure of merit
    	write(unit=lun,fmt="(/,2(a,i2),/)") "         SOLUTION NUMBER ",i," / ",nsol
    	do k=init_line(j),end_line(j)
     	   write(unit=lun,fmt="(a)") trim(file_lines%line(k))
   	  end do
    end do
    close(unit=lun)
    call cpu_time(fin)
    write(unit=*,fmt="(a,f10.6)")  " => Total CPU time : ",fin-ini
  End Program reorder

