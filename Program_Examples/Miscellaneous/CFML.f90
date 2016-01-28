!!
!! Extracted from CrysFML
!!
 Module CFML_String_Utilities  !Just the procedures used in the program ILL_Pubs_WoS

    implicit none

    private

    !---- List of public functions ----!
    public :: Pack_String, String_Count, L_case

    !---- List of public subroutines ----!
    public :: lcase, Get_Separator_Pos, SString_Replace, Get_Substring_Positions

    !---- Definitions ----!

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!----
    !!---- Character Function L_Case(Text) Result (Mtext)
    !!----    character (len=*), intent(in) :: text   !  In -> String: "InPUT Line"
    !!----    character (len=len(text))     :: mtex   ! Out -> String: "input line"
    !!----
    !!----    Conversion to lower case, text is not modified
    !!----
    !!---- Update: February - 2005
    !!
    Function L_Case(Text) Result (Mtext)
       !---- Argument ----!
       character (len=*), intent(in) :: text
       character (len=len(text))     :: mtext

       !---- Local variables ----!
       integer, parameter :: inc = ICHAR("A") - ICHAR("a")
       integer            :: leng, pos

       mtext=text
       leng=len_trim(mtext)
       do pos=1,leng
          if (mtext(pos:pos) >= "A" .and. mtext(pos:pos) <= "Z")           &
              mtext(pos:pos) = CHAR ( ICHAR(mtext(pos:pos)) - inc )
       end do

       return
    End Function L_Case

    !!----
    !!---- Character Function Pack_String(String) Result (Strp)
    !!----    character (len=*), intent(in) :: String
    !!----    character (len=*)             :: Strp
    !!----
    !!----    Pack a string: the function provides a string without empty spaces
    !!----
    !!---- Update: February - 2005
    !!
    Function Pack_String(String) Result (Strp)
       !---- Argument ----!
       character (len=*), intent(in)    :: string
       character (len=len_trim(string)) :: strp

       !---- Local variables ----!
       integer ::  i,n

       n=0
       strp=" "
       do i=1,len(string)
          if (string(i:i) /= " ") then
             n=n+1
             strp(n:n)=string(i:i)
          end if
       end do

       return
    End Function Pack_String

    !!----
    !!---- Function String_Count(string,substr) result(coun)
    !!----    character(len=*), intent(in) :: string
    !!----    character(len=*), intent(in) :: substr
    !!----    integer                      :: coun
    !!----
    !!----  Function counting the number of times a substring appears in a string
    !!----
    !!---- Updated: May - 2014
    !!
    Function String_Count(string,substr) result(coun)
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: substr
      integer                      :: coun
      ! --- Local variables ---!
      character(len=len_trim(string)) :: cut_string
      integer :: i,lstr
      coun=0
      lstr=len_trim(substr)-1
      cut_string=string
      do
        i=index(cut_string,trim(substr))
        if (i == 0) exit
        coun=coun+1
        cut_string=cut_string(i+lstr:)
      end do
      return
    End Function String_Count


    !---------------------!
    !---- Subroutines ----!
    !---------------------!


    !!----
    !!---- Subroutine Get_Separator_Pos(line,car,pos,ncar)
    !!----   character(len=*),      intent(in)  :: line  ! In -> Input String
    !!----   character(len=1),      intent(in)  :: car   ! In -> Separator character
    !!----   integer, dimension(:), intent(out) :: pos   ! Out -> Vector with positions of "car" in "Line"
    !!----   integer,               intent(out) :: ncar  ! Out -> Number of appearance of "car" in "Line"
    !!----
    !!----    Determines the positions of the separator character "car" in string "Line" and generates
    !!----    the vector Pos containing the positions. The number of times the character "car" appears
    !!----    In "Line" is stored in "ncar".
    !!----    The separator "car" is not counted within substrings of "Line" that are written within
    !!----    quotes. The following example illustrates the functionning of the subroutine
    !!----
    !!----       !       12345678901234567890123456789012345678901234567890
    !!----        line =' 23, "List, of, authors", this book, year=1989'
    !!----
    !!----    A call like:  call Get_Separator_Pos(line,',',pos,ncar) provides
    !!----    ncar= 3
    !!----    pos= (/ 4, 25, 36, 0, ..../)
    !!----
    !!---- Updated: December 2009
    !!
    Subroutine Get_Separator_Pos(line,car,pos,ncar)
      character(len=*),      intent(in)  :: line
      character(len=1),      intent(in)  :: car
      integer, dimension(:), intent(out) :: pos
      integer,               intent(out) :: ncar
      integer :: i,j,k

      ncar=0
      j=0
      do i=1,len_trim(line)
        j=j+1
        if(line(j:j) == '"') then  !A chains of characters is found, advance up the the next "
          do k=1,len_trim(line)    !the character "car" is ignored if it is within " "
            j=j+1
            if(line(j:j) /= '"') cycle
            exit
          end do
        end if
        if(line(j:j) == car) then
          ncar=ncar+1
          pos(ncar)=j
        end if
      end do
      return
    End Subroutine Get_Separator_Pos

    !!----
    !!---- Subroutine Get_Substring_Positions(string,substr,pos,nsubs)
    !!----   character(len=*),      intent(in)  :: string   ! In -> Input String
    !!----   character(len=*),      intent(in)  :: substr   ! In -> Substring
    !!----   integer, dimension(:), intent(out) :: pos      ! Out -> Vector with positions of the firs character of "substr" in "String"
    !!----   integer,               intent(out) :: nsubs    ! Out -> Number of appearance of "substr" in "String"
    !!----
    !!----    Determines the positions of the substring "substr" in "String" and generates
    !!----    the vector Pos containing the positions of the first character of "substr" in "String".
    !!----    The number of times the "substr" appears in "String" is stored in "nsubs".
    !!----
    !!----     Updated: May 2014

    Subroutine Get_Substring_Positions(string,substr,pos,nsubs)
      character(len=*),      intent(in)  :: string
      character(len=*),      intent(in)  :: substr
      integer, dimension(:), intent(out) :: pos
      integer,               intent(out) :: nsubs
      integer :: i,j,lsubs

      nsubs=0
      lsubs=len_trim(substr)
      j=0
      do i=1,len_trim(string)
        j=j+1
        if(string(j:j+lsubs-1) == trim(substr)) then
          nsubs=nsubs+1
          pos(nsubs)=j
        end if
      end do
      return
    End Subroutine Get_Substring_Positions

    !!----
    !!---- Subroutine Lcase(Line)
    !!----    character(len=*), intent(in out) :: Line
    !!----
    !!----    Conversion to lower case. Line is modified
    !!----
    !!---- Updated: February - 2005
    !!
    Subroutine Lcase(line)
       !---- Argument ----!
       character (len=*), intent(in out) :: line

       line=l_case(line)

       return
    End Subroutine Lcase
    !!----
    !!----
    !!---- Subroutine SString_Replace(string, substr, rep_string,warning)
    !!----    character(len=*), intent(in out) :: string
    !!----    character(len=*), intent(in)     :: substr
    !!----    character(len=*), intent(in)     :: rep_string
    !!----    character(len=*), intent(out)    :: warning
    !!----
    !!----    Subroutine to replace a substring (substr) by another one (rep_string)
    !!----    within a given string (string). The original string is modified on output.
    !!----    If len_trim(warning) /= 0, one of the substrings will not be complete,
    !!----    it works as a warning or error condition without interrupting the
    !!----    procedure.
    !!----
    !!---- Updated: May - 2014
    !!
    Subroutine SString_Replace(string, substr, rep_string,warning)
      character(len=*), intent(in out) :: string
      character(len=*), intent(in)     :: substr
      character(len=*), intent(in)     :: rep_string
      character(len=*), intent(out)    :: warning
      ! --- Local variables ---!
      integer                                      :: i,j,lstr,ncount,nsubs,d,dmax
      integer,            dimension(:),allocatable :: pos
      character(len=1024),dimension(:),allocatable :: splitted_string

      lstr=len(substr)
      warning=" "
      i=index(rep_string,substr)
      if(i /= 0) then !Check if the substring to be replaced is contained in the replacing string
         !In such case the alternative short code doesn't work ... we have to use the longer analysis below
         ncount=String_Count(string,trim(substr))+1
         allocate(pos(ncount))
         allocate(splitted_string(ncount))
         call Get_Substring_Positions(string,substr,pos,nsubs)
         dmax=0
         do i=2,nsubs
           d=pos(i)-pos(i-1)
           if(d > dmax) dmax=d
         end do
         if(dmax > 1024) write(unit=warning,fmt="(a)") " => Warning! ... string too long to be fetch into the splitted_string"
         !Construct the splitted string
         j=1
         splitted_string(j)=string(1:pos(j)-1)
         do
           j=j+1
           if(j > nsubs) exit
           splitted_string(j)=string(pos(j-1)+lstr:pos(j)-1)
         end do
         splitted_string(ncount)=string(pos(nsubs)+lstr:)
         !Construct now the full string
         string=""
         do i=1,nsubs
           string=trim(string)//trim(splitted_string(i))//rep_string
         end do
         string=trim(string)//trim(splitted_string(ncount))

      else  !The following short code works easily when substr is not contained in rep_string

         do
           i=index(string,substr)
           if (i == 0) exit
           string=string(1:i-1)//rep_string//trim(string(i+lstr:))
         end do

      end if
      return
    End Subroutine SString_Replace

 End Module CFML_String_Utilities
