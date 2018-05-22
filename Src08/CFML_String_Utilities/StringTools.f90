!!----
!!---- SUBMODULE CFML_String_Utilities
!!----
!!----
!!
Submodule (CFML_String_Utilities) StrTools
   !---- Parameters ----!
   implicit none
   
 Contains
   !!----
   !!---- FUNCTION EQUAL_SETS_TEXT
   !!----
   !!----    Determine if two sets of text lines are equal irrespective of the
   !!----    order of the lines. 
   !!----
   !!----    The function is true if the two sets of text
   !!----    have the same lines in whatever order.  Two lines are equal only
   !!----    if they have the same length and all their component characters
   !!----    are equal and in the same order.
   !!----
   !!---- Update: February - 2005
   !!
   Module Pure Function Equal_Sets_Text(str1,n1,str2,n2) result(Equal_sets_texto)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in) :: str1              ! Vector of String
      character(len=*), dimension(:), intent(in) :: str2              ! Vector of String
      integer,                        intent(in) :: n1                ! Number of lines on Text1
      integer,                        intent(in) :: n2                ! Number of lines on str2
      logical                                    :: Equal_sets_texto  !
       
      !---- Local variables ----!
      integer :: i,j
      logical :: info

      Equal_sets_texto=.false.

      if (n1 /= n2) return
      if (len(str1) /= len(str2)) return

      do i=1,n1
         info=.false.
         do j=1,n2
            if (str1(i) == str2(j)) then
               info=.true.
               exit
            end if
         end do
         if (.not. info) return
      end do

      Equal_sets_texto=.true.

      return
   End Function Equal_Sets_Text


   !!----
   !!---- FUNCTION L_CASE
   !!----
   !!----    Conversion to lower case, text is not modified
   !!----
   !!---- Update: February - 2005
   !!
   Module Pure Function L_Case(Str) Result (LStr)
      !---- Argument ----!
      character (len=*), intent(in) :: Str    ! Input String
      character (len=len(str))      :: LStr   ! lower case of Text

      !---- Local variables ----!
      integer, parameter :: INC = ICHAR("A") - ICHAR("a")
      integer            :: leng, pos

      lstr=str
      leng=len_trim(lstr)
      do pos=1,leng
         if (lstr(pos:pos) >= "A" .and. lstr(pos:pos) <= "Z")           &
             lstr(pos:pos) = CHAR ( ICHAR(lstr(pos:pos)) - INC )
      end do

      return
   End Function L_Case
   
   !!----
   !!---- FUNCTION PACK_STRINg(Str) Result (Strp)
   !!----
   !!----    Pack a string: the function provides a string without empty spaces
   !!----
   !!---- Update: February - 2005
   !!
   Module Pure Function Pack_String(Str) Result (Strp)
      !---- Argument ----!
      character (len=*), intent(in) :: str    ! Input String
      character (len=len_trim(str)) :: strp   ! Output string

      !---- Local variables ----!
      integer ::  i,n

      n=0
      strp=" "
      do i=1,len(str)
         if (str(i:i) /= " ") then
            n=n+1
            strp(n:n)=str(i:i)
         end if
      end do

      return
   End Function Pack_String

   !!----
   !!---- FUNCTION STRING_COUNT
   !!----
   !!----  counting the number of times a substring appears in a string
   !!----
   !!---- Updated: May - 2014
   !!
   Module Pure Function String_Count(str,substr) result(coun)
      !---- Arguments ----!
      character(len=*), intent(in) :: str       ! Input String
      character(len=*), intent(in) :: substr    ! Substring model
      integer                      :: coun      ! Number

      !---- Local variables ---!
      character(len=len_trim(str)) :: cut_string
      integer                      :: i,lstr

      !> Init
      coun=0

      lstr=len_trim(substr)-1
      cut_string=str

      do
         i=index(cut_string,trim(substr))
         if (i == 0) exit

         coun=coun+1
         cut_string=cut_string(i+lstr:)
      end do

      return
   End Function String_Count

   !!----
   !!---- Function Strip_String(string, to_strip) Result(striped_string)
   !!----
   !!---- Strip a string from a particular word
   !!----
   !!---- Update: January - 2010
   !!
   Module Pure Function Strip_String(str, to_strip) Result(striped_string)
      !---- Arguments----!
      character (len = *), intent(in) :: str            ! Input string
      character (len = *), intent(in) :: to_strip       ! Pattern
      character (len = len_trim(str)) :: striped_string

      !---- Local variables ----!
      integer :: i

      striped_string=trim(str)

      i=index(str,trim(to_strip),back=.true.)
      if (i > 0) striped_string=str(1:i-1)

      return
   End Function Strip_String
   
   !!----
   !!---- FUNCTION U_CASE
   !!----
   !!----    Conversion to upper case, text is not modified
   !!----
   !!---- Update: February - 2005
   !!
   Module Pure Function U_Case(Str) Result (UStr)
      !---- Argument ----!
      character (len=*), intent(in) :: Str   ! Input string
      character (len=len(Str))      :: UStr  ! Upper conversion

      !---- Local variables ----!
      integer, parameter :: INC = ICHAR("A") - ICHAR("a")
      integer            :: leng, pos

      UStr=Str
      leng=len_trim(UStr)
      do pos=1,leng
         if (UStr(pos:pos) >= "a" .and. UStr(pos:pos) <= "z")           &
             UStr(pos:pos) = CHAR ( ICHAR(UStr(pos:pos)) + INC )
      end do

      return
   End Function U_Case  
   
   !!----
   !!---- SUBROUTINE CUTST
   !!----
   !!----    Removes the first word of the input String.
   !!----    Provides (optionally) a string with the first word.
   !!----
   !!---- Update: February - 2005
   !!
   Module Subroutine Cutst(Str1,nlong1,Str2,nlong2)
      !---- Argument ----!
      character (len=*),           intent(in out) :: Str1     ! Input string / Out: string without the first word
      character (len=*), optional, intent(   out) :: Str2     ! The first word of String on Input
      integer,           optional, intent(   out) :: nlong1   ! Give the length of Str1 on Output
      integer,           optional, intent(   out) :: nlong2   ! Give the length of Str2 on Output

      !---- Local variables ----!
      integer  :: k,iniz1

      !---- Initializing variables ----!
      if (present(nlong1)) nlong1=0
      if (present(nlong2)) nlong2=0

      !---- Initializing to blank the directive ----!
      if (present(str2)) str2=" "

      !---- Elimination of possible blanks on the left ----!
      str1=adjustl(str1)
      if (len_trim(str1) <= 0) return

      k=len(str1)
      iniz1=index(str1," ")

      if (k ==1) then
         if (present(str2)) str2=str1
         if (present(nlong2)) nlong2=1
         str1=" "
      else
         if (iniz1 > 0) then
            if (present(str2))  str2=str1(1:iniz1-1)
            if (present(nlong2)) nlong2=len_trim(str1(1:iniz1-1))
            str1=str1(iniz1:)
         else
            if (present(str2))  str2=str1
            if (present(nlong2)) nlong2=len_trim(str1)
            str1=" "
         end if
      end if

      str1=adjustl(str1)
      if(present(nlong1)) nlong1=len_trim(str1)

      return
   End Subroutine Cutst
   
   !!----
   !!---- SUBROUTINE GET_DATETIME
   !!----
   !!---- Define a string containing the information about Date and Time
   !!----
   !!---- Update: January 2017
   !!
   Module Subroutine Get_DateTime(Str)
      !---- Argument ----!
      character(len=*), intent(out) :: Str  ! String containing the Date and Time

      !---- Local Variables ----!
      character (len=10) :: dat
      character (len=10) :: tim

      call date_and_time(date=dat,time=tim)

      Str="Date: "//dat(7:8)//"/"//dat(5:6)//"/"//dat(1:4)//      &
          "  Time: "//tim(1:2)//":"//tim(3:4)//":"//tim(5:10)

      return
   End Subroutine Get_DateTime

   !!----
   !!---- Subroutine Get_Dirname(Str, Directory)
   !!----
   !!---- Obtain the Path from a string
   !!----
   !!---- Update: January - 2010
   !!
   Module Subroutine Get_Dirname(Str,Directory)
      !---- Argument ----!
      Character (Len=*), Intent (In)  :: Str          ! String containing Path + Filename
      Character (Len=*), Intent (Out) :: Directory    ! Path

      !---- Local Variables ----!
      Integer :: I

      Directory=" "
      I = Index(str, Ops_Sep, Back = .True.)
      If (I > 0) Then
         Directory = str(1:I)
      End If

      Return
   End Subroutine Get_Dirname

   !!----
   !!---- SUBROUTINE GET_EXTENSION
   !!----
   !!---- Obtaining the extension of the Filename
   !!----
   !!---- Written: December - 2012
   !!
   Module Subroutine Get_Extension(filename, extension, dotted)
      !---- Arguments ----!
      character(len=*),  intent(in)  :: filename     ! Input filename
      character(len=*),  intent(out) :: extension    ! Extension of the file
      logical, optional, intent(in)  :: dotted       ! If True, the extension will be returned with a dot

      !---- Local Variables ----!
      integer :: idx
      logical :: dot

      !> Search for the last dot.
      idx = index(filename, '.', back=.true.)

      !> If no dot was found in the filename, then the file has no extension.
      if (idx == 0) then
          extension = " "
      else
          !> Handle the optional dotted argument.
          dot=.true.
          if (present(dotted)) dot=dotted

          if (.not. dot) idx = idx + 1

          ! The extension is set.
          extension = filename(idx:)
      end if

      return
   End Subroutine Get_Extension
   
   !!----
   !!---- Subroutine Get_Filename(Str, Filename)
   !!----
   !!---- Get Filename from a String containing in general Path+Filename
   !!----
   !!---- Update: June - 2011
   !!
   Module Subroutine Get_Filename(Str,Filename)
      !---- Argument ----!
      Character (Len=*), Intent (In)  :: Str       ! String containing Path + Filename
      Character (Len=*), Intent (Out) :: Filename  ! Filename

      !---- Local Variables ----!
      Integer :: i

      i = Index(str, OPS_SEP, Back = .True.)

      if (i > 0) then
          Filename = Str(I+1:)
      else
          Filename = trim(Str)
      end If

      Return
   End Subroutine Get_Filename
   
   !!----
   !!---- SUBROUTINE GET_SEPARATOR_POS
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
   !!---- Update: December 2009
   !!
   Module Subroutine Get_Separator_Pos(Str,car,pos,ncar)
      !---- Arguments ----!
      character(len=*),      intent(in)  :: Str   ! Inout String
      character(len=1),      intent(in)  :: car   ! Separator character
      integer, dimension(:), intent(out) :: pos   ! Vector with positions of "sep" in "Line"
      integer,               intent(out) :: ncar  ! Number of appearance of "sep" in "Line"
     
      !---- Local variables ----!
      integer :: i,j,k

      !> init
      ncar=0
      pos=0
      
      j=0
      do i=1,len_trim(str)
         j=j+1
         if (str(j:j) == '"') then  !A chains of characters is found, advance up the the next "
            do k=1,len_trim(str)    !the character "car" is ignored if it is within " "
               j=j+1
               if (str(j:j) /= '"') cycle
               exit
            end do
         end if
       
         if (str(j:j) == car) then
            ncar=ncar+1
            pos(ncar)=j
         end if
      end do
     
      return
   End Subroutine Get_Separator_Pos

   !!----
   !!---- SUBROUTINE GET_SUBSTRING_POSITIONS
   !!----
   !!----    Determines the positions of the substring "substr" in "String" and generates
   !!----    the vector Pos containing the positions of the first character of "substr" in "String".
   !!----    The number of times the "substr" appears in "String" is stored in "nsubs".
   !!----
   !!----     Updated: May 2014
   Module Subroutine Get_Substring_Positions(str,substr,pos,nsubs)
      !---- Arguments ----!
      character(len=*),      intent(in)  :: str         ! In -> Input String                                                         
      character(len=*),      intent(in)  :: substr      ! In -> Substring                                                            
      integer, dimension(:), intent(out) :: pos         ! Out -> Vector with positions of the firs character of "substr" in "String" 
      integer,               intent(out) :: nsubs       ! Out -> Number of appearance of "substr" in "String"                        
     
      !---- Local Variables ----!
      integer :: i,j,lsubs

      !> init
      nsubs=0
      pos=0
      
      lsubs=len_trim(substr)
      j=0
      do i=1,len_trim(str)
         j=j+1
         if (str(j:j+lsubs-1) == trim(substr)) then
            nsubs=nsubs+1
            pos(nsubs)=j
         end if
      end do
     
      return
   End Subroutine Get_Substring_Positions
   
   !!----
   !!---- SUBROUTINE GET_WORD
   !!----
   !!----    Determines the number of words (Ic) in the string "Line" and generates a
   !!----    character vector "Dire" with separated words.
   !!----    Control of errors is possible by inquiring the global variables Err_CFML
   !!----    and Err_CFML_Mess. The last modification allows to treat strings between
   !!----    quotes as a single word.
   !!----
   !!---- Update: July - 2011
   !!
   Module Subroutine Get_Word(Str,dire,ic)
      !---- Argument ----!
      character (len=*),                 intent ( in) :: Str   ! Input string
      character (len=*), dimension(:),   intent (out) :: dire  ! Vector of Words
      integer,                           intent (out) :: ic    ! Number of words

      !---- Local variables ----!
      character (len=len(Str))  :: line1,line2
      integer                   :: nlong2
      integer                   :: ndim, j

      !> Init
      call Clear_error()
      ic=0
      ndim=size(dire)
      line1=Str

      do
         line1=adjustl(line1)
         if(line1(1:1) == '"') then
            j=index(line1(2:),'"')
            if( j > 0) then
              line2=line1(2:j)
              nlong2=len_trim(line2)
              line1 = line1(j+2:)
            else
              Err_CFML=.true.
              Err_CFML_Flag=2
              Err_CFML_Msg="Non balanced quotes!"
              exit
            end if
         else
            call cutst(line1,str2=line2,nlong2=nlong2)
         end if
         if (nlong2 == 0) exit
         ic=ic+1
         if (ic > ndim) then
            Err_CFML=.true.
            Err_CFML_Flag=2
            Err_CFML_Msg="Dimension of DIRE exceeded"
            exit
         end if
         dire(ic)=line2(:nlong2)
      end do

      return
   End Subroutine Get_Word
   
   !!----
   !!---- SUBROUTINE LCASE(LINE)
   !!----
   !!----    Conversion to lower case. Line is modified
   !!----
   !!---- Update: February - 2005
   !!
   Module Subroutine Lcase(line)
      !---- Argument ----!
      character (len=*), intent(in out) :: line

      line=l_case(line)

      return
   End Subroutine Lcase
   
   !!----
   !!---- SUBROUTINE SUBSTRING_REPLACE
   !!----
   !!----    Subroutine to replace a substring (substr) by another one (rep_string)
   !!----    within a given string (string). The original string is modified on output.
   !!----    If len_trim(warning) /= 0, one of the substrings will not be complete,
   !!----    it works as a warning or error condition without interrupting the
   !!----    procedure.
   !!----
   !!---- Updated: May - 2014
   !!
   Module Subroutine SubString_Replace(string, substr, repstr, warning)
      !---- Arguments ----!
      character(len=*), intent(in out) :: string   ! Input/output string
      character(len=*), intent(in)     :: substr   ! Subtring to be replaced
      character(len=*), intent(in)     :: repstr   ! String for add
      character(len=*), intent(out)    :: warning  ! Message

      ! --- Local variables ---!
      integer                                      :: i,j,lstr,ncount,nsubs,d,dmax
      integer,            dimension(:),allocatable :: pos
      character(len=1024),dimension(:),allocatable :: splitted_string

      lstr=len(substr)
      warning=" "

      i=index(repstr,substr)
      if (i /= 0) then
         !> Check if the substring to be replaced is contained in the replacing string
         !> In such case the alternative short code doesn't work ... we have to use the longer analysis below
         ncount=String_Count(string,trim(substr))+1
         allocate(pos(ncount))
         allocate(splitted_string(ncount))
         call Get_Substring_Positions(string,substr,pos,nsubs)

         dmax=0
         do i=2,nsubs
            d=pos(i)-pos(i-1)
            if (d > dmax) dmax=d
         end do
         if (dmax > 1024) write(unit=warning,fmt="(a)") " => Warning! ... string too long to be fetch into the splitted_string"

         !>Construct the splitted string
         j=1
         splitted_string(j)=string(1:pos(j)-1)
         do
            j=j+1
            if (j > nsubs) exit
            splitted_string(j)=string(pos(j-1)+lstr:pos(j)-1)
         end do
         splitted_string(ncount)=string(pos(nsubs)+lstr:)

         !> Construct now the full string
         string=""
         do i=1,nsubs
            string=trim(string)//trim(splitted_string(i))//repstr
         end do
         string=trim(string)//trim(splitted_string(ncount))

      else
         !> The following short code works easily when substr is not contained in repstr
         do
            i=index(string,substr)
            if (i == 0) exit
            string=string(1:i-1)//repstr//trim(string(i+lstr:))
         end do
      end if

      return
   End Subroutine SubString_Replace
   
   !!----
   !!---- SUBROUTINE UCASE
   !!----
   !!----    Conversion to upper case. Line is modified
   !!----
   !!---- Update: February - 2005
   !!
   Module Subroutine Ucase(Str)
      !---- Argument ----!
      character (len=*), intent(in out) :: Str

      str=u_case(str)

      return
   End Subroutine Ucase
   
End Submodule StrTools   