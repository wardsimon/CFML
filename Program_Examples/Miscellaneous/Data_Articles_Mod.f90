!!---- Module containing types and procedures to generate ISI Web strings for searching articles
!!---- Author: Juan Rodriguez-Carvajal (ILL)

  Module Data_Articles_Mod
   !Use CFML_String_Utilities, only: Pack_String, get_separator_pos, SString_Replace
   Use CFML_procedures, only: Pack_String, get_separator_pos, SString_Replace
   implicit none
   private
   public :: ISI_string
   character(len=1), parameter, public :: tab=achar(9), line_feed=achar(10)
   character(len=1), parameter, public :: comma=","
   integer, public :: i_str !Logical unit for writing strange characters not currently handled

   Type, public :: article
     character(len=50)  :: WOS=" "
     character(len=50)  :: ISBN=" "
     character(len=20)  :: Numb=" "
     character(len=256) :: Authors=" "
     character(len=256) :: Title=" "
     character(len=256) :: Journal=" "
     character(len=10)  :: Volume=" "
     character(len=20)  :: Pages=" "
     integer            :: year=0
     integer            :: citations=0
     character(len=6)   :: cyear=" "
     character(len=256) :: Book_Authors=" "
     character(len=256) :: Book_Title=" "
     character(len=256) :: Book_Pages=" "
     character(len=256) :: DOI=" "
     character(len=10),dimension(5)  :: instrument=" "
     integer            :: typ = 0 !=1 article =0 report
   End Type article

   Type(article), public, dimension(:), allocatable :: articles

   contains

   !  ASCII 1-127  ichar("a")=97 - ichar("z")=122    ichar("A")= 65 - ichar("Z")=90
   !
   !
   subroutine ISI_string(artic,ISI_str,code)
     type(article),     intent(in) :: artic
     character(len=*),  intent(out):: ISI_str
     character(len=*),  intent(in) :: code
     !--- Local variables ---!
     character(len=512) :: authors, mtitle
     character(len=20),dimension(20)  :: tit_words
     character(len=60),dimension(20)  :: author
     character(len=60)  :: autnam,author_1,author_2
     character(len=8)   :: initials
     integer :: i,j,k,l,nau, ncar,nc,nw,naut
     integer, dimension(60) :: pos
     integer, dimension(10) :: pb
     ISI_Str=" "

     !Title of the article
     mtitle=artic%title
     i=index(mtitle,'"')     !supress quotes
     if(i /= 0) then
       mtitle(i:i) = " "
       i=index(mtitle,'"')
       mtitle(i:i) = " "
       mtitle=adjustl(mtitle)
     end if

     if(index(code,"A") /= 0 .or. len_trim(code)== 0) then
        !Eliminate articles, parenthesis anb booleans from the title to simplify the search
        !call purge_bool_par_art(mtitle)
        authors=artic%authors
        i=index(authors,'"')     !supress quotes
        if(i /= 0) then
          authors(i:i) = " "
          i=index(authors,'"')
          authors(i:i) = " "
          authors=adjustl(authors)
        end if
        i=0
        do         !Suppress dots and compact initials of authors B.E.F. => BEF
          i=i+1
          if(authors(i:i) == ".") then
            authors=authors(1:i-1)//authors(i+1:)
            i=i-1
          end if
          if(i > len_trim(authors)) exit
        end do
        !Treat the name of the authors having blanks
        call get_separator_pos(trim(authors),comma,pos,ncar)
        nau=ncar+1
        if(nau > 1) then
          author(nau)=authors(pos(ncar)+1:)
          j=0
          do i=1,ncar
            author(i)= authors(j+1:pos(i)-1)
            j=pos(i)
          end do
        else
          author(1)=authors
        end if
        authors=" "
        naut=min(3,nau) ! Limit the output to 3 authors
        if(code == "AA" .or. code == "AY" .or. code == "YA" .or. artic%year==0) naut=nau
        do i=1,naut
            k=index(trim(author(i))," ",back=.true.)
            initials=author(i)(k+1:)
            autnam=adjustl(author(i)(1:k-1))
            call get_separator_pos(trim(autnam)," ",pb,nc)
            author_1=" "
            author_2=" "
            if (nc > 0) then
              author_1=trim(Pack_String(autnam))//" "//trim(initials)
              author_2=  trim(autnam(pb(nc)+1:))//" "//trim(initials)//autnam(1:1)
              author(i) = "("//trim(author_1)//" OR "//trim(author_2)//")"
            else
              call get_separator_pos(trim(autnam),"-",pb,nc)
              if (nc > 0) then
                author_1=trim(Pack_String(autnam))//" "//trim(initials)
                author_2=  autnam(1:pb(nc)-1)//trim(autnam(pb(nc)+1:))//" "//trim(initials)
                author(i) = "("//trim(author_1)//" OR "//trim(author_2)//")"
              end if
            end if
            if(len_trim(author(i)) > 0) then
              authors=trim(authors)//" "//trim(author(i))//" AND"
            end if
        end do
        i=Len_Trim(authors)-2
        if(authors(i:i+2) == "AND") authors=authors(1:i-1)
     end if

     Select Case(trim(code))
       case("AA")   !only authors
         write(unit=ISI_str,fmt="(a,i5,a)") "(AU=("//trim(authors)//"))"

       case("AY","YA")  !authors and year
         if(artic%year /= 0) then
           write(unit=ISI_str,fmt="(a,i5,a)") "(AU=("//trim(authors)//") AND PY=",artic%year,")"
         !else
         !  write(unit=ISI_str,fmt="(a)") "(AU=("//trim(authors)//"))"
         end if

       case("TA","AT")  !authors and title
         write(unit=ISI_str,fmt="(a,a)") "(AU=("//trim(authors)//")",' AND TI="'//trim(mtitle)//'")'

       case("TY","YT") !Title and year
         if(artic%year /= 0) then
           write(unit=ISI_str,fmt="(a,i5,a)") "(PY=",artic%year,' AND TI="'//trim(mtitle)//'")'
         !else
         !  write(unit=ISI_str,fmt="(a)") '(TI="'//trim(mtitle)//'")'
         end if

       case("IA","AI") !ISBN and autors
         write(unit=ISI_str,fmt="(a,a)") "(AU=("//trim(authors)//")",' AND IS="'//trim(artic%ISBN)//'")'

       case default  !authors, year and title
         if(artic%year /= 0) then
           write(unit=ISI_str,fmt="(a,i5,a)") "(AU=("//trim(authors)//") AND PY=",artic%year,' AND TI="'//trim(mtitle)//'")'
         else
           write(unit=ISI_str,fmt="(a,a)") "(AU=("//trim(authors)//")",' AND TI="'//trim(mtitle)//'")'
         end if
     End Select

     if(len_trim(ISI_str) > 1) call Replace_n_Search_nonascii(ISI_str)
     return
   end subroutine ISI_string

   Subroutine Replace_n_Search_nonascii(string)
     character(len=*), intent(in out) :: string
     !--- Local variables
     character(len=60)  :: warn
     character(len=30)  :: stchar
     integer :: i
     !Convert characters that need a special encoding (two or three bytes) to equivalent ascii
     !(single byte) character.
     call SString_Replace(string,"&lt;","<",warn)
     call SString_Replace(string,"&gt;",">",warn)
     call SString_Replace(string,"&#039;"," ",warn)
     call SString_Replace(string,"&#034;"," ",warn)
     call SString_Replace(string,"#"," ",warn)
     call SString_Replace(string,"ê","e",warn)
     call SString_Replace(string,"é","e",warn)
     call SString_Replace(string,"’"," ",warn)
     call SString_Replace(string,"á","a",warn)
     call SString_Replace(string,"í","i",warn)
     call SString_Replace(string,"ö","o",warn)
     call SString_Replace(string,"ó","o",warn)
     call SString_Replace(string,"ñ","n",warn)
     call SString_Replace(string,"ç","c",warn)
     call SString_Replace(string,"?","-",warn)
     call SString_Replace(string,"å","a",warn)
     call SString_Replace(string,"è","e",warn)
     call SString_Replace(string,"à","a",warn)
     call SString_Replace(string,"ï","i",warn)
     call SString_Replace(string,"ú","u",warn)
     call SString_Replace(string,"ü","u",warn)
     call SString_Replace(string,"ò","o",warn)
     call SString_Replace(string,"û","u",warn)

     !Now search for non-handled non-ascii characters
     do i=1,len_trim(string)
      if(iachar(string(i:i)) > 127) then
        stchar=" "
        if(iachar(string(i+1:i+1)) > 127) then
          stchar=string(i:i+1)
        else
          stchar=string(max(1,i-14):min(i+15,len_trim(string)))
        end if
        write(unit=i_str,fmt="(a,i5,a)") &
        " -> Strange character: "//string(i:i)//" <- Value of iachar: ",iachar(string(i:i)),"    "//trim(stchar)
      end if
     end do
     return
   End Subroutine Replace_n_Search_nonascii

   Subroutine purge_bool_par_art(string)
     character(len=*), intent(in out) :: string
     !--- Local variables
     character(len=60)  :: warn
     character(len=20), dimension(:), allocatable :: words
     integer :: i,j,nw
     integer,dimension(180) :: pos

     If(string(1:2) == "A ") string=adjustl(string(2:))
     !call SString_Replace(string,"(","",warn)
     !call SString_Replace(string,")","",warn)
     call SString_Replace(string," and "," ",warn)
     call SString_Replace(string," or "," ",warn)
     call SString_Replace(string,"The ","",warn)
     call SString_Replace(string," the "," ",warn)
     call SString_Replace(string," a "," ",warn)
     call SString_Replace(string," of "," ",warn)
     call SString_Replace(string,"="," ",warn)
     call SString_Replace(string," from "," ",warn)
     call SString_Replace(string,"."," ",warn)
     call SString_Replace(string,","," ",warn)
     call SString_Replace(string," in "," ",warn)
     call SString_Replace(string," near "," ",warn)
     string=adjustl(string)

     !---Attempt to eliminate formulae ?
     call get_separator_pos(trim(string)," ",pos,j)
     nw=j+1
     allocate(words(nw))
     words(1)=string(1:pos(1)-1)
     do i=2,nw-1
       words(i)=string(pos(i-1)+1:pos(i)-1)
     end do
     words(nw)=string(pos(nw-1)+1:)
     string=" "
     do i=1,nw
       if(Contains_a_number(words(i))) cycle
       if(len_trim(words(i)) < 4) cycle
       string=trim(string)//" "//trim(words(i))
     end do
     i=index(string,"(")
     j=index(string,")")
     if(j == 0 .and. i > 0) string(i:i) =""
     if(i == 0 .and. j > 0) string(j:j) =""
     return
   End Subroutine purge_bool_par_art

   Function Contains_a_number(string) result(ok)
     character(len=*), intent(in) :: string
     logical :: ok
     character(len=*),parameter,dimension(10) :: numbers=(/"0","1","2","3","4","5","6","7","8","9"/)
     integer :: i
     ok=.false.
     do i=1,10
       if(index(string,numbers(i)) /= 0) then
         ok=.true.
         exit
       end if
     end do
   End Function Contains_a_number


  End Module Data_Articles_Mod
