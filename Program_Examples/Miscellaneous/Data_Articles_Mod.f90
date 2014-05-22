!!---- Module containing types and procedures to generate ISI Web strings for searching articles
!!---- Author: Juan Rodriguez-Carvajal (ILL)

  Module Data_Articles_Mod
   Use CFML_String_Utilities, only: Pack_String, get_separator_pos, SString_Replace
   implicit none
   private
   public :: ISI_string, ISI_string_title
   character(len=1), parameter, public :: tab=achar(9), line_feed=achar(10)
   character(len=1), parameter, public :: comma=","

   Type, public :: article
     character(len=20)  :: Numb=" "
     character(len=256) :: Authors=" "
     character(len=256) :: Title=" "
     character(len=256) :: Journal=" "
     character(len=10)  :: Volume=" "
     character(len=20)  :: Pages=" "
     integer            :: year=0
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
   subroutine ISI_string(artic,ISI_str)
     type(article),    intent(in) :: artic
     character(len=*), intent(out):: ISI_str
     !--- Local variables ---!
     character(len=512) :: authors, mtitle
     character(len=20),dimension(20)  :: tit_words
     character(len=60),dimension(20)  :: author
     character(len=60)  :: autnam,author_1,author_2
     character(len=8)   :: initials
     integer :: i,j,k,l,nau, ncar,nc,nw               ! 1   2   3   4   5   6   7   8   9  10  11  12  13
     character(len=1),dimension(13):: list_non_ascii=(/"å","ä","á","é","è","í","ï","ö","ó","ü","ù","ú","ñ"/)
     character(len=1),dimension(13):: list_____ascii=(/"a","a","a","e","e","i","i","o","o","u","u","u","n"/)
     integer, dimension(60) :: pos
     integer, dimension(10) :: pb

     ISI_Str=" "
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
     do i=1,nau
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
         end if
         if(len_trim(author(i)) > 0) then
           authors=trim(authors)//" "//trim(author(i))//" AND"
         end if
     end do
     i=Len_Trim(authors)-2
     if(authors(i:i+2) == "AND") authors=authors(1:i-1)

     !Title of the article
     mtitle=artic%title
     i=index(mtitle,'"')     !supress quotes
     if(i /= 0) then
       mtitle(i:i) = " "
       i=index(mtitle,'"')
       mtitle(i:i) = " "
       mtitle=adjustl(mtitle)
     end if
     call get_separator_pos(trim(mtitle)," ",pos,ncar)
     nw=ncar+1
     tit_words(nw)=mtitle(pos(ncar)+1:)
     j=0
     do i=1,ncar
       tit_words(i)= mtitle(j+1:pos(i)-1)
       j=pos(i)
     end do
     mtitle=" "
     do i=1,nw
         if(len_trim(tit_words(i)) > 4) then
           mtitle=trim(mtitle)//" "//trim(tit_words(i))//" AND"
         end if
     end do
     i=Len_Trim(mtitle)-2
     if(mtitle(i:i+2) == "AND") mtitle=mtitle(1:i-1)

     ISI_str=" "
     if( nau > 2) then
       write(unit=ISI_str,fmt="(a,i5,a)") "(AU=("//trim(authors)//") AND PY=",artic%year,")"
     else
       write(unit=ISI_str,fmt="(a,i5,a)") "(AU=("//trim(authors)//") AND PY=",artic%year," AND TI=("//trim(mtitle)//"))"
     end if

     !Now replace non-ascii characters by the equivalent ascii ones
     do i=1,len_trim(ISI_str)
      do j=1,13
        if(List_non_ascii(j) == ISI_str(i:i)) then
          ISI_str(i:i)=list_____ascii(j)
          exit
        end if
      end do
     end do
   end subroutine ISI_string

   subroutine ISI_string_title(artic,ISI_str)
     type(article),    intent(in) :: artic
     character(len=*), intent(out):: ISI_str
     !--- Local variables ---!
     character(len=512) :: authors, mtitle
     character(len=20),dimension(20)  :: tit_words
     character(len=60),dimension(20)  :: author
     character(len=60)  :: autnam,author_1,author_2,warn
     character(len=8)   :: initials
     integer :: i,j,l,nau, ncar,nc,nw               ! 1   2   3   4   5   6   7   8   9  10  11  12  13
     character(len=1),dimension(13):: list_non_ascii=(/"å","ä","á","é","è","í","ï","ö","ó","ü","ù","ú","ñ"/)
     character(len=1),dimension(13):: list_____ascii=(/"a","a","a","e","e","i","i","o","o","u","u","u","n"/)
     integer, dimension(60) :: pos
     integer, dimension(10) :: pb

     authors=artic%authors
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
     ISI_str=" "
     if( artic%year == 0) then
       ISI_str='(TI="'//trim(mtitle)//'")'
     else
       write(unit=ISI_str,fmt="(a,i5,a)") "(PY=",artic%year,' AND TI="'//trim(mtitle)//'")'
     end if
     ISI_str=trim(ISI_str)
     !Now replace non-ascii characters by the equivalent ascii ones
     do i=1,len_trim(ISI_str)
      do j=1,13
        if(List_non_ascii(j) == ISI_str(i:i)) then
          ISI_str(i:i)=list_____ascii(j)
          exit
        end if
      end do
      if(iachar(ISI_str(i:i)) > 127) then
        !write(unit=*,fmt="(a,i5)") " -> Strange character: "//ISI_str(i:i)//" <- Value of iachar: ",iachar(ISI_str(i:i))
        ISI_str(i:i)=" "
      end if
     end do
     call SString_Replace(ISI_str,"&lt;"," ",warn)
     call SString_Replace(ISI_str,"&gt;"," ",warn)
     call SString_Replace(ISI_str,"&#039;"," ",warn)
     call SString_Replace(ISI_str,"&#034;"," ",warn)
     call SString_Replace(ISI_str,"#"," ",warn)
     return
   End subroutine ISI_string_title

  End Module Data_Articles_Mod