!!---- Program to generate ISI Web strings for searching articles starting with the
!!---- information obtained from the ILL database saved in "Tabulated" format.
!!---- Author: Juan Rodriguez-Carvajal (ILL)
  Program ILL_Pubs_WoS
   use CFML_String_Utilities, only: String_Count, SString_Replace,lcase
   use Data_Articles_Mod
   implicit none
   integer, parameter      :: max_chars=204800
   character(len=1024)     :: line, ISI_Str
   character(len=max_chars):: DOI_Str,Title_Str
   character(len=80)       :: fileinf,chain,name_jour

   integer :: ier,n,i,j, nart, ncar, npub, n_doi,nlog,n_title
   integer :: narg,iart=1
   logical :: esta, doi_only=.false., inc_auth=.false.
   character(len=12)   :: nam_inst=" "
   character(len=16)   :: third_arg=" "
   character(len=180)  :: file_inst
   integer             :: num_art



   write(unit=*,fmt="(/a)") "-------------------------------------------------------------------------------------"
   write(unit=*,fmt="( a)") " Program to generate a long string for Advanced Search in the Web of Science "
   write(unit=*,fmt="( a)") "     It needs a text file from the ILL database in Tabulated format."
   write(unit=*,fmt="( a)") "                       May 2014, JRC -ILL"
   write(unit=*,fmt="( a)") " Usage: -> ILL_Pubs_WoS  my_input_file  my_output_file [doi_only or include_authors] "
   write(unit=*,fmt="(a/)") "-------------------------------------------------------------------------------------"
   narg= COMMAND_ARGUMENT_COUNT()

   if( narg <= 1) then

       write(unit=*,fmt="(a)")    " => The program 'ILL_Pubs_WoS' should be invoked at least with two arguments as indicated in the banner "
       write(unit=*,fmt="(a)")    "    The two arguments are: "
       write(unit=*,fmt="(a)")    "     1: the name of the saved file from the ILL database "
       write(unit=*,fmt="(a)")    "     2: the name of the file in which the search string will be written "
       write(unit=*,fmt="(a)")    "    "
       write(unit=*,fmt="(a)")    " => An additional optional arguments is the type of output: "
       write(unit=*,fmt="(a)")    "    The values of the additional argument are: doi_only or include_authors "
       stop

   else

     call GET_COMMAND_ARGUMENT(1, fileinf)
     call GET_COMMAND_ARGUMENT(2, file_inst)
     if(narg > 2) then
       call GET_COMMAND_ARGUMENT(3, third_arg)
       call lcase(third_arg)
       if(trim(third_arg) == "doi_only")        doi_only=.true.
       if(trim(third_arg) == "include_authors") inc_auth=.true.
     end if

   end if
   inquire(file=fileinf,exist=esta)
   if(.not. esta) then
     write(unit=*,fmt="(a)")    " => The input file: "//trim(fileinf)//" doesn't exist!"
     stop
   end if

   num_art=0
   open(unit=iart,file=fileinf,status="old",action="read", position="rewind")
   nart=0
   !Determine the number of papers in the file
   do
     read(unit=iart,fmt="(a)",iostat=ier)  line
     if(ier /= 0) exit
     if(len_trim(line) == 0) cycle
     if(line(1:6) =="Number") nart=nart+1
   end do
   rewind(unit=iart)
   allocate(articles(nart))

   n=0
   n_doi=0; n_title=0

   do
     read(unit=iart,fmt="(a)",iostat=ier)  line
     if(ier /= 0) then
       write(unit=*,fmt="(a,i7,a)") " => Number of Articles read: ",n," -> Last article: "//trim(articles(n)%Numb)
       exit
     end if
     if(len_trim(line) == 0) cycle
     j=index(line,tab)
     Select Case(line(1:j-1))
       Case("Number")
         n=n+1
         articles(n)%Numb= line(j+1:)       !<  1  Numero
         articles(n)%instrument(1)= trim(Nam_inst)
       Case("Author")
         articles(n)%Authors= line(j+1:)    !<  2  Auteurs
       Case("Title")
         articles(n)%Title= line(j+1:)      !<  3  Titre
       Case("Journal title")
         articles(n)%Journal= line(j+1:)    !<  4  Journal
       Case("Volume")
         articles(n)%Volume= line(j+1:)     !<  5  Volume
       Case("Pages")
         articles(n)%Pages= line(j+1:)      !<  6  Pages
       Case("DOI")
         i=index(line,"Keyword")
         if( i /= 0) then
            articles(n)%DOI= line(j+1:i-1)        !<  6  DOI
         else
            articles(n)%DOI= line(j+1:)        !<  6  DOI
         end if
       Case("Year")
         read (unit=line(j+1:),fmt=*,iostat=ier) articles(n)%year  !< 7 Annee
         if(ier /= 0) then
           write(unit=*,fmt="(a,a)")  " => Error reading the year in the article: "//trim(articles(n)%Numb), &
           " ... The document "//trim(articles(n)%Numb)//" is discarded"
           n=n-1
           cycle
         end if
       Case Default
         cycle
     End Select
     articles(n)%typ=1
   end do
   close(unit=iart)
   nart=n
   open(newunit=i_str,file="strange_chars.inf",status="replace",action="write")

   DOI_Str=" "; Title_Str=" "
   npub=0
   do j=1,nart
      if(articles(j)%year == 0) then
           write(unit=*,fmt="(a,a)")  " => No year available in the article: "//trim(articles(j)%Numb), &
           " ... The document "//trim(articles(j)%Numb)//" is discarded"
           cycle
      end if
      if(len_trim(articles(j)%DOI) == 0) then
        if(doi_only) cycle
        n_title=n_title+1
        if(inc_auth) then
           call ISI_string(articles(j),ISI_str,"include_authors")
        else
           call ISI_string(articles(j),ISI_str)
        end if
        if(n_title == 1) then
          Title_Str=trim(ISI_str)
        else
          nlog=len_trim(ISI_str)
          ncar=len_trim(Title_Str)
          !write(unit=*,fmt="(3i9)")  nlog, ncar,nlog+ncar
          if(nlog+ncar <= max_chars) then
             Title_Str=trim(Title_Str)//" OR"//line_feed//trim(ISI_str)
          else
             exit
          end if
        end if
      else
        n_doi=n_doi+1
        ISI_str="DO=("//trim(articles(j)%DOI)//")"
        if(n_doi == 1) then
          DOI_Str=trim(ISI_str)
        else
          nlog=len_trim(ISI_str)
          ncar=len_trim(DOI_Str)
          if(nlog+ncar <= max_chars) then
             DOI_Str=trim(DOI_Str)//" OR"//line_feed//trim(ISI_str)
          else
             exit
          end if
        end if
      end if
      npub=npub+1
   end do
   open(unit=iart,file=trim(file_inst),status="replace",action="write")
     if(doi_only) then
        write(unit=iart,fmt="(a)")  trim(DOI_str)
     else
        write(unit=iart,fmt="(a)")  trim(DOI_str)//" OR"//line_feed//trim(Title_Str)
     end if
   close(unit=iart)

   write(unit=*,fmt="(a,i6)") " => Number of articles saved   : ",npub
   write(unit=*,fmt="(a,i6,tr2,f6.2,a)") " => Number of articles with DOI: ",n_doi,100.0*real(n_doi)/real(nart),"%"
   write(unit=*,fmt="(a)") " => String to paste in Advanced Search of WoS in file: "//trim(file_inst)

   stop
  End Program ILL_Pubs_WoS