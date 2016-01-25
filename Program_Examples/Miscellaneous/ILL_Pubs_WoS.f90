!!---- Program to generate ISI Web strings for searching articles starting with the
!!---- information obtained from the ILL database saved in "Tabulated" format.
!!---- Author: Juan Rodriguez-Carvajal (ILL)
  Program ILL_Pubs_WoS
   use CFML_String_Utilities, only: lcase
   use Data_Articles_Mod
   implicit none
   integer, parameter      :: max_chars=204800
   character(len=1024)     :: line, ISI_Str
   character(len=max_chars):: DOI_Str,Title_Str,WOS_Str,ISBN_Str
   character(len=80)       :: fileinf,chain,name_jour

   integer :: ier,n,i,j, nart, ncar, npub, n_doi=0,nlog,n_title, n_wos=0,n_isbn=0
   integer :: narg,iart=1
   logical :: esta, doi_only=.false., inc_code=.false., non_doi=.false.
   character(len=2)    :: CODE=" "
   character(len=12)   :: nam_inst=" "
   character(len=16)   :: third_arg=" "
   character(len=180)  :: file_inst
   integer             :: num_art


   open(newunit=i_str,file="ILL_Pubs_WoS.log",status="replace",action="write")

   write(unit=*,fmt="(/a)") "-------------------------------------------------------------------------------------"
   write(unit=*,fmt="( a)") " Program to generate a long string for Advanced Search in the Web of Science "
   write(unit=*,fmt="( a)") "     It needs a text file from the ILL database in Tabulated format."
   write(unit=*,fmt="( a)") "                       January 2016, JRC -ILL"
   write(unit=*,fmt="( a)") " Usage: -> ILL_Pubs_WoS input_file moutput_file [doi_only, non_doi or CODE] "
   write(unit=*,fmt="(a/)") "-------------------------------------------------------------------------------------"

   write(unit=i_str,fmt="(/a)") "-------------------------------------------------------------------------------------"
   write(unit=i_str,fmt="( a)") " Program to generate a long string for Advanced Search in the Web of Science "
   write(unit=i_str,fmt="( a)") "     It needs a text file from the ILL database in Tabulated format."
   write(unit=i_str,fmt="( a)") "                   January 2016, JRC -ILL"
   write(unit=i_str,fmt="( a)") " Usage: -> ILL_Pubs_WoS input_file moutput_file [doi_only, non_doi or CODE] "
   write(unit=i_str,fmt="(a/)") "-------------------------------------------------------------------------------------"

   narg= COMMAND_ARGUMENT_COUNT()

   if( narg <= 1) then

     write(unit=*,fmt="(a)") " => The program 'ILL_Pubs_WoS' should be invoked at least with two arguments as indicated above"
     write(unit=*,fmt="(a)") "    The two arguments are: "
     write(unit=*,fmt="(a)") "     1: the name of the saved file from the ILL database "
     write(unit=*,fmt="(a)") "     2: the name of the file in which the search string will be written "
     write(unit=*,fmt="(a)") "    "
     write(unit=*,fmt="(a)") " => An additional optional argument concerns the type of output: "
     write(unit=*,fmt="(a)") "    The values of the additional argument are: doi_only, non_doi or CODE "
     write(unit=*,fmt="(a)") "    In which CODE may be: AA (only authors), AY (authors + year), "
     write(unit=*,fmt="(a)") "                          TY (title+year),   TA (authors + title)  "
     write(unit=*,fmt="(a)") "    The codes are commutative, e.g. AY is the same as YA  "
     write(unit=*,fmt="(a)") "    If no code is provided it is equivalent to: authors + year + title"
     write(unit=*,fmt="(a)") "    The '+' sign corresponds to a logical 'AND'"
     stop

   else

     call GET_COMMAND_ARGUMENT(1, fileinf)
     call GET_COMMAND_ARGUMENT(2, file_inst)
     if(narg > 2) then
       call GET_COMMAND_ARGUMENT(3, third_arg)
       CODE=trim(third_arg)
       call lcase(third_arg)
       if(trim(third_arg) == "doi_only")  then
         doi_only=.true.
         code=" "
       else if(trim(third_arg) == "non_doi") then
         non_doi =.true.
         code=" "
       else
         inc_code=.true.
         write(unit=*,fmt="(a)") " => CODE: "//CODE
       end if
     else
       code=" "      
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
   n_doi=0; n_title=0; n_isbn=0; n_wos=0

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
       Case("ISBN")
         articles(n)%ISBN= line(j+1:)       !<  5  Journal
       Case("Volume")
         articles(n)%Volume= line(j+1:)     !<  6  Volume
       Case("Pages")
         articles(n)%Pages= line(j+1:)      !<  7  Pages
       Case("WoS number")
            articles(n)%WOS= adjustl(line(j+1:))  !< 8 WOS
       Case("DOI","DOI added")
         i=index(line,"Keyword")
         if( i /= 0) then
            articles(n)%DOI= line(j+1:i-1)     !<  9  DOI
         else
            articles(n)%DOI= line(j+1:)        !<  9  DOI
         end if
       Case("Year")
         read (unit=line(j+1:),fmt=*,iostat=ier) articles(n)%year  !< 10 Annee
         if(ier /= 0) then
           if(len_trim(articles(n)%DOI) == 0 .and. len_trim(articles(n)%WOS) == 0) then
             write(unit=*,fmt="(a,a)")  " => Error reading the year in the article: "//trim(articles(n)%Numb), &
             " ... The document "//trim(articles(n)%Numb)//" is discarded"
             write(unit=i_str,fmt="(a,a)")  " => Error reading the year in the article: "//trim(articles(n)%Numb), &
             " ... The document "//trim(articles(n)%Numb)//" is discarded"
             n=n-1
             cycle
           else
             articles(n)%year=0
           end if
         end if
       Case Default
         cycle
     End Select
     articles(n)%typ=1
   end do
   close(unit=iart)
   nart=n

   DOI_Str=" "; Title_Str=" ";WOS_Str=" "; ISBN_Str=" "

   do j=1,nart
      !if(articles(j)%year == 0 .and. len_trim(articles(j)%DOI) == 0  &
      !   .and. len_trim(articles(j)%WOS) == 0  .and. len_trim(articles(j)%ISBN) == 0) then
      !     write(unit=*,fmt="(a,a)")  " => No year/DOI/WOS/ISBN available in the article: "//trim(articles(j)%Numb), &
      !     " ... The document "//trim(articles(j)%Numb)//" is discarded"
      !     write(unit=i_str,fmt="(a,a)")  " => No year/DOI/WOS/ISBN available in the article: "//trim(articles(j)%Numb), &
      !     " ... The document "//trim(articles(j)%Numb)//" is discarded"
      !     cycle
      !end if

      !Article with WOS
      if(len_trim(articles(j)%WOS) /= 0) then   !WOS accession number is provided (prioritary)
        if(non_doi) cycle
        n_wos=n_wos+1
        ISI_str="UT=("//trim(articles(j)%WOS)//")"
        if(n_wos == 1) then
          WOS_Str=trim(ISI_str)
        else
          nlog=len_trim(ISI_str)
          ncar=len_trim(WOS_Str)
          if(nlog+ncar <= max_chars) then
             WOS_Str=trim(WOS_Str)//" OR"//line_feed//trim(ISI_str)
          else
             exit
          end if
        end if

      !Article with DOI
      else if (len_trim(articles(j)%DOI) /= 0) then  !DOI is provided
        if(non_doi) cycle
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

      !Article with ISBN
      else if (len_trim(articles(j)%ISBN) /= 0) then  !ISBN is provided
        if(non_doi .or. doi_only) cycle
        call ISI_string(articles(j),ISI_str,"IA")
        n_isbn=n_isbn+1
        if(n_isbn == 1) then
          ISBN_Str=trim(ISI_str)
        else
          nlog=len_trim(ISI_str)
          ncar=len_trim(ISBN_Str)
          if(nlog+ncar <= max_chars) then
             ISBN_Str=trim(ISBN_Str)//" OR"//line_feed//trim(ISI_str)
          else
             exit
          end if
        end if
      !Article without DOI,WOS and ISBN: need to construct a string with authors, title, year
      !This part should be improved
      !if(len_trim(articles(j)%DOI) == 0 .and. len_trim(articles(j)%WOS) == 0 .and. articles(j)%year /= 0) then
      else
        if(doi_only) cycle
        n_title=n_title+1
        if(inc_code .or. non_doi) then
          call ISI_string(articles(j),ISI_str,CODE)
        else
           call ISI_string(articles(j),ISI_str," ")
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
      end if
      
   end do
   
   npub=n_doi+n_wos+n_isbn+n_title
   
   open(unit=iart,file=trim(file_inst),status="replace",action="write")
     if(doi_only) then
        write(unit=iart,fmt="(a)")  trim(WOS_str)//" OR"//line_feed//trim(DOI_str)//" OR"//line_feed//trim(ISBN_str)
     else if(non_doi) then
        write(unit=iart,fmt="(a)")  trim(Title_Str)
     else
        write(unit=iart,fmt="(a)")  &
        trim(WOS_str)//" OR"//line_feed//trim(DOI_str)//" OR"//line_feed//trim(ISBN_str)//" OR"//line_feed//trim(Title_Str)
     end if
   close(unit=iart)

   write(unit=*,fmt="(/a,i6)")                 " => Number of articles treated             : ",nart
   write(unit=*,fmt="(a,i6)")                  " => Number of articles saved               : ",npub
   write(unit=i_str,fmt="(/a,i6)")             " => Number of articles treated             : ",nart
   write(unit=i_str,fmt="(a,i6)")              " => Number of articles saved               : ",npub

   if(non_doi) then
     write(unit=*,fmt="(a)")                   " => The saved articles have no DOI in the ILL database "
     write(unit=i_str,fmt="(a)")               " => The saved articles have no DOI in the ILL database "
   else
     write(unit=*,fmt="(a,i6,tr2,f6.2,a)")     " => Number of articles with DOI            : ",n_doi,100.0*real(n_doi)/real(nart),"%"
     write(unit=i_str,fmt="(a,i6,tr2,f6.2,a)") " => Number of articles with DOI            : ",n_doi,100.0*real(n_doi)/real(nart),"%"
     write(unit=*,fmt="(a,i6,tr2,f6.2,a)")     " => Number of articles with WOS            : ",n_wos,100.0*real(n_wos)/real(nart),"%"
     write(unit=i_str,fmt="(a,i6,tr2,f6.2,a)") " => Number of articles with WOS            : ",n_wos,100.0*real(n_wos)/real(nart),"%"
   end if
   if(n_isbn > 0) then
     write(unit=*,fmt="(a,i6,tr2,f6.2,a)")     " => Number of articles with ISBN           : ",n_isbn,100.0*real(n_isbn)/real(nart),"%"
     write(unit=i_str,fmt="(a,i6,tr2,f6.2,a)") " => Number of articles with ISBN           : ",n_isbn,100.0*real(n_isbn)/real(nart),"%"
   end if
   if(n_title > 0) then
     write(unit=*,fmt="(a,i6,tr2,f6.2,a)")     " => Number of articles without DOI/WOS/ISBN: ",n_title,100.0*real(n_title)/real(nart),"%"
     write(unit=i_str,fmt="(a,i6,tr2,f6.2,a)") " => Number of articles without DOI/WOS/ISBN: ",n_title,100.0*real(n_title)/real(nart),"%"
   end if
   write(unit=*,fmt="(a)") " => String to paste in Advanced Search of WoS in file: "//trim(file_inst)
   write(unit=*,fmt="(a)") " => Log file: ILL_Pubs_WoS.log "
   write(unit=i_str,fmt="(a)") " => String to paste in Advanced Search of WoS in file: "//trim(file_inst)
   close(unit=i_str)

   stop
  End Program ILL_Pubs_WoS
