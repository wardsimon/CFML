PROGRAM DIFFaX_to_FAULTS

    implicit none                             ! turn off the implicit type of the variables

    !Declare Global Constants
      integer, parameter :: ip = 5,&          ! ip    -  standard input device unit number
                            op = 6,&          ! op    -  standard output device unit number
                            f_in = 10,&       ! f_in  -  input file unit number
                            f_out = 20        ! f_out -  output file unit number
      integer, parameter :: maxcharD = 200    ! maximum number of character per line in the DIFFaX input file
      integer, parameter :: maxcharF = 132    ! maximum number of character per line in the FAULTS input file
      character(len=10), parameter :: ext = ".flts" ! extension of the FAULTS input file created by the program

    !Declare Global Variables
      character(len=1)                                 :: keyw
      character(len=20)                                :: infile          ! name of the DIFFaX input file provided by the user
      character(len=20)                                :: shortname       ! name of the input file without extension
      character(len=30)                                :: filename        ! name of the FAULTS input file created by the program
      character(len=maxcharD),dimension(:),allocatable :: dfile           ! list of lines of the DIFFaX input file
      character(len=maxcharD),dimension(:),allocatable :: tfile           ! temporary list of lines of the DIFFaX in put file, with splitted lines of comment
      integer,dimension(:),allocatable                 :: tlineD          ! remember the line number in original DIFFaX input file for each element in tfile
      character(len=maxcharF),dimension(:),allocatable :: ffile           ! list of lines of the FAULTS input file
      integer                                          :: maxlinesD       ! maximum number of lines of the DIFFaX input file
      integer                                          :: maxlinesF       ! maximum number of lines of the FAULTS input file
      integer                                          :: nblinesD        ! number of lines in DIFFaX input file
      integer                                          :: nblinesT        ! number of lines in tfile
      integer                                          :: nlfinalF        ! number of lines in ffile
      integer                                          :: S1,S2,S3,S4,S5  ! line nb where each section begins :
                                                                    ! 1-INSTRUMENTAL 2-STRUCTURAL 3-LAYER 4-STACKING 5-TRANSITIONS
      integer                                          :: nlayers         ! number of layers
      integer                                          :: narg            ! Number of command line arguments

    !BEGINNING OF THE PROGRAM
      !---- Arguments on the command line ----!
      narg=command_argument_count()

      infile=" "
      if (narg > 0) then
         call get_command_argument(1,infile)
      end if

      CALL Salute()
      CALL Read_diffax(infile)               ! Copy DIFFaX input file into array dfile
      CALL Split_comments()            ! Split lines which begin with a comment,
                                       !               fait commencer toutes lignes de commentaires par un point d'exclamation,
                                       !               cadre toutes les lignes à gauche
      CALL Index_sections()            ! Index the different sections S1,S2,S3,S4,S5
      CALL Extract_data()              ! analyse the DIFFaX information
      CALL Write_flts()                ! create the FAULTS input file

      write(unit=*,fmt="(/,a)") " => Normal end, please press <Enter> to finish ..."
      read(unit=*,fmt="(a)") keyw

    !END OF THE PROGRAM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains

      SUBROUTINE Salute()

        write(op,"(a)") "                                                                             "
        write(op,"(a)") "_____________________________________________________________________________"
        write(op,"(a)") "                                                                             "
        write(op,"(a)") "                           ____ DIFFaX2FAULTS ____                           "
        write(op,"(a)") "                                                                             "
        write(op,"(a)") "_____________________________________________________________________________"
        write(op,"(a)") "_____________________________________________________________________________"
        write(op,"(a)") "                                                                             "
        write(op,"(a)") "         A tool to convert DIFFaX input files to FAULTS input files          "
        write(op,"(a)") "                                                                             "
        write(op,"(a)") "                                                                             "
        write(op,"(a)") "            Authors of DIFFaX:                                               "
        write(op,"(a)") "              M. M. J. Treacy , M. W. Deem & J. M. Newsam                    "
        write(op,"(a)") "                                                                             "
        write(op,"(a)") "            Authors of FAULTS:                                               "
        write(op,"(a)") "             Montse Casas-Cabanas     (CIC Energigune)                       "
        write(op,"(a)") "             Jokin Rikarte            (CIC Energigune)                       "
        write(op,"(a)") "             Marine Reynaud           (CIC Energigune)                       "
        write(op,"(a)") "             Juan Rodriguez-Carvajal  (Institut Laue-Langevin)               "
        write(op,"(a)") "                                                                             "
        write(op,"(a)") "                                                                             "
        write(op,"(a)") "                             [version: Jan. 2015]                            "
        write(op,"(a)") "_____________________________________________________________________________"
        write(op,"(a)") "                                                                             "

      END SUBROUTINE Salute

      SUBROUTINE Read_diffax(filen)
        character(len=*), intent(in) :: filen
        ! Local Variables
          integer             :: ier         ! controls : 0 = successful, >0 error
          integer             :: i,k         ! counters
          character(LEN=maxcharD)  :: line
          character(LEN=100)  :: linenb
          logical             :: error

        if(len_trim(filen) == 0) then
          write(UNIT=op,FMT="(a)",advance="no") " => Enter the full name of the DIFFaX input file (with extension): "
          read (UNIT=ip,FMT="(a)") infile
          if(len_trim(infile) == 0) then
               write(unit=*,fmt="(/,a)") " => No file provided, please press <Enter> to finish ..."
               read(unit=*,fmt="(a)") keyw
               stop
          end if
        else
          infile=filen
        end if

        OPEN (unit=f_in, file=trim(infile), status='old', action='read', position='rewind', iostat=ier)   ! Old = the DIFFaX input file should exist
                                                                                           ! ier = control :  0 = successful, >0 error
            if (ier /= 0) then
                write(UNIT=op,FMT="(2a)") " => Error opening the DIFFaX input file ", trim(infile)
            else
                write(UNIT=op,FMT="(3a)") " "
                write(UNIT=op,FMT="(3a)") " => The DIFFaX input file ", trim(infile), " was successfully opened"
            end if

            ! Count lines in infile


            nblinesD=0
            error=.false.
            linenb="  "

            do
                read(unit=f_in, fmt="(a)", iostat=ier) line
                nblinesD=nblinesD+1
                if (ier /= 0) then
                     exit
                end if
            end do

            write(unit=op,fmt="(a,i3,a)") " => The DIFFaX input file contains ", nblinesD, " lines"

            maxlinesD=nblinesD
            maxlinesF=2*maxcharD

            ! Allocate arrays

            if (allocated (dfile)) deallocate(dfile)                ! list of lines of the DIFFaX input file
            allocate(dfile(maxlinesD))

            if (allocated (tfile)) deallocate(tfile)                ! temporary list of lines of the DIFFaX in put file, with splitted lines of comment
            allocate(tfile(maxlinesF))

            if (allocated (tlineD)) deallocate(tlineD)              ! remember the line number in original DIFFaX input file for each element in tfile
            allocate(tlineD(maxlinesF))

            if (allocated (ffile)) deallocate(ffile)                ! list of lines of the FAULTS input file
            allocate(ffile(maxlinesF))

      !character(maxcharD),dimension(maxlinesD)   :: dfile          ! list of lines of the DIFFaX input file
      !character(maxcharD),dimension(maxlinesF)   :: tfile          ! temporary list of lines of the DIFFaX in put file, with splitted lines of comment
      !integer,dimension(maxlinesF)               :: tlineD         ! remember the line number in original DIFFaX input file for each element in tfile
      !character(maxcharF),dimension(maxlinesF)   :: ffile          ! list of lines of the FAULTS input file



            ! Read infile

            rewind(unit=f_in)                        ! go back to the beginning of the file infile
            dfile="  "

            do i=1,maxlinesD
                read(unit=f_in, fmt="(a)", iostat=ier) line
                !write(unit=op, fmt="(a)") dfile(i)
                k=len_trim(adjustl(line))
                if(k > maxcharF) then
                    write(UNIT=op, FMT="(a)") " => Caution: Some information may be lost during the conversion: "
                    write(UNIT=op, FMT="(a,i3,a,i3,a)") "    The length of line ", i, &
                                                               " is greater than ", maxcharF, " characters:"
                    write(UNIT=op, FMT="(a,a50,a)") "    ``", line, "...''"
                end if
                dfile(i)=adjustl(line(1:maxcharF))
                if (ier /= 0) then
                     exit
                     !error=.true.
                     !linenb = trim(linenb) // ", " // trim(adjustl(xstring))
                end if
            end do

           ! if (error) then
           !     write(UNIT=op,FMT="(2a)") " => Errors reading the DIFFaX input file ", infile
           !     write(UNIT=op,FMT="(3a)") "    Line(s)", trim(linenb), " was(were) not read correctly"
           ! else
           !     write(UNIT=op,FMT="(3a)") " => The DIFFaX input file ", trim(infile), " was successfully read"
           !     write(unit=op,fmt="(a,i3)") " => The DIFFaX input file contains ", nblinesD, " lines"
           ! end if


        CLOSE (f_in)

        !create the file name of the FAULTS input file
         k=index(infile,".")
         shortname=infile(1:(k-1))
         filename = trim(shortname)//ext


      END SUBROUTINE Read_diffax

      SUBROUTINE Split_comments()

        ! Local Variables
          integer             :: i,j,k,c      ! counters
          character(maxcharD) :: texto,com,ligne

        tfile = "  "
        tlineD = 0
        j=0

        do i=1,nblinesD
           j=j+1
           texto=adjustl(dfile(i))                    ! cadrage a gauche de la ligne i
           do c=1,10    ! max 10 comments
                if(texto(1:1) == "{") then
                    k=index(texto,"}")
                    if(k == 0) then
                        write(UNIT=op,FMT="(a,i3,a)") " => Caution: Closing brace is missing on line ", i, &
                                                      ". This line is treated as a comment."
                        ligne="! " // texto
                        tfile(j) = ligne(1:maxcharF)
                        tlineD(j)= i
                        exit
                    else
                        com= "! " // texto(1:k)
                        texto=adjustl(texto(k+1:))
                        if(texto(1:1)==" ") then
                            ligne=trim(com)
                            tfile(j)=ligne(1:maxcharF)
                            tlineD(j)= i
                            exit
                        else if(texto(1:1)=="!") then
                            ligne= trim(com) // " " // trim(texto(2:))
                            tfile(j)=ligne(1:maxcharF)
                            tlineD(j)= i
                            exit
                        else if(texto(1:1)=="{") then
                            ligne=trim(com)
                            tfile(j)=ligne(1:maxcharF)
                            tlineD(j)= i
                            j=j+1
                            cycle
                        else
                            ligne=trim(com)
                            tfile(j)=ligne(1:maxcharF)
                            tlineD(j)= i
                            j=j+1
                            tfile(j)=texto(1:maxcharD)
                            tlineD(j)= i
                            exit
                        end if
                    end if
                else
                    tfile(j)=texto(1:maxcharD)
                    tlineD(j)= i
                end if
            end do
        end do

        nblinesT=j

      END SUBROUTINE Split_comments

      SUBROUTINE Index_sections()

        ! Local Variables
          integer             :: i,k          ! counters
          character(maxcharD) :: texto,key

        ! Initialize variables
        k=0

        do i=1,nblinesT
            texto=adjustl(tfile(i))                    ! cadrage a gauche de la ligne i
            if (texto(1:1)==" " .or. texto(1:1)=="!") then
                cycle
            else
                k=index(texto," ")
                key=U_Case(texto(1:k-1))
                if(key=="INSTRUMENTAL") then
                    S1=i
                else if(key=="STRUCTURAL") then
                    S2=i
                else if(key=="LAYER") then
                    texto=adjustl(texto(k:))
                    k=index(texto," ")
                    key=texto(1:k-1)
                    if(key=="1") then
                        S3=i
                    else
                        cycle
                    end if
                else if(key=="STACKING") then
                    S4=i
                else if(key=="TRANSITIONS") then
                    S5=i
                else
                    cycle
                end if
            end if
        end do

      END SUBROUTINE Index_sections

      SUBROUTINE Extract_data()

        ! Local Variables
          !integer            :: ier         ! control : 0 = successful, >0 error
          integer             :: i,j !,k       ! counters
          character(maxcharD) :: ligne

        ! Initialize variables
          i=0
          j=3
          ffile=" "


        ! TITLE section
         ffile(1) = "  TITLE"
         ffile(2) = "  " // shortname
         ffile(3) = "  "

         if(S1>1) then
            do i=1,S1-1
                j=j+1
                ligne=tfile(i)
                ffile(j)=ligne(1:maxcharF)
            end do
         end if

        ! INSTRUMENTAL section
            CALL instrumental(j)

        ! STRUCTURAL section
            CALL structural(j)

        ! LAYER section
            CALL layers(j)

        ! STACKING section
            CALL stacking(j)

        ! TRANSITIONS section
            CALL transitions(j)

        ! CALCULATION section
            CALL calculation(j)

        nlfinalF=j                              ! Final number of lines in the FAULTS input file

        if (nlfinalF>maxlinesF) then
            write(unit=op,fmt="(a)") " => DIFFaX2FAULTS has reached the maximum number of lines permitted for" &
                                        // "the FAULTS input file. Please check that no information was lost."
        end if

      END SUBROUTINE Extract_data

      SUBROUTINE instrumental(j)

            ! Argument Intent
             integer,             intent(inout) :: j         ! counters: line nb in DIFFaX and FAULTS respectively

            ! Local Variables
             integer             :: i,k,k1,k2,n                ! counters
             character(maxcharD) :: texte,key,ligne            ! manipulation de chaine de caracteres
             character(10)       :: A,B,C,D                    ! lecture de valeurs sur une meme ligne
             integer             :: ier                        ! control : 0 = successful, >0 error

            k=0
            k1=0
            k2=0
            n=0
            texte=" "
            key=" "
            ligne=" "


            do i=S1,S2-1
               j=j+1
               texte=adjustl(tfile(i))
               if(texte(1:1)=="!" .or. texte(1:1)==" ") then
                    ffile(j)=texte(1:maxcharF)
                    cycle
               else
                    n=n+1
                    k=index(texte," ")
                    key=U_Case(texte(1:k-1))
                    texte=adjustl(texte(k:))
                    if(n==1) then
                        if(texte==" ") then                      ! No comments after "INSTRUMENTAL"
                            ffile(j) = "  INSTRUMENTAL AND SIZE BROADENING"
                        else                                     ! Comment(s) present after "INSTRUMENTAL"
                            ligne = "  INSTRUMENTAL AND SIZE BROADENING    ! " // texte
                            ffile(j) = ligne(1:maxcharF)
                        end if
                    else if(n==2) then
                        ffile(j) = " !Type of radiation"
                        j=j+1
                        if(key=="X-RAY" .or. key=="NEUTRON" &
                           .or. key=="ELECTRON") then
                            if(texte==" ") then                         ! No comments after keyword
                                ligne = "  Radiation   " // trim(key)
                                ffile(j) = ligne(1:maxcharF)
                            else                                        ! Comment(s) present after "INSTRUMENTAL"
                                ligne = "  Radiation   " // trim(key) // "    " // texte
                                ffile(j) = ligne(1:maxcharF)
                            end if
                        else
                            ligne = " ! ******** Error: keyword " // key // " not recognized"
                            ffile(j) = ligne
                            write(unit=op,fmt="(3a,i4)") " => Error: keyword ", trim(key), " not recognized at line # ", tlineD(i)
                        end if
                    else if(n==3) then
                        ffile(j) = " !             Lambda1       Lambda2      ratio"
                        j=j+1
                        if(texte==" ") then                                                  ! No comments
                            ligne = "  Wavelength   " // trim(key) // "         0.0000       0.0000"
                            ffile(j) = ligne(1:maxcharF)
                        else                                                                 ! Comment(s) present
                            ligne = "  Wavelength   " // trim(key) // "         0.0000   0.0000    " // texte
                            ffile(j) = ligne(1:maxcharF)
                        end if
                    else if(n==4) then
                        ffile(j) = " !Instrumental aberrations      zero      sycos     sysin"
                        j=j+1
                        ffile(j) = "  Aberrations                   0.0000    0.0000    0.0000"
                        j=j+1
                        ffile(j) = "                                0.00      0.00      0.00"
                        j=j+1
                        ffile(j) = " !instr. broadening       u           v           w           x" &
                                   // "         Dg         Dl"
                        j=j+1
                        if(key=="GAUSSIAN" .or. key=="LORENTZIAN" &
                               .or. key=="PSEUDO-VOIGT") then
                            read(unit=texte,fmt=*,iostat=ier) A, B, C, D
                            ligne = "  "// trim(key) //"            "// trim(A) //"       "// trim(B) &
                                   //"       "// trim(C) //"       "// " 0.1  " //"       "// "8000.0" &
                                   //"       "// "8000.0"
                            k=index(texte,"T")
                            k1=index(texte,"!")
                            k2=index(texte,"{")
                            if(k/=0 .or. k1/=0 .or. k2/=0) then                                                    ! TRIM keyword or comment(s) present
                                texte = texte(k:)
                                ligne = trim(ligne) //"       "// trim(texte)
                                ffile(j) = ligne(1:maxcharF)
                                j=j+1
                                ffile(j) = "                          0.00        0.00        0.00        " &
                                           // "0.00        0.00        0.00          "
                            else                                                             ! No comment or TRIM keyword present
                                ffile(j) = ligne(1:maxcharF)
                                j=j+1
                                ffile(j) = "                          0.00        0.00        0.00        " &
                                           // "0.00        0.00        0.00          "
                            end if
                        else if(key=="NONE") then
                            ligne = "  PSEUDO-VOIGT             0.01       0.00       0.00       0.01       " &
                                    // "8000.0       8000.0     !TRIM"
                            k1=index(texte,"!")
                            k2=index(texte,"{")
                            if(k1/=0 .or. k2/=0) then                                                    ! comment(s) present
                                texte = texte(k:)
                                ligne = trim(ligne) //"       "// trim(texte)
                                ffile(j) = ligne(1:maxcharF)
                                j=j+1
                                ffile(j) = "                          0.00        0.00        0.00        " &
                                           // "0.00        0.00        0.00          "
                            else                                                             ! No comment or TRIM keyword present
                                ffile(j) = ligne(1:maxcharF)
                                j=j+1
                                ffile(j) = "                          0.00        0.00        0.00        " &
                                           // "0.00        0.00        0.00          "
                            end if
                        else
                            ligne = " ! ******** Error: keyword " // key // " not recognized"
                            ffile(j) = ligne
                            write(unit=op,fmt="(3a,i4)") " => Error: keyword ", trim(key), " not recognized at line # ", tlineD(i)
                        end if
                    else
                           write(unit=op,fmt="(a,i4)") " => Error in INSTRUMENTAL section at line ", tlineD(i)
                           write(unit=ligne,fmt="(2a,i4)") " !Error in the INSTRUMENTAL section " &
                                                           // "of the DIFFAX input file at line", tlineD(i)
                           ffile(j) = ligne(1:maxcharF)
                    end if
               end if
            end do

        j=j+1

         return
      END SUBROUTINE instrumental

      SUBROUTINE structural(j)

            ! Argument Intent
             integer,             intent(inout) :: j         ! counters: line nb in DIFFaX and FAULTS respectively

            ! Local Variables
             integer             :: i,k,k1,n                   ! counters
             character(maxcharD) :: texte,key,ligne            ! manipulation de chaines de caracteres
             integer             :: ier                        ! control : 0 = successful, >0 error
             character(25)       :: A, B                       ! Lwidth values

            k=0
            k1=0
            n=0
            texte=" "
            key=" "
            ligne=" "
            A=" "
            B=" "

           do i=S2,S3-1
               j=j+1
               texte=adjustl(tfile(i))
               if(texte(1:1)=="!" .or. texte(1:1)==" ") then
                    ffile(j)=texte(1:maxcharF)
                    cycle
               else
                    n=n+1
                    k=index(texte," ")
                    key=U_Case(texte(1:k-1))
                    texte=adjustl(texte(k:))
                    if(n==1) then
                        if(texte==" ") then                      ! No comments after "INSTRUMENTAL"
                            ffile(j) = "  STRUCTURAL"
                        else                                     ! Comment(s) present after "INSTRUMENTAL"
                            ligne = "  STRUCTURAL    ! " // texte
                            ffile(j) = ligne(1:maxcharF)
                        end if
                        j=j+1
                        ffile(j+1) = "!!Layers to be stacked to calculate the average cell (for Bragg positions in PRF file)"
                        ffile(j+2) = "! CalcAverCell 'nb_of_layers'"
                        ffile(j+3) = "! 'sequence of layers'"
                        ffile(j+4) = "!!Average cell parameters (for Bragg positions in PRF file)"
                        ffile(j+5) = "! AverCell 'a b c alpha beta gamma'"
                        ffile(j+6) = "!!Average Space Group (for Bragg positions in PRF file)"
                        ffile(j+7) = "! SpGR 'space_group'"
                        ffile(j+8) = "!!FullProf Studio commands"
                        ffile(j+9) = "! FST_CMD 'FullProf Studio Commands'"
                        j=j+10
                    else if (n==2) then
                        ffile(j) = " !         a        b        c        gamma"
                        j=j+1
                        ligne = "  Cell     " // trim(key) // "   " // texte
                        ffile(j) = ligne(1:maxcharF)
                        j=j+1
                        ffile(j) = "           0.00     0.00     0.00     0.00"
                    else if (n==3) then
                        ffile(j) = " !Laue symmetry"
                        j=j+1
                        ligne = "  Symm     " // trim(key) // "   " // texte
                        ffile(j) = ligne(1:maxcharF)
                    else if (n==4) then
                        ffile(j) = " !Number of layer types"
                        j=j+1
                        ligne = "  NLAYERS     " // trim(key) // "   " // texte
                        ffile(j) = ligne(1:maxcharF)
                        read(unit=key,fmt="(i3)")  nlayers
                    else if (n==5) then
                        ffile(j) = " !Layer width"
                        j=j+1
                        texte=adjustl(tfile(i))
                        k=index(texte," ")
                        key=U_Case(texte(1:k-1))
                        if(key=="INFINITE") then
                            ligne = "  Lwidth     " // texte
                            ffile(j) = ligne(1:maxcharF)
                        else
                            k=index(texte,"!")                 ! check if comments are present after Lwidth values
                            k1=index(texte,"{")
                            key=" "
                            if (k/=0 .and. (k<=k1  .or. k1==0)) then
                                key=texte(1:k-1)
                                texte=texte(k:)
                            else if(k1/=0 .and. (k>k1 .or. k==0)) then
                                key=texte(1:k1-1)
                                texte=texte(k1:)
                            else
                                key=texte
                                texte=" "
                            end if
                            read(unit=key,fmt=*,iostat=ier) A, B
                            if(ier/=0) then
                                read(unit=key,fmt=*,iostat=ier) A
                                B = A
                            end if
                            ligne = "  Lwidth      " // trim(A) // "      "// trim(B) // "       " // texte
                            ffile(j) = ligne(1:maxcharF)
                            j=j+1
                            ffile(j) = "              0.00     0.00"
                        end if
                    else
                        write(unit=op,fmt="(a,i4)") " => Error in STRUCTURAL section at line ", tlineD(i)
                        write(unit=ligne,fmt="(a,i4)") "! ******** Error in the STRUCTURAL section " &
                                                       // "of the DIFFAX input file at line ", tlineD(i)
                        ffile(j) = ligne(1:maxcharF)
                    end if
                end if
            end do

            if(n<5) then                                ! if there was not any Lwidth in the DIFFaX input file
                ffile(j) = " !Layer width"
                j=j+1
                ffile(j) = "!  ******** The layers' width was not defined in the DIFFaX input file. It is set to INFINITE."
                j=j+1
                ffile(j) = "  Lwidth   INFINITE"
            end if

        j=j+1

         return
      END SUBROUTINE structural

      SUBROUTINE layers(j)

            ! Argument Intent
             integer,             intent(inout) :: j         ! counters: line nb in DIFFaX and FAULTS respectively

            ! Local Variables
             integer                :: i,k,k1,n                     ! counter
             character(maxcharD)    :: texte,key,ligne            ! manipulation de chaines de caracteres
             integer, dimension(50) :: L                          ! indexation de lignes
             integer                :: nactlay                    ! nb de layers comptées
             character(4)           :: atom, atom1, atom2         ! for manipulation of atoms' names

            k=0
            k1=0
            n=0
            texte=" "
            key=" "
            ligne=" "
            L=0

            do i=S3,S4-1                                        ! index lines containing the keyword LAYER
                k=index(tfile(i),"LAYER")
                if (k/=0) then
                    n=n+1
                    L(n)=i
                else
                    cycle
                end if
            end do

            L(n+1)=S4
            nactlay=n
            n=0

           if (nactlay /= nlayers) then
                write(unit=op,fmt="(a,i3,2a,i3,a)") " => Error: ", nactlay, " layers were found in the LAYER section",&
                                                    "    while", nlayers, "were declared in the STRUCTURAL section"
           else
                write(unit=op,fmt="(a,i3,a)") " => LAYERS section:", nactlay, " layers were found in the LAYER section", &
                                              "    while", nlayers, " layers were declared in the STRUCTURAL section"
           end if

           i=S3-1
           do
               i=i+1
               if(i>S4-1) exit
               j=j+1
               texte=adjustl(tfile(i))
               if(texte(1:1)=="!" .or. texte(1:1)==" ") then
                    ffile(j)=texte(1:maxcharF)
                    cycle
               else
                    k=index(texte," ")
                    key=U_Case(texte(1:k-1))
                    texte=adjustl(texte(k:))
                    if(key == "LAYER") then
                        n=n+1
                        k1=index(texte, "=")
                        if(k1 /= 0) then
                            ligne = "  " // trim(key) // " " // trim(texte)
                            ffile(j) = ligne(1:maxcharF)
                     !      cycle
                        else
                            ligne = "  " // trim(key) // " " // trim(texte)
                            ffile(j) = ligne(1:maxcharF)
                            do i=L(n)+1,L(n+1)-1
                                j=j+1
                                texte = adjustl(tfile(i))
                                if(texte(1:1)=="!" .or. texte(1:1)==" ") &
                                    then
                                    ffile(j)=texte(1:maxcharF)
                                    cycle
                                else
                                    k=index(texte," ")
                                    key=U_Case(texte(1:k-1))
                                    if(key == "NONE" .or. key == "CENTROSYMMETRIC") &
                                        then
                                        ffile(j) = " !Layer symmetry"
                                        j=j+1
                                        ligne = "  LSYM   " // key
                                        ffile(j) = ligne(1:maxcharF)
                                        j=j+1
                                        ffile(j) = " !Atom name   Number      x        y        z        Biso        Occ "
                                        cycle
                                    else
                                        atom = texte(1:4)
                                        texte = adjustl(texte(5:))
                                        k = index(atom," ")
                                        if(k/=0) then
                                            atom1 = atom(1:k-1)
                                            atom2 = atom(k+1:4)
                                            atom  = trim(atom1) // trim(atom2) // "  "
                                        end if
                                        ligne = "  Atom  " // atom // "  " // texte
                                        ffile(j) = ligne(1:maxcharF)
                                        j=j+1
                                        ffile(j) = "                          0.00     0.00" &
                                                // "     0.00     0.00     0.00"
                                    end if
                                end if
                            end do
                            i=i-1
                        end if
                    else
                        ligne = " ! ******** Error: keyword " // key // " not recognized"
                        ffile(j) = ligne
                        write(unit=op,fmt="(3a,i4)") " => Error: keyword ", trim(key), " not recognized at line # ", tlineD(i)
                    end if
                end if
            end do

        j=j+1

        return
      END SUBROUTINE layers

      SUBROUTINE stacking(j)

            ! Argument Intent
             integer,             intent(inout) :: j           ! counters: line nb in DIFFaX and FAULTS respectively

            ! Local Variables
             integer             :: i,k,k1,n                   ! counters
             character(maxcharD) :: texte,key,ligne            ! manipulation de chaines de caracteres

            k=0
            k1=0
            n=0
            texte=" "
            key=" "
            ligne=" "

            i=S4-1

           do
               i=i+1
               if(i>=S5) return
               j=j+1
               texte=adjustl(tfile(i))
               if(texte(1:1)=="!" .or. texte(1:1)==" ") then
                    ffile(j)=texte(1:maxcharF)
                    cycle
               else
                    k=index(texte," ")
                    key=U_Case(texte(1:k-1))
                    texte=adjustl(texte(k:))
                    if(key=="STACKING") then
                        ligne = "  STACKING     " // texte
                        ffile(j) = ligne(1:maxcharF)
                    else if(key=="EXPLICIT") then
                        ffile(j) = " !Stacking type"
                        j=j+1
                        ligne = "  EXPLICIT     " // texte
                        ffile(j) = ligne(1:maxcharF)
                        i=i+1
                        j=j+1
                        texte=adjustl(tfile(i))
                        k=index(texte," ")
                        key=U_Case(texte(1:k-1))
                        if (key == "RANDOM") then
                            texte=texte(k+1:)
                            ffile(j) = " !             Number of layers"
                            j=j+1
                            ligne = "  RANDOM       " // texte
                            ffile(j) = ligne(1:maxcharF)
                        else                            ! specific sequence of layers
                            ffile(j) = " SPECIFIC"
                            j=j+1
                            do
                                texte=adjustl(tfile(i))
                                ffile(j)=texte(1:maxcharF)
                                i=i+1
                                j=j+1
                                if(i>=S5) return
                            end do
                        end if
                    else if(key=="RECURSIVE") then
                        ffile(j) = " !Stacking type"
                        j=j+1
                        ligne = "  RECURSIVE    " // texte
                        ffile(j) =  ligne(1:maxcharF)
                        j=j+1
                        ffile(j) = " !Number of layers"
                        i=i+1
                        j=j+1
                        texte=adjustl(tfile(i))
                        k=index(texte," ")
                        key=U_Case(texte(1:k-1))
                        if(key=="INFINITE") then
                            ffile(j) = "  " // texte(1:maxcharF-2)
                        else
                            ffile(j) = "  " // texte(1:maxcharF-2)
                            j=j+1
                            ffile(j) = "   0.00"
                        end if
                    end if
                end if
            end do

        j=j+1

        return
      END SUBROUTINE stacking

      SUBROUTINE transitions(j)

            ! Argument Intent
             integer,             intent(inout) :: j           ! counters: line nb in DIFFaX and FAULTS respectively

            ! Local Variables
             integer             :: i,k,k1,n,l                 ! counter
             character(maxcharD) :: texte,key,ligne,more       ! manipulation de chaines de caracteres
             character(25)       :: A, B, C                    ! Lwidth values

            k=0
            k1=0
            n=0
            l=1
            texte=" "
            key=" "
            ligne=" "
            more=" "
            A=" "
            B=" "
            C=" "

           do i=S5,nblinesT
               j=j+1
               texte=adjustl(tfile(i))
               if(texte(1:1)=="!" .or. texte(1:1)==" ") then
                    ffile(j)=texte(1:maxcharF)
                    cycle
               else
                    k=index(texte," ")
                    key=U_Case(texte(1:k-1))
                    texte=adjustl(texte(k:))
                    if(key=="TRANSITIONS") then
                        ligne = "  TRANSITIONS     " // texte
                        ffile(j) = ligne(1:maxcharF)
                        j=j+1
                        ffile(j) = " !LT   alpha      tx       ty      tz"
                        j=j+1
                        ffile(j) = " !FW      C11     C22      C33     " &
                                   // "C12      C23      C31"
                    else
                        n=n+1
                        if(n>nlayers) then
                            l=l+1
                            n=1
                            if(l>nlayers) cycle
                        end if
                        k=index(texte," ")
                        A=texte(1:k-1)
                        texte=adjustl(texte(k:))
                        k=index(texte," ")
                        B=texte(1:k-1)
                        texte=adjustl(texte(k:))
                        k=index(texte," ")
                        C=texte(1:k-1)
                        texte=adjustl(texte(k:))
                        k=index(texte,"!")
                        k1=index(texte,"{")
                        if(k==1 .or. k1==1) then
                            more=" "
                        else if(k>1 .and. (k<=k1 .or. k1==0)) then
                            more=texte(1:k-1)
                            texte=texte(k:)
                        else if(k1>1 .and. (k>k1 .or. k==0)) then
                            more=texte(1:k1-1)
                            texte=texte(k1:)
                        else  ! k==0 .and. k1==0
                            more=texte
                            texte=" "
                        end if
                        write(unit=ligne,fmt="(a,i2,a,i2)") " !Layer ", l, " to layer ", n
                        ffile(j) = ligne(1:maxcharF)
                        j=j+1
                        ligne = "  LT  " // trim(key) // "    " // trim(A) // "    " &
                                // trim(B) // "    " // trim(C) // "    " // trim(texte)
                        ffile(j) = ligne(1:maxcharF)
                        j=j+1
                        ffile(j) = "      0.00      0.00      0.00      0.00"
                        j=j+1
                        if(more /= " ") then
                            ligne = "  FW  " // trim(more)
                        else
                            ligne = "  FW  0.00000   0.00000   0.00000" &
                                   // "   0.00000   0.00000   0.00000"
                        end if
                        ffile(j) = ligne(1:maxcharF)
                        j=j+1
                        ffile(j) = "      0.00      0.00      0.00" &
                                   // "      0.00      0.00      0.00"
                    end if
                end if
            end do

        j=j+1

        return
      END SUBROUTINE transitions

      SUBROUTINE calculation(j)

            ! Argument Intent
             integer,             intent(inout) :: j         ! counters: line nb in ffile

            ffile(j+1) = "  CALCULATION"
            ffile(j+2) = "  Simulation"
            ffile(j+3) = " !Type        th2_min  th2_max  d_theta"
            ffile(j+4) = "  Powder      5.0      120.0    0.02"
            ffile(j+5) = "! Replace_Files"
            j=j+6

            ffile(j+1) = "! Other available options for the CALCULATION section (to be filed by the user): "
            j=j+2

            ffile(j+1) = "! CALCULATION"
            ffile(j+2) = "! Simulation"
            ffile(j+3) = "!!Type        i_plane  l_upper  loglin  brightness"
            ffile(j+4) = "! SADP        1        4.0      0       10.0"
            ffile(j+5) = "! Replace_Files"
            j=j+6

            ffile(j+1) = "! CALCULATION"
            ffile(j+2) = "! LMA"
            ffile(j+3) = "! Corrmax   30"
            ffile(j+4) = "! Maxfun    2400"
            ffile(j+5) = "! Tol       0.1E-04"
            ffile(j+6) = "! Nprint    0"
            ffile(j+7) = "! Replace_Files"
            j=j+8

            ffile(j+1) =  "! EXPERIMENTAL"
            ffile(j+2) =  "!!Filename                    Scale factor     code"
            ffile(j+3) =  "! FILE             CFILE.dat  1.0              0.00"
            ffile(j+4) =  "! Excluded_Regions 0"
            ffile(j+5) =  "! FFORMAT          free  "
            ffile(j+6) =  "!!Linear interpolation"
            ffile(j+7) =  "! Bgrinter    4.bgr"
            ffile(j+8) =  "!!Polynomial  Number of coefficients "
            ffile(j+9) =  "! Bgrcheb     2"
            ffile(j+10) = "!!Polynomial coefficients "
            ffile(j+11) = "! 1.00000     0.20000"
            ffile(j+12) = "! 0.0         0.0  "
            ffile(j+13) = "!!Number of pattern backgrounds"
            ffile(j+14) = "! Bgrnum  1 "
            ffile(j+15) = "!!Pattern file     Filename    Scale factor   code"
            ffile(j+16) = "! Bgrpatt           CFILE.sub   1.0000         0.00    !CFILE.hkl"
            j=j+17

        return
      END SUBROUTINE calculation

      SUBROUTINE Write_flts()

        ! Local Variables
          integer             :: ier,ok      ! control : 0 = successful, >0 error
          integer             :: j           ! counters

        !open (create) the FAULTS input file
        OPEN (unit=f_out, file=filename, status='replace', action='write', position='rewind', iostat=ier)

            if (ier /= 0) then
                write(UNIT=op,FMT="(2a)") " => Error creating the FAULTS input file ", filename
            else
                write(UNIT=op,FMT="(3a)") " => The FAULTS input file ", trim(filename) , " was successfully created"
                write(UNIT=op,FMT="(3a)") " => Writing the FAULTS input file ", trim(filename), " ..."
            end if

            ok=0
            ier=0

            do j=1,nlfinalF
                write(unit=f_out,fmt="(a)",iostat=ier) ffile(j)
                !write(unit=op,fmt="(a)") tfile(j)
                ok=ok+ier
            end do

            write(UNIT=f_out, FMT="(a)") " !End of the FAULTS input file"
            !write(UNIT=f_out, FMT="(5i3)")  S1,S2,S3,S4,S5

        CLOSE (f_out)

            if (ok /= 0) then
                write(UNIT=op,FMT="(2a)") " => Error writing in the FAULTS input file ", filename
            else
                write(UNIT=op,FMT="(3a)") " => The FAULTS input file ", trim(filename) , " was successfully written"
                write(UNIT=op,FMT="(3a)") " => Normal end of DIFFaX_to_FAULTS"
            end if


      END SUBROUTINE Write_flts                    !!! CFML_String8Util.f90

      SUBROUTINE Ucase(line)
        !---- Argument ----!
        character (len=*), intent(in out) :: line

        line=u_case(line)
        return
      End SUBROUTINE Ucase

      FUNCTION U_Case(Text) Result (Mtext)         !!! CFML_String8Util.f90
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
      End FUNCTION U_Case

END PROGRAM DIFFaX_to_FAULTS
