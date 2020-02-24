  Module DataRed_Mod
    Use CFML_GlobalDeps
    Use CFML_gSpaceGroups
    Use CFML_Metrics
    Use CFML_Strings, only: File_Type, Reading_File
    Use CFML_Reflections
    Use Twin_Mod
    Implicit None
    !Everything is public
    integer, parameter :: nref=400000, inp=1, iou=2, ihkl=3, irej=4, iin=10
    integer, dimension(nref)             :: idomain
    integer, dimension(nref)             :: itreat, iord, nequv, ini, fin, numor, ivk, icod, warn
    integer, dimension(nmax)             :: contr
    integer, dimension(:),   allocatable :: hkl
    real(kind=cp),    dimension(3,nref)  :: h
    real(kind=cp),    dimension(nref)    :: intav, sigmav, sigstat,lambda_laue, &
                                            tbar,twotet,absorpt
    real(kind=cp),    dimension(3)       :: h1,h2,h3
    real(kind=cp),    dimension(3)       :: kv
    real(kind=cp),    dimension(3,48)    :: kvec
    real(kind=cp),    dimension(4)       :: angles
    real(kind=cp),    dimension(256)     :: weight  ! A maximum of 256 reflections equivalent to one of them can be treated
    real(kind=cp),    dimension(3,nmax)  :: hkln

    character(len=512)          :: filein, fileout, filered
    character(len=256)          :: line, cmdline, title
    character(len=50)           :: forma
    character(len=20)           :: line1,wmess
    character(len=1)            :: ans
    !integer, dimension(:,:), allocatable :: hkls !indices for superspace case

    Type, public :: Conditions_Type
      logical :: Friedel    =.true.
      logical :: prop       =.false.
      logical :: file_given =.false.
      logical :: cell_given =.false.
      logical :: twinned    =.false.
      logical :: wave_given =.false.
      logical :: powder     =.false.
      logical :: domain     =.false.
      logical :: transf_ind =.false.
      logical :: hkl_real   =.false.
      logical :: eps_given  =.false.
      logical :: statistics =.false.
      logical :: scale_given=.false.
      logical :: absent     =.false.
      logical :: lauetof    =.false.
      character(len=:), allocatable :: fileout
      integer                       :: hkl_type
      real(kind=cp), dimension(3,3) :: transhkl
      real(kind=cp), dimension(3)   :: cel,ang
      real(kind=cp)                 :: epsg, wavel
      real(kind=cp)                 :: warning=0.25  ! 25% error for warning equivalent reflections
    End Type Conditions_Type

    Type, extends(Refl_Type)    :: ObsRef
      real(kind=cp),dimension(3):: hr
      real(kind=cp)             :: intens
      real(kind=cp)             :: sigma
      real(kind=cp)             :: twtheta
      real(kind=cp)             :: omega
      real(kind=cp)             :: chi
      real(kind=cp)             :: phi
    End Type ObsRef

    Type(Refl_Type), dimension(:), allocatable :: Reflex
    Type(SPG_Type)              :: SpG, twSpG
    Type(SuperSpaceGroup_Type)  :: SSpG
    Type(Cell_G_Type)           :: celda
    !Type(Group_k_Type)          :: Gk

    character(len=*),parameter,dimension(0:1) :: warn_mess=["                      ",  &
                                                             " <- Dubious reflection"]
    real    :: total, sig, suma, suman, sumaw, sumanw, Rint, Rwint, aver_sig, &
               wavel,sigg, int_rej, epsg, delt, warning, scal_fact, q2, aver_int
    integer :: i,j,k, ier, nr=0, iv, ns, rej, ival,  nkv, nn, &
               lenf, nin, hkl_type, ivp, nk, nequiv, L, drej, Lmin,i_refout, a,b

    contains

    Subroutine Read_DataRed_File(cfl,cond,twins)
      type(File_Type),       intent(in)  :: cfl
      type(Conditions_Type), intent(out) :: cond
      type(Twin_Type),       intent(out) :: twins

      !--- Local variables ---!
      integer :: i,j
      character(len=:), allocatable :: keyw
      character(len=:), allocatable :: line

!  TITLE   YMn2O5 Phase comm
!  INPFIL  ymnm25K.col
!  OUTFIL  ymnm25K
!  WAVE    2.34
!  CELL    7.2570 8.4918 5.6780 90 90 90
!  SPGR    P -1
!  EPSIL   0.0001
!  KVEC    0.5 0 0.25
!  HKL_T   5


      do i=1,cfl%nlines
         line=adjustl(cfl%line(i)%str)
         if(len_trim(line) == 0) cycle
         if(line(1:1) == "!" .or. line(1:1) == "#") cycle
         j=index(line," ")
         if(j == 0) then
           keyw=u_case(trim(line))
         else
           keyw=u_case(trim(line(1:j-1)))
         end if


         Select Case(keyw)

            Case("OUTFIL")
              cond%fileout=adjustl(line(j:))
              prop=.true.

            Case("KVEC")
              read(unit=line(6:),fmt=*,iostat=ier) kv
              if(ier /= 0) then
                write(unit=*,fmt="(a)") " => Error reading the propagation vector"
                stop
              end if
              prop=.true.

            Case("TRANS")
              cond%transf_ind=.true.
              read(unit=line(j:),fmt=*)  ((cond%transhkl(i,j),j=1,3),i=1,3)

            Case("SCALE")
              scale_given=.true.
              read(unit=line(j:),fmt=*)  cond%scal_fact

            Case("STATI")
              cond%statistics=.true.

            Case("SPGR")
              line1=adjustl(line(6:))

            Case("EPSIL")
              read(unit=line(j:),fmt=*)  cond%epsg
               cond%eps_given=.true.

            Case("NFRDL")
               cond%Friedel=.false.

            Case("CELL")
               read(unit=line(j:),fmt=*)  cel, ang
               call Set_Crystal_Cell(cel,ang,Celda)
               cond%cell_given=.true.

            Case("WAVE")
               read(unit=line(j:),fmt=*)  cond%wavel
               cond%wave_given=.true.

            Case("HKL_T")
               read(unit=line(j:),fmt=*)  cond%hkl_type

            Case("WARN")
               read(unit=line(j:),fmt=*)  cond%warning

            Case("POWDER")
               cond%powder=.true.

            Case("TWIN")
              call read_twinlaw(iin)
              twinned=.true.
              if(iubm) lambda=wavel

            Case("DOMAIN")
              domain=.true.

         End Select

      end do

    End Subroutine Read_DataRed_File

    Subroutine Header_Output()
       write(unit=iou,fmt="(a)") "       ==============================="
       write(unit=iou,fmt="(a)") "       DATA REDUCTION PROGRAM: DataRed"
       write(unit=iou,fmt="(a)") "       ==============================="
       write(unit=iou,fmt="(a)") "          JRC-ILL version:30-2-2020"
       write(unit=iou,fmt="(a)") "       "

       write(unit=iou,fmt="(a,a)")   " Input        Control  file: ", trim(filered)
       write(unit=iou,fmt="(a,a)")   " Input    Reflections  file: ", trim(filehkl)
       write(unit=iou,fmt="(a,a)")   " Averaged Reflections  file: ", trim(fileout)//".int"
       write(unit=iou,fmt="(a,a//)") " Rejected Reflections  file: ", trim(fileout)//".rej"
       write(unit=iou,fmt="(a,a)")   " General       Output  file: ", trim(fileout)//".out"

       Select Case(hkl_type)
         Case(0)       !Shelx-like input file (3i4,2f8.2)
            write(unit=iou,fmt="(a/)") " Format of the reflections file =>  ShelX-like format h k l intens sigma (3i4,2f8.2)"
         case(1)
            write(unit=iou,fmt="(a/)") " Format of the reflections file =>  FREE format:  h k l intens sigma "
         case(2)
            write(unit=iou,fmt="(a)")  " Data from COLL5 (two lines per reflection) *.fsq           "
            write(unit=iou,fmt="(a)")  " Format of the reflections file =>  (i6,3f7.3,2f10.2,3f8.2) "
            write(unit=iou,fmt="(a/)") " For reading the items: Numor, h k l Int Sigma gamma nu phi"
         case(3)
            write(unit=iou,fmt="(a)")  " Data in FullProf *.int format           "
            write(unit=iou,fmt="(a)")  " Format of the reflections read in the file "
            write(unit=iou,fmt="(a/)") " For reading the items:  h k l Int Sigma angles(1:4)"
         case(4)
            write(unit=iou,fmt="(a)")  " Data from COLL5 (1 line per reflection) *.col (hkl-integer)      "
            write(unit=iou,fmt="(a)")  " Format of the reflections file =>  (i6,3i4,2f10.2,4f8.2) "
            write(unit=iou,fmt="(a/)") " For reading the items: Numor, h k l Int Sigma theta omega chi phi"
         case(5)
            write(unit=iou,fmt="(a)")  " Data from COLL5 (1 line per reflection) *.col (hkl-real)"
            write(unit=iou,fmt="(a)")  " Format of the reflections file =>  "//trim(forma)
            write(unit=iou,fmt="(a/)") " For reading the items: Numor, h k l Int Sigma theta omega chi phi"
         case(6)
            write(unit=iou,fmt="(a)")  " Data from COLL5 (1 line per reflection) *.col (hkl-real) 6T2-LLB"
            write(unit=iou,fmt="(a)")  " Format of the reflections file =>  (i4,3f6.2,2f10.2,4f8.2) "
            write(unit=iou,fmt="(a/)") " For reading the items: Numor, h k l Int Sigma theta omega chi phi"
         case(7)
            write(unit=iou,fmt="(a)")  " Data from SXD in format adequate for FullProf"
            write(unit=iou,fmt="(a,a)")" Format of the reflections file => " , trim(forma)
            if(prop) then
                write(unit=iou,fmt="(a/)") " For reading:  h k l ivk Fsqr s(Fsqr) Cod  Lambda Twotheta Absorpt. Tbar "
            else
                write(unit=iou,fmt="(a/)") " For reading:  h k l Fsqr s(Fsqr) Cod  Lambda Twotheta Absorpt. Tbar "
            end if
         case(8)
            write(unit=iou,fmt="(a)")  " Data from D3 file (1 line per reflection) *.bpbres (hkl-real) D3-ILL"
            write(unit=iou,fmt="(a)")  " (Warning: if the orientation of the crytal is not in a  "
            write(unit=iou,fmt="(a)")  "  symmetry direction, use 'P 1' as space group)"
            write(unit=iou,fmt="(a)")  " Format of the Flipping Ratio reflections file =>  (i8,6f8.3,2f10.2) "
            write(unit=iou,fmt="(a/)") " For reading the items: Numor, h k l omega gamma nu R Sigma(R) "
         case(9)
            write(unit=iou,fmt="(a)")  " Data from Modified COLL5 (1 line per reflection) *.col (hkl-real)"
            write(unit=iou,fmt="(a)")  " Format of the reflections file =>  (i6,3f10.4,2f10.2,4f8.2) "
            write(unit=iou,fmt="(a/)") " For reading the items: Numor, h k l Int Sigma theta omega chi phi"
         case(10)
            write(unit=iou,fmt="(a)")  " Data from SHELX HKLF5-format (domains have been treated outside DataRed)"
            write(unit=iou,fmt="(a)")  " Format of the reflections file =>  (3i4,2f8.0,i4) "
            write(unit=iou,fmt="(a/)") " For reading the items: h k l Int Sigma domain_code"
       End Select


       write(unit=iou,fmt="(a,i6)") " => Number of reflections read: ",nr

       if(scale_given) then
         write(unit=iou,fmt="(a,f8.4,a)") " => A scale factor ",scal_fact," has been applied to intensities and sigmas "
         intens(1:nr)=intens(1:nr)*scal_fact
         sigma(1:nr)= sigma(1:nr)*scal_fact
       end if

       if(hkl_real) then
          write(unit=*,fmt="(a)")   " => hkl indices will be treated as real numbers "
          write(unit=iou,fmt="(a)") " => hkl indices will be treated as real numbers "
       end if

       if(prop) write(unit=iou,fmt="(a,3f8.4)") " => Input propagation vector: ", kv

       if(hkl_type /= 10) then
         if(statistics) then
            write(unit=*,fmt="(a)")   &
            " => Statistical errors are considered for sigmas of average intensisites (propagation error formula)"
            write(unit=iou,fmt="(a)") &
            " => Statistical errors are considered for sigmas of average intensisites (propagation error formula)"
         else
            write(unit=*,fmt="(a)")   &
            " => Statistics is NOT considered for sigmas of average intensities: exp. variance weighted with 1/sigmas^2"
            write(unit=iou,fmt="(a)") &
            " => Statistics is NOT considered for sigmas of average intensities: exp. variance weighted with 1/sigmas^2"
         end if
       end if

       if(transf_ind) then
          write(unit=iou,fmt="(a)") " => Input indices trasformed with matrix: "
          write(unit=iou,fmt="(a,3f7.2,a)") "         (Hnew)     (",transhkl(1,:)," )   (Hold)"
          write(unit=iou,fmt="(a,3f7.2,a)") "         (Knew)  =  (",transhkl(2,:)," )   (Kold)"
          write(unit=iou,fmt="(a,3f7.2,a)") "         (Lnew)     (",transhkl(3,:)," )   (Lold)"

          if(prop) then
            kv=matmul(transhkl,kv)
            write(unit=iou,fmt="(a,3f8.4)") " => Transformed propagation vector: ", kv
          end if
       end if

    End Subroutine Header_Output

  End Module DataRed_Mod
