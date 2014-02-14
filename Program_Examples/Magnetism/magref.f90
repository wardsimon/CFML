Program MagRef
 use CFML_crystallographic_symmetry,only: space_group_type, Write_SpaceGroup
 use CFML_Atom_TypeDef,             only: Atom_List_Type, Write_Atom_List, MAtom_list_Type
 use CFML_crystal_metrics,          only: Crystal_Cell_Type, Write_Crystal_Cell
 use CFML_Reflections_Utilities,    only: Hkl_s
 use CFML_IO_Formats,               only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
 use CFML_Propagation_vectors,      only: K_Equiv_Minus_K
 use CFML_Magnetic_Symmetry
 use CFML_Magnetic_Structure_Factors

 implicit none

 type (file_list_type)       :: fich_cfl
 type (Space_Group_Type)     :: SpG
 type (MagSymm_k_Type)       :: MGp
 type ( Atom_list_Type)      :: A
 type (MAtom_list_Type)      :: Am
 type (Crystal_Cell_Type)    :: Cell
 type (MagH_Type)            :: Mh

 character(len=256)          :: filcod     !Name of the input file
 character(len=1)            :: sig
 real                        :: sn,sf2
 integer                     :: Num, lun=1, ier,i,j,m,ih,ik,il,iv, n_ini,n_end
 complex                     :: fc

 integer                     :: narg,iargc
 Logical                     :: esta, arggiven=.false.
      !---- Arguments on the command line ----!
      narg=COMMAND_ARGUMENT_COUNT()

      if(narg > 0) then
              call GET_COMMAND_ARGUMENT(1,filcod)
              arggiven=.true.
      end if

      write(unit=*,fmt="(/,/,6(a,/))")                                                 &
           "              ------ P r o g r a m     M a g R e f  ------"               , &
           "                    ---- Version 0.2 April-2005 ----"                    , &
           "    ******************************************************************"  , &
           "    * Calculates magnetic structure factors, and magnetic interaction*"  , &
           "    * vectors from magnetic structures by reading a *.CFL file       *"  , &
           "    ******************************************************************"  , &
           "                    (JRC- April 2005, testing stage )"
    write(unit=*,fmt=*) " "

     if(.not. arggiven) then
       write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.cfl (give xx): "
       read(unit=*,fmt="(a)") filcod
       if(len_trim(filcod) == 0) stop
     end if

     open(unit=lun,file=trim(filcod)//".cal", status="replace",action="write")
     write(unit=lun,fmt="(/,/,6(a,/))")                                                &
           "              ------ P r o g r a m     M a g R e f  ------"               , &
           "                    ---- Version 0.2 April-2005 ----"                    , &
           "    ******************************************************************"  , &
           "    * Calculates magnetic structure factors, and magnetic interaction*"  , &
           "    * vectors from magnetic structures by reading a *.CFL file       *"  , &
           "    ******************************************************************"  , &
           "                    (JRC- April 2005, testing stage )"

     inquire(file=trim(filcod)//".cfl",exist=esta)
     if( .not. esta) then
       write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cfl does'nt exist!"
       stop
     end if
     call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)

     If(err_form) then
       write(unit=*,fmt="(a)") trim(err_form_mess)
     else

     !Test of the type "file_list_type" by writing at the end of the file
     write(unit=lun,fmt="(/,a)") "    =========================="
     write(unit=lun,fmt="( a )") "    Text of the input CFL file"
     write(unit=lun,fmt="(a,/)") "    =========================="
     write(unit=*,fmt="(/,a)") "    =========================="
     write(unit=*,fmt="( a )") "    Text of the input CFL file"
     write(unit=*,fmt="(a,/)") "    =========================="
     do i=1,fich_cfl%nlines
       write(unit=lun,fmt="(a,i5,a)") " Line:",i,"  "//fich_cfl%line(i)
       write(unit=*,fmt="(a,i5,a)") " Line:",i,"  "//fich_cfl%line(i)
     end do


       call Write_Crystal_Cell(Cell,lun)
       call Write_SpaceGroup(SpG,lun,full=.true.)
       call Write_Atom_List(A,lun=lun)

       n_ini=1
       n_end=fich_cfl%nlines
       call Readn_Set_Magnetic_Structure(fich_cfl,n_ini,n_end,MGp,Am)
       if(err_MagSym) then
         write(unit=*,fmt="(a)") "   "//err_MagSym_Mess
         stop
       end if
       call Write_Magnetic_Structure(lun,MGp,Am)
    !!---- TYPE :: MAGSYMM_K_TYPE
    !!--..
    !!---- Type, Public :: MagSymm_k_Type
    !!----    character(len=31)                   :: MagModel   ! Name to characterize the magnetic symmetry
    !!----    character(len=1)                    :: Latt       ! Symbol of the crystallographic lattice
    !!----    integer                             :: nmsym      ! Number of magnetic operators per crystallographic operator
    !!----    integer                             :: centred    ! =0 centric centre not at origin, =1 acentric, =2 centric (-1 at origin)
    !!----    integer                             :: mcentred   ! =1 Anti/a-centric Magnetic symmetry, = 2 centric magnetic symmetry
    !!----    integer                             :: nkv        ! Number of independent propagation vectors
    !!----    real(kind=sp),dimension(3,12)       :: kvec       ! Propagation vectors
    !!----    integer                             :: Num_Lat    ! Number of centring lattice vectors
    !!----    real(kind=sp), dimension(3,4)       :: Ltr        ! Centring translations
    !!----    integer                             :: Numops     ! Reduced number of crystallographic Symm. Op.
    !!----    integer                             :: Multip     ! General multiplicity of the space group
    !!----    character(len=40),   dimension(48)  :: SymopSymb  ! Alphanumeric Symbols for SYMM
    !!----    type( Sym_Oper_Type),dimension(48)  :: SymOp      ! Crystallographic symmetry operators
    !!----    character(len=40),   dimension(48,8):: MSymopSymb ! Alphanumeric Symbols for MSYMM
    !!----    type(MSym_Oper_Type),dimension(48,8):: MSymOp     ! Magnetic symmetry operators
    !!---- End Type MagSymm_k_Type
    !!----
    !!----  MAGH_TYPE
    !!----     Type, Public  :: MagH_Type
    !!----        logical              :: keqv_minus  !True if k equivalent to -k
    !!----        integer              :: mult    !Multiplicity of the reflection (useful for powder calculations)
    !!----        integer              :: num_k   !number of the propagation vector vk
    !!----        real                 :: signp   !+1 for -vk   and -1 for +vk
    !!----        real                 :: s       !sinTheta/Lambda
    !!----        real                 :: sqMiV   !Square of the Magnetic Interaction vector
    !!----        real, dimension(3)   :: H       ! H +/- k
    !!----        complex,dimension(3) :: MsF     !magnetic structure factor
    !!----        complex,dimension(3) :: MiV     !magnetic interaction vector
    !!----     End Type HR_Type
    !!----
    !!----    Define the scatering vector vector  H+k and the sign -1 for H+k and +1 for H-k.
    !!----    Includes the magnetic interaction vector MiV = Mper = M

      do

         write(unit=*,fmt="(a,i2,a)",advance="no") &
         " => Enter a magnetic reflections as 4 integers -> (h,k,l,m)=H+sign(m)*k(abs(m)): "
         read(unit=*,fmt=*) ih,ik,il,m
         if( m == 0) exit
         !construct partially the object Mh
         j=sign(1,m)
         sig="+"
         if( j == -1) sig="-"
         Mh%signp=real(-j)  ! sign "+" for H-k and "-" for H+k
         iv=abs(m)
         Mh%num_k=iv
         Mh%h= real((/ih,ik,il/)) - Mh%signp*MGp%kvec(:,iv)
         Mh%s = hkl_s(Mh%h,Cell)
         Mh%keqv_minus=K_Equiv_Minus_K(MGp%kvec(:,iv),MGp%latt)
         !Calculate magnetic structure factor and magnetic interaction vector
         call Calc_Magnetic_StrF_MiV(Cell,MGp,Am,Mh)
         write(unit=lun,fmt="(/,a,3i4,a,3f8.4,a)") "  Reflection: (",ih,ik,il,") "//sig//" (",MGp%kvec(:,iv),")"
         write(unit=lun,fmt="(a,3f8.4,a)")         "              (",Mh%h,")"
         write(unit=*,  fmt="(/,a,3i4,a,3f8.4,a)") "  Reflection: (",ih,ik,il,") "//sig//" (",MGp%kvec(:,iv),")"
         write(unit=*,  fmt="(a,3f8.4,a)")         "              (",Mh%h,")"
         write(unit=lun,fmt="(a,2(3f8.4,a))") "  Magnetic Structure   Factor : (",real(Mh%MsF),")+i(",aimag(Mh%MsF),") "
         write(unit=lun,fmt="(a,2(3f8.4,a))") "  Magnetic Interaction Vector : (",real(Mh%MiV),")+i(",aimag(Mh%MiV),") "
         write(unit=*,  fmt="(a,2(3f8.4,a))") "  Magnetic Structure   Factor : (",real(Mh%MsF),")+i(",aimag(Mh%MsF),") "
         write(unit=*,  fmt="(a,2(3f8.4,a))") "  Magnetic Interaction Vector : (",real(Mh%MiV),")+i(",aimag(Mh%MiV),") "
         write(unit=lun,fmt="(a,f12.5 )")     "  Square of Mag. int.  Vector : ",Mh%sqMiV
         write(unit=*  ,fmt="(a,f12.5 )")     "  Square of Mag. int.  Vector : ",Mh%sqMiV
       end do


       write(unit=*,fmt="(a)") " Normal End of program: MagRef "
       write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".cal"
     end if

     close(unit=lun)
End Program MagRef

