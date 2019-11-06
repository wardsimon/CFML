    Module read_data
       use CFML_GlobalDeps,               only : sp , cp
       use CFML_String_Utilities,         only : number_lines , reading_lines , findfmt ,iErr_fmt, getword, err_string, &
                                                 err_string_mess, getnum, U_case, cutst
       use CFML_Optimization_General,     only : Opt_Conditions_Type
       use CFML_LSQ_TypeDef,              only : LSQ_Conditions_type, LSQ_State_Vector_Type
       use CFML_Crystal_Metrics,          only : Crystal_Cell_Type, Set_Crystal_Cell
       use CFML_Reflections_utilities,    only : Reflect_Type, get_maxnumref, HKL_uni
       use CFML_Crystallographic_Symmetry,only : Space_Group_Type, set_SpaceGroup, ERR_Symm_Mess, &
                                                 Err_Symm
       use CFML_IO_Formats               ,only : Read_CIF_Cell, Read_CIF_HM, Read_CIF_Atom
       use CFML_Atom_TypeDef             ,only : Atom_List_Type, Allocate_Atom_List, Atom_Equiv_List_Type, &
                                                 Set_Atom_Equiv_List, Atom_Equiv_Type

       use diffax_mod
       implicit none

       !private
       !public subroutines
       !public:: read_structure_file, length, Read_Bgr_patterns
       Type(Crystal_Cell_Type),                        public :: acell
       Type(Reflect_Type),  dimension(:), allocatable, public :: reflections
       character(len=20),                              public :: spgsymb
       Type(Space_Group_Type),                         public :: SpG

       integer, parameter :: max_bgr_num=15 !maximum number of background patterns
       integer, parameter :: max_n_cheb =24 !maximum number of Chebychev coefficients
       integer            :: LGBLT=0        !number of global parameters
       integer(kind=2), dimension(17+max_bgr_num+max_n_cheb) :: Lglb
       integer(kind=2), dimension(5,max_a,max_l)             :: Latom
       integer(kind=2), dimension(10,max_l,max_l)            :: Ltrans

       real,            dimension(17+max_bgr_num+max_n_cheb) :: ref_glb
       real,            dimension(5,max_a,max_l)             :: ref_atom
       real,            dimension(10,max_l,max_l)            :: ref_trans

       real,            dimension(17+max_bgr_num+max_n_cheb) :: val_glb
       real,            dimension(5,max_a,max_l)             :: val_atom
       real,            dimension(10,max_l,max_l)            :: val_trans
       real,            dimension(max_bgr_num)               :: secnd_phases
       real                                                  :: area_f

!
       character(len=*),dimension(17),parameter :: names_glb = [ "Scale_Factor","zero_shift  ", &
                                                                 "sycos       ","sysin       ", &
                                                                 "U           ","V           ", &
                                                                 "W           ","x           ", &
                                                                 "Dg          ","Dl          ", &
                                                                 "cell_a      ","cell_b      ", &
                                                                 "cell_c      ","cell_gamma  ", &
                                                                 "diameter_a  ","diameter_b  ", &
                                                                 "num_layers  "]
!      Definition of crys_2d_type

       Type,Public :: crys_2d_type
         integer    :: rad_type                     !radiation type   >  rad_type
         real       :: lambda , lambda2 , ratio     !radiation wavelength    > lambda
         integer    :: broad                        !type of instrumental broadening>  blurring
         real       :: cell_a                       !cell parameter a   >cell_a
         real       :: cell_b                       !cell parameter b    >cell_b
         real       :: cell_c                       !cell parameter c   >cell_c
         real       :: cell_gamma                   !gamma   >cell_gamma
         real       :: layer_a = 0.0, layer_b = 0.0 !layer characteristic widths (optional)
         real       :: zero_shift = 0.0
         real       :: sycos = 0.0
         real       :: sysin = 0.0
         real       :: p_u = 0.0                    !u  > pv_u
         real       :: p_v = 0.0                    !v  > pv_v
         real       :: p_w = 0.0                    !w  > pv_w
         real       :: fwhm = 0.0                   !fwmh ( cannot be negative)  >fwhm
         real       :: p_gamma = 0.0                !pseudo-voigt gamma    > pv_gamma
         real       :: p_x                          ! x > pv_x
         real       :: p_dg, p_dl                   ! gaussian and lorentzian average volumetric size
         real       :: tolerance                    !d-> Maximum deviation that intensities can have
         real       :: bckg_level = 100.0           ! Background level for powder simulations
         character(len=132)  :: sym         !Laue symmetry  >pnt_grp
         character(len=132)  :: calctype    !type of calculation
         real                :: patscal
         integer             :: n_typ        !number of layer types    > n_layers
         real                :: l_cnt = 0.0  !number of layers in stacking section
         integer             :: SymGrpNo
         integer             :: n_cycles=0
         integer             :: n_actual       !number of unique layers
         integer             :: npar = 0       !total number of parameters to be refined
         integer             :: num_temp=0     ! maximum number of temperatures in SAN
         integer             :: init_config    ! flag to randomly initialaze SAN configuration if = 1.
         integer             :: yy             !max number of atoms per layer
         integer             :: n_seq          !number of semirandom sequences
         integer             :: num_bgrpatt
         integer             :: cheb_nump      !number of parameters in cheb polynomial
         logical             :: recrsv         !layer sequencing recursive
         logical             :: xplcit         !layer sequencing explicit
         logical             :: finite_width   !layer width in plane ab
         logical             :: inf_thick      !infinite stacking
         logical             :: trm            !trim keyword: to supress peak at origin >trim_origin
         logical             :: bgrpatt, bgrinter, bgrcheb
         character(len=4),  dimension(:,:),allocatable :: a_name    !atom name (4characters wide)>atom_l o a_name?
         character(len=20), dimension(:),  allocatable :: bfilepat  !filename for bgrpatt
         character(len=20), dimension(:),  allocatable :: bfilehkl  !filename for reflection indices and positions associated to a
                                                                    !background pattern
         integer, dimension(:),     allocatable  :: centro              !layer symmetry : none or centrosymmetric  > only_real
         integer, dimension(:),     allocatable  :: cod
         integer, dimension(:,:),   allocatable  :: a_num               !atom number  > a_number
         integer, dimension(:)  ,   allocatable  :: l_seq               !d-> array containing the explicitly defined sequence of layers.
         integer, dimension(:)  ,   allocatable  :: l_actual            !Contains the layer number that layer i is structurally identical to.
         integer, dimension(:)  ,   allocatable  :: l_n_atoms           !number of atoms in each layer
         integer, dimension(:),     allocatable  :: p                   ! refinable parameter index
         real,    dimension(:),     allocatable  :: mult                ! refinable parameter multiplicator
         real,    dimension(:,:,:), allocatable  :: a_pos               ! xyz coordinates of the atoms   >a_pos
         real,    dimension(:,:),   allocatable  :: a_B                 !debye waller factors
         real,    dimension(:,:),   allocatable  :: a_occup             !occupation factor (between 0 and 1)
         real,    dimension(:)  ,   allocatable  :: high_atom, low_atom !lower and upper atom positions
         real,    dimension(:,:),   allocatable  :: l_alpha             !l_alpha
         real,    dimension(:,:,:), allocatable  :: l_r                 !transitions vector
         real,    dimension(:,:),   allocatable  :: r_b11  , r_B22 , r_B33 , r_B12 , r_B23, r_b13
         real,    dimension(:),     allocatable  :: chebp
         real,    dimension(:),     allocatable  :: bscalpat !scale factor for bgrpatt
         real,    dimension(:)  ,   allocatable  :: list     !vector containing all the parameters that can be optimized
         real,    dimension(:)  ,   allocatable  :: Pv_refi  !vector containing the free parameters to optimize (restrictions taken into account)

         character(len=132)               :: ttl          !title
         logical,dimension(:),allocatable :: fundamental  !array which defines if each layer is fundamental or not
         logical                          :: randm, semirandm, spcfc !defines what of the explicit cases is
         integer,dimension(:),allocatable :: original     !
         integer                          :: fls, lls, otls, stls  !data of the sequence in the semirandom case

       End Type crys_2d_type

       integer, parameter ,private       :: i_data=30
       logical,            public, save  :: Err_crys=.false.
       character(len=120), public, save  :: Err_crys_mess=" "

       character(len=132), dimension(:),allocatable :: tfile     !List of lines of the input file (UPPER CASE)
       Integer               :: numberl   !number of lines in the input file
       Integer               :: np        !to count number of parameters with refinement codes = npar
       Integer, dimension(8) :: sect_indx=0 !Indices (line numbers) of the sevent sections:
                                            !1:TITLE, 2:INSTRUMENTAL, 3:STRUCTURAL, 4:LAYER
                                            !5:STACKING , 6:TRANSITIONS, 7:CALCULATION,
                                            !8:EXPERIMENTAL
      integer, parameter,public :: max_fst=40
      integer,public            :: num_fst !number of fp_studio commands
      type(crys_2d_type),              save,  public  :: crys
      type(Opt_Conditions_Type),       save,  public  :: opti
      type(LSQ_Conditions_type),       save,  public  :: cond
      type(LSQ_State_Vector_Type),     save,  public  :: Vs
      real, dimension (3),                    public  :: aver_cell, aver_ang
      integer, dimension(max_avercell),       public  :: lay_avcell          !sequence of layers to be stacked for calculating the average cell
      character(len=180), dimension(max_fst), public  :: fst_cmd             !Commands for FP_STUDIO
      character(len=100),                     public  :: filenam
      integer,                                public  :: nlay_avcell         !number of layers needed for calculating the average cell
      Real, dimension(:,:),     allocatable,  public  :: bgr_patt
      Real, dimension(:,:),     allocatable,  public  :: bgr_hkl_pos
      integer, dimension(:,:,:),allocatable,  public  :: bgr_hkl_ind
      integer, dimension(:),    allocatable,  public  :: bgr_hkl_nref
      logical,                                public  :: calculate_aberrations=.false.
      logical, save,                          public  :: aver_cellgiven=.false.,aver_spggiven=.false., aver_cellcalc=.false.
      logical, save,                          public  :: fst_given=.false.

   contains

    Subroutine Read_Bgr_patterns()
      integer :: i,j,k,i_pat,nlin,maxref
      real    :: ai,st,af
      character(len=132) :: line
      logical :: done
      done=.false.
      nlin=0; maxref=0
      do i=1,crys%num_bgrpatt
        if(len_trim(crys%bfilehkl(i)) /= 0) then
          call number_lines(trim(crys%bfilehkl(i)),nlin)
          if(nlin > maxref) maxref=nlin
          done=.true.
        end if
      end do
      if(done) then
        allocate(bgr_hkl_pos(maxref,crys%num_bgrpatt))
        allocate(bgr_hkl_ind(3,maxref,crys%num_bgrpatt))
        allocate(bgr_hkl_nref(crys%num_bgrpatt))
        bgr_hkl_nref=0
        bgr_hkl_pos=0.0
        bgr_hkl_ind=0
     end if
      do i=1,crys%num_bgrpatt
         open(newunit=i_pat,file=trim(crys%bfilepat(i)),status="old",action="read",position="rewind")
         read(unit=i_pat,fmt=*) ai,st,af
         read(unit=i_pat,fmt=*) bgr_patt(:,i)
         close(unit=i_pat)
         if(.not. done) cycle
         if(len_trim(crys%bfilehkl(i)) /= 0) then
           call number_lines(trim(crys%bfilehkl(i)),nlin)
           bgr_hkl_nref(i)=nlin-3
           open(newunit=i_pat,file=trim(crys%bfilehkl(i)),status="old",action="read",position="rewind")
           read(unit=i_pat,fmt="(a)") line
           read(unit=i_pat,fmt="(a)") line
           read(unit=i_pat,fmt="(a)") line
           do j=1, bgr_hkl_nref(i)
             read(unit=i_pat,fmt=*) bgr_hkl_ind(:,j,i),k,ai,bgr_hkl_pos(j,i)
           end do
           close(unit=i_pat)
         end if
      end do
    End Subroutine Read_Bgr_patterns

    Subroutine Set_TFile(namef,logi)
      character(len=*),  intent(in)  :: namef
      logical,           intent(out) :: logi
      ! integer :: i
      call number_lines(namef, numberl)                            ! count number of lines in file

      logi=.true.
      if (numberl == 0) then
        err_crys=.true.
        err_crys_mess="ERROR :  The file "//trim(namef)//" contains nothing"
        logi=.false.
        return
      else
        if (allocated (tfile)) deallocate (tfile)
        allocate (tfile(numberl))
        tfile=" "
        call reading_lines(namef, numberl, tfile)   ! we 'load' the file in tfile so we can close the unit
      end if
      call Set_Sections_index()
      return
    End Subroutine Set_TFile

    Subroutine read_fraction (line, real_num)

        character(len=20), intent(in)  :: line        !string of numbers to convert to real
        real              ,intent(out) :: real_num    ! string of real characters

        integer                                       ::  j , ierr
        real                                          :: numerator, denominator

        j=0
        j=index(line,"/")
        if(j == 0) then                  !it is a real number
          read ( line, fmt=*,iostat=ierr) real_num
          if(ierr /= 0) then
            ERR_String= .true.
            ERR_String_Mess=" Error reading a real number! "
            return
          end if
        else         !it was expressed as a ratio. delete the slash, '/'.
          read(unit =line(1:j-1),fmt = *,iostat=ierr) numerator
          if(ierr /= 0) then
            ERR_String= .true.
            ERR_String_Mess=" Error reading the numerator of a fraction! -> "//line(1:j-1)
            return
          end if
          read(unit =line(j+1:),fmt = *,iostat=ierr) denominator  ! now contains two arguments, numerator and denominator.
          if(ierr /= 0) then
            ERR_String= .true.
            ERR_String_Mess=" Error reading the denominator of a fraction! -> "//trim(line(j+1:))
            return
          end if
          real_num = numerator/ denominator
        end if
    End subroutine read_fraction

    Subroutine Set_Sections_index()
      integer  :: k, l
      character(len=132) :: txt

      k=0
      !r=0   !counts n_actual
      l=0   !counts n_layer
      global:   do
        k=k+1; if(k > numberl) exit
        txt=u_case(adjustl(tfile(k)))
        if(txt(1:5) == "TITLE") then
          sect_indx(1) = k
          do
            k=k+1; if(k > numberl) exit
            txt=u_case(adjustl(tfile(k)))
            if(txt(1:12) == "INSTRUMENTAL") then
              sect_indx(2) = k
              do
                k=k+1; if(k > numberl) exit
                txt=u_case(adjustl(tfile(k)))
                if(txt(1:10) == "STRUCTURAL") then
                  sect_indx(3) = k
                  do
                    k=k+1; if(k > numberl) exit
                    txt=u_case(adjustl(tfile(k)))
                    if(txt(1:5) == "LAYER") then
                       sect_indx(4) = k
                       !n_actual=r+1
                       n_layers=l+1
                       do
                         k=k+1; if(k > numberl) exit
                         txt=u_case(adjustl(tfile(k)))
                         if(txt(1:8) == "STACKING") then
                            sect_indx(5) = k
                            do
                              k=k+1; if(k > numberl) exit
                              txt=u_case(adjustl(tfile(k)))
                              if(txt(1:11) == "TRANSITIONS") then
                                 sect_indx(6) = k
                                 do
                                   k=k+1; if(k > numberl) exit
                                   txt=u_case(adjustl(tfile(k)))
                                   if(txt(1:11) == "CALCULATION") then
                                      sect_indx(7) = k
                                      do
                                        k=k+1; if(k > numberl) exit
                                        txt=u_case(adjustl(tfile(k)))
                                        if(txt(1:12) == "EXPERIMENTAL") then
                                           sect_indx(8) = k
                                           exit global
                                     end if
                                   end do
                                 end if
                               end do
                             end if
                           end do
                         end if
                       end do
                    end if
                  end do
                end if
              end do
            end if
          end do
        end if
      end do global


    End Subroutine Set_Sections_index


 !1:TITLE, 2:INSTRUMENTAL, 3:STRUCTURAL, 4:LAYER
 !5:STACKING , 6:TRANSITIONS, 7:CALCULATION,
 !8:EXPERIMENTAL

    Subroutine Read_TITLE(logi)
      logical, intent(out) :: logi

      integer :: i,i1,i2 !,ier,k
      character(len=132) :: txt=" "

      logical :: ok_title

      logi=.true.
      i1=sect_indx(1)
      i2=sect_indx(2)-1

      i=i1

      ok_title=.false.

      do
        i=i+1
        if(i > i2) exit
        txt=adjustl(tfile(i))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        crys%ttl = txt
        ok_title=.true.
      end do

    End Subroutine Read_TITLE

    Subroutine Read_INSTRUMENTAL(logi)
      logical, intent(out) :: logi

      integer :: i,i1,i2,k, ier !,l
      character(len=132) :: txt
      character(len=25)  :: keyv
      logical :: ok_rad, ok_wave, ok_uvw, ok_abe
      !real  :: ab

      logi=.true.
      i1=sect_indx(2)
      i2=sect_indx(3)-1

      i=i1

      ok_rad=.false.; ok_wave=.false.; ok_uvw=.false. ; ok_abe=.true.
      !Initialise codes to zero as well as values
       vs%code=0; vs%pv=0.0; vs%spv=0.0; vs%nampar=" "; vs%dpv=0.0; vs%np=0

      np=0
      do

        i=i+1; if(i > i2) exit
        txt=u_case(adjustl(tfile(i)))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        keyv=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(keyv)

          Case("RADIATION")

                if (index(txt , 'X-RAY')/=0 ) then
                  crys%rad_type = 0
                else if (index(txt , 'NEUTRON')/=0) then
                  crys%rad_type = 1
                else if (index(txt , 'ELECTRON')/=0) then
                  crys%rad_type = 2
                else
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading radiation type"
                  logi = .false.
                  return
                end if
                ok_rad=.true.

          Case("WAVELENGTH")

            read(unit=txt,fmt=*, iostat=ier) crys%lambda, crys%lambda2, crys%ratio
              if(ier /= 0 ) then
                read(unit=txt,fmt=*, iostat=ier) crys%lambda, crys%lambda2
                if(ier /= 0 ) then
                  read(unit=txt,fmt=*, iostat=ier) crys%lambda
                  if(ier/=0) then
                     Err_crys=.true.
                     Err_crys_mess="ERROR reading wavelength instruction"
                     logi=.false.
                     return
                  end if
                else
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading ratio"
                  logi=.false.
                  return
                end if
              end if
            ok_wave=.true.

          Case("ABERRATIONS")

             read(unit=txt,fmt=*, iostat=ier) crys%zero_shift, crys%sycos, crys%sysin
             val_glb(2:4)=[crys%zero_shift, crys%sycos, crys%sysin]
             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading instrumental aberrations"
               logi=.false.
               return
             end if
             if(sum(abs(val_glb(2:4))) > 0.0001)  calculate_aberrations=.true.
             i=i+1
             txt=adjustl(tfile(i))
             !Reading refinement codes
             read(unit=txt,fmt=*, iostat=ier) ref_glb(2:4)   !codes for zero_shift,sycos,sysin
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the instrumental aberrations"
               logi=.false.
               return
             end if
             ok_abe=.true.
             if(sum(abs(ref_glb(2:4))) > 0.1)  calculate_aberrations=.true.


          Case("PSEUDO-VOIGT")
             crys%broad=ps_vgt
             read(unit=txt,fmt=*, iostat=ier) crys%p_u, crys%p_v, crys%p_w, crys%p_x, crys%p_dg, crys%p_dl
             val_glb(5:10)=[crys%p_u, crys%p_v, crys%p_w, crys%p_x, crys%p_dg, crys%p_dl]
             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading Pseudo-Voigt instruction"
               logi=.false.
               return
             end if

             if(index(txt,"TRIM") /= 0)   crys%trm=.true.

             i=i+1
             txt=adjustl(tfile(i))
             !Reading refinement codes
             read(unit=txt,fmt=*, iostat=ier) ref_glb(5:10) !codes for U,V,W,X,Dg,DL
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the Pseudo-Voigt instruction"
               logi=.false.
               return
             end if
             ok_uvw=.true.

          Case("GAUSSIAN")

             crys%broad= pv_gss

             read(unit=txt,fmt=*, iostat=ier) crys%p_u, crys%p_v, crys%p_w, crys%p_dg
             crys%p_dl=100000
             val_glb(5:10)=[crys%p_u, crys%p_v, crys%p_w, crys%p_x, crys%p_dg, crys%p_dl]
             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading Gaussian/Lorentzian instruction"
               logi=.false.
               return
             end if

             if(index(txt,"TRIM") /= 0) crys%trm=.true.

             i=i+1
             txt=adjustl(tfile(i))
             !Reading refinement codes
             read(unit=txt,fmt=*, iostat=ier) ref_glb(5:9) !codes for U,V,W,X,Dg
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the Gaussian/Lorentzian instruction"
               logi=.false.
               return
             end if
             ok_uvw=.true.

          Case("LORENTZIAN")

             crys%broad= pv_lrn

             !read(unit=txt,fmt=*, iostat=ier) crys%p_u, crys%p_v, crys%p_w, crys%p_dl
             read(unit=txt,fmt=*, iostat=ier) ref_glb(5:10)
             crys%p_dg=100000
             val_glb(5:10)=[crys%p_u, crys%p_v, crys%p_w, crys%p_x, crys%p_dg, crys%p_dl]
             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading Gaussian/Lorentzian instruction"
               logi=.false.
               return
             end if

             if(index(txt,"TRIM") /= 0) crys%trm=.true.

             i=i+1
             txt=adjustl(tfile(i))
             !Reading refinement codes
             read(unit=txt,fmt=*, iostat=ier) ref_glb(5:8),ref_glb(10)
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the Gaussian/Lorentzian instruction"
               logi=.false.
               return
             end if
            ok_uvw=.true.

          Case Default

             cycle

        End Select
      end do

      if(ok_rad .and. ok_wave .and. ok_uvw .and. ok_abe) then
        return
      else
        Err_crys=.true.
        Err_crys_mess="ERROR: Radiation type, wavelength or UVW parameters missing!"
        logi=.false.
      end if

    End Subroutine Read_INSTRUMENTAL

    Subroutine Read_STRUCTURAL(logi)
      logical, intent(out) :: logi

      integer :: i,i1,i2,j,k, ier !,l
      character(len=132) :: txt
      character(len=25)  :: keyv
      logical :: ok_cell, ok_symm, ok_nlayers, ok_lwidth
      !real  :: ab

      logi=.true.
      i1=sect_indx(3)
      i2=sect_indx(4)-1
      num_fst=0
      i=i1

      crys%finite_width=.false.

      ok_cell=.false.; ok_symm=.false.; ok_nlayers=.false. ; ok_lwidth=.false.

      do

        i=i+1; if(i > i2) exit
        txt=u_case(adjustl(tfile(i)))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        keyv=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(keyv)

          Case("SPGR")
                txt=adjustl(tfile(i))
                spgsymb=adjustl(txt(5:))  !read symbol of space group
                j=index(spgsymb,"!")
                if(j /= 0) spgsymb=trim(spgsymb(1:j-1))
                call Set_SpaceGroup(spgsymb, SpG)
                if (Err_Symm)then
                  Err_crys=.true.
                  Err_crys_mess="Error in space group symbol: "//trim(spgsymb)
                  logi = .false.
                  return
                end if
                aver_spggiven=.true.

          Case("AVERCELL")
                read(unit=txt,fmt=*,iostat=ier) aver_cell, aver_ang   !read average cell parameters
                if (ier /= 0)then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading average cell parameters"
                  logi = .false.
                  return
                end if
                aver_cellgiven=.true.
                call Set_Crystal_Cell(aver_cell, aver_ang, aCell)

         Case("FST_CMD")
                num_fst=num_fst+1
                if(index(txt,"SEQ") /= 0) fst_given=.true.
                txt=adjustl(tfile(i))
                fst_cmd(num_fst)= adjustl(txt(8:)) !store the fst command

          Case("CALCAVERCELL")
                write(unit=*,fmt="(a)") " => An average cell will be calculated at the end of the refinement."
                read(unit=txt, fmt=*, iostat=ier) nlay_avcell          !read number of layers to be read
                if (ier /= 0)then
                    Err_crys=.true.
                    Err_crys_mess="ERROR reading the number of transition vectors for calculating the average cell"
                    logi = .false.
                    return
                end if
                !write(unit=*,fmt="(a)") "Number of vectors to be read for calculating the average cell: "
                !write(unit=*,fmt=*)  nvect_avcell
                do
                    i=i+1
                    txt=adjustl(tfile(i))
                    if(txt(1:1) == "!")cycle
                    exit
                end do
                read (unit=txt, fmt=*, iostat=ier) lay_avcell(1:nlay_avcell)   !read layers for calculating the average cell
                lay_avcell(nlay_avcell+1)=lay_avcell(1)                        !to find the last vector needed for the average cell,
                                                                               !we suppose that the last layer is the same as the first one
                if (ier /= 0)then
                    Err_crys=.true.
                    Err_crys_mess="ERROR reading the transition vectors for calculating the average cell"
                    logi = .false.
                    return
                end if
                !write(unit=*,fmt="(a)") "The vectors needed for the average cell are"
                !do j = 1,nvect_avcell
                !        write(unit=*,fmt="(a)") vect_avcell(j)
                !end do
                aver_cellcalc=.true.



          Case("CELL")
                read(unit=txt,fmt=*,iostat=ier) crys%cell_a, crys%cell_b, crys%cell_c, crys%cell_gamma   !read cell parameters
                val_glb(11:14)=[crys%cell_a, crys%cell_b, crys%cell_c, crys%cell_gamma ]
                if (ier /= 0)then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading cell parameters"
                  logi = .false.
                  return
                end if
                crys%cell_gamma = crys%cell_gamma * deg2rad

             i=i+1
             txt=adjustl(tfile(i))
             !Reading refinement codes
             read(unit=txt,fmt=*, iostat=ier) ref_glb(11:14)
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the cell parameters"
               logi=.false.
               return
             end if
             ok_cell=.true.

          Case("SYMM")

             j=index(txt, " ")
             if( j < 2) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the SYMM instruction"
               logi=.false.
               return
             end if
             crys%sym=txt(1:j-1)

             read(unit=txt(j:) ,fmt=*, iostat = ier) crys%tolerance
             if(ier /= 0 ) crys%tolerance = 1

             crys%tolerance = crys%tolerance * eps2   ! convert from a percentage

             if (crys%sym == "-1") then
               Crys%Symgrpno = 1
             else if ( crys%sym == "2/M(1)") then
               Crys%Symgrpno = 2
             else if ( crys%sym == " 2/M(2)") then
               Crys%Symgrpno = 3
             else if ( crys%sym == "MMM")  then
               Crys%Symgrpno = 4
             else if ( crys%sym == "-3") then
               Crys%Symgrpno = 5
             else if ( crys%sym == "-3M") then
               Crys%Symgrpno = 6
             else if ( crys%sym == "4/M")  then
               Crys%Symgrpno = 7
             else if ( crys%sym == "4/MMM")  then
               Crys%Symgrpno = 8
             elseif ( crys%sym == "6/M") then
               Crys%Symgrpno = 9
             else if ( crys%sym == "6/MMM")  then
               Crys%Symgrpno = 10
             else if ( crys%sym == "AXIAL")  then
               Crys%Symgrpno = 11
             else if ( crys%sym == "UNKNOWN")  then
               Crys%Symgrpno = -1
             end if

             ok_symm=.true.

          Case("NLAYERS")

             read(unit=txt,fmt=*, iostat=ier)  crys%n_typ
             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading number of layer types"
               logi=.false.
               return
             end if
             if(crys%n_typ==0) then
              Err_crys=.true.
              Err_crys_mess="ERROR reading number of layer types. It can not be zero."
              logi=.false.
              return
             end if
             ok_nlayers=.true.

          Case ("LWIDTH")

             if (index(txt, 'INFINITE') /= 0) then
               ok_lwidth=.true.
               cycle
             else if (index(txt, 'INFINITE') == 0) then
               read(unit=txt,fmt=*, iostat=ier)  crys%layer_a, crys%layer_b
               val_glb(15:16)=[crys%layer_a, crys%layer_b]
               crys%finite_width=.true.

               if (crys%layer_a == 0 .or. crys%layer_b == 0) then
                 Err_crys=.true.
                 Err_crys_mess="ERROR reading layer width parameters: width along a or b cannot be zero"
                 logi=.false.
                 return
               end if
               if(ier /= 0 ) then
                 Err_crys=.true.
                 Err_crys_mess="ERROR reading layer width parameters"
                 logi=.false.
                 return
               end if

               i=i+1
               txt=adjustl(tfile(i))
              ! Reading refinement codes
               !k=index(txt,"(")
               !l=index(txt,")")

               read(unit=txt,fmt=*, iostat=ier) ref_glb(15:16)  !Codes for Wa,Wb

               if(ier /= 0) then
                 Err_crys=.true.
                 Err_crys_mess="ERROR reading the refinement codes of layer width parameters"
                 logi=.false.
                 return
               end if
               ok_lwidth=.true.
             else
                Err_crys=.true.
                Err_crys_mess="ERROR: No parameter in the Layer width instruction: this must be given!"
                logi=.false.
                return
             end if



          Case Default
             Write(unit=*,fmt="(a)") " => WARNING! Unknown keyword: "//trim(keyv)
             cycle

          End Select
      end do

      if(ok_cell .and. ok_symm .and. ok_nlayers .and. ok_lwidth) then
        return
      else
        Err_crys=.true.
        Err_crys_mess="ERROR: cell parameters, Laue symmetry, number of layer types, or layers width missing!"
        logi=.false.
      end if

    End Subroutine Read_STRUCTURAL

    Subroutine Read_LAYER(logi)
      logical, intent(out) :: logi

      integer :: i,i1,i2 !,i3
      integer :: k, ier, a, a1, a2, r, tmp, j, m, l, nitem
      character(len=20), dimension(30) :: citem
      character(len=132) :: txt, layer
      character(len=25)  :: keyv

      logical :: ok_lsym, ok_atom
      integer,          dimension(:),     allocatable  :: d !counts n\BA of atoms in unique layer
      !real  :: ab

      if (allocated (crys%fundamental)) deallocate(crys%fundamental)
      allocate(crys%fundamental(crys%n_typ))
      crys%fundamental=.false.

      if (allocated (fundamental)) deallocate(fundamental)
      allocate(fundamental(crys%n_typ))
      fundamental=.false.

      if (allocated (crys%original)) deallocate(crys%original)
      allocate(crys%original(crys%n_typ))
      crys%original=0

      if (allocated (original)) deallocate(original)
      allocate(original(crys%n_typ))
      original=0

      if (allocated (d)) deallocate(d)
      allocate(d(max_a))
      d=0

      logi=.true.
      i1=sect_indx(4)
      i2=sect_indx(5)-1

      i=i1-1

      ok_lsym=.false.; ok_atom=.false.

       a = 0      ! counts l_n_atoms
       r = 0      ! counts n_actual
       a1= 0      ! reads layer number
       a2= 0      ! reads layer number
       !yy=0       ! counts n\BA of atoms in layer type
       j=0
       l=0
       m=0
      do
        i=i+1; if(i > i2) exit
        txt=u_case(adjustl(tfile(i)))

        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        keyv=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(keyv)

          Case ("LAYER")
            if (index(txt,'=') /= 0) then        ! search for '=' sign
              j = j+1
              read (unit = txt, fmt =*, iostat = ier) a1, layer, a2
              crys%fundamental(j)=.false.
              crys%original(j)= a2
              if(ier /= 0) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading layer equivalencies"
                logi=.false.
                return
              end if
              if (a2 >= a1 .OR. a2 < 1) then
                write (*,*) "Layer", a1, " cannot be equal to layer" ,a2, "."
                logi = .false.
                return
              else
                crys%l_actual(j)  = crys%l_actual(a2)
              end if
            else
              j=j+1
              r=r+1
              crys%l_actual(j) = r
              crys%fundamental(j)=.true.
            end if
            crys%n_actual = r
            n_layers=j

          Case("LSYM")
            if(INDEX(txt, 'CENTROSYMMETRIC') == 1) then
              crys%centro(r) = CENTRO
            else if (INDEX(txt, 'NONE') == 1) then
              crys%centro(r) = NONE
            else
                Err_crys=.true.
                Err_crys_mess="ERROR: error reading layer symmetry"
                logi=.false.
                return
            end if

                ok_lsym=.true.

          Case ("ATOM")
            ok_atom=.false.
            a=a+1

            d(r)=d(r)+1

            call getword(txt, citem, nitem)
            if(nitem<7) then
              write(unit=*,fmt="(3a,i3)")" => ERROR reading atomic parameters. Parameter missing in atom ", &
                                             trim(citem(2))," from layer ", r
              logi=.false.
              return
            end if

            do m=1,3

              call read_fraction(citem(2+m), crys%a_pos(m, d(r),r))
              val_atom(m, d(r),r)=crys%a_pos(m, d(r),r)
              if(ERR_String) then
                write(unit=*,fmt="(a)") trim(ERR_String_Mess)
                logi=.false.
                return
              end if
            end do
            crys%a_name(d(r),r)=citem(1)
            read(unit=citem(2),fmt=*,iostat=ier) crys%a_num (d(r),r)
            read(unit=citem(6),fmt=*,iostat=ier) crys%a_B (d(r),r)
            val_atom(4,d(r),r)=crys%a_B (d(r),r)

            call read_fraction(citem(7), crys%a_occup (d(r),r))
            !read(unit=citem(7),fmt=*,iostat=ier) crys%a_occup (d(r),r)
            val_atom(5,d(r),r)=crys%a_occup (d(r),r)
             if(crys%a_occup(d(r),r)<0 .or. crys%a_occup(d(r),r)>1) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading atomic parameters. Occupation must be between 0 and 1"
               logi=.false.
               return
             end if
            if(ier /= 0) then
              Err_crys=.true.
              Err_crys_mess="ERROR reading atomic parameters"
              logi=.false.
              return
            end if

            tmp = crys%a_pos(3,a,r)                          !to asign lower and upper positions ****
            IF(tmp > crys%high_atom(r)) crys%high_atom(r) = tmp
            IF(tmp < crys%low_atom(r))   crys%low_atom(r) = tmp

            i=i+1
            txt=adjustl(tfile(i))

            !reading refinement codes
            read(unit=txt,fmt=*, iostat=ier) ref_atom(1:5,d(r),r)
            if(ier /= 0) then
                   Err_crys=.true.
                   Err_crys_mess="ERROR reading atomic parameters refinement codes"
                   logi=.false.
                   return
            end if
            crys%l_n_atoms(r) = d(r)
            ok_atom=.true.
            crys%n_actual=r
           !if (d(r) > yy) yy=d(r)

          Case Default

             cycle

          End Select

       end do

       if(abs(crys%n_typ - n_layers) >= 1.0) then                      !NIK
            Err_crys=.true.                                            !
            Err_crys_mess="ERROR with the number of layer types."      !
            logi=.false.                                               !
            return                                                     !
       end if                                                          !

       ! make sure we have genuine lowest and highest atoms in each layer.
       DO  i = 1, crys%n_actual
         IF(l_symmetry(i) == centro) THEN
           IF(-crys%low_atom(i) > crys%high_atom(i)) THEN
              crys%high_atom(i) = -crys%low_atom(i)
           ELSE IF(-crys%low_atom(i) < crys%high_atom(i)) THEN
                   crys%low_atom(i) = -crys%high_atom(i)
           END IF
         END IF
       END DO
       crys%yy=maxval(crys%l_n_atoms)
       !  write( *,*) "ok_lsym .and. ok_atom" , ok_lsym , ok_atom
       if(ok_lsym .and. ok_atom) then
         return
       else
         Err_crys=.true.
         Err_crys_mess="ERROR: layer symmetry or layer atomic parameters missing!"
         logi=.false.
       end if

    End Subroutine Read_LAYER

    Subroutine Read_STACKING(logi)
    logical, intent(out) :: logi

      integer :: i,i1,i2,k, ier,j, j1, j2, l1, l2, r !,l
      integer,dimension(132)                         :: Inte    !  Vector of integer numbers
      integer                                        :: n_int  !number of integers per line
      real(kind=sp),      dimension(132)             :: rel !  Vector of integer numbers
      character(len=132) :: txt
      character(len=25)  :: keyv, seq
      logical :: ok_explicit, ok_recursive
      !real  :: ab

      logi=.true.
      i1=sect_indx(5)
      i2=sect_indx(6)-1

      i=i1

      ok_explicit=.false.; ok_recursive=.false.

      crys%randm=.false.; crys%semirandm=.false.; crys%spcfc=.false.

      do

        i=i+1; if(i > i2) exit
        txt=u_case(adjustl(tfile(i)))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        keyv=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(keyv)
          Case ("EXPLICIT")
            crys%xplcit = .true.
            crys%l_cnt = 0.0
            i=i+1
            txt=adjustl(tfile(i))
            txt=u_case(adjustl(tfile(i)))
            if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) then
              i=i+1
              txt=u_case(adjustl(tfile(i)))
            end if
            if (index(txt , 'SPECIFIC')/=0 ) then
              crys%spcfc=.true.
              crys%inf_thick = .false.
              read (unit = txt, fmt = *, iostat = ier) lstype
              do
                i=i+1
                if (i>i2) exit
                txt=adjustl(tfile(i))
                call getnum(txt, rel,Inte, n_int)    ! to know number of integers per line
                l1=crys%l_cnt+1
                l2=crys%l_cnt+n_int
                read(unit=txt,fmt=*, iostat=ier)(crys%l_seq(r), r=l1,l2)  ! to read stacking sequence

                r= r-1
                crys%l_cnt = r
              end do

            else if (index(txt , 'SEMIRANDOM')/=0) then
              crys%semirandm=.true.
              read (unit = txt, fmt = *, iostat = ier) lstype, crys%l_cnt
                if   (crys%l_cnt==0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR :  Number of layers is missing"
                  logi = .false.
                  return
                end if
                val_glb(17) = crys%l_cnt
                ref_glb(17) = 0.0
                crys%inf_thick = .false.
                rndm = .true.
                crys%n_seq=0
                i=i+1
                txt=u_case(adjustl(tfile(i)))
                if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
                if (index(txt , 'SEQ')==1) then
                  read (unit = txt, fmt = *, iostat = ier) seq, j1, j2, l1, l2
                  crys%fls=j1
                  crys%lls=j2
                  crys%otls=l1
                  crys%stls=l2
                  crys%n_seq=crys%n_seq+1
                  j=j1
                  do
                    crys%l_seq(j)= l1
                    if (j==j2)then
                      exit
                    end if
                    j=j+1
                    crys%l_seq(j) = l2
                    if (j==j2)then
                      exit
                    end if
                    j=j+1
                  end do
                else
                  !Err_crys=.true.
                  !Err_crys_mess="ERROR : Layer explicit sequences missing"
                  write(*,"(a)") " => ERROR : Layer explicit sequences missing"
                  logi = .false.
                  return
                end if

            else if (index(txt , 'RANDOM')/=0) then
              crys%randm=.true.
              read (unit = txt, fmt = *, iostat = ier) lstype, crys%l_cnt

                if   (crys%l_cnt==0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR :  Number of layers is missing"
                  logi = .false.
                  return
                end if
                crys%inf_thick = .false.
                rndm = .true.
                val_glb(17) = crys%l_cnt
                ref_glb(17) = 0.0

            else
              Err_crys=.true.
              Err_crys_mess="ERROR reading explicit layer sequence"
              logi = .false.
              return
            end if
            ok_explicit=.true.

          Case("RECURSIVE")
            crys%recrsv = .true.
            i=i+1
            txt=u_case(adjustl(tfile(i)))
            if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) then
              i=i+1
              txt=u_case(adjustl(tfile(i)))
            end if
            if (index(txt , 'INFINITE')==1) then
              crys%inf_thick = .true.
            else
              read(unit=txt ,fmt= * , iostat = ier)crys%l_cnt
              val_glb(17) = crys%l_cnt
              if (crys%l_cnt >= 1023.0) then
                crys%inf_thick = .true.
              else
                crys%inf_thick = .false.
              end if
              i=i+1
              txt=adjustl(tfile(i))
              !Reading refinement codes
              read(unit=txt,fmt=*, iostat=ier) ref_glb(17)
              if(ier /= 0) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading the refinement codes of the number of layers in the crystal"
                logi=.false.
                return
              end if

            end if
            ok_recursive=.true.

          Case Default
            cycle
          End Select
     end do

     if(ok_recursive  .or. ok_explicit) then
        return
     else
        Err_crys=.true.
        Err_crys_mess="ERROR: layer stacking parameters missing!"
        logi=.false.
     end if

    End Subroutine Read_STACKING

    Subroutine Read_TRANSITIONS(logi)
      logical, intent(out) :: logi
      integer :: i,i1,i2, k, ier, j, l, nitem, m
      real, parameter :: eps =1.0E-4
      real  :: suma !,ab
      character(len=20), dimension(30) :: citem
      character(len=132) :: txt
      character(len=25)  :: keyv
      logical :: ok_lt

      logi=.true.
      ok_lt=.false.
      i1=sect_indx(6)
      i2=sect_indx(7)-1

      i=i1+1
      l=1
      j=0
      m=0

        do
          if(i > i2) exit
          txt=u_case(adjustl(tfile(i)))
          if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) then
           i=i+1
           cycle
          end if
          k=index(txt," ")
          keyv=txt(1:k-1)
          txt=adjustl(txt(k+1:))


        SELECT CASE (keyv)

          CASE ("LT")
            if (j==n_layers) then
                l=l+1
                j=0
            end if

            j=j+1

            call getword(txt, citem, nitem)
            do m=1, 3
              call read_fraction(citem(1+m), crys%l_r (m,j,l))
              val_trans(m+1,j,l)=crys%l_r (m,j,l)
              if(ERR_String) then
                write(unit=*,fmt="(a)") " => ERROR reading transition vectors coordinates", trim(ERR_String_Mess)
                logi=.false.
                return
              end if
            end do

            call read_fraction(citem(1), crys%l_alpha (j,l))
           !read(unit=citem(1),fmt=*,iostat=ier) crys%l_alpha (j,l)
            val_trans(1,j,l) = crys%l_alpha (j,l)
            if(ERR_String) then
                write(unit=*,fmt="(a)") " => ERROR reading layer probabilities", trim(ERR_String_Mess)
                logi=.false.
                return
            end if
           !if(ier /= 0) then
           !       Err_crys=.true.
           !       Err_crys_mess="ERROR reading layer probabilities"
           !       logi=.false.
           !       return
           !end if

            i=i+1
            txt=adjustl(tfile(i))

            read (unit = txt, fmt =*, iostat = ier) ref_trans(1:4, j, l)


            i=i+1

          CASE ("FW")
            read (unit = txt, fmt =*, iostat = ier) crys%r_b11(j,l) , crys%r_b22(j,l) , crys%r_b33(j,l) , &
                                                    crys%r_b12(j,l) , crys%r_b13(j,l) , crys%r_b23(j,l)
            val_trans(5:10, j,l)=[crys%r_b11(j,l) , crys%r_b22(j,l) , crys%r_b33(j,l) , &
                                 crys%r_b12(j,l) , crys%r_b13(j,l) , crys%r_b23(j,l)]

            i=i+1
            txt=adjustl(tfile(i))

            read(unit=txt,fmt=*, iostat=ier) ref_trans(5:10, j, l)

            if(ier /= 0) then
              Err_crys=.true.
              Err_crys_mess="ERROR reading fats waller parameters"
              logi=.false.
              return
            end if
            i=i+1

          CASE ("TABLE")
            if      (u_case(trim(txt)) == "A") then
                l = 1
            else if (u_case(trim(txt)) == "B") then
                l = 2
            else if (u_case(trim(txt)) == "C") then
                l = 3
            else if (u_case(trim(txt)) == "ALPHA") then
                l = 4
            else if (u_case(trim(txt)) == "FW11") then
                l = 5
            else if (u_case(trim(txt)) == "FW22") then
                l = 6
            else if (u_case(trim(txt)) == "FW33") then
                l = 7
            else if (u_case(trim(txt)) == "FW12") then
                l = 8
            else if (u_case(trim(txt)) == "FW13") then
                l = 9
            else if (u_case(trim(txt)) == "FW23") then
                l = 10
            end if
            table(l) = .true.
            m=1 !counts "FROM"-layer
            do while (m <= n_layers)
                i=i+1
                txt=adjustl(tfile(i))
                !check if there is a comment ! {} or blank field
                if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) then
                    !do nothink
                else
                    call getword(txt, citem, nitem)
                    !check the number of elements within the row
                    if (nitem /= n_layers) then
                        Err_crys=.true.
                        Err_crys_mess = "ERROR reading table of probabilities. Number of columns is not equal to n_layers."
                        logi=.false.
                        return
                    end if
                    !read elements within the row
                    do j = 1, n_layers   !counts "TO"-layer
                        if (l == 4) then
                            call read_fraction(citem(j), crys%l_alpha(j,m))
                            val_trans(1,j,m) = crys%l_alpha(j,m)
                        else if (l < 4) then
                            call read_fraction(citem(j), crys%l_r (l,j,m))
                            val_trans(l+1,j,m) = crys%l_r(l,j,m)
                        else if (l == 5) then
                            call read_fraction(citem(j), crys%r_b11(j,m))
                            val_trans(l,j,m) = crys%r_b11(j,m)
                        else if (l == 6) then
                            call read_fraction(citem(j), crys%r_b22(j,m))
                            val_trans(l,j,m) = crys%r_b22(j,m)
                        else if (l == 7) then
                            call read_fraction(citem(j), crys%r_b33(j,m))
                            val_trans(l,j,m) = crys%r_b33(j,m)
                        else if (l == 8) then
                            call read_fraction(citem(j), crys%r_b12(j,m))
                            val_trans(l,j,m) = crys%r_b12(j,m)
                        else if (l == 9) then
                            call read_fraction(citem(j), crys%r_b13(j,m))
                            val_trans(l,j,m) = crys%r_b13(j,m)
                        else if (l == 10) then
                            call read_fraction(citem(j), crys%r_b23(j,m))
                            val_trans(l,j,m) = crys%r_b23(j,m)
                        end if
                    end do
                    m=m+1
                end if
            end do !do while(m_t <= n_layers)
            i=i+1
            l=1
            j=0
          !end of CASE ("TABLE")

          CASE DEFAULT
            cycle
          end select


        end do


        DO  l = 1, crys%n_typ        !check if l_apha sums one for each layer
          suma = 0.0
          DO  j = 1, crys%n_typ
            suma = suma + crys%l_alpha(j,l)
          END DO

          IF(ABS(suma - 1.0) > eps) THEN
            write(op,"(a,i3,a)") " => Stacking probabilities from LAYER ",l," do not sum  1."
            write(op,"(a,f10.5)")" => The sum is ", suma
            logi = .false.
          end if
        END DO
        ok_lt=.true.
        if(ok_lt) then
          return
        else
          Err_crys=.true.
          Err_crys_mess="ERROR: layer transition parameters missing!"
          logi=.false.
        end if

    End Subroutine Read_TRANSITIONS

    Subroutine Read_CALCULATION(logi)
    logical, intent(out) :: logi

      integer :: i,i1,i2,k, ier,m !,n,j
      character(len=132) :: txt
      character(len=25)  :: keyv
      logical :: ok_sim, ok_opt, ok_eps, ok_acc, ok_iout, ok_mxfun, ok_lma, ok_boxp, ok_maxfun, &
                 ok_range, ok_corrmax, ok_tol, ok_nprint, ok_sadp, ok_streak

      logi=.true.
      i1=sect_indx(7)
      if (sect_indx(8) == 0) then
        i2= numberl
      else
        i2=sect_indx(8)-1 !if simulation there is no experimental section
      end if

      i=i1


      ok_sim=.false.; ok_opt=.false.; ok_eps=.false.; ok_acc=.false.; ok_iout=.false.; ok_mxfun=.false.
      ok_lma=.false.; ok_boxp=.false.; ok_maxfun=.false.; ok_corrmax=.false.; ok_tol=.false. ; ok_nprint=.false.
      ok_range=.false.; ok_sadp=.false.


      do
        i=i+1; if(i > i2) exit
        txt=u_case(adjustl(tfile(i)))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle

        k=index(txt," ")
        keyv=txt(1:k-1)

        Select Case(keyv)

          Case ('REPLACE_FILES')  !Tune the opening output files to replace previous files
             replace_files = .true.

          Case ('SIMULATION')
            opt = 0
            ok_sim = .true.

          Case ("STREAK")
            funct_num = 1
            num_streak = 1
            allocate(h_streak(num_streak))
            allocate(k_streak(num_streak))
            allocate(l0_streak(num_streak))
            allocate(l1_streak(num_streak))
            allocate(dl_streak(num_streak))
            allocate(streak_flags(num_streak + 1))
            streak_flags(1) = 1
            read (unit=txt(k:), fmt=*, iostat=ier) adapt_quad, h_streak(1), k_streak(1), l0_streak(1), l1_streak(1), dl_streak(1) !there must be some variables
            if (ier /= 0 ) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading STREAK information"
                logi=.false.
                return
            end if
            n_high = nint((l1_streak(1)-l0_streak(1))/dl_streak(1)+1.0_dp)
            streak_flags(2) = n_high
            ok_streak =.true.

          Case ("UNBROADEN")
            unbroaden = .true.

          Case ("LOGARITHM")
            logarithm = .true.

          Case ("POWDER")
            k=index(txt," ")
            txt=adjustl(txt(k+1:))
            m=index(txt,"BCKG_LEVEL")
            if(m /= 0) then
              read(unit=txt(m+10:),fmt=*,iostat=ier) crys%bckg_level
              if(ier /= 0 ) crys%bckg_level=100.0
            end if
            m=index(txt,"SCALEF")
            if( m /= 0) then
              read(unit=txt(m+6:),fmt=*,iostat=ier) crys%patscal
              if(ier /= 0 ) then
                 crys%patscal=1.0
              end if
              txt=txt(1:m-1)
            else
              crys%patscal=1.0
            end if

            funct_num = 3
              read (unit = txt, fmt = *, iostat=ier) th2_min, th2_max, d_theta
                if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading 2theta ranges"
                  logi=.false.
                  return
                end if
                IF(th2_min < zero) THEN
                  WRITE(op,"(a)") ' => ERROR: 2theta min is negative.'
                  logi = .false.
                  RETURN
                END IF
                IF(th2_min >= one_eighty) THEN
                  WRITE(op,"(a)") ' => ERROR: 2theta min exceeds 180 degrees.'
                  logi = .false.
                  RETURN
                END IF
                IF(th2_max <= zero) THEN
                  WRITE(op,"(a)") ' => ERROR: Negative (or zero) value for 2theta max.'
                  logi = .false.
                  RETURN
                END IF
                IF(th2_max > one_eighty) THEN
                  WRITE(op,"(a)") ' => ERROR: 2theta max exceeds 180 degrees.'
                  logi = .false.
                  RETURN
                END IF
                ! if th2_max = 180, reduce it slightly to keep us out of trouble
                IF(th2_max == one_eighty) th2_max = th2_max - eps4
                IF(d_theta <= zero) THEN
                  WRITE(op,"(a)") ' => ERROR: Negative (or zero) value for 2theta increment.'
                  logi = .false.
                  RETURN
                END IF
                IF(th2_min >= th2_max) THEN
                  WRITE(op,"(a)") ' => ERROR: 2theta min is larger than 2theta max.'
                  logi = .false.
                  RETURN
                END IF
                IF((th2_max - th2_min) < d_theta) THEN
                  WRITE(op,"(a)") ' => ERROR: 2theta increment is larger than 2theta max - 2theta min.'
                  logi = .false.
                  RETURN
                END IF

              !--- Calculation of n_high (maximum index of the array for pattern calculation)
              thmin=th2_min
              thmax=th2_max
              step_2th=d_theta
              n_high = nint((th2_max-th2_min)/d_theta+1.0_dp)
              write(*,"(a,3f7.2 )") " => 2Theta max, 2Theta min, stepsize:",  th2_min, th2_max, d_theta
              th2_min = th2_min * deg2rad
              th2_max = th2_max * deg2rad
              d_theta = half * deg2rad * d_theta
              ok_range=.true.

          Case ("SADP")
              k=index(txt," ")
              txt=adjustl(txt(k+1:))
              funct_num = 4
            	read (unit = txt, fmt = *, iostat=ier) i_plane, l_upper, loglin, brightness

            	if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading SADP parameters"
                  logi=.false.
                  return
              end if
            	if(i_plane < 1 .or. i_plane > 4) then
                WRITE(op,"(a)") ' => ERROR: Illegal reciprocal plane choice'
                logi = .false.
                return
              end if
              if(loglin < 0 .or. loglin > 1) then
                WRITE(op,"(a)") ' => ERROR: Illegal intensity scaling type.'
                logi = .false.
                return
              end if
              if(brightness <= ZERO) then
                WRITE(op,"(a)") ' => ERROR: Illegal value for brightness. Must be positive'
            	  logi = .false.
                return
              end if
              ok_sadp = .true.


          Case ("LMA")
            opt=4
            txt=adjustl(txt(k+1:))
            Cond%constr=.false.
            Cond%reached=.false.
            Cond%corrmax=50.0
            Cond%nfev=0
            Cond%njev=0
            Cond%icyc=0
            Cond%npvar=0
            Cond%iw=0
            Cond%tol=1.0e-5
            Cond%percent=0.0

          Case ("CORRMAX")
            txt=adjustl(txt(k+1:))
            read(unit=txt,fmt=*, iostat=ier) Cond%corrmax
            if(ier /= 0 ) then
              Err_crys=.true.
              Err_crys_mess="ERROR reading corrmax"
              logi=.false.
              return
            end if
            ok_corrmax=.true.

          Case ("MAXFUN")
            txt=adjustl(txt(k+1:))
            read(unit=txt,fmt=*, iostat=ier) Cond%icyc
            if(ier /= 0 ) then
              Err_crys=.true.
              Err_crys_mess="ERROR reading maxfun"
              logi=.false.
              return
            end if
            ok_maxfun=.true.

          Case ("TOL")
             txt=adjustl(txt(k+1:))
             read(unit=txt,fmt=*, iostat=ier) Cond%tol
             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading tolerance"
               logi=.false.
               return
             end if
             ok_tol=.true.

          Case ("NPRINT")
             txt=adjustl(txt(k+1:))
             read(unit=txt,fmt=*, iostat=ier) Cond%nprint
             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading nprint"
               logi=.false.
               return
             end if
             ok_nprint=.true.

            if (Cond%constr .and. ok_corrmax .and. ok_tol .and. ok_maxfun) then
              ok_lma=.true.
            elseif(ok_corrmax .and. ok_tol .and. ok_maxfun .and. ok_nprint .and. .not. Cond%constr ) then
              ok_lma=.true.
            end if


          Case ('LOCAL_OPTIMIZER')
            opt = 3
            opti%iquad = 1
            txt=adjustl(txt(k+1:))
            read(unit = txt, fmt = *, iostat=ier) opti%method
              if(ier /= 0 ) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading optimization method"
                logi=.false.
                return
              end if
              ok_opt=.true.

          Case ("MXFUN")
            txt=adjustl(txt(k+1:))
            read(unit = txt, fmt = *, iostat=ier) opti%mxfun
            if(ier /= 0 ) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading max no of function evaluations. Recommended 20 * no of parameters to refine"
                logi=.false.
                return
            end if
            ok_mxfun=.true.

          Case ("EPS")
            txt=adjustl(txt(k+1:))
            read(unit = txt, fmt = *, iostat=ier) opti%eps
            if(ier /= 0 ) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading stopping criterion"
                logi=.false.
                return
            end if
            ok_eps=.true.

          Case ("IOUT")
            txt=adjustl(txt(k+1:))
            read(unit = txt, fmt = *, iostat=ier) opti%iout
            if(ier /= 0 ) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading creation of output file instruction"
                logi=.false.
                return
            end if
            ok_iout=.true.

          Case ("ACC")
            txt=adjustl(txt(k+1:))
            read(unit = txt, fmt = *, iostat=ier) opti%acc
            if(ier /= 0 ) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading stopping accuracy"
                logi=.false.
                return
            end if
            ok_acc=.true.

          Case Default
            cycle
          End Select
      end do


     if (ok_sim .and. (ok_range .or. ok_sadp)) then
       return
     else if (ok_opt .and. ok_acc .and. ok_mxfun .and. ok_eps .and. ok_iout) then
       return
     else if (ok_lma) then
       return
     else if (ok_streak) then
       return
     else
        Err_crys=.true.
        Err_crys_mess="ERROR:  calculation parameters missing!"
        logi=.false.
     end if


    End Subroutine Read_CALCULATION

    Subroutine Read_EXPERIMENTAL(logi)
    logical, intent(out) :: logi


      integer :: i,i1,i2,k,j, ier, m
      character(len=132) :: txt
      character(len=25)  :: keyv
      logical :: ok_file, ok_fformat, ok_bgrnum, ok_bgrinter, ok_bgrcheb, ok_bgrpatt, ok_expstreak

      logi=.true.
      i1=sect_indx(8)
      i2=numberl

      i=i1
      m=0

      ok_file=.false.; ok_fformat=.false.; ok_bgrnum=.false.
      ok_bgrinter=.false.; ok_bgrcheb=.false.; ok_bgrpatt=.false.
      ok_expstreak=.false.
      crys%bgrpatt=.false.; crys%bgrinter=.false.; crys%bgrcheb=.false.
      do

        i=i+1; if(i > i2) exit
        txt=u_case(adjustl(tfile(i)))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        keyv=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(keyv)


        Case("FILE")   !Here read tfile that conserves the true file name

            read(unit=tfile(i)(k+1:),fmt=*, iostat=ier)   dfile, crys%patscal, ref_glb(1)
            val_glb(1)=crys%patscal
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading file instruction"
                  logi=.false.
                  return
              end if

              ok_file=.true.

          Case ('REPLACE_FILES')  !Tune the opening output files to replace previous files
             replace_files = .true.

          Case("EXCLUDED_REGIONS")
              read(unit=txt,fmt=*, iostat=ier) nexcrg
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading number of excluded region"
                  logi=.false.
                  return
              end if
              do j=1,nexcrg
                 i=i+1; if(i > i2) exit
                 txt=adjustl(tfile(i))
                 read(unit=txt,fmt=*, iostat=ier)   alow(j),ahigh(j)
                 if(ier /= 0 ) then
                     Err_crys=.true.
                     write(unit=Err_crys_mess,fmt="(a,i3)"),"ERROR reading the excluded region number ",j
                     logi=.false.
                     return
                 end if
              end do

          Case("FFORMAT")
              read(unit=txt,fmt=*,iostat=ier) fmode
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading file format instruction"
                  logi=.false.
                  return
              end if
              ok_fformat=.true.

          Case("EXPSTREAK")
              read(unit=txt,fmt=*,iostat=ier) num_streak, crys%patscal, ref_glb(1)
              val_glb(1)=crys%patscal
              allocate (file_streak(num_streak))
              allocate (h_streak(num_streak))
              allocate (k_streak(num_streak))
              allocate (l1_streak(num_streak))
              allocate (l0_streak(num_streak))
              allocate (dl_streak(num_streak))
              allocate (streak_flags(num_streak + 1))
              !crys%patscal = 1.0
              j=1   !count streak files
              do while (j <= num_streak)
                  i=i+1
                  !txt=u_case(adjustl(tfile(i)))
                  txt=adjustl(tfile(i)) !Conserve the original name of file_streaks (important for Linux and MacOS)
                  if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) then
                      !do nothink
                  else
                      read (unit = txt,fmt=*,iostat=ier) file_streak(j), h_streak(j), k_streak(j)
                      j=j+1
                  end if
              end do    ! end of do while(j <= num_streak)
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading experimental streak file instruction"
                  logi=.false.
                  return
              end if
              streakOrPowder = .true.
              ok_expstreak = .true.
              ok_file = .true.
              adapt_quad = 1

          Case("BGRINTER")
            crys%bgrinter=.true.
            mode = "INTERPOLATION"
            read(unit=tfile(i)(k+1:),fmt=*, iostat=ier) background_file
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading background linear interpolation instruction"
                  logi=.false.
                  return
              end if
              ok_bgrinter=.true.

          Case("BGRCHEB")
              crys%bgrcheb=.true.
              read(unit=txt,fmt=*,iostat=ier)  crys%cheb_nump
              do
                i=i+1
                txt=adjustl(tfile(i))
                if(txt(1:1) == "!") cycle
                exit
              end do
              read (unit=txt,fmt=*,iostat=ier) crys%chebp(1:crys%cheb_nump)
              val_glb(18:17+crys%cheb_nump)=crys%chebp(1:crys%cheb_nump)
              do
                i=i+1
                txt=adjustl(tfile(i))
                if(txt(1:1) == "!") cycle
                exit
              end do
              read(unit=txt,fmt=*, iostat=ier)ref_glb(18:17+crys%cheb_nump)
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading background Chebychev instruction"
                  logi=.false.
                  return
              end if
              ok_bgrcheb=.true.

          Case("BGRNUM")
              read(unit=txt,fmt=*,iostat=ier)  crys%num_bgrpatt
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading number of background patterns "
                  logi=.false.
                  return
              end if
              ok_bgrnum=.true.

          Case("BGRPATT")
            crys%bgrpatt=.true.
            m=m+1
            write(*,"(a,i3)") " => Reading background pattern #", m
            !First try to read the name of the hkl file
            read(unit=tfile(i)(k+1:),fmt=*, iostat=ier)crys%bfilepat(m), crys%bscalpat(m), &
                                                       ref_glb(17+crys%cheb_nump+m),crys%bfilehkl(m)
            if(ier == 0) then
                write(*,"(4a)") " => The background pattern and hkl files have been successfully read: ", &
                                trim(crys%bfilepat(m))," ", trim(crys%bfilehkl(m))
            else  !if not re-read only the background file name the scale factor and the code
                write(*,"(a,i3)") " => No hkl file was found for background pattern #", m
            read(unit=tfile(i)(k+1:),fmt=*, iostat=ier)crys%bfilepat(m), crys%bscalpat(m), &
                                                       ref_glb(17+crys%cheb_nump+m)
                if(ier == 0) then
                    write(*,"(2a)") " => The background pattern file has been successfully read: ", &
                                trim(crys%bfilepat(m))
                end if
            end if
            if(ier /= 0 ) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading background pattern instruction"
                logi=.false.
                return
            end if
            val_glb(17+crys%cheb_nump+m)= crys%bscalpat(m)
              ok_bgrpatt=.true.

          Case Default

             cycle

        End Select
        crys%num_bgrpatt=m
      end do

      LGBLT=17+crys%cheb_nump+crys%num_bgrpatt

        if(ok_file .and. ok_fformat .and. (ok_expstreak .or. ok_bgrinter .or. ok_bgrcheb .or. (ok_bgrnum .and. ok_bgrpatt) )) then
        return
      else
        Err_crys=.true.
        Err_crys_mess="ERROR: file details, file format, background file name or background calculation missing!"
        logi=.false.
      end if

    End Subroutine Read_EXPERIMENTAL

    Subroutine Set_Crys()

        crys%trm = .false.
        crys%finite_width = .false.

        crys%recrsv = .false.
        crys%xplcit = .false.
        crys%inf_thick = .false.

        crys%yy = max_a

        if (allocated (crys%a_name)) deallocate(crys%a_name)
        allocate(crys%a_name(max_a,max_l))
        crys%a_name=" "

        if (allocated (crys%a_num)) deallocate(crys%a_num)
        allocate(crys%a_num(max_a,max_l))
        crys%a_num=0

        if (allocated (crys%l_seq)) deallocate(crys%l_seq)
        allocate(crys%l_seq(xp_max))
        crys%l_seq = 0

        if (allocated (crys%l_actual)) deallocate(crys%l_actual)
        allocate(crys%l_actual(max_l))
        crys%l_actual=0

        if (allocated (crys%a_pos)) deallocate(crys%a_pos)
        allocate(crys%a_pos(3,max_a,max_l))
        crys%a_pos=0


        if (allocated (crys%a_B)) deallocate(crys%a_B)
        allocate(crys%a_B(max_a,max_l))
        crys%a_B=0

        if (allocated (crys%centro)) deallocate(crys%centro)
        allocate(crys%centro(max_l))
        crys%centro=0

        if (allocated (crys%high_atom)) deallocate(crys%high_atom)
        allocate(crys%high_atom(max_l))
        crys%high_atom=0

        if (allocated (crys%low_atom)) deallocate(crys%low_atom)
        allocate(crys%low_atom(max_l))
        crys%low_atom=0

        if (allocated (crys%l_n_atoms)) deallocate(crys%l_n_atoms)
        allocate(crys%l_n_atoms(max_l))
        crys%l_n_atoms=0

        if (allocated (crys%a_occup)) deallocate(crys%a_occup)
        allocate(crys%a_occup(max_a,max_l))
        crys%a_occup=0

        if (allocated (crys%l_alpha)) deallocate(crys%l_alpha)
        allocate(crys%l_alpha(max_a,max_l))
        crys%l_alpha=0


        if (allocated (crys%l_r)) deallocate(crys%l_r)
        allocate(crys%l_r(3,max_a,max_l))
        crys%l_r=0


        if (allocated (crys%r_b11)) deallocate(crys%r_b11)
        allocate(crys%r_b11(max_a,max_l))
        crys%r_b11=0

        if (allocated (crys%r_b22)) deallocate(crys%r_b22)
        allocate(crys%r_b22(max_a,max_l))
        crys%r_b22=0

        if (allocated (crys%r_b33)) deallocate(crys%r_b33)
        allocate(crys%r_b33(max_a,max_l))
        crys%r_b33=0

        if (allocated (crys%r_b12)) deallocate(crys%r_b12)
        allocate(crys%r_b12(max_a,max_l))
        crys%r_b12=0

        if (allocated (crys%r_b23)) deallocate(crys%r_b23)
        allocate(crys%r_b23(max_a,max_l))
        crys%r_b23=0

        if (allocated (crys%r_b13)) deallocate(crys%r_b13)
        allocate(crys%r_b13(max_a,max_l))
        crys%r_b13=0

        if (allocated (crys%chebp)) deallocate(crys%chebp)
        allocate(crys%chebp(max_n_cheb))
        crys%chebp=0


        if (allocated (crys%bfilepat)) deallocate(crys%bfilepat)
        allocate(crys%bfilepat(max_bgr_num))
        crys%bfilepat=" "

        if (allocated (crys%bfilehkl)) deallocate(crys%bfilehkl)
        allocate(crys%bfilehkl(max_bgr_num))
        crys%bfilehkl=" "

        if (allocated (crys%bscalpat)) deallocate(crys%bscalpat)
        allocate(crys%bscalpat(max_bgr_num))
        crys%bscalpat=0


        if(allocated (crys%list)) deallocate(crys%list)
        allocate(crys%list(max_npar))
        crys%list=0

        if(allocated (crys%cod)) deallocate(crys%cod)
        allocate(crys%cod(max_npar))
        crys%cod=0

        if(allocated (crys%p)) deallocate(crys%p)
        allocate(crys%p(max_npar))
        crys%p=0

        if(allocated (crys%mult)) deallocate(crys%mult)
        allocate(crys%mult(max_npar))
        crys%mult=0

        if(allocated (crys% Pv_refi)) deallocate(crys% Pv_refi)
        allocate(crys% Pv_refi(max_npar))
        crys% Pv_refi=0

        return
    End Subroutine Set_Crys
!______________________________________________________________________________________________________

     Function length(string) result(leng)    !from diffax

        character(len=*), intent(in) :: string
        integer :: leng
        integer :: i
        i=index(string," ")
        if( i==0) then
          leng=len_trim(string)
        else
          leng=i-1
        end if
      End Function length
!______________________________________________________________________________________________________

      Subroutine read_structure_file (namef,  logi)
        character(len=*), intent(in ):: namef
        logical,          intent(out):: logi
        !Local variables
        logical :: esta
        integer :: a,b,l,j, i, Lcode_max
        !integer :: n, ier

        logi = .true.
        call Set_Crys()
        b = 0
        a = 0
        l = 0
        j = 0

        lglb=0;  ref_glb=0.0
        latom=0; ref_atom=0.0
        ltrans=0; ref_trans=0.0



        inquire(file=namef,exist=esta)
        if( .not. esta) then
          Err_crys=.true.
          Err_crys_mess="ERROR :  The file "//trim(namef)//" doesn't exist"
          logi=.false.
          return
        !else
        !  open(unit=i_data,file=trim(namef),status="old",action="read",position="rewind",iostat=ier)
        !  if(ier /= 0) then
        !    Err_crys=.true.
        !    Err_crys_mess="ERROR opening the file "//trim(namef)
        !    logi=.false.
        !    return
        !  end if
        end if

        call Set_TFile(namef,logi)

         if(.not. logi) return

         call Read_TITLE(logi)
         if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
         end if

         call Read_INSTRUMENTAL(logi)
         if(.not. logi) then
           write(*,"(a)")  " => "//Err_crys_mess
           return
         end if

         call Read_STRUCTURAL(logi)
         if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
         end if

         call Read_LAYER(logi)
         if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
         end if

         call Read_STACKING (logi)
         if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
         end if

        call Read_TRANSITIONS (logi)
         if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
         end if

         call Read_CALCULATION (logi)

         if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
         end if


        if (opt /= 0) then !not necessary for simulation

          call Read_EXPERIMENTAL(logi)

          if(.not. logi) then
            write(*,"(a)")  " => "//Err_crys_mess
            return
          end if
        end if

        call Treat_codes(Lcode_max)
        call Modify_codes(Lcode_max)
        call Update_all(Lcode_max)
        crys%npar=np

        if (opt == 4) then         !construction of some optimization variables  (LMA)
          if (np==0) then                                                        !
            write(*,"(a)") " => WARNING: No parameters are being refined!"  !NIK
          end if                                                                 !
          opti%npar = Lcode_max
          Cond%npvar=opti%npar
          vs%pv(1:Lcode_max)= crys%Pv_refi(1:Lcode_max)
          do i=1,crys%npar
           vs%code(crys%p(i)) = 1
          end do
          vs%np= opti%npar
        end if

        return
      End subroutine  read_structure_file

      Subroutine Treat_codes(Lcode_max)
       integer, intent (out)     :: Lcode_max
       integer                   :: k, i, j, iyy

        Lcode_max=0
        do k=1, LGBLT
           Lglb(k) =  int(abs(ref_glb(k)/10.0))  !ordinal
           iyy= Lglb(k)
           if(iyy > Lcode_max) Lcode_max=iyy
           ref_glb(k) = ((abs(ref_glb(k))-(10.0*REAL(Lglb(k))))*SIGN(1.0,(ref_glb(k))))  !multiplicador
        end do

        do k=1, crys%n_actual
          do i=1, crys%l_n_atoms(k)
            do j=1,5
              Latom(j,i,k) = int(abs(ref_atom(j,i,k)/10.0))  !ordinal
              iyy= Latom(j,i,k)
              if(iyy > Lcode_max) Lcode_max=iyy
              ref_atom(j,i,k) = ((abs(ref_atom(j,i,k))-(10.0*REAL(Latom(j,i,k))))*SIGN(1.0,(ref_atom(j,i,k))))  !multiplicador
            end do
          end do
        end do

        do k=1, crys%n_typ
          do i=1, crys%n_typ
            do j=1, 10
              Ltrans(j, i, k) = int(abs(ref_trans(j,i,k)/10.0))  !ordinal
              iyy= Ltrans(j, i, k)
              if(iyy > Lcode_max) Lcode_max=iyy
              ref_trans(j, i, k) = ((abs(ref_trans(j,i,k))-(10.0*REAL(Ltrans(j, i, k))))*SIGN(1.0,(ref_trans(j,i,k))))  !multiplicador
            end do
          end do
        end do

      end Subroutine Treat_codes
!
      Subroutine Restore_Codes()
        !---- Local Variables ----!
        integer :: i, j, k

        do i=1,LGBLT
          ref_glb(i)=sign(1.0_cp,ref_glb(i))*(real(10*lglb(i))+abs(ref_glb(i)))
        end do
        do k=1, crys%n_actual
          do i=1, crys%l_n_atoms(k)
            do j=1,5
              ref_atom(j,i,k) =sign(1.0_cp,ref_atom(j,i,k))*(real(10*latom(j,i,k))+abs(ref_atom(j,i,k)))
            end do
          end do
        end do

        do k=1, crys%n_typ
          do i=1, crys%n_typ
            do j=1, 10
              ref_trans(j,i,k) =sign(1.0_cp,ref_trans(j,i,k))*(real(10*ltrans(j,i,k))+abs(ref_trans(j,i,k)))
            end do
          end do
        end do

        return
      End Subroutine Restore_Codes

     Subroutine Modify_codes(Lcode_max)
       integer, intent (in out) :: Lcode_max
       integer                  :: k, l, i, j, n_given,  &
                                   n_att, ndisp, nn, ni, ndispm, &
                                   maxs !,n_pat, kk, iof, nd
       integer, dimension(Lcode_max) :: dispo
       real(kind=cp), parameter      :: e=0.001
       !
       !  Check correlated parameters and already used codes
       !
       ndisp=0
       n_given=0
       dispo(:)=0
       ! First Pass
        Do L=1, Lcode_max
          ni=0
          do k=1, LGBLT
            IF (L == Lglb(k)) ni=ni+1
          end do

          do k=1, crys%n_actual
            do i=1, crys%l_n_atoms(k)
              do j=1,5
                IF (L == Latom(j,i,k)) ni=ni+1
              end do
            end do
          end do

          do k=1, crys%n_typ
            do i=1, crys%n_typ
              do j=1, 10
                IF (L == Ltrans(j, i, k)) ni=ni+1
              end do
            end do
          end do

          if(ni==0) then
            ndisp=ndisp+1  !number of disponible codes
            dispo(ndisp)=L  !disponible code number
          else
            n_given=n_given+1  !number of given codes
          end if
!
       End Do !=1,Lcode_max
!
       !if(ndisp == 0) return  !all codes have been attributed

        !
        ! Attributing numbers to parameters (not already attributed) with multipliers
        ! different from zero. First the attribution is taken from the vector disp() and
        ! continued, after finishing the disponible codes, from Lcode_max+1, ...
        ! If after attributing the codes ndisp /=0, then a displacement of all parameters
        ! is done and the maximum number of parameters to be refined is diminished by
        ! ndisp
        !
        ndispm=ndisp !number of disponible codes before attributing code numbers
        ni=0
        nn=Lcode_max
        n_att=0
       ! Test Global Parameters
         do j =1,LGBLT
           if (abs(ref_glb(j)) > e .and. lglb(j) == 0 ) then
             if(abs(ref_glb(j)) > 1.001) then
                ref_glb(j)=sign(1.0_cp,ref_glb(j))*(abs(ref_glb(j))-1.0)
             end if
             if(ndisp==0) then
               nn=nn+1
               lglb(j)=nn
               n_att=n_att+1
             else
               ni=ni+1
               lglb(j)=dispo(ni)
               n_att=n_att+1
               ndisp=ndisp-1
             end if
           end if
         end do

         ! Test Atom Parameters

         do k=1, crys%n_actual
           do i=1, crys%l_n_atoms(k)
             do j=1,5
                if(abs(ref_atom(j,i,k)) > e .and. Latom(j,i,k) == 0 ) then
                  if(abs(ref_atom(j,i,k)) > 1.001) then
                     ref_atom(j,i,k)=sign(1.0_cp,ref_atom(j,i,k))*(abs(ref_atom(j,i,k))-1.0)
                  end if
                  if(ndisp==0) then
                    nn=nn+1
                    Latom(j,i,k)=nn
                    n_att=n_att+1
                  else
                    ni=ni+1
                    Latom(j,i,k)=dispo(ni)
                    n_att=n_att+1
                    ndisp=ndisp-1
                  end if
                end if
             end do
           end do
         end do

         ! Test Transition Parameters

         do k=1, crys%n_typ
           do i=1, crys%n_typ
             do j=1, 10
                if(abs(ref_trans(j,i,k)) > e .and. Ltrans(j,i,k) == 0 ) then
                  if(abs(ref_trans(j,i,k)) > 1.001) then
                     ref_trans(j,i,k)=sign(1.0_cp,ref_trans(j,i,k))*(abs(ref_trans(j,i,k))-1.0)
                  end if
                  if(ndisp==0) then
                    nn=nn+1
                    Ltrans(j,i,k)=nn
                    n_att=n_att+1
                  else
                    ni=ni+1
                    Ltrans(j,i,k)=dispo(ni)
                    n_att=n_att+1
                    ndisp=ndisp-1
                  end if
                end if
             end do
           end do
         end do

         if(ndisp == 0) then
           !maxs=nn   !this is an error!!!!
           if(nn > Lcode_max) Lcode_max=nn
           return  !all parameters have been attributed
         end if

         maxs=n_given+n_att    !Number of refined parameters (given + attributed)

         ! Third Pass ndisp /=0 => Displacement of codes needed to avoid holes in the matrix.

       n_att=ndispm-ndisp+1


       DO L=Lcode_max, maxs+1,-1
        nn=0
         do j =1,LGBLT
           if (L == lglb(j)) then
             lglb(j)=dispo(n_att)
             nn=nn+1
           end if
         end do

         do k=1, crys%n_actual
           do i=1, crys%l_n_atoms(k)
             do j=1,5
                if(L == Latom(j,i,k)) then
                  Latom(j,i,k)=dispo(n_att)
                  nn=nn+1
                end if
             end do
           end do
         end do

         do k=1, crys%n_typ
           do i=1, crys%n_typ
             do j=1, 10
                if(L == Ltrans(j,i,k)) then
                  Ltrans(j,i,k)=dispo(n_att)
                  nn=nn+1
                end if
             end do
           end do
         end do

         if(nn == 0) cycle
         n_att=n_att+1
         if(n_att > ndispm) exit
       END DO !i=Lcode_max,maxs+1,-1
       Lcode_max=maxs

       return
     End Subroutine Modify_codes

     Subroutine Update_all(Lcode_Max)
     integer, intent (in out) :: Lcode_max
     integer                  :: L,i,j,k


     np=0
     DO L=1,Lcode_max

       do j =1,LGBLT
         if (L == lglb(j)) then
           np=np+1
           crys%list(np)=val_glb(j)
           crys%p(np)=lglb(j)
           crys%Pv_refi(crys%p(np))=val_glb(j)
           crys%mult(np)=ref_glb(j)
           if(j <= 17) then
              namepar(np) = names_glb(j)
           else if (j> 17 .and. j <= 17+crys%cheb_nump ) then
              write(namepar(np),"(a,i2.2)") "ChebCoeff_",j-17
           else if (j> 17+crys%cheb_nump .and. j <= LGBLT ) then
              write(namepar(np),"(a,i2.2)") "Bkg_Scale_",j-(17+crys%cheb_nump)
           end if
           vs%nampar(crys%p(np))=namepar(np)
         end if
       end do

       do k=1, crys%n_actual
         do i=1, crys%l_n_atoms(k)
           do j=1,5
              if(L == Latom(j,i,k)) then
                np=np+1
                crys%list(np)=val_atom(j,i,k)
                crys%p(np)=Latom(j,i,k)
                crys%Pv_refi(crys%p(np))=val_atom(j,i,k)
                crys%mult(np)=ref_atom(j,i,k)
                Select Case(j)
                  Case(1) ; write(namepar(np),"(a,2i2.2)") 'pos_x',i,k
                  Case(2) ; write(namepar(np),"(a,2i2.2)") 'pos_y',i,k
                  Case(3) ; write(namepar(np),"(a,2i2.2)") 'pos_z',i,k
                  Case(4) ; write(namepar(np),"(a,2i2.2)") 'Biso',i,k
                  case(5) ; write(namepar(np),"(a,2i2.2)") 'Occ',i,k
                End Select
                vs%nampar(crys%p(np))=namepar(np)
              end if
           end do
         end do
       end do

       do k=1, crys%n_typ
         do i=1, crys%n_typ
           do j=1, 10
              if(L == Ltrans(j,i,k)) then
                np=np+1
                crys%list(np)=val_trans(j,i,k)
                crys%p(np)=Ltrans(j,i,k)
                crys%Pv_refi(crys%p(np))=val_trans(j,i,k)
                crys%mult(np)=ref_trans(j,i,k)
                Select Case (j)
                   Case(1) ; write(namepar(np),"(a,2i2.2)")"alpha",k,i
                   Case(2) ; write(namepar(np),"(a,2i2.2)")"tx"   ,k,i
                   Case(3) ; write(namepar(np),"(a,2i2.2)")"ty"   ,k,i
                   Case(4) ; write(namepar(np),"(a,2i2.2)")"tz"   ,k,i
                   Case(5) ; write(namepar(np),"(a,i2.2)") "FW_11",k,i
                   Case(6) ; write(namepar(np),"(a,i2.2)") "FW_22",k,i
                   Case(7) ; write(namepar(np),"(a,i2.2)") "FW_33",k,i
                   Case(8) ; write(namepar(np),"(a,i2.2)") "FW_12",k,i
                   Case(9) ; write(namepar(np),"(a,i2.2)") "FW_31",k,i
                   Case(10); write(namepar(np),"(a,i2.2)") "FW_23",k,i
                End Select
                vs%nampar(crys%p(np))=namepar(np)
              end if
           end do
         end do
       end do

     END DO !i=Lcode_max,maxs+1,-1


     End Subroutine Update_all

        ! TemplateFromCIF
        ! Aim: write simple .flst file form data stored in fileName.cif
        ! Usage: to faults program one should provide a fileName.cif with additional parameter, which corresponds to
        !        number of layers.
        ! Remember that: 1) structure described in .cif file should have "P1" spacegroup.
        !                2) condition to check whether the atom relates to layers is:
        !                   (myListAtom%atom(i)%x(3) < thickness*k) .AND. (myListAtom%atom(i)%x(3) >= thickness*(k-1))
        !                   which means that plane with higher z-coordinte is not included into that layer,
        !                   but included into the next one

        subroutine TemplateFromCIF (fileName, input_n_layers, thickness)

            !---Input/Output---!
            character(len=*), intent(in)    :: fileName
            integer,          intent(in)    :: input_n_layers
            real,             intent(in)    :: thickness

            !---Local variables. Primitives---!
            integer                                         :: i,j,k, numLines, ini, numOfAtom
            !integer                                         :: i_atom
            integer, dimension(:),allocatable               :: n_zero_atom
            real,dimension(6)                               :: Cellda
            !real                                            :: delta
            character(len=100),dimension(:),allocatable     :: arrayOfStrings
            character(len=100)                              :: SG !name of SpaceGroup
            !logical                                         :: match = .false.

            !---Local variables. Types---!
            Type(Atom_List_Type)                             ::  myListAtom
            Type(Crystal_Cell_Type)                          ::  myCell
            Type(Space_Group_Type)                           ::  mySG
            Type(Atom_Equiv_List_Type)                       ::  myELA
            !Type(Atom_Equiv_Type), dimension(input_n_layers) ::  atom

            !---Subroutine---!

            !read cif file and initialize structure variables
            call Set_Crys()
            call Number_Lines(fileName,numLines)
            allocate(arrayOfStrings(numLines))
            call Reading_Lines(fileName,numLines,arrayOfStrings)
            ini=1
            call Read_CIF_Cell(arrayOfStrings,ini,numLines,Cellda)
            call Set_Crystal_Cell(Cellda(1:3),Cellda(4:6),myCell)
            ini=1
            call Read_CIF_Atom(arrayOfStrings,ini,numLines,numOfAtom,myListAtom)
            ini=1
            call Read_CIF_HM(arrayOfStrings,ini,numLines,SG)                            !find and obtain spacegroup name from the array of strings
            call Set_SpaceGroup(SG,mySG)
            call Set_Atom_Equiv_List(mySG,myCell,myListAtom,myELA)
            write (unit = crys%ttl, fmt=*) "Template from CIF file: ", fileName

            !default parameters
            crys%rad_type = 0
            crys%lambda = 1.540560
            crys%lambda2 = 1.544390
            crys%ratio = 0.5
            crys%cell_a = Cellda(1)
            crys%cell_b = Cellda(2)
            crys%cell_c = Cellda(3)
            crys%cell_gamma = Cellda(6)*deg2rad     !needed for correct write into .flts file
            crys%broad = PS_VGT
            crys%p_u = 0.004133
            crys%p_v = -0.007618
            crys%p_w = 0.006255
            crys%p_x = 0.018961
            crys%p_dg = 50000
            crys%p_dl = 50000
            crys%patscal = 1.0
            crys%trm = .true.
            crys%bckg_level = 100.0
            crys%sym = "-1"
            crys%SymGrpNo = 1
            crys%recrsv = .true.
            crys%inf_thick = .true.

            crys%n_typ = input_n_layers
            if (allocated (crys%fundamental)) deallocate(crys%fundamental)
            allocate(crys%fundamental(crys%n_typ))
            if (allocated (crys%original)) deallocate(crys%original)
            allocate(crys%original(crys%n_typ))
            do i = 1, input_n_layers
                crys%fundamental(i)=.true.
                crys%original(i) = i
            end do
            if (allocated (fundamental)) deallocate(fundamental)
            allocate(fundamental(crys%n_typ))
            if (allocated (original)) deallocate(original)
            allocate(original(crys%n_typ))
            do i = 1, input_n_layers
                fundamental(i)=.true.
                original(i) = i
            end do
            if (allocated (n_zero_atom)) deallocate(n_zero_atom)
            allocate(n_zero_atom(crys%n_typ))
            !initialize all numbers of atoms equal to zero
            do k = 1, MAX_L
                crys%l_n_atoms = 0
            end do
            !now iterate over all atoms in atom list and find that which references to some layer
            do i = 1, myListAtom%natoms
                do k = 1, input_n_layers
                    if ( (myListAtom%atom(i)%x(3) < thickness*k) .AND. (myListAtom%atom(i)%x(3) >= thickness*(k-1)) ) then
                        crys%l_n_atoms(k) = crys%l_n_atoms(k) + 1
                        crys%a_name  (  crys%l_n_atoms(k),k) = myListAtom%atom(i)%ChemSymb
                        crys%a_num   (  crys%l_n_atoms(k),k) = crys%l_n_atoms(k)
                        crys%a_pos   (1,crys%l_n_atoms(k),k) = myListAtom%atom(i)%x(1)
                        crys%a_pos   (2,crys%l_n_atoms(k),k) = myListAtom%atom(i)%x(2)
                        crys%a_pos   (3,crys%l_n_atoms(k),k) = myListAtom%atom(i)%x(3)
                        crys%a_occup (  crys%l_n_atoms(k),k) = myListAtom%atom(i)%occ
                    end if
                end do
            end do
            !find zero atom coordinates, zero atoms is determined as atom with lowest z-coordinate within layer
            do k = 1, input_n_layers
                n_zero_atom(k) = 1
            end do
            do k = 1, input_n_layers
                do j = 1, crys%l_n_atoms(k)
                    if ( &!(crys%a_pos(1,j,k) <= crys%a_pos(1,n_zero_atom(k),k)) .AND. &
                         !(crys%a_pos(2,j,k) <= crys%a_pos(2,n_zero_atom(k),k)) .AND. &
                         (crys%a_pos(3,j,k) <= crys%a_pos(3,n_zero_atom(k),k)) ) then
                        n_zero_atom(k) = j
                    end if
                end do
            end do
            !as far as coordinates of atoms are not transformed to ones in zeroatoms reference system
            !we need to calculate transitions vector between the layers
            do i = 1, input_n_layers-1
                crys%l_r(1,i+1,i) = crys%a_pos(1,n_zero_atom(i+1),i+1) - crys%a_pos(1,n_zero_atom(i),i)
                crys%l_r(2,i+1,i) = crys%a_pos(2,n_zero_atom(i+1),i+1) - crys%a_pos(2,n_zero_atom(i),i)
                crys%l_r(3,i+1,i) = crys%a_pos(3,n_zero_atom(i+1),i+1) - crys%a_pos(3,n_zero_atom(i),i)
                crys%l_alpha(i+1,i) = 1.0
            end do
            crys%l_r(1,1,input_n_layers) = crys%a_pos(1,n_zero_atom(1             ),1             ) -    &
                                           crys%a_pos(1,n_zero_atom(input_n_layers),input_n_layers)
            crys%l_r(2,1,input_n_layers) = crys%a_pos(2,n_zero_atom(1             ),1             ) -    &
                                           crys%a_pos(2,n_zero_atom(input_n_layers),input_n_layers)
            crys%l_r(3,1,input_n_layers) = crys%a_pos(3,n_zero_atom(1             ),1             )    &
                                         - crys%a_pos(3,n_zero_atom(input_n_layers),input_n_layers) + one
            crys%l_alpha(1,input_n_layers) = 1.0
            !this is needed to write a correct coordinates to .flts file
            do k = 1, input_n_layers
                do i = 1, crys%l_n_atoms(k)
                    crys%a_pos(1, i, k) = crys%a_pos(1, i, k)*pi2
                    crys%a_pos(2, i, k) = crys%a_pos(2, i, k)*pi2
                    crys%a_pos(3, i, k) = crys%a_pos(3, i, k)*pi2
                end do
            end do
            !now we can transform the coordinates of atoms to the reference of system zero atom
            do k = 1, input_n_layers
                do i = 1, crys%l_n_atoms(k)
                    crys%a_pos(1, i, k) = crys%a_pos(1, i, k) - crys%a_pos(1,n_zero_atom(k),k)
                    crys%a_pos(2, i, k) = crys%a_pos(2, i, k) - crys%a_pos(2,n_zero_atom(k),k)
                    crys%a_pos(3, i, k) = crys%a_pos(3, i, k) - crys%a_pos(3,n_zero_atom(k),k)
                end do
            end do

            !by default put a powder calculation command
            opt = 0
            funct_num = 3
            thmin = 10
            thmax = 120
            step_2th = 0.01
            !ask program to write transition parameters in form of a table
            do i=1,4
                table(i) = .true.
            end do
        end subroutine TemplateFromCIF
    End Module read_data
