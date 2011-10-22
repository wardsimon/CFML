    Module read_data
       use CFML_GlobalDeps,           only : sp , cp
       use CFML_String_Utilities,     only : number_lines , reading_lines , findfmt ,iErr_fmt, getword, err_string, &
                                             err_string_mess, getnum, Ucase
       use CFML_Optimization_General, only : Nelder_Mead_Simplex,  Opt_Conditions_Type
       use CFML_Simulated_Annealing
       use diffax_mod
       implicit none

       private
       !public subroutines
       public:: read_structure_file, new_getfil , length, prune, rdrnge, choice

!      Declaration of diffraction_pattern_type

       type,public :: crys_2d_type
       integer                                          :: rad_type                     !radiation type   >  rad_type
       real                                             :: lambda , lambda2 , ratio     !radiation wavelength    > lambda
       integer                                          :: broad                        !type of instrumental broadening>  blurring
       real                                             :: cell_a                       !cell parameter a   >cell_a
       real                                             :: cell_b                       !cell parameter b    >cell_b
       real                                             :: cell_c                       !cell parameter c   >cell_c
       real                                             :: cell_gamma                   !gamma   >cell_gamma
       real                                             :: layer_a = 0.0, layer_b = 0.0 !layer characteristic widths (optional)
       real                                             :: ref_layer_a, ref_layer_b, rang_layer_a, rang_layer_b
       real                                             :: p_u = 0.0                    !u  > pv_u
       real                                             :: p_v = 0.0                    !v  > pv_v
       real                                             :: p_w = 0.0                    !w  > pv_w
       real                                             :: fwhm = 0.0                   !fwmh ( cannot be negative)  >fwhm
       real                                             :: p_gamma = 0.0                !pseudo-voigt gamma    > pv_gamma
       real                                             :: p_x                          ! x > pv_x
       real                                             :: p_dg, p_dl                   ! gaussian and lorentzian average volumetric size
       real                                             :: ref_p_u, ref_p_v, ref_p_w, ref_p_x, ref_p_dg, ref_p_dl
       real                                             :: rang_p_v, rang_p_u, rang_p_w, rang_p_x, rang_p_dg, rang_p_dl
       real                                             :: tolerance                    !d-> Maximum deviation that intensities can have
       real                                             :: rang_cell_a ,rang_cell_b,rang_cell_c,rang_cell_gamma
       real                                             :: t_ini=0.0, anneal, accept        ! In SAN: initial temperature, kirkpatrick factor of annealing, minimum number of accepted configurations
       character (len=132)                              :: sym                          !Laue symmetry  >pnt_grp
       character (len=132)                              :: calctype                     !type of calculation
       real                                             :: ref_cell_a,ref_cell_b, ref_cell_c, ref_cell_gamma  ! index of refinement: if zero, no refinement, if one, refinement to be done
       integer                                          :: n_typ                        !number of layer types    >n_layers
       real                                             :: l_cnt = 0                    !number of layers in stacking section
       real                                             :: ref_l_cnt
       real                                             :: rang_l_cnt
       integer                                          :: SymGrpNo
       integer                                          :: n_cycles=0
       integer                                          :: n_actual                     !number of unique layers
       integer                                          :: npar = 0                     !total number of parameters to be refined
       integer                                          :: num_temp=0                     ! maximum number of temperatures in SAN
       integer                                          :: init_config                  ! flag to randomly initialaze SAN configuration if = 1.
       logical                                          :: recrsv                       !layer sequencing recursive
       logical                                          :: xplcit                       !layer sequencing explicit
       logical                                          :: finite_width                 !layer width in plane ab
       logical                                          :: inf_thick                    !infinit stacking
       logical                                          :: trm                          !trim keyword: to supress peak at origin >trim_origin
       character(len=4), dimension(:,:),   allocatable  :: a_name                       !atom name (4characters wide)>atom_l o a_name?
       integer,          dimension(:),     allocatable  :: centro                       !layer simetry : none or centrosymmetric  >only_real
       integer,          dimension(:),     allocatable  :: cod
       integer,          dimension(:,:),   allocatable  :: a_num                        !atom number  > a_number
       real,             dimension(:,:,:), allocatable  :: ref_a_pos
       integer,          dimension(:)  ,   allocatable  :: l_seq                        !d-> array containing the explicitly defined sequence of layers.
       integer,          dimension(:)  ,   allocatable  :: l_actual                     !Contains the layer number that layer i is structurally identical to.
       integer,          dimension(:)  ,   allocatable  :: l_n_atoms                    !number of atoms in each layer
       integer,          dimension(:),     allocatable  :: p                            ! refinable parameter index
       real,             dimension(:),     allocatable  :: mult                         ! refinable parameter multiplicator
       real,             dimension(:,:),   allocatable  :: ref_l_alpha                  !index of refinement of l_alpha
       real,             dimension(:,:,:), allocatable  :: ref_l_r                      ! refinement index of transitions vector
       real,             dimension(:,:,:), allocatable  :: rang_a_pos
       real,             dimension(:,:,:), allocatable  :: a_pos                        ! xyz coordinates of the atoms   >a_pos
       real,             dimension(:,:),   allocatable  :: a_B                          !debye waller factors
       real,             dimension(:,:),   allocatable  :: ref_a_B                      !index of refinement of debye waller factors
       real,             dimension(:,:),   allocatable  :: rang_a_B                     !range of refinement of debye waller factors
       real,             dimension(:,:),   allocatable  :: a_occup                      !occupation factor (between 0 and 1)
       real,             dimension(:)  ,   allocatable  :: high_atom, low_atom          !lower and upper atom positions
       real,             dimension(:,:),   allocatable  :: l_alpha                      !l_alpha
       real,             dimension(:,:),   allocatable  :: rang_l_alpha                 !range of refinement of l_alpha
       real,             dimension(:,:,:), allocatable  :: l_r                          !transitions vector
       real,             dimension(:,:,:), allocatable  :: rang_l_r                     !range of refinement of transitions vector
       real,             dimension(:,:),   allocatable  :: r_b11  , r_B22 , r_B33 , r_B12 , r_B23, r_B31
       real,             dimension(:)  ,   allocatable  ::  list                       !vector containing all the parameters to optimize
       real,             dimension(:),     allocatable  :: vlim1                        !Low-limit value of parameters
       real,             dimension(:),     allocatable  :: vlim2                        !Low-limit value of parameters
       end  type crys_2d_type

       integer, parameter ,private                      :: i_data=30
       logical,            public, save                 :: Err_crys=.false.
       character(len=120), public, save                 :: Err_crys_mess=" "

   contains


!_____________________________________________________________________________________________________
    Subroutine new_getfil(stfile, l_crys, l_opti,l_san,vecsan, olg)


       character(len=31) ,             intent (in out)      :: stfile
       type (crys_2d_type),            intent (   out)      :: l_crys
       type (Opt_Conditions_Type),     intent (   out)      :: l_opti
       type (SimAnn_Conditions_type),  intent (   out)      :: l_san
       type (State_Vector_Type),       intent (   out)      :: vecsan
       logical                   ,     intent (   out)      :: olg

        write(unit=op,fmt="(a)",advance="no") ' => Enter the complete name of the structure input file: '
        read(unit= *,fmt="(a)") stfile

        WRITE(op,*) " => Looking for scattering factor data file '",  sfname(:),"'"
        OPEN(UNIT = sf, FILE = sfname)
        WRITE(op,*) " => Opening scattering factor data file '",  sfname(:),"'"

        call read_structure_file    (stfile, l_crys, l_opti, l_san,vecsan, olg)
        if (err_crys) print*, trim(err_crys_mess)

        return
    End subroutine
!______________________________________________________________________________________________________

      INTEGER*4 FUNCTION choice(flag, list, n)          !from diffax_read


      IMPLICIT NONE
      CHARACTER (LEN=*), INTENT(IN)        :: flag
      CHARACTER (LEN=80), INTENT(IN)       :: list(n)
      INTEGER*4, INTENT(IN)                :: n

      INTEGER*4 i, j1, j2

      i = 1
      j1 = length(flag)

      10 j2 = length(list(i))
! see if the string contained in list(i) is identical to that in flag
      IF(j1 == j2 .AND. INDEX(flag, list(i)(1:j2)) == 1) GO TO 20
      i = i + 1
      IF(i <= n) GO TO 10

      20 IF(i > n) i = -1
      choice = i

      RETURN
      END FUNCTION choice
!______________________________________________________________________________________________________

      LOGICAL FUNCTION rdrnge()      !from diffax_read

      rdrnge = .false.

      123 WRITE(op,"(a)") ' Enter angular range:'
      WRITE(op,"(a)") '   2theta min, 2theta max, 2theta increment : '
      READ(cntrl,*,ERR=123) th2_min, th2_max, d_theta
      IF(cfile) WRITE(op,"(1X, g11.4, 2(3X, g11.4))") th2_min, th2_max, d_theta

      IF(th2_min < zero) THEN
        WRITE(op,"(a)") ' ERROR: 2theta min is negative.'
        RETURN
      END IF

      IF(th2_min >= one_eighty) THEN
        WRITE(op,"(a)") ' ERROR: 2theta min exceeds 180 degrees.'
        RETURN
      END IF

      IF(th2_max <= zero) THEN
        WRITE(op,"(a)") ' ERROR: Negative (or zero) value for 2theta min.'
        RETURN
      END IF

      IF(th2_max > one_eighty) THEN
        WRITE(op,"(a)") ' ERROR: 2theta max exceeds 180 degrees.'
        RETURN
      END IF

! if th2_max = 180, reduce it slightly to keep us out of trouble
      IF(th2_max == one_eighty) th2_max = th2_max - eps4

      IF(d_theta <= zero) THEN
        WRITE(op,"(a)") ' ERROR: Negative (or zero) value for 2theta increment.'
        RETURN
      END IF

      IF(th2_min >= th2_max) THEN
        WRITE(op,"(a)") ' ERROR: 2theta min is larger than 2theta max.'
        RETURN
      END IF

      IF((th2_max - th2_min) < d_theta) THEN
        WRITE(op,"(a)") ' ERROR: 2theta increment is larger than 2theta max - 2theta min.'
        RETURN
      END IF

      th2_min = th2_min * deg2rad
      th2_max = th2_max * deg2rad
      d_theta = half * deg2rad * d_theta

      rdrnge = .true.
      RETURN

      END FUNCTION rdrnge

!______________________________________________________________________________________________________

      INTEGER*4 FUNCTION prune(line)          !from diffax_read


      CHARACTER (LEN=*), INTENT(IN OUT)        :: line


      INTEGER*4 lin_len, i

      prune = 0

      lin_len = LEN(line)
      i = lin_len
      10 IF(i > 0) THEN
        IF(line(i:i) == ' ') THEN
          IF(i > 1) THEN
            i = i - 1
            GO TO 10
          END IF
        END IF
      END IF

      IF(i > 0) prune = i

      RETURN
      END FUNCTION prune
!______________________________________________________________________________________________________

     Function length(string) result(leng)    !from diffax_read

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

        Subroutine read_structure_file (namef, crys, opti,opsan,sanvec, logi)

        character(len=*)    ,          intent(in    )     :: namef
        type (crys_2d_type ),          intent(   out)     :: crys
        type (Opt_Conditions_Type),    intent(   out)     :: opti
        type (State_Vector_Type),      intent(   out)     :: sanvec
        type (SimAnn_Conditions_type), intent(   out)     :: opsan
        logical                   ,    intent(   out)     :: logi



        logical                                        :: esta
        character (len=132)                            :: txt ,   broad , layer , random, semirandom, seq
        integer                                        :: numberl, i1,i2    !number of lines
        integer                                        :: ier,i, j,  a1, a2, k , l , z ,m ,r, pos1, pos2, iflag , &
                                                          y, q, l1, l2, j1, j2, yy
        integer,dimension(132)                         :: Inte    !  Vector of integer numbers
        integer                                        :: n_int  !number of integers per line
        real(kind=sp),      dimension(132)             :: rel !  Vector of integer numbers
        real,               dimension(132)             :: num_real
        character(len=132), dimension(:),allocatable   :: tfile
        integer, parameter                             :: eps =1.0E-4 , nrp=80
        integer, parameter                             :: read_out = 30
        character(len=132), dimension(15)              :: word  = " "! Out -> Vector of Words
        integer                                        :: n_word   ! Out -> Number of words
        real                                           :: sum , tmp , ab
        character (len = 4)                            :: trima
        character (len = 80 )                          :: list (1:12)
        real              , dimension(80)              :: label


        crys%trm = .false.
        crys%finite_width = .false.

        crys%recrsv = .false.
        crys%xplcit = .false.
        crys%inf_thick = .false.

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

        if (allocated (crys%ref_a_pos)) deallocate(crys%ref_a_pos)
        allocate(crys%ref_a_pos(3,max_a,max_l))
        crys%ref_a_pos=0

        if (allocated (crys%rang_a_pos)) deallocate(crys%rang_a_pos)
        allocate(crys%rang_a_pos(3,max_a,max_l))
        crys%rang_a_pos=0

        if (allocated (crys%a_B)) deallocate(crys%a_B)
        allocate(crys%a_B(max_a,max_l))
        crys%a_B=0

        if (allocated (crys%ref_a_B)) deallocate(crys%ref_a_B)
        allocate(crys%ref_a_B(max_a,max_l))
        crys%ref_a_B=0

        if (allocated (crys%rang_a_B)) deallocate(crys%rang_a_B)
        allocate(crys%rang_a_B(max_a,max_l))
        crys%rang_a_B=0

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

        if (allocated (crys%ref_l_alpha)) deallocate(crys%ref_l_alpha)
        allocate(crys%ref_l_alpha(max_a,max_l))
        crys%ref_l_alpha=0

        if (allocated (crys%rang_l_alpha)) deallocate(crys%rang_l_alpha)
        allocate(crys%rang_l_alpha(max_a,max_l))
        crys%rang_l_alpha=0

        if (allocated (crys%l_r)) deallocate(crys%l_r)
        allocate(crys%l_r(3,max_a,max_l))
        crys%l_r=0

        if (allocated (crys%ref_l_r)) deallocate(crys%ref_l_r)
        allocate(crys%ref_l_r(3,max_a,max_l))
        crys%ref_l_r=0

        if (allocated (crys%rang_l_r)) deallocate(crys%rang_l_r)
        allocate(crys%rang_l_r(3,max_a,max_l))
        crys%rang_l_r=0

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

        if (allocated (crys%r_b31)) deallocate(crys%r_b31)
        allocate(crys%r_b31(max_a,max_l))
        crys%r_b31=0

        if(allocated (crys%vlim1)) deallocate(crys%vlim1)
        allocate(crys%vlim1(nrp))
        crys%vlim1=0

        if(allocated (crys%vlim2)) deallocate(crys%vlim2)
        allocate(crys%vlim2(nrp))
        crys%vlim2=0

        if(allocated (crys%list)) deallocate(crys%list)
        allocate(crys%list(nrp))
        crys%list=0

        if(allocated (crys%cod)) deallocate(crys%cod)
        allocate(crys%cod(nrp))
        crys%cod=0

        if(allocated (crys%p)) deallocate(crys%p)
        allocate(crys%p(nrp))
        crys%p=0

        if(allocated (crys%mult)) deallocate(crys%mult)
        allocate(crys%mult(nrp))
        crys%mult=0

        yy = 0

        inquire(file=namef,exist=esta)
        if( .not. esta) then
          Err_crys=.true.
          Err_crys_mess="ERROR :  The file "//trim(namef)//" doesn't exist"
          logi = .false.
          return
        else
          open(unit=i_data,file=trim(namef),status="old",action="read",position="rewind",iostat=ier)
          if(ier /= 0) then
            Err_crys=.true.
            Err_crys_mess="ERROR opening the file "//trim(namef)
            logi = .false.
            return
          end if
        end if

        call number_lines(namef, numberl)                            ! count number of lines in file
        if (numberl==0) then
          err_crys=.true.
          err_crys_mess="ERROR :  The file "//trim(namef)//" contains nothing"
          logi = .false.
          return
        else
         if (allocated (tfile)) deallocate (tfile)
         allocate (tfile(numberl))
         tfile=" "
        call reading_lines(namef, numberl, tfile)                    ! we 'charge' the file in tfile  so we can close the unit
        end if
        close(unit=i_data)

        do i = 1, numberl                 ! To read in case insensitive mode
          call Ucase(tfile(i))
        end do

!    begin to read
        logi = .true.

        Do l = 1, numberl

             txt=adjustl(tfile(l))

               if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,"{") == 1 ) cycle   !skip comments and blank lines

               if (INDEX(txt,'INSTRUMENTAL') ==1) then                                    !search for sections

                    z=l+1

                  DO

                    txt = adjustl (tfile (z))

                    if(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                      z=z+1
                      cycle
                    else if (index(txt, 'STRUCTURAL')==1  .or. index(txt,'LAYER')==1 .or. index(txt,'STACKING')==1 .or. &
                             index(txt, 'TRANSITIONS')==1 .or. index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                      exit
                    else                                              !read radiation type
                      if (index(txt , 'X-RAY')/=0 ) then
                          crys%rad_type = 0
                          z=z+1
                      elseif (index(txt , 'NEUTRON')/=0) then
                          crys%rad_type = 1
                          z=z+1
                      elseif (index(txt , 'ELECTRON')/=0) then
                          crys%rad_type = 2
                          z=z+1
                      else
                          Err_crys=.true.
                          Err_crys_mess="ERROR reading radiation type"
                          logi = .false.
                          return
                      end if

                        do                                                            !read lambda
                        txt = adjustl(tfile(z))

                        if( index(txt,'{')==1 .or. len_trim(txt) == 0 .or. index(txt,"!") == 1) then
                                z=z+1
                                cycle
                        elseif (index(txt, 'STRUCTURAL')==1  .or. index(txt,'LAYER')==1 .or. &
                                index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 .or. &
                                index(txt,'CALCULATION')==1  .or. index(txt,'EXPERIMENTAL')==1) then
                                exit
                        else
                                call getword (txt, word, n_word)
                                 if (n_word == 3) then
                                   read (unit = txt, fmt = *, iostat = ier) crys%lambda    , crys%lambda2, crys%ratio
                                    if (ier /= 0)then
                                        Err_crys=.true.
                                        Err_crys_mess="ERROR reading lambda"
                                        logi = .false.
                                        return
                                    end if
                                      z = z + 1
                                 else
                                   read (unit = txt, fmt = *, iostat = ier) crys%lambda
                                    if (ier /= 0)then
                                        Err_crys=.true.
                                        Err_crys_mess="ERROR reading lambda"
                                        logi = .false.
                                        return
                                    end if
                                      z = z + 1
                                 end if
                                do

                                txt = adjustl(tfile(z))
                                if (INDEX (txt, 'NONE') /=0) then    !read instrumental data
                                       crys%broad = none
                                        z = z+1
                                elseif (index(txt, 'STRUCTURAL')==1  .or. index(txt,'LAYER')==1 .or. &
                                        index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 .or. &
                                        index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                    exit
                                else if( index(txt,'{')==1 .or. len_trim(txt) == 0 .or. index(txt,"!") == 1) then
                                    z=z+1
                                    cycle
                                else if (INDEX (txt, 'GAUSSIAN') /=0) then                     !we need to know number of parameters to distinguish them
                                      z = z+1
                                    if ((index(txt,'{') /= 0)) then
                                       j=index(txt,'{')-1
                                       if( j > 1) then
                                         call getword (txt(1:j), word, n_word)
                                       else
                                         n_word=0
                                       end if
                                    else
                                       call getword (txt, word, n_word)
                                    end if
                                    if (iErr_fmt/=0) then
                                       Err_string=.true.
                                       Err_string_mess="Error reading number of words"
                                       logi = .false.
                                       return
                                    elseif (n_word == 2) then
                                       crys%broad = pv_gss
                                       read(unit=txt ,fmt=*, iostat = ier) broad,  crys%fwhm
                                               if ( crys%fwhm <= 0)  crys%broad = none            ! if 0, equivalent to NONE
                                               if ( crys%fwhm < 0) then
                                                  Err_crys=.true.
                                                  Err_crys_mess="ERROR :  fwhm cannot be negative"
                                                  logi = .false.
                                                  return
                                               end if

                                    elseif (n_word == 4) then
                                       crys%broad = gauss

                                       read(unit=txt ,fmt=*, iostat = ier) broad, crys%p_u,  crys%p_v,  crys%p_w

                                    elseif (n_word == 5) then
                                       crys%broad = gauss

                                       read(unit=txt ,fmt=*, iostat = ier) broad,  crys%p_u,  crys%p_v,  crys%p_w, trima
                                       crys%trm = .true.

                                    else
                                       Err_crys=.true.
                                       Err_crys_mess="ERROR reading instrumental broadening"
                                       logi = .false.
                                       return
                                    end if
                                else if (INDEX (txt, 'LORENTZIAN') /=0) then
                                          z = z+1
                                    if ((index(txt,'{')/= 0)) then
                                       j=index(txt,'{')-1
                                       if( j > 1) then
                                         call getword (txt(1:j), word, n_word)
                                       else
                                         n_word=0
                                       end if
                                    else
                                       call getword (txt, word, n_word)
                                    end if
                                    if (err_string) then
                                       Err_crys=.true.
                                       Err_crys_mess="ERROR reading number of words"
                                       logi = .false.
                                       return
                                    elseif (n_word == 2) then
                                       crys%broad = pv_lrn

                                       read(unit=txt ,fmt=*, iostat = ier) broad,  crys%fwhm
                                       if ( crys%fwhm <= eps )  crys%broad = none
                                       if ( crys%fwhm < 0) then
                                              Err_crys=.true.
                                              Err_crys_mess="ERROR :  fwhm cannot be negative"
                                              logi = .false.
                                              return
                                       end if

                                    elseif (n_word == 4) then
                                       crys%broad = lorenz
                                       read(unit=txt ,fmt=*, iostat = ier) broad,  crys%p_u,  crys%p_v,  crys%p_w

                                    elseif (n_word == 5) then
                                       crys%broad = lorenz
                                       read(unit=txt ,fmt=*, iostat = ier) broad,  crys%p_u,  crys%p_v,  crys%p_w,trima
                                       crys%trm = .true.

                                    else
                                       Err_crys=.true.
                                       Err_crys_mess="ERROR reading instrumental broadening"
                                       logi = .false.
                                       return
                                    end if
                                elseif (INDEX (txt, 'PSEUDO-VOIGT')/=0) then

                                    crys%broad = ps_vgt
                                    if ((index(txt,'{')/= 0)) then
                                       j=index(txt,'{')-1
                                       if( j > 1) then
                                         call getword (txt(1:j), word, n_word)
                                       else
                                         n_word=0
                                       end if
                                    else
                                       call getword (txt, word, n_word)
                                    end if

                                    if (Err_string) then
                                       err_crys = .true.
                                       Err_crys_mess="ERROR reading number of words"
                                       logi = .false.
                                       return
                                    elseif (n_word == 7) then
                                       read(unit=txt ,fmt=*, iostat = ier) broad,  crys%p_u,  crys%p_v,  crys%p_w,  crys%p_x  , &
                                                                           crys%p_dg, crys%p_dl

                                           z = z + 1

                                             do
                                                 txt = adjustl(tfile (z))

                                                if (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                                    z = z +1
                                                     cycle

                                                ELSE IF ( index(txt,'STRUCTURAL')==1 .OR. index(txt,'LAYER')==1 .or. &
                                                          index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 .or. &
                                                          index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                                     exit

                                                else

                                                  y = index (txt, '(')    ! check if there are parenthesis and eliminate them
                                                    if (y /= 0 ) then
                                                      q = index (txt, ')')
                                                      txt (y:y) = ' '
                                                      txt (q:q) = ' '
                                                    end if
                                                   if (index (txt, '{')/= 0) then        ! comments present
                                                      j=index(txt,'{')-1
                                                      if( j > 1) then
                                                        call getword (txt(1:j), word, n_word)
                                                      else
                                                        n_word=0
                                                      end if
                                                   else
                                                      call getword (txt, word, n_word)
                                                   end if

                                                   if (n_word == 6 ) then

                                                        read (unit = txt, fmt = *) crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, &
                                                                                   crys%ref_p_x,  crys%ref_p_dg,  crys%ref_p_dl

                                                        if (  crys%ref_p_u /= 0 .or. crys%ref_p_v /= 0 .or. crys%ref_p_w /= 0 &
                                                              .or. crys%ref_p_x /= 0 .or. crys%ref_p_dg /= 0   &
                                                              .or. crys%ref_p_dl /= 0  ) then

                                                            Err_crys=.true.
                                                            Err_crys_mess=&
                                                            "ERROR :  Range of refinement missing for pseudo_voigt parameters "
                                                            logi = .false.
                                                            return
                                                         end if
                                                          crys%rang_p_u= 0.0
                                                          crys%rang_p_v= 0.0
                                                          crys%rang_p_w= 0.0
                                                          crys%rang_p_x= 0.0
                                                          crys%rang_p_dg= 0.0
                                                          crys%rang_p_dl= 0.0

                                                          z = z +1



                                                   elseif (n_word == 12) then
                                                      read (unit = txt, fmt = *)crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w,    &
                                                                                crys%ref_p_x, crys%ref_p_dg,  crys%ref_p_dl

                                                      if (crys%ref_p_u /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_u
                                                         write(unit=namepar(crys%npar),fmt="(a)")'u'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_u)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_u))-(10.0*  &
                                                                             REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_u ))
                                                         crys%vlim1 (crys%npar) = crys%p_u- crys%rang_p_u
                                                         if (crys%vlim1(crys%npar) .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_u + crys%rang_p_u
                                                      end if
                                                      if (crys%ref_p_v /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_v
                                                         write(unit=namepar(crys%npar),fmt="(a)")'v'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_v)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_v))-(10.0*  &
                                                              REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_v ))
                                                         crys%vlim1 (crys%npar) = crys%p_v- crys%rang_p_v
                                                         crys%vlim2 (crys%npar) = crys%p_v + crys%rang_p_v
                                                         if (crys%vlim2(crys%npar)  .GT. 0 ) crys%vlim2(crys%npar) = 0
                                                      end if
                                                      if (crys%ref_p_w /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_w
                                                         write(unit=namepar(crys%npar),fmt="(a)")'w'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_w)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_w))-(10.0* &
                                                               REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_w ))
                                                         crys%vlim1 (crys%npar) = crys%p_w- crys%rang_p_w
                                                         if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_w + crys%rang_p_w
                                                      end if
                                                      if (crys%ref_p_x /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_x
                                                         write(unit=namepar(crys%npar),fmt="(a)")'x'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_x)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_x))-(10.0* &
                                                                REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_x ))
                                                         crys%vlim1 (crys%npar) = crys%p_x- crys%rang_p_x
                                                         if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0

                                                         crys%vlim2 (crys%npar) = crys%p_x + crys%rang_p_x
                                                      end if
                                                      if (crys%ref_p_dg /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_dg
                                                         write(unit=namepar(crys%npar),fmt="(a)")'Dg'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_dg)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_dg))-(10.0* &
                                                                 REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_dg ))
                                                         crys%vlim1 (crys%npar) = crys%p_dg- crys%rang_p_dg
                                                         if (crys%vlim1(crys%npar)   .LT. 0 ) then
                                                                crys%vlim1(crys%npar) = 0

                                                         end if
                                                         crys%vlim2 (crys%npar) = crys%p_dg + crys%rang_p_dg
                                                      end if
                                                      if (crys%ref_p_dl /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_dl
                                                         write(unit=namepar(crys%npar),fmt="(a)")'Dl'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_dl)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_dl))-(10.0* &
                                                                 REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_dl ))
                                                         crys%vlim1 (crys%npar) = crys%p_dl - crys%rang_p_dl
                                                         if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_dl + crys%rang_p_dl
                                                      end if
                                                     if (crys%ref_p_u == 0) crys%rang_p_u = 0.0
                                                     if (crys%ref_p_v == 0) crys%rang_p_v = 0.0
                                                     if (crys%ref_p_w == 0) crys%rang_p_w = 0.0
                                                     if (crys%ref_p_x == 0) crys%rang_p_x = 0.0
                                                     if (crys%ref_p_dg == 0) crys%rang_p_dg = 0.0
                                                     if (crys%ref_p_dl == 0) crys%rang_p_dl = 0.0

                                                      z = z +1


                                                   else
                                                      Err_crys=.true.
                                                      Err_crys_mess="ERROR reading pseudo_voigt refinement parameters"
                                                      logi = .false.
                                                      return
                                                   end if

                                                end if


                                            end do


                                    elseif (n_word == 8) then
                                       read(unit=txt ,fmt=*, iostat = ier) broad,  crys%p_u,  crys%p_v,  crys%p_w,  crys%p_x , &
                                                                           crys%p_dg, crys%p_dl  , trima

                                       crys%trm = .true.
                                        z = z + 1

                                             do
                                                 txt = adjustl(tfile (z))

                                                if (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                                    z = z +1
                                                     cycle

                                                ELSE IF ( index(txt,'STRUCTURAL')==1 .or. index(txt,'LAYER')==1 .or. &
                                                          index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 .or. &
                                                          index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                                     exit

                                                else

                                                  y = index (txt, '(')    ! check if there are parenthesis and eliminate them
                                                    if (y /= 0 ) then
                                                      q = index (txt, ')')
                                                      txt (y:y) = ' '
                                                      txt (q:q) = ' '
                                                    end if
                                                   if (index (txt, '{')/= 0) then        ! comments present
                                                      j=index(txt,'{')-1
                                                      if( j > 1) then
                                                        call getword (txt(1:j), word, n_word)
                                                      else
                                                        n_word=0
                                                      end if
                                                    else
                                                      call getword (txt, word, n_word)
                                                   end if

                                                   if (n_word == 6 ) then

                                                        read (unit = txt, fmt = *) crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w,  &
                                                                                   crys%ref_p_x,  crys%ref_p_dg,  crys%ref_p_dl

                                                        if (  crys%ref_p_u /= 0 .or. crys%ref_p_v /= 0 .or. crys%ref_p_w /= 0 &
                                                              .or. crys%ref_p_x /= 0 .or. crys%ref_p_dg /= 0 &
                                                              .or. crys%ref_p_dl /= 0  ) then
                                                            Err_crys=.true.
                                                            Err_crys_mess=&
                                                             "ERROR :  Range of refinement missing for pseudo_voigt parameters "
                                                            logi = .false.
                                                            return
                                                         end if
                                                          crys%rang_p_u= 0.0
                                                          crys%rang_p_v= 0.0
                                                          crys%rang_p_w= 0.0
                                                          crys%rang_p_x= 0.0
                                                          crys%rang_p_dg= 0.0
                                                          crys%rang_p_dl= 0.0

                                                          z = z +1

                                                   elseif (n_word == 12) then
                                                      read (unit = txt, fmt = *) crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w,  &
                                                      crys%ref_p_x,  crys%ref_p_dg,  crys%ref_p_dl ,  crys%rang_p_u, &
                                                      crys%rang_p_v, crys%rang_p_w, crys%rang_p_x, crys%rang_p_dg, crys%rang_p_dl


                                                      if (crys%ref_p_u /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_u
                                                         write(unit=namepar(crys%npar),fmt="(a)")'u'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_u)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_u))-(10.0* &
                                                                 REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_u ))
                                                         crys%vlim1 (crys%npar) = crys%p_u- crys%rang_p_u
                                                         if (crys%vlim1(crys%npar) .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_u + crys%rang_p_u
                                                      end if
                                                      if (crys%ref_p_v /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_v
                                                         write(unit=namepar(crys%npar),fmt="(a)")'v'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_v)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_v))-(10.0* &
                                                                 REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_v ))
                                                         crys%vlim1 (crys%npar) = crys%p_v- crys%rang_p_v
                                                         if (crys%vlim1(crys%npar) .GT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_v + crys%rang_p_v
                                                      end if
                                                      if (crys%ref_p_w /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_w
                                                         write(unit=namepar(crys%npar),fmt="(a)")'w'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_w)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_w))-(10.0* &
                                                                 REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_w ))
                                                         crys%vlim1 (crys%npar) = crys%p_w- crys%rang_p_w
                                                         if (crys%vlim1(crys%npar) .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_w + crys%rang_p_w
                                                      end if
                                                      if (crys%ref_p_x /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_x
                                                         write(unit=namepar(crys%npar),fmt="(a)")'x'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_x)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_x))-(10.0* &
                                                                 REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_x ))
                                                         crys%vlim1 (crys%npar) = crys%p_x- crys%rang_p_x
                                                         if (crys%vlim1(crys%npar) .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_x + crys%rang_p_x
                                                      end if
                                                      if (crys%ref_p_dg /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_dg
                                                         write(unit=namepar(crys%npar),fmt="(a)")'Dg'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_dg)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_dg))-(10.0* &
                                                                 REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_dg ))
                                                         crys%vlim1 (crys%npar) = crys%p_dg- crys%rang_p_dg
                                                         if (crys%vlim1(crys%npar) .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_dg + crys%rang_p_dg
                                                      end if
                                                      if (crys%ref_p_dl /= 0) then
                                                         crys%npar = crys%npar + 1        ! to count npar
                                                         crys%list (crys%npar) =crys%p_dl
                                                         write(unit=namepar(crys%npar),fmt="(a)")'Dl'
                                                         crys%cod(crys%npar) = 1
                                                         ab =  (abs(crys%ref_p_dl)/10.0)
                                                         crys%p(crys%npar)= int(ab)
                                                         crys%mult(crys%npar) = ((abs(crys%ref_p_dl))-(10.0* &
                                                                 REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_p_dl ))
                                                         crys%vlim1 (crys%npar) = crys%p_dl - crys%rang_p_dl
                                                         if (crys%vlim1(crys%npar) .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                         crys%vlim2 (crys%npar) = crys%p_dl + crys%rang_p_dl
                                                      end if
                                                     if (crys%ref_p_u == 0) crys%rang_p_u = 0.0
                                                     if (crys%ref_p_v == 0) crys%rang_p_v = 0.0
                                                     if (crys%ref_p_w == 0) crys%rang_p_w = 0.0
                                                     if (crys%ref_p_x == 0) crys%rang_p_x = 0.0
                                                     if (crys%ref_p_dg == 0) crys%rang_p_dg = 0.0
                                                     if (crys%ref_p_dl == 0) crys%rang_p_dl = 0.0

                                                      z = z +1


                                                   else
                                                      Err_crys=.true.
                                                      Err_crys_mess="ERROR reading pseudo_voigt refinement parameters"
                                                      logi = .false.
                                                      return
                                                   end if

                                                end if


                                            end do

                                    else
                                       Err_crys=.true.
                                       Err_crys_mess="ERROR reading instrumental broadening"
                                       logi = .false.
                                       return
                                    end if

                                elseif( len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                       exit

                                else
                                    Err_crys=.true.
                                    Err_crys_mess="ERROR reading instrumental parameters"
                                    logi = .false.
                                    return
                                end if
                                end do

                        end if
                        end do

                    end if
                  END DO


               elseif (INDEX(txt,'STRUCTURAL') ==1 ) then      !structural section

                  z = l+1
                  do

                    txt = adjustl (tfile (z))

                    if(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                z=z+1
                                cycle
                    elseif ( index(txt,'LAYER')==1 .or. index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 .or. &
                     index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                exit
                    else

                      read(unit=txt ,fmt=*, iostat = ier) crys%cell_a, crys%cell_b, crys%cell_c, crys%cell_gamma   !read cell parameters
                        if (ier /= 0)then
                          Err_crys=.true.
                          Err_crys_mess="ERROR reading cell parameters"
                          logi = .false.
                          return
                        end if
                        z = z +1
                        crys%cell_gamma = crys%cell_gamma * deg2rad

                        DO

                          txt = adjustl (tfile (z))
                          IF(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                z=z+1
                                cycle
                          ELSEIF ( index(txt,'LAYER')==1 .or. index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 &
                                .or. index(txt,'CALCULATION')==1  .or. index(txt,'EXPERIMENTAL')==1) then
                                exit

                          ELSE
                            k = index (txt, '(')
                            if (k /= 0 ) then              !eliminate parenthesis
                               m = index (txt, ')')
                               txt (k:k) = ' '
                               txt (m:m) = ' '
                            end if
                            if (index (txt, '{')/= 0) then        ! comments present
                               j=index(txt,'{')-1
                               if( j > 1) then
                                 call getword (txt(1:j), word, n_word)
                               else
                                 n_word=0
                               end if
                            else
                               call getword (txt, word, n_word)
                            end if

                            if (n_word == 4 ) then

                               read(unit=txt ,fmt=*, iostat = ier) crys%ref_cell_a, crys%ref_cell_b, crys%ref_cell_c, &
                                                                   crys%ref_cell_gamma   !read cell parameters

                               if (crys%ref_cell_a /= 0.0 .or. crys%ref_cell_b /= 0.0 .or. crys%ref_cell_c /= 0.0   &
                                   .or. crys%ref_cell_gamma /= 0.0 )then
                                   Err_crys=.true.
                                   Err_crys_mess="ERROR :  Range of refinement missing for cell parameters"
                                   logi = .false.
                                   return
                               end if
                                 z = z +1
                               crys%rang_cell_a = 0.0
                               crys%rang_cell_b = 0.0
                               crys%rang_cell_c = 0.0
                               crys%rang_cell_gamma = 0.0

                            elseif (n_word == 8) then

                              read(unit=txt ,fmt=*, iostat = ier) crys%ref_cell_a, crys%ref_cell_b, crys%ref_cell_c, &
                              crys%ref_cell_gamma, crys%rang_cell_a, crys%rang_cell_b, crys%rang_cell_c, crys%rang_cell_gamma   !read cell parameters
                              if (ier /= 0)then
                                  Err_crys=.true.
                                  Err_crys_mess="ERROR reading cell refinement parameters"
                                  logi = .false.
                                  return
                              end if
                                 z = z + 1


                              if (crys%ref_cell_a /= 0 ) then
                                 crys%npar = crys%npar + 1       !to count npar
                                 crys%list (crys%npar) = crys%cell_a
                                 crys%cod(crys%npar) = 1
                                 ab =  (abs(crys%ref_cell_a)/10.0)
                                 crys%p(crys%npar)= INT(ab)
                                 crys%mult(crys%npar) = ((abs(crys%ref_cell_a))-(10.0* &
                                         REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_cell_a))
                                 namepar(crys%npar) = 'cell_a'
                                 crys%vlim1 (crys%npar) = crys%cell_a - crys%rang_cell_a
                                 if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0
                                 crys%vlim2 (crys%npar) = crys%cell_a + crys%rang_cell_a
                              end if
                              if (crys%ref_cell_b /= 0 ) then
                                 crys%npar = crys%npar + 1
                                 crys%list (crys%npar) = crys%cell_b
                                 namepar(crys%npar) = 'cell_b'
                                 crys%cod(crys%npar) = 1
                                 ab =  (abs(crys%ref_cell_b)/10.0)
                                 crys%p(crys%npar)= int(ab)
                                 crys%mult(crys%npar) = ((abs(crys%ref_cell_b))-(10.0* &
                                         REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_cell_b))
                                 crys%vlim1 (crys%npar) = crys%cell_b - crys%rang_cell_b
                                 if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0
                                 crys%vlim2 (crys%npar) = crys%cell_b + crys%rang_cell_b
                              end if
                              if (crys%ref_cell_c /= 0 ) then
                                 crys%npar = crys%npar + 1
                                 crys%list (crys%npar) = crys%cell_c
                                 namepar(crys%npar) = 'cell_c'
                                 crys%cod(crys%npar) =1
                                 ab =  (abs(crys%ref_cell_c)/10.0)
                                 crys%p(crys%npar)= int(ab)
                                 crys%mult(crys%npar) = ((abs(crys%ref_cell_c))-(10.0* &
                                         REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_cell_c))
                                 crys%vlim1 (crys%npar) = crys%cell_c - crys%rang_cell_c
                                 if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0
                                 crys%vlim2 (crys%npar) = crys%cell_c + crys%rang_cell_c
                              end if
                              if (crys%ref_cell_gamma /= 0 ) then
                                 crys%npar = crys%npar + 1
                                 crys%list (crys%npar) = crys%cell_gamma
                                 namepar(crys%npar) = 'cell_gamma'
                                 crys%cod(crys%npar) = 1
                                 ab =  (abs(crys%ref_cell_gamma)/10.0)
                                 crys%p(crys%npar)= int(ab)
                                 crys%mult(crys%npar) = ((abs(crys%ref_cell_gamma))-(10.0* &
                                         REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_cell_gamma))
                                 crys%vlim1 (crys%npar) = (crys%cell_gamma)- (crys%rang_cell_gamma * deg2rad)
                                 if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0
                                 crys%vlim2 (crys%npar) = crys%cell_gamma + (crys%rang_cell_gamma * deg2rad)
                              end if

                              if (crys%ref_cell_a == 0 ) crys%rang_cell_a = 0.0
                              if (crys%ref_cell_b == 0 ) crys%rang_cell_b = 0.0
                              if (crys%ref_cell_c == 0 ) crys%rang_cell_c = 0.0
                              if (crys%ref_cell_gamma == 0 ) crys%rang_cell_gamma = 0.0

                            else
                                Err_crys=.true.
                                Err_crys_mess="ERROR reading cell refinement parameters"
                                logi = .false.
                                return
                            end if

                            DO

                             txt = adjustl (tfile (z))

                             IF(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                z=z+1
                                cycle
                             ELSE IF ( index(txt,'LAYER')==1 .or. index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 &
                                  .or. index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                exit

                             ELSE

                                if (index (txt, '{')/= 0) then        ! comments present
                                   j=index(txt,'{')-1
                                   if( j > 1) then
                                     call getword (txt(1:j), word, n_word)
                                   else
                                     n_word=0
                                   end if
                                else
                                   call getword (txt, word, n_word)
                                end if
                                if (n_word == 2) then
                                   read(unit=txt ,fmt=*, iostat = ier) crys%sym  , crys%tolerance           ! read symmetry
                                   if (ier /= 0)then
                                      Err_crys=.true.
                                      Err_crys_mess="ERROR reading symmetry"
                                      logi = .false.
                                      return
                                   end if
                                   z=z+1
                                else

                                   read(unit=txt ,fmt=*, iostat = ier) crys%sym         ! read symmetry
                                   if (ier /= 0)then
                                     Err_crys=.true.
                                     Err_crys_mess="ERROR reading symmetry"
                                     logi = .false.
                                     return
                                   end if
                                   z=z+1
                                   crys%tolerance = 1
                                end if

                                list(1)  = '-1 '
                                list(2)  = '2/M(1) '
                                list(3)  = '2/M(2) '
                                list(4)  = 'MMM '
                                list(5)  = '-3 '
                                list(6)  = '-3M '
                                list(7)  = '4/M '
                                list(8)  = '4/MMM '
                                list(9)  = '6/M '
                                list(10) = '6/MMM '
                                list(11) = 'AXIAL '
                                list(12) = 'UNKNOWN '

                                iflag = choice(crys%sym, list, 12)
                                crys%sym  = list(iflag)

                                IF(iflag >= 1 .AND. iflag <= 11) THEN
                                     crys%symgrpno = iflag
                                ELSE IF(iflag == 12) THEN
                                    ! a value of UNKNOWN = -1 alerts OPTIMZ that user entered 'UNKNOWN'
                                     crys%symgrpno = -1
                                END IF

                                crys%tolerance = crys%tolerance * eps2   ! convert from a percentage

                                do
                                  txt = adjustl (tfile (z))

                                  IF(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                    z=z+1
                                    cycle
                                  ELSE IF ( index(txt,'LAYER')==1 .or. index(txt,'STACKING')==1 .or. &
                                            index(txt,'TRANSITIONS')==1 .or. index(txt,'CALCULATION')==1 .or. &
                                            index(txt,'EXPERIMENTAL')==1) then
                                    exit

                                  ELSE
                                      read(unit=txt ,fmt=*, iostat = ier) crys%n_typ                  !read number of  layers
                                         if (ier /= 0)then
                                           Err_crys=.true.
                                           Err_crys_mess="ERROR reading number of layers"
                                           logi = .false.
                                           return
                                         end if
                                         z = z + 1

                                     DO  i = 1, crys%n_typ        ! initialization of upper and lower positions
                                        crys%high_atom(i) = zero
                                        crys%low_atom(i)  = zero
                                     END DO

                                     do
                                         txt = adjustl (tfile (z))
                                         IF(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                              z=z+1
                                              cycle
                                         ELSE IF ( index(txt,'LAYER')==1 .or. index(txt,'STACKING')==1 .or. &
                                                   index(txt,'TRANSITIONS')==1 .or. index(txt,'CALCULATION')==1 .or. &
                                                   index(txt,'EXPERIMENTAL')==1) then
                                              exit

                                         ELSE
                                              if ((index(txt,'{')/= 0)) then                               !need to know number of words to distinguish between diameter of width along a
                                                 j=index(txt,'{')-1
                                                 if( j > 1) then
                                                   call getword (txt(1:j), word, n_word)
                                                 else
                                                   n_word=0
                                                 end if
                                              else
                                                 call getword (txt, word, n_word)
                                              end if
                                              if (iErr_fmt/=0) then
                                                 Err_string=.true.
                                                 Err_string_mess="Error reading number of words"
                                                 logi = .false.
                                                 return
                                              elseif (n_word == 1) then
                                                  if (INDEX (txt, 'INFINITE') /=0) then
                                                          crys%finite_width = .false.
                                                          z = z +1
                                                  elseif ( INDEX (txt, 'LAYER') /=0) then            ! if next section, infinite is considered
                                                          crys%finite_width = .false.
                                                          backspace 0
                                                  ! elseif (len_trim(txt) == 0 ) then
                                                  !         crys%finite_width = .false.
                                                  else
                                                     crys%finite_width = .true.
                                                     backspace 0
                                                     read (unit =tfile, fmt = *, iostat = ier) crys%layer_a
                                                     crys%layer_b = crys%layer_a
                                                     z = z + 1

                                                     do
                                                        txt = adjustl(tfile (z))

                                                        if (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. &
                                                            index(txt,"!") == 1) then
                                                            z = z +1
                                                             cycle

                                                        ELSE IF ( index(txt,'LAYER')==1 .or. &
                                                                  index(txt,'STACKING')==1 .or. &
                                                                  index(txt,'TRANSITIONS')==1 .or. &
                                                                  index(txt,'CALCULATION')==1 .or. &
                                                                  index(txt,'EXPERIMENTAL')==1) then
                                                             exit

                                                        else

                                                          y = index (txt, '(')    ! check if there are parenthesis and eliminate them
                                                            if (y /= 0 ) then
                                                              q = index (txt, ')')
                                                              txt (y:y) = ' '
                                                              txt (q:q) = ' '
                                                           end if
                                                           if (index (txt, '{')/= 0) then        ! comments present
                                                             j=index(txt,'{')-1
                                                             if( j > 1) then
                                                               call getword (txt(1:j), word, n_word)
                                                             else
                                                               n_word=0
                                                             end if
                                                           else
                                                              call getword (txt, word, n_word)
                                                           end if

                                                           if (n_word == 1 ) then

                                                               read (unit = txt, fmt = *) crys%ref_layer_a
                                                                       crys%ref_layer_b = crys%ref_layer_a
                                                               if (  crys%ref_layer_a /= 0  ) then
                                                                   Err_crys=.true.
                                                                   Err_crys_mess=&
                                                                   "ERROR :  Range of refinement missing for layer diameter "
                                                                   logi = .false.
                                                                   return
                                                                end if
                                                                 crys%rang_layer_a= 0.0
                                                                 crys%rang_layer_b= 0.0
                                                                 z = z +1
                                                           elseif (n_word == 2) then
                                                              read (unit = txt, fmt = *) crys%ref_layer_a,crys%rang_layer_a
                                                               crys%ref_layer_b = crys%ref_layer_a
                                                              if (crys%ref_layer_a /= 0) then
                                                                 crys%npar = crys%npar + 1        ! to count npar
                                                                 crys%list (crys%npar) =crys%layer_a
                                                                 write(unit=namepar(crys%npar),fmt="(a)")'diameter_a'
                                                                 crys%cod(crys%npar) = 1
                                                                 ab =  (abs(crys%ref_layer_a)/10.0)
                                                                 crys%p(crys%npar)= int(ab)
                                                                 crys%mult(crys%npar) = ((abs(crys%ref_layer_a ))-(10.0* &
                                                                         REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_layer_a ))
                                                                 crys%vlim1 (crys%npar) = crys%layer_a - crys%rang_layer_a
                                                                 crys%vlim2 (crys%npar) = crys%layer_a + crys%rang_layer_a

                                                                 crys%npar = crys%npar + 1
                                                                 crys%list (crys%npar) =crys%layer_b
                                                                 write(unit=namepar(crys%npar),fmt="(a)")'diameter_b'
                                                                 crys%cod(crys%npar) = 1
                                                                 ab =  (abs(crys%ref_layer_a)/10.0)
                                                                 crys%p(crys%npar)= int(ab)
                                                                 crys%mult(crys%npar) = ((abs(crys%ref_layer_a ))-(10.0* &
                                                                         REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_layer_a ))
                                                                 crys%vlim1 (crys%npar) = crys%layer_a - crys%rang_layer_a
                                                                 crys%vlim2 (crys%npar) = crys%layer_a + crys%rang_layer_a
                                                             end if



                                                             if (crys%ref_layer_a == 0) then
                                                                    crys%rang_layer_a  = 0.0
                                                                    crys%rang_layer_b  = 0.0
                                                             end if

                                                              z = z +1


                                                           else
                                                              Err_crys=.true.
                                                             Err_crys_mess= &
                                                              "ERROR reading layer diameter refinement parameters"
                                                              logi = .false.
                                                              return
                                                           end if

                                                        end if


                                                     end do
                                                  end if

                                              elseif (n_word == 2) then
                                                      crys%finite_width = .true.

                                                      read(unit=txt ,fmt=*, iostat = ier) crys%layer_a, crys%layer_b

                                                      if (ier /= 0 ) then
                                                             Err_crys=.true.
                                                             Err_crys_mess="ERROR reading a profile  file: end of file1"
                                                             logi = .false.
                                                             return
                                                      end if
                                                      z = z + 1

!*******************************************************************************************************************
                                                      do
                                                          txt = adjustl(tfile (z))

                                                         if (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. &
                                                             index(txt,"!") == 1) then
                                                             z = z +1
                                                              cycle

                                                         ELSE IF ( index(txt,'LAYER')==1 .or. &
                                                                   index(txt,'STACKING')==1 .or. &
                                                                   index(txt,'TRANSITIONS')==1 .or. &
                                                                   index(txt,'CALCULATION')==1 .or. &
                                                                   index(txt,'EXPERIMENTAL')==1) then
                                                              exit

                                                         else

                                                           y = index (txt, '(')    ! check if there are parenthesis and eliminate them
                                                             if (y /= 0 ) then
                                                               q = index (txt, ')')
                                                               txt (y:y) = ' '
                                                               txt (q:q) = ' '
                                                             end if
                                                            if (index (txt, '{')/= 0) then        ! comments present
                                                              j=index(txt,'{')-1
                                                              if( j > 1) then
                                                                call getword (txt(1:j), word, n_word)
                                                              else
                                                                n_word=0
                                                              end if
                                                            else
                                                               call getword (txt, word, n_word)
                                                            end if

                                                           if (n_word == 2 ) then

                                                                read (unit = txt, fmt = *) crys%ref_layer_a , &
                                                                                           crys%ref_layer_b

                                                                if (  crys%ref_layer_a /= 0  ) then
                                                                    Err_crys=.true.
                                                                    Err_crys_mess= &
                                                                     "ERROR :  Range of refinement missing for layer dimensions "
                                                                    logi = .false.
                                                                    return
                                                                 end if
                                                                  crys%rang_layer_a= 0.0
                                                                  crys%rang_layer_b= 0.0
                                                                  z = z +1
                                                           else if (n_word == 4) then
                                                             read (unit = txt, fmt = *) crys%ref_layer_a,&
                                                             crys%ref_layer_b, crys%rang_layer_a , crys%rang_layer_b

                                                             if (crys%ref_layer_a /= 0) then
                                                               crys%npar = crys%npar + 1        ! to count npar
                                                               crys%list (crys%npar) =crys%layer_a
                                                               write(unit=namepar(crys%npar),fmt="(a)")'diameter_a'
                                                               crys%cod(crys%npar) = 1
                                                               ab =  (abs(crys%ref_layer_a)/10.0)
                                                               crys%p(crys%npar)= int(ab)
                                                               crys%mult(crys%npar) = ((abs(crys%ref_layer_a ))-(10.0* &
                                                                       REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_layer_a ))
                                                               crys%vlim1 (crys%npar) = crys%layer_a - crys%rang_layer_a
                                                               crys%vlim2 (crys%npar) = crys%layer_a + crys%rang_layer_a
                                                             end if
                                                             if (crys%ref_layer_a /= 0) then
                                                              crys%npar = crys%npar + 1
                                                              crys%list (crys%npar) =crys%layer_b
                                                              write(unit=namepar(crys%npar),fmt="(a)")'diameter_b'
                                                              crys%cod(crys%npar) = 1
                                                              ab =  (abs(crys%ref_layer_b)/10.0)
                                                              crys%p(crys%npar)= int(ab)
                                                              crys%mult(crys%npar) = ((abs(crys%ref_layer_b ))-(10.0* &
                                                                     REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_layer_b ))
                                                              crys%vlim1 (crys%npar) = crys%layer_b - crys%rang_layer_b
                                                              crys%vlim2 (crys%npar) = crys%layer_b + crys%rang_layer_b
                                                             end if



                                                             if (crys%ref_layer_a == 0) crys%rang_layer_a  = 0.0
                                                             if (crys%ref_layer_b == 0) crys%rang_layer_b  = 0.0

                                                                     z = z +1


                                                           else
                                                               Err_crys=.true.
                                                               Err_crys_mess="ERROR reading layer dimensions refinement parameters"
                                                               logi = .false.
                                                               return
                                                           end if
                                                         end if


                                                      end do
!*************************************************************************************************************

                                             else
                                                      Err_crys=.true.
                                                      Err_crys_mess="ERROR reading a profile  file: end of file2"
                                                      logi = .false.
                                                      return
                                             end if

                                         end if

                                     end do

                                  END IF
                                end do
                             END IF

                            END DO

                          END IF

                        END DO
                    end if
                  end do

               elseif  (INDEX(txt,'LAYER') == 1 ) then      !Cannot be /=0 because layer can appear in comments!!
                           k = 0
                           m = 0
                           r = 0      ! counts n_actual
                           j = 0      ! counts l_actual
                   do
                         txt = adjustl(tfile(z))

                         if (index(txt,'=') /= 0) then        ! search for '=' sign
                            j = j+1
                            z = z +1
                            read (unit = txt, fmt =*, iostat = ier) layer, a1, layer, a2
                            if (a2 > a1 .OR. a2 < 1) then
                               write (*,*) "Layer", a1, " cannot be equal to layer" ,a2, "."
                               logi = .false.
                               return
                            else
                              crys%l_actual (j) = crys%l_actual (a2)
                            end if
                            cycle

                         elseif (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                             z = z +1
                             cycle

                         elseif (index(txt, 'STRUCTURAL')== 1 .or.  index(txt,'STACKING')==1 .or. &
                                 index(txt,'TRANSITIONS')==1 .or. index(txt,'CALCULATION')==1 .or. &
                                 index(txt,'EXPERIMENTAL')==1) then

                             exit
                         else
                             if  (INDEX(txt,'LAYER') == 1 ) then
                                        z = z + 1
                                        cycle
                             end if
                             j = j+1
                             r = r + 1
                             crys%l_actual(j) = j
                             crys%n_actual = r

                             if    (INDEX(txt, 'CENTROSYMMETRIC') == 1) then
                                crys%centro(r) = CENTRO
                             else
                                crys%centro(r) = NONE
                             end if
                             z = z +1
                             yy=k
                             k = 0      !to count  number atoms

                             do

                                txt = adjustl(tfile(z))

                                if (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                    z = z + 1
                                    cycle
                                elseif ( index (txt, '{') == 1 ) then
                                    z = z + 1
                                    cycle
                                elseif (index(txt,'LAYER')==1 .or. index(txt, 'STRUCTURAL')== 1 .or.  &
                                        index(txt,'STACKING')==1 .or. index(txt,'CALCULATION')==1  .or. &
                                        index(txt,'TRANSITIONS')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                    exit
                                else
                                    if (index(txt,'{')/= 0) then         ! if line contents comments
                                            k = k+1
                                      read (unit = txt, fmt = '(a4i5)')   crys%a_name(k,j), crys%a_num(k,j)
                                            pos1 = index(txt,'.')         ! localize real and/or fractional numbers
                                            pos2 = index(txt,'/')

                                      if (pos1 /= 0 .and. ( pos1<=pos2 .or. pos2 == 0 )) then    ! first number is decimal
                                         i1= index(txt,'.')-1
                                         i2= index(txt,'{')-1
                                         if( i1 >= 1 .and. i2 >=1) then
                                           call getword (txt(i1:i2), word, n_word)
                                         else
                                           n_word=0
                                           word=txt
                                         end if
                                         call read_fraction (word,n_word, num_real)             ! we convert fractions to decimals using subroutine read_fraction

                                            crys%a_pos(1, k,j)  = num_real(1)
                                            crys%a_pos(2, k,j)  = num_real(2)
                                            crys%a_pos(3, k,j)  = num_real(3)
                                            crys%a_B (k,j)      = num_real(4)
                                            crys%a_occup(k,j)   = num_real(5)
                                            z = z + 1

                                      else if ((pos2 /= 0 .and. pos2<=pos1) .or. (pos1 == 0 .and. pos2/=0)) then  ! first number is fractional
                                         i1= index(txt,'/')-1
                                         i2= index(txt,'{')-1
                                         if( i1 >= 1 .and. i2 >=1) then
                                           call getword (txt(i1:i2), word, n_word)
                                         else
                                           n_word=0
                                           word=txt
                                         end if
                                         call read_fraction (word,n_word, num_real)

                                            crys%a_pos(1, k,j)  = num_real(1)
                                            crys%a_pos(2, k,j)  = num_real(2)
                                            crys%a_pos(3, k,j)  = num_real(3)
                                            crys%a_B (k,j)      = num_real(4)
                                            crys%a_occup(k,j)   = num_real(5)
                                            z = z + 1

                                      else
                                         Err_crys=.true.
                                         Err_crys_mess="ERROR :  Atomic positions must be real numbers"
                                         logi = .false.
                                         return
                                      end if

                                         tmp = crys%a_pos(3,k,j)                          !to asign lower and upper postions ****
                                         IF(tmp > crys%high_atom(j)) crys%high_atom(j) = tmp
                                         IF(tmp < crys%low_atom(j))   crys%low_atom(j) = tmp

                                         do
                                            txt = adjustl(tfile (z))

                                            if (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then


                                                    z = z +1
                                                   cycle

                                            elseif ( index(txt,'LAYER')==1 .or. index(txt, 'STRUCTURAL')== 1 .or.  &
                                                     index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 .or. &
                                                     index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                                exit

                                            else

                                               y = index (txt, '(')    ! check if there are parenthesis and eliminate them
                                               if (y /= 0 ) then
                                                  q = index (txt, ')')
                                                  txt (y:y) = ' '
                                                  txt (q:q) = ' '
                                               end if
                                               if (index (txt, '{')/= 0) then        ! comments present
                                                  j=index(txt,'{')-1
                                                  if( j > 1) then
                                                    call getword (txt(1:j), word, n_word)
                                                  else
                                                    n_word=0
                                                  end if
                                               else
                                                  call getword (txt, word, n_word)
                                               end if

                                               if (n_word == 4 ) then

                                                  read (unit = txt, fmt = *) crys%ref_a_pos(1, k,j),crys%ref_a_pos(2, k,j),&
                                                                             crys%ref_a_pos(3, k,j) , crys%ref_a_B(k,j)
                                                  if (  crys%ref_a_pos(1, k,j)/= 0 .or.  crys%ref_a_pos(2, k,j)/=0 .or.    &
                                                        crys%ref_a_pos(3, k,j)/= 0 ) then
                                                           Err_crys=.true.
                                                           Err_crys_mess=&
                                                           "ERROR :  Range of refinement missing for atomic positions "
                                                           logi = .false.
                                                           return
                                                  end if
                                                  crys%rang_a_pos(1, k,j)= 0.0
                                                  crys%rang_a_pos(2, k,j)= 0.0
                                                  crys%rang_a_pos(3, k,j)= 0.0
                                                  crys%rang_a_B(k,j) = 0.0
                                                  z = z +1
                                                  exit
                                               elseif (n_word == 8) then
                                                  read (unit = txt, fmt = *) crys%ref_a_pos(1, k,j),crys%ref_a_pos(2, k,j),&
                                                                             crys%ref_a_pos(3, k,j),crys%ref_a_B(k,j), &
                                                                             crys%rang_a_pos(1, k,j), crys%rang_a_pos(2, k,j), &
                                                                             crys%rang_a_pos(3, k,j) , crys%rang_a_B(k,j)

                                                  if (crys%ref_a_pos(1, k,j)/= 0) then
                                                           crys%npar = crys%npar + 1        ! to count npar
                                                           crys%list (crys%npar) =crys%a_pos(1, k,j)
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'pos_x',k,j   !att: solo para i y j de 1 digito!!!!
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_a_pos(1, k,j))/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_a_pos(1, k,j)))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_a_pos(1, k,j)))
                                                           crys%vlim1 (crys%npar) = crys%a_pos(1, k,j) - crys%rang_a_pos(1, k,j)
                                                           crys%vlim2 (crys%npar) = crys%a_pos(1, k,j) + crys%rang_a_pos(1, k,j)
                                                  end if
                                                  if (crys%ref_a_pos(2, k,j)/= 0) then
                                                           crys%npar = crys%npar + 1
                                                           crys%list (crys%npar) = crys%a_pos(2, k,j)
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'pos_y',k,j
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_a_pos(2, k,j))/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_a_pos(2, k,j)))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_a_pos(2, k,j)))
                                                           crys%vlim1 (crys%npar) = crys%a_pos(2, k,j) - crys%rang_a_pos(2, k,j)
                                                           crys%vlim2 (crys%npar) = crys%a_pos(2, k,j) + crys%rang_a_pos(2, k,j)
                                                  end if
                                                  if (crys%ref_a_pos(3, k,j)/= 0) then
                                                           crys%npar = crys%npar + 1
                                                           crys%list (crys%npar) = crys%a_pos(3, k,j)
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'pos_z',k,j
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_a_pos(3, k,j))/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_a_pos(3, k,j)))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_a_pos(3, k,j)))
                                                           crys%vlim1 (crys%npar) = crys%a_pos(3, k,j) - crys%rang_a_pos(3, k,j)
                                                           crys%vlim2 (crys%npar) = crys%a_pos(3, k,j) + crys%rang_a_pos(3, k,j)
                                                  end if
                                                  if (crys%ref_a_B( k,j)/= 0) then
                                                           crys%npar = crys%npar + 1
                                                           crys%list (crys%npar) = crys%a_B( k,j)
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'Biso',k,j
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_a_B(k,j))/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_a_B( k,j)))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_a_B( k,j)))
                                                           crys%vlim1 (crys%npar) = crys%a_B( k,j) - crys%rang_a_B( k,j)
                                                           if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                           crys%vlim2 (crys%npar) = crys%a_B( k,j) + crys%rang_a_B( k,j)
                                                  end if

                                                  if (crys%ref_a_pos(1, k,j)== 0) crys%rang_a_pos (1,k,j) = 0.0
                                                  if (crys%ref_a_pos(2, k,j)== 0) crys%rang_a_pos (2,k,j) = 0.0
                                                  if (crys%ref_a_pos(3, k,j)== 0) crys%rang_a_pos (3,k,j) = 0.0
                                                  if (crys%ref_a_B( k,j)== 0) crys%rang_a_B (k,j) = 0.0
                                                  z = z +1

                                                  exit
                                               elseif ( n_word == 7 ) then  ! if we are in the atomic positions line, exit loop

                                                  exit

                                               else
                                                  Err_crys=.true.
                                                  Err_crys_mess="ERROR reading atomic positions refinement parameters"
                                                  logi = .false.
                                                  return
                                               end if
                                            end if


                                         end do


                                    else             ! line doesn't contain comments
                                            k = k+1
                                       read (unit = txt, fmt = '(a4i5)')   crys%a_name(k,j), crys%a_num(k,j)

                                            pos1 = index(txt,'.')
                                            pos2 = index(txt,'/')

                                       if (pos1 /= 0 .and. ( pos1<=pos2 .or. pos2 == 0 )) then
                                         i1= index(txt,'.')-1
                                         if( i1 >= 1 ) then
                                           call getword (txt(i1:), word, n_word)
                                         else
                                           n_word=0
                                           word=txt
                                         end if

                                         call read_fraction (word,n_word, num_real)

                                            crys%a_pos(1, k,j)  = num_real(1)
                                            crys%a_pos(2, k,j)  = num_real(2)
                                            crys%a_pos(3, k,j)  = num_real(3)
                                            crys%a_B (k,j)      = num_real(4)
                                            crys%a_occup(k,j)   = num_real(5)
                                            z = z + 1

                                       elseif ((pos2 /= 0 .and. pos2<=pos1) .or. (pos1 == 0 .and. pos2/=0))  then
                                            i1= index(txt,'/')-1
                                            if( i1 >= 1 ) then
                                              call getword (txt(i1:), word, n_word)
                                            else
                                              n_word=0
                                              word=txt
                                            end if

                                            call read_fraction (word,n_word, num_real)

                                            crys%a_pos(1, k,j)  = num_real(1)
                                            crys%a_pos(2, k,j)  = num_real(2)
                                            crys%a_pos(3, k,j)  = num_real(3)
                                            crys%a_B (k,j)      = num_real(4)
                                            crys%a_occup(k,j)   = num_real(5)
                                            z = z + 1

                                       else

                                         Err_crys=.true.
                                         Err_crys_mess="ERROR :  Atomic positions must be real numbers"
                                         logi = .false.
                                         return

                                       end if

                                         tmp = crys%a_pos(3,k,j)                          !to asign lower and upper postions ***
                                         IF(tmp > crys%high_atom(j)) crys%high_atom(j) = tmp
                                         IF(tmp < crys%low_atom(j))   crys%low_atom(j) = tmp

                                         do
                                            txt = adjustl(tfile (z))
                                            if (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                                   z = z +1
                                                   cycle

                                            elseif ( index(txt,'LAYER')==1 .or. index(txt, 'STRUCTURAL')== 1 .or.  &
                                                     index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 .or. &
                                                     index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                                exit

                                            else
                                                   y = index (txt, '(')    ! check if there are parenthesis and eliminate them
                                                   if (y /= 0 ) then
                                                        q = index (txt, ')')
                                                        txt (y:y) = ' '
                                                        txt (q:q) = ' '
                                                   end if

                                                   if (index (txt, '{')/= 0) then        ! comments present
                                                         j=index(txt,'{')-1
                                                         if( j > 1) then
                                                           call getword (txt(1:j), word, n_word)
                                                         else
                                                           n_word=0
                                                         end if
                                                   else
                                                        call getword (txt, word, n_word)
                                                   end if
                                                   if (n_word == 4 ) then

                                                        read (unit = txt, fmt = *) crys%ref_a_pos(1, k,j),crys%ref_a_pos(2, k,j),&
                                                                                   crys%ref_a_pos(3, k,j), crys%ref_a_B(k,j)
                                                        if (  crys%ref_a_pos(1, k,j)/= 0 .or.  crys%ref_a_pos(2, k,j)/=0 .or.  &
                                                              crys%ref_a_pos(3, k,j)/= 0 ) then
                                                           Err_crys=.true.
                                                           Err_crys_mess=&
                                                           "ERROR :  Range of refinement missing for atomic positions "
                                                           logi = .false.
                                                           return
                                                        end if
                                                        crys%rang_a_pos(1, k,j)= 0.0
                                                        crys%rang_a_pos(2, k,j)= 0.0
                                                        crys%rang_a_pos(3, k,j)= 0.0
                                                        crys%rang_a_B(k,j) = 0.0
                                                        z = z + 1
                                                        exit
                                                   elseif (n_word == 8) then
                                                        read (unit =txt,fmt=*) crys%ref_a_pos(1, k,j),crys%ref_a_pos(2, k,j),&
                                                                               crys%ref_a_pos(3, k,j), crys%ref_a_B(k,j), &
                                                                               crys%rang_a_pos(1, k,j), crys%rang_a_pos(2, k,j),&
                                                                               crys%rang_a_pos(3, k,j), crys%rang_a_B(k,j)

                                                        if (crys%ref_a_pos(1, k,j)/= 0) then
                                                           crys%npar = crys%npar + 1        ! to count npar
                                                           crys%list (crys%npar) =crys%a_pos(1, k,j)
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'pos_x',k,j
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_a_pos(1, k,j))/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_a_pos(1, k,j)))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_a_pos(1, k,j)))
                                                           crys%vlim1 (crys%npar) = crys%a_pos(1, k,j) - crys%rang_a_pos(1, k,j)
                                                           crys%vlim2 (crys%npar) = crys%a_pos(1, k,j) + crys%rang_a_pos(1, k,j)
                                                        end if
                                                        if (crys%ref_a_pos(2, k,j)/= 0) then
                                                           crys%npar = crys%npar + 1
                                                           crys%list (crys%npar) =crys%a_pos(2, k,j)
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'pos_y',k,j
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_a_pos(2, k,j))/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_a_pos(2, k,j)))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_a_pos(2, k,j)))
                                                           crys%vlim1 (crys%npar) = crys%a_pos(2, k,j) - crys%rang_a_pos(2, k,j)
                                                           crys%vlim2 (crys%npar) = crys%a_pos(2, k,j) + crys%rang_a_pos(2, k,j)
                                                        end if
                                                        if (crys%ref_a_pos(3, k,j)/= 0) then
                                                           crys%npar = crys%npar + 1
                                                           crys%list (crys%npar) = crys%a_pos(3, k,j)
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'pos_z',k,j
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_a_pos(3, k,j))/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_a_pos(3, k,j)))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_a_pos(3, k,j)))
                                                           crys%vlim1 (crys%npar) = crys%a_pos(3, k,j) - crys%rang_a_pos(3, k,j)
                                                           crys%vlim2 (crys%npar) = crys%a_pos(3, k,j) + crys%rang_a_pos(3, k,j)
                                                        end if
                                                        if (crys%ref_a_B( k,j)/= 0) then
                                                           crys%npar = crys%npar + 1
                                                           crys%list (crys%npar) = crys%a_B( k,j)
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'Biso',k,j
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_a_B(k,j))/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_a_B( k,j)))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_a_B( k,j)))
                                                           crys%vlim1 (crys%npar) = crys%a_B( k,j) - crys%rang_a_B( k,j)
                                                           if (crys%vlim1(crys%npar)   .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                           crys%vlim2 (crys%npar) = crys%a_B( k,j) + crys%rang_a_B( k,j)
                                                        end if
                                                        if (crys%ref_a_pos(1, k,j)== 0) crys%rang_a_pos (1,k,j) = 0.0
                                                        if (crys%ref_a_pos(2, k,j)== 0) crys%rang_a_pos (2,k,j) = 0.0
                                                        if (crys%ref_a_pos(3, k,j)== 0) crys%rang_a_pos (3,k,j) = 0.0
                                                        if (crys%ref_a_B( k,j)== 0) crys%rang_a_B(k,j) = 0.0
                                                        z = z + 1
                                                        exit
                                                   elseif ( n_word == 7 ) then
                                                        exit


                                                   else
                                                        Err_crys=.true.
                                                        Err_crys_mess="ERROR reading atomic positions refinement parameters"
                                                        logi = .false.
                                                        return
                                                   end if


                                            end if
                                         end do
                                    end if
                                end if
                                crys%l_n_atoms(j) = k

                                if (k>yy) yy=k
                             end do
                         end if
                   end do
!********************************************************************************************************************


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

!********************************************************************************************************************

               elseif  (INDEX(txt,'STACKING') == 1 ) then        !stacking section


                   z = z +1
                 DO
                    txt = adjustl (tfile (z))

                    if(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                z=z+1
                                cycle
                    elseif ( index(txt,'TRANSITIONS')==1 .or. index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1  ) then
                                exit

                    elseif (INDEX(txt , 'EXPLICIT')==1) then
                            crys%xplcit = .true.
                            z=z+1
                            crys%l_cnt = 0
                            DO

                             txt = adjustl (tfile (z))
                             if(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                z=z+1
                                cycle
                             elseif ( index(txt,'TRANSITIONS')/=0 .or. index(txt,'CALCULATION')==1 .or. &
                                      index(txt,'EXPERIMENTAL')==1 ) then
                                exit
                             elseif (index(txt , 'RANDOM')==1) then
                                   read (unit = txt, fmt = *, iostat = ier) random, crys%l_cnt
                                        if   (crys%l_cnt==0) then
                                            Err_crys=.true.
                                            Err_crys_mess="ERROR :  Number of layers is missing"
                                            logi = .false.
                                            return
                                        end if
                                   crys%inf_thick = .false.
                                   rndm = .true.
                                   z=z+1
                             elseif (index(txt , 'SEMIRANDOM')==1) then
                                   read (unit = txt, fmt = *, iostat = ier) semirandom, crys%l_cnt
                                        if   (crys%l_cnt==0) then
                                            Err_crys=.true.
                                            Err_crys_mess="ERROR :  Number of layers is missing"
                                            logi = .false.
                                            return
                                        end if
                                   crys%inf_thick = .false.
                                   rndm = .true.
                                   z=z+1
                                   do
                                     txt = adjustl (tfile (z))
                                     if(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                        z=z+1
                                        cycle
                                     elseif ( index(txt,'TRANSITIONS')/=0 .or. index(txt,'CALCULATION')==1 .or. &
                                              index(txt,'EXPERIMENTAL')==1 ) then
                                        if (crys%l_seq(1) == 0) then
                                          write(*,*) "ERROR : Layer explicit sequences missing"
                                          logi = .false.
                                          return
                                        else
                                          exit
                                        end if
                                        elseif (index(txt , 'SEQ')==1) then
                                        read (unit = txt, fmt = *, iostat = ier) seq, j1, j2, l1, l2
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
                                     end if
                                        z=z+1
                                        cycle
                                   end do
                             else
                                  do
                                   txt = adjustl (tfile (z))

                                   if(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                         z=z+1
                                         cycle
                                   elseif ( index(txt,'TRANSITIONS')/=0 .or. index(txt,'CALCULATION')==1 .or. &
                                            index(txt,'EXPERIMENTAL')==1 ) then
                                         exit
                                   else

                                         crys%inf_thick = .false.
                                         call getnum(txt, rel,inte, n_int)    ! to know number of integers per line
                                         i1=crys%l_cnt+1
                                         i2=crys%l_cnt+n_int
                                         read(unit=txt,fmt=*, iostat=ier)(crys%l_seq(r), r=i1,i2)  ! to read stacking sequence
                                         r= r-1
                                         crys%l_cnt = r
                                         z=z+1

                                   end if
                                  end do
                              end if

                            END DO

                    elseif (INDEX(txt , 'RECURSIVE')==1) then

                                crys%recrsv = .true.
                                z=z+1
                           DO

                             txt = adjustl (tfile (z))

                             if(len_trim(txt) == 0 .or. index(txt,'{')==1 .or. index(txt,"!") == 1) then
                                z=z+1
                                cycle
                             elseif ( index(txt,'TRANSITIONS')/=0 .or. index(txt,'CALCULATION')==1 .or. &
                                      index(txt,'EXPERIMENTAL')==1) then
                                exit
                             elseif (index(txt , 'INFINITE')==1) then
                                crys%inf_thick = .true.
                                z=z+1
                             else
                                read(unit=txt ,fmt= * , iostat = ier)crys%l_cnt

                                crys%inf_thick = .false.
                                z=z+1
!//////////////////////////////////////////////////////////////////////////////
                                         do
                                                   txt = adjustl(tfile (z))

                                                  if (len_trim(txt) == 0 .or. index(txt,'{')==1 .or. &
                                                      index(txt,"!") == 1) then
                                                      z = z +1
                                                       cycle

                                                  elseif ( index(txt,'LAYER')==1 .or. index(txt, 'STRUCTURAL')== 1 .or. &
                                                           index(txt,'STACKING')==1 .or. index(txt,'TRANSITIONS')==1 .or. &
                                                           index(txt,'CALCULATION')==1 .or. index(txt,'EXPERIMENTAL')==1) then
                                                       exit

                                                  else
                                                      y = index (txt, '(')    ! check if there are parenthesis and eliminate them
                                                      if (y /= 0 ) then
                                                        q = index (txt, ')')
                                                        txt (y:y) = ' '
                                                        txt (q:q) = ' '
                                                      end if
                                                      if (index (txt, '{')/= 0) then        ! comments present
                                                        j=index(txt,'{')-1
                                                        if( j > 1) then
                                                          call getword (txt(1:j), word, n_word)
                                                        else
                                                          n_word=0
                                                        end if
                                                      else
                                                        call getword (txt, word, n_word)
                                                      end if
                                                      if (n_word == 1 ) then

                                                        read (unit = txt, fmt = *) crys%ref_l_cnt
                                                        if (  crys%ref_l_cnt/= 0 ) then
                                                           Err_crys=.true.
                                                           Err_crys_mess=&
                                                           "ERROR :  Range of refinement missing for number of layers in crystal "
                                                           logi = .false.
                                                           return
                                                        end if
                                                        crys%rang_l_cnt= 0.0
                                                        z = z + 1
                                                      elseif (n_word == 2) then
                                                        read (unit = txt, fmt = *) crys%ref_l_cnt, crys%rang_l_cnt

                                                        if (crys%ref_l_cnt/= 0) then
                                                           crys%npar = crys%npar + 1        ! to count npar
                                                           crys%list (crys%npar) =crys%l_cnt
                                                           write(unit=namepar(crys%npar),fmt="(a,2i1)")'num_layers'
                                                           crys%cod(crys%npar) = 1
                                                           ab =  (abs(crys%ref_l_cnt)/10.0)
                                                           crys%p(crys%npar)= int(ab)
                                                           crys%mult(crys%npar) = ((abs(crys%ref_l_cnt))-(10.0* &
                                                                   REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_l_cnt))
                                                           crys%vlim1 (crys%npar) = crys%l_cnt - crys%rang_l_cnt
                                                           crys%vlim2 (crys%npar) = crys%l_cnt + crys%rang_l_cnt
                                                        end if
                                                        if (crys%ref_l_cnt== 0) crys%rang_l_cnt = 0.0
                                                        z = z + 1
                                                      end if
                                                  end if

                                          end do

                             end if


                           END DO

                    else
                         Err_crys=.true.
                         Err_crys_mess="ERROR reading layer sequence"
                         logi = .false.
                         return
                    end if

                    end do

               elseif  (INDEX(txt,'TRANSITIONS') == 1 ) then         ! transitions section

                      z = z + 1
                      do
                           txt = adjustl (tfile (z))
                           if (len_trim (txt) == 0 .or. index(txt,'{')==1 .or. INDEX(txt,'!')==1 ) then
                                z = z+1
                                cycle
                           elseif (INDEX(txt,'CALCULATION') == 1 ) then

                                exit
                           else
                                i = 1
                               do
                                      txt = adjustl (tfile (z))
                                  if  (len_trim (txt) == 0 .or. INDEX(txt,'{')==1 .or. INDEX(txt,'!')==1  ) then
                                        z=z+1
                                        cycle
                                  elseif (INDEX(txt,'CALCULATION') == 1 ) then
                                        exit
                                  else
                                      j=1
                                    do

                                       txt = adjustl (tfile (z))

                                       if  (len_trim (txt) == 0 .or. INDEX(txt,'{')==1 .or. INDEX(txt,'!')==1 ) then
                                            z=z+1
                                            cycle
                                       elseif (INDEX(txt,'CALCULATION') == 1 ) then
                                            exit
                                       elseif (j .GT. crys%n_typ ) then
                                            exit
                                       else

                                                   k = index (txt,'(')             !  if Fats Waller parameters are present we eliminate parenthesis

                                                if (k /=  0) then
                                                    txt (k:k) = ' '
                                                    m = index (txt,')')
                                                    txt (m:m) = ' '

                                                end if

                                                if (index (txt, '{')/= 0) then        ! comments present
                                                   j=index(txt,'{')-1
                                                   if( j > 1) then
                                                     call getword (txt(1:j), word, n_word)
                                                   else
                                                     n_word=0
                                                   end if
                                                else
                                                   call getword (txt, word, n_word)
                                                end if


                                                if (n_word == 4) then               ! only probabilities and stacking vectors

                                                   call read_fraction (word,n_word, num_real)

                                                   crys%l_alpha (j,i) = num_real(1)
                                                   crys%l_r (1,j,i)   = num_real(2)
                                                   crys%l_r (2,j,i)   = num_real(3)
                                                   crys%l_r (3,j,i)   = num_real(4)

                                                   z=z+1



                                                 elseif (n_word == 10) then                          ! Fats Waller present

                                                   call read_fraction (word,n_word, num_real)

                                                   crys%l_alpha (j,i) = num_real(1)
                                                   crys%l_r (1,j,i)   = num_real(2)
                                                   crys%l_r (2,j,i)   = num_real(3)
                                                   crys%l_r (3,j,i)   = num_real(4)
                                                   crys%r_b11 (j,i)        = num_real(5)
                                                   crys%r_b22 (j,i)        = num_real(6)
                                                   crys%r_b33 (j,i)        = num_real(7)
                                                   crys%r_b12 (j,i)        = num_real(8)
                                                   crys%r_b31 (j,i)        = num_real(9)
                                                   crys%r_b23 (j,i)        = num_real(10)

                                                   z=z+1


                                                else
                                                   Err_crys=.true.
                                                   Err_crys_mess="ERROR reading layer transitions"
                                                   logi = .false.
                                                   return
                                                end if

                                                  DO
                                                    txt = adjustl (tfile (z))


                                                      IF  (len_trim (txt) == 0 .or. INDEX(txt,'{')==1 .or. INDEX(txt,'!')==1 ) then
                                                                 z=z+1
                                                                 cycle
                                                      ELSEIF (INDEX(txt,'CALCULATION') == 1 ) then
                                                                 exit

                                                      ELSE

                                                            txt = adjustl (tfile (z))
                                                            k = index (txt, '(')    ! check if there are parenthesis and eliminate them
                                                             if (k /= 0 ) then
                                                               m = index (txt, ')')
                                                               txt (k:k) = ' '
                                                               txt (m:m) = ' '
                                                             end if

                                                             if (index (txt, '{')/= 0) then        ! comments present
                                                             	 j=index(txt,'{')-1
                                                             	 if( j > 1) then
                                                                   call getword (txt(1:j), word, n_word)
                                                                 else
                                                                   n_word=0
                                                                 end if
                                                             else
                                                                 call getword (txt, word, n_word)
                                                             end if
                                                             if (n_word == 4 ) then

                                                              read (unit = txt, fmt = *) crys%ref_l_alpha (j,i), &
                                                                crys%ref_l_r (1,j,i),crys%ref_l_r (2,j,i),crys%ref_l_r (3,j,i)
                                                              if (crys%ref_l_alpha(j,i)/=0 .or. crys%ref_l_r(1,j,i)/=0 .or. &
                                                                  crys%ref_l_r(2,j,i)/=0 .or. crys%ref_l_r(3,j,i)/=0) then
                                                                  Err_crys=.true.
                                                                  Err_crys_mess=&
                                                                  "ERROR :  Range of refinement missing for transition parameters"
                                                                  logi = .false.
                                                                  return
                                                              end if
                                                              crys%rang_l_alpha (j,i)= 0.0
                                                              crys%rang_l_r (1,j,i)  = 0.0
                                                              crys%rang_l_r (2,j,i)  = 0.0
                                                              crys%rang_l_r (3,j,i)  = 0.0
                                                               z = z + 1

                                                                 exit

                                                             elseif (n_word == 8) then

                                                              read (unit = txt, fmt = *) crys%ref_l_alpha (j,i), &
                                                                crys%ref_l_r (1,j,i),crys%ref_l_r (2,j,i),         &
                                                                crys%ref_l_r (3,j,i), crys%rang_l_alpha (j,i),     &
                                                                crys%rang_l_r(1,j,i), crys%rang_l_r(2,j,i) ,&
                                                                crys%rang_l_r(3,j,i)
                                                                z = z + 1
                                                              if (crys%ref_l_alpha (j,i) /= 0)  then
                                                                crys%npar = crys%npar + 1    !to count npar
                                                                crys%list (crys%npar) = crys%l_alpha (j,i)
                                                                write(unit=namepar(crys%npar),fmt="(a,2i1)")'alpha',i,j
                                                                crys%cod(crys%npar) = 1
                                                                ab =  (abs(crys%ref_l_alpha (j,i))/10.0)
                                                                crys%p(crys%npar)= int(ab)
                                                                crys%mult(crys%npar) = ((abs(crys%ref_l_alpha (j,i)))-(10.0* &
                                                                     REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_l_alpha (j,i)))
                                                                crys%vlim1(crys%npar)=crys%l_alpha(j,i)- crys%rang_l_alpha (j,i)
                                                                if (crys%vlim1(crys%npar)  .LT. 0 ) crys%vlim1(crys%npar) = 0
                                                                crys%vlim2(crys%npar) = crys%l_alpha(j,i)+crys%rang_l_alpha(j,i)
                                                                if (crys%vlim2(crys%npar)  > 1 ) crys%vlim2(crys%npar) = 1

                                                              end if
                                                              if (crys%ref_l_r(1,j,i) /= 0 ) then
                                                                  crys%npar = crys%npar + 1
                                                                  crys%list (crys%npar) = crys%l_r(1,j,i)
                                                                  write(unit=namepar(crys%npar),fmt="(a,2i1)")'tx',i,j
                                                                  crys%cod(crys%npar) = 1
                                                                  ab =  (abs(crys%ref_l_r(1,j,i))/10.0)
                                                                  crys%p(crys%npar)= int(ab)
                                                                  crys%mult(crys%npar) = ((abs(crys%ref_l_r(1,j,i)))-(10.0* &
                                                                          REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_l_r(1,j,i)))
                                                                  crys%vlim1(crys%npar) = crys%l_r (1,j,i) - crys%rang_l_r (1,j,i)
                                                                  crys%vlim2(crys%npar) = crys%l_r (1,j,i) + crys%rang_l_r (1,j,i)
                                                              end if
                                                              if (crys%ref_l_r(2,j,i) /= 0 ) then
                                                                  crys%npar = crys%npar + 1
                                                                  crys%list (crys%npar) = crys%l_r(2,j,i)
                                                                  write(unit=namepar(crys%npar),fmt="(a,2i1)")'ty',i,j
                                                                  crys%cod(crys%npar) = 1
                                                                  ab =  (abs(crys%ref_l_r(2,j,i))/10.0)
                                                                  crys%p(crys%npar)= int(ab)
                                                                  crys%mult(crys%npar) = ((abs(crys%ref_l_r(2,j,i)))-(10.0* &
                                                                  REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_l_r(2,j,i)))
                                                                  crys%vlim1(crys%npar) = crys%l_r (2,j,i) - crys%rang_l_r (2,j,i)
                                                                  crys%vlim2(crys%npar) = crys%l_r (2,j,i) + crys%rang_l_r (2,j,i)
                                                              end if
                                                              if (crys%ref_l_r(3,j,i) /= 0 )  then
                                                                  crys%npar = crys%npar + 1
                                                                  crys%list (crys%npar) = crys%l_r(3,j,i)
                                                                  write(unit=namepar(crys%npar),fmt="(a,2i1)")'tz',i,j
                                                                  crys%cod(crys%npar) = 1
                                                                  ab =  (abs(crys%ref_l_r(3,j,i))/10.0)
                                                                  crys%p(crys%npar)= int(ab)
                                                                  crys%mult(crys%npar) = ((abs(crys%ref_l_r(3,j,i)))-(10.0* &
                                                                  REAL(crys%p(crys%npar))))*SIGN(1.0,(crys%ref_l_r(3,j,i)))
                                                                  crys%vlim1(crys%npar) = crys%l_r (3,j,i) - crys%rang_l_r (3,j,i)
                                                                  crys%vlim2(crys%npar) = crys%l_r (3,j,i) + crys%rang_l_r (3,j,i)
                                                              end if

                                                              if (crys%ref_l_alpha (j,i) == 0)  crys%rang_l_alpha(j,i) = 0.0
                                                              if (crys%ref_l_r(1,j,i) == 0 )  crys%rang_l_r(1,j,i) = 0.0
                                                              if (crys%ref_l_r(2,j,i) == 0 )  crys%rang_l_r(2,j,i) = 0.0
                                                              if (crys%ref_l_r(3,j,i) == 0 )  crys%rang_l_r(3,j,i) = 0.0
                                                              exit
                                                             else
                                                                Err_crys=.true.
                                                                Err_crys_mess=&
                                                                "ERROR reading stacking transitions refinement parameters"
                                                                logi = .false.
                                                                return
                                                             end if


                                                      END IF
                                                  END DO


                                                     txt = adjustl(tfile(z))

                                       end if

                                       j = j + 1

                                    end do

                                 end if
                                     i = i+1
                               end do

                           end if
                      end do

                           DO  i = 1, crys%n_typ        !check if l_apha sums one for each layer
                                sum = 0
                                DO  j = 1, crys%n_typ

                                    sum = sum + crys%l_alpha(j,i)
                                END DO

                                IF(ABS(sum - 1) > eps) THEN
                                 write(op,*) "Stacking probabilities from LAYER ",i," do not sum to 1."
                                 write (*,*)  'the sum is'  , sum
                                 logi = .false.
                                end if
                           END DO

                           do

                            if (INDEX(txt,'CALCULATION') == 1 ) then

                                exit
                            else
                               z = z +1
                               txt = adjustl(tfile(z))

                            end if
                           end do

               elseif  (INDEX(txt,'CALCULATION') == 1 ) then

                        z=l+1
                 do

                    txt = adjustl (tfile (z))

                    if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                      z = z+1
                      cycle
                    else if (index (txt, 'EXPERIMENTAL' ) == 1 .or. z > numberl) then

                      exit

                    else
                      read (unit = txt, fmt = *, iostat=ier) crys%calctype


                    select case (crys%calctype)

                      case ('SIMULATION')


                        opt = 0
                        z = z +1

                        do

                          txt = adjustl (tfile (z))

                          if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                                  z = z+1
                                  cycle
                          else if ( z > numberl ) then

                                  exit
                          else

                            read (unit = txt, fmt = *, iostat=ier) th2_min, th2_max, d_theta
                            z = z +1
                               IF(th2_min < zero) THEN
                                   WRITE(op,"(a)") ' ERROR: 2theta min is negative.'
                                   logi = .false.
                                   RETURN
                               END IF

                               IF(th2_min >= one_eighty) THEN
                                   WRITE(op,"(a)") ' ERROR: 2theta min exceeds 180 degrees.'
                                   logi = .false.
                                   RETURN
                               END IF

                               IF(th2_max <= zero) THEN
                                   WRITE(op,"(a)") ' ERROR: Negative (or zero) value for 2theta min.'
                                   logi = .false.
                                   RETURN
                               END IF

                               IF(th2_max > one_eighty) THEN
                                   WRITE(op,"(a)") ' ERROR: 2theta max exceeds 180 degrees.'
                                   logi = .false.
                                   RETURN
                               END IF

                             ! if th2_max = 180, reduce it slightly to keep us out of trouble
                               IF(th2_max == one_eighty) th2_max = th2_max - eps4

                               IF(d_theta <= zero) THEN
                                   WRITE(op,"(a)") ' ERROR: Negative (or zero) value for 2theta increment.'
                                   logi = .false.
                                   RETURN
                               END IF

                               IF(th2_min >= th2_max) THEN
                                   WRITE(op,"(a)") ' ERROR: 2theta min is larger than 2theta max.'
                                   logi = .false.
                                   RETURN
                               END IF

                               IF((th2_max - th2_min) < d_theta) THEN
                                    WRITE(op,"(a)") ' ERROR: 2theta increment is larger than 2theta max - 2theta min.'
                                    logi = .false.
                                    RETURN
                               END IF

                              th2_min = th2_min * deg2rad
                              th2_max = th2_max * deg2rad
                              d_theta = half * deg2rad * d_theta

                          end if
                        end do


                      case ('SIMPLEX')


                        opt = 2
                        z = z +1
                        opti%iquad = 1

                       do

                         txt = adjustl (tfile (z))

                         if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then

                             z=z+1
                             cycle
                         else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                             exit
                         else
                             read(unit = txt, fmt = *, iostat=ier) opti%mxfun
                             z = z +1
                             txt = adjustl (tfile (z))

                             do

                                 if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or.index(txt,'{')==1) then
                                   z=z+1
                                   txt = adjustl (tfile (z))
                                   cycle
                                 else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                                   exit
                                 else
                                   read(unit = txt, fmt = *, iostat=ier) opti%eps
                                   z = z + 1
                                   txt = adjustl (tfile (z))
                                   do

                                     if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                                       z=z+1
                                       txt = adjustl (tfile (z))
                                       cycle
                                     else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                                       exit
                                     else

                                       read(unit = txt, fmt = *, iostat=ier) opti%iout
                                       z = z + 1
                                       txt = adjustl (tfile (z))

                                       do

                                          if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                                              z=z+1
                                              txt = adjustl (tfile (z))
                                              cycle
                                          else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                                              exit
                                          else

                                              read(unit = txt, fmt = *, iostat=ier) opti%acc
                                              z = z + 1
                                              txt = adjustl (tfile (z))


                                          end if
                                       end do
                                     end if
                                   end do
                                 end if
                             end do
                         end if
                       end do



!///////////////////////CONVERSION TO SIMPLEX VARIABLES///////////////////////////////////////////

                    vector(1:nrp)  = crys%list (1:nrp)
                    !   code(1:nrp) = crys%cod(1:nrp)
                            numpar = crys%npar
                    !      nm_cycl  = crys%n_cycles
!////////////////////////////////////////////////////////////////////////////////////////////////
                        n_plex = maxval(crys%p)
                        opti%loops = 2 * n_plex
                        opti%iquad=1


                        do i = 1, numpar
                           label(crys%p(i)) = 0
                        end do

                        do i=1, numpar                    !construction of SIMPLEX 'in' vectors

                          if (label (crys%p(i))  == 0 ) then        !not to overwrite in config
                             v_plex(crys%p(i)) = crys%list(i)
                           !  nampar(pnum(i)) = namepar(i)
                           !  vlim1(pnum(i)) = crys%vlim1(i)
                           !  vlim2(pnum(i)) = crys%vlim2(i)
                          end if
                          label (crys%p(i)) = 1
                        end do

                      case ('SAN')
                       opt = 1
                       z = z +1

                       do

                         txt = adjustl (tfile (z))

                         if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then

                             z=z+1
                             cycle
                         else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                             exit
                         else

                             read(unit = txt, fmt = *, iostat=ier) crys%n_cycles
                             z = z +1
                             txt = adjustl (tfile (z))

                             do

                                 if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or.index(txt,'{')==1) then
                                   z=z+1
                                   txt = adjustl (tfile (z))
                                   cycle
                                 else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                                   exit
                                 else
                                   read(unit = txt, fmt = *, iostat=ier) crys%num_temp
                                   z = z + 1
                                   txt = adjustl (tfile (z))
                                   do

                                     if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                                       z=z+1
                                       txt = adjustl (tfile (z))
                                       cycle
                                     else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                                       exit
                                     else

                                       read(unit = txt, fmt = *, iostat=ier) crys%t_ini
                                       z = z + 1
                                       txt = adjustl (tfile (z))
                                       do

                                          if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                                              z=z+1
                                              txt = adjustl (tfile (z))
                                              cycle
                                          else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                                              exit
                                          else

                                              read(unit = txt, fmt = *, iostat=ier) crys%anneal
                                              z = z + 1
                                              txt = adjustl (tfile (z))
                                              do

                                                if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                                                  z=z+1
                                                  txt = adjustl (tfile (z))
                                                  cycle
                                                else if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
                                                  exit
                                                else

                                                  read(unit = txt, fmt = *, iostat=ier) crys%accept
                                                  z = z + 1
                                                  txt = adjustl (tfile (z))
                                                  do

                                                     if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then

                                                        z=z+1
                                                        txt = adjustl (tfile (z))
                                                        cycle
                                                     else if ( index (txt, 'EXPERIMENTAL') == 1 ) then
                                                        exit
                                                     else

                                                        read(unit = txt, fmt = *, iostat=ier) crys%init_config
                                                        z = z + 1
                                                        txt = adjustl (tfile (z))
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

                       opsan%num_therm = 0
                       opsan%nalgor = 0
                       sanvec%bound(1:crys%npar) = 0

!///////////////////////CONVERSION TO SAN VARIABLES//////////////////////////
                    sanvec%state(1:nrp)  = crys%list (1:nrp)
                      sanvec%code(1:nrp) = crys%cod(1:nrp)
                             sanvec%npar = crys%npar
                        opsan%num_temps  = crys%num_temp
                             opsan%t_ini = crys%t_ini
                            opsan%anneal = crys%anneal
                           opsan%accept  = crys%accept
                        opsan%initconfig = crys%init_config
                          opsan%nm_cycl  = crys%n_cycles

!////////////////////////////////////////////////////////////////////////////////////////////////

                         npar = maxval(crys%p)
                         numpar = crys%npar
                         label(crys%p(1:numpar)) = 0




                        do i=1, numpar                    !construction of SAN 'in' vectors
                          if(label (crys%p(i)) == 0 ) then   !not to overwrite in config
                             sanvec%config(crys%p(i)) = crys%list(i)
                             sanvec%nampar(crys%p(i)) = namepar(i)
                             sanvec%low(crys%p(i)) = crys%vlim1(i)
                             sanvec%high(crys%p(i)) = crys%vlim2(i)
                          end if
                          label (crys%p(i)) = 1
                        end do

                      case default
                        write(op,*) 'calculation type not valid'
                        logi = .false.
                    end select

                    end if
                 end do

               else if(INDEX(txt,'EXPERIMENTAL') == 1 ) then
                 z=l+1
                 do
                    txt = adjustl (tfile (z))
                    if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                      z = z+1
                      cycle
                    else if ( z > numberl) then
                      exit
                    else
                      read (unit = txt, fmt = *, iostat=ier) dfile , th2_min, th2_max, d_theta
                      z = z +1

                      if (th2_min /= 0 .and.  th2_max /= 0  .and. d_theta /= 0) then
                        th2_min = th2_min * deg2rad
                        th2_max = th2_max * deg2rad
                        d_theta = half * deg2rad * d_theta
                      end if

                      do
                         txt = adjustl (tfile (z))
                         if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                           z = z+1
                           cycle
                         elseif ( z > numberl) then
                           exit
                         else
                           read (unit = txt, fmt = *, iostat=ier) fmode
                           z = z + 1
                           do
                             txt = adjustl (tfile (z))
                             if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                                   z = z+1
                                   cycle
                             elseif ( z > numberl) then
                                   exit
                             else
                               read (unit = txt, fmt = *, iostat=ier)  background_file
                               z = z +1
                               do
                                  txt = adjustl (tfile (z))
                                  if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) then
                                    z = z+1
                                     cycle
                                  else if ( z > numberl) then
                                    exit
                                  else
                                    read (unit = txt, fmt = *, iostat=ier) mode
                                    z = z + 1
                                  end if
                               end do
                             end if
                           end do
                         end if
                      end do
                    end if
                 end do
               else
                    if(index(txt,"!") == 1  .or. len_trim(txt) == 0 .or. index(txt,'{')==1) cycle
               end if
        End do


!        end  of reading section

!///////////////////////CONVERSION TO DIFFAX VARIABLES////////////////

                              rad_type = crys%rad_type
                              lambda   = crys%lambda
                              lambda2  = crys%lambda2
                                 ratio = crys%ratio
                              blurring = crys%broad      ! Diffax utiliza varias variables, primero asigna blurring a: none, gauss, pv_gss, lorenz, pv_lrn , ps_vgt; y luego asigna cada una de estas variables a los enteros  : 0,1,4,2,5,y 3 respectivamente. Todas estas variables estan definidas como enteros
                              cell_a   = crys%cell_a
                              cell_b   = crys%cell_b
                              cell_c   = crys%cell_c
                            cell_gamma = crys%cell_gamma
                             pnt_grp   = crys%sym        ! hay otra variable: symgrpno
                            l_symmetry = crys%centro     ! mismo problema k en blurring: 1o asigna l_symmetry a none o centro y luego asigna estos a 0 y 1 respectivamente
           a_name(1:yy,1:crys%n_typ)   = crys%a_name(1:yy,1:crys%n_typ)
           a_number(1:yy,1:crys%n_typ) = crys%a_num(1:yy,1:crys%n_typ)
          a_pos(1:3,1:yy,1:crys%n_typ) = crys%a_pos(1:3,1:yy,1:crys%n_typ)
             a_b(1:yy,1:crys%n_typ)    = crys%a_B(1:yy,1:crys%n_typ)
            a_occup(1:yy,1:crys%n_typ) = crys%a_occup(1:yy,1:crys%n_typ)
                              recrsv   = crys%recrsv
                              xplcit   = crys%xplcit
                          finite_width = crys%finite_width
                            inf_thick  = crys%inf_thick
  l_alpha(1:crys%n_typ,1:crys%n_typ)   = crys%l_alpha(1:crys%n_typ,1:crys%n_typ)
  l_r(1:3,1:crys%n_typ,1:crys%n_typ)   = crys%l_r(1:3,1:crys%n_typ,1:crys%n_typ)
                                Wa     = crys%layer_a
                                Wb     = crys%layer_b
                               pv_u    = crys%p_u
                               pv_v    = crys%p_v
                               pv_w    = crys%p_w
                               pv_x    = crys%p_x
                               pv_dg   = crys%p_dg
                               pv_dl   = crys%p_dl
                               fwhm    = crys%FWHM
                              pv_gamma = crys%p_gamma
                           trim_origin = crys%trm
                            n_actual   = crys%n_actual
                              l_cnt    = crys%l_cnt
                       l_seq(1:xp_max) = crys%l_seq(1:xp_max)
                l_actual(1:crys%n_typ) = crys%l_actual(1:crys%n_typ)
                           SymGrpNo    = crys%SymGrpNo
              r_b11(1:yy,1:crys%n_typ) = crys%r_b11 (1:yy,1:crys%n_typ)
              r_b22(1:yy,1:crys%n_typ) = crys%r_b22 (1:yy,1:crys%n_typ)
              r_b33(1:yy,1:crys%n_typ) = crys%r_b33 (1:yy,1:crys%n_typ)
              r_b12(1:yy,1:crys%n_typ) = crys%r_b12 (1:yy,1:crys%n_typ)
              r_b31(1:yy,1:crys%n_typ) = crys%r_b31 (1:yy,1:crys%n_typ)
              r_b23(1:yy,1:crys%n_typ) = crys%r_b23 (1:yy,1:crys%n_typ)
              l_n_atoms(1:crys%n_typ)  = crys%l_n_atoms(1:crys%n_typ)
               low_atom(1:crys%n_typ)  = crys%low_atom(1:crys%n_typ)
              high_atom(1:crys%n_typ)  = crys%high_atom(1:crys%n_typ)
                             tolerance = crys%tolerance
                                  pnum = crys%p(1:numpar)
                                  mult = crys%mult(1:numpar)
                              n_layers = crys%n_typ



        open (unit=read_out, file='read.out', status='replace', action='write')
        write(read_out, FMT = "(A30,1X,3(F15.4,1X))")   "lambda , lambda2 , ratio",  lambda , lambda2 , ratio
        write(read_out, FMT = "(A30,1X, 6(F15.4,1X))")   "instrumental parameters", pv_u,  pv_v, pv_w,  pv_x , pv_dg, pv_dl
        write(read_out, FMT = "(A30,1X,4(F15.6,1X))") "cell parameters", cell_a, cell_b, cell_c, cell_gamma
        write(read_out, FMT = "(A30,1X,I4)") "number of layer types", n_layers
        write(read_out, FMT = "(A30,1X,I4)") "number of unique layers", n_actual
        write(read_out, FMT = "(20X,A70)") "layer number   Atom name   number   x      y      z         Biso      Occ "
        write(read_out, FMT = "(20X,A70)") "-----------------------------------------------------------------------------"
        do i=1, n_actual
          do j=1,l_n_atoms(i)
            write(*,*) "n_actual, l_n_atoms",i,j
            write(read_out, FMT = "(25X, I4, 4X, A4, 4X, I4, 4X,3(F5.4, 4X), 2(F5.4, 4X))") i, a_name(i,j), a_number(i,j), a_pos(1:3,i,j), a_b(i,j), a_occup(i,j)
          end do
        end do
        close (unit=read_out)
        write(*,*) "n_actual",  n_actual
        do i=1, n_layers
          do j=1, n_layers
             write(*,*) " i,j,l_alpha(i,j)" ,i,j,l_alpha(i,j)
          end do
        end do
        return
        End subroutine  read_structure_file
!_____________________________________________________________________________________________________________

        Subroutine read_fraction (line, n_comp, real_num)

        character(len=*),dimension(:), intent(in out) :: line        !string of numbers to convert to real
        integer,                       intent(in    ) :: n_comp      !number of components in line
        real            ,dimension(:), intent(   out) :: real_num    ! string of real characters

        integer                                       :: i, j
        real                                          :: numerator, denominator

        do i = 1, n_comp
            j = index (line(i), '/')
            if(j == 0) then                 !!it was already a number
              read ( line(i), fmt=*) real_num(i)
            else     !!it was expressed as a ratio. delete the slash, '/'.
              line(i)(j:j) = ' '
              read(unit =line(i),fmt = *) numerator, denominator  !! now contains two arguments, numerator and denominator.
              real_num(i) = numerator/ denominator
            end if
        end do

        End subroutine read_fraction

    End module read_data
