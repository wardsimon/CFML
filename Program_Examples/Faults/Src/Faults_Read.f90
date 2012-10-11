    Module read_data
       use CFML_GlobalDeps,           only : sp , cp
       use CFML_String_Utilities,     only : number_lines , reading_lines , findfmt ,iErr_fmt, getword, err_string, &
                                             err_string_mess, getnum, Ucase, cutst
       use CFML_Optimization_General, only : Opt_Conditions_Type
       use CFML_LSQ_TypeDef,          only : LSQ_Conditions_type, LSQ_State_Vector_Type
       use CFML_Simulated_Annealing
       use diffax_mod
       implicit none

       private
       !public subroutines
       public:: read_structure_file, new_getfil , length

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
       real                                             :: zero_shift = 0.0
       real                                             :: sycos = 0.0
       real                                             :: sysin = 0.0
       real                                             :: p_u = 0.0                    !u  > pv_u
       real                                             :: p_v = 0.0                    !v  > pv_v
       real                                             :: p_w = 0.0                    !w  > pv_w
       real                                             :: fwhm = 0.0                   !fwmh ( cannot be negative)  >fwhm
       real                                             :: p_gamma = 0.0                !pseudo-voigt gamma    > pv_gamma
       real                                             :: p_x                          ! x > pv_x
       real                                             :: p_dg, p_dl                   ! gaussian and lorentzian average volumetric size
       real                                             :: ref_p_u, ref_p_v, ref_p_w, ref_p_x, ref_p_dg, ref_p_dl
       real                                             :: ref_zero_shift, ref_sycos, ref_sysin
       real                                             :: rang_p_v, rang_p_u, rang_p_w, rang_p_x, rang_p_dg, rang_p_dl
       real                                             :: rang_zero_shift, rang_sycos, rang_sysin
       real                                             :: tolerance                    !d-> Maximum deviation that intensities can have
       real                                             :: rang_cell_a ,rang_cell_b,rang_cell_c,rang_cell_gamma
       real                                             :: t_ini=0.0, anneal, accept        ! In SAN: initial temperature, kirkpatrick factor of annealing, minimum number of accepted configurations
       character (len=132)                              :: sym                          !Laue symmetry  >pnt_grp
       character (len=132)                              :: calctype                     !type of calculation
       real                                             :: ref_cell_a,ref_cell_b, ref_cell_c, ref_cell_gamma  ! index of refinement: if zero, no refinement, if one, refinement to be done
       integer                                          :: n_typ                        !number of layer types    >n_layers
       real                                             :: l_cnt = 0.0                  !number of layers in stacking section
       real                                             :: ref_l_cnt
       real                                             :: rang_l_cnt
       integer                                          :: SymGrpNo
       integer                                          :: n_cycles=0
       integer                                          :: n_actual                     !number of unique layers
       integer                                          :: npar = 0                     !total number of parameters to be refined
       integer                                          :: num_temp=0                   ! maximum number of temperatures in SAN
       integer                                          :: init_config                  ! flag to randomly initialaze SAN configuration if = 1.
       integer                                          :: yy                           !max number of atoms per layer
       logical                                          :: recrsv                       !layer sequencing recursive
       logical                                          :: xplcit                       !layer sequencing explicit
       logical                                          :: finite_width                 !layer width in plane ab
       logical                                          :: inf_thick                    !infinit stacking
       logical                                          :: trm                          !trim keyword: to supress peak at origin >trim_origin
       character(len=4), dimension(:,:),   allocatable  :: a_name                       !atom name (4characters wide)>atom_l o a_name?
       integer,          dimension(:),     allocatable  :: centro                       !layer symmetry : none or centrosymmetric  >only_real
       integer,          dimension(:),     allocatable  :: cod
       integer,          dimension(:,:),   allocatable  :: a_num                        !atom number  > a_number
       real,             dimension(:,:,:), allocatable  :: ref_a_pos
       integer,          dimension(:)  ,   allocatable  :: l_seq                        !d-> array containing the explicitly defined sequence of layers.
       integer,          dimension(:)  ,   allocatable  :: l_actual                     !Contains the layer number that layer i is structurally identical to.
       integer,          dimension(:)  ,   allocatable  :: l_n_atoms                    !number of atoms in each layer
       integer,          dimension(:),     allocatable  :: p                            ! refinable parameter index
       real,             dimension(:),     allocatable  :: mult                         ! refinable parameter multiplicator
       real,             dimension(:,:),   allocatable  :: ref_l_alpha                  ! index of refinement of l_alpha
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
       real,             dimension(:)  ,   allocatable  :: list                         !vector containing all the parameters to optimize
       real,             dimension(:),     allocatable  :: vlim1                        !Low-limit value of parameters
       real,             dimension(:),     allocatable  :: vlim2                        !Low-limit value of parameters
       real,             dimension(:)  ,   allocatable  :: Pv_refi                      !vector containing the free parameters to optimize (restrictions taken into account)
       end  type crys_2d_type

       integer, parameter ,private                      :: i_data=30
       logical,            public, save                 :: Err_crys=.false.
       character(len=120), public, save                 :: Err_crys_mess=" "

       character(len=132), dimension(:),allocatable     :: tfile     !List of lines of the input file (UPPER CASE)
       Integer                                          :: numberl   !number of lines in the input file
       Integer                                          :: np        !to count number of parameters to be refined =npar
       Integer, dimension(7)                            :: sect_indx=0 !Indices (line numbers) of the sevent sections:
                                                                       !1:INSTRUMENTAL, 2:STRUCTURAL, 3:LAYER
                                                                       !4:STACKING , 5:TRANSITIONS, 6:CALCULATION,
                                                                       !7:EXPERIMENTAL

      type (crys_2d_type),            save,  public  :: crys
      type (Opt_Conditions_Type),     save,  public  :: opti
      type (SimAnn_Conditions_type ), save,  public  :: opsan
      type (LSQ_Conditions_type),     save,  public  :: cond
      type (LSQ_State_Vector_Type),   save,  public  :: Vs

      !integer, parameter ::  nrp=80

   contains

    Subroutine Set_TFile(namef,logi)
      character(len=*),  intent(in)  :: namef
      logical,           intent(out) :: logi
      integer :: i
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
        call reading_lines(namef, numberl, tfile)   ! we 'charge' the file in tfile  so we can close the unit
      end if
      close(unit=i_data)
      do i = 1, numberl                 ! To read in case insensitive mode
        call Ucase(tfile(i))
      end do
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
        txt=adjustl(tfile(k))
        if(txt(1:12) == "INSTRUMENTAL") then
          sect_indx(1) = k
          do
            k=k+1; if(k > numberl) exit
            txt=adjustl(tfile(k))
            if(txt(1:10) == "STRUCTURAL") then
              sect_indx(2) = k
              do
                k=k+1; if(k > numberl) exit
                txt=adjustl(tfile(k))
                if(txt(1:5) == "LAYER") then
                  sect_indx(3) = k
                  !n_actual=r+1
                  n_layers=l+1
                  do
                    k=k+1; if(k > numberl) exit
                    txt=adjustl(tfile(k))
                    if(txt(1:8) == "STACKING") then
                       sect_indx(4) = k
                       do
                         k=k+1; if(k > numberl) exit
                         txt=adjustl(tfile(k))
                         if(txt(1:11) == "TRANSITIONS") then
                            sect_indx(5) = k
                            do
                              k=k+1; if(k > numberl) exit
                              txt=adjustl(tfile(k))
                              if(txt(1:11) == "CALCULATION") then
                                 sect_indx(6) = k
                                 do
                                   k=k+1; if(k > numberl) exit
                                   txt=adjustl(tfile(k))
                                   if(txt(1:12) == "EXPERIMENTAL") then
                                      sect_indx(7) = k
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
      end do global

      !write(*,"(a,7i6)") " Indices: ", sect_indx

    End Subroutine Set_Sections_index


 !1:INSTRUMENTAL, 2:STRUCTURAL, 3:LAYER
 !4:STACKING , 5:TRANSITIONS, 6:CALCULATION,
 !7:EXPERIMENTAL
    Subroutine Read_INSTRUMENTAL(logi)
      logical, intent(out) :: logi

      integer :: i,i1,i2,k, ier, l
      character(len=132) :: txt
      character(len=25)  :: key
      logical :: ok_rad, ok_wave, ok_uvw, ok_abe
      real  :: ab

      logi=.true.
      i1=sect_indx(1)
      i2=sect_indx(2)-1

      i=i1

      ok_rad=.false.; ok_wave=.false.; ok_uvw=.false. ; ok_abe=.true.
      vs%np = 0
      np=0
      do

        i=i+1; if(i > i2) exit
        txt=adjustl(tfile(i))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        key=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(key)

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
                read(unit=txt,fmt=*, iostat=ier) crys%lambda
                if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading wavelength instruction"
                  logi=.false.
                  return
                end if
             end if
             ok_wave=.true.

          Case("ABERRATIONS")

             read(unit=txt,fmt=*, iostat=ier) crys%zero_shift, crys%sycos, crys%sysin

             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading instrumental aberrations"
               logi=.false.
               return
             end if

             i=i+1
             txt=adjustl(tfile(i))
             !Reading refinement codes
             read(unit=txt,fmt=*, iostat=ier) crys%ref_zero_shift, crys%ref_sycos,  crys%ref_sysin
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the instrumental aberrations"
               logi=.false.
               return
             end if

             !Reading range of parameters
             k=index(txt,"(")
             l=index(txt,")")
             if(k /= 0) then
                read(unit=txt(k+1:l-1),fmt=*, iostat=ier) crys%rang_zero_shift,crys%rang_sycos, crys%rang_sysin
                if(ier /= 0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading the ranges of the instrumental aberrations"
                  logi=.false.
                  return
                end if
             else
               Err_crys=.true.
               Err_crys_mess="ERROR: No parameter ranges in the instrumental aberrations: these must be given!"
               logi=.false.
               return
             end if
             !Creating refinement vectors
             if (abs(crys%ref_zero_shift) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%zero_shift
               write(unit=namepar(np),fmt="(a)")'zero_shift'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_zero_shift)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_zero_shift))-(10.0*  &
                                   REAL(crys%p(np))))*SIGN(1.0,(crys%ref_zero_shift ))
               crys%vlim1(crys%p(np)) = crys%zero_shift- crys%rang_zero_shift
               if (crys%vlim1(crys%p(np)) .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%zero_shift + crys%rang_zero_shift
               crys%Pv_refi(crys%p(np)) = crys%zero_shift

             end if
             if (abs(crys%ref_sycos) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%sycos
               write(unit=namepar(np),fmt="(a)")'sycos'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_sycos)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_sycos))-(10.0*  &
               REAL(crys%p(np))))*SIGN(1.0,(crys%ref_sycos ))
               crys%vlim1(crys%p(np)) = crys%sycos- crys%rang_sycos
               crys%vlim2(crys%p(np)) = crys%sycos + crys%rang_sycos
               if (crys%vlim2(crys%p(np))  .GT. 0 ) crys%vlim2(crys%p(np)) = 0
               crys%Pv_refi(crys%p(np))=crys%sycos
             end if
             if (abs(crys%ref_sysin) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%sysin
               write(unit=namepar(np),fmt="(a)")'sysin'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_sysin)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_sysin))-(10.0* &
                     REAL(crys%p(np))))*SIGN(1.0,(crys%ref_sysin ))
               crys%vlim1(crys%p(np)) = crys%sysin- crys%rang_sysin
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%sysin + crys%rang_sysin
               crys%Pv_refi(crys%p(np))= crys%sysin
             end if

             ok_abe=.true.


          Case("PSEUDO-VOIGT")
             crys%broad=ps_vgt
             read(unit=txt,fmt=*, iostat=ier) crys%p_u, crys%p_v, crys%p_w, crys%p_x, crys%p_dg, crys%p_dl
             if(ier /= 0 ) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading Pseudo-Voigt instruction"
               logi=.false.
               return
             end if

             if(index(txt,"TRIM") /= 0) crys%trm=.true.

             i=i+1
             txt=adjustl(tfile(i))
             !Reading refinement codes
             read(unit=txt,fmt=*, iostat=ier) crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, &
                                              crys%ref_p_x,  crys%ref_p_dg,  crys%ref_p_dl
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the Pseudo-Voigt instruction"
               logi=.false.
               return
             end if

             !Reading range of parameters
             k=index(txt,"(")
             l=index(txt,")")
             if(k /= 0) then
                read(unit=txt(k+1:l-1),fmt=*, iostat=ier) crys%rang_p_u,crys%rang_p_v, crys%rang_p_w, &
                                                       crys%rang_p_x, crys%rang_p_dg, crys%rang_p_dl
                if(ier /= 0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading the ranges of the Pseudo-Voigt parameters"
                  logi=.false.
                  return
                end if

             else
               Err_crys=.true.
               Err_crys_mess="ERROR: No parameter ranges in the Pseudo-Voigt instruction: these must be given!"
               logi=.false.
               return
             end if
             !Creating refinement vectors
             if (abs(crys%ref_p_u) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_u
               write(unit=namepar(np),fmt="(a)")'u'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_u)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_u))-(10.0*  &
                                   REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_u ))
               crys%vlim1(crys%p(np)) = crys%p_u- crys%rang_p_u
               if (crys%vlim1(crys%p(np)) .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_u + crys%rang_p_u
               crys%Pv_refi(crys%p(np)) = crys%p_u

             end if
             if (abs(crys%ref_p_v) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_v
               write(unit=namepar(np),fmt="(a)")'v'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_v)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_v))-(10.0*  &
               REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_v ))
               crys%vlim1(crys%p(np)) = crys%p_v- crys%rang_p_v
               crys%vlim2(crys%p(np)) = crys%p_v + crys%rang_p_v
               if (crys%vlim2(crys%p(np))  .GT. 0 ) crys%vlim2(crys%p(np)) = 0
               crys%Pv_refi(crys%p(np))=crys%p_v
             end if
             if (abs(crys%ref_p_w) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_w
               write(unit=namepar(np),fmt="(a)")'w'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_w)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_w))-(10.0* &
                     REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_w ))
               crys%vlim1(crys%p(np)) = crys%p_w- crys%rang_p_w
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_w + crys%rang_p_w
               crys%Pv_refi(crys%p(np))= crys%p_w
             end if
             if (abs(crys%ref_p_x) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_x
               write(unit=namepar(np),fmt="(a)")'x'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_x)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_x))-(10.0* &
                      REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_x ))
               crys%vlim1(crys%p(np)) = crys%p_x- crys%rang_p_x
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_x + crys%rang_p_x
               crys%Pv_refi(crys%p(np))=crys%p_x
             end if
             if (abs(crys%ref_p_dg) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_dg
               write(unit=namepar(np),fmt="(a)")'Dg'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_dg)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_dg))-(10.0* &
                 REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_dg ))
               crys%vlim1(crys%p(np)) = crys%p_dg- crys%rang_p_dg
               if (crys%vlim1(crys%p(np))   .LT. 0 ) then
                 crys%vlim1(crys%p(np)) = 0
               end if
               crys%vlim2(crys%p(np)) = crys%p_dg + crys%rang_p_dg
               crys%Pv_refi(crys%p(np))=crys%p_dg
             end if
             if (abs(crys%ref_p_dl) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_dl
               write(unit=namepar(np),fmt="(a)")'Dl'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_dl)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_dl))-(10.0* &
                 REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_dl ))
               crys%vlim1(crys%p(np)) = crys%p_dl - crys%rang_p_dl
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_dl + crys%rang_p_dl
               crys%Pv_refi(crys%p(np)) =crys%p_dl
             end if

             if (abs(crys%ref_p_u) > 0.0) crys%rang_p_u = 0.0
             if (abs(crys%ref_p_v) > 0.0) crys%rang_p_v = 0.0
             if (abs(crys%ref_p_w) > 0.0) crys%rang_p_w = 0.0
             if (abs(crys%ref_p_x) > 0.0) crys%rang_p_x = 0.0
             if (abs(crys%ref_p_dg) > 0.0) crys%rang_p_dg = 0.0
             if (abs(crys%ref_p_dl) > 0.0) crys%rang_p_dl = 0.0

             ok_uvw=.true.

          Case("GAUSSIAN")

             crys%broad= pv_gss

             read(unit=txt,fmt=*, iostat=ier) crys%p_u, crys%p_v, crys%p_w, crys%p_dg
             crys%p_dl=100000
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
             read(unit=txt,fmt=*, iostat=ier) crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, crys%ref_p_dg
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the Gaussian/Lorentzian instruction"
               logi=.false.
               return
             end if

             !Reading range of parameters
             k=index(txt,"(")
             l=index(txt,")")
             if(k /= 0) then
                read(unit=txt(k+1:l-1),fmt=*, iostat=ier) crys%rang_p_u,crys%rang_p_v, crys%rang_p_w, crys%rang_p_dg

                if(ier /= 0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading the ranges of the Gaussian/Lorentzian parameters"
                  logi=.false.
                  return
                end if
             else
               Err_crys=.true.
               Err_crys_mess="ERROR: No parameter ranges in the Gaussian/Lorentzian instruction: these must be given!"
               logi=.false.
               return
             end if
             !Creating refinement vectors
             if (abs(crys%ref_p_u) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_u
               write(unit=namepar(np),fmt="(a)")'u'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_u)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_u))-(10.0*  &
                                   REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_u ))
               crys%vlim1(crys%p(np)) = crys%p_u- crys%rang_p_u
               if (crys%vlim1(crys%p(np)) .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_u + crys%rang_p_u
               crys%Pv_refi(crys%p(np)) = crys%p_u

             end if
             if (abs(crys%ref_p_v) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_v
               write(unit=namepar(np),fmt="(a)")'v'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_v)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_v))-(10.0*  &
               REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_v ))
               crys%vlim1(crys%p(np)) = crys%p_v- crys%rang_p_v
               crys%vlim2(crys%p(np)) = crys%p_v + crys%rang_p_v
               if (crys%vlim2(crys%p(np))  .GT. 0 ) crys%vlim2(crys%p(np)) = 0
               crys%Pv_refi(crys%p(np))=crys%p_v
             end if
             if (abs(crys%ref_p_w) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_w
               write(unit=namepar(np),fmt="(a)")'w'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_w)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_w))-(10.0* &
                     REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_w ))
               crys%vlim1(crys%p(np)) = crys%p_w- crys%rang_p_w
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_w + crys%rang_p_w
               crys%Pv_refi(crys%p(np))= crys%p_w
             end if
             if (abs(crys%ref_p_dg) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_dg
               write(unit=namepar(np),fmt="(a)")'Dg'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_dg)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_dg))-(10.0* &
                 REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_dg ))
               crys%vlim1(crys%p(np)) = crys%p_dg- crys%rang_p_dg
               if (crys%vlim1(crys%p(np))   .LT. 0 ) then
                 crys%vlim1(crys%p(np)) = 0
               end if
               crys%vlim2(crys%p(np)) = crys%p_dg + crys%rang_p_dg
               crys%Pv_refi(crys%p(np))=crys%p_dg
             end if

             if (abs(crys%ref_p_u) > 0.0) crys%rang_p_u = 0.0
             if (abs(crys%ref_p_v) > 0.0) crys%rang_p_v = 0.0
             if (abs(crys%ref_p_w) > 0.0) crys%rang_p_w = 0.0
             if (abs(crys%ref_p_dg) > 0.0) crys%rang_p_dg = 0.0


             ok_uvw=.true.

          Case("LORENTZIAN")

             crys%broad= pv_lrn

             read(unit=txt,fmt=*, iostat=ier) crys%p_u, crys%p_v, crys%p_w, crys%p_dl
             crys%p_dg=100000
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
             read(unit=txt,fmt=*, iostat=ier) crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, crys%ref_p_dl
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the Gaussian/Lorentzian instruction"
               logi=.false.
               return
             end if

             !Reading range of parameters
             k=index(txt,"(")
             l=index(txt,")")
             if(k /= 0) then
                read(unit=txt(k+1:l-1),fmt=*, iostat=ier) crys%rang_p_u,crys%rang_p_v, crys%rang_p_w , crys%rang_p_dl

                if(ier /= 0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading the ranges of the Gaussian/Lorentzian parameters"
                  logi=.false.
                  return
                end if
             else
               Err_crys=.true.
               Err_crys_mess="ERROR: No parameter ranges in the Gaussian/Lorentzian instruction: these must be given!"
               logi=.false.
               return
             end if
             !Creating refinement vectors
             if (abs(crys%ref_p_u) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_u
               write(unit=namepar(np),fmt="(a)")'u'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_u)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_u))-(10.0*  &
                                   REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_u ))
               crys%vlim1(crys%p(np)) = crys%p_u- crys%rang_p_u
               if (crys%vlim1(crys%p(np)) .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_u + crys%rang_p_u
               crys%Pv_refi(crys%p(np)) = crys%p_u

             end if
             if (abs(crys%ref_p_v) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_v
               write(unit=namepar(np),fmt="(a)")'v'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_v)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_v))-(10.0*  &
               REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_v ))
               crys%vlim1(crys%p(np)) = crys%p_v- crys%rang_p_v
               crys%vlim2(crys%p(np)) = crys%p_v + crys%rang_p_v
               if (crys%vlim2(crys%p(np))  .GT. 0 ) crys%vlim2(crys%p(np)) = 0
               crys%Pv_refi(crys%p(np))=crys%p_v
             end if
             if (abs(crys%ref_p_w) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_w
               write(unit=namepar(np),fmt="(a)")'w'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_w)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_w))-(10.0* &
                     REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_w ))
               crys%vlim1(crys%p(np)) = crys%p_w- crys%rang_p_w
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_w + crys%rang_p_w
               crys%Pv_refi(crys%p(np))= crys%p_w
             end if
             if (abs(crys%ref_p_dl) > 0.0) then
               np = np + 1        ! to count npar
               crys%list (np) =crys%p_dl
               write(unit=namepar(np),fmt="(a)")'Dl'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_p_dl)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_p_dl))-(10.0* &
                 REAL(crys%p(np))))*SIGN(1.0,(crys%ref_p_dl ))
               write(*,*)  crys%p(np) , crys%p_dl ,   crys%rang_p_dl
               crys%vlim1(crys%p(np)) = crys%p_dl - crys%rang_p_dl
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%p_dl + crys%rang_p_dl
               crys%Pv_refi(crys%p(np)) =crys%p_dl
             end if

             if (abs(crys%ref_p_u) > 0.0) crys%rang_p_u = 0.0
             if (abs(crys%ref_p_v) > 0.0) crys%rang_p_v = 0.0
             if (abs(crys%ref_p_w) > 0.0) crys%rang_p_w = 0.0
             if (abs(crys%ref_p_dl) > 0.0) crys%rang_p_dl = 0.0

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

      integer :: i,i1,i2,j,k, ier, l
      character(len=132) :: txt
      character(len=25)  :: key
      logical :: ok_cell, ok_symm, ok_nlayers, ok_lwidth
      real  :: ab

      logi=.true.
      i1=sect_indx(2)
      i2=sect_indx(3)-1

      i=i1

      crys%finite_width=.false.

      ok_cell=.false.; ok_symm=.false.; ok_nlayers=.false. ; ok_lwidth=.false.

      do

        i=i+1; if(i > i2) exit
        txt=adjustl(tfile(i))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        key=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(key)

          Case("CELL")
                read(unit=txt,fmt=*,iostat=ier) crys%cell_a, crys%cell_b, crys%cell_c, crys%cell_gamma   !read cell parameters
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
             read(unit=txt,fmt=*, iostat=ier)  crys%ref_cell_a, crys%ref_cell_b, crys%ref_cell_c, &
                                               crys%ref_cell_gamma   !read refinement codes of cell parameters
             if(ier /= 0) then
               Err_crys=.true.
               Err_crys_mess="ERROR reading the refinement codes of the cell parameters"
               logi=.false.
               return
             end if

             !Reading range of parameters
             k=index(txt,"(")
             l=index(txt,")")
             if(k /= 0) then
               read(unit=txt(k+1:l-1),fmt=*, iostat=ier) crys%rang_cell_a, crys%rang_cell_b, crys%rang_cell_c, crys%rang_cell_gamma
               if(ier /= 0) then
                 Err_crys=.true.
                 Err_crys_mess="ERROR reading the ranges of the cell parameters"
                 logi=.false.
                 return
               end if
             else
               Err_crys=.true.
               Err_crys_mess="ERROR: No parameter ranges in the cell instruction: these must be given!"
               logi=.false.
               return
             end if

             !Creating refinement vectors
             if (abs(crys%ref_cell_a) > 0.0 ) then
               np = np + 1       !to count npar
               crys%list (np) = crys%cell_a
               crys%cod(np) = 1
               ab =  (abs(crys%ref_cell_a)/10.0)
               crys%p(np)= INT(ab)
               crys%mult(np) = ((abs(crys%ref_cell_a))-(10.0* &
                       REAL(crys%p(np))))*SIGN(1.0,(crys%ref_cell_a))
               namepar(np) = 'cell_a'
               crys%vlim1(crys%p(np)) = crys%cell_a - crys%rang_cell_a
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%cell_a + crys%rang_cell_a
               crys%Pv_refi(crys%p(np)) =crys%cell_a
             end if
             if (abs(crys%ref_cell_b) > 0.0 ) then
               np = np + 1
               crys%list (np) = crys%cell_b
               namepar(np) = 'cell_b'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_cell_b)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_cell_b))-(10.0* &
                       REAL(crys%p(np))))*SIGN(1.0,(crys%ref_cell_b))
               crys%vlim1(crys%p(np)) = crys%cell_b - crys%rang_cell_b
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%cell_b + crys%rang_cell_b
               crys%Pv_refi(crys%p(np)) =crys%cell_b
             end if
             if (abs(crys%ref_cell_c) > 0.0 ) then
               np = np + 1
               crys%list (np) = crys%cell_c
               namepar(np) = 'cell_c'
               crys%cod(np) =1
               ab =  (abs(crys%ref_cell_c)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_cell_c))-(10.0* &
                       REAL(crys%p(np))))*SIGN(1.0,(crys%ref_cell_c))
               crys%vlim1(crys%p(np)) = crys%cell_c - crys%rang_cell_c
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%cell_c + crys%rang_cell_c
               crys%Pv_refi(crys%p(np)) =crys%cell_c
             end if
             if (abs(crys%ref_cell_gamma) > 0.0 ) then
               np = np + 1
               crys%list (np) = crys%cell_gamma
               namepar(np) = 'cell_gamma'
               crys%cod(np) = 1
               ab =  (abs(crys%ref_cell_gamma)/10.0)
               crys%p(np)= int(ab)
               crys%mult(np) = ((abs(crys%ref_cell_gamma))-(10.0* &
                       REAL(crys%p(np))))*SIGN(1.0,(crys%ref_cell_gamma))
               crys%vlim1(crys%p(np)) = (crys%cell_gamma)- (crys%rang_cell_gamma * deg2rad)
               if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
               crys%vlim2(crys%p(np)) = crys%cell_gamma + (crys%rang_cell_gamma * deg2rad)
               crys%Pv_refi(crys%p(np)) =crys%cell_gamma
             end if

             if (abs(crys%ref_cell_a) > 0.0 ) crys%rang_cell_a = 0.0
             if (abs(crys%ref_cell_b) > 0.0 ) crys%rang_cell_b = 0.0
             if (abs(crys%ref_cell_c) > 0.0 ) crys%rang_cell_c = 0.0
             if (abs(crys%ref_cell_gamma) > 0.0 ) crys%rang_cell_gamma = 0.0

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
             elseif ( crys%sym == "2/M(1)") then
               Crys%Symgrpno = 2
             elseif ( crys%sym == " 2/M(2)") then
               Crys%Symgrpno = 3
             elseif ( crys%sym == "MMM")  then
               Crys%Symgrpno = 4
             elseif ( crys%sym == "-3") then
               Crys%Symgrpno = 5
             elseif ( crys%sym == "-3M") then
               Crys%Symgrpno = 6
             elseif ( crys%sym == "4/M")  then
               Crys%Symgrpno = 7
             elseif ( crys%sym == "4/MMM")  then
               Crys%Symgrpno = 8
             elseif ( crys%sym == "6/M") then
               Crys%Symgrpno = 9
             elseif ( crys%sym == "6/MMM")  then
               Crys%Symgrpno = 10
             elseif ( crys%sym == "AXIAL")  then
               Crys%Symgrpno = 11
             elseif ( crys%sym == "UNKNOWN")  then
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

             ok_nlayers=.true.

          Case ("LWIDTH")

             if (INDEX (txt, 'INFINITE') /= 0) then
               return
             elseif  (INDEX (txt, 'INFINITE') == 0) then
               read(unit=txt,fmt=*, iostat=ier)  crys%layer_a, crys%layer_b
               crys%finite_width=.true.
               if (crys%layer_b==0) then
                  crys%layer_b= crys%layer_a
                  ier=0
               end if

               if(ier /= 0 ) then
                 Err_crys=.true.
                 Err_crys_mess="ERROR reading layer width parameters"
                 logi=.false.
                 return
               end if
               if (crys%layer_b==0) then
                  crys%layer_b= crys%layer_a
               end if
               i=i+1
               txt=adjustl(tfile(i))
              ! Reading refinement codes
               k=index(txt,"(")
               l=index(txt,")")
                 read(unit=txt(:k),fmt=*, iostat=ier) crys%ref_layer_a, crys%ref_layer_b
                   if (crys%ref_layer_b==0) then
                     crys%ref_layer_b = crys%ref_layer_a
                     ier=0
                   end if
                   if(ier /= 0) then
                     Err_crys=.true.
                     Err_crys_mess="ERROR reading the refinement codes of layer width parameters"
                     logi=.false.
                     return
                   end if

                   ! Reading range of parameters
                   if(k /= 0) then
                     read(unit=txt(k+1:l-1),fmt=*, iostat=ier) crys%rang_layer_a, crys%rang_layer_b

                     if (abs(crys%ref_layer_b) > 0.0) then
                       crys%rang_layer_b = crys%rang_layer_a
                       ier=0
                     end if
                     if(ier /= 0) then
                       Err_crys=.true.
                       Err_crys_mess="ERROR reading the ranges of the Pseudo-Voigt parameters"
                       logi=.false.
                       return
                     end if
                   else
                     Err_crys=.true.
                     Err_crys_mess="ERROR: No parameter ranges in the Layer width instruction: this must be given!"
                     logi=.false.
                     return
                   end if

                   !Creating refinement vectors
                   if (abs(crys%ref_layer_a) > 0.0) then
                     np = np + 1        ! to count npar
                     crys%list (np) =crys%layer_a
                     write(unit=namepar(np),fmt="(a)")'diameter_a'
                     crys%cod(np) = 1
                     ab =  (abs(crys%ref_layer_a)/10.0)
                     crys%p(np)= int(ab)
                     crys%mult(np) = ((abs(crys%ref_layer_a ))-(10.0* &
                       REAL(crys%p(np))))*SIGN(1.0,(crys%ref_layer_a ))
                     crys%vlim1(crys%p(np)) = crys%layer_a - crys%rang_layer_a
                     crys%vlim2(crys%p(np)) = crys%layer_a + crys%rang_layer_a
                     np = np + 1
                     crys%list (np) =crys%layer_b
                     crys%Pv_refi(crys%p(np)) =crys%layer_a
                     write(unit=namepar(np),fmt="(a)")'diameter_b'
                     crys%cod(np) = 1
                     ab =  (abs(crys%ref_layer_a)/10.0)
                     crys%p(np)= int(ab)
                     crys%mult(np) = ((abs(crys%ref_layer_a ))-(10.0* &
                       REAL(crys%p(np))))*SIGN(1.0,(crys%ref_layer_a ))
                     crys%vlim1(crys%p(np)) = crys%layer_a - crys%rang_layer_a
                     crys%vlim2(crys%p(np)) = crys%layer_a + crys%rang_layer_a
                     crys%Pv_refi(crys%p(np)) =crys%layer_b
                   end if

                   if (abs(crys%ref_layer_a) > 0.0) then
                     crys%rang_layer_a  = 0.0
                     crys%rang_layer_b  = 0.0
                   end if
             else
                Err_crys=.true.
                Err_crys_mess="ERROR: No parameter in the Layer width instruction: this must be given!"
                logi=.false.
                return
             end if

             ok_lwidth=.true.

          Case Default

             cycle

          End Select
      end do

      if(ok_cell .and. ok_symm .and. ok_nlayers .and. ok_lwidth) then
        return
      else
        Err_crys=.true.
        Err_crys_mess="ERROR: cell parameters, Laue symmetry or number of layer types missing!"
        logi=.false.
      end if

    End Subroutine Read_STRUCTURAL

    Subroutine Read_LAYER(logi)
      logical, intent(out) :: logi

      integer :: i,i1,i2,i3, k, ier, a, a1, a2, r, tmp, j, m, l, nitem
      character(len=20), dimension(30) :: citem
      character(len=132) :: txt, layer
      character(len=25)  :: key

      logical :: ok_lsym, ok_atom
      integer,          dimension(:),     allocatable  :: d !counts nº of atoms in unique layer
      real  :: ab

      if (allocated (d)) deallocate(d)
      allocate(d(max_a))
      d=0

      logi=.true.
      i1=sect_indx(3)
      i2=sect_indx(4)-1

      i=i1-1

      ok_lsym=.false.; ok_atom=.false.

       a = 0      ! counts l_n_atoms
       r = 0      ! counts n_actual
       a1= 0      ! reads layer number
       a2= 0      ! reads layer number
       !yy=0       ! counts nº of atoms in layer type
       j=0
       l=0
       m=0
      do
        i=i+1; if(i > i2) exit
        txt=adjustl(tfile(i))

        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        key=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(key)

          Case ("LAYER")
            if (index(txt,'=') /= 0) then        ! search for '=' sign
              j = j+1

              read (unit = txt, fmt =*, iostat = ier) a1, layer, a2
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
                crys%l_actual(j) = crys%l_actual(a2)

              end if
            else
              j=j+1
              r=r+1
              crys%l_actual(j) = r

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
            do m=1, 3

              call read_fraction(citem(2+m), crys%a_pos(m, d(r),r))
              if(ERR_String) then
                write(unit=*,fmt="(a)") trim(ERR_String_Mess)
                logi=.false.
                return
              end if
            end do
            crys%a_name(d(r),r)=citem(1)
            read(unit=citem(2),fmt=*,iostat=ier) crys%a_num (d(r),r)
            read(unit=citem(6),fmt=*,iostat=ier) crys%a_B (d(r),r)
            read(unit=citem(7),fmt=*,iostat=ier) crys%a_occup (d(r),r)

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
            read (unit = txt, fmt =*, iostat=ier) crys%ref_a_pos(1, d(r),r), crys%ref_a_pos(2, d(r),r), &
                                                  crys%ref_a_pos(3, d(r),r), crys%ref_a_B(d(r),r)
            if(ier /= 0) then
                   Err_crys=.true.
                   Err_crys_mess="ERROR reading atomic parameters refinement codes"
                   logi=.false.
                   return
            end if

            !reading range of atomic parameters
            k=index(txt,"(")
            l=index(txt,")")
            if (k/=0) then
              read (unit = txt(k+1:l-1), fmt =*, iostat=ier) crys%rang_a_pos(1,d(r),r), crys%rang_a_pos(2,d(r),r), &
                                                             crys%rang_a_pos(3,d(r),r) , crys%rang_a_B(d(r),r)
                if(ier /= 0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading the ranges of the atomic positions"
                  logi=.false.
                  return
                end if
            else
               Err_crys=.true.
               Err_crys_mess="ERROR: No parameter ranges in the atomic instruction: these must be given!"
               logi=.false.
               return
            end if

            if (abs(crys%ref_a_pos(1, d(r),r)) > 0.0) then
              np = np + 1        ! to count npar
              crys%list (np) =crys%a_pos(1, d(r),r)
              write(unit=namepar(np),fmt="(a,2i1)")'pos_x',d(r),r   !att: solo para i y j de 1 digito!!!!
              crys%cod(np) = 1
              ab =  (abs(crys%ref_a_pos(1, d(r),r))/10.0)
              crys%p(np)= int(ab)
              crys%mult(np) = ((abs(crys%ref_a_pos(1, d(r),r)))-(10.0* &
              REAL(crys%p(np))))*SIGN(1.0,(crys%ref_a_pos(1,d(r),r)))
              crys%vlim1(crys%p(np)) = crys%a_pos(1, d(r),r) - crys%rang_a_pos(1, d(r),r)
              crys%vlim2(crys%p(np)) = crys%a_pos(1, d(r),r) + crys%rang_a_pos(1, d(r),r)
              crys%Pv_refi(crys%p(np)) =crys%a_pos(1, d(r),r)
            end if
            if (abs(crys%ref_a_pos(2,d(r),r)) > 0.0) then
              np = np + 1
              crys%list (np) = crys%a_pos(2, d(r),r)
              write(unit=namepar(np),fmt="(a,2i1)")'pos_y',d(r),r
              crys%cod(np) = 1
              ab =  (abs(crys%ref_a_pos(2, d(r),r))/10.0)
              crys%p(np)= int(ab)
              crys%mult(np) = ((abs(crys%ref_a_pos(2, d(r),r)))-(10.0* &
                REAL(crys%p(np))))*SIGN(1.0,(crys%ref_a_pos(2, d(r),r)))
              crys%vlim1(crys%p(np)) = crys%a_pos(2, d(r),r) - crys%rang_a_pos(2, d(r),r)
              crys%vlim2(crys%p(np)) = crys%a_pos(2, d(r),r) + crys%rang_a_pos(2, d(r),r)
              crys%Pv_refi(crys%p(np)) =crys%a_pos(2, d(r),r)
            end if
            if (abs(crys%ref_a_pos(3, d(r),r)) > 0.0) then
              np = np + 1
              crys%list (np) = crys%a_pos(3, d(r),r)
              write(unit=namepar(np),fmt="(a,2i1)")'pos_z',d(r),r
              crys%cod(np) = 1
              ab =  (abs(crys%ref_a_pos(3, d(r),r))/10.0)
              crys%p(np)= int(ab)
              crys%mult(np) = ((abs(crys%ref_a_pos(3, d(r),r)))-(10.0* &
                REAL(crys%p(np))))*SIGN(1.0,(crys%ref_a_pos(3, k,r)))
              crys%vlim1(crys%p(np)) = crys%a_pos(3, d(r),r) - crys%rang_a_pos(3, d(r),r)
              crys%vlim2(crys%p(np)) = crys%a_pos(3, d(r),r) + crys%rang_a_pos(3, d(r),r)
              crys%Pv_refi(crys%p(np)) =crys%a_pos(3, d(r),r)
            end if
            if (abs(crys%ref_a_B( d(r),r)) > 0.0) then
              np = np + 1
              crys%list (np) = crys%a_B( d(r),r)
              write(unit=namepar(np),fmt="(a,2i1)")'Biso',d(r),r
              crys%cod(np) = 1
              ab =  (abs(crys%ref_a_B(d(r),r))/10.0)
              crys%p(np)= int(ab)
              crys%mult(np) = ((abs(crys%ref_a_B( d(r),r)))-(10.0* &
                REAL(crys%p(np))))*SIGN(1.0,(crys%ref_a_B( d(r),r)))
              crys%vlim1(crys%p(np)) = crys%a_B( d(r),r) - crys%rang_a_B( d(r),r)
              if (crys%vlim1(crys%p(np))   .LT. 0 ) crys%vlim1(crys%p(np)) = 0
              crys%vlim2(crys%p(np)) = crys%a_B( d(r),r) + crys%rang_a_B( d(r),r)
              crys%Pv_refi(crys%p(np)) =crys%a_B( d(r),r)
            end if
            if (abs(crys%ref_a_pos(1, d(r),r)) > 0.0) crys%rang_a_pos (1,d(r),r) = 0.0
            if (abs(crys%ref_a_pos(2, d(r),r)) > 0.0) crys%rang_a_pos (2,d(r),r) = 0.0
            if (abs(crys%ref_a_pos(3, d(r),r)) > 0.0) crys%rang_a_pos (3,d(r),r) = 0.0
            if (abs(crys%ref_a_B( d(r),r)) > 0.0) crys%rang_a_B (d(r),r) = 0.0

            crys%l_n_atoms(r) = d(r)
            ok_atom=.true.
            crys%n_actual=r
           !if (d(r) > yy) yy=d(r)

          Case Default

             cycle

          End Select

       end do

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

      integer :: i,i1,i2,k, ier,j, j1, j2, l1, l2, r, l
      integer,dimension(132)                         :: Inte    !  Vector of integer numbers
      integer                                        :: n_int  !number of integers per line
      real(kind=sp),      dimension(132)             :: rel !  Vector of integer numbers
      character(len=132) :: txt
      character(len=25)  :: key, seq
      logical :: ok_explicit, ok_recursive
      real  :: ab

      logi=.true.
      i1=sect_indx(4)
      i2=sect_indx(5)-1

      i=i1

      ok_explicit=.false.; ok_recursive=.false.

      do

        i=i+1; if(i > i2) exit
        txt=adjustl(tfile(i))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        key=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(key)
          Case ("EXPLICIT")
            crys%xplcit = .true.
            crys%l_cnt = 0.0
            i=i+1
            txt=adjustl(tfile(i))
            if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) then
              i=i+1
              txt=adjustl(tfile(i))
            end if
            if (index(txt , 'SPECIFIC')/=0 ) then
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


            else if (index(txt , 'RANDOM')/=0) then
              read (unit = txt, fmt = *, iostat = ier) lstype, crys%l_cnt
                if   (crys%l_cnt==0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR :  Number of layers is missing"
                  logi = .false.
                  return
                end if
                  crys%inf_thick = .false.
                  rndm = .true.

            else if (index(txt , 'SEMIRANDOM')/=0) then
              read (unit = txt, fmt = *, iostat = ier) lstype, crys%l_cnt
                if   (crys%l_cnt==0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR :  Number of layers is missing"
                  logi = .false.
                  return
                end if
                crys%inf_thick = .false.
                rndm = .true.
                i=i+1
                txt=adjustl(tfile(i))
                if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
                if (index(txt , 'SEQ')==1) then
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
                else
                  write(*,*) "ERROR : Layer explicit sequences missing"
                  logi = .false.
                  return
                end if
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
            txt=adjustl(tfile(i))
            if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) then
              i=i+1
              txt=adjustl(tfile(i))
            end if
            if (index(txt , 'INFINITE')==1) then
              crys%inf_thick = .true.
            else
              read(unit=txt ,fmt= * , iostat = ier)crys%l_cnt
              crys%inf_thick = .false.
              i=i+1
              txt=adjustl(tfile(i))
              !Reading refinement codes
              read(unit=txt,fmt=*, iostat=ier) crys%ref_l_cnt
              if(ier /= 0) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading the refinement codes of the number of layers in the crystal"
                logi=.false.
                return
              end if

              !Reading range of parameters
              k=index(txt,"(")
              l=index(txt,")")
              if(k /= 0) then
                read(unit=txt(k+1:l-1),fmt=*, iostat=ier) crys%rang_l_cnt
                if(ier /= 0) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading the ranges of the of the number of layers in the crystal"
                  logi=.false.
                  return
                end if
              else
                Err_crys=.true.
                Err_crys_mess="ERROR: No parameter ranges in the number of layers in the crystal: these must be given!"
                logi=.false.
                return
              end if
              !Creating refinement vectors
              if (abs(crys%ref_l_cnt) > 0.0) then
                np = np + 1        ! to count npar
                crys%list (np) =crys%l_cnt
                write(unit=namepar(np),fmt="(a,2i1)")'num_layers'
                crys%cod(np) = 1
                ab =  (abs(crys%ref_l_cnt)/10.0)
                crys%p(np)= int(ab)
                crys%mult(np) = ((abs(crys%ref_l_cnt))-(10.0* &
                     REAL(crys%p(np))))*SIGN(1.0,(crys%ref_l_cnt))
                crys%vlim1(crys%p(np)) = crys%l_cnt - crys%rang_l_cnt
                crys%vlim2(crys%p(np)) = crys%l_cnt + crys%rang_l_cnt
                crys%Pv_refi(crys%p(np)) = crys%l_cnt
              end if
              if (abs(crys%ref_l_cnt) > 0.0) crys%rang_l_cnt = 0.0
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
      integer, parameter :: eps =1.0E-4
      real  :: ab, suma
      character(len=20), dimension(30) :: citem
      character(len=132) :: txt
      character(len=25)  :: key
      logical :: ok_lt

      logi=.true.
      ok_lt=.false.
      i1=sect_indx(5)
      i2=sect_indx(6)-1

      i=i1+1
      l=1
      j=0
      m=0

        do
          if(i > i2) exit
          txt=adjustl(tfile(i))
          if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) then
           i=i+1
           cycle
          end if
          k=index(txt," ")
          key=txt(1:k-1)
          txt=adjustl(txt(k+1:))

          SELECT CASE (key)

          CASE ("LT")
            j=j+1

            call getword(txt, citem, nitem)
            do m=1, 3
              call read_fraction(citem(1+m), crys%l_r (m,j,l))
              if(ERR_String) then
                write(unit=*,fmt="(a)") trim(ERR_String_Mess)
                logi=.false.
                return
              end if
            end do

            read(unit=citem(1),fmt=*,iostat=ier) crys%l_alpha (j,l)
            if(ier /= 0) then
                   Err_crys=.true.
                   Err_crys_mess="ERROR reading atomic parameters"
                   logi=.false.
                   return
            end if

            i=i+1
            txt=adjustl(tfile(i))

            read (unit = txt, fmt =*, iostat = ier) crys%ref_l_alpha (j,l), crys%ref_l_r (1,j,l), &
                                                    crys%ref_l_r (2,j,l), crys%ref_l_r (3,j,l)

            k=index(txt,"(")
            m=index(txt,")")
            if(k /= 0) then
              read(unit=txt(k+1:m-1),fmt=*, iostat=ier) crys%rang_l_alpha (j,l), crys%rang_l_r(1,j,l), &
                                                        crys%rang_l_r(2,j,l) , crys%rang_l_r(3,j,l)
            end if
            if (abs(crys%ref_l_alpha (j,l)) > 0.0)  then
              np = np + 1    !to count npar
              crys%list (np) = crys%l_alpha (j,l)
              write(unit=namepar(np),fmt="(a,2i1)")'alpha',l,j
              crys%cod(np) = 1
              ab =  (abs(crys%ref_l_alpha (j,l))/10.0)
              crys%p(np)= int(ab)
              crys%mult(np) = ((abs(crys%ref_l_alpha (j,l)))-(10.0* &
                REAL(crys%p(np))))*SIGN(1.0,(crys%ref_l_alpha (j,l)))
              crys%vlim1(crys%p(np))=crys%l_alpha(j,l)- crys%rang_l_alpha (j,l)
              if (crys%vlim1(crys%p(np))  .LT. 0 ) crys%vlim1(crys%p(np)) = 0
              crys%vlim2(crys%p(np)) = crys%l_alpha(j,l)+crys%rang_l_alpha(j,l)
              if (crys%vlim2(crys%p(np))  > 1 ) crys%vlim2(crys%p(np)) = 1
              crys%Pv_refi(crys%p(np)) = crys%l_alpha (j,l)
            end if
            if (abs(crys%ref_l_r(1,j,l)) > 0.0 ) then
              np = np + 1
              crys%list (np) = crys%l_r(1,j,l)
              write(unit=namepar(np),fmt="(a,2i1)")'tx',l,j
              crys%cod(np) = 1
              ab =  (abs(crys%ref_l_r(1,j,l))/10.0)
              crys%p(np)= int(ab)
              crys%mult(np) = ((abs(crys%ref_l_r(1,j,l)))-(10.0* &
                REAL(crys%p(np))))*SIGN(1.0,(crys%ref_l_r(1,j,l)))
              crys%vlim1(crys%p(np)) = crys%l_r (1,j,l) - crys%rang_l_r (1,j,l)
              crys%vlim2(crys%p(np)) = crys%l_r (1,j,l) + crys%rang_l_r (1,j,l)
              crys%Pv_refi(crys%p(np)) = crys%l_r(1,j,l)
            end if
            if (abs(crys%ref_l_r(2,j,l)) > 0.0 ) then
              np = np + 1
              crys%list (np) = crys%l_r(2,j,l)
              write(unit=namepar(np),fmt="(a,2i1)")'ty',l,j
              crys%cod(np) = 1
              ab =  (abs(crys%ref_l_r(2,j,l))/10.0)
              crys%p(np)= int(ab)
              crys%mult(np) = ((abs(crys%ref_l_r(2,j,l)))-(10.0* &
                REAL(crys%p(np))))*SIGN(1.0,(crys%ref_l_r(2,j,l)))
              crys%vlim1(crys%p(np)) = crys%l_r (2,j,l) - crys%rang_l_r (2,j,l)
              crys%vlim2(crys%p(np)) = crys%l_r (2,j,l) + crys%rang_l_r (2,j,l)
              crys%Pv_refi(crys%p(np)) = crys%l_r(2,j,l)
            end if
            if (abs(crys%ref_l_r(3,j,l)) > 0.0 )  then
              np = np + 1
              crys%list (np) = crys%l_r(3,j,l)
              write(unit=namepar(np),fmt="(a,2i1)")'tz',l,j
              crys%cod(np) = 1
              ab =  (abs(crys%ref_l_r(3,j,l))/10.0)
              crys%p(np)= int(ab)
              crys%mult(np) = ((abs(crys%ref_l_r(3,j,l)))-(10.0* &
                REAL(crys%p(np))))*SIGN(1.0,(crys%ref_l_r(3,j,l)))
              crys%vlim1(crys%p(np)) = crys%l_r (3,j,l) - crys%rang_l_r (3,j,l)
              crys%vlim2(crys%p(np)) = crys%l_r (3,j,l) + crys%rang_l_r (3,j,l)
              crys%Pv_refi(crys%p(np)) = crys%l_r(3,j,l)
            end if

            if (abs(crys%ref_l_alpha (j,l)) > 0.0)  crys%rang_l_alpha(j,l) = 0.0
            if (abs(crys%ref_l_r(1,j,l)) > 0.0 )  crys%rang_l_r(1,j,l) = 0.0
            if (abs(crys%ref_l_r(2,j,l)) > 0.0)  crys%rang_l_r(2,j,l) = 0.0
            if (abs(crys%ref_l_r(3,j,l)) > 0.0 )  crys%rang_l_r(3,j,l) = 0.0

            if (j==n_layers) then
                l=l+1
                j=0
            end if
            i=i+1
            ok_lt=.true.


          CASE ("FW")
            read (unit = txt, fmt =*, iostat = ier) crys%r_b11 (j,l) , crys%r_b22 (j,l) , crys%r_b33 (j,l) , &
                                                    crys%r_b12 (j,l) , crys%r_b31 (j,l) , crys%r_b23 (j,l)
            if(ier /= 0) then
              Err_crys=.true.
              Err_crys_mess="ERROR reading fats waller parameters"
              logi=.false.
              return
            end if
            i=i+1

          CASE DEFAULT
            cycle
          end select
        end do

        DO  l = 1, crys%n_typ        !check if l_apha sums one for each layer
          suma = 0
          DO  j = 1, crys%n_typ
            suma = suma + crys%l_alpha(j,l)
          END DO

          IF(ABS(suma - 1.0) > eps) THEN
            write(op,*) "Stacking probabilities from LAYER ",l," do not sum to 1."
            write (*,*)  'the sum is'  , suma
            logi = .false.
          end if
        END DO
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

      integer :: i,i1,i2,k, ier,n, j
      character(len=132) :: txt
      character(len=25)  :: key
      logical :: ok_sim, ok_opt, ok_eps, ok_acc, ok_iout, ok_mxfun, ok_lmq, ok_boxp, ok_maxfun, ok_corrmax, ok_tol, ok_nprint

      logi=.true.
      i1=sect_indx(6)
      if (sect_indx(7) == 0) then
        i2= numberl
      else
        i2=sect_indx(7)-1 !if simulation there is no experimental section
      end if

      i=i1

      ok_sim=.false.; ok_opt=.false.; ok_eps=.false.; ok_acc=.false.; ok_iout=.false.; ok_mxfun=.false.
      ok_lmq=.false.; ok_boxp=.false.; ok_maxfun=.false.; ok_corrmax=.false.; ok_tol=.false. ; ok_nprint=.false.

      do
        i=i+1; if(i > i2) exit
        txt=adjustl(tfile(i))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle

        k=index(txt," ")
        key=txt(1:k-1)

        Select Case(key)

          Case ('SIMULATION')

                  opt = 0
                  txt=adjustl(txt(k+1:))
                  read (unit = txt, fmt = *, iostat=ier) th2_min, th2_max, d_theta
                  if(ier /= 0 ) then
                    Err_crys=.true.
                    Err_crys_mess="ERROR reading 2theta ranges"
                    logi=.false.
                    return
                  end if

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
                  ok_sim = .true.

          Case ("LMQ")
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

          Case ("BOXP")
            txt=adjustl(txt(k+1:))

            read(unit=txt,fmt=*, iostat=ier) Cond%percent
              if(ier /= 0 ) then
                Err_crys=.true.
                Err_crys_mess="ERROR reading boxp parameter"
                logi=.false.
                return
              end if
              if( Cond%percent <= 0.0 .or. Cond%percent >= 300.0)  then
                Cond%constr=.false.
              else
                Cond%constr=.true.
              end if
            ok_boxp=.true.
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
      !Initialise codes to zero
            vs%code=0

            if (Cond%constr .and. ok_corrmax .and. ok_tol .and. ok_maxfun) then
              ok_lmq=.true.
            elseif(ok_corrmax .and. ok_tol .and. ok_maxfun .and. ok_nprint .and. .not. Cond%constr ) then
              ok_lmq=.true.
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

     !    Case ("SIMPLEX")  !TO BE CHECKED /ELIMINATED
     !
     !        i=i+1; if(i > i2) exit
     !        txt=adjustl(tfile(i))
     !        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
     !        k=index(txt," ")
     !        key=txt(1:k-1)
     !        txt=adjustl(txt(k+1:))
     !
     !
     !              read(unit = txt, fmt = *, iostat=ier) opti%mxfun
     !              do
     !                z = z+1; if(z > numberl) exit
     !                txt = adjustl (tfile (z))
     !                if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or.index(txt,'{')==1) cycle
     !                if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
     !                  z=z-1
     !                  exit
     !                else
     !                  read(unit = txt, fmt = *, iostat=ier) opti%eps
     !                  do
     !                    z = z+1; if(z > numberl) exit
     !                    txt = adjustl (tfile (z))
     !                    if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) cycle
     !                    if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
     !                      z=z-1
     !                      exit
     !                    else
     !                      read(unit = txt, fmt = *, iostat=ier) opti%iout
     !                      do
     !                        z = z+1; if(z > numberl) exit
     !                        txt = adjustl (tfile (z))
     !                        if(index(txt,"!") == 1 .or. len_trim(txt) == 0 .or. index(txt,'{')==1) cycle
     !                        if ( index (txt, 'EXPERIMENTAL' ) == 1 ) then
     !                          z=z-1
     !                          exit
     !                        else
     !                          read(unit = txt, fmt = *, iostat=ier) opti%acc
     !                        end if
     !                      end do
     !                    end if
     !                  end do
     !                end if
     !              end do
     !
     !
!////!////////////////CONVERSION TO SIMPLEX VARIABLES///////////////////////////////////////////
     !
     !            !vector(1:nrp)  = crys%list (1:nrp)
     !            !   code(1:nrp) = crys%cod(1:nrp)
     !                    numpar = np
     !            !      nm_cycl  = crys%n_cycles
!////!/////////////////////////////////////////////////////////////////////////////////////////
     !            n_plex = maxval(crys%p)
     !            opti%loops = 2 * n_plex
     !            opti%iquad=1
     !
     !            opti%npar = n_plex
     !
     !            do i = 1, numpar
     !               label(crys%p(i)) = 0
     !            end do
     !
     !            do i=1, numpar                    !construction of SIMPLEX 'in' vectors
     !
     !              if (label (crys%p(i))  == 0 ) then        !not to overwrite in config
     !               !  v_plex(crys%p(i)) = crys%Pv_refi(i)
     !
     !               !  nampar(pnum(i)) = namepar(i)
     !                 sanvec%low(pnum(i)) = crys%vlim1(i)
     !                 sanvec%high(pnum(i)) = crys%vlim2(i)
     !
     !              end if
     !              label (crys%p(i)) = 1
     !            end do
     !

          Case Default
            cycle
          End Select
      end do

     if(ok_sim) then
       return
     else if (ok_opt .and. ok_acc .and. ok_mxfun .and. ok_eps .and. ok_iout) then
       return
     else if (ok_lmq) then
       return
     else
        Err_crys=.true.
        Err_crys_mess="ERROR:  calculation parameters missing!"
        logi=.false.
     end if

    End Subroutine Read_CALCULATION

    Subroutine Read_EXPERIMENTAL(logi)
    logical, intent(out) :: logi

      integer :: i,i1,i2,k,j, ier
      character(len=132) :: txt
      character(len=25)  :: key
      logical :: ok_file, ok_fformat, ok_bgr, ok_bcalc

      logi=.true.
      i1=sect_indx(7)
      i2=numberl

      i=i1

      ok_file=.false.; ok_fformat=.false.; ok_bgr=.false. ; ok_bcalc=.false.

      do

        i=i+1; if(i > i2) exit
        txt=adjustl(tfile(i))
        if( len_trim(txt) == 0 .or. txt(1:1) == "!" .or. txt(1:1) == "{" ) cycle
        k=index(txt," ")
        key=txt(1:k-1)
        txt=adjustl(txt(k+1:))

        Select Case(key)


        Case("FILE")

            read(unit=txt,fmt=*, iostat=ier)   dfile !, th2_min, th2_max, d_theta
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading file instruction"
                  logi=.false.
                  return
              end if
              if (th2_min /= 0 .and.  th2_max /= 0  .and. d_theta /= 0) then
                th2_min = th2_min * deg2rad
                th2_max = th2_max * deg2rad
                d_theta = half * deg2rad * d_theta
              end if
              ok_file=.true.

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

          Case("BGR")
            read(unit=txt,fmt=*,iostat=ier)  background_file
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading background instruction"
                  logi=.false.
                  return
              end if
              ok_bgr=.true.

          Case("BCALC")
            read (unit = txt, fmt = *, iostat=ier) mode
              if(ier /= 0 ) then
                  Err_crys=.true.
                  Err_crys_mess="ERROR reading background calculation instruction"
                  logi=.false.
                  return
              end if
              ok_bcalc=.true.


          Case Default

             cycle

        End Select
      end do

      if(ok_file .and. ok_fformat .and. ok_bgr .and. ok_bcalc) then
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
        allocate(crys%vlim1(max_npar))
        crys%vlim1=0

        if(allocated (crys%vlim2)) deallocate(crys%vlim2)
        allocate(crys%vlim2(max_npar))
        crys%vlim2=0

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


!_____________________________________________________________________________________________________
    Subroutine new_getfil(stfile,vecsan, olg)


       character(len=31) ,             intent (in out)      :: stfile
       type (State_Vector_Type),       intent (   out)      :: vecsan
       logical                   ,     intent (   out)      :: olg

        write(unit=op,fmt="(a)",advance="no") ' => Enter the complete name of the structure input file: '
        read(unit= *,fmt="(a)") stfile

        !WRITE(op,fmt=*) "=> Looking for scattering factor data file '",  sfname(:),"'"
        OPEN(UNIT = sf, FILE = sfname)
        !WRITE(op,fmt=*) "=> Opening scattering factor data file '",  sfname(:),"'"

        call read_structure_file(stfile, vecsan, olg)

        if (err_crys) then
          print*, trim(err_crys_mess)
        else
          write(op, fmt=*) "=> Structure input file read in"
        end if
        return
    End subroutine

!______________________________________________________________________________________________________

 !    INTEGER*4 FUNCTION choice(flag, list, n)          !from diffax
 !
 !
 !    IMPLICIT NONE
 !    CHARACTER (LEN=*), INTENT(IN)        :: flag
 !    CHARACTER (LEN=80), INTENT(IN)       :: list(n)
 !    INTEGER*4, INTENT(IN)                :: n
 !
 !    INTEGER*4 i, j1, j2
 !
 !    i = 1
 !    j1 = length(flag)
 !
 !    10 j2 = length(list(i))
!!see if the string contained in list(i) is identical to that in flag
 !    IF(j1 == j2 .AND. INDEX(flag, list(i)(1:j2)) == 1) GO TO 20
 !    i = i + 1
 !    IF(i <= n) GO TO 10
 !
 !    20 IF(i > n) i = -1
 !    choice = i
 !
 !    RETURN
 !    END FUNCTION choice
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

        Subroutine read_structure_file (namef, sanvec, logi)

        character(len=*)    ,          intent(in    )     :: namef
        type (State_Vector_Type),      intent(   out)     :: sanvec
        logical                   ,    intent(   out)     :: logi



        logical                                        :: esta
        integer                                        :: ier, a,b,l,j, i
        !real,               dimension(80)              :: label

        logi = .true.

        call Set_Crys()

        b = 0
        a = 0
        l = 0
        j = 0

        inquire(file=namef,exist=esta)
        if( .not. esta) then
          Err_crys=.true.
          Err_crys_mess="ERROR :  The file "//trim(namef)//" doesn't exist"
          logi=.false.
          return
        else
          open(unit=i_data,file=trim(namef),status="old",action="read",position="rewind",iostat=ier)
          if(ier /= 0) then
            Err_crys=.true.
            Err_crys_mess="ERROR opening the file "//trim(namef)
            logi=.false.
            return
          end if
        end if

        call Set_TFile(namef,logi)
        if(.not. logi) return

        call Read_INSTRUMENTAL(logi)


        if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
        else
          write(*,"(a,i2)")     " Radiation:            ", crys%rad_type
          write(*,"(a,3f10.4)") " Lambda:               ", crys%lambda, crys%lambda2, crys%ratio
          write(*,"(a,6f10.2)") " peak-width parameters:", crys%p_u, crys%p_v, crys%p_w, crys%p_x, crys%p_dg, crys%p_dl
          write(*,"(a,6f10.2)") " peak-width codes:     ", crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, &
                                                           crys%ref_p_x,  crys%ref_p_dg,  crys%ref_p_dl
          write(*,"(a,6f10.2)") " peak-width ranges:    ", crys%rang_p_u,crys%rang_p_v, crys%rang_p_w, &
                                                           crys%rang_p_x, crys%rang_p_dg, crys%rang_p_dl
        end if
        call Read_STRUCTURAL(logi)
        if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
        else
          write(*,"(a,4f10.4)") " Cell parameters:       ", crys%cell_a, crys%cell_b, crys%cell_c, crys%cell_gamma
          write(*,"(a,4f10.2)") " cell parameter codes:  ", crys%ref_cell_a, crys%ref_cell_b, crys%ref_cell_c, &
                                                            crys%ref_cell_gamma
          write(*,"(a,4f10.2)") " cell parameter ranges: ",    crys%rang_cell_a, crys%rang_cell_b, crys%rang_cell_c, &
                                                            crys%rang_cell_gamma
          write(*,*)            " Laue symmetry:         ", crys%sym
          write(*,"(a,2f10.2)")     " layer width parameters:", crys%layer_a, crys%layer_b
          write(*,"(a,6f10.2)") " layer width codes:     ", crys%ref_layer_a, crys%ref_layer_b
          write(*,"(a,6f10.2)") " layer width ranges:    ", crys%rang_layer_a, crys%rang_layer_b
        end if

        call Read_LAYER(logi)

        if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
        else
            b=1
            a=1
          do b=1, crys%n_actual
              write(*,*)            " Layer symmetry:         ", crys%centro(b)
            do a=1, crys%l_n_atoms(b)
              write(*,*) "a,b,crys%a_num(a,b)"  , a,b,crys%a_num(a,b)
              write(*,"(2a,i4, 5f10.2)") " Atomic parameters:       ", crys%a_name(a,b), crys%a_num(a,b), crys%a_pos(1, a,b), &
                                           crys%a_pos(2, a,b), crys%a_pos(3, a,b), crys%a_B (a,b), crys%a_occup(a,b)
              write(*,"(a,4f10.2)")      " Atomic parameter codes:  ", crys%ref_a_pos(1, a,b), crys%ref_a_pos(2, a,b), &
                                           crys%ref_a_pos(3, a,b), crys%ref_a_B(a,b)
              write(*,"(a,4f10.2)")      " atomic parameter ranges: ", crys%rang_a_pos(1, a,b), crys%rang_a_pos(2, a,b), &
                                           crys%rang_a_pos(3, a,b), crys%rang_a_B(a,b)
            end do
          end do
           write(*,*) "n_layers", n_layers
        end if

        call Read_STACKING (logi)
        if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
        else
          if (crys%xplcit) then
            write(*, "(2a)") "Stacking type: EXPLICIT, ", lstype
            write(*, "(a, f5.2)") "Number of layers:", crys%l_cnt
            !a = 1
            !do a=1, int(crys%l_cnt)
                if (crys%l_seq(a) /=0) then
                  write(*,*) "Sequence:", crys%l_seq(1:crys%l_cnt)
                end if
            !end do
          else
             write(*, "(2a)") "Stacking type: RECURSIVE"
             if (crys%inf_thick) then
               write (*, "(a)") "Number of layers: INFINITE"
             else
               write (*, "(a, f5.2)") "Number of layers:", crys%l_cnt
             end if
           end if
        end if

        call Read_TRANSITIONS (logi)

        if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
        else
          l=1
          j=1
          do l=1, n_layers
            do j=1, n_layers
              write(*, "(a,i2, a, i2,a)") "Layer transition parameters between layers ", l, " and ", j, ":"
              write(*, "(a,4f10.2)")  "Stacking probability and stacking vector:", crys%l_alpha (j,l), crys%l_r (1,j,l), &
                                      crys%l_r (2,j,l), crys%l_r (3,j,l)
              write(*, "(a, 4f10.2)") "Refinement codes:", crys%ref_l_alpha (j,l), crys%ref_l_r (1,j,l),crys%ref_l_r (2,j,l), &
                                      crys%ref_l_r (3,j,l)
              write(*, "(a, 4f10.2)") "Refinement ranges:",crys%rang_l_alpha (j,l), crys%rang_l_r(1,j,l), crys%rang_l_r(2,j,l) ,&
                                      crys%rang_l_r(3,j,l)
              write(*, "(a, 6f10.2)") "Fats Waller parameters:",crys%r_b11 (j,l) , crys%r_b22 (j,l) , crys%r_b33 (j,l) , &
                                      crys%r_b12 (j,l) , crys%r_b31 (j,l) , crys%r_b23 (j,l)
            end do
          end do
        end if


        call Read_CALCULATION (logi)

        crys%npar=np

        if(.not. logi) then
          write(*,"(a)")  " => "//Err_crys_mess
          return
        else

          if (opt==3) then        !construction of some optimization variables(DFP)
            numpar = crys%npar
            n_plex = maxval(crys%p)
            opti%loops = 2 * n_plex
            opti%iquad=1
            opti%npar = n_plex
            do i = 1, numpar
              sanvec%low(crys%p(i)) = crys%vlim1(i)
              sanvec%high(crys%p(i)) = crys%vlim2(i)
            end do
          end if
          if (opt==4) then         !construction of some optimization variables  (LMQ)

            opti%npar =maxval(crys%p)
            Cond%npvar=opti%npar
            numpar = crys%npar
            vs%pv(1:opti%npar)= crys%Pv_refi(1:opti%npar)
            vs%code(1:opti%npar) = 1
            vs%np= opti%npar
          end if


          write(*,"(a,i2)") " Type of calculation (0=simulation, 3=local_optimizer, 4=LMQ): ",  opt
          if (opt == 3) then
            write(*,"(2a)") " Method of calculation: ", opti%method
            write(*,"(a,i4 )") " Maximum number of function evaluations: ", opti%mxfun
            write(*,"(a,f5.2 )") "Stopping criterion: ", opti%eps
            write(*,"(a,i2 )") "Output file indicator: ", opti%iout
            write(*,"(a,f12.9 )") "Accuracy: ", opti%acc
          elseif (opt == 4) then
            write(*,"(a,f5.2, 2i4, f12.9, i4)") "percent, corrmax, maxfun, tol, nprint: ", Cond%percent, Cond%corrmax, &
                                             Cond%icyc, cond%tol, cond%nprint
          else
            write(*,"(a,3f7.2 )") "2Theta max, 2Theta min, stepsize:",  th2_min, th2_max, d_theta
            th2_min = th2_min * deg2rad
            th2_max = th2_max * deg2rad
            d_theta = half * deg2rad * d_theta
          end if
        end if

        if (opt /= 0) then !not necessary for simulation

          call Read_EXPERIMENTAL(logi)

          if(.not. logi) then
            write(*,"(a)")  " => "//Err_crys_mess
            return
          else
            if (th2_min == 0 .and.  th2_max == 0 .and. d_theta == 0) then
              write(*, "(2a)") "File name:", dfile
            else
              write(*,"(2a,3f8.5)") " File parameters:", dfile , th2_min, th2_max, d_theta
            end if
            write(*, "(2a)") "File format: ", fmode
            write(*, "(2a)") "Background file name: ", background_file
            write(*, "(2a)") "Background calculation type: ", mode
          end if
        end if

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
           a_name(1:crys%yy,1:crys%n_typ)   = crys%a_name(1:crys%yy,1:crys%n_typ)
           a_number(1:crys%yy,1:crys%n_typ) = crys%a_num(1:crys%yy,1:crys%n_typ)
          a_pos(1:3,1:crys%yy,1:crys%n_typ) = crys%a_pos(1:3,1:crys%yy,1:crys%n_typ)
             a_b(1:crys%yy,1:crys%n_typ)    = crys%a_B(1:crys%yy,1:crys%n_typ)
            a_occup(1:crys%yy,1:crys%n_typ) = crys%a_occup(1:crys%yy,1:crys%n_typ)
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
                              pv_gamma = crys%p_gamma
                           trim_origin = crys%trm
                            n_actual   = crys%n_actual
                              l_cnt    = crys%l_cnt
                       l_seq(1:xp_max) = crys%l_seq(1:xp_max)
                l_actual(1:crys%n_typ) = crys%l_actual(1:crys%n_typ)
                           SymGrpNo    = crys%SymGrpNo
              r_b11(1:crys%yy,1:crys%n_typ) = crys%r_b11 (1:crys%yy,1:crys%n_typ)
              r_b22(1:crys%yy,1:crys%n_typ) = crys%r_b22 (1:crys%yy,1:crys%n_typ)
              r_b33(1:crys%yy,1:crys%n_typ) = crys%r_b33 (1:crys%yy,1:crys%n_typ)
              r_b12(1:crys%yy,1:crys%n_typ) = crys%r_b12 (1:crys%yy,1:crys%n_typ)
              r_b31(1:crys%yy,1:crys%n_typ) = crys%r_b31 (1:crys%yy,1:crys%n_typ)
              r_b23(1:crys%yy,1:crys%n_typ) = crys%r_b23 (1:crys%yy,1:crys%n_typ)
              l_n_atoms(1:crys%n_typ)  = crys%l_n_atoms(1:crys%n_typ)
               low_atom(1:crys%n_typ)  = crys%low_atom(1:crys%n_typ)
              high_atom(1:crys%n_typ)  = crys%high_atom(1:crys%n_typ)
                             tolerance = crys%tolerance
                        pnum(1:numpar) = crys%p(1:numpar)
                        mult(1:numpar) = crys%mult(1:numpar)
                              n_layers = crys%n_typ

        return
        End subroutine  read_structure_file

    End module read_data
