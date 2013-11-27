
  Module Dif_compl

      use diffax_mod
      use CFML_GlobalDeps,            only : sp
      use CFML_Diffraction_Patterns , only : diffraction_pattern_type
      use CFML_Optimization_General,  only : Opt_Conditions_Type
      use read_data,                  only : opti, crys_2d_type, crys, cond
      use CFML_LSQ_TypeDef,           only : LSQ_State_Vector_Type

      implicit none

      private

      !public subroutines
      public :: scale_factor_lmq, Write_Prf, Write_ftls, Faults2diffax, vs2faults, Var_assign

      contains
!________________________________________________________________________________________________________________________

       Subroutine Var_assign(state)

    real, dimension(:), intent(in) :: state
    !local variables
    integer :: i, j, k, a, b

          !////////CALCULATED PATTERN VARIABLES ASSIGNMENT////////////////////////////////////

      do i=1, crys%npar

        if (namepar(i) ==  'u')          pv_u   = state(i)
        if (namepar(i) ==  'v')          pv_v   = state(i)
        if (namepar(i) ==  'w')          pv_w   = state(i)
        if (namepar(i) ==  'x')          pv_x   = state(i)
        if (namepar(i) ==  'Dg')         pv_dg  = state(i)
        if (namepar(i) ==  'Dl')         pv_dl  = state(i)
        if (namepar(i) ==  'cell_a')     cell_a = state(i)
        if (namepar(i) ==  'cell_b')     cell_b = state(i)
        if (namepar(i) ==  'cell_c')     cell_c = state(i)
        if (namepar(i) ==  'num_layers') l_cnt  = state(i)
        if (namepar(i) ==  "diameter_a") Wa = state(i)
        if (namepar(i) ==  "diameter_b") Wa = state(i)
        if (namepar(i) ==  'zero_shift') crys%zero_shift  = state(i)
        if (namepar(i) ==  'sycos')      crys%sycos  = state(i)
        if (namepar(i) ==  'sysin')      crys%sysin  = state(i)
        if (index( namepar(i) ,'cell_gamma') == 1)   cell_gamma = state(i)

        do j=1, n_layers
          do k=1, n_atoms
            if (index (namepar(i) , 'pos_x' )== 1)     then
                read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
                a_pos(1,a,b)  = state(i) * pi2                         !need to invert conversion done by routine nmcoor (diffax_calc)
            end if
            if (index (namepar(i) ,'pos_y' )== 1)    then
                read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
                a_pos(2,a,b)  = state(i) * pi2
            end if
            if (index (namepar(i), 'pos_z' ) == 1 )   then
                read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
                a_pos(3,a,b)  = state(i) * pi2
            end if
            if (index (namepar(i),'Biso')==1) then
                read (unit = namepar(i)(5:6), fmt = "(2i1)" ) a,b
                a_b(a,b)  = state(i)

            end if
            if (index( namepar(i) ,  'alpha' ) == 1)    then
                read (unit = namepar(i)(6:7), fmt = "(2i1)" ) b,a
                l_alpha(a,b)  = state(i)

            end if
            if (index (namepar(i), 'tx' )== 1 )    then
                read (unit = namepar(i)(3:4), fmt = "(2i1)" ) b,a
                l_r(1,a,b)  = state(i)
            end if
            if (index (namepar(i), 'ty' )== 1 )    then
                read (unit = namepar(i)(3:4), fmt = "(2i1)" ) b,a
                l_r(2,a,b)  = state(i)
            end if
            if (index (namepar(i), 'tz' ) == 1)     then
                read (unit = namepar(i)(3:4), fmt = "(2i1)" ) b,a
                l_r(3,a,b)  = state(i)
            end if
          end do
        end do
      end do

    End subroutine Var_assign

    Subroutine vs2faults(vs, crys)

       type(crys_2d_type),     intent(in out)   :: crys
       type(LSQ_State_Vector_type), Intent(In)  :: Vs

       !local variables
       integer        :: n, j, i
       real, dimension(max_npar):: shift

       do i=1, opti%npar
         shift(1:opti%npar) = vs%pv(1:opti%npar) - crys%Pv_refi(1:opti%npar)
       end do

       do i = 1, crys%npar
             crys%list(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
       end do

       call Var_assign(crys%list)

 !     crys%n_typ     = n
 !     crys%yy        = j
 !     crys%rad_type  = rad_type
 !     crys%lambda    = lambda
 !     crys%lambda2   = lambda2
 !     crys%ratio     = ratio
 !     crys%broad     = blurring
 !     crys%cell_a    = cell_a
 !     crys%cell_b    = cell_b
 !     crys%cell_c    = cell_c
 !     crys%cell_gamma= cell_gamma
 !     crys%sym       = pnt_grp
 !     crys%centro    = l_symmetry
 !     crys%a_name(1:j,1:n)     = a_name(1:j,1:n)
 !     crys%a_num(1:j,1:n)      = a_number(1:j,1:n)
 !     crys%a_pos(1:3,1:j,1:n)  = a_pos(1:3,1:j,1:n)
 !     crys%a_B(1:j,1:n)        = a_b(1:j,1:n)
 !     crys%a_occup(1:j,1:n)    = a_occup(1:j,1:n)
 !     crys%recrsv              = recrsv
 !     crys%xplcit              = xplcit
 !     crys%finite_width        = finite_width
 !     crys%inf_thick           = inf_thick
 !     crys%l_alpha(1:n,1:n)    = l_alpha(1:n,1:n)
 !     crys%l_r(1:3,1:n,1:n)    = l_r(1:3,1:n,1:n)
 !     crys%layer_a             = Wa
 !     crys%layer_b             = Wb
 !     crys%p_u                 = pv_u
 !     crys%p_v                 = pv_v
 !     crys%p_w                 = pv_w
 !     crys%p_x                 = pv_x
 !     crys%p_dg                = pv_dg
 !     crys%p_dl                = pv_dl
 !     crys%p_gamma             = pv_gamma
 !     crys%trm                 = trim_origin
 !     crys%n_actual            = n_actual
 !     crys%l_cnt               = l_cnt
 !     crys%l_seq(1:xp_max)     = l_seq(1:xp_max)
 !     crys%l_actual(1:n)       = l_actual(1:n)
 !     crys%SymGrpNo            = SymGrpNo
 !     crys%r_b11 (1:j,1:n)     = r_b11 (1:j,1:n)
 !     crys%r_b22 (1:j,1:n)     = r_b22 (1:j,1:n)
 !     crys%r_b33 (1:j,1:n)     = r_b33 (1:j,1:n)
 !     crys%r_b12 (1:j,1:n)     = r_b12 (1:j,1:n)
 !     crys%r_b31 (1:j,1:n)     = r_b31 (1:j,1:n)
 !     crys%r_b23 (1:j,1:n)     = r_b23 (1:j,1:n)
 !     crys%l_n_atoms(1:n)      = l_n_atoms(1:n)
 !      crys%low_atom(1:n)      = low_atom(1:n)
 !     crys%high_atom(1:n)      = high_atom(1:n)
 !     crys%tolerance           = tolerance
 !     crys%mult(1:crys%npar)   = mult(1:crys%npar)
 !     crys%n_typ               = n_layers
 !     crys%ttl                 = ttl
 !     crys%fundamental         = fundamental
 !     crys%original            = original
 !     crys%randm               = randm
 !     crys%semirandm           = semirandm
 !     crys%spcfc               = spcfc
 !     crys%fls                 = fls
 !     crys%lls                 = lls
 !     crys%otls                = otls
 !     crys%stls                = stls


       End subroutine vs2faults

       Subroutine Faults2diffax(crys)

       Type(crys_2d_type),     intent(in) :: crys
       integer        :: n, j

                   n  = crys%n_typ
                   j  = crys%yy
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
          a_name(1:j,1:n)   = crys%a_name(1:j,1:n)
          a_number(1:j,1:n) = crys%a_num(1:j,1:n)
         a_pos(1:3,1:j,1:n) = crys%a_pos(1:3,1:j,1:n)
            a_b(1:j,1:n)    = crys%a_B(1:j,1:n)
           a_occup(1:j,1:n) = crys%a_occup(1:j,1:n)
                   recrsv   = crys%recrsv
                   xplcit   = crys%xplcit
               finite_width = crys%finite_width
                 inf_thick  = crys%inf_thick
         l_alpha(1:n,1:n)   = crys%l_alpha(1:n,1:n)
         l_r(1:3,1:n,1:n)   = crys%l_r(1:3,1:n,1:n)
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
              l_actual(1:n) = crys%l_actual(1:n)
                 SymGrpNo   = crys%SymGrpNo
             r_b11(1:j,1:n) = crys%r_b11 (1:j,1:n)
             r_b22(1:j,1:n) = crys%r_b22 (1:j,1:n)
             r_b33(1:j,1:n) = crys%r_b33 (1:j,1:n)
             r_b12(1:j,1:n) = crys%r_b12 (1:j,1:n)
             r_b31(1:j,1:n) = crys%r_b31 (1:j,1:n)
             r_b23(1:j,1:n) = crys%r_b23 (1:j,1:n)
             l_n_atoms(1:n) = crys%l_n_atoms(1:n)
              low_atom(1:n) =  crys%low_atom(1:n)
             high_atom(1:n) = crys%high_atom(1:n)
             tolerance      = crys%tolerance
          mult(1:crys%npar) = crys%mult(1:crys%npar)
                   n_layers = crys%n_typ
                        ttl = crys%ttl
                fundamental = crys%fundamental
                   original = crys%original
                      randm = crys%randm
                  semirandm = crys%semirandm
                      spcfc = crys%spcfc
                       fls  = crys%fls
                       lls  = crys%lls
                       otls = crys%otls
                       stls = crys%stls

       End subroutine Faults2diffax

       Subroutine Write_ftls(crys,i_ftls, vs)
          !-----------------------------------------------
          !   D u m m y   A r g u m e n t s
          !-----------------------------------------------
          Type(crys_2d_type),     intent(in  out) :: crys
          integer,                intent(in)      :: i_ftls
          type(LSQ_State_Vector_type), Intent(In) :: Vs


          integer                            :: a,b,c, j, l,i ,n, print_width
          CHARACTER(LEN=80)                  :: list(2)

          call vs2faults(vs, crys)

          write(i_ftls,"(a)")          "TITLE"
          write(i_ftls,"(a)")          crys%ttl

          write(i_ftls,"(a)")              "  "
          write(i_ftls,"(a)")          "INSTRUMENTAL  AND  SIZE  BROADENING"
          write(i_ftls,"(a)")          "!type of radiation"
          if (crys%rad_type == 0 ) then
            write(i_ftls,"(a)")     " Radiation      X-RAY"
          elseif (crys%rad_type == 1 ) then
            write(i_ftls,"(a)")     " Radiation      NEUTRON"
          else
            write(i_ftls,"(a)")     " Radiation      ELECTRON"
          end if
          write(i_ftls,"(a)")          "!               lambda1    lambda2   ratio "
          write(i_ftls,"(a,3f10.4)") " Wavelength  ", crys%lambda , crys%lambda2 , crys%ratio

          write(i_ftls,"(a)")          "!instrumental aberrations      zero      sycos     sysin"
          write(i_ftls,"(a,3f10.4)") " Aberrations               ", crys%zero_shift, crys%sycos, crys%sysin
          write(i_ftls,"(tr25,3f10.2)") crys%ref_zero_shift, crys%ref_sycos,  crys%ref_sysin

          if (crys%broad == ps_vgt .and. crys%trm) then
            write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dg         Dl"
            write(i_ftls,"(a,4f12.6,2f11.2, a)") " PSEUDO-VOIGT   ", pv_u, pv_v, pv_w, pv_x, pv_dg,pv_dl, " TRIM"
            write(i_ftls,"(tr16,4f12.2,2f11.2)")  crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, crys%ref_p_x,  crys%ref_p_dg, &
                                           crys%ref_p_dl
          elseif (crys%broad == ps_vgt .and. .not. crys%trm ) then
            write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dg         Dl"
            write(i_ftls,"(a,4f12.6,2f11.2)")   " PSEUDO-VOIGT   ", pv_u, pv_v, pv_w, pv_x, pv_dg, pv_dl
            write(i_ftls,"(tr16,4f12.2,2f11.2)")  crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, crys%ref_p_x,  crys%ref_p_dg, &
                                                  crys%ref_p_dl
          elseif (crys%broad == pv_gss .and. crys%trm) then
            write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dg"
            write(i_ftls,"(a,4f12.6,f11.2, a)") "    GAUSSIAN    ", pv_u, pv_v, pv_w, pv_x, pv_dg, "TRIM"
            write(i_ftls,"(tr16,4f12.2,f11.2)")  crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, crys%ref_p_x,  crys%ref_p_dg
          elseif (crys%broad == pv_gss .and. .not. crys%trm ) then
            write(i_ftls,"(a)")  "!instr. broadening       u           v           w           x         Dg"
            write(i_ftls,"(a,4f12.6,1f11.2)") "    GAUSSIAN    ", pv_u, pv_v, pv_w, pv_x, pv_dg
            write(i_ftls,"(tr16,4f12.2,f11.2)")  crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, crys%ref_p_x,  crys%ref_p_dg
          elseif (crys%broad == pv_lrn .and. crys%trm ) then
            write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dl"
            write(i_ftls,"(a,4f12.6,1f11.2 a)") "   LORENTZIAN   ", pv_u, pv_v, pv_w, pv_x, pv_dl, "TRIM"
            write(i_ftls,"(tr16,4f12.2,f11.2)")  crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, crys%ref_p_x,  &
                                                 crys%ref_p_dl
          elseif   (crys%broad==pv_lrn .and. .not. crys%trm) then
            write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dl"
            write(i_ftls,"(a,4f12.6,1f11.2)") "   LORENTZIAN   ", pv_u, pv_v, pv_w, pv_x, pv_dl
            write(i_ftls,"(tr16,4f12.2,f11.2)")  crys%ref_p_u, crys%ref_p_v,  crys%ref_p_w, crys%ref_p_x,  &
                                                 crys%ref_p_dl
          else
            write(*,*) "ERROR writing *.ftls file: Problem with instrumental parameters!"
            return
          end if


          write(i_ftls,"(a)")              "  "
          write(i_ftls,"(a)")          " STRUCTURAL  "
          write(i_ftls,"(a)")          " !         a            b            c          gamma "
          write(i_ftls,"(a,3f12.6,f12.2)")   " CELL  ", cell_a, cell_b, cell_c, cell_gamma*rad2deg
          write(i_ftls,"(tr7,4f12.2)")  crys%ref_cell_a, crys%ref_cell_b, crys%ref_cell_c, crys%ref_cell_gamma
          write(i_ftls,"(a)")        "!Laue symmetry"
          write(i_ftls,*)            "SYMM  ", pnt_grp
          write(i_ftls,"(a)")        "!number of layer types"
          write(i_ftls,"(a,i2)")            " NLAYERS  ", n_layers
          write(i_ftls,"(a)")       "!layer width"
          if (crys%finite_width) then
            write(i_ftls,"(a,2f10.2)") " LWIDTH",   Wa, Wb
            write(i_ftls,"(2f10.2)")    crys%ref_layer_a, crys%ref_layer_b
          else
            write(i_ftls,"(a)")        " LWIDTH  INFINITE"
          end if

          b=1
          a=1
          c=0
          do b=1, n_layers
            write(i_ftls,"(a)")              "  "
            if (fundamental(b) .eqv. .false.) then
              write(i_ftls,"(a,i2,a,i2)") "LAYER",b," =", original(b)
            else
              c=c+1
              write(i_ftls,"(a, i2)")  " LAYER", b
              list(1) = 'NONE '
              list(2) = 'CENTROSYMMETRIC '
             !WRITE(dmp,100) 'symmetry = ', list(l_symmetry(i)+1)
              write(i_ftls,"(a)")        "!Layer symmetry   "
              write(i_ftls,"(2a)")      " LSYM   ", list(l_symmetry(c)+1)
              do a=1, crys%l_n_atoms(c)
                write(i_ftls,"(2a)")      "!Atom name   number   x         y         z         Biso      Occ  "
                write(i_ftls,"(2a,i9, 5f10.5)") " ATOM ", a_name(a,c), a_number(a,c), a_pos(1, a,c)/pi2, &
                                           a_pos(2, a,c)/pi2,a_pos(3, a,c)/pi2, a_B (a,c), a_occup(a,c)
                write(i_ftls,"(tr19,5f10.2)") crys%ref_a_pos(1, a,c), crys%ref_a_pos(2, a,c), &
                                                  crys%ref_a_pos(3, a,c), crys%ref_a_B(a,c), crys%ref_a_occup(a,c)
              end do
            end if
          end do
          write(i_ftls,"(a)")              "  "
          write(i_ftls,"(a)")          " STACKING"
          write(i_ftls,"(a)")     "!stacking type"
          if (crys%xplcit) then
            write(i_ftls, "(a)") " EXPLICIT "
            !if (rndm) then
            !   write(i_ftls, " (a,f5.2)") "RANDOM ",l_cnt
           ! else
            !   write(i_ftls,"(a)") lstype
               if (semirandm) then
                 write(i_ftls,"(a)", advance="no") " SEMIRANDOM "
                 write(i_ftls,"(a)") "!number of layers"
                 write(i_ftls,"(f6.1)") l_cnt
                 write(i_ftls,"(a,i5,i5,i5,i5)") " SEQ", fls, lls, otls, stls
                ! i=1
                ! do i=1,crys%n_seq
                   !write(i_ftls,"(a)") "SEQ"             !----------------TO BE FINISHED
                ! end do
               elseif(spcfc) then
                 write(i_ftls,"(a)") " SPECIFIC"
                 write(i_ftls, "(25i2)") crys%l_seq(1:crys%l_cnt)
               elseif(randm) then
                 write(i_ftls,"(a)", advance="no") " RANDOM "
                 write(i_ftls,"(a)") "!number of layers"
                 write(i_ftls,"(f6.1)") l_cnt
               else
                 write(i_ftls,"(a)") "ez du ondo egin"
               end if
           ! end if
            !a = 1
            !do a=1, int(crys%l_cnt)
            !    if (crys%l_seq(a) /=0) then
            !      write(i_ftls,*)  crys%l_seq(1:crys%l_cnt)
            !    end if
            !end do                                      !_______________________________
          else
             write(i_ftls, "(a)") " RECURSIVE"
             write(i_ftls,"(a)") "!number of layers"
             if (crys%inf_thick) then
               write (i_ftls, "(a)") " INFINITE"
             else
               write (i_ftls, "( f6.1)") l_cnt
               write (i_ftls, "( f6.1,a)")  crys%ref_l_cnt
             end if
           end if
          write(i_ftls,"(a)")              "  "
          write(i_ftls,"(a)")          " TRANSITIONS"
          l=1
          j=1
          do l=1, n_layers
            write(i_ftls,"(a)")              "  "
            write(i_ftls,"(a)")              "  "
            do j=1, n_layers
              write(i_ftls, "(a,i2, a, i2)") "!layer ", l, " to layer ", j
              write(i_ftls, "(a, 4f10.6)")  "LT ",  l_alpha (j,l), l_r (1,j,l), l_r (2,j,l), l_r (3,j,l)


              write(i_ftls, "(f13.6,3f10.2)")  crys%ref_l_alpha (j,l), crys%ref_l_r (1,j,l),crys%ref_l_r (2,j,l), &
                                                      crys%ref_l_r (3,j,l)
              write(i_ftls, "(a, 6f10.2)") "FW ",r_b11 (j,l) , r_b22 (j,l) , r_b33 (j,l) , &
                                      r_b12 (j,l) ,r_b31 (j,l) , r_b23 (j,l)
            end do
          end do
          write(i_ftls,"(a)")              "  "
          write(i_ftls,"(a)")     " CALCULATION  "
          if (opt == 0) then
            write(i_ftls,"(a)")          " SIMULATION"
            write(i_ftls,"(3f10.4)")     th2_min, th2_max, d_theta
          elseif (opt == 3) then
            write(i_ftls,"(2a)")          " LOCAL_OPTIMIZER   ", opti%method
            write(i_ftls,"(a,i4)")          " MXFUN  ", opti%mxfun
            write(i_ftls,"(a,f10.4)")          " EPS  ", opti%eps
            write(i_ftls,"(a, i2)")          " IOUT  ", opti%iout
            write(i_ftls,"(a,f10.7)")          " ACC  ", opti%acc
          elseif (opt == 4) then
            write(i_ftls,"(a)")          " LMQ"
            if (Cond%constr ) write(i_ftls,"(a,f5.2)")          " BOXP    " , Cond%percent
            write(i_ftls,"(a,i4)")    " CORRMAX    ", cond%corrmax
            write(i_ftls,"(a,i4)")    " MAXFUN     ", cond%icyc
            write(i_ftls,"(a,g14.6)")    " TOL     ", cond%tol
            write(i_ftls,"(a,i2)")    " Nprint     ", cond%nprint
          else
            write(*,*) "ERROR writing *.ftls file: Problem with calculation section"
            return
          end if

          if(opt == 3 .or. opt == 4) then
            write(i_ftls,"(a)")              "  "
            write(i_ftls,"(a)")          " EXPERIMENTAL"
            write(i_ftls,"(2a)")         " FILE  ", dfile
            if (nexcrg /= 0) then
              write(i_ftls,"(a, i2)")    " EXCLUDED_REGIONS  ",  nexcrg
              do i=1,nexcrg
                write(i_ftls,"(2f10.4)")  alow(i),ahigh(i)
              end do

            end if
            write(i_ftls,"(2a)")         " FFORMAT  ",fmode
            write(i_ftls,"(2a)")         " BGR    ",background_file
            write(i_ftls,"(2a)")         " BCALC    ",mode
          end if

          return

       End Subroutine Write_ftls


       Subroutine Write_Prf(diff_pat,i_prf)
          !-----------------------------------------------
          !   D u m m y   A r g u m e n t s
          !-----------------------------------------------
          Type(Diffraction_Pattern_Type), intent(in) :: diff_pat
          integer,                        intent(in) :: i_prf
          !-----------------------------------------------
          !   L o c a l   V a r i a b l e s
          !-----------------------------------------------
          integer ::  i, j, iposr, ihkl, irc, nvk

          real :: twtet, dd, scl,yymi,yyma
          character (len=1)   :: tb
          character (len=50)  :: forma1,forma2
          !character (len=200) :: cell_sp_string
          !-----------------------------------------------
          !check for very high values of intensities and rescal everything in such a case
          ! scl: scale factor scl=1.0 (normal ymax < 1e6, 0.1 multiplier)
          yyma=diff_pat%ymax
          scl=1.0
          do
            if(yyma < 1.0e6) exit !on exit we have the appropriate value of scl
            scl=scl*0.1
            yyma=yyma*scl
          end do
          yymi=diff_pat%ymin*scl
          tb=CHAR(9)

          if(yyma < 100.0) then
           forma1='(f12.4,4(a,f8.4))'
          else if(yyma < 1000.0) then
           forma1='(f12.4,4(a,f8.3))'
          else if(yyma < 10000.0) then
           forma1='(f12.4,4(a,f8.2))'
          else if(yyma < 100000.0) then
           forma1='(f12.4,4(a,f8.1))'
          else
           forma1='(f12.4,4(a,f8.0))'
          end if
          !cell_sp_string=" "
          !write(unit=cell_sp_string,fmt="(a,3f10.5,3f10.4,a)")"  CELL: ",cellp(1)%cell(:),cellp(1)%ang(:),"   SPGR: "//symb(1)
          write(i_prf,'(A)') trim(diff_pat%title)

          write(i_prf,'(i3,i7,5f12.5,i5)') 1,diff_pat%npts,lambda,lambda2,0.0,0.0,0.0,0

          nvk=0
          WRITE(i_prf,'(17I5)') n_hkl, nvk , nexcrg

          do  j=1,nexcrg
            write(i_prf,'(2f14.5)')alow(j),ahigh(j)
          end do

          WRITE(i_prf,'(15a)')' 2Theta',tb,'Yobs',tb,'Ycal',tb,  &
                'Yobs-Ycal',tb,'Backg',tb,'Posr',tb,'(hkl)',tb,'K'

          do  i=1,diff_pat%npts
            twtet=diff_pat%x(i)
            dd=(diff_pat%y(i)-diff_pat%ycalc(i))*scl
            do  j=1,nexcrg
              if( twtet >= alow(j) .AND. twtet <= ahigh(j) ) dd=0.0
            end do
            WRITE(i_prf,forma1) twtet,tb,diff_pat%y(i)*scl,tb,diff_pat%ycalc(i)*scl,tb,  &
                dd-yyma/4.0,tb,diff_pat%bgr(i)*scl-yymi/4.0
          end do

          !Reflections
          iposr=0
          irc=1
          ihkl=0
          DO i=1,n_hkl
            WRITE(i_prf,'(f12.4,9a,i8,a,3i3,a,2i3)')  &
                dos_theta(i),tb,'        ',tb,'        ',tb,'        ',  &
                tb,'        ',tb,iposr, tb//'(',hkl_list(:,i),')'//tb,ihkl,irc
          END DO

          RETURN
       End Subroutine Write_Prf

       Subroutine scale_factor_lmq(pat, fvec, chi2,r)
         type (diffraction_pattern_type), intent(in out):: pat
         Real (Kind=cp),Dimension(:),     Intent(in Out):: fvec
         real,                            Intent(   out):: chi2, r
         real                                           :: a ,  c
         integer                                        :: punts=0
         integer                                        :: i,j

         a=0
         punts=0
         do_a: Do j = 1, pat%npts
           do i=1,nexcrg
            if(pat%x(j) >= alow(i) .and. pat%x(j) <= ahigh(i)) cycle do_a
           end do
            punts=punts+1
            pat%ycalc(j)  = brd_spc(j)
            a = a + pat%ycalc(j)
         End do do_a

         c=0
         do_c:Do  j = 1, pat%npts
           do i=1,nexcrg
            if(pat%x(j) >= alow(i) .and. pat%x(j) <= ahigh(i)) cycle do_c
           end do
           c = c + (pat%y(j) - pat%bgr(j))
         End do do_c

         pat%scal = c/a

         Do j = 1, pat%npts
          pat%ycalc(j) = pat%scal * pat%ycalc(j)+ pat%bgr(j)
         End do
         call calc_par_lmq(pat, punts, r,fvec, chi2)

         return
       End subroutine scale_factor_lmq

       Subroutine calc_par_lmq (pat, punts, r, fvec, chi2)
         type (diffraction_pattern_type), intent(in    ) :: pat
         integer,                         intent(in    ) :: punts
         real,                            intent(   out) :: r
         Real (Kind=cp),Dimension(:),     Intent(in out) :: fvec
         real,                            Intent(   out) :: chi2
         !----
         real                                            :: a,b,c,delta
         integer                                         :: j,i

         a=0.0
         b=0.0
         do_a: do j=1, pat%npts
           do i=1,nexcrg
            if(pat%x(j) >= alow(i) .and. pat%x(j) <= ahigh(i)) cycle do_a
           end do
           a= a + pat%y(j)
           b= b + abs(pat%y(j) - pat%ycalc(j))
         end do do_a

         r =  b/a *100.0
         c=0.0
         do_c: do j=1, pat%npts
           do i=1,nexcrg
            if(pat%x(j) >= alow(i) .and. pat%x(j) <= ahigh(i)) cycle do_c
           end do
           fvec(j)=  (pat%y(j) - pat%ycalc(j))/sqrt(pat%sigma(j))
           c = c + fvec(j)*fvec(j)
         end do do_c

         chi2= c/(punts-opti%npar)
         return
       End subroutine calc_par_lmq

  End module Dif_compl
!________________________________________________________________________________________________________________
  Module dif_ref
    use CFML_GlobalDeps,            only : sp , cp
    use CFML_String_Utilities,      only : number_lines , reading_lines ,  init_findfmt, findfmt ,iErr_fmt, getword, &
                                           err_string, err_string_mess, getnum, Ucase

    use CFML_Crystal_Metrics,       only : Set_Crystal_Cell, Crystal_Cell_Type
    use CFML_Diffraction_patterns , only : diffraction_pattern_type
    use CFML_Optimization_LSQ,      only : Levenberg_Marquardt_Fit
    use CFML_LSQ_TypeDef,           only : LSQ_Conditions_type, LSQ_State_Vector_Type
    use CFML_Math_General,          only : spline, splint, sind, cosd
    use diffax_mod
    use read_data,                  only : crys, read_structure_file, length,   opti , cond, vs
    use diffax_calc,                only : salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, optimz,point,  &
                                           gospec, gostrk, gointr,gosadp, getfnm, nmcoor
    use Dif_compl,                  only : scale_factor_lmq, Write_Prf, Faults2diffax, vs2faults, Var_assign

    implicit none

    public  ::  Cost_LMQ, apply_aberrations
    type (diffraction_pattern_type), save :: difpat
    integer, parameter, public            :: iout = 25,i_out=20


    contains


!! Subroutine Levenberg_Marquard_Fit(Model_Functn, m, c, Vs, chi2, infout,residuals)      Cost_LMQ, Nop, Cond, Vs, chi2, texte
!!--..            Integer,                     Intent(In)      :: m        !Number of observations
!!--..            type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
!!--..            type(LSQ_State_Vector_type), Intent(In Out)  :: Vs       !State vector
!!--..            Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
!!--..            character(len=*),            Intent(out)     :: infout   !Information about the refinement (min length 256)
!!--..            Real (Kind=cp), dimension(:),optional, intent(out) :: residuals
!!--..         End Subroutine

!!Interface No_Fderivatives
!!--..           Subroutine Model_Functn(m, n, x, fvec, iflag)             !Model Function subroutine
!!--..             Use CFML_GlobalDeps, Only: cp
!!--..             Integer,                       Intent(In)    :: m, n    !Number of observations and free parameters
!!--..             Real (Kind=cp),Dimension(:),   Intent(In)    :: x       !Array with the values of free parameters: x(1:n)
!!--..             Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec    !Array of residuals fvec=(y-yc)/sig : fvec(1:m)
!!--..             Integer,                       Intent(In Out):: iflag   !If iflag=1 calculate only fvec without changing fjac
!!--..           End Subroutine Model_Functn                               !If iflag=2 calculate only fjac keeping fvec fixed
!!--..         End Interface No_Fderivatives

    Subroutine Cost_LMQ(m,npar,v,fvec,iflag)                 !Levenberg Marquardt
      Integer,                       Intent(In)    :: m      !is the number of observations
      Integer,                       Intent(In)    :: npar   !is the number of free parameters
      Real (Kind=cp),Dimension(:),   Intent(In)    :: v      !List of free parameters values
      Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec   !Residuals
      Integer,                       Intent(In Out):: iflag  !=0 for printing, 1: Full calculation, 2: Calculation of derivatives
                                                             ! If modified to a negative number the algorithm stops.
      !local variables
      logical                  :: ok
      integer                  :: j ,i, k, a, b
      real, dimension(max_npar):: shift, state
 !     real, dimension(m),save  :: svdfvec
      real                     :: rf
      real, save               :: chi2
      integer, save            :: iter=0

      shift(1:npar) = v(1:npar) - vector(1:npar)

      Select Case(iflag)

        Case(1)  !Calculation of fvec and updating completely the parameters

           fvec=0.0
           write(i_out,"(a)")" --------FCOST-------"
           do i = 1, crys%npar
             state(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
             write(i_Out,"(a,i3,a,4f14.5,i5)") " State(",i,"):"//namepar(i),state(i), crys%list(i) ,  mult(i), shift(crys%p(i)),crys%p(i)
           end do
           crys%list(:) = state(:)
           vector(1:npar) = v(1:npar) !vector upload
           crys%Pv_refi(1:npar) = v(1:npar)
           call Pattern_Calculation(state,ok)
           if(.not. ok) then
             print*, "Error calculating spectrum, please check input parameters"
           else
             call scale_factor_lmq(difpat, fvec, chi2,rf)
           end if
           numcal = numcal + 1  !Counter for optimz, detun, etc
           iter = iter + 1
           write(*,"(a,i4,2(a,f14.4))")  " => Iteration ",iter,"   R-Factor = ",rf,"   Chi2 = ",chi2
           write(i_out,"(a,i4,2(a,f14.4))")  " => Iteration ",iter,"   R-Factor = ",rf,"   Chi2 = ",chi2

        Case(2)  !Calculation of numerical derivatives

           do i = 1, crys%npar
             state(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
           end do
           call Pattern_Calculation(state,ok)
           if(.not. ok) then
             print*, "Error calculating spectrum, please check input parameters"
           else
             call scale_factor_lmq(difpat, fvec, chi2, rf)
           end if

        Case(0)  !Printing
           write(*,"(a,i4,a,f14.4)")  " => Iteration ",iter,"   Chi2 = ",chi2
           do i=1, crys%npar
             write(*,"(a,f14.5)")  "  ->  "//namepar(i), state(i)
           end do

      End Select

      ok = .true.

      IF(cfile) CLOSE(UNIT = cntrl)
      return

    End subroutine Cost_LMQ

    Subroutine Pattern_Calculation(state,ok)
      real, dimension(:), intent(in) :: state
      logical,            intent(out):: ok


      call Var_assign(state)
      ok = .true.
      if ((conv_d == 1 .or.  numcal== 0) .and. ok ) ok = get_g()
      if ((conv_d == 1 .or.  numcal== 0  .or. conv_e==1) .and. (ok .AND. rndm)) ok = getlay()
      if ((conv_c == 1 .or.  numcal== 0) .and. ok ) CALL sphcst()
      if ( numcal == 0 .and. ok ) CALL detun()
      if ((conv_b == 1 .or.  conv_c == 1 .or. conv_d == 1 .or. conv_e==1   .or. &
           numcal == 0 .or.  conv_f==1)  .and. ok ) CALL optimz(infile, ok)

      IF(.NOT. ok) then
        IF(cfile) CLOSE(UNIT = cntrl)
        return
      END IF
      CALL gospec(infile,outfile,ok)
      return

      call Apply_Aberrations()

    End Subroutine Pattern_Calculation

    Subroutine Apply_Aberrations()
      !--- Modifies brd_spc by spline interpolation after applying zero-shift
      !--- displacement and transparency (Bragg-Brentano)
      real(kind=cp), dimension(n_high) :: true_2th, broad_spect, der2v
      integer       :: i
      real(kind=cp) :: aux,tt,shift,ycal,t

      do i=1,n_high
         true_2th(i)=thmin+real(i-1,kind=cp)*step_2th
         broad_spect(i) = brd_spc(i)  !needed because brd_spc is double precision
      end do
      aux=1.0e+35

      call spline(true_2th,broad_spect,n_high,aux,aux,der2v)

      ! Shift(2Theta) = zero + SyCos * cos(Theta) + SySin * sin(2Theta)     !Bragg-Brentano
      ! Shift(2Theta) = zero + SyCos * cos(2Theta) + SySin * sin(2Theta)    ! Debye-Scherrer
      do i=1,difpat%npts
        tt=difpat%x(i)
        t=tt
        if(i_geom == 0) t=t*0.5 !Bragg Brentano (SyCos: displacement, Sysin: transparency)
        shift=crys%zero_shift + crys%sycos * cosd(t) + crys%sysin * sind(tt)
        tt=tt-shift
        call splint(difpat%x,broad_spect,der2v,difpat%npts,tt,ycal)
        difpat%ycalc(i)=ycal
        brd_spc(i)=ycal
      end do
      return
    End Subroutine Apply_Aberrations

   End module dif_ref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   PROGRAM FAULTS

     use CFML_GlobalDeps,              only : sp , cp
     use CFML_String_Utilities,        only : number_lines , reading_lines ,  init_findfmt, findfmt ,iErr_fmt, &
                                              getword, err_string, err_string_mess, getnum, Ucase, lcase
     use CFML_Diffraction_patterns,    only : read_pattern , diffraction_pattern_type , err_diffpatt, err_diffpatt_mess,  &
                                              read_background_file
     use CFML_Optimization_General,    only : Opt_Conditions_Type, Local_Optimize
     use CFML_Crystal_Metrics,         only : Set_Crystal_Cell, Crystal_Cell_Type
     use CFML_Optimization_LSQ,        only : Levenberg_Marquardt_Fit,Info_LSQ_LM
     use CFML_LSQ_TypeDef,             only : LSQ_Conditions_type, LSQ_State_Vector_Type
     use diffax_mod
     use read_data,                    only : read_structure_file, length,   &
                                              crys, opti, cond, Vs, Err_crys, Err_crys_mess
     use diffax_calc ,                 only : salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, optimz,point,  &
                                              gospec, gostrk, gointr,gosadp, chk_sym, get_sym, overlp, nmcoor , getfnm
     use Dif_compl,                    only : scale_factor_lmq, Write_Prf, write_ftls, Faults2diffax, vs2faults, Var_assign
     use dif_ref

     implicit none


      real                    :: rpl,  theta , ymax, ymini , ymin ,deg
      LOGICAL                 :: ok, ending , gol, p_ok, arggiven=.false.
      INTEGER                 :: i ,n ,j, l , ier , fn_menu,a,b,c ,aa,bb,cc, e, narg
      character(len=100)      :: pfile, bfile , bmode
      character(len=100)      :: pmode, filenam
      character(len=10)       :: time, date
      Real (Kind=cp)          :: chi2     !final Chi2
      character(len=3000)     :: infout   !Information about the refinement (min length 256)

      pi = four * ATAN(one)
      pi2 = two * pi
      deg2rad = pi / one_eighty
      rad2deg = one_eighty / pi

      ending = .false.
      ok = .true.
      sfname = 'data.sfc'
      cntrl=ip


      !---- Arguments on the command line ----!
      narg=command_argument_count()

      if (narg > 0) then
         call get_command_argument(1,infile)
         arggiven=.true.
      end if

      CALL salute()

      if(.not. arggiven) then
        write(unit=op,fmt="(a)",advance="no") ' => Enter the complete name of the structure input file: '
        read(unit= *,fmt="(a)") infile
      end if
      !WRITE(op,fmt=*) "=> Looking for scattering factor data file '",  sfname(:),"'"
      OPEN(UNIT = sf, FILE = sfname)
      !WRITE(op,fmt=*) "=> Opening scattering factor data file '",  sfname(:),"'"

      filenam = trim(infile(1:(index(infile,'.')-1)))

      open(unit=i_out, file=trim(filenam)//".out",status="replace",action="write")
      CALL salute(i_out)


      call read_structure_file(infile, gol)

      if (err_crys) then
        write(unit=*,fmt="(a)") " ERROR in "//trim(infile)//": "//trim(err_crys_mess)
        stop
      else
        write(op, fmt=*) "=> Structure input file read in"
      end if

      call faults2diffax(crys)

      IF(gol) then
        ok = sfc()
      else
        GO TO 999
      end if


      do i = 1, crys%npar              ! to avoid repetitions
             if (index (namepar(i) , 'pos' ) == 1 )  conv_a = 1
             if (index (namepar(i) , 't' ) == 1 .or. index (namepar(i) , 'alpha' ) == 1 ) conv_b = 1
             if (index (namepar(i) , 'cell' ) == 1 ) conv_c = 1
             if (index (namepar(i) , 'alpha' ) == 1 ) conv_d = 1
             if (index (namepar(i) , 'num' ) == 1) conv_e = 1
             if (index (namepar(i) , 'Biso')==1) conv_f=1
             if (index (namepar(i) , 't' ) == 1 .or. index (namepar(i) , 'alpha' ) == 1  .or. &
                 index (namepar(i) , 'num' ) == 1 .or. index (namepar(i) , 'Biso')==1 .or.       &
                 index (namepar(i) , 'pos')==1 .or. index (namepar(i) , 'cell')==1 ) conv_g = 1 !used in diffax_calc not to recalculate the whole spectra if only instrumental parameters are refined
      end do

      !Simulation or optimization

      if(opt /= 0) then
         !Reading experimental pattern and scattering factors:
         write(*,"(a)") " => Reading Pattern file="//trim(dfile)
         call  read_pattern(dfile,difpat,fmode )
           if(Err_diffpatt) then
             print*, trim(err_diffpatt_mess)
           else
             if (th2_min <= 0.0000001 .and.  th2_max <= 0.0000001  .and. d_theta <= 0.0000001) then
                thmin    =  difpat%xmin    ! if not specified in input file, take the values from the observed pattern
                thmax    =  difpat%xmax
                step_2th =  difpat%step
                n_high = nint((thmax-thmin)/step_2th+1.0_cp)
                th2_min = thmin * deg2rad
                th2_max = thmax * deg2rad
                d_theta = half * deg2rad * step_2th
             end if
           end if
         write(*,"(a)") " => Reading Background file="//trim(background_file)
         call read_background_file(background_file, mode ,difpat)
               if(Err_diffpatt) print*, trim(err_diffpatt_mess)

      end if

      !Operations before starting :
      check_sym = .false.
      IF(symgrpno == UNKNOWN) THEN
        symgrpno = get_sym(ok)
        IF(.NOT. ok) GO TO 999
        WRITE(op,"(a)") 'Diffraction point symmetry is '//pnt_grp
        IF(symgrpno /= 1) THEN
          WRITE(op,"(a,i6)") '  to within a tolerance of one part in ',  &
              nint(one / tolerance)
        END IF
      ELSE
        check_sym = .true.
        CALL chk_sym(ok)
        IF(.NOT. ok) GO TO 999
      END IF
      IF(ok) ok = get_g()
      IF(ok .AND. rndm) ok = getlay()
      IF(ok) CALL sphcst()
      IF(ok) CALL detun()
      IF(.NOT. ok) GO TO 999
      ! See if there are any optimizations we can do
      IF(ok) CALL optimz(infile, ok)
      !write(*,"(a)") " => Calling dump file: "//trim(filenam)//".dmp"
      call dump(infile, i_out, p_ok) !Writing read control file
      call overlp()
      call nmcoor ()

       Select  case (opt)

          Case (0)

          ! What type of intensity output does the user want?
           10  IF(ok) THEN
              ! WRITE(op,"(a)") 'Enter function number: '
              ! WRITE(op,"(a)") '0 POINT, 1 STREAK, 2 INTEGRATE, 3 POWDER PATTERN, 4 SADP : '
              ! READ(cntrl,*,END=999) n
                n=3 !Powder pattern
             END IF

            ! Do what the user asked for.
             IF(ok) THEN

                IF(n == 0) THEN
                   CALL point(ok)
                ELSE IF(n == 1) THEN
                   CALL gostrk(infile,outfile,ok)
                ELSE IF(n == 2) THEN
                   CALL gointr(ok)
                ELSE IF(n == 3) THEN
                   write(*,*) "calculating powder diffraction pattern"
                   CALL gospec(infile,outfile,ok)

                     Do j = 1, n_high
                         ycalcdef(j) = brd_spc(j)
                     end do

                     CALL getfnm(filenam, outfile, '.dat', ok)

                     OPEN(UNIT = iout, FILE = outfile, STATUS = 'new')
                        write(unit = iout,fmt = *)'!', outfile
                        write(unit = iout,fmt = '(3f12.4)')thmin, step_2th,thmax
                     !  theta = thmin +(j-1)*d_theta
                        write(unit = iout,fmt = '(8f12.2)') ( ycalcdef(j), j=1, n_high    )

                    CLOSE(UNIT = iout)
                    ok = .true.
                ELSE IF(n == 4) THEN
                   CALL gosadp(infile,outfile,ok)
                ELSE
                   WRITE(op,"(a)") ' Unknown function type.'
                END IF

              END IF


              IF(ok .AND. n /= 3) THEN
              96   WRITE(op,"(a)") ' => Enter 1 to return to function menu.'
                   READ(cntrl,*,ERR=96,END=999) fn_menu
                   IF(fn_menu == 1) GO TO 10
              END IF


          Case (4) !LMQ

            chi2o = 1.0E10                         !initialization of agreement factor
          !  rpl = 0
            do i=1, opti%npar                       !creation of the step sizes
              vector(i) = crys%Pv_refi(i)
            end do

            call Levenberg_Marquardt_Fit(cost_LMQ, difpat%npts, cond, Vs, chi2, infout)

            !Output the final list of refined parameters
            Call Info_LSQ_LM(Chi2,op,Cond,Vs)
            Call Info_LSQ_LM(Chi2,i_out,Cond,Vs)

            write(*,"(a)") " => "// trim(infout)

            CALL getfnm(filenam, outfile, '.prf', ok)
              if (ok) then
                OPEN(UNIT = iout, FILE = outfile, STATUS = 'replace')
                call Write_Prf(difpat,iout)
              else
                write(*,*) 'The outfile .prf cannot be created'
              end if
              !CALL getfnm(trim(filenam)//"_new", outfile, '.ftls', ok)
             if (ok) then
               filenam=trim(filenam)
               i=index(filenam,"_new")
               if(i /= 0) filenam=filenam(1:i-1)
               OPEN(UNIT = i_flts, FILE = trim(filenam)//"_new.flts", STATUS = 'replace',action="write")
               call Write_ftls(crys,i_flts, vs)
             else
               write(*,*) 'The outfile .flts cannot be created'
             end if


!-------------------------------------------------------------------------------------------------------------------------------

            Case default

                print*,"problems reading mode "


          End select


      ending = .true.

      999 IF(cfile) CLOSE(UNIT = cntrl)
      IF(ok .AND. ending) THEN
        WRITE(op,"(a)") ' => FAULTS ended normally....'
      ELSE
        WRITE(op,"(a)") ' => FAULTS was terminated abnormally!'
      END IF

   END PROGRAM FAULTS


    Subroutine Write_FST(fst_file,v,cost)
       Use CFML_GlobalDeps, only: Cp
       character(len=*),     intent(in):: fst_file
       real(kind=cp),dimension(:),    intent(in):: v
       real(kind=cp),                 intent(in):: cost
       return
    End Subroutine Write_FST

