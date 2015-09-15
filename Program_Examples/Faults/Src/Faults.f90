
  Module Dif_compl

      use diffax_mod
      use CFML_GlobalDeps,                only : sp,dp, to_deg
      use CFML_Diffraction_Patterns ,     only : diffraction_pattern_type
      use CFML_Optimization_General,      only : Opt_Conditions_Type
      use CFML_Math_General,              only : sind, asind, Euclidean_norm, acosd
      use CFML_String_Utilities,          only : pack_string,Get_Separator_Pos,u_case,l_case
      use read_data
      use CFML_LSQ_TypeDef,               only : LSQ_State_Vector_Type
      use CFML_Crystal_Metrics,           only : Crystal_Cell_Type, Set_Crystal_Cell
      use CFML_Crystallographic_Symmetry, only : Space_Group_Type, Set_SpaceGroup
      use CFML_IO_Formats,                only : Write_Cif_Template, Write_CFL
      use CFML_Atom_TypeDef,              only : Atom_List_Type, Allocate_Atom_List


      implicit none

      private

      !public subroutines
      public :: Write_Prf, Write_ftls, Faults2diffax, vs2faults, Var_assign, calc_avcell,Write_FST_VESTA

      contains
!________________________________________________________________________________________________________________________
!   This subroutine copy the "state" or "cys%List" vector to DIFFAX variables
    Subroutine Var_assign(state)

      real, dimension(:), intent(in) :: state
      !local variables
      integer :: i, j, k, a, b

          !////////CALCULATED PATTERN VARIABLES ASSIGNMENT////////////////////////////////////

      do i=1, crys%npar
        if (namepar(i) ==  'Scale_Factor') crys%patscal = state(i)
        if (namepar(i) ==  'U')          pv_u   = state(i)
        if (namepar(i) ==  'V')          pv_v   = state(i)
        if (namepar(i) ==  'W')          pv_w   = state(i)
        if (namepar(i) ==  'x')          pv_x   = state(i)
        if (namepar(i) ==  'Dg')         pv_dg  = state(i)
        if (namepar(i) ==  'Dl')         pv_dl  = state(i)
        if (namepar(i) ==  'cell_a')     cell_a = state(i)
        if (namepar(i) ==  'cell_b')     cell_b = state(i)
        if (namepar(i) ==  'cell_c')     cell_c = state(i)
        if (namepar(i) ==  'num_layers') l_cnt  = state(i)
        if (namepar(i) ==  "diameter_a") Wa = state(i)
        if (namepar(i) ==  "diameter_b") Wb = state(i)
        if (namepar(i) ==  'zero_shift') crys%zero_shift  = state(i)
        if (namepar(i) ==  'sycos')      crys%sycos  = state(i)
        if (namepar(i) ==  'sysin')      crys%sysin  = state(i)
        if (index( namepar(i) ,'cell_gamma') == 1)   cell_gamma = state(i)*deg2rad
        if (index (namepar(i), 'ChebCoeff_' ) == 1)     then
            read (unit = namepar(i)(11:12), fmt = "(i2)" ) a
            crys%chebp(a)  = state(i)
        end if
        if (index (namepar(i), 'Bkg_Scale_' ) == 1)     then
            read (unit = namepar(i)(11:12), fmt = "(i2)" ) a
            crys%bscalpat(a)  = state(i)
        end if

        do j=1, n_layers
          do k=1, n_atoms
            if (index (namepar(i) , 'pos_x' )== 1)     then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) a,b
                a_pos(1,a,b)  = state(i) * pi2       !need to invert conversion done by routine nmcoor (diffax_calc)
            end if
            if (index (namepar(i) ,'pos_y' )== 1)    then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) a,b
                a_pos(2,a,b)  = state(i) * pi2
            end if
            if (index (namepar(i), 'pos_z' ) == 1 )   then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) a,b
                a_pos(3,a,b)  = state(i) * pi2
            end if
            if (index (namepar(i),'Biso')==1) then
                read (unit = namepar(i)(5:8), fmt = "(2i2)" ) a,b
                a_b(a,b)  = state(i)
            end if
            if (index (namepar(i),'Occ')==1) then
                read (unit = namepar(i)(4:7), fmt = "(2i2)" ) a,b
                a_occup(a,b)  = state(i)
            end if
            if (index( namepar(i) ,  'alpha' ) == 1)    then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) b,a
                l_alpha(a,b)  = state(i)

            end if
            if (index (namepar(i), 'tx' )== 1 )    then
                read (unit = namepar(i)(3:6), fmt = "(2i2)" ) b,a
                l_r(1,a,b)  = state(i)
            end if
            if (index (namepar(i), 'ty' )== 1 )    then
                read (unit = namepar(i)(3:6), fmt = "(2i2)" ) b,a
                l_r(2,a,b)  = state(i)
            end if
            if (index (namepar(i), 'tz' ) == 1)     then
                read (unit = namepar(i)(3:6), fmt = "(2i2)" ) b,a
                l_r(3,a,b)  = state(i)
            end if
            if (index (namepar(i), 'FW_11' ) == 1)     then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) b,a
                r_b11(a,b)  = state(i)
            end if
            if (index (namepar(i), 'FW_22' ) == 1)     then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) b,a
                r_b11(a,b)  = state(i)
            end if
            if (index (namepar(i), 'FW_33' ) == 1)     then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) b,a
                r_b11(a,b)  = state(i)
            end if
            if (index (namepar(i), 'FW_12' ) == 1)     then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) b,a
                r_b11(a,b)  = state(i)
            end if
            if (index (namepar(i), 'FW_31' ) == 1)     then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) b,a
                r_b11(a,b)  = state(i)
            end if
            if (index (namepar(i), 'FW_23' ) == 1)     then
                read (unit = namepar(i)(6:9), fmt = "(2i2)" ) b,a
                r_b11(a,b)  = state(i)
            end if
          end do
        end do
      end do

    End subroutine Var_assign

    Subroutine vs2faults(vs, crys)
       !This subroutine calculates the final values of crys before writing the new flts-file
       type(crys_2d_type),          intent(in out) :: crys
       type(LSQ_State_Vector_type), Intent(In)     :: Vs

       !local variables
       integer                  :: n, j, i
       real, dimension(max_npar):: shift

       shift(1:opti%npar) = vs%pv(1:opti%npar) - crys%Pv_refi(1:opti%npar)

       do i = 1, crys%npar
           crys%list(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
       end do

       call Var_assign(crys%list)

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
             r_b13(1:j,1:n) = crys%r_b13 (1:j,1:n)
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

     subroutine calc_avcell(nlay_avcell,lay_avcell,verif)
      integer , intent(in)       ::    nlay_avcell    !number of layers needed for calculating the average cell
      integer , dimension(:), intent(in)       ::    lay_avcell     !sequence of layers to be stacked for calculating the average cell
      logical , intent(out)      ::    verif


      real, dimension(3)   :: c_avcell               !c vector of the average cell
      real, dimension(3)   :: c_avnorm               !to calculate the norm of c
      real, dimension(3)   :: unita, unitb
      real, dimension(3,3) :: Metric_Tensor
      character(len=132)   :: txt,key
      integer              :: i,j,k,l,m, ier

              verif=.true.
              unita=(/ 1, 0, 0 /)
              unitb=(/ 0, 1, 0 /)
              c_avcell=0
              Metric_tensor= 0.0
              Metric_Tensor(1,1)=cell_a * cell_a
              Metric_Tensor(2,2)=cell_b * cell_b
              Metric_Tensor(3,3)=cell_c * cell_c
              Metric_Tensor(1,2)=cell_a * cell_b *cos(cell_gamma)
              Metric_Tensor(2,1)=Metric_Tensor(1,2)

              do i=1,nlay_avcell
                  l=lay_avcell(i)                        !layer origine
                  j=lay_avcell(i+1)                      !layer destination

   ! write(unit=*,fmt="(a,i2,a,i2,a,i2)") "For average cell, vector ", i, " starts from layer ", l, "goes to layer", j

                  c_avcell=c_avcell+l_r(:,j,l)     !calculate the vector C of the average cell in the FAULTS cell

   ! write(unit=*,fmt="(a,f12.5,f12.5,f12.5)") "The average c vector is now ",  c_avcell(1), c_avcell(2), c_avcell(3)

               end do

              write(unit=*,fmt="(a,f12.5,f12.5,f12.5)") " => The average c vector is ",  c_avcell

              c_avnorm=(/ c_avcell(1)*cell_a, c_avcell(2)*cell_b, c_avcell(3)*cell_c /)

   ! write(unit=*,fmt="(a,f12.5,f12.5,f12.5)") "(for norm calculation) ",  c_avnorm(1), c_avnorm(2), c_avnorm(3)

   ! write(unit=*,fmt="(a)") "Now calculate the cell parameters of the average cell"

                                                       !now calculate the cell parameters of the average cell
             aver_cell(1)=cell_a                       !a from DIFFAX variable

   ! write(unit=*,fmt="(a,f12.5)") "a=", aver_cell(1)

             aver_cell(2)=cell_b                       !b from DIFFAX variable

   ! write(unit=*,fmt="(a,f12.5)") "b=", aver_cell(2)

             !aver_cell(3)=Euclidean_norm(3,c_avnorm)   !c
             !aver_cell(3)=sqrt(dot_product(c_avnorm,c_avnorm)+ 2.0*c_avnorm(1)*c_avnorm(2)*cos(cell_gamma))
             aver_cell(3)=sqrt(dot_product(c_avcell,matmul(Metric_Tensor,c_avcell)))

   ! write(unit=*,fmt="(a,f12.5)") "c=", aver_cell(3)

             aver_ang(1)= &
             !acosd(dot_product(c_avcell,unitb)/(Euclidean_norm(3,c_avcell)))    !alpha
             acosd(dot_product(unitb,Matmul(Metric_Tensor,c_avcell))/(aver_cell(2)*aver_cell(3)))    !alpha


   ! write(unit=*,fmt="(a,f12.5)") "alpha=", aver_ang(1)

             aver_ang(2)= &
             !acosd(dot_product(c_avcell,unita)/(Euclidean_norm(3,c_avcell)))    !beta
             acosd(dot_product(unita,Matmul(Metric_Tensor,c_avcell))/(aver_cell(1)*aver_cell(3)))    !alpha

   ! write(unit=*,fmt="(a,f12.5)") "beta=", aver_ang(2)

             aver_ang(3)=cell_gamma*rad2deg            !gamma from DIFFAX variable

   ! write(unit=*,fmt="(a,f12.5)") "gamma=", aver_ang(3)

         return
    end subroutine calc_avcell

    Subroutine Write_ftls(crys,i_ftls, vs)
       !-----------------------------------------------
       !   D u m m y   A r g u m e n t s
       !-----------------------------------------------
       Type(crys_2d_type),          intent(in  out) :: crys
       integer,                     intent(in)      :: i_ftls
       type(LSQ_State_Vector_type), intent(in)      :: Vs


       integer                            :: a,b,c, j, l,i ,n, print_width,m
       CHARACTER(LEN=80), dimension(2)    :: list

      !call vs2faults(vs, crys)
       call Restore_Codes() !the ref_* variables are now the original codes

       !Here we use the DIFFAX variables or the type Crys
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
       write(i_ftls,"(tr25,3f10.2)") ref_glb(2:4)

       if (crys%broad == ps_vgt .and. crys%trm) then
         write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dg         Dl"
         write(i_ftls,"(a,4f12.6,2f11.2, a)") " Pseudo-Voigt   ", pv_u, pv_v, pv_w, pv_x, pv_dg,pv_dl, " TRIM"
         write(i_ftls,"(tr16,4f12.2,2f11.2)")  ref_glb(5:10)
       elseif (crys%broad == ps_vgt .and. .not. crys%trm ) then
         write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dg         Dl"
         write(i_ftls,"(a,4f12.6,2f11.2)")   " Pseudo-Voigt   ", pv_u, pv_v, pv_w, pv_x, pv_dg, pv_dl
         write(i_ftls,"(tr16,4f12.2,2f11.2)")  ref_glb(5:10)
       elseif (crys%broad == pv_gss .and. crys%trm) then
         write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dg"
         write(i_ftls,"(a,4f12.6,f11.2, a)") "    Gaussian    ", pv_u, pv_v, pv_w, pv_x, pv_dg, "TRIM"
         write(i_ftls,"(tr16,4f12.2,f11.2)")  ref_glb(5:9)
       elseif (crys%broad == pv_gss .and. .not. crys%trm ) then
         write(i_ftls,"(a)")  "!instr. broadening       u           v           w           x         Dg"
         write(i_ftls,"(a,4f12.6,1f11.2)") "    Gaussian    ", pv_u, pv_v, pv_w, pv_x, pv_dg
         write(i_ftls,"(tr16,4f12.2,f11.2)")  ref_glb(5:9)
       elseif (crys%broad == pv_lrn .and. crys%trm ) then
         write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dl"
         write(i_ftls,"(a,4f12.6,1f11.2 a)") "   Lorentzian   ", pv_u, pv_v, pv_w, pv_x, pv_dl, "TRIM"
         write(i_ftls,"(tr16,4f12.2,f11.2)")  ref_glb(5:8), ref_glb(10)
       elseif   (crys%broad==pv_lrn .and. .not. crys%trm) then
         write(i_ftls,"(a)") "!instr. broadening       u           v           w           x         Dl"
         write(i_ftls,"(a,4f12.6,1f11.2)") "   Lorentzian   ", pv_u, pv_v, pv_w, pv_x, pv_dl
         write(i_ftls,"(tr16,4f12.2,f11.2)")  ref_glb(5:8), ref_glb(10)
       else
         write(*,*) "ERROR writing *.ftls file: Problem with instrumental parameters!"
         return
       end if


       write(i_ftls,"(a)")              "  "
       write(i_ftls,"(a)")          " STRUCTURAL  "
       if(aver_cellcalc) then
          write(i_ftls,"(/,a)")        "!    Layers to be stacked to calculate the average cell "
          write(i_ftls,"(a,i2)")       "CalcAverCell ", nlay_avcell
          write(i_ftls,"(25i3)")           lay_avcell(1:nlay_avcell)
       end if
       if(aver_cellgiven) then
          write(i_ftls,"(/,a)")        "!    Average cell parameters (for Bragg positions in PRF file) "
          write(i_ftls,"(a,6f12.5)")   "AverCell ", aver_cell, aver_ang
       end if
       if(aver_spggiven) then
          write(i_ftls,"(/,a)") "!    Average Space Group (for Bragg positions in PRF file) "
          write(i_ftls,"(a,/)")   "SpGR "//trim(spgsymb)
       end if

       if(fst_given) then
          write(i_ftls,"(/,a,/)") "!    FullProf Studio commands "
          do i=1,num_fst
            write(i_ftls,"(a)") "FST_CMD  "//trim(fst_cmd(i))
          end do
          write(i_ftls,"(a)") " "
       end if

       write(i_ftls,"(a)")          "!         a            b            c          gamma "
       write(i_ftls,"(a,3f12.6,f12.2)")   " Cell  ", cell_a, cell_b, cell_c, cell_gamma*rad2deg
       write(i_ftls,"(tr7,4f12.2)")  ref_glb(11:14)
       write(i_ftls,"(a)")        "!Laue symmetry"
       write(i_ftls,*)            "Symm  ", pnt_grp
       write(i_ftls,"(a)")        "!number of layer types"
       write(i_ftls,"(a,i2)")     " NLAYERS  ", n_layers
       write(i_ftls,"(a)")        "!layer width"
       if (crys%finite_width) then
         write(i_ftls,"(a,2f10.2)") " Lwidth",   Wa, Wb
         write(i_ftls,"(2f10.2)")    ref_glb(15:16)
       else
         write(i_ftls,"(a)")        " Lwidth  Infinite"
       end if

       b=1
       a=1
       c=0
       do b=1, n_layers
         write(i_ftls,"(a)")              "  "
         if (fundamental(b) .eqv. .false.) then
           write(i_ftls,"(a,i2,a,i2)") " LAYER",b," =", original(b)
         else
           c=c+1
           write(i_ftls,"(a, i2)")  " LAYER", b
           list(1) = 'None '
           list(2) = 'CentroSymmetric '
          !WRITE(dmp,100) 'symmetry = ', list(l_symmetry(i)+1)
           write(i_ftls,"(a)")        "!Layer symmetry   "
           write(i_ftls,"(2a)")      " LSYM   ", list(l_symmetry(c)+1)
           write(i_ftls,"(/a)")      "!Atom name   number   x         y         z         Biso      Occ  "
           do a=1, crys%l_n_atoms(c)
             write(i_ftls,"(2a,i9, 5f10.5)") " Atom ", a_name(a,c), a_number(a,c), a_pos(1, a,c)/pi2, &
                                        a_pos(2, a,c)/pi2,a_pos(3, a,c)/pi2, a_B (a,c), a_occup(a,c)
             write(i_ftls,"(tr19,5f10.2)") ref_atom(1:5,a,c)
           end do
         end if
       end do
       write(i_ftls,"(a)")              "  "
       write(i_ftls,"(a)")          " STACKING"
       write(i_ftls,"(a)")     "!stacking type"
       if (crys%xplcit) then
         write(i_ftls, "(a)") " EXPLICIT "
            if (semirandm) then
              write(i_ftls,"(a)") "!number of layers"
              write(i_ftls,"(a)", advance="no") " SEMIRANDOM "
              write(i_ftls,"(f6.1)") l_cnt
              write(i_ftls,"(a,i5,i5,i5,i5)") " SEQ", fls, lls, otls, stls
            elseif(spcfc) then
              write(i_ftls,"(a)") " SPECIFIC"
              write(i_ftls, "(25i2)") crys%l_seq(1:nint(crys%l_cnt))
            elseif(randm) then
             write(i_ftls,"(a)") "!number of layers"
              write(i_ftls,"(a)", advance="no") " RANDOM "
              write(i_ftls,"(f6.1)") l_cnt
            end if
       else
          write(i_ftls, "(a)") " RECURSIVE"
          write(i_ftls,"(a)") "!number of layers"
          if (crys%inf_thick) then
            write (i_ftls, "(a)") " INFINITE"
          else
            write (i_ftls, "( f6.1)") l_cnt
            write (i_ftls, "( f6.1,a)")  ref_glb(17)
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
           write(i_ftls, "(a, 4f12.6)")  "LT ",  l_alpha (j,l), l_r (1,j,l), l_r (2,j,l), l_r (3,j,l)
           write(i_ftls, "(f15.2,3f12.2)")  ref_trans(1:4, j, l)
           write(i_ftls, "(a, 6f12.6)") "FW ",r_b11 (j,l) , r_b22 (j,l) , r_b33 (j,l) , &
                                              r_b12 (j,l) , r_b13 (j,l) , r_b23 (j,l)
           write(i_ftls, "(f15.2, 5f12.2)") ref_trans(5:10, j, l)
         end do
       end do
       write(i_ftls,"(a)")              "  "
       write(i_ftls,"(a)")     " CALCULATION  "
       if (opt == 0) then
         write(i_ftls,"(a)")          " SIMULATION"
         if (funct_num == 3) then
           write(i_ftls,"(a)")          " ! Range of powder pattern:th2_min, th2_max, d_theta;   Scale_Factor and Background Level"
           write(i_ftls,"(a, 3f10.4,a,g14.5,a,f10.4)")  "POWDER", th2_min, th2_max, d_theta,"    ScaleF ",crys%patscal, "   Bckg_Level ",crys%Bckg_Level
         else
         	 write(i_ftls,"(a)")          " !Selected Area Diffraction Pattern: i_plane, l_upper, loglin, brightness"
         	 write(i_ftls,"(a, i2, f10.4, i2, f10.4)")  "SADP", i_plane, l_upper, loglin, brightness
         end if
       elseif (opt == 3) then
         write(i_ftls,"(2a)")          " LOCAL_OPTIMIZER   ", opti%method
         write(i_ftls,"(a,i4)")          " MXFUN  ", opti%mxfun
         write(i_ftls,"(a,f10.4)")          " EPS  ", opti%eps
         write(i_ftls,"(a, i2)")          " IOUT  ", opti%iout
         write(i_ftls,"(a,f10.7)")          " ACC  ", opti%acc
       elseif (opt == 4) then
         write(i_ftls,"(a)")          " LMA"
         if (Cond%constr) write(i_ftls,"(a,f5.2)")          " BOXP    " , Cond%percent
         write(i_ftls,"(a,i4)")    " Corrmax    ", cond%corrmax
         write(i_ftls,"(a,i4)")    " Maxfun     ", cond%icyc
         write(i_ftls,"(a,g14.6)")    " Tol     ", cond%tol
         write(i_ftls,"(a,i2)")    " Nprint     ", cond%nprint
       else
         write(*,*) "ERROR writing *.ftls file: Problem with calculation section"
         return
       end if
       if(replace_files) write(i_ftls,"(a)") " Replace_Files"

       if(opt == 3 .or. opt == 4) then
         write(i_ftls,"(a)")              "  "
         write(i_ftls,"(a)")          " EXPERIMENTAL"
         write(i_ftls,"(a)")          "!Filename                    Scale factor     code"
         write(i_ftls,"(2a, g14.6,f9.2)")         " File  ", dfile, crys%patscal, ref_glb(1)
         if (nexcrg /= 0) then
           write(i_ftls,"(a, i2)")    " Excluded_Regions  ",  nexcrg
           do i=1,nexcrg
             write(i_ftls,"(2f10.4)")  alow(i),ahigh(i)
           end do

         end if
         write(i_ftls,"(2a)")         " Fformat  ",fmode
         if (crys%bgrinter) then
           write(i_ftls,"(a)")"!Linear interpolation"
           write(i_ftls,"(2a)")         " BgrInter    ", background_file
         end if
         if (crys%bgrcheb) then
           write(i_ftls,"(a)") "!Polynomial  Number of coefficients"
           write(i_ftls,"(a, i2)")         " BgrCheb          ", crys%cheb_nump
           write(i_ftls,"(a)") "!Chebychev Polynomial coefficients"
           write(i_ftls,"(24f12.5)")  crys%chebp(1:crys%cheb_nump)
           write(i_ftls,"(24f12.2)")  ref_glb(18:17+crys%cheb_nump)
         end if

         write(i_ftls,"(a)")        "!Number of background patterns"  !NIK
         write(i_ftls,"(a, i2)")    " BgrNum ",  crys%num_bgrpatt     !

         if (crys%bgrpatt) then
           if(any(len_trim(crys%bfilehkl) /= 0)) then
             write(i_ftls,"(a)")                  "!Pattern file           Filename      Scale factor     code        hkl-file"
           else
             write(i_ftls,"(a)")                  "!Pattern file           Filename      Scale factor     code"
           end if
           do i=1, crys%num_bgrpatt
             write(i_ftls,"(2a, g18.5,f10.2,tr4,a)")  " BgrPatt    ",adjustr(crys%bfilepat(i)), crys%bscalpat(i), &
                                                    ref_glb(17+crys%cheb_nump+i), trim(crys%bfilehkl(i))
           end do
         end if

       end if

       return

    End Subroutine Write_ftls

    Subroutine Write_FST_VESTA()
      Type(Space_Group_Type) :: SpGr
      Type(Atom_List_Type)   :: Atm
      Type(Crystal_Cell_Type):: celld
      integer :: i,j,k,l,lact,n,m, nseq, nlayers,nl,nat,ncar,ier
      integer, dimension(1000)      :: seq
      integer, dimension(20)        :: pos
      real(kind=cp), dimension(3)   :: cvect,stck_vect,celax,celang
      real(kind=dp), dimension(6)   :: cell_fst
      real(kind=cp)                 :: caxis !Perpendicular to the plane of the layers
      real(kind=cp),     dimension(:,:),allocatable :: xyz
      real(kind=cp),     dimension(:),  allocatable :: biso,occ
      character(len=20), dimension(:),  allocatable :: atnam
      character(len=2)  :: elem
      character(len=16) :: aux
      character(len=132):: box_cmd
      logical :: box_given, mod_seq,ok

      cell_fst=(/cell_a,cell_b,cell_c,90.0_dp,90.0_dp,cell_gamma*to_deg/)
      !Writing file for FullProf Studio
      write(unit=i_fst,fmt="(a)") "! "
      write(unit=i_fst,fmt="(a)") "!File Generated by the program FAULTS"
      write(unit=i_fst,fmt="(a)") "!Title: "//trim(ttl)
      write(unit=i_fst,fmt="(a)") "SPACEG  P 1"
      !Calculate the effective unit cell depending on the number of layers to be diplayed
      !Analysis of the fst_cmd lines
      box_given=.false.
      mod_seq=.false.
      nl=0
      stck_vect=(/0.0,0.0,0.25/)
      do i=1,num_fst
        aux=u_case(fst_cmd(i))
        j=index(aux,"SEQ")
        if( j /= 0) then
          read(unit=fst_cmd(i)(j+3:),fmt=*,iostat=ier) nlayers
          if(nl+nlayers > 1000) exit
          read(unit=fst_cmd(i)(j+3:),fmt=*,iostat=ier) nlayers,seq(nl+1:nl+nlayers)
          if(ier /= 0) then
            write(unit=*,fmt="(a)") " => Error reading the sequence of layers! No FST file will be generated ..."
            return
          end if
          nl=nl+nlayers
        end if
        if(index(aux,"BOX") /= 0) then
          box_cmd=fst_cmd(i)
          box_given=.true.
        end if
        j=index(aux,"STACK_VECT")
        if( j /= 0) then
          read(unit=fst_cmd(i)(j+10:),fmt=*,iostat=ier) stck_vect
          if(ier /= 0) then
            stck_vect=(/0.0,0.0,0.25/)
            write(unit=*,fmt="(a)") " => Error reading the initial stacking vector for FST..."
            write(unit=*,fmt="(a)") "    Taking the default value: [0.0, 0.0, 0.25] "
          end if
        end if
      end do
      nlayers=nl !Total number of layers read in the file
      !Checking if there are impossible sequences and modify the command subsequently
      do
        do i=2,nlayers
          L=seq(i-1) !current number of layer         !(origin layer)
          if(L == 0) cycle
          j=seq(i)                                    !(destination layer)
          !if(crys%l_alpha (j,L) < 0.0001) then !impossible sequences
          if(l_alpha (j,L) < 0.0001) then !impossible sequences
            !write(*,*) "from", L, "to", j
            !write(*,*) "crys%l_alpha ", crys%l_alpha(j,L)   ! this is the value before refinement
            !write(*,*) "l_alpha ", l_alpha(j,L)             ! this is the values after refinement
            seq(i)=0
            mod_seq=.true.
            exit
          end if
        end do
        if(i < nlayers) then
          do j=i+1,nlayers
            seq(j-1)=seq(j)
          end do
          nlayers=nlayers-1
        else if (i==nlayers .and. seq(i)==0) then          !new
          nlayers=nlayers-1                                !new
          exit                                             !new
        else
          exit
        end if
      end do
      !Rewrite the SEQ command
      if(mod_seq) then
        write(unit=*,fmt="(a)")    " => The sequence has been modified by eliminating impossible stacking "// &
                                   "(alpha < 0.0001)! Check your FST_COMMAND SEQ lines!"
        write(unit=*,fmt="(a,i3)") " => The final number of layers is: ",nlayers
        write(unit=*,fmt="(a)")    " => The effective sequence is the following:"
        aux="(a,   i3)"
        write(unit=aux(4:6),fmt="(i3)") nlayers
        write(unit=*,fmt=aux)    "    SEQ: ",seq(1:nlayers)
      end if
      !Calculate the total number of atoms and construct a label for each of them
      nat=0
      do i=1,nlayers
        j=seq(i) !current number of layer
        lact=crys%l_actual(j)
        nat=nat+(crys%centro(lact)+1)*crys%l_n_atoms(lact) !Adding number of atoms in layer j
      end do
      write(unit=*,fmt="(a,3i5)") " => Number of layers and atoms for FP_Studio: ",nlayers,nat
      allocate(xyz(3,nat),atnam(nat))
      call Allocate_Atom_List(nat,Atm,ok) !For CIF file
      ok=.not. ok
      xyz=0.0; atnam=" "
      !Atoms position of the first given layer (for the moment in the given basis: cell_a,cell_b,cell_c)
      j=seq(1)
      lact=crys%l_actual(j)
      k=0
      cvect=stck_vect !Start with a non zero vector to include all the layers within the new unit cell
      do i=1,crys%l_n_atoms(lact)
        k=k+1
        !xyz(:,k) =crys%a_pos(:,i,lact)+cvect
        xyz(:,k) =a_pos(:,i,lact)/pi2+cvect
        !elem=crys%a_name(i,lact)(1:2)
        elem=a_name(i,lact)(1:2)
        if(index("0123456789",elem(2:2)) /= 0) elem(2:2)=" "
        aux=" "
        if(ok) then
          Atm%atom(k)%ChemSymb=elem
          Atm%atom(k)%Occ=a_occup(i,lact)
          Atm%atom(k)%biso =a_b(i,lact)
        end if
        write(unit=aux,fmt="(3(a,i3))") elem,i,"_",j,"_",1
        aux=Pack_String(aux)
        !atnam(k)=trim(aux)//"  "//crys%a_name(i,j)
        atnam(k)=aux//elem
        if(crys%centro(lact) == centro) then
          k=k+1
          !xyz(:,k) = -crys%a_pos(:,i,lact)+cvect
          xyz(:,k) = -a_pos(:,i,lact)/pi2+cvect
          if(ok) then
            Atm%atom(k)%ChemSymb=elem
            Atm%atom(k)%Occ=a_occup(i,lact)
            Atm%atom(k)%biso =a_b(i,lact)
          end if
          m=index(atnam(k-1)," ")
          atnam(k)=atnam(k-1)(1:m-1)//"c"//atnam(k-1)(m+1:)
        end if
      end do

      !Calculate the total displacement vector by adding the stacking vectors
      !and the atoms positions of each layer
      do i=2,nlayers
        L=seq(i-1)
        j=seq(i)
        lact=crys%l_actual(j)
        cvect=cvect+l_r(:,j,L)
        !write(*,"(a,2i4,2(a,3f10.4))") " Layers: ",L,j," Vector: ",l_r(:,j,L), " Sum: ",cvect
        !Atom positions of each layer taking into account the stacking vector
        do n=1,crys%l_n_atoms(lact)
          k=k+1
          !xyz(:,k) = crys%a_pos(:,n,lact) + cvect
          xyz(:,k) = a_pos(:,n,lact)/pi2 + cvect
          !elem=crys%a_name(n,lact)(1:2)
          elem=a_name(n,lact)(1:2)
          if(index("0123456789",elem(2:2)) /= 0) elem(2:2)=" "
          aux=" "
          if(ok) then
            Atm%atom(k)%ChemSymb=elem
            Atm%atom(k)%Occ=a_occup(n,lact)
            Atm%atom(k)%biso =a_b(n,lact)
          end if
          write(unit=aux,fmt="(3(a,i3))") elem,n,"_",j,"_",i
          aux=Pack_String(aux)
          atnam(k)=aux//elem
          if(crys%centro(lact) == centro) then
            k=k+1
            !xyz(:,k) = -crys%a_pos(:,n,lact) + cvect
            xyz(:,k) = -a_pos(:,n,lact)/pi2 + cvect
            if(ok) then
              Atm%atom(k)%ChemSymb=elem
              Atm%atom(k)%Occ=a_occup(n,lact)
              Atm%atom(k)%biso =a_b(n,lact)
            end if
            m=index(atnam(k-1)," ")
            atnam(k)=atnam(k-1)(1:m-1)//"c"//atnam(k-1)(m+1:)
          end if
        end do
      end do
      cvect=cvect+stck_vect !add another stacking vector to include all layers within the new cell
      cell_fst(3)=cvect(3)*cell_c
      xyz(1:2,1:nat)=mod(xyz(1:2,1:nat)+10.0_cp,1.0_cp)
      xyz(3,1:nat)=xyz(3,1:nat)/cvect(3) !Correcting the z-fractional coordinate to the new cell
      write(unit=i_fst,fmt="(a,6f12.5)") "CELL ",cell_fst
      if(box_given) then
        write(unit=i_fst,fmt="(a)") trim(box_cmd)
      else
        write(unit=i_fst,fmt="(a)") "box -1.30 1.30  -1.30 1.30 -0.1 1.1"
      end if
      !Writing atoms
      do i=1,nat
        write(unit=i_fst,fmt="(a,3f12.5)") "Atom "//atnam(i),xyz(:,i)
      end do
      write(unit=i_fst,fmt="(a)") " "
      ! Writing CIF file for VESTA and VESTA file
      if(ok) then
        celax=cell_fst(1:3); celang=cell_fst(4:6)
        Call Set_Crystal_Cell(celax,celang, Celld)
        Call Set_SpaceGroup("P 1",SpGr)
        do i=1,nat
          Atm%atom(i)%x=xyz(:,i)
          m=index(atnam(i)," ")
          Atm%atom(i)%Lab=atnam(i)(1:m-1)
          Atm%atom(i)%Mult=1
        end do
        call Write_Cif_Template(trim(filenam)//"_flts.cif",2,"! File Generated by the program FAULTS: "//trim(ttl),Celld,SpGr,Atm)
        call Write_Vesta_File()
      end if
      !Writing the rest of FST commands
      do i=1,num_fst
        j=index(L_case(fst_cmd(i)),"seq")
        if( j /= 0) cycle
        j=index(L_case(fst_cmd(i)),"box")
        if( j /= 0) cycle
        !Put a command per line
        call Get_Separator_Pos(fst_cmd(i),";",pos,ncar)
        if(ncar == 0) then
          write(unit=i_fst,fmt="(a)") trim(fst_cmd(i))
        else
          k=1
          do L=1,ncar
            write(unit=i_fst,fmt="(a)") fst_cmd(i)(k:pos(L)-1)
            k=pos(L)+1
          end do
          write(unit=i_fst,fmt="(a)") fst_cmd(i)(k:)
        end if
      end do
      ! the unit is closed in the calling program
      return
      contains
        Subroutine Write_Vesta_File()

           integer                                     :: lun, i, j, n,np, pol,k,cent
           character(len=2)                            :: elem
           character(len=80)                           :: box_cmd,aux
           character(len=2), dimension(:), allocatable :: poly
           character(len=80),dimension(:), allocatable :: cmd_bond
           character(len=80)  :: cmd

           open(newunit=lun,file=trim(filenam)//"_flts.vesta",action="write",status="replace")
           write(unit=lun,fmt="(a)") "#VESTA_FORMAT_VERSION 3.1.9"
           write(unit=lun,fmt="(a)") " "
           write(unit=lun,fmt="(a)") " "
           write(unit=lun,fmt="(a)") "CRYSTAL"
           write(unit=lun,fmt="(a)") "TITLE"
           write(unit=lun,fmt="(a)") "  "//trim(ttl)
           write(unit=lun,fmt="(a)") " "
           write(unit=lun,fmt="(a)") "IMPORT_STRUCTURE"
           write(unit=lun,fmt="(a)") trim(filenam)//"_flts.cif"
           write(unit=lun,fmt="(a)") " "
           write(unit=lun,fmt="(a)") "BOUND"
           box_cmd=" "
           do i=1,num_fst
             aux=u_case(fst_cmd(i))
             j=index(aux,"BOX")
             if(j /= 0) then
               box_cmd=fst_cmd(i)(j+3:)
               exit
             end if
           end do
           if(len_trim(box_cmd) == 0) then
              write(unit=lun,fmt="(a)") "    -0.1      1.1      -0.1      1.1      -0.1     1.1 "
           else
              write(unit=lun,fmt="(a)") trim(box_cmd)
           end if
           write(unit=lun,fmt="(a)")"  0   0   0   0  0"

           write(unit=lun,fmt="(a)") "SBOND"

           n=0; np=0
           pol=0
           allocate(cmd_bond(Atm%natoms), poly(Atm%natoms))
           cmd_bond=" "
           poly=" "
           do i=1,num_fst
               aux=u_case(fst_cmd(i))
               j=index(aux,"CONN")
               if(j /= 0) then
                 n=n+1
                 cmd=adjustl(fst_cmd(i)(j+4:))
                 j=index(cmd," ")
                 elem=cmd(1:j-1)
                 cmd=adjustl(cmd(j:))
                 j=index(cmd," ")
                 cmd_bond(n)=elem//"  "//cmd(1:j)//"  "//cmd(j+1:)
                 cycle
               end if
               j=index(aux,"POLY")
               if(j /= 0)  then
                  pol=1
                  np=np+1
                  poly(np)=adjustl(fst_cmd(i)(j+4:))
               end if
           end do
           do i=1,n
             j=index(cmd_bond(i)," ")
             elem=cmd_bond(i)(1:j-1)
             cent=0
             do k=1,np
               if(trim(poly(k)) == trim(elem)) then
                  cent=1
                  exit
               end if
             end do
             write(unit=lun,fmt="(i3,a,5i3)") i,"  "//trim(cmd_bond(i)),0,1,cent,0,1
           end do
           write(unit=lun,fmt="(a)")    "  0 0 0 0"
           write(unit=lun,fmt="(a)")    " "
           write(unit=lun,fmt="(a)")    "STYLE"
           write(unit=lun,fmt="(a)")    "MODEL   2  1  0"
           write(unit=lun,fmt="(a)")    "SURFS   0  1  1"
           write(unit=lun,fmt="(a)")    "SECTS  96  0"
           write(unit=lun,fmt="(a,i1)") "POLYS  ",pol
           call flush(lun)
           close(unit=lun)
        End Subroutine Write_Vesta_File

    End Subroutine Write_FST_VESTA


    Subroutine Write_Prf(diff_pat,i_prf)
       !-----------------------------------------------
       !   D u m m y   A r g u m e n t s
       !-----------------------------------------------
       Type(Diffraction_Pattern_Type), intent(in) :: diff_pat
       integer,                        intent(in) :: i_prf
       !-----------------------------------------------
       !   L o c a l   V a r i a b l e s
       !-----------------------------------------------
       integer ::  i, j,k, iposr, ihkl, irc, nvk,nphase,ideltr,MaxNumRef
       logical :: verif

       real :: twtet, dd, scl,yymi,yyma,t1,t2,sintlm
       character (len=1)   :: tb
       character (len=50)  :: forma1,forma2
       integer, dimension(max_bgr_num+1) :: nref
       integer, dimension(max_bgr_num+1) :: ibgr
       logical :: tick_marks
       !character (len=200) :: cell_sp_string
       !-----------------------------------------------
       !check for very high values of intensities and rescal everything in such a case
       ! scl: scale factor scl=1.0 (normal ymax < 1e6, 0.1 multiplier)
       t1=maxval(diff_pat%y)
       t2=maxval(diff_pat%ycalc)
       yyma=max (t1,t2)
       scl=1.0
       nref=0
       ibgr =0
       nphase=1
       tick_marks=.true.

       !Calculate average Bragg positions in case the average cell and space group has been given

       if(aver_cellcalc) then
         write(unit=*,fmt="(a)") " => Calculating average cell ... "
           call calc_avcell(nlay_avcell,lay_avcell,verif)
           if (.not. verif)then
                write(unit=*,fmt="(a)") " => ERROR calculating the average cell"
                write(unit=*,fmt="(a)") " => No tick marks for pseudo-reflections in PRF file"
                tick_marks=.false.
           else
                call Set_Crystal_Cell(aver_cell, aver_ang, aCell)
                write(unit=*,fmt="(a)") " => The average cell is:"
                write(unit=*,fmt="(a,f12.5,a,f12.5,a,f12.5)") &
                     "    a=    ", aCell%cell(1), "    b=   ",aCell%cell(2),"   c=    ", aCell%cell(3)
                write(unit=*,fmt="(a,f12.5,a,f12.5,a,f12.5)") &
                     "    alpha=", aCell%ang(1), "    beta=",aCell%ang(2), "   gamma=", aCell%ang(3)
           end if
       end if

       !write(unit=*,fmt="(a)") " => Calculating Bragg positions ... "
       n_hkl=0
       if((aver_cellgiven .or. aver_cellcalc) .and. aver_spggiven .and. tick_marks) then
          write(unit=*,fmt="(a)") " => Writting Bragg positions ... "
          sintlm=sind(0.5*(thmax+2.0))/lambda
          MaxNumRef = get_maxnumref(sintlm,aCell%CellVol,mult=spg%Multip)
          if (allocated(reflections)) then
            deallocate(reflections)
            allocate (reflections(MaxNumRef))
          else
            allocate (reflections(MaxNumRef))
          end if

          !> Procedure to calculation all reflections
          !call HKL_GEN(cell,grp_espacial,.true.,val1,val2,num,reflections) !Not ordered
          call HKL_UNI(acell,spg,.true.,0.0,sintlm,"s",n_hkl,reflections) !Ordered
          if(allocated(dos_theta))   deallocate(dos_theta)
          allocate(dos_theta(n_hkl))
          if(allocated(hkl_list))   deallocate(hkl_list)
          allocate(hkl_list(3,n_hkl))
          dos_theta=0.0; hkl_list=0
          do i=1,n_hkl
            dos_theta(i)=2.0*asind(reflections(i)%S*lambda)+crys%zero_shift
            hkl_list(:,i)=reflections(i)%h
          end do
       end if
       nref(1)=n_hkl
       ibgr(1)=1
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
       !Calculate the number of phases taking into account the number of background patterns
       !having an accompanying file with reflections
       if(crys%bgrpatt) then
         do i=1,crys%num_bgrpatt
           if(len_trim(crys%bfilehkl(i)) /= 0) then
            nphase=nphase+1
            nref(nphase)=bgr_hkl_nref(i)
            ibgr(nphase) = i
          end if
         end do
       end if
       write(i_prf,'(A)') trim(diff_pat%title)

       write(i_prf,'(i3,i7,5f12.5,i5)') nphase,diff_pat%npts,lambda,lambda2,0.0,0.0,0.0,0

       nvk=0
       if(nphase <= 8) then
         write(i_prf,'(17i6)')(Nref(i),i=1,nphase), (nvk,i=1,nphase) , nexcrg
       else
         write(i_prf,'(16i6)')(Nref(i),i=1,8), (nvk,i=1,8)
         write(i_prf,'(17i6)')(Nref(i),i=9,nphase), (nvk,i=9,nphase), nexcrg
       end if

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
       ideltr=INT(yyma/16)
       DO i=1,n_hkl
         WRITE(i_prf,'(f12.4,9a,i8,a,3i3,a,2i3)')  &
             dos_theta(i),tb,'        ',tb,'        ',tb,'        ',  &
             tb,'        ',tb,iposr, tb//'(',hkl_list(:,i),')'//tb,ihkl,irc
       END DO
       do i=2,nphase
         iposr=-(i-1)*ideltr
         k=ibgr(i)
         do j=1,nref(i)
           twtet=bgr_hkl_pos(j,k)+crys%zero_shift
           write(i_prf,'(f12.4,9a,i8,a,3i3,a,2i3)')  &
             twtet,tb,'        ',tb,'        ',tb,'        ',  &
             tb,'        ',tb,iposr, tb//'(',bgr_hkl_ind(:,j,k),')'//tb,ihkl,irc
         end do
       end do

       RETURN
    End Subroutine Write_Prf


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
    use CFML_Math_General,          only : spline, splint, sind, cosd, asind
    use diffax_mod
    use read_data,                  only : crys, read_structure_file, length,   opti , cond, vs, &
                                           bgr_patt, calculate_Aberrations
    use diffax_calc,                only : salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, optimz,point,  &
                                           gospec, gostrk, gointr,gosadp, getfnm, nmcoor
    use Dif_compl,                  only : Write_Prf, Faults2diffax, vs2faults, Var_assign

    implicit none

    public  ::  Cost_LMA, apply_aberrations
    type (diffraction_pattern_type), save :: difpat
    integer, parameter, public            :: iout = 25,i_out=20


    contains


!! Subroutine Levenberg_Marquard_Fit(Model_Functn, m, c, Vs, chi2, infout,residuals)      Cost_LMA, Nop, Cond, Vs, chi2, texte
!!--..            Integer,                     Intent(In)      :: m        !Number of observations
!!--..            type(LSQ_conditions_type),   Intent(In Out)  :: c        !Conditions of refinement
!!--..            type(LSQ_State_Vector_type), Intent(In Out)  :: Vs       !State vector
!!--..            Real (Kind=cp),              Intent(out)     :: chi2     !final Chi2
!!--..            character(len=*),            Intent(out)     :: infout   !Information about the refinement (min length 256)
!!--..            Real (Kind=cp), dimension(:),optional, intent(out) :: residuals
!! End Subroutine

!!Interface No_derivatives
!!--..           Subroutine Model_Functn(m, n, x, fvec, iflag)             !Model Function subroutine
!!--..             Use CFML_GlobalDeps, Only: cp
!!--..             Integer,                       Intent(In)    :: m, n    !Number of observations and free parameters
!!--..             Real (Kind=cp),Dimension(:),   Intent(In)    :: x       !Array with the values of free parameters: x(1:n)
!!--..             Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec    !Array of residuals fvec=(y-yc)/sig : fvec(1:m)
!!--..             Integer,                       Intent(In Out):: iflag   !If iflag=1 calculate only fvec without changing fjac
!!--..           End Subroutine Model_Functn                               !If iflag=2 calculate only fjac keeping fvec fixed
!!--..         End Interface No_Fderivatives

    Subroutine Cost_LMA(m,npar,v,fvec,iflag)                 !Levenberg Marquardt
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
      real                     :: rf
      real, save               :: chi2,chiold=1.0e35
      integer, save            :: iter=-1

      shift(1:npar) = v(1:npar) - vector(1:npar)

      Select Case(iflag)

        Case(1)  !Calculation of fvec and updating completely the parameters

           fvec=0.0
           do i = 1, crys%npar
             state(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
           end do
          !crys%list(:) = state(:)
           vector(1:npar) = v(1:npar) !vector upload
           crys%Pv_refi(1:npar) = v(1:npar)
           call Pattern_Calculation(state,ok)
           if(.not. ok) then
             write(*,"(a)") " => Error calculating spectrum, please check the input parameters"
           else
             call calc_fullpat_lma(difpat, fvec, chi2,rf)
           end if
           numcal = numcal + 1  !Counter for optimz, detun, etc
           if(chi2 < chiold) then
              chiold=chi2
              iter = iter + 1
              chiold=chi2
              write(i_out,"(a)")" ------------------------------------------------------------------------------------------------"
              write(i_out,"(a)")"              Parameter Name           New_Value     Old_Value    Multiplier         Shift ParNum"
              write(i_out,"(a)")" ------------------------------------------------------------------------------------------------"
              do i = 1, crys%npar
                 write(i_Out,"(a,i3,a,4f14.5,i5)") " State(",i,"):  "//namepar(i),state(i), crys%list(i), mult(i), &
                                                  shift(crys%p(i)),crys%p(i)
                 if(Cond%nprint < 0) &
                 write(*,"(a,i3,a,4f14.5,i5)") " State(",i,"):  "//namepar(i),state(i), crys%list(i), mult(i), &
                                                  shift(crys%p(i)),crys%p(i)
              end do
              write(*,"(a,i4,2(a,f12.5))")  " => Iteration ",iter,"   R-Factor = ",rf,"   Chi2 = ",chi2
              write(i_out,"(a,i4,2(a,f12.5))")  " => Iteration ",iter,"   R-Factor = ",rf,"   Chi2 = ",chi2
           end if
           crys%list(:) = state(:)
           !svdfvec(1:m)=fvec(1:m)

        Case(2)  !Calculation of numerical derivatives

           do i = 1, crys%npar
             state(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
           end do
           call Pattern_Calculation(state,ok)
           if(.not. ok) then
             write(*,"(a)") " => Error calculating spectrum, please check the input parameters"
           else
             call calc_fullpat_lma(difpat, fvec, chi2, rf)
           end if
           !fvec(1:m)=svdfvec(1:m)

        Case(0)  !Printing
           write(*,"(a,i4,a,f14.4)")  " => Iteration ",iter,"   Chi2 = ",chi2
           do i=1, crys%npar
             write(*,"(a,f14.5)")  "  ->  "//namepar(i), state(i)
           end do

      End Select

      ok = .true.

      IF(cfile) CLOSE(UNIT = cntrl)
      return

    End subroutine Cost_LMA

    Subroutine Bck_Chebychev(bk)
      real, dimension(:), intent(out) :: bk
      real    :: x,c,thx,rj
      integer :: i,j

      bk=0.0
      do i=1,difpat%npts
        x=difpat%x(i)
        bk(i)=crys%chebp(1)
        thx=2.0*(x-0.5*(difpat%x(difpat%npts)+difpat%x(1)))/(difpat%x(difpat%npts)-difpat%x(i))
        if(thx < -1.0) thx=-1.0
        if(thx >  1.0) thx= 1.0
        c=acos(thx)
        do j=2,crys%cheb_nump
          rj=real(j)
          bk(i)=bk(i)+crys%chebp(j)*cos(rj*c)
        end do
      end do
    End Subroutine Bck_Chebychev


    Subroutine calc_fullpat_lma(pat, fvec, chi2,r)
      type (diffraction_pattern_type), intent(in out):: pat
      Real (Kind=cp),Dimension(:),     Intent(in Out):: fvec
      real,                            Intent(   out):: chi2, r
      ! Local variables
      real                                           :: a,b,c,delta
      integer                                        :: punts
      integer                                        :: i,j
      real, dimension(pat%npts)                      :: bk

      pat%scal = crys%patscal
      Do j = 1, pat%npts
       pat%ycalc(j)  = brd_spc(j)
       pat%ycalc(j)  = pat%scal * pat%ycalc(j)+ pat%bgr(j)
      End do

      !Adding Chebychev polynomial background
      if(crys%bgrcheb) then
        call Bck_Chebychev(bk)
        pat%ycalc=pat%ycalc+bk
      end if

      if(crys%num_bgrpatt > 0) then !Adding contributions of background patterns
        do i=1,crys%num_bgrpatt
          pat%ycalc=pat%ycalc+crys%bscalpat(i)*bgr_patt(:,i)
        end do
      end if

      a=0.0; b=0.0; c=0.0
      punts=0
      do_a: do j=1, pat%npts
        do i=1,nexcrg
         if(pat%x(j) >= alow(i) .and. pat%x(j) <= ahigh(i)) cycle do_a
        end do
        punts=punts+1
        a= a + pat%y(j)
        delta= pat%y(j) - pat%ycalc(j)
        b= b + abs(delta)
        fvec(j)=  delta/sqrt(pat%sigma(j))
        c = c + fvec(j)*fvec(j)
      end do do_a
      r =  b/a *100.0
      !chi2= c/(punts-opti%npar)
      chi2= c/(pat%npts-opti%npar)  !Use all points => this makes Chi2 smaller but consistent with
      return                        !the calculation in CrysFML
    End subroutine calc_fullpat_lma

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
      if(calculate_aberrations) call Apply_Aberrations()
      return

    End Subroutine Pattern_Calculation

    !Subroutine Apply_Aberrations(pat)
    Subroutine Apply_Aberrations()
    !type (diffraction_pattern_type), intent(in out):: pat
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
        call splint(true_2th,broad_spect,der2v,n_high,tt,ycal)
        brd_spc(i)=ycal
      end do
      return
    End Subroutine Apply_Aberrations

   End module dif_ref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   PROGRAM FAULTS

     use CFML_GlobalDeps,              only : sp , cp, OPS_SEP
     use CFML_String_Utilities,        only : number_lines , reading_lines ,  init_findfmt, findfmt ,iErr_fmt, &
                                              getword, err_string, err_string_mess, getnum, Ucase, lcase
     use CFML_Diffraction_patterns,    only : read_pattern , diffraction_pattern_type , err_diffpatt, err_diffpatt_mess,  &
                                              read_background_file
     use CFML_Optimization_General,    only : Opt_Conditions_Type, Local_Optimize
     use CFML_Crystal_Metrics,         only : Set_Crystal_Cell, Crystal_Cell_Type
     use CFML_Optimization_LSQ,        only : Levenberg_Marquardt_Fit,Info_LSQ_LM
     use CFML_LSQ_TypeDef,             only : LSQ_Conditions_type, LSQ_State_Vector_Type
     use CFML_Random_Generators,       only : random_poisson
     use diffax_mod
     use read_data,                    only : filenam, read_structure_file, length, bgr_patt, Read_Bgr_patterns,  &
                                              crys, opti, cond, Vs, Err_crys, Err_crys_mess,fst_given
     use diffax_calc ,                 only : salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, optimz,point,  &
                                              gospec, gostrk, gointr,gosadp, chk_sym, get_sym, overlp, nmcoor , getfnm
     use Dif_compl,                    only : Write_Prf, write_ftls, Faults2diffax, vs2faults, Var_assign, Write_FST_VESTA
     use dif_ref

     implicit none


      real                    :: rpl,  theta , ymax, ymini , ymin ,deg, tini,tfin
      LOGICAL                 :: ok, ending , gol, p_ok, arggiven=.false.
      INTEGER                 :: i ,n ,j, l , ier , fn_menu,a,b,c ,aa,bb,cc, e, narg
      character(len=100)      :: pfile, bfile , bmode
      character(len=100)      :: pmode
      character(len=256)      :: path_name
      character(len=10)       :: time, date
      character(len=1)        :: keyw
      Real (Kind=cp)          :: chi2     !final Chi2
      character(len=3000)     :: infout   !Information about the refinement (min length 256)

      pi = four * ATAN(one)
      pi2 = two * pi
      deg2rad = pi / one_eighty
      rad2deg = one_eighty / pi

      ending = .false.
      ok = .true.

      call get_environment_variable('FULLPROF',path_name)
      i=len_trim(Path_name)
      if( i == 0) then
         sfname = 'data.sfc'
      else
        if(path_name(i:i) == OPS_SEP) then
           sfname=trim(path_name)//'data.sfc'
        else
           sfname=trim(path_name)//OPS_SEP//'data.sfc'
        end if
      end if

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
      call cpu_time(tini)
      !WRITE(op,fmt=*) "=> Looking for scattering factor data file '",  sfname(:),"'"
      OPEN(UNIT = sf, FILE = sfname, status="old",action="read",position="rewind")
      !WRITE(op,fmt=*) "=> Opening scattering factor data file '",  sfname(:),"'"

      filenam = trim(infile(1:(index(infile,'.')-1)))
      if (index(infile,'_new')>0) then
      filenam = trim(filenam(1:(index(infile,'_new')-1)))
      end if

     if(replace_files) then
       Call getfnm(filenam,outfile, '.out', ok,replace_files)
     else
       Call getfnm(filenam,outfile, '.out', ok)            !check for any existing .out file, if any choose another filename not to overwrite the existing one.
     end if
     open(unit=i_out, file=trim(outfile),status="replace",action="write")
    ! open(unit=i_out, file=trim(filenam)//".out",status="replace",action="write")
      CALL salute(i_out)

      call read_structure_file(infile, gol)

      if (err_crys) then
        write(unit=*,fmt="(a)") " => ERROR in "//trim(infile)//": "//trim(err_crys_mess)
        call Close_Faults()
      else
        write(op, fmt="(a)") " => Structure input file read in"
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
         if (crys%bgrinter) then
           write(*,"(a)") " => Reading Background file="//trim(background_file)
           call read_background_file(background_file, mode ,difpat)
           if(Err_diffpatt) print*, trim(err_diffpatt_mess)
           if(crys%num_bgrpatt > 0) then
             allocate(bgr_patt(difpat%npts,crys%num_bgrpatt))
             call Read_Bgr_patterns()
           end if
         end if
      end if

      !Operations before starting :
      check_sym = .false.
      IF(symgrpno == UNKNOWN) THEN
        symgrpno = get_sym(ok)
        IF(.NOT. ok) GO TO 999
        WRITE(op,"(a)") ' => Diffraction point symmetry is '//pnt_grp
        IF(symgrpno /= 1) THEN
          WRITE(op,"(a,i6)") '    to within a tolerance of one part in ',  &
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

            WRITE(op,"(a)") ' => Start simulation'
          ! What type of intensity output does the user want?


                IF(funct_num == 3) THEN
                   write(unit=*,fmt="(a)") " => Calculating powder diffraction pattern"
                   CALL gospec(infile,outfile,ok)

                     Do j = 1, n_high
                         ycalcdef(j) = crys%patscal*brd_spc(j)+crys%bckg_level
                         call random_poisson(ycalcdef(j),i)
                         ycalcdef(j)=real(i)
                     end do

                     if(replace_files) then
                       Call getfnm(filenam,outfile, '.dat', ok,replace_files)
                     else
                       CALL getfnm(filenam, outfile, '.dat', ok)
                     end if

                     OPEN(UNIT = iout, FILE = outfile, STATUS = 'replace')
                        write(unit = iout,fmt = *)'!', outfile
                        write(unit = iout,fmt = '(3f12.4)')thmin, step_2th,thmax
                     !  theta = thmin +(j-1)*d_theta
                        write(unit = iout,fmt = '(8f12.2)') ( ycalcdef(j), j=1, n_high    )

                    CLOSE(UNIT = iout)
                    ok = .true.
                ELSE IF(funct_num == 4) THEN
                   CALL gosadp(infile,outfile,ok)
                ELSE
                   WRITE(op,"(a)") ' => Unknown function type.'
                END IF




        !     IF(ok .AND. n /= 3) THEN
        !     96   WRITE(op,"(a)") ' => Enter 1 to return to function menu.'
        !          READ(cntrl,*,ERR=96,END=999) fn_menu
        !          IF(fn_menu == 1) GO TO 10
        !     END IF

!
          Case (4) !LMA

            WRITE(op,"(a)") ' => Start LMA refinement'

            chi2o = 1.0E10                         !initialization of agreement factor
          !  rpl = 0
            do i=1, opti%npar                       !creation of the step sizes
              vector(i) = crys%Pv_refi(i)
            end do

            call Levenberg_Marquardt_Fit(cost_LMA, difpat%npts, cond, Vs, chi2, infout)

            !Output the final list of refined parameters
            Call Info_LSQ_LM(Chi2,op,Cond,Vs)
            Call Info_LSQ_LM(Chi2,i_out,Cond,Vs)

            write(*,"(a)") " => "// trim(infout)


       call vs2faults(vs, crys)
      !call Restore_Codes() !the ref_* variables are now the original codes


            if(replace_files) then
              Call getfnm(filenam,outfile, '.prf', ok,replace_files)
            else
              Call getfnm(filenam, outfile, '.prf', ok)
            end if
             if (ok) then
               OPEN(UNIT = iout, FILE = outfile, STATUS = 'replace')
               call Write_Prf(difpat,iout)
             else
               write(*,"(a)") ' => The output file .prf cannot be created'
             end if
             !CALL getfnm(trim(filenam)//"_new", outfile, '.ftls', ok)
             if (ok) then
               filenam=trim(filenam)
               i=index(filenam,"_new")
               if(i /= 0) filenam=filenam(1:i-1)
               OPEN(UNIT = i_flts, FILE = trim(filenam)//"_new.flts", STATUS = 'replace',action="write")
               call Write_ftls(crys,i_flts, vs)
             else
               write(*,"(a)") ' => The output file .flts cannot be created'
             end if


!-------------------------------------------------------------------------------------------------------------------------------

            Case default

                write(*,"(a)") " => Problems reading mode "


          End select

      if(fst_given) then  !Generate file for FullProf Studio
          if(replace_files) then
              Call getfnm(filenam,outfile, '.fst', ok,replace_files)
          else
              Call getfnm(filenam,outfile, '.fst', ok)
          end if
        OPEN(UNIT = i_fst, FILE = trim(outfile), STATUS = 'replace',action="write")
        call Write_FST_VESTA()
        close(unit=i_fst)
      end if
      ending = .true.

      999 IF(cfile) CLOSE(UNIT = cntrl)
      IF(ok .AND. ending) THEN
        WRITE(op,"(a)") ' => FAULTS ended normally....'
      ELSE
        WRITE(op,"(a)") ' => FAULTS was terminated abnormally!'
      END IF
      call cpu_time(tfin)
      tini=tfin-tini
      tini=tini/60.0
      tfin=int(tini)
      tini=(tini-tfin)*60.0
      write(op,"(a,i4,a,f8.4,a)") " => Total CPU-time: ",int(tfin)," minutes and ",tini," seconds"
      write(i_out,"(a,i4,a,f8.4,a)") " => Total CPU-time: ",int(tfin)," minutes and ",tini," seconds"
      call Close_Faults()

   END PROGRAM FAULTS

