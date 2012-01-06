
  Module Diff_ref

        use diffax_mod
        use CFML_GlobalDeps,            only : sp
        use CFML_Diffraction_Patterns , only : diffraction_pattern_type

        implicit none

        private

        !public subroutines
        public :: scale_factor

        contains
!________________________________________________________________________________________________________________________

       Subroutine scale_factor (pat,r)
         type (diffraction_pattern_type)  , intent (in out) :: pat
         real                             , intent (   out) :: r
         real                                               :: a ,  c
         real                                               :: thmin, thmax  , thmin_o
         integer                                            :: n_high   , punts=0
         integer                                            ::  j

         n_high = INT(half*(th2_max-th2_min)/d_theta) + 1    !rad

         a=0
         Do j = 1, n_high
            pat%ycalc(j)  = brd_spc(j)
            a = a + pat%ycalc(j)
         End do

         thmin=th2_min * rad2deg        !deg
         thmax=th2_max * rad2deg
         thmin_o = pat%xmin * deg2rad  !rad

         if( th2_min /= thmin_o ) then
           punts= INT(half*(th2_min-thmin_o)/d_theta) ! to calculate Rp and chi2 without the excluded regions
         end if

         c=0
         Do  j = 1, pat%npts
           if (pat%x(j) < thmin) cycle
           if (pat%x(j) > thmax) exit
           c = c + (pat%y(j) - pat%bgr(j))
         End do

         pat%scal = c/a

         Do j = 1, n_high
          pat%ycalc(j) = pat%scal * pat%ycalc(j)+ pat%bgr(j)
         End do

         call calc_par(pat, n_high,punts, r)
         return
      End subroutine scale_factor

!____________________________________________________________________________________________________________________________

    Subroutine calc_par (pat, n_high, punts, r)
      type (diffraction_pattern_type), intent(in    ) :: pat
      integer                        , intent(in    ) :: n_high, punts
      real                           , intent(   out) :: r
      real                                            :: chi2
      real                                            :: a,b,c
      integer                                         :: j

      a=0.0
      b=0.0
      do j=1, n_high
        a= a + pat%y(j+punts)
        b= b + (abs(pat%y(j+punts) - pat%ycalc(j)))
      end do
        r =  b/a *100
        c=0.0
      do j=1, n_high
        c = c  + ((1/(pat%sigma(j+punts)**2))*(pat%y(j+punts)-pat%ycalc(j))**2 )
      end do
        chi2= c/n_high *100
      write (*,*) r, chi2   , n_high
      return
    End subroutine calc_par

  End module Diff_ref
!________________________________________________________________________________________________________________
   Module dif_ref
     use CFML_GlobalDeps,            only : sp , cp
     use CFML_String_Utilities,      only : number_lines , reading_lines ,  init_findfmt, findfmt ,iErr_fmt, getword, &
                                            err_string, err_string_mess, getnum, Ucase
     use CFML_Simulated_Annealing
     use CFML_Crystal_Metrics,       only : Set_Crystal_Cell, Crystal_Cell_Type
     use CFML_Diffraction_patterns , only : diffraction_pattern_type
     use diffax_mod
     use read_data,                  only : crys_2d_type, read_structure_file, length, rdrnge, choice
     use diffax_calc,                only : salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, optimz,point,  &
                                            gospec, gostrk, gointr,gosadp, getfnm, nmcoor
     use Diff_ref,                   only : scale_factor


    implicit none

    public  :: Cost3 , F_cost

    type (State_Vector_Type), public  :: st
    type (diffraction_pattern_type)   :: difpat

    contains


    Subroutine  Cost3(vref,rp)      !Simulated annealing

          real,   dimension (:), intent (in    ) :: vref
          real,                  intent (   out) :: rp
          logical                                :: ok
          integer                                :: n , j ,i,k  , label , a, b
          real, dimension(80)                    :: shift, state, menor=1 , multi,tar
          character                              :: cons


          do i= 1, npar                          !shift calculation

                shift(i) = vref(i) - gen(i)

          end do

        !*******RESTRICTIONS*******

         do i = 1, numpar

             if (index (namepar(i) , 'alpha' ) == 1 )  then       !To avoid negative values of alpha
                ! read (unit = namepar(i)(6:7), fmt = "(2i1)" ) b,a
                ! if  (l_alpha(a,b) .le. menor(b) ) then
                !        menor(b) = l_alpha(a, b)
                !        multi(b) = mult(i)
                ! end if
                 state(i) = vector(i) +  mult(i) * shift(pnum(i))
                 if (state(i) < zero ) then
                        write(*,*) 'Attention, shift was higher than alpha:  new shift applied'
                        shift(pnum(i)) = shift(pnum(i))/2
                        if (state(i) < zero ) then
                        shift(pnum(i)) = -shift(pnum(i))
                        end if
                 end if

             end if
          end do
        !!!!!!!!!!!!!!!!!!!

          do i=1, numpar   !assignment of new values and more restrictions

            state(i) = vector(i) +  mult(i) * shift(pnum(i))

            if (index (namepar(i) , 'Biso' ) == 1 .and. state(i) < 0 )   state(i) = (-1.0) * state(i)  !Biso only >0
            if (index (namepar(i) , 'v' ) == 1 .and. state(i) < 0 )      state(i) = (-1.0) * state(i)  !v only <0
            if (index (namepar(i) , 'Dg') == 1 .or. index (namepar(i) , 'Dl') == 1 .or.  index (namepar(i) , 'u') == 1 .or.  &
                index (namepar(i) , 'w') == 1 .or. index (namepar(i) , 'x') == 1   ) then  !only >0

              if ( state(i) < 0)  state(i) =  -state(i)

            end if

          End do

         vector(:) = state(:)

         !!!!!!!!!!!!!!!!!!!!!!!!!!
          tar(:)=0                                 !generation of gen
           do i=1, numpar
             if (tar(pnum(i)) == 0 ) gen(pnum(i)) = state(i)
             tar(pnum(i)) = 1
           end do


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !////////CALCULATED PATTERN VARIABLES ASSIGNMENT////////////////////////////////////

          do i=1, numpar

            write(*,*) namepar(i) ,  state(i)
            if (namepar(i) ==  'u')    pv_u  = state(i)
            if (namepar(i) ==  'v')    pv_v  = state(i)
            if (namepar(i) ==  'w')    pv_w  = state(i)
            if (namepar(i) ==  'x')    pv_x  = state(i)
            if (namepar(i) ==  'Dg')   pv_dg = state(i)
            if (namepar(i) ==  'Dl')   pv_dl = state(i)
            if (namepar(i) ==  'cell_a')    cell_a = state(i)
            if (namepar(i) ==  'cell_b')    cell_b = state(i)
            if (namepar(i) ==  'cell_c')    cell_c = state(i)
            if (namepar(i) ==  'num_layers')  l_cnt= state(i)
            if (index( namepar(i) ,'cell_gamma') == 1)   cell_gamma = state(i)

            do j=1, n_layers
              do k=1, n_atoms
                if (index (namepar(i) , 'pos_x' )== 1)     then
                    read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
                    a_pos(1,a,b)  = state(i) * pi2
                end if
                if (index (namepar(i) ,'pos_y' )== 1)    then
                    read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
                    a_pos(2,a,b)  = state(i) * pi2
                end if
                if (index (namepar(i), 'pos_z' ) == 1 )   then
                    read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
                    a_pos(3,a,b)  = state(i) * pi2

                end if
                if (index (namepar(i) ,'Biso')== 1) then
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
                    read (unit = namepar(i)(3:4), fmt = "(2i1)" )b,a
                    l_r(2,a,b)  = state(i)
                end if
                if (index (namepar(i), 'tz' ) == 1)     then
                    read (unit = namepar(i)(3:4), fmt = "(2i1)" ) b,a
                    l_r(3,a,b)  = state(i)
                end if
              end do
            end do
          end do

          !//////////////////////////////////////////////////////////////////////////////////////

          ok = .true.

          if ((conv_d == 1 .or. numcal== 0) .and. ok) ok = get_g()
          if ((conv_d == 1  .or.  numcal== 0 .or. conv_e==1) .and. (ok .and. rndm) ) ok = getlay()
          if ((conv_c == 1 .or. numcal== 0) .and. ok ) CALL sphcst()
          if (numcal == 0 .and. ok) CALL detun()
          if ((conv_b == 1 .or. conv_c == 1 .or. conv_d == 1  .or. numcal == 0 .or. &
              conv_e == 1 .or. conv_f == 1) .and. ok) CALL optimz(infile, ok)
          IF(.NOT.ok) then
           IF(cfile) CLOSE(UNIT = cntrl)
           return
          END IF
          CALL gospec(infile,outfile,ok)
          call scale_factor (difpat,rp)

          numcal = numcal + 1
          write(*,*) ' => Calculated Rp    :   ' , rp
          write(*,*) ' => Best Rp up to now:   ' , rpo

          if (rp .LT. rpo ) then                  !To keep calculated intensity for the best value of rp
            rpo = rp
            statok(1:numpar) = state( 1:numpar)
            write(*,*) 'writing calculated best pattern up to now. Rp :       ', rpo
            do j = 1, n_high
              ycalcdef(j) = difpat%ycalc(j)
            end do
          end if

          ok = .true.
          IF(cfile) CLOSE(UNIT = cntrl)
          return
    End subroutine Cost3

  ! Subroutine Nelder_Mead_Simplex(Model_Functn, Nop, P, Step, Var, Func, C, Ipr)
  !    !---- Arguments ----!
  !    integer,                      intent(in)      :: nop
  !    real(kind=cp), dimension(:),  intent(in out)  :: p, step
  !    real(kind=cp), dimension(:),  intent(out)     :: var
  !    real(kind=cp),                intent(out)     :: func
  !    type(opt_conditions_Type),    intent(in out)  :: c
  !    integer, optional,            intent(in)      :: Ipr
  !
  !    Interface
  !       Subroutine Model_Functn(n,x,f,g)
  !          use CFML_GlobalDeps,  only: cp
  !          integer,                             intent(in) :: n
  !          real(kind=cp),dimension(:),          intent(in) :: x
  !          real(kind=cp),                       intent(out):: f
  !          real(kind=cp),dimension(:),optional, intent(out):: g
  !       End Subroutine Model_Functn
  !    End Interface




    Subroutine  F_cost(n_plex,v,rplex,g)      !SIMPLEX
    use CFML_GlobalDeps,  only: cp
    integer,                        intent (in    ) :: n_plex
    real(kind=cp),  dimension(:),   intent (in    ) :: v
    real(kind=cp),                  intent (   out) :: rplex
    real(kind=cp),dimension(:),optional, intent(out):: g

    logical                 :: ok
    integer                 :: n , j ,i,k  , label , a, b
    real, dimension(80)     :: shift, state, menor=1, multi, tar
    character               :: cons

          do i= 1, n_plex
                shift(i) = v(i) - gen(i)
          end do

           !*******RESTRICTIONS*******

          do i = 1, numpar
             if (index (namepar(i) , 'alpha' ) == 1 )  then                     !To avoid negative values of alpha
                ! read (unit = namepar(i)(6:7), fmt = "(2i1)" ) b,a
                ! if  (l_alpha(a,b) .le. menor(b) ) then
                !        menor(b) = l_alpha(a, b)
                !        multi(b) = mult(i)
                ! end if
                 state(i) = vector(i) +  mult(i) * shift(pnum(i))
                 if(state(i) < zero .or. state(i) > 1) then
                    write(*,*) 'Attention, shift was higher than alpha:  new shift applied'
                    shift(pnum(i)) = shift(pnum(i))/2
                    state(i) = vector(i) +  mult(i) * shift(pnum(i))
                    if (state(i) < zero .or. state(i) > zero) then
                           shift(pnum(i)) = - shift (pnum(i))
                    end if
                 end if

             end if
          end do
        !!!!!!!!!!!!!!!!!!!
          do i=1, numpar                                                            !assignment of new values and more restrictions

            state(i) = vector(i) +  mult(i) * shift(pnum(i))
           ! write(*,*) state(i), vector(i), mult(i), shift(pnum(i))
            if (index (namepar(i),'Biso' ) == 1 .and. state(i) .lt. 0 )   state(i) = (-1.0) * state(i)  !Biso only >0
            if (index (namepar(i),'v' ) == 1 .and. state(i) .gt. 0 )   state(i) = (-1.0) * state(i)  !v only <0
            if (index (namepar(i),'Dg') == 1 .or. index (namepar(i), 'Dl') == 1 .or. &
                index (namepar(i),'u')  == 1 .or. index (namepar(i), 'w') == 1  .or. &
                index (namepar(i), 'x') == 1 ) then  !only >0
                if ( state(i) < 0.0)  state(i) = (-1.0) * state(i)
            end if
          End do


          vector(:) = state(:)

          do i=1, numpar
            write(*,*) namepar(i), state(i)
          end do


          tar(:)=0                                 !generation of gen
          do i=1, numpar
                if (tar(pnum(i)) == 0 ) then
                  gen(pnum(i)) = state(i)
                 ! write(*,*) 'gen', gen(pnum(i)), state(i)
                end if
               tar(pnum(i)) = 1
           end do
          ! write(*,*) 'prova', gen(1:n_plex)

          !////////CALCULATED PATTERN VARIABLES ASSIGNMENT////////////////////////////////////

          do i=1, numpar

            if (namepar(i) ==  'u')    pv_u  = state(i)
            if (namepar(i) ==  'v')    pv_v  = state(i)
            if (namepar(i) ==  'w')    pv_w  = state(i)
            if (namepar(i) ==  'x')    pv_x  = state(i)
            if (namepar(i) ==  'Dg')   pv_dg = state(i)
            if (namepar(i) ==  'Dl')   pv_dl = state(i)
            if (namepar(i) ==  'cell_a')    cell_a = state(i)
            if (namepar(i) ==  'cell_b')    cell_b = state(i)
            if (namepar(i) ==  'cell_c')    cell_c = state(i)
            if (namepar(i) ==  'num_layers')  l_cnt = state(i)
            if (index( namepar(i) ,'cell_gamma') == 1)   cell_gamma = state(i)

            do j=1, n_layers
              do k=1, n_atoms
                if (index (namepar(i) , 'pos_x' )== 1)     then
                    read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
                    a_pos(1,a,b)  = state(i) * pi2
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


          !//////////////////////////////////////////////////////////////////////////////////////

          ok = .true.
          if ((conv_d == 1 .or. numcal== 0) .and. ok ) ok = get_g()
          if ((conv_d == 1  .or.  numcal== 0 .or. conv_e==1) .and. (ok .AND. rndm)) ok = getlay()
          if ((conv_c == 1 .or. numcal== 0) .and. ok ) CALL sphcst()
          if ( numcal == 0 .and. ok ) CALL detun()
          if ((conv_b == 1 .or. conv_c == 1 .or. conv_d == 1 .or. conv_e==1   .or. &
              numcal == 0 .or. conv_f==1)  .and. ok ) CALL optimz(infile, ok)

          IF(.NOT. ok) then
            IF(cfile) CLOSE(UNIT = cntrl)
            return
          END IF
          CALL gospec(infile,outfile,ok)
          call scale_factor(difpat,rplex)
          numcal = numcal + 1
          write(*,*) ' => Calculated Rp    :   ' , rplex
          write(*,*) ' => Best Rp up to now:   ' , rpo
      !    write(*,*) var_plex(:)

          if (rplex < rpo ) then                  !To keep calculated intensity for the best value of rplex
            rpo = rplex
            statok(1:numpar) = state( 1:numpar)

            write(*,*)  ' => Writing the best calculated pattern up to now. Rp : ', rpo
            do j = 1, n_high
              ycalcdef(j) = difpat%ycalc(j)
            end do
            do j=1, l_cnt
              l_seqdef(j) = l_seq(j)
            end do
          end if
          ok = .true.
          IF(cfile) CLOSE(UNIT = cntrl)
          return
    End subroutine F_cost

   End module dif_ref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   PROGRAM FAULTS

     use CFML_GlobalDeps,              only : sp , cp
     use CFML_String_Utilities,        only : number_lines , reading_lines ,  init_findfmt, findfmt ,iErr_fmt, &
                                              getword, err_string, err_string_mess, getnum, Ucase
     use CFML_Diffraction_patterns,    only : read_pattern , diffraction_pattern_type , err_diffpatt, err_diffpatt_mess,  &
                                              read_background_file
     use CFML_Simulated_Annealing
     use CFML_Optimization_General,    only : Nelder_Mead_Simplex,  Opt_Conditions_Type
     use CFML_Crystal_Metrics,         only : Set_Crystal_Cell, Crystal_Cell_Type
     use diffax_mod
     use read_data,                    only : crys_2d_type, new_getfil, read_structure_file, length , rdrnge, choice
     use diffax_calc ,                 only : salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, optimz,point,  &
                                              gospec, gostrk, gointr,gosadp, chk_sym, get_sym, overlp, nmcoor , getfnm
     use Diff_ref ,                    only : scale_factor
     use dif_ref,                      only : cost3 , difpat , F_cost, st

     implicit none

      type (crys_2d_type)                                     :: crys
      type (Opt_Conditions_Type)                              :: opt_c
      type (SimAnn_Conditions_type )                          :: san
     !type(State_Vector_Type)                                 :: sv

      real                                                    :: rpl,  theta ,thmin, thmax , thmin_o  , ymax, ymini , ymin ,deg
      LOGICAL                                                 :: ok, ending , gol
      INTEGER*4                                               :: i ,n ,j, l , ier , fn_menu,a,b,c ,aa,bb,cc, p_resta
      character(len=100)                                      :: pfile, bfile , bmode
      character(len=100)                                      :: pmode, filenam
      integer, parameter                                      :: out = 25
      real, dimension (:), allocatable                        :: resta
      character(len=10)                                       :: time, date

      pi = four * ATAN(one)
      pi2 = two * pi
      deg2rad = pi / one_eighty
      rad2deg = one_eighty / pi

      ending = .false.
      ok = .true.

      CALL salute()


      sfname = 'data.sfc'
      cntrl=ip

      call new_getfil(infile, crys, opt_c, san,st,   gol)                           ! parametros para calculo

      san%Cost_function_name="R-factor"

      IF(gol) then
        ok = sfc()
      else
        GO TO 999
      end if

      do i = 1, numpar              ! to avoid repetitions
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

       Select  case (opt)

          Case (0)
              !Algunas operaciones antes de empezar:
              check_sym = .false.
              IF(symgrpno == UNKNOWN) THEN
                symgrpno = get_sym(ok)
                IF(.NOT.ok) GO TO 999
                WRITE(op,200) 'Diffraction point symmetry is ',pnt_grp
                IF(symgrpno /= 1) THEN
                  WRITE(op,201) '  to within a tolerance of one part in ',  &
                      nint(one / tolerance)
                END IF
              ELSE
                check_sym = .true.
                CALL chk_sym(ok)
                IF(.NOT.ok) GO TO 999
              END IF

              n_high = INT(half*(th2_max-th2_min)/d_theta) + 1

              call overlp()
              call nmcoor ()

              if(allocated(difpat%ycalc) ) deallocate(difpat%ycalc)
              allocate(difpat%ycalc(n_high))

              IF(ok) ok = get_g()
              IF(ok .AND. rndm) ok = getlay()
              IF(ok) CALL sphcst()

              IF(ok) CALL detun()
              IF(.NOT.ok) GO TO 999

! See if there are any optimizations we can do
              IF(ok) CALL optimz(infile, ok)

! What type of intensity output does the user want?
           10  IF(ok) THEN
              ! WRITE(op,100) 'Enter function number: '
              ! WRITE(op,100) '0 POINT, 1 STREAK, 2 INTEGRATE, 3 POWDER PATTERN, 4 SADP : '
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

                   CALL gospec(infile,outfile,ok)

                       n_high = INT(half*(th2_max-th2_min)/d_theta) + 1

                     Do j = 1, n_high
                         ycalcdef(j) = brd_spc(j)
                     end do

                     thmax=th2_max * rad2deg
                     thmin=th2_min * rad2deg
                     d_theta = rad2deg * 2 *d_theta

                     filenam = trim(infile(1:(index(infile,'.')-1)))
                   !  write(*,*) filenam
                     CALL getfnm(filenam, outfile, '.dat', ok)
                !        write(*,*) outfile

                     OPEN(UNIT = out, FILE = outfile, STATUS = 'new')
                        write(unit = out,fmt = *)'!', outfile
                        write(unit = out,fmt = '(3f12.2)')thmin, d_theta,thmax
                     !  theta = thmin +(j-1)*d_theta
                        write(unit = out,fmt = '(8f12.2)') ( ycalcdef(j), j=1, n_high    )

                    CLOSE(UNIT = out)
                    ok = .true.
                ELSE IF(n == 4) THEN
                   CALL gosadp(infile,outfile,ok)
                ELSE
                   WRITE(op,100) 'Unknown function type.'
                END IF

              END IF


              IF(ok .AND. n /= 3) THEN
              96   WRITE(op,100) 'Enter 1 to return to function menu.'
                   READ(cntrl,*,ERR=96,END=999) fn_menu
                   IF(fn_menu == 1) GO TO 10
              END IF

          Case (1)     !SIMULATED ANNEALING

              !Lectura del  pattern experimental y de los scattering factors:


              call  read_pattern (dfile,difpat,fmode )
                 if(Err_diffpatt) then
                    print*, trim(err_diffpatt_mess)
                 else
                    if (th2_min == 0 .and.  th2_max == 0  .and. d_theta == 0) then
                       th2_min =  difpat%xmin    ! if not specified in input file, take the values from the observed pattern
                       th2_max =  difpat%xmax
                       d_theta =  difpat%step
                       th2_min = th2_min * deg2rad
                       th2_max = th2_max * deg2rad
                       d_theta = half * deg2rad * d_theta
                    end if
                 end if
              call read_background_file(background_file, mode ,difpat)
                    if(Err_diffpatt) print*, trim(err_diffpatt_mess)

              !Fin de lectura
              if(allocated (resta)) deallocate(resta)
              allocate(resta(difpat%npts))
              !Algunas operaciones antes de empezar:
              check_sym = .false.
              IF(symgrpno == UNKNOWN) THEN
                symgrpno = get_sym(ok)
                IF(.NOT.ok) GO TO 999
                WRITE(op,200) 'Diffraction point symmetry is ',pnt_grp
                IF(symgrpno /= 1) THEN
                  WRITE(op,201) '  to within a tolerance of one part in ',  &
                      nint(one / tolerance)
                END IF
              ELSE
                check_sym = .true.
                CALL chk_sym(ok)
                IF(.NOT.ok) GO TO 999
              END IF

              n_high = INT(half*(th2_max-th2_min)/d_theta) + 1


              if(allocated(difpat%ycalc) ) deallocate(difpat%ycalc)
              allocate(difpat%ycalc(n_high))
              filenam = trim(infile(1:(index(infile,'.')-1)))

              open(unit=san_out,file=trim(filenam)//".out", status="replace",action="write")

              gen(1:st%npar) = st%config(1:st%npar)
              difpat%step = difpat%step * 5.0    !we enlarge the step in order to accelerate the calculation of the theoretical pattern
              vector(1:numpar) = st%state(1:numpar)
              call overlp()
              call nmcoor ()

              rpo = 1000                         !initialization of agreement factor


             !Call Set_SimAnn_StateV (sv%npar,Con,Bounds,namepar,sv%config,sv)
            !  Call SimAnneal_gen(san_out, cost3)
              Call Simanneal_Gen(cost3,san,st,san_out)


              write(unit=san_out,fmt="(/,a,/,f15.4)")  " => Final configuration (for file *.san).  Rp: ", rpo
              write(unit=san_out,fmt="(/,a)") &
                    "  NUM           Value           Name"

              do i = 1, numpar
                    write(unit=san_out,fmt="(i5,f15.4,tr10, a8)") i, statok(i) , namepar(i)
              end do

              thmin=th2_min * rad2deg
              thmax=th2_max * rad2deg
              ymax = maxval(difpat%y)
              ymini= -0.2 * ymax
              ymin = ymini - 0.5* ymax
              ymax = ymax + 0.5*ymax
              CALL getfnm(filenam, outfile, '.pgf', ok)
              OPEN(UNIT = out, FILE = outfile, STATUS = 'new')
              call DATE_AND_TIME(date, time)

                WRITE(out,'(a)')      '# .PGF (WinPLOTR Graphics file) created by FullProf:'
                WRITE (out, '(12a)')  '# ' , date(7:8),'-',date(5:6),'-',date(1:4), '  at ',&
                                            time(1:2),':',time(3:4),':',time(5:6)
                WRITE(out,'(a)')      '#'
                WRITE(out,'(a)') "# X SPACE:           1  0"
                WRITE(out,'(a)') "# MAIN LEGEND TEXT:  "//trim(filenam)
                WRITE(out,'(a)') "# X LEGEND TEXT   : 2Theta (degrees)"
                WRITE(out,'(a)') "# Y LEGEND TEXT   : Diffracted-Intensity (arb.units)"
                WRITE(out,'(a,4f14.6,2i4)') "# XMIN XMAX: " ,   thmin, thmax, thmin, thmax,1,1
                WRITE(out,'(a,4f14.6,2i4)') "# YMIN YMAX: " ,  ymin ,ymax,ymin,ymax,1,1
                WRITE(out,'(a)') "# X AND Y GRADUATIONS:   6  8  5  5"
                WRITE(out,'(a)') "# WRITE TEXT (X grad., Y grad. , Yneg. grad. , file_name):   1  1  1  1"
                WRITE(out,'(a)') "# GRID (X and Y):            0  0"
                WRITE(out,'(a)') "# FRAME FEATURES:           0.70    3    3    1    4    3"
                WRITE(out,'(a)') "# DRAW ERROR BARRS       :  N"
                WRITE(out,'(a)') "# MAIN TITLE COLOR       :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# X LEGEND COLOR         :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# Y LEGEND COLOR         :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# X GRADUATIONS COLOR    :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# Y GRADUATIONS COLOR    :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# BACKGROUND SCREEN COLOR:  RGB(240,202,166)"
                WRITE(out,'(a)') "# BACKGROUND TEXT COLOR  :  RGB(255,  0,  0)"
                WRITE(out,'(a)') "# BACKGROUND PLOT COLOR  :  RGB(255,255,255)"
                WRITE(out,'(a)') "# PLOT FRAME COLOR       :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# NUMBER OF PATTERNS:            4          "
                WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
                WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',1
                write(out,'(a,a)')          '#        FILE NAME  : ', " Observed "
                write(out,'(a,a)')          '#            TITLE  : ', " Yobs(res) "
                write(out,'(a,i10)')        '#  NUMBER OF POINTS : ',difpat%npts
                write(out,'(a)')            '#            MARKER : 4'       !open circles
                write(out,'(a,F6.1)')       '#              SIZE : 1.5'     !size 1
                write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0,  0,255)' !Red
                write(out,'(a,i6)')         '#             STYLE : 0'       !Points non continuous line
                write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
                write(out,'(a,i6)')         '#        DATA: X Y  '
                do i=1,difpat%npts
                  write(unit=out,fmt="(f10.6,3f14.6)") difpat%x(i), difpat%y(i)
                end do
                WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
                WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',2
                write(out,'(a,a)')          '#        FILE NAME  : ', " Calculated "
                write(out,'(a,a)')          '#            TITLE  : ', " Ycal(res) "
                write(out,'(a,i10)')        '#  NUMBER OF POINTS : ', n_high
                write(out,'(a)')            '#            MARKER : 4'       !open circles
                write(out,'(a,F6.1)')       '#              SIZE : 0.0'     !size
                write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0, 0,0)' !Red
                write(out,'(a,i6)')         '#             STYLE : 1'       !Continuous line
                write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
                write(out,'(a,i6)')         '#        DATA: X Y  '
                do i=1,n_high
                  theta = thmin +(i-1)*d_theta*rad2deg * 2
                  write(unit = out,fmt = "(f10.6,3f14.6)") theta,  ycalcdef(i)
                end do
                p_resta= n_high-1
                WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
                WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',3
                write(out,'(a,a)')          '#        FILE NAME  : ', " Difference"
                write(out,'(a,a)')          '#            TITLE  : ', " Yobs-Ycal "
                write(out,'(a,i10)')        '#  NUMBER OF POINTS : ', p_resta
                write(out,'(a)')            '#            MARKER : 4'       !open circles
                write(out,'(a,F6.1)')       '#              SIZE : 0.0'     !size 1
                write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  255,  0,0)' !?
                write(out,'(a,i6)')         '#             STYLE : 1'       !Points non continuous line
                write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
                        write(out,'(a,i6)')         '#        DATA: X Y  '
                thmin_o = difpat%xmin * deg2rad      !rad

                if    ( th2_min/= thmin_o ) then
                  punts= INT(half*(th2_min-thmin_o)/d_theta)
                end if

                do i=1,n_high-1
                     resta(i) =(difpat%y(i+punts) - ycalcdef(i)) + ymini
                     write(unit=out,fmt="(2f15.7)") difpat%x(i+punts) , resta(i)
                end do
                WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
                WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',4
                write(out,'(a,a)')          '#        FILE NAME  : ', " Bragg_position "
                write(out,'(a,a)')          '#            TITLE  : ', " Bragg_position "
                write(out,'(a,i10)')        '#  NUMBER OF POINTS : '  , d_punt
                write(out,'(a)')            '#            MARKER : 8'       !open circles
                write(out,'(a,F6.1)')       '#              SIZE : 3.0'     !size 0
                write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0,  128,  0) ' !black
                write(out,'(a,i6)')         '#             STYLE : 0'       !Continuous line
                write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
                write(out,'(a,i6)')         '#        DATA: X Y  '

                   deg = ymini * 0.25
                   aa=h_min
                   Do a=h_min, h_max
                      bb=k_min
                     do b=k_min, k_max
                        cc= zero
                       do c=1, 30
                         if (dos_theta(aa,bb,cc) /= zero) then
                            write(unit= out,fmt = "(2F15.7, 5x, a, 3i3, a, i3)") dos_theta(aa, bb, cc)  , deg,  "(", aa,bb,cc,")", 1
                         end if
                         cc=cc+1
                       End do
                        bb=bb+1
                     End do
                      aa=aa+1
                  End do

                write(out,"(a)") "# END OF FILE "
                close(unit=out)




          Case (2) !SIMPLEX

              !Lectura del  pattern experimental y de los scattering factors:

              call  read_pattern(dfile,difpat,fmode )
                if(Err_diffpatt) then
                  print*, trim(err_diffpatt_mess)
                else
                  if (th2_min == 0 .and.  th2_max == 0  .and. d_theta == 0) then
                     th2_min =  difpat%xmin    ! if not specified in input file, take the values from the observed pattern
                     th2_max =  difpat%xmax
                     d_theta =  difpat%step
                     th2_min = th2_min * deg2rad
                     th2_max = th2_max * deg2rad
                     d_theta = half * deg2rad * d_theta
                  end if
                end if
              call read_background_file(background_file, mode ,difpat)
                    if(Err_diffpatt) print*, trim(err_diffpatt_mess)

              !Fin de lectura

              if(allocated (resta)) deallocate(resta)
              allocate(resta(difpat%npts))

              !Algunas operaciones antes de empezar:
              check_sym = .false.
              IF(symgrpno == UNKNOWN) THEN
                symgrpno = get_sym(ok)
                IF(.NOT.ok) GO TO 999
                WRITE(op,200) 'Diffraction point symmetry is ',pnt_grp
                IF(symgrpno /= 1) THEN
                  WRITE(op,201) '  to within a tolerance of one part in ',  &
                      nint(one / tolerance)
                END IF
              ELSE
                check_sym = .true.
                CALL chk_sym(ok)
                IF(.NOT.ok) GO TO 999
              END IF

              n_high = INT(half*(th2_max-th2_min)/d_theta) + 1


              if(allocated(difpat%ycalc) ) deallocate(difpat%ycalc)
              allocate(difpat%ycalc(n_high))
              filenam = trim(infile(1:(index(infile,'.')-1)))

              gen(1:n_plex) = v_plex(1:n_plex)
              call overlp()
              call nmcoor ()

              rpo = 1000                         !initialization of agreement factor
              rpl = 0
              do i=1, n_plex                       !creation of the step sizes
                steplex(i) = 0.05 * v_plex(i)
              end do
              open (unit=23, file='nelder_mess.out', status='replace', action='write')
              call  Nelder_Mead_Simplex( F_cost,n_plex ,v_plex(1:n_plex) , &
                                         steplex(1:n_plex), var_plex(1:n_plex), rpl, opt_c, ipr=23)
              write(*,*)'Rp', rpo
              write(*,*) '________________________________________________'
              write(*,'(3a)') ' Parameter     refined value     sigma'
              write(*,*) '________________________________________________'
              do i = 1, numpar
                    write(*,*)  namepar(i)  ,statok(i) , var_plex(i)
              end do

              thmin=th2_min * rad2deg
              thmax=th2_max * rad2deg
              ymax = maxval(difpat%y)
              ymini= -0.2 * ymax
              ymin = ymini - 0.5* ymax
              ymax = ymax + 0.5*ymax

              CALL getfnm(filenam, outfile, '.pgf', ok)
              if (ok) then
                OPEN(UNIT = out, FILE = outfile, STATUS = 'replace')
                call DATE_AND_TIME(date, time)

                WRITE(out,'(a)')      '# .PGF (WinPLOTR Graphics file) created by FullProf:'
                WRITE (out, '(12a)')  '# ' , date(7:8),'-',date(5:6),'-',date(1:4), '  at ',&
                                            time(1:2),':',time(3:4),':',time(5:6)
                WRITE(out,'(a)')      '#'
                WRITE(out,'(a)') "# X SPACE:           1  0"
                WRITE(out,'(a)') "# MAIN LEGEND TEXT:  "//trim(filenam)
                WRITE(out,'(a)') "# X LEGEND TEXT   : 2Theta (degrees)"
                WRITE(out,'(a)') "# Y LEGEND TEXT   : Diffracted-Intensity (arb.units)"
                WRITE(out,'(a,4f14.6,2i4)') "# XMIN XMAX: " ,   thmin, thmax, thmin, thmax,1,1
                WRITE(out,'(a,4f14.6,2i4)') "# YMIN YMAX: " ,  ymin ,ymax,ymin,ymax,1,1
                WRITE(out,'(a)') "# X AND Y GRADUATIONS:   6  8  5  5"
                WRITE(out,'(a)') "# WRITE TEXT (X grad., Y grad. , Yneg. grad. , file_name):   1  1  1  1"
                WRITE(out,'(a)') "# GRID (X and Y):            0  0"
                WRITE(out,'(a)') "# FRAME FEATURES:           0.70    3    3    1    4    3"
                WRITE(out,'(a)') "# DRAW ERROR BARRS       :  N"
                WRITE(out,'(a)') "# MAIN TITLE COLOR       :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# X LEGEND COLOR         :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# Y LEGEND COLOR         :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# X GRADUATIONS COLOR    :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# Y GRADUATIONS COLOR    :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# BACKGROUND SCREEN COLOR:  RGB(240,202,166)"
                WRITE(out,'(a)') "# BACKGROUND TEXT COLOR  :  RGB(255,  0,  0)"
                WRITE(out,'(a)') "# BACKGROUND PLOT COLOR  :  RGB(255,255,255)"
                WRITE(out,'(a)') "# PLOT FRAME COLOR       :  RGB(  0,  0,  0)"
                WRITE(out,'(a)') "# NUMBER OF PATTERNS:            4          "
                WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
                WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',1
                write(out,'(a,a)')          '#        FILE NAME  : ', " Observed "
                write(out,'(a,a)')          '#            TITLE  : ', " Yobs(res) "
                write(out,'(a,i10)')        '#  NUMBER OF POINTS : ',difpat%npts
                write(out,'(a)')            '#            MARKER : 4'       !open circles
                write(out,'(a,F6.1)')       '#              SIZE : 1.5'     !size 1
                write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0,  0,255)' !Red
                write(out,'(a,i6)')         '#             STYLE : 0'       !Points non continuous line
                write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
                write(out,'(a,i6)')         '#        DATA: X Y  '
                do i=1,difpat%npts
                  write(unit=out,fmt="(f10.6,3f14.6)") difpat%x(i), difpat%y(i)
                end do
                WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
                WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',2
                write(out,'(a,a)')          '#        FILE NAME  : ', " Calculated "
                write(out,'(a,a)')          '#            TITLE  : ', " Ycal(res) "
                write(out,'(a,i10)')        '#  NUMBER OF POINTS : ', n_high
                write(out,'(a)')            '#            MARKER : 4'       !open circles
                write(out,'(a,F6.1)')       '#              SIZE : 0.0'     !size
                write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0, 0,0)' !Red
                write(out,'(a,i6)')         '#             STYLE : 1'       !Continuous line
                write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
                write(out,'(a,i6)')         '#        DATA: X Y  '
                do i=1,n_high
                  theta = thmin +(i-1)*d_theta *rad2deg * 2
                  write(unit = out,fmt = "(f10.6,3f14.6)") theta,  ycalcdef(i)
                end do
                p_resta= n_high-1
                WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
                WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',3
                write(out,'(a,a)')          '#        FILE NAME  : ', " Difference"
                write(out,'(a,a)')          '#            TITLE  : ', " Yobs-Ycal "
                write(out,'(a,i10)')        '#  NUMBER OF POINTS : ', p_resta
                write(out,'(a)')            '#            MARKER : 4'       !open circles
                write(out,'(a,F6.1)')       '#              SIZE : 0.0'     !size 1
                write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  255,  0,0)' !?
                write(out,'(a,i6)')         '#             STYLE : 1'       !Points non continuous line
                write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
                        write(out,'(a,i6)')         '#        DATA: X Y  '
                thmin_o = difpat%xmin * deg2rad      !rad

                if    ( th2_min/= thmin_o ) then
                  punts= INT(half*(th2_min-thmin_o)/d_theta)
                end if

                do i=1,n_high-1
                     resta(i) =(difpat%y(i+punts) - difpat%ycalc(i))+ymini
                     write(unit=out,fmt="(2f15.7)") difpat%x(i+punts) , resta(i)!, difpat%y(i+punts)    ,difpat%ycalc(i), punts, ymini
                end do

                WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
                WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',4
                write(out,'(a,a)')          '#        FILE NAME  : ', " Bragg_position "
                write(out,'(a,a)')          '#            TITLE  : ', " Bragg_position "
                write(out,'(a,i10)')        '#  NUMBER OF POINTS : '  , d_punt
                write(out,'(a)')            '#            MARKER : 8'       !open circles
                write(out,'(a,F6.1)')       '#              SIZE : 3.0'     !size 0
                write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0,  128,  0) ' !black
                write(out,'(a,i6)')         '#             STYLE : 0'       !Continuous line
                write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
                write(out,'(a,i6)')         '#        DATA: X Y  '

               else
                write(*,*) 'The outfile cannot be created'
               end if
                     deg = ymini * 0.25
                     aa=h_min
                   Do a=h_min, h_max
                      bb=k_min
                     do b=k_min, k_max
                        cc= zero
                       do c=1, 30
                         if (dos_theta(aa,bb,cc) /= zero) then
                            write(unit= out,fmt = "(2F15.7, 5x, a, 3i3, a, i3)") dos_theta(aa, bb, cc)  , deg,  "(", aa,bb,cc,")", 1
                         end if
                         cc=cc+1
                       End do
                        bb=bb+1
                     End do
                      aa=aa+1
                  End do

                write(out,"(a)") "# END OF FILE "
                close(unit=out)


          Case default

                print*,"problems reading mode "


          End select



      ending = .true.

      close(unit = san_out)
      999 IF(cfile) CLOSE(UNIT = cntrl)
      IF(ok .AND. ending) THEN
        WRITE(op,100) ' => FAULTS ended normally....'
      ELSE
        WRITE(op,100) ' => FAULTS was terminated abnormally!'
      END IF
      100 FORMAT(1X, a)
      101 FORMAT(1X, i3)
      200 FORMAT(1X, 2A)
      201 FORMAT(1X, a, i6)

   END PROGRAM FAULTS


    Subroutine Write_FST(fst_file,v,cost)
       Use CFML_GlobalDeps, only: Cp
       character(len=*),     intent(in):: fst_file
       real(kind=cp),dimension(:),    intent(in):: v
       real(kind=cp),                 intent(in):: cost
       return
    End Subroutine Write_FST
