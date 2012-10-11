
  Module Diff_ref

      use diffax_mod
      use CFML_GlobalDeps,            only : sp
      use CFML_Diffraction_Patterns , only : diffraction_pattern_type
      use CFML_Optimization_General,  only : Opt_Conditions_Type
      use read_data,                  only : opti

      implicit none

      private

      !public subroutines
      public :: scale_factor, scale_factor_lmq, Write_Prf

      contains
!________________________________________________________________________________________________________________________

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

       Subroutine scale_factor_lmq(pat, fvec, chi2)
         type (diffraction_pattern_type)  , intent (in out) :: pat
         Real (Kind=cp),Dimension(:),         Intent(   Out):: fvec
         real                                  ,Intent( out):: chi2
         real                                               :: a ,  c, r
         real                                               :: thmin, thmax  , thmin_o
         integer                                            :: n_high, punts=0
         integer                                            :: i,j

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
         integer                        , intent(in    ) :: punts
         real                           , intent(   out) :: r
         Real (Kind=cp),Dimension(:),     Intent(   out) :: fvec
         real,                            Intent(   out) :: chi2
         !----
         real                                            :: a,b,c
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
         write (*,*) r, chi2 , punts-opti%npar
         return
       End subroutine calc_par_lmq

       Subroutine scale_factor (pat,r, chi2)
         type (diffraction_pattern_type)  , intent (in out) :: pat
         real                             , intent (   out) :: r
         real                             , intent (   out) :: chi2
         real                                               :: a ,  c
         real                                               :: thmin, thmax  , thmin_o
         integer                                            :: n_high, punts=0
         integer                                            :: i,j

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

         call calc_par(pat, punts, r, chi2)
         return
      End subroutine scale_factor

!____________________________________________________________________________________________________________________________

    Subroutine calc_par (pat, punts, r, chi2)
      type (diffraction_pattern_type), intent(in    ) :: pat
      integer                        , intent(in    ) :: punts
      real                           , intent(   out) :: r
      real                           , intent(   out) :: chi2
      real                                            :: a,b,c
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
        c = c + ((pat%y(j) - pat%ycalc(j))**2/pat%sigma(j))
      end do do_c

      chi2= c/(punts-opti%npar)
      write (*,*) r, chi2 , punts-opti%npar
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
    use CFML_Optimization_LSQ,      only : Levenberg_Marquardt_Fit
    use CFML_LSQ_TypeDef,           only : LSQ_Conditions_type, LSQ_State_Vector_Type
    use diffax_mod
    use read_data,                  only : crys, read_structure_file, length,   opti , cond, vs
    use diffax_calc,                only : salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, optimz,point,  &
                                           gospec, gostrk, gointr,gosadp, getfnm, nmcoor
    use Diff_ref,                   only : scale_factor, scale_factor_lmq, Write_Prf

    implicit none

    public  :: F_cost, Cost_LMQ !Cost3
    type (diffraction_pattern_type),  save         :: difpat
    type (State_Vector_Type), public, save         :: st

    contains


!   Subroutine  Cost3(vref,rp)      !Simulated annealing
!
!         real,   dimension (:), intent (in    ) :: vref
!         real,                  intent (   out) :: rp
!         logical                                :: ok
!         integer                                :: n , j ,i,k  , label , a, b
!         real, dimension(80)                    :: shift, state, menor=1 , multi,tar
!         character                              :: cons
!
!
!         do i= 1, numpar                          !shift calculation
!
!               shift(i) = vref(i) - gen(i)
!
!         end do
!
!       !*******RESTRICTIONS*******
!
!        do i = 1, numpar
!
!            if (index (namepar(i) , 'alpha' ) == 1 )  then       !To avoid negative values of alpha
!               ! read (unit = namepar(i)(6:7), fmt = "(2i1)" ) b,a
!               ! if  (l_alpha(a,b) .le. menor(b) ) then
!               !        menor(b) = l_alpha(a, b)
!               !        multi(b) = mult(i)
!               ! end if
!                state(i) = vector(i) +  mult(i) * shift(pnum(i))
!                if (state(i) < zero ) then
!                       write(*,*) 'Attention, shift was higher than alpha:  new shift applied'
!                       shift(pnum(i)) = shift(pnum(i))/2
!                       if (state(i) < zero ) then
!                       shift(pnum(i)) = -shift(pnum(i))
!                       end if
!                end if
!
!            end if
!        end do
!
!
!        do i=1, numpar   !assignment of new values and more restrictions
!
!           state(i) = vector(i) +  mult(i) * shift(pnum(i))
!
!           if (index (namepar(i) , 'Biso' ) == 1 .and. state(i) < 0 )   state(i) = (-1.0) * state(i)  !Biso only >0
!           if (index (namepar(i) , 'v' ) == 1 .and. state(i) < 0 )      state(i) = (-1.0) * state(i)  !v only <0
!           if (index (namepar(i) , 'Dg') == 1 .or. index (namepar(i) , 'Dl') == 1 .or.  index (namepar(i) , 'u') == 1 .or.  &
!               index (namepar(i) , 'w') == 1 .or. index (namepar(i) , 'x') == 1   ) then  !only >0
!
!             if ( state(i) < 0)  state(i) =  -state(i)
!
!           end if
!
!        End do
!
!        vector(:) = state(:)
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!
!         tar(:)=0                                 !generation of gen
!          do i=1, numpar
!            if (tar(pnum(i)) == 0 ) gen(pnum(i)) = state(i)
!            tar(pnum(i)) = 1
!          end do
!
!
!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         !////////CALCULATED PATTERN VARIABLES ASSIGNMENT////////////////////////////////////
!
!         do i=1, numpar
!
!           write(*,*) namepar(i) ,  state(i)
!           if (namepar(i) ==  'u')    pv_u  = state(i)
!           if (namepar(i) ==  'v')    pv_v  = state(i)
!           if (namepar(i) ==  'w')    pv_w  = state(i)
!           if (namepar(i) ==  'x')    pv_x  = state(i)
!           if (namepar(i) ==  'Dg')   pv_dg = state(i)
!           if (namepar(i) ==  'Dl')   pv_dl = state(i)
!           if (namepar(i) ==  'cell_a')    cell_a = state(i)
!           if (namepar(i) ==  'cell_b')    cell_b = state(i)
!           if (namepar(i) ==  'cell_c')    cell_c = state(i)
!           if (namepar(i) ==  'num_layers')  l_cnt= state(i)
!           if (index( namepar(i) ,'cell_gamma') == 1)   cell_gamma = state(i)
!
!           do j=1, n_layers
!             do k=1, n_atoms
!               if (index (namepar(i) , 'pos_x' )== 1)     then
!                   read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
!                   a_pos(1,a,b)  = state(i) * pi2
!               end if
!               if (index (namepar(i) ,'pos_y' )== 1)    then
!                   read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
!                   a_pos(2,a,b)  = state(i) * pi2
!               end if
!               if (index (namepar(i), 'pos_z' ) == 1 )   then
!                   read (unit = namepar(i)(6:7), fmt = "(2i1)" ) a,b
!                   a_pos(3,a,b)  = state(i) * pi2
!
!               end if
!               if (index (namepar(i) ,'Biso')== 1) then
!                   read (unit = namepar(i)(5:6), fmt = "(2i1)" ) a,b
!                   a_b(a,b)  = state(i)
!               end if
!               if (index( namepar(i) ,  'alpha' ) == 1)    then
!                   read (unit = namepar(i)(6:7), fmt = "(2i1)" ) b,a
!                   l_alpha(a,b)  = state(i)
!
!               end if
!               if (index (namepar(i), 'tx' )== 1 )    then
!                   read (unit = namepar(i)(3:4), fmt = "(2i1)" ) b,a
!                   l_r(1,a,b)  = state(i)
!               end if
!               if (index (namepar(i), 'ty' )== 1 )    then
!                   read (unit = namepar(i)(3:4), fmt = "(2i1)" )b,a
!                   l_r(2,a,b)  = state(i)
!               end if
!               if (index (namepar(i), 'tz' ) == 1)     then
!                   read (unit = namepar(i)(3:4), fmt = "(2i1)" ) b,a
!                   l_r(3,a,b)  = state(i)
!               end if
!             end do
!           end do
!         end do
!
!
!
!         ok = .true.
!
!         if ((conv_d == 1 .or. numcal== 0) .and. ok) ok = get_g()
!         if ((conv_d == 1  .or.  numcal== 0 .or. conv_e==1) .and. (ok .and. rndm) ) ok = getlay()
!         if ((conv_c == 1 .or. numcal== 0) .and. ok ) CALL sphcst()
!         if (numcal == 0 .and. ok) CALL detun()
!         if ((conv_b == 1 .or. conv_c == 1 .or. conv_d == 1  .or. numcal == 0 .or. &
!             conv_e == 1 .or. conv_f == 1) .and. ok) CALL optimz(infile, ok)
!         IF(.NOT.ok) then
!          IF(cfile) CLOSE(UNIT = cntrl)
!          return
!         END IF
!         CALL gospec(infile,outfile,ok)
!         call scale_factor (difpat,rp)
!
!         numcal = numcal + 1
!         write(*,*) ' => Calculated Rp    :   ' , rp
!         write(*,*) ' => Best Rp up to now:   ' , rpo
!
!         if (rp .LT. rpo ) then                  !To keep calculated intensity for the best value of rp
!           rpo = rp
!           statok(1:numpar) = state( 1:numpar)
!           write(*,*) 'writing calculated best pattern up to now. Rp :       ', rpo
!           do j = 1, n_high
!             ycalcdef(j) = difpat%ycalc(j)
!           end do
!         end if
!
!         ok = .true.
!         IF(cfile) CLOSE(UNIT = cntrl)
!         return
!   End subroutine Cost3

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

    Subroutine Cost_LMQ(m,npar,v,fvec,iflag)   !Levenberg Marquardt
      Integer,                       Intent(In)    :: m !is the number of observations (Num_spots)
      Integer,                       Intent(In)    :: npar !is the number of free parameters
      Real (Kind=cp),Dimension(:),   Intent(In)    :: v !List of free parameters values
      Real (Kind=cp),Dimension(:),   Intent(In Out):: fvec   !Residuals Num_spots
      Integer,                       Intent(In Out):: iflag  !=0 for printing, 1: Full calculation, 2: Calculation of derivatives
                                                             ! If modified to a negative number the algorithm stops.
      !local variables
      logical                  :: ok
      integer                  :: j ,i, k, a, b
      real, dimension(max_npar):: shift, state
      real                     :: chi2
      integer, save            :: iter=0

      fvec=0.0
      !chi2=chi2o
      write(*,*)"--------FCOST-------"

      do i= 1, opti%npar
         shift(i) = v(i) - vector(i)
      end do

      do i = 1, crys%npar
         state(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
      end do


      !******* state(i) corrections *******

      do i = 1, crys%npar
         if (index (vs%nampar(i) , 'alpha' ) == 1 .and. (state(i) < zero .or. state(i) > 1)) then
                write(*,*) 'Attention, shift was higher/lower than accepted values for alpha:  new shift applied'
                write(*,*) "alpha before" , vs%nampar(i), state(i), shift(vs%code(i)) , crys%vlim1(crys%p(i)) ,&
                            crys%vlim2(vs%code(i))
                shift(crys%p(i)) = shift(crys%p(i))/2
                state(i) = crys%list(i) +  mult(i) * shift(vs%code(i))
                write(*,*) "alpha after" , vs%nampar(i), state(i), shift(vs%code(i))
                do j=1,crys%npar
                  write(*,*) "j", j, vs%nampar(j)
                  if (index (vs%nampar(j) , 'alpha' ) == 1 .and.  vs%code(i) == vs%code(j)) &
                      state(j) = crys%list(j) +  mult(j) * shift(vs%code(i))
                  write(*,*) "other alpha" , vs%nampar(j), state(j), shift(vs%code(i))
                end do

                if (state(i) < zero .or. state(i) > zero) then
                  shift(vs%code(i)) = - shift (vs%code(i))
                  state(i) = crys%list(i) +  mult(i) * shift(vs%code(i))
                  write(*,*) "alpha again" , vs%nampar(i), state(i), shift(vs%code(i))
                end if
         else if (index (vs%nampar(i),'Biso' ) == 1 .and. state(i) .lt. 0 ) then
                  state(i) = (-1.0) * state(i)  !Biso only >0
         else if (index (vs%nampar(i),'v' ) == 1 .and. state(i) .gt. 0 ) then
                  state(i) = (-1.0) * state(i)  !v only <0
         else if (index (vs%nampar(i),'Dg') == 1 .or. index (vs%nampar(i), 'Dl') == 1 .or. index (vs%nampar(i),'u')  == 1 .or. &
                  index (vs%nampar(i), 'w') == 1 .or. index (vs%nampar(i), 'x') == 1 ) then  !only >0
                   state(i) = (-1.0) * state(i)
         else
            cycle
         end if
      End do

      !update  (only if iflag = 1)
      if(iflag == 1) then
         crys%list(:) = state(:)
         do i=1, opti%npar                 !vector upload
           vector(i) = v(i)
         end do
      end if

      !////////CALCULATED PATTERN VARIABLES ASSIGNMENT////////////////////////////////////

      do i=1, crys%npar
        if (vs%nampar(i) ==  'u')    pv_u  = state(i)
        if (vs%nampar(i) ==  'v')    pv_v  = state(i)
        if (vs%nampar(i) ==  'w')    pv_w  = state(i)
        if (vs%nampar(i) ==  'x')    pv_x  = state(i)
        if (vs%nampar(i) ==  'Dg')   pv_dg = state(i)
        if (vs%nampar(i) ==  'Dl')   pv_dl = state(i)
        if (vs%nampar(i) ==  'cell_a')    cell_a = state(i)
        if (vs%nampar(i) ==  'cell_b')    cell_b = state(i)
        if (vs%nampar(i) ==  'cell_c')    cell_c = state(i)
        if (vs%nampar(i) ==  'num_layers')  l_cnt = state(i)
        if (index( vs%nampar(i) ,'cell_gamma') == 1)   cell_gamma = state(i)

        do j=1, n_layers
          do k=1, n_atoms
            if (index (vs%nampar(i) , 'pos_x' )== 1)     then
                read (unit = vs%nampar(i)(6:7), fmt = "(2i1)" ) a,b
                a_pos(1,a,b)  = state(i) * pi2
            end if
            if (index (vs%nampar(i) ,'pos_y' )== 1)    then
                read (unit = vs%nampar(i)(6:7), fmt = "(2i1)" ) a,b
                a_pos(2,a,b)  = state(i) * pi2
            end if
            if (index (vs%nampar(i), 'pos_z' ) == 1 )   then
                read (unit = vs%nampar(i)(6:7), fmt = "(2i1)" ) a,b
                a_pos(3,a,b)  = state(i) * pi2
            end if
            if (index (vs%nampar(i),'Biso')==1) then
                read (unit = vs%nampar(i)(5:6), fmt = "(2i1)" ) a,b
                a_b(a,b)  = state(i)

            end if
            if (index( vs%nampar(i) ,  'alpha' ) == 1)    then
                read (unit = vs%nampar(i)(6:7), fmt = "(2i1)" ) b,a
                l_alpha(a,b)  = state(i)

            end if
            if (index (vs%nampar(i), 'tx' )== 1 )    then
                read (unit = vs%nampar(i)(3:4), fmt = "(2i1)" ) b,a
                l_r(1,a,b)  = state(i)
            end if
            if (index (vs%nampar(i), 'ty' )== 1 )    then
                read (unit = vs%nampar(i)(3:4), fmt = "(2i1)" ) b,a
                l_r(2,a,b)  = state(i)
            end if
            if (index (vs%nampar(i), 'tz' ) == 1)     then
                read (unit = vs%nampar(i)(3:4), fmt = "(2i1)" ) b,a
                l_r(3,a,b)  = state(i)
            end if
          end do
        end do
      end do

      ok = .true.
      if ((conv_d == 1 .or.  numcal == 0) .and. ok ) ok = get_g()
      if ((conv_d == 1 .or.  numcal == 0  .or.  conv_e==1) .and. (ok .AND. rndm)) ok = getlay()
      if ((conv_c == 1 .or.  numcal == 0) .and. ok ) CALL sphcst()
      if ( numcal == 0 .and. ok ) CALL detun()
      if ((conv_b == 1 .or.  conv_c == 1  .or. conv_d == 1 .or. conv_e==1   .or. &
           numcal == 0 .or.  conv_f ==1)  .and. ok ) CALL optimz(infile, ok)

      IF(.NOT. ok) then
        IF(cfile) CLOSE(UNIT = cntrl)
        return
      END IF
      CALL gospec(infile,outfile,ok)

      if(.not. ok) then
        print*, "Error calculating spectrum, please check input parameters"
      else
        call scale_factor_lmq(difpat,fvec, chi2)
      end if
      numcal = numcal + 1

      if(iflag == 0) then
        iter = iter + cond%nprint
        write(*,"(a,i4)")  " => Iteration number ",iter
        do i=1, crys%npar
          write(*,"(a,f14.5)")  "  ->  "//vs%nampar(i), state(i)
        end do
        write(*,*) ' => Calculated Chi2    :   ' , chi2
        write(*,*) ' => Best Chi2 up to now:   ' , chi2o
      end if


      if (chi2 < chi2o ) then                  !To keep calculated intensity for the best value of rplex
        chi2o = chi2
        statok(1:crys%npar) = state( 1:crys%npar)
        write(*,*)  ' => Writing the best calculated pattern up to now. Chi2 : ', chi2o
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

    End subroutine Cost_LMQ
!--------------------------------------------------------------------------------------------------------------------------------------------------

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
      integer                 :: j ,i, k, a, b
      real, dimension(300)    :: shift, state
      real                    :: chi2

      write(*,*)"--------FCOST-------"

       !--- to avoid warnings
      if(present(g)) g=0.0
      do i= 1, opti%npar
            shift(i) = v(i) - vector(i)
      end do

      do i = 1, crys%npar
         state(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
      end do

      !******* state(i) corrections *******

      do i = 1, crys%npar
         if (index (namepar(i) , 'alpha' ) == 1 .and. (state(i) < zero .or. state(i) > 1)) then
                write(*,*) 'Attention, shift was higher/lower than accepted values for alpha:  new shift applied'
                write(*,*) "alpha before" , namepar(i), state(i), shift(crys%p(i)) , crys%vlim1(crys%p(i)) ,&
                            crys%vlim2(crys%p(i))
                shift(crys%p(i)) = shift(crys%p(i))/2
                state(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
                write(*,*) "alpha after" , namepar(i), state(i), shift(crys%p(i))
                do j=1,crys%npar
                  write(*,*) "j", j, namepar(j)
                  if (index (namepar(j) , 'alpha' ) == 1 .and.  crys%p(i) == crys%p(j)) &
                      state(j) = crys%list(j) +  mult(j) * shift(crys%p(i))
                  write(*,*) "other alpha" , namepar(j), state(j), shift(crys%p(i))
                end do

                if (state(i) < zero .or. state(i) > zero) then
                  shift(crys%p(i)) = - shift (crys%p(i))
                  state(i) = crys%list(i) +  mult(i) * shift(crys%p(i))
                  write(*,*) "alpha again" , namepar(i), state(i), shift(crys%p(i))
                end if
         else if (index (namepar(i),'Biso' ) == 1 .and. state(i) .lt. 0 ) then
                  state(i) = (-1.0) * state(i)  !Biso only >0
         else if (index (namepar(i),'v' ) == 1 .and. state(i) .gt. 0 ) then
                  state(i) = (-1.0) * state(i)  !v only <0
         else if (index (namepar(i),'Dg') == 1 .or. index (namepar(i), 'Dl') == 1 .or. index (namepar(i),'u')  == 1 .or. &
                  index (namepar(i), 'w') == 1  .or. index (namepar(i), 'x') == 1 ) then  !only >0
                   state(i) = (-1.0) * state(i)
         else
            cycle
         end if
      End do

      !update

      crys%list(:) = state(:)

      do i=1, opti%npar
             vector(i) = v(i)
      end do

      do i=1, crys%npar
        write(*,*)  namepar(i), state(i)
      end do


      !////////CALCULATED PATTERN VARIABLES ASSIGNMENT////////////////////////////////////

      do i=1, crys%npar

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
      if(.not. ok) then
        print*, "Error calculating spectrum, please check input parameters"
      else
        call scale_factor(difpat,rplex, chi2)
      end if
      numcal = numcal + 1
      write(*,*) ' => Calculated Rp    :   ' , rplex
      write(*,*) ' => Best Rp up to now:   ' , rpo
      write(*,*) var_plex(:)

      if (rplex < rpo ) then                  !To keep calculated intensity for the best value of rplex
        rpo = rplex
        statok(1:crys%npar) = state( 1:crys%npar)

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
                                              getword, err_string, err_string_mess, getnum, Ucase, lcase
     use CFML_Diffraction_patterns,    only : read_pattern , diffraction_pattern_type , err_diffpatt, err_diffpatt_mess,  &
                                              read_background_file
     use CFML_Simulated_Annealing
     use CFML_Optimization_General,    only : Nelder_Mead_Simplex,  Opt_Conditions_Type, Local_Optimize
     use CFML_Crystal_Metrics,         only : Set_Crystal_Cell, Crystal_Cell_Type
     use CFML_Optimization_LSQ,        only : Levenberg_Marquardt_Fit
     use CFML_LSQ_TypeDef,             only : LSQ_Conditions_type, LSQ_State_Vector_Type
     use diffax_mod
     use read_data,                    only : new_getfil, read_structure_file, length,   &
                                              crys, opti, opsan, cond, Vs
     use diffax_calc ,                 only : salute , sfc, get_g, get_alpha, getlay , sphcst, dump, detun, optimz,point,  &
                                              gospec, gostrk, gointr,gosadp, chk_sym, get_sym, overlp, nmcoor , getfnm
     use Diff_ref ,                    only : scale_factor, scale_factor_lmq, Write_Prf
     use dif_ref,                      only :  difpat , F_cost, st, cost_LMQ !, cost3

     implicit none


      real                                                    :: rpl,  theta ,thmin, thmax , thmin_o  , ymax, ymini , ymin ,deg
      LOGICAL                                                 :: ok, ending , gol, p_ok
      INTEGER*4                                               :: i ,n ,j, l , ier , fn_menu,a,b,c ,aa,bb,cc, p_resta ,e
      character(len=100)                                      :: pfile, bfile , bmode
      character(len=100)                                      :: pmode, filenam
      integer, parameter                                      :: out = 25
      real, dimension (:), allocatable                        :: resta
      character(len=10)                                       :: time, date
      Real (Kind=cp)                                          :: chi2     !final Chi2
      character(len=3000)                                     :: infout   !Information about the refinement (min length 256)

      pi = four * ATAN(one)
      pi2 = two * pi
      deg2rad = pi / one_eighty
      rad2deg = one_eighty / pi

      ending = .false.
      ok = .true.

      CALL salute()


      sfname = 'data.sfc'
      cntrl=ip

      call new_getfil(infile, st,   gol)                           ! parametros para calculo
      opsan%Cost_function_name="R-factor"

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
                IF(.NOT. ok) GO TO 999
                WRITE(op,200) 'Diffraction point symmetry is ',pnt_grp
                IF(symgrpno /= 1) THEN
                  WRITE(op,201) '  to within a tolerance of one part in ',  &
                      nint(one / tolerance)
                END IF
              ELSE
                check_sym = .true.
                CALL chk_sym(ok)
                IF(.NOT. ok) GO TO 999
              END IF

              n_high = INT(half*(th2_max-th2_min)/d_theta) + 1

              call overlp()

              if(allocated(difpat%ycalc) ) deallocate(difpat%ycalc)
              allocate(difpat%ycalc(n_high))

              IF(ok) ok = get_g()
              IF(ok .AND. rndm) ok = getlay()
              IF(ok) CALL sphcst()
              IF(ok) CALL detun()
              IF(.NOT. ok) GO TO 999
! See if there are any optimizations we can do
              IF(ok) CALL optimz(infile, ok)

! What type of intensity output does the user want?
           10  IF(ok) THEN
              ! WRITE(op,100) 'Enter function number: '
              ! WRITE(op,100) '0 POINT, 1 STREAK, 2 INTEGRATE, 3 POWDER PATTERN, 4 SADP : '
              ! READ(cntrl,*,END=999) n
                n=3 !Powder pattern
             END IF

             call dump (infile, p_ok)
             call nmcoor ()
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

                     n_high = nINT(half*(th2_max-th2_min)/d_theta + 1.0)

                     Do j = 1, n_high
                         ycalcdef(j) = brd_spc(j)
                         !write(*,*) ycalcdef(j)
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
                        write(unit = out,fmt = '(3f12.4)')thmin, d_theta,thmax
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

!!!       Case (1)     !SIMULATED ANNEALING
!!!
!!!           !Lectura del  pattern experimental y de los scattering factors:
!!!
!!!
!!!           call  read_pattern (dfile,difpat,fmode )
!!!              if(Err_diffpatt) then
!!!                 print*, trim(err_diffpatt_mess)
!!!              else
!!!                 if (th2_min == 0 .and.  th2_max == 0  .and. d_theta == 0) then
!!!                    th2_min =  difpat%xmin    ! if not specified in input file, take the values from the observed pattern
!!!                    th2_max =  difpat%xmax
!!!                    d_theta =  difpat%step
!!!                    th2_min = th2_min * deg2rad
!!!                    th2_max = th2_max * deg2rad
!!!                    d_theta = half * deg2rad * d_theta
!!!                 end if
!!!              end if
!!!           call read_background_file(background_file, mode ,difpat)
!!!                 if(Err_diffpatt) print*, trim(err_diffpatt_mess)
!!!
!!!           !Fin de lectura
!!!           if(allocated (resta)) deallocate(resta)
!!!           allocate(resta(difpat%npts))
!!!           !Algunas operaciones antes de empezar:
!!!           check_sym = .false.
!!!           IF(symgrpno == UNKNOWN) THEN
!!!             symgrpno = get_sym(ok)
!!!             IF(.NOT.ok) GO TO 999
!!!             WRITE(op,200) 'Diffraction point symmetry is ',pnt_grp
!!!             IF(symgrpno /= 1) THEN
!!!               WRITE(op,201) '  to within a tolerance of one part in ',  &
!!!                   nint(one / tolerance)
!!!             END IF
!!!           ELSE
!!!             check_sym = .true.
!!!             CALL chk_sym(ok)
!!!             IF(.NOT.ok) GO TO 999
!!!           END IF
!!!
!!!           n_high = INT(half*(th2_max-th2_min)/d_theta) + 1
!!!
!!!
!!!           if(allocated(difpat%ycalc) ) deallocate(difpat%ycalc)
!!!           allocate(difpat%ycalc(n_high))
!!!           filenam = trim(infile(1:(index(infile,'.')-1)))
!!!
!!!           open(unit=san_out,file=trim(filenam)//".out", status="replace",action="write")
!!!
!!!           gen(1:st%npar) = st%config(1:st%npar)
!!!           difpat%step = difpat%step * 5.0    !we enlarge the step in order to accelerate the calculation of the theoretical pattern
!!!           vector(1:c) = st%state(1:numpar)
!!!
!!!
!!!           IF(ok) ok = get_g()
!!!           call dump (infile, p_ok)
!!!           call overlp()
!!!           call nmcoor ()
!!!           rpo = 1000                         !initialization of agreement factor
!!!
!!!
!!!          !Call Set_SimAnn_StateV (sv%npar,Con,Bounds,namepar,sv%config,sv)
!!!         !  Call SimAnneal_gen(san_out, cost3)
!!! !!----    type(SimAnn_Conditions_type),intent(in out)  :: san
!!! !!----    type(State_Vector_Type),     intent(in out)  :: st
!!!           Call Simanneal_Gen(cost3,opsan,st,san_out)
!!!
!!!
!!!           write(unit=san_out,fmt="(/,a,/,f15.4)")  " => Final configuration (for file *.san).  Rp: ", rpo
!!!           write(unit=san_out,fmt="(/,a)") &
!!!                 "  NUM           Value           Name"
!!!
!!!           do i = 1, numpar
!!!                 write(unit=san_out,fmt="(i5,f15.4,tr10, a8)") i, statok(i) , namepar(i)
!!!           end do
!!!
!!!           thmin=th2_min * rad2deg
!!!           thmax=th2_max * rad2deg
!!!           ymax = maxval(difpat%y)
!!!           ymini= -0.2 * ymax
!!!           ymin = ymini - 0.5* ymax
!!!           ymax = ymax + 0.5*ymax
!!!           CALL getfnm(filenam, outfile, '.pgf', ok)
!!!           OPEN(UNIT = out, FILE = outfile, STATUS = 'new')
!!!           call DATE_AND_TIME(date, time)
!!!
!!!             WRITE(out,'(a)')      '# .PGF (WinPLOTR Graphics file) created by FullProf:'
!!!             WRITE (out, '(12a)')  '# ' , date(7:8),'-',date(5:6),'-',date(1:4), '  at ',&
!!!                                         time(1:2),':',time(3:4),':',time(5:6)
!!!             WRITE(out,'(a)')      '#'
!!!             WRITE(out,'(a)') "# X SPACE:           1  0"
!!!             WRITE(out,'(a)') "# MAIN LEGEND TEXT:  "//trim(filenam)
!!!             WRITE(out,'(a)') "# X LEGEND TEXT   : 2Theta (degrees)"
!!!             WRITE(out,'(a)') "# Y LEGEND TEXT   : Diffracted-Intensity (arb.units)"
!!!             WRITE(out,'(a,4f14.6,2i4)') "# XMIN XMAX: " ,   thmin, thmax, thmin, thmax,1,1
!!!             WRITE(out,'(a,4f14.6,2i4)') "# YMIN YMAX: " ,  ymin ,ymax,ymin,ymax,1,1
!!!             WRITE(out,'(a)') "# X AND Y GRADUATIONS:   6  8  5  5"
!!!             WRITE(out,'(a)') "# WRITE TEXT (X grad., Y grad. , Yneg. grad. , file_name):   1  1  1  1"
!!!             WRITE(out,'(a)') "# GRID (X and Y):            0  0"
!!!             WRITE(out,'(a)') "# FRAME FEATURES:           0.70    3    3    1    4    3"
!!!             WRITE(out,'(a)') "# DRAW ERROR BARRS       :  N"
!!!             WRITE(out,'(a)') "# MAIN TITLE COLOR       :  RGB(  0,  0,  0)"
!!!             WRITE(out,'(a)') "# X LEGEND COLOR         :  RGB(  0,  0,  0)"
!!!             WRITE(out,'(a)') "# Y LEGEND COLOR         :  RGB(  0,  0,  0)"
!!!             WRITE(out,'(a)') "# X GRADUATIONS COLOR    :  RGB(  0,  0,  0)"
!!!             WRITE(out,'(a)') "# Y GRADUATIONS COLOR    :  RGB(  0,  0,  0)"
!!!             WRITE(out,'(a)') "# BACKGROUND SCREEN COLOR:  RGB(240,202,166)"
!!!             WRITE(out,'(a)') "# BACKGROUND TEXT COLOR  :  RGB(255,  0,  0)"
!!!             WRITE(out,'(a)') "# BACKGROUND PLOT COLOR  :  RGB(255,255,255)"
!!!             WRITE(out,'(a)') "# PLOT FRAME COLOR       :  RGB(  0,  0,  0)"
!!!             WRITE(out,'(a)') "# NUMBER OF PATTERNS:            4          "
!!!             WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
!!!             WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',1
!!!             write(out,'(a,a)')          '#        FILE NAME  : ', " Observed "
!!!             write(out,'(a,a)')          '#            TITLE  : ', " Yobs(res) "
!!!             write(out,'(a,i10)')        '#  NUMBER OF POINTS : ',difpat%npts
!!!             write(out,'(a)')            '#            MARKER : 4'       !open circles
!!!             write(out,'(a,F6.1)')       '#              SIZE : 1.5'     !size 1
!!!             write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0,  0,255)' !Red
!!!             write(out,'(a,i6)')         '#             STYLE : 0'       !Points non continuous line
!!!             write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
!!!             write(out,'(a,i6)')         '#        DATA: X Y  '
!!!             do i=1,difpat%npts
!!!               write(unit=out,fmt="(f10.6,3f14.6)") difpat%x(i), difpat%y(i)
!!!             end do
!!!             WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
!!!             WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',2
!!!             write(out,'(a,a)')          '#        FILE NAME  : ', " Calculated "
!!!             write(out,'(a,a)')          '#            TITLE  : ', " Ycal(res) "
!!!             write(out,'(a,i10)')        '#  NUMBER OF POINTS : ', n_high
!!!             write(out,'(a)')            '#            MARKER : 4'       !open circles
!!!             write(out,'(a,F6.1)')       '#              SIZE : 0.0'     !size
!!!             write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0, 0,0)' !Red
!!!             write(out,'(a,i6)')         '#             STYLE : 1'       !Continuous line
!!!             write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
!!!             write(out,'(a,i6)')         '#        DATA: X Y  '
!!!             do i=1,n_high
!!!               theta = thmin +(i-1)*d_theta*rad2deg * 2
!!!               write(unit = out,fmt = "(f10.6,3f14.6)") theta,  ycalcdef(i)
!!!             end do
!!!             p_resta= n_high-1
!!!             WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
!!!             WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',3
!!!             write(out,'(a,a)')          '#        FILE NAME  : ', " Difference"
!!!             write(out,'(a,a)')          '#            TITLE  : ', " Yobs-Ycal "
!!!             write(out,'(a,i10)')        '#  NUMBER OF POINTS : ', p_resta
!!!             write(out,'(a)')            '#            MARKER : 4'       !open circles
!!!             write(out,'(a,F6.1)')       '#              SIZE : 0.0'     !size 1
!!!             write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  255,  0,0)' !?
!!!             write(out,'(a,i6)')         '#             STYLE : 1'       !Points non continuous line
!!!             write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
!!!                     write(out,'(a,i6)')         '#        DATA: X Y  '
!!!             thmin_o = difpat%xmin * deg2rad      !rad
!!!
!!!             if    ( th2_min/= thmin_o ) then
!!!               punts= INT(half*(th2_min-thmin_o)/d_theta)
!!!             end if
!!!
!!!             do i=1,n_high-1
!!!                  resta(i) =(difpat%y(i+punts) - ycalcdef(i)) + ymini
!!!                  write(unit=out,fmt="(2f15.7)") difpat%x(i+punts) , resta(i)
!!!             end do
!!!             WRITE(out,'(a1,128a1)')     '#',('-',i=1,128)
!!!             WRITE(out,'(a,i6)')         '# >>>>>>>> PATTERN #: ',4
!!!             write(out,'(a,a)')          '#        FILE NAME  : ', " Bragg_position "
!!!             write(out,'(a,a)')          '#            TITLE  : ', " Bragg_position "
!!!             write(out,'(a,i10)')        '#  NUMBER OF POINTS : '  , d_punt
!!!             write(out,'(a)')            '#            MARKER : 8'       !open circles
!!!             write(out,'(a,F6.1)')       '#              SIZE : 3.0'     !size 0
!!!             write(out,'(a,a16)')        '#          RGBCOLOR : RGB(  0,  128,  0) ' !black
!!!             write(out,'(a,i6)')         '#             STYLE : 0'       !Continuous line
!!!             write(out,'(a,i6)')         '#         PEN WIDTH : 1'       !current_pen_width
!!!             write(out,'(a,i6)')         '#        DATA: X Y  '
!!!
!!!                deg = ymini * 0.25
!!!                aa=h_min
!!!                Do a=h_min, h_max
!!!                   bb=k_min
!!!                  do b=k_min, k_max
!!!                     cc= zero
!!!                    do c=1, 30
!!!                      if (dos_theta(aa,bb,cc) /= zero) then
!!!                         write(unit= out,fmt = "(2F15.7, 5x, a, 3i3, a, i3)") dos_theta(aa, bb, cc)  , deg,  "(", aa,bb,cc,")", 1
!!!                      end if
!!!                      cc=cc+1
!!!                    End do
!!!                     bb=bb+1
!!!                  End do
!!!                   aa=aa+1
!!!               End do
!!!
!!!             write(out,"(a)") "# END OF FILE "
!!!             close(unit=out)
!!!
!!!


          Case (2) !SIMPLEX

              !Lectura del  pattern experimental y de los scattering factors:
              write(*,"(a)") " => Reading Pattern file="//trim(dfile)
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
              write(*,"(a)") " => Reading Background file="//trim(background_file)
              call read_background_file(background_file, mode ,difpat)
                    if(Err_diffpatt) print*, trim(err_diffpatt_mess)

              !Fin de lectura

              if(allocated (resta)) deallocate(resta)
              allocate(resta(difpat%npts))

              !Algunas operaciones antes de empezar:
              check_sym = .false.
              IF(symgrpno == UNKNOWN) THEN
                symgrpno = get_sym(ok)
                IF(.NOT. ok) GO TO 999
                WRITE(op,200) 'Diffraction point symmetry is ',pnt_grp
                IF(symgrpno /= 1) THEN
                  WRITE(op,201) '  to within a tolerance of one part in ',  &
                      nint(one / tolerance)
                END IF
              ELSE
                check_sym = .true.
                CALL chk_sym(ok)
                IF(.NOT. ok) GO TO 999
              END IF

              n_high = INT(half*(th2_max-th2_min)/d_theta) + 1


              if(allocated(difpat%ycalc) ) deallocate(difpat%ycalc)
              allocate(difpat%ycalc(n_high))
              filenam = trim(infile(1:(index(infile,'.')-1)))

             ! gen(1:n_plex) = v_plex(1:n_plex)


              IF(ok) ok = get_g()
              write(*,"(a)") " => Calling dump file="//trim(infile)
              call dump(infile, p_ok)

              call overlp()
              call nmcoor ()

              rpo = 1000                         !initialization of agreement factor
             ! rpl = 0
              do i=1, opti%npar                       !creation of the step sizes
                steplex(i) = 0.2 * crys%NP_refi(i)
                vector(i) = crys%NP_refi(i)
              end do
              open (unit=23, file='nelder_mess.out', status='replace', action='write')
              !write(*,*) "npar", opti%npar, numpar
              call  Nelder_Mead_Simplex( F_cost,opti%npar  ,crys%NP_refi(1:opti%npar) , &
                                         steplex(1:opti%npar), var_plex(1:opti%npar), rpl, opti, ipr=23)
              write(*,*)'Rp', rpo
              write(*,*) '_____________________________________'
              write(*,'(3a)') ' Parameter     refined value    '
              write(*,*) '_____________________________________'
              do i = 1, numpar
                    write(*,*)  namepar(i)  ,statok(i)
              end do


              CALL getfnm(filenam, outfile, '.prf', ok)
              if (ok) then
                OPEN(UNIT = out, FILE = outfile, STATUS = 'replace')
                call Write_Prf(difpat,out)
              else
                write(*,*) 'The outfile cannot be created'
              end if
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

          Case (3) !Local optimizer

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

             ! gen(1:n_plex) = v_plex(1:n_plex)
              IF(ok) ok = get_g()
              call dump (infile, p_ok)
              call overlp()
              write(*,*) "overlp"
              call nmcoor ()
              write(*,*) "opti%npar", opti%npar

              chi2o = 1000                         !initialization of agreement factor
            !  rpl = 0
              do i=1, opti%npar                       !creation of the step sizes
                vector(i) = crys%NP_refi(i)
                write(*,*) "crys%vlim1(i)", crys%vlim1(i)
              end do

              open (unit=23, file='local_optimizer.out', status='replace', action='write')
              if (opti%method == "DFP_NO-DERIVATIVES" .or. opti%method == "LOCAL_RANDOM" .or. opti%method == "UNIRANDOM") then
                call Lcase(opti%method )

                call Local_Optimize( F_cost,crys%NP_refi(1:opti%npar) ,  rpl, opti, mini=crys%vlim1(1:opti%npar),&
                                    maxi=crys%vlim2(1:opti%npar),ipr=23  )
              else
                write(*,*) "simplex to be added"
              end if
              write(*,*)'Rp', rpo
              write(*,*) '______________________________________'
              write(*,'(3a)') ' Parameter     refined value     '
              write(*,*) '______________________________________'
              do i = 1, numpar
                    write(*,*)  namepar(i)  ,statok(i)
              end do

              CALL getfnm(filenam, outfile, '.prf', ok)
              if (ok) then
                OPEN(UNIT = out, FILE = outfile, STATUS = 'replace')
                call Write_Prf(difpat,out)
              else
                write(*,*) 'The outfile cannot be created'
              end if
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
          Case (4) !LMQ

            !Lectura del  pattern experimental y de los scattering factors:
            write(*,"(a)") " => Reading Pattern file="//trim(dfile)
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
            write(*,"(a)") " => Reading Background file="//trim(background_file)
            call read_background_file(background_file, mode ,difpat)
                  if(Err_diffpatt) print*, trim(err_diffpatt_mess)

            !Fin de lectura

            if(allocated (resta)) deallocate(resta)
            allocate(resta(difpat%npts))

            !Algunas operaciones antes de empezar:
            check_sym = .false.
            IF(symgrpno == UNKNOWN) THEN
              symgrpno = get_sym(ok)
              IF(.NOT. ok) GO TO 999
              WRITE(op,200) 'Diffraction point symmetry is ',pnt_grp
              IF(symgrpno /= 1) THEN
                WRITE(op,201) '  to within a tolerance of one part in ',  &
                    nint(one / tolerance)
              END IF
            ELSE
              check_sym = .true.
              CALL chk_sym(ok)
              IF(.NOT. ok) GO TO 999
            END IF

            n_high = INT(half*(th2_max-th2_min)/d_theta) + 1


            if(allocated(difpat%ycalc) ) deallocate(difpat%ycalc)
            allocate(difpat%ycalc(n_high))
            filenam = trim(infile(1:(index(infile,'.')-1)))

           ! gen(1:n_plex) = v_plex(1:n_plex)


            IF(ok) ok = get_g()
            write(*,"(a)") " => Calling dump file="//trim(infile)
            call dump(infile, p_ok)

            call overlp()
            call nmcoor ()

            chi2o = 1.0E10                         !initialization of agreement factor
          !  rpl = 0
            do i=1, opti%npar                       !creation of the step sizes
              vector(i) = crys%NP_refi(i)
            end do
            open (unit=23, file='nelder_mess.out', status='replace', action='write')

            call Levenberg_Marquardt_Fit(cost_LMQ, difpat%npts, cond, Vs, chi2, infout)
            write(*,*) "infout" // trim(infout)
            !write(*,*)'Rp', rpo
            write(*,*) '_____________________________________'
            write(*,'(3a)') ' Parameter     refined value    '
            write(*,*) '_____________________________________'
            do i = 1, numpar
               write(*,*)  namepar(i)  ,statok(i)
            end do

            CALL getfnm(filenam, outfile, '.prf', ok)
              if (ok) then
                OPEN(UNIT = out, FILE = outfile, STATUS = 'replace')
                call Write_Prf(difpat,out)
              else
                write(*,*) 'The outfile cannot be created'
              end if

!-------------------------------------------------------------------------------------------------------------------------------

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

