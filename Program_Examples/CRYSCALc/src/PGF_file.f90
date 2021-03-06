!     Last change:  TR   24 Jan 2007    3:04 pm
!subroutine create_PGF_file(X,Y, h,k,l, npts, X_legend )
subroutine create_PGF_file(pgf_file, X,Y, string, npts, X_legend )
 !USE cryscalc_module, only : debug_proc
 use atome_module
 implicit none
  character (len=*), intent(in)          :: pgf_file
  REAL,    INTENT(IN), DIMENSION(npts)   :: X, Y
  !integer, intent(in), dimension(npts)   :: h,k,l
  CHARACTER(LEN=*), DIMENSION(npts)     :: string
  INTEGER, INTENT(IN)                    :: npts
  CHARACTER(LEN=*), INTENT(IN)           :: X_legend
  !local variables:
  INTEGER                                :: i, i1, long
  INTEGER                                :: PGF_style
  character (len=64)                     :: titre
  REAL                                   :: Xmin, Xmax, Ymin, Ymax
  LOGICAL                                :: FCF_plot, FCF_plot_stl

  !if(debug_proc%level_2)  call write_debug_proc_level(2, "create_PGF_file")

 ! creation d'un fichier format PGF (pour WinPLOTR): F2 = f(sinTheta/wave)

 FCF_plot     = .false.
 FCF_plot_stl = .false.
 long = len_trim(X_legend)

 ! entete

 open(unit=11, file=trim(pgf_file))
 WRITE(11, '(a)') '# .PGF (winPLOTR Graphics file) created by CRYSCALC:'
 WRITE(11, '(a)') '#'

 IF(long ==8) then
  if(X_legend(1:8) == 'sinTheta') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: F2 = f(sinTheta/lambda)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : sin(Theta)/lambda (A-1)'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : F2'
   titre =  'F2_=_f(sinTheta/lambda)'
   PGF_style = 0
  END IF
 ELSEIF(long ==5) then
  IF(X_legend(1:5) == 'd_hkl') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: F2 = f(d_hkl)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : d_hkl(A)'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : F2'
   titre =  'F2_=_f(d_hkl)'
   PGF_style = 0
  END IF
 ELSEIF(long ==5) THEN
  IF(X_legend(1:5) == 'theta') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: F2 = f(Theta)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Theta(deg)'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : F2'
   titre =  'F2_=_f(Theta)'
   PGF_style = 0
  END IF
 ELSEIF(long ==3) then
  if(X_legend(1:3) == 'fcf') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: F2c = f(F2o)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : F2_obs'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : F2_calc'
   titre =  'F2_c = f(F2_o)'
   PGF_style = 0
   FCF_plot  = .true.
  end if
 ELSEIF(long ==7) then
  if(X_legend(1:7) == 'fcf_stl') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: F2c - F2o = f(sinTheta/lambda)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : sinTheta/lambda'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : F2_calc - F2_obs'
   titre =  'F2_c - F2_o = f(SinTheta/lambda)'
   PGF_style = 0
   FCF_plot_stl  = .true.
  end if


 ELSEIF(long ==4) then
  if(X_legend(1:4) == 'bcoh') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Neutron coherent scattering length'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : bcoh (10-12 cm)'
   titre =  'bcoh_=_f(element)'
   PGF_style = 0
  end if
 ELSEIF(long ==6) then
  if(X_legend(1:6) == 'sedinc') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Neutron incoherent cross section'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : sed_inc (barns)'
   titre =  'sed_inc_=_f(element)'
   PGF_style = 1
  end if
 ELSEIF(long ==3) then
  if(X_legend(1:3) == 'sea') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Neutron absorption cross section'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : sea (barns)'
   titre =  'sea_=_f(element)'
   PGF_style = 1
  end if
 ELSEIF(long ==11) then
  if(X_legend(1:11) == 'RE_data_mod') then
   WRITE(11,'(a)')      '# DATA FROM:  Atomic data and nuclear data tables 44, 191-207 (1990)'
   WRITE(11,'(a)')      '#             Resonance effects in neutron scattering lengths or rare-earth nuclides'
   WRITE(11,'(a)')      '#             J.E. Lynn and P.A. Seeger, Los Alamos National Laboratory'
   WRITE(11,'(2a)')     '# MAIN LEGEND TEXT: Neutron coherent scattering length for ', trim(current_RE_label)
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Energy (ev)'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : |b| (10^-12 cm)'
   titre =  'mod_b = f(energy)'
   PGF_style = 1
  END IF
 ELSEIF(long ==10) then
  if(X_legend(1:10) == 'RE_data_Re') then
   WRITE(11,'(a)')      '# DATA FROM:  Atomic data and nuclear data tables 44, 191-207 (1990)'
   WRITE(11,'(a)')      '#             Resonance effects in neutron scattering lengths or rare-earth nuclides'
   WRITE(11,'(a)')      '#             J.E. Lynn and P.A. Seeger, Los Alamos National Laboratory'
   WRITE(11,'(2a)')     '# MAIN LEGEND TEXT: Neutron coherent scattering length for ', trim(current_RE_label)
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Energy (ev)'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : Re(b) (10^-12 cm)'
   titre =  'Re_b = f(energy)'
   PGF_style = 1
 elseif(X_legend(1:10) == 'RE_data_Im') then
   WRITE(11,'(a)')      '# DATA FROM:  Atomic data and nuclear data tables 44, 191-207 (1990)'
   WRITE(11,'(a)')      '#             Resonance effects in neutron scattering lengths or rare-earth nuclides'
   WRITE(11,'(a)')      '#             J.E. Lynn and P.A. Seeger, Los Alamos National Laboratory'
   WRITE(11,'(2a)')     '# MAIN LEGEND TEXT: Neutron coherent scattering length for ', trim(current_RE_label)
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Energy (ev)'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : Im(b) (10^-12 cm)'
   titre =  'Im_b = f(energy)'
   PGF_style = 1
  end if
 ELSEIF(X_legend(1:2) == 'Ag') then
  IF(X_legend(3:) == '_tics') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Total interaction cross section (Ag radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : barns'
   titre =  'TICS_X_Ag = f(element)'
   PGF_style = 1
  ELSEIF(long ==3 .and. X_legend(3:) == '_cam') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Mass attenuation coefficient (Ag radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : cm2/g'
   titre =  'MAC_X_Ag = f(element)'
   PGF_style = 1
  endif

 ELSEIF(X_legend(1:2) == 'Mo') then
  IF(X_legend(3:) == '_tics') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Total interaction cross section (Mo radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : barns'
   titre =  'TICS_X_Mo = f(element)'
   PGF_style = 1
  ELSEIF(long ==3 .and. X_legend(3:) == '_cam') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Mass attenuation coefficient (Mo radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : cm2/g'
   titre =  'MAC_X_Mo = f(element)'
   PGF_style = 1
  endif

 ELSEIF(X_legend(1:2) == 'Cu') then
  PGF_style = 1
  IF(X_legend(3:) == '_tics') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Total interaction cross section (Cu radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : barns'
   titre =  'TICS_X_Cu = f(element)'
  ELSEIF(long ==3 .and. X_legend(3:) == '_cam') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Mass attenuation coefficient (Cu radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : cm2/g'
   titre =  'MAC_X_Cu = f(element)'
  endif

 ELSEIF(X_legend(1:2) == 'Ni') then
  PGF_style = 1
  IF(long > 3 .and. X_legend(3:) == '_tics') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Total interaction cross section (Ni radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : barns'
   titre =  'TICS_X_Ni = f(element)'
  ELSEIF(long >3 .and. X_legend(3:) == '_cam') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Mass attenuation coefficient (Ni radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : cm2/g'
   titre =  'MAC_X_Ni = f(element)'
  endif

 ELSEIF(X_legend(1:2) == 'Co') then
  PGF_style = 1
  IF(X_legend(3:) == '_tics') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Total interaction cross section (Co radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : barns'
   titre =  'TICS_X_Co = f(element)'
  ELSEIF(X_legend(3:) == '_cam') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Mass attenuation coefficient (Co radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : cm2/g'
   titre =  'MAC_X_Co = f(element)'
  endif

 ELSEIF(X_legend(1:2) == 'Fe') then
  PGF_style = 1
  IF(X_legend(3:) == '_tics') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Total interaction cross section (Fe radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : barns'
   titre =  'TICS_X_Fe = f(element)'
  ELSEIF(X_legend(3:) == '_cam') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Mass attenuation coefficient (Fe radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : cm2/g'
   titre =  'MAC_X_Fe = f(element)'
  endif

 ELSEIF(X_legend(1:2) == 'Cr') then
  PGF_style = 1
  IF(X_legend(3:) == '_tics') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Total interaction cross section (Cr radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : barns'
   titre =  'TICS_X_Cr = f(element)'
  ELSEIF(X_legend(3:) == '_cam') then
   WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Mass attenuation coefficient (Cr radiation)'
   WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
   WRITE(11,'(a)')      '# Y LEGEND TEXT   : cm2/g'
   titre =  'MAC_X_Cr = f(element)'
  endif

 ELSEIF(long ==8 .and. X_legend(1:8) == 'density') then
  WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Atomic density'
  WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
  WRITE(11,'(a)')      '# Y LEGEND TEXT   : density'
  titre =  'atomic_density = f(element)'
  PGF_style = 1

 ELSEIF(long ==6 .and. X_legend(1:6) == 'radius') then
  WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Atomic radius'
  WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
  WRITE(11,'(a)')      '# Y LEGEND TEXT   : radius (A)'
  titre =  'atomic_radius = f(element)'
  PGF_style = 1

 ELSEIF(long ==6 .and. X_legend(1:6) == 'weight') then
  WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Atomic weight'
  WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
  WRITE(11,'(a)')      '# Y LEGEND TEXT   : weight'
  titre =  'atomic_weight = f(element)'
  PGF_style = 1

 ELSEIF(long ==2 .and. X_legend(1:2) == 'f0') then
  i1 = INDEX(X_legend, ':')
  WRITE(11,'(a)')      '# MAIN LEGEND TEXT: X ray scattering factor for '//X_legend(i1+1:)
  WRITE(11,'(a)')      '# X LEGEND TEXT   : sin(Theta)/lambda (A-1)'
  WRITE(11,'(a)')      '# Y LEGEND TEXT   : f0'
  WRITE(titre,'(a)')   TRIM(X_legend(i1+1:))
  PGF_style = 1

 endif


 Xmin = minval(X(1:npts))
 Xmax = MAXVAL(X(1:npts))

 Ymin = minval(Y(1:npts))
 Ymax = maxval(Y(1:npts))


 !write(11, '(a,4(x,F4.2),a)') '#   XMIN XMAX     : ', 0., Xmax, 0., Xmax, ' 1 1'
 if(.not. FCF_plot) then
 write(11, '(a,2(1x,F6.2),a)') '#   XMIN XMAX     : ', 0., real(int((10*Xmax+1)))/10
 else
 write(11, '(a,2(1x,F15.5))')  '#   YMIN YMAX     : ', Xmin, Xmax
 endif
 write(11, '(a,2(1x,F15.5))')  '#   YMIN YMAX     : ', Ymin, Ymax

 WRITE(11,'(a,i6)')     '# NUMBER OF PATTERNS: 1'
 WRITE(11,'(a1,70a1)')  '#',('-',i=1,70)
 WRITE(11,'(a,i6)')     '# >>>>>>>> PATTERN #: 1'
 write(11,'(2a)')       '#             TITLE : ',trim(titre)
 write(11,'(a,i8)')     '#  NUMBER OF POINTS : ', npts
 write(11,'(a)')        '#            MARKER : 4'
 write(11,'(a)')        '#              SIZE : 1.5'
 write(11,'(a,I1)')     '#             STYLE : ', PGF_style
 write(11,'(a)')        '#   DATA: X Y COMM'

 if(.not. FCF_plot) then
  do i=1, npts
   write(11,'(F10.6,1x,F15.5,5x,a)') X(i), Y(i), string(i)
  END do
 else
  do i=1, npts
   write(11,'(F15.5,1x,F15.5,5x,a)') X(i), Y(i), string(i)
  END do
 end if
 WRITE(11,'(a)') '# END OF FILE'


 CLOSE(UNIT=11)
 return

end subroutine create_PGF_file
!--------------------------------------------------------------------------------------------------

subroutine create_PGF_file_multi(PGF_unit, X,Y, string, npts)
 implicit none
  INTEGER, INTENT(IN)                    :: PGF_unit
  REAL,    INTENT(IN), DIMENSION(npts)   :: X, Y
  CHARACTER(LEN=*),    DIMENSION(npts)   :: string
  INTEGER, INTENT(IN)                    :: npts
  !local variables:
  INTEGER                                :: i

  do i=1, npts
   write(PGF_unit,'(F10.6,1x,F15.5,5x,a)') X(i), Y(i), string(i)
  END do

 RETURN
end subroutine create_PGF_file_multi
