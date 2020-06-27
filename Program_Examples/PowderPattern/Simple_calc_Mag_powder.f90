 !!----
 !!---- PROGRAM:  Simple_Calculation_of_Powder_Patterns
 !!----
 !!---- Simple program for calculating powder patterns by reading a CIF or a CFL file
 !!----
 !!---- If a CIF file is read there is no control on powder pattern parameters. A standard
 !!---- pattern corresponding to D2B with lambda=1.594 is calculated. The parameters of this
 !!---- powder pattern are described by the components of the object PPC (powder pattern constants)
 !!---- Default Powder Diffraction Pattern
 !!----    stlmax=0.6; PPC%Title="Default D2B-Powder Pattern"
 !!----    PPC%U=0.009687; PPC%V=-0.035561; PPC%W=0.049284; PPC%LAMBDA=1.594; PPC%X=0.008905
 !!----    PPC%Thmin=1.00; PPC%step=0.05;  PPC%Thmax= int(2.0*asind(stlmax*1.594)); PPC%job=1
 !!----    PPC%Ls=1900.0; nf=30; PPC%bkg=50.0
 !!----    powfile="powder_pattern.dat"
 !!----
 !!---- If a CFL is provided the directives for calculating the powder pattern must be provided
 !!---- before the list of the atoms. The following directives can be given:
 !!----
 !!----   TITLE  whatever set of characters identifying the job
 !!----   CELL   a b c alpha beta gamma (four real values, cell parameters)
 !!----   SHUB   SpgSymbol[:a,b,c;1/2,0,0]     Standard Shubnikov symbo and setting change
 !!----   UVWX   u v w x    (four real values defining the resolution of the instrument TCH Pseudo-Voigt)
 !!----   LAMBDA  lambda    (1 real value, wavelength)
 !!----   JOBTYPE   Neutrons
 !!----   PATTERN  Tmin  step Tmax  (intial angle, step and final angle for calculations of powder pattern)
 !!----   SIZE_L   Lorentzian_Size   (one real in angstroms)
 !!----   POWFILE  name_of_powder_diffraction_data_file
 !!----   BACKGD   bkg   (1 real, level of background in counts)
 !!----
 !!----   ATOM  Label  ScatSymb   x    y    z    Biso   Occ   Moment:  mx  my  mz
 !!----
 !!---- The number of ATOM directives is not limited. The items are 2 strings five reals,
 !!---- the string "Moment:" and three reals corresponding to the components of the magnetic moment
 !!----
 !!---- Remember that for a magnetic atom the scattering form-factor name should be of the
 !!---- form MXXV or JMMV, where XX is the chemical element and V is the valence
 !!----
 !!---- The program uses CrysFML and a module called Gen_Powder_Pattern where the subroutine for
 !!---- calculating the powder diffraction pattern is stored
 !!----
 Module Gen_Mag_Powder_Pattern
    !---- Use Modules ----!
    use CFML_GlobalDeps,           only: to_Deg,cp
    use CFML_Math_General,         only: asind,locate
    use CFML_Reflections_Utilities,only: Reflect_List_Type
    use CFML_Structure_Factors,    only: Strf_List_Type
    use CFML_Diffraction_Patterns, only: Diffraction_Pattern_Type, Allocate_Diffraction_Pattern
    use CFML_PowderProfiles_CW,    only: PseudoVoigt

    !---- Variables ----!
    implicit none

    private

    public  :: calc_powder_pattern, Write_PRF
    private :: TCH

    Type, public :: PowPat_CW_Conditions
       character(len=140) :: title
       integer :: job
       real    :: Lambda, U, V, W, X, Ls
       real    :: Thmin, Thmax, step
       real    :: scalef,bkg
    End Type PowPat_CW_Conditions

 Contains
    !!----
    !!---- Pure Subroutine TCH(Hg,Hl,Fwhm,Eta)
    !!----
    !!---- Calculation of eta and FWHM of the pV-function for the
    !!---- T-C-H representation.
    !!----
    !!
    Pure Subroutine TCH(Hg,Hl,Fwhm,Eta)
       !---- Arguments ----!
       real, intent(in)  :: hg
       real, intent(in)  :: hl
       real, intent(out) :: fwhm
       real, intent(out) :: eta

       !---- Variables ----!
       real, parameter :: o1= 2.69269, o2=2.42843, o3=4.47163, o4= 0.07842
       real, parameter :: e1= 1.36603, e2=0.47719, e3=0.11116
       real            :: ctl, tlr

       ! There is no exception handling because it is supposed to be
       ! perfomed before calling TCH
       ctl=hg**5.0+o1*hg**4.0*hl+o2*hg**3.0*hl**2.0+o3*hg**2.0*hl**3.0+  &
           o4*hg*hl**4.0+hl**5.0
       fwhm=ctl**0.2
       tlr = hl/fwhm
       eta = max(1.0e-06, e1*tlr-e2*tlr*tlr+e3*tlr**3.0)  !eta

       Return
    End Subroutine TCH

    Subroutine Calc_Powder_Pattern(Ppc,Hkl,Stf,Pat)
       !---- Argument ----!
       Type(PowPat_CW_Conditions),     intent(in)  :: PPC
       Type(Reflect_List_Type),        intent(in)  :: hkl
       Type(Strf_List_Type),           intent(in)  :: stf
       Type(Diffraction_Pattern_Type), intent(out) :: Pat

       !--- Local Variables ----!
       integer :: i,j,npts,i1,i2
       real    :: step,Intens,Bragg,Hl,Hg, ss,cs,tt,th1,th2,LorentzF, Y,eta,fwhm,chw

       npts=(PPC%Thmax-PPC%Thmin)/PPC%step + 1.02
       call Allocate_Diffraction_Pattern(Pat,npts)
       Pat%Title=adjustl(Trim(PPC%title))
       i=len_trim(Pat%Title)
       write(unit=Pat%Title(i+2:),fmt="(a,f7.4,f7.1,a)") " => Lambda & Lorentzian Size: ", &
                  PPC%Lambda,PPC%Ls," (angstroms)"

       Pat%scat_var="2-Theta"
       Pat%instr="Default D2B-Powder Pattern"
       Pat%xmin= PPC%Thmin
       Pat%xmax= PPC%Thmax
       Pat%ymin= 0.0
       Pat%ymax=0.0
       Pat%scal=1.0
       Pat%monitor=0.0
       Pat%step=PPC%step
       Pat%Tsamp=300.0
       Pat%Tset=300.0
       Pat%npts=npts
       Pat%ct_step=.true.
       Pat%conv=(/PPC%Lambda,PPC%Lambda,0.0,0.0,0.0/)
       chw=15.0
       do i=1,npts
          Pat%x(i)=Pat%xmin+real(i-1)*Pat%step
       end do

       Y= to_deg*PPC%Lambda/PPC%Ls
       do i=1,hkl%nref
          ss=PPC%Lambda*hkl%ref(i)%S
          cs=sqrt(abs(1.0-ss*ss))
          tt=ss/cs
          LorentzF=0.5/(ss*ss*cs)
          Bragg=2.0*asind(ss)
          HG=sqrt(tt*(PPC%U*tt+PPC%V)+PPC%W)
          HL=PPC%X*tt + Y/cs
          call TCH(hg,hl,fwhm,eta)
          Select Case(nint(eta*10.0))
             Case(:2)
                chw=25.0
             case(3:5)
                chw=45.0
             case(6:7)
                chw=60.0
             case(8:)
                chw=90.0
          End Select

          th1=Bragg-chw*fwhm
          th2=Bragg+chw*fwhm
          i1=Locate(Pat%x,npts,th1)
          i2=Locate(Pat%x,npts,th2)
          i1=max(i1,1)
          i2=min(i2,npts)
          Intens= LorentzF *hkl%ref(i)%mult * (stf%strf(i)%sqNuc+stf%strf(i)%sqMiv) * PPC%Scalef
          do j=i1,i2
             Pat%ycalc(j)=Pat%ycalc(j)+ PseudoVoigt( Pat%x(j)-Bragg, (/fwhm,eta /) ) * Intens
          end do
          Pat%ymin= minval(Pat%ycalc)
          Pat%ymax= maxval(Pat%ycalc)
       end do

       return
    End Subroutine Calc_Powder_Pattern

    Subroutine Write_PRF(fileprf,lambda,Pat,hkl)
       character(len=*),               intent(in)  :: fileprf
       real(kind=cp),                  intent(in)  :: lambda
       Type(Diffraction_Pattern_Type), intent(in)  :: Pat
       Type(Reflect_List_Type),        intent(in)  :: hkl

       integer :: i,j,i_prf
       character(len=*),parameter :: tb=char(9)
       character (len=50) :: forma1,forma2
       real :: ymax,scl
       open(newunit=i_prf,file=trim(fileprf),status="replace",action="write")

       ymax=Pat%ymax
       scl=1.0
       do
         if(ymax < 1.0e6) exit !on exit we have the appropriate value of scl
         scl=scl*0.1
         ymax=ymax*scl
       end do
       if(ymax < 100.0) then
        forma1='(f12.4,4(a,f8.4))'
       else if(ymax < 1000.0) then
        forma1='(f12.4,4(a,f8.3))'
       else if(ymax < 10000.0) then
        forma1='(f12.4,4(a,f8.2))'
       else if(ymax < 100000.0) then
        forma1='(f12.4,4(a,f8.1))'
       else
        forma1='(f12.4,4(a,f8.0))'
       end if
       write(i_prf,'(a)') trim(Pat%Title)  !//"  CELL:"
                                !  N_phases, N_points, Lamda1,Lambda2,zero,shift1,shif2,Ixunit
       write(i_prf,'(I3,I7,5f12.5,i5)')1,Pat%npts,lambda,lambda,0.0,0.0,0.0,0
       write(i_prf,'(17i6)')  hkl%Nref, 0 , 0
       write(i_prf,'(15a)')' 2Theta',tb,'Yobs',tb,'Ycal',tb,  &
        'Yobs-Ycal',tb,'Backg',tb,'Posr',tb,'(hkl)',tb,'K'
       !dd=(y(i,n_pat)-yc(i,n_pat))*scl
       do  i=1,Pat%npts
         write(i_prf,forma1) Pat%x(i),tb,Pat%ycalc(i)*scl,tb,Pat%ycalc(i)*scl,tb, -ymax/4.0,tb,0.0
       end do
       !Writing reflections
       do j=1,hkl%Nref
         !iposr=-(k-1)*ideltr
           write(i_prf,'(f12.4,9a,i8,a,3i3,a,2i3)')  2.0*asind(hkl%Ref(j)%s*Lambda), &
               tb,'        ',tb,'        ',tb,'        ',  &
               tb,'        ',tb,0, tb//'(',hkl%Ref(j)%h,')'//tb,hkl%Ref(j)%imag,1
       end do
      close(unit=i_prf)
    End Subroutine Write_PRF

  End Module Gen_Mag_Powder_Pattern

  !!----
  !!----  Program Calculation_of_Powder_Patterns
  !!----
  !!----
  !!---- Update: June - 2009
  !!
  Program Simple_Calculation_of_Powder_Patterns
     !---- Use Modules ----!
     use CFML_Math_General,              only: asind,sind
     use CFML_Atom_TypeDef,              only: Atom_List_Type, Allocate_Atom_List,Write_Atom_List
     use CFML_Crystal_Metrics,           only: Crystal_Cell_type, set_Crystal_Cell,write_crystal_cell
     use CFML_string_utilities,          only: u_case
     use CFML_Reflections_Utilities,     only: Reflect_List_Type
     Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Set_SpaceGroup, &
                                               Write_Magnetic_Space_Group,Magnetic_Space_Group_Type
     Use CFML_Structure_Factors,         only: Write_Structure_Factors, Magnetic_Structure_Factors,&
                                               ERR_SFac,ERR_SFac_Mess
     use CFML_Diffraction_Patterns,      only: Diffraction_Pattern_Type, &
                                               Allocate_Diffraction_Pattern
     use CFML_IO_Formats,                only: Readn_set_Xtal_Structure,err_form_mess,err_form, &
                                               file_list_type, Get_moment_ctr
     use CFML_Structure_Factors,         only: Strf_List_Type

     use Gen_Mag_Powder_Pattern

     !---- Variables ----!
     implicit none

     integer                :: i,j,k,l,m,n,maxnumref,ier,mult,nf,codini=0
     integer                :: lun=1,lp=2
     real, dimension(3)     :: ad,ang,x,fr,codes=(/11.0,21.0,31.0/)
     character(len=1)       :: ans,outa
     integer, dimension(3)  :: ncel,h
     real                   :: stlmax,tini,tfin,sn,sf2,tim,ftim=1.0,box
     character(len=132)     :: line,powfile,filcod,prf_file
     character(len=3)       :: mode
     character(len=8)       :: units="seconds",radiation
     character(len=4),dimension(:),allocatable :: ch
     logical                :: full=.true.

     Type(Crystal_Cell_type)        :: cell
     Type(Magnetic_Space_Group_Type):: SpG
     Type(Atom_List_Type)           :: A
     Type(Reflect_List_Type)        :: hkl
     Type(Strf_List_Type)           :: stf
     Type(Diffraction_Pattern_Type) :: Pat
     Type(PowPat_CW_Conditions)     :: PPC
     Type(file_list_type)           :: fich_cfl
     integer                        :: narg
     Logical                        :: esta, arggiven=.false., fail


     !---- Arguments on the command line ----!
     narg=command_argument_count()

     if (narg > 0) then
        call get_command_argument(1,filcod)
        arggiven=.true.
     end if


     write(unit=*,fmt="(/,/,6(a,/))")                                                     &
          "            ------ PROGRAM SIMPLE POWDER PATTERN CALCULATION  ------"        , &
          "                    ---- Version 0.2 January-2020----"                         , &
          "    **********************************************************************"  , &
          "    * Calculates powder diffraction pattern from a *.CFL or a *.CIF file *"  , &
          "    **********************************************************************"  , &
          "                          (JRC- January-2020 )"
     write(unit=*,fmt=*) " "

     if (.not. arggiven) then
        write(unit=*,fmt="(a)", advance='no') " => Code of the file xx.cif(cfl) (give xx): "
        read(unit=*,fmt="(a)") filcod
        if(len_trim(filcod) == 0) stop
     end if


     inquire(file=trim(filcod)//".mcif",exist=esta)
     if (esta) then
        call Readn_set_Xtal_Structure(trim(filcod)//".mcif",Cell,SpG,A,Mode="CIF")
     else
        inquire(file=trim(filcod)//".cfl",exist=esta)
        if ( .not. esta) then
           write(unit=*,fmt="(a)") " File: "//trim(filcod)//".mcif (or .cfl) does'nt exist!"
           stop
        end if
        call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)
        mode="CFL"
     end if

     if (err_form) then
        write(unit=*,fmt="(a)") trim(err_form_mess)
     else
       open(unit=lun,file=trim(filcod)//".powder", status="replace",action="write")
       write(unit=lun,fmt="(/,/,6(a,/))")                                                 &
          "            ------ PROGRAM SIMPLE POWDER PATTERN CALCULATION  ------"        , &
          "                    ---- Version 0.2 January-2020----"                         , &
          "    **********************************************************************"  , &
          "    * Calculates powder diffraction pattern from a *.CFL or a *.CIF file *"  , &
          "    **********************************************************************"  , &
          "                          (JRC- January-2020 )"
       ! Calculate a default Powder Diffraction Pattern
       stlmax=0.6; PPC%Title="Powder Pattern of Structure provided in: "//trim(filcod)//".cfl"
       PPC%U=0.0002; PPC%V=-0.0002; PPC%W=0.012; PPC%LAMBDA=1.54056; PPC%X=0.0015
       PPC%Thmin=1.00; PPC%step=0.05;  PPC%Thmax= int(2.0*asind(stlmax*1.54056)); PPC%job=0
       PPC%Ls=1900.0;  nf=30; PPC%bkg=50.0
       powfile=trim(filcod)//".xys"
       units=" seconds"
       tim=0.0

       ! Write initial structure information in the .powder file
       call Write_Crystal_Cell(Cell,lun)
       call Write_Magnetic_Space_Group(SpG,lun,full)
       !Get information on moment constraints and modify the list of atoms accordingly
       write(unit=lun,fmt="(/,a,/)") " => Symmetry constraints in magnetic moments:"
       do i=1,A%natoms
         if(A%Atom(i)%moment < 0.001) cycle !Skip non-magnetic atoms
         call Get_moment_ctr(A%Atom(i)%X,A%Atom(i)%M_xyz,Spg,codini,codes)
         write(unit=lun,fmt="(a12,3(a,3f10.4))") "     "//A%Atom(i)%Lab," Pos:",A%Atom(i)%X," Mom:",A%Atom(i)%M_xyz," Codes:",codes
       end do
       call Write_Atom_List(A,level=2,lun=lun)

       ! Look for calculation conditions in the CFL file that are provided before the list of atoms
       if (mode == "CFL") then
          do i=1,fich_cfl%nlines
             line=adjustl(u_case(fich_cfl%line(i)))
             !if(line(1:4) == "ATOM") exit

             Select Case(Trim(line(1:7)))
                Case("TITLE")
                   PPC%title= adjustl(fich_cfl%line(i)(6:))

                Case("UVWX")
                   read(unit=line(7:),fmt=*,iostat=ier) PPC%U,PPC%V,PPC%W,PPC%X
                   if (ier /= 0) then
                      PPC%U=0.02; PPC%V=-0.02; PPC%W=0.12; PPC%X=0.0015
                   end if

                Case("LAMBDA")
                   read(unit=line(7:),fmt=*,iostat=ier) PPC%LAMBDA
                   if(ier /= 0) PPC%LAMBDA=1.56

                Case("BACKGD")
                   read(unit=line(7:),fmt=*,iostat=ier) PPC%bkg
                   if(ier /= 0) PPC%bkg=20.0

                Case("JOBTYPE")
                   radiation = adjustl(fich_cfl%line(i)(8:))
                   if (radiation(1:1) == "N" .or. radiation(1:1) == "n") then
                      PPC%job=1
                      write(*,*) " => Neutrons Job ...."
                   else
                      PPC%job=0
                      write(*,*) " => X-rays Job ...."
                   end if

                Case("PATTERN")
                   read(unit=line(8:),fmt=*,iostat=ier) PPC%Thmin, PPC%step,  PPC%Thmax
                   if (ier /= 0) then
                      PPC%Thmin=1.00; PPC%step=0.05;  PPC%Thmax= 120.0
                   end if

                Case("SIZE_L")
                   read(unit=line(8:),fmt=*,iostat=ier) PPC%Ls
                   if (ier /= 0) then
                      PPC%Ls=1900.0
                   end if
                   write(*,*) " => Lorentzian size: ", PPC%Ls

                Case("POWFILE")
                   powfile=adjustl(fich_cfl%line(i)(8:))
                   if (len_trim(powfile) == 0) then
                      powfile="powder_pattern.xys"
                   end if
                   write(*,*) " => Powder pattern file: "//trim(powfile)

             End Select
          end do

       End if

       ! Calculate sinTheta/Lambda max from 2Thetamax     PPC%Thmax= int(2.0*asind(stlmax*1.56))
       stlmax=sind(min((PPC%Thmax+10.0)*0.5,90.0))/PPC%lambda

     end if !if error

     write(unit=lun,fmt="(/,a)")  " => CALCULATION OF NEUTRON POWDER DIFFRACTION PATTERN"
     !PPC%title=Trim(PPC%title)
     write(unit=lun,fmt="(  a,4f10.5)")  " => Resolution parameters UVWX: ",PPC%U,PPC%V,PPC%W,PPC%X
     write(unit=lun,fmt="(  a, f10.5)")  " => Lambda: ",PPC%lambda
     write(unit=lun,fmt="(  a, f10.5,a)")" => Background level: ",PPC%bkg," counts"
     write(unit=lun,fmt="(  a,2f10.2)")  " => Lorentzian size: ",PPC%Ls
     write(unit=lun,fmt="(  a,3f10.5)")  " => 2Theta range and step: ",PPC%Thmin,PPC%step,PPC%Thmax
     write(unit=lun,fmt="(  a,3f10.5)")  " => Maximum sin(Theta)/Lambda (for generating reflections): ",stlmax

     ! Now calculate structure factors
     write(*,*) " => Calculating structure factors ..."
     call cpu_time(tini)
     call Magnetic_Structure_Factors(Cell,A,SpG,stlmax,hkl,Stf,lun)
     if (ERR_SFac) then
         write(*,*) " => Error in calculations of Structure Factors"
         write(*,*) " => "//trim(ERR_SFac_Mess)
         stop
     end if
     call cpu_time(tfin)
     tim=tim+ tfin-tini
     write(unit=*,fmt="(a,f15.3,a)")   "  => CPU-time used for Structure_Factors: ",(tfin-tini)*ftim,units
     write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time used for Structure_Factors: ",(tfin-tini)*ftim,units
     call Write_Structure_Factors(lun,hkl,stf,full)
     write(*,"(a,i8)") "  => Total number of generated reflections is ",hkl%nref
     write(unit=lun,fmt="(a,i9)") " => Total number of generated reflections is ",hkl%nref

     call cpu_time(tini)
     PPC%Scalef=cell%RCellVol
     call Calc_powder_pattern(PPC,hkl,Stf,Pat)
     call cpu_time(tfin)
     tim=tim+ tfin-tini
     write(*,*) " => CPU-time used for Calc_powder_pattern: ",(tfin-tini)*ftim,units
     write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time used for Calc_powder_pattern: ",(tfin-tini)*ftim,units
     write(*,*) " => CPU-time for all calculations: ",tim*ftim,units
     write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time for all calculations: ",tim*ftim,units
     write(*,*) " => Writing powder pattern file: "//trim(powfile)
     open(unit=lp,file=trim(powfile),status="replace",action="write")
       write(unit=lp,fmt="(a)") "XYDATA"
       write(unit=lp,fmt="(a,f12.5)") "TITLE "//trim(Pat%Title)
       write(unit=lp,fmt="(a,5f12.5)") "UVWX_Size_L: ", PPC%U, PPC%V, PPC%W,PPC%X,PPC%Ls
       write(unit=lp,fmt="(a)") "FILE: "//trim(powfile)
       write(unit=lp,fmt="(a)") "TEMP    273.0   273.0"
       write(unit=lp,fmt="(a)") "INTER   1.0000  1.0000  0 0.00000 <- internal multipliers for X, Y-Sigma, Interpol, StepIn"
       write(unit=lp,fmt="(a,3f10.4)")"! Angular range (min,step,max):", Pat%xmin,Pat%step,Pat%xmax
       write(unit=lp,fmt="(a)")"!     2theta             Y       "
       do i=1,Pat%npts
         write(unit=lp,fmt="(2f16.4)") Pat%x(i),Pat%ycalc(i)+PPC%bkg
       end do
     close(unit=lun)
     if(len_trim(powfile) == 0) then
       prf_file=trim(filcod)//".prf"
     else
       i=index(powfile,".",back=.true.)
       prf_file=powfile(1:i)//"prf"
     end if
     write(*,*) " => Writing PRF file: "//trim(prf_file)
     call Write_PRF(prf_file,PPc%Lambda,Pat,hkl)
     stop

  End Program Simple_Calculation_of_Powder_Patterns

