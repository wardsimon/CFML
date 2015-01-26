 !!----
 !!---- PROGRAM:  Calculation_of_Powder_Patterns
 !!----
 !!---- Program for calculating powder patterns by reading a CIF or a CFL file (extension of a simpler
 !!---- program Simple_Calculation_of_Powder_Patterns
 !!---- The program can be used for generating a super-cell and perturb the atom positions
 !!---- to see the effect on the powder diffraction pattern.
 !!----
 !!---- If a CIF file is read there is no control on powder pattern parameters. A standard
 !!---- pattern corresponding to CuKa X-ray radiation is calculated. The parameters of this
 !!---- powder pattern are described by the components of the object PPC (powder pattern constants)
 !!---- Default Powder Diffraction Pattern
 !!----    stlmax=0.6; PPC%Title="Default Powder Pattern"
 !!----    PPC%U=0.0002; PPC%V=-0.0002; PPC%W=0.012; PPC%LAMBDA=1.54056; PPC%X=0.0015
 !!----    PPC%Thmin=1.00; PPC%step=0.05;  PPC%Thmax= int(2.0*asind(stlmax*1.54056)); PPC%job=0
 !!----    PPC%Ls=1900.0; PPC%Gs=  1900.0; nf=30; PPC%bkg=50.0
 !!----    powfile="powder_pattern.dat"
 !!----
 !!---- If a CFL is provided the directives for calculating the powder pattern must be provided
 !!---- before the list of the atoms. The following directives can be given:
 !!----
 !!----   TITLE  whatever set of characters identifying the job
 !!----   CELL   a b c alpha beta gamma (four real values, cell parameters)
 !!----   SPGR   SpgSymbol    (Symbol (HM/Hall) or number of the space group)
 !!----   UVWX   u v w x    (four real values defining the resolution of the instrument TCH Pseudo-Voigt)
 !!----   LAMBDA  lambda    (1 real value, wavelength)
 !!----   JOBTYPE   X-rays  (or Neutrons, if the first character is not N or n, "X-rays" is the default)
 !!----   REFLIST   (if given, the program outputs the list of structure factors)
 !!----   PATTERN  Tmin  step Tmax  (intial angle, step and final angle for calculations of powder pattern)
 !!----   SIZE_LG  Lorentzian_Size  Gaussian_Size  (two reals in angstroms)
 !!----   POWFILE  name_of_powder_diffraction_data_file
 !!----   BACKGD   bkg   (1 real, level of background in counts)
 !!----   SUPCEL    n1 n2 n3   outa   (SuperCell calculatio, 3 integer plus 1 character (outa=y or outa=n)
 !!----             If outa is "y" the program outputs the list of all generated atoms
 !!----   RANDOMB  box  (The atoms are randomly generating in a box around the ideal positions of "box"
 !!----                  angstroms side)
 !!----
 !!----   ATOM  Label  ChemSymb   x    y    z    Biso   Occ
 !!----
 !!---- The number of ATOM directives is not limited. The items are 2 strings and five reals.
 !!----
 !!---- The program uses CrysFML and a module called Atoms_Generation where the subroutine for
 !!---- calculating the atoms and generate the powder diffraction pattern are stored
 !!----
 Module Gen_Powder_Pattern_and_Atoms
    !---- Use Modules ----!
    use CFML_GlobalDeps,           only: to_Deg
    use CFML_Math_General,         only: asind,locate
    use CFML_Atom_TypeDef,         only: Atom_Type, Atom_List_Type,Allocate_Atom_List
    use CFML_Reflections_Utilities,only: Reflection_List_Type
    use CFML_Diffraction_Patterns, only: Diffraction_Pattern_Type, Allocate_Diffraction_Pattern
    use CFML_PowderProfiles_CW,    only: TCH_pVoigt,PseudoVoigt

    !---- Variables ----!
    implicit none

    private

    public :: gen_powder_pattern, gen_atoms

    Type, public :: PowPat_CW_Conditions
       character(len=140) :: title
       integer :: job
       real    :: Lambda, U, V, W, X, Ls, Gs
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

    !!----
    !!---- Subroutine Gen_Powder_Pattern(Ppc,Hkl,Pat)
    !!----
    !!
    Subroutine Gen_Powder_Pattern(Ppc,Hkl,Pat)
       !---- Argument ----!
       Type(PowPat_CW_Conditions),     intent(in)  :: PPC
       Type(Reflection_List_Type),     intent(in)  :: hkl
       Type(Diffraction_Pattern_Type), intent(out) :: Pat

       !--- Local Variables ----!
       integer :: i,j,npts,i1,i2
       real    :: step,Intens,Bragg,Hg,Hl, ss,cs,tt,th1,th2,LorentzF, IG,Y,eta,fwhm,chw

       npts=(PPC%Thmax-PPC%Thmin)/PPC%step + 1.02
       call Allocate_Diffraction_Pattern(Pat,npts)
       Pat%Title=adjustl(Trim(PPC%title))
       i=len_trim(Pat%Title)
       write(unit=Pat%Title(i+2:),fmt="(a,f7.4,2f7.1)") " => lambda,Ls,Gs: ", &
                                                        PPC%Lambda,PPC%Ls,PPC%Gs

       Pat%scat_var="2-Theta"
       Pat%instr="Calculated Pattern"
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

       IG=to_deg*PPC%Lambda/PPC%Gs
       Y= to_deg*PPC%Lambda/PPC%Ls
       do i=1,hkl%nref
          ss=PPC%Lambda*hkl%ref(i)%S
          cs=sqrt(abs(1.0-ss*ss))
          tt=ss/cs
          LorentzF=0.5/(ss*ss*cs)
          Bragg=2.0*asind(ss)
          HG=sqrt(tt*(PPC%U*tt+PPC%V)+PPC%W  + IG*IG/(cs*cs))
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
          Intens= LorentzF *hkl%ref(i)%mult * hkl%ref(i)%Fc**2 * PPC%Scalef
          do j=i1,i2
             Pat%ycalc(j)=Pat%ycalc(j)+ PseudoVoigt( Pat%x(j)-Bragg, (/fwhm,eta /) ) * Intens
          end do
       end do

       return
    End Subroutine Gen_Powder_Pattern

    !!----
    !!---- Subroutine Gen_Atoms(Chem,Xyz,Biso,Nat,Ncells,Fr,Atoms)
    !!----
    !!----
    !!---- Update: Apr - 2009
    !!
    Subroutine Gen_Atoms(Chem,Xyz,Biso,Nat,Ncells,Fr,Atoms)
       !---- Argument ----!
       character(len=*),dimension(:), intent(in) :: chem
       real,          dimension(:,:), intent(in) :: xyz
       real,          dimension(  :), intent(in) :: Biso
       integer,                       intent(in) :: nat
       integer,       dimension(3),   intent(in) :: ncells
       real,          dimension(3),   intent(in) :: fr
       Type(Atom_List_Type),         intent(out) :: atoms

       !---- Local variables ----!
       integer                                   :: i,j,k,m,n,natoms
       real, dimension(3)                        :: tr,u,ton
       real,          dimension(:,:),allocatable :: pos
       integer,       dimension(  :),allocatable :: ptr

       ! Here we generate all atoms
       natoms=nat*ncells(1)*ncells(2)*ncells(3)
       call Allocate_Atom_List(natoms,Atoms)
       ton=real(ncells)
       ton=1.0/ton
       call RANDOM_SEED()
       m=0
       do i=0,ncells(1)-1
          do j=0,ncells(2)-1
             do k=0,ncells(3)-1
                do n=1,nat
                   m=m+1
                   tr=(/real(i),real(j),real(k)/)
                   call random_number(u)
                   u=(2.0*u-1.0)*fr !box

                   write(unit=Atoms%atom(m)%Lab,fmt="(a,3i2.2)") trim(chem(n))//"_",i,j,k
                   Atoms%atom(m)%ChemSymb = chem(n)
                   Atoms%atom(m)%SfacSymb = chem(n)
                   Atoms%atom(m)%X        = (tr(:) + xyz(:,n) + u(:))*ton(:)
                   Atoms%atom(m)%Occ      = 1.0
                   Atoms%atom(m)%Biso     = Biso(n)
                end do
             end do
          end do
       end do
       Atoms%natoms=m

       return
    End Subroutine Gen_Atoms

  End Module Gen_Powder_Pattern_and_Atoms

  !!----
  !!----  Program Gen_Powder
  !!----
  !!----
  !!---- Update: April - 2009
  !!
  Program Calculation_of_Powder_Patterns
     !---- Use Modules ----!
     !use f2kcli
     use CFML_Math_General,              only: asind,sind
     use CFML_Atom_TypeDef,              only: Atom_List_Type, Allocate_Atom_List,Write_Atom_List
     use CFML_Crystal_Metrics,           only: Crystal_Cell_type, set_Crystal_Cell,write_crystal_cell
     use CFML_string_utilities,          only: u_case
     use CFML_Reflections_Utilities,     only: Reflection_List_Type,Hkl_uni,get_maxnumref
     Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Set_SpaceGroup, Get_Orbit, &
                                               Write_SpaceGroup
     Use CFML_Structure_Factors,         only: Write_Structure_Factors, Calc_hkl_StrFactor,&
                                               Init_Calc_hkl_StrFactors,ERR_SFac,ERR_SFac_Mess
     use CFML_Diffraction_Patterns,      only: Diffraction_Pattern_Type, &
                                               Allocate_Diffraction_Pattern
     use CFML_IO_Formats,                only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type

     use Gen_Powder_Pattern_and_Atoms

     !---- Variables ----!
     implicit none

     integer                :: i,j,k,l,m,n,maxnumref,ier,mult,nf
     integer                :: n1,n2,n3,nat,natoms,lun=1,lp=2
     real, dimension(3)     :: ad,ang,x,fr
     character(len=1)       :: ans,outa
     integer, dimension(3)  :: ncel,h
     real, dimension(3,192) :: orb
     real                   :: stlmax,tini,tfin,sn,sf2,tim,ftim,box
     character(len=132)     :: line,powfile,filcod
     character(len=3)       :: mode
     character(len=8)       :: units,radiation
     character(len=35)      :: time
     character(len=4),dimension(:),allocatable :: ch
     real, dimension(:,:),allocatable  :: xyz
     real, dimension(:),  allocatable  :: biso

     Type(Crystal_Cell_type)        :: cell
     Type(Space_Group_Type)         :: SpG
     Type(Atom_List_Type)           :: Atm,A
     Type(Reflection_List_Type)     :: hkl
     Type(Diffraction_Pattern_Type) :: Pat
     Type(PowPat_CW_Conditions)     :: PPC
     Type(file_list_type)           :: fich_cfl
     integer                        :: narg
     Logical                        :: esta, arggiven=.false., supcel_given=.false., out_atoms=.false.,&
                                       reflist=.false.,fail

     !---- Arguments on the command line ----!
     narg=command_argument_count()

     if (narg > 0) then
        call get_command_argument(1,filcod)
        arggiven=.true.
        i=index(filcod,".",back=.true.)
        if(i /= 0) filcod=filcod(1:i-1)
     end if


     write(unit=*,fmt="(/,/,6(a,/))")                                                     &
          "              ------ PROGRAM POWDER PATTERN CALCULATION  ------"             , &
          "                     ---- Version 0.1 April-2009----"                        , &
          "    **********************************************************************"  , &
          "    * Calculates powder diffraction pattern from a *.CFL or a *.CIF file *"  , &
          "    **********************************************************************"  , &
          "                          (JRC- April 2009 )"
     write(unit=*,fmt=*) " "

     if (.not. arggiven) then
        write(unit=*,fmt="(a)", advance='no') " => Code of the file xx.cif(cfl) (give xx): "
        read(unit=*,fmt="(a)") filcod
        if(len_trim(filcod) == 0) stop
        i=index(filcod,".",back=.true.)
        if(i /= 0) filcod=filcod(1:i-1)
     end if


     inquire(file=trim(filcod)//".cif",exist=esta)
     if (esta) then
        Mode="CIF"
        call Readn_set_Xtal_Structure(trim(filcod)//".cif",Cell,SpG,A,Mode="CIF")
     else
        inquire(file=trim(filcod)//".cfl",exist=esta)
        if ( .not. esta) then
           write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cif (or .cfl) does'nt exist!"
           stop
        end if
        Mode="CFL"
        call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)
     end if

    if (err_form) then
       write(unit=*,fmt="(a)") trim(err_form_mess)
    else
       open(unit=lun,file=trim(filcod)//".powder", status="replace",action="write")
       write(unit=lun,fmt="(/,/,6(a,/))")                                                 &
          "              ------ PROGRAM: POWDER PATTERN CALCULATION  ------"             , &
          "                     ---- Version 0.1 April-2009----"                        , &
          "    ***********************************************************************"  , &
          "    * Calculates powder diffraction patterns from a *.CFL or a *.CIF file *"  , &
          "    ***********************************************************************"  , &
          "                            (JRC- April 2009 )"
       ! Calculate a default Powder Diffraction Pattern
       stlmax=0.6; PPC%Title="Default Powder Pattern"
       PPC%U=0.0002; PPC%V=-0.0002; PPC%W=0.012; PPC%LAMBDA=1.54056; PPC%X=0.0015
       PPC%Thmin=1.00; PPC%step=0.05;  PPC%Thmax= int(2.0*asind(stlmax*1.54056)); PPC%job=0
       PPC%Ls=1900.0; PPC%Gs=  1900.0; nf=30; PPC%bkg=50.0
       powfile="powder_pattern.dat"
       units=" seconds"
       tim=0.0

       ! Write initial structure information in the .powder file
       call Write_Crystal_Cell(Cell,lun)
       call Write_SpaceGroup(SpG,lun)
       call Write_Atom_List(A,lun=lun)

       ! Look for calculation conditions in the CFL file that are provided before the list of atoms
       if (mode == "CFL") then
          do i=1,fich_cfl%nlines
             line=adjustl(u_case(fich_cfl%line(i)))
             if(line(1:4) == "ATOM") exit

             Select Case(Trim(line(1:7)))
                Case("TITLE")
                   PPC%title= adjustl(fich_cfl%line(i)(6:))

                Case("SUPCEL")
                   read(unit=line(7:),fmt=*,iostat=ier) n1,n2,n3,outa
                   if (ier == 0) then
                      ncel=(/n1,n2,n3/)
                      supcel_given=.true.
                      write(*,*) " => Supercell calculations"
                      if (outa == "Y") then
                         out_atoms=.true.
                         write(*,*) "    with output of generated atoms"
                      end if
                   end if

                Case("UVWX")
                   read(unit=line(7:),fmt=*,iostat=ier) PPC%U,PPC%V,PPC%W,PPC%X
                   if (ier /= 0) then
                      PPC%U=0.02; PPC%V=-0.02; PPC%W=0.12; PPC%X=0.0015
                   end if

                Case("REFLIST")
                   reflist=.true.
                   write(*,*) " => Reflection list will be output"

                Case("RANDOMB")
                   read(unit=line(8:),fmt=*,iostat=ier)  box
                   if(ier /= 0) box=0.5

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

                Case("SIZE_LG")
                   read(unit=line(8:),fmt=*,iostat=ier) PPC%Ls, PPC%Gs
                   if (ier /= 0) then
                      PPC%Ls=1900.0; PPC%Gs=  1900.0
                   end if
                   write(*,*) " => Lorentzian and Gaussian sizes: ", PPC%Ls, PPC%Gs

                Case("POWFILE")
                   powfile=adjustl(fich_cfl%line(i)(8:))
                   if (len_trim(powfile) == 0) then
                      powfile="powder_pattern.dat"
                   end if
                   write(*,*) " => Powder pattern file: "//trim(powfile)

             End Select
          end do

          if (supcel_given) then
             !Determine the box for random generation in fractional units
             fr=0.5*box/Cell%cell   !The 1/2 factor comes because box is the side in angstroms around the position
                                    !and the generation of atoms are in the two senses of each direction

             !Construct the supercell
             ad=Cell%cell*real(ncel); ang=Cell%ang
             call set_Crystal_Cell(ad,ang,cell)

             write(unit=lun,fmt="(/,a)")    " -- SUPERCELL Calculations, atoms randomly generated around ideal positions"
             write(unit=lun,fmt="(a,i4,a)") " -- Box side for atom generation 1/",nf," in fractional units of the subcell"
             write(unit=lun,fmt="(a,i3,2(a,i3))") "    SuperCell : ",ncel(1),"x",ncel(2),"x",ncel(3)
             write(unit=lun,fmt="(a,3f9.5)") "    Box Values: ",fr
             call write_crystal_cell(cell,lun)

             !Save coordinates of the read atoms
             nat=A%natoms*SpG%multip
             allocate(ch(nat),biso(nat),xyz(3,nat))
             m=0
             do i=1,A%natoms
                x=A%atom(i)%x
                call Get_Orbit(x,Spg,Mult,orb)
                do j=1,mult
                   m=m+1
                   xyz(:,m)=orb(:,j)
                   biso(m)=A%atom(i)%Biso
                   ch(m)=A%atom(i)%ChemSymb
                end do
             end do
             nat=m
             natoms=n1*n2*n3*nat
             write(*,*) " => Box Values: ",fr
             write(*,"(a,i8,a)") "  => A total number of ",natoms," atoms will be generated"   !", continue? "
             !read(*,*) ans
             !if(.not. (ans == 'y'.or. ans=='Y') ) stop
             call set_SpaceGroup("P 1",SpG)    !Space group is set to P 1  before structure factor calculations
             tim=0.0
             call cpu_time(tini)
             call gen_atoms(ch,xyz,Biso,nat,ncel,fr,Atm)
             if (out_atoms) then
                write(unit=lun,fmt="(a)")  " => List of atoms positions generated in the super-cell"
                do i=1,Atm%Natoms
                   write(unit=lun,fmt="(a,tr2,3f12.6,f12.4,f12.5)")  "Atom "//Atm%Atom(i)%Lab//Atm%Atom(i)%SfacSymb,   &
                                  Atm%Atom(i)%x,Atm%Atom(i)%biso,Atm%Atom(i)%occ
                end do
                call flush(lun)
             end if
             call cpu_time(tfin)
             tim=tim+ tfin-tini
             call factim(tfin-tini,time)
             write(*,*) " => CPU-time used for gen_atoms: ",(tfin-tini)*ftim,units
          else  !Normal calculation
             call Allocate_Atom_List(A%natoms,Atm,fail)
             if (fail) then
                write(*,"(a,i8,a)") "  => Error allocating the atom list for ", A%natoms," atoms"
                stop
             else
                Atm=A
             end if
          end if

       else ! CIF file no interaction
          call Allocate_Atom_List(A%natoms,Atm,fail)
          if (fail) then
             write(*,"(a,i8,a)") "  => Error allocating the atom list for ", A%natoms," atoms"
             stop
          else
             Atm=A
          end if
       End if

       ! Calculate sinTheta/Lambda max from 2Thetamax     PPC%Thmax= int(2.0*asind(stlmax*1.56))
       stlmax=sind(min((PPC%Thmax+10.0)*0.5,90.0))/PPC%lambda

    end if !if error

    if (PPC%job == 0) then      !X-rays
       write(unit=lun,fmt="(/,a)")  " => CALCULATION OF X-RAY POWDER DIFFRACTION PATTERN "
       PPC%title=Trim(PPC%title)//"; X-RAYS: "
    else
       write(unit=lun,fmt="(/,a)")  " => CALCULATION OF NEUTRON POWDER DIFFRACTION PATTERN"
       PPC%title=Trim(PPC%title)//"; NEUTRONS: "
    end if
    write(unit=lun,fmt="(  a,4f10.5)")  " => Resolution parameters UVWX: ",PPC%U,PPC%V,PPC%W,PPC%X
    write(unit=lun,fmt="(  a, f10.5)")  " => Lambda: ",PPC%lambda
    write(unit=lun,fmt="(  a, f10.5,a)")" => Background level: ",PPC%bkg," counts"
    write(unit=lun,fmt="(  a,2f10.2)")  " => Lorentzian and Gaussian sizes: ",PPC%Ls, PPC%Gs
    write(unit=lun,fmt="(  a,3f10.5)")  " => 2Theta range and step: ",PPC%Thmin,PPC%step,PPC%Thmax
    write(unit=lun,fmt="(  a,3f10.5)")  " => Maximum sin(Theta)/Lambda (for generating reflections): ",stlmax

    ! Now calculate a powder diffraction pattern
    ! First generate reflections and calculate structure factors
    Mult=2*SpG%NumOps
    if (supcel_given) Mult=1
       MaxNumRef = get_maxnumref(stlmax,Cell%CellVol,mult=Mult)
       call cpu_time(tini)
       call Hkl_Uni(Cell,SpG,.true.,0.0,stlmax,"s",MaxNumRef,hkl)
       call cpu_time(tfin)
       tim=tim+ tfin-tini
       call factim(tfin-tini,time)
       write(*,"(2(a,i8))") "  => Total number of effectively generated reflections is ",hkl%nref, &
                            "  Starting estimate was: ",MaxNumRef
       write(*,*) " => CPU-time used for Hkl_Uni: ",(tfin-tini)*ftim,units
       write(unit=lun,fmt="(a,i9)") " => Total number of generated reflections is ",hkl%nref
       write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time used for Hkl_Uni: ",(tfin-tini)*ftim,units

       call cpu_time(tini)
       if (PPC%job == 1) then      !Neutrons
          call Init_Calc_hkl_StrFactors(Atm,mode="NUC",lun=lun)
       else if(PPC%job == 0) then !Xrays
          call Init_Calc_hkl_StrFactors(Atm,mode="XRA",lambda=PPC%lambda,lun=lun)
       end if
       call cpu_time(tfin)
       tim=tim+ tfin-tini
       call factim(tfin-tini,time)
       write(*,*) " => CPU-time used for Init_Structure_Factors: ",(tfin-tini)*ftim,units
       write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time used for Init_Structure_Factors: ",(tfin-tini)*ftim,units
       if (ERR_SFac) then
          write(*,*) " => Error in initialization of Structure Factors"
          write(*,*) " => "//trim(ERR_SFac_Mess)
          stop
       end if

       call cpu_time(tini)
       write(*,*) " => Calculating structure factors ... please be patient!"
       if (PPC%job == 1) then      !Neutrons
          do i=1, hkl%nref
             sn=hkl%ref(i)%s * hkl%ref(i)%s
             h= hkl%ref(i)%h
             call Calc_hkl_StrFactor("P","N",h,sn,Atm,SpG,sf2)
             hkl%ref(i)%Fc=sqrt(sf2)
          end do
       else if(PPC%job == 0) then !X-rays
          do i=1, hkl%nref
             sn=hkl%ref(i)%s * hkl%ref(i)%s
             h= hkl%ref(i)%h
             call Calc_hkl_StrFactor("P","X",h,sn,Atm,SpG,sf2)
             hkl%ref(i)%Fc=sqrt(sf2)
          end do
       end if
       call cpu_time(tfin)
       tim=tim+ tfin-tini
       call factim(tfin-tini,time)
       write(*,*) " => CPU-time used for Structure_Factors: ",(tfin-tini)*ftim,units
       write(*,*) " => CPU-time used for Structure_Factors: "//trim(time)
       write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time used for Structure_Factors: ",(tfin-tini)*ftim,units
       write(unit=lun,fmt="(a)") " => CPU-time used for Structure_Factors: "//trim(time)

       if (ERR_SFac) then
          write(*,*) " => Error in calculations of Structure Factors"
          write(*,*) " => "//trim(ERR_SFac_Mess)
          stop
      end if
      if (reflist) then
         if (radiation(1:1) == "N") then
            call Write_Structure_Factors(lun,hkl,mode="NUC")
         else
            call Write_Structure_Factors(lun,hkl,mode="XRA")
         end if
      end if

      call cpu_time(tini)
      PPC%Scalef=cell%RCellVol
      call gen_powder_pattern(PPC,hkl,Pat)
      call cpu_time(tfin)
      tim=tim+ tfin-tini
      call factim(tfin-tini,time)
      write(*,*) " => CPU-time used for gen_powder_pattern: ",(tfin-tini)*ftim,units
      write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time used for gen_powder_pattern: ",(tfin-tini)*ftim,units
      call factim(tim,time)
      write(*,*) " => CPU-time for all calculations: ",tim*ftim,units
      write(*,*) " => CPU-time for all calculations: "//trim(time)
      write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time for all calculations: ",tim*ftim,units
      write(unit=lun,fmt="(a)") " => CPU-time for all calculations: "//trim(time)
      open(unit=lp,file=trim(powfile),status="replace",action="write")
      write(unit=lp,fmt="(a)") "!"//trim(Pat%title)
      write(unit=lp,fmt="(3f10.4)") Pat%xmin,Pat%step,Pat%xmax
      write(unit=lp,fmt="(8f16.4)") Pat%ycalc+PPC%bkg
      ! Alternative two-column output (comment two previous lines and uncomment the following ones)
      !write(unit=lp,fmt="(a)") "! twotheta  countrate "
      !do i=1,Pat%npts
      !    write(unit=lp,fmt="(2f16.4)") Pat%xmin+i*Pat%step,Pat%ycalc(i)+PPC%bkg
      !end do
      close(unit=lun)
      stop

 contains
     !!----
     !!----
     !!----
     !!----
     !!---- Update:
     !!
     Subroutine Factim(Timum,Time)
        !---- Arguments ----!
        real,             intent(in) :: timum
        character(len=35),intent(out):: time

        !---- Variables ----!
        integer :: seconds, minutes, hours
        integer :: scnds

        scnds=nint(timum)
        minutes=scnds/60
        seconds=mod(scnds,60)
        hours=minutes/60
        minutes=mod(minutes,60)
        write(unit=time,fmt="(3(i3,a))") hours," hours,",minutes," minutes,",seconds," seconds"
        Select Case(scnds)
           Case(60:3600)
              ftim=1.0/60.0
              units=" minutes"
           Case(3601:)
              ftim=1.0/60.0/60.0
              units=" hours"
           Case Default
              ftim=1.0
              units=" seconds"
        End Select

        return
     End Subroutine Factim
  End Program Calculation_of_Powder_Patterns

