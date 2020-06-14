module CFML_API_PowderSimulation

use CFML_GlobalDeps,  only: to_Deg
use CFML_Math_General,         only: asind,locate
use CFML_Atom_TypeDef,         only: Atom_Type, Atom_List_Type,Allocate_Atom_List
use CFML_Reflections_Utilities,only: Reflection_List_Type
use CFML_Diffraction_Patterns, only: Diffraction_Pattern_Type, Allocate_Diffraction_Pattern
use CFML_PowderProfiles_CW,    only: TCH_pVoigt,PseudoVoigt


implicit none

    Type :: PowPat_CW_Conditions
       character(len=140) :: title
       integer :: job
       real    :: Lambda, U, V, W, X, Ls, Gs
       real    :: Thmin, Thmax, step
       real    :: scalef,bkg
    End Type PowPat_CW_Conditions

CONTAINS

subroutine compute_powder_pattern(input_filename, mode, PowPat_Conditions, Pat)
     !---- Use Modules ----!
     use CFML_Math_General,              only: asind,sind
     use CFML_Atom_TypeDef,              only: Atom_List_Type, Allocate_Atom_List,Write_Atom_List
     use CFML_Crystal_Metrics,           only: Crystal_Cell_type, set_Crystal_Cell,write_crystal_cell
     use CFML_string_utilities,          only: u_case
     use CFML_Reflections_Utilities,     only: Reflection_List_Type,Hkl_uni,get_maxnumref
     Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Set_SpaceGroup, &
                                               Write_SpaceGroup
     Use CFML_Structure_Factors,         only: Write_Structure_Factors, Structure_Factors,&
                                               Init_Structure_Factors,ERR_SFac,ERR_SFac_Mess
     use CFML_Diffraction_Patterns,      only: Diffraction_Pattern_Type, &
                                               Allocate_Diffraction_Pattern
     use CFML_IO_Formats,                only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type

     !---- Variables ----!
     implicit none
     ! Input/Output
     character(len=*) :: input_filename
     character(len=3)       :: mode
     Type(PowPat_CW_Conditions)     :: PowPat_Conditions
     Type(Diffraction_Pattern_Type) :: Pat

     !
     integer                :: i,j,k,l,m,n,maxnumref,ier,mult,nf
     integer                :: lun=1,lp=2
     real, dimension(3)     :: ad,ang,x,fr
     character(len=1)       :: ans,outa
     integer, dimension(3)  :: ncel,h
     real                   :: stlmax,tini,tfin,sn,sf2,tim,ftim=1.0,box
     character(len=132)     :: line
     character(len=8)       :: units="seconds",radiation
     character(len=4),dimension(:),allocatable :: ch
     !
     Type(Crystal_Cell_type)        :: cell
     Type(Space_Group_Type)         :: SpG
     Type(Atom_List_Type)           :: A
     Type(Reflection_List_Type)     :: hkl
     Type(file_list_type)           :: fich_cfl
     integer                        :: narg
     Logical                        :: esta, arggiven=.false., fail

     if (Mode .eq. "CIF") then
        call Readn_set_Xtal_Structure(trim(input_filename),Cell,SpG,A,Mode="CIF")
     else
        fich_cfl%nlines=0
        call Readn_set_Xtal_Structure(trim(input_filename),Cell,SpG,A,Mode="CFL",file_list=fich_cfl)
        if (Mode .eq. "CFL") then
            call Read_Power_Pattern_Simulation_conditions(fich_cfl,PowPat_Conditions)
        end if
     end if

    if (err_form) then
        write(*,*) trim(err_form_mess)
    end if

    ! Calculate sinTheta/Lambda max from 2Thetamax     PowPat_Conditions%Thmax= int(2.0*asind(stlmax*1.56))
    stlmax=sind(min((PowPat_Conditions%Thmax+10.0)*0.5,90.0))/PowPat_Conditions%lambda

    ! Now calculate a powder diffraction pattern
    ! First generate reflections and calculate structure factors
    Mult=2*SpG%NumOps
    MaxNumRef = get_maxnumref(stlmax,Cell%CellVol,mult=Mult)
    call Hkl_Uni(Cell,SpG,.true.,0.0,stlmax,"s",MaxNumRef,hkl)

    if (PowPat_Conditions%job == 1) then      !Neutrons
       call Init_Structure_Factors(hkl,A,Spg,mode="NUC",lun=lun)
    else if(PowPat_Conditions%job == 0) then !Xrays
       call Init_Structure_Factors(hkl,A,Spg,mode="XRA",lambda=PowPat_Conditions%lambda,lun=lun)
    end if

    if (PowPat_Conditions%job == 1) then      !Neutrons
       call Structure_Factors(A,SpG,hkl,mode="NUC")
    else if(PowPat_Conditions%job == 0) then !X-rays
       call Structure_Factors(A,SpG,hkl,mode="XRA",lambda=PowPat_Conditions%lambda)
    end if

    if (ERR_SFac) then
       write(*,*) " => Error in calculations of Structure Factors"
       write(*,*) " => "//trim(ERR_SFac_Mess)
       stop
    end if

    if (PowPat_Conditions%job == 1) then
       call Write_Structure_Factors(lun,hkl,mode="NUC")
    else
       call Write_Structure_Factors(lun,hkl,mode="XRA")
    end if

    PowPat_Conditions%Scalef=cell%RCellVol
    call Calc_powder_pattern(PowPat_Conditions,hkl,Pat)

end subroutine


subroutine Read_Power_Pattern_Simulation_conditions(fich_cfl,PowPat_Conditions)
     use CFML_IO_Formats,                only: file_list_type
     use CFML_string_utilities,          only: u_case

     !---- Variables ----!
     !
     implicit none
     Type(file_list_type)           :: fich_cfl
     Type(PowPat_CW_Conditions)     :: PowPat_Conditions

     character(len=8)       :: radiation
     integer                :: ier

     !
     integer :: i
     character(len=132)     :: line

     do i=1,fich_cfl%nlines
         line=adjustl(u_case(fich_cfl%line(i)))
         if(line(1:4) == "ATOM") exit

         Select Case(Trim(line(1:7)))
            Case("TITLE")
               PowPat_Conditions%title= adjustl(fich_cfl%line(i)(6:))

            Case("UVWX")
               read(unit=line(7:),fmt=*,iostat=ier) PowPat_Conditions%U,PowPat_Conditions%V,PowPat_Conditions%W,PowPat_Conditions%X
               if (ier /= 0) then
                  PowPat_Conditions%U=0.02; PowPat_Conditions%V=-0.02; PowPat_Conditions%W=0.12; PowPat_Conditions%X=0.0015
               end if

            Case("LAMBDA")
               read(unit=line(7:),fmt=*,iostat=ier) PowPat_Conditions%LAMBDA
               if(ier /= 0) PowPat_Conditions%LAMBDA=1.56

            Case("BACKGD")
               read(unit=line(7:),fmt=*,iostat=ier) PowPat_Conditions%bkg
               if(ier /= 0) PowPat_Conditions%bkg=20.0

            Case("JOBTYPE")
               radiation = adjustl(fich_cfl%line(i)(8:))
               if (radiation(1:1) == "N" .or. radiation(1:1) == "n") then
                  PowPat_Conditions%job=1
               else
                  PowPat_Conditions%job=0
               end if

            Case("PATTERN")
               read(unit=line(8:),fmt=*,iostat=ier) PowPat_Conditions%Thmin, PowPat_Conditions%step,  PowPat_Conditions%Thmax
               if (ier /= 0) then
                  PowPat_Conditions%Thmin=1.00; PowPat_Conditions%step=0.05;  PowPat_Conditions%Thmax= 120.0
               end if

            Case("SIZE_LG")
               read(unit=line(8:),fmt=*,iostat=ier) PowPat_Conditions%Ls
               if (ier /= 0) then
                  PowPat_Conditions%Ls=1900.0
               end if

            Case("POWFILE")
               ! No output file

         End Select
     end do

end subroutine


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

Subroutine Calc_Powder_Pattern(PowPat_Conditions,Hkl,Pat)
   !---- Argument ----!
   Type(PowPat_CW_Conditions),     intent(in)  :: PowPat_Conditions
   Type(Reflection_List_Type),     intent(in)  :: hkl
   Type(Diffraction_Pattern_Type), intent(out) :: Pat

   !--- Local Variables ----!
   integer :: i,j,npts,i1,i2
   real    :: step,Intens,Bragg,Hl,Hg, ss,cs,tt,th1,th2,LorentzF, Y,eta,fwhm,chw

   npts=(PowPat_Conditions%Thmax-PowPat_Conditions%Thmin)/PowPat_Conditions%step + 1.02
   call Allocate_Diffraction_Pattern(Pat,npts)
   Pat%Title=adjustl(Trim(PowPat_Conditions%title))
   i=len_trim(Pat%Title)
   write(unit=Pat%Title(i+2:),fmt="(a,f7.4,f7.1)") " => lambda,Ls: ", &
              PowPat_Conditions%Lambda,PowPat_Conditions%Ls

   Pat%scat_var="2-Theta"
   Pat%instr="Calculated Pattern"
   Pat%xmin= PowPat_Conditions%Thmin
   Pat%xmax= PowPat_Conditions%Thmax
   Pat%ymin= 0.0
   Pat%ymax=0.0
   Pat%scal=1.0
   Pat%monitor=0.0
   Pat%step=PowPat_Conditions%step
   Pat%Tsamp=300.0
   Pat%Tset=300.0
   Pat%npts=npts
   Pat%ct_step=.true.
   Pat%conv=(/PowPat_Conditions%Lambda,PowPat_Conditions%Lambda,0.0,0.0,0.0/)
   chw=15.0
   do i=1,npts
      Pat%x(i)=Pat%xmin+real(i-1)*Pat%step
   end do

   Y= to_deg*PowPat_Conditions%Lambda/PowPat_Conditions%Ls
   do i=1,hkl%nref
      ss=PowPat_Conditions%Lambda*hkl%ref(i)%S
      cs=sqrt(abs(1.0-ss*ss))
      tt=ss/cs
      LorentzF=0.5/(ss*ss*cs)
      Bragg=2.0*asind(ss)
      HG=sqrt(tt*(PowPat_Conditions%U*tt+PowPat_Conditions%V)+PowPat_Conditions%W)
      HL=PowPat_Conditions%X*tt + Y/cs
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
      Intens= LorentzF *hkl%ref(i)%mult * hkl%ref(i)%Fc**2 * PowPat_Conditions%Scalef
      do j=i1,i2
         Pat%ycalc(j)=Pat%ycalc(j)+ PseudoVoigt( Pat%x(j)-Bragg, (/fwhm,eta /) ) * Intens
      end do
   end do

   return
End Subroutine Calc_Powder_Pattern

end module
