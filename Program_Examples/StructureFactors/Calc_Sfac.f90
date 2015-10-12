!!----
!!---- Program: CALC_SFAC
!!----          Example of simple program using CFML
!!----
!!---- Author: Juan Rodriguez-Carvajal
!!---- Revision: November-2008
!!
Program Calc_Structure_Factors
   !---- Use Modules ----!
  !use f2kcli  !Comment for Fortran 2003 compliant (Environment module) compilers
   use CFML_crystallographic_symmetry, only: space_group_type, Write_SpaceGroup
   use CFML_Atom_TypeDef,              only: Atom_List_Type, Write_Atom_List
   use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Write_Crystal_Cell
   use CFML_Reflections_Utilities,     only: Reflection_List_Type, Hkl_Uni, get_maxnumref
   use CFML_IO_Formats,                only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
   use CFML_Structure_Factors,         only: Structure_Factors, Write_Structure_Factors, &
                                             Init_Structure_Factors,Calc_StrFactor
   use CFML_String_Utilities,          only: u_case

   !---- Variables ----!
   implicit none

   type (file_list_type)       :: fich_cfl
   type (space_group_type)     :: SpG
   type (Atom_list_Type)       :: A
   type (Crystal_Cell_Type)    :: Cell
   type (Reflection_List_Type) :: hkl

   character(len=256)          :: filcod     !Name of the input file
   character(len=132)          :: line
   character(len=15)           :: sinthlamb  !String with stlmax (2nd cmdline argument)
   real                        :: stlmax     !Maximum Sin(Theta)/Lambda
   real                        :: sn,sf2, Lambda
   integer                     :: MaxNumRef, Num, lun=1, ier,i
   complex                     :: fc

   integer                     :: narg
   Logical                     :: esta, arggiven=.false.,sthlgiven=.false.

   !---- Arguments on the command line ----!
   narg=command_argument_count()

   if (narg > 0) then
      call get_command_argument(1,filcod)
      arggiven=.true.
   end if

   if (narg > 1) then
      call get_command_argument(2,sinthlamb)
      read(unit=sinthlamb,fmt=*,iostat=ier) stlmax
      if (ier == 0) sthlgiven=.true.
   end if

   write(unit=*,fmt="(/,/,6(a,/))")                                                  &
        "            ------ PROGRAM STRUCTURE FACTORS ------"                      , &
        "              ---- Version 0.2 November-2008----"                         , &
        "    *******************************************************************"  , &
        "    * Calculates structure factors reading a *.CFL or a *.CIF file    *"  , &
        "    *******************************************************************"  , &
        "                      (JRC- November 2008 )"
   write(unit=*,fmt=*) " "

   if (.not. arggiven) then
      write(unit=*,fmt="(a)", advance='no') " => Code of the file xx.cif(cfl) (give xx): "
      read(unit=*,fmt="(a)") filcod
      if(len_trim(filcod) == 0) stop
   end if
   if (.not. sthlgiven) then
      write(unit=*,fmt="(a)", advance='no') " => Maximum sinTheta/Lambda: "
      read(unit=*,fmt=*) stlmax
   end if

   open(unit=lun,file=trim(filcod)//".sfa", status="replace",action="write")
   write(unit=lun,fmt="(/,/,6(a,/))")                                                  &
          "            ------ PROGRAM STRUCTURE FACTORS ------"                      , &
          "              ---- Version 0.2 November-2008----"                         , &
          "    *******************************************************************"  , &
          "    * Calculates structure factors reading a *.CFL or a *.CIF file    *"  , &
          "    *******************************************************************"  , &
          "                      (JRC- November 2008 )"

   inquire(file=trim(filcod)//".cif",exist=esta)
   if (esta) then
      call Readn_set_Xtal_Structure(trim(filcod)//".cif",Cell,SpG,A,Mode="CIF")
   else
      inquire(file=trim(filcod)//".cfl",exist=esta)
      if ( .not. esta) then
         write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cif (or .cfl) does'nt exist!"
         stop
      end if
      call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)
   end if

   if (err_form) then
      write(unit=*,fmt="(a)") trim(err_form_mess)
   else
      call Write_Crystal_Cell(Cell,lun)
      call Write_SpaceGroup(SpG,lun)
      call Write_Atom_List(A,lun=lun)
      !Look for wavelength in CFL file
      lambda=0.70926 !Mo kalpha (used only for x-rays)
       do i=1,fich_cfl%nlines
         line=adjustl(fich_cfl%line(i))
         if(U_Case(line(1:6)) == "LAMBDA") then
           read(unit=line(7:),fmt=*,iostat=ier) lambda
           if(ier /= 0) lambda=0.70926
         end if
       end do

      MaxNumRef = get_maxnumref(stlmax,Cell%CellVol,mult=SpG%NumOps)

      call Hkl_Uni(Cell,Spg,.true.,0.0,stlmax,"s",MaxNumRef,hkl)

      !> Calculation for neutron scattering
      call Init_Structure_Factors(hkl,A,Spg,mode="NUC",lun=lun)
      call Structure_Factors(A,SpG,hkl,mode="NUC")
      call Write_Structure_Factors(lun,hkl,mode="NUC")

      !> Test of another structure factor subroutine
      write(unit=lun,fmt="(/,a,/)") "   H   K   L   Mult  SinTh/Lda    |Fc|       Phase        F-Real      F-Imag      Num"
      do i=1, hkl%nref
         sn=hkl%ref(i)%s * hkl%ref(i)%s
         call Calc_StrFactor("P","N",i,sn,A,Spg,sf2,fc=fc)
         write(unit=lun,fmt="(3i4,i5,5f12.5,i8,f12.5)") hkl%ref(i)%h, hkl%ref(i)%mult, &
              hkl%ref(i)%S, hkl%ref(i)%Fc, hkl%ref(i)%Phase, real(fc), aimag(fc), i, sqrt(sf2)
      end do

      !> Calculation for X-rays assume Mo-kalpha if Lambda not given
      call Init_Structure_Factors(hkl,A,Spg,lambda=lambda,lun=lun)
      call Structure_Factors(A,SpG,hkl)
      call Write_Structure_Factors(lun,hkl)

      !> Calculation for Electron Diffraction
      call Init_Structure_Factors(hkl,A,Spg,Mode="ELE",lun=lun)
      call Structure_Factors(A,SpG,hkl,Mode="ELE")
      call Write_Structure_Factors(lun,hkl,Mode="ELE")

      write(unit=*,fmt="(a)") " Normal End of: PROGRAM STRUCTURE FACTORS "
      write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".sfa"
   end if

   !> Test of the type "file_list_type" by writing at the end of the file
   !do i=1,fich_cfl%nlines
   !   write(unit=lun,fmt="(a,i5,a)") " Line:",i,"  "//fich_cfl%line(i)
   !end do

   close(unit=lun)
End Program Calc_Structure_Factors

