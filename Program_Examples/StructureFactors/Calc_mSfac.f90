!!----
!!---- Program: Calc_Magnetic_Structure_Factors
!!----          Example of simple program using CFML
!!----
!!---- Author: Juan Rodriguez-Carvajal
!!---- Revision: January 2020
!!
Program Calc_Magnetic_Structure_Factors
   !---- Use Modules ----!
  !use f2kcli  !Comment for Fortran 2003 compliant (Environment module) compilers
   use CFML_crystallographic_symmetry, only: magnetic_space_group_type,&
                                             Write_Magnetic_Space_Group 
   use CFML_Atom_TypeDef,              only: Atom_List_Type, Write_Atom_List
   use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Write_Crystal_Cell
   use CFML_Reflections_Utilities,     only: Reflect_List_Type
   use CFML_IO_Formats,                only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
   use CFML_Structure_Factors,         only: Magnetic_Structure_Factors, Write_Structure_Factors, & 
                                             Strf_List_Type 
   use CFML_String_Utilities,          only: u_case

   !---- Variables ----!
   implicit none

   type (file_list_type)            :: fich_cfl
   type (magnetic_space_group_type) :: SpG
   type (Atom_list_Type)            :: A
   type (Crystal_Cell_Type)         :: Cell
   type (Reflect_List_Type)         :: hkl
   type (Strf_List_Type)            :: Stf

   character(len=256)          :: filcod     !Name of the input file
   character(len=132)          :: line
   character(len=15)           :: sinthlamb  !String with stlmax (2nd cmdline argument)
   real                        :: stlmax     !Maximum Sin(Theta)/Lambda
   real                        :: Lambda
   integer                     :: lun=1, ier,i
   integer                     :: narg
   Logical                     :: esta, arggiven=.false.,sthlgiven=.false., full=.true.

   !---- Arguments on the command line ----!
   narg=command_argument_count()

   if (narg > 0) then
      call get_command_argument(1,filcod)
      i=index(filcod,".")
      if( i/= 0) filcod=filcod(1:i-1)
      arggiven=.true.
   end if

   if (narg > 1) then
      call get_command_argument(2,sinthlamb)
      read(unit=sinthlamb,fmt=*,iostat=ier) stlmax
      if (ier == 0) sthlgiven=.true.
   end if

   write(unit=*,fmt="(/,/,6(a,/))")                                                  &
        "            ------ PROGRAM Magnetic STRUCTURE FACTORS ------"             , &
        "                 ---- Version 0.0 January 2020----"                       , &
        "    *******************************************************************"  , &
        "    * Calculates structure factors reading a *.CFL or an *.mCIF file  *"  , &
        "    *******************************************************************"  , &
        "                      (JRC- January 2020 )"
   write(unit=*,fmt=*) " "

   if (.not. arggiven) then
      write(unit=*,fmt="(a)", advance='no') " => Code of the file xx.mcif(cfl) (give xx): "
      read(unit=*,fmt="(a)") filcod
      if(len_trim(filcod) == 0) stop
   end if
   if (.not. sthlgiven) then
      write(unit=*,fmt="(a)", advance='no') " => Maximum sinTheta/Lambda: "
      read(unit=*,fmt=*) stlmax
   end if

   open(unit=lun,file=trim(filcod)//".sfa", status="replace",action="write")
   write(unit=lun,fmt="(/,/,6(a,/))")                                                  &
          "            ------ PROGRAM Magnetic STRUCTURE FACTORS ------"             , &
          "              ---- Version 0.0 January-2020----"                          , &
          "    *******************************************************************"  , &
          "    * Calculates structure factors reading a *.CFL or a *.mCIF file   *"  , &
          "    *******************************************************************"  , &
          "                      (JRC- January-2020 )"

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
   end if

   if (err_form) then
      write(unit=*,fmt="(a)") trim(err_form_mess)
   else
      call Write_Crystal_Cell(Cell,lun)
      call Write_Magnetic_Space_Group(SpG,lun,full)
      call Write_Atom_List(A,level=2,lun=lun)
      !Look for wavelength in CFL file
      lambda=0.70926 !Mo kalpha (used only for x-rays)
       do i=1,fich_cfl%nlines
         line=adjustl(fich_cfl%line(i))
         if(U_Case(line(1:6)) == "LAMBDA") then
           read(unit=line(7:),fmt=*,iostat=ier) lambda
           if(ier /= 0) lambda=0.70926
         end if
       end do

      call Magnetic_Structure_Factors(Cell,A,SpG,stlmax,hkl,Stf,lun)
      call Write_Structure_Factors(lun,hkl,stf,full)


      write(unit=*,fmt="(a)") " Normal End of: PROGRAM Magnetic STRUCTURE FACTORS "
      write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".sfa"
   end if

    close(unit=lun)
End Program Calc_Magnetic_Structure_Factors

