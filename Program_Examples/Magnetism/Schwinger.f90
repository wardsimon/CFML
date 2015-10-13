Program Schwinger
 use CFML_GlobalDeps
 use CFML_Math_3D,                  only: cross_product
 use CFML_crystallographic_symmetry,only: space_group_type, Write_SpaceGroup
 use CFML_Atom_TypeDef,             only: Atom_List_Type, Write_Atom_List
 use CFML_crystal_metrics,          only: Crystal_Cell_Type, Write_Crystal_Cell
 use CFML_Reflections_Utilities,    only: Hkl_s
 use CFML_String_Utilities,         only: l_case,u_case
 use CFML_IO_Formats,               only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
 use CFML_Scattering_Chemical_Tables

 implicit none

 type (file_list_type)        :: fich_cfl
 type (Space_Group_Type)      :: SpG
 type ( Atom_list_Type)       :: A
 type (Crystal_Cell_Type)     :: Cell

 character(len=256)           :: line,filcod,mess
 real(kind=cp)                :: sn,s2,theta,lambda,flip_right,flip_left,up,down
 real(kind=cp), dimension(3)  :: hn
 integer                      :: lun=1, ier,i,j,ih,ik,il, n_ini, n_end
 complex(kind=cp)             :: fn,fx,fe,fsru,fsrd,fslu,fsld
 real(kind=cp), dimension(6)  :: extc
 real(kind=cp), dimension(3,3):: ub

 integer                     :: narg,iext,ity
 Logical                     :: esta, arggiven=.false.,left=.false.,ext=.false.,ubgiven=.false.,ok

 Type :: Scattering_Species_Type
    integer                                        :: Num_Species
    character(len=6),    dimension(:), allocatable :: Symb
    real(kind=cp),       dimension(:), allocatable :: br,bi
    real(kind=cp),       dimension(:), allocatable :: delta_fp,delta_fpp
    type(Xray_Form_Type),dimension(:), allocatable :: Xcoef
 End Type Scattering_Species_Type
 Type(Scattering_Species_Type) :: Scattf, add_Scat

 real(kind=dp), parameter :: pn  =  0.2695420113693928312
 real(kind=dp), parameter :: schw= -0.00014699

      !---- Arguments on the command line ----!
      narg=COMMAND_ARGUMENT_COUNT()

      if(narg > 0) then
              call GET_COMMAND_ARGUMENT(1,filcod)
              arggiven=.true.
      end if

      write(unit=*,fmt="(/,/,7(a,/))")                                                 &
           "              ------ P r o g r a m   S c h w i n g e r  ------"          , &
           "                  ---- Version 0.0 October-2015 ----"                    , &
           "    ******************************************************************"  , &
           "    *    Calculates the amplitude of the Schwinger Scattering for     *"  , &
           "    * individual reflections. The structure is read from a *.CFL file *"  , &
           "    ******************************************************************"  , &
           "                            (JRC- October 2015)"
    write(unit=*,fmt=*) " "

     if(.not. arggiven) then
       write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.cfl (give xx): "
       read(unit=*,fmt="(a)") filcod
       if(len_trim(filcod) == 0) stop
     end if
     i=index(filcod,".",back=.true.)
     if( i /= 0) filcod=filcod(1:i-1)

     open(unit=lun,file=trim(filcod)//".cal", status="replace",action="write")
     write(unit=lun,fmt="(/,/,7(a,/))")                                                &
           "              ------ P r o g r a m   S c h w i n g e r  ------"          , &
           "                  ---- Version 0.0 October-2015 ----"                    , &
           "    ******************************************************************"  , &
           "    *    Calculates the amplitude of the Schwinger Scattering for     *"  , &
           "    * individual reflections. The structure is read from a *.CFL file *"  , &
           "    ******************************************************************"  , &
           "                            (JRC- October 2015)"

     inquire(file=trim(filcod)//".cfl",exist=esta)
     if( .not. esta) then
       write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cfl does'nt exist!"
       stop
     end if
     call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)

     If(err_form) then
       write(unit=*,fmt="(a)") trim(err_form_mess)
     else

       !Test of the type "file_list_type" by writing at the end of the file
       write(unit=lun,fmt="(/,a)") "    =========================="
       write(unit=lun,fmt="( a )") "    Text of the input CFL file"
       write(unit=lun,fmt="(a,/)") "    =========================="
       do i=1,fich_cfl%nlines
         write(unit=lun,fmt="(a,i5,a)") " Line:",i,"  "//fich_cfl%line(i)
       end do


       call Write_Crystal_Cell(Cell,lun)
       call Write_SpaceGroup(SpG,lun,full=.true.)
       call Write_Atom_List(A,lun=lun)
       lambda=0.8
       do i=1,fich_cfl%nlines
         line=adjustl(fich_cfl%line(i))
         if(U_Case(line(1:6)) == "LAMBDA") then
           read(unit=line(7:),fmt=*,iostat=ier) lambda
           if(ier /= 0) lambda=0.8
         end if
         if(U_Case(line(1:4)) == "EXTI") then
           read(unit=line(5:),fmt=*,iostat=ier) extc, iext
           if(ier /= 0) then
              read(unit=line(5:),fmt=*,iostat=ier) extc(1)
              if(ier /= 0) extc(:) = 0.0
              iext=1
           end if
           if(any(extc > 0.0001)) ext=.true.
         end if
       end do

       call Additional_Scattering_Factors(fich_cfl,add_Scat,ok,mess)
       if(.not. ok) then
         write(unit=*,fmt="(a)") " => "//trim(mess)
         stop
       end if

       !  Set nuclear, x-ray, neutron and magnetic form factors coefficients for all the atoms
       !  in the structure
       if(add_Scat%Num_Species > 0) then
         call Set_Form_Factors(A,Scattf,ok,mess,lambda,lun,add_Scat)
       else
         call Set_Form_Factors(A,Scattf,ok,mess,lambda,lun)
       end if
       if(.not. ok) then
         write(unit=*,fmt="(a)") " => "//trim(mess)
         stop
       end if
       do
         write(unit=*,fmt="(a)",advance="no")  " => Enter a reflection as 3 integers -> (h,k,l): "
         read(unit=*,fmt=*,iostat=ier) ih,ik,il
         if(ier /= 0) cycle
         if(abs(ih)+abs(ik)+abs(il) == 0 ) exit
         hn=real([ih,ik,il])
         sn = hkl_s(hn,Cell)
         theta=lambda*sn
         theta=asin(theta)
         s2=sn*sn
         call Calc_General_StrFactor(hn,s2,A,SpG,Scattf,fn,fx,fe)
         fsru=Schwinger_Amplitude(hn,[0.0,0.0, 1.0],theta,fe,Left=.false.)
         fslu=Schwinger_Amplitude(hn,[0.0,0.0, 1.0],theta,fe,Left=.true.)
         fsrd=Schwinger_Amplitude(hn,[0.0,0.0,-1.0],theta,fe,Left=.false.)
         fsld=Schwinger_Amplitude(hn,[0.0,0.0,-1.0],theta,fe,Left=.true.)
         up= (fn+fsru)*conjg(fn+fsru); down= (fn+fsrd)*conjg(fn+fsrd)
         if(up < 1.0e-6 .and. down < 1.0e-6) then
            flip_right=1.0
         else
            flip_right=up/down
         end if
         up=(fn+fslu)*conjg(fn+fslu); down= (fn+fsld)*conjg(fn+fsld)
         if(up < 1.0e-6 .and. down < 1.0e-6) then
            flip_left=1.0
         else
            flip_left=up/down
         end if
         write(unit=*,fmt="(/,a,3i4,a,f8.5)") "  Reflection: (",ih,ik,il,") sinTheta/Lambda=",sn
         write(unit=*,fmt="(a,2(f10.5,a))")   "  Nuclear     Structure   Factor : (",real(fn),  " ) + i (",aimag(fn),  " ) "
         write(unit=*,fmt="(a,2(f10.5,a))")   "  X-rays      Structure   Factor : (",real(fx),  " ) + i (",aimag(fx),  " ) "
         write(unit=*,fmt="(a,2(f10.5,a))")   "  Electrostatic Structure Factor : (",real(fe),  " ) + i (",aimag(fe),  " ) "
         write(unit=*,fmt="(a,2(f10.5,a))")   "  Schwinger Amplitude  right-up  : (",real(fsru)," ) + i (",aimag(fsru)," ) "
         write(unit=*,fmt="(a,2(f10.5,a))")   "  Schwinger Amplitude  right-down: (",real(fsrd)," ) + i (",aimag(fsrd)," ) "
         write(unit=*,fmt="(a,2(f10.5,a))")   "  Schwinger Amplitude  left-up   : (",real(fslu)," ) + i (",aimag(fslu)," ) "
         write(unit=*,fmt="(a,2(f10.5,a))")   "  Schwinger Amplitude  left-down : (",real(fsld)," ) + i (",aimag(fsld)," ) "
         write(unit=*,fmt="(a, f12.6)")      "            Flipping ratio right : ", flip_right
         write(unit=*,fmt="(a, f12.6)")      "            Flipping ratio left  : ", flip_left

         write(unit=lun,fmt="(/,a,3i4,a,f8.5)") "  Reflection: (",ih,ik,il,") sinTheta/Lambda=",sn
         write(unit=lun,fmt="(a,2(f10.5,a))")   "  Nuclear     Structure   Factor : (",real(fn),  " ) +i (",aimag(fn),  " ) "
         write(unit=lun,fmt="(a,2(f10.5,a))")   "  X-rays      Structure   Factor : (",real(fx),  " ) +i (",aimag(fx),  " ) "
         write(unit=lun,fmt="(a,2(f10.5,a))")   "  Electrostatic Structure Factor : (",real(fe),  " ) +i (",aimag(fe),  " ) "
         write(unit=lun,fmt="(a,2(f10.5,a))")   "  Schwinger Amplitude  right-up  : (",real(fsru)," ) +i (",aimag(fsru)," ) "
         write(unit=lun,fmt="(a,2(f10.5,a))")   "  Schwinger Amplitude  right-down: (",real(fsrd)," ) +i (",aimag(fsrd)," ) "
         write(unit=lun,fmt="(a,2(f10.5,a))")   "  Schwinger Amplitude  left-up   : (",real(fslu)," ) +i (",aimag(fslu)," ) "
         write(unit=lun,fmt="(a,2(f10.5,a))")   "  Schwinger Amplitude  leftt-down: (",real(fsld)," ) +i (",aimag(fsld)," ) "
         write(unit=lun,fmt="(a, f12.6)")      "            Flipping ratio right : ", flip_right
         write(unit=lun,fmt="(a, f12.6)")      "            Flipping ratio left  : ", flip_left
       end do
       write(unit=*,fmt="(a)") " Normal End of program: Schwinger "
       write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".cal"
     end if

     close(unit=lun)

    contains

    Subroutine Additional_Scattering_Factors(fil,add_Scatt,ok,mess)
      Type(File_List_Type),          intent(in)  :: fil
      Type(Scattering_Species_Type), intent(out) :: add_Scatt
      logical,                       intent(out) :: ok
      character(len=*),              intent(out) :: mess
      !Local variables
      integer, parameter :: N_add = 20
      character(len=132) :: line
      character(len=4), dimension(N_add)   :: names
      real(kind=cp),    dimension(N_add)   :: b_real,b_imag
      real(kind=cp),    dimension(N_add)   :: d_fp,d_fpp,cc
      real(kind=cp),    dimension(4,N_add) :: ac,bc
      integer :: i,nsp,j
      ok=.true.
      mess=" "
      b_real=0.0; b_imag=0.0; d_fp=0.0; d_fpp=0.0; cc=0.0
      ac=0.0; bc=0.0
      nsp=0
      do i=1,fil%nlines
        line=adjustl(fil%line(i))
        if(U_case(line(1:2)) == "B_") then
          nsp=nsp+1
          j=index(line," ")
          names(nsp)=line(3:j-1)
          read(unit=line(j:),fmt=*,iostat=ier) b_real(nsp),b_imag(nsp)
          if(ier /= 0) then
            ok=.false.
            mess="Error reading scattering length on line containing: "//trim(line)
            return
          end if
        end if
        if(U_case(line(1:5)) == "DELT_") then
          nsp=nsp+1
          j=index(line," ")
          names(nsp)=line(6:j-1)
          read(unit=line(j:),fmt=*,iostat=ier) d_fp(nsp),d_fpp(nsp)
          if(ier /= 0) then
            ok=.false.
            mess="Error reading anomalous scattering terms on line containing: "//trim(line)
            return
          end if
        end if
        if(U_case(line(1:7)) == "XCOEFF_") then
          nsp=nsp+1
          j=index(line," ")
          names(nsp)=line(8:j-1)
          read(unit=line(j:),fmt=*,iostat=ier) ac(:,nsp),bc(:,nsp),cc(nsp)
          if(ier /= 0) then
            ok=.false.
            mess="Error reading X-ray scattering coefficients on line containing: "//trim(line)
            return
          end if
        end if
      end do
      if(nsp > N_add) then
        ok=.false.
        write(unit=mess,fmt="(a,i3,a)") "The number of additional scattering factors is limited to ",N_add," !!!"
        return
      end if
      if(nsp > 0) then
        call Allocate_Scattering_Species(nsp,add_Scatt)
        do i=1,nsp
          add_Scatt%Symb(i)      = names(i)
          add_Scatt%br(i)        = b_real(i)
          add_Scatt%bi(i)        = b_imag(i)
          add_Scatt%delta_fp(i)  = d_fp(i)
          add_Scatt%delta_fpp(i) = d_fpp(i)
          add_Scatt%Xcoef(i)%Symb= names(i)
          add_Scatt%Xcoef(i)%a   = ac(:,i)
          add_Scatt%Xcoef(i)%b   = bc(:,i)
          add_Scatt%Xcoef(i)%c   = cc(i)
        end do
      else
        add_Scatt%Num_species=0
      end if

    End Subroutine Additional_Scattering_Factors

    Subroutine Set_Form_Factors(Atm,Scf,ok,mess,lambda,lun,Add_Scatt)
      type(Atom_List_Type),                             intent(in out):: Atm
      Type(Scattering_Species_Type),                    intent(out)   :: Scf
      logical,                                          intent(out)   :: ok
      character(len=*),                                 intent(out)   :: mess
      real(kind=cp),                optional,           intent(in)    :: lambda
      integer,                      optional,           intent(in)    :: lun
      Type(Scattering_Species_Type),optional,           intent(in)    :: Add_Scatt

      !---- Local variables ----!
      character(len=4)                        :: symbcar
      integer                                 :: i,k,n,m,L
      character(len=4), dimension(atm%natoms) :: symb
      character(len=4), dimension(atm%natoms) :: elem
      real(kind=cp),    dimension(atm%natoms) :: bs
      real(kind=cp)                           :: b,dmin,d
      character(len=*), parameter             :: digpm="0123456789+-"
      logical                                 :: found

      call set_chem_info()
       !---- Getting Fermi Lengths of atoms ----!
       symb="    "; Elem="    "
       bs=0.0
       n=0
       ok=.true.
       do i=1,atm%natoms
          symbcar=u_case(atm%atom(i)%chemsymb)
          call Get_ChemSymb(symbcar, atm%atom(i)%chemsymb, atm%atom(i)%Z) !Getting the atomic number
          call Get_Fermi_Length(symbcar,b)                                !equal to the charge of the nuclei
          if (abs(b) < 0.0001) then
             ok=.false.
             Mess="The Fermi Length of Species "//symbcar//" was not found"
             return
          else
             if(any(Elem == symbcar)) cycle
             n=n+1
             bs(n) = b
             Elem(n)=u_case(atm%atom(i)%chemsymb)
             symb(n)=atm%atom(i)%SfacSymb
          end if
       end do
       call Remove_chem_info()

       Call Allocate_Scattering_Species(n,Scf)

       Scf%br=bs(1:n)

       if(present(lambda)) then
         call Set_Delta_Fp_Fpp()

          !---- Select wavelength (by default is CuKalpha1: k=5 in the list) ----!
          dmin=1000.0
          do i=1,5
             d=abs(lambda-Xray_Wavelengths(i)%Kalfa(1))
             if (d < dmin) then
                dmin=d
                k=i        !Selection of the index for fp and fpp lists
             end if
          end do

          !---- Found Species on Anomalous_ScFac ----!
          do i=1,Scf%Num_species
             symbcar=l_case(Elem(i))
             do j=1,Num_Delta_Fp
                if (symbcar /= Anomalous_ScFac(j)%Symb) cycle
                Scf%delta_fp(i)=Anomalous_ScFac(j)%fp(k)
                Scf%delta_fpp(i)=Anomalous_ScFac(j)%fpp(k)
                exit
             end do
          end do
         call Remove_Delta_Fp_Fpp()
       end if

       call Set_Xray_Form()

       !Look for X-ray scattering coefficients
       do i=1,Scf%Num_species
          symbcar=l_case(symb(i)) !Scattering factors
          k=index(symbcar,"+")
          j=index(symbcar,"-")
          if(k == 0 .and. j == 0) then !simple element or magnetic form Factor
            if(len_trim(symbcar) > 2) then !magnetic form Factor -> use the chemical symbol
               symbcar=l_case(Elem(i))
            end if
          end if
          found=.false.
          do j=1,Num_Xray_Form
             if (symbcar /= Xray_form(j)%Symb) cycle
             Scf%xcoef(i)=Xray_form(j)
             found=.true.
             exit
          end do
          if(.not. found) then
            ok=.false.
            mess="Error: X-ray scattering form factor coefficients not found for "//symbcar
            return
          end if
       end do
       call Remove_Xray_Form()

       m=0
       do i=1,atm%natoms
         symbcar=u_case(atm%atom(i)%chemsymb)
         if(u_case(trim(atm%atom(i)%SfacSymb)) == "MPOL") m=m+1
         do j=1,Scf%Num_species
           if(symbcar == Elem(j)) then
             atm%atom(i)%ind(1)=j
             Scf%Xcoef(j)%Z=atm%atom(i)%Z
             Scf%Symb(j)=atm%atom(i)%SfacSymb
             exit
           end if
         end do
       end do

       if(present(Add_Scatt)) then
         if(Add_Scatt%Num_Species > 0) then
           do i=1,Add_Scatt%Num_Species
             do j=1,Scf%Num_species
                if(Scf%Symb(j) == Add_Scatt%Symb(i)) then
                  if(abs(Add_Scatt%br(i))> 0.00001)  Scf%br(j)=Add_Scatt%br(i)
                  if(abs(Add_Scatt%bi(i))> 0.00001)  Scf%bi(j)=Add_Scatt%bi(i)
                  if(abs(Add_Scatt%delta_fp(i))> 0.00001)  Scf%delta_fp(j) =Add_Scatt%delta_fp(i)
                  if(abs(Add_Scatt%delta_fpp(i))> 0.00001) Scf%delta_fpp(j)=Add_Scatt%delta_fpp(i)
                  if(abs(sum(Add_Scatt%Xcoef(i)%a))> 0.00001) then
                    Scf%Xcoef(j)%a=Add_Scatt%Xcoef(i)%a
                    Scf%Xcoef(j)%b=Add_Scatt%Xcoef(i)%b
                    Scf%Xcoef(j)%c=Add_Scatt%Xcoef(i)%c
                  end if
                end if
             end do
           end do
         end if
       end if

       !---- Printing Information ----!
       if (present(lun)) then
          if(present(lambda)) then
            write(unit=lun,fmt="(/,a,f10.6,a)")  "  WAVELENGTH: ",lambda," Angstroms"
          else
            write(unit=lun,fmt="(/,a)")  "  WAVELENGTH NOT PROVIDED! "
          end if
          write(unit=lun,fmt="(/,a)")  "  INFORMATION FROM TABULATED NEUTRON SCATTERING FACTORS"
          write(unit=lun,fmt="(a,/)")  "  ==================================================="
          write(unit=lun,fmt="(a)")    "  FERMI LENGTHS "
          write(unit=lun,fmt="(a,i3)") "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)")  "   Atom     Fermi Length (Br,Bi)[10^(-12) cm]      Atomic Number"
          do k=1,Scf%Num_species
             write(unit=lun,fmt="(a,2F10.6,tr20,i8)")  "     "//Scf%Symb(k), Scf%br(k), Scf%bi(k), Scf%Xcoef(k)%Z
          end do
          write(unit=lun,fmt="(/,/)")
          write(unit=lun,fmt="(/,a)")  "  INFORMATION FROM TABULATED X-RAY SCATTERING FACTORS"
          write(unit=lun,fmt="(a,/)")  "  ==================================================="
          write(unit=lun,fmt="(/,a,/)")    "   ATOMIC SCATTERING FACTOR COEFFICIENTS: {A(i),B(i),I=1,4},C  Dfp  Dfpp "
          write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)") &
               " Atom(ScF)  Atom(ScFX)           a1       b1       a2       b2       a3       b3       a4       b4        c      Dfp     Dfpp"
          do k=1,Scf%Num_species
             write(unit=lun,fmt="(a,11F9.5)")    &
                            "  "//Scf%Symb(k)//"        "//Scf%Xcoef(k)%Symb//"        ", &
                           (Scf%xcoef(k)%a(L),Scf%xcoef(k)%b(L), L=1,4), Scf%xcoef(k)%c, &
                           Scf%delta_fp(k), Scf%delta_fpp(k)
          end do
          write(unit=lun,fmt="(/,/)")
       end if

    End Subroutine Set_Form_Factors

    Subroutine Allocate_Scattering_Species(n,Scf)
      integer,                       intent(in)  :: n
      type(Scattering_Species_Type), intent(out) :: Scf
      integer :: i
      Scf%Num_Species=n
      allocate(Scf%br(n),Scf%bi(n),Scf%delta_fp(n),Scf%delta_fpp(n),Scf%symb(n))
      Scf%br=0.0; Scf%bi=0.0; Scf%delta_fp=0.0; Scf%delta_fpp=0.0; Scf%symb= " "
      allocate(Scf%Xcoef(n))
      do i=1,Scf%Num_Species
        Scf%Xcoef(i)%Z=0
        Scf%Xcoef(i)%a=0.0
        Scf%Xcoef(i)%b=0.0
        Scf%Xcoef(i)%c=0.0
      end do
      return
    End Subroutine Allocate_Scattering_Species

    Subroutine Mpol_Ffactor(coeff,ac,sn,mff)
       real(kind=cp), dimension(4), intent(in) :: coeff
       real(kind=cp), dimension(:), intent(in) :: ac
       real(kind=cp), intent(in) :: sn
       real(kind=cp), intent(out):: mff

       ! Local variables
       integer :: i
       real(kind=cp)    :: j0,j2,j4,j6
       j0=ac(7); j2=ac(14); j4=ac(21); j6=ac(28)
       do i=1,5,2
         j0=j0+ac(i)*EXP(-ac(i+1)*sn)      !<j0>      value for Q=H+k
         j2=j2+ac(i+7)*EXP(-ac(i+8)*sn)    !<j2>/sn value for Q=H+k
         j4=j4+ac(i+14)*EXP(-ac(i+15)*sn)  !<j4>/sn value for Q=H+k
         j6=j6+ac(i+21)*EXP(-ac(i+22)*sn)  !<j6>/sn value for Q=H+k
       end do
       j2=j2*sn  !<j2> value for Q=H+k
       j4=j4*sn  !<j4> value for Q=H+k
       j6=j6*sn  !<j6> value for Q=H+k
       mff= pn*coeff(1)*j0 + coeff(2)*j2 + coeff(3)*j4 + coeff(4)*j6
       return
    End Subroutine Mpol_Ffactor

    Function Schwinger_Amplitude(hn,pol,theta,fe,Left,UB)  result(Schwinger)
      real(kind=cp), dimension(3),            intent(in) :: hn,pol
      real(kind=cp),                          intent(in) :: theta !in radians
      complex(kind=cp),                       intent(in) :: fe
      logical,                      optional, intent(in) :: left
      real(kind=cp), dimension(3,3),optional, intent(in) :: UB
      complex(kind=cp)                                   :: Schwinger

      real(kind=cp),dimension(3) :: uvect, hc

      if(present(left)) then
        if(left) then
          uvect=(/0.0,0.0,-1.0/)
        else
          uvect=(/0.0,0.0,1.0/)
        end if
      else if(present(UB)) then
        hc=matmul(UB,hn)
        uvect=cross_product((/0.0,1.0,0.0/),hc)
        uvect=-uvect/sqrt(dot_product(uvect,uvect))
      else
        Schwinger=0.0
        return
      end if
      !change of sign of "i" because the convention used by crystallographers is Q=Kf-Ki
      Schwinger= schw*cmplx(0.0,-1.0)*dot_product(pol,uvect)*fe/abs(tan(theta))
    End Function Schwinger_Amplitude

    Subroutine Calc_General_StrFactor(hn,sn,Atm,Grp,Scf,fn,fx,fe)
       !---- Arguments ----!
       real(kind=cp),dimension(3),         intent(in) :: hn
       real(kind=cp),                      intent(in) :: sn !(sinTheta/Lambda)**2
       type(atom_list_type),               intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp
       type(Scattering_Species_Type),      intent(in) :: Scf
       complex,                            intent(out):: fn,fx,fe

       !---- Local Variables ----!
       integer                               :: i,j,k,m
       real(kind=cp)                         :: arg,anis,scosr,ssinr,b,s
       real(kind=cp)                         :: a1,a3,b1,b3,av,bv,nffr,nffi        !fn
       real(kind=cp)                         :: xa1,xa3,xb1,xb3,xav,xbv,xffr,xffi  !fx
       real(kind=cp)                         :: ea1,ea3,eb1,eb3,eav,ebv,effr       !fe
       real(kind=cp),dimension(3)            :: h
       real(kind=cp),dimension(6)            :: beta

       !--- Initialising local variables
       a1 =0.0;  a3 =0.0
       b1 =0.0;  b3 =0.0
       ea1=0.0;  ea3=0.0
       eb1=0.0;  eb3=0.0
       xa1=0.0;  xa3=0.0
       xb1=0.0;  xb3=0.0

       do i=1,Atm%natoms
          arg=0.0
          scosr=0.0
          ssinr=0.0
          do k=1,grp%NumOps
             h=matmul(hn,grp%Symop(k)%Rot)
             arg=tpi*(dot_product(h,Atm%atom(i)%x)+ dot_product(hn,grp%Symop(k)%tr))
             anis=1.0
             if(Atm%atom(i)%thtype == "aniso") then
               beta=Atm%atom(i)%u(1:6)
               anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                    +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
               anis=exp(-anis)
             end if
             scosr=scosr+COS(arg)*anis  !Real part of geometrical structure factor for the current atom
             ssinr=ssinr+SIN(arg)*anis  !Imaginary part of geometrical structure factor for the current atom
          end do ! symmetry

          b= atm%atom(i)%occ * exp(-atm%atom(i)%biso*sn)
          !Calculation of scattering factors
          j=atm%atom(i)%ind(1)  !pointer to the form factor coefficients

          nffr = Scf%br(j)*b
          nffi = Scf%bi(j)*b

          xffr=Scf%Xcoef(j)%c
          do k=1,4
            xffr=xffr+Scf%Xcoef(j)%a(k)*exp(-Scf%Xcoef(j)%b(k)*sn)
          end do

          effr = (real(Scf%Xcoef(j)%Z)-xffr)*b  !<- Here delta_fp is not used ....
          !write(*,"(tr4,a,i4,3f10.4)")  Atm%atom(i)%chemsymb//": Z, xffr, b, effr", Scf%Xcoef(j)%Z, xffr, b, effr
          xffr = (xffr+Scf%delta_fp(j))*b       ! (f0+Deltaf')*OCC*Tiso
          xffi = Scf%delta_fpp(j)*b             !     Deltaf" *OCC*Tiso

          a1 = a1 + nffr*scosr  ! F=A+iB: components of A  and B (ai,bi)
          b1 = b1 + nffi*scosr  ! a2,b2,a4,b4 are components for anisotropic form factors
          a3 = a3 + nffi*ssinr  ! they are not used here
          b3 = b3 + nffr*ssinr  ! For general case: av = a1-a2-a3-a4, bv = b1-b2+b3+b4

          xa1 = xa1 + xffr*scosr  ! F=A+iB: components of A  and B (ai,bi)
          xb1 = xb1 + xffi*scosr  ! a2,b2,a4,b4 are components for anisotropic form factors
          xa3 = xa3 + xffi*ssinr  ! they are not used here
          xb3 = xb3 + xffr*ssinr  ! For general case: av = a1-a2-a3-a4, bv = b1-b2+b3+b4

          ea1 = ea1 + effr*scosr  ! No anomalous imaginary component is used here
          eb3 = eb3 + effr*ssinr  ! there is no anomalous scattering neutron + electrons
          !write(*,"(tr4,a,tr4,3f10.4)") Atm%atom(i)%chemsymb//": effr,   ea1, eb3",  effr, ea1, eb3

       end do ! Atoms

       av = a1-a3    !real part of the Nuclear structure factor
       bv = b1+b3    !imaginary part of the Nuclear structure factor
       fn=cmplx(av,bv) * Grp%Centred * Grp%NumLat

       xav = xa1-xa3    !real part of the X-rays structure factor
       xbv = xb1+xb3    !imaginary part of the X-rays structure factor
       fx=cmplx(xav,xbv) * Grp%Centred * Grp%NumLat

       fe=cmplx(ea1,eb3) * Grp%Centred * Grp%NumLat

    End Subroutine Calc_General_StrFactor


End Program Schwinger

