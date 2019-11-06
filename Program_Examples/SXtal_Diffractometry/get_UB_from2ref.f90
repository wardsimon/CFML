!!----
!!---- Program: get_UB_from2ref
!!----          Example of simple program using CFML, testing CFML_Geometry_SXTAL subroutines
!!----
!!---- Author: Juan Rodriguez-Carvajal
!!---- Revision: Nov-2019
!!
Program get_UB_from2ref
   !---- Use Modules ----!
   use CFML_GlobalDeps,      only: cp
   use CFML_Crystal_Metrics, only: Crystal_Cell_Type,set_crystal_cell,Write_Crystal_Cell
   use CFML_Geometry_SXTAL,  only: z1frfc,GenUB,Normal_Beam_Angles,Get_Angs_NB
   use CFML_String_Utilities,only: L_case

   !---- Variables ----!
   implicit none

   character(len=80)          :: name_file
   character(len=132)         :: line
   type (Crystal_Cell_Type)   :: cell
   integer                    :: i,narg , ier
   real,dimension(3)          :: celda, angulo,h1,h2,z1,z2
   real,dimension(4)          :: ang1,ang2
   real,dimension(3,3)        :: UB
   real                       :: lambda


   !---- Initializing ----!
   name_file=" "
   narg=COMMAND_ARGUMENT_COUNT()
   write(unit=*,fmt="(/,a)") "    ==============================================="
   write(unit=*,fmt="(a)")   "        Calculation of UB from two reflections     "
   write(unit=*,fmt="(a)")   "           Calculation of motors given UB          "
   write(unit=*,fmt="(a)")   "    ==============================================="
   write(unit=*,fmt="(a)")   "  A CFL file with cell parameters, wavelength and instrument"
   write(unit=*,fmt="(a)")   "  file is needed as well as the instrument file *.geom"
   if(narg > 0) then
      call GET_COMMAND_ARGUMENT(1,name_file)
   Else
      write(unit=*,fmt="(a)",advance="no") " => Enter the name of the input CFL-file: "
      read(unit=*,fmt="(a)") name_file
   end if
   if (len_trim(name_file)==0) Call Abort_Program()

   call Read_CFL()  !Internal subroutine to read the input file
   call set_crystal_cell(celda,angulo,cell)


   !---- Procedure ----!
   do
      !> Set the range for HKL calculation
      write(unit=*,fmt=*) " "
      write(unit=*,fmt="(a)",advance="no") " => Give hkl (3 reals) and motors (4 reals) of the  first reflection : "
      read(unit=*,fmt=*,iostat=ier) h1,ang1
      if(ier /= 0) cycle
      if (sum(abs(h1)) < 0.001) exit
      write(unit=*,fmt="(a)",advance="no") " => Give hkl (3 reals) and motors (4 reals) of the second reflection : "
      read(unit=*,fmt=*,iostat=ier) h2,ang2
      call z1frfc(lambda,ang1(1),ang1(2),ang1(3),ang1(4),z1)
      call z1frfc(lambda,ang2(1),ang2(2),ang2(3),ang2(4),z2)
      call GenUB(cell%BL_M,h1,h2,z1,z2,UB,ier)
      if(ier /= 0) then
        write(unit=*,fmt="(a)") " => UB-matrix cannot be calculated with the provided reflections!"
        cycle
      else
        write(unit=*,fmt="(a)") " => UB-matrix:"
        do i=1,3
          write(unit=*,fmt="(3f14.5)") UB(i,:)
        end do

      end if

      !Calculation of motors for other reflections
      do
        write(unit=*,fmt="(a)",advance="no") " => Give hkl for calculating motors: "
        read(unit=*,fmt=*,iostat=ier) h1
        if(ier /= 0) exit
        if(sum(abs(h1)) < 0.001) exit
        i=1
        call Normal_Beam_Angles(lambda,ub,h1,i,ang1,ier)
        if(ier /= 0) then
          write(unit=*,fmt="(a)") " => Normal-beam angles cannot be calculated"
          cycle
        else !ANBCAL(1:4) -> GAMMA, OMEGA(NB), NU, THETA
          write(unit=*,fmt="(a,4f14.5)") " => Normal-beam: Gamma, omega, nu, theta = ",ang1
        end if

        z1=matmul(UB,h1)
        ang1(4) = asind(0.5*lambda*sqrt(dot_product(z1,z1)))
        call Get_Angs_NB(lambda,z1,ang1(1),ang1(2),ang1(3),ier)
        if(ier /= 0) then
          write(unit=*,fmt="(a)") " => Normal-beam angles cannot be calculated"
          cycle
        else !ANBCAL(1:4) -> GAMMA, OMEGA(NB), NU, THETA
          write(unit=*,fmt="(a,4f14.5)") " => Normal-beam: Gamma, omega, nu, theta = ",ang1
        end if

        !Four circle angles
        call calang(h1,ang1(1),ang1(2),ang1(3),ang1(4),ier,lambda,ub,1)
        if(ier /= 0) then
          write(unit=*,fmt="(a)") " => Four circle angles cannot be calculated for Geom=1"
        else !Gamma,omega,chi, phi
          write(unit=*,fmt="(a,4f14.5)") " => Four circle (Geom=1): Gamma, omega, chi, phi = ",ang1(1:4)
        end if
        call calang(h1,ang1(1),ang1(2),ang1(3),ang1(4),ier,lambda,ub,2)
        if(ier /= 0) then
          write(unit=*,fmt="(a)") " => Four circle angles cannot be calculated for Geom=2"
        else !Gamma,omega,chi, phi
          write(unit=*,fmt="(a,4f14.5)") " => Four circle (Geom=2): Gamma, omega, chi, phi = ",ang1(1:4)
        end if
        call calang(h1,ang1(1),ang1(2),ang1(3),ang1(4),ier,lambda,ub,3)
        if(ier /= 0) then
          write(unit=*,fmt="(a)") " => Four circle angles cannot be calculated for Geom=3"
        else !Gamma,omega,chi, phi
          write(unit=*,fmt="(a,4f14.5)") " => Four circle (Geom=3): Gamma, omega, chi, phi = ",ang1(1:4)
        end if
        call calang(h1,ang1(1),ang1(2),ang1(3),ang1(4),ier,lambda,ub,4)
        if(ier /= 0) then
          write(unit=*,fmt="(a)") " => Four circle angles cannot be calculated for Geom=4"
          cycle
        else !Gamma,omega,chi, phi
          write(unit=*,fmt="(a,4f14.5)") " => Four circle (Geom=4): Gamma, omega, chi, phi = ",ang1(1:4)
        end if
      end do
   end do

   write(unit=*,fmt="(a)") " => Program finished ... "

   contains

     Subroutine Read_CFL()
       logical :: cell_read, lambda_read
       character(len=132) :: lline
       open(unit=1,file=trim(name_file),status="old",action="read",position="rewind",iostat=ier)
       cell_read=.false.; lambda_read=.false.
       if(ier /= 0) then
         write(unit=*,fmt="(a)") " => Error opening the file: "//trim(name_file)
         Call Abort_Program()
       end if
       do
         read(unit=1,fmt="(a)",iostat=ier) line
         if(ier /= 0) exit
         if(len_trim(line) == 0) cycle
         line=adjustl(line)
         if(line(1:1) == "!" .or. line(1:1) == "#" ) cycle
         lline=adjustl(L_case(line))


         if(lline(1:6) == "lambda" ) then
           read(unit=line(7:),fmt=*,iostat=ier) lambda !wavelength
           if(ier /= 0) then
             write(unit=*,fmt="(a)") " => Error reading the wavelength "
             Call Abort_Program()
           end if
           lambda_read=.true.
           cycle
         end if

         if(lline(1:4) == "cell" ) then
           read(unit=line(6:),fmt=*,iostat=ier) celda,angulo !cell parameters
           if(ier /= 0) then
             write(unit=*,fmt="(a)") " => Error reading the cell parameters "
             Call Abort_Program()
           end if
           cell_read=.true.
           cycle
         end if
       end do

       if(.not. cell_read) then
        write(unit=*,fmt="(a)") " => No unit cell has been given!"
        Call Abort_Program()
       end if
       if(.not. lambda_read) then
        write(unit=*,fmt="(a)") " => No Lambda has been given!"
        Call Abort_Program()
       end if
       close(unit=1)
       return
     End Subroutine Read_CFL

     Subroutine Abort_Program()
       write(unit=*,fmt="(a)") " => Program aborted!"
       stop
     End Subroutine Abort_Program

     !!---- Subroutine calang(h,tteta,om,ch,ph,ierr,wave,ub,igeom)
     !!---- This is a copy of the subroutine existing in CFML_Geometry_SXTAL module
     !!---- in which the dependency on instrument has been removed.
     Subroutine calang(h,tteta,om,ch,ph,ierr,wave,ub,igeom)
        !---- Arguments ----!
        real(kind=cp),Dimension(3),    Intent(In) :: h
        real(kind=cp),                 Intent(Out):: tteta,om,ch,ph
        Integer,                       Intent(Out):: ierr
        real(kind=cp),                 intent(in) :: wave
        real(kind=cp), dimension(3,3), intent(in) :: ub
        integer,                       intent(in) :: igeom

        !--- Local Variables ---!
        real(kind=cp)                 :: ttmax,ttmin,sint,chmax,diff,ds,theta
        real(kind=cp), Dimension(3)   :: z1

        if(igeom == 2) then  !high-chi
           chmax= 200.0
        else
           chmax= 170.0
        end if

        ttmin= 0.0
        ttmax= 180.0
        ierr=0

        z1=Matmul(ub,h)
        sint=0.5*wave*Sqrt(Dot_Product(z1,z1))
        If (abs(sint) > 1.0) Then
           ierr=1
           Return
        End If
        theta =asind(sint)
        om=theta       !Theta = omega in bisecting geometry
        tteta=2.0*theta

        If (tteta < ttmin .or. tteta > ttmax) Then
           ierr=1
           Return
        End If

        ds = Sqrt(z1(1)*z1(1)+z1(2)*z1(2))
        Select Case(igeom)
           Case(1,2)    !1:Bisecting Geometry  (PSI=0), 2:Bisecting Geometry - High CHI
              ph = atan2d(z1(2),z1(1))  !Eqs 38 BL
              ch = atan2d (z1(3),ds)
              If(igeom == 2) Then
                ch=180.0-ch  !Eqs 40 BL
                ph=180.0+ph
              End If

              If(ph >  180.0) ph=ph-360.0
              If(ph < -180.0) ph=ph+360.0
              If(ch >  180.0) ch=ch-360.0
              If(ch < -180.0) ch=ch+360.0
              If (chmax > 180.0) Then
                 diff=chmax-180.0
                 If (ch < -180.0 + diff .AND. ch >= -180.0)  ch=ch+360.0
              End If

           Case(3)    !3:  Geometry NORMAL BEAM
              ch=-90.0
              ph=-atan2d(z1(1),z1(2))
              om=theta+90.0-atan2d(z1(3),ds)

           Case(4)    !4:  Geometry parallel (PSI=90) ..D15-D16  (PSI=90)
              ch=+90.0
              ph=atan2d(z1(1),-z1(2))     ! Eqs 41 BL (with different definition of omega
              om=theta+atan2d(-ds,z1(3))  ! Omega=Theta+atan2d(-ds,z1(3))

        End Select

        Return
     End Subroutine calang

End Program get_UB_from2ref
