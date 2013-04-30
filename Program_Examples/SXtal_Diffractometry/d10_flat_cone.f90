  !!---- Below the program using the FlatCone modules -------
  !!---- Program written by JRC the 28 February 2013 (Just 45 years after an earthquake at Sevilla!)
  !!----
  !!---- The program is just to test the continuity of rotation angles
  !!---- when turning around a zone axis put perpendicular to the detector plane
  !!---- It detects the integer reflections apearing in the detector during
  !!---- the Psi-scan.

  Program D10_Flat_Cone
   use CFML_GlobalDeps,      only: cp
   use CFML_Math_General,    only: asind,sind,cosd,atan2d
   use CFML_Math_3D,         only: Invert=>Invert_A, Cross_Product
   use CFML_Crystal_Metrics, only: Crystal_Cell_Type, Cart_Vector,Rot_Matrix, &
                                   ERR_Crys_Mess,ERR_Crys, set_crystal_cell
   use CFML_String_Utilities,only: L_case
   use CFML_Geometry_Calc,   only: Set_Rotation_Matrix
   use CFML_Geometry_SXTAL,  only: z4frgn, z1frz4, Get_FlatCone_Angles_D10

   !Adding new utilities for calculating the structure factor of reflections
   use CFML_crystallographic_symmetry,only: space_group_type
   use CFML_Atom_TypeDef,             only: Atom_List_Type
   use CFML_IO_Formats,               only: Readn_set_Xtal_Structure,err_form_mess,err_form
   use CFML_Structure_Factors,        only: Calc_hkl_StrFactor,Init_Calc_hkl_StrFactors
   use CFML_Reflections_Utilities,    only: Hkl_Absent

   implicit none

   Character(len=150)            :: cfl_file, fileout, line
   integer                       :: i,j,k,n,ier,narg, npoints,nref
   real(kind=cp)                 :: Lambda, ruvw, dstar,omega,chi,phi,psi,alpha, &
                                    gamma,nu,mu,step,psi1,psi2, sn, fsq
   real(kind=cp), dimension(3)   :: h1,h2,uvw,r1,z1,z2,z3,z4,dL
   real(kind=cp), dimension(3,3) :: U_Mat,UB,UB_inv,UB_invt, R_psi, Gibbs,Rot, &
                                    M_Omega,M_Chi,M_Phi,M_set
   real(kind=cp), dimension(2,3) :: limits
   type (Crystal_Cell_Type)      :: Cell
   type (space_group_type)       :: SpG
   type (Atom_list_Type)         :: Atm
   logical                       :: ok,new_ref
   character(len=80)             :: mess
   integer,       parameter      :: maxang=5000, nfil=128, maxref=500
   logical,       dimension(maxang) :: inlim
   real(kind=cp), dimension(maxang) :: om,ch,ph,ps
   real(kind=cp), dimension(maxref) :: om_keep, ch_keep, ph_keep, stl2
   integer,       dimension(maxref) :: pix
   ! Flat cone detector
   real(kind=cp), parameter                :: span_ang=120.0_cp
   real(kind=cp), dimension(3,nfil,maxang) :: q_hkl
   real(kind=cp), dimension(3,nfil)        :: zfc
   integer,       dimension(3,maxref)      :: hkl
   integer,       dimension(3)             :: i_uvw, ihkl
   logical                                 :: UBM_read=.false.

   write(unit=*,fmt="(/a)")"---------------------------------------------"
   write(unit=*,fmt="(a)") "  PROGRAM to simulate FlatCone motions on D10"
   write(unit=*,fmt="(a)") "  28 February 2013, JRC"
   write(unit=*,fmt="(a/)")"---------------------------------------------"

   narg=COMMAND_ARGUMENT_COUNT()
   if(narg > 0) then
      call GET_COMMAND_ARGUMENT(1,cfl_file)
   Else
      write(unit=*,fmt="(a)") " => Error! An input *.cfl file should be provided!"
      Call Abort_Program()
   end if

   ! Read the cell parameters (and eventually the space group and atoms), lambda and
   ! the orientation of the crystal from the input file my_name.cfl.
   ! The derived types Cell, SpG, Atm are set in that subroutine

   call Read_CFL()  !Internal subroutine to read the input file
   write(unit=*,fmt="(a)") " => CFL-file read successfully!"

   !Calculate the UB-matrix from given orientation angles if no UB-matrix has been read
   if(.not. UBM_Read) UB=Matmul(U_Mat,Cell%BL_m)
   UB_inv=Invert(UB)
   UB_invt=transpose(UB_inv)

   !Initialize module of structure factors
   if(Atm%Natoms > 0) then
     call Init_Calc_hkl_StrFactors(Atm,"NUC",lambda)
   end if

   do
     write(unit=*,fmt="(/a)",advance="no") &
     " => Enter the reciprocal plane (two reflections) or the zone axis to be explored (<cr> to stop): "
     read(unit=*,fmt="(a)") line
     if(len_trim(line) == 0) exit
     read(unit=line,fmt=*,iostat=ier) h1,h2
     if(ier /= 0) then   !Try to read directly the zone axis
       read(unit=line,fmt=*,iostat=ier) uvw
     else
       uvw=Cross_Product(h1,h2)
     end if
     i_uvw=nint(uvw)
     z1=Cart_Vector("BLD",uvw,Cell)   !Coordinates in the Busing-Levy Crystal system
     ruvw=sqrt(dot_Product(z1,z1))
     z1=matmul(UB_invt,uvw)             !Cartesian coordinates in reciprocal space (L-system)
     r1=z1/sqrt(dot_product(z1,z1))     !Unitary vector along the zone-axis when all angles are set to zero
     dstar=1.0/ruvw
     write(unit=*,fmt="(a,3f8.4,a,f8.4,a)") " => The zone axis is: (",uvw,") of module: ",ruvw," angstroms"
     write(unit=*,fmt="(a)",advance="no") " => Enter the level number to be explored: "
     read(unit=*,fmt=*,iostat=ier) n
     if(ier /= 0) exit
     !Calculate the mu-angle
     mu=n*dstar*lambda
     if(abs(mu) > 1.0_cp) then
       write(unit=*,fmt="(a)") " => Level not accessible !!!"
       cycle
     end if
     mu=asind(mu)
     dL=(/0.0_cp, sind(mu), cosd(mu)/)
     write(unit=*,fmt="(a,3f10.5,a)") " => The unitary vector r1 is: (",r1,")"
     write(unit=*,fmt="(a,3f10.5,a)") " => The unitary vector dL is: (",dL,")"
     write(unit=*,fmt="(a,f8.4,a)")   " => The mu-angle  is: ",mu," degrees"

     !We can calculate here the angles (gamma,nu) of each pixel in the detector
     !using the expressions:
     !      sin(nu)=-sin(alpha)sin(mu)
     !   tan(gamma)=cos(alpha)/cos(mu)/sin(alpha)
     ! See document for the definition of alpha.
     !
     ! The corresponding (fixed for given mu) diffraction vectors z4 can be obtained by
     ! using the subroutine z4frgn(wave,ga,nu,z4) for each pixel (get z4 from gamma and nu)
     !
     step=span_ang/real(nfil-1)
     write(unit=*,fmt="(/,a)") "     Pixel    Alpha    Gamma      Nu                z4(x)     Z4(y)     Z4(z)"
     do i=1,nfil
      alpha= 90.0-step*real(i-1)  !Assume that the first pixed is at alpha=90
      nu=asind(-sind(alpha)*sind(mu))
      gamma=atan2d(cosd(alpha),cosd(mu)*sind(alpha))
      call z4frgn(lambda,gamma,nu,zfc(:,i))
      write(*,"(i10,3f9.4,tr10,3f10.5)")  i, alpha,gamma,nu, zfc(:,i)
     end do

     ! Calculate the setting angles omega, chi, phi

     call Get_FlatCone_Angles_D10(r1,mu,psi1,psi2,npoints,limits,ps,om,ch,ph,inlim, Gibbs)

     write(unit=*,fmt="(/a)") " ---------------------------------------"
     write(unit=*,fmt="( a)") "     Psi      Omega      Chi       Phi       z4=Rot.r1 (Should be always the same)"
     write(unit=*,fmt="( a)") " ---------------------------------------"

     ok=.false.
     nref=0
     do i=1,npoints
        if(.not. inlim(i)) cycle  ! Output only the accessible range
        ok=.true.
        psi=ps(i)
        R_psi= Rot_Matrix(dL,psi)
        !Total rotation matrix
        Rot=Matmul(R_psi,Gibbs)  !This matrix changes a vector of the L-system when
                                 !all angles are zero to a position in which the zone
                                 !axis is put parallel to the Z-axis of detector system and
                                 !then a rotation of angle Psi is done around the zone axis
        !Calculate here the list reciprocal coordinates of the recorded data on the detector
        !From of each z4 of the pixels get z1, using z1frz4(z4,om,ch,ph,z1),
        !get finally hkl or each pixel from hkl=inv(UB)z1

        !Test of the matrix               Done and everything seems OK
        !call Chi_mat(ch(i),M_Chi)
        !call Phi_mat(ph(i),M_Phi)
        !call Phi_mat(om(i),M_Omega)
        !M_set=Matmul(M_Omega,Matmul(M_Chi,M_Phi))
        !M_set=Matrix_OmegaChiPhi(om(i),ch(i),ph(i),"D")
        !z3=Matmul(M_set,r1)
        !z3=Matmul(transpose(M_set),r1)

        do j=1,nfil
          call z1frz4(zfc(:,j),om(i),ch(i),ph(i),z1)
          sn=0.25*dot_product(z1,z1)
          z1=Matmul(UB_inv,z1)
          z2=abs(z1-real(nint(z1)))
          if(sum(z2) < 0.05) then
            ihkl=nint(z1)
            !Test if it is a forbidden reflection
            if(.not. Hkl_Absent(ihkl,SpG)) then
              new_ref=.true.
              do k=1,nref
                if(sum(ihkl-hkl(:,k)) == 0) then
                  new_ref=.false.
                  exit
                end if
              end do
              if(new_ref) then
                nref=nref+1
                hkl(:,nref)=ihkl
                om_keep(nref)=om(i)
                ch_keep(nref)=ch(i)
                ph_keep(nref)=ph(i)
                pix(nref)    = j
                stl2(nref)   = sn
              end if
            end if
          end if
          q_hkl(:,j,i)=z1
        end do
        z2=Matmul(Rot,r1)  !This should be always equal to dL
        write(unit=*,fmt="(4f10.4,tr5,a,3f8.4,a)") Psi, om(i),ch(i),ph(i),"(",z2,")"
       ! write(unit=*,fmt="(4f10.4,tr5,4(a,3f8.4))") Psi, om(i),ch(i),ph(i),"(",z2, &
       ! ")  -> Q-ini=(",q_hkl(:,1,i),")   Q-cent=(",q_hkl(:,nfil/2,i),")   Q-fin=(",q_hkl(:,nfil,i)
     end do
     if(.not. ok)  write(unit=*,fmt="(a)")  " => Reciprocal plane not accessible!"
     if(nref > 0) then
       write(unit=*,fmt="(/,a,3i4,a,i3)")        " => List of accessible reflections around zone axis [",i_uvw,"] level number ",n
       write(unit=*,fmt="(a,f8.4,a)")   " => The mu-angle  is: ",mu," degrees"
       write(unit=*,fmt="(/,a)") " Nref             H   K   L    h.uvw        Omega     Chi     Phi      Pixel       F^2"
       if(Atm%Natoms > 0) then
         do i=1,nref
           call Calc_hkl_StrFactor("SXtal","Neutrons",hkl(:,i),stl2(i),Atm,SpG,fsq)
           write(unit=*,fmt="(i5,tr10,3i4,tr5,i4,tr5,3f9.3,i8,f12.4)") i, hkl(:,i), dot_Product(hkl(:,i),i_uvw), &
                 om_keep(i),ch_keep(i),ph_keep(i),pix(i),fsq
         end do
       else
         do i=1,nref
           write(unit=*,fmt="(i5,tr10,3i4,tr5,i4,tr5,3f9.3,i8)") i, hkl(:,i), dot_Product(hkl(:,i),i_uvw), &
                 om_keep(i),ch_keep(i),ph_keep(i),pix(i)
         end do
       end if
     else
       write(unit=*,fmt="(a)") " => No integer reflection accessible"
     end if
   end do
   write(unit=*,fmt="(a)")  " => Program terminated normally ..."

   contains

     Subroutine Abort_Program()
       write(unit=*,fmt="(a)") " => Program aborted!"
       stop
     End Subroutine Abort_Program

     Subroutine Read_CFL()
       logical :: cell_read, lambda_read, orient_read

       !With this call Cell,SpG and Atm are set
       call Readn_set_Xtal_Structure(trim(cfl_file),Cell,SpG,Atm,Mode="CFL")
       if(err_form) then
        write(unit=*,fmt="(a)") "  "//trim(err_form_mess)
        Call Abort_Program()
       end if

       ! Re-open the CFL file to read the rest of relevant items
       open(unit=1,file=trim(cfl_file),status="old",action="read",position="rewind",iostat=ier)
       if(ier /= 0) then
         write(unit=*,fmt="(a)") " => Error opening the file: "//trim(cfl_file)
         stop
       end if
       cell_read=.false.; lambda_read=.false.; orient_read=.false.
       limits(1,1) = -15.0_cp;  limits(2,1) = 55.0_cp   !default omega limits
       limits(1,2) = -15.0_cp;  limits(2,2) = 26.0_cp   !default chi limits
       limits(1,3) =-179.99_cp; limits(2,3) = 180.0_cp  !default phi limits
       psi1=0.0; psi2=360.0; npoints=360  !default psi-range

       do
         read(unit=1,fmt="(a)",iostat=ier) line
         if(ier /= 0) exit
         if(len_trim(line) == 0) cycle
         line=adjustl(line)
         if(line(1:1) == "!" .or. line(1:1) == "#" ) cycle
         line=adjustl(L_case(line))

         !if(line(1:4) == "cell" ) then
         !  read(unit=line(5:),fmt=*,iostat=ier) h1,h2 !six numbers h1=(a,b,c), h2=(alpha,beta,gamma)
         !  if(ier /= 0) then
         !    write(unit=*,fmt="(a)") " => Error reading the cell parameters "
         !    Call Abort_Program()
         !  end if
         !  Call Set_Crystal_Cell(h1,h2,Cell) !Setting the derived type Cell
         !  if(ERR_Crys) then
         !    write(unit=*,fmt="(a)") " => Error: "//trim(ERR_Crys_Mess)
         !    Call Abort_Program()
         !  else
         !    cell_read=.true.
         !  end if
         !  cycle
         !end if

         if(line(1:9) == "psi_range" ) then
           read(unit=line(10:),fmt=*,iostat=ier) psi1,psi2,npoints !Three angles h1=(phi_x,phi_y,phi_z)  to construct the U-matrix
           if(ier /= 0) then
             psi1=0.0; psi2=360.0; npoints=360
           end if
           cycle
         end if

         if(line(1:6) == "orient" ) then
           read(unit=line(7:),fmt=*,iostat=ier) h1 !Three angles h1=(phi_x,phi_y,phi_z)  to construct the U-matrix
           if(ier /= 0) then
             write(unit=*,fmt="(a)") " => Error reading the orientation angles "
             Call Abort_Program()
           end if
           Call Set_Rotation_Matrix(h1,U_Mat) !Setting the U-Matrix
           orient_read=.true.
           cycle
         end if

         if(line(1:6) == "ub_mat" ) then
           do i=1,3
             read(unit=1,fmt=*,iostat=ier) UB(i,:)
           end do
           if(ier /= 0) then
             write(unit=*,fmt="(a)") " => Error reading the UB-matrix "
             Call Abort_Program()
           end if
           ubm_read=.true.
           cycle
         end if

         if(line(1:6) == "lambda" ) then
           read(unit=line(7:),fmt=*,iostat=ier) lambda !wavelength
           if(ier /= 0) then
             write(unit=*,fmt="(a)") " => Error reading the wavelength "
             Call Abort_Program()
           end if
           lambda_read=.true.
           cycle
         end if

         if(line(1:12) == "omega_limits" ) then
           read(unit=line(13:),fmt=*,iostat=ier) limits(1,1),limits(2,1)
           if(ier /= 0) then
             write(unit=*,fmt="(a)") " => Error reading the omega limits "
             Call Abort_Program()
           end if
           cycle
         end if

         if(line(1:10) == "chi_limits" ) then
           read(unit=line(11:),fmt=*,iostat=ier) limits(1,2),limits(2,2)
           if(ier /= 0) then
             write(unit=*,fmt="(a)") " => Error reading the chi limits "
             Call Abort_Program()
           end if
           cycle
         end if

       end do

       !if(.not. cell_read) then
       ! write(unit=*,fmt="(a)") " => No unit cell has been given!"
       ! Call Abort_Program()
       !end if
       if(.not. lambda_read) then
        write(unit=*,fmt="(a)") " => No Lambda has been given!"
        Call Abort_Program()
       end if
       if(.not. orient_read .and. .not. ubm_read) then
        write(unit=*,fmt="(a)") " => No orientation has been given!"
        Call Abort_Program()
       end if
       close(unit=1)
       return
     End Subroutine Read_CFL

  End Program D10_Flat_Cone