   module moment_mod

     implicit none

     integer, parameter ::  natom=50,nv=24,ncel=1000
     real,    parameter ::  rad=57.29577951,tpi=6.28318531, pi=tpi/2.0
     integer            ::  ifou, nat, nvk, nce
     integer, parameter ::  inp=1,iou=2,iom=3    !logical units
     integer            ::  ier                  !iostat control value

     real, dimension(5,natom)     :: x
     real, dimension(3,nv,natom)  :: rr, ri
     real, dimension(nv,natom)    :: ph
     real, dimension(nv)          :: phk

     real, dimension(3,natom,ncel):: xcry,rmo,xcar,rmc
     real, dimension(  natom,ncel):: rmod

     real, dimension(nv,natom)    :: phko
     real, dimension(3,nv)        :: pvk
     real, dimension(nv)          :: phkop


      character(len=4), dimension(natom) ::  nam
      character(len=4), dimension(2)     ::  nfst
      character(len=20)                  ::  symb
      character(len=80)                  ::  title

      real :: asq,bsq,csq,ab,ac,bc, astar,bstar,cstar,castar, &
              cbstar,ccstar,sastar, sbstar,scstar,a,b,c,ca,cb,cc

      real, dimension(3,3) :: ortd, orti

      logical              ::  phasgiven,momencal
      real, dimension(3,4) ::  ltr
      integer :: nlat,inlat

      character(len=*), parameter, dimension(7)::  latt =             &
       (/ "  P: { 000 }                                       ",  &
          "  A: { 000;  0  1/2 1/2 }+                         ",  &
          "  B: { 000; 1/2  0  1/2 }+                         ",  &
          "  C: { 000; 1/2 1/2  0  }+                         ",  &
          "  I: { 000; 1/2 1/2 1/2 }+                         ",  &
          "  R: { 000; 2/3 1/3 1/3; 1/3 2/3 2/3  }+           ",  &
          "  F: { 000;  0  1/2 1/2; 1/2  0  1/2; 1/2 1/2  0 }+"/)


   contains

      subroutine latsym(symb)
!---------------------------------------------------------------------
!     SYMB is the space group symbol in International notation, only the
!     first character is taken into account for lattice identification.
!     This subroutine constructs the lattice centring vectors.
!---------------------------------------------------------------------
         character (len=*), intent(in)   :: symb
         character (len=1)               :: lat
         lat=symb(1:1)
         nlat=1
         ltr(1,1)=0.0
         ltr(2,1)=0.0
         ltr(3,1)=0.0
         if(lat == "P" .OR. lat == "p") then
           nlat=1
           inlat=1
         else if(lat == "A" .OR. lat == "a") then
           nlat=2
           inlat=2
           ltr(1,2)=0.0
           ltr(2,2)=0.5
           ltr(3,2)=0.5
         else if(lat == "B" .OR. lat == "b") then
           nlat=2
           inlat=3
           ltr(1,2)=0.5
           ltr(2,2)=0.0
           ltr(3,2)=0.5
         else if(lat == "C" .OR. lat == "c") then
           nlat=2
           inlat=4
           ltr(1,2)=0.5
           ltr(2,2)=0.5
           ltr(3,2)=0.0
         else if(lat == "I" .OR. lat == "i") then
           nlat=2
           inlat=5
           ltr(1,2)=0.5
           ltr(2,2)=0.5
           ltr(3,2)=0.5
         else if(lat == "R" .OR. lat == "r") then
           nlat=3
           inlat=6
           ltr(1,2)=2.0/3.0
           ltr(2,2)=1.0/3.0
           ltr(3,2)=1.0/3.0
           ltr(1,3)=1.0/3.0
           ltr(2,3)=2.0/3.0
           ltr(3,3)=2.0/3.0
         else if(lat == "F" .OR. lat == "f") then
           nlat=4
           inlat=7
           ltr(1,2)=0.5
           ltr(2,2)=0.5
           ltr(3,2)=0.0
           ltr(1,3)=0.5
           ltr(2,3)=0.0
           ltr(3,3)=0.5
           ltr(1,4)=0.0
           ltr(2,4)=0.5
           ltr(3,4)=0.5
         end if
         return
      end subroutine latsym

      subroutine build_orth()
!---------------------------------------------------------
! This subroutine calculates the reciprocal cell and the matrix for
! conversions to the orthonormal reference system.
!---------------------------------------------------------
         real :: sa,sb, sc
         if(ca-1.0 > 0.0) ca=cos(ca/rad)
         if(cb-1.0 > 0.0) cb=cos(cb/rad)
         if(cc-1.0 > 0.0) cc=cos(cc/rad)
         sa=sqrt(1.0-ca*ca)
         sb=sqrt(1.0-cb*cb)
         sc=sqrt(1.0-cc*cc)
         castar=(cb*cc-ca)/(sb*sc)
         cbstar=(ca*cc-cb)/(sa*sc)
         ccstar=(ca*cb-cc)/(sa*sb)
         sastar=sqrt (1.0-castar*castar)
         sbstar=sqrt (1.0-cbstar*cbstar)
         scstar=sqrt (1.0-ccstar*ccstar)
         astar=1.0/(a*sb*scstar)
         bstar=1.0/(b*sa*scstar)
         cstar=1.0/(c*sa*sbstar)
         asq=a*a
         bsq=b*b
         csq=c*c
         ab=a*b*cc
         ac=a*c*cb
         bc=b*c*ca

!  Building the crystallographic to orthonormal system matrix
!  The orthonormal system is defined having the axes as:
!  X axis coincides with the crystallographic A
!  Y axis belongs to the plane A,B
!  Z axis is perpendicular to the plane A,B

         ortd(1,1)=a
         ortd(1,2)=0.0
         ortd(1,3)=0.0
         ortd(2,1)=b*cc
         ortd(2,2)=b*sc
         ortd(2,3)=0.0
         ortd(3,1)=c*cb
         ortd(3,2)=-c*sb*castar
         ortd(3,3)=1/cstar


         orti(1,1)=1.0/ortd(1,1)
         orti(1,2)=0.0
         orti(1,3)=0.0
         orti(2,1)=-ortd(2,1)/(ortd(1,1)*ortd(2,2))
         orti(2,2)=1.0/ortd(2,2)
         orti(2,3)=0.0
         orti(3,1)=ortd(2,1)*ortd(3,2)-ortd(3,1)*ortd(2,2)
         orti(3,1)=cstar*orti(3,1)/(ortd(1,1)*ortd(2,2))
         orti(3,2)=-cstar*ortd(3,2)/ortd(2,2)
         orti(3,3)=cstar
         return
      end subroutine build_orth

      subroutine coordinates(t,n)
      !Constructs the arrays xcry and xcar containing the coordinates of the
      !atoms in cell numbered "n" corresponding to the lattice vector "t"
         real, dimension(3), intent(in) :: t
         integer,            intent(in) :: n
         real, dimension(3) ::  xc,xo
         integer :: i,j
         xc=0.0
         xo=0.0
         do i=1,nat
           xc(:)=x(1:3,i)+t(:)
           xcry(:,i,n)=xc(:)
           xo = matmul (xc,ortd) !Convert to cartesian
           xcar(:,i,n)=xo(:)
         end do
         return
      end subroutine coordinates

      subroutine spherical(xo,ss,theta,phi)
      !Calculates the spherical coordinates corresponding to the
      !Cartesian vector xo. It is valid for coordinates and moment values.
      !The angles theta and phi are output in degrees.
         real, intent(in)                         :: xo(3)
         real, intent(in out)                     :: ss
         real, intent(in out)                     :: theta
         real, intent(in out)                     :: phi
         ss=sqrt(dot_product(xo,xo))
         if(ss > 0.0) then
           theta=xo(3)/ss
           if(abs(theta) > 1.0) then
             theta=sign(1.0,theta)
           end if
           theta=acos(theta)
           if(abs(theta) < 1.0e-5 .or. abs(theta-pi) < 1.0e-5) then
             phi=0.0
           else
             phi=atan2(xo(2),xo(1))
           end if
         else
           theta=0.0
           phi=0.0
         end if
         theta=theta*rad
         phi=phi*rad
         return
      end subroutine spherical

      subroutine  momento(t,n,ia,cm,r)
      !This subroutine gives the value of the moment amplitude
      !of the atom "ia" in cell "n" corresponding to translation "t".
      !It provides also the cartesian components of the magnetic moment.
         real, dimension(3),intent(in)   :: t
         integer,           intent(in)   :: n
         integer,           intent(in)   :: ia
         real              ,intent(out)  :: r
         real, dimension(3),intent(out)  :: cm
         real, dimension(3) :: xc
         call magmom_calc(t,n,ia)
         xc(1)=rmo(1,ia,n)/a
         xc(2)=rmo(2,ia,n)/b
         xc(3)=rmo(3,ia,n)/c
         cm = matmul (xc,ortd) !Convert to cartesian
         r=sqrt(dot_product(cm,cm))
         return
      end subroutine momento



      subroutine magmom_calc(t,n,ia)
      !   Calculation of the magnetic moments of all the atoms located in a cell
      !   at a translation T from the zero-cell at (000). Crystallographic
      !   components stored in vector RMo(j,I,n)(j=1,2,3,I:atom index, n:cell index).
      !   If ia<=0 all atoms are calculated (from 1 to nat), otherwise only the moment
      !   of the atoms corresponding to the given index "ia" is calculated.
         real, dimension(3), intent(in)  :: t
         integer,            intent(in)  :: n
         integer,            intent(in)  :: ia

         integer :: i1,i2,i,j,iv
         real    :: arg,cost,sint

         if(ia <= 0) then    !Select the calculation for all atoms or a single atom
           i1=1
           i2=nat
         else
           i1=ia
           i2=ia
         end if

         do i=i1,i2                        ! Loop over atoms
           do j=1,3                        ! Loop over components x,y,z
             rmo(j,i,n)=0.0
             do iv=1,nvk                   ! Loop over propagation vectors
               arg=dot_product(t,pvk(:,iv)) ! Scalar product: k.R(l)
               if(ifou == 2) arg=arg+dot_product(x(1:3,i),pvk(:,iv)) ! Scalar products: k.r(i)
               cost=cos(tpi*( arg+ph(iv,i)+phk(iv) ) )
               sint=sin(tpi*( arg+ph(iv,i)+phk(iv) ) )
               rmo(j,i,n)=rmo(j,i,n)+rr(j,iv,i)*cost+ri(j,iv,i)*sint
             end do                        ! End loop over propagation vectors
           end do                          ! End loop over components x,y,z
         end do                            ! End loop over atoms
         return
      end subroutine magmom_calc

      subroutine phase_give(ival)
         integer, intent (in) :: ival
         integer              :: iphas, iv
         if(ival == 2) phasgiven=.FALSE.
         if(nvk >= 1 .and. .not. phasgiven) then
           iphas=nvk
           write(unit=*,fmt="(a,i3,a)")  &
               "  => Please enter ",iphas," phases of Fourier coefficients "
           write(unit=*,fmt="(a)") "     for the different k-vectors (in fractions of 2*Pi)"
           phk(1)=0.0
           do iv=1,iphas
             do
               write(unit=*,fmt="(a,i3,a)",advance="no") "  -> Enter phase",iv,": "
               read(unit=*,fmt=*,iostat=ier) phk(iv)
               if(ier == 0) exit
             end do
           end do
           phasgiven=.TRUE.
         end if
         return
      end subroutine phase_give

      subroutine phase_search()
         real,    dimension(3)      :: t,cm
         integer, dimension(3)      :: ic
         integer, dimension(natom)  :: inda
         real, dimension(nv,2)      :: pin
         real    :: diffo, rmino,rmaxo, rmax, rmin, dif, ran, rm
         integer :: n, iv, j, ik, i, iato, natin, np, i1,i2,i3, jl, nstep

         diffo= 99999999.0
         rmino= diffo
         call random_number(ran)
         rmaxo=-ran/diffo

         ic=1

         do iv=1,nvk                       !Search the maximum indices
           do j=1,3                              !for an approximant cell
             if(pvk(j,iv) /= 0) then
               ik=nint(abs(1.0/pvk(j,iv)))+1
               if(ik > ic(j)) ic(j)=ik
             end if
           end do
         end do
         n=1
         write(unit=iou,fmt="(/,70a1)") ("*",i=1,70)
         write(unit=iou,fmt="(10a1,a)") ("-",i=1,10),"   SEARCH FOR PHASES PHI(k)"
         write(unit=iou,fmt="(70a1,/)") ("*",i=1,70)
         write(unit=iou,fmt="(a,3f5.1,a)")  &
             " => Domain for cells search: from (0.0 0.0 0.0) to (", (real(ic(j)),j=1,3),")"
         do
           write(unit=*,fmt="(a)") " =>  Give the number of the testing atom "
           write(unit=*,fmt="(a)",advance="no") "    (If =0 you will be asked for a list): "
           read(unit=*,fmt=*,iostat=ier) iato
           if( ier == 0) exit
         end do
         if(iato == 0) then
           do
             write(unit=*,fmt="(a)",advance="no") " => How many atoms should be included?: "
             read(unit=*,fmt=*,iostat=ier) natin
             if( ier == 0) exit
           end do
           do i=1,natin
             do
               write(unit=*,fmt="(a,i2,a)",advance="no") " -> Enter number of the atom (",i,"): "
               read(unit=*,fmt=*,iostat=ier) inda(i)
               if( ier == 0) exit
             end do
           end do
         else
           natin=1
           inda(1)=iato
         end if
         write(unit=iou,fmt="(a,i3/a,18i3/26x,18i3/26x,18i3)")  &
             " => Number of atoms to be tested: ",natin,  &
             " => Indices of the atoms: ",(inda(i),i=1,natin)
         do
           write(unit=*,fmt="(a)",advance="no") " => Give the number of random trials: "
           read(unit=*,fmt=*,iostat=ier) nstep
           if( ier == 0) exit
         end do
         write(unit=iou,fmt="(/,a,i6,/)") " => Number of random trials: ",nstep
         write(unit=*,fmt="(a,/)") " => Give the intervals for phase search:"
         do iv=1,nvk
           do
             write(unit=*,fmt="(a,i2,a)",advance="no")  &
                 " -> Phase interval (2 reals modulo 2pi) for vector(",iv,"): "
             read(unit=*,fmt=*,iostat=ier) pin(iv,1),pin(iv,2)
             if( ier /= 0) cycle
             if(pin(iv,1) > pin(iv,2)) then
               write(unit=*,fmt="(a)") " => WARNING!: ",pin(iv,1),&
                              " must be lower than ",pin(iv,2),""
               cycle
             end if
             exit
           end do
         end do

         do np=0,nstep
           do iv=1,nvk
             call random_number(ran)
             phk(iv)=pin(iv,1)+(pin(iv,2)-pin(iv,1))*ran
           end do
           rmax=0.0
           rmin=99999.0
           do i1=0,ic(1)
             do i2=0,ic(2)
               do i3=0,ic(3)
                 do jl=1,nlat
                   t(:)=real((/i1,i2,i3/))+ ltr(:,jl)
                   do i=1,natin
                     iato=inda(i)
                      call momento(t,n,iato,cm,rm)
                     if(rm > rmax) then
                       rmax=rm
                     end if
                     if(rm < rmin) then
                       rmin=rm
                     end if
                   end do
                 end do
               end do
             end do
           end do
           dif=rmax-rmin
           if(dif < diffo) then
             diffo=dif
             rmino=rmin
             rmaxo=rmax
             write(unit=*,fmt="(a,2f8.3,a,5f8.4)")  &
                 " =>Moments: ", rmino,rmaxo," =>Phases: ",(phk(iv),iv=1,nvk)
             write(unit=iou,fmt="(a,2f8.3,a,5f8.4)")  &
                 " =>Moments: ", rmino,rmaxo," =>Phases: ",(phk(iv),iv=1,nvk)
             do iv=1,nvk
               phkop(iv)=phk(iv)
             end do
           end if
         end do
         write(unit=iou,fmt="(/a,2f8.3/a/)")  &
             " => Optimal moment fluctuation range: ", rmino,rmaxo,  &
             "    For the following phases: "
         do iv=1,nvk
           phk(iv)=phkop(iv)
           write(unit=iou,fmt="(a,i1,a,f8.4,a,f7.2,a)")  &
               " -> Optimal additional phase for Sk",iv,": ",  &
               phkop(iv)," (modulo 2pi) ->", phkop(iv)*rad*tpi," degrees"
           write(unit=iou,fmt="(a,2f8.4,a)")  &
               "    (Search within the range: ",(pin(iv,j),j=1,2)," ) "
           write(unit=*,fmt=  "(a,i1,a,f8.4,a,f10.2,a)")  &
               " -> Optimal additional phase for Sk",iv,": ",  &
               phkop(iv)," (modulo 2pi) ->", phkop(iv)*rad*tpi," degrees"
           write(unit=*,fmt="(a,2f8.4,a)")  &
               "    (Search within the range: ",(pin(iv,j),j=1,2)," ) "
         end do
         return
      end subroutine phase_search

   end module moment_mod


!---------------------------------------------------------------------------------


   module moment_inout

      use moment_mod
      use CFML_String_Utilities, only: FindFmt,  Init_FindFmt, ierr_fmt, mess_findfmt

      implicit none

     ! private
     ! public :: phase_bidon,


   contains



      subroutine input_data()
         character (len=1)  :: ans
         character (len=16) :: filout
         character (len=16) :: filin
         logical :: leer
         integer :: j

         leer=.FALSE.
         write(unit=*,fmt="(a)",advance="no")" => Interactive input (i=def) or read file (r): "
         read(unit=*,fmt="(a)") ans
         if(ans == "r" .OR. ans == "R") then
           leer=.TRUE.
           do
             write(unit=*,fmt="(a)",advance="no") " => Name of the  input file: "
             read(unit=*,fmt="(a)") filin
             open(unit=inp,file=filin,status="old",action="read",iostat=ier)
             if(ier /= 0) then
               write(unit=*,fmt="(a)") " => Problem reading the file: ",trim(filin)
               stop
             end if
             exit
           end do
         end if
         write(unit=*,fmt="(a)",advance="no") " => Name of the output file: "
         read(unit=*,fmt="(a)") filout
         open(unit=iou, file=filout,status="replace",action="write")
         write(unit=iou,fmt="(/,/,a,79a1)") " ",("-",j=1,79)
         write(unit=iou,fmt="(a)") "                         ----------------------"
         write(unit=iou,fmt="(a)") "                         --- Program MOMENT ---"
         write(unit=iou,fmt="(a)") "                         ----------------------"
         write(unit=iou,fmt="(a)") "                            (JRC-LLB,Sept-94)"
         write(unit=iou,fmt="(a,79a1,/,/)") " ",("-",j=1,79)
         if(leer) then
           read(unit=inp,fmt="(tr4,a)") title
           write(unit=iou,fmt="(a,a,/)")" => Data from input file: ",filin
           write(unit=iou,fmt="(a,a)") "    Title:",title
           write(unit=iou,fmt="(a,/)")  "    ------"
           call readfile()
           close(unit=inp)
         else
           call interactive_input()
         end if
         call latsym(symb)
         write(unit=iou,fmt="(/,a,a,/)")  " => Lattice type:",latt(inlat)
         return
      end subroutine input_data

      subroutine interactive_input()
         integer :: i,j,iv
         nfst(1)="S(k)"
         nfst(2)="T(k)"
         write(unit=*,fmt="(a)",advance="no")" => Give a title for the job: "
         read(unit=*,fmt="(a)") title
         write(unit=iou,fmt="(a)")" => Data given interactively"
         write(unit=iou,fmt="(a,a)") "    Title:",title
         write(unit=iou,fmt="(a/)")  "    ------"
!-------------------------------------------
!        Enter the unit cell
!-------------------------------------------
         do
           write(unit=*,fmt="(a)",advance="no")" => Unit cell parameters: "
           read(unit=*,fmt=*,iostat=ier)a,b,c,ca,cb,cc
           if(ier == 0) exit
         end do
         write(unit=iou,fmt="(a,6f10.4)")" => Unit cell :",a,b,c,ca,cb,cc
         call build_orth()
!-------------------------------------------
!        Enter the space group symbol
!-------------------------------------------
         do
           write(unit=*,fmt="(a)",advance="no") " => Space group: "
           read(unit=*,fmt="(a)",iostat=ier) symb
           if(ier == 0) exit
         end do
         write(unit=iou,fmt="(a,a)")" => Space group: ",symb
!-------------------------------------------
         do
           write(unit=*,fmt="(a)",advance="no")  &
             " => Type of Fourier coefficient S(k) -> 1 or T(k) -> 2: "
           read(unit=*,fmt=*,iostat=ier) ifou
           if(ier == 0) exit
         end do
         write(unit=iou,fmt="(a,a)")" => Type of Fourier Coefficients: ",nfst(ifou)
         write(unit=iou,fmt="(a)") " => For atom j of the reference cell (000) is given by: "
         if(ifou == 1) then
           write(unit=iou,fmt="(a)") "   S(k,j)= 0.5 {R(k,j) + i I(k,j)} exp{-2pi Phase(k,j)}"
         else
           write(unit=iou,fmt="(a)")  &
               "   T(k,j)= 0.5 {R(k,j) + i I(k,j)} exp{-2pi(k.r(j)+Phase(k,j))}"
         end if

         write(unit=iou,fmt="(a)") " => The magnetic moments are calculated using the formula: "
         if(ifou == 1) then
           write(unit=iou,fmt="(a/a/a)")  &
               "     M(L,j) = Sum{k}[ R(k,j) cos{2pi*(k.R(L)+Phase(k,j)}+ ",  &
               "                     +I(k,j) sin{2pi*(k.R(L)+Phase(k,j)}] ",  &
               "   Where R(L) is the lattice vector corresponding to cell L"
         else
           write(unit=iou,fmt="(a/a/a)")  &
               "     M(L,j) = Sum{k}[ R(k,j) cos{2pi*(k.R(L,j)+Phase(k,j)}+ ",  &
               "                     +I(k,j) sin{2pi*(k.R(L,j)+Phase(k,j)}] ",  &
               "   Where R(L,j) is the vector position of atom j in cell L"
         end if
!------------------------------------------------
!        Enter the number of propagation vectors
!------------------------------------------------
         do
           write(unit=*,fmt="(a)")   " => Number of propagation vectors"
           write(unit=*,fmt="(a)",advance="no") "      (k & -k are counted as one): "
           read(unit=*,fmt=*,iostat=ier)  nvk
           if(ier == 0) exit
         end do
         if(nvk < 0) nvk=abs(nvk)
         write(unit=iou,fmt="(a,i5)")" => Number of propagation vectors: ",nvk
!-------------------------------------------
!        Enter the propagation vectors
!-------------------------------------------
         do i=1,nvk
           do
             write(unit=*,fmt="(a,i1,a)",advance="no") " => Propagation vector (",i,"): "
             read(unit=*,fmt=*,iostat=ier) pvk(:,i)
             if(ier == 0) exit
           end do
           write(unit=iou,fmt="(a,i2,a,3f9.4)")"    Vector (",i,"): ", pvk(:,i)
         end do
!--------------------------------------------
!        Read the atom positions
!-------------------------------------------
         do
           write(unit=*,fmt="(a)",advance="no") " => Number of atoms: "
           read(unit=*,fmt=*,iostat=ier) nat
           if(ier == 0) exit
         end do
!--------------------------------------------------------------------------
         do i=1,nat
!--------------------------------------------------------------------------
           write(unit=*,fmt="(a,i1,a)",advance="no") " -> Name of atom (",i,"): "
           read(unit=*,fmt="(a)") nam(i)
           do
             write(unit=*,fmt="(a,i1,a)",advance="no") " -> Coordinates  (",i,")  (x,y,z): "
             read(unit=*,fmt=*,iostat=ier) x(1:3,i)
             if(ier == 0) exit
           end do
           write(unit=*,fmt="(a)")  &
               " -> Real & Imaginary parts and Phase of Fourier coefficients:"
           do iv=1,nvk
             do
                write(unit=*,fmt="(a,i1,a,i1,a)",advance="no")  &
                 " -> Rx,Ry,Rz  Ix,Iy,Iz  Phase(",i,")  Vector(",iv,"): "
                read(unit=*,fmt=*,iostat=ier) rr(:,iv,i),ri(:,iv,i),ph(iv,i)
                if(ier == 0) exit
             end do
           end do
!------------------------------------------------------------------
         end do
!------------------------------------------------------------------
         write(unit=iou,fmt="(/,a)")  &
             " => Atoms and Fourier coefficients of magnetic moments"
         write(unit=iou,fmt="(/,a)") "     Atom       X         Y         Z    "
         do i=1,nat
           write(unit=iou,fmt="(a,a4,a,3f10.5)") " -> ",nam(i),": ", x(1:3,i)
         end do
         write(unit=iou,fmt="(/,a)")  &
             "     Atom     Rx      Ry      Rz      Ix      Iy      Iz   Phase"
         do i=1,nat
           write(unit=iou,fmt="(a,a4,a,6f8.3,f8.5)")  &
               "    ",nam(i),": ",rr(:,1,i),ri(:,1,i),ph(1,i)
           do iv=2,nvk
             write(unit=iou,fmt="(a,6f8.3,f8.5)")  &
                 "          ",rr(:,iv,i),ri(:,iv,i),ph(iv,i)
           end do
         end do
         return
      end subroutine interactive_input

      subroutine cells_enter()
         real, dimension(3) :: t
         character(len=1)   :: ans
         integer :: j

         do
            nce=nce+1
            do
               write(unit=*,fmt="(a)",advance="no") " => Enter the origin of the Cell (m,n,p): "
               read(unit=*,fmt=*,iostat=ier) (t(j),j=1,3)
               if(ier == 0) exit
            end do
            call magmom_calc(t,nce,0)
            call write_mom(t)
            write(unit=*,fmt="(/,a)",advance="no") " -> More cells (y/n)?: "
            read(unit=*,fmt="(a)") ans
            if(ans == "y" .or. ans == "Y") cycle
            exit
         end do
         return
      end subroutine cells_enter

      subroutine angle_calc()
         real, dimension(3) :: t1, t2, cm1, cm2
         character(len=1)   :: ans
         integer :: iato1, iato2
         real    :: rm1, rm2, s
         do
           do
             write(unit=*,fmt="(a)",advance="no") " => Give two atom numbers: "
             read(unit=*,fmt=*,iostat=ier) iato1,iato2
             if(ier == 0) exit
           end do
           do
             write(unit=*,fmt="(a)",advance="no") " => Cell of the  first atom (m1,n1,p1): "
             read(unit=*,fmt=*,iostat=ier) t1
             if(ier == 0) exit
           end do
           do
             write(unit=*,fmt="(a)",advance="no") " => Cell of the second atom (m2,n2,p2): "
             read(unit=*,fmt=*,iostat=ier) t2
             if(ier == 0) exit
           end do
           call momento(t1,1,iato1,cm1,rm1)
           call momento(t2,2,iato2,cm2,rm2)
           write(unit=*,fmt="(a)")
           write(unit=*,fmt="(a,f8.3,a,3f8.3,a)")  " => Moment(Atom 1): ",rm1, " [",cm1," ]"
           write(unit=*,fmt="(a,f8.3,a,3f8.3,a)")  " => Moment(Atom 2): ",rm2, " [",cm2," ]"
           if(rm1 < 1.0e-08 .or. rm2 < 1.0e-08) then
             write(unit=*,fmt="(a)") " => One of the atoms has a nul moment!"
           end if
           s=dot_product(cm1,cm2)
           s=s/rm1/rm2
           if(abs(s) > 1.0) s=sign(1.0,s)
           s=acos(s)*rad
           write(unit=*,fmt="(a,f8.3,a)") " => Angle: ",s, " degrees"
           write(unit=*,fmt="(/,a)",advance="no") " -> More angles (y/n)?: "
           read(unit=*,fmt="(a)") ans
           if(ans == "y" .or. ans == "Y") cycle
           exit
         end do
         return
      end subroutine angle_calc

      subroutine phase_bidon()
         character(len=1)   :: ans
         real, dimension(3) ::  t,cm
         integer :: i,j,iv,n,np, iato, nstep
         real    :: rm1,rm2, ro, epsil, rm, delt
         n=1
         t=0.0
         do
           do
             write(unit=*,fmt="(a)",advance="no") " => Give the number of the testing atom: "
             read(unit=*,fmt=*,iostat=ier) iato
             if(ier == 0) exit
           end do
           do
             write(unit=*,fmt="(a)",advance="no") " => Give Limits for Moment (m1<m2): "
             read(unit=*,fmt=*,iostat=ier) rm1,rm2
             if(ier == 0 .AND. rm1 <= rm2) exit
           end do
           ro=0.5*(rm1+rm2)
           do
             write(unit=*,fmt="(a)",advance="no") " => Give the number of steps: "
             read(unit=*,fmt=*,iostat=ier) nstep
             if(ier == 0) exit
           end do
           do iv=2,nvk
             phko(iv,iato)=0.0
             epsil=abs(rm2-rm1)
             do np=0,nstep
               phk(iv)=0.5*real(np)/real(nstep)
               call momento(t,n,iato,cm,rm)
               if(rm >= rm1 .AND. rm <= rm2) then
                 delt=abs(ro-rm)
                 if(delt < epsil) then
                   epsil=delt
                   phko(iv,iato)=phk(iv)
                 end if
                 write(unit=*,fmt="(a,3f10.4)")" Phases&Mom: ",(phk(j),j=2,nvk),rm
               end if
             end do
           end do
           write(unit=*,fmt="(/,a)",advance="no") " -> More atoms (y/n)?: "
           read(unit=*,fmt="(a)") ans
           if(ans == "y" .OR. ans == "Y") cycle
           exit
         end do
         write(unit=*,fmt="(/,a,/)")"  => Optimun phases for all the atoms:"
         do i=1,nat
           write(unit=*,fmt="(a,a,4f10.5)") " -> Atom: ",nam(i),(phko(iv,i),iv=2,nvk)
           write(unit=*,fmt="(a,4f10.5)") "              ",(phko(iv,i)*rad*tpi,iv=2,nvk)
         end do
         return
      end subroutine phase_bidon

!-------------------------------------------------------------------------

      subroutine readfile()
         character (len=80)  :: fmtfields
         character (len=132) :: fmtformat
         character (len=132) :: aline
         character (len=9)   :: key,lab
         integer :: i,j,iv,ive
         nfst(1)="S(k)"
         nfst(2)="T(k)"
         !-----------------------------
         !        Read the unit cell
         !-----------------------------

         call init_findfmt(1)
         do
           fmtfields = "9ffffff"
           call findfmt(inp,aline,fmtfields,fmtformat)
           if(ierr_fmt == 13)  cycle
           exit
         end do
         if (ierr_fmt /= 0)  call finalize()
         read(unit=aline,fmt=fmtformat)  key,a,b,c,ca,cb,cc
         if(key(1:4) /= "CELL") then
           write(unit=*,fmt="(a)") " => Warning!, the CELL card could be in error."
         end if
         write(unit=*  ,fmt="(a,6f10.4)")" => Unit cell :",a,b,c,ca,cb,cc
         write(unit=iou,fmt="(a,6f10.4)")" => Unit cell :",a,b,c,ca,cb,cc
         call build_orth()
!-------------------------------------
!        Read the space group symbol
!-------------------------------------
         read(unit=inp,fmt="(a)")aline
         if(aline(1:4) /= "SPGR") then
           write(unit=*,fmt="(a)") " => Warning!, the SPGR card could be in error: "//trim(aline)
         end if
         aline =adjustl(aline(6:))
         if(len_trim(aline) == 0)  aline="P-1"
         symb=trim(aline)
         write(unit=*  ,fmt="(a,a)")" => Space group: ",symb
         write(unit=iou,fmt="(a,a)")" => Space group: ",symb
!--------------------------------------------------------
!        Read the type of fourier coefficients in input
!--------------------------------------------------------
         do
           fmtfields = "9i"
           call findfmt(inp,aline,fmtfields,fmtformat)
           if(ierr_fmt == 13)  cycle
           exit
         end do
         if (ierr_fmt /= 0)  call finalize()
         read(unit=aline,fmt=fmtformat) key,ifou
         if(abs(ifou) == 1) ifou=1
         if(abs(ifou) == 5) ifou=2
         if(key(1:4) /= "JBTV") then
           write(unit=*,fmt="(a)")  " => Warning!, the JBTV card could be in error: "//trim(aline)
         end if
         write(unit=*  ,fmt="(a,a)")" => Type of Fourier Coefficients: ",nfst(ifou)
         write(unit=iou,fmt="(a,a)")" => Type of Fourier Coefficients: ",nfst(ifou)
         write(unit=iou,fmt="(a)") " => For atom j of the reference cell (000) is given by: "
         if(ifou == 1) then
           write(unit=iou,fmt="(a)") "   S(k,j)= 0.5 {R(k,j) + i I(k,j)} exp{-2pi Phase(k,j)}"
         else
           write(unit=iou,fmt="(a)")  &
               "   T(k,j)= 0.5 {R(k,j) + i I(k,j)} exp{-2pi(k.r(j)+Phase(k,j))}"
         end if

         write(unit=iou,fmt="(a)") " => The magnetic moments are calculated using the formula: "
         if(ifou == 1) then
           write(unit=iou,fmt="(a/a/a)")  &
               "     M(L,j) = Sum{k}[ R(k,j) cos{2pi*(k.R(L)+Phase(k,j)}+ ",  &
               "                     +I(k,j) sin{2pi*(k.R(L)+Phase(k,j)}] ",  &
               "     Where R(L) is the lattice vector corresponding to cell L"
         else
           write(unit=iou,fmt="(a/a/a)")  &
               "     M(L,j) = Sum{k}[ R(k,j) cos{2pi*(k.R(L,j)+Phase(k,j)}+ ",  &
               "                     +I(k,j) sin{2pi*(k.R(L,j)+Phase(k,j)}] ",  &
               "     Where R(L,j) is the vector position of atom j in cell L"
         end if

!-----------------------------------------------
!        Read the number of propagation vectors
!-----------------------------------------------
         do
           fmtfields = "9i"
           call findfmt(inp,aline,fmtfields,fmtformat)
           if(ierr_fmt == 13)  cycle
           exit
         end do
         if (ierr_fmt /= 0)   call finalize()
         read(unit=aline,fmt=fmtformat)  key,nvk
         write(unit=*  ,fmt="(a,i5)")" => Number of propagation vectors: ",nvk
         write(unit=iou,fmt="(a,i5)")" => Number of propagation vectors: ",nvk
!--------------------------------------
!        Read the propagation vectors
!--------------------------------------
         do i=1,nvk
           do
             fmtfields = "9fff"
             call findfmt(inp,aline,fmtfields,fmtformat)
             if(ierr_fmt == 13)  cycle
             exit
           end do
           if (ierr_fmt /= 0)   call finalize()
           read(unit=aline,fmt=fmtformat)  key,pvk(:,i)
           if(key(1:3) /= "Vk_") then
             write(unit=*,fmt="(a)") " => Warning!, the number of vectors k could be in error: "//trim(aline)
           end if
           write(unit=*  ,fmt="(a,i2,a,3f9.4)")"    Vector (",i,"): ", pvk(:,i)
           write(unit=iou,fmt="(a,i2,a,3f9.4)")"    Vector (",i,"): ", pvk(:,i)
         end do
!----------------------------------
!        Read the atom positions
!----------------------------------
         nat=0

         do i=1,natom
           do
             fmtfields = "94fffff4"
             call findfmt(inp,aline,fmtfields,fmtformat)
             if(ierr_fmt == 13)  cycle
             exit
           end do
           if(ierr_fmt == -1) exit
           if (ierr_fmt /= 0)  call finalize()
           nat=nat+1
           read(unit=aline,fmt=fmtformat) key,nam(nat),x(:,nat),lab
           if(key(1:4) /= "ATOM") then
             write(unit=*,fmt="(a)") " => Warning!, the ATOM card could be in error "//trim(aline)
           end if
           do iv=1,nvk
             do
               fmtfields = "4ifffffff"
               call findfmt(inp,aline,fmtfields,fmtformat)
               if(ierr_fmt == 13)  cycle
               exit
             end do
             if (ierr_fmt /= 0)   call finalize()
             read(unit=aline,fmt=fmtformat)key,ive, rr(:,iv,nat),ri(:,iv,nat), ph(iv,nat)
             if(key(1:4) /= "MK_P") then
               write(unit=*,fmt="(a)") " => Warning!, the MK_P card could be in error "//trim(aline)
             end if
             if(abs(ive) /= iv .and. ive /= 0)  write(unit=*,fmt="(a)")  &
                 " => Warning!, the order of Fourier coefficient could be in error "
           end do
         end do  !loop_atom

         write(unit=iou,fmt="(/a)") " => Atoms and Fourier coefficients of magnetic moments"
         write(unit=iou,fmt="(/a)") "     Atom       X         Y         Z         Biso    Occ"
         do i=1,nat
           write(unit=iou,fmt="(a,a4,a,5f10.5)") " -> ",nam(i),": ", x(:,i)
         end do

         write(unit=iou,fmt="(/a)") "     Atom     Rx      Ry      Rz      Ix      Iy      Iz   Phase"
         do i=1,nat
           write(unit=iou,fmt="(a,a4,a,6f8.3,f8.5)")  "    ",nam(i),": ",rr(:,1,i),ri(:,1,i),ph(1,i)
           do iv=2,nvk
             write(unit=iou,fmt="(a,6f8.3,f8.5)")  "          ",rr(:,iv,i),ri(:,iv,i),ph(iv,i)
           end do
         end do
         return
      end subroutine readfile

      subroutine finalize()
         integer :: i
         if(ierr_fmt == -1) then
            write(unit=*,fmt="(a)")  " => End of input file "
         else
            do i=1,mess_findfmt%nlines
             write(unit=*,fmt="(a)")  mess_findfmt%txt(i)
            end do
         end if
         stop
      end subroutine finalize



      subroutine write_mom(t)
         real, intent(in)   :: t(3)
         real, dimension(3) :: xc,xo
         integer :: i,j,iv
         real    :: phi,theta,ss

         write(unit=*,fmt="(/a,3(f5.1,a)/)")  &
             "   ATOMS and MAGNETIC MOMENTS in cell: ( ",t(1), "," , t(2), "," ,t(3), " )"
         write(unit=iou,fmt="(//80a1)") ("-",j=1,80)
         write(unit=iou,fmt="(a,3(f5.1,a)/)")  &
             "   ATOMS and MAGNETIC MOMENTS in cell: ( ",t(1), "," , t(2), "," ,t(3), " )"
         do iv=1,nvk
           write(unit=iou,fmt="(a,i1,a,f8.4,a,f7.2,a)")  " -> Additional phase for Sk",iv,": ",  &
               phk(iv)," (modulo 2pi) ->", phk(iv)*rad*tpi," degrees"
           write(unit=*,fmt=  "(a,i1,a,f8.4,a,f10.2,a)") " -> Additional phase for Sk",iv,": ",  &
               phk(iv)," (modulo 2pi) ->", phk(iv)*rad*tpi," degrees"
         end do

         write(unit=*,fmt="(/,a)")    &
             "     Atom       X         Y         Z        XO        YO       ZO"
         write(unit=iou,fmt="(/,a)")  &
             "     Atom       X         Y         Z        XO        YO       ZO"
         call coordinates(t,nce)
         do i=1,nat
           write(unit=*,fmt="(a,a4,a,6f10.5)")  &
               " -> ",nam(i),": ", xcry(:,i,nce), xcar(:,i,nce)
           write(unit=iou,fmt="(a,a4,a,6f10.5)")  &
               " -> ",nam(i),": ", xcry(:,i,nce), xcar(:,i,nce)
         end do

!        Calculation of the modules of the magnetic moments during writing

         write(unit=*,fmt="(/a)")  &
             " Atom        Mx      My       Mz       MxO     MyO       MzO     Modul     Phi    Theta"
         write(unit=iou,fmt="(/a)")  &
             " Atom        Mx      My       Mz       MxO     MyO       MzO     Modul     Phi    Theta"
         do i=1,nat
          xc(1)=rmo(1,i,nce)/a
           xc(2)=rmo(2,i,nce)/b
           xc(3)=rmo(3,i,nce)/c
           xo = matmul (xc,ortd) !Convert to cartesian
           call spherical(xo,ss,theta,phi)
           rmc(:,i,nce)=xo(:)
           rmod(i,nce)=ss
           write(unit=*,fmt="(a,a4,a,9f9.3)")  &
               " ",nam(i),": ",rmo(:,i,nce),xo(:),rmod(i,nce) ,phi,theta
           write(unit=iou,fmt="(a,a4,a,9f9.3)")  &
               " ",nam(i),": ",rmo(:,i,nce),xo(:),rmod(i,nce) ,phi,theta
         end do
         return
      end subroutine write_mom
!-------------------------------------------------------------------------

      subroutine show_input(ilo)
        integer, intent(in) :: ilo
         character (len=1)   :: ans
         integer :: i,j,iv

         write(unit=ilo,fmt="(/a/)") "   Information about LATTICE and PROPAGATION VECTORS"
         write(unit=*,fmt="(a/4x,3f9.4,3f10.4)")  &
             " => Cell a,       b,       c,  cos(alpha),cos(beta),cos(gamma):",  &
              a,b,c,ca,cb,cc

         write(unit=ilo,fmt="(a,a)")     " => Space group: ",symb
         write(unit=ilo,fmt="(a,i5)")" => Number of propagation vectors: ",nvk
         if(nvk < 0) write(ilo,"(a/)")  "    Vectors -k are taken into account automatically  "
         do i=1,nvk
           write(unit=ilo,fmt="(a,i2,a,3f9.4)")"    Vector (",i,"): ",pvk(:,i)
         end do
         if(ilo == 6) then
           write(unit=*,fmt="(a)") " => Press <cr> to continue ..."
           read(unit=*,fmt="(a)") ans
         end if
         write(unit=ilo,fmt="(/a/)") "   Information about ATOMS and FOURIER COEFFICIENTS "
         write(unit=ilo,fmt="(/a)") "     Atom       X         Y         Z         Biso    Occ"
         do i=1,nat
           write(unit=ilo,fmt="(a,a4,a,5f10.5)") " -> ",nam(i),": ", x(:,i)
         end do
         write(unit=ilo,fmt="(/a)")  &
             "     Atom     Rx      Ry      Rz      Ix      Iy      Iz   Phase"
         do i=1,nat
           write(unit=ilo,fmt="(a,a4,a,6f8.3,f8.5)")  "    ",nam(i),": ",rr(:,1,i),ri(:,1,i),ph(1,i)
           do iv=2,nvk
             write(unit=ilo,fmt="(a,6f8.3,f8.5)")  "          ",rr(:,iv,i),ri(:,iv,i),ph(iv,i)
           end do
         end do
         if(ilo == 6) then
           write(unit=*,fmt="(a)") " => Press <cr> to continue ..."
           read(unit=*,fmt="(a)") ans
         end if
         write(unit=ilo,fmt="(/a)") " => Additional phase for each vector K:"
         do iv=1,nvk
           write(unit=ilo,fmt="(a,i2,a,f8.4,a,f10.4,a)")  "     Vector(",iv,") -> Phase: ",phk(iv),"(modulo 2pi)",  &
                                                          phk(iv)*rad*tpi," degrees"
         end do
         return
      end subroutine show_input

      subroutine write_file()
         character (len=1) :: ans
         character (len=16) :: filmom
         logical :: cartes
         real,    dimension(3) :: t(3),cm(3)
         integer, dimension(3) :: ic1(3),ic2(3)
         integer :: i,j,i1,i2,i3,jl,ia,n

         cartes=.FALSE.
         write(unit=*,fmt="(a)",advance="no") " => Name of the output file: "
         read(unit=*,fmt="(a)") filmom
         write(unit=*,fmt="(a)",advance="no") " => Cartesian (C) or Crystallographic framework (X=def): "
         read(unit=*,fmt="(a)") ans
         if(ans == "c" .OR. ans == "C") cartes=.TRUE.
         open(unit=iom, file=filmom,status="unknown")
         call show_input(iom)
         write(unit=iom,fmt="(/a)")  &
             "Atom       Cell           x       y       z       Mx      My      Mz    Modul"
         do
             write(unit=*,fmt="(a,/,a,/,a)",advance="no")  &
             " => Give integer indices for the range of unit cells",  &
             "    (e.g.:  x1 x2    y1 y2    z1 z2 ",  &
             "            -2  2    -1  3    -1  2 ,for example...): "
             read(unit=*,fmt=*,iostat=ier) (ic1(i),ic2(i), i=1,3)
             if(ier == 0) exit
         end do
         n=0
         do i1=ic1(1),ic2(1)
           do i2=ic1(2),ic2(2)
             do i3=ic1(3),ic2(3)
               do jl=1,nlat
                 t(:)=real((/i1,i2,i3/))+ ltr(:,jl)
                 n =n +1
                 call coordinates(t,n)
                 do ia=1,nat
                   call momento(t,n,ia,cm,rmod(ia,n))
                   if(cartes) then
                     write(unit=iom,fmt="(a4,1x,3f5.1,1x,3f8.2,1x,3f8.3,1x,f6.2)")  &
                           nam(ia), t(:), xcar(:,ia,n), cm(:), rmod(ia,n)
                   else
                     write(unit=iom,fmt="(a4,1x,3f5.1,1x,3f8.4,1x,3f8.3,1x,f6.2)")  &
                           nam(ia),t(:), xcry(:,ia,n), rmo(:,ia,n), rmod(ia,n)
                   end if
                 end do
               end do
             end do
           end do
         end do
         write(unit=iom,fmt="(/,a,i5)")" => Total number of cells: ",n
         close(unit=iom)
         return
      end subroutine write_file
!----------------------------------------------------------------
      subroutine write_mag3d()
         character (len=2),dimension(10)    :: spec
         real,    dimension(3) :: t(3),cm(3)
         integer :: i,j,jl,ia,n,nspec,ns
         real    :: phi,theta,ss

         open(unit=iom,file="mag3d.cry",status="unknown")
         write(unit=iom,fmt="(a,a)")"N",title
         write(unit=iom,fmt="(a,3f8.4,3f8.2)") "C", a,b,c,acos(ca)*rad,acos(cb)*rad,acos(cc)*rad
         write(unit=iom,fmt="(a)")"S x,y,z"
         n=0
         nspec=1
         spec(1)=nam(1)(1:2)
      at:do  ia=2,nat
           do ns=1,nspec
             if(nam(ia)(1:2) == spec(ns)) cycle at
           end do
           nspec=nspec+1
           spec(nspec)=nam(ia)(1:2)
         end do at
!        Writing Atom cards
         do jl=1,nlat
           t(:)=ltr(:,jl)
           call coordinates(t,n)
           do ia=1,nat
             write(unit=iom,fmt="(a,a4,1x,3f8.4,a)") "A ",nam(ia),xcry(:,ia,n),"  0.000"
           end do
         end do
!          Writing F cards
         do ns=1,nspec
           write(unit=iom,fmt="(a,a,a)")"F ",spec(ns)," 1 0.9"
           write(unit=iom,fmt="(a,a,a)")"F ",spec(ns)//"M ","3 0. 1. 0.05 .9 .8 .01"
         end do
!          Writing Q-form cards
!             Q FeM FORM Fe1
     ato:do  ia=1,nat
           do ns=1,nspec
             if(nam(ia)(1:2) == spec(ns)) then
               write(unit=iom,fmt="(a,a,a)")"Q ",spec(ns)//"M FORM ",nam(ia)
               cycle ato
             end if
           end do
         end do ato
!          Writing Q-PROP card
!             Q PROP 0 0 0.5
         do i=1,nvk
           write(unit=iom,fmt="(a,3f9.4)")"Q PROP ", pvk(:,i)
         end do
!         Writing Q-MU and Q-SDIR cards
!            Q Fe1  MU 1.2        (in Bohr magnetons)
!            Q Fe1 SDIR  50 -30   (theta,phi in degrees)
         do ia=1,nat
           call momento(t,n,ia,cm,ss)
           call spherical(cm,ss,theta,phi)
           write(unit=iom,fmt="(a,a4,a,f5.2)")  "Q ",nam(ia)," MU  ", ss
           write(unit=iom,fmt="(a,a4,a,2f8.2)") "Q ",nam(ia)," SDIR  ", theta,phi
         end do
!                         Writing X cards, examples:
!        X PERS .5 0.766 0.6428
!        X ARRO .333 0.1667 .0833 .02 1.2
!        X TRAN  1  0  0  0  1  0  0  0  1
!        X SYMB Fe1 .25 0.25
!        X CM/A 1

         write(unit=iom,fmt="(a)") "X PERS .5 0.766 0.6428"
         write(unit=iom,fmt="(a)") "X ARRO .333 0.1667 .0833 .02 1.2"
         write(unit=iom,fmt="(a)") "X TRAN  1  0  0  0  1  0  0  0  1"
         write(unit=iom,fmt="(a)") "X CM/A 1"
         do ia=1,nat
           write(unit=iom,fmt="(a,a,a)")"X SYMB ",nam(ia), " 0.25  0.25"
         end do
         close(unit=iom)
         return
      end subroutine write_mag3d

   end module moment_inout

!-------------------------------------------------------------------------


      program moment
      use moment_mod
      use moment_inout

!     Initialization of some variables
      phasgiven=.FALSE.
      nce=0
      iin=0
!
      write(unit=*,fmt="(a)") "   ------------------------------------------"
      write(unit=*,fmt="(a)") "   -        --- Program MOMENT ---          -"
      write(unit=*,fmt="(a)") "   -           (JRC-LLB,Sept-94)            -"
      write(unit=*,fmt="(a)") "   ------------------------------------------"
      write(unit=*,fmt="(a)") " "
      do
       do
          write(unit=*,fmt="(/,a,/)") "      OPTIONS: (enter the option number)  "
          write(unit=*,fmt="(a)")     "      0: Read file or enter input data"
          write(unit=*,fmt="(a)")     "      1: Show the input data"
          write(unit=*,fmt="(a)")     "      2: Change Phase angle between Fourier Coeff."
          write(unit=*,fmt="(a)")     "      3: Calculate magnetic moments in several cells"
          write(unit=*,fmt="(a)")     "      4: Write a file for MAG3D"
          write(unit=*,fmt="(a)")     "      5: Angles between selected magnetic moments"
          write(unit=*,fmt="(a)")     "      6: Search phase between Fourier Coeff.(general)"
          write(unit=*,fmt="(a)")     "      7: Write a file with magnetic moments for plot"
          write(unit=*,fmt="(a)")     "     20: Stop calculations"
          write(unit=*,fmt="(a)")     "  "
          write(unit=*,fmt="(a)",advance="no") " =>  Option: "
          read(unit=*,fmt=*,iostat=ier) ival
          if (ier == 0) exit
       end do

      if(iin == 0 .and. ival < 20 .and. ival /= 0) then
        write(unit=*,fmt="(a)") " => PLEASE, select first option 0 !"
        cycle
       end if

       select case (ival)
         case(  0)
            call input_data()
            iin=1
         case(  1)
            call show_input(6)
         case(  2)
            call phase_give(ival)
         case(  3)
            call phase_give(ival)
            call cells_enter()
         case(  4)
!           Call Phase_bidon()
            call write_mag3d()
         case(  5)
            call angle_calc()
         case(  6)
            call phase_search()
         case(  7)
            call write_file()
         case( 20)
            stop
       end select
      end do

      end program moment
!-------------------------------------------------------------------------

