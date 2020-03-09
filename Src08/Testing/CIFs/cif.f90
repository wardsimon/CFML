!!----
!!----
!!----
!!----
 Program Test_CIF_CFL
    !---- Use Modules ----!
    use CFML_Globaldeps
    !use CFML_Symmetry_Tables
    use CFML_Metrics
    use CFML_Maths
    use CFML_gSpaceGroups
    use CFML_Rational
    use CFML_IOForm
    use CFML_Atoms
    use CFML_Strings,only: File_type,Set_Symb_From_Mat !pack_string
    implicit none

    character(len=256)                  :: fname
    character(len=256)                  :: setting
    character(len=25)                   :: forma
    character(len=5)                    :: aux
    type(Cell_G_Type)                   :: Cell,Celln
    !type(Spg_Type)                      :: Grp
    type(AtList_Type)                   :: Atm
    type(File_type)                     :: flist
    type(SuperSpaceGroup_Type)          :: Grp
    character (len=*), dimension(26),parameter   :: &
    cdd=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r', &
         's','t','u','v','w','x','y','z']
    !type(rational),   dimension(:,:),allocatable :: Mat
    !character(len=40),dimension(:,:),allocatable :: matrix
    integer :: i, j, L,k, d,Dd,nsg, ind, indexg, num_group, ier,mult,codini
    real(kind=cp) :: start, fin
    real(kind=cp), dimension(3,256) :: orb,morb
    real(kind=cp), dimension(3)     :: codes=1.0
    integer, dimension(256) :: ptr

    call Set_Eps_Math(0.0002_cp)
    do

       write(*,'(/,a)',advance='no') " => Introduce the name of the CFL file: "
       read(*,"(a)") fname
       if(len_trim(fname) == 0) exit
       call CPU_TIME(start)
       call Readn_Set_Xtal_Structure(fname,Cell,Grp,Atm,"MAtm_std","CFL")!,file_list=flist) !,Iphase,Job_Info,file_list,CFrame)
       if(Err_CFML%Ierr == 0) then
          !write(*,"(/,a,/)")  " => Content of the CFL-file: "//flist%Fname
          !do i=1,flist%nlines
          !   write(*,"(i6,a)") i,"    "//flist%line(i)%Str
          !end do
          call Write_Crystal_Cell(Cell)
          if(len_trim(Grp%setting) /= 0) then
            write(*,"(/,a)") " => Transformed Cell"
            if(Grp%D > 4) then
              i=index(Grp%setting,"d")
              setting=Grp%setting(1:d-2)//";0,0,0"
            else
              setting=Grp%setting
            end if
            call Change_Setting_Cell(Cell,setting,Celln)
            call Write_Crystal_Cell(Celln)
          end if
          call Write_SpaceGroup_Info(Grp)
          !if(Atm%natoms > 0) call Write_Info_Atom_List(Atm)

          if(Atm%natoms > 0) then
             write(*,"(a,i5)") "  Number of atoms:",Atm%natoms
             call Write_Atom_List(Atm)
             !Calculate all atoms in the unit cell
             do i=1,Atm%natoms
               codini=1
               call Get_moment_ctr(Atm%Atom(i)%x,Atm%Atom(i)%moment,Grp,codini,codes) !,Ipr=6)
               !write(*,"(a,2(i3,3f10.5))") "Moment: ",i,Atm%Atom(i)%moment,codini,codes
               call Get_Orbit(Atm%Atom(i)%x,Atm%Atom(i)%moment,Grp,Mult,orb,morb,ptr)
               write(*,"(a)") " => Orbit of atom: "//Atm%Atom(i)%Lab
               do j=1,Mult
                   write(*,"(i5,3f10.5,tr8,3f10.5,i8)") j,orb(:,j),morb(:,j),ptr(j)
               end do
             end do
          end if
       else
          write(*,'(/,a)') " => ERROR: "//trim(Err_CFML%Msg)
       end if
       call CPU_TIME(fin)
       write(*,"(/,a,f12.3,a)") "CPU_TIME for this calculation: ",fin-start," seconds"
    end do

    contains
    !!--++
    !!--++  Subroutine Get_Orbit(x,mom,Spg,Mult,orb,morb,ptr)
    !!--++     !---- Arguments ----!
    !!--++     real(kind=cp), dimension(3),    intent (in) :: x
    !!--++     type(Magnetic_Space_Group_type),intent (in) :: spg
    !!--++     integer,                        intent(out) :: mult
    !!--++     real(kind=cp),dimension(:,:),   intent(out) :: orb
    !!--++     integer,dimension(:),optional,  intent(out) :: ptr
    !!--++
    !!--++     Calculates the orbit and mutiplicity of an atom position x
    !!--++     in a crystal structure described by a Shubnikov Group
    !!--++     The values of the moments along the orbit is also calculated
    !!--++
    !!--++     Updated: January 2020
    !!--++
    Subroutine Get_Orbit(x,mom,Spg,Mult,orb,morb,ptr)
       !---- Arguments ----!
       real(kind=cp), dimension(:),    intent (in) :: x,mom
       class(SpG_type),                intent (in) :: spg
       integer,                        intent(out) :: mult
       real(kind=cp),dimension(:,:),   intent(out) :: orb,morb
       integer,dimension(:),optional,  intent(out) :: ptr

       !---- Local variables ----!
       integer                                :: j, nt
       real(kind=cp), dimension(3)            :: xx,v,mmom,w
       real(kind=cp), dimension(3,3)          :: rmat
       character(len=1)                       :: laty

       laty="P"
       mult=1
       orb(:,1)=x(:)
       morb(:,1)=mom(:)
       if(present(ptr)) ptr(mult) = 1
       do_ext: do j=2,Spg%Multip
          rmat=Spg%Op(j)%Mat(1:3,1:3)
          xx=Apply_OP(Spg%Op(j),x)
          xx=modulo_lat(xx)
          do nt=1,mult
             v=orb(:,nt)-xx(:)
             if(Spg%Num_lat > 0) then
               if (is_Lattice_vec(v,Spg%Lat_tr)) cycle do_ext
             else
               if (Zbelong(v)) cycle do_ext
             end if
          end do
          mmom=Spg%Op(j)%dt*Spg%Op(j)%time_inv*matmul(rmat,mom)
          mult=mult+1
          orb(:,mult)=xx(:)
          morb(:,mult)=mmom(:)
          if(present(ptr)) ptr(mult) = j   !Pointer to symmetry operator
       end do do_ext
       return
    End Subroutine Get_Orbit

   !!
   !!----  Subroutine get_moment_ctr(xnr,TFourier,Spg,codini,codes,ord,ss,att,Ipr)
   !!----     real(kind=cp), dimension(3),            intent(in    ) :: xnr    !Atom position (fractional coordinates)
   !!----     Complex(kind=cp), dimension(3),         intent(in out) :: TFourier !Fourier coefficients at position xnr
   !!----     class(SuperSpaceGroup_Type),            intent(in)     :: Spg    !Super Space Group
   !!----     Integer,                                intent(in out) :: codini !Last attributed parameter
   !!----     real(kind=cp), dimension(3),            intent(in out) :: codes  !codewords for positions
   !!----     integer,                       optional,intent(in)     :: ord
   !!----     integer, dimension(:),         optional,intent(in)     :: ss
   !!----     real(kind=cp), dimension(:,:), optional,intent(in)     :: att
   !!----     integer,                       optional,intent(in)     :: Ipr
   !!----
   !!----  Subroutine to get the appropriate constraints in the refinement codes of
   !!----  magnetic moment parameters.
   !!----  Algorithm based in the Wigner theorem.
   !!----  The vector Mom = Sum { R Moment} displays the symmetry constraints to be
   !!----  applied to the magnetic moments. The sum runs over all magnetic
   !!----  matrices of the stabilizer of the particular atom position in the given
   !!----  space group.
   !!----
   !!----   Updated: March 2020
   !!----
   !!
   Subroutine Get_TFourier_ctr(xnr,kinf,TFourier,codes,SpG,codini,ord,ss,att,Ipr)
       real(kind=cp), dimension(3),            intent(in)     :: xnr
       type(kvect_info_type),                  intent(in)     :: kinf
       real(kind=cp), dimension(:),            intent(in out) :: TFourier
       real(kind=cp), dimension(:),            intent(in out) :: codes
       class(SpG_Type),                        intent(in)     :: SpG
       Integer,                                intent(in out) :: codini
       integer,                       optional,intent(in)     :: ord
       integer, dimension(:),         optional,intent(in)     :: ss
       real(kind=cp), dimension(:,:), optional,intent(in)     :: att
       integer,                       optional,intent(in)     :: Ipr

       ! Local variables
       integer,          dimension(kinf%nk)    :: brack_m  ! hs=(H,[m])
       real(kind=cp),    dimension(spg%d-1)    :: ts  ! ts=(t,tI)
       integer,          dimension(kinf%nk)    :: mE  ![m].E
       integer,          dimension(kinf%nk,3)  :: Mx  !M
       integer,          dimension(3,3)        :: g, magm   !g, magm= delta * det(g) * g
       integer,      dimension(kinf%nk,kinf%nk):: ep  !E
       character(len=1), dimension(26)         :: codd
       integer                                 :: i,j,order,iq,ir,d,n,k,iqt,ig,ip,iss
       real(kind=cp)                           :: suma,alpha
       integer,           dimension(48)        :: ss_ptr
       real(kind=cp),     dimension(3,48)      :: atr
       real(kind=cp),     dimension(6*kinf%nq) :: cod,multi,val
       real(kind=cp),     dimension(3)         :: x
       real(kind=dp),     dimension(6,6)       :: subm
       real(kind=dp),  dimension(6*kinf%nq,6*kinf%nq) :: Ctr,sCtr
       real(kind=cp),     parameter                   :: epss=0.001_cp
       real(kind=cp), dimension(6*kinf%nq)            ::TFourL,sTF


       !Test if all codes are given ... in such a case the user constraints
       !are prevalent

       d=SpG%d-1
       suma=0.0
       !iq=0
       n=0
       do k=1,kinf%nq
         do j=1,6
            n=n+1
            suma=suma+abs(codes(n))
         end do
       end do

       !if(suma < epss .or. iq > 0 ) return  !No refinement is required
       if(suma < epss ) return  !No refinement is required

       x=xnr
       where(x < 0.0) x=x+1.0
       where(x > 1.0) x=x-1.0

       if(present(ord) .and. present(ss) .and. present(att)) then
         order=ord
         ss_ptr(1:order) = ss(1:ord)
         atr(:,1:order)  = att(:,1:ord)
       else
         call get_stabilizer(x,SpG,order,ss_ptr,atr)
         if(present(ipr)) Write(unit=ipr,fmt="(a,i3)") " => Superspace stabilizer without identity, order:",order
       end if


       TFourL=TFourier
       sCtr=0.0_cp
       sTF=0.0
       if(order > 1) then
         n=6*kinf%nq
         !do ir=1,spg%multip
         do ig=1,order
           ir=ss_ptr(ig)
               g(:,:) = SpG%Op(ir)%Mat(1:3,1:3)   !                          /  g    0   t  \
                ts(:) = SpG%Op(ir)%Mat(1:d,d+1)   !   Superspace operator:  |  Mx   ep   tI |    !ts=(t,tI)
              Mx(:,:) = SpG%Op(ir)%Mat(4:d,1:3)   !                         \   0    0   1 /
              Ep(:,:) = SpG%Op(ir)%Mat(4:d,4:d)
            magm(:,:) = g(:,:)*SpG%Op(ir)%time_inv*SpG%Op(ir)%dt
            !Identify the q_coefficients of mM to select the proper Tfourier
            !Construction of the Cr-matrix
            Ctr=0.0
            do ip=1,kinf%nq
              brack_m(:)=kinf%q_coeff(:,ip) !nk-components vector
              mE=matmul(brack_m,Ep)   ![m].Ep  nk-components vector
              do iqt=1,kinf%nq
                iq=iqt
                if(equal_vector(mE,kinf%q_coeff(:,iqt))) then
                  iss=1
                  exit
                end if
                if(equal_vector(mE,-kinf%q_coeff(:,iqt))) then
                  iss=-1
                  exit
                end if
              end do !iq
              alpha=tpi*dot_product(brack_m,matmul(Mx,x)+ts(4:))
              subm(1:3,1:3)=cos(alpha)*magm
              subm(4:6,1:3)=sin(alpha)*magm
              subm(1:3,4:6)=-iss*sin(alpha)*magm
              subm(4:6,4:6)= iss*cos(alpha)*magm
              !Block  ip,iq
              i=(ip-1)*6; j=(iq-1)*6
              Ctr(i+1:i+6,j+1:j+6)=subm
              Tfourl(i+1:i+6)=matmul(subm,TFourier(j+1:j+6))
              sTF(i+1:i+6)=sTF(i+1:i+6)+TfourL(i+1:i+6)
            end do !ip
            sCtr=sCtr+Ctr !Adding constraint matrices for each operator of stabilizer
            if(present(ipr)) then
              do ip=1,kinf%nq
                 i=(ip-1)*6
                 write(unit=ipr,fmt='(a,i2,a,t20,a,t55,6f14.4,4i3)') '     Operator ',ir,": ",trim(Spg%Symb_Op(ir)),Tfourl(i+1:i+6),brack_m(:)
              end do
            end if
         end do  !ig operators
         sCtr=sCtr/order
         sTF=sTF/order
         if(present(ipr)) then
           do ip=1,kinf%nq
              i=(ip-1)*6
              write(unit=ipr,fmt='(a,6f14.4)')     '     Sum of TFour: ',sTf(i+1:i+6)
           end do
         end if
         call Get_Refinement_Codes(n,sTF,sCtr,iss,multi,codd,TFourL)
         cod=0.0
         do j=1,n
           if(codd(j) /= "0") then
             do i=1,iss
               if(codd(j) == cdd(i)) then
                 cod(j)=codini+i
                 exit
               end if
             end do
           end if
         end do
         TFourier=TFourL
         codes=0.0
         do j=1,n
           if(abs(multi(j)) > epss)  codes(j) = sign(1.0_cp, multi(j))*(abs(cod(j))*10.0_cp + abs(multi(j)) )
         end do
         codini=codini+iss
         if(present(Ipr)) then
           Write(Ipr,"(a,i4)")      " Number of free parameters: ",iss
           Write(Ipr,"(a,24F14.6)") " Basic Values: ",val(1:iss)
           write(Ipr,"(a,24f14.6)") " Multipliers: ",(multi(j), j=1,n)
           write(Ipr,"(28a)")       " String with free parameters: ( ",(codd(j)//", ",j=1,n-1),codd(n)//" )"
           write(Ipr,"(a,24i6)")    " Resulting integer codes: ", nint(cod(1:n))
           write(Ipr,"(a,24f14.6)") " Final codes: ",codes(1:n)
           write(Ipr,"(a,24f14.6)") " Constrained Fourier Coefficients: ",TFourier
         end if

       else !No restrictions
         codd(1:n)=cdd(1:n)
         multi(1:n)=1.0_cp
         do j=1,n
           cod(j)=codini+j
           codes(j) = (abs(cod(j))*10.0_cp + abs(multi(j)))
         end do
         codini=codini+n
         if(present(Ipr)) then
           write(Ipr,"(a,24f10.6)") " General position, no constraints in Fourier Coefficients "
           write(Ipr,"(28a)")       " String with free parameters: ( ",(codd(j)//", ",j=1,n-1),codd(n)//" )"
           write(Ipr,"(a,24i6)")    " Resulting integer codes: ", nint(cod(1:n))
           write(Ipr,"(a,24f14.6)") " Final codes: ",codes(1:n)
           write(Ipr,"(a,24f14.6)") " Constrained Fourier Coefficients: ",TFourier
         end if

       end if
       return
   End Subroutine Get_TFourier_ctr

    Subroutine Get_moment_ctr(xnr,moment,Spg,codini,codes,ord,ss,att,Ipr)
       real(kind=cp), dimension(3),            intent(in)     :: xnr
       real(kind=cp), dimension(:),            intent(in out) :: moment
       class(SpG_type),                         intent(in)   :: Spg
       Integer,                                intent(in out) :: codini
       real(kind=cp), dimension(:),            intent(in out) :: codes
       integer,                       optional,intent(in)     :: ord
       integer, dimension(:),         optional,intent(in)     :: ss
       real(kind=cp), dimension(:,:), optional,intent(in)     :: att
       integer,                       optional,intent(in)     :: Ipr

       ! Local variables
       character(len=1),  dimension(3)   :: codd
       character(len=:), allocatable     :: mag
       integer                           :: i,j,order,n,ig,iss
       real(kind=cp)                     :: suma
       integer,           dimension(48)  :: ss_ptr
       real(kind=cp),     dimension(3,48):: atr
       real(kind=cp),     dimension(3)   :: cod,multi
       real(kind=cp),     dimension(3)   :: x
       real(kind=cp),     dimension(3,3) :: magm  !g, magm= delta * det(g) * g
       real(kind=dp),     dimension(3,3) :: sCtr
       real(kind=cp),     dimension(3)   :: momentL,TotMom


       !Test if all codes are given ... in such a case the user constraints
       !are prevalent

       suma=0.0_cp
       !iq=0
       n=3 !Real moments -> three components
       do j=1,3
          suma=suma+abs(codes(j))
       end do

       if(suma < epss ) return  !No refinement is required
       if(present(Ipr)) then
         write(Ipr,"(/,a)")         " => Calculation of symmetry constraints for magnetic moments "
       end if
       x=xnr
       !where(x < 0.0) x=x+1.0
       !where(x > 1.0) x=x-1.0

       if(present(ord) .and. present(ss) .and. present(att)) then
         order=ord
         ss_ptr(1:order) = ss(1:ord)
         atr(:,1:order)  = att(:,1:ord)
       else
         call get_stabilizer(x,SpG,order,ss_ptr,atr)
         if(present(ipr)) Write(unit=ipr,fmt="(a,i3)") " => Stabilizer without identity, order:",order
       end if

       momentL=moment
       sCtr=0.0_cp
       if(order > 1) then
         do ig=1,order
           j=ss_ptr(ig)
           magm(:,:) = real(Spg%Op(j)%Mat(1:3,1:3))*Spg%Op(j)%dt*Spg%Op(j)%time_inv
           mag=Set_Symb_From_Mat(magm,["u","v","w"])
           sCtr=sCtr+magm !Adding constraint matrices for each operator of stabilizer
           if(present(ipr)) then
             write(unit=ipr,fmt='(a,i2,a,t20,a,t55,a,t75,9f8.4)') '     Operator ',ig,": ",trim(Spg%Symb_Op(j)), &
              trim(mag), sCtr
           end if
         end do  !ig operators
         sCtr=sCtr/order
         suma=sum(abs(sCtr))
         !write(*,"(a,f10.4,a,i3)") " suma:",suma, "Mag_Type:", spg%mag_type
         if(suma < epss .or. spg%mag_type == 2) then !This corresponds to a grey point group
            moment=0.0_cp
            codes=0.0_cp
            if(present(Ipr)) then
              write(Ipr,"(a)")         " Grey point group or symmetry paramagnetic site: the attached moment is zero "
              write(Ipr,"(a,24f14.6)") " Final codes: ",codes(1:n)
              write(Ipr,"(a,24f14.6)") " Constrained moment: ",moment
            end if
            return
         end if
         TotMom=matmul(sCtr,momentL)
         call Get_Refinement_Codes(n,TotMom,sCtr,iss,multi,codd,momentL)
         cod=0.0
         do j=1,n
           if(codd(j) /= "0") then
             do i=1,iss
               if(codd(j) == cdd(i)) then
                 cod(j)=codini+i
                 exit
               end if
             end do
           end if
         end do
         moment=momentL
         codes=0.0
         do j=1,n
           if(abs(multi(j)) > epss)  codes(j) = sign(1.0_cp, multi(j))*(abs(cod(j))*10.0_cp + abs(multi(j)) )
         end do
         codini=codini+iss
         if(present(Ipr)) then
           Write(unit=Ipr,fmt="(a,i4)")       " Number of free parameters: ",iss
           write(unit=Ipr,fmt="(a,3f14.6)")   " Multipliers: ",(multi(j), j=1,n)
           write(unit=Ipr,fmt="(28a)")        " String with free parameters: ( ",(codd(j)//", ",j=1,n-1),codd(n)//" )"
           write(unit=Ipr,fmt="(a,3i6)")      " Resulting integer codes: ", nint(cod(1:n))
           write(unit=Ipr,fmt="(a,3f14.6)")   " Final codes: ",codes(1:n)
           write(unit=Ipr,fmt="(a,3f14.6)")   " Constrained Moment: ",moment
         end if

       else !No restrictions

         codd(1:n)=cdd(1:n)
         multi(1:n)=1.0_cp
         do j=1,n
           cod(j)=codini+j
           codes(j) = (abs(cod(j))*10.0_cp + abs(multi(j)))
         end do
         codini=codini+n
         if(present(Ipr)) then
           write(unit=Ipr,fmt="(a)")         " General position, no constraints in moment "
           write(unit=Ipr,fmt="(28a)")       " String with free parameters: ( ",(codd(j)//", ",j=1,n-1),codd(n)//" )"
           write(unit=Ipr,fmt="(a,24i6)")    " Resulting integer codes: ", nint(cod(1:n))
           write(unit=Ipr,fmt="(a,24f14.6)") " Final codes: ",codes(1:n)
           write(unit=Ipr,fmt="(a,24f14.6)") " Constrained moment: ",moment
         end if

       end if
       return
    End Subroutine Get_moment_ctr

    Subroutine Get_Refinement_Codes(n,vect_val,Ctr,iss,multi,codd,vect_out)
      integer,                       intent(in)    :: n !dimension of the vector and the matrix
      real(kind=cp), dimension(:),   intent(in)    :: vect_val
      real(kind=dp), dimension(:,:), intent(in out):: Ctr
      integer,                       intent(out)   :: iss
      real(kind=cp), dimension(:),   intent(out)   :: multi
      character(len=*), dimension(:),intent(out)   :: codd
      real(kind=cp), dimension(:),   intent(out)   :: vect_out
      !--- Local variables ---!
      real(kind=cp), dimension(n)   :: val
      integer,       dimension(n)   :: pti
      real(kind=dp), dimension(n,n) :: zv
      integer                       :: i,j,k,kval,ier,ip
      real(kind=dp)                 :: zmi
      real(kind=dp), dimension(n)   :: Wr, Wi
      logical,       dimension(n)   :: done

      !Diagonalize the matrix and pickup the lambda=1 eigenvalues
      !The corresponding eigenvector contains the constraints of all moment components
      !Calling the general diagonalization subroutine from EisPack
      call Diagonalize_RGen(n,Ctr,wr,wi,.true.,zv)
      iss=0
      pti=0
      kval=0
      do i=1,n
        if(abs(wr(i)-1.0_dp) < epss .and. abs(wi(i)) < epss) then
          iss=iss+1   !Number of eigenvalues = 1 => number of free parameters
          pti(iss)=i !This points to the eigenvectors with eigenvalue equal to 1.
          zmi=1.0e6 !normalize the eigenvectors so that the minimum (non-zero value) is 1.
          j=1
          do k=1,n
            if(abs(zv(k,i)) < epss) cycle
            if(abs(zv(k,i)) < zmi) then
              zmi=abs(zv(k,i))
              kval=k  !This is the basis value
              j=nint(sign(1.0_dp,zv(k,i)))
            end if
          end do
          zv(:,i)=j*zv(:,i)/zmi  !This provides directly the multipliers for a single lambda=1 eigenvalue
          val(iss)=vect_val(kval) !This is the basis value to construct the new Moment
        end if
      end do
      codd="0"
      vect_out=0.0
      multi=0.0
      done=.false.
      where(abs(vect_val) < epss) done=.true.
      Select Case(iss)
        case(1)
          vect_out(:)=val(1)*zv(:,pti(1))
          where(abs(vect_out) > epss)  codd(:)=cdd(1)
          multi(:)=zv(:,pti(1))
        case(2:)
          ip=0
          do i=1,n
            if(.not. done(i)) then
              if(abs(vect_val(i)) > epss) then
                ip=ip+1
                codd(i)=cdd(ip)
                multi(i)=1.0
                vect_out(i)=vect_val(i)
                done(i)=.true.
                do j=i+1,n
                  if(.not. done(j)) then
                    if(abs(vect_val(i)-vect_val(j)) < epss) then
                      codd(j)=cdd(ip)
                      multi(j)=1.0
                      vect_out(j)=vect_val(i)
                      done(j)=.true.
                    else if(abs(vect_val(i)+vect_val(j)) < epss) then
                      codd(j)=cdd(ip)
                      multi(j)=-1.0
                      vect_out(j)=-vect_val(i)
                      done(j)=.true.
                    end if
                  end if
                end do
              end if
            end if
          end do
      End Select
    End Subroutine Get_Refinement_Codes

End Program Test_CIF_CFL