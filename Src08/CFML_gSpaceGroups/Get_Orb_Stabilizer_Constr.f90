!!----
SubModule (CFML_gSpaceGroups) Stabilizer_Constraints
   implicit none

    character (len=*), dimension(26),parameter   :: &
    cdd=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r', &
         's','t','u','v','w','x','y','z']
   Contains

   !!----
   !!---- GET_STABILIZER
   !!----    Subroutine to obtain the list of symmetry operator of a space group that leaves
   !!----    invariant an atomic position. This subroutine provides a pointer to the symmetry
   !!----    operators of the site point group and the additional translation with respect to
   !!----    the canonical representant.
   !!----
   !!---- 13/06/2019
   !!
   Module Subroutine Get_Stabilizer(X, Spg,Order,Ptr,Atr)
      !---- Arguments ----!
      real(kind=cp), dimension(3),  intent (in)  :: x     ! real space position (fractional coordinates)
      class(Spg_Type),              intent (in)  :: Spg   ! Space group
      integer,                      intent(out)  :: order ! Number of sym.op. keeping invariant the position x
      integer, dimension(:),        intent(out)  :: ptr   ! Array pointing to the symmetry operators numbers
                                                          ! of the stabilizer of x
      real(kind=cp), dimension(:,:),intent(out)  :: atr   ! Associated additional translation to the symmetry operator

      !---- Local variables ----!
      real(kind=cp), dimension(3)    :: xx, tr
      integer                        :: j,n1,n2,n3

      !> Init
      order = 1              ! Identity belongs always to the stabilizer
      ptr   = 0; ptr(1)= 1
      atr   = 0.0_cp

      do n1=-1,1
         do n2=-1,1
            do n3=-1,1
               tr=real([n1, n2, n3])
               do j=2,Spg%multip
                  xx=Apply_OP(Spg%Op(j),x) + tr - x
                  if (sum(abs(xx)) > 2.0 * EPSS) cycle
                  order=order+1
                  ptr(order)=j
                  atr(:,order)=tr
               end do
            end do
         end do
      end do

   End Subroutine Get_Stabilizer

   Module Subroutine Get_Orbit(x,mom,Spg,Mult,orb,morb,ptr,convl)
      !---- Arguments ----!
      real(kind=cp), dimension(:),                intent (in) :: x,mom
      class(SpG_Type),                            intent (in) :: spg
      integer,                                    intent(out) :: mult
      real(kind=cp),dimension(:,:), allocatable,  intent(out) :: orb,morb
      integer, dimension(:),allocatable, optional,intent(out) :: ptr
      logical,                           optional,intent(in)  :: convl

      !---- Local variables ----!
      integer                                          :: i, j, nt,d
      real(kind=cp), dimension(Spg%d)                  :: xs,xsp
      real(kind=cp), dimension(Spg%d)                  :: ms,msp
      real(kind=cp), dimension(Spg%d-1)                :: v
      real(kind=cp), dimension(Spg%d,Spg%d,Spg%multip) :: Om
      logical                                          :: conv

      conv=.true.
      if(present(convl)) conv=convl
      d=Spg%d-1
      xs(Spg%d)=1.0
      ms(Spg%d)=1.0
      allocate(orb(d,Spg%multip),morb(d,Spg%multip),ptr(Spg%multip))
      orb=0.0; morb=0.0

      Select Type(SpG)

        type is (SuperSpaceGroup_Type)
           Om=SpG%Om
           xs(1:3)=x   !Extend the position and moment to superspace
           ms(1:3)=mom
           do i=1,SpG%nk
             xs(3+i)=dot_product(x,SpG%kv(:,i))
             ms(3+i)=dot_product(mom,SpG%kv(:,i))
           end do
        class default
           do i=1,Spg%Multip
             Om(:,:,i)=Spg%Op(i)%Mat
           end do
           xs(1:d)=x
           ms(1:d)=mom
      End Select

      mult=1
      orb(:,1)=xs(1:d)
      morb(:,1)=ms(1:d)
      if(present(ptr)) ptr(mult) = 1

      do_ext: do j=2,Spg%Multip
         xsp=matmul(Om(:,:,j),xs)
         xsp=modulo_lat(xsp)
         do nt=1,mult
            v=orb(1:d,nt)-xsp(1:d)
            if(Spg%Num_lat > 0 .and. .not. conv) then
              if (is_Lattice_vec(v,Spg%Lat_tr)) cycle do_ext
            else
              if (Zbelong(v)) cycle do_ext
            end if
         end do
         msp(1:d)=Spg%Op(j)%dt*Spg%Op(j)%time_inv*matmul(Om(1:d,1:d,j),ms(1:d))
         mult=mult+1
         orb(:,mult)=xsp(1:d)
         morb(:,mult)=msp(1:d)
         if(present(ptr)) ptr(mult) = j   !Pointer to symmetry operator
      end do do_ext

   End Subroutine Get_Orbit

   Module Subroutine Get_moment_ctr(xnr,moment,Spg,codini,codes,ord,ss,att,Ipr,ctr_code)
      real(kind=cp), dimension(3),            intent(in)     :: xnr
      real(kind=cp), dimension(:),            intent(in out) :: moment
      class(SpG_type),                        intent(in)     :: Spg
      Integer,                                intent(in out) :: codini
      real(kind=cp), dimension(:),            intent(in out) :: codes
      integer,                       optional,intent(in)     :: ord
      integer, dimension(:),         optional,intent(in)     :: ss
      real(kind=cp), dimension(:,:), optional,intent(in)     :: att
      integer,                       optional,intent(in)     :: Ipr
      character(len=*),              optional,intent(out)    :: ctr_code

      ! Local variables
      character(len=1),  dimension(3)   :: codd
      character(len=15), dimension(3)   :: St_Cod
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
      if(present(ctr_code)) ctr_code="(0,0,0)"
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

      if(present(ctr_code)) then
          do j=1,3
            St_cod(j)=" "
            if(abs(multi(j)) < 0.00001) then
              St_cod(j) = "0"
            else
              if(multi(j) == 1.0_cp) then
                 St_cod(j)=codd(j)
              else if(multi(j) == -1.0_cp) then
                 St_cod(j)="-"//codd(j)
              else
                 write(unit=St_cod(j),fmt="(f10.5,a)")  multi(j),"*"//codd(j)
                 St_cod(j)=adjustl(St_cod(j))
              end if
            end if
          end do
         write(unit=ctr_code,fmt="(28a)") " ( ",(St_cod(j)//", ",j=1,n-1),St_cod(j)//" )"
         ctr_code=pack_string(ctr_code)
      end if
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
     integer                       :: i,j,k,kval,ip !,ier
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
         zv(1:n,i)=j*zv(1:n,i)/zmi  !This provides directly the multipliers for a single lambda=1 eigenvalue
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
         vect_out(1:n)=val(1)*zv(1:n,pti(1))
         where(abs(vect_out) > epss)  codd(:)=cdd(1)
         multi(1:n)=zv(1:n,pti(1))
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

   !!
   !!----  Subroutine get_moment_ctr(xnr,TFourier,Spg,codini,codes,ord,ss,att,Ipr,ctr_code)
   !!----     real(kind=cp), dimension(3),            intent(in    ) :: xnr    !Atom position (fractional coordinates)
   !!----     Complex(kind=cp), dimension(3),         intent(in out) :: TFourier !Fourier coefficients at position xnr
   !!----     class(SuperSpaceGroup_Type),            intent(in)     :: Spg    !Super Space Group
   !!----     Integer,                                intent(in out) :: codini !Last attributed parameter
   !!----     real(kind=cp), dimension(3),            intent(in out) :: codes  !codewords for positions
   !!----     integer,                       optional,intent(in)     :: ord
   !!----     integer, dimension(:),         optional,intent(in)     :: ss
   !!----     real(kind=cp), dimension(:,:), optional,intent(in)     :: att
   !!----     integer,                       optional,intent(in)     :: Ipr
   !!----     character(len=*),              optional,intent(out)    :: ctr_code
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
   Module Subroutine Get_TFourier_ctr(xnr,TFourier,codes,SpG,codini,mode,ord,ss,att,Ipr,ctr_code)
      real(kind=cp), dimension(3),            intent(in)     :: xnr
      real(kind=cp), dimension(:,:),          intent(in out) :: TFourier
      real(kind=cp), dimension(:,:),          intent(in out) :: codes
      class(SuperSpaceGroup_Type),            intent(in)     :: SpG
      Integer,                                intent(in out) :: codini
      character(len=*),                       intent(in)     :: mode
      integer,                       optional,intent(in)     :: ord
      integer, dimension(:),         optional,intent(in)     :: ss
      real(kind=cp), dimension(:,:), optional,intent(in)     :: att
      integer,                       optional,intent(in)     :: Ipr
      character(len=*),dimension(:), optional,intent(out)    :: ctr_code

      ! Local variables
      integer,          dimension(SpG%nk)     :: brack_m  ! hs=(H,[m])
      real(kind=cp),    dimension(spg%d-1)    :: ts  ! ts=(t,tI)
      integer,          dimension(SpG%nk)     :: mE  ![m].E
      integer,          dimension(SpG%nk,3)   :: Mx  !M
      integer,          dimension(3,3)        :: g, magm   !g, magm= delta * det(g) * g
      integer,      dimension(SpG%nk,SpG%nk)  :: ep  !E
      character(len=1),  dimension(26)        :: codd
      character(len=15), dimension(26)        :: st_cod
      integer                                 :: i,j,order,iq,ir,d,n,k,iqt,ig,ip,iss,nq
      real(kind=cp)                           :: suma,alpha
      integer,           dimension(48)        :: ss_ptr
      real(kind=cp),     dimension(3,48)      :: atr
      real(kind=cp),     dimension(6*SpG%nq)  :: cod,multi,val
      real(kind=cp),     dimension(3)         :: x
      real(kind=dp),     dimension(6,6)       :: subm
      real(kind=dp),  dimension(6*SpG%nq,6*SpG%nq) :: Ctr,sCtr
      real(kind=cp),     parameter                 :: epss=0.001_cp
      real(kind=cp), dimension(6*SpG%nq)           ::TFourL,sTF


      !Test if all codes are given ... in such a case the user constraints
      !are prevalent

      d=SpG%d-1
      suma=0.0
      !iq=0
      n=0
      nq=size(TFourier,dim=2)
      do k=1,nq
        do j=1,6
           n=n+1
           suma=suma+abs(codes(j,k))
        end do
      end do

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

      do ip=1,nq
        i=(ip-1)*6
        TFourL(i+1:i+6)=TFourier(:,ip)
      end do

      sCtr=0.0_cp
      sTF=0.0
      if(order > 1) then
        n=6*nq

        do ig=1,order
          ir=ss_ptr(ig)
              g(:,:) = SpG%Op(ir)%Mat(1:3,1:3)   !                          /  g    0   t  \
               ts(:) = SpG%Op(ir)%Mat(1:d,d+1)   !   Superspace operator:  |  Mx   ep   tI |    !ts=(t,tI)
             Mx(:,:) = SpG%Op(ir)%Mat(4:d,1:3)   !                         \   0    0   1 /
             Ep(:,:) = SpG%Op(ir)%Mat(4:d,4:d)
             magm(:,:) = g(:,:)
             if(mode(1:1) == "M")  magm(:,:) = magm(:,:)*SpG%Op(ir)%time_inv*SpG%Op(ir)%dt
           !Identify the q_coefficients of mM to select the proper Tfourier
           !Construction of the Cr-matrix
           Ctr=0.0
           do ip=1,nq
             brack_m(:)=SpG%q_coeff(:,ip) !nk-components vector
             mE=matmul(brack_m,Ep)   ![m].Ep  nk-components vector
             do iqt=1,nq
               iq=iqt
               if(equal_vector(mE,SpG%q_coeff(:,iqt))) then
                 iss=1
                 exit
               end if
               if(equal_vector(mE,-SpG%q_coeff(:,iqt))) then
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
             Tfourl(i+1:i+6)=matmul(subm,TFourier(:,iq))
             sTF(i+1:i+6)=sTF(i+1:i+6)+TfourL(i+1:i+6)
           end do !ip
           sCtr=sCtr+Ctr !Adding constraint matrices for each operator of stabilizer
           if(present(ipr)) then
             do ip=1,nq
                i=(ip-1)*6
                write(unit=ipr,fmt='(a,i2,a,t20,a,t55,6f14.4,4i3)') '     Operator ',ir,": ",trim(Spg%Symb_Op(ir)),Tfourl(i+1:i+6),brack_m(:)
             end do
           end if
        end do  !ig operators
        sCtr=sCtr/order
        sTF=sTF/order
        if(present(ipr)) then
          do ip=1,nq
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
        do ip=1,nq
          i=(ip-1)*6
          TFourier(:,ip)=TFourL(i+1:i+6)
        end do

        codes=0.0; j=0
        do ip=1,nq
          do i=1,6
            j=j+1
            if(abs(multi(j)) > epss)  codes(i,ip) = sign(1.0_cp, multi(j))*(abs(cod(j))*10.0_cp + abs(multi(j)) )
          end do
        end do
        codini=codini+iss
        if(present(Ipr)) then
          Write(Ipr,"(a,i4)")      " Number of free parameters: ",iss
          Write(Ipr,"(a,24F14.6)") " Basic Values: ",val(1:iss)
          write(Ipr,"(a,24f14.6)") " Multipliers: ",(multi(j), j=1,n)
          write(Ipr,"(28a)")       " String with free parameters: ( ",(codd(j)//", ",j=1,n-1),codd(n)//" )"
          write(Ipr,"(a,24i6)")    " Resulting integer codes: ", nint(cod(1:n))
          do ip=1,nq
            write(Ipr,"(a,i2, 6f14.6)") " Constrained Fourier Coefficients: ",ip,TFourier(:,ip)
            write(Ipr,"(a,tr2,6f14.6)") " Final codes: ",codes(:,ip)
          end do
        end if

      else !No restrictions

        codd(1:n)=cdd(1:n)
        multi(1:n)=1.0_cp
        j=0
        do ip=1,nq
          do i=1,6
            j=j+1
            cod(j)=codini+j
            codes(i,ip) = abs(cod(j))*10.0_cp + abs(multi(j))
          end do
        end do
        codini=codini+n
        if(present(Ipr)) then
          write(Ipr,"(a,24f10.6)") " General position, no constraints in Fourier Coefficients "
          write(Ipr,"(28a)")       " String with free parameters: ( ",(codd(j)//", ",j=1,n-1),codd(n)//" )"
          write(Ipr,"(a,24i6)")    " Resulting integer codes: ", nint(cod(1:n))
          do ip=1,nq
            write(Ipr,"(a,i2, 6f14.6)") " Constrained Fourier Coefficients: ",ip,TFourier(:,ip)
            write(Ipr,"(a,tr2,6f14.6)") " Final codes: ",codes(:,ip)
          end do
        end if

      end if

      if(present(ctr_code)) then
        St_cod="0"
        n=0
        do ip=1,nq
          do j=1,6
            n=n+1
            St_cod(n)=" "
            if(abs(multi(n)) < 0.00001) then
              St_cod(n) = " 0 "
            else
              if(multi(n) == 1.0_cp) then
                 St_cod(n)=codd(n)
              else if(multi(n) == -1.0_cp) then
                 St_cod(n)="-"//codd(n)
              else
                 write(unit=St_cod(n),fmt="(f10.5,a)")  multi(n),"*"//codd(n)
                 St_cod(n)=adjustl(St_cod(n))
              end if
            end if
          end do
        end do
        n=0
        do ip=1,nq
          write(unit=ctr_code(ip),fmt="(8a)") " ( ",(St_cod(n+j)//", ",j=1,5),St_cod(n+6)//" )"
          ctr_code(ip)=pack_string(ctr_code(ip))
          n=n+6
        end do
      end if

   End Subroutine Get_TFourier_ctr

End SubModule Stabilizer_Constraints