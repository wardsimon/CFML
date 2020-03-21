!!----
!!----
!!----
SubModule (CFML_Reflections) Generation_of_general_Reflections
   Contains
   !!----
   !!---- Gener_Reflections
   !!----    Calculate unique reflections between two values of
   !!----    sin_theta/lambda.  The output is not ordered.
   !!----
   !!---- 21/06/2019
   !!
   Module Subroutine Gener_Reflections(Cell,Sintlmax,Mag,Num_Ref,Reflex,SpG,kinfo,order,powder,mag_only,Friedel)
      !---- Arguments ----!
      class(Cell_G_Type),                          intent(in)     :: Cell
      real(kind=cp),                               intent(in)     :: Sintlmax
      logical,                                     intent(in)     :: Mag
      integer,                                     intent(out)    :: Num_ref
      class(Refl_Type), dimension(:), allocatable, intent(out)    :: Reflex
      class(Spg_Type) ,              optional,     intent(in)     :: SpG
      type(kvect_info_type),         optional,     intent(in)     :: Kinfo
      character(len=*),              optional,     intent(in)     :: Order
      logical,                       optional,     intent(in)     :: Powder
      logical,                       optional,     intent(in)     :: Mag_only
      logical,                       optional,     intent(in)     :: Friedel

      !---- Local variables ----!
      real(kind=cp)         :: epsr=0.00001, delt=0.0001
      real(kind=cp)         :: sval,max_s !,vmin,vmax
      integer               :: h,k,l,hmin,kmin,lmin,hmax,kmax,lmax, maxref,i,j,indp,indj, &
                               maxpos, mp, iprev,Dd, nf, ia, i0, nk, nharm,n
      integer,      dimension(:),   allocatable :: hh,kk,nulo
      integer,      dimension(:,:), allocatable :: hkl,hklm
      integer,      dimension(:),   allocatable :: indx,indtyp,ind,ini,fin,itreat

      real(kind=cp),dimension(:),   allocatable :: sv,sm
      logical                                   :: kvect,ordering,magg,Frd

      !> Init
      Dd=3
      ordering=.false.
      kvect=.false.
      magg=.false.
      Frd=.true.
      if(present(Friedel)) Frd=Friedel
      if (present(mag_only)) magg=mag_only
      if (present(order) .or. present(powder)) ordering=.true.
      if (present(kinfo)) then
         nk=kinfo%nk
         nharm=kinfo%nq
         kvect=.true.
      end if
      if (kvect) Dd=3+nk ! total dimension of the reciprocal space

      hmax=nint(Cell%cell(1)*2.0*sintlmax+1.0)
      kmax=nint(Cell%cell(2)*2.0*sintlmax+1.0)
      lmax=nint(Cell%cell(3)*2.0*sintlmax+1.0)
      hmin=-hmax; kmin=-kmax; lmin= -lmax
      maxref= (2*hmax+1)*(2*kmax+1)*(2*lmax+1)
      if (kvect) then
         do k=1,nk
            maxref=maxref*2*nharm
         end do
         if (present(kinfo)) then
            max_s=maxval(kinfo%sintlim)
         else
            max_s=sintlmax
         end if
      end if

      allocate(hkl(Dd,maxref), indx(maxref), indtyp(maxref), ind(maxref), sv(maxref))
      allocate(hh(Dd), kk(Dd), nulo(Dd))
      nulo=0
      indtyp=0
      num_ref=0

      !> Generation of fundamental reflections
      i0=0
      ext_do: do h=hmin,hmax
         do k=kmin,kmax
            do l=lmin,lmax
               hh=0
               hh(1:3)=[h,k,l]
               sval=H_S(hh,Cell)
               if (sval > sintlmax) cycle

               num_ref=num_ref+1
               if (num_ref > maxref) then
                  num_ref=maxref
                  exit ext_do
               end if

               if (sval < epsr) i0=num_ref !localization of 0 0 0  reflection
               sv(num_ref)=sval
               hkl(:,num_ref)=hh
            end do
         end do
      end do ext_do

      !> Generation of satellites
      !> The generated satellites corresponds to those obtained from the list
      !> of +/-kinfo%q_qcoeff
      nf=num_ref
      if (kvect) then
         do_ex: do n=1,kinfo%nq
            do i=1,nf
               hh=hkl(:,i)
               do ia=-1,1,2
                  hh(4:3+nk)=ia*kinfo%q_coeff(1:nk,n)
                  sval=H_S(hh, Cell, nk, kinfo%kv)
                  if (sval > max_s) cycle

                  num_ref=num_ref+1
                  if (num_ref > maxref) then
                     num_ref=maxref
                     exit do_ex
                  end if

                  sv(num_ref)=sval
                  hkl(:,num_ref)=hh
               end do !ia
            end do  !i
         end do do_ex
      end if

      !> elimination of the reflection (0000..)
      do i=i0+1,num_ref
         sv(i-1)=sv(i)
         hkl(:,i-1)=hkl(:,i)
      end do
      num_ref=num_ref-1

      !> Determination of reflection character and extinctions (Lattice + others)
      n=0
      do i=1,num_ref
         hh=hkl(:,i)
         mp=0
         if (present(SpG)) then
            if (SpG%Num_Lat /= 0) then
               if (H_Latt_Absent(hh,SpG%Lat_tr,SpG%Num_Lat)) cycle
            end if
            if (H_Absent(hh,SpG)) then
               if (SpG%Mag_type /= 2 .and. Mag) then
                  if (mH_Absent(hh,SpG)) then
                     cycle
                  else
                     mp=1   !pure magnetic
                  end if
               else
                  cycle
               end if
            else
               if (SpG%Mag_type /= 2 .and. Mag) then
                  if (mH_Absent(hh,SpG)) then
                     mp=0  !pure nuclear
                  else
                     mp=2
                  end if
               end if
            end if
         end if

         n=n+1
         hkl(:,n)=hkl(:,i)
         sv(n)=sv(i)
         indtyp(n)=mp
      end do
      num_ref=n

      if (ordering) then
         indx=sort(sv,num_ref)
         allocate(hklm(Dd,num_ref),sm(num_ref))
         do i=1,num_ref
            j=indx(i)
            hklm(:,i)=hkl(:,j)
            sm(i)=sv(j)
            ind(i)=indtyp(j)
         end do

         !> contains now the type of reflection in the proper order
         indtyp(1:num_ref)=ind(1:num_ref)
         hkl(:,1:num_ref)=hklm(:,1:num_ref)
         sv(1:num_ref)=sm(1:num_ref)
         if (present(SpG) .and. present(powder)) then
            deallocate(hkl,sv,indx,ind)
            allocate(ini(num_ref),fin(num_ref),itreat(num_ref))
            itreat=0; ini=0; fin=0
            indp=0
            do i=1,num_ref       !Loop over all reflections
               if (itreat(i) == 0) then   !If not yet treated do the following
                  hh(:)=hklm(:,i)
                  indp=indp+1  !update the number of independent reflections
                  itreat(i)=i  !Make this reflection treated
                  ini(indp)=i  !put pointers for initial and final equivalent reflections
                  fin(indp)=i
                  do j=i+1,num_ref  !look for equivalent reflections to the current (i) in the list
                     if (abs(sm(i)-sm(j)) > delt) exit
                     kk=hklm(:,j)
                     if (h_equiv(hh,kk,SpG,Frd)) then ! if  hh eqv kk (Friedel law according to Frd)
                        itreat(j) = i                 ! add kk to the list equivalent to i
                        fin(indp) = j
                     end if
                  end do
               end if !itreat
            end do

            !> Selection of the most convenient independent reflections
            allocate(hkl(Dd,indp),sv(indp),ind(indp))
            do i=1,indp
               maxpos=0
               indj=ini(i)
               iprev=itreat(indj)
               do j=ini(i),fin(i)
                  if (iprev /= itreat(j)) cycle
                  hh=hklm(:,j)
                  mp=count(hh > 0)
                  if (mp > maxpos) then
                     indj=j
                     maxpos=mp
                  end if
               end do !j
               hkl(:,i)=hklm(:,indj)
               if (hkl(1,i) < 0) hkl(:,i)=-hkl(:,i)
               sv(i)=sm(indj)
               ind(i)=indtyp(indj)
            end do
            indtyp(1:indp)=ind(1:indp)
            num_ref=indp
         end if  !SpG and Powder

      end if !present "order"

      !> Final assignments
      if (allocated(reflex)) deallocate(reflex)
      allocate(reflex(num_ref))

      do i=1,num_ref
         allocate(reflex(i)%h(Dd)) !needed in f95
         reflex(i)%h      = hkl(:,i)
         if(dd > 3) then
           kk               = abs(hkl(4:3+nk,i))
           reflex(i)%pcoeff = 0 !Fundamental reflections point to the Fourier coefficient [00...]
           do_n: do n=1,kinfo%nq
              do k=1,nk
                 if (equal_vector(kk(1:nk),abs(kinfo%q_coeff(1:nk,n))))  then
                    reflex(i)%pcoeff=n
                    exit do_n
                 end if
              end do
           end do do_n
         end if
         reflex(i)%s      = sv(i)
         if (present(SpG)) then
            reflex(i)%mult = h_mult(reflex(i)%h,SpG,Frd)
         else
           reflex(i)%mult = 1
         end if
         reflex(i)%imag = indtyp(i)
      end do
   End Subroutine Gener_Reflections

End SubModule Generation_of_general_Reflections