  Module DataRed_Treat_Reflections_Mod
    Use CFML_GlobalDeps
    Use CFML_Rational
    Use CFML_Maths,   only: zbelong,epss,sort,set_eps_math
    Use CFML_Metrics, only: Cell_Type, Cell_G_Type
    Use CFML_gSpaceGroups
    Use CFML_Reflections
    Use CFML_Propagation_Vectors

    Use Datared_Mod
    Use DataRed_rnw_Reflections_Mod
    Use Twin_Mod

    implicit none
    private
    public :: Treat_Reflections

    Interface Treat_Reflections
      module procedure Treat_Reflections_Conv
      module procedure Treat_Reflections_nonConv
    End Interface Treat_Reflections

   contains

    Subroutine Treat_Reflections_nonConv(mode,R,cond,cell,SpG,kinfo,tw,lun)
      character(Len=*),     intent(in)     :: mode
      type(Reflection_List),intent(in out) :: R
      type(Conditions_Type),intent(in out) :: cond
      type(Cell_G_type),    intent(in)     :: Cell
      class(SpG_type),      intent(in)     :: SpG
      type(kvect_Info_Type),intent(in)     :: kinfo
      type(Twin_Type),      intent(in)     :: tw
      integer,optional,     intent(in)     :: lun
    End Subroutine Treat_Reflections_nonConv

    Subroutine Treat_Reflections_Conv(R,cond,cell,SpG,kinfo,Gk,tw,lun)
      type(Reflection_List),intent(in out) :: R
      type(Conditions_Type),intent(in out) :: cond
      type(Cell_G_type),    intent(in)     :: Cell
      class(SpG_type),      intent(in)     :: SpG
      type(kvect_Info_Type),intent(in)     :: kinfo
      type(Group_k_Type),   intent(in)     :: Gk
      type(Twin_Type),      intent(in)     :: tw
      integer,optional,     intent(in)     :: lun

      ! Local Variables

      character(len=20)                    :: line1,wmess
      integer, dimension(R%nref)           :: itreat, iord, nequv, ini, fin, numor, ivk, icod, warn
      integer, dimension(nmax)             :: contr
      integer, dimension(3)                :: hkl
      real(kind=cp),    dimension(R%nref)  :: intav, sigmav, sigstat,lambda_laue, tbar,twotet,absorpt
      integer,          dimension(3,R%nref):: hkls
      real(kind=cp),    dimension(4)       :: angles
      real(kind=cp),    dimension(3)       :: h1,h2,h3
      real(kind=cp),    dimension(256)     :: weight  ! A maximum of 256 reflections equivalent to one of them can be treated
      real(kind=cp),    dimension(3,nmax)  :: hkln
      character(len=*),parameter,dimension(0:1) :: warn_mess=["                      ",  &
                                                               " <- Dubious reflection"]
      real    :: total, sig, suma, suman, sumaw, sumanw, Rint, Rwint, aver_sig, &
                 wavel,sigg, int_rej, epsg, delt,  scal_fact, q2, aver_int
      integer :: i,j,k, iou, ier, iv, ns, rej, ival,  nkv, nn, ihkl, irej,&
                 lenf, nin, ivp, nk, nequiv, L, drej, Lmin,i_refout
      logical :: absent,twin_acc

      warn=0
      iou=6
      if(present(lun)) iou=lun

      open(newunit=ihkl,file=trim(cond%fileout)//".int",status="replace",action="write")
      open(newunit=irej,file=trim(cond%fileout)//".rej",status="replace",action="write")

      if(cond%hkl_type == 10) then !The only data reduction applied concerns the symmetry ABSENCES
        write(unit=ihkl,fmt="(a)") cond%title
        write(unit=ihkl,fmt="(a)") "(3i4,2F14.4,i5,4f8.2)"
        write(unit=ihkl,fmt="(f9.5,a)") cond%wavel,"   0   1"
        rej=0
        do i=1,R%nref
          wmess=" "
          if(R%Ref(i)%idomain < 0) then
            R%Ref(i)%intens = -1.0
            R%Ref(i)%sigma  =  1.0
          end if
          if(h_absent(R%Ref(i)%h,SpG)) then
            rej=rej+1
            wmess="  <-  Forbidden reflection"
            write(unit=irej,fmt="(3i4,2F12.3,i5,f8.2,a)") &
            R%Ref(i)%h,R%Ref(i)%intens,R%Ref(i)%sigma,R%Ref(i)%twtheta,wmess
          end if
          write(unit=ihkl,fmt="(3i4,2F14.4,i5,4f8.2,a)") &
          R%Ref(i)%h,R%Ref(i)%intens,R%Ref(i)%sigma,abs(R%Ref(i)%idomain), R%Ref(i)%twtheta, 0.0, 0.0, 0.0,wmess
        end do
        write(unit=*,fmt="(/,a,i10)")  " => Number of reflections read               : ", R%nref
        write(unit=*,fmt="(a,i10)")    " => Number of rejected (absences) reflections: ", rej
        write(unit=*,fmt="(a)")        " => No averaging/merging reflections for HKL_TYPE=10 (HKLF5 Shelx) "
        write(unit=iou,fmt="(/,a,i10)")" => Number of reflections read               : ", R%nref
        write(unit=iou,fmt="(a,i10)")  " => Number of rejected (absences) reflections: ", rej
        write(unit=iou,fmt="(a)")      " => No averaging/merging reflections for HKL_TYPE=10 (HKLF5 Shelx) "
        return
      end if
      !
      ! All other cases for conventional case with no use of magnetic groups
      !
      !  First loop over reflections
      !
      nin=0
      itreat(:)=0
      ini(:)=0
      fin(:)=0
      total=0.0
      if(cond%transf_ind) then
        write(unit=iou,fmt="(/,a)") "    DETAILED TREATMENT OF REFLECTIONS USING TRANSFORMED INDICES AS GIVEN ABOVE"
        write(unit=iou,fmt="(a,/)") "    =========================================================================="
      else
        write(unit=iou,fmt="(/,a)") "    DETAILED TREATMENT OF REFLECTIONS USING INDICES AS IN DATA COLLECTION"
        write(unit=iou,fmt="(a,/)") "    ====================================================================="
      end if
      if(cond%prop) then
       write(unit=iou,fmt="(a)")  &
       "  No      hr      kr      lr     NEqv   Ini   Fin      <Inten>       Sigma      h   k   l   iv       Numor    2Theta"
      else if(cond%hkl_real) then
       write(unit=iou,fmt="(a)")  &
       "  No      hr      kr      lr    NEqv   Ini   Fin     <Inten>      Sigma       Numor    2Theta"
      else
       write(unit=iou,fmt="(a)")  &
       "      No   h   k   l  NEqv   Ini   Fin     <Inten>       Sigma       Numor    2Theta"
      end if

      rej=0
      i_refout=max(R%nref/60,100)

      !----------------------------------------------------
      do i=1,R%nref        !Loop over all measured reflections
          if(mod(i,i_refout) == 0) write(unit=*,fmt="(a,i10)") "=> Reflection:",i
          if(itreat(i) == 0) then   !If not yet treated do the following
            absent=.false.
            h1(:)=R%Ref(i)%hr

            if(cond%prop) then   !Propagation vector is given
             ivp=0
             hkl(:) = nint( h1 )
             h3(:)=hkl(:)-h1(:)
              !Verify first if the reflection is fundamental or is a satellite
              if(Zbelong(h3) .and. SpG%Num_Lat == 0) then
                absent=h_absent(nint(h1),SpG)    !Fundamental reflections
                ivp=0
              else
                do k=1, Gk%nk
                  h3(:)=  h1-Gk%stark(:,k)
                  if(Zbelong(h3)) then
                    ivp=k
                    if(Gk%k_equiv_minusk) then
                      hkl(:) = int( h3 )
                    else
                      hkl(:) = nint( h3(:))
                    end if
                    exit
                  end if
                end do
              end if
            else         !No Propagation vector is given
              h1=nint(h1)
              absent=h_absent(nint(h1),SpG)
            end if

            if(absent) then  !reject absent reflections
              rej=rej+1
              if(cond%prop .or. cond%hkl_real) then
               write(unit=irej,fmt="(3f8.4,2F12.3,4f8.2,i12)") &
               R%Ref(i)%hr,R%Ref(i)%intens,R%Ref(i)%sigma,R%Ref(i)%twtheta, R%Ref(i)%omega, R%Ref(i)%chi, R%Ref(i)%phi, R%Ref(i)%numor
              else
               write(unit=irej,fmt="(3i4,2F12.3,4f8.2,i12)") &
               R%Ref(i)%h,R%Ref(i)%intens,R%Ref(i)%sigma,R%Ref(i)%twtheta, R%Ref(i)%omega, R%Ref(i)%chi, R%Ref(i)%phi, R%Ref(i)%numor
              end if
              cycle
            end if
            nin=nin+1    !update the numer of independent reflections
            itreat(i)=i  !Make this reflection treated
            !  suma=R%Ref(i)%intens
            sig =1.0/R%Ref(i)%sigma**2
            ini(nin)=i   !put pointers for initial and final equivalent reflections
            fin(nin)=i
            nequv(nin)=1 !One reflection for the moment equivalent to itself
            if(cond%prop) then
              ivk(nin)=ivp
              hkls(:,nin)=hkl(:)
              write(unit=iou,fmt="(i8,3f8.4,a,2f12.3,a,3i4,i5,i12,f10.4)")&
                      i,h1,"                  ",R%Ref(i)%intens,R%Ref(i)%sigma,"   ",hkl,ivp,R%Ref(i)%numor,R%Ref(i)%twtheta
            else if(cond%hkl_real) then
              write(unit=iou,fmt="(i8,3f8.4,a,2f12.3,i12,f10.4,a)") &
                     i,h1,"                  ",R%Ref(i)%intens,R%Ref(i)%sigma,R%Ref(i)%numor,R%Ref(i)%twtheta,trim(warn_mess(warn(i)))
            else
              write(unit=iou,fmt="(i8,3i4,a,2f12.3,i12,f10.4)") &
                     i,nint(h1),"                  ",R%Ref(i)%intens,R%Ref(i)%sigma,R%Ref(i)%numor,R%Ref(i)%twtheta
            end if
            !-----------------------------------------------------------
            do j=i+1,R%nref       !look for equivalent reflections to the current (i) in the list
              if(abs(R%Ref(i)%s-R%Ref(j)%s) > 0.0002) exit
               if(cond%hkl_type == 7 .and. abs(R%Ref(i)%Lambda_laue-R%Ref(j)%Lambda_laue) > 0.1) cycle
               h2=R%Ref(j)%hr
               if(cond%prop) then
                 iv=0
                 hkl(:) = nint( h2 )
                  do k=1, Gk%nk
                    h3(:)=  h2-Gk%stark(:,k)
                    if(Zbelong(h3)) then
                     iv=k
                     hkl(:) = nint( h3(:))
                     exit
                    end if
                  end do
                 if(iv /= ivp) cycle

                 if(iv == 0) then

                    if(h_equiv(nint(h1),nint(h2),SpG,cond%Friedel)) then
                     itreat(j) = i
                     nequv(nin)=nequv(nin)+1
                     sig=sig + 1.0/R%Ref(j)%sigma**2
                     fin(nin)=j
                     write(unit=iou,fmt="(i8,3f8.4,a,2f12.3,a,3i4,i5,i12,f10.4)") &
                           j,h2,"                  ",R%Ref(j)%intens,R%Ref(j)%sigma,"   ",hkl,iv,R%Ref(j)%numor,R%Ref(j)%twtheta
                    end if

                 else

                    if(hk_equiv(h1,h2,Gk,cond%Friedel)) then
                     itreat(j) = i
                     nequv(nin)=nequv(nin)+1
                     sig=sig + 1.0/R%Ref(j)%sigma**2
                     fin(nin)=j
                     write(unit=iou,fmt="(i8,3f8.4,a,2f12.3,a,3i4,i5,i12,f10.4)") &
                           j,h2,"                  ",R%Ref(j)%intens,R%Ref(j)%sigma,"   ",hkl,iv,R%Ref(j)%numor,R%Ref(j)%twtheta
                    end if

                 end if

               else if (cond%hkl_real)  then
                 if(h_equiv(nint(h1),nint(h2),SpG,cond%Friedel)) then
                  itreat(j) = i
                  nequv(nin)=nequv(nin)+1
                  sig=sig + 1.0/R%Ref(j)%sigma**2
                  fin(nin)=j
                  write(unit=iou,fmt="(i8,3f8.4,a,2f12.3,i12,f10.4)") &
                        j,h2,"                  ",R%Ref(j)%intens,R%Ref(j)%sigma,R%Ref(j)%numor,R%Ref(j)%twtheta
                 end if
               else
                 if(h_equiv(nint(h1),nint(h2),Spg,cond%Friedel)) then ! if h1 eqv h2
                  itreat(j) = i                                 ! add h2 to the list equivalent to i
                  nequv(nin)=nequv(nin)+1                       ! update the number of equivalents
                  sig=sig + 1.0/R%Ref(j)%sigma**2
                  fin(nin)=j
                  write(unit=iou,fmt="(i8,3i4,a,2f12.3,i12,f10.4)")  &
                        j,nint(h2),"                  ",R%Ref(j)%intens,R%Ref(j)%sigma,R%Ref(j)%numor,R%Ref(j)%twtheta
                 end if
               end if


            end do

            ns=0
            do j=ini(nin),fin(nin)
              if(itreat(j) == i) then
                ns=ns+1
                weight(ns)=(1.0/R%Ref(j)%sigma**2)/sig
              end if
            end do

            suma=0.0
            ns=0
            do j=ini(nin),fin(nin)
              if(itreat(j) == i) then
                ns=ns+1
                suma=suma+weight(ns)*R%Ref(j)%intens
              end if
            end do


            intav(nin)=suma

            suma=0.0
            ns=0
            do j=ini(nin),fin(nin)
              if(itreat(j) == i) then
                ns=ns+1
                delt= intav(nin)-R%Ref(j)%intens
                if(abs(delt)/intav(nin) > cond%warning) warn(nin)=1
                suma=suma+weight(ns)*delt*delt
              end if
            end do

            sigmav(nin)=sqrt(suma)

            ! Average sigma using propagation error formula:
            !  I = Sum(v)/n => var(I)= 1/n^2 Sum(var(v)) => s(I)= sqrt(Sum(var(v))/n
            sigstat(nin)=0.0
            nn=0
            sigg=0.0
            do j=ini(nin),fin(nin)
              if(itreat(j) == i) then
                nn=nn+1
                sigg=sigg+R%Ref(j)%sigma
                sigstat(nin)=sigstat(nin)+R%Ref(j)%sigma*R%Ref(j)%sigma
              end if
            end do
            sigstat(nin) = sqrt(sigstat(nin))/real(nn)
            !
            sigg=sigg/real(nn)
            if(sigmav(nin) < sigg) sigmav(nin) = sigg
            ! Use statistical errors instead of experimental variance
            if(cond%statistics) sigmav(nin)=sigstat(nin)

            total=total+sigmav(nin)

            if(cond%prop .or. cond%hkl_real ) then
             write(unit=iou,fmt="(a,3f8.4,3i7,2f12.3,a/)")"   =>", h1,nequv(nin),ini(nin),fin(nin),intav(nin),sigmav(nin),&
                                                          trim(warn_mess(warn(nin)))
            else
            !hkl(:) = nint(h1(:))
            !ha=asu_hkl(hkl,grp_espacial)
             write(unit=iou,fmt="(a,3i4,i4,2i7,2f12.3,a/)")"   =>   ",nint(h1),nequv(nin),ini(nin),fin(nin),intav(nin),sigmav(nin),&
                                                           trim(warn_mess(warn(nin)))
            end if
          end if !itreat
      end do
      !
      ! Second loop over reflections to calculate R-int
      !
      ns=0
      suma =0.0
      suman=0.0
      sumaw =0.0
      sumanw=0.0
      do i=1,nin
        k=ini(i)
        if(nequv(i) < 2 ) cycle
        sig=0.0
        nn=0
        do j=ini(i),fin(i)
          if(itreat(j) == k) then
            nn=nn+1
            sig=sig+1.0/R%Ref(j)%sigma**2
          end if
        end do
        sig=1.0/sig
        do j=ini(i),fin(i)
          if(itreat(j) == k) then
            ns=ns+1
            suma=suma+abs(intav(i)-R%Ref(j)%intens)
            suman=suman+R%Ref(j)%intens
            sumaw=sumaw+ sig*((intav(i)-R%Ref(j)%intens)**2)/R%Ref(j)%sigma**2
            sumanw=sumanw+sig*(R%Ref(j)%intens/R%Ref(j)%sigma)**2
          end if
        end do
      end do
      Rint = 100.0*suma/max(1.0,suman)
      Rwint= 100.0*sqrt(sumaw/max(1.0,sumanw))
      aver_sig= total/real(nin)
      aver_int=suman/real(nin)
      nequiv=ns
      !
      !  Writing the list of averaged reflections
      !

      i=len_trim(cond%forma)
      cond%forma=cond%forma(1:i-1)//",a)"
      drej=0
      int_rej=0.0

      if(cond%powder) then

         write(unit=ihkl,fmt="(a)") cond%title
         write(unit=ihkl,fmt="(a)") "Intensities to be read from FullProf JBT=-3, IRF=2"
         do i=1,nin
           j=ini(i)
           hkl(:)=  Get_Asymm_Unit_H(R%Ref(j)%h,SpG)
           k=H_MULT(hkl,SpG,cond%Friedel)
           sigg=sqrt(intav(i))
           write(unit=ihkl,fmt="(3i4,i6,2F14.4)") hkl(:),k,sigg,0.5*sigmav(i)/sigg
         end do

      else  !Single X-tals

       if(cond%prop) then

          if(cond%hkl_type /= 7) then
            write(unit=ihkl,fmt="(a)") cond%title
            write(unit=ihkl,fmt="(a)") "(3i4,i5,2F14.4,i5,4f8.2)"
            if(cond%hkl_type /= 8) then
              write(unit=ihkl,fmt="(f9.5,a)") cond%wavel,"   0   0"
            else
              write(unit=ihkl,fmt="(f9.5,a)") cond%wavel,"   2   0"
            end if
          end if

          if(cond%domain) then
             if(.not. Gk%k_equiv_minusk) then
               nk= 2
               h1= Gk%stark(:,1)
               h2=-Gk%stark(:,1)
               if(cond%hkl_type /= 7) then
                 write(unit=ihkl,fmt="(i6,a)") nk, "  ! Number of Propagation vectors (Domain averaged: only k,-k)"
                 write(unit=ihkl,fmt="(i4,3f10.4)") 1, h1
                 write(unit=ihkl,fmt="(i4,3f10.4)") 2, h2
               end if
               nk=Gk%nk/2
             else
               h1=Gk%stark(:,1)
               write(unit=ihkl,fmt="(i6,a)") 1, "  ! Number of Propagation vectors (Domain averaged for all k in the star)"
               write(unit=ihkl,fmt="(i4,3f10.4)") 1, h1
               nk=Gk%nk
             end if


             do i=1,nin
               if(ivk(i) > nk) then  !belongs to -k (this not happens if k eqv -k)
                 ivk(i)=2
               else if(ivk(i) > 0) then   !belongs to  k
                 ivk(i)=1
               end if
             end do

             itreat(:) = 0
             ival=0
             do i=1,nin
              k=ini(i)
              if(itreat(i) == 0) then
               ns = 1
               suma = intav(i)
               sig  = sigmav(i)
               suman=0.0
                do j=i+1,nin
                    if(abs(R%Ref(k)%s-R%Ref(ini(j))%s) > 0.0002) exit
                    if( any(hkls(:,i) /= hkls(:,j)) ) cycle
                    if( ivk(i) /= ivk(j) ) cycle
                    ns=ns+1
                    suman=suman + abs(intav(i)-intav(j))
                    suma=suma + intav(j)
                    sig = sig + sigmav(j)
                    itreat(j)=j
                end do
                itreat(i)=i
                suma=real(nk)*suma/real(ns)
                suman=real(nk)*suman/real(ns)
                suman=100.0*suman/suma
                sig=real(nk)*sig/real(ns)
                ival=ival+1
                if(cond%hkl_type /= 7) then
                  write(unit=ihkl,fmt="(3i4,i5,2f14.4,i5,4f8.2,a,f8.2,i5,a)")    &
                  hkls(:,i),ivk(i), suma, sig,1, R%Ref(k)%twtheta, R%Ref(k)%omega, R%Ref(k)%chi, R%Ref(k)%phi,"     ", suman, ns,trim(warn_mess(warn(k)))
                else
                  write(unit=ihkl,fmt=cond%forma)    &
                  hkls(:,i),ivk(i), suma, sig,R%Ref(k)%icod, R%Ref(k)%Lambda_Laue, R%Ref(k)%twtheta,R%Ref(k)%absorpt,R%Ref(k)%tbar,trim(warn_mess(warn(k)))
                end if
              end if !itreat
             end do
             write(unit=*,fmt="(a,i6)") " => Number of domain-averaged reflections    : ",ival

          else
             write(unit=ihkl,fmt="(i6,a)") Gk%nk, "  ! Number of Propagation vectors (star_k)"
             do i=1, Gk%nk
               write(unit=ihkl,fmt="(i4,3f10.4)") i, Gk%stark(:,i)
             end do
             do i=1,nin
               j=ini(i)
               if(cond%hkl_type /= 7) then
                 write(unit=ihkl,fmt="(3i4,i5,2F14.4,i5,4f8.2,a)") &
                 hkls(:,i),ivk(i), intav(i),sigmav(i),1, R%Ref(j)%twtheta, R%Ref(j)%omega, R%Ref(j)%chi, R%Ref(j)%phi,trim(warn_mess(warn(i)))
               else
                  write(unit=ihkl,fmt=cond%forma)    &
                  hkls(:,i),ivk(i), intav(i),sigmav(i),R%Ref(k)%icod, R%Ref(k)%Lambda_Laue, R%Ref(k)%twtheta,R%Ref(k)%absorpt,R%Ref(k)%tbar,trim(warn_mess(warn(i)))
               end if
             end do

          end if

       else  !no propagation vectors

          if(cond%hkl_type /= 7) then
            write(unit=ihkl,fmt="(a)") cond%title
            if(cond%hkl_type == 8) then
              write(unit=ihkl,fmt="(a)") "(3i4,2f10.6,i5,f8.5,3f8.2,a)"
            else
              write(unit=ihkl,fmt="(a)") "(3i4,2F14.4,i5,4f8.2)"
            end if
          end if

          if(cond%twinned) then

            if(cond%hkl_type /= 7) write(unit=ihkl,fmt="(f9.5,a)") cond%wavel,"   0   1"
            write(unit=*,fmt="(/,a)")   " => Indices have been transformed according to TWIN commands "
            write(unit=iou,fmt="(/,a)") " => Indices have been transformed according to TWIN commands "
            do i=1,nin
              j=ini(i)
              angles(1)= R%Ref(j)%twtheta
              angles(2)= R%Ref(j)%omega
              angles(3)= R%Ref(j)%chi
              angles(4)= R%Ref(j)%phi
              call get_domain_contrib(R%Ref(j)%hr,tw,hkln,contr,angles,SpG,Cell)
              twin_acc=.false.
              Lmin=0
              do L=1,tw%nmat
               if(contr(L) == 1) then  !Handles the possibility that the initial reflection has not
                 Lmin=L                !integer indices
                 twin_acc=.true.
                 exit
               end if
              end do
              do L=tw%nmat,Lmin+1,-1
                if(contr(L) == 1) then
                   write(unit=ihkl,fmt="(3i4,2F14.4,i5)") nint(hkln(:,L)),-1.0,0.0,L
                   twin_acc=.true.
                end if
              end do
              if(.not. twin_acc .or. Lmin == 0) then
                 if(drej == 0) write(unit=irej,fmt="(/,a,/)") " => LIST OF REJECTED REFLECTIONS FOR DOMAIN CONTRIBUTION:"
                 int_rej=int_rej+ intav(i)
                 drej=drej+1
                 write(unit=irej,fmt="(/,a,3f6.2,2f10.2,a)") &
                 " => Reflection:", R%Ref(j)%hr, intav(i),sigmav(i), " has been rejected for the present TWIN LAW"
                 do L=1,tw%nmat
                  if(contr(L) == 0) then
                   write(unit=irej,fmt="(a,i3,a,3f6.2)")   &
                   "                Domain:",L,"    -> ",hkln(:,L)
                  else
                   write(unit=irej,fmt="(a,i3,a,3f6.2,a)") &
                   "                Domain:",L,"    -> ",hkln(:,L)," Forbidden in "//trim(tw%twin_SpG)
                  end if
                 end do
              end if
              if(Lmin > 0) write(unit=ihkl,fmt="(3i4,2f14.4,i5,4f8.2,a)") &
                           nint(hkln(:,Lmin)),intav(i),sigmav(i),Lmin, angles(:),trim(warn_mess(warn(i)))
            end do
             if(drej > 0) int_rej=int_rej/real(drej)

          else

            if(cond%hkl_type /= 7) then
              if(cond%hkl_type /= 8) then
                write(unit=ihkl,fmt="(f9.5,a)") cond%wavel,"   0   0"
              else                            !23456789    2       0      0.9400    0.9400    0      0       0
                write(unit=ihkl,fmt="(a)")   "! Lambda   itypd    ipow    Polarp    Polarm  UB_mat  Ext  read_strf"
                write(unit=ihkl,fmt="(f9.5,a)") cond%wavel,"    2       0      0.9400    0.9400    0      0       0"
              end if
            end if

            do i=1,nin
              j=ini(i)
              hkl(:)=  nint(R%Ref(j)%hr)
              if(cond%hkl_type /= 7) then
                if(cond%hkl_type == 8) then
                  !                 gamma    omega    nu
                  !call z1frnb(cond%wavel,R%Ref(j)%phi, R%Ref(j)%omega,R%Ref(j)%chi,h1)
                  h1=z1frnb(cond%wavel,R%Ref(j)%phi, R%Ref(j)%omega,R%Ref(j)%chi)
                  q2=1.0 - h1(3)*h1(3)/dot_product(h1,h1)
                  write(unit=ihkl,fmt="(3i4,2f10.6,i5,f8.5,3f8.2,a)") hkl(:),intav(i),sigmav(i),1, q2, 0.0, &
                                                                 0.0, 0.0,trim(warn_mess(warn(i)))
                else
                  hkl=Get_Asymm_Unit_H(hkl,SpG)
                  write(unit=ihkl,fmt="(3i4,2f14.4,i5,4f8.2,a)") hkl(:),intav(i),sigmav(i),1, R%Ref(j)%twtheta, R%Ref(j)%omega, &
                                                                 R%Ref(j)%chi, R%Ref(j)%phi,trim(warn_mess(warn(i)))
                end if
              else
                  write(unit=ihkl,fmt=cond%forma)    &
                  hkl(:), intav(i),sigmav(i),icod(j), R%Ref(j)%Lambda_Laue, R%Ref(j)%twtheta,R%Ref(j)%absorpt,R%Ref(j)%tbar,trim(warn_mess(warn(i)))
              end if
            end do

          end if

       end if

      end if

      write(unit=iou,fmt="(/,a,i10)")" => Number of reflections read               : ", R%nref
      write(unit=iou,fmt="(a,i10)")  " => Number of valid independent   reflections: ", nin
      write(unit=iou,fmt="(a,i10)")  " => Number of obs. with equival.  reflections: ", ns
      write(unit=iou,fmt="(a,i10)")  " => Number of rejected (absences) reflections: ", rej
      write(unit=iou,fmt="(a,f10.2)")" => R-internal for equivalent reflections (%): ", Rint
      write(unit=iou,fmt="(a,f10.2)")" => R-weighted for equivalent reflections (%): ", Rwint
      write(unit=iou,fmt="(a,f10.2)")" => Average Intensity of reflections         : ", aver_int
      write(unit=iou,fmt="(a,f10.2)")" => Average sigma for equivalent reflections : ", aver_sig


      write(unit=*,fmt="(/,a,i10)")" => Number of reflections read               : ", R%nref
      write(unit=*,fmt="(a,i10)")  " => Number of valid independent   reflections: ", nin
      write(unit=*,fmt="(a,i10)")  " => Number of obs. with equival.  reflections: ", nequiv
      write(unit=*,fmt="(a,i10)")  " => Number of rejected (absences) reflections: ", rej
      write(unit=*,fmt="(a,f10.2)")" => R-internal for equivalent reflections (%): ", Rint
      write(unit=*,fmt="(a,f10.2)")" => R-weighted for equivalent reflections (%): ", Rwint
      write(unit=*,fmt="(a,f10.2)")" => Average Intensity of reflections         : ", aver_int
      write(unit=*,fmt="(a,f10.2)")" => Average sigma for equivalent reflections : ", aver_sig
      if(cond%twinned) then
         write(unit=iou,fmt="(a,i10)")    " => Number of domain rejected reflections    : ", drej
         write(unit=*  ,fmt="(a,i10)")    " => Number of domain rejected reflections    : ", drej
         write(unit=iou,fmt="(a,f10.2)")  " => Average intensity of domain rejected ref.: ", int_rej
         write(unit=*,  fmt="(a,f10.2)")  " => Average intensity of domain rejected ref.: ", int_rej
      end if


    End Subroutine Treat_Reflections_Conv

  End Module DataRed_Treat_Reflections_Mod