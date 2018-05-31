  Module CFML_FlipR_Mod

      use CFML_GlobalDeps,            only: cp,dp,pi,tpi
      use CFML_Math_General,          only: factorial
      !use fullprof_params,            only: nat_p, xl, ac, ptr, Spgr, asfa, bsfa, irf, &
      !                                      jfou,jtyp, ipow, maxs, n_t, npsc, dfp, dfpp, &
      !                                      isap, icryg, deriv, phasen, nscat, a, lp,  &
      !                                      ntyp, refs, atext,MaxPOINT,MaxREFLX,bk,ggi,cci,  &
      !                                      cmat,  gam, sp_pointers, den_nlm, &
      !                                       polarp, polarm, symb, &
      !                                      read_strf, icod, par
      !
      USE CFML_Spherical_Harmonics,   only: Real_Spher_Harm_Ucvec,Real_Spher_HarmCharge_Ucvec, Int_Slater_Bessel
      USE CFML_Reflections_Utilities, only: unit_cart_hkl
      USE CFML_Crystal_Metrics,       only: cart_vector
      USE CFML_Extinction_Correction, only: Correct_FlippingRatios
      USE CFML_IO_Messages ,          only: Error_Message

      implicit none
      private
      Public  :: Flipr, simpleFlipr, Schwinger_Amplitude
      !Changed to public for accessing directly to the form-factors by external procedures
      Public  :: Mpol_Ffactor,Mag_Ffactor
      Private :: Dipo_Ffactor, Disp_Ffactor, Quad_Ffactor,&
                 Hxap_Ffactor, Int_Slater_Bessel
      real(kind=dp), parameter :: pn  =  0.2695420113693928312
      real(kind=dp), parameter :: schw= -0.00014699
      logical           :: err_flipr=.false.
      character(len=80) :: err_flipr_mess=" "

    contains

    Function Schwinger_Amplitude(hn,pol,cellp,Atm,Grp,Left,UB) Result(Schwinger)
      real(kind=cp), dimension(3),            intent(in) :: hn,pol
      type(Crystal_Cell_Type),                intent(in) :: cellp
      type(Atom_List_Type),                   intent(in) :: Atm
      type(Space_Group_Type),                 intent(in) :: Grp
      logical,                      optional, intent(in) :: left
      real(kind=cp), dimension(3,3),optional, intent(in) :: UB
      complex(kind=cp)                                   :: Schwinger

      real(kind=cp),dimension(3) :: uvect, hc
      real(kind=cp) :: theta

      if(present(left)) then
        if(left) then
          uvect=(/0.0,0.0,1.0/)
        else
          uvect=(/0.0,0.0,-1.0/)
        end if
      else if(present(UB)) then
        hc=matmul(UB,hn)
        uvect=cross_product((/0.0,1.0,0.0/),hc)
        uvect=uvect/sqrt(dot_product(uvect,uvect))
      else
        Schwinger=0.0
        return
      end if

    End Function Schwinger_Amplitude

    Subroutine Calc_Electrostatic_StrFactor(hn,sn,Atm,Grp,fc)
       !---- Arguments ----!
       real(kind=cp),dimension(3),         intent(in) :: hn
       real(kind=cp),                      intent(in) :: sn !(sinTheta/Lambda)**2
       type(atom_list_type),               intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp
       complex,                            intent(out):: fc

       !---- Local Variables ----!
       integer                               :: i,j,k,m
       real(kind=cp)                         :: arg,anis,scosr,ssinr,b,s
       real(kind=cp)                         :: a1,a3,b1,b3,av,bv,f
       real(kind=cp),dimension(3)            :: h
       real(kind=cp),dimension(6)            :: beta
       real(kind=cp),dimension(Atm%natoms)   :: otr,oti,afpxn

       !--- Initialising local variables
       a1=0.0;  a3=0.0
       b1=0.0;  b3=0.0
       otr=0.0
       oti=0.0

       do i=1,Atm%natoms
          arg=0.0
          scosr=0.0
          ssinr=0.0
          do k=1,grp%NumOps
             h=matmul(hn,grp%Symop(k)%Rot)
             arg=tpi*(dot_product(h,Atm%atom(i)%x)+grp%Symop(k)%tr)
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

          b= atm%atom(i)%occ * exp(-atm%atom(i)%biso*s*s)
          otr(i) = afpxn(i)*b     ! (f0+Deltaf')*OCC*Tiso
          oti(i) =  afpp(i)*b     !     Deltaf" *OCC*Tiso
          a1 = a1 + otr(i)*scosr  ! F=A+iB: components of A  and B (ai,bi)
          b1 = b1 + oti(i)*scosr  ! a2,b2,a4,b4 are components for anisotropic form factors
          a3 = a3 + oti(i)*ssinr  ! they are not used here
          b3 = b3 + otr(i)*ssinr  ! For general case: av = a1-a2-a3-a4, bv = b1-b2+b3+b4

       end do ! Atoms

       av = a1-a3    !real part of the structure factor
       bv = b1+b3    !imaginary part of the structure factor
       fc=cmplx(av,bv)

    End Subroutine Calc_Electrostatic_StrFactor

    Subroutine Flipr(mode,iext,extc,polarp,polarm,hn,sn,cellp,Atm,Grp,flr,deriv,fn,fmag)
      !  Standard calculation of flipping ratios, including the complex
      !  form-factors and magnetic form factors. The calculation of special form
      !  factors is performed in subroutine form_factor. That of magnetic form
      !  factors is performed in subroutines mpol_ffactor and/or mag_ffactor
      character(len=1),                   intent(in) :: mode !Powder "P" or single Xtal "X"
      integer,                            intent(in) :: iext
      real(kind=cp), dimension(6),        intent(in) :: extc
      real(kind=cp), dimension(3),        intent(in) :: hn
      real(kind=cp),                      intent(in) :: sn,polarp, polarm
      type(Crystal_Cell_Type),            intent(in) :: cellp
      type(Atom_List_Type),               intent(in) :: Atm
      type(Space_Group_Type),             intent(in) :: Grp
      real(kind=cp),                      intent(out):: flr
      real(kind=cp),dimension(:),optional,intent(out):: deriv
      complex,                   optional,intent(out):: fn,fmag
      !-----------------------------------------------
      !   L o c a l   P a r a m e t e r s
      !-----------------------------------------------
      integer, parameter :: nat_pl=nat_p     ! Normally nat_pL=nat_p
      real(kind=cp), parameter :: EPSIL=1.E-06
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      integer :: i, j, k, n, ni, ii, irr, mni, i_sc
      real(kind=cp), dimension(nat_pl)   :: xi
      real(kind=cp), dimension(Atm%natoms)        :: otr, oti, frc, frs, fic, fis
      real(kind=cp), dimension(Atm%natoms)        :: motr, moti, mfrc, mfrs, mfic, mfis !magnetic part
      real(kind=cp), dimension(3)                 :: t
      real(kind=cp), dimension(Atm%Nvar)            :: mdxir,mdxii !magnetic part
      real(kind=cp), dimension(Atm%natoms,Atm%Nvar) :: drc,drs,dic,dis,dr,di
      real(kind=cp), dimension(3)        :: h
      real(kind=cp), dimension(3,3)      :: sm
      real(kind=cp) :: fr, fi, ffr, ffi, cosr, cosi, sinr, sini, scosr, scosi, &
                       ssinr, ssini, a1, a2, a3, a4, b1, b2, b3, b4, temp, scal
      real(kind=cp) :: av,bv, x, yy, arg, arg2, exparg, der, q
      real(kind=cp) :: mav, mbv          ! real(kind=cp) and imaginary part of the magnetic structure factor.
      real(kind=cp) :: mfr,mfi,mffr,mffi !magnetic form factors
      real(kind=cp) :: ma1, ma2, ma3, ma4, mb1, mb2, mb3, mb4
      real(kind=cp) :: mda, mdb, mda1, mda2, mda3, mda4, mdb1, mdb2, mdb3, mdb4
      real(kind=cp) :: mcosr, mcosi, msinr, msini, mscosr, mscosi, &
                       mssinr, mssini
      real(kind=cp) :: MAG,NUC                ! Magnetic and Nuclear squared structure factors
      real(kind=cp) :: AAlpha, Iplus, Iminus  ! Intermediate variables used in the FR calculus.
      real(kind=cp) :: DIpA, DIpB, DImA, DImB ! Intermediate variables used calculus of FR derivatives.
      Logical       :: magnetic
      ! real(kind=cp), dimension(3) :: hl  !! h refered to the local axes of the magnetic atom
      real(kind=cp), dimension(3,3) :: MatRot  !! Rotation matrix relating the cartesian reference system
                                               !! defined by the lattice to the cartesian reference system
                                               !! of the atoms in the asymetric unit.
      !--------------------------------------------------------------------------
      !   L o c a l   V a r i a b l e s  R e l a t e d  t o  E x t i n c t i o n
      !--------------------------------------------------------------------------
      real(kind=cp) :: DIpyp, DIpym, DIpypm,&  ! Intermediate variables used calculus of FR derivatives.
                       DImyp, DImym, DImypm    ! with respect to extinction parameters
      real(kind=cp) :: dflrmav, dflrmbv        !! Derivatives of flipping ratios with respect to mav and mbv
      real(kind=cp) :: dflryp, dflrym, dflrypm !! Derivatives of flipping ratios with respect to
                                               !! extinction coefficients
      real(kind=cp) :: yp, ym, ypm, ppp, ppm, mpp, mpm  !! Extinction correction + polar-corrections
      real(kind=cp), dimension(6) :: dyp, dym, dypm     !! Derivatives of yp, ym, and ypm extinction
                                                        !! coef. with respect to input parameters
      real(kind=cp), dimension(6) ::dymag               !! stores the derivatives of yp, ym, and ypm on mav and mbv
      !-----------------------------------------------
      !  SN: (sintheta/Lambda)**2
      n=Atm%natoms
      drc(1:n,1:nat_pl)=0.0
      drs(1:n,1:nat_pl)=0.0
      dic(1:n,1:nat_pl)=0.0
      dis(1:n,1:nat_pl)=0.0
      dr(1:n,1:nat_pl) =0.0
      di(1:n,1:nat_pl) =0.0
!      dxir(:)=0.0 ; dxii(:)=0.0
      mdxir(:)=0.0 ; mdxii(:)=0.0
      a1=0.0    ;   ma1=0.0
      a2=0.0    ;   ma2=0.0
      a3=0.0    ;   ma3=0.0
      a4=0.0    ;   ma4=0.0
      b1=0.0    ;   mb1=0.0
      b2=0.0    ;   mb2=0.0
      b3=0.0    ;   mb3=0.0
      b4=0.0    ;   mb4=0.0
      mav=0.0
      mbv=0.0


      !  LOOP OVER ATOMS N ATOMS OF PHASE IPH

      !---------------------------
      DO i=1,n           !loop over atoms
      !---------------------------
        ipiof=i+iof
        magnetic=.false.

        ffr= 0.0                  !
        ffi= 0.0                  !
        mfr= 0.0                  ! We start reseting all form factors
        mfi= 0.0                  !
        mffr= 0.0                 !
        mffi= 0.0                 !

        IF(ntyp(ipiof) == 'DIPO' .or. ntyp(ipiof) == 'QUAD' .or. ntyp(ipiof) == 'HXAP' .or. &
           ntyp(ipiof) == 'OCTU' .or. ntyp(ipiof) == 'DISP' .or. ntyp(ipiof) == 'MULT' ) THEN
                magnetic=.true.  !Magnetic atoms are given with special form factors
                mffr= 1.0        ! Isotropic magnetic form factor is just set to 1
                mffr= 0.0
        END IF

        xi(1:nat_pl)=xl(ipiof,1:nat_pl)

        temp=EXP(-xl(ipiof,4)*sn)   !exp{-Bi (sintheta/Lambda)^2}

        ni=ptr(ipiof,1,n_pat)

        IF(ni > 0) THEN                ! Nuclear contribution
            ffi=dfpp(ni,n_pat)         ! Imaginary Fermi length
            ffr=dfp(ni,n_pat)          ! Fermi Length
        ELSE                           ! Nuclear and Magnetic contribution
            ffi=dfpp(-ni,n_pat)        ! Imaginary Fermi length
            ffr=dfp(-ni,n_pat)         ! Fermi Length
            if(ntyp(ipiof) == 'MPOL') then
               mfr= 1.0                  ! Especial magnetic form factor is just set to 1
               mfi= 0.0
               mni=-ptr(ipiof,2,n_pat)
               !This is a magnetic atom with simple form factor (spin+orbital contribution): mu {<j0>+C2<j2>+C4<j4>+C6<j6>}
               !Nevertheless, magnetic is kept to be .false. in order to avoid the calculation of multipolar terms
               !in subroutine mag_ffactor called below if magnetic=.true. .
               CALL mpol_ffactor(n_pat,mni,ipiof,sn,mffr,mffi,mdxir,mdxii)
               if (.not.(maxs == 0 .OR. abs(icryg) >= 2) ) then
                  if (ABS(mffr) > epsil) then
                    dr(i,:)= mdxir/mffr
                  end if
                  if (ABS(mffi) > epsil) then
                    di(i,:)= mdxii/mffi
                  end if
               end if
            else
               MatRot=den_nlm(sp_pointers(ipiof))%mat
               mffr= 1.0    ! Isotropic magnetic form factor is just set to 1
               mffi= 0.0
               mni=-ptr(ipiof,2,n_pat)
            end if
        END IF

        ! Loop over symmetry operators
        ! Set indices for applying selected symmetry operators
        ! Identity is always applied and must be the first sym.op.

        scosr=0.0 ; mscosr=0.0        !
        scosi=0.0 ; mscosi=0.0        !  We just reset all the needed quantities.
        ssinr=0.0 ; mssinr=0.0        !
        ssini=0.0 ; mssini=0.0        !

        !+++++++++++++++++++++++++
        DO  irr=1,irl   !Loop over symmetry operators
        !+++++++++++++++++++++++++
          sm(:,:)=Spgr(iph)%Symop(irr)%Rot(:,:)
             t(:)=SpGr(iph)%Symop(irr)%tr(:)
          x=0.0
          DO  ii=1,3
            x=t(ii)*hn(ii)+x
            yy=0.0
            DO  j=1,3
              yy= hn(j)*sm(j,ii)+yy
            END DO
            h(ii)=yy
          END DO
          arg=x
          DO j=1,3
            arg=h(j)*xi(j)+arg
          END DO
          arg=tpi*arg
          arg2=h(1)*h(1)*xi(13)+h(2)*h(2)*xi(14)+ h(3)*h(3)*xi(15)+        &
               2.0*h(1)*h(2)*xi(16)+ 2.0*h(1)*h(3)*xi(17)+2.0*h(2)*h(3)*xi(18)
          exparg=EXP(-arg2)
            fr=1.0
            fi=0.0
          if(ni < 0) then  !special form factors
            !Nuclear contribution fr=Fermi length, fi=0.0
            !Magnetic contribution with special form factors (d-orbitals, multipoles, etc)
            if(magnetic) then

               h= unit_cart_hkl(hn,cellp(iph))     ! transform (h,k,l) to direct (notice we are using hn!!!)
                                                    ! cell in cartesian coordinates
               h= matmul(MatRot,h)                  ! Applies the rotation to local cartesian
                                                    ! coordinates
               h= matmul(cellp(iph)%Orth_Cr_cel,h)  ! Passes the result to direct
                                                    ! cell coordinates
               h= matmul(sm,h)                      ! Applies the current symetry operation
                                                    ! to the result
               h= matmul(cellp(iph)%Cr_Orth_cel,h)  ! Returns it back to cartesian coordinates

               if(maxs == 0 .OR. abs(icryg) >= 2 ) then
                 CALL mag_ffactor(n_pat,ipiof,sn,h,mfr,mfi)
               else
                 CALL mag_ffactor(n_pat,ipiof,sn,h,mfr,mfi,mdxir,mdxii)
               end if
            end if
          end if
          !Nuclear contribution if read_strf == 0
          if (read_strf == 0) then
             cosr=COS(arg)*exparg*fr   !fr*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
             cosi=COS(arg)*exparg*fi   !fi*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
             sinr=SIN(arg)*exparg*fr   !fr*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
             sini=SIN(arg)*exparg*fi   !fi*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})

             scosr=scosr+cosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)
             scosi=scosi+cosi          !FIC= SIG fi(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)

            if(icent == 1) then
             ssinr=ssinr+sinr          !FRS= SIG fr(j,s)sin{2pi(hT Rs rj+ts)}*Ta(s)
             ssini=ssini+sini          !FIS= SIG fi(j,s)sin{2pi(hT Rs rj+ts)}*Ta(s)
            end if
          end if

         !Magnetic contribution

          mcosr=COS(arg)*exparg*mfr   !fr*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
          mcosi=COS(arg)*exparg*mfi   !fi*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
          msinr=SIN(arg)*exparg*mfr   !fr*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
          msini=SIN(arg)*exparg*mfi   !fi*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})

          mscosr=mscosr+mcosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)
          mscosi=mscosi+mcosi          !FIC= SIG fi(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)

         IF(icent == 1) then
          mssinr=mssinr+msinr          !FRS= SIG fr(j,s)sin{2pi(hT Rs rj+ts)}*Ta(s)
          mssini=mssini+msini          !FIS= SIG fi(j,s)sin{2pi(hT Rs rj+ts)}*Ta(s)
         END IF

         ! Calculation of some quantities useful for derivatives
         IF(.not.(maxs == 0 .OR. abs(icryg) >= 2 )) then

         ! For the calculation of derivatives of MA1 and MB1

          IF(ni < 0 .AND. ABS(mfr) > epsil ) THEN
            DO j=19,nat_pl
              drc(i,j)= drc(i,j)+mcosr*mdxir(j)/mfr
            END DO
          END IF
          ! For the calculation of derivatives of MA4 and MB4
          IF(ni < 0 .AND. ABS(mfi) > epsil) THEN
            DO j=19,nat_pl
              dic(i,j)= dic(i,j)+mcosi*mdxii(j)/mfi
            END DO
          END IF

          IF(icent == 1) THEN

          ! For the calculation of derivatives of MA3 and MB3

            IF(ni < 0 .AND. ABS(mfr) > epsil) THEN
              DO j=19,nat_pl
                drs(i,j)= drs(i,j)+msinr*mdxir(j)/mfr
              END DO
            END IF
            ! For the calculation of derivatives of MA2 and MB2
            IF(ni < 0 .AND. ABS(mfi) > epsil ) THEN
              DO j=19,nat_pl
                dis(i,j)= dis(i,j)+msini*mdxii(j)/mfi
              END DO
            END IF

          END IF
         END IF     ! derivative calculations?

        !+++++++++++++++++++++++++
        END DO     !over symmetry operators  !  END LOOP SYMM.OP.
        !+++++++++++++++++++++++++

        if (read_strf==0) then
           frc(i)=scosr    !Components of geometrical struture factor of atom i
           fis(i)=ssini    !all components contribute to real(kind=cp) and imaginary parts
           frs(i)=ssinr    !in the most general case.
           fic(i)=scosi

           otr(i)=ffr*xl(ipiof,5)*temp     ! (f0+Deltaf')*OCC*Tiso
           oti(i)=ffi*xl(ipiof,5)*temp     !     Deltaf" *OCC*Tiso

           !-----CALCULATE A AND B OF F
           a1 = a1 + otr(i)*frc(i)    ! components of A  and B
           a4 = a4 + oti(i)*fic(i)    ! A(h) = a1 - a2 - a3 - a4
           b1 = b1 + oti(i)*frc(i)    ! B(h) = b1 - b2 + b3 + b4
           b4 = b4 + otr(i)*fic(i)    !

           if(icent == 1) then
             a2 = a2 + otr(i)*fis(i)
             a3 = a3 + oti(i)*frs(i)
             b2 = b2 + oti(i)*fis(i)
             b3 = b3 + otr(i)*frs(i)
           end if
        end if

        mfrc(i)=mscosr    !Components of geometrical struture factor of atom i
        mfis(i)=mssini    !all components contribute to real(kind=cp) and imaginary parts
        mfrs(i)=mssinr    !in the most general case.
        mfic(i)=mscosi

        motr(i)=mffr*xl(ipiof,5)*temp     ! mu*fr(mag)*OCC*Tiso
        moti(i)=mffi*xl(ipiof,5)*temp     ! mu*fi(mag) *OCC*Tiso


        !-----CALCULATE Am AND Bm OF Fm
        ma1 = ma1 + motr(i)*mfrc(i)    ! components of A  and B
        ma4 = ma4 + moti(i)*mfic(i)    ! Am(h) = ma1 - ma2 - ma3 - ma4
        mb1 = mb1 + moti(i)*mfrc(i)    ! Bm(h) = mb1 - mb2 + mb3 + mb4
        mb4 = mb4 + motr(i)*mfic(i)    !

        IF(icent == 1) THEN
          ma2 = ma2 + motr(i)*mfis(i)
          ma3 = ma3 + moti(i)*mfrs(i)
          mb2 = mb2 + moti(i)*mfis(i)
          mb3 = mb3 + motr(i)*mfrs(i)
        END IF

      !---------------------------
      END DO  !over atoms
      !---------------------------

        av=a1-a2-a3-a4        !real part of the nuclear structure factor
        bv=b1-b2+b3+b4        !imaginary part of the nuclear structure factor
        mav=ma1-ma2-ma3-ma4   !real part of the magnetic structure factor
        mbv=mb1-mb2+mb3+mb4   !imaginary part of the magnetic structure factor

        ! Correction with scale factor

        select case(abs(icod(nn,n_pat)))
          case (1)
            i_sc= 1
            scal= par(iph,i_sc,n_pat)
          case (2)
            i_sc= 2
            scal= par(iph,i_sc,n_pat)
          case (3)
            i_sc= 3
            scal= par(iph,i_sc,n_pat)
          case (4)
            i_sc= 4
            scal= par(iph,i_sc,n_pat)
          case (5)
            i_sc= 5
            scal= par(iph,i_sc,n_pat)
          case (6)
            i_sc= 12
            scal= par(iph,i_sc,n_pat)
          case default
            i_sc= 0
            scal= 0.0
            call Error_Message('  Warning!!, Integer identifying scale factor larger than 6. A zero is assumed.')
        end select

        mav= scal*mav
        mbv= scal*mbv

        q= gam(nn,n_pat)
        AAlpha= q * ( av*mav + bv*mbv )
        MAG= mav * mav + mbv * mbv
        NUC=  av *  av +  bv *  bv

        ggi(nn,n_pat) =icent*av     !the different components of the nuclear and magnetic
        cci(nn,n_pat) =icent*bv     !structure factors
        asfa(nn,n_pat)=icent*mav    !Store the magnetic structure factor
        bsfa(nn,n_pat)=icent*mbv    !The value of the structure factor corresponds to a primitive cell

        ! Extinction correction

        if (iext /= 0) then
           call Correct_FlippingRatios(iext,Lambda,q,extc,sn,hn,Av,Bv,mav,mbv,yp,ym,ypm,dyp,dym,dypm,dymag)

           ppp=((1.0+polarp)*yp+(1.0-polarp)*ym)*0.5 !+Pp
           ppm=((1.0+polarp)*yp-(1.0-polarp)*ym)*0.5 !+Pm
           mpp=((1.0-polarm)*yp+(1.0+polarm)*ym)*0.5 !-Pp
           mpm=((1.0-polarm)*yp-(1.0+polarm)*ym)*0.5 !-Pm

           Iplus=  (NUC+q*q*MAG)*ppp  + 2.0*q*(av*mav+ bv*mbv)*ppm  +(1.0-q)*q*MAG*ypm
           Iminus= (NUC+q*q*MAG)*mpp  + 2.0*q*(av*mav+ bv*mbv)*mpm  +(1.0-q)*q*MAG*ypm

           if(.not. (maxs == 0 .OR. abs(icryg) >= 2)) then
             DIpA=  2.0*q*q*mav*ppp  + 2.0*q* av*ppm + 2.0*(1.0-q)*q*mav*ypm
             DIpB=  2.0*q*q*mbv*ppp  + 2.0*q* bv*ppm + 2.0*(1.0-q)*q*mbv*ypm
             DImA=  2.0*q*q*mav*mpp  + 2.0*q* av*mpm + 2.0*(1.0-q)*q*mav*ypm
             DImB=  2.0*q*q*mbv*mpp  + 2.0*q* bv*mpm + 2.0*(1.0-q)*q*mbv*ypm

             DIpyp=  (NUC+q*q*MAG+2.0*q*(av*mav+ bv*mbv))*(1.0+polarp)/2.0
             DIpym=  (NUC+q*q*MAG-2.0*q*(av*mav+ bv*mbv))*(1.0-polarp)/2.0
             DIpypm= (1.0-q)*q*MAG

             DImyp=  (NUC+q*q*MAG+2.0*q*(av*mav+ bv*mbv))*(1.0-polarm)/2.0
             DImym=  (NUC+q*q*MAG-2.0*q*(av*mav+ bv*mbv))*(1.0+polarm)/2.0
             DImypm= (1.0-q)*q*MAG

             DIpA= DIpA+ DIpyp*dymag(1)+ DIpym*dymag(3)+ DIpypm*dymag(5)
             DIpB= DIpB+ DIpyp*dymag(2)+ DIpym*dymag(4)+ DIpypm*dymag(6)

             DImA= DImA+ DImyp*dymag(1)+ DImym*dymag(3)+ DImypm*dymag(5)
             DImB= DImB+ DImyp*dymag(2)+ DImym*dymag(4)+ DImypm*dymag(6)
           end if
        else

           Iplus=  NUC + 2.0*polarp* AAlpha+ q*MAG
           Iminus= NUC - 2.0*polarm* AAlpha+ q*MAG

           ! Derivative with respect to the real and imaginary magnetic part mav and mbv
           if(.not. (maxs == 0 .OR. abs(icryg) >= 2)) then
             DIpA= 2.0*q*(mav +polarp*av)
             DIpB= 2.0*q*(mbv +polarp*bv)
             DImA= 2.0*q*(mav -polarm*av)
             DImB= 2.0*q*(mbv -polarm*bv)
             ! No derivatives with respect to extinction here.
             DIpyp= 0.0
             DIpym= 0.0
             DIpypm= 0.0

             DImyp= 0.0
             DImym= 0.0
             DImypm= 0.0
           end if
        end if

        ! The following lines are valid if input data are flipping ratios
        ! Flipping ratio is calculated and stored in flr
        flr=Iplus/Iminus
        ! Derivatives of the flipping ratio with respect to mav and mbv
        if(.not. (maxs == 0 .OR. abs(icryg) >= 2)) then
          dflrmav= (DIpA-flr*DImA)/Iminus
          dflrmbv= (DIpB-flr*DImB)/Iminus
          dflryp= (DIpyp-flr*DImyp)/Iminus
          dflrym= (DIpym-flr*DImym)/Iminus
          dflrypm= (DIpypm-flr*DImypm)/Iminus
        end if


        ! The following lines are valid if input data are asymmetries
        ! Asymmetry is calculated and stored in flr
        !        flr= (Iplus-Iminus)/(Iplus+Iminus)
        !        dflrmav= 2.0*(Iminus*DIpA-Iplus*DImA)/(Iplus+Iminus)**2
        !        dflrmbv= 2.0*(Iminus*DIpB-Iplus*DImB)/(Iplus+Iminus)**2
        !        dflryp= 2.0*(Iminus*DIpyp-Iplus*DImyp)/(Iplus+Iminus)**2
        !        dflrym= 2.0*(Iminus*DIpym-Iplus*DImym)/(Iplus+Iminus)**2
        !        dflrypm= 2.0*(Iminus*DIpypm-Iplus*DImypm)/(Iplus+Iminus)**2

       !-----Calculate the phase of the magnetic structure factor in radians
      if(jfou(n_pat) /= 0) then
        phasen(nn,n_pat)=0.0
        if(icent == 1) then
          if (abs(mag) < 0.00001) then
            arg=0.0
          else
            arg= mbv/sqrt(mag)
          end if
          if(abs(arg) > 1.0) arg=sign(1.0_cp,arg)
          if (abs(mag) >= 0.000001) phasen(nn,n_pat)=asin(arg)
          if (mav < 0.0) phasen(nn,n_pat)=pi-phasen(nn,n_pat)
          if (phasen(nn,n_pat) < 0.0) phasen(nn,n_pat)=phasen(nn,n_pat)+tpi
        else
          if(mav < 0.0) phasen(nn,n_pat) = pi
        end if
      end if

      IF(maxs == 0 .OR. abs(icryg) >= 2) RETURN   !
      !----------------------------------------------------------------------------
      !-----CALCULATE DERIVATIVES of Flipping Ratios
      !     Only parameters that affect the magnetic form factor can be refined.
      !     So, only for such parameters, FR derivatives are calculated.
      !     All the other parameter's derivatives are set to zero.
      !----------------------------------------------------------------------------

      if (i_sc /= 0) then                    ! Scale factors
        k=lpar(iph,i_sc,n_pat)
        if (k /= 0) then
           if ( scal > 0.0001) then
              der= mav/scal*dflrmav+ mbv/scal*dflrmbv
              deriv(k)= sign(1.0_cp,apar(iph,i_sc,n_pat))*der+deriv(k)
           end if
        end if
      end if

      do j=13,18                 ! Extinction correction coefficients
        k=lpar(iph,j,n_pat)
        if(k /= 0) then
          der= dflryp*dyp(j) +dflrym*dym(j) +dflrypm*dypm(j)
          deriv(k) = sign(1.0_cp,apar(iph,j,n_pat))*der+deriv(k)
        end if
      end do

      DO i=1,n           !Loop over atoms
        ipiof=i+iof
        !--------------------Form-factor parameters
        IF(n_t(ipiof) >= 4) THEN

          DO j=19,nat_pl
           k=lp(ipiof,j)
           IF(k /= 0) THEN
             mda1= motr(i)*dr(i,j)*mfrc(i)+motr(i)*drc(i,j)
             mda2= motr(i)*dr(i,j)*mfis(i)+motr(i)*dis(i,j)
             mda3= moti(i)*di(i,j)*mfrs(i)+moti(i)*drs(i,j)
             mda4= moti(i)*di(i,j)*mfic(i)+moti(i)*dic(i,j)

             mdb1= moti(i)*di(i,j)*mfrc(i)+moti(i)*drc(i,j)
             mdb2= moti(i)*di(i,j)*mfis(i)+moti(i)*dis(i,j)
             mdb3= motr(i)*dr(i,j)*mfrs(i)+motr(i)*drs(i,j)
             mdb4= motr(i)*dr(i,j)*mfic(i)+motr(i)*dic(i,j)

             mda= mda1 -mda2 -mda3 -mda4
             mdb= mdb1 -mdb2 +mdb3 +mdb4

             der= dflrmav*mda +dflrmbv*mdb
             deriv(k) = sign(1.0_cp,a(ipiof,j))*der+deriv(k)

           END IF
          END DO
        !-------------------
        END IF
        !-------------------

      END DO !I=1,N over atoms

      RETURN
    End Subroutine Flipr

    Subroutine Mag_Ffactor(n_pat,ipiof,sn,h,mffr,mffi,mdxir,mdxii)
       !-----------------------------------------------
       !   D u m m y   A r g u m e n t s
       !-----------------------------------------------
       integer,                              intent(in)   :: n_pat
       integer,                              intent(in)   :: ipiof
       real(kind=cp),                        intent(in)   :: sn
       real(kind=cp), dimension(3),          intent(in)   :: h
       real(kind=cp),                        intent(out)  :: mffr
       real(kind=cp),                        intent(out)  :: mffi
       real(kind=cp), dimension(:), optional,intent(out)  :: mdxir
       real(kind=cp), dimension(:), optional,intent(out)  :: mdxii
       !-----------------------------------------------
       !   L o c a l   V a r i a b l e s
       !-----------------------------------------------
       logical :: rdx
       integer :: i
       rdx= .false.
          if(present(mdxir)) then
             mdxir(:)=0.0
          end if
          if(present(mdxii)) then
             mdxii(:)=0.0
          end if
       Select Case (ntyp(ipiof))
           case ('DIPO')
                if (present(mdxir).and. present(mdxii)) then
                   call dipo_ffactor(n_pat,ipiof,sn,h,mffr,mffi,mdxir,mdxii)
                else
                   call dipo_ffactor(n_pat,ipiof,sn,h,mffr,mffi)
                end if
                rdx= .true.
           case ('QUAD')
                if (present(mdxir).and. present(mdxii)) then
                   call quad_ffactor(n_pat,ipiof,sn,h,mffr,mffi,mdxir,mdxii)
                else
                   call quad_ffactor(n_pat,ipiof,sn,h,mffr,mffi)
                end if
                rdx= .true.
           case ('HXAP')
                if (present(mdxir).and. present(mdxii)) then
                   call hxap_ffactor(n_pat,ipiof,sn,h,mffr,mffi,mdxir,mdxii)
                else
                   call hxap_ffactor(n_pat,ipiof,sn,h,mffr,mffi)
                end if
                rdx= .true.
           case ('DISP')
                if (present(mdxir).and. present(mdxii)) then
                   call disp_ffactor(ipiof,sn,h,mffr,mffi,mdxir,mdxii)
                else
                   call disp_ffactor(ipiof,sn,h,mffr,mffi)
                end if
                rdx= .true.
           case ('MULT')
                if (present(mdxir).and. present(mdxii)) then
                   call mult_ffactor(ipiof,sn,h,mffr,mffi,mdxir,mdxii)
                else
                   call mult_ffactor(ipiof,sn,h,mffr,mffi)
                end if
                rdx= .true.
        End Select

        if (rdx) then
           do i=19, 32
              if (lp(ipiof,i) == 0) then
                  if (present(mdxir)) then
                     mdxir (i) = 0.0
                  end if
                  if (present(mdxii)) then
                     mdxii (i) = 0.0
                  end if
              end if
           end do
        end if
        mffr= pn* mffr
        mffi= pn* mffi
        if (present(mdxir)) then
           mdxir= pn*mdxir
        end if
        if (present(mdxii)) then
           mdxii= pn*mdxii
        end if
        return
    End Subroutine Mag_Ffactor

    SubroutiNe Mpol_Ffactor(n_pat,ni,ipiof,sn,mffr,mffi,mdxir,mdxii)
       integer,                              intent(in) :: n_pat
       integer,                              intent(in) :: ni
       integer,                              intent(in) :: ipiof
       real(kind=cp),                        intent(in) :: sn
       real(kind=cp),                        intent(out):: mffr
       real(kind=cp),                        intent(out):: mffi
       real(kind=cp), dimension(:), optional,intent(out):: mdxir
       real(kind=cp), dimension(:), optional,intent(out):: mdxii
       ! Local variables
       integer :: i
       real(kind=cp)    :: j0,j2,j4,j6
          j0=ac(7,ni,n_pat)
          j2=ac(14,ni,n_pat)
          j4=ac(21,ni,n_pat)
          j6=ac(28,ni,n_pat)
          DO i=1,5,2
            j0=j0+ac(i,ni,n_pat)*EXP(-ac(i+1,ni,n_pat)*sn)      !<j0>      value for Q=H+k
            j2=j2+ac(i+7,ni,n_pat)*EXP(-ac(i+8,ni,n_pat)*sn)    !<j2>/sn value for Q=H+k
            j4=j4+ac(i+14,ni,n_pat)*EXP(-ac(i+15,ni,n_pat)*sn)  !<j4>/sn value for Q=H+k
            j6=j6+ac(i+21,ni,n_pat)*EXP(-ac(i+22,ni,n_pat)*sn)  !<j6>/sn value for Q=H+k
          END DO
            j2=j2*sn  !<j2> value for Q=H+k
            j4=j4*sn  !<j4> value for Q=H+k
            j6=j6*sn  !<j6> value for Q=H+k
            mffr= pn*( xl(ipiof,19)*j0 + xl(ipiof,20)*j2 + xl(ipiof,21)*j4 + xl(ipiof,22)*j6)
            mffi=0.0
          if(present(mdxir)) then
            mdxir(:) =0.0
            IF(lp(ipiof,19) /= 0) mdxir(19) = pn*j0
            IF(lp(ipiof,20) /= 0) mdxir(20) = pn*j2
            IF(lp(ipiof,21) /= 0) mdxir(21) = pn*j4
            IF(lp(ipiof,22) /= 0) mdxir(22) = pn*j6
          end if
          if(present(mdxii)) then
               mdxii(:) =0.0
          end if
       return
    End Subroutine Mpol_Ffactor

    Subroutine Dipo_ffactor(n_pat,ipiof,sn,huc,mffr,mffi,mdxir,mdxii)
       integer, intent(in)                                :: n_pat
       integer, intent(in)                                :: ipiof
       real(kind=cp), intent(in)                          :: sn
       real(kind=cp), intent(in) , dimension(3)           :: huc
       real(kind=cp), intent(out)                         :: mffr
       real(kind=cp), intent(out)                         :: mffi
       real(kind=cp), optional,intent(out), dimension(:)  :: mdxir
       real(kind=cp), optional,intent(out), dimension(:)  :: mdxii
       !-----------------------------------------------
       !   L o c a l   P a r a m e t e r s
       !-----------------------------------------------

       real(kind=cp), parameter :: D2R= 0.01745329251994
       integer, parameter :: ZPAR= 25
       integer, parameter :: ORBPAR= 20
       integer, parameter :: MAGPAR= 19

       !
       !   The following lines of code have been automatically generated
       !

       real(kind=cp), parameter ::                                                              &
                         AA01= 0.2185097040, AA02= -0.2185097040, AA03= -0.1261566430,  &
                         AA04= 0.2523132860, AA05=  0.2820947770
       !
       !   The preceding lines of code have been automatically generated
       !
       !-----------------------------------------------------------------------------------------------------
       !   L o c a l   V a r i a b l e s :   R a d i a l   P a r t   o f   t h e   W a v e   F u n c t i o n
       !-----------------------------------------------------------------------------------------------------
       real(kind=cp) :: Zeff          !   exponent and power in the Slater-Type function if Rad_Model >= 2
       integer :: Rad_Model  !   R(r)~ r**n * exp(-chi*r)         chi= Zeff/n     (n=Rad_Model)
       !-----------------------------------------------------------------------------------------------------
       !   L o c a l   V a r i a b l e s :   A n g u l a r   P a r t   o f   t h e   W a v e   F u n c t i o n
       !-----------------------------------------------------------------------------------------------------
       real(kind=cp), dimension(0:0)      :: Devpl0   ! Term with l=0 comming from p-terms
       real(kind=cp), dimension(0:0,1:5)  :: DDevpl0  ! Derivatives of the preceding term with respect
                                                      ! to normalized input parameters
       real(kind=cp), dimension(-2:2)     :: Devpl2   ! Term with l=2 comming from p-terms
       real(kind=cp), dimension(-2:2,1:5) :: DDevpl2  ! Derivatives of the preceding term with respect
                                                      ! to normalized input parameters
       real(kind=cp) :: a, b, c, ph1, ph2
       real(kind=cp) :: norm, mag
       real(kind=cp), dimension (3,3)     :: normcorr  !   Term for the correction to the derivatives
       real(kind=cp), dimension (3)       :: deraux
       !
       !
       !      \phi = R_p(r) \left\{ a*y_1^{1 +}+b \exp\{i ph1\}*y_1^{1 -} +c \exp\{i ph2\}*y_1^{0 +} \right\}
       !
       !      \rho= ?\phi?^2= R(r)^2*\sum_l={0,2} \sum_{m=-l}^l DevdlX(m) *y_l^{abs(m),sign(1,m)}(\hat{r})   (X=0,2)
       !
       !      fm(\vec{h})= 4\pi \sum_l={0,2} (i)**l*\int_{0}^{\infty} r**2*R(r)^2*j_l(2\pi h r)*dr *
       !                        \sum_{m=-l}^l DevdlX(m) *y_l^{abs(m),sign(1,m)}(\hat{h})                     (X=0,2)
       !
       !-----------------------------------------------
       !   L o c a l   V a r i a b l e s
       !-----------------------------------------------
       real(kind=cp) ::  YY, summ
       real(kind=cp), dimension(5) ::  summd
       real(kind=cp), dimension (0:1) :: IRJ, DIRJ
       real(kind=cp):: hm    !modulus of q=2*pi*|h|
       integer      :: l, m, p

       hm=2.0*tpi*sqrt(sn)
       Rad_Model = den_nlm(sp_pointers(ipiof))%nl(2)

       DDevpl0= 0.0
       DDevpl2= 0.0
       xl(ipiof, ORBPAR+ 2:ORBPAR+ 4:2)= MOD(xl(ipiof, ORBPAR+ 2:ORBPAR+ 4:2), 360.0_cp)

       Zeff= xl(ipiof, ZPAR)
       a  = xl(ipiof, ORBPAR)
       b  = xl(ipiof, ORBPAR+1)
       ph1= xl(ipiof, ORBPAR+2) * D2R
       c  = xl(ipiof, ORBPAR+3)
       ph2= xl(ipiof, ORBPAR+4) * D2R
       mag= xl(ipiof, MAGPAR)
       norm= sqrt(a**2+b**2+c**2)

       ! If norm = 0 We have nothing to do !!
       if (norm <= 1.e-5) then
          call Error_Message("Zero norm of an orbital, nothing is calculated", 6)
         mffr= 0.0
          mffi= 0.0
          if (present(mdxir)) then
             mdxir= 0.0
          end if
          if (present(mdxii)) then
             mdxii= 0.0
          end if
          return
       end if
       a= a/norm
       b= b/norm
       c= c/norm

       ! If not, we store the correction to derivatives due to normalization

          normcorr= RESHAPE ((/ 1.0-a**2, -a*b, -a*c, &
                                -b*a, 1.0-b**2, -b*c, &
                                -c*a, -c*b, 1.0-c**2 /),(/3,3/))/norm
      !
      !   The following lines of code have been automatically generated
      !
      Devpl2( 2)= +AA01 *a*a +AA02 *b*b
      Devpl2( 1)= +AA01 *2.0*a*c*COS(   -ph2)
      Devpl2( 0)= +AA03 *a*a +AA03 *b*b +AA04 *c*c
      Devpl2(-1)= +AA01 *2.0*b*c*COS(ph1-ph2)
      Devpl2(-2)= +AA01 *2.0*a*b*COS(   -ph1)
      Devpl0( 0)= +AA05 *a*a +AA05 *b*b +AA05 *c*c

      DDevpl2( 2, 1)= +AA01 *2.0*a
      DDevpl2( 2, 2)= +AA02 *2.0*b
      DDevpl2( 1, 1)= +AA01 *2.0*c*COS(   -ph2)
      DDevpl2( 1, 4)= +AA01 *2.0*a*COS(   -ph2)
      DDevpl2( 1, 5)= +AA01 *2.0*a*c*SIN(   -ph2)*D2R
      DDevpl2( 0, 1)= +AA03 *2.0*a
      DDevpl2( 0, 2)= +AA03 *2.0*b
      DDevpl2( 0, 4)= +AA04 *2.0*c
      DDevpl2(-1, 2)= +AA01 *2.0*c*COS(ph1-ph2)
      DDevpl2(-1, 3)= -AA01 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevpl2(-1, 4)= +AA01 *2.0*b*COS(ph1-ph2)
      DDevpl2(-1, 5)= +AA01 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevpl2(-2, 1)= +AA01 *2.0*b*COS(   -ph1)
      DDevpl2(-2, 2)= +AA01 *2.0*a*COS(   -ph1)
      DDevpl2(-2, 3)= +AA01 *2.0*a*b*SIN(   -ph1)*D2R

      DDevpl0( 0, 1)= +AA05 *2.0*a
      DDevpl0( 0, 2)= +AA05 *2.0*b
      DDevpl0( 0, 4)= +AA05 *2.0*c

       !
       !   The preceding lines of code have been automatically generated
       !
       if ( 0 <= Rad_Model .and. Rad_Model <= 1 ) then
          call IRJ_ITB(n_pat,-ptr(ipiof,2,n_pat),sn,IRJ,DIRJ)
       else if (Rad_Model < 0) then
          call IRJ_Slater(-Rad_Model, 2, Zeff, hm, IRJ, DIRJ)
       else
          call IRJ_Hydrog(Rad_Model, 1, 2, Zeff, hm, IRJ, DIRJ)
       end if
       mffr=0.0
       mffi=0.0

       ! Calculation for l=2
       l=2

       summ= 0.0
       summd= 0.0

       do m=-l,l
         p= sign(1,m)
         YY= Real_Spher_Harm_Ucvec(l,ABS(m),p,huc)
         summ= summ+ YY*Devpl2(m)
         summd= summd+ YY*DDevpl2(m,:)
       end do

       mffr= mffr- IRJ(l/2)*summ                                   !   The minus sign comes from the term i**l with l=2.

       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+4)= mdxir(ORBPAR:ORBPAR+4)-IRJ(l/2)*summd  !   The minus sign comes from the term i**l with l=2.
          mdxir(ZPAR)= mdxir(ZPAR)- DIRJ(l/2)*summ                     !   The minus sign comes from the term i**l with l=2.
       end if

       ! Calculation for l=0
       l=0

       YY= 1.0/SQRT(4.0*PI)
       summ = YY*Devpl0(0)
       mffr= mffr+ IRJ(l/2)*summ
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+4)= mdxir(ORBPAR:ORBPAR+4)+ DIRJ(l/2)*YY*DDevpl0(0,1:5)
          mdxir(ZPAR)= mdxir(ZPAR) +DIRJ(l/2)*summ
       end if

        !   At this point the true mffr is 4*PI*MAG times the calculated mffr,
        ! and the derivatives are calculated with respect to normalized parameters but no with respect to
        ! input parameters
        !   We now correct the derivatives for normalization, scale (MAG), and 4*PI factor. The derivative
        ! with respect to the scale is also calculated.

       if (PRESENT(mdxir)) then
          mdxir(MAGPAR)= mffr*4.0*PI
       end if
       mffr= mffr*4.0*PI*MAG
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+4)= mdxir(ORBPAR:ORBPAR+4) *4.0*PI*MAG
          mdxir(ZPAR)= mdxir(ZPAR) *4.0*PI*MAG
          deraux= matmul(normcorr,(/mdxir(ORBPAR), mdxir(ORBPAR+1:ORBPAR+3:2)/))
          mdxir(ORBPAR)= deraux(1)
          mdxir(ORBPAR+1:ORBPAR+3:2)= deraux(2:3)
       end if

       return

    End Subroutine dipo_ffactor

    Subroutine Disp_ffactor(ipiof,sn,huc,mffr,mffi,mdxir,mdxii)
       integer,       intent(in)                     :: ipiof
       real(kind=cp), intent(in)                     :: sn
       real(kind=cp), intent(in) , dimension(3)      :: huc
       real(kind=cp), intent(out)                    :: mffr
       real(kind=cp), intent(out)                    :: mffi
       real(kind=cp), optional,intent(out), dimension(:)  :: mdxir
       real(kind=cp), optional,intent(out), dimension(:)  :: mdxii
       !-----------------------------------------------
       !   L o c a l   P a r a m e t e r s
       !-----------------------------------------------

       real(kind=cp), parameter :: D2R= 0.01745329251994
       ! real(kind=cp), parameter :: a0= 0.529177249
       integer, parameter :: ZPAR= 27
       integer, parameter :: ORBPAR= 20
       integer, parameter :: MAGPAR= 19

       !
       !   The following lines of code have been automatically generated
       !

      real(kind=cp), parameter ::                                                              &
                         AA01= 0.2185097040, AA02= -0.2185097040, AA03= -0.1261566430,  &
                         AA04= 0.2523132860, AA05=  0.2820947770
       !
       !   The preceding lines of code have been automatically generated
       !
       !-----------------------------------------------------------------------------------------------------
       !   L o c a l   V a r i a b l e s :   R a d i a l   P a r t   o f   t h e   W a v e   F u n c t i o n
       !-----------------------------------------------------------------------------------------------------
       real(kind=cp) :: Zeffs, Zeffp          !   exponent and power in the Slater-Type function
       integer :: ns, np             !   R_s(r)~ r**ns * exp(-chis*r)
       integer :: Rad_Model          !   R_p(r)~ r**np * exp(-chip*r)
       !-----------------------------------------------------------------------------------------------------
       !   L o c a l   V a r i a b l e s :   A n g u l a r   P a r t   o f   t h e   W a v e   F u n c t i o n
       !-----------------------------------------------------------------------------------------------------
       real(kind=cp), dimension(0:0)       :: Devpl0    !   Term with l=0 comming from p-terms
       real(kind=cp), dimension (0:0,1:7)  :: DDevpl0   !  Derivatives of the preceding term with respect
                                                        !  to normalized input parameters
       real(kind=cp), dimension (-2:2)     :: Devpl2        !   Term with l=2 comming from p-terms
       real(kind=cp), dimension (-2:2,1:7) :: DDevpl2   !  Derivatives of the preceding term with respect
                                                        !  to normalized input parameters
       real(kind=cp), dimension(0:0)       :: Devsl0          !   Term with l=0 comming from s-terms
       real(kind=cp), dimension (0:0,1:7)  :: DDevsl0   !  Derivatives of the preceding term with respect
                                                        !  to normalized input parameters
       real(kind=cp), dimension(-1:1)      :: Devspl1   !   Term with l=1 comming from sp-terms
       real(kind=cp), dimension (-1:1,1:7) :: DDevspl1  !  Derivatives of the preceding term with respect
                                                        !  to normalized input parameters
       real(kind=cp) :: a, b, c, d, ph1, ph2, ph3
       real(kind=cp) :: norm, mag
       real(kind=cp), dimension (4,4) :: normcorr  !   Term with the correction to the derivatives
       real(kind=cp), dimension (4)   :: deraux    !   due to normallization
       !
       !
       !      phi = R_p(r)*\left\{ a*y_2^{1 +}+b \exp\{i ph1\}*y_2^{1 -} +c \exp\{i ph2\}*y_2^{0 +} \right\} +
       !            R_s(r)*\left\{ d \exp\{i ph3\}*y_0^{0 +} \right\}
       !
       !      \rho= ?\phi?^2= R_s(r)^2*Devsl0(0)*y_l^{0 +} + R_s(r)R_p \sum_{m=-1}^1 Devspl1 y_1^{abs(m),sign(1,m)}
       !                     +R_p(r)^2 \sum_l={0,2} \sum_{m=-l}^l Devdl?(m) *y_l^{abs(m),sign(1,m)}
       !
       !
       !-----------------------------------------------
       !   L o c a l   V a r i a b l e s
       !-----------------------------------------------
       real(kind=cp) ::  YY, summ
       real(kind=cp), dimension(7) ::  summd
       real(kind=cp), dimension(0:1) :: IRJ, DIRJ
       real(kind=cp):: hm     ! modulus of q=2*pi*|h|
       integer :: m, p, l

       Rad_Model = den_nlm(sp_pointers(ipiof))%nl(3)
       hm= 2.0*tpi*sqrt(sn)
       np= abs(den_nlm(sp_pointers(ipiof))%nl(1))
       ns= abs(den_nlm(sp_pointers(ipiof))%nl(2))
       DDevsl0= 0.0
       DDevpl0= 0.0
       DDevpl2= 0.0
       DDevspl1= 0.0
       Zeffp= xl(ipiof, ZPAR)
       Zeffs= xl(ipiof, ZPAR+1)
       xl(ipiof, ORBPAR+ 2:ORBPAR+ 6:2)= MOD(xl(ipiof, ORBPAR+ 2:ORBPAR+ 6:2), 360.0_cp)

       a  = xl(ipiof, ORBPAR)
       b  = xl(ipiof, ORBPAR+1)
       ph1= xl(ipiof, ORBPAR+2) * D2R
       c  = xl(ipiof, ORBPAR+3)
       ph2= xl(ipiof, ORBPAR+4) * D2R
       d  = xl(ipiof, ORBPAR+5)
       ph3= xl(ipiof, ORBPAR+6) * D2R
       mag= xl(ipiof, MAGPAR)
       norm= sqrt(a**2+b**2+c**2+d**2)

       ! If norm = 0 We have nothing to do !!
       if (norm <= 1.e-5) then
          call Error_Message("Zero norm of an orbital, nothing is calculated", 6)
          mffr= 0.0
          mffi= 0.0
          if (PRESENT(mdxir)) then
             mdxir= 0.0
          end if
          if (PRESENT(mdxii)) then
             mdxii= 0.0
          end if
          return
       end if
       a= a/norm
       b= b/norm
       c= c/norm
       d= d/norm

       ! If not, we store the correction to derivatives due to normalization

       normcorr= RESHAPE ((/ 1.0-a**2, -a*b, -a*c, -a*d, &
                             -b*a, 1.0-b**2, -b*c, -b*d, &
                             -c*a, -c*b, 1.0-c**2, -c*d, &
                             -d*a, -d*b, -d*c, 1.0-d**2 /),(/4,4/))/norm
       !
       !   The following lines of code have been automatically generated
       !
      Devpl2( 2) = +AA01 *a*a +AA02 *b*b
      Devpl2( 1) = +AA01 *2.0*a*c*COS(   -ph2)
      Devpl2( 0) = +AA03 *a*a +AA03 *b*b +AA04 *c*c
      Devpl2(-1) = +AA01 *2.0*b*c*COS(ph1-ph2)
      Devpl2(-2) = +AA01 *2.0*a*b*COS(   -ph1)
      Devpl0( 0) = +AA05 *a*a +AA05 *b*b +AA05 *c*c
      Devspl1( 1)= +AA05 *2.0*a*d*COS(   -ph3)
      Devspl1( 0)= +AA05 *2.0*c*d*COS(ph2-ph3)
      Devspl1(-1)= +AA05 *2.0*b*d*COS(ph1-ph3)
      Devsl0( 0) = +AA05 *d*d

      DDevpl2( 2, 1)= +AA01 *2.0*a
      DDevpl2( 2, 2)= +AA02 *2.0*b
      DDevpl2( 1, 1)= +AA01 *2.0*c*COS(   -ph2)
      DDevpl2( 1, 4)= +AA01 *2.0*a*COS(   -ph2)
      DDevpl2( 1, 5)= +AA01 *2.0*a*c*SIN( -ph2)*D2R
      DDevpl2( 0, 1)= +AA03 *2.0*a
      DDevpl2( 0, 2)= +AA03 *2.0*b
      DDevpl2( 0, 4)= +AA04 *2.0*c
      DDevpl2(-1, 2)= +AA01 *2.0*c*COS(ph1-ph2)
      DDevpl2(-1, 3)= -AA01 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevpl2(-1, 4)= +AA01 *2.0*b*COS(ph1-ph2)
      DDevpl2(-1, 5)= +AA01 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevpl2(-2, 1)= +AA01 *2.0*b*COS(   -ph1)
      DDevpl2(-2, 2)= +AA01 *2.0*a*COS(   -ph1)
      DDevpl2(-2, 3)= +AA01 *2.0*a*b*SIN(   -ph1)*D2R

      DDevpl0( 0, 1)= +AA05 *2.0*a
      DDevpl0( 0, 2)= +AA05 *2.0*b
      DDevpl0( 0, 4)= +AA05 *2.0*c

      DDevspl1( 1, 1)= +AA05 *2.0*d*COS(   -ph3)
      DDevspl1( 1, 6)= +AA05 *2.0*a*COS(   -ph3)
      DDevspl1( 1, 7)= +AA05 *2.0*a*d*SIN(   -ph3)*D2R

      DDevspl1( 0, 4)= +AA05 *2.0*d*COS(ph2-ph3)
      DDevspl1( 0, 5)= -AA05 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevspl1( 0, 6)= +AA05 *2.0*c*COS(ph2-ph3)
      DDevspl1( 0, 7)= +AA05 *2.0*c*d*SIN(ph2-ph3)*D2R

      DDevspl1(-1, 2)= +AA05 *2.0*d*COS(ph1-ph3)
      DDevspl1(-1, 3)= -AA05 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevspl1(-1, 6)= +AA05 *2.0*b*COS(ph1-ph3)
      DDevspl1(-1, 7)= +AA05 *2.0*b*d*SIN(ph1-ph3)*D2R

      DDevsl0( 0, 6)= +AA05 *2.0*d
       !
       !   The preceding lines of code have been automatically generated
       !
       ! if ((0<= Rad_Model).and.(Rad_Model<= 1)) then
       !    call IRJ_ITB(n_pat,-ptr(ipiof,2,n_pat),sn,IRJ,DIRJ)
       ! else if (Rad_Model < 0) then
       !    call IRJ_Slater(-Rad_Model, 2, Zeff, 2*pi*hm, IRJ, DIRJ)
       ! else
       !    call IRJ_Hydrog(Rad_Model, 1, 2, Zeff, 2*pi*hm, IRJ, DIRJ)
       ! end if
       if (Rad_Model < 0) then
          call IRJ_Slater(np, 2, Zeffp, hm, IRJ, DIRJ)
       else
          call IRJ_Hydrog(np,1,2,Zeffp, hm, IRJ, DIRJ)
       end if

       mffr=0.0
       mffi=0.0

       ! Calculation for l=2 (p terms)

       l=2

       summ=0.0
       summd=0.0
       do m=-l,l
         p= sign(1,m)
         YY= Real_Spher_Harm_Ucvec(l,ABS(m),p,huc)
         summ= summ+ YY*Devpl2(m)
         summd(1:5)= summd(1:5)+ YY*DDevpl2(m,1:5)
         if (PRESENT(mdxir)) then
         end if
       end do

       mffr= mffr- IRJ(l/2)*summ         !   The minus sign comes from the term i**l with l=2.

       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+4)= mdxir(ORBPAR:ORBPAR+4)-IRJ(l/2)*summd(1:5)  !   The minus sign comes from the term i**l with l=2.
          mdxir(ZPAR)= mdxir(ZPAR)- DIRJ(l/2)*summ !   The minus sign comes from the term i**l with l=2.
       end if

       ! Calculation for l=0 (p terms)

       l=0

       YY= 1.0/SQRT(4*PI)
       summ = YY*Devpl0(0)
       mffr= mffr+ IRJ(l/2)*summ
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+4)= mdxir(ORBPAR:ORBPAR+4)+IRJ(l/2)*YY*DDevpl0(0,1:5)
          mdxir(ZPAR)= mdxir(ZPAR) +DIRJ(l/2)*summ
       end if

       ! Calculation for l=0 (s terms)

       if (Rad_Model < 0) then
          call IRJ_Slater(ns, 0, Zeffs, hm, IRJ, DIRJ)
       else
          call IRJ_Hydrog(ns,0,0,Zeffs, hm, IRJ, DIRJ)
       end if

       l=0

       YY= 1.0/SQRT(4*PI)
       summ = YY*Devsl0(0)
       mffr= mffr+ IRJ(l/2)*summ
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR+5:ORBPAR+6)= mdxir(ORBPAR+5:ORBPAR+6)+IRJ(l/2)*YY*DDevsl0(0,6:7)
          mdxir(ZPAR+1)= mdxir(ZPAR+1) +DIRJ(l/2)*summ
       end if

       !   At this point the true mffr is 4*PI*MAG times the calculated mffr,
       ! and the derivatives are calculated with respect to normalized parameters but no with respect to
       ! input parameters
       !   We now correct the derivatives for normalization, scale (MAG), and 4*PI factor. The derivative
       ! with respect to the scale is also calculated.


       if (PRESENT(mdxir)) then
          mdxir(MAGPAR)= mffr*4.0*PI
       end if
       mffr= mffr*4.0*PI*MAG
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+6)= mdxir(ORBPAR:ORBPAR+6) *4.0*PI*MAG
          mdxir(ZPAR:ZPAR+1)= mdxir(ZPAR:ZPAR+1) *4.0*PI*MAG
          deraux= matmul(normcorr,(/mdxir(ORBPAR), mdxir(ORBPAR+1:ORBPAR+5:2)/))
          mdxir(ORBPAR)= deraux(1)
          mdxir(ORBPAR+1:ORBPAR+5:2)= deraux(2:4)
       end if

       ! Calculation for l=1 (s-p terms)

       l=1
       if (Rad_Model < 0) then
          YY= 2*Zeffs**(ns+1.5)*Zeffp**(np+1.5)/(factorial(2*ns+2)*factorial(2*np+2))**0.5
          IRJ(0)= YY* Int_Slater_Bessel(ns+np, 1, Zeffs+Zeffp, hm)
          if (Zeffp > 0.0001) then
             DIRJ(0)= (np+1.5)* IRJ(0)/Zeffp +    &
                      YY* Int_Slater_Bessel(ns+np+1, 1, Zeffs+Zeffp, hm)
          else
             DIRJ(0)= 0.0
          end if
          if (Zeffs > 0.0001) then
             DIRJ(1)= (ns+1.5)* IRJ(0)/Zeffs +    &
                      YY* Int_Slater_Bessel(ns+np+1, 1, Zeffs+Zeffp, hm)
          else
             DIRJ(1)= 0.0
          end if
       else
          call IRJ_SP(ns,Zeffp,Zeffs, hm,IRJ,DIRJ)
       end if

       summ=0
       do m=-l,l
         p= sign(1,m)
         YY= Real_Spher_Harm_Ucvec(l,ABS(m),p,huc)
         summ= summ+ YY*Devspl1(m)
         if (PRESENT(mdxii)) then
            mdxii(ORBPAR:ORBPAR+6)= mdxii(ORBPAR:ORBPAR+6)+IRJ(0)*YY*DDevspl1(m,1:7)
         end if
       end do

       mffi= mffi+ IRJ(0)*summ

       if (PRESENT(mdxii)) then
          mdxii(ZPAR:ZPAR+1)= mdxii(ZPAR:ZPAR+1)+ DIRJ*summ
       end if

       !   At this point the true mffi is 4*PI*MAG times the calculated mffi,
       ! and the derivatives are calculated with respect to normalized parameters but no with respect to
       ! input parameters
       !   We now correct the derivatives for normalization, scale (MAG), and 4*PI factor. The derivative
       !   with respect to the scale is also calculated.
       if (PRESENT(mdxii)) then
          mdxii(MAGPAR)= mffi*4.0*PI
       end if

       mffi= 4.0*PI* mffi*MAG
       if (PRESENT(mdxii)) then
          mdxii(ORBPAR:ORBPAR+6)= mdxii(ORBPAR:ORBPAR+6) *4.0*PI*MAG
          mdxii(ZPAR:ZPAR+1)= mdxii(ZPAR:ZPAR+1) *4.0*PI*MAG
          deraux= matmul(normcorr,(/mdxii(ORBPAR), mdxii(ORBPAR+1:ORBPAR+5:2)/))
          mdxii(ORBPAR)= deraux(1)
          mdxii(ORBPAR+1:ORBPAR+5:2)= deraux(2:4)
       end if
       return
    End Subroutine disp_ffactor

    Subroutine Quad_ffactor(n_pat,ipiof,sn,huc,mffr,mffi,mdxir,mdxii)
       integer, intent(in)                                :: n_pat
       integer, intent(in)                                :: ipiof
       real(kind=cp), intent(in)                          :: sn
       real(kind=cp), intent(in) , dimension(3)           :: huc
       real(kind=cp), intent(out)                         :: mffr
       real(kind=cp), intent(out)                         :: mffi
       real(kind=cp), optional,intent(out), dimension(:)  :: mdxir
       real(kind=cp), optional,intent(out), dimension(:)  :: mdxii
       !-----------------------------------------------
       !   L o c a l   P a r a m e t e r s
       !-----------------------------------------------

       real(kind=cp), parameter :: D2R= 0.01745329251994
       integer, parameter :: ZPAR= 29
       integer, parameter :: ORBPAR= 20
       integer, parameter :: MAGPAR= 19

       !
       !   The following lines of code have been automatically generated
       !
      real(kind=cp), parameter ::                                                       &
                         AA01=  0.2384136320, AA02= -0.2384136320, AA03= 0.1685838850,  &
                         AA04= -0.1685838850, AA05=  0.1802237630, AA06=-0.1802237630,  &
                         AA07=  0.1560783540, AA08= -0.0637187213, AA09= 0.2207281290,  &
                         AA10=  0.0402992591, AA11= -0.1611970370, AA12= 0.2417955550,  &
                         AA13=  0.0637187213, AA14= -0.1560783540, AA15= 0.0901118815,  &
                         AA16=  0.2820947770
       !
       !   The preceding lines of code have been automatically generated
       !
       !-----------------------------------------------------------------------------------------------------
       !   L o c a l   V a r i a b l e s :   R a d i a l   P a r t   o f   t h e   W a v e   F u n c t i o n
       !-----------------------------------------------------------------------------------------------------
       real(kind=cp) :: Zeff          !   Exponent and power in the Slater-Type function if Rad_Model >= 3
       integer       ::  Rad_Model !   R(r)~ r**n * exp(-chi*r)    chi= Zeff/n (n=Rad_Model)
       !-----------------------------------------------------------------------------------------------------
       !   L o c a l   V a r i a b l e s :   A n g u l a r   P a r t   o f   t h e   W a v e   F u n c t i o n
       !-----------------------------------------------------------------------------------------------------
       real(kind=cp), dimension(0:0)       :: Devdl0   !   Term with l=0 comming from d-terms
       real(kind=cp), dimension (0:0,1:9)  :: DDevdl0  !  Derivatives of the preceding term with respect
                                                       !  to normalized input parameters
       real(kind=cp), dimension (-2:2)     :: Devdl2   !   Term with l=2 comming from d-terms
       real(kind=cp), dimension (-2:2,1:9) :: DDevdl2  !  Derivatives of the preceding term with respect
                                                       !  to normalized input parameters
       real(kind=cp), dimension (-4:4)     :: Devdl4   !   Term with l=4 comming from d-terms
       real(kind=cp), dimension (-4:4,1:9) :: DDevdl4  !  Derivatives of the preceding term with respect
                                                       !  to normalized input parameters
       real(kind=cp) :: a, b, c, d, e, ph1, ph2, ph3, ph4
       real(kind=cp) :: norm, mag
       real(kind=cp), dimension (5,5) :: normcorr      !   Term with the correction to the derivatives
       real(kind=cp), dimension (5)   :: deraux        !   due to normalization
       !
       !
       !      \phi = R_p(r) \left\{ a*y_2^{1 +}+b \exp\{i ph1\}*y_2^{1 -} +c \exp\{i ph2\}*y_2^{0 +} \right\}
       !
       !      \rho= ?\phi?^2= R(r)^2*\sum_l={0,2,4} \sum_{m=-l}^l DevdlX(m) *y_l^{abs(m),sign(1,m)}(\hat{r})   (X=0,2,4)
       !
       !      fm(\vec{h})= 4\pi \sum_l={0,2,4} (i)**l*\int_{0}^{\infty} r**2*R(r)^2*j_l(2\pi h r)*dr *
       !                        \sum_{m=-l}^l DevdlX(m) *y_l^{abs(m),sign(1,m)}(\hat{h})                     (X=0,2,4)
       !
       !-----------------------------------------------
       !   L o c a l   V a r i a b l e s
       !-----------------------------------------------
       real(kind=cp)                 :: YY, summ
       real(kind=cp), dimension(9)   :: summd
       real(kind=cp), dimension(0:2) :: IRJ, DIRJ
       real(kind=cp):: hm     ! modulus of q=2*pi*|h|
       integer :: l, m, p

       hm=2.0*tpi*sqrt(sn)  !4*pi*sintheta/lambda
       Rad_Model= den_nlm(sp_pointers(ipiof))%nl(2)
       DDevdl0= 0.0
       DDevdl2= 0.0
       DDevdl4= 0.0
       if (present(mdxir)) then
          mdxir(ZPAR)=0.0
          mdxir(ORBPAR:ORBPAR+8)=0.0
       end if
       if (present(mdxii)) then
          mdxii(ZPAR)=0.0
          mdxii(ORBPAR:ORBPAR+8)=0.0
       end if
       xl(ipiof, ORBPAR+ 2:ORBPAR+ 8:2)= MOD(xl(ipiof, ORBPAR+ 2:ORBPAR+ 8:2), 360.0_cp)

       Zeff= xl(ipiof, ZPAR)
       a  = xl(ipiof, ORBPAR)
       b  = xl(ipiof, ORBPAR+1)
       ph1= xl(ipiof, ORBPAR+2) * D2R
       c  = xl(ipiof, ORBPAR+3)
       ph2= xl(ipiof, ORBPAR+4) * D2R
       d  = xl(ipiof, ORBPAR+5)
       ph3= xl(ipiof, ORBPAR+6) * D2R
       e  = xl(ipiof, ORBPAR+7)
       ph4= xl(ipiof, ORBPAR+8) * D2R
       mag= xl(ipiof, MAGPAR)
       norm= sqrt(a**2+b**2+c**2+d**2+e**2)

         !! If norm = 0 We have no things to do !!

       if (norm <= 1.e-5) then
          call Error_Message("Zero norm of an orbital, nothing is calculated", 6)
          mffr= 0.0
          mffi= 0.0
          if (PRESENT(mdxir)) then
             mdxir= 0.0
          end if
          if (PRESENT(mdxii)) then
             mdxii= 0.0
          end if
          return
       end if

       ! If not, we store the correction to derivatives due to normalization

       a= a/norm
       b= b/norm
       c= c/norm
       d= d/norm
       e= e/norm
       normcorr= RESHAPE ((/ 1.0-a**2, -a*b, -a*c, -a*d, -a*e, &
                             -b*a, 1.0-b**2, -b*c, -b*d, -b*e, &
                             -c*a, -c*b, 1.0-c**2, -c*d, -c*e, &
                             -d*a, -d*b, -d*c, 1.0-d**2, -d*e, &
                             -e*a, -e*b, -e*c, -e*d, 1.0-e**2 /),(/5,5/))/norm
       !
       !   The following lines of code have been automatically generated
       !
      Devdl4( 4)= +AA01 *a*a +AA02 *b*b
      DDevdl4( 4, 1)= +AA01 *2.0*a
      DDevdl4( 4, 2)= +AA02 *2.0*b
      Devdl4( 3)= +AA03 *2.0*a*c*COS(   -ph2) +AA04 *2.0*b*d*COS(ph1-ph3)
      DDevdl4( 3, 1)= +AA03 *2.0*c*COS(   -ph2)
      DDevdl4( 3, 2)= +AA04 *2.0*d*COS(ph1-ph3)
      DDevdl4( 3, 3)= -AA04 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevdl4( 3, 4)= +AA03 *2.0*a*COS(   -ph2)
      DDevdl4( 3, 5)= +AA03 *2.0*a*c*SIN(   -ph2)*D2R
      DDevdl4( 3, 6)= +AA04 *2.0*b*COS(ph1-ph3)
      DDevdl4( 3, 7)= +AA04 *2.0*b*d*SIN(ph1-ph3)*D2R
      Devdl4( 2)= +AA05 *c*c +AA06 *d*d +AA07 *2.0*a*e*COS(   -ph4)
      DDevdl4( 2, 1)= +AA07 *2.0*e*COS(   -ph4)
      DDevdl4( 2, 4)= +AA05 *2.0*c
      DDevdl4( 2, 6)= +AA06 *2.0*d
      DDevdl4( 2, 8)= +AA07 *2.0*a*COS(   -ph4)
      DDevdl4( 2, 9)= +AA07 *2.0*a*e*SIN(   -ph4)*D2R
      Devdl4( 1)= +AA08 *2.0*a*c*COS(   -ph2) +AA08 *2.0*b*d*COS(ph1-ph3) +AA09 *2.0*c*e*COS(ph2-ph4)
      DDevdl4( 1, 1)= +AA08 *2.0*c*COS(   -ph2)
      DDevdl4( 1, 2)= +AA08 *2.0*d*COS(ph1-ph3)
      DDevdl4( 1, 3)= -AA08 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevdl4( 1, 4)= +AA08 *2.0*a*COS(   -ph2) +AA09 *2.0*e*COS(ph2-ph4)
      DDevdl4( 1, 5)= +AA08 *2.0*a*c*SIN(   -ph2)*D2R -AA09 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevdl4( 1, 6)= +AA08 *2.0*b*COS(ph1-ph3)
      DDevdl4( 1, 7)= +AA08 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevdl4( 1, 8)= +AA09 *2.0*c*COS(ph2-ph4)
      DDevdl4( 1, 9)= +AA09 *2.0*c*e*SIN(ph2-ph4)*D2R
      Devdl4( 0)= +AA10 *a*a +AA10 *b*b +AA11 *c*c +AA11 *d*d +AA12 *e*e
      DDevdl4( 0, 1)= +AA10 *2.0*a
      DDevdl4( 0, 2)= +AA10 *2.0*b
      DDevdl4( 0, 4)= +AA11 *2.0*c
      DDevdl4( 0, 6)= +AA11 *2.0*d
      DDevdl4( 0, 8)= +AA12 *2.0*e
      Devdl4(-1)= +AA13 *2.0*a*d*COS(   -ph3) +AA08 *2.0*b*c*COS(ph1-ph2) +AA09 *2.0*d*e*COS(ph3-ph4)
      DDevdl4(-1, 1)= +AA13 *2.0*d*COS(   -ph3)
      DDevdl4(-1, 2)= +AA08 *2.0*c*COS(ph1-ph2)
      DDevdl4(-1, 3)= -AA08 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevdl4(-1, 4)= +AA08 *2.0*b*COS(ph1-ph2)
      DDevdl4(-1, 5)= +AA08 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevdl4(-1, 6)= +AA13 *2.0*a*COS(   -ph3) +AA09 *2.0*e*COS(ph3-ph4)
      DDevdl4(-1, 7)= +AA13 *2.0*a*d*SIN(   -ph3)*D2R -AA09 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevdl4(-1, 8)= +AA09 *2.0*d*COS(ph3-ph4)
      DDevdl4(-1, 9)= +AA09 *2.0*d*e*SIN(ph3-ph4)*D2R
      Devdl4(-2)= +AA07 *2.0*b*e*COS(ph1-ph4) +AA05 *2.0*c*d*COS(ph2-ph3)
      DDevdl4(-2, 2)= +AA07 *2.0*e*COS(ph1-ph4)
      DDevdl4(-2, 3)= -AA07 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevdl4(-2, 4)= +AA05 *2.0*d*COS(ph2-ph3)
      DDevdl4(-2, 5)= -AA05 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevdl4(-2, 6)= +AA05 *2.0*c*COS(ph2-ph3)
      DDevdl4(-2, 7)= +AA05 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevdl4(-2, 8)= +AA07 *2.0*b*COS(ph1-ph4)
      DDevdl4(-2, 9)= +AA07 *2.0*b*e*SIN(ph1-ph4)*D2R
      Devdl4(-3)= +AA03 *2.0*a*d*COS(   -ph3) +AA03 *2.0*b*c*COS(ph1-ph2)
      DDevdl4(-3, 1)= +AA03 *2.0*d*COS(   -ph3)
      DDevdl4(-3, 2)= +AA03 *2.0*c*COS(ph1-ph2)
      DDevdl4(-3, 3)= -AA03 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevdl4(-3, 4)= +AA03 *2.0*b*COS(ph1-ph2)
      DDevdl4(-3, 5)= +AA03 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevdl4(-3, 6)= +AA03 *2.0*a*COS(   -ph3)
      DDevdl4(-3, 7)= +AA03 *2.0*a*d*SIN(   -ph3)*D2R
      Devdl4(-4)= +AA01 *2.0*a*b*COS(   -ph1)
      DDevdl4(-4, 1)= +AA01 *2.0*b*COS(   -ph1)
      DDevdl4(-4, 2)= +AA01 *2.0*a*COS(   -ph1)
      DDevdl4(-4, 3)= +AA01 *2.0*a*b*SIN(   -ph1)*D2R


      Devdl2( 2)= +AA07 *c*c +AA14 *d*d +AA06 *2.0*a*e*COS(   -ph4)
      DDevdl2( 2, 1)= +AA06 *2.0*e*COS(   -ph4)
      DDevdl2( 2, 4)= +AA07 *2.0*c
      DDevdl2( 2, 6)= +AA14 *2.0*d
      DDevdl2( 2, 8)= +AA06 *2.0*a*COS(   -ph4)
      DDevdl2( 2, 9)= +AA06 *2.0*a*e*SIN(   -ph4)*D2R
      Devdl2( 1)= +AA07 *2.0*a*c*COS(   -ph2) +AA07 *2.0*b*d*COS(ph1-ph3) +AA15 *2.0*c*e*COS(ph2-ph4)
      DDevdl2( 1, 1)= +AA07 *2.0*c*COS(   -ph2)
      DDevdl2( 1, 2)= +AA07 *2.0*d*COS(ph1-ph3)
      DDevdl2( 1, 3)= -AA07 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevdl2( 1, 4)= +AA07 *2.0*a*COS(   -ph2) +AA15 *2.0*e*COS(ph2-ph4)
      DDevdl2( 1, 5)= +AA07 *2.0*a*c*SIN(   -ph2)*D2R -AA15 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevdl2( 1, 6)= +AA07 *2.0*b*COS(ph1-ph3)
      DDevdl2( 1, 7)= +AA07 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevdl2( 1, 8)= +AA15 *2.0*c*COS(ph2-ph4)
      DDevdl2( 1, 9)= +AA15 *2.0*c*e*SIN(ph2-ph4)*D2R
      Devdl2( 0)= +AA06 *a*a +AA06 *b*b +AA15 *c*c +AA15 *d*d +AA05 *e*e
      DDevdl2( 0, 1)= +AA06 *2.0*a
      DDevdl2( 0, 2)= +AA06 *2.0*b
      DDevdl2( 0, 4)= +AA15 *2.0*c
      DDevdl2( 0, 6)= +AA15 *2.0*d
      DDevdl2( 0, 8)= +AA05 *2.0*e
      Devdl2(-1)= +AA14 *2.0*a*d*COS(   -ph3) +AA07 *2.0*b*c*COS(ph1-ph2) +AA15 *2.0*d*e*COS(ph3-ph4)
      DDevdl2(-1, 1)= +AA14 *2.0*d*COS(   -ph3)
      DDevdl2(-1, 2)= +AA07 *2.0*c*COS(ph1-ph2)
      DDevdl2(-1, 3)= -AA07 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevdl2(-1, 4)= +AA07 *2.0*b*COS(ph1-ph2)
      DDevdl2(-1, 5)= +AA07 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevdl2(-1, 6)= +AA14 *2.0*a*COS(   -ph3) +AA15 *2.0*e*COS(ph3-ph4)
      DDevdl2(-1, 7)= +AA14 *2.0*a*d*SIN(   -ph3)*D2R -AA15 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevdl2(-1, 8)= +AA15 *2.0*d*COS(ph3-ph4)
      DDevdl2(-1, 9)= +AA15 *2.0*d*e*SIN(ph3-ph4)*D2R
      Devdl2(-2)= +AA06 *2.0*b*e*COS(ph1-ph4) +AA07 *2.0*c*d*COS(ph2-ph3)
      DDevdl2(-2, 2)= +AA06 *2.0*e*COS(ph1-ph4)
      DDevdl2(-2, 3)= -AA06 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevdl2(-2, 4)= +AA07 *2.0*d*COS(ph2-ph3)
      DDevdl2(-2, 5)= -AA07 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevdl2(-2, 6)= +AA07 *2.0*c*COS(ph2-ph3)
      DDevdl2(-2, 7)= +AA07 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevdl2(-2, 8)= +AA06 *2.0*b*COS(ph1-ph4)
      DDevdl2(-2, 9)= +AA06 *2.0*b*e*SIN(ph1-ph4)*D2R


      Devdl0( 0)= +AA16 *a*a +AA16 *b*b +AA16 *c*c +AA16 *d*d +AA16 *e*e
      DDevdl0( 0, 1)= +AA16 *2.0*a
      DDevdl0( 0, 2)= +AA16 *2.0*b
      DDevdl0( 0, 4)= +AA16 *2.0*c
      DDevdl0( 0, 6)= +AA16 *2.0*d
      DDevdl0( 0, 8)= +AA16 *2.0*e
       !
       !   The preceding lines of code have been automatically generated
       !

       if ((0 <= Rad_Model) .and. (Rad_Model <= 2)) then
          call IRJ_ITB(n_pat,-ptr(ipiof,2,n_pat),sn,IRJ,DIRJ)
       else if (Rad_Model < 0) then
          call IRJ_Slater(-Rad_Model, 4, Zeff, hm, IRJ, DIRJ)
       else
          call IRJ_Hydrog(Rad_Model,2,4,Zeff, hm, IRJ, DIRJ)
       end if
       mffr=0.0
       mffi=0.0

       ! Calculation for l=4
       l=4

       summ= 0.0
       summd= 0.0
       do m=-l,l
         p= sign(1,m)
         YY= Real_Spher_Harm_Ucvec(l,ABS(m),p,huc)
         summ= summ+ YY*Devdl4(m)
         summd= summd + YY*DDevdl4(m,:)
       end do

       mffr= mffr+ IRJ(l/2)*summ

       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+8)= mdxir(ORBPAR:ORBPAR+8)+IRJ(l/2)*summd
          mdxir(ZPAR)= mdxir(ZPAR)+ DIRJ(l/2)*summ
       end if

       ! Calculation for l=2
       l=2

       summ= 0.0
       summd= 0.0
       do m=-l,l
         p= sign(1,m)
         YY= Real_Spher_Harm_Ucvec(l,ABS(m),p,huc)
         summ= summ+ YY*Devdl2(m)
         summd= summd+ YY*DDevdl2(m,:)
       end do

       mffr= mffr- IRJ(l/2)*summ                !   The minus sign comes from the term i**l with l=2.

       if (PRESENT(mdxir)) then
         mdxir(ORBPAR:ORBPAR+8)= mdxir(ORBPAR:ORBPAR+8)-IRJ(l/2)*summd  !   The minus sign comes from the term i**l with l=2.
         mdxir(ZPAR)= mdxir(ZPAR)- DIRJ(l/2)*summ   !   The minus sign comes from the term i**l with l=2.
       end if

       ! Calculation for l=0
       l=0

       YY= 1.0/SQRT(4*PI)
       summ = YY*Devdl0(0)
       mffr= mffr+ IRJ(l/2)*summ
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+8)= mdxir(ORBPAR:ORBPAR+8)+ DIRJ(l/2)*YY*DDevdl0(0,1:9)
          mdxir(ZPAR)= mdxir(ZPAR) +DIRJ(l/2)*summ
       end if

       ! At this point the true mffr and mffi are 4*PI*MAG times the calculated
       ! mffr and mffi, and the derivatives are calculated with respect to normalized
       ! parameters but no with respect to input parameters
       ! We now correct the derivatives for normalization, scale (MAG), and 4*PI factor.
       ! The derivative with respect to the scale is also calculated.

       if (PRESENT(mdxir)) then
          mdxir(MAGPAR)= mffr*4.0*PI
       end if
       mffr= mffr*4.0*PI*MAG
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+8)= mdxir(ORBPAR:ORBPAR+8) *4.0*PI*MAG
          mdxir(ZPAR)= mdxir(ZPAR) *4.0*PI*MAG
          deraux= matmul(normcorr,(/mdxir(ORBPAR), mdxir(ORBPAR+1:ORBPAR+7:2)/))
          mdxir(ORBPAR)= deraux(1)
          mdxir(ORBPAR+1:ORBPAR+7:2)= deraux(2:5)
       end if
       return
    End Subroutine quad_ffactor

    Subroutine Hxap_Ffactor(n_pat,ipiof,sn,huc,mffr,mffi,mdxir,mdxii)
       integer,                              intent(in)  :: n_pat
       integer,                              intent(in)  :: ipiof
       real(kind=cp),                        intent(in)  :: sn
       real(kind=cp), dimension(3) ,         intent(in)  :: huc
       real(kind=cp),                        intent(out) :: mffr
       real(kind=cp),                        intent(out) :: mffi
       real(kind=cp), dimension(:), optional,intent(out) :: mdxir
       real(kind=cp), dimension(:), optional,intent(out) :: mdxii
       !-----------------------------------------------
       !   L o c a l   P a r a m e t e r s
       !-----------------------------------------------

       real(kind=cp), parameter :: D2R= 0.01745329251994
       integer, parameter :: ZPAR= 33
       integer, parameter :: ORBPAR= 20
       integer, parameter :: MAGPAR= 19

       !
       !   The following lines of code have been automatically generated
       !

      real(kind=cp), PARAMETER ::                                                              &
                         AA01=  0.2548005880, AA02= -0.2548005880, AA03=  0.1801712210,  &
                         AA04= -0.1801712210, AA05=  0.1881827120, AA06= -0.1881827120,  &
                         AA07=  0.1214714200, AA08= -0.1214714200, AA09=  0.1086473390,  &
                         AA10=  0.1629710050, AA11= -0.1629710050, AA12=  0.1717865170,  &
                         AA13= -0.1717865170, AA14= -0.0443550907, AA15=  0.1774203480,  &
                         AA16=  0.0221775454, AA17= -0.0858932585, AA18=  0.2217754420,  &
                         AA19= -0.0118543962, AA20=  0.0711263791, AA21= -0.1778159290,  &
                         AA22=  0.2370879200, AA23= -0.0221775454, AA24=  0.0858932585,  &
                         AA25=  0.0443550907, AA26=  0.1517177520, AA27= -0.1517177520,  &
                         AA28= -0.1175200640, AA29=  0.1175200640, AA30= -0.2035507110,  &
                         AA31=  0.0678502396, AA32= -0.0678502396, AA33=  0.1146878450,  &
                         AA34= -0.1146878450, AA35=  0.1332552280, AA36= -0.0444184095,  &
                         AA37= -0.0993225873, AA38=  0.1025799290, AA39=  0.0993225798,  &
                         AA40=  0.0769349411, AA41= -0.1795148700, AA42=  0.0256449804,  &
                         AA43=  0.1538698820, AA44= -0.1025799290, AA45= -0.1332552280,  &
                         AA46=  0.1456731260, AA47= -0.1456731260, AA48= -0.0940315947,  &
                         AA49= -0.1880631890, AA50=  0.1486770060, AA51=  0.1151647120,  &
                         AA52=  0.0594708063, AA53= -0.2102610470, AA54=  0.1261566130,  &
                         AA55=  0.1682088380, AA56= -0.1486770060, AA57= -0.1151647120,  &
                         AA58=  0.0940315947, AA59=  0.2820948060
       !
       !   The preceding lines of code have been automatically generated
       !
       !-----------------------------------------------------------------------------------------------------
       !   L o c a l   V a r i a b l e s :   R a d i a l   P a r t   o f   t h e   W a v e   F u n c t i o n
       !-----------------------------------------------------------------------------------------------------
       real(kind=cp) :: Zeff          !   exponent and power in the Slater-Type function   if Rad_Model >= 4
       integer :: Rad_Model  !   R(r)~ r**n * exp(-chi*r)         chi= Zeff/n    (n=Rad_Model)
       !-----------------------------------------------------------------------------------------------------
       !   L o c a l   V a r i a b l e s :   A n g u l a r   P a r t   o f   t h e   W a v e   F u n c t i o n
       !-----------------------------------------------------------------------------------------------------
       real(kind=cp), dimension (0:0)       :: Devfl0  !   Term with l=0 comming from d-terms
       real(kind=cp), dimension (0:0,1:13)  :: DDevfl0 !  Derivatives of the preceding term with respect
                                                       !  to normalized input parameters
       real(kind=cp), dimension (-2:2)      :: Devfl2  !   Term with l=2 comming from d-terms
       real(kind=cp), dimension (-2:2,1:13) :: DDevfl2 !  Derivatives of the preceding term with respect
                                                       !  to normalized input parameters
       real(kind=cp), dimension (-4:4)      :: Devfl4  !   Term with l=4 comming from d-terms
       real(kind=cp), dimension (-4:4,1:13) :: DDevfl4 !  Derivatives of the preceding term with respect
                                                       !  to normalized input parameters
       real(kind=cp), dimension (-6:6)      :: Devfl6  !   Term with l=4 comming from d-terms
       real(kind=cp), dimension (-6:6,1:13) :: DDevfl6 !  Derivatives of the preceding term with respect
                                                       !  to normalized input parameters
       real(kind=cp) :: a, b, c, d, e, f, g, ph1, ph2, ph3, ph4, ph5, ph6
       real(kind=cp) :: norm, mag
       real(kind=cp), dimension (7,7) :: normcorr  !   Term with the correction to the derivatives
       real(kind=cp), dimension (7)   :: deraux        !   due to normalization
       !
       !
       !      \phi = R_f(r) \left\{ a*y_3^{3 +}+b \exp\{i ph1\}*y_3^{3 -} +c \exp\{i ph2\}*y_3^{2 +} +
       !                            d \exp\{i ph3\}*y_3^{2 -}+e \exp\{i ph4\}*y_3^{1 +}+
       !                            f \exp\{i ph5\}*y_3^{1 -} + g \exp\{i ph6\}*y_3^{0 +} \right\}
       !
       !      \rho= ?\phi?^2= R(r)^2*\sum_l={0,2,4,6} \sum_{m=-l}^l DevdlX(m) *y_l^{abs(m),sign(1,m)}(\hat{r})   (X=0,2,4,6)
       !
       !      fm(\vec{h})= 4\pi \sum_l={0,2,4,6} (i)**l*\int_{0}^{\infty} r**2*R(r)^2*j_l(2\pi h r)*dr *
       !                        \sum_{m=-l}^l DevdlX(m) *y_l^{abs(m),sign(1,m)}(\hat{h})                     (X=0,2,4,6)
       !
       !-----------------------------------------------
       !   L o c a l   V a r i a b l e s
       !-----------------------------------------------
       real(kind=cp) ::  YY, summ
       real(kind=cp), dimension(13) ::  summd
       real(kind=cp), dimension (0:3) :: IRJ, DIRJ
       real(kind=cp):: hm     ! modulus of q=2*pi*|h|
       integer      :: l, m, p

       hm=2.0*tpi*sqrt(sn)

       Rad_Model= den_nlm(sp_pointers(ipiof))%nl(2)

       DDevfl0= 0.0
       DDevfl2= 0.0
       DDevfl4= 0.0
       DDevfl6= 0.0
       ! Phases modulo 360 degrees
       xl(ipiof, ORBPAR+ 2:ORBPAR+12:2)= MOD(xl(ipiof, ORBPAR+ 2:ORBPAR+12:2), 360.0_cp)
       ! Effective charge
       Zeff= xl(ipiof, ZPAR)
       ! Read all the parameters a=A33+, b=A33-, ph1=Phi(33-), c=A32+.......
       ! and transform phases in radians
       a  = xl(ipiof, ORBPAR)
       b  = xl(ipiof, ORBPAR+1)
       ph1= xl(ipiof, ORBPAR+2) * D2R
       c  = xl(ipiof, ORBPAR+3)
       ph2= xl(ipiof, ORBPAR+4) * D2R
       d  = xl(ipiof, ORBPAR+5)
       ph3= xl(ipiof, ORBPAR+6) * D2R
       e  = xl(ipiof, ORBPAR+7)
       ph4= xl(ipiof, ORBPAR+8) * D2R
       f  = xl(ipiof, ORBPAR+9)
       ph5= xl(ipiof, ORBPAR+10) * D2R
       g  = xl(ipiof, ORBPAR+11)
       ph6= xl(ipiof, ORBPAR+12) * D2R
       mag= xl(ipiof, MAGPAR)
       !
       norm= sqrt(a**2+b**2+c**2+d**2+e**2+f**2+g**2)

       !! If norm = 0 We have nothing to do !!
       if (norm <= 1.e-5) then
          call Error_Message("Zero norm of an orbital, nothing is calculated", 6)
          mffr= 0.0
          mffi= 0.0
          if (PRESENT(mdxir)) then
             mdxir= 0.0
          end if
          if (PRESENT(mdxii)) then
             mdxii= 0.0
          end if
          return
       end if

       !! If not, we store the correction to derivatives due to normalization

       a= a/norm
       b= b/norm
       c= c/norm
       d= d/norm
       e= e/norm
       f= f/norm
       g= g/norm
       normcorr= RESHAPE ((/ 1.0-a**2, -a*b, -a*c, -a*d, -a*e, -a*f, -a*g, &
                             -b*a, 1.0-b**2, -b*c, -b*d, -b*e, -b*f, -b*g, &
                             -c*a, -c*b, 1.0-c**2, -c*d, -c*e, -c*f, -c*g, &
                             -d*a, -d*b, -d*c, 1.0-d**2, -d*e, -d*f, -d*g, &
                             -e*a, -e*b, -e*c, -e*d, 1.0-e**2, -e*f, -e*g, &
                             -f*a, -f*b, -f*c, -f*d, -f*e, 1.0-f**2, -f*g, &
                             -g*a, -g*b, -g*c, -g*d, -g*e, -g*f, 1.0-e**2 /),(/7,7/))/norm
      !
      !   The following lines of code have been automatically generated
      !
      Devfl6( 6)= +AA01 *a*a +AA02 *b*b
      DDevfl6( 6, 1)= +AA01 *2.0*a
      DDevfl6( 6, 2)= +AA02 *2.0*b
      Devfl6( 5)= +AA03 *2.0*a*c*COS(   -ph2) +AA04 *2.0*b*d*COS(ph1-ph3)
      DDevfl6( 5, 1)= +AA03 *2.0*c*COS(   -ph2)
      DDevfl6( 5, 2)= +AA04 *2.0*d*COS(ph1-ph3)
      DDevfl6( 5, 3)= -AA04 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevfl6( 5, 4)= +AA03 *2.0*a*COS(   -ph2)
      DDevfl6( 5, 5)= +AA03 *2.0*a*c*SIN(   -ph2)*D2R
      DDevfl6( 5, 6)= +AA04 *2.0*b*COS(ph1-ph3)
      DDevfl6( 5, 7)= +AA04 *2.0*b*d*SIN(ph1-ph3)*D2R
      Devfl6( 4)= +AA05 *c*c +AA06 *d*d +AA07 *2.0*a*e*COS(   -ph4) +AA08 *2.0*b*f*COS(ph1-ph5)
      DDevfl6( 4, 1)= +AA07 *2.0*e*COS(   -ph4)
      DDevfl6( 4, 2)= +AA08 *2.0*f*COS(ph1-ph5)
      DDevfl6( 4, 3)= -AA08 *2.0*b*f*SIN(ph1-ph5)*D2R
      DDevfl6( 4, 4)= +AA05 *2.0*c
      DDevfl6( 4, 6)= +AA06 *2.0*d
      DDevfl6( 4, 8)= +AA07 *2.0*a*COS(   -ph4)
      DDevfl6( 4, 9)= +AA07 *2.0*a*e*SIN(   -ph4)*D2R
      DDevfl6( 4,10)= +AA08 *2.0*b*COS(ph1-ph5)
      DDevfl6( 4,11)= +AA08 *2.0*b*f*SIN(ph1-ph5)*D2R
      Devfl6( 3)= +AA09 *2.0*a*g*COS(   -ph6) +AA10 *2.0*c*e*COS(ph2-ph4) +AA11 *2.0*d*f*COS(ph3-ph5)
      DDevfl6( 3, 1)= +AA09 *2.0*g*COS(   -ph6)
      DDevfl6( 3, 4)= +AA10 *2.0*e*COS(ph2-ph4)
      DDevfl6( 3, 5)= -AA10 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevfl6( 3, 6)= +AA11 *2.0*f*COS(ph3-ph5)
      DDevfl6( 3, 7)= -AA11 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl6( 3, 8)= +AA10 *2.0*c*COS(ph2-ph4)
      DDevfl6( 3, 9)= +AA10 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevfl6( 3,10)= +AA11 *2.0*d*COS(ph3-ph5)
      DDevfl6( 3,11)= +AA11 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl6( 3,12)= +AA09 *2.0*a*COS(   -ph6)
      DDevfl6( 3,13)= +AA09 *2.0*a*g*SIN(   -ph6)*D2R
      Devfl6( 2)= +AA12 *e*e +AA13 *f*f +AA14 *2.0*a*e*COS(   -ph4) +AA14 *2.0*b*f*COS(ph1-ph5) +AA15 *2.0*c*g*COS(ph2-ph6)
      DDevfl6( 2, 1)= +AA14 *2.0*e*COS(   -ph4)
      DDevfl6( 2, 2)= +AA14 *2.0*f*COS(ph1-ph5)
      DDevfl6( 2, 3)= -AA14 *2.0*b*f*SIN(ph1-ph5)*D2R
      DDevfl6( 2, 4)= +AA15 *2.0*g*COS(ph2-ph6)
      DDevfl6( 2, 5)= -AA15 *2.0*c*g*SIN(ph2-ph6)*D2R
      DDevfl6( 2, 8)= +AA12 *2.0*e +AA14 *2.0*a*COS(   -ph4)
      DDevfl6( 2, 9)= +AA14 *2.0*a*e*SIN(   -ph4)*D2R
      DDevfl6( 2,10)= +AA13 *2.0*f +AA14 *2.0*b*COS(ph1-ph5)
      DDevfl6( 2,11)= +AA14 *2.0*b*f*SIN(ph1-ph5)*D2R
      DDevfl6( 2,12)= +AA15 *2.0*c*COS(ph2-ph6)
      DDevfl6( 2,13)= +AA15 *2.0*c*g*SIN(ph2-ph6)*D2R
      Devfl6( 1)= +AA16 *2.0*a*c*COS(   -ph2) +AA16 *2.0*b*d*COS(ph1-ph3) +AA17 *2.0*c*e*COS(ph2-ph4) +&
                  AA17 *2.0*d*f*COS(ph3-ph5) +AA18 *2.0*e*g*COS(ph4-ph6)
      DDevfl6( 1, 1)= +AA16 *2.0*c*COS(   -ph2)
      DDevfl6( 1, 2)= +AA16 *2.0*d*COS(ph1-ph3)
      DDevfl6( 1, 3)= -AA16 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevfl6( 1, 4)= +AA16 *2.0*a*COS(   -ph2) +AA17 *2.0*e*COS(ph2-ph4)
      DDevfl6( 1, 5)= +AA16 *2.0*a*c*SIN(   -ph2)*D2R -AA17 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevfl6( 1, 6)= +AA16 *2.0*b*COS(ph1-ph3) +AA17 *2.0*f*COS(ph3-ph5)
      DDevfl6( 1, 7)= +AA16 *2.0*b*d*SIN(ph1-ph3)*D2R -AA17 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl6( 1, 8)= +AA17 *2.0*c*COS(ph2-ph4) +AA18 *2.0*g*COS(ph4-ph6)
      DDevfl6( 1, 9)= +AA17 *2.0*c*e*SIN(ph2-ph4)*D2R -AA18 *2.0*e*g*SIN(ph4-ph6)*D2R
      DDevfl6( 1,10)= +AA17 *2.0*d*COS(ph3-ph5)
      DDevfl6( 1,11)= +AA17 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl6( 1,12)= +AA18 *2.0*e*COS(ph4-ph6)
      DDevfl6( 1,13)= +AA18 *2.0*e*g*SIN(ph4-ph6)*D2R
      Devfl6( 0)= +AA19 *a*a +AA19 *b*b +AA20 *c*c +AA20 *d*d +AA21 *e*e +AA21 *f*f +AA22 *g*g
      DDevfl6( 0, 1)= +AA19 *2.0*a
      DDevfl6( 0, 2)= +AA19 *2.0*b
      DDevfl6( 0, 4)= +AA20 *2.0*c
      DDevfl6( 0, 6)= +AA20 *2.0*d
      DDevfl6( 0, 8)= +AA21 *2.0*e
      DDevfl6( 0,10)= +AA21 *2.0*f
      DDevfl6( 0,12)= +AA22 *2.0*g
      Devfl6(-1)= +AA23 *2.0*a*d*COS(   -ph3) +AA16 *2.0*b*c*COS(ph1-ph2) +AA24 *2.0*c*f*COS(ph2-ph5) +&
                  AA17 *2.0*d*e*COS(ph3-ph4) +AA18 *2.0*f*g*COS(ph5-ph6)
      DDevfl6(-1, 1)= +AA23 *2.0*d*COS(   -ph3)
      DDevfl6(-1, 2)= +AA16 *2.0*c*COS(ph1-ph2)
      DDevfl6(-1, 3)= -AA16 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevfl6(-1, 4)= +AA16 *2.0*b*COS(ph1-ph2) +AA24 *2.0*f*COS(ph2-ph5)
      DDevfl6(-1, 5)= +AA16 *2.0*b*c*SIN(ph1-ph2)*D2R -AA24 *2.0*c*f*SIN(ph2-ph5)*D2R
      DDevfl6(-1, 6)= +AA23 *2.0*a*COS(   -ph3) +AA17 *2.0*e*COS(ph3-ph4)
      DDevfl6(-1, 7)= +AA23 *2.0*a*d*SIN(   -ph3)*D2R -AA17 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl6(-1, 8)= +AA17 *2.0*d*COS(ph3-ph4)
      DDevfl6(-1, 9)= +AA17 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl6(-1,10)= +AA24 *2.0*c*COS(ph2-ph5) +AA18 *2.0*g*COS(ph5-ph6)
      DDevfl6(-1,11)= +AA24 *2.0*c*f*SIN(ph2-ph5)*D2R -AA18 *2.0*f*g*SIN(ph5-ph6)*D2R
      DDevfl6(-1,12)= +AA18 *2.0*f*COS(ph5-ph6)
      DDevfl6(-1,13)= +AA18 *2.0*f*g*SIN(ph5-ph6)*D2R
      Devfl6(-2)= +AA25 *2.0*a*f*COS(   -ph5) +AA14 *2.0*b*e*COS(ph1-ph4) +AA15 *2.0*d*g*COS(ph3-ph6) +AA12 *2.0*e*f*COS(ph4-ph5)
      DDevfl6(-2, 1)= +AA25 *2.0*f*COS(   -ph5)
      DDevfl6(-2, 2)= +AA14 *2.0*e*COS(ph1-ph4)
      DDevfl6(-2, 3)= -AA14 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevfl6(-2, 6)= +AA15 *2.0*g*COS(ph3-ph6)
      DDevfl6(-2, 7)= -AA15 *2.0*d*g*SIN(ph3-ph6)*D2R
      DDevfl6(-2, 8)= +AA14 *2.0*b*COS(ph1-ph4) +AA12 *2.0*f*COS(ph4-ph5)
      DDevfl6(-2, 9)= +AA14 *2.0*b*e*SIN(ph1-ph4)*D2R -AA12 *2.0*e*f*SIN(ph4-ph5)*D2R
      DDevfl6(-2,10)= +AA25 *2.0*a*COS(   -ph5) +AA12 *2.0*e*COS(ph4-ph5)
      DDevfl6(-2,11)= +AA25 *2.0*a*f*SIN(   -ph5)*D2R +AA12 *2.0*e*f*SIN(ph4-ph5)*D2R
      DDevfl6(-2,12)= +AA15 *2.0*d*COS(ph3-ph6)
      DDevfl6(-2,13)= +AA15 *2.0*d*g*SIN(ph3-ph6)*D2R
      Devfl6(-3)= +AA09 *2.0*b*g*COS(ph1-ph6) +AA10 *2.0*c*f*COS(ph2-ph5) +AA10 *2.0*d*e*COS(ph3-ph4)
      DDevfl6(-3, 2)= +AA09 *2.0*g*COS(ph1-ph6)
      DDevfl6(-3, 3)= -AA09 *2.0*b*g*SIN(ph1-ph6)*D2R
      DDevfl6(-3, 4)= +AA10 *2.0*f*COS(ph2-ph5)
      DDevfl6(-3, 5)= -AA10 *2.0*c*f*SIN(ph2-ph5)*D2R
      DDevfl6(-3, 6)= +AA10 *2.0*e*COS(ph3-ph4)
      DDevfl6(-3, 7)= -AA10 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl6(-3, 8)= +AA10 *2.0*d*COS(ph3-ph4)
      DDevfl6(-3, 9)= +AA10 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl6(-3,10)= +AA10 *2.0*c*COS(ph2-ph5)
      DDevfl6(-3,11)= +AA10 *2.0*c*f*SIN(ph2-ph5)*D2R
      DDevfl6(-3,12)= +AA09 *2.0*b*COS(ph1-ph6)
      DDevfl6(-3,13)= +AA09 *2.0*b*g*SIN(ph1-ph6)*D2R
      Devfl6(-4)= +AA07 *2.0*a*f*COS(   -ph5) +AA07 *2.0*b*e*COS(ph1-ph4) +AA05 *2.0*c*d*COS(ph2-ph3)
      DDevfl6(-4, 1)= +AA07 *2.0*f*COS(   -ph5)
      DDevfl6(-4, 2)= +AA07 *2.0*e*COS(ph1-ph4)
      DDevfl6(-4, 3)= -AA07 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevfl6(-4, 4)= +AA05 *2.0*d*COS(ph2-ph3)
      DDevfl6(-4, 5)= -AA05 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevfl6(-4, 6)= +AA05 *2.0*c*COS(ph2-ph3)
      DDevfl6(-4, 7)= +AA05 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevfl6(-4, 8)= +AA07 *2.0*b*COS(ph1-ph4)
      DDevfl6(-4, 9)= +AA07 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevfl6(-4,10)= +AA07 *2.0*a*COS(   -ph5)
      DDevfl6(-4,11)= +AA07 *2.0*a*f*SIN(   -ph5)*D2R
      Devfl6(-5)= +AA03 *2.0*a*d*COS(   -ph3) +AA03 *2.0*b*c*COS(ph1-ph2)
      DDevfl6(-5, 1)= +AA03 *2.0*d*COS(   -ph3)
      DDevfl6(-5, 2)= +AA03 *2.0*c*COS(ph1-ph2)
      DDevfl6(-5, 3)= -AA03 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevfl6(-5, 4)= +AA03 *2.0*b*COS(ph1-ph2)
      DDevfl6(-5, 5)= +AA03 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevfl6(-5, 6)= +AA03 *2.0*a*COS(   -ph3)
      DDevfl6(-5, 7)= +AA03 *2.0*a*d*SIN(   -ph3)*D2R
      Devfl6(-6)= +AA01 *2.0*a*b*COS(   -ph1)
      DDevfl6(-6, 1)= +AA01 *2.0*b*COS(   -ph1)
      DDevfl6(-6, 2)= +AA01 *2.0*a*COS(   -ph1)
      DDevfl6(-6, 3)= +AA01 *2.0*a*b*SIN(   -ph1)*D2R


      Devfl4( 4)= +AA26 *c*c +AA27 *d*d +AA28 *2.0*a*e*COS(   -ph4) +AA29 *2.0*b*f*COS(ph1-ph5)
      DDevfl4( 4, 1)= +AA28 *2.0*e*COS(   -ph4)
      DDevfl4( 4, 2)= +AA29 *2.0*f*COS(ph1-ph5)
      DDevfl4( 4, 3)= -AA29 *2.0*b*f*SIN(ph1-ph5)*D2R
      DDevfl4( 4, 4)= +AA26 *2.0*c
      DDevfl4( 4, 6)= +AA27 *2.0*d
      DDevfl4( 4, 8)= +AA28 *2.0*a*COS(   -ph4)
      DDevfl4( 4, 9)= +AA28 *2.0*a*e*SIN(   -ph4)*D2R
      DDevfl4( 4,10)= +AA29 *2.0*b*COS(ph1-ph5)
      DDevfl4( 4,11)= +AA29 *2.0*b*f*SIN(ph1-ph5)*D2R
      Devfl4( 3)= +AA30 *2.0*a*g*COS(   -ph6) +AA31 *2.0*c*e*COS(ph2-ph4) +AA32 *2.0*d*f*COS(ph3-ph5)
      DDevfl4( 3, 1)= +AA30 *2.0*g*COS(   -ph6)
      DDevfl4( 3, 4)= +AA31 *2.0*e*COS(ph2-ph4)
      DDevfl4( 3, 5)= -AA31 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevfl4( 3, 6)= +AA32 *2.0*f*COS(ph3-ph5)
      DDevfl4( 3, 7)= -AA32 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl4( 3, 8)= +AA31 *2.0*c*COS(ph2-ph4)
      DDevfl4( 3, 9)= +AA31 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevfl4( 3,10)= +AA32 *2.0*d*COS(ph3-ph5)
      DDevfl4( 3,11)= +AA32 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl4( 3,12)= +AA30 *2.0*a*COS(   -ph6)
      DDevfl4( 3,13)= +AA30 *2.0*a*g*SIN(   -ph6)*D2R
      Devfl4( 2)= +AA33 *e*e +AA34 *f*f +AA35 *2.0*a*e*COS(   -ph4) +AA35 *2.0*b*f*COS(ph1-ph5) +AA36 *2.0*c*g*COS(ph2-ph6)
      DDevfl4( 2, 1)= +AA35 *2.0*e*COS(   -ph4)
      DDevfl4( 2, 2)= +AA35 *2.0*f*COS(ph1-ph5)
      DDevfl4( 2, 3)= -AA35 *2.0*b*f*SIN(ph1-ph5)*D2R
      DDevfl4( 2, 4)= +AA36 *2.0*g*COS(ph2-ph6)
      DDevfl4( 2, 5)= -AA36 *2.0*c*g*SIN(ph2-ph6)*D2R
      DDevfl4( 2, 8)= +AA33 *2.0*e +AA35 *2.0*a*COS(   -ph4)
      DDevfl4( 2, 9)= +AA35 *2.0*a*e*SIN(   -ph4)*D2R
      DDevfl4( 2,10)= +AA34 *2.0*f +AA35 *2.0*b*COS(ph1-ph5)
      DDevfl4( 2,11)= +AA35 *2.0*b*f*SIN(ph1-ph5)*D2R
      DDevfl4( 2,12)= +AA36 *2.0*c*COS(ph2-ph6)
      DDevfl4( 2,13)= +AA36 *2.0*c*g*SIN(ph2-ph6)*D2R
      Devfl4( 1)= +AA37 *2.0*a*c*COS(   -ph2) +AA37 *2.0*b*d*COS(ph1-ph3) +AA38 *2.0*c*e*COS(ph2-ph4) +&
                  AA38 *2.0*d*f*COS(ph3-ph5) +AA39 *2.0*e*g*COS(ph4-ph6)
      DDevfl4( 1, 1)= +AA37 *2.0*c*COS(   -ph2)
      DDevfl4( 1, 2)= +AA37 *2.0*d*COS(ph1-ph3)
      DDevfl4( 1, 3)= -AA37 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevfl4( 1, 4)= +AA37 *2.0*a*COS(   -ph2) +AA38 *2.0*e*COS(ph2-ph4)
      DDevfl4( 1, 5)= +AA37 *2.0*a*c*SIN(   -ph2)*D2R -AA38 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevfl4( 1, 6)= +AA37 *2.0*b*COS(ph1-ph3) +AA38 *2.0*f*COS(ph3-ph5)
      DDevfl4( 1, 7)= +AA37 *2.0*b*d*SIN(ph1-ph3)*D2R -AA38 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl4( 1, 8)= +AA38 *2.0*c*COS(ph2-ph4) +AA39 *2.0*g*COS(ph4-ph6)
      DDevfl4( 1, 9)= +AA38 *2.0*c*e*SIN(ph2-ph4)*D2R -AA39 *2.0*e*g*SIN(ph4-ph6)*D2R
      DDevfl4( 1,10)= +AA38 *2.0*d*COS(ph3-ph5)
      DDevfl4( 1,11)= +AA38 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl4( 1,12)= +AA39 *2.0*e*COS(ph4-ph6)
      DDevfl4( 1,13)= +AA39 *2.0*e*g*SIN(ph4-ph6)*D2R
      Devfl4( 0)= +AA40 *a*a +AA40 *b*b +AA41 *c*c +AA41 *d*d +AA42 *e*e +AA42 *f*f +AA43 *g*g
      DDevfl4( 0, 1)= +AA40 *2.0*a
      DDevfl4( 0, 2)= +AA40 *2.0*b
      DDevfl4( 0, 4)= +AA41 *2.0*c
      DDevfl4( 0, 6)= +AA41 *2.0*d
      DDevfl4( 0, 8)= +AA42 *2.0*e
      DDevfl4( 0,10)= +AA42 *2.0*f
      DDevfl4( 0,12)= +AA43 *2.0*g
      Devfl4(-1)= +AA39 *2.0*a*d*COS(   -ph3) +AA37 *2.0*b*c*COS(ph1-ph2) +AA44 *2.0*c*f*COS(ph2-ph5) +&
                  AA38 *2.0*d*e*COS(ph3-ph4) +AA39 *2.0*f*g*COS(ph5-ph6)
      DDevfl4(-1, 1)= +AA39 *2.0*d*COS(   -ph3)
      DDevfl4(-1, 2)= +AA37 *2.0*c*COS(ph1-ph2)
      DDevfl4(-1, 3)= -AA37 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevfl4(-1, 4)= +AA37 *2.0*b*COS(ph1-ph2) +AA44 *2.0*f*COS(ph2-ph5)
      DDevfl4(-1, 5)= +AA37 *2.0*b*c*SIN(ph1-ph2)*D2R -AA44 *2.0*c*f*SIN(ph2-ph5)*D2R
      DDevfl4(-1, 6)= +AA39 *2.0*a*COS(   -ph3) +AA38 *2.0*e*COS(ph3-ph4)
      DDevfl4(-1, 7)= +AA39 *2.0*a*d*SIN(   -ph3)*D2R -AA38 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl4(-1, 8)= +AA38 *2.0*d*COS(ph3-ph4)
      DDevfl4(-1, 9)= +AA38 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl4(-1,10)= +AA44 *2.0*c*COS(ph2-ph5) +AA39 *2.0*g*COS(ph5-ph6)
      DDevfl4(-1,11)= +AA44 *2.0*c*f*SIN(ph2-ph5)*D2R -AA39 *2.0*f*g*SIN(ph5-ph6)*D2R
      DDevfl4(-1,12)= +AA39 *2.0*f*COS(ph5-ph6)
      DDevfl4(-1,13)= +AA39 *2.0*f*g*SIN(ph5-ph6)*D2R
      Devfl4(-2)= +AA45 *2.0*a*f*COS(   -ph5) +AA35 *2.0*b*e*COS(ph1-ph4) +AA36 *2.0*d*g*COS(ph3-ph6) +AA33 *2.0*e*f*COS(ph4-ph5)
      DDevfl4(-2, 1)= +AA45 *2.0*f*COS(   -ph5)
      DDevfl4(-2, 2)= +AA35 *2.0*e*COS(ph1-ph4)
      DDevfl4(-2, 3)= -AA35 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevfl4(-2, 6)= +AA36 *2.0*g*COS(ph3-ph6)
      DDevfl4(-2, 7)= -AA36 *2.0*d*g*SIN(ph3-ph6)*D2R
      DDevfl4(-2, 8)= +AA35 *2.0*b*COS(ph1-ph4) +AA33 *2.0*f*COS(ph4-ph5)
      DDevfl4(-2, 9)= +AA35 *2.0*b*e*SIN(ph1-ph4)*D2R -AA33 *2.0*e*f*SIN(ph4-ph5)*D2R
      DDevfl4(-2,10)= +AA45 *2.0*a*COS(   -ph5) +AA33 *2.0*e*COS(ph4-ph5)
      DDevfl4(-2,11)= +AA45 *2.0*a*f*SIN(   -ph5)*D2R +AA33 *2.0*e*f*SIN(ph4-ph5)*D2R
      DDevfl4(-2,12)= +AA36 *2.0*d*COS(ph3-ph6)
      DDevfl4(-2,13)= +AA36 *2.0*d*g*SIN(ph3-ph6)*D2R
      Devfl4(-3)= +AA30 *2.0*b*g*COS(ph1-ph6) +AA31 *2.0*c*f*COS(ph2-ph5) +AA31 *2.0*d*e*COS(ph3-ph4)
      DDevfl4(-3, 2)= +AA30 *2.0*g*COS(ph1-ph6)
      DDevfl4(-3, 3)= -AA30 *2.0*b*g*SIN(ph1-ph6)*D2R
      DDevfl4(-3, 4)= +AA31 *2.0*f*COS(ph2-ph5)
      DDevfl4(-3, 5)= -AA31 *2.0*c*f*SIN(ph2-ph5)*D2R
      DDevfl4(-3, 6)= +AA31 *2.0*e*COS(ph3-ph4)
      DDevfl4(-3, 7)= -AA31 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl4(-3, 8)= +AA31 *2.0*d*COS(ph3-ph4)
      DDevfl4(-3, 9)= +AA31 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl4(-3,10)= +AA31 *2.0*c*COS(ph2-ph5)
      DDevfl4(-3,11)= +AA31 *2.0*c*f*SIN(ph2-ph5)*D2R
      DDevfl4(-3,12)= +AA30 *2.0*b*COS(ph1-ph6)
      DDevfl4(-3,13)= +AA30 *2.0*b*g*SIN(ph1-ph6)*D2R
      Devfl4(-4)= +AA28 *2.0*a*f*COS(   -ph5) +AA28 *2.0*b*e*COS(ph1-ph4) +AA26 *2.0*c*d*COS(ph2-ph3)
      DDevfl4(-4, 1)= +AA28 *2.0*f*COS(   -ph5)
      DDevfl4(-4, 2)= +AA28 *2.0*e*COS(ph1-ph4)
      DDevfl4(-4, 3)= -AA28 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevfl4(-4, 4)= +AA26 *2.0*d*COS(ph2-ph3)
      DDevfl4(-4, 5)= -AA26 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevfl4(-4, 6)= +AA26 *2.0*c*COS(ph2-ph3)
      DDevfl4(-4, 7)= +AA26 *2.0*c*d*SIN(ph2-ph3)*D2R
      DDevfl4(-4, 8)= +AA28 *2.0*b*COS(ph1-ph4)
      DDevfl4(-4, 9)= +AA28 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevfl4(-4,10)= +AA28 *2.0*a*COS(   -ph5)
      DDevfl4(-4,11)= +AA28 *2.0*a*f*SIN(   -ph5)*D2R


      Devfl2( 2)= +AA46 *e*e +AA47 *f*f +AA48 *2.0*a*e*COS(   -ph4) +AA48 *2.0*b*f*COS(ph1-ph5) +AA49 *2.0*c*g*COS(ph2-ph6)
      DDevfl2( 2, 1)= +AA48 *2.0*e*COS(   -ph4)
      DDevfl2( 2, 2)= +AA48 *2.0*f*COS(ph1-ph5)
      DDevfl2( 2, 3)= -AA48 *2.0*b*f*SIN(ph1-ph5)*D2R
      DDevfl2( 2, 4)= +AA49 *2.0*g*COS(ph2-ph6)
      DDevfl2( 2, 5)= -AA49 *2.0*c*g*SIN(ph2-ph6)*D2R
      DDevfl2( 2, 8)= +AA46 *2.0*e +AA48 *2.0*a*COS(   -ph4)
      DDevfl2( 2, 9)= +AA48 *2.0*a*e*SIN(   -ph4)*D2R
      DDevfl2( 2,10)= +AA47 *2.0*f +AA48 *2.0*b*COS(ph1-ph5)
      DDevfl2( 2,11)= +AA48 *2.0*b*f*SIN(ph1-ph5)*D2R
      DDevfl2( 2,12)= +AA49 *2.0*c*COS(ph2-ph6)
      DDevfl2( 2,13)= +AA49 *2.0*c*g*SIN(ph2-ph6)*D2R
      Devfl2( 1)= +AA50 *2.0*a*c*COS(   -ph2) +AA50 *2.0*b*d*COS(ph1-ph3) +AA51 *2.0*c*e*COS(ph2-ph4) +&
                  AA51 *2.0*d*f*COS(ph3-ph5) +AA52 *2.0*e*g*COS(ph4-ph6)
      DDevfl2( 1, 1)= +AA50 *2.0*c*COS(   -ph2)
      DDevfl2( 1, 2)= +AA50 *2.0*d*COS(ph1-ph3)
      DDevfl2( 1, 3)= -AA50 *2.0*b*d*SIN(ph1-ph3)*D2R
      DDevfl2( 1, 4)= +AA50 *2.0*a*COS(   -ph2) +AA51 *2.0*e*COS(ph2-ph4)
      DDevfl2( 1, 5)= +AA50 *2.0*a*c*SIN(   -ph2)*D2R -AA51 *2.0*c*e*SIN(ph2-ph4)*D2R
      DDevfl2( 1, 6)= +AA50 *2.0*b*COS(ph1-ph3) +AA51 *2.0*f*COS(ph3-ph5)
      DDevfl2( 1, 7)= +AA50 *2.0*b*d*SIN(ph1-ph3)*D2R -AA51 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl2( 1, 8)= +AA51 *2.0*c*COS(ph2-ph4) +AA52 *2.0*g*COS(ph4-ph6)
      DDevfl2( 1, 9)= +AA51 *2.0*c*e*SIN(ph2-ph4)*D2R -AA52 *2.0*e*g*SIN(ph4-ph6)*D2R
      DDevfl2( 1,10)= +AA51 *2.0*d*COS(ph3-ph5)
      DDevfl2( 1,11)= +AA51 *2.0*d*f*SIN(ph3-ph5)*D2R
      DDevfl2( 1,12)= +AA52 *2.0*e*COS(ph4-ph6)
      DDevfl2( 1,13)= +AA52 *2.0*e*g*SIN(ph4-ph6)*D2R
      Devfl2( 0)= +AA53 *a*a +AA53 *b*b +AA54 *e*e +AA54 *f*f +AA55 *g*g
      DDevfl2( 0, 1)= +AA53 *2.0*a
      DDevfl2( 0, 2)= +AA53 *2.0*b
      DDevfl2( 0, 8)= +AA54 *2.0*e
      DDevfl2( 0,10)= +AA54 *2.0*f
      DDevfl2( 0,12)= +AA55 *2.0*g
      Devfl2(-1)= +AA56 *2.0*a*d*COS(   -ph3) +AA50 *2.0*b*c*COS(ph1-ph2) +AA57 *2.0*c*f*COS(ph2-ph5) +&
                  AA51 *2.0*d*e*COS(ph3-ph4) +AA52 *2.0*f*g*COS(ph5-ph6)
      DDevfl2(-1, 1)= +AA56 *2.0*d*COS(   -ph3)
      DDevfl2(-1, 2)= +AA50 *2.0*c*COS(ph1-ph2)
      DDevfl2(-1, 3)= -AA50 *2.0*b*c*SIN(ph1-ph2)*D2R
      DDevfl2(-1, 4)= +AA50 *2.0*b*COS(ph1-ph2) +AA57 *2.0*f*COS(ph2-ph5)
      DDevfl2(-1, 5)= +AA50 *2.0*b*c*SIN(ph1-ph2)*D2R -AA57 *2.0*c*f*SIN(ph2-ph5)*D2R
      DDevfl2(-1, 6)= +AA56 *2.0*a*COS(   -ph3) +AA51 *2.0*e*COS(ph3-ph4)
      DDevfl2(-1, 7)= +AA56 *2.0*a*d*SIN(   -ph3)*D2R -AA51 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl2(-1, 8)= +AA51 *2.0*d*COS(ph3-ph4)
      DDevfl2(-1, 9)= +AA51 *2.0*d*e*SIN(ph3-ph4)*D2R
      DDevfl2(-1,10)= +AA57 *2.0*c*COS(ph2-ph5) +AA52 *2.0*g*COS(ph5-ph6)
      DDevfl2(-1,11)= +AA57 *2.0*c*f*SIN(ph2-ph5)*D2R -AA52 *2.0*f*g*SIN(ph5-ph6)*D2R
      DDevfl2(-1,12)= +AA52 *2.0*f*COS(ph5-ph6)
      DDevfl2(-1,13)= +AA52 *2.0*f*g*SIN(ph5-ph6)*D2R
      Devfl2(-2)= +AA58 *2.0*a*f*COS(   -ph5) +AA48 *2.0*b*e*COS(ph1-ph4) +AA49 *2.0*d*g*COS(ph3-ph6) +AA46 *2.0*e*f*COS(ph4-ph5)
      DDevfl2(-2, 1)= +AA58 *2.0*f*COS(   -ph5)
      DDevfl2(-2, 2)= +AA48 *2.0*e*COS(ph1-ph4)
      DDevfl2(-2, 3)= -AA48 *2.0*b*e*SIN(ph1-ph4)*D2R
      DDevfl2(-2, 6)= +AA49 *2.0*g*COS(ph3-ph6)
      DDevfl2(-2, 7)= -AA49 *2.0*d*g*SIN(ph3-ph6)*D2R
      DDevfl2(-2, 8)= +AA48 *2.0*b*COS(ph1-ph4) +AA46 *2.0*f*COS(ph4-ph5)
      DDevfl2(-2, 9)= +AA48 *2.0*b*e*SIN(ph1-ph4)*D2R -AA46 *2.0*e*f*SIN(ph4-ph5)*D2R
      DDevfl2(-2,10)= +AA58 *2.0*a*COS(   -ph5) +AA46 *2.0*e*COS(ph4-ph5)
      DDevfl2(-2,11)= +AA58 *2.0*a*f*SIN(   -ph5)*D2R +AA46 *2.0*e*f*SIN(ph4-ph5)*D2R
      DDevfl2(-2,12)= +AA49 *2.0*d*COS(ph3-ph6)
      DDevfl2(-2,13)= +AA49 *2.0*d*g*SIN(ph3-ph6)*D2R


      Devfl0( 0)= +AA59 *a*a +AA59 *b*b +AA59 *c*c +AA59 *d*d +AA59 *e*e +AA59 *f*f +AA59 *g*g
      DDevfl0( 0, 1)= +AA59 *2.0*a
      DDevfl0( 0, 2)= +AA59 *2.0*b
      DDevfl0( 0, 4)= +AA59 *2.0*c
      DDevfl0( 0, 6)= +AA59 *2.0*d
      DDevfl0( 0, 8)= +AA59 *2.0*e
      DDevfl0( 0,10)= +AA59 *2.0*f
      DDevfl0( 0,12)= +AA59 *2.0*g
       !
       !   The preceding lines of code have been automatically generated
       !

       if ((0<= Rad_Model).and.(Rad_Model<= 3)) then
          call IRJ_ITB(n_pat,-ptr(ipiof,2,n_pat),sn,IRJ,DIRJ)
       else if (Rad_Model < 0) then
          call IRJ_Slater(-Rad_Model, 6, Zeff, hm, IRJ, DIRJ)
       else
          call IRJ_Hydrog(Rad_Model,3,6,Zeff, hm, IRJ, DIRJ)
       end if

       mffr=0.0
       mffi=0.0

       if(PRESENT(mdxir)) then
         mdxir(ORBPAR:ORBPAR+12)= 0
         mdxir(ZPAR)= 0
       end if

       ! Calculation for l=6
       l=6
       summ=  0.0
       summd= 0.0
       do m=-l,l
         p= sign(1,m)
         YY= Real_Spher_Harm_Ucvec(l,ABS(m),p,huc)
         summ= summ+ YY*Devfl6(m)
         summd= summd+ YY*DDevfl6(m,:)
       end do

       mffr= mffr- IRJ(l/2)*summ               !   The minus sign comes from the term i**l with l=6.

       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+12)= mdxir(ORBPAR:ORBPAR+12)-IRJ(l/2)*summd  !   The minus sign comes from the term i**l with l=6.
          mdxir(ZPAR)= mdxir(ZPAR)- DIRJ(l/2)*summ !   The minus sign comes from the term i**l with l=6.
       end if
       ! Calculation for l=4
       l=4

       summ= 0.0
       summd= 0.0
       do m=-l,l
         p= sign(1,m)
         YY= Real_Spher_Harm_Ucvec(l,ABS(m),p,huc)
         summ= summ+ YY*Devfl4(m)
         summd= summd+ YY*DDevfl4(m,:)
       end do

       mffr= mffr+ IRJ(l/2)*summ

       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+12)= mdxir(ORBPAR:ORBPAR+12)+IRJ(l/2)*summd
          mdxir(ZPAR)= mdxir(ZPAR)+ DIRJ(l/2)*summ
       end if

       ! Calculation for l=2
       l=2

       summ= 0.0
       summd= 0.0
       do m=-l,l
         p= sign(1,m)
         YY= Real_Spher_Harm_Ucvec(l,ABS(m),p,huc)
         summ= summ+ YY*Devfl2(m)
         summd= summd+ YY*DDevfl2(m,:)
       end do

       mffr= mffr- IRJ(l/2)*summ               !   The minus sign comes from the term i**l with l=2.

       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+12)= mdxir(ORBPAR:ORBPAR+12)-IRJ(l/2)*summd !   The minus sign comes from the term i**l with l=2.
          mdxir(ZPAR)= mdxir(ZPAR)- DIRJ(l/2)*summ !   The minus sign comes from the term i**l with l=2.
       end if

       ! Calculation for l=0
       l=0

       YY= 1.0/SQRT(4*PI)
       summ = YY*Devfl0(0)
       mffr= mffr+ IRJ(l/2)*summ
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+12)= mdxir(ORBPAR:ORBPAR+12)+ DIRJ(l/2)*YY*DDevfl0(0,1:13)
          mdxir(ZPAR)= mdxir(ZPAR) +DIRJ(l/2)*summ
       end if

       ! At this point the true mffr and mffi are 4*PI*MAG times the calculated
       ! mffr and mffi, and the derivatives are calculated with respect to normalized
       ! parameters but no with respect to input parameters
       ! We now correct the derivatives for normalization, scale (MAG), and 4*PI factor.
       ! The derivative with respect to the scale is also calculated.

       if (PRESENT(mdxir)) then
          mdxir(MAGPAR)= mffr*4.0*PI
       end if
       mffr= mffr*4.0*PI*MAG
       if (PRESENT(mdxir)) then
          mdxir(ORBPAR:ORBPAR+12)= mdxir(ORBPAR:ORBPAR+12) *4.0*PI*MAG
          mdxir(ZPAR)= mdxir(ZPAR) *4.0*PI* MAG
          deraux= matmul(normcorr,(/mdxir(ORBPAR), mdxir(ORBPAR+1:ORBPAR+11:2)/))
          mdxir(ORBPAR)= deraux(1)
          mdxir(ORBPAR+1:ORBPAR+11:2)= deraux(2:7)
       end if

       return
    End Subroutine Hxap_Ffactor

    Subroutine Mult_Ffactor(ipiof,sn,huc,mffr,mffi,mdxir,mdxii)
       integer,                              intent(in)   :: ipiof
       real(kind=cp),                        intent(in)   :: sn
       real(kind=cp), dimension(3),          intent(in)   :: huc
       real(kind=cp),                        intent(out)  :: mffr
       real(kind=cp),                        intent(out)  :: mffi
       real(kind=cp), dimension(:), optional,intent(out)  :: mdxir
       real(kind=cp), dimension(:), optional,intent(out)  :: mdxii
       !-----------------------------------------------
       !   L o c a l   P a r a m e t e r s
       !-----------------------------------------------

       real(kind=cp), parameter :: a0= 0.529177249 !Bohr radius in angstroms
       integer, parameter :: ORBPARS1= 19 ! l=0 m=0
       integer, parameter :: CHIPARS1= 20
       integer, parameter :: ORBPARS2= 21 ! l=0 m=0
       integer, parameter :: CHIPARS2= 22
       integer, parameter ::  ORBPARP= 23 ! l=1 m=-1,0,1 (23, 24, 25)
       integer, parameter ::  CHIPARP= 26
       integer, parameter ::  ORBPARD= 27 ! l=2 m=-2,-1,0,1,2 (27, 28, 29, 30, 31)
       integer, parameter ::  CHIPARD= 32
       integer, parameter ::  ORBPARF= 33 ! l=3 m=-3,-2,-1,0,1,2,3 (33, 34, 35, 36, 37, 38, 39)
       integer, parameter ::  CHIPARF= 40
       integer, parameter ::  ORBPARG= 41 ! l=4 m=-4,...,4 (41, 42, 43, 44, 45, 46, 47, 48, 49)
       integer, parameter ::  CHIPARG= 50

       !-----------------------------------------------
       !   L o c a l   V a r i a b l e s
       !-----------------------------------------------
       real(kind=cp) :: hm     ! modulus of q=2*pi*|h|
       integer :: l, m, p
       integer :: n, ns
       real(kind=cp) :: chis1, chis2, chip, chid, chif, chig,chi
       real(kind=cp) :: ISB, DISB, YY, summ
       real(kind=cp), dimension( 0:0) :: Devl0_1, Devl0_2
       real(kind=cp), dimension(-1:1) :: Devl1
       real(kind=cp), dimension(-2:2) :: Devl2
       real(kind=cp), dimension(-3:3) :: Devl3
       real(kind=cp), dimension(-4:4) :: Devl4
       chig=  xl(ipiof, CHIPARG)/a0
       devl4= xl(ipiof, ORBPARG:ORBPARG+8)
       chif=  xl(ipiof, CHIPARF)/a0
       devl3= xl(ipiof, ORBPARF:ORBPARF+6)
       chid=  xl(ipiof, CHIPARD)/a0
       devl2= xl(ipiof, ORBPARD:ORBPARD+4)
       chip=  xl(ipiof, CHIPARP)/a0
       devl1= xl(ipiof, ORBPARP:ORBPARP+2)
       chis2=   xl(ipiof, CHIPARS2)/a0
       devl0_2= xl(ipiof, ORBPARS2)
       chis1=   xl(ipiof, CHIPARS1)/a0
       devl0_1= xl(ipiof, ORBPARS1)
       ns=sp_pointers(ipiof)

       mffr= 0.0
       mffi= 0.0

       if (PRESENT(mdxir)) then
          mdxir(19:50)= 0.0
       end if

       if (PRESENT(mdxii)) then
          mdxii(19:50)= 0.0
       end if

       hm=4.0*pi*sqrt(sn)  !(Q in angstroms -1), sn=(sin(th)/lambda)^2

       if (chig > 0) then             ! l=4 terms (i)**l =1
         l=4
         chi= chig
         n= den_nlm(ns)%nl(6) !ng : exponent of the Slater part
         ISB = chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n,l,chi,hm)
         if (present(mdxir)) DISB= ISB/chi*(n+3)- chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n+1,l,chi,hm)
         summ= 0.0
         do m=-l,l
            p=sign(1,m)
            YY= Real_Spher_HarmCharge_Ucvec(l,abs(m),p,huc)
            summ= summ +YY* Devl4(m)
            if (present(mdxir)) then
               mdxir(ORBPARG+l+m)= mdxir(ORBPARG+l+m)+ ISB*YY
            end if
         end do
         mffr= mffr+ ISB*summ
         if (present(mdxir)) then
            mdxir(CHIPARG)= mdxir(CHIPARG)+ DISB*summ/a0
         end if
       end if

       if (chif > 0) then             ! l=3 terms (i)**l = -i
         l=3
         chi= chif
         n= den_nlm(ns)%nl(5) !nf
         ISB = chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n,l,chi,hm)
         if (present(mdxii)) DISB= ISB/chi*(n+3)- chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n+1,l,chi,hm)
         summ= 0.0
         do m=-l,l
            p=sign(1,m)
            YY= Real_Spher_HarmCharge_Ucvec(l,abs(m),p,huc)
            ! write (i_out,'("Real_Spher_Harm*ISB: ",2(i4),f14.7)') l,m,YY*ISB
            summ= summ +YY* Devl3(m)
            if (present(mdxii)) then
               mdxii(ORBPARF+l+m)= mdxii(ORBPARF+l+m)- ISB*YY   ! (i)**l = -i
            end if
         end do
         mffi= mffi- ISB*summ                                   ! (i)**l = -i
         if (present(mdxii)) then
            mdxii(CHIPARF)= mdxii(CHIPARF)- DISB*summ/a0        ! (i)**l = -i
         end if
       end if

       if (chid > 0) then             ! l=2 terms (i)**l =-1
         l=2
         chi= chid
         n= den_nlm(ns)%nl(4) !nd
         ISB = chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n,l,chi,hm)
         if (present(mdxir)) DISB= ISB/chi*(n+3)- chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n+1,l,chi,hm)
         summ= 0.0
         do m=-l,l
            p=sign(1,m)
            YY= Real_Spher_HarmCharge_Ucvec(l,abs(m),p,huc)
            summ= summ +YY* Devl2(m)
            if (present(mdxir)) then
               mdxir(ORBPARD+l+m)= mdxir(ORBPARD+l+m)- ISB*YY    ! (i)**l = -1
            end if
         end do
         mffr= mffr- ISB*summ                                    ! (i)**l = -1
         if (present(mdxir)) then
            mdxir(CHIPARD)= mdxir(CHIPARD)- DISB*summ/a0         ! (i)**l = -1
         end if
       end if

       if (chip > 0) then             ! l=1 terms (i)**l =  i
         l=1
         chi= chip
         n= den_nlm(ns)%nl(3) !np
         ISB = chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n,l,chi,hm)
         if (present(mdxii)) DISB= ISB/chi*(n+3)- chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n+1,l,chi,hm)
         summ= 0.0
         do m=-l,l
            p=sign(1,m)
            YY= Real_Spher_HarmCharge_Ucvec(l,abs(m),p,huc)
            summ= summ +YY* Devl1(m)
            if (present(mdxii)) then
               mdxii(ORBPARP+l+m)= mdxii(ORBPARP+l+m)+ ISB*YY    ! (i)**l = i
            end if
         end do
         mffi= mffi+ ISB*summ
         if (present(mdxii)) then
            mdxii(CHIPARP)= mdxii(CHIPARP)+ DISB*summ/a0         ! (i)**l = i
         end if
       end if

       if (chis2 > 0) then             ! l=0 terms (i)**l =1
         l=0
         chi= chis2
         n=den_nlm(ns)%nl(2) !ns2
         ISB = chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n,l,chi,hm)
         if (present(mdxir)) DISB= ISB/chi*(n+3)- chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n+1,l,chi,hm)
         summ= 0.0
         do m=-l,l
            p=sign(1,m)
            YY= Real_Spher_HarmCharge_Ucvec(l,abs(m),p,huc)
            summ= summ +YY* Devl0_2(m)
            if (present(mdxir)) then
               mdxir(ORBPARS2+l+m)= mdxir(ORBPARS2+l+m)+ ISB*YY
            end if
         end do
         mffr= mffr+ ISB*summ
         if (present(mdxir)) then
            mdxir(CHIPARS2)= mdxir(CHIPARS2)+ DISB*summ/a0
         end if
       end if

      if (chis1 > 0) then             ! l=0 terms (i)**l =1
         l=0
         chi= chis1
         n=den_nlm(ns)%nl(1)      !ns1
         ISB = chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n,l,chi,hm)
         if (present(mdxir)) DISB= ISB/chi*(n+3)- chi**(n+3)/factorial(n+2)*Int_Slater_Bessel(n+1,l,chi,hm)
         summ= 0.0
         do m=-l,l
            p=sign(1,m)
            YY= Real_Spher_HarmCharge_Ucvec(l,abs(m),p,huc)
            summ= summ +YY* Devl0_1(m)
            if (present(mdxir)) then
               mdxir(ORBPARS1+l+m)= mdxir(ORBPARS1+l+m)+ ISB*YY
            end if
         end do
         mffr= mffr+ ISB*summ
         if (present(mdxir)) then
            mdxir(CHIPARS1)= mdxir(CHIPARS1)+ DISB*summ/a0
         end if
       end if
       mffr= 4.0*PI*mffr
       mffi= 4.0*PI*mffi

       if (present(mdxir)) then
          mdxir(19:50)= 4.0*PI*mdxir(19:50)
       end if

       if (present(mdxii)) then
          mdxii(19:50)= 4.0*PI*mdxii(19:50)
       end if
       return
    End Subroutine Mult_Ffactor

    Subroutine IRJ_Hydrog(n,l,ll,Zeff,s,IRJ,DIRJ)
        ! Returns the integrals \int_0^{\infty} R_{n,l}^2(Zeff;r) j_{k}(r s) r^2 dr   for k=0, 2, 4, ... ll
        ! and their derivatives with respect to Zeff (IRJ and DIRJ vectors respectively).
        ! R_{n,l}(Zeff;r) is an hidrogen-like orbital with effective atomic number Zeff. n >= l+1
        integer, intent(in)                         :: n, l, ll
        real(kind=cp), intent (in)                  :: Zeff, s
        real(kind=cp), dimension (0:), intent (out) :: IRJ, DIRJ
        !-------------------------------------
        !   L o c a l   V a r i a b l e s
        !-------------------------------------
        real(kind=cp), parameter :: a0= 0.529177249
        real(kind=cp) :: chi
        integer :: i, j, k
        real(kind=cp) :: A, B
        chi= 2.0*Zeff/real(n)/a0

        IRJ= 0.0
        DIRJ= 0.0

        if (n < l+1 ) then
           ! Warning message
           err_flipr=.true.
           err_flipr_mess=" Invalid n in hydrogen-like radial function"
           return
        end if
        if (chi <= 0.0) then
           err_flipr=.true.
           err_flipr_mess="Zeff smaller or equal zero in hydrogen-like radial function"
           return
        end if
        A=factorial(n+l)* factorial(n-l-1)/ 2.0/ n
        do i=0,2*(n-l-1)
           B=0.0
           do j=max(0,i-(n-l-1)),min(i,(n-l-1))
              B= B+ 1.0/(1.0*factorial(i-j)*factorial(j)*factorial(i-j+2*l+1)* factorial(j+2*l+1)* &
                         factorial(n-l-1-i+j)* factorial(n-l-1-j))
           end do
           B=B*(-2)**i
           do k=0, ll/2
              IRJ(k)= IRJ(k)+ B*chi**(2*l+3+i)*Int_Slater_Bessel(2*l+i,2*k,chi,s)
              DIRJ(k)= DIRJ(k)+B*((2*l+3+i)*chi**(2*l+2+i)*Int_Slater_Bessel(2*l+i,2*k,chi,s)- &
                               chi**(2*l+3+i)*Int_Slater_Bessel(2*l+1+i,2*k,chi,s))
           end do
        end do
           IRJ= IRJ*A
           DIRJ= DIRJ*A*2.0/n/a0
    End Subroutine IRJ_Hydrog

    Subroutine IRJ_ITB(n_pat,ni,sn,IRJ,DIRJ)
       ! Returns the integrals \int_0^{\infty} R_{n,l}^2(r) j_{k}(r s) r^2 dr for k=0, 2, 4, 6
       ! approached to the expression given in the International Tables of Crystallography Vol. B.
       INTEGER, INTENT(IN)                        :: n_pat
       INTEGER, INTENT(IN)                        :: ni
       real(kind=cp), INTENT(IN)                  :: sn
       real(kind=cp), INTENT(OUT), dimension(0:)  :: IRJ
       real(kind=cp), INTENT(OUT), dimension(0:)  :: DIRJ
       ! Local variables
       integer :: i,j
       DIRJ= 0.0
       IRJ= 0.0
       do j=0,min(size(IRJ),3)
          IRJ(j)= ac((j+1)*7,ni,n_pat)                ! j=l/2
!          j0=ac(7,ni,n_pat)
!          j2=ac(14,ni,n_pat)
!          j4=ac(21,ni,n_pat)
!          j6=ac(28,ni,n_pat)
          DO i=1,5,2
            IRJ(j)= IRJ(j)+ ac(j*7+i,ni,n_pat)*EXP(-ac(7*j+i+1,ni,n_pat)*sn)
!            j0=j0+ac(i,ni,n_pat)*EXP(-ac(i+1,ni,n_pat)*sn)      !<j0> value for Q=H+k
!            j2=j2+ac(i+7,ni,n_pat)*EXP(-ac(i+8,ni,n_pat)*sn)    !<j2> value for Q=H+k
!            j4=j4+ac(i+14,ni,n_pat)*EXP(-ac(i+15,ni,n_pat)*sn)  !<j4> value for Q=H+k
!            j6=j6+ac(i+21,ni,n_pat)*EXP(-ac(i+22,ni,n_pat)*sn)  !<j6> value for Q=H+k
          END DO
       end do
       IRJ(1:)=IRJ(1:)*sn
!      IRJ(0)= IRJ(0)/sn
       return
    End Subroutine IRJ_ITB

    Subroutine IRJ_SP(n,Zeffp,Zeffs,s,IRJ,DIRJ)
        ! Returns the integral  \int_0^{\infty} R_{n,0}(Zeff_s;r)R_{n,1}(Zeff_p;r) j_{1}(r s) r^2 dr
        ! and their derivatives with respect to Zeff_p and Zeff_s (IRJ(0) and DIRJ(0:1)respectively).
        ! R_{n,l}(Zeff;r) is an hidrogen-like orbital with effective atomic number Zeff. n >= l+1
        integer,                       intent(in)   :: n
        real(kind=cp),                 intent (in)  :: Zeffp, Zeffs, s
        real(kind=cp), dimension (0:), intent (out) :: IRJ, DIRJ
        !-------------------------------------
        !   L o c a l   V a r i a b l e s
        !-------------------------------------
        real(kind=cp), parameter :: a0= 0.529177249
        real(kind=cp) :: chis, chip, ISB0, ISB1, A, B, sum0, sump, sums
        integer :: i, j
        chis= Zeffs/real(n)/a0
        chip= Zeffp/real(n)/a0

        IRJ= 0.0
        DIRJ= 0.0

        if (n < 2 ) then
           ! Warning message
           err_flipr=.true.
           err_flipr_mess=" Invalid n in hydrogen-like radial function "
           return
        end if
        if ((chis <= 0.0).or.(chip <= 0.0)) then
           err_flipr=.true.
           err_flipr_mess="Zeff smaller or equal zero in hydrogen-like radial function"
           return
        end if
        A=sqrt( 1.0*factorial(n+1)* factorial(n)* factorial(n-1)*factorial(n-2))/ n *8.0
        do i=0,2*n-3
           ISB0= Int_Slater_Bessel(i+1,1,chis+chip,s)
           ISB1= Int_Slater_Bessel(i+2,1,chis+chip,s)
           sum0= 0.0
           sums= 0.0
           sump= 0.0
           do j=max(0,i-(n-1)),min(i,n-2)
              B= 1.0*(-2)**i/(1.0*factorial(i-j)*factorial(i-j+1)*factorial(n-i+j-1)*factorial(j)*factorial(j+3)*factorial(n-2-j))
              sum0= sum0 +B* chis**(i-j+1.5)* chip**(j+2.5)
              sums= sums +B* (i-j+1.5)*chis**(i-j+0.5)* chip**(j+2.5)
              sump= sump +B* chis**(i-j+1.5)*(j+2.5)* chip**(j+1.5)
           end do
           IRJ(0)= IRJ(0)+ sum0*ISB0
           DIRJ(0)= DIRJ(0)+ sump*ISB0 -sum0*ISB1
           DIRJ(1)= DIRJ(1)+ sums*ISB0 -sum0*ISB1
        end do
        IRJ= IRJ*A
        DIRJ= DIRJ*A/n/a0

    End Subroutine IRJ_SP

    Subroutine IRJ_Slater(n,ll,chi,s,IRJ,DIRJ)
        ! Returns the integrals
        ! \int_0^{\infty} [chi**(2n+3)/(2n+2)! r**2n exp(-2*chi*r)]^2 j_{k}(r s) r^2 dr
        ! for k=0, 2, 4, ... ll
        ! and their derivatives with respect to chi (IRJ and DIRJ vectors respectively).
        integer,                       intent(in)  :: n, ll
        real(kind=cp),                 intent (in) :: chi, s
        real(kind=cp), dimension (0:), intent (out):: IRJ, DIRJ
        !-------------------------------------
        !   L o c a l   V a r i a b l e s
        !-------------------------------------
        integer :: i

        IRJ= 0.0
        DIRJ= 0.0

        do i=0,ll/2
           IRJ(i)= (2*chi)**(2*n+3)/factorial(2*n+2)*Int_Slater_Bessel(2*n,2*i,2*chi,s)
           DIRJ(i)= (2*n+3)/chi* IRJ(i) -&
                 2*(2*chi)**(2*n+3)/factorial(2*n+2)*Int_Slater_Bessel(2*n+1,2*i,2*chi,s)
        end do

    End Subroutine IRJ_Slater



end module CFML_FlipR_Mod
