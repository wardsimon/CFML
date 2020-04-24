SubModule (CFML_Profiles) PRF_031
   Contains

    !!----
    !!---- FUNCTION TCH_PVOIGT
    !!----
    !!----
    !!---- Update: April - 2009
    !!
    Pure Module Function TCH_pVoigt(X,Par) Result (Pv_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: pv_val

       !--- Local variables ---!
       real(kind=cp)            :: Hg,Hl,eta,H,x2,ag,bg,al,bl,lor,gauss,r
       real(kind=cp), parameter :: o1= 2.69269, o2=2.42843, o3=4.47163, o4= 0.07842
       real(kind=cp), parameter :: e1= 1.36603, e2=0.47719, e3=0.11116

       Hg=par(1)
       Hl=par(2)

       !> Calculate H and eta from Tomson-Cox-Hasting formula
       H=hg**5+o1*hg**4*hl+o2*hg**3*hl**2+o3*hg**2*hl**3+o4*hg*hl**4+hl**5
       H=abs(H)**0.2_cp
       r = hl/H                       !HL/FWHM
       eta = max( 1.0e-06_cp, r*(e1 -(e2 + e3*r)*r) )  !eta
       x2=x*x
       ag= 0.93943727869965133377234032841018_cp/H
       bg= 2.7725887222397812376689284858327_cp/(H*H)
       al= 0.63661977236758134307553505349006_cp/H
       bl= 4.0_cp/(H*H)
       gauss = ag* exp(-bg*x2)
       lor   = al/(1.0_cp+bl*x2)
       pv_val = eta*lor + (1.0_cp - eta)*gauss

       return
    End Function TCH_pVoigt

    !!----
    !!---- SUBROUTINE TCH_PVOIGT_DER
    !!----
    !!----  Pv_Val is the value of the function for the argument x
    !!----  Par=(/HG,HL/)  and Dpar(1:3)=(/derx,derHG,derHL/)
    !!----
    !!---- Update: April - 2009
    !!
    Pure Module Subroutine TCH_pVoigt_Der(X,Par,Pv_Val,dPar)
       !---- Arguments ----!
       real(kind=cp),                       intent(in) :: x
       real(kind=cp), dimension(:),         intent(in) :: par
       real(kind=cp),                       intent(out):: pv_val
       real(kind=cp), optional,dimension(:),intent(out):: dpar

       !--- Local variables ---!
       real(kind=cp), parameter :: o1= 2.69269_cp, o2=2.42843_cp, o3=4.47163_cp, o4= 0.07842_cp
       real(kind=cp), parameter :: e1= 1.36603_cp, e2=0.47719_cp, e3=0.11116_cp
       real(kind=cp) :: Hg,Hl,eta,H,x2,ag,bg,al,bl,lor,gauss, invH,invH2,r
       real(kind=cp) :: derH,derEta,derx,dlorH,dgaussH,lorp,gaussp, &
                        dhdhg,dhdhl,deta,detag,detal,derHg,derHl

       Hg=par(1)
       Hl=par(2)
       !Calculate H and eta from Tomson-Cox-Hasting formula
       H=hg**5+o1*hg**4*hl+o2*hg**3*hl**2+o3*hg**2*hl**3+o4*hg*hl**4+hl**5
       H=abs(H)**0.2_cp
       r = hl/H                       !HL/FWHM
       eta = max( 1.0e-06_cp, r*(e1 -(e2 + e3*r)*r) )  !eta
       x2=x*x
       invH=1.0_cp/H
       invH2=invH*invH
       ag= 0.93943727869965133377234032841018_cp*invH
       bg= 2.7725887222397812376689284858327_cp*invH2
       al= 0.63661977236758134307553505349006_cp*invH
       bl= 4.0_cp*invH2
       gauss = ag* exp(-bg*x2)
       lor   = al/(1.0_cp+bl*x2)
       pv_val = eta*lor + (1.0_cp - eta)*gauss

       if(present(dpar)) then
          dhdhg = 0.2_cp/H**4*(5.0_cp*hg**4+ 4.0_cp*o1* hg**3*hl+  &
                  3.0_cp*o2*hg*hg*hl*hl + 2.0_cp*o3*hg*hl**3 + o4*hl**4)
          dhdhl = 0.2_cp/H**4*(o1*hg**4+ 2.0_cp*o2*hg**3*hl+  &
                  3.0_cp*o3*hg*hg*hl*hl + 4.0_cp*o4*hg*hl**3 + 5.0_cp*hl**4)
           deta = e1-(2.0_cp*e2 - 3.0_cp*e3*r)*r  !derivative of ETA w.r.t. r
          detag = -r*deta*dhdhg*invH
          detal = (1.0_cp-r*dhdhl)*deta*invH

          derEta= lor-gauss  !Eta
          lorp = -2.0_cp *lor*lor*bl*x/al  !x
          gaussp = -2.0_cp * gauss * bg * x  !x
          derx=eta*lorp+(1.0-eta)*gaussp  !x

          dlorH= (2.0_cp*bl*lor*x2/al -1.0_cp)*invH*lor
          dgaussH= (2.0_cp*bg*x2-1.0_cp)*invH*gauss
          derH=eta*dlorH + (1.0_cp-eta) * dgaussH

          derHG= derH * dhdhg + derEta * detag  !Chain rule
          derHL= derH * dhdhl + derEta * detal
          dpar(1:3)=(/derx,derHG,derHL/)
       end if

       return
    End Subroutine TCH_pVoigt_Der

End SubModule PRF_031
