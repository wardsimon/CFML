Program MagRef
 use CFML_GlobalDeps
 use CFML_Math_3D,                  only: Veclength
 use CFML_crystallographic_symmetry,only: space_group_type, Write_SpaceGroup, get_stabilizer
 use CFML_Atom_TypeDef,             only: Atom_List_Type, Write_Atom_List, MAtom_list_Type, &
                                          Get_Atom_2nd_Tensor_Ctr
 use CFML_crystal_metrics,          only: Crystal_Cell_Type, Write_Crystal_Cell
 use CFML_Reflections_Utilities,    only: Hkl_s
 use CFML_IO_Formats,               only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
 use CFML_Propagation_vectors,      only: K_Equiv_Minus_K
 use CFML_Magnetic_Symmetry
 use CFML_Magnetic_Structure_Factors
 use CFML_Geometry_SXTAL,           only: z1frnb

 implicit none

 type (file_list_type)       :: fich_cfl
 type (Space_Group_Type)     :: SpG
 type (MagSymm_k_Type)       :: MGp
 type ( Atom_list_Type)      :: A
 type (MAtom_list_Type)      :: Am
 type (Crystal_Cell_Type)    :: Cell
 type (MagH_Type)            :: Mh

 character(len=256)           :: line,filcod,fileflip     !Name of the input file
 character(len=1)             :: sig
 real(kind=cp)                :: sn,sf2,lambda,polarup,polardown
 real(kind=cp), dimension(3)  :: u_vect
 integer                      :: Num, lun=1, ier,i,j,m,ih,ik,il,iv, n_ini, n_end, fliptyp
 complex(kind=cp)             :: fc
 real(kind=cp), dimension(6)  :: extc
 real(kind=cp), dimension(3,3):: ub

 integer                     :: narg,i_flip,iext,ity
 Logical                     :: esta, arggiven=.false.,flipp=.false.,ext=.false.,ubgiven=.false.
      !---- Arguments on the command line ----!
      narg=COMMAND_ARGUMENT_COUNT()

      if(narg > 0) then
              call GET_COMMAND_ARGUMENT(1,filcod)
              arggiven=.true.
      end if

      write(unit=*,fmt="(/,/,7(a,/))")                                                 &
           "              ------ P r o g r a m     M a g R e f  ------"               , &
           "                    ---- Version 0.3 June-2014 ----"                    , &
           "    ******************************************************************"  , &
           "    * Calculates magnetic structure factors, and magnetic interaction*"  , &
           "    * vectors from magnetic structures by reading a *.CFL file       *"  , &
           "    * Calculates also Induced Magnetic Structure Factor Tensors      *"  , &
           "    ******************************************************************"  , &
           "                            (JRC- June 2014)"
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
           "              ------ P r o g r a m     M a g R e f  ------"               , &
           "                    ---- Version 0.3 June-2014 ----"                    , &
           "    ******************************************************************"  , &
           "    * Calculates magnetic structure factors, and magnetic interaction*"  , &
           "    * vectors from magnetic structures by reading a *.CFL file       *"  , &
           "    * Calculates also Induced Magnetic Structure Factor Tensors      *"  , &
           "    ******************************************************************"  , &
           "                            (JRC- June 2014)"

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
     write(unit=*,fmt="(/,a)") "    =========================="
     write(unit=*,fmt="( a )") "    Text of the input CFL file"
     write(unit=*,fmt="(a,/)") "    =========================="
     do i=1,fich_cfl%nlines
       write(unit=lun,fmt="(a,i5,a)") " Line:",i,"  "//fich_cfl%line(i)
       write(unit=*,fmt="(a,i5,a)") " Line:",i,"  "//fich_cfl%line(i)
     end do


       call Write_Crystal_Cell(Cell,lun)
       call Write_SpaceGroup(SpG,lun,full=.true.)
       call Write_Atom_List(A,lun=lun)
       !Search for flipping ratio files in order to extract phased observed magnetic structure factors
       lambda=0.8
       fliptyp=1
       polarup=1.0 ; polardown=1.0
       do i=1,fich_cfl%nlines
         line=adjustl(fich_cfl%line(i))
         if(U_Case(line(1:4) == "FLIP") then
           fileflip=adjustl(line(5:))
           flipp=.true.
         end if
         if(U_Case(line(1:6) == "LAMBDA") then
           read(unit=line(7:),fmt=*,iostat=ier) lambda
           if(ier /= 0) lambda=0.8
         end if
         if(U_Case(line(1:7) == "POLARUD") then
           read(unit=line(8:),fmt=*,iostat=ier) polarup,polardown
           if(ier /= 0) then
             polarup=1.0 ; polardown=1.0
           end if
         end if
         if(U_Case(line(1:4) == "EXTI") then
           read(unit=line(5:),fmt=*,iostat=ier) extc, iext
           if(ier /= 0) then
              read(unit=line(5:),fmt=*,iostat=ier) extc(1)
              if(ier /= 0) extc(:) = 0.0
              iext=1
           end if
           if(any(extc > 0.0001)) ext=.true.
         end if
         if(flipp .and. ext) exit
       end do

       !If flipp is true, extract only observed magnetic structure factors
       if(flipp) then
         if(SpG%centred /= 2) then
           write(unit=*,fmt="(a)") " => The given space group has no inversion centre at the origin!"
           write(unit=*,fmt="(a)") " => Phased magnetic structure factors cannot be obtained!"
           stop
         end if
         open(newunit=i_flip,file=fileflip,status="old",action="read",iostat=ier,position="rewind")
         if(ier/=0) then
           write(unit=*,fmt="(a)") " => File: "//trim(fileflip)//" not found!"
           stop
         end if
         !Detect the type of file *.int or other
         if(index(fileflip,".int") /= 0) fliptyp=2
         if(fliptyp == 2) then

           read(unit=i_flip,fmt="(a)") line
           read(unit=i_flip,fmt="(a)") line
           read(unit=i_flip,fmt="(a)") line
           read(unit=line,fmt=*) lambda,ity, i, polarup,polardown,j
           if(ity /= 2) then
             write(unit=*,fmt="(a)") " => The file: "//trim(fileflip)//" doesn't contain flipping ratios!"
             stop
           end if
           if(j /= 0) then
             read(unit=i_flip,fmt=*,iostat=ier) ub(1,:)
             read(unit=i_flip,fmt=*,iostat=ier) ub(2,:)
             read(unit=i_flip,fmt=*,iostat=ier) ub(3,:)
             if(ier == 0) ubgiven=.true.
           end if
           !Now the position for reading is the good one for getting hkl,flip and sflip
         end if
         call readn_set_MsF()
         stop
       end if

       n_ini=1
       n_end=fich_cfl%nlines
       call Readn_Set_Magnetic_Structure(fich_cfl,n_ini,n_end,MGp,Am)
       if(err_MagSym) then
         write(unit=*,fmt="(a)") "   "//err_MagSym_Mess
         stop
       end if
       if(Am%suscept) then
          call Write_Magnetic_Structure(lun,MGp,Am,cell=Cell)
       else
          call Write_Magnetic_Structure(lun,MGp,Am)
       end if

      do
         if(Am%suscept) then
           write(unit=*,fmt="(a)",advance="no") &
           " => Enter a reflections as 3 integers -> (h,k,l) : "
           read(unit=*,fmt=*) ih,ik,il
           m=1
         else
           write(unit=*,fmt="(a)",advance="no") &
           " => Enter a magnetic reflection as 4 integers -> (h,k,l,m)=H+sign(m)*k(abs(m)): "
           read(unit=*,fmt=*) ih,ik,il,m
         end if
         if( m == 0 .or. abs(ih)+abs(ik)+abs(il) == 0 ) exit
         !construct partially the object Mh
         j=sign(1,m)
         sig="+"
         if( j == -1) sig="-"
         Mh%signp=real(-j)  ! sign "+" for H-k and "-" for H+k
         iv=abs(m)
         Mh%num_k=iv
         Mh%h= real((/ih,ik,il/)) - Mh%signp*MGp%kvec(:,iv)
         Mh%s = hkl_s(Mh%h,Cell)
         Mh%keqv_minus=K_Equiv_Minus_K(MGp%kvec(:,iv),MGp%latt)

         if(Am%suscept) then
            !write(unit=*,fmt="(a,i2,a)",advance="no") &
            !" => Enter the strength(in Tesla) and direction of applied magnetic field: "
            !read(unit=*,fmt=*) Am%MagField, Am%dir_MField
            call Calc_Induced_Sk(cell,SpG,Am%MagField,Am%dir_MField,Am,6)
            call Calc_Magnetic_Strf_Tensor(SpG,Am,Mh)
            write(unit=lun,fmt="(/,a,3i4,a)")  "  Reflection: (",ih,ik,il,") "
            write(unit=*,  fmt="(/,a,3i4,a)")  "  Reflection: (",ih,ik,il,") "
            write(unit=lun,fmt="(a)")          "  Real part of Tensorial Magnetic Structure Factor: "
            write(unit=*,  fmt="(a)")          "  Real part of Tensorial Magnetic Structure Factor: "
            do i=1,3
              write(unit=lun,fmt="(a,3f12.5)") "       ",real(Mh%TMsF(i,:))
              write(unit=*  ,fmt="(a,3f12.5)") "       ",real(Mh%TMsF(i,:))
            end do
            write(unit=lun,fmt="(a)")          "  Imaginary part of Tensorial Magnetic Structure Factor: "
            write(unit=*  ,fmt="(a)")          "  Imaginary part of Tensorial Magnetic Structure Factor: "
            do i=1,3
              write(unit=lun,fmt="(a,3f12.5)") "       ",aimag(Mh%TMsF(i,:))
              write(unit=*  ,fmt="(a,3f12.5)") "       ",aimag(Mh%TMsF(i,:))
            end do
            call Calc_Induced_MsF_MiV(cell,Am%MagField,Am%dir_MField,Mh)
            write(unit=lun,fmt="(a,2(3f8.4,a))") "  Magnetic Structure   Factor : (",real(Mh%MsF),")+i(",aimag(Mh%MsF),") "
            write(unit=lun,fmt="(a,2(3f8.4,a))") "  Magnetic Interaction Vector : (",real(Mh%MiV),")+i(",aimag(Mh%MiV),") "
            write(unit=*,  fmt="(a,2(3f8.4,a))") "  Magnetic Structure   Factor : (",real(Mh%MsF),")+i(",aimag(Mh%MsF),") "
            write(unit=*,  fmt="(a,2(3f8.4,a))") "  Magnetic Interaction Vector : (",real(Mh%MiV),")+i(",aimag(Mh%MiV),") "
            write(unit=lun,fmt="(a,f12.5 )")     "  Square of Mag. int.  Vector : ",Mh%sqMiV
            write(unit=*  ,fmt="(a,f12.5 )")     "  Square of Mag. int.  Vector : ",Mh%sqMiV

         else

            !Calculate magnetic structure factor and magnetic interaction vector
            call Calc_Magnetic_StrF_MiV(Cell,MGp,Am,Mh)
            write(unit=lun,fmt="(/,a,3i4,a,3f8.4,a)") "  Reflection: (",ih,ik,il,") "//sig//" (",MGp%kvec(:,iv),")"
            write(unit=lun,fmt="(a,3f8.4,a)")         "              (",Mh%h,")"
            write(unit=*,  fmt="(/,a,3i4,a,3f8.4,a)") "  Reflection: (",ih,ik,il,") "//sig//" (",MGp%kvec(:,iv),")"
            write(unit=*,  fmt="(a,3f8.4,a)")         "              (",Mh%h,")"
            write(unit=lun,fmt="(a,2(3f8.4,a))") "  Magnetic Structure   Factor : (",real(Mh%MsF),")+i(",aimag(Mh%MsF),") "
            write(unit=lun,fmt="(a,2(3f8.4,a))") "  Magnetic Interaction Vector : (",real(Mh%MiV),")+i(",aimag(Mh%MiV),") "
            write(unit=*,  fmt="(a,2(3f8.4,a))") "  Magnetic Structure   Factor : (",real(Mh%MsF),")+i(",aimag(Mh%MsF),") "
            write(unit=*,  fmt="(a,2(3f8.4,a))") "  Magnetic Interaction Vector : (",real(Mh%MiV),")+i(",aimag(Mh%MiV),") "
            write(unit=lun,fmt="(a,f12.5 )")     "  Square of Mag. int.  Vector : ",Mh%sqMiV
            write(unit=*  ,fmt="(a,f12.5 )")     "  Square of Mag. int.  Vector : ",Mh%sqMiV
         end if
       end do


       write(unit=*,fmt="(a)") " Normal End of program: MagRef "
       write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".cal"
     end if

     close(unit=lun)

    contains



    !!---- Subroutine Calc_Induced_MsF_MiV(cell,MField,dir_MField,Mh)
    !!----    !---- Arguments ----!
    !!----    type(Crystal_Cell_type),    intent(in)     :: Cell
    !!----    real(kind=cp),              intent(in)     :: MField
    !!----    real(kind=cp),dimension(3), intent(in)     :: dir_MField
    !!----    type(MagH_Type),            intent(in out) :: Mh
    !!----
    !!----  This subroutine completes the object Mh when the tensorial magnetic structure
    !!----  factor has been previously calculated.
    !!----
    !!----  Created: June 2014 (JRC)
    !!----
    Subroutine Calc_Induced_MsF_MiV(cell,MField,dir_MField,Mh)
       !---- Arguments ----!
       type(Crystal_Cell_type),    intent(in)     :: Cell
       real(kind=cp),              intent(in)     :: MField
       real(kind=cp),dimension(3), intent(in)     :: dir_MField
       type(MagH_Type),            intent(in out) :: Mh
       !--- Local variables ---!
       real(kind=cp)                  :: s
       real(kind=cp),    dimension(3) :: u_vec,er,ed
       complex(kind=cp), dimension(3) :: Mc, MiV

       !
       u_vect=MField * dir_MField / Veclength(Cell%Cr_Orth_cel,dir_MField)
       Mh%MsF=matmul(Mh%TMsF,u_vect)  !Magnetic structure factor from TMsF tensor
       !---- Calculation of the Magnetic Interaction vector ----!
       s  = 2.0*Mh%s            !1/d=r*, M = M// + Mp   => Mp = M - M// = M - (M.e)e
       er = Mh%h/s              !unitary vector referred to the reciprocal basis
       ed = matmul(cell%GR,er)  !  "        "       "             direct    "
       Mc  = Mh%MsF / Cell%cell                !Magnetic structure factor in basis {a,b,c}
       MiV = Mc - dot_product(er,Mc) * ed      !Magnetic interaction vector in basis {a,b,c}
       Mh%MiV  =  MiV * Cell%cell              !Magnetic interaction vector in basis {e1,e2,e3}
       Mh%MiVC = matmul(Cell%Cr_Orth_cel,MiV)  !Magnetic interaction vector in Cartesian components
       Mh%sqMiV= dot_product(Mh%MiVC, Mh%MiVC)
       return
    End Subroutine Calc_Induced_MsF_MiV

    Subroutine readn_set_MsF()
      use CFML_Structure_Factors,       only: Calc_hkl_StrFactor
      use CFML_Extinction_Corrections,  only: Correct_FlippingRatios
      real(kind=cp), dimension(3) :: h,hc
      real(kind=cp)     :: q, flipping_ratio, sigma, omeg,gam,nu,  &
                           sn, AN,BN,AM,BM,yp,ym,ypm, Valmub=0.2695 !0.5*(gamma_N * r0)
      character(len=30) :: fmm
      complex(kind=cp)  :: fc
       !Reading flipping ratio file

       do
          if(fliptyp == 1) then
            read(unit=i_flip,fmt="(i8,6f8.3,2f10.2)",iostat=ier)  numor,h(:), omeg, gam, nu,  flipping_ratio, sigma
          else
            read(unit=i_flip,fmt=*,iostat=ier)  h(:),flipping_ratio, sigma,i,q
          end if
          if(ier /= 0) exit
          !Calculation of q=sin^2(alpha) if needed
          if(ubgiven) then
              hc=matmul(ub,h)
              q = 1.0-hc(3)*hc(3)/dot_product(hc,hc)  !sin^2(alpha) = 1 - (h(3)/|h|)^2
          else
              if(fliptyp == 1) then
                call z1frnb(lambda,gam,omeg,nu,hc)
                q=1.0-hc(3)*hc(3)/dot_product(hc,hc)
              end if
          end if
          sn=hkl_s(h,cell)
          sn=sn*sn
          call Calc_hkl_StrFactor("S","N",h,sn,A,SpG,sf2,fc=fc)
          AN=real(fc)
          BN=aimag(fc)
          if(ext) then
            Call Correct_FlippingRatios(iext,Lambda,q,extc,sn,hkl,AN,BN,AM,BM,yp,ym,ypm,dyp,dym,dypm,dymag)
          end if
          !Now solve the equation to get Fm
          !Approximations:
          !ppp=1/2((1+p+)yp+(1-p+)ym)    ppm=1/2((1+p+)yp-(1-p+)ym)
          !mpp=1/2((1-p-)yp+(1+p-)ym)    mpm=1/2((1-p-)yp-(1+p-)ym)
          !General expression for flipping ratio
          !R=I+/I-
          !I+= (NN*+ MM*q2)ppp + 2(AN.AM + BN.BM)q.ppm + MM*ypm q(1-q)
          !I-= (NN*+ MM*q2)mpp + 2(AN.AM + BN.BM)q.mpm + MM*ypm q(1-q)
          ! Particular case for a centrosymmetric structure M*=M, N*=N, B=0 AN=N,AM=M
          !I+= (N2+ M2.q2)ppp + 2.N.M.q.ppm + M2 ypm q(1-q)
          !I-= (N2+ M2.q2)mpp + 2.N.M.q.mpm + M2 ypm q(1-q)
          ! Case of no extinction: yp=1,ym=1, ypm=1
          ! ppp=1, mpp=1, ppm=p+  mpm=-p-
          !I+= (N2+ M2.q2) + 2.N.M.q.p+ + M2 q(1-q)
          !I-= (N2+ M2.q2) - 2.N.M.q.p- + M2 q(1-q)


          !     (1+ g2.q2)ppp + 2.g.q.ppm + g2 ypm q(1-q)
          !R= ------------------------------------------------
          !     (1+ g2.q2)mpp + 2.g.q.mpm + g2 ypm q(1-q)
          !
          !
            sf2=sf2*ys

       end do
      write(unit=lun,fmt="(/,a,/)") "   H   K   L   Mult  SinTh/Lda    |Fc|       Phase        F-Real      F-Imag      Num"
      do i=1, hkl%nref
         sn=hkl%ref(i)%s * hkl%ref(i)%s
         call Calc_StrFactor("P","N",i,sn,A,Spg,sf2,fc=fc)
         write(unit=lun,fmt="(3i4,i5,5f12.5,i8,f12.5)") hkl%ref(i)%h, hkl%ref(i)%mult, &
              hkl%ref(i)%S, hkl%ref(i)%Fc, hkl%ref(i)%Phase, real(fc), aimag(fc), i, sqrt(sf2)
      end do
    End Subroutine readn_set_MsF

End Program MagRef

