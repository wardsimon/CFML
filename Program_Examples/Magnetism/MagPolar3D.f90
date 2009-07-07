Program MagPolar3D
 use CFML_crystallographic_symmetry,only: space_group_type, Write_SpaceGroup
 use CFML_Atom_TypeDef,             only: Atom_List_Type, Write_Atom_List, MAtom_list_Type
 use CFML_crystal_metrics,          only: Crystal_Cell_Type, Write_Crystal_Cell
 use CFML_Reflections_Utilities,    only: Hkl_s
 use CFML_IO_Formats,               only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
 use CFML_Propagation_vectors,      only: K_Equiv_Minus_K
 use CFML_Magnetic_Symmetry
 use CFML_Magnetic_Structure_Factors
 use CFML_Polarimetry

 implicit none

 type (file_list_type)       :: fich_cfl
 type (Space_Group_Type)     :: SpG
 type (MagSymm_k_Type)       :: MGp
 type (Atom_list_Type)       :: A
 type (MAtom_list_Type)      :: Am
 type (Crystal_Cell_Type)    :: Cell
 type (MagHD_Type)           :: Mh
 type (Magnetic_domain_type) :: Mag_Dom
 Type (Polar_Info_type)      :: polari

 character(len=256)          :: filcod     !Name of the input file
 character(len=1)            :: sig
 real                        :: sn,sf2
 real, dimension(3)          :: vpl
 integer                     :: Num, lun=1, ier,i,j,m,ih,ik,il,iv, n_ini,n_end, &
                                ich, nd, nch
 complex                     :: fc

 integer                     :: narg,iargc
 Logical                     :: esta, arggiven=.false.

      !---- Arguments on the command line ----!
      !---- Arguments on the command line ----!
      narg=COMMAND_ARGUMENT_COUNT()

      if(narg > 0) then
              call GET_COMMAND_ARGUMENT(1,filcod)
              arggiven=.true.
      end if

     write(unit=*,fmt="(/,/,8(a,/))")                                                  &
           "              ------ P r o g r a m    M a g P o l a r (DOMAINS) ------"  , &
           "                    ---- Version 0.1 November-2006 ----"                 , &
           "    ******************************************************************"  , &
           "    * Calculates magnetic structure factors, and magnetic interaction*"  , &
           "    * vectors from magnetic structures by reading a *.CFL file       *"  , &
           "    * Magnetic S-domains and chirality domains are considered        *"  , &
           "    ******************************************************************"  , &
           "                    (JRC- November 2006, testing stage )"
    write(unit=*,fmt=*) " "

     if(.not. arggiven) then
       write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.cfl (give xx): "
       read(unit=*,fmt="(a)") filcod
       if(len_trim(filcod) == 0) stop
     end if

     open(unit=lun,file=trim(filcod)//".cal", status="replace",action="write")
     write(unit=lun,fmt="(/,/,8(a,/))")                                                &
           "              ------ P r o g r a m   M a g P o l a r  (DOMAINS) ------"  , &
           "                    ---- Version 0.1 November-2006 ----"                 , &
           "    ******************************************************************"  , &
           "    * Calculates magnetic structure factors, and magnetic interaction*"  , &
           "    * vectors from magnetic structures by reading a *.CFL file       *"  , &
           "    * Magnetic S-domains and chirality domains are considered        *"  , &
           "    ******************************************************************"  , &
           "                    (JRC- November 2006, testing stage )"

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

       n_ini=1
       n_end=fich_cfl%nlines
       call Readn_Set_Magnetic_Structure(fich_cfl,n_ini,n_end,MGp,Am,Mag_dom=Mag_dom)
       if(err_MagSym) then
         write(unit=*,fmt="(a)") "   "//err_MagSym_Mess
         stop
       end if
       if(Mag_dom%nd == 0) then
          Mag_dom%nd=1
          Mag_Dom%chir=.false.
          do i=1,3
           Mag_Dom%DMat(i,i,1)=1.0
          end do
          Mag_Dom%pop(1,1)= 1.0
       end if
       call Write_Magnetic_Structure(lun,MGp,Am,Mag_dom=Mag_dom)
       nch=1
       if(Mag_Dom%chir) nch=2

      do

         write(unit=*,fmt="(/a,i2,a)",advance="no") &
                            " => Enter a magnetic reflections as 4 integers -> (h,k,l,m)=H+sign(m)*k(abs(m)): "
         read(unit=*,fmt=*) ih,ik,il,m
         if( m == 0) exit
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

         !Calculate magnetic structure factor and magnetic interaction vector
         call Calc_Magnetic_StrF_MiV_Dom(Cell,MGp,Am,Mag_Dom,Mh)

         write(unit=lun,fmt="(/,a,3i4,a,3f8.4,a)") "  Reflection: (",ih,ik,il,") "//sig//" (",MGp%kvec(:,iv),")"
         write(unit=lun,fmt="(a,3f8.4,a)")         "              (",Mh%h,")"
         write(unit=*,  fmt="(/,a,3i4,a,3f8.4,a)") "  Reflection: (",ih,ik,il,") "//sig//" (",MGp%kvec(:,iv),")"
         write(unit=*,  fmt="(a,3f8.4,a)")         "              (",Mh%h,")"

         do nd=1,Mag_Dom%nd
           do ich=1,nch
              write(unit=lun,fmt="(/2(a,i3))")      "  => Magnetic Domain #",nd,"  Chirality Domain #",ich
              write(unit=lun,fmt="(a,2(3f8.4,a))") "  Magnetic Structure   Factor : (",real(Mh%MsF(:,ich,nd)),")+i(",&
                                                      aimag(Mh%MsF(:,ich,nd)),") "
              write(unit=lun,fmt="(a,2(3f8.4,a))") "  Magnetic Interaction Vector : (",real(Mh%MiV(:,ich,nd)),")+i(",&
                                                      aimag(Mh%MiV(:,ich,nd)),") "
              write(unit=*,  fmt="(/2(a,i3))")      "  => Magnetic Domain #",nd,"  Chirality Domain #",ich
              write(unit=*,  fmt="(a,2(3f8.4,a))") "  Magnetic Structure   Factor : (",real(Mh%MsF(:,ich,nd)),")+i(",&
                                                      aimag(Mh%MsF(:,ich,nd)),") "
              write(unit=*,  fmt="(a,2(3f8.4,a))") "  Magnetic Interaction Vector : (",real(Mh%MiV(:,ich,nd)),")+i(",&
                                                      aimag(Mh%MiV(:,ich,nd)),") "
           end do
         end do

         write(unit=lun,fmt="(/a,2(3f8.4,a))") "  Average Magnetic Interaction Vector : (",real(Mh%AMiV(:)),")+i(",&
                                                  aimag(Mh%AMiV(:)),") "
         write(unit=lun,fmt="(a,f12.5 )")      "  Square of Average Mag. int.  Vector : ",Mh%sqAMiV
         write(unit=lun,fmt="(a,f12.5 )")      "  Average Square of Mag. int.  Vectors: ",Mh%sqMiV
         write(unit=*,  fmt="(/a,2(3f8.4,a))") "  Average Magnetic Interaction Vector : (",real(Mh%AMiV(:)),")+i(",&
                                                  aimag(Mh%AMiV(:)),") "
         write(unit=*,  fmt="(a,f12.5 )")      "  Square of Average Mag. int.  Vector : ",Mh%sqAMiV
         write(unit=*,  fmt="(a,f12.5 )")      "  Average Square of Mag. int.  Vectors: ",Mh%sqMiV
         do
          write(unit=*,  fmt="(/a)") "    Polarisation Matrix calculation: "
          write(unit=*,  fmt="( a)",advance="no") " => Enter another reciprocal lattice vector in the horizontal plane: "
          read(unit=*,fmt=*,iostat=ier) vpl
          if(ier /= 0) vpl=(/0.0,0.0,0.0/)
          if(dot_product(vpl,vpl) < 0.00001) then
            write(unit=*,  fmt="(a)") "    The given vector has module = 0, please provide a non-zero vector!"
            cycle
          end if
          exit
         end do
         call Set_Polar_Info(Cell, Mh%h, vpl, 1.0, cmplx(0.0,0.0), Mh%AMiV, Polari)
         call Write_Polar_Info(Polari, info="p")
         call Write_Polar_line(Polari)
       end do


       write(unit=*,fmt="(a)") " Normal End of program: MagPolar3D"
       write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".cal"
     end if

     close(unit=lun)
     stop
End Program MagPolar3D
