!!----
!!----
!!----
SubModule (CFML_IOForm) IO_CFL
   !---- Variables ----!
   implicit none

 Contains
   !!----
   !!---- READ_CFL_ATOM
   !!----    Subroutine to read atoms parameters
   !!----
   !!----         ATOM   Label  ChemSymb   x y z B Occ Us or Moment
   !!----
   !!----     For charge obtained from Label: Label[+/-][number]
   !!----
   !!---- 07/05/2020
   !!
   Module Subroutine Read_CFL_Atoms(cfl, AtmList, Type_Atm, d)
      !---- Arguments ----!
      type(File_Type),      intent(in)     :: cfl     ! Containing information
      Type(AtList_Type),    intent(out)    :: AtmList
      character(len=*),     intent(in)     :: Type_Atm
      integer,              intent(in)     :: d

      !---- Local variables -----!
      character(len=:),   allocatable   :: line,mom_comp
      character(len=:),   allocatable   :: dire
      integer                           :: i, j, na, npos, n_oc, n_mc,n_dc,n_uc

      !> Init
      call clear_error()
      call Allocate_Atom_List(0, Atmlist, Type_Atm, d)

      !> Calculate number of Atoms
      na=0
      mom_comp=" "
      do i=1,cfl%nlines
          line=adjustl(cfl%line(i)%str)
          if (len_trim(line) == 0) cycle
          if (line(1:1) == "!" .or. line(1:1) == "#") cycle

          if (u_case(line(1:12)) == "ATM_MOM_COMP") then
             j=index(line,"!")
             if ( j /= 0) then
                mom_comp=adjustl(line(13:j-1))
             else
                mom_comp=adjustl(line(13:))
             end if
          end if
          npos=index(u_case(line),'ATOM')
          if (npos > 0) na=na+1
      end do
      if (na == 0) return             ! No atoms in the lines

      !> Allocate List
      call Allocate_Atom_List(na, Atmlist, Type_Atm, d)
      if (len_trim(mom_comp) > 2) Atmlist%mcomp=mom_comp

      na=0
      do i=1,cfl%nlines
         line=adjustl(cfl%line(i)%str)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == "!" .or. line(1:1) == "#") cycle

         !> Truncate line from symbols: # and !
         npos=index(line,'!')
         if (npos > 0) line=line(:npos-1)

         !> Eliminate Tabs
         do
            npos=index(line,TAB)
            if (npos == 0) exit
            line(npos:npos)=' '
         end do

         !> ATOM Directive
         dire=adjustl(u_case(line(1:4)))
         if (trim(dire) /= "ATOM") cycle

         na=na+1
         call read_atom(line, Atmlist%atom(na))
         Atmlist%atom(na)%UType="B"
         Atmlist%atom(na)%ThType="ISO"

         !> Trial to read anisotropic thermal and
         !> magnetic moment parameters
         j=i
         do
            j=j+1
            if ( j < cfl%nlines ) then
               line=adjustl(cfl%line(j)%str)
               if (u_case(line(1:4)) == "ATOM") exit

               npos=index(line," ")
               if (npos <= 1) then
                  err_CFML%Ierr=1
                  err_CFML%Msg="Read_CFL_Atoms: Error in line "//trim(line)
                  return
               end if

               select case (u_case(line(1:npos-1)))
                  case ("MOMENT")
                     call read_Moment(line,Atmlist%atom(na))

                  case ("U_IJ")
                     call read_UTherms(line,Atmlist%atom(na))
                     Atmlist%atom(na)%UType= "U"
                     Atmlist%atom(na)%ThType= "ANI"

                  case ("B_IJ")
                     call read_UTherms(line,Atmlist%atom(na))
                     Atmlist%atom(na)%UType="B"
                     Atmlist%atom(na)%ThType="ANI"

                  case ("BETA")
                     call read_UTherms(line,Atmlist%atom(na))
                     Atmlist%atom(na)%UType= "BETA"
                     Atmlist%atom(na)%ThType="ANI"
               end select
               if (Err_CFML%Ierr /= 0) return

            else
               exit
            end if
         end do

         select type(at => Atmlist%atom)
            class is(MAtm_Std_Type)
               n_oc=0; n_mc=0; n_dc=0; n_uc=0
               j=i
               do
                  j=j+1
                  if ( j < cfl%nlines ) then
                     line=adjustl(cfl%line(j)%str)
                     if (u_case(line(1:4)) == "ATOM") exit

                     npos=index(line," ")
                     if (npos <= 1) then
                        err_CFML%Ierr=1
                        err_CFML%Msg="Read_CFL_Atoms: Error in line "//trim(line)
                        return
                     end if

                     select case (u_case(line(1:npos-1)))
                        case ("O_CS")
                           n_oc=n_oc+1
                           call read_modulation_amplitudes(line,At(na),"O_CS",n_oc)

                        case ("M_CS")
                           n_mc=n_mc+1
                           call read_modulation_amplitudes(line,At(na),"M_CS",n_mc)

                        case ("D_CS")
                           n_dc=n_dc+1
                           call read_modulation_amplitudes(line,At(na),"D_CS",n_dc)

                        case ("U_CS")
                           n_uc=n_uc+1
                           call read_modulation_amplitudes(line,At(na),"U_CS",n_uc)
                     end select
                     if (Err_CFML%Ierr /= 0) return
                  else
                     exit
                  end if
               end do
               At(na)%n_oc=n_oc
               At(na)%n_mc=n_mc
               At(na)%n_dc=n_dc
               At(na)%n_uc=n_uc
         end select
      end do
      AtmList%natoms=na

   End Subroutine Read_CFL_Atoms

   !!----
   !!---- READ_CFL_CELL
   !!----
   !!----    Obtaining Cell Parameter from CFL Format
   !!----
   !!---- 07/05/2020
   !!
   Module Subroutine Read_CFL_Cell(cfl, Cell, CFrame)
      !---- Arguments ----!
      type(File_Type),            intent(in)     :: cfl     ! Containing information
      class(Cell_Type),           intent(out)    :: Cell    ! Cell object
      character(len=*), optional, intent( in)    :: CFrame

      !---- Local variables -----!
      integer                              :: i, iv,n_ini,n_end
      real(kind=cp), dimension (6)         :: vcell, std
      character(len=132),dimension(1)      :: lines

      !> Init
      call clear_error()
      if (cfl%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: 0 lines "
         return
      end if

      !> Search: CELL
      do i=1,cfl%nlines
         lines(1)=adjustl(u_case(cfl%line(i)%str))
         if (lines(1)(1:4) == "CELL") exit
         lines(1)=" "
      end do

      if (len_trim(lines(1)) == 0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: Instruction 'CELL' not provided "
         return
      end if

      !> Eliminate Tabs
      do
         iv=index(lines(1),TAB)
         if (iv == 0) exit
         lines(1)(iv:iv)=' '
      end do

      n_ini=1; n_end=1
      vcell=0.0
      std=0.0
      call Read_Key_ValueSTD(lines, n_ini, n_end,"CELL", vcell, std, iv)
      if (iv /= 6) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: Problems reading cell parameters!"
         return
      end if
      if (present(CFrame)) then
         call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, CarType=CFrame, Vscell=std(1:3), Vsang=std(4:6))
      else
         call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, Vscell=std(1:3), Vsang=std(4:6))
      end if

   End Subroutine Read_CFL_Cell

   !!----
   !!---- READ_CFL_KVECTORS
   !!----
   !!---- Read K-vectors information
   !!----
   !!---- 07/05/2020
   !!
   Module Subroutine Read_CFL_KVectors(cfl, Kvec)
      !---- Arguments ----!
      type(File_Type),         intent(in)     :: cfl
      type(kvect_info_Type),   intent(out)    :: Kvec

      !---- Local Variables ----!
      integer                      :: i,j,ier,nk,nq,iv
      character(len=:),allocatable :: uline,line

      !> Init
      call clear_error()
      if (cfl%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: 0 lines "
         return
      end if

      nk=0; nq=0
      do i=1,cfl%nlines
         line=adjustl(cfl%line(i)%str)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == "!" .or. line(1:1) == "#") cycle

         !> Eliminate Tabs
         do
            iv=index(line,TAB)
            if (iv == 0) exit
            line(iv:iv)=' '
         end do

         j=index(line,"!")
         if(j /= 0) line=line(:j-1)

         j=index(line," ")
         if ( j == 0) then
            uline=u_case(line)
         else
            uline=u_case(line(:j-1))
         end if
         line=adjustl(line(j+1:))

         select case(trim(uline))
            case("NQVECT","NKVECT","NKVEC","NQVEC")
               read(unit=line,fmt=*,iostat=ier) Kvec%nk, Kvec%nq
               if (ier /= 0) then
                  Err_CFML%Ierr=1
                  Err_CFML%Msg="Error reading the number of k-vectors and/or number of Q-coefficients"
                  return
               end if
               allocate(Kvec%kv(3,Kvec%nk),Kvec%q_coeff(Kvec%nk,Kvec%nq))
               allocate(Kvec%nharm(Kvec%nk),Kvec%sintlim(Kvec%nk))
               Kvec%kv=0.0_cp; Kvec%q_coeff=1; Kvec%nharm=1; Kvec%sintlim=1.0

            case("QVECT","KVECT","KVEC","QVEC")
               if (Kvec%nk > 0) then
                  nk=nk+1
                  read(unit=line,fmt=*,iostat=ier) Kvec%kv(:,nk)
                  if (ier /= 0) then
                     Err_CFML%Ierr=1
                     write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the k-vector #",nk
                     return
                  end if
               end if

           case("NHARM")
              if (Kvec%nk > 0) then
                 read(unit=line,fmt=*,iostat=ier) Kvec%nharm(1:Kvec%nk)
                 if (ier /= 0) then
                    Err_CFML%Ierr=1
                    Err_CFML%Msg = "Error reading the nk harmonics !"
                    return
                 end if
              end if

           case("SINTL")
              if (Kvec%nk > 0) then
                 read(unit=line,fmt=*,iostat=ier) Kvec%sintlim(1:Kvec%nk)
                 if (ier /= 0) then
                    Err_CFML%Ierr=1
                    Err_CFML%Msg = "Error reading the maximum sinTheta/Lambda for harmonics!"
                    return
                 end if
              end if

           case("Q_COEFF")
              nq=nq+1
              read(unit=line,fmt=*,iostat=ier) Kvec%q_coeff(1:Kvec%nk,nq)
              if (ier /= 0) then
                 Err_CFML%Ierr=1
                 write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the Q-coefficent # ",nq
                 return
              end if
         end select
      end do

      if (Kvec%nk /= nk) then
         Err_CFML%Ierr=1
         write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of k-vectors,",Kvec%nk, ", does not correspond with the prescribed number: ",nk
         return
      end if

      if (Kvec%nq /= nq) then
         Err_CFML%Ierr=1
         write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of expected Q-coefficients,",Kvec%nq, ", does not correspond with number of read Q-coefficients ",nq
         return
      end if

   End Subroutine Read_CFL_KVectors

   !!----
   !!---- READ_CFL_SPG
   !!----
   !!---- Read Space group information ina CFL file
   !!----
   !!---- 07/05/2020
   !!
   Module Subroutine Read_CFL_SpG(cfl, SpG, xyz_type)
       !---- Arguments ----!
       Type(File_Type),                 intent(in)     :: cfl
       class(SpG_Type),                 intent(out)    :: SpG
       character(len=*), optional,      intent(in)     :: xyz_type

       !--- Local Variables ---!
       integer                           :: i,j,ngen,nk,nq,iv,ier
       character(len=:),     allocatable :: line,uline,setting,strcode
       character(len=40), dimension(192) :: gen
       logical                           :: change_setting

       !> Init
       call clear_error()

       !> Look for the appropriate keywords to construct the space group:
       !> Crystallographic, Shubnikov, or superspace
       ngen=0
       setting=" "
       change_setting=.false.

       strcode="xyz"
       if (present(xyz_type)) strcode=trim(xyz_type)

       do i=1,cfl%nlines
          line=adjustl(cfl%line(i)%str)
          if (len_trim(line) == 0) cycle
          if (line(1:1) == "!" .or. line(1:1) == "#") cycle

          !> Eliminate Tabs
          do
             iv=index(line,TAB)
             if (iv == 0) exit
             line(iv:iv)=' '
          end do

          j=index(line,"!")
          if (j /= 0) line=line(:j-1)

          j=index(line,"::")
          if (j /= 0) then
             setting=trim(adjustl(line(j+2:)))
             if (len_trim(setting) /= 0) change_setting=.true.
             line=line(:j-1)
          end if

          j=index(line," ")
          uline=u_case(line(:j-1))

          line=adjustl(line(j+1:))
          select case(trim(uline))
             case("HALL","SPGR","SPACEG")
                call Set_SpaceGroup(line, SpG)
                exit

             case("SHUB")
                call Set_SpaceGroup(line,"SHUBN",SpG)
                exit

             case("SSG","SUPER","SSPG")
                call Set_SpaceGroup(line,"SUPER",SpG, strcode)
                exit

             case("GENLIST","GENERATORS","LIST")
                call Set_SpaceGroup(line,SpG)
                exit

             case("GEN","SYMM")
                ngen=ngen+1
                gen(ngen)=line
          end select
       end do
       if (ngen > 0) call Set_SpaceGroup("  ",SpG,ngen,gen)
       if (Err_CFML%Ierr == 1) return

       if (change_setting) then
          if (strcode == "xyz")  then
             call Change_Setting_SpaceG(setting, SpG)
          else
             call Change_Setting_SpaceG(setting, SpG,strcode)
          end if
       end if
       if (Err_CFML%Ierr == 1) return

       !> Now read q-vectors and other items if the class of SpG is SuperSpaceGroup_Type
       Select Type (SpG)
          Class is (SuperSpaceGroup_Type)
             if (allocated(SpG%kv))      deallocate (SpG%kv)
             if (allocated(SpG%nharm))   deallocate (SpG%nharm)
             if (allocated(SpG%sintlim)) deallocate (SpG%sintlim)
             if (allocated(SpG%Om))      deallocate (SpG%Om)
             if (allocated(SpG%q_coeff)) deallocate (SpG%q_coeff)

             allocate(SpG%Om(SpG%D,Spg%D,SpG%Multip))
             do i=1,SpG%Multip
                SpG%Om(:,:,i)=real(SpG%Op(i)%Mat(:,:))
             end do

             nk=0; nq=0
             do i=1,cfl%nlines
                line=adjustl(cfl%line(i)%str)
                if (len_trim(line) == 0) cycle
                if (line(1:1) == "!" .or. line(1:1) == "#") cycle

                j=index(line,"!")
                if (j /= 0) line=line(:j-1)

                j=index(line," ")
                uline=u_case(line(:j-1))

                line=adjustl(line(j+1:))
                select case(trim(uline))
                   case("NQVECT","NKVECT")
                      read(unit=line,fmt=*,iostat=ier) Spg%nk, Spg%nq
                      if (ier /= 0) then
                         Err_CFML%Ierr=1
                         Err_CFML%Msg="Error reading the number of k-vectors and/or number of Q-coefficients"
                         return
                      end if
                      allocate(Spg%kv(3,Spg%nk),Spg%q_coeff(Spg%nk,Spg%nq))
                      allocate(Spg%nharm(Spg%nk),Spg%sintlim(Spg%nk))
                      SpG%kv=0.0_cp; SpG%q_coeff=1; Spg%nharm=1; Spg%sintlim=1.0

                   case("QVECT","KVECT")
                      if (Spg%nk > 0) then
                         nk=nk+1
                         read(unit=line,fmt=*,iostat=ier) Spg%kv(:,nk)
                         if (ier /= 0) then
                            Err_CFML%Ierr=1
                            write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the k-vector #",nk
                            return
                         end if
                      end if

                   case("NHARM")
                      if (Spg%nk > 0) then
                         read(unit=line,fmt=*,iostat=ier) Spg%nharm(1:Spg%nk)
                         if (ier /= 0) then
                            Err_CFML%Ierr=1
                            Err_CFML%Msg = "Error reading the nk harmonics !"
                            return
                         end if
                      end if

                   case("SINTL")
                      if (Spg%nk > 0) then
                         read(unit=line,fmt=*,iostat=ier) Spg%sintlim(1:Spg%nk)
                         if (ier /= 0) then
                            Err_CFML%Ierr=1
                            Err_CFML%Msg = "Error reading the maximum sinTheta/Lambda for harmonics!"
                            return
                         end if
                      end if

                   case("Q_COEFF")
                      nq=nq+1
                      read(unit=line,fmt=*,iostat=ier) Spg%q_coeff(1:Spg%nk,nq)
                      if (ier /= 0) then
                         Err_CFML%Ierr=1
                         write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the Q-coefficent # ",nq
                         return
                      end if

                end select
             end do

             if (Spg%nk /= (Spg%D-4)) then
                Err_CFML%Ierr=1
                write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of k-vectors,",Spg%nk, ", does not correspond with the additional dimensions of the group ",Spg%D-4
                return
             end if

             if (Spg%nq /= nq) then
                Err_CFML%Ierr=1
                write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of expected Q-coefficients,",Spg%nq, ", does not correspond with number of read Q-coefficients ",nq
                return
             end if

       End Select

   End Subroutine Read_CFL_SpG

   !!----
   !!---- WRITE_CFL_ATOMS
   !!----
   !!----    Write the atoms in the asymmetric unit for a CFL file
   !!----
   !!---- 08/05/2020
   !!
   Module Subroutine Write_CFL_Atoms(AtmList, Lun, Cell)
      !---- Arguments ----!
      Type(AtList_Type),            intent(in) :: AtmList
      integer,            optional, intent(in) :: Lun
      class(Cell_G_Type), optional, intent(in) :: Cell

      !---- Local Variables ----!
      character(len=36)              :: forma,fom
      character(len=30),dimension(6) :: text
      real(kind=cp), dimension(6)    :: u,bet,sb
      integer                        :: i, j, iunit, leng, maxl,ish

      !> Unit
      iunit=6
      if (present(lun)) iunit=lun

      if (AtmList%natoms == 0) then
         write (unit=iunit,fmt="(a)") " There aren't atoms defined!"
         return
      end if

      !> Determine the maximum length of the atom labels
      maxl=0
      do i=1,AtmList%natoms
         leng=len_trim(atmList%atom(i)%lab)
         if (leng > maxl) maxl=leng
      end do
      maxl=max(maxl,4)+1

      !> define format fom
      fom   ="(a,tr  ,a)"
      forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"

      ish=maxl-4
      select case(ish)
         case(:9)
            write(unit=fom(6:6),fmt="(i1)") ish

         case(10:)
            write(unit=fom(6:7),fmt="(i2)") ish
      end select

      select case(maxl)
         case(:9)
            write(unit=forma(5:5),fmt="(i1)") maxl
         case(10:)
            write(unit=forma(5:6),fmt="(i2)") maxl
      end select

      write (unit=iunit,fmt=fom) "!     ", &
            "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"

      select type (at => AtmList%atom)
         class is (Atm_Std_Type)
            do i=1,AtmList%natoms

               do j=1,3
                  text(j)=string_NumStd(at(i)%x(j), at(i)%x_std(j))
               end do
               text(4)=string_NumStd(at(i)%U_iso, at(i)%U_iso_std)
               text(5)=string_NumStd(at(i)%Occ, at(i)%Occ_std)

               write (unit=iunit,fmt=forma) "Atom   ",trim(at(i)%lab),at(i)%ChemSymb, &
                     (text(j),j=1,5), at(i)%moment, at(i)%charge, " # "//trim(at(i)%AtmInfo)

               select case (l_case(at(i)%ThType))
                  case ('ani')
                     select case (l_case(at(i)%UType))
                        case ('beta')
                           u=at(i)%u(1:6)
                           if (present(cell)) then
                              bet=Get_U_from_Betas(u,cell)
                              sb=Get_U_from_Betas(at(i)%u_std,cell)
                           else
                              bet=u
                               sb=at(i)%u_std(1:6)
                           end if
                           do j=1,6
                              text(j)=string_NumStd(bet(j), sb(j))
                           end do
                           if (present(cell)) then
                              write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
                           else
                              write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                           end if

                        case ('u_ij')
                           u=at(i)%u(1:6)
                           if (present(cell)) then
                              bet=Get_Betas_from_U(u,cell)
                              sb=Get_Betas_from_U(at(i)%u_std,cell)
                           else
                              bet=u
                               sb=at(i)%u_std(1:6)
                           end if
                           do j=1,6
                              text(j)=string_NumStd(bet(j), sb(j))
                           end do
                           if (present(cell)) then
                              write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
                           else
                              write (unit=iunit,fmt="(a,tr1,6a14)") "U_ij  ", text
                           end if

                     end select
               end select
            end do
      end select

   End Subroutine Write_CFL_Atoms

   !!----
   !!---- Write_CFL_File
   !!----
   !!----    Write a CFL file
   !!----
   !!---- 08/05/2020
   !!
   Module Subroutine Write_CFL_File(Lun,Cell, SpG, Atm, Title)
      !---- Arguments ----!
      integer,                     intent(in)    :: lun
      class(Cell_G_Type),          intent(in)    :: Cell
      class(SpG_Type),             intent(in)    :: SpG
      Type(AtList_Type), optional, intent(in)    :: Atm
      character(len=*),  optional, intent(in)    :: Title

      !----- Local variables -----!
      integer                         :: j
      real(kind=cp), dimension(6)     :: a,sa
      character(len=30), dimension(6) :: text

      !> Title
      if (present(title)) write(unit=lun,fmt="(a)") "TITLE "//trim(title)

      write(unit=lun,fmt='(a)')" "
      write(unit=lun,fmt='(a)') "!  Automatically generated CFL file (Write_CFL)"
      write(unit=lun,fmt='(a)')" "

      !> Cell
      a(1:3)=Cell%Cell
      a(4:6)=Cell%ang
      sa(1:3)=Cell%scell
      sa(4:6)=Cell%sang
      do j=1,6
         text(j)=string_NumStd(a(j), sa(j))
      end do
      write(unit=lun,fmt="(a)") "!         a               b               c            alpha           beta            gamma"
      write(unit=lun,fmt="(a,6a16)") "Cell ",text
      write(unit=lun,fmt='(a)')" "

      !> Space group
      write(unit=lun,fmt="(a,i3)")"!     Space Group # ",SpG%NumSpg
      write(unit=lun,fmt="(a,a)") "Spgr  ",SpG%spg_symb
      write(unit=lun,fmt='(a)')" "

      !> Atoms
      if (present(Atm)) then
         call Write_CFL_Atoms(Atm,Lun,cell)
         write(unit=lun,fmt='(a)')" "
      end if

   End Subroutine Write_CFL_File

End SubModule IO_CFL