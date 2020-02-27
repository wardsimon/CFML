!!----
!!----
!!----
SubModule (CFML_IOForm) IOF_CFL
   implicit none
   Contains

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Class (SPG_Type),             intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++    Type(Job_Info_type), optional,intent(out)  :: Job_Info
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CFL File
    !!--++
    !!--++ Update: February - 2020
    !!
    Module Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,Type_Atm,CFrame,NPhase,Job_Info,xyz_type)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       class (Cell_Type),            intent(out)  :: Cell
       class (SPG_Type),             intent(out)  :: SpG
       Type (AtList_Type),           intent(out)  :: A
       character(len=*),             intent(in)   :: Type_Atm
       character(len=*),    optional,intent(in)   :: CFrame
       Integer,             optional,intent(in)   :: Nphase
       Type(Job_Info_type), optional,intent(out)  :: Job_Info
       character(len=*),    optional,intent(in)   :: xyz_type

       !---- Local variables ----!
       character(len=132)               :: line
       integer                          :: i, nauas, ndata, iph, n_ini,n_end,k
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip
       real(kind=cp),dimension(:),allocatable:: vet

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!
       iph=1
       if (present(nphase)) iph=nphase
       if (present(Job_Info)) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Get_Job_Info(file_dat,n_ini,n_end,Job_info)
       end if

       !---- Reading Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if(present(CFrame)) then
         call Read_CFL_Cell(file_dat,n_ini,n_end,Cell,CFrame) !Read and construct Cell
       else
         call Read_CFL_Cell(file_dat,n_ini,n_end,Cell) !Read and construct Cell
       end if
       if (err_CFML%Ierr /= 0) return

       !---- Reading Space Group Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       call Read_CFL_SpG(file_dat,n_ini,n_end,SpG)
       if (err_CFML%Ierr /= 0) return

       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       !---- Calculating number of Atoms in the Phase ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (l_case(line(1:4)) == "atom")  nauas=nauas+1
       end do

       if (nauas > 0) then
          call Read_CFL_Atoms(file_dat,n_ini,n_end,A,Type_Atm,SpG%D-1)
          if (err_CFML%Ierr /= 0) return
          if(allocated(vet)) deallocate(vet)
          Select Type (SpG)
            type is (SuperSpaceGroup_Type)
              allocate(vet(SpG%D-1))
              do i=1,A%natoms
                 vet(1:3)=A%atom(i)%x
                 do k=1,Spg%nk
                   vet(3+k)=dot_product(vet(1:3),SpG%kv(:,k))
                 end do
                 A%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
                 if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)
                 Select Type (at => A%atom(i))
                   Class is (MAtm_Std_Type)
                     At%Xs=vet
                 End Select
              end do
            class Default
              allocate(vet(3))
              do i=1,A%natoms
                 vet=A%atom(i)%x
                 A%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
                 if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)
              end do
          End Select

          !Select Type (Cell)   !This piece of code should be executed before structure factor calculations
          !  Type is(Cell_G_Type)
          !     do i=1,A%natoms
          !        if (A%atom(i)%thtype == "ani") then  !For Structure factor calculations we use betas
          !           select case (trim(A%atom(i)%Utype))
          !              case ("U")
          !                 A%atom(i)%u(1:6) =  Get_Betas_from_U(A%atom(i)%u(1:6),Cell)
          !              case ("B")
          !                 A%atom(i)%u(1:6) =  Get_Betas_from_B(A%atom(i)%u(1:6),Cell)
          !           end select
          !           A%atom(i)%Utype="beta"
          !        end if
          !     end do
          !End Select
       end if
    End Subroutine Readn_Set_XTal_CFL
   !!----
   !!---- READ_CFL_ATOM
   !!----    Subroutine to read atoms parameters
   !!----
   !!----    Format: ATOM   Label  ChemSymb   x y z B Occ Us or Moment
   !!----
   !!----     For charge obtained from Label: Label[+/-][number]
   !!----
   !!----
   !!---- 25/06/2019
   !!
   Module Subroutine Read_CFL_Atoms(lines,n_ini, n_end, At_List,Type_Atm,d)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in)     :: lines
      integer,                        intent(in out) :: n_ini
      integer,                        intent(in)     :: n_end
      Type (AtList_Type),             intent(out)    :: At_List
      character(len=*),               intent(in)     :: Type_Atm
      integer,                        intent(in)     :: d

      !---- Local variables -----!
      character(len=12),  parameter     :: DIGPM="0123456789+-"
      character(len=:),   allocatable   :: line,mom_comp
      character(len=180), dimension(1)  :: extline
      character(len=:),   allocatable   :: dire
      character(len=:),   allocatable   :: label
      character(len=2)                  :: Chemsymb
      integer                           :: i, j, n, na, npos, nlong, iv,n_oc, n_mc,n_dc,n_uc
      real(kind=cp), dimension (20)     :: vet1, vet2


      !> Init
      call clear_error()
      call Allocate_Atom_List(0, At_list, Type_Atm,d)

      !> Calculate number of Atoms
      na=0
      mom_comp=" "
      do i=n_ini, n_end
         !> Checks
         if (len_trim(lines(i)) <=0) cycle                 ! Empty line
         line=adjustl(lines(i))
         if (line(1:1) =='!' .or. line(1:1) =='#') cycle   ! Comment line
         if(u_case(line(1:12)) == "ATM_MOM_COMP") then
           j=index(line,"!")
           if( j /= 0) then
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
      call Allocate_Atom_List(na, At_list,Type_Atm,d)
      if(len_trim(mom_comp) > 2) At_list%mcomp=mom_comp
      write(*,"(a)") " Moment components: "//At_list%mcomp
      !> Read Information
      na=0

      do i=n_ini, n_end
         !> Checks
         if (len_trim(lines(i)) <=0) cycle                 ! Empty line
         line=adjustl(lines(i))
         if (line(1:1) =='!' .or. line(1:1) =='#') cycle   ! Comment line

         !> Truncate line from symbols: # and !
         npos=index(line,'!')                             ! Cut line according to !
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
          call read_atom(line,At_list%atom(na))
          At_list%atom(na)%UType="B"
          At_list%atom(na)%ThType="ISO"
          !---- Trial to read anisotropic thermal and magnetic moment parameters ----!
          j=i
          do
            j=j+1
            if( j < size(lines) ) then
               line=adjustl(lines(j))
               if(u_case(line(1:4)) == "ATOM") exit
               npos=index(line," ")
               select case (u_case(line(1:npos-1)))
                 case ("MOMENT")
                    call Read_Moment(line,At_list%atom(na))
                 case ("U_IJ")
                    call read_uvals(line,At_list%atom(na),"U_IJ")
                    At_list%atom(na)%UType= "U"
                    At_list%atom(na)%ThType= "ANI"
                 case ("B_IJ")
                    call read_uvals(line,At_list%atom(na),"B_IJ")
                    At_list%atom(na)%UType="B"
                    At_list%atom(na)%ThType="ANI"
                 case ("BETA")
                    call read_uvals(line,At_list%atom(na),"BETA")
                    At_list%atom(na)%UType= "BETA"
                    At_list%atom(na)%ThType="ANI"
               end select
               if(Err_CFML%Ierr /= 0) exit
            else
              exit
            end if
          end do

          Select Type(at => At_list%atom)
           class is(MAtm_Std_Type)
             n_oc=0; n_mc=0; n_dc=0; n_uc=0
             j=i
             do
               j=j+1
               if( j < size(lines) ) then
                  line=adjustl(lines(j))
                  if(u_case(line(1:4)) == "ATOM") exit
                  npos=index(line," ")
                  !write(*,"(i3,a)") na, " -> "//At(na)%lab//u_case(line(1:npos-1))//" "//lines(j)(npos:)
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
                  if(Err_CFML%Ierr /= 0) exit
               else
                 exit
               end if
             end do
             At(na)%n_oc=n_oc
             At(na)%n_mc=n_mc
             At(na)%n_dc=n_dc
             At(na)%n_uc=n_uc
          End Select
      end do
      At_List%natoms=na

   End Subroutine Read_CFL_Atoms

    !!----
    !!---- Subroutine Read_Atom(Line,Atomo)
    !!----    character(len=*), intent(in out ) :: line    !  In -> Input String with ATOM directive
    !!----    Type (Atom_Type), intent(out)     :: Atomo   ! Out -> Parameters on variable Atomo
    !!----
    !!----    Subroutine to read the atom parameters from a given "line"
    !!----    it construct the object Atomo of type Atom.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Atom(line,Atom)
       !---- Arguments ----!
       character(len=*), intent(in out ) :: line
       class(Atm_Type),  intent(out)     :: Atom

       !---- Local variables -----!
       integer                           :: iv, nlong1,n,ier,q
       real(kind=cp), dimension (10)     :: vet1
       real(kind=cp), dimension (10)     :: vet2
       character(len=4)                  :: dire
       character(len=40)                 :: magmom
       character(len=:),   allocatable   :: label
       character(len=132), dimension(1)  :: filevar
       character(len=*), parameter       :: digpm="0123456789+-"

       !---- Init ----!
       call clear_error()
       q=0
       label="               "
       iv=index(line,"#")
       if(iv /= 0) then
        atom%AtmInfo=line(iv+1:)
        line=line(1:iv-1)
       end if
       iv=index(line,"Moment:") !Attemp to read magnetic moment components in the same line
       magmom=" "
       if(iv /= 0) then
         magmom=line(iv+7:)  !magmon should contain magnetic moment
         line=line(1:iv-1)   !Line after removing "Moment:" and information
       end if
       call Cut_String(line,nlong1,dire)
       if (u_case(dire) /= "ATOM") then
          err_CFML%Ierr=1
          err_CFML%Msg=" Error reading the ATOM keyword"
          return
       end if

       !---- Atom Label ----!
       call Cut_String(line,nlong1,label)
       atom%lab=label

       !---- Atom Type (Chemical symbol & Scattering Factor) ----!
       call Cut_String(line,nlong1,label)

       if((U_case(label(1:1)) == "M" .or. U_case(label(1:1)) == "J" ) & !Magnetic atom
          .and. index(digpm(1:10),label(4:4)) /= 0) then
          atom%ChemSymb=U_case(label(2:2))//L_case(label(3:3))
          atom%Magnetic=.true.
       else
          n=index(digpm,label(2:2))
          if (n /=0) then
            atom%ChemSymb=U_case(label(1:1))
          else
            atom%ChemSymb=U_case(label(1:1))//L_case(label(2:2))
          end if
       end if
       atom%SfacSymb=label(1:4)

       !---- Parameters ----!
       filevar(1)="atm "//trim(line)
       n=1
       call Read_Key_ValueSTD(filevar,n,n,"atm",vet1,vet2,iv)
       if (iv <= 0) then
          err_CFML%Ierr=1
          err_CFML%Msg="Error reading parameters of atom:"//atom%lab
          return
       end if

       !---- Coordinates  ----!
       if (iv < 3) then
          err_CFML%Ierr=1
          err_CFML%Msg="Error reading Coordinates of atom:"//atom%lab
          return
       end if

       atom%x(:)=vet1(1:3)             !---- Coordinates ----!
       if (iv > 3) atom%U_iso=vet1(4)  !---- Biso ----!
       if (iv > 4) atom%occ=vet1(5)    !---- Occ ----!
       if (iv > 5) atom%mom=vet1(6)    !---- Module of moment ----!
       if (iv > 6) atom%charge=vet1(7) !---- Charge ----!
       Select Type (Atom)
         class is (Atm_Std_Type)
           atom%x_std(:)=vet2(1:3)
           if (iv > 3) atom%U_iso_std=vet2(4)
           if (iv > 4) atom%occ_std=vet2(5)
       End Select

       !Attempt to get the oxidation state from "Atom%Label"
       label=atom%lab
       if(abs(atom%charge) < epsv) then
         iv=index(label,"+")
         Select Case(iv)
           Case(0) !No + sign
             n=index(label,"-")
             Select Case(n)
               Case(2) !Element with a single character symbol F-1
                  read(unit=label(3:),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
               Case(3) !Element in the form: F1- or Br-1
                  read(unit=label(2:2),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) then
                        read(unit=label(4:4),fmt="(i1)",iostat=ier)  q
                        if (ier /= 0) q=0
                  end if
               Case(4) !Element in the form: Br1-
                  read(unit=label(3:3),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
             End Select
             q=-q   !anions
           Case(2) !Element with a single character symbol C+4
                  read(unit=label(3:),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
           Case(3) !Element in the form: C4+ or Fe+3
                  read(unit=label(2:2),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) then
                        read(unit=label(4:4),fmt="(i1)",iostat=ier)  q
                        if (ier /= 0) q=0
                  end if
           Case(4) !Element in the form: Fe3+
                  read(unit=label(3:3),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
         End Select
         atom%charge=real(q)
       end if

       !Now read the components of the magnetic moment if given in the same line
       if(len_trim(magmom) /= 0) then
         filevar(1)="mom "//trim(magmom)
         call Read_Key_ValueSTD(filevar,n,n,"mom",vet1,vet2,iv)

         !---- Moment components  ----!
         if (iv < 3) then
            err_CFML%Ierr=1
            err_CFML%Msg="Error reading magnetic moment components of atom:"//atom%lab
            return
         end if
         atom%moment(:)=vet1(1:3)
         if(atom%mom < 0.001) atom%mom=maxval(abs(atom%moment(:)))
         Select Type (Atom)
           class is (Atm_Std_Type)
             atom%moment_std(:)=vet2(1:3)
         End Select
       end if

       return
    End Subroutine Read_Atom

    !!----
    !!---- Subroutine Read_Uvals(Line,Atomo,Ulabel)
    !!----    character(len=*),  intent(in out)  :: line      !  In -> String
    !!----    Type (Atom_Type),  intent(in out)  :: Atomo     !  In -> Atomo variable
    !!----                                                      Out ->
    !!----    character(len=4),  intent(in)      :: ulabel    !  In -> u_ij, b_ij, beta
    !!----
    !!----    Subroutine to read the anisotropic thermal parameters from a given Line
    !!----    it complets the object Atomo of type Atom.
    !!----    Assumes the string Line has been read from a file and
    !!----    starts with one of the words (u_ij, b_ij or beta), that is removed before reading
    !!----    the values of the parameters.
    !!----
    !!---- Update: February - 2020
    !!
    Subroutine Read_Uvals(Line,Atom,Ulabel)
       !---- Arguments ----!
       character(len=*),  intent(in )     :: line
       class(Atm_Type),   intent(in out)  :: Atom
       character(len=*),  intent(in)      :: ulabel
       !---- Local variables -----!
       character(len=len(line)),dimension(1):: line2
       real(kind=cp), dimension (6)         :: vet1,vet2
       integer                              :: iv,n

       call clear_error()

       line2(1)=line
       n=1
       call Cut_String(line2(1))
       line2(1)=trim(ulabel)//" "//line2(1)(1:len(line2(1))-5)  !this form of writing is to avoid gfortran warning -Wstring-overflow
       call Read_Key_ValueSTD(line2,n,n,trim(ulabel),vet1,vet2,iv)

        if (iv /= 6) then
          ERR_CFML%Ierr=1
          ERR_CFML%Msg="  Error reading the anisotropic thermal parameters of atom:"//atom%lab
          return
       end if
       atom%U(1:6)=vet1(1:6)
       Select Type(Atom)
         class is(Atm_std_Type)
          atom%U_std(1:6)=vet2(1:6)
       End Select

    End Subroutine Read_Uvals

    Subroutine Read_Moment(Line,Atom)
       !---- Arguments ----!
       character(len=*),  intent(in )     :: line
       class(Atm_Type),   intent(in out)  :: Atom

       !---- Local variables -----!
       character(len=len(line)),dimension(1):: line2
       real(kind=cp), dimension (6)         :: vet1,vet2
       integer                              :: iv,n

       call clear_error()

       line2(1)=line
       n=1
       call Cut_String(line2(1))
       !line2(1)="Moment "//line2(1)(1:len(line2(1))-5)  !this form of writing is to avoid gfortran warning -Wstring-overflow
       line2(1)="Moment "//line2(1)  !this form of writing is to avoid gfortran warning -Wstring-overflow
       write(*,"(a)")
       call Read_Key_ValueSTD(line2,n,n,"Moment",vet1,vet2,iv)

        if (iv /= 3) then
          ERR_CFML%Ierr=1
          ERR_CFML%Msg="  Error reading the Magnetic Moment of atom:"//atom%lab
          return
       end if
       atom%moment(1:3)=vet1(1:3)

       Select Type(at => Atom)
         type is(Atm_Std_Type)
           at%moment_std(1:3)=vet2(1:3)
         class is(MAtm_Std_Type)
           at%moment_std(1:3)=vet2(1:3)
       End Select

    End Subroutine Read_Moment

    Subroutine Read_Modulation_Amplitudes(Line,Atom,Ulabel,nt)
       !---- Arguments ----!
       character(len=*),    intent(in )     :: line
       class(MAtm_Std_Type),intent(in out)  :: Atom
       character(len=*),    intent(in)      :: ulabel
       integer,             intent(in)      :: nt  !number of the amplitude
       !---- Local variables -----!
       character(len=len(line)),dimension(1):: line2
       real(kind=cp), dimension (12)        :: vet1,vet2
       integer                              :: iv,n,nq,ier

       call clear_error()

       line2(1)=line
       n=1
       call Cut_String(line2(1))
       read(unit=line2(1),fmt=*,iostat=ier) nq
       if(ier /= 0) then
         Err_CFML%Ierr=1
         Err_CFML%Msg="Error reading the Q_coeff number in modulation amplitude line: "//trim(Ulabel)
         return
       end if
       call Cut_String(line2(1)) !Strip the first word, here "nq"
       !Reconstruct the line without "nq"
       line2(1)=trim(Ulabel)//" "//line2(1)(1:len(line2(1))-5)  !this form of writing is to avoid gfortran warning -Wstring-overflow
       call Read_Key_ValueSTD(line2,n,n,trim(Ulabel),vet1,vet2,iv)

       Select Case(trim(Ulabel))
         Case("O_CS")
            if (iv /= 2) then
              ERR_CFML%Ierr=1
              ERR_CFML%Msg="  Error reading the Modulation Amplitudes: "//trim(Ulabel)//" of atom:"//atom%lab
              return
            end if
            Atom%poc_q(nt) = nq
            Atom%Ocs(1:2,nt) = vet1(1:2)
         Case("M_CS")
            if (iv /= 6) then
              ERR_CFML%Ierr=1
              ERR_CFML%Msg="  Error reading the Modulation Amplitudes: "//trim(Ulabel)//" of atom:"//atom%lab
              return
            end if
            Atom%pmc_q(nt) = nq
            Atom%Mcs(1:6,nt) = vet1(1:6)
         Case("D_CS")
            if (iv /= 6) then
              ERR_CFML%Ierr=1
              ERR_CFML%Msg="  Error reading the Modulation Amplitudes: "//trim(Ulabel)//" of atom:"//atom%lab
              return
            end if
            Atom%pdc_q(nt) = nq
            Atom%Dcs(1:6,nt) = vet1(1:6)
         Case("U_CS")
            if (iv /= 12) then
              ERR_CFML%Ierr=1
              ERR_CFML%Msg="  Error reading the Modulation Amplitudes: "//trim(Ulabel)//" of atom:"//atom%lab
              return
            end if
            Atom%puc_q(nt) = nq
            Atom%Ucs(1:12,nt) = vet1(1:12)
       End Select

    End Subroutine Read_Modulation_Amplitudes

   !!----
   !!---- READ_CFL_CELL
   !!----    Obtaining Cell Parameter from CFL Format
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_CFL_Cell_Lines(lines, n_ini, n_end, Cell, CFrame)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      class(Cell_Type),                intent(out)    :: Cell    ! Cell object
      Character(len=*), optional,      intent( in)    :: CFrame
      !---- Local variables -----!
      integer                              :: iv
      real(kind=cp), dimension (6)         :: vcell, std

      !> Init
      call clear_error()

      call Read_Key_ValueSTD(lines,n_ini,n_end,"CELL",vcell,std,iv)
      if (iv /= 6) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: Problems reading cell parameters!"
         return
      end if
      if(present(CFrame)) then
         call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, CarType=CFrame,Vscell=std(1:3), Vsang=std(4:6))
      else
         call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, Vscell=std(1:3), Vsang=std(4:6))
      end if
   End Subroutine Read_CFL_Cell_Lines

    Module Subroutine Read_CFL_Cell_FileTyp(cfl, Cell, CFrame)
      !---- Arguments ----!
      Type(File_Type),            intent(in)     :: cfl   ! Containing information
      class(Cell_Type),           intent(out)    :: Cell    ! Cell object
      Character(len=*), optional, intent( in)    :: CFrame
      !---- Local variables -----!
      integer                              :: i, iv,n_ini,n_end
      real(kind=cp), dimension (6)         :: vcell, std
      character(len=132),dimension(1)      :: lines

      !> Init
      call clear_error()

      n_ini=1; n_end=1
      lines(1)=" "
      do i=1,cfl%nlines
        lines(1)=adjustl(u_case(cfl%line(i)%str))
        if(lines(1)(1:4) == "CELL") exit
      end do
      if(len_trim(lines(1)) == 0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: Instruction 'CELL' not provided "
         return
      end if
      call Read_Key_ValueSTD(lines,n_ini,n_end,"CELL",vcell,std,iv)
      if (iv /= 6) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: Problems reading cell parameters!"
         return
      end if
      if(present(CFrame)) then
         call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, CarType=CFrame,Vscell=std(1:3), Vsang=std(4:6))
      else
         call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, Vscell=std(1:3), Vsang=std(4:6))
      end if
   End Subroutine Read_CFL_Cell_FileTyp

   !!----
    !!---- Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
    !!----   character(len=*), dimension(:), intent( in) :: file_dat     !Lines of text (content of a file)
    !!----   integer,                        intent( in) :: i_ini,i_end  !Lines to explore
    !!----   type(job_info_type),            intent(out) :: Job_info     !Object to be constructed here
    !!----
    !!----
    !!----    Constructor of the object Job_info. The arrary of strings file_dat
    !!----    have to be provided as input. It contains lines corresponding to the
    !!----    input control file. The analysis of the command lines is not given here.
    !!----
    !!---- Update: February - 2005, February -2020
    !!
    Module Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
       !---- Arguments ----!
       character(len=*), dimension(:), intent( in) :: file_dat
       integer,                        intent( in) :: i_ini,i_end
       type(job_info_type),            intent(out) :: Job_info

       !---- Local Variables ----!
       integer                           :: i,nphas, ncmd,n_pat,ier, j
       integer, dimension(i_end-i_ini+1) :: ip,ic,ipt
       real(kind=cp)                     :: a1,a2,a3,a4,a5
       character(len=120)                :: line, fmtfields, fmtformat

       !--- Initialize FindFMT
       call Init_FindFMT(i_ini)
       nphas=0
       ncmd=0
       n_pat=0
       ip=i_end
       ic=0
       ipt=0
       Job_info%title=" General Job: CrysFML"
       Job_info%Num_Patterns=1

       do i=i_ini,i_end
          line=u_case(adjustl(file_dat(i)))
          if (line(1:5) == "TITLE") Job_info%title=line(7:)
          if (line(1:5) == "NPATT") then
             read(unit=line(7:), fmt=*,iostat=ier) Job_info%Num_Patterns
             if (ier /= 0) Job_info%Num_Patterns=1
          end if
          if (line(1:6) == "PHASE_") then
             nphas=nphas+1
             ip(nphas)=i
          end if
          if (line(1:4) == "CMDL") then
             ncmd=ncmd+1
             ic(ncmd)=i
          end if
          if (line(1:5) == "PATT_") then
             n_pat=n_pat+1
             ipt(n_pat)=i
          end if
       end do

       if (nphas == 0) then
          nphas=1
          ip(nphas)=0
       end if
       if (n_pat == 0) then
          n_pat=1
          ipt(n_pat) = 0
       end if

       if (Job_info%Num_Patterns /= n_pat) Job_info%Num_Patterns = n_pat
       Job_info%Num_Phases=nphas
       Job_info%Num_Cmd=ncmd

       if (allocated(Job_Info%Patt_typ)) deallocate(Job_Info%Patt_typ)
       allocate(Job_Info%Patt_typ(n_pat))

       if (allocated(Job_Info%Phas_nam)) deallocate(Job_Info%Phas_nam)
       allocate(Job_Info%Phas_nam(nphas))

       if (allocated(Job_Info%range_stl)) deallocate(Job_Info%range_stl)
       allocate(Job_Info%range_stl(n_pat))

       if (allocated(Job_Info%range_q)) deallocate(Job_Info%range_q)
       allocate(Job_Info%range_q(n_pat))

       if (allocated(Job_Info%range_d)) deallocate(Job_Info%range_d)
       allocate(Job_Info%range_d(n_pat))

       if (allocated(Job_Info%range_2theta)) deallocate(Job_Info%range_2theta)
       allocate(Job_Info%range_2theta(n_pat))

       if (allocated(Job_Info%range_energy)) deallocate(Job_Info%range_energy)
       allocate(Job_Info%range_energy(n_pat))

       if (allocated(Job_Info%range_tof)) deallocate(Job_Info%range_tof)
       allocate(Job_Info%range_tof(n_pat))

       if (allocated(Job_Info%lambda)) deallocate(Job_Info%lambda)
       allocate(Job_Info%lambda(n_pat))

       if (allocated(Job_Info%ratio)) deallocate(Job_Info%ratio)
       allocate(Job_Info%ratio(n_pat))

       if (allocated(Job_Info%dtt1)) deallocate(Job_Info%dtt1)
       allocate(Job_Info%dtt1(n_pat))

       if (allocated(Job_Info%dtt2)) deallocate(Job_Info%dtt2)
       allocate(Job_Info%dtt2(n_pat))

       !---- Initialize all variables
       Job_Info%Patt_typ    =" "
       Job_Info%Phas_nam    =" "
       Job_Info%range_stl%mina=0.0
       Job_Info%range_stl%maxb=0.0
       Job_Info%range_q%mina=0.0
       Job_Info%range_q%maxb=0.0
       Job_Info%range_d%mina=0.0
       Job_Info%range_d%maxb=0.0
       Job_Info%range_2theta%mina=0.0
       Job_Info%range_2theta%maxb=0.0
       Job_Info%range_Energy%mina=0.0
       Job_Info%range_Energy%maxb=0.0
       Job_Info%range_tof%mina=0.0
       Job_Info%range_tof%maxb=0.0
       Job_Info%Lambda%mina=0.0
       Job_Info%Lambda%maxb=0.0
       Job_Info%ratio = 0.0
       Job_Info%dtt1 = 0.0
       Job_Info%dtt2 = 0.0
       if (ncmd > 0) then
          if (allocated(Job_Info%cmd)) deallocate(Job_Info%cmd)
          allocate(Job_Info%cmd(ncmd))
          Job_Info%cmd=" "
       end if

       !---- Fill the different fields of Job_Info
       !---- Start with patterns
       fmtfields = "9fffff"

       !---- First asks if there is a PATT_ card, if not a standard is taken
       if (ipt(1) /= 0) then
          do n_pat=1, Job_info%Num_Patterns
             i=ipt(n_pat)
             line=u_case(adjustl(file_dat(i)))
             line=line(8:)
             call findfmt(0,line,fmtfields,fmtformat)
             read(unit=line,fmt=fmtformat) Job_Info%Patt_typ(n_pat), a1,a2,a3,a4,a5
             if (Err_CFML%Ierr /= 0) return
             line=u_case(Job_Info%Patt_typ(n_pat))

             select case(line(1:9))
                case("XRAY_2THE","NEUT_2THE","XRAY_SXTA","NEUT_SXTA")
                   if ( a1 <= 0.000001) a1=1.5405
                   if ( a2 <= 0.000001) then
                      a2=a1
                      a3=0.0
                   end if
                   if (a5 <= a4) a5=120.0
                   Job_Info%Lambda(n_pat)%mina=a1
                   Job_Info%Lambda(n_pat)%maxb=a2
                   Job_Info%ratio(n_pat)=a3
                   Job_Info%range_2theta(n_pat)%mina=a4
                   Job_Info%range_2theta(n_pat)%maxb=a5
                   a4=sind(0.5*a4)/a1
                   a5=sind(0.5*a5)/a2
                   Job_Info%range_stl(n_pat)%mina=a4
                   Job_Info%range_stl(n_pat)%maxb=a5
                   Job_Info%range_q(n_pat)%mina=a4*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
                   Job_Info%range_d(n_pat)%mina=0.5/a5
                   Job_Info%range_d(n_pat)%maxb=0.5/a4

                case("NEUT_TOF ")
                   if (a1 <= 0.000001) a1=1000.0
                   if (a4 <= a3) a4=2.0*abs(a3)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=a2
                   Job_Info%range_tof(n_pat)%mina=a3
                   Job_Info%range_tof(n_pat)%maxb=a4
                   Job_Info%range_d(n_pat)%mina=0.5*(-1.0+sqrt(1.0+4.0*a2*a3/a1/a1))
                   Job_Info%range_d(n_pat)%maxb=0.5*(-1.0+sqrt(1.0+4.0*a2*a4/a1/a1))
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

                case("XRAY_ENER")
                   if (a1 <= 0.000001) a1=12.4 !(=hc(keV.Angstr.)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=0.0
                   Job_Info%range_energy(n_pat)%mina=a3
                   Job_Info%range_energy(n_pat)%maxb=a4
                   if (a3 <= 0.00001) a3=0.01
                   if (a4 <= 0.00001) a4=2.00
                   Job_Info%range_d(n_pat)%mina=a1/a4
                   Job_Info%range_d(n_pat)%maxb=a1/a3
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

             end select
          end do

       else
          n_pat=1
          a1=1.5405
          a2=a1
          a3=0.0
          a4=0.0
          a5=120.0
          Job_Info%Patt_typ(n_pat)="XRAY_2THE"
          Job_Info%Lambda(n_pat)%mina=a1
          Job_Info%Lambda(n_pat)%maxb=a2
          Job_Info%ratio(n_pat)=a3
          Job_Info%range_2theta(n_pat)%mina=a4
          Job_Info%range_2theta(n_pat)%maxb=a5
          a4=sind(0.5*a4)/a1
          a5=sind(0.5*a5)/a2
          Job_Info%range_stl(n_pat)%mina=a4
          Job_Info%range_stl(n_pat)%maxb=a5
          Job_Info%range_q(n_pat)%mina=a4*4.0*pi
          Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
          Job_Info%range_d(n_pat)%mina=0.5/a5
          Job_Info%range_d(n_pat)%maxb=0.5/a4
       end if

       !---- Phase names
       if (ip(1) /= 0) then
          do i=1,nphas
             j=ip(i)
             line=adjustl(file_dat(j))
             Job_Info%Phas_nam(i)=line(8:)
          end do
       else
          Job_Info%Phas_nam(1)= Job_info%title
       end if

       !---- Command Lines, stored but not analysed here
       do i=1,ncmd
          j=ic(i)
          line=adjustl(file_dat(j))
          Job_Info%cmd(i)=line(8:)
       end do

       return
    End Subroutine Get_Job_Info

    Module Subroutine Read_CFL_SpG_lines(lines, n_ini, n_end, SpG, xyz_type)
       !---- Arguments ----!
       character(len=*), dimension(:),  intent(in)     :: lines
       integer,                         intent(in out) :: n_ini
       integer,                         intent(in)     :: n_end
       class(SpG_Type),                 intent(out)    :: SpG
       character(len=*), optional,      intent(in)     :: xyz_type

       !--- Local Variables ---!
       integer :: i,j,k,ngen,nk,nq
       character(len=:),     allocatable :: line,uline,setting,strcode
       character(len=40), dimension(192) :: gen
       logical :: change_setting
       integer :: ier

       call clear_error()
       !Look for the appropriate keywords to construct the space group: Crystallographic, Shubnikov, or superspace
       ngen=0
       setting=" "
       change_setting=.false.
       strcode="xyz"
       if(present(xyz_type)) strcode=xyz_type
       !write(*,"(a,2i5)") " n_ini,n_end:",n_ini,n_end
       do i=n_ini,n_end
         line=trim(adjustl(lines(i)))
         if(len_trim(line) == 0) cycle
         if(line(1:1) == "!" .or. line(1:1) == "#") cycle
         j=index(line,"!")
         if(j /= 0) line=line(1:j-1)
         j=index(line,"::")
         if(j /= 0) then
           setting=trim(adjustl(line(j+2:)))
           if(len_trim(setting) /= 0) change_setting=.true.
           line=line(1:j-1)
         end if
        ! write(*,"(a)") "line: "//line
         j=index(line," ")
         uline=u_case(line(1:j-1))
         !write(*,"(a)") "uline: "//uline
         line=adjustl(line(j+1:))
         !write(*,"(a/)") "line: "//line
         Select Case(uline)
           Case("HALL","SPGR","SPACEG")
              call Set_SpaceGroup(line,SpG)
           Case("SHUB")
              call Set_SpaceGroup(line,"SHUBN",SpG)
           Case("SSG","SUPER","SSPG")
              call Set_SpaceGroup(line,"SUPER",SpG,strcode)
           Case("GENLIST","GENERATORS","LIST")
              call Set_SpaceGroup(line,SpG)
           Case("GEN","SYMM")
              ngen=ngen+1
              gen(ngen)=line
         End Select
       end do

       if(ngen > 0) then
          !write(*,"(10a)") (trim(gen(i))//";", i=1,ngen)
           call Set_SpaceGroup("  ",SpG,ngen,gen)
       end if

       if(Err_CFML%Ierr == 1) return
       if(change_setting) then
         if(strcode == "xyz")  then
           call Change_Setting_SpaceG(setting, SpG)
         else
           call Change_Setting_SpaceG(setting, SpG,strcode)
         end if
       end if
       if(Err_CFML%Ierr == 1) return

       !Now read q-vectors and other items if the class of SpG is SuperSpaceGroup_Type
       Select Type (SpG)
         Class is (SuperSpaceGroup_Type)
           if(allocated(SpG%kv)) deallocate (SpG%kv)
           if(allocated(SpG%nharm)) deallocate (SpG%nharm)
           if(allocated(SpG%sintlim)) deallocate (SpG%sintlim)
           if(allocated(SpG%Om)) deallocate (SpG%Om)
           if(allocated(SpG%q_coeff)) deallocate (SpG%q_coeff)
           allocate(SpG%Om(SpG%D,Spg%D,SpG%Multip))
           do i=1,SpG%Multip
              SpG%Om(:,:,i)=real(SpG%Op(i)%Mat(:,:))
           end do
           nk=0; nq=0
           do i=n_ini,n_end
             line=trim(adjustl(lines(i)))
             if(len_trim(line) == 0) cycle
             if(line(1:1) == "!" .or. line(1:1) == "#") cycle
             j=index(line,"!")
             if(j /= 0) line=line(1:j-1)
             j=index(line," ")
             uline=u_case(line(1:j-1))
             line=adjustl(line(j+1:))
             Select Case(uline)
               Case("NQVECT","NKVECT")
                  Read(unit=line,fmt=*,iostat=ier) Spg%nk, Spg%nq
                  if(ier /= 0) then
                    Err_CFML%Ierr=1
                    Err_CFML%Msg="Error reading the number of k-vectors and/or number of Q-coefficients"
                    return
                  end if
                  allocate(Spg%kv(3,Spg%nk),Spg%q_coeff(Spg%nk,Spg%nq))
                  allocate(Spg%nharm(Spg%nk),Spg%sintlim(Spg%nk))
                  SpG%kv=0.0_cp; SpG%q_coeff=1; Spg%nharm=1; Spg%sintlim=1.0
               Case("QVECT","KVECT")
                  if(Spg%nk > 0) then
                    nk=nk+1
                    Read(unit=line,fmt=*,iostat=ier) Spg%kv(:,nk)
                    if(ier /= 0) then
                      Err_CFML%Ierr=1
                      write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the k-vector #",nk
                      return
                    end if
                  end if
               Case("NHARM")
                  if(Spg%nk > 0) then
                    Read(unit=line,fmt=*,iostat=ier) Spg%nharm(1:Spg%nk)
                    if(ier /= 0) then
                      Err_CFML%Ierr=1
                      Err_CFML%Msg = "Error reading the nk harmonics !"
                      return
                    end if
                  end if
               Case("SINTL")
                  if(Spg%nk > 0) then
                    Read(unit=line,fmt=*,iostat=ier) Spg%sintlim(1:Spg%nk)
                    if(ier /= 0) then
                      Err_CFML%Ierr=1
                      Err_CFML%Msg = "Error reading the maximum sinTheta/Lambda for harmonics!"
                      return
                    end if
                  end if
               Case("Q_COEFF")
                  nq=nq+1
                  Read(unit=line,fmt=*,iostat=ier) Spg%q_coeff(1:Spg%nk,nq)
                  if(ier /= 0) then
                    Err_CFML%Ierr=1
                    write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the Q-coefficent # ",nq
                    return
                  end if

             End Select
           end do
           if(Spg%nk /= (Spg%D-4)) then
              Err_CFML%Ierr=1
              write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of k-vectors,",Spg%nk, ", does not correspond with the additional dimensions of the group ",Spg%D-4
              return
           end if
           if(Spg%nq /= nq) then
              Err_CFML%Ierr=1
              write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of expected Q-coefficients,",Spg%nq, ", does not correspond with number of read Q-coefficients ",nq
              return
           end if

       End Select

    End Subroutine Read_CFL_SpG_lines

    Module Subroutine Read_CFL_SpG_FileTyp(cfl, SpG, xyz_type)
       !---- Arguments ----!
       Type(File_Type),                 intent(in)     :: cfl
       class(SpG_Type),                 intent(out)    :: SpG
       character(len=*), optional,      intent(in)     :: xyz_type

       !--- Local Variables ---!
       integer :: i,j,k,ngen,nk,nq
       character(len=:),     allocatable :: line,uline,setting,strcode
       character(len=40), dimension(192) :: gen
       logical :: change_setting
       integer :: ier

       call clear_error()
       !Look for the appropriate keywords to construct the space group: Crystallographic, Shubnikov, or superspace
       ngen=0
       setting=" "
       change_setting=.false.
       strcode="xyz"
       if(present(xyz_type)) strcode=xyz_type
       !write(*,"(a,2i5)") " n_ini,n_end:",n_ini,n_end
       do i=1,cfl%nlines
         line=adjustl(cfl%line(i)%str)
         if(len_trim(line) == 0) cycle
         if(line(1:1) == "!" .or. line(1:1) == "#") cycle
         j=index(line,"!")
         if(j /= 0) line=line(1:j-1)
         j=index(line,"::")
         if(j /= 0) then
           setting=trim(adjustl(line(j+2:)))
           if(len_trim(setting) /= 0) change_setting=.true.
           line=line(1:j-1)
         end if
        ! write(*,"(a)") "line: "//line
         j=index(line," ")
         uline=u_case(line(1:j-1))
         !write(*,"(a)") "uline: "//uline
         line=adjustl(line(j+1:))
         !write(*,"(a/)") "line: "//line
         Select Case(uline)
           Case("HALL","SPGR","SPACEG")
              call Set_SpaceGroup(line,SpG)
           Case("SHUB")
              call Set_SpaceGroup(line,"SHUBN",SpG)
           Case("SSG","SUPER","SSPG")
              call Set_SpaceGroup(line,"SUPER",SpG,strcode)
           Case("GENLIST","GENERATORS","LIST")
              call Set_SpaceGroup(line,SpG)
           Case("GEN","SYMM")
              ngen=ngen+1
              gen(ngen)=line
         End Select
       end do

       if(ngen > 0) then
          !write(*,"(10a)") (trim(gen(i))//";", i=1,ngen)
           call Set_SpaceGroup("  ",SpG,ngen,gen)
       end if

       if(Err_CFML%Ierr == 1) return
       if(change_setting) then
         if(strcode == "xyz")  then
           call Change_Setting_SpaceG(setting, SpG)
         else
           call Change_Setting_SpaceG(setting, SpG,strcode)
         end if
       end if
       if(Err_CFML%Ierr == 1) return

       !Now read q-vectors and other items if the class of SpG is SuperSpaceGroup_Type
       Select Type (SpG)
         Class is (SuperSpaceGroup_Type)
           if(allocated(SpG%kv)) deallocate (SpG%kv)
           if(allocated(SpG%nharm)) deallocate (SpG%nharm)
           if(allocated(SpG%sintlim)) deallocate (SpG%sintlim)
           if(allocated(SpG%Om)) deallocate (SpG%Om)
           if(allocated(SpG%q_coeff)) deallocate (SpG%q_coeff)
           allocate(SpG%Om(SpG%D,Spg%D,SpG%Multip))
           do i=1,SpG%Multip
              SpG%Om(:,:,i)=real(SpG%Op(i)%Mat(:,:))
           end do
           nk=0; nq=0
           do i=1,cfl%nlines
             line=adjustl(cfl%line(i)%str)
             if(len_trim(line) == 0) cycle
             if(line(1:1) == "!" .or. line(1:1) == "#") cycle
             j=index(line,"!")
             if(j /= 0) line=line(1:j-1)
             j=index(line," ")
             uline=u_case(line(1:j-1))
             line=adjustl(line(j+1:))
             Select Case(uline)
               Case("NQVECT","NKVECT")
                  Read(unit=line,fmt=*,iostat=ier) Spg%nk, Spg%nq
                  if(ier /= 0) then
                    Err_CFML%Ierr=1
                    Err_CFML%Msg="Error reading the number of k-vectors and/or number of Q-coefficients"
                    return
                  end if
                  allocate(Spg%kv(3,Spg%nk),Spg%q_coeff(Spg%nk,Spg%nq))
                  allocate(Spg%nharm(Spg%nk),Spg%sintlim(Spg%nk))
                  SpG%kv=0.0_cp; SpG%q_coeff=1; Spg%nharm=1; Spg%sintlim=1.0
               Case("QVECT","KVECT")
                  if(Spg%nk > 0) then
                    nk=nk+1
                    Read(unit=line,fmt=*,iostat=ier) Spg%kv(:,nk)
                    if(ier /= 0) then
                      Err_CFML%Ierr=1
                      write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the k-vector #",nk
                      return
                    end if
                  end if
               Case("NHARM")
                  if(Spg%nk > 0) then
                    Read(unit=line,fmt=*,iostat=ier) Spg%nharm(1:Spg%nk)
                    if(ier /= 0) then
                      Err_CFML%Ierr=1
                      Err_CFML%Msg = "Error reading the nk harmonics !"
                      return
                    end if
                  end if
               Case("SINTL")
                  if(Spg%nk > 0) then
                    Read(unit=line,fmt=*,iostat=ier) Spg%sintlim(1:Spg%nk)
                    if(ier /= 0) then
                      Err_CFML%Ierr=1
                      Err_CFML%Msg = "Error reading the maximum sinTheta/Lambda for harmonics!"
                      return
                    end if
                  end if
               Case("Q_COEFF")
                  nq=nq+1
                  Read(unit=line,fmt=*,iostat=ier) Spg%q_coeff(1:Spg%nk,nq)
                  if(ier /= 0) then
                    Err_CFML%Ierr=1
                    write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the Q-coefficent # ",nq
                    return
                  end if

             End Select
           end do
           if(Spg%nk /= (Spg%D-4)) then
              Err_CFML%Ierr=1
              write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of k-vectors,",Spg%nk, ", does not correspond with the additional dimensions of the group ",Spg%D-4
              return
           end if
           if(Spg%nq /= nq) then
              Err_CFML%Ierr=1
              write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of expected Q-coefficients,",Spg%nq, ", does not correspond with number of read Q-coefficients ",nq
              return
           end if

       End Select

    End Subroutine Read_CFL_SpG_FileTyp

    Module Subroutine Read_kinfo(cfl, kinfo)
       !---- Arguments ----!
       type(File_Type),                 intent(in)     :: cfl     ! File_type object
       type(kvect_info_Type),           intent(out)    :: kinfo
       !
       integer :: i,j,ier,nk,nq
       character(len=:),allocatable :: uline,line

       nk=0; nq=0
       do i=1,cfl%nlines
         line=adjustl(cfl%line(i)%str)
         if(len_trim(line) == 0) cycle
         if(line(1:1) == "!" .or. line(1:1) == "#") cycle
         j=index(line,"!")
         if(j /= 0) line=line(1:j-1)
         j=index(line," ")
         if( j == 0) then
           uline=u_case(line)
         else
           uline=u_case(line(1:j-1))
         end if
         line=adjustl(line(j+1:))

         Select Case(uline)

           Case("NQVECT","NKVECT","NKVEC")
              Read(unit=line,fmt=*,iostat=ier) kinfo%nk, kinfo%nq
              if(ier /= 0) then
                Err_CFML%Ierr=1
                Err_CFML%Msg="Error reading the number of k-vectors and/or number of Q-coefficients"
                return
              end if
              allocate(kinfo%kv(3,kinfo%nk),kinfo%q_coeff(kinfo%nk,kinfo%nq))
              allocate(kinfo%nharm(kinfo%nk),kinfo%sintlim(kinfo%nk))
              kinfo%kv=0.0_cp; kinfo%q_coeff=1; kinfo%nharm=1; kinfo%sintlim=1.0

           Case("QVECT","KVECT","KVEC")
              if(kinfo%nk > 0) then
                nk=nk+1
                Read(unit=line,fmt=*,iostat=ier) kinfo%kv(:,nk)
                if(ier /= 0) then
                  Err_CFML%Ierr=1
                  write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the k-vector #",nk
                  return
                end if
              end if

           Case("NHARM")
              if(kinfo%nk > 0) then
                Read(unit=line,fmt=*,iostat=ier) kinfo%nharm(1:kinfo%nk)
                if(ier /= 0) then
                  Err_CFML%Ierr=1
                  Err_CFML%Msg = "Error reading the nk harmonics !"
                  return
                end if
              end if

           Case("SINTL")
              if(kinfo%nk > 0) then
                Read(unit=line,fmt=*,iostat=ier) kinfo%sintlim(1:kinfo%nk)
                if(ier /= 0) then
                  Err_CFML%Ierr=1
                  Err_CFML%Msg = "Error reading the maximum sinTheta/Lambda for harmonics!"
                  return
                end if
              end if

           Case("Q_COEFF")
              nq=nq+1
              Read(unit=line,fmt=*,iostat=ier) kinfo%q_coeff(1:kinfo%nk,nq)
              if(ier /= 0) then
                Err_CFML%Ierr=1
                write(unit=Err_CFML%Msg,fmt="(a,i2)") "Error reading the Q-coefficent # ",nq
                return
              end if
         End Select
       end do

       if(kinfo%nk /= nk) then
          Err_CFML%Ierr=1
          write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of k-vectors,",kinfo%nk, ", does not correspond with the prescribed number: ",nk
          return
       end if

       if(kinfo%nq /= nq) then
          Err_CFML%Ierr=1
          write(unit=Err_CFML%Msg,fmt="(2(a,i2))") "The number of expected Q-coefficients,",kinfo%nq, ", does not correspond with number of read Q-coefficients ",nq
          return
       end if

     End Subroutine Read_kinfo
End SubModule IOF_CFL