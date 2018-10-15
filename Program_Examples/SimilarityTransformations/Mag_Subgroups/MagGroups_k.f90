 Module transformations
    use CFML_GlobalDeps
    Use CFML_Math_General,     only: Modulo_Lat, equal_matrix, Zbelong
    Use CFML_Math_3D,          only: Determ_A, matrix_inverse
    Use CFML_Crystallographic_Symmetry, only: Set_SpaceGroup, Space_Group_Type, Get_Crystal_System, &
                                              Get_SymSymb, Get_Centring_Vectors, SpGr_Equal,Get_Hallsymb_From_Gener,&
                                              NS_Space_Group_Type,Ltr,Get_GenSymb_from_Gener, &
                                              Lattice_trans, Copy_NS_SpG_To_SpG, Latsym, inlat
    Use CFML_String_Utilities, only: l_case, pack_string, Frac_Trans_2Dig
    Use CFML_Crystal_Metrics
    Use CFML_Symmetry_Tables, only : Latt, Sys_cry
    implicit none
    private
    public :: Setting_Change_NonConv, get_C_Subgroups
    logical,public     ::      err_trans=.false.
    character(len=120) ::      ERR_transmess= " "
    real(kind=cp), parameter, private :: eps_symm  = 0.0002_cp

 contains

    Subroutine Get_C_SubGroups(SpG,SubG,nsg,point)
       !---- Arguments ----!
       type (Space_Group_Type) ,             intent( in) :: SpG
       type (Space_Group_Type) ,dimension(:),intent(out) :: SubG
       integer,                              intent(out) :: nsg
       logical, dimension(:,:), optional,    intent(out) :: point
       !--- Local variables ---!
       integer                            :: i,L,j,k, nc, maxg,ng , nla, i1,i2,nop
       character (len=30), dimension(192) :: gen
       logical                            :: newg, cen_added

       maxg=size(SubG)
       !---- Remove first the generators of centring translations ----!
       ng=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings
       if (SpG%centred /= 1) nop=nop*2
       nla=ng
       nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
       L=0
       !---- Determine first the triclinic subgroups
       cen_added=.false.
       do
           L=L+1
           newg=.true.
           call set_spacegroup(" ",SubG(L),gen,ng,"gen")
           do j=1,L-1
              if (SpGr_Equal(SubG(L), SubG(j))) then
                 newg=.false.
                 exit
              end if
           end do
           if (newg) then
              call get_HallSymb_from_gener(SubG(L))
           else
              L=L-1
           end if
           if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
              ng=ng+1
              gen(ng)=SpG%SymopSymb(nc)
              cen_added=.true.
           else
              exit
           end if
       end do

       !---- Determine first the groups with only one rotational generator
       do i=2,nop
          ng=nla+1
          gen(ng) = SpG%SymopSymb(i)
          cen_added=.false.
          do
             L=L+1
             if (L > maxg) then
                nsg=maxg
                return
             end if
             newg=.true.
             call set_spacegroup(" ",SubG(L),gen,ng,"gen")
             do j=1,L-1
                if (SpGr_Equal(SubG(L), SubG(j))) then
                   newg=.false.
                   exit
                end if
             end do
             if (newg) then
                call get_HallSymb_from_gener(SubG(L))
             else
                L=L-1
             end if
             if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
                ng=ng+1
                gen(ng)=SpG%SymopSymb(nc)
                cen_added=.true.
             else
                exit
             end if
          end do
       end do

       !---- Determine now the groups with two rotational generator ----!

       do i1=2,nop-1
          gen(nla+1) = SpG%SymopSymb(i1)
          do i2 = i1+1,nop
             gen(nla+2) = SpG%SymopSymb(i2)
             ng=nla+2
             cen_added=.false.
             do
                L=L+1
                if (L > maxg) then
                   nsg=maxg
                   return
                end if
                newg=.true.
                call set_spacegroup(" ",SubG(L),gen,ng,"gen")
                do j=1,L-1
                   if (SpGr_Equal(SubG(L), SubG(j))) then
                      newg=.false.
                      exit
                   end if
                end do
                if (newg) then
                   call get_HallSymb_from_gener(SubG(L))
                else
                   L=L-1
                end if
                if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
                   ng=ng+1
                   gen(ng)=SpG%SymopSymb(nc)
                   cen_added=.true.
                else
                   exit
                end if
             end do
          end do
       end do
       nsg=L
       if(present(point)) then
         point=.false.
         do j=1,nsg
           L=1
           do i=1,SpG%multip
              do k=L,SubG(j)%multip
               if(SubG(j)%SymopSymb(k) == SpG%SymopSymb(i)) then
                  point(i,j) = .true.
                  L=k+1
                  exit
               end if
              end do
           end do
         end do
       end if

       return
    End Subroutine Get_C_SubGroups


     Subroutine Setting_Change_NonConv(Mat,Orig,Spg,Spgn,Matkind,purge)
       !---- Arguments ----!
       real(kind=cp), dimension(3,3),intent(in) :: Mat
       real(kind=cp), dimension(3),  intent(in) :: Orig
       type (Space_Group_Type),      intent(in) :: SpG
       type (NS_Space_Group_Type),  intent(out) :: SpGn
       character (len=*), optional,  intent(in) :: Matkind
       logical,           optional,  intent(in) :: purge

       !--- Local variables ---!
       integer                 :: ifail, i, j, k, L, im, nc, m, ngm,n,ngen,nlat,isyst
       real(kind=cp)           :: det,rmin,rmax
       character(len=40)       :: transla
       character(len=1)        :: LatSymb,Cryst
       real(kind=cp), dimension (3,3), parameter :: e = reshape ((/1.0,0.0,0.0,  &
                                                                   0.0,1.0,0.0,  &
                                                                   0.0,0.0,1.0/),(/3,3/))
       real(kind=cp), dimension (3,192)    :: newlat = 0.0 !big enough number of centring tranlations
       real(kind=cp), dimension (3,3)      :: S, Sinv, rot, rotn  !S is the ITC matrix P.
       integer,       dimension (3,3)      :: nulo,irot
       real(kind=cp), dimension (  3)      :: tr, trn, v
       logical                             :: lattl,change_only_origin
       character(len=80)                   :: symbsg
       character(len=60),dimension(15)     :: gen
       character(len=180)                  :: setting
       real(kind=cp),  dimension(3,3,Spg%Multip) :: sm
       real(kind=cp),  dimension(3,Spg%Multip)   :: tm

       !call Init_Err_Symm()
       change_only_origin=.false.
       nulo=0
       !---- Up to here all "conventionnal" translational generators have been obtained
       !---- Set the minimum and maximum admissible component of translations
       select case (SpG%CrystalSys)
          Case("Triclinic")
             rmin=0.0
             rmax=1.0
          Case("Monoclinic")
             rmin=0.5
             rmax=0.5
          Case("Orthorhombic")
             rmin=0.5
             rmax=0.5
             if (SpG%SPG_lat == "F") then
                rmin=0.25
                rmax=0.75
             end if
          Case("Tetragonal")
             rmin=0.25
             rmax=0.75
          Case("Rhombohedral","Hexagonal","Trigonal")
             rmin=1.0/6.0
             rmax=5.0/6.0
          Case("Cubic")
             rmin=0.25
             rmax=0.75
          Case default
             rmin=0.5
             rmax=0.5
       end select

       call get_setting_info(Mat,orig,setting,matkind)
       symbsg=Pack_String(SpG%spg_symb)
       if (present(matkind)) then
          if (matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
             S=Mat
          else
             S=transpose(Mat)
          end if
       else
          S=transpose(Mat)
       end if
       setting = trim(setting)//" det:"
       if(equal_matrix(S,e,3)) change_only_origin=.true.
       det=determ_a(Mat)
       i=len_trim(setting)
       write(unit=setting(i+2:),fmt="(f6.2)") det
       !write(unit=*,fmt="(a)") " => Setting Symbol: "//trim(setting)
       call matrix_inverse(S,Sinv,ifail)
       if (ifail /= 0) then
          err_trans=.true.
          ERR_transmess= "Inversion Matrix Failed on: Setting_Change_NonConv"
          return
       end if

       L=0
       if (SpG%NumLat > 1) then  !Original lattice is centered
          do i=2,SpG%NumLat      !Transform the centring vectors to the new lattice
             v=Modulo_Lat(matmul(Sinv,SpG%Latt_trans(:,i)))
             if (sum(v) < eps_symm) cycle
             L=L+1
             newlat(:,L)=v
          end do
       end if
       do i=1,3  !Test the basis vectors of the original setting
         rot(:,i)=Modulo_Lat(Sinv(:,i))
         if (sum(rot(:,i)) < eps_symm) cycle
         L=L+1
         newlat(:,L)=rot(:,i)
       end do

       if (det > 1 ) then  !The new lattice is centred
          im=nint(det)-1         !Determine the new lattice translations
          ngm=L+im
          doi: do i=0,im
             v(1) = i
             do j=0,im
                v(2) = j
                do k=0,im
                   v(3) = k
                   if (nint(sum(v)) == 0) cycle
                   tr=Modulo_Lat(matmul(Sinv,v))
                   if (sum(tr) < eps_symm) cycle
                   lattl =.true.
                   do m=1,L
                      if (sum(abs(tr-newlat(:,m))) < eps_symm) then
                         lattl =.false.
                         exit
                      end if
                   end do
                   if (lattl) then ! new lattice translation
                      L=L+1
                      newlat(:,L) = tr(:)
                      if (L == ngm) exit doi
                   end if
                end do !k
             end do !j
          end do doi !i
       end if

       call get_centring_vectors(L,newlat,LatSymb)  !Complete the centring vectors
       nlat=L+1
       !Now we have L centring translations
       call LatSym(LatSymb,L,newlat)  !provides the value of the global variable inlat: index of the type of lattice
       SpGn%SPG_lat      = LatSymb

       !---- Change of symmetry operator under a change of basis and origin
       !----  A'= M A,  origin O =>  X'=inv(Mt)(X-O)
       !----  Symmetry operator C = (R,T)  -> C' = (R',T')
       !----   R' = inv(Mt) R Mt                 ITC:    R'= inv(P) R P
       !----   T' = inv(Mt) (T -(E-R)O)                  T'= inv(P) (T-(E-R)O)
       sm=0.0
       tm=0.0
       sm(:,:,1)=SpG%SymOp(1)%Rot
       tm(:,1)=SpG%SymOp(1)%tr
       n=1
       do_i:do i=2,SpG%NumOps
          Rot=SpG%SymOp(i)%rot
          Rotn=matmul(matmul(Sinv,Rot),S)
          do k=n,1,-1
            if(equal_matrix(Rotn,sm(:,:,k),3))  cycle do_i
          end do
          tr=SpG%SymOp(i)%tr
          trn=matmul(Sinv,tr-matmul(e-Rot,orig))
          if(present(purge)) then
            irot=abs(nint(rotn))
            if ( any(irot > 1) ) cycle    !Conserve only the conventional forms  |aij|=1,0
            if (.not. Zbelong(Rotn)) cycle
            if ( any((trn < rmin .and. trn > 0.0) .or. trn > rmax) ) cycle
          end if
          n=n+1
          sm(:,:,N)=Rotn
          tm(:,n)=Modulo_Lat(trn)
       end do do_i

       SpGn%Centred=SpG%Centred
       SpGn%Centre_coord=SpG%Centre_coord
       if (SpG%Centred /= 1) then !the space group is centro-symmetric
          nc=SpG%NumOps+1
          Rot=SpG%SymOp(nc)%rot
          tr=SpG%SymOp(nc)%tr
          trn=matmul(Sinv,tr-matmul(e-Rot,orig)) ! matmul(Sinv,tr-2*orig)
          trn= Modulo_Lat(trn)
          if(sum(abs(trn)) > 3.0*eps_symm) then
            SpGn%Centred=0
            SpGn%Centre_coord=0.5*trn
          end if
       end if

       !Do another thing we conserve the transformations and generate ourself the new group
       !The new multiplicity is n * i * nlat
       i=1
       if(SpGn%Centred /= 1) i=2
       SpGn%multip= n * i * nlat  !nlat=L+1

       allocate(SpGn%SymOp(SpGn%multip), SpGn%SymOpSymb(SpGn%multip))
       SpGn%NumOps=n
       do i=1,SpGn%NumOps
         SpGn%SymOp(i)%Rot=sm(:,:,i)
         SpGn%SymOp(i)%tr=tm(:,i)
       end do

       allocate(SpGn%Latt_trans(3,nlat))
       SpGn%NumLat    = nlat
       SpGn%Latt_trans= Ltr(:,1:nlat)
       if(purge) then
         call Get_Crystal_System(n,nint(sm),isyst,Cryst)
         SpGn%CrystalSys   = Sys_cry(isyst)
         select case (SpGn%CrystalSys)
            Case("Triclinic")
               SpGn%SPG_latsy    = "a"//LatSymb
            Case("Monoclinic")
               SpGn%SPG_latsy    = "m"//LatSymb
            Case("Orthorhombic")
               SpGn%SPG_latsy    = "o"//LatSymb
            Case("Tetragonal")
               SpGn%SPG_latsy    = "t"//LatSymb
            Case("Rhombohedral")
               SpGn%SPG_latsy    = "h"//LatSymb
            Case("Hexagonal")
               SpGn%SPG_latsy    = "h"//LatSymb
            Case("Trigonal")
               SpGn%SPG_latsy    = "h"//LatSymb
            Case("Cubic")
               SpGn%SPG_latsy    = "c"//LatSymb
            Case default
               SpGn%SPG_latsy    = "a"//LatSymb
         end select
       else
         SpGn%CrystalSys   = SpG%CrystalSys
         SpGn%NumSpg=SpG%NumSpg
         SpGn%SPG_latsy    = SpG%SPG_latsy(1:1)//LatSymb
       end if

       SpGn%PG           = SpG%PG
       SpGn%Laue         = SpG%laue
       SpGn%Num_gen      = SpG%Num_gen
       SpGn%SG_setting   = setting
       SpGn%Bravais      = Latt(inlat)

       Select Case (SpGn%Centred)
           Case(0,2)
             call Frac_Trans_2Dig(SpGn%Centre_coord,transla)
             SpGn%centre="Centric, -1 at "//trim(transla)
           Case Default
             SpGn%centre="Acentric"
       End Select
       m=SpGn%Numops
       if (SpGn%centred /= 1) then
          do i=1,SpGn%Numops
             m=m+1
             SpGn%Symop(m)%Rot(:,:) = -SpGn%Symop(i)%Rot(:,:)
             SpGn%Symop(m)%tr(:)    =  modulo_lat(-SpGn%Symop(i)%tr(:)+2.0*SpGn%Centre_coord)
          end do
       end if
       ngm=m
       if (SpGn%NumLat > 1) then
          do L=2,SpGn%NumLat
             do i=1,ngm
                m=m+1
                trn=SpGn%Symop(i)%tr(:) + SpGn%Latt_trans(:,L)
                SpGn%Symop(m)%Rot(:,:) = SpGn%Symop(i)%Rot(:,:)
                SpGn%Symop(m)%tr(:)    = modulo_lat(trn)
             end do
          end do
       end if
       do i=1,SpGn%multip
          call Get_SymSymb(SpGn%Symop(i)%Rot(:,:), &
                           SpGn%Symop(i)%tr(:)   , &
                           SpGn%SymopSymb(i))
       end do
       !Try to assign a Hall symbol to the space group in the new setting
       !If the hall symbol has been found and the symbol exists in the table the H-M symbol is also set.
       SpGn%hall="From:"//trim(SpG%hall)
       SpGn%spg_symb="From:"//trim(SpG%spg_symb)
       if(change_only_origin) then
         SpGn%spg_symb=trim(symbsg)
       else
         if(SpGn%NumSpg == 0) then
            SpGn%spg_symb="From:"//trim(symbsg)
         end if
       end if
       !Generate a general symbol, first select generators as a function of SpGn$Numops
       n=SpGn%NumOps
       Select Case(SpGn%centred)
         Case(0)
           gen(1)=SpGn%SPG_lat
           gen(2)=SpGn%SymopSymb(n+1)
           ngen=2
         Case(1)
           gen(1)=SpGn%SPG_lat
           ngen=1
         Case(2)
           gen(1)="-"//SpGn%SPG_lat
           ngen=1
       End Select

       Select Case(n)
         case(1:3)
           ngen=ngen+1
           gen(ngen)=SpGn%SymopSymb(2)
         case(4:)
           ngen=ngen+1
           gen(ngen)=SpGn%SymopSymb(2)
           ngen=ngen+1
           gen(ngen)=SpGn%SymopSymb(3)
       End Select
       !
       call Get_GenSymb_from_Gener(gen,ngen,SpGn%ghall)

       return
    End Subroutine Setting_Change_NonConv

    Subroutine Get_Setting_Info(Mat,orig,setting,matkind)
       !---- Arguments ----!
       real(kind=cp), dimension (3,3),intent( in)    :: Mat
       real(kind=cp), dimension (  3),intent( in)    :: orig
       character (len=*),             intent(out)    :: setting
       character (len=*), optional,   intent( in)    :: matkind

       !---- local variables ----!
       real(kind=cp), dimension (  3), parameter  :: nul = (/ 0.0, 0.0, 0.0/)
       real(kind=cp), dimension (3,3)  :: S
       character (len=22)     :: tro
       integer                :: i

       if (present(matkind)) then
          if (matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
             S=transpose(Mat)
          else
             S=Mat
          end if
       else
          S=Mat
       end if

       call Get_SymSymb(S,nul,setting)
       i=index(setting,",")
       setting="a'="//setting(1:i)//" b'="//setting(i+1:)
       i=index(setting,",",back=.true.)
       setting=setting(1:i)//" c'="//setting(i+1:)
       do i=1,len_trim(setting)
          if (setting(i:i) == "x")  setting(i:i) = "a"
          if (setting(i:i) == "y")  setting(i:i) = "b"
          if (setting(i:i) == "z")  setting(i:i) = "c"
       end do

       call Frac_Trans_2Dig(Orig,tro)
       i=len_trim(setting)
       setting(i+2:)=" -> Origin: "//trim(tro)

       return
    End Subroutine Get_Setting_Info

 End Module transformations



 Program MagGroups_k
!!---------------------------------------------------------------------------------
!!--- Program MagGroups_k
!!--- Purpose: Generate maximal magnetic groups corresponding to a conmensurate propagation vector.
!!---          Input: Space Group, cell parameters, atom positions, k-vector
!!---
!!--- (C) Created by JRC April 2015. Based in Similar
!!--- Author: Juan Rodriguez-Carvajal (Institut Laue-Langevin)
!!---------------------------------------------------------------------------------
    use CFML_GlobalDeps
    use transformations
    Use CFML_Math_General,     only: equal_matrix
    Use CFML_Crystallographic_Symmetry, only: Set_SpaceGroup, Write_SpaceGroup, Space_Group_Type, &
                                              get_stabilizer, Get_multip_pos, Get_SymSymb, &
                                              get_orbit,applyso, get_T_SubGroups,NS_Space_Group_Type,&
                                              Lattice_trans, similar_transf_SG, Copy_NS_SpG_To_SpG, &
                                              Setting_Change,Get_Laue_PG
    Use CFML_String_Utilities, only: l_case, number_lines, reading_lines, pack_string
    Use CFML_Crystal_Metrics
    Use CFML_Atom_TypeDef,     only:  Atom_list_Type, Allocate_Atom_list, Write_Atom_List, &
                                      Set_Atom_Equiv_List,Atom_Equiv_List_Type

    Use CFML_Geometry_Calc,    only: point_list_type, get_transf_list, allocate_point_list, &
                                     deallocate_point_list,set_orbits_inlist,Set_New_AsymUnit, &
                                     err_geom_Mess,err_geom
    Use CFML_IO_formats
    Use CFML_Propagation_Vectors

    Use CFML_Math_3D, only: set_eps, determ_a

    Implicit None

    type(file_list_type) :: file_dat
    integer                                 :: nlines, nlong, n_ini,n_end
    type (Crystal_Cell_Type)                :: Cell, Cell_n
    type (Space_Group_Type)                 :: SpaceGroup,SpaceGroup_n, SubGk
    type (NS_Space_Group_Type)              :: SpGn
    type (Group_k_Type)                     :: Gk
    type (Space_Group_Type),dimension (48)  :: SubGroup,SubG_C
    type (Atom_list_Type)                   :: A, A_n     !List of atoms in the asymmetric unit
    type (Atom_Equiv_List_Type)             :: Ate,Ate_n  !List of all atoms in the cell
    type (Point_list_type)                  :: pl, pl_n
    logical, dimension(192,48)              :: pointr,pointc
    character(len=1)         :: ans
    character(len=5)         :: so_ord
    character(len=12)        :: nam
    character(len=20)        :: spp, sppg,symb !symbol of space group
    character(len=80)        :: line , title, cmdline
    character(len=256)       :: filcod,outfil,texto
    integer, parameter       :: lun1=1,lun2=6,lun=2
    integer                  :: i, j, numops, ier, ln, nauas, natc, iid,  nmag, len_cmdline, &
                                lenf, lr, nsg, nat, i1,i2,nsgc
    integer                  :: l,ng, indx, k, order,mulg, norbi, n, m, ifail, indice,i_cfl,mult
    real                     :: seconds, End_time, start_time, rminutes, hours, det,occ
    real,    dimension(3,3)  :: trans,identity      !matrix transforming the cells
    real,    dimension(3  )  :: orig       !origin of the transformed cell in old cell
    real,    dimension(3  )  :: xp         !auxiliary 3D-vector
    real,    dimension(3  )  :: kv,magcel  !Propagation vector and cell multiplication along the axes
    real,    dimension(6  )  :: cel        !cell parameters
    logical                  :: iprin, k_given, full_given, esta
    integer, dimension (192) :: spg_ptr, ptr, norb, sbg_ptr
    real,   dimension(3,192) :: xo !orbits matrix
    character(len=3)         :: tinv
    integer                  :: narg
    integer, dimension(3)    :: kinv

    narg=COMMAND_ARGUMENT_COUNT()
    len_cmdline=0
    if(narg > 0) then
            call GET_COMMAND_ARGUMENT(1,cmdline)
            len_cmdline=len_trim(cmdline)
            outfil=" "
            if(narg > 1) then
              call GET_COMMAND_ARGUMENT(1,outfil)
              i=index(outfil,".")
              if(i /= 0) outfil=outfil(1:i-1)
            end if
    end if

    call set_eps(0.001)   !Needed for atom position comparisons

    write(unit=*,fmt="(/,/,6(a,/))")                                                 &
    "                      ------ PROGRAM MagGroups_k  ------",                         &
    "                       ---- Version 1.0 April-2015 ----",                          &
    "  *******************************************************************************",  &
    "  * Determination of magnetic space groups compatible with a propagation vector *",  &
    "  *******************************************************************************",  &
    "                           (ILL JRC- April 2015 )"

      if(len_cmdline /=0) then
        lenf=index(cmdline,".")-1
        if(lenf <= 0) then
          filcod=cmdline
        else
          filcod=cmdline(1:lenf)
        end if
        if(len_trim(outfil) == 0) outfil=filcod
        ln=len_trim(filcod)
        lr=len_trim(outfil)
      else
        write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.cfl (give xx): "
        read(unit=*,fmt="(a)") filcod
        ln=len_trim(filcod)
        write(unit=*,fmt="(3a)",advance="no") " => Code of the output file (.mag) ( <cr>= ",filcod(1:ln),") :"
        read(unit=*,fmt="(a)") outfil
        lr=len_trim(outfil)
        IF(lr == 0) THEN
          outfil=filcod
          lr=len_trim(outfil)
        END IF
      end if

     inquire(file=trim(filcod)//".cfl",exist=esta)
     if( .not. esta) then
       write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cfl doesn't exist!"
       stop
     end if
     call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpaceGroup,A,Mode="CFL",file_list=file_dat)
     If(err_form) then
       write(unit=*,fmt="(a)") trim(err_form_mess)
       stop
     end if

    open(unit=lun,file=outfil(1:lr)//".magr",status="replace",action="write")

    !Writing titles and content of the input file
    write(unit=lun,fmt="(/,/,6(a,/))")                                                  &
    "                      ------ PROGRAM MagGroups_k  ------",                         &
    "                       ---- Version 1.0 April-2015 ----",                          &
    "  *******************************************************************************",  &
    "  * Determination of magnetic space groups compatible with a propagation vector *",  &
    "  *******************************************************************************",  &
    "                           (ILL JRC- April 2015 )"
    write(unit=lun,fmt="(/,a,a,/)")" => Content of the input file: ",filcod(1:ln)//".cfl"

 !  Second reading and start calculations
 !  -------------------------------------
    CALL CPU_TIME(start_time)

    n_ini=1
    n_end=file_dat%nlines
    orig=0.0

     n_ini=1
     indice=2
     k_given=.false.

    do i=1,n_end
      line=l_case(file_dat%line(i))
      j=index(line,"kv")
      if(j /= 0) then
        line=adjustl(line(j:))
        j=index(line," ")
        read(unit=line(j:),fmt=*) kv
        k_given=.true.
        exit
      end if
    end do
    magcel = 1.0
    identity=reshape([ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 ],[3,3])
    if(k_given) then
       where(abs(kv) > 0.00001)
         kinv=nint(1.0/kv)
       elsewhere
         kinv=1
       end where
       magcel=kinv
       trans=reshape([ magcel(1),0.0,0.0, 0.0,magcel(2),0.0, 0.0,0.0,magcel(3) ],[3,3])
    else
       kv=0.0
       trans=identity
    end if
    ! Print the whole information on cell, space group and atoms
    call Write_SpaceGroup(SpaceGroup,lun,full=.true.)
    call Write_Crystal_Cell(Cell,lun)
    call Write_Atom_List(A,lun=lun)
    !Calculate all atoms within a unit cell -> Ate
    call Set_Atom_Equiv_List(SpaceGroup,cell,A,Ate,lun)
    !Construct and write the propagation vector group
    call K_Star(kv,Spacegroup,Gk)
    call Write_Group_k(Gk,lun)
    !Start working with the propagation vector group, get first the space group
    !from the generators of Gk
    call set_Gk(Gk,SubGk)
    if(SubGk%Spg_Symb /= Spacegroup%Spg_Symb) then
      write(unit=lun,fmt="(/,a,/,a,/,a/)")"    ------------------------------------------------------------------", &
                                          " => Writing the Propagation Vector Group as a conventional Space Group:", &
                                          "    ------------------------------------------------------------------"
      call Write_SpaceGroup(SubGk,lun,full=.true.)
    else
      write(unit=lun,fmt="(/,a,/)")" => The propagation vector group is the same as the original group"
    end if

    if(.not. equal_matrix(identity,trans,3)) then
       write(unit=lun,fmt="(/,a,/,a,/)")" => Change of Space group setting according to the transformation:", &
                           "    --------------------------------------------------------------"
       write(unit=lun,fmt="(a)") "                 Matrix M (A'= M A)              Origin"
       do i=1,3
          write(unit=lun,fmt="(a,(3f8.4,a,f10.5))") "            ",trans(i,1:3), "          ", orig(i)
       end do
       det=determ_a(trans)
       det=abs(det)
       write(unit=lun,fmt="(/,a,/)")" => Maximal subgroup of the original space group according to the cell transformation:"
       !call Similar_Transf_SG(trans,orig,SubGk,SpaceGroup_n)
       call Setting_Change_NonConv(trans,Orig,SubGk,Spgn,purge=.true.)
       call Copy_NS_SpG_To_SpG(Spgn,SpaceGroup_n)
       call Get_Laue_PG(SpaceGroup_n, SpaceGroup_n%Laue, SpaceGroup_n%PG)
       call Write_SpaceGroup(SpaceGroup_n,lun, .true.)
       call Change_Setting_Cell(Cell,trans,Cell_n)

       write(unit=lun,fmt="(/,a,/,a,/)")" => Change of Unit cell according to the transformation:", &
                                        "    ----------------------------------------------------"
       write(unit=lun,fmt="(a)")        "                Matrix M (A'= M A)               Origin"
       do i=1,3
          write(unit=lun,fmt="(a,(3f8.4,a,f10.5))") "            ",trans(i,1:3), "          ", orig(i)
       end do
       call Write_Crystal_Cell(Cell_n,lun)

       ! It is needed to increase the asymmetric unit because some atoms may be lost
       ! after decreasind the number of symmetry operators in the new subgroup.
       ! call change_setting_atoms(Cell_n,A,trans,orig,A_n)  <= not adequate in this context
       write(unit=*,fmt="(a)")" => Creating The new Asymmetric unit"
       call Set_new_AsymUnit(SpaceGroup_n,Ate,trans,orig,A_n,debug="D")
       if(err_geom) then
         write(unit=*,fmt="(a)") trim(err_geom_Mess)
       end if
       ! Write the atoms the in new asymmetric unit
       call Write_Atom_List(A_n,lun=lun)
       write(unit=lun,fmt="(/,a,/,a,/)")" => List of all atoms in the new cell", &
                                        "    ---------------------------------"
       call Set_Atom_Equiv_List(SpaceGroup_n,cell_n,A_n,Ate_n,lun)

       call allocate_point_list(SpaceGroup%multip,pl,ier)  !Allocate the point list for original cell
       if(ier /= 0) then
          write(unit=*,fmt="(a)")  " => Error allocating  ", SpaceGroup%multip," points for PL"
       end if
       if(det >= 1.0) then
         nat=SpaceGroup%multip*det+1
       else
         nat=SpaceGroup%multip+1
       end if
       call Allocate_Point_List(nat,pl_n,ier)   !allocate the new point list in the new cell
       if(ier /= 0) then
          write(unit=*,fmt="(a)")  " => Error allocating  ", nat," points for PL_n"
       end if
    else
       SpaceGroup_n=SubGk
    end if
    !Determine all translationengleiche subgroups of the new space group
    call get_T_SubGroups(SpaceGroup_n,SubGroup,nsg,pointr)
    call get_C_Subgroups(SpaceGroup_n,SubG_C,nsgc,pointc)
    write(unit=lun,fmt="(/,/,a,a,a,/)") " => LIST of Translationengleiche Subgroups of index 2 for: ",&
                                             SpaceGroup_n%Spg_Symb,SpaceGroup_n%hall
    do i=1,nsg
      indx=SpaceGroup_n%Multip/SubGroup(i)%multip
      if(indx /= indice) cycle
      ng=SubGroup(i)%numops
      write(unit=lun,fmt="(4a,i2,30a)") " => ", SubGroup(i)%Spg_Symb, SubGroup(i)%hall,&
        " Index: [",indx,"]   ->  { ", ( trim(SubGroup(i)%SymopSymb(l))//" ; ",l=1,ng-1),&
        trim(SubGroup(i)%SymopSymb(ng))," }    ", trim(SubGroup(i)%centre)
    end do

    write(unit=lun,fmt="(/,/,a,a,a,/)") " => LIST of Classengleiche Subgroups of index 2 for: ",&
                                             SpaceGroup_n%Spg_Symb,SpaceGroup_n%hall
    do i=1,nsgc
      indx=SpaceGroup_n%Multip/SubG_C(i)%multip
      if(indx /= indice) cycle
      ng=SubG_C(i)%numops
      write(unit=lun,fmt="(4a,i2,30a)") " => ", SubG_C(i)%Spg_Symb, SubG_C(i)%hall,&
        " Index: [",indx,"]   ->  { ", ( trim(SubG_C(i)%SymopSymb(l))//" ; ",l=1,ng-1),&
        trim(SubG_C(i)%SymopSymb(ng))," }    ", trim(SubG_C(i)%centre)
    end do

    !Write the possible magnetic groups
    write(unit=lun,fmt="(3(/,a))") "---------------------------------",&
                                   "  Maximal Magnetic Space Groups",&
                                   "---------------------------------"
    do j=1,nsg
      indx=SpaceGroup_n%Multip/SubGroup(j)%multip
      if( indx /= indice) cycle
      write(unit=lun,fmt="(a)") "Magnetic Group derived from: "//trim(SpaceGroup_n%Spg_Symb)//" and subgroup: "//trim(SubGroup(j)%hall)
      do i=1,SpaceGroup_n%Multip
        tinv=",+1"
        if(.not. pointr(i,j)) tinv=",-1"
        write(unit=lun,fmt="(i4,a)") i,"  "//trim(SpaceGroup_n%SymOpSymb(i))//tinv
      end do
    end do

    do j=1,nsgc
      indx=SpaceGroup_n%Multip/SubG_C(j)%multip
      if( indx /= indice) cycle
      write(unit=lun,fmt="(a)") "Magnetic Group derived from: "//trim(SpaceGroup_n%Spg_Symb)//" and subgroup: "//trim(SubG_C(j)%hall)
      do i=1,SpaceGroup_n%Multip
        tinv=",+1"
        if(.not. pointc(i,j)) tinv=",-1"
        write(unit=lun,fmt="(i4,a)") i,"  "//trim(SpaceGroup_n%SymOpSymb(i))//tinv
      end do
    end do

    write(unit=*,fmt="(a,a)") " => Results in file: ",  trim(outfil)//".magr"
    stop
 End Program MagGroups_k
