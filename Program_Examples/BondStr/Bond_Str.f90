!!----
!!---- Program Bond_STR
!!----
!!---- Using CrysFML versions above 4.01
!!----
Program Bond_Str
   !---- Use Modules ----!
   !use f2kcli  !To be used with Lahey Fortran by removing the comment
   use CFML_GlobalDeps,                  only: Cp
   use CFML_String_Utilities,            only: u_case, pack_string, setnum_std
   use CFML_Math_General,                only: sort,Set_Epsg, Set_Epsg_default,Modulo_Lat, Equal_Vector
   use CFML_Math_3D,                     only: Determ_A, Cross_Product, Polyhedron_Volume, Get_Spheric_Coord, &
                                               Get_Cart_From_Spher, Get_Centroid_Coord, err_Math3D, err_Math3D_Mess
   use CFML_Crystal_Metrics,             only: Write_Crystal_Cell, Crystal_Cell_Type
   use CFML_Crystallographic_Symmetry,   only: ApplySo,Write_SpaceGroup,Space_Group_Type
   use CFML_Atom_TypeDef,                only: Write_Atom_List, Atom_list_Type, Allocate_Atom_List, Deallocate_Atom_List,AtList1_ExtenCell_AtList2
   use CFML_Geometry_Calc,               only: Calc_Dist_Angle_Sigma, Coord_Info, Err_Geom, Err_Geom_Mess, Distance, &
                                               Allocate_Coordination_Type, Set_TDist_Coordination
   use CFML_IO_Formats,                  only: Readn_set_Xtal_Structure,ERR_Form_Mess,err_form,file_list_type,Write_CFL
   Use CFML_Scattering_Chemical_Tables
   Use CFML_Export_VTK
   use CFML_BVS_Energy_Calc

   !---- Variables ----!
   Implicit None

   type (Space_Group_Type)         :: SpGr
   type (Crystal_Cell_Type)        :: Cell
   type (Atom_list_Type)           :: A
   type (Atoms_Conf_list_Type)     :: Ac
   type (File_List_Type)           :: Fich_cfl

   character(len=256)              :: filcod
   character(len=256)              :: line,title
   character(len=80),dimension(15) :: bvparm,bvelparm
   character(len=30)               :: ElectCon
   character(len=4)                :: atname,atm
   character(len=2)                :: Chem

   integer                         :: lun=1, ier,i, lr,ln, i_cfl=2, i_cons=3
   integer                         :: narg,iargc,maxc
   integer                         :: nc,n_bvpar,n_bvelpar    ,n1,n2,j
   integer                         :: ndimx,ndimy,ndimz

   logical                         :: esta, arggiven=.false.,sout=.false.,cif=.false.
   logical                         :: read_bvparm=.false., restr=.false., bvs_calc=.true.
   logical                         :: vdist=.false.,read_bvelparm=.false.
   logical                         :: map_calc=.false., soft=.false., bvel_calc=.false.
   real                            :: ttol=20.0,dmax,dangl, gii
   real                            :: drmax,delta,qval,tini,tfin,qn
    !!----
    !!---- Character(len=22), dimension(:), allocatable, public :: Electronic_Configuration
    !!----
    !!----  Electronic configurations of the elements. Alloctate to Num_Chem_Info = 108
    !!----  The indices correspond to the atomic number Z
    !!----
    !!----  Created: January 2015
    !!
    Character(len=30), dimension(:), allocatable :: Electronic_Configuration

   ! Arguments on the command line
   lr=0
   ln=0
   !narg=iargc()
   narg=Command_Argument_Count()

   if (narg > 0) then
      !call getarg(1,filcod)
      call GET_COMMAND_ARGUMENT(1,filcod)
      i=index(filcod,".cfl")
      if(i /= 0) filcod=filcod(1:i-1)
      i=index(filcod,".cif")
      if(i /= 0) filcod=filcod(1:i-1)
      arggiven=.true.
   end if

   write(unit=*,fmt="(/,/,6(a,/))")                                                        &
           "                      =============================="                        , &
           "                      ====== PROGRAM BOND_STR ======"                        , &
           "                      =============================="                        , &
           "    ***********************************************************************" , &
           "    * Distances, angles and Bond-Valence Sums from  *.cfl or *.cif files  *" , &
           "    ***********************************************************************" , &
           "                     (JRC - ILL, version: December 2014 )"
   write(unit=*,fmt=*) " "
   call set_epsg(0.001_cp)  !needed for well controlling the calculation of multiplicities
   if(.not. arggiven) then
      write(unit=*,fmt="(a)") " => Code of the file xx.cfl(cif) (give xx): "
      read(unit=*,fmt="(a)") filcod
      if(len_trim(filcod) == 0) stop
   end if

   inquire(file=trim(filcod)//".cfl",exist=esta)
   call cpu_time(tini)
   if(esta) then
      call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpGr,A,Mode="CFL",file_list=fich_cfl)
      cif=.false.
   else
      inquire(file=trim(filcod)//".cif",exist=esta)
      if(.not. esta) then
        write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cfl (or .cif) does'nt exist!"
        stop
      end if
      call Readn_set_Xtal_Structure(trim(filcod)//".cif",Cell,SpGr,A,Mode="CIF",file_list=fich_cfl)
      cif=.true.
   end if

   if (err_form) then
      write(unit=*,fmt="(a)") trim(ERR_Form_Mess)
      write(unit=*,fmt="(/,a)") " => PROGRAM BOND_STR finished in error!"
      stop
   end if

   open(unit=lun,file=trim(filcod)//".bvs", status="replace",action="write")
   write(unit=lun,fmt="(/,/,6(a,/))")                                                  &
        "                      =============================="                        , &
        "                      ====== PROGRAM BOND_STR ======"                        , &
        "                      =============================="                        , &
        "    ***********************************************************************" , &
        "    * Distances, angles and Bond-Valence Sums from  *.cfl or *.cif files  *" , &
        "    ***********************************************************************" , &
        "                    (JRC - ILL, version: February 2010 )"
   write(unit=lun,fmt="(a,/)") " "
   dmax=3.2
   dangl=0.0
   title=" "
   delta=0.0

   if (cif) then
      write(unit=lun,fmt="(a,/)") " => Data obtained from CIF file: "//trim(filcod)//".cif"

      open(unit=i_cfl,file="CFL_file.cfl",status="replace",action="write")
      write(unit=i_cfl,fmt="(a)") "Title  CFL-file generated from CIF file: "//trim(filcod)//".cif"
      call Write_CFL(i_cfl,Cell,SpGr,A)
      close(unit=i_cfl)

      write(unit=*,fmt="(a)") " => A CFL-file has been generated from CIF -> CFL_file.cfl"
      write(unit=*,fmt="(a)") "    This file may be used to add instructions for BOND_STR"

   else
      write(unit=lun,fmt="(a,/)") " => Content of the input file: "
      do i=1,fich_cfl%nlines
         write(unit=lun,fmt="(tr10,a)") fich_cfl%line(i)
      end do
   end if

   write(unit=lun,fmt="(a,/)") " "

   call Write_Crystal_Cell(Cell,lun)
   call Write_SpaceGroup(SpGr,lun,full=.true.)
   call Write_Atom_List(A,lun=lun)

   ! Check if Bond-Valence calculation are possible
   do i=1,A%natoms
     if(abs(A%atom(i)%charge) <= 0.001)  then
       bvs_calc=.false.
       exit
     end if
   end do

   ! Look for specific instructions for BOND_STR in CFL-file
   ttol=30.0
   if (cif) then
      sout=.true.
   else
      n_bvpar=0; n_bvelpar=0
      bvparm=" "; bvelparm=" "
      do i=1,fich_cfl%nlines
         line=adjustl(u_case(fich_cfl%line(i)))

         if (line(1:4) == "TITL") then
            title=line(6:)
            cycle
         end if

         if (line(1:5) == "BVPAR") then
            read_bvparm=.true.
            n_bvpar=n_bvpar+1
            nc=index(line,' ')
            bvparm(n_bvpar)=adjustl(line(nc+1:))
            cycle
         end if

         if (line(1:7) == "BVELPAR") then
            read_bvelparm=.true.
            n_bvelpar=n_bvelpar+1
            nc=index(line,' ')
            bvelparm(n_bvpar)=adjustl(line(nc+1:))
            cycle
         end if

         if (line(1:4) == "DMAX") then
            read(unit=line(5:),fmt=*,iostat=ier) dmax,dangl
            if (ier /= 0) then
               dmax=3.2
               dangl=0.0
            end if
            cycle
         end if

         if (line(1:7) == "SOFTBVS") then
            soft=.true.
            cycle
         end if

         if (line(1:3) == "TOL") then
            read(unit=line(4:),fmt=*,iostat=ier) ttol
            if (ier /= 0) then
               ttol=30.0
            end if
            if(ttol < 0.001) ttol=30.0
            cycle
         end if

         if (line(1:6) == "DISTAN")  sout=.true.

         if (line(1:5) == "RESTR")  restr=.true.

         !--- JGP ----!
         if (line(1:4) == "MAP " .or. line(1:4) == "BVEL") then

            read(unit=line(5:),fmt=*,iostat=ier) ndimx,ndimy,ndimz,atname,drmax,delta
            if (ier /= 0) then
               read(unit=line(5:),fmt=*,iostat=ier) ndimx,ndimy,ndimz,atname,drmax
               if (ier /= 0) then
                  ndimx=32
                  ndimy=32
                  ndimz=32
                  atname=" "
                  drmax=4.0
               end if
            end if

            atm=u_case(atname)
            call Get_Chemsymb(atm,chem)
            chem=u_case(chem)
            nc=index(atname,"+")

            if(nc > 0) then
              read(atname(nc+1:),*,iostat=ier) qval
            else
              nc=index(atname,"-")
              if(nc > 0) read(atname(nc:),*,iostat=ier) qval
            end if

            if(line(1:3) == "MAP") then
              map_calc=.true.
            else
              bvel_calc=.true.
            end if
            cycle

         end if

         if (line(1:5) == "VDIST") vdist=.true.

      end do
   end if !cif

   ! Distances, Angles, Restraints Calculations
   if (restr) then
      call Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spgr, A, lun, i_cons)
   else if(sout) then
      call Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spgr, A, lun)
   else
      call Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spgr, A)
   end if
   if (err_geom) then
      write(unit=lun,fmt="(/,a)") ERR_Geom_Mess
      write(unit=*,  fmt="(/,a)") ERR_Geom_Mess
      stop
   end if

   ! Pair distribution bonds
   call Calc_PDB(A,Dmax,lun)

   if(bvs_calc .or. vdist .or. map_calc .or. bvel_calc) then
     if(map_calc .or. bvel_calc ) then
        call Allocate_Atoms_Conf_List(A%natoms+1,Ac)
        Ac%atom(A%natoms+1)%Lab=atname
        Ac%atom(A%natoms+1)%ChemSymb=chem
        Ac%atom(A%natoms+1)%SfacSymb=atname
        Ac%atom(A%natoms+1)%charge=qval
     else
        call Allocate_Atoms_Conf_List(A%natoms,Ac)
     end if
     Ac%atom(1:A%natoms)=A%atom

     if(bvel_calc) then
        call Species_on_List(Ac,SpGr%Multip,ttol,.true.) !using covalent radius
     else
        call Species_on_List(Ac,SpGr%Multip,ttol)
     end if
     if (err_conf) then
        write(unit=*,fmt="(/,a)")   "    Species on list: "//trim(ERR_Conf_Mess)
        write(unit=*,fmt="(a,/)")   " => PROGRAM BOND_STR finished in error!"
        write(unit=lun,fmt="(/,a)") "    Species on list: "//trim(ERR_Conf_Mess)
        write(unit=lun,fmt="(a)")   " => PROGRAM BOND_STR finished in error!"
        stop
     end if
   end if

   ! Distortion calculations based on Centroids
   if (vdist) call Calc_Distortion_IVTON(Ac, Spgr, Cell,0,lun)

   ! BVS Calculations
   if (bvs_calc) then
      if(bvel_calc) then
         if(read_bvelparm) then
            call Set_Table_BVEL_Params(Ac,N_bvelpar,bvelparm)
         else
            call Set_Table_BVEL_Params(Ac)
         end if


         write(unit=lun,fmt="(a,/,a,/)")  &
         " Bond-Valence Energy parameters (D0,Rmin,alpha) for Morse Potential:  D0*{exp(alpha(dmin-d))-1}^2-1", &
               "   (data read from internal table or provided by the user)"
         do n1=1,Ac%N_Cations
            do j=1,Ac%N_Anions
               n2=Ac%N_Cations+j
               write(unit=lun,fmt="(2(a,i3,a,a4),/,3(a,f9.5),/,3(a,f9.5),a)")           &
                     "   Type",n1,": ",Ac%Species(n1)," with type",n2,": ",Ac%Species(n2), &
                     "    D0  =",Table_Dzero(n1,n2),"       Rmin =",Table_Rmin(n1,n2),   &
                     "  Alpha =",Table_Alpha(n1,n2),"  Av. Coord.=",Table_Avcoor(n1,n2), &
                     "    R0  =",Table_Rzero(n1,n2),"   R-cutoff =",Table_Rcutoff(n1,n2),&
                     "   => Reference: "//trim(references(Table_ref(n1,n2)))
               write(unit=lun,fmt="(2(a,a,a,f6.3,a),/)") &
                     "   Cation (Eff. radius): ",Ac%Species(n1),"(",Ac%Radius(n1),")   ",  &
                     "   Anion  (Eff. radius): ",Ac%Species(n2),"(",Ac%Radius(n2),")"
            end do
         end do
      else
         ! Setting Tables for B and D0
         if (read_bvparm) then
           if(soft) then
             call Set_Table_d0_b(Ac,n_bvpar,bvparm,soft)  !The table completes with extra data
           else
             call Set_Table_d0_b(Ac,n_bvpar,bvparm)  !The table completes with extra data
           end if
         else
           if(soft) then
             call Set_Table_d0_b(Ac,soft=soft)  !The table completes with extra data
           else
             call Set_Table_d0_b(Ac)
           end if
         end if
      end if


      if (err_conf) then
         write(unit=*,fmt="(/,a)")   "    Setting Table: "//trim(ERR_Conf_Mess)
         write(unit=*,fmt="(a,/)")   " => PROGRAM BOND_STR finished in error!"
         write(unit=lun,fmt="(/,a)") "    Setting Table: "//trim(ERR_Conf_Mess)
         write(unit=lun,fmt="(a)")   " => PROGRAM BOND_STR finished in error!"
         stop
      end if

      !write(unit=lun,fmt="(/,a)") " => Electronic Configuration of the different species"
      !do i=1,Ac%natoms
      !  Call Get_Electronic_Configuration(Ac%atom(i)%ChemSymb,Ac%atom(i)%charge,ElectCon,qn)
      !  write(unit=lun,fmt="(2(a,f8.3))") "   "//Ac%atom(i)%ChemSymb,Ac%atom(i)%charge," : "//trim(ElectCon),qn
      !end do
      !stop

      ! BVS
      ! Maps of BVS or BVEL
      if (map_calc) then

         if(delta > 0.001) then
           call Calc_Map_BVS(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,delta)
           !call Calc_Map_BVS_fast(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,delta)
         else
           call Calc_Map_BVS(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax)
         end if

      else if(bvel_calc) then

         call Calc_Map_BVEL(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,delta)

      else

         call Calc_BVS(Ac,lun,filecod=filcod)

      end if

   else
      write(unit=lun,fmt="(/,a)")" => Bond-Valence calculations imposible: ionic charges must be provided!"
   end if

   if (err_conf) then
      write(unit=lun,fmt="(/,a)") ERR_Conf_Mess
      write(unit=*,  fmt="(/,a)") ERR_Conf_Mess
   end if

   if (.not. cif .and. bvs_calc .and. .not. map_calc .and. .not. bvel_calc ) then
      call system("type "//trim(filcod)//"_sum.bvs")
   end if
   write(unit=*,fmt="(/,a)")   " => Normal End of: PROGRAM BOND_STR "
   write(unit=lun,fmt="(/,a)") " => Normal End of: PROGRAM BOND_STR "
   write(unit=*,fmt="(a)")     " => Global results in File: "//trim(filcod)//".bvs"
   if(bvs_calc) then
     if (map_calc) then
       write(unit=*,fmt="(a)")     " => Bond Valence Map in File: "//trim(filcod)//".map"
       write(unit=*,fmt="(a)")     " => Binary VESTA File: "//trim(filcod)//"_bvs.pgrid"
     else if (bvel_calc) then
       write(unit=*,fmt="(a)")     " => Bond Valence Energy Landscape in File: "//trim(filcod)//"_bvel.map"
       write(unit=*,fmt="(a)")     " => Binary BVEL VESTA File: "//trim(filcod)//"_bvel.pgrid"
     else
       write(unit=*,fmt="(a)")     " => Summary of BVS in File: "//trim(filcod)//"_sum.bvs"
     end if
   end if

   call cpu_time(tfin)
   write(unit=*,fmt="(a,f8.4,a)")     " => CPU-time: ",tfin-tini," seconds"
   write(unit=lun,fmt="(a,f8.4,a)")   " => CPU-time: ",tfin-tini," seconds"
   close(unit=lun)
   call Deallocate_Atoms_Conf_List(Ac)

 Contains
   !!----
   !!----
   !!----
   !!
   Subroutine Calc_Distortion_IVTON(At, Spg, Cell, Model_Vi,Iunit)
      !---- Arguments ----!
      type (Atoms_Conf_list_Type), intent(in)  :: At
      type (Space_Group_Type),     intent(in)  :: Spg
      type (Crystal_Cell_Type),    intent(in)  :: Cell
      integer,                     intent(in)  :: Model_Vi
      integer, optional,           intent(in)  :: Iunit

      !---- Local Variables ----!
      real(kind=cp), dimension(3)     :: Atc,tr, cent, bari, cent_fr, bari_fr,te
      real(kind=cp), dimension(3,12)  :: Atv
      real(kind=cp)                   :: Vp, rm, srm, spher, delta, eccent, aspher, vecc, &
                                         vaspher, vspher, Vi, Vdis, dd, q1,q2
      integer                         :: lun,i,j,k,n,cn,ecn,l1,l2,sig1,sig2
      character(len=40)               :: car

      ! Init control
      if (At%natoms <=0) return

      lun=6
      if (present(iunit)) lun=iunit

      write(unit=lun,fmt="(/,a)")    "  -----------------------------------------------"
      write(unit=lun,fmt="(a)")      "  ---- Distortion for Coordination Polyhedra ----"
      write(unit=lun,fmt="(a)")      "  -----------------------------------------------"
      write(unit=lun,fmt="(a)") " "


      do i=1,At%natoms
         ! Coordination number
         cn=coord_info%coord_num(i)
         if (cn < 3) cycle
         l1=At%Atom(i)%ind(1)      !--------
         q1=At%Atom(i)%charge
         sig1=SIGN(1.0_cp,q1)
         ! Convert to cartesian coordinates
         atc=matmul(cell%cr_orth_cel,At%atom(i)%x)
         ! Printing Info
         write(unit=lun, fmt='(a,5x,3f10.4,5x,3f10.4)') ' Central Atom (Fract. & Cartesian Coord.): '//trim(At%atom(i)%Lab), &
                                                          At%atom(i)%x,atc
         write(unit=lun, fmt='(a, i3)') ' Atoms around for which distances are calculated: ', cn
         write(unit=lun, fmt='(a)') &
              ' Label         Fract. Coord.          Symm. Operator             Traslation                 '// &
              '    Ext. Coord.               Cartesian Coord.           Distance'

         ecn=0
         do j=1,cn
            q2=At%Atom(coord_info%n_cooatm(j,i))%charge
            sig2=SIGN(1.0_cp,q2)
            if(sig1 == sig2) cycle
            l2=At%Atom(coord_info%n_cooatm(j,i))%ind(1)
            dd=coord_info%dist(j,i)
            if (dd > (At%Radius(l1)+At%Radius(l2))*(1.0+0.01*At%tol)) cycle
            ecn=ecn+1
            k=coord_info%n_cooatm(j,i)
            n=coord_info%n_sym(j,i)
            tr=coord_info%tr_coo(:,j,i)
            te=ApplySO(Spg%SymOp(n),At%atom(k)%x)+tr
            Atv(:,ecn)=te
            Atv(:,ecn)=matmul(cell%cr_orth_cel,Atv(:,ecn))
            write(unit=lun, fmt='(3x,a,t10,3f8.4,t38,a,t57,3f9.4,t90,3f8.4,5x,3f8.4,5x,f10.4)') &
                  At%atom(k)%lab, At%atom(k)%x, trim(Spg%SymOpSymb(n)), tr, te,Atv(:,ecn), coord_info%dist(j,i)
         end do
         write(unit=lun, fmt='(a, i3)') ' Coordination number: ', ecn

         ! Polehedra volume (Vp)
         vp= Polyhedron_Volume(ecn,Atv,Atc)
         write(unit=lun, fmt='(a, f8.3)') ' Polyhedron volume: ', vp

         ! Centroid
         call Get_Centroid_Coord(ecn,Atv,cent,bari)
         if(err_Math3D) then
          write(unit=*,fmt="(t10,a)") "=>"//trim(err_Math3D_mess)
          write(unit=lun,fmt="(t10,a)") "=>"//trim(err_Math3D_mess)
          cycle
         end if
         cent_fr=matmul(cell%Orth_Cr_cel,cent)
         bari_fr=matmul(cell%Orth_Cr_cel,bari)

         ! Average distance (Rs)
         rm=0.0
         do j=1,ecn
            rm=rm+distance(cent,Atv(1:3,j))
         end do
         rm=rm/real(ecn)

         ! Standard deviation of Rs
         srm=0.0
         do j=1,ecn
            srm=srm+(distance(cent,Atv(1:3,j)) - rm)**2
         end do
         srm=sqrt(srm/real(ecn-1))

         ! Displacement Central Atom (Delta)
         delta=distance(cent,Atc)

         ! Other geometric parameters
         eccent=delta/rm                   ! Eccentrincity
         aspher=srm/rm                     ! Asphericity
         spher=1.0-aspher                  ! Sphericity

         vecc=1.0 - ((rm-delta)/rm)**3     ! Volume eccentrincity
         vaspher=3.0*aspher                ! Volume asphericity
         vspher=1.0 - vaspher              ! Volume sphericity

         ! Ideal volume (Vi): Maximum-volume polyhedra
         if (model_vi == 0) then
            select case (ecn)
               case (4) ! Tetrahedron
                  vi=(8.0/(9.0*sqrt(3.0)))*(rm**3)

               case (5) ! trigonal bipyramid
                  vi=(sqrt(3.0)/2.0)*(rm**3)

               case (6) ! Octahedron
                  vi=(4.0/3.0)*(rm**3)

               case (7) ! Pentagonal bipyramid
                  vi=1.58509*(rm**3)

               case (8) ! Bisdisphenoid
                  vi=1.81571*(rm**3)

               case (9) ! Tricapped trigonal prism
                  vi=2.0437*(rm**3)

               case (10)
                  vi=2.24*(rm**3)

               case (11)
                  vi=2.4*(rm**3)

               case (12) ! Icosahedron
                  vi=2.5362*(rm**3)

               case default
                  vi=vp
            end select
         else
            ! Special cases for Vi calculations

         end if

         ! Volume discrepancy
         Vdis=100.0*((vi - vp)/vi)

         write(unit=lun, fmt='(a,3f10.4)') ' Centroid   coordinates: ',cent_fr
         write(unit=lun, fmt='(a,3f10.4)') ' Baricenter coordinates: ',bari_fr
         call setnum_std(rm,srm,car)
         write(unit=lun, fmt='(a)')       ' Average distance from Centroid: '//trim(car)
         write(unit=lun, fmt='(a,f10.3)') ' Distance of the central atom to Centroid: ',delta
         write(unit=lun, fmt='(a,f10.3)') ' Linear eccentricity: ', eccent
         write(unit=lun, fmt='(a,f10.3)') ' Linear  asphericity: ', aspher
         write(unit=lun, fmt='(a,f10.3)') ' Linear   sphericity: ', spher

         write(unit=lun, fmt='(a,f10.3)') ' Volume of  coordination Polyhedron: ',vp
         write(unit=lun, fmt='(a,f10.3)') ' Volume of the idealized Polyhedron: ',vi
         write(unit=lun, fmt='(a,f10.3)') ' Volume   discrepancy factor (v%)  : ',vdis

         write(unit=lun, fmt='(a,f10.3)') ' Volume eccentricity: ',vecc
         write(unit=lun, fmt='(a,f10.3)') ' Volume  asphericity: ',vaspher
         write(unit=lun, fmt='(a,f10.3)') ' Volume   sphericity: ',vspher
         write(unit=lun, fmt='(a)') ' '
      end do

      return
   End Subroutine Calc_Distortion_IVTON

   !!----
   !!----
   !!----
   Subroutine Calc_PDB(A,Dmax,lun)
      !---- Arguments ----!
      type (Atom_list_Type), intent(in)         :: A
      real,                  intent(in)         :: Dmax
      integer,               intent(in)         :: lun

      !---- Local variables ----!
      integer, parameter                        :: Max_NSpecies=20
      character(len=2), dimension(Max_NSpecies) :: species
      character(len=2)                          :: car

      integer                                   :: n_spec
      integer, dimension(:,:),allocatable       :: bond_spec
      integer                                   :: i,j,n,n_c,n1,n2

      integer, parameter                        :: Max_dis=8000
      integer                                   :: n_dis,n_ini
      integer, dimension(Max_Dis)               :: indx
      real,dimension(Max_Dis)                   :: disbond
      real                                      :: dis_ini

      type (Atom_list_Type)                     :: Ac

      !---- Init control ----!
      if (A%natoms <=0) return

      write(unit=lun,fmt="(/,a)")          "  ---------------------------------------------"
      write(unit=lun,fmt="(a,f6.3,a)")     "  {--- BONDS DISTRIBUTIONS (UP TO ",dmax, " ) ---}"
      write(unit=lun,fmt="(a)")            "  ---------------------------------------------"

      !---- Make a copy of A ----!
      call Allocate_Atom_List(A%natoms,Ac)
      Ac%Atom=A%atom

      !---- Calculate the number of different species in the List ----!
      n_spec=0
      species=" "
      at1:do n=1,Ac%natoms
         if (n_spec > 0) then
            car=u_case(Ac%Atom(n)%ChemSymb)
            do i=1,n_spec
               if (car == species(i)) cycle at1
            end do
         end if
         if (n_spec == Max_NSpecies) then
            write(unit=lun,fmt="(/,a)")  "Overflow the number of different species the program can use"
            call Deallocate_atom_list(Ac)
            return
         end if
         n_spec=n_spec+1
         species(n_spec)=u_case(Ac%Atom(n)%ChemSymb)
      end do at1

      !---- Assign the index in the Atom list respect to species ----!
      do n=1,Ac%natoms
         do i=1,n_spec
            if (u_case(Ac%Atom(n)%ChemSymb) == species(i)) then
               Ac%atom(n)%ind(1)=i
               exit
            end if
         end do
      end do

      !---- Creating Table for bonds between species ----!
      allocate(bond_spec(n_spec,n_spec))
      bond_spec=0

      do n=1,coord_info%natoms
         n1=Ac%atom(n)%ind(1)
         n_c=coord_info%coord_num(n)
         do i=1,n_c
            j=coord_info%n_cooatm(i,n)
            n2=Ac%atom(j)%ind(1)
            bond_spec(n1,n2)=bond_spec(n1,n2)+1
            bond_spec(n2,n1)=bond_spec(n2,n1)+1
         end do
      end do

      !---- Calculate the Distributions Information ----!
      do n1=1,n_spec
         do n2=1,n_spec
            if (bond_spec(n1,n2) == 0) cycle
           !if (n1 > n2) cycle !This is an error!!!!!

            n_dis=0
            indx=0
            disbond=0.0
            do n=1,coord_info%natoms
               if (Ac%atom(n)%ind(1) /= n1) cycle
               n_c=coord_info%coord_num(n)
               do i=1,n_c
                  j=coord_info%n_cooatm(i,n)
                  if (Ac%atom(j)%ind(1) /= n2) cycle
                  if (n_dis+1 > Max_Dis) then
                     write(unit=lun,fmt="(a,i6/)")  "Overflow the number of distances the program can handle: ",max_dis
                     deallocate(bond_spec)
                     call Deallocate_atom_list(Ac)
                     return
                  end if
                  n_dis=n_dis+1
                  disbond(n_dis)=coord_info%dist(i,n)
               end do
            end do

            call sort(disbond,n_dis,indx)

            write(unit=lun,fmt='(/,a,i4)') " Bond type: "//species(n1)//" - "//species(n2)//"      Num: ",n_dis
            write(unit=lun,fmt='(a)')   "  Num. Bonds         Distance"

            n_ini=1
            dis_ini=disbond(indx(n_ini))
            do i=2,n_dis
               if (abs(disbond(indx(i))-dis_ini) <= 0.001) then
                  if (i < n_dis) cycle
                  write(unit=lun,fmt='(3x,i3,12x,f12.4)') i-n_ini+1,disbond(indx(n_ini))
               else
                  write(unit=lun,fmt='(3x,i3,12x,f12.4)') i-n_ini,disbond(indx(n_ini))
                  n_ini=i
                  dis_ini=disbond(indx(i))
                  if (n_ini == n_dis) write(unit=lun,fmt='(3x,i3,12x,f12.4)') 1,dis_ini
               end if
            end do

         end do
      end do

      deallocate(bond_spec)
      call Deallocate_atom_list(Ac)

      return
   End Subroutine Calc_PDB

    Subroutine Calc_Map_BVSn(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta)
       !---- Arguments ----!
       type (Atoms_Conf_List_type), intent(in) :: A
       type (Space_Group_Type),     intent(in) :: SpG
       Type (Crystal_Cell_Type),    intent(in) :: Cell
       character(len=*),            intent(in) :: Filecod
       integer,                     intent(in) :: ndimx
       integer,                     intent(in) :: ndimy
       integer,                     intent(in) :: ndimz
       character(len=*),            intent(in) :: atname  !Species of the probe atom e.g. LI+1
       real(kind=cp),               intent(in) :: drmax
       real(kind=cp), optional,     intent(in) :: delta

       !---- Local variables ----!
       character(len=4)                             :: car,atm,chem
       integer                                      :: i,j,k,n,n1,n2,np,L,nL
       integer                                      :: nx1,nx2,ny1,ny2,nz1,nz2,nn
       integer                                      :: i1,j1,k1
       integer                                      :: jbvs
       integer,      dimension(:),     allocatable  :: ind
       real(kind=cp)                                :: rx1,ry1,rz1,qval,dif,q1,q2,sig1,sig2,rep
       real(kind=cp)                                :: sbvs, dd, occ, radius,sig,dmin
       real(kind=cp), dimension(3)                  :: pto,pta,step,extend
       real(kind=cp), dimension(:,:,:), allocatable :: map_bvs
       real(kind=cp), dimension(:,:),   allocatable :: peaks
       real(kind=cp), dimension(:),     allocatable :: VD_peaks
       type (Atom_list_Type)                        :: At1, At2
       logical                                      :: new_peak,anion,cation

       !---- Initial conditions ----!
       if (A%natoms <= 0) return
       if (ndimx <= 0 .or. ndimy <= 0 .or. ndimz <= 0) return

       !---- Preparing the Variables and calculations ----!
       call Allocate_Atom_List(A%natoms,At1)
       At1%atom=A%atom  !Atom list A%Atom(:)%ind(1) contains the pointer to the species
                        !coming from A (Atoms_Conf_List_type, set in the main program)
       atm=u_case(atname)
       n1=0
       do i=1,At1%natoms  !This exclude all atoms with the same species as Atname
          n2=A%Atom(i)%ind(1)
          if (trim(u_case(A%species(n2))) == trim(atm)) then
             At1%atom(i)%active=.false.
             q1=At1%atom(i)%charge
             sig1=SIGN(1.0_cp,q1)
             n1=n2
             radius=A%radius(n1)
             car=A%Species(n1)
          end if
       end do
       if (n1 ==0) then
          err_conf=.true.
          ERR_Conf_Mess="The point atom "//atname//" is not in the Species Table"
          return
       end if

       call AtList1_ExtenCell_AtList2(Spg,At1,At2,.true.)
       call Deallocate_atom_list(At1)
       !check that all species are well set in the list
       do n=1,At2%natoms
           n2=At2%Atom(n)%ind(1)
           if (n2 ==0) then
              err_conf=.true.
              ERR_Conf_Mess="The atom "//trim(At2%Atom(n)%lab)//" is not in the Species Table"
              return
           end if
       end do

       allocate(map_bvs(ndimx,ndimy,ndimz))
       map_bvs=0.0

       if(.not. present(delta)) then
          allocate(peaks(4,ndimx)) ! A maximum of ndimx peaks will be stored
          peaks=0.0
          allocate(VD_peaks(ndimx)) ! Vector holding the differences of bond-valences
          VD_peaks=0.0
          allocate(ind(ndimx))      ! Vector pointer for ordering the peaks
          ind=0
       end if

       ! Determination of the valence expected for a good position
       i=index(car,"+")
       cation=.false.
       anion=.false.
       if( i /= 0) then
         read(unit=car(i+1:),fmt=*) qval
         cation=.true.
       else
         i=index(car,"-")
         read(unit=car(i+1:),fmt=*) qval
         anion=.true.
       end if
       step=(/ 1.0/real(ndimx),  1.0/real(ndimy), 1.0/real(ndimz) /)
       extend=(/ drmax/cell%cell(1),  drmax/cell%cell(2), drmax/cell%cell(3) /)

       np=0
       !---- Map Calculation ----!
       do k=1,ndimz
          pto(3)=(k-1)*step(3)
          do j=1,ndimy
             pto(2)=(j-1)*step(2)
             do i=1,ndimx
                pto(1)=(i-1)*step(1)

                rx1=pto(1)-extend(1)
                if (rx1 <= 0.0) then
                    nx1=int(rx1)-1
                else
                    nx1=int(rx1)
                end if
                nx2=int(pto(1)+extend(1))

                ry1=pto(2)-extend(2)
                if (ry1 <= 0.0) then
                    ny1=int(ry1)-1
                else
                    ny1=int(ry1)
                end if
                ny2=int(pto(2)+extend(2))

                rz1=pto(3)-extend(3)
                if (rz1 <= 0.0) then
                    nz1=int(rz1)-1
                else
                    nz1=int(rz1)
                end if
                nz2=int(pto(3)+extend(3))

                sbvs=0.0
                rep=0.0
                do n=1,At2%natoms
                   n2=At2%Atom(n)%ind(1)
                   !if((cation .and. n2 <= A%N_cations)  .or. (anion .and. n2 > A%N_cations) ) cycle
                   q2=At2%Atom(n)%charge
                   sig2= SIGN(1.0_cp,q2)
                   sig=(radius+A%radius(n2))*0.99
                   do k1=nz1,nz2
                      do j1=ny1,ny2
                         do i1=nx1,nx2
                            pta=At2%Atom(n)%x+real((/i1,j1,k1/))
                            occ=At2%Atom(n)%VarF(1)
                            dd=Distance(pto,pta,Cell)
                            if (dd > drmax) cycle
                            if (sig1 == sig2) then
                               rep=rep + (sig/dd)**18
                               cycle
                            end if
                            sbvs=sbvs+occ*exp((table_d0(n1,n2)-dd)/table_b(n1,n2))
                         end do
                      end do
                   end do
                end do
                dif=abs(sbvs-qval)
                if(present(delta)) then
                  if(dif > delta .or. rep > 0.01) sbvs=0.0_cp
                else
                  !Algorithm for selecting the optimum positions
                  if(dif < 0.15 .and. np < ndimx) then
                    new_peak=.true.
                    do L=1,np
                       dd=Distance(peaks(1:3,L),pto,Cell)
                       if( dd < 0.8 ) then
                         new_peak=.false.
                         nL=L
                         exit
                       end if
                    end do
                    if(new_peak) then
                      np=np+1
                      peaks(1:3,np)= pto
                      peaks(4,np)  = sbvs
                      VD_peaks(np) = dif
                    else  !now compare with the peak stored at nL and interchange them if sbvs is more favourable
                      if( dif < abs(qval-peaks(4,nL)) ) then
                        peaks(1:3,nL) = pto
                        peaks(  4,nL) = sbvs
                        VD_peaks(nL)  = dif
                      end if
                    end if
                  end if
                end if
                !end of peaks construction
                map_bvs(i,j,k)=sbvs
             end do
          end do
       end do
       !Sorting the favourable peak positions for inserting the species of atom Atname

       if(.not. present(delta)) call sort(VD_peaks,np,ind)
       !---- Export a File ----!
       !call Get_LogUnit(jbvs)
       open(newunit=jbvs,file=trim(filecod)//".map",status="replace",action="write")

       write (unit=jbvs, fmt='(a)') "BVS Map Calculations using Bond_STR Program"
       write (unit=jbvs, fmt='(a)') "BVS Map for species "//trim(car)
       write (unit=jbvs, fmt='(a,3f12.4,3x,3f8.3)') "CELL ",cell%cell,cell%ang
       write (unit=jbvs, fmt='(a)')     "SPGR  "//trim(SpG%spg_symb)
       write (unit=jbvs, fmt='(a)')     "! List ot atoms  "
       do i=1,A%natoms
         write(unit=jbvs, fmt='(a,t20,5f12.5)')"Atom "//trim(A%Atom(i)%lab)//"  "//A%Atom(i)%SfacSymb, &
                                                A%Atom(i)%x, A%Atom(i)%biso, A%Atom(i)%occ
       end do
       if(.not. present(delta)) then
         write (unit=jbvs, fmt='(a)')     "! List ot favourable positions for inserting the species "//trim(car)
         do i=1,np
           j=ind(i)
           write(unit=jbvs, fmt='(a,i4,a,3f10.5,a,f8.3)')"#",i,"  Position: (",peaks(1:3,j),")  Valence: ",peaks(4,j)
         end do
       end if
       write (unit=jbvs, fmt='(a)')     "! Grid values: ndimx,ndimy,ndimz [0,ndimx-1] [0,ndimy-1] [0,ndimz-1]  "
       write (unit=jbvs, fmt='(a,9i6)') "GRID ",ndimx,ndimy,ndimz,0,ndimx-1,0,ndimy-1,0,ndimz-1
       write (unit=jbvs, fmt='(a)')     "DENSITY_MAP"
       write (unit=jbvs,fmt='(8g12.5)') map_bvs
       close(unit=jbvs)
       call write_grid_VESTA(map_bvs,cell,"Bond Valence Map",trim(filecod)//"_bvs","P")

       !---- End Procedure ----!
       if (allocated(map_bvs)) deallocate(map_bvs)
       call deallocate_atom_list(At2)

       return
    End Subroutine Calc_Map_BVSn


    Subroutine Calc_Map_BVS_fast(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta)
       !---- Arguments ----!
       type (Atoms_Conf_List_type), intent(in) :: A
       type (Space_Group_Type),     intent(in) :: SpG
       Type (Crystal_Cell_Type),    intent(in) :: Cell
       character(len=*),            intent(in) :: Filecod
       integer,                     intent(in) :: ndimx
       integer,                     intent(in) :: ndimy
       integer,                     intent(in) :: ndimz
       character(len=*),            intent(in) :: atname  !Species of the probe atom e.g. LI+1
       real(kind=cp),               intent(in) :: drmax
       real(kind=cp),               intent(in) :: delta

       !---- Local variables ----!
       character(len=4)                             :: car,atm,chem
       integer                                      :: i,j,k,n,n1,n2,np,L,nL
       integer                                      :: nx1,nx2,ny1,ny2,nz1,nz2,nn
       integer                                      :: i1,j1,k1
       integer                                      :: jbvs
       integer,      dimension(:),     allocatable  :: ind
       real(kind=cp)                                :: rx1,ry1,rz1,qval,dif,q1,q2,sig1,sig2,rep
       real(kind=cp)                                :: sbvs, dd, occ, radius,sig
       real(kind=cp), dimension(3)                  :: pto,pta,step,extend,qd
       real(kind=cp), dimension(:,:,:), allocatable :: map_bvs
       real(kind=cp), dimension(:,:),   allocatable :: peaks
       real(kind=cp), dimension(:),     allocatable :: VD_peaks
       type (Atom_list_Type)                        :: At1, At2
       logical                                      :: new_peak,anion,cation

       !---- Initial conditions ----!
       if (A%natoms <= 0) return
       if (ndimx <= 0 .or. ndimy <= 0 .or. ndimz <= 0) return

       !---- Preparing the Variables and calculations ----!
       call Allocate_Atom_List(A%natoms,At1)
       At1%atom=A%atom  !Atom list A%Atom(:)%ind(1) contains the pointer to the species
                        !coming from A (Atoms_Conf_List_type, set in the main program)
       atm=u_case(atname)
       do i=1,At1%natoms  !This exclude all atoms with the same species as Atname
          n2=A%Atom(i)%ind(1)
          if (trim(u_case(A%species(n2))) == trim(atm)) then
             At1%atom(i)%active=.false.
             q1=At1%atom(i)%charge
             sig1=SIGN(1.0_cp,q1)
          end if
       end do
       qd(:)=1.0/cell%rcell(:)

       extend=(/ drmax/cell%cell(1),  drmax/cell%cell(2), drmax/cell%cell(3) /)

       call AList2AList_inbox(Spg,At1,extend,At2)
       call Deallocate_atom_list(At1)
       !check that all species are well set in the list
       do n=1,At2%natoms
           n2=At2%Atom(n)%ind(1)
           !write(*,"(a,i4,a)") " => Atom #",n,"  "//trim(At2%Atom(n)%lab)//"  Species: "//A%Species(n2)
           if (n2 ==0) then
              err_conf=.true.
              ERR_Conf_Mess="The atom "//trim(At2%Atom(n)%lab)//" is not in the Species Table"
              return
           end if
       end do



       allocate(map_bvs(ndimx,ndimy,ndimz))
       map_bvs=0.0

       n1=0
       call Get_Chemsymb(atm,chem)
       chem=u_case(chem)
       do i=1,A%n_spec
          car=A%Species(i)
          if (index(car,trim(chem)) /= 0) then
             n1=i
             radius=A%radius(i)
             exit
          end if
       end do
       if (n1 ==0) then
          err_conf=.true.
          ERR_Conf_Mess="The point atom "//atname//" is not in the Species Table"
          return
       end if
       ! Determination of the valence expected for a good position
       i=index(car,"+")
       cation=.false.
       anion=.false.
       if( i /= 0) then
         read(unit=car(i+1:),fmt=*) qval
         cation=.true.
       else
         i=index(car,"-")
         read(unit=car(i+1:),fmt=*) qval
         anion=.true.
       end if
       step=(/ 1.0/real(ndimx),  1.0/real(ndimy), 1.0/real(ndimz) /)

       np=0
       !---- Map Calculation ----!
       do k=1,ndimz
          pto(3)=(k-1)*step(3)
          do j=1,ndimy
             pto(2)=(j-1)*step(2)
             do i=1,ndimx
                pto(1)=(i-1)*step(1)
                sbvs=0.0
                rep=0.0
                do_n: do n=1,At2%natoms
                   pta=At2%Atom(n)%x
                   do L=1,3
                      if (abs(pto(l)-pta(l))*qd(l) > drmax) cycle do_n
                   end do
                   n2=At2%Atom(n)%ind(1)
                   q2=At2%Atom(n)%charge
                   sig2= SIGN(1.0_cp,q2)
                   sig=(radius+A%radius(n2))*0.99
                   occ=At2%Atom(n)%VarF(1)
                   dd=Distance(pto,pta,Cell)
                   if (dd > drmax) cycle do_n
                   if (sig1 == sig2) then
                      rep=rep + (sig/dd)**18
                      cycle do_n
                   end if
                   sbvs=sbvs+occ*exp((table_d0(n1,n2)-dd)/table_b(n1,n2))
                end do do_n
                !write(*,"(a,3f8.4,3(a,f14.4))") " => Position: ",pto,"  BVS: ",sbvs,"  Rep: ",rep,"  dmin:",dmin,"  atom: "//At2%atom(nn)%lab
                dif=abs(sbvs-qval)
                if(dif > delta .or. rep > 0.01) sbvs=0.0_cp
                map_bvs(i,j,k)=sbvs
             end do
          end do
       end do
       !---- Export a File ----!
       open(newunit=jbvs,file=trim(filecod)//".map",status="replace",action="write")
       write (unit=jbvs, fmt='(a)') "BVS Map Calculations using Bond_STR Program"
       write (unit=jbvs, fmt='(a)') "BVS Map for species "//trim(car)
       write (unit=jbvs, fmt='(a,3f12.4,3x,3f8.3)') "CELL ",cell%cell,cell%ang
       write (unit=jbvs, fmt='(a)')     "SPGR  "//trim(SpG%spg_symb)
       write (unit=jbvs, fmt='(a)')     "! List ot atoms  "
       do i=1,A%natoms
         write(unit=jbvs, fmt='(a,t20,5f12.5)')"Atom "//trim(A%Atom(i)%lab)//"  "//A%Atom(i)%SfacSymb, &
                                                A%Atom(i)%x, A%Atom(i)%biso, A%Atom(i)%occ
       end do
       write (unit=jbvs, fmt='(a)')     "! Grid values: ndimx,ndimy,ndimz [0,ndimx-1] [0,ndimy-1] [0,ndimz-1]  "
       write (unit=jbvs, fmt='(a,9i6)') "GRID ",ndimx,ndimy,ndimz,0,ndimx-1,0,ndimy-1,0,ndimz-1
       write (unit=jbvs, fmt='(a)')     "DENSITY_MAP"
       write (unit=jbvs,fmt='(8g12.5)') map_bvs
       close(unit=jbvs)
       call write_grid_VESTA(map_bvs,cell,"Bond Valence Map",trim(filecod)//"_bvs","P")

       !---- End Procedure ----!
       if (allocated(map_bvs)) deallocate(map_bvs)
       call deallocate_atom_list(At2)

       return
    End Subroutine Calc_Map_BVS_fast

    Pure Function BVEL(d,d0,rmin,alpha,occ) result(BVEL_value)
      real, intent(in) :: d,d0,rmin,alpha,occ
      real             :: BVEL_value

      BVEL_value=occ*d0*(exp(alpha*(rmin-d))-1.0)**2-1.0

    End Function BVEL

    Subroutine Calc_Map_BVEL(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta)
       !---- Arguments ----!
       type (Atoms_Conf_List_type), intent(in) :: A
       type (Space_Group_Type),     intent(in) :: SpG
       Type (Crystal_Cell_Type),    intent(in) :: Cell
       character(len=*),            intent(in) :: Filecod
       integer,                     intent(in) :: ndimx
       integer,                     intent(in) :: ndimy
       integer,                     intent(in) :: ndimz
       character(len=*),            intent(in) :: atname  !Species of the probe atom e.g. LI+1
       real(kind=cp),               intent(in) :: drmax
       real(kind=cp),               intent(in) :: delta

       !---- Local variables ----!
       character(len=4)                             :: car,atm,chem
       integer                                      :: i,j,k,n,n1,n2,np,L,nL
       integer                                      :: nx1,nx2,ny1,ny2,nz1,nz2,nn
       integer                                      :: i1,j1,k1,sig1,sig2
       integer                                      :: jbvs
       integer,      dimension(:),     allocatable  :: ind
       real(kind=cp)                                :: rx1,ry1,rz1,qval,dif,q1,q2,rep
       real(kind=cp)                                :: sbvs, dd, occ, radius,sig,dmin, rho, coor,&
                                                       rzero,rcutoff,dzero, alpha,c_rep,c_atr
       real(kind=cp), dimension(3)                  :: pto,pta,step,extend
       real(kind=cp), dimension(:,:,:), allocatable :: map_bvs
       real(kind=cp), dimension(:,:),   allocatable :: peaks
       real(kind=cp), dimension(:),     allocatable :: VD_peaks
       type (Atom_list_Type)                        :: At1, At2
       logical                                      :: new_peak,anion,cation
       !Coulomb constant (1/4pie0) to get energies in eV, use electron charges and distances in angstroms
       real(kind=cp), parameter :: ke=14.399644850445155254866066


       !---- Initial conditions ----!
       if (A%natoms <= 0) return
       if (ndimx <= 0 .or. ndimy <= 0 .or. ndimz <= 0) return

       !---- Preparing the Variables and calculations ----!
       call Allocate_Atom_List(A%natoms,At1)
       At1%atom=A%atom  !Atom list A%Atom(:)%ind(1) contains the pointer to the species
                        !coming from A (Atoms_Conf_List_type, set in the main program)
       atm=u_case(atname)
       n1=0
       do i=1,At1%natoms  !This exclude all atoms with the same species as Atname
          n2=A%Atom(i)%ind(1)
          if (trim(u_case(A%species(n2))) == trim(atm)) then
             At1%atom(i)%active=.false.
             q1=At1%atom(i)%charge
             sig1=nint(SIGN(1.0_cp,q1))
             n1=n2
             radius=A%radius(n1)
             car=A%Species(n1)
          end if
       end do
       if (n1 ==0) then
          err_conf=.true.
          ERR_Conf_Mess="The point atom "//atname//" is not in the Species Table"
          return
       end if

       call AtList1_ExtenCell_AtList2(Spg,At1,At2,.true.)
       call Deallocate_atom_list(At1)
       !check that all species are well set in the list
       do n=1,At2%natoms
           n2=At2%Atom(n)%ind(1)
           if (n2 ==0) then
              err_conf=.true.
              ERR_Conf_Mess="The atom "//trim(At2%Atom(n)%lab)//" is not in the Species Table"
              return
           end if
       end do

       allocate(map_bvs(ndimx,ndimy,ndimz))
       map_bvs=0.0
       ! Determination of the valence expected for a good position
       i=index(car,"+")
       cation=.false.
       anion=.false.
       if( i /= 0) then
         read(unit=car(i+1:),fmt=*) qval
         cation=.true.
       else
         i=index(car,"-")
         read(unit=car(i+1:),fmt=*) qval
         anion=.true.
       end if
       step=(/ 1.0/real(ndimx),  1.0/real(ndimy), 1.0/real(ndimz) /)
       extend=(/ drmax/cell%cell(1),  drmax/cell%cell(2), drmax/cell%cell(3) /)

       np=0
       !---- Map Calculation ----!
       do k=1,ndimz
          pto(3)=(k-1)*step(3)
          do j=1,ndimy
             pto(2)=(j-1)*step(2)
             do i=1,ndimx
                pto(1)=(i-1)*step(1)

                rx1=pto(1)-extend(1)
                if (rx1 <= 0.0) then
                    nx1=int(rx1)-1
                else
                    nx1=int(rx1)
                end if
                nx2=int(pto(1)+extend(1))

                ry1=pto(2)-extend(2)
                if (ry1 <= 0.0) then
                    ny1=int(ry1)-1
                else
                    ny1=int(ry1)
                end if
                ny2=int(pto(2)+extend(2))

                rz1=pto(3)-extend(3)
                if (rz1 <= 0.0) then
                    nz1=int(rz1)-1
                else
                    nz1=int(rz1)
                end if
                nz2=int(pto(3)+extend(3))

                sbvs=0.0
                rep=0.0
                do n=1,At2%natoms
                   n2=At2%Atom(n)%ind(1)
                   q2=At2%Atom(n)%charge
                   sig2= nint(SIGN(1.0_cp,q2))
                   rho=(radius+A%radius(n2))*0.74
                   !coor=Table_Avcoor(n1,n2)
                   !rzero=Table_Rzero(n1,n2)
                   !rcutoff=Table_Rcutoff(n1,n2)
                   dzero=Table_Dzero(n1,n2)
                    dmin=Table_Rmin(n1,n2)
                   alpha=Table_Alpha(n1,n2)
                   occ=At2%Atom(n)%VarF(1)
                   c_rep=occ*q1*q2
                   c_atr=occ*dzero
                   do k1=nz1,nz2
                      do j1=ny1,ny2
                         do i1=nx1,nx2
                            pta=At2%Atom(n)%x+real((/i1,j1,k1/))
                            dd=Distance(pto,pta,Cell)
                            if (dd > drmax) cycle
                            if (sig1 == sig2) then
                                rep=rep + c_rep*(erfc(dd/rho)-erfc(drmax/rho))/dd
                            else
                               sbvs=sbvs+ c_atr*(exp(alpha*(dmin-dd))-1.0)**2-1.0
                            end if
                         end do
                      end do
                   end do

                end do
                !Multiply the repulsion term by the Coulomb constant to convert to eV
                map_bvs(i,j,k)=sbvs+ke*rep
             end do
          end do
       end do

       !---- Export a File ----!
       !call Get_LogUnit(jbvs)
       open(newunit=jbvs,file=trim(filecod)//".map",status="replace",action="write")

       write (unit=jbvs, fmt='(a)') "BVEL Map Calculations using Bond_STR Program"
       write (unit=jbvs, fmt='(a)') "BVEL Map for species "//trim(car)
       write (unit=jbvs, fmt='(a,3f12.4,3x,3f8.3)') "CELL ",cell%cell,cell%ang
       write (unit=jbvs, fmt='(a)')     "SPGR  "//trim(SpG%spg_symb)
       write(unit=jbvs,fmt="(a,/,a,/)")  &
       " Bond-Valence Energy parameters (D0,Rmin,alpha) for Morse Potential:  D0*{exp(alpha(dmin-d))-1}^2-1", &
             "   (data read from internal table or provided by the user)"
       do n1=1,A%N_Cations
          do j=1,A%N_Anions
             n2=A%N_Cations+j
             write(unit=jbvs,fmt="(2(a,i3,a,a4),/,3(a,f9.5),/,3(a,f9.5),a)")           &
                   "   Type",n1,": ",A%Species(n1)," with type",n2,": ",A%Species(n2), &
                   "    D0  =",Table_Dzero(n1,n2),"       Rmin =",Table_Rmin(n1,n2),   &
                   "  Alpha =",Table_Alpha(n1,n2),"  Av. Coord.=",Table_Avcoor(n1,n2), &
                   "    R0  =",Table_Rzero(n1,n2),"   R-cutoff =",Table_Rcutoff(n1,n2),&
                   "   => Reference: "//trim(references(Table_ref(n1,n2)))
             write(unit=jbvs,fmt="(2(a,a,a,f6.3,a),/)") &
                   "   Cation (Eff. radius): ",A%Species(n1),"(",A%Radius(n1),")   ",  &
                   "   Anion  (Eff. radius): ",A%Species(n2),"(",A%Radius(n2),")"
          end do
       end do
       write (unit=jbvs, fmt='(a)')     "! List ot atoms  "
       do i=1,A%natoms
         write(unit=jbvs, fmt='(a,t20,5f12.5)')"Atom "//trim(A%Atom(i)%lab)//"  "//A%Atom(i)%SfacSymb, &
                                                A%Atom(i)%x, A%Atom(i)%biso, A%Atom(i)%occ
       end do
       write (unit=jbvs, fmt='(a)')     "! Grid values: ndimx,ndimy,ndimz [0,ndimx-1] [0,ndimy-1] [0,ndimz-1]  "
       write (unit=jbvs, fmt='(a,9i6)') "GRID ",ndimx,ndimy,ndimz,0,ndimx-1,0,ndimy-1,0,ndimz-1
       write (unit=jbvs, fmt='(a)')     "DENSITY_MAP"
       write (unit=jbvs,fmt='(10g14.5)') map_bvs
       close(unit=jbvs)
       call write_grid_VESTA(map_bvs,cell,"Bond Valence Energy Landscape",trim(filecod)//"_bvel","P")
       return
    End Subroutine Calc_Map_BVEL

    !!----
    !!---- Subroutine Atlist1_Extencell_Atlist2(Spg,A,B,Conven)
    !!----    type(Space_Group_Type), intent(in)     :: SpG       !  In -> Space Group Information
    !!----    type(atom_list_type),  intent(in)      :: A         !  In -> Atom List (asymmetric unit)
    !!----    type(atom_list_type),  intent(out)     :: B         ! Out -> Atoms in unit cell
    !!----    logical,                intent(in)     :: conven    !  In -> .true. for using the whole conventional unit cell
    !!----
    !!----    Subroutine to generate atoms in the primitive (conven=.false.) or the conventional
    !!----    unit cell (conven=.true.), Excluding atoms with A%atom(:)%active=.false.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine AList2AList_inbox(Spg,A,extend,C)
       !---- Arguments ----!
       type(Space_Group_Type),     intent(in)     :: SpG
       type(atom_list_type),       intent(in)     :: A
       real(kind=cp),dimension(3), intent(in)     :: extend
       type(atom_list_type),       intent(out)    :: C

       !---- Local Variables ----!
       type(atom_list_type)                     :: b
       real(kind=cp),dimension(3)               :: xo,xx,tr
       integer,      dimension(3)               :: imin,imax
       real(kind=cp),dimension(:,:),allocatable :: u
       integer                                  :: k,j,l,nt,npeq,n,i1,i2,i3
       character(len=4)                         :: fmm

       npeq=SpG%multip
       !---- Init proccess ----!
       imax=ceiling(extend)
       imin=-imax
       imax=imax+1
       n=1
       do i=1,3
         n=n*(imax(i)-imin(i)+1)
       end do
       call allocate_atom_list(n*npeq*A%natoms,b)  !n number of tested unit cells
       allocate(u(3,n*npeq))
       n=0
       do k=1,A%natoms
          if (.not. A%atom(k)%active) cycle
          l=1
          n=n+1
          B%Atom(n)=A%Atom(k)
          xo    = modulo_lat(A%atom(k)%x)
          u(:,l)= xo
          B%Atom(n)%x=xo
          do i1=imin(1),imax(1)
            do i2=imin(2),imax(2)
              do i3=imin(3),imax(3)
                tr=real([i1,i2,i3])
                do_eq:do j=2,npeq
                   xx=ApplySO(SpG%SymOp(j),xo)
                   xx=modulo_lat(xx)+tr
                   if(.not. inbox(xx,extend)) cycle
                   do nt=1,l
                      if (equal_vector(u(:,nt),xx,3)) then
                         B%atom(n-(l-nt))%occ=B%atom(n-(l-nt))%occ+A%atom(k)%occ
                         cycle do_eq
                      end if
                   end do
                   l=l+1
                   u(:,l)=xx(:)
                   n=n+1
                   select case (l)
                      case(:9)
                         write(unit=fmm,fmt="(i1)") l
                      case(10:99)
                         write(unit=fmm,fmt="(i2)") l
                      case(100:999)
                         write(unit=fmm,fmt="(i3)") l
                   end select
                   B%Atom(n)=A%Atom(k)
                   B%Atom(n)%lab      =trim(A%Atom(k)%lab)//"_"//adjustl(fmm)
                   B%Atom(n)%x        =xx
                   B%Atom(n)%active   =.true.
                   B%Atom(n)%Mult     = 1.0
                end do do_eq
              end do
            end do
          end do
       end do

       B%natoms=n

       call allocate_atom_list(n,C)

       C%natoms=n
       C%atom(1:n)=B%atom(1:n)

       call deallocate_atom_list(B)

       return

    End Subroutine AList2AList_inbox

    logical function inbox(x,ext)
      real(kind=cp),dimension(3),intent(in) :: x,ext
      inbox=.true.
      do i=1,3
        if(x(i) <= -ext(i) .or. x(i) > ext(i)+1.0) then
          inbox=.false.
          exit
        end if
      end do
    end function inbox

    !!----
    !!---- Subroutine Get_Electronic_Configuration (nam,charge,ElectConf,Qn)
    !!----    character(len=*), intent (in) :: nam        !Chemical Symbol
    !!----    integer,          intent (in) :: charge     !Ionic charge
    !!----    character(len=*), intent(out) :: ElectConf  !Electronic configuration
    !!----    real(kind=cp),    intent(out) :: Qn         !Principal quantum number n of the last occupied shell
    !!----
    !!----    Provides the electronic configuration from the chemical symbol of the element
    !!----    and the valence as an integer. In case of problems the returned electronic configuration is empty.
    !!----
    !!---- Created: January - 2015
    !!
    Subroutine Get_Electronic_Configuration (nam,rcharge,ElectConf,Qn)
       character(len=*), intent (in) :: nam        !Chemical Symbol
       real(kind=cp),    intent (in) :: rcharge    !Ionic charge
       character(len=*), intent(out) :: ElectConf  !Electronic configuration
       real(kind=cp),    intent(out) :: Qn         !Principal quantum number n of the last occupied shell

       !---- Local variables ----!
       character(len=1) :: sh
       character(len=2) :: atm_car
       character(len=6) :: shell
       integer          :: i,Z,le,j,ne,nu,qq, charge

       ElectConf=" "
       Qn=1.0
       atm_car=u_case(nam(1:2))
       z=0
       charge=nint(rcharge)
       if (atm_car(2:2) > "Z" .or. atm_car(2:2) < "A") atm_car(2:2)=" "
       if (.not. allocated(chem_info) ) call set_chem_info()
       do i=1,Num_Chem_Info
          if (index(atm_car,chem_info(i)%Symb) /=0) then
              Z=Chem_Info(i)%Z
              exit
          end if
       end do
       if(z == 0) then
        ElectConf="Error: Z=0"
        return
       end if
       if(.not. allocated(Electronic_Configuration)) call Set_Electronic_Configuration()
       ElectConf=Electronic_Configuration(z)
       le=len_trim(ElectConf)
       i=index(trim(ElectConf)," ",back=.true.)+1
       read(unit=ElectConf(i:i),fmt=*) qq !provisory principal quantum number
       Qn=qq
       read(unit=ElectConf(i+2:),fmt=*) ne
       sh=ElectConf(i+1:i+1)

       Select Case (charge)
         Case(0)
            return
         Case(1:)  !Cation: remove electrons according to the charge
              nu=ne-charge
              Select Case(nu)
                 case(1:)
                   write(unit=ElectConf(i+2:),fmt="(i1)") nu
                 Case(0)
                   ElectConf=ElectConf(1:i-1)
                   i=index(trim(ElectConf)," ",back=.true.)+1
                   read(unit=ElectConf(i:i),fmt=*) Qn
                 Case(:-1)
              End Select
              Select Case (sh)
                case("s")
                  if(nu > 2) then
                    ElectConf(i+2:i+2)="2"
                    ne=nu-2
                    write(unit=shell,fmt="(i1,a,i1)") qq,"p",ne
                    ElectConf=trim(ElectConf)//" "//trim(shell)
                  else
                    write(unit=ElectConf(i+2:),fmt="(i1)") nu
                  end if
                case("p")
                  if(nu > 6) then
                    ElectConf(i+2:i+2)="6"
                    ne=nu-6
                    write(unit=shell,fmt="(i1,a,i1)") qq,"d",ne
                    ElectConf=trim(ElectConf)//" "//trim(shell)
                  else
                    write(unit=ElectConf(i+2:i+2),fmt="(i1)") nu
                  end if
              End select

         Case(:-1) !Anion : Add electrons according to the charge
              nu=ne-charge

              Select Case (sh)
                case("s")
                  if(nu > 2) then
                    ElectConf(i+2:i+2)="2"
                    ne=nu-2
                    write(unit=shell,fmt="(i1,a,i1)") qq,"p",ne
                    ElectConf=trim(ElectConf)//" "//trim(shell)
                  else
                    write(unit=ElectConf(i+2:i+2),fmt="(i1)") nu
                  end if
                case("p")
                  if(nu > 6) then
                    ElectConf(i+2:i+2)="6"
                    ne=nu-6
                    write(unit=shell,fmt="(i1,a,i1)") qq,"d",ne
                    ElectConf=trim(ElectConf)//" "//trim(shell)
                  else
                    write(unit=ElectConf(i+2:i+2),fmt="(i1)") nu
                  end if
              End select
       End Select
       return
    End Subroutine Get_Electronic_Configuration

    !!----
    !!---- Subroutine Set_Electronic_Configuration()
    !!----
    !!----    Electronic configuration of the elements
    !!----
    !!---- Created: January - 2015
    !!
    Subroutine Set_Electronic_Configuration()
       if(.not. allocated(Electronic_Configuration)) then
         allocate(Electronic_Configuration(Num_Chem_Info))
       else
         return
       end if
       Electronic_Configuration(1)   = "1s1"                         !    "H "
       Electronic_Configuration(2)   = "1s2"                         !    "He"
       Electronic_Configuration(3)   = "1s2 2s1"                     !    "Li"
       Electronic_Configuration(4)   = "1s2 2s2"                     !    "Be"
       Electronic_Configuration(5)   = "1s2 2s2 2p1"                 !    "B "
       Electronic_Configuration(6)   = "1s2 2s2 2p2"                 !    "C "
       Electronic_Configuration(7)   = "1s2 2s2 2p3"                 !    "N "
       Electronic_Configuration(8)   = "1s2 2s2 2p4"                 !    "O "
       Electronic_Configuration(9)   = "1s2 2s2 2p5"                 !    "F "
       Electronic_Configuration(10)  = "1s2 2s2 2p6"                 !    "Ne"
       Electronic_Configuration(11)  = "[Ne] 3s1"                    !    "Na"
       Electronic_Configuration(12)  = "[Ne] 3s2"                    !    "Mg"
       Electronic_Configuration(13)  = "[Ne] 3s2 3p1"                !    "Al"
       Electronic_Configuration(14)  = "[Ne] 3s2 3p2"                !    "Si"
       Electronic_Configuration(15)  = "[Ne] 3s2 3p3"                !    "P "
       Electronic_Configuration(16)  = "[Ne] 3s2 3p4"                !    "S "
       Electronic_Configuration(17)  = "[Ne] 3s2 3p5"                !    "Cl"
       Electronic_Configuration(18)  = "[Ne] 3s2 3p6"                !    "Ar"
       Electronic_Configuration(19)  = "[Ar] 4s1"                    !    "K "
       Electronic_Configuration(20)  = "[Ar] 4s2"                    !    "Ca"
       Electronic_Configuration(21)  = "[Ar] 3d1 4s2"                !    "Sc"
       Electronic_Configuration(22)  = "[Ar] 3d2 4s2"                !    "Ti"
       Electronic_Configuration(23)  = "[Ar] 3d3 4s2"                !    "V "
       Electronic_Configuration(24)  = "[Ar] 3d5 4s1"                !    "Cr"
       Electronic_Configuration(25)  = "[Ar] 3d5 4s2"                !    "Mn"
       Electronic_Configuration(26)  = "[Ar] 3d6 4s2"                !    "Fe"
       Electronic_Configuration(27)  = "[Ar] 3d7 4s2"                !    "Co"
       Electronic_Configuration(28)  = "[Ar] 3d8 4s2"                !    "Ni"
       Electronic_Configuration(29)  = "[Ar] 3d10 4s1"               !    "Cu"
       Electronic_Configuration(30)  = "[Ar] 3d10 4s2"               !    "Zn"
       Electronic_Configuration(31)  = "[Ar] 3d10 4s2 4p1"           !    "Ga"
       Electronic_Configuration(32)  = "[Ar] 3d10 4s2 4p2"           !    "Ge"
       Electronic_Configuration(33)  = "[Ar] 3d10 4s2 4p3"           !    "As"
       Electronic_Configuration(34)  = "[Ar] 3d10 4s2 4p4"           !    "Se"
       Electronic_Configuration(35)  = "[Ar] 3d10 4s2 4p5"           !    "Br"
       Electronic_Configuration(36)  = "[Ar] 3d10 4s2 4p6"           !    "Kr"
       Electronic_Configuration(37)  = "[Kr] 5s1"                    !    "Rb"
       Electronic_Configuration(38)  = "[Kr] 5s2"                    !    "Sr"
       Electronic_Configuration(39)  = "[Kr] 4d1 5s2"                !    "Y "
       Electronic_Configuration(40)  = "[Kr] 4d2 5s2"                !    "Zr"
       Electronic_Configuration(41)  = "[Kr] 4d4 5s1"                !    "Nb"
       Electronic_Configuration(42)  = "[Kr] 4d5 5s1"                !    "Mo"
       Electronic_Configuration(43)  = "[Kr] 4d5 5s2"                !    "Tc"
       Electronic_Configuration(44)  = "[Kr] 4d7 5s1"                !    "Ru"
       Electronic_Configuration(45)  = "[Kr] 4d8 5s1"                !    "Rh"
       Electronic_Configuration(46)  = "[Kr] 4d10"                   !    "Pd"
       Electronic_Configuration(47)  = "[Kr] 4d10 5s1"               !    "Ag"
       Electronic_Configuration(48)  = "[Kr] 4d10 5s2"               !    "Cd"
       Electronic_Configuration(49)  = "[Kr] 4d10 5s2 5p1"           !    "In"
       Electronic_Configuration(50)  = "[Kr] 4d10 5s2 5p2"           !    "Sn"
       Electronic_Configuration(51)  = "[Kr] 4d10 5s2 5p3"           !    "Sb"
       Electronic_Configuration(52)  = "[Kr] 4d10 5s2 5p4"           !    "Te"
       Electronic_Configuration(53)  = "[Kr] 4d10 5s2 5p5"           !    "I "
       Electronic_Configuration(54)  = "[Kr] 4d10 5s2 5p6"           !    "Xe"
       Electronic_Configuration(55)  = "[Xe] 6s1"                    !    "Cs"
       Electronic_Configuration(56)  = "[Xe] 6s2"                    !    "Ba"
       Electronic_Configuration(57)  = "[Xe] 5d1 6s2"                !    "La"
       Electronic_Configuration(58)  = "[Xe] 4f1 5d1 6s2"            !    "Ce"
       Electronic_Configuration(59)  = "[Xe] 4f3 5d0 6s2"            !    "Pr"
       Electronic_Configuration(60)  = "[Xe] 4f4 5d0 6s2"            !    "Nd"
       Electronic_Configuration(61)  = "[Xe] 4f5 5d0 6s2"            !    "Pm"
       Electronic_Configuration(62)  = "[Xe] 4f6 5d0 6s2"            !    "Sm"
       Electronic_Configuration(63)  = "[Xe] 4f7 5d0 6s2"            !    "Eu"
       Electronic_Configuration(64)  = "[Xe] 4f7 5d1 6s2"            !    "Gd"
       Electronic_Configuration(65)  = "[Xe] 4f9 5d0 6s2"            !    "Tb"
       Electronic_Configuration(66)  = "[Xe] 4f10 5d0 6s2"           !    "Dy"
       Electronic_Configuration(67)  = "[Xe] 4f11 5d0 6s2"           !    "Ho"
       Electronic_Configuration(68)  = "[Xe] 4f12 5d0 6s2"           !    "Er"
       Electronic_Configuration(69)  = "[Xe] 4f13 5d0 6s2"           !    "Tm"
       Electronic_Configuration(70)  = "[Xe] 4f14 5d0 6s2"           !    "Yb"
       Electronic_Configuration(71)  = "[Xe] 4f14 5d1 6s2"           !    "Lu"
       Electronic_Configuration(72)  = "[Xe] 4f14 5d2 6s2"           !    "Hf"
       Electronic_Configuration(73)  = "[Xe] 4f14 5d3 6s2"           !    "Ta"
       Electronic_Configuration(74)  = "[Xe] 4f14 5d4 6s2"           !    "W "
       Electronic_Configuration(75)  = "[Xe] 4f14 5d5 6s2"           !    "Re"
       Electronic_Configuration(76)  = "[Xe] 4f14 5d6 6s2"           !    "Os"
       Electronic_Configuration(77)  = "[Xe] 4f14 5d7 6s2"           !    "Ir"
       Electronic_Configuration(78)  = "[Xe] 4f14 5d9 6s1"           !    "Pt"
       Electronic_Configuration(79)  = "[Xe] 4f14 5d10 6s1"          !    "Au"
       Electronic_Configuration(80)  = "[Xe] 4f14 5d10 6s2"          !    "Hg"
       Electronic_Configuration(81)  = "[Xe] 4f14 5d10 6s2 6p1"      !    "Tl"
       Electronic_Configuration(82)  = "[Xe] 4f14 5d10 6s2 6p2"      !    "Pb"
       Electronic_Configuration(83)  = "[Xe] 4f14 5d10 6s2 6p3"      !    "Bi"
       Electronic_Configuration(84)  = "[Xe] 4f14 5d10 6s2 6p4"      !    "Po"
       Electronic_Configuration(85)  = "[Xe] 4f14 5d10 6s2 6p5"      !    "At"
       Electronic_Configuration(86)  = "[Xe] 4f14 5d10 6s2 6p6"      !    "Rn"
       Electronic_Configuration(87)  = "[Rn] 7s1"                    !    "Fr"
       Electronic_Configuration(88)  = "[Rn] 7s2"                    !    "Ra"
       Electronic_Configuration(89)  = "[Rn] 6d2 7s1"                !    "Ac"
       Electronic_Configuration(90)  = "[Rn] 5f0 6d2 7s2"            !    "Th"
       Electronic_Configuration(91)  = "[Rn] 5f2 6d1 7s2"            !    "Pa"
       Electronic_Configuration(92)  = "[Rn] 5f3 6d1 7s2"            !    "U "
       Electronic_Configuration(93)  = "[Rn] 5f4 6d1 7s2"            !    "Np"
       Electronic_Configuration(94)  = "[Rn] 5f6 6d0 7s2"            !    "Pu"
       Electronic_Configuration(95)  = "[Rn] 5f7 6d0 7s2"            !    "Am"
       Electronic_Configuration(96)  = "[Rn] 5f7 6d1 7s2"            !    "Cm"
       Electronic_Configuration(97)  = "[Rn] 5f8 6d1 7s2"            !    "Bk"
       Electronic_Configuration(98)  = "[Rn] 5f10 6d0 7s2"           !    "Cf"
       Electronic_Configuration(99)  = "[Rn] 5f11 6d0 7s2"           !    "Es"
       Electronic_Configuration(100) = "[Rn] 5f12 6d0 7s2"           !    "Fm"
       Electronic_Configuration(101) = "[Rn] 5f13 + 6d0 7s2"         !    "Md"
       Electronic_Configuration(102) = "[Rn] 5f14 + 6d0 7s2"         !    "No"
       Electronic_Configuration(103) = "[Rn] 5f14 + 6d1 7s2"         !    "Lw"
       Electronic_Configuration(104) = "[Rn] 5f14 + 6d2 7s2"         !    "Rf"
       Electronic_Configuration(105) = "[Rn] 5f14 + 6d3 7s2"         !    "Db"
       Electronic_Configuration(106) = "[Rn] 5f14 + 6d4 7s2"         !    "Sg"
       Electronic_Configuration(107) = "[Rn] 5f14 + 6d5 7s2"         !    "Bh"
       Electronic_Configuration(108) = "[Rn] 5f14 + 6d6 7s2"         !    "Hs"
    End Subroutine Set_Electronic_Configuration





End Program Bond_Str

