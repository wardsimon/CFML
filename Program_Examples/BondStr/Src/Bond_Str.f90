!!----
!!---- Program Bond_STR
!!----
!!---- Using CrysFML versions above 4.01
!!----
Program Bond_Str
   !---- Use Modules ----!
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
   use CFML_IO_Formats,                  only: Readn_set_Xtal_Structure,ERR_Form_Mess,err_form,file_list_type,   &
                                               Write_Cif_Template,Write_CFL
   Use CFML_Scattering_Chemical_Tables
   Use CFML_Export_VTK
   use CFML_BVSpar
   use CFML_BVS_Energy_Calc
   use CFML_Percolation

   !---- Variables ----!
   Implicit None

   type (Space_Group_Type)         :: SpGr
   type (Crystal_Cell_Type)        :: Cell
   type (Atom_list_Type)           :: A
   type (Atoms_Conf_list_Type)     :: Ac
   type (File_List_Type)           :: Fich_cfl

   character(len=256)              :: filcod,restr_file
   character(len=256)              :: line,title,lineor
   character(len=80),dimension(15) :: bvparm,bvelparm,fst_cmd
   character(len=4)                :: atname,atm
   character(len=2)                :: Chem

   integer                         :: lun=1, ier,i, lr,ln, i_cfl=2, i_cons=3
   integer                         :: narg,npix
   integer                         :: nc,n_bvpar,n_bvelpar,n1,n2,j,n_fst
   integer                         :: ndimx,ndimy,ndimz

   logical                         :: esta, arggiven=.false.,sout=.false.,cif=.false.,out_cif=.false.
   logical                         :: read_bvparm=.false., restr=.false., bvs_calc=.true., percolation=.false.
   logical                         :: vdist=.false.,read_bvelparm=.false.,rest_file=.false.
   logical                         :: map_calc=.false., soft=.false., bvel_calc=.false., outp=.false., outf=.false.
   real                            :: ttol=20.0,dmax,dangl, rdmax,ramin
   real                            :: drmax,delta,qval,tini,tfin,qn,qp,vol 
   
   Real(kind=cp) Emin,Eini,Eend,dE,dE_ini
   Real(kind=cp), Dimension(3) :: E_percol,E_percol_aux
   Real(kind=cp), Dimension(:,:,:), Allocatable :: bvel_map

   ! Arguments on the command line
   lr=0
   ln=0
   narg=Command_Argument_Count()

   if (narg > 0) then
      call GET_COMMAND_ARGUMENT(1,filcod)
      i=index(filcod,".cfl")
      if(i /= 0) filcod=filcod(1:i-1)
      i=index(filcod,".cif")
      if(i /= 0) filcod=filcod(1:i-1)
      arggiven=.true.
   end if

   write(unit=*,fmt="(/,/,8(a,/))")                                                        &
           "                      =============================="                        , &
           "                      ====== PROGRAM BOND_STR ======"                        , &
           "                      =============================="                        , &
           "    ***********************************************************************" , &
           "    * Distances, angles and Bond-Valence Sums from  *.cfl or *.cif files  *" , &
           "    *     Calculation of BVS and Bond-Valence Energy Landscape maps       *" , &
           "    ***********************************************************************" , &
           "                     (JRC - ILL, version: June 2017 )"
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
      out_cif=.true.
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
   write(unit=lun,fmt="(/,/,8(a,/))")                                                  &
           "                      =============================="                        , &
           "                      ====== PROGRAM BOND_STR ======"                        , &
           "                      =============================="                        , &
           "    ***********************************************************************" , &
           "    * Distances, angles and Bond-Valence Sums from  *.cfl or *.cif files  *" , &
           "    *     Calculation of BVS and Bond-Valence Energy Landscape maps       *" , &
           "    ***********************************************************************" , &
           "                     (JRC - ILL, version: June 2017 )"
   write(unit=lun,fmt="(a,/)") " "
   dmax=3.2
   dangl=0.0
   title=" "
   delta=0.0
   rdmax=2.5
   ramin=45.0

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
   !Calculation of the total charge in the asymmetric unit
   if(bvs_calc) then
     qp=0.0; qn=0.0
     do i=1,A%natoms
       if(A%atom(i)%charge > 0.0) then
        qp=qp+A%atom(i)%charge*A%atom(i)%occ
       else
        qn=qn+A%atom(i)%charge*A%atom(i)%occ
       end if
     end do
     write(unit=lun,fmt="(/,a,f10.4)") " => Weighted positive charge of the asymmetric unit: ",qp
     write(unit=lun,fmt="(  a,f10.4)") " => Weighted negative charge of the asymmetric unit: ",qn
     if(abs(qp+qn) > 0.1) then
       write(unit=lun,fmt="(  a,f10.4)") " => WARNING!! Electro-neutrality is not satisfied (q > 0.1), total charge: ",qn+qp
     else
       write(unit=lun,fmt="(  a,f10.4)") " => Electro-neutrality satisfied (q <= 0.1), total charge: ",qn+qp
     end if
   end if

   ! Look for specific instructions for BOND_STR in CFL-file
   ttol=30.0
   if (cif) then
      sout=.true.
   else
      n_bvpar=0; n_bvelpar=0; n_fst=0
      bvparm=" "; bvelparm=" "
      do i=1,fich_cfl%nlines
         lineor=adjustl(fich_cfl%line(i))
         line=u_case(lineor)

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

         if (line(1:7) == "FST_CMD") then
            n_fst=n_fst+1
            nc=index(lineor,' ')
            fst_cmd(n_fst)=adjustl(lineor(nc+1:))
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

         if (line(1:6) == "RESTDA") then
            read(unit=line(7:),fmt=*,iostat=ier) rdmax,ramin
            if (ier /= 0) then
               rdmax=2.5
               ramin=45.0
            end if
            cycle
         end if

         if(line(1:10) == "RESTR_FILE") then
          rest_file=.true.
          nc=index(lineor,' ')
          restr_file=adjustl(lineor(nc+1:))
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

         if (line(1:6) == "OUTMAP")  then
            outp=.true.
            if(index(line,"GFOU") /= 0) outf=.true.
         end if
         
         if (line(1:11) == "PERCOLATION")  then
            percolation=.true.
         end if
         
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
   if( .not. soft .and. map_calc ) soft=.true.
   if(out_cif) then !Generate a CIF file from the CFL file in order to write a complete VESTA file
     call Write_Cif_Template(trim(filcod)//"_gen.cif",2,title,Cell,SpGr,A)
   end if

   ! Distances, Angles, Restraints Calculations
   if (restr) then
      if(rest_file) then
        write(*,*) trim(restr_file)//"   => rdmax,ramin: ",rdmax,ramin
        call Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spgr, A, lun, i_cons,filrest=restr_file,rdmax=rdmax,ramin=ramin)
      else
        call Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spgr, A, lun, i_cons,rdmax=rdmax,ramin=ramin)
      end if
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
   if (vdist) call Calc_PDB(A,Dmax,lun)

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
        !if(soft) then
          call Species_on_List(Ac,SpGr%Multip,ttol,.true.,soft) !using softBVS for determining Morse parameters
        !else
        !  call Species_on_List(Ac,SpGr%Multip,ttol,.true.) !using covalent radius and table for oxides
        !end if
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
         write(unit=lun,fmt="(/,a)") "    ---------------------------------------------------"
         write(unit=lun,fmt="(a)")   "       Calculation of Bond-Valence Energy Landscape Map"
         write(unit=lun,fmt="(a,/)") "    ---------------------------------------------------"
         if(soft) then
             if(read_bvelparm) then
                call Set_Table_BVEL_Params(Ac,N_bvelpar,bvelparm,soft)
             else
                call Set_Table_BVEL_Params(Ac,soft=soft)
             end if
         else
             if(read_bvelparm) then
                call Set_Table_BVEL_Params(Ac,N_bvelpar,bvelparm)
             else
                call Set_Table_BVEL_Params(Ac)
             end if
         end if
         if(.not. err_conf) then
            write(unit=lun,fmt="(/a,/,a,/)")  &
            " Bond-Valence Energy parameters (D0,Rmin,alpha) for Morse Potential:  D0*[{exp(alpha(dmin-d))-1}^2-1]", &
                  "   (data read from internal table, provided by the user or calculated from softBVS parameters)"

            if(soft) then
              write(unit=lun,fmt="(a,/)") " Morse parameters D0 and Rmin calculated from softBVS parameters"
            end if
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
         end if
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

      ! BVS
      ! Maps of BVS or BVEL
      if (map_calc) then

         write(unit=lun,fmt="(/,a)") "    ----------------------------------------------------------"
         write(unit=lun,fmt="(a)")   "       Calculation of Bond-Valence Map with softBVS parameters"
         write(unit=lun,fmt="(a,/)") "    ----------------------------------------------------------"
         write (unit=lun, fmt='(/,a,f10.4,a)')" => Global distance cutoff:",drmax," angstroms"
         if(delta > 0.001) then
           call Calc_Map_BVS(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,delta,vol)
           write (unit=lun, fmt='(a,f10.4,a)')  " => Delta (Output for valence +/- delta):",delta," valence units"
           write (unit=lun, fmt='(a,f10.4,a)')  " => Available volume for ion mobility in the unit cell:",vol," angstroms^3"
           write (unit=lun, fmt='(a,f10.2,a)')  " => Volume  fraction for ion mobility in the unit cell:",vol/Cell%CellVol*100.0, " %"
         else
           call Calc_Map_BVS(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax)
         end if

      else if(bvel_calc) then
         Write(unit=*,fmt="(/,a)") " => Calculation of Bond-Valence Energy Landscape Map (it can take some minutes) ...."
         write (unit=lun, fmt='(/,a,f10.4,a)')" => Global distance cutoff:",drmax," angstroms"
         if(delta > 0.01) then
           if(outp) then
           	 if(percolation) then
               call Calc_Map_BVEL(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,delta,vol,emin,npix,outf,bvel_map)
           	 else
               call Calc_Map_BVEL(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,delta,vol,emin,npix,outf)
             end if
           else
           	 if(percolation) then
               call Calc_Map_BVEL(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,delta,vol,emin,npix,bvel_map=bvel_map)
           	 else
               call Calc_Map_BVEL(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,delta,vol,emin,npix)
             end if
           end if
           write (unit=lun, fmt='(/,a,f10.4,a)')" => Value of Delta (for volume calculation) :",delta," eV"
           write (unit=lun, fmt='(a,f10.4,a)')  " => Available volume for ion mobility in the unit cell:",vol," angstroms^3"
           write (unit=lun, fmt='(a,f10.2,a)')  " => Volume  fraction for ion mobility in the unit cell:",vol/Cell%CellVol*100.0, " %"
           write (unit=lun, fmt='(a,f10.4)')    " => Minum Energy (in eV):", emin
           write (unit=lun, fmt='(a,i8)')       " => Number of pixels with Emin < Energy < Emin+Delta: ",npix
           
           if(percolation) then           	
           	 call Percolation_Calc()
           end if
        
         else
           if(outp) then
             call Calc_Map_BVEL(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax,outp=outf)
           else
             call Calc_Map_BVEL(Ac,Spgr,Cell,trim(filcod),ndimx,ndimy,ndimz,atname,drmax)
           end if
         end if
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
   write(unit=*,fmt="(a)")     " => Global results in File: "//trim(filcod)//".bvs"
   if(bvs_calc) then
     if (map_calc) then
       write(unit=*,fmt="(a)")     " => Bond Valence Map in File: "//trim(filcod)//".map"
       write(unit=lun,fmt="(a)")   " => Bond Valence Map in File: "//trim(filcod)//".map"
     else if (bvel_calc) then
       write(unit=*,fmt="(a)")     " => Bond Valence Energy Landscape in File: "//trim(filcod)//"_bvel.map"
       write(unit=lun,fmt="(a)")   " => Bond Valence Energy Landscape in File: "//trim(filcod)//"_bvel.map"
     else
       write(unit=*,fmt="(a)")     " => Summary of BVS in File: "//trim(filcod)//"_sum.bvs"
       write(unit=lun,fmt="(a)")   " => Summary of BVS in File: "//trim(filcod)//"_sum.bvs"
     end if
   end if
   if (map_calc  .or. bvel_calc) then
     call Write_Vesta_File()
     write(unit=*,fmt="(a)")     " => VESTA File: "//trim(filcod)//"_str.vesta"
     write(unit=lun,fmt="(a)")   " => VESTA File: "//trim(filcod)//"_str.vesta"
   end if
   call cpu_time(tfin)
   write(unit=lun,fmt="(/,a)")        " => Normal End of: PROGRAM BOND_STR "
   write(unit=*,fmt="(a,f8.4,a)")     " => CPU-time: ",tfin-tini," seconds"
   write(unit=lun,fmt="(a,f8.4,a)")   " => CPU-time: ",tfin-tini," seconds"
   close(unit=lun)
   call Deallocate_Atoms_Conf_List(Ac)

 Contains

   Subroutine Percolation_Calc()
   	  real :: Eminim
      Character(len=*), Dimension(3), parameter :: axis=["a","b","c"]
      bvel_map    = bvel_map - Emin
      Eminim      = 0.0
      Eini        = 0.0
      Eend        = 3.0
      dE          = 0.5
      E_percol(:) = -1.
      
      ! Percolation analysis
      
      Write(unit=*,fmt="(/,a)") " => Calculation of percolation energies (it can take some minutes) ...."
      Write(unit=lun,fmt="(/,a)") " => Computing first estimation of percolation energies (it can take some minutes) ...."
      
      Call Percol_Analysis(bvel_map,Eminim,Eini,Eend,dE,E_percol)
      
      Do i = 1 , 3
         If (E_percol(i) > 0.0) Then
            Write(unit=lun,fmt="(tr4,3a,f6.2,a)") "Percolation along ", axis(i), ": Yes, Percolation energy: ", E_percol(i), " eV"
         Else
            Write(unit=lun,fmt="(tr4,3a)") "Percolation along ", axis(i), ": No"
         End If
      End Do
      
      Write(unit=lun,fmt="(/,a)") " => Refining energies...."
      
      dE_ini = dE
      Do i = 1 , 3
         dE = dE_ini
         If (E_percol(i) > 0.0) Then
            Write(unit=lun,fmt="(tr4,2a)") "axis ", axis(i)
            Eini = E_percol(i) - dE
            dE   = 0.1
            Eend = E_percol(i) + 0.01
            Write(unit=lun,fmt="(tr6,a,2(f5.2,a))") "Searching percolation between ",Eini," and ",Eend," eV"           
            Call Percol_Analysis(bvel_map,Eminim,Eini,Eend,dE,E_percol_aux,axis=i)
            E_percol(i) = E_percol_aux(i)
            Write(unit=lun,fmt="(tr8,a,f6.2,a)") "Percolation energy above Emin: ", E_percol(i), " eV"
            Eini = E_percol(i) - dE
            dE   = 0.01
            Eend = E_percol(i) + 0.01
            Write(unit=lun,fmt="(tr6,a,2(f5.2,a))") "Searching percolation between ",Eini," and ",Eend," eV"
            Call Percol_Analysis(bvel_map,Eminim,Eini,Eend,dE,E_percol_aux,axis=i)
            E_percol(i) = E_percol_aux(i)
            Write(unit=lun,fmt="(tr8,a,2(f6.2,a),/)") "Percolation energy above Emin: ", E_percol(i), " eV,   Isosurface for VESTA: ", E_percol(i)+Emin," eV"
         End If
      End Do
   	
   End Subroutine Percolation_Calc
   
   
   Subroutine Write_Vesta_File()

      integer                                     :: lun, i, j, n,np, pol,k,cent
      character(len=2)                            :: elem
      character(len=80)                           :: box_cmd,aux
      character(len=2), dimension(:), allocatable :: poly
      character(len=80),dimension(:), allocatable :: cmd_bond
      character(len=80)  :: cmd

      open(newunit=lun,file=trim(filcod)//"_str.vesta",action="write",status="replace")
      write(unit=lun,fmt="(a)") "#VESTA_FORMAT_VERSION 3.1.9"
      write(unit=lun,fmt="(a)") " "
      write(unit=lun,fmt="(a)") " "
      write(unit=lun,fmt="(a)") "CRYSTAL"
      write(unit=lun,fmt="(a)") "TITLE"
      write(unit=lun,fmt="(a)") "  "//trim(title(5:))
      write(unit=lun,fmt="(a)") " "
      write(unit=lun,fmt="(a)") "IMPORT_STRUCTURE"
      write(unit=lun,fmt="(a)") trim(filcod)//"_gen.cif"
      write(unit=lun,fmt="(a)") " "
      write(unit=lun,fmt="(a)") "IMPORT_DENSITY 1"
      if (map_calc) then
        write(unit=lun,fmt="(a)") "+1.000000E+000 "//trim(filcod)//"_bvs.pgrid"
      else if (bvel_calc) then
        write(unit=lun,fmt="(a)") "+1.000000E+000 "//trim(filcod)//"_bvel.pgrid"
      end if
      write(unit=lun,fmt="(a)") " "
      write(unit=lun,fmt="(a)") "BOUND"
      box_cmd=" "
      do i=1,n_fst
        aux=u_case(fst_cmd(i))
        j=index(aux,"BOX")
        if(j /= 0) then
          box_cmd=fst_cmd(i)(j+3:)
          exit
        end if
      end do
      if(len_trim(box_cmd) == 0) then
         write(unit=lun,fmt="(a)") "    -0.1      1.1      -0.1      1.1      -0.1     1.1 "
      else
         write(unit=lun,fmt="(a)") trim(box_cmd)
      end if
      write(unit=lun,fmt="(a)")"  0   0   0   0  0"

      write(unit=lun,fmt="(a)") "SBOND"

      n=0; np=0
      pol=0
      allocate(cmd_bond(A%natoms), poly(A%natoms))
      cmd_bond=" "
      poly=" "
      do i=1,n_fst
          aux=u_case(fst_cmd(i))
          j=index(aux,"CONN")
          if(j /= 0) then
            n=n+1
            cmd=adjustl(fst_cmd(i)(j+4:))
            j=index(cmd," ")
            elem=cmd(1:j-1)
            cmd=adjustl(cmd(j:))
            j=index(cmd," ")
            cmd_bond(n)=elem//"  "//cmd(1:j)//"  "//cmd(j+1:)
            cycle
          end if
          j=index(aux,"POLY")
          if(j /= 0)  then
             pol=1
             np=np+1
             poly(np)=adjustl(fst_cmd(i)(j+4:))
          end if
      end do
      do i=1,n
        j=index(cmd_bond(i)," ")
        elem=cmd_bond(i)(1:j-1)
        cent=0
        do k=1,np
          if(trim(poly(k)) == trim(elem)) then
             cent=1
             exit
          end if
        end do
        write(unit=lun,fmt="(i3,a,5i3)") i,"  "//trim(cmd_bond(i)),0,1,cent,0,1
      end do
      write(unit=lun,fmt="(a)")    "  0 0 0 0"
      write(unit=lun,fmt="(a)")    " "
      write(unit=lun,fmt="(a)")    "STYLE"
      write(unit=lun,fmt="(a)")    "MODEL   2  1  0"
      write(unit=lun,fmt="(a)")    "SURFS   0  1  1"
      write(unit=lun,fmt="(a)")    "SECTS  96  0"
      write(unit=lun,fmt="(a,i1)") "POLYS  ",pol
      call flush(lun)
      close(unit=lun)
   End Subroutine Write_Vesta_File
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

End Program Bond_Str

