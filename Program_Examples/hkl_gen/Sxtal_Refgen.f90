Module Ref_Gen

    use CFML_IO_Formats,      only: file_list_type
    use CFML_Math_General,    only: sind,cosd,acosd,asind,Co_Prime_Vector
    use CFML_ILL_Instrm_data, only: Err_ILLdata_Mess, Err_ILLdata, &
                                    SXTAL_Orient_type, Current_Orient, diffractometer_type, &
                                    Current_Instrm, Read_Current_Instrm, Update_Current_Instrm_UB,&
                                    Set_default_Instrument,Write_Current_Instrm_data
    use CFML_String_Utilities,only: l_case, Get_LogUnit
    use CFML_crystal_metrics, only: Crystal_Cell_Type,Set_Crystal_Cell, Write_Crystal_Cell,Zone_Axis_type, &
                                    Get_basis_from_uvw
    use CFML_Geometry_SXTAL

    implicit none

    private
    public :: read_sxtal_geom,  calc_angles

 contains

    Subroutine read_sxtal_geom(ipr,cell,file_dat,ok,mess,ord,spher,schwinger,hlim,sgiven,smin,smax,opgiven,iop_ord)
        !---- Arguments ----!
        integer,                intent( in)    :: ipr
        type(Crystal_Cell_Type),intent(in out) :: cell
        Type(file_list_type),   intent( in)    :: file_dat
        logical,                intent(out)    :: ok
        character(len=*),       intent(out)    :: mess
        integer,dimension(3),   intent(out)    :: ord
        logical,                intent(out)    :: spher
        logical,                intent(out)    :: schwinger
        integer,dimension(3,2), intent(out)    :: hlim
        logical, optional,      intent(out)    :: sgiven,opgiven
        real,    optional,      intent(out)    :: smin,smax
        integer, optional,      intent(out)    :: iop_ord
        !---- Local variables ----!
        character(len=132)   :: line, file_inst
        real                 :: wave,wav,tmin,tmax,omeg,s1,s2,ang12,tet1,tet2
        real, dimension(3,3) :: ub
        integer, dimension(3):: uvw
        real,    dimension(3):: h1,h2,hv
        real, dimension(6)   :: dcel,incel
        integer              :: i,ier,j, n, lun, i_def, igeom
        logical              :: esta, ub_read, wave_read, inst_read, trang_read
        type(Zone_Axis_Type) :: zone_axis

        ok=.false.
        ub_read=.false.
        wave_read=.false.
        inst_read=.false.
        spher=.false.
        schwinger=.false.
        trang_read=.false.
        igeom=0
        mess= " "
        ord=(/3,2,1/)
        incel(1:3)=cell%cell
        incel(4:6)=cell%ang
        ub=reshape((/1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0/),(/3,3/))
        hlim=0
        if(present(opgiven)) opgiven=.false. !Initialise to false, put to .true. if ANGOR is provided

        do j=1,file_dat%nlines
            line=adjustl(file_dat%line(j))
            line=l_case(line)
            if (line(1:1) ==" ") cycle
            if (line(1:1) =="!") cycle
            i=index(line,"!")
            if( i /= 0) line=line(1:i-1)
            i=index(line," ")-1
            select case (line(1:i))

                case("srang")
                    if(present(sgiven) .and. present(smin) .and. present(smax)) then
                        read(unit=line(6:),fmt=*,iostat=ier) smin,smax
                        if(ier /= 0) then
                            mess="Error reading Smin and Smax in CFL file: "//trim(file_inst)
                            sgiven=.false.
                            return
                        else
                            sgiven=.true.
                        end if
                    end if

                case("trang")
                    if(present(sgiven) .and. present(smin) .and. present(smax) .and. wave_read) then
                        read(unit=line(6:),fmt=*,iostat=ier) tmin,tmax
                        if(ier /= 0) then
                            mess="Error reading 2theta_min and 2Theta_max in CFL file: "//trim(file_inst)
                            sgiven=.false.
                            return
                        else
                            smin=sind(0.5*tmin)/wave
                            smax=sind(0.5*tmax)/wave
                            sgiven=.true.
                            trang_read=.true.
                        end if
                    end if

                case("hlim")
                    ! hmin hmax  kmin kmax  lmin lmax
                    read(unit=line(6:),fmt=*,iostat=ier) ((hlim(i,n),n=1,2),i=1,3)
                    if(ier /= 0 ) then
                        mess="Error reading hkl-limits in CFL file: "//trim(file_inst)
                        return
                    end if

                case("instr")
                    line=file_dat%line(j)
                    i=index(line,"!")
                    if(i /= 0) line=line(1:i-1)
                    i=index(trim(line)," ",back=.true.)
                    file_inst=trim(adjustl(line(i:)))
                    inquire(file=trim(file_inst), exist=esta)
                    if(.not. esta) then
                        mess="File "//trim(file_inst)//" not found!"
                        return
                    end if
                    call Read_Current_Instrm(trim(file_inst))
                    ub=Current_Orient%ub
                    inst_read=.true.

                case("ubmat")
                    do i=1,3
                        read(unit=file_dat%line(j+i),fmt=*,iostat=ier) ub(i,:)
                        if(ier /= 0 ) then
                            mess="Error reading the orientation matrix in CFL file: "//trim(file_inst)
                            return
                        end if
                    end do
                    ub_read=.true.

                case("geom")
                    read(unit=line(6:),fmt=*,iostat=ier) igeom
                    if(ier /= 0 ) then
                        mess="Error reading the diffraction geometry in CFL file: "//trim(file_inst)
                        return
                    end if

                case("wave")
                    read(unit=line(6:),fmt=*,iostat=ier) wave
                    if(ier /= 0 ) then
                        mess="Error reading the wavelength in CFL file: "//trim(file_inst)
                        return
                    end if
                    wave_read=.true.

                case("order")
                    read(unit=line(6:),fmt=*,iostat=ier) ord
                    if(ier /= 0 ) then
                        mess="Error reading hkl-ordering in CFL file: "//trim(file_inst)
                        return
                    end if

                case("angor")

                    if(present(iop_ord)) then
                        iop_ord=0
                        i=index(line,"2theta")
                        if(i /= 0) iop_ord=1
                        i=index(line,"gamma")
                        if(i /= 0) iop_ord=1
                        i=index(line,"omega")
                        if(i /= 0) iop_ord=2
                        i=index(line,"chi")
                        if(i /= 0) iop_ord=3
                        i=index(line,"nu")
                        if(i /= 0) iop_ord=3
                        i=index(line,"phi")
                        if(i /= 0) iop_ord=4

                        if(present(opgiven)) opgiven=.true.

                    end if

                case("spher")
                    spher=.true.

                case("schwinger")
                    schwinger=.true.

                case("orient_hh")
                    if(.not. wave_read .and. .not. inst_read) then
                      write(unit=*,fmt="(a)") " => Wavelength shoud be provided!"
                      return
                    else if(.not. wave_read) then
                      wave=Current_Orient%wave
                    end if
                    read(unit=line(10:),fmt=*) h1,h2,omeg
                    hv=cross_product(h1,h2)
                    call Co_Prime_Vector(nint(hv),uvw)
                    call Get_basis_from_uvw(0.5,uvw,cell,Zone_Axis,ok)
                    if(ok) then
                      Call Get_UB_from_uvw_hkl_omega(wave,Cell,Zone_Axis,h1,omeg,UB,ok,mess)
                      if(.not. ok) then
                        write(unit=*,fmt="(a)") " => UB cannot be calculated from ORIENT_HH instruction "
                        write(unit=*,fmt="(a)") " => "//trim(mess)
                        return
                      end if
                    else
                      write(unit=*,fmt="(a)") " => UB cannot be calculated from ORIENT_HH instruction "
                      write(unit=*,fmt="(a,3i4,tr3,3i4,a,f7.2,a,3i4)") " => Problem with reflections:",nint(h1),nint(h2), " with omega(h1)=",omeg," and calculated [uvw]:",uvw
                      return
                    end if
                    ub_read=.true.

                case("orient_vh")
                    if(.not.wave_read .and. .not. inst_read) then
                      write(unit=*,fmt="(a)") " => Wavelength shoud be provided!"
                      return
                    else if(.not. wave_read) then
                      wave=Current_Orient%wave
                    end if

                    read(unit=line(10:),fmt=*) uvw,h1,omeg
                    call Get_basis_from_uvw(0.5,uvw,cell,Zone_Axis,ok)
                    if(ok) then
                      Call Get_UB_from_uvw_hkl_omega(wave,Cell,Zone_Axis,h1,omeg,UB,ok,mess)
                      if(.not. ok) then
                        write(unit=*,fmt="(a)") " => UB cannot be calculated from ORIENTH instruction "
                        write(unit=*,fmt="(a)") " => "//trim(mess)
                        return
                      end if
                    else
                      write(unit=*,fmt="(a)") " => UB cannot be calculated from ORIENTH instruction "
                      write(unit=*,fmt="(a,6f7.3,a,f7.2)") " => Problem with Zone axis and reflection: ",uvw,h1, " with omega(h1)=",omeg
                      return
                    end if
                    ub_read=.true.
            end select
        end do

        if(igeom /= 0 .and. inst_read) Current_Instrm%igeom=igeom !overrides the instrument geometry

        inquire(file="ubfrom.raf",exist=esta)  !checking if "ubfrom.raf" is present in the current directory
                                               !If present, its content is prioritary w.r.t. the provided values
                                               !in the CFL file or instrument file.
        if(.not. ub_read  .or. esta) then      !Read UB-matrix from file ubfrom.raf
            if(.not. esta) then
               if(inst_read) then
                 ub=Current_Orient%ub
               else
                 mess="No UB-matrix in instrument/CFL files, no ubfrom.raf available!  "
                 return
               end if
            else  !Read everything from ubfrom.raf
               call get_logunit(i_def)
               open(unit=i_def, file="ubfrom.raf", status="old",action="read",position="rewind")
               do i=1,3
                 read(unit=i_def,fmt="(a)",iostat=ier) line
                 read(unit=line(1:50),fmt=*,iostat=ier) ub(i,:)
                 if(ier /= 0 ) then
                     mess="Error reading the orientation matrix in file: ubfrom.raf "
                     return
                 end if
               end do
               ub_read=.true.
               read(unit=i_def,fmt=*,iostat=ier) wav  !Read wavelength from file ubfrom.raf
               if(ier == 0 ) then
                 wave=wav
                 wave_read=.true.
               end if
               read(unit=i_def,fmt=*,iostat=ier) dcel !Read unit cell from file ubfrom.raf
               if(ier == 0 ) then
                 if(sum(abs(dcel-incel)) > 0.001) then  !Output information is the unit cell in ubfrom.raf
                   incel=dcel                           !is different from the values provided in CFL file
                   !Reconstruct the unit cell type
                   call Set_Crystal_Cell(incel(1:3),incel(4:6),Cell,"A")
                   write(unit=ipr,fmt="(/,/,a/)")  "  => New unit cell read from ubfrom.raf !"
                   call Write_Crystal_Cell(Cell,ipr)
                 end if
               end if
               if(trang_read) then !re-calculate smin,smax for the new wavelength is 2theta-range has been provided
                  smin=sind(0.5*tmin)/wave
                  smax=sind(0.5*tmax)/wave
               end if
               close(unit=i_def)
            end if
        end if

        if(.not. inst_read ) then
            file_inst="default_instrument.geom"
            call Set_default_Instrument()
            write(unit=*,fmt="(a)") " => No Instrument file has been provided!"
            write(unit=*,fmt="(a)") " => A default instrument file (default_instrm.geom) has been written"
            call get_logunit(i_def)
            open(unit=i_def, file=trim(file_inst), status="replace",action="write")
            call Write_Current_Instrm_data(i_def,trim(file_inst))
            call flush(i_def)
            close(unit=i_def)
        end if

        if(ub_read .and. wave_read) call Update_Current_Instrm_UB(trim(file_inst),UB,wave)

        !Check that the UB-matrix is consistent with the cell parameters
        call cell_fr_UB(ub,ipr,dcel)

        if( sum(abs(dcel-incel)) > 1.0) then
            mess="Incompatible UB-matrix and cell parameters in file "//trim(file_inst)
            return
        end if
        if(.not. Err_ILLdata) ok=.true.
    End Subroutine read_sxtal_geom

    Subroutine check_limits(ang,ok_ang)
        real, dimension(:),    intent (in) :: ang
        logical, dimension(:), intent(out) :: ok_ang
        integer :: i

        ok_ang=.false.
        do i=1,Current_Instrm%num_ang
            if(ang(i) >=  Current_Instrm%ang_Limits(i,1) .and.  ang(i) <=  Current_Instrm%ang_Limits(i,2)) ok_ang(i) =.true.
        end do
        return
    End subroutine check_limits

    Subroutine calc_angles(geom,sig,hr,ang,comment)
        integer,            intent (in)    :: geom
        integer,            intent (in out):: sig
        real, dimension(3), intent (in)    :: hr
        real, dimension(4), intent (out)   :: ang
        character(len=*),   intent (out)   :: comment
        !---- Local variables ----!
        real                  :: tteta,om,ch,ph
        integer               :: ierr
        logical, dimension(4) :: ok_ang
        logical               :: ok_om, ok_ch,ok_ph, ok_tteta
        equivalence   (ok_tteta,ok_ang(1)),( ok_om,ok_ang(2)), (ok_ch,ok_ang(3)), (ok_ph,ok_ang(4))

        comment=" "

        Select Case (geom)

            Case(1,2,4)
                call calang(hr,tteta,om,ch,ph,ierr)
                if(ierr /= 0) comment="Outside 2Theta limits"
                ang= (/tteta,om,ch,ph/)
                call check_limits(ang,ok_ang)
                if(.not. ok_om)    comment="Outside Omega limits"
                if(.not. ok_ch)    comment="Outside Chi limits"
                if(.not. ok_ph)    comment="Outside Phi limits"

            Case(3,-3)
                if(geom == 3) then
                    call Normal_Beam_angles(Current_Orient%wave,Current_Orient%ub,hr,sig,ang,ierr)
                else
                    call Normal_Beam_angles(Current_Orient%wave,Current_Orient%ub,hr,sig,ang,ierr,nusign=-1)
                end if
                if(ierr /= 0) then
                    if(ierr == 1) comment="Outside resolution sphere"
                    if(ierr == 2) comment="Blind zone for Nu"
                    if(ierr == 3) comment="Blind zone for Gamma"
                else
                    call check_limits(ang,ok_ang)
                    if(.not. ok_tteta) comment="Outside Gamma limits"
                    if(.not. ok_om)    comment="Outside Omega limits"
                    if(.not. ok_ch)    comment="Outside Nu limits"
                end if

        End Select

        return

    End subroutine calc_angles

End Module Ref_Gen

Program Sxtal_Ref_Gen

    use CFML_GlobalDeps,               only: cp,dp,pi
    use CFML_Math_general,             only: sort
    use CFML_Math_3D,                  only: cross_product
    use CFML_crystallographic_symmetry,only: space_group_type, Write_SpaceGroup, Set_SpaceGroup
    use CFML_Atom_TypeDef,             only: Atom_List_Type, Write_Atom_List,MAtom_list_Type
    use CFML_crystal_metrics,          only: Crystal_Cell_Type, Write_Crystal_Cell
    use CFML_Reflections_Utilities,    only: Reflection_List_Type, Hkl_Uni, Hkl_Gen_Sxtal, get_maxnumref, Hkl_Equiv_List
    use CFML_IO_Formats,               only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
    use CFML_Structure_Factors,        only: Structure_Factors, Write_Structure_Factors, &
                                             Init_Structure_Factors,Calc_General_StrFactor,&
                                             Scattering_Species_Type, Allocate_Scattering_Species, &
                                             Additional_Scattering_Factors, Set_Form_Factors
    use CFML_ILL_Instrm_data,          only: Err_ILLdata_Mess, Err_ILLdata, &
                                             SXTAL_Orient_type, Current_Orient, diffractometer_type, &
                                             Current_Instrm, Write_Current_Instrm_data
    use CFML_Geometry_SXTAL
    use CFML_Magnetic_Symmetry
    use CFML_Magnetic_Structure_Factors

    use Ref_Gen

    implicit none

    type (file_list_type)               :: fich_cfl
    type (space_group_type)             :: SpG,SpGT
    type (Atom_list_Type)               :: A
    type (Crystal_Cell_Type)            :: Cell
    type (Reflection_List_Type)         :: hkl

    type (MagSymm_k_Type)               :: MGp
    type (MAtom_list_Type)              :: Am
    type (MagH_Type)                    :: Mh
    type (MagH_List_Type)               :: Mhkl

    character(len=256)                  :: filcod     !Name of the input file
    character(len=256)                  :: mess       !Message after reading the CFL file for ref_gen
    character(len=15)                   :: sinthlamb  !String with stlmax (2nd cmdline argument)
    character(len=30)                   :: comment
    character(len=1)                    :: keyv
    real                                :: stlmin,stlmax     !Minimum and Maximum Sin(Theta)/Lambda
    real                                :: tteta, ss, dspc, SqMiV
    real,    dimension(3)               :: hr, z1, vk
    real,    dimension(3,3)             :: ub
    real                                :: sn,s2,theta,flip_right,flip_left,up,down,Nuc
    real                                :: start, tend
    integer, dimension(3)               :: h, ord
    integer, dimension(3,2)             :: hlim
    integer, dimension(3,48)            :: hlist
    real, dimension(4)                  :: ang
    real, dimension(:,:),allocatable    :: angles
    real, dimension(:,:),allocatable    :: reflx
    integer, dimension(:),  allocatable :: ind
    real,    dimension(:),  allocatable :: fst
    integer                             :: MaxNumRef, Num, lun=1, ier,i,j, ierr,i_hkl=2, n, iop
    integer                             :: narg, mul, sig, n_ini, n_end, nm, mu, nv
    logical                             :: esta, arggiven=.false.,sthlgiven=.false., ok=.false., schwinger=.false.,&
                                           spher=.false.,lim, iop_given=.false., mag_structure=.false.
    complex                             :: fn,fx,fe,fsru,fsrd,fslu,fsld
    Type(Scattering_Species_Type)       :: Scattf, add_Scat
    real(kind=dp), parameter            :: schw= -0.00014699



    !---- Arguments on the command line ----!
    narg=COMMAND_ARGUMENT_COUNT()
    stlmin=0.0

    if(narg > 0) then
        call GET_COMMAND_ARGUMENT(1,filcod)
        arggiven=.true.
        i=index(filcod,'.cfl',back=.true.)
        if( i /= 0) filcod=filcod(1:i-1)
    end if

    if(narg > 1) then
        call GET_COMMAND_ARGUMENT(2,sinthlamb)
        read(unit=sinthlamb,fmt=*,iostat=ier) stlmax
        if(ier == 0) sthlgiven=.true.
    end if

    if(narg > 2) then
        call GET_COMMAND_ARGUMENT(3,sinthlamb)
        read(unit=sinthlamb,fmt=*,iostat=ier) stlmin
        if(ier /= 0) stlmin=0.0
        if(stlmin > stlmax) then
            tteta=stlmin
            stlmin=stlmax
            stlmax=tteta
        end if
    end if

    if(narg > 3) then
        call GET_COMMAND_ARGUMENT(4,sinthlamb)
        read(unit=sinthlamb,fmt=*,iostat=ier) iop
        if(ier /= 0) iop=0
        iop_given=.true.
    end if

    write(unit=*,fmt="(/,/,7(a,/))")                                                  &
          "                      ------ PROGRAM SXTAL_REFGEN ------"                 , &
          "                      ---- Version  0.3 April-2013  ----"                 , &
          "    *******************************************************************"  , &
          "    * Generates single crystal reflections and nuc. structure factors *"  , &
          "    * (if atoms are given ) reading a *.CFL file (4C & NB geometries) *"  , &
          "    *******************************************************************"  , &
          "                     (JRC- April-2007, updated April 2013 )"
    write(unit=*,fmt=*) " "

    if(.not. arggiven) then
        write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.cfl (give xx): "
        read(unit=*,fmt="(a)") filcod
        if(len_trim(filcod) == 0) stop
    end if

    open(unit=lun,file=trim(filcod)//".sfa", status="replace",action="write")

    write(unit=lun,fmt="(/,/,7(a,/))")                                                &
          "                      ------ PROGRAM SXTAL_REFGEN ------"                 , &
          "                      ---- Version  0.3 April-2013  ----"                 , &
          "    *******************************************************************"  , &
          "    * Generates single crystal reflections and nuc. structure factors *"  , &
          "    * (if atoms are given ) reading a *.CFL file (4C & NB geometries) *"  , &
          "    *******************************************************************"  , &
          "                     (JRC- April-2007, updated April 2013 )"

    inquire(file=trim(filcod)//".cfl",exist=esta)
    if( .not. esta) then
        write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cfl doesn't exist!"
        stop
    end if

    call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)

    call cpu_time(start)

    If(err_form) then

        write(unit=*,fmt="(a)") "  "//trim(err_form_mess)

    else

        call Write_Crystal_Cell(Cell,lun)
        call Write_SpaceGroup(SpG,lun)
        if(A%natoms > 0) call Write_Atom_List(A,lun=lun)


        n_ini=1
        n_end=fich_cfl%nlines
        call Readn_Set_Magnetic_Structure(fich_cfl,n_ini,n_end,MGp,Am)
        if(err_MagSym) then
            write(unit=*,fmt="(a)") "   "//err_MagSym_Mess
            mag_structure=.false.
        else
            mag_structure=.true.
            !if(Am%natoms > 0) then
            call Write_Magnetic_Structure(lun,MGp,Am)
            !end if
        end if

        !re-read the input file searching for INSTRM, UBM, GEOM, WAVE, ORDER, SPHER items
        if(.not. sthlgiven) then
            call read_sxtal_geom(lun,cell,fich_cfl,ok,mess,ord,spher,schwinger,hlim, &
                                 sgiven=sthlgiven,smin=stlmin,smax=stlmax,opgiven=iop_given,iop_ord=iop)
        else
            if(.not. iop_given) then
                call read_sxtal_geom(lun,cell,fich_cfl,ok,mess,ord,spher,schwinger,hlim,  &
                                     opgiven=iop_given,iop_ord=iop)
            else
                call read_sxtal_geom(lun,cell,fich_cfl,ok,mess,ord,spher,schwinger,hlim)
            end if
        end if

        if(.not. ok) then
            write(unit=*,fmt="(a)") "     Mess: "//trim(mess)
            if(Err_ILLdata) write(unit=*,fmt="(a)") " ILL_data: "//trim(Err_ILLdata_Mess)
            stop
        end if

        if(.not. sthlgiven) then
            write(unit=*,fmt="(a)",advance="no") " => Minimum and Maximum sinTheta/Lambda: "
            read(unit=*,fmt=*) stlmin,stlmax
        end if
        MaxNumRef = get_maxnumref(stlmax,Cell%CellVol,mult=SpG%Multip)
        if(spher) MaxNumRef=MaxNumRef*Spg%NumOps*max(Spg%Centred,1)
        !write(unit=*,fmt="(a,i6,a)") " Generate at least ",MaxNumRef," reflections"
        lim=.false.
        if(sum(abs(hlim)) /= 0) lim=.true.
        if(spher) then
            if(lim) then
                call Hkl_gen_sxtal(cell,SpG,stlmin,stlmax,MaxNumRef,hkl,ord,hlim)
            else
                call Hkl_gen_sxtal(cell,SpG,stlmin,stlmax,MaxNumRef,hkl,ord)
            end if
        else
            call Hkl_Uni(cell,SpG,.true.,stlmin,stlmax,"s",MaxNumRef,hkl, no_order=.true.)
        end if

        !Calculation of Structure factors for neutron scattering
        if(A%natoms /= 0) then
            call Init_Structure_Factors(hkl,A,Spg,mode="NUC",lun=lun)
            call Structure_Factors(A,SpG,hkl,mode="NUC")
        end if

        !Allocation of memory and initialization of arrays
        if(spher) then
            n=hkl%nref
        else
            n=hkl%nref*Spg%NumOps*max(Spg%Centred,1)
        end if
        if(allocated(angles)) deallocate (angles)
        allocate(angles(4,n))
        angles=0.0
        if(allocated(reflx)) deallocate (reflx)
        allocate(reflx(3,n))
        reflx=0.0
        if(allocated(fst)) deallocate (fst)
        allocate(fst(n))
        fst=0.0
        if(allocated(ind)) deallocate(ind)
        allocate(ind(n))

        ind=(/(i,i=1,n)/)

        sig=1
        if(Current_instrm%BL_frame == "z-down") sig=-1

        open(unit=i_hkl, file=trim(filcod)//".hkl", status="replace",action="write")
        call Write_Current_Instrm_data(lun)

        if(.not. iop_given) then
            write(unit=*,fmt="(a)") " => Order of reflections: "
            write(unit=*,fmt="(a)") "              0: As set on CFL file (hkl-running) "
            write(unit=*,fmt="(a)") "              1: Ascending order on 2Theta "
            write(unit=*,fmt="(a)") "              2: Ascending order on Omega "
            write(unit=*,fmt="(a)") "              3: Ascending order on Chi/Nu "
            write(unit=*,fmt="(a)") "              4: Ascending order on  Phi "
            write(unit=*,fmt="(a)",advance="no") " => Enter option: "
            read(unit=*,fmt=*,iostat=ier) iop
            if(ier /= 0) iop=0
        end if
        n=0

        write(unit=lun,fmt="(/,a,2f8.5,a/)")" => Reflection generation between sinTheta/Lambda min and max: ", &
                                                 stlmin,stlmax, " angstroms^-1"

        if(Current_instrm%igeom == 3 .or. Current_instrm%igeom == -3) then
            write(unit=lun,fmt="(/,/,a)")"----------------------------------------------------------------------------"
            write(unit=lun,fmt="(a)")    "   h   k   l       F2        Re(F)       Im(F)     Gamma     Omega      Nu  "
            write(unit=lun,fmt="(a)")    "----------------------------------------------------------------------------"
        else
            write(unit=lun,fmt="(/,/,a)")"---------------------------------------------------------------------------------------"
            write(unit=lun,fmt="(    a)")"   h   k   l       F2        Re(F)       Im(F)     2Theta     Omega     Chi       Phi  "
            write(unit=lun,fmt="(    a)")"---------------------------------------------------------------------------------------"
        end if

        if(spher) then
            do i=1,hkl%nref
                hr=hkl%Ref(i)%h
                call calc_angles(Current_Instrm%igeom,sig,hr,ang,comment)

                Select Case (Current_Instrm%igeom)

                    Case(1,2,4)
                        write(unit=lun,fmt="(3i4,3f12.5,4f10.3,tr2,a)") hkl%Ref(i)%h,hkl%ref(i)%Fc,hkl%ref(i)%A,hkl%ref(i)%B,&
                                                                        ang,trim(comment)

                    Case(3,-3)
                        write(unit=lun,fmt="(3i4,3f12.5,3f10.3,tr2,a)") hkl%Ref(i)%h,hkl%ref(i)%Fc,hkl%ref(i)%A,hkl%ref(i)%B,&
                                                                        ang(1:3),trim(comment)

                End Select

                if(len_trim(comment) == 0) then
                    n=n+1
                    angles(:,n)  = ang(:)
                    reflx(:,n)   = hkl%Ref(i)%h
                    fst(n)       = hkl%ref(i)%Fc
                end if

            end do

        else

            do i=1,hkl%nref
                h=hkl%Ref(i)%h
                call Hkl_Equiv_List(h,SpG,.true.,Mul,Hlist)

                do j=1,mul
                    hr=Hlist(:,j)
                    !z1=matmul(Current_Orient%UB,hr)
                    call calc_angles(Current_Instrm%igeom,sig,hr,ang,comment)

                    Select Case (Current_Instrm%igeom)
                        Case(1,2,4)
                            write(unit=lun,fmt="(3i4,3f12.5,4f10.3,tr2,a)") hlist(:,j),hkl%ref(i)%Fc,hkl%ref(i)%A,hkl%ref(i)%B,&
                                                                            ang,trim(comment)
                        Case(3,-3)
                            write(unit=lun,fmt="(3i4,3f12.5,3f10.3,tr2,a)") hkl%Ref(i)%h,hkl%ref(i)%Fc,hkl%ref(i)%A,hkl%ref(i)%B,&
                                                                            ang(1:3),trim(comment)
                    End Select

                    if(len_trim(comment) == 0) then
                        n=n+1
                        angles(:,n)  = ang
                        reflx(:,n)   = hlist(:,j)
                        fst(n)       = hkl%ref(i)%Fc
                    end if

                end do

            end do
        end if
        !write final output file with a selected ordering according to a motor
        if(iop /= 0) then
            call sort(angles(iop,:),n,ind)
        end if
        do i=1,n
            j=ind(i)
            write(unit=i_hkl,fmt="(3i4,f12.5,4f10.3)") nint(reflx(:,j)),fst(j), angles(:,j)
        end do
        close(unit=i_hkl)

        ! Generation of magnetic satellites if needed
        nm=0
        if(mag_structure) then
            !First generate the full set of reflections according to the lattice of the
            !space group.
            comment=SpG%SPG_lat//" -1"
            call Set_SpaceGroup(comment,SpGT)
            if(lim) then
                call Hkl_gen_sxtal(cell,SpGT,stlmin,stlmax,MaxNumRef,hkl,ord,hlim)
            else
                call Hkl_gen_sxtal(cell,SpGT,stlmin,stlmax,MaxNumRef,hkl,ord)
            end if
            open(unit=i_hkl, file=trim(filcod)//".mhkl", status="replace",action="write")
            call Gen_Satellites(Cell,MGp,stlmax,Mhkl,hkl=hkl)  !Mhkl magnetic satellites
            write(unit=lun,fmt="(/,/,a)") " "

            if(Am%natoms > 0) then
                call Init_Mag_Structure_Factors(Mhkl,Am,MGp,lun)
                if(err_msfac) then
                    write(unit=*,fmt="(a)")  "  "//err_msfac_mess
                    stop
                end if
                call Mag_Structure_Factors(Cell,Am,MGp,Mhkl)
                call Calc_Mag_Interaction_Vector(Mhkl,Cell)  !in {e1,e2,e3} basis (complete Mhkl)
            end if

            nm=Mhkl%nref
            if(allocated(angles)) deallocate (angles)
            allocate(angles(4,nm))
            angles=0.0
            if(allocated(reflx)) deallocate (reflx)
            allocate(reflx(3,nm))
            reflx=0.0
            if(allocated(fst)) deallocate (fst)
            allocate(fst(nm))
            fst=0.0
            if(allocated(ind)) deallocate(ind)
            allocate(ind(nm))

            ind=(/(i,i=1,nm)/)

            write(unit=lun,fmt="(/,a)")   "    LIST OF REFLECTIONS AND MAGNETIC STRUCTURE FACTORS"
            write(unit=lun,fmt="(a,/)")   "    =================================================="
            write(unit=lun,fmt="(a,i6)") " => Total number of reflections  : ",Mhkl%Nref
            write(unit=lun,fmt="(a,i6)") " => Number of propagation vectors: ",MGp%nkv
            do i=1,MGp%nkv
                write(unit=lun,fmt="(a,i2,a,3f8.4,a)") " => Propagation vectors #",i," = (",MGp%kvec(:,i)," )"
            end do

            if(Current_instrm%igeom == 3 .or. Current_instrm%igeom == -3) then
                write(unit=lun,fmt="(/,/,a)")"-------------------------------------------------------------------------"// &
                                             "-------------------------------------------------------------------------"//&
                                             "-------------------------------------------------------------------------"
                write(unit=lun,fmt="(a)") &
                "    Hr      Kr      Lr       H   K   L   nvk   Mult   dspc      |MiV|^2     Mrx      Mry      Mrz      "// &
                "Mix      Miy      Miz     MiVrx    MiVry    MiVrz    MiVix    MiViy    MiViz    Gamma     Omega      Nu  "

                write(unit=lun,fmt="(a)")"-------------------------------------------------------------------------"// &
                                         "-------------------------------------------------------------------------"//&
                                         "-------------------------------------------------------------------------"
            else
                write(unit=lun,fmt="(/,/,a)")"-------------------------------------------------------------------------"// &
                                             "-------------------------------------------------------------------------"//&
                                            "--------------------------------------------------------------------------------------"
                write(unit=lun,fmt="(a)") &
                "    Hr      Kr      Lr       H   K   L   nvk   Mult   dspc      |MiV|^2     Mrx      Mry      Mrz      "// &
                "Mix      Miy      Miz     MiVrx    MiVry    MiVrz    MiVix    MiViy    MiViz    2Theta     Omega     Chi       Phi"
                write(unit=lun,fmt="(    a)")"-------------------------------------------------------------------------"// &
                                             "-------------------------------------------------------------------------"//&
                                            "--------------------------------------------------------------------------------------"
            end if

            nm=0
            do i=1,Mhkl%nref
                hr=Mhkl%Mh(i)%h
                ss=Mhkl%Mh(i)%signp
                mul=Mhkl%Mh(i)%mult
                nv =Mhkl%Mh(i)%Num_k
                vk=MGp%kvec(:,nv)
                h=nint(hr+sig*vk)
                sqMiV= Mhkl%Mh(i)%sqMiV
                dspc=0.5/Mhkl%Mh(i)%S

                call calc_angles(Current_Instrm%igeom,sig,hr,ang,comment)

                Select Case (Current_Instrm%igeom)
                    Case(1,2,4)
                        write(unit=lun,fmt="(3f8.3,tr2,3i4,i5,i6,f9.4,f13.5,12f9.4,4f10.3,tr2,a)") hr,h, -sig*nv, mul, dspc,sqMiV, &
                        real(Mhkl%Mh(i)%MsF),aimag(Mhkl%Mh(i)%MsF), real(Mhkl%Mh(i)%MiV),aimag(Mhkl%Mh(i)%MiV), &
                        ang,trim(comment)
                    Case(3,-3)
                        write(unit=lun,fmt="(3f8.3,tr2,3i4,i5,i6,f9.4,f13.5,12f9.4,3f10.3,tr2,a)") hr,h, -sig*nv, mul, dspc,sqMiV, &
                        real(Mhkl%Mh(i)%MsF),aimag(Mhkl%Mh(i)%MsF), real(Mhkl%Mh(i)%MiV),aimag(Mhkl%Mh(i)%MiV), &
                        ang(1:3),trim(comment)
                End Select

                if(len_trim(comment) == 0) then
                    nm=nm+1
                    angles(:,nm)  = ang
                    reflx(:,nm)   = hr
                    fst(nm)       = Mhkl%Mh(i)%sqMiV
                end if
            end do
            !write final output file with a selected ordering according to a motor
            if(iop /= 0) then
                call sort(angles(iop,:),nm,ind)
            end if
            do i=1,nm
                j=ind(i)
                write(unit=i_hkl,fmt="(3f9.4,f12.5,4f10.3)") reflx(:,j),fst(j), angles(:,j)
            end do
        end if

        if(Schwinger .and. A%natoms /= 0) then

            call Additional_Scattering_Factors(fich_cfl,add_Scat,ok,mess)
            if(.not. ok) then
              write(unit=*,fmt="(a)") " => "//trim(mess)
              stop
            end if
            !  Set nuclear, x-ray, neutron and magnetic form factors coefficients for all the atoms
            !  in the structure
            if(add_Scat%Num_Species > 0) then
              call Set_Form_Factors(A,Scattf,ok,mess,Current_Orient%wave,lun,add_Scat)
            else
              call Set_Form_Factors(A,Scattf,ok,mess,Current_Orient%wave,lun)
            end if
            if(.not. ok) then
              write(unit=*,fmt="(a)") " => "//trim(mess)
              stop
            end if
           write(unit=lun,fmt="(/,a)") &
           "   H   K   L   sinT/L  Intensity     FNr       FNi           FXr       FXi           FEr       FEi      Flip_left  Flip_right"// &
                                    "  SRUr      SRUi          SRDr      SRDi          SLUr      SLUi          SLDr      SLDi"
           n=0
           angles=0.0
           ind=(/(i,i=1,hkl%NRef)/)
           do i=1,hkl%NRef
              hr=real(hkl%ref(i)%h)
              sn = hkl%ref(i)%s
              theta=Current_Orient%wave*sn
              theta=asin(theta)
              s2=sn*sn
              call Calc_General_StrFactor(hr,s2,A,SpG,Scattf,fn,fx,fe)
              Nuc=fn*conjg(fn)
              fsru=Schwinger_Amplitude(hr,[0.0,0.0, 1.0],theta,fe,Left=.false.)
              fslu=Schwinger_Amplitude(hr,[0.0,0.0, 1.0],theta,fe,Left=.true.)
              fsrd=Schwinger_Amplitude(hr,[0.0,0.0,-1.0],theta,fe,Left=.false.)
              fsld=Schwinger_Amplitude(hr,[0.0,0.0,-1.0],theta,fe,Left=.true.)
              up= (fn+fsru)*conjg(fn+fsru); down= (fn+fsrd)*conjg(fn+fsrd)
              if(up < 1.0e-6 .and. down < 1.0e-6) then
                 flip_right=1.0
              else
                 flip_right=up/down
              end if
              up=(fn+fslu)*conjg(fn+fslu); down= (fn+fsld)*conjg(fn+fsld)
              if(up < 1.0e-6 .and. down < 1.0e-6) then
                 flip_left=1.0
              else
                 flip_left=up/down
              end if

              call calc_angles(Current_Instrm%igeom,sig,hr,ang,comment)

              Select Case (Current_Instrm%igeom)

                  Case(1,2,4)
                      write(unit=lun,fmt="(3i4,f9.5,f10.4,3(2f10.4,tr4),2f10.5,4(2f10.6,tr4),4f10.3,tr2,a)") &
                            hkl%Ref(i)%h, hkl%Ref(i)%s,Nuc,fn,fx,fe,flip_left,flip_right,&
                            fsru,fsrd,fslu,fsld,ang,trim(comment)

                  Case(3,-3)
                      write(unit=lun,fmt="(3i4,f9.5,f10.4,3(2f10.4,tr4),2f10.5,4(2f10.6,tr4),3f10.3,tr2,a)") &
                            hkl%Ref(i)%h, hkl%Ref(i)%s,Nuc,fn,fx,fe,flip_left,flip_right,&
                            fsru,fsrd,fslu,fsld,ang(1:3),trim(comment)
              End Select

              if(len_trim(comment) == 0) then
                  n=n+1
                  angles(1:4,n)  = (/flip_right,flip_right,flip_right,flip_right/)
                  reflx(:,n)   = hr
                  fst(n)       = Nuc
              end if

           end do
            !write final output file with a selected ordering according maximum right flipping ratio

            call sort(angles(1,:),n,ind)

            do i=n,1,-1
                j=ind(i)
                write(unit=i_hkl,fmt="(3f9.4,f12.5,4f10.3)") reflx(:,j),fst(j), angles(:,j)
            end do

        end if

        write(unit=*,fmt="(a)")    " Normal End of: PROGRAM SXTAL_REFGEN "
        write(unit=*,fmt="(a,i5)") " Number of accessible nuclear  reflections: ",n
        if(mag_structure) write(unit=*,fmt="(a,i5)") " Number of accessible magnetic reflections: ",nm
        write(unit=*,fmt="(a)")                      " Results in Files: "//trim(filcod)//".sfa  and "//trim(filcod)//".hkl"
        if(mag_structure) write(unit=*,fmt="(a)")    " Magnetic Satellites in: "//trim(filcod)//".mhkl"
    end if

    close(unit=lun)

    call cpu_time(tend)

    write(unit=*,fmt="(/,a,f10.2,a)")  "  CPU-Time: ", tend-start," seconds"
    write(unit=*,fmt="(/,a)") " => Press <enter> to finish "
    read(unit=*,fmt="(a)") keyv

    stop

  contains

    Function Schwinger_Amplitude(hn,pol,theta,fe,Left,UB)  result(Schwinger)
      real(kind=cp), dimension(3),            intent(in) :: hn,pol
      real(kind=cp),                          intent(in) :: theta !in radians
      complex(kind=cp),                       intent(in) :: fe
      logical,                      optional, intent(in) :: left
      real(kind=cp), dimension(3,3),optional, intent(in) :: UB
      complex(kind=cp)                                   :: Schwinger

      real(kind=cp),dimension(3) :: uvect, hc

      if(present(left)) then
        if(left) then
          uvect=(/0.0,0.0,-1.0/)
        else
          uvect=(/0.0,0.0,1.0/)
        end if
      else if(present(UB)) then
        hc=matmul(UB,hn)
        uvect=cross_product((/0.0,1.0,0.0/),hc)
        uvect=-uvect/sqrt(dot_product(uvect,uvect))
      else
        Schwinger=0.0
        return
      end if
      !change of sign of "i" because the convention used by crystallographers is Q=Kf-Ki
      Schwinger= schw*cmplx(0.0,-1.0)*dot_product(pol,uvect)*fe/abs(tan(theta))
    End Function Schwinger_Amplitude



End Program Sxtal_Ref_Gen
