Program Optimizing_structures
   !---- Use Modules ----!
   use crystallographic_symmetry, only: space_group_type, Write_SpaceGroup
   use string_utilities,          only: u_case
   use Atom_Module,               only: Atom_List_Type, Write_Atom_List
   use crystal_types,             only: Crystal_Cell_Type, Write_Crystal_Cell
   use Reflections_Utilities,     only: Reflection_List_Type
   use IO_Formats,                only: Readn_set_Xtal_Structure,err_mess_form,err_form,&
                                        file_list_type, File2File_List
   use Configuration_Calculations,only: calc_BVS
   use Structure_Factor_Module,   only: Structure_Factors, Write_Structure_Factors, &
                                        Init_Structure_Factors, err_sfac,err_mess_sfac
   use Refinement_Codes,          only: NP_Max,NP_Refi,V_BCon,V_Bounds,V_Name,V_Vec,&
                                        Allocate_VParam,Init_RefCodes, Read_RefCodes_File, &
                                        Write_Info_RefCodes, Err_RefCodes, Err_Mess_RefCodes, &
                                        allocate_restparam, Write_restraints_ObsCalc
   use Optimization_SAN,          only: SimAnn_Conditions_type, state_Vector_Type, Multistate_Vector_Type, &
                                        err_mess_SAN,err_SAN, Simanneal_Gen,Set_SimAnn_Cond, &
                                        Set_SimAnn_StateV,Write_SimAnn_Cond, Write_SimAnn_StateV, &
                                        Write_SimAnn_MStateV, SimAnneal_MultiConf,Set_SimAnn_MStateV, &
                                        Sann_Opt_MultiConf
   use observed_reflections,      only: Read_observations,Observation_Type,Observation_List_Type, &
                                        Write_ObsCalc_SFactors,err_observ,err_mess_observ, SumGobs,&
                                        wavel_int, Write_FoFc_Powder
   use cost_functions,            only: Cell,A,Ac,SpG,hkl,Oh,Icost,wcost,Err_cost,Err_Mess_cost, &
                                        General_Cost_function, readn_set_costfunctpars, write_costfunctpars, &
                                        Write_FinalCost,wavel,diff_mode
   !---- Local Variables ----!
   implicit none
   type (file_list_type)              :: fich_cfl,fich_rest
   type (SimAnn_Conditions_type),save :: c
   type (state_Vector_Type)           :: vs
   type (Multistate_Vector_Type)      :: mvs


   character(len=256)                 :: filcod     !Name of the input file
   character(len=256)                 :: filhkl     !Name of the hkl-file
   character(len=256)                 :: filrest    !Name of restraints file
   character(len=256)                 :: line       !Text line
   character(len=256)                 :: fst_cmd    !Commands for FP_Studio
   integer                            :: MaxNumRef, Num, lun=1, ier,i,j,k, i_hkl=2, n
   real                               :: start,fin, mindspc, maxsintl
   integer                            :: narg,iargc
   Logical                            :: esta, arggiven=.false.,sthlgiven=.false., &
                                         fst_out=.false., local_opt=.false., rest_file=.false.

    !---- Arguments on the command line ----!
    narg=iargc()

    if(narg > 0) then
            call getarg(1,filcod)
            arggiven=.true.
    end if

    write (unit=*,fmt="(/,/,6(a,/))")                                                  &
            "            ------ PROGRAM FOR OPTIMIZING X-TAL STRUCTURES ------"            , &
            "                    ---- Version 0.2 April-2007----"                          , &
            "    *******************************************************************************"  , &
            "    * Optimizes a X-tal structure reading integrated intensities and a *.CFL file *"  , &
            "    *******************************************************************************"  , &
            "                             (JRC- April 2007 )"
   write (unit=*,fmt=*) " "

   if(.not. arggiven) then
     write(unit=*,fmt="(a)") " => Code of the file xx.cfl (give xx): "
     read(unit=*,fmt="(a)") filcod
     if(len_trim(filcod) == 0) stop
   end if

   open(unit=lun,file=trim(filcod)//".out", status="replace",action="write")
    write(unit=lun,fmt="(/,/,6(a,/))")                                                  &
            "            ------ PROGRAM FOR OPTIMIZING X-TAL STRUCTURES ------"            , &
            "                    ---- Version 0.2 April-2007----"                          , &
            "    *******************************************************************************"  , &
            "    * Optimizes a X-tal structure reading integrated intensities and a *.CFL file *"  , &
            "    *******************************************************************************"  , &
            "                             (JRC- April 2007 )"

   inquire(file=trim(filcod)//".cfl",exist=esta)
   if( .not. esta) then
     write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cif (or .cfl) does'nt exist!"
     stop
   end if

   call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)

   If(err_form) then
     write(unit=*,fmt="(a)") trim(err_mess_form)
   else
     call Write_Crystal_Cell(Cell,lun)
     call Write_SpaceGroup(SpG,lun)
     call Write_Atom_List(A,lun=lun)

     np_max= A%natoms*11  ! x,y,z,biso,occ,betas
     call allocate_vparam(np_max)

    !Read what cost functions will be minimized
     Call Readn_Set_CostFunctPars(fich_cfl)
     if(Err_cost) then
       write(unit=*,fmt="(a)") trim(Err_Mess_cost)
       stop
     end if
     Call Write_CostFunctPars(lun)

     !Allocate objects for restraints
     if(icost(2) == 1 .or. icost(3) == 1 .or. icost(4) == 1) then
       !Look first for restraints file
          do i=1,fich_cfl%nlines
            line=adjustl(u_case(fich_cfl%line(i)))
            if(line(1:10) == "RESTR_FILE") then
              j=index(line,"!")
              if(j /= 0) line=trim(line(1:j-1))
              filrest= adjustl(line(12:))
              inquire(file=trim(filrest),exist=esta)
              call File2File_List(filrest,fich_rest)
              If(err_form) then
                 write(unit=*,fmt="(a)") trim(err_mess_form)
                 stop
              else
                 rest_file=.true.
              end if
              exit
            end if
          end do
        if(rest_file) then
           call allocate_restparam(fich_rest)
        else
           call allocate_restparam(fich_cfl)
        end if
     end if

     !Determine the refinement codes from vary/fix instructions
     call Init_RefCodes(A)
     call Read_RefCodes_File(fich_cfl,1,fich_cfl%nlines,A,Spg)
     if(Err_RefCodes) then
       write(unit=*,fmt="(a)") trim(Err_Mess_RefCodes)
     end if
     if(rest_file) then
        call Read_RefCodes_File(fich_rest,1,fich_rest%nlines,A,Spg)
        if(Err_RefCodes) then
          write(unit=*,fmt="(a)") trim(Err_Mess_RefCodes)
        end if
     end if
     call Write_Info_RefCodes(A,Spg,lun)

     !--- Look for FP_Studio commands
     fst_cmd=" "
     do i=1,fich_cfl%nlines
       line=adjustl(u_case(fich_cfl%line(i)))
       if(line(1:7) == "FST_CMD") then
         j=index(line,"!")
         if(j /= 0) line=trim(line(1:j-1))
         fst_cmd= adjustl(line(9:))
         fst_out=.true.
         exit
       end if
     end do

     !--- Look for Local_Opt commands
     do i=1,fich_cfl%nlines
       line=adjustl(u_case(fich_cfl%line(i)))
       if(line(1:9) == "LOCAL_OPT") then
         Local_Opt=.true.
         exit
       end if
     end do

     if(Icost(1) == 1 .or. Icost(7) == 1) then
          ! Reading observed structure factors squared and construct hkl if Icost(1) =1
          ! First detect the name of the hkl file
          esta=.false.
          do i=1,fich_cfl%nlines
            line=adjustl(u_case(fich_cfl%line(i)))
            if(line(1:7) == "HKL-OBS") then
              j=index(line,"!")
              if(j /= 0) line=trim(line(1:j-1))
              filhkl= adjustl(line(9:))
              inquire(file=trim(filhkl),exist=esta)
              exit
            end if
          end do

          if(.not. esta) then
            write(unit=*,fmt="(a)") " => No hkl-file (or wrong name!) has been provided! "
            Icost(1)=0; wcost(1)=0.0
            Icost(7)=0; wcost(7)=0.0
          end if

          mindspc=0.0       !Read minimun d-spacing to be considered in the list
          do i=1,fich_cfl%nlines
            line=u_case(fich_cfl%line(i))
            if(line(1:12) == "MIN-DSPACING") then
              read(unit=line(14:),fmt=*,iostat=ier) mindspc
              if(ier /= 0) mindspc=0.0
              exit
            end if
          end do
          if(Icost(7) == 1) then
             call read_observations(trim(filhkl),Cell,SpG,.true.,Oh,hkl)
          else
             call read_observations(trim(filhkl),Cell,SpG,.true.,hkl)
          end if
          if(err_observ) then
             write(unit=*,fmt="(a)") " => Error reading observations"
             write(unit=*,fmt="(a)") trim(err_mess_observ)
            stop
          end if

          if(Icost(7) == 1) then
            wavel=wavel_int
          end if

        ! Change the number of reflections to be considered according to the
        ! resolution (mindspc) read above. Change also SumGobs
         if( mindspc > 0.001) then
            n=0
            maxsintl=0.5/mindspc
            maxsintl=maxsintl*1.01

            do i=1,hkl%Nref
              n=n+1
              if(hkl%Ref(i)%s > maxsintl) then
                hkl%Nref=n
                exit
              end if
            end do

            if(Icost(7) == 1) then
              do i=1,Oh%Nobs
                  j=Oh%Ob(i)%ncont
                  k=Oh%Ob(i)%p(j)
                  if(k > n) then
                    Oh%Nobs=i-1
                    exit
                  end if
              end do
            end if

            SumGobs=0.0
            if(Icost(7) == 1) then
              do i=1,Oh%nobs
                SumGobs=SumGobs+Oh%Ob(i)%Gobs
              end do
            else
              do i=1,n
                SumGobs=SumGobs+abs(hkl%Ref(i)%Fo)
              end do
            end if

            write(unit= * ,fmt="(/,a,i5,/,a,f8.4)") " => Number of reflections to consider: ",n,&
                                                    "    Resolution: ",mindspc
            write(unit=lun,fmt="(/,a,i5,/,a,f8.4)") " => Number of reflections to consider: ",n,&
                                                    "    Resolution: ",mindspc
            if(Icost(7) == 1) then
              write(unit= * ,fmt="(a,i5)") " => Number of Observations to consider: ",Oh%nobs
              write(unit=lun,fmt="(a,i5)") " => Number of Observations to consider: ",Oh%nobs
            end if
         end if

        !Set up Structure factor tables and initial values
         call Init_Structure_Factors(hkl,A,Spg,mode=diff_mode,lambda=wavel,lun=lun)
         call Structure_Factors(A,SpG,hkl,mode=diff_mode,lambda=wavel)
          if(err_sfac) then
             write(unit=*,fmt="(a)") " => Error in Structure factor calculations"
             write(unit=*,fmt="(a)") trim(err_mess_sfac)
            stop
          end if
     End if !Icost(1) == 1

    !Set up the simulated annealing conditions
     call Set_SimAnn_Cond(fich_cfl,c)
      if(err_SAN) then
         write(unit=*,fmt="(a)") " => Error setting Simulated Annealing conditions"
         write(unit=*,fmt="(a)") trim(err_mess_SAN)
        stop
      end if
     call Write_SimAnn_Cond(lun,c)

     call cpu_time(start)

     if(c%num_conf == 1) then
        !Set up the Simulated Annealing State Vector
         call Set_SimAnn_StateV(NP_Refi,V_BCon,V_Bounds,V_Name,V_Vec,vs)
          if(err_SAN) then
             write(unit=*,fmt="(a)") " => Error setting Simulated Annealing State Vector"
             write(unit=*,fmt="(a)") trim(err_mess_SAN)
            stop
          end if
         call Write_SimAnn_StateV(lun,vs,"INITIAL STATE")
         if(fst_out) then
           call Simanneal_Gen(General_Cost_function,c,vs,lun,fst="simann.fst  "//trim(fst_cmd))
         else
           call Simanneal_Gen(General_Cost_function,c,vs,lun)
         end if
         call Write_SimAnn_StateV(lun,vs,"FINAL STATE")

     else   ! Multi-configuration Simulated Annealing

         call Set_SimAnn_MStateV(NP_Refi,c%num_conf,V_BCon,V_Bounds,V_Name,V_Vec,mvs)
          if(err_SAN) then
             write(unit=*,fmt="(a)") " => Error setting Simulated Annealing State Vector"
             write(unit=*,fmt="(a)") trim(err_mess_SAN)
            stop
          end if
         call Write_SimAnn_MStateV(lun,mvs,"INITIAL STATE")

         if (Local_Opt) then
           if(fst_out) then
              call SAnn_Opt_MultiConf(General_Cost_function,c,mvs,lun,fst="simann.fst  "//trim(fst_cmd))
           else
              call SAnn_Opt_MultiConf(General_Cost_function,c,mvs,lun)
           end if
         else
           if(fst_out) then
              call SimAnneal_MultiConf(General_Cost_function,c,mvs,lun,fst="simann.fst  "//trim(fst_cmd))
           else
              call SimAnneal_MultiConf(General_Cost_function,c,mvs,lun)
           end if
         end if

         call Write_SimAnn_MStateV(lun,mvs,"FINAL STATE",1)
     end if
     call Write_FinalCost(lun)
     call Write_Atom_List(A,lun=lun)

     if(Icost(1) == 1) call Write_ObsCalc_SFactors(lun,hkl,mode=diff_mode)
     if(Icost(7) == 1) call Write_FoFc_Powder(lun,Oh,hkl,mode=diff_mode)
     if(icost(2) == 1 .or. icost(3) == 1 .or. icost(4) == 1) then
        call Write_restraints_ObsCalc(A,lun)
     end if
     if(icost(5) == 1 .or. icost(6) == 1) then
       call Calc_BVS(Ac,lun)
     end if
     write(unit=*,fmt="(a)") " Normal End of: PROGRAM FOR OPTIMIZING X-TAL STRUCTURES "
     write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".out"
     call cpu_time(fin)
     write(unit=*,fmt="(a,f10.2,a)")  "  CPU-Time: ", fin-start," seconds"
     write(unit=*,fmt="(a,f10.2,a)")  "  CPU-Time: ", (fin-start)/60.0," minutes"
     write(unit=lun,fmt="(/,a,f10.2,a)")  "  CPU-Time: ", fin-start," seconds"
   end if

   close(unit=lun)
   stop
End Program Optimizing_structures

    Subroutine Write_FST(fst_file,v,cost)
      use String_Utilities,  only:  get_logunit
      Use Refinement_Codes,  only:  VState_to_AtomsPar
      use cost_functions,    only:  Cell,A,SpG

       character(len=*),     intent(in):: fst_file
       real,dimension(:),    intent(in):: v
       real,                 intent(in):: cost
       !----- Local variables -----!
       integer :: lun,i,nc, ier
       character(len=132)                 :: file_fst,fst_cmd
       character(len=30), dimension(10)   :: cmds

       i=index(fst_file,".fst")
       file_fst=fst_file(1:i+3)
       fst_cmd=adjustl(fst_file(i+4:))
       nc=0

       do
         i=index(fst_cmd,";")
         if(i /= 0) then
           nc=nc+1
           cmds(nc)=fst_cmd(1:i-1)
           fst_cmd=fst_cmd(i+1:)
         else
           nc=nc+1
           cmds(nc)=fst_cmd
           exit
         end if
       end do

       !Update the atom parameters
       call VState_to_AtomsPar(A,mode="V")

       call get_logunit(lun)
       open(unit=lun, file=trim(file_fst), status="replace",action="write",position="rewind",iostat=ier)
       if(ier == 0) then
        write(unit=lun,fmt="(a)") "TITLE  FST-file generated with Write_FST"
        write(unit=lun,fmt="(a)") "SPACEG "//trim(Spg%SPG_Symb)
        write(unit=lun,fmt="(a,3f12.5,3f8.3)") "CELL ",cell%cell, cell%ang
        write(unit=lun,fmt="(a)") "BOX   -0.20  1.20   -0.20  1.20    -0.20  1.20 "
        do i=1,A%natoms
          write(unit=lun,fmt="(a,a,3f12.5)")"Atom "//A%atom(i)%lab, A%atom(i)%chemsymb, A%atom(i)%x
        end do
        do i=1,nc
          write(unit=lun,fmt="(a)") trim(cmds(i))
        end do
        write(unit=lun,fmt="(a)") "!"
        call flush(lun)
        close(unit=lun)
       end if
       return
    End Subroutine Write_FST


