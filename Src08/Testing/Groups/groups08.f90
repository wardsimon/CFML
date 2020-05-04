 Module Get_gSpG
   use CFML_GlobalDeps
   use CFML_gSpaceGroups
   use CFML_Rational
   use CFML_Maths
   implicit none
   private
   Public :: Set_gSpG
   integer,          dimension(0:2), parameter :: CENT=[2,1,2]
   character(len=1), dimension(10),  parameter :: XYZ=["x","y","z","t","u","v","w","p","q","r"]
   contains
     Subroutine Set_gSpG(str,SpG,Mode,Setting)
       character(len=*),          intent(in) :: str
       Class(SpG_Type),           intent(out):: SpG
       Character(len=*),          intent(in) :: Mode
       Character(len=*),optional, intent(in) :: Setting
       ! --- Local Variables ---!

       Select Case (trim(Mode))
         Case("SHUBN","SUPER")
           if(present(setting)) then
             call Set_SpaceGroup(Str,Mode,SpG,Setting=Setting)
           else
             call Set_SpaceGroup(Str,Mode,SpG)
           end if
         Case default
           call Set_SpaceGroup(Str,SpG)
           if(present(setting)) then
             call Change_Setting_SpaceG(setting,SpG)
             !call Change_Set_SpaceG(setting,SpG)
           end if
       End Select
     End Subroutine Set_gSpG

   Subroutine Change_Set_SpaceG(setting, SpaceG,xyz_type)
      !---- Arguments ----!
      character(len=*),           intent(in )    :: setting
      class(spg_type),            intent(in out) :: SpaceG
      character(len=*), optional, intent(in )    :: xyz_type

      !---- Local variables ----!
      Type(rational), dimension(SpaceG%D,SpaceG%D)     :: Pmat,invPmat
      Type(rational), dimension(SpaceG%D-1,SpaceG%D-1) :: rot,roti,identd
      Type(rational), dimension(SpaceG%D-1)            :: v
      Type(rational), dimension(:,:),allocatable       :: newLat
      Type(rational) :: det
      integer :: i,j,k,l,n,m,Npos,d,Dd,im
      character(len=6) :: Strcode
      character(len=80), dimension(:),allocatable :: gen_lat
      type(spg_type)   :: SpG, P1_t !Auxiliary space groups
      logical          :: centring

      Dd=SpaceG%D
      d=Dd-1
      call Get_Mat_From_Symb(setting, Pmat)
      if(err_CFML%Ierr /= 0) return
      rot=Pmat(1:d,1:d)
      det=Rational_Determ(rot)
      if(det < 0_LI ) then
         err_CFML%Ierr=1
         err_CFML%Msg ="The determinant of the transformation matrix should be positive"
         return
      end if
      L=(SpaceG%Num_Lat+d)*nint(det)
      allocate(newLat(SpaceG%D-1,L))
      newLat=0//1
      invPmat=Rational_Inverse_Matrix(Pmat)
      SpaceG%setting=trim(setting)
      SpaceG%mat2std=Get_Symb_From_Mat(invPmat,StrCode="abc" )
      centring=.false.
      Strcode="xyz"
      if(present(xyz_type)) Strcode=trim(xyz_type)
      roti=Rational_Inverse_Matrix(rot)
      call rational_identity_matrix(identd)

      L=0
      if (SpaceG%Num_Lat > 0) then  !Original lattice is centered
         do_i:do i=1,SpaceG%Num_Lat      !Transform the centring vectors to the new lattice
            v=Rational_Modulo_Lat(matmul(roti,SpaceG%Lat_tr(:,i)))
            if (sum(v) == 0_LI) cycle
            do j=1,L
               if(Rational_Equal(v,newlat(:,j))) cycle do_i
            end do
            L=L+1
            newlat(:,L)=v
            !write(*,"(i8,a,10a)") L," -> ",(trim(rational_string(v(j)))//"  ",j=1,d)
         end do do_i
      end if

      do_i2:do i=1,d  !Test also the basis vectors of the original setting
        v=Rational_Modulo_Lat(roti(1:d,i))
        if (sum(v) == 0_LI) cycle
            do j=1,L
               if(Rational_Equal(v,newlat(:,j))) cycle do_i2
            end do
        L=L+1
        newlat(:,L)=v
        !write(*,"(i8,a,10a)") L," -> ",(trim(rational_string(v(j)))//"  ",j=1,d)
      end do do_i2

      if(L > 0) then !Generate the group P1 with the primitive set of lattice centring in order to get its total number
        allocate(gen_lat(L))
        do i=1,L
          gen_lat(i)=" "
          do j=1,d
            gen_lat(i)=trim(gen_lat(i))//XYZ(j)//"+"//trim(rational_string(newlat(j,i)))//","
          end do
          n=len_trim(gen_lat(i))
          gen_lat(i)(n+1:n+1)="1"
        end do
        call group_constructor(gen_lat,P1_t)
        L=P1_t%Num_Lat
        do i=1,L
          newlat(:,i)=P1_t%Lat_tr(:,i)
        end do
      end if

      SpG%Num_Lat=L
      if(SpG%Num_Lat > 0) centring=.true.

      Npos=SpaceG%Numops*cent(SpaceG%centred)
      SpG%Multip=Npos*(1+SpG%Num_Lat)
      Spg%NumOps=SpaceG%NumOps

      call Allocate_SpaceGroup(Dd, SpG%Multip, SpG)

      if(SpG%Num_Lat > 0) then
        if(allocated(SpG%Lat_tr)) deallocate(SpG%Lat_tr)
        allocate(SpG%Lat_tr(d,SpG%Num_Lat))
        SpG%Lat_tr(:,:)=newlat(:,1:L)
      end if

      do i=1,Npos          !Transform the first Npos operators
        SpG%Op(i)%Mat=matmul(matmul(invPmat,SpaceG%Op(i)%Mat),Pmat)
        SpG%Op(i)%Mat(1:d,Dd)=Rational_Modulo_Lat(SpG%Op(i)%Mat(1:d,Dd))
        SpG%Op(i)%time_inv=SpaceG%Op(i)%time_inv
        SpG%Op(i)%dt=SpaceG%Op(i)%dt
        SpG%Symb_Op(i)=Get_Symb_from_Mat(SpG%Op(i)%Mat, Strcode,SpG%Op(i)%time_inv)
      end do

      if(centring) then
        n=Npos
        do L=1,SpG%Num_Lat
           do k=1,Npos
             SpG%Op(k+n)%Mat=SpG%Op(k)%Mat             !Same matrix as the first NumOps operators
             SpG%Op(k+n)%time_inv=SpG%Op(k)%time_inv   !Same time inversion
             SpG%Op(k+n)%Mat(1:d,Dd)=Rational_Modulo_Lat(SpG%Op(k)%Mat(1:d,Dd)+SpG%Lat_tr(:,L)) !Different translation
             SpG%Op(k+n)%dt=SpG%Op(k)%dt
             SpG%Symb_Op(k+n)=Get_Symb_from_Mat(SpG%Op(k+n)%Mat, Strcode,SpG%Op(k+n)%time_inv)
           end do
           n=n+Npos
        end do
      end if

      !Reallocate the array components of the initial space group
      call Allocate_SpaceGroup(Dd, SpG%Multip, SpaceG)
      SpaceG%Multip=SpG%Multip
      SpaceG%Op=SpG%Op             !Copy the array components keeping the names,labels, etc.
      SpaceG%Symb_Op=SpG%Symb_Op
      SpaceG%Num_Lat=Spg%Num_lat
      SpaceG%Num_aLat=SpaceG%Num_aLat*(1+SpG%Num_Lat)

      if(SpaceG%Num_Lat > 0) then
         if(allocated(SpaceG%Lat_tr)) deallocate(SpaceG%Lat_tr)
         allocate(SpaceG%Lat_tr(d,SpaceG%Num_Lat))     ! Centring vectors
         SpaceG%Lat_tr=SpG%Lat_tr                      ! Direct assignment of vectors
      end if

      if(SpaceG%Num_aLat > 0) then
         if(allocated(SpaceG%aLat_tr)) deallocate(SpaceG%aLat_tr)
         allocate(SpaceG%aLat_tr(d,SpaceG%Num_aLat))     ! Anti-translation vectors
         m=0
         if(SpaceG%Mag_Type /= 2) then
           do k=1,SpaceG%multip
             if(rational_equal(SpaceG%Op(k)%Mat(1:d,1:d),identd) .and. SpaceG%Op(k)%time_inv == -1) then
               m=m+1
               SpaceG%aLat_tr(:,m) = SpaceG%Op(k)%Mat(1:d,Dd)
             end if
           end do
         end if
         SpaceG%Num_aLat=m
      end if
      !Transform the symmetry operators of the list of generators !TO DO

   End Subroutine Change_Set_SpaceG


 End Module Get_gSpG

!!----
!!----
!!----
!!----
 Program Test_Groups
    !---- Use Modules ----!
    use CFML_Globaldeps, only: cp, err_cfml
    use CFML_Symmetry_Tables
    use CFML_gSpaceGroups
    use Get_gSpG

    character(len=256)                  :: generatorList
    character(len=180)                  :: setting
    !character(len=25)                   :: forma="(i5,tr2,a,   i4,a,i8)"
    character(len=25)                   :: mode
    character(len=5)                    :: aux
    type(Spg_Type)                      :: Grp
    type(Spg_Type), dimension(512)      :: sGrp
    !integer, dimension(:,:),allocatable :: table
    !integer, dimension(:,:),allocatable :: G
    !integer, dimension(:),  allocatable :: ord
    integer :: i, j, L, nsg, indexg, num_group, ier !, ind
    real(kind=cp) :: start, fin
    logical :: set_given, sup_given, datb_given, full, shub_given

    !> Init
    call Set_Conditions_NumOp_EPS(4096) !Maximum admissible multiplicity

    do
       set_given=.false.; sup_given=.false.; datb_given=.false.; shub_given=.false.
       write(*,'(/,a,/)') " => Examples of input to the next question:"
       write(*,'(a)') "1236   :: a,b,c;0,0,0      shub        <-- Shubnikov group type 1236 in standard setting"
       write(*,'(a)') "123    :: a,-c,b;1/2,0,0   cryst       <-- Space group type 123 in non-standard setting"
       write(*,'(a)') "Pn'ma                                  <-- Shubnikov group type Pn'ma in standard setting"
       write(*,'(a)') "Pn'ma  :: -c,b,a;0,0,0     shub        <-- Shubnikov group type Pn'ma in the non-standard setting Pcmn'"
       write(*,'(a)') "B 2 C B                                <-- space #41 of standard symbol Aba2"
       write(*,'(a)') "13232  :: sup                          <-- SuperSpace group 13232 in standard setting"
       write(*,'(a)') "221    :: a2,-a1,a3,a4;0,0,1/2,0  sup  <-- SuperSpace group #221 in non-standard setting"
       write(*,'(a)') "Pnma(0,0,g)ss0  :: sup                 <-- SuperSpace group of standard symbol"
       write(*,'(a)') "x,-y,z,t,1;x,y,z,t+1/2,-1              <-- The generators of a Magnetic SuperSpace group"
       write(*,'(a)') "x1,-x2,x3,x4,1;x1,x2,x3,x4+1/2,-1      <-- The generators of a Magnetic SuperSpace group"
       write(*,'(/,a)',advance='no') " => Introduce generators, HM or Hall symbol, or the number of a standard group as indicated above: "
       read(*,'(a)') generatorList
       if (len_trim(generatorList) == 0) exit

       !> Determine if it is a number
       aux=" "
       read(unit=generatorList,fmt=*,iostat=ier) num_group
       if (ier == 0) then
          !Now determine if setting is appearing in which case it uses data bases
          i=index(generatorList,"::")
          if(i /= 0) then
             j=index(generatorList,"sup")  !Superspace
             if(j /= 0) then
               sup_given=.true.
               if(index(generatorList,"a1") /= 0) then
                 setting= adjustl(generatorlist(i+2:j-1))
                 set_given=.true.
               end if
             else
               j=index(generatorList,"shub")  !Shubnikov group
               if(j /= 0) then
                 shub_given=.true.
                 setting= adjustl(generatorlist(i+2:j-1))
                 if(len_trim(setting) /= 0) set_given=.true.
               else
                 setting= adjustl(generatorlist(i+2:))
                 if(len_trim(setting) /= 0) set_given=.true.
               end if
             end if
             generatorlist=generatorlist(1:i-1)
          end if
       else
          i=index(generatorList,"::")
          if(i /= 0) then
             j=index(generatorList,"sup")  !Superspace
             if(j /= 0) then
               sup_given=.true.
               if(index(generatorList,"a1") /= 0) then
                 setting= adjustl(generatorlist(i+2:j-1))
                 set_given=.true.
               end if
             else
               j=index(generatorList,"shub")  !Shubnikov group
               if(j /= 0) then
                 shub_given=.true.
                 setting= adjustl(generatorlist(i+2:j-1))
                 if(len_trim(setting) /= 0) set_given=.true.
               else
                 setting= adjustl(generatorlist(i+2:))
                 if(len_trim(setting) /= 0) set_given=.true.
               end if
             end if
             generatorlist=generatorlist(1:i-1)
          end if
          i=index(generatorList,"'")
          j=index(generatorList,"_")
          if(i /= 0 .or. j /= 0) datb_given=.true.
       end if

       call CPU_TIME(start)
       if(sup_given) then
          if(set_given) then
             call Set_gSpG(generatorList,Grp,"SUPER",Setting)
          else
             call Set_gSpG(generatorList,Grp,"SUPER")
          end if
       else
         if(shub_given) then
            if(set_given ) then
               call Set_gSpG(generatorList,Grp,"SHUBN",Setting)
            else if (datb_given) then
               call Set_gSpG(generatorList,Grp,"SHUBN")
            end if
         else
            if(set_given ) then
               call Set_gSpG(generatorList,Grp,"GEN",Setting)
            else
               call Set_gSpG(generatorList,Grp,"GEN")
            end if
         end if
       end if
       if (Err_CFML%Ierr /= 0) then
          write(*,'(/,4x,a)') trim(Err_CFML%Msg)
          cycle
       else
          call Write_SpaceGroup_Info(Grp)
       end if

       full=.false.
       do
          write(*,'(/,a)',advance='no') "Introduce the index of subgroups (if = 0, no restriction, if < 0 no calculation): "
          read(*,"(a)") mode
          if(index(mode,"full") /= 0) full=.true.
          read(mode,*,iostat=ier) indexg
          if(ier == 0) exit
       end do
       !> Testing Get_subgroups_cosets
       if (indexg == 0) then
         if(full) then
            call get_subgroups_full(Grp,sGrp,nsg)
         else
            call get_subgroups(Grp,sGrp,nsg)
         end if
       else if(indexg > 0) then
         if(full) then
            call get_subgroups_full(Grp,sGrp,nsg,indexg)
         else
          call get_subgroups(Grp,sGrp,nsg,indexg)
         end if
       else
          cycle
       end if
       if (Err_CFML%Ierr /= 0) then
          write(*,'(/,4x,a)') trim(Err_CFML%Msg)
          cycle
       end if

       if (nsg > 0) Then
          do L=1,nsg
             write(*,"(/2(a,i3))") "  SUB-GROUP NUMBER #",L, " of index: ",Grp%multip/sGrp(L)%multip
             call Identify_Group(sGrp(L)) !.false.
             if (Err_CFML%Ierr /= 0) then
                write(*,'(/,4x,"=> Error in the identification of the group: ",a)') trim(Err_CFML%Msg)
             end if
             call Write_SpaceGroup_Info(sGrp(L))
          end do
       end if

       call CPU_TIME(fin)
       write(*,"(a,f12.3,a)") "CPU_TIME for this calculation: ",fin-start," seconds"
    end do

End Program Test_Groups