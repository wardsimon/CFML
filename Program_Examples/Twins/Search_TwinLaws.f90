!!----
!!----   PROGRAM:  Search_TwinLaws
!!----   Author: J. Rodriguez-Carvajal (ILL) January 2011, Updated: December 2011
!!----
!!----   This program uses the method of A. Santoro (Acta Cryst A30, 224 (1974)) for
!!----   getting possible twin laws from special metrics.
!!----
!!----   The program looks for solutions of the equation (18) of the cited paper:
!!----                    delta = Bm  GD (Bm)t  - GD
!!----   A solution is accepted if Sum(Abs(delta)) < eps = tol*Sum(Abs(GD))/9
!!----   The matrices Bm are generated by systematically scanning integers values.
!!----   Bm are rational matrices of the form Nij/ni. For details see the cited article.
!!----
!!----   The CPU time depends strongly on the extend of the range for integers for generating
!!----   the rational Bm matrices. The number of solutions depends also on the tolerance
!!----   and the special metric relations between the cell parameters.
!!----   The highest admissible integer is 5.
!!----
!!----   INPUT FILE (My_File.cfl)
!!----   ------------------------
!!----   The input data are simply the unit cell parameters, the symbol of the space group,
!!----   the keyword TOL (for tolerance) followed by its value. If the value is greater
!!----   than 1 it is suposed to be given in percentage and a division by 100 is performed.
!!----   The integer indices defining the range for search are also to be given.
!!----
!!----   The input file should have the extension ".cfl" and an example of its content
!!----   is given below:
!!----
!!----  ---- start of the CFL file in line below
!!----Title  Test of the program Search_TwinLaws
!!----cell   5.486  5.486  11.000 90.0 90.0 90.0
!!----Spgr  I 4/m m m
!!----tol  0.03                     !tolerance for eq (18)     eps=tol*Sum(Abs(GD))/9
!!----indices    -3 3  0 3  -4 4    ! maximum -5 5
!!----  ---- end of the CFL file
!!----
!!----   For running the program from a terminal, one has to invoke it followed by the
!!----   name of the CFL file. (e.g.  twins my_cfl_file.cfl)
!!----
!!----
!!----
!!----   OUTPUT FILES (My_File_twins.cfl and My_File.twins)
!!----   --------------------------------------------------
!!----   Two output files are generated by Search_TWIN_LAWS, they are named according to
!!----   the base-name of the initial CFL file. For the input file My_File.cfl, the two
!!----   files My_File_twins.cfl and My_File.twins are generated. The new CFL file may
!!----   be used by other programs using information about twins (e.g. Esmeralda) for
!!----   predicting peak positions.
!!----
!!----   An example of the output CFL file corresponding to the previous one is given
!!----   below:
!!----
!!----  ---- start of the output CFL file in line below
!!----   Title   CFL-file generated by Search_TWIN_LAWS from test_twin_search.cfl
!!----   Title  Test of the program Search_TwinLaws
!!----   cell   5.486  5.486  11.000 90.0 90.0 90.0
!!----   Spgr  I 4/m m m
!!----   tol  0.03                     !tolerance for eq (18)     eps=tol*Sum(Abs(GD))/9
!!----   indices    -3 3  0 3  -4 4    ! maximum -5 5
!!----
!!----   TWIN_nam  from test_twin_search.cfl
!!----   TWIN_typ  4
!!----   TWIN_rot  -0.6718 -1.5030 -1.3333      131.8103
!!----   TWIN_rot   0.0002  1.4142  0.0000      180.0000
!!----   TWIN_rot  -1.8804 -0.3740 -1.5000      120.0000
!!----   TWIN_rot   0.0000  0.4987 -2.0000       90.0000
!!----   TWIN_rot   1.0000  1.0000  0.0000      180.0000
!!----   TWIN_rot   2.3791  0.3740 -0.5000      120.0000
!!----   TWIN_rot   2.0017  1.3367  0.6667      109.4712
!!----   TWIN_rot   1.8355  1.5030  0.6667       70.5288
!!----     ---- end of the output CFL file
!!----
!!---- The file My_File.twins is the full output files from the program and contains
!!---- self-explanatory information.
!!----
!!---- The twin laws are given as a rotation axis in Cartesian components with respect
!!---- to the crystal and a rotation angle. The Cartesian system used in this program
!!---- corresponds to that defined by the following transformations:
!!----         (a)   (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )  (i)
!!----         (b) = (         0         ,     b sinalpha      , b cosalpha)  (j)
!!----         (c)   (         0         ,         0           , c         )  (k)
!!---- The rotation, corresponding to the given angle, of the reference lattice with respect
!!---- to the given axis produces a quasi-coincidence with the original one.
!!----
!!----
  Program Search_TwinLaws
    use CFML_Crystal_Metrics,  only: Crystal_Cell_Type, Err_Crys, Err_Crys_Mess, Init_err_crys,  &
                                     Change_Setting_Cell,Set_Crystal_Cell, Write_Crystal_Cell,   &
                                     get_primitive_cell
    use CFML_String_utilities, only: l_case
    use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                              Write_SpaceGroup

    use CFML_Math_3D,       only: Invert_A
    use CFML_Geometry_Calc, only: Get_Anglen_Axis_From_RotMat
    use CFML_Math_General,  only: acosd, Equal_Matrix, Trace, Co_prime_Vector
    implicit none
    integer, parameter          :: Max_Sol=2000
    real,    dimension(3)       :: cel1,ang1,axc
    integer, dimension(3,3)     :: Nu,nB
    integer, dimension(3)       :: axis,ax
    integer, dimension(3,3,48)  :: sym !Symmetry operators of the point group
                                       !in a primitive basis
    real, dimension(3,3,Max_Sol):: Bsol !Independent solutions
    real, dimension(3,3)        :: base,Rot,Wi,W,GD,DELTA,Bm,Bmi,Dt,Lm,Lmi
    character(len=1)            :: lat_type, key
    type (Crystal_Cell_Type)    :: cell, pcell
    type (space_group_type)     :: SpG
    logical                     :: arggiven,centred, cell_given=.false., spg_given=.false.
    character(len=20)           :: SpG_symb
    character(len=256)          :: fileinp, fileout, new_cfl, line,line_ini
    integer                     :: i,j,lun=1,lout=2,i_cfl=3,i1,i2,i3,i4,i5,i6,i7,i8,i9,n1,n2,n3
    integer                     :: ia1,ia2,ib1,ib2,ic1,ic2, in1,in2,narg,nop,n,ier
    integer(kind=8)             :: iratio,im,nn
    real                        :: tol,tolf,sm,start,fin,remain,angle


    narg=command_argument_count()
    if(narg /= 0) then
      call get_command_argument(1,fileinp)
      arggiven=.true.
      if(index(fileinp,".cfl") == 0) fileinp=trim(fileinp)//".cfl"
      open(unit=lun,file=trim(fileinp),status="old",action="read",position="rewind",iostat=ier)
    end if

    write(unit=*,fmt="(/,/,6(a,/))")                                                          &
         "                  ------ PROGRAM: Search_TwinLAWS ------                   "  , &
         " **************************************************************************"  , &
         " * This program uses the method of A. Santoro (Acta Cryst A30, 224 (1974))*"  , &
         " * for getting possible twin laws from special metrics of the crystal.    *"  , &
         " **************************************************************************"  , &
         "                         (JRC- ILL, December-2011)                         "
    write(unit=*,fmt=*) " "

    if(ier /=0 .or. narg == 0) then
    ! reading data
      do
        write(unit=*,fmt="(a)",advance="no") " => Enter the code of the input CFL file: "
        read(unit=*,fmt="(a)") fileinp
        fileinp=trim(fileinp)//".cfl"
        open(unit=lun,file=trim(fileinp),status="old",action="read",position="rewind",iostat=ier)
        if(ier /= 0) cycle
        exit
      end do
    end if

    fileout=fileinp(1:len_trim(fileinp)-4)//".twins"
    new_cfl=fileinp(1:len_trim(fileinp)-4)//"_twin.cfl"
    open(unit=lout,file=trim(fileout),status="replace",action="write")
    open(unit=i_cfl,file=trim(new_cfl),status="replace",action="write")

    write(unit=i_cfl,fmt="(a)") "Title   CFL-file generated by Search_TWIN_LAWS from "//trim(fileinp)
    ! Initialisation of indices
     ia1=-5; ia2=5; ib1=0; ib2=5; ic1=-5; ic2=5; in1=1; in2=5
     tol=0.05 !Percentage of the total sum of absolute values of GD

    !Reading the CFL file
    do
      read(unit=lun,fmt="(a)",iostat=ier) line_ini
      if(ier /= 0) exit
      write(unit=i_cfl,fmt="(a)") trim(line_ini)  !Rewrites in the new CFL file
      line=l_case(line_ini)
      i=index(line,"cell")
      if( i /= 0) then
         read(unit=line(5:),fmt=*) cel1,ang1
         cell_given=.true.
      end if
      i=index(line,"spgr")
      if( i /= 0) then
         SpG_Symb=adjustl(line_ini(5:))
         lat_type=SpG_Symb(1:1)
         spg_given=.true.
      end if
      i=index(line,"tol")
      if( i /= 0) then
        read(unit=line(4:),fmt=*) tol
        if(tol > 1.0) tol=tol*0.01
      end if
      i=index(line,"indices")
      if( i /= 0) then
        read(unit=line(8:),fmt=*) ia1,ia2,ib1,ib2,ic1,ic2
        if(ia1 < -5) ia1 = -5;  if(ia2 > 5) ia2 = 5
        if(ib1 < -5) ib1 = -5;  if(ib2 > 5) ib2 = 5
        if(ic1 < -5) ic1 = -5;  if(ic2 > 5) ic2 = 5
      end if
    end do
    !End reading CFL file

    if(.not. cell_given) then
      write(unit=*,fmt="(a)") " => Error no cell given in file: "//trim(fileinp)
      stop
    end if
    if(.not. spg_given) then
      write(unit=*,fmt="(a)") " => Error no space group given in file: "//trim(fileinp)
      stop
    end if

    write(unit=i_cfl,fmt="(/,a)") "TWIN_nam  from "//trim(fileinp)
    write(unit=i_cfl,fmt="(a)")   "TWIN_typ  4 "
    centred=.true.
    if(lat_type == "P" .or. lat_type == "p") centred=.false.

    call Set_Crystal_Cell(cel1,ang1,Cell)
    call Set_SpaceGroup(SpG_Symb,SpG)

    write(unit=lout,fmt="(a)") "                  -------------------------------"
    write(unit=lout,fmt="(a)") "                     PROGRAM: Search_TwinLAWS    "
    write(unit=lout,fmt="(a)") "                  -------------------------------"
    write(unit=lout,fmt="(/,4(a,/))")                                                     &
         " **************************************************************************"  , &
         " * This program uses the method of A. Santoro (Acta Cryst A30, 224 (1974))*"  , &
         " * for getting possible twin laws from special metrics of the crystal.    *"  , &
         " **************************************************************************"
    write(unit=lout,fmt="(a)") "                  Author: J. Rodriguez-Carvajal (ILL)"
    write(unit=lout,fmt="(a)") "            Version 0.1, January 2011 (updated December 2001)"
    write(unit=lout,fmt="(a)") "         ----------------------------------------------------"
    write(unit=lout,fmt="(/,a)") "  => INPUT CELL DATA"
    call Write_Crystal_Cell(Cell,Lout)

    Lm=Transpose(Cell%Orth_Cr_cel)   !L matrices as defined by Santoro
    Lmi=Transpose(Cell%Cr_Orth_cel)

    !TBL=Matmul(transpose(Lm),Matmul(Cell%BL_M,Cell%GD)) !Matrix to pass a column vector referred to
    !Not needed TBL=I in our convention !!!!             !the conventional Cartesian basis to the
                                                         !Busing-Levy reference system
    call Get_Primitive_Cell(lat_type,cell,pcell,Wi)
    Wi=reshape( (/1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0/), (/3,3/) )
    if(centred) then
      W=invert_a(Wi)
      write(unit=lout,fmt="(/,a)") "  => Input subcell data in primitive setting"
      write(unit=lout,fmt="(/,a)") "  => Transformation matrix to primitive cell and inverse:"
      do i=1,3
         write(unit=lout,fmt="(2(a,3f10.4))")"                 ",Wi(i,:),"                 ",W(i,:)
      end do
      call Write_Crystal_Cell(pCell,Lout)
    else
      W=Wi
    end if

    call Write_SpaceGroup(SpG,Lout)
    ! A'  = M  A  => C = (R,T), C'= (R',T') => R'  = inv(Mt) R Mt  ; T' = inv(Mt) T
    ! In our case A' is primitive so M=Wi
    nop=SpG%numops !number of symmetry operators excluding lattice centrings
    if (SpG%centred /= 1) nop=nop*2
    write(unit=Lout,fmt="(/,a)") " => Point Symmetry operators referred to the primitive basis"
    do i=1,nop
       sym(:,:,i) = Nint(Matmul(Matmul(transpose(W),SpG%Symop(i)%rot),transpose(Wi)))
       write(unit=Lout,fmt="(a,i3,a,3(3i3,a))") " =>  SymOp # :",i, "    (",sym(1,:,i)," /",sym(2,:,i)," /",sym(3,:,i)," )"
    end do

    base=transpose(pCell%Cr_Orth_cel)   !Provides a matrix with rows equal to the basis vectors in cartesian components
    GD=pCell%GD  !metric tensor of the primitive cell
    tolf=tol*Sum(abs(GD))/9.0

    iratio= (ia2-ia1+1)*(ia2-ia1+1)*(ia2-ia1+1)*(ib2-ib1+1)*(ib2-ib1+1)*(ib2-ib1+1)*(ic2-ic1+1)*(ic2-ic1+1)*(ic2-ic1+1)
    iratio=iratio*(in2-in1+1)*(in2-in1+1)*(in2-in1+1)
    write(unit=*,fmt="(a,i10)") " => The maximum number of test calculations is ",iratio
    im=iratio/400
    n=0
    call cpu_time(start)
 dox: do i1=ia1,ia2                      !         |i1  i4  i7|
     do i2=ib1,ib2                       !    Nu = |i2  i5  i8|
      do i3=ic1,ic2                      !         |i3  i6  i9|
       do i4=ia1,ia2
        do i5=ib1,ib2                    !         |n1   0   0|
         do i6=ic1,ic2                   !    n  = | 0  n2   0|
          do i7=ia1,ia2                  !         | 0   0  n3|
           do i8=ib1,ib2
            do i9=ic1,ic2
             j=i1*i5*i9+i4*i8*i3+i2*i6*i7-i3*i5*i7-i8*i6*i1-i2*i4*i9     !determinant (much faster than calling determ_A)
             Nu=reshape((/i1,i2,i3,i4,i5,i6,i7,i8,i9/),(/3,3/))
             do n1=in1,in2
              do n2=in1,in2
               doi: do n3=in1,in2
                nn=nn+1
                if(j/(n1*n2*n3) /= 1) cycle !Determinant of B should be 1
                nB=reshape((/ i1/n1,i2/n2,i3/n3,i4/n1,i5/n2,i6/n3,i7/n1,i8/n2,i9/n3 /),(/3,3/))
                Bm(1,:)=real(Nu(1,:))/real(n1)
                Bm(2,:)=real(Nu(2,:))/real(n2)
                Bm(3,:)=real(Nu(3,:))/real(n3)
                do i=1,n !n is the current number of solutions
                  if(Equal_Matrix(Bm, Bsol(:,:,i),3)) cycle doi
                end do

                if(mod(nn,im) == 0)  then
                  write(unit=*,fmt="(a,i12,a,f12.3)") "  Status: ",iratio-nn," tests remaining  ->  SM: ",sm
                  call cpu_time(fin)
                  remain=real(iratio-nn)*(fin-start)/real(nn)
                  write(*,"(a,i6,a)")"  Approximate remaining CPU-Time: ",nint(remain)," seconds"
                end if

                !Fundamental equation for Twinning condition: delta approx. zero
                delta=matmul(Bm,matmul(GD,transpose(Bm)))-GD
                sm=sum(abs(delta))
                if(sm > tolf) cycle

                ! Test of symmetry operations
                do i=1,nop
                  if(Equal_Matrix(nB,sym(:,:,i),3)) cycle doi
                end do
                ! A possible twin law has been found
                n=n+1
                if(n > Max_Sol) exit dox
                Bsol(:,:,n)=Bm
                write(unit=lout,fmt="(/,a,i3)") " => TWIN law number : ",n
                write(unit=lout,fmt="(a)") &
                " =>      Nij       ni       Nint(B)                    B                                Delta:"
                write(unit=lout,fmt="(a,3i3,a,i1,a,3i3,2(5x,3f10.5))")  &
                "     ",Nu(1,:),"      ",n1,"     ",nB(1,:),Bm(1,:),delta(1,:)
                write(unit=lout,fmt="(a,3i3,a,i1,a,3i3,2(5x,3f10.5))")  &
                "     ",Nu(2,:),"      ",n2,"     ",nB(2,:),Bm(2,:),delta(2,:)
                write(unit=lout,fmt="(a,3i3,a,i1,a,3i3,2(5x,3f10.5))")  &
                "     ",Nu(3,:),"      ",n3,"     ",nB(3,:),Bm(3,:),delta(3,:)

                !Determination of the twin axis and angle
                Bmi=invert_A(Bm)
                Dt=Matmul(W,matmul(Bmi,Wi))
                Rot=transpose( Matmul( Lm,matmul(Dt,Lmi) )  )
                Call Get_Anglen_Axis_From_RotMat(Rot,axc,angle)

                write(unit=i_cfl,fmt="(a,3f8.4,f14.4)") "TWIN_rot ",axc,angle
                write(unit=lout,fmt="(a,3f10.4,a)")  &
                " => Cartesian axis of the twin Law:   [ ",axc," ]"
                ax=Nint(100.0*Matmul(Transpose(Lm),axc))
                Call Co_Prime_vector(ax,axis)
                write(unit=lout,fmt="(a,3i4,a,f10.4)")  &
                " => Axis and angle of the twin Law:   [ ",axis," ]    Alpha = ",angle

                write(unit=lout,fmt="(a)") "    ----------------------------------------------------------------------------------"
                write(unit=*,fmt="(a,i5,a,f7.3,a,3i4,a,f10.4)") " => New solution: ",n," SM -> ",sm,&
                "   ->  Axis and angle of the twin Law:   [ ",axis," ]    Alpha = ",angle

               end do doi   !n3
              end do    !n2
             end do    !n1
            end do    !i9
           end do     !i8
          end do      !i7
         end do       !i6
        end do        !i5
       end do         !i4
      end do          !i3
     end do           !i2
    end do  dox       !i1
    call cpu_time(fin)
    write(unit=*,fmt="(/,a,f10.2,a)")  "  CPU-Time: ", fin-start," seconds"
    write(unit=*,fmt="(/,a)") " => Press <enter> to finish "
    read(unit=*,fmt="(a)") key

    stop
  End Program Search_TwinLaws
