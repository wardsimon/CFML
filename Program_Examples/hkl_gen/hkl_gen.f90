!!----
!!---- Program: HKL_GEN
!!----          Example of simple program using CFML
!!----
!!---- Author: Juan Rodriguez-Carvajal
!!---- Revision: Nov-2008
!!
Program Test_HKL_GEN
   !---- Use Modules ----!
   use CFML_Math_General,               only: asind
   use CFML_Crystallographic_Symmetry,  only: Space_Group_Type, set_SpaceGroup, &
                                              Write_SpaceGroup
   use CFML_Crystal_Metrics,            only: Crystal_Cell_Type,set_crystal_cell,Write_Crystal_Cell
   use CFML_Reflections_utilities,      only: Reflect_Type, get_maxnumref, HKL_uni !HKL_gen

   !---- Variables ----!
   implicit none

   character(len=1)           :: car
   character(len=12)          :: name_file
   character(len=20)          :: spgr
   character(len=35)          :: texto
   type (Space_Group_Type)    :: grp_espacial
   type (Crystal_Cell_Type)   :: cell
   integer                    :: i,num , ier, MaxNumRef
   real,dimension(3)          :: celda, angulo
   real                       :: val1,val2,lambda,angle_2theta
   logical                    :: info

   type (Reflect_Type),allocatable, dimension(:) :: reflections

   !---- Initializing ----!
   info=.true.
   num=0

   !---- Procedure ----!
   do
      write(unit=*,fmt="(/,a)") " ==================================="
      write(unit=*,fmt="(a)")   "     List of unique reflections     "
      write(unit=*,fmt="(a)")   " ==================================="
      write(unit=*,fmt="(a)")   " "

      if (info) then
         write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
         read(unit=*,fmt="(a)") spgr
         if (len_trim(spgr)==0) exit

         !> Set the SpaceGroup Information
         call set_spacegroup(spgr,grp_espacial)

         !> Load Cell Parameters according to the Crystal System
         !> defined by the Space Group
         Select Case(grp_espacial%CrystalSys)
            case("Triclinic")
               write(unit=*,fmt="(a)",advance="no") " => Cell parameters (a,b,c,alpha,beta,gamma): "
               read(unit=*,fmt=*,iostat=ier) celda(1),celda(2),celda(3),angulo(1),angulo(2),angulo(3)

            case("Monoclinic")
               write(unit=*,fmt="(a)",advance="no") " => Cell parameters (a,b,c,beta): "
               angulo(1)=90.0
               angulo(3)=90.0
               read(unit=*,fmt=*,iostat=ier) celda(1),celda(2),celda(3),angulo(2)

            case("Orthorhombic")
               write(unit=*,fmt="(a)",advance="no") " => Cell parameters (a,b,c): "
               angulo(1:3)=90.0
               read(unit=*,fmt=*,iostat=ier) celda(1),celda(2),celda(3)

            case("Tetragonal")
               write(unit=*,fmt="(a)",advance="no") " => Cell parameters (a,c): "
               angulo(1:3)=90.0
               read(unit=*,fmt=*,iostat=ier) celda(1),celda(3)
               celda(2)=celda(1)

            case("Rhombohedral","Hexagonal")
               write(unit=*,fmt="(a)",advance="no") " => Cell parameters (a,c): "
               angulo(1:2)=90.0
               angulo(3)=120.0
               read(unit=*,fmt=*,iostat=ier) celda(1),celda(3)
               celda(2)=celda(1)

            case("Cubic")
               write(unit=*,fmt="(a)",advance="no") " => Cell parameter (a): "
               angulo(1:3)=90.0
               read(unit=*,fmt=*,iostat=ier) celda(1)
               celda(2)=celda(1)
               celda(3)=celda(1)
         End Select
         if( ier /= 0 ) cycle

         !> Set the Cell parameters
         call set_crystal_cell(celda,angulo,cell)

         !> Wavelength
         info=.false.
         write(unit=*,fmt="(a)",advance="no") " => Give the wavelength: "
         read(unit=*,fmt=*) lambda
         !Resolution sphere d*(max)=2.0/Lambda => maximum admissible sinTheta/Lambda =1/Lambda
      end if

      !> Set the range for HKL calculation
      write(unit=*,fmt=*) " "
      write(unit=*,fmt="(a)",advance="no") " => Interval in Sin_Theta/Lambda (2 reals, if val1 < 0 => stops): "
      read(unit=*,fmt=*,iostat=ier) val1,val2
      if(ier /= 0) cycle
      if (val1 < 0.0) exit
      if(val2 > 1.0/Lambda ) then
        val2 = 1.0/Lambda
        write(unit=*,fmt="(a,f8.4)")  " => The maximum Sin_Theta/Lambda available for the current wavelength is: ",val2
        write(unit=*,fmt="(a,2f8.4,a)") " => Interval in Sin_Theta/Lambda is then changed to: (",val1,val2,")"
      end if

      texto = "  1/Angtroms (Sin_Theta/Lambda)"
      car="s"

      !> Allocating the vectors
      MaxNumRef = get_maxnumref(val2+0.2,Cell%CellVol,mult=grp_espacial%Multip)
      write(unit=*,fmt="(a,i10)") " => Maximum number of reflections: ", MaxNumRef
      if (allocated(reflections)) then
        deallocate(reflections)
        allocate (reflections(MaxNumRef))
      else
        allocate (reflections(MaxNumRef))
      end if

      !> Procedure to calculation the all reflections
      !call HKL_GEN(cell,grp_espacial,.true.,val1,val2,num,reflections) !Not ordered
      call HKL_UNI(cell,grp_espacial,.true.,val1,val2,car,num,reflections) !Ordered

      !> Output Information
      write(unit=*,fmt="(a)",advance="no") " => Name of the output file: "
      read(unit=*,fmt="(a)") name_file

      open(unit=1,file=trim(name_file),status="REPLACE",action="WRITE")
      write(unit=1,fmt="(a)")   "    =========================="
      write(unit=1,fmt="(a)")   "    LIST OF UNIQUE REFLECTIONS"
      write(unit=1,fmt="(a,/)") "    =========================="

      !> Write the SpaceGroup information
      call Write_SpaceGroup(grp_espacial,1)

      !> Write the Cell parameters information
      call Write_Crystal_Cell(Cell,1)

      !> Write the Reflections
      write(unit=1,fmt="(/,a,2f8.4,a,/)") " => List of reflections within: ",val1,val2,texto
      do i=1,num
         angle_2theta=2.0*asind(reflections(i)%S*lambda)
         write(unit=1,fmt="(3i4,i5,2f10.5,i8)") reflections(i)%h, reflections(i)%mult, &
                                                reflections(i)%S, angle_2theta,i
      end do
      close(unit=1)
   end do

   write(unit=*,fmt="(a)") " => Program finished ... "
End Program Test_HKL_GEN
