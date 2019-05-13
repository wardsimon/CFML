!!----
!!----
!!----
SubModule (CFML_gSpaceGroups) Spg_050
   Contains
   
   !!----
   !!---- GET_POINTGROUP_NUM
   !!----
   !!----    Obtain the ordinal number corresponding to the Point Group
   !!----    symbol according to Point_Group array. Zero if Error is present
   !!----
   !!---- 11/05/2019 
   !!
   Module Function Get_PointGroup_Num(Str_PG) Result(N)
      !---- Arguments ----!
      character(len=*), intent (in) :: Str_PG   ! String containing the PG information
      integer                       :: N        ! Return value on the vector POINT_GROUP 

      !---- Local Variables ----!
      integer                       :: i
      character(len=:), allocatable :: pg

      N=0
      pg=adjustl(Str_PG)

      do i=1,42
         if (pg(1:5) == POINT_GROUP(i)) then
            N=i
            exit
         end if
      end do

      !> return previous numbers for m3 and m3m
      select case (N)
         case (40) ! m3 -> m-3
            N=36

         case (41) ! m3m -> m-3m
            Ni=39

         case (42) ! -3m1 -> -3m
            N=23
      end select
      
   End Function Get_PointGroup_Num
   
   !!----
   !!---- GET_POINTGROUP_STR
   !!----
   !!----    Obtain the string for the Point Group. Blank if error
   !!----
   !!---- 11/05/2019 
   !!
   Module Function Get_PointGroup_Str(N) Result(Str_PG)
      !---- Arguments ----!
      integer,          intent( in) :: N
      character(len=:), allocatable :: Str_PG

      Str_PG="  "
      if (Ni < 1 .or. N > 42) return

      Str_PG=POINT_GROUP(N)
   End Function Get_PointGroup_Str
    
End SubModule Spg_050   