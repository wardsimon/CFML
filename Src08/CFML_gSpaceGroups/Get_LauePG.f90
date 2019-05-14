!!----
!!----
!!----
SubModule (CFML_gSpaceGroups) Spg_055
   Contains
   
   !!----
   !!---- GET_LAUE_PG
   !!----    Subroutine to get the information of Laue and Point Group
   !!----    Acta Cryst. A55, (1999) 383-395.
   !!----
   !!---- 13/05/2019 
   !!
   Module Subroutine Get_Laue_PG(Ops, nops, Centro, Laue, Pg)
      !---- Arguments ----!
      type(Symm_Oper_Type), dimension(:), intent(in) :: Ops    ! Reduced operators (Numops)
      integer,                            intent(in) :: NOps   ! Numops
      logical,                            intent(in) :: Centro ! .True. if centred consideration
      character(len=*),                   intent(out):: Laue   ! Laue string 
      character(len=*),                   intent(out):: Pg     ! Point group string          

      !---- Local variables ----!
      integer :: nrot_1, nrot_1b
      integer :: nrot_2, nrot_2b
      integer :: nrot_3, nrot_3b
      integer :: nrot_4, nrot_4b
      integer :: nrot_6, nrot_6b
      integer :: i,n_m,ndet,ind

      !> Init
      Laue=" "
      Pg=" "
      if (Nops ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Get_Laue_PG@GSPACEGROUPS: The symmetry operators is zero!"
         return
      end if  

      nrot_1  = 0
      nrot_2  = 0
      nrot_3  = 0
      nrot_4  = 0
      nrot_6  = 0
      nrot_1b = 0
      nrot_2b = 0
      nrot_3b = 0
      nrot_4b = 0
      nrot_6b = 0
      n_m = 0

      do i=1, Nops
         ndet= Get_Rotation_Order(Ops(i)%Mat(1:3,1:3))
         select case (ndet)
            case (-6)
               nrot_6b=nrot_6b +1
            case (-4)
               nrot_4b=nrot_4b +1
            case (-3)
               nrot_3b=nrot_3b +1
            case (-2)
               nrot_2b=nrot_2b +1
            case (-1)
               nrot_1b=nrot_1b +1
            case ( 1)
               nrot_1 =nrot_1  +1
            case ( 2)
               nrot_2 =nrot_2  +1
            case ( 3)
               nrot_3 =nrot_3  +1
            case ( 4)
               nrot_4 =nrot_4  +1
            case ( 6)
               nrot_6 =nrot_6  +1
            case default
               err_CFML%Ierr=1
               err_CFML%Msg="Get_Laue_PG@GSPACEGROUPS: Problems in Rotation order calculation!"
               return
         end select
      end do

      n_m = nrot_1  + nrot_2  + nrot_3  + nrot_4  + nrot_6  + &
            nrot_1b + nrot_2b + nrot_3b + nrot_4b + nrot_6b

      !> Cubic 
      if ( (nrot_3 + nrot_3b == 8) ) then
         select case (n_m)
            case (12)
               if (.not. centro) then
                  pg="23"
               else
                  pg="m-3"
               end if
               laue="m-3"

            case (24)
               if (.not. centro) then
                  if (nrot_4  == 6) pg="432"
                  if (nrot_4b == 6) pg="-43m"
               else
                  pg="m-3m" 
               end if  
               laue="m-3m"
         end select
         return
      end if   

      !> Hexagonal
      if ( (nrot_6 + nrot_6b == 2) ) then
         select case (n_m)
            case (6)
               if (.not. centro) then
                  if (nrot_6  == 2) pg="6"
                  if (nrot_6b == 2) pg="-6"
               else
                  pg="6/m"
               end if  
               laue="6/m"

            case (12)
               if (.not. centro) then
                  if (nrot_6 == 2) then
                     if (nrot_2  == 7) pg="622"
                     if (nrot_2b == 6) pg="6mm"
                     
                  else if (nrot_6b == 2) then
                     pg="-6m2" 
                  end if    
               else
                  pg="6/mmm"
               end if 
               laue="6/mmm"
         end select
         return
      end if   

      !> Trigonal
      if ( (nrot_3 + nrot_3b == 2) ) then
         select case (n_m)
            case (3)
               if (.not. centro) then
                  pg="3"
               else
                  pg="-3"
               end if
               laue="-3"

            case (6)
               if (.not. centro) then
                  if (nrot_2 == 3) pg="32"
                  if (nrot_2b== 3) pg="3m"
               else
                  pg="-3m"
               end if  
               laue="-3m"
         end select
         return
      end if   

      !> Tetragonal
      if ( (nrot_4 + nrot_4b == 2) ) then
         select case (n_m)
            case (4)
               if (.not. centro) then
                  if (nrot_4  == 2 ) pg="4"
                  if (nrot_4b == 2 ) pg="-4"
               else
                  pg="4/m"
               end if
               laue="4/m"

            case (8)
               if (.not. centro) then
                  if (nrot_4 == 2) then
                     if (nrot_2 == 5)  pg="422"
                     if (nrot_2b == 4) pg="4mmm"
                  
                  else if (nrot_4b == 2) then
                     pg="-4m2"
                  end if
               else
                  pg="4/mmm"
               end if
               laue="4/mmm"
         end select
         return
      end if
               
      !> Orthorhombic 
      if ( (nrot_2 + nrot_2b == 3) ) then
         if (.not. centro) then
            if (nrot_2  == 3 ) pg="222"
            if (nrot_2b == 2 ) pg="mm2"
         else 
            pg="mmm"
         end if
         laue="mmm"
         return
      end if   

      !> Monoclinic
      if ( (nrot_2 + nrot_2b == 1)  ) then
         if (.not. centro) then
            if (nrot_2  == 1 ) pg="2"
            if (nrot_2b == 1 ) pg="m"
         else
            pg="2/m"
         end if
         laue="2/m"
         return
      end if   

      !> Triclinic 
      if (n_m == 1) then
         if (.not. centro) then
            pg="1"
         else
            pg="-1"
         end if
         laue="-1"
      end if  
   End Subroutine Get_Laue_PG
   
   !!----
   !!---- GET_LAUE_NUM
   !!----
   !!----    Obtain the ordinal number corresponding to the Laue-Class
   !!----    symbol according to Laue_Class array. Zero if error is present
   !!----
   !!---- 11/05/2019 
   !!
   Module Function Get_Laue_Num(Str_Laue) Result(N)
      !---- Arguments ----!
      character(len=*), intent (in) :: Str_Laue
      integer                       :: N

      !---- Local Variables ----!
      integer                       :: i
      character(len=:), allocatable :: laue

      N=0
      laue=adjustl(Str_laue)

      do i=1,16
         if (laue(1:5) == LAUE_CLASS(i)) then
            N=i
            exit
         end if
      end do
      
      if (N==15) N=13
      if (N==16) N=14

   End Function Get_Laue_Num
   
   !!----
   !!---- GET_LAUE_STR
   !!----    Obtain the string for the Laue-Class. Blank if error
   !!----
   !!---- 11/05/2019 
   !!
   Module Function Get_Laue_Str(N) Result(Str_Laue)
      !---- Arguments ----!
      integer,          intent( in) :: N
      character(len=:), allocatable :: Str_Laue

      Str_Laue="  " 
      if (N < 1 .or. N > 16) return

      str_laue=LAUE_CLASS(N)
   End Function Get_Laue_Str
   
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
    
End SubModule Spg_055  