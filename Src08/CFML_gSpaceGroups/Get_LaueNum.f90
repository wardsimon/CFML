!!----
!!----
!!----
!!----
SubModule (CFML_GSpaceGroups) Spg_051
   Contains
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
    
End SubModule Spg_051   