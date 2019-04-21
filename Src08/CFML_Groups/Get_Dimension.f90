!!----
!!----
!!----
!!
SubModule (CFML_Groups) CFML_GRP_003
   Contains
   !!----
   !!---- GET_DIMENSION_GENER
   !!----
   !!---- 19/04/2019
   !!
   Module Function Get_Dimension_Gener(Symb) Result(d)
      !---- Arguments ----! 
      character(len=*), intent(in) :: Symb
      integer                      :: d
      
      !---- Local Variables ----!
      integer, dimension(10) :: Pos
      integer                :: np,i,j
      
      !> Init
      d=0
      
      call Get_Separator_Pos(symb, ",", pos, np)
      if (np <= 0) return
      
      !> Verify if time reversal symbol is provided
      !> reading an integer supposed to be invt
      read(unit=symb(pos(np)+1:),fmt=*, iostat=i) j 

      if (i == 0) then
         d=np+1     !> time_inv provided  
      else
         d=np+2     !>time_inv not-provided
      end if
      
      return
   End Function Get_Dimension_Gener

End SubModule CFML_GRP_003   
   
