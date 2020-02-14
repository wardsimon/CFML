!!----
!!----
!!----
!!
SubModule (CFML_gSpaceGroups) SPG_016
   Contains

   !!----
   !!---- Get_Generators_from_Str
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Get_Generators_from_Str(StrGen, d, gen, ngen)
      !---- Arguments ----!
      character(len=*),                            intent(in)  :: StrGen
      integer,                                     intent(out) :: d
      character(len=*), dimension(:), allocatable, intent(out) :: gen
      integer,                                     intent(out) :: ngen

      !--- Local variables ---!
      character(len=:), allocatable :: symbol, ListGen
      integer, dimension(20) :: Pos
      integer                :: i,j,k,np
      logical                :: timerev_provided

      !> Init
      ngen=0

      !> Determine the dimension from the generatorList (first operator)
      ListGen=trim(StrGen)//"  "

      i=index(ListGen,";")
      if (i == 0) then
         !> there is only one generator and ";" is not given
         Symbol=ListGen
         ngen=1
      else
         Symbol=ListGen(1:i-1)
      end if
      d=Get_Dimension_SymmOp(Symbol)

      !> List of Generators
      call Get_Separator_Pos(ListGen,";",pos,np)

      !> Verify if there is a final ";" not followed by a generator
      if (np == 0) then !there is only one generator and ";" is not given
         np=1
         Pos(np)=len_trim(ListGen)+2 !add artificially a position for ";"
      end if
      if (len_trim(ListGen(Pos(np)+1:)) == 0) then
         ngen=np  !final ;
      else
         ngen=np+1
      end if
      allocate(gen(ngen))

      !> Obtaining generators
      j=1
      do i=1, np
         k=pos(i)
         gen(i)=ListGen(j:k-1)
         j = k + 1
      end do
      if (ngen > np) then
         gen(ngen)=ListGen(j:)
         i=index(gen(ngen),";")
         if (i /= 0) gen(ngen)(i:i) = " "
      end if

      !> Time reverse
      timerev_provided=.false.
      do i=1,ngen
         call Get_Separator_Pos(gen(i),",",pos,np)
         if (np < d-1) cycle
         timerev_provided=.true.
      end do

      !> Add time inversion in those operators that have not been Read
      if (timerev_provided) then
         do i=1,ngen
            call Get_Separator_Pos(gen(i),",",pos,np)
            if (np < d-1) gen(i)=trim(gen(i))//",1"
         end do
      end if
   End Subroutine Get_Generators_from_Str

End SubModule SPG_016

