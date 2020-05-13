
Submodule (CFML_Atoms) Atm_Symmetry_Constraints
  contains
     Module Subroutine Check_Symmetry_Constraints(SpG,Atm)
       class(SpG_Type),   intent(in)     :: SpG
       type(AtList_Type), intent(in out) :: Atm

       !--- Local variables ---!
       integer :: i,codini
       real(kind=cp), dimension(3)   :: codes
       real(kind=cp), dimension(6,8) :: codeT
       codini=1
       do i=1,Atm%natoms
          if(Atm%Atom(i)%Magnetic) then
            codes=1.0
            call Get_moment_ctr(Atm%Atom(i)%x,Atm%Atom(i)%moment,SpG,codini,codes)
          end if
       end do
       Select Type (SpG)
         type is (SuperSpaceGroup_Type)
            do i=1,Atm%natoms
              Select Type(at => Atm%Atom(i))
                class is (MAtm_Std_Type)
                  if(at%n_mc > 0) then
                    codeT=1.0
                    call Get_TFourier_ctr(at%x,at%Mcs(:,1:at%n_mc),codeT(:,1:at%n_mc),SpG,codini,"M")
                  end if
                  if(at%n_dc > 0) then
                    codeT=1.0
                    call Get_TFourier_ctr(at%x,at%Dcs(:,1:at%n_dc),codeT(:,1:at%n_dc),SpG,codini,"D")
                  end if
              End Select
            end do
       End Select
       Atm%symm_checked=.true.
     End Subroutine Check_Symmetry_Constraints

End Submodule Atm_Symmetry_Constraints
