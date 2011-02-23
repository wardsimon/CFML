!
! The painfull (and better) way to call CFML functions from C:
!   Have access to CFML types and functions directly from C, using a binding wrapper
!

module BindCFML
  use ISO_C_BINDING

Type, public :: F_Crystal_Cell_Type
     real(kind=C_FLOAT),dimension(3)   :: cell, ang
     real(kind=C_FLOAT),dimension(3)   :: cell_std, ang_std
     real(kind=C_FLOAT),dimension(3)   :: rcell, rang
     real(kind=C_FLOAT),dimension(3,3) :: GD,GR
     real(kind=C_FLOAT),dimension(3,3) :: Cr_Orth_cel
     real(kind=C_FLOAT),dimension(3,3) :: Orth_Cr_cel
     real(kind=C_FLOAT),dimension(3,3) :: BL_M
     real(kind=C_FLOAT),dimension(3,3) :: BL_Minv
     real(kind=C_FLOAT)                :: CellVol
     real(kind=C_FLOAT)                :: RCellVol
     character (len=1)                 :: CartType
 End Type F_Crystal_Cell_Type


  private :: CToF_Crystal_Cell_Type, FToC_Crystal_Cell_Type
  public  :: BindC_COMPLAINS_THAT
contains

!-----
!----- Private subroutines
!-----


  subroutine CToF_Crystal_Cell_Type(f_in, f_out)
    use CFML_Crystal_Metrics,           only: Crystal_Cell_Type !, Write_Crystal_Cell
    implicit none
    type(F_Crystal_Cell_Type), intent(in)                      :: f_in
    type(Crystal_Cell_Type), intent(out)                       :: f_out

    f_out%cell = f_in%cell
    f_out%ang = f_in%ang
    f_out%cell_std = f_in%cell_std
    f_out%ang_std = f_in%ang_std
    f_out%rcell = f_in%rcell
    f_out%rang = f_in%rang
    f_out%GD = f_in%GD
    f_out%GR = f_in%GR
    f_out%Cr_Orth_cel = f_in%Cr_Orth_cel
    f_out%Orth_Cr_cel = f_in%Orth_Cr_cel
    f_out%BL_M = f_in%BL_M
    f_out%BL_Minv = f_in%BL_Minv

    return
  end subroutine CToF_Crystal_Cell_Type

  subroutine FToC_Crystal_Cell_Type(f_in, f_out)
    use CFML_Crystal_Metrics,           only: Crystal_Cell_Type !, Write_Crystal_Cell
    implicit none
    type(Crystal_Cell_Type), intent(in)                          :: f_in
    type(F_Crystal_Cell_Type), intent(out)                       :: f_out

    f_out%cell = f_in%cell
    f_out%ang = f_in%ang
    f_out%cell_std = f_in%cell_std
    f_out%ang_std = f_in%ang_std
    f_out%rcell = f_in%rcell
    f_out%rang = f_in%rang
    f_out%GD = f_in%GD
    f_out%GR = f_in%GR
    f_out%Cr_Orth_cel = f_in%Cr_Orth_cel
    f_out%Orth_Cr_cel = f_in%Orth_Cr_cel
    f_out%BL_M = f_in%BL_M
    f_out%BL_Minv = f_in%BL_Minv

    return
  end subroutine FToC_Crystal_Cell_Type


!-----
!----- Public subroutines
!-----


end module BindCFML
