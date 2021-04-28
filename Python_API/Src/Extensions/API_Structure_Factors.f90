! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/API_Structure_Factors.f90
! @brief     CFML Structure Factors Fortran binding
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_Structure_Factors

  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use CFML_GlobalDeps,                 only: Cp
  use CFML_Structure_Factors,          only: &
       Structure_Factors, &
       Write_Structure_Factors

  use API_Crystallographic_Symmetry, only: &
       Space_Group_Type_p, &
       get_space_group_type_from_arg

  use API_Atom_TypeDef, only: &
       Atom_list_type_p, &
       get_atom_list_type_from_arg

  use API_Reflections_Utilities, only: &
      Reflection_List_type_p, &
      get_reflection_list_from_arg

  use API_IO_Formats, only: &
       Job_info_type_p, &
       get_job_info_type_from_arg

  implicit none

contains

!!$  function structure_factors_structure_factors(self_ptr, args_ptr) result(r) bind(c)
!!$
!!$    type(c_ptr), value :: self_ptr
!!$    type(c_ptr), value :: args_ptr
!!$    type(c_ptr)        :: r
!!$    type(tuple)        :: args
!!$    type(dict)         :: retval
!!$
!!$    integer            :: num_args
!!$    integer            :: ierror
!!$    integer            :: ii
!!$
!!$    type(list)   :: index_obj
!!$    type(object) :: arg_obj
!!$
!!$    type(Atom_list_type_p)          :: atom_list_p
!!$    type(Space_Group_type_p)        :: spg_p
!!$    type(Reflection_List_type_p)    :: reflection_list_p
!!$
!!$    character(len=3)       :: mode
!!$    character(len=16)      :: pattern
!!$    real(kind=cp)          :: lambda
!!$    
!!$    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
!!$    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
!!$    call unsafe_cast_from_c_ptr(args, args_ptr)
!!$    ! Check if the arguments are OK
!!$    ierror = args%len(num_args)
!!$    
!!$    if (num_args /= 3) then
!!$       call raise_exception(TypeError, "structure_factors_structure_factors expects exactly 3 arguments")
!!$       !@atom_list @space_group.as_fortran_object(), @job_info, @reflection_list
!!$       call args%destroy
!!$       return
!!$    endif
!!$
!!$    call get_atom_list_type_from_arg(args, atom_list_p, 0)
!!$
!!$    call get_space_group_type_from_arg(args, spg_p, 1)
!!$
!!$    call get_reflection_list_from_arg(args, reflection_list_p, 2)
!!$    
!!$    call Structure_Factors(atom_list_p%p, spg_p%p, reflection_list_p%p)
!!$
!!$    ierror = dict_create(retval)
!!$    r = retval%get_c_ptr()
!!$
!!$  end function structure_factors_structure_factors

  function structure_factors_structure_factors(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr)        :: r
    type(tuple)        :: args
    type(dict)         :: retval

    integer            :: num_args
    integer            :: ierror
    integer            :: ii

    type(list)   :: index_obj
    type(object) :: arg_obj

    type(Atom_list_type_p)          :: atom_list_p
    type(Space_Group_type_p)        :: spg_p
    type(Reflection_List_type_p)    :: reflection_list_p
    type(job_info_type_p)           :: job_p

    character(len=3)       :: mode
    character(len=16)      :: pattern
    real(kind=cp)          :: lambda
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    
    if (num_args /= 4) then
       call raise_exception(TypeError, "structure_factors_structure_factors expects exactly 4 arguments")
       !@atom_list @space_group, @reflection_list,  @job_info
       call args%destroy
       return
    endif

    call get_atom_list_type_from_arg(args, atom_list_p, 0)

    call get_space_group_type_from_arg(args, spg_p, 1)

    call get_reflection_list_from_arg(args, reflection_list_p, 2)

    call get_job_info_type_from_arg(args, job_p, 3)
    

    select case (job_p%p%patt_typ(1))
    case ("XRAY_2THE", "XRAY_SXTAL", "XRAY_ENER")
       mode = "NUC"
       call Structure_Factors(atom_list_p%p, spg_p%p, reflection_list_p%p, mode) 
    case("NEUT_2THE", "NEUT_SXTAL", "NEUT_TOF" )
       mode = "XRA"
       lambda = job_p%p%lambda(1)%mina
       call Structure_Factors(atom_list_p%p, spg_p%p, reflection_list_p%p, mode, lambda) 
    case default
       write(*,*) 'Default calculation'
       call Structure_Factors(atom_list_p%p, spg_p%p, reflection_list_p%p)
    end select

       

       !!----
    !!---- Subroutine Structure_Factors(Atm,Grp,Reflex,Mode,lambda)
    !!----    type(atom_list_type),               intent(in)     :: Atm    !List of atoms
    !!----    type(space_group_type),             intent(in)     :: Grp    !Space group
    !!----    type(reflection_list_type),         intent(in out) :: Reflex !It is completed on output
    !!----    character(len=*), optional,         intent(in)     :: Mode   !"NUC","ELE" for neutrons, electrons else: XRays
    !!----    real(kind=cp), optional,            intent(in)     :: lambda !Needed for Xrays
    !!----
    !!----    Calculate the Structure Factors from a list of Atoms
    !!----    and a set of reflections. A call to Init_Structure_Factors
    !!----    is a pre-requisite for using this subroutine. In any case
    !!----    the subroutine calls Init_Structure_Factors if SF_initialized=.false.
    !!----
    !!---- Update: February - 2005
    !!
    !call Structure_Factors(atom_list_p%p, spg_p%p, reflection_list_p%p) !, mode, lambda)

    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function structure_factors_structure_factors

  function structure_factors_write_structure_factors(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Reflection_List_type_p)    :: reflection_list_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "write_structure_factors expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_reflection_list_from_arg(args, reflection_list_p)

    !
    call Write_Structure_Factors(6, reflection_list_p%p)

    !
    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function structure_factors_write_structure_factors

end module API_Structure_Factors
