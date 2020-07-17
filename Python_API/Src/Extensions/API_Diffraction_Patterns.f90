! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/API_Diffraction_Patterns.f90
! @brief     Diffraction Patterns utilities based on CrysFML
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_Diffraction_Patterns

  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use CFML_Diffraction_Patterns, only: &
       Diffraction_Pattern_Type

  use CFML_API_Calc_Powder_Pattern, only: &
       Calc_Powder_Pattern, &
       Powder_Pattern_Simulation_Conditions_type

  use API_reflections_utilities, only: &
       Reflection_List_type_p, &
       get_reflection_list_from_arg

  implicit none

  type Diffraction_Pattern_type_p
     type(Diffraction_Pattern_type), pointer :: p
  end type Diffraction_Pattern_type_p

  type Powder_Pattern_Simulation_Conditions_type_p
     type(Powder_Pattern_Simulation_Conditions_type), pointer :: p
  end type Powder_Pattern_Simulation_Conditions_type_p

contains

  subroutine get_powpat_from_args(args, pow_pat_sim_c_p, indx)
    type(tuple)                            :: args
    type(Powder_Pattern_Simulation_Conditions_type_p), intent(out)  :: pow_pat_sim_c_p
    integer, optional                      :: indx

    type(object) :: arg_obj
    type(dict) :: arg_dict

    integer :: ierror
    integer :: ii
    type(object) :: t

    if (present(indx)) then
       ierror = args%getitem(arg_obj, indx)
    else
       ierror = args%getitem(arg_obj, 0)
    endif

    ierror = cast(arg_dict, arg_obj)

    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%job   , "job")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%Lambda, "lambda")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%U     , "u_resolution")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%V     , "v_resolution")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%W     , "w_resolution")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%X     , "x_resolution")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%Ls    , "lorentzian_size")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%Thmin , "theta_min")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%Thmax , "theta_max")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%step  , "theta_step")
    ierror = arg_dict%getitem(pow_pat_sim_c_p%p%bkg   , "background")

  end subroutine get_powpat_from_args


  subroutine get_diffraction_pattern_type_from_arg(args, diffraction_pattern_type_pointer, indx)
    type(tuple)                            :: args
    type(Diffraction_Pattern_type_p), intent(out)  :: diffraction_pattern_type_pointer
    integer, optional                      :: indx

    type(object) :: arg_obj
    type(list) :: arg_list
    integer :: diffraction_pattern_type_p12(12)

    integer :: ierror
    integer :: ii
    type(object) :: t

    if (present(indx)) then
       ierror = args%getitem(arg_obj, indx)
    else
       ierror = args%getitem(arg_obj, 0)
    endif

    ierror = cast(arg_list, arg_obj)
    do ii=1,12
       ierror = arg_list%getitem(t, ii-1)
       ierror = cast(diffraction_pattern_type_p12(ii), t)
    enddo
    diffraction_pattern_type_pointer = transfer(diffraction_pattern_type_p12, diffraction_pattern_type_pointer)

  end subroutine get_diffraction_pattern_type_from_arg


  function diffraction_patterns_compute_powder_pattern(self_ptr, args_ptr) result(r) bind(c)

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

    type(Reflection_List_type_p)     :: hkl_p
    type(Diffraction_Pattern_type_p) :: pattern_p
    type(Powder_Pattern_Simulation_Conditions_type_p) :: powpat_p
    integer                          :: pattern_p12(12)

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)

    if (num_args /= 3) then
       call raise_exception(TypeError, "Calc_Powder_Pattern expects exactly 3 arguments")
       !Powpat_conditions, Reflectionlist, Scalef
       call args%destroy
       return
    endif

    allocate(powpat_p%p)
    call get_powpat_from_args(args, powpat_p, 0)
    call get_reflection_list_from_arg(args, hkl_p, 1)
    ierror = args%getitem(arg_obj, 2)
    ierror = cast_nonstrict(powpat_p%p%scalef, arg_obj)
    allocate(pattern_p%p)
    call Calc_Powder_Pattern(powpat_p%p,hkl_p%p,pattern_p%p)
    deallocate(powpat_p%p)

    pattern_p12 = transfer(pattern_p,pattern_p12)
    ierror = list_create(index_obj)
    do ii=1,12
       ierror = index_obj%append(pattern_p12(ii))
    end do
    ierror = dict_create(retval)
    ierror = retval%setitem("address", index_obj)

    r = retval%get_c_ptr()

  end function diffraction_patterns_compute_powder_pattern

  function diffraction_patterns_get_x(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Diffraction_Pattern_type_p) :: diffraction_pattern_type_pointer

    type(ndarray) :: x_array

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "diffraction_patterns_get_x expects exactly 1 argument")
       call args%destroy
       return
    endif

    ! Doing boring stuff
    call get_diffraction_pattern_type_from_arg(args, diffraction_pattern_type_pointer)
    ierror = ndarray_create(x_array, diffraction_pattern_type_pointer%p%x)
    ierror = dict_create(retval)
    ierror = retval%setitem("x", x_array)
    r = retval%get_c_ptr()

  end function diffraction_patterns_get_x

  function diffraction_patterns_get_y(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Diffraction_Pattern_type_p) :: diffraction_pattern_type_pointer

    type(ndarray) :: y_array

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "diffraction_patterns_get_y expects exactly 1 argument")
       call args%destroy
       return
    endif

    ! Doing boring stuff
    call get_diffraction_pattern_type_from_arg(args, diffraction_pattern_type_pointer)
    ierror = ndarray_create(y_array, diffraction_pattern_type_pointer%p%ycalc)
    ierror = dict_create(retval)
    ierror = retval%setitem("y", y_array)
    r = retval%get_c_ptr()

  end function diffraction_patterns_get_y


end module API_Diffraction_Patterns
 
