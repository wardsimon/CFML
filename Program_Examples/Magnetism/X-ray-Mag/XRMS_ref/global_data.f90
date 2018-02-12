module global_data

   use CFML_crystal_metrics,              only: crystal_cell_type
   use CFML_crystallographic_symmetry,    only: space_group_type
   use CFML_atom_typedef,                 only: atom_list_type, matom_list_type
   use CFML_Magnetic_Symmetry,            only: magsymm_k_type
   use CFML_IO_Formats,                   only: file_list_type
   use CFML_LSQ_TypeDef

   use XRMS_types

   implicit none

   !  Structure variables
   type(crystal_cell_type)          :: unit_cell
   type(space_group_type)           :: space_gr
   type(atom_list_type)             :: atom_list
   type(matom_list_type)            :: magn_atom_list
   type(magsymm_k_type)             :: magn_model
   type(magh_list_type_xrms)        :: sample_refl_list
   type(file_list_type)             :: str_file

   real                                :: scale_fact = 100.0
   real,parameter                      :: lambda = 1.627517 ! In Angstroms
   real,parameter                      :: tau = 0.4258 ! Incommensurate component of vec(k)
   real,dimension(3)                   :: specular_axis
   real,dimension(3,3)                 :: matrix_mom_star
   !  The array matrix_mom_star transforms
   !  the coordinates of the magnetic moment
   !  from the frame of basis (e_1, e_2, e_3) = (vec(a)/a, vec(b)/b, vec(c)/c)
   !  to the reference frame of the reciprocal space.
   !  The coordinates of m_n (magnetic moment) transforms acording to M:
   !  m_n_star_i = M_ij*m_n_j
   
   type(LSQ_data_type)                 :: dat
   type(LSQ_state_vector_type)         :: vs
   type(LSQ_conditions_type),save      :: cond
   
   !  File variables
   character(len=50)                   :: filename, dataname, outname
   character(len=:),allocatable        :: codename
   integer                             :: luo=7

   contains

end module global_data
