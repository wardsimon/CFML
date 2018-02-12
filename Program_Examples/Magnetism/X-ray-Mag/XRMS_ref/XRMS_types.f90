module XRMS_types

   use CFML_GlobalDeps,          only : cp
   
   implicit none
   
   !
   !  Type :: MAGH_TYPE_XRMS
   !
   type, public :: MagH_Type_XRMS
      real(kind=cp),dimension(3)       :: H           ! H +/- k
      real(kind=cp)                    :: Fobs2       ! Squared observed magn structure factor
      real(kind=cp)                    :: SFobs2      ! Sigma of F2obs
      real(kind=cp),dimension(3)       :: kf          ! Scattered vector (in lab frame)
      real(kind=cp),dimension(3,3)     :: Tmatrix     ! Transformation matrix (recipr -> lab frame)
   end type MagH_Type_XRMS

   !
   !  Type :: MAGH_LIST_TYPE_XRMS
   !
   type, public :: MagH_List_Type_XRMS
      integer                                            :: NRef
      type(magh_type_xrms),dimension(:),allocatable      :: MH
   end type MagH_List_Type_XRMS
   
   contains
   
end module XRMS_types
