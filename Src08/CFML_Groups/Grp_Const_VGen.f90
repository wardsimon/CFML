!!----
!!----
!!----
SubModule (CFML_Groups) CFML_GRP_009
   Contains
   
   !!----
   !!---- GROUP_CONSTRUCTOR_GEN
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Group_Constructor_GenV(GenV, Grp, StrCode)
      !---- Arguments ----!
      character(len=*),dimension(:),intent(in)     :: GenV
      class(Spg_Type),              intent(in out) :: Grp
      character(len=*),optional,    intent(in)     :: StrCode
       
      !--- Local variables ---!
      character(len=40),    dimension(:),  allocatable :: gen
      type(Symm_Oper_Type), dimension(:),  allocatable :: Op
      type(rational),       dimension(:),  allocatable :: centre_coord,anticentre_coord
      type(rational),       dimension(:,:),allocatable :: Lat_tr, aLat_tr
      type(rational),       dimension(:,:),allocatable :: Mat
       
      integer :: d,i,ngen,invt,multip,centred,anticentred,Numops,num_lat,num_alat,mag_type

      !> Init
      call Clear_Error()
      
      !> Initializes Grp
      call Initialize_Group(Grp)
      
      call Check_Generators(GenV,gen)
      if (Err_CFML%Ierr /= 0) return
      
      d=get_dimension_gener(gen(1))
      ngen = size(gen)
      do i=1,ngen
         Grp%generators_list=trim(Grp%generators_list)//trim(gen(i))//";"
      end do
      Grp%generators_list=Grp%generators_list(1:len_trim(Grp%generators_list)-1)
      
      allocate(Op(maxnum_op))
      do i=1,maxnum_op
         call Allocate_Operator(d,Op(i))
      end do
      allocate(Mat(d,d))

      !>Construct the list of the generators on top of Op. The identity is always the first operator
      do i=1,ngen
         call Get_Mat_From_Symb(gen(i),Mat,invt)
         if (Err_CFML%Ierr /= 0) return
         Op(i+1)%Mat=Mat
         Op(i+1)%time_inv=invt
      end do
      ngen=ngen+1
      
      !> Construct the raw Group
      call Get_Group_From_Gener(ngen,Op,multip)
      if (Err_CFML%Ierr /= 0) return

      ! Allocate provisionally to Multip the lattice translations and anti-Translations
      allocate(Lat_tr(d-1,multip), aLat_tr(d-1,multip))
      allocate(centre_coord(d-1),anticentre_coord(d-1))

      call Reorder_Operators(multip, Op, centred, centre_coord, anticentred, anticentre_coord, &
                             Numops, num_lat, num_alat, Lat_tr, aLat_tr, mag_type)
      if (Err_CFML%Ierr /= 0) return

      Grp%multip=multip
      Grp%d=d
      if (allocated(Grp%Op)) deallocate(Grp%Op)
      call Allocate_Operators(d,multip,Grp%Op)
      Grp%Op(1:multip)=Op(1:multip)

      if (allocated(Grp%Symb_Op)) Deallocate(Grp%Symb_Op)
      allocate(Grp%Symb_Op(multip))
      do i=1,multip
         Grp%Symb_Op(i)=trim(Get_Symb_from_Oper(Op(i)))
      end do

      if (num_lat > 0) then
         if (allocated(Grp%Lat_tr)) Deallocate(Grp%Lat_tr)
         allocate(Grp%Lat_tr(1:d-1,1:Num_Lat))
      end if
      
      if (num_alat > 0) then
         if (allocated(Grp%aLat_tr)) Deallocate(Grp%aLat_tr)
         allocate(Grp%aLat_tr(1:d-1,1:Num_aLat))
      end if
      
      if (allocated(Grp%centre_coord)) Deallocate(Grp%centre_coord)
      if (allocated(Grp%anticentre_coord)) Deallocate(Grp%anticentre_coord)
      allocate(Grp%centre_coord(1:d-1))
      allocate(Grp%anticentre_coord(1:d-1))
      Grp%Numops           = Numops
      Grp%centred          = centred
      Grp%anticentred      = anticentred
      Grp%mag_type         = mag_type
      Grp%num_lat          = num_lat
      Grp%num_alat         = num_alat
      Grp%centre_coord     = centre_coord
      Grp%anticentre_coord = anticentre_coord
      if (num_lat  > 0)  Grp%Lat_tr = Lat_tr(1:d-1,1:Num_Lat)
      if (num_alat > 0)  Grp%aLat_tr=aLat_tr(1:d-1,1:Num_aLat)
      
      if (present(StrCode)) then
         if (trim(StrCode) /= 'xyz') then
            do i=1,Grp%Multip
               Grp%Symb_Op(i)=Get_Symb_from_Oper(Grp%Op(i),StrCode)
            end do
         end if
      end if 
       
      return
   End Subroutine Group_Constructor_GenV
   
End SubModule CFML_GRP_009   