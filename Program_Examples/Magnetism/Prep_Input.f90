Module prep_input
 !---- Use Modules ----! 
 use CFML_GlobalDeps,             only: sp, cp, eps
 use CFML_Atom_TypeDef,           only: mAtom_list_Type, allocate_mAtom_list
 use CFML_String_Utilities,       only: l_case
 use CFML_Reflections_Utilities,  only: Hkl_s
 use CFML_IO_Formats,             only: file_list_type
 use CFML_Propagation_vectors,    only: K_Equiv_Minus_K
 use CFML_crystal_Metrics,        only: Crystal_Cell_Type
 use CFML_Keywords_Code_Parser
 use CFML_Magnetic_Symmetry 
 use CFML_Magnetic_Structure_Factors

 implicit none

   !---- List of public subroutines ----!
 public :: Readn_Set_sqMiV_Data, Readn_Set_Job, Calc_sqMiV_Data, Copy_mAtom_List, Copy_MagDom_List, &
           Write_ObsCalc, Mag_Dom, mA

 integer, save               :: lun=1,lan=2,lin=3
 integer, save               :: Nf2
 real, public                :: Scale

   !---- Known used types are decleared belove the new types
   
    !!----
    !!---- TYPE :: mOBSERVATION_TYPE
    !!--..
    !!---- Type, public :: mObservation_Type
    !!----    integer                 :: ncont ! number of contributing reflections
    !!-----   integer,dimension(3,12) :: Hd    ! Indices of contributing reflections
    !!-----   integer,dimension(  12) :: icod  ! Pointer to the number of domain
    !!-----   integer,dimension(  12) :: p     ! Pointer to the list of independent reflections
    !!----    real(kind=sp)           :: Gobs  ! Observed  Magnetic Intensity
    !!----    real(kind=sp)           :: SGobs ! Sigma of  Gobs
    !!---- End Type mObservation_Type
    !!----
    !!---- Created: January - 2012
    !!
    Type, public :: mObservation_Type
       integer                                  :: ncont ! number of contributing reflections
       real(kind=sp),allocatable,dimension(:,:) :: Hd    ! (3,ncont) Indices of contributing reflections
       integer,allocatable,dimension(:)         :: icod  ! (ncont)   Pointer to the number of domain
       integer,allocatable,dimension(:)         :: p     ! (ncont)   Pointer to the list of independent reflections
       real(kind=sp)                            :: Gobs  ! Observed  Magnetic Intensity
       real(kind=sp)                            :: sGobs ! Sigma of  Gobs
       real(kind=sp)                            :: wGobs ! Weight of Gobs, usually 1/sGobs^2
    End Type mObservation_Type

    !!----
    !!---- TYPE :: mOBSERVATION_LIST_TYPE
    !!--..
    !!---- Type, public :: mObservation_List_Type
    !!----    integer                                          :: Nobs  ! Number of Observations
    !!----    type(mobservation_type),allocatable,dimension(:) :: Ob    ! Observation (F2, etc...)
    !!---- End Type mObservation_List_Type
    !!----
    !!---- Created: January - 2012
    !!
    Type, public :: mObservation_List_Type
       integer                                          :: Nobs  ! Number of Observations
       type(mObservation_type),allocatable,dimension(:) :: Ob    ! Observation (F2, etc...)
    End Type mObservation_List_Type

 type (file_list_type)             :: file_cfl
 type (MagSymm_k_Type)             :: MGp
 type (MagHD_Type)                 :: Mh
 type (MagHD_List_Type)            :: Mhlist
 type (mObservation_Type)          :: Ob
 type (mObservation_List_Type)     :: Oblist
 type (Crystal_Cell_Type)          :: Cell
 type (mAtom_List_Type),save       :: mA,mA_clone
 type (Magnetic_Domain_type),save  :: Mag_Dom,Mag_dom_clone
 
contains

!---- Subroutines ----!
!******************************************!
    Subroutine Readn_Set_sqMiV_Data(file_cfl)
!******************************************!
!    reads integrated intensities for the list of reflections given in *.int
       !---- Arguments ----!
       type(file_list_type),                intent (in)     :: file_cfl

       !---- Local variables ----!
       character(len=256)          :: title,formt
       character(len=256)          :: line,fildat
       character(len=1)            :: sig
       integer                     :: ier,num_k,m,ih,ik,il,ncont,iv,i,j,ityp,ipow
       real,   dimension(3)        :: vk,v_k
       real                        :: Gobs,SGobs,lambda
       logical                     :: esta
    
          esta=.false.
          do i=1,file_cfl%nlines
            line=adjustl(l_case(file_cfl%line(i)))
            if(line(1:10) == "f2mag_data") then
              line=file_cfl%line(i)
              j=index(line,"!")
              if(j /= 0) line=trim(line(1:j-1))
              fildat= adjustl(line(12:))
              inquire(file=trim(fildat),exist=esta)
              exit
            end if
          end do

          if(.not. esta) then
            write(unit=*,fmt="(a)") " => No .int has been provided! "
          end if

     open(unit=lin,file=trim(fildat), status="old",action="read")
! temporarily programmed for the special case k,-k, i.e. m=1,2     
     read(unit=lin,fmt=*) title
     read(unit=lin,fmt=*) formt
     read(unit=lin,fmt=*) lambda,ityp,ipow
     read(unit=lin,fmt=*) num_k

     do i=1,num_k
      read(unit=lin,fmt=*) m,vk
     enddo

     Nf2=0
     do !just to get Nf2
      read(unit=lin,fmt=*,end=1) ih,ik,il,m,Gobs,SGobs,ncont
      Nf2=Nf2+1
     enddo

1    rewind(unit=lin)

     if(allocated(Oblist%Ob)) deallocate(Oblist%Ob)
     allocate(Oblist%Ob(Nf2))
     
     if(allocated(Mhlist%Mh)) deallocate(Mhlist%Mh)
     allocate(Mhlist%Mh(Nf2))
     Mhlist%Nref= Nf2

! temporarily programmed for twin coefficient ncont=1
     do i=1,Nf2
      allocate(Oblist%Ob(i)%Hd(1:3,1))
     enddo
     Oblist%Nobs= Nf2
     
     read(unit=lin,fmt=*) title
     read(unit=lin,fmt=*) formt
     read(unit=lin,fmt=*) lambda,ityp,ipow
     read(unit=lin,fmt=*) num_k

     do i=1,num_k
      if(i.eq.1) then
       read(unit=lin,fmt=*) m,vk
      else
       read(unit=lin,fmt=*) m,v_k
       if(maxval(abs(vk+v_k)).gt.eps) then
        write(unit=*,fmt="(a)") " => Only k and -k are handled presently! "
        stop
       endif
      endif
     enddo

     do i=1,Nf2
      read(unit=lin,fmt=*,end=2) ih,ik,il,m,Oblist%Ob(i)%Gobs,Oblist%Ob(i)%SGobs,Oblist%Ob(i)%ncont
      if(m==1) then
       Oblist%Ob(i)%Hd(:,1)=real((/ih,ik,il/)) + vk(:)
       j=1
      endif
      if(m==2) then
       Oblist%Ob(i)%Hd(:,1)=real((/ih,ik,il/)) + v_k(:)
       j=-1
      endif

      !construct partially the object Mh
         if(j== 1) sig='+'
         if(j==-1) sig='-'
         Mh%signp=real(-j)  ! sign "+" for H-k and "-" for H+k
         iv=abs(j)
         Mh%num_k=iv
         Mh%h= real((/ih,ik,il/)) - Mh%signp*MGp%kvec(:,iv)
         Mh%s = hkl_s(Mh%h,Cell)
         Mh%keqv_minus=K_Equiv_Minus_K(MGp%kvec(:,iv),MGp%latt)
         MhList%Mh(i)=Mh
     enddo
       Oblist%Ob(:)%wGobs=1.d0/max([(Oblist%Ob(i)%sGobs**2,i=1,Nf2)],eps)

2    close(lin)     

    End Subroutine Readn_Set_sqMiV_Data
    
!******************************************!
    Subroutine Readn_Set_Job(file_cfl)
!******************************************!
       !---- Arguments ----!
       Type(file_list_type),   intent( in)  :: file_cfl

       !---- Local variables ----!
       character(len=132)   :: line
       integer              :: i,j

       i=0
       do j=1,file_cfl%nlines
          line=adjustl(file_cfl%line(j))
          line=l_case(line)
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle
          i=index(line,"!")
          if( i /= 0) line=trim(line(1:i-1))

          select case (line(1:4))
          case ("opti")    !Optimization
              i=index(line,"f2mag")
              if(i /= 0) then
               call Calc_sqMiV_Data(lun)
              endif
          end select
       end do
    End Subroutine Readn_Set_Job 

!******************************************!
    Subroutine Calc_sqMiV_Data(lun)
!******************************************!
! gets magn str Factor and Miv from Calc_Magnetic_StrF_MiV
! 
       !---- Argument ----!
       integer, optional,           intent(in) :: lun

       !---- Local variables ----!
       integer                     :: j,n
       real                        :: cost
       character(len=1)            :: mode

       if (present(lun)) write(lun,'(a)') '   hkl                      Imag   sImag   Icalc'
       n=Oblist%Nobs       
       do j=1,Nf2 ! This is loop over reflections with Int

         !Calculate magnetic structure factor and magnetic interaction vector
         ! as mode='y' MiV w.r.t. cartesian crystallographic frame
         mode='y'
         Mh=MhList%Mh(j)
!         call Calc_Magnetic_StrF_MiV(Cell,MGp,mA,Mh) ! was withput domains 
         call Calc_Magnetic_StrF_MiV_Dom(Cell,MGp,mA,Mag_Dom,Mh) !CFML_Msfac

!----    without mode components with respect to direct cell system {e1,e2,e3}
!----    with mode components with respect to the cartesian frame
         MhList%Mh(j)%sqMiV=Mh%sqMiV
       end do !end loop over reflections

!---     Cost_sqMiV(cost,Scale)
       if (present(lun)) then
          Scale=sum( [(MhList%Mh(j)%sqMiV*Oblist%Ob(j)%Gobs *Oblist%Ob(j)%wGobs,j=1,n)])/ &
                sum( [(MhList%Mh(j)%sqMiV**2 *Oblist%Ob(j)%wGobs,j=1,n)] )

         do j=1,Nf2 ! This is loop over reflections with Int
          write(unit=lun,fmt="(6f8.4)") Oblist%Ob(j)%Hd(:,1),&
                                        Oblist%Ob(j)%Gobs,Oblist%Ob(j)%SGobs,Scale*MhList%Mh(j)%sqMiV
         end do !end loop over reflections
         cost=sum( ([(Oblist%Ob(j)%wGobs* (Oblist%Ob(j)%Gobs-Scale*MhList%Mh(j)%sqMiV)**2,j=1,n)]) )/(n-NP_Refi)

         write(unit=lun,fmt="(/,a)")    "=================================="
         write(unit=lun,fmt="(a,f12.0)") 'Initial Cost(F2mag): Sum(|Obs-Scale*Calc|^2) / Sum(sigma^2) ) / (Nobs-Nref) =',cost
         write(unit=lun,fmt="(a,/)")    "=================================="

       end if

    End Subroutine Calc_sqMiV_Data
   
!******************************************!
    !!----
    !!---- Subroutine Copy_mAtom_List(mA, mA_clone)
    !!----    type(mAtom_list_type),   intent(in)  :: mA      !  In -> Atom list
    !!----    type(mAtom_list_type),   intent(out) :: mA_clone!   Out -> Atom list
    !!----
    !!----    Subroutine to copy an atom list to another one
    !!----
    !!---- Created: January - 2012
    !!
!******************************************!
    Subroutine Copy_mAtom_List(mA, mA_clone)
!******************************************!
       !---- Arguments ----!
       type(mAtom_list_type),   intent(in)  :: mA  
       type(mAtom_list_type),   intent(out) :: mA_clone
       !---- Local variables ----!
       integer                    :: n

       n=mA%natoms
       call Allocate_mAtom_List(n,mA_clone)
       mA_clone%atom(1:n)=mA%atom(1:n)
       return
    End Subroutine Copy_mAtom_List

!******************************************!
    !!----
    !!---- Subroutine Copy_MagDom_List(Mag_dom, Mag_dom_clone)
    !!----    type(mAtom_list_type),   intent(in)  :: Mag_dom      !  In -> Atom list
    !!----    type(mAtom_list_type),   intent(out) :: Mag_dom_clone!   Out -> Atom list
    !!----
    !!----    Subroutine to copy an mag domain list to another one
    !!----
    !!---- Created: January - 2012
    !!
!******************************************!
    Subroutine Copy_MagDom_List(Mag_dom, Mag_dom_clone)
!******************************************!
       !---- Arguments ----!
       type(Magnetic_Domain_type),   intent(in)  :: Mag_dom  
       type(Magnetic_Domain_type),   intent(out) :: Mag_dom_clone
       !---- Local Variables ----!

       
       Mag_dom_clone%nd=Mag_dom%nd
       Mag_dom_clone%Chir=Mag_dom%Chir
       Mag_dom_clone%DMat=Mag_dom%DMat
       Mag_dom_clone%Pop=Mag_dom%Pop
       Mag_dom_clone%Mpop=Mag_dom%Mpop
       Mag_dom_clone%Lpop=Mag_dom%Lpop
       Mag_dom_clone%Lab=Mag_dom%Lab

       return
    End Subroutine Copy_MagDom_List

!******************************************!
    Subroutine Write_ObsCalc(file_cfl)
!******************************************!

       !---- Arguments ----!
       type(file_list_type),                intent (in)     :: file_cfl
       !---- Local variables ----!
       character(len=132)   :: line
       integer              :: i,j

        i=0
        do j=1,file_cfl%nlines
          line=adjustl(file_cfl%line(j))
          line=l_case(line)
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle
          i=index(line,"!")

          if( i /= 0) line=trim(line(1:i-1))

          select case (line(1:4))
          case ("opti")    !Optimization

          i=index(line,"f2mag")
           if(i /= 0) then
            write(lun,'(a)') '     hkl                      Iobs      Ical'
            do i=1,Nf2
             write(unit=lun,fmt='(3(f8.4),2x,2(f12.2))') Oblist%Ob(i)%Hd(:,1),Oblist%Ob(i)%Gobs,Scale*Mhlist%Mh(i)%sqMiV
            enddo
           endif

          end select
        enddo
   
   End Subroutine Write_ObsCalc
   
end module prep_input