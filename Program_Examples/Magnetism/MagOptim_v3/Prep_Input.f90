 Module prep_input
 !---- Use Modules ----!
 use CFML_GlobalDeps,             only: sp, cp, eps
 use CFML_Atom_TypeDef,           only: mAtom_list_Type, allocate_mAtom_list
 use CFML_String_Utilities,       only: L_Case, U_Case, Getword
 use CFML_Reflections_Utilities,  only: Hkl_s
 use CFML_IO_Formats,             only: file_list_type
 use CFML_Propagation_vectors,    only: K_Equiv_Minus_K
 use CFML_crystal_Metrics,        only: Crystal_Cell_Type
 use CFML_Keywords_Code_Parser
 use CFML_Magnetic_Symmetry
 use CFML_Magnetic_Structure_Factors
 use CFML_Polarimetry

 implicit none

   !---- List of public subroutines ----!
 public :: Readn_Set_Data, Readn_Set_Job, Copy_mAtom_List, Copy_MagDom_List, Write_ObsCalc, MagDom_to_Dataset

 integer, save   :: lun=1,lan=2,lin=3
 integer, save   :: Nobs,Nset,Nref,Nf2,Npol
 real, public    :: Scalef,cost,costPol,costF2

   !---- OLD types are declared below the NEW types
   !---- NEW types
   !!----    mOBSERVATION_TYPE
   !!----    mOBSERVATION_LIST_TYPE
   !!----    Multidata_type  !NEW 02/12

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
       real(kind=sp)                            :: sGobs ! Sigma of  Gobs,
       real(kind=sp)                            :: wGobs ! Weight of Gobs, 1/sGobs^2
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

    !!
    !!---- TYPE :: Multidata_type
    !!---- integer                                              :: Nset      ! Number of Datasets
    !!---- integer,allocatable, dimension(:)                    :: Nobs      ! Number of Observations in the Dataset
    !!---- character(len=128),allocatable, dimension(:)         :: Datfile   ! Name of the File with the Dataset
    !!---- logical,allocatable, dimension(:)                    :: SNP       !True if the polarisation Dataset
    !!---- type(Magnetic_Domain_type),allocatable, dimension(:) :: MagDom    ! Domains for the Dataset
    !!
    !!----
    !!---- Created: February - 2012
    !!
    Type, public :: Multidata_type
       integer                                              :: Nset      ! Number of Datasets
       integer,allocatable, dimension(:)                    :: Nobs      ! Number of Observations in the Dataset
       character(len=128),allocatable, dimension(:)         :: Datfile   ! Name of the File with the Dataset
       logical,allocatable, dimension(:)                    :: SNP       ! True if the polarisation Dataset Blume-Maleyev type
       type(Magnetic_Domain_type),allocatable, dimension(:) :: MagDom    ! Domains for the Dataset
    End Type Multidata_type

    !!----
    !!----  MAGHD_MULTILIST_TYPE
    !!----
    Type, Public  :: MagHD_MultiList_Type
       integer                                   :: Nset
       Type(MagHD_List_Type),allocatable, dimension(:) :: Mhlist
    End Type MagHD_MultiList_Type
    !!----
    !!---- Created: February - 2012 OZ

 type (file_list_type)                  :: file_cfl
 type (MagSymm_k_Type)                  :: MGp
 type (MagHD_Type)                      :: Mh
 type (MagHD_List_Type)                 :: Mhlist
 type (mObservation_Type)               :: Ob
 type (mObservation_List_Type)          :: Oblist
 type (Crystal_Cell_Type)               :: Cell
 type (mAtom_List_Type),public          :: mA,mA_clone
 type (Magnetic_Domain_type),public,save:: Mag_Dom,AllMag_Dom,Mag_dom_clone
 type (Multidata_type),save             :: Multidata
 type (MagHD_MultiList_Type)            :: MhMultilist
 type (Polar_calc_type)                 :: Polari
 type (Polar_Calc_List_type)            :: Polarilist
 type (Polar_CalcMulti_List_type),public  :: PolariMultilist
 type (Polar_calc_sVs_type)             :: PolarisVs
 type (Polar_Calc_sVs_List_type)        :: PolarisVslist
 type (Polar_CalcMulti_sVs_List_type),public:: PolariMultisVslist
 type (Polar_obs_type)                  :: Polaro
 type (Polar_Obs_List_type)             :: Polarolist
 type (Polar_ObsMulti_List_type),public    :: PolaroMultilist

!---- List of public arrays
! public :: Mag_Dom, mA, PolaroMultilist, PolariMultilist, PolariMultisVslist


contains

!---- Subroutines ----!
!******************************************!
    Subroutine Readn_Set_Data(file_cfl)
!******************************************!
!     Analyses and Reads observations given in different files
!     Nobs - total number of reflections

      !---- Arguments ----!
      type(file_list_type),   intent(in) :: file_cfl

      !---- Local variables ----!
      character(len=50), dimension(25):: label
      character(len=256)              :: title,formt
      character(len=256)              :: lowline
      character(len=1)                :: sig
      integer                         :: ier,ic,num_k,m,ih,ik,il,ncont,iv,i,j,ityp,ipow,ip,iset,&
                                         iobs,ind,nd
      integer                         :: n_ini, n_end
      real,   dimension(3)            :: vk,v_k
      real                            :: Gobs,SGobs,lambda
      logical                         :: kvect_begin=.true.

      Nset=0
      Npol=0

      !---
      !--- Read the wave vector as in Readn_Set_Magnetic_Structure
      !---
      do i=1,file_cfl%nlines
        lowline=l_case(adjustl(file_cfl%line(i)))
        if (lowline(1:13) == "mag_structure") then
           n_ini=i
        end if
        if (lowline(1:7) =="end_mag" ) then
           n_end=i
           exit
        end if
      end do

      num_k=0
      i=n_ini
      do
        i=i+1
        if(i >= n_end) exit
        lowline=l_case(adjustl(file_cfl%line(i)))
        ! Read propagation vectors
        if (lowline(1:5) == "kvect" .and. kvect_begin) then
          num_k=num_k+1
          read(unit=lowline(6:),fmt=*,iostat=ier) MGp%kvec(:,num_k)
          if (ier /= 0) then
             err_magsym=.true.
             ERR_MagSym_Mess=" Error reading propagation vectors"
             return
          end if
          do !repeat reading until continuous KVECT lines are exhausted
             i=i+1
             lowline=adjustl(l_case(file_cfl%line(i)))
             if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
             if (lowline(1:5) == "kvect") then
                num_k=num_k+1
                read(unit=lowline(6:),fmt=*,iostat=ier) MGp%kvec(:,num_k)
                if (ier /= 0) then
                   err_magsym=.true.
                   ERR_MagSym_Mess=" Error reading propagation vectors"
                   return
                end if
             else
                i=i-1
                kvect_begin=.false.
                exit
             end if
          end do
          cycle
        end if
      end do
      MGp%nkv   =num_k

      !---
      !--- Analyse how many data are provided in cfl file
      !---
      do i=1,file_cfl%nlines
        lowline=l_case(adjustl(file_cfl%line(i)))
        if(lowline(1:10) == "f2mag_data") then
         Nset=Nset+1
        end if

        if(lowline(1:12) == "cryopad_data") then
         Nset=Nset+1
         Npol=Npol+1
        endif

        if(lowline(1:10) == "mupad_data") then
         Nset=Nset+1
         Npol=Npol+1
        endif
      end do

      if(allocated(Multidata%datfile)) deallocate(Multidata%datfile)
      allocate(Multidata%datfile(Nset))
      Multidata%Nset= Nset

      if(allocated(Multidata%MagDom)) deallocate(Multidata%MagDom)
      allocate(Multidata%MagDom(Nset))
      if(allocated(Multidata%Nobs)) deallocate(Multidata%Nobs)
      allocate(Multidata%Nobs(Nset))
      if(allocated(Multidata%SNP)) deallocate(Multidata%SNP)
      allocate(Multidata%SNP(Nset))

      if(allocated(PolaroMultilist%Polarolist)) deallocate(PolaroMultilist%Polarolist)
      allocate(PolaroMultilist%Polarolist(Npol))

      if(allocated(MhMultilist%Mhlist)) deallocate(MhMultilist%Mhlist)
      allocate(MhMultilist%Mhlist(Nset))
      MhMultilist%Nset= Nset

      Nset=0
      !---
      !--- Write Multidata%datfile,Multidata%MagDom,etc
      !---
      do i=1,file_cfl%nlines
        lowline=l_case(adjustl(file_cfl%line(i)))

        if(lowline(1:10) == "f2mag_data") then
          Nset=Nset+1
          lowline=adjustl(lowline(12:))
          call getword(lowline,label,ic)
          Multidata%datfile(Nset)=l_case(label(1))
          nd=ic-1 !number of domains read from f2mag_data line

          !--- Here I assume that two chiral domains are given in Mag_data
          if(nd >= 2 .and. modulo(nd,2) == 0) then
            Multidata%MagDom(Nset)%nd=nd/2
            ind=1
            do j=2,nd,2
              Multidata%MagDom(Nset)%Lab(1,ind)=l_case(label(j))
              Multidata%MagDom(Nset)%Lab(2,ind)=l_case(label(j+1))
              ind=ind+1
            end do
          else
            write(unit=*,fmt="(a)") " => All domains should be given in magdom and _data! "
            stop
          end if

          call Countref_Int_file(Multidata%datfile(Nset),Nf2)
          Multidata%Nset=Nset
          Multidata%Nobs(Nset)=Nf2
          Multidata%SNP(Nset)=.false.
        end if

        if(lowline(1:12) == "cryopad_data") then
          Nset=Nset+1
          lowline=adjustl(lowline(14:))
          call getword(lowline,label,ic)
          Multidata%datfile(Nset)=label(1)
          nd=ic-1 !number of domains read from cryopad_data line

          !--- Here I assume that two chiral domains are given in Mag_data
          if(nd >= 2 .and. modulo(nd,2) ==0) then
            Multidata%MagDom(Nset)%nd=nd/2
            ind=1
            do j=2,nd,2
              Multidata%MagDom(Nset)%Lab(1,ind)=l_case(label(j))
              Multidata%MagDom(Nset)%Lab(2,ind)=l_case(label(j+1))
              ind=ind+1
            end do
          else
           write(unit=*,fmt="(a)") " => All domains should be given in magdom and _data! "
           stop
          end if

          call Countref_Pol_file(Multidata%datfile(Nset),Npol)
          Multidata%Nset=Nset
          Multidata%Nobs(Nset)=Npol
          Multidata%SNP(Nset)=.true.
        end if

        if(lowline(1:10) == "mupad_data") then
          Nset=Nset+1
          lowline=adjustl(lowline(12:))
          call getword(lowline,label,ic)
          Multidata%datfile(Nset)=label(1)
          nd=ic-1 !number of domains read from mupad_data line

          !--- Here I assume that two chiral domains are given in Mag_data
          if(nd >= 2 .and. modulo(nd,2) == 0) then
            Multidata%MagDom(Nset)%nd=nd/2
            ind=1
            do j=2,nd,2
              Multidata%MagDom(Nset)%Lab(1,ind)=l_case(label(j))
              Multidata%MagDom(Nset)%Lab(2,ind)=l_case(label(j+1))
              ind=ind+1
            end do
          else
            write(unit=*,fmt="(a)") " => All domains should be given in magdom and _data! "
            stop
          end if
          call Countref_Pol_file(Multidata%datfile(Nset),Npol)
          Multidata%Nset=Nset
          Multidata%Nobs(Nset)=Npol
          Multidata%SNP(Nset)=.true.
        end if

      end do

      do iset=1,Nset !Loop on all data

        if(Multidata%SNP(iset)) then
          !---
          !--- Polarisation data
          !---

          PolaroMultilist%Nset= iset
          open(unit=lan,file=trim(Multidata%datfile(iset)),status="old",action="read")

          Npol=Multidata%Nobs(iset)
          if(allocated(PolaroMultilist%Polarolist(iset)%Polaro)) deallocate(PolaroMultilist%Polarolist(iset)%Polaro)
          allocate(PolaroMultilist%Polarolist(iset)%Polaro(Npol))
          PolaroMultilist%Polarolist(iset)%Nref= Npol

          if(allocated(MhMultilist%Mhlist(iset)%Mh)) deallocate(MhMultilist%Mhlist(iset)%Mh)
          allocate(MhMultilist%Mhlist(iset)%Mh(Npol))
          MhMultilist%Mhlist(iset)%Nref= Npol

          do iobs=1,Npol !Loop on polarised observations
            read(unit=lan,fmt=*,iostat=ier) ih,ik,il,m,PolaroMultilist%Polarolist(iset)%Polaro(iobs)%SPV, &
                     PolaroMultilist%Polarolist(iset)%Polaro(iobs)%P  !2 vectors in scattering plane, Pin
            if(ier /= 0) exit
            read(unit=lan,fmt=*) PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij
            read(unit=lan,fmt=*) PolaroMultilist%Polarolist(iset)%Polaro(iobs)%soPij

            ! calculate weight Polarolist%Polaro(i)%woPij 1/sigma**2
            PolaroMultilist%Polarolist(iset)%Polaro(iobs)%woPij= 1.d0 / max(eps, &
            PolaroMultilist%Polarolist(iset)%Polaro(iobs)%soPij**2)

            !construct partially the object Mh
            j=sign(1,m)
            sig="+"
            if( j == -1) sig="-"
            Mh%signp=real(-j)  ! sign "+" for H-k and "-" for H+k
            iv=abs(m)

            Mh%num_k=iv
            Mh%h= real((/ih,ik,il/)) - Mh%signp*MGp%kvec(:,iv)
            Mh%s = hkl_s(Mh%h,Cell)
            Mh%keqv_minus=K_Equiv_Minus_K(MGp%kvec(:,iv),MGp%latt)
            MhMultilist%Mhlist(iset)%Mh(iobs)=Mh

            PolaroMultilist%Polarolist(iset)%Polaro(iobs)%H=MhMultilist%Mhlist(iset)%Mh(iobs)%h
          end do !Loop on observations

          close(lan)

        else

          !---
          ! Integrated intensity data
          !---

          open(unit=lin,file=trim(Multidata%datfile(iset)),status="old",action="read")
          Nf2=Multidata%Nobs(iset)

           ! temporarily programmed for the special case k,-k
          read(unit=lin,fmt=*) title
          read(unit=lin,fmt=*) formt
          read(unit=lin,fmt=*) lambda,ityp,ipow
          read(unit=lin,fmt=*) num_k

          do i=1,num_k
            if(i.eq.1) then
              read(unit=lin,fmt=*) m,vk
            else
              read(unit=lin,fmt=*) m,v_k
            end if
          end do

          if(allocated(Oblist%Ob)) deallocate(Oblist%Ob)
          allocate(Oblist%Ob(Nf2))

          ! temporarily programmed for ncont=1
          do iobs=1,Nf2
            allocate(Oblist%Ob(iobs)%Hd(1:3,1))
          end do
          Oblist%Nobs= Nf2

          !Multidata%Nobs(iset) should be the same as Oblist%Nobs
          if(allocated(MhMultilist%Mhlist(iset)%Mh)) deallocate(MhMultilist%Mhlist(iset)%Mh)
          allocate(MhMultilist%Mhlist(iset)%Mh(Nf2))
          MhMultilist%Mhlist(iset)%Nref= Nf2

          do iobs=1,Nf2
            read(unit=lin,fmt=*,end=3) ih,ik,il,m,Oblist%Ob(iobs)%Gobs,Oblist%Ob(iobs)%SGobs,ncont

            if(m == 1) then
              Oblist%Ob(iobs)%Hd(:,1)=real((/ih,ik,il/)) + vk(:)
              j=1; sig="+"
            endif
            if(m == 2) then
              Oblist%Ob(iobs)%Hd(:,1)=real((/ih,ik,il/)) + v_k(:)
              j=-1; sig="-"
            endif
            !construct partially the object Mh
            Mh%signp=real(-j)  ! sign "+" for H-k and "-" for H+k
            iv=abs(j)
            MGp%kvec(:,iv)=vk(:)
            Mh%num_k=iv
            Mh%h= real((/ih,ik,il/)) - Mh%signp*MGp%kvec(:,iv)
            Mh%s = hkl_s(Mh%h,Cell)
            Mh%keqv_minus=K_Equiv_Minus_K(MGp%kvec(:,iv),MGp%latt)
            MhMultilist%Mhlist(iset)%Mh(iobs)=Mh

          end do

!         calculate weight Oblist%Ob(i)%wGobs 1/sGobs
          Oblist%Ob(:)%wGobs=1.d0/max([(Oblist%Ob(i)%sGobs**2,i=1,Oblist%Nobs)],eps)

           3    close(lin)
        end if !if .not. SNP(iset)
      end do !Loop on all data

      Nobs=sum(MultiData%Nobs(1:Nset))

    End Subroutine Readn_Set_Data

!******************************************!
    Subroutine Countref_Pol_file(datfile,Npol)
!******************************************!
! reads file with polarization data (counts entries)
!
       !---- Argument ----!
      character(len=*),intent(in) :: datfile
      integer,intent(out)         :: Npol

       !---- Local variables ----!
      integer                     :: ih,ik,il,m
      logical             :: pesta, argtaken=.false.
      real                :: vpl, pol


     inquire(file=trim(datfile),exist=pesta)
     if( .not. pesta) then
       write(unit=*,fmt="(a)") " File: "//trim(datfile)//" doesn't exist!"
       stop
     end if
     open(unit=lan,file=trim(datfile), status="old",action="read")
     Npol=0

     do !just to get Npol
      read(unit=lan,fmt=*,end=1) ih,ik,il,m,vpl,pol  !2 vectors in scattering plane, Pin
      read(unit=lan,fmt=*) Polaro%oPij
      read(unit=lan,fmt=*) Polaro%soPij
      Npol=Npol+1
     enddo

1    close(lan)

    End Subroutine Countref_Pol_File

!******************************************!
    Subroutine Countref_Int_File(datfile,Nf2)
!******************************************!
! reads file with intensity data (counts entries)
!
       !---- Argument ----!
      character(len=*),intent(in) :: datfile
      integer,intent(out)         :: Nf2

       !---- Local variables ----!
      character(len=256)          :: title,formt,lambda
      character(len=1)            :: sig
      integer                     :: ier,num_k,m,ih,ik,il,ncont,i,ityp,ipow
      real,   dimension(3)        :: vk,v_k
      real                        :: Gobs,SGobs
      logical                     :: pesta, argtaken=.false.

     inquire(file=trim(datfile),exist=pesta)
     if( .not. pesta) then
       write(unit=*,fmt="(a)") " File: "//trim(datfile)//" does'nt exist!"
       stop
     end if
     open(unit=lin,file=trim(datfile), status="old",action="read")

! temporarily programmed for the special case k,-k, i.e. m=1,2
     read(unit=lin,fmt=*) title
     read(unit=lin,fmt=*) formt
     read(unit=lin,fmt=*) lambda,ityp,ipow
     read(unit=lin,fmt=*) num_k

     if(num_k /= MGp%nkv*2) then
        write(unit=*,fmt="(a)") " => Mismatch between number of wavevectors in .cfl and .int files! "
        stop
     end if


     do i=1,num_k
      if(i.eq.1) then
       read(unit=lin,fmt=*) m,vk
      else
       read(unit=lin,fmt=*) m,v_k
       if(maxval(abs(vk+v_k)).gt.eps) then
        write(unit=*,fmt="(a)") " => Only k and -k are handled presently! "
        stop
       end if
      end if
     end do

     Nf2=0
     do !just to get Nf2
      read(unit=lin,fmt=*,end=1) ih,ik,il,m,Gobs,SGobs,ncont
      Nf2=Nf2+1
     enddo

1    close(lin)

    End Subroutine Countref_Int_File

!******************************************!
    Subroutine Readn_Set_Job(file_cfl)
!******************************************!
       !---- Arguments ----!
       Type(file_list_type),   intent( in)  :: file_cfl

       !---- Local variables ----!
       character(len=132)   :: line
       integer              :: i,j,iset,i1,i2,i3,i4

       i=0; i1=0; i2=0; i3=0; i4=0

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

              i3=index(line,"f2mag_cryopad")
               if(i3 /= 0) then
                do iset=1,Nset
                 if(Multidata%SNP(iset)) then
                  if(allocated(PolariMultilist%Polarilist)) deallocate(PolariMultilist%Polarilist)
                  allocate(PolariMultilist%Polarilist(iset))
                  call Calc_Polar_Dom_Data(iset,lun)
                 else
                  call Calc_sqMiV_Data(iset,lun)
                 end if
                end do
               end if

              i4=index(line,"f2mag_mupad")
               if(i4 /= 0) then
                do iset=1,Nset
                  if(Multidata%SNP(iset)) then
                    if(allocated(PolariMultisVslist%PolarisVslist)) deallocate(PolariMultisVslist%PolarisVslist)
                    allocate(PolariMultisVslist%PolarisVslist(iset))
                    call Calc_Polar_CrSec_Data(iset,lun)
                  else
                    call Calc_sqMiV_Data(iset,lun)
                  end if
                end do
               end if

               if(i3 == 0 .and. i4 == 0) then
                call Calc_sqMiV_Data(iset,lun)
               end if

             else

              i1=index(line,"cryopad")
              if(i1 /= 0) then
               do iset=1,Nset
                 if(Multidata%SNP(iset)) then
                   if(allocated(PolariMultilist%Polarilist)) deallocate(PolariMultilist%Polarilist)
                   allocate(PolariMultilist%Polarilist(iset))
                   call Calc_Polar_Dom_Data(iset,lun)
                 end if
               end do
              end if

              i2=index(line,"mupad")
              if(i2 /= 0) then
               do iset=1,Nset
                 if(Multidata%SNP(iset)) then
                    if(allocated(PolariMultisVslist%PolarisVslist)) deallocate(PolariMultisVslist%PolarisVslist)
                    allocate(PolariMultisVslist%PolarisVslist(iset))
                    call Calc_Polar_CrSec_Data(iset,lun)
                 end if
               end do
              end if

             end if

          end select
       end do
    End Subroutine Readn_Set_Job

!******************************************!
    Subroutine Calc_sqMiV_Data(iset,lun)
!******************************************!
! gets magn str Factor and Miv from Calc_Magnetic_StrF_MiV
! lun is used only at 0-cycle !!
       !---- Argument ----!
       integer, optional,           intent(in) :: iset,lun

       !---- Local variables ----!
       integer                     :: j

       Nf2=MhMultilist%MhList(iset)%Nref
       Mag_Dom=Multidata%MagDom(iset)
       !write(*,*) " Mag_Dom - Occ:", Mag_Dom%MPop(1:2,1)
       do j=1,Nf2 ! This is loop over reflections with Int

         !Calculate magnetic structure factor and magnetic interaction vector
         ! as mode='Car' MiV w.r.t. cartesian crystallographic frame

         Mh=MhMultilist%MhList(iset)%Mh(j)
         call Calc_Magnetic_StrF_MiV_Dom(Cell,MGp,mA,Mag_Dom,Mh) !CFML_Msfac, sums domains

!----    without components with respect to direct cell system {e1,e2,e3}
!----    with components with respect to the cartesian frame
         !write(*,*) "  h k l   M^2 ",Mh%h,Mh%sqMiV
         MhMultilist%MhList(iset)%Mh(j)%sqMiV=Mh%sqMiV

       end do !end loop over reflections

!---     Cost_sqMiV(cost,Scalef)
       if (present(lun)) then
          Scalef=sum( [(MhMultilist%MhList(iset)%Mh(j)%sqMiV*Oblist%Ob(j)%Gobs *Oblist%Ob(j)%wGobs,j=1,Nf2)])/ &
                 sum( [(MhMultilist%MhList(iset)%Mh(j)%sqMiV**2 *Oblist%Ob(j)%wGobs,j=1,Nf2)] )
         cost=sum( ([(Oblist%Ob(j)%wGobs* (Oblist%Ob(j)%Gobs-Scalef*MhMultilist%MhList(iset)%Mh(j)%sqMiV)**2, &
                                                                                 j=1,Nf2)]) ) /Nf2 !-NP_Refi)

         write(unit=lun,fmt="(/,a)")    "================================================================================"
         write(unit=lun,fmt="(a,f12.0)") 'Initial partial Cost(F2mag): Sum({|IObs-Scale*ICalc|/Sigma}^2)) / Nobs =',cost
         write(unit=lun,fmt="(a,/)")    "================================================================================"
         write(lun,'(a,i2,1x,a)') 'iset= ',iset, MultiData%datfile(iset)
         write(lun,'(a)') '     h      k       l          Imag        sImag      Icalc'

         do j=1,Nf2 ! This is loop over reflections with Int
          write(unit=lun,fmt="(3f8.4, 3f12.4)") Oblist%Ob(j)%Hd(:,1),&
                                        Oblist%Ob(j)%Gobs,Oblist%Ob(j)%SGobs,MhMultilist%MhList(iset)%Mh(j)%sqMiV
!          write(unit=lun,fmt="(3f8.4, 3f12.4)") MhMultilist%MhList(iset)%Mh(j)%H,&
!                                        Oblist%Ob(j)%Gobs,Oblist%Ob(j)%SGobs,MhMultilist%MhList(iset)%Mh(j)%sqMiV
         end do !end loop over reflections

       end if

    End Subroutine Calc_sqMiV_Data

!******************************************!
    Subroutine Calc_Polar_Dom_Data(iset,lun)
!******************************************!
    ! gets polarization matrices from Calc_Polar_Dom (via nc, mxyz,tc), domains are considered in Calc_Polar_Dom
    ! lun is used only at 0-cycle !!
      !---- Argument ----!
      integer, optional, intent(in) :: iset,lun

      !---- Local variables ----!
      integer                     :: j,ich,nch,nd
      real(kind=cp)               :: Pin
      real(kind=cp),dimension(3)  :: SPV

      Npol=Multidata%Nobs(iset)
      Mag_dom=Multidata%MagDom(iset)

      if(allocated(PolariMultilist%Polarilist(iset)%Polari)) deallocate(PolariMultilist%Polarilist(iset)%Polari)
      allocate(PolariMultilist%Polarilist(iset)%Polari(Npol))
      PolariMultilist%Polarilist(iset)%Nref= Npol

      do j=1,Npol ! Loop over hkl observations
         SPV=PolaroMultilist%Polarolist(iset)%Polaro(j)%SPV
         Pin=PolaroMultilist%Polarolist(iset)%Polaro(j)%P

         !Calculate magnetic structure factor and magnetic interaction vector
         ! as mode='Car' MiV w.r.t. cartesian crystallographic frame

         Mh=MhMultilist%MhList(iset)%Mh(j)
         call Calc_Magnetic_StrF_MiV_Dom(Cell,MGp,mA,Mag_Dom,Mh) !CFML_Msfac
         !call Calc_Polar_Dom(Cell,Mh%h,SPV,Pin,cmplx(0.0,0.0),Mag_Dom,Mh,PolariMultilist%Polarilist(iset)%Polari(j),B_Q="BM") !CFML_Polar
         call Calc_Polar_Dom(Cell,Mh%h,SPV,Pin,cmplx(0.0,0.0),Mag_Dom,Mh,PolariMultilist%Polarilist(iset)%Polari(j)) !CFML_Polar
         !without components with respect to direct cell system {e1,e2,e3}
         !with components with respect to the cartesian frame
         MhMultilist%MhList(iset)%Mh(j)%sqMiV=Mh%sqMiV

      end do !end loop over hkl observations

      !---     Cost_cryopad(cost)
      if (present(lun)) then

        cost=0.0
        do j=1,MultiData%Nobs(iset) !loop over observations
          cost =  cost + sum(PolaroMultilist%Polarolist(iset)%Polaro(j)%woPij * &
                 ( (PolariMultilist%Polarilist(iset)%Polari(j)%Pij - &
                    PolaroMultilist%Polarolist(iset)%Polaro(j)%oPij)**2))
        end do !end loop over observations

        cost=cost/(9*MultiData%Nobs(iset)) !-NP_Refi)

        write(unit=lun,fmt="(/,a)")    "==========================================================================="
        write(unit=lun,fmt="(a,f12.0)") 'Initial partial-Cost(cryopad): Sum({|PObs-PCalc|/sigma}^2) / Nobs =',cost
        write(unit=lun,fmt="(a,/)")    "==========================================================================="
        write(lun,'(a,i2,1x,a)') 'iset= ',iset, MultiData%datfile(iset)
        write(lun,'(a)') '    Pobs                   Pcalc'

        do j=1,MultiData%Nobs(iset) !loop over observations
          write(unit=lun,fmt='(3(f10.6))')    PolaroMultilist%Polarolist(iset)%Polaro(j)%H
          write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(j)%oPij(:,1), &
          PolariMultilist%Polarilist(iset)%Polari(j)%Pij(:,1)
          write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(j)%oPij(:,2), &
          PolariMultilist%Polarilist(iset)%Polari(j)%Pij(:,2)
          write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(j)%oPij(:,3), &
          PolariMultilist%Polarilist(iset)%Polari(j)%Pij(:,3)
        end do !end loop over observations

      end if

    End Subroutine Calc_Polar_Dom_Data

!******************************************!
    Subroutine Calc_Polar_CrSec_Data(iset,lun)
!******************************************!
      ! gets polarized cross sections from Calc_Polar_CrSec (via sVs)
      !
      !---- Argument ----!
      integer, optional,           intent(in) :: iset,lun

      !---- Local variables ----!
      integer                     :: i,j,iobs,nd,ich,nch
      real(kind=cp)               :: Pin
      real(kind=cp),dimension(3)  :: SPV
      real(kind=cp),dimension(3,3):: Ipp,Ipm,Imp,Imm

      Npol=Multidata%Nobs(iset)
      Mag_dom=Multidata%MagDom(iset)

      if(allocated(PolariMultisVslist%PolarisVslist(iset)%PolarisVs)) deallocate(PolariMultisVslist%PolarisVslist(iset)%PolarisVs)
      allocate(PolariMultisVslist%PolarisVslist(iset)%PolarisVs(Npol))
      PolariMultisVslist%PolarisVslist(iset)%Nref= Npol

      do iobs=1,Npol ! Loop over hkl observations
        SPV=PolaroMultilist%Polarolist(iset)%Polaro(iobs)%SPV
        Pin=PolaroMultilist%Polarolist(iset)%Polaro(iobs)%P
        nch=1
        if(Mag_Dom%chir) nch=2

        do nd=1,Mag_Dom%nd !loop over S-domains
          do ich=1,nch !loop over chiral domains
            !Calculate magnetic structure factor and magnetic interaction vector

            Mh=MhMultilist%MhList(iset)%Mh(iobs)
            call Calc_Magnetic_StrF_MiV_Dom(Cell,MGp,mA,Mag_Dom,Mh) !CFML_Msfac
            call Calc_Polar_CrSec(Cell,Mh%h,SPV,Pin,cmplx(0.0,0.0),Mag_Dom,Mh,Ipp,Ipm,Imp,Imm) !CFML_Polar
            !write(*,*) Mh%h
            !without components with respect to direct cell system {e1,e2,e3}
            !with components with respect to the cartesian frame
            MhMultilist%MhList(iset)%Mh(iobs)%sqMiV=Mh%sqMiV
          end do !loop over chiral domains
        end do !loop over S-domains

        if(Pin > 0.d0) then
          do i=1,3
            do j=1,3
              if((Ipp(i,j)+Imp(i,j)).le.eps) then
                PolariMultisVslist%PolarisVslist(iset)%PolarisVs(iobs)%Pij(i,j) = 0.d0
              else
                PolariMultisVslist%PolarisVslist(iset)%PolarisVs(iobs)%Pij(i,j) = & ! calculate using the same structure parameters AND
                            (Ipp(i,j)-Imp(i,j))/(Ipp(i,j)+Imp(i,j))    ! the set of populations for the iset-th dataset
              end if
            end do
          end do
        end if

        if(Pin < 0.d0) then
          do i=1,3
            do j=1,3
              if((Imm(i,j)+Ipm(i,j)).le.eps) then
                PolariMultisVslist%PolarisVslist(iset)%PolarisVs(iobs)%Pij(i,j) = 0.d0
              else
                PolariMultisVslist%PolarisVslist(iset)%PolarisVs(iobs)%Pij(i,j) = &
                              -(Imm(i,j)-Ipm(i,j))/(Imm(i,j)+Ipm(i,j))
              end if
            end do
          end do
        end if

      end do !end loop over hkl observations

      !---   Cost_mupad(cost)
      if (present(lun)) then

        cost=0.0
        do j=1,MultiData%Nobs(iset) !loop over observations
          cost =  cost + sum(PolaroMultilist%Polarolist(iset)%Polaro(j)%woPij * &
                  ( (PolariMultisVslist%PolarisVslist(iset)%PolarisVs(j)%Pij - &
                     PolaroMultilist%Polarolist(iset)%Polaro(j)%oPij)**2))
        end do !end loop over observations
        cost=cost/(9*MultiData%Nobs(iset)) !-NP_Refi)

        write(unit=lun,fmt="(/,a)")    "==========================================================================="
        write(unit=lun,fmt="(a,f12.0)") 'Initial partial-Cost(mupad): Sum({|PObs-PCalc|/sigma}^2) / Nobs =',cost
        write(unit=lun,fmt="(a,/)")    "==========================================================================="
        write(lun,'(a,i2,1x,a)') 'iset= ',iset, MultiData%datfile(iset)
        write(lun,'(a)') '    Pobs                   Pcalc'

        do j=1,MultiData%Nobs(iset) !loop over observations

            write(unit=lun,fmt='(3(f10.6))')    PolaroMultilist%Polarolist(iset)%Polaro(j)%H
            write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(j)%oPij(:,1), &
            PolariMultisVslist%PolarisVslist(iset)%PolarisVs(j)%Pij(:,1)
            write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(j)%oPij(:,2), &
            PolariMultisVslist%Polarisvslist(iset)%PolarisVs(j)%Pij(:,2)
            write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(j)%oPij(:,3), &
            PolariMultisVslist%PolarisVslist(iset)%PolarisVs(j)%Pij(:,3)

        end do !end loop over observations

      end if

    End Subroutine Calc_Polar_CrSec_Data

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
      type(file_list_type),     intent (in)     :: file_cfl
      !---- Local variables ----!
      character(len=132)   :: line
      integer              :: i,j,n,iset,iobs,i1,i2,i3,i4

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
           i3=index(line,"f2mag_cryopad")
           if(i3 /= 0) then
             do iset=1,MultiData%Nset
                write(lun,'(a,i2,1x,a)') 'iset= ',iset, MultiData%datfile(iset)
                if(Multidata%SNP(iset)) then
                  write(lun,'(a)') '    Pobs                   Pcalc'
                  do iobs=1,MultiData%Nobs(iset) !loop over observations

                   write(unit=lun,fmt='(3(f10.6))')    PolaroMultilist%Polarolist(iset)%Polaro(iobs)%H
                   write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,1), &
                   PolariMultilist%Polarilist(iset)%Polari(iobs)%Pij(:,1)
                   write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,2), &
                   PolariMultilist%Polarilist(iset)%Polari(iobs)%Pij(:,2)
                   write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,3), &
                   PolariMultilist%Polarilist(iset)%Polari(iobs)%Pij(:,3)

                  end do !end loop over observations
                else
                  write(lun,'(a)') '     h      k       l          Imag        sImag      Icalc'
                  Nf2=MhMultilist%MhList(iset)%Nref
                  do n=1,Nf2 ! This is loop over reflections with Int
                    write(unit=lun,fmt="(3f8.4, 3f12.4)") Oblist%Ob(n)%Hd(:,1),Oblist%Ob(n)%Gobs, &
                                     Oblist%Ob(n)%SGobs,Scalef*MhMultilist%MhList(iset)%Mh(n)%sqMiV
                  end do
                end if
             end do
           end if

           i4=index(line,"f2mag_mupad")
           if(i4 /= 0) then
             do iset=1,MultiData%Nset
                write(lun,'(a,i2,1x,a)') 'iset= ',iset, MultiData%datfile(iset)
                if(Multidata%SNP(iset)) then
                  write(lun,'(a)') '    Pobs                   Pcalc'
                  do iobs=1,MultiData%Nobs(iset) !loop over observations

                    write(unit=lun,fmt='(3(f10.6))')    PolaroMultilist%Polarolist(iset)%Polaro(iobs)%H
                    write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,1), &
                    PolariMultisVslist%PolarisVslist(iset)%PolarisVs(iobs)%Pij(:,1)
                    write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,2), &
                    PolariMultisVslist%Polarisvslist(iset)%PolarisVs(iobs)%Pij(:,2)
                    write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,3), &
                    PolariMultisVslist%PolarisVslist(iset)%PolarisVs(iobs)%Pij(:,3)

                  end do !end loop over observations
                else
                 write(lun,'(a)') '     h      k       l          Imag        sImag      Icalc'
                 Nf2=MhMultilist%MhList(iset)%Nref
                 do n=1,Nf2 ! This is loop over reflections with Int
                  write(unit=lun,fmt="(3f8.4, 3f12.4)") Oblist%Ob(n)%Hd(:,1),Oblist%Ob(n)%Gobs, &
                                    Oblist%Ob(n)%SGobs,Scalef*MhMultilist%MhList(iset)%Mh(n)%sqMiV
                 end do
                end if
             end do
           end if

           if(i3 == 0 .and. i4 == 0) then !f2mag pure case
             do iset=1,MultiData%Nset
               write(lun,'(a,i2,1x,a)') 'iset= ',iset, MultiData%datfile(iset)
               write(lun,'(a)') '     h      k       l          Imag        sImag      Icalc'
               Nf2=MhMultilist%MhList(iset)%Nref
               do n=1,Nf2 ! This is loop over reflections with Int
                 write(unit=lun,fmt="(3f8.4, 3f12.4)") Oblist%Ob(n)%Hd(:,1),Oblist%Ob(n)%Gobs, &
                                    Oblist%Ob(n)%SGobs,Scalef*MhMultilist%MhList(iset)%Mh(n)%sqMiV
               end do
             end do
           end if

         else

           i1=index(line,"cryopad")
           if(i1 /= 0) then
             do iset=1,MultiData%Nset
                if(Multidata%SNP(iset)) then
                  write(lun,'(a,i2,1x,a)') 'iset= ',iset, MultiData%datfile(iset)
                  write(lun,'(a)') '    Pobs                   Pcalc'
                  do iobs=1,MultiData%Nobs(iset) !loop over observations

                   write(unit=lun,fmt='(3(f10.6))')    PolaroMultilist%Polarolist(iset)%Polaro(iobs)%H
                   write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,1), &
                   PolariMultilist%Polarilist(iset)%Polari(iobs)%Pij(:,1)
                   write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,2), &
                   PolariMultilist%Polarilist(iset)%Polari(iobs)%Pij(:,2)
                   write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,3), &
                   PolariMultilist%Polarilist(iset)%Polari(iobs)%Pij(:,3)

                  end do !end loop over observations
                end if
             end do
           end if

           i2=index(line,"mupad")
           if(i2 /= 0) then
           do iset=1,MultiData%Nset
              if(Multidata%SNP(iset)) then
                write(lun,'(a,i2,1x,a)') 'iset= ',iset, MultiData%datfile(iset)
                write(lun,'(a)') '    Pobs                   Pcalc'
                do iobs=1,MultiData%Nobs(iset) !loop over observations

                  write(unit=lun,fmt='(3(f10.6))')    PolaroMultilist%Polarolist(iset)%Polaro(iobs)%H
                  write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,1), &
                  PolariMultisVslist%PolarisVslist(iset)%PolarisVs(iobs)%Pij(:,1)
                  write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,2), &
                  PolariMultisVslist%Polarisvslist(iset)%PolarisVs(iobs)%Pij(:,2)
                  write(unit=lun,fmt='(2(1x,3f8.4))') PolaroMultilist%Polarolist(iset)%Polaro(iobs)%oPij(:,3), &
                  PolariMultisVslist%PolarisVslist(iset)%PolarisVs(iobs)%Pij(:,3)

                end do !end loop over observations
              end if
           end do
           end if
         end if

        end select
      end do

   End Subroutine Write_ObsCalc

!******************************************!
   Subroutine MagDom_to_Dataset(Mag_Dom)
!******************************************!
    !!--++ writes Mag_Dom (read from *cfl file or updated by VState_to_AtomsPar after the SA cycle)
    !!--++ into Multidata%MagDom(iset)
    !!--++ Created: February - 2012
    !!--++ Works when both chiral domains given

     !---- Arguments ----!
     Type(Magnetic_Domain_type),   intent(in out)  :: Mag_Dom

     !---- Local variables ----!
     integer              :: iset,i1,j1,i2,j2,idom,dpast,ndomset
     integer              :: place(2)
     real                 :: domsum
     character(len=10)    :: lab

     !----If not all domains are given
      if(Sum(Multidata%MagDom(1:Nset)%nd) /= Mag_dom%nd) then
         write(unit=*,fmt="(a)") " => Mismatch between domains in Mag_Structure and Mag_data! "
         stop
      end if

     do iset=1,Nset

       do i1=1,Mag_Dom%nd
         do j1=1,2 !ich1
           lab=Mag_dom%lab(j1,i1)
           idom = 0

           dom:do i2=1,Multidata%MagDom(iset)%nd
             do j2=1,2 !ich2
               if(Multidata%MagDom(iset)%lab(j2,i2) == lab) then
                 idom=1
                 place=[j2,i2]
                 exit dom
               end if
             end do
           end do dom

           if (idom == 1) then
             j2=place(1); i2=place(2)
             Multidata%MagDom(iset)%DMat(:,:,i2)=Mag_dom%DMat(:,:,i1)
             Multidata%MagDom(iset)%Chir=Mag_dom%Chir
             Multidata%MagDom(iset)%pop(j2,i2) =Mag_dom%pop(j1,i1)
             Multidata%MagDom(iset)%Mpop(j2,i2)=Mag_dom%Mpop(j1,i1)
             Multidata%MagDom(iset)%Lpop(j2,i2)=Mag_dom%Lpop(j1,i1)
           endif

         end do
       end do

     end do !on iset

     !---- Normalization to Sum_of_domains_onedataset=1
     dpast=0
     do iset=1,Nset
       ndomset=Multidata%MagDom(iset)%nd
       domsum=sum(Multidata%MagDom(iset)%pop)

       do i1=1,ndomset
         do j1=1,2 !ich1
           Multidata%MagDom(iset)%pop(j1,i1)=Multidata%MagDom(iset)%pop(j1,i1)/domsum
           Mag_dom%pop(j1,(i1+dpast))=Multidata%MagDom(iset)%pop(j1,i1)
         end do
       end do
       dpast=dpast+ndomset
     end do !on iset

   End Subroutine MagDom_to_Dataset

 End module prep_input