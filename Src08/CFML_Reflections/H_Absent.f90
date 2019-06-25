!!----
!!----
!!----
SubModule (CFML_Reflections) RFL_002
   Contains
   !!----
   !!---- H_ABSENT
   !!----   .True. if the reflection is absent for this Space group
   !!----
   !!---- 20/06/2019
   !!
   Module Function H_Absent(H, SpG) Result(Info)
      !---- Arguments ----!
      integer, dimension(:), intent (in) :: H
      class(SpG_Type),       intent (in) :: SpG
      logical                            :: info

      !---- Local Variables ----!
      integer                            :: i, Dd
      integer, dimension(size(h))        :: k
      integer,dimension(size(h),size(h)) :: Mat
      real(kind=cp)                      :: r1,r2
      real(kind=cp),dimension(size(h))   :: tr

      !> Init
      info=.false.
       
      Dd=size(h)
      do i=1,SpG%Multip
         Mat=SpG%Op(i)%Mat(1:Dd,1:Dd)
         
         k=matmul(h, Mat)
         if (h_equal(h,k)) then
            tr=SpG%Op(i)%Mat(1:Dd, Dd+1)
            r1=dot_product(tr,real(h,kind=cp))
            r2=nint(r1)
            if (abs(r1-r2) > EPS_REF) then
               info=.true.
               exit
            end if
         end if
      end do
   End Function H_Absent
   
   !!----
   !!---- MH_ABSENT
   !!----   .True. if the magnetic reflection is absent for this Space group
   !!----
   !!---- 20/06/2019 
   !!
   Module Function mH_Absent(H,SpG) Result(Info)
      !---- Arguments ----!
      integer, dimension(:), intent (in) :: H
      class(SpG_Type),       intent (in) :: SpG
      logical                            :: info

      !---- Local Variables ----!
      integer,       dimension(size(h))        :: k
      integer,       dimension(size(h),size(h)):: Mat
      real(kind=cp), dimension(size(h))        :: tr
      integer                                  :: i,n_id,Dd
      real(kind=cp)                            :: r1 

      !> Init
      info=.false.
      Dd=size(h)
      
      select case(SpG%mag_type)
         case(1,3)
            n_id=0
            do i=1,SpG%Multip
               Mat=SpG%Op(i)%Mat(1:Dd,1:Dd)
               tr=SpG%Op(i)%Mat(1:Dd,Dd+1)
               k =matmul(h, Mat)
               if (h_equal(h,k)) then
                  r1=cos(TPI*dot_product(Tr,real(h)))
                  n_id=n_id + Trace(Mat*nint(r1))
               end if
            end do
            if (n_id == 0) info=.true.
            
         case(2)
            info=.true.
            
         case(4)
            info=.true.
            do i=1,SpG%Num_aLat
               tr=SpG%aLat_tr(:,i)
               r1=cos(TPI*dot_product(tr,real(h)))
               if (nint(r1) == -1) info=.false.
            end do
      end Select
   End Function mH_Absent
   
   !!----
   !!---- H_LAT_ABSENT
   !!----  .True. if reflection is absent for Lattice conditions
   !!----
   !!---- 20/06/2019
   !!
   Module Function H_Latt_Absent(H, Latt, n) Result(info)
      !---- Arguments ----!
      integer,        dimension(:),  intent (in) :: h
      type(rational), dimension(:,:),intent (in) :: Latt
      integer,                       intent (in) :: n
      logical                                    :: info

      !---- Local Variables ----!
      integer                          :: k,i
      real(kind=cp),dimension(size(h)) :: Lat
      real(kind=cp)                    :: r1,r2

      !> Init 
      info=.false.
      
      do i=1,n
         Lat=Latt(:,i)
         r1=dot_product(Lat,real(h,kind=cp))
         r2=nint(r1)
         k=nint(2.0_cp*r1)
         if (mod(k,2) /= 0) info=.true.
         exit
      end do
   End Function H_Latt_Absent
    
End SubModule RFL_002   