!!----
!!----
!!----
SubModule (CFML_Reflections) RFL_013
   Contains
   
   !!----
   !!---- WRITE_INFO_REFLIST
   !!----    Write information about the Reflection List
   !!----
   !!---- 24/06/2019 
   !!
   Module Subroutine Write_Info_RefList(Reflex, Iunit, Mode)
      !---- Arguments ----!
      type(RefList_Type),         intent(in) :: Reflex
      integer,          optional, intent(in) :: Iunit
      character(len=*), optional, intent(in) :: Mode

      !---- Local variables ----!
      integer, dimension(:,:), allocatable :: hh
      integer               :: i,n,lun
      integer               :: hmax,kmax,lmax,hmin,kmin,lmin
      real(kind=cp)         :: delta

      !> Init
      lun=6
      if (present(iunit)) lun=iunit

      if (present(mode)) then
         select case (l_case(mode(1:3)))
            case("neu")
               write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(NEUTRONS)"
               write(unit=lun,fmt="(a)")     "    ==================================================="
            case default
               write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(X-RAYS)"
               write(unit=lun,fmt="(a)")     "    ================================================="
         end select

      else
         write(unit=lun,fmt="(a)")   "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(X-RAYS)"
         write(unit=lun,fmt="(a)")   "    ================================================="
      end if
      
      if (reflex%nref <=0) then
         write(unit=lun,fmt="(/,a)")   "There aren't reflections on the List!"
         return
      end if   

      n=reflex%nref
      
      !> This part only works with Intel and not
      !> with GFortran 8.1
      
      !hmax=maxval(reflex%ref(1:n)%h(1))
      !kmax=maxval(reflex%ref(1:n)%h(2))
      !lmax=maxval(reflex%ref(1:n)%h(3))
      
      !hmin=minval(reflex%ref(1:n)%h(1))
      !kmin=minval(reflex%ref(1:n)%h(2))
      !lmin=minval(reflex%ref(1:n)%h(3))
      
      if (allocated(hh)) deallocate(hh)
      allocate(hh(n,3))
      do i=1,n
         hh(i,:)=reflex%ref(i)%h(1:3)
      end do
      hmax=maxval(hh(1:n,1))
      kmax=maxval(hh(1:n,2))
      lmax=maxval(hh(1:n,3))
      
      hmin=minval(hh(1:n,1))
      kmin=minval(hh(1:n,2))
      lmin=minval(hh(1:n,3))
            

      write(unit=lun,fmt="(/,a,/)") "   H   K   L   Mult  SinTh/Lda    |Fobs|     sig(Fobs)      |Fc|       Delta"
      select type (r => reflex%ref)
         type is (refl_type)
            do i=1,n
               write(unit=lun,fmt="(3i4,i5,f12.5)") r(i)%h, r(i)%mult, r(i)%S
            end do
            
         type is (srefl_type)
            do i=1,n
               delta=r(i)%Fo-r(i)%Fc
               write(unit=lun,fmt="(3i4,i5,5f12.5)") r(i)%h, r(i)%mult, r(i)%S,   &
                                                     r(i)%Fo,r(i)%SFo, r(i)%Fc, delta
            end do
            
         type is (mrefl_type)
            do i=1,n
            end do
      end select
                  
      write(unit=lun,fmt="(a)") " "
      write(unit=lun,fmt="(a)") " "
      write(unit=lun,fmt="(a,i6)") " => Number of Reflections: ", reflex%nref
      write(unit=lun,fmt="(a,i4,tr3,a,i4,tr3,a,i4)") " => H_max: ",hmax," K_max: ",kmax," L_max: ",lmax
      write(unit=lun,fmt="(a,i4,tr3,a,i4,tr3,a,i4)") " => H_min: ",hmin," K_min: ",kmin," L_min: ",lmin
      write(unit=lun,fmt="(a)") " "

      return
   End Subroutine Write_Info_RefList
    
End SubModule RFL_013   