!!----
!!----
!!----
SubModule (CFML_Reflections) RFL_Init_Reflist
   implicit none
   Contains

   !!----
   !!---- INITIALIZE_REFLIST()
   !!----    Initialize the Reflection List Variable
   !!----
   !!---- 24/06/2019
   !!
   Module Subroutine Initialize_RefList(N, Reflex)
      !---- Arguments ----!
      integer,             intent(in)    :: N
      type(RefList_Type),  intent(in out) :: Reflex

      !---- Local Variables ----!
      integer :: i

      if (allocated(reflex%ref)) deallocate(reflex%ref)
      select case (n)
         case (0)
             reflex%Nref=0

         case (1:)
            reflex%Nref=n
            allocate(reflex%ref(n))

            associate (r => reflex%ref)
               select type (r)
                 class is (Refl_Type)
                    do i=1,n
                       r(i)%h     =0
                       r(i)%mult  =0
                       r(i)%s     =0.0_cp
                       r(i)%Imag  =0
                       r(i)%Pcoeff=0
                    end do
               end select

               select type (r)
                 class is (SRefl_Type)
                    do i=1,n
                       r(i)%Fo   =0.0_cp
                       r(i)%Fc   =0.0_cp
                       r(i)%sFo  =0.0_cp
                       r(i)%Phase=0.0_cp
                       r(i)%A    =0.0_cp
                       r(i)%B    =0.0_cp
                       r(i)%W    =1.0_cp
                    end do
               end select

               select type (r)
                 type is (MRefl_Type)
                    do i=1,n
                       r(i)%smIvo=0.0_cp
                       r(i)%msF=cmplx(0.0_cp,0.0_cp)
                       r(i)%mIv=cmplx(0.0_cp,0.0_cp)
                    end do
               end select
            end associate

      end select

   End Subroutine Initialize_RefList

End SubModule RFL_Init_Reflist