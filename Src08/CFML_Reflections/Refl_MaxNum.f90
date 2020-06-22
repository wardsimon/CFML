!!----
!!----
!!----
SubModule (CFML_Reflections) RFL_Maximum_Number_of_Reflections
   Contains

   !!----
   !!---- GET_MAXNUMREF
   !!----    Provides un upper limit of the expected maximum number of
   !!----    reflections up to SinTLMax for a volume VolCell of the
   !!----    primitive cell.
   !!----    If SinTLMin is given, the result is the number of reflections
   !!----    in the interval (SinTLMin,SinTLMax).
   !!----    If Mult is provided the result is divided by half this multiplicity
   !!----    so we obtain an estimation of the expected mumber of unique reflections.
   !!----
   !!---- 21/06/2019
   !!
   Module Function Get_MaxNumRef(SinTLMax, VolCell, SinTLMin, Mult) Result(numref)
      !---- Arguments ----!
      real(kind=cp),           intent(in) :: SinTLMax    ! Maximum sinTheta/Lambda
      real(kind=cp),           intent(in) :: VolCell     ! Direct Cell Volume
      real(kind=cp), optional, intent(in) :: SinTLMin    ! Minimum sinTheta/Lambda
      integer,       optional, intent(in) :: Mult        ! General Multiplicity
      integer                             :: numref

      !---- Local Variables ----!
      real(kind=cp) :: r3

      r3=8.0_cp * SinTLMax * SinTLMax * SinTLMax * 1.05_cp

      if (present(SinTLMin)) r3= r3 - 8.0_cp*SinTLMin * SinTLMin * SinTLMin

      numref=nint(4.0_cp*PI*r3*VolCell / 3.0_cp)

      !> The factor 2 is given because, for high symmetry, sometimes the obtained number
      !> is not enough for allocating the real number of reflections
      if (present(Mult)) numref=2 * numref/max(1,Mult)
   End Function Get_MaxNumRef
End SubModule RFL_Maximum_Number_of_Reflections