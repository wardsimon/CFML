Module Module_Test
   !---- Use Modules ----!
   
   implicit none
   
   !---- Variables /Definitions /Parameters
   
 Contains
 
    !!--++
    !!--++ Subroutine Recip(A,Ang,Ar,Angr,Vol,Volr)
    !!--++    real(kind=cp), dimension(3), intent(in ) :: a        !  In -> a,b,c
    !!--++    real(kind=cp), dimension(3), intent(in ) :: ang      !  In -> alpha,beta,gamma
    !!--++    real(kind=cp), dimension(3), intent(out) :: ar       !  In -> a*,b*,c*
    !!--++    real(kind=cp), dimension(3), intent(out) :: angr     !  In -> alpha*,beta*,gamma*
    !!--++    real(kind=cp),               intent(out) :: vol      ! Out -> Vol
    !!--++    real(kind=cp),               intent(out) :: volr     ! Out -> Vol*
    !!--++
    !!--++    (PRIVATE)
    !!--++    Calculates the reciprocal lattice vectors and cell volume
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Recip(A,Ang,Ar,Angr,Vol,Volr)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in ) :: a,ang
       real(kind=cp), dimension(3), intent(out) :: ar,angr
       real(kind=cp),               intent(out) :: vol,volr

       !---- Local Variables ----!
       integer        :: i
       real(kind=cp)  :: s,p,cose

       p=1.0
       s=1.0
       do i=1,3
          cose=cosd(ang(i))
          p=p*cose
          s=s-cose*cose
       end do
       vol=sqrt(abs(s+2.0*p))

       do i=1,3
          vol=vol*a(i)
       end do
       volr=1.0/vol

       ar(1)=a(2)*a(3)*sind(ang(1))/vol
       ar(2)=a(3)*a(1)*sind(ang(2))/vol
       ar(3)=a(1)*a(2)*sind(ang(3))/vol
       angr(1)=(cosd(ang(2))*cosd(ang(3))-cosd(ang(1)))/(sind(ang(2))*sind(ang(3)))
       angr(2)=(cosd(ang(1))*cosd(ang(3))-cosd(ang(2)))/(sind(ang(1))*sind(ang(3)))
       angr(3)=(cosd(ang(2))*cosd(ang(1))-cosd(ang(3)))/(sind(ang(2))*sind(ang(1)))
       do i=1,3
          angr(i)=acosd(angr(i))
       end do

       return
    End Subroutine Recip
    
End Module_Test

Program Testing
    !---- Use modules ----!
    Use CFML_Crystal_Metrics
    Use Module_Test
    
    implicit none
    
    class(CrysCell_Type) , allocatable :: celda
    
    
    celda=cryscell_type([3.0, 2.0, 1.0], [90.0, 107.12, 90.0],[0.1,0.1,0.2],[0.1,0.5,0.1] ,0.0,0.0)
    
    celda%vol=Volume_Cell(celda%cell, celda%ang)
    celda%svol=SigmaVolume(celda)
    print*," Volumen de la Celda: ",celda%vol
    print*," Sigma: ",celda%svol

    

End Program Testing
