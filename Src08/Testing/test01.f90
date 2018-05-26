Module Module_Test
   !---- Use Modules ----!
   Use CFML_GlobalDeps
   Use CFML_Crystal_Metrics
   Use CFML_String_Utilities, only: U_Case
   
   implicit none
   
   !---- Variables /Definitions /Parameters
   
 Contains
    !!----
    !!---- SUBROUTINE GET_CRYST_FAMILY
    !!----
    !!----    Obtain the Crystal Family, Symbol and System from cell parameters
    !!----
    !!---- Update: May - 2005
    !!----
    Subroutine Get_Cryst_Family(Cell, Str_Family, Str_Symbol, Str_System)
       !---- Arguments ----!
       class(CrysCell_Type),   intent(in ) :: Cell
       character(len=*),       intent(out) :: Str_Family
       character(len=*),       intent(out) :: Str_Symbol
       character(len=*),       intent(out) :: Str_System

       !---- Local variables ----!
       integer, dimension(3) :: icodp, icoda
       integer               :: n1,n2

       !> Init
       Str_Family=" "
       Str_Symbol=" "
       Str_System=" "

       icodp=0; icoda=0

       !> Codification 

       !> a
       icodp(1)=1

       !> b 
       if (abs(cell%cell(2)-cell%cell(1)) <= eps) then
          icodp(2)=icodp(1)
       else
          icodp(2)=2
       end if

       !> c
       if (abs(cell%cell(3)-cell%cell(1)) <= eps) then
          icodp(3)=icodp(1)
       else
          icodp(3)=3
       end if

       !> alpha
       icoda(1)=1

       !> beta
       if (abs(cell%ang(2)-cell%ang(1)) <= eps) then
          icoda(2)=icoda(1)
       else
          icoda(2)=2
       end if

       !>gamma 
       if (abs(cell%ang(3)-cell%ang(1)) <= eps) then
          icoda(3)=icoda(1)
       else
          icoda(3)=3
       end if

       n1=count(icoda==icoda(1))  ! angles
       n2=count(icodp==icodp(1))  ! parameters
       
       select case (n1)
          case (1) ! All angles are differents
             Str_Family="Triclinic"
             Str_Symbol ="a"
             Str_System ="Triclinic"

          case (2) ! two angles are equal
             if (icoda(1) == icoda(2)) then
                if (abs(cell%ang(3)-120.0) <= eps) then
                   if (icodp(1)==icodp(2)) then
                      !> Hexagonal 
                      Str_Family="Hexagonal"
                      Str_Symbol ="h"
                      Str_System ="Hexagonal"
                   end if
                else
                   !> Monoclinic 
                   Str_Family="Monoclinic"
                   Str_Symbol ="m"
                   Str_System ="Monoclinic"
                end if

             else
                !> Monoclic b-unique setting
                if (abs(cell%ang(1)-90.0) <= eps) then
                   Str_Family="Monoclinic"
                   Str_Symbol ="m"
                   Str_System ="Monoclinic"
                end if
             end if

          case (3) ! all angles are equals
             if (abs(cell%ang(1) - 90.000) <= eps) then
                select case (n2)
                   case (1) 
                      !> Orthorhombic 
                      Str_Family="Orthorhombic"
                      Str_Symbol ="o"
                      Str_System ="Orthorhombic"

                   case (2)
                      !> Tetragonal 
                      if (icodp(1)==icodp(2)) then
                         Str_Family="Tetragonal"
                         Str_Symbol ="t"
                         Str_System ="Tetragonal"
                      end if   

                   case (3)
                      !> Cubic
                      Str_Family="Cubic"
                      Str_Symbol ="c"
                      Str_System ="Cubic"
                end select

             else
                if (n2 == 3) then
                   !> Hexagonal with rhombohedral axes 
                   Str_Family="Hexagonal"
                   Str_Symbol ="h"
                   Str_System ="Trigonal"
                end if
             end if

       end select ! n
       
       !> Error?
       if (len_trim(Str_Family) <= 0) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          err_CFML%Msg=" Error obtaining the Crystal Family. Please, check the cell parameters!"
       end if 

       return
    End Subroutine Get_Cryst_Family
    
    
    
 
End Module Module_Test

Program Testing
    !---- Use modules ----!
    Use Module_Test
    
    implicit none
    
    real, dimension (3) :: abc, ang,sabc,sang
    character(len=40)   :: car1, car2, car3
    class(CrysCell_Type) , allocatable :: celda, celda2
    
    !> Init values (We must initialize the class variable!!!!)
    if (.not. allocated(celda)) then
       print*,' ALLOCATING VARIABLE !!!!'
       celda=crysCell_type([0.0,0.0,0.0],[0.0,0.0,0.0], &
                           [0.0,0.0,0.0],[0.0,0.0,0.0],0.0,0.0)
    end if                       
                        
    celda2=cryscell_m_type([0.0,0.0,0.0],[0.0,0.0,0.0], &
                           [0.0,0.0,0.0],[0.0,0.0,0.0], &
                            0.0,0.0,                    & 
                           [0.0,0.0,0.0],[0.0,0.0,0.0], 0.0, &                   
                           reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,3]), &
                           reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,3]), &
                           reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,3]), &
                           reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,3]), &
                           reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,3]), &
                           reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,3]), " ")
                           
    abc =[12.2134, 17.2091, 8.0012]
    ang =[84.23, 107.12, 95.88]
    sabc=[0.0003,0.0004,0.0002]
    sang=[0.02,0.03,0.01]
    
    !> celda 1
    call Set_crystal_Cell(abc,ang,celda," ",sabc, sang)
    call write_crystal_cell(celda)
    call Get_Cryst_Family(celda,car1, car2, car3)
    write(unit=*,fmt='(a)') "Family: "//trim(car1)//"   Symbol: "//trim(car2)//"     System: "//trim(car3)
    
    !> celda 2
    ang =[90.0, 90.0, 90.0]
    sang=[0.0,0.0,0.0]
    call Set_crystal_Cell(abc,ang,celda2," ",sabc, sang)
    call write_crystal_cell(celda2)
    call Get_Cryst_Family(celda2,car1, car2, car3)
    write(unit=*,fmt='(a)') "Family: "//trim(car1)//"   Symbol: "//trim(car2)//"     System: "//trim(car3)
    
    

End Program Testing
