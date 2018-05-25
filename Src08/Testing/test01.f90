Module Module_Test
   !---- Use Modules ----!
   Use CFML_GlobalDeps
   Use CFML_Crystal_Metrics
   
   implicit none
   
   !---- Variables /Definitions /Parameters
   
 Contains
 
    
 
End Module Module_Test

Program Testing
    !---- Use modules ----!
    Use Module_Test
    
    implicit none
    
    real, dimension (3) :: abc, ang,sabc,sang
    class(CrysCell_Type) , allocatable :: celda
    
    !> Init values
    celda=crysCell_type([0.0,0.0,0.0],[0.0,0.0,0.0], &
                        [0.0,0.0,0.0],[0.0,0.0,0.0],0.0,0.0)
    
    abc =[12.2134, 17.2091, 8.0012]
    ang =[84.23, 107.12, 95.88]
    sabc=[0.0003,0.0004,0.0002]
    sang=[0.02,0.03,0.01]
    
    
    print*,abc
    print*,ang
    print*,sabc
    print*,sang
    print*,'Llamando a Set_Crystal_Cell'
    
    call Set_crystal_Cell(abc,ang,celda," ",sabc, sang)
    print*,'Llamando a Write_Crystal_Cell'
    call write_crystal_cell(celda)

End Program Testing
