!!----
!!----
!!----
!!----
SubModule (CFML_SpaceG) SPG_039
   Contains 

    !!----
    !!---- IDENTIFY_CRYSTALLOGRAPHIC_PG
    !!----
    !!----  Determines the crystallographic point group of the group G.
    !!----
    !!---- 24/04/2019 
    !!
    Module Subroutine Identify_Crystallographic_PG(G)
        !---- Arguments ----!
        type(spg_type), intent(in out) :: G

        !---- Local variables ----!
        integer                                            :: i,j,n,d,t
        integer                                            :: numops,nRepSymOp
        logical                                            :: selected
        type(rational)                                     :: det,tr
        integer,               dimension(6)                :: nRot
        type(Symm_Oper_Type),  dimension(:),   allocatable :: repSymOp ! representative operations
        integer,               dimension(:,:), allocatable :: idd

        !> Init
        call Clear_Error()
        
        nRot(:)  = 0  ! number of selected rotations of order 1,2,...,6
        G%pg     = ""
        numops = G%multip / (G%num_lat+G%num_alat+1)
        if (G%centred /= 1 .or. G%anticentred /= 1) numops = numops / 2
        allocate(repSymOp(numops))
        
        !> Get the rotations of the representative matrices
        nRepSymOp = 0
        do i = 1 , G%Multip
           !if (Rational_Equal(G%op(i)%Mat(1:3,1:3),identity)) cycle
           det = rational_determ(G%op(i)%Mat(1:3,1:3))
           if (mod(det%numerator,det%denominator) /= 0) then
              Err_CFML%Ierr = 1
              Err_CFML%Msg = "Identify_Crystallographic_PG@SPACEG: "// &
                             "Determinant is not an integer."
              return
           end if
           d = det%numerator / det%denominator
           selected = .false.
           do j = 1 , nRepSymOp
              if (Rational_Equal(G%op(i)%Mat(1:3,1:3),repSymOp(j)%Mat(1:3,1:3)) .or. &
                 Rational_Equal(G%op(i)%Mat(1:3,1:3),-repSymOp(j)%Mat(1:3,1:3))) then
                 selected = .true.
                 exit
              end if
           end do
           if (.not. selected) then
              nRepSymOp = nRepSymOp + 1
              if (nRepSymOp > numops) then
                 Err_CFML%Ierr = 1
                 Err_CFML%Msg = "Identify_Crystallographic_PG@SPACEG: nRepSymOp > numops"
                 return
              end if
              repSymOp(nRepSymOp) = G%op(i)
              tr  = rational_trace(G%op(i)%Mat(1:3,1:3))
              if (mod(tr%numerator,tr%denominator) /= 0) then
                 Err_CFML%Ierr = 1
                 Err_CFML%Msg = "Identify_Crystallographic_PG@SPACEG: Trace is not an integer."
                 return
              end if
              t = tr%numerator / tr%denominator
              select case (abs(t))
                 case (0)
                    nRot(3) = nRot(3) + 1
                   
                 case(1)
                    if (d * t ==  1) then
                       nRot(4) = nRot(4) + 1
                    else
                       nRot(2) = nRot(2) + 1
                    end if
                   
                 case(2)
                    nRot(6) = nRot(6) + 1
                   
                 case(3)
                    nRot(1) = nRot(1) + 1
              end select
           end if
        end do
        
        !> Get the point group
        allocate(idd(numops,2))
        if (nRot(3) == 8) then ! Cubic
           if (nRepSymOp == 12) then
              if (G%Centred == 1 .and. G%anticentred == 1) then
                 G%PG = "23"
              else
                 G%PG = "m-3"
              end if
           
           else if (nRepSymOp == 24) then
              if (G%Centred == 1 .and. G%anticentred == 1) then
                 call Get_Rotations(repSymOp,nRepSymOp,4,n,idd)
                 if (n == 6 .and. idd(1,2) == 1) then
                    G%PG = "432"
                 else if (n == 6 .and. idd(1,2) == -1) then
                    G%PG = "-43m"
                 end if
              else
                 G%PG = "m-3m"
              end if
           end if
        
        else if (nRot(6) == 2) then ! Hexagonal
           if (nRepSymOp == 6) then
              if (G%Centred == 1 .and. G%anticentred == 1) then
                 call Get_Rotations(repSymOp,nRepSymOp,6,n,idd)
                 if (idd(1,2) == 1) then
                    G%PG = "6"
                 else
                    G%PG = "-6"
                 end if
              else
                 G%PG = "6/m"
              end if
           
           else if (nRepSymOp == 12) then
              if (G%Centred == 1 .and. G%anticentred == 1) then
                 call Get_Rotations(repSymOp,nRepSymOp,6,n,idd)
                 if (idd(1,2) == 1) then
                    call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
                    if (n == 7 .and. (idd(1,2) == -1 .or. idd(2,2) == -1)) then
                       G%PG = "6mm"
                    else
                       G%PG = "622"
                    end if
                 
                 else if (idd(1,2) == -1) then
                    G%PG = "-6m2"
                 end if
                
              else
                 G%PG = "6/mmm"
              end if
           end if
        
        else if (nRot(3) == 2) then ! Trigonal
           if (nRepSymOp == 3) then
              if (G%Centred == 1 .and. G%anticentred == 1) then
                 G%PG = "3"
              else
                 G%PG = "-3"
              end if
            
           else if (nRepSymOp == 6) then
              if (G%Centred == 1 .and. G%anticentred == 1) then
                 call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
                 if (n == 3 .and. idd(1,2) == 1) then
                    G%PG = "32"
                 else if (n == 3 .and. idd(1,2) == -1) then
                    G%PG = "3m"
                 end if
              else
                 G%PG = "-3m"
              end if
           end if
        
        else if (nRot(4) == 2) then ! Tetragonal
           if (nRepSymOp == 4) then
              if (G%Centred == 1 .and. G%anticentred == 1) then
                 call Get_Rotations(repSymOp,nRepSymOp,4,n,idd)
                 if (n == 2 .and. idd(1,2) == 1) then
                    G%PG = "4"
                 else if (n == 2 .and. idd(1,2) == -1) then
                    G%PG = "-4"
                 end if
              else
                 G%PG = "4/m"
              end if
           
           else if (nRepSymOp == 8) then
              if (G%Centred == 1 .and. G%anticentred == 1) then
                 call Get_Rotations(repSymOp,nRepSymOp,4,n,idd)
                 if (n == 2 .and. idd(1,2) == 1) then
                    call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
                    if (n == 5 .and. (idd(1,2) == -1 .or. idd(2,2) == -1)) then
                       G%PG = "4mm"
                    else
                       G%PG = "422"
                    end if
                 
                 else if (n == 2 .and. idd(1,2) == -1) then
                    G%PG = "-4m2"
                 end if
              
              else
                 G%PG = "4/mmm"
              end if
           end if
        
        else if (nRot(2) == 3) then ! Orthorhombic
           if (G%Centred == 1 .and. G%anticentred == 1) then
              call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
              if (n == 3 .and. idd(1,2) == -1 .or. idd(2,2) == -1) then
                 G%PG = "mm2"
              else
                 G%PG = "222"
              end if
           
           else
               G%PG = "mmm"
           end if
        
        else if (nRot(2) == 1) then ! Monoclinic
           if (G%Centred == 1 .and. G%anticentred == 1) then
              call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
              if (idd(1,2) == 1) then
                 G%PG = "2"
              else
                 G%PG = "m"
              end if
            
            else
               G%PG = "2/m"
            end if
        
        else
           if (G%Centred == 1 .and. G%anticentred == 1) then ! Triclinic
              G%PG = "1"
           else
              G%PG = "-1"
           end if
        end if

        if (G%PG == "") then
           Err_CFML%Ierr = 1
           Err_CFML%Msg = "Identify_Crystallographic_PG@SPACEG: Unable to identify point group"
        end if

    End Subroutine Identify_Crystallographic_PG
    
End Submodule SPG_039    