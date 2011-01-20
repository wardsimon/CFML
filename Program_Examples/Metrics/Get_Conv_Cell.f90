  !!---- Program Get_Conventional_Unit_Cells
  !!----
  !!----  Author: Juan Rodriguez-Carvajal, ILL (2008)
  !!----
  !!---- This program needs as input unit cell parameters and provides the conventional
  !!---- unit cell parameters and the trasformation matrix between the input cell and
  !!---- the conventional cell(s). The program is similar to REDUC from Yvon Le Page but
  !!---- it has been developed from the scratch using new procedures that are in CrysFML.
  !!---- These procedures are based on well known articles of the literature about
  !!---- reduced and conventional cells. In particular the search of two-fold axes is
  !!---- based on the paper by Yvon Le Page in J.Appl.Cryst. 15, 255 (1982) and the search
  !!---- of the Niggli cell is based on the algorithm developed by
  !!---- I. Krivy and B. Gruber, Acta Cryst A32, 297 (1976).
  !!----
  !!----           (Api) = transfm (Aic)      (AN)  = trans (Api)
  !!----           (AN)  = trans . transfm   (Aic)  = Prod (Aic)
  !!----                    prod=matmul(trans,transfm)
  !!----           (Acc) = Tr (AN) = Tr Prod (Aic) = Tr . trans . transfm (Aic)

  Program Get_Conventional_Unit_Cells
    use CFML_String_Utilities, only: Ucase
    use CFML_Crystal_Metrics,  only: Crystal_Cell_Type, Err_Crys, ERR_Crys_Mess, &
                                     Set_Crystal_Cell,Get_Cryst_Family, Niggli_Cell,     &
                                     Write_Crystal_Cell, Get_Primitive_Cell,Get_two_fold_axes,  &
                                     Twofold_Axes_Type, Get_Conventional_Cell, Change_Setting_Cell
    use CFML_Math_3D,          only: determ_A, determ_V, invert_A
    use CFML_Math_General,     only: sort

    implicit none

    !
    !  Examples of cells for testing the program:                       Input           Output
    !  4.000 4.472  4.583   79.03  64.13  64.15                            P       !Pseudo-tetragonal I
    !  4.583 4.472  4.000   64.15  64.13  79.03                            P       !Pseudo-tetragonal I
    !  5.211 5.222  5.209   59.80  60.32  59.85                            P       !Pseudo-cubic F
    !  2.04230261 2.51256847 2.66664577 114.966057 108.821213 97.3560257   P       !Monoclinic C
    !  1.86333036 2.08734274 2.78944445 81.4266052 70.4884567 63.4908752   P       !Orthorhombic F
    !  1.863 2.087  2.789   81.2   70.30  63.29                            P       !Pseudo-Orthorhombic F
    !  5.211 5.212  5.209   59.90  60.12  59.85                            F       !Pseudo-Rhombohedral
    !  5.234 9.065  11.91   89.91  90.01  90.03                            C       !Pseudo-Hexagonal

    real, dimension(3,3)    :: transfm,trans,prod,finalm,mat
    integer, dimension(3,3) :: tr
    real, dimension(6)      :: ad
    real, dimension(12)     :: del
    integer, dimension(12)  :: ind
    character(len=1)        :: ls,ans,mon
    character(len=80)       :: message
    character(len=11)       :: metr
    logical                 :: ok,cell_trans
    integer                 :: i,j,p,n,ier, det, ntwot, nold , io=6 !Standard output
    real                    :: rmi,rma, tol, told
    Type(Twofold_Axes_Type) :: twofold,twf
    Type(Twofold_Axes_Type), dimension(12) :: otwf   !Individual two-fold axes for searching
                                                     !monoclinic cells
    type(Crystal_Cell_Type) :: cellc, cellp, celln, cell,cellt

      del=0.0
      write(unit=*,fmt="(/,a)") " --------------------------------------------"
      write(unit=*,fmt="(a)")   "     Program: Get_Conventional_Unit_Cells    "
      write(unit=*,fmt="(a)")   "  Author: J. Rodriguez-Carvajal, ILL(2008)   "
      write(unit=*,fmt="(a)")   "            Updated January 2011             "
      write(unit=*,fmt="(a,/)") " --------------------------------------------"
      write(unit=*,fmt="(a)",advance="no") " => Output on the screen (<cr>) or in file (f) : "
      read(unit=*,fmt="(a)") ans
      if(ans == "f") then
        io=1
        open(unit=io,file="conventional_cells.out",status="replace",action="write")
        write(unit=io,fmt="(/,a)") " --------------------------------------------"
        write(unit=io,fmt="(a)")   "     Program: Get_Conventional_Unit_Cells    "
        write(unit=io,fmt="(a)")   "  Author: J. Rodriguez-Carvajal, ILL(2008)   "
        write(unit=io,fmt="(a)")   "            Updated January 2011             "
        write(unit=io,fmt="(a,/)") " --------------------------------------------"
      end if
      do
        write(unit=*,fmt="(a)",advance="no") " => Enter the cell parameters: "
        read(unit=*,fmt=*,iostat=ier) ad
        if( ier == 0) exit
      end do
      Call Set_Crystal_Cell(ad(1:3),ad(4:6), Cellc)
      if(Err_Crys) then
        write(unit=*,fmt="(a)") " => "//trim(ERR_Crys_Mess)
        write(unit=*,fmt="(a)") " => The provided unit cell parameters contain some error! Abnormal end!"
        stop
      end if
      do
        write(unit=*,fmt="(a)",advance="no") " => Enter the centring type (P,A,B,C,I,R,F): "
        read(unit=*,fmt=*,iostat=ier) ls
        if( ier == 0) exit
      end do
      call ucase(ls)
      if(index("PABCIRF",ls) == 0) ls="P"

      call Get_Primitive_Cell(ls,Cellc,cellp,transfm)
      call Niggli_Cell(Cellp,celln=celln,trans=trans)
      if(Err_Crys) then
       write(unit=*,fmt="(a)") " => "//trim(ERR_Crys_Mess)
       stop
      end if

      tol=3.0
      write(unit=*,fmt="(a)",advance="no") " => Enter angular tolerance in degrees (<cr> = 3 deg.): "
      read(unit=*,fmt="(a)") message
      if(len_trim(message) /= 0) then
        read(unit=message,fmt=*,iostat=ier) tol
        if(ier /= 0) tol=3.0
        message=" "
      end if
      told=0.2
      write(unit=*,fmt="(a)",advance="no") " => Enter distance tolerance in angstroms (<cr> = 0.2 A.): "
      read(unit=*,fmt="(a)") message
      if(len_trim(message) /= 0) then
        read(unit=message,fmt=*,iostat=ier) told
        if(ier /= 0) told=0.2
        message=" "
      end if
      call Get_two_fold_axes(celln,tol,twofold)
      write(unit=io,fmt="(/,a,6f10.4,a)")"             Input Cell (Aic): ",ad,"  Centring: "//ls
      write(unit=io,fmt="(a,6f10.4)")    "   Primitive input Cell (Api): ",Cellp%cell,Cellp%ang
      write(unit=io,fmt="(a,6f10.4)")    "            Niggli Cell  (AN): ",celln%cell,celln%ang
      write(unit=io,fmt="(/,a)")  "             (Aic) formal column matrix containing the input cell vectors: "
      write(unit=io,fmt="(a)")    "             (Api) formal column matrix containing the primitive cell vectors: "
      write(unit=io,fmt="(a)")    "              (AN) formal column matrix containing the Niggli cell vectors: "
      write(unit=io,fmt="(a,/)")  "             (Acc) formal column matrix containing the conventional cell vectors: "
      write(unit=io,fmt="(a)")   "                            (Api) = M (Aic)                          (AN) = N (Api)"
      do i=1,3
        write(unit=io,fmt="(a,3f12.6,tr5,3f12.6)") "     TransF:    ",transfm(i,:),trans(i,:)
      end do
      prod=matmul(trans,transfm)
      metr= "    Pseudo-"
      if( twofold%ntwo > 0) then
        write(unit=io,fmt="(a)")        " => Two-fold axes (indices in the Niggli cell)"
        write(unit=io,fmt="(a)")        " =>       Direct       Reciprocal    Dot    Cross      Length "
        rma=-100.0
        do i=1,twofold%ntwo
         write(unit=io,fmt="(a,3i4,tr4,3i4,i6,f10.3,f12.5)")  "     ", &
            twofold%dtwofold(:,i),twofold%rtwofold(:,i),twofold%dot(i),twofold%cross(i),twofold%maxes(i)
         if( twofold%cross(i) > rma) rma= twofold%cross(i)
         del(i)=twofold%cross(i)
        end do
        call sort(del,twofold%ntwo,ind)

        ntwot = twofold%ntwo
        do i=1,twofold%ntwo
          j=ind(i)
          del(twofold%ntwo-i+1)=twofold%cross(j)  !reorder by decreasing errors
          otwf(twofold%ntwo-i+1)%ntwo=1
          otwf(twofold%ntwo-i+1)%tol=tol
          otwf(twofold%ntwo-i+1)%caxes(:,1)=twofold%caxes(:,j)
          otwf(twofold%ntwo-i+1)%dtwofold(:,1)=twofold%dtwofold(:,j)
          otwf(twofold%ntwo-i+1)%rtwofold(:,1)=twofold%rtwofold(:,j)
          otwf(twofold%ntwo-i+1)%dot(1)=twofold%dot(j)
          otwf(twofold%ntwo-i+1)%cross(1)=twofold%cross(j)
          otwf(twofold%ntwo-i+1)%maxes(1)=twofold%maxes(j)
          otwf(twofold%ntwo-i+1)%a(:)=twofold%a(:)
          otwf(twofold%ntwo-i+1)%b(:)=twofold%b(:)
          otwf(twofold%ntwo-i+1)%c(:)=twofold%c(:)
        end do

        !!!------ FIRST: test with all twofold axes found
        call Get_Conventional_Cell(twofold,Cell,tr,message,ok,told)
        det=determ_A(tr)
        !Here Tr is the matrix transforming the Niggli cell to the conventional cell
        if(ok) then
          if(rma < 0.1)  metr="Metrically "
          if(abs(det) > 0) then
            call Write_New_Cell()
            call Write_New_Monoc_Cell()
          end if
          !call Write_Crystal_Cell(Cell)

         !!!------ SECOND: select smaller blocks of twofold axes
          if( rma >= 0.1) then
            p=twofold%ntwo+1
            del(1)=del(1)-0.06
            nold=0
            do j=1,p-1       !Loop decreasing the tolerance for selecting better 2-fold axes
              rmi=del(j)+0.05
              if(rmi < 0.1) then
                rmi=0.1
                metr="Metrically "
              else
                metr="    Pseudo-"
              end if
              n=0
              do i=1,twofold%ntwo
                if(twofold%cross(i) > rmi) cycle
                n=n+1
                twf%caxes(:,n)    = twofold%caxes(:,i)
                twf%dtwofold(:,n) = twofold%dtwofold(:,i)
                twf%rtwofold(:,n) = twofold%rtwofold(:,i)
                twf%dot(n)        = twofold%dot(i)
                twf%cross(n)      = twofold%cross(i)
                twf%maxes(n)      = twofold%maxes(i)
              end do

              if( n == nold ) then
                 cycle
              else
                 nold=n
              end if
              if( n >= 1 ) then
                twf%ntwo=n
                twf%tol=twofold%tol
                twf%a=twofold%a
                twf%b=twofold%b
                twf%c=twofold%c
                call Get_Conventional_Cell(twf,Cell,tr,message,ok)
                det=determ_A(tr)
                if(ok) then
                  call Write_New_Cell()
                  call Write_New_Monoc_Cell()
                  write(unit=io, fmt="(/,a,i3,a,f10.4)")  " => Two-fold axes:  ",n, "  Angular Discrepancy (deg)", rmi
                  write(unit=io,fmt="(a)")           " =>       Direct       Reciprocal    Dot    Cross      Length "
                  do i=1,n
                    write(unit=io,fmt="(a,3i4,tr4,3i4,i6,f10.3,f12.5)")  "     ",&
                       twf%dtwofold(:,i),twf%rtwofold(:,i), twf%dot(i),twf%cross(i),twf%maxes(i)
                  end do
                !else
                !  write(unit=io, fmt="(  a)")   " => No lattice results: "//trim(message)
                end if
              end if
              if(n <= 1) exit
              twofold=twf   !copy back on Twofold
            end do
          end if !rma >0.1
        else
          write(unit=io,fmt="(a)") " => An unexpected error has occurred: "//trim(message)
          write(unit=io,fmt="(a)") " => Change the angular/distance tolerance to obtain proper two-fold axes."
        end if !ok
      else
        !The Niggli cell is accepted as triclinic cell
        write(unit=io,fmt="(a,3f10.5,3f9.3)") "  Cell (Triclinic, Niggli Cell): ",Celln%cell,Celln%ang
        write(unit=io,fmt="(a,3f12.6,a)")     "                                   /",prod(1,:),"\"
        write(unit=io,fmt="(a,3f12.6,a)")     "     Final Tranformation Matrix:  | ",prod(2,:)," |"
        write(unit=io,fmt="(a,3f12.6,a)")     "                                   \",prod(3,:),"/"
      end if

     !!!------ THIRD: Determine all monoclinic cells considering one-by-one all the two-fold axes
     write(unit=*,fmt="(/,a,i3,a)")       " => There are ",ntwot," two-fold axes"
     write(unit=*,fmt="(a)",advance="no") " => Do you want to display the monoclinic cells? (<cr>=no): "
     read(unit=*,fmt="(a)") mon
     if(mon == "y" .or. mon == "Y") then
       write(unit=io, fmt="(/,/,/,130a)") " ",("*",i=1,120)
       write(unit=io,fmt="(  a,/,a)")   &
       " => MONOCLINIC CELLS: In this section each two-fold axis is considered as a unique b-axis of M-cells",&
       "    The metrics of the monoclinic cells may be of higher symmetry, e.g. beta=90.0"
       write(unit=io, fmt="(130a,/)") " ",("*",i=1,120)
       do j=1,ntwot
         call Get_Conventional_Cell(otwf(j),Cell,tr,message,ok)
         det=determ_A(tr)
         rma=otwf(j)%cross(1)
         if(rma < 0.1) then
           metr="Metrically "
         else
           metr="    Pseudo-"
         end if
         if(abs(det) > 0) then
          call Write_New_Cell()
          call Write_New_Monoc_Cell()
         end if
         write(unit=io, fmt="(/,a,f10.4)")  " => Two-fold axis along monoclinic b-axis,  Angular Discrepancy (deg)", rma
         write(unit=io,fmt="(a)")           " =>       Direct       Reciprocal    Dot    Cross      Length "
         write(unit=io,fmt="(a,3i4,tr4,3i4,i6,f10.3,f12.5)")  "     ",&
         otwf(j)%dtwofold(:,1),otwf(j)%rtwofold(:,1), otwf(j)%dot(1),otwf(j)%cross(1),otwf(j)%maxes(1)
       end do
     end if
     if(ans == "f") write(unit=*,fmt="(a)") " => Normal End of the program, results in file: conventional_cells.out"
     stop

  Contains

    Subroutine Write_New_Cell()
      write(unit=io, fmt="(/,130a)") " ",("-",i=1,120)
      write(unit=io,fmt="(  a)")   &
      " => The new Cell is "//metr//trim(message)//"  and the transformation matrix from then Niggli cell is:"
      finalm=matmul(real(tr),prod)
      write(unit=io, fmt="(130a,/)") " ",("-",i=1,120)
      write(unit=io,fmt="(a,i3,2i4,a)")     "                         /",tr(1,:), " \"
      write(unit=io,fmt="(a,i3,2i4,a,i4)") "  (Acc) = Tr (AN);  Tr: | ",tr(2,:), "  |   Determinant: ",det
      write(unit=io,fmt="(a,i3,2i4,a,/)")   "                         \",tr(3,:), " /"
      write(unit=io,fmt="(a,3f10.5,3f9.3)") "  Conventional Cell: ",Cell%cell,Cell%ang
      write(unit=io,fmt="(a,3f12.6,a)")     "                                   /",finalm(1,:),"\"
      write(unit=io,fmt="(a,3f12.6,a)")     "     Final Tranformation Matrix:  | ",finalm(2,:)," |       (Acc) = Ftr (Aic)"
      write(unit=io,fmt="(a,3f12.6,a)")     "                                   \",finalm(3,:),"/"
      return
    End Subroutine Write_New_Cell

    Subroutine Write_New_Monoc_Cell()
      cell_trans=.false.
      if(trim(message) == "Monoclinic, A-centred cell") then
         mat=reshape( (/0.0,0.0,1.0, 0.0,-1.0,0.0, 1.0,0.0,0.0/),(/3,3/))
         call Change_Setting_Cell(Cell,Mat,Cellt)
         cell_trans=.true.
      else if(trim(message) == "Monoclinic, I-centred cell") then
         mat=reshape( (/1.0,0.0,1.0, 0.0,1.0,0.0, 0.0,0.0,1.0/),(/3,3/))
         call Change_Setting_Cell(Cell,Mat,Cellt)
         if(Cellt%ang(2) < 80.0) then
           mat=reshape( (/1.0,0.0,-1.0, 0.0, 1.0,0.0, 0.0,0.0, 1.0/),(/3,3/))
           call Change_Setting_Cell(Cell,Mat,Cellt)
         end if
         cell_trans=.true.
      end if
      if(cell_trans) then
        finalm=matmul(Mat,finalm)
        write(unit=io,fmt="(/,a,3f10.5,3f9.3)") "  Equivalent C-centred Cell: ",Cellt%cell,Cellt%ang
        write(unit=io,fmt="(a,3f12.6,a)")     "                                   /",finalm(1,:),"\"
        write(unit=io,fmt="(a,3f12.6,a)")     "     Final Tranformation Matrix:  | ",finalm(2,:)," |       (Acc) = Ftr (Aic)"
        write(unit=io,fmt="(a,3f12.6,a)")     "                                   \",finalm(3,:),"/"
      end if
      return
    End Subroutine Write_New_Monoc_Cell

  End Program Get_Conventional_Unit_Cells
