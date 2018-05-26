Submodule (CFML_Crystal_Metrics) Niggli_Cell 

 Contains
    !!--++
    !!--++ SUBROUTINE NIGGLI_CELL_ABC
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the Niggli cell when the input is the list of cell parameters
    !!--++    provided as a 6D vector. Calls the subroutine Niggli_Cell_Nigglimat for
    !!--++    the effective calculations
    !!--++
    !!--++ Update: October - 2008
    !!
    Module Subroutine Niggli_Cell_ABC(CellV,Niggli_Point,Cell,Trans)    
       !---- Arguments ----!
       real(kind=cp),dimension(6),             intent(in out) :: CellV         ! Cell parameters in a vector
       real(kind=cp),dimension(5),   optional, intent(out)    :: Niggli_Point  ! Niggli points
       class(CrysCell_Type),         optional, intent(out)    :: cell          ! Cell Object
       real(kind=cp), dimension(3,3),optional, intent(out)    :: trans         ! Transformation matrix

       !---- Local variables ----!
       real(kind=cp), dimension(2,3)    :: n_mat
       type(CrysCell_M_Type)            :: celda

       n_mat(1,1)=cellV(1)*cellV(1)
       n_mat(1,2)=cellV(2)*cellV(2)
       n_mat(1,3)=cellV(3)*cellV(3)
       n_mat(2,1)=cellV(2)*cellV(3)*cosd(cellV(4))
       n_mat(2,2)=cellV(1)*cellV(3)*cosd(cellV(5))
       n_mat(2,3)=cellV(1)*cellV(2)*cosd(cellV(6))

       if (present(Niggli_Point)) then
          if (present(trans)) then
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda,trans)
          else
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda)
          end if
       else if(present(trans)) then
          call Niggli_Cell_nigglimat(n_mat,cell=celda,trans=trans)
       else
          call Niggli_Cell_nigglimat(n_mat,cell=celda)
       end if
       if (Err_CFML%state) return
       
       !> Reconstruct the new cell (Niggli Cell)
       CellV(1) = sqrt(n_mat(1,1))
       CellV(2) = sqrt(n_mat(1,2))
       CellV(3) = sqrt(n_mat(1,3))
       CellV(4) = acosd(n_mat(2,1)/(cellV(2)*cellV(3)))
       CellV(5) = acosd(n_mat(2,2)/(cellV(1)*cellV(3)))
       CellV(6) = acosd(n_mat(2,3)/(cellV(1)*cellV(2)))
       
       !> The exactly type of cell is unknown
       if (present(cell)) then
          call Set_Crystal_Cell(cellV(1:3),cellV(4:6), Cell)
       end if 

       return
    End Subroutine Niggli_Cell_abc
    
    !!--++
    !!--++ SUBROUTINE NIGGLI_CELL_NIGGLIMAT
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the Niggli cell when the input is the Niggli Matrix (part of the metrics)
    !!--++    of a primitive cell. Applies the scalar algorithm of
    !!--++    I. Krivy and B. Gruber, Acta Cryst A32, 297 (1976)
    !!--++    If Trans is present, Celln should also be present.
    !!--++
    !!--++ Update: January - 2011
    !!
    Module Subroutine Niggli_Cell_Nigglimat(N_Mat,Niggli_Point,Cell,Trans)    
       !---- Arguments ----!
       real(kind=cp),dimension(2,3),              intent(in out) :: n_mat
       real(kind=cp),dimension(5),      optional, intent(out)    :: Niggli_Point
       class(CrysCell_Type),            optional, intent(out)    :: cell
       real(kind=cp), dimension(3,3),   optional, intent(out)    :: trans

       !--- Local variables ---!
       integer, parameter            :: NUMITER=100
       type(CrysCell_M_Type)         :: Cell1,Cell2
       real(kind=cp)                 :: A,B,C,u,v,w,eps
       real(kind=cp), dimension(3,3) :: trm
       real(kind=cp), dimension(3)   :: cel,ang
       integer                       :: iu,iv,iw, ncount ! ncount is the counter no more that Numiter=100
                                                         ! iterations are permitted. In case of exhausting
                                                         ! the iteration Err_Crys=.true. but the current
                                                         ! cell is output anyway
       real(kind=cp),parameter        :: epr=0.0001      ! Relative epsilon
       logical                        :: ok

       
       !> Init
       call clear_error()
       
       ! N is a Niggli cell of L if  (i) it is as Buerger cell of L and
       !                            (ii) |90-alpha| + |90-beta| + |90-gamma| -> maximum
       !                  / a.a  b.b  c.c \       /  s11  s22  s33 \
       !   Niggli matrix  |               |   =   |                |
       !                  \ b.c  a.c  a.b /       \  s23  s13  s12 /
       !
       ! I. Krivy and B. Gruber, Acta Cryst A32, 297 (1976)
       ! Krivy-Gruber algorithms safely implemented (suggestion of Ralf Grosse-Kunsleve)
       ! R.W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, Acta Cryst A60, 1-6 (2004)
       ! Epsilon: e
       !    x < y -> x < y-e;    x > y -> y < x-e
       !   x <= y -> .not. y < x-e;   x >= y -> .not. x < y-e
       !   x == y -> .not. (x < y-e .or. y < x-e)
       !
       A=n_mat(1,1)
       B=n_mat(1,2)
       C=n_mat(1,3)
       u=2.0*n_mat(2,1)
       v=2.0*n_mat(2,2)
       w=2.0*n_mat(2,3)
       eps=epr*(A*B*C)**(1.0/6.0)
       
       ncount=0
       ok=.true.
       if (present(trans)) then
          !> Construct the input cell Cellp from its Niggli parameters
          cel(1) = sqrt(A)
          cel(2) = sqrt(B)
          cel(3) = sqrt(C)
          ang(1) = acosd(u/(cel(2)*cel(3)*2.0))
          ang(2) = acosd(v/(cel(1)*cel(3)*2.0))
          ang(3) = acosd(w/(cel(1)*cel(2)*2.0))
          call Set_Crystal_Cell(cel,ang, Cell1)
       end if

       do
          ncount=ncount+1
          if (ncount > NUMITER) then
             ok=.false.
             exit
          end if

          !---- if(A > B .or. ( A == B  .and. abs(u) > abs(v)) ) then  ! A1
          if (B < A-eps .or. ( .not.( A < B-eps .or. B < A-eps)  .and. abs(v) < abs(u)-eps ) ) then  ! A1
             call swap(A,B)
             call swap(u,v)
          end if

          !---- if(B > C .or. ( B == C .and. abs(v) > abs(w)) ) then  ! A2
          if (C < B-eps .or. ( .not.( C < B-eps .or. B < C-eps) .and. abs(w) < abs(v)-eps) ) then  ! A2
             call swap(B,C)
             call swap(v,w)
             cycle
          end if

          !---- if (u*v*w > 0.0) then                                 ! A3
          iu=1; iv=1; iw=1
          if ( u < -eps) iu=-1
          if ( v < -eps) iv=-1
          if ( w < -eps) iw=-1
          if (abs(u) < eps) iu=0
          if (abs(v) < eps) iv=0
          if (abs(w) < eps) iw=0
          if (iu*iv*iw > 0) then                                      ! A3
             u=abs(u)
             v=abs(v)
             w=abs(w)
          else                                                        ! A4
             u=-abs(u)
             v=-abs(v)
             w=-abs(w)
          end if

          !---- if( abs(u) > B .or. ( u == B .and. 2.0*v < w) .or. ( u == -B .and. w < 0.0)) then  ! A5
          if ( B < abs(u)-eps  .or. ( .not.(u < B-eps .or. B < u-eps) .and. 2.0*v < w-eps) .or. &
             ( .not.(u < -B-eps .or. -B < u-eps) .and. w < -eps)) then  ! A5
             iu=1; if( u < -eps) iu=-1
             C = B+C - u * iu
             v =  v  - w * iu
             u = u - 2.0*B*iu
             cycle
          end if

          !---- if( abs(v) > A .or. ( v == A .and. 2.0*u < w) .or. ( v == -A .and. w < 0.0)) then  ! A6
          if ( A < abs(v)-eps .or. (.not. (v < A-eps .or. A < v-eps) .and. 2.0*u < w-eps) .or. &
             ( .not.( v < -A-eps .or. -A < v-eps) .and. w < -eps)) then  ! A6
             iv=1; if( v < -eps) iv=-1
             C = A+C - v * iv
             u =  u  - w * iv
             v = v - 2.0*A*iv
             cycle
          end if

          !---- if( abs(w) > A .or. ( w == A .and. 2.0*u < v) .or. ( w == -A .and. v < 0.0)) then  ! A7
          if ( A < abs(w)-eps .or. ( .not. (w < A-eps .or. A < w-eps) .and. 2.0*u < v-eps) .or. &
             ( .not. (w < -A-eps .or. -A < w-eps) .and. v < -eps)) then  ! A7
             iw=1; if( w < -eps) iw=-1
             B = A+B - w * iw
             u =  u  - v * iw
             w = w - 2.0*A*iw
             cycle
          end if

          !---- if(u+v+w+A+B < 0.0 .or. (u+v+w+A+B == 0.0 .and. 2.0*(A+v)+w > 0.0 )) then  ! A8
          if (u+v+w+A+B < -eps .or. ( abs(u+v+w+A+B) < eps .and. 2.0*(A+v)+w > eps )) then  ! A8
             C=A+B+C+u+v+w
             u=2.0*B+u+w
             v=2.0*A+v+w
             cycle
          end if
          exit
       end do

       !> Reconstruct the new Niggli matrix
       n_mat(1,1)=A; n_mat(1,2)=B; n_mat(1,3)=C
       n_mat(2,1)=0.5*u; n_mat(2,2)=0.5*v; n_mat(2,3)=0.5*w

       if (.not. ok) Then
          Err_CFML%state=.true.
          Err_CFML%Flag=2         
          Err_CFML%Msg=" The limit of iterations in Niggli_Cell_NiggliMat has been reached!"
          return
       end if

       if (present(Niggli_point)) then
          Niggli_point(1)= A/C
          Niggli_point(2)= B/C
          Niggli_point(3)= u/C
          Niggli_point(4)= v/C
          Niggli_point(5)= w/C
       end if

       if (present(cell)) then
          !Reconstruct the new cell (Niggli Cell)
          cel(1) = sqrt(A)
          cel(2) = sqrt(B)
          cel(3) = sqrt(C)
          ang(1) = acosd(u/(cel(2)*cel(3)*2.0))
          ang(2) = acosd(v/(cel(1)*cel(3)*2.0))
          ang(3) = acosd(w/(cel(1)*cel(2)*2.0))
          call Set_Crystal_Cell(cel,ang, Cell2)
          
          if (present(trans)) then
             call Get_Transfm_Matrix(cell1,cell2,trm)
             if (.not. err_CFML%state) then
                trans=trm
             else
                trans=identity
             end if
          end if
       end if

       !> The type of the cell is unknown
       if (present(cell)) then
          cel=cell2%cell
          ang=cell2%ang 
          call Set_Crystal_Cell(cel,ang, Cell)
       end if

       return
    End Subroutine Niggli_Cell_nigglimat
    
    !!--++
    !!--++ SUBROUTINE NIGGLI_CELL_PARAMS
    !!--++
    !!--++    (OVERLOAD)
    !!--++     Calculates the Niggli cell when the input is the list of cell parameters
    !!--++     provided as six scalars.
    !!--++     Calls the subroutine Niggli_Cell_Nigglimat for the effective calculations
    !!--++
    !!--++ Update: October - 2008
    !!
    Module Subroutine Niggli_Cell_Params(A,B,C,Al,Be,Ga,Niggli_Point,Cell,Trans)
       !---- Arguments ----!
       real(kind=cp),                           intent (in out)  :: a,b,c,al,be,ga
       real(kind=cp),dimension(5),    optional, intent(out)      :: Niggli_Point
       class(CrysCell_Type),          optional, intent(out)      :: cell
       real(kind=cp), dimension(3,3), optional, intent(out)      :: trans

       !--- Local variables ---!
       type(CrysCell_M_Type)            :: celda
       real(kind=cp), dimension(2,3)    :: n_mat

       !> Init
       call clear_error()
       
       if ( al+be < ga+1.0  .or. al+ga < be+1.0 .or. be+ga < al+1.0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" The provided angles cannot set a unit cell!"
          return
       end if

       call Set_Crystal_Cell((/a,b,c/),(/al,be,ga/), Celda)
       if (Err_CFML%state) return

       n_mat(1,1)=Celda%GD(1,1); n_mat(1,2)=Celda%GD(2,2); n_mat(1,3)=Celda%GD(3,3)
       n_mat(2,1)=Celda%GD(2,3); n_mat(2,2)=Celda%GD(1,3); n_mat(2,3)=Celda%GD(1,2)

       if (present(Niggli_Point)) then
          if (present(trans)) then
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda,trans)
          else
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda)
          end if
       else if(present(trans)) then
          call Niggli_Cell_nigglimat(n_mat,cell=celda,trans=trans)
       else
          call Niggli_Cell_nigglimat(n_mat,cell=celda)
       end if
       if (Err_CFML%state) return

       !> Final values
       a=celda%cell(1); b=celda%cell(2); c=celda%cell(3)
       al=celda%ang(1); be=celda%ang(2); ga=celda%ang(3)
       
       !> we don't know which type is cell variable
       if (present(cell)) then
          call Set_Crystal_Cell([a,b,c],[al,be,ga], Cell)
       end if

       return
    End Subroutine Niggli_Cell_Params

    !!--++
    !!--++ SUBROUTINE NIGGLI_CELL_TYPE
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the Niggli cell when the input is an object of type Crystal_Cell_Type
    !!--++    Calls the subroutine Niggli_Cell_Nigglimat for the effective calculations
    !!--++
    !!--++ Update: October - 2008
    !!
    Module Subroutine Niggli_Cell_Type(Cell,Niggli_Point,Celln,Trans)
       !---- Arguments ----!
       class(CrysCell_M_Type),                  intent(in out ) :: cell
       real(kind=cp),dimension(5),    optional, intent(out)     :: Niggli_Point
       class(CrysCell_Type),          optional, intent(out)     :: celln
       real(kind=cp), dimension(3,3), optional, intent(out)     :: trans

       !--- Local variables ---!
       type(CrysCell_M_Type)           :: celda
       real(kind=cp), dimension(2,3)   :: n_mat

       !> Init
       call clear_error()
       
       call Set_Crystal_Cell(cell%cell,cell%ang,celda)
       n_mat(1,1)=Celda%GD(1,1); n_mat(1,2)=Celda%GD(2,2); n_mat(1,3)=Celda%GD(3,3)
       n_mat(2,1)=Celda%GD(2,3); n_mat(2,2)=Celda%GD(1,3); n_mat(2,3)=Celda%GD(1,2)

       if (present(Niggli_Point)) then
          if (present(trans)) then
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda,trans)
          else
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda)
          end if
          
       else if(present(trans)) then
          call Niggli_Cell_nigglimat(n_mat,cell=celda,trans=trans)
       else
          call Niggli_Cell_nigglimat(n_mat,cell=celda)
       end if
       if (Err_CFML%state) return

       call Set_Crystal_Cell(celda%cell,celda%ang,cell)
       
       !> Problem with the type of celln
       if (present(celln)) then
          call Set_Crystal_Cell(celda%cell,celda%ang,celln)
       end if

       return
    End Subroutine Niggli_Cell_Type

    !!--++
    !!--++ SUBROUTINE NIGGLI_CELL_VECT
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the Niggli cell when the input is given as three vectors
    !!--++    in Cartesian components. A test of linear indenpendency is performed.
    !!--++    Calls the subroutine Niggli_Cell_Nigglimat for the effective calculations
    !!--++
    !!--++ Update: October - 2008
    !!
    Module Subroutine Niggli_Cell_Vect(A,B,C,Niggli_Point,Cell,Trans)
       !---- Arguments ----!
       real(kind=cp),dimension(3),                intent(in)     :: a,b,c
       real(kind=cp),dimension(5),      optional, intent(out)    :: Niggli_Point
       class(CrysCell_Type),            optional, intent(out)    :: cell
       real(kind=cp), dimension(3,3),   optional, intent(out)    :: trans

       !--- Local variables ---!
       type(CrysCell_M_Type)         :: celda
       real(kind=cp), dimension(2,3) :: n_mat
       real(kind=cp)                 :: det

       !> Init
       call clear_error()
       
       det=determ_Vec(a,b,c)
       if (abs(det) < 0.0001) then
          Err_CFML%state=.true.
          Err_CFML%flag=2
          ERR_CFML%Msg=" The three input vectors are nor linearly independent!"
          return
       end if
       n_mat(1,1)=dot_product(a,a); n_mat(1,2)=dot_product(b,b); n_mat(1,3)=dot_product(c,c)
       n_mat(2,1)=dot_product(b,c); n_mat(2,2)=dot_product(a,c); n_mat(2,3)=dot_product(a,b)

       if (present(Niggli_Point)) then
          if (present(trans)) then
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda,trans)
          else
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda)
          end if
          
       else if(present(trans)) then
          call Niggli_Cell_nigglimat(n_mat,cell=celda,trans=trans)
          
       else
          call Niggli_Cell_nigglimat(n_mat,cell=celda)
       end if
       if (Err_CFML%state) return
       
       !> We don't know which type of cell was used.
       if (present(cell)) then
          call Set_Crystal_Cell(celda%cell,celda%ang, Cell) 
       end if 

       return
    End Subroutine Niggli_Cell_Vect

End Submodule Niggli_Cell 