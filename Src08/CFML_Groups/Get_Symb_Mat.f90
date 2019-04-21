!!----
!!----
!!----
!!
SubModule (CFML_Groups) CFML_GRP_012
   Contains
   
   !!----
   !!---- GET_SYMB_FROM_MAT
   !!----
   !!---- 19/04/19
   !!
   Module Function Get_Symb_from_Mat(Mat, Strcode, Invt) Result(Symb)
      !---- Arguments ----!
      type(rational), dimension(:,:), intent(in) :: Mat
      character(len=*), optional,     intent(in) :: strcode
      integer,          optional,     intent(in) :: invt
      character(len=80)                          :: symb

      !---- Local Variables ----!
      character(len=3),dimension(10)    :: x_typ
      character(len=15)                 :: car
      character(len=40)                 :: translation
      character(len=15),dimension(10)   :: sym
      integer                           :: i,j,Dd,d,k
      logical                           :: abc_type

      Dd=size(Mat,dim=1)
      d=Dd-1
      x_typ=xyz
      abc_type=.false.
      if (present(strcode)) then
         Select Case (trim(strcode))
            case("xyz")
               x_typ=xyz
               
            case("x1x2x3")
               x_typ=x1x2x3
               
            case("abc")
               x_typ=abc
               abc_type=.true.
          
            case Default
          	   x_typ=xyz
         End Select
      end if
      
      !> Main
      symb=" "
      translation=" "
      do i=1, d
         sym(i)=" "
         do j=1,d
            if (Mat(i,j) == 1_LI) then
               sym(i) = trim(sym(i))//"+"//trim(x_typ(j))
            
            else if(Mat(i,j) == -1_LI) then
               sym(i) =  trim(sym(i))//"-"//trim(x_typ(j))
             
            else if(Mat(i,j) /= 0_LI) then
               car=adjustl(Rational_String(Mat(i,j)))
               k=index(car,"/")
               if (k /= 0) then
                  if (car(1:1) == "1") then
                     car=trim(x_typ(j))//car(k:)
                     
                  else if(car(1:2) == "-1") then
                     car="-"//trim(x_typ(j))//car(k:)
                     
                  else
                     car=car(1:k-1)//trim(x_typ(j))//car(k:)
                  end if
               
               else
                  car=trim(car)//trim(x_typ(j))
               end if
               
               if (Mat(i,j) > 0_LI) car="+"//trim(car)
               sym(i)=trim(sym(i))//pack_string(car)
            end if
         end do

         !> Write here the translational part for each component
         if (Mat(i,Dd) /= 0_LI) then
            car=adjustl(Rational_String(Mat(i,Dd)))
            if (abc_type) then
               translation=trim(translation)//","//trim(car)
            else
               if (car(1:1) == "-") then
                  sym(i)=trim(sym(i))//trim(car)
               else
                  sym(i)=trim(sym(i))//"+"//trim(car)
               end if
            end if
         
         else
            if (abc_type) translation=trim(translation)//",0"
         end if
         
         sym(i)=adjustl(sym(i))
         if (sym(i)(1:1) == "+")  then
            sym(i)(1:1) = " "
            sym(i)=adjustl(sym(i))
         end if
         sym(i)=pack_string(sym(i))
      end do
       
      symb=sym(1)
      do i=2,d
         symb=trim(symb)//","//trim(sym(i))
      end do
      
      if (abc_type)then
         symb=trim(symb)//";"//trim(translation(2:))
       
      else
         if (present(invt)) then
            write(unit=car,fmt="(i2)") invt !Rational_String(Mat(Dd,Dd))
            car=adjustl(car)
            symb=trim(symb)//","//trim(car)
         end if
      end if
      
      return
   End Function Get_Symb_from_Mat

End SubModule CFML_GRP_012   
   
