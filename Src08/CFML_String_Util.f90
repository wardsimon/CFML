!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2018  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----               Ross Angel         (University of Pavia) 
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!----
!!---- MODULE: CFML_String_Utilities
!!----   INFO: Manipulation of strings with alfanumeric characters
!!----
!!----
!!
 Module CFML_String_Utilities
    !---- Use Modules ----!
    use CFML_GlobalDeps,   only: cp, ops_sep, err_cfml,err_cfml_flag,err_cfml_msg, clear_error
    use CFML_Math_General, only: Negligible, Zbelong
    use ieee_arithmetic,   only: ieee_is_nan,ieee_is_finite

    implicit none

    private

    !---- List of public functions ----!
    public :: Equal_Sets_Text, L_Case, Pack_String, U_Case, Strip_String, String_Count

    !---- List of public subroutines ----!
    public :: Cutst, Format_Real, Frac_Trans_1Dig, Frac_Trans_2Dig, Get_DateTime, Get_Dirname, &
              Get_Extension, Get_Filename, Get_Fraction_1Dig, Get_Fraction_2Dig,               &
              Get_Mat_From_Symb, Get_Separator_Pos, Get_Substring_Positions, Get_Symb_From_Mat,& 
              Get_Transf, Get_Vec_From_String, Get_num, Get_Numstd, get_word, Inc_LineNum,     &
              Init_FindFmt, Lcase, Number_lines, NumCol_from_NumFmt, Read_Key_Str,             &
              Read_Key_strVal, Read_Key_Value, Read_Key_ValueSTD, Reading_Lines, Set_numstd,   &
              Substring_Replace, Ucase, FindFMT


    !> Parameters
    character (len=*), parameter :: DIGIT         ="0123456789.-"   ! Character parameter for numbers
        
    !> TYPES
    !!----
    !!---- TYPE :: ERR_TEXT_TYPE
    !!--..
    !!
    Type, Public :: Err_Text_Type
       integer                           :: nlines=0
       character (len=132), dimension(5) :: txt=" "
    End Type Err_Text_Type

    !> Variables
    integer,              public :: iErr_fmt      ! Error code value (should be normally = 0) for findFMT procedure
    Type (Err_Text_Type), public :: Mess_FindFMT  ! Text composed (findFMT) 

    Interface
       Module Pure Function Equal_Sets_Text(str1,n1,str2,n2) result(Equal_sets_texto)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in) :: str1              ! Vector of String
          character(len=*), dimension(:), intent(in) :: str2              ! Vector of String
          integer,                        intent(in) :: n1                ! Number of lines on Text1
          integer,                        intent(in) :: n2                ! Number of lines on str2
          logical                                    :: Equal_sets_texto  !
       End Function Equal_Sets_Text
       
       Module Pure Function L_Case(Str) Result (LStr)
          !---- Argument ----!
          character (len=*), intent(in) :: Str    ! Input String
          character (len=len(str))      :: LStr   ! lower case of Text
       End Function L_Case   
       
       Module Pure Function Pack_String(Str) Result (Strp)
          !---- Argument ----!
          character (len=*), intent(in) :: str    ! Input String
          character (len=len_trim(str)) :: strp   ! Output string
       End Function Pack_String   
       
       Module Pure Function String_Count(str,substr) result(coun)
          !---- Arguments ----!
          character(len=*), intent(in) :: str       ! Input String
          character(len=*), intent(in) :: substr    ! Substring model
          integer                      :: coun      ! Number
       End Function String_Count   
       
       Module Pure Function Strip_String(str, to_strip) Result(striped_string)
          !---- Arguments----!
          character (len = *), intent(in) :: str            ! Input string
          character (len = *), intent(in) :: to_strip       ! Pattern
          character (len = len_trim(str)) :: striped_string
       End Function Strip_String   
       
       Module Pure Function U_Case(Str) Result (UStr)
          !---- Argument ----!
          character (len=*), intent(in) :: Str   ! Input string
          character (len=len(Str))      :: UStr  ! Upper conversion
       End Function U_Case   
        
       Module Subroutine BuildFMT(iFld,nCar,nStr,FMTstring)
          !---- Arguments ----!
          Integer,           intent(in    ) ::   iFld       ! Format type
          Integer,           intent(in out) ::   nCar       ! integer/real field: number of characters in field
                                                            ! character field: number of characters to skip before A field
          Integer,           intent(in out) ::   nStr       ! current character number in FMTstring
          Character (len=*) ,intent(in out) ::   FMTstring  ! FORTRAN format string
       End Subroutine BuildFMT
       
       Module Subroutine FindFmt(IUnit,aLine,FMTfields,FMTstring,idebug)
          !---- Arguments ----!
          Integer ,           intent(in    ) ::  IUnit      ! Logical unit number
          Character (len=*) , intent(in out) ::  aLine      ! character string to be decoded
          Character (len=*) , intent(in    ) ::  FMTfields  ! description of the format fields (e.g. IIFIF)
          Character (len=*) , intent(   out) ::  FMTstring  ! format of the line (e.g. (I5,I1,F8.0,I4,F7.0,)
          Integer ,optional,  intent(in    ) ::  idebug     ! Logical unit number for writing the input file
       End Subroutine FindFmt
       
       Module Subroutine FindFMT_Err(aLine,nC_L)
          !---- Arguments ----!
          Character(len=*), intent(in) ::   aLine      ! Current data line                   
          Integer,         intent (in) ::   nC_L       ! location of last character treated  
       End Subroutine FindFMT_Err
       
       Module Subroutine Inc_LineNum(line_n)
          !---- Argument ----!
          integer, intent(in) :: line_n
       End Subroutine Inc_LineNum  
       
       Module Subroutine Init_FindFMT(nline)
          !---- Arguments ----!
          integer, optional, intent(in) :: nline
       End Subroutine Init_FindFMT   
       
       Module Subroutine SGetFTMfield(GetFTMfield,FMTfields,nFld,nFldMax)
          !---- Arguments ----!
          Integer ,          intent(out)    ::  GetFTMfield
          Character (len=*) ,intent( in)    ::  FMTfields        !  -> format descriptor
          Integer ,          intent(in out) ::  nFld             ! <-> current field in format descriptor
          Integer ,          intent( in)    ::  nFldMax          !  -> max. number of fields in format descriptor
       End Subroutine SGetFTMfield
      
       Module Subroutine TreatMCharField(iFld,aLine,L_Line,nC_L,nC_X)
          !---- Arguments ----!
          Integer,           intent(in out)  :: iFld      ! <-> "A" format size (1 to 9)
          Character (len=*), intent(in)      :: aLine     !  -> data line to be analysed
          Integer,           intent(in)      :: L_Line    !  -> true length of data Line
          Integer,           intent(in out)  :: nC_L      ! <-> current character in data line
          Integer,           intent(out)     :: nC_X      ! <-  number of characters in X format field (now nx -> trn)
       End Subroutine TreatMCharField
 
       Module Subroutine TreatNumerField(iFld,aLine,L_Line,nC_L,nCar)
          !---- Arguments ----!
          Integer ,          intent( in)    ::  iFld   ! field type
          Character (len=*), intent(in out) ::  aLine  ! data line
          Integer ,          intent( in)    ::  L_Line ! true length of the data line
          Integer ,          intent(in out) ::  nC_L   ! counts characters in data line
          Integer ,          intent(in out) ::  nCar   ! counts characters in format field
       End Subroutine TreatNumerField
       
       Module Subroutine Cutst(Str1,nlong1,Str2,nlong2)
          !---- Argument ----!
          character (len=*),           intent(in out) :: Str1     ! Input string / Out: string without the first word
          character (len=*), optional, intent(   out) :: Str2     ! The first word of String on Input
          integer,           optional, intent(   out) :: nlong1   ! Give the length of Str1 on Output
          integer,           optional, intent(   out) :: nlong2   ! Give the length of Str2 on Output
       End Subroutine Cutst
       
       Module Subroutine Get_DateTime(Str)
          !---- Argument ----!
          character(len=*), intent(out) :: Str  ! String containing the Date and Time
       End Subroutine Get_DateTime
       
       Module Subroutine Get_Dirname(Str,Directory)
          !---- Argument ----!
          Character (Len=*), Intent (In)  :: Str          ! String containing Path + Filename
          Character (Len=*), Intent (Out) :: Directory    ! Path
       End Subroutine Get_Dirname
       
       Module Subroutine Get_Extension(filename, extension, dotted)
          !---- Arguments ----!
          character(len=*),  intent(in)  :: filename     ! Input filename
          character(len=*),  intent(out) :: extension    ! Extension of the file
          logical, optional, intent(in)  :: dotted       ! If True, the extension will be returned with a dot
       End Subroutine Get_Extension
       
       Module Subroutine Get_Filename(Str,Filename)
          !---- Argument ----!
          Character (Len=*), Intent (In)  :: Str       ! String containing Path + Filename
          Character (Len=*), Intent (Out) :: Filename  ! Filename
       End Subroutine Get_Filename  
       
       Module Subroutine Get_Separator_Pos(Str,car,pos,ncar)
          !---- Arguments ----!
          character(len=*),      intent(in)  :: Str   ! Inout String
          character(len=1),      intent(in)  :: car   ! Separator character
          integer, dimension(:), intent(out) :: pos   ! Vector with positions of "sep" in "Line"
          integer,               intent(out) :: ncar  ! Number of appearance of "sep" in "Line" 
       End Subroutine Get_Separator_Pos 
       
       Module Subroutine Get_Substring_Positions(str,substr,pos,nsubs)
          !---- Arguments ----!
          character(len=*),      intent(in)  :: str         ! In -> Input String                                                         
          character(len=*),      intent(in)  :: substr      ! In -> Substring                                                            
          integer, dimension(:), intent(out) :: pos         ! Out -> Vector with positions of the firs character of "substr" in "String" 
          integer,               intent(out) :: nsubs       ! Out -> Number of appearance of "substr" in "String"                        
       End Subroutine Get_Substring_Positions 
       
       Module Subroutine Lcase(line)
          !---- Argument ----!
          character (len=*), intent(in out) :: line
       End Subroutine Lcase 
       
       Module Subroutine SubString_Replace(string, substr, repstr, warning)
          !---- Arguments ----!
          character(len=*), intent(in out) :: string   ! Input/output string
          character(len=*), intent(in)     :: substr   ! Subtring to be replaced
          character(len=*), intent(in)     :: repstr   ! String for add
          character(len=*), intent(out)    :: warning  ! Message 
       End Subroutine SubString_Replace 
       
       Module Subroutine Ucase(Str)
          !---- Argument ----!
          character (len=*), intent(in out) :: Str  
       End Subroutine Ucase  
       
       Module Subroutine Format_Real(Val,W,fmtcar) 
          !---- Arguments ----!
          real,             intent(in)  :: val        ! value to be output
          integer,          intent(in)  :: w          ! Width
          character(len=*), intent(out) :: fmtcar 
       End Subroutine Format_Real
       
       Module Subroutine Frac_Trans_1Dig(Vec,Str)
          !---- Argument ----!
          real(kind=cp), dimension(3), intent( in)   :: Vec  ! Vector
          character (len=* ),          intent(out)   :: Str  ! String with conversion to fractional
       End Subroutine Frac_Trans_1Dig 
       
       Module Subroutine Frac_Trans_2Dig(Vec,Str)
          !---- Argument ----!
          real(kind=cp), dimension(3), intent( in) :: Vec   ! Vector
          character (len=* ),          intent(out) :: Str   ! String with conversion to fractional  
       End Subroutine Frac_Trans_2Dig
       
       Module Subroutine Get_Fraction_1Dig(V,Str)
          !---- Argument ----!
          real(kind=cp),    intent( in) :: V   !  Real value
          character(len=*), intent(out) :: Str !  Fracction in character form
       End Subroutine Get_Fraction_1Dig 
       
       Module Subroutine Get_Fraction_2Dig(V,Str)
          !---- Argument ----!
          real(kind=cp),    intent( in) :: V   !  Real value
          character(len=*), intent(out) :: Str !  Fracction in character form
       End Subroutine Get_Fraction_2Dig 
       
       Module Subroutine Get_Mat_From_Symb(Symb,Mat,cod)
          !---- Arguments ----!
          character(len=*),                intent(in)  :: Symb   ! String                               
          real(kind=cp),dimension(3,3),    intent(out) :: Mat    ! (/"u","v","w"/) or (/"x","y","z"/)   
          character(len=1), dimension(3),  intent(in)  :: cod    ! Output matrix  
       End Subroutine Get_Mat_From_Symb 
       
       Module Subroutine Get_Vec_From_String(Str,Cod, Vec)
          !---- Arguments ----!
          character(len=*),                intent(in)  :: str   ! Input string
          character(len=1), dimension(3),  intent(in)  :: cod   ! Code
          real(kind=cp),dimension(3),      intent(out) :: vec   ! Vector
       End Subroutine Get_Vec_From_String
       
       Module Subroutine Get_Symb_From_Mat(Mat,Symb,cod)
          !---- Arguments ----! 
          real(kind=cp),dimension(3,3),    intent(in)  :: Mat    ! Array
          character(len=*),                intent(out) :: Symb   ! Symbol
          character(len=1), dimension(3),  intent(in)  :: cod    ! Codes
       End Subroutine Get_Symb_From_Mat
       
       Module Subroutine Get_Transf(str,mat,v,cod)
          !---- Arguments ----!
          character(len=*),                          intent(in)  :: str      ! Input string
          real(kind=cp),dimension(3,3),              intent(out) :: mat      ! Matrix      
          real(kind=cp),dimension(3),                intent(out) :: v        ! Vector      
          character(len=1), dimension(4), optional,  intent(in)  :: cod      ! Code
       End Subroutine Get_Transf  
       
       Module Subroutine Get_Num(Str,vet,ivet,iv)
          !---- Argument ----!
          character (len=*),          intent ( in) :: Str   ! Input String to convert
          real(kind=cp), dimension(:),intent (out) :: vet   ! Vector of real numbers
          integer, dimension(:),      intent (out) :: ivet  ! Vector of integer numbers
          integer,                    intent (out) :: iv    ! Number of numbers in Vet/Ivet
       End Subroutine Get_Num
       
       Module Subroutine Get_NumStd(Str, value, std, ic)
          !----Arguments ----!
          character(len=*),             intent( in) :: Str     ! Input String
          real(kind=cp), dimension(:),  intent(out) :: value   ! Vector of values with real numbers
          real(kind=cp), dimension(:),  intent(out) :: std     ! Vector of standard deviation values
          integer,                      intent(out) :: ic      ! Number of components of vector Value
       End Subroutine Get_NumStd
       
       Module Subroutine Get_Word(Str,dire,ic)
          !---- Argument ----!
          character (len=*),                 intent ( in) :: Str   ! Input string
          character (len=*), dimension(:),   intent (out) :: dire  ! Vector of Words
          integer,                           intent (out) :: ic    ! Number of words
       End Subroutine Get_Word
       
       Module Subroutine NumCol_from_NumFmt(Str,n_col)
          !---- Argument ----!
          character (len=*), intent(in)  :: Str    ! Input format string
          Integer,           intent(out) :: n_col  ! Integer number of columns
       End Subroutine NumCol_from_NumFmt
       
       Module Subroutine Read_Fract(str,value)
          !---- Arguments ----!
          Character(len=*), intent(in) :: str     ! Input String
          real(kind=cp),    intent(out):: value   ! Value
       End Subroutine Read_Fract
       
       Module Subroutine Set_NumStd(Value, Std, Str)
          !---- Argument ----!
          real(kind=cp),   intent(in)  :: Value    ! Value
          real(kind=cp),   intent(in)  :: Std      ! Standard deviation
          character(len=*),intent(out) :: Str      ! String containing the information
       End Subroutine Set_NumStd
       
       Module Subroutine Read_Key_Str(filevar,nline_ini,nline_end,keyword,string,comment)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)      :: filevar     ! Input vector of String
          integer,                        intent(in out)  :: nline_ini   ! Pointer to initial position to search
                                                                         ! Out -> Pointer to final position in search
          integer,                        intent(in)      :: nline_end   ! Pointer to final position to search
          character(len=*),               intent(in)      :: keyword     ! Word to search
          character(len=*),               intent(out)     :: string      ! Rest of the input string
          character(len=1), optional,     intent(in)      :: comment     ! Character that define a comment line
       End Subroutine Read_Key_Str
       
       Module Subroutine Read_Key_StrVal(filevar,nline_ini,nline_end,keyword,string,vet,ivet,iv,comment)
          !---- Arguments ----!
          character(len=*), dimension(:),           intent(in)      :: filevar       !  In -> Input vector of String
          integer,                                  intent(in out)  :: nline_ini     !  In -> Pointer to initial position to search
                                                                                     ! Out -> Pointer to final position in search
          integer,                                  intent(in)      :: nline_end     !  In -> Pointer to final position to search
          character(len=*),                         intent(in)      :: keyword       !  In -> Word to search
          character(len=*),                         intent(out)     :: string        ! Out -> Rest of the input string
          real(kind=cp),dimension(:),     optional, intent(out)     :: vet           ! Out -> Vector for real numbers
          integer,dimension(:),           optional, intent(out)     :: ivet          ! Out -> Vector for integer numbers
          integer,                        optional, intent(out)     :: iv            ! Out -> Number of numbers
          character(len=1),               optional, intent(in)      :: comment       ! Character that define a comment line
       End Subroutine Read_Key_StrVal
       
       Module Subroutine Read_Key_Value(filevar,nline_ini,nline_end,keyword,vet,ivet,iv,comment,line_key)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)     :: filevar            !  In -> Input vector of String
          integer,                        intent(in out) :: nline_ini          !  In -> Pointer to initial position to search
                                                                               ! Out -> Pointer to final position in search
          integer,                        intent(in)     :: nline_end          !  In -> Pointer to final position to search
          character(len=*),               intent(in)     :: keyword            !  In -> Word to search
          real(kind=cp),dimension(:),     intent(out)    :: vet                ! Out -> Vector for real numbers
          integer,dimension(:),           intent(out)    :: ivet               ! Out -> Vector for integer numbers
          integer,                        intent(out)    :: iv                 ! Out -> Number of components
          character(len=1),     optional, intent(in)     :: comment            ! Consider the character passed in comment as a comment to skip the line
          character(len=*),     optional, intent(out)    :: line_key           ! Out -> Cut line where keyword is read
       End Subroutine Read_Key_Value 
       
       Module Subroutine Read_Key_ValueSTD(filevar,nline_ini,nline_end,keyword,vet1,vet2,iv,comment)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: filevar         !  In -> Input vector of String
          integer,                         intent(in out) :: nline_ini       !  In -> Pointer to initial position to search
                                                                             ! Out -> Pointer to final position in search
          integer,                         intent(in)     :: nline_end       !  In -> Pointer to final position to search
          character(len=*),                intent(in)     :: keyword         !  In -> Word to search
          real(kind=cp),dimension(:),      intent(out)    :: vet1            ! Out -> Vector of real numbers
          real(kind=cp),dimension(:),      intent(out)    :: vet2            ! Out -> Vector of standard deviations
          integer,                         intent(out)    :: iv              ! Out -> Number of components
          character(len=1),      optional, intent(in)     :: comment         ! Consider the character passed in comment as a comment to skip the line
       End Subroutine Read_Key_ValueSTD
          
    End Interface 
     

 Contains

    !!----
    !!---- SUBROUTINE NUMBER_LINES
    !!----
    !!----    Return the number of lines contained in a file.
    !!----    The file will be opened and closed before returning to the calling unit.
    !!----    If 'input_string' is present, return the number of lines until 'input_string' is founded
    !!----    as first string in the line
    !!----    (example : input_string =='END' : avoid Q peaks in a SHELX file)
    !!----
    !!---- Update: March-2014
    !!
    Subroutine Number_Lines(filename,n, input_string)
       !---- Arguments ----!
       character(len=*),           intent(in)  :: filename      ! Nane of File
       integer,                    intent(out) :: n             ! Number of lines
       character(len=*), optional, intent(in)  :: input_string  ! String to Exit

       !---- Local Variables ----!
       logical            :: info,opn
       integer            :: lun,cond,olun
       character (len=256):: read_line
       integer            :: long

       !---- Init ----!
       n=0
       info=.false.
       cond=0

       long=0
       if (present(input_string)) long = len_trim(input_string)

       !---- Exist filename ? ----!
       inquire (file=trim(filename),exist=info)
       if (.not. info) return

       !> Check if the file is already opened
       inquire(file=trim(filename),opened=opn, number=olun)   
       if (opn) then
          rewind(olun)
          lun=olun
       else
          !> Open file
          open(newunit=lun,file=trim(filename), status="old",action="read", position="rewind")
       end if   

       !---- Counting lines ----!
       do
          read(unit=lun,fmt="(a)",iostat=cond) read_line
          if (cond /= 0) exit
          read_line=adjustl(read_line)

          if (present(input_string)) then
             if (u_case(read_line(1:long)) == u_case(input_string(1:long))) exit
          end if
          n=n+1
       end do

       if (.not. opn) then
          close(unit=lun)
       else
          rewind(unit=lun)
       end if

       return
    End Subroutine Number_Lines
    
    !!----
    !!---- SUBROUTINE READING_LINES
    !!----
    !!----    Read nlines of the file and put the information on Filevar. The file is opened to read the
    !!----    lines and closed before returning to the calling unit.
    !!----
    !!---- Update: February - 2005, March-2018 (reintroducing the "opened" inquire,JRC)
    !!
    Subroutine Reading_Lines(filename,nlines,filevar)
       !---- Arguments ----!
       character(len=*),               intent( in) :: filename   ! Filename
       integer,                        intent( in) :: nlines     ! Number of lines to be read
       character(len=*), dimension(:), intent(out) :: filevar    ! String vector

       !---- Local Variables ----!
       logical :: info,opn
       integer :: lun,i,olun

       !---- Init ----!
       call clear_error()
       info=.false.

       !---- Exist filename ? ----!
       inquire (file=trim(filename),exist=info)
       if (.not. info) then
          Err_CFML=.true.
          Err_CFML_Flag=2
          ERR_cfml_Msg="The file"//trim(filename)//" does not exist "
          return
       end if

       inquire(file=trim(filename),opened=opn, number=olun)   !Check if the file is already opened
       if (opn) then
          rewind(olun)
          lun=olun
       else
          open(newunit=lun,file=filename, status="old",action="read", position="rewind")
       end if

       !---- Reading... ----!
       do i=1,nlines
          read(unit=lun,fmt="(a)") filevar(i)
       end do

       if (.not. opn) close(unit=lun)

       return
    End Subroutine Reading_Lines
    
 End Module CFML_String_Utilities
