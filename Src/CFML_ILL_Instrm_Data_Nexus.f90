!-------------------------------------------------------
!---- Crystallographic Fortran Modules Library (CrysFML)
!-------------------------------------------------------
!---- The CrysFML project is distributed under LGPL. In agreement with the
!---- Intergovernmental Convention of the ILL, this software cannot be used
!---- in military applications.
!----
!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!----
!---- Authors: Juan Rodriguez-Carvajal (ILL)
!----          Javier Gonzalez-Platas  (ULL)
!----          Nebil Ayape Katcho      (ILL)
!----
!---- Contributors: Laurent Chapon     (ILL)
!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!----               Tierry Roisnel     (CDIFX,Rennes France)
!----               Eric Pellegrini    (ILL)
!----
!---- This library is free software; you can redistribute it and/or
!---- modify it under the terms of the GNU Lesser General Public
!---- License as published by the Free Software Foundation; either
!---- version 3.0 of the License, or (at your option) any later version.
!----
!---- This library is distributed in the hope that it will be useful,
!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!---- Lesser General Public License for more details.
!----
!---- You should have received a copy of the GNU Lesser General Public
!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!----
!----
!---- MODULE: CFML_ILL_Instrm_Data
!----   INFO: Subroutines related to Instrument information from ILL
!----
Module CFML_ILL_Instrm_Data_Nexus
    use, intrinsic :: iso_c_binding
    implicit none

#ifdef USE_HDF

    interface
        subroutine c_read_init_nxs (string_filename) bind ( c )
            use iso_c_binding, only: c_char
            character ( kind=c_char  ):: string_filename(*)
        end subroutine c_read_init_nxs
    end interface

    interface
        subroutine print_c(string) bind(C, name="print_C")
            use iso_c_binding, only: c_char
            character(kind=c_char) :: string(*)
        end subroutine print_c
    end interface

    interface
        subroutine c_read_header_nxs ( ) bind ( c )
            use iso_c_binding
        end subroutine c_read_header_nxs
    end interface

    interface
        subroutine c_get_header_numor_nxs (n) bind ( c )
            use iso_c_binding
            integer ( c_int ) :: n
        end subroutine c_get_header_numor_nxs
    end interface

    interface
        subroutine c_get_header_scantype_nxs (n)  bind ( c )
            use iso_c_binding
            integer ( c_int ) ::  n
        end subroutine c_get_header_scantype_nxs
    end interface

    interface
        Subroutine c_get_header_instr_name_nxs(fname,size1) bind ( c )
            use, intrinsic :: iso_c_binding
            character(kind=c_char)  :: fname(*)
            integer ( c_int ) ::size1
        End Subroutine c_get_header_instr_name_nxs
    end interface

    interface
        Subroutine c_get_header_subt_nxs(fname,size1) bind ( c )
            use, intrinsic :: iso_c_binding
            character(kind=c_char) :: fname(*)
            integer size1
        End Subroutine c_get_header_subt_nxs
    end interface

    interface
        subroutine c_read_integer_bloc_nxs () bind ( c )
            use iso_c_binding
        end subroutine c_read_integer_bloc_nxs
    end interface

    interface
        subroutine c_get_integer_bloc_nxs(tabi,n2) bind ( c )
            use iso_c_binding
            integer ( c_int ) :: tabi(*)
            integer ( c_int ), VALUE:: n2
        end subroutine c_get_integer_bloc_nxs
    end interface

    interface
        subroutine c_read_float_bloc_nxs () bind ( c )
            use iso_c_binding
        end subroutine c_read_float_bloc_nxs
    end interface

    interface
        subroutine c_get_float_bloc_nxs(tabf,n2) bind ( c )
            use iso_c_binding
            real ( c_float ) :: tabf(*)
            integer ( c_int ), VALUE:: n2
        end subroutine c_get_float_bloc_nxs
    end interface

    interface
        subroutine c_read_data_bloc_param_nxs () bind ( c )
            use iso_c_binding
        end subroutine c_read_data_bloc_param_nxs
    end interface

    interface
        subroutine c_get_data_bloc_time_nxs(tabf,n2) bind ( c )
            use iso_c_binding
            real ( c_float ) :: tabf(*)
            integer ( c_int ), VALUE:: n2
        end subroutine c_get_data_bloc_time_nxs
    end interface

    interface
        subroutine c_get_data_bloc_moni_nxs(tabf,n2) bind ( c )
            use iso_c_binding
            real ( c_float ) :: tabf(*)
            integer ( c_int ), VALUE:: n2
        end subroutine c_get_data_bloc_moni_nxs
    end interface

    interface
        subroutine c_get_data_bloc_total_count_nxs(tabf,n2) bind ( c )
            use iso_c_binding
            real ( c_float ) :: tabf(*)
            integer ( c_int ), VALUE:: n2
        end subroutine c_get_data_bloc_total_count_nxs
    end interface

    interface
        subroutine c_get_data_bloc_angle1_nxs(tabf,n2) bind ( c )
            use iso_c_binding
            real ( c_float ) :: tabf(*)
            integer ( c_int ), VALUE:: n2
        end subroutine c_get_data_bloc_angle1_nxs
    end interface

    interface
        subroutine c_get_data_bloc_data_full_nxs(tabf,n2) bind ( c )
            use iso_c_binding
            integer ( c_int ) :: tabf(*)
            integer ( c_int ), VALUE:: n2
        end subroutine c_get_data_bloc_data_full_nxs
    end interface

    private

    !---- Public Subroutines ----!
    public :: Read_Numor_D19_NXS2,fra1,&
        Read_Init_NXS, Read_Header_NXS, get_Header_Numor_NXS,&
        get_Header_ScanType_NXS,get_Header_Instr_Name_NXS,get_Header_SubT_NXS,&
        Read_IntegerBloc_NXS,get_IntegerBloc_NXS,&
        Read_FloatBloc_NXS,get_FloatBloc_NXS,&
        Read_DataBlocParam_NXS,get_DataBlocTime_NXS,get_DataBlocMoni_NXS,&
        get_DataBlocTotalCount_NXS,get_DataBlocAngle1_NXS,&
        get_DataBlocDataFull_NXS


  !---- Private Subroutines ----!
  !private::


!    interface
!       subroutine c_read_header_nxs() bind ( c )
!         use iso_c_binding
!       end subroutine c_read_header_nxs
!     end interface

Contains

    Subroutine fra1

    End Subroutine

    !!----
    !!---- Subroutine Read_Numor_D19(filename,n,frames)
    !!----    character(len=*)               , intent(in)    :: filename ! The input numor
    !!----    type(SXTAL_numor_type)         , intent(inout) :: n        ! The output numor structure
    !!----    integer, optional, dimension(:), intent(in)    :: frames   ! The frames to include in the numor structure
    !!----
    !!---- Subroutine to read a Numor of D19 Instrument at ILL
    !!----
    !!---- Counts: 640 x 256 = 163840
    !!----
    !!---- Update: 14/03/2011
    !!
    Subroutine Read_Numor_D19_NXS2

    End Subroutine Read_Numor_D19_NXS2

    Subroutine Read_Init_NXS(filename)
        use, intrinsic :: iso_c_binding
        implicit none
        character(len=*)               , intent(in)    :: filename
        !write (*,*) " ---------------   FRA21  " , filename
        call c_read_init_nxs(filename//".nxs"//C_NULL_CHAR)
    End Subroutine Read_Init_NXS


    Subroutine Read_Header_NXS
        use, intrinsic :: iso_c_binding
        implicit none
        call c_read_header_nxs()
    End Subroutine Read_Header_NXS

    subroutine get_Header_Numor_NXS (n)
        implicit none
        integer n
        call c_get_header_numor_nxs(n)
    end subroutine get_Header_Numor_NXS

    subroutine get_Header_ScanType_NXS (n)
        implicit none
        integer n
        call c_get_header_scantype_nxs(n)
    end subroutine get_Header_ScanType_NXS

    Subroutine get_Header_Instr_Name_NXS(fname,size1)
        use, intrinsic :: iso_c_binding
        implicit none
        character(len=*)               , intent(in)    :: fname
        integer size1
        call c_get_header_instr_name_nxs(fname,size1)
    End Subroutine get_Header_Instr_Name_NXS


    Subroutine get_Header_SubT_NXS(fname,size1)
        use, intrinsic :: iso_c_binding
        implicit none
        character(len=*)               , intent(in)    :: fname
        integer size1
        call c_get_header_subt_nxs(fname,size1)
    End Subroutine get_Header_SubT_NXS


    Subroutine Read_IntegerBloc_NXS
        use, intrinsic :: iso_c_binding
        implicit none
        call c_read_integer_bloc_nxs()
    End Subroutine Read_IntegerBloc_NXS

    Subroutine get_IntegerBloc_NXS(tabint ,size1 )
        use, intrinsic :: iso_c_binding
        implicit none
        integer , DIMENSION( : ) , intent(out)  :: tabint(:)
        integer size1,i
        call c_get_integer_bloc_nxs(tabint,size1)
    !      do i=1,size1
    !         write (*,*) ' cpp --> ',tabint(i)
    !      end do
    End Subroutine get_IntegerBloc_NXS


    Subroutine Read_FloatBloc_NXS
        use, intrinsic :: iso_c_binding
        implicit none
        call c_read_float_bloc_nxs()
    End Subroutine Read_FloatBloc_NXS

    Subroutine get_FloatBloc_NXS(tabf ,size1 )
        use, intrinsic :: iso_c_binding
        implicit none
        real , DIMENSION( : ) , intent(out)  :: tabf(:)
        integer size1,i
        call c_get_float_bloc_nxs(tabf,size1)
    !      do i=1,size1
    !         write (*,*) ' cpp --> ',tabf(i)
    !      end do
    End Subroutine get_FloatBloc_NXS


    Subroutine Read_DataBlocParam_NXS
        use, intrinsic :: iso_c_binding
        implicit none
        call c_read_data_bloc_param_nxs()
    End Subroutine Read_DataBlocParam_NXS

    Subroutine get_DataBlocTime_NXS(tabf ,size1 )
        use, intrinsic :: iso_c_binding
        implicit none
        real , DIMENSION( : ) , intent(out)  :: tabf(:)
        integer size1,i
        call c_get_data_bloc_time_nxs(tabf,size1)
    !      do i=1,size1
    !         write (*,*) ' cpp --> ',tabf(i)
    !      end do
    End Subroutine get_DataBlocTime_NXS


    Subroutine get_DataBlocMoni_NXS(tabf ,size1 )
        use, intrinsic :: iso_c_binding
        implicit none
        real , DIMENSION( : ) , intent(out)  :: tabf(:)
        integer size1,i
        call c_get_data_bloc_moni_nxs(tabf,size1)
    !      do i=1,size1
    !         write (*,*) ' cpp --> ',tabf(i)
    !      end do
    End Subroutine get_DataBlocMoni_NXS

    Subroutine get_DataBlocTotalCount_NXS(tabf ,size1 )
        use, intrinsic :: iso_c_binding
        implicit none
        real , DIMENSION( : ) , intent(out)  :: tabf(:)
        integer size1,i
        call c_get_data_bloc_total_count_nxs(tabf,size1)
    !      do i=1,size1
    !         write (*,*) ' cpp --> ',tabf(i)
    !      end do
    End Subroutine get_DataBlocTotalCount_NXS


    Subroutine get_DataBlocAngle1_NXS(tabf ,size1 )
        use, intrinsic :: iso_c_binding
        implicit none
        real , DIMENSION( : ) , intent(out)  :: tabf(:)
        integer size1,i
        call c_get_data_bloc_angle1_nxs(tabf,size1)
    !      do i=1,size1
    !         write (*,*) ' cpp --> ',tabf(i)
    !      end do
    End Subroutine get_DataBlocAngle1_NXS

    Subroutine get_DataBlocDataFull_NXS(tabf ,size1 )
        use, intrinsic :: iso_c_binding
        implicit none
        integer , DIMENSION(:) , intent(out)  :: tabf(:)
        integer size1,i
        call c_get_data_bloc_data_full_nxs(tabf,size1)
    !      do i=1,size1
    !         write (*,*) ' cpp --> ',tabf(i)
    !      end do
    End Subroutine get_DataBlocDataFull_NXS

#endif

End Module CFML_ILL_Instrm_Data_Nexus
