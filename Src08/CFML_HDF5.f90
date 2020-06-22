!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2019  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_HDF5
!!----         Module for working with HDF5 files.
!!----
!!
Module CFML_HDF5
    !---- Use Modules ----!
    use HDF5
    use CFML_ILL_Instrm_Data, ONLY: SXtal_Numor_type,Initialize_Numor,&
                                    err_ILLData,err_ILLData_mess
    use CFML_Strings,         ONLY: l_case

    implicit none

    private

    public: Read_Nexus_D19

    contains

    Interface
        Module Subroutine Read_Nexus_D19(filename,numor,counts)
            !---- Arguments ----!
            character(len=*),                       intent(in)  :: filename
            type(SXTAL_NUMOR_type),                 intent(out) :: numor
            integer, dimension(:,:,:), allocatable, intent(out) :: counts
        End Subroutine Read_Nexus_D19
    End Interface

End Module CFML_HDF5