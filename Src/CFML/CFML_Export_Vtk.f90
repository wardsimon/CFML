!!----
!!---- Copyleft(C) 1999-2011,              Version: 5.0
!!---- Laurent Chapon
!!----
!!---- MODULE: CFML_Export_Vtk
!!----   INFO: Module to export volumetric data format to VTK format
!!----
!!---- HISTORY
!!----    Update: 23/12/2011
!!----
!!----
!!----
Module CFML_Export_VTK

use CFML_Globaldeps, Only: cp, dp
use CFML_Crystal_metrics, Only: Crystal_Cell_Type

implicit none

contains
!!----
!!---- Subroutine array3D_2_vts(Z,filename)
!!----    real(kind=dp), dimension(:,:,:), intent(in) :: Z			! 3D Array to export
!!----    Character(len=*),                intent(in) :: filename ! Name of the vtk file
!!----	Export a 3D array representing for example a density
!!----  into a vts file (vtk file format). This format correponds to the new xml format of vtk
!!----  as opposed to the legacy ascii format. vts corresponds to a structured grid
!!----  and the data is compressed using raw compression.
!!----  This routine should be used for pointData, i.e. data corresponding to a distribution:
!!----  Each point is defined by a
!!---- Update: December - 2011
!!
subroutine array3D_pointdata_2_vts(Z,cell,xmin,xmax,ymin,ymax,zmin,zmax,filename)
	real(kind=cp), dimension(:,:,:), intent(in) :: Z			! 3D Array to export
	type(Crystal_Cell_Type), intent(in)         :: cell   ! The unit-cell
	real(kind=cp), intent(in)                   :: xmin   ! Limits of x
	real(kind=cp), intent(in)                   :: xmax   ! Limits of x
	real(kind=cp), intent(in)                   :: ymin   ! Limits of y
	real(kind=cp), intent(in)                   :: ymax   ! Limits of y
	real(kind=cp), intent(in)                   :: zmin   ! Limits of z
	real(kind=cp), intent(in)                   :: zmax   ! Limits of z
	Character(len=*),                intent(in) :: filename ! Name of the vtk file
	! Local variables
	integer :: nx, ny, nz
	integer :: vtk_id, ier
	logical :: i_error
	integer :: i,j,k
	integer :: nbytes, nbits
	character(len=256) :: tpc, tamere
	integer :: offset
	real(kind=cp) ,dimension(3) :: pos

	! Check that the file is valid
	open(unit=vtk_id,file=trim(filename)//'.vts',form='unformatted',status='unknown',iostat=ier)
	if (ier/=0) then
		i_error=.true.
	end if
	! Check the extent of the array (C style 0:n-1)
	nx=size(Z,1)-1
	ny=size(Z,2)-1
	nz=size(Z,3)-1
	! Determine the number of bits used
	nbytes=sizeof(Z(1,1,1))
	nbits=8*nbytes
	write(tpc,*) nbits
	! Calculate the offset between points and Z
	offset=nbytes*(nx+1)*(ny+1)*(nz+1)+4
	write(tamere,*) offset
	! Write the header
  write(vtk_id) trim('<?xml version="1.0"?>')
	write(vtk_id) trim('<VTKFile byte_order="LittleEndian" version="0.1" type="StructuredGrid"> \n')
	write(vtk_id) '<StructuredGrid WholeExtent="', 0, nx, 0, ny, 0, nz, '">'
	write(vtk_id) '<Piece Extent="', 0, nx, 0, ny, 0, nz,'">'
	! Points represent the
	write(vtk_id) '<Points>'
	write(vtk_id) '<DataArray NumberOfComponents="3" offset="0" type="Float'//trim(tpc)//'" Name="points" format="appended"/>'
	write(vtk_id) '</Points>'
	write(vtk_id) '<PointData scalars="density">'
	write(vtk_id) '<DataArray NumberOfComponents="1" offset="'//trim(tamere)//'" type="Float'//trim(tpc)//'" Name="density" format="appended"/>'
	write(vtk_id) '</PointData>'
	write(vtk_id) '</Piece>'
	write(vtk_id) '</StructuredGrid>'
  !
	write(vtk_id) '<AppendedData encoding="raw">'
	! This is where the data is written
	do i=0,nx
		do j=0,ny
			do k=0,nz
				pos=(/xmin+(xmax-xmin)*i/nx,ymin+(ymax-ymin)*j/ny,zmin+(zmax-zmin)*k/nz/)
				pos=Matmul(cell%Cr_Orth_cel,pos)
				write(vtk_id) pos
				write(vtk_id) Z(i,j,k)
			end do
		end do
	end do
	! Finished to write the data, close the xml
	write(vtk_id) '</AppendedData>'
  write(vtk_id) '</VTKFile>'
  ! Close the file
  close(unit=vtk_id)

end subroutine array3D_pointdata_2_vts

end Module CFML_Export_Vtk


