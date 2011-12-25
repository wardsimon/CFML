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

use CFML_Globaldeps, Only: dp, dp
use CFML_Crystal_metrics, Only: Crystal_Cell_Type
use CFML_String_Utilities, Only: Pack_String

implicit none

character(1), parameter:: end_line=char(10)

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
	real(kind=dp), dimension(:,:,:), intent(in) :: Z      ! 3D Array to export
	type(Crystal_Cell_Type), intent(in)         :: cell   ! The unit-cell
	real(kind=dp), intent(in)                   :: xmin   ! Limits of x
	real(kind=dp), intent(in)                   :: xmax   ! Limits of x
	real(kind=dp), intent(in)                   :: ymin   ! Limits of y
	real(kind=dp), intent(in)                   :: ymax   ! Limits of y
	real(kind=dp), intent(in)                   :: zmin   ! Limits of z
	real(kind=dp), intent(in)                   :: zmax   ! Limits of z
	Character(len=*),                intent(in) :: filename ! Name of the vtk file
	! Local variables
	integer :: nx, ny, nz
	integer :: vtk_id, ier
	logical :: i_error
	integer :: i,j,k
	integer :: nbytes, nbits
	character(len=48)   :: nbit_c, offset_c, nx_c, ny_c, nz_c
	character(len=1000) :: header
	integer :: offset
	real(kind=dp), dimension(3) :: pos
    integer :: ntot ! Total number of points
    integer :: ntot_bytes


	! Open the file and check that the file is valid
	open(unit=vtk_id,file=trim(filename)//'.vts', &
         form       = 'UNFORMATTED',  &
         access     = 'SEQUENTIAL',   &
         action     = 'WRITE',        &
         convert    = 'LITTLE_ENDIAN',&
         recordtype = 'STREAM',       &
         buffered   = 'YES',          &
         iostat=ier)
	 if (ier/=0) then
		i_error=.true.
	 end if

	 ! Get the extent in each direction. Will be written as C array [0,n-1]
	 nx=size(Z,1)-1
     ny=size(Z,2)-1
     nz=size(Z,3)-1

     ! Total number of points
     ntot=(nx+1)*(ny+1)*(nz+1)

     ! Determine the number of bits used
     nbytes=sizeof(Z(1,1,1))
     nbits=8*nbytes
     !
     ntot_bytes=ntot*nbytes*3
     ! Calculate the offset between points and Z
     offset=nbytes*ntot*3+4
     !
     ! ----------- Here starts the header
	 write(unit=vtk_id) '<?xml version="1.0"?>'//end_line
     write(unit=vtk_id) '<VTKFile byte_order="LittleEndian" version="0.1" type="StructuredGrid">'//end_line
     write(unit=nx_c,fmt=*) nx
     write(unit=ny_c,fmt=*) ny
     write(unit=nz_c,fmt=*) nz
     write(unit=header,fmt='(A)') '<StructuredGrid WholeExtent="0 '//trim(Pack_String(nx_c))//' 0 '//trim(Pack_String(ny_c))//' 0 '//trim(Pack_String(nz_c))//'">'
     write(unit=vtk_id) trim(header)//end_line
     write(unit=header,fmt='(A)') '<Piece Extent="0 '//trim(Pack_String(nx_c))//' 0 '//trim(Pack_String(ny_c))//' 0 '//trim(Pack_String(nx_c))//'">'
     write(unit=vtk_id) trim(header)//end_line
     write(unit=vtk_id) '<Points>'//end_line
     write(unit=nbit_c,fmt=*) nbits
     write(unit=header,fmt='(A)') '<DataArray NumberOfComponents="3" offset="0" type="Float'//trim(Pack_String(nbit_c))//'" Name="points" format="appended"/>'
     write(unit=vtk_id) trim(header)//end_line
     write(unit=vtk_id) '</Points>'//end_line
     write(unit=vtk_id) '<PointData scalars="density">'//end_line
     write(unit=offset_c,fmt=*) offset
     write(unit=header,fmt='(A)') '<DataArray NumberOfComponents="1" offset="'//trim(Pack_string(offset_c))//'" type="Float'//trim(Pack_String(nbit_c))//'" Name="density" format="appended"/>'
     write(unit=vtk_id) trim(header)//end_line
     write(unit=vtk_id) '</PointData>'//end_line
     write(unit=vtk_id) '</Piece>'//end_line
     write(unit=vtk_id) '</StructuredGrid>'//end_line
     write(unit=vtk_id) '<AppendedData encoding="raw">'//end_line
     write(unit=vtk_id) '_'
    !------------- Here finishes the header
    ! This is where the data is written

    ! First we write the number of bytes for the points
     write(unit=vtk_id) ntot*nbytes*3
    ! Then write the coordinates of the points
	 do i=0,nx
        do j=0,ny
            do k=0,nz
                pos=(/xmin+(xmax-xmin)*i/nx,ymin+(ymax-ymin)*j/ny,zmin+(zmax-zmin)*k/nz/)
                pos=Matmul(cell%Cr_Orth_cel,pos)
                write(vtk_id) pos
            end do
          end do
        end do
     ! write the number of bytes for the point data
       write(unit=vtk_id) ntot*nbytes
     ! write the point data
      do i=0,nx
        do j=0,ny
            do k=0,nz
                write(vtk_id) Z(i+1,j+1,k+1)
            end do
        end do
     end do
  ! Finished to write the data, close the xml
  write(vtk_id) end_line
  write(vtk_id) '</AppendedData>'//end_line
  write(vtk_id) '</VTKFile>'//end_line
  ! Close the file
  close(unit=vtk_id)

end subroutine array3D_pointdata_2_vts


end Module CFML_Export_Vtk


