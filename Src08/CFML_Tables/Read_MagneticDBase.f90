!!----
!!----
!!----
!!----
SubModule (CFML_Magnetic_Database) MagDB_002
   Contains
   !!----
   !!---- READ_MAGNETIC_DATA
   !!----
   !!---- Read data about magnetic space groups
   !!---- input data from magnetic_table.dat
   !!----
   !!---- 24/04/2019 
   !!
   Module Subroutine Read_Magnetic_Data()
      !---- Local Variables ----!
      integer            :: i,j,k,n,m,i_mag,ier
      character(len=512) :: fullprof_suite, database
    
      !> Init
      call clear_error()
       
      !> open data file
      call GET_ENVIRONMENT_VARIABLE("FULLPROF",fullprof_suite)
      n=len_trim(fullprof_suite)
      if (n == 0) then
         err_CFML%IErr=1 
         write(unit=err_cfml%msg,fmt="(a)") " => The FULLPROF environment variable is not defined! "//newline// &
                                            "    This is needed for localizing the data base: magnetic_data.txt"//newline// &
                                            "    that should be within the %FULLPROF%/Databases directory"
         return
      end if   

      if (fullprof_suite(n:n) /= OPS_SEP) then
         database=trim(fullprof_suite)//OPS_SEP//"Databases"//OPS_SEP//'magnetic_data.txt'
      else
         database=trim(fullprof_suite)//"Databases"//OPS_SEP//'magnetic_data.txt'
      end if

      !> Open 
      open(newunit=i_mag,File=trim(database),status="old",action="read",position="rewind",iostat=ier)
      if ( ier /= 0) then
         err_CFML%IErr=1  
         err_CFML%Msg="Read_Magnetic_Data@SPACEG: Problem opening the data base: "//trim(database)
         return
      end if
    
      if (.not. Magnetic_DBase_allocated) call allocate_Magnetic_DBase()
      
      !> read nonhexagonal point operators
      do i=1,48
         Read(i_mag,*) n,point_op_label(i),point_op_xyz(i),  &
                       ((point_op_matrix(k,j,i),j=1,3),k=1,3)
         if (n /= i) then
            err_CFML%IErr=1 
            err_CFML%Msg="Read_Magnetic_Data@SPACEG: Error in numbering of nonhexagonal point operators"
            return
         end if
      end do
      
      !> read hexagonal point operators
      do i=1, 24
         Read(i_mag,*) n,point_op_hex_label(i), point_op_hex_xyz(i),  &
                       ((point_op_hex_matrix(k,j,i),j=1,3),k=1,3)
         if (n /= i) then
            err_CFML%IErr=1 
            err_CFML%Msg="Read_Mag_Data@SPACEG: Error in numbering of hexagonal point operators"
            return
         end if
      end do
    
      !> read data for each magnetic space group
      do i=1,1651
         Read(i_mag,*) (nlabelparts_bns(j,i),j=1,2),nlabel_bns(i),  &
                       spacegroup_label_bns(i),(nlabelparts_og(j,i),j=1,3),  &
                       nlabel_og(i),spacegroup_label_og(i)
      
         Read(i_mag,*) magtype(i)
         If (magtype(i) == 4) Then
            Read(i_mag,*) ((bnsog_point_op(j,k,i),j=1,3),k=1,3),  &
                          (bnsog_origin(j,i),j=1,3),bnsog_origin_denom(i)
         End If
      
         Read(i_mag,*) ops_count(i)
         Read(i_mag,*) (ops_bns_point_op(j,i),(ops_bns_trans(k,j,i),k=1,3),  &
                       ops_bns_trans_denom(j,i),ops_bns_timeinv(j,i), j=1,ops_count(i))
         Read(i_mag,*) lattice_bns_vectors_count(i)
         Read(i_mag,*) ((lattice_bns_vectors(k,j,i),k=1,3),  &
                       lattice_bns_vectors_denom(j,i), j=1,lattice_bns_vectors_count(i))
         Read(i_mag,*) wyckoff_site_count(i)
         Do j=1,wyckoff_site_count(i)
            Read(i_mag,*) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
            Do k=1,wyckoff_pos_count(j,i)
               Read(i_mag,*) (wyckoff_bns_fract(m,k,j,i),m=1,3),  &
                             wyckoff_bns_fract_denom(k,j,i),  &
                             ((wyckoff_bns_xyz(m,n,k,j,i),m=1,3),n=1,3),  &
                             ((wyckoff_bns_mag(m,n,k,j,i),m=1,3),n=1,3)
            End Do
         End Do
      
         If (magtype(i) == 4) Then
            Read(i_mag,*) ops_count(i)
            Read(i_mag,*) (ops_og_point_op(j,i),(ops_og_trans(k,j,i),k=1,3),  &
                          ops_og_trans_denom(j,i),ops_og_timeinv(j,i), j=1,ops_count(i))
            Read(i_mag,*) lattice_og_vectors_count(i)
            Read(i_mag,*) ((lattice_og_vectors(k,j,i),k=1,3),  &
                          lattice_og_vectors_denom(j,i), j=1,lattice_og_vectors_count(i))
            Read(i_mag,*) wyckoff_site_count(i)
            Do j=1,wyckoff_site_count(i)
               Read(i_mag,*) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
               Do k=1,wyckoff_pos_count(j,i)
                  Read(i_mag,*) (wyckoff_og_fract(m,k,j,i),m=1,3),  &
                                wyckoff_og_fract_denom(k,j,i),              &
                                ((wyckoff_og_xyz(m,n,k,j,i),m=1,3),n=1,3),    &
                                ((wyckoff_og_mag(m,n,k,j,i),m=1,3),n=1,3)
               End Do
            End Do
         End If
      End Do
    
      !> close data file
      Close(i_mag)
   End Subroutine Read_Magnetic_Data
   
   !!----
   !!---- READ_MAGNETIC_BINARY
   !!----
   !!----     read data about magnetic space groups
   !!----      input data from magnetic_table.dat
   !!----
   !!---- 24/04/2019 
   !!
   Module Subroutine Read_Magnetic_Binary()
      !---- Local Variables ----!
      integer :: i,j,k,n,m,i_mag

      !> open data file
      Open(newunit=i_mag,File='magnetic_data.bin',status="old",action="read",form="unformatted",access="stream")
    
      !>For the old Lahey compiler use this
      !Open(unit=i_mag,File='magnetic_data.bin',status="old",action="read",form="unformatted",access="transparent") ! For Lahey
    
      !> read nonhexagonal point operators
      do i=1,48
         Read(i_mag)n,point_op_label(i),point_op_xyz(i),  &
                    ((point_op_matrix(k,j,i),j=1,3),k=1,3)
      end do
    
      !> read hexagonal point operators
      do i=1,24
         Read(i_mag) n,point_op_hex_label(i), point_op_hex_xyz(i),  &
                     ((point_op_hex_matrix(k,j,i),j=1,3),k=1,3)
      end do
    
      !> read data for each magnetic space group
      do i=1,1651
         Read(i_mag) (nlabelparts_bns(j,i),j=1,2),nlabel_bns(i),  &
                     spacegroup_label_bns(i),(nlabelparts_og(j,i),j=1,3),  &
                     nlabel_og(i),spacegroup_label_og(i)
         Read(i_mag) magtype(i)
         If (magtype(i) == 4) Then
            Read(i_mag) ((bnsog_point_op(j,k,i),j=1,3),k=1,3),  &
                        (bnsog_origin(j,i),j=1,3),bnsog_origin_denom(i)
         End If
         Read(i_mag) ops_count(i)
         Read(i_mag) (ops_bns_point_op(j,i),(ops_bns_trans(k,j,i),k=1,3),  &
                     ops_bns_trans_denom(j,i),ops_bns_timeinv(j,i), j=1,ops_count(i))
         Read(i_mag) lattice_bns_vectors_count(i)
         Read(i_mag) ((lattice_bns_vectors(k,j,i),k=1,3),  &
                     lattice_bns_vectors_denom(j,i), j=1,lattice_bns_vectors_count(i))
         Read(i_mag) wyckoff_site_count(i)
         Do j=1,wyckoff_site_count(i)
            Read(i_mag) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
            Do k=1,wyckoff_pos_count(j,i)
               Read(i_mag) (wyckoff_bns_fract(m,k,j,i),m=1,3),  &
                           wyckoff_bns_fract_denom(k,j,i),  &
                           ((wyckoff_bns_xyz(m,n,k,j,i),m=1,3),n=1,3),  &
                           ((wyckoff_bns_mag(m,n,k,j,i),m=1,3),n=1,3)
            End Do
         End Do
         If (magtype(i) == 4) Then
            Read(i_mag) ops_count(i)
            Read(i_mag) (ops_og_point_op(j,i),(ops_og_trans(k,j,i),k=1,3),  &
                        ops_og_trans_denom(j,i),ops_og_timeinv(j,i), j=1,ops_count(i))
            Read(i_mag) lattice_og_vectors_count(i)
            Read(i_mag) ((lattice_og_vectors(k,j,i),k=1,3),  &
                        lattice_og_vectors_denom(j,i), j=1,lattice_og_vectors_count(i))
            Read(i_mag) wyckoff_site_count(i)
            Do j=1,wyckoff_site_count(i)
               Read(i_mag) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
               Do k=1,wyckoff_pos_count(j,i)
                  Read(i_mag) (wyckoff_og_fract(m,k,j,i),m=1,3),  &
                              wyckoff_og_fract_denom(k,j,i),              &
                              ((wyckoff_og_xyz(m,n,k,j,i),m=1,3),n=1,3),    &
                              ((wyckoff_og_mag(m,n,k,j,i),m=1,3),n=1,3)
               End Do
            End Do
         End If
      End Do
    
      !> close data file
      Close(i_mag)
   End Subroutine Read_Magnetic_Binary
   
   
End SubModule MagDB_002  