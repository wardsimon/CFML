    Program read_ssg_database
      Use CFML_GlobalDeps
      Use CFML_gSpaceGroups
      Use CFML_SuperSpace_Database
      Use CFML_IOForm, only : Read_CFL_SpG,Read_CFL_Cell,Read_kinfo
      use CFML_Rational

      implicit none

      integer :: iclass,nmod,i,j,k,m,multip 
      character(len=50)  :: str,forma
      !type(SuperSpaceGroup_Type) :: SSpaceGroup

      !call Read_SSG_DBase()
      !if(Err_CFML%Ierr /= 0) then
      !  write(*,"(a)") "   !!! "//trim(Err_CFML%Msg)//" !!!"
      !  stop
      !end if
      
	  do
        write(*,"(a)",advance="no") " => Enter the number (or the symbol) of the SSG: "
        read(*,"(a)") str
        !if(m <= 0) exit
		if(len_trim(str) == 0) exit
        !if(m > 16697) then
        !  write(*,"(a)") " => There are only 16697 superspace groups in the database! "
        !  cycle
        !end if
		
		call Read_single_SSG(str,m)
		if(err_CFML%Ierr /= 0) then 
		  write(*,"(a)") trim(err_CFML%Msg)
		  cycle
		end if
        write(*,"(4(a,i5))") " Order number:",m, " Group number:",igroup_number(m), " Bravais class:",igroup_class(m), "  Basic Group #:",igroup_spacegroup(m)
        iclass=igroup_class(m)
        nmod=iclass_nmod(iclass)
        write(*,"(a,tr4,a)") "Numeric Label: "//group_nlabel(m), "    Group Label: "//group_label(m)
        write(*,"(a,i3)") " Number of operators:", igroup_nops(m)
        do k=1,igroup_nops(m)
          write(*,"(a,i3)") " Operator #",k
          forma="(  i4)"
          write(forma(2:3),"(i2)") nmod+4
          do i=1,nmod+4
            write(*,forma) igroup_ops(i,1:nmod+4,k,m)
            !write(*,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))
          end do
        end do

      end do

    End Program read_ssg_database