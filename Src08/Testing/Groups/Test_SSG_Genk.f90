    Program Test_ssg_Genk
      Use CFML_GlobalDeps
      Use CFML_gSpaceGroups
      Use CFML_gSpaceGroups
      Use CFML_Propagation_Vectors, only: K_Star, Write_Group_K, Set_Gk, Group_k_Type
      Use CFML_SuperSpace_Database
      Use CFML_IOForm, only : Read_CFL_SpG,Read_CFL_Cell,Read_kinfo
      use CFML_Rational

      implicit none

      integer                         :: nkv,i,j,k,m,multip
      character(len=50)               :: str,forma
      real(kind=cp), dimension(3,4)   :: kv
      type(Group_k_Type),dimension(4) :: Gk
      type(SpG_Type)                  :: SpG
      type(SpG_Type),dimension(4)     :: Grpk
      logical                         :: ext=.true.
      !call Read_SSG_DBase()
      !if(Err_CFML%Ierr /= 0) then
      !  write(*,"(a)") "   !!! "//trim(Err_CFML%Msg)//" !!!"
      !  stop
      !end if

	    do
        write(*,"(a)",advance="no") " => Enter the number (or the symbol) of a space group: "
        read(*,"(a)") str
	    	if(len_trim(str) == 0) exit
        write(*,"(a)",advance="no") " => Enter the number of propagation vectors: "
        read(*,*) nkv
        do i=1,nkv
           write(*,"(a)",advance="no") " => Enter a propagation vector: "
           read(*,*) kv(:,i)
        end do

		    call Set_SpaceGroup(str,SpG)
		    call Write_SpaceGroup_Info(SpG)

		    do i=1,nkv
		      call K_Star(kv(:,i),SpG,Gk(i),ext)
		      call Set_Gk(Gk(i),Grpk(i),ext)
		      write(*,"(//a,3f10.5/)") " => GROUP OF THE PROPAGATION VECTOR: ",kv(:,i)
		      call Write_Group_K(Gk(i))
		      write(*,"(//a/)") " => FULL PROPAGATION VECTOR SPACE GROUP INCLUDING -k (Extended Little Group): "
		      call Write_SpaceGroup_Info(Grpk(i))
		    end do
		    !Determine now the possible superspace groups supposing that the propagation vectors are purely magnetic
		    !That is magnetic superspace groups of type 4. Just adding the operation 1'(0 0 0 1/2) and adding in the
		    !addiditonal dimensions a phase of different types depending on the order of operators
        !Take up to three generators of the Extended Little Group
        !Determine the order
      end do

    contains
      Subroutine Get_Symbol_SSG_from_Operators()
        !Translations along extra coordinates allowed for the different order of operators
        !  t =  0  1/2  1/3  -1/3   1/4  -1/4   1/6   -1/6
        !Symb   0   s    t   -t      q    -q     h     -h
      End Subroutine Get_Symbol_SSG_from_Operators
    End Program Test_ssg_Genk