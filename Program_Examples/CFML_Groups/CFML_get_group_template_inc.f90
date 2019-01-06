! Template algorithm for generating an arbitrary group characterized by matrices of
! whatever kind and dimensions
      max_op=size(Op)
      n=size(Op(1)%Mat,dim=1)
      call Allocate_Operator(n,Opt)
      done=.false.
      done(1,:) = .true.
      done(:,1) = .true.
      tb(1,:) = [(i,i=1,max_op)] !It is supposed that the first generator is the identity
      tb(:,1) = [(i,i=1,max_op)]
      nt=ngen
      Err_group=.false.
      Err_group_mess=" "
      !Ensure that determinants of generators are calculated
      do i=1,ngen
        Op(i)%dt=rdet(Op(i)%Mat(1:3,1:3))
      end do

      do_ext:do
        n=nt
        do i=1,n
          do_j:do j=1,n
            if(done(i,j)) cycle
            Opt=Op(i)*Op(j)
            do k=1,nt
              if(Opt == Op(k)) then
                tb(i,j)=k
                done(i,j)=.true.
                cycle do_j
              end if
            end do
            done(i,j)=.true.
            nt=nt+1
            if(nt > max_op) then
              nt=nt-1
              exit do_ext
            end if
            tb(i,j)=nt
            Op(nt)=Opt
          end do do_j
        end do
        if ( n == nt) exit do_ext
      end do do_ext

      if(any(done(1:nt,1:nt) .eqv. .false. ) ) then
      	Err_group=.true.
      	Err_group_mess="Table of SSG operators not exhausted! Increase the expected order of the group!"
      end if
      if(nt == max_op) then
        Err_group=.true.
      	write(Err_group_mess,"(a,i5,a)") "Max_Order (",max_op,") reached! The provided generators may not form a group!"
      end if

      multip=nt
      if(present(table)) then
        allocate(Table(multip,multip))
        Table=tb
      end if