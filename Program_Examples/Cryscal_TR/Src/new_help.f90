program replace
 implicit none
 character(len=256) :: read_line
 integer            :: i1, i2, ier, n, case_numor

 open(unit=1, file='case.txt')
 open(unit=2, file='help_line.txt')

 write(*,*) '  > Enter CASE numor: '
 read(*,*) case_numor
 n=0
 do
  read(unit=1, '(a)', iostat=ier) read_line
  if(ier/=0) exit
  read_line = adjustl(read_line)
  i1 = index(read_line,'call write_info(')
  if (i1/=0) then
   i1 = index(read_line, '(')
   i2 = index(read_line, ')')

   n=n+1
   if(case_numor < 10) then
    write(2, '(a,i1,a,i2,2a)') '  HELP_line(',case_numor,',',n, ") = ",read_line(i1+1:i2-1)
   elseif(case_numor < 100) then
    write(2, '(a,i2,a,i2,2a)') '  HELP_line(',case_numor,',',n, ") = ",read_line(i1+1:i2-1)
   endif
  endif

 end do

 if(case_numor < 10) then
  write(2, '(a,i1,a,i2)') '  HELP_lines_nb(',case_numor,') = ', n
 elseif(case_numor < 100) then
  write(2, '(a,i2,a,i2)') '  HELP_lines_nb(',case_numor,') = ', n
 endif

end program replace
