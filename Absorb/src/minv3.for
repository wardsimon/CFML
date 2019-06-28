C -------------------------------------------------------------
C
C SUBROUTINE INVERT3 FROM LWF
C GENERALISED FEB 2003 INTO A REAL*4 AND REAL*8 ROUTINES MINV3 AND MINV38
C 
C NEEDED BECAUSE LU DECOMPOSITION DOES NOT HANDLE DET=0 VERY WELL.
C
C KEEP THE INVERT3 SUBROUTINE FOR REFERENCE!!
C
      subroutine INVERT3(a,B,d)
c     INVERTS a 3x3 matrix IN A AND PUTS (A-1)T IN B
C     on return
c     d is the value of the determinant of the matrix. if d is less than
c     1.0e-6 no replacement takes place.
      real*8 a(9),p(9),B(3,3),D
      d = 0.
      ii = 0
      io = 0
      jo = 3
      ko = 6
      do 1 ic = 1,3
      j = 2
      k = 3
      do 2 i = 1,3
      ii = ii+1
      p(ii) = a(jo+j)*a(ko+k) - a(jo+k)*a(ko+j)
      j = k
    2 k = i
      d = d + p(ii)*a(ii)
      jo = ko
      ko = io
    1 io = io + 3
      if (Dabs(d) .lt. 1.0e-7) go to 4
      do  i = 1,3
	do j=1,3
	  k= 3*(i-1)+j
	  B(i,j) = p(k)/d
	enddo
      enddo
    4 return
      end
C
C
	subroutine MINV3(a,B,DMIN,d,IERR)
C
C SUBROUTINE TO INVERT A 3X3 MATRIX IN REAL*4
C using the minv38 routine
C
C IF DETERMINANT IS LESS THAN DMIN, IERR=1 AND NO INVERSION IS PERFORMED
C
	REAL A(3,3), B(3,3), DMIN,D
	REAL,PARAMETER::DOUBLE=SELECTED_REAL_KIND(p=13,r=200)

	REAL(KIND=DOUBLE):: A8(3,3), B8(3,3), DMIN8,D8
C
C COPY A TO A8 
	DO I=1,3
		DO J=1,3
			A8(I,J)=A(I,J)
		ENDDO
	ENDDO
	DMIN8=DMIN
C
	CALL MINV38(A8,B8,DMIN8,D8,IERR)
C
	D=D8
	DO I=1,3
		DO J=1,3
			B(I,J)=B8(I,J)
		ENDDO
	ENDDO
	RETURN
	END	

      subroutine MINV38(a,B,DMIN,d,IERR)
c     INVERTS a 3x3 matrix IN A AND PUTS (A-1) IN B
C     on return
c     d is the value of the determinant of the matrix. if d is less than
c     DMIN no replacement takes place.
      REAL,PARAMETER::DOUBLE=SELECTED_REAL_KIND(p=13,r=200)
	REAL(KIND=DOUBLE):: a(9),p(9),B(3,3),D,DMIN
	INTEGER IERR
	IERR=0
      d = 0.
C
C CLEAR RESULT MATRIX 
		DO I=1,3
			DO J=1,3
				B(I,J)=0.0D0
			ENDDO
		ENDDO
      ii = 0
      io = 0
      jo = 3
      ko = 6
      do 1 ic = 1,3
      j = 2
      k = 3
      do 2 i = 1,3
      ii = ii+1
      p(ii) = a(jo+j)*a(ko+k) - a(jo+k)*a(ko+j)
      j = k
    2 k = i
      d = d + p(ii)*a(ii)
      jo = ko
      ko = io
    1 io = io + 3
      if (Dabs(d) .lt. DMIN)THEN
		D=0.0
		IERR=1
		RETURN
	ENDIF
      do  i = 1,3
	do j=1,3
	  k= 3*(i-1)+j
	  B(i,j) = p(k)/d
	enddo
      enddo
      return
      end