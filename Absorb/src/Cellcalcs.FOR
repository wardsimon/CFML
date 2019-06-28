C -------------------------------------------------------------
C
C FUNCTION RETURNS UNIT-CELL VOLUME FROM CELL PARAMS
C
	FUNCTION GETCELLVOL(A,IRAD)
	DIMENSION A(6),C(3)

	RAD=ATAN(1.0)/45.0
C
C IRAD=1 FOR ANGLES IN RADIANS, 0 FOR ANGLES IN DEGREES
	FACT=1.0
	IF(IRAD .EQ. 0)FACT=RAD
	   
	DO I=1,3
		C(I)=COS(A(3+I)*FACT)
	ENDDO
C
C
	SUM= 1.0 +2.0*C(1)*C(2)*C(3)
	DO I=1,3
		SUM=SUM-C(I)*C(I)
	ENDDO
	GETCELLVOL=A(1)*A(2)*A(3)*SQRT(SUM)
	RETURN
	END
C ****************************************************************
C RECIP CELL FROM UB
C
	SUBROUTINE RECIP_CELL
	INCLUDE 'unitcell.inc'
	CHARACTER*256 ETEXT
	rad=45.0/atan(1.0)
c
c get recip metric tensor
	ubtranspose=transpose(ub)
	gstar=matmul(ubtranspose,ub)
c
c check for symmetric
	if(abs(gstar(1,2)-gstar(2,1)) .gt. 0.0001)goto99
	if(abs(gstar(1,3)-gstar(3,1)) .gt. 0.0001)goto99
	if(abs(gstar(3,2)-gstar(2,3)) .gt. 0.0001)goto99
c
c get recip lattice axes
	do i=1,3
		if(gstar(i,i) .le. 0)goto99
		recell(i)=sqrt(gstar(i,i))
	enddo
c
c get angles in degrees
	recell(4)=acos(gstar(2,3)/recell(2)/recell(3))*rad
	recell(5)=acos(gstar(1,3)/recell(1)/recell(3))*rad
	recell(6)=acos(gstar(1,2)/recell(1)/recell(2))*rad
c
c get (trans(UB))-1 for shelx dcosines
	ubti=transpose(ub)
	CALL MINV3(UBTI,UBTI,1.0E-7,DET,IERR)
	IF(IERR .EQ. 0)RETURN
C
C 
	WRITE(ETEXT,11)UBTI
11	FORMAT('INVALID UB: (UBTRANSPOSE)-1= ',3(10x,3f10.6,':'))
	CALL ERROR(ETEXT)
	RETURN
c  errors here 
c
99	write(ETEXT,10)gstar
10	format('  INVALID RECIP METRIC TENSOR: ',3(10x,3f10.6,':'))
	CALL ERROR(ETEXT)
	RETURN
	end
c
C ****************************************************************
C REAL CELL FROM RECIP METRIC TENSOR GSTAR
C
	SUBROUTINE REAL_CELL
	INCLUDE 'unitcell.inc'
		include 'files.inc'
	REAL G(3,3),GT(3,3)
	CHARACTER*256 ETEXT
	rad=45.0/atan(1.0)
c
c to avoid problems with big cells do test and rescaling before inversion of G*
      small = 1.0
      do i=1,3
          if(gstar(i,i) .lt. tiny(0.0))goto 95        ! diagonal element zero
          if(gstar(i,i) .lt. small)small=gstar(i,i)
      enddo
      gt=gstar/small      ! rescaled version of gstar
c
c get real metric tensor
	CALL MINV3(GT,GT,1.0E-9,DET,IERR)   
	if(IERR .NE. 0)then
		write(ETEXT,98)det
98		format('   Determinant of Scaled Recip Metric tensor = ',f10.4)
c write gstar to log file
c
		goto 90
      else
          g=gt/small                      ! get metric tensor scaled
		aub(7)=1.0/sqrt(det)/small**1.5			! VOLUME
	endif
c
c check for symmetric
	if(abs(g(1,2)-g(2,1)) .gt. 0.0001)goto99
	if(abs(g(1,3)-g(3,1)) .gt. 0.0001)goto99
	if(abs(g(3,2)-g(2,3)) .gt. 0.0001)goto99
c
C get real lattice axes
	do i=1,3
		if(g(i,i) .le. 0)goto99
		aub(i)=sqrt(g(i,i))
	enddo
c
c get angles in degrees
	AUB(4)=acos(g(2,3)/aub(2)/aub(3))*rad
	AUB(5)=acos(g(1,3)/aub(1)/aub(3))*rad
	AUB(6)=acos(g(1,2)/aub(1)/aub(2))*rad
c
	return
c  errors here 
c
95    write(ETEXT,12)
12    format('INVALID RECIPROCAL METRIC TENSOR: ')
      GOTO 90
      
99	write(ETEXT,10)g
10    format('INVALID METRIC TENSOR: ',3(10x,3f10.6,';'))
c
c write out gstar if problem BEFORE CALLING ERROR TO CLOSE PROGRAM!
90	write(iprt,11)
11	format(/' The Recip metric tensor is:'/)
	do i=1,3
		write(iprt,20)(gstar(i,j),j=1,3)
20		format(10x,3f10.6)
	enddo
	CALL ERROR(ETEXT)
	RETURN
	end
c
c*********************************************************************
C
	SUBROUTINE ROTATE_UB
C
C SUBROUTINE ROTATES UB ON INPUT TO B-L AXES ON BASIS OF BLAXES CARD
C
	include 'unitcell.inc'
	INTEGER IT(3,3)
	CHARACTER*256 ETEXT
C
C NEED IBL TO TELL US ORIENTATION OF UB W.R.T. BUSING-LEVY AXES
C
C AFTER CALLING THIS SUBROUTINE UB IS THE ORIENTATION MATRIX ON BL AXES
C SO INPUT UB IS LOST.
C
	IBLTEST=0
	DO I=1,3
		IBLTEST=IBLTEST+IABS(IBL(I))
	ENDDO
C
	IF(IBLTEST .NE. 6)THEN
		WRITE(ETEXT,20)IBLTEST
20		FORMAT('ERROR IN EXP FILE: '
	1			' BLAXES CARD ENTRIES ADD TO ',I1,
     2	5X,'AXIAL SYSTEM ASSUMED BUSING-LEVY')
		DO I=1,3
			IBL(I)=I	! BUSING-LEVY AXES
		ENDDO
		CALL WARNING(ETEXT)
		RETURN			! NO ROTATION REQUIRED		
	ENDIF
C
C ALL PRESENT AND CORRECT: ROTATE UB
C
	DO I=1,3
		DO J=1,3
			IT(I,J)=0.
		ENDDO
	ENDDO
	DO I=1,3
		IT(IABS(IBL(I)),I)=SIGN(1,IBL(I))
	ENDDO
C
	UB=MATMUL(IT,UB)
C
	RETURN
	END
C
C
	SUBROUTINE DC_TOSHELX(DC)
C
C SUBROUTINE USES THE UB MATRIX AND RECIPROCAL LATTICE PARAMETERS
C TO CONVERT BL DIRECTION COSINES INTO SHELX COSINES W.R.T. UNIT-CELL AXES
C 
C NOTE: THE AXIAL ORIENTATION OF THE UB IS B-L WHICH IS THE SAME AS
C       THE INTERNAL DCOSINES ON THE PHI-AXIS SYSTEM SO NO WORRIES ABOUT
C       THE ORIENTATION OF COORDINATE SYSTEM
C
	INCLUDE 'unitcell.inc'
	include 'reflection.inc'
	REAL DC(6),DCI(3),DCD(3)
C
	IF(IUB .NE. 3)THEN
		DO I=3,6
			DC(I)=0.
		ENDDO
		DC(1)=1.
		DC(2)=1.
		RETURN
	ENDIF
			
	DCI=MATMUL(UBTRANSPOSE,COSP)
	DCD=MATMUL(UBTRANSPOSE,COSD)
C
	DO I=1,3
		ID=2*I
		IP=ID-1
		DC(IP)=DCI(I)/RECELL(I)
		DC(ID)=DCD(I)/RECELL(I)
	ENDDO
	RETURN
	END
