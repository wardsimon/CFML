C -------------------------------------------------------------
C ---- HighPressure software
C -------------------------------------------------------------
C                          ABSORB 7.0
C ABSORPTION CORRECTION PROGRAM FOR CRYSTALS IN DAC'S
C INCLUDING DAC ABSORPTION, GASKET SHADOWING WITH PARTIAL ABSN
C AND CRYSTAL ABSORBTION CORRECTIONS. 
C
C WRITTEN APRIL 1992 R.J. ANGEL, UCL MINERAL PHYSICS.
C MODIFIED SEPT 1993 R.J. ANGEL - INCLUSION OF GASKET ABSORBTION
C
C ***** TRNS AND DTRNS CALCULATIONS CHECKED OK AGAINST CYLINDER DATA
C       IN ITABLES VOL. C, by RJA SEPT 94
C
C ***** STTH IN SUBROUTINE ITOF HAD ABS FUNCTION ADDED - RJA 31-JAN-95
C
C ***** MORE USER INTERACTION ADDED FOR DAC DIMENSIONS - RJA MAY 1996
C
C ***** PORTED TO LAHEY PC FORTRAN JULY 1998 RJA
C
C ***** MODIFIED TO WRITE BETA=0.0 ON OUTPUT IF NO ABSORPTION CALCULATED
C        RJA 11 NOV 99
C ***** new absorption coeffs inserted for Be and diamond 
C        RJA 9-Dec-1999
C
C Correction TO GENERAL EXPRESSION FOR MONOCHROMATOR POLARISATION  20-MAR-2000
c
C
C *****VERSION 5.0 STARTED SEPT 2001 RJA
C		NEW EXPANDED INPUT FILE FORMAT
C		SWITCH TO BUSING-LEVY COORDINATE SYSTEM 
C         (ONLY AFFECTS DIRN COSINES AND DAC CORRECTION ALGEBRA)
C		SPLIT UP CODE INTO SEVERAL FILES
C
C
C *****VERSION5.1 STARTED APRIL 2002 RJA
C          RESTRUCTURED READING OF CONTROL INFO INTO A CRYSTAL FILE AND
C          MAKE USE OF INSTURMENT PARAMETER FILE
C
C *****VERSION 5.2 AUGUST 2002 RJA
C		RESTRUCTURING OF PROGRAM
C
C *****VERSION 5.3 OCTOBER 2002 RJA
C		CLEAN UP OF BUGS, ITOF RESTRUCTURED, ADDED STUFF
C         JAN 2003 -
C            REORGANISED COMMON BLOCKS
C            SWITCHED TO FREE FORMAT READ OF EXP FILE
C         APRIL 2003 REAL*8 CHANGED TO KIND(8) 
C                    COMPILED WITH /FPE=0
C
C *****VERSION 5.4  JUNE/JULY 2003 -
C              ADDED TWO ANVILS
C			 RECONFIGURED PSII, PSID CALCS
C			 REORGANISED I/O 
C
C *****VERSION 6 BUILT DIRECT FROM V5.4 FOR DIALOG/WINDOWS
C
C *****version 6.1: BUG FIXES TO 6 AND TRAPS FOR ILLEGAL COMBINATIONS OF INPUT
C
C ******VERSION 7  NEW STAND-ALONE VERSION TO BETTER INTERFACE WITH COMMERCIAL
C       SOFTWARE. THIS RUNS FROM THE CONFIG FILE
C          
C *****************************************************************
C BASED IN PART UPON ABSORB (BURNHAM, 1966; MODIFIED BY LW FINGER)*
C *****************************************************************

C
C  THIS is the main program for the ABSORB dialog
C
C

	INCLUDE 'FLAGS.INC'			! TO INITIALISE FLAGS
	INCLUDE 'FILES.INC'			! TO INITIALISE FORMAT FLAGS
	INCLUDE 'UNITCELL.INC'		! NEEDED TO CHECK FOR UB
	CHARACTER STARS*25
	STARS='*************************'
C
C Version information
C                   123456789012      
      version_date='15-JUL-2017'      !not to exceed 12 characters
      version_no=7.21                 !f5.2, so only 2 dp
C
C INITIALISE FLAGS
	IFATAL=0			! FATAL ERROR FLAG
	IBUILD=0			! CRYSTAL MODEL BUILD FLAG
	IDFORM_IN=0			! INPUT FORMAT - NOT SET
	IDFORM_OUT=0		! OUTPUT FORMAT - NOT SET
	IABS_RUN=1			! DEFAULT IS SHOW MESSAGES
C
C      
C SET THE LOCK for safety first
	CALL LOCK(1,ierr)
C DELETE THE ERROR LOG: 
	CALL DELETE_ERROR_LOG
      if(ierr .ne. 0)then
          WRITE(ESTRING,'(''Absorb lock file exists:'',
     1          '' another version of the program may be running'')')
          call ERROR(estring)   ! cannot use error because the error log may be in use by another program
          goto 9998       ! do not use close program because that uses lock file
      endif    

 
C
C SET THE CONFIGURATION
C
	CALL SET_CONFIG
C
C OPEN THE CONFIG FILE AND GET THE FILE NAMES: 
	CALL READ_CONFIG_FILE(IERR)
C
	IF(IERR .NE. 0)THEN
C WARNING IN CONFIG FILE ABOUT FORMAT NUMBERS
		IF(IERR .EQ. 3 .or. IERR .EQ. 4)THEN
			CALL FIX_FILEFORMAT(IERR)			! WILL CLEAR IERR IF SUCCESSFUL
			IF(IERR .EQ. 0)GOTO 10
C
C OFOLLOWING ARE  ERRORS READING CONFIG FILE - 
		ELSEIF(IERR .EQ. -1)THEN
			WRITE(ESTRING,'(''Cannot open Absorb_config.dat file'')')
		ELSEIF(IERR .EQ. 1)THEN
			WRITE(ESTRING,'(''Absorb_config.dat file not found'')')
		ELSEIF(IERR .EQ. 2)THEN
			WRITE(ESTRING,'(''Error reading Absorb_config.dat file'')')
		ENDIF
		CALL ERROR(ESTRING)
		GOTO 9999
	ENDIF
C	IF(IABS_RUN .NE. 0)CALL CREATEMESSAGEBOX
C 
C NOW OPEN THE FILES
C
10	CALL OPEN_ABSORB_FILES(IERR)
	IF(IERR .NE. 0)then
          estring=''
          select case(ierr)
          case(1)
              estring=trim('Error opening exp file')
          case(2)
              estring=trim('Error opening print file')
          case(3)
              estring=trim('Error opening input data file')
          case(4)
              estring=trim('Error opening output data file')
          case default
              estring=trim('Unknown error while opening files')
          end select
          call error(estring)
          goto 9999
      endif
      

C WRITE HEADER TO PRINT FILE
C
	CALL WRITE_HEADER(iprt)
C
C   READ DATA FROM EXP FILE
C
	CALL READ_CRYSTAL_DATA
	IF(IFATAL .EQ. 1)GOTO 9999		! ERROR BOXED POPPED BY READ ROUTINES
C
C DO ANY NECESSARY PRE-PROCESSING
	CALL PRE_PROCESS
C
C
	IF(IGRID .NE. 0)CALL OPEN_GRID
	
C
C SET UP CRYSTAL MODEL
C
30	CALL BUILD_MODEL
	IF(IFATAL .EQ. 1)GOTO 9999
	IF(IBUILD .NE. 1)GOTO 9999
C
C WRITE DATA FILE NAMES TO PRINT FILE
C
	WRITE(IPRT,23)STARS,STARS
23	FORMAT(/A25,'****** DATA FILES*********',A25/)
C
	WRITE(IPRT,24)FNAME(3)(1:Ilen(3))
24	FORMAT(/5X,'INPUT DATA FILE: ',A)
	WRITE(IPRT,28)FNAME(4)(1:Ilen(4))
28	FORMAT(/5X,'OUTPUT DATA FILE: ',A)
C
	WRITE(IPRT,29)STARS,STARS
29	FORMAT(/A25,'**** REFLECTION DATA *****',a25/)
C
C CHECK THAT WE HAVE THE NECESSARY INFORMATION FOR INPUT FILE FORMAT:
C
C IF SHELX, NEED UB MATRIX AND 
C  IBL TO TELL US ORIENTATION OF DCOSINES W.R.T. BUSING-LEVY AXES
C
	IF(IDFORM_IN .eq. 2 .OR. IDFORM_IN .EQ. 3)THEN
C
C FIRST CHECK FOR UB MATRIX PRESENT
		IF(IUB .NE. 3)THEN
			WRITE(ESTRING,25)
25			FORMAT('UB MISSING FROM EXP FILE - NEEDED ',
	1				'TO CONVERT DIRECTION COSINES IN HKL FILE')
			CALL ERROR(ESTRING)
			GOTO 9999
		ENDIF
C
C NOW CHECK FOR BLAXES CARD
		IBLTEST=0
		DO I=1,3
			IBLTEST=IBLTEST+IABS(IBL(I))
		ENDDO
		IF(IBLTEST .NE. 6)THEN
			WRITE(ESTRING,40)
40			FORMAT('NEED BLAXES CARD TO INTERPRET DIRECTION COSINES ',
	1'ON REFLECTION DATA FILE, BUT BLAXES HAS INVALID INFORMATION OR ',
     1'IS ABSENT')
			CALL ERROR(ESTRING)
			GOTO 9999
		ENDIF
	ENDIF
C
C  COULD POP A MESSAGE BOX HERE TO SAY DATA BEING PROCESSED
	CALL PROCESS_DATA 
	IF(IFATAL .EQ. 1)GOTO 9999
C
C CLEAN UP
C
C	retlog= DlgSet(dlg, idc_message,
C	1'PROCESSING COMPLETE WRITING CIF FILE',DLG_TITLE)

	CALL WRITE_CIF_ABSORB
	IF(IFATAL .EQ. 1)GOTO 9999

C	retlog= DlgSet(dlg, idc_message,
C	1'ABSORB FINISHED SUCCESSFULLY',DLG_TITLE)

C SUCCESSFUL TERMINATION
C
C ALL EXITS THROUGH HERE
9999	CALL CLOSE_PROGRAM
9998	END

C
C****************************************************************************
C

	SUBROUTINE OPEN_ABSORB_FILES(IERR)
c
c THIS IS TO BE USED ONLY FOR ABSORB STANDALONE
C
C IT RESETS ISEL TO 0 IF THERE IS AN ERROR
C
	CHARACTER*12 F(4)
	CHARACTER*20 TITLE
	character*512 Xdir,sname,ext
	INTEGER I
	INCLUDE 'FILES.INC'
	character*512 :: file = ""C
	character*(*),parameter :: filter = 
	1"All Files"C//"*.*"C//""C
	DATA F/'experiment','print','input data','output data'/
C
c INITIALISATION
1	IERR=0			! WILL BE USED TO INDICATE ERROR FILE
	isel=0

C
C EXP FILE
	OPEN(UNIT=icnt,FILE=FNAME(1)(1:ILEN(1)),
	1			ACTION='READ',STATUS='OLD',err=10)
      isel(1)=1
C
C PRINT FILE
      OPEN(UNIT=IPRT,FILE=FNAME(2)(1:ILEN(2)),STATUS='UNKNOWN',
	1		ACCESS='APPEND',CARRIAGECONTROL='FORTRAN',ERR=20)
      isel(2)=1
C
C INPUT DATA FILE
      OPEN(UNIT=IIN,FILE=FNAME(3)(1:ILEN(3)),
	1			STATUS='OLD',ACTION='READ',ERR=30)
      isel(3)=1
C
C OUTPUT DATAFILE
      OPEN(UNIT=IOUT,FILE=FNAME(4)(1:ILEN(4)),STATUS='REPLACE',ERR=40)
      isel(4)=1
C
C HERE IF SUCCESS
	RETURN
C
C ERRORS
10	IERR=1
	RETURN
20	IERR=2
	RETURN
30	IERR=3
	RETURN
40	IERR=4
	RETURN
	END
C
C****************************************************************************
C
	SUBROUTINE CLOSE_PROGRAM
	INCLUDE 'FILES.INC'

C
C CLOSE MESSAGEBOX
C	IF(IABS_RUN .NE. 0)CALL CLOSEMESSAGEBOX
C LAST ACTION IS TO DELETE THE LOCK FILE
	CALL LOCK(0,ierr)
      if(ierr .ne. 0)then
          write(estring,'(''Lock file deleted during Absorb run: '',
     1    ''your results may be corrupted'')')
          call error(estring)      
      endif    
C
C CLOSE FILES
C
	CLOSE(UNIT=ICNT)
	CLOSE(UNIT=IPRT)
	CLOSE(UNIT=IIN)
	CLOSE(UNIT=IOUT)
	CLOSE(UNIT=IFERR)

      
      RETURN
	END
C
C*****************************************************************************
C
	SUBROUTINE DELETE_ERROR_LOG
	USE KERNEL32
	USE DFLIB
	CHARACTER*255 NAME
C
	ILEN=GETMODULEFILENAME(NULL,NAME,255)
	IF(ILEN .EQ. 0)RETURN
C
	I=INDEX(NAME,'\',BACK=.TRUE.)			! BUILD EXE DIR NAME
	NAME=NAME(1:I)//'absorb_error.log'
	I=LEN_TRIM(NAME)
	
	OPEN(UNIT=90,FILE=NAME(1:I),status='unknown')
	CLOSE(UNIT=90,STATUS='DELETE')
C
	RETURN
	END


C
C*****************************************************************************
C
	SUBROUTINE LOCK(ILOCK,IERR)
	USE KERNEL32
	USE DFLIB
      use ifport
	CHARACTER*255 NAME
      integer ierr,i
      save name,i

C
	INTEGER ILOCK
C
C
      ierr=0
	IF(ILOCK .EQ. 1)THEN
C
		ILEN=GETMODULEFILENAME(NULL,NAME,255)
		IF(ILEN .EQ. 0)RETURN
C
		I=INDEX(NAME,'\',BACK=.TRUE.)			! BUILD EXE DIR NAME
		NAME=NAME(1:I)//'absorb.lck'
		I=LEN_TRIM(NAME)
C
c		OPEN(UNIT=99,FILE=NAME(1:I),status='unknown',share='DENYRW')
c documentation says share=denywr is default
c might want to add error handling later
      open(unit=99,file=name(1:i),status='new',err=10,share='denynone')
c
c here if ok
          close(unit=99,status='keep')          
      ELSE
          open(unit=99,file=name(1:i),status='old',err=10)
		CLOSE(UNIT=99,STATUS='DELETE')          
	ENDIF
	RETURN
c errors
10    ierr=1
      return
      
      
	END

C*****************************************************************************
       block data filestreams
        include 'files.inc'
	data icnt,iprt,iin,iout,IPAR/3,4,1,2,9/
        end
C
c***************************************************************************
	SUBROUTINE OPEN_GRID
	INCLUDE 'FILES.INC'
	CHARACTER*512 NAME	
C
C OPENS GRID LISTING IN SAME DIRECTORY AS PRINT FILE
C
	I=LEN_TRIM(DIR(2))
	NAME=DIR(2)(1:I)//'ABSORB_GRID.LST'

	OPEN(UNIT=15,FILE=NAME,STATUS='REPLACE',
	1	CARRIAGECONTROL='FORTRAN',err=20)
	RETURN
cc here if absorb_grid.lst locked
20	WRITE(Estring,'('' ABSORB_GRID.LST file cannot be opened,'',
	1		''no grid listing will be written to this file'')')
	call WARNING(Estring)
	IGRID=0		! STOP GRID WRITING IN FILE NOT OPENED (EG SET READ ONLY)

	RETURN
	END



C******************************************************************************
C
C  UTILITY ERROR AND WRITING ROUTINES
C
C******************************************************************************
C



	SUBROUTINE MESSAGE(TEXT)
	INCLUDE 'FILES.INC'
	INCLUDE 'FLAGS.INC'
	character*256 text
C
C

C
C NOW POP UP ERROR MESSAGE
	IF(IABS_RUN .EQ. 1)CALL MESSAGEBOX(TEXT)
	RETURN
	END	


	SUBROUTINE WARNING(TEXT)
	INCLUDE 'FILES.INC'
	INCLUDE 'FLAGS.INC'
	character*256 text

C  WRITE ERROR TO PRINT FILE
C
C
	L=LEN_TRIM(TEXT)
	IF(L .GT. 100)L=100
	IF(ISEL(2) .EQ. 1)THEN
	WRITE(IPRT,10)TEXT(1:L)
10	FORMAT(//'*****WARNING:'/5X,A//'*****PROGRAM WILL CONTINUE'/)
	ENDIF
C
C WRITE ERROR TO ERROR_LOG
	CALL OPEN_ERROR_LOG
	WRITE(IFERR,10)TEXT(1:L)
      CLOSE(UNIT=IFERR)

C
C NOW POP UP ERROR MESSAGE
	IF(IABS_RUN .EQ. 1)CALL WARNINGBOX(TEXT)
	RETURN
	END
C
	SUBROUTINE ERROR(TEXT)
	INCLUDE 'FILES.INC'
	INCLUDE 'FLAGS.INC'
	character*256 text
C
C
	I=LEN_TRIM(TEXT)
	IF(I .GT. 100)I=100
C
C  WRITE ERROR TO PRINT FILE
C
	IF(ISEL(2) .EQ. 1)THEN
		WRITE(IPRT,10)TEXT(1:I)
10	FORMAT(//'*****FATAL ERROR:'/5X,A//'*****PROGRAM WILL HALT'/)
	ENDIF
C
C WRITE ERROR TO ERROR_LOG
	CALL OPEN_ERROR_LOG
	WRITE(IFERR,10)TEXT(1:I)
      CLOSE(UNIT=IFERR)

C
C SET FATAL ERROR FLAG IN COMMON
	IFATAL=1

C
C NOW POP UP ERROR MESSAGE
	IF(IABS_RUN .EQ. 1)CALL ERRORBOX(TEXT)
	RETURN
	END
C
C**********************************************************************************************************
C
	SUBROUTINE OPEN_ERROR_LOG
	USE KERNEL32
	USE DFLIB
	CHARACTER*255 NAME
	CHARACTER*256 TEXT
	INTEGER NERR,I
	INCLUDE 'FILES.INC'
	SAVE NERR,NAME,I

C
C  
C
	IF(NERR .EQ. 0)THEN
		I=GETMODULEFILENAME(NULL,NAME,255)
		IF(I .EQ. 0)RETURN
C
		I=INDEX(NAME,'\',BACK=.TRUE.)			! BUILD EXE DIR NAME
		NAME=NAME(1:I)//'absorb_error.log'
		I=LEN_TRIM(NAME)
		OPEN(UNIT=IFERR,FILE=NAME(1:I),status='NEW',
	1		CARRIAGECONTROL='FORTRAN')
C
C WRITE THE HEADER
C
		CALL WRITE_HEADER(IFERR)
      ELSE
		OPEN(UNIT=IFERR,FILE=NAME(1:I),status='OLD',access='APPEND',
	1		CARRIAGECONTROL='FORTRAN')
          
	ENDIF
	NERR=NERR+1
	
C
	RETURN
	END
C
C**********************************************************************************************************
C

	SUBROUTINE PROCESS_INDICATOR(i)
	INTEGER I
	RETURN
	END
C
C**********************************************************************************************************
C
	SUBROUTINE FIX_FILEFORMAT(IERR)
C
C THIS TRIES TO FIX PROBLEMS IF AN INVALID FILE FORMAT HAS BEEN GIVEN IN ABSORB_CONFIG,DAT
C
	INCLUDE 'FILES.INC'
C
C ON ENTRY IERR IS THE VALUE FROM READ_CONFIG
C
	IF(IERR .EQ. 3)THEN
		IF(INDEX(FNAME(3),'hkl') .ne. 0)THEN
			WRITE(ESTRING,200)IDFORM_IN
200			FORMAT('INPUT FORMAT NUMBER ',I5,
	1			' NOT RECOGNISED - SET TO HKL FILE')
			CALL WARNING(ESTRING)
			IDFORM_IN=2	
			IERR=0	
		ELSE
			WRITE(ESTRING,210)IDFORM_IN
210			FORMAT('INPUT FORMAT NUMBER ',I5,
	1			' NOT RECOGNISED')
			CALL ERROR(ESTRING)
			RETURN
		ENDIF
	ELSE	
C
		IF(INDEX(FNAME(4),'hkl') .ne. 0)THEN
			WRITE(ESTRING,220)IDFORM_OUT
220			FORMAT('OUTPUT FORMAT NUMBER ',I5,
	1			' NOT RECOGNISED - SET TO HKL FILE')
			CALL WARNING(ESTRING)
			IDFORM_OUT=4		
			IERR=0
		ELSE
			WRITE(ESTRING,230)IDFORM_OUT
230			FORMAT('OUTPUT FORMAT NUMBER ',I5,
	1			' NOT RECOGNISED')
			CALL ERROR(ESTRING)
		ENDIF
	ENDIF	
	RETURN
	END
