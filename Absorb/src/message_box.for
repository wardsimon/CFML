C**************************************************************************
C 
C   SUBROUTINES TO POP BOXES
C
C
C**************************************************************************



C**************************************************************************


	SUBROUTINE ERRORBOX(TEXT)
c    
      USE DFLOGM
	IMPLICIT NONE
	INCLUDE 'RESOURCE.FD'
	logical retlog
	INTEGER I,RETINT
	character*256 text
      TYPE (dialog) DLG,dlg_error
	INTEGER CONTROL_NAME,CALLBACKTYPE
	COMMON/TESTDLG/dlg, control_name, callbacktype
	


	RETLOG=DLGINIT(IDD_ERROR,DLG_ERROR)	! INITIALISE BOX
C
C
	I=LEN_TRIM(TEXT)
C
C WRITE STUFF INTO INFO BOX
	RETLOG=DLGSET(DLG_ERROR,idC_ERROR1,text(1:I),DLG_TITLE)
C Activate the modal dialog.
      retint = DlgModal( dlg_ERROR )
C
C     Release dialog resources.
      CALL DlgUninit( dlg_ERROR )

c
c

	RETURN
	END
C
	SUBROUTINE WARNINGBOX(TEXT)
c    
      USE DFLOGM
	IMPLICIT NONE
	INCLUDE '..\code\FILES.INC'
	INCLUDE 'RESOURCE.FD'
	logical retlog
	INTEGER L,RETINT
	character*256 text
      TYPE (dialog) dlg_WARN
	INTEGER CONTROL_NAME,CALLBACKTYPE
	
	L=LEN_TRIM(TEXT)

	RETLOG=DLGINIT(IDD_WARNING,DLG_WARN)	! INITIALISE BOX
C
C WRITE STUFF INTO INFO BOX
	RETLOG=DLGSET(DLG_WARN,idC_WARNTEXT,text(1:L),DLG_TITLE)
C Activate the modal dialog.
      retint = DlgModal( dlg_WARN )
C
C     Release dialog resources.
      CALL DlgUninit( dlg_WARN )
	RETURN
	END
C
C**********************************************************************
C
C MESSAGE BOX IS DIFFRENT - NEEDS TO REMAIN UP...THESE DO NOT YET WORK!!!
C
	SUBROUTINE CREATEMESSAGEBOX
c    
      USE DFLOGM
	IMPLICIT NONE
	INCLUDE 'RESOURCE.FD'
	logical retlog
	INTEGER I,RETINT
	character*256 text
      TYPE (dialog) DLG,dlg_MESSAGE
	INTEGER CONTROL_NAME,CALLBACKTYPE
	COMMON/TESTDLG/dlg, control_name, callbacktype
	


	RETLOG=DLGINIT(IDD_MESSAGE,DLG_MESSAGE)	! INITIALISE BOX
C
C
      retint = DlgModELESS( dlg_message )

c
	RETURN
	END
	SUBROUTINE MESSAGEBOX(TEXT)
c    
      USE DFLOGM
	IMPLICIT NONE
	INCLUDE 'RESOURCE.FD'
	logical retlog
	INTEGER I,RETINT
	character*256 text
      TYPE (dialog) DLG,dlg_MESSAGE
	INTEGER CONTROL_NAME,CALLBACKTYPE
	COMMON/TESTDLG/dlg, control_name, callbacktype
	
C
	I=LEN_TRIM(TEXT)
C
C WRITE STUFF INTO INFO BOX
	RETLOG=DLGSET(DLG_MESSAGE,IDC_MESSAGE,text(1:I),DLG_TITLE)
      

C Activate the modal dialog.
      retint = DlgModal( dlg_MESSAGE)
C
C     Release dialog resources.
      CALL DlgUninit( dlg_MESSAGE)

      RETURN
	END
	SUBROUTINE CLOSEMESSAGEBOX
c    
      USE DFLOGM
	IMPLICIT NONE
	INCLUDE 'RESOURCE.FD'
	logical retlog
	INTEGER I,RETINT
      TYPE (dialog) DLG,dlg_MESSAGE
	INTEGER CONTROL_NAME,CALLBACKTYPE
	COMMON/TESTDLG/dlg, control_name, callbacktype
C
C     Release dialog resources...
      CALL DlgUninit( dlg_MESSAGE )
c
	RETURN
	END



C**********************************************************************


	SUBROUTINE GET_FORMAT(IN)
      USE DFLOGM
	IMPLICIT NONE
	INCLUDE '..\code\FILES.INC'
	INCLUDE 'RESOURCE.FD'
	logical retlog,pushed_state
	INTEGER INNN,IN,J,JS,RETINT
	EXTERNAL RADIOFORMAT
	CHARACTER*256 TEXT
      TYPE (dialog) dlg_format
	INTEGER CONTROL_NAME,CALLBACKTYPE
	COMMON/LOCALIN/INNN

	INNN=IN
C
C IN=3 FOR INPUT DATAFILE, 4 FOR OUTPUT DATA FILE	
C
	RETLOG=DLGINIT(IDD_format,DLG_format)	! INITIALISE BOX
C
C
C WRITE FILE NAME INTO INFO BOX
	J=INDEX(FNAME(IN),'   ')-1
	JS=1
	IF(J .GT. 40)JS=J-39
	RETLOG=DLGSET(DLG_FORMAT,idC_WARNFILENAME,FNAME(in)(JS:J),DLG_TITLE)
C
C SET CORRECT FORMAT ON RADIO BUTTON
C
	IF(IN .EQ. 3)THEN
	RETLOG=DLGSET(DLG_FORMAT,idC_FORMATRADIO1,'WinIntegrStp',DLG_TITLE)
	ELSE
	RETLOG=DLGSET(DLG_FORMAT,idC_FORMATRADIO1,
	1					'Rfine abs/avg',DLG_TITLE)
	ENDIF
C
c set incommensurate button not active
	retlog = DLGSET (dlg_FORMAT, IDC_incomm2, .false.)


c
c set the button callbacks
c
        retlog = DlgSetSub( dlg_format, IDC_formatRadio1, radioformat)
        retlog = DlgSetSub( dlg_format, IDC_formatRadio2, radioformat)
c


C Activate the modal dialog.
      retint = DlgModal( dlg_format )



C     Release dialog resources.
      CALL DlgUninit( dlg_format )
c
c

	RETURN
	END
c
	subroutine radioformat( dlg_FORMAT, control_name, callbacktype )
      USE DFLOGM
	IMPLICIT NONE
	INCLUDE '..\code\FILES.INC'
	INCLUDE 'RESOURCE.FD'
	logical retlog,PUSHED_STATE,CHECKED_STATE
	INTEGER IN,CONTROL_NAME,CALLBACKTYPE
	character*16 text
      TYPE (dialog) dlg_FORMAT
	COMMON/LOCALIN/IN
C
C read the buttons: changed from test on retlog to pushed_state 3 feb 2012
C
	retlog = DLGGET (dlg_format, IDC_formatRADIO1, pushed_state)
	if(pushed_state)then
		if(in .eq. 3)then
			idform_in=1
		else
			idform_out=1
		endif
	else
		retlog = DLGGET (dlg_format, IDC_formatRADIO2, pushed_state)
		if(pushed_state)then
			if(in .eq. 3)then
				idform_in=2
			else
				idform_out=2
			endif
		else
			if(in .eq. 3)then
				idform_in=3
			else
				idform_out=3
			endif
		endif
	endif
c
c read the incomm button
	retlog = DLGGET (dlg_format, IDC_incomm1, CHECKED_state)
	if(checked_state .eq. .true.)then
		icflag=1
	else
		icflag=0
	endif

	RETURN
	END