C -------------------------------------------------------------
C
C ********************************************************************************
C
C  INTENSITY.FOR HOLDS ROUTINES FOR MODIFYING THE INCIDENT BEAM INTENSITY
C  FOR POINT TO POINT ACROSS GAUSSIAN GRIDS
C
C  DOES NOT DEPEND ON THE CHOICE OF MODEL DESCRIPTION (SPHERE, HKL ETC)
C
C
C  CREATED AND ADDED TO PROJECT 24 FEB 2012
C
C ********************************************************************************
C
	FUNCTION GET_INTENSITY(X)
c
C
	REAL X(3)
C
C RETURNS THE INCIDENT BEAM INTENSITY AT GRID COORDS X,Y,Z
C
	GET_INTENSITY=1.0
	RETURN
	END
C
C******************************************************************************
	SUBROUTINE SET_CONFIG
	INCLUDE 'FLAGS.INC'

C	
C  SETS ICONFIG=0 TO INDICATE BASE CONFIGURATTON REST OF PROGRAM
C
	ICONFIG=0
C
	RETURN
	END
C******************************************************************************

C
C  THIS SUBROUTINE CAN DO ANY PRE-PROCESSING THAT IS NECESSARY 
C  AFTER THE EXP FILE HAS BEEN READ AND BEFORE THE CRYSTAL MODEL IS BUILT
	SUBROUTINE PRE_PROCESS
	INCLUDE 'FLAGS.INC'
C  PRE_PROCESS SETS ICONFIG=0 TO SET NORMAL CONFIGURATION
	ICONFIG=0

	RETURN
	END
C
	SUBROUTINE READ_SPECIALS(STRING)
C
	CHARACTER*80 STRING
	RETURN
	END
C
C******************************************************************************
C
	SUBROUTINE WRITE_CRYSTAL_SPECIAL
C
C DUMMY
C
	RETURN
	END
C
C******************************************************************************
C
	SUBROUTINE WRITE_REFLECTION_SPECIALS
C
C DUMMY
C
	RETURN
	END
