C 
C -------------------------------------------------------------
	SUBROUTINE CALCMU(CHEMISTRY)		
C	
C SUBROUTINE TO INTERPRET CHEMISTRY STRING AND RETURN ABSCO IN Micron-1 
C 
C CHEMISTRY SHOULD BE UPPER CASE ELEMENT SYMBOLS PLUS NUMBERS: SI1 O2
C  CAN ALSO HAVE Z= TO SET NO FORM UNITS/CELL
C 
C   IABSCO=0 IF FAIL
C	
	INCLUDE 'CRYSTALMODEL.INC'		! NEEDS UNITCELL VOL (CVOL), AND SETS ABSCO
	INCLUDE 'UNITCELL.INC'			! NEEDS WAVELENGTH
	include 'chemistry.inc'			! SAVE CHEMISTRY INFO
C
	INTEGER ISING(92)	! FLAG TO SAY SINGLE CHARACTER IN ELEMENT SYMBOL
C
	REAL AtomicWeight(92)	! ATOMIC WEIGHTS
	REAL MASSABS(6,92)		! MASS ABSORPTION COEFFS BY WAVELENGTH: Ag added 6-April-2012 as #6
C
	CHARACTER CHEMISTRY*60, CHEM*60, ETEXT*256
C
C*******************ELEMENT SYMBOLS*************************
      data ElementString(1:46)/
     1'H HELIBEB C N O F NENAMGALSIP S CLARK CASCTIV '/
      data ElementString(47:92)/
     1'CRMNFECONICUZNGAGEASSEBRKRRBSRY ZRNBMOTCRURHPD'/
      data ElementString(93:138)/
     1'AGCDINSNSBTEI XECSBALACEPRNDPMSMEUGDTBDYHOERTM'/
      data ElementString(139:184)/
     1'YBLUHFTAW REOSIRPTAUHGTLPBBIPOATRNFRRAACTHPAU '/
C***********************************************************
      data AtomicWeight /
     1      1.007,4.000,6.939,9.012,10.81,12.01,14.01,16.00,
     2      19.00,20.18,22.99,24.31,26.98,28.09,30.97,32.06,
     3      35.45,39.95,39.10,40.08,44.96,47.90,50.94,52.00,
     4      54.94,55.85,58.93,58.71,63.54,65.37,69.72,72.59,
     5      74.92,78.96,79.91,83.80,85.47,87.62,88.90,91.22,
     6      92.91,95.94,98.00,101.10,102.91,106.4,107.87,112.40,
     7      114.82,118.69,121.75,127.60,126.90,131.30,132.91,137.34,
     8      138.91,140.12,140.91,144.24,147.00,150.35,152.00,
     9      157.25,158.92,162.50,164.93,167.26,168.93,173.04,
     1      174.97,178.49,180.95,183.85,186.2,190.2,192.2,
     2      195.09,196.97,200.59,204.37,207.19,208.98,210.0,
     3      210.0,222.0,223.0,226.0,227.0,232.04,231.0,238.03/
C*************************************************************
C
C Mass attenuation coefficients from I Tables Vol C 1992
c  Values are in units of cm2.g-1
c
      data (MassAbs(1,I),I=1,92) /
     1      0.412,0.498,1.30,3.44,7.59,15.0,24.7,37.8,    !{Cr Ka}
     2      51.5,74.1,94.9,126.,155.,196.,230.,281.,
     3      316.,342.,421.,490.,516.,590.,74.7,86.8,
     4      97.5,113.,124.,144.,153.,171.,183.,199.,
     5      219.,234.,260.,277.,303.,328.,358.,386.,
     6      416.,442.,474.,501.,536.,563.,602.,626.,
     7      663.,691.,723.,740.,796.,721.,760.,570.,
     8      225.,238.,238.,251.,294.,279.,309.,298.,
     9      332.,325.,347.,352.,386.,387.,431.,425.,
     1      432.,457.,501.,499.,520.,541.,551.,541.,
     2      597.,643.,666.,691.,680.,734.,758.,743.,
     3      739.,768.,738.,766./
      data (MassAbs(2,I),I=1,92) /
     1      0.400,0.381,0.839,2.09,4.55,8.99,14.9,22.8,    !{Fe Ka}
     2      31.3,45.2,58.2,77.8,95.9,122.,144.,177.,
     3      200.,218.,270.,314.,332.,358.,399.,492.,
     4      61.6,71.0,78.5,91.3,96.8,108.,116.,127.,
     5      139.,149.,165.,176.,193.,210.,229.,247.,
     6      267.,284.,305.,323.,346.,363.,389.,405.,
     7      428.,447.,471.,483.,522.,540.,569.,586.,
     8      618.,561.,448.,455.,194.,204.,203.,195.,
     9      219.,214.,228.,232.,253.,251.,280.,277.,
     1      283.,301.,327.,327.,340.,357.,361.,339.,
     2      403.,420.,434.,452.,444.,477.,493.,487.,
     3      530.,485.,482.,528./
      data (MassAbs(3,I),I=1,92)/
     1      0.397,0.343,0.693,1.67,3.59,7.07,11.7,18.0,    !{Co Ka-bar}
     2      24.7,35.8,46.2,61.9,76.4,97.8,115.,142.,
     3      161.,176.,218.,255.,269.,291.,325.,408.,
     4      393.,57.2,63.2,73.5,78.0,87.1,93.4,102.,
     5      112.,120.,133.,142.,156.,170.,185.,200.,
     6      216.,230.,247.,262.,280.,295.,316.,329.,
     7      349.,364.,383.,394.,425.,440.,465.,480.,
     8      507.,535.,565.,505.,400.,176.,419.,161.,
     9      180.,176.,187.,191.,206.,206.,229.,227.,
     1      231.,246.,268.,268.,278.,276.,295.,273.,
     2      331.,343.,355.,370.,363.,392.,403.,398.,
     3      461.,406.,394.,420./

      data (MassAbs(4,I),I=1,92)/
     1      .391,.292,.500,1.11,2.31,4.51,7.44,11.5,      !{Cu Ka-bar}
     2      15.8,22.9,29.7,40.0,49.6,63.7,75.5,93.3,
     3      106.,116.,145.,170.,180.,200.,219.,247.,
     4      270.,302.,321.,48.8,51.8,57.9,62.1,67.9,
     5      74.7,80.0,89.0,95.2,104.,113.,124.,139.,
     6      145.,154.,166.,176.,189.,199.,213.,222.,
     7      236.,247.,259.,267.,288.,299.,317.,325.,
     8      348.,368.,390.,404.,426.,434.,434.,403.,
     9      321.,362.,129.,132.,140.,142.,156.,155.,
     1      158.,168.,187.,184.,191.,188.,201.,188.,
     2      226.,235.,244.,254.,248.,267.,277.,273.,
     3      317.,306.,271.,288./

      data (MassAbs(5,I),I=1,92)/
     1      0.373,0.202,0.198,0.256,0.368,0.576,0.845,1.22,    !{Mo Kalpha bar}
     2      1.63,2.35,3.03,4.09,5.11,6.64,7.970,9.99,
     3      11.5,12.8,16.20,19.30,20.8,23.4,26.0,29.9,
     4      33.1,37.6,41.0,46.9,49.1,54.0,57.0,61.2,
     5      66.1,69.5,75.6,79.30,85.10,90.6,97.0,16.30,
     6      17.7,18.8,20.4,21.7,23.3,24.7,26.5,27.8,
     7      29.5,31.0,32.7,33.8,36.7,38.2,40.7,42.3,
     8      44.9,47.7,50.7,53.0,56.3,57.8,60.9,62.6,
     9      65.8,68.3,71.3,74.4,77.9,80.4,84.0,86.9,
     1      90.4,93.8,97.4,100.0,104.0,107.0,112.0,115.0,
     2      118.0,122.0,126.0,132.0,117.0,108.0,87.0,88.0,
     3      90.8,96.5,101.0,102.0/
c
c each line is one block of Table 4.3.4.3 , page 230- of ITables vol C
c
      data (MassAbs(6,I),I=1,92)/
     1	0.367,0.193,0.179,0.209,0.267,0.374,0.503,0.685,	!{Ag Kalpha bar}
	2	0.879,1.23,1.56,2.09,2.59,3.35,4.01,5.02,
	3	5.79,6.46,8.19,9.79,10.6,11.9,13.3,15.4,
	4	17.0,19.4,21.2,24.4,25.6,28.2,29.8,32.1,
	5	34.8,36.8,40.3,42.5,45.9,49.1,52.9,55.9,
	6	59.8,72.0,66.0,11.4,12.3,13.0,14.0,14.6,
	7	15.6,16.4,17.3,17.9,19.4,20.2,21.5,22.4,
	8	23.8,25.3,26.9,28.1,29.9,30.8,32.4,33.4,
	9	35.1,36.5,38.1,39.9,41.8,43.2,45.1,46.7,
	1	48.7,50.5,52.5,54.1,56.3,58.3,60.7,62.6,
	1	64.5,66.6,69.1,72.3,75.1,72.1,77.0,79.3,
	2	82.4,83.9,87.8,88.6/
C
C
C***********************************************************************
C 
C SETUP
C
C   SETTING FLAGS FOR SINGLE ELEMENT FLAGS
C
	DO I=1,92
		ISING(I)=0
	ENDDO
	ISING(1)=1		! H      
	ISING(5)=1		! B			
	ISING(6)=1		! C			
	ISING(7)=1		! N	
	ISING(8)=1		! O	
	ISING(9)=1		! F	
	ISING(15)=1		! P			
	ISING(16)=1		! S			
	ISING(19)=1		! K
	ISING(23)=1		! V
	ISING(39)=1		! Y
	ISING(53)=1		! I
	ISING(74)=1		! W
	ISING(92)=1		! U
C
C  GET WAVELENGTH

	IF(IWAVE .EQ. 0)THEN
		WRITE(ETEXT,10)
10		FORMAT('NO WAVELENGTH SELECTED: ',
	1		   ' ABS COEFF CANNOT BE CALCULATED')
		CALL WARNING(ETEXT)
		IABSCO=0
		ABSCO=0.
		RETURN
	ENDIF
C
C LOCAL COPY OF INPUT STRING
	CHEM=CHEMISTRY				
C
C FIRST LOCATE AND INTERPRET Z
C
	Z=1
	Js=INDEX(CHEM,'Z=')
      if(js == 0)then
          js=index(CHEM,'Z =')     ! added 22-Dec-2014
          j=js+3
      else
          j=js+2
      endif
      do
            if(chem(j:j) /=' ')exit
            j=j+1     !js is now the loc(Z) and j where the number should start
      enddo
      
	IF(J .GT. 2)THEN
		K=INDEX(CHEM(J:60),' ')-1		! LOOKING FOR TERMINATOR
		L=INDEX(CHEM(J:60),',')-1
		JSIZ=MAX(K,L)		
		IF(JSIZ .NE. 0)THEN
			READ(CHEM(J:J+JSIZ),*)Z
			CHEM(JS:J+JSIZ)='          '	!CLEAR
		ENDIF
	ENDIF
C
C NOW SEARCH THE STRING BY LOOPING OVER THE ELEMENT SYMBOLS FOR DOUBLE SYMBOLS
C
	DO I=1,92
		NE(I)=0.			!INITIALISATION
		IF(ISING(I) .EQ. 0)THEN
			J=INDEX(CHEM,ELEMENTSTRING(2*I-1:2*I))+2
			IF(J .GT. 2)THEN
				L=INDEX(CHEM(J:60),',')-1 		! LOOKING FOR TERMINATOR as comma
				IF(L .LT. 1)L=60				! NO COMMA AT ALL 
C
C							! INSTEAD OF PREVIOUS LOOK FOR BLANK, LOOK FOR NEXT CHAR
C
				K=J-1
20				K=K+1
				IF(ICHAR(CHEM(K:K)) .GT. 64 .AND. 
	1				ICHAR(CHEM(K:K)) .LT. 91)GOTO 30	! character
				IF(K .LT. 60)GOTO 20
30				K=K-J-1			! K NOW POINTS TO 1 BEFORE NEXT CHAR
				JSIZ=MIN(K,L)		
				NE(I)=1.0				!NO NUMBER MEANS 1
				IF(JSIZ .NE. 0)THEN
					READ(CHEM(J:J+JSIZ),*,END=40)NE(I)
				ENDIF
40				CHEM(J-2:J+JSIZ)='  '		!CLEAR SYMBOL AND NUMBER
			ENDIF
		ENDIF
	ENDDO
C
C REPEAT LOOP FOR SINGLE CHARACTER ELEMENT SYMBOLS
	DO I=1,92
		IF(ISING(I) .EQ. 1)THEN
			J=INDEX(CHEM,ELEMENTSTRING(2*I-1:2*I-1))+1
			IF(J .GT. 1)THEN
				L=INDEX(CHEM(J:60),',')-1 		! LOOKING FOR TERMINATOR as comma
				IF(L .LT. 1)L=60				! NO COMMA AT ALL 
C
C							! INSTEAD OF PREVIOUS LOOK FOR BLANK, LOOK FOR NEXT CHAR
C
				K=J-1
120				K=K+1
				IF(ICHAR(CHEM(K:K)) .GT. 64 .AND. 
	1				ICHAR(CHEM(K:K)) .LT. 91)GOTO 130	! character
				IF(K .LT. 60)GOTO 120
130				K=K-J-1			! K NOW POINTS TO 1 BEFORE NEXT CHAR
				JSIZ=MIN(K,L)		
				NE(I)=1.0				!NO NUMBER MEANS 1
				IF(JSIZ .NE. 0)THEN
					READ(CHEM(J:J+JSIZ),*,END=140)NE(I)
				ENDIF
140				CHEM(J-1:J+JSIZ)='  '		!CLEAR SYMBOL AND NUMBER

			ENDIF
		ENDIF
	ENDDO
C
C SEARCH STRING FOR NON-NUMERICAL JUNK
	DO I=1,60
		IF(CHEM(I:I) .NE. ' ')THEN
		IF(ICHAR(CHEM(I:I)) .LT. 48 .OR. ICHAR(CHEM(I:I)) .GT. 57)THEN
		 IF(CHEM(I:I) .NE. ',' .AND. CHEM(I:I) .NE. '.')THEN
		WRITE(ETEXT,200)
200	FORMAT('UNINTERPRETED CHARACTER ON CELL CONTENTS CARD')
		CALL WARNING(ETEXT)
		ENDIF
		ENDIF
	ENDIF
	ENDDO
C
C LOOP OVER ELEMENT LIST AND ACCUMULATE SUM FOR TOTAL WEIGHT
C
	WEIGHTSUM=0.
      DO I= 1,92
	  IF(NE(I) .GT. 0.00001)THEN
          NoElements = NoElements + 1
          WeightSum = WeightSum + NE(I)* AtomicWeight(I)
        END IF
	ENDDO
C
C NOW THE SUMS TO GET MU
	SUMMUL=0.
     	DO I= 1,92
	  IF(NE(I) .NE. 0.00001)THEN	  
		SUMMUL = SUMMUL + NE(I)* AtomicWeight(I) 
	1				* MassAbs(IWave,I)
        END IF
	ENDDO
      Density = WeightSum*Z/(0.6023*CVol)
      ABSCO = SumMuL * Density/WEIGHTSUM

	IF(ABSCO .LT. 0.01)THEN
		ABSCO=0.
		IABSCO=0
		WRITE(ETEXT,202)
202		FORMAT('ABSCO CALCULATED NEGATIVE OR ZERO',
	1			' - SET TO ZERO')
		CALL WARNING(ETEXT)
	ELSE
		ABSCO=ABSCO/10000.		! RESCALE TO MICRON-1
		IABSCO=1
	ENDIF
	RETURN
	END
