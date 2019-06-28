C 
C -------------------------------------------------------------
      SUBROUTINE DATE(STRING)
      CHARACTER*11 STRING		!CHANGED TO 11 FOR Y2K
      CHARACTER*3 MONTH(12)
      CHARACTER*10 TTIME,DDATE,ZONE,TEMP
      INTEGER IVAL(8)
C
      DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL',
     1           'AUG','SEP','OCT','NOV','DEC'/
C
      CALL DATE_AND_TIME(DDATE,TTIME,ZONE,IVAL)
C
C Keep compiler quiet
      temp = ddate
      temp = ttime
      temp = zone
      WRITE(STRING,10)IVAL(3),MONTH(IVAL(2)),IVAL(1)
10    FORMAT(I2,'-',A3,'-',I4)
      RETURN
      END
C
C
      SUBROUTINE TIME(STRING)
      CHARACTER*8 STRING
      CHARACTER*10 TTIME,DDATE,ZONE,temp
      INTEGER IVAL(8)
C
      CALL DATE_AND_TIME(DDATE,TTIME,ZONE,IVAL)
C Keep compiler quiet
      temp = ddate
      temp = ttime
      temp = zone
C
      WRITE(STRING,10)IVAL(5),IVAL(6)
10    FORMAT(I2.2,':',I2.2)
      RETURN
      END