C
C    FILES.INC HAS THE UNIT NUMBERS AND FORMAT FOR I/O
C
C  FNAME HOLDS THE FULL PATH+FILE NAME OF THE FILES, LENGTH ILEN
C  DIR HOLDS THE PATH NAME WITH THE FINAL '/'
C  SHORTNAME HOLDS THE FILENAME WITHOUT EXTENSION
C
      integer icnt,iprt,iin,iout,igrid,ipar,idform_in,idform_out
      integer isel(4),nindex,icflag,ILEN(4),IABS_RUN,iferr
      character*512 fname(4),dir(4),shortname(4)
      character*256 estring
	  character*12 version_date
	  real version_no
	  common/files/icnt,iprt,iin,iout,igrid,ipar,iferr,idform_in,
     1   idform_out,isel,nindex,icflag,iabs_run,ILEN,fname,shortname,dir
      common/version/version_no,version_date

