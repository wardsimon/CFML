c
c     information for setting up crystal grid
c

       REAL XPT(32),YPT(32,32),ZPT(32,32,32),XW(32),YW(32,32),
     1       ZW(32,32,32)
	   real gauss(6,32),gweight(6,32)

	   COMMON/DEFGRID/XPT,YPT,ZPT,XW,YW,ZW
	   COMMON/GCOEFFS/GAUSS,GWEIGHT
