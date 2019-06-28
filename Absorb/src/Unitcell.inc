C     
C     UNITCELL.INC HAS THE CELL, UB MATRIX AND THE GONIOMETER INFO
C
C
	  real recell(7),gstar(3,3),ub(3,3),ubti(3,3),ubtranspose(3,3),
     1 xlamb,AUB(7),ACELL(7)
	  integer iub,ibl(3),parity(4),icell,iwave
	  common/xxx/iub,icell,ub,recell,gstar,aub,acell,ubti,ubtranspose
	  common/gonio/parity,ibl,xlamb,iwave