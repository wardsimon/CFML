!
!     CRYSTALMODEL.INC 
!      Holds the common blocks with the info about the model of the crystal
!      but not the default grid stuff
!
	  character*80 title
	  common/exptitle/title

      real cvol,volam,fact,eps,ngrid(3),tthmin,tthmax,scale_i
      integer iabsco,ishad,npt(3),ichem,ISKIP_REJECT
      integer inegsig,inegint,imodelorigin
      common/control/cvol,volam,fact,eps,iabsco,ishad,iskip_reject,
     1       ichem,npt,ngrid,tthmin,tthmax,scale_i,inegsig,inegint,
     2       imodelorigin


	  real xmin(3),xmax(3),xori(3),b(100,4),xvol,absco,xtal_rad,
     1 numer(32,32,32,100),corner(200,3),hkld(100,4),model_origin(3)

	  integer i1,ncorn,iedge(200,2),nedge,imodelform


      common/model/i1,b,xmin,xmax,xvol,numer,corner,ncorn,iedge,
     1           nedge,absco,XTAL_RAD,imodelform,xori,hkld,model_origin

