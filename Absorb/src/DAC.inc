C     DAC.INC INCLUDES DATA ON THE DAC AND THE GASKET
C     
C
      REAL psiabs(3,100),abspsi_par(2),tdiamond(2),tberyl(2),
     1mu_dia(2),mu_be(2),mu_media,psimax(2),phizero,tcrystal

      
	  integer idac,ipsiNEG,NTERMS_PSI,nobsc,npsireg,IMEDIA,ixanvil
      integer ismall


	  common/dac/idac,ismall,ipsiNEG,NTERMS_PSI,psiabs,ABSPSI_PAR,
     1          tdiamond,tberyl,mu_dia,mu_be,mu_media,psimax,
     2          nobsc,npsireg,IMEDIA,phizero,ixanvil,tcrystal


      real tg,rgasket,gasketmu,thresho
      integer iopaque

	  common/gasket/iopaque,tg,rgasket,gasketmu,thresho
C
C
      REAL TPHIDAC(3,3),TDACPHI(3,3)
      COMMON/DACPHI/TPHIDAC,TDACPHI

