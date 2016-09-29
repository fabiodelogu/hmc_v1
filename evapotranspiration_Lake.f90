!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!     Subroutine Evapotraspiration Lake
!***********************************************************************************

	SUBROUTINE evapotranspiration_Lake(iRows,iCols,d)
	implicit none
	
	include 'DeclarationH.f90'

	CHARACTER*12  xstring,sVariable
	integer i,j
	integer iRows,iCols,iD
	real*8 d     
	real*8 dEvlake
      


	!Aggiorno il volume dei laghi a monte delle dighe
	DO iD=1,iNdam
		i=a2dXYDam(iD,2)
		j=a2dXYDam(iD,1)
        !Equazione riempimento dei laghi
		IF(a1dCodeDam(iD).gt.0)THEN 
		!Lago distribuito
			dEvlake=SUM(SUM(a2dEvapot,DIM=1,MASK=a2iChoice.eq.a1dCodeDam(iD))) !DIM=1 columns
			!write(*,*)dEvlake,a1dDamVolume(iD)
			a1dDamVolume(iD)=a1dDamVolume(iD)-dEvlake/1000*a2dAreaCell(i,j) !in m^3
		ENDIF
		!write(*,*)a1dDamVolume(iD)
		if(a1dDamVolume(iD).lt.0.0)a1dDamVolume(iD)=0.0
		dEvlake=0
	ENDDO

	!Aggiorno il volume dei laghi senza dighe
	DO iD=1,iNlake
		i=a2dXYLake(iD,2)
		j=a2dXYLake(iD,1)
        IF(a1dCodeLake(iD).gt.0)THEN 
		!Lago distribuito
			dEvlake=SUM(SUM(a2dEvapot,DIM=1,MASK=a2iChoice.EQ.a1dCodeLake(iD)))
			a1dVlake(iD)=a1dVlake(iD)-dEvlake/1000*a2dAreaCell(i,j) !in m^3
		ENDIF
		if(a1dVlake(iD).lt.0.0)a1dVlake(iD)=0.0
	END DO


	RETURN

	END
      
!***********************************************************************************

