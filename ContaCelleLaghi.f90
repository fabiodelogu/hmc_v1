!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!------------------------------------------------------------------
!	Conta le celle dei diversi laghi a monte delle dighe
!------------------------------------------------------------------

 SUBROUTINE ContaCelleLaghi(iRows,iCols)
	include 'DeclarationH.f90'


	integer*4 id
    integer*4 name_lenght
	REAL*8 dDintegr  !Intervallo di integrazione del routing in secondi
	real*8 a2dConta(iRows,iCols)

!------------------------------------------------------------------
!   Calcolo il numero di celle corrispondente ai laghi a monte della diga
	do id=1,iNdam

		a2dConta=0
		WHERE (a2iChoice.eq.a1dCodeDam(iD).AND.a2dDem.GT.0.0)
			a2dConta=1
		ENDWHERE
		a1dNumCodeDam(iD)=SUM(SUM(a2dConta,DIM=1,MASK=a2dConta.EQ.1)) !Numero di celle del lago
	enddo
!------------------------------------------------------------------
!   Calcolo il numero di celle corrispondente ai laghi
	do id=1,iNlake
		
		a2dConta=0
		WHERE (a2iChoice.eq.a1dCodeLake(iD).AND.a2dDem.GT.0.0)
			a2dConta=1
		ENDWHERE
		a1dNumCodeLake(iD)=SUM(SUM(a2dConta,DIM=1,MASK=a2dConta.EQ.1)) !Numero di celle del lago
	enddo
	RETURN
	END

