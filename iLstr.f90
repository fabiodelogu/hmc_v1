!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!   Funzione per il calcolo della lunhezza di una stringa
!***********************************************************************************

	integer function iLStr(a)
	character*500 a

	last=1
	
	do while (a(last:last).ne.' ')
	last=last+1
    enddo

	iLStr=last-1
	
	RETURN
	end 


