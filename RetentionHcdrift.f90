!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************   
!	Subroutine RetentionHdrift
!	Calcola la ritenzione iniziale dovuta alla vegetazione
!***********************************************************************************
SUBROUTINE RetentionHdrift(iRows,iCols)
!
implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
REAL*8 a2dMaxRetention(iRows,iCols),rainp(iRows,iCols),retp(iRows,iCols),err(iRows,iCols)


IF(dFlagLai.eq.0)THEN
	a2dMaxRetention = (0.038*a2dS+0.4909)*1
ELSE
	!a2dMaxRetention = a2dLai*0.4
	WHERE (a2dDem.GT.0.0)
		a2dMaxRetention = 0.95 + 0.5*a2dLai-0.006*a2dLai**2
	ENDWHERE
	WHERE (a2dMaxRetention.GT.50.0)a2dMaxRetention =0.0
ENDIF


retp=a2dRetention
rainp=a2dRain

WHERE ((a2dRetention+a2dRain).LE.a2dMaxRetention.AND.a2dDem.GT.0.0)
	a2dRetention=a2dRetention+a2dRain
	a2dRain=0.0
ENDWHERE	

WHERE((a2dRetention+a2dRain).GT.a2dMaxRetention.AND.a2dDem.GT.0.0)
	a2dRain=(a2dRetention+a2dRain)-a2dMaxRetention
	a2dRetention=a2dMaxRetention
ENDWHERE

err=rainp-(a2dRetention-retp)-a2dRain

a2dAreeWT=a2dRain	

WHERE (a2dRetention.lt.0.0)
	a2dRetention=0.0
ENDWHERE


END