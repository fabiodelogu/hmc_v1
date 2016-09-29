!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
! Check Input Meteorological Map:
!		- Rain 
!		- Temperature 
!		- Relative Humidity 
!		- Wind
!		- Short Wave Radiation
!***********************************************************************************  
 
subroutine CheckMeteoMapBinary(xstring)

implicit none
include 'DeclarationH.f90'
real*8 rm,tm,km,um,wm
CHARACTER*12  xstring

!rm=SUM(SUM(a2dRain,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!tm=SUM(SUM(a2dTemp,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!km=SUM(SUM(a2dK,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!um=SUM(SUM(a2dUm,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!wm=SUM(SUM(a2dW,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!888 FORMAT (1000(A4,1x,f9.2))
!write(*,888)'Rm:',rm,' Tm:',tm,' Km:',km,' Um:',um,' Wm:',wm

WHERE (a2dDem.gt.0.0)
	WHERE(a2dRain.gt.850.or.a2dRain.lt.0.0)
		a2dRain=0
	ENDWHERE
	WHERE(a2dTemp.lt.-70)
		a2dTemp=-70
	ENDWHERE
	WHERE(a2dTemp.gt.60)
		a2dTemp=60
	ENDWHERE
	WHERE(a2dK.lt.0)
		a2dK=0
	ENDWHERE
    WHERE(a2dK.gt.1412)
		a2dK=1412
	ENDWHERE
    WHERE(a2dK.lt.0)
		a2dK=0
	ENDWHERE
	WHERE(a2dW.gt.80)
		a2dW=80
	ENDWHERE
	WHERE(a2dW.lt.0)
		a2dW=0
	ENDWHERE
	WHERE(a2dUm.gt.1)
		a2dUm=1
	ENDWHERE
	WHERE(a2dUm.lt.0)
		a2dUm=0
	ENDWHERE
ENDWHERE


return
end subroutine


