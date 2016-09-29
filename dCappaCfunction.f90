!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

REAL*8 FUNCTION dCappaCfunction(dH,dDth,dIf)

implicit none
include 'DeclarationH.f90' ! Dichiarazioni per variabili comuni
REAL*8 dCappaCnew,dH,dDth,dCappaMa,dIf

if(dH.lt.0.0)dH=0.000001
dCappaCfunction=dCappaV*0.0+0.1+dCappaC*(dIf**0.5)*dH**dBc

dCappaMa=3600/dDth*0.7

if(dCappaCfunction.gt.dCappaMa)then
	dCappaCfunction=dCappaMa
endif

END 