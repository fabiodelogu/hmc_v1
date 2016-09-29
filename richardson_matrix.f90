!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!*****************************************************************************************
	subroutine richardson_matrix(iRows,iCols,Ta,Patm,Rb)

	implicit none	
	include 'DeclarationH.f90'

!-----------------------------------------------------------------------------------------
!	Dichiarazioni delle variabili
!   Per il calcolo del numero di Richardson (bulk) le uscite saranno:
!	Pt temperatura potenziale 
!	Rb numero di richardson (bulk)
!-----------------------------------------------------------------------------------------
	
	integer iRows,iCols
	real*8 Ta(iRows,iCols),Patm(iRows,iCols)
	real*8 PT(iRows,iCols),PT0(iRows,iCols),Rb(iRows,iCols)
	real*8 zrif,ga,Rd,cp
	parameter (zrif=3.0, cp=1004.0, ga=9.81, Rd=287.0)	! Altezza di riferimento per il vento [m]	
														! Calore specifico a pressione costante [J/kg/a2dK]
														! Accelerazione di gravit� [m/a2dS^2]	
														! Rd costante dei gas per l'aria [J/kg a2dK]
Rb=-0.9 !Inizializzazione
WHERE (a2dCTime.gt.0.0.and.a2dW.gt.0.0)
      PT=Ta*(1000.0/Patm)**(Rd/cp)
      PT0=a2dTS*(1000.0/Patm)**(Rd/cp)
      Rb=(ga/PT)*(PT-PT0)*zrif/(a2dW**2) 
ELSEWHERE(a2dCTime.gt.0.0)
     PT=Ta*(1000.0/Patm)**(Rd/cp)
     PT0=a2dTS*(1000.0/Patm)**(Rd/cp)
     Rb=(ga/PT)*(PT-PT0)*zrif/(0.1**2)
ENDWHERE

! qua nel programma di francesca Rb (Ric) vengono numeri compresi fra 0 e -1 verificare un p�!

	return
	end subroutine

