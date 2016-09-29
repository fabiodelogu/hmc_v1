!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine read_land_data
!	Legge il file sBasin.info.txt
!***********************************************************************************  
 
subroutine SaltaMeteoData(iRows,iCols,a1dLatP,a1dLonP,a1dZP,a1dLatT,a1dLonT,a1dZT,a1dLatK,a1dLonK,a1dZK,a1dLatW, &
			a1dLonW,a1dZW,a1dLatUm,a1dLonUm,a1dZUm)

implicit none
include 'DeclarationH.f90'

INTEGER i,iLapseChoice,iRows,iCols,iNStazS
REAL*8 dMax,dMin,matrice_var(iRows,iCols)

REAL*8 a1dLatP(iNStazP),a1dLonP(iNStazP),a1dZP(iNStazP)
REAL*8 a1dLatT(iNStazT),a1dLonT(iNStazT),a1dZT(iNStazT)
REAL*8 a1dLatK(iNStazK),a1dLonK(iNStazK),a1dZK(iNStazK)
REAL*8 a1dLatW(iNStazW),a1dLonW(iNStazW),a1dZW(iNStazW)
REAL*8 a1dLatUm(iNStazUm),a1dLonUm(iNStazUm),a1dZUm(iNStazUm)

REAL*8 a1dP1(iNStazP),a1dT1(iNStazT),a1dK1(iNStazK),a1dW1(iNStazW),a1dUm1(iNStazUm)
REAL*8, ALLOCATABLE :: a1dLatS(:),a1dLonS(:),a1dZS(:)
!-----------------------------------------------------------------------------------------  
!   Lettura file di input meteo per Interpolazione
!-----------------------------------------------------------------------------------------  
READ(21,*)(a1dP1(i),i=1,iNStazP)	!RAINFALL IN mm
READ(22,*)(a1dT1(i),i=1,iNStazT)	!TEMPERATURE IN °C
READ(23,*)(a1dK1(i),i=1,iNStazK) 	!LA RADIAZIONE E' in a2dW/m^2
READ(24,*)(a1dW1(i),i=1,iNStazW)	!WIND SPEED IN m a2dS^-1
READ(25,*)(a1dUm1(i),i=1,iNStazUm)	!RELATIVE HUMIDITY % nei file di input l'umidità è già divisa per 100!	 
a1dUm1=a1dUm1/100.0

!-----------------------------------------------------------------------------------------  
!   Interpolazione dati micrometeorologici (chiamata alla subroutine di interpolazione dei
!	dati ed assegnazione dei valori massimo e minimo all'interno dei quali fare 
!	l'interpolazione dMin e dMax e scelta dell'interpolazione da fare o come regressione 
!	sulle quote (iLapseChoice=1), come adiabatico (iLapseChoice=2), solo IDW (iLapseChoice=3))
!-----------------------------------------------------------------------------------------  


return
end subroutine


