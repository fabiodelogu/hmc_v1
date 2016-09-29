!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Reding Meteorological Data series
!		- Rain
!		- Temperature
!		- Relative Umidity 
!		- Wind
!		- Short Wave Radiation
!***********************************************************************************  
 
subroutine ReadMeteoData(a1dLatP,a1dLonP,a1dZP,a1dLatT,a1dLonT,a1dZT,a1dLatK,a1dLonK,a1dZK,a1dLatW,a1dLonW,a1dZW, &
			a1dLatUm,a1dLonUm,a1dZUm)

implicit none
include 'DeclarationH.f90'

character*500 rainfile,temperaturefile,radiationfile,windfile,umidityfile
character*50 file
REAL*8 a1dLatP(iNStazP),a1dLonP(iNStazP),a1dZP(iNStazP)
REAL*8 a1dLatT(iNStazT),a1dLonT(iNStazT),a1dZT(iNStazT)
REAL*8 a1dLatK(iNStazK),a1dLonK(iNStazK),a1dZK(iNStazK)
REAL*8 a1dLatW(iNStazW),a1dLonW(iNStazW),a1dZW(iNStazW)
REAL*8 a1dLatUm(iNStazUm),a1dLonUm(iNStazUm),a1dZUm(iNStazUm)
INTEGER I,J
INTEGER*4 iLStr,iPathLenght


iPathLenght = iLStr(sPathRain)
rainfile	= sPathRain(1:iPathLenght)//'time_series_rain_'//sBasin(1:iNameLenght)//'.txt'

iPathLenght     = iLStr(sPathTemp)
temperaturefile = sPathTemp(1:iPathLenght)//'time_series_temperature_'//sBasin(1:iNameLenght)//'.txt'

iPathLenght =  iLStr(sPathRad)
radiationfile	= sPathRad(1:iPathLenght)//'time_series_radiation_'//sBasin(1:iNameLenght)//'.txt'

iPathLenght = iLStr(sPathWind)
windfile		=sPathWind(1:iPathLenght)//'time_series_wind_'//sBasin(1:iNameLenght)//'.txt'

iPathLenght = iLStr(sPathUmid)
umidityfile		=sPathUmid(1:iPathLenght)//'time_series_humidity_'//sBasin(1:iNameLenght)//'.txt'

!-----------------------------------------------------------------------------------
!   OPENing Input File
!-----------------------------------------------------------------------------------

!	Apertura unità per file dati da stazioni di misura a terra (unit=21..25)     	
open(821,file=rainfile,status='old')
open(822,file=temperaturefile,status='old')
open(823,file=radiationfile,status='old')
open(824,file=windfile,status='old')
open(825,file=umidityfile,status='old')


!-----------------------------------------------------------------------------------
!	Lettura dei file dei dati in ingresso per le stazioni di misura (lon,lat,quota)
!-----------------------------------------------------------------------------------
READ(821,*)(a1dLatP(j),j=1,iNStazP)  
READ(821,*)(a1dLonP(j),j=1,iNStazP) 
READ(821,*)(a1dZP(j),j=1,iNStazP) 

READ(822,*)(a1dLatT(i),i=1,iNStazT)
READ(822,*)(a1dLonT(i),i=1,iNStazT)			
READ(822,*)(a1dZT(i),i=1,iNStazT)

READ(823,*)(a1dLatK(i),i=1,iNStazK)
READ(823,*)(a1dLonK(i),i=1,iNStazK)
READ(823,*)(a1dZK(i),i=1,iNStazK)

READ(824,*)(a1dLatW(i),i=1,iNStazW)		
READ(824,*)(a1dLonW(i),i=1,iNStazW)
READ(824,*)(a1dZW(i),i=1,iNStazW)

READ(825,*)(a1dLatUm(i),i=1,iNStazUm)
READ(825,*)(a1dLonUm(i),i=1,iNStazUm)		
READ(825,*)(a1dZUm(i),i=1,iNStazUm)

return
end subroutine


