!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
! Delete Input Meteorological Map:
!		- Rain 
!		- Temperature 
!		- Relative Umidity 
!		- Wind
!		- Short Wave Radiation
!***********************************************************************************  
 
subroutine DeleteMeteoMapBinary(xstring,iRows,iCols)

implicit none
include 'DeclarationH.f90'

CHARACTER*12  xstring,ext
CHARACTER*500  file,rainfile,temperaturefile,radiationfile,windfile,umidityfile,maskfile,laifile,line,laifilegz
CHARACTER*200 sCheck
INTEGER i,j,ios,is,iTLength,iiw,iw
INTEGER iRows,iCols,iFhum
INTEGER*4 iLStr,iPathLenght,temp(iRows,iCols)
REAL*8 a2dMaskSnow(iRows,iCols),dTW
INTEGER anno,mese,giorno,ora,flag, min,iFgzip

if(dMeteoGroup.eq.3)THEN !Grouping aaaamm
	iPathLenght = iLStr(sPathRain)
	rainfile	= sPathRain(1:iPathLenght)//xstring(1:4)//sBar//xstring(5:6)//sBar//'Rain_'//xstring

	iPathLenght = iLStr(sPathRain)
	maskfile	= sPathRain(1:iPathLenght)//xstring(1:4)//sBar//xstring(5:6)//sBar//'SnowMask_'//xstring

	iPathLenght = iLStr(sPathTemp)
	temperaturefile =sPathTemp(1:iPathLenght)//xstring(1:4)//sBar//xstring(5:6)//sBar//'Temperature_'//xstring

	iPathLenght = iLStr(sPathRad)
	radiationfile	=sPathRad(1:iPathLenght)//xstring(1:4)//sBar//xstring(5:6)//sBar//'Radiation_'//xstring
	!radiationfile	=sPathRad(1:iPathLenght)//'Rad_'//xstring

	iPathLenght = iLStr(sPathWind)
	windfile		=sPathWind(1:iPathLenght)//xstring(1:4)//sBar//xstring(5:6)//sBar//'Wind_'//xstring

	iPathLenght = iLStr(sPathUmid)
	umidityfile		=sPathUmid(1:iPathLenght)//xstring(1:4)//sBar//xstring(5:6)//sBar//'RelUmid_'//xstring
endif

if(dMeteoGroup.eq.2)THEN !Grouping aaaa
	iPathLenght = iLStr(sPathRain)
	rainfile	= sPathRain(1:iPathLenght)//xstring(1:4)//sBar//'Rain_'//xstring

	iPathLenght = iLStr(sPathRain)
	maskfile	= sPathRain(1:iPathLenght)//xstring(1:4)//sBar//'SnowMask_'//xstring

	iPathLenght = iLStr(sPathTemp)
	temperaturefile =sPathTemp(1:iPathLenght)//xstring(1:4)//sBar//'Temperature_'//xstring

	iPathLenght = iLStr(sPathRad)
	radiationfile	=sPathRad(1:iPathLenght)//xstring(1:4)//sBar//'Radiation_'//xstring
	!radiationfile	=sPathRad(1:iPathLenght)//'Rad_'//xstring

	iPathLenght = iLStr(sPathWind)
	windfile		=sPathWind(1:iPathLenght)//xstring(1:4)//sBar//'Wind_'//xstring

	iPathLenght = iLStr(sPathUmid)
	umidityfile		=sPathUmid(1:iPathLenght)//xstring(1:4)//sBar//'RelUmid_'//xstring
endif

if(dMeteoGroup.eq.1)THEN !Grouping one dir
	iPathLenght = iLStr(sPathRain)
	rainfile	= sPathRain(1:iPathLenght)//'Rain_'//xstring

	iPathLenght = iLStr(sPathRain)
	maskfile	= sPathRain(1:iPathLenght)//'SnowMask_'//xstring

	iPathLenght = iLStr(sPathTemp)
	temperaturefile =sPathTemp(1:iPathLenght)//'Temperature_'//xstring

	iPathLenght = iLStr(sPathRad)
	radiationfile	=sPathRad(1:iPathLenght)//'Radiation_'//xstring
	
	iPathLenght = iLStr(sPathWind)
	windfile		=sPathWind(1:iPathLenght)//'Wind_'//xstring

	iPathLenght = iLStr(sPathUmid)
	umidityfile		=sPathUmid(1:iPathLenght)//'RelUmid_'//xstring
endif
!-----------------------------------------------------------------------------------
!   Opening Input Compressed (for Windows) File
!-----------------------------------------------------------------------------------

!call system('7z e -tgzip '//rainfile//'.gz > uscita.txt')
!open(21,file=rainfile,status='old')

temp=-8888

!Lettura Piogge
iTLength=iLStr(rainfile)
sCheck=''
sCheck(1:iTLength+3)=rainfile(1:iTLength)//'.gz'
dTW=dMinWait/dWin !Number of attempts waiting for meteo dat
iiw=1
528 open(22,file=sCheck,status='old',err=228)
close(22)
if(iLinux.eq.10)then
!****Linux****
	line='rm -r  '//rainfile(1:iTLength)//'.gz'
else
	line='erase '//rainfile(1:iTLength)//'.gz'
endif
call system(line)

!Lettura Temperature
iTLength=iLStr(temperaturefile)
sCheck=''
sCheck(1:iTLength+3)=temperaturefile(1:iTLength)//'.gz'
open(22,file=sCheck,status='old',err=228)
close(22)
if(iLinux.eq.10)then
	line='rm -r '//temperaturefile(1:iTLength)//'.gz'
else
	line='erase '//temperaturefile(1:iTLength)//'.gz'	
endif
call system(line)

!Lettura Radiazioni
iTLength=iLStr(radiationfile)
sCheck=''
sCheck(1:iTLength+3)=radiationfile(1:iTLength)//'.gz'
open(22,file=sCheck,status='old',err=228)
close(22)
!line='C:\Programmi\7-Zip\7z.exe e -tgzip '//radiationfile(1:iTLength)//'.gz > lozip.xt'
if(iLinux.eq.10)then
!****Linux****
	line='rm -r '//radiationfile(1:iTLength)//'.gz'
else
	line='erase '//radiationfile(1:iTLength)//'.gz'
endif
call system(line)

!Lettura Vento
iTLength=iLStr(windfile)
sCheck=''
sCheck(1:iTLength+3)=windfile(1:iTLength)//'.gz'
open(22,file=sCheck,status='old',err=228)
close(22)

!line='C:\Programmi\7-Zip\7z.exe e -tgzip '//windfile(1:iTLength)//'.gz > lozip.xt'
if(iLinux.eq.10)then
!****Linux****
	line='rm -r '//windfile(1:iTLength)//'.gz'
else
	line='erase '//windfile(1:iTLength)//'.gz'
endif
call system(line)

!Lettura Umidità
iTLength=iLStr(umidityfile)
iFhum=1 !Se 1 nome file RelUmid
open(22,file=umidityfile(1:iTLength)//'.gz',status='old',form='unformatted', &
          access='direct',recl=iRows*4,iostat=ios,err=225)
close(22)
if(1.eq.2)then !Prova con altro nome
225	iPathLenght = iLStr(sPathUmid)
	umidityfile		=sPathUmid(1:iPathLenght)//'Humidity_'//xstring
    if(dMeteoGroup.eq.2)then
		umidityfile		=sPathUmid(1:iPathLenght)//xstring(1:4)//sBar//'Humidity_'//xstring
	elseif(dMeteoGroup.eq.3)then
		umidityfile		=sPathUmid(1:iPathLenght)//xstring(1:4)//sBar//xstring(5:6)//sBar//'Humidity_'//xstring
	endif

	iTLength=iLStr(umidityfile)
	sCheck=''
	sCheck(1:iTLength+3)=umidityfile(1:iTLength)//'.gz'
	open(22,file=sCheck,status='old',err=228)
	close(22)
    	
	iFhum=0
endif

!line='C:\Programmi\7-Zip\7z.exe e -tgzip '//umidityfile(1:iTLength)//'.gz > lozip.xt'
if(iLinux.eq.10)then
!****Linux****
	line='rm -r '//umidityfile(1:iTLength)//'.gz'
else
	line='erase '//umidityfile(1:iTLength)//'.gz'
endif
call system(line)

if(1.eq.2)then
228 write(*,*)'File ',sCheck,' not found'
	write(*,*)'Program will be terminated '
	stop
endif

return
end subroutine


