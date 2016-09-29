!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
! Read Input Meteorological Map:
!		- Rain 
!		- Temperature 
!		- Relative Umidity 
!		- Wind
!		- Short Wave Radiation
!***********************************************************************************  
 
subroutine ReadMeteoMapBinary(xstring,iRows,iCols)

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
528 open(22,file=sCheck,status='old',err=428)
close(22)
GO TO 328
428 do iw=iiw,dWin
		write(*,*)'Wait...',iw,' of ',int(dWin)
		CALL SLEEP (int(dTW)) 
		iiw=iiw+1
		GO TO 528
	enddo
GO TO 228

328 if(iLinux.eq.10)then
!****Linux****
	line='gunzip -c '//rainfile(1:iTLength)//'.gz > Rain_'//xstring
else
	line='7z.exe e -tgzip '//rainfile(1:iTLength)//'.gz > lozip.xt'
endif
call system(line)
open(22,file='Rain_'//xstring,status='old',form='unformatted', &
          access='direct',recl=iRows*4,iostat=ios)
do j=1,iCols
   read(22,rec=j)(temp(iRows-i+1,j),i=1,iRows)
end do

a2dRain=real(temp)/dRescFct
temp=0
close(22)
if(iLinux.eq.10)then
	call system('rm -r Rain_*')
else
	call system('erase Rain_*')	
endif
WHERE (a2dDem.LE.-10.0)
	a2dRain=-9999
ENDWHERE
temp=0
 !leggo anche la maschera della neve

if(1.eq.1)then
	iTLength=iLStr(maskfile)
	a2dMaskSnow=0.0
	open(22,file=maskfile(1:iTLength)//'.gz',status='old',form='unformatted', &
			  access='direct',recl=iRows*4,iostat=ios,err=224)
	close(22)

	if(iLinux.eq.10)then
	!*** Linux ****
		line='gunzip -c '//maskfile(1:iTLength)//'.gz > SnowMask_'//xstring
	else
		line='7z.exe e -tgzip '//maskfile(1:iTLength)//'.gz > lozip.xt'
	endif


	call system(line)
	open(22,file='SnowMask_'//xstring,status='old',form='unformatted', &
			  access='direct',recl=iRows*4,iostat=ios,err=224)
	do j=1,iCols
	   read(22,rec=j)(temp(iRows-i+1,j),i=1,iRows)
	end do
	close(22)
	a2dMaskSnow=real(temp)/dRescFct
	WHERE(a2dMaskSnow.eq.1)
		a2dEvapot=-1
	ENDWHERE
	if(iLinux.eq.10)then
		call system('rm -r SnowMask_*')	
	else
		call system('erase SnowMask_*')	
	endif
endif
224 continue

!Lettura Temperature
iTLength=iLStr(temperaturefile)
sCheck=''
sCheck(1:iTLength+3)=temperaturefile(1:iTLength)//'.gz'
open(22,file=sCheck,status='old',err=228)
close(22)
if(iLinux.eq.10)then
	line='gunzip -c '//temperaturefile(1:iTLength)//'.gz > Temperature_'//xstring
else
	line='7z.exe e -tgzip '//temperaturefile(1:iTLength)//'.gz > lozip.xt'	
endif
call system(line)
open(22,file='Temperature_'//xstring,status='old',form='unformatted', &
          access='direct',recl=iRows*4,iostat=ios)
do j=1,iCols
   read(22,rec=j)(temp(iRows-i+1,j),i=1,iRows)
end do
close(22)
a2dTemp=real(temp)/dRescFct
WHERE (a2dDem.LE.0.0)
	a2dTemp=-9999
ENDWHERE

temp=0
if(iLinux.eq.10)then
	call system('rm -r Temperature_*')
else
	call system('erase Temperature_*')
endif

!Lettura Radiazioni
iTLength=iLStr(radiationfile)
sCheck=''
sCheck(1:iTLength+3)=radiationfile(1:iTLength)//'.gz'
open(22,file=sCheck,status='old',err=228)
close(22)
!line='C:\Programmi\7-Zip\7z.exe e -tgzip '//radiationfile(1:iTLength)//'.gz > lozip.xt'
if(iLinux.eq.10)then
!****Linux****
	line='gunzip -c '//radiationfile(1:iTLength)//'.gz > Radiation_'//xstring
else
	line='7z.exe e -tgzip '//radiationfile(1:iTLength)//'.gz > lozip.xt'
endif
call system(line)
!Orba
open(22,file='Radiation_'//xstring,status='old',form='unformatted', &
!open(22,file='Rad_'//xstring,status='old',form='unformatted', &
          access='direct',recl=iRows*4,iostat=ios,err=222)
do j=1,iCols
   read(22,rec=j)(temp(iRows-i+1,j),i=1,iRows)
end do
close(22)
a2dK=real(temp)/dRescFct
WHERE (a2dDem.LE.0.0)
	a2dK=-9999
ENDWHERE
temp=0
if(1.eq.2)then !Se non trovo la radiazione metto 0
222 a2dK=0.0
endif
if(iLinux.eq.10)then
	call system('rm -r Radiation_*')
else
	call system('erase Radiation_*')
endif
!call system('erase Rad_*')

!Lettura Vento
iTLength=iLStr(windfile)
sCheck=''
sCheck(1:iTLength+3)=windfile(1:iTLength)//'.gz'
open(22,file=sCheck,status='old',err=228)
close(22)

!line='C:\Programmi\7-Zip\7z.exe e -tgzip '//windfile(1:iTLength)//'.gz > lozip.xt'
if(iLinux.eq.10)then
!****Linux****
	line='gunzip -c '//windfile(1:iTLength)//'.gz > Wind_'//xstring
else
	line='7z.exe e -tgzip '//windfile(1:iTLength)//'.gz > lozip.xt'
endif
call system(line)
open(22,file='Wind_'//xstring,status='old',form='unformatted', &
          access='direct',recl=iRows*4,iostat=ios)
do j=1,iCols
   read(22,rec=j)(temp(iRows-i+1,j),i=1,iRows)
end do
close(22)
a2dW=real(temp)/dRescFct
WHERE (a2dDem.LE.0.0)
	a2dW=-9999
ENDWHERE
temp=0
if(iLinux.eq.10)then
	call system('rm -r Wind_*')
else
	call system('erase Wind_*')	
endif

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
	line='gunzip -c '//umidityfile(1:iTLength)//'.gz > RelUmid_'//xstring
	iFhum=1 !In Linux apro il file decompresso con lo stesso nome RelUmid
else
	line='7z.exe e -tgzip '//umidityfile(1:iTLength)//'.gz > lozip.xt'
endif
call system(line)
if(iFhum.eq.1)then
	open(22,file='RelUmid_'//xstring,status='old',form='unformatted', &
			  access='direct',recl=iRows*4,iostat=ios)
else

	open(22,file='Humidity_'//xstring,status='old',form='unformatted', &
			  access='direct',recl=iRows*4,iostat=ios)
endif

do j=1,iCols
   read(22,rec=j)(temp(iRows-i+1,j),i=1,iRows)
end do
close(22)
a2dUm=real(temp)/(dRescFct*100)
!a2dUm=real(temp)/10
WHERE (a2dDem.LE.0.0)
	a2dUm=-9999
ENDWHERE
temp=0
if(iLinux.eq.10)then
	call system('rm -r RelUmid_*')
else
	call system('erase RelUmid_*')
	call system('erase Humidity_*')
endif
!-----------------------------------------------------------------------------------
!	Read LAI
!-----------------------------------------------------------------------------------

READ( xstring(1:4), '(i4)' )  anno
READ( xstring(5:6), '(i2)' )  mese
READ( xstring(7:8), '(i2)' )  giorno
READ( xstring(9:10), '(i2)' ) ora
READ( xstring(11:12), '(i2)') min

if(ora.eq.0)THEN
	if(giorno.eq.1.or.giorno.eq.15)THEN
		iPathLenght = iLStr(sPathLai)
		laifile		=sPathLai(1:iPathLenght)//'LAI_'//xstring(1:8)//'0000'
		iTLength=iLStr(laifile)
		iFgzip=1 !If 1 the file is not gzipped
		!Check if life file is gzipped
		open(22,file=laifile(1:iTLength)//'.gz',status='old',form='unformatted', &
				  access='direct',recl=iRows*4,iostat=ios,err=233)
		close(22)

		if(iLinux.eq.10)then
		!****Linux****
			line='gunzip -c '//laifile(1:iTLength)//'.gz > LAI_'//xstring(1:8)//'0000'
		else
			line='7z.exe e -tgzip '//laifile(1:iTLength)//'.gz > lozip.xt'
		endif
		call system(line)
		iFgzip=0

233		open(22,file=laifile(1:iTLength),status='old',form='unformatted', &
				  access='direct',recl=iRows*4,iostat=ios,err=223)
		close(22)

		!line='7z.exe e -tgzip '//laifile(1:iTLength)//'.gz > lozip.xt'
	!	call system(line)
		!open(22,file='LAI_'//xstring(1:8)//'0000',status='old',form='unformatted', &
		open(22,file=laifile(1:iTLength),status='old',form='unformatted', &
				  access='direct',recl=iRows*4,iostat=ios)
		do j=1,iCols
		   read(22,rec=j)(temp(iRows-i+1,j),i=1,iRows)
		end do
		close(22)
		
		a2dLai=real(temp)/dRescFct
		dFlagLai=1
		write(*,*)'Uso il LAI '//xstring(1:8)
		if(iFgzip.eq.0)then !If the file is gzipped, delete the unzipped version
			if(iLinux.eq.10)then
				call system('rm -r LAI_*')
			else
				call system('erase LAI_*')				
			endif
		endif
		if(1.eq.2)then !No LAI file available
223			dFlagLai=0	
			write(*,*)'Non uso il LAI'
		endif
		

	endif
endif



if(1.eq.2)then
228 write(*,*)'File ',sCheck,' not found'
	write(*,*)'Program will be terminated '
	stop
endif

!-----------------------------------------------------------------------------------
! Delete txt input files
!-----------------------------------------------------------------------------------
!call system('del '//rainfile//'.txt') !delate txt file
!call system('del '//temperaturefile//'.txt')
!call system('del '//radiationfile//'.txt')
!call system('del '//windfile//'.txt')
!call system('del '//umidityfile//'.txt')

return
end subroutine


