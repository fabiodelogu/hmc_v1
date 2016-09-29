!***********************************************************************************
! Read Input Meteorological Map:
!		- Rain 
!
!***********************************************************************************  
 
subroutine ReadRainMapBinary(xstring,iRows,iCols)

implicit none
include 'DeclarationH.f90'

CHARACTER*12  xstring,ext
CHARACTER*500  file,rainfile,sRainout,maskfile,laifile,line
INTEGER i,j,ios,is,iTLength
INTEGER iRows,iCols
INTEGER*4 iLStr,iPathLenght,temp(iRows,iCols)
REAL*8 a2dMaskSnow(iRows,iCols)
INTEGER anno,mese,giorno,ora,flag, min

iPathLenght = iLStr(sPathRain)
rainfile	= sPathRain(1:iPathLenght)//'Rain_'//xstring
sRainout = 'Rain_'//xstring

iPathLenght = iLStr(sPathRain)
maskfile	= sPathRain(1:iPathLenght)//'SnowMask_'//xstring

!-----------------------------------------------------------------------------------
!   Opening Input Compressed (for Windows) File
!-----------------------------------------------------------------------------------

!call system('7z e -tgzip '//rainfile//'.gz > uscita.txt')
!open(21,file=rainfile,status='old')

temp=-8888

!Lettura Piogge
iTLength=iLStr(rainfile)
call ReadMapBinary(iRows,iCols,rainfile,sRainout,a2dRain,10,iLinux)

WHERE (a2dDem.LE.-10.0)
	a2dRain=-9999
ENDWHERE
temp=0

 !leggo la maschera della neve

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
	a2dMaskSnow=real(temp)/10
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



return
end subroutine


