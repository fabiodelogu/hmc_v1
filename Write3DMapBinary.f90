!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************


!***********************************************************************************
! WRITE OUTPUT Meteorological Map:
! INPUT
! xstring: data aaaammggHHMM
! iRows, iCols: numero righe e colonne
! sVariable: nome della variabile, usata per ricostruire il nome del file da leggere
! sPathDati: percorso mappe input
! a2dMap: mappa da scrivere su file
!***********************************************************************************  
 
subroutine Write3DMapBinary(xstring,iRows,iCols,iTime,sVariable,sPathDati,a2dMap,iMultiply,iTgzip)

implicit none

CHARACTER*12  xstring
CHARACTER*500  file,sFileName,line
CHARACTER(*) sVariable
CHARACTER*500 sPathDati
INTEGER i,j,k,ios,is,iLength
INTEGER iRows,iCols,iTime
INTEGER iTgzip !1 se file salvato con Windows, 10 se con Linux
INTEGER*4 iLStr,iPathLenght,iPathLenght2
REAL*8 a2dMap(iRows,iCols,iTime)
INTEGER*4 a2iMap(iRows,iCols,iTime),iMultiply
LOGICAL :: file_exist
!INTEGER*4 temp(iRows,iCols)

iPathLenght = iLStr(sPathDati)
iPathLenght2 = iLStr(sVariable)
sFileName	= sPathDati(1:iPathLenght)//sVariable(1:iPathLenght2)//'_'//xstring


!-----------------------------------------------------------------------------------
!   Opening Input Compressed (for Windows) File
!-----------------------------------------------------------------------------------
a2iMap=0.0
a2iMap=nint(a2dMap*real(iMultiply))
!Lunghezza Mappa Output
iLength=iLStr(sFileName)

open(22,file=sFileName(1:iLength),form='unformatted', &
          access='direct',recl=iRows*4,iostat=ios)
inquire(file=sFileName,exist=file_exist)
if(file_exist)then
else
	write(*,*)'Cannot create ',sFileName,' '
	write(*,*)'Program will be terminated '
	stop
endif

!temp=int(a2dMap)
do k=1,iTime
	do j=1,iCols
		write(22,rec=j)(a2iMap(iRows-i+1,j,k),i=1,iRows)
	end do
enddo
close(22)

!****Linux***
!line='gzip -f '//sFileName(1:iLength)//' > LogZipWriteMap.txt'

!Windows
!line='C:\Programmi\7-Zip\7z.exe a -tzip '//sFileName(1:iLength)//'.gz  '//sFileName(1:iLength)//'> LogZipWriteMap.txt'
if(iTgzip.eq.10)then
!****Linux***
	line='gzip -f '//sFileName(1:iLength)//' > LogZipWriteMap.txt'
else
	line='7z.exe a -tzip '//sFileName(1:iLength)//'.gz  '//sFileName(1:iLength)//'> LogZipWriteMap.txt'
endif

call system(line)
if(iTgzip.eq.10)then
	!line='rm -r '//sFileName(1:iLength)
else
	line='erase '//sFileName(1:iLength)
	call system(line)
endif


return
end subroutine

