!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
! READ OUTPUT Meteorological Map:
! INPUT
! xstring: data aaaammggHHMM
! iRows, iCols: numero righe e colonne
! sVariable: nome della variabile, usata per ricostruire il nome del file da leggere
! sPathDati: percorso mappe input
! a2dMap: mappa da scrivere su file
!***********************************************************************************  
 
subroutine ReadMapBinary(iRows,iCols,sPathDati,sNameOut,a2dMap,iDivide,iTgzip)

implicit none

CHARACTER*12  xstring
CHARACTER*500  file,sFileName,line,sVariable
CHARACTER*500 sPathDati,sNameOut
INTEGER i,j,ios,is,iLength,iLo
INTEGER iTgzip !1 se file salvato con tar gzip
INTEGER iRows,iCols
INTEGER*4 iLStr,iPathLenght,iPathLenght2
REAL*8 a2dMap(iRows,iCols)
INTEGER*4 a2iMap(iRows,iCols),iDivide
!INTEGER*4 temp(iRows,iCols)

iPathLenght = iLStr(sPathDati)
sFileName	= sPathDati(1:iPathLenght)
iLo = iLStr(sNameOut)

a2dMap=0.0
!-----------------------------------------------------------------------------------
!   Opening Input Compressed (for Windows) File
!-----------------------------------------------------------------------------------

open(22,file=sFileName(1:iPathLenght)//'.gz',form='unformatted', &
          access='direct',recl=iRows*4,iostat=ios,status='old',err=324)
close(22)
if(iTgzip.eq.1)then
	line='7z.exe e -tgzip '//sFileName(1:iPathLenght)//'.gz > lozip.xt'
elseif(iTgzip.eq.0)then
	line='7z.exe e '//sFileName(1:iPathLenght)//'.gz > lozip.xt'
elseif(iTgzip.eq.10)then
!****Linux****
	line='gunzip -c '//sFileName(1:iPathLenght)//'.gz > '//sNameOut(1:iLo)
end if


call system(line)

iLength=iLStr(sFileName)

open(22,file=sNameOut(1:iLo),form='unformatted', &
          access='direct',recl=iRows*4,iostat=ios)

!temp=int(a2dMap)
a2iMap=0.0
do j=1,iCols
   read(22,rec=j)(a2iMap(iRows-i+1,j),i=1,iRows)
end do
close(22)

a2dMap=a2iMap/real(iDivide)

!Linux
!line='gzip -f '//sFileName(1:iLength)//' > LogZipWriteMap.txt'
!Windows
!line='C:\Programmi\7-Zip\7z.exe a -tzip '//sFileName(1:iLength)//'.gz  '//sFileName(1:iLength)//'> LogZipWriteMap.txt'
!line='7z.exe a -tzip '//sFileName(1:iLength)//'.gz  '//sFileName(1:iLength)//'> LogZipWriteMap.txt'

if(iTgzip.eq.10)then
	line='rm -r '//sNameOut(1:iLo)
else
	line='erase '//sNameOut(1:iLo)
endif
call system(line)

if(1.eq.2)then
	324 write(*,*)'File ',sFileName(1:iPathLenght)//'.gz not found'
	write(*,*)'Program will be terminated '
	stop
endif
return
end subroutine


