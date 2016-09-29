
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
 
subroutine WriteEsriFile(iRows,iCols,sVar,tmp,iWtype)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols,iWtype
integer i,j
character*50 file,sVar
INTEGER*4 iLStr,iPathLenght,iVarLenght
CHARACTER*12  xstring,a
REAL*8 tmp(iRows,iCols)

iPathLenght = iLStr(sPathLandData)
iVarLenght = iLStr(sVar)
!700 FORMAT (1000(f9.2,1x))
700 FORMAT (1000(f9.2,1x))
701 FORMAT (1000(f12.6,1x))
!-----------------------------------------------------------------------------------
!   Writing matrix V
!-----------------------------------------------------------------------------------

open(unit=1,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.'//sVar(1:iVarLenght)//'.txt',status='unknown')

write(1,*)'ncols	',iCols
write(1,*)'nrows	',iRows
write(1,*)'xllcorner	',dXDemLon
write(1,*)'yllcorner	',dXDemLat
write(1,*)'cellsize  ',dDemPasLat
write(1,*)'NODATA_value  -9999'
do i=1,iRows,1
	if(iWtype.eq.1)then
		write(1,700) (tmp(iRows-i+1,j),j=1,iCols)
	endif
	if(iWtype.eq.2)then
		write(1,701) (tmp(iRows-i+1,j),j=1,iCols)
	endif
enddo
close(1)

return
end subroutine