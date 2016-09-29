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
 
subroutine ReadSectionsDam(iRows,iCols)

implicit none
include 'DeclarationH.f90'

integer iSS,j,iRows,iCols
character*50 file
INTEGER*4 iLStr,iPathLenght


iPathLenght = iLStr(sPathLandData)

!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing Porosity file (unit=6) 
!-----------------------------------------------------------------------------------


open(unit=8,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.sections.txt',status='old')

do iSS=1,iNumBasins
	read(8,*)a2dXYsections(iSS,2),a2dXYsections(iSS,1)
	a2dXYsections(iSS,2)=iRows-a2dXYsections(iSS,2)+1 !!Vale per VDA
end do

close(8)
return
end subroutine