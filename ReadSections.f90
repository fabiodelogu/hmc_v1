!***********************************************************************************
!	Subroutine read_land_data
!	Legge il file sBasin.info.txt
!***********************************************************************************  
 
subroutine ReadSections(iRows,iCols)

implicit none
include 'DeclarationH.f90'

integer iS,j,iRows,iCols
character*50 file
INTEGER*4 iLStr,iPathLenght


iPathLenght = iLStr(sPathLandData)

!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing Porosity file (unit=6) 
!-----------------------------------------------------------------------------------


open(unit=8,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.sections.txt',status='old')

do iS=1,iNumBasins
	read(8,*)a2dXYsections(iS,2),a2dXYsections(iS,1)
	a2dXYsections(iS,2)=iRows-a2dXYsections(iS,2)+1 !!Vale per VDA
end do

close(8)
return
end subroutine