!***********************************************************************************
!	Subroutine read_land_data
!	Legge il file sBasin.info.txt
!***********************************************************************************  
 
subroutine WriteStateMatrix(iRows,iCols,xstring)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,j
character*50 file
INTEGER*4 iLStr,iPathLenght
CHARACTER*12  xstring
REAL*8 tmp(iRows,iCols)

iPathLenght = iLStr(sPathLandData)
700 FORMAT (1000(f9.2,1x))

!-----------------------------------------------------------------------------------
!   Writing matrix V
!-----------------------------------------------------------------------------------

open(unit=1,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//xstring//'.V.txt',status='unknown')
do i=1,iRows,1
	write(1,700) (a2dV(i,j),j=1,iCols) 
enddo
close(1)
!-----------------------------------------------------------------------------------
!   Writing matrix Retention
!-----------------------------------------------------------------------------------

open(unit=1,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//xstring//'.Retention.txt',status='unknown')
do i=1,iRows,1
	write(1,700) (a2dRetention(i,j),j=1,iCols) 
enddo
close(1)
!-----------------------------------------------------------------------------------
!   Writing matrix Hydro Level
!-----------------------------------------------------------------------------------

open(unit=1,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//xstring//'.Hydro.txt',status='unknown')
do i=1,iRows,1
	write(1,700) (a2dHydro(i,j),j=1,iCols) 
enddo
close(1)

!-----------------------------------------------------------------------------------
!   Writing matrix Routing
!-----------------------------------------------------------------------------------
701 FORMAT (1000(f9.4,x))
open(unit=1,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//xstring//'.Routing.txt',status='unknown')
do i=1,iRows,1
	write(1,701) (a2dRouting(i,j),j=1,iCols) 
enddo
close(1)

!-----------------------------------------------------------------------------------
!   Writing matrix Water Table Level
!-----------------------------------------------------------------------------------
tmp=0
WHERE(a2dDem.gt.0.0)
	tmp=(a2dDem-a2dVwt)*1000
ENDWHERE
702 FORMAT (1000(f11.5,1x))
open(unit=1,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//xstring//'.Vw.txt',status='unknown')
do i=1,iRows,1
	write(1,702) (tmp(i,j),j=1,iCols) 
enddo
close(1)

return
end subroutine