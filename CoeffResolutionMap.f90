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
 
subroutine CoeffResolutionMap(iRows,iCols)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,j
character*50 file,sVar
INTEGER*4 iLStr,iPathLenght
REAL*8 dNCellMin,dACellMin,dCoffMean
REAL*8 dMaxW,dMinW
REAL*8, ALLOCATABLE :: a2dArea(:,:),a2dSecWidth(:,:)

ALLOCATE  (a2dArea(iRows,iCols),a2dSecWidth(iRows,iCols))
iPathLenght = iLStr(sPathLandData)

a2dArea=0.0
a2dSecWidth=-9999
!-----------------------------------------------------------------------------------
!   Opening, Reading and Closing a2dDem file (unit=2) 
!-----------------------------------------------------------------------------------
a2dCoeffResol=0.0
dRateResol=0.5 !
open(unit=2,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.area.txt',status='old',err=514)
do i=1,6
	read(2,*) 
enddo
do i=1,iRows,1
	read(2,*) (a2dArea(iRows-i+1,j),j=1,iCols)
enddo
write(*,*)'*****************'
write(*,*)'Read file of Aree'
write(*,*)'*****************'

close(2)

!Coeff to regulate subsurface and deep flow
dNCellMin=MINVAL(MINVAL(a2dArea,dim=1,mask=a2iMask.gt.0.0.and.a2dArea.gt.0.0))
dACellMin=MINVAL(MINVAL(a2dAreaCell,dim=1,mask=a2iMask.gt.0.0.and.a2dAreaCell.gt.0.0))
write(*,*)dNCellMin,dACellMin
IF(dNCellMin.gt.1)THEN !Area matrix and not cells number matrix
	write(*,*)'Conversion from area to cell number'
	WHERE(a2dDem.gt.0.and.a2iChoice.gt.0)
		a2dArea=a2dArea/a2dAreaCell*1000000 !Approximate formula
	ENDWHERE
ENDIF


!Section width
WHERE(a2dDem.gt.0.and.a2iChoice.gt.0)
	a2dSecWidth=0.01*(a2dArea*a2dAreaCell/1000000)**0.4*1000 !width in m
ENDWHERE

dMaxW=MAXVAL(MAXVAL(a2dSecWidth,dim=1,mask=a2iChoice.gt.0.0))
dMinW=MINVAL(MINVAL(a2dSecWidth,dim=1,mask=a2iChoice.eq.1.and.a2dDem.gt.0))

write(*,*),'Max and Min Width: ',dMaxW,dMinW

if(1.eq.0)then
	514 continue	
endif

a2dCoeffResol=1.0
WHERE(a2dDem.gt.0.and.a2iChoice.ge.0)
	a2dCoeffResol=dexp(-dsqrt(a2dAreaCell)*0.0007)
ENDWHERE
WHERE(a2dCoeffResol.lt.0)
	a2dCoeffResol=0.0
ENDWHERE

write(*,*),'Coeff Res Mean: ',a2dCoeffResol((iRows/2),(iCols/2))

! Eliminate the border values
a2dCoeffResol(1,:)=-9999
a2dCoeffResol(:,1)=-9999
a2dCoeffResol(:,iCols)=-9999
a2dCoeffResol(iRows,:)=-9999

a2dSecWidth(1,:)=-9999
a2dSecWidth(:,1)=-9999
a2dSecWidth(:,iCols)=-9999
a2dSecWidth(iRows,:)=-9999

DEALLOCATE  (a2dArea,a2dSecWidth)


return
end subroutine