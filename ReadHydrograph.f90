!***********************************************************************************
!	Reding Meteorological Data series
!		- Rain
!		- Temperature
!		- Relative Umidity 
!		- Wind
!		- Short Wave Radiation
!***********************************************************************************  
 
subroutine ReadHydrograph(dTstart,d)

implicit none
include 'DeclarationH.f90'

character*500 rainfile,temperaturefile,radiationfile,windfile,umidityfile
character*50 file

INTEGER i,j,nb
INTEGER*4 iLStr,iNameLength
REAL*8 t,dTstart,d

iNameLength = iLStr(sBasin)
i=1
OPEN(210,file=sBasin(1:iNameLenght)//'IdrogrammiHyperCdrift.out',status='old',form='formatted',access='sequential')

604 READ(210,*,end=603) t,(a2dQsections(nb,i),nb=1,iNumBasins)
i=i+1
GO TO 604
603 CLOSE(210)

dTstart=nint(t*3600/d)

return
end subroutine


