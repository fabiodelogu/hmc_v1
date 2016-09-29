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
 
subroutine ReadLandData(iRows,iCols)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,j
character*50 file,sVar
INTEGER*4 iLStr,iPathLenght,a2iTmp(iRows,iCols)
REAL*8 dDD,dem_min,aa

dDD=50
iPathLenght = iLStr(sPathLandData)

!-----------------------------------------------------------------------------------
!   Opening, Reading and Closing a2dDem file (unit=2) 
!-----------------------------------------------------------------------------------
open(unit=2,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.dem.txt',status='old')
do i=1,6
	read(2,*) 
enddo
do i=1,iRows,1
	read(2,*) (a2dDem(iRows-i+1,j),j=1,iCols)
enddo
close(2)
!Controllo che il dem non sia <= 0 ed eventualmente lo alzo
dem_min=minval(minval(a2dDem,DIM = 1,MASK=a2dDem.gt.-1000),DIM=1)
write(*,*)'Dem lower value: ',dem_min
IF(dem_min<0)THEN 
	dem_min=abs(dem_min)+0.1
	WHERE(a2dDem.gt.-1000)
		a2dDem=a2dDem+dem_min
	ENDWHERE
ENDIF

!-----------------------------------------------------------------------------------
!   Opening, Reading and Closing Pointers file (unit=3) 
!-----------------------------------------------------------------------------------
open(unit=3,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.pnt.txt',status='old')
do i=1,6
	read(3,*) 
enddo	
do i=1,iRows  
	read(3,*)(a2iPun(iRows-i+1,j),j=1,iCols)
enddo   
close(3)
!-----------------------------------------------------------------
!     FLOW DIRECTIONS
!      
!                                        1 7---------8--------9
!										   !         !        !
!                                          !         !        !
!     0 NO FLOW   1....9 YES               !         !        !
!     I0=(MAT-1)/3-1              Asse J 0 4---------0--------6
!                                          !         !        !
!     J0=MAT-5-3*I0                        !         !        !
!                                          !         !        !
!                                       -1 1---------2--------3
!                                         -1         0        1
!                                                  Asse I
!		PUNTATORI vers. ARCMap
!			   32 64 128
!			   16  0  1
!				8  4  2
!-----------------------------------------------------------------    
! Convert ARCMap to Drift Format
IF((maxval(maxval(a2iPun,DIM = 1),DIM = 1)).gt.10)THEN
	a2iTmp=a2iPun

	WHERE(a2iTmp.eq.32)
		a2iPun=7 !Ok
	ENDWHERE
	WHERE(a2iTmp.eq.64)
		a2iPun=8 !Ok
	ENDWHERE
	WHERE(a2iTmp.eq.128)
		a2iPun=9 !Ok
	ENDWHERE
	WHERE(a2iTmp.eq.1)
		a2iPun=6 !Ok
	ENDWHERE
	WHERE(a2iTmp.eq.2)
		a2iPun=3	!Ok
	ENDWHERE
	WHERE(a2iTmp.eq.4)
		a2iPun=2	!Ok
	ENDWHERE
	WHERE(a2iTmp.eq.8)
		a2iPun=1 !Ok
	ENDWHERE
	WHERE(a2iTmp.eq.16)
		a2iPun=4
	ENDWHERE


	
ENDIF


!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing CN file (unit=4) 
!-----------------------------------------------------------------------------------
open(unit=4,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.cn.txt',status='old')
do i=1,6
	read(4,*) 
enddo	
do i=1,iRows
	read(4,*)(a2dCurNum(iRows+1-i,j),j=1,iCols)
enddo
close(4)

!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing Choice file (unit=5) 
!-----------------------------------------------------------------------------------
open(unit=5,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.choice.txt',status='old')
do i=1,6
	read(5,*) 
enddo	
do i=1,iRows
	read(5,*)(a2iChoice(iRows+1-i,j),j=1,iCols)
enddo
close(5)

open(unit=5,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.rainmul.txt',status='old',err=514)
do i=1,6
	read(5,*) 
enddo	
do i=1,iRows
	read(5,*)(a2dRmul(iRows+1-i,j),j=1,iCols)
enddo
close(5)

514 continue
!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing Cell area file if needed (unit=5) 
!-----------------------------------------------------------------------------------

IF(dCelLat.lt.0.0.and.dCelLon.lt.0.0)THEN
	open(unit=5,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.areacell.txt',status='old',err=512)
	do i=1,6
		read(5,*) 
	enddo	
	do i=1,iRows
		read(5,*)(a2dAreaCell(iRows+1-i,j),j=1,iCols)
	enddo
	close(5)

	!Compute mean cell dimension (meteres)
	dCelLat=SUM(SUM(a2dAreaCell,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	WHERE(a2dDem.gt.0)
		a2iMask=1 !Pongo la Mask uguale al dem perchè non lo ho
	ENDWHERE
	dCelLat=dCelLat/float(sum(sum(a2iMask,dim=1,mask=a2iMask.gt.0.0)))
	dCelLat=sqrt(dCelLat)*1000 !Da Km a M
	dCelLon=dCelLat
	WHERE(a2dDem.gt.0)
		a2dAreaCell=a2dAreaCell*1000000 !da Km2 a m2
	ENDWHERE
	IF(1.eq.2)then
512		write(*,*)'raster file with cells dimensions not found. Exit from the program'
		STOP
	ENDIF

ENDIF

! Eliminate the border values
a2iPun(1,:)=-9999
a2iPun(:,1)=-9999
a2iPun(:,iCols)=-9999
a2iPun(iRows,:)=-9999

a2iChoice(1,:)=-9999
a2iChoice(:,1)=-9999
a2iChoice(:,iCols)=-9999
a2iChoice(iRows,:)=-9999

a2dAreaCell(1,:)=-9999
a2dAreaCell(:,1)=-9999
a2dAreaCell(:,iCols)=-9999
a2dAreaCell(iRows,:)=-9999

a2dDem(1,:)=-9999
a2dDem(:,1)=-9999
a2dDem(:,iCols)=-9999
a2dDem(iRows,:)=-9999

a2dCurNum(1,:)=-9999
a2dCurNum(:,1)=-9999
a2dCurNum(:,iCols)=-9999
a2dCurNum(iRows,:)=-9999


!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing Alpha and Beta file (unit=5) 
!-----------------------------------------------------------------------------------

	open(unit=7,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.alpha.txt',status='old',err=511)
	do i=1,6
		read(7,*) 
	enddo	
	do i=1,iRows
		read(7,*)(a2dAlpha(iRows+1-i,j),j=1,iCols)
	enddo
	close(7)

	open(unit=7,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.beta.txt',status='old',err=511)
	do i=1,6
		read(7,*) 
	enddo	
	do i=1,iRows
		read(7,*)(a2dBeta(iRows+1-i,j),j=1,iCols)
	enddo
	close(7)

	!Se non ho le pendenze beta e alpha le calcolo in base al DEM
    if(1.eq.2)then
511		CALL wt_alpha(iRows,iCols,dDD)
		sVar='alpha'
		CALL WriteEsriFile(iRows,iCols,sVar,a2dAlpha,2)
		sVar='beta'
		CALL WriteEsriFile(iRows,iCols,sVar,a2dBeta,2)
	endif

	!Check a2dAlpha a2dBeta are not 0
    WHERE(a2dDem.gt.0.0.AND.a2dAlpha.le.0.0)
		a2dAlpha=0.00001
	ENDWHERE
	WHERE(a2dDem.gt.0.0.AND.a2dBeta.le.0.0)
		a2dBeta=0.00001
	ENDWHERE



!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing Ct and Cf file (unit=5) 
!-----------------------------------------------------------------------------------

	open(unit=7,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.ct.txt',status='old',err=518)
	do i=1,6
		read(7,*) 
	enddo	
	do i=1,iRows
		read(7,*)(a2dCt(iRows+1-i,j),j=1,iCols)
	enddo
	close(7)

	!If ct file is not available uses the one from command line
    if(1.eq.2)then
518		write(*,*)'**********************'
		write(*,*)'ct FILE not Found. Continuum uses ct from command line'
		write(*,*)'**********************'
		write(*,*)''
	endif

	open(unit=7,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.cf.txt',status='old',err=519)
	do i=1,6
		read(7,*) 
	enddo	
	do i=1,iRows
		read(7,*)(a2dCf(iRows+1-i,j),j=1,iCols)
	enddo
	close(7)

	!If ct file is not available uses the one from command line
    if(1.eq.2)then
519		write(*,*)'**********************'
		write(*,*)'cf FILE not Found. Continuum uses cf from command line'
		write(*,*)'**********************'
		write(*,*)''
	endif

!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing Uc and Uh file (unit=5) 
!-----------------------------------------------------------------------------------

	open(unit=7,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.uc.txt',status='old',err=520)
	do i=1,6
		read(7,*) 
	enddo	
	do i=1,iRows
		read(7,*)(a2dCappaC(iRows+1-i,j),j=1,iCols)
	enddo
	close(7)

	!If ct file is not available uses the one from command line
    if(1.eq.2)then
520		write(*,*)'**********************'
		write(*,*)'uc FILE not Found. Continuum uses uc from command line'
		write(*,*)'**********************'
		write(*,*)''
	endif

	open(unit=7,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.uh.txt',status='old',err=521)
	do i=1,6
		read(7,*) 
	enddo	
	do i=1,iRows
		read(7,*)(a2dCappaV(iRows+1-i,j),j=1,iCols)
	enddo
	close(7)

	!If ct file is not available uses the one from command line
    if(1.eq.2)then
521		write(*,*)'**********************'
		write(*,*)'uh FILE not Found. Continuum uses uh from command line'
		write(*,*)'**********************'
		write(*,*)''
	endif

!-----------------------------------------------------------------------------------
where(a2iMask.eq.-9999)
a2dCTime=-9999.0
!a2dDem=-9999.0
!a2iPun=-9999.0
a2dCurNum=-9999.0
!a2iChoice=-9999.0
endwhere
!-----------------------------------------------------------------------------------



return
end subroutine