!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine ReadInfoDataDam
!	Legge il file con le informazioni delle sezioni se ho le dighe: basinDem.info.txt
!***********************************************************************************  
 
subroutine ReadInfoS3M(iRows,iCols)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols,i,j
CHARACTER*500 sPathNature
INTEGER*4 iLStr,iPathLenght
REAL*8 dDintegr !Dt integrazione Hcdrift
!Opening file info snow
iFlagSnow=1

OPEN(4,file=sFileInfoSnow,status='old',err=513)
  
!Screen message
write(*,*)
write(*,*)'File info for SNOW found! S3M WILL be used'
write(*,*)

!Reading S3M.info.txt
READ(4,*)
READ(4,*)dTrif		 ! Temperatura di riferimento
READ(4,*)
READ(4,*)dRoS0		 ! densità neve fresca
READ(4,*)
READ(4,*)dRoMax		 ! densità massima neve 
READ(4,*)
READ(4,*)dThres		! soglia cumulata di Precip giornaliera
READ(4,*)
READ(4,*)dExpRoLow		! Esponente compattazione densita, assenza melting
READ(4,*)
READ(4,*)dExpRoHigh		! Esponente compattazione densita, con melting elevato
READ(4,*)
READ(4,*)snodata	 ! nodata value
READ(4,*)
READ(4,*)sPathNature	 ! nome del path Carta della Natura
READ(4,*)
READ(4,*)(sMcStag(i),i=1,4) !Coeff. Melting Stagionali
CLOSE(4)
!-----------------------------------------------------------------------------------
!	Opening, Reading and Closing Nature Map for Glaciers (unit=5) 
!-----------------------------------------------------------------------------------
iPathLenght = iLStr(sPathNature)

open(unit=5,file=sPathNature(1:iPathLenght)//sBasin(1:iNameLenght)//'.nature.txt',status='old',err=512)
do i=1,6
	read(5,*) 
enddo	
do i=1,iRows
	read(5,*)(a2dNature(iRows+1-i,j),j=1,iCols)
enddo
close(5)

! Non ho file per funzionamento parte neve
if(1.eq.0)then 
513	iFlagSnow=0
write(*,*)
write(*,*)'File info for SNOW NOT found! S3M will NOT be used'
write(*,*)
endif

512 return
end subroutine