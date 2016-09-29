!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine lettura turbinati
!***********************************************************************************  
 
subroutine ReadTurbinatiDam(iRows,iCols,dSimLength,d)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i,j,iFlagInterp,it,iit,iiit,iShift,in,iRank
CHARACTER*500 sF,sFileTurbinati,sName
CHARACTER*12 xstring
CHARACTER*50 sNmTmp
INTEGER*4 iLStr,iPathLenght,iSezLenght
REAL*8 dSimLength,dTemp,d !Lunghezza simulazione, variabile temp, dt delle turbinate
REAL*8 dHyTemp(int(dSimLength)) !vettore idrogramma temporaneo
REAL*8 dTurb(int(dSimLength))
iPathLenght = iLStr(sPathTurb)

dTurb=0

DO it=1,iNcentr 
	iit=0
	dHyTemp=0.0
	sName=a1sNameTurbinate(it)
	iSezLenght = iLStr(sName)
	open(31,file=sPathTurb(1:iPathLenght)//'plantdata_'//sName(1:iSezLenght)//'.txt',status='old',err=12)
11	iit=iit+1
	if(iit.gt.dSimLength)Then !Se ho più dati della lunghezza della simulazione mi fermo
		go to 111
	endif

	!read(31,*,end=111)dTemp
	read(31,'(a50)',end=111)sNmTmp
	do in=1,50
		if(sNmTmp(in:in).eq.'N')then
			go to 11
		endif
	enddo
	read(sNmTmp,*)dTemp
	!Trasformo in mm/h
	i=a2dXYCen(it,2)
	j=a2dXYCen(it,1)
	if(dTemp.ge.0.0)then
		!dHyTemp(iit)=dTemp*1000*d/(dCelLat*dCelLon)
		dHyTemp(iit)=dTemp*1000*d/a2dAreaCell(i,j) 
		!dHyTemp(iit)=50*1000*d/(dCelLat*dCelLon) !test
	endif
	go to 11
111	CLOSE(31)

	!Smoothing della la turbinata
	iRank = SIZE (dHyTemp,dim=1)
	CALL SmoothSerie(dHyTemp,iRank,2)
	
	DO iiit=1,iit-1
		a1dIdroTurbinate(it,iiit)=dHyTemp(iiit)
		!dTurb(iiit)=dTurb(iiit)+dHyTemp(iiit)*(dCelLat*dCelLon)/(1000*d)
		dTurb(iiit)=dTurb(iiit)+dHyTemp(iiit)*a2dAreaCell(i,j)/(1000*d) 
	ENDDO

	IF(1.eq.2)Then !Se non c'è il file del turbinato lo pongo a 0.0
12		a1dIdroTurbinate(it,:)=-9999
	ENDIF



ENDDO




return
end subroutine