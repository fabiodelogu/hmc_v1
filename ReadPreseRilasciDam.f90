!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine lettura prese rilasci
!***********************************************************************************  
 
subroutine ReadPreseRilasciDam(iRows,iCols,dSimLength,d)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i,j,iFlagInterp,it,iit,iiit,iShift,in,iPr,iRank
CHARACTER*500 sF,sFileTurbinati,sName
CHARACTER*12 xstring
CHARACTER*500 sNmTmp,sTemp1,sTemp2
INTEGER*4 iLStr,iPathLenght,iSezLenght,iL1,iL2
REAL*8 dSimLength,dTemp,d !Lunghezza simulazione, variabile temp, dt delle turbinate
REAL*8 dHyTemp(int(dSimLength)) !vettore idrogramma temporaneo
REAL*8 dTurb(int(dSimLength))
iPathLenght = iLStr(sPathTurb)

dTurb=0

DO it=1,iNril 
	iit=0
	dHyTemp=0.0
	sName=a1sNameRilasci(it)
	iSezLenght = iLStr(sName)
	open(31,file=sPathTurb(1:iPathLenght)//'plantdata_'//sName(1:iSezLenght)//'.txt',status='old',err=13)
12	iit=iit+1
	if(iit.gt.dSimLength)Then !Se ho più dati della lunghezza della simulazione mi fermo
		go to 112
	endif

	!read(31,*,end=111)dTemp
	read(31,'(a50)',end=112)sNmTmp
	do in=1,50
		if(sNmTmp(in:in).eq.'N')then
			go to 12
		endif
	enddo
	read(sNmTmp,*)dTemp

	!Coordinate rilasci
	i=a2dXYRilascio(it,2)
	j=a2dXYRilascio(it,1)
	if(dTemp.ge.0.0)then
		!Lascio in m^3/s
		dHyTemp(iit)=dTemp*1.0 !levare il *0.0
		!Trasformo in mm/h
		!a1dIdroRila(it,iit)=dTemp*1000*d/(dCelLat*dCelLon)*1.0 !levare il *0.0
		a1dIdroRila(it,iit)=dTemp*1000*d/a2dAreaCell(i,j) 
		!dHyTemp(iit)=50*1000*d/(dCelLat*dCelLon) !test
	endif
	go to 12
112	CLOSE(31)
	
	IF(1.eq.2)Then !Se non c'è il file del turbinato lo pongo a 0.0		
13		a1dIdroRila(it,:)=0.0
		!Se ho il rilascio massimo disponibile lo uso all'80%
		if(a1dQmaxRil(it).gt.0.0)then
			!Lunghezza turbinata
			iit = int(dSimLength)
			!Coordinate rilasci
			i=a2dXYRilascio(it,2)
			j=a2dXYRilascio(it,1)
			dHyTemp(:)=a1dQmaxRil(it)*0.8
			a1dIdroRila(it,:)=a1dQmaxRil(it)*1000*d/a2dAreaCell(i,j)*0.8
		endif
		write(*,*)'Plant data ',sName(1:iSezLenght),' not found'
	ENDIF
	!Smoothing della la turbinata
	iRank = SIZE (a1dIdroRila,dim=2)
	CALL SmoothSerie(a1dIdroRila(it,:),iRank,3)
	iRank = SIZE (dHyTemp,dim=1)
	CALL SmoothSerie(dHyTemp,iRank,3)
	
	!!Associo alle prese i corrispondenti rilasci pesati e shiftati
	DO iPr=1,iNprese
	!Shift temporale
	sTemp1=a1sNameRilasci(it)
	sTemp2=a1sNamePrese(iPr)
	iL1=iLStr(sTemp1)
	iL2=iLStr(sTemp2)
	if(sTemp1(1:iL1).eq.sTemp2(1:iL2))then
		iShift=nint(a1dTcorrprese(iPr)*60/d)
		DO iiit=1,iit-1
			
			if(iiit-iShift.ge.1)then
				a1dIdroPrese(iPr,iiit-iShift)=dHyTemp(iiit)*a1dPesoPresa(iPr) !in m^3/s
			endif
			!dTurb(iiit+iShift)=dTurb(iiit+iShift)+dHyTemp(iiit)*(dCelLat*dCelLon)/(1000*d)
		ENDDO
		!Gli ultimi iShift step ipotizzo l'ultimo valore di turbinata costante
		DO iiit=iit-1-iShift,iit-1
			a1dIdroPrese(iPr,iiit)=dHyTemp(iit-1)*a1dPesoPresa(iPr)
		ENDDO
	ENDIF
	ENDDO




ENDDO




return
end subroutine