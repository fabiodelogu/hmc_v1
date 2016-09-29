!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine read_info_data
!	Legge il file con le info delle dighe con i laghi a monte: sBasinDamDigheInfo.txt
!***********************************************************************************  
 
subroutine ReadInfoDigheLaghi(iRows,iCols)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i,iFlagInterp
CHARACTER*50 info,text
CHARACTER*12 xstring
CHARACTER*50 sFileInvVol
INTEGER*4 iLStr,iPathLenght,inC,inCt,ii
REAL*8 dDintegr !Dt integrazione Hcdrift

!Opening info dighe
OPEN(20,file=sFileInfoDighe,status='old')
  
 
READ(20,*)
READ(20,*)
READ(20,*)
inCt=1
DO i=1,iNdam
		!---nome sezioni di chiusura
	    READ(20,*)
		!----coord i e j della sezione di chiusura o dell'opera di sbarramento
		READ(20,*)a2dXYDam(i,2),a2dXYDam(i,1)
		a2dXYDam(i,2)=iRows-a2dXYDam(i,2)+1
		IF(a2dXYDam(i,2).le.0)THEN
			write(*,*)'*************'
			write(*,*)'Diga numero ',i,' numero righe negative; Esco'
			write(*,*)'*************'
			STOP
		ENDIF
		!----numero di centrali a valle della diga
		READ(20,*)inC
		!----Codice della diga sulla matrice choiche
		READ(20,*)a1dCodeDam(i)
		!----Volume massimo contenibile dalla diga
		READ(20,*)a1dVdamMax(i)
		!----Volume iniziale contenibile in diga (se non letto da file)
		READ(20,*)a1dDamVolume(i)
		!----Portata critica dallo scarico laterale o altro tipo
		READ(20,*)a1dQ_sLC(i)
		!----Lunghezza equivalente scarico superficiale della diga
		READ(20,*)adL(i)
		!----Altezza massima di invaso
		READ(20,*)a1dHmax(i)
		!----Coefficiente serbatoio lineare equivalente allo scarico superficiale della diga
		READ(20,*)a1dCoefDighe(i)
		!----path e nome file curve invaso volume
		READ(20,*)sFileInvVol
		!----Leggo la curva invaso volume
		CALL ReadCurveInvasoVolumiDam(sFileInvVol,i)
	 		
		DO ii=inCt,inC+inCt-1
		 !----nome della centrale
	     READ(20,*)a1sNameTurbinate(ii)
		 !----coordinate della centrale
	     READ(20,*) a2dXYCen(ii,2),a2dXYCen(ii,1)
         a2dXYCen(ii,2)=iRows-a2dXYCen(ii,2)+1
		 !----tempo di corrivazione diga-centrale (dismesso)
	     READ(20,*) a1dTcorrturb(ii)
		 !----portata max impianto
	     READ(20,*)a1dQturbMax(ii)
		 !----flag
	     READ(20,*) 
		 a1iDigaCentrale(ii)=i !Diga corrispondente alla centrale
		 !write(*,*)a1sNameTurbinate(ii),' ',a1dQturbMax(ii)
		 !----flag se conosco o no le portate turbinate
		 !CALL checkturbinate(pathevento,filename_centrale(i,ii),i,ii)
		 
		 !IF(flag_turb(i,ii).eq.1)n_turbinate=n_turbinate+1
		ENDDO
		inCt=inCt+inC
		
	    READ(20,*)
ENDDO
CLOSE(20)



return
end subroutine