!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine read_info_data
!	Legge il file con le info dei punti di immissione con eventuale bypass
!***********************************************************************************  
 
subroutine ReadInfoJoin(iRows,iCols)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i
CHARACTER*50 info,text
CHARACTER*12 xstring
INTEGER*4 iLStr,iPathLenght

!Opening info join
OPEN(20,file=sFileInfoJoin,status='old')
  
 
READ(20,*)

DO i=1,iNjoin 
		READ(20,*)
		!---nome di riferimento dell'unione dei corsi d'acqua
	    READ(20,*)
		!----coord i e j della sezione del corso d'acqua master
		READ(20,*)a2dXYMain(i,2),a2dXYMain(i,1)
		a2dXYMain(i,2)=iRows-a2dXYMain(i,2)+1
		!----coord i e j della sezione del corso d'acqua immissario che riceve acqua dal master
		READ(20,*)a2dXYImm(i,2),a2dXYImm(i,1)
		a2dXYImm(i,2)=iRows-a2dXYImm(i,2)+1
		!----coord i e j della sezione del corso d'acqua immissario dove arriva l'acqua dal master
		READ(20,*)a2dXYOut(i,2),a2dXYOut(i,1)
		a2dXYOut(i,2)=iRows-a2dXYOut(i,2)+1
		!----Soglia oltre la quale innescare lo svaso nell'immissario
		READ(20,*)a1dThreshLiv(i)
		  
ENDDO
CLOSE(20)



return
end subroutine