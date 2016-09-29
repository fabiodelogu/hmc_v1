!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine SaveStateDigheLaghi
!	Salva lo stato dei laghi a monte delle dighe e dei laghi
!***********************************************************************************  
 
subroutine SaveStateDigheLaghi(iRows,iCols,xstring)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i,iFlagInterp
CHARACTER*50 info,text
CHARACTER*12 xstring
INTEGER*4 iLStr,iPathLenght,inC,inCt,ii
REAL*8 dDintegr !Dt integrazione Hcdrift
CHARACTER*500 sPathState

iPathLenght = iLStr(sPathResults)
700 FORMAT (1000(f9.2,1x))

sPathState=sPathResults(1:iPathLenght)//'StateStorage'//sBasin(1:iNameLenght)//xstring//'.txt'

!Opening info dighe
OPEN(20,file=sPathState,status='unknown')
  
 
DO i=1,iNdam
		!----coord i e j della sezione di chiusura o dell'opera di sbarramento
		WRITE(20,*)a2dXYDam(i,2),a2dXYDam(i,1)
		!----Codice della diga sulla matrice choiche
		WRITE(20,*)a1dCodeDam(i)
		!----Volume massimo contenibile dalla diga
		WRITE(20,*)a1dVdamMax(i)
		!----Volume iniziale contenibile in diga (se non letto da file)
		WRITE(20,*)a1dDamVolume(i)
		!---nome sezioni di chiusura
	    WRITE(20,*)

ENDDO
DO i=1,iNlake
		!---nome sezioni di chiusura
	    WRITE(20,*)
		!----coord i e j della sezione di chiusura o dell'opera di sbarramento
		WRITE(20,*)a2dXYLake(i,2),a2dXYLake(i,1)
		!----Codice del lago sulla matrice choich
		WRITE(20,*)a1dCodeLake(i)
		!----Volume minimo per avere portata
		WRITE(20,*)a1dVlakeMin(i)
		!----Volume iniziale del lago
		WRITE(20,*)a1dVlake(i)

ENDDO
CLOSE(20)



return
end subroutine