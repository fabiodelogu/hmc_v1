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
 
subroutine ReadStateDigheLaghi(iRows,iCols,xstring)

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
OPEN(20,file=sPathState,status='old',err=58)
  
 
DO i=1,iNdam
		!----coord i e j della sezione di chiusura o dell'opera di sbarramento
		READ(20,*)
		!----Codice della diga sulla matrice choiche
		READ(20,*)
		!----Volume massimo contenibile dalla diga
		READ(20,*)
		!----Volume iniziale contenibile in diga (se non letto da file)
		READ(20,*)a1dDamVolume(i)
		!---spazio
	    READ(20,*)

ENDDO
DO i=1,iNlake
		!---nome sezioni di chiusura
	    READ(20,*)
		!----coord i e j della sezione di chiusura o dell'opera di sbarramento
		READ(20,*)
		!----Codice del lago sulla matrice choich
		READ(20,*)
		!----Volume minimo per avere portata
		READ(20,*)
		!----Volume iniziale del lago
		READ(20,*)a1dVlake(i)

ENDDO
CLOSE(20)



58 return
end subroutine