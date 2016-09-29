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
 
subroutine ReadInfoLaghi(iRows,iCols)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i
CHARACTER*50 info,text
CHARACTER*12 xstring
INTEGER*4 iLStr,iPathLenght
REAL*8 dDintegr !Dt integrazione Hcdrift

!Opening info dighe
OPEN(20,file=sFileInfoLaghi,status='old')
  
 
READ(20,*)

DO i=1,iNlake
		READ(20,*)
		!---nome sezioni di chiusura
	    READ(20,*)
		!----coord i e j della sezione di chiusura o dell'opera di sbarramento
		READ(20,*)a2dXYLake(i,2),a2dXYLake(i,1)
		a2dXYLake(i,2)=iRows-a2dXYLake(i,2)+1
		IF(a2dXYLake(i,2).le.0)THEN
			write(*,*)'*************'
			write(*,*)'Lago numero ',i,' numero righe negative; Esco'
			write(*,*)'*************'
			STOP
		ENDIF
		!----Codice del lago sulla matrice choich
		READ(20,*)a1dCodeLake(i)
		!----Volume minimo per avere portata
		READ(20,*)a1dVlakeMin(i)
		!----Volume iniziale del lago
		READ(20,*)a1dVlake(i)
		!----Costante di svuotamento del lago in 1/h
		READ(20,*)a1dCostLaghi(i)
	  
ENDDO
CLOSE(20)



return
end subroutine