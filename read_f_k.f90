!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!     Subroutine READf_k
!     Read a2dK and f
!***********************************************************************************    

SUBROUTINE ReadFK(iRows,iCols)

IMPLICIT NONE
INCLUDE 'DeclarationH.f90'

INTEGER i,j,n
INTEGER iRows,iCols
REAL*8 kcn(100),fcn(100),a2dCostF2(iRows,iCols)
INTEGER*4 iLStr,iPathLenght
!***********************************************************************************

iPathLenght = iLStr(sPathLandData)

OPEN (unit=11,file=sPathLandData(1:iPathLenght)//'valori_fo_noIa_AMC2.txt',status='old',err=888)   !ATT: quando cambi f_k cambia anche CN
DO i=1,100
	READ(11,*)kcn(i),fcn(i) !la prima colonna si può levare (e variare anche il file di conseguneza)
END DO
CLOSE (11)
!***********************************************************************************
 
IF(1.eq.0)THEN
888 write(*,*)' '
	write(*,*)'File CN to fo conversion not found!'
	write(*,*)'Use standard CN to fo conversion'
	write(*,*)' '
	fcn=(/603.3,572.9,543.8,515.9,489.2,463.6,439.2,415.9,393.6,372.3, &
		352.0,332.7,314.3,296.8,280.2,264.4,249.4,235.2,221.8,209.1,197.1, &
		185.8,175.1,165.1,155.6,146.8,138.4,130.7,123.4,116.6,110.3,104.4, &
		98.9,93.8,89.1,84.7,80.7,77.0,73.7,70.5,67.7,65.1,62.8,60.6,58.7,56.9, &
		55.3,53.9,52.6,51.4,50.4,49.4,48.6,47.8,47.1,46.4,45.9,45.3,44.8,44.3, &
		43.8,43.3,42.8,42.3,41.8,41.2,40.7,40.1,39.5,38.8,38.1,37.3,36.5,35.7, &
		34.7,33.8,32.8,31.7,30.5,29.4,28.1,26.8,25.5,24.1,22.6,21.2,19.7,18.1, &
		16.5,14.9,13.3,11.7,10.0,8.4,6.8,5.2,3.6,2.1,0.6,0.1/)

ENDIF
    

FORALL(i=1:iRows,j=1:iCols,a2dCon(i,j).lt.1.or.a2dCon(i,j).gt.99) a2dCon(i,j)=1
		
FORALL(i=1:iRows,j=1:iCols)
	a2dCostF(i,j)=fcn(a2dCon(i,j))
END FORALL


!***********************************************************************************
! Preparo le costanti del metodo di Horton

!Calcolo f1
a2dCostF1=a2dCf*a2dCostF 
!Esponente di Horton per I>g
a2dCostChFix=((1-a2dCt)*a2dCostF+a2dCt*a2dCostF1)/((1-a2dCt)*a2dS)
!Termine correttivo di Horton per Ct num. 1
a2dC1=a2dCostF1*a2dCt/(1-a2dCt)
!Termine correttivo di Horton per Ct num. 2
a2dF2=a2dCostF1/(1-a2dCt)


RETURN
END

!***********************************************************************************
