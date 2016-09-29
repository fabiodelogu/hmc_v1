!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!     Subroutine Storage
!     Calculate a2dS for the CN
!***********************************************************************************

	SUBROUTINE storage(iRows,iCols)
    
	implicit none
	include 'DeclarationH.f90'
	
	integer iRows,iCols
	integer i,j

	!Initialization
	a2dS=500
	a2dCon=1

	forall(i=1:iRows,j=1:iCols,a2dCurNum(i,j).gt.0.and.a2dCurNum(i,j).le.100) 
		a2dCurNum(i,j)=a2dCurNum(i,j)!-9
		a2dS(i,j)=(1000.0/a2dCurNum(i,j)-10)*25.4	
		a2dCon(i,j)=int(a2dCurNum(i,j))
	end forall


	WHERE (a2dDem.gt.0.0.and.a2dS.lt.1.0)
		a2dS=1.0
		a2dVwt=0.0
	ENDWHERE


	WHERE (a2dS.lt.0.0)
		a2dS=0.0
	ENDWHERE

   
!-----------------------------------------------------------------------------------
!	Chiamata alla Subroutine Readf_k
!	Legge a2dK(Cn),f(Cn) e crea la matrice per il sBasin dei kcn e fcn
!-----------------------------------------------------------------------------------
	call ReadFK(iRows,iCols)

	return
    END

!***********************************************************************************

