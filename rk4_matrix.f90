!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!*****************************************************************************************
!	Subroutine per un Runge - Kutta di quarto ordine 
!*****************************************************************************************

! Changes
! 12/03/2014: Insert dsqrt and dexp instead of sqrt and exp

	subroutine rk4_matrix(iRows,iCols,t,tstep,TD,P_it,Ch,lambda,beta,Rn,UmR,WR,Ta,Trif,ea,rhoa,PP,Ts_up)

	implicit none
    include 'DeclarationH.f90'
!-----------------------------------------------------------------------------------------
!	Dichiarazioni delle variabili
!-----------------------------------------------------------------------------------------

	integer iRows,iCols
    integer i, n

	real*8 h, t, tstep, Trif
	real*8 Ts_up(iRows,iCols)
    real*8 k1(iRows,iCols), k2(iRows,iCols), k3(iRows,iCols), k4(iRows,iCols)
	real*8 Rn(iRows,iCols), WR(iRows,iCols), Ta(iRows,iCols), P_it(iRows,iCols), rhoa(iRows,iCols), &
			TD(iRows,iCols), UmR(iRows,iCols), PP(iRows,iCols)
	real*8 Ch(iRows,iCols),lambda(iRows,iCols),beta(iRows,iCols),ea(iRows,iCols)
	real*8 omega,pi,cp
	real*8 dTsmax !delta max for applying Runge Kutta

	parameter(pi=3.14, cp=1004.0, omega=1.0/(60.0*60.0*24.0))	! Pi greco [-], Calore specifico a pressione costante [J/kg/a2dK], Lunghezza del giorno [a2dS]
	dTsmax=2
	h=tstep/2.0

	WHERE (a2dCTime.gt.0.0)

	Ts_up=0.0
	
	Ts_up=a2dTS
	k1 = tstep *( 2*dsqrt(pi*omega)/P_it*(Rn-rhoa*cp*Ch*WR*(Ts_up-Ta)-rhoa*lambda*Ch*WR*beta*(0.611*dexp(17.3*(Ts_up-Trif)/ &
			(237.3+Ts_up-Trif))-ea)/PP*0.622)-2*3.14*omega*(Ts_up-TD))
	WHERE (k1.gt.dTsmax)
		k1 = dTsmax
	ENDWHERE
	
	Ts_up=a2dTS+k1/2.0
	k2 = tstep *( 2*dsqrt(pi*omega)/P_it*(Rn-rhoa*cp*Ch*WR*(Ts_up-Ta)-rhoa*lambda*Ch*WR*beta*(0.611*dexp(17.3*(Ts_up-Trif)/ &
			(237.3+Ts_up-Trif))-ea)/PP*0.622)-2*3.14*omega*(Ts_up-TD))
	
	WHERE (k2.gt.dTsmax)
		k2 = dTsmax
	ENDWHERE
	Ts_up=a2dTS+k2/2.0
	k3 = tstep *( 2*dsqrt(pi*omega)/P_it*(Rn-rhoa*cp*Ch*WR*(Ts_up-Ta)-rhoa*lambda*Ch*WR*beta*(0.611*dexp(17.3*(Ts_up-Trif)/ &
			(237.3+Ts_up-Trif))-ea)/PP*0.622)-2*3.14*omega*(Ts_up-TD))
	WHERE (k3.gt.dTsmax/2)
		k3 = dTsmax/2
	ENDWHERE

	Ts_up=a2dTS+k3
	k4 = tstep *( 2*dsqrt(pi*omega)/P_it*(Rn-rhoa*cp*Ch*WR*(Ts_up-Ta)-rhoa*lambda*Ch*WR*beta*(0.611*dexp(17.3*(Ts_up-Trif)/ &
			(237.3+Ts_up-Trif))-ea)/PP*0.622)-2*3.14*omega*(Ts_up-TD))
	
    WHERE (k4.gt.dTsmax)
		k4 = dTsmax
	ENDWHERE

	Ts_up = a2dTS + (k1 + (2.*(k2 + k3)) + k4)/6.0
	
	ENDWHERE

	return
	end subroutine