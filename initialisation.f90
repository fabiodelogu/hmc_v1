!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!     Subroutine Initialisation
!***********************************************************************************  
 
SUBROUTINE initialisation
implicit none
include 'DeclarationH.f90'
integer i,j
real*8 dem_max

!	Inizializzazione del volume iniziale (a2dCon 0.2 abbiamo imposto che almeno per un EF pari a 0.45)
!	qui si è scelto un valore indicativo non quello risultante dall'inversione dell'arctg che identifica
!	la forma di EF. La forma di EF al momento è da rimettere a posto!


	a3dTemp24=0.0
	a3dTMarked=0.0
	 				
	a2dCumRain=0.0
	a2dEsf=0.0
	vol_sot=0.0
	a2dRetention=0.001


	a2dRain=0.0
	a2dEvapot=0.0
    a2dTs=0.0

	dem_max=maxval(maxval(a2dDem,DIM = 1),DIM = 1)

    a2dV=dCPI/2*a2dS+dCPI/2*a2dS*(dem_max-a2dDem)/dem_max
	a2dVwt=a2dVwt*(dem_max-a2dDem)/dem_max
	a2dHydro=0.0*(dem_max-a2dDem)/dem_max
	WHERE (a2dV.lt.0.0)
		a2dV=0.0
		a2dVwt=0.0
    ENDWHERE
    WHERE (a2iChoice.lt.1)
		a2dHydro=0.0
    ENDWHERE
    WHERE (a2iChoice.eq.1)
		a2dV=a2dS*0.98
    ENDWHERE
RETURN
END

!***********************************************************************************
