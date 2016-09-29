!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************


!***********************************************************************************
!     Subroutine Thermal Inertia
!***********************************************************************************

! Changes
! 12/03/2014: Insert dsqrt and dexp instead of sqrt and exp

	subroutine thermal_inertia_MATRIX(iRows,iCols,VV,P_it)

	implicit none
	include 'DeclarationH.f90'

	integer iRows,iCols
	real*8 VV(iRows,iCols), P_it(iRows,iCols), C(iRows,iCols)
	real*8 rhot, cpt, rhow, cpw, rhod
	real*8 kq, kw, ko, ks, kdry, ksat
	real*8 n
	real*8 a2dQ
	real*8 Ke(iRows,iCols), kter(iRows,iCols),VVtmp(iRows,iCols)

!-----------------------------------------------------------------------------------
!	Costanti utilizzate nella subroutine
!-----------------------------------------------------------------------------------

!	Densità del terreno [kg m^-3]
rhot=2700
!	Calore Specifico [J kg^-1 a2dK^-1]
cpt=733
!	Densità dell'acqua [kg m^-3]
rhow=1000
!	Calore specifico dell'acqua [J kg^-1 a2dK^-1]
cpw=4186
!	Conduttività termica del quarzo [a2dW m^-1 a2dK^-1]
kq=7.7
!	Conduttività termica dell'acqua [a2dW m^-1 a2dK^-1]
kw=0.57
!	Conduttività termica degli altri minerali [a2dW m^-1 a2dK^-1]
!ko=2 !Casentino
ko=4 !Orba

!-----------------------------------------------------------------------------------
!	Porosità del terreno [%]
!n=0.4
!	Contenuto di quarzo [%]
!a2dQ=0.4 !Casentino
a2dQ=0.5 !Orba
!-----------------------------------------------------------------------------------
C=0.0
VVtmp=0.0
VVtmp=VV
WHERE (a2dCTime.gt.0.0.and.VV.gt.10.0)
	VVtmp=1.0
ENDWHERE
WHERE (a2dCTime.gt.0.0.and.VV.lt.0.0)
	VVtmp=0.0
ENDWHERE



WHERE (a2dCTime.gt.0.0)
!	Heat Capacity [J K^-1 m^-3]
	C=(1-dPorosity)*rhot*cpt+dPorosity*VVtmp*rhow*cpw
ENDWHERE

!	Densità secca [kg m^-3]
rhod=(1-dPorosity)*rhot

!	Conduttività termica secca [a2dW m^-1 a2dK^-1]
kdry=(0.135*rhod+64.7)/(rhot-0.947*rhod)

!	Conduttività termica nei solidi [a2dW m^-1 a2dK^-1]	
ks=kq**(a2dQ)*ko**(1-a2dQ)

!	Conduttività termica satura
ksat=ks**(1-dPorosity)*kw**(dPorosity)

WHERE (a2dCTime.gt.0.0)

!	Kersten number (funzione solo del grado di saturazione VV per terreni sottili)
	WHERE(VVtmp.ge.0.1) 
		Ke=log10(VVtmp)+1
	ELSEWHERE
		Ke=log10(0.1)+1
	ENDWHERE

!	Conduttività termica [W m^-1 K^-1]
	kter=Ke*(ksat-kdry)+kdry

!	Inerzia Termica [J m^-2 K S^-(1/2)]
	P_it=dsqrt(C*kter)

!write(*,*) Ke, kter, P_it

ENDWHERE

return
end subroutine
