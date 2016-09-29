!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************


!------------------------------------------------------------------
!	restituisce: a1dQ_sv=portata dagli organi di scarico delle dighe
!------------------------------------------------------------------

 SUBROUTINE SvasoDamVolume(dVdam,dVdamMax,iNdighe,dDintegr,a1dQ_sv)
	include 'DeclarationH.f90'


	integer*4 iNdighe,nn,i,ii
    integer*4 name_lenght
	REAL*8 dDintegr  !Intervallo di integrazione del routing in secondi
	real*8 dVdam(iNdighe),dVdamMax(iNdighe),a1dQ_sv(iNdighe)
	real*8 dTV !Percentuale di volume oltre il quale inizia lo svaso dagli organi di scarico
!------------------------------------------------------------------

!---Inizializzazione a1dQ_sv (portata da svasare a valle della diga)
    a1dQ_sv=0
!------------------------------------------------------------------
	dTV=0.90
!   Calcolo il volume svasato se supero il 0.95V max contenibile in diga
	do i=1,iNdighe
		if(a1dCoefDighe(i).gt.0.0.and.dVdam(i).gt.dTV*dVdamMax(i))then
			ii=a2dXYDam(i,2)
			j=a2dXYDam(i,1)
			a1dQ_sv(i)=(dVdam(i)-dTV*dVdamMax(i))*a1dCoefDighe(i) !In m3/s
			!Confronto con Q critica scarico laterale
			if(a1dQ_sv(i).gt.a1dQ_sLC(i))then
				a1dQ_sv(i)=a1dQ_sLC(i)
			endif

			!if(i.eq.1)then
			!	write(*,*)a1dQ_sv(i),dVdam(i)
			!endif

			dVdam(i)=dVdam(i)-a1dQ_sv(i)*dDintegr !Aggiorno il volume della diga m3
			a1dQ_sv(i)=a1dQ_sv(i)*1000*dDintegr/(a2dAreaCell(ii,j)) !In mm
			!Controllo che il volume non sia negativo
			if(dVdam(i).lt.0.0)then
				dVdam(i)=0.0
				a1dQ_sv(i)=0.0
			endif
		endif
	enddo

	RETURN
	END

