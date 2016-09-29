!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************


!------------------------------------------------------------------
!	Gestisce la dinamica dei laghi come serbatoi
!------------------------------------------------------------------

 SUBROUTINE TankLaghi(dVlake,dVlakeMin,iNlaghi,dDintegr,dQ_laghi,dCostLaghi)
	include 'DeclarationH.f90'


	integer*4 iNlaghi,nn,iL,i,j,ii,jj,iii,jjj
    integer*4 name_lenght
	REAL*8 dDintegr  !Intervallo di integrazione del routing in secondi
	real*8 dVlake(iNlaghi),dVlakeMin(iNlaghi),dQ_laghi(iNlaghi),dCostLaghi(iNlaghi)

!------------------------------------------------------------------

!---Inizializzazione a1dQ_sv (portata da svasare a valle della diga)
    a1dQ_sv=0
!------------------------------------------------------------------
!   Calcolo il volume svasato se supero il V max contenibile in diga
	do iL=1,iNlaghi

		if(dVlake(iL).gt.dVlakeMin(iL))then
			i=a2dXYLake(iL,2)
			j=a2dXYLake(iL,1)
			dQ_laghi(iL)=(dVlake(iL)-dVlakeMin(iL))*(dCostLaghi(iL)/3600)/(a2dAreaCell(i,j))*3600*1000 !In mm/h, portata in uscita dal lago
			dVlake(iL)=dVlake(iL)-(dVlake(iL)-dVlakeMin(iL))*(dCostLaghi(iL)/3600)*dDintegr !Svuoto il lago
			
			ii=int((a2iPun(i,j)-1)/3)-1
			jj=a2iPun(i,j)-5-3*ii
			iii=i+ii
			jjj=j+jj
			!Equazione riempimento dei laghi
			IF(a1dCodeLake(iL).gt.0)THEN 
				WHERE (a2iChoice.eq.a1dCodeLake(iL).AND.a2dDem.GT.0.0)
					!Converto il volume in un livello medio
					a2dHydro=dVlake(iL)/(a1dCodeLake(iL)*a2dAreaCell(i,j))*1000 !in mm               
				ENDWHERE
			ENDIF
			!Sommo il contributo dei laghi ad a2dDeepFlow, che contribuisce ad a2dIntensity
			!nella subroutine HortonMatrixHcdrift
			a2dDeepFlow(iii,jjj)=a2dDeepFlow(iii,jjj)+dQ_laghi(iL)*dDintegr/3600 !mm nel dt
		endif
		
	enddo

	RETURN
	END

