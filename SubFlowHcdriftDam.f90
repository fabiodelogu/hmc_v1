!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************   
!	Subroutine SubFlowHcdrift
!	Calcola la variazione del contenuto d'acqua di ogni cella del sBasin
!***********************************************************************************

	SUBROUTINE SubFlowHcdriftDam(iRows,iCols,d,dDintegr,dTotHypodFlow,dTotDeepFlow,dTotEsf)
	implicit none

	include 'DeclarationH.f90'
	
	integer iRows,iCols
	real*8 dTotHypodFlow,dTotDeepFlow,Kalpha_zero_mean
	real*8 d,f(iRows,iCols),dTotEsf,tmp(iRows,iCols)
	real*8 Vloss_matrix(iRows,iCols),a2dVtmp(iRows,iCols)
	real*8 dDintegr !Dt del routing
	real*8 dDth !Dt su cui devo integrare il deflusso ipodermico
	real*8 dRate
	integer*4 i,j,ii,jj,iii,jjj,ia,ja,oldj,oldi,iD


!Pongo dDth=dDintegr, ma potrei differenziare canale e versante	
	dDth=dDintegr


	dTotEsf=0.0
	dTotHypodFlow=0.0
	a2dVtmp=0.0
	f=0.0
    !dTotDeepFlow=0.0

	!Riempio la diga con il Volume di Horton della cella diga
	!Pongo poi a 0 il V(t) nelle celle corrispondenti alle dighe
	DO iD=1,iNdam
		i=a2dXYDam(iD,2)
		j=a2dXYDam(iD,1)
		a1dDamVolume(iD)=a1dDamVolume(iD)+a2dV(i,j)*(a2dAreaCell(i,j))/1000
		!write(*,*)a2dV(i,j),' ',a2dV(i,j)*(a2dAreaCell(i,j))/1000
		a2dV(i,j)=0.0
	ENDDO

	WHERE (a2dDem.gt.0.0.and.a2dV.lt.0.0) a2dV=0.0
    WHERE (a2dDem.gt.0.0.and.a2dVLoss.lt.0.0) a2dVLoss=0.0

    DO j=1,iCols
		DO i=1,iRows
		  
			IF (a2dDem(i,j).gt.0.0) then
			
!	Calcola il volume di uscita dalla cella nei due casi: a2dV > o < di a2dS;
!	Qsup è la portata che esce dalla parte superiore della cella e si 
!	aggiunge al deflusso superficiale               
!	Il contatore punta alla cella successiva(controllare se vale per 
!   l'ultima cella)	  

				ii=int((a2iPun(i,j)-1)/3)-1
				jj=a2iPun(i,j)-5-3*ii
				iii=i+ii
				jjj=j+jj

!-------
!   SENZA Water Table
!	Riempimento della cella successiva
				IF(iFlagDeepFlow.eq.0)THEN
					IF(iii.ge.1.and.jjj.ge.1) a2dVtmp(iii,jjj)=a2dVtmp(iii,jjj)+vol_sot(i,j)
!-------
				ELSE
!	Riempimento della cella successiva
					!dRate=(1-a2dCoeffResol(i,j))*dsin(a2dBeta(i,j))+a2dCoeffResol(i,j)
					dRate=dsin(a2dBeta(i,j))			
					if(dRate.gt.1)dRate=0.99
					if(dRate.lt.dRateMin)then
						dRate=dRateMin
					endif
					if(a2dVwtMax(i,j).eq.0)dRate=1 !Celle dove non ho falda

					!IF(iii.ge.1.and.jjj.ge.1) a2dV(iii,jjj)=a2dV(iii,jjj)+vol_sot(i,j)*(1-a2dKAlphaZero(i,j))*dcos(a2dAlpha(i,j))**2
					IF(iii.ge.1.and.jjj.ge.1) a2dVtmp(iii,jjj)=a2dVtmp(iii,jjj)+vol_sot(i,j)*dRate
!	Volume perso verso gli strati profondi
					!IF(iii.ge.1.and.jjj.ge.1) a2dVLoss(iii,jjj)=a2dVLoss(iii,jjj)+vol_sot(i,j)*a2dKAlphaZero(i,j)+vol_sot(i,j)*(1-a2dKAlphaZero(i,j))*dsin(a2dAlpha(i,j))**2
					IF(iii.ge.1.and.jjj.ge.1) a2dVLoss(iii,jjj)=a2dVLoss(iii,jjj)+vol_sot(i,j)*(1-dRate)
				ENDIF
				!if(iii.ge.1.and.jjj.ge.1)then
				!	f(i,j)=dRate
				!endif
				
			ENDIF
		    
		END DO
	END DO

	!aggiorno il V
	a2dV=a2dV+a2dVtmp

	a2dEsf=0.0
	!Drains to surface along channeled network because problems of spatial resolution
	!WHERE (a2dDem.gt.0.0.and.a2dV.le.a2dS.and.a2iChoice.eq.1)
	!		a2dEsf=a2dCoeffResol*(a2dV/a2dS)/(1000.0*3600) !in m/sec
	!		a2dV=a2dV-a2dCoeffResol*(a2dV/a2dS)*dDth/3600
	!ENDWHERE
	!Exfiltration calculation
	WHERE (a2dDem.gt.0.0.and.a2dV.gt.a2dS)
			a2dEsf=(a2dV-a2dS)/(1000.0*dDth) !in m/sec
			a2dV=1.0*a2dS 
	ENDWHERE

    !Kalpha_zero_mean=(SUM(SUM(a2dKAlphaZero,DIM=1,MASK=a2dDem.GT.0.0)))/dBasinArea
	!dTotEsf=SUM(SUM(a2dEsf,DIM=1,MASK=a2dDem.gt.0.0.and.a2dV.gt.a2dS*(1-dCt+dCf*dCt)/(1-dCt+dCf+dCf*dCt))) !DIM=1 columns
	!dTotHypodFlow=SUM(SUM(a2dV,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	!dTotDeepFlow=SUM(SUM(a2dVLoss,DIM=1,MASK=a2dDem.GT.0.0)) 


	!Da levare, Riempio la diga con il deflusso ipodermico della cella diga
	!DO iD=1,iNdam
	!	i=a2dXYDam(iD,2)
	!	j=a2dXYDam(iD,1)
		!Sommo tutto il V perchè è solo il contributo delle celle a monte, V è annullato 
		!ad inizio subroutine
		!a1dDamVolume(iD)=a1dDamVolume(iD)+a2dV(i,j)*(a2dAreaCell(i,j))/1000
	!ENDDO
	i=a2dXYsections(1,2)
	j=a2dXYsections(1,1)
    dRate=dsin(a2dBeta(i,j))
	dVHypod=dVHypod+vol_sot(i,j)*dRate/1000*(a2dAreaCell(i,j))

	return
	END

!***********************************************************************************
