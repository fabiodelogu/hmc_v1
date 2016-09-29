!***********************************************************************************   
!	Subroutine SubFlowHcdrift
!	Calcola la variazione del contenuto d'acqua di ogni cella del sBasin
!***********************************************************************************

	SUBROUTINE SubFlowHcdrift(iRows,iCols,d,dDintegr,dTotHypodFlow,dTotDeepFlow,dTotEsf)
	implicit none

	include 'DeclarationH.f90'
	
	integer iRows,iCols
	real*8 dTotHypodFlow,dTotDeepFlow,Kalpha_zero_mean
	real*8 d,f(iRows,iCols),dTotEsf
	real*8 Vloss_matrix(iRows,iCols),a2dVtmp(iRows,iCols)
	real*8 dDintegr !Dt del routing
	real*8 Vprec(iRows,iCols),err(iRows,iCols)
	real*8 dDth !Dt su cui devo integrare il deflusso ipodermico
	real*8 dRate
	integer*4 i,j,ii,jj,iii,jjj,ia,ja,oldj,oldi


!Pongo dDth=dDintegr, ma potrei differenziare canale e versante	
	dDth=dDintegr


	dTotEsf=0.0
	dTotHypodFlow=0.0
    !dTotDeepFlow=0.0
	a2dVtmp=0.0
	f=0.0

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
					dRate=dsin(a2dBeta(i,j))
					if(dRate.gt.1)dRate=0.99
					if(dRate.lt.0.1)dRate=0.1
					if(a2dVwtMax(i,j).eq.0)dRate=1 !Celle dove non ho falda
					IF(iii.ge.1.and.jjj.ge.1) a2dVtmp(iii,jjj)=a2dVtmp(iii,jjj)+vol_sot(i,j)*dRate

!	Volume perso verso gli strati profondi
					IF(iii.ge.1.and.jjj.ge.1) a2dVLoss(iii,jjj)=a2dVLoss(iii,jjj)+vol_sot(i,j)*(1-dRate)
					if(iii.ge.1.and.jjj.ge.1)then
						f(i,j)=dRate
					endif
				
				ENDIF
				
				
			ENDIF
		    
		END DO
	END DO

	!aggiorno il V
	a2dV=a2dV+a2dVtmp
	a2dEsf=0.0
	!err=a2dV-Vprec-a2dVtmp
	!Vprec=a2dV

	WHERE (a2dDem.gt.0.0.and.a2dV.gt.a2dS)

			!a2dEsf=(a2dV-a2dS*(1-dCt+dCf*dCt)/(1-dCt+dCf+dCf*dCt))/(1000.0*dDth) !*d/3600	!in m/sec
			a2dEsf=(a2dV-a2dS)/(1000.0*dDth) !in m/sec
			!f=(a2dV-a2dS*(1-dCt+dCf*dCt)/(1-dCt+dCf+dCf*dCt))
					
			!a2dV=1.0*a2dS*(1-dCt+dCf*dCt)/(1-dCt+dCf+dCf*dCt)
			a2dV=1.0*a2dS 

	ENDWHERE
    i=a2dXYsections(1,2)
	j=a2dXYsections(1,1)
    dRate=dsin(a2dBeta(i,j))
	dVHypod=dVHypod+vol_sot(i,j)*dRate/1000*(a2dAreaCell(i,j))

	!err=a2dV+a2dEsf*1000.0*dDth-Vprec

    !Kalpha_zero_mean=(SUM(SUM(a2dKAlphaZero,DIM=1,MASK=a2dDem.GT.0.0)))/dBasinArea
	!dTotEsf=SUM(SUM(a2dEsf,DIM=1,MASK=a2dDem.gt.0.0.and.a2dV.gt.a2dS*(1-dCt+dCf*dCt)/(1-dCt+dCf+dCf*dCt))) !DIM=1 columns
	!dTotHypodFlow=SUM(SUM(a2dV,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	!dTotDeepFlow=SUM(SUM(a2dVLoss,DIM=1,MASK=a2dDem.GT.0.0)) 

	!write(996,*)Kalpha_zero_mean

	return
	END

!***********************************************************************************
