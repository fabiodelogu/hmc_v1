!***********************************************************************************   
!	Subroutine convolutionHyperCdrift
!	Esegue le equazioni di Hcdrift e la convoluzione
!***********************************************************************************
subroutine convolutionHyperCdriftDam(iRows,iCols,d,dDintegr,dDQ,dThydro,a2dVolumIn, &
			a2dVolumInNet,dTot,dTotIa,dVolTot,dTotHypodFlow,dTotDeepFlow,dTotEsf)

!Esegue la convoluzione del Drift continuo completamente distribuito

implicit none ! Nessuna dichiarazione implicita di variabili
include 'DeclarationH.f90' ! Dichiarazioni per variabili comuni

integer*4 iCountHyper !Contatore tempo di integrazione
integer*4 iTmax !numero step integrazione fine
integer*4 iRows,iCols !Dimensioni dello squadro dem
integer*4 iTtemp !Contatore interno alla Subroutine
real*8 dThydro !Istante temporale di riferimento per l'idrogramma
real*8 iCountT !Contatore tempo simulazione
real*8 dDintegr !Dt del routing
real*8 dDQ !Dt della portata
real*8 d !Dt delle variabili meteorologiche
real*8 a2dQ(iRows,iCols) !Portata istantanea su ogni cella
real*8 a2dVolumIn,a2dVolumInNet,dVolTot,dTotIa,dTot !Input di HortonMatrixHcdrift
real*8 dTotHypodFlow,dTotDeepFlow,Kalpha_zero_mean,dTotEsf !Input di SubFlowHcdrift

REAL*8 a2dCappaCact(iRows,iCols) !Velocità canale attuale
REAL*8 a2dHydroprec(iRows,iCols)
REAL*8 a2dQtmp(iRows,iCols) !Portata in uscita da SurfaceRouting
integer h,mm,m,i,j,iTq,iS
!Inizializzo la portata
a2dQ=0.0
!Numero di step di integrazione se fDintegr è il dt  di integrazione fine in sec
!e d l'intervallo di campionamento delle variabili meteo
iTmax=int(d/dDintegr)
dErr=0
iTtemp=int(dThydro)


!Se c'è il deep flow lo inizializzo a 0 perchè devo sommare per iTmax steps
IF(iFlagDeepFlow.eq.1)THEN
	a2dVLoss=0.0
ENDIF
	!******************************************************************
	!Blocco di Integrazione delle routine con passo d'integrazione fine
	!******************************************************************

	iTq=1
    dVtot2=SUM(SUM(a2dV,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	DO iCountHyper=1,iTmax
	
!-----------------------------------------------------------------------------------
!	Infiltration/Runoff computation by using Modified Horton method
!-----------------------------------------------------------------------------------   
		CALL HortonMatrixHcdrift(iRows,iCols,d,dDintegr,a2dVolumIn,a2dVolumInNet,dTot,dTotIa,dVolTot)
!-----------------------------------------------------------------------------------
!	Hypodermic Flow
!-----------------------------------------------------------------------------------  
		a2dEsf=0.0
		CALL SubFlowHcdriftDam(iRows,iCols,d,dDintegr,dTotHypodFlow,dTotDeepFlow,dTotEsf)
!-----------------------------------------------------------------------------------
!	Surface riuting computation by tank model
!----------------------------------------------------------------------------------- 
		a2dQtmp=0.0
		CALL SurfaceRoutingDam(iRows,iCols,dDintegr,dThydro,a2dQtmp) 
	

		!Portata nelle Celle Canale in mm/s
		IF(iFlagVcVar.eq.1)THEN !cCappaC variabile
			a2dCappaCact=0.0
			!Uso la pendenza
			
			WHERE(a2iChoice.eq.1.AND.a2dCTime.GT.0.0)
				!a2dCappaCact=0.1+dCappaC*(dtan(a2dAlpha)**0.5)*a2dHydroprec**dBc
				a2dQ=a2dQ+a2dQtmp
			ENDWHERE
			
		ELSE
			WHERE(a2iChoice.eq.1.AND.a2dCTime.GT.0.0)
				a2dQ=a2dQ+a2dHydroprec*dCappaC/3600	
			ENDWHERE
		ENDIF

		
		!Portata nelle Celle Versante in mm/s
		WHERE(a2iChoice.eq.0.AND.a2dCTime.GT.0.0)
			a2dQ=a2dQ+a2dQtmp
		ENDWHERE
		!Controllo che l'intervallo di integrazione della Portata sia raggiunto 
		IF(real(iTq*dDintegr).eq.dDQ)then
			!Media su fDQ e porto da mm/s (la costante è in 1/h) a m^3/s.
			!a2dHydro è il livello in mm
			WHERE(a2dCTime.GT.0.0)
				a2dQ=a2dQ/(real(iTq)*1000)*a2dAreaCell
			ENDWHERE


			!Calcolo la portata nelle sezioni di controllo, in a2dXYsections ci sono le coordinate
			DO iS=1,iNsec
				i=a2dXYsections(iS,2)
				j=a2dXYsections(iS,1)
				a2dQsections(iS,iTtemp)=a2dQ(i,j)
				
			ENDDO
			!Incremento il contatore del tempo perchè vado con lo stesso passo di t esterno
			iTtemp=iTtemp+1
			!Annullo la matrice della portata e il contatore su dDQ	
			a2dQ=0.0
			iTq=0
		ENDIF	
		iTq=iTq+1
	ENDDO


dHyTot=SUM(SUM(a2dHydro,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
dRoutTot=SUM(SUM(a2dRouting,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
dVtot=SUM(SUM(a2dV,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
!dTotIa=dTotIa+SUM(SUM(a2dRetention,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
dTotIa=SUM(SUM(a2dRetention,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
dcS=SUM(SUM(a2dS,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
dTot=dTot+SUM(SUM(vol_sot,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
!Se voglio cumulare su tutto il periodo metto dTotDeepFlow=dTotDeepFlow+...
dErr=dErr/real(iTmax)
RETURN
END
