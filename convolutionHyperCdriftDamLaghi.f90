!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************


!***********************************************************************************   
!	Subroutine convolutionHyperCdrift
!	Esegue le equazioni di Hcdrift e la convoluzione
!***********************************************************************************
subroutine convolutionHyperCdriftDamLaghi(iRows,iCols,d,dDintegr,dDQ,dThydro,a2dVolumIn, &
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

real*8 a2dVdt(iRows,iCols) !Velocità e Livello per il Dt
real*8 dVcmax !Velocità massima
real*8 dDintegrAct !Dt del routing attuale (variabile)
real*8 dRainMax !I pioggia max [mm/h]
REAL*8 dDtref !Dt minimo di riferimento (s)
REAL*8 dTmax,dDemPassM
REAL*8 dDtrefRatio !Fattore divisorio per limitare il passo di integrazione
REAL*8 dDintegrMin !Passo di integrazione minimo
REAL*8 dRate

INTEGER iFlagVariable; !Se 1 uso step di integrazione variabile se o 0 



integer h,mm,m,i,j,iTq,iS


!Inizializzo la portata
a2dQ=0.0
!Flag dt variabile 1 variabile, 0 costante
iFlagVariable=1
! Se dt variabile ne seguito viene cambiato
dDintegrAct=dDintegrAct
!Intialization of maximum ratio between d/dDintegr
dDtrefRatio=3
IF(iFlagVariable.eq.1)THEN
	!Inizializzo il dt di riferimento
	dDtref=6
	dDintegrMin=25 ! s
	dDemPassM=sqrt(dCelLat*dCelLon) !Passo medio del DEM
	IF(dDemPassM<=150)THEN !m
		dDtref=5 ! s
		dDtrefRatio=3
		dDintegrMin=5 ! s

	ENDIF
	IF(dDemPassM<10)THEN !m
		dDtref=1 ! s
		dDtrefRatio=3
		dDintegrMin=1 ! s

	ENDIF
	IF(dDemPassM>1000)THEN !m
		dDtref=60 ! s
		dDtrefRatio=2
		dDintegrMin=600 ! s

	ENDIF
	!Inizializzo la velocità sui canali
	a2dVdt=0.0
	!***********************************************************************************************
	!Valutazione del passo di integrazione variabile
	! 
	WHERE (a2iChoice.eq.1.AND.a2dDem.GT.0.0)		
		a2dCappaCact=a2dCappaC*(dtan(a2dBeta)**0.5)*a2dHydro**dBc
		!a2dVdt=a2dHydro*a2dCappaCact/(1000*3600)*dCelLat*dCelLon !m^3/s
		a2dVdt=a2dHydro*a2dCappaCact/(1000*3600)*a2dAreaCell !m^3/s
	ENDWHERE
	WHERE (a2iChoice.eq.0.AND.a2dDem.GT.0.0)
		a2dVdt=a2dHydro*a2dCappaV/(1000*3600)*a2dAreaCell !m^3/s
	ENDWHERE
	!Stima velocità
	WHERE (a2dHydro.gt.0.AND.a2dDem.GT.0.0) !Controllo che il livello sia >0 e calcolo la Vel.
		a2dVdt=a2dVdt/(a2dHydro/1000*dDemPassM) !m/s
	ELSEWHERE
		a2dVdt=0
	ENDWHERE
	dRainMax=maxval(maxval(a2dRain,DIM = 1),DIM = 1)*3600.0/d !Intensità massima
	dVcmax=maxval(maxval(a2dVdt,DIM = 1),DIM = 1)
	IF(dVcmax.le.0.1)dVcmax=0.1
	dDintegrAct=dDemPassM/dVcmax*0.6 !Prima stima passo temporale

	!write(*,*)'Rain max: ',dRainMax
	IF(dDintegrAct.gt.d/dDtrefRatio)dDintegrAct=d/dDtrefRatio !Se più grande di d/3
	IF(dRainMax.gt.1.and.dDintegrAct.gt.dDintegr)THEN  !Se piove inizio a integrare fine
		dDintegrAct=(dDintegr)*dexp(-dRainMax/3.0)+dDintegr
	ENDIF

	IF(dDintegrAct.lt.dDintegrMin)dDintegrAct=dDintegrMin !Se troppo piccolo impongo soglia minima

	!write(*,*)'Intevallo di integrazione: ',dDintegrAct,'s, Vel max=',dVcmax,'m/s'

	dDintegrAct=d/int(d/dDintegrAct)
	dDintegrAct=int(dDintegrAct/dDtref)*dDtref

	!Controllo di non avere aumenti troppo repentine del Dt
	if(dDintegrAct.gt.dDintegrPrec+dDtref)dDintegrAct=dDintegrPrec+dDtref

	!Aggiornamento della variabile a2dRouting:l'intensità non deve cambiare da un passo all'altro 
	a2dRouting=a2dRouting/dDintegrPrec*dDintegrAct

	dDintegrPrec=dDintegrAct

!***********************************************************************************************

ENDIF


!Numero di step di integrazione se fDintegr è il dt  di integrazione fine in sec
!e d l'intervallo di campionamento delle variabili meteo
!iTmax=int(d/dDintegr)
iTmax=int(d/dDintegrAct)
dTmax=d/dDintegrAct
!Se dTmax non è intero devo fare il complemento a d (in sec) 
IF(dTmax.eq.iTmax)THEN
	dTmax=0
ELSE
	dTmax=(dTmax-real(iTmax))*dDintegrAct
ENDIF

dErr=0
iTtemp=1

!write(*,*)'Intevallo di integrazione finale: ',dDintegrAct

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
		CALL HortonMatrixHcdrift(iRows,iCols,d,dDintegrAct,a2dVolumIn,a2dVolumInNet,dTot,dTotIa,dVolTot)
!-----------------------------------------------------------------------------------
!	Hypodermic Flow
!-----------------------------------------------------------------------------------  
		a2dEsf=0.0
		CALL SubFlowHcdriftDam(iRows,iCols,d,dDintegrAct,dTotHypodFlow,dTotDeepFlow,dTotEsf)
!-----------------------------------------------------------------------------------
!	Surface riuting computation by tank model
!----------------------------------------------------------------------------------- 
		a2dQtmp=0.0
		CALL SurfaceRoutingDamLaghi(iRows,iCols,dDintegrAct,dThydro,a2dQtmp,d) 
	

		!Portata nelle Celle Canale in mm/s
		IF(iFlagVcVar.eq.1)THEN !cCappaC variabile
			a2dCappaCact=0.0
			!Uso la pendenza
			WHERE(a2iChoice.eq.1.AND.a2dCTime.GT.0.0)				
				a2dQ=a2dQ+a2dQtmp
			ENDWHERE
			
		ELSE
			WHERE(a2iChoice.eq.1.AND.a2dCTime.GT.0.0)
				a2dQ=a2dQ+a2dHydroprec*a2dCappaC/3600	
			ENDWHERE
		ENDIF

		
		!Portata nelle Celle Versante in mm/s
		WHERE(a2iChoice.eq.0.AND.a2dCTime.GT.0.0)
			a2dQ=a2dQ+a2dQtmp
		ENDWHERE
		!Controllo che l'intervallo di integrazione della Portata sia raggiunto 
		IF(real(iTq)*dDintegrAct.ge.dDQ-dTmax*1.001)then !1.001 per motivi di approssimazione numerica
			!Media su fDQ e porto da mm/s (la costante è in 1/h) a m^3/s.
			!a2dHydro è il livello in mm
			WHERE(a2dCTime.GT.0.0)
				a2dQ=a2dQ/(real(iTq)*1000)*a2dAreaCell
				a2dQmap=a2dQ
			ENDWHERE


			!Calcolo la portata nelle sezioni di controllo, in a2dXYsections ci sono le coordinate
			DO iS=1,iNsec
				i=a2dXYsections(iS,2)
				j=a2dXYsections(iS,1)
				a2dQsections(iS,iTtemp)=a2dQ(i,j)
				!Deflusso subsuperificiale
				dRate=dsin(a2dBeta(i,j))
				if(dRate.gt.1)dRate=0.99
				if(dRate.lt.0)dRate=0.1
				if(a2dVwtMax(i,j).eq.0)dRate=1 !Celle dove non ho falda
				
				
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
