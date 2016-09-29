!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************



SUBROUTINE SurfaceRoutingDamLaghi(iRows,iCols,dDintegr,dThydro,a2dQtmp,d)
implicit none
include 'DeclarationH.f90' ! Dichiarazioni per variabili comuni

INTEGER iRows,iCols
REAL*8 a2dIntensityIn(iRows,iCols) !Intensità in ingrsso alla cella superficiale
REAL*8 dCappa !Costante svuotamento superficiale temporanea
REAL*8 dDintegr  !Intervallo di integrazione del routing in secondi
REAL*8 dDth !Dt su cui devo integrare il deflusso ipodermico
INTEGER*4 i,j,ii,jj,iii,jjj,ia,ja,oldj,oldi,iD,iP,iIte
integer*4 iTtemp !Contatore interno alla Subroutine
real*8 dThydro !Istante temporale di riferimento per l'idrogramma
REAL*8 a2dCappaCact(iRows,iCols) !actual friction coefficiento on channels
REAL*8 a2dUhact(iRows,iCols) !actual constant coefficiento on hillslopes
REAL*8 dPend
REAL*8 dQril !Portata in corrispondenza dei rilasci
REAL*8 dEPS,Hprec(iNprese),dDH
REAL*8 a2dQtmp(iRows,iCols),a2dQout(iRows,iCols)
REAL*8 dR1,dR2,dRm,a2dHydroprec(iRows,iCols)
REAL*8 dVlake !Volume temporaneo che entra nel lago
REAL*8 a1dQ_sv(iNdam)	!Portata in uscita fagli organi di scarico (scarico laterale..)
REAL*8 a1dQ_rilaTmp(iNcentr) !portata temporanea dei rilasci dall dighe
REAL*8 dH1,dH2,dHrate !Livelli nei 2 rami della confluenza dei fiumi
REAL*8 dH1f,dQt !H1 fittizia
REAL*8 matnan(iRows,iCols),nansum,dRoutPrec,dDHPrec
REAL*8 d !Dt delle variabili meteo
integer*4 iShift,iRank

iTtemp=int(dThydro)
!Pongo dDth=dDintegr, ma potrei differenziare canale e versante	
dDth=dDintegr

!Velocità Massima in cnalae per evitare problemi numerici
dCappaMax=3600/dDth*0.5
!
a2dIntensityIn=0.0
a2dCappaCact=0.0
a2dQout=0	!Portata in uscita da una cella
a2dQtmp=0	!Portata in volume in uscita da una cella
a1dQ_sv=0   !Portata svasata dagli organi di scarico della diga
!Check sul tirante
WHERE (a2dHydro.LT.0.0)
	a2dHydro=0.0000001
ENDWHERE
WHERE (a2dHydro.GT.100000.0)
	a2dHydro=0.0000001
ENDWHERE


a2dHydroprec=a2dHydro
WHERE(a2dIntensity.lt.0.0)
	a2dIntensity=0.0
ENDWHERE
WHERE(a2dRouting.lt.0.0)
	a2dRouting=0.0
ENDWHERE
WHERE (a2dDem.GT.0.0)
	!Calcolo l'input alla cella superficiale dato dall'esfiltrazione più il runoff (mm/h)
	a2dIntensityIn=a2dIntensity+(a2dEsf)*1000.0*3600.0+(1-a2dCoeffResol)*a2dRouting/dDth*3600.0 !CONTROLLA LA CONVERSIONE

ENDWHERE

!Aggiungo le intensità dovute alle turbinature degli impianti Dighe
DO iD=1,iNcentr
    i=a2dXYCen(iD,2)
	j=a2dXYCen(iD,1)
	!Calcolo le intensità e levo il volume dalle dighe corrispondenti
	if(a1dIdroTurbinate(iD,iTtemp+1).ge.0)then !Se ho le turbinate (si presume che i tempi siano relativi all'immissione)
		a2dIntensityIn(i,j)=a2dIntensityIn(i,j)+a1dIdroTurbinate(iD,iTtemp+1)
		a1dDamVolume(a1iDigaCentrale(iD))=a1dDamVolume(a1iDigaCentrale(iD)) &
					-a1dIdroTurbinate(iD,iTtemp+1)*(a2dAreaCell(i,j))/(1000*3600)*dDth !m^3
	else
		!!Metto la Qmax in mm/h se Vdiga>Vmax se no proporzionale alla potenza 
		a1dQ_rilaTmp(iD)=a1dQturbMax(iD)
		if(a1dDamVolume(a1iDigaCentrale(iD)).lt.a1dVdamMax(a1iDigaCentrale(iD)))then
			a1dQ_rilaTmp(iD)=a1dQturbMax(iD)*(a1dDamVolume(a1iDigaCentrale(iD))/a1dVdamMax(a1iDigaCentrale(iD)))**6
		endif
		if(a1dQ_rilaTmp(iD).lt.0.)a1dQ_rilaTmp(iD)=0.0
		!write(*,*)a1dQ_rilaTmp(iD),a1dQturbMax(iD)
		a2dIntensityIn(i,j)=a2dIntensityIn(i,j)+a1dQ_rilaTmp(iD)*1000*3600/(a2dAreaCell(i,j))
		a1dDamVolume(a1iDigaCentrale(iD))=a1dDamVolume(a1iDigaCentrale(iD)) &
					-a1dQ_rilaTmp(iD)*dDth
	endif
	!controllo che il volume sia >0
	if(a1dDamVolume(a1iDigaCentrale(iD)).lt.0.0)a1dDamVolume(a1iDigaCentrale(iD))=0.0
ENDDO

!Aggiungo le intensità dovute alle turbinature degli impianti Prese-Rilasci
DO iD=1,iNril
    i=a2dXYRilascio(iD,2)
	j=a2dXYRilascio(iD,1)
	if(a1dIdroRila(iD,iTtemp+1).lt.0.)a1dIdroRila(iD,iTtemp+1)=0.0
	a2dIntensityIn(i,j)=a2dIntensityIn(i,j)+a1dIdroRila(iD,iTtemp+1)
ENDDO

!**********************************************************************************
!Equazione superficiale per versanti
a2dUhact=0.0
a2dUhact=a2dCappaV

WHERE (a2iChoice.eq.0.AND.a2dUhact.gt.dCappaMax) !numerical check
	a2dUhact=dCappaMax
ENDWHERE
WHERE (a2iChoice.eq.0.AND.a2dDem.GT.0.0)
    !Eulero diretto
	a2dHydro=a2dHydro+a2dIntensityIn*dDth/3600 &
              -a2dHydro*a2dUhact*dDth/3600
	a2dQout=a2dHydroprec*a2dUhact*dDth/3600
ENDWHERE
!matnan = float(ISNAN (a2dHydro))
WHERE(a2dHydro.lt.0.0)
	matnan=1.0
ENDWHERE
!nansum=SUM(SUM(matnan,DIM=1,MASK=a2dDem.GT.0.0))
!nansum=MINVAL(MINVAL(a2dHydro,DIM=1,MASK=a2dDem.GT.0.0))


!**********************************************************************************
!Equazione superficiale per canali con Cost funzione di H o costante
IF(iFlagVcVar.eq.1)THEN !Se falso Costante come nei versanti (in disuso)
	WHERE (a2iChoice.eq.1.AND.a2dDem.GT.0.0)
		
		a2dCappaCact=a2dCappaC*(dtan(a2dBeta)**0.5)*a2dHydro**dBc
		WHERE (a2dCappaCact.gt.dCappaMax)
			a2dCappaCact=dCappaMax
		ENDWHERE
		a2dQout=a2dHydro*a2dCappaCact*dDth/3600
		!Equazione serbatoio superficiale (Input: runoff(ha anche il routing) e esfiltrazione)
		!Eulero diretto
		a2dHydro=a2dHydro+a2dIntensityIn*dDth/3600  &
				-a2dHydro*a2dCappaCact*dDth/3600
	ENDWHERE

	!***Metodo Trapezi 1 Iterazione
		
	WHERE (a2iChoice.eq.1.AND.a2dDem.GT.0.0)
		a2dQout=a2dCappaC*(dtan(a2dBeta)**0.5)*(0.5*a2dHydroprec**(1+dBc)+0.5*a2dHydro**(1+dBc))*dDth/3600
		!a2dQout=a2dCappaCact*(dtan(a2dBeta)**0.5)*(0.5*a2dHydroprec**(1+dBc)+0.5*a2dHydro**(1+dBc))*dDth/3600

		WHERE (a2dQout.gt.(a2dHydroprec+a2dIntensityIn*dDth/3600)*0.7)
			a2dQout=(a2dHydroprec+a2dIntensityIn*dDth/3600)*0.7
		ENDWHERE
		a2dHydro=a2dHydroprec+a2dIntensityIn*dDth/3600  &
				-a2dQout
	ENDWHERE
	!***Fine metodo Trapezi

ELSE
	
ENDIF


a2dRouting=0.0 !Inizializzo con 0, lo butto in Horton (CONTROLLARE)

!Aggiorno il volume dei laghi a monte delle dighe
DO iD=1,iNdam
    i=a2dXYDam(iD,2)
	j=a2dXYDam(iD,1)
	!L'eventuale aggiornamento di livello lo metto nella diga
	IF(a1dCodeDam(iD).gt.0)THEN 
		!Lago distribuito
		dVlake=SUM(SUM(a2dIntensityIn,DIM=1,MASK=a2iChoice.EQ.a1dCodeDam(iD)))
		a1dDamVolume(iD)=a1dDamVolume(iD)+dVlake*dDth/(3600*1000)*(a2dAreaCell(i,j))!in m^3

	ELSE
		!Lago puntiforme
		a1dDamVolume(iD)=a1dDamVolume(iD)+a2dQout(i,j)/1000*(a2dAreaCell(i,j))		
	ENDIF
	a2dQout(i,j)=0.0
ENDDO

!Aggiorno il volume dei laghi senza dighe
DO iD=1,iNlake
	i=a2dXYLake(iD,2)
	j=a2dXYLake(iD,1)
	IF(a1dCodeLake(iD).gt.0)THEN 
		!Lago distribuito
		dVlake=SUM(SUM(a2dIntensityIn,DIM=1,MASK=a2iChoice.EQ.a1dCodeLake(iD)))
		a1dVlake(iD)=a1dVlake(iD)+dVlake*dDth/(3600*1000)*(a2dAreaCell(i,j)) !in m^3
	ELSE
		!Lago puntiforme
		a1dVlake(iD)=a1dVlake(iD)+a2dQout(i,j)/1000*(a2dAreaCell(i,j))
	ENDIF
END DO

WHERE (a2dHydro.LT.0.0)
	a2dHydro=0.0000001
ENDWHERE

!Calcolo la porzione di acqua che va nella cella successiva
!Essa verrà sommata alla pioggia nella subroutine di Horton
!l'istante successivo
DO i=1,iRows
    DO j=1,iCols
   
		IF (a2dDem(i,j).gt.0)then

			 
			!il contatore punta alla cella successiva(controllare se vale per 
			!l'ultima cella)	  

	        ii=int((a2iPun(i,j)-1)/3)-1
			jj=a2iPun(i,j)-5-3*ii
	        iii=i+ii
			jjj=j+jj

			!Indici non accettabili
	        IF(iii.ge.1.and.jjj.ge.1)then

				!integrazione del routing in mm/passo_integrazione_del_routing
				!L'acqua viene mandata nella cella successiva e utilizzata nella Subrotine
				!Horton
				dRm=a2dQout(i,j) !Trapezi
				a2dRouting(iii,jjj)=a2dRouting(iii,jjj)+dRm  ![mm]
				a2dQtmp(i,j)=dRm/dDth			

			endif
			

		endif

		
  444	continue	 
	
		  
	END DO
END DO

!!Tolgo dal routing la turbinata
DO iP=1,iNprese
	dPend=1
	i=a2dXYPresa(iP,2)
	j=a2dXYPresa(iP,1)
	!transform in in mm
	dDH=a1dIdroPrese(iP,iTtemp+1)/(a2dAreaCell(i,j))*dDth*1000 !In mm
	!write(*,*)'Presa in mm',dDH
	if(1.eq.1)then
		ii=int((a2iPun(i,j)-1)/3)-1
		jj=a2iPun(i,j)-5-3*ii
		iii=i+ii
		jjj=j+jj
		dRoutPrec=0.0
		!Indici non accettabili
		IF(iii.ge.1.and.jjj.ge.1)then
			dRoutPrec=a2dRouting(iii,jjj)			
			!Lascio il 2% della portata in alveo
			if(dDH.gt.0.98*dRoutPrec)then
				dDHPrec=(0.98*dRoutPrec)/dDH
				dDH=0.98*dRoutPrec				
				!Cambio i rilasci
				if(iNprese.eq.iNril)then
					iShift=nint(a1dTcorrprese(iP)*60/d)
					iRank = SIZE (a1dIdroRila,dim=2)
					if(iTtemp+iShift.lt.iRank)then
						a1dIdroRila(iP,iTtemp+iShift)=dDHPrec*a1dIdroRila(iP,iTtemp+iShift)
					endif
				endif
			endif
			
			a2dRouting(iii,jjj)=a2dRouting(iii,jjj)-dDH

			!Deflusso minimo vitale
			dRoutPrec=a2dRouting(iii,jjj)*a2dAreaCell(i,j)/(dDth*1000) !m3/s
			if(dRoutPrec.lt.a1dQminEco(iP).and.dDH.gt.0.0)then
				
				dDHPrec=dDH
				dDH=dDH+a2dRouting(iii,jjj)-a1dQminEco(iP)/(a2dAreaCell(i,j))*dDth*1000 !In mm
				if(dDH.lt.0.0)dDH=0.0
				!Cambio i rilasci
				if(iNprese.eq.iNril)then
					iShift=nint(a1dTcorrprese(iP)*60/d)
					iRank = SIZE (a1dIdroRila,dim=2)
					if(iTtemp+iShift.lt.iRank)then
						a1dIdroRila(iP,iTtemp+iShift)=dDH/dDHPrec*a1dIdroRila(iP,iTtemp+iShift)
					endif
				endif
			endif
			if(a2dRouting(iii,jjj).lt.0.0)then
				a2dRouting(iii,jjj)=0.0
			endif
		ENDIF
	endif
ENDDO

!Calcolo contributi dell'unione tra due fiumi con risalita
DO iD=1,iNjoin
	i=a2dXYMain(iD,2)
	j=a2dXYMain(iD,1)
	dH2=a2dHydro(i,j)
	ia=a2dXYImm(iD,2)
	ja=a2dXYImm(iD,1)
	dH1=a2dHydro(ia,ja)
	if(dH1.lt.dH2.and.dH2.gt.0)then
		dQt=a2dQtmp(ia,ja)/1000*a2dAreaCell(i,j) !m3/s
		dH1f=dsqrt(dH1**2+1000*1000*2*(dQt**2)/(a2dAreaCell(i,j)*9.8*dH1/1000)) !Derivato dall'equazione delle spinte
		
		!write(*,*)dH1f,dH2,dH1,dH1f/dH2
		if(dH1f/dH2.lt.a1dThreshLiv(iD))then
			ii=a2dXYOut(iD,2)
			jj=a2dXYOut(iD,1)
			dHrate=dH1f/dH2
			if(dHrate.lt.0.5)dHrate=0.5
			
			!Master
			a2dRouting(ii,jj)=a2dRouting(ii,jj)+a2dRouting(i,j)*(1-dHrate) !Routing dove immetto la derivazione del Master
			a2dRouting(i,j)=a2dRouting(i,j)*dHrate !Routing che rimane in alveo
			!Immissario
			a2dRouting(ii,jj)=a2dRouting(ii,jj)+a2dRouting(ia,ja)*(1-dHrate) !Routing dove immetto la derivazione dell'immissario
			a2dRouting(ia,ja)=a2dRouting(ia,ja)*dHrate !Routing che rimane in alveo
		endif
	endif

ENDDO
!!Svaso dagli organi di scarico (scarichi laterali...)
CALL SvasoDamVolume(a1dDamVolume,a1dVdamMax,iNdam,dDintegr,a1dQ_sv)
!!Portata a valle della diga
DO iD=1,iNdam
    i=a2dXYDam(iD,2)
	j=a2dXYDam(iD,1)
	ii=int((a2iPun(i,j)-1)/3)-1
	jj=a2iPun(i,j)-5-3*ii
	iii=i+ii
	jjj=j+jj
	a2dRouting(iii,jjj)=a2dRouting(iii,jjj)+a1dQ_sv(iD) !Uscita da organi di scarico subito a valle della diga
	!Equazione riempimento dei laghi
	WHERE (a2iChoice.eq.a1dCodeDam(iD).AND.a2dDem.GT.0.0)
		!Converto il volume in un livello medio
		a2dHydro=a1dDamVolume(iD)/(a1dNumCodeDam(iD)*a2dAreaCell(i,j))*1000 !in mm               
	ENDWHERE
ENDDO


return
END
