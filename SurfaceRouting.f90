


SUBROUTINE SurfaceRouting(iRows,iCols,dDintegr,a2dQtmp)
implicit none
include 'DeclarationH.f90' ! Dichiarazioni per variabili comuni

INTEGER iRows,iCols
REAL*8 a2dIntensityIn(iRows,iCols) !Intensità in ingrsso alla cella superficiale
REAL*8 dCappa !Costante svuotamento superficiale temporanea
REAL*8 dDintegr  !Intervallo di integrazione del routing in secondi
REAL*8 dDth !Dt su cui devo integrare il deflusso ipodermico
INTEGER*4 i,j,ii,jj,iii,jjj,ia,ja,oldj,oldi,iIt

REAL*8 a2dCappaCact(iRows,iCols) !Velocità canale attuale
REAL*8 dCappaCfunction !Funzione di calcolo della Costante Canale variabile
REAL*8 dPend
REAL*8 a2dQtmp(iRows,iCols),a2dQout(iRows,iCols)
REAL*8 a2dHydroprec(iRows,iCols),dHv(iRows,iCols),err(iRows,iCols)

REAL*8 dR1,dR2,dRm,dCappa2
!Pongo dDth=dDintegr, ma potrei differenziare canale e versante	
dDth=dDintegr

!Velocità Massima in cnalae per evitare problemi numerici
dCappaMax=3600/dDth*0.7

a2dIntensityIn=0.0
a2dCappaCact=0.0
a2dQout=0	!Portata in uscita da una cella
a2dQtmp=0	!Portata in volume in uscita da una cella
!Check sul tirante
WHERE (a2dHydro.LT.0.0.OR.a2dDem.LT.0.0)
	a2dHydro=0.00000001
ENDWHERE


a2dHydroprec=a2dHydro

WHERE (a2dDem.GT.0.0)
	!Calcolo l'input alla cella superficiale dato dall'esfiltrazione più il runoff (mm/h)
	a2dIntensityIn=a2dIntensity+(a2dEsf)*1000*3600 !CONTROLLA LA CONVERSIONE
ENDWHERE


!Equazione superficiale per versanti
WHERE (a2iChoice.eq.0.AND.a2dDem.GT.0.0)
	!Equazione serbatoio superficiale (Input: runoff(ha anche il routing) e esfiltrazione)
	!a2dHydro=a2dIntensityIn/(dCappaV*1)*(1-dexp(-dCappaV*dDth/3600))+ &
    !            a2dHydro*dexp(-dCappaV*dDth/3600)
	!Eulero
	!a2dHydro=a2dIntensityIn*dDth/3600+ &
    !           a2dHydro*dexp(-dCappaV*dDth/3600)
	!Eulero diretto
	a2dHydro=a2dHydro+a2dIntensityIn*dDth/3600 &
              -a2dHydro*dCappaV*dDth/3600
	a2dQout=a2dHydroprec*dCappaV*dDth/3600
ENDWHERE

!Equazione superficiale per canali con Cost funzione di H o costante
IF(iFlagVcVar.eq.1)THEN !Se falso Costante come nei versanti (in disuso)
	
	WHERE (a2iChoice.eq.1.AND.a2dDem.GT.0.0)
	
		a2dCappaCact=0.01+dCappaC*(dtan(a2dBeta)**0.5)*a2dHydro**dBc
		WHERE (a2dCappaCact.gt.dCappaMax)
			a2dCappaCact=dCappaMax
		ENDWHERE
		a2dQout=a2dHydro*a2dCappaCact*dDth/3600
		!Equazione serbatoio superficiale (Input: runoff(ha anche il routing) e esfiltrazione)
		!Eulero diretto
		a2dHydro=a2dHydro+a2dIntensityIn*dDth/3600  &
				-a2dHydro*a2dCappaCact*dDth/3600
		
		!a2dQout=a2dHydro*a2dCappaCact*dDth/3600

	ENDWHERE
	!***Metodo Eulero implicito
	DO iIt=1,0
		WHERE (a2iChoice.eq.1.AND.a2dDem.GT.0.0)
			a2dCappaCact=0.1+dCappaC*(dtan(a2dBeta)**0.5)*a2dHydro**dBc
			WHERE (a2dCappaCact.gt.dCappaMax)
				a2dCappaCact=dCappaMax
			ENDWHERE
			a2dHydro=(2*a2dHydroprec+a2dIntensityIn*dDth/3600  &
					-a2dHydro*a2dCappaCact*dDth/3600)*0.5
		ENDWHERE
	ENDDO


		!***Metodo Trapezi 1 Iterazione
		
	WHERE (a2iChoice.eq.1.AND.a2dDem.GT.0.0)
		a2dQout=dCappaC*(dtan(a2dBeta)**0.5)*(0.5*a2dHydroprec**(1+dBc)+0.5*a2dHydro**(1+dBc))*dDth/3600
		WHERE (a2dQout.gt.(a2dHydroprec+a2dIntensityIn*dDth/3600)*0.7)
			a2dQout=(a2dHydroprec+a2dIntensityIn*dDth/3600)*0.7
		ENDWHERE

		a2dHydro=a2dHydroprec+a2dIntensityIn*dDth/3600  &
				-a2dQout
	ENDWHERE

		!***Fine metodo Trapezi

	
ELSE
	!Nel caso Vc costante
	
ENDIF

a2dRouting=0.0 !Inizializzo con 0, lo butto in Horton (CONTROLLARE)
!a2dRoutingV2C=0 !DA LEVARE


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
				!dRm=(dCappa2/3600)*dDth*a2dHydro(i,j)
				!dRm=(dCappa/3600)*dDth*a2dHydroprec(i,j)
				!dRm=((dCappa/3600)*dDth*a2dHydroprec(i,j)+(dCappa2/3600)*dDth*a2dHydro(i,j))*0.5
				!dRm=a2dHydroprec(i,j)*(1-dexp(-dCappa*dDth/3600))+ &
				!		a2dIntensityIn(i,j)*dDth/3600-a2dIntensityIn(i,j)/(dCappa*1)*(1-dexp(-dCappa*dDth/3600))

				!dRm=a2dHydroprec(i,j)*(1-dexp(-dCappa*dDth/3600)) !Eulero diretto
				
				dRm=a2dQout(i,j) !Trapezi
				a2dRouting(iii,jjj)=a2dRouting(iii,jjj)+dRm  ![mm]
				a2dQtmp(i,j)=dRm/dDth
				!a2dRoutingV2C(i,j)=dRm
				
			endif
			

		endif

		
  444	continue	 
	
		  
	END DO
END DO

!a2dHydro=a2dHydro+err
!WHERE (a2dHydro.LT.0.0.OR.a2dDem.LT.0.0)
!	a2dHydro=0.000001
!ENDWHERE
!err=a2dHydroprec+dHv-a2dRoutingV2C-a2dHydro

return
END
