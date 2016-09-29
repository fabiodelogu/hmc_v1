


SUBROUTINE SurfaceRoutingDam(iRows,iCols,dDintegr,dThydro,a2dQtmp)
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
REAL*8 a2dCappaCact(iRows,iCols) !Velocità canale attuale
REAL*8 dCappaCfunction !Funzione di calcolo della Costante Canale variabile
REAL*8 dPend
REAL*8 dQril !Portata in corrispondenza dei rilasci
REAL*8 dEPS,Hprec(iNprese),dDH,dNcell
REAL*8 a2dQtmp(iRows,iCols),a2dQout(iRows,iCols)
REAL*8 dR1,dR2,dRm,a2dHydroprec(iRows,iCols)

iTtemp=int(dThydro)
!Pongo dDth=dDintegr, ma potrei differenziare canale e versante	
dDth=dDintegr

!Velocità Massima in cnalae per evitare problemi numerici
dCappaMax=3600/dDth*0.7

a2dIntensityIn=0.0
a2dCappaCact=0.0
a2dQout=0	!Portata in uscita da una cella
a2dQtmp=0	!Portata in volume in uscita da una cella

!Check sul tirante
WHERE (a2dHydro.LT.0.0)
	a2dHydro=0.0000001
ENDWHERE


a2dHydroprec=a2dHydro

WHERE (a2dDem.GT.0.0)
	!Calcolo l'input alla cella superficiale dato dall'esfiltrazione più il runoff (mm/h)
	a2dIntensityIn=a2dIntensity+(a2dEsf)*1000*3600 !CONTROLLA LA CONVERSIONE
ENDWHERE

!Aggiungo le intensità dovute alle turbinature degli impianti Dighe
DO iD=1,iNcentr
    i=a2dXYCen(iD,2)
	j=a2dXYCen(iD,1)
	if(a1dIdroTurbinate(iD,iTtemp+1).ge.0.0)then !Se ho -9999 non faccio nulla
		a2dIntensityIn(i,j)=a2dIntensityIn(i,j)+a1dIdroTurbinate(iD,iTtemp+1)
	endif
ENDDO

!Aggiungo le intensità dovute alle turbinature degli impianti Prese-Rilasci
DO iD=1,iNril
    i=a2dXYRilascio(iD,2)
	j=a2dXYRilascio(iD,1)
	a2dIntensityIn(i,j)=a2dIntensityIn(i,j)+a1dIdroRila(iD,iTtemp+1)
ENDDO

!Equazione superficiale per versanti
WHERE (a2iChoice.eq.0.AND.a2dDem.GT.0.0)
   !Eulero diretto
	a2dHydro=a2dHydro+a2dIntensityIn*dDth/3600 &
              -a2dHydro*dCappaV*dDth/3600
	a2dQout=a2dHydroprec*dCappaV*dDth/3600
ENDWHERE

!Equazione superficiale per canali con Cost funzione di H o costante
IF(iFlagVcVar.eq.1)THEN !Se falso Costante come nei versanti (in disuso)
	WHERE (a2iChoice.eq.1.AND.a2dDem.GT.0.0)
		
		a2dCappaCact=0.1+dCappaC*(dtan(a2dBeta)**0.5)*a2dHydro**dBc
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
		a2dQout=dCappaC*(dtan(a2dBeta)**0.5)*(0.5*a2dHydroprec**(1+dBc)+0.5*a2dHydro**(1+dBc))*dDth/3600
		WHERE (a2dQout.gt.(a2dHydroprec+a2dIntensityIn*dDth/3600)*0.7)
			a2dQout=(a2dHydroprec+a2dIntensityIn*dDth/3600)*0.7
		ENDWHERE
		a2dHydro=a2dHydroprec+a2dIntensityIn*dDth/3600  &
				-a2dQout
	ENDWHERE

		!***Fine metodo Trapezi

ELSE
	
ENDIF




!do iP=1,1
!	i=a2dXYPresa(iP,2)
!	j=a2dXYPresa(iP,1)
!	write(*,*)a2dHydro(i,j)
!enddo


a2dRouting=0.0 !Inizializzo con 0, lo butto in Horton (CONTROLLARE)


!Pongo a 0 il livello nelle celle corrispondenti alle dighe
DO iD=1,iNdam
    i=a2dXYDam(iD,2)
	j=a2dXYDam(iD,1)
	!L'eventuale aggiornamento di livello lo metto nella diga
	a1dDamVolume(iD)=a1dDamVolume(iD)+a2dQout(i,j)/1000*(a2dAreaCell(i,j))
	a2dQout(i,j)=0.0
ENDDO


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
dNcell=1

DO iP=1,iNprese
	dPend=1
	i=a2dXYPresa(iP,2)
	j=a2dXYPresa(iP,1)
	dDH=a1dIdroPrese(iP,iTtemp+1)/(a2dAreaCell(i,j))*dDth*1000/dNcell
	if(1.eq.1)then
		!Suddivido la portata della presa su dNcell celle
		DO iIte=1,dNcell
			ii=int((a2iPun(i,j)-1)/3)-1
			jj=a2iPun(i,j)-5-3*ii
			iii=i+ii
			jjj=j+jj
			!Indici non accettabili
			IF(iii.ge.1.and.jjj.ge.1)then
				a2dRouting(iii,jjj)=a2dRouting(iii,jjj)-dDH
				if(a2dRouting(iii,jjj).lt.0.0)then
					a2dRouting(iii,jjj)=0.0
				endif
			ENDIF
			i=iii
			j=jjj
		ENDDO
	endif
ENDDO

!!Riempimento volume della diga con acqua a monte
DO iD=1,iNdam
    i=a2dXYDam(iD,2)
	j=a2dXYDam(iD,1)
	a1dDamVolume(iD)=a1dDamVolume(iD)+a2dRouting(i,j)/1000*(a2dAreaCell(i,j))
	a2dRouting(i,j)=0.0
ENDDO


return
END
