!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!------------------------------------------------------------------
!	restituisce: Q_out=idrogramma di uscita dall'opera
!------------------------------------------------------------------

!Funzione che calcola il coefficiente di svaso per serbatoio lineare equivalente,
!In base a livelli e volumi (misurati o stimati)

 SUBROUTINE CoeffDamEquiv(dVdam,dVdamMax,iNdighe)
	include 'DeclarationH.f90'


	integer*4 iNdighe,nn,i,ii,iRank,iFlagD
    integer*4 name_lenght,iMaxC
	real*8 dVdam(iNdighe),dVdamMax(iNdighe),dQtmp
	real*8 a1dH(iNdighe)
	real*8 dTV !Percentuale di volume oltre il quale inizia lo svaso dagli organi di scarico
	real*8 dQsprec !Portata scarico superficiale step precedente

	dTV=0.90
	a1dH=0
	
	do i=1,iNdighe
		dQsprec=0.0
		if(dVdam(i).gt.dVdamMax(i)*dTV)then !Considero il 93% volume totale
			!Check sulle lunghezze (forse da levare)
			if(adL(i).le.0)then
				adL(i)=40 !m
			endif
			!Portata scarico superficiale step precedente
			dQsprec=a1dCoefDighe(i)*(dVdam(i)-dTV*dVdamMax(i))

			!Interpolazione lineare
			iRank = SIZE (a1d_Volume,dim=2)
			!write(*,*)iRank
			iFlagD=0
			do ii=iRank,2,-1

				if(a1d_Volume(i,ii).gt.dVdam(i).and.a1d_Volume(i,ii-1).lt.dVdam(i))then
					a1dH(i)=a1d_Level(i,ii-1)+(a1d_Level(i,ii)-a1d_Level(i,ii-1))* &
					(dVdam(i)-a1d_Volume(i,ii-1))/(a1d_Volume(i,ii)-a1d_Volume(i,ii-1))
					iFlagD=1
				endif

				if(a1d_Volume(i,ii).lt.0.0)then !Non ho la curva invaso-volumi
					iFlagD=2
				endif
				
			enddo
			
			if(iFlagD.eq.0)then !Sono fuori dalla curva invaso volume
				iMaxC=MAXLOC(a1d_Level(i,:),1) 
				write(*,*)iMaxC
				a1dH(i)=a1d_Level(i,iMaxC)+(a1d_Level(i,iMaxC)-a1d_Level(i,iMaxC-1))* &
					(dVdam(i)-a1d_Volume(i,iMaxC))/(a1d_Volume(i,iMaxC)-a1d_Volume(i,iMaxC-1))
			endif


			if(iFlagD.eq.2)then !Non ho la curva invaso volume uso relazione lineare
				a1dH(i)=a1dHmax(i)*dVdam(i)/dVdamMax(i)
			endif
			!write(*,*)'Altezza',a1dH(i),' flag',iFlagD
			if(a1dH(i).gt.a1dHmax(i)-3)then !Controllo che H sia > Hmax - 3 m e inizio a svasare acqua
				dQtmp=0.385*adL(i)*((2*9.81)**0.5)*(a1dH(i)-(a1dHmax(i)-3))**1.5 !In m3/s
				!write(*,*)'Q svaso prima',dQtmp
				if(dQtmp.gt.a1dQ_sLC(i))dQtmp=a1dQ_sLC(i) !Controllo che non sia > Q massima smaltibile
				if(a1dH(i).gt.a1dHmax(i))then !Controllo che H sia > Hmax 
					dQtmp=(dVdam(i)-dVdamMax(i))/3600 !Butto via tutto il volume 					
					!Se è minore dello scarico massimo, utilizzo lo scarico massimo
					dQtmp=max(dQtmp,a1dQ_sLC(i));
				endif
				
			else
				dQtmp=0
			endif
						
			!Coefficiente medio del corrispondente serbatoio lineare
			if(dQtmp.gt.0)then
				a1dCoefDighe(i)=(0.5*dQtmp+0.5*dQsprec)/(dVdam(i)-dTV*dVdamMax(i)) !In 1/s
			else
				a1dCoefDighe(i)=0.000 ! In 1/s	
			endif
			!write(*,*)'Volume ',dVdam(i)
			!write(*,*)'diga ',i,'Portata ',dQtmp,'H ',a1dH(i),'H max ',a1dHmax(i),'Coeff. ',a1dCoefDighe(i)
		endif
	enddo

	RETURN
	END

