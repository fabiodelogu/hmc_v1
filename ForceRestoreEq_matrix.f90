!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!
!	Soubroutines ForceRestoreEq
!	Risoluzione Numerica Eq. Differenziale force restore equation del primo ordine
!	a2dCon il metodo di Runge Kutta
!	Restituisce l'evapotraspirazione 'a2dEvapot'
!
!	a2dTS(j,i)=a2dTS(j,i)+[2*sqrt(3.14*omega)*(Rn-H-LE)/P_it-2*3.14*omega*(a2dTS(j,i)-TD)]*dt
!
!***********************************************************************************

! Changes
! 12/03/2014: Insert dsqrt and dexp instead of sqrt and exp


subroutine ForceRestoreEq_matrix(iRows,iCols,iStart,d,t_periodo,mm_day,h_day,dDayStep,xstring)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer	i,j,jj,iStart
integer ii,ks,imask1

real*8 d,dDayStep,t_periodo,mm_day,h_day

real*8 tstep,t,RRR,Ta24(iRows,iCols),Ta12(iRows,iCols)

real*8 Ta_gc,EF(iRows,iCols), EFg(iRows,iCols)
real*8 lambda(iRows,iCols), rhow(iRows,iCols), ea(iRows,iCols), Rn(iRows,iCols)
real*8 epsilon_a(iRows,iCols), rhoa(iRows,iCols), Rb(iRows,iCols), Stab_psi(iRows,iCols)
real*8 psi, Chn, sigma, cp, ab, esSTAR(iRows,iCols)
real*8 costK_FR1, Trif, pi, epsilon_s
real*8 rho_cost, g_cost
real*8 P_it(iRows,iCols), VV(iRows,iCols),Ta(iRows,iCols),beta(iRows,iCols),Ch(iRows,iCols)
real*8 omega

real*8 k_c(iRows,iCols),k_b(iRows,iCols),lima,limb,ctwp(iRows,iCols),k_c1(iRows,iCols)
real*8 k_b1(iRows,iCols)

real*8 TD(iRows,iCols),es,dTsMean

real*8 Ts_update(iRows,iCols),a2dDiffLST(iRows,iCols)

real*8 H(iRows,iCols),LE(iRows,iCols),G(iRows,iCols)

real*8 Ts_pixel_med(iNCoord),Ta_pixel_med(iNCoord),TD_pixel_med(iNCoord),VV_pixel_med(iNCoord)
real*8 Rn_pixel_med(iNCoord),LE_pixel_med(iNCoord),H_pixel_med(iNCoord),G_pixel_med(iNCoord)

real*8 TDm,LSTm,T12m,T24m,Tam

! DECIDERE SE MANTENERE IL TUTTO PER DIVERSI BACINI...	
real*8 Ts_m(iNumBasins),Ta_m(iNumBasins),H_m(iNumBasins),LE_m(iNumBasins),LE_pm(iNumBasins)
real*8 ET_m(iNumBasins),EF_m(iNumBasins),Rn_m(iNumBasins),G_m(iNumBasins)
real*8 VV_m(iNumBasins),W_m(iNumBasins),beta_m(iNumBasins),Ch_m(iNumBasins)
real*8 rhoa_m(iNumBasins),rhow_m(iNumBasins),lambda_m(iNumBasins),TD_m(iNumBasins)
real*8 ea_m(iNumBasins),esSTAR_m(iNumBasins),UmR_m(iNumBasins),PP_m(iNumBasins)
real*8 P_it_m(iNumBasins)
real*8 Ta12_m(iNumBasins),Ta24_m(iNumBasins)
real*8 dLSTmax

INTEGER*4 iLStr,iPathLenght,iTLength,temp(iRows,iCols),ios,iTstep
INTEGER anno,mese,giorno,ora,flag, minute
CHARACTER*12  xstring,sVariable
CHARACTER*500 sChfile

!-----------------------------------------------------------------------------------		
!	Calcolo delle quantità utili alla determinazione di a2dTS
!-----------------------------------------------------------------------------------
	
!	CONSTANTS USED IN THE PROGRAM	
!	Temperatura di riferimento [a2dK]
!	Pigreco [-]
!	Calore specifico a pressione costante [J/kg/a2dK]
!	delta LST max °/hour
	
parameter (Trif=273.15, pi=3.14, cp=1004.0, dLSTmax=12)
			
		
!	CONSTANTS FOR Ch COMPUTATION
!	Valore di psi
!parameter (RRR=-7.6) !Valore minimo, Chn=0.0005
!parameter (RRR=-4.6) !Valore massimo, Chn=0.01
!parameter (RRR=-5.85)


!	CONSTANTS FOR Rn COMPUTATION
!	Emissività del suolo [-]
!	Costante Stefan-Boltzmann [a2dW/m^2 a2dK]
!	Valore di albedo

parameter (epsilon_s=0.96, sigma=0.00000005576, ab=0.23)
	
!	CONSTANTS FOR dTs/dt COMPUTATION 	
!	Lunghezza del giorno [a2dS]

parameter (omega=1.0/(60*60*24), lima=0.1, limb=0.9)		



IF(xstring(5:6) == '12'.OR.xstring(5:6) == '01'.OR.xstring(5:6) == '02')THEN
	RRR=-7.3
ENDIF
IF(xstring(5:6) == '03'.OR.xstring(5:6) == '04'.OR.xstring(5:6) == '05')THEN
	RRR=-5.8
ENDIF
IF(xstring(5:6) == '06'.OR.xstring(5:6) == '07'.OR.xstring(5:6) == '08'.OR.xstring(5:6) == '09')THEN
	RRR=-4.8
ENDIF
IF(xstring(5:6) == '10'.OR.xstring(5:6) == '11')THEN
	RRR=-5.9
endif
!RRR=-5.9 se CH costante

psi=log(2.0)				
ctwp=0.4*a2dCt


! int(d) = dStep delle forzanti misurate [a2dS] 

! Discretizzazione del tempo di ogni dStep [a2dS], d è diviso 4 step per
! Runge Kutta
!tstep=int(int(d))	   !450.0
tstep=MIN(int(int(d)),900)
write(*,*)'LST tim step [sec]',tstep

	
!-----------------------------------------------------------------------------------------	
!	Calcolo costanti per la relazione a2dCon beta

k_b=(limb-lima)/(a2dCt-ctwp)    !coefficente angolare retta
k_c=lima-(limb-lima)/(a2dCt-ctwp)*ctwp

k_b1=(1-limb)/(1-a2dCt)    !coefficente angolare retta
k_c1=1-k_b1
!-----------------------------------------------------------------------------------------

EF=0.0
H=0.0
LE=0.0
G=0.0
!a2dEvapot=0.0
Stab_psi=0.0
Ch=0.0
Ta12=0.0
Ta24=0.0

lambda=0.0
rhow=1.0
ea=0.0
Rn=0.0
epsilon_a=0.0
esSTAR=0.0
P_it=0.0
VV=0.0
beta=0.0
temp=0

imask1 = 1

Chn=dexp(RRR)! Ch neutro

!*********************************************************************************************
!Tentativo di lettura del Ch
READ( xstring(1:4), '(i4)' )  anno
READ( xstring(5:6), '(i2)' )  mese
READ( xstring(7:8), '(i2)' )  giorno
READ( xstring(9:10), '(i2)' ) ora
READ( xstring(11:12), '(i2)') minute


if(ora.eq.0)THEN
	if(giorno.eq.1)THEN !I file di CH sono riferiti alle 00 del 1 di ogni mese
		sPathCh=sPathLai
		iPathLenght = iLStr(sPathCh)
		sChfile		=sPathCh(1:iPathLenght)//'Chn_'//xstring(1:8)//'0000'
		iTLength = iLStr(sChfile)
		open(22,file=sChfile(1:iTLength),status='old',form='unformatted', &
				  access='direct',recl=iRows*4,iostat=ios,err=225)
		do j=1,iCols
		   read(22,rec=j)(temp(iRows-i+1,j),i=1,iRows)
		end do
		close(22)

		a2dChn=real(temp)/1000
		iFlagCh=1 !pongo = 1 perchè ho il file dei Ch
		write(*,*)'Uso il Chn '//xstring(1:8)
	else
225		continue
		iFlagCh=0 !Se rimane 0 non ho i file esterni di Ch
	endif
endif

!*********************************************************************************************

!IF(ks.gt.2*dDayStep/24+1*dDayStep/24)ks=1           ! dShift pari a 2 ore ho bisogno di un matrice che contenga 3 ore precedenti

ks=int(dShift*dDayStep/24+1)

!a2dPres atmospheric pressure
!a2dTemp Temperatura dell'aria [°C]
!a2dK	Radiazione solare incidente [a2dW/m^2]
!a2dW	Velocità del vento [m/a2dS]
!a2dUm Umidità [%]
		
WHERE (a2dCTime.gt.0.0)
	a2dPres=101.3*((293-0.0065*a2dDem)/293)**5.26	![kPa]
	VV=a2dV/a2dS										!VV Volume invasato in ogni singola cella [%]
ENDWHERE


! Controllo per dati sgarri e assegnazione dei valori precedenti

WHERE (a2dW.LT.0.0.and.a2dCTime.gt.0.0) a2dW=a2dWPrec
WHERE (a2dTemp.LT.-60.0.and.a2dCTime.gt.0.0) a2dTemp=a2dTempPrec
WHERE (a2dUm.LT.0.0.and.a2dCTime.gt.0.0) a2dUm=a2dUmPrec
WHERE (a2dK.LT.0.0.and.a2dCTime.gt.0.0) a2dK=a2dKPrec
		
WHERE (a2dCTime.gt.0.0)
	a2dWPrec=a2dW
	a2dTempPrec=a2dTemp
	a2dUmPrec=a2dUm
	a2dKPrec=a2dK
ENDWHERE

!Inerzia termica [J m^-2 a2dK a2dS^-(1/2)]
CALL thermal_inertia_matrix(iRows,iCols,VV,P_it)	

!-----------------------------------------------------------------------------------------	
!EQUATION INTEGRATION
!-----------------------------------------------------------------------------------------	

! dichiarare Ta nelle sub se necessario

Ta=-9999
WHERE (a2dCTime.gt.0.0)
	Ta=a2dTemp+Trif  ! T aria in [a2dK]
ENDWHERE

!-----------------------------------------------------------------------------------------	
!media sulle finestre di 12 e 24 ore da otimizzare eventualmente subroutine per il calcolo della tdeep... VEDIAMO alla FINE
		

FORALL(i=1:iRows,j=1:iCols,ii=2:int(dDayStep),a2dCTime(i,j).gt.0.0)
	 a3dTemp24(i,j,ii-1)=a3dTemp24(i,j,ii)
ENDFORALL

FORALL(i=1:iRows,j=1:iCols,a2dCTime(i,j).gt.0.0)
	 a3dTemp24(i,j,int(dDayStep))=Ta(i,j)
ENDFORALL


FORALL(i=1:iRows,j=1:iCols,a2dCTime(i,j).gt.0.0)
  Ta24(i,j)=sum(a3dTemp24(i,j,1:int(dDayStep)))/(dDayStep) !media sulle 24 ore precedenti
  Ta12(i,j)=sum(a3dTemp24(i,j,int(dDayStep/2+1):int(dDayStep)))/int(dDayStep/2.) !media sulle 12 ore precedenti
ENDFORALL

!-----------------------------------------------------------------------------------------	


!Initial condition of LSTwhen starting a new run, no LST status is available

FORALL(i=1:iRows,j=1:iCols,ii=2:ks,a2dCTime(i,j).gt.0.0)
	 a3dTMarked(i,j,ii-1)=a3dTMarked(i,j,ii)
ENDFORALL

IF(iStart.eq.0)THEN	
	dTsMean=sum(sum(a2dTs,DIM = 1,MASK=a2dDem.gt.0),DIM=1)/dBasinArea
	IF(dTsMean.eq.0)then					
		WHERE(a2dCTime.gt.0.0)a2dTS=Ta+1.0 !in [a2dK]   all'inizio del run (istante iniziale)
	endif									
ENDIF
				
IF(mm_day.lt.1.and.iFlagStateVar.eq.0)THEN							!mm_day = contatore del giorno che parte da zero
	FORALL(i=1:iRows,j=1:iCols,a2dCTime(i,j).gt.0.0)
		a3dTMarked(i,j,1:ks)=17.3+Trif	    !media su 24 ore del giorno prima
	ENDFORALL
ELSE
	FORALL(i=1:iRows,j=1:iCols,a2dCTime(i,j).gt.0.0)
		a3dTMarked(i,j,ks)=Ta24(i,j)+(Ta12(i,j)-Ta24(i,j))/exp(1.0)  
	ENDFORALL
ENDIF


WHERE(a2dCTime.gt.0.0)
!-----------------------------------------------------------------------------------------	
!	latent heat of vaporization [J/kg]
	lambda=(2.5-2.36*0.001*(a2dTemp))*1000000  
!-----------------------------------------------------------------------------------------
!	densità dell'acqua (row) in kg/m^3
	rhow=1000.0-0.019549*dabs(a2dTemp-3.98)**1.68
!-----------------------------------------------------------------------------------------
!	pressione di vapore (ea) in kPa - Ta [°C] e a2dUm [%] il risultato è in kPa
	ea=a2dUm*0.611*dexp(17.3*a2dTemp/(237.3+a2dTemp))
!-----------------------------------------------------------------------------------------
!	Calcolo della emissività effettiva dell'atmosfera [%]
!	Unità di misura: ea [kPa] *10 = [millibars]
	epsilon_a=0.740+0.0049*ea*10.0
!-----------------------------------------------------------------------------------------
!	densità dell'aria (costante dei gas per l'aria R=0.288) [kg/m^3]
!	a2dPres [kPa] Ta [a2dK]
	rhoa=a2dPres/(Ta*0.288)
!-----------------------------------------------------------------------------------------
ENDWHERE

beta=0.0
WHERE (VV.lt.ctwp.and.a2dCTime.gt.0.0)
	beta=lima
ELSEWHERE(VV.ge.ctwp.and.VV.le.dCt.and.a2dCTime.gt.0.0)
	beta=k_b*VV+k_c
ELSEWHERE(a2dCTime.gt.0.0)
	beta=k_b1*VV+k_c1
ENDWHERE

!-----------------------------------------------------------------------------------------
!	Chiamata della subroutine per il calcolo del numero di Richardson METTERE IL TUTTO IN FORMA MATRICIALE
!-----------------------------------------------------------------------------------------

CALL richardson_matrix(iRows,iCols,Ta,a2dPres,Rb)

!-----------------------------------------------------------------------------------------
!	Calcolo di Ch QUI Rb è MATRICE QUINDI ANCHE Ch E Stab_psi
!-----------------------------------------------------------------------------------------
!	Il valore di Stab_psi è sempre compreso fra 1 e 3 in una formulazioni di questo tipo


WHERE(Rb.le.0.0.and.a2dCTime.gt.0.0)	
	Stab_psi=1+dexp(psi)*(1-dexp(10*Rb))
ELSEWHERE(a2dCTime.gt.0.0)
	Stab_psi=1.0
ENDWHERE
IF(iFlagCh.eq.0.0)THEN
	!Ch se non ho Chn in lettura (es. da Achab)
	Ch=Chn*Stab_psi
ELSE
	!Ch se ho Chn esterno (es. da Achab)
	Ch=a2dChn*Stab_psi
ENDIF

!-----------------------------------------------------------------------------------------
!	Calcolo della radiazione netta [a2dW/m^2] Rn
!-----------------------------------------------------------------------------------------
!	Unità di misura: sigma [a2dW/m^2°K^4], epsilon_a [%], epsilon_g [%], ab [-], Ta [°K], a2dTS [°K]
!					 a2dK [a2dW/m^2], 
WHERE(a2dCTime.gt.0.0)
	Rn=a2dK*(1.0-ab)+sigma*epsilon_a*Ta**4-sigma*epsilon_s*a2dTS**4
ENDWHERE

!-----------------------------------------------------------------------------------------
!	Ciclo do su int(d) dell'algoritmo di Runge - Kutta
!-----------------------------------------------------------------------------------------
TD=0
FORALL(i=1:iRows,j=1:iCols,a2dCTime(i,j).gt.0.0)
	TD(i,j)=a3dTMarked(i,j,ks-dShift*int(dDayStep/24)) 
ENDFORALL


!Check LST e Tdeep
!LSTm=SUM(SUM(a2dTS,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!TDm=SUM(SUM(TD,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!T12m=SUM(SUM(Ta12,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!T24m=SUM(SUM(Ta24,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!889 FORMAT (1000(A6,1x,f9.2))
!write(*,889)'LSTm:',LSTm,' TDm:',TDm,' Ta12m:',T12m,' Ta24m:',T24m,'Tam: ',Tam


Ts_update=0.0
!Versione vecchia
!DO jj = tstep,int(d),tstep 
!   CALL rk4_matrix(iRows,iCols,dble(jj),tstep,TD,P_it,Ch,lambda,beta,Rn,a2dUm,a2dW,Ta,Trif,ea,rhoa,a2dPres,Ts_update)
!ENDDO

DO jj = 1,int(d)/int(tstep)   
	CALL rk4_matrix(iRows,iCols,dble(jj),tstep,TD,P_it,Ch,lambda,beta,Rn,a2dUm,a2dW,Ta,Trif,ea,rhoa,a2dPres,Ts_update)
	!write(*,*)'Step ',jj,' dim:',tstep
ENDDO

!Check on LST Temperature

a2dDiffLST=0.0
WHERE(a2dCTime.gt.0.0)
	a2dDiffLST=dabs(Ts_update-a2dTS) 	
ENDWHERE
WHERE(a2dCTime.gt.0.0)
	WHERE(a2dDiffLST.gt.dLSTmax*d/3600.0) !Wheret dLST is too large pu old LST
		Ts_update=a2dTS	
	ENDWHERE

	WHERE(Ts_update.lt.273.15-70.0) !Wheret LST is too small put old LST
		Ts_update=a2dTS	
	ENDWHERE

	WHERE(Ts_update.gt.273.15+70.0) !70°C as maximum soil temperature
		Ts_update=a2dTS	
	ENDWHERE
	!Last cchek
	WHERE(Ts_update.lt.273.15-70.0.or.Ts_update.gt.273.15+70.0) !Wheret LST is too large  ot too small put standard LST
		Ts_update=273.15+20	
	ENDWHERE
ENDWHERE
!-----------------------------------------------------------------------------------------
!	Calcolo dei flussi e dell'evapotraspirazione
!-----------------------------------------------------------------------------------------		



WHERE(a2dCTime.gt.0.0)
	a2dTS=Ts_update  
	esSTAR=0.611*dexp(17.3*(a2dTS-Trif)/(237.3+a2dTS-Trif))
	H=rhoa*cp*Ch*a2dW*(a2dTS-Ta)
	LE=rhoa*lambda*Ch*a2dW*beta*(esSTAR-ea)/a2dPres*0.622	
	G=H+LE-Rn
ENDWHERE



WHERE (LE.gt.0.0.AND.a2dCTime.gt.0.0.and.a2dEvapot.ge.0) 
	EF=LE/(LE+H)
	!Evaporative Fraction media giornaliera
	EFg=EFg+EF
	a2dEvapot=LE/(rhow*lambda)*1000*d  ![mm]   
ELSEWHERE(a2dCTime.gt.0.0)
	a2dEvapot=0.0  ![mm]
ENDWHERE

!Rescaling because of spatial resolution
!WHERE (a2dCTime.gt.0.0.and.a2iChoice.ge.0)
!	a2dEvapot=a2dEvapot*(0.7*dexp(-dsqrt(a2dAreaCell)/10.0)+0.3)
!ENDWHERE
!Ts media     	  
dLST=SUM(SUM(a2dTS,DIM=1,MASK=a2dDem.GT.0.0))
iStart=1

RETURN
END SUBROUTINE

!*****************************************************************************************





