!***********************************************************************************
!	 Hyper C-DRiFt:
!
!	- Modified Horton method for infiltration 
!		(S. Gabellani, F. Silvestro, R. Rudari G. Boni. "General calibration methodology for infiltration scheme in flash flood modeling".
!		 Nat. Hazards Earth Syst. Sci., 8, 1317-1327, 2008);
!	- Runoff Routing;
!	- Hypodermic flow Routing; 
!	- Complete Energy Balance using Force Restore equation to compute soil temperature;
!	- Tdeep filter;
!	- WaterTable  (working...).
!***********************************************************************************
PROGRAM MainContinuum

IMPLICIT NONE 
INCLUDE 'DeclarationH.f90' ! Common variables declaration

!-----------------------------------------------------------------------------------
!	Dichiarazione dei TARGET dei vettori o delle matrici dinamiche comuni
!-----------------------------------------------------------------------------------
!	Matrici Tempoinvarianti
REAL*8,  ALLOCATABLE, TARGET :: a2dACTime(:,:),a2dADem(:,:),a2dAS(:,:),a2dACurNum(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dACostF(:,:),a2dACostK(:,:),a2dACostF1(:,:),a2dACostChFix(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dAAlpha(:,:),a2dAKAlpha(:,:),a2dAWaterTable(:,:)	
REAL*8,  ALLOCATABLE, TARGET :: a2dALimSupMatrix(:,:),a2dAKAlphaZero(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dAAreeWT(:,:),a2dADiffMaxMatrix(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dAC1(:,:),a2dAF2(:,:)
INTEGER, ALLOCATABLE, TARGET :: a2iAMask(:,:),a2dACon(:,:),a2iAPun(:,:),a2iAChoice(:,:)
INTEGER,  ALLOCATABLE, TARGET :: a2iAYLat(:,:), a2iAXLon(:,:)

! Matrici di Prova per Sublflow
REAL*8,  ALLOCATABLE, TARGET :: a2dAVwt(:,:),a2dADeepFlow(:,:),a2dAVwtMax(:,:),a2dABeta(:,:)

!	Matrici Tempovarianti
REAL*8,  ALLOCATABLE, TARGET :: Arain(:,:),Acumrain(:,:),Acumrain_prec(:,:)
REAL*8,  ALLOCATABLE, TARGET :: ATemp(:,:),AK(:,:),AW(:,:),AUm(:,:),APres(:,:)
REAL*8,  ALLOCATABLE, TARGET :: Asat_m(:,:),Aesf(:,:),Avol_sot(:,:),AV(:,:),AVloss(:,:)
REAL*8,  ALLOCATABLE, TARGET :: Aevapot(:,:),Aintensity(:,:),Aritenzione(:,:)
REAL*8,  ALLOCATABLE, TARGET :: AEF_mat(:,:),AW_prec(:,:),ATT_prec(:,:),AUm_prec(:,:),AK_prec(:,:)
REAL*8,  ALLOCATABLE, TARGET :: AT_segnato(:,:,:),ATa_vet24(:,:,:),ATs(:,:)
REAL*8,  ALLOCATABLE, TARGET :: Avloss_aree_wt_t(:,:),Awatertable_update_matrix(:,:)
REAL*8,  ALLOCATABLE, TARGET :: Avolume_zero_matrix(:,:),Avolume_update_matrix(:,:)
REAL*8,  ALLOCATABLE, TARGET :: Avloss_aree_wt_t_cum(:,:),Avloss_cum(:,:)
INTEGER, ALLOCATABLE, TARGET :: Anpixel(:)

REAL*8,  ALLOCATABLE, TARGET :: a2dAHydro(:,:),a2dARouting(:,:),a2dARoutingV2C(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dAQsections(:,:),a2dAQmap(:,:),a2dALai(:,:),a2dAChn(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dAAreaCell(:,:)

INTEGER,  ALLOCATABLE, TARGET :: a2dAXYsections(:,:)

!	Matrici Neve
REAL*8,  ALLOCATABLE, TARGET :: a2dANature(:,:),a2dASnowFall(:,:),a2dAExp(:,:),a2dAAge(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dASWE(:,:),a2dAAlbedo(:,:),a2dARoS0(:,:),a2dARoS(:,:),a2dAMeanDayTemp(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dAMeltingDayCum(:,:),a2dAMelting(:,:),a2dAsMc(:,:),a2dAsMcGlac(:,:)

!-----------------------------------------------------------------------------------
!	Dichiarazioni Locali
!-----------------------------------------------------------------------------------
!	Dichiarazione dei vettori o delle matrici dinamiche del main del programma		
REAL*8 a2dVolumIn(:),a2dVolumInNet(:)
REAL*8 a2dQSotT(:,:),a2dQT(:,:),a2dEvTMat(:,:)
REAL*8 a1dLatP(:),a1dLonP(:),a1dLatT(:),a1dLonT(:),a1dLatK(:),a1dLonK(:),a1dLatW(:),a1dLonW(:),a1dLatUm(:),a1dLonUm(:)
REAL*8 a1dZP(:),a1dZT(:),a1dZK(:),a1dZW(:),a1dZUm(:)


ALLOCATABLE a2dVolumIn,a2dVolumInNet
ALLOCATABLE a2dQT,a2dQSotT,a2dEvTMat
ALLOCATABLE a1dLatP,a1dLonP,a1dLatT,a1dLonT,a1dLatK,a1dLonK,a1dLatW,a1dLonW,a1dLatUm,a1dLonUm
ALLOCATABLE a1dZP,a1dZT,a1dZK,a1dZW,a1dZUm

!	Dichiarazione delle variabili utilizzate del main del programma
INTEGER i,j,n,nn,nb,tt
INTEGER iRows,iCols,iStart,iStepRain,iClass,iLVect
INTEGER*4 iLStr

REAL*8 dSimLength,d,dd,dIst,dStep,dDayStep,t,h,m,mm,dTc
REAL*8 dCTimeMax,dVolumInPrec,dVolumInNetPrec,r
REAL*8 dIniz,dFin,dVolTot,dTot,dTotEsf,dTotIa,dTotHypodFlow,dTotDeepFlow

REAL*8 dETDayMean
REAL*8 rain_m_ist,tot_vloss_cum

REAL*8, ALLOCATABLE :: dSatBac(:),dSatFC(:),dEvTTot(:)

REAL*8, ALLOCATABLE :: matrice_var(:,:),rain2file(:,:)

REAL*8 diff_quantile,tot_vloss_aree_wt_t_cum

INTEGER iFlagInterp,iFlagGrid
CHARACTER*12  xstring,xstring_up,xstring_do,sDataIni,sVariabile
CHARACTER*80 sBuffer

real*8 dDintegr !Dt del routing
real*8 dTstart	!Tempo di inizio della simulazione
real*8 dTo! !Vale 1 per ripresa run e 0 per per run dall'inizio periodo
real*8 dSt1,dSt2,RTC,dTime_spent !Variabili per calcolo dei tempi
real*8 dHourSD,dMonthSD !Hour and Month saving data !Hour saving data



!***********************************************************************************
!	Argomenti da riga di comando
!***********************************************************************************
dKsatRatio=1
!Leggo costante di svuotamento canale
CALL GETARG(1,sBuffer)
read(sBuffer,*)dCappaC
!Leggo costante di svuotamento versante
CALL GETARG(2,sBuffer)
read(sBuffer,*)dCappaV
!Leggo ct, capillarità Horton
CALL GETARG(3,sBuffer)
read(sBuffer,*)dCt
!Leggo cf svuotamento serbatoio Horotn
CALL GETARG(4,sBuffer)
read(sBuffer,*)dCf
!Leggo il nome del bacino
CALL GETARG(5,sBuffer)
sBasin=sBuffer
!Leggo coefficiente di umidità iniziale
CALL GETARG(6,sBuffer)
read(sBuffer,*)dCPI
!Leggo la Hmax della WT
CALL GETARG(7,sBuffer)
read(sBuffer,*)dHbr
!Leggo il rapporto tra Ksaturo vertical e orizzontale Ksh/Ksv
CALL GETARG(8,sBuffer)
read(sBuffer,*)dKsatRatio
!Leggo la massima pendenza per cui ho sottosuolo
CALL GETARG(9,sBuffer)
read(sBuffer,*)dSlopeMax

if(dSlopeMax<40)then
	dSlopeMax=70 !Valore standard
endif


!***********************************************************************************
!***********************************************************************************

!-----------------------------------------------------------------------------------
!   Dichiarazione del nome del bacino idrografico analizzato e del nome del run con
!	cui identificare i files in uscita 
!sBasin='casentino'	
sRunName='_Continuum_'
 
iNameLenght = iLStr(sBasin)
iNameLenghtRun = iLStr(sRunName)

!-----------------------------------------------------------------------------------
!   Reading Information File
!-----------------------------------------------------------------------------------
CALL ReadInfoData(iRows,iCols,iFlagInterp,xstring,dDintegr,sDataIni)
!iFlagVcVar=1 !Se 1 dCappaC variable se no costante METTERE IN INFO!

!	Calcolo del passo temporale in secondi dei dati in ingresso
iStepRain=dHours/dMinutes*60
dd=3600*dHours
d=60*dMinutes
dIst=dx
dStep=d/dIst
dCTimeMax=0
dDayStep=24*60/dMinutes

iNsec=iNumBasins !Numero sezioni lette su file


!-----------------------------------------------------------------------------------
!   Reading Information Snow File
!-----------------------------------------------------------------------------------
CALL ReadInfoS3M(iRows,iCols)

!-----------------------------------------------------------------------------------
!	Allocazione dei TARGET dei vettori o delle matrici comuni
!-----------------------------------------------------------------------------------
!	Matrici Tempoinvarianti
ALLOCATE(a2iAMask(iRows,iCols),a2dACon(iRows,iCols))
ALLOCATE(a2dACTime(iRows,iCols),a2dADem(iRows,iCols),a2iAPun(iRows,iCols),a2dAS(iRows,iCols),a2dACurNum(iRows,iCols))
ALLOCATE(a2iAChoice(iRows,iCols))
ALLOCATE(a2dACostF(iRows,iCols),a2dACostK(iRows,iCols),a2dACostF1(iRows,iCols),a2dACostChFix(iRows,iCols))
ALLOCATE(a2dAAlpha(iRows,iCols),a2dAKAlpha(iRows,iCols),a2dAWaterTable(iRows,iCols))
ALLOCATE(a2dALimSupMatrix(iRows,iCols),a2dAKAlphaZero(iRows,iCols))
ALLOCATE(a2dAAreeWT(iRows,iCols),a2dADiffMaxMatrix(iRows,iCols))
ALLOCATE(a2dAC1(iRows,iCols),a2dAF2(iRows,iCols))

ALLOCATE(a2dAVwt(iRows,iCols),a2dADeepFlow(iRows,iCols),a2dAVwtMax(iRows,iCols),a2dABeta(iRows,iCols))

a2iMask			      => a2iAMask	         
a2dCon				  => a2dACon		         	       
a2dCTime			  => a2dACTime	         
a2dDem				  => a2dADem		         
a2iPun			      => a2iAPun	  
a2dS				  => a2dAS	
a2dCurNum			  => a2dACurNum
a2iChoice			  => a2iAChoice  
a2dCostF			  => a2dACostF	         
a2dCostK			  => a2dACostK	         
a2dCostF1			  => a2dACostF1
a2dCostChFix		  => a2dACostChFix        
a2dAlpha			  => a2dAAlpha                               
a2dKAlpha		   	  => a2dAKAlpha  
a2dWaterTable		  => a2dAWaterTable 
a2dAreeWT			  => a2dAAreeWT

a2dKAlphaZero		  => a2dAKAlphaZero

a2dDiffMaxMatrix      => a2dADiffMaxMatrix	

a2dC1				=> a2dAC1
a2dF2				=> a2dAF2

a2dVwt				=> a2dAVwt

a2dDeepFlow		=> a2dADeepFlow
a2dVwtMax		=> a2dAVwtMax
a2dBeta		=> a2dABeta

a2dBeta=0.0
a2dVwtMax=0.0
a2dVwt=0.0
a2dDeepFlow=0.0
a2iMask=0
!	Matrici Tempovarianti
ALLOCATE(Arain(iRows,iCols),Acumrain(iRows,iCols),Acumrain_prec(iRows,iCols))
ALLOCATE(ATemp(iRows,iCols),AK(iRows,iCols),AW(iRows,iCols),AUm(iRows,iCols),APres(iRows,iCols))
ALLOCATE(Asat_m(iRows,iCols),Aesf(iRows,iCols),Avol_sot(iRows,iCols),AV(iRows,iCols),AVloss(iRows,iCols))
ALLOCATE(Aevapot(iRows,iCols),Aintensity(iRows,iCols),Aritenzione(iRows,iCols))
ALLOCATE(AEF_mat(iRows,iCols),AW_prec(iRows,iCols),ATT_prec(iRows,iCols),AUm_prec(iRows,iCols),AK_prec(iRows,iCols))
ALLOCATE(AT_segnato(iRows,iCols,dShift*INT(60/dMinutes)+1),ATa_vet24(iRows,iCols,INT(dDayStep)), ATs(iRows,iCols))
ALLOCATE(Avloss_aree_wt_t(iRows,iCols))
ALLOCATE(Avolume_zero_matrix(iRows,iCols),Avolume_update_matrix(iRows,iCols),Awatertable_update_matrix(iRows,iCols))
ALLOCATE(Avloss_aree_wt_t_cum(iRows,iCols),Avloss_cum(iRows,iCols))
ALLOCATE(Anpixel(iNCoord))


ALLOCATE(a2dAHydro(iRows,iCols),a2dARouting(iRows,iCols),a2dARoutingV2C(iRows,iCols))
ALLOCATE(a2dAXYsections(iNsec,2),a2dAQmap(iRows,iCols),a2dALai(iRows,iCols),a2dAChn(iRows,iCols))
ALLOCATE(a2dAAreaCell(iRows,iCols))

a2dRain			    => Arain				 
a2dCumRain			=> Acumrain			
a2dCumRainPrec	    => Acumrain_prec                               	                             
a2dTemp			    => ATemp	         
a2dK				=> AK				       
a2dW				=> AW				       
a2dUm				=> AUm			         
a2dPres			    => APres		         
a2dSatM			    => Asat_m		       
a2dEsf				=> Aesf			       
vol_sot			    => Avol_sot	       
a2dV				=> AV	
a2dVLoss			=> AVloss				       
a2dEvapot			=> Aevapot		     
a2dIntensity		=> Aintensity	   
a2dRetention		=> Aritenzione    
a2dEF		    	=> AEF_mat		     
a2dWPrec			=> AW_prec		     
a2dTempPrec			=> ATT_prec	
a2dUmPrec			=> AUm_prec	   
a2dKPrec			=> AK_prec
a3dTMarked		    => AT_segnato	   
a3dTemp24		    => ATa_vet24
a2dTS				=> ATs	  
a2dVLossAreeWTStep  => Avloss_aree_wt_t
!Inizializzo la Tsuolo
a2dTs=0.0

a2dWaterTableUpdateMatrix	=> Awatertable_update_matrix
a2dVolumeZeroMatrix			=> Avolume_zero_matrix
a2dVolumeUpDateMatrix		=> Avolume_update_matrix
a2dVLossAreeWTStepCum		=> Avloss_aree_wt_t_cum
a2dVLossCum					=> Avloss_cum

a2dHydro			=> a2dAHydro
a2dRouting			=> a2dARouting
a2dRoutingV2C		=> a2dARoutingV2C
!Inizializzo le variabili per il routing superficiale
a2dHydro=0
a2dRouting=0.0
a2dRoutingV2C=0.0


a2dXYsections		=> a2dAXYsections
!Inizializzo il vettore con le sezione
a2dXYsections=0.0

a2dQmap		=>	a2dAQmap
!Inizializzo la mappa di portate
a2dQmap=0.0

a2dLai		=> a2dALai
!Inizializzo la mappa di LAI e il flag Lai
a2dLai=0.0
dFlagLai=0 !0 se non ho mappa LAI, 1 se si

!Inizializzo il Ch neutro
a2dChn	=>	a2dAChn
a2dChn=0.0
iFlagCh=0 ! 0 se non ho chb esterno

!Inizializzo la mappa con le aree delle celle
a2dAreaCell	=>	a2dAAreaCell
a2dAreaCell=0.0
IF(dCelLat.gt.0.0.and.dCelLon.gt.0.0)THEN
	a2dAreaCell=(dCelLat*dCelLon)	
ENDIF

!	Matrici Snow
ALLOCATE(a2dANature(iRows,iCols),a2dASnowFall(iRows,iCols),a2dAExp(iRows,iCols),a2dAAge(iRows,iCols))
ALLOCATE(a2dASWE(iRows,iCols),a2dAAlbedo(iRows,iCols),a2dARoS0(iRows,iCols),a2dARoS(iRows,iCols),a2dAMeanDayTemp(iRows,iCols))
ALLOCATE(a2dAMeltingDayCum(iRows,iCols),a2dAMelting(iRows,iCols),a2dAsMc(iRows,iCols),a2dAsMcGlac(iRows,iCols))

a2dNature	=>	a2dANature
a2dSnowFall	=>	a2dASnowFall
a2dExp	=> a2dAExp
a2dAge	=>	a2dAAge	
a2dSWE	=>	a2dASWE
a2dAlbedo	=>	a2dAAlbedo
a2dRoS0	=>	a2dARoS0
a2dRoS	=>	a2dARoS
a2dMeanDayTemp	=>	a2dAMeanDayTemp
a2dMeltingDayCum	=>	a2dAMeltingDayCum
a2dMelting	=>	a2dAMelting
a2dsMc	=>	a2dAsMc
a2dsMcGlac	=>	a2dAsMcGlac

a2dNature=0.0
a2dSnowFall=0.0
a2dExp=0.0
a2dAge=0.0
a2dSWE=0.0
a2dAlbedo=0.0
a2dRoS0=0.0
a2dRoS=0.0
a2dMeanDayTemp=0.0
a2dMeltingDayCum=0.0
a2dMelting=0.0
a2dsMc=0.0
a2dsMcGlac=0.0

!-----------------------------------------------------------------------------------
!   Chiamata alla Subroutine READ_land_data
!   Lettura dei file territoriali del sBasin
!-----------------------------------------------------------------------------------
dTstart=0

CALL ReadLandData(iRows,iCols)

!Change units the surface flow parameters in the one used in the code
if(dCappaV.lt.0.05)then
	dCappaV=dCappaV*3600 !from 1/s to 1/h
	dCappaC=dCappaC*(3600*1000)/(sqrt(dCelLat*dCelLon)*1000**(dBc+1)) !from m^0.5/s to 1/(h*mm^0.5)	
endif

!-----------------------------------------------------------------------------------
!   Generate the conversion from meteo grid and dem grid
!-----------------------------------------------------------------------------------
!Check if meteo and dem grid are different


iFlagGrid=0 !1 the Grids are different
IF((dXDemLat.ne.dXCNLat).or.(dXDemLon.ne.dXCNLon).or.(dDemPasLon.ne.dCNPasLon))THEN 
	ALLOCATE(a2iAYLat(iRows,iCols),a2iAXLon(iRows,iCols))
	
	a2iYLat			      => a2iAYLat	         
	a2iXLon				  => a2iAXLon	
	CALL Create_Meteo_GRID(iRows,iCols,iRowsMeteo,iColsMeteo)
	iFlagGrid=1
ELSE
	iRowsMeteo=iRows
	iColsMeteo=iCols
ENDIF



!-----------------------------------------------------------------------------------
!Scrittura Tc indicativo stimato dell'area di interesse
dTc=IDNINT(0.27*sqrt(0.6*float(iRows)*float(iCols)*dCelLat*dCelLon/1000000)+0.25)
WRITE(*,*)'Tc=',dTc,'[Hour]'



dCTimeMax=dTc*3600 !48 ore dopo la fine della pioggia
a2dCTime=a2dDem !Pongo il Ctime uguale al dem perchè non lo ho
WHERE(a2dDem.gt.0)
	a2iMask=1 !Pongo la Mask uguale al dem perchè non lo ho
ENDWHERE
!-----------------------------------------------------------------------------------
!   Reading cell coordinate of the sections to use in Continuum
!-----------------------------------------------------------------------------------
CALL ReadSections(iRows,iCols)


iClass=dCTimeMax/dStep+1 !%Da mettere nella forma operativa
!	Lunghezza in dStep dell'idrogramma di portata
iLVect=iClass+dd/dStep
!	Lunghezza in dStep della simulazione	
dSimLength=dd/dStep-1	

ALLOCATE (a2dVolumIn(iNumBasins),a2dVolumInNet(iNumBasins))
ALLOCATE (a2dQT(iNumBasins,0:iLVect))
ALLOCATE (a2dQSotT(iNumBasins,0:iLVect),a2dEvTMat(10,0:iLVect))
ALLOCATE (a1dLatP(iNStazP),a1dLonP(iNStazP),a1dLatT(iNStazT),a1dLonT(iNStazT),a1dLatK(iNStazK),a1dLonK(iNStazK))
ALLOCATE (a1dLatW(iNStazW),a1dLonW(iNStazW),a1dLatUm(iNStazUm),a1dLonUm(iNStazUm))
ALLOCATE (a1dZP(iNStazP),a1dZT(iNStazT),a1dZK(iNStazK),a1dZW(iNStazW),a1dZUm(iNStazUm))

ALLOCATE(a2dAQsections(iNsec,int(dIst)))

a2dQsections => a2dAQsections
!Inizializzo il vettore delle portate
a2dQsections =0.0

ALLOCATE (dSatBac(iNumBasins),dSatFC(iNumBasins),dEvTTot(iNumBasins))

!-------------------------------------------------------------------------------------
! Se riprendo un run interrotto, leggo le variabili di stato e l'idrogramma per capire
! a che step sono. Solo con HyperCdrift.
!-------------------------------------------------------------------------------------
dTo=0

if(iFlagTypeConv.eq.10)THEN
	iLinux=10 !Linux
ELSE
	iLinux=0 !Windows
ENDIF

if(iFlagStateVar.eq.1)THEN 
	
	CALL ReadStateMatrix(iRows,iCols,sDataIni)
	!CALL ReadStateMatrixBinary(iRows,iCols,sDataIni) !Da usare nella versione operativa
	dTo=1
ENDIF


!-----------------------------------------------------------------------------------
!   Lettura dei file dei dati micrometeorologici se ho le serie temporali
!-----------------------------------------------------------------------------------
IF(iFlagInterp==1)THEN
	CALL ReadMeteoData(a1dLatP,a1dLonP,a1dZP,a1dLatT,a1dLonT,a1dZT,a1dLatK,a1dLonK,a1dZK,a1dLatW,a1dLonW,a1dZW,a1dLatUm,a1dLonUm,a1dZUm)
ENDIF
!-----------------------------------------------------------------------------------
!   Open Output Files
!-----------------------------------------------------------------------------------
!CALL open_output_files
!-----------------------------------------------------------------------------------
!	Computing total catchment area in number of cells
dBasinArea=sum(sum(a2iMask,dim=1,mask=a2iMask.gt.0.0)) !DIM=1 columns
!Scrittura dell'area del sBasin idrografico in km^2
WRITE(*,'(A11,1x,f12.1,1x,A7)')'AREA BACINO',dBasinArea*dCelLat*dCelLon/1000000,'[km^2]'

!-----------------------------------------------------------------------------------
!	Chiamata alla Subroutine Storage (dove si calcola a2dS a partire dal CN)
!-----------------------------------------------------------------------------------
CALL storage(iRows,iCols)

DEALLOCATE(a2dACurNum)

!----------------------------------------------------------------------------------- 
!	Intervallo dei parametri sui quali si fa la calibrazione (il terzo argomento 
!	indica il passo a2dCon cui vengono presi i parametri)      
!----------------------------------------------------------------------------------- 

!	Scrittura del file di controllo di inizio simulazione
OPEN(unit=500,file='INIZIO.out')
WRITE(500,*)'Simulation Start'
CLOSE(500)

nn=0


			
dVolumInPrec=0.0
dVolumInNetPrec=0.0
a2dVolumIn(1)=0.0
a2dVolumInNet(1)=0.0
	
		
		!File dove vengono scritti gli idrogrammi
OPEN(210,file=sBasin(1:iNameLenght)//'HydrographContinuum.out',status='unknown', &
		form='formatted',access='sequential')
		!File dove vengono scritti i volumi
OPEN(unit=111,file=sBasin(1:iNameLenght)//'_volumes.out')
!File dove vengono scritte le TS a scala di bacino (solo sperimentale)
OPEN(unit=311,file=sBasin(1:iNameLenght)//'_Ts.out') 
		

!-----------------------------------------------------------------------------------
!	Chiamata alla Subroutine Initialisation
!-----------------------------------------------------------------------------------
if(iFlagStateVar.ne.1)THEN !La eseguo solo se non ho già le variabili di di stato
	CALL initialisation
ENDIF
!-----------------------------------------------------------------------------------
!	Chiamata alla Subroutine wt_bedrock
!-----------------------------------------------------------------------------------
IF(iFlagDeepFlow.eq.1)THEN
	CALL wt_bedrock(iRows,iCols)
ENDIF
!-----------------------------------------------------------------------------------	  
	   
nn=nn+1

WRITE(6,'(A4,I4,1x,A4,f3.1,1x,A4,f4.2)')'run',nn,' ct=',dCt,' cf=',dCf

	dHyTot=SUM(SUM(a2dHydro,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	dRoutTot=SUM(SUM(a2dRouting,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	dVtot=SUM(SUM(a2dV,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	!dTotIa=dTotIa+SUM(SUM(a2dRetention,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	dTotIa=SUM(SUM(a2dRetention,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	dcS=SUM(SUM(a2dS,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	dTot=dTot+SUM(SUM(vol_sot,DIM=1,MASK=a2dDem.GT.0.0)) !DIM=1 columns
	dTotDeepFlow=SUM(SUM((a2dVwt-a2dVwtMax)*1000,DIM=1,MASK=a2dDem.GT.0.0)) !dTotDeepFlow+SUM(SUM(a2dVLoss,DIM=1,MASK=a2dDem.GT.0.0)) 

WRITE(111,'(50(f32.10,1x))')  &
	a2dVolumIn(1), &      
	dVtot*dCelLat*dCelLon/1000, &
	dTot*dCelLat*dCelLon/1000,	&											
	dEvTTot(1)*dCelLat*dCelLon/1000, &										
	dRoutTot*dCelLat*dCelLon/1000, &		
	dTotIa*dCelLat*dCelLon/1000, &
	dHyTot*dCelLat*dCelLon/1000, &
	dQtot, &
	dVHypod, &
	dVWtot, &
	dVtot/dBasinArea, &
	dVtot2/dBasinArea, &
	dcS/dBasinArea, &
	a2dVolumIn(1)/(dBasinArea*dCelLat*dCelLon)*1000, &
	dTotDeepFlow*dCelLat*dCelLon/1000, &
	dErr

h=0.0
iStart=0
		
!Se riprende un run interrotto salto per gli step temporali necessari i dati in lettura
!dTstart va messo da codice, valutare se levare questa parte
DO t=1,dTstart/dx
	IF(iFlagInterp==1)THEN
		CALL SaltaMeteoData(iRows,iCols,a1dLatP,a1dLonP,a1dZP,a1dLatT,a1dLonT,a1dZT,a1dLatK,a1dLonK,a1dZK,a1dLatW,a1dLonW,a1dZW,a1dLatUm,a1dLonUm,a1dZUm)
	ENDIF
	dTo=1
ENDDO

!-----------------------------------------------------------------------------------------
!   	Inizia ciclo sul tempo, dura fino alla fine della Simualzione
!-----------------------------------------------------------------------------------------  
dDintegrPrec=dDintegr !Inizializzo il dt precedente
DO t=dTo,dSimLength,1

!     Contatore dei giorni			
	mm=int(t/(dDayStep*dIst))
			
!     Contatore dello dStep del giorno partendo dalle 00.00
	h=h+1 
	m=t/dIst

	IF (m.eq.(t/dIst)) then
        IF (mm.eq.t/(dDayStep*dIst)) WRITE(*,'(A5,1x,f5.0)')'day=',mm
            WRITE(*,*)'h=',m

!-----------------------------------------------------------------------------------------
!   Metorological data - Interpolation or Reading 
!-----------------------------------------------------------------------------------------  
		IF(iFlagInterp==1)THEN
			! Meto data interpolation
			CALL MeteoDataInterp(iRows,iCols,a1dLatP,a1dLonP,a1dZP,a1dLatT,a1dLonT,a1dZT,a1dLatK,a1dLonK,a1dZK,a1dLatW,a1dLonW,a1dZW,a1dLatUm,a1dLonUm,a1dZUm)
			! Updating Date		    
			CALL DeterminaData_stepParam(xstring,xstring_up,int(dMinutes))
			xstring_do=xstring
			xstring=xstring_up
	
		ELSE
		    ! Reading Meteorological Input Maps
            IF(iFlagGrid.eq.1)THEN	
				!Meteo and Dem grids different 
				CALL ReadMeteoMapBinaryGeoloc(xstring,iRows,iCols)
			ELSE
				!Meteo and Dem grids equal
				CALL ReadMeteoMapBinary(xstring,iRows,iCols)
			ENDIF
			! Updating Date
		    CALL DeterminaData_stepParam(xstring,xstring_up,int(dMinutes))
			
			xstring_do=xstring
			xstring=xstring_up
		ENDIF

					
		WHERE(a2dRain.lt.0.0.and.a2dDem.GT.0.0)a2dRain=0.0
		a2dVolumIn(1)=a2dVolumIn(1)+SUM(SUM(a2dRain*a2dAreaCell,DIM=1,MASK=a2dDem.GT.0.0))/1000 !DIM=1 columns

!-----------------------------------------------------------------------------------------
!   SNOW  MODULE
!-----------------------------------------------------------------------------------------  
		!If "basin"InfoSnow.txt file is not found iFlagSnow is set to 0
		IF(iFlagSnow.eq.1)THEN
			CALL SRaM(t,xstring_do,d,iRows,iCols)
		ENDIF
										
!-----------------------------------------------------------------------------------
!	Energy Balance 
!-----------------------------------------------------------------------------------				    			                
		CALL ForceRestoreEq_matrix(iRows,iCols,iStart,d,t,mm,h,dDayStep,xstring_do)
		dLST=dLST/dBasinArea
		WRITE(311,'(f6.2)') dLST
		dLST=-9999
!-----------------------------------------------------------------------------------
!	Evapotraspiration computation
!-----------------------------------------------------------------------------------
				
	
		CALL evapotranspiration_FR_matrix(iRows,iCols,d,dEvTTot,xstring_do)
		dETDayMean=dETDayMean+dEvTTot(1)/dBasinArea
		
					
		IF(h.eq.dDayStep)THEN
			WRITE(103,'(F8.2,1x)')dETDayMean
			dETDayMean=0.0
			h=0
		ENDIF

		a2dAreeWT=a2dRain
!-----------------------------------------------------------------------------------
!	Retention computation
!-----------------------------------------------------------------------------------

		CALL RetentionHdrift(iRows,iCols)
!***********************************************************************************
!-----------------------------------------------------------------------------------
!	Entro nella covoluzione di Continuum (All'interno anche Horton e deflusso Ipodermico
!-----------------------------------------------------------------------------------

					
		dSt1 =  RTC()
		CALL convolutionHyperCdrift(iRows,iCols,d,dDintegr,dStep,t,a2dVolumIn, &
				a2dVolumInNet,dTot,dTotIa,dVolTot,dTotHypodFlow,dTotDeepFlow,dTotEsf)
		dSt2 = RTC( )
		dTime_spent =dSt2-dSt1 
		WRITE(*,*)'tempo di calcolo: ',dTime_spent,' secondi'
		!Scrivo le portate
		DO tt=1,dx
			WRITE(210,600) t+((tt-1)*d/(dIst*3600)),(a2dQsections(nb,tt),nb=1,iNumBasins)
			WRITE(*,*)'tempo ',t+((tt-1)*d/(dIst*3600)),'Q=',a2dQsections(1,tt)
		ENDDO

						


!-----------------------------------------------------------------------------------
!	Chiamata alla Subroutine wt_subflow
!	Ripartizione della parte di subflow in direzione della falda profonda
!----------------------------------------------------------------------------------- 
		IF(iFlagDeepFlow.eq.1)THEN
			!Subroutine che calcola il deflusso profondo semplificato
			CALL QDeepFlow_al2(iRows,iCols,d,dTotDeepFlow)
		ENDIF

		!Scrivo la mappa di Portate attuale
		sVariabile='portata'
		!CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariabile,sPathResults,a2dQmap,100)
		!Scrivo la mappa di Ts attuale
		sVariabile='Ts'
		CALL WriteMeteoMapBinary(xstring_do,iRows,iCols,sVariabile,sPathResults,a2dTs,10,iLinux)

!-----------------------------------------------------------------------------------------  
		IF(iFlagStateSave.eq.0)THEN !Each hour
				CALL WriteStateMatrixBinary(iRows,iCols,xstring_do)
				IF(iNlake+iNdam.ge.1)THEN
					CALL SaveStateDigheLaghi(iRows,iCols,xstring_do) !Stati dei laghi
				ENDIF
				
			ELSEIF(iFlagStateSave.eq.1)THEN !Once a day
				read(xstring_do(9:10),*)dHourSD
				IF (dHourSD.eq.dHourState)then
					CALL WriteStateMatrixBinary(iRows,iCols,xstring_do)
					IF(iNlake+iNdam.ge.1)THEN
						CALL SaveStateDigheLaghi(iRows,iCols,xstring_do) !Stati dei laghi
					ENDIF

				ENDIF
			ELSEIF(iFlagStateSave.eq.2)THEN !Once a month
				read(xstring_do(9:10),*)dHourSD
				read(xstring_do(7:8),*)dMonthSD
				IF (dHourSD.eq.dHourState.and.dMonthSD.eq.dMonthState)then
					CALL WriteStateMatrixBinary(iRows,iCols,xstring_do)
					IF(iNlake+iNdam.ge.1)THEN
						CALL SaveStateDigheLaghi(iRows,iCols,xstring_do) !Stati dei laghi
					ENDIF

				ENDIF
			ENDIF
	
           IF(iFlagOutSave.eq.0)THEN !Each hour
				!Scrivo la mappa di Portate attuale
				sVariabile='Q'
				CALL WriteMeteoMapBinary(xstring_do,iRows,iCols,sVariabile,sPathResults,a2dQmap,100,iLinux)
				!Scrivo la mappa di Ts attuale
				sVariabile='Ts'
				CALL WriteMeteoMapBinary(xstring_do,iRows,iCols,sVariabile,sPathResults,a2dTs,int(dRescFct),iLinux)

				
			ELSEIF(iFlagOutSave.eq.1)THEN !Once a day
				read(xstring_do(9:10),*)dHourSD
				IF (dHourSD.eq.dHourOut)then
					!Scrivo la mappa di Portate attuale
					sVariabile='Q'
					CALL WriteMeteoMapBinary(xstring_do,iRows,iCols,sVariabile,sPathResults,a2dQmap,100,iLinux)
					!Scrivo la mappa di Ts attuale
					sVariabile='Ts'
					CALL WriteMeteoMapBinary(xstring_do,iRows,iCols,sVariabile,sPathResults,a2dTs,int(dRescFct),iLinux)

					
				ENDIF
			ELSEIF(iFlagOutSave.eq.2)THEN !Once a month
				read(xstring_do(9:10),*)dHourOut
				read(xstring_do(7:8),*)dMonthOut
				IF (dHourSD.eq.dHourOut.and.dMonthSD.eq.dMonthOut)then
					!Scrivo la mappa di Portate attuale
					sVariabile='Q'
					CALL WriteMeteoMapBinary(xstring_do,iRows,iCols,sVariabile,sPathResults,a2dQmap,100,iLinux)
					!Scrivo la mappa di Ts attuale
					sVariabile='Ts'
					CALL WriteMeteoMapBinary(xstring_do,iRows,iCols,sVariabile,sPathResults,a2dTs,int(dRescFct),iLinux)

				ENDIF
			ENDIF
			
!-----------------------------------------------------------------------------------
!	SCRITTURA FILE DI BILANCIO
!
!	Scrittura dei file di bilancio
!	La variabile a2dVolumIn indica la precipitazione, mentre a2dVolumInNet indica il ruscellamento; sia la precipitazione sia il ruscellamento non sono azzerate
!	ad ogni giro e di conseguenza vengono calcolate ad un certo istante togliendo il valore della variabile all'istante precedente. La variabile dTot indica il 
!	volume del subflow, mentre dEvTTot indica l'evapotraspirazione potenziale al cui interno vengono espresse le due componenti evapotraspirazione da terreno
!	ed evaporazione dalle piante. Il volume dTotHypodFlow è il volume invasato in un certo istante, mentre la quantità dTotIa indica il volume d'acqua trattenuto dalle 
!	piante. Queste variabili vengono invece azzerate ad ogni istante temporale diversamente a quello che viene fatto per le prime due variabili espresse. 
!	Il bilancio di massa viene fatto come:
!	P - R - EVTtot - ESFiltrazione (non calcolata qui) = Vtot + Ia
!	o analogamente in termini di variabili all'interno del programma come:
!	a2dVolumIn - a2dVolumInNet - dEvTTot - dTot = dTotHypodFlow + dTotIa
!	In generale il volume compreso in Vtot e Ia è trascurabile sul lungo periodo
!	rispetto a quello legato agli altri termini.
!	Nel file 101 vengono scritte quindi le seguenti quantità: (in ordine)
!	precipitazione  [mm] 
!	ruscellamento   [mm]
!	volume istantaneamente trasferito come sub-flow
!	evpotrasp totale [mm]
!	volume medio invasato in quell'istante  [mm]
!   dTotIa
!	volume medio perso verso lo strato profondo in quell'istante  [mm]
!	evaporazione da a2dRetention vegetativa  [mm]
!-----------------------------------------------------------------------------------
			
		dQtot=dQtot+a2dQsections(1,1)*d
		dVWtot=dVWtot+a2dAlpha(a2dXYsections(1,2),a2dXYsections(1,1))*a2dCostF1(a2dXYsections(1,2),a2dXYsections(1,1))*d/(3600*1000)*dCelLat*dCelLon
		WRITE(111,'(5000(f32.10,1x))')  &
    		a2dVolumIn(1), &      
    		dVtot*dCelLat*dCelLon/1000, &
    		dTot*dCelLat*dCelLon/1000,	&											
    		dEvTTot(1)*dCelLat*dCelLon/1000, &										
    		dRoutTot*dCelLat*dCelLon/1000, &		
    		dTotIa*dCelLat*dCelLon/1000, &
			dHyTot*dCelLat*dCelLon/1000, &
			dQtot, &
			dVHypod, &
			dVWtot, &
			dVtot/dBasinArea, &
			dVtot2/dBasinArea, &
			dcS/dBasinArea, &
			a2dVolumIn(1)/(dBasinArea*dCelLat*dCelLon)*1000, &
			dTotDeepFlow*dCelLat*dCelLon/1000, &
			dErr

!-----------------------------------------------------------------------------------

			dVolumInPrec=a2dVolumIn(1)
			dVolumInNetPrec=a2dVolumInNet(1)

!-----------------------------------------------------------------------------------
!     Fine if sull passo temporale
!-----------------------------------------------------------------------------------

	ENDIF

! Inizializzo a 0 la matrice dell'esfiltrazione per l'istante successivo
	a2dEsf=0.0

!-----------------------------------------------------------------------------------
!	Fine loop sulla durata della simulazione Ciclo sul tempo t
!-----------------------------------------------------------------------------------
      
ENDDO
!-----------------------------------------------------------------------------------
!	Ciclo per calcolo in Hcdrift dell'idrogramma per tc dopo dSimLength
!-----------------------------------------------------------------------------------
DO t=dSimLength+1,iLVect
	a2dRain=0.0 !Pongo la pioggia a 0 (i dati meteo sono finiti)
	CALL convolutionHyperCdrift(iRows,iCols,d,dDintegr,dStep,t,a2dVolumIn, &
				a2dVolumInNet,dTot,dTotIa,dVolTot,dTotHypodFlow,dTotDeepFlow,dTotEsf)
	!Scrivo le portate
			
	WRITE(210,600) t+((tt-1)*d/(dIst*3600)),(a2dQsections(nb,tt),nb=1,iNumBasins)
	WRITE(*,*)'tempo ',t+((tt-1)*d/(dIst*3600)),'Q=',a2dQsections(1,tt)

	IF(iFlagDeepFlow.eq.1)THEN
	    !Subroutine che calcola il deflusso profondo semplificato
		CALL QDeepFlow_al2(iRows,iCols,d,dTotDeepFlow)
	ENDIF
	dQtot=dQtot+a2dQsections(1,1)*d
	dVWtot=dVWtot+a2dAlpha(a2dXYsections(1,2),a2dXYsections(1,1))*a2dCostF1(a2dXYsections(1,2),a2dXYsections(1,1))*d/(3600*1000)*dCelLat*dCelLon*dKsatRatio										   

	WRITE(111,'(50(f32.10,1x))')  &
  				a2dVolumIn(1), &      
   				dVtot*dCelLat*dCelLon/1000, &
   				dTot*dCelLat*dCelLon/1000,	&											
   				dEvTTot(1)*dCelLat*dCelLon/1000, &										
   				dRoutTot*dCelLat*dCelLon/1000, &		
   				dTotIa*dCelLat*dCelLon/1000, &
				dHyTot*dCelLat*dCelLon/1000, &
				dQtot, &
				dVHypod, &
				dVWtot, &
				dVtot/dBasinArea, &
				dVtot2/dBasinArea, &
				dcS/dBasinArea, &
				a2dVolumIn(1)/(dBasinArea*dCelLat*dCelLon)*1000, &
				dTotDeepFlow*dCelLat*dCelLon/1000, &
				dErr
ENDDO


 
iStart=0

!-----------------------------------------------------------------------------------
!	Apertura e Scrittura dei file per idrogramma, subflow e ruscellamento sul bacino (runoff)
!-----------------------------------------------------------------------------------
!	Definizione dei format di Scrittura
600 FORMAT (f8.1,600(f13.2))
601 FORMAT (600(f9.2,1x))


!	Scrittura del file di controllo di fine simulazione
OPEN(unit=500,file='FINE.out')
WRITE(500,*)'Fine della Simulazione'
CLOSE(500)

!-----------------------------------------------------------------------------------
!   De-allocazione delle matrici e dei vettori
!-----------------------------------------------------------------------------------
DEALLOCATE(a2iAMask,a2dACon)
DEALLOCATE(a2dAS,a2dACTime,a2dADem,a2iAPun,a2iAChoice)
DEALLOCATE(a2dACostF,a2dACostK,a2dACostF1,a2dACostChFix)
DEALLOCATE(a2dAAlpha,a2dAKAlpha,a2dAWaterTable)

DEALLOCATE(a2dAVwt,a2dADeepFlow,a2dAVwtMax,a2dABeta)

DEALLOCATE(a2dAC1,a2dAF2)

DEALLOCATE(Arain,Acumrain,Acumrain_prec)
DEALLOCATE(ATemp,AK,AW,AUm,APres)
DEALLOCATE(Asat_m,Aesf,Avol_sot,AV,AVloss)
DEALLOCATE(Aevapot,Aintensity,Aritenzione)
DEALLOCATE(AEF_mat,AW_prec,ATT_prec,AUm_prec,AK_prec)
DEALLOCATE(AT_segnato,ATa_vet24,ATs)



DEALLOCATE(a2dAHydro,a2dARouting,a2dARoutingV2C)
DEALLOCATE(a2dAXYsections,a2dAQmap,a2dALai,a2dAChn)
DEALLOCATE(a2dAAreaCell)

DEALLOCATE(dSatBac,dSatFC,dEvTTot)
DEALLOCATE(a2dAQsections)

DEALLOCATE (a1dLatP,a1dLonP,a1dLatT,a1dLonT,a1dLatK,a1dLonK,a1dLatW,a1dLonW,a1dLatUm,a1dLonUm)
DEALLOCATE (a1dZP,a1dZT,a1dZK,a1dZW,a1dZUm)
DEALLOCATE (a2dVolumIn,a2dVolumInNet)

!-----------------------------------------------------------------------------------
!   End Main
!-----------------------------------------------------------------------------------
END PROGRAM

!-----------------------------------------------------------------------------------



