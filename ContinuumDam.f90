!***********************************************************************************
!	 Continuum:
!
!	- Modified Horton method for infiltration
!	- Runoff Routing;
!	- Subsurface flow Routing;
!	- Complete Energy Balance using Force Restore equation to compute soil temperature;
!	- Tdeep filter;
!	- WaterTable  and Deep Flow routing.
!***********************************************************************************

PROGRAM MainContinuumDam

IMPLICIT NONE 
INCLUDE 'DeclarationH.f90' ! Common variables declaration

!-----------------------------------------------------------------------------------
!	Dichiarazione dei TARGET dei vettori o delle matrici dinamiche comuni
!-----------------------------------------------------------------------------------
!	Matrici Tempoinvarianti
REAL*8,  ALLOCATABLE, TARGET :: a2dACTime(:,:),a2dADem(:,:),a2dAS(:,:),a2dACurNum(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dACostF(:,:),a2dACostK(:,:),a2dACostF1(:,:),a2dACostChFix(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a2dAAlpha(:,:),a2dAKAlpha(:,:),a2dAWaterTable(:,:)	
REAL*8,  ALLOCATABLE, TARGET :: a2dAKAlphaZero(:,:)
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

! Matrici utilizzate per la gestione delle Dighe, Rilasci, ..etc
REAL*8,  ALLOCATABLE, TARGET :: a1dADamVolume(:),a1dAIdroTurbinate(:,:)
INTEGER,  ALLOCATABLE, TARGET :: a2dAXYDam(:,:),a2dAXYCen(:,:)
CHARACTER*500, ALLOCATABLE, TARGET :: a1sANameTurbinate(:),a1sANameRilasci(:),a1sANamePrese(:)
REAL*8,  ALLOCATABLE, TARGET :: a1dAQturbMax(:),a1dATcorrturb(:)
INTEGER,  ALLOCATABLE, TARGET :: a2dAXYPresa(:,:),a2dAXYRilascio(:,:)
REAL*8,  ALLOCATABLE, TARGET :: a1dATcorrprese(:),a1dAIdroPrese(:,:),a1dAIdroRila(:,:),a1dAPesoPresa(:)

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
INTEGER i,j,n,nn,nb,tt,iD
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
CHARACTER*12  xstring,xstring_up,xstring_do,sDataIni
CHARACTER*80 sBuffer

real*8 dDintegr !Dt del routing
real*8 dTstart	!Tempo di inizio della simulazione
real*8 dTo! !Vale 1 per ripresa run e 0 per per run dall'inizio periodo
real*8 dSt1,dSt2,RTC,dTime_spent !Variabili per calcolo dei tempi
real*8 dHourSD !Hour saving data
!***********************************************************************************
!	Argomenti da riga di comando
!***********************************************************************************
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
CALL ReadInfoDataDam(iRows,iCols,iFlagInterp,xstring,dDintegr,sDataIni)
!iFlagVcVar=1 !Se 1 dCappaC variable se no costante METTERE IN INFO!

!	Calcolo del passo temporale in secondi dei dati in ingresso
iStepRain=dHours/dMinutes*60
dd=3600*dHours
d=60*dMinutes
dIst=dx
dStep=d/dIst
dCTimeMax=0
dDayStep=24*60/dMinutes

!	Dimensioni della matrice che conserva le coord dei pixel contenuti nel pixel SAT
iNsec=iNumBasins !Numero sezioni lette su file
!-----------------------------------------------------------------------------------
!	Allocazione dei TARGET dei vettori o delle matrici comuni
!-----------------------------------------------------------------------------------
!	Matrici Tempoinvarianti
ALLOCATE(a2iAMask(iRows,iCols),a2dACon(iRows,iCols))
ALLOCATE(a2dACTime(iRows,iCols),a2dADem(iRows,iCols),a2iAPun(iRows,iCols),a2dAS(iRows,iCols),a2dACurNum(iRows,iCols))
ALLOCATE(a2iAChoice(iRows,iCols))
ALLOCATE(a2dACostF(iRows,iCols),a2dACostK(iRows,iCols),a2dACostF1(iRows,iCols),a2dACostChFix(iRows,iCols))
ALLOCATE(a2dAAlpha(iRows,iCols),a2dAKAlpha(iRows,iCols),a2dAWaterTable(iRows,iCols))
ALLOCATE(a2dAKAlphaZero(iRows,iCols))
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
a2dHydro=0.000001
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
	a2dAreaCell=(dCelLat*dCelLon);	
ENDIF

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
IF((dXDemLat.ne.dXCNLat).or.(dXDemLon.ne.dXCNLon))THEN 
	ALLOCATE(a2iAYLat(iRows,iCols),a2iAXLon(iRows,iCols))
	
	a2iYLat			      => a2iAYLat	         
	a2iXLon				  => a2iAXLon	
	CALL Create_Meteo_GRID(iRows,iCols,iRowsMeteo,iColsMeteo)
	iFlagGrid=1
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
!   Reading sections Rows Cols in Hypercddrift
!-----------------------------------------------------------------------------------
CALL ReadSectionsDam(iRows,iCols)
!	Massimo tempo di corrivazione 

iClass=0 !dCTimeMax/dStep+1 !%Da mettere nella forma operativa
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

ALLOCATE(a2dAQsections(iNsec,0:iLVect))

a2dQsections => a2dAQsections
!Inizializzo il vettore delle portate
a2dQsections =0.0

ALLOCATE (dSatBac(iNumBasins),dSatFC(iNumBasins),dEvTTot(iNumBasins))
!-------------------------------------------------------------------------------------
! Se riprendo un run interrotto, leggo le variabili di stato e l'idrogramma per capire
! a che step sono. Solo con HyperCdrift.
!-------------------------------------------------------------------------------------
dTo=0
if(iFlagStateVar.eq.1)THEN 
	!CALL ReadHydrograph(dTstart,d)
	!CALL ReadStateMatrix(iRows,iCols,sDataIni)
	CALL ReadStateMatrixBinary(iRows,iCols,sDataIni) !Da usare nella versione operativa
	dTo=1
ENDIF

if(iFlagTypeConv.eq.10)THEN
	iLinux=10 !Linux
ELSE
	iLinux=0 !Windows
ENDIF

!-----------------------------------------------------------------------------------
!   Lettura dei file dei dati micrometeorologici 
!-----------------------------------------------------------------------------------
IF(iFlagInterp==1)THEN
	CALL ReadMeteoData(a1dLatP,a1dLonP,a1dZP,a1dLatT,a1dLonT,a1dZT,a1dLatK,a1dLonK,a1dZK,a1dLatW,a1dLonW,a1dZW,a1dLatUm,a1dLonUm,a1dZUm)
ENDIF

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
!	leggo informazioni su numero di dighe e di centrali
!-----------------------------------------------------------------------------------


!Opening info dighe
OPEN(20,file=sFileInfoDighe,status='old')
READ(20,*)iNdam !Numero di dighe
READ(20,*)iNcentr !Numero centrali 
CLOSE(20)

ALLOCATE(a1dADamVolume(iNdam),a2dAXYDam(iNdam,2),a2dAXYCen(iNcentr,2),a1dAIdroTurbinate(iNcentr,int(dSimLength)+10))
ALLOCATE(a1sANameTurbinate(iNcentr),a1dAQturbMax(iNcentr),a1dATcorrturb(iNcentr))
!Inizializzo matrici relativi alle dighe
a1dDamVolume	=>	a1dADamVolume
a2dXYDam		=>	a2dAXYDam
a2dXYCen		=>	a2dAXYCen
a1dIdroTurbinate	=>	a1dAIdroTurbinate
a1sNameTurbinate	=>	a1sANameTurbinate
a1dQturbMax		=>	a1dAQturbMax
a1dTcorrturb	=>	a1dATcorrturb


a1dDamVolume=0.0
a2dXYDam=0.0
a2dXYCen=0.0
a1dIdroTurbinate=0.0 !*1000*d/(dCelLat*dCelLon)
a1dQturbMax=0.0
a1dTcorrturb=0.0

! Legge tutte le informazioni dall' infoDighe
CALL ReadInfoDighe(iRows,iCols)
!Leggo la serie dei turbinati
CALL ReadTurbinatiDam(iRows,iCols,dSimLength,d)


!-----------------------------------------------------------------------------------
!	leggo informazioni su numero di prese e rilasci
!-----------------------------------------------------------------------------------


OPEN(20,file=sFileInfoPrese,status='old')
READ(20,*)iNprese !Numero di prese
READ(20,*)iNril !Numero di rilasci
CLOSE(20)


ALLOCATE(a2dAXYPresa(iNprese,2),a2dAXYRilascio(iNril ,2))
ALLOCATE(a1dATcorrprese(iNprese),a1dAIdroPrese(iNprese,int(dSimLength)+10),a1dAIdroRila(iNril,int(dSimLength)+10),a1dAPesoPresa(iNprese))
ALLOCATE(a1sANameRilasci(iNril),a1sANamePrese(iNprese))


a1sNameRilasci	=>	a1sANameRilasci
a1sNamePrese	=>	a1sANamePrese
a2dXYPresa	=>	a2dAXYPresa
a2dXYRilascio	=>	a2dAXYRilascio
a1dTcorrprese	=>	a1dATcorrprese
a1dIdroPrese	=>	a1dAIdroPrese
a1dIdroRila	=>	a1dAIdroRila
a1dPesoPresa	=>	a1dAPesoPresa

a2dXYPresa=0.0
a2dXYRilascio=0.0
a1dTcorrprese=0.0
a1dIdroPrese=0.0
a1dIdroRila=0.0
a1dPesoPresa=0.0

! Legge tutte le info dell'infoPrese
CALL ReadInfoPrese(iRows,iCols)
!Legge la serie dei turbinati dai rilasci
CALL ReadPreseRilasciDam(iRows,iCols,dSimLength,d)
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
	
		h=0.0
		iStart=0
	    !Se riprende un run interrotto salto per gli step temporali necessari i dati in lettura
		!dTstart va messo da codice, valutare se levare questa parte

		dTo=1
		DO t=1,dTstart/dx
			IF(iFlagInterp==1)THEN
				CALL SaltaMeteoData(iRows,iCols,a1dLatP,a1dLonP,a1dZP,a1dLatT,a1dLonT,a1dZT,a1dLatK,a1dLonK,a1dZK,a1dLatW,a1dLonW,a1dZW,a1dLatUm,a1dLonUm,a1dZUm)
			ENDIF
			dTo=1
		ENDDO

		!Ciclo sul tempo
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
						! Meteo data interpolation
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

!-----------------------------------------------------------------------------------
!	Energy Balance 
!-----------------------------------------------------------------------------------				    			                
					CALL ForceRestoreEq_matrix(iRows,iCols,iStart,d,t,mm,h,dDayStep,xstring_do)

!-----------------------------------------------------------------------------------
!	Evapotraspiration computation
!-----------------------------------------------------------------------------------

					
					CALL evapotranspiration_FR_matrix(iRows,iCols,d,dEvTTot,xstring_do)
					dETDayMean=dETDayMean+dEvTTot(1)/dBasinArea
										
					IF(h.eq.17.0*dDayStep/24.0)THEN
						WRITE(103,'(F7.4,1x)')dETDayMean
						dETDayMean=0.0
					ENDIF
!-----------------------------------------------------------------------------------
!	Retention computation
!-----------------------------------------------------------------------------------

					CALL RetentionHdrift(iRows,iCols)
!***********************************************************************************
!-----------------------------------------------------------------------------------
!	Entro nella covoluzione di HyperCdrift (All'interno anche Horton e deflusso Ipodermico
!-----------------------------------------------------------------------------------

					dSt1 =  RTC()
					CALL convolutionHyperCdriftDam(iRows,iCols,d,dDintegr,dStep,t,a2dVolumIn, &
						a2dVolumInNet,dTot,dTotIa,dVolTot,dTotHypodFlow,dTotDeepFlow,dTotEsf)
					dSt2 = RTC( )
					dTime_spent =dSt2-dSt1 
					WRITE(*,*)'tempo di calcolo: ',dTime_spent,' secondi'
					!Scrivo le portate
					DO tt=1,dx
						WRITE(210,600) ((t+tt-1)*d/(dIst*3600)),(a2dQsections(nb,int(t)+tt-1),nb=1,iNumBasins)
						WRITE(*,*)'tempo ',t+tt-1,'Q=',a2dQsections(1,int(t)+tt-1)
					ENDDO

!-----------------------------------------------------------------------------------
!	Chiamata alla Subroutine wt_subflow
!	Ripartizione della parte di subflow in direzione della falda profonda
!----------------------------------------------------------------------------------- 
					IF(iFlagDeepFlow.eq.1)THEN
					    !Subroutine che calcola il deflusso profondo semplificato
						CALL QDeepFlow_al2(iRows,iCols,d,dTotDeepFlow)

					ENDIF
!-----------------------------------------------------------------------------------------  

					read(xstring_do(9:10),*)dHourSD
					IF (dHourSD.eq.0.or.dHourSD.eq.12)then
						CALL WriteStateMatrixBinary(iRows,iCols,xstring_do)
						!CALL WriteStateMatrix(iRows,iCols,xstring_do)
						h=0.0				
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
					dQtot=dQtot+a2dQsections(1,int(t))*d
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
					(a1dDamVolume(iD),iD=1,iNdam), & !Volume invasato nelle dighe
					dErr


!-----------------------------------------------------------------------------------

					dVolumInPrec=a2dVolumIn(1)
					dVolumInNetPrec=a2dVolumInNet(1)


!-----------------------------------------------------------------------------------
!     Fine if sull passo temporale
!-----------------------------------------------------------------------------------
				!DEALLOCATE(rain2file) !??? DA LEVARE?

	ENDIF

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
			
	WRITE(210,600) (t*d/(dIst*3600)),(a2dQsections(nb,int(t)),nb=1,iNumBasins)
	WRITE(*,*)'tempo ',t,'Q=',a2dQsections(1,int(t))
	IF(iFlagDeepFlow.eq.1)THEN
	    !Subroutine che calcola il deflusso profondo semplificato
		CALL QDeepFlow_al2(iRows,iCols,d,dTotDeepFlow)
	ELSE
		a2dDeepFlow=0.0
	ENDIF
ENDDO



!-----------------------------------------------------------------------------------
!	Apertura e Scrittura dei file per idrogramma, subflow e ruscellamento sul bacino (runoff)
!-----------------------------------------------------------------------------------
!	Definizione dei format di Scrittura
600 FORMAT (f8.1,600(f13.2))
601 FORMAT (600(f9.2,1x))

!	Scrittura del file dei parametri
!DO i=1,3
!	WRITE(109,'(600(f9.2,1x))')(par(j,i),j=1,nn)
!ENDDO


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


DEALLOCATE(a1dADamVolume,a2dAXYDam,a2dAXYCen,a1dAIdroTurbinate)
DEALLOCATE(a1sANameTurbinate,a1sANameRilasci,a1sANamePrese)
DEALLOCATE(a1dAQturbMax,a1dATcorrturb)
DEALLOCATE(a2dAXYPresa,a2dAXYRilascio)
DEALLOCATE(a1dATcorrprese,a1dAIdroPrese,a1dAIdroRila,a1dAPesoPresa)
!-----------------------------------------------------------------------------------
!   End Main
!-----------------------------------------------------------------------------------
END PROGRAM

!-----------------------------------------------------------------------------------



