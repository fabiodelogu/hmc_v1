!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!-----------------------------------------------------------------------------------
!	Dichiarazione dei parameter comuni al programma
!-----------------------------------------------------------------------------------
	integer iNumBasins,dNSimul,dShift
	integer iNCoord
	real*8 dNull,dCPI,dPorosity,dBc,fFlagRouting,dRateMin

	parameter (dShift=2)			! dShift in ore per la determinazione di Tdeep
	parameter (dNull=-9999.0)	! Dichiarazione del dNull 
	parameter (dNSimul=700)		! Numero delle simulazioni in funzione delle coppie di parametri
	parameter (fFlagRouting=1)	!Flag routing 1:uso pendenza, 0 no
    parameter (iNCoord=1)       !Orba
	parameter (dPorosity=0.4)	! Porosità media del terreno [%]
	parameter (dBc=0.5)	! Esponente formula dCappaCt=dCappaC*h^dBc
	parameter (dRateMin=0.1) !rapporto minimo del HypodermicFlow
!-----------------------------------------------------------------------------------
!	Dichiarazione delle variabili utilizzate in comune al programma
!-----------------------------------------------------------------------------------

	integer iNStazP,iNStazT,iNStazW,iNStazUm,iNStazK,iLinux
	integer*4 iNameLenght,iNameLenghtRun
	real*8 dHours,dMinutes,dx
	real*8 dCelLat,dCelLon,dDemPasLat,dDemPasLon,dXDemLat,dXDemLon, & 
				dCNPasLat,dCNPasLon,dXCNLat,dXCNLon
	real*8 dSatPas
	real*8 dBasinArea,dFlagLai
	character*250 sBasin,sRunName
	character*500 sPathRain,sPathTemp,sPathRad,sPathUmid,sPathWind,sPathLandData,sPathResults,sPathLai,sPathCh
	real*8 dDintegrPrec,dLST
	integer iRowsMeteo,iColsMeteo
	character*1 sBar
	real*8 dRescFct,dMeteoGroup,dMonthState,dHourState,dMonthOut,dHourOut,dFreqMeteo,dHourM
	integer iFlagStateSave,iFlagOutSave
	real*8 dMinWait,dWin,dTmaxFor

    common /com/ iNStazP,iNStazT,iNStazW,iNStazUm,iNStazK,iNumBasins,dCPI,iLinux
	common /com/ iNameLenght,iNameLenghtRun
	common /com/ dHours,dMinutes,dx
	common /com/ dCelLat,dCelLon,dDemPasLat,dDemPasLon,dXDemLat,dXDemLon, &
				       dCNPasLat,dCNPasLon,dXCNLat,dXCNLon
	common /com/ dSatPas
    common /com/ dBasinArea
	common /com/ sBasin,sRunName
	common /com/ sPathRain,sPathTemp,sPathRad,sPathUmid,sPathWind,sPathLandData,sPathResults,sPathLai,sPathCh
	common /com/ dFlagLai
	common /com/ dDintegrPrec,dLST
	common /com/ iRowsMeteo,iColsMeteo
	common /com/ sBar
	common /com/ dRescFct,dMeteoGroup,dMonthState,dHourState,dMonthOut,dHourOut,dFreqMeteo,dHourM
	common /com/ iFlagStateSave,iFlagOutSave
	common /com/ dMinWait,dWin,dTmaxFor
!-----------------------------------------------------------------------------------
!	Dichiarazione dei Target dei vettori o delle matrici comuni al programma
!-----------------------------------------------------------------------------------
!	Matrici Tempoinvarianti
	integer, pointer ::  a2iMask(:,:),a2dCon(:,:),a2iPun(:,:),a2iChoice(:,:)
	real*8, pointer  ::  a2dCTime(:,:),a2dDem(:,:),a2dS(:,:),a2dCurNum(:,:)
	real*8, pointer  ::  a2dCostF(:,:),a2dCostK(:,:),a2dCostF1(:,:),a2dCostChFix(:,:)
	real*8, pointer  ::  a2dAlpha(:,:),a2dKAlpha(:,:),a2dWaterTable(:,:)
	real*8, pointer  ::  a2dKAlphaZero(:,:)
	real*8, pointer  ::  a2dAreeWT(:,:),a2dDiffMaxMatrix(:,:)
	real*8, pointer  ::  a2dCt(:,:),a2dCf(:,:),a2dCappaC(:,:),a2dCappaV(:,:)
	integer, pointer  ::  a2iYLat(:,:), a2iXLon(:,:)

	common /com/  a2iMask,a2dCon,a2iPun,a2iChoice
	common /com/  a2dCTime,a2dDem,a2dS,a2dCurNum
	common /com/  a2dCostF,a2dCostK,a2dCostF1,a2dCostChFix
	common /com/  a2dAlpha,a2dKAlpha,a2dWaterTable
	common /com/  a2dKAlphaZero
	common /com/  a2dAreeWT,a2dDiffMaxMatrix
	common /com/  a2dCt,a2dCf,a2dCappaC,a2dCappaV
	common /com/  a2iYLat,a2iXLon

!-----------------------------------------------------------------------------------
!	Matrici Tempovarianti
	real*8, pointer  ::  a2dRain(:,:),a2dCumRain(:,:),a2dCumRainPrec(:,:)
	real*8, pointer  ::  a2dTemp(:,:),a2dK(:,:),a2dW(:,:),a2dUm(:,:),a2dPres(:,:)
	real*8, pointer  ::  a2dSatM(:,:),a2dEsf(:,:),vol_sot(:,:),a2dV(:,:),a2dVLoss(:,:)
	real*8, pointer  ::  a2dEvapot(:,:),a2dIntensity(:,:),a2dRetention(:,:)
	real*8, pointer  ::  a2dEF(:,:),a2dWPrec(:,:),a2dTempPrec(:,:),a2dUmPrec(:,:),a2dKPrec(:,:)
	real*8, pointer  ::  a3dTMarked(:,:,:),a3dTemp24(:,:,:),a2dTS(:,:)
	real*8, pointer  ::  a2dVLossAreeWTStep(:,:),a2dWaterTableUpdateMatrix(:,:)
	real*8, pointer  ::  a2dVolumeZeroMatrix(:,:),a2dVolumeUpDateMatrix(:,:)
	real*8, pointer  ::  a2dVLossAreeWTStepCum(:,:),a2dVLossCum(:,:)
	real*8, pointer	 ::  a2dC1(:,:),a2dF2(:,:)
	real*8, pointer	 ::  a2dEvtCum(:,:)

	integer, pointer ::  a2iNPixel(:)

	common /com/  a2dRain,a2dCumRain,a2dCumRainPrec
	common /com/  a2dTemp,a2dK,a2dW,a2dUm,a2dPres
	common /com/  a2dSatM,a2dEsf,vol_sot,a2dV,a2dVLoss
	common /com/  a2dEvapot,a2dIntensity,a2dRetention
	common /com/  a2dEF,a2dWPrec,a2dTempPrec,a2dUmPrec,a2dKPrec
	common /com/  a3dTMarked,a3dTemp24,a2dTS
	common /com/  a2iNPixel
	common /com/  a2dVLossAreeWTStep,a2dWaterTableUpdateMatrix
	common /com/  a2dVolumeZeroMatrix,a2dVolumeUpDateMatrix
	common /com/  a2dVLossAreeWTStepCum,a2dVLossCum
	common /com/  a2dC1,a2dF2
	common /com/  a2dEvtCum
!-----------------------------------------------------------------------------------
!	Matrici Continuum
	real*8, pointer :: a2dHydro(:,:),a2dRouting(:,:),a2dQsections(:,:),a2dRoutingV2C(:,:)
	real*8, pointer :: a2dVwt(:,:),a2dDeepFlow(:,:),a2dQmap(:,:),a2dLai(:,:),a2dChn(:,:)
	real*8, pointer :: a2dVwtMax(:,:),a2dBeta(:,:)
	real*8, pointer :: a2dAreaCell(:,:)
	integer, pointer :: a2dXYsections(:,:)
	real*8, pointer :: a2dQsubssections(:,:)

	real*8 dCappaC,dCappaV,dCt,dCf,dHbr,dKsatRatio,dSlopeMax 

	integer*4 iNsec,iFlagTypeConv,iFlagDeepFlow,iFlagStateVar,iFlagVcVar,iFlagCh

	real*8 dHyTot,dRoutTot,dVtot,dQtot,dcS,dcV,dErr,dVtot2,dCappaMax

	common /com/  a2dHydro,a2dRouting,a2dXYsections,a2dQsections,a2dRoutingV2C,dKsatRatio
	common /com/  a2dVwt,a2dDeepFlow,a2dQmap,a2dLai,a2dChn,a2dVwtMax,dHbr,a2dBeta
	common /com/  a2dAreaCell
	common /com/  a2dQsubssections

	common /com/  dCappaC,dCappaV,dCt,dCf,iNsec,iFlagTypeConv,iFlagDeepFlow,iFlagStateVar,iFlagVcVar,iFlagCh

	common /com/ dHyTot,dRoutTot,dVtot,dQtot,dcS,dcV,dErr,dVtot2,dCappaMax,dSlopeMax 
!-----------------------------------------------------------------------------------	
!	Variabili per le dighe

	integer, pointer :: a2dXYDam(:,:),a2dXYCen(:,:),a2dXYPresa(:,:),a2dXYRilascio(:,:),a1iDigaCentrale(:)
	real*8, pointer :: a1dDamVolume(:)
	real*8, pointer :: a1dIdroTurbinate(:,:),a1dQturbMax(:),a1dTcorrturb(:),a1dIdroPrese(:,:),a1dIdroRila(:,:)
	real*8, pointer :: a1dQmaxRil(:),a1dQminEco(:)
	real*8, pointer :: a1dTcorrprese(:),a1dPesoPresa(:)
	integer*4 iNdam,iNcentr,iNprese,iNril,iFlagPrese
	real*8 dVHypod,dVWtot

	character*500,pointer ::  a1sNameTurbinate(:),a1sNameRilasci(:),a1sNamePrese(:)
	character*500 sFileInfoDighe,sPathTurb,sFileInfoPrese

	common /com/  a2dXYDam,a1dDamVolume,a1dIdroTurbinate,iNdam,sFileInfoDighe,a1sNameTurbinate
	common /com/  a2dXYPresa,a2dXYRilascio,a1dIdroPrese,a1dIdroRila,sFileInfoPrese,a1sNameRilasci,a1sNamePrese
	common /com/  a1dQmaxRil,a1dQminEco
	common /com/  a2dXYCen,iNcentr,sPathTurb,a1dQturbMax,a1dTcorrturb,iNprese,iNril,a1dTcorrprese,a1dPesoPresa,a1iDigaCentrale
	common /com/  dVHypod,dVWtot,iFlagPrese

!-----------------------------------------------------------------------------------	
!	Variabili per le dighe con lago a monte

	real*8, pointer :: a1dVdamMax(:),a1dCodeDam(:),a1dNumCodeDam(:),a1dQ_sLC(:),a1dCostDighe(:)
	real*8, pointer :: a1dVlake(:),a1dVlakeMin(:),a1dQ_laghi(:),a1dCostLaghi(:),a1dCodeLake(:),a1dNumCodeLake(:)
	integer*4 iNlake 
	integer, pointer :: a2dXYLake(:,:)
	character*500 sFileInfoLaghi
	real*8, pointer :: a1d_Level(:,:),a1d_Volume(:,:),a1dHmax(:),adL(:),a1dCoefDighe(:)

	common /com/  a1dVdamMax,a1dCodeDam,a1dNumCodeDam,a1dQ_sLC,a1dCostDighe
	common /com/  a1dVlake,a1dVlakeMin,a1dQ_laghi,a1dCostLaghi,iNlake,a2dXYLake,sFileInfoLaghi,a1dCodeLake,a1dNumCodeLake

	common /com/  a1d_Level,a1d_Volume,a1dHmax,adL,a1dCoefDighe

!-----------------------------------------------------------------------------------	
!	Variabili per unione fiumi e risalita corrente
	integer, pointer :: a2dXYMain(:,:),a2dXYImm(:,:),a2dXYOut(:,:)
	real*8, pointer :: a1dThreshLiv(:)
	integer*4 iNjoin
	character*500 sFileInfoJoin
	common /com/  a2dXYMain,a2dXYImm,a2dXYOut,iNjoin,a1dThreshLiv,sFileInfoJoin

!-----------------------------------------------------------------------------------	
!	Variabili per neve
	real*8 dTrif,dRoS0,dRoMax,dThres,dExpRoLow,dExpRoHigh,snodata,sMcStag(4),iFlagSnow
	real*8, pointer :: a2dNature(:,:),a2dSnowFall(:,:),a2dExp(:,:),a2dAge(:,:)
	real*8, pointer :: a2dSWE(:,:),a2dAlbedo(:,:),a2dRoS0(:,:),a2dRoS(:,:),a2dMeanDayTemp(:,:)
	real*8, pointer :: a2dMeltingDayCum(:,:),a2dMelting(:,:),a2dsMc(:,:),a2dsMcGlac(:,:)
	character*500 sFileInfoSnow
	common /com/ dTrif,dRoS0,dRoMax,dThres,dExpRoLow,dExpRoHigh,snodata,sMcStag,iFlagSnow
	common /com/ sFileInfoSnow
	common /com/ a2dNature,a2dSnowFall,a2dExp,a2dAge
	common /com/ a2dSWE,a2dAlbedo,a2dRoS0,a2dRoS,a2dMeanDayTemp
    common /com/ a2dMeltingDayCum,a2dMelting,a2dsMc,a2dsMcGlac

!
	real*8, pointer :: a2dRmul(:,:)
    common /com/ a2dRmul

!-----------------------------------------------------------------------------------	
!	Variables to account od spatial scale problems
	real*8, pointer ::  a2dCoeffResol(:,:)
	real*8 dRateResol
	common /com/ a2dCoeffResol,dRateResol