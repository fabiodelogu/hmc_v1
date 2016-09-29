!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine ReadInfoDataDam
!	Legge il file con le informazioni delle sezioni se ho le dighe: basinDem.info.txt
!***********************************************************************************  
 
subroutine ReadInfoDataDam(iRows,iCols,iFlagInterp,xstring,dDintegr,sDataIni)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i,iFlagInterp
CHARACTER*50 info,text
CHARACTER*12 xstring,sDataIni,sTest
INTEGER*4 iLStr,iPathLenght
REAL*8 dDintegr !Dt integrazione Hcdrift
!Opening bacino.info.txt
info=sBasin(1:iNameLenght)//'Dam.info.txt' 
OPEN(20,file=info,status='old',err=543)
READ(20,*)sTest 
!Test if info file is in Ini version. If yes use parse file reader
IF(sTest=='[]')THEN
	CLOSE(20)
	write(*,*)'Info File of Ini Type'
	CALL ReadInfoDataDamIni(iRows,iCols,iFlagInterp,xstring,dDintegr,sDataIni)
	GO TO 544

ENDIF

!Reading bacino.info.txt
READ(20,*) iFlagTypeConv			!Type of Operative System (1=Windows, 10=Linux)
READ(20,*)
READ(20,*) iFlagDeepFlow			!Computing DeepFlow (1=yes, 0=no)
READ(20,*)
READ(20,*) iFlagStateVar			!Restart a run(1=yes, 0=no)
READ(20,*)
READ(20,*) dMeteoGroup				!meteo data grouping (1, 2, 3)
READ(20,*)
READ(20,*) iNumBasins				!Read Sections (basin) number
READ(20,*)
READ(20,*) dXDemLat, dXDemLon		!Lower Left DEM angle coordinate
READ(20,*)
READ(20,*) dXCNLat, dXCNLon			!Lower Left Meteo angle coordinate
READ(20,*)
READ(20,*) dDemPasLat, dDemPasLon   !DEM Lat and Lon cellsizes
READ(20,*)
READ(20,*) dCNPasLat, dCNPasLon		!Meteo Lat and Lon cellsizes
IF((dXDemLat.ne.dXCNLat).or.(dXDemLon.ne.dXCNLon).or.(dDemPasLon.ne.dCNPasLon))THEN		!If the meteo and dem grids are different
	READ(20,*)
	READ(20,*)iRowsMeteo,iColsMeteo   !Meteo dimensions Rows, Cols
ENDIF
READ(20,*)
READ(20,*) dCelLat, dCelLon			!Lat and Lon cellsizes in meters
READ(20,*)
READ(20,*)dHours					!Rainfall Duration (hours)
READ(20,*)
READ(20,*)dMinutes					!Temporal dStep of Rainfall Matrix (minutes)
READ(20,*)
READ(20,*)dx						!Number of Integral Steps Between Two Rainfall Matrix (dx)
READ(20,*)
READ(20,*)dDintegr					!Integration Step for Continuum (seconds)
READ(20,*)
READ(20,*)dRescFct					! Meteo data rescaling factor 
									!	sono contenuti nel pixel del satellite
READ(20,*)
READ(20,*)iFlagInterp               !Interpolation flag
									!	if Flag is equal 1 reads Meto Data Series 
									!	and performs IDW interpolation
!IF(iFlagInterp==1)THEN
IF(1==1)THEN
	READ(20,*)
	READ(20,*)iNStazP				!N° of Rain gauges
	READ(20,*)
	READ(20,*)iNStazT				!N° of Thermometers
	READ(20,*)
	READ(20,*)iNStazK				!N° of Radiometers
	READ(20,*)
	READ(20,*)iNStazW				!N° of Anemometers
	READ(20,*)
	READ(20,*)iNStazUm				!N° of Hygrometers
ENDIF

READ(20,*)
READ(20,'(A12)')xstring				! Run starting data 
READ(20,*)
READ(20,'(A12)')sDataIni			! Run Initial Condition data 
READ(20,*)
READ(20, '(A500)')sPathRain			!Rainfall Data/Map Path  
READ(20,*)
READ(20, '(A500)')sPathTemp			!Temperature Data/Map Path
READ(20,*)
READ(20, '(A500)')sPathRad			!Radiation Data/Map Path
READ(20,*)
READ(20, '(A500)')sPathWind			!Wind Data/Map Path
READ(20,*)
READ(20, '(A500)')sPathUmid			!Relatvie Umidity Data/Map Path
READ(20,*)
READ(20, '(A500)')sPathLandData		!Land Data Path	

READ(20,*)
READ(20, '(A500)')sPathResults		!Risultati
READ(20,*)
READ(20, '(A500)')sPathLai			!LAI Data/Map Path
READ(20,*)
READ(20, '(A500)')sFileInfoDighe	!File Info dighe
READ(20,*)
READ(20, '(A500)')sFileInfoPrese	!File Info Prese e Rilasci
READ(20,*)
READ(20, '(A500)')sPathTurb			!Path dei turbinati
READ(20,*)
READ(20, '(A500)')sFileInfoLaghi	!File info Lake
READ(20,*)
READ(20, '(A500)')sFileInfoJoin		!File info confluenze
READ(20,*)
READ(20, '(A500)')sFileInfoSnow  !File info Snow Characteristics
READ(20,*)
READ(20,*)iFlagStateSave,dMonthState,dHourState !Code for state variables saving
READ(20,*)
READ(20,*)iFlagOutSave,dMonthOut,dHourOut !Code for state variables saving
READ(20,*)
READ(20,*)dFreqMeteo !Freqeuncy of meteo data maps
dFreqMeteo=dFreqMeteo/60
CLOSE(20)

!Minuti di attesa
dMinWait=120
!Check on meteo data rescaling factor
544 if(dRescFct.eq.10.or.dRescFct.eq.100.or.dRescFct.eq.1000)then
else
	dRescFct=10 !default
endif
dMinWait=dMinWait*60
dWin=dMinWait/10

iFlagVcVar=1 !dCappaC variabile (1=yes, 0=no)
!  CTime header for extracting matrix dimension iRows & iCols 
iPathLenght = iLStr(sPathLandData)
OPEN(unit=1,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.dem.txt',status='old')
READ(1,*)text,iCols
READ(1,*)text,iRows
CLOSE(1)

IF(1.eq.0)THEN
543 write(*,*)'Info file not found, the program will be terminated'
	stop
ENDIF



return
end subroutine