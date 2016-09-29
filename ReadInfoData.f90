!***********************************************************************************
!	Subroutine read_info_data
!	Legge il file sBasin.info.txt
!***********************************************************************************  
 
subroutine ReadInfoData(iRows,iCols,iFlagInterp,xstring,dDintegr,sDataIni)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i,iFlagInterp
CHARACTER*50 info,text
CHARACTER*12 xstring,sDataIni
INTEGER*4 iLStr,iPathLenght
REAL*8 dDintegr !Dt integrazione Hcdrift
!Opening bacino.info.txt
info=sBasin(1:iNameLenght)//'.info.txt' 
OPEN(20,file=info,status='old')
  
!Reading bacino.info.txt
READ(20,*)
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
READ(20,*) dXCNLat, dXCNLon			!Lower Left CN angle coordinate
READ(20,*)
READ(20,*) dDemPasLat, dDemPasLon   !DEM Lat and Lon cellsizes
READ(20,*)
READ(20,*) dCNPasLat, dCNPasLon		!CN Lat and Lon cellsizes
IF((dXDemLat.ne.dXCNLat).or.(dXDemLon.ne.dXCNLon))THEN		!If the meteo and dem grids are different
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
READ(20,*)dRescFct					!Passo delle immagini satellitari per individare i pixel del sBasin che 
									!	sono contenuti nel pixel del satellite
READ(20,*)
READ(20,*)iFlagInterp               !Interpolation flag
									!	if Flag is equal 1 reads Meto Data Series 
									!	and performs IDW interpolation
!IF(iFlagInterp==1)THEN
IF(1==1)THEN
	READ(20,*)
	READ(20,*)iNStazP		!N° of Rain gauges
	READ(20,*)
	READ(20,*)iNStazT		!N° of Thermometers
	READ(20,*)
	READ(20,*)iNStazK		!N° of Radiometers
	READ(20,*)
	READ(20,*)iNStazW		!N° of Anemometers
	READ(20,*)
	READ(20,*)iNStazUm		!N° of Hygrometers
ENDIF

READ(20,*)
READ(20,'(A12)')xstring		! Run starting data 
READ(20,*)
READ(20,'(A12)')sDataIni		! Run Initial Condition data 
READ(20,*)
READ(20, '(A500)')sPathRain		!Rainfall Data/Map Path  
READ(20,*)
READ(20, '(A500)')sPathTemp		!Temperature Data/Map Path
READ(20,*)
READ(20, '(A500)')sPathRad		!Radiation Data/Map Path
READ(20,*)
READ(20, '(A500)')sPathWind		!Wind Data/Map Path
READ(20,*)
READ(20, '(A500)')sPathUmid		!Relatvie Umidity Data/Map Path

READ(20,*)
READ(20, '(A500)')sPathLandData !Land Data Path

READ(20,*)
READ(20, '(A500)')sPathResults  !Risultati
READ(20,*)
READ(20, '(A500)')sPathLai  !LAI Data/Map Path
READ(20,*)
READ(20, '(A500)')sFileInfoSnow  !File info Snow Characteristics
READ(20,*)
READ(20,*)iFlagStateSave,dMonthState,dHourState !Code for state variables saving
READ(20,*)
READ(20,*)iFlagOutSave,dMonthOut,dHourOut !Code for state variables saving

CLOSE(20)

!Check on meteo data rescaling factor
if(dRescFct.eq.10.or.dRescFct.eq.100.or.dRescFct.eq.1000)then
else
	dRescFct=10 !default
endif
iFlagVcVar=1 !dCappaC variabile (1=yes, 0=no)

!  CTime header for extracting matrix dimension iRows & iCols 
iPathLenght = iLStr(sPathLandData)
OPEN(unit=1,file=sPathLandData(1:iPathLenght)//sBasin(1:iNameLenght)//'.dem.txt',status='old')
READ(1,*)text,iCols
READ(1,*)text,iRows
CLOSE(1)

return
end subroutine