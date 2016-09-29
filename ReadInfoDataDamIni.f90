subroutine ReadInfoDataDamIni(iRows,iCols,iFlagInterp,xstring,dDintegr,sDataIni)

    IMPLICIT NONE
	INCLUDE 'DeclarationH.f90' ! 

    INTEGER iRows,iCols
    CHARACTER str*500
	CHARACTER*12 xstring,sDataIni
	CHARACTER*50 info,text
    external setIniFilename, getValue
    logical error, simulate
    real pi
	REAL*8 dDintegr
	INTEGER*4 iFlagInterp

	!Opening bacino.info.txt
	info=sBasin(1:iNameLenght)//'Dam.info.txt' 
    call setIniFilename(info)

    !Operationa System 
    call getValue('', 'SO', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(i4)') iFlagTypeConv
        write(*,'(a,i4)') 'SO = ', iFlagTypeConv
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'SO'
		stop
    endif
	!Computing deep Flow (1=yes, 0=no)
    call getValue('', 'DeepFlow', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(i4)') iFlagDeepFlow
        write(*,'(a,i4)') 'DeepFlow = ', iFlagDeepFlow
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'DeepFlow'
		stop
    endif
	!Restart a run (1=yes, 0=no)
	call getValue('', 'RR', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(i4)') iFlagStateVar
        write(*,'(a,i4)') 'RR = ', iFlagStateVar
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'RR'
		stop
    endif
	!Meteo data grouping (1=one dir; 2=year,aaaa ; 3=year-month, aaaamm)
	call getValue('', 'MeteoGroup', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dMeteoGroup
        write(*,'(a,f5.2)') 'MeteoGroup = ', dMeteoGroup
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'MeteoGroup'
		stop
    endif
	!Number of sections
	call getValue('', 'iSectionNumber', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) iNumBasins
        write(*,'(a,i4)') 'iSectionNumber= ', iNumBasins
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'iSectionNumber'
		stop
    endif
	!Lower Left DEM angle coordinate	
    call getValue('', 'GeoDEM', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *)dXDemLat, dXDemLon
        write(*,'(a,f9.6,f9.6)') 'GeoDEM= ', dXDemLat, dXDemLon
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'GeoDEM'
		stop
    endif

    !Lower Left Meteo Data angle coordinate
	call getValue('', 'GeoMETEO', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *)dXCNLat, dXCNLon	
        write(*,'(a,f9.6,f9.6)') 'GeoMETEO= ', dXCNLat, dXCNLon	
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'GeoMETEO'
		stop
    endif
	if (dXCNLat.gt.dXDemLat) then
		write(*,*) 'LatLL_METEO > LatLL_DEM. Stop'
		stop
	endif
	if (dXCNLon.gt.dXDemLon) then
		write(*,*) 'LonLL_METEO > LonLL_DEM. Stop'
		stop
	endif
	!DEM Lat and Lon cellsizes
	call getValue('', 'passDem', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *)dDemPasLat, dDemPasLon 	
        write(*,'(a,f10.7,f10.7)') 'passDem= ', dDemPasLat, dDemPasLon 
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'passDem'
		stop
    endif
    !Meteo Data Lat and Lon cellsizes	
    call getValue('', 'passMeteo', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *)dCNPasLat, dCNPasLon	 	
        write(*,'(a,f10.7,f10.7)') 'passMeteo= ', dCNPasLat, dCNPasLon	 
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'passMeteo'
		stop
    endif
	

	IF((dXDemLat.ne.dXCNLat).or.(dXDemLon.ne.dXCNLon).or.(dDemPasLon.ne.dCNPasLon))THEN		!If the meteo and dem grids are different
	!Meteodata num rows and cols file (only if Dem and Meteo grids are different)
		call getValue('', 'DimMeteo', str, error)
		if (error .eqv. .false.) then
			write(*,'(a,a)') 'str = ', trim(str)
			read(str, *)iRowsMeteo,iColsMeteo	 	
			write(*,'(a,i4,i4)') 'DimMeteo= ', iRowsMeteo,iColsMeteo	 
		else
			write(*,'(a)') 'Error: section or keyword not found. Stop'
			write(*,*) 'DimMeteo'
			stop
		endif

	ENDIF
	!Lat and Lon cellsizes in meters
    call getValue('', 'DimDemMeters', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *)dCelLat, dCelLon		 	
        write(*,'(a,f10.2,f10.2)') 'DimDemMeters= ', dCelLat, dCelLon		 
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'DimDemMeters'
		stop
    endif
    !Rainfall Duration (hours)
	call getValue('', 'Duration', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dHours
        write(*,'(a,f15.1)') 'Duration= ', dHours
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'Duration'
		stop
    endif
	!Temporal dStep of Rainfall Matrix (minutes)																									
    call getValue('', 'iStepDelta', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dMinutes
        write(*,'(a,f15.1)') 'iStepDelta= ', dMinutes
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'iStepDelta'
		stop
    endif

	!Number of Integral Steps Between Two Rainfall Matrix (dx)
	call getValue('', 'intStep', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dx
        write(*,'(a,f6.1)') 'intStep= ', dx
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'intStep'
		stop
    endif

	!Temporal dStep of Continuum integration (seconds)
	call getValue('', 'ConvInt', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dDintegr
        write(*,'(a,f6.1)') 'ConvInt= ', dDintegr
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'ConvInt'
		stop
    endif																								
	!Meteo Data Rescaling factor (permitted: 10 or 100 or 1000)
	call getValue('', 'iScaleFactor', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dRescFct	
        write(*,'(a,f6.1)') 'iScaleFactor= ', dRescFct	
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'iScaleFactor'
		stop
    endif	
	!Interpolation flag - if Flag is equal 1  reads Meto Data Series and performs IDW interpolation, else put it equal 0
	
	call getValue('', 'IntFlag', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) iFlagInterp	
        write(*,'(a,i4)') 'IntFlag= ', iFlagInterp	
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'IntFlag'
		stop
    endif	
	!N° of Rain gauges
	call getValue('', 'Nraing', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) iNStazP	
        write(*,'(a,i4)') 'Nraing= ', iNStazP	
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'Nraing'
		stop
    endif	
	!N° of Thermometers
	call getValue('', 'Ntmpng', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) iNStazT	
        write(*,'(a,i4)') 'Ntmpng= ', iNStazT	
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'Ntmpng'
		stop
    endif																	
    !N° of Radiometers
	call getValue('', 'Nradmg', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) iNStazK	
        write(*,'(a,i4)') 'Nradmg= ', iNStazK	
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'Nradmg'
		stop
    endif																	

    !N° of Anemometers	
	call getValue('', 'Nwindg', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) iNStazW	
        write(*,'(a,i4)') 'Nwindg= ', iNStazW	
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'Nwindg'
		stop
    endif		
	!N° of Hydrometers	
	call getValue('', 'Nhygrm', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) iNStazUm	
        write(*,'(a,i4)') 'Nhygrm= ', iNStazUm	
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'Nhygrm'
		stop
    endif		


	!Run starting data	
	call getValue('', 'sDateStart', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a12)') xstring	
        write(*,'(a,a)') 'sDateStart= ', xstring		
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sDateStart'
		stop
    endif		

	!Run initial conditions	
	call getValue('', 'sDateStatus', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a12)') sDataIni	
        write(*,'(a,a)') 'sDateStatus= ', sDataIni		
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sDateStatus'
		stop
    endif	
    
	!Rainfall Data/Map Path			
	call getValue('', 'sPathDataTimeRain', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)') sPathRain	
        write(*,'(a,a)') 'sPathDataTimeRain= ', sPathRain		
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sPathDataTimeRain'
		stop
    endif

	!Temperature Data/Map Path	
	call getValue('', 'sPathDataTimeTemp', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)') sPathTemp	
        write(*,'(a,a)') 'sPathDataTimeTemp= ', sPathTemp		
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sPathDataTimeTemp'
		stop
    endif
	
	!Radiation Data/Map Path	
	call getValue('', 'sPathDataTimeRadi', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)') sPathRad	
        write(*,'(a,a)') 'sPathDataTimeRadi= ', sPathRad		
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sPathDataTimeRadi'
		stop
    endif	

	!Wind Data/Map Path	
	call getValue('', 'sPathDataTimeWind', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)') sPathWind	
        write(*,'(a,a)') 'sPathDataTimeWind= ', sPathWind		
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sPathDataTimeWind'
		stop
    endif	
	
	!Rel Humidity Data/Map Path	
	call getValue('', 'sPathDataTimeRelu', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)') sPathUmid	
        write(*,'(a,a)') 'sPathDataTimeRelu= ', sPathUmid		
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sPathDataTimeRelu'
		stop
    endif	
	
		

    !Land Data Path
    call getValue('', 'sPathDataStaticLandData', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sPathLandData
        write(*,'(a,a)') 'sPathLandData = ', trim(sPathLandData)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sPathLandData'
		stop
    endif

	!Path 2D output
    call getValue('', 'sOut2D', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sPathResults
        write(*,'(a,a)') 'sOut2D = ', trim(sPathResults)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sOut2D'
		stop
    endif

	!Quasi static Data/Map Path (LAI)
    call getValue('', 'sQuasiStatic', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sPathLai
        write(*,'(a,a)') 'sQuasiStatic = ', trim(sPathLai)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sQuasiStatic'
		stop
    endif

    !File info dams
    call getValue('', 'sInfoDams', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sFileInfoDighe
        write(*,'(a,a)') 'sInfoDams = ', trim(sFileInfoDighe)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sInfoDams'
		stop
    endif

	!File info outlet-release
    call getValue('', 'sInfoOutlet', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sFileInfoPrese
        write(*,'(a,a)') 'sInfoOutlet = ', trim(sFileInfoPrese)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sInfoOutlet'
		stop
    endif

	!Path hydraulic plants series 
    call getValue('', 'sPlantsSeries', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sPathTurb	
        write(*,'(a,a)') 'sPlantsSeries = ', trim(sPathTurb	)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sPlantsSeries'
		stop
    endif

	!File info Lakes
    call getValue('', 'sInfoLakes', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sFileInfoLaghi	
        write(*,'(a,a)') 'sInfoLakes= ', trim(sFileInfoLaghi	)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sInfoLakes'
		stop
    endif

	!File info join tributaries
    call getValue('', 'sInfoJoin', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sFileInfoJoin		
        write(*,'(a,a)') 'sInfoJoin= ', trim(sFileInfoJoin)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sInfoJoin'
		stop
    endif

    !Snow info file
    call getValue('', 'sInfoSnow', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, '(a)')sFileInfoSnow		
        write(*,'(a,a)') 'sInfoSnow= ', trim(sFileInfoSnow)
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'sInfoSnow'
		stop
    endif

	!State save modality
    call getValue('', 'dFlagS', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *)iFlagStateSave,dMonthState,dHourState 	 	
        write(*,'(a,i4,x,f10.7,x,f10.7)') 'dFlagS= ', iFlagStateSave,dMonthState,dHourState 	 
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'dFlagS'
		stop
    endif

	!Output save modality
    call getValue('', 'dFlagO', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *)iFlagOutSave,dMonthOut,dHourOut 	 	
        write(*,'(a,i4,x,f10.7,x,f10.7)') 'dFlagO= ', iFlagOutSave,dMonthOut,dHourOut 	 
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'dFlagO'
		stop
    endif

	!Frequency of MeteoData(hours)
	call getValue('', 'iStepFrequency', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dFreqMeteo
        write(*,'(a,f15.1)') 'iStepFrequency= ', dFreqMeteo
    else
		write(*,'(a)') 'Error: section or keyword not found. Stop'
		write(*,*) 'iStepFrequency'
		stop
    endif
	dFreqMeteo=dFreqMeteo/60
	!Number of minutes waiting for Meteo data
	call getValue('', 'dMinWait', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dMinWait
        write(*,'(a,f15.1)') 'dMinWait= ', dMinWait
    else
		write(*,'(a)') 'Error: section or keyword not found. Continue With'
		write(*,*) 'dMinWait=0'
		dMinWait=0
    endif
	
    !Number of hours of routing after last observation
	call getValue('', 'dTmaxFor', str, error)
    if (error .eqv. .false.) then
        write(*,'(a,a)') 'str = ', trim(str)
        read(str, *) dTmaxFor
        write(*,'(a,f15.1)') 'dTmaxFor= ', dTmaxFor
    else
		write(*,'(a)') 'Error: section or keyword not found. Continue With'
		write(*,*) 'dTmaxFor=-9999; estimated corrivation time used'
		dTmaxFor=-9999
    endif
return
end subroutine