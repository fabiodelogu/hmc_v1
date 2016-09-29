!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!     Subroutine Evapotraspiration
!***********************************************************************************

	SUBROUTINE evapotranspiration_FR_matrix(iRows,iCols,d,dEvTTot,xstring)
	implicit none
	
	include 'DeclarationH.f90'

	CHARACTER*12  xstring
	CHARACTER*500  sVariable
	integer i,j
	integer iRows,iCols
	real*8 d     
	real*8 dEvTTot,evt_att(iRows,iCols),Vprec(iRows,iCols),Rprec(iRows,iCols),err(iRows,iCols)
    real*8 dHourSD,dMonthSD !Hour and Month saving data
	 
	!dEvTTot=0.0
	evt_att=0.0

	!Vprec=a2dV
	!Rprec=a2dRetention

	

	WHERE (a2dCTime.gt.0.0)	 
		evt_att=a2dEvapot
	ENDWHERE

	WHERE (a2dRetention.GT.0.0.AND.a2dRetention.GT.evt_att.AND.a2dCTime.gt.0.0.AND.a2iChoice.le.1)
  		   
			a2dRetention=a2dRetention-a2dEvapot

	ELSEWHERE(a2dRetention.GT.0.0.AND.a2dRetention.lt.evt_att.AND.a2iChoice.le.1)
				
 			WHERE (a2dV.GE.(evt_att-a2dRetention))
				a2dV=a2dV-(evt_att-a2dRetention)
			ELSEWHERE 
				evt_att=a2dRetention+a2dV
				a2dV=0.0	
			ENDWHERE
			
			a2dRetention=0.0

	ELSEWHERE(a2iChoice.le.1)  !cioè a2dRetention.eq.0.0 ma non sui laghi

			WHERE(a2dV.GE.evt_att)
				a2dV=a2dV-evt_att !tolgo evt da a2dV solo quando "non piove" cioè a2dRetention=0 quando piove l'evaporazione dal suolo è trascurabile                    
			ELSEWHERE
				evt_att=a2dV
				a2dV=0.0					       
			ENDWHERE                    
				
	ENDWHERE

	a2dEvtCum=a2dEvtCum+evt_att
									                    
	!err	= Vprec+Rprec-evt_att-a2dV-a2dRetention
			
	dEvTTot=dEvTTot+SUM(SUM(evt_att,DIM=1,MASK=a2dCTime.GT.0.0)) !DIM=1 columns

	!Scrivo la mappa di EVT attuale
    IF(iFlagOutSave.eq.0)THEN !Each hour
		!Scrivo la mappa di Evt
		sVariable='evt'
		CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathResults,evt_att,int(dRescFct),iLinux)				
	ELSEIF(iFlagOutSave.eq.1)THEN !Once a day
		read(xstring(9:10),*)dHourSD
		
		IF (dHourSD.eq.dHourOut.and.dMonthOut.ne.1)then		
			sVariable='evt'
			CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathResults,evt_att,int(dRescFct),iLinux)				
		ENDIF
		IF (dHourSD.eq.dHourOut.and.dMonthOut.eq.1)then		
			sVariable='evtDaily'
			CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathResults,a2dEvtCum,int(dRescFct),iLinux)				
			a2dEvtCum=0.0
		ENDIF

	ELSEIF(iFlagOutSave.eq.2)THEN !Once a month
		read(xstring(9:10),*)dHourSD
		read(xstring(7:8),*)dMonthSD
		IF (dHourSD.eq.dHourOut.and.dMonthSD.eq.dMonthOut)then
			sVariable='evt'
			CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathResults,evt_att,int(dRescFct),iLinux)
		ENDIF
	ENDIF
			

	a2dEvapot=0.0 !Inizializzo l'evapotraspirazione
	RETURN

	END
      
!***********************************************************************************

