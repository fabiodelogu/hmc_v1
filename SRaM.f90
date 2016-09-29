!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!-----------------------------------------------------------------------------------
! UTILIZZA OSSERVAZIONI DI TEMPERATURA, RADIAZIONE AD ONDE CORTE INCIDENTE,
! PRECIPITAZIONE
!-----------------------------------------------------------------------------------


subroutine SRaM(itime,sData,dtSec,iRows,iCols)

implicit none
include 'DeclarationH.f90'
integer*4 iSwitchNetRad,iSwitchLW,iHours
integer iRows,iCols
CHARACTER*12 sData
real*8  dDt,dtSec,dtMin,dtHour,itime
real*8  a2dLWRad(iRows,iCols)
real*8 Mm,Cm
iSwitchNetRad=0 !Calcolo SW rad con Albedo
iSwitchLW=0 !Calcolo LW rad con T
iHours=0

dtHour=3600.0/dtSec
dtMin=dtSec/60.0
write(*,*) dtHour
dDt=24            !se lo esprimo in secondo dDt/Tau Tau=86400 s
dDt=dDt*3600/86400

!Da levare!!!!
if(itime.lt.0)then
	WHERE (a2dDem.GE.500.0)
		a2dTemp=-2.0
		a2dRain=5.0
	ENDWHERE
endif

     
WHERE(a2dDem.GE.0.0)
	a2dMeanDayTemp=a2dMeanDayTemp+a2dTemp/(dtHour*24)
ENDWHERE
!---Cumulata su dt di Precipitazione solida
WHERE(a2dDem.GE.0.0.AND.a2dTemp.LT.dTrif.AND.a2dRain.GT.0.0)
	a2dSnowFall=a2dRain !a2dSnowFall+a2dRain
ENDWHERE

    
WHERE(a2dDem.GE.0.0.AND.a2dMeltingDayCum.LT.3.0)
    a2dExp=dExpRoLow !No Melting 0.033=1/30 -> raggiunge dRoMax in 30 giorni
ELSEWHERE
    a2dExp=dExpRoHigh !Si Melting   0.2=1/5 -> raggiunge dRoMax in 5 giorni
ENDWHERE

IF(itime.GT.1.AND.sData(9:10) == '00')THEN
	WHERE(a2dDem.GE.0.0.AND.a2dSnowFall.GT.3.0)
        a2dAge=0
    ELSEWHERE(a2dDem.GE.0.0.AND.a2dSWE.GT.0.0)
        a2dAge=a2dAge+1
    ENDWHERE

    WHERE(a2dDem.GE.0.0.AND.a2dSWE.EQ.0.0)
        a2dAge=0
    ENDWHERE

        
    WHERE(a2dDem.GE.0.0.AND.a2dMeanDayTemp.GT.0.0)
        a2dAlbedo=0.4+0.44*dexp(-a2dAge*0.12)
    ENDWHERE
    WHERE(a2dDem.GE.0.0.AND.a2dMeanDayTemp.LE.0.0)
        a2dAlbedo=0.4+0.44*dexp(-a2dAge*0.05)
    ENDWHERE

    
ENDIF


     !WHERE(a2dDem.GE.0.0.AND.a2dSWE.GT.0.0.AND.a2dSnowFall.GT.dThres)
WHERE(a2dDem.GE.0.0.AND.a2dSnowFall.GT.dThres)      
    a2dRoS0=67.9+51.3*dexp(a2dTemp/2.6)
ELSEWHERE 
    a2dRoS0=0
ENDWHERE

WHERE(a2dDem.GE.0.0.AND.a2dRoS0.GT.200)
    a2dRoS0=200
ENDWHERE

WHERE(a2dDem.GE.0.0.AND.a2dSWE.GT.0.0.AND.a2dSnowFall.GT.dThres)
    !a1dRoS(t)=(a1dSWE(t-1)+a1dP(t))/(a1dP(t)/dRo0+(a1dSWE(t-1)/a1dRoS(t-1))); 
     a2dRoS=a2dSWE/(a2dSnowFall/a2dRoS0+(a2dSWE-a2dSnowFall)/a2dRoS)		
ELSEWHERE(a2dDem.GE.0.0.AND.a2dSWE.GT.0.0)				
     a2dRoS=a2dRoS 												
ENDWHERE

WHERE(a2dDem.GE.0.0.AND.a2dSWE.EQ.0.0.AND.a2dSnowFall.GT.dThres)
     a2dRoS=a2dRoS0
										
ENDWHERE

a2dSnowFall=0.0

!---A new snow parameterization for the Meteo-France climate model
!Part I: validation in stand-alone experiments
!H. Douville, J.-F. Royer, J.-F. Mahfouf Climate Dynamics (1995) 12:21-35
!a1dRoS(t)=(a1dRoS(t-1)-dRomax)*exp(-.033*Dt)+dRomax

IF(itime.GT.1.AND.sData(9:10) == '00')THEN  
	WHERE(a2dDem.GE.0.0.AND.a2dSWE.GT.0.0)			
		a2dRoS=(a2dRoS-dRoMax)*dexp(-a2dExp*dDt)+dRoMax 		
    ENDWHERE
ENDIF


IF(itime.GT.1.AND.sData(9:10) == '00')THEN 
	a2dMeanDayTemp=0.0
ENDIF
IF(itime.GT.1.AND.sData(9:10) == '00')THEN  		
	a2dMeltingDayCum=0.0		
ENDIF

!---Coefficiente di melting 
IF(dXCNLat.GT.-10)THEN !Controllo in che emisfero sono
	IF(sData(5:6) == '12'.OR.sData(5:6) == '01'.OR.sData(5:6) == '02')THEN
		a2dsMc=sMcStag(1)
	ENDIF
	IF(sData(5:6) == '03'.OR.sData(5:6) == '04'.OR.sData(5:6) == '05')THEN
		a2dsMc=sMcStag(2)
	ENDIF
	IF(sData(5:6) == '06'.OR.sData(5:6) == '07'.OR.sData(5:6) == '08'.OR.sData(5:6) == '09')THEN
		a2dsMc=sMcStag(3)
	ENDIF
	IF(sData(5:6) == '10'.OR.sData(5:6) == '11')THEN
		a2dsMc=sMcStag(4)
	ENDIF
ELSE
	IF(sData(5:6) == '12'.OR.sData(5:6) == '01'.OR.sData(5:6) == '02')THEN
		a2dsMc=sMcStag(3)
	ENDIF
	IF(sData(5:6) == '03'.OR.sData(5:6) == '04'.OR.sData(5:6) == '05')THEN
		a2dsMc=sMcStag(4)
	ENDIF
	IF(sData(5:6) == '06'.OR.sData(5:6) == '07'.OR.sData(5:6) == '08'.OR.sData(5:6) == '09')THEN
		a2dsMc=sMcStag(1)
	ENDIF
	IF(sData(5:6) == '10'.OR.sData(5:6) == '11')THEN
		a2dsMc=sMcStag(2)
	ENDIF
ENDIF

     
WHERE(a2dsMc.LT.0.0)
    a2dsMc=0
ENDWHERE

WHERE(a2dDem.GT.0.0.AND.a2dSWE.GT.0.0.AND.a2dNature.EQ.51)
	a2dsMcGlac=a2dsMc
ELSEWHERE(a2dDem.GT.0.0.AND.a2dSWE.EQ.0.0.AND.a2dNature.EQ.51)
    a2dsMcGlac=3.0*a2dsMc
ENDWHERE



     !---AGGIORNA LA MAPPA PIOGGIA/NEVE
 a2dLWRad=0.0
CALL Snow(int(dtMin),iRows,iCols,iHours,dTrif,a2dsMc,a2dsMcGlac,a2dDem,a2dTemp,a2dRain,a2dK,a2dLWRad,a2dAlbedo,a2dNature, &
			iSwitchNetRad,iSwitchLW,a2dMelting,a2dSWE)
!---Cumulata giornaliera di Melting !DA VERIFICARE Land sulla T!!!
WHERE(a2dDem.GE.0.0.AND.a2dMelting.GT.0.0)
    a2dMeltingDayCum=a2dMeltingDayCum+a2dMelting
ENDWHERE


!---GHIACCIAI sui ghiacciai, individuati da Carta della Natura,se SWE è 
!inferiore a 100 mm (che a seconda della densità della neve sono 
! 10-25 cm) impongo uno spessore fisso
!WHERE(a2dCdN.EQ.51.AND.a2dSWE.LT.100)
!	a2dSWE=100
!ENDWHERE

!Mm=SUM(SUM(a2dMelting,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!Cm=SUM(SUM(a2dSWE,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea
!899 FORMAT (1000(A6,1x,f9.2))
!write(*,*)'**** Neve*****'
!write(*,899)'Mm:',Mm,' Cm:',Cm 
!write(*,*)'**************'


WHERE(a2dDem.LT.0.0)
   a2dRain=0.0
   a2dSWE=-9999.0
   a2dSnowFall=-9999.0
ENDWHERE
     
!Vento, Radiazione e Temperatura pari a 0 sotto la neve
!Il bilancio di energia fa tendere la LST a 0
WHERE(a2dDem.GE.0.0.AND.a2dSWE.GT.2.0)
	a2dK=0.0
    a2dTemp=0.0
	a2dW=0.1
ENDWHERE
     
IF(sData(9:10) == '00')THEN  
    iHours=0
ENDIF

!Check on rainfall
WHERE(a2dRain.LT.0.0)
   a2dRain=0.0
ENDWHERE
     


  !----------------------------------------------------------------------------------

END subroutine










