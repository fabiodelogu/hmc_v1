!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!     Subroutine Horton
!***********************************************************************************

subroutine HortonMatrixHcdrift(iRows,iCols,d,dDintegr,a2dVolumIn,a2dVolumInNet,dTot,dTotIa,dVolTot)
    
implicit none	
include 'DeclarationH.f90'

integer iRows,iCols
integer i,j
real*8 d,dTot,b
real*8 dDintegr !Dt del routing
real*8 dDth !Dt su cui devo integrare Horton
real*8 a2dG(iRows,iCols)
real*8 a2dCh(iRows,iCols)
real*8 a2dErr(iRows,iCols),a2dVtemp(iRows,iCols)


!GUARDARE SE NECESSARIO TENERE a2dVolumIn MULTIPLI E CONTROLLARE PASSAGGI SUBS

real*8 a2dVolumIn,a2dVolumInNet,dVolTot,dTotIa

!Pongo dDth=dDintegr, ma potrei differenziare canale e versante	
dDth=dDintegr

dTot=0.0
dTotIa=0.0
dVolTot=0.0

!a2dS in mm and a2dRain a2dIntensity in mm/hr    


a2dErr=0.0 !Errore bilancio finale
a2dVtemp=a2dV !V passo precedente
   
WHERE (a2dRain.lt.0.0.or.a2dRain.gt.845.0) a2dRain=0.0
		 
!a2dVolumIn=a2dVolumIn+SUM(SUM(a2dRain,DIM=1,MASK=a2dDem.GT.0.0))*dCelLat*dCelLon/1000 !DIM=1 columns

!Equazione del filtro di Horton		       
a2dG=a2dCostF-(a2dCostF-a2dCostF1)/a2dS*a2dV
b=dCf

!a2dCh=((1-dCt)*a2dCostF+dCt*a2dCostF1)/((1-dCt)*a2dS) !Da levare
a2dCh=a2dCostChFix

WHERE (a2dV.lt.a2dCt*a2dS.AND.a2dDem.gt.0.0)
			vol_sot=0.0
ELSEWHERE(a2dDem.gt.0.0)
			vol_sot=a2dF2*(a2dV-a2dCt*a2dS)/a2dS*dDth/3600  ! *dDth/3600 perchè a2dCostF1 è in mm/h ma lavoro in mm/dDth
ENDWHERE
!-----------------------------------------------------------------------------------------------------------------------------------
!CAMBIARE A2DRain se no viene annullata, serve una variabile temporanea
!Oppure porta fuori ritenzione


!L'ingresso ad Horton è dato dalla pioggia più il routing delle celle a monte		
a2dIntensity=a2dRain*3600.0/d+a2dCoeffResol*a2dRouting/dDth*3600.0+a2dDeepFlow*3600.0/d
!a2dIntensity=a2dRain*3600.0/d+1*a2dRouting/dDth*3600.0+a2dDeepFlow*3600.0/d

! a2dIntensity=0
WHERE (a2dIntensity.EQ.0.0.AND.a2dDem.GT.0.0.AND.a2dV.ge.a2dCt*a2dS)    
	!a2dV=a2dC1/a2dF2*a2dS+(a2dV-a2dC1/a2dF2*a2dS)*exp(-a2dF2/a2dS*dDth/3600)
	a2dV=a2dV-a2dF2*(a2dV-a2dCt*a2dS)/a2dS*dDth/3600	
ENDWHERE
				
! 0<a2dIntensity<=g					
WHERE (a2dIntensity.GT.0.0.and.a2dIntensity.LE.a2dG.AND.a2dDem.GT.0.0)
						
	WHERE (a2dV.lt.a2dCt*a2dS.AND.a2dDem.GT.0.0)
		a2dV=a2dV+a2dRain/d*dDth+a2dCoeffResol*a2dRouting+a2dDeepFlow/d*dDth
		a2dIntensity=0							
	ELSEWHERE(a2dDem.GT.0.0)
		a2dV=a2dV+a2dIntensity*dDth/3600-vol_sot

		!a2dV=(a2dIntensity+a2dC1)/a2dF2*a2dS+(a2dV-(a2dIntensity+a2dC1)/a2dF2*a2dS)*exp(-a2dF2/a2dS*dDth/3600)
		a2dIntensity=0
	ENDWHERE
	

		
ENDWHERE

!a2dIntensity>g
WHERE (a2dIntensity.GT.a2dG.AND.a2dDem.GT.0.0)						
            
	!WHERE (a2dV.lt.dCt*a2dS.AND.a2dDem.GT.0.0) a2dCh=a2dCostF/a2dS

	!a2dV=a2dS*(1-exp(-a2dCh*dDth/3600))+a2dV*exp(-a2dCh*dDth/3600)
	WHERE (a2dV.lt.a2dCt*a2dS.AND.a2dDem.GT.0.0)
		a2dV=a2dV+a2dG*dDth/3600
	ELSEWHERE (a2dDem.GT.0.0)
		a2dV=a2dV+a2dG*dDth/3600-vol_sot
	ENDWHERE
	a2dIntensity=a2dIntensity-a2dG

ENDWHERE			

WHERE (a2dV.GT.a2dS.AND.a2dDem.GT.0.0)
	a2dIntensity=a2dIntensity+(a2dV-a2dS)/dDth*3600
	a2dV=a2dS			
ENDWHERE	

!a2dRain=a2dIntensity/(1000*3600.)

!-----------------------------------------------------------------------------------------------------------------------------------

WHERE (a2dRain.lt.0.0.AND.a2dDem.GT.0.0) a2dRain=0.0 !controllo 


!a2dAreeWT=a2dV+vol_sot+a2dIntensity*dDth/3600-(a2dDeepFlow/d*dDth+a2dRain/d*dDth+a2dRouting)			


!Calcolo l'errore di bilancio e lo suddivido tra V e Ruscellamento
WHERE (a2dDem.GT.0.0)	
	a2dErr=a2dDeepFlow/d*dDth+a2dRain/d*dDth+a2dCoeffResol*a2dRouting-a2dIntensity*dDth/3600-(a2dV-a2dVtemp+vol_sot)
	!a2dIntensity=a2dIntensity+a2dErr*3600/(dDth*2)
	!a2dV=a2dV+a2dErr
	!WHERE (a2dIntensity.lt.0.0)
	!	a2dIntensity=0.0
	!ENDWHERE
	WHERE (a2dV.lt.0.0)
		a2dV=0.0
	ENDWHERE
	!a2dErr=a2dRt-a2dIntensity*dDth/3600-(a2dV-a2dVtemp+vol_sot)
	!a2dErr=a2dDeepFlow/d*dDth+a2dRain/d*dDth+a2dRouting-a2dIntensity*dDth/3600-(a2dV-a2dVtemp+vol_sot)

ENDWHERE

dErr=dErr+SUM(SUM(a2dErr,DIM=1,MASK=a2dDem.GT.0.0))/dBasinArea !DIM=1 columns


RETURN
END

