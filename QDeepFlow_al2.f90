!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

	SUBROUTINE QDeepFlow_al2(iRows,iCols,d,dTotDeepFlow)
	implicit none

	include 'DeclarationH.f90'
		
	integer iRows,iCols,i,j,ii,jj,iii,jjj
	real*8 dTotHypodFlow,dTotDeepFlow,Kalpha_zero_mean,iNgr,dHt,dHm
	real*8 d,dVwt,gamma,a2dVwtemp(iRows,iCols),err(iRows,iCols),a2dDarcy(iRows,iCols)
	real*8 dDintegr !Dt del routing
	real*8 dB,dA
	gamma=0.0013
	a2dDeepFlow=0.0
	a2dDarcy=0.0
	a2dVwtemp=0.0
	!dB=SUM(SUM(a2dVwt,DIM=1,MASK=a2dDem.GT.0.0))+SUM(SUM(a2dVLoss,DIM=1,MASK=a2dDem.GT.0.0))

	WHERE(a2dDem.GT.0.0)
			a2dVwtemp=a2dVwt+a2dVLoss/1000
	ENDWHERE
 

	!dA=dVwt+SUM(SUM(a2dVwtemp,DIM=1,MASK=a2dDem.GT.0.0))
	DO j=3,iCols-1
		DO i=3,iRows-1
		iNgr=0
		dHt=0
		dHm=0.0
		ii=int((a2iPun(i,j)-1)/3)-1
		jj=a2iPun(i,j)-5-3*ii
		iii=i+ii
		jjj=j+jj
		IF (a2dDem(i,j).gt.0.0.and.a2dDem(iii,jjj).gt.0.0) then
			DO ii=i-1,i+1
				DO jj=j-1,j+1
					if(a2dDem(ii,jj).gt.0.0.and.(ii.ne.i.and.jj.ne.j))then
						if((a2dVwt(i,j)-a2dVwt(ii,jj)).gt.0)then
							
							dHt=dHt+(a2dVwt(i,j)-a2dVwt(ii,jj))
							iNgr=iNgr+1
						endif
					endif

				ENDDO
			ENDDO
			IF(iNgr.gt.0)THEN
				dHm=dHt/iNgr
			
				a2dDarcy(i,j)=dHm/sqrt(a2dAreaCell(i,j))*a2dCostF1(i,j)*d/3600*dKsatRatio
				if(a2dDarcy(i,j).gt.(a2dVwt(i,j)-a2dVwtMax(i,j))*1000)a2dDarcy(i,j)=(a2dVwt(i,j)-a2dVwtMax(i,j))*1000

				DO ii=i-1,i+1
					DO jj=j-1,j+1
						if(a2dDem(ii,jj).gt.0.0.and.(ii.ne.i.and.jj.ne.j))then
							if((a2dVwt(i,j)-a2dVwt(ii,jj)).gt.0)then
								a2dVwtemp(ii,jj)=a2dVwtemp(ii,jj)+a2dDarcy(i,j)*(a2dVwt(i,j)-a2dVwt(ii,jj))/(dHt*1000)
								dHm=dHm
							endif
						endif

					ENDDO
				ENDDO
			ENDIF

		ENDIF
		!Cella di chiusura
		IF (a2dDem(i,j).gt.0.0.and.a2dDem(iii,jjj).lt.0.0) then
			a2dDarcy(i,j)=a2dAlpha(i,j)*a2dCostF1(i,j)*d/(3600*1000)*dKsatRatio
			if(a2dDarcy(i,j).gt.(a2dVwt(i,j)-a2dVwtMax(i,j)))a2dDarcy(i,j)=(a2dVwt(i,j)-a2dVwtMax(i,j))

		ENDIF
	
		  
		ENDDO
	ENDDO

	WHERE(a2dDem.GT.0.0)
			a2dVwtemp=a2dVwtemp-a2dDarcy/1000
	ENDWHERE
	dVwt=a2dDarcy(a2dXYsections(1,2),a2dXYsections(1,1))
	
	
	!Se si fa interagire la WT con la superficie
	WHERE(a2dDem.GT.0.0.and.a2dVwtemp.gt.a2dDem)
		a2dDeepFlow=(a2dVwtemp-a2dDem)*d/3600*1000
		a2dVwtemp=a2dDem
	ENDWHERE
	!Aggiorni i Vwt
	WHERE(a2dDem.GT.0.0)
		a2dVwt=a2dVwtemp
	ENDWHERE

	!----
	WHERE(a2dDem.GT.0.0)
		a2dV=a2dV+a2dDeepFlow
	ENDWHERE
	a2dDeepFlow=0
	WHERE(a2dDem.GT.0.0.and.a2dV.gt.a2dS)
		a2dDeepFlow=a2dV-a2dS
		a2dV=a2dS
	ENDWHERE
	!----
	!dA=dVwt+SUM(SUM(a2dVwt,DIM=1,MASK=a2dDem.GT.0.0))+SUM(SUM(a2dDeepFlow,DIM=1,MASK=a2dDem.GT.0.0))
	
	
	dTotDeepFlow=SUM(SUM((a2dVwt-a2dVwtMax)*1000,DIM=1,MASK=a2dDem.GT.0.0)) !dTotDeepFlow+SUM(SUM(a2dVLoss,DIM=1,MASK=a2dDem.GT.0.0)) 

	return
	END
