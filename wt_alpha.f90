!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************


subroutine wt_alpha(iRows,iCols,DD)
implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,ii,iii,j,jj,jjj,kkk,ttt
integer a,b,perc_tot

real*8 pi,DD,diff,fn,fMean,fNumPen
real*8 pdistance(iRows,iCols),LDD(iRows,iCols),mask_perc_tot(iRows,iCols)
real*8 diff_DD(iRows,iCols),pend(iRows,iCols),pend2(iRows,iCols),pend3(iRows,iCols)
real*8 dBmin,dBmax,dem_max
real*8 dDistanceT

!------------------------------------------------------------------

LDD=0.0
pdistance=0.0
diff_DD=-9999.0
a2dAlpha=-9999.0
mask_perc_tot=-9999.0
pend=0

dDistanceT=500
IF(dCelLat.ge.100.and.dCelLat.lt.1000)THEN
     dDistanceT=2000
ENDIF
IF(dCelLat.ge.5000.and.dCelLat.lt.20000)THEN
     dDistanceT=30000
ENDIF

WHERE(a2dDem.le.0.and.a2dDem.gt.-1000)
	a2dDem=0.2
ENDWHERE
!Effettuo la correzione del DEM per non avere Laghi
CALL Laghi(iRows,iCols)

perc_tot=0
!Il ciclo lascia righe e colonne di bordo perchè se il DEM
!arriva sul bordo ho problemi numerici
DO i=1,iRows
	DO j=1,iCols 
				
		a=i
		b=j


		IF(a2dDem(i,j).gt.0.0)THEN
			perc_tot=perc_tot+1
			fNumPen=0
			DO WHILE((a2dDem(a,b).gt.0.0).and.diff_DD(a,b).eq.-9999)

				IF((a.gt.0.and.a.le.iRows).and.(b.gt.0.and.b.le.iCols))THEN
				!IF((a.gt.1.and.a.le.iRows-1).and.(b.gt.1.and.b.le.iCols-1))THEN	
					iii=a+(INT((a2iPun(a,b)-1)/3)-1)
					jjj=b+a2iPun(a,b)-5-3*(INT((a2iPun(a,b)-1)/3)-1)
					LDD(a,b)=SQRT(((a-iii)*dCelLat)**2+((b-jjj)*dCelLon)**2)
					IF(iii.lt.1.or.jjj.lt.1)THEN
						EXIT
					ENDIF
					!write(*,*)jjj,iii,' ',a,b
					pdistance(a,b)=LDD(a,b)					
					diff_DD(a,b)=a2dDem(a,b)-a2dDem(iii,jjj)
					mask_perc_tot(a,b)=perc_tot

					!Pendenza media sui canali
					if(datan2(diff_DD(a,b),LDD(a,b)).gt.0.0)then
						fNumPen=fNumPen+1
						pend(a,b)=pend(a,b)+datan2(diff_DD(a,b),LDD(a,b))
					endif
					DO WHILE(a2dDem(a,b)-a2dDem(iii,jjj).le.DD.AND.(iii.gt.0.and.iii.le.iRows).and.(jjj.gt.0.and.jjj.le.iCols) &
								.and.a2dDem(iii,jjj).gt.0.0.and.LDD(a,b).lt.dDistanceT)	
					
						mask_perc_tot(a,b)=perc_tot
						diff_DD(a,b)=a2dDem(a,b)-a2dDem(iii,jjj)
						ii=iii+(INT((a2iPun(iii,jjj)-1)/3)-1)
						jj=jjj+a2iPun(iii,jjj)-5-3*(INT((a2iPun(iii,jjj)-1)/3)-1)	
										
						IF(a2dDem(a,b)-a2dDem(ii,jj).le.DD.and.(ii.gt.0.and.ii.le.iRows).and.(jj.gt.0.and.jj.le.iCols))THEN	
							LDD(a,b)=LDD(a,b)+SQRT(((ii-iii)*dCelLat)**2+((jj-jjj)*dCelLon)**2)
							!Pendenza media sui canali
							if(datan2(diff_DD(a,b),LDD(a,b)).gt.0.0)then
								if(a2iChoice(a,b).eq.1)then
									fNumPen=fNumPen+1
									pend(a,b)=pend(a,b)+datan2(diff_DD(a,b),LDD(a,b))
								endif
								if(a2iChoice(a,b).eq.0.and.LDD(a,b).lt.500)then
									fNumPen=fNumPen+1
									pend(a,b)=pend(a,b)+datan2(diff_DD(a,b),LDD(a,b))
								endif
							endif
						ENDIF
									
						iii=ii
						jjj=jj	
						
						!write(*,*)jjj,iii

						IF(diff_DD(iii,jjj).ne.-9999)THEN
							DO WHILE(a2dDem(a,b)-a2dDem(iii,jjj).le.DD.AND.(iii.gt.0.and.iii.le.iRows).and. &
										(jjj.gt.0.and.jjj.le.iCols).and.a2dDem(iii,jjj).gt.0.0.and.LDD(a,b).lt.dDistanceT)	
					    
								mask_perc_tot(a,b)=perc_tot					
								diff_DD(a,b)=a2dDem(a,b)-a2dDem(iii,jjj)
								ii=iii+(INT((a2iPun(iii,jjj)-1)/3)-1)
								jj=jjj+a2iPun(iii,jjj)-5-3*(INT((a2iPun(iii,jjj)-1)/3)-1)	
									
								IF(a2dDem(a,b)-a2dDem(ii,jj).le.DD.and.(ii.gt.0.and.ii.le.iRows).and.(jj.gt.0.and.jj.le.iCols))THEN	
									LDD(a,b)=LDD(a,b)+SQRT(((ii-iii)*dCelLat)**2+((jj-jjj)*dCelLon)**2)
									!Pendenza media sui canali
									if(datan2(diff_DD(a,b),LDD(a,b)).gt.0.0)then
										if(a2iChoice(a,b).eq.1)then
											fNumPen=fNumPen+1
											pend(a,b)=pend(a,b)+datan2(diff_DD(a,b),LDD(a,b))
										endif
										if(a2iChoice(a,b).eq.0.and.LDD(a,b).lt.500)then
											fNumPen=fNumPen+1
											pend(a,b)=pend(a,b)+datan2(diff_DD(a,b),LDD(a,b))
										endif
									endif											
								ENDIF		
								iii=ii
								jjj=jj

							ENDDO									
						ENDIF

					ENDDO					
					
					if(fNumPen.gt.0.0)then	
						pend(a,b)=pend(a,b)/fNumPen
					endif
										
					a2dAlpha(a,b)=datan2(DD,LDD(a,b))  !Angolo in radianti
					if(diff_DD(a,b).lt.0.9.or.diff_DD(a,b).gt.500)then
						diff_DD(a,b)=0.9
					endif
					if(diff_DD(a,b).lt.1.and.LDD(a,b).lt.4*dCelLat)then
						LDD(a,b)=4*dCelLat
					endif

					a2dAlpha(a,b)=datan2(diff_DD(a,b),LDD(a,b))
					    	
					ii=a+(INT((a2iPun(a,b)-1)/3)-1)
					jj=b+a2iPun(a,b)-5-3*(INT((a2iPun(a,b)-1)/3)-1)
					IF(a2dDem(ii,jj).gt.0.0)THEN
						a=ii
						b=jj
						fNumPen=0
	
					ELSE
						EXIT !esce ma conserva gli indici della fine percorso svolto
					ENDIF

				ENDIF
		     
			ENDDO !FINE DI UN PERCORSO COMPLETO SEGUENDO I PUNTATORI		
		
			ii=a+(INT((a2iPun(a,b)-1)/3)-1)
			jj=b+a2iPun(a,b)-5-3*(INT((a2iPun(a,b)-1)/3)-1)
		       
		ENDIF
        
	ENDDO
ENDDO

a2dBeta=pend


WHERE(a2iChoice.lt.1)
	pend=0
ENDWHERE
pend2=pend
pend=0

!Smoothing della pendenza sui canali
DO i=1,iRows 
	DO j=1,iCols
		if(a2iChoice(i,j).eq.1)then
			fn=0
			DO ii=i-1,i+1
				DO jj=j-1,j+1
					if(pend2(ii,jj).gt.0.0)then
						fn=fn+1
						pend(i,j)=pend(i,j)+pend2(ii,jj)
					endif
				ENDDO
			ENDDO
			if(fn.lt.1)fn=1
			pend(i,j)=pend(i,j)/fn
			if(LDD(i,j).le.4*dCelLat.and.diff_DD(i,j).lt.2)then
				pend(i,j)=a2dAlpha(i,j)
			endif
			if(pend(i,j).gt.0.0)then
				a2dAlpha(i,j)=pend(i,j)
				a2dBeta(i,j)=pend(i,j)
			endif
		endif
		
	ENDDO
ENDDO

!Controllo la pendenza minima che non sia sotto una certa soglia
dBmin=minval(minval(pend,DIM = 1,MASK=pend.gt.0),DIM=1)
dBmin=max(dBmin,0.00001)

WHERE(a2dDem.gt.0.and.a2dBeta.eq.0)
	a2dBeta=a2dAlpha
ENDWHERE

WHERE(a2dDem.gt.0.and.a2dBeta.lt.dBmin)
	a2dBeta=dBmin
ENDWHERE

return
end subroutine