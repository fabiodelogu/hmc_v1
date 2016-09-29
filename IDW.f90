!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

    subroutine IDW(iRows,iCols,nstaz,latSt,longSt,zSt,val,Vmin,Vmax,lapse,matrice_var)

	implicit none
	include 'DeclarationH.f90'

	integer iRows,iCols
	integer i,j,ii,nstaz,label,conta,indici(nstaz),nstazmax
	real*8 lapse,nx,ny
	
	real*8 Vmin,Vmax
	real*8 dist(nstaz)
    real*8 Sd,val(nstaz)
    real*8 latSt(nstaz),longSt(nstaz),zSt(nstaz),zStazioni(nstaz)
	real*8 matrice_var(iRows,iCols),dDmax,dEsp

	Sd=0.0
	
	matrice_var=0.0
	
	nstazmax=4 !Numero massimo di stazioni per interpolazione in un punto
	if(nstazmax.gt.nstaz)then
		nstazmax=nstaz
	endif

	if(Vmax.eq.15)then !Se è pioggiaa. Occhio se variano i range!!!
		dDmax=1000/dCelLat*25
	else
		dDmax=1000/dCelLat*25
	endif
	do i=1,nstaz
		if(val(i).ge.Vmin.and.val(i).le.Vmax)then
			val(i)=val(i)-lapse*zSt(i)
		else
			val(i)= val(i)
		endif
	enddo


	do j=1,iCols
		do i=1,iRows
			Sd=0.0
			dist=99999.0
			if(a2dDem(i,j).ge.0.0.AND.a2dCTime(i,j).gt.0.0)then	 
	     
				do ii=1,nstaz
					nx=(longSt(ii)-dXDemLon)/(dDemPasLon)+1
					ny=(latSt(ii)-dXDemLat)/(dDemPasLat)+1
					if(val(ii).ge.Vmin.and.val(ii).le.Vmax) then
						dist(ii)=sqrt(dble(j-nx)**2+dble(i-ny)**2)
						if(dist(ii).le.0.0001) dist(ii)=1. 
						!if(dist(ii).le.dDmax )then
						!	!dEsp=1.5+1*val(ii)/100
						!	matrice_var(i,j)=matrice_var(i,j)+val(ii)/(dist(ii)**2)
						!	Sd=Sd+1/(dist(ii)**2)
						!endif
					endif
				enddo
				!Ordino le stazione
				CALL indexx(nstaz,dist,indici)
				do ii=1,nstazmax
					if(val(indici(ii)).ge.Vmin.and.val(indici(ii)).le.Vmax) then
						if(dist(indici(ii)).le.dDmax )then
						dEsp=2 !3-exp(-val(indici(ii))/100)
						matrice_var(i,j)=matrice_var(i,j)+val(indici(ii))/(dist(indici(ii))**dEsp)
						Sd=Sd+1/(dist(indici(ii))**dEsp)
						endif
					endif
				enddo
				if(Sd.eq.0.0)then
					Sd=1.0
					matrice_var(i,j)=0.0
				endif
				matrice_var(i,j)=matrice_var(i,j)/Sd+(lapse*a2dDem(i,j))

			else

				matrice_var(i,j)=-9999	

			endif

			Sd=0
		 
		enddo
	enddo

	return
	end subroutine
