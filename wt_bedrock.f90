!***********************************************************************************   
!	Subroutine wt_bedrock
!	Calcola i Vwm per ogni pixel ed eventualmente inizializza la WT
!***********************************************************************************

!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

subroutine wt_bedrock(iRows,iCols)
implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,ii,iii,j,jj,jjj,kkk,ttt
integer a,b,perc_tot,di(iRows,iCols)

real*8 dBmin,dBmax,dem_max,dBextreme,dHbrmin,dem_min,a2Profiles(4,50),fMin,fCan,fPen,fOv

dBmin=10
dBmax=0
di=0
fMin=0
fCan=0
fPen=0
fOv=0

!Da inserire le funzioni di Fortran min e max
DO i=1,iRows 
	DO j=1,iCols
		if(a2dDem(i,j).gt.0)then
			if(a2dAlpha(i,j).lt.dBmin.and.a2dAlpha(i,j).gt.0.and.a2dAlpha(i,j).gt.0.001)dBmin=a2dAlpha(i,j)
			if(a2dAlpha(i,j).gt.dBmax.and.a2dAlpha(i,j).gt.0)dBmax=a2dAlpha(i,j)
		endif
		
	ENDDO
ENDDO

!dBmin=min(min(a2dAlpha,DIM=1,MASK=a2dDem.GT.0.0))
dem_max=maxval(maxval(a2dDem,DIM = 1),DIM = 1)
dem_min=minval(minval(a2dDem,DIM = 1,MASK=a2dDem.gt.0),DIM=1)

dHbrmin=10 !Valore minimo in mm per non mettere proprio a 0
if(dHbr<dHbrmin)dHbrmin=dHbr*0.9
!IMPORTANTE,  definisce l'angolo massimo oltre il quale non ho strato di suolo profondo
if(dSlopeMax.lt.0.01)dSlopeMax =0.01
dBextreme=3.14/180*dSlopeMax !70 orba 40 Vda
if(dBmax.gt.dBextreme)then
	dBmax=dBextreme
	dHbrmin=0
endif

WHERE(a2dDem.gt.0)
	a2dVwtMax=dHbr*(1-(dtan(a2dAlpha)-dtan(dBmin))/(dtan(dBmax)-dtan(dBmin))*(1-dHbrmin/dHbr))
ENDWHERE
WHERE(a2dDem.gt.0.and.a2dVwtMax.gt.dHbr)
	a2dVwtMax=dHbr
ENDWHERE
WHERE(a2dDem.gt.0.and.a2dVwtMax.lt.0)
	a2dVwtMax=0
ENDWHERE
WHERE(a2dDem.gt.0.and.a2dVwtMax.ge.0)
	a2dVwtMax=a2dDem-a2dVwtMax/1000
ENDWHERE
WHERE(a2dDem.gt.0.and.a2dVwtMax.lt.0.0)
	a2dVwtMax=a2dDem
ENDWHERE

if(iFlagStateVar.eq.0)then !Non ho lo stato iniziale della WT
	open(22,file='WT_initial.txt',status='old')
	read(22,*)fMin
	read(22,*)fCan
	read(22,*)fPen
	read(22,*)fOv
	close(22)    
	WHERE(a2dDem.gt.0)
		a2dVwt=(dHbr-fMin)*((dtan(a2dAlpha)-dtan(dBmin))/(dtan(dBmax)-dtan(dBmin)))+fMin !wet

		a2dVwt=a2dVwt/1000

		di=a2dVwt*1000
	ENDWHERE

	WHERE(a2dDem.gt.0)
		!a2dVwt=a2dVwt+a2dVwtMax !Ver 1
		a2dVwt=a2dDem-a2dVwt
	ENDWHERE
	!Inizializzazione aggiunta
	WHERE(a2dAlpha.gt.fPen)a2dVwt=a2dVwtMax+fOv/1000 !Riempimento WT in mm per pendenze elevate
	WHERE(a2iChoice.eq.1)a2dVwt=a2dDem-fCan/1000 !Riempimento WT in mm sotto i canali
else	!Leggo da file lo stato iniziale della WT
	WHERE(a2dDem.gt.0)
		!a2dVwt=a2dVwt+a2dVwtMax !Ver 1
		a2dVwt=a2dDem-a2dVwt/1000
	ENDWHERE
endif
!Chek limite superiore
WHERE(a2dVwt.gt.a2dDem)
	a2dVwt=a2dDem
ENDWHERE
!Chek limite inferiore
WHERE(a2dVwt.lt.a2dVwtMax)
	a2dVwt=a2dVwtMax
ENDWHERE

a=93
b=247
ttt=1
DO WHILE((a2iChoice(a,b).eq.0))

		IF((a.gt.0.and.a.le.iRows).and.(b.gt.0.and.b.le.iCols))THEN
					
			ii=int((a2iPun(a,b)-1)/3)-1
			jj=a2iPun(a,b)-5-3*ii
			iii=a+ii
			jjj=b+jj
			a2Profiles(1,ttt)=a2dDem(iii,jjj)
			a2Profiles(2,ttt)=a2dVwt(iii,jjj)
			a2Profiles(3,ttt)=a2dVwtMax(iii,jjj)
			a2Profiles(4,ttt)=a2dAlpha(iii,jjj)
			ttt=ttt+1
		ENDIF
		a=iii
		b=jjj

ENDDO

di=(a2dDem-a2dVwt)*1000
di=(a2dVwt-a2dVwtMax)*1000
return
end subroutine