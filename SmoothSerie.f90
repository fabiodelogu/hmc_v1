!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine lettura turbinati
!***********************************************************************************  
 
subroutine SmoothSerie(a1dDataT,iRank,iNsmooth)

implicit none
include 'DeclarationH.f90'

INTEGER i,iFlagInterp,it,iit,iiit,iRank,in,iNsmooth
CHARACTER*500 sF,sFileTurbinati,sName
CHARACTER*12 xstring
CHARACTER*50 sNmTmp
REAL*8 dSimLength,dNum,dMean,dNumT,dVdif
REAL*8 a1dDataT(iRank),a1dDataTmp(iRank),dSumB,dSumA

a1dDataTmp=0
a1dDataTmp=a1dDataT

iRank = SIZE (a1dDataT,dim=1)
dSumB = SUM(a1dDataT, DIM = 1)
dNumT=0
dVdif=0

DO it=1+iNsmooth,iRank-iNsmooth-1
	dNum=0
	dMean=0
	if(a1dDataT(it).gt.0.0)dNumT=dNumT+1

	DO iit=it-iNsmooth,it+iNsmooth
		if(a1dDataT(iit).ge.0.0)then
			dMean=dMean+a1dDataT(iit)
			dNum=dNum+1
		endif
	ENDDO
	if(dNum.gt.0)then
		dMean=dMean/dNum
		a1dDataTmp(it)=dMean
	endif
ENDDO

a1dDataT=a1dDataTmp
dSumA = SUM(a1dDataT, DIM = 1)
if(dNumT.gt.0)then
	dVdif=(dSumB-dSumA)/dNumT
endif
!write(*,*)(dSumA-dSumB)/dSumB,dVdif,dSumA,dSumB

if(dVdif.lt.0.0)then !Se la differenza è positiva
	WHERE(a1dDataT.gt.abs(dVdif))
		a1dDataT=a1dDataT+dVdif
	ENDWHERE
else
	WHERE(a1dDataT.gt.0)
		a1dDataT=a1dDataT+dVdif
	ENDWHERE
endif

return
end subroutine