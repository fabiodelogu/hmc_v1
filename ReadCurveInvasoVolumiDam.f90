!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine lettura curve invaso-volume
!***********************************************************************************  
 
subroutine ReadCurveInvasoVolumiDam(sFileInvVol,iDiga)

implicit none
include 'DeclarationH.f90'

INTEGER iFlagInterp,i,iit,iDiga,itl,itll
CHARACTER*50 sFileInvVol
INTEGER*4 iLStr,iPathLenght
REAL*8 dVtemp,dHtemp,dRate

iPathLenght = iLStr(sPathTurb)



!Inserisco i dati in  modo che i valori più alti di livello e volume risultino alla coordinata
!massima della matrice. Sotto al valore inferiore di validità inserisco H e V minimi

do i=1,2 !Il primo giro calcola la lunghezza della serie di dati
    iit=0
	a1d_Level(iDiga,:)=a1d_Level(iDiga,1)
	a1d_Volume(iDiga,:)=a1d_Volume(iDiga,1)
	open(32,file=sPathTurb(1:iPathLenght)//sFileInvVol,status='old',err=14)
	read(32,*)
	read(32,*)
	read(32,*)
13	iit=iit+1
	
	read(32,*,end=113)dHtemp,dVtemp
	if(i.eq.2)then
		!Se troppo lunga la serie la campiono
		if(iit/itll.GT.0.AND.iit/itll.LE.SIZE(a1d_Volume,dim=2))then
			a1d_Level(iDiga,iit/itll)=dHtemp
			a1d_Volume(iDiga,iit/itll)=dVtemp
		endif
	endif

	go to 13
113	CLOSE(32)
    itl=iit-1
	itll=itl/SIZE(a1d_Volume,dim=2)
	itll=AINT(real(itll))+1
	if(itll.lt.1)itll=1
	write(*,*)'turb length and dl',itl,itll
	

enddo
	

IF(1.eq.2)Then !Se non c'è il file curva invaso volume la pongo -9999, nei calcoli considero andamenteo lineare
14	a1d_Level(iDiga,:)=-9999
	a1d_Volume(iDiga,:)=-9999
	write(*,*)'non trovo ',sFileInvVol,', Curvalineare'    
ENDIF


return
end subroutine