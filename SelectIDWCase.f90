!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!*******Interpolazione dati micrometeorologici

SUBROUTINE SelectIDWCase(iRows,iCols,nstaz,latSt,longSt,zSt,val,Vmin,Vmax,label,matrice_var)

IMPLICIT NONE
INCLUDE 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i,j,ii,nstaz,label,conta
REAL*8 ar,br,rr_int,ad_lapse,nx,ny
	
REAL*8 Vmin,Vmax
REAL*8 dist(nstaz)
REAL*8 val(nstaz),valori(nstaz)
REAL*8 latSt(nstaz),longSt(nstaz),zSt(nstaz),zStazioni(nstaz)
REAL*8 matrice_var(iRows,iCols)
    
parameter(ad_lapse=-0.006)

br=0.0
ar=0.0
rr_int=0.0
valori=0.0
zStazioni=0.0
conta=0
	
!-----------------------------------------------------------------------------------
!  label = 1 , temperature station linear regression + IDW
!  label = 2 , temperature station adiabatic lapse rate +IDW
!  label = else , IDW
!-----------------------------------------------------------------------------------

DO i=1,nstaz
	IF(val(i).ge.Vmin.and.val(i).le.Vmax)then
 		conta=conta+1	    
 		valori(conta)=val(i)
 		zStazioni(conta)=zSt(i)
 	ENDIF
ENDDO

IF(conta.gt.0)THEN

	SELECT CASE (label)

	CASE (1)

! SE NON CI SONO STAZIONI SUFF. APPLICA IL LAPSE RATE ADIABATICO
	IF(conta.gt.1)THEN
		CALL LinRegr(conta,zStazioni,valori,ar,br,rr_int)
		CALL IDW(iRows,iCols,nstaz,latSt,longSt,zSt,val,Vmin,Vmax,br,matrice_var)	
	ELSE 
		CALL IDW(iRows,iCols,nstaz,latSt,longSt,zSt,val,Vmin,Vmax,ad_lapse,matrice_var)	
	ENDIF

	CASE (2)

	CALL IDW(iRows,iCols,nstaz,latSt,longSt,zSt,val,Vmin,Vmax,ad_lapse,matrice_var)

	CASE DEFAULT

	CALL IDW(iRows,iCols,nstaz,latSt,longSt,zSt,val,Vmin,Vmax,br,matrice_var)
    
	END SELECT

ELSE
	
	!matrice_var=-9999 !Se non ho dati in questo modo salta tutto!
	matrice_var=0.0

ENDIF

	
	
!-----------------------------------------------------------------------------------

RETURN
END SUBROUTINE


